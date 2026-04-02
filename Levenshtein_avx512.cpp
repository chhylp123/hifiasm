#include "Levenshtein_distance.h"
#include <immintrin.h>

#define init_simd_ed(PSA, PNA, THRE, ABS_DIAG, R_ERR, R_PE, SI, TN, CUT, BD, I, MM, PEQ_MM, LZ, IBD) {\
	(R_ERR)[(SI)] = INT32_MAX; (R_PE)[(SI)] = -1; (IBD)[(SI)] = ((THRE)<<1) - (ABS_DIAG)[(SI)];\
	if(((PNA)[(SI)] <= (TN) + (CUT)) && ((TN) <= (PNA)[(SI)] + (CUT))) {\
		(BD) = (((THRE)<<1)+1)-(ABS_DIAG)[(SI)]; (BD) = (((BD)<=(PNA)[(SI)])?(BD):(PNA)[(SI)]); (LZ) |= (((__mmask8)1u) << (SI));\
		for ((I) = 0, (MM) = (((Word)1)<<((ABS_DIAG)[(SI)])); (I) < (BD); (I)++) {\
			(PEQ_MM)[seq_nt4_table[(uint8_t)(PSA)[(SI)][(I)]]][(SI)] |= (MM); (MM) <<= 1;\
		}\
	}\
}

#define ed_core_64x8(PEQz, VPz, VNz, Xz, D0z, HNz, HPz) {	\
	/**(X) = (Peq)|(VN);**/\
	(Xz) = _mm512_or_si512((PEQz), (VNz));\
	/**(D0) = (((VP) + ((X)&(VP))) ^ (VP)) | (X);**/\
	(D0z) = _mm512_or_si512(_mm512_xor_si512(_mm512_add_epi64((VPz), _mm512_and_si512((Xz), (VPz))), (VPz)), (Xz));\
	/**(HN) = (VP)&(D0);**/\
	(HNz) = _mm512_and_si512((VPz), (D0z));\
	/**(HP) = (VN) | ~((VP) | (D0));**/\
	(HPz) = _mm512_or_si512((VNz), _mm512_andnot_si512(_mm512_or_si512((VPz), (D0z)), _mm512_set1_epi64(-1)));\
	/**(X) = (D0) >> 1;**/\
	(Xz)  = _mm512_srli_epi64((D0z), 1);\
	/**(VN) = (X)&(HP);**/\
	(VNz) = _mm512_and_si512((Xz), (HPz));\
	/**(VP) = (HN) | ~((X) | (HP));**/\
	(VPz) = _mm512_or_si512((HNz), _mm512_andnot_si512(_mm512_or_si512((Xz), (HPz)), _mm512_set1_epi64(-1)));\
}

#define ed_core_upx8(PEQz, PSA, PNA, IBD, HT, CC, MMK, SI) {	\
	if((HT) & (((__mmask8)1u) << (SI))) {\
		(IBD)[(SI)]++;\
		if((IBD)[(SI)] < (PNA)[(SI)]) {\
			(CC) = seq_nt4_table[(uint8_t)(PSA)[(SI)][(IBD)[(SI)]]];\
			if((CC) < 4) (PEQz)[(CC)] = _mm512_or_si512((PEQz)[(CC)], (MMK)[(SI)]);\
		}\
	}\
}

#define ed_tail_upx8(HT, SI, ST, AI, PNA, ABS_DIAG, K, ERR_MM, VP_MM, VN_MM, THRE, R_ERR, R_PE, BD, I) {\
	if((HT) & (((__mmask8)1u) << (SI))) {\
		(ST)[(SI)] -= (ABS_DIAG)[(SI)]; (AI)[(SI)] += (PNA)[(SI)] + (ABS_DIAG)[(SI)];\
		for ((K)[(SI)] = 0; (ST)[(SI)] < 0 && (K)[(SI)] < (AI)[(SI)]; (K)[(SI)]++, (ST)[(SI)]++) {\
			(ERR_MM)[(SI)] += ((VP_MM)[(SI)]&(1ULL)); (VP_MM)[(SI)]>>=1;\
			(ERR_MM)[(SI)] -= ((VN_MM)[(SI)]&(1ULL)); (VN_MM)[(SI)]>>=1;\
		}\
		if (((ERR_MM)[(SI)] <= (THRE)) && ((ERR_MM)[(SI)] <= (R_ERR)[(SI)])) {\
			(R_ERR)[(SI)] = (ERR_MM)[(SI)]; (R_PE)[(SI)] = (ST)[(SI)];\
		}\
		(ST)[(SI)] -= (K)[(SI)]; (BD)++; (I) = (SI);\
	}\
}

#define ed_tail_ck8(MBEST, SI, R_PE, ST, K, THRE, UGE_MM, ERR_MM, AI, HT) {\
	(K)[(SI)]++;\
	if((MBEST) & (((__mmask8)1u) << (SI))) {\
		(R_PE)[(SI)] = (ST)[(SI)] + (K)[(SI)];\
	}\
	if((K)[(SI)] >= (AI)[(SI)]) (HT) &= ~(((__mmask8)1u) << (SI));\
	if((K)[(SI)] == (THRE)) (UGE_MM)[(SI)] = (ERR_MM)[(SI)];\
}


void ed_band_cal_semi_64_w_absent_diag_avx8(char **psa, int32_t *pna, char *tstr, int32_t tn, int32_t thre, int32_t *abs_diag_a, int64_t *r_err, int64_t *r_pe)
{
	// r_err[0] = r_err[1] = r_err[2] = r_err[3] = r_err[4] = r_err[5] = r_err[6] = r_err[7] = thre+1;
	// r_pe[0] = r_pe[1] = r_pe[2] = r_pe[3] = r_pe[4] = r_pe[5] = r_pe[6] = r_pe[7] = -1;
	/**   
	ed_band_cal_semi_64_w_absent_diag_avx4(psa, pna, tstr, tn, thre, abs_diag_a, r_err, r_pe);
	ed_band_cal_semi_64_w_absent_diag_avx4(psa + 4, pna + 4, tstr, tn, thre, abs_diag_a + 4, r_err + 4, r_pe + 4);
	return;
	**/

	
	Word mm, Peq_mm[5][AVX_GS] = {{0}}, *VN_mm = NULL, *VP_mm = NULL, c = 0; __m512i Peq[5], VP, VN, X, D0, HN, HP, lone, E, C, bestE, bestPE, curPE, cutPE, threPE, ugE, mmk[AVX_GS]; 
	__mmask8 lz = ((__mmask8)0u), ht = (((__mmask8)1u)<<AVX_GS)-1, mtf, mt, mbest; int32_t bd, ibd[AVX_GS], i, last_high = (thre<<1), tn0 = tn - 1, cut = thre+last_high;

	lone = _mm512_set1_epi64(1);
	VP = _mm512_setzero_si512();

	VN_mm = Peq_mm[0];
	VN_mm[0] = (((Word)1)<<(abs_diag_a[0]))-1; VN_mm[1] = (((Word)1)<<(abs_diag_a[1]))-1; VN_mm[2] = (((Word)1)<<(abs_diag_a[2]))-1; VN_mm[3] = (((Word)1)<<(abs_diag_a[3]))-1; 
	VN_mm[4] = (((Word)1)<<(abs_diag_a[4]))-1; VN_mm[5] = (((Word)1)<<(abs_diag_a[5]))-1; VN_mm[6] = (((Word)1)<<(abs_diag_a[6]))-1; VN_mm[7] = (((Word)1)<<(abs_diag_a[7]))-1; 
	VN = _mm512_loadu_si512(VN_mm);

	VN_mm[0] = abs_diag_a[0]; VN_mm[1] = abs_diag_a[1]; VN_mm[2] = abs_diag_a[2]; VN_mm[3] = abs_diag_a[3]; 
	VN_mm[4] = abs_diag_a[4]; VN_mm[5] = abs_diag_a[5]; VN_mm[6] = abs_diag_a[6]; VN_mm[7] = abs_diag_a[7];
	E = _mm512_loadu_si512(VN_mm);

	memset(VN_mm, 0, (sizeof((*VN_mm))*AVX_GS)); VN_mm = NULL;///reset

	init_simd_ed(psa, pna, thre, abs_diag_a, r_err, r_pe, 0, tn, cut, bd, i, mm, Peq_mm, lz, ibd);
	init_simd_ed(psa, pna, thre, abs_diag_a, r_err, r_pe, 1, tn, cut, bd, i, mm, Peq_mm, lz, ibd);
	init_simd_ed(psa, pna, thre, abs_diag_a, r_err, r_pe, 2, tn, cut, bd, i, mm, Peq_mm, lz, ibd);
	init_simd_ed(psa, pna, thre, abs_diag_a, r_err, r_pe, 3, tn, cut, bd, i, mm, Peq_mm, lz, ibd);
	init_simd_ed(psa, pna, thre, abs_diag_a, r_err, r_pe, 4, tn, cut, bd, i, mm, Peq_mm, lz, ibd);
	init_simd_ed(psa, pna, thre, abs_diag_a, r_err, r_pe, 5, tn, cut, bd, i, mm, Peq_mm, lz, ibd);
	init_simd_ed(psa, pna, thre, abs_diag_a, r_err, r_pe, 6, tn, cut, bd, i, mm, Peq_mm, lz, ibd);
	init_simd_ed(psa, pna, thre, abs_diag_a, r_err, r_pe, 7, tn, cut, bd, i, mm, Peq_mm, lz, ibd);

	ht &= lz;
	if(ht == 0) return;

	Peq[0] = _mm512_loadu_si512(Peq_mm[0]); 
	Peq[1] = _mm512_loadu_si512(Peq_mm[1]); 
	Peq[2] = _mm512_loadu_si512(Peq_mm[2]); 
	Peq[3] = _mm512_loadu_si512(Peq_mm[3]); 
	Peq[4] = _mm512_setzero_si512();

	C = _mm512_set1_epi64(cut);

	mm = ((Word)1 << (thre<<1));///for the incoming char/last char**
	mmk[0] = _mm512_mask_set1_epi64(VP, 1, mm);
	mmk[1] = _mm512_mask_set1_epi64(VP, 2, mm);
	mmk[2] = _mm512_mask_set1_epi64(VP, 4, mm);
	mmk[3] = _mm512_mask_set1_epi64(VP, 8, mm);
	mmk[4] = _mm512_mask_set1_epi64(VP, 16, mm);
	mmk[5] = _mm512_mask_set1_epi64(VP, 32, mm);
	mmk[6] = _mm512_mask_set1_epi64(VP, 64, mm);
	mmk[7] = _mm512_mask_set1_epi64(VP, 128, mm);

	i = 0; 

	while (i < tn0) {
		ed_core_64x8(Peq[seq_nt4_table[(uint8_t)tstr[i]]], VP, VN, X, D0, HN, HP);
		E = _mm512_add_epi64(_mm512_xor_si512(lone, _mm512_and_si512(D0, lone)), E);
		ht = _mm512_cmple_epi64_mask(E, C);
		ht &= lz;
		if(ht == 0) return;

		Peq[0] = _mm512_srli_epi64(Peq[0], 1);
		Peq[1] = _mm512_srli_epi64(Peq[1], 1);
		Peq[2] = _mm512_srli_epi64(Peq[2], 1);
		Peq[3] = _mm512_srli_epi64(Peq[3], 1);
		i++; ///c = 4;

		ed_core_upx8(Peq, psa, pna, ibd, ht, c, mmk, 0);
		ed_core_upx8(Peq, psa, pna, ibd, ht, c, mmk, 1);
		ed_core_upx8(Peq, psa, pna, ibd, ht, c, mmk, 2);
		ed_core_upx8(Peq, psa, pna, ibd, ht, c, mmk, 3);
		ed_core_upx8(Peq, psa, pna, ibd, ht, c, mmk, 4);
		ed_core_upx8(Peq, psa, pna, ibd, ht, c, mmk, 5);
		ed_core_upx8(Peq, psa, pna, ibd, ht, c, mmk, 6);
		ed_core_upx8(Peq, psa, pna, ibd, ht, c, mmk, 7);
	}

	ed_core_64x8(Peq[seq_nt4_table[(uint8_t)tstr[i]]], VP, VN, X, D0, HN, HP);
	E = _mm512_add_epi64(_mm512_xor_si512(lone, _mm512_and_si512(D0, lone)), E);
	ht = _mm512_cmple_epi64_mask(E, C);
	ht &= lz;
	if(ht == 0) return;


	// site = tn - 1 - abs_diag;/**up bound**/
	// ai = pn - tn + abs_diag; /**in most cases, ai = (thre<<1)**/
	int32_t st[AVX_GS] = {tn-1, tn-1, tn-1, tn-1, tn-1, tn-1, tn-1, tn-1}; 
	int32_t ai[AVX_GS] = {-tn, -tn, -tn, -tn, -tn, -tn, -tn, -tn};
	int32_t k[AVX_GS] = {0}; i = -1; bd = 0;
	int64_t err_mm[AVX_GS], uge_mm[AVX_GS] = {INT32_MAX, INT32_MAX, INT32_MAX, INT32_MAX, INT32_MAX, INT32_MAX, INT32_MAX, INT32_MAX}; 
	VN_mm = Peq_mm[0]; VP_mm = Peq_mm[1];
	_mm512_storeu_si512(VN_mm, VN); _mm512_storeu_si512(VP_mm, VP); _mm512_storeu_si512(err_mm, E);

	ed_tail_upx8(ht, 0, st, ai, pna, abs_diag_a, k, err_mm, VP_mm, VN_mm, thre, r_err, r_pe, bd, i);
	ed_tail_upx8(ht, 1, st, ai, pna, abs_diag_a, k, err_mm, VP_mm, VN_mm, thre, r_err, r_pe, bd, i);
	ed_tail_upx8(ht, 2, st, ai, pna, abs_diag_a, k, err_mm, VP_mm, VN_mm, thre, r_err, r_pe, bd, i);
	ed_tail_upx8(ht, 3, st, ai, pna, abs_diag_a, k, err_mm, VP_mm, VN_mm, thre, r_err, r_pe, bd, i);
	ed_tail_upx8(ht, 4, st, ai, pna, abs_diag_a, k, err_mm, VP_mm, VN_mm, thre, r_err, r_pe, bd, i);
	ed_tail_upx8(ht, 5, st, ai, pna, abs_diag_a, k, err_mm, VP_mm, VN_mm, thre, r_err, r_pe, bd, i);
	ed_tail_upx8(ht, 6, st, ai, pna, abs_diag_a, k, err_mm, VP_mm, VN_mm, thre, r_err, r_pe, bd, i);
	ed_tail_upx8(ht, 7, st, ai, pna, abs_diag_a, k, err_mm, VP_mm, VN_mm, thre, r_err, r_pe, bd, i);

	if(bd <= 0) return;

	if(bd > 1) {
		VN = _mm512_loadu_si512(VN_mm); VP = _mm512_loadu_si512(VP_mm); E = _mm512_loadu_si512(err_mm); i = 0;

		bestE = _mm512_loadu_si512(r_err); ///threE = _mm512_set1_epi64(thre);
		if(k[0] >= ai[0]) ht &= ((__mmask8)(255-1));
		if(k[1] >= ai[1]) ht &= ((__mmask8)(255-2));
		if(k[2] >= ai[2]) ht &= ((__mmask8)(255-4));
		if(k[3] >= ai[3]) ht &= ((__mmask8)(255-8));
		if(k[4] >= ai[4]) ht &= ((__mmask8)(255-16));
		if(k[5] >= ai[5]) ht &= ((__mmask8)(255-32));
		if(k[6] >= ai[6]) ht &= ((__mmask8)(255-64));
		if(k[7] >= ai[7]) ht &= ((__mmask8)(255-128));

		err_mm[0] = r_pe[0]; err_mm[1] = r_pe[1]; err_mm[2] = r_pe[2]; err_mm[3] = r_pe[3]; 
		err_mm[4] = r_pe[4]; err_mm[5] = r_pe[5]; err_mm[6] = r_pe[6]; err_mm[7] = r_pe[7]; 
		bestPE = _mm512_loadu_si512(err_mm);
		err_mm[0] = st[0] + k[0]; err_mm[1] = st[1] + k[1]; err_mm[2] = st[2] + k[2]; err_mm[3] = st[3] + k[3];
		err_mm[4] = st[4] + k[4]; err_mm[5] = st[5] + k[5]; err_mm[6] = st[6] + k[6]; err_mm[7] = st[7] + k[7];
		curPE = _mm512_loadu_si512(err_mm);	
		err_mm[0] = st[0] + thre; err_mm[1] = st[1] + thre; err_mm[2] = st[2] + thre; err_mm[3] = st[3] + thre;
		err_mm[4] = st[4] + thre; err_mm[5] = st[5] + thre; err_mm[6] = st[6] + thre; err_mm[7] = st[7] + thre;
		threPE = _mm512_loadu_si512(err_mm);
		err_mm[0] = st[0] + ai[0]; err_mm[1] = st[1] + ai[1]; err_mm[2] = st[2] + ai[2]; err_mm[3] = st[3] + ai[3];
		err_mm[4] = st[4] + ai[4]; err_mm[5] = st[5] + ai[5]; err_mm[6] = st[6] + ai[6]; err_mm[7] = st[7] + ai[7];
		cutPE = _mm512_loadu_si512(err_mm);

		ugE = _mm512_loadu_si512(uge_mm);

		// mtf = _mm512_cmpge_epi64_mask(curPE, threPE) | ((__mmask8)(~ht));
		mtf = _mm512_cmpge_epi64_mask(curPE, threPE);

		while ((ht != 0) && ((mtf|((__mmask8)(~ht))) != (__mmask8)255)) {
			E = _mm512_add_epi64(E, _mm512_and_si512(VP, lone)); VP = _mm512_srli_epi64(VP, 1);
			E = _mm512_sub_epi64(E, _mm512_and_si512(VN, lone)); VN = _mm512_srli_epi64(VN, 1);
			// i++;
			
			curPE = _mm512_add_epi64(curPE, lone);
			ht &= _mm512_cmple_epi64_mask(curPE, cutPE);
			if (ht == 0) break;

			mbest = _mm512_cmple_epi64_mask(E, bestE) & ht;
			
			bestE = _mm512_mask_mov_epi64(bestE, mbest, E);
			bestPE = _mm512_mask_mov_epi64(bestPE, mbest, curPE);

			mt = _mm512_cmpeq_epi64_mask(curPE, threPE);
			ugE = _mm512_mask_mov_epi64(ugE, mt&ht, E);

			mtf |= mt;
	
			// if(mbest && i < thre) _mm512_storeu_si512(err_mm, E);

			// ed_tail_ck8(mbest, 0, r_pe, st, k, thre, uge_mm, err_mm, ai, ht);
			// ed_tail_ck8(mbest, 1, r_pe, st, k, thre, uge_mm, err_mm, ai, ht); 
			// ed_tail_ck8(mbest, 2, r_pe, st, k, thre, uge_mm, err_mm, ai, ht);
			// ed_tail_ck8(mbest, 3, r_pe, st, k, thre, uge_mm, err_mm, ai, ht);
			// ed_tail_ck8(mbest, 4, r_pe, st, k, thre, uge_mm, err_mm, ai, ht);
			// ed_tail_ck8(mbest, 5, r_pe, st, k, thre, uge_mm, err_mm, ai, ht);
			// ed_tail_ck8(mbest, 6, r_pe, st, k, thre, uge_mm, err_mm, ai, ht);
			// ed_tail_ck8(mbest, 7, r_pe, st, k, thre, uge_mm, err_mm, ai, ht);
		}
		

		while (ht != 0) {
			E = _mm512_add_epi64(E, _mm512_and_si512(VP, lone)); VP = _mm512_srli_epi64(VP, 1);
			E = _mm512_sub_epi64(E, _mm512_and_si512(VN, lone)); VN = _mm512_srli_epi64(VN, 1);
			// i++;
			
			curPE = _mm512_add_epi64(curPE, lone);
			ht &= _mm512_cmple_epi64_mask(curPE, cutPE);
			if (ht == 0) break;

			mbest = _mm512_cmple_epi64_mask(E, bestE) & ht;
			
			bestE = _mm512_mask_mov_epi64(bestE, mbest, E);
			bestPE = _mm512_mask_mov_epi64(bestPE, mbest, curPE);
		}

		cutPE = _mm512_set1_epi64(thre);
		ht = _mm512_cmpgt_epi64_mask(bestE, cutPE);
		bestE = _mm512_mask_set1_epi64(bestE, ht, INT32_MAX);
		bestPE = _mm512_mask_set1_epi64(bestPE, ht, -1);

		ht = _mm512_cmple_epi64_mask(ugE, cutPE) & _mm512_cmpeq_epi64_mask(ugE, bestE);
		bestPE = _mm512_mask_mov_epi64(bestPE, ht, threPE);

		_mm512_storeu_si512(r_err, bestE);
		_mm512_storeu_si512(r_pe, bestPE);
		_mm512_storeu_si512(uge_mm, ugE);
	} else {///bd == 1
		while (k[i] < ai[i]) {
			err_mm[i] += (VP_mm[i]&(1ULL)); VP_mm[i]>>=1; 
			err_mm[i] -= (VN_mm[i]&(1ULL)); VN_mm[i]>>=1;
			++k[i];
			if ((err_mm[i] <= thre) && (err_mm[i] <= r_err[i])) {
				r_err[i] = err_mm[i]; r_pe[i] = st[i] + k[i];
			}
			if(k[i] == thre) uge_mm[i] = err_mm[i];
		}
		if((uge_mm[i] <= thre) && (uge_mm[i] == r_err[i])) r_pe[i] = st[i] + thre;
	}
}
