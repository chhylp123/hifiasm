#include "Levenshtein_distance.h"
#include <immintrin.h>


#define init_simd_ed4(PSA, PNA, THRE, ABS_DIAG, R_ERR, R_PE, SI, TN, CUT, BD, I, MM, PEQ_MM, LZ, IBD) {\
	(R_ERR)[(SI)] = INT32_MAX; (R_PE)[(SI)] = -1; (IBD)[(SI)] = ((THRE)<<1) - (ABS_DIAG)[(SI)];\
	if(((PNA)[(SI)] <= (TN) + (CUT)) && ((TN) <= (PNA)[(SI)] + (CUT))) {\
		(BD) = (((THRE)<<1)+1)-(ABS_DIAG)[(SI)]; (BD) = (((BD)<=(PNA)[(SI)])?(BD):(PNA)[(SI)]); (LZ) |= (((int32_t)1u) << (SI));\
		for ((I) = 0, (MM) = (((Word)1)<<((ABS_DIAG)[(SI)])); (I) < (BD); (I)++) {\
			(PEQ_MM)[seq_nt4_table[(uint8_t)(PSA)[(SI)][(I)]]][(SI)] |= (MM); (MM) <<= 1;\
		}\
	}\
}

#define ed_core_64x4(PEQz, VPz, VNz, Xz, D0z, HNz, HPz) {	\
	/**(X) = (Peq)|(VN);**/\
	(Xz) = _mm256_or_si256((PEQz), (VNz)); \
	/**(D0) = (((VP) + ((X)&(VP))) ^ (VP)) | (X);**/\
	(D0z) = _mm256_or_si256(_mm256_xor_si256(_mm256_add_epi64((VPz), _mm256_and_si256((Xz), (VPz))), (VPz)), (Xz)); \
	/**(HN) = (VP)&(D0);**/\
	(HNz) = _mm256_and_si256((VPz), (D0z)); \
	/**(HP) = (VN) | ~((VP) | (D0));**/\
	(HPz) = _mm256_or_si256((VNz), _mm256_andnot_si256(_mm256_or_si256((VPz), (D0z)), _mm256_set1_epi64x(-1))); \
	/**(X) = (D0) >> 1;**/\
	(Xz) = _mm256_srli_epi64((D0z), 1); \
	/**(VN) = (X)&(HP);**/\
	(VNz) = _mm256_and_si256((Xz), (HPz)); \
	/**(VP) = (HN) | ~((X) | (HP));**/\
	(VPz) = _mm256_or_si256((HNz), _mm256_andnot_si256(_mm256_or_si256((Xz), (HPz)), _mm256_set1_epi64x(-1))); \
}

#define ed_core_upx4(PEQz, PSA, PNA, IBD, HT, CC, MMK, SI) {	\
	if((HT) & (((int32_t)1u) << (SI))) {\
		(IBD)[(SI)]++;\
		if((IBD)[(SI)] < (PNA)[(SI)]) {\
			(CC) = seq_nt4_table[(uint8_t)(PSA)[(SI)][(IBD)[(SI)]]];\
			if((CC) < 4) (PEQz)[(CC)] = _mm256_or_si256((PEQz)[(CC)], (MMK)[(SI)]);\
		}\
	}\
}

#define ed_tail_upx4(HT, SI, ST, AI, PNA, ABS_DIAG, K, ERR_MM, VP_MM, VN_MM, THRE, R_ERR, R_PE, BD, I) {\
	if((HT) & (((int32_t)1u) << (SI))) {\
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

#define ED_TAIL_LANE(K) do { \
    if ((ht) & (((int32_t)1u) << (K))) { \
        st = tn - 1 - abs_diag_a[(K)]; \
        ai = pna[(K)] - tn + abs_diag_a[(K)]; \
        for (i = 0, uge = INT64_MAX; st < 0 && i < ai; i++, st++) { \
            err_mm[(K)] += ((VP_mm[(K)] >> i) & 1ULL); \
            err_mm[(K)] -= ((VN_mm[(K)] >> i) & 1ULL); \
        } \
        if ((err_mm[(K)] <= thre) && (err_mm[(K)] <= r_err[(K)])) { \
            r_err[(K)] = err_mm[(K)]; \
            r_pe[(K)] = st; \
        } \
        st -= i; \
        while (i < ai) { \
            err_mm[(K)] += ((VP_mm[(K)] >> i) & 1ULL); \
            err_mm[(K)] -= ((VN_mm[(K)] >> i) & 1ULL); \
            ++i; \
            if ((err_mm[(K)] <= thre) && (err_mm[(K)] <= r_err[(K)])) { \
                r_err[(K)] = err_mm[(K)]; \
                r_pe[(K)] = st + i; \
            } \
            if (i == thre) uge = err_mm[(K)]; \
        } \
        if ((uge <= thre) && (uge == r_err[(K)])) r_pe[(K)] = st + thre; \
    } \
} while (0)

void ed_band_cal_semi_64_w_absent_diag_avx4(char **psa, int32_t *pna, char *tstr, int32_t tn, int32_t thre, int32_t *abs_diag_a, int64_t *r_err, int64_t *r_pe)
{
	// r_err[0] = r_err[1] = r_err[2] = r_err[3] = r_err[4] = r_err[5] = r_err[6] = r_err[7] = thre+1;
	// r_pe[0] = r_pe[1] = r_pe[2] = r_pe[3] = r_pe[4] = r_pe[5] = r_pe[6] = r_pe[7] = -1;   
	
	Word mm, Peq_mm[5][AVX_GS2] = {{0}}, *VN_mm = NULL, *VP_mm = NULL, c = 0; __m256i Peq[5], VP, VN, X, D0, HN, HP, lone, E, C, mmk[AVX_GS2]; 
	int32_t lz = 0, ht = (((int32_t)1u)<<AVX_GS2)-1; int32_t bd, ibd[AVX_GS2], i, last_high = (thre<<1), tn0 = tn - 1, cut = thre+last_high;

	lone = _mm256_set1_epi64x(1);
	VP = _mm256_setzero_si256();

	VN_mm = Peq_mm[0];
	VN_mm[0] = (((Word)1)<<(abs_diag_a[0]))-1; VN_mm[1] = (((Word)1)<<(abs_diag_a[1]))-1; VN_mm[2] = (((Word)1)<<(abs_diag_a[2]))-1; VN_mm[3] = (((Word)1)<<(abs_diag_a[3]))-1; 
	VN = _mm256_loadu_si256((const __m256i *)VN_mm);

	VN_mm[0] = abs_diag_a[0]; VN_mm[1] = abs_diag_a[1]; VN_mm[2] = abs_diag_a[2]; VN_mm[3] = abs_diag_a[3]; 
	E = _mm256_loadu_si256((const __m256i *)VN_mm);

	memset(VN_mm, 0, (sizeof((*VN_mm))*AVX_GS2)); VN_mm = NULL;///reset

	init_simd_ed4(psa, pna, thre, abs_diag_a, r_err, r_pe, 0, tn, cut, bd, i, mm, Peq_mm, lz, ibd);
	init_simd_ed4(psa, pna, thre, abs_diag_a, r_err, r_pe, 1, tn, cut, bd, i, mm, Peq_mm, lz, ibd);
	init_simd_ed4(psa, pna, thre, abs_diag_a, r_err, r_pe, 2, tn, cut, bd, i, mm, Peq_mm, lz, ibd);
	init_simd_ed4(psa, pna, thre, abs_diag_a, r_err, r_pe, 3, tn, cut, bd, i, mm, Peq_mm, lz, ibd);

	ht &= lz;
	if(ht == 0) return;

	Peq[0] = _mm256_loadu_si256((const __m256i *)Peq_mm[0]); 
	Peq[1] = _mm256_loadu_si256((const __m256i *)Peq_mm[1]); 
	Peq[2] = _mm256_loadu_si256((const __m256i *)Peq_mm[2]); 
	Peq[3] = _mm256_loadu_si256((const __m256i *)Peq_mm[3]); 
	Peq[4] = _mm256_setzero_si256();
	
	C = _mm256_set1_epi64x(cut + 1);///_mm512_set1_epi64x(cut);

	mm = ((Word)1 << (thre<<1));///for the incoming char/last char**
	mmk[0] = _mm256_set_epi64x(0, 0, 0, mm);
	mmk[1] = _mm256_set_epi64x(0, 0, mm, 0);
	mmk[2] = _mm256_set_epi64x(0, mm, 0, 0);
	mmk[3] = _mm256_set_epi64x(mm, 0, 0, 0);

	i = 0; 

	while (i < tn0) {
		ed_core_64x4(Peq[seq_nt4_table[(uint8_t)tstr[i]]], VP, VN, X, D0, HN, HP);
		E = _mm256_add_epi64(_mm256_xor_si256(lone, _mm256_and_si256(D0, lone)), E);
		// ht = _mm512_cmple_epi64_mask(E, C);
		ht =_mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpgt_epi64(C, E)));
		ht &= lz;
		if(ht == 0) return;

		Peq[0] = _mm256_srli_epi64(Peq[0], 1);
		Peq[1] = _mm256_srli_epi64(Peq[1], 1);
		Peq[2] = _mm256_srli_epi64(Peq[2], 1);
		Peq[3] = _mm256_srli_epi64(Peq[3], 1);
		i++; ///c = 4;

		ed_core_upx4(Peq, psa, pna, ibd, ht, c, mmk, 0);
		ed_core_upx4(Peq, psa, pna, ibd, ht, c, mmk, 1);
		ed_core_upx4(Peq, psa, pna, ibd, ht, c, mmk, 2);
		ed_core_upx4(Peq, psa, pna, ibd, ht, c, mmk, 3);
	}

	ed_core_64x4(Peq[seq_nt4_table[(uint8_t)tstr[i]]], VP, VN, X, D0, HN, HP);
	E = _mm256_add_epi64(_mm256_xor_si256(lone, _mm256_and_si256(D0, lone)), E);
	// ht = _mm512_cmple_epi64_mask(E, C);
	ht =_mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpgt_epi64(C, E)));
	ht &= lz;
	if(ht == 0) return;

	int32_t st, ai; int64_t err_mm[AVX_GS2], uge;
	VN_mm = Peq_mm[0]; VP_mm = Peq_mm[1];
	_mm256_storeu_si256((__m256i *)VN_mm, VN); _mm256_storeu_si256((__m256i *)VP_mm, VP); _mm256_storeu_si256((__m256i *)err_mm, E);

	ED_TAIL_LANE(0);
	ED_TAIL_LANE(1);
	ED_TAIL_LANE(2);
	ED_TAIL_LANE(3);
}
