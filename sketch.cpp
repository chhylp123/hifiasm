#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "kvec.h"
#include "htab.h"
#include "ksort.h"
#include "Correct.h"

#define MAX_HIGH_OCC     8   // TODO: don't hard code if we need to tune this parameter
#define MAX_MAX_HIGH_OCC 16
#define GMC(a, x,y,xn) ((a)[(x)*(xn)+(y)])
#define GL(x, i) ((int64_t)((uint32_t)((x).a[(i)])))
#define A_M(p, i) ((i) >= 0 && (p).a[(i)].rid > 0)

void debug_refine(ha_mz1_t *ma, uint64_t *mmt, int32_t sn, int32_t n, int32_t m, int32_t end)
{
	uint64_t ks = end;
	int64_t t = 0, i, k, sp = -1, ep = -1, ovlp, tot = mmt[end]&0xffffffff, nt = 0;;
	while (ks != 0xffffffff)
	{
		i = ks/m; k = ks%m;
		ks = mmt[ks]>>32;
		if(ks == 0xffffffff || (int32_t)(ks/m) == (i-1))
		{
			t++;
			ovlp = ((MIN(ep, (int64_t)ma[k].pos) >= MAX(sp, (int64_t)(ma[k].pos+1-ma[k].span)))? 
						MIN(ep, (int64_t)ma[k].pos) - MAX(sp, (int64_t)(ma[k].pos+1-ma[k].span)) + 1:0);
			if(ovlp != 0) fprintf(stderr, "ERROR-OVLP\n");
			if(sp == -1 || sp > (ma[k].pos+1-ma[k].span)) sp = ma[k].pos+1-ma[k].span;
			if(ep == -1 || ep < ma[k].pos) ep = ma[k].pos;
			nt += (ma[k].rid);
		}
	}
	if(t != sn) fprintf(stderr, "ERROR-TN, t: %ld, sn: %d\n", t, sn);
	if(nt != tot) fprintf(stderr, "ERROR-TOT, nt: %ld, tot: %ld\n", nt, tot);
}

void dbg_boundary(ha_mz1_v *p, st_mt_t *mt, int32_t w, int32_t k, int32_t tot_l)
{
	if(tot_l < w + k -1) return;
	int32_t i, m, n = p->n, s, a;
	for (i = 0; i < n; i++){
		if(GL(*mt, i) >= w+k-1){
			for (m = s = a = 0; m <= i; m++){
				if(!A_M(*p, m)) continue;
				if(GL(*mt, m) <= w+k-1){
					a++;
					if(mt->a[m]&0x100000000) s++;
				}
			}
			if(a > 0 && s == 0){
				fprintf(stderr, "\nERROR1, s: %d, n: %d, tot_l: %d, end_l: %ld\n", s, n, tot_l, GL(*mt, i));
				for (m = s = a = 0; m <= i; m++){
					if(!A_M(*p, m)) continue;
					if(GL(*mt, m) <= w+k-1){
						fprintf(stderr, "lp: %ld\n", GL(*mt, m));
						a++;
						if(mt->a[m]&0x100000000) s++;
					}
				}
				
			} 
			
			break;
		}
	}
	if(i == n){
		for (m = s = a = 0; m < n; m++){
			if(!A_M(*p, m)) continue;
			if(GL(*mt, m) <= w+k-1){
				a++;
				if(mt->a[m]&0x100000000) s++;
			}
		}
		if(a > 0 && s == 0) fprintf(stderr, "ERROR2\n");
	}

	for (i = n-1; i >= 0; i--)
	{
		if (GL(*mt, i) + w <= tot_l + 1) {
			for (m = i, s = a = 0; m < n; m++){
				if(!A_M(*p, m)) continue;
				if(GL(*mt, m) + w >= tot_l + 1){
					a++;
					if(mt->a[m]&0x100000000) s++;
				} 				
			}
			if(a > 0 && s == 0) {
				fprintf(stderr, "\nERROR3, s: %d, n: %d, tot_l: %d, end_l: %ld\n", s, n, tot_l, GL(*mt, i));
				for (m = i, s = a = 0; m < n; m++){
					if(!A_M(*p, m)) continue;
					if(GL(*mt, m) + w >= tot_l + 1){
						fprintf(stderr, "lp: %ld\n", GL(*mt, m));
						a++;
						if(mt->a[m]&0x100000000) s++;
					} 				
				}
			}
			
			break;
		}
	}
	
	if(i < 0){
		for (m = s = a = 0; m < n; m++){
			if(!A_M(*p, m)) continue;
			if(GL(*mt, m) + w >= tot_l + 1){
				a++;
				if(mt->a[m]&0x100000000) s++;
			} 				
		}
		if(a > 0 && s == 0) fprintf(stderr, "ERROR4\n");
	}
}

void debug_pl(const char *str, int len, int w, int k, int is_hpc, ha_mz1_v *p, const void *hf, st_mt_t *mt)
{
	int i, l, dbi, dbcnt = 0, kmer_span = 0;
	tiny_queue_t tq;
	memset(&tq, 0, sizeof(tiny_queue_t));
	uint64_t shift1 = k - 1, mask = (1ULL<<k) - 1, kmer[4] = {0,0,0,0};
	
    for (i = l = dbi = 0; i < len; ++i) {
        int c = seq_nt4_table[(uint8_t)str[i]];
        if (c < 4) { // not an ambiguous base
            int z;
            if (is_hpc) {
                int skip_len = 1;
                if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
                    for (skip_len = 2; i + skip_len < len; ++skip_len)
                        if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
                            break;
                    i += skip_len - 1; // put $i at the end of the current homopolymer run
                }
                tq_push(&tq, skip_len);
                kmer_span += skip_len;
                ///how many bases that are covered by this HPC k-mer
                ///kmer_span includes at most k HPC elements
                if (tq.count > k) kmer_span -= tq_shift(&tq);
            } else kmer_span = l + 1 < k? l + 1 : k; 
            ///kmer_span should be used for HPC k-mer
            ///non-HPC k-mer, kmer_span should be k
            ///kmer_span is used to calculate anchor pos on reverse complementary strand

            kmer[0] = (kmer[0] << 1 | (c&1))  & mask;                  // forward k-mer
            kmer[1] = (kmer[1] << 1 | (c>>1)) & mask;
            kmer[2] = kmer[2] >> 1 | (uint64_t)(1 - (c&1))  << shift1; // reverse k-mer
            kmer[3] = kmer[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift1;
            if (kmer[1] == kmer[3]) continue; // skip "symmetric k-mers" as we don't know it strand
            z = kmer[1] < kmer[3]? 0 : 1; // strand
            ++l;

            if (l >= k && kmer_span < 256) {
                uint64_t y;
                int32_t cnt;
                y = yak_hash64_64(kmer[z<<1|0]) + yak_hash64_64(kmer[z<<1|1]);
                cnt = hf? ha_ft_cnt(hf, y) : 0;

				for (dbi = 0; dbi < mt->n; dbi++)
				{
					if(p->a[dbi].x == y && p->a[dbi].rid == cnt && p->a[dbi].pos == i && p->a[dbi].rev == z && p->a[dbi].span == kmer_span)
					{
						if(l != (int)mt->a[dbi]) fprintf(stderr, "ERROR\n");
						dbcnt++;
					}
				}
            }
        } else l = 0, tq.count = tq.front = 0, kmer_span = 0;
	}

	if(dbcnt != mt->n) fprintf(stderr, "ERROR\n");
	if(mt->n != (int)p->n) fprintf(stderr, "ERROR\n");
	for (dbi = 1; dbi < mt->n; dbi++)
	{
		if(p->a[dbi].pos <= p->a[dbi-1].pos || (int)mt->a[dbi] <= (int)mt->a[dbi-1])
		{
			fprintf(stderr, "ERROR\n");
		}
	}
}

static inline int mz1_mzcmp(const ha_mz1_t *a, const ha_mz1_t *b){return a->rid < b->rid? -1 : a->rid > b->rid? 1 : ((a->x > b->x) - (a->x < b->x));}
#define mz1_mz_lt(a, b) (mz1_mzcmp(&(a), &(b)) < 0)
KSORT_INIT(mz1_mz, ha_mz1_t, mz1_mz_lt)

static inline int mz2_mzcmp(const ha_mzl_t *a, const ha_mzl_t *b){return a->rid < b->rid? -1 : a->rid > b->rid? 1 : ((a->x > b->x) - (a->x < b->x));}
#define mz2_mz_lt(a, b) (mz2_mzcmp(&(a), &(b)) < 0)
KSORT_INIT(mz2_mz, ha_mzl_t, mz2_mz_lt)


#define HA_SC_INIT(sf, HType, VType, RidBits, PosBits)\
inline void sf##_hf_select(VType *p, int32_t si, int32_t ei, int32_t n, int32_t len, int32_t sample_dist, HType *b, int32_t force)\
{\
	if(ei - si <= 1) return;\
	int32_t ps = si < 0? 0 : p->a[si].pos;\
	int32_t pe = ei == n? len : p->a[ei].pos;\
	int32_t j, k, st = si + 1, en = ei;\
	int32_t max_high_occ = (int32_t)((double)(pe - ps) / sample_dist + .499);\
	if (max_high_occ > MAX_MAX_HIGH_OCC)\
		max_high_occ = MAX_MAX_HIGH_OCC;\
	for (j = st, k = 0; j < en && k < max_high_occ; ++j, ++k)\
		b[k] = p->a[j], b[k].pos = j; /** b[].pos keeps the index in p->a[]**/\
	ks_heapmake_##sf##_mz(k, b); /** initialize the binomial heap**/\
	for (; j < en; ++j) { /** if there are more, choose top max_high_occ**/\
		if (sf##_mz_lt(p->a[j], b[0])) { /** then update the heap**/\
			b[0] = p->a[j], b[0].pos = j;\
			ks_heapdown_##sf##_mz(0, k, b);\
		}\
	}\
	/**ks_heapsort_mz(k, b); // sorting is not needed for now**/\
	for (j = 0; j < k; ++j)\
		if (b[j].rid < pe - ps || force)\
			p->a[b[j].pos].rid = 0;\
}\
static inline int sf##_mzcmp_l(const VType *p, int32_t ai, int32_t bi)\
{\
	if(ai >= 0 && bi >= 0){\
		HType *a = &(p->a[ai]), *b = &(p->a[bi]);\
		if(a->rid > 0 && b->rid > 0) return sf##_mzcmp(a, b);\
		return (a->rid == 0) - (b->rid == 0);\
	}\
	return (ai < 0) - (bi < 0);\
}\
int32_t sf##_qfw(VType *p, st_mt_t *mt, int32_t n, int32_t tot_l, int32_t ws, int32_t i, int32_t *mi)\
{\
	int32_t m, si;\
	for (si = i, (*mi) = -1; i < n; i++){\
		if(GL(*mt, i) >= ws || (i+1 < n && GL(*mt, i) < ws && GL(*mt, i+1) > ws) || \
				(i+1 == n && tot_l >= ws && GL(*mt, i) < ws)){\
			for (m = si; m <= i; m++){\
				if(!A_M(*p, m)) continue;\
				if(sf##_mzcmp_l(p, *mi, m) >= 0) (*mi) = m;\
			}\
			if((*mi) >= 0 && A_M(*p, *mi)){\
				for (m = si; m <= i; m++){\
					if(!A_M(*p, m)) continue;\
					if(sf##_mzcmp_l(p, *mi, m) == 0) mt->a[m] |= 0x100000000;\
				}\
			}\
			break;\
		}\
	}\
	return i;\
}\
static void sf##_select_mz_h(VType *p, st_mt_t *mt, int len, int sample_dist, int32_t w, int32_t k, int32_t tot_l)\
{ /**for high-occ minimizers, choose up to max_high_occ in each high-occ streak**/\
	int32_t i, mi = -1, si, last0 = -1, n = (int32_t)p->n, m = 0, ws = w + k - 1;\
	if (n == 0) return;\
	assert((int64_t)(n) < (int64_t)((((uint64_t)1)<<PosBits)));\
	for (i = m = 0, last0 = -1; i <= n; ++i) {\
		if (i == n || p->a[i].rid == 0) {\
			if (i - last0 > 1) {\
				int32_t ps = last0 < 0? 0 : p->a[last0].pos;\
				int32_t pe = i == n? len : p->a[i].pos;\
				if(((int32_t)((double)(pe - ps) / sample_dist + .499)) > 0){\
					last0 = -2;\
					m++;\
					break;\
				}\
			}\
			last0 = i;\
		}\
	}\
	if (m == 0) return; /**no high-frequency k-mers; do nothing**/\
	if(last0 >= -1) goto sf##_ff;\
	i = 0;\
	i = sf##_qfw(p, mt, n, tot_l, ws, i, &mi);\
	if(i == n) goto sf##_ff;\
	for (si = 0, i++; i < n; i++){\
		for (; si < i; si++){\
			if(GL(*mt, si) + w > GL(*mt, i)) break;\
		}\
		/**a new minimum; then write the old min**/\
		if(sf##_mzcmp_l(p, i, mi) <= 0) {\
			if(A_M(*p, mi)) mt->a[mi] |= 0x100000000;\
			mi = i;\
		}/**old min has moved outside the window**/\
		else if(si > mi){\
			if(A_M(*p, mi)) mt->a[mi] |= 0x100000000;\
			for (m = si, mi = -1; m <= i; m++){\
				if(sf##_mzcmp_l(p, mi, m) >= 0) mi = m;\
            }\
			if(A_M(*p, mi)){\
				for (m = si; m <= i; m++){\
					if(!A_M(*p, m)) continue;\
					if(sf##_mzcmp_l(p, mi, m) == 0) mt->a[m] |= 0x100000000;\
				}\
			}\
		}\
	}\
	if(A_M(*p, mi)) mt->a[mi] |= 0x100000000;\
	for (i = n - 1; si < n && GL(*mt, si) + w <= tot_l + 1; si++){\
		if(si > mi){\
			if(A_M(*p, mi)) mt->a[mi] |= 0x100000000;\
			for (m = si, mi = -1; m <= i; m++){\
				if(sf##_mzcmp_l(p, mi, m) >= 0) mi = m;\
            }\
			if(A_M(*p, mi)){\
				for (m = si; m <= i; m++){\
					if(!A_M(*p, m)) continue;\
					if(sf##_mzcmp_l(p, mi, m) == 0) mt->a[m] |= 0x100000000;\
				}\
			}\
		}\
	}\
	/**dbg_boundary(p, mt, w, k, tot_l);**/\
	HType b[MAX_MAX_HIGH_OCC];\
	for (i = 0, last0 = -1; i <= n; ++i) {\
		if (i == n || p->a[i].rid == 0) {\
			if (i - last0 > 1) {\
				int32_t ps = last0 < 0? 0 : p->a[last0].pos;\
				int32_t pe = i == n? len : p->a[i].pos;\
				if(((int32_t)((double)(pe - ps) / sample_dist + .499)) > 0){\
					for (m = last0 + 1, mi = 0; m < i; ++m){\
						if(mt->a[m]&0x100000000) p->a[m].rid = 0, mi++;\
					}\
					if(mi == 0) sf##_hf_select(p, last0, i, n, len, sample_dist, b, 0);\
				}\
			}\
			last0 = i;\
		}\
	}\
	sf##_ff:\
	for (i = n = 0; i < (int32_t)p->n; ++i) /**squeeze out filtered minimizers**/\
		if (p->a[i].rid == 0)\
			p->a[n++] = p->a[i];\
	p->n = n;\
}\
void sf##_refine_select(VType *mz, int32_t sidx, int32_t eidx, int32_t sn, int32_t min_freq, st_mt_t *mm, int32_t *rsi, int32_t *rei)\
{\
	int32_t n = sn, m = eidx + 1 - sidx, i, k, t, mk=-1;\
	uint64_t ix, kx, ks;\
	kv_resize(uint64_t, *mm, mm->n+n*m);\
	HType *ma = mz->a + sidx;\
	uint64_t *mmt = mm->a + mm->n;\
	/**fprintf(stderr, "[M::%s::] ==> +n: %d, m: %d, sn: %d, sidx: %d, eidx: %d\n", __func__, n, m, sn, sidx, eidx);**/\
	for (i = 0; i < n; i++) /**how many selected minimizers**/\
	{\
		for (k = 0, mk = -1; k < m; k++) /**how many minimizers in total**/\
		{\
			if((int32_t)(ma[k].rid)<min_freq) continue;\
			ks = ma[k].pos + 1 - ma[k].span; t = -1;\
			if(i > 0){\
				for (t = k-1; t >= 0 && (ma[t].pos >= ks||(int32_t)(ma[t].rid)<min_freq); t--);\
			}\
			ix = (i <= 0?0:(t<0?0xffffffff:(GMC(mmt, i-1,t,m)&0xffffffff)));\
			if(ix < 0xffffffff) ix += (ma[k].rid);\
			kx = (mk < 0?0xffffffff:(GMC(mmt, i, mk,m)&0xffffffff));\
			ks = MIN(ix, kx);\
			/**fprintf(stderr, "ks: %lu, i: %d (n-%d), k: %d (m-%d), ix: %lu, kx: %lu, t: %d, mk: %d\n", ks, i, n, k, m, ix, kx, t, mk);**/\
			if((ks&0xffffffff) == 0xffffffff) ks |= ((uint64_t)0xffffffff)<<32;\
			else if(ks == ix) ks |= (uint64_t)(i>0?(i-1)*m+t:0xffffffff)<<32;\
			else if(ks == kx) ks |= (uint64_t)(mk>=0?i*m+mk:0xffffffff)<<32;\
			GMC(mmt, i,k,m) = ks;\
			mk = k;\
		}\
	}\
	/**fprintf(stderr, "[M::%s::] ==> ++n: %d, m: %d, sn: %d, sidx: %d, eidx: %d\n", __func__, n, m, sn, sidx, eidx);**/\
	ks = (n-1)*m + mk; ix = (uint64_t)-1; kx = 0;\
	while (ks != 0xffffffff)\
	{\
		i = ks/m; k = ks%m;\
		ks = mmt[ks]>>32;\
		/**fprintf(stderr, "i: %d, k: %d, ks: %lu\n", i, k, ks);**/\
		if(ks == 0xffffffff || (int32_t)(ks/m) == (i-1)){\
			mm->a[sidx+k] = 1;\
			ix = MIN((uint64_t)k, ix); kx = MAX((uint64_t)k, kx);\
		}\
	}\
	/**debug_refine(ma, mmt, sn, n, m, (n-1)*m + mk);**/\
	if(rsi) (*rsi) = ix + sidx;\
	if(rei) (*rei) = kx + sidx;\
}\
void sf##_refine_sketch(VType *p, ha_pt_t *pt, int32_t rlen, int32_t dp_min_len, float er, int32_t min_freq, st_mt_t *mt)\
{\
	/**fprintf(stderr, "[M::%s::] ==> #########10#########, rlen: %d\n", __func__, rlen);**/\
	int32_t i, n = p->n, bd, len = MIN(rlen, dp_min_len), sublen, cnt, ei, li, ri;\
	int32_t sn =  len*er + 1;\
	kv_resize(uint64_t, *mt, (int64_t)p->n);\
	mt->n = p->n; memset(mt->a, 0, sizeof(uint64_t)*p->n);\
	for (i = 0; i < n; i++) p->a[i].rid = ha_pt_cnt(pt, p->a[i].x);\
	for (i = cnt = 0, bd = -1, ei = -1; i < n; i++){\
		if((int32_t)(p->a[i].rid)<min_freq) continue;\
		sublen = p->a[i].pos + 1;\
		if(sublen > len) break;\
		else ei = i;\
		if((int32_t)(p->a[i].pos + 1 - p->a[i].span) > bd){\
			bd = p->a[i].pos;\
			cnt++;\
		}\
	}\
	/**fprintf(stderr, "[M::%s::] ==> +cnt: %d, sn: %d, ei: %d, n: %d\n", __func__, cnt, sn, ei, n);**/\
	if(cnt >= sn) sf##_refine_select(p, 0, ei, sn, min_freq, mt, NULL, &li);\
	else{\
		li = i-1;\
		for (i = 0; i <= li; i++) mt->a[i] = 1;\
	}\
	if(len < rlen){\
		for (i = n-1, cnt = 0, bd = rlen+1, ei = -1; i >= 0; i--){\
			if((int32_t)(p->a[i].rid)<min_freq) continue;\
			sublen = rlen - (p->a[i].pos + 1 - p->a[i].span);\
			if(sublen > len) break;\
			else ei = i;\
			if((int32_t)(p->a[i].pos) < bd){\
				bd = p->a[i].pos + 1 - p->a[i].span;\
				cnt++;\
			}\
		}\
		/**fprintf(stderr, "[M::%s::] ==> -cnt: %d, sn: %d, ei: %d, n: %d\n", __func__, cnt, sn, ei, n);**/\
		if(cnt >= sn) sf##_refine_select(p, ei, n-1, sn, min_freq, mt, &ri, NULL);\
		else {\
			ri = i+1;\
			for (i = ri; i <= n-1; i++) mt->a[i] = 1;\
		}\
		/**fprintf(stderr, "[M::%s::] ==> --cnt: %d, sn: %d, ei: %d, n: %d\n", __func__, cnt, sn, ei, n);**/\
		if(ri - li >= 2){\
			li++; ri--;\
			sn =  (p->a[ri].pos - p->a[li].pos + p->a[li].span)*er + 1;\
			for (i = li, cnt = 0, bd = -1; i <= ri; i++){\
				if((int32_t)(p->a[i].rid)<min_freq) continue;\
				if((int32_t)(p->a[i].pos + 1 - p->a[i].span) > bd){\
					bd = p->a[i].pos;\
					cnt++;\
					if(cnt >= sn) break;\
				}\
			}\
			if(cnt >= sn) sf##_refine_select(p, li, ri, sn, min_freq, mt, NULL, NULL);\
			else for (i = li; i <= ri; i++) mt->a[i] = 1;\
		}\
	}\
	/**fprintf(stderr, "[M::%s::] ==> #########20#########, p->n: %u, n: %d\n", __func__, p->n, n);**/\
	for (i = sn = 0; i < n; i++){\
		if(mt->a[i]){\
			p->a[sn] = p->a[i];\
			sn++;\
		}\
	}\
	/**if(p->n != sn) fprintf(stderr, "[M::%s::] ==> #########21#########, p->n: %u, sn: %d\n", __func__, p->n, sn);**/\
	p->n = sn;\
}\
/**\
 * Find symmetric (w,k)-minimizers on a DNA sequence\
 *\
 * @param str    DNA sequence\
 * @param len    length of $str\
 * @param w      find a minimizer for every $w consecutive k-mers\
 * @param k      k-mer size\
 * @param rid    reference ID; will be copied to the output $p array\
 * @param is_hpc homopolymer-compressed or not\
 * @param p      minimizers\
 */\
void sf##_ha_sketch(const char *str, int len, int w, int k, uint32_t rid, int is_hpc, VType *p, const void *hf, int sample_dist, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, ha_pt_t *pt, int min_freq, int32_t dp_min_len, float dp_e, st_mt_t *mt, int32_t ws, int32_t is_unique)\
{   /**in default, w = 51, k = 51, is_hpc = 1**/\
    extern void *ha_ct_table;\
    static const HType dummy = { UINT64_MAX, (((uint64_t)1)<<RidBits) - 1, 0, 0, 0};\
    uint64_t shift1 = k - 1, mask = (1ULL<<k) - 1, kmer[4] = {0,0,0,0};\
    int i, j, l, tl = 0, buf_pos, min_pos, kmer_span = 0;\
    HType buf[256], min = dummy;\
    uint32_t buf_p[256], min_s = (uint32_t)-1;\
    tiny_queue_t tq;\
    assert(len > 0 && (int64_t)(len) < (int64_t)((((uint64_t)1)<<PosBits)) && (int64_t)(rid) < (int64_t)((((uint64_t)1)<<RidBits)) && (w > 0 && w < 256) && (k > 0 && k <= 63));\
    if (dbg_ct != NULL) dbg_ct->a.n = 0;\
    if (k_flag != NULL) {\
        kv_resize(uint8_t, k_flag->a, (uint64_t)len);\
        k_flag->a.n = len;\
        memset(k_flag->a.a, 0, k_flag->a.n);\
    }\
    memset(buf, 0xff, w * sizeof(HType));\
    memset(&tq, 0, sizeof(tiny_queue_t));\
    /**len/w is the evaluated minimizer numbers**/\
    kv_resize(HType, *p, p->n + len/w);\
    kv_resize(uint64_t, *mt, (int64_t)p->m); mt->n = p->n;\
    for (i = l = tl = buf_pos = min_pos = 0; i < len; ++i) {\
        int c = seq_nt4_table[(uint8_t)str[i]];\
        HType info = dummy;\
        if (c < 4) { /**not an ambiguous base**/\
            int z;\
            if (is_hpc) {\
                int skip_len = 1;\
                if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {\
                    for (skip_len = 2; i + skip_len < len; ++skip_len)\
                        if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)\
                            break;\
                    i += skip_len - 1; /**put $i at the end of the current homopolymer run**/\
                }\
                tq_push(&tq, skip_len);\
                kmer_span += skip_len;\
                /**how many bases that are covered by this HPC k-mer\
                kmer_span includes at most k HPC elements**/\
                if (tq.count > k) kmer_span -= tq_shift(&tq);\
            } else kmer_span = l + 1 < k? l + 1 : k;\
            /**kmer_span should be used for HPC k-mer\
            non-HPC k-mer, kmer_span should be k\
            kmer_span is used to calculate anchor pos on reverse complementary strand**/\
            if (k_flag != NULL) k_flag->a.a[i] = 1;/**lable all useful base, which are not ignored by HPC**/\
            kmer[0] = (kmer[0] << 1 | (c&1))  & mask;/**forward k-mer**/\
            kmer[1] = (kmer[1] << 1 | (c>>1)) & mask;\
            kmer[2] = kmer[2] >> 1 | (uint64_t)(1 - (c&1))  << shift1; /**reverse k-mer**/\
            kmer[3] = kmer[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift1;\
            if (kmer[1] == kmer[3]) continue; /** skip "symmetric k-mers" as we don't know it strand**/\
            z = kmer[1] < kmer[3]? 0 : 1; /** strand**/\
            ++l; tl++;\
            if (l >= k && kmer_span < 256) {\
                uint64_t y;\
                int32_t cnt, filtered;\
                y = yak_hash64_64(kmer[z<<1|0]) + yak_hash64_64(kmer[z<<1|1]);\
                cnt = hf? ha_ft_cnt(hf, y) : 0;\
				filtered = (cnt >= 1<<28);\
                if(is_unique && (!filtered)) {\
					filtered = (cnt == 0);\
					cnt = (cnt == 1? 0:cnt);\
				}\
                if (dbg_ct != NULL) kv_push(uint64_t, dbg_ct->a, ((((uint64_t)(query_ct_index(ha_ct_table, y))<<1)|filtered)<<32)|(uint64_t)(i));\
                if (!filtered) info.x = y, info.rid = cnt, info.pos = i, info.rev = z, info.span = kmer_span; /** initially ha_mz1_t::rid keeps the k-mer count**/\
                if (k_flag != NULL) k_flag->a.a[i]++;\
                if (k_flag != NULL && filtered > 0) k_flag->a.a[i]++;\
            }\
        } else l = 0, tq.count = tq.front = 0, kmer_span = 0;\
        buf[buf_pos] = info; /**need to do this here as appropriate buf_pos and buf[buf_pos] are needed below**/\
        buf_p[buf_pos] = l;\
        if (l == w + k - 1 && min.x != UINT64_MAX) { /**special case for the first window - because identical k-mers are not stored yet**/\
			for (j = buf_pos + 1; j < w; ++j){\
                if (sf##_mzcmp(&min, &buf[j]) == 0 && buf[j].pos != min.pos){\
                    kv_push(HType, *p, buf[j]); kv_push(uint64_t, *mt, buf_p[j]);\
                }\
            }\
            for (j = 0; j < buf_pos; ++j){\
                if (sf##_mzcmp(&min, &buf[j]) == 0 && buf[j].pos != min.pos){\
                    kv_push(HType, *p, buf[j]); kv_push(uint64_t, *mt, buf_p[j]);\
                }\
            }\
        }\
        /**\
         * There are three cases:\
         * 1. info.x <= min.x, means info is a new minimizer\
         * 2. info.x > min.x, info is not a new minimizer\
         *    (1) buf_pos != min_pos, do nothing\
         *    (2) buf_pos == min_pos, means current minimizer has moved outside the window\
         * **/\
        /**three cases: 1.**/\
        if (sf##_mzcmp(&min, &info) >= 0) { /**a new minimum; then write the old min**/\
            if (l >= w + k && min.x != UINT64_MAX){\
                kv_push(HType, *p, min); kv_push(uint64_t, *mt, min_s);\
            }\
            min = info, min_pos = buf_pos, min_s = buf_p[buf_pos];\
        } else if (buf_pos == min_pos) { /**old min has moved outside the window**/\
            if (l >= w + k - 1 && min.x != UINT64_MAX){\
                kv_push(HType, *p, min); kv_push(uint64_t, *mt, min_s);\
            }\
            /**buf_pos == min_pos, means current minimizer has moved outside the window\
            so for now we need to find a new minimizer at the current window (w k-mers)**/\
            for (j = buf_pos + 1, min = dummy; j < w; ++j) /**the two loops are necessary when there are identical k-mers**/\
                if (sf##_mzcmp(&min, &buf[j]) >= 0) min = buf[j], min_pos = j, min_s = buf_p[j]; /** >= is important s.t. min is always the closest k-mer**/\
            for (j = 0; j <= buf_pos; ++j)\
                if (sf##_mzcmp(&min, &buf[j]) >= 0) min = buf[j], min_pos = j, min_s = buf_p[j];\
            if (l >= w + k - 1 && min.x != UINT64_MAX) { /**write identical k-mers**/\
                for (j = buf_pos + 1; j < w; ++j) /**these two loops make sure the output is sorted**/\
                    if (sf##_mzcmp(&min, &buf[j]) == 0 && min.pos != buf[j].pos){\
                        kv_push(HType, *p, buf[j]); kv_push(uint64_t, *mt, buf_p[j]);\
                    }\
                for (j = 0; j <= buf_pos; ++j)\
                    if (sf##_mzcmp(&min, &buf[j]) == 0 && min.pos != buf[j].pos){\
                        kv_push(HType, *p, buf[j]); kv_push(uint64_t, *mt, buf_p[j]);\
                    }\
            }\
        }\
        if (++buf_pos == w) buf_pos = 0;\
    }\
    if (min.x != UINT64_MAX){\
        kv_push(HType, *p, min); kv_push(uint64_t, *mt, min_s);\
    }\
	/**debug_pl(str, len, w, k, is_hpc, p, hf, mt);**/\
    if (sample_dist > w) sf##_select_mz_h(p, mt, len, sample_dist, ws, k, tl);\
    if (dp_min_len > 0 && pt && mt) sf##_refine_sketch(p, pt, len, dp_min_len, dp_e, min_freq, mt);\
	for (i = 0; i < (int)p->n; ++i) /**populate .rid as this was keeping counts**/\
		p->a[i].rid = rid;\
}

HA_SC_INIT(mz1, ha_mz1_t, ha_mz1_v, 28, 27)
HA_SC_INIT(mz2, ha_mzl_t, ha_mzl_v, 31, 32)