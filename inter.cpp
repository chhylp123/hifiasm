#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <zlib.h>
#include "kseq.h" // FASTA/Q parser
#include "kavl.h"
#include "khash.h"
#include "kalloc.h"
#include "kthread.h"
#include "inter.h"
#include "Overlaps.h"
#include "CommandLines.h"
#include "htab.h"
#include "Hash_Table.h"
KSEQ_INIT(gzFile, gzread)

#define MG_SEED_IGNORE     (1ULL<<41)
#define MG_SEED_TANDEM     (1ULL<<42)
#define MG_SEED_KEPT       (1ULL<<43)

#define MG_MAX_SEG        255
#define MG_SEED_SEG_SHIFT  48
#define MG_SEED_SEG_MASK   (0xffULL<<(MG_SEED_SEG_SHIFT))
#define mg_seg_id(a) ((int32_t)(((a).y&MG_SEED_SEG_MASK) >> MG_SEED_SEG_SHIFT))

#define MG_SEED_WT_SHIFT   56
#define MG_MAX_SHORT_K  15

#define MG_SHORT_K_EXT 10000 ///1000 in minigraph


#define generic_key(x) (x)
KRADIX_SORT_INIT(gfa64, uint64_t, generic_key, 8)

struct mg_tbuf_s {
	void *km;
	int frag_gap;
};
typedef struct mg_tbuf_s mg_tbuf_t;


mg_tbuf_t *mg_tbuf_init(void)
{
	mg_tbuf_t *b;
	b = (mg_tbuf_t*)calloc(1, sizeof(mg_tbuf_t));
	b->km = km_init();
	return b;
}

void mg_tbuf_destroy(mg_tbuf_t *b)
{
	if (b == 0) return;
	if (b->km) km_destroy(b->km);
	free(b);
}

void *mg_tbuf_get_km(mg_tbuf_t *b)
{
	return b->km;
}

typedef struct {
	int w, k, bw, max_gap, is_HPC, hap_n, occ_weight, max_gap_pre;
	int max_lc_skip, max_lc_iter, min_lc_cnt, min_lc_score;
    float chn_pen_gap;
} mg_idxopt_t;

typedef struct {
    ha_abufl_t *abl;
    st_mt_t sp;
    Candidates_list clist;
    overlap_region_alloc olist;
} ma_ov_buf_t;

typedef struct { // global data structure for kt_pipeline()
    const void *ha_flt_tab;
    const ha_pt_t *ha_idx;
    const mg_idxopt_t *opt;
    const ma_ug_t *ug;
	kseq_t *ks;
    int64_t chunk_size;
    uint64_t n_thread;
    uint64_t total_base;
    uint64_t total_pair;
} uldat_t;

typedef struct { // data structure for each step in kt_pipeline()
    const mg_idxopt_t *opt;
    const void *ha_flt_tab;
    const ha_pt_t *ha_idx;
    const ma_ug_t *ug;
	int n, m, sum_len;
	uint64_t *len, id;
	char **seq;
    ha_mzl_v *mzs;
    st_mt_t *sps;
    mg_tbuf_t **buf;
} utepdat_t;

typedef struct {
	///off: start idx in mg128_t * a[];
	///cnt: how many eles in this chain
	///a[off, off+cnt) saves the eles in this chain
	int32_t off, cnt:31, inner_pre:1;
	///ref_id|rev
	uint32_t v;
	///chain in ref: [rs, re)
	///chain in query: [qs, qe)
	int32_t rs, re, qs, qe;
	///score: chain score
	int32_t score, dist_pre;
	uint32_t hash_pre;
} mg_lchain_t;

typedef struct {
	uint32_t v, d;
	int32_t pre;
} mg_pathv_t;

KHASH_MAP_INIT_INT(sp, sp_topk_t)
KHASH_MAP_INIT_INT(sp2, uint64_t)

///mg128_t->y: weight(8)seg_id(8)flag(8)span(8)pos(32)
///mg128_t->x: rid(31)rev(1)pos(33); keep reference
typedef struct { uint64_t x, y; } mg128_t;
#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, mg128_t, sort_key_128x, 8) 
void radix_sort_128x(mg128_t *beg, mg128_t *end);

typedef struct {
	uint32_t n; ///length of candidate list 
    uint64_t q_span:31, rev:1, q_pos:32;
	uint32_t qid:16, weight:15, is_tandem:1;
	const ha_idxposl_t *cr; ///candidate list
} mg_match_t;

// shortest path
typedef struct {
	// input
	///(lj_ref_id)|(lj_ref_rev^1)
	uint32_t v;
	///target_dist should like the overlap length in string graph
	///it should be used to evaluate if the identified path is close to real path/alignment
	int32_t target_dist;
	uint32_t target_hash;
	///inner: if li and lj are at the same ref id
	///meta: j
	uint32_t meta:30, check_hash:1, inner:1;
	/**
	 * There are two cases:
	 * (1) lj->qs************lj->qe
	 * 				li->qs************li->qe
	 * (2) lj->qs************lj->qe
	 * 									li->qs************li->qe
	 * qlen = li->qs - lj->qe;///might be negative
	 * **/
	int32_t qlen;
	// output
	uint32_t n_path:31, is_0:1;///I guess n_path is how many path from src to dest
	int32_t path_end;///looks like an idx to alignment
	int32_t dist, mlen;
	uint32_t hash;
	// aux
	uint64_t srt_key;
} mg_path_dst_t;

typedef struct {
	uint32_t srt;
	int32_t i;
} gc_frag_t;

///I think this structure is just used for iteration
///iterate each ref id, instead of alignment id
typedef struct sp_node_s {
	uint64_t di; // dist<<32 | node_id in avl tree(doesn't matter too much)
	uint32_t v;///ref_id|rev
	int32_t pre;
	uint32_t hash;
	int32_t is_0;
	KAVL_HEAD(struct sp_node_s) head;
} sp_node_t, *sp_node_p;

typedef struct {
	int32_t k, mlen;//k: number of walks from src to this node
	int32_t qs, qe;
	sp_node_t *p[MG_MAX_SHORT_K]; // this forms a max-heap
} sp_topk_t;

#define gc_frag_key(p) ((p).srt)
KRADIX_SORT_INIT(gc, gc_frag_t, gc_frag_key, 4)

#define dst_key(p) ((p).srt_key)
KRADIX_SORT_INIT(dst, mg_path_dst_t, dst_key, 8)

#define sp_node_cmp(a, b) (((a)->di > (b)->di) - ((a)->di < (b)->di))
KAVL_INIT(sp, sp_node_t, head, sp_node_cmp)

#define sp_node_lt(a, b) ((a)->di < (b)->di)
KSORT_INIT(sp, sp_node_p, sp_node_lt)

void init_mg_opt(mg_idxopt_t *opt, int is_HPC, int k, int w, int hap_n)
{
    opt->k = k; 
    opt->w = w;
    opt->hap_n = hap_n;
    opt->is_HPC = is_HPC;
    opt->bw = 2000;
    opt->max_gap = 5000;
    opt->occ_weight = 20;
    opt->max_gap_pre = 10000;///1000 in minigraph
    opt->max_lc_iter = 10000;
    opt->chn_pen_gap = 0.19;///using minimap2's value
	opt->max_lc_skip = 25;// mo->max_gc_skip = 25;
	opt->max_lc_iter = 10000;
	opt->min_lc_cnt = 2; 
	opt->min_lc_score = 30;
}

void uidx_build(ma_ug_t *ug, mg_idxopt_t *opt)
{
    int flag = asm_opt.flag;
    asm_opt.flag |= HA_F_NO_HPC;
    ha_flt_tab = ha_ft_ug_gen(&asm_opt, &(ug->u), opt->is_HPC, opt->k, opt->w, 1, opt->hap_n*10);
    ha_idx = ha_pt_ug_gen(&asm_opt, ha_flt_tab, &(ug->u), opt->is_HPC, opt->k, opt->w, 1);
    asm_opt.flag = flag;
}

void uidx_destory()
{
    ha_ft_destroy(ha_flt_tab); 
    ha_pt_destroy(ha_idx);
}


///only use non-repetitive minimizers
static mg_match_t *collect_matches(void *km, int *_n_m, int max_occ, const void *ha_flt_tab, const ha_pt_t *ha_idx, int check_unique, const ha_mzl_v *mv, int64_t *n_a, int *rep_len, int *n_mini_pos, int32_t **mini_pos)
{
	int rep_st = 0, rep_en = 0, n_m, tn, tw;
	size_t i;
	mg_match_t *m;
	*n_mini_pos = 0;
    KMALLOC(km, *mini_pos, mv->n);///mv->n how many minimizers in query
	m = (mg_match_t*)kmalloc(km, mv->n * sizeof(mg_match_t));
	for (i = 0, n_m = 0, *rep_len = 0, *n_a = 0; i < mv->n; ++i) {
		const ha_idxposl_t *cr;
		ha_mzl_t *z = &mv->a[i];
		cr = ha_ptl_get(ha_idx, z->x, &tn);        
        tw = ha_ft_cnt(ha_flt_tab, z->x);
		if (tw > max_occ) { ///the frequency of repetitive regions; ignore those minimizers
			int en = z->pos + 1, st = en - z->span;//[st, en)
			if (st > rep_en) { ///just record the length of repetive regions
				*rep_len += rep_en - rep_st;
				rep_st = st, rep_en = en;
			} else rep_en = en;
		} else {
			mg_match_t *q = &m[n_m++];
			q->q_pos = z->pos, q->q_span = z->span, q->rev = z->rev, q->cr = cr, q->n = tn, q->qid = 0;
			q->is_tandem = 0, q->weight = 255;
            if(check_unique && tw != 1) q->is_tandem = 1, q->weight = 15;
			*n_a += q->n;///how many candidates 
			(*mini_pos)[(*n_mini_pos)++] = z->pos;///minimizer offset in query
		}
	}
	*rep_len += rep_en - rep_st; ///the length of repetitive regions 
	*_n_m = n_m;
	return m;
}

mg128_t *collect_seed_hits(void *km, const mg_idxopt_t *opt, int max_occ, const void *ha_flt_tab, const ha_pt_t *ha_idx, 
const ma_ug_t *ug, const ha_mzl_v *mv, int64_t *n_a, int *rep_len, int *n_mini_pos, int32_t **mini_pos)
{
    int i, n_m;
    mg128_t *a = NULL;
    mg_match_t *m = collect_matches(km, &n_m, max_occ, ha_flt_tab, ha_idx, 1, mv, n_a, rep_len, n_mini_pos, mini_pos);
    a = (mg128_t*)kmalloc(km, *n_a * sizeof(mg128_t));///n_a: how many available candidates in total 
    for (i = 0, *n_a = 0; i < n_m; ++i) {///n_m: how many available seeds, instead of candidates
		mg_match_t *q = &m[i];
		const ha_idxposl_t *r = q->cr;
		uint32_t k;
		for (k = 0; k < q->n; ++k) {///q->n: number of candidates belonging to seed m[i]
			mg128_t *p;
			p = &a[(*n_a)++];///pick up a slot for one candidate
			if (r[k].rev == q->rev) // forward strand
				p->x = (uint64_t)(r[k].rid)<<33|r[k].pos; ///reference: rid(31)|rev(1)|pos(32)                
			else // reverse strand
				p->x = (uint64_t)(r[k].rid)<<33 | 1ULL<<32 | (ug->g->seq[r[k].rid].len - (r[k].pos + 1 - r[k].span) - 1);
			p->y = (uint64_t)q->q_span << 32 | q->q_pos;
			p->y |= (uint64_t)q->qid << MG_SEED_SEG_SHIFT;
			if (q->is_tandem) p->y |= MG_SEED_TANDEM;
			p->y |= (uint64_t)q->weight << MG_SEED_WT_SHIFT; 
			///p->y: weight(8)seg_id(8)flag(8)span(8)pos(32)
			///p->x: rid(31)rev(1)pos(33); keep reference
		}
	}
    kfree(km, m);
	radix_sort_128x(a, a + (*n_a));
	return a;
}

///r is 1000 in default
///remove isolated hits, whic are not close enough to others
static int64_t flt_anchors(int64_t n_a, mg128_t *a, int32_t r)
{
	int64_t i, j;
	for (i = 0; i < n_a; ++i) {
		for (j = i - 1; j >= 0; --j) {
			/**
			 * a is sorted by x
			 * a[].x: ref_id(31)rev(1)r_pos(32)
			 * a[].y: weight(8)query_id(8)flag(8)span(8)q_pos(32)
			 **/
			int32_t dq;
			int64_t dr = a[i].x - a[j].x;///a is sorted by x
			if (dr > r) break;///if two candidates coming from differnt unitigs, dr would be extremly large
			dq = (int32_t)a[i].y - (int32_t)a[j].y;
			if (dq > r || dq < 0) continue;
			a[j].y |= MG_SEED_KEPT;
			a[i].y |= MG_SEED_KEPT;
			break;
		}
	}
	for (i = n_a - 1; i >= 0; --i) {
		if (a[i].y & MG_SEED_KEPT) continue;
		for (j = i + 1; j < n_a; ++j) {
			int32_t dq;
			int64_t dr = a[j].x - a[i].x;
			if (dr > r) break;
			dq = (int32_t)a[j].y - (int32_t)a[i].y;
			if (dq > r || dq < 0) continue;
			a[j].y |= MG_SEED_KEPT;
			a[i].y |= MG_SEED_KEPT;
			break;
		}
	}
	for (i = j = 0; i < n_a; ++i)
		if (a[i].y & MG_SEED_KEPT)
			a[j++] = a[i];
	return j;
}

static inline float mg_log2(float x) // NB: this doesn't work when x<2
{
	union { float f; uint32_t i; } z = { x };
	float log_2 = ((z.i >> 23) & 255) - 128;
	z.i &= ~(255 << 23);
	z.i += 127 << 23;
	log_2 += (-0.34484843f * z.f + 2.02466578f) * z.f - 0.67487759f;
	return log_2;
}

inline int32_t normal_sc(uint64_t w, int32_t sc)
{
    if(w < 255){
        int32_t tmp = (int)(0.00392156862745098 * w * sc); // 0.00392... = 1/255
        sc = tmp > 1? tmp : 1;
    }
    return sc;
}
// ai[].x: ref_id(31)rev(1)r_pos(32)
// ai[].y: weight(8)query_id(8)flag(8)span(8)q_pos(32)
// comput_sc(&a[i], &a[j], max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, is_cdna, n_segs);
static inline int32_t comput_sc(const mg128_t *ai, const mg128_t *aj, int32_t max_dist_x, int32_t max_dist_y, int32_t bw, float chn_pen_gap)
{	
	int32_t dq = (int32_t)ai->y - (int32_t)aj->y, dr = (int32_t)ai->x - (int32_t)aj->x, dd, dg, q_span, sc;
	///ai and aj has already been sorted by x
	///which means ai->x >= aj->x
	if (dq <= 0 || dq > max_dist_x) return INT32_MIN;
    if (dr <= 0 || dr > max_dist_y) return INT32_MIN;
	dd = dr > dq? dr - dq : dq - dr; ///indel, dd is always >= 0
	if (dd > bw) return INT32_MIN;
	dg = dr < dq? dr : dq;///MIN(dr, dq)
	q_span = aj->y>>32&0xff;///query span; should be ai->y>>32&0xff, is it a bug?
    sc = normal_sc(aj->y>>MG_SEED_WT_SHIFT, (q_span<dg?q_span:dg)); ///positive part of sc
    ///dd: there are indels
    ///dg > q_span: there are some bases that are not covered between ai and aj
    ///it is if (dd || dg > q_span) in minigraph
	if (dd) {
		float lin_pen, log_pen;
		lin_pen = chn_pen_gap * (float)dd;
		log_pen = dd >= 2? mg_log2(dd) : 0.0f; // mg_log2() only works for dd>=2
        sc -= (int)(lin_pen + log_pen);
	}
	return sc;
}


///p[]: id of last
///f[]: the score ending at i, not always the peak
///v[]: keeps the peak score up to i;
///t[]: id of next
///min_cnt = 2; min_sc = 30; extra_u = 0
///u = mg_chain_backtrack(n, f, p, v, t, min_cnt, min_sc, 0, &n_u, &n_v);
uint64_t *mg_chain_backtrack(void *km, int64_t n, const int32_t *f, const int64_t *p, int32_t *v, int32_t *t, int32_t min_cnt, int32_t min_sc, int32_t extra_u, int32_t *n_u_, int32_t *n_v_)
{
	mg128_t *z;
	uint64_t *u;
	int64_t i, k, n_z, n_v;
	int32_t n_u;
	// v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
	*n_u_ = *n_v_ = 0;
	for (i = 0, n_z = 0; i < n; ++i) // precompute n_z
		if (f[i] >= min_sc) ++n_z;
	if (n_z == 0) return 0;
	KMALLOC(km, z, n_z);
	for (i = 0, k = 0; i < n; ++i) // populate z[]
		if (f[i] >= min_sc) z[k].x = f[i], z[k++].y = i;
	radix_sort_128x(z, z + n_z);///sort by score

	memset(t, 0, n * 4);///t is a buffer
	///from the largest to the smallest
	for (k = n_z - 1, n_v = n_u = 0; k >= 0; --k) { // precompute n_u
		int64_t n_v0 = n_v;
		int32_t sc;
		///note t[i] == 0 is not used to find local alignment
		///say if we have already found a long chain, then the secondary might be able to merged to the long chain
		///t[i] == 0 is used to find those chains
		for (i = z[k].y; i >= 0 && t[i] == 0; i = p[i])
			++n_v, t[i] = 1;
		sc = i < 0? z[k].x : (int32_t)z[k].x - f[i];
		if (sc >= min_sc && n_v > n_v0 && n_v - n_v0 >= min_cnt)
			++n_u;///how many chains, including primary chains and non-primary chains
		else n_v = n_v0;
	}
	KMALLOC(km, u, n_u + extra_u);
	memset(t, 0, n * 4);
	for (k = n_z - 1, n_v = n_u = 0; k >= 0; --k) { // populate u[]
		int64_t n_v0 = n_v;
		int32_t sc;
		for (i = z[k].y; i >= 0 && t[i] == 0; i = p[i])
			v[n_v++] = i, t[i] = 1;
		sc = i < 0? z[k].x : (int32_t)z[k].x - f[i];
		if (sc >= min_sc && n_v > n_v0 && n_v - n_v0 >= min_cnt)
			u[n_u++] = (uint64_t)sc << 32 | (n_v - n_v0);
		else n_v = n_v0;
	}
	kfree(km, z);
	assert(n_v < INT32_MAX);
	*n_u_ = n_u, *n_v_ = n_v;
	return u;
}

//u[]: sc|occ of chains
//v[]: idx of each element
static mg128_t *compact_a(void *km, int32_t n_u, uint64_t *u, int32_t n_v, int32_t *v, mg128_t *a)
{
	mg128_t *b, *w;
	uint64_t *u2;
	int64_t i, j, k;

	// write the result to b[]
	KMALLOC(km, b, n_v);
	for (i = 0, k = 0; i < n_u; ++i) {
		int32_t k0 = k, ni = (int32_t)u[i];
		for (j = 0; j < ni; ++j)
			b[k++] = a[v[k0 + (ni - j - 1)]];///write all elements of a chain together
	}
	kfree(km, v);

	// sort u[] and a[] by the target position, such that adjacent chains may be joined
	KMALLOC(km, w, n_u);
	for (i = k = 0; i < n_u; ++i) {///n_u: how many chains
		///x: ref_id(31)rev(1)r_pos(32)
		w[i].x = b[k].x, w[i].y = (uint64_t)k<<32|i;
		k += (int32_t)u[i];
	}
	radix_sort_128x(w, w + n_u);///sort by ref_id(31)rev(1)r_pos(32); r_pos is the start pos of chain
	KMALLOC(km, u2, n_u);
	for (i = k = 0; i < n_u; ++i) {///note merge chain; just place close chains together
		///j is chain id; n is how many elements in j-th chain
		int32_t j = (int32_t)w[i].y, n = (int32_t)u[j];
		u2[i] = u[j];
		memcpy(&a[k], &b[w[i].y>>32], n * sizeof(mg128_t));
		k += n;
	}
	memcpy(u, u2, n_u * 8);
	memcpy(b, a, k * sizeof(mg128_t)); // write _a_ to _b_ and deallocate _a_ because _a_ is oversized, sometimes a lot
	kfree(km, a); kfree(km, w); kfree(km, u2);
	return b;
}

/* Input:
 *   a[].x: ref_id(31)rev(1)r_pos(32)
 *   a[].y: weight(8)query_id(8)flag(8)span(8)q_pos(32)
 *   n: length of a[]
 * Output:
 *   n_u: #chains
 *   u[]: score<<32 | #anchors (sum of lower 32 bits of u[] is the returned length of a[])
 * input a[] is deallocated on return
 */
///is_cdna is is_splice
mg128_t *mg_lchain_dp(int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, int min_cnt, int min_sc, float chn_pen_gap, 
					  int64_t n, mg128_t *a, int *n_u_, uint64_t **_u, void *km)
{ // TODO: make sure this works when n has more than 32 bits
	int32_t *f, *t, *v, n_u, n_v;
	int64_t *p, i, j, max_ii, st = 0;
	uint64_t *u;

	if (_u) *_u = 0, *n_u_ = 0;
	if (n == 0 || a == 0) return 0;
	KMALLOC(km, p, n);///id of last cell
	KMALLOC(km, f, n);///f[] is the score ending at i, not always the peak
	KMALLOC(km, v, n);///v[] keeps the peak score up to i;
	KCALLOC(km, t, n);///t doesn't matter too much; it is mainly used to accelrate the iteration

	// a[].x: ref_id(31)rev(1)r_pos(32)
 	// a[].y: weight(8)query_id(8)flag(8)span(8)q_pos(32)
	// fill the score and backtrack arrays
	for (i = st = 0, max_ii = -1; i < n; ++i) {
		int64_t max_j = -1, end_j;
        ///max_f -> minimizer span length in query, which is the initial score
		int32_t max_f = normal_sc(a[i].y>>MG_SEED_WT_SHIFT, a[i].y>>32&0xff), n_skip = 0;
		///until we are at the same rid, same direction, and the coordinates are close enough
		while (st < i && (a[i].x>>32 != a[st].x>>32 || a[i].x > a[st].x + max_dist_x)) ++st;
		///max_iter = 10000 in default, which means dp can go back to up to 10000 cells
		if (i - st > max_iter) st = i - max_iter; 
		for (j = i - 1; j >= st; --j) {
			int32_t sc;
			sc = comput_sc(&a[i], &a[j], max_dist_x, max_dist_y, bw, chn_pen_gap);
			if (sc == INT32_MIN) continue;
			sc += f[j];
			if (sc > max_f) {
				max_f = sc, max_j = j;
				if (n_skip > 0) --n_skip;
			} else if (t[j] == (int32_t)i) {///note we scan j backwards; we don't need to update t[] for each i
				if (++n_skip > max_skip)
					break;
			}
			if (p[j] >= 0) t[p[j]] = i;
		}
		end_j = j;///end_j might be > 0
		///if not close enough, select a new max
		if (max_ii < 0 || (int64_t)(a[i].x - a[max_ii].x) > (int64_t)max_dist_x) {///select a new max
			int32_t max = INT32_MIN;
			max_ii = -1;
			for (j = i - 1; j >= st; --j)
				if (max < f[j]) max = f[j], max_ii = j;
		}
		///note: it will happen when `max_ii` < `end_j`; 
		///iteration is terminated at `end_j` mostly because of `max_skip` and `max_iter`
		if (max_ii >= 0 && max_ii < end_j) {
			int32_t tmp;
			tmp = comput_sc(&a[i], &a[max_ii], max_dist_x, max_dist_y, bw, chn_pen_gap);
			if (tmp != INT32_MIN && max_f < tmp + f[max_ii])
				max_f = tmp + f[max_ii], max_j = max_ii;
		}
		// v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
		f[i] = max_f, p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; 
		if (max_ii < 0 || ((int64_t)(a[i].x - a[max_ii].x) <= (int64_t)max_dist_x && f[max_ii] < f[i]))
			max_ii = i;
	}
	///after mg_chain_backtrack, the results are saved in u and v;
	u = mg_chain_backtrack(km, n, f, p, v, t, min_cnt, min_sc, 0, &n_u, &n_v);
	*n_u_ = n_u, *_u = u; // NB: note that u[] may not be sorted by score here
	kfree(km, p); kfree(km, f); kfree(km, t);
	if (n_u == 0) {
		kfree(km, a); kfree(km, v);
		return 0;
	}
	//u[]: sc|occ of chains
	//v[]: idx of each element
	return compact_a(km, n_u, u, n_v, v, a);
}

///qlen: query length
///u[]: sc|occ of chains
///a[]: candidate list
mg_lchain_t *mg_lchain_gen(void *km, int qlen, int n_u, uint64_t *u, mg128_t *a)
{
	mg128_t *z;
	mg_lchain_t *r;
	int i, k;

	if (n_u == 0) return 0;
	KCALLOC(km, r, n_u); KMALLOC(km, z, n_u);
	// u[] is sorted by query position
	for (i = k = 0; i < n_u; ++i) {
		/**
		 * a[].x: ref_id(31)rev(1)r_pos(32)
		 * a[].y: weight(8)query_id(8)flag(8)span(8)q_pos(32)
		**/
		///u[]: sc(32)occ(32)
		int32_t qs = (int32_t)a[k].y + 1 - (a[k].y>>32 & 0xff);
		z[i].x = (uint64_t)qs << 32 | u[i] >> 32;
		z[i].y = (uint64_t)k << 32 | (int32_t)u[i];
		k += (int32_t)u[i];
	}
	radix_sort_128x(z, z + n_u);//sort by qs|sc

	// populate r[]
	for (i = 0; i < n_u; ++i) {
		mg_lchain_t *ri = &r[i];
		/**
		 * z[].x: query start pos| chain score
		 * z[].y: idx in a[] | chain occ
		 * a[].x: ref_id(31)rev(1)r_pos(32)
		 * a[].y: weight(8)query_id(8)flag(8)span(8)q_pos(32)
		 * **/
		int32_t k = z[i].y >> 32, q_span = a[k].y >> 32 & 0xff;
		ri->off = k;
		ri->cnt = (int32_t)z[i].y;
		ri->score = (uint32_t)z[i].x;
		ri->v = a[k].x >> 32;///ref_id|rev
		ri->rs = (int32_t)a[k].x + 1 > q_span? (int32_t)a[k].x + 1 - q_span : 0; // for HPC k-mer
		ri->qs = z[i].x >> 32;
		ri->re = (int32_t)a[k + ri->cnt - 1].x + 1;
		ri->qe = (int32_t)a[k + ri->cnt - 1].y + 1;
	}
	kfree(km, z);
	return r;
}

static int32_t get_mini_idx(const mg128_t *a, int32_t n, const int32_t *mini_pos)
{
	int32_t x, L = 0, R = n - 1;
	x = (int32_t)a->y;
	while (L <= R) { // binary search
		int32_t m = ((uint64_t)L + R) >> 1;
		int32_t y = mini_pos[m];
		if (y < x) L = m + 1;
		else if (y > x) R = m - 1;
		else return m;
	}
	return -1;
}

/* Before:
 *   a[].x: ref_id(31)rev(1)r_pos(32)
 *   a[].y: weight(8)query_id(8)flag(8)span(8)q_pos(32)
 * After:
 *   a[].x: idx_in_minimizer_arr(32)r_pos(32)
 *   a[].y: weight(8)query_id(8)flag(8)span(8)q_pos(32)
 */
void mg_update_anchors(int32_t n_a, mg128_t *a, int32_t n, const int32_t *mini_pos)
{
	int32_t st, j, k;
	if (n_a <= 0) return;
	st = get_mini_idx(&a[0], n, mini_pos);
	assert(st >= 0);
	for (k = 0, j = st; j < n && k < n_a; ++j)
		if ((int32_t)a[k].y == mini_pos[j])
			a[k].x = (uint64_t)j << 32 | (a[k].x & 0xffffffffU), ++k;
	assert(k == n_a);
}

static int32_t find_max(int32_t n, const gc_frag_t *gf, uint32_t x)
{
	int32_t s = 0, e = n;
	if (n == 0) return -1;
	if (gf[n-1].srt < x) return n - 1;
	if (gf[0].srt >= x) return -1;
	while (e > s) { // TODO: finish this block
		int32_t m = s + (e - s) / 2;
		if (gf[m].srt >= x) e = m;
		else s = m + 1;
	}
	assert(s == e);
	return s;
}
///target_dist should like the overlap length in string graph
///it should be used to evaluate if the identified path is close to real path/alignment
static int32_t mg_target_dist(const asg_t *g, const mg_lchain_t *l0, const mg_lchain_t *l1)
{
	/**
	 case 1: l0->qs************l0->qe
	 					l1->qs************l1->qe
	 case 2: l0->qs************l0->qe
	                                       l1->qs************l1->qe

	*****l0->rs************l0->re*****
				                                     ****l1->rs************l1->re**
	 * **/
	///min_dist = l1->rs + (g->seg[l0->v>>1].len - l0->re); 
	// below equals (l1->qs - l0->qe) - min_dist + g->seg[l1->v>>1].len; see mg_gchain1_dp() for the calculation of min_dist
	//(l1->qs - l0->qe) is the gap in query, min_dist is the gap in reference
	return (l1->qs - l0->qe) - (g->seq[l0->v>>1].len - l0->re) + (g->seq[l1->v>>1].len - l1->rs);
	// when l0->v == l1->v, the above becomes (l1->qs - l0->qe) - (l1->rs - l0->re), which is what we want
}

static inline sp_node_t *gen_sp_node(void *km, uint32_t v, int32_t d, int32_t id)
{
	sp_node_t *p;
	KMALLOC(km, p, 1);
	p->v = v, p->di = (uint64_t)d<<32 | id, p->pre = -1, p->is_0 = 1;
	return p;
}

///max_dist is like the overlap length in string graph
///mg_shortest_k(km, g, li->v^1, n_dst, dst, max_dist_g + (g->seg[li->v>>1].len - li->rs), MG_MAX_SHORT_K, 0, 0, 1, 0);
///the end position of qs is li->qs; dst[]->->qlen indicate the region that need to be checked in bases
mg_pathv_t *mg_shortest_k(void *km0, const asg_t *g, uint32_t src, int32_t n_dst, mg_path_dst_t *dst, int32_t max_dist, int32_t max_k, int32_t ql, const char *qs, int is_rev, int32_t *n_pathv)
{
	sp_node_t *p, *root = 0, **out;
	sp_topk_t *q;
	khash_t(sp) *h;
	khash_t(sp2) *h2;
	void *km;
	khint_t k;
	int absent;
	int32_t i, j, n_done, n_found, n_seeds = 0;
	uint32_t id, n_out, m_out;
	int8_t *dst_done;
	mg_pathv_t *ret = 0;
	uint64_t *dst_group, *seeds = 0;
	/** //path
	void *h_seeds = 0; 
	mg128_v mini = {0,0,0}; 
	**/

	if (n_pathv) *n_pathv = 0;///for us, n_pathv = NULL
	if (n_dst <= 0) return 0;///n_dst: how many candidate nodes
	for (i = 0; i < n_dst; ++i) { // initialize
		mg_path_dst_t *t = &dst[i];
		///if src and dest are at the same ref id, there are already one path
		if (t->inner)///if two chains are at the same ref id
			t->dist = 0, t->n_path = 1, t->path_end = -1;
		else
			t->dist = -1, t->n_path = 0, t->path_end = -1;
	}
	if (max_k > MG_MAX_SHORT_K) max_k = MG_MAX_SHORT_K;
	km = km_init2(km0, 0x4000);

	/** //path
	///for the first time, we just check th reachability without sequence (qs);
	///but for the second round, we need to check sequence
	///qs is the sequence between two minimizers
	if (ql > 0 && qs) { // build the seed hash table for the query
		mg_sketch(km, qs, ql, MG_SHORT_KW, MG_SHORT_KK, 0, &mini);
		// mini->a[].x = hash_key<<8 | kmerSpan
		// mini->a[].y = rid<<32 | lastPos<<1 | strand
		if (is_rev)///is_rev = 1;
			for (i = 0; i < mini.n; ++i)///reverse qs[0, ql) to qs(ql, 0]
				mini.a[i].y = (ql - (((int32_t)mini.a[i].y>>1) + 1 - MG_SHORT_KK) - 1) << 1 | ((mini.a[i].y&1)^1);
		///h_seeds is the ordinary hash index
		h_seeds = mg_idx_a2h(km, mini.n, mini.a, 0, &seeds, &n_seeds);
		///h_seeds+seeds+n_seeds ----> hash index of qs[0, ql)
	}
	**/

	KCALLOC(km, dst_done, n_dst);
	KMALLOC(km, dst_group, n_dst);
	// multiple dst[] may have the same dst[].v. We need to group them first.
	// in other words, one ref id may have multiple dst alignment chains
	for (i = 0; i < n_dst; ++i) 
		dst_group[i] = (uint64_t)dst[i].v<<32 | i;
	radix_sort_gfa64(dst_group, dst_group + n_dst);

	h2 = kh_init2(sp2, km); // this hash table keeps all destinations from the same ref id
	kh_resize(sp2, h2, n_dst * 2);
	///please note that one contig in ref may have multiple alignment chains
	///so h2 is a index that helps us to query it
	///key(h2) = ref id; value(h2) = start_idx | occ
	for (i = 1, j = 0; i <= n_dst; ++i) {
		if (i == n_dst || dst_group[i]>>32 != dst_group[j]>>32) {
			k = kh_put(sp2, h2, dst_group[j]>>32, &absent);
			kh_val(h2, k) = (uint64_t)j << 32 | (i - j);
			assert(absent);
			j = i;
		}
	}

	

	h = kh_init2(sp, km); // this hash table keeps visited vertices; path to each visited vertice
	kh_resize(sp, h, 16);

	m_out = 16, n_out = 0;///16 is just the initial size
	KMALLOC(km, out, m_out);
	
	/**
	typedef struct {
		int32_t k, mlen;//k: number of walks from src to this node
		int32_t qs, qe;
		sp_node_t *p[MG_MAX_SHORT_K]; // this forms a max-heap; all path
	} sp_topk_t;
	**/
	id = 0;
	p = gen_sp_node(km, src, 0, id++);///just malloc a node for src; the distance is 0
	p->hash = __ac_Wang_hash(src);
	kavl_insert(sp, &root, p, 0);///should be avl tree

	k = kh_put(sp, h, src, &absent);///here is a hash table
	q = &kh_val(h, k);
	q->k = 1, q->p[0] = p, q->mlen = 0, q->qs = q->qe = -1;

	n_done = 0;
	///the key of avl tree: #define sp_node_cmp(a, b) (((a)->di > (b)->di) - ((a)->di < (b)->di))
	///the higher bits of (*)->di is distance to src node
	///so the key of avl tree is distance
	while (kavl_size(head, root) > 0) {///thr first root is src
		int32_t i, nv;
		asg_arc_t *av;
		sp_node_t *r;
		///note that one node might be visited multiple times if there are circles
		///delete the first node
		r = kavl_erase_first(sp, &root); // take out the closest vertex in the heap (as a binary tree)
		//fprintf(stderr, "XX\t%d\t%d\t%d\t%c%s[%d]\t%d\n", n_out, kavl_size(head, root), n_finished, "><"[(r->v&1)^1], g->seg[r->v>>1].name, r->v, (int32_t)(r->di>>32));
		if (n_out == m_out) KEXPAND(km, out, m_out);
		///higher 32 bits might be the distance to root node
		// lower 32 bits now for position in the out[] array
		r->di = r->di>>32<<32 | n_out; ///n_out is just the id in out
		///so one node id in graph might be saved multiple times in avl tree and out[]
		out[n_out++] = r;///out[0] = src

		///r->v is the dst vertex id
		///sometimes k==kh_end(h2). Some nodes are found by graph travesal but not in linear chain alignment
		k = kh_get(sp2, h2, r->v);
		// we have reached one dst vertex
		// note that one dst vertex may have multipe alignment chains
		// we can visit some nodes in graph which are not reachable during chaining
		// h2 is used to determine if one node is reachable or not
		if (k != kh_end(h2)) { 
			///node r->v might be visited multiple times
			int32_t j, dist = r->di>>32, off = kh_val(h2, k) >> 32, cnt = (int32_t)kh_val(h2, k);
			//src can reach ref id r->v; there might be not only one alignment chain in r->v
			//so we need to scan all of them
			for (j = 0; j < cnt; ++j) {
				mg_path_dst_t *t = &dst[(int32_t)dst_group[off + j]];
				int32_t done = 0;
				///the src and dest are at the same ref id, say we directly find the shortest path
				if (t->inner) {
					done = 1;
				} else {
					int32_t mlen = 0, copy = 0;
					///in the first round, we just check reachability without sequence
					///so h_seeds = NULL; we can assume mlen = 0
					/** //path
					mlen = h_seeds? path_mlen(out, n_out - 1, h, t->qlen) : 0;
					**/
					//if (mg_dbg_flag & MG_DBG_GC1) fprintf(stderr, "  src=%c%s[%d],qlen=%d\tdst=%c%s[%d]\ttarget_distx=%d,target_hash=%x\tdistx=%d,mlen=%d,hash=%x\n", "><"[src&1], g->seg[src>>1].name, src, ql, "><"[t->v&1], g->seg[t->v>>1].name, t->v, t->target_dist - g->seg[src>>1].len, t->target_hash, dist - g->seg[src>>1].len, mlen, r->hash);
					// note: t indicates a linear alignmnet, instead of a node in graph
					///target_dist should be the distance on query
					if (t->n_path == 0) { // keep the shortest path
						copy = 1;
					} else if (t->target_dist >= 0) { // we have a target distance; choose the closest
						if (dist == t->target_dist && t->check_hash && r->hash == t->target_hash) { // we found the target path
							copy = 1, done = 1;
						} else {
							int32_t d0 = t->dist, d1 = dist;
							d0 = d0 > t->target_dist? d0 - t->target_dist : t->target_dist - d0;
							d1 = d1 > t->target_dist? d1 - t->target_dist : t->target_dist - d1;
							///if the new distance (d1) is smaller than the old distance (d0), update the results
							///in other words, the length of new path should be closer to t->target_dist
							if (d1 - mlen/2 < d0 - t->mlen/2) copy = 1;
						}
					}
					if (copy) {
						t->path_end = n_out - 1, t->dist = dist, t->hash = r->hash, t->mlen = mlen, t->is_0 = r->is_0;
						if (t->target_dist >= 0) {
							///src is from li from li to lj, so the dis is generally increased
							///target_dist should be the distance on query
							if (dist == t->target_dist && t->check_hash && r->hash == t->target_hash) done = 1;
							else if (dist > t->target_dist + MG_SHORT_K_EXT) done = 1;
						}
					}
					++t->n_path;///we found a path to the alignment t
					if (t->n_path >= max_k) done = 1;
				}
				if (dst_done[off + j] == 0 && done)
					dst_done[off + j] = 1, ++n_done;
			}
			///if all alignments have been settle down
			///pre-end; accelerate the loop
			if (n_done == n_dst) break;
		}

		///below is used to push new nodes to avl tree for iteration
		nv = asg_arc_n(g, r->v);
		av = asg_arc_a(g, r->v);
		for (i = 0; i < nv; ++i) { // visit all neighbors
			asg_arc_t *ai = &av[i];
			///v_lv is the (dest_length - overlap_length); it is a normal path length in string graph
			///ai->v_lv is the path length from r->v to ai->w
			///(r->di>>32)
			int32_t d = (r->di>>32) + (uint32_t)ai->ul;
			if (d > max_dist) continue; // don't probe vertices too far away
			///ai->w is the dest ref id; we insert a new ref id, instead of an alignment chain
			k = kh_put(sp, h, ai->v, &absent);///one node might be visited multiple times
			q = &kh_val(h, k);
			if (absent) { // a new vertex visited
				///q->k: number of walks from src to ai->w
				q->k = 0, q->qs = q->qe = -1; q->mlen = 0;
				///h_seeds = NULL; so q->mlen = 0
				/** //path
				q->mlen = h_seeds && d + gfa_arc_lw(g, *ai) <= max_dist? node_mlen(km, g, ai->w, &mini, h_seeds, n_seeds, seeds, &q->qs, &q->qe) : 0;
				**/
				//if (ql && qs) fprintf(stderr, "ql=%d,src=%d\tv=%c%s[%d],n_seeds=%d,mlen=%d\n", ql, src, "><"[ai->w&1], g->seg[ai->w>>1].name, ai->w, n_seeds, q->mlen);
			}
			///if there are less than <max_k> walks from src to ai->w, directly add
			///if there are more, keep the smallest <max_k> walks
			if (q->k < max_k) { // enough room: add to the heap
				p = gen_sp_node(km, ai->v, d, id++);
				p->pre = n_out - 1;///the parent node of this one 
				p->hash = r->hash + __ac_Wang_hash(ai->v);
				p->is_0 = r->is_0;
				/** //path
				if (ai->rank > 0) p->is_0 = 0;
				**/
				kavl_insert(sp, &root, p, 0);
				q->p[q->k++] = p;
				ks_heapup_sp(q->k, q->p);///adjust heap by distance
			} else if (q->p[0]->di>>32 > d) { // shorter than the longest path so far: replace the longest
				p = kavl_erase(sp, &root, q->p[0], 0);
				if (p) {
					p->di = (uint64_t)d<<32 | (id++);
					p->pre = n_out - 1;
					p->hash = r->hash + __ac_Wang_hash(ai->v);
					p->is_0 = r->is_0;
					/** //path
					if (ai->rank > 0) p->is_0 = 0;
					**/
					kavl_insert(sp, &root, p, 0);
					ks_heapdown_sp(0, q->k, q->p);
				} else {
					fprintf(stderr, "Warning: logical bug in gfa_shortest_k(): q->k=%d,q->p[0]->{d,i}={%d,%d},d=%d,src=%u,max_dist=%d,n_dst=%d\n", q->k, (int32_t)(q->p[0]->di>>32), (int32_t)q->p[0]->di, d, src, max_dist, n_dst);
					km_destroy(km);
					return 0;
				}
			} // else: the path is longer than all the existing paths ended at ai->w
		}
	}

	kfree(km, dst_group);
	kfree(km, dst_done);
	kh_destroy(sp, h);
	mg_idx_hfree(h_seeds);
	kfree(km, seeds);
	kfree(km, mini.a);
	// NB: AVL nodes are not deallocated. When km==0, they are memory leaks.

	for (i = 0, n_found = 0; i < n_dst; ++i)
		if (dst[i].n_path > 0) ++n_found;///n_path might be larger than 16
	///we can assume n_pathv = NULL for now 
	if (n_found > 0 && n_pathv) { // then generate the backtrack array
		int32_t n, *trans;
		KCALLOC(km, trans, n_out); // used to squeeze unused elements in out[]
		for (i = 0; i < n_dst; ++i) { // mark dst vertices with a target distance
			mg_path_dst_t *t = &dst[i];
			if (t->n_path > 0 && t->target_dist >= 0 && t->path_end >= 0)
				trans[(int32_t)out[t->path_end]->di] = 1;
		}
		for (i = 0; i < n_out; ++i) { // mark dst vertices without a target distance
			k = kh_get(sp2, h2, out[i]->v);
			if (k != kh_end(h2)) { // TODO: check if this is correct!
				int32_t off = kh_val(h2, k)>>32, cnt = (int32_t)kh_val(h2, k);
				for (j = off; j < off + cnt; ++j)
					if (dst[j].target_dist < 0)
						trans[i] = 1;
			}
		}
		for (i = n_out - 1; i >= 0; --i) // mark all predecessors
			if (trans[i] && out[i]->pre >= 0)
				trans[out[i]->pre] = 1;
		for (i = n = 0; i < n_out; ++i) // generate coordinate translations
			if (trans[i]) trans[i] = n++;
			else trans[i] = -1;

		*n_pathv = n;
		KMALLOC(km0, ret, n);
		for (i = 0; i < n_out; ++i) { // generate the backtrack array
			mg_pathv_t *p;
			if (trans[i] < 0) continue;
			p = &ret[trans[i]];
			p->v = out[i]->v, p->d = out[i]->di >> 32;
			p->pre = out[i]->pre < 0? out[i]->pre : trans[out[i]->pre];
		}
		for (i = 0; i < n_dst; ++i) // translate "path_end"
			if (dst[i].path_end >= 0)
				dst[i].path_end = trans[dst[i].path_end];
	}

	km_destroy(km);
	return ret;
}


int32_t mg_gchain1_dp(void *km, const ma_ug_t *ug, int32_t *n_lc_, mg_lchain_t *lc, int32_t qlen, int32_t max_dist_g, int32_t max_dist_q, int32_t bw, int32_t max_skip,
					  int32_t ref_bonus, float chn_pen_gap, float chn_pen_skip, float mask_level, int32_t max_gc_seq_ext, const char *qseq, const mg128_t *an, uint64_t **u_)
{
	int32_t i, j, k, m_dst, n_dst, n_ext, n_u, n_v, n_lc = *n_lc_;
	int32_t *f, *v, *t;
	int64_t *p;
	uint64_t *u;
	mg_path_dst_t *dst;
	gc_frag_t *a;
	mg_lchain_t *swap;
	char *qs;
	asg_t *g = ug->g;

	*u_ = 0;
	if (n_lc == 0) return 0;

	KMALLOC(km, a, n_lc);
	///n_lc how many linear chains
	for (i = n_ext = 0; i < n_lc; ++i) { // a[] is a view of frag[]; for sorting
		mg_lchain_t *r = &lc[i];
		gc_frag_t *ai = &a[i];
		int32_t is_isolated = 0, min_end_dist_g;
		r->dist_pre = -1;///indicate parent in graph chain
		min_end_dist_g = g->seq[r->v>>1].len - r->re;///r->v: ref_id|rev
		if (r->rs < min_end_dist_g) min_end_dist_g = r->rs;
		if (min_end_dist_g > max_dist_g) is_isolated = 1; // if too far from segment ends
		else if (min_end_dist_g>>3 > r->score) is_isolated = 1; // if the lchain too small relative to distance to the segment ends
		ai->srt = (uint32_t)is_isolated<<31 | r->qe;
		ai->i = i;
		if (!is_isolated) ++n_ext;
	}
	///if the alignment is too far from segment ends, which means it cannot contribute to graph alignment
	if (n_ext < 2) { // no graph chaining needed; early return
		kfree(km, a);
		KMALLOC(km, u, n_lc);
		for (i = 0; i < n_lc; ++i)
			u[i] = (uint64_t)lc[i].score<<32 | 1;
		*u_ = u;
		return n_lc;
	}
	radix_sort_gc(a, a + n_lc);///sort by: is_isolated(1):qe

	KMALLOC(km, v, n_lc);
	KMALLOC(km, f, n_ext);
	KMALLOC(km, p, n_ext);
	KCALLOC(km, t, n_ext);
	// KMALLOC(km, qs, max_dist_q + 1);//for 

	m_dst = n_dst = 0, dst = 0;
	///n_ext is number of linear chains that might be included in graph chains
	///sorted by the positions in query; sorted by qe of each chain
	for (i = 0; i < n_ext; ++i) { // core loop
		gc_frag_t *ai = &a[i];
		mg_lchain_t *li = &lc[ai->i];///linear chain; sorted by qe, i.e. end position in query
		///note segi is query id, instead of ref id; it is not such useful
		/**
		 * a[].x: idx_in_minimizer_arr(32)r_pos(32)
		 * a[].y: weight(8)query_id(8)flag(8)span(8)q_pos(32)
		**/
		{ // collect end points potentially reachable from _i_
			int32_t x = li->qs + bw, n_skip = 0;
			if (x > qlen) x = qlen;
			///collect alignments that can be reachable from the left side
			///that is, a[x].qe <= x
			x = find_max(i, a, x);
			n_dst = 0;
			for (j = x; j >= 0; --j) { // collect potential destination vertices
				gc_frag_t *aj = &a[j];
				//potential chains that might be overlapped with the left side of li
				mg_lchain_t *lj = &lc[aj->i];
				mg_path_dst_t *q;
				int32_t target_dist, dq;
				///lj->qs >= li->qs && lj->qe <= li->qs, so lj is contained
				if (lj->qs >= li->qs) continue; // lj is contained in li on the query coordinate
				///lj->qs************lj->qe
				///			li->qs************li->qe
				if (lj->qe > li->qs) { // test overlap on the query
					int o = lj->qe - li->qs;
					///mask_level = 0.5, if overlap is too long
					///note here is the overlap in query, so too long overlaps might be wrong
					if (o > (lj->qe - lj->qs) * mask_level || o > (li->qe - li->qs) * mask_level)
						continue;
				}
				dq = li->qs - lj->qe;///dq might be smaller than 0
				if (dq > max_dist_q) break; // if query gap too large, stop
				///The above filter chains like:
				///1. lj is contained in li
				///2. the overlap between li and lj is too large
				///3. li and lj are too far
				/**
				 lj->qs************lj->qe
				 			li->qs************li->qe

				*****lj->rs************lj->re*****
				                                     ****li->rs************li->re**
				**/
				///above we have checked gap/overlap in query 
				///then we need to check gap/overlap in reference
				if (li->v != lj->v) { // the two linear chains are on two different refs
					// minimal graph gap; the real graph gap might be larger
					int32_t min_dist = li->rs + (g->seq[lj->v>>1].len - lj->re); 
					if (min_dist > max_dist_g) continue; // graph gap too large
					//note here min_dist - (lj->qs - li->qe) > bw is important
					//min_dist is always larger than 0, (lj->qs - li->qe) might be negative 
					if (min_dist - bw > li->qs - lj->qe) continue; ///note seg* is the query id, instead of ref id
					target_dist = mg_target_dist(g, lj, li);
					if (target_dist < 0) continue; // this may happen if the query overlap is far too large
				} else if (lj->rs >= li->rs || lj->re >= li->re) { // not colinear
					continue;
				} else {///li->v == lj->v and colinear; at the same ref id
						/**
						 	case 1: lj->qs************lj->qe
												li->qs************li->qe
							case 2: lj->qs************lj->qe
																li->qs************li->qe

							*****lj->rs************lj->re*****
																			****li->rs************li->re**
						* **/
					///w is indel, w is always positive
					int32_t dr = li->rs - lj->re, w = dr > dq? dr - dq : dq - dr;
					///note that l*->v is the ref id, while seg* is the query id
					if (w > bw) continue; // test bandwidth
					if (dr > max_dist_g || dr < -max_dist_g) continue;
					if (lj->re > li->rs) { // test overlap on the graph segment
						int o = lj->re - li->rs;
						if (o > (lj->re - lj->rs) * mask_level || o > (li->re - li->rs) * mask_level)
							continue;
					}
					target_dist = mg_target_dist(g, lj, li);
				}
				if (n_dst == m_dst) KEXPAND(km, dst, m_dst); // TODO: watch out the quadratic behavior!
				q = &dst[n_dst++];///q saves information for i->j
				memset(q, 0, sizeof(mg_path_dst_t));
				///note v is (rid:rev), so two alignment chains might be at the same ref id with different directions
				q->inner = (li->v == lj->v);
				q->v = lj->v^1;///must be v^1 instead of v
				q->meta = j;
				///lj->qs************lj->qe
				///			li->qs************li->qe
				q->qlen = li->qs - lj->qe;///might be negative
				q->target_dist = target_dist;///cannot understand the target_dist
				q->target_hash = 0;
				q->check_hash = 0;
				if (t[j] == i) {///this pre-cut is weird; attention
					if (++n_skip > max_skip)
						break;
				}
				if (p[j] >= 0) t[p[j]] = i;
			}
		}
		///the above saves all linear chains that might be reached to the left side of chain i
		///all those chains are saved to dst<n_dst>
		{ // confirm reach-ability
			int32_t k;
			// test reach-ability without sequences
			/**
			*****lj->rs************lj->re*****
				                                     ****li->rs************li->re***
			(g->seg[li->v>>1].len - li->rs) ----> is like the node length in string graph
	 		**/
			mg_shortest_k(km, g, li->v^1, n_dst, dst, max_dist_g + (g->seq[li->v>>1].len - li->rs), MG_MAX_SHORT_K, 0, 0, 1, 0);
			// remove unreachable destinations
			for (j = k = 0; j < n_dst; ++j) {
				mg_path_dst_t *dj = &dst[j];
				int32_t sc;
				if (dj->n_path == 0) continue; // unreachable
				sc = cal_sc(dj, li, lc, an, a, f, bw, ref_bonus, chn_pen_gap);
				if (sc == INT32_MIN) continue; // out of band
				if (sc + li->score < 0) continue; // negative score and too low
				dst[k] = dst[j];
				dst[k++].srt_key = INT64_MAX/2 - (int64_t)sc; // sort in the descending order
			}
			n_dst = k;
			if (n_dst > 0) {
				radix_sort_dst(dst, dst + n_dst);
				// discard weaker chains if the best score is much larger (assuming base-level heuristic won't lift it to the top chain)
				// dst[0].srt_key has the largest score
				for (j = 1; j < n_dst; ++j) 
					if (dst[j].srt_key - dst[0].srt_key > li->score)//discard chains with too small weight
						break;
				n_dst = j;
				if (n_dst > max_gc_seq_ext) n_dst = max_gc_seq_ext; // discard weaker chains
			}
		}
		if (n_dst > 0) { // find paths with sequences
			int32_t min_qs = li->qs;
			for (j = 0; j < n_dst; ++j) {
				const mg_lchain_t *lj;
				assert(dst[j].n_path > 0);
				///a[]->srt = (uint32_t)is_isolated<<31 | r->qe;
				///a[]->i = i;
				lj = &lc[a[dst[j].meta].i];
				if (lj->qe < min_qs) min_qs = lj->qe;
			}
			///qs keeps the sequence at the gap between the li and lj in query
			memcpy(qs, &qseq[min_qs], li->qs - min_qs);
			mg_shortest_k(km, g, li->v^1, n_dst, dst, max_dist_g + (g->seg[li->v>>1].len - li->rs), MG_MAX_SHORT_K, li->qs - min_qs, qs, 1, 0);
			if (mg_dbg_flag & MG_DBG_GC1) fprintf(stderr, "[src:%d] q_intv=[%d,%d), src=%c%s[%d], n_dst=%d, max_dist=%d, min_qs=%d, lc_score=%d\n", ai->i, li->qs, li->qe, "><"[(li->v&1)^1], g->seg[li->v>>1].name, li->v^1, n_dst, max_dist_g + (g->seg[li->v>>1].len - li->rs), min_qs, li->score);
		}
		{ // DP
			int32_t max_f = li->score, max_j = -1, max_d = -1, max_inner = 0;
			uint32_t max_hash = 0;
			for (j = 0; j < n_dst; ++j) {
				mg_path_dst_t *dj = &dst[j];
				int32_t sc;
				sc = cal_sc(dj, li, lc, an, a, f, bw, ref_bonus, chn_pen_gap);
				if (sc == INT32_MIN) continue;
				if (mg_dbg_flag & MG_DBG_GC1) {
					mg_lchain_t *lj = &lc[a[dj->meta].i];
					fprintf(stderr, "  [dst:%d] dst=%c%s[%d], n_path=%d, target=%d, opt_dist=%d, score=%d, q_intv=[%d,%d), g_intv=[%d,%d)\n", dj->meta, "><"[dj->v&1], g->seg[dj->v>>1].name, dj->v, dj->n_path, dj->target_dist - g->seg[li->v>>1].len, dj->dist - g->seg[li->v>>1].len, sc, lj->qs, lj->qe, lj->rs, lj->re);
				}
				if (sc > max_f) max_f = sc, max_j = dj->meta, max_d = dj->dist, max_hash = dj->hash, max_inner = dj->inner;
			}
			f[i] = max_f, p[i] = max_j;
			li->dist_pre = max_d;
			li->hash_pre = max_hash;
			li->inner_pre = max_inner;
			v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f;
			if (mg_dbg_flag & MG_DBG_GC1) fprintf(stderr, " [opt:%d] opt=%d, max_f=%d\n", ai->i, max_j, max_f);
		}
	}
	kfree(km, dst);
	kfree(km, qs);
	if (mg_dbg_flag & MG_DBG_GC1) {
		int32_t mmax_f = 0, mmax_i = -1;
		for (i = 0; i < n_ext; ++i) if (f[i] > mmax_f) mmax_f = f[i], mmax_i = i;
		i = mmax_i; while (i >= 0) { fprintf(stderr, "[best] i=%d, seg=%s, max_f=%d, chn_pen_gap=%f\n", a[i].i, g->seg[lc[a[i].i].v>>1].name, f[i], chn_pen_gap); i = p[i]; }
	}
	///n_ext: number of useful chains
	///n_lc - n_ext: number of isoated chains
	u = mg_chain_backtrack(km, n_ext, f, p, v, t, 0, 0, n_lc - n_ext, &n_u, &n_v);
	kfree(km, f); kfree(km, p); kfree(km, t);
	///store the extra isoated chains
	for (i = 0; i < n_lc - n_ext; ++i) {
		u[n_u++] = (uint64_t)lc[a[n_ext + i].i].score << 32 | 1;
		v[n_v++] = n_ext + i;
	}

	///reorganize lc;
	KMALLOC(km, swap, n_v);
	for (i = 0, k = 0; i < n_u; ++i) {
		int32_t k0 = k, ni = (int32_t)u[i];
		for (j = 0; j < ni; ++j)
			swap[k++] = lc[a[v[k0 + (ni - j - 1)]].i];
	}
	assert(k == n_v);
	memcpy(lc, swap, n_v * sizeof(mg_lchain_t));
	*n_lc_ = n_v;
	*u_ = u;

	kfree(km, a);
	kfree(km, swap);
	kfree(km, v);
	return n_u;
}


void mg_map_frag(const void *ha_flt_tab, const ha_pt_t *ha_idx, const ma_ug_t *ug, const uint32_t qid, const int qlen, const char *qseq, ha_mzl_v *mz, 
st_mt_t *sp, mg_tbuf_t *b, int32_t w, int32_t k, int32_t hpc, int32_t mz_sd, int32_t mz_rewin, const mg_idxopt_t *opt)
{
    mg128_t *a = NULL;
    int64_t n_a;
    int32_t *mini_pos;
    int i, rep_len, n_mini_pos, n_lc, max_chain_gap_qry, max_chain_gap_ref;
    uint64_t *u;
	mg_lchain_t *lc;
    mz->n = 0;
    mz2_ha_sketch(qseq, qlen, w, k, 0, hpc, mz, ha_flt_tab, mz_sd, NULL, NULL, NULL, -1, -1, -1, sp, mz_rewin, 1);
    ///a[]->y: weight(8)seg_id(8)flag(8)span(8)pos(32);--->query
    ///a[]->x: rid(31)rev(1)rpos(33);--->reference
    a = collect_seed_hits(b->km, opt, opt->hap_n, ha_flt_tab, ha_idx, ug, mz, &n_a, &rep_len, &n_mini_pos, &mini_pos);
    if (opt->max_gap_pre > 0 && opt->max_gap_pre * 2 < opt->max_gap) n_a = flt_anchors(n_a, a, opt->max_gap_pre);
    max_chain_gap_qry = max_chain_gap_ref = opt->max_gap;
    if (n_a == 0) {
        if(a) kfree(b->km, a);
        a = 0, n_lc = 0, u = 0;
    } else {
		a = mg_lchain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_lc_skip, opt->max_lc_iter, 
		opt->min_lc_cnt, opt->min_lc_score, opt->chn_pen_gap, n_a, a, &n_lc, &u, b->km);
    }

	if (n_lc) {///n_lc is how many chain we found
		lc = mg_lchain_gen(b->km, qlen, n_lc, u, a);
		for (i = 0; i < n_lc; ++i)///update a[] since ref_id|rev has already been saved to lc[].v
			mg_update_anchors(lc[i].cnt, &a[lc[i].off], n_mini_pos, mini_pos);///update a[].x
	} else lc = 0;
	kfree(b->km, mini_pos); kfree(b->km, u);

	/**
	 * up to here, a[] has been changed
	 * a[].x: idx_in_minimizer_arr(32)r_pos(32)
	 * a[].y: weight(8)query_id(8)flag(8)span(8)q_pos(32)
	**/
	
}

static void worker_for_ul_alignment(void *data, long i, int tid) // callback for kt_for()
{
    utepdat_t *s = (utepdat_t*)data;
    // ma_ov_buf_t *b = s->mo[tid];
    mg_map_frag(s->ha_flt_tab, s->ha_idx, s->ug, s->id+i, s->len[i], s->seq[i], &(s->mzs[tid]), &(s->sps[tid]), s->buf[tid], s->opt->w, s->opt->k, 
    s->opt->is_HPC, asm_opt.mz_sample_dist, asm_opt.mz_rewin, s->opt);

}

static void *worker_ul_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
    uldat_t *p = (uldat_t*)data;
    ///uint64_t total_base = 0, total_pair = 0;
    if (step == 0) { // step 1: read a block of sequences
        int ret;
        uint64_t l;
		utepdat_t *s;
		CALLOC(s, 1);
        s->ha_flt_tab = p->ha_flt_tab; s->ha_idx = p->ha_idx; s->id = p->total_pair; s->opt = p->opt; s->ug = p->ug;
        while ((ret = kseq_read(p->ks)) >= 0) 
        {
            if (p->ks->seq.l < (uint64_t)p->opt->k) continue;
            if (s->n == s->m) {
                s->m = s->m < 16? 16 : s->m + (s->n>>1);
                REALLOC(s->len, s->m);
                REALLOC(s->seq, s->m);
            }
            l = p->ks->seq.l;
            MALLOC(s->seq[s->n], l);
            s->sum_len += l;
            memcpy(s->seq[s->n], p->ks->seq.s, l);
            s->len[s->n++] = l;
            if (s->sum_len >= p->chunk_size) break;            
        }
        p->total_pair += s->n;
        if (s->sum_len == 0) free(s);
		else return s;
    }
    else if (step == 1) { // step 2: alignment
		uint64_t i;
        utepdat_t *s = (utepdat_t*)in;
        CALLOC(s->mzs, p->n_thread);
        CALLOC(s->sps, p->n_thread);

		s->buf = (mg_tbuf_t**)calloc(p->n_thread, sizeof(mg_tbuf_t*));
		for (i = 0; i < p->n_thread; ++i) s->buf[i] = mg_tbuf_init();
        /**
        CALLOC(s->pos, s->n);
        **/
        kt_for(p->n_thread, worker_for_ul_alignment, s, s->n);
        for (i = 0; i < (uint64_t)s->n; ++i) {
            free(s->seq[i]);
            p->total_base += s->len[i];
        }
        free(s->seq); free(s->len);

		for (i = 0; i < p->n_thread; ++i) mg_tbuf_destroy(s->buf[i]);
		free(s->buf);
		return s;
    }
    else if (step == 2) { // step 3: dump
        utepdat_t *s = (utepdat_t*)in;
        /**
        int i;
        for (i = 0; i < s->n; ++i) {
            // if(s->pos[i].a == NULL) continue;
            // kv_push(pe_hit_hap, p->hits, s->pos[i]);
            if(s->pos[i].s == (uint64_t)-1) continue;
            kv_push(pe_hit, p->hits.a, s->pos[i]);
        }
        free(s->pos);
        **/
        free(s);
    }
    return 0;
}

int alignment_ul_pipeline(uldat_t* sl, const enzyme *fn)
{
    double index_time = yak_realtime();
    int i;
    for (i = 0; i < fn->n; i++){
        gzFile fp;
        if ((fp = gzopen(fn->a[i], "r")) == 0) return 0;
        sl->ks = kseq_init(fp);
        kt_pipeline(3, worker_ul_pipeline, sl, 3);
        kseq_destroy(sl->ks);
        gzclose(fp);
    }
    fprintf(stderr, "[M::%s::%.3f] ==> Qualification\n", __func__, yak_realtime()-index_time);
    return 1;
}

int ul_align(mg_idxopt_t *opt, const enzyme *fn, void *ha_flt_tab, ha_pt_t *ha_idx, ma_ug_t *ug)
{
    uldat_t sl; memset(&sl, 0, sizeof(sl));
    sl.ha_flt_tab = ha_flt_tab;
    sl.ha_idx = ha_idx;
    sl.opt = opt;   
    sl.chunk_size = 20000000;
    sl.n_thread = asm_opt.thread_num;
    sl.ug = ug;
    alignment_ul_pipeline(&sl, fn);
    return 1;
}

void ul_resolve(ma_ug_t *ug, int hap_n)
{
    mg_idxopt_t opt;
    init_mg_opt(&opt, 0, 19, 10, hap_n);
    uidx_build(ug, &opt);
    ul_align(&opt, asm_opt.ar, ha_flt_tab, ha_idx, ug);




    uidx_destory();
}