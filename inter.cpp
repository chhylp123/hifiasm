#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <zlib.h>
#include <math.h>
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
#include "Correct.h"
#include "Process_Read.h"
#include "Assembly.h"
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
	int w, k, bw, max_gap, is_HPC, hap_n, occ_weight, max_gap_pre, max_gc_seq_ext, seed;
	int max_lc_skip, max_lc_iter, min_lc_cnt, min_lc_score, max_gc_skip, ref_bonus;
	int min_gc_cnt, min_gc_score, sub_diff, best_n;
    float chn_pen_gap, mask_level, pri_ratio;
} mg_idxopt_t;

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

typedef struct {
	int32_t qs, qe, rs, re;
	uint32_t v;
} mg_coor_t;

///mg128_t->y: weight(8)seg_id(8)flag(8)span(8)pos(32)
///mg128_t->x: rid(31)rev(1)pos(33); keep reference
typedef struct { uint64_t x, y; } mg128_t;
#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, mg128_t, sort_key_128x, 8) 
void radix_sort_128x(mg128_t *beg, mg128_t *end);

typedef struct {
	int32_t off, cnt;
	uint32_t v;
	int32_t score;
} mg_llchain_t;

typedef struct {
	int32_t id, parent;
	int32_t off, cnt;
	int32_t n_anchor, score;
	int32_t qs, qe;
	int32_t plen, ps, pe;
	int32_t blen, mlen;
	float div;
	uint32_t hash;
	int32_t subsc, n_sub;
	uint32_t mapq:8, flt:1, dummy:23;
} mg_gchain_t;

typedef struct {
	size_t n,m;
    uint64_t *a, tl;
	kvec_t(char) cc;
} mg_dbn_t;

typedef struct {
	int32_t cnt;
	uint32_t v;
	int32_t score;
	uint32_t qs, qe, ts, te;
} mg_lres_t;

typedef struct {
    int32_t n_gc, n_lc;
    mg_gchain_t *gc;///g_chain; idx in l_chains
    mg_lres_t *lc;///l_chain
    uint64_t qid, qlen;
} mg_gres_t;

typedef struct {
	size_t n,m;
    mg_gres_t *a;
	uint64_t total_base;
    uint64_t total_pair;
} mg_gres_a;

typedef struct { // global data structure for kt_pipeline()
    const void *ha_flt_tab;
    const ha_pt_t *ha_idx;
    const mg_idxopt_t *opt;
    const ma_ug_t *ug;
	const asg_t *rg;
	const ug_opt_t *uopt;
	kseq_t *ks;
    int64_t chunk_size;
    uint64_t n_thread;
    uint64_t total_base;
    uint64_t total_pair;
	mg_gres_a hits;
	mg_dbn_t nn;
} uldat_t;

typedef struct {
	uint64_t asm_size;
    uint64_t asm_cov;
} mul_ov_t;

///three levels: 
///level-0: minimizers
///level-1: linear chains
///level-2: g chains
///gc[] saves the idx in lc[], lc saves the idx in a[]
typedef struct {
	void *km;
	int32_t n_gc, n_lc, n_a, rep_len;
	mg_gchain_t *gc;///g_chain; idx in l_chains
	mg_llchain_t *lc;///l_chain
	mg128_t *a; // minimizer positions; see comments above mg_update_anchors() for details
	uint64_t qid, qlen;
} mg_gchains_t;

typedef struct {
	uint32_t n; ///length of candidate list 
    uint64_t q_span:31, rev:1, q_pos:32;
	uint32_t qid:16, weight:15, is_tandem:1;
	const ha_idxposl_t *cr; ///candidate list
} mg_match_t;

typedef struct {
	uint64_t qse, rse, gld;
} lc_srt_t;

#define lc_srt_key(p) ((p).qse)
KRADIX_SORT_INIT(lc_srt, lc_srt_t, lc_srt_key, member_size(lc_srt_t, qse))

typedef struct {
	uint64_t x, e;
	int32_t d;
	uint32_t id;
} eg_srt_t;

#define eg_srt_x_key(p) ((p).x)
KRADIX_SORT_INIT(eg_srt_x, eg_srt_t, eg_srt_x_key, member_size(eg_srt_t, x))
#define eg_srt_d_key(p) ((p).d)
KRADIX_SORT_INIT(eg_srt_d, eg_srt_t, eg_srt_d_key, member_size(eg_srt_t, d))



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
	int32_t qlen/**, so**/;
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
	uint32_t hash;///hash is path hash, instead of node hash
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

KHASH_MAP_INIT_INT(sp, sp_topk_t)
KHASH_MAP_INIT_INT(sp2, uint64_t)


typedef struct { // data structure for each step in kt_pipeline()
    const mg_idxopt_t *opt;
    const void *ha_flt_tab;
    const ha_pt_t *ha_idx;
    const ma_ug_t *ug;
	const asg_t *rg;
	const ug_opt_t *uopt;
	int n, m, sum_len;
	uint64_t *len, id;
	char **seq;
    ha_mzl_v *mzs;///useless
    st_mt_t *sps;///useless
	mg_gchains_t **gcs;///useless
    mg_tbuf_t **buf;///useless
	ha_ovec_buf_t **hab;
} utepdat_t;

void init_mg_opt(mg_idxopt_t *opt, int is_HPC, int k, int w, int hap_n)
{
    opt->k = k; 
    opt->w = w;
    opt->hap_n = hap_n;
    opt->is_HPC = is_HPC;
    opt->bw = 10000;///2000 in minigraph
    opt->max_gap = 500000;///5000 in minigraph
    opt->occ_weight = 20;
    opt->max_gap_pre = 10000;///1000 in minigraph
    opt->max_lc_iter = 10000;
    opt->chn_pen_gap = 0.19;///using minimap2's value
	opt->max_lc_skip = 25;// mo->max_gc_skip = 25;
	opt->max_lc_iter = 10000;
	opt->min_lc_cnt = 2; 
	opt->min_lc_score = 30;
	opt->max_gc_skip = 25;
	opt->ref_bonus = 0;
	opt->mask_level = 0.5f;
	opt->max_gc_seq_ext = 5;
	opt->seed = 11;
	opt->min_gc_cnt = 3, opt->min_gc_score = 50;
	opt->sub_diff = 6;
	opt->best_n = 5;
	opt->pri_ratio = 0.8f;
}

void uidx_l_build(ma_ug_t *ug, mg_idxopt_t *opt, int cutoff)
{
    ha_flt_tab = ha_ft_ul_gen(&asm_opt, &(ug->u), opt->k, opt->w, cutoff);
    ha_idx = ha_pt_ul_gen(&asm_opt, ha_flt_tab, &(ug->u), opt->k, opt->w, cutoff);	
	fprintf(stderr, "[M::%s] Index has been built.\n", __func__);
}

void uidx_build(ma_ug_t *ug, mg_idxopt_t *opt)
{
    int flag = asm_opt.flag;
    asm_opt.flag |= HA_F_NO_HPC;
    ha_flt_tab = ha_ft_ug_gen(&asm_opt, &(ug->u), opt->is_HPC, opt->k, opt->w, 1, opt->hap_n*5);
    ha_idx = ha_pt_ug_gen(&asm_opt, ha_flt_tab, &(ug->u), opt->is_HPC, opt->k, opt->w, 1);
    asm_opt.flag = flag;
	fprintf(stderr, "[M::%s] Index has been built.\n", __func__);
}

void uidx_destory()
{
    ha_ft_destroy(ha_flt_tab); 
    ha_pt_destroy(ha_idx);
	ha_flt_tab = NULL; ha_idx = NULL;
}

void mg_gres_a_des(mg_gres_a *p)
{
	uint64_t i = 0;
	for (i = 0; i < p->n; i++){
		free(p->a[i].lc); free(p->a[i].gc);
	}
	free(p->a);
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
		if ((tw > max_occ) || (check_unique && tw != 1)) { ///the frequency of repetitive regions; ignore those minimizers
			int en = z->pos + 1, st = en - z->span;//[st, en)
			if (st > rep_en) { ///just record the length of repetive regions
				*rep_len += rep_en - rep_st;
				rep_st = st, rep_en = en;
			} else rep_en = en;
		} else {
			mg_match_t *q = &m[n_m++];
			q->q_pos = z->pos, q->q_span = z->span, q->rev = z->rev, q->cr = cr, q->n = tn, q->qid = 0;
			q->is_tandem = 0, q->weight = 255;
            if(check_unique && tw != 1) q->is_tandem = 1, q->weight = 1;
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
int64_t flt_anchors(int64_t n_a, mg128_t *a, int32_t r)
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


void extend_coordinates(mg_lchain_t *ri, int64_t qlen, int64_t rlen)
{
	int64_t qs, qe, rs, re, qtail, rtail;
	qs = ri->qs; qe = ri->qe - 1; rs = ri->rs; re = ri->re - 1;
	if(ri->v&1) {
		rs = rlen - ri->re; re = rlen - ri->rs - 1;
	}
	
	if(qs <= rs) {
        rs -= qs; qs = 0;
    } else {
    	qs -= rs; rs = 0;
    }

	qtail = qlen - qe - 1; rtail = rlen - re - 1;
    if(qtail <= rtail) {
        qe = qlen - 1; re += qtail;
    }
    else
    {
		re = rlen - 1; qe += rtail; 
    }

	ri->qs = qs; ri->qe = qe + 1; ri->rs = rs; ri->re = re + 1;
	if(ri->v&1) {
		ri->rs = rlen - re - 1; ri->re = rlen - rs;
	}
}
///qlen: query length
///u[]: sc|occ of chains
///a[]: candidate list
mg_lchain_t *mg_lchain_gen(void *km, int qlen, int n_u, uint64_t *u, mg128_t *a, const ma_ug_t *ug)
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
		// fprintf(stderr, "+0+\tA\tutg%.6d%c\t%c\tqs:%u\tqe:%u\tql:%d\tts:%u\tte:%u\ttl:%u\n", 
        //         (ri->v>>1)+1, "lc"[ug->u.a[ri->v>>1].circ], "+-"[ri->v&1], ri->qs, ri->qe, qlen, ri->rs, ri->re, ug->u.a[ri->v>>1].len);
		// extend_coordinates(ri, qlen, ug->u.a[ri->v>>1].len);
		// fprintf(stderr, "-0-\tA\tutg%.6d%c\t%c\tqs:%u\tqe:%u\tql:%d\tts:%u\tte:%u\ttl:%u\n", 
        //         (ri->v>>1)+1, "lc"[ug->u.a[ri->v>>1].circ], "+-"[ri->v&1], ri->qs, ri->qe, qlen, ri->rs, ri->re, ug->u.a[ri->v>>1].len);
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
///the end position of qs is li->qs; dst[]->->qlen indicate the region that need to be checked in bases
mg_pathv_t *mg_shortest_k(void *km0, const asg_t *g, uint32_t src, int32_t n_dst, mg_path_dst_t *dst, int32_t max_dist, int32_t max_k, /** //pathint32_t ql, const char *qs, int is_rev, **/int32_t *n_pathv)
{
	sp_node_t *p, *root = 0, **out;
	sp_topk_t *q;
	khash_t(sp) *h;
	khash_t(sp2) *h2;
	void *km;
	khint_t k;
	int absent;
	int32_t i, j, n_done, n_found;
	uint32_t id, n_out, m_out;
	int8_t *dst_done;
	mg_pathv_t *ret = 0;
	uint64_t *dst_group;
	/** //path
	int32_t n_seeds = 0;
	uint64_t *seeds = 0;
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

	h2 = kh_init2(sp2, km); // (h2+dst_group) keeps all destinations from the same ref id
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

	

	h = kh_init2(sp, km); // h keeps visited vertices; path to each visited vertice
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
	p->hash = __ac_Wang_hash(src);///hash is path hash, instead of node hash
	kavl_insert(sp, &root, p, 0);///should be avl tree
	///each src corresponds to one node in the hash table <h>, but corresponds to <MG_MAX_SHORT_K> node in the AVL tree <root> 
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
		///note that one (sp_node_t->v) might be visited multiple times if there are circles
		///so there might be multipe nodes with the same (sp_node_t->v)
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
							///src is from li from li to lj, so the dis is generally increased; dijkstra algorithm
							///target_dist should be the distance on query
							if (dist == t->target_dist && t->check_hash && r->hash == t->target_hash) done = 1;
							else if ((dist > t->target_dist + MG_SHORT_K_EXT) && (dist > (t->target_dist>>4))) done = 1;
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
			// h keeps visited vertices; path to each visited vertice
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
			} else if ((int32_t)(q->p[0]->di>>32) > d) { // shorter than the longest path so far: replace the longest
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
	/** //path
	mg_idx_hfree(h_seeds);
	kfree(km, seeds);
	kfree(km, mini.a);
	**/
	// NB: AVL nodes are not deallocated. When km==0, they are memory leaks.

	for (i = 0, n_found = 0; i < n_dst; ++i)
		if (dst[i].n_path > 0) ++n_found;///n_path might be larger than 16
	///we can assume n_pathv = NULL for now 
	if (n_found > 0 && n_pathv) { // then generate the backtrack array
		int32_t n, *trans;
		///n_out: how many times that nodes in graph have been visited 
		///note one node might be visited multiples times
		KCALLOC(km, trans, n_out); // used to squeeze unused elements in out[]
		///n_dst: number of alignment chains
		for (i = 0; i < n_dst; ++i) { // mark dst vertices with a target distance
			mg_path_dst_t *t = &dst[i];
			if (t->n_path > 0 && t->target_dist >= 0 && t->path_end >= 0)
				trans[(int32_t)out[t->path_end]->di] = 1;///(int32_t)out[]->di: traverse track corresponds to the alignment chain dst[]
		}
		for (i = 0; (uint32_t)i < n_out; ++i) { // mark dst vertices without a target distance
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
		for (i = n = 0; (uint32_t)i < n_out; ++i) // generate coordinate translations
			if (trans[i]) trans[i] = n++;
			else trans[i] = -1;

		*n_pathv = n;
		KMALLOC(km0, ret, n);
		for (i = 0; (uint32_t)i < n_out; ++i) { // generate the backtrack array
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

static inline int32_t cal_sc(const mg_path_dst_t *dj, const mg_lchain_t *li, const mg_lchain_t *lc, const mg128_t *an, const gc_frag_t *a, const int32_t *f,
							 int bw, int ref_bonus, float chn_pen_gap)
{
	const mg_lchain_t *lj;
	int32_t gap, sc;
	float lin_pen, log_pen;
	if (dj->n_path == 0) return INT32_MIN;
	gap = dj->dist - dj->target_dist;
	lj = &lc[a[dj->meta].i];
	if (gap < 0) gap = -gap;
	if (gap > bw) return INT32_MIN;
	if (lj->qe <= li->qs) sc = li->score;
	else sc = (int32_t)((double)(li->qe - lj->qe) / (li->qe - li->qs) * li->score + .499); // dealing with overlap on query
	//sc += dj->mlen; // TODO: is this line the right thing to do?
	if (dj->is_0) sc += ref_bonus;
	lin_pen = chn_pen_gap * (float)gap;
	log_pen = gap >= 2? mg_log2(gap) : 0.0f;
	sc -= (int32_t)(lin_pen + log_pen);
	sc += f[dj->meta];
	return sc;
}

void transfor_icoord(const int64_t iqs, const int64_t iqe, const int64_t irs, const int64_t ire, const uint8_t rev, 
const int64_t qlen, const int64_t rlen, int32_t *r_qs, int32_t *r_qe, int32_t *r_rs, int32_t *r_re)
{
    int64_t qs, qe, rs, re, qtail, rtail;
    qs = iqs; qe = iqe - 1; rs = irs; re = ire - 1;
    if(rev) {
        rs = rlen - ire; re = rlen - irs - 1;
    }
    
    if(qs <= rs) {
        rs -= qs; qs = 0;
    } else {
        qs -= rs; rs = 0;
    }

    qtail = qlen - qe - 1; rtail = rlen - re - 1;
    if(qtail <= rtail) {
        qe = qlen - 1; re += qtail;
    }
    else
    {
        re = rlen - 1; qe += rtail; 
    }

	if(r_qs) (*r_qs) = qs; if(r_qe) (*r_qe) = qe + 1;
	if(r_rs) (*r_rs) = rs; if(r_re) (*r_re) = re + 1;
    if(rev) {
        if(r_rs) (*r_rs) = rlen - re - 1; 
		if(r_re) (*r_re) = rlen - rs;
    }
}

void transfor_coord(mg_lchain_t *ri, const int64_t qlen, const int64_t rlen, 
int32_t *r_qs, int32_t *r_qe, int32_t *r_rs, int32_t *r_re)
{
    int64_t qs, qe, rs, re, qtail, rtail;
    qs = ri->qs; qe = ri->qe - 1; rs = ri->rs; re = ri->re - 1;
    if(ri->v&1) {
        rs = rlen - ri->re; re = rlen - ri->rs - 1;
    }
    
    if(qs <= rs) {
        rs -= qs; qs = 0;
    } else {
        qs -= rs; rs = 0;
    }

    qtail = qlen - qe - 1; rtail = rlen - re - 1;
    if(qtail <= rtail) {
        qe = qlen - 1; re += qtail;
    }
    else
    {
        re = rlen - 1; qe += rtail; 
    }

	if(r_qs) (*r_qs) = qs; if(r_qe) (*r_qe) = qe + 1;
	if(r_rs) (*r_rs) = rs; if(r_re) (*r_re) = re + 1;
    if(ri->v&1) {
        if(r_rs) (*r_rs) = rlen - re - 1; 
		if(r_re) (*r_re) = rlen - rs;
    }
}

int64_t get_nn_ov(const uint32_t v, const uint32_t w, const asg_t *g)
{
	uint32_t i;
	uint32_t nv = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v), *p = NULL;
	for (i = 0; i < nv; i++) {
		if(av[i].del) continue;
		if(av[i].v == w) {
			// o -= av[i].ol;
			p = &(av[i]);
			break;
		}
	}
	return p?p->ol:0;
}

int64_t get_lchain_ovlp(mg_lchain_t *lp, mg_lchain_t *la, const asg_t *g, const int64_t qlen, const ma_ug_t *ug)
{
	int64_t o = lp->qe - la->qs, oj;
	uint32_t v = la->v^1, w = lp->v^1, i;
	int32_t pqe, aqs;
	if(o <= 0) return 0;
	if(v == w) return o;
	transfor_coord(lp, qlen, ug->u.a[lp->v>>1].len, NULL, &pqe, NULL, NULL);
	transfor_coord(la, qlen, ug->u.a[la->v>>1].len, &aqs, NULL, NULL, NULL);
	uint32_t nv = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v), *p = NULL;
	for (i = 0; i < nv; i++) {
		if(av[i].del) continue;
		if(av[i].v == w) {
			// o -= av[i].ol;
			p = &(av[i]);
			break;
		}
	}
	oj = o;
	if(p) oj = pqe - aqs - p->ol;
	if(o > oj) o = oj;
	if(o < 0) o = 0;
	return o;
}



int64_t get_lchain_gap(mg_lchain_t *lp, mg_lchain_t *la, const asg_t *g, const int64_t qlen, const ma_ug_t *ug, int32_t double_ol)
{
	int64_t gg = la->qs - lp->qe, ggj;
	uint32_t v = la->v^1, w = lp->v^1, i;
	int32_t aqs, pqe;
	if(double_ol == 0 && gg >= 0) return gg;
	if(v == w) return gg;
	transfor_coord(la, qlen, ug->u.a[la->v>>1].len, &aqs, NULL, NULL, NULL);
	transfor_coord(lp, qlen, ug->u.a[lp->v>>1].len, NULL, &pqe, NULL, NULL);
	uint32_t nv = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v), *p = NULL;;
	for (i = 0; i < nv; i++) {
		if(av[i].del) continue;
		if(av[i].v == w) {
			// gg += av[i].ol;
			p = &(av[i]);
			break;
		}
	}
	ggj = gg;
	if(p) ggj = aqs - pqe + p->ol + (double_ol?p->ol:0);
	// if(gg < ggj) gg = ggj;
	// return gg;
	return ggj;
}

int64_t max_ovlp(const asg_t *g, uint32_t v)
{
	uint32_t i, nv = asg_arc_n(g, v), o = 0;
	asg_arc_t *av = asg_arc_a(g, v);
	for (i = 0; i < nv; i++) {
		if(av[i].del) continue;
		if(o < av[i].ol) o = av[i].ol;
	}
	return o;
}

int64_t specific_ovlp(const ma_ug_t *ug, const ug_opt_t *uopt, const uint32_t v, const uint32_t w)
{
	if(ug->u.a[v>>1].circ || ug->u.a[w>>1].circ) return 0;
	uint32_t rv, rw, r = 0, i;
	const ma_hit_t_alloc *x = NULL;
	asg_arc_t t; memset(&t, 0, sizeof(t));
	if(v&1) rv = ug->u.a[v>>1].start^1;
	else rv = ug->u.a[v>>1].end^1;
	if(w&1) rw = ug->u.a[v>>1].end;
	else rw = ug->u.a[v>>1].start;
	x = &(uopt->sources[rv>>1]);
	for (i = 0; i < x->length; i++) {
		if(Get_tn(x->buffer[i]) == (rw>>1)) {
			r = ma_hit2arc(&(x->buffer[i]), uopt->coverage_cut[rv>>1].e - uopt->coverage_cut[rv>>1].s, 
			uopt->coverage_cut[rw>>1].e - uopt->coverage_cut[rw>>1].s, uopt->max_hang, asm_opt.max_hang_rate, uopt->min_ovlp, &t);
			if(r < 0) return 0;
			if((t.ul>>32)!=rv || t.v!=rw) return 0;
			return t.ol;
		}
	}
	return 0;
}

void extend_lchain(mg_lchain_t *lc, int32_t n_lc, int32_t qlen, const ma_ug_t *ug)
{
	int32_t i;
	for (i = 0; i < n_lc; ++i) {
		extend_coordinates(&lc[i], qlen, ug->u.a[lc[i].v>>1].len);
	}
}

void compress_lchain(mg_lchain_t *lc, int32_t n_lc, int32_t qlen, const ma_ug_t *ug, const mg128_t *a)
{
	int32_t i, k, q_span;
	mg_lchain_t *ri = NULL;
	for (i = 0; i < n_lc; ++i) {
		ri = &lc[i];
		k = ri->off;
		ri->rs = (int32_t)a[k].x + 1 > q_span? (int32_t)a[k].x + 1 - q_span : 0; // for HPC k-mer
        ri->qs = (int32_t)a[k].y + 1 - (a[k].y>>32 & 0xff);
        ri->re = (int32_t)a[k + ri->cnt - 1].x + 1;
        ri->qe = (int32_t)a[k + ri->cnt - 1].y + 1;
	}
}

void print_gchain(gc_frag_t *a, const int64_t *p, mg_lchain_t *lc, const int64_t nlc, const ma_ug_t *ug, int32_t qlen)
{
	int64_t k, i;
	gc_frag_t *ai = NULL;
    mg_lchain_t *li = NULL;
	for (k = 0; k < nlc; k++) {
		fprintf(stderr, "\n");
		for (i = k; i >= 0; i = p[i]) {
			ai = &a[i]; li = &lc[ai->i];
			fprintf(stderr, "*\tXXXXXX\tutg%.6d%c\t%c\tqs:%u\tqe:%u\tql:%d\tts:%u\tte:%u\ttl:%u\n", 
					(li->v>>1)+1, "lc"[ug->u.a[li->v>>1].circ], "+-"[li->v&1], li->qs, li->qe, qlen, li->rs, li->re, ug->u.a[li->v>>1].len);
		}
	}
}

// void extend_graph_coordnates(const ma_ug_t *ug, const ug_opt_t *uopt, const mg_lchain_t *lp, const mg_lchain_t *la, 
// mg_coor_t *gp, mg_coor_t *ga, int32_t *go, int32_t *gg)
// {
// 	int32_t so = specific_ovlp(ug, uopt, lp->v^1, la->v^1);
// }

int32_t mg_gchain1_dp(void *km, const ma_ug_t *ug, const asg_t *rg, int32_t *n_lc_, mg_lchain_t *lc, int32_t qlen, int32_t max_dist_g, int32_t max_dist_q, int32_t bw, int32_t max_skip,
                      int32_t ref_bonus, float chn_pen_gap, float mask_level, int32_t max_gc_seq_ext, const ug_opt_t *uopt, const mg128_t *an, uint64_t **u_)
{
    int32_t i, j, k, m_dst, n_dst, n_ext, n_u, n_v, n_lc = *n_lc_, rrs, rre;
    int32_t *f, *v, *t, li_qs, li_qe, li_rs, li_re, lj_qs, lj_qe, lj_rs, lj_re;
    int64_t *p;
    uint64_t *u;
    mg_path_dst_t *dst;
    gc_frag_t *a;
    mg_lchain_t *swap;
    // char *qs;
    asg_t *g = ug->g;

    *u_ = 0;
    if (n_lc == 0) return 0;

	// extend_lchain(lc, n_lc, qlen, ug);
    KMALLOC(km, a, n_lc);
    ///n_lc how many linear chains
    for (i = n_ext = 0; i < n_lc; ++i) { // a[] is a view of frag[]; for sorting
        mg_lchain_t *r = &lc[i];
        gc_frag_t *ai = &a[i];
        int32_t is_isolated = 0, min_end_dist_g;
		transfor_coord(r, qlen, ug->u.a[r->v>>1].len, NULL, NULL, &rrs, &rre);
        r->dist_pre = -1;///indicate parent in graph chain
        min_end_dist_g = g->seq[r->v>>1].len - rre;///r->v: ref_id|rev
        if (rrs < min_end_dist_g) min_end_dist_g = rrs;
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
		// compress_lchain(lc, n_lc, qlen, ug, an);
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
		int32_t mm_ovlp = max_ovlp(ug->g, li->v^1);
		transfor_coord(li, qlen, ug->u.a[li->v>>1].len, &li_qs, &li_qe, &li_rs, &li_re);
        ///note segi is query id, instead of ref id; it is not such useful
        /**
         * a[].x: idx_in_minimizer_arr(32)r_pos(32)
         * a[].y: weight(8)query_id(8)flag(8)span(8)q_pos(32)
        **/
        { // collect end points potentially reachable from _i_
            int32_t x = li->qs + bw + mm_ovlp, n_skip = 0;
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
                int32_t target_dist, dq/**, so = specific_ovlp(ug, uopt, li->v^1, lj->v^1)**/;
				transfor_coord(lj, qlen, ug->u.a[lj->v>>1].len, &lj_qs, &lj_qe, &lj_rs, &lj_re);
                ///lj->qs >= li->qs && lj->qe <= li->qs, so lj is contained
                if (lj->qs >= li->qs) continue; // lj is contained in li on the query coordinate
				/**
				  * doesn't work for overlap graph
                if (lj_qe > li_qs) { // test overlap on the query
                    int o = lj_qe - li_qs - so;///get_lchain_ovlp(lj, li, ug->g, qlen, ug); 
					///mask_level = 0.5, if overlap is too long
                    ///note here is the overlap in query, so too long overlaps might be wrong
                    if (o > (lj->qe - lj->qs) * mask_level || o > (li->qe - li->qs) * mask_level)
                        continue;
                }
				**/
                dq = li_qs - lj_qe;///dq might be smaller than 0
                if (dq > max_dist_q) break; // if query gap too large, stop
                ///The above filter chains like:
                ///1. lj is contained in li
                ///2. the overlap between li and lj is too large
                ///3. li and lj are too far
                ///above we have checked gap/overlap in query 
                ///then we need to check gap/overlap in reference
                if (li->v != lj->v) { // the two linear chains are on two different refs
                    // minimal graph gap; the real graph gap might be larger
                    int32_t min_dist = li_rs + (g->seq[lj->v>>1].len - lj_re); 
                    if (min_dist > max_dist_g) continue; // graph gap too large
                    //note here min_dist - (lj->qs - li->qe) > bw is important
                    //min_dist is always larger than 0, (lj->qs - li->qe) might be negative 
					/**
				  	  * doesn't work for overlap graph
					min_dist -= so;
					if (min_dist - bw > li->qs - lj->qe) continue; ///note seg* is the query id, instead of ref id
					**/
					target_dist = mg_target_dist(g, lj, li);
                    if (target_dist < 0) continue; // this may happen if the query overlap is far too large
                } else if (lj->rs >= li->rs || lj->re >= li->re) { // not colinear
                    continue;
                } else {///li->v == lj->v and colinear; at the same ref id
                    ///w is indel, w is always positive
                    int32_t dr = li->rs - lj->re, dq = li->qs - lj->qe, w = dr > dq? dr - dq : dq - dr;
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
                q->qlen = li->qs - lj->qe;///might be negative
				/**
				  * doesn't work for overlap graph
				q->so = 0;
				if(li->v != lj->v && lj->qe > li->qs) {
					lj_qe = lj->qe; li_qs = li->qs + g->seq[lj->v>>1].len - so;
					q->so = lj_qe - li_qs;
					if(q->so < 0) q->so = 0;
				}
				**/
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
            // (g->seg[li->v>>1].len - li->rs) ----> is like the node length in string graph
            mg_shortest_k(km, g, li->v^1, n_dst, dst, max_dist_g + (g->seq[li->v>>1].len - li->rs), MG_MAX_SHORT_K, /**0, 0, 1,**/ 0);
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
                    if ((int64_t)(dst[j].srt_key - dst[0].srt_key) > li->score)//discard chains with too small weight
                        break;
                n_dst = j;
                if (n_dst > max_gc_seq_ext) n_dst = max_gc_seq_ext; // discard weaker chains
            }
        }
		/** //path
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
        }**/
        { // DP
            int32_t max_f = li->score, max_j = -1, max_d = -1, max_inner = 0;
            uint32_t max_hash = 0;
            for (j = 0; j < n_dst; ++j) {
                mg_path_dst_t *dj = &dst[j];
                int32_t sc;
                sc = cal_sc(dj, li, lc, an, a, f, bw, ref_bonus, chn_pen_gap);
                if (sc == INT32_MIN) continue;
                if (sc > max_f) max_f = sc, max_j = dj->meta, max_d = dj->dist, max_hash = dj->hash, max_inner = dj->inner;
            }
            f[i] = max_f, p[i] = max_j;
            li->dist_pre = max_d;
            li->hash_pre = max_hash;
            li->inner_pre = max_inner;
            v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f;
        }
    }
    kfree(km, dst);

	// print_gchain(a, p, lc, n_ext, ug, qlen);

    // kfree(km, qs);
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
	// compress_lchain(lc, *n_lc_, qlen, ug, an);

    kfree(km, a);
    kfree(km, swap);
    kfree(km, v);
    return n_u;
}


static inline void copy_lchain(mg_llchain_t *q, const mg_lchain_t *p, int32_t *n_a, mg128_t *a_new, const mg128_t *a_old)
{
	q->cnt = p->cnt, q->v = p->v, q->score = p->score;
	memcpy(&a_new[*n_a], &a_old[p->off], q->cnt * sizeof(mg128_t));
	q->off = *n_a;
	(*n_a) += q->cnt;
}

void mg_gchain_extra(const asg_t *g, mg_gchains_t *gs)
{
	int32_t i, j, k;
	for (i = 0; i < gs->n_gc; ++i) { // iterate over gchains
		mg_gchain_t *p = &gs->gc[i];
		const mg_llchain_t *q;
		const mg128_t *last_a;
		int32_t q_span, rest_pl, tmp, n_mini;

		p->qs = p->qe = p->ps = p->pe = -1, p->plen = p->blen = p->mlen = 0, p->div = -1.0f;
		if (p->cnt == 0) continue;
		///some linear chains in middle might be [].cnt == 0
		///but for the first and the last linear chains, [].cnt > 0
		assert(gs->lc[p->off].cnt > 0 && gs->lc[p->off + p->cnt - 1].cnt > 0); // first and last lchains can't be empty
		q = &gs->lc[p->off];
		q_span = (int32_t)(gs->a[q->off].y>>32&0xff);
		/**
		 * a[].x: idx_in_minimizer_arr(32)r_pos(32)
 		 * a[].y: weight(8)query_id(8)flag(8)span(8)q_pos(32)
		 * **/
		p->qs = (int32_t)gs->a[q->off].y + 1 - q_span;///calculated by the first lchain
		p->ps = (int32_t)gs->a[q->off].x + 1 - q_span;///calculated by the first lchain
		tmp = (int32_t)(gs->a[q->off].x>>32);
		assert(p->qs >= 0 && p->ps >= 0);
		q = &gs->lc[p->off + p->cnt - 1];///last lchain
		p->qe = (int32_t)gs->a[q->off + q->cnt - 1].y + 1;
		p->pe = g->seq[q->v>>1].len - (int32_t)gs->a[q->off + q->cnt - 1].x - 1; // this is temporary
		n_mini = (int32_t)(gs->a[q->off + q->cnt - 1].x>>32) - tmp + 1;
		assert(p->n_anchor > 0);

		rest_pl = 0; // this value is never used if the first lchain is not empty (which should always be true)
		last_a = &gs->a[gs->lc[p->off].off];///first minizers in the first linear chain
		for (j = 0; j < p->cnt; ++j) { // iterate over lchains
			const mg_llchain_t *q = &gs->lc[p->off + j];
			int32_t vlen = g->seq[q->v>>1].len;///node length in graph
			p->plen += vlen;
			for (k = 0; k < q->cnt; ++k) { // iterate over anchors
				const mg128_t *r = &gs->a[q->off + k];
				int32_t pl, ql = (int32_t)r->y - (int32_t)last_a->y;
				int32_t span = (int32_t)(r->y>>32&0xff);
				if (j == 0 && k == 0) { // the first anchor on the first lchain
					pl = ql = span;
				} else if (j > 0 && k == 0) { // the first anchor but not on the first lchain
					pl = (int32_t)r->x + 1 + rest_pl;
				} else {
					pl = (int32_t)r->x - (int32_t)last_a->x;
				}
				if (ql < 0) ql = -ql, n_mini += (int32_t)(last_a->x>>32) - (int32_t)(r->x>>32); // dealing with overlapping query at junctions
				p->blen += pl > ql? pl : ql;
				p->mlen += pl > span && ql > span? span : pl < ql? pl : ql;
				last_a = r;
			}
			if (q->cnt == 0) rest_pl += vlen;
			else rest_pl = vlen - (int32_t)gs->a[q->off + q->cnt - 1].x - 1;
		}
		p->pe = p->plen - p->pe;
		assert(p->pe >= p->ps);
		// here n_mini >= p->n_anchor should stand almost all the time
		p->div = n_mini >= p->n_anchor? log((double)n_mini / p->n_anchor) / q_span : log((double)p->n_anchor / n_mini) / q_span;
	}
}

// reorder gcs->a[] and gcs->lc[] such that they are in the same order as gcs->gc[]
void mg_gchain_restore_order(void *km, mg_gchains_t *gcs)
{
	int32_t i, n_a, n_lc;
	mg_llchain_t *lc;
	mg128_t *a;
	KMALLOC(km, lc, gcs->n_lc);
	KMALLOC(km, a, gcs->n_a);
	n_a = n_lc = 0;
	for (i = 0; i < gcs->n_gc; ++i) {
		mg_gchain_t *gc = &gcs->gc[i];
		assert(gc->cnt > 0);
		memcpy(&lc[n_lc], &gcs->lc[gc->off], gc->cnt * sizeof(mg_llchain_t));
		memcpy(&a[n_a], &gcs->a[gcs->lc[gc->off].off], gc->n_anchor * sizeof(mg128_t));
		n_lc += gc->cnt, n_a += gc->n_anchor;
	}
	memcpy(gcs->lc, lc, gcs->n_lc * sizeof(mg_llchain_t));
	memcpy(gcs->a, a, gcs->n_a * sizeof(mg128_t));
	kfree(km, lc); kfree(km, a);
}

// sort chains by score
void mg_gchain_sort_by_score(void *km, mg_gchains_t *gcs)
{
	mg128_t *z;
	mg_gchain_t *gc;
	int32_t i;
	KMALLOC(km, z, gcs->n_gc);
	KMALLOC(km, gc, gcs->n_gc);
	for (i = 0; i < gcs->n_gc; ++i)
		z[i].x = (uint64_t)gcs->gc[i].score << 32 | gcs->gc[i].hash, z[i].y = i;
	radix_sort_128x(z, z + gcs->n_gc);
	for (i = gcs->n_gc - 1; i >= 0; --i)
		gc[gcs->n_gc - 1 - i] = gcs->gc[z[i].y];
	memcpy(gcs->gc, gc, gcs->n_gc * sizeof(mg_gchain_t));
	kfree(km, z); kfree(km, gc);
	mg_gchain_restore_order(km, gcs); // this put gcs in the proper order
}

///u[]: sc|occ of chains
///a[]: candidate list
///gcs[0] = mg_gchain_gen(0, b->km, gi->g, n_gc, u, lc, a, hash, opt->min_gc_cnt, opt->min_gc_score);
// TODO: if frequent malloc() is a concern, filter first and then generate gchains; or generate gchains in thread-local pool and then move to global malloc()
mg_gchains_t *mg_gchain_gen(void *km_dst, void *km, const asg_t *g, int32_t n_u, const uint64_t *u, const mg_lchain_t *lc, const mg128_t *a,
							uint32_t hash, int32_t min_gc_cnt, int32_t min_gc_score)
{
	mg_gchains_t *gc;
	mg_llchain_t *tmp;
	int32_t i, j, k, st, n_g, n_a, s_tmp, n_tmp, m_tmp;
	KCALLOC(km_dst, gc, 1);
	// count the number of gchains and remaining anchors
	// filter out low-quality g_chains
	for (i = 0, st = 0, n_g = n_a = 0; i < n_u; ++i) {
		///nui: how many linear chaisn in i-th g_chain
		int32_t m = 0, nui = (int32_t)u[i];
		for (j = 0; j < nui; ++j) m += lc[st + j].cnt; // m is the number of anchors in this gchain
		if (m >= min_gc_cnt && (int64_t)(u[i]>>32) >= min_gc_score)
			++n_g, n_a += m;
		st += nui;
	}
	if (n_g == 0) return gc;

	// preallocate
	gc->km = km_dst;
	gc->n_gc = n_g, gc->n_a = n_a;
	KCALLOC(km_dst, gc->gc, n_g);///all graph chains
	KMALLOC(km_dst, gc->a, n_a);///all anchors, aka minimizers

	// core loop
	tmp = 0; s_tmp = n_tmp = m_tmp = 0;
	for (i = k = 0, st = 0, n_a = 0; i < n_u; ++i) {
		int32_t n_a0 = n_a, m = 0, nui = (int32_t)u[i]; ///nui: how many linear chaisn in i-th g_chain
		for (j = 0; j < nui; ++j) m += lc[st + j].cnt; ///how many minizers in i-th g_chain
		if (m >= min_gc_cnt && (int64_t)(u[i]>>32) >= min_gc_score) {
			mg_llchain_t *q;
			uint32_t h = hash;

			gc->gc[k].score = u[i]>>32; ///chain score
			gc->gc[k].off = s_tmp; ///all minimizers of k-th chain: gc->a[gc->gc[k].off, )

			for (j = 0; j < nui; ++j) {///how many linear chains
				const mg_lchain_t *p = &lc[st + j];
				h += __ac_Wang_hash(p->qs) + __ac_Wang_hash(p->re) + __ac_Wang_hash(p->v);
			}
			gc->gc[k].hash = __ac_Wang_hash(h);///hash key for the k-th graph chain

			if (n_tmp == m_tmp) KEXPAND(km, tmp, m_tmp);
			// copy the first lchain to gc->a[] and tmp[] (aka, gc->lc[])
			// for the first lchain, it is easy and we just copy all its anchors
			copy_lchain(&tmp[n_tmp++], &lc[st], &n_a, gc->a, a); 
			///0-th lchain has been stored
			///process the remaining chains
			for (j = 1; j < nui; ++j) { 
				const mg_lchain_t *l0 = &lc[st + j - 1], *l1 = &lc[st + j];
				if (!l1->inner_pre) { // bridging two segments; if l0 and l1 are at different reference
					int32_t s, n_pathv;
					mg_path_dst_t dst;
					mg_pathv_t *p;
					memset(&dst, 0, sizeof(mg_path_dst_t));
					dst.v = l0->v ^ 1;
					assert(l1->dist_pre >= 0);
					dst.target_dist = l1->dist_pre;
					dst.target_hash = l1->hash_pre;///hash value of the whole path
					dst.check_hash = 1;
					p = mg_shortest_k(km, g, l1->v^1, 1, &dst, dst.target_dist, MG_MAX_SHORT_K, &n_pathv);
					if (n_pathv == 0 || dst.target_hash != dst.hash)
						fprintf(stderr, "%c[%d] -> %c[%d], dist=%d, target_dist=%d\n", "><"[(l1->v^1)&1], l1->v^1, "><"[(l0->v^1)&1], l0->v^1, dst.dist, dst.target_dist);
					assert(n_pathv > 0);
					assert(dst.target_hash == dst.hash);
					for (s = n_pathv - 2; s >= 1; --s) { // path found in a backward way, so we need to reverse it
						if (n_tmp == m_tmp) KEXPAND(km, tmp, m_tmp);
						q = &tmp[n_tmp++];
						q->off = q->cnt = q->score = 0;
						q->v = p[s].v^1; // when reversing a path, we also need to flip the orientation
					}
					kfree(km, p);
					if (n_tmp == m_tmp) KEXPAND(km, tmp, m_tmp);
					copy_lchain(&tmp[n_tmp++], l1, &n_a, gc->a, a);
				}
				else { // if both of them are at the same linear chain, just merge them
					#if 1
					int32_t k;
					mg_llchain_t *t = &tmp[n_tmp - 1];//the last lchain, have alread done
					assert(l0->v == l1->v);
					// a[].x: ref_id(31)rev(1)r_pos(32)
 					// a[].y: weight(8)query_id(8)flag(8)span(8)q_pos(32)
					for (k = 0; k < l1->cnt; ++k) {
						const mg128_t *ak = &a[l1->off + k];
						if ((int32_t)ak->x > l0->re && (int32_t)ak->y > l0->qe)//find colinear anchors
							break;
					}
					assert(k < l1->cnt);
					t->cnt += l1->cnt - k, t->score += l1->score;
					memcpy(&gc->a[n_a], &a[l1->off + k], (l1->cnt - k) * sizeof(mg128_t));
					n_a += l1->cnt - k;
					#else // don't use this block; for debugging only
					if (n_tmp == m_tmp) KEXPAND(km, tmp, m_tmp);
					copy_lchain(&tmp[n_tmp++], l1, &n_a, gc->a, a);
					#endif
				}
			}
			gc->gc[k].cnt = n_tmp - s_tmp;
			gc->gc[k].n_anchor = n_a - n_a0;
			++k, s_tmp = n_tmp;
		}
		st += nui;//nui: how many linear chains in this gchain
	}
	assert(n_a <= gc->n_a);

	gc->n_a = n_a;
	gc->n_lc = n_tmp;
	KMALLOC(km_dst, gc->lc, n_tmp);
	memcpy(gc->lc, tmp, n_tmp * sizeof(mg_llchain_t));
	kfree(km, tmp);

	mg_gchain_extra(g, gc);
	mg_gchain_sort_by_score(km, gc);
	return gc;
}


// set r[].{id,parent,subsc}, ASSUMING r[] is sorted by score
// mg_gchain_set_parent(b->km, opt->mask_level, gcs[0]->n_gc, gcs[0]->gc, opt->sub_diff, 0);
void mg_gchain_set_parent(void *km, float mask_level, int n, mg_gchain_t *r, int sub_diff, int hard_mask_level)
{
	int i, j, k, *w;
	uint64_t *cov;
	if (n <= 0) return;
	for (i = 0; i < n; ++i) r[i].id = i;
	cov = (uint64_t*)kmalloc(km, n * sizeof(uint64_t));
	w = (int*)kmalloc(km, n * sizeof(int));
	w[0] = 0, r[0].parent = 0;///the first gchain is a primary hits; since all gchains have already been sorted by scores
	for (i = 1, k = 1; i < n; ++i) {///start from the 1-th chain, instead of the 0-th chain
		mg_gchain_t *ri = &r[i];
		int si = ri->qs, ei = ri->qe, n_cov = 0, uncov_len = 0;
		if (hard_mask_level) goto skip_uncov;
		for (j = 0; j < k; ++j) { // traverse existing primary hits to find overlapping hits
			mg_gchain_t *rp = &r[w[j]];
			int sj = rp->qs, ej = rp->qe;
			if (ej <= si || sj >= ei) continue;///no overlaps
			if (sj < si) sj = si;///MAX(si, sj)
			if (ej > ei) ej = ei;///MIN(ei, ej)
			cov[n_cov++] = (uint64_t)sj<<32 | ej;///overlap coordinates
		}
		if (n_cov == 0) {
			goto set_parent_test; // no overlapping primary hits; then i is a new primary hit
		} else if (n_cov > 0) { // there are overlapping primary hits; find the length not covered by existing primary hits
			int j, x = si;
			radix_sort_gfa64(cov, cov + n_cov);
			for (j = 0; j < n_cov; ++j) {
				if ((int)(cov[j]>>32) > x) uncov_len += (cov[j]>>32) - x;
				x = (int32_t)cov[j] > x? (int32_t)cov[j] : x;
			}
			if (ei > x) uncov_len += ei - x;
		}
skip_uncov:
		for (j = 0; j < k; ++j) { // traverse existing primary hits again
			mg_gchain_t *rp = &r[w[j]];
			int sj = rp->qs, ej = rp->qe, min, max, ol;
			if (ej <= si || sj >= ei) continue; // no overlap
			min = ej - sj < ei - si? ej - sj : ei - si;///chain length
			max = ej - sj > ei - si? ej - sj : ei - si;///chain length
			ol = si < sj? (ei < sj? 0 : ei < ej? ei - sj : ej - sj) : (ej < si? 0 : ej < ei? ej - si : ei - si); // overlap length; TODO: this can be simplified
			if ((float)ol / min - (float)uncov_len / max > mask_level) {
				int cnt_sub = 0;
				ri->parent = rp->parent;
				rp->subsc = rp->subsc > ri->score? rp->subsc : ri->score;
				if (ri->cnt >= rp->cnt) cnt_sub = 1;
				if (cnt_sub) ++rp->n_sub;
				break;
			}
		}
set_parent_test:
		if (j == k) w[k++] = i, ri->parent = i, ri->n_sub = 0;
	}
	kfree(km, cov);
	kfree(km, w);
}

// set r[].flt, i.e. mark weak suboptimal chains as filtered
int mg_gchain_flt_sub(float pri_ratio, int min_diff, int best_n, int n, mg_gchain_t *r)
{
	if (pri_ratio > 0.0f && n > 0) {
		int i, k, n_2nd = 0;
		for (i = k = 0; i < n; ++i) {
			int p = r[i].parent;
			if (p == i) { // primary
				r[i].flt = 0, ++k;
			} else if ((r[i].score >= r[p].score * pri_ratio || r[i].score + min_diff >= r[p].score) && n_2nd < best_n) {
				if (!(r[i].qs == r[p].qs && r[i].qe == r[p].qe && r[i].ps == r[p].ps && r[i].pe == r[p].pe)) // not identical hits; TODO: check path as well
					r[i].flt = 0, ++n_2nd, ++k;
				else r[i].flt = 1;
			} else r[i].flt = 1;
		}
		return k;
	}
	return n;
}

// recompute gcs->gc[].{off,n_anchor} and gcs->lc[].off, ASSUMING they are properly ordered (see mg_gchain_restore_order)
void mg_gchain_restore_offset(mg_gchains_t *gcs)
{
	int32_t i, j, n_a, n_lc;
	for (i = 0, n_a = n_lc = 0; i < gcs->n_gc; ++i) {
		mg_gchain_t *gc = &gcs->gc[i];
		gc->off = n_lc;
		for (j = 0, gc->n_anchor = 0; j < gc->cnt; ++j) {
			mg_llchain_t *lc = &gcs->lc[n_lc + j];
			lc->off = n_a;
			n_a += lc->cnt;
			gc->n_anchor += lc->cnt;
		}
		n_lc += gc->cnt;
	}
	assert(n_lc == gcs->n_lc && n_a == gcs->n_a);
}

// hard drop filtered chains, ASSUMING gcs is properly ordered
void mg_gchain_drop_flt(void *km, mg_gchains_t *gcs)
{
	int32_t i, n_gc, n_lc, n_a, n_lc0, n_a0, *o2n;
	if (gcs->n_gc == 0) return;
	KMALLOC(km, o2n, gcs->n_gc);
	for (i = 0, n_gc = 0; i < gcs->n_gc; ++i) {
		mg_gchain_t *r = &gcs->gc[i];
		o2n[i] = -1;
		if (r->flt || r->cnt == 0) continue;
		o2n[i] = n_gc++;
	}
	n_gc = n_lc = n_a = 0;
	n_lc0 = n_a0 = 0;
	for (i = 0; i < gcs->n_gc; ++i) {
		mg_gchain_t *r = &gcs->gc[i];
		if (o2n[i] >= 0) {
			memmove(&gcs->a[n_a], &gcs->a[n_a0], r->n_anchor * sizeof(mg128_t));
			memmove(&gcs->lc[n_lc], &gcs->lc[n_lc0], r->cnt * sizeof(mg_llchain_t));
			gcs->gc[n_gc] = *r;
			gcs->gc[n_gc].id = n_gc;
			gcs->gc[n_gc].parent = o2n[gcs->gc[n_gc].parent];
			++n_gc, n_lc += r->cnt, n_a += r->n_anchor;
		}
		n_lc0 += r->cnt, n_a0 += r->n_anchor;
	}
	assert(n_lc0 == gcs->n_lc && n_a0 == gcs->n_a);
	kfree(km, o2n);
	gcs->n_gc = n_gc, gcs->n_lc = n_lc, gcs->n_a = n_a;
	if (n_a != n_a0) {
		KREALLOC(gcs->km, gcs->a, gcs->n_a);
		KREALLOC(gcs->km, gcs->lc, gcs->n_lc);
		KREALLOC(gcs->km, gcs->gc, gcs->n_gc);
	}
	mg_gchain_restore_offset(gcs);
}

// estimate mapping quality
///mg_gchain_set_mapq(b->km, gcs, qlen, mz->n, opt->min_gc_score);
void mg_gchain_set_mapq(void *km, mg_gchains_t *gcs, int qlen, int max_mini, int min_gc_score)
{
	static const float q_coef = 40.0f;
	int64_t sum_sc = 0;
	float uniq_ratio, r_sc, r_cnt;
	int i, t_sc, t_cnt;
	if (gcs == 0 || gcs->n_gc == 0) return;
	t_sc = qlen < 100? qlen : 100;
	t_cnt = max_mini < 10? max_mini : 10;
	if (t_cnt < 5) t_cnt = 5;
	r_sc = 1.0 / t_sc;
	r_cnt = 1.0 / t_cnt;
	for (i = 0; i < gcs->n_gc; ++i)
		if (gcs->gc[i].parent == gcs->gc[i].id)
			sum_sc += gcs->gc[i].score;///primary chain
	uniq_ratio = (float)sum_sc / (sum_sc + gcs->rep_len);
	for (i = 0; i < gcs->n_gc; ++i) {
		mg_gchain_t *r = &gcs->gc[i];
		if (r->parent == r->id) {///primary chain
			int mapq, subsc;
			float pen_s1 = (r->score > t_sc? 1.0f : r->score * r_sc) * uniq_ratio;
			float x, pen_cm = r->n_anchor > t_cnt? 1.0f : r->n_anchor * r_cnt;
			pen_cm = pen_s1 < pen_cm? pen_s1 : pen_cm;
			subsc = r->subsc > min_gc_score? r->subsc : min_gc_score;
			x = (float)subsc / r->score;
			mapq = (int)(pen_cm * q_coef * (1.0f - x) * logf(r->score));
			mapq -= (int)(4.343f * logf(r->n_sub + 1) + .499f);
			mapq = mapq > 0? mapq : 0;
			if (r->score > subsc && mapq == 0) mapq = 1;
			r->mapq = mapq < 60? mapq : 60;
		} else r->mapq = 0;
	}
}

void mg_map_frag(const void *ha_flt_tab, const ha_pt_t *ha_idx, const ma_ug_t *ug, const asg_t *rg, const uint32_t qid, const int qlen, const char *qseq, ha_mzl_v *mz, 
st_mt_t *sp, mg_tbuf_t *b, int32_t w, int32_t k, int32_t hpc, int32_t mz_sd, int32_t mz_rewin, const mg_idxopt_t *opt, const ug_opt_t *uopt, mg_gchains_t **gcs)
{
    mg128_t *a = NULL;
    int64_t n_a;
    int32_t *mini_pos;
    int i, rep_len, n_mini_pos, n_lc, max_chain_gap_qry, max_chain_gap_ref, n_gc;
	uint32_t hash;
    uint64_t *u;
	mg_lchain_t *lc;
	km_stat_t kmst;
	(*gcs) = NULL;

	hash  = qid;
    hash ^= __ac_Wang_hash(qlen) + __ac_Wang_hash(opt->seed);
    hash  = __ac_Wang_hash(hash);

    mz->n = 0;
    mz2_ha_sketch(qseq, qlen, w, k, 0, hpc, mz, ha_flt_tab, mz_sd, NULL, NULL, NULL, -1, -1, -1, sp, mz_rewin, 1);
	///a[]->y: weight(8)seg_id(8)flag(8)span(8)pos(32);--->query
    ///a[]->x: rid(31)rev(1)rpos(33);--->reference
    a = collect_seed_hits(b->km, opt, 1/**opt->hap_n**/, ha_flt_tab, ha_idx, ug, mz, &n_a, &rep_len, &n_mini_pos, &mini_pos);
	/**
	// might be recover
	if (opt->max_gap_pre > 0 && opt->max_gap_pre * 2 < opt->max_gap) n_a = flt_anchors(n_a, a, opt->max_gap_pre);
    max_chain_gap_qry = max_chain_gap_ref = opt->max_gap;
	**/
	max_chain_gap_qry = max_chain_gap_ref = qlen*2;
    if (n_a == 0) {
        if(a) kfree(b->km, a);
        a = 0, n_lc = 0, u = 0;
    } else {
		a = mg_lchain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_lc_skip, opt->max_lc_iter, 
		opt->min_lc_cnt, opt->min_lc_score, opt->chn_pen_gap, n_a, a, &n_lc, &u, b->km);
    }
	
	if (n_lc) {///n_lc is how many chain we found
		lc = mg_lchain_gen(b->km, qlen, n_lc, u, a, ug);
		for (i = 0; i < n_lc; ++i)///update a[] since ref_id|rev has already been saved to lc[].v
			mg_update_anchors(lc[i].cnt, &a[lc[i].off], n_mini_pos, mini_pos);///update a[].x
	} else lc = 0;
	kfree(b->km, mini_pos); kfree(b->km, u);
	// fprintf(stderr, "++0++qid: %u, qlen: %d, n_a: %ld, n_lc: %d\n", qid, qlen, n_a, n_lc);
	/**
	 * up to here, a[] has been changed
	 * a[].x: idx_in_minimizer_arr(32)r_pos(32)
	 * a[].y: weight(8)query_id(8)flag(8)span(8)q_pos(32)
	**/
	// for (i = 0; i < n_lc; i++) {
	// 	mg_lchain_t *ri = &lc[i];
	// 	fprintf(stderr, "+0)))))))))))))))))))))))))))+\tA\tutg%.6d%c\t%c\tqs:%u\tqe:%u\tql:%d\tts:%u\tte:%u\ttl:%u\n", 
    //             (ri->v>>1)+1, "lc"[ug->u.a[ri->v>>1].circ], "+-"[ri->v&1], ri->qs, ri->qe, qlen, ri->rs, ri->re, ug->u.a[ri->v>>1].len);
	// }
	max_chain_gap_qry = max_chain_gap_ref = opt->max_gap;
	n_gc = mg_gchain1_dp(b->km, ug, rg, &n_lc, lc, qlen, max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_gc_skip, opt->ref_bonus,
						opt->chn_pen_gap, opt->mask_level, opt->max_gc_seq_ext, uopt, a, &u);
	// for (i = 0; i < n_lc; i++) {
	// 	mg_lchain_t *ri = &lc[i];
	// 	fprintf(stderr, "-0-\tA\tutg%.6d%c\t%c\tqs:%u\tqe:%u\tql:%d\tts:%u\tte:%u\ttl:%u\n", 
    //             (ri->v>>1)+1, "lc"[ug->u.a[ri->v>>1].circ], "+-"[ri->v&1], ri->qs, ri->qe, qlen, ri->rs, ri->re, ug->u.a[ri->v>>1].len);
	// }
	(*gcs) = mg_gchain_gen(0, b->km, ug->g, n_gc, u, lc, a, hash, opt->min_gc_cnt, opt->min_gc_score);
	(*gcs)->rep_len = rep_len; (*gcs)->qid = qid; (*gcs)->qlen = qlen;
	kfree(b->km, a);
	kfree(b->km, lc);
	kfree(b->km, u);

	mg_gchain_set_parent(b->km, opt->mask_level, (*gcs)->n_gc, (*gcs)->gc, opt->sub_diff, 0);
	mg_gchain_flt_sub(opt->pri_ratio, k * 2, opt->best_n, (*gcs)->n_gc, (*gcs)->gc);
	mg_gchain_drop_flt(b->km, (*gcs));
	mg_gchain_set_mapq(b->km, (*gcs), qlen, mz->n, opt->min_gc_score);

	if (b->km) {
		km_stat(b->km, &kmst);
		if (kmst.n_blocks != kmst.n_cores) {
			fprintf(stderr, "[E::%s] memory leak at %u\n", __func__, qid);
			abort();
		}
		if (kmst.largest > 1U<<28) {
			km_destroy(b->km);
			b->km = km_init();
		}
	}
	// fprintf(stderr, "++6++qid: %u, (*gcs)->n_gc: %d\n", qid, (*gcs)->n_gc);

}

static void worker_for_ul_alignment(void *data, long i, int tid) // callback for kt_for()
{
    utepdat_t *s = (utepdat_t*)data;
    mg_map_frag(s->ha_flt_tab, s->ha_idx, s->ug, s->rg, s->id+i, s->len[i], s->seq[i], &(s->mzs[tid]), &(s->sps[tid]), s->buf[tid], s->opt->w, s->opt->k, 
    s->opt->is_HPC, asm_opt.mz_sample_dist, asm_opt.mz_rewin, s->opt, s->uopt, &(s->gcs[i]));
}

void dump_gaf(mg_gres_a *hits, const mg_gchains_t *gs, uint32_t only_p)
{
	if (gs == NULL || gs->n_gc == 0 || gs->n_lc == 0) return;
	uint64_t i, j;
	int64_t q_span;
	mg_gres_t *p = NULL;
	kv_pushp(mg_gres_t, *hits, &p); memset(p, 0, sizeof(*p));
	p->n_gc = 0; p->n_lc = 0; p->qid = gs->qid; p->qlen = gs->qlen;
	// p->n_gc = gs->n_gc; p->n_lc = gs->n_lc; p->qid = gs->qid; p->qlen = gs->qlen;
	// MALLOC(p->gc, p->n_gc); memcpy(p->gc, gs->gc, p->n_gc);
	for (i = 0; i < (uint64_t)gs->n_gc; ++i) {
		const mg_gchain_t *t = &gs->gc[i];///one of the gchain
		if(only_p && t->id != t->parent) continue;
		if (t->cnt == 0) continue;
		p->n_gc++; p->n_lc += t->cnt;
	}
	if (p->n_gc == 0) {
		hits->n--;
		return;
	}
	MALLOC(p->gc, p->n_gc); MALLOC(p->lc, p->n_lc); 
	p->n_gc = p->n_lc = 0;
	for (i = 0; i < (uint64_t)gs->n_gc; ++i) {
		const mg_gchain_t *t = &gs->gc[i];///one of the gchain
		if(only_p && t->id != t->parent) continue;
		if (t->cnt == 0) continue;
		p->gc[p->n_gc] = *t; p->gc[p->n_gc].off = p->n_lc;
		for (j = 0; j < (uint64_t)t->cnt; ++j) {
			const mg_llchain_t *q = &gs->lc[t->off + j];
			p->lc[p->n_lc+j].cnt = q->cnt;
			p->lc[p->n_lc+j].score = q->score;
			p->lc[p->n_lc+j].v = q->v;
			if(q->cnt) {
				q_span = (int32_t)(gs->a[q->off].y>>32&0xff);
				p->lc[p->n_lc+j].qs = (int32_t)gs->a[q->off].y + 1 - q_span;///calculated by the first lchain
				p->lc[p->n_lc+j].ts = (int32_t)gs->a[q->off].x + 1 - q_span;///calculated by the first lchain
				p->lc[p->n_lc+j].qe = (int32_t)gs->a[q->off + q->cnt - 1].y + 1;
				p->lc[p->n_lc+j].te = (int32_t)gs->a[q->off + q->cnt - 1].x + 1;
			} else {
				p->lc[p->n_lc+j].qs = p->lc[p->n_lc+j].qe = p->lc[p->n_lc+j].ts = p->lc[p->n_lc+j].te = (uint32_t)-1;
			}
			
			// mg_sprintf_lite(s, "%c%s", "><"[q->v&1], g->seg[q->v>>1].name);
		}
		p->n_gc++; p->n_lc += t->cnt;
	}
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
        s->ha_flt_tab = p->ha_flt_tab; s->ha_idx = p->ha_idx; s->id = p->total_pair; 
		s->opt = p->opt; s->ug = p->ug; s->uopt = p->uopt; s->rg = p->rg;
        while ((ret = kseq_read(p->ks)) >= 0) 
        {
            if (p->ks->seq.l < (uint64_t)p->opt->k) continue;
            if (s->n == s->m) {
                s->m = s->m < 16? 16 : s->m + (s->n>>1);
                REALLOC(s->len, s->m);
                REALLOC(s->seq, s->m);
            }
			/**if(asm_opt.flag & HA_F_VERBOSE_GFA)**/ {
				kv_push(uint64_t, p->nn, p->ks->name.l+p->nn.tl);
				kv_resize(char, p->nn.cc, p->ks->name.l+p->nn.tl);
				memcpy(p->nn.cc.a+p->nn.tl, p->ks->name.s, p->ks->name.l);
				p->nn.tl += p->ks->name.l;
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
		CALLOC(s->gcs, s->n);

		s->buf = (mg_tbuf_t**)calloc(p->n_thread, sizeof(mg_tbuf_t*));
		for (i = 0; i < p->n_thread; ++i) s->buf[i] = mg_tbuf_init();
        
        kt_for(p->n_thread, worker_for_ul_alignment, s, s->n);
        for (i = 0; i < (uint64_t)s->n; ++i) {
            free(s->seq[i]);
            p->total_base += s->len[i];
        }
        free(s->seq); free(s->len); 
		
		for (i = 0; i < p->n_thread; ++i) {
			mg_tbuf_destroy(s->buf[i]);
			free(s->mzs[i].a); free(s->sps[i].a);
		}
		
		free(s->buf); free(s->mzs); free(s->sps);
		return s;
    }
    else if (step == 2) { // step 3: dump
        utepdat_t *s = (utepdat_t*)in;
		uint64_t i;
        for (i = 0; i < (uint64_t)s->n; ++i) {
            // if(s->pos[i].s == (uint64_t)-1) continue;
            // kv_push(pe_hit, p->hits.a, s->pos[i]);
			if(!s->gcs[i]) continue;
			dump_gaf(&(p->hits), s->gcs[i], 1);
			free(s->gcs[i]->gc); free(s->gcs[i]->a); free(s->gcs[i]->lc); free(s->gcs[i]);
        }
        free(s->gcs);
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
	sl->hits.total_base = sl->total_base;
	sl->hits.total_pair = sl->total_pair;
    fprintf(stderr, "[M::%s::%.3f] ==> Qualification\n", __func__, yak_realtime()-index_time);
    return 1;
}

static void *worker_ul_scall_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
    uldat_t *p = (uldat_t*)data;
    ///uint64_t total_base = 0, total_pair = 0;
    if (step == 0) { // step 1: read a block of sequences
        int ret;
        uint64_t l;
		utepdat_t *s;
		CALLOC(s, 1);
        s->ha_flt_tab = p->ha_flt_tab; s->ha_idx = p->ha_idx; s->id = p->total_pair; 
		s->opt = p->opt; s->ug = p->ug; s->uopt = p->uopt; s->rg = p->rg;
        while ((ret = kseq_read(p->ks)) >= 0) 
        {
            if (p->ks->seq.l < (uint64_t)p->opt->k) continue;
            if (s->n == s->m) {
                s->m = s->m < 16? 16 : s->m + (s->n>>1);
                REALLOC(s->len, s->m);
                REALLOC(s->seq, s->m);
            }
			
			append_ul_t(&UL_INF, NULL, p->ks->name.s, p->ks->name.l, NULL, 0, NULL, 0);		
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
        utepdat_t *s = (utepdat_t*)in;

		uint64_t i;
		CALLOC(s->hab, p->n_thread);
		for (i = 0; i < p->n_thread; ++i) s->hab[i] = ha_ovec_init(0, 0, 1);
		///debug
		/**
		uint64_t i;
        CALLOC(s->mzs, p->n_thread);
        CALLOC(s->sps, p->n_thread);
		CALLOC(s->gcs, s->n);

		s->buf = (mg_tbuf_t**)calloc(p->n_thread, sizeof(mg_tbuf_t*));
		for (i = 0; i < p->n_thread; ++i) s->buf[i] = mg_tbuf_init();
        
        kt_for(p->n_thread, worker_for_ul_alignment, s, s->n);
        for (i = 0; i < (uint64_t)s->n; ++i) {
            free(s->seq[i]);
            p->total_base += s->len[i];
        }
        free(s->seq); free(s->len); 
		
		for (i = 0; i < p->n_thread; ++i) {
			mg_tbuf_destroy(s->buf[i]);
			free(s->mzs[i].a); free(s->sps[i].a);
		}
		**/
		for (i = 0; i < p->n_thread; ++i) ha_ovec_destroy(s->hab[i]);
		free(s->hab);
		// free(s->buf); free(s->mzs); free(s->sps);
		return s;
    }
    else if (step == 2) { // step 3: dump
        utepdat_t *s = (utepdat_t*)in;
		uint64_t i, rid;
        for (i = 0; i < (uint64_t)s->n; ++i) {
			///debug
			/**
            // if(s->pos[i].s == (uint64_t)-1) continue;
            // kv_push(pe_hit, p->hits.a, s->pos[i]);
			if(!s->gcs[i]) continue;
			dump_gaf(&(p->hits), s->gcs[i], 1);
			free(s->gcs[i]->gc); free(s->gcs[i]->a); free(s->gcs[i]->lc); free(s->gcs[i]);
			**/
			rid = s->id + i;
			append_ul_t(&UL_INF, &rid, NULL, 0, s->seq[i], s->len[i], NULL, 0);	
			// fprintf(stderr, "%.*s\n", (int)s->len[i], s->seq[i]);	
			free(s->seq[i]); p->total_base += s->len[i];
        }
		///debug
		/**
        free(s->gcs);
		**/
        free(s->len); free(s->seq); free(s);
    }
    return 0;
}

int scall_ul_pipeline(uldat_t* sl, const enzyme *fn)
{
    double index_time = yak_realtime();
    int i;

	init_all_ul_t(&UL_INF, &R_INF);
    for (i = 0; i < fn->n; i++){
        gzFile fp;
        if ((fp = gzopen(fn->a[i], "r")) == 0) return 0;
        sl->ks = kseq_init(fp);
        kt_pipeline(3, worker_ul_scall_pipeline, sl, 3);
        kseq_destroy(sl->ks);
        gzclose(fp);
    }
	sl->hits.total_base = sl->total_base;
	sl->hits.total_pair = sl->total_pair;
    fprintf(stderr, "[M::%s::%.3f] ==> Qualification\n", __func__, yak_realtime()-index_time);
	fprintf(stderr, "[M::%s::%.3f] ==> # reads: %lu, # bases: %lu\n", __func__, yak_realtime()-index_time, 
							UL_INF.n, sl->total_base);
	
    return 1;
}


int print_ul_rs(all_ul_t *U_INF)
{
	uint32_t i;
	ul_vec_t *p = NULL;
	UC_Read ur;
	init_UC_Read(&ur);
	for (i = 0; i < U_INF->n; i++) {
		p = &(U_INF->a[i]);
		retrieve_ul_t(&ur, NULL, U_INF, i, 0, 0, -1);
		fprintf(stderr, ">%s\n", p->n_n);
		fprintf(stderr, "%.*s\n", (int)ur.length, ur.seq);		
	}

	destory_UC_Read(&ur);
	return 1;
}

inline void get_ulname(mg_dbn_t *name, int32_t rid, char **rn, int32_t *rl)
{
	(*rn) = name->cc.a + (rid>0?name->a[rid-1]:0);
	(*rl) = name->a[rid] - (rid>0?name->a[rid-1]:0);
}

void print_gaf(const ma_ug_t *ug, mg_gres_a *hits, mg_dbn_t *name)
{
	uint64_t i, q;
	int32_t k, nl, m;
	char *nn; mg_gchain_t *gc; mg_lres_t *lc; 
	for (i = 0; i < hits->n; i++) {
		q = hits->a[i].qid;
		nn = name->cc.a + (q>0?name->a[q-1]:0);
		nl = name->a[q] - (q>0?name->a[q-1]:0);
		for (k = 0; k < hits->a[i].n_gc; k++) {
			gc = &(hits->a[i].gc[k]);
			fprintf(stderr, "S\t%.*s\tq:id:%lu\tl:n:%d\n", nl, nn, q, gc->cnt);
			for (m = 0; m < gc->cnt; m++) {
				lc = &(hits->a[i].lc[gc->off + m]);
				fprintf(stderr, "*\tA\tutg%.6d%c\t%c\tqs:%u\tqe:%u\tql:%lu\tts:%u\tte:%u\ttl:%u\tcnt:%d\n", 
				(lc->v>>1)+1, "lc"[ug->u.a[lc->v>>1].circ], "+-"[lc->v&1], lc->qs, lc->qe, hits->a[i].qlen, lc->ts, lc->te, ug->u.a[lc->v>>1].len, lc->cnt);
			}
		}
	}
}

void write_ul_hits(mg_gres_a *hits, mg_dbn_t *nn, const char *fn)
{
	char *buf = (char*)calloc(strlen(fn) + 25, 1);
    sprintf(buf, "%s.ul.aln.bin", fn);
    FILE* fp = fopen(buf, "w");
	uint32_t i;

    fwrite(&hits->n, sizeof(hits->n), 1, fp);
	for (i = 0; i < hits->n; i++) {
		fwrite(&hits->a[i].qid, sizeof(hits->a[i].qid), 1, fp);
		fwrite(&hits->a[i].qlen, sizeof(hits->a[i].qlen), 1, fp);
		fwrite(&hits->a[i].n_gc, sizeof(hits->a[i].n_gc), 1, fp);
		fwrite(&hits->a[i].n_lc, sizeof(hits->a[i].n_lc), 1, fp);
		fwrite(hits->a[i].gc, sizeof(mg_gchain_t), hits->a[i].n_gc, fp);
		fwrite(hits->a[i].lc, sizeof(mg_lres_t), hits->a[i].n_lc, fp);
	}
    // fwrite(hits->a, sizeof(mg_gres_t), hits->n, fp);
	fwrite(&hits->total_pair, sizeof(hits->total_pair), 1, fp);
	fwrite(&hits->total_base, sizeof(hits->total_base), 1, fp);

	fwrite(&(nn->n), sizeof(nn->n), 1, fp);
	fwrite(nn->a, sizeof(uint64_t), nn->n, fp);
	fwrite(&(nn->tl), sizeof(nn->tl), 1, fp);
	fwrite(&(nn->cc.n), sizeof(nn->cc.n), 1, fp);
	fwrite(nn->cc.a, sizeof(char), nn->cc.n, fp);
    // write_dbug(ug, fp);

    fclose(fp);
	fprintf(stderr, "[M::%s::] ==> UL alignments have been written\n", __func__);
    free(buf);
}

int load_ul_hits(mg_gres_a *hits, mg_dbn_t *nn, const char *fn)
{
    uint64_t flag = 0;
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
    sprintf(buf, "%s.ul.aln.bin", fn);

    FILE* fp = NULL; 
    fp = fopen(buf, "r"); 
    if(!fp) {
        free(buf);
        return 0;
    }
    uint32_t i;

    kv_init(*hits);
    flag += fread(&hits->n, sizeof(hits->n), 1, fp);
    hits->m = hits->n; MALLOC(hits->a, hits->n);
	for (i = 0; i < hits->n; i++) {
		flag += fread(&hits->a[i].qid, sizeof(hits->a[i].qid), 1, fp);
		flag += fread(&hits->a[i].qlen, sizeof(hits->a[i].qlen), 1, fp);
		flag += fread(&hits->a[i].n_gc, sizeof(hits->a[i].n_gc), 1, fp);
		flag += fread(&hits->a[i].n_lc, sizeof(hits->a[i].n_lc), 1, fp);
		MALLOC(hits->a[i].gc, hits->a[i].n_gc); MALLOC(hits->a[i].lc, hits->a[i].n_lc);
		flag += fread(hits->a[i].gc, sizeof(mg_gchain_t), hits->a[i].n_gc, fp);
		flag += fread(hits->a[i].lc, sizeof(mg_lres_t), hits->a[i].n_lc, fp);
	}
    // flag += fread(hits->a, sizeof(mg_gres_t), hits->n, fp);
	flag += fread(&hits->total_pair, sizeof(hits->total_pair), 1, fp);
	flag += fread(&hits->total_base, sizeof(hits->total_base), 1, fp);

	memset(nn, 0, sizeof(*nn));
	flag += fread(&(nn->n), sizeof(nn->n), 1, fp);
	nn->m = nn->n; MALLOC(nn->a, nn->n);
	flag += fread(nn->a, sizeof(uint64_t), nn->n, fp);
	flag += fread(&(nn->tl), sizeof(nn->tl), 1, fp);
	flag += fread(&(nn->cc.n), sizeof(nn->cc.n), 1, fp);
	nn->cc.m = nn->cc.n; MALLOC(nn->cc.a, nn->cc.n);
	flag += fread(nn->cc.a, sizeof(char), nn->cc.n, fp);

    free(buf);

    // if(!test_dbug(ug, fp))
    // {
    //     free(hits->a.a);
    //     kv_init(hits->a);
    //     fclose(fp);
    //     fprintf(stderr, "[M::%s::] ==> Renew Hi-C linkages\n", __func__);
    //     return 0;
    // }

    fclose(fp);
    fprintf(stderr, "[M::%s::] ==> UL alignments have been loaded\n", __func__);
    return 1;
}

void get_asm_cov(ma_ug_t *ug, uint64_t ul_base, mul_ov_t *aov)
{
	int64_t ss = asm_opt.hg_size;
	if(ss < 0) {
		uint64_t i, k, an;
		int64_t sp;
		asg_t *g = ug->g;
		asg_arc_t *av = NULL;
		for (i = 0, ss = 0; i < g->n_seq; i++) {
			sp = g->seq[i].len; av = asg_arc_a(g, i); an = asg_arc_n(g, i);
			for (k = 0; k < an; k++) {
				if(av[k].del) continue;
				if((av[k].v) < i) {
					sp -= ((int64_t)av[k].ol);
				}
			}
			if(sp > 0) ss += sp;
		}		
	} else {
		ss *= asm_opt.polyploidy;
	}

	if(ss <= 0) ss = 1;
	
	aov->asm_cov = ul_base/ss; aov->asm_size = ss;
	fprintf(stderr, "[M::%s::] ==> asm_cov: %lu, asm_size: %lu\n", __func__, aov->asm_cov, aov->asm_size);
}

int32_t spec_ovlp_occ(eg_srt_t *a, int32_t a_n, int32_t st, int32_t vv, int32_t c_thres)
{
	int32_t i, dst = a[st].d, occ = 1;
	if(occ >= c_thres) return 1;
	for (i = st + 1; i < a_n; i++) {
		if(a[i].id == a[st].id) continue;
		if(a[i].d - dst <= vv) {
			occ++;
			if(occ >= c_thres) return 1; 
		}
	}

	for (i = st - 1; i >= 0; i--) {
		if(a[i].id == a[st].id) continue;
		if(dst - a[i].d <= vv) {
			occ++;
			if(occ >= c_thres) return 1; 
		}
	}	
	return 0;
}


int32_t get_spec_ovlp_occ(eg_srt_t *a, int32_t a_n, int32_t st, int32_t vv, int32_t c_thres, int32_t *s, int32_t *e, kvec_t_u64_warp *res)
{
	int32_t i, dst = a[st].d, occ = 1, pp;
	(*s) = (*e) = st; res->a.n = 0;
	for (i = st + 1; i < a_n; i++) {
		if(a[i].d - dst <= vv) {
			(*e) = i; 
			if(a[i].id == a[st].id) continue;
			occ++; kv_push(uint64_t, res->a, (((uint64_t)(a[i].id))<<32)|i);
		} else {
			break;
		}
	}

	for (i = st - 1; i >= 0; i--) {
		if(dst - a[i].d <= vv) {
			(*s) = i; 
			if(a[i].id == a[st].id) continue;
			occ++; kv_push(uint64_t, res->a, (((uint64_t)(a[i].id))<<32)|i);
		} else {
			break;
		}
	}	
	if(occ >= c_thres) {
		radix_sort_gfa64(res->a.a, res->a.a + res->a.n);
		for (i = 0, pp = -1, occ = 0; i < (int32_t)res->a.n; i++) {
			if((int32_t)(res->a.a[i]>>32) != pp) {
				pp = (res->a.a[i]>>32);
				res->a.a[occ] = res->a.a[i];
				occ++;
			}
		}
		res->a.n = occ;
		if(occ >= c_thres) return occ;
		return 0;
	}
	else {
		return 0;
	}
}

void clean_ul_g(asg_t *xg)
{
	uint32_t n_vtx = xg->n_seq * 2, v, i, nv, ie = 0, ike = 0;
	asg_arc_t *av = NULL;
	uint8_t* bs_flag = NULL; CALLOC(bs_flag, n_vtx);
	buf_t b; memset(&b, 0, sizeof(buf_t)); b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    uint64_t max_dist = get_bub_pop_max_dist_advance(xg, &b);
	for (v = 0; v < xg->n_seq; v++) xg->seq[v].c = 0;
	for (v = 0; v < n_vtx; ++v) {
		if(bs_flag[v] != 0) continue;
		if (asg_arc_n(xg, v) < 2 || xg->seq[v>>1].del) continue;
		if(asg_bub_pop1_primary_trio(xg, NULL, v, max_dist, &b, (uint32_t)-1, (uint32_t)-1, 0, NULL, NULL, NULL, 0, 0, NULL)) {
			//beg is v, end is b.S.a[0]
			//note b.b include end, does not include beg
			for (i = 0; i < b.b.n; i++) {
				if(b.b.a[i]==v || b.b.a[i]==b.S.a[0]) continue;
				bs_flag[b.b.a[i]] = bs_flag[b.b.a[i]^1] = 1;
			}
			bs_flag[v] = 2; bs_flag[b.S.a[0]^1] = 3;
		}
	}

	for (v = 0; v < n_vtx; ++v) {
		if(bs_flag[v] != 0) continue;
		nv = asg_arc_n(xg, v);
		if (nv >= 2) {
			av = asg_arc_a(xg, v);
			for (i = 0; i < nv; ++i){
				if (av[i].ol == 0) {
					av[i].del = 1;
					asg_arc_del(xg, av[i].v^1, (av[i].ul>>32)^1, 1);
					// fprintf(stderr, "---q0-utg%.6d%c, q1-utg%.6d%c\n", 
					// (int32_t)((av[i].ul>>33)+1), "lc"[ug->u.a[av[i].ul>>33].circ], 
					// (int32_t)((av[i].v)>>1)+1, "lc"[ug->u.a[av[i].v].circ]);
				}

				// fprintf(stderr, "xxxx-nv: %u, q0-utg%.6d%c, q1-utg%.6d%c\n", nv, 
				// 	(int32_t)((av[i].ul>>33)+1), "lc"[ug->u.a[av[i].ul>>33].circ], 
				// 	(int32_t)((av[i].v)>>1)+1, "lc"[ug->u.a[av[i].v].circ]);
			}
		}
	}

	for (i = 0; i < xg->n_arc; i++) {
		if(xg->arc[i].ol == 0) {
			ie++;
			if(!xg->arc[i].del) ike++;
		}
	}
	
	fprintf(stderr, "[M::%s::] ==> # fill gaps: %u, # keep gaps: %u\n", __func__, ie, ike);
	free(bs_flag); free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
}

// int32_t max_cluster(int32_t mmi, double vv, int32_t min_off, eg_srt_t *a, int32_t a_n, int32_t st, int32_t st_occ, int32_t *s, int32_t *e, kvec_t_u64_warp *res) 
// {
// 	int32_t i, k, iocc, ovlp;
// 	for (i = st, iocc = 0; i < k; i++) {
// 		ovlp = (a[i].d > mmi? a[i].d - mmi: mmi - a[i].d) * vv;
// 		if(ovlp < min_off) ovlp = min_off;
// 		// fprintf(stderr, "i-%lu, ovlp: %d, td.a[i].d: %d, qid: %u\n", i, ovlp, td.a[i].d, td.a[i].id);
// 		// if(spec_ovlp_occ(td.a + l, k-l, i - l, ovlp, c_thres)) break;
// 		iocc = get_spec_ovlp_occ(td.a + l, k-l, i - l, ovlp, c_thres, &is, &ie, &tidx);
// 		if(iocc >= c_thres) break;
// 	}
// }

void get_ul_g(mul_ov_t *aov, mg_gres_a *hits, ma_ug_t *ug, const asg_t *rg, 
double cov_thres, double vv, int32_t min_off, int32_t min_read_ovlp)
{
	int64_t c_thres = (aov->asm_cov*cov_thres)>2?(aov->asm_cov*cov_thres):2;
	uint64_t i, k, l, m, v0, v1, r0, r1;
	int32_t qs, qe, rs, re, qs0, qe0, qs1, qe1, ovlp, mmi, nngc2 = 0, is, ie, iocc, m_iocc, max_i;
	mg_gres_t *p = NULL;
	mg_gchain_t *gc = NULL, *gc0, *gc1;
	mg_lres_t *lf = NULL, *ll = NULL;
	asg_t *xg = copy_read_graph(ug->g);
	asg_arc_t *pe = NULL;
	kvec_t(lc_srt_t) tt; kv_init(tt); lc_srt_t *pt = NULL;
	kvec_t(eg_srt_t) td; kv_init(td); eg_srt_t *pd = NULL;
	kvec_t_u64_warp tidx; kv_init(tidx.a);
	///for debug
	kvec_t(eg_srt_t) dbg_vw_srt; kv_init(dbg_vw_srt);
	for (i = 0; i < hits->n; i++) {
		// fprintf(stderr, "+i+: %lu\n",i);
		p = &(hits->a[i]); tt.n = 0;
		// fprintf(stderr, "-i-: %lu\n",i);
		if(p->n_gc < 2) continue;
		nngc2++;
		// fprintf(stderr, "\nsis: %lu, p->n_gc: %d\n",i,p->n_gc);
		for (k = 0; k < (uint64_t)p->n_gc; k++) {
			gc = &(p->gc[k]);
			assert(gc->cnt > 0);
			lf = &(p->lc[gc->off]); ll = gc->cnt>1?&(p->lc[gc->off+gc->cnt-1]):NULL;
			assert(lf->qs != (uint32_t)-1);
			if(ll) assert(ll->qs != (uint32_t)-1);

			transfor_icoord(lf->qs, lf->qe, lf->ts, lf->te, lf->v&1, p->qlen, ug->g->seq[lf->v>>1].len, 
			&qs, ll?NULL:&qe, &rs, ll?NULL:&re);
			if(ll) {
				transfor_icoord(ll->qs, ll->qe, ll->ts, ll->te, ll->v&1, p->qlen, ug->g->seq[ll->v>>1].len, 
				NULL, &qe, NULL, &re);
			} else {
				ll = lf;
			}
			if(qe - qs < min_read_ovlp || re - rs < min_read_ovlp) continue;
			kv_pushp(lc_srt_t, tt, &pt); 
			pt->qse = qs; pt->qse <<= 32; pt->qse |= qe;
			pt->rse = rs; pt->rse <<= 32; pt->rse |= re;
			pt->gld = i; pt->gld <<= 32; pt->gld |= k;
			// fprintf(stderr, ">>>>k: %lu, qs: %d, qe: %d, qs-utg%.6d%c, qe-utg%.6d%c\n", k, qs, qe, 
			// (int32_t)((lf->v>>1)+1), "lc"[ug->u.a[lf->v>>1].circ], 
			// (int32_t)((ll->v>>1)+1), "lc"[ug->u.a[ll->v>>1].circ]);
			// fprintf(stderr, "lf_qs: %u, lf_qe: %u, lf_ts: %u, lf_te: %u\n", lf->qs, lf->qe, lf->ts, lf->te);
			// fprintf(stderr, "ll_qs: %u, ll_qe: %u, ll_ts: %u, ll_te: %u\n", ll->qs, ll->qe, ll->ts, ll->te);
		}
		// fprintf(stderr, "eie: %lu\n",i);
		radix_sort_lc_srt(tt.a, tt.a + tt.n);
		for (k = 0; k < tt.n; k++) {
			for (m = k + 1; m < tt.n; m++) {
				gc0 = &(p->gc[(uint32_t)(tt.a[k].gld)]); 
				v0 = p->lc[gc0->off+gc0->cnt-1].v;
				gc1 = &(p->gc[(uint32_t)(tt.a[m].gld)]);
				v1 = p->lc[gc1->off].v;
				if((v0>>1) == (v1>>1)) continue;

				qs0 = tt.a[k].qse>>32; qe0 = (uint32_t)(tt.a[k].qse);
				qs1 = tt.a[m].qse>>32; qe1 = (uint32_t)(tt.a[m].qse);
				// fprintf(stderr, "++++k: %lu, qs0: %d, qe0: %d, qs1: %d, qe1: %d, q0-utg%.6d%c, q1-utg%.6d%c\n", 
				// k, qs0, qe0, qs1, qe1, (int32_t)((v0>>1)+1), "lc"[ug->u.a[v0>>1].circ], (int32_t)((v1>>1)+1), "lc"[ug->u.a[v1>>1].circ]);
				if(qs1 <= qs0 && qe1 >= qe0) continue;///contain
				if(qs0 <= qs1 && qe0 >= qe1) continue;///contain
				if(ug->u.a[v0>>1].circ || ug->u.a[v1>>1].circ) continue;
				ovlp = ((MIN((qe0), (qe1)) > MAX((qs0), (qs1)))? MIN((qe0), (qe1)) - MAX((qs0), (qs1)):0);
				r0 = v0&1?(ug->u.a[v0>>1].start>>1):(ug->u.a[v0>>1].end>>1);
				r1 = v1&1?(ug->u.a[v1>>1].end>>1):(ug->u.a[v1>>1].start>>1);
				// fprintf(stderr, "----k: %lu, ovlp: %d\n", k, ovlp);
				// if((ovlp == 0) || (ovlp <= ((qe0 - qs0)*vv) && ovlp <= ((qe1 - qs1)*vv)) || 
				//  (asg_arc_n(ug->g, v0) == 0 && asg_arc_n(ug->g, v1^1) == 0)) {
				if(/**(asg_arc_n(ug->g, v0) == 0 && asg_arc_n(ug->g, v1^1) == 0) 
												&& **/(ovlp < (int32_t)(MIN(rg->seq[r0].len, rg->seq[r1].len)))) {
					kv_pushp(eg_srt_t, td, &pd);
					pd->d = MAX((qs0), (qs1)) - MIN((qe0), (qe1));
					pd->x = v0<v1?((v0<<32)|v1):(((v1^1)<<32)|(v0^1));
					pd->id = p->qid; 
					pd->e = (uint32_t)(tt.a[k].gld);
					pd->e <<= 32; pd->e |= (uint32_t)(tt.a[m].gld);
				}
			}
		}
	}
	fprintf(stderr, "td.n: %d\n", (int)td.n);
	radix_sort_eg_srt_x(td.a, td.a + td.n);
	for (k = 1, l = 0; k <= td.n; ++k) 
    {   
        if (k == td.n || td.a[k].x != td.a[l].x) 
        {
			if(k - l >= (uint64_t)c_thres) {
				for (i = l+1, mmi = l; i < k; i++) {
					if(td.a[mmi].d > td.a[i].d) mmi = i;
				}
				mmi = td.a[mmi].d < 0? -td.a[mmi].d:0;
				if(mmi != 0) {
					for (i = l; i < k; i++) td.a[i].d += mmi;
				}


				radix_sort_eg_srt_d(td.a + l, td.a + k);
				for (i = l, iocc = 0, tidx.a.n = 0; i < k; i++) {
					ovlp = (td.a[i].d > mmi? td.a[i].d - mmi: mmi - td.a[i].d) * vv;
					if(ovlp < min_off) ovlp = min_off;
					// fprintf(stderr, "i-%lu, ovlp: %d, td.a[i].d: %d, qid: %u\n", i, ovlp, td.a[i].d, td.a[i].id);
					// if(spec_ovlp_occ(td.a + l, k-l, i - l, ovlp, c_thres)) break;
					iocc = get_spec_ovlp_occ(td.a + l, k-l, i - l, ovlp, c_thres, &is, &ie, &tidx);
					// fprintf(stderr, "c_thres-%ld, iocc-%d\n", c_thres, iocc);
					if(iocc >= c_thres) break;
				}

				

				if(i < k) {
					m_iocc = iocc; max_i = i;
					for (i = ie + 1; i < k; i++) {
						iocc = get_spec_ovlp_occ(td.a + l, k-l, i - l, ovlp, m_iocc, &is, &ie, &tidx);
						if(iocc > m_iocc) m_iocc = iocc, max_i = i;
						i = ie + l;
					}
					///for debug
					kv_pushp(eg_srt_t, dbg_vw_srt, &pd);
					pd->x = m_iocc; pd->e = td.a[l].x;
					
					v0 = (uint32_t)td.a[l].x; v1 = td.a[l].x>>32;
					pe = asg_arc_pushp(xg);
					pe->del = 0; pe->strong = 0; pe->el = 0; pe->no_l_indel = 0; pe->ol = 0;
					pe->v = v0; pe->ul = v1<<32; pe->ul += xg->seq[v1>>1].len;

					v0 = (td.a[l].x>>32)^1; v1 = ((uint32_t)td.a[l].x)^1;
					pe = asg_arc_pushp(xg);
					pe->del = 0; pe->strong = 0; pe->el = 0; pe->no_l_indel = 0; pe->ol = 0;
					pe->v = v0; pe->ul = v1<<32; pe->ul += xg->seq[v1>>1].len;

					// fprintf(stderr, "++++q0-utg%.6d%c, q1-utg%.6d%c, k-l: %lu, c_thres: %ld, flag: %u\n", 
					// (int32_t)((td.a[l].x>>33)+1), "lc"[ug->u.a[td.a[l].x>>33].circ], 
					// (int32_t)(((uint32_t)td.a[l].x)>>1)+1, "lc"[ug->u.a[(((uint32_t)td.a[l].x)>>1)].circ], k-l, c_thres,
					// (asg_arc_n(ug->g, ((uint32_t)td.a[l].x)^1) == 0 && asg_arc_n(ug->g, (td.a[l].x>>32)) == 0));
				}
			}
            l = k;
        }
    }

	xg->is_srt = 0; xg->idx = 0; free(xg->idx); 
    asg_cleanup(xg);
	clean_ul_g(xg);

	///for debug
	fprintf(stderr, "[M::%s::] ==> nngc2: %d\n", __func__, nngc2);
	radix_sort_eg_srt_x(dbg_vw_srt.a, dbg_vw_srt.a + dbg_vw_srt.n);
	for (max_i = (int32_t)dbg_vw_srt.n - 1; max_i >= 0; --max_i) {
		pd = &(dbg_vw_srt.a[max_i]);
		fprintf(stderr, "++++q0-utg%.6d%c, q1-utg%.6d%c, occ: %lu, c_thres: %ld, flag: %u\n", 
					(int32_t)((pd->e>>33)+1), "lc"[ug->u.a[pd->e>>33].circ], 
					(int32_t)(((uint32_t)pd->e)>>1)+1, "lc"[ug->u.a[(((uint32_t)pd->e)>>1)].circ], pd->x, c_thres,
					(asg_arc_n(ug->g, ((uint32_t)pd->e)^1) == 0 && asg_arc_n(ug->g, (pd->e>>32)) == 0));
	}

	
	kv_destroy(tt); kv_destroy(td); kv_destroy(tidx.a); kv_destroy(dbg_vw_srt);
	asg_destroy(xg);
} 

int ul_align(mg_idxopt_t *opt, const ug_opt_t *uopt, const asg_t *rg, const enzyme *fn, void *ha_flt_tab, ha_pt_t *ha_idx, ma_ug_t *ug)
{
    uldat_t sl; memset(&sl, 0, sizeof(sl));
    sl.ha_flt_tab = ha_flt_tab;
    sl.ha_idx = ha_idx;
    sl.opt = opt;   
    sl.chunk_size = 200000000;
    sl.n_thread = asm_opt.thread_num;
    sl.ug = ug;
	sl.rg = rg;
	sl.uopt = uopt;
	if(!load_ul_hits(&sl.hits, &sl.nn, asm_opt.output_file_name)) {
    	alignment_ul_pipeline(&sl, fn);
		write_ul_hits(&sl.hits, &sl.nn, asm_opt.output_file_name);
	}

	mul_ov_t aov; memset(&aov, 0, sizeof(aov));
	get_asm_cov(ug, sl.hits.total_base, &aov);
	fprintf(stderr, "[M::%s::] ==> total_pair: %lu, total_base: %lu, n: %d\n", 
										__func__, sl.hits.total_pair, sl.hits.total_base, (int32_t)sl.hits.n);



	get_ul_g(&aov, &sl.hits, ug, rg, 0.51, 0.1, 500, 1000);

	

	// print_gaf(ug, &(sl.hits), &(sl.nn));
	mg_gres_a_des(&(sl.hits)); free(sl.nn.a); free(sl.nn.cc.a);
    return 1;
}

void ul_resolve(ma_ug_t *ug, const asg_t *rg, const ug_opt_t *uopt, int hap_n)
{
	fprintf(stderr, "[M::%s::] ==> UL\n", __func__);
    mg_idxopt_t opt;
    init_mg_opt(&opt, 0, 19, 10, hap_n);
	int exist = (asm_opt.load_index_from_disk? uidx_load(&ha_flt_tab, &ha_idx, asm_opt.output_file_name) : 0);
    if(exist == 0) uidx_build(ug, &opt);
	if(exist == 0) uidx_write(ha_flt_tab, ha_idx, asm_opt.output_file_name);
    ul_align(&opt, uopt, rg, asm_opt.ar, ha_flt_tab, ha_idx, ug);
    uidx_destory();
}

int ul_v_call(mg_idxopt_t *opt, const ug_opt_t *uopt, const enzyme *fn, void *ha_flt_tab, ha_pt_t *ha_idx, ma_ug_t *ug)
{
    uldat_t sl; memset(&sl, 0, sizeof(sl));
    sl.ha_flt_tab = ha_flt_tab;
    sl.ha_idx = ha_idx;
    sl.opt = opt;   
    sl.chunk_size = 100000000;
    sl.n_thread = asm_opt.thread_num;
    sl.ug = ug;
	sl.uopt = uopt;
	scall_ul_pipeline(&sl, fn);
	// print_ul_rs(&UL_INF);
	debug_retrieve_rc_sub(&UL_INF, &R_INF, &(ug->u), 100);
	// if(!load_ul_hits(&sl.hits, &sl.nn, asm_opt.output_file_name)) {
    // 	scall_ul_pipeline(&sl, fn);
	// 	write_ul_hits(&sl.hits, &sl.nn, asm_opt.output_file_name);
	// }

	return 1;
}

ma_ug_t *dedup_HiFis(ma_hit_t_alloc* src, int64_t min_ovlp, int64_t max_hang, int64_t gap_fuzz)
{
	uint64_t i, k, qn, tn, n_read = R_INF.total_reads;
	int32_t r; asg_arc_t t, *p = NULL; 
	asg_t *rg = asg_init();

	rg->m_seq = rg->n_seq = n_read; MALLOC(rg->seq, rg->m_seq);
	for (i = 0; i < n_read; ++i) rg->seq[i].len = Get_READ_LENGTH(R_INF, i);
	
	for (i = 0; i < n_read; i++) {
		if(rg->seq[i].del) continue;
		for (k = 0; k < src[i].length; k++) {
			if(!src[i].buffer[k].el) continue;
			qn = Get_qn(src[i].buffer[k]); tn = Get_tn(src[i].buffer[k]);
			if(rg->seq[qn].del || rg->seq[tn].del) continue;
			if((Get_qe(src[i].buffer[k]) - Get_qs(src[i].buffer[k])) < min_ovlp) continue;
			if((Get_te(src[i].buffer[k]) - Get_ts(src[i].buffer[k])) < min_ovlp) continue;
			r = ma_hit2arc(&(src[i].buffer[k]), rg->seq[qn].len, rg->seq[tn].len, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
			if (r == MA_HT_QCONT) {
				rg->seq[qn].del = 1;
			} else if(r == MA_HT_TCONT) {
				rg->seq[tn].del = 1;
			}
			if(rg->seq[i].del) break;
		}
	}

	for (i = 0; i < n_read; i++) {
		if(rg->seq[i].del) continue;
		for (k = 0; k < src[i].length; k++) {
			if(!src[i].buffer[k].el) continue;
			qn = Get_qn(src[i].buffer[k]); tn = Get_tn(src[i].buffer[k]);
			if(rg->seq[qn].del || rg->seq[tn].del) continue;
			if((Get_qe(src[i].buffer[k]) - Get_qs(src[i].buffer[k])) < min_ovlp) continue;
			if((Get_te(src[i].buffer[k]) - Get_ts(src[i].buffer[k])) < min_ovlp) continue;
			r = ma_hit2arc(&(src[i].buffer[k]), rg->seq[qn].len, rg->seq[tn].len, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
			if (r >= 0) {
                p = asg_arc_pushp(rg);
                *p = t;
            }
		}
	}

	asg_cleanup(rg); asg_symm(rg);
	asg_arc_del_trans(rg, gap_fuzz);
	ma_ug_t *ug = NULL;
	ug = ma_ug_gen(rg);
	asg_destroy(rg);

	for (i = k = 0; i < ug->u.n; i++) k += ug->u.a[i].len;
	fprintf(stderr, "[M::%s::] # unitigs: %lu, # bases: %lu\n", __func__, (uint64_t)ug->u.n, k);
	return ug;
}

void ul_load(const ug_opt_t *uopt)
{
	fprintf(stderr, "[M::%s::] ==> UL\n", __func__);
	mg_idxopt_t opt;
	ma_ug_t *ug = dedup_HiFis(uopt->sources, uopt->min_ovlp, uopt->max_hang, uopt->gap_fuzz);
	int cutoff = asm_opt.hom_cov * asm_opt.high_factor;
	init_aux_table();
	init_mg_opt(&opt, !(asm_opt.flag&HA_F_NO_HPC), 19, 10, cutoff);
	/**
	int exist = (asm_opt.load_index_from_disk? uidx_load(&ha_flt_tab, &ha_idx, asm_opt.output_file_name) : 0);
    if(exist == 0) uidx_l_build(ug, &opt, cutoff);
	if(exist == 0) uidx_write(ha_flt_tab, ha_idx, asm_opt.output_file_name);
	**/
	ul_v_call(&opt, uopt, asm_opt.ar, ha_flt_tab, ha_idx, ug);
	ma_ug_destroy(ug); destory_all_ul_t(&UL_INF);
}