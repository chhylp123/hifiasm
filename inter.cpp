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

void ha_get_ul_candidates_interface(ha_abufl_t *ab, int64_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k, const ul_idx_t *uref, overlap_region_alloc *overlap_list, overlap_region_alloc *overlap_list_hp, Candidates_list *cl, double bw_thres, 
								 int max_n_chain, int keep_whole_chain, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* chain_idx, overlap_region* f_cigar, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t high_occ, void *km);
#define G_CHAIN_BW 16//128
#define FLANK_M (0x7fffU)
#define P_CHAIN_COV 0.985
#define P_FRAGEMENT_CHAIN_COV 0.20
#define P_FRAGEMENT_PRIMARY_CHAIN_COV 0.70
#define P_FRAGEMENT_PRIMARY_SECOND_COV 0.25
#define P_CHAIN_SCORE 0.6
#define G_CHAIN_GAP 0.1
#define UG_SKIP 5
#define RG_SKIP 25
#define G_CHAIN_TRANS_RATE 0.25
#define G_CHAIN_TRANS_WEIGHT -1
#define G_CHAIN_INDEL 128
#define W_CHN_PEN_GAP 0.1
#define N_GCHAIN_RATE 0.04
#define PRIMARY_UL_CHAIN_MIN 75000

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

#define GC_OFFSET_RATE 0.0001
#define GC_OFFSET_POS 8

#define SEC_LEN_DIF 0.03
#define REA_ALIGN_CUTOFF 32
#define CHUNK_SIZE 1000000000

#define generic_key(x) (x)
KRADIX_SORT_INIT(gfa64, uint64_t, generic_key, 8)

#define ul_ov_srt_qe_key(p) ((p).qe)
KRADIX_SORT_INIT(ul_ov_srt_qe, ul_ov_t, ul_ov_srt_qe_key, member_size(ul_ov_t, qe))

#define ul_ov_srt_qs_key(p) ((p).qs)
KRADIX_SORT_INIT(ul_ov_srt_qs, ul_ov_t, ul_ov_srt_qs_key, member_size(ul_ov_t, qs))

#define ul_ov_srt_tn_key(p) ((p).tn)
KRADIX_SORT_INIT(ul_ov_srt_tn, ul_ov_t, ul_ov_srt_tn_key, member_size(ul_ov_t, tn))

#define ul_ov_srt_qn_key(p) ((p).qn)
KRADIX_SORT_INIT(ul_ov_srt_qn, ul_ov_t, ul_ov_srt_qn_key, member_size(ul_ov_t, qn))

#define utg_ct_t_x_key(p) ((p).x)
KRADIX_SORT_INIT(utg_ct_t_x_srt, utg_ct_t, utg_ct_t_x_key, member_size(utg_ct_t, x))

#define utg_ct_t_s_key(p) ((p).s)
KRADIX_SORT_INIT(utg_ct_t_s_srt, utg_ct_t, utg_ct_t_s_key, member_size(utg_ct_t, s))

#define hap_ev_cov_key(x) ((x).cov)
KRADIX_SORT_INIT(hap_ev_cov_srt, haplotype_evdience, hap_ev_cov_key, member_size(haplotype_evdience, cov))

#define uc_block_t_qe_key(x) ((x).qe)
KRADIX_SORT_INIT(uc_block_t_qe_srt, uc_block_t, uc_block_t_qe_key, member_size(uc_block_t, qe));

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
	///base-alignment
	double bw_thres, diff_ec_ul, diff_ec_ul_low; int max_n_chain, ec_ul_round;
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

#define mg_pathv_t_v_srt_key(x) ((x).v)
KRADIX_SORT_INIT(mg_pathv_t_v_srt, mg_pathv_t, mg_pathv_t_v_srt_key, member_size(mg_pathv_t, v))
#define mg_pathv_t_d_srt_key(x) ((x).d)
KRADIX_SORT_INIT(mg_pathv_t_d_srt, mg_pathv_t, mg_pathv_t_d_srt_key, member_size(mg_pathv_t, d))


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

typedef struct {
	FILE *fp;
	ul_vec_t u;
	uint64_t flag;
} ucr_file_t;

typedef struct { // global data structure for kt_pipeline()
    const void *ha_flt_tab;
    const ha_pt_t *ha_idx;
    const mg_idxopt_t *opt;
    const ma_ug_t *ug;
	const asg_t *rg;
	const ug_opt_t *uopt;
	const ul_idx_t *uu;
	ucr_file_t *ucr_s;
	kseq_t *ks;
    int64_t chunk_size;
    uint64_t n_thread;
    uint64_t total_base;
    uint64_t total_pair;
	mg_gres_a hits;
	mg_dbn_t nn;
	uint64_t num_bases, num_corrected_bases, num_recorrected_bases;
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

typedef struct {
	kv_ul_ov_t lo;
	kv_ul_ov_t tk;
	kvec_t_u64_warp srt;
}glchain_t;

typedef struct {
	mg_lchain_t *a;
	size_t n, m;
}vec_mg_lchain_t;

typedef struct {
	mg_path_dst_t *a;
	size_t n, m;
}vec_mg_path_dst_t;

typedef struct {
	sp_node_t **a;
	size_t n, m;
}vec_sp_node_t;
typedef struct {
	mg_pathv_t *a;
	size_t n, m;
}vec_mg_pathv_t;

typedef struct {
	vec_mg_lchain_t l;
	vec_mg_lchain_t swap;
	vec_mg_path_dst_t dst;
	vec_sp_node_t out;
	vec_mg_pathv_t path;
	kvec_t(uint64_t) v;
	kvec_t(int64_t) f;
	st_mt_t dst_done;
}gdpchain_t;

typedef struct { // data structure for each step in kt_pipeline()
    const mg_idxopt_t *opt;
    const void *ha_flt_tab;
    const ha_pt_t *ha_idx;
    const ma_ug_t *ug;
	const asg_t *rg;
	const ug_opt_t *uopt;
	const ul_idx_t *uu;
	int n, m, sum_len;
	uint64_t *len, id;
	char **seq;
    ha_mzl_v *mzs;///useless
    st_mt_t *sps;///useless
	mg_gchains_t **gcs;///useless
    mg_tbuf_t **buf;///useless
	ha_ovec_buf_t **hab;
	glchain_t *ll;
	gdpchain_t *gdp;
	// glchain_t *sec_ll;
	uint64_t num_bases, num_corrected_bases, num_recorrected_bases;
} utepdat_t;


void hc_glchain_destroy(glchain_t *b)
{
	if (!b) return;
	kv_destroy(b->lo); kv_destroy(b->tk); kv_destroy(b->srt.a);
}


void hc_gdpchain_destroy(gdpchain_t *b)
{
	if (!b) return;
	kv_destroy(b->l); kv_destroy(b->swap); kv_destroy(b->dst); kv_destroy(b->out);
	kv_destroy(b->path); kv_destroy(b->v); kv_destroy(b->f); kv_destroy(b->dst_done);
}

void init_mg_opt(mg_idxopt_t *opt, int is_HPC, int k, int w, int hap_n, int max_n_chain, double bw_thres, 
double diff_ec_ul, double diff_ec_ul_low, int ec_ul_round)
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
	opt->max_n_chain = max_n_chain;
	opt->bw_thres = bw_thres;
	opt->diff_ec_ul = diff_ec_ul;
	opt->diff_ec_ul_low = diff_ec_ul_low;
	opt->ec_ul_round = ec_ul_round;
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
///t[]: used for buffer
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
        ///max_f -> score of minimizer
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
			if (p[j] >= 0) t[p[j]] = i;//p[]: prefix idx; means there is a chain longer than 2
		}
		end_j = j;///end_j might be > 0; just the end idx of backwards
		///if not close enough, select a new max
		///max_ii is just used to rescue best-score in case best-score appears before end_j
		if (max_ii < 0 || (int64_t)(a[i].x - a[max_ii].x) > (int64_t)max_dist_x) {///select a new max
			int32_t max = INT32_MIN;
			max_ii = -1;
			for (j = i - 1; j >= st; --j)
				if (max < f[j]) max = f[j], max_ii = j;
		}
		///note: it will happen when `max_ii` < `end_j`; 
		///iteration is terminated at `end_j` mostly because of `max_skip` and `max_iter`
		///max_ii is just used to rescue best-score in case best-score appears before end_j
		if (max_ii >= 0 && max_ii < end_j) {
			int32_t tmp;
			tmp = comput_sc(&a[i], &a[max_ii], max_dist_x, max_dist_y, bw, chn_pen_gap);
			if (tmp != INT32_MIN && max_f < tmp + f[max_ii])
				max_f = tmp + f[max_ii], max_j = max_ii;
		}
		// v[] keeps the peak score up to i (as score might decerase); f[] is the score ending at i, not always the peak
		f[i] = max_f, p[i] = max_j;//p[]: prefix idx
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
	//u[]: sc|occ of chains; chain is mostly sorted by the score; at least the first chain has the largest score
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
	///dst is how many candidates
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
	///for normal graph traversal, one node just has one parental node; here each node has at most 16 parental nodes
	q->k = 1, q->p[0] = p, q->mlen = 0, q->qs = q->qe = -1;

	n_done = 0;
	///the key of avl tree: #define sp_node_cmp(a, b) (((a)->di > (b)->di) - ((a)->di < (b)->di))
	///the higher bits of (*)->di is distance to src node
	///so the key of avl tree is distance
	///in avl tree <root>, one node might be saved multipe times
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
				mg_path_dst_t *t = &dst[(int32_t)dst_group[off + j]];///t is a linear alignment at r->v
				int32_t done = 0;
				///the src and dest are at the same ref id, say we directly find the shortest path
				if (t->inner) {//usually the first node, which is same to src
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
					if (t->n_path == 0) { // means this alignment has never been visited before; keep the shortest path anyway
						copy = 1;
					 // we have a target distance; choose the closest; 
					 // there is already several paths reaching the linear alignment <t>
					} else if (t->target_dist >= 0) { 
						// we found the target path; hash is the path hash including multiple nodes, instead of node hash 
						if (dist == t->target_dist && t->check_hash && r->hash == t->target_hash) { 
							copy = 1, done = 1;
						} else {
							int32_t d0 = t->dist, d1 = dist;
							d0 = d0 > t->target_dist? d0 - t->target_dist : t->target_dist - d0;
							d1 = d1 > t->target_dist? d1 - t->target_dist : t->target_dist - d1;
							///if the new distance (d1) is smaller than the old distance (d0), update the results
							///the length of new path should be closer to t->target_dist
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

int64_t max_ovlp_src(const ug_opt_t *uopt, uint32_t v)
{
	ma_hit_t_alloc* src = uopt->sources;
    int64_t min_ovlp = uopt->min_ovlp, max_hang = uopt->max_hang;
	uint32_t i, qn, tn, o = 0, x = v>>1; asg_arc_t e;

	for (i = 0; i < src[x].length; i++) {
        qn = Get_qn(src[x].buffer[i]); tn = Get_tn(src[x].buffer[i]);
        if(ma_hit2arc(&(src[x].buffer[i]), Get_READ_LENGTH(R_INF, qn), Get_READ_LENGTH(R_INF, tn), 
														max_hang, asm_opt.max_hang_rate, min_ovlp, &e) < 0) {
			continue;
		}
		if((e.ul>>32) != v) continue;
		if(o < e.ol) o = e.ol;
    }

	return o;
}

int64_t specific_ovlp(const ma_ug_t *ug, const ug_opt_t *uopt, const uint32_t v, const uint32_t w)
{
	if(ug->u.a[v>>1].circ || ug->u.a[w>>1].circ) return 0;
	uint32_t rv, rw, i; int32_t r;
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
    ///n_lc how many linear chains; just filter some linear chains
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
    mz2_ha_sketch(qseq, qlen, w, k, 0, hpc, mz, ha_flt_tab, mz_sd, NULL, NULL, NULL, -1, -1, -1, sp, mz_rewin, 1, NULL);
	///a[]->y: weight(8)seg_id(8)flag(8)span(8)pos(32);--->query
    ///a[]->x: rid(31)rev(1)rpos(33);--->reference
    a = collect_seed_hits(b->km, opt, 1/**opt->hap_n**/, ha_flt_tab, ha_idx, ug, mz, &n_a, &rep_len, &n_mini_pos, &mini_pos);
	/**
	// might be recover
	if (opt->max_gap_pre > 0 && opt->max_gap_pre * 2 < opt->max_gap) n_a = flt_anchors(n_a, a, opt->max_gap_pre);
    max_chain_gap_qry = max_chain_gap_ref = opt->max_gap;
	**/
	max_chain_gap_qry = max_chain_gap_ref = qlen*2;
    if (n_a == 0) {//no matched minimizer
        if(a) kfree(b->km, a);
        a = 0, n_lc = 0, u = 0;
    } else {
		a = mg_lchain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_lc_skip, opt->max_lc_iter, 
		opt->min_lc_cnt, opt->min_lc_score, opt->chn_pen_gap, n_a, a, &n_lc, &u, b->km);
    }
	
	if (n_lc) {///n_lc is how many linear chain we found
		lc = mg_lchain_gen(b->km, qlen, n_lc, u, a, ug);//lc->the status of each chain; u->idx of each chain;
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

uint32_t overlap_statistics(overlap_region_alloc* olist, ma_ug_t *ug, int64_t *tt, uint8_t mm)
{
	uint32_t k, sp = (uint32_t)-1, ep = (uint32_t)-1, l = 0;
	for (k = 0; k < olist->length; k++) {
		/**
		if(b->olist.list[k].y_id != 38) continue;
		**/		
		/**
		
		for (z = 0, te = ta = tua = 0; z < b->olist.list[k].w_list_length; z++) {
			if(b->olist.list[k].w_list[z].y_end != -1) {
				te += b->olist.list[k].w_list[z].error;
				ta += b->olist.list[k].w_list[z].x_end + 1 - b->olist.list[k].w_list[z].x_start;
				fprintf(stderr, "x->[%lu, %lu), y->[%d, %d), e->%d\n", b->olist.list[k].w_list[z].x_start, b->olist.list[k].w_list[z].x_end+1, 
					b->olist.list[k].w_list[z].y_start, b->olist.list[k].w_list[z].y_end+1, b->olist.list[k].w_list[z].error);
			}
			else {
				tua += b->olist.list[k].w_list[z].x_end + 1 - b->olist.list[k].w_list[z].x_start;
			}
		}
		fprintf(stderr, "[M::utg%.6d%c::is_match:%u] x->[%u, %u); y->[%u, %u), ualigned->%u, e_rate->%f\n", b->olist.list[k].y_id+1, "lc"[s->ug->u.a[b->olist.list[k].y_id].circ], 
		b->olist.list[k].is_match == 1, b->olist.list[k].x_pos_s, b->olist.list[k].x_pos_e+1, b->olist.list[k].y_pos_s, b->olist.list[k].y_pos_e+1, tua, (float)te/(float)ta);
		**/
		if(tt){
			uint32_t z;
			for (z = 0; z < olist->list[k].w_list_length; z++) {
				if(olist->list[k].w_list[z].y_end != -1) {
					if(tt) *tt += olist->list[k].w_list[z].x_end+1-olist->list[k].w_list[z].x_start;
				}
			}
		}
		if(olist->list[k].is_match == mm) {
			if(sp == (uint32_t)-1 || ep < olist->list[k].x_pos_s) {
				if(sp != (uint32_t)-1) l += ep + 1 - sp;
				sp = olist->list[k].x_pos_s;
				ep = olist->list[k].x_pos_e;
			} else {
				ep = MAX(ep, olist->list[k].x_pos_e);
			}
			if(ug) {
				fprintf(stderr, "[M::utg%.6d%c::is_match->%u] rev->%u, x->[%u, %u), y->[%u, %u)\n", (int)olist->list[k].y_id+1, "lc"[ug->u.a[olist->list[k].y_id].circ], olist->list[k].is_match, 
				olist->list[k].y_pos_strand, olist->list[k].x_pos_s, olist->list[k].x_pos_e+1, olist->list[k].y_pos_s, olist->list[k].y_pos_e+1);
			}
		}
	}
	if(sp != (uint32_t)-1) l += ep + 1 - sp;
	return l;
}
/**
void replace_ul(overlap_region_alloc* olist, Correct_dumy* dumy, haplotype_evdience_alloc* hap, const ul_idx_t *uu)
{
	int64_t k, z, n = 0, c_qs, c_qe, c_ts, c_te, c_rev, p_qs, p_qe, p_te, p_ts, p_rev;
	uint64_t *sc = NULL, *track = NULL;
	overlap_region *c = NULL, *p = NULL;
	dumy->length = 0;
	for (k = 0; k < olist->length; k++) {///has already sorted by x_pos_e
		if(olist->list[k].is_match!=1) continue;
		dumy->overlapID[dumy->length] = (uint64_t)-1;
		dumy->overlapID[dumy->length] <<= 32;
		dumy->overlapID[dumy->length] |= k;
		dumy->length++;
	}

	kv_resize(uint64_t, hap->snp_srt, dumy->length);
	hap->snp_srt.n = dumy->length;
	memset(hap->snp_srt.a, 0, hap->snp_srt.n*sizeof(uint64_t));

	sc = dumy->overlapID; track = hap->snp_srt.a; n = dumy->length;
	for (k = 0; k < n; k++) {
		c = &(olist->list[(uint32_t)track[k]]);
		for (z = k-1; z >= 0; z--) {
			p = &(olist->list[(uint32_t)track[z]]);
		}
	}
}
**/

uint64_t update_ava_het_site(haplotype_evdience_alloc *h, uint64_t oid, uint64_t *beg, uint64_t *end, uint64_t is_srt)
{
	uint64_t k, l, i, occ = 0, n = h->length, need_srt = 0; SnpStats *s = NULL;
	haplotype_evdience tt;
	l = beg? (*beg):0; if(end) (*end) = n; if(beg) (*beg) = n;
	if(l < n && h->list[l].overlapID > oid){
		if(end) (*end) = l;
		return 0;
	}
	for (k = l + 1; k <= n; ++k) {  
		if(h->list[l].overlapID > oid) {
			if(end) (*end) = l;
			break;
		}
        if (k == n || h->list[k].overlapID != h->list[l].overlapID) {
			if(h->list[l].overlapID == oid) {
				for (i = l; i < k; i++) {
					if(h->list[i].type!=1) continue; 
					s = &(h->snp_stat.a[h->list[i].overlapSite]);
					if(s->score == 1 && (!(s->occ_0 < 2 || s->occ_1 < 2))) {
						if(l+occ != i) {
							tt = h->list[l+occ];
							h->list[l+occ] = h->list[i];
							h->list[i] = tt;
						}
						if((occ>0) && (h->list[l+occ].cov<h->list[l+occ-1].cov)) need_srt = 1;
						occ++;
					}
				}
				if(beg) (*beg) = l;
				if(end) (*end) = k;
				break;
			}
            l = k;
        }
    }
	// if(oid == 160) {
	// 	fprintf(stderr, "###[M::%s] l:%lu, occ:%lu\n", __func__, l, occ);
	// 	for (k = l; k < l + occ; k++) {
	// 		fprintf(stderr, "h->list[%lu]:%u\n", k, h->list[k].cov);
	// 	}
	// }
	if(occ && is_srt && need_srt) {
		radix_sort_hap_ev_cov_srt(h->list+l, h->list+l+occ);
	}

	return occ;
}


uint64_t gl_chain_gen(overlap_region_alloc* olist, const ul_idx_t *uref, kv_ul_ov_t *res, uint32_t rec_trans, haplotype_evdience_alloc *hap, void *km)
{
	uint64_t k, o2 = 0, si = 0, ei = 0; ul_ov_t *p = NULL;
	res->n = 0;
	for (k = 0; k < olist->length; k++) {
		if(olist->list[k].is_match==2) o2++;
		if((!rec_trans) && olist->list[k].is_match!=1) continue;
		if(rec_trans && olist->list[k].is_match!=1 && olist->list[k].is_match!=2) continue;
		kv_pushp_km(km, ul_ov_t, *res, &p);
		p->qn = k/**olist->list[k].x_id**/; p->qs = olist->list[k].x_pos_s; p->qe = olist->list[k].x_pos_e+1;
		p->tn = olist->list[k].y_id; p->sec = 0; p->rev = olist->list[k].y_pos_strand;
		p->el = (olist->list[k].is_match==1?1:0);
		if(p->rev) {
			p->ts = uref->ug->u.a[p->tn].len - (olist->list[k].y_pos_e+1); 
			p->te = uref->ug->u.a[p->tn].len - olist->list[k].y_pos_s;
		} else {
			p->ts = olist->list[k].y_pos_s; 
			p->te = olist->list[k].y_pos_e+1;
		}
		if(olist->list[k].is_match==2) {
			p->sec = update_ava_het_site(hap, k, &si, &ei, 1);
			assert(p->sec > 0);
			si = ei;
		}
	}
	return o2;
}

int32_t find_ul_ov_max(int32_t n, const ul_ov_t *a, uint32_t x)
{
    int32_t s = 0, e = n;
    if (n == 0) return -1;
    if (a[n-1].qe < x) return n - 1;
    if (a[0].qe >= x) return -1;
    while (e > s) { // TODO: finish this block
        int32_t m = s + (e - s) / 2;
        if (a[m].qe >= x) e = m;
        else s = m + 1;
    }
    assert(s == e);
    return s;
}



int64_t get_ecov(const ul_idx_t *uref, ul_ov_t *lv, ul_ov_t *lw, int64_t qlen, int64_t bw, double diff_ec_ul)
{
	int64_t dis_q = lv->qe - lw->qe, dis_t = 0, dif, mm;
	uint32_t i, v = ((lv->tn<<1)|lv->rev)^1, w = ((lw->tn<<1)|lw->rev)^1;
	const asg_t *g = uref->ug->g;
	uint32_t nv = asg_arc_n(g, v);
	asg_arc_t *av = asg_arc_a(g, v);
	for (i = 0; i < nv; i++) {
        if(av[i].del || av[i].v != w) continue;
		dis_t = ((uint32_t)av[i].ul);
		dis_t -= (lv->rev?lv->ts:g->seq[v>>1].len-lv->te);
		break;
    }
	
	dif = (dis_q>dis_t? dis_q-dis_t:dis_t-dis_q);
	mm = MAX(dis_q, dis_t); mm *= diff_ec_ul; if(mm < bw) mm = bw;
	// if((v>>1) == 1163 && (w>>1) == 1168) fprintf(stderr, ">>>>>>dis_q:%ld, dis_t:%ld, dif:%ld, mm:%ld\n", dis_q, dis_t, dif, mm);
	if(dif <= mm) return 1;
	return 0;
}

int64_t gl_exact_chain(kv_ul_ov_t *res, kv_ul_ov_t *ex, const ul_idx_t *uref, int64_t bw, double diff_ec_ul, 
int64_t qlen, uint64_t *srt, uint64_t *idx, uint64_t *track, void *km)
{
	// fprintf(stderr, "*****************\n");
	uint32_t li_v, lj_v;
	int64_t mm_ovlp, x, i, j, k, sc, csc, mm_sc, mm_idx;
	ul_ov_t *li = NULL, *lj = NULL;
	const asg_t *g = uref->ug->g;
	radix_sort_ul_ov_srt_qe(res->a, res->a + res->n);
	for (i = 0; i < (int64_t)res->n; ++i) {
		li = &(res->a[i]); li_v = (li->tn<<1)|li->rev;
		mm_ovlp = max_ovlp(g, li_v^1);
		x = (li->qs + mm_ovlp)*diff_ec_ul; 
		if(x < bw) x = bw; 
		x += li->qs + mm_ovlp;
		if (x > qlen+1) x = qlen+1;
		x = find_ul_ov_max(i, res->a, x); 
		csc = retrieve_u_cov_region(uref, li->tn, 0, li->ts, li->te, NULL);
		mm_sc = csc; mm_idx = -1;
		// fprintf(stderr, "---i:%ld, csc:%ld, li->tn:%u, li->ts:%u, li->te:%u\n", i, csc, li->tn, li->ts, li->te);
        for (j = x; j >= 0; --j) { // collect potential destination vertices
			lj = &(res->a[j]); lj_v = (lj->tn<<1)|lj->rev;
			// if(lj->qs >= li->qs) continue; // lj is contained in li on the query coordinate
			if(li_v != lj_v && get_ecov(uref, li, lj, qlen, bw, diff_ec_ul)) {
				sc = csc + (track[j]>>32);
				if(sc > mm_sc) mm_sc = sc, mm_idx = j;
			}
		}
		// 4294967295L
		track[i] = mm_sc; track[i] <<= 32;
		track[i] |= (mm_idx>=0?mm_idx:((uint64_t)0x7FFFFFFF));
		srt[i] = mm_sc; srt[i] <<= 32; srt[i] |= i;
		// fprintf(stderr, "+++i:%ld, mm_idx:%ld, mm_sc:%ld\n", i, mm_idx, mm_sc);
		// fprintf(stderr, "[M::utg%.6d%c] qs->%u; qe->%u\n\n", li->tn+1, "lc"[uref->ug->u.a[li->tn].circ], li->qs, li->qe);
	}
	
	int64_t n_v, n_u, n_v0; 
	radix_sort_gfa64(srt, srt+res->n); ex->n = res->n;
	for (k = (int64_t)res->n-1, n_v = n_u = 0; k >= 0; --k) {
		// fprintf(stderr, "\nk:%ld\n", k);
		n_v0 = n_v;
		for (i = (uint32_t)srt[k]; i >= 0 && (track[i]&((uint64_t)0x80000000)) == 0;) {
			ex->a[n_v++] = res->a[i]; track[i] |= ((uint64_t)0x80000000);
			// fprintf(stderr, "+i:%ld, ", i);
			// fprintf(stderr, "[M::utg%.6d%c] qs->%u; qe->%u\n", res->a[i].tn+1, "lc"[uref->ug->u.a[res->a[i].tn].circ], res->a[i].qs, res->a[i].qe);
			if((track[i]&((uint64_t)0x7FFFFFFF)) == ((uint64_t)0x7FFFFFFF)) i = -1;
			else i = track[i]&((uint64_t)0x7FFFFFFF);
			// if(i>=(int64_t)res->n) fprintf(stderr, "ERROR->i:%ld, res->n:%d, n_v:%ld, qlen:%ld\n", i, (int32_t)res->n, n_v, qlen);
			// fprintf(stderr, "next_i:%ld\n", i);
			// i = (olist->list[i].y_id == (uint32_t)-1?-1:olist->list[i].y_id);
			// fprintf(stderr, "-i:%ld\n", i);
		}
		if(n_v0 == n_v) continue;
		///keep the whole score; do not cut score like minigraph
		// sc = (i<0?(srt[k]>>32):((srt[k]>>32)-olist->list[i].x_id));
		sc = srt[k]>>32;
		idx[n_u++] = ((uint64_t)sc<<32)|(n_v-n_v0);
	}
	// if(n_v != (int64_t)res->n) {
	// 	fprintf(stderr, "\nERROR->n_v:%ld, res->n:%d, qlen:%ld\n", n_v, (int32_t)res->n, qlen);
	// 	for (k = 0; k < (int64_t)res->n; k++) {
	// 		fprintf(stderr, "(%ld)srt-sc:%lu, srt-i:%u\n", k, srt[k]>>32, (uint32_t)srt[k]);
	// 	}

	// 	for (k = 0; k < (int64_t)res->n; k++) {
	// 		fprintf(stderr, "(%ld)track-sc:%lu, track-pi:%lu\n", k, track[k]>>32, track[k]&((uint64_t)0x7FFFFFFF));
	// 	}
	// }
	// if(n_v && ex->a[0].qn == 6) {
	// 	for (k = 0, n_v = n_v0 = 0; k < n_u; k++) {
	// 		n_v0 = n_v; n_v += (uint32_t)idx[k];
	// 		fprintf(stderr, "\n");
	// 		for (i = n_v0; i < n_v; i++) {
	// 			fprintf(stderr, "[%u, %u]\n", ex->a[i].qs, ex->a->qe);
	// 		}
	// 	}
	// }

	for (k = 0, n_v = n_v0 = 0; k < n_u; k++) {
		n_v0 = n_v; n_v += (uint32_t)idx[k];
		res->a[k].qn = idx[k]>>32;
		res->a[k].ts = n_v0; res->a[k].te = n_v;
		res->a[k].qs = ex->a[n_v-1].qs;
		res->a[k].qe = ex->a[n_v0].qe;
	}
	res->n = n_u;
	return res->n;
}

uint64_t get_het_site(haplotype_evdience_alloc *hap, uint32_t oid)
{
	uint64_t k, l, i, occ = 0; SnpStats *s = NULL;
	for (k = 1, l = 0; k <= hap->length; ++k) {  
        if (k == hap->length || hap->list[k].overlapID != hap->list[l].overlapID) {
            if(hap->list[l].overlapID != oid) {
                l = k;
                continue;
            }
            for (i = l; i < k; i++) {
                if(hap->list[i].type!=1) continue; 
                s = &(hap->snp_stat.a[hap->list[i].overlapSite]);
                if(s->score == 1 && (!(s->occ_0 < 2 || s->occ_1 < 2))) {
					occ++;
                }
            }
            l = k;
        }
    }

	return (occ&((uint64_t)0x3FFFFFFF));
}


int64_t get_chain_x(overlap_region* ot, int64_t q)
{
	int64_t x, y, off, i, lx = -1, ly = -1; Fake_Cigar* o = &(ot->f_cigar);
	x = get_fake_gap_pos(o, o->length - 1); 
	off = get_fake_gap_shift(o, o->length - 1);
	y = x - ot->x_pos_s + ot->y_pos_s + off;
	// if(ot->x_id == 98 && (q == 6681 || q == 6990)) fprintf(stderr, "o->length->%u, q->%ld, y->%ld, x->%ld\n", o->length, q, y, x);
	if(y == q) return x;

	for (i = 0; i < (int64_t)o->length; i++){
		x = get_fake_gap_pos(o, i); off = get_fake_gap_shift(o, i);
		y = x - ot->x_pos_s + ot->y_pos_s + off;
		// if(ot->x_id == 98 && (q == 6681 || q == 6990)) fprintf(stderr, "+i->%ld, q->%ld, y->%ld, x->%ld\n", i, q, y, x);
		if(q < y) {
			lx = x; ly = y;
			break;
		}
	}

	if(i == 0 || i == (int64_t)o->length) {
        fprintf(stderr, "ERROR at %s:%d, x_id->%u, y_id->%u, q->%ld, i->%ld, yi_s->%u, yi_e->%u\n", 
		__FILE__, __LINE__, ot->x_id, ot->y_id, q, i, ot->y_pos_s, ot->y_pos_e);
        exit(0);
    }

	x = get_fake_gap_pos(o, i-1); off = get_fake_gap_shift(o, i-1);
	y = x - ot->x_pos_s + ot->y_pos_s + off;
	y = (((double)(q - y))/((double)(ly - y)))*((double)(lx -x)) + x;
	// y = q - y + x;
	if(y < ot->x_pos_s) y = ot->x_pos_s;
	if(y > ot->x_pos_e) y = ot->x_pos_e;
	return y;
}



//[s, e]
double es_win_err(overlap_region* o, int64_t winLen, int64_t s, int64_t e)
{
	int64_t si, ei, os, k, tErr = 0, tLen = 0, minE, maxS, ov;
	os = (o->x_pos_s/winLen)*winLen;
	si = (s-os)/winLen; ei = (e-os)/winLen;
	for (k = si+1; k <= ei-1; k++) {
		tLen += o->w_list[k].x_end+1-o->w_list[k].x_start;
		if(o->w_list[k].y_end != -1) {
			tErr += o->w_list[k].error;
		} else {
			tErr += o->w_list[k].x_end+1-o->w_list[k].x_start;
		}
	}

	k = si;
	maxS = MAX(s, (int64_t)(o->w_list[k].x_start)); minE = MIN(e, (int64_t)(o->w_list[k].x_end)) + 1;
	ov = minE > maxS? minE - maxS:0; 
	if(ov == 0) {
		fprintf(stderr, "WARNNING-1, o->w_list_length->%u, o->x_id->%u, s->%ld, e->%ld, w_list_s->%lu, w_list_e->%lu, winLen->%ld, o->x_pos_s->%u, o->x_pos_e->%u, si->%ld, flag->%d\n", 
		o->w_list_length, o->x_id, s, e, o->w_list[k].x_start, o->w_list[k].x_end, winLen, o->x_pos_s, o->x_pos_e, si, o->w_list[k].y_end);
	}
	tLen += ov/**o->w_list[k].x_end+1-o->w_list[k].x_start**/;
	if(o->w_list[k].y_end != -1) {
		tErr += (ov*o->w_list[k].error)/(o->w_list[k].x_end+1-o->w_list[k].x_start);
	} else {
		tErr += ov/**o->w_list[k].x_end+1-o->w_list[k].x_start**/;
	}

	k = ei;
	maxS = MAX(s, (int64_t)(o->w_list[k].x_start)); minE = MIN(e, (int64_t)(o->w_list[k].x_end)) + 1;
	ov = minE > maxS? minE - maxS:0; 
	if(ov == 0) {
		fprintf(stderr, "WARNNING-2, o->w_list_length->%u, o->x_id->%u, s->%ld, e->%ld, w_list_s->%lu, w_list_e->%lu, winLen->%ld, o->x_pos_s->%u, o->x_pos_e->%u, ei->%ld, flag->%d\n", 
		o->w_list_length, o->x_id, s, e, o->w_list[k].x_start, o->w_list[k].x_end, winLen, o->x_pos_s, o->x_pos_e, ei, o->w_list[k].y_end);
	}
	tLen += ov/**o->w_list[k].x_end+1-o->w_list[k].x_start**/;
	if(o->w_list[k].y_end != -1) {
		tErr += (ov*o->w_list[k].error)/(o->w_list[k].x_end+1-o->w_list[k].x_start);
	} else {
		tErr += ov/**o->w_list[k].x_end+1-o->w_list[k].x_start**/;
	}

	return ((double)tErr)/((double)tLen);
}

int64_t gen_contain_chain(const ul_idx_t *uref, utg_ct_t *p, overlap_region* o, kv_ul_ov_t *chains, double diff_ec_ul, int64_t winLen, void *km)
{
	int64_t y_s, y_e, y_bs, y_be, x_s, x_e, q_s, q_e;
	if(o->y_pos_strand) {
		y_s = uref->ug->u.a[o->y_id].len - p->e;
		y_e = uref->ug->u.a[o->y_id].len - p->s - 1;
	} else {
		y_s = p->s; y_e = p->e - 1;
	}


	y_s = MAX(y_s, (int64_t)o->y_pos_s); y_e = MIN(y_e, (int64_t)o->y_pos_e);
	if(y_s > y_e) return 0;
	x_s = get_chain_x(o, y_s); x_e = get_chain_x(o, y_e) + 1;
	if(x_s >= x_e) fprintf(stderr, "+++y_s->%ld, y_e->%ld, x_s->%ld, x_e->%ld\n", y_s, y_e, x_s, x_e);
	if(o->y_pos_strand) {
		y_bs = uref->ug->u.a[o->y_id].len - (y_e+1);
		y_be = uref->ug->u.a[o->y_id].len - y_s;
	} else {
		y_bs = y_s; y_be = y_e + 1;
	}

	q_s = 0; q_e = p->e - p->s;
	if(p->x&1) {
		q_s += (p->e - y_be);
		q_e -= (y_bs - p->s);
	} else {
		q_s += (y_bs - p->s);
		q_e -= (p->e - y_be); 
	}
	// if(q_s < 0 || q_e < 0 || q_s >= (int64_t)(p->e - p->s) || q_e > (int64_t)(p->e - p->s)) fprintf(stderr, "ERROR\n");
	if(winLen > 0 && diff_ec_ul > 0 && es_win_err(o, winLen, x_s, x_e-1) > diff_ec_ul) return 0;

	ul_ov_t *x = NULL;
	kv_pushp_km(km, ul_ov_t, *chains, &x);
	x->qn = o->x_id; x->qs = x_s; x->qe = x_e;

	/**x->tn = p->x>>1;**/x->tn = (uint32_t)(0x80000000); x->tn |= (p->x>>1);
	x->ts = q_s; x->te = q_e; x->el = 1;x->sec = 0; x->rev = ((o->y_pos_strand == (p->x&1))?0:1);
	// if(x->qn == 0 /**&& ((x->tn<<1)>>1) == 302**/) {
	// 	/**if(o->x_id == 0 && (o->y_id == 46 || o->y_id == 48))**/ {
	// 		// fprintf(stderr, "\nUL[%u]\t%u\t%u\t%c\tUTG[%u]\t%u\t%u\n", o->x_id, o->x_pos_s, o->x_pos_e,
	// 		// 				"+-"[o->y_pos_strand], o->y_id, o->y_pos_s, o->y_pos_e);
	// 		// fprintf(stderr, "Contain[%u]\t%c\ts[%u]\te[%u]\n", p->x>>1, "+-"[p->x&1], p->s, p->e);				
	// 		fprintf(stderr, "U[%u]\t%u\t%u\t%c\tR[%u]\t%u\t%u\tUid[%u]\n", x->qn, x->qs, x->qe, 
	// 											"+-"[x->rev], ((x->tn<<1)>>1), x->ts, x->te, o->y_id);
	// 	}									
	// }
	
	return 1;
}

int64_t debug_utg_ct_t(const ul_idx_t *uref, overlap_region* o, utg_ct_t *ct_a, int64_t ct_n, ma_utg_t *u, utg_ct_t *z, haplotype_evdience *he_a, int64_t he_n)
{
	int64_t k, i, l, rs, re, ss, m = 0;
	utg_ct_t *p = NULL;
	if(ct_a && ct_n) {
		for (i = 0; i < ct_n; i++) {
			p = &(ct_a[i]);
			for (k = 0; k < he_n; k++) {
				ss = o->y_pos_strand?uref->ug->u.a[o->y_id].len - he_a[k].cov - 1:he_a[k].cov;
				if(ss >= p->s && ss < p->e) break;
			}
			if(k < he_n) m++;
		}
	}
	if(u) {
		for (i = l = 0; i < u->n; i++) {
			rs = l; re = l + Get_READ_LENGTH(R_INF, (u->a[i]>>33));
			l += (uint32_t)u->a[i];
			for (k = 0; k < he_n; k++) {
				ss = o->y_pos_strand?uref->ug->u.a[o->y_id].len - he_a[k].cov - 1:he_a[k].cov;
				if(ss >= rs && ss < re) break;
			}
			if(k < he_n) m++;
		}
	}

	if(z) {
		for (k = 0; k < he_n; k++) {
			ss = o->y_pos_strand?uref->ug->u.a[o->y_id].len - he_a[k].cov - 1:he_a[k].cov;
			if(ss >= z->s && ss < z->e) break;
		}
		if(k < he_n) m++;
	}
	return m;
}

int64_t get_het_occ(haplotype_evdience *he_a, int64_t he_n, int64_t c_k, int64_t ylen, utg_ct_t *p, int64_t rev)
{
	int64_t k, occ = 0, ss;
	if(!rev) {
		for (k = c_k; k >= 0; k--) {
			if(he_a[k].cov >= p->s && he_a[k].cov < p->e) {
				occ++;
			} else {
				break;
			}
		}

		for (k = c_k+1; k < he_n; k++) {
			if(he_a[k].cov >= p->s && he_a[k].cov < p->e) {
				occ++;
			} else {
				break;
			}
		}
	} else {
		for (k = c_k; k >= 0; k--) {
			ss = ylen - he_a[k].cov - 1;
			if(ss >= p->s && ss < p->e) {
				occ++;
			} else {
				break;
			}
		}

		for (k = c_k+1; k < he_n; k++) {
			ss = ylen - he_a[k].cov - 1;
			if(ss >= p->s && ss < p->e) {
				occ++;
			} else {
				break;
			}
		}
	}

	assert(occ);
	return occ;
}


int64_t rescue_contain_ul_chains(const ul_idx_t *uref, overlap_region* o, haplotype_evdience *he_a, int64_t he_n, utg_ct_t *ct_a, int64_t ct_n, 
kv_ul_ov_t *chains, double diff_ec_ul, int64_t winLen, int64_t rescue_trans, void *km)
{
	int64_t i, k, ss, ff, t0 = 0;
	uint64_t ys, ye;
	utg_ct_t *p = NULL; 
	// if(o->x_id == 0) {
	// 	fprintf(stderr, "\no->y_id->%u\n", o->y_id);
	// 	for (i = 0; i < ct_n; i++) {
	// 		p = &(ct_a[i]);
	// 		fprintf(stderr, "***rid->%u, rev->%u, s->%u, e->%u\n", p->x>>1, p->x&1, p->s, p->e);
	// 	}
	// }

	if(o->y_pos_strand == 0){
		ys = o->y_pos_s; ye = o->y_pos_e + 1;
		for (i = k = 0; i < ct_n; i++) {
			p = &(ct_a[i]);
			if(p->e <= ys) continue; 
			if(p->s >= ye) break;
			ff = 1;
			if(he_a && he_n > 0) {
				for (; k < he_n; k++) {
					if(he_a[k].cov >= p->s && he_a[k].cov < p->e) {
						ff = 0;
						break;
					}
					if(he_a[k].cov >= p->e) break;
				}
			}
			// if(ff == debug_utg_ct_t(uref, o, p, he_a, he_n)) fprintf(stderr, "ERROR\n");
			if(ff) {
				///push ovlp
				t0 += gen_contain_chain(uref, p, o, chains, diff_ec_ul, winLen, km);
			} else if(rescue_trans) {
                if(gen_contain_chain(uref, p, o, chains, diff_ec_ul, winLen, km)){
                    t0++; chains->a[chains->n-1].el = 0;
					chains->a[chains->n-1].sec = get_het_occ(he_a, he_n, k, uref->ug->u.a[o->y_id].len, p, o->y_pos_strand);
                }
            }
			// if(!ff)	t0++;	
		}
		
	} else {
		ys = uref->ug->u.a[o->y_id].len - (o->y_pos_e+1); 
		ye = uref->ug->u.a[o->y_id].len - o->y_pos_s;
		for (i = 0, k = he_n - 1; i < ct_n; i++) {
			p = &(ct_a[i]); 
			if(p->e <= ys) continue; 
			if(p->s >= ye) break;
			ff = 1;
			if(he_a && he_n > 0) {
				for (; k >= 0; k--) {
					ss = uref->ug->u.a[o->y_id].len - he_a[k].cov - 1;
					if(ss >= p->s && ss < p->e) {
						ff = 0;
						break;
					}
					if(ss >= p->e) break;
				}
			}
			// if(ff == debug_utg_ct_t(uref, o, p, he_a, he_n)) fprintf(stderr, "ERROR\n");
			if(ff) {
				///push ovlp
				t0 += gen_contain_chain(uref, p, o, chains, diff_ec_ul, winLen, km);
			} else if(rescue_trans) {
                if(gen_contain_chain(uref, p, o, chains, diff_ec_ul, winLen, km)){
                    t0++; chains->a[chains->n-1].el = 0;
					chains->a[chains->n-1].sec = get_het_occ(he_a, he_n, k, uref->ug->u.a[o->y_id].len, p, o->y_pos_strand);
                }
            }
			// if(!ff)	t0++;
		}
	}
	// if(debug_utg_ct_t(uref, o, ct_a, ct_n, he_a, he_n)!=t0) fprintf(stderr, "ERROR\n");
	// fprintf(stderr, "***[M::%s] o->y_id:%u, chains->n:%u\n", __func__, o->y_id, (uint32_t)chains->n);
	return t0;
}


int64_t rescue_trans_ul_chains(const ul_idx_t *uref, overlap_region* o, haplotype_evdience *he_a, int64_t he_n, ma_utg_t *u, 
kv_ul_ov_t *chains, double diff_ec_ul, int64_t winLen, int64_t rescue_trans, uint64_t *cis_occ, void *km)
{
	uint64_t ys, ye, i, l;
	int64_t k, ff, ss, t0 = 0;
	utg_ct_t p; 
	if(cis_occ) (*cis_occ) = 0;

	if(o->y_pos_strand == 0) {
		ys = o->y_pos_s; ye = o->y_pos_e + 1;
		for (i = k = l = 0; i < u->n; i++) {
			p.x = u->a[i]>>32; p.s = l; p.e = l + Get_READ_LENGTH(R_INF, (u->a[i]>>33));
			l += (uint32_t)u->a[i];
			if(p.e <= ys) continue; 
			if(p.s >= ye) break;
			ff = 1;
			if(he_a && he_n) {
				for (; k < he_n; k++) {
					if(he_a[k].cov >= p.s && he_a[k].cov < p.e) {
						ff = 0;
						break;
					}
					if(he_a[k].cov >= p.e) break;
				}
			}
			// if(ff == debug_utg_ct_t(uref, o, 0, 0, 0, &p, he_a, he_n)) fprintf(stderr, "ERROR\n");
			if(ff) {
				///push ovlp
				t0 += gen_contain_chain(uref, &p, o, chains, diff_ec_ul, winLen, km);
			} else if(rescue_trans) {
				if(gen_contain_chain(uref, &p, o, chains, diff_ec_ul, winLen, km)){
					t0++; chains->a[chains->n-1].el = 0; if(cis_occ) (*cis_occ)++;
					chains->a[chains->n-1].sec = get_het_occ(he_a, he_n, k, uref->ug->u.a[o->y_id].len, &p, o->y_pos_strand);
				}
			}
			// if(!ff)	t0++;	
		}
	} else {
		ys = uref->ug->u.a[o->y_id].len - (o->y_pos_e+1); 
        ye = uref->ug->u.a[o->y_id].len - o->y_pos_s;
        for (i = l = 0, k = he_n - 1; i < u->n; i++) {
            p.x = u->a[i]>>32; p.s = l; p.e = l + Get_READ_LENGTH(R_INF, (u->a[i]>>33));
			l += (uint32_t)u->a[i];
            if(p.e <= ys) continue; 
            if(p.s >= ye) break;
			ff = 1;
			if(he_a && he_n) {
				for (; k >= 0; k--) {
					ss = uref->ug->u.a[o->y_id].len - he_a[k].cov - 1;
					if(ss >= p.s && ss < p.e) {
						ff = 0;
						break;
					}
					if(ss >= p.e) break;
				}
			}
            // if(ff == debug_utg_ct_t(uref, o, 0, 0, 0, &p, he_a, he_n)) fprintf(stderr, "ERROR\n");
            if(ff) {
                ///push ovlp
                t0 += gen_contain_chain(uref, &p, o, chains, diff_ec_ul, winLen, km);
            } else if(rescue_trans) {
				if(gen_contain_chain(uref, &p, o, chains, diff_ec_ul, winLen, km)){
					t0++; chains->a[chains->n-1].el = 0; if(cis_occ) (*cis_occ)++;
					chains->a[chains->n-1].sec = get_het_occ(he_a, he_n, k, uref->ug->u.a[o->y_id].len, &p, o->y_pos_strand);
				}
			}
            // if(!ff)  t0++;
        }
	}
	// if(debug_utg_ct_t(uref, o, NULL, 0, u, he_a, he_n)!=t0) fprintf(stderr, "ERROR\n");
	// fprintf(stderr, "t0->%ld\n", t0);
    return t0;
}

int64_t dedup_sort_ul_ov_t(ul_ov_t *a, int64_t a_n)
{
	int64_t k, l, z, r, i, qo, to; float rr = 0.9;
	for (k = 1, l = i = 0; k <= a_n; k++) {
		if(k == a_n || a[k].tn != a[l].tn) {///remove the duplicated contained alignments
			for (z = l; z < k; z++) {
				for (r = i-1; r >= 0 && a[r].tn == a[z].tn; r--){
					/**
                	if(a[z].qn == a[r].qn && a[z].qs == a[r].qs && a[z].qe == a[r].qe && 
					   a[z].tn == a[r].tn && a[z].ts == a[r].ts && a[z].te == a[r].te && 
					   a[z].sec == a[r].sec && a[z].el == a[r].el && a[z].rev == a[r].rev) {
						   break;
					}
					**/
					if(a[z].qn == a[r].qn && a[z].tn == a[r].tn && a[z].rev == a[r].rev) {
						qo = ((MIN(a[z].qe, a[r].qe) > MAX(a[z].qs, a[r].qs))? 
															MIN(a[z].qe, a[r].qe) - MAX(a[z].qs, a[r].qs):0);
						to = ((MIN(a[z].te, a[r].te) > MAX(a[z].ts, a[r].ts))?
															MIN(a[z].te, a[r].te) - MAX(a[z].ts, a[r].ts):0);
						if(qo >= ((a[r].qe - a[r].qs)*rr) && qo >= ((a[z].qe - a[z].qs)*rr) && 
						   to >= ((a[r].te - a[r].ts)*rr) && to >= ((a[z].te - a[z].ts)*rr)) {
							break;
						}
					}
				}
				if(r >= 0 && a[r].tn == a[z].tn) {
					if(a[z].el) a[r].el = 1;
					continue;
				}
				a[i++] = a[z];
			}
			l = k;
		}
	}
	return i;
}


uint32_t check_contain_pair(const ug_opt_t *uopt, uint32_t x, uint32_t y, uint32_t check_el)
{
	ma_hit_t_alloc* src = uopt->sources;
	int64_t min_ovlp = uopt->min_ovlp;
	int64_t max_hang = uopt->max_hang;
	uint64_t z, qn, tn; int32_t r = 1; asg_arc_t e;
	for (z = 0; z < src[x].length; z++) {
		if(check_el && (!src[x].buffer[z].el)) continue;
		qn = Get_qn(src[x].buffer[z]); tn = Get_tn(src[x].buffer[z]);
		if(tn != y) continue;
		r = ma_hit2arc(&(src[x].buffer[z]), Get_READ_LENGTH(R_INF, qn), Get_READ_LENGTH(R_INF, tn), max_hang, asm_opt.max_hang_rate, min_ovlp, &e);
		if(r == MA_HT_QCONT || r == MA_HT_TCONT) break;
	}

	if(z < src[x].length) return 1;
	return 0;
}

void debug_contain_ovlps(ul_ov_t *a, uint64_t a_n, const ug_opt_t *uopt)
{
	uint64_t k, i, f; ul_ov_t *z = NULL, *w = NULL;
	for (k = 0; k < a_n; k++) {
		z = &(a[k]);
		if(!(z->tn&((uint32_t)(0x80000000)))) continue;
		for (i = 0, f = z->qn; i < a_n; i++) {
			w = &(a[i]);
			if(i == k) continue;
			if(w->tn&((uint32_t)(0x80000000))) continue;
			if(z->qs >= w->qs && z->qe <= w->qe && check_contain_pair(uopt, (z->tn<<1)>>1, w->tn, 1)) {
				f = (uint32_t)-1;
				break;
			}
		}
		if(z->qn != f) fprintf(stderr, "ERROR\n");
	}
}

ma_hit_t* query_ovlp_src(const ug_opt_t *uopt, uint32_t v, uint32_t w, int64_t o, double diff_ec_ul, uint32_t *ol)
{
	ma_hit_t_alloc* src = uopt->sources;
    int64_t min_ovlp = uopt->min_ovlp;
    int64_t max_hang = uopt->max_hang, d, l, max_l;
    uint64_t z, qn, tn, x = v>>1; int32_t r = 1; asg_arc_t e;
	l = (o>=0?o:-o); //l *= diff_ec_ul;
	if(l <= 0) return NULL;
    for (z = 0; z < src[x].length; z++) {
        qn = Get_qn(src[x].buffer[z]); tn = Get_tn(src[x].buffer[z]);
        if(tn != (w>>1)) continue;
        r = ma_hit2arc(&(src[x].buffer[z]), Get_READ_LENGTH(R_INF, qn), Get_READ_LENGTH(R_INF, tn), max_hang, asm_opt.max_hang_rate, min_ovlp, &e);
        if(r < 0) continue;
		if((e.ul>>32) != v || e.v != w) continue;
		// if(v == 56 && w == 25) fprintf(stderr, "+xxxx, o:%ld, e.ol:%u\n", o, e.ol);
		// if(v == 25 && w == 56) fprintf(stderr, "-xxxx, o:%ld, e.ol:%u\n", o, e.ol);
		d = (o>=e.ol?o-e.ol:e.ol-o);
		max_l = MAX(l, e.ol);
		if(d <= (max_l*diff_ec_ul)) {
			if(ol) (*ol) = e.ol;
			return &(src[x].buffer[z]);
		}
    }
	return NULL;
}

int64_t infer_rovlp(ul_ov_t *li, ul_ov_t *lj, uc_block_t *bi, uc_block_t *bj, All_reads *ridx, ma_ug_t *ug)
{
	int64_t in, is, ie, irev, iqs, iqe, jn, js, je, jrev, jqs, jqe, ir, jr, ts, te, max_s, min_e, s_shift, e_shift;
	
	if(li) {
		in = ug?ug->u.a[li->tn].len:Get_READ_LENGTH(R_INF, li->tn); 
		is = li->ts; ie = li->te; irev = li->rev; iqs = li->qs; iqe = li->qe;
	} else if(bi) {
		in = ug?ug->u.a[bi->hid].len:Get_READ_LENGTH(R_INF, bi->hid); 
		is = bi->ts; ie = bi->te; irev = bi->rev; iqs = bi->qs; iqe = bi->qe;
	} else {
		return 0;
	}

	if(lj) {
		jn = ug?ug->u.a[lj->tn].len:Get_READ_LENGTH(R_INF, lj->tn); 
		js = lj->ts; je = lj->te; jrev = lj->rev; jqs = lj->qs; jqe = lj->qe;
	} else if(bj) {
		jn = ug?ug->u.a[bj->hid].len:Get_READ_LENGTH(R_INF, bj->hid); 
		js = bj->ts; je = bj->te; jrev = bj->rev; jqs = bj->qs; jqe = bj->qe;
	} else {
		return 0;
	}

	max_s = MAX(iqs, jqs); min_e = MIN(iqe, jqe);
	if(min_e <= max_s) return 0;
	s_shift = get_offset_adjust(max_s - iqs, iqe-iqs, ie-is);
	e_shift = get_offset_adjust(iqe - min_e, iqe-iqs, ie-is);
	if(irev) {
		ts = s_shift; s_shift = e_shift; e_shift = ts;
	}
	is += s_shift; ie-= e_shift;

	// if(li && lj && li->tn == 324 && lj->tn == 319 && li->qs == 63841) {
	// 	fprintf(stderr, "+++in:%ld, is:%ld, ie:%ld, irev:%ld, jn:%ld, js:%ld, je:%ld, jrev:%ld\n", in, is, ie, irev, jn, js, je, jrev);
	// }

	s_shift = get_offset_adjust(max_s - jqs, jqe-jqs, je-js);
    e_shift = get_offset_adjust(jqe - min_e, jqe-jqs, je-js);
    if(jrev) {
        ts = s_shift; s_shift = e_shift; e_shift = ts;
    }
    js += s_shift; je-= e_shift;

	if(irev) {
		ts = in - ie; te = in - is;
		is = ts; ie = te;
	}

	if(jrev) {
		ts = jn - je; te = jn - js;
		js = ts; je = te;
	}
	
	// if(li && lj && li->tn == 324 && lj->tn == 319 && li->qs == 63841) {
	// 	fprintf(stderr, "---in:%ld, is:%ld, ie:%ld, irev:%ld, jn:%ld, js:%ld, je:%ld, jrev:%ld\n", in, is, ie, irev, jn, js, je, jrev);
	// }

	if(is <= js) {
		js -= is; is = 0;
	} else {
		is -= js; js = 0;
	}

	ir = in - ie; jr = jn - je;

	if(ir <= jr){
        ie = in; je += ir;        
    }
    else {
		je = jn; ie += jr;
    }

	ir = ie - is; jr = je - js;
	return MAX(ir, jr);
}

void debug_infer_read_ovlp(const ug_opt_t *uopt, double diff_ec_ul, ul_ov_t *li, ul_ov_t *lj, ma_utg_t *u, 
uint32_t i_idx, uint32_t j_idx, All_reads *ridx, ma_ug_t *ug)
{
	uint32_t li_v, lj_v; ma_hit_t *t = NULL;
	li_v = (((uint32_t)(li->tn))<<1)|((uint32_t)(li->rev));
	lj_v = (((uint32_t)(lj->tn))<<1)|((uint32_t)(lj->rev));
	if(lj->qe <= li->qs || li_v == lj_v) fprintf(stderr, "ERROR-1\n");
    t = query_ovlp_src(uopt, li_v^1, lj_v^1, infer_rovlp(li, lj, NULL, NULL, ridx, ug), diff_ec_ul, NULL);
	// ((int64_t)(lj->qe))-((int64_t)(li->qs))
	if(!t /**&& (li_v^1) == 648 && (lj_v^1) == 638 && li->qs == 63841**/) {
		fprintf(stderr, "ERROR-2, li_v^1->%u, li->qs->%u, li->qe->%u, lj_v^1->%u, lj->qs->%u, lj->qe->%u, infer_rovlp->%ld\n", 
		li_v^1, li->qs, li->qe, lj_v^1, lj->qs, lj->qe, infer_rovlp(li, lj, NULL, NULL, ridx, ug));
	}
}

uint64_t infer_read_ovlp(const ul_idx_t *uref, overlap_region_alloc* olist, kv_ul_ov_t *in, kv_ul_ov_t *res, double diff_ec_ul, int64_t winLen, const ug_opt_t *uopt, ul_contain *ct, void *km)
{
	uint64_t t, k, l, t_0, pb, cut = res->n, c_occ = 0;; 
	ma_ug_t *ug = uref->ug;
	ma_utg_t *u = NULL; 
	overlap_region* o = NULL;
	ul_ov_t *z = NULL;
	utg_ct_t p;
	// res->n = 0;
	for (t = 0; t < in->n; t++) {
		if(!(in->a[t].tn&(uint32_t)(0x80000000))) {///uid
			u = &(ug->u.a[in->a[t].tn]); o = &(olist->list[in->a[t].qn]);
			assert(o->y_id == in->a[t].tn);
			for (k = l = 0, pb = res->n+2; k < u->n; k++) {
				p.x = u->a[k]>>32; 
				p.s = l; p.e = l + Get_READ_LENGTH(R_INF, (u->a[k]>>33));
				l += (uint32_t)u->a[k];
				if(p.e <= in->a[t].ts) continue;
				if(p.s >= in->a[t].te) break;
				t_0 = gen_contain_chain(uref, &p, o, res, -1, -1, km);	
				assert(t_0 > 0);
				// if(t_0 == 0) {
				// 	fprintf(stderr, "ERROR-2, o->x_id:%u, o->y_id:%u, k:%lu, u->n:%lu, p.s:%u, p.e:%u, ts:%u, te:%u, rev:%u\n", 
				// 	o->x_id, o->y_id, k, (uint64_t)u->n, p.s, p.e, in->a[t].ts, in->a[t].te, in->a[t].rev);
				// }
				res->a[res->n-1].el = in->a[t].el;
				res->a[res->n-1].sec = in->a[t].sec;
				res->a[res->n-1].tn <<= 1; 
				res->a[res->n-1].tn >>= 1;
				res->a[res->n-1].qn = o->x_id;
				if(res->n >= pb) {
					if(in->a[t].rev == 0) {
						// if(res->a[res->n-1].qs > res->a[res->n-2].qe) fprintf(stderr, "ERROR-3\n");
						if(!(res->a[res->n-2].qs<=res->a[res->n-1].qs && res->a[res->n-1].qs <= res->a[res->n-2].qe 
							&& res->a[res->n-2].qe <= res->a[res->n-1].qe)) {
								fprintf(stderr, "ERROR-3\n");
						}
						// if(res->a[res->n-1].qs == res->a[res->n-2].qe) {
						// 	if(res->a[res->n-2].qe < in->a[t].qe) res->a[res->n-2].qe++;
						// 	else if(res->a[res->n-1].qs > 0) res->a[res->n-1].qs--;
						// }
					} else {
						if(!(res->a[res->n-1].qs<=res->a[res->n-2].qs && res->a[res->n-2].qs <= res->a[res->n-1].qe 
							&& res->a[res->n-1].qe <= res->a[res->n-2].qe)) {
								fprintf(stderr, "ERROR-4\n");
						}
					}
					
					// debug_infer_read_ovlp(uopt, diff_ec_ul, 
					// 	in->a[t].rev?&(res->a[res->n-2]):&(res->a[res->n-1]), 
					// 			in->a[t].rev?&(res->a[res->n-1]):&(res->a[res->n-2]), u, k, k-1);
				}
			}
		} else {///rid
			kv_push(ul_ov_t, *res, in->a[t]);
			// res->a[res->n-1].tn <<= 1; 
			// res->a[res->n-1].tn >>= 1;
			res->a[res->n-1].qn = o->x_id;
			c_occ++;
		}
	}
	if(res->n != cut) {
		radix_sort_ul_ov_srt_qe(res->a + cut, res->a + res->n);
		if(c_occ) {
			int64_t ci, cn = cut;
			for (k = cut; k < res->n; k++) {
				z = &(res->a[k]);
				if(z->tn&((uint32_t)(0x80000000))) continue;
				if(ct->is_c.a[z->tn] == 0) continue;

				for (ci = k+1; ci < (int64_t)(res->n); ci++) {
					if(res->a[ci].qe > z->qe) break;
					if(res->a[ci].qn == (uint32_t)-1) continue;
					if(!(res->a[ci].tn&((uint32_t)(0x80000000)))) continue;
					if(z->qs <= res->a[ci].qs && z->qe >= res->a[ci].qe) {
						if(check_contain_pair(uopt, (res->a[ci].tn<<1)>>1, z->tn,  1)) {
							res->a[ci].qn = (uint32_t)-1;
						}
					}
				}

				for (ci = k-1; ci >= cn; ci--) {
					if(res->a[ci].qe <= z->qs) break;
					if(res->a[ci].qn == (uint32_t)-1) continue;
					if(!(res->a[ci].tn&((uint32_t)(0x80000000)))) continue;
					if(z->qs <= res->a[ci].qs && z->qe >= res->a[ci].qe) {
						if(check_contain_pair(uopt, (res->a[ci].tn<<1)>>1, z->tn,  1)) {
							res->a[ci].qn = (uint32_t)-1;
						}
					}
				}
			}

			// debug_contain_ovlps(res->a+cut, res->n-cut, uopt);

			for (k = l = cut; k < res->n; k++) {
				if(res->a[k].qn == (uint32_t)-1) continue;
				if(k != l) {
					res->a[l] = res->a[k];
				}
				res->a[l].tn <<= 1; res->a[l].tn >>= 1;
				++l;
			}
			res->n = l;
		}
	}

	return res->n - cut;
}


int64_t gl_chain_refine(overlap_region_alloc* olist, Correct_dumy* dumy, haplotype_evdience_alloc *hap, glchain_t *ll, const ul_idx_t *uref, double diff_ec_ul, int64_t winLen, int64_t qlen, const ug_opt_t *uopt, void *km)
{
	ll->tk.n = ll->lo.n = 0;
	kv_ul_ov_t *idx = &(ll->lo);
	ul_contain *ct = uref->ct;
	gl_chain_gen(olist, uref, idx, 0, hap, km);
	if(idx->n == 0) return 0;
	kv_resize_km(km, ul_ov_t, ll->tk, idx->n);
	kv_resize_km(km, uint64_t, ll->srt.a, idx->n);
	kv_resize_km(km, uint64_t, hap->snp_srt, idx->n);
 	if(gl_exact_chain(idx, &(ll->tk), uref, G_CHAIN_BW, diff_ec_ul, qlen, dumy->overlapID, ll->srt.a.a, hap->snp_srt.a, km)) {
		kv_ul_ov_t *chains = &(ll->tk); ul_ov_t *p = NULL; 
		uint64_t k, z, ff, s, e, an, cn, sft = 50, si = 0, ei = 0, resc = 0, chains_pl = chains->n;
		radix_sort_ul_ov_srt_qs(idx->a, idx->a + idx->n);
		for (k = 0; k < olist->length; k++) {
			if(olist->list[k].is_match!=2) continue;

			s = olist->list[k].x_pos_s; e = olist->list[k].x_pos_e+1;
			for (z = ff = 0; z < idx->n; z++) {
				if((s+sft) >= idx->a[z].qs && e <= (idx->a[z].qe+sft)) {
					ff = 1;
					break;
				}
				if(idx->a[z].qs >= (e+sft)) break;
			}
			cn = ((uint32_t)(ct->idx.a[olist->list[k].y_id]));
			if(ff && cn==0) continue;
			an = update_ava_het_site(hap, k, &si, &ei, cn);
			// if(an != get_het_site(hap, k)) fprintf(stderr, "an->%lu, get_het_site->%lu\n", an, get_het_site(hap, k));
			
			if(cn > 0 && an > 0) {
				resc += rescue_contain_ul_chains(uref, &(olist->list[k]), hap->list+si, an, 
				ct->rids.a + ((ct->idx.a[olist->list[k].y_id])>>32), cn, chains, diff_ec_ul, winLen, 0, km);
			}

			if(ff == 0) {
				kv_pushp_km(km, ul_ov_t, *chains, &p);
				p->qn = k/**olist->list[k].x_id**/; p->qs = olist->list[k].x_pos_s; p->qe = olist->list[k].x_pos_e+1;
				p->tn = olist->list[k].y_id; p->el = 0; p->sec = (an&((uint64_t)0x3FFFFFFF))/**get_het_site(hap, k)**/; 
				p->rev = olist->list[k].y_pos_strand;
				if(p->rev) {
					p->ts = uref->ug->u.a[p->tn].len - (olist->list[k].y_pos_e+1); 
					p->te = uref->ug->u.a[p->tn].len - olist->list[k].y_pos_s;
				} else {
					p->ts = olist->list[k].y_pos_s; 
					p->te = olist->list[k].y_pos_e+1;
				}
				
			}
			si = ei;
		}

		if(resc > 0) {
			radix_sort_ul_ov_srt_tn(chains->a + chains_pl, chains->a + chains->n);
			ff = dedup_sort_ul_ov_t(chains->a + chains->n - resc, resc);
			chains->n = chains->n - resc + ff;
		}
	}

	if(ll->tk.n > 0) infer_read_ovlp(uref, olist, &(ll->tk), &(ll->lo), diff_ec_ul, winLen, uopt, ct, km);
	else ll->lo.n = 0;
	return 1;
}


/**
void fill_edge_weight(ul_ov_t *a, int64_t a_n, const ug_opt_t *uopt, int64_t bw, double diff_ec_ul, int64_t qlen)
{
	uint32_t li_v, lj_v;
	ul_ov_t *li = NULL, *lj = NULL;
	int64_t mm_ovlp, x, i, j, o;
	ma_hit_t *t = NULL;
	for (i = 0; i < a_n; i++) {
		li = &(a[i]); li_v = (li->tn<<1)|li->rev;
		mm_ovlp = max_ovlp_src(uopt, li_v^1);
		x = (li->qs + mm_ovlp)*diff_ec_ul; 
        if(x < bw) x = bw; 
        x += li->qs + mm_ovlp;
        if (x > qlen+1) x = qlen+1;
		x = find_ul_ov_max(i, a, x); 
		for (j = x; j >= 0; --j) { // collect potential destination vertices
            lj = &(a[j]); lj_v = (lj->tn<<1)|lj->rev;
			if(lj->qe <= li->qs) break;
			if(li_v == lj_v) continue;
			t = query_ovlp_src(uopt, li_v, lj_v, ((int64_t)(lj->qe))-((int64_t)(li->qs)), diff_ec_ul);
			if(t) {
				t->bl;
			}
        }
	}
	

}
**/

int64_t get_ecov_adv_back(const ul_idx_t *uref, const ug_opt_t *uopt, uint32_t v, uint32_t w, int64_t bw, double diff_ec_ul, int64_t dq, uint32_t *is_contain)
{
	int64_t dt = -1, dif, mm; if(is_contain) (*is_contain) = 0;
	const asg_t *g = uref?uref->ug->g:NULL;
	uint32_t nv, i; asg_arc_t *av = NULL;
	if(g) {
		nv = asg_arc_n(g, v); av = asg_arc_a(g, v);
		for (i = 0; i < nv; i++) {
			if(av[i].del || av[i].v != w) continue;
			dt = av[i].ol;
			break;
		}
	}

	if(dt < 0 && uopt) {
		ma_hit_t_alloc* src = uopt->sources;
		int64_t min_ovlp = uopt->min_ovlp;
		int64_t max_hang = uopt->max_hang;
		uint64_t z, qn, tn, x = v>>1; int32_t r = 1; asg_arc_t e;
		for (z = 0; z < src[x].length; z++) {
			qn = Get_qn(src[x].buffer[z]); tn = Get_tn(src[x].buffer[z]);
			if(tn != (w>>1)) continue;
			r = ma_hit2arc(&(src[x].buffer[z]), Get_READ_LENGTH(R_INF, qn), Get_READ_LENGTH(R_INF, tn), max_hang, asm_opt.max_hang_rate, min_ovlp, &e);
			if(r < 0) {
				if(r == MA_HT_QCONT || r == MA_HT_TCONT) {
					if(src[x].buffer[z].rev == ((uint32_t)(v^w))) {
						dt = Get_qe(src[x].buffer[z]) - Get_qs(src[x].buffer[z]);
						if(dt < Get_te(src[x].buffer[z]) - Get_ts(src[x].buffer[z])) {
							dt = Get_te(src[x].buffer[z]) - Get_ts(src[x].buffer[z]);
						}
						if(is_contain) (*is_contain) = 1;
						break;
					}
				}
				continue;
			}
			if((e.ul>>32) != v || e.v != w) continue;
			dt = e.ol;
			break;
		}
	}
	if(dt < 0) return 0;
	dif = (dq>dt? dq-dt:dt-dq);
	mm = MAX(dq, dt); mm *= diff_ec_ul; if(mm < bw) mm = bw;
	// if((v>>1) == 1163 && (w>>1) == 1168) fprintf(stderr, ">>>>>>dis_q:%ld, dis_t:%ld, dif:%ld, mm:%ld\n", dis_q, dis_t, dif, mm);
	if(dif <= mm) return 1;
	return 0;
}

ma_hit_t *get_ug_edge_src(ma_ug_t *ug, ma_hit_t_alloc *src, int64_t max_hang, int64_t min_ovlp, uint32_t uv, uint32_t uw)
{
	if(ug->u.a[uv>>1].circ || ug->u.a[uw>>1].circ) return NULL;
	uint32_t v, w, k, qn, tn; int32_t r; asg_arc_t t;
    v = ((uv&1)?(ug->u.a[uv>>1].start^1):(ug->u.a[uv>>1].end^1));
    w = ((uw&1)?(ug->u.a[uw>>1].end):(ug->u.a[uw>>1].start));
    ma_hit_t_alloc *x = &(src[v>>1]);

	for (k = 0; k < x->length; k++) {
        qn = Get_qn(x->buffer[k]);
        tn = Get_tn(x->buffer[k]);
        if(qn == (v>>1) && tn == (w>>1)) {
            r = ma_hit2arc(&(x->buffer[k]), Get_READ_LENGTH(R_INF, v>>1), Get_READ_LENGTH(R_INF, w>>1),
			max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
            if(r < 0) continue;
            if((t.ul>>32)!=v || t.v!=w) continue;
			return &(x->buffer[k]);
        }
    }

	return NULL;
}

///mode: 0->ug; 1->read
int64_t get_ecov_adv(const ul_idx_t *uref, const ug_opt_t *uopt, uint32_t v, uint32_t w, int64_t bw, double diff_ec_ul, int64_t dq, uint64_t mode, int64_t *contain_off)
{
	int64_t dt = -1, dif, mm; (*contain_off) = 0;
	uint32_t nv, i; asg_arc_t *av = NULL; ma_hit_t *x = NULL;
	if(!mode) {
		const asg_t *g = uref?uref->ug->g:NULL;
		nv = asg_arc_n(g, v); av = asg_arc_a(g, v);
		for (i = 0; i < nv; i++) {
			if(av[i].del || av[i].v != w) continue;
			dt = av[i].ol; (*contain_off) = av[i].ou;
			// if(v==1772 && w==1769) fprintf(stderr, "+++v:%u, w:%u, ou:%u\n", v, w, av[i].ou);
			// if((v>>1) == 3012 && (w>>1) == 3011) fprintf(stderr, "******************\n");
			if(av[i].ou >= OU_MASK) {
				x = get_ug_edge_src(uref->ug, uopt->sources, uopt->max_hang, uopt->min_ovlp, 
                                                    						av[i].ul>>32, av[i].v);
				(*contain_off) = x->cc;
				// if(v==1772 && w==1769) fprintf(stderr, "---v:%u, w:%u, cc:%u\n", v, w, x->cc);
			}
			break;
		}
	}else {
		ma_hit_t_alloc* src = uopt->sources;
		int64_t min_ovlp = uopt->min_ovlp;
		int64_t max_hang = uopt->max_hang;
		uint64_t z, qn, tn, x = v>>1; int32_t r = 1; asg_arc_t e;
		for (z = 0; z < src[x].length; z++) {
			qn = Get_qn(src[x].buffer[z]); tn = Get_tn(src[x].buffer[z]);
			if(tn != (w>>1)) continue;
			r = ma_hit2arc(&(src[x].buffer[z]), Get_READ_LENGTH(R_INF, qn), Get_READ_LENGTH(R_INF, tn), max_hang, asm_opt.max_hang_rate, min_ovlp, &e);
			if(r < 0) continue;
			if((e.ul>>32) != v || e.v != w) continue;
			dt = e.ol; (*contain_off) = src[x].buffer[z].cc;
			break;
		}
	}
	if(dt < 0) return 0;
	dif = (dq>dt? dq-dt:dt-dq);
	mm = MAX(dq, dt); mm *= diff_ec_ul; if(mm < bw) mm = bw;
	// if((v>>1) == 1163 && (w>>1) == 1168) fprintf(stderr, ">>>>>>dis_q:%ld, dis_t:%ld, dif:%ld, mm:%ld\n", dis_q, dis_t, dif, mm);
	if(dif <= mm) return 1;
	return 0;
}

void get_rr_tse(const ul_idx_t *uref, ul_ov_t *li, uint32_t *ts, uint32_t *te, uint32_t *tl)
{
	(*tl) = uref?uref->ug->g->seq[li->tn].len:Get_READ_LENGTH(R_INF, li->tn); 
	if(!(li->rev)) {
		(*ts) = li->ts; (*te) = li->te;
	} else {
		(*ts) = (*tl) - li->te; (*te) = (*tl) - li->ts;
	}
}
/**
uint32_t checkM(uint32_t v, uint32_t l, const ul_idx_t *uref, const asg_t *g, uint32_t in, uint32_t its, 
uint32_t iqs, uint32_t iqe, int64_t bw, double diff_ec_ul, int64_t qlen, ul_ov_t *a, uint32_t a_n) 
{
	int64_t vl = uref?uref->ug->g->seq[v>>1].len:Get_READ_LENGTH(R_INF, (v>>1)), t_dis, q_dis, kcs, mm_ovlp, x, k; 
	t_dis = ((int64_t)(l + vl)) - ((int64_t)(in - its));
	kcs = iqs; kcs -= t_dis; if(kcs < 0) kcs = 0;
	uint32_t nv = asg_arc_n(g, v), i, lk_v, kts, kte, kn; 
	asg_arc_t *av = asg_arc_a(g, v), *p = NULL; mm_ovlp = -1; ul_ov_t *lk;
	for (i = 0, p = NULL; i < nv; i++) {
		if(av[i].del) continue;
		if((int32_t)(av[i].ol) > mm_ovlp) {
			p = &(av[i]); mm_ovlp = av[i].ol;
		}
	}
	if(!p) return 0;

	x = (kcs + mm_ovlp)*diff_ec_ul; 
	if(x < bw) x = bw; 
	x += kcs + mm_ovlp;
	if (x > qlen+1) x = qlen+1;
	x = find_ul_ov_max(a_n, a, x); 
	for (k = x; k >= 0; --k) {
		lk = &(a[k]); lk_v = ((lk->tn<<1)|lk->rev)^1;
		if(lk->qe <= kcs) break;//evan this pair has a overlap, its length will be very small; just ignore
    	if(lk->qs >= kcs) continue; // lk is contained in li on the query coordinate

		get_rr_tse(uref, lk, &kts, &kte, &kn);
		///t_dis and q_dis might be < 0
		t_dis = ((int64_t)(l + kn - kte)) - ((int64_t)(in - its));
		q_dis = ((int64_t)(iqs)) - ((int64_t)(lk->qe));
	}

}

void best_path_ext(const ul_idx_t *uref, const ug_opt_t *uopt, int64_t g_gap, ul_ov_t *a, uint32_t a_n,
int64_t bw, double diff_ec_ul, uint64_t *track, ul_ov_t *li)
{
	if(a_n <= 0) return;
	const asg_t *g = uref?uref->ug->g:NULL; asg_arc_t *av = NULL, *p = NULL; ul_ov_t *lk;
	uint32_t nv, i, v, io, in, its, ite, kn, kts, kte; int64_t mm, l, max_dist, k, t_dis, q_dis; 

	get_rr_tse(uref, li, &its, &ite, &in);
	io = in - ite; max_dist = g_gap + in - its;

	if(g) {
		v = (((li->tn<<1)|li->rev)^1); mm = 1; l = 0;
		while (mm >= 0) {

			



			nv = asg_arc_n(g, v); av = asg_arc_a(g, v); mm = -1;
			for (i = 0, p = NULL; i < nv; i++) {
				if(av[i].del) continue;
				if((int32_t)(av[i].ol) > mm) {
					p = &(av[i]); mm = av[i].ol;
				}
			}
			if(p) {
				l += (uint32_t)(p->ul); 
				if(l > max_dist) break;
				for (k = a_n-1; k >= 0; k--) {
					lk = &(a[k]);
					if((lk->qe+g_gap) <= li->qs) break; 
					if(p->v == (((lk->tn<<1)|lk->rev)^1)) {
						///check if lk can be directly reachedc from li
						if(lk->qe > li->qs && (track[k]&((uint64_t)0x80000000))) {
							;
						}
						get_rr_tse(uref, lk, &kts, &kte, &kn);
						///t_dis and q_dis might be < 0
						t_dis = ((int64_t)(l + kn - kte)) - ((int64_t)(in - its));
						q_dis = ((int64_t)(li->qs)) - ((int64_t)(lk->qe));
					}
				}
				v = p->v;
			}
		}
	}
}


int64_t gl_chain_advance(kv_ul_ov_t *res, ul_ov_t *ex, const ul_idx_t *uref, const ug_opt_t *uopt, int64_t bw, 
double diff_ec_ul, int64_t qlen, int64_t max_skip, uint64_t *srt, uint64_t *idx, uint64_t *track, float trav_rate, void *km)
{
	uint32_t li_v, lj_v, rev_n, gapLen = (trav_rate>0?trav_rate*qlen:0);
	int64_t mm_ovlp, x, i, j, k, sc, csc, mm_sc, mm_idx, qo, n_skip, n_all;
	ul_ov_t *li = NULL, *lj = NULL, rev_t;
	radix_sort_ul_ov_srt_qe(res->a, res->a + res->n);
	for (i = 0; i < (int64_t)res->n; ++i) {
		li = &(res->a[i]); li_v = (li->tn<<1)|li->rev;
		// if(!(li->el)) continue;
		mm_ovlp = uref?max_ovlp(uref->ug->g, li_v^1):max_ovlp_src(uopt, li_v^1);
		x = (li->qs + mm_ovlp)*diff_ec_ul; 
		if(x < bw) x = bw; 
		x += li->qs + mm_ovlp;
		if (x > qlen+1) x = qlen+1;
		x = find_ul_ov_max(i, res->a, x); 
		if(li->el) {
			csc = uref?retrieve_u_cov_region(uref, li->tn, 0, li->ts, li->te, NULL):li->te-li->ts;
		} else {
			csc = -1;///for cis overlap, the csc should be >1000; so -1 for trans overlaps should be fine
		}
		mm_sc = csc; mm_idx = -1; n_skip = n_all = 0;
		for (j = x; j >= 0; --j) { // collect potential destination vertices
			lj = &(res->a[j]); lj_v = (lj->tn<<1)|lj->rev;
			// if((lj->qe+gapLen) <= li->qs) break; 
			if(lj->qe <= li->qs) break;//evan this pair has a overlap, its length will be very small; just ignore
			if(lj->qs >= li->qs) continue; // lj is contained in li on the query coordinate
			qo = infer_rovlp(li, lj, NULL, NULL); ///overlap length in query (UL read)
			if(li_v != lj_v && get_ecov_adv(uref, uopt, li_v^1, lj_v^1, bw, diff_ec_ul, qo)) {
				sc = csc + (track[j]>>32);
				if(sc > mm_sc) mm_sc = sc, mm_idx = j;
				if(res->a[j].sec == i && res->a[j].el) n_skip++;
				if((track[j]&((uint64_t)0x7FFFFFFF)) != ((uint64_t)0x7FFFFFFF)) {
					res->a[(track[j]&((uint64_t)0x7FFFFFFF))].sec = i;
				}
				track[j] |= ((uint64_t)0x80000000); 
			} else {
				if(track[j]&((uint64_t)0x80000000)) track[j] -= ((uint64_t)0x80000000);
			}
			n_all++;
		}

		if(n_all > max_skip) n_all = max_skip;
		else n_all -= 2; //allow one mismatch; note here must be -2

		if(li->el && (mm_idx<0 || n_skip<n_all)) {///do not extend for trans overlaps
			gapLen = (trav_rate>0?trav_rate*qlen*li->el:0);
			if(gapLen > 0) {
				// if((lj->qe+gapLen) <= li->qs) break; 
				///since graph traversal just has one path, so this step might be quite easy
				
				best_path_ext(uref, uopt, gapLen, res->a, x+1, bw, diff_ec_ul, li);
			}
		}
		// 4294967295L
		track[i] = mm_sc; track[i] <<= 32;
		track[i] |= (mm_idx>=0?mm_idx:((uint64_t)0x7FFFFFFF));
		srt[i] = mm_sc; srt[i] <<= 32; srt[i] |= i;
		// fprintf(stderr, "+++i:%ld, mm_idx:%ld, mm_sc:%ld\n", i, mm_idx, mm_sc);
		// fprintf(stderr, "[M::utg%.6d%c] qs->%u; qe->%u\n\n", li->tn+1, "lc"[uref->ug->u.a[li->tn].circ], li->qs, li->qe);
	}

	for (i = 0; i < (int64_t)res->n; ++i) {
		if(track[i]&((uint64_t)0x80000000)) track[i] -= ((uint64_t)0x80000000);
	}
	
	int64_t n_v, n_u, n_v0; 
	radix_sort_gfa64(srt, srt+res->n); //ex->n = res->n;
	for (k = (int64_t)res->n-1, n_v = n_u = 0; k >= 0; --k) {
		// fprintf(stderr, "\nk:%ld\n", k);
		n_v0 = n_v;
		for (i = (uint32_t)srt[k]; i >= 0 && (track[i]&((uint64_t)0x80000000)) == 0;) {
				ex[n_v++] = res->a[i]; track[i] |= ((uint64_t)0x80000000);
				// fprintf(stderr, "+i:%ld, ", i);
				// fprintf(stderr, "[M::utg%.6d%c] qs->%u; qe->%u\n", res->a[i].tn+1, "lc"[uref->ug->u.a[res->a[i].tn].circ], res->a[i].qs, res->a[i].qe);
				if((track[i]&((uint64_t)0x7FFFFFFF)) == ((uint64_t)0x7FFFFFFF)) i = -1;
				else i = track[i]&((uint64_t)0x7FFFFFFF);
		}
		if(n_v0 == n_v) continue;
		///keep the whole score; do not cut score like minigraph
		sc = (i<0?(srt[k]>>32):((srt[k]>>32)-(track[i]>>32)));
		// sc = srt[k]>>32;
		idx[n_u++] = ((uint64_t)sc<<32)|(n_v-n_v0);
	}

	for (k = 0, n_v = n_v0 = 0; k < n_u; k++) {
		n_v0 = n_v; n_v += (uint32_t)idx[k];
		res->a[k].qn = idx[k]>>32;
		res->a[k].ts = n_v0; res->a[k].te = n_v;
		
		rev_n = ((uint32_t)idx[k])>>1; 
		///we need to consider contained reads; so determining qs is not such easy
		res->a[k].qs = (uint32_t)-1; res->a[k].qe = ex[n_v0].qe;
		for (i = 0; i < rev_n; i++) {
			rev_t = ex[n_v0+i]; 
			ex[n_v0+i] = ex[n_v0+rev_n-i-1];
			ex[n_v0+rev_n-i-1] = rev_t;
			if(res->a[k].qs > ex[n_v0+i].qs) res->a[k].qs = ex[n_v0+i].qs;
			if(res->a[k].qs > ex[n_v0+rev_n-i-1].qs) res->a[k].qs = ex[n_v0+rev_n-i-1].qs;
		}
		if(i < ((uint32_t)idx[k]) && res->a[k].qs < ex[n_v0+i].qs) {
			res->a[k].qs = ex[n_v0+i].qs;
		}
	}
	res->n = n_u;
	return res->n;
}
**/

int64_t determine_containment_chain(const ug_opt_t *uopt, uint64_t *track, uint64_t *flag, kv_ul_ov_t *res, int32_t nc, int64_t *nsc, int64_t mm_idx, int64_t bw, double diff_ec_ul, uint32_t el, All_reads *ridx, ma_ug_t *ug)
{
	int64_t i, k, pk, ak, e, off = 128, qo, tt = 0, ii; ul_ov_t *li = NULL, *lk = NULL;
	uint32_t li_v, lk_v, is_c;
	if(nc<=0) return 0;
	for (k = tt = ak = 0; k < nc; k++) {
		if(!(flag[res->a[k].sec]&((uint64_t)0x80000000))) {
			flag[res->a[k].sec] |= ((uint64_t)0x80000000); tt++;
		} else {
			res->a[ak++].sec = res->a[k].sec;
		}
	}
	if(tt==nc) return nsc[0] - nsc[1];
	assert(ak>0);
	e = res->a[res->a[ak-1].sec].qe;//e is the smallest qe
	for (i = mm_idx, pk = 0; i >= 0;) {///i++, li->qe--
		li = &(res->a[i]); li_v = (li->tn<<1)|li->rev;
		if((track[i]&((uint64_t)0x7FFFFFFF)) == ((uint64_t)0x7FFFFFFF)) i = -1;
		else i = track[i]&((uint64_t)0x7FFFFFFF);
		if(li->qe + off < e || tt == nc) break;//128 is used to tolerate indels; 
		for (k = pk, ii = 0; k < ak; k++) {///k++, lk->qe--
			if(res->a[k].sec == ((uint32_t)0x3FFFFFFF)) continue;
			lk = &(res->a[res->a[k].sec]); lk_v = (lk->tn<<1)|lk->rev;
			if(li->qe + off >= lk->qe) {
				if(ii == 0) pk = k;
				if(li->qs <= lk->qs + off) {
					qo = infer_rovlp(li, lk, NULL, NULL, ridx, ug); ///overlap length in query (UL read)
					if(li_v != lk_v && get_ecov_adv_back(NULL, uopt, li_v^1, lk_v^1, bw, diff_ec_ul, qo, &is_c)) {
						if(is_c) {
							tt++; res->a[k].sec = ((uint32_t)0x3FFFFFFF);
							if(el) nsc[0] -= ((int64_t)(lk->te-lk->ts));
							else nsc[!(lk->el)] -= ((int64_t)(lk->te-lk->ts));
						}
					}
				}
				ii = 1;
			}
		}
	}

	assert(nsc[0]>=0 && nsc[1]>=0);
	return nsc[0] - nsc[1];
}

uint64_t push_sc_pre(int64_t mm_sc, int64_t mm_idx)
{
	uint32_t sc = (mm_sc>=0?(((uint32_t)(mm_sc))|((uint32_t)(0x80000000))):((uint32_t)(-mm_sc)));
	uint64_t x = sc; x <<= 32; x |= (mm_idx>=0?mm_idx:((uint64_t)0x7FFFFFFF));
	return x;
}

int64_t pop_sc(uint64_t x)
{
	int64_t sc;
	x >>= 32;
	if(x&((uint64_t)(0x80000000))) {
		sc = x - ((uint64_t)(0x80000000));
	} else {
		sc = x; sc *= -1;
	}
	return sc;
}

int64_t pop_pre(uint64_t x)
{
	if((x&((uint64_t)0x7FFFFFFF)) == ((uint64_t)0x7FFFFFFF)) return -1;
	else return (x&((uint64_t)0x7FFFFFFF));
}

///mode: 0->ug; 1->read
int64_t gl_chain_advance(kv_ul_ov_t *res, ul_ov_t *ex, const ul_idx_t *uref, const ug_opt_t *uopt, int64_t bw, 
double diff_ec_ul, int64_t qlen, int64_t max_skip, uint64_t *srt, uint64_t *idx, uint64_t *track, int64_t trans_sc, 
uint64_t mode, All_reads *ridx, ma_ug_t *ug, int64_t debug_i, void *km)
{
	// fprintf(stderr, "\n+++[M::%s] res->n:%u\n", __func__, (uint32_t)res->n);
	if(res->n == 0) return 0;
	uint32_t li_v, lj_v, rev_n;
	int64_t mm_ovlp, x, i, j, k, sc, csc, mm_sc, mm_idx, qo, share, n_el = 0;
	ul_ov_t *li = NULL, *lj = NULL, rev_t;
	radix_sort_ul_ov_srt_qe(res->a, res->a + res->n);
	for (i = 1, j = 0; i <= (int64_t)res->n; i++) {
        if (i == (int64_t)res->n || res->a[i].qe != res->a[j].qe) {
            if(i - j > 1) {
                radix_sort_ul_ov_srt_qs(res->a+j, res->a+i);
            }
            j = i;
        }
    }
	for (i = 0; i < (int64_t)res->n; ++i) {
		li = &(res->a[i]); li_v = (li->tn<<1)|li->rev;
		mm_ovlp = mode?max_ovlp_src(uopt, li_v^1):max_ovlp(uref->ug->g, li_v^1);		
		x = (li->qs + mm_ovlp)*diff_ec_ul; 
		if(x < bw) x = bw; 
		x += li->qs + mm_ovlp;
		if (x > qlen+1) x = qlen+1;
		x = find_ul_ov_max(i, res->a, x+G_CHAIN_INDEL); 
		if(li->el) csc = mode?retrieve_r_cov_region(uref, li->tn, 0, li->ts, li->te, NULL):retrieve_u_cov_region(uref, li->tn, 0, li->ts, li->te, NULL);
		else csc = (trans_sc*li->sec); //trans overlaps
		mm_sc = csc; mm_idx = -1; 
		// if(i == 37 || i == 36 || i == 35 || i == 32) fprintf(stderr, "*i:%ld, x:%ld, mm_sc:%ld\n", i, x, mm_sc);
		for (j = x; j >= 0; --j) { // collect potential destination vertices
			lj = &(res->a[j]); lj_v = (lj->tn<<1)|lj->rev;
			// if((lj->qe+gapLen) <= li->qs) break; 
			if(lj->qe+G_CHAIN_INDEL <= li->qs) break;//even this pair has a overlap, its length will be very small; just ignore
			// if(lj->qs >= li->qs+G_CHAIN_INDEL) continue; // lj is contained in li on the query coordinate; 128 for indel offset
			if(lj->qs >= li->qs) continue;
			qo = infer_rovlp(li, lj, NULL, NULL, ridx, ug); ///overlap length in query (UL read)
			// if(i == 37 || i == 36 || i == 35 || i == 32) fprintf(stderr, ">i:%ld, j:%ld, qo:%ld\n", i, j, qo);
			if(li_v != lj_v && get_ecov_adv(uref, uopt, li_v^1, lj_v^1, bw, diff_ec_ul, qo, mode, &share)) {
				// if(i == 37 || i == 36 || i == 35 || i == 32) fprintf(stderr, "#i:%ld, j:%ld, share:%ld\n", i, j, share);
				sc = csc + pop_sc(track[j]);
				// if((!mode)&&i==11&&j==10) {
				// 	fprintf(stderr,"+share:%ld, i:%ld, j:%ld, li_v^1:%u, lj_v^1:%u\n", 
				// 															share, i, j, li_v^1, lj_v^1);
				// }
				// if((mode&&i==21&&j==20) || (mode&&i==22&&j==21) || (mode&&i==23&&j==22)) {
				// 	fprintf(stderr,"-share:%ld, i:%ld, j:%ld, li_v^1:%u, lj_v^1:%u\n", 
				// 															share, i, j, li_v^1, lj_v^1);
				// }
				if(li->el && lj->el) sc -= (share>=csc?csc:share);///csc must be larger than 0
				// if((!li->el) && (!lj->el)) sc -= ((share>=o_csc?o_csc:share)*(-trans_scl));
				if(sc > mm_sc) mm_sc = sc, mm_idx = j;
			}
		}

		track[i] = push_sc_pre(mm_sc, mm_idx);
		srt[i] = track[i]>>32; srt[i] <<= 32; srt[i] |= i;
		n_el += li->el;

		// if(mode) {
		// 	fprintf(stderr, "[M::%.*s] i:%ld, li->el:%u, li->score:%ld (raw_sc:%u), mm_idx:%ld, mm_sc:%ld, q[%u, %u), t[%u, %u), rev:%c\n", 
		// 		(int32_t)Get_NAME_LENGTH(R_INF, li->tn), Get_NAME(R_INF, li->tn), i, li->el, csc, li->te - li->ts, mm_idx, mm_sc, li->qs, li->qe, li->ts, li->te, "+-"[li->rev]);
		// } else {
		// 	fprintf(stderr, "[M::utg%.6u%c] i:%ld, li->el:%u, li->score:%ld (raw_sc:%u), mm_idx:%ld, mm_sc:%ld, q[%u, %u), t[%u, %u), rev:%c\n", 
		// 		li->tn+1, "lc"[uref->ug->u.a[li->tn].circ], i, li->el, csc, li->te - li->ts, mm_idx, mm_sc, li->qs, li->qe, li->ts, li->te, "+-"[li->rev]);
		// }
		// if(!mode) {
		// 	fprintf(stderr, "[M::utg%.6d%c] qs->%u; qe->%u\n", li->tn+1, "lc"[uref->ug->u.a[li->tn].circ], li->qs, li->qe);
		// }
	}
	
	int64_t n_v, n_u, n_v0, le, lnv; 
	radix_sort_gfa64(srt, srt+res->n); 
	for (k = (int64_t)res->n-1, n_v = n_u = 0; k >= 0; --k) {
		n_v0 = n_v; i = (uint32_t)srt[k];
		if(res->a[i].el) { ///chain must start from cis alignments
			for (le = -1; i >= 0 && (track[i]&((uint64_t)0x80000000)) == 0;) {
				if(res->a[i].el) {
					le = -1;
				}else if(n_v>n_v0 && ex[n_v-1].el) {
					le = i; lnv = n_v;///cut the cis alignments in the end
				}
				ex[n_v++] = res->a[i]; track[i] |= ((uint64_t)0x80000000);
				// fprintf(stderr, "[M::utg%.6d%c] qs->%u; qe->%u\n", res->a[i].tn+1, "lc"[uref->ug->u.a[res->a[i].tn].circ], res->a[i].qs, res->a[i].qe);
				i = pop_pre(track[i]);
			}
		}
		if(n_v0 == n_v) continue;
		if(le >= 0) {
			i = le; n_v = lnv;
		}
		if(n_v0 == n_v) continue;
		// fprintf(stderr, "[++chain::] beg_idx->%u, end_idx->%ld, le->%ld, chain_n->%ld\n", (uint32_t)srt[k], i, le, n_v - n_v0);
		///keep the whole score; do not cut score like minigraph
		// sc = pop_sc(srt[k]);
		sc = (i<0?(pop_sc(srt[k])):(pop_sc(srt[k])-pop_sc(track[i])));
		// fprintf(stderr, "++[M::%s] k:%ld, n_v0:%ld, n_v:%ld, le:%ld, sc:%ld, beg:%u, end:%ld, p_score:%ld, cut_score:%ld\n", 
		// __func__, k, n_v0, n_v, le, sc, (uint32_t)srt[k], i, pop_sc(srt[k]), i<0?0:pop_sc(track[i]));
		if(sc /**<=**/< 0) {///sc might be 0, if the UL alignment cannot cover the whole overlap between two HiFi reads
			n_v = n_v0;
			continue;
		}
		// idx[n_u++] = push_sc_pre(sc, n_v-n_v0);	
		idx[n_u++] = ((uint64_t)sc<<32)|(n_v-n_v0);	
	}
	// fprintf(stderr, "[M::%s] n_u:%ld, n_v:%ld\n", __func__, n_u, n_v);
	for (k = 0, n_v = n_v0 = 0; k < n_u; k++) {
		n_v0 = n_v; n_v += (uint32_t)idx[k];
		// fprintf(stderr, "[M::%s] k:%ld, n_v0:%ld, n_v:%ld\n", __func__, k, n_v0, n_v);
		res->a[k].qn = idx[k]>>32;//score
		res->a[k].ts = n_v0; res->a[k].te = n_v;///idx
		
		rev_n = ((uint32_t)idx[k])>>1; 
		///we need to consider contained reads; so determining qs is not such easy
		res->a[k].qs = (uint32_t)-1; res->a[k].qe = ex[n_v0].qe;
		for (i = 0; i < rev_n; i++) {
			rev_t = ex[n_v0+i]; ex[n_v0+i] = ex[n_v-i-1]; ex[n_v-i-1] = rev_t;
			
			if(res->a[k].qs > ex[n_v0+i].qs) res->a[k].qs = ex[n_v0+i].qs;
			if(res->a[k].qs > ex[n_v-i-1].qs) res->a[k].qs = ex[n_v-i-1].qs;
			n_el -= ex[n_v0+i].el; n_el -= ex[n_v-i-1].el;
			ex[n_v0+i].sec = ex[n_v-i-1].sec = SEC_MODE; 
		}
		if(((uint32_t)idx[k])&1) {
			if(res->a[k].qs > ex[n_v0+i].qs) res->a[k].qs = ex[n_v0+i].qs;
			n_el -= ex[n_v0+i].el; ex[n_v0+i].sec = SEC_MODE; 
		}
		assert(ex[n_v0].el && ex[n_v-1].el);
		// fprintf(stderr, "[M::%s] k:%ld, qs:%u, qe:%u, chain_occ:%u, chain_score:%u\n", __func__, k, 
		// 					res->a[k].qs, res->a[k].qe, res->a[k].te - res->a[k].ts, res->a[k].qn);
	}
	// if(n_el) {
	// 	fprintf(stderr, "[M::%s] debug_i->%ld, n_el->%ld, n_u->%ld, n_v->%ld\n", __func__, debug_i, n_el, n_u, n_v);
	// }
	assert(n_el == 0);
	res->n = n_u;
	radix_sort_ul_ov_srt_qn(res->a, res->a + res->n);//sort by score
	// fprintf(stderr, "---[M::%s] n_u:%ld, n_v:%ld\n", __func__, n_u, n_v);
	return n_v;
}


int64_t gl_chain_advance_back(kv_ul_ov_t *res, ul_ov_t *ex, const ul_idx_t *uref, const ug_opt_t *uopt, int64_t bw, 
double diff_ec_ul, int64_t qlen, int64_t max_skip, uint64_t *srt, uint64_t *idx, uint64_t *track, float trans_allow, 
All_reads *ridx, ma_ug_t *ug, void *km)
{
	uint32_t li_v, lj_v, rev_n, is_c, nc, s_nc;
	int64_t mm_ovlp, x, i, j, k, sc, csc, mm_sc, mm_idx, qo, trans_scl = (int64_t)(((float)(1))/trans_allow), nsc[2];
	ul_ov_t *li = NULL, *lj = NULL, rev_t;
	radix_sort_ul_ov_srt_qe(res->a, res->a + res->n);
	for (i = s_nc = 0; i < (int64_t)res->n; ++i) {
		li = &(res->a[i]); li_v = (li->tn<<1)|li->rev;
		mm_ovlp = uref?max_ovlp(uref->ug->g, li_v^1):max_ovlp_src(uopt, li_v^1);
		x = (li->qs + mm_ovlp)*diff_ec_ul; 
		if(x < bw) x = bw; 
		x += li->qs + mm_ovlp;
		if (x > qlen+1) x = qlen+1;
		x = find_ul_ov_max(i, res->a, x); 
		csc = uref?retrieve_u_cov_region(uref, li->tn, 0, li->ts, li->te, NULL):li->te-li->ts;
		if(!(li->el)) {
			csc *= -trans_scl; //trans overlaps
			if(csc >= 0) csc = -1;
		}
		mm_sc = csc; mm_idx = -1; nc = nsc[0] = nsc[1] = 0;
		for (j = x; j >= 0; --j) { // collect potential destination vertices
			lj = &(res->a[j]); lj_v = (lj->tn<<1)|lj->rev;
			// if((lj->qe+gapLen) <= li->qs) break; 
			if(lj->qe <= li->qs) break;//even this pair has a overlap, its length will be very small; just ignore
			// if(lj->qs >= li->qs) continue; // lj is contained in li on the query coordinate
			qo = infer_rovlp(li, lj, NULL, NULL, ridx, ug); ///overlap length in query (UL read)
			if(li_v != lj_v && get_ecov_adv_back(uref, uopt, li_v^1, lj_v^1, bw, diff_ec_ul, qo, &is_c)) {
				if(!is_c) {
					sc = csc + pop_sc(track[j]);
					if(sc > mm_sc) mm_sc = sc, mm_idx = j;
				} else if(!uref) {///with uref, retrieve_u_cov_region has already consider contained reads
					res->a[nc++].sec = j; 
					if(li->el) nsc[0] += lj->te-lj->ts;
					else nsc[!(lj->el)] += lj->te-lj->ts;
				}
			}
		}
		
		if(nc && (!uref) && mm_idx>=0) {///deal with containments
			mm_sc += determine_containment_chain(uopt, track, srt, res, nc, nsc, mm_idx, bw, diff_ec_ul, li->el, ridx, ug);
			s_nc++;
		}

		track[i] = push_sc_pre(mm_sc, mm_idx);
		srt[i] = track[i]>>32; srt[i] <<= 32; srt[i] |= i;
		// fprintf(stderr, "[M::utg%.6d%c] qs->%u; qe->%u\n\n", li->tn+1, "lc"[uref->ug->u.a[li->tn].circ], li->qs, li->qe);
	}

	if(s_nc) {
		for (i = 0; i < (int64_t)res->n; ++i) {
			if(srt[i]&((uint64_t)0x80000000)) srt[i]-=((uint64_t)0x80000000);
		}
	}
	
	int64_t n_v, n_u, n_v0, le, lnv; 
	radix_sort_gfa64(srt, srt+res->n); //ex->n = res->n;
	for (k = (int64_t)res->n-1, n_v = n_u = 0; k >= 0; --k) {
		n_v0 = n_v; i = (uint32_t)srt[k];
		if(i>=0 && (res->a[i].el)) { ///chain must start from cis alignments
			for (le = -1; i >= 0 && (track[i]&((uint64_t)0x80000000)) == 0;) {
				if(res->a[i].el) {
					le = -1;
				}else if(n_v>n_v0 && ex[n_v-1].el) {
					le = i; lnv = n_v;///cut the cis alignments in the end
				}
				ex[n_v++] = res->a[i]; track[i] |= ((uint64_t)0x80000000);
				// fprintf(stderr, "[M::utg%.6d%c] qs->%u; qe->%u\n", res->a[i].tn+1, "lc"[uref->ug->u.a[res->a[i].tn].circ], res->a[i].qs, res->a[i].qe);
				i = pop_pre(track[i]);
			}
		}
		if(n_v0 == n_v) continue;
		if(le >= 0) {
			i = le; n_v = lnv;
		}
		if(n_v0 == n_v) continue;
		///keep the whole score; do not cut score like minigraph
		// sc = pop_sc(srt[k]);
		sc = (i<0?(pop_sc(srt[k])):(pop_sc(srt[k])-pop_sc(track[i])));
		if(sc < 0) {
			n_v = n_v0;
			continue;
		}
		// idx[n_u++] = push_sc_pre(sc, n_v-n_v0);	
		idx[n_u++] = ((uint64_t)sc<<32)|(n_v-n_v0);	
	}

	for (k = 0, n_v = n_v0 = 0; k < n_u; k++) {
		n_v0 = n_v; n_v += (uint32_t)idx[k];
		res->a[k].qn = idx[k]>>32;
		res->a[k].ts = n_v0; res->a[k].te = n_v;
		
		rev_n = ((uint32_t)idx[k])>>1; 
		///we need to consider contained reads; so determining qs is not such easy
		res->a[k].qs = (uint32_t)-1; res->a[k].qe = ex[n_v0].qe;
		for (i = 0; i < rev_n; i++) {
			rev_t = ex[n_v0+i]; 
			ex[n_v0+i] = ex[n_v0+rev_n-i-1];
			ex[n_v0+rev_n-i-1] = rev_t;
			if(res->a[k].qs > ex[n_v0+i].qs) res->a[k].qs = ex[n_v0+i].qs;
			if(res->a[k].qs > ex[n_v0+rev_n-i-1].qs) res->a[k].qs = ex[n_v0+rev_n-i-1].qs;
		}
		if(i < ((uint32_t)idx[k]) && res->a[k].qs < ex[n_v0+i].qs) {
			res->a[k].qs = ex[n_v0+i].qs;
		}
	}
	res->n = n_u;
	radix_sort_ul_ov_srt_qn(res->a, res->a + res->n);
	return n_v;
}

uint32_t check_trans_rate(ul_ov_t *a, int64_t a_n, float trans_thres)
{
	uint32_t sp = (uint32_t)-1, ep = (uint32_t)-1, tts = (uint32_t)-1, tte = 0, el = 0, iel = 0; 
	int64_t k;
	for (k = a_n-1; k >= 0; k--) {
		if(a[k].qs < tts) tts = a[k].qs;
		if(a[k].qe > tte) tte = a[k].qe;
		if(!(a[k].el)) continue;
		if(sp == (uint32_t)-1 || a[k].qe <= sp) {
			if(sp != (uint32_t)-1) el += ep - sp;
            sp = a[k].qs;
            ep = a[k].qe;
		} else {
			sp = MIN(sp, a[k].qs);
		}
	}
	if(sp != (uint32_t)-1) el += ep - sp;
	iel = (tte - tts) - el;
	// fprintf(stderr, "[M::%s] el:%u, iel:%u\n", __func__, el, iel);
	if((iel == 0) || (iel <= ((tte - tts)*trans_thres))) return 1;
	return 0;
}

uint32_t ff_chain(kv_ul_ov_t *idx, int64_t qlen, float cov_rate, float trans_thres, ul_ov_t *a, 
overlap_region_alloc* olist, haplotype_evdience_alloc *hap, const ul_idx_t *uref, double diff_ec_ul, int64_t winLen, 
void *km)
{
	if(idx->n <= 0) return 0;
	ul_ov_t *m = &(idx->a[idx->n-1]); //largest chain
	// fprintf(stderr, "[M::%s] m->score:%u, m->qs:%u, m->qe:%u, chain_n:%u\n", __func__, m->qn, m->qs, m->qe, m->te-m->ts);
	if((m->qe-m->qs) <= (qlen*cov_rate)) return 0;
	if(check_trans_rate(a+m->ts, m->te-m->ts, trans_thres)) return 1;
	if(olist && hap && uref) {
		int64_t idx_n = idx->n, z, i, het_n, resc_tk = 0, f = 0; 
		uint64_t si; ma_utg_t *u = NULL;
		for (z = m->ts; z < m->te; z++) {
			if(a[z].el) {
				kv_push_km(km, ul_ov_t, *idx, a[z]);
			} else {
				i = a[z].qn; si = 0;
				het_n = update_ava_het_site(hap, i, &si, NULL, 1);
				assert(het_n > 0 && olist->list[i].is_match == 2);
				u = &(uref->ug->u.a[olist->list[i].y_id]);
				if(u->n > 1) {
					resc_tk += rescue_trans_ul_chains(uref, &(olist->list[i]), hap->list+si, het_n, u, 
            		idx, diff_ec_ul, winLen, 0, NULL, km);
				}
			}
		}

		if(resc_tk) {
			radix_sort_ul_ov_srt_qe(idx->a+idx_n, idx->a+idx->n);
			f = check_trans_rate(idx->a+idx_n, idx->n-idx_n, trans_thres);
		}
		idx->n = idx_n;
		return f;
	}
	return 0;
}

void dump_chain(kv_ul_ov_t *des, ul_ov_t *src, ul_ov_t *chain, void *km)
{
	///note: dump results to <des> may change <chain>, so we should save <chain> in advance
	uint64_t beg = chain->ts, occ = chain->te - chain->ts;
	kv_resize_km(km, ul_ov_t, *des, occ); des->n = occ;
	memcpy(des->a, src + beg, occ*sizeof((*src)));
}

int64_t dedup_sort_contains(ul_ov_t *a, int64_t a_n, ul_contain *ct, const ug_opt_t *uopt)
{
	int64_t k, l, ci; ul_ov_t *z = NULL;
	for (k = 0; k < a_n; k++) {
		z = &(a[k]);
		if(z->tn&((uint32_t)(0x80000000))) continue;///contained alignment
		if(ct->is_c.a[z->tn] == 0) continue;

		for (ci = k+1; ci < a_n; ci++) {
			if(a[ci].qe > z->qe + G_CHAIN_INDEL) break;///128 is for indel
			if(a[ci].qn == (uint32_t)-1) continue;
			if(!(a[ci].tn&((uint32_t)(0x80000000)))) continue;
			if(z->qs <= a[ci].qs + G_CHAIN_INDEL && z->qe + G_CHAIN_INDEL >= a[ci].qe) {
				if(check_contain_pair(uopt, (a[ci].tn<<1)>>1, z->tn,  1)) {
					a[ci].qn = (uint32_t)-1;
				}
			}
		}

		for (ci = k-1; ci >= 0; ci--) {
			if(a[ci].qe + G_CHAIN_INDEL <= z->qs) break;
			if(a[ci].qn == (uint32_t)-1) continue;
			if(!(a[ci].tn&((uint32_t)(0x80000000)))) continue;
			if(z->qs <= a[ci].qs + G_CHAIN_INDEL && z->qe + G_CHAIN_INDEL >= a[ci].qe) {
				if(check_contain_pair(uopt, (a[ci].tn<<1)>>1, z->tn,  1)) {
					a[ci].qn = (uint32_t)-1;
				}
			}
		}
	}
	for (k = l = 0; k < a_n; k++) {
		if(a[k].qn == (uint32_t)-1) continue;
		if(k != l) a[l] = a[k];
		a[l].tn <<= 1; a[l].tn >>= 1;
		++l;
	}
	return l;
}

void ins_merge_ul_ov(kv_ul_ov_t *idx, int64_t idx_s, int64_t idx_e, ul_ov_t q)
{	
	int64_t k, ii, s = -1, e = -1, ovlp = 0, qs = q.qs, qe = q.qe;
	for (k = idx_s, ii = -1; k < idx_e; k++) {
		if(ii == -1 && q.qs > idx->a[k].qs) ii = k;
		if(((int64_t)(idx->a[k].qs)) >= e) {
			if(s >= 0 && e >= 0) {
				ovlp += ((MIN(e, qe) > MAX(s, qs))?(MIN(e, qe) - MAX(s, qs)):0);
			}
			s = idx->a[k].qs; e = idx->a[k].qe;
		} else {
			if(e < ((int64_t)(idx->a[k].qe))) e = idx->a[k].qe;
		}
	}

	if(s >= 0 && e >= 0) {
		ovlp += ((MIN(e, qe) > MAX(s, qs))?(MIN(e, qe) - MAX(s, qs)):0);
	}
}

void dump_all_chain(kv_ul_ov_t *idx, kv_ul_ov_t *ax, int64_t ax_new_occ, int64_t qlen, float primary_cov_rate, float primary_score_rate) {
	if(idx->n <= 0) return;
	ul_ov_t *m = &(idx->a[idx->n-1]); //largest chain
	ul_ov_t *a = ax->a + ax->n; int64_t k, i, z, l, idx_n = idx->n; uint64_t ovlp;
	if((m->qe-m->qs) > (qlen*primary_cov_rate)) { ///found a primary chain
		for (k = m->ts, l = 0; k < m->te; k++) {
			a[l] = a[k]; a[l].tn |= ((uint32_t)(0x80000000));
			l++;
		}
		ax->n += l;
	} else {
		for (k = idx_n-1; k >= 0; k--) {
			for (i = idx_n-1; i > k; i--) {
				if(idx->a[i].qn == (uint32_t)-1) continue;///just remove totally contained alignments
				if(idx->a[k].qn > idx->a[i].qn*primary_score_rate) continue;///consider score
				ovlp = ((MIN(idx->a[k].qe, idx->a[i].qe) > MAX(idx->a[k].qs, idx->a[i].qs))?
								(MIN(idx->a[k].qe, idx->a[i].qe) - MAX(idx->a[k].qs, idx->a[i].qs)):0);
				if(ovlp > ((idx->a[k].qe-idx->a[k].qs)*primary_cov_rate)) {
					for (z = idx->a[k].ts; z < idx->a[k].te; z++) a[z].el = 1;
					idx->a[k].qn = (uint32_t)-1;
					break;
				}
			}
			// ins_merge_ul_ov(idx, idx_n, idx->n, idx->a[k]);
		}
		for (k = 0, l = 0; k < ax_new_occ; k++) {
			if(a[k].el) continue;
			a[l] = a[k]; l++;
		}
		radix_sort_ul_ov_srt_qe(a, a + l);
		ax->n += l;
	}
}

int64_t dump_all_chain_simple(kv_ul_ov_t *idx, kv_ul_ov_t *ax, int64_t ax_new_occ, int64_t qlen, 
float primary_cov_rate, float fragement_cov_rate, float primary_fragment_cov_rate, 
float primary_fragment_second_score_rate, float trans_thres, uint64_t mini_primary_fragment_len) 
{
	if(idx->n <= 0) return 0;
	ul_ov_t *m = &(idx->a[idx->n-1]); //largest chain
	ul_ov_t *a = ax->a + ax->n; int64_t k, z, l, idx_n = idx->n, ovlp, om, ok, ff = 0;
	// fprintf(stderr, "[M::%s] m->score:%u, m->qs:%u, m->qe:%u, chain_n:%u\n", __func__, m->qn, m->qs, m->qe, m->te-m->ts);
	if(((m->qe-m->qs) > (qlen*primary_cov_rate)) && 
						(check_trans_rate(a+m->ts, m->te-m->ts, trans_thres))) { ///found a primary chain
		for (k = m->ts, l = 0; k < m->te; k++) {
			a[l] = a[k]; a[l].tn |= ((uint32_t)(0x80000000)); a[l].el = 1;
			l++;
		}
		ax->n += l; ff = 1;
	} else {
		if((((m->qe-m->qs) > (qlen*primary_fragment_cov_rate)) || ((m->qe - m->qs) > mini_primary_fragment_len)) 
		&& (check_trans_rate(a+m->ts, m->te-m->ts, trans_thres))) {
			om = m->qe - m->qs;
			for (k = 0; k < idx_n-1; k++) {
				ovlp = ((MIN((m->qe), (idx->a[k].qe)) > MAX((m->qs), (idx->a[k].qs)))? 
										(MIN((m->qe), (idx->a[k].qe)) - MAX((m->qs), (idx->a[k].qs))):0);
				if(ovlp == 0) continue;
				ok = idx->a[k].qe - idx->a[k].qs;
				if(ok > om ) ok = om;
				if((ovlp > ok*0.1/**0.25**/) && idx->a[k].qn > (m->qn*primary_fragment_second_score_rate)) break;
			}

			if(k >= idx_n-1) {
				for (k = m->ts; k < m->te; k++) {
					// if(a[k].el) a[k].tn |= ((uint32_t)(0x80000000));
					a[k].tn |= ((uint32_t)(0x80000000));
				}
				ff = 1;
			}
		}

		radix_sort_ul_ov_srt_qe(idx->a, idx->a + idx->n);
		for (k = 0; k < idx_n; k++) {
			if(k < idx_n-1 && idx->a[k].qe > idx->a[k+1].qs) break;//not one chain
			if((idx->a[k].qe - idx->a[k].qs) > (qlen*fragement_cov_rate)) {///large enough fragements
				if(!check_trans_rate(a+idx->a[k].ts, idx->a[k].te-idx->a[k].ts, trans_thres)) break;
			}
		}

		if(k == idx_n) {///only if there is a clear chain (with holes)
			for (k = 0; k < idx_n; k++) {
				if((idx->a[k].qe - idx->a[k].qs) <= (qlen*fragement_cov_rate)) continue;
				for (z = idx->a[k].ts; z < idx->a[k].te; z++) {
					// if(a[z].el) a[z].tn |= ((uint32_t)(0x80000000));
					a[z].tn |= ((uint32_t)(0x80000000));
				}
				ff = 1;
			}
		} 
		/**
		for (k = 0, l = 0; k < ax_new_occ; k++) {
			if(!(a[k].el)) continue;
			a[l] = a[k]; l++;
		}

		radix_sort_ul_ov_srt_qe(a, a + l);
		ax->n += l;
		**/
		for (k = 0, l = 0; k < ax_new_occ; k++) {
			if(a[k].el || (a[k].tn&((uint32_t)(0x80000000)))) {
				a[l] = a[k]; l++;
			}
		}
		ax_new_occ = l;

		radix_sort_ul_ov_srt_qe(a, a + ax_new_occ);
		for (k = 1, l = 0; k <= ax_new_occ; k++) {
			if (k == ax_new_occ || a[k].qe != a[l].qe) {
				if(k - l > 1) radix_sort_ul_ov_srt_qs(a+l, a+k);
				l = k;
			}
		}
		ax->n += ax_new_occ;
	}
	return ff;
}

void save_tmp_chains(ul_ov_t *idx_a, uint64_t idx_n, uint64_t *idx_buf_0, uint64_t *idx_buf_1, ul_ov_t *cc_a, uint64_t cc_n, uint64_t *cc_buf)
{
	uint64_t k;
	for (k = 0; k < idx_n; k++) ;
}

void debug_reverse_chain(ul_ov_t *a, int64_t a_n)
{
	int64_t rev_n = a_n>>1, i; ul_ov_t rev_t;
	for (i = 0; i < rev_n; i++) {
		rev_t = a[i]; a[i] = a[a_n-i-1]; a[a_n-i-1] = rev_t;
	}
}

uint32_t quick_primary_assgin(const ul_idx_t *uref, const ug_opt_t *uopt, int64_t bw, double diff_ec_ul, ul_ov_t *a, int64_t a_n)
{
	int64_t k, qo, share, f = 1; ul_ov_t *li, *lk; uint32_t li_v, lk_v;
	for (k = a_n - 1, li = NULL; k >= 0; k--) {
		lk = &(a[k]); lk_v = (lk->tn<<1)|lk->rev; lk->sec = SEC_MODE; 
		if(!(lk->tn&((uint32_t)(0x80000000)))) continue;
		if(li && lk->qe > li->qs) { ///lk is overlapped with li
			li->tn <<= 1; li->tn >>= 1; lk->tn <<= 1; lk->tn >>= 1;
			qo = infer_rovlp(li, lk, NULL, NULL, NULL, NULL);
			li->tn |= ((uint32_t)(0x80000000)); lk->tn |= ((uint32_t)(0x80000000));
			if(qo && li_v != lk_v && get_ecov_adv(uref, uopt, li_v^1, lk_v^1, bw, diff_ec_ul, qo, 1, &share)) {
				li->sec = k; f++; ///the end of a chain is a cis overlap
			} else {
				f = 0;
				break;
			}
		}
		li = lk; li_v = lk_v; 
	}
	return f;
}

void assgin_primary_chains(ul_ov_t *a, int64_t a_n, int64_t is_srt, const ul_idx_t *uref, const ug_opt_t *uopt, 
uint64_t *track, uint64_t *srt, int64_t bw, double diff_ec_ul, int64_t qlen, int64_t is_ungap)
{
	if(a_n == 0) return;
	int64_t i, j, k, mm_ovlp, x, csc, mm_sc, mm_idx, qo, share, sc; 
	uint32_t li_v, lj_v; ul_ov_t *li = NULL, *lj = NULL;
	if(is_srt) {
		radix_sort_ul_ov_srt_qe(a, a + a_n);
		for (i = 1, j = 0; i <= a_n; i++) {
			if (i == a_n || a[i].qe != a[j].qe) {
				if(i - j > 1) {
					radix_sort_ul_ov_srt_qs(a+j, a+i);
				}
				j = i;
			}
		}
	}

	if(uref && uopt && track && srt) {
		if(quick_primary_assgin(uref, uopt, bw, diff_ec_ul, a, a_n) == 0) {
			for (i = 0; i < a_n; ++i) {
				li = &(a[i]); li_v = (li->tn<<1)|li->rev; li->sec = SEC_MODE; 
				mm_sc = csc = 1; mm_idx = -1; 
				if(li->tn&((uint32_t)(0x80000000))) {
					mm_ovlp = max_ovlp_src(uopt, li_v^1);
					x = (li->qs + mm_ovlp)*diff_ec_ul; 
					if(x < bw) x = bw; 
					x += li->qs + mm_ovlp;
					if (x > qlen+1) x = qlen+1;
					x = find_ul_ov_max(i, a, x+G_CHAIN_INDEL);
					
					for (j = x; j >= 0; --j) {
						lj = &(a[j]); lj_v = (lj->tn<<1)|lj->rev;
						if(lj->qe+G_CHAIN_INDEL <= li->qs) break;//even this pair has a overlap, its length will be very small; just ignore
						if(is_ungap && (lj->qs >= li->qs+G_CHAIN_INDEL)) continue; // lj is contained in li on the query coordinate; 128 for indel offset
						if(!(lj->tn&((uint32_t)(0x80000000)))) continue;
						li->tn <<= 1; li->tn >>= 1; lj->tn <<= 1; lj->tn >>= 1;
						qo = infer_rovlp(li, lj, NULL, NULL, NULL, NULL);
						li->tn |= ((uint32_t)(0x80000000)); lj->tn |= ((uint32_t)(0x80000000));
						if(qo && li_v != lj_v && get_ecov_adv(uref, uopt, li_v^1, lj_v^1, bw, diff_ec_ul, qo, 1, &share)) {
							sc = csc + pop_sc(track[j]);
							if(sc > mm_sc) mm_sc = sc, mm_idx = j;
						}
					}
				}
				track[i] = push_sc_pre(mm_sc, mm_idx);
				srt[i] = track[i]>>32; srt[i] <<= 32; srt[i] |= i;
				// fprintf(stderr, "[M::] i->%ld; mm_idx->%ld\n", i, mm_idx);
			}

			radix_sort_gfa64(srt, srt+a_n); 
			for (k = a_n-1; k >= 0; --k) {
				i = (uint32_t)srt[k];
				// if(i < 0 || i >= a_n) fprintf(stderr, "sbsbsbsbsbsb, k->%ld, i->%ld, a_n->%ld\n", k, i, a_n);
				if(a[i].el && (a[i].tn&((uint32_t)(0x80000000)))) {
					for (; i >= 0 && (track[i]&((uint64_t)0x80000000)) == 0;) {
						track[i] |= ((uint64_t)0x80000000); j = i;
						i = pop_pre(track[i]); 
						if(i >= 0) {
							// if(i == j) fprintf(stderr, "sb\n");
							a[j].sec = i; 
						}
					}
				}
			}
		}
	}	
}

int64_t gl_chain_refine_advance(overlap_region_alloc* olist, Correct_dumy* dumy, haplotype_evdience_alloc *hap, glchain_t *ll, st_mt_t *sps, const ul_idx_t *uref, double diff_ec_ul, int64_t winLen, int64_t qlen, const ug_opt_t *uopt, 
int64_t debug_i, void *km)
{
	// ll->tk.n = ll->lo.n = 0;
	kv_ul_ov_t *idx = &(ll->lo);
	ul_contain *ct = uref->ct;
	uint64_t o2 = gl_chain_gen(olist, uref, idx, 0, hap, km);
	if(idx->n == 0) return 0;
	// fprintf(stderr, "[M::%s] qlen:%ld, idx->n:%u\n", __func__, qlen, (uint32_t)idx->n);
	uint64_t k, an, cn, si = 0, ei = 0, resc = 0, resc_tk = 0, tk_pl = 0, f = 0, occ = 0, cis_occ = 0, t_cis = 0;
	ma_utg_t *u = NULL; overlap_region *o = NULL;

	kv_resize_km(km, uint64_t, ll->srt.a, idx->n);
    kv_resize_km(km, uint64_t, hap->snp_srt, idx->n);
	kv_resize_km(km, ul_ov_t, ll->tk, ll->tk.n+idx->n);
	///note: there are three rounds of gl_chain_advance()
	///the first two rounds could reuse dumy->overlapID. But for the last round, dumy->overlapID is not long enough 
	occ = gl_chain_advance(idx, ll->tk.a+ll->tk.n, uref, uopt, G_CHAIN_BW, diff_ec_ul, qlen, UG_SKIP, dumy->overlapID, ll->srt.a.a, hap->snp_srt.a, G_CHAIN_TRANS_WEIGHT, 0, NULL, uref->ug, debug_i, km);
	if(occ) {
		if(ff_chain(idx, qlen, P_CHAIN_COV, G_CHAIN_TRANS_RATE, ll->tk.a+ll->tk.n, NULL, NULL, NULL, diff_ec_ul, winLen, km)) {
			f = 1; //dump_chain(idx, ll->tk.a+ll->tk.n, &(idx->a[idx->n-1]), km); 
			for (k = idx->a[idx->n-1].ts; k < idx->a[idx->n-1].te; k++) {
				olist->list[ll->tk.a[ll->tk.n+k].qn].x_pos_strand = 1;
			}
		} else if(o2) {///means there are trans overlaps
			gl_chain_gen(olist, uref, idx, 1, hap, km);
			kv_resize_km(km, uint64_t, ll->srt.a, idx->n);
    		kv_resize_km(km, uint64_t, hap->snp_srt, idx->n);
			kv_resize_km(km, ul_ov_t, ll->tk, ll->tk.n+idx->n);
			///chain all U-matches
			occ = gl_chain_advance(idx, ll->tk.a+ll->tk.n, uref, uopt, G_CHAIN_BW, diff_ec_ul, qlen, UG_SKIP, dumy->overlapID, ll->srt.a.a, hap->snp_srt.a, G_CHAIN_TRANS_WEIGHT, 0, NULL, uref->ug, debug_i, km);
			if(ff_chain(idx, qlen, P_CHAIN_COV, G_CHAIN_TRANS_RATE, ll->tk.a+ll->tk.n, olist, hap, uref, diff_ec_ul, winLen, km)) {
				f = 1; //dump_chain(idx, ll->tk.a+ll->tk.n, &(idx->a[idx->n-1]), km); 
				for (k = idx->a[idx->n-1].ts; k < idx->a[idx->n-1].te; k++) {
					olist->list[ll->tk.a[ll->tk.n+k].qn].x_pos_strand = 1;
				}
			}
		}
	}

	if(!f) {///if f == 1, only dump primary chain; otherwise dump all chains
		///we can save all data to buffer like ll->srt.a in advance; in case we don't need third round of chaining
		///means no trans overlaps, no need to do third round of chaining
		// if(!o2) {
		// 	;
		// }
		/**
		uint64_t z;
		for (k = 0; k < idx->n; k++) {
			if(check_trans_rate(ll->tk.a+ll->tk.n+idx->a[k].ts, idx->a[k].te-idx->a[k].ts, G_CHAIN_TRANS_RATE)) {
				for (z = idx->a[k].ts; z < idx->a[k].te; z++) {
					olist->list[ll->tk.a[ll->tk.n+z].qn].x_pos_strand = 1;
				}
			} else {
				for (z = idx->a[k].ts; z < idx->a[k].te; z++) {
					if(!(ll->tk.a[ll->tk.n+z].el)) continue;
					olist->list[ll->tk.a[ll->tk.n+z].qn].x_pos_strand = 1;
				}
			}	
		}
		**/
		for (k = 0; k < occ; k++) {
			olist->list[ll->tk.a[ll->tk.n+k].qn].x_pos_strand = 1;
		}
		// kv_resize_km(km, ul_ov_t, *idx, occ); idx->n = occ;
		// memcpy(idx->a, ll->tk.a+ll->tk.n, occ*sizeof((*(idx->a))));
	}
	// for (k = 0; k < idx->n; k++) olist->list[idx->a[k].qn].x_pos_strand = 1;
		
	for (k = 0, idx->n = 0, tk_pl = ll->tk.n, t_cis = 0; k < olist->length; k++) {
		o = &(olist->list[k]);
		///if f == 1, no matter 
		if(o->x_pos_strand && (f || o->is_match == 1)){
			u = &(uref->ug->u.a[o->y_id]);///overlaped reads
			resc_tk += rescue_trans_ul_chains(uref, o, NULL, 0, u, &(ll->tk), -1, -1, 0, NULL, km);
		} else if((!f) && o->is_match == 2) {
			// fprintf(stderr, "###[M::%s] # k:%lu, # o->y_id:%u\n", __func__, k, o->y_id);
			an = update_ava_het_site(hap, k, &si, &ei, 1);
			// if(an != get_het_site(hap, k)) fprintf(stderr, "an->%lu, get_het_site->%lu\n", an, get_het_site(hap, k));
			assert(an > 0);

			cn = ((uint32_t)(ct->idx.a[o->y_id]));		
			if(cn > 0) {
				resc += rescue_contain_ul_chains(uref, o, hap->list+si, an, 
				ct->rids.a + ((ct->idx.a[o->y_id])>>32), cn, idx, diff_ec_ul, winLen, 0, km);
			}

			u = &(uref->ug->u.a[o->y_id]);
			if(u->n > 1 || o->x_pos_strand) {///no redundant items here
				resc_tk += rescue_trans_ul_chains(uref, o, hap->list+si, an, u, 
				&(ll->tk), diff_ec_ul, winLen, o->x_pos_strand, &cis_occ, km);
				t_cis += cis_occ;
			}

			si = ei;
		} 
	}

	assert(ll->tk.n == resc_tk+tk_pl);
	assert(idx->n == resc);
	if(f) assert(resc==0);

		
	if(!f) {///dedup contained alignments
		if(idx->n) {///if some contained alignments have been rescued
			radix_sort_ul_ov_srt_tn(idx->a, idx->a + idx->n);
			idx->n = dedup_sort_ul_ov_t(idx->a, idx->n);///different trans alignments may have the same contained alignment
		}
		resc = idx->n;
		// fprintf(stderr, "***[M::%s] # contain:%lu, # non-contain:%lu\n", __func__, resc, (uint64_t)(ll->tk.n-tk_pl));
		for (k = tk_pl; k < ll->tk.n; k++) {///dump all non-contained reads
			kv_push_km(km, ul_ov_t, *idx, ll->tk.a[k]);
			if(idx->a[idx->n-1].tn&((uint32_t)(0x80000000))) {
				idx->a[idx->n-1].tn -= ((uint32_t)(0x80000000));
			}
		}
		ll->tk.n = tk_pl;

		radix_sort_ul_ov_srt_qe(idx->a, idx->a + idx->n);
		if(resc) {///need to dedup contained alignment again
			// fprintf(stderr, "[M::%s] idx->n:%lu, resc:%lu\n", __func__, (uint64_t)idx->n, resc);
			idx->n = dedup_sort_contains(idx->a, idx->n, ct, uopt);
		}
		
		///note: need sps for third round of gl_chain_advance() as dumy->overlapID might be not long enough
		kv_resize_km(km, uint64_t, *sps, idx->n);
		kv_resize_km(km, uint64_t, ll->srt.a, idx->n);
		kv_resize_km(km, uint64_t, hap->snp_srt, idx->n);
		kv_resize_km(km, ul_ov_t, ll->tk, ll->tk.n+idx->n);
		occ = gl_chain_advance(idx, ll->tk.a+ll->tk.n, uref, uopt, G_CHAIN_BW, diff_ec_ul, qlen, UG_SKIP, sps->a, ll->srt.a.a, hap->snp_srt.a, G_CHAIN_TRANS_WEIGHT, 1, &R_INF, NULL, debug_i, km);
		// fprintf(stderr, "***[M::%s] ll->tk.n:%u, occ:%lu\n", __func__, (uint32_t)ll->tk.n, occ);
		f = dump_all_chain_simple(idx, &(ll->tk), occ, qlen, P_CHAIN_COV, P_FRAGEMENT_CHAIN_COV, 
		P_FRAGEMENT_PRIMARY_CHAIN_COV, 0.1/**P_FRAGEMENT_PRIMARY_SECOND_COV**/, G_CHAIN_TRANS_RATE, PRIMARY_UL_CHAIN_MIN);
		// fprintf(stderr, ">>>[M::%s] ll->tk.n:%u\n", __func__, (uint32_t)ll->tk.n);
		// dump_all_chain(idx, &(ll->tk), occ, qlen, P_CHAIN_COV, P_CHAIN_SCORE);
	} else {
		///for primary chain, each element x: (x->tn & (uint32_t)(0x80000000))
		assgin_primary_chains(ll->tk.a+tk_pl, ll->tk.n-tk_pl, 1, NULL, NULL, NULL, NULL, G_CHAIN_BW, diff_ec_ul, qlen, 1);
		// radix_sort_ul_ov_srt_qe(ll->tk.a+tk_pl, ll->tk.a+ll->tk.n);
	} 

	///if f == 0, results have already been sorted by qe|qs
	if(f) {
		kv_resize_km(km, uint64_t, ll->srt.a, ll->tk.n-tk_pl); kv_resize_km(km, uint64_t, hap->snp_srt, ll->tk.n-tk_pl);
		assgin_primary_chains(ll->tk.a+tk_pl, ll->tk.n-tk_pl, 0, uref, uopt, ll->srt.a.a, hap->snp_srt.a, G_CHAIN_BW, diff_ec_ul, qlen, 1);
	}
	// debug_reverse_chain(ll->tk.a+tk_pl, ll->tk.n-tk_pl);
	/**
	if(idx->n > 0) {
		an = infer_read_ovlp(uref, olist, idx	, &(ll->tk), diff_ec_ul, winLen, uopt, ct, km);
		// if(an) fill_edge_weight(ll->tk.a+ll->tk.n-an, an, uopt, G_CHAIN_BW, diff_ec_ul, qlen);
	}
	**/
	return 1;
}

int64_t g_adjacent_dis(const asg_t *g, uint32_t v, uint32_t w)
{
	uint32_t nv, i; asg_arc_t *av = NULL;
	nv = asg_arc_n(g, v); av = asg_arc_a(g, v);
	for (i = 0; i < nv; i++) {
        if(av[i].del || av[i].v != w) continue;
		return (uint32_t)av[i].ul;
	}
	return -1;
}

void get_r_offset(ma_ug_t *ug, mg_lchain_t *x, int64_t *rs, int64_t *re, int64_t *qs, int64_t *qe)
{
	if(qs) *qs = x->qs; if(qe) *qe = x->qe;
	if(!(x->score&1)) {
		if(rs) *rs = x->rs + x->off;
		if(re) *re = x->re + x->off;
	} else {
		if(rs) *rs = x->off + ug->g->seq[x->score>>1].len - x->re;
		if(re) *re = x->off + ug->g->seq[x->score>>1].len - x->rs;
	}
}

void get_u_offset(ma_ug_t *ug, mg_lchain_t *x, int64_t *rs, int64_t *re, int64_t *qs, int64_t *qe)
{
	if(qs) *qs = x->qs; if(qe) *qe = x->qe;
	if(!(x->v&1)) {
		if(rs) *rs = x->rs + x->off;
		if(re) *re = x->re + x->off;
	} else {
		if(rs) *rs = x->off + ug->g->seq[x->v>>1].len - x->re;
		if(re) *re = x->off + ug->g->seq[x->v>>1].len - x->rs;
	}
}

void l2g_chain(const ul_idx_t *uref, kv_ul_ov_t *lidx, vec_mg_lchain_t *res)
{
	uint64_t k; 
	res->n = 0; kv_resize(mg_lchain_t, *res, lidx->n); res->n = lidx->n;
	for (k = 0; k < lidx->n; k++) {
		memset(&(res->a[k]), 0, sizeof(res->a[k]));
		res->a[k].v = (lidx->a[k].tn<<1)|(lidx->a[k].rev);
		res->a[k].off = lidx->a[k].qn; res->a[k].score = lidx->a[k].sec;
		res->a[k].qs = lidx->a[k].qs; res->a[k].qe = lidx->a[k].qe;
		res->a[k].rs = lidx->a[k].ts; res->a[k].re = lidx->a[k].te; 
		if(lidx->a[k].el) {
			res->a[k].score = retrieve_u_cov_region(uref, lidx->a[k].tn, 0, lidx->a[k].ts, lidx->a[k].te, NULL);
		}
	}
}

int64_t l2g_res_chain(ma_ug_t *ug, ul_ov_t *a, uint64_t a_n, vec_mg_lchain_t *gchains, double diff_rate)
{
	if(a_n <= 0) return 0;
	uint64_t k, m; int64_t l, rs, re, qs, qe, dq, dr, dif, mm; a_n++; asg_t *g = ug->g;
	gchains->n = 0; kv_resize(mg_lchain_t, *gchains, a_n); gchains->n = a_n;
	memset(&(gchains->a[0]), 0, sizeof(gchains->a[0])); 
	gchains->a[0].cnt = a_n - 1; gchains->a[0].v = (uint32_t)-1;
	for (k = 1, m = 0, l = 0; k < a_n; k++, m++) {
		memset(&(gchains->a[k]), 0, sizeof(gchains->a[k]));
		gchains->a[k].v = (a[m].tn<<1)|(a[m].rev); gchains->a[k].dist_pre = -1;
		gchains->a[k].off = a[m].qn; gchains->a[k].score = a[m].sec;
		gchains->a[k].qs = a[m].qs; gchains->a[k].qe = a[m].qe;
		gchains->a[k].rs = a[m].ts; gchains->a[k].re = a[m].te; 
		if(k > 1) {
			gchains->a[k-1].dist_pre = g_adjacent_dis(g, gchains->a[k].v^1, gchains->a[k-1].v^1);
			assert(gchains->a[k-1].dist_pre >= 0);
			l += g->seq[gchains->a[k-1].v>>1].len + gchains->a[k-1].dist_pre;
		}
	}

	if(diff_rate < 0) return 1;

	mg_lchain_t s = gchains->a[1], e = gchains->a[a_n-1];
	s.off = 0; 
	l -= (int64_t)g->seq[gchains->a[a_n-1].v>>1].len;
	if(l < 0) l = 0; 
	e.off = l;

	get_u_offset(ug, &s, &rs, NULL, &qs, NULL); get_u_offset(ug, &e, NULL, &re, NULL, &qe);
	dq = qe - qs; dr = re - rs;
	dif = (dq>dr? dq-dr:dr-dq);
	mm = MAX(dq, dr); mm *= diff_rate; 
	if(dif <= mm) return 1;
	return 0;
}

int64_t check_elen_gchain(ul_ov_t *a, int64_t a_n, float trans_thres)
{
	uint32_t sp_e, ep_e, ts, te, tl = 0, el = 0, iel = 0; 
	sp_e = ep_e = ts = te = (uint32_t)-1;
	int64_t k;
	for (k = a_n-1; k >= 0; k--) {
		if(ts == (uint32_t)-1 || a[k].qe <= ts) {
			if(ts != (uint32_t)-1) tl += te - ts;
            ts = a[k].qs; te = a[k].qe;
		} else {
			ts = MIN(ts, a[k].qs);
		}
		if(!(a[k].el)) continue;

		if(sp_e == (uint32_t)-1 || a[k].qe <= sp_e) {
			if(sp_e != (uint32_t)-1) el += ep_e - sp_e;
            sp_e = a[k].qs;
            ep_e = a[k].qe;
		} else {
			sp_e = MIN(sp_e, a[k].qs);
		}
	}
	if(ts != (uint32_t)-1) tl += te - ts;
	if(sp_e != (uint32_t)-1) el += ep_e - sp_e;

	iel = tl - el;
	// fprintf(stderr, "[M::%s] el:%u, iel:%u\n", __func__, el, iel);
	if((iel == 0) || (iel <= (tl*trans_thres))) return 1;
	return 0;
}

int64_t ds_check_vec_mg_lchain_t(mg_lchain_t *a, int64_t a_n, kv_ul_ov_t *buf, overlap_region_alloc* olist, 
haplotype_evdience_alloc *hap, const ul_idx_t *uref, double diff_ec_ul, int64_t winLen, float trans_thres)
{
	if(a_n <= 0) return 0;
	int64_t k, resc_tk, het_n; ul_ov_t *p; uint64_t si; ma_utg_t *u = NULL; buf->n = 0;
	for (k = 0; k < a_n; k++) {
		if(a[k].off < 0) continue;
		kv_pushp(ul_ov_t, *buf, &p);
		p->qs = a[k].qs; p->qe = a[k].qe; p->el = 1;
		if(a[k].score < 0) p->el = 0;
	}
	if(buf->n <= 0) return 0;
	if(check_elen_gchain(buf->a, buf->n, trans_thres)) return 1;

	buf->n = 0; resc_tk = 0;
	for (k = 0; k < a_n; k++) {
		if(a[k].off < 0) continue;
		if(a[k].score >= 0) {
			kv_pushp(ul_ov_t, *buf, &p);
			p->qs = a[k].qs; p->qe = a[k].qe; p->el = 1;
		} else {
			si = 0;
			het_n = update_ava_het_site(hap, a[k].off, &si, NULL, 1);
			assert(het_n > 0 && olist->list[a[k].off].is_match == 2);
			u = &(uref->ug->u.a[olist->list[a[k].off].y_id]);
			if(u->n > 1) {
				resc_tk += rescue_trans_ul_chains(uref, &(olist->list[a[k].off]), hap->list+si, het_n, u, 
				buf, diff_ec_ul, winLen, 0, NULL, NULL);
			}
		}		
	}

	if(resc_tk) {
		radix_sort_ul_ov_srt_qe(buf->a, buf->a+buf->n);
		if(check_elen_gchain(buf->a, buf->n, trans_thres)) return 1;
	}

	return 0;
}

int64_t check_trans_rate_gap(vec_mg_lchain_t *uc, kv_ul_ov_t *buf, overlap_region_alloc* olist, 
haplotype_evdience_alloc *hap, const ul_idx_t *uref, double diff_ec_ul, int64_t winLen, float trans_thres)
{
	if(uc->n <= 0) return 0;
	int64_t k, m, ucn = uc->n, k_cnt, z; mg_lchain_t *ix; buf->n = 0; 
	for (k = m = 0; k < ucn; k += k_cnt) {
		ix = &(uc->a[k]); assert(ix->v == (uint32_t)-1); k_cnt = ix->cnt + 1;
		if(ds_check_vec_mg_lchain_t(uc->a + k + 1, ix->cnt, buf, olist, hap, uref, diff_ec_ul, winLen, trans_thres)) {
			if(m == k) {
				m += k_cnt; 
			} else {
				for (z = 0; z < k_cnt; z++) uc->a[m++] = uc->a[k+z];
			}
		}
	}
	uc->n = m;
	if(uc->n) return 1;
	return 0;
}

int64_t hc_gchain1_dp(void *km, const ul_idx_t *uref, const ma_ug_t *ug, vec_mg_lchain_t *lc, vec_mg_lchain_t *sw, vec_mg_path_dst_t *dst, vec_sp_node_t *out, vec_mg_pathv_t *path,
int64_t qlen, const ug_opt_t *uopt, int64_t bw, double diff_thre, double ng_diff_thre, uint64_t *srt, st_mt_t *bf, int64_t *f, uint64_t *p, uint64_t *v);
uint32_t gen_max_gchain_adv(void *km, const ul_idx_t *uref, int64_t ulid, st_mt_t *idx, vec_mg_lchain_t *e, kv_ul_ov_t *raw_idx, kv_ul_ov_t *raw_chn,
int64_t qlen, float primary_cov_rate, float primary_fragment_cov_rate, float primary_fragment_second_score_rate, uint64_t mini_primary_fragment_len, 
const asg_t *g, st_mt_t *dst_done, vec_sp_node_t *out, vec_mg_pathv_t *res, uint64_t *b, vec_mg_lchain_t *gchains);



void debug_intermediate_chain(ma_ug_t *ug, mg_lchain_t *a, int64_t a_n, int64_t is_uovlp, int64_t debug_i)
{
	int64_t k; mg_lchain_t *p, *c;
	int64_t prs, pre, pqs, pqe, crs, cre, cqs, cqe;
	int64_t tot = 0, fal = 0;
	for (k = a_n-1; k >= 0; k--) {
		c = &(a[k]); 
		if(c->v!=(uint32_t)-1) {
			fprintf(stderr, "+(%ld) [M::utg%.6u%c::%c] qs:%d, qe:%d, rs:%d, re:%d, off:%d\n", k, (c->v>>1)+1, "lc"[ug->u.a[c->v>>1].circ], "+-"[c->v&1],
			c->qs, c->qe, c->rs, c->re, c->off);
		}
		
		
		if(c->hash_pre != (uint32_t)-1) {
			p = &(a[c->hash_pre]); tot++;
			if(is_uovlp) {
				get_u_offset(ug, p, &prs, &pre, &pqs, &pqe);
			} else {
				get_r_offset(ug, p, &prs, &pre, &pqs, &pqe);
			}
			if(is_uovlp) {
				get_u_offset(ug, c, &crs, &cre, &cqs, &cqe);
			} else {
				get_r_offset(ug, c, &crs, &cre, &cqs, &cqe);
			}

			fprintf(stderr, "cur (%ld) [M::utg%.6u%c::%c] qs:%ld, qe:%ld, rs:%ld, re:%ld, pidx:%u\n", 
			k, (c->v>>1)+1, "lc"[ug->u.a[c->v>>1].circ], "+-"[c->v&1], cqs, cqe, crs, cre, c->hash_pre);

			fprintf(stderr, "pre (%ld) [M::utg%.6u%c::%c] qs:%ld, qe:%ld, rs:%ld, re:%ld, pidx:%u\n", 
			k, (p->v>>1)+1, "lc"[ug->u.a[p->v>>1].circ], "+-"[p->v&1], pqs, pqe, prs, pre, p->hash_pre);
			
			if((!(prs<=crs&&pre<=cre&&pqs<=cqs&&pqe<=cqe)) || (!(prs<=pre&&pqs<=pqe&&crs<=cre&&cqs<=cqe))) {
				fal++;
				// fprintf(stderr, "[M::%s::a_n->%ld, k->%ld]p->dist_pre:%d, c->dist_pre:%d, c->pre_idx:%u\n", __func__, a_n, k, p->dist_pre, c->dist_pre, c->hash_pre);
				// fprintf(stderr, "[M::%s::]prs->%ld, pre->%ld, pqs->%ld, pqe->%ld, crs->%ld, cre->%ld, cqs->%ld, cqe->%ld\n", 
				// __func__, prs, pre, pqs, pqe, crs, cre, cqs, cqe);
			}
			assert(prs<=pre&&pqs<=pqe&&crs<=cre&&cqs<=cqe);
			// assert(prs<=crs&&pre<=cre&&pqs<=cqs&&pqe<=cqe);
		}
	}	
	if(fal) fprintf(stderr, "[M::%s::tot->%ld, fal->%ld] ulid->%ld\n", __func__, tot, fal, debug_i);
}

int64_t adjust_utg_chain_qoffset(uint64_t *r_srt, int64_t *r_pos, int64_t *q_pos, int64_t r_off)
{
	if(r_off < ((int64_t)(r_srt[0]>>32))) {
		// fprintf(stderr, "r_off:%ld, r_srt[0]:%ld\n", r_off, ((int64_t)(r_srt[0]>>32))); 
		return q_pos[(uint32_t)r_srt[0]]; ///have small chance
	}
	if(r_off >= ((int64_t)(r_srt[3]>>32))) {
		// fprintf(stderr, "r_off:%ld, r_srt[3]:%ld\n", r_off, ((int64_t)(r_srt[3]>>32))); 
		return q_pos[(uint32_t)r_srt[3]]; ///have small chance
	}
	int64_t k, rdis, qdis;
	for (k = 0; k < 3; k++) {
		if(r_off >= ((int64_t)(r_srt[k]>>32)) && r_off < ((int64_t)(r_srt[k+1]>>32))) break;
	}
	// if(k >= 3) {
	// 	if(rs >= ((int64_t)(r_srt[2]>>32)) && rs <= ((int64_t)(r_srt[3]>>32))) k = 2;
	// }
	if(k >= 3) {
		for (k = 0; k < 3; k++) {
			if(r_off >= ((int64_t)(r_srt[k]>>32)) && r_off <= ((int64_t)(r_srt[k+1]>>32))) break;
		}
	}
	assert(k < 3);

	rdis = r_pos[(uint32_t)r_srt[k+1]] - r_pos[(uint32_t)r_srt[k]]; 
	qdis = q_pos[(uint32_t)r_srt[k+1]] - q_pos[(uint32_t)r_srt[k]]; 
	return q_pos[(uint32_t)r_srt[k]] + get_offset_adjust(r_off-r_pos[(uint32_t)r_srt[k]], rdis, qdis);
}

void update_uovlp_chain_qse(ma_ug_t *ug, int64_t sidx, int64_t eidx, mg_lchain_t *a, int64_t a_n)
{
	// fprintf(stderr, "******[M::%s::] sidx:%ld, eidx:%ld\n", __func__, sidx, eidx);
	if(eidx - sidx <= 1) return;
	///for ug chains, sidx >= 0 && eidx < a_n
	assert(sidx>=0 && eidx<a_n);
	
	int64_t r_pos[4], q_pos[4];  uint64_t r_srt[4];
	int64_t i, rs, re;//qs or qe might be -1, while rs and re should >= 0
	if(sidx >= 0) {
		get_u_offset(ug, &(a[sidx]), &r_pos[0], &r_pos[1], &q_pos[0], &q_pos[1]);
	} else {
		get_u_offset(ug, &(a[0]), &r_pos[0], &r_pos[1], &q_pos[0], &q_pos[1]);
	}

	if(eidx < a_n) {
		get_u_offset(ug, &(a[eidx]), &r_pos[2], &r_pos[3], &q_pos[2], &q_pos[3]);
	} else {
		get_u_offset(ug, &(a[a_n-1]), &r_pos[2], &r_pos[3], &q_pos[2], &q_pos[3]);
	}
	// assert((left_q[0] >= 0 && left_q[1] >= 0) || (right_q[0] >= 0 && right_q[1] >= 0)); ///assert(re >= rs);
	///for ug chains, sidx >= 0 && eidx < a_n
	assert(q_pos[0] >= 0 && q_pos[1] >= 0 && q_pos[2] >= 0 && q_pos[3] >= 0);
	r_srt[0] = r_pos[0]; r_srt[0] <<= 32;
	r_srt[1] = r_pos[1]; r_srt[1] <<= 32; r_srt[1] += 1;
	r_srt[2] = r_pos[2]; r_srt[2] <<= 32; r_srt[2] += 2;
	r_srt[3] = r_pos[3]; r_srt[3] <<= 32; r_srt[3] += 3;
	radix_sort_gfa64(r_srt, r_srt + 4);
	// assert(r_pos[(uint32_t)r_srt[0]] == (r_srt[0]>>32));
	// assert(r_pos[(uint32_t)r_srt[1]] == (r_srt[1]>>32));
	// assert(r_pos[(uint32_t)r_srt[2]] == (r_srt[2]>>32));
	// assert(r_pos[(uint32_t)r_srt[3]] == (r_srt[3]>>32));

	for (i = sidx+1; i < eidx; i++) {

		get_u_offset(ug, &(a[i]), &rs, &re, NULL, NULL);
		a[i].qs = adjust_utg_chain_qoffset(r_srt, r_pos, q_pos, rs);
		a[i].qe = adjust_utg_chain_qoffset(r_srt, r_pos, q_pos, re);
		
		if(i + 1 < eidx) {
			q_pos[0] = a[i].qs; q_pos[1] = a[i].qe; r_pos[0] = rs; r_pos[1] = re;
			r_srt[0] = r_pos[0]; r_srt[0] <<= 32;
			r_srt[1] = r_pos[1]; r_srt[1] <<= 32; r_srt[1] += 1;
			r_srt[2] = r_pos[2]; r_srt[2] <<= 32; r_srt[2] += 2;
			r_srt[3] = r_pos[3]; r_srt[3] <<= 32; r_srt[3] += 3;
			radix_sort_gfa64(r_srt, r_srt + 4);
			// assert(r_pos[(uint32_t)r_srt[0]] == (r_srt[0]>>32));
			// assert(r_pos[(uint32_t)r_srt[1]] == (r_srt[1]>>32));
			// assert(r_pos[(uint32_t)r_srt[2]] == (r_srt[2]>>32));
			// assert(r_pos[(uint32_t)r_srt[3]] == (r_srt[3]>>32));
		}
	}
	


	// if(right_q[0] < 0 || right_q[1] < 0) {
	// 	for (i = sidx+1; i < eidx; i++) {
	// 		get_u_offset(ug, &(a[i]), &rs, &re, NULL, NULL);
	// 		a[i].qs = left_q[0] + get_offset_adjust(rs - left_r[0], left_r[1]-left_r[0], left_q[1]-left_q[0]);
	// 		a[i].qe = left_q[1] + (re - left_r[1]);
	// 		left_q[0] = a[i].qs; left_q[1] = a[i].qe;
	// 		left_r[0] = rs; left_r[1] = re;
	// 	}
	// }

	// if(left_q[0] < 0 || left_q[1] < 0) {
	// 	for (i = eidx-1; i > sidx; i--) {
	// 		get_u_offset(ug, &(a[i]), &rs, &re, NULL, NULL);
	// 		a[i].qe = right_q[1] - get_offset_adjust(right_r[1]-re, right_r[1]-right_r[0], right_q[1]-right_q[0]);
	// 		a[i].qs = right_q[0] - (right_r[0]-rs);
	// 		right_q[0] = a[i].qs; right_q[1] = a[i].qe;
	// 		right_r[0] = rs; right_r[1] = re;
	// 	}
	// }

	// fprintf(stderr, "******[M::%s::] right_q[0]:%ld, right_q[1]:%ld\n", __func__, right_q[0], right_q[1]);
}


void debug_update_uovlp_chain_qse(ma_ug_t *ug, mg_lchain_t *a, int64_t a_n, int64_t ulid)
{
	if(a_n <= 2) return;
	int64_t r0_s, r0_e, q0_s, q0_e, k;
	int64_t r1_s, r1_e, q1_s, q1_e;
	get_u_offset(ug, &(a[0]), &r0_s, &r0_e, &q0_s, &q0_e);
	get_u_offset(ug, &(a[a_n-1]), &r1_s, &r1_e, &q1_s, &q1_e);
	/**if(r0_s <= r1_s && r0_e <= r1_e && q0_s <= q1_s && q0_e <= q1_e)**/ {
		for (k = 1; k + 1 < a_n; k++) a[k].qs = a[k].qe = -1;
		update_uovlp_chain_qse(ug, 0, a_n-1, a, a_n);
		for (k = 1; k < a_n; k++) {
			get_u_offset(ug, a + k - 1, &r0_s, &r0_e, &q0_s, &q0_e);
			get_u_offset(ug, a + k, &r1_s, &r1_e, &q1_s, &q1_e);
			// if(!(r0_s <= r1_s && r0_e <= r1_e && q0_s <= q1_s && q0_e <= q1_e)) {
			// 	fprintf(stderr, "ulid:%ld, r0_s:%ld, r1_s:%ld, r0_e:%ld, r1_e:%ld, q0_s:%ld, q1_s:%ld, q0_e:%ld, q1_e:%ld\n", ulid, 
			// 	r0_s, r1_s, r0_e, r1_e, q0_s, q1_s, q0_e, q1_e);
			// }
			// assert(q0_s <= q1_s && q0_e <= q1_e && q0_s <= q0_e && q1_s <= q1_e);
			if(r0_s <= r1_s) assert(q0_s <= q1_s);
			if(r0_e <= r1_e) assert(q0_e <= q1_e);
			// assert(r0_s <= r1_s && r0_e <= r1_e && q0_s <= q1_s && q0_e <= q1_e);
			assert(q0_s>=0 && q0_e>=0 && q1_s>=0 && q1_e>=0);
		}
	}
}

void fill_unaligned_alignments(ma_ug_t *ug, mg_lchain_t *a, int64_t a_n, int64_t offset, int64_t ulid)
{
	// fprintf(stderr, "a_n->%ld, offset->%ld\n", a_n, offset);
	if(a_n == 0) return;
	int64_t k, l;
	for (k = 0, l = ug->g->seq[a[0].v>>1].len; k < a_n; k++) {
		l -= ug->g->seq[a[k].v>>1].len;
		// fprintf(stderr, "k->%ld, l->%ld [M::utg%.6u%c::%c]\n", k, l, 
		// (a[k].v>>1)+1, "lc"[ug->u.a[a[k].v>>1].circ], "+-"[a[k].v&1]);

		if(a[k].off < 0) a[k].qs = a[k].qe = -1; 
		a[k].off = l; a[k].hash_pre = (uint32_t)-1;
		if(k > 0) a[k].hash_pre = offset + k - 1;

		l += ug->g->seq[a[k].v>>1].len + a[k].dist_pre;
		
	}

	// debug_update_uovlp_chain_qse(ug, a, a_n, ulid);

	for (l = -1, k = 0; k <= a_n; k++) {
		if(k == a_n || a[k].qs >= 0) { ///a[k] and a[l] are anchors
			if(k-l>1) update_uovlp_chain_qse(ug, l, k, a, a_n);
			l = k;
		}
	}
}

void update_ul_vec_t_ug(const ul_idx_t *uref, ul_vec_t *rch, vec_mg_lchain_t *uc, int64_t ulid)
{
	int64_t k, ucn = uc->n, a_n, m, l; ma_ug_t *ug = uref->ug; mg_lchain_t *ix, *a; uc_block_t *z;
	for (k = 0, a_n = 0; k < ucn; k += ix->cnt + 1) {
		ix = &(uc->a[k]); assert(ix->v == (uint32_t)-1); ix->hash_pre = (uint32_t)-1; ix->off = -1;
		fill_unaligned_alignments(ug, uc->a + k + 1, ix->cnt, k + 1, ulid); a_n += ix->cnt;
	}

	///up to now, given a <x> in swap
	///x->ts and x->te are the coordinates in unitig
	///x->qs and x->qe are the coordinates in UL
	///x->dist_pre is the idx of this chain at rch 
	// debug_intermediate_chain(uref->ug, uc->a, uc->n, 1, debug_i);
	// dd_ul_vec_t(uref, swap->a, swap->n, rch);
	rch->bb.n = 0; kv_resize(uc_block_t, rch->bb, (uint64_t)a_n);
	for (k = 0; k < ucn; k += ix->cnt + 1) {
		ix = &(uc->a[k]); assert(ix->v == (uint32_t)-1); ix->hash_pre = (uint32_t)-1; ix->off = -1;
		a = uc->a + k + 1; a_n = ix->cnt;
		for (m = 0; m < a_n; m++) {
			kv_pushp(uc_block_t, rch->bb, &z); 
			z->hid = (a[m].v>>1); z->rev = (a[m].v&1); 
			z->pchain = 1; z->base = 0; z->el = 1;
			z->qs = a[m].qs; z->qe = a[m].qe;
			z->te = a[m].re; z->ts = a[m].rs;
			z->pidx = k + 1 + m; z->pdis = z->aidx = (uint32_t)-1;
		}
	}

	a_n = rch->bb.n; a = uc->a;
	radix_sort_uc_block_t_qe_srt(rch->bb.a, rch->bb.a + rch->bb.n); 
	for (k = 0, l = m = -1; k < a_n; k++) {
		a[rch->bb.a[k].pidx].off = k;
		if(m < (int64_t)rch->bb.a[k].pidx) {
			m = rch->bb.a[k].pidx; l = k;
		}
	}

	for (k = 0; k < a_n; k++) {
		if(a[rch->bb.a[k].pidx].hash_pre == (uint32_t)-1) {
			rch->bb.a[k].pidx = rch->bb.a[k].pdis = rch->bb.a[k].aidx = (uint32_t)-1;
			continue;
		}
		m = rch->bb.a[k].pidx;
		rch->bb.a[k].pidx = a[a[m].hash_pre].off; 
		rch->bb.a[k].pdis = a[a[m].hash_pre].dist_pre;
		rch->bb.a[rch->bb.a[k].pidx].aidx = k;
	}
	// fprintf(stderr, "+ulid->%ld\n", ulid);
	uint32_t sp = (uint32_t)-1, ep = (uint32_t)-1; k = l/**a_n - 1**/;///start from the max chain
	for (l = 0; k >= 0; ) {
		if(sp == (uint32_t)-1 || rch->bb.a[k].qe <= sp) {
            if(sp != (uint32_t)-1) l += ep - sp;
            sp = rch->bb.a[k].qs; ep = rch->bb.a[k].qe;
        } else {
            sp = MIN(sp, rch->bb.a[k].qs);
        }
		if(rch->bb.a[k].pidx == (uint32_t)-1) k = -1;
		else k = rch->bb.a[k].pidx;
	}
	// fprintf(stderr, "-ulid->%ld\n", ulid);
	rch->dd = 0;
	if(sp != (uint32_t)-1) l += ep - sp;
	l = (int64_t)rch->rlen - l;
	if(l == 0) {
		rch->dd = 1;
	} else if(l < ((int64_t)rch->rlen)*0.001) {
		rch->dd = 2;
	}
}

void print_raw_chains(vec_mg_lchain_t *uc, int64_t ulid)
{
	if(uc->n <= 0) return;
	int64_t k, m, ucn = uc->n, k_cnt, a_n; mg_lchain_t *ix, *a; 
	for (k = m = 0; k < ucn; k += k_cnt) {
		ix = &(uc->a[k]); assert(ix->v == (uint32_t)-1); k_cnt = ix->cnt + 1;
		a = uc->a + k + 1; a_n = ix->cnt;
		for (m = 0; m < a_n; m++) {
			if(a[m].off < 0) break;
		}
		if(m < a_n) fprintf(stderr, "ulid->%ld\n", ulid);
	}
}

void hc_shortest_k(void *km0, const asg_t *g, uint32_t src, int32_t n_dst, mg_path_dst_t *dst, int32_t max_dist, int32_t max_k, 
st_mt_t *dst_done, uint64_t *dst_group, vec_sp_node_t *out, vec_mg_pathv_t *res, uint64_t first_src_ban, 
uint64_t detect_mul_way, float len_dif);

void debug_ul_vec_t_chain(void *km, const asg_t *g, ul_vec_t *rch, st_mt_t *dst_done, vec_sp_node_t *out)
{
	if(rch->dd == 0) return;
	uint64_t k, i, v, w, nv; int64_t dd; asg_arc_t *av; 
	mg_path_dst_t dst; uint64_t dst_group; uc_block_t *a = rch->bb.a;
	for (k = 0; k < rch->bb.n; k++) {
		if(a[k].pidx != (uint32_t)-1) assert(a[a[k].pidx].aidx == k);
		if(a[k].aidx != (uint32_t)-1) assert(a[a[k].aidx].pidx == k);
		if(a[k].pidx == (uint32_t)-1) continue;

		v = (a[k].hid<<1)|a[k].rev; v ^= 1;
		w = (a[a[k].pidx].hid<<1)|a[a[k].pidx].rev; w ^= 1;
		
		nv = asg_arc_n(g, v); av = asg_arc_a(g, v);
		for (i = 0; i < nv; i++) {
			if(av[i].v == w) break;
		}
		// if(i >= nv) {
		// 	// fprintf(stderr, "[M::%s::]\n", __func__);
		// 	fprintf(stderr, "[M::%s::]\tutg%.6dl(%c)\t->\tutg%.6dl(%c)\n", __func__, (int32_t)(v>>1)+1, "+-"[v&1], (int32_t)(w>>1)+1, "+-"[w&1]);
		// }
		if(i < nv) {
			dd = (int64_t)((uint32_t)(av[i].ul));
		} else {
			memset(&dst, 0, sizeof(dst));
            dst.v = w; 
            dst.target_dist = a[k].pdis;
            dst.target_hash = 0; dst.check_hash = 0;
			hc_shortest_k(km, g, v, 1, &dst, dst.target_dist, MG_MAX_SHORT_K, dst_done, &dst_group, out, NULL, 1, 0, 0);
			dd = dst.dist;
		}
		if(a[k].pdis != dd) {
			fprintf(stderr, "[M::%s::]\tutg%.6dl(%c)\t->\tutg%.6dl(%c)\tdist_pre:%d\td:%ld\n", __func__, (int32_t)(v>>1)+1, "+-"[v&1], 
			(int32_t)(w>>1)+1, "+-"[w&1], a[k].pdis, dd);
		}
	}
}

int64_t gl_chain_refine_advance_combine(mg_tbuf_t *b, ul_vec_t *rch, overlap_region_alloc* olist, Correct_dumy* dumy, haplotype_evdience_alloc *hap, st_mt_t *sps, glchain_t *ll, gdpchain_t *gdp, const ul_idx_t *uref, double diff_ec_ul, int64_t winLen, int64_t qlen, const ug_opt_t *uopt, 
int64_t debug_i, void *km)
{
	ll->tk.n = ll->lo.n = 0;
	kv_ul_ov_t *idx = &(ll->lo);
	uint64_t o2 = gl_chain_gen(olist, uref, idx, 0, hap, km);
	if(idx->n == 0) return 0;
	// fprintf(stderr, "[M::%s] qlen:%ld, idx->n:%u\n", __func__, qlen, (uint32_t)idx->n);
	int64_t max_idx, occ = 0, f = 0;

	kv_resize_km(km, uint64_t, ll->srt.a, idx->n);
    kv_resize_km(km, uint64_t, hap->snp_srt, idx->n);
	kv_resize_km(km, ul_ov_t, ll->tk, idx->n);
	///chain exact U-matches
	occ = gl_chain_advance(idx, ll->tk.a, uref, uopt, G_CHAIN_BW, /**diff_ec_ul**/N_GCHAIN_RATE, qlen, UG_SKIP, dumy->overlapID, ll->srt.a.a, hap->snp_srt.a, G_CHAIN_TRANS_WEIGHT, 0, NULL, uref->ug, debug_i, km);
	if(occ) {
		if(ff_chain(idx, qlen, 0.99/**P_CHAIN_COV**/, G_CHAIN_TRANS_RATE, ll->tk.a, NULL, NULL, NULL, diff_ec_ul, winLen, km)) {
			f = l2g_res_chain(uref->ug, ll->tk.a+idx->a[idx->n-1].ts, idx->a[idx->n-1].te-idx->a[idx->n-1].ts, &(gdp->swap), -1/**N_GCHAIN_RATE**/);
		} else if(o2) {///means there are trans overlaps
			gl_chain_gen(olist, uref, idx, 1, hap, km);
			kv_resize_km(km, uint64_t, ll->srt.a, idx->n);
    		kv_resize_km(km, uint64_t, hap->snp_srt, idx->n);
			kv_resize_km(km, ul_ov_t, ll->tk, idx->n);
			///chain all U-matches
			occ = gl_chain_advance(idx, ll->tk.a, uref, uopt, G_CHAIN_BW, /**diff_ec_ul**/N_GCHAIN_RATE, qlen, UG_SKIP, dumy->overlapID, ll->srt.a.a, hap->snp_srt.a, G_CHAIN_TRANS_WEIGHT, 0, NULL, uref->ug, debug_i, km);
			if(ff_chain(idx, qlen, 0.99/**P_CHAIN_COV**/, G_CHAIN_TRANS_RATE, ll->tk.a, olist, hap, uref, diff_ec_ul, winLen, km)) {
				f = l2g_res_chain(uref->ug, ll->tk.a+idx->a[idx->n-1].ts, idx->a[idx->n-1].te-idx->a[idx->n-1].ts, &(gdp->swap), -1/**N_GCHAIN_RATE**/);
			}
		}
	}

	if(!f) {
		gl_chain_gen(olist, uref, idx, 1, hap, km);
		l2g_chain(uref, idx, &(gdp->l)); ll->tk.n = 0;
		///buffer
		kv_resize(uint64_t, ll->srt.a, gdp->l.n); kv_resize(uint64_t, hap->snp_srt, gdp->l.n);
		kv_resize(uint64_t, gdp->v, gdp->l.n); kv_resize(int64_t, gdp->f, gdp->l.n);
		max_idx = hc_gchain1_dp(b->km, uref, uref->ug, &(gdp->l), &(gdp->swap), &(gdp->dst), &(gdp->out), &(gdp->path), rch->rlen,
		uopt, G_CHAIN_BW, diff_ec_ul, -1, ll->srt.a.a, sps, gdp->f.a, hap->snp_srt.a, gdp->v.a);
		if(max_idx >= 0 && gen_max_gchain_adv(b->km, uref, debug_i, sps, &(gdp->l), &(ll->tk), NULL, rch->rlen, P_CHAIN_COV, 0.3/**P_FRAGEMENT_PRIMARY_CHAIN_COV**/, 
			0.1/**P_FRAGEMENT_PRIMARY_SECOND_COV**/, PRIMARY_UL_CHAIN_MIN, uref->ug->g, &(gdp->dst_done), &(gdp->out), &(gdp->path), ll->srt.a.a, &(gdp->swap))) {
			// print_raw_chains(&(gdp->swap), debug_i);
			f = check_trans_rate_gap(&(gdp->swap), &(ll->tk), olist, hap, uref, diff_ec_ul, winLen, G_CHAIN_TRANS_RATE);
		}
	}

	if(f) update_ul_vec_t_ug(uref, rch, &(gdp->swap), debug_i);
	// debug_ul_vec_t_chain(km, uref->ug->g, rch, &(gdp->dst_done), &(gdp->out));
	return 1;
}

uint64_t kv_ul_ov_t_statistics(kv_ul_ov_t *olist, uint64_t qn, int64_t *occ)
{
	int64_t k, l = 0;
	uint32_t sp = (uint32_t)-1, ep = (uint32_t)-1;
	for (k = olist->n-1; k >= 0 && olist->a[k].qn == qn; k--) {
		if(!(olist->a[k].el)) continue;
		if(!(olist->a[k].tn&((uint32_t)(0x80000000)))) continue;
		if(sp == (uint32_t)-1 || olist->a[k].qe <= sp) {
			if(sp != (uint32_t)-1) l += ep - sp;
			sp = olist->a[k].qs;
			ep = olist->a[k].qe;
		} else {
			sp = MIN(sp, olist->a[k].qs);
		}
		(*occ)++;
	}
	if(sp != (uint32_t)-1) l += ep - sp;
	return l;
}

static void worker_for_ul_scall_alignment(void *data, long i, int tid) // callback for kt_for()
{
    utepdat_t *s = (utepdat_t*)data;
	ha_ovec_buf_t *b = s->hab[tid];
	glchain_t *bl = &(s->ll[tid]);
	int64_t /**rid = s->id+i,**/ winLen = MIN((((double)THRESHOLD_MAX_SIZE)/s->opt->diff_ec_ul), WINDOW);
	uint64_t align = 0;
	int fully_cov, abnormal;
	void *km = s->buf?(s->buf[tid]?s->buf[tid]->km:NULL):NULL;
	// if(s->id+i!=601) return;
	// fprintf(stderr, "[M::%s] rid:%ld\n", __func__, s->id+i);
	// if (memcmp(UL_INF.nid.a[s->id+i].a, "d0aab024-b3a7-40fb-83cc-22c3d6d951f8", UL_INF.nid.a[s->id+i].n-1)) return;
	// fprintf(stderr, "[M::%s::] ==> len: %lu\n", __func__, s->len[i]);
	ha_get_ul_candidates_interface(b->abl, i, s->seq[i], s->len[i], s->opt->w, s->opt->k, s->uu, &b->olist, &b->olist_hp, &b->clist, s->opt->bw_thres, 
		s->opt->max_n_chain, 1, &(b->k_flag), &b->r_buf, &(b->tmp_region), NULL, &(b->sp), asm_opt.hom_cov, km);
	
	clear_Cigar_record(&b->cigar1);
	clear_Round2_alignment(&b->round2);
	// return;
	// b->num_correct_base += overlap_statistics(&b->olist, NULL, 0);	

	b->self_read.seq = s->seq[i]; b->self_read.length = s->len[i]; b->self_read.size = 0;
	correct_ul_overlap(&b->olist, s->uu, &b->self_read, &b->correct, &b->ovlp_read, &b->POA_Graph, &b->DAGCon,
			&b->cigar1, &b->hap, &b->round2, 0, 1, &fully_cov, &abnormal, s->opt->diff_ec_ul, winLen, km);

	// uint64_t k;
	// for (k = 0; k < b->olist.length; k++) {
	// 	if(b->olist.list[k].is_match == 1) b->num_correct_base += b->olist.list[k].x_pos_e+1-b->olist.list[k].x_pos_s;
	// 	if(b->olist.list[k].is_match == 2) b->num_recorrect_base += b->olist.list[k].x_pos_e+1-b->olist.list[k].x_pos_s;
	// }


	// gl_chain_refine(&b->olist, &b->correct, &b->hap, bl, s->uu, s->opt->diff_ec_ul, winLen, s->len[i], km);
	gl_chain_refine_advance(&b->olist, &b->correct, &b->hap, bl, &(s->sps[tid]), s->uu, s->opt->diff_ec_ul, winLen, s->len[i], s->uopt, s->id+i, km);
	// return;
	// b->num_read_base += b->self_read.length;
	// b->num_correct_base += b->correct.corrected_base;
	// b->num_recorrect_base += b->round2.dumy.corrected_base;
	memset(&b->self_read, 0, sizeof(b->self_read));
	align = kv_ul_ov_t_statistics(&(bl->tk), i, &(b->num_recorrect_base));
	if(align == s->len[i]) {
		free(s->seq[i]); s->seq[i] = NULL;
	}
	b->num_correct_base += align;

	// uint64_t k;
	// b->num_read_base += overlap_statistics(&b->olist, NULL, NULL, 1);
	// for (k = 0; k < bl->tk.n; k++) {
	// 	if(bl->tk.a[k].sec == 0) b->num_correct_base += bl->tk.a[k].qe - bl->tk.a[k].qs;
	// 	if(bl->tk.a[k].sec > 0) b->num_recorrect_base += bl->tk.a[k].qe - bl->tk.a[k].qs;
	// }
	// for (k = 0; k < bl->lo.n; k++) {
	// 	b->num_read_base += bl->lo.a[k].qe - bl->lo.a[k].qs;
	// }

	// uint32_t l1 = overlap_statistics(&b->olist, s->uu->ug, 1), l2 = overlap_statistics(&b->olist, s->uu->ug, 2);
    // 
	// if(l1 == 0 && l2 > 0) fprintf(stderr, "[M::%s::%lu::no_match]\n", UL_INF.nid.a[s->id+i].a, s->len[i]);
	// fprintf(stderr, "[M::%s::%lu::] l1->%u; l2->%u\n", UL_INF.nid.a[s->id+i].a, s->len[i], l1, l2);
	if(km) {
        destory_overlap_region_alloc_buf(km, &b->olist, 1); 
        destory_Correct_dumy_buf(km, &b->correct, 1);
        destoryHaplotypeEvdience_buf(km, &b->hap, 1);
    }
}

static void worker_for_ul_rescall_alignment(void *data, long i, int tid) // callback for kt_for()
{
    utepdat_t *s = (utepdat_t*)data;
	ha_ovec_buf_t *b = s->hab[tid];
	glchain_t *bl = &(s->ll[tid]);
	int64_t /**rid = s->id+i,**/ winLen = MIN((((double)THRESHOLD_MAX_SIZE)/s->opt->diff_ec_ul), WINDOW);
	// uint64_t align = 0;
	int fully_cov, abnormal;
	// void *km = s->buf?(s->buf[tid]?s->buf[tid]->km:NULL):NULL;
	// if(s->id+i!=873) return;
	// fprintf(stderr, "\n[M::%s] rid:%ld, len:%lu\n", __func__, s->id+i, s->len[i]);
	// if (memcmp(UL_INF.nid.a[s->id+i].a, "d0aab024-b3a7-40fb-83cc-22c3d6d951f8", UL_INF.nid.a[s->id+i].n-1)) return;
	// fprintf(stderr, "[M::%s::] ==> len: %lu\n", __func__, s->len[i]);
	ha_get_ul_candidates_interface(b->abl, i, s->seq[i], s->len[i], s->opt->w, s->opt->k, s->uu, &b->olist, &b->olist_hp, &b->clist, s->opt->bw_thres, 
		s->opt->max_n_chain, 1, &(b->k_flag), &b->r_buf, &(b->tmp_region), NULL, &(b->sp), 1, NULL);
	
	clear_Cigar_record(&b->cigar1);
	clear_Round2_alignment(&b->round2);
	// return;
	// b->num_correct_base += overlap_statistics(&b->olist, NULL, 0);	

	b->self_read.seq = s->seq[i]; b->self_read.length = s->len[i]; b->self_read.size = 0;
	correct_ul_overlap(&b->olist, s->uu, &b->self_read, &b->correct, &b->ovlp_read, &b->POA_Graph, &b->DAGCon,
			&b->cigar1, &b->hap, &b->round2, 0, 1, &fully_cov, &abnormal, s->opt->diff_ec_ul, winLen, NULL);

	// uint64_t k;
	// for (k = 0; k < b->olist.length; k++) {
	// 	if(b->olist.list[k].is_match == 1) b->num_correct_base += b->olist.list[k].x_pos_e+1-b->olist.list[k].x_pos_s;
	// 	if(b->olist.list[k].is_match == 2) b->num_recorrect_base += b->olist.list[k].x_pos_e+1-b->olist.list[k].x_pos_s;
	// }


	// gl_chain_refine(&b->olist, &b->correct, &b->hap, bl, s->uu, s->opt->diff_ec_ul, winLen, s->len[i], km);
	gl_chain_refine_advance_combine(s->buf[tid], &(UL_INF.a[i]), &b->olist, &b->correct, &b->hap, &(s->sps[tid]), bl, &(s->gdp[tid]), s->uu, s->opt->diff_ec_ul, winLen, s->len[i], s->uopt, s->id+i, NULL);
	// return;
	// b->num_read_base += b->self_read.length;
	// b->num_correct_base += b->correct.corrected_base;
	// b->num_recorrect_base += b->round2.dumy.corrected_base;
	memset(&b->self_read, 0, sizeof(b->self_read));
	if(UL_INF.a[i].dd) {
		free(s->seq[i]); s->seq[i] = NULL; b->num_correct_base++;
	}
	s->hab[tid]->num_read_base++;

	// align = kv_ul_ov_t_statistics(&(bl->tk), i, &(b->num_recorrect_base));
	// if(align == s->len[i]) {
	// 	free(s->seq[i]); s->seq[i] = NULL;
	// }
	// b->num_correct_base += align;

	// uint64_t k;
	// b->num_read_base += overlap_statistics(&b->olist, NULL, NULL, 1);
	// for (k = 0; k < bl->tk.n; k++) {
	// 	if(bl->tk.a[k].sec == 0) b->num_correct_base += bl->tk.a[k].qe - bl->tk.a[k].qs;
	// 	if(bl->tk.a[k].sec > 0) b->num_recorrect_base += bl->tk.a[k].qe - bl->tk.a[k].qs;
	// }
	// for (k = 0; k < bl->lo.n; k++) {
	// 	b->num_read_base += bl->lo.a[k].qe - bl->lo.a[k].qs;
	// }

	// uint32_t l1 = overlap_statistics(&b->olist, s->uu->ug, 1), l2 = overlap_statistics(&b->olist, s->uu->ug, 2);
    // 
	// if(l1 == 0 && l2 > 0) fprintf(stderr, "[M::%s::%lu::no_match]\n", UL_INF.nid.a[s->id+i].a, s->len[i]);
	// fprintf(stderr, "[M::%s::%lu::] l1->%u; l2->%u\n", UL_INF.nid.a[s->id+i].a, s->len[i], l1, l2);
	// fprintf(stderr, "[M::%s::rid->%ld] done\n", __func__, s->id+i);
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

void push_uc_block_t(const ug_opt_t *uopt, kv_ul_ov_t *z, char **seq, uint64_t *len, uint64_t b_id)
{
	uint64_t k, l, rid;
	for (k = 1, l = 0; k <= z->n; k++) {
        if(k == z->n || z->a[k].qn != z->a[l].qn) {
			rid = b_id + z->a[l].qn;
			// fprintf(stderr, "rid->%lu, b_id->%lu, l->%lu, z->a[l].qn->%u, len[z->a[l].qn]->%lu, seq[z->a[l].qn]->%u\n", rid, b_id, l, z->a[l].qn, len[z->a[l].qn], seq[z->a[l].qn]?1:0);
			append_ul_t(&UL_INF, &rid, NULL, 0, seq[z->a[l].qn], len[z->a[l].qn], z->a + l, k - l, P_CHAIN_COV, uopt, 0);
			// append_ul_t_back(&UL_INF, &rid, NULL, 0, seq[z->a[l].qn], len[z->a[l].qn], z->a + l, k - l, P_CHAIN_COV);
			l = k;
		}
	}
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
		s->opt = p->opt; s->uu = p->uu; s->uopt = p->uopt; s->rg = p->rg;
        while ((ret = kseq_read(p->ks)) >= 0) 
        {
            if (p->ks->seq.l < (uint64_t)p->opt->k) continue;
            if (s->n == s->m) {
                s->m = s->m < 16? 16 : s->m + (s->n>>1);
                REALLOC(s->len, s->m);
                REALLOC(s->seq, s->m);
            }
			
			append_ul_t(&UL_INF, NULL, p->ks->name.s, p->ks->name.l, NULL, 0, NULL, 0, P_CHAIN_COV, s->uopt, 0);		
            l = p->ks->seq.l;
            MALLOC(s->seq[s->n], l);
            s->sum_len += l;
            memcpy(s->seq[s->n], p->ks->seq.s, l);
			// fprintf(stderr, "s->n->%d, l->%lu\n", s->n, l);
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
		CALLOC(s->hab, p->n_thread); CALLOC(s->ll, p->n_thread); CALLOC(s->sps, p->n_thread); 
		// CALLOC(s->buf, p->n_thread);
		for (i = 0; i < p->n_thread; ++i) {
			// s->buf[i] = mg_tbuf_init();
			// s->hab[i] = ha_ovec_buf_init(s->buf[i]->km, 0, 0, 1);
			// s->buf[i] = NULL;
			// s->hab[i] = ha_ovec_buf_init(NULL, 0, 0, 1);
			s->hab[i] = ha_ovec_init(0, 0, 1);
		}
		fprintf(stderr, "[M::%s::Start] ==> s->id: %lu, s->n:% d\n", __func__, s->id, s->n);
		kt_for(p->n_thread, worker_for_ul_scall_alignment, s, s->n);
		fprintf(stderr, "[M::%s::Done] ==> s->id: %lu, s->n:% d\n", __func__, s->id, s->n);
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
		for (i = 0; i < p->n_thread; ++i) {
			s->num_bases += s->hab[i]->num_read_base;
			s->num_corrected_bases += s->hab[i]->num_correct_base;
			s->num_recorrected_bases += s->hab[i]->num_recorrect_base;
			// mg_tbuf_destroy(s->buf[i]);
			ha_ovec_destroy(s->hab[i]); kv_destroy(s->sps[i]);
			free(s->ll[i].lo.a); /**free(s->ll[i].tk.a);**/ free(s->ll[i].srt.a.a);
		}
		free(s->hab); free(s->sps); /**free(s->ll);**/ // free(s->buf);
		//free(s->mzs); free(s->sps);
		return s;
    }
    else if (step == 2) { // step 3: dump
        utepdat_t *s = (utepdat_t*)in;
		uint64_t i, rid;
		p->num_bases += s->num_bases;
		p->num_corrected_bases += s->num_corrected_bases;
		p->num_recorrected_bases += s->num_recorrected_bases;
		fprintf(stderr, "[M::%s::dump_start] ==> s->id: %lu, s->n:% d\n", __func__, s->id, s->n);
		for (i = 0; i < p->n_thread; ++i) {
			push_uc_block_t(s->uopt, &(s->ll[i].tk), s->seq, s->len, s->id);
			free(s->ll[i].tk.a);
		}
		for (i = 0; i < (uint64_t)s->n; ++i) {
			rid = s->id + i;
			if(UL_INF.n > rid && UL_INF.a[rid].rlen != s->len[i]) {
				append_ul_t(&UL_INF, &rid, NULL, 0, s->seq[i], s->len[i], NULL, 0, P_CHAIN_COV, s->uopt, 0);
			}
			free(s->seq[i]);
		}
		fprintf(stderr, "[M::%s::dump_done] ==> s->id: %lu, s->n:% d\n", __func__, s->id, s->n);
		/**
        for (i = 0; i < (uint64_t)s->n; ++i) {
			///debug
			
            // if(s->pos[i].s == (uint64_t)-1) continue;
            // kv_push(pe_hit, p->hits.a, s->pos[i]);
			// if(!s->gcs[i]) continue;
			// dump_gaf(&(p->hits), s->gcs[i], 1);
			// free(s->gcs[i]->gc); free(s->gcs[i]->a); free(s->gcs[i]->lc); free(s->gcs[i]);
			
			rid = s->id + i;
			append_ul_t(&UL_INF, &rid, NULL, 0, s->seq[i], s->len[i], NULL, 0, P_CHAIN_COV);	
			// fprintf(stderr, "%.*s\n", (int)s->len[i], s->seq[i]);	
			free(s->seq[i]); p->total_base += s->len[i];
        }
		**/
		///debug
		/**
        free(s->gcs);
		**/
		
        free(s->ll); free(s->len); free(s->seq); free(s);
    }
    return 0;
}


static void *worker_ul_rescall_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
    uldat_t *p = (uldat_t*)data;
    ///uint64_t total_base = 0, total_pair = 0;
    if (step == 0) { // step 1: read a block of sequences
        int ret;
        uint64_t l;
		utepdat_t *s;
		CALLOC(s, 1);
        s->ha_flt_tab = p->ha_flt_tab; s->ha_idx = p->ha_idx; s->id = p->total_pair; 
		s->opt = p->opt; s->uu = p->uu; s->uopt = p->uopt; s->rg = p->rg;
        while ((ret = kseq_read(p->ks)) >= 0) 
        {
            if (p->ks->seq.l < (uint64_t)p->opt->k) continue;
            if (s->n == s->m) {
                s->m = s->m < 16? 16 : s->m + (s->n>>1);
                REALLOC(s->len, s->m);
                REALLOC(s->seq, s->m);
            }
			
			// append_ul_t(&UL_INF, NULL, p->ks->name.s, p->ks->name.l, NULL, 0, NULL, 0, P_CHAIN_COV, s->uopt);		
            l = p->ks->seq.l;
            MALLOC(s->seq[s->n], l);
            s->sum_len += l;
            memcpy(s->seq[s->n], p->ks->seq.s, l);
			// fprintf(stderr, "s->n->%d, l->%lu\n", s->n, l);
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
		CALLOC(s->hab, p->n_thread); CALLOC(s->ll, p->n_thread); CALLOC(s->buf, p->n_thread); 
		CALLOC(s->gdp, p->n_thread); CALLOC(s->mzs, p->n_thread); CALLOC(s->sps, p->n_thread); 

		// CALLOC(s->buf, p->n_thread);
		for (i = 0; i < p->n_thread; ++i) {
			s->hab[i] = ha_ovec_init(0, 0, 1); s->buf[i] = mg_tbuf_init();
		}
		kt_for(p->n_thread, worker_for_ul_rescall_alignment, s, s->n);
		for (i = 0; i < p->n_thread; ++i) {
			p->num_bases += s->hab[i]->num_read_base;
			p->num_corrected_bases += s->hab[i]->num_correct_base;
			// s->num_recorrected_bases += s->hab[i]->num_recorrect_base;
			ha_ovec_destroy(s->hab[i]); hc_glchain_destroy(&(s->ll[i]));
			mg_tbuf_destroy(s->buf[i]); hc_gdpchain_destroy(&(s->gdp[i])); 
			kv_destroy(s->mzs[i]); kv_destroy(s->sps[i]); //free(s->seq[i]);
		}
		free(s->hab); free(s->ll); ///free(s->len); free(s->seq); 
		free(s->buf); free(s->gdp); free(s->mzs); free(s->sps); ///free(s);
		return s;
    } else if (step == 2) { // step 3: dump
		utepdat_t *s = (utepdat_t*)in; int64_t i, rid;
		for (i = 0; i < s->n; ++i) {
			rid = s->id + i;
			if(UL_INF.a[rid].dd == 0 && p->ucr_s && p->ucr_s->flag == 1) {
				write_compress_base_disk(p->ucr_s->fp, rid, s->seq[i], s->len[i], &(p->ucr_s->u));
			}
			free(s->seq[i]);
		}
		free(s->len); free(s->seq); free(s);
	}
    return 0;
}

int32_t init_ucr_file_t(uldat_t *sl, char* file, uint64_t mode)
{
	if(mode == 1 || mode == 2) {
		char *gfa_name = (char*)malloc(strlen(file)+25);
		sprintf(gfa_name, "%s.uidx.ucr.bin", file);
		CALLOC(sl->ucr_s, 1); 
		sl->ucr_s->flag = mode;
		sl->ucr_s->fp = fopen(gfa_name, mode==1?"w":"r");
		if (!(sl->ucr_s->fp)) {
			free(gfa_name);
			return 0;
		}
		free(gfa_name);
		return 1;
	}	
	return 0;
}


void destory_ucr_file_t(uldat_t *sl)
{
	if(sl->ucr_s) {
		free(sl->ucr_s->u.r_base.a);
		free(sl->ucr_s->u.bb.a);
		free(sl->ucr_s->u.N_site.a);

		fclose(sl->ucr_s->fp);
		free(sl->ucr_s);
		sl->ucr_s = NULL;
	}
}


void debug_sl_compress_base_disk_0(uldat_t *sl, char* gfa_name)
{
	int32_t ret, rid = 0, sr_0, sr_1; gzFile fp;
	init_ucr_file_t(sl, gfa_name, 1);
	fp = gzopen(gfa_name, "r"); assert(fp);
    sl->ks = kseq_init(fp);
	while ((ret = kseq_read(sl->ks)) >= 0) {
		write_compress_base_disk(sl->ucr_s->fp, rid, sl->ks->seq.s, sl->ks->seq.l, &(sl->ucr_s->u));
		rid++;   
	}
	kseq_destroy(sl->ks);
    gzclose(fp);

	destory_ucr_file_t(sl);


	uint64_t ulid; uint32_t ulen; kvec_t(char) des; kv_init(des);
	init_ucr_file_t(sl, gfa_name, 2); rid = 0;

	fp = gzopen(gfa_name, "r"); assert(fp);
    sl->ks = kseq_init(fp);
	while (1) {
		sr_0 = kseq_read(sl->ks); des.n = 0; 
		if(sr_0 >= 0) sr_0 = 1;
		else sr_0 = 0;
		if(sr_0 == 0) sl->ks->seq.l = 1;
		kv_resize(char, des, sl->ks->seq.l);
		sr_1 = load_compress_base_disk(sl->ucr_s->fp, &ulid, des.a, &ulen, &(sl->ucr_s->u));
		if(sr_0 != sr_1) fprintf(stderr, "[M::%s::] rid->%d, sr_0->%d, sr_1->%d\n", __func__, rid, sr_0, sr_1);
		assert(sr_0 == sr_1);
		if(sr_0 == 0 || sr_1 == 0) break;
		// if(rid != (int64_t)ulid) fprintf(stderr, "[M::%s::] rid->%d, ulid->%lu, sr_0->%d, sr_1->%d\n", __func__, rid, ulid, sr_0, sr_1);
		assert(rid == (int64_t)ulid);
		assert(sl->ks->seq.l == ulen);
		assert(memcmp(sl->ks->seq.s, des.a, ulen) == 0);
		rid++;   
	}
	kseq_destroy(sl->ks); kv_destroy(des);
    gzclose(fp);

	destory_ucr_file_t(sl);
	fprintf(stderr, "[M::%s::] ==> Have checked %d UL reads\n", __func__, rid);
	
}

void debug_sl_compress_base_disk_0(uldat_t *sl, const enzyme *fn)
{
	int32_t i;
	for (i = 0; i < fn->n; i++) debug_sl_compress_base_disk_0(sl, fn->a[i]);
	exit(1);
}



utg_rid_dt *get_r_ug_region(utg_rid_t *idx, uint64_t *n, uint64_t rid)
{
	(*n) = idx->idx[rid+1] - idx->idx[rid];
	return (*n)?idx->p.a + idx->idx[rid]:NULL;
}

void rov2uov(uint64_t rid, const ul_idx_t *uref, utg_rid_dt *ru_map, uc_block_t *rovlp, ul_ov_t *res, uint32_t adjust_rev)
{
	uint64_t ori = ru_map->u&1, ts, te;
	if(!ori) {
		ts = rovlp->ts; te = rovlp->te;
	} else {
		ts = uref->r_ug->rg->seq[rid].len - rovlp->te;
		te = uref->r_ug->rg->seq[rid].len - rovlp->ts;
	}
	ts += ru_map->off; te += ru_map->off;
	memset(res, 0, sizeof(*res));
	res->qn = 0; res->qs = rovlp->qs; res->qe = rovlp->qe; 
	res->tn = ru_map->u>>1; res->ts = ts; res->te = te; 
	res->el = rovlp->el; res->rev = (rovlp->rev == ori?0:1);
	if(adjust_rev && res->rev) {///for linear chaining
		res->ts = uref->ug->g->seq[res->tn].len - te;
		res->te = uref->ug->g->seq[res->tn].len - ts;
	}
}

void print_ul_ov_t(ul_ov_t *xs, const char* cmd)
{
	fprintf(stderr, "%s\t%s\t%u\t%u\t%c\t%.*s\t%u\t%u\n", cmd, UL_INF.nid.a[xs->qn].a, xs->qs, xs->qe, 
	"+-"[xs->rev], (int)Get_NAME_LENGTH(R_INF, ((xs->tn<<1)>>1)), Get_NAME(R_INF, ((xs->tn<<1)>>1)), xs->ts, xs->te);
}

void gl_rg2ug_gen(ul_vec_t *r_cl, kv_ul_ov_t *u_cl, const ul_idx_t *uref, uint64_t is_el, uint64_t n_pchain)
{
	uint64_t k, a_k, a_n; uc_block_t *z; utg_rid_dt *a; ul_ov_t *p;
	u_cl->n = 0;
	for (k = 0; k < r_cl->bb.n; k++) {
		z = &(r_cl->bb.a[k]);
		if(z->base) continue;
		if(is_el && (!(z->el))) continue;
		if(z->pchain == n_pchain) continue;
		a = get_r_ug_region(uref->r_ug, &a_n, z->hid);
		if(!a) continue;
		for (a_k = 0; a_k < a_n; a_k++) {
			kv_pushp(ul_ov_t, *u_cl, &p);
			// fprintf(stderr, "\n+[M::%s::] %u\t%u\t%c\t%.*s(%u)\t%u\t%u\n", __func__, z->qs, z->qe, "+-"[z->rev], 
			// (int)Get_NAME_LENGTH(R_INF, z->hid), Get_NAME(R_INF, z->hid), (uint32_t)Get_READ_LENGTH(R_INF, z->hid), z->ts, z->te);
			// fprintf(stderr, "*[M::%s::] utg%.6d%c(%u)\t%c\t%u\n", __func__, 
			// (int32_t)(a[a_k].u>>1)+1, "lc"[uref->ug->u.a[a[a_k].u>>1].circ], uref->ug->u.a[a[a_k].u>>1].len, 
			// "+-"[a[a_k].u&1], a[a_k].off);
			rov2uov(z->hid, uref, &(a[a_k]), z, p, 1);
			p->el = 1; p->tn <<= 1; p->tn |= p->rev; p->qn = k/**uref->r_ug->idx[z->hid] + a_k**/;//for linear chain
			// if(k == 2) {
			// 	fprintf(stderr, "[M::%s::] p->ts:%u, p->te:%u, z->ts:%u, z->te:%u, a[a_k].off:%u\n", __func__, p->ts, p->te, z->ts, z->te, a[a_k].off);
			// }
			// fprintf(stderr, "-[M::%s::] %u\t%u\t%c\tutg%.6d%c(%u)\t%u\t%u\n", __func__, p->qs, p->qe, "+-"[p->rev], 
			// 				(int32_t)(p->tn>>1)+1, "lc"[uref->ug->u.a[p->tn>>1].circ], uref->ug->u.a[p->tn>>1].len, p->ts, p->te);
		}
	}
}

void adjust_rev_tse(ul_ov_t *x, int64_t tlen, int64_t *ts, int64_t *te)
{
	*ts = x->ts; *te = x->te;
	if(x->rev) {
		*ts = tlen - x->te; *te = tlen - x->ts;
	}
}

uint64_t get_add_cov_score(const ul_idx_t *uref, int64_t ps, int64_t pe, int64_t cs, int64_t ce, int64_t uid, int64_t *cov_i)
{
	int64_t os = MAX(ps, cs), oe = MIN(pe, ce);
	int64_t ovlp = ((oe > os)? (oe - os):0);
	// fprintf(stderr, "ovlp:%ld, os:%ld, oe:%ld, ps:%ld, pe:%ld, cs:%ld, ce::%ld\n", ovlp, os, oe, ps, pe, cs, ce);
	if(ovlp > 0) {
		return (os>cs?retrieve_u_cov_region(uref, uid, 0, cs, os, cov_i):0) + 
						(ce>oe?retrieve_u_cov_region(uref, uid, 0, oe, ce, cov_i):0);
	} 

	return retrieve_u_cov_region(uref, uid, 0, cs, ce, cov_i);
}

uint64_t linear_chain_dp(ul_ov_t *ch, int64_t ch_n, ul_ov_t *sv, const ul_idx_t *uref, const ug_opt_t *uopt, int64_t bw, 
double diff_ec_ul, int64_t qlen, int64_t max_skip, uint64_t *idx, uint64_t *track, ma_ug_t *ug, int64_t chain_offset)
{	///all in[].el must be 1
	if(ch_n == 0) return 0;
    int64_t /**mm_ovlp, x,**/ i, j, k, sc, csc, mm_sc, mm_idx, its, ite, jts, jte, cov_i, dq, dt, dd, mm;
    ul_ov_t *li = NULL, *lj = NULL;
    radix_sort_ul_ov_srt_qe(ch, ch + ch_n);
	for (i = 1, j = 0; i <= ch_n; i++) {
		// if(i < ch_n) {
		// 	li = &(ch[i]); 
		// 	fprintf(stderr, "##(%ld) %u\t%u\t%c\tutg%.6d%c(%u)\t%u\t%u\tmm_idx:%ld\tmm_sc:%ld\n", i, li->qs, li->qe, "+-"[li->rev], 
        // 	(int32_t)(li->tn)+1, "lc"[uref->ug->u.a[li->tn].circ], uref->ug->u.a[li->tn].len, li->ts, li->te, mm_idx, mm_sc);
		// }
        if (i == ch_n || ch[i].qe != ch[j].qe) {
            if(i - j > 1) {
                radix_sort_ul_ov_srt_qs(ch+j, ch+i);
            }
            j = i;
        }
    }

	// fprintf(stderr, "[M::%s::] ch_n:%ld\n", __func__, ch_n);
	for (i = 0; i < ch_n; ++i) {
		li = &(ch[i]); 
		// mm_ovlp = max_ovlp_src(uopt, ((li->tn<<1)|li->rev)^1);  
		// x = (li->qs + mm_ovlp)*diff_ec_ul; 
        // if(x < bw) x = bw; 
        // x += li->qs + mm_ovlp;
        // if (x > qlen+1) x = qlen+1;
        // x = find_ul_ov_max(i, ch, x+G_CHAIN_INDEL); 
		adjust_rev_tse(li, ug->g->seq[li->tn].len, &its, &ite);
		cov_i = 0;
		csc = retrieve_u_cov_region(uref, li->tn, 0, its, ite, &cov_i);
		mm_sc = csc; mm_idx = -1; 
		for (j = i-1/**x**/; j >= 0; --j) {
			lj = &(ch[j]); 
			// fprintf(stderr, "<0>\n");
			if(lj->qs <= li->qs && lj->qe <= li->qe && lj->ts <= li->ts && lj->te <= li->te) {///co-linear
				assert(li->tn == lj->tn && li->rev == lj->rev);
				// fprintf(stderr, "<1>\n");
				if(lj->qs == li->qs && lj->qe == li->qe && lj->ts == li->ts && lj->te == li->te) continue;
				dq = li->qe - lj->qs; dt = li->te - lj->ts;
				dd = (dq>dt? dq-dt:dt-dq);
				mm = MAX(dq, dt); mm *= diff_ec_ul; if(mm < bw) mm = bw;
				// fprintf(stderr, "+++i->%ld, j->%ld, dd->%ld, mm->%ld\n", i, j, dd, mm);
				if(dd <= mm) {///pass distance checking
					adjust_rev_tse(lj, ug->g->seq[lj->tn].len, &jts, &jte);
					sc = get_add_cov_score(uref, jts, jte, its, ite, li->tn, &cov_i) + pop_sc(track[j]);
					
					if((sc > mm_sc) || (sc == mm_sc && mm_idx == -1)) { ///must be >= instead of >
						mm_sc = sc, mm_idx = j;
					}
					// fprintf(stderr, "<i->%ld, its:%ld, ite:%ld>, <j->%ld, jts:%ld, jte:%ld> sc:%ld, pop_sc(track[j]):%ld, csc:%ld\n", 
					// i, its, ite, j, jts, jte, sc, pop_sc(track[j]), csc);
				}
			}
		}

		// fprintf(stderr, "##(%ld) %u\t%u\t%c\tutg%.6d%c(%u)\t%u\t%u\tmm_idx:%ld\tmm_sc:%ld\n", i, li->qs, li->qe, "+-"[li->rev], 
        // (int32_t)(li->tn)+1, "lc"[uref->ug->u.a[li->tn].circ], uref->ug->u.a[li->tn].len, li->ts, li->te, mm_idx, mm_sc);
		track[i] = push_sc_pre(mm_sc, mm_idx); 
		li->sec = (mm_idx<0?0x3FFFFFFF:i-mm_idx); sv[i] = *li;
	}

	int64_t n_u;
	for (k = ch_n-1, n_u = 0; k >= 0; --k) {
		if(track[k]&((uint64_t)0x80000000)) continue;
		i = k; ch[n_u]=sv[i]; sc = pop_sc(track[i]);
		for (;i>=0;) {
			track[i] |= ((uint64_t)0x80000000);
			if(sv[i].qs < ch[n_u].qs) ch[n_u].qs = sv[i].qs;
			if(sv[i].ts < ch[n_u].ts) ch[n_u].ts = sv[i].ts;
			if(sv[i].qe > ch[n_u].qe) ch[n_u].qe = sv[i].qe;
			if(sv[i].te > ch[n_u].te) ch[n_u].te = sv[i].te;
			// ch[n_u].qn = i;//start idx of read alignment in chain
			i = pop_pre(track[i]);
		}
		adjust_rev_tse(&(ch[n_u]), ug->g->seq[ch[n_u].tn].len, &its, &ite);
		ch[n_u].ts = its; ch[n_u].te = ite; ch[n_u].sec = (sc>0x3FFFFFFF?0x3FFFFFFF:sc);
		// ch[n_u].qn += chain_offset; //start idx of read alignment in chain
		// ch[n_u].tn = k + chain_offset; //end idx of read alignment in chain
		ch[n_u].qn = k + chain_offset; //end idx of read alignment in chain
		n_u++;
	}
	for (i = 0; i < ch_n; ++i) {
		adjust_rev_tse(&(sv[i]), ug->g->seq[sv[i].tn].len, &its, &ite);
		sv[i].ts = its; sv[i].te = ite; 
		k = pop_pre(track[i]); sv[i].tn = k>=0?k+chain_offset:(uint32_t)-1;
	}
	return n_u;
}

void gen_linear_chains(kv_ul_ov_t *res, kv_ul_ov_t *buf, const ul_idx_t *uref, const ug_opt_t *uopt, int64_t bw, 
double diff_ec_ul, int64_t qlen, int64_t max_skip, glchain_t *bufg, st_mt_t *bufs)
{
	uint64_t k, l, z, an, m;
	radix_sort_ul_ov_srt_tn(res->a, res->a + res->n);
	kv_resize(ul_ov_t, *buf, res->n); buf->n = res->n;
	for (k = 1, l = m = 0; k <= res->n; k++) {
        if(k == res->n || res->a[k].tn != res->a[l].tn) {///qn <- (tn|rev)
			kv_resize(uint64_t, bufg->srt.a, k-l); 
			kv_resize(uint64_t, *bufs, k-l); 
			for (z = l; z < k; z++) res->a[z].tn>>=1;
			// fprintf(stderr, "\n*[M::%s::] %c\tutg%.6d%c(%u)\tocc:[%lu, %lu)\n", __func__, "+-"[res->a[l].rev], (int32_t)(res->a[l].tn)+1, 
			// "lc"[uref->ug->u.a[res->a[l].tn].circ], uref->ug->u.a[res->a[l].tn].len, l, k);
			an = l + linear_chain_dp(res->a+l, k-l, buf->a+l, uref, uopt, bw, diff_ec_ul, qlen, max_skip, bufg->srt.a.a, bufs->a, uref->ug, l);
			for (z = l; z < an; z++) res->a[m++] = res->a[z];
			// fprintf(stderr, "#occ:[%lu, %lu)\n", l, an);
			l = k;
		}
	}
	res->n = m;
}

void gen_end_coord(ul_ov_t *z, int64_t qlen, int64_t tlen, int64_t *r_qs, int64_t *r_qe, int64_t *r_ts, int64_t *r_te)
{
    int64_t qs, qe, ts, te, qtail, ttail;
    qs = z->qs; qe = z->qe; ts = z->ts; te = z->te;
    if(z->rev) {
        ts = tlen - z->te; te = tlen - z->ts;
    }
    
    if(qs <= ts) {
        ts -= qs; qs = 0;
    } else {
        qs -= ts; ts = 0;
    }

    qtail = qlen - qe; ttail = tlen - te;
    if(qtail <= ttail) {
        qe = qlen; te += qtail;
    }
    else {
        te = tlen; qe += ttail; 
    }

	if(r_qs) (*r_qs) = qs; if(r_qe) (*r_qe) = qe;
	if(r_ts) (*r_ts) = ts; if(r_te) (*r_te) = te;
    if(z->rev) {
        if(r_ts) (*r_ts) = tlen - te; 
		if(r_te) (*r_te) = tlen - ts;
    }
}

uint32_t is_end_check(uint32_t v, ul_ov_t *z, asg_t *g)
{
	if(v&1) {
		if(z->ts==0) return 1;
	} else {
		if(z->te==g->seq[v>>1].len) return 1;
	}
	
	return 0;
}

int64_t simple_g_chain_dp(kv_ul_ov_t *in, ul_ov_t *buf, const ul_idx_t *uref, const ug_opt_t *uopt, int64_t bw, 
double diff_ec_ul, int64_t qlen, int64_t max_skip, uint64_t *srt, uint64_t *idx, uint64_t *track)
{
	if(in->n == 0) return 0;
	uint32_t ai_v, aj_v, rev_n; ma_ug_t *ug = uref->ug;
    int64_t mm_ovlp, x, i, j, k, sc, csc, mm_sc, mm_idx, qo, share, in_n = in->n;
	int64_t iqs, iqe, its, ite, i_end, j_end;
    ul_ov_t *ai, *aj, *e_ai, *e_aj, rev_t;
	for (i = 0; i < in_n; i++) {
		gen_end_coord(&(in->a[i]), qlen, ug->u.a[in->a[i].tn].len, NULL, &iqe, NULL, NULL);
		srt[i] = iqe; srt[i] <<= 32; srt[i] |= (uint64_t)i;
	}
	radix_sort_gfa64(srt, srt+in_n); 
	for (i = 0; i < in_n; i++) buf[i] = in->a[(uint32_t)srt[i]];
	memcpy(in->a, buf, in_n *sizeof((*buf)));///all alignments have been sorted by the real end-qe

	for (i = 0; i < in_n; ++i) {
		ai = &(in->a[i]); ai_v = (ai->tn<<1)|ai->rev; e_ai = &(buf[i]); i_end = 0;
		gen_end_coord(ai, qlen, ug->u.a[ai->tn].len, &iqs, &iqe, &its, &ite);
		e_ai->qs = iqs; e_ai->qe = iqe; e_ai->ts = its; e_ai->te = ite;
		mm_ovlp = max_ovlp(uref->ug->g, ai_v^1);
		x = (e_ai->qs + mm_ovlp)*diff_ec_ul; 
        if(x < bw) x = bw; 
        x += e_ai->qs + mm_ovlp;
        if (x > qlen+1) x = qlen+1;
		x = find_ul_ov_max(i, buf, x+G_CHAIN_INDEL); 
		i_end = is_end_check(ai_v^1, ai, uref->ug->g);
		csc = mm_sc = e_ai->sec; mm_idx = -1; 
		for (j = x; j >= 0; --j) { // collect potential destination vertices
            aj = &(in->a[j]); aj_v = (aj->tn<<1)|aj->rev; e_aj = &(buf[i]); j_end = 0;
			if(e_aj->qe+G_CHAIN_INDEL <= e_ai->qs) break;//even this pair has a overlap, its length will be very small; just ignore
            if(e_aj->qs >= e_ai->qs+G_CHAIN_INDEL) continue; // lj is contained in li on the query coordinate; 128 for indel offset
			qo = infer_rovlp(e_ai, e_aj, NULL, NULL, NULL, ug); ///overlap length in query (UL read)
			if(ai_v != aj_v && get_ecov_adv(uref, uopt, ai_v^1, aj_v^1, bw, diff_ec_ul, qo, 0, &share)) {
				sc = csc + pop_sc(track[j]); j_end = is_end_check(aj_v, aj, uref->ug->g);
				if(i_end && j_end) sc -= (share>=csc?csc:share);
				if(sc > mm_sc) mm_sc = sc, mm_idx = j;
			}
		}
		track[i] = push_sc_pre(mm_sc, mm_idx);
        srt[i] = track[i]>>32; srt[i] <<= 32; srt[i] |= i;
	}

	int64_t n_v, n_u, n_v0; 
    radix_sort_gfa64(srt, srt+in_n); 
	for (k = in_n-1, n_v = n_u = 0; k >= 0; --k) {
		n_v0 = n_v;
		for (i = (uint32_t)srt[k]; i >= 0 && (track[i]&((uint64_t)0x80000000)) == 0;){
			buf[n_v] = in->a[i]; 
			gen_end_coord(&(buf[n_v]), qlen, ug->u.a[buf[n_v].tn].len, &iqs, &iqe, &its, &ite);
			buf[n_v].qs = iqs; buf[n_v].qe = iqe; buf[n_v].ts = its; buf[n_v].te = ite;

			track[i] |= ((uint64_t)0x80000000);
            i = pop_pre(track[i]);
			n_v++;
		}
		if(n_v0 == n_v) continue;
		sc = (i<0?(pop_sc(srt[k])):(pop_sc(srt[k])-pop_sc(track[i])));
		idx[n_u++] = ((uint64_t)sc<<32)|(n_v-n_v0); 
	}

	for (k = 0, n_v = n_v0 = 0; k < n_u; k++) {
        n_v0 = n_v; n_v += (uint32_t)idx[k];
        in->a[k].qn = idx[k]>>32;//score
        in->a[k].ts = n_v0; in->a[k].te = n_v;///idx
        
        rev_n = ((uint32_t)idx[k])>>1; 
        ///we need to consider contained reads; so determining qs is not such easy
        in->a[k].qs = (uint32_t)-1; in->a[k].qe = buf[n_v0].qe;
        for (i = 0; i < rev_n; i++) {
            rev_t = buf[n_v0+i]; buf[n_v0+i] = buf[n_v-i-1]; buf[n_v-i-1] = rev_t;
            
            if(in->a[k].qs > buf[n_v0+i].qs) in->a[k].qs = buf[n_v0+i].qs;
            if(in->a[k].qs > buf[n_v-i-1].qs) in->a[k].qs = buf[n_v-i-1].qs;
        }
        if(((uint32_t)idx[k])&1) {
            if(in->a[k].qs > buf[n_v0+i].qs) in->a[k].qs = buf[n_v0+i].qs;
        }
        // fprintf(stderr, "[M::%s] k:%ld, qs:%u, qe:%u, chain_occ:%u, chain_score:%u\n", __func__, k, 
        //                  res->a[k].qs, res->a[k].qe, res->a[k].te - res->a[k].ts, res->a[k].qn);
    }

	in->n = n_u;
    radix_sort_ul_ov_srt_qn(in->a, in->a + in->n);//sort by score
	return n_v;
}

/**
uint32_t uov2rov(const ul_idx_t *uref, ul_ov_t *r_al, ul_ov_t *ul_al, ul_ov_t *res)
{
	int64_t y_s, y_e, y_bs, y_be, x_s, x_e, q_s, q_e, s_shift, e_shift;
	y_s = MAX(r_al->ts, ul_al->ts); y_e = MIN(r_al->te, ul_al->te);
	if(y_s > y_e) return 0;
	res->tn = r_al->qn; res->ts = y_s; res->te = y_e; res->el = 1; res->rev = r_al->rev; res->sec = 0;
	s_shift = get_offset_adjust(y_s-r_al->ts, r_al->te-r_al->ts, r_al->qe-r_al->qs);
	e_shift = get_offset_adjust(r_al->te-y_e, r_al->te-r_al->ts, r_al->qe-r_al->qs);
	if(r_al->rev) {
		y_s = s_shift; s_shift = e_shift; e_shift = y_s;
	}
	res->qn = 0; res->qs = r_al->qs+s_shift; res->qe = r_al->qe-e_shift;
	return 1;
}


void update_ul_vec_t()
{

}

void ug2rg_gen(ul_ov_t *a, int64_t an, ul_vec_t *qn, const ul_idx_t *uref, ul_vec_t *rch)
{
	ul_ov_t *ot, p, res; uint64_t i, l, m; 
	ma_utg_t *u; uc_block_t *b; int64_t z, ff, iqs, iqe, its, ite; 


	for (z = 0; z < an; z++) {
		gen_end_coord(&(a[z]), rch->rlen, uref->ug->u.a[a[z].tn].len, &iqs, &iqe, NULL, NULL);
		o = &(a[z]); u = &(uref->ug->u.a[o->tn]);
		for (i = l = 0; i < u->n; i++) {
			p.tn = o->tn; p.rev = (u->a[i]>>32)&1; p.qn = u->a[i]>>33;///tn is unitig, qn is HiFi read
			p.qs = 0; p.qe = Get_READ_LENGTH(R_INF, (u->a[i]>>33));
			p.ts = l; p.te = l + Get_READ_LENGTH(R_INF, (u->a[i]>>33));
			l += (uint32_t)u->a[i];
			if(p.te <= o->ts) continue; 
			if(p.ts >= o->te) break;
			ff = uov2rov(uref, &p, o, &res);
			assert(ff);
			if(ff) {
				kv_pushp(uc_block_t, rch->bb, &b);
				b->hid = res.tn; b->rev = res.rev; b->base = 0; b->el = res.el;
				b->pchain = 1; b->qs = res.qs; b->qe = res.qe; b->ts = res.ts; b->te = res.te;
			}
		}
	}
}
**/



void extend_end_coord(mg_lchain_t *li, ul_ov_t *ui, const int64_t qlen, const int64_t rlen, int64_t *r_qs, int64_t *r_qe, int64_t *r_rs, int64_t *r_re)
{
    int64_t qs = 0, qe = 0, rs = 0, re = 0, rev = 0, qtail = 0, rtail = 0;
	if(li) {
		qs = li->qs; qe = li->qe; rs = li->rs; re = li->re; rev = li->v&1; 
		if(rev) {
			rs = rlen - li->re; re = rlen - li->rs;
		}
	}
	if(ui) {
		qs = ui->qs; qe = ui->qe; rs = ui->ts; re = ui->te; rev = ui->rev; 
		if(rev) {
        rs = rlen - ui->te; re = rlen - ui->ts;
    	}
	}
    
    
    
    if(qs <= rs) {
        rs -= qs; qs = 0;
    } else {
        qs -= rs; rs = 0;
    }

    qtail = qlen - qe; rtail = rlen - re;
    if(qtail <= rtail) {
        qe = qlen; re += qtail;
    }
    else {
        re = rlen; qe += rtail; 
    }

	if(r_qs) (*r_qs) = qs; if(r_qe) (*r_qe) = qe;
	if(r_rs) (*r_rs) = rs; if(r_re) (*r_re) = re;
    if(rev) {
        if(r_rs) (*r_rs) = rlen - re; 
		if(r_re) (*r_re) = rlen - rs;
    }
}

void dump_linear_chain(asg_t *g, kv_ul_ov_t *lidx, kv_ul_ov_t *autom, vec_mg_lchain_t *res, int64_t qlen)
{
	uint64_t k; int64_t iqs, iqe, its, ite;
	res->n = 0; kv_resize(mg_lchain_t, *res, lidx->n); res->n = lidx->n;
	for (k = 0; k < lidx->n; k++) {
		memset(&(res->a[k]), 0, sizeof(res->a[k]));
		// res->a[k].v = (autom->a[lidx->a[k].tn].tn<<1)|lidx->a[k].rev;
		res->a[k].v = (lidx->a[k].tn<<1)|(lidx->a[k].rev);
		///.off -> idx of original chain; cnt -> score of the chain
		res->a[k].off = k; res->a[k].score = lidx->a[k].sec;
		res->a[k].qs = lidx->a[k].qs; res->a[k].qe = lidx->a[k].qe;
		res->a[k].rs = lidx->a[k].ts; res->a[k].re = lidx->a[k].te; 
		extend_end_coord(&(res->a[k]), NULL, qlen, g->seq[res->a[k].v>>1].len, &iqs, &iqe, &its, &ite);
		res->a[k].qs = iqs; res->a[k].qe = iqe; res->a[k].rs = its; res->a[k].re = ite;

		// fprintf(stderr, "chain_id:%d\t%u\t%u\t%c\tutg%.6dl(%u)\t%u\t%u\n", 
		// res->a[k].off, res->a[k].qs, res->a[k].qe, "+-"[res->a[k].v&1], (int32_t)(res->a[k].v>>1)+1, 
		// g->seq[res->a[k].v>>1].len, res->a[k].rs, res->a[k].re);
	}
}

int64_t find_mg_lchain_max(int64_t n, const mg_lchain_t *a, int32_t x)
{
    int64_t s = 0, e = n;
    if (n == 0) return -1;
    if (a[n-1].qe < x) return n - 1;
    if (a[0].qe >= x) return -1;
    while (e > s) { // TODO: finish this block
        int64_t m = s + (e - s) / 2;
        if (a[m].qe >= x) e = m;
        else s = m + 1;
    }
    assert(s == e);
    return s;
}


int64_t hc_target_len(asg_t *g, mg_lchain_t *s, mg_lchain_t *e)
{
	// int64_t ql = s->qe - e->qe, tp, tm;
	// if((s->v^1)&1) tp = g->seq[s->v>>1].len - s->re;
	// else tp = s->rs;

	// if((e->v^1)&1) tm = g->seq[e->v>>1].len - e->re;
	// else tm = e->rs;
	int64_t ql = (int64_t)s->qs - (int64_t)e->qe, tp, tm;
	int64_t sts, ete;
	sts = (s->v&1)?g->seq[s->v>>1].len-s->re:s->rs; tp = g->seq[s->v>>1].len - sts;
	ete = (e->v&1)?g->seq[e->v>>1].len-e->rs:e->re; tm = g->seq[e->v>>1].len - ete;
	// fprintf(stderr, "[M::%s::] ql:%ld, tp:%ld, tm:%ld, sts:%ld, ete:%ld\n", __func__, ql, tp, tm, sts, ete);
	return ql + tp - tm;
}

inline int32_t cal_gchain_sc(const mg_path_dst_t *dj, const mg_lchain_t *li, const mg_lchain_t *lc, int64_t *f, int64_t b_w, float diff_thre, float chn_pen_gap)
{
	// const mg_lchain_t *lj;
	int32_t gap, sc;
	float lin_pen, log_pen;
	if (dj->n_path == 0) return INT32_MIN;
	gap = dj->dist - dj->target_dist;
	// lj = &lc[dj->meta];
	if (gap < 0) gap = -gap;
	if ((gap > ((dj->target_dist)*diff_thre)) && (gap > b_w)) return INT32_MIN;
	// if (lj->qe <= li->qs) sc = li->score;
	// else sc = (int32_t)((double)(li->qe - lj->qe) / (li->qe - li->qs) * li->score + .499); // dealing with overlap on query
	sc = li->score;
	//sc += dj->mlen; // TODO: is this line the right thing to do?
	// if (dj->is_0) sc += ref_bonus;
	lin_pen = chn_pen_gap * (float)gap;
	log_pen = gap >= 2? mg_log2(gap) : 0.0f;
	sc -= (int32_t)(lin_pen + log_pen);
	sc += f[dj->meta];
	return sc;
}

void set_trans_arr(uint64_t *trans, vec_sp_node_t *out, int64_t idx)
{
	int64_t i;
	for (i = idx; i >=0; ) {
		// if(out->a[idx]->v == 65) {
		// 	fprintf(stderr, "+[M::%s::] out->a[%ld]->v:%u\n", __func__, i, out->a[i]->v);
		// }
		trans[i]++;
		i = out->a[i]->pre;
	}
}

int64_t select_mul_way_nodes(mg_pathv_t *a, int64_t a_n, vec_sp_node_t *out, float len_dif, int32_t m_pathn, uint64_t *flag)
{
	int64_t i, k, pd, kd, occ, tt = 0; uint32_t pp; mg_pathv_t *p = NULL;
	if(a_n > 1) radix_sort_mg_pathv_t_d_srt(a, a + a_n);
	for (i = 0; i < a_n; i++) {
		p = &(a[i]); pp = (p->d<<1)>>1; pd = out->a[p->pre]->di>>32; occ = 1;
		
		if(p->d&0x80000000) break;
		for (k = 0; k < a_n; k++) {
			if(k == i) continue; ///same alignment
			if(((a[k].d<<1)>>1) == pp) continue; //same path
			if(a[k].v!=p->v) continue;
			kd = out->a[a[k].pre]->di>>32;
			if(kd >= (pd*(1-len_dif)) && kd <= (pd*(1+len_dif))) occ++;
		}
		if(occ >= m_pathn) {
			flag[p->pre] = 1; tt++;
		}
	}

	return tt;
}

int32_t phase_mul_ways(vec_mg_pathv_t *res, st_mt_t *dst_done, vec_sp_node_t *out, int32_t n_dst, mg_path_dst_t *dst, int32_t max_k, float len_dif)
{
	int64_t i, j, z, zl, n = 0, n_mpath, kk_p, od, res_n = res->n, pid; mg_pathv_t *h;
	dst_done->n = 0; kv_resize(uint64_t, *dst_done, out->n); 
	uint64_t *trans = dst_done->a; memset(dst_done->a, 0, out->n*sizeof(*(dst_done->a)));

	n_mpath = 0;
	for (i = 0; i < n_dst; ++i) { // mark dst vertices with a target distance
		mg_path_dst_t *t = &dst[i];
		if (t->n_path > 0 && t->target_dist >= 0 && t->path_end >= 0){
			assert((int32_t)(out->a[t->path_end]->di>>32) == t->target_dist);
			if(t->n_path >= max_k) {
				t->n_path = 0;
				for (z = zl = t->path_end; z >= 0;) {
					zl = z;
					z = out->a[z]->pre;
				}
				n += 2; trans[t->path_end] = trans[zl] = 1;
			} else {
				kk_p = 0;
				for (j = t->path_end; j < (int32_t)out->n; j++) {
					od = out->a[j]->di>>32;
					if(od >= (t->target_dist*(1-len_dif)) && 
							od <= (t->target_dist*(1+len_dif))) {
						if(out->a[j]->v == t->v) kk_p++;
					} else {
						break;
					}
				}
				for (j = t->path_end-1; j >=0; j--) {
					od = out->a[j]->di>>32;
					if(od >= (t->target_dist*(1-len_dif)) && 
							od <= (t->target_dist*(1+len_dif))) {
						if(out->a[j]->v == t->v) kk_p++;
					} else {
						break;
					}
				}
				assert(kk_p > 0 && kk_p <= t->n_path);
				n_mpath += kk_p;
			}
		}
	}

	// if(detect_mul_way && src == 74) {
	// 	fprintf(stderr, "+[M::%s::] src:%u, dst:%u, n_mpath:%d\n", __func__, src, dst[0].v, n_mpath);
	// }

	if(n_mpath > 1) {
		for (i = 0, pid = 0; i < n_dst; ++i) { // mark dst vertices with a target distance
			mg_path_dst_t *t = &dst[i];
			if (t->n_path > 0 && t->target_dist >= 0 && t->path_end >= 0){
				assert((int32_t)(out->a[t->path_end]->di>>32) == t->target_dist);
				assert(t->n_path < max_k);
				for (j = t->path_end; j < (int32_t)out->n; j++) {
					od = out->a[j]->di>>32;
					if(od >= (t->target_dist*(1-len_dif)) && 
							od <= (t->target_dist*(1+len_dif))) {
						if(out->a[j]->v == t->v) {
							for (z = j; z >= 0;) {
								kv_pushp(mg_pathv_t, *res, &h);
								h->v = out->a[z]->v; h->pre = z;
								h->d = pid; if(j!=t->path_end) h->d |= 0x80000000;
								z = out->a[z]->pre;
							}
							pid++;
						}
					} else {
						break;
					}
				}
				for (j = t->path_end-1; j >=0; j--) {
					od = out->a[j]->di>>32;
					if(od >= (t->target_dist*(1-len_dif)) && 
							od <= (t->target_dist*(1+len_dif))) {
						if(out->a[j]->v == t->v) {
							for (z = j; z >= 0;) {
								kv_pushp(mg_pathv_t, *res, &h);
								h->v = out->a[z]->v; h->pre = z;
								h->d = pid; if(j!=t->path_end) h->d |= 0x80000000;
								z = out->a[z]->pre;
							}
							pid++;
						}
					} else {
						break;
					}
				}
			}
		}

		radix_sort_mg_pathv_t_v_srt(res->a + res_n, res->a + res->n);
		for (i = res_n+1, j = res_n/**, n = 0**/; i <= (int64_t)res->n; ++i) {
			if (i == (int64_t)res->n || res->a[i].v != res->a[j].v) {
				n += select_mul_way_nodes(res->a + j, i - j, out, len_dif, n_mpath, trans);
				j = i;
			}
		}

		res->n = res_n;
	}


	if(n > 0) {//found some nodes
		for (i = n = 0; (uint32_t)i < out->n; ++i) { // generate coordinate translations
			if (trans[i]) {
				trans[i] = n++;
			} else {
				trans[i] = (uint32_t)-1;
			}
		}

		kv_resize(mg_pathv_t, *res, res->n + n); //res->n += n;
		for (i = 0; (uint32_t)i < out->n; ++i) { // generate the backtrack array
			mg_pathv_t *p;
			if (trans[i] == (uint32_t)-1) continue;
			p = &res->a[trans[i]+res->n];
			p->v = out->a[i]->v, p->d = out->a[i]->di >> 32;
			if(out->a[i]->pre < 0) {
				p->pre = out->a[i]->pre;
			} else {
				if(trans[out->a[i]->pre] == (uint32_t)-1) p->pre = -2;
				else p->pre = trans[out->a[i]->pre];
			}
		}

		res->n += n;
		for (i = 0; i < n_dst; ++i) // translate "path_end"
			if (dst[i].path_end >= 0)
				dst[i].path_end = trans[dst[i].path_end];
	}

	return n_mpath;
}


///max_dist is like the overlap length in string graph
///first_src_ban do not allow co-linear chain at the same node
void hc_shortest_k(void *km0, const asg_t *g, uint32_t src, int32_t n_dst, mg_path_dst_t *dst, int32_t max_dist, int32_t max_k, 
st_mt_t *dst_done, uint64_t *dst_group, vec_sp_node_t *out, vec_mg_pathv_t *res, uint64_t first_src_ban, 
uint64_t detect_mul_way, float len_dif)
{
    sp_node_t *p, *root = 0;
    sp_topk_t *q;
    khash_t(sp) *h;///
    khash_t(sp2) *h2;///alignment->vertice index
    void *km;
    khint_t k;
    int absent;
    int32_t i, j, n_done, n_found;
    uint32_t id;

    // if (res) res->n = 0;///for us, n_pathv = NULL
    if (n_dst <= 0) return;///n_dst: how many candidate vertices
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

    // multiple dst[] may have the same dst[].v. We need to group them first.
    // in other words, one ref id may have multiple dst alignment chains
	dst_done->n = 0; kv_resize(uint64_t, *dst_done, (uint64_t)n_dst);
    for (i = 0; i < n_dst; ++i) {
		dst_group[i] = ((((uint64_t)dst[i].v)<<32)|((uint64_t)i));
		dst_done->a[i] = 0;
	}
        
    radix_sort_gfa64(dst_group, dst_group + n_dst);

    h2 = kh_init2(sp2, km); // (h2+dst_group) keeps all destinations from the same ref id
    kh_resize(sp2, h2, n_dst * 2);
    ///please note that one contig in ref may have multiple alignment chains
    ///so h2 is a index that helps us to query it
    ///key(h2) = ref id; value(h2) = start_idx | occ
    for (i = 1, j = 0; i <= n_dst; ++i) {
        if (i == n_dst || dst_group[i]>>32 != dst_group[j]>>32) {
            k = kh_put(sp2, h2, dst_group[j]>>32, &absent);
            kh_val(h2, k) = (((uint64_t)j)<<32)|((uint64_t)(i-j));
            assert(absent);
            j = i;
        }
    }

    h = kh_init2(sp, km); // h keeps visited vertices; path to each visited vertice
    kh_resize(sp, h, 16);

	out->n = 0; kv_resize(sp_node_t*, *out, 16); ///16 is just the initial size    
    id = 0;
    p = gen_sp_node(km, src, 0, id++);///just malloc a node for src; the distance is 0
    p->hash = __ac_Wang_hash(src);///hash is path hash, instead of node hash
    kavl_insert(sp, &root, p, 0);///should be avl tree; p is a node at avl-tree

    ///each cell in the hash table <h> corresponds to one node in the graph
	///each cell in the AVL tree <root> is a path, corresponds to <MG_MAX_SHORT_K> node in the graph
    k = kh_put(sp, h, src, &absent);
    q = &kh_val(h, k);
    ///for normal graph traversal, one node just has one parental node; here each node has at most 16 parental nodes
    q->k = 1, q->p[0] = p, q->mlen = 0, q->qs = q->qe = -1;

    n_done = 0; first_src_ban = first_src_ban?0:1;
    ///the key of avl tree: #define sp_node_cmp(a, b) (((a)->di > (b)->di) - ((a)->di < (b)->di))
    ///the higher bits of (*)->di is distance to src node
    ///so the key of avl tree is distance
    ///in avl tree <root>, one node might be saved multipe times
    while (kavl_size(head, root) > 0) {///thr first root is src
        int32_t i, nv;
        asg_arc_t *av;
        sp_node_t *r;
        ///note that one node in the graph (sp_node_t->v) might be visited multiple times if there are circles
        ///so there might be multipe cells in the avl-tree with the same (sp_node_t->v)
        ///delete the first cell
        r = kavl_erase_first(sp, &root); // take out the closest vertex in the heap (as a binary tree)
        //fprintf(stderr, "XX\t%d\t%d\t%d\t%c%s[%d]\t%d\n", n_out, kavl_size(head, root), n_finished, "><"[(r->v&1)^1], g->seg[r->v>>1].name, r->v, (int32_t)(r->di>>32));
        
		///higher 32 bits might be the distance to root node
        // lower 32 bits now for position in the out[] array
		///r->pre keep the pre-node in the path; follow the pre it is able to recover the whole path
        r->di = ((r->di>>32)<<32)|((uint64_t)out->n); ///n_out is just the id in out
        ///so one node id in graph might be saved multiple times in avl tree and out[]
		kv_push(sp_node_t*, *out, r);

        ///r->v is the dst vertex id
        ///sometimes k==kh_end(h2). Some nodes are found by graph travesal but not in linear chain alignment
        k = kh_get(sp2, h2, r->v);
        // we have reached one dst vertex
        // note that one dst vertex may have multipe alignment chains
        // we can visit some nodes in graph which are not reachable during chaining
        // h2 is used to determine if one node is reachable or not
		// if(src == 2844) {
		// 	fprintf(stderr, "******src->%u, dst->%u, max_dist->%d\n", src, r->v, max_dist);
		// }
        if (k != kh_end(h2) && first_src_ban) { 
            ///node r->v might be visited multiple times
            int32_t j, dist = r->di>>32, off = kh_val(h2, k) >> 32, cnt = (int32_t)kh_val(h2, k);
			// if(src == 2844) {
			// 	fprintf(stderr, "----src->%u, dst->%u, max_dist->%d, cnt->%d\n", src, r->v, max_dist, cnt);
			// }
            //src can reach ref id r->v; there might be not only one alignment chain in r->v
            //so we need to scan all of them
            for (j = 0; j < cnt; ++j) {
                mg_path_dst_t *t = &dst[(int32_t)dst_group[off + j]];///t is a linear alignment at r->v
                int32_t done = 0;
				// if((src>>1) == 51) {
				// 	fprintf(stderr, "###src->%u, dst->%u, max_dist->%d, dist:%d\n", src, r->v, max_dist, dist);
				// }
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
				   // means this alignment has never been visited before; keep it anyway
				   // note here is the alignment, instead of node
					
				    // if(src == 2844) {
					// 	fprintf(stderr, ">>src->%u, dst->%u, target_dist->%d, dist->%d, max_dist->%d\n", 
					// 	src, r->v, t->target_dist, dist, max_dist);
					// }
                    if (t->n_path == 0) { 
                        copy = 1;
                     // we have a target distance; choose the closest; 
                     // there is already several paths reaching the linear alignment <t>
                    } else if (t->target_dist >= 0) { 
                        // we found the target path; hash is the path hash including multiple nodes, instead of node hash 
                        if (dist == t->target_dist && t->check_hash && r->hash == t->target_hash) { 
                            copy = 1, done = 1;
                        } else {
                            int32_t d0 = t->dist, d1 = dist;
                            d0 = d0 > t->target_dist? d0 - t->target_dist : t->target_dist - d0;
                            d1 = d1 > t->target_dist? d1 - t->target_dist : t->target_dist - d1;
                            ///if the new distance (d1) is smaller than the old distance (d0), update the results
                            ///the length of new path should be closer to t->target_dist
                            if (d1 - mlen/2 < d0 - t->mlen/2) copy = 1;
                        }
                    }
                    if (copy) {
                        t->path_end = out->n-1, t->dist = dist, t->hash = r->hash, t->mlen = mlen, t->is_0 = r->is_0;
                        if (t->target_dist >= 0) {
                            ///src is from li from li to lj, so the dis is generally increased; dijkstra algorithm
                            ///target_dist should be the distance on query
                            if (dist == t->target_dist && t->check_hash && r->hash == t->target_hash) {
								done = 1;
							} else if ((dist > t->target_dist + MG_SHORT_K_EXT) && (dist > (t->target_dist>>4)) && 
								(dist > (t->target_dist*1.25))) {
								done = 1;
							}
                        }
                    }
                    ++t->n_path;///we found a path to the alignment t
                    if (t->n_path >= max_k) done = 1;
                }
                if (detect_mul_way == 0 && dst_done->a[off + j] == 0 && done)
                    dst_done->a[off + j] = 1, ++n_done;
            }
            ///if all alignments have been settle down
            ///pre-end; accelerate the loop
            if (n_done == n_dst) break;
        }
		first_src_ban = 1;
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
                p->pre = out->n - 1;///the parent node of this one 
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
                    p->pre = out->n - 1;
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
                    return;
                }
            } // else: the path is longer than all the existing paths ended at ai->w
        }
    }
    kh_destroy(sp, h);
    // NB: AVL nodes are not deallocated. When km==0, they are memory leaks.

    for (i = 0, n_found = 0; i < n_dst; ++i)
        if (dst[i].n_path > 0) ++n_found;///n_path might be larger than 16
    ///we can assume n_pathv = NULL for now 
    if (n_found > 0 && res) { // then generate the backtrack array
        int32_t n, n_mpath = 1; dst_done->n = 0; kv_resize(uint64_t, *dst_done, out->n); 
		uint64_t *trans = dst_done->a; memset(dst_done->a, 0, out->n*sizeof(*(dst_done->a)));

		if(detect_mul_way) {
			n_mpath = phase_mul_ways(res, dst_done, out, n_dst, dst, max_k, len_dif);
		}
		
		if(n_mpath == 1) {
			// KCALLOC(km, trans, n_out); // used to squeeze unused elements in out[]
			///n_out: how many times that nodes in graph have been visited 
			///note one node might be visited multiples times
			///n_dst: number of alignment chains
			for (i = 0; i < n_dst; ++i) { // mark dst vertices with a target distance
				mg_path_dst_t *t = &dst[i];
				if (t->n_path > 0 && t->target_dist >= 0 && t->path_end >= 0)
					trans[(uint32_t)out->a[t->path_end]->di] = 1;///(int32_t)out[]->di: traverse track corresponds to the alignment chain dst[]
			}
			// for (i = 0; (uint32_t)i < out->n; ++i) { // mark dst vertices without a target distance
			//     k = kh_get(sp2, h2, out->a[i]->v);
			//     if (k != kh_end(h2)) { // TODO: check if this is correct!
			//         int32_t off = kh_val(h2, k)>>32, cnt = (int32_t)kh_val(h2, k);
			//         for (j = off; j < off + cnt; ++j)
			//             if (dst[j].target_dist < 0)
			//                 trans[i] = 1;
			//     }
			// }
			for (i = (int32_t)(out->n) - 1; i >= 0; --i) // mark all predecessors
				if (trans[i] && out->a[i]->pre >= 0)
					trans[out->a[i]->pre] = 1;
			for (i = n = 0; (uint32_t)i < out->n; ++i) // generate coordinate translations
				if (trans[i]) trans[i] = n++;
				else trans[i] = (uint32_t)-1;
			
			kv_resize(mg_pathv_t, *res, res->n + n); //res->n += n;
			for (i = 0; (uint32_t)i < out->n; ++i) { // generate the backtrack array
				mg_pathv_t *p;
				if (trans[i] == (uint32_t)-1) continue;
				p = &res->a[trans[i]+res->n];
				p->v = out->a[i]->v, p->d = out->a[i]->di >> 32;
				p->pre = out->a[i]->pre < 0? out->a[i]->pre:trans[out->a[i]->pre];
			}
			res->n += n;
			for (i = 0; i < n_dst; ++i) // translate "path_end"
				if (dst[i].path_end >= 0)
					dst[i].path_end = trans[dst[i].path_end];
		}
    }

    km_destroy(km);
}

///p[]: id of last
///f[]: the score ending at i, not always the peak
///v[]: keeps the peak score up to i;
///t[]: used for buffer
///min_cnt = 2; min_sc = 30; extra_u = 0
///u = mg_chain_backtrack(n, f, p, v, t, min_cnt, min_sc, 0, &n_u, &n_v);
int64_t hc_chain_backtrack(int64_t n, const int64_t *f, const uint64_t *p, uint64_t *srt, uint64_t *u, uint64_t *v, 
int64_t *n_u_, int64_t *n_v_)
{
	if(n_u_) *n_u_ = 0; if(n_v_) *n_v_ = 0;
	int64_t i, k, n_v, n_srt, n_v0, n_u, sc;
	if (n == 0) return 0;
	// v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
	*n_u_ = *n_v_ = 0;	
	for (i = 0, k = 0; i < n; ++i) {
		if(f[i] >= 0) {
			srt[k] = (uint64_t)f[i]; srt[k] <<= 32; srt[k] |= (((uint64_t)i)<<1); k++;
		}
	}
	n_srt = k;
	radix_sort_gfa64(srt, srt + n_srt); ///sort by score

	///from the largest to the smallest
	for (k = n_srt-1, n_v = n_u = 0; k >= 0; --k) { // precompute n_u
		n_v0 = n_v;
		for (i = ((uint32_t)srt[k])>>1; i >= 0 && (srt[i]&1) == 0; i = (p[i]==(uint64_t)-1?-1:p[i])) {
			v[n_v++] = i; srt[i] |= 1;
		}
		if(n_v <= n_v0) continue;
		sc = i < 0? srt[k]>>32: (int64_t)(srt[k]>>32)-f[i]; 
		u[n_u++] = (((uint64_t)sc)<<32) | ((uint64_t)(n_v-n_v0));
	}

	if(n_u_) *n_u_ = n_u; if(n_v_) *n_v_ = n_v;
	return n_u;
}

void set_ul_ov_t_by_mg_lchain_t(ul_ov_t *u, mg_lchain_t *l)
{
	u->tn = l->v>>1; u->rev = (l->v&1); 
	u->ts = l->rs; u->te = l->re; 
	u->qs = l->qs; u->qe = l->qe;
}

uint64_t primary_chain_check(uint64_t *idx, int64_t idx_n, mg_lchain_t *a)
{
	if(idx_n <= 0) return 0;
	ul_ov_t m; memset(&m, 0, sizeof(m)); int64_t m_sc = -1, i, a_n; uint64_t s_idx, e_idx, ovlp, novlp;
	for (i = a_n = 0; i < idx_n; ++i) {
        if(((int64_t)(idx[i]>>32)) > m_sc) {
            m_sc = ((int64_t)(idx[i]>>32)); m.qn = i;
            m.ts = a_n; m.te = a_n + ((uint32_t)idx[i]);
            m.qs = a[m.ts].qs; m.qe = a[m.te-1].qe;
        }
        a_n += ((uint32_t)idx[i]);
    }
	assert(a[m.ts].qs<=a[m.te-1].qs && a[m.te-1].qe>=a[m.ts].qe);

	for (i = a_n = 0; i < idx_n; ++i) {
		s_idx = a[a_n].qs; a_n += ((uint32_t)idx[i]); e_idx = a[a_n-1].qe;
		if(i == m.qn) continue;
		ovlp = ((MIN(m.qe, e_idx) > MAX(m.qs, s_idx))? (MIN(m.qe, e_idx) - MAX(m.qs, s_idx)):0);
		novlp = (e_idx - s_idx) - ovlp;
		if(novlp > ((m.qe-m.qs)*GC_OFFSET_RATE) && novlp > GC_OFFSET_POS) break;
	}
	if(i >= idx_n) return 1;
	return 0;
}



int64_t hc_gchain1_dp(void *km, const ul_idx_t *uref, const ma_ug_t *ug, vec_mg_lchain_t *lc, vec_mg_lchain_t *sw, vec_mg_path_dst_t *dst, vec_sp_node_t *out, vec_mg_pathv_t *path,
int64_t qlen, const ug_opt_t *uopt, int64_t bw, double diff_thre, double ng_diff_thre, uint64_t *srt, st_mt_t *bf, int64_t *f, uint64_t *p, uint64_t *v)
{
	bf->n = 0;
	if(lc->n == 0) return 0;
	int64_t i, j, lc_n = lc->n, n_ext, mm_ovlp, target_dist, max_target_dist, x, m_idx, m_sc, qo, sc; 
	int64_t max_f, max_j = -1, max_d = -1, max_inner = 0, share; uint32_t max_hash = 0; int64_t k, k0, n_u, n_v, ni;
	mg_lchain_t *r, *li, *lj; mg_path_dst_t *q; asg_t *g = ug->g; uint64_t isolated, *u; ul_ov_t ui, uj;
	for (i = n_ext = 0; i < lc_n; i++) {
		r = &lc->a[i]; r->dist_pre = -1; isolated = 0;///dist_pre -> parent in graph chain
		if((r->re < g->seq[r->v>>1].len) && (r->rs > 0)) isolated = 1;///UL contained in one vertice
		if (!isolated) ++n_ext;
		srt[i] = r->qe; srt[i] <<= 32; srt[i] |= (uint64_t)i; srt[i] |= (isolated<<63);
	}
	radix_sort_gfa64(srt, srt+lc_n); 
	for (i = 1, j = 0; i <= lc_n; i++) {
		if (i == lc_n || (srt[i]>>32) != (srt[j]>>32)) {
			if(i - j > 1) {
				for (x = j; x < i; x++) {
					srt[x] <<= 32; srt[x] >>= 32; srt[x] |= ((uint64_t)lc->a[(uint32_t)srt[x]].qs)<<32;
				}
				radix_sort_gfa64(srt+j, srt+i);
			}
			j = i;
		}
	}

	kv_resize(mg_lchain_t, *sw, (uint64_t)lc_n); sw->n = lc_n;
	for (i = 0; i < lc_n; i++) sw->a[i] = lc->a[(uint32_t)srt[i]];
	memcpy(lc->a, sw->a, lc_n *sizeof((*(lc->a))));
	// fprintf(stderr, "[M::%s::] n_ext:%ld, lc_n:%ld\n", __func__, n_ext, lc_n);

	if(ng_diff_thre >= 0) {
		//first non-gap chain
		for (i = 0; i < n_ext; ++i) { // core loop
			li = &lc->a[i]; set_ul_ov_t_by_mg_lchain_t(&ui, li);
			mm_ovlp = max_ovlp(g, li->v^1);
			x = (li->qs + mm_ovlp)*diff_thre; if(x < bw) x = bw; 
			x += li->qs + mm_ovlp;
			if (x > qlen+1) x = qlen+1;
			x = find_mg_lchain_max(i, lc->a, x+G_CHAIN_INDEL);
			// fprintf(stderr, "\nli->(%ld)\tutg%.6d%c(%u)\tqs:%u\tqe:%u\t%c\trs:%u\tre:%u\tsrc:%u\tscore:%d, x:%ld\n", 
			// 	i, (int32_t)(li->v>>1)+1, "lc"[ug->u.a[li->v>>1].circ], ug->u.a[li->v>>1].len,
			// 	li->qs, li->qe, "+-"[li->v&1], li->rs, li->re, li->v^1, li->score, x);
			max_f = li->score, max_j = -1;
			// collect potential destination vertices
			for (j = x; j >= 0; --j) { 
				lj = &lc->a[j]; ///extend_end_coord(lj, qlen, g->seq[lj->v>>1].len, &jqs, &jqe, &jrs, &jre);
				//even this pair has a overlap, its length will be very small; just ignore; only for non-gapped chains
				if(lj->qe+G_CHAIN_INDEL <= li->qs) break;
				if(lj->qs >= li->qs/**+G_CHAIN_INDEL**/) continue;
				set_ul_ov_t_by_mg_lchain_t(&uj, lj);
				qo = infer_rovlp(&ui, &uj, NULL, NULL, NULL, (ma_ug_t *)ug); ///overlap length in query (UL read)
				if(li->v!=lj->v && get_ecov_adv(uref, uopt, li->v^1, lj->v^1, bw, ng_diff_thre, qo, 0, &share)) {
					sc = li->score + f[j];
					if(sc > max_f) {
						max_f = sc; max_j = j;
					}
				}
			}

			f[i] = max_f, p[i] = max_j<0?(uint64_t)-1:max_j;
			li->dist_pre = max_j<0?-1:g_adjacent_dis(g, li->v^1, lc->a[max_j].v^1); li->inner_pre = 0;
			li->hash_pre = max_j<0?0:(__ac_Wang_hash((li->v^1))+__ac_Wang_hash((lc->a[max_j].v^1))); 
			// fprintf(stderr, "i->%ld, utg%.6d%c->utg%.6d%c, max_f:%ld\n", i, (int32_t)(li->v>>1)+1, "lc"[ug->u.a[li->v>>1].circ], 
			// max_j<0?0:(int32_t)(lc->a[max_j].v>>1)+1, max_j<0?'*':"lc"[ug->u.a[lc->a[max_j].v>>1].circ], max_f);
		}

		kv_resize(uint64_t, *bf, (uint64_t)lc_n); u = bf->a;
		hc_chain_backtrack(n_ext, f, p, srt, u, v, &n_u, &n_v);
		for (i = 0; i < lc_n - n_ext; ++i) {
			u[n_u++] = (((uint64_t)lc->a[n_ext + i].score)<<32) | 1;
			v[n_v++] = n_ext + i;
		}

		sw->n = 0; kv_resize(mg_lchain_t, *sw, (uint64_t)n_v); m_idx = m_sc = -1; bf->n = 0;
		for (i = 0, k = 0; i < n_u; ++i) {
			k0 = k, ni = (int32_t)u[i];
			for (j = 0; j < ni; ++j) {
				sw->a[k++] = lc->a[v[k0 + (ni - j - 1)]];
			}
			if(m_idx < 0 || m_sc < ((int64_t)(u[i]>>32))) {
				m_idx = i; m_sc = ((int64_t)(u[i]>>32));
			}
		}
		assert(k == n_v); bf->n = n_u;

		if(primary_chain_check(u, n_u, sw->a)) {
			memcpy(lc->a, sw->a, n_v*sizeof(mg_lchain_t));
			return m_idx;
		} 
	}
	bf->n = 0;

	///then gapped-chaining
	for (i = 0; i < n_ext; ++i) { // core loop
		li = &lc->a[i]; 
		mm_ovlp = max_ovlp(g, li->v^1);
		x = (li->qs + mm_ovlp)*diff_thre; if(x < bw) x = bw; 
        x += li->qs + mm_ovlp;
        if (x > qlen+1) x = qlen+1;
		x = find_mg_lchain_max(i, lc->a, x+G_CHAIN_INDEL);
		// fprintf(stderr, "\nli->(%ld)\tutg%.6d%c(%u)\tqs:%u\tqe:%u\t%c\trs:%u\tre:%u\tsrc:%u\tscore:%d, x:%ld\n", 
		// 	i, (int32_t)(li->v>>1)+1, "lc"[ug->u.a[li->v>>1].circ], ug->u.a[li->v>>1].len,
		// 	li->qs, li->qe, "+-"[li->v&1], li->rs, li->re, li->v^1, li->score, x);

		// collect potential destination vertices
		for (dst->n = 0, max_target_dist= -1, j = x; j >= 0; --j) { 
			lj = &lc->a[j]; ///extend_end_coord(lj, qlen, g->seq[lj->v>>1].len, &jqs, &jqe, &jrs, &jre);
			//lj contained in li; actually in circle, this might happen; need to deal with it later
			if(lj->qs >= li->qs/**+G_CHAIN_INDEL**/) continue;
			///if there is a circle, the two linear chains might be at the same vertice
			target_dist = hc_target_len(g, li, lj);
			// fprintf(stderr, "j:%ld, target_dist:%ld\n", j, target_dist);
			if(target_dist < 0) continue;
			kv_pushp(mg_path_dst_t, *dst, &q);
			memset(q, 0, sizeof(*q));
			q->inner = 0;//we set q->inner = 0 to allow circles
			q->v = lj->v^1;///must be v^1 instead of v
			q->meta = j;
			///lj->qs************lj->qe
			///			li->qs************li->qe
			q->qlen = li->qs - lj->qe;///might be negative; this is the region that need to be checked in base-level
			q->target_dist = target_dist;///cannot understand the target_dist
			q->target_hash = 0;
			q->check_hash = 0;
			if(max_target_dist < target_dist) max_target_dist = target_dist;
			///not sure how to use this cut-off
			// if (t[j] == i) {
			// 	if (++n_skip > max_skip)
			// 		break;
			// }
			// if (p[j] >= 0) t[p[j]] = i;
			// if((li->v>>1)==10 && ((lj->v>>1)==15||(lj->v>>1)==14)) max_target_dist = 100000;
			// fprintf(stderr, "+++lj->(%ld)\tutg%.6d%c(%u)\t%u\t%u\t%c\ttarget_dist:%d\n", 
			// j, (int32_t)(lj->v>>1)+1, "lc"[ug->u.a[lj->v>>1].circ], ug->u.a[lj->v>>1].len,
			// lj->qs, lj->qe, "+-"[lj->v&1], q->target_dist);
		}

		// confirm reach-ability
		max_f = li->score, max_j = -1, max_d = -1, max_inner = 0; max_hash = 0;
		if(dst->n) {
			max_target_dist *= (1+diff_thre); if(max_target_dist < bw) max_target_dist = bw;
			hc_shortest_k(km, g, li->v^1, dst->n, dst->a, max_target_dist, MG_MAX_SHORT_K, bf, srt, out, NULL, 1, 0, 0);
			// remove unreachable destinations
			//TODO: check sequence identity
            for (j = 0; j < (int64_t)dst->n; ++j) {
                mg_path_dst_t *dj = &dst->a[j];
                if (dj->n_path == 0) continue; // unreachable
                sc = cal_gchain_sc(dj, li, lc->a, f, bw, diff_thre, W_CHN_PEN_GAP);		

				// fprintf(stderr, "---dj->(%ld)\tutg%.6d%c(%u)\tsc:%d\tmax_f:%ld\ttarget_dist:%d\tdj->dist:%d\n", 
				// j, (int32_t)(dj->v>>1)+1, "lc"[ug->u.a[dj->v>>1].circ], ug->u.a[dj->v>>1].len, sc, max_f, dj->target_dist, dj->dist);	
                if (sc == INT32_MIN) continue; // out of band
				// fprintf(stderr, "+max_f->%d, max_j->%d\n", max_f, max_j);
				if (sc < 0) continue;// negative score
				// fprintf(stderr, "++max_f->%d, max_j->%d\n", max_f, max_j);
				if (sc > max_f) {
					max_f = sc, max_j = dj->meta, max_d = dj->dist, max_hash = dj->hash, max_inner = dj->inner;
					// fprintf(stderr, "+++max_f->%d, max_j->%d\n", max_f, max_j);
				}
            }
		}

		f[i] = max_f, p[i] = max_j<0?(uint64_t)-1:max_j;
		li->dist_pre = max_d;
		li->hash_pre = max_hash;
		li->inner_pre = max_inner;
		// fprintf(stderr, "i->%ld, utg%.6d%c->utg%.6d%c, max_f:%ld\n", i, (int32_t)(li->v>>1)+1, "lc"[ug->u.a[li->v>>1].circ], 
		// max_j<0?0:(int32_t)(lc->a[max_j].v>>1)+1, max_j<0?'*':"lc"[ug->u.a[lc->a[max_j].v>>1].circ], max_f);
	}

	kv_resize(uint64_t, *bf, (uint64_t)lc_n); u = bf->a;
	hc_chain_backtrack(n_ext, f, p, srt, u, v, &n_u, &n_v);
	for (i = 0; i < lc_n - n_ext; ++i) {
        u[n_u++] = (((uint64_t)lc->a[n_ext + i].score)<<32) | 1;
        v[n_v++] = n_ext + i;
    }

	sw->n = 0; kv_resize(mg_lchain_t, *sw, (uint64_t)n_v); m_idx = m_sc = -1; bf->n = 0;
	for (i = 0, k = 0; i < n_u; ++i) {
        k0 = k, ni = (int32_t)u[i];
        for (j = 0; j < ni; ++j) {
			sw->a[k++] = lc->a[v[k0 + (ni - j - 1)]];
		}
		if(m_idx < 0 || m_sc < ((int64_t)(u[i]>>32))) {
			m_idx = i; m_sc = ((int64_t)(u[i]>>32));
		}
    }
    assert(k == n_v); bf->n = n_u;
	memcpy(lc->a, sw->a, n_v*sizeof(mg_lchain_t));
	return m_idx;
}


void debug_gchain(void *km, const asg_t *g, mg_lchain_t *a, uint64_t n, st_mt_t *dst_done, vec_sp_node_t *out)
{
	uint64_t k, i, v, w, nv; int64_t dd; asg_arc_t *av; mg_path_dst_t dst; uint64_t dst_group;
	for (k = 1; k < n; k++) {
		v = a[k].v^1; w = a[k-1].v^1;
		nv = asg_arc_n(g, v); av = asg_arc_a(g, v);
		for (i = 0; i < nv; i++) {
			if(av[i].v == w) break;
		}
		// if(i >= nv) {
		// 	// fprintf(stderr, "[M::%s::]\n", __func__);
		// 	fprintf(stderr, "[M::%s::]\tutg%.6dl(%c)\t->\tutg%.6dl(%c)\n", __func__, (int32_t)(v>>1)+1, "+-"[v&1], (int32_t)(w>>1)+1, "+-"[w&1]);
		// }
		if(i < nv) {
			dd = (int64_t)((uint32_t)(av[i].ul));
		} else {
			memset(&dst, 0, sizeof(dst));
            dst.v = w; 
            dst.target_dist = a[k-1].dist_pre;
            dst.target_hash = 0; dst.check_hash = 0;
			hc_shortest_k(km, g, v, 1, &dst, dst.target_dist, MG_MAX_SHORT_K, dst_done, &dst_group, out, NULL, 1, 0, 0);
			dd = dst.dist;
		}
		if(a[k-1].dist_pre != dd) {
			fprintf(stderr, "[M::%s::]\tutg%.6dl(%c)\t->\tutg%.6dl(%c)\tdist_pre:%d\td:%ld\n", __func__, (int32_t)(v>>1)+1, "+-"[v&1], 
			(int32_t)(w>>1)+1, "+-"[w&1], a[k-1].dist_pre, dd);
		}
	}
}

void debug_gchain2(const asg_t *g, mg_pathv_t *a, uint64_t n)
{
	uint64_t k, i, v, w, nv; asg_arc_t *av;
	for (k = 1; k < n; k++) {
		v = a[k-1].v; w = a[k].v;
		nv = asg_arc_n(g, v); av = asg_arc_a(g, v);
		for (i = 0; i < nv; i++) {
			if(av[i].v == w) break;
		}
		if(i >= nv) {
			// fprintf(stderr, "[M::%s::]\n", __func__);
			fprintf(stderr, "[M::%s::]\tutg%.6dl(%c)\t->\tutg%.6dl(%c)\n", __func__, (int32_t)(v>>1)+1, "+-"[v&1], (int32_t)(w>>1)+1, "+-"[w&1]);
		}
	}
}


void reverse_track(mg_pathv_t *a, uint64_t a_n)
{
	int64_t k, hn = (a_n>>1); mg_pathv_t z;
	for (k = a_n-1; k >= 1; k--) a[k].d -= a[k-1].d;

	for (k = 0; k < hn; k++) {
		z = a[k]; a[k] = a[a_n - k - 1]; a[a_n - k - 1] = z;
		a[k].v ^= 1; a[a_n - k - 1].v ^= 1;
	}
	if(a_n&1) a[k].v ^= 1;
}

void dbg_print(mg_pathv_t *a, int64_t a_n)
{
	int64_t k;
	for (k = 0; k < a_n; k++) {
		if(a[k].pre == -2) break;
	}
	
	if(k < a_n) {
		fprintf(stderr, "+[M::%s::] src:utg%.6dl(v:%u), dst:utg%.6dl(v:%u)\n", __func__,
		(int32_t)(a[0].v>>1)+1, a[0].v, (int32_t)(a[a_n-1].v>>1)+1, a[a_n-1].v);
		for (k = 0; k < a_n; k++) {
			fprintf(stderr, "-[M::%s::] utg%.6dl(v:%u), pre:%d, d:%u\n", __func__, (int32_t)(a[k].v>>1)+1, a[k].v, 
			a[k].pre, a[k].d);
		}
	}
	
}

uint32_t gen_gchain_track(void *km, mg_lchain_t *a, int64_t a_n, const asg_t *g, 
st_mt_t *dst_done, vec_sp_node_t *out, vec_mg_pathv_t *res)
{
	int64_t k, p_n/**, trav_occ = 0**/; mg_lchain_t *l0, *l1; mg_path_dst_t dst; uint64_t dst_group; mg_pathv_t *p;
	res->n = 0; kv_pushp(mg_pathv_t, *res, &p); p->v = (uint32_t)-1; p->d = 0; p->pre = 0;
	for (k = 1; k < a_n; k++) {
		l0 = a + k - 1; l1 = a + k;
		assert(!l1->inner_pre); assert(l1->dist_pre >= 0);
		memset(&dst, 0, sizeof(dst));
		dst.v = l0->v^1; 
		assert(l1->dist_pre >= 0);
		dst.target_dist = l1->dist_pre;
		dst.target_hash = l1->hash_pre;
		dst.check_hash = 1; p_n = res->n;
		if((dst.target_hash != (__ac_Wang_hash((l1->v^1))+__ac_Wang_hash(dst.v))) || 
													(g_adjacent_dis(g, l1->v^1, dst.v) != dst.target_dist)) {
			hc_shortest_k(km, g, l1->v^1, 1, &dst, dst.target_dist*(1+SEC_LEN_DIF), MG_MAX_SHORT_K, dst_done, &dst_group, out, res, 1, 1, SEC_LEN_DIF);
			// debug_gchain2(g, res->a + p_n, res->n - p_n);
			// fprintf(stderr, "[M::%s::n->%ld]\tutg%.6dl(%c)\t->\tutg%.6dl(%c)\n", __func__, res->n - p_n, 
			// (int32_t)(l0->v>>1)+1, "+-"[l0->v&1], (int32_t)(l1->v>>1)+1, "+-"[l1->v&1]);
			
			// fprintf(stderr, "\n-[M::%s::res->n->%u::p_n->%ld] utg%.6dl(v:%u) -> utg%.6dl(v:%u)\n", 
			// 			__func__, (uint32_t)res->n, p_n, (int32_t)(l1->v>>1)+1, l1->v, (int32_t)(l0->v>>1)+1, l0->v);
			// dbg_print(res->a + p_n, res->n - p_n);
			assert(res->n - p_n > 1); assert(dst.target_hash == dst.hash); 
			res->a[p_n-1].d = res->a[res->n-1].d - res->a[res->n-2].d; res->n--; 
			reverse_track(res->a + p_n, res->n - p_n); res->n--;///reomve l1 from res
			// trav_occ++;
		} else {
			res->a[p_n-1].d = dst.target_dist;
		}

		kv_pushp(mg_pathv_t, *res, &p); p->v = (uint32_t)-1; p->pre = k; p->d = 0;
		
	}
	// fprintf(stderr, "[M::%s::]\ta_n:%ld\ttrav_occ:%ld\n", __func__, a_n, trav_occ);
	return res->n;
}

void print_chain(mg_lchain_t *a, uint32_t a_n)
{
	uint32_t k;
	for (k = 0; k < a_n; k++) {
		if(a[k].off!=-1) {
			fprintf(stderr, "%u\t%u\t%c\tutg%.6dl\t%u\t%u\n", 
					a[k].qs, a[k].qe, "+-"[a[k].v&1], (int32_t)(a[k].v>>1)+1, a[k].rs, a[k].re);
		} else {
			fprintf(stderr, "*\t*\t%c\tutg%.6dl\t*\t*\n", "+-"[a[k].v&1], (int32_t)(a[k].v>>1)+1);
		}
	}
}

void update_exist_chain(const ul_idx_t *uref, ul_ov_t *ch, uint64_t *idx, int64_t idx_n, int64_t tid, int64_t bw, 
double diff_ec_ul, mg_lchain_t *res)
{
	int64_t i, j, cov_i, i_qs, i_qe, i_ts, i_te, j_qs, j_qe, j_ts, j_te, dq, dt, dd, mm, sc = 0; 
	int64_t tlen = uref->ug->g->seq[tid].len; memset(res, 0, sizeof(*res));
	ul_ov_t *li, *lj;
	if(idx_n <= 0) return;
	i = idx_n - 1; 
	res->qs = (ch[idx[i]].qs<<1)>>1; res->qe = ch[idx[i]].qe;
	res->rs = ch[idx[i]].ts; res->re = ch[idx[i]].te;
	for (; i >= 0; i--) {
		li = &(ch[idx[i]]); 
		i_qs = (li->qs<<1)>>1; i_qe = li->qe; 
		i_ts = (li->rev?(tlen-li->te):(li->ts));
		i_te = (li->rev?(tlen-li->ts):(li->te));

		if((int64_t)((li->qs<<1)>>1) < res->qs) res->qs = ((li->qs<<1)>>1); 
		if((int64_t)li->ts < res->rs) res->rs = li->ts;
		if((int64_t)li->qe > res->qe) res->qe = li->qe;
		if((int64_t)li->te > res->re) res->re = li->te;

		cov_i = 0; j = i + 1;
		if(j < idx_n) {
			lj = &(ch[idx[j]]); 
			j_qs = (lj->qs<<1)>>1; j_qe = lj->qe; 
			j_ts = (lj->rev?(tlen-lj->te):(lj->ts));
			j_te = (lj->rev?(tlen-lj->ts):(lj->te));
			if(j_qs <= i_qs && j_qe <= i_qe && j_ts <= i_ts && j_te <= i_te) {///co-linear
				assert(li->rev == lj->rev);
				if(j_qs == i_qs && j_qe == i_qe && j_ts == i_ts && j_te == i_te) continue;
				dq = i_qe - j_qs; dt = i_te - j_ts;
                dd = (dq>dt? dq-dt:dt-dq);
                mm = MAX(dq, dt); mm *= diff_ec_ul; if(mm < bw) mm = bw;
				if(dd <= mm) {///pass distance checking
                    sc += get_add_cov_score(uref, lj->ts, lj->te, li->ts, li->te, tid, &cov_i);
                }
			}
		} else {
			sc += retrieve_u_cov_region(uref, tid, 0, li->ts, li->te, &cov_i);
		}
	}
	res->score = (sc>0x3FFFFFFF?0x3FFFFFFF:sc);
}

void debug_ll_chains(const ul_idx_t *uref, uint64_t *ix, int64_t ix_n, int64_t p_sidx, int64_t p_eidx, mg_lchain_t *chain_a,
kv_ul_ov_t *raw_idx, kv_ul_ov_t *raw_chn, uint64_t *b, int64_t bw, double diff_ec_ul, int64_t qlen)
{
	int64_t k, i, z, a_n, ss, ee, b_n; uint64_t qs, qe, ts, te; uint32_t mk = 0x80000000; mg_lchain_t nn;
	int64_t iqs, iqe, its, ite, tsc;
	
	for (k = p_sidx; k < p_eidx; k++) {
		i = raw_idx->a[chain_a[k].off].qn; 
		qs = raw_chn->a[i].qs; qe = raw_chn->a[i].qe;
		ts = raw_chn->a[i].ts; te = raw_chn->a[i].te;
		// fprintf(stderr, "++++++++++++[M::%s::idx:%ld]\n", __func__, i);
		for (;i>=0;) {
			// fprintf(stderr, "--[M::%s::i->%ld]\n", __func__, i);
			if(raw_chn->a[i].qs < qs) qs = raw_chn->a[i].qs;
        	if(raw_chn->a[i].ts < ts) ts = raw_chn->a[i].ts;
        	if(raw_chn->a[i].qe > qe) qe = raw_chn->a[i].qe;
        	if(raw_chn->a[i].te > te) te = raw_chn->a[i].te;
			raw_chn->a[i].qs |= mk;

			if(raw_chn->a[i].tn == (uint32_t)-1) i = -1;
			else i = raw_chn->a[i].tn;
		}
		assert(raw_idx->a[chain_a[k].off].qs == qs && raw_idx->a[chain_a[k].off].qe == qe && 
				raw_idx->a[chain_a[k].off].ts == ts && raw_idx->a[chain_a[k].off].te == te);
	}


	///dedup
	for (z = a_n = 0; z < ix_n; ++z) {
		ss = a_n; ee = a_n + ((uint32_t)ix[z]); tsc = 0;
		/**if(ss != p_sidx || ee != p_eidx)**/ {
			for (k = ss; k < ee; k++) {
				i = raw_idx->a[chain_a[k].off].qn; b_n = 0;
				qs = ((raw_chn->a[i].qs<<1)>>1); qe = raw_chn->a[i].qe;
				ts = raw_chn->a[i].ts; te = raw_chn->a[i].te;
				// fprintf(stderr, "++++++++++++[M::%s::idx:%ld]\n", __func__, i);
				for (;i>=0;) {
					// fprintf(stderr, "--[M::%s::i->%ld]\n", __func__, i);
					if(((raw_chn->a[i].qs<<1)>>1) < qs) qs = ((raw_chn->a[i].qs<<1)>>1);
					if(raw_chn->a[i].ts < ts) ts = raw_chn->a[i].ts;
					if(raw_chn->a[i].qe > qe) qe = raw_chn->a[i].qe;
					if(raw_chn->a[i].te > te) te = raw_chn->a[i].te;
					// raw_chn->a[i].qs |= mk;
					//update here!!!!!!!
					// if(!(raw_chn->a[i].qs&mk)) b[b_n++] = i;
					b[b_n++] = i;
					if(raw_chn->a[i].tn == (uint32_t)-1) i = -1;
					else i = raw_chn->a[i].tn;
				}
				assert(raw_idx->a[chain_a[k].off].qs == qs && raw_idx->a[chain_a[k].off].qe == qe && 
						raw_idx->a[chain_a[k].off].ts == ts && raw_idx->a[chain_a[k].off].te == te);
				
				update_exist_chain(uref, raw_chn->a, b, b_n, raw_idx->a[chain_a[k].off].tn, bw, diff_ec_ul, &nn);
				nn.v = (raw_idx->a[chain_a[k].off].tn<<1)|raw_idx->a[chain_a[k].off].rev;
				extend_end_coord(&nn, NULL, qlen, uref->ug->g->seq[raw_idx->a[chain_a[k].off].tn].len, 
				&iqs, &iqe, &its, &ite);
				nn.qs = iqs; nn.qe = iqe; nn.rs = its; nn.re = ite;
				if(!(chain_a[k].score == nn.score && chain_a[k].qs == nn.qs && chain_a[k].qe == nn.qe 
															&& chain_a[k].rs == nn.rs && chain_a[k].re == nn.re)){
					fprintf(stderr, "[M::%s::] chain_a[k].score->%d, nn.score->%d\n", __func__, chain_a[k].score, nn.score);
					fprintf(stderr, "[M::%s::] chain_a[k].qs->%d, nn.qs->%d, chain_a[k].qe->%d, nn.qe->%d, chain_a[k].rs->%d, nn.rs->%d, chain_a[k].re->%d, nn.re->%d\n", __func__, 
					chain_a[k].qs, nn.qs, chain_a[k].qe, nn.qe, chain_a[k].rs, nn.rs, chain_a[k].re, nn.re);
				}
				assert(chain_a[k].score == nn.score && chain_a[k].qs == nn.qs && chain_a[k].qe == nn.qe 
															&& chain_a[k].rs == nn.rs && chain_a[k].re == nn.re);
				tsc += nn.score;
			}
			assert(tsc >= ((int64_t)(ix[z]>>32)));
		}

		a_n += ((uint32_t)ix[z]);
	}



	for (k = p_sidx; k < p_eidx; k++) {
		i = raw_idx->a[chain_a[k].off].qn; 
		qs = ((raw_chn->a[i].qs<<1)>>1); qe = raw_chn->a[i].qe;
		ts = raw_chn->a[i].ts; te = raw_chn->a[i].te;
		for (;i>=0;) {
			// fprintf(stderr, "--[M::%s::i->%ld]\n", __func__, i);
			if(raw_chn->a[i].qs&mk) raw_chn->a[i].qs -= mk;
			if(raw_chn->a[i].qs < qs) qs = raw_chn->a[i].qs;
        	if(raw_chn->a[i].ts < ts) ts = raw_chn->a[i].ts;
        	if(raw_chn->a[i].qe > qe) qe = raw_chn->a[i].qe;
        	if(raw_chn->a[i].te > te) te = raw_chn->a[i].te;
			

			if(raw_chn->a[i].tn == (uint32_t)-1) i = -1;
			else i = raw_chn->a[i].tn;
		}
		assert(raw_idx->a[chain_a[k].off].qs == qs && raw_idx->a[chain_a[k].off].qe == qe && 
				raw_idx->a[chain_a[k].off].ts == ts && raw_idx->a[chain_a[k].off].te == te);
	}
}

void dedup_second_chain(const ul_idx_t *uref, uint64_t *ix, int64_t ix_n, int64_t p_sidx, int64_t p_eidx, mg_lchain_t *chain_a,
kv_ul_ov_t *raw_idx, kv_ul_ov_t *raw_chn, uint64_t *b, int64_t bw, double diff_ec_ul, int64_t qlen)
{
	int64_t k, i, z, a_n, ss, ee, b_n; uint64_t qs, qe, ts, te; uint32_t mk = 0x80000000; mg_lchain_t nn;
	int64_t iqs, iqe, its, ite, tsc;
	
	for (k = p_sidx; k < p_eidx; k++) {
		i = raw_idx->a[chain_a[k].off].qn; 
		qs = raw_chn->a[i].qs; qe = raw_chn->a[i].qe;
		ts = raw_chn->a[i].ts; te = raw_chn->a[i].te;
		// fprintf(stderr, "++++++++++++[M::%s::idx:%ld]\n", __func__, i);
		for (;i>=0;) {
			// fprintf(stderr, "--[M::%s::i->%ld]\n", __func__, i);
			if(raw_chn->a[i].qs < qs) qs = raw_chn->a[i].qs;
        	if(raw_chn->a[i].ts < ts) ts = raw_chn->a[i].ts;
        	if(raw_chn->a[i].qe > qe) qe = raw_chn->a[i].qe;
        	if(raw_chn->a[i].te > te) te = raw_chn->a[i].te;
			raw_chn->a[i].qs |= mk;

			if(raw_chn->a[i].tn == (uint32_t)-1) i = -1;
			else i = raw_chn->a[i].tn;
		}
		assert(raw_idx->a[chain_a[k].off].qs == qs && raw_idx->a[chain_a[k].off].qe == qe && 
				raw_idx->a[chain_a[k].off].ts == ts && raw_idx->a[chain_a[k].off].te == te);
	}


	///dedup
	for (z = a_n = 0; z < ix_n; ++z) {
		ss = a_n; ee = a_n + ((uint32_t)ix[z]); tsc = 0;
		if(ss != p_sidx || ee != p_eidx) {
			for (k = ss; k < ee; k++) {
				i = raw_idx->a[chain_a[k].off].qn; b_n = 0;
				qs = ((raw_chn->a[i].qs<<1)>>1); qe = raw_chn->a[i].qe;
				ts = raw_chn->a[i].ts; te = raw_chn->a[i].te;
				// fprintf(stderr, "++++++++++++[M::%s::idx:%ld]\n", __func__, i);
				for (;i>=0;) {
					// fprintf(stderr, "--[M::%s::i->%ld]\n", __func__, i);
					if(((raw_chn->a[i].qs<<1)>>1) < qs) qs = ((raw_chn->a[i].qs<<1)>>1);
					if(raw_chn->a[i].ts < ts) ts = raw_chn->a[i].ts;
					if(raw_chn->a[i].qe > qe) qe = raw_chn->a[i].qe;
					if(raw_chn->a[i].te > te) te = raw_chn->a[i].te;
					// raw_chn->a[i].qs |= mk;
					//update here!!!!!!!
					if(!(raw_chn->a[i].qs&mk)) b[b_n++] = i;
					if(raw_chn->a[i].tn == (uint32_t)-1) i = -1;
					else i = raw_chn->a[i].tn;
				}
				assert(raw_idx->a[chain_a[k].off].qs == qs && raw_idx->a[chain_a[k].off].qe == qe && 
						raw_idx->a[chain_a[k].off].ts == ts && raw_idx->a[chain_a[k].off].te == te);
				
				update_exist_chain(uref, raw_chn->a, b, b_n, raw_idx->a[chain_a[k].off].tn, bw, diff_ec_ul, &nn);
				nn.v = (raw_idx->a[chain_a[k].off].tn<<1)|raw_idx->a[chain_a[k].off].rev;
				extend_end_coord(&nn, NULL, qlen, uref->ug->g->seq[raw_idx->a[chain_a[k].off].tn].len, &iqs, &iqe, &its, &ite);
				nn.qs = iqs; nn.qe = iqe; nn.rs = its; nn.re = ite;
				assert(chain_a[k].score >= nn.score);
				tsc += (chain_a[k].score - nn.score);
			}
		}
		///TODO: also update qs, qe
		tsc = ((int64_t)(ix[z]>>32)) - tsc; if(tsc < 0) tsc = 0;
		ix[z] <<= 32; ix[z] >>= 32; ix[z] |= ((uint64_t)tsc)<<32;

		a_n += ((uint32_t)ix[z]);
	}



	for (k = p_sidx; k < p_eidx; k++) {
		i = raw_idx->a[chain_a[k].off].qn; 
		qs = ((raw_chn->a[i].qs<<1)>>1); qe = raw_chn->a[i].qe;
		ts = raw_chn->a[i].ts; te = raw_chn->a[i].te;
		for (;i>=0;) {
			// fprintf(stderr, "--[M::%s::i->%ld]\n", __func__, i);
			if(raw_chn->a[i].qs&mk) raw_chn->a[i].qs -= mk;
			if(raw_chn->a[i].qs < qs) qs = raw_chn->a[i].qs;
        	if(raw_chn->a[i].ts < ts) ts = raw_chn->a[i].ts;
        	if(raw_chn->a[i].qe > qe) qe = raw_chn->a[i].qe;
        	if(raw_chn->a[i].te > te) te = raw_chn->a[i].te;
			

			if(raw_chn->a[i].tn == (uint32_t)-1) i = -1;
			else i = raw_chn->a[i].tn;
		}
		assert(raw_idx->a[chain_a[k].off].qs == qs && raw_idx->a[chain_a[k].off].qe == qe && 
				raw_idx->a[chain_a[k].off].ts == ts && raw_idx->a[chain_a[k].off].te == te);
	}
}

uint32_t gen_max_gchain(void *km, const ul_idx_t *uref, int64_t ulid, st_mt_t *idx, vec_mg_lchain_t *e, kv_ul_ov_t *raw_idx, kv_ul_ov_t *raw_chn,
int64_t qlen, float primary_cov_rate, float primary_fragment_cov_rate, float primary_fragment_second_score_rate, const asg_t *g, st_mt_t *dst_done, 
vec_sp_node_t *out, vec_mg_pathv_t *res, uint64_t *b, int64_t bw, double diff_ec_ul)
{
	if(idx->n <= 0) return 0;
	int64_t a_n, idx_n = idx->n, i, m_sc = 0, is_done = 0; uint64_t s_idx, e_idx, om, ok, ovlp, novlp; 
	ul_ov_t m; memset(&m, 0, sizeof(m)); m_sc = -1; mg_lchain_t *a = e->a;
	for (i = a_n = 0; i < idx_n; ++i) {
		if(((int64_t)(idx->a[i]>>32)) > m_sc) {
			m_sc = ((int64_t)(idx->a[i]>>32)); m.qn = i;
			m.ts = a_n; m.te = a_n + ((uint32_t)idx->a[i]);
			m.qs = a[m.ts].qs; m.qe = a[m.te-1].qe;
		}
		// dedup_second_chain(NULL, 0, a_n, a_n + ((uint32_t)idx->a[i]), a, raw_idx, raw_chn);
		a_n += ((uint32_t)idx->a[i]);
	}
	assert(a[m.ts].qs<=a[m.te-1].qs && a[m.te-1].qe>=a[m.ts].qe);
	// print_chain(a + m.ts, m.te - m.ts);
	///for debug
	// dedup_second_chain(uref, idx->a, idx_n, m.ts, m.te, a, raw_idx, raw_chn, b, bw, diff_ec_ul);
	// debug_ll_chains(uref, idx->a, idx_n, m.ts, m.te, a, raw_idx, raw_chn, b, bw, diff_ec_ul, qlen);

	if((m.qe - m.qs) > (qlen*primary_cov_rate)) is_done = 1;
	if(is_done == 0) {
		// for (i = a_n = 0; i < idx_n; ++i) {
		// 	s_idx = a[a_n].qs; a_n += ((uint32_t)idx->a[i]); e_idx = a[a_n-1].qe;
		// 	if(i == m.qn) continue;
		// 	if(s_idx < m.qs || e_idx < m.qs || s_idx > m.qe || e_idx > m.qe) break;
		// }
		// if(i >= idx_n) is_done = 2;///no alignment that is on the left or the right side of the primary chain
		for (i = a_n = 0; i < idx_n; ++i) {
			s_idx = a[a_n].qs; a_n += ((uint32_t)idx->a[i]); e_idx = a[a_n-1].qe;
			if(i == m.qn) continue;
			ovlp = ((MIN(m.qe, e_idx) > MAX(m.qs, s_idx))? (MIN(m.qe, e_idx) - MAX(m.qs, s_idx)):0);
			novlp = (e_idx - s_idx) - ovlp;
			if(novlp > ((m.qe-m.qs)*GC_OFFSET_RATE) && novlp > GC_OFFSET_POS) break;
		}
		if(i >= idx_n) is_done = 2;///no alignment that is on the left or the right side of the primary chain
	}

	if(is_done == 0) {
		if((m.qe - m.qs) > (qlen*primary_fragment_cov_rate)) {
			dedup_second_chain(uref, idx->a, idx_n, m.ts, m.te, a, raw_idx, raw_chn, b, bw, diff_ec_ul, qlen);

			om = m.qe - m.qs;
			for (i = a_n = 0; i < idx_n; ++i) {
				s_idx = a[a_n].qs; a_n += ((uint32_t)idx->a[i]); e_idx = a[a_n-1].qe;
				if(i == m.qn) continue;
				ovlp = ((MIN(m.qe, e_idx) > MAX(m.qs, s_idx))? (MIN(m.qe, e_idx) - MAX(m.qs, s_idx)):0);
				if(ovlp == 0) continue;
				ok = e_idx - s_idx;
				if(ok > om) ok = om;
				if((ovlp > ok*0.1) && ((int64_t)(idx->a[i]>>32)) > (m_sc*primary_fragment_second_score_rate)) break;
			}
			if(i >= idx_n) is_done = 3;
		}
	}

	if(is_done && gen_gchain_track(km, a + m.ts, m.te - m.ts, g, dst_done, out, res)) {///try to find a path
		for (i = m.ts, e->n = 0; i < (int64_t)m.te; i++) a[e->n++] = a[i];
		// fprintf(stderr, "--[M::%s::id->%ld] [%u, %u), res->n:%lu\n", __func__, ulid, m.qs, m.qe, (uint64_t)res->n);
		kv_resize(mg_lchain_t, *e, res->n); a = e->a;
		for (i = ((int64_t)res->n)-1; i >= 0; i--) {
			
			if(res->a[i].v == (uint32_t)-1) {
				a[i] = a[res->a[i].pre]; a[i].dist_pre = res->a[i].d;
				// fprintf(stderr, "ulid:%ld\t%u\t%u\t%c\tutg%.6dl\t%u\t%u\tdist_pre:%d\n", ulid, a[i].qs, a[i].qe, "+-"[a[i].v&1], (int32_t)(a[i].v>>1)+1, a[i].rs, a[i].re, a[i].dist_pre);
			}
			else {
				a[i].v = res->a[i].v; a[i].off = -1; a[i].dist_pre = res->a[i].d;
				// fprintf(stderr, "ulid:%ld\t*\t*\t%c\tutg%.6dl\t*\t*\tdist_pre:%d\n", ulid, "+-"[a[i].v&1], (int32_t)(a[i].v>>1)+1, a[i].dist_pre);
			}
		}
		e->n = res->n;

		


		// debug_gchain(km, g, e->a, e->n, dst_done, out);


		return 1;
	}

	return 0;
}


void update_exist_chain_adv(const ul_idx_t *uref, ul_ov_t *ch, uint64_t *idx, int64_t idx_n, int64_t tid, mg_lchain_t *res)
{
	int64_t i, j, cov_i, sc = 0; memset(res, 0, sizeof(*res));
	ul_ov_t *li, *lj;
	if(idx_n <= 0) return;
	i = 0; 
	res->qs = (ch[idx[i]].qs<<1)>>1; res->qe = ch[idx[i]].qe;
	res->rs = ch[idx[i]].ts; res->re = ch[idx[i]].te;
	for (i = 0; i < idx_n; i++) {
		li = &(ch[idx[i]]); 

		if((int64_t)((li->qs<<1)>>1) < res->qs) res->qs = ((li->qs<<1)>>1); 
		if((int64_t)li->ts < res->rs) res->rs = li->ts;
		if((int64_t)li->qe > res->qe) res->qe = li->qe;
		if((int64_t)li->te > res->re) res->re = li->te;

		cov_i = 0; j = i - 1;
		if(j >= 0) {
			lj = &(ch[idx[j]]); 
			sc += get_add_cov_score(uref, lj->ts, lj->te, li->ts, li->te, tid, &cov_i);
		} else {
			sc += retrieve_u_cov_region(uref, tid, 0, li->ts, li->te, &cov_i);
		}
	}
	res->score = (sc>0x3FFFFFFF?0x3FFFFFFF:sc);
}


void dedup_second_chain_adv(const ul_idx_t *uref, ul_ov_t *gb, int64_t gb_n, mg_lchain_t *chain_a,
kv_ul_ov_t *raw_idx, kv_ul_ov_t *raw_chn, uint64_t *b, int64_t qlen, int64_t ulid)
{
	int64_t k, i, z, ss, ee, b_n, n_s; uint64_t qs, qe, ts, te; uint32_t mk = 0x80000000/**, pi**/; mg_lchain_t nn;
	int64_t iqs, iqe, its, ite, tsc;

	for (z = gb_n - 1; z >= 0; z--) {///start from the best chain
		ss = gb[z].ts; ee = gb[z].te; tsc = 0; gb[z].qs = qlen; gb[z].qe = 0;
		for (k = ss; k < ee; k++) {
			i = raw_idx->a[chain_a[k].off].qn; b_n = 0;
			qs = ((raw_chn->a[i].qs<<1)>>1); qe = raw_chn->a[i].qe;
			ts = raw_chn->a[i].ts; te = raw_chn->a[i].te;
			for (n_s = 0;i>=0;) {
				if(((raw_chn->a[i].qs<<1)>>1) < qs) qs = ((raw_chn->a[i].qs<<1)>>1);
				if(raw_chn->a[i].ts < ts) ts = raw_chn->a[i].ts;
				if(raw_chn->a[i].qe > qe) qe = raw_chn->a[i].qe;
				if(raw_chn->a[i].te > te) te = raw_chn->a[i].te;

				if(!(raw_chn->a[i].qs&mk)) b[b_n++] = i;
				else n_s++;

				raw_chn->a[i].qs |= mk;
				if(raw_chn->a[i].tn == (uint32_t)-1) i = -1;
				else i = raw_chn->a[i].tn;
			}
			assert(b_n > 0);
			assert(raw_idx->a[chain_a[k].off].qs == qs && raw_idx->a[chain_a[k].off].qe == qe && 
					raw_idx->a[chain_a[k].off].ts == ts && raw_idx->a[chain_a[k].off].te == te);
			// if(!(raw_idx->a[chain_a[k].off].qs == qs && raw_idx->a[chain_a[k].off].qe == qe && 
			// 		raw_idx->a[chain_a[k].off].ts == ts && raw_idx->a[chain_a[k].off].te == te)) {
			// 	fprintf(stderr, "\n[M::%s::ulid->%ld******] raw_idx_offset:%d, qs:%lu, qe:%lu, ts:%lu, te:%lu, raw_idx->qs:%u, raw_idx->qe:%u, raw_idx->ts:%u, raw_idx->te:%u\n", 
			// 	__func__, ulid, chain_a[k].off, qs, qe, ts, te, 
			// 	raw_idx->a[chain_a[k].off].qs, raw_idx->a[chain_a[k].off].qe, 
			// 	raw_idx->a[chain_a[k].off].ts, raw_idx->a[chain_a[k].off].te); 

			// 	for (i = raw_idx->a[chain_a[k].off].qn;i>=0;) {
			// 		if(((raw_chn->a[i].qs<<1)>>1) < qs) qs = ((raw_chn->a[i].qs<<1)>>1);
			// 		if(raw_chn->a[i].ts < ts) ts = raw_chn->a[i].ts;
			// 		if(raw_chn->a[i].qe > qe) qe = raw_chn->a[i].qe;
			// 		if(raw_chn->a[i].te > te) te = raw_chn->a[i].te;
			// 		fprintf(stderr, "[M::%s->pieces (%ld)] qs->%u, qe->%u, ts->%u, te->%u\n", __func__, i, 
            // 				((raw_chn->a[i].qs<<1)>>1), raw_chn->a[i].qe, raw_chn->a[i].ts, raw_chn->a[i].te);
			// 		if(raw_chn->a[i].tn == (uint32_t)-1) i = -1;
			// 		else i = raw_chn->a[i].tn;
			// 	}
			// }
			
			update_exist_chain_adv(uref, raw_chn->a, b, b_n, raw_idx->a[chain_a[k].off].tn, &nn);
			nn.v = (raw_idx->a[chain_a[k].off].tn<<1)|raw_idx->a[chain_a[k].off].rev;
			extend_end_coord(&nn, NULL, qlen, uref->ug->g->seq[raw_idx->a[chain_a[k].off].tn].len, &iqs, &iqe, &its, &ite);
			nn.qs = iqs; nn.qe = iqe; nn.rs = its; nn.re = ite;
			assert(chain_a[k].score >= nn.score);
			tsc += (chain_a[k].score - nn.score);

			if(n_s) {
				// raw_idx->a[chain_a[k].off].qn = b[b_n-1];
				// for (i = 0, pi = (uint32_t)-1; i < b_n; i++) {
				// 	raw_chn->a[b[i]].tn = pi; pi = b[i];
				// }
				///don't update chain_a[k] as it will be used for taceback in the next step
				chain_a[k].score = nn.score; 
				chain_a[k].qs = nn.qs; chain_a[k].qe = nn.qe; 
				chain_a[k].rs = nn.rs; chain_a[k].re = nn.re; 
				// raw_idx->a[chain_a[k].off].sec = nn.score; 
				// raw_idx->a[chain_a[k].off].qs = nn.qs; 
				// raw_idx->a[chain_a[k].off].qe = nn.qe; 
				// raw_idx->a[chain_a[k].off].ts = nn.rs; 
				// raw_idx->a[chain_a[k].off].te = nn.re;
			} else {
				// if(!(chain_a[k].score == nn.score && chain_a[k].qs == nn.qs && chain_a[k].qe == nn.qe && chain_a[k].rs == nn.rs && chain_a[k].re == nn.re)) {
				// 	fprintf(stderr, "++++[M::%s::k->%ld] chain_a[k].score->%d, nn.score->%d, chain_a[k].qs->%d, nn.qs->%d, chain_a[k].qe->%d, nn.qe->%d, chain_a[k].rs->%d, nn.rs->%d, chain_a[k].re->%d, nn.re->%d\n", __func__, 
				// 	k, chain_a[k].score, nn.score, chain_a[k].qs, nn.qs, chain_a[k].qe, nn.qe, chain_a[k].rs, nn.rs, chain_a[k].re, nn.re);
				// }
				assert(chain_a[k].score == nn.score && chain_a[k].qs == nn.qs && chain_a[k].qe == nn.qe && chain_a[k].rs == nn.rs && chain_a[k].re == nn.re);
			}

			if((int32_t)gb[z].qs > nn.qs) gb[z].qs = nn.qs;///update qs and qe
			if((int32_t)gb[z].qe < nn.qe) gb[z].qe = nn.qe;
		}
		tsc = (int64_t)(gb[z].qn) - tsc; if(tsc < 0) tsc = 0; gb[z].qn = tsc;///update score
		// if(gb[z].qe <= gb[z].qs) {
		// 	fprintf(stderr, "++++[M::%s::] gb[%ld].qe->%u, gb[%ld].qs->%u\n", __func__, z, gb[z].qe, z, gb[z].qs);
		// }
		assert(gb[z].qe > gb[z].qs);
		// fprintf(stderr, "[M::%s::z->%ld] score->%u, qs->%u, qe->%u, occ->%u\n", 
		// __func__, z, gb[z].qn, gb[z].qs, gb[z].qe, gb[z].te-gb[z].ts); 
	}

	for (i = 0; i < (int64_t)raw_chn->n; i++){
		if(raw_chn->a[i].qs&mk) raw_chn->a[i].qs -= mk;
	}
	radix_sort_ul_ov_srt_qn(gb, gb + gb_n);//sort by scores
}

uint32_t gen_max_gchain_adv(void *km, const ul_idx_t *uref, int64_t ulid, st_mt_t *idx, vec_mg_lchain_t *e, kv_ul_ov_t *raw_idx, kv_ul_ov_t *raw_chn,
int64_t qlen, float primary_cov_rate, float primary_fragment_cov_rate, float primary_fragment_second_score_rate, uint64_t mini_primary_fragment_len, 
const asg_t *g, st_mt_t *dst_done, vec_sp_node_t *out, vec_mg_pathv_t *res, uint64_t *b, vec_mg_lchain_t *gchains)
{
	gchains->n = 0;
	if(idx->n <= 0) return 0;
	int64_t a_n, idx_n = idx->n, i, k, is_done = 0, n_mchain = 0; uint64_t om, ok, ovlp, novlp; 
	ul_ov_t *m = NULL, *p = NULL; mg_lchain_t *a = e->a, *g_item; int64_t raw_idx_n = raw_idx->n;
	ul_ov_t *gb = NULL; int64_t gb_n = 0;
	for (i = a_n = 0; i < idx_n; ++i) {
		kv_pushp(ul_ov_t, *raw_idx, &p);
		p->qn = ((int64_t)(idx->a[i]>>32));//score
		p->ts = a_n; p->te = a_n + ((uint32_t)idx->a[i]);
		p->qs = a[p->ts].qs; p->qe = a[p->te-1].qe; p->tn = 0;//(tn = 1) -> normal; (t = 0) -> duplicated chain
		a_n += ((uint32_t)idx->a[i]);
		// fprintf(stderr, "[M::%s::i->%ld] score->%u, qs->%u, qe->%u, occ->%u\n", __func__, i, p->qn, p->qs, p->qe, p->te-p->ts); 
	}
	gb = raw_idx->a + raw_idx_n; gb_n = raw_idx->n - raw_idx_n;
	radix_sort_ul_ov_srt_qn(gb, gb + gb_n);//sort by scores
	m = &(gb[gb_n-1]);///max chain
	// assert(a[m.ts].qs<=a[m.te-1].qs && a[m.te-1].qe>=a[m.ts].qe);
	// print_chain(a + m.ts, m.te - m.ts);
	///for debug
	// dedup_second_chain(uref, idx->a, idx_n, m.ts, m.te, a, raw_idx, raw_chn, b, bw, diff_ec_ul);
	// debug_ll_chains(uref, idx->a, idx_n, m.ts, m.te, a, raw_idx, raw_chn, b, bw, diff_ec_ul, qlen);

	if((m->qe - m->qs) > (qlen*primary_cov_rate)) {
		is_done = 1; m->tn = 1; n_mchain = 1;
	}

	if(is_done == 0) {
		for (i = gb_n - 2; i >= 0; i--) {///from the second best chain
			p = &(gb[i]);
			ovlp = ((MIN(m->qe, p->qe) > MAX(m->qs, p->qs))? (MIN(m->qe, p->qe) - MAX(m->qs, p->qs)):0);
			novlp = (p->qe - p->qs) - ovlp;
			if(novlp > ((m->qe-m->qs)*GC_OFFSET_RATE) && novlp > GC_OFFSET_POS) break;
		}
		if(i < 0) { ///no alignment that is on the left or the right side of the primary chain
			is_done = 2; m->tn = 1; n_mchain = 1;
		}
	}

	if(is_done == 0) {
		if(((m->qe - m->qs) > (qlen*primary_fragment_cov_rate)) || ((m->qe - m->qs) > mini_primary_fragment_len)) {
			// dedup_second_chain(uref, idx->a, idx_n, m.ts, m.te, a, raw_idx, raw_chn, b, bw, diff_ec_ul, qlen);
			if(raw_chn && raw_idx) dedup_second_chain_adv(uref, gb, gb_n, a, raw_idx, raw_chn, b, qlen, ulid);
			for (k = gb_n-1, n_mchain = 0; k >= 0; k--) {
				m = &(gb[k]);///max chain
				// fprintf(stderr, "++[M::%s::k->%ld] score->%u, qs->%u, qe->%u\n", __func__, k, m->qn, m->qs, m->qe); 
				if(((m->qe - m->qs) <= (qlen*primary_fragment_cov_rate)) && 
															((m->qe - m->qs) <= mini_primary_fragment_len)) break;
				om = m->qe - m->qs;
				for (i = gb_n-1; i >= 0; i--) {
					if(i == k) continue;
					p = &(gb[i]);
					// fprintf(stderr, "--[M::%s::i->%ld] score->%u, qs->%u, qe->%u\n", __func__, i, p->qn, p->qs, p->qe); 
					ovlp = ((MIN(m->qe, p->qe) > MAX(m->qs, p->qs))? (MIN(m->qe, p->qe) - MAX(m->qs, p->qs)):0);
					if(ovlp == 0) continue;
					ok = p->qe - p->qs;
					if(p->tn == 1 && ((ovlp > ok*0.1) || (ovlp > om*0.1))) break;
					if(ok > om) ok = om;
					if((ovlp > ok*0.1) && (p->qn > (m->qn*primary_fragment_second_score_rate))) break;
				}
				if(i < 0) {
					is_done = 3; m->tn = 1; n_mchain++;
				} else {
					 break;
				 }
			}
		}
	}

	if(is_done) {
		gchains->n = 0;
		for (k = gb_n - n_mchain; k < gb_n; k++) {
			kv_pushp(mg_lchain_t, *gchains, &g_item);
			g_item->v = (uint32_t)-1; 
			g_item->qs = gb[k].qs; g_item->qe = gb[k].qe; 
			g_item->rs = gb[k].ts; g_item->re = gb[k].te; 
			g_item->cnt = g_item->off = 0;
			gen_gchain_track(km, a + g_item->rs, g_item->re - g_item->rs, g, dst_done, out, res);
			kv_resize(mg_lchain_t, *gchains, gchains->n + res->n); ///a = gchains->a + gchains->n;
			for (i = 0, g_item = &(gchains->a[gchains->n-1]); i < ((int64_t)res->n); i++) {
				if(res->a[i].v == (uint32_t)-1) {
					gchains->a[i+gchains->n] = a[res->a[i].pre + g_item->rs]; 
					gchains->a[i+gchains->n].dist_pre = res->a[i].d;

					// fprintf(stderr, "+[M::%s::]\tutg%.6dl(%c)\n", __func__, 
					// (int32_t)(gchains->a[i+gchains->n].v>>1)+1, "+-"[gchains->a[i+gchains->n].v&1]);
				} else {
					gchains->a[i+gchains->n].v = res->a[i].v; 
					gchains->a[i+gchains->n].off = -1; 
					gchains->a[i+gchains->n].dist_pre = res->a[i].d;
					///the nodes detected by the graph chaining should be fully covered
					gchains->a[i+gchains->n].rs = 0;
					gchains->a[i+gchains->n].re = uref->ug->g->seq[res->a[i].v>>1].len;
					// fprintf(stderr, "aaaaaaa, ulid->%ld\n", ulid);
					// fprintf(stderr, "-[M::%s::]\tutg%.6dl(%c)\n", __func__, 
					// (int32_t)(gchains->a[i+gchains->n].v>>1)+1, "+-"[gchains->a[i+gchains->n].v&1]);
				}
			}
			g_item->cnt = res->n;
			gchains->n += res->n;
			// fprintf(stderr, "sbsbsbsb, ulid->%ld\n", ulid);
			// debug_gchain(km, g, gchains->a + gchains->n - res->n, res->n, dst_done, out);
		}
	}
	raw_idx->n = raw_idx_n;
	return n_mchain;
}

int64_t extract_rovlp_by_ug(utg_ct_t *p, mg_lchain_t* o, vec_mg_lchain_t *chains, int64_t tOff)
{
	int64_t rs, re;
	rs = MAX((int32_t)p->s, o->rs); re = MIN((int32_t)p->e, o->re);
	if(rs > re) return 0;
	mg_lchain_t *x = NULL;
	kv_pushp(mg_lchain_t, *chains, &x); memset(x, 0, sizeof(*x));
	x->v = (p->x>>1)<<1; x->v += (((o->v&1) == (p->x&1))?0:1); x->rs = rs; x->re = re; x->off = tOff; 
	x->hash_pre = (uint32_t)-1; x->dist_pre = -1; x->qs = x->qe = -1;
	return 1;
}




void update_existing_anchors(ul_vec_t *rch, ma_ug_t *ug, ma_utg_t *u, vec_mg_lchain_t *res, int64_t res_n0, 
mg_lchain_t *uo, kv_ul_ov_t *raw_idx, kv_ul_ov_t *raw_chn)
{
	int64_t z = -1, m = res->n-1, midx, mdif, ovlp, novlp, mbeg, left[2], right[2]; 
	uint64_t zv;
	if(uo->off >= 0) z = raw_idx->a[uo->off].qn;

	for (; z >= 0;) {
		zv = ((rch->bb.a[raw_chn->a[z].qn].hid<<1)+rch->bb.a[raw_chn->a[z].qn].rev);
		for (midx = mdif = -1, mbeg = m; m >= res_n0; m--) {
			if(zv==res->a[m].v) {
				if(((int32_t)raw_chn->a[z].ts == res->a[m].rs && (int32_t)raw_chn->a[z].te == res->a[m].re)) {
					midx = m; mdif = 0; 
					break;
				} else {
					ovlp = ((MIN((int32_t)raw_chn->a[z].te, res->a[m].re) > MAX((int32_t)raw_chn->a[z].ts, res->a[m].rs))? 
								(MIN((int32_t)raw_chn->a[z].te, res->a[m].re) - MAX((int32_t)raw_chn->a[z].ts, res->a[m].rs)):0);
					novlp = (raw_chn->a[z].te - raw_chn->a[z].ts - ovlp) + (res->a[m].re - res->a[m].rs - ovlp);
					if(midx==-1 || mdif>novlp) {
						midx = m; mdif = novlp;
					}
				}
			}
		}
		if(mdif != 0) {
			for (m = res->n-1; m > mbeg; m--) {
				if(zv == res->a[m].v) {
					if(((int32_t)raw_chn->a[z].ts == res->a[m].rs && (int32_t)raw_chn->a[z].te == res->a[m].re)) {
						midx = m; mdif = 0; 
						break;
					} else {
						ovlp = ((MIN((int32_t)raw_chn->a[z].te, res->a[m].re) > MAX((int32_t)raw_chn->a[z].ts, res->a[m].rs))? 
									(MIN((int32_t)raw_chn->a[z].te, res->a[m].re) - MAX((int32_t)raw_chn->a[z].ts, res->a[m].rs)):0);
						novlp = (raw_chn->a[z].te - raw_chn->a[z].ts - ovlp) + (res->a[m].re - res->a[m].rs - ovlp);
						if(midx==-1 || mdif>novlp) {
							midx = m; mdif = novlp;
						}
					}
				}
			}
		}
		m = midx; 
		// if(m < 0) fprintf(stderr, ">>>>>>[M::%s::] z->%ld\n", __func__, z);
		assert(m >= res_n0);
		res->a[m].qs = raw_chn->a[z].qs; res->a[m].qe = raw_chn->a[z].qe; 
		res->a[m].dist_pre = raw_chn->a[z].qn;///the idx of this chain at rch 
		// fprintf(stderr, "%c, res_n0->%ld, m->%ld, raw_idx->%u, qs->%d, qe->%d, ts->%d, te->%d, mdif->%ld\n", 
		// "+-"[(uo->v&1)], res_n0, m, raw_chn->a[z].qn, raw_chn->a[z].qs, raw_chn->a[z].qe, raw_chn->a[z].ts, raw_chn->a[z].te, mdif);
		


		left[0] = left[1] = -1; right[0] = right[1] = u->len+1;
		if(m > 0) get_r_offset(ug, &(res->a[m-1]), &left[0], &left[1], NULL, NULL);
		if(m + 1 < (int64_t)res->n) get_r_offset(ug, &(res->a[m+1]), &right[0], &right[1], NULL, NULL);

		// if(!(uo->v&1)) {///forward
		// 	if(m > 0) {
		// 		left[0] = a[m-1].rs; left[1] = a[m-1].re;
		// 	}
		// 	if(m + 1 < (int64_t)a_n) {
		// 		right[0] = a[m+1].rs; right[1] = a[m+1].re;
		// 	}
		// } else {//reverse
		// 	if(m > 0) {
		// 		right[0] = a[m-1].rs; right[1] = a[m-1].re;
		// 	}
		// 	if(m + 1 < (int64_t)a_n) { 
		// 		left[0] = a[m+1].rs; left[1] = a[m+1].re;
		// 	}
		// }
		///otherwise a[m] is not co-linear with a[m-1] and a[m+1]
		if(raw_chn->a[z].ts>=left[0]&&raw_chn->a[z].ts<=right[0]
						&&raw_chn->a[z].te>=left[1]&&raw_chn->a[z].te<=right[1]) {
            res->a[m].rs = raw_chn->a[z].ts; res->a[m].re = raw_chn->a[z].te; 
		}

		if(raw_chn->a[z].tn == (uint32_t)-1) z = -1;
        else z = raw_chn->a[z].tn;
	}
}

void gl_ug2rg_gen(ul_vec_t *rch, ma_ug_t *ug, mg_lchain_t *uo, vec_mg_lchain_t *res, int64_t tOff, 
kv_ul_ov_t *raw_idx, kv_ul_ov_t *raw_chn)
{
	// fprintf(stderr, "\n[M::%s::] uo->qs:%d, uo->qe:%d\n", __func__, uo->qs, uo->qe);
	///uo is a unitig alignment
	ma_utg_t *u = &(ug->u.a[uo->v>>1]); int64_t res_n0 = res->n;
	uint64_t rs = uo->rs, re = uo->re, i, l; utg_ct_t p; 

	for (i = l = 0; i < u->n; i++) {
		p.x = u->a[i]>>32; p.s = l; p.e = l + Get_READ_LENGTH(R_INF, (u->a[i]>>33));
		l += (uint32_t)u->a[i];
		if(p.e <= rs) continue; 
		if(p.s >= re) break;
		// fprintf(stderr, "[M::%s::]rs->%lu, re->%lu, p.s->%u, p.e->%u\n", __func__, rs, re, p.s, p.e); 
		assert(extract_rovlp_by_ug(&p, uo, res, tOff));
		res->a[res->n-1].score = uo->v; res->a[res->n-1].cnt = i;
		///for res->a[res->n-1]
		///ts and te are the coordinates in unitig (res->a[res->n-1].score>>1), instead of HiFi read (res->a[res->n-1].v>>1)
		///qs and qe are the coordinates in UL, 
	}

	mg_lchain_t *a = res->a + res_n0, t; uint64_t a_n = res->n - res_n0; 
	if(uo->v&1) {
		for (i = 0; i < (a_n>>1); i++) {
			t = a[i]; a[i] = a[a_n-i-1]; a[a_n-i-1] = t;
		}
	}

	update_existing_anchors(rch, ug, u, res, res_n0, uo, raw_idx, raw_chn);
}


void update_rovlp_chain_qse(ma_ug_t *ug, int64_t sidx, int64_t eidx, mg_lchain_t *a, int64_t a_n)
{
	if(eidx - sidx <= 1) return;
	assert(sidx>=0||eidx<a_n);
	// fprintf(stderr, "******[M::%s::] sidx:%ld, eidx:%ld\n", __func__, sidx, eidx);
	int64_t left_r[2], right_r[2], left_q[2], right_q[2]; 
	int64_t i, rs, re;//qs or qe might be -1, while rs and re should >= 0
	if(sidx >= 0) {
		get_r_offset(ug, &(a[sidx]), &left_r[0], &left_r[1], &left_q[0], &left_q[1]);
	} else {
		get_r_offset(ug, &(a[0]), &left_r[0], &left_r[1], &left_q[0], &left_q[1]);
	}

	if(eidx < a_n) {
		get_r_offset(ug, &(a[eidx]), &right_r[0], &right_r[1], &right_q[0], &right_q[1]);
	} else {
		get_r_offset(ug, &(a[a_n-1]), &right_r[0], &right_r[1], &right_q[0], &right_q[1]);
	}
	assert((left_q[0] >= 0 && left_q[1] >= 0) || (right_q[0] >= 0 && right_q[1] >= 0)); ///assert(re >= rs);
	// fprintf(stderr, "##[M::%s::] right_q[0]:%ld, right_q[1]:%ld, left_q[0]:%ld, left_q[1]:%ld\n", 
	// __func__, right_q[0], right_q[1], left_q[0], left_q[1]);
	// fprintf(stderr, "##[M::%s::] right_r[0]:%ld, right_r[1]:%ld, left_r[0]:%ld, left_r[1]:%ld\n", 
	// __func__, right_r[0], right_r[1], left_r[0], left_r[1]);
	// if(left_q[0] >= 0 && left_q[1] >= 0 && right_q[0] >= 0 && right_q[1] >= 0) {
	// 	rlen[0] = (right_r[0] - left_r[0]); rlen[1] = (right_r[1] - left_r[1]); 
	// 	qlen[0] = (right_q[0] - left_q[0]); qlen[1] = (right_q[1] - left_q[1]); 
	// 	for (i = sidx+1; i < eidx; i++) {
	// 		get_r_offset(ug, &(a[i]), &rs, &re, NULL, NULL);
	// 		a[i].qs = left_q[0] + get_offset_adjust((rs - left_r[0]), rlen[0], qlen[0]);
	// 		a[i].qe = left_q[1] + get_offset_adjust((re - left_r[1]), rlen[1], qlen[1]);
	// 	}
	// }

	if(left_q[0] >= 0 && left_q[1] >= 0 && right_q[0] >= 0 && right_q[1] >= 0) {
		// fprintf(stderr, "+++sidx:%ld+++ left_qs:%ld, left_qe:%ld, left_rs:%ld, left_re:%ld\n", 
		// sidx, left_q[0], left_q[1], left_r[0], left_r[1]);
		// fprintf(stderr, "---eidx:%ld--- right_qs:%ld, right_qe:%ld, right_rs:%ld, right_re:%ld\n", 
		// eidx, right_q[0], right_q[1], right_r[0], right_r[1]);

		for (i = sidx+1; i < eidx; i++) {
			get_r_offset(ug, &(a[i]), &rs, &re, NULL, NULL);
			// a[i].qs = left_q[0] + get_offset_adjust((rs - left_r[0]), rlen[0], qlen[0]);
			///a[i].qs>=left_q[0] && a[i].qs<left_q[0]
			// a[i].qs = left_q[0] + get_offset_adjust(rs - left_r[0], left_r[1]-left_r[0], left_q[1]-left_q[0]);
			a[i].qs = left_q[0] + get_offset_adjust(rs - left_r[0], right_r[0]-left_r[0], right_q[0]-left_q[0]);
			// a[i].qe = left_q[1] + get_offset_adjust((re - left_r[1]), rlen[1], qlen[1]);
			a[i].qe = left_q[1] + get_offset_adjust((re - left_r[1]), right_r[1] - left_r[1], right_q[1] - left_q[1]);
			
			// fprintf(stderr, ">>>i:%ld<<< a[i].qs:%u, a[i].qe:%u, rs:%ld, re:%ld\n", i, a[i].qs, a[i].qe, rs, re);

			left_q[0] = a[i].qs; left_q[1] = a[i].qe;
			left_r[0] = rs; left_r[1] = re;
		}
	}

	if(right_q[0] < 0 || right_q[1] < 0) {
		for (i = sidx+1; i < eidx; i++) {
			get_r_offset(ug, &(a[i]), &rs, &re, NULL, NULL);
			///a[i].qs>=left_q[0] && a[i].qs<left_q[0]
			a[i].qs = left_q[0] + get_offset_adjust(rs - left_r[0], left_r[1]-left_r[0], left_q[1]-left_q[0]);
			///a[i].qe>=left_q[1]
			a[i].qe = left_q[1] + (re - left_r[1]);
			left_q[0] = a[i].qs; left_q[1] = a[i].qe;
			left_r[0] = rs; left_r[1] = re;
		}
	}

	if(left_q[0] < 0 || left_q[1] < 0) {
		for (i = eidx-1; i > sidx; i--) {
			get_r_offset(ug, &(a[i]), &rs, &re, NULL, NULL);
			a[i].qe = right_q[1] - get_offset_adjust(right_r[1]-re, right_r[1]-right_r[0], right_q[1]-right_q[0]);
			a[i].qs = right_q[0] - (right_r[0]-rs);
			right_q[0] = a[i].qs; right_q[1] = a[i].qe;
			right_r[0] = rs; right_r[1] = re;
		}
	}
	// if(left_q[0] < 0) left_q[0] = right_q[0] - (right_r[0] - left_r[0]);
	// if(left_q[1] < 0) left_q[1] = right_q[1] - (right_r[1] - left_r[1]);
	// if(right_q[0] < 0 || right_q[1] < 0) {
	// 	right_q[0] = left_q[0] + (right_r[0] - left_r[0]);
	// 	right_q[1] = left_q[1] + (right_r[1] - left_r[1]);
	// }	
	
	// fprintf(stderr, "******[M::%s::] right_q[0]:%ld, right_q[1]:%ld\n", __func__, right_q[0], right_q[1]);
}

void gen_rovlp_chain_by_ul(ul_vec_t *rch, const ul_idx_t *uref, kv_ul_ov_t *raw_idx, kv_ul_ov_t *raw_chn, mg_lchain_t *a, int64_t a_n, vec_mg_lchain_t *res)
{
	if(a_n == 0) return;
	int64_t k, l, res_n0 = res->n, tt = 0; ma_ug_t *ug = uref->ug;
	// fprintf(stderr, "\n[M::%s::a_n->%ld]\n", __func__, a_n);
	///a[0, a_n) is a gchain of untigs
	for (k = 0, l = ug->g->seq[a[0].v>>1].len; k < a_n; k++) {
		l -= ug->g->seq[a[k].v>>1].len;
		gl_ug2rg_gen(rch, ug, &(a[k]), res, l, raw_idx, raw_chn);
		l += ug->g->seq[a[k].v>>1].len + a[k].dist_pre;
	}
	mg_lchain_t *x = res->a + res_n0; int64_t x_n = res->n - res_n0;
	for (l = -1, k = 0; k <= x_n; k++) {
		if(k < x_n) {
			if(k > 0) x[k].hash_pre = k-1+res_n0;
			else x[k].hash_pre = (uint32_t)-1;
			if(x[k].qs >= 0) tt++;
		}
		if(k == x_n || x[k].qs >=0) { ///x[k] and x[l] are anchors
			if(k-l>1) update_rovlp_chain_qse(ug, l, k, x, x_n);
			l = k;
		}
	}
	assert(tt > 0);
}



int64_t convert_mg_lchain_t(utg_ct_t *p, mg_lchain_t *o)
{
	int64_t rs = p->s, re = p->e;
	rs = MAX(rs, o->rs); re = MIN(re, o->re);
	assert(rs < re);	
	if(!(p->x&1)) {
		o->rs = rs-p->s; o->re = re-p->s;
	} else {
		o->rs = p->e-re; o->re = p->e-rs;
	}
	return 1;
}

void renew_mg_lchains(ma_ug_t *ug, mg_lchain_t *a, int64_t a_n)
{
	if (a_n <= 0) return;
	uint32_t rev = (a[0].score&1); ma_utg_t *u = &(ug->u.a[a[0].score>>1]);
	uint64_t i, l; int64_t k; utg_ct_t p; 

	if(!rev) {
		for (i = l = 0; i < u->n; i++) {
			if(i == (uint64_t)a[0].cnt) break;
			l += (uint32_t)u->a[i];
		}
		assert(i < u->n);
		for (k = 0; k < a_n; k++, i++) {
			assert(a[k].cnt == (int64_t)i && (a[k].v>>1) == (u->a[i]>>33)); 
			p.x = u->a[i]>>32; p.s = l; p.e = l + Get_READ_LENGTH(R_INF, (u->a[i]>>33));
			convert_mg_lchain_t(&p, &a[k]);
			l += (uint32_t)u->a[i];
		}		
	} else {
		for (i = l = 0; i < u->n; i++) {
			if(i == (uint64_t)a[a_n-1].cnt) break;
			l += (uint32_t)u->a[i];
		}
		assert(i < u->n);
		for (k = a_n-1; k >= 0; k--, i++) {
			assert(a[k].cnt == (int64_t)i && (a[k].v>>1) == (u->a[i]>>33)); 
			p.x = u->a[i]>>32; p.s = l; p.e = l + Get_READ_LENGTH(R_INF, (u->a[i]>>33));
			convert_mg_lchain_t(&p, &a[k]);
			l += (uint32_t)u->a[i];
		}
	}
}


int64_t g_adjacent_dis_mul(const asg_t *g, ma_hit_t_alloc *src, int64_t max_hang, int64_t min_ovlp, uint32_t v, uint32_t w)
{
	uint32_t i;
	if(g) {
		uint32_t nv; asg_arc_t *av = NULL;
		nv = asg_arc_n(g, v); av = asg_arc_a(g, v);
		for (i = 0; i < nv; i++) {
			if(av[i].del || av[i].v != w) continue;
			return (uint32_t)av[i].ul;
		}
	}

	if(src) {
		ma_hit_t_alloc *x = &(src[v>>1]); uint32_t qn, tn; 
		int32_t r; asg_arc_t e;
		for (i = 0; i < x->length; i++) {
			qn = Get_qn(x->buffer[i]);
			tn = Get_tn(x->buffer[i]);
			if(qn == (v>>1) && tn == (w>>1)) {
				r = ma_hit2arc(&(x->buffer[i]), Get_READ_LENGTH(R_INF, qn), Get_READ_LENGTH(R_INF, tn),
					max_hang, asm_opt.max_hang_rate, min_ovlp, &e);
				if(r < 0) continue;
				if((e.ul>>32) != v || e.v != w) continue;
				return (uint32_t)e.ul;
			}
		}
	}
	
	return -1;
}


void dd_ul_vec_t(const ul_idx_t *uref, mg_lchain_t *a, int64_t a_n, ul_vec_t *rch)
{
	int64_t k, l, ovlp, novlp, tt, rs, re; uint64_t i; uc_block_t *z; mg_lchain_t *p, *c;
	for (l = 0, k = 1; k <= a_n; k++) {
		if(k == a_n || a[k].score != a[l].score) { ///x[k] and x[l] come from the same unitig
			renew_mg_lchains(uref->ug, a + l, k - l);
			l = k;
		}
	}

	
	for (i = 0; i < rch->bb.n; i++) {
		rch->bb.a[i].pidx = 0xfffffffe;
		rch->bb.a[i].aidx = rch->bb.a[i].pdis = (uint32_t)-1;
	}

	for (k = 0; k < a_n; k++) {
		if(a[k].dist_pre >= 0) {///not a new alignment
			z = &(rch->bb.a[a[k].dist_pre]); 
			assert(z->hid == (a[k].v>>1) && z->rev == (a[k].v&1));
			if((int64_t)z->qs == a[k].qs && (int64_t)z->qe == a[k].qe && (int64_t)z->ts == a[k].rs && (int64_t)z->te == a[k].re) {
				z->pidx = k; z->pchain = 1;
			} else {
				ovlp = novlp = tt = 0;
				ovlp = ((MIN((int64_t)z->qe, a[k].qe) > MAX((int64_t)z->qs, a[k].qs))? (MIN((int64_t)z->qe, a[k].qe)-MAX((int64_t)z->qs, a[k].qs)):0);
				tt += ovlp;
				novlp += (a[k].qe - a[k].qs - ovlp) + (z->qe - z->qs - ovlp);

				ovlp = ((MIN((int64_t)z->te, a[k].re) > MAX((int64_t)z->ts, a[k].rs))? (MIN((int64_t)z->te, a[k].re)-MAX((int64_t)z->ts, a[k].rs)):0);
				tt += ovlp;
				novlp += (a[k].re - a[k].rs - ovlp) + (z->te - z->ts - ovlp);
				if(novlp > 8 || novlp > (tt*0.01)) {
					kv_pushp(uc_block_t, rch->bb, &z); 
					z->hid = (a[k].v>>1); z->rev = (a[k].v&1); 
					z->pchain = 2; z->base = 0; z->el = 1;
					z->qs = a[k].qs; z->qe = a[k].qe;
					z->te = a[k].re; z->ts = a[k].rs;
					z->pidx = k;
				} else {
					z->qs = a[k].qs; z->qe = a[k].qe;
					z->te = a[k].re; z->ts = a[k].rs;
					z->pidx = k; z->pchain = 1;
				}
			}
		} else {
			kv_pushp(uc_block_t, rch->bb, &z); 
			z->hid = (a[k].v>>1); z->rev = (a[k].v&1); 
			z->pchain = 2; z->base = 0; z->el = 1;
			z->qs = a[k].qs; z->qe = a[k].qe;
			z->te = a[k].re; z->ts = a[k].rs;
			z->pidx = k;
		}
	}

	for (i = k = 0; i < rch->bb.n; i++) {
		if(rch->bb.a[i].pidx == 0xfffffffe && (rch->bb.a[i].pchain != 1 || rch->bb.a[i].pchain != 0)) continue;
		rch->bb.a[k++] = rch->bb.a[i];
	}
	rch->bb.n = k;

	radix_sort_uc_block_t_qe_srt(rch->bb.a, rch->bb.a + rch->bb.n);
	for (i = 0; i < rch->bb.n; i++) {
		if(rch->bb.a[i].pidx == 0xfffffffe) {
			rch->bb.a[i].pidx = (uint32_t)-1;
		} else {
			a[rch->bb.a[i].pidx].dist_pre = i;
		}
	}
	// fprintf(stderr, "\n[M::%s::]\n", __func__);
	for (i = 0, k = -1; i < rch->bb.n; i++) {
		if(rch->bb.a[i].pidx == (uint32_t)-1) continue;
		if(k < 0) k = i;///in case there is only one UL-to-HiFi alignment
		if(a[rch->bb.a[i].pidx].hash_pre == (uint32_t)-1) {
			rch->bb.a[i].pidx = (uint32_t)-1;
			continue;
		}
		c = &(a[rch->bb.a[i].pidx]); p = &(a[a[rch->bb.a[i].pidx].hash_pre]);
		rch->bb.a[i].pidx = a[a[rch->bb.a[i].pidx].hash_pre].dist_pre; 
		rch->bb.a[rch->bb.a[i].pidx].aidx = i;
		tt = g_adjacent_dis_mul(uref->r_ug->rg, NULL, -1, -1, 
			((rch->bb.a[i].hid<<1)|((uint32_t)rch->bb.a[i].rev))^1, 
			((rch->bb.a[rch->bb.a[i].pidx].hid<<1)|((uint32_t)rch->bb.a[rch->bb.a[i].pidx].rev))^1);
		if(tt >= 0) {
			rch->bb.a[i].pdis = tt;
			// get_r_offset(uref->ug, p, NULL, &rs, NULL, NULL);
			// get_r_offset(uref->ug, c, NULL, &re, NULL, NULL);
			// fprintf(stderr, "+i->%lu: dis->%u, record_dis->%ld\n", i, rch->bb.a[i].pdis, re-rs);
		} else {
			rs = p->off + uref->ug->g->seq[p->score>>1].len;
			re = c->off + uref->ug->g->seq[c->score>>1].len;
			if(re >= rs) rch->bb.a[i].pdis = re - rs;
			else rch->bb.a[i].pdis = (uint32_t)-1;
			// fprintf(stderr, "-i->%lu: dis->%u\n", i, rch->bb.a[i].pdis);
		}

		k = i;
		// if(rch->bb.a[i].base || rch->bb.a[i].pchain == 0 || rch->bb.a[i].el == 0) {
		// 	fprintf(stderr, "+++(%lu) base:%u, pchain:%u, el:%u\n", 
		// 	i, rch->bb.a[i].base, rch->bb.a[i].pchain, rch->bb.a[i].el);
		// }
	}

	///make sure if this UL read has been done
	uint32_t sp = (uint32_t)-1, ep = (uint32_t)-1;
	for (l = 0 ; k >= 0; ) {
		// fprintf(stderr, "k->%ld\n", k);
		if(sp == (uint32_t)-1 || rch->bb.a[k].qe <= sp) {
            if(sp != (uint32_t)-1) l += ep - sp;
            sp = rch->bb.a[k].qs; ep = rch->bb.a[k].qe;
        } else {
            sp = MIN(sp, rch->bb.a[k].qs);
        }
		if(rch->bb.a[k].pidx == (uint32_t)-1) k = -1;
		else k = rch->bb.a[k].pidx;
	}

	if(sp != (uint32_t)-1) l += ep - sp;
	if(l == (int64_t)rch->rlen) rch->dd = 1;

	// for (i = 0; i < rch->bb.n; i++) {
	// 	if(rch->bb.a[i].pidx == (uint32_t)-1) continue;
	// 	if(rch->bb.a[i].base || rch->bb.a[i].pchain == 0 || rch->bb.a[i].el == 0) {
	// 		fprintf(stderr, "(%lu) base:%u, pchain:%u, el:%u\n", 
	// 		i, rch->bb.a[i].base, rch->bb.a[i].pchain, rch->bb.a[i].el);
	// 	}
	// 	assert((!(rch->bb.a[i].base)) && (rch->bb.a[i].pchain) && (rch->bb.a[i].el));
	// 	assert((!(rch->bb.a[rch->bb.a[i].pidx].base)) && (rch->bb.a[rch->bb.a[i].pidx].pchain) 
	// 														&& (rch->bb.a[rch->bb.a[i].pidx].el));
	// }
	// int64_t exact = 0, inexact = 0; uc_block_t *z;
	// for (k = 0; k < a_n; k++) {
	// 	// fprintf(stderr, "(%ld) a->qs:%d, a->qe:%d, a->rs:%d, a->re:%d\n", k, a[k].qs, a[k].qe, a[k].rs, a[k].re);
	// 	if(a[k].dist_pre < 0) continue;
	// 	z = &(rch->bb.a[a[k].dist_pre]);
	// 	assert(z->hid == (a[k].v>>1) && z->rev == (a[k].v&1));
	// 	if((int64_t)z->qs == a[k].qs && (int64_t)z->qe == a[k].qe && (int64_t)z->ts == a[k].rs && (int64_t)z->te == a[k].re) {
	// 		exact++;
	// 	} else {
	// 		inexact++;
	// 		// fprintf(stderr, "+z->qs:%u, z->qe:%u, z->ts:%u, z->te:%u\n", z->qs, z->qe, z->ts, z->te);
	// 		// fprintf(stderr, "-a->qs:%d, a->qe:%d, a->rs:%d, a->re:%d\n\n", a[k].qs, a[k].qe, a[k].rs, a[k].re);
	// 	}
	// }
	// fprintf(stderr, "[M::%s::exact->%ld, inexact->%ld]\n", __func__, exact, inexact);
}

void update_ul_vec_t(const ul_idx_t *uref, kv_ul_ov_t *raw_idx, kv_ul_ov_t *raw_chn, ul_vec_t *rch, 
vec_mg_lchain_t *uc, vec_mg_lchain_t *swap)
{
	int64_t k, ucn = uc->n; mg_lchain_t *ix; 
	for (k = 0, swap->n = 0; k < ucn; k += ix->cnt + 1) {
		ix = &(uc->a[k]); assert(ix->v == (uint32_t)-1);
		// fprintf(stderr, "\n[M::%s::ucn->%ld, k->%ld, kcnt->%d]\n", __func__, ucn, k, ix->cnt);
		gen_rovlp_chain_by_ul(rch, uref, raw_idx, raw_chn, uc->a + k + 1, ix->cnt, swap);
	}

	///up to now, given a <x> in swap
	///x->ts and x->te are the coordinates in unitig (x->score>>1), instead of HiFi read (x->v>>1)
	///x->qs and x->qe are the coordinates in UL
	///x->dist_pre is the idx of this chain at rch 
	// debug_intermediate_chain(uref->ug, swap->a, swap->n);
	dd_ul_vec_t(uref, swap->a, swap->n, rch);
}

void print_ru_raw_chains(kv_ul_ov_t *raw_idx, kv_ul_ov_t *raw_chn, vec_mg_lchain_t *gch, ul_vec_t *rch, ma_ug_t *ug)
{
	int64_t k, i; uint64_t ts, te, qs, qe;
	for (k = 0; k < (int64_t)gch->n; k++) {
		i = raw_idx->a[gch->a[k].off].qn;
		qs = ((raw_chn->a[i].qs<<1)>>1); qe = raw_chn->a[i].qe;
		fprintf(stderr, "\n[M::%s->overall chain (%ld)] utg%.6d%c(%c), qs->%u, qe->%u, qlen->%u, ts->%u, te->%u, tlen->%u\n", __func__, k, 
		(int32_t)(gch->a[k].v>>1)+1, "lc"[ug->u.a[(gch->a[k].v>>1)].circ], "+-"[(gch->a[k].v&1)],
		raw_idx->a[gch->a[k].off].qs, raw_idx->a[gch->a[k].off].qe, rch->rlen,
		raw_idx->a[gch->a[k].off].ts, raw_idx->a[gch->a[k].off].te, ug->u.a[(gch->a[k].v>>1)].len);
		ts = raw_chn->a[i].ts; te = raw_chn->a[i].te;
		for (;i>=0;) {
			if(((raw_chn->a[i].qs<<1)>>1) < qs) qs = ((raw_chn->a[i].qs<<1)>>1);
			if(raw_chn->a[i].ts < ts) ts = raw_chn->a[i].ts;
			if(raw_chn->a[i].qe > qe) qe = raw_chn->a[i].qe;
			if(raw_chn->a[i].te > te) te = raw_chn->a[i].te;

			fprintf(stderr, "[M::%s->chain pieces (%ld)] qs->%u, qe->%u, ts->%u, te->%u\n", __func__, i, 
			((raw_chn->a[i].qs<<1)>>1), raw_chn->a[i].qe, raw_chn->a[i].ts, raw_chn->a[i].te);

			if(raw_chn->a[i].tn == (uint32_t)-1) i = -1;
			else i = raw_chn->a[i].tn;
		}

	}
	
}

///sps and hap are just vector for uint64_t; used for buffer
uint32_t direct_gchain(mg_tbuf_t *b, ul_vec_t *rch, glchain_t *ll, gdpchain_t *gdp, st_mt_t *sps, haplotype_evdience_alloc *hap, const ul_idx_t *uref, const ug_opt_t *uopt,
int64_t bw, double diff_ec_ul, int64_t max_skip, int64_t ulid)
{
	// if(ulid != 7768/** && ulid != 44522**/) return 0;
	kv_ul_ov_t *idx = &(ll->lo), *init = &(ll->tk); int64_t max_idx;
	idx->n = init->n = 0;
	gl_rg2ug_gen(rch, idx, uref, 1, 2);
	if(idx->n == 0) return 0;
	///generate linear chains
	gen_linear_chains(idx, init, uref, uopt, bw, diff_ec_ul, rch->rlen, max_skip, ll, sps);
	assert(idx->n);
	if(idx->n == 0) return 0;
	// fprintf(stderr, "\n++[M::%s::%.*s(id:%ld), len:%u] idx->n:%lu\n", __func__, UL_INF.nid.a[ulid].n, UL_INF.nid.a[ulid].a,
	// ulid, rch->rlen, (uint64_t)idx->n);

	dump_linear_chain(uref->ug->g, idx, init, &(gdp->l), rch->rlen);
	// fprintf(stderr, "\n+++[M::%s::id->%ld, len->%u] idx->n:%lu\n", __func__, ulid, rch->rlen, (uint64_t)idx->n);
	// kv_resize(uint64_t, ll->srt.a, idx->n); kv_resize(uint64_t, hap->snp_srt, idx->n); kv_resize(uint64_t, gdp->v, idx->n);
	// occ = gl_chain_advance(&(gdp->l), &(gdp->swap), uref, uopt, G_CHAIN_BW, diff_ec_ul, qlen, UG_SKIP, dumy->overlapID, ll->srt.a.a, hap->snp_srt.a, G_CHAIN_TRANS_WEIGHT, 0, NULL, uref->ug, debug_i, km);
	// print_ru_raw_chains(idx, init,  &(gdp->l), rch, uref->ug);

	///buffer
	kv_resize(uint64_t, ll->srt.a, gdp->l.n); kv_resize(uint64_t, hap->snp_srt, gdp->l.n);
	kv_resize(uint64_t, gdp->v, gdp->l.n); kv_resize(int64_t, gdp->f, gdp->l.n);
	max_idx = hc_gchain1_dp(b->km, uref, uref->ug, &(gdp->l), &(gdp->swap), &(gdp->dst), &(gdp->out), &(gdp->path), rch->rlen,
	uopt, bw, diff_ec_ul, N_GCHAIN_RATE, ll->srt.a.a, sps, gdp->f.a, hap->snp_srt.a, gdp->v.a);
	// fprintf(stderr, "++++[M::%s::id->%ld, len->%u] gdp->l.n:%lu\n", __func__, ulid, rch->rlen, (uint64_t)gdp->l.n);
	// fprintf(stderr, "+[M::%s::] gdp->l.n:%lu\n", __func__, (uint64_t)gdp->l.n);
	//sps has the chain idx; gdp->l has the chain
	if(max_idx >= 0 && gen_max_gchain_adv(b->km, uref, ulid, sps, &(gdp->l), idx, init, rch->rlen, P_CHAIN_COV, 0.3/**P_FRAGEMENT_PRIMARY_CHAIN_COV**/, 
			0.1/**P_FRAGEMENT_PRIMARY_SECOND_COV**/, PRIMARY_UL_CHAIN_MIN, uref->ug->g, &(gdp->dst_done), &(gdp->out), &(gdp->path), ll->srt.a.a, &(gdp->swap))) {
	// if(max_idx >= 0 && gen_max_gchain(b->km, uref, ulid, sps, &(gdp->l), idx, init, rch->rlen, P_CHAIN_COV, 0.3/**P_FRAGEMENT_PRIMARY_CHAIN_COV**/, 
	// 		0.1/**P_FRAGEMENT_PRIMARY_SECOND_COV**/, uref->ug->g, &(gdp->dst_done), &(gdp->out), &(gdp->path), ll->srt.a.a, bw, diff_ec_ul)) {
		// update_ul_vec_t(rch, &(gdp->l), uref);
		// fprintf(stderr, "\n++[M::%s::(id:%ld), len:%u]\n", __func__, ulid, rch->rlen);
		update_ul_vec_t(uref, idx, init, rch, &(gdp->swap), &(gdp->l));
		// __ac_X31_hash_string("hehe");

		return (rch->dd == 1?1:0);
	// } else {
	// 	// uint64_t i;
	// 	fprintf(stderr, "unsuccess->[M::%s::id->%ld, len->%u] gdp->l.n:%lu\n", __func__, ulid, rch->rlen, (uint64_t)gdp->l.n);
	// 	for (i = 0; i < gdp->l.n; ++i) {
	// 		fprintf(stderr, "(%lu)\t%u\t%u\t%c\tutg%.6d%c(%u)\t%u\t%u\tsrc:%u\tscore:%d\n", 
	// 			i, gdp->l.a[i].qs, gdp->l.a[i].qe, "+-"[gdp->l.a[i].v&1], (int32_t)(gdp->l.a[i].v>>1)+1, "lc"[uref->ug->u.a[gdp->l.a[i].v>>1].circ], uref->ug->u.a[gdp->l.a[i].v>>1].len,
	// 			gdp->l.a[i].rs, gdp->l.a[i].re, gdp->l.a[i].v^1, gdp->l.a[i].score);
	// 	}
	}

	

	// occ = gl_chain_advance(idx, ll->tk.a+ll->tk.n, uref, uopt, G_CHAIN_BW, diff_ec_ul, qlen, UG_SKIP, dumy->overlapID, ll->srt.a.a, hap->snp_srt.a, G_CHAIN_TRANS_WEIGHT, 0, NULL, uref->ug, debug_i, km);

	// simple_g_chain_dp(idx, buf->a, uref, uopt, bw, diff_ec_ul, rch->rlen, max_skip, ll->srt.a.a, hap->snp_srt.a, sps->a);
	// if(check_extension_end(idx, rch->rlen, buf->a)) {
	// 	// ug2rg_gen(idx->a[idx->n-1].qs, idx->a[idx->n-1].qe, buf->a + idx->a[idx->n-1].ts, idx->a[idx->n-1].te - idx->a[idx->n-1].ts, rch);
	// } else {///need graph chaining

	// }

	return 0;
}


static void worker_for_ul_gchains_alignment(void *data, long i, int tid)
{
	ul_vec_t *p = &(UL_INF.a[i]);
	if(p->dd == 1) return; //fully aligned
	if(p->bb.n == 1 && p->bb.a[0].base) return;///no alignment
	if(p->bb.n == 0) return;///no alignment
	utepdat_t *s = (utepdat_t*)data;
	s->hab[tid]->num_read_base++;
	s->hab[tid]->num_correct_base += direct_gchain(s->buf[tid], p, &(s->ll[tid]), &(s->gdp[tid]), &(s->sps[tid]), &(s->hab[tid]->hap), s->uu, s->uopt, G_CHAIN_BW, s->opt->diff_ec_ul, UG_SKIP, i);
	// gl_chain_refine_advance(&b->olist, &b->correct, &b->hap, bl, s->uu, s->opt->diff_ec_ul, winLen, s->len[i], s->uopt, s->id+i, km);
}

uint64_t work_ul_gchains(uldat_t *sl)
{
	utepdat_t s; uint64_t i; memset(&s, 0, sizeof(s));
	s.id = 0; s.opt = sl->opt; s.ug = sl->ug; s.uopt = sl->uopt; s.rg = sl->rg; s.uu = sl->uu; 
	CALLOC(s.hab, sl->n_thread); CALLOC(s.buf, sl->n_thread); CALLOC(s.ll, sl->n_thread); 
	CALLOC(s.gdp, sl->n_thread); CALLOC(s.mzs, sl->n_thread); CALLOC(s.sps, sl->n_thread); 

	for (i = 0; i < sl->n_thread; ++i) {
		s.hab[i] = ha_ovec_init(0, 0, 1); s.buf[i] = mg_tbuf_init();
	}

	kt_for(sl->n_thread, worker_for_ul_gchains_alignment, &s, UL_INF.n);

	for (i = 0; i < sl->n_thread; ++i) {
		s.sum_len += s.hab[i]->num_read_base; s.n += s.hab[i]->num_correct_base;
		ha_ovec_destroy(s.hab[i]); mg_tbuf_destroy(s.buf[i]); hc_glchain_destroy(&(s.ll[i]));
		hc_gdpchain_destroy(&(s.gdp[i])); kv_destroy(s.mzs[i]); kv_destroy(s.sps[i]);
	}

	free(s.hab); free(s.buf); free(s.ll); free(s.gdp); free(s.mzs); free(s.sps);
	fprintf(stderr, "[M::%s::] # try:%d, # done:%d\n", __func__, s.sum_len, s.n); 
	return s.n;
}

void print_ul_ovlps(all_ul_t *x, int32_t prt_ovlp)
{
	uint64_t k, i, ucov_occ = 0, cov_occ = 0, ucov_len = 0, cov_len = 0, unaligned_len = 0, unaligned_occ = 0, aligned_occ = 0;
	ul_vec_t *p = NULL; nid_t *z = NULL; uc_block_t *m = NULL;
	for (k = 0; k < x->n; k++) {
		z = &(x->nid.a[k]);
		p = &(x->a[k]);
		fprintf(stderr, "S\t%.*s\tq:id:%lu\tl:%u\tdd:%d\n", (int32_t)z->n, z->a, k, p->rlen, 
		((p->bb.n == 1&&p->bb.a[0].base)||(p->bb.n==0))?-1:(int32_t)p->dd);
		if(prt_ovlp) {
			for (i = 0; i < p->bb.n; i++) {
				m = &(p->bb.a[i]);
				if(m->base) {
					ucov_occ++;
					ucov_len += (m->qe-(m->hid&FLANK_M)) - (m->qs+((m->hid>>15)&FLANK_M));
					fprintf(stderr, "B\t%.*s\t%u\t%u\t%u\n", 
					(int32_t)z->n, z->a, p->rlen, (m->qs+((m->hid>>15)&FLANK_M)), (m->qe-(m->hid&FLANK_M)));
				} else {
					fprintf(stderr, "A\t%.*s\t%u\t%u\t%u\t%c\t%.*s\t%u\t%u\t%u\n", 
					(int32_t)z->n, z->a, p->rlen, m->qs, m->qe, "+-"[m->rev], 
					(int32_t)Get_NAME_LENGTH(R_INF, m->hid), Get_NAME(R_INF, m->hid), 
					(uint32_t)Get_READ_LENGTH(R_INF, m->hid), m->ts, m->te);
					if(m->el) cov_occ++;
				}
			}
		}
		if((p->bb.n == 1 && p->bb.a[0].base)||(p->bb.n == 0)) {
			unaligned_len += p->rlen; unaligned_occ++;
		} else {
			aligned_occ++;
		}
		cov_len += p->rlen;
	}
	cov_len -= ucov_len;
	fprintf(stderr, "[M::%s::] ==>aligned_occ:%lu, unaligned_occ:%lu\n", __func__, aligned_occ, unaligned_occ);
	fprintf(stderr, "[M::%s::] ==>cov_len:%lu, ucov_len:%lu, unaligned_len:%lu\n", 
																	__func__, cov_len, ucov_len-unaligned_len, unaligned_len);
}

void print_all_ul_t_stat(all_ul_t *x)
{
	uint64_t k, i, ucov_occ = 0, cov_occ = 0, ucov_len = 0, cov_len = 0;
	ul_vec_t *p = NULL;
	for (k = 0; k < x->n; k++) {
		p = &(x->a[k]);
		for (i = 0; i < p->bb.n; i++) {
			if(p->bb.a[i].base/**.hid&x->mm**/) {
				ucov_occ++;
				ucov_len += (p->bb.a[i].qe-(p->bb.a[i].hid&FLANK_M)) - 
					(p->bb.a[i].qs+((p->bb.a[i].hid>>15)&FLANK_M));
			} else {
				cov_occ++;
			}
		}
		cov_len += p->rlen;
	}
	cov_len -= ucov_len;
	fprintf(stderr, "[M::%s::] ==>cov_occ:%lu, ucov_occ:%lu\n", __func__, cov_occ, ucov_occ);
	fprintf(stderr, "[M::%s::] ==>cov_len:%lu, ucov_len:%lu\n", __func__, cov_len, ucov_len);
}

void print_ovlp_src_bl_stat(all_ul_t *x, const ug_opt_t *uopt)
{
	uint64_t k, z, tc, ta;
	ma_hit_t_alloc* src = uopt->sources;
	for (k = tc = ta = 0; k < R_INF.total_reads; k++) {
		if(x->ridx.idx.a[k+1] - x->ridx.idx.a[k] == 0) continue;
		tc++;
		for (z = 0; z < src[k].length; z++) {
			if(src[k].buffer[z].bl) {
				ta++;
				break;
			}
		}
	}
	
	fprintf(stderr, "[M::%s::] ==> # HiFi reads:%lu, # covered HiFi reads:%lu, # chained HiFi reads:%lu\n", 
	__func__, R_INF.total_reads, tc, ta);

	uint64_t tt[4] = {0};
	for (k = 0; k < x->n; k++) tt[x->a[k].dd]++;

	fprintf(stderr, "[M::%s::] ==> # passed UL reads:%lu, # fully corrected UL reads:%lu, # almost fully corrected UL reads:%lu, # UL reads have primary chains:%lu\n", 
	__func__, tt[0]+tt[1]+tt[2]+tt[3], tt[1], tt[2], tt[3]);
}

void gen_ul_vec_rid_t(all_ul_t *x)
{
	ul_vec_rid_t *ridx = &(x->ridx);
	uint64_t k, i, l, m, *a, a_n; ul_vec_t *p = NULL;
	ridx->idx.n = ridx->idx.m = R_INF.total_reads + 1; CALLOC(ridx->idx.a, ridx->idx.n);

	for (k = 0; k < x->n; k++) {
		p = &(x->a[k]);
		for (i = 0; i < p->bb.n; i++) {
			if(p->bb.a[i].base) continue;
			ridx->idx.a[p->bb.a[i].hid]++;
		}
	}

	for (k = l = 0; k < ridx->idx.n; k++) {
		m = ridx->idx.a[k];
		ridx->idx.a[k] = l;
		l += m;
	}

	ridx->occ.n = ridx->occ.m = l; MALLOC(ridx->occ.a, ridx->occ.n);
	for (k = 0; k < R_INF.total_reads; k++) {
		a = ridx->occ.a + ridx->idx.a[k];
		a_n = ridx->idx.a[k+1] - ridx->idx.a[k];
		if(a_n) a[a_n-1] = 0;
	}

	for (k = 0; k < x->n; k++) {
		p = &(x->a[k]);
		for (i = 0; i < p->bb.n; i++) {
			if(p->bb.a[i].base) continue;
			a = ridx->occ.a + ridx->idx.a[p->bb.a[i].hid];
			a_n = ridx->idx.a[p->bb.a[i].hid+1] - ridx->idx.a[p->bb.a[i].hid];
			if(a_n) {
				if(a[a_n-1] == a_n-1) a[a_n-1] = (k<<32)|i;
				else a[a[a_n-1]++] = (k<<32)|i;
			}
		}
	}
	
}

int32_t find_ul_block_max_reverse(int32_t n, const uc_block_t *a, uint32_t x)
{
    int32_t s = 0, e = n;
    if (n == 0) return n;
    if (a[0].qe < x) return 0;///max qe
    if (a[n-1].qe >= x) return n;///min qe

    while (e > s) { // TODO: finish this block
        int32_t m = s + (e - s) / 2;
        // if (a[m].qe >= x) e = m;
        // else s = m + 1;
		if (a[m].qe > x) s = m + 1;
		else e = m;
    }
    assert(s == e);
    return s;
}

int32_t find_ul_block_max(int32_t n, const uc_block_t *a, uint32_t x)
{
    int32_t s = 0, e = n;
    if (n == 0) return -1;
    if (a[n-1].qe < x) return n - 1;
    if (a[0].qe >= x) return -1;
    while (e > s) { // TODO: finish this block
        int32_t m = s + (e - s) / 2;
        if (a[m].qe >= x) e = m;
        else s = m + 1;
    }
    assert(s == e);
    return s;
}
/**
void determine_connective(all_ul_t *m, const ug_opt_t *uopt, int64_t bw, double diff_ec_ul, ul_vec_t *p, uint32_t ii, uint64_t rid)
{
	if((p->bb.a[ii].base) || p->bb.a[ii].hid != rid) fprintf(stderr, "ERROR\n");
	if(p->bb.n <= ii + 1) return;
	uint32_t li_v, lk_v, k, ol; int64_t mm_ovlp, x; 
	uc_block_t *li = NULL, *lk = NULL;
	ma_hit_t *t = NULL;
	li = &(p->bb.a[ii]); li_v = (((uint32_t)(li->hid))<<1)|((uint32_t)(li->rev));
	mm_ovlp = max_ovlp_src(uopt, li_v^1);
	x = (li->qs + mm_ovlp)*diff_ec_ul; 
	if(x < bw) x = bw; 
	x += li->qs + mm_ovlp;
	if (x > p->rlen+1) x = p->rlen+1;
	x = find_ul_block_max_rev(p->bb.n - ii - 1, p->bb.a + ii + 1, x) + ii + 1;
	for (k = x; k < p->bb.n; ++k) { // collect potential destination vertices
		lk = &(p->bb.a[k]); lk_v = (((uint32_t)(lk->hid))<<1)|((uint32_t)(lk->rev));
		if(lk->qe <= li->qs) break;//even this pair has a overlap, its length will be very small; just ignore
		if((li_v == lk_v) || (lk->base)) continue;
		// if(li->qs <= 0) continue;///means the UL read does not longer than the overlap between li and lk
		// if(lk->qs <= 0) continue;//the UL read should be cover the whole HiFi reads li and lk
		if(((li->te - li->ts)*1.05) < Get_READ_LENGTH(R_INF, li->hid)) continue;
		if(((lk->te - lk->ts)*1.05) < Get_READ_LENGTH(R_INF, lk->hid)) continue;
		x = infer_rovlp(NULL, NULL, li, lk, &R_INF, NULL);
		t = query_ovlp_src(uopt, li_v^1, lk_v^1, x, diff_ec_ul, &ol);
		if(t) {
			// sum = t->bl + ol;
			// t->bl = (sum & 0x7fffffffU);
			t->bl++;
		}
	}
}
**/

///note: we only label reliable chains
void determine_connective_adv(all_ul_t *m, const ug_opt_t *uopt, int64_t bw, double diff_ec_ul, ul_vec_t *p, uint32_t ii, uint64_t rid)
{
	assert((!p->bb.a[ii].base)&&(p->bb.a[ii].hid == rid));
	if(ii <= 0) return;
	if(!(p->bb.a[ii].pchain)) return; ///not a primary chain
	if(!(p->bb.a[ii].el)) return; ///not a cis alignment
	uint32_t li_v, lk_v, ol; int64_t mm_ovlp, k, x; 
	uc_block_t *li = NULL, *lk = NULL;
	ma_hit_t *t = NULL;
	li = &(p->bb.a[ii]); li_v = (((uint32_t)(li->hid))<<1)|((uint32_t)(li->rev));
	mm_ovlp = max_ovlp_src(uopt, li_v^1);
	x = (li->qs + mm_ovlp)*diff_ec_ul; 
	if(x < bw) x = bw; 
	x += li->qs + mm_ovlp;
	if (x > p->rlen+1) x = p->rlen+1;
	x = find_ul_block_max(ii, p->bb.a, x+G_CHAIN_INDEL);
	for (k = x; k >= 0; --k) { // collect potential destination vertices
		lk = &(p->bb.a[k]); lk_v = (((uint32_t)(lk->hid))<<1)|((uint32_t)(lk->rev));
		if(lk->qe+G_CHAIN_INDEL <= li->qs) break;//even this pair has a overlap, its length will be very small; just ignore
		if(lk->base || (!(lk->pchain)) || (!(lk->el))) continue; 
		if(li_v == lk_v) continue;
		// if(li->qs <= 0) continue;///means the UL read does not longer than the overlap between li and lk
		// if(lk->qs <= 0) continue;//the UL read should be cover the whole HiFi reads li and lk
		// if(((li->te - li->ts)*1.05) < Get_READ_LENGTH(R_INF, li->hid)) continue;
		// if(((lk->te - lk->ts)*1.05) < Get_READ_LENGTH(R_INF, lk->hid)) continue;
		if((li->te - li->ts) < Get_READ_LENGTH(R_INF, li->hid)) continue;
		if((lk->te - lk->ts) < Get_READ_LENGTH(R_INF, lk->hid)) continue;
		x = /**((int64_t)(lk->qe))-((int64_t)(li->qs))**/infer_rovlp(NULL, NULL, li, lk, &R_INF, NULL);
		t = query_ovlp_src(uopt, li_v^1, lk_v^1, x, diff_ec_ul, &ol);
		if(t) {
			// sum = t->bl + ol;
			// t->bl = (sum & 0x7fffffffU);
			t->bl++;
		}
	}
}

void determine_connective_backtrack(all_ul_t *m, const ug_opt_t *uopt, ul_vec_t *p, uint32_t ii, uint64_t rid)
{
	assert((!p->bb.a[ii].base)&&(p->bb.a[ii].hid == rid)&&(p->bb.a[ii].el));
	if(ii <= 0) return;
	if(!(p->bb.a[ii].pchain)) return; ///not a primary chain
	if(p->bb.a[ii].pidx == (uint32_t)-1) return; ///not connected
	uint32_t li_v, lk_v, z, qn, tn; int32_t r; uc_block_t *li = NULL, *lk = NULL; asg_arc_t t; 
	li = &(p->bb.a[ii]); li_v = (((uint32_t)(li->hid))<<1)|((uint32_t)(li->rev)); li_v^=1;
	if((li->te - li->ts) < Get_READ_LENGTH(R_INF, li->hid)) return;
	ma_hit_t_alloc *x = &(uopt->sources[li_v>>1]);
    int64_t min_ovlp = uopt->min_ovlp;
    int64_t max_hang = uopt->max_hang;
	// if(rid == 4217) fprintf(stderr, "rid: %lu, p->bb.n: %u\n", rid, p->bb.n);

	for (lk = &(p->bb.a[li->pidx]); lk; ) {
		lk_v = (((uint32_t)(lk->hid))<<1)|((uint32_t)(lk->rev)); lk_v^=1;
		assert((!(lk->base)) && (lk->pchain) && (lk->el));
		if((lk->te - lk->ts) >= Get_READ_LENGTH(R_INF, lk->hid)) {
			for (z = 0; z < x->length; z++) {
				qn = Get_qn(x->buffer[z]);
				tn = Get_tn(x->buffer[z]);
				if(qn == (li_v>>1) && tn == (lk_v>>1)) {
					r = ma_hit2arc(&(x->buffer[z]), Get_READ_LENGTH(R_INF, li_v>>1), Get_READ_LENGTH(R_INF, lk_v>>1),
						max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
					if(r < 0) continue;
					if((t.ul>>32) != li_v || t.v != lk_v) continue;
					break;
				}
			}
			if(z < x->length) x->buffer[z].bl++;
		}

		lk = ((lk->pidx==(uint32_t)-1)?NULL:&(p->bb.a[lk->pidx]));
	}
}

static void update_ovlp_src(void *data, long i, int tid) // callback for kt_for()
{
	uldat_t *sl = (uldat_t *)data;
	ma_hit_t_alloc* src = sl->uopt->sources;
	uint64_t z, k, *a, a_n;
	for (z = 0; z < src[i].length; z++) src[i].buffer[z].bl = 0;
	a = UL_INF.ridx.occ.a + UL_INF.ridx.idx.a[i];
	a_n = UL_INF.ridx.idx.a[i+1] - UL_INF.ridx.idx.a[i];
	// fprintf(stderr, "++i->%d, a_n->%lu++\n", i, a_n);
	for (k = 0; k < a_n; k++) {
		///note: we only label reliable chains
		// if(i == 4217) fprintf(stderr, "\nul_id: %lu\n", a[k]>>32);
		determine_connective_backtrack(&UL_INF, sl->uopt, &(UL_INF.a[a[k]>>32]), (uint32_t)(a[k]), i);
		// determine_connective_adv(&UL_INF, sl->uopt, G_CHAIN_BW, sl->opt->diff_ec_ul, &(UL_INF.a[a[k]>>32]), (uint32_t)(a[k]), i);
		// determine_connective(&UL_INF, sl->uopt, G_CHAIN_BW, sl->opt->diff_ec_ul, 
		// 										&(UL_INF.a[a[k]>>32]), (uint32_t)(a[k]), i);
	}
	// fprintf(stderr, "--i->%d, a_n->%lu--\n", i, a_n);
}

uint64_t* get_hifi2ul_list(all_ul_t *x, uint64_t hid, uint64_t* a_n)
{
	(*a_n) = x->ridx.idx.a[hid+1] - x->ridx.idx.a[hid];
	return x->ridx.occ.a + x->ridx.idx.a[hid];
}


static void update_ovlp_src_bl(void *data, long i, int tid)
{
	uldat_t *sl = (uldat_t *)data;
	ma_hit_t_alloc* src = sl->uopt->sources;
	uint64_t z, sum; uint32_t qn, tn; int32_t idx;
	for (z = 0; z < src[i].length; z++) {
		qn = Get_qn(src[i].buffer[z]);
		tn = Get_tn(src[i].buffer[z]);
		if(qn > tn) continue;
		idx = get_specific_overlap(&(src[tn]), tn, qn);
		assert(idx != -1);
		sum = src[i].buffer[z].bl + src[tn].buffer[idx].bl;
		src[i].buffer[z].bl = src[tn].buffer[idx].bl = sum/**(sum&0x7fffffffU)**/;
	}
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
	fprintf(stderr, "[M::%s::] ==> # reads: %lu, # bases: %lu\n", __func__, UL_INF.n, sl->total_base);
	fprintf(stderr, "[M::%s::] ==> # bases: %lu; # corrected bases: %lu; # recorrected bases: %lu\n", 
	__func__, sl->num_bases, sl->num_corrected_bases, sl->num_recorrected_bases);
	gen_ul_vec_rid_t(&UL_INF);
    return 1;
}


int rescall_ul_pipeline(uldat_t* sl, const enzyme *fn)
{
    double index_time = yak_realtime();
    int32_t i; uint32_t k;
	///clean UL_INF
	for (k = 0; k < UL_INF.n; k++) {
		free(UL_INF.a[k].bb.a); free(UL_INF.a[k].N_site.a); free(UL_INF.a[k].r_base.a); 
		memset(&(UL_INF.a[k]), 0, sizeof(UL_INF.a[k]));
	}

    for (i = 0; i < fn->n; i++){
        gzFile fp;
        if ((fp = gzopen(fn->a[i], "r")) == 0) return 0;
        sl->ks = kseq_init(fp);
        kt_pipeline(3, worker_ul_rescall_pipeline, sl, 3);
        kseq_destroy(sl->ks);
        gzclose(fp);
    }
	sl->hits.total_base = sl->total_base;
	sl->hits.total_pair = sl->total_pair;
    fprintf(stderr, "[M::%s::%.3f] ==> Qualification\n", __func__, yak_realtime()-index_time);
	fprintf(stderr, "[M::%s::] ==> # reads: %lu, # bases: %lu, # fully corrected reads: %lu\n", 
	__func__, UL_INF.n, sl->total_base, sl->num_corrected_bases);
	// fprintf(stderr, "[M::%s::] ==> # bases: %lu; # corrected bases: %lu; # recorrected bases: %lu\n", 
	// __func__, sl->num_bases, sl->num_corrected_bases, sl->num_recorrected_bases);
	// gen_ul_vec_rid_t(&UL_INF);
    return 1;
}

int print_ul_rs(all_ul_t *U_INF)
{
	uint32_t i;
	UC_Read ur;
	init_UC_Read(&ur);
	for (i = 0; i < U_INF->n; i++) {
		retrieve_ul_t(&ur, NULL, U_INF, i, 0, 0, -1);
		fprintf(stderr, ">%s\n", U_INF->nid.a[i].a);
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
    init_mg_opt(&opt, 0, 19, 10, hap_n, 0, 0, 0.05, asm_opt.ul_error_rate_low, asm_opt.ul_ec_round);
	int exist = (asm_opt.load_index_from_disk? uidx_load(&ha_flt_tab, &ha_idx, asm_opt.output_file_name, NULL) : 0);
    if(exist == 0) uidx_build(ug, &opt);
	if(exist == 0) uidx_write(ha_flt_tab, ha_idx, asm_opt.output_file_name, NULL);
    ul_align(&opt, uopt, rg, asm_opt.ar, ha_flt_tab, ha_idx, ug);
    uidx_destory();
}

void ul_v_call(uldat_t *sl, const enzyme *fn)
{
	scall_ul_pipeline(sl, fn);
	// UL_INF;
	// print_ul_rs(&UL_INF);
	// debug_retrieve_rc_sub(uopt, &UL_INF, &R_INF, (ul_idx_t *)sl.uu, 100);
	// if(!load_ul_hits(&sl.hits, &sl.nn, asm_opt.output_file_name)) {
    // 	scall_ul_pipeline(&sl, fn);
	// 	write_ul_hits(&sl.hits, &sl.nn, asm_opt.output_file_name);
	// }
}

void ul_v_recall(uldat_t *sl, const enzyme *fn)
{
	rescall_ul_pipeline(sl, fn);
	// UL_INF;
	// print_ul_rs(&UL_INF);
	// debug_retrieve_rc_sub(uopt, &UL_INF, &R_INF, (ul_idx_t *)sl.uu, 100);
	// if(!load_ul_hits(&sl.hits, &sl.nn, asm_opt.output_file_name)) {
    // 	scall_ul_pipeline(&sl, fn);
	// 	write_ul_hits(&sl.hits, &sl.nn, asm_opt.output_file_name);
	// }
}

void print_dedup_HiFis_seq(ma_ug_t *ug)
{
	uint64_t i;
	ma_utg_t *p = NULL;
	for (i = 0; i < ug->u.n; i++) {
		p = &(ug->u.a[i]);
		CALLOC(p->s, p->len+1);
		retrieve_u_seq(NULL, p->s, p, 0, 0, -1, NULL);
		p->s[p->len] = '\0';
	}

	FILE* output_file = fopen("dedup_HiFis_seq.gfa", "w");
    ma_ug_print(ug, NULL, NULL, NULL, NULL, "utg", output_file);
    fclose(output_file);

	output_file = fopen("dedup_HiFis_seq.noseq.gfa", "w");
	ma_ug_print_simple(ug, NULL, NULL, NULL, NULL, "utg", output_file);
	fclose(output_file);
	exit(1);
}

void push_coverage_track(ucov_t *cc, ul_contain *ct, uint64_t uid, ma_utg_t *u, asg_t *rg, ma_hit_t_alloc* src, int64_t min_ovlp, int64_t max_hang, uint64_t is_el, uint64_t is_del)
{
	uint64_t k, l, z, dp, ct_n; utg_ct_t *ct_a = NULL; 
	cc->idx[uid] = cc->interval.n;
	ct_n = ((uint32_t)(ct->idx.a[uid])); ct_a = ct->rids.a + ((ct->idx.a[uid])>>32);
	for (z = 0; z < ct_n; z++) {
		kv_push(uint64_t, cc->interval, ct_a[z].s<<1);
		kv_push(uint64_t, cc->interval, (ct_a[z].e<<1)|1);
	}
	for (k = l = 0; k < u->n; k++) {
		kv_push(uint64_t, cc->interval, l<<1);
		kv_push(uint64_t, cc->interval, ((l + Get_READ_LENGTH(R_INF, u->a[k]>>33))<<1)|1);
		/**
		i = u->a[k]>>33;///rid
		for (z = 0; z < src[i].length; z++) {
            if(is_el && (!src[i].buffer[z].el)) continue;
			if(is_del && (!src[i].buffer[z].del)) continue;
            qn = Get_qn(src[i].buffer[z]); tn = Get_tn(src[i].buffer[z]);
            if(!rg->seq[tn].del) continue;
            if((Get_qe(src[i].buffer[z]) - Get_qs(src[i].buffer[z])) < min_ovlp) continue;
            if((Get_te(src[i].buffer[z]) - Get_ts(src[i].buffer[z])) < min_ovlp) continue;
            r = ma_hit2arc(&(src[i].buffer[z]), rg->seq[qn].len, rg->seq[tn].len, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
			if(r != MA_HT_TCONT) continue;///tn is contained
			ori = (u->a[k]>>32)&1;
			if(ori == 0) {
				qs = Get_qs(src[i].buffer[z]); qe = Get_qe(src[i].buffer[z]);
			} else {
				qs = (Get_READ_LENGTH(R_INF, i)) - Get_qe(src[i].buffer[z]);
				qe = (Get_READ_LENGTH(R_INF, i)) - Get_qs(src[i].buffer[z]);
			}
			kv_push(uint64_t, cc->interval, (l+qs)<<1);
			kv_push(uint64_t, cc->interval, ((l+qe)<<1)|1);
        }
		**/
        l += (uint32_t)u->a[k];
    }
	cc->idx[uid+1] = cc->interval.n;


	radix_sort_gfa64(cc->interval.a+cc->idx[uid], cc->interval.a+cc->interval.n);
	for (k = cc->idx[uid], dp = 0; k < cc->interval.n; ++k) {
        ///if a[j] is qe
        if (cc->interval.a[k]&1) --dp;
        else ++dp;
		l = cc->interval.a[k]>>1; l <<= 32; l += dp;
		cc->interval.a[k] = l;
    }
}


uint32_t check_if_fully_contain(uint32_t sid, uint32_t lid, uint32_t ori, uint8_t *rset, asg_t *rg, 
ma_hit_t_alloc* src, int64_t min_ovlp, int64_t max_hang, int64_t gap_fuzz)
{
	uint32_t rid, k, qn, tn, ff = 1; int32_t r; asg_arc_t t;
	return 1;

	rid = lid;
	for (k = 0; k < src[rid].length; k++) {
		if(!src[rid].buffer[k].el) continue;
		qn = Get_qn(src[rid].buffer[k]); tn = Get_tn(src[rid].buffer[k]);
		if(rg->seq[qn].del || rg->seq[tn].del) continue;
		if((Get_qe(src[rid].buffer[k]) - Get_qs(src[rid].buffer[k])) < min_ovlp) continue;
		if((Get_te(src[rid].buffer[k]) - Get_ts(src[rid].buffer[k])) < min_ovlp) continue;
		r = ma_hit2arc(&(src[rid].buffer[k]), rg->seq[qn].len, rg->seq[tn].len, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
		if(r < 0) continue;
		rset[t.v] = ((t.ul>>32)&1)+1;
	}



	rid = sid;
	for (k = 0; k < src[rid].length; k++) {
		if(!src[rid].buffer[k].el) continue;
		qn = Get_qn(src[rid].buffer[k]); tn = Get_tn(src[rid].buffer[k]);
		if(rg->seq[qn].del || rg->seq[tn].del) continue;
		if((Get_qe(src[rid].buffer[k]) - Get_qs(src[rid].buffer[k])) < min_ovlp) continue;
		if((Get_te(src[rid].buffer[k]) - Get_ts(src[rid].buffer[k])) < min_ovlp) continue;
		r = ma_hit2arc(&(src[rid].buffer[k]), rg->seq[qn].len, rg->seq[tn].len, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
		if(r < 0) continue;
		if(rset[t.v] != ((((t.ul>>32)&1)^ori)+1)) {
			ff = 0; 
			break;
		}
	}


	rid = lid;
	for (k = 0; k < src[rid].length; k++) {
		if(!src[rid].buffer[k].el) continue;
		qn = Get_qn(src[rid].buffer[k]); tn = Get_tn(src[rid].buffer[k]);
		if(rg->seq[qn].del || rg->seq[tn].del) continue;
		if((Get_qe(src[rid].buffer[k]) - Get_qs(src[rid].buffer[k])) < min_ovlp) continue;
		if((Get_te(src[rid].buffer[k]) - Get_ts(src[rid].buffer[k])) < min_ovlp) continue;
		r = ma_hit2arc(&(src[rid].buffer[k]), rg->seq[qn].len, rg->seq[tn].len, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
		if(r < 0) continue;
		rset[t.v] = 0;
	}

	return ff;
}

ul_contain *ul_contain_gen(ma_ug_t *ug, asg_t *rg, ma_hit_t_alloc* src, int64_t min_ovlp, int64_t max_hang, uint64_t is_el, uint64_t is_del)
{
	uint64_t k, l, i, z, t, qn, tn, ori, qs, qe, ovlp, o_z, o_r, o_o;
	ul_contain *p = NULL; ma_utg_t *u = NULL;
	utg_ct_t *m = NULL;
	int32_t r; asg_arc_t e;
	CALLOC(p, 1);
	p->idx.n = p->idx.m = ug->u.n; CALLOC(p->idx.a, p->idx.n);
	p->is_c.n = rg->n_seq; CALLOC(p->is_c.a, p->is_c.n);
	for (t = 0; t < ug->u.n; t++) {
		u = &(ug->u.a[t]);
		p->idx.a[t] = p->rids.n; p->idx.a[t] <<= 32;

		for (k = l = 0; k < u->n; k++) {
			i = u->a[k]>>33;///rid

			for (z = 0; z < src[i].length; z++) {
				if(is_el && (!src[i].buffer[z].el)) continue;
				if(is_del && (!src[i].buffer[z].del)) continue;
				qn = Get_qn(src[i].buffer[z]); tn = Get_tn(src[i].buffer[z]);
				if(!rg->seq[tn].del) continue;
				if((Get_qe(src[i].buffer[z]) - Get_qs(src[i].buffer[z])) < min_ovlp) continue;
				if((Get_te(src[i].buffer[z]) - Get_ts(src[i].buffer[z])) < min_ovlp) continue;
				r = ma_hit2arc(&(src[i].buffer[z]), rg->seq[qn].len, rg->seq[tn].len, max_hang, asm_opt.max_hang_rate, min_ovlp, &e);
            	if(r != MA_HT_TCONT) continue;///tn is contained
				p->is_c.a[qn] = 1;
				ori = (u->a[k]>>32)&1;
				if(ori == 0) {
                	qs = Get_qs(src[i].buffer[z]); qe = Get_qe(src[i].buffer[z]);
				} else {
					qs = (Get_READ_LENGTH(R_INF, i)) - Get_qe(src[i].buffer[z]);
					qe = (Get_READ_LENGTH(R_INF, i)) - Get_qs(src[i].buffer[z]);
				}
				qs += l; qe += l;
				kv_pushp(utg_ct_t, p->rids, &m);
				m->x = tn; m->x <<= 1; m->x |= (src[i].buffer[z].rev == ori?0:1);
				m->s = qs; m->e = qe/** + 1**/;
			}

			l += (uint32_t)u->a[k];
		}

		radix_sort_utg_ct_t_x_srt(p->rids.a + (p->idx.a[t]>>32), p->rids.a + p->rids.n);
		/**
		for (k = (p->idx.a[t]>>32) + 1, l = i = (p->idx.a[t]>>32); k <= p->rids.n; ++k) { 
			if (k == p->rids.n || (p->rids.a[k].x>>1) != (p->rids.a[l].x>>1)) {
				p->rids.a[i] = p->rids.a[l]; i++;
				l = k;
			}
		}
		**/
		for (k = (p->idx.a[t]>>32) + 1, l = i = (p->idx.a[t]>>32); k <= p->rids.n; ++k) { 
			if (k == p->rids.n || p->rids.a[k].x != p->rids.a[l].x) {
				for (z = l; z < k; z++) {
					for (r = (int64_t)i-1; r >= 0 && p->rids.a[r].x == p->rids.a[z].x; r--) {
						ovlp = ((MIN(p->rids.a[z].e, p->rids.a[r].e) > MAX(p->rids.a[z].s, p->rids.a[r].s))?
									(MIN(p->rids.a[z].e, p->rids.a[r].e) - MAX(p->rids.a[z].s, p->rids.a[r].s)):0);
						if(ovlp) {
							o_z = p->rids.a[z].e - p->rids.a[z].s; 
							o_r = p->rids.a[r].e - p->rids.a[r].s;
							o_o = MIN(o_z, o_r);
							if((ovlp <= o_o*1.05) && (ovlp >= o_o*0.95)) break;
						}
						// if(p->rids.a[z].s == p->rids.a[r].s && p->rids.a[z].e == p->rids.a[r].e) break;
					}
					if(r >= 0 && p->rids.a[r].x == p->rids.a[z].x) continue;
					p->rids.a[i++] = p->rids.a[z];
				}
				l = k;
			}
		}
		p->rids.n = i;

		radix_sort_utg_ct_t_s_srt(p->rids.a + (p->idx.a[t]>>32), p->rids.a + p->rids.n);
		p->idx.a[t] |= (p->rids.n - (p->idx.a[t]>>32));
	}
	
	// fprintf(stderr, "p->rids.n:%u, p->idx.n:%u\n", (uint32_t)p->rids.n, (uint32_t)p->idx.n);

	return p;
}


void debug_append_inexact_edges(ma_ug_t *ug, const ug_opt_t *uopt) {
	uint32_t n_asymm = 0, n_disconnect = 0, z, v, w, k, nv; asg_arc_t *av = NULL;
	for (z = 0; z < ug->g->n_arc; ++z) {
		if(ug->g->arc[z].del) continue;
		if(!get_ug_edge_src(ug, uopt->sources, uopt->max_hang, uopt->min_ovlp, 
													ug->g->arc[z].ul>>32, ug->g->arc[z].v)) {
			n_disconnect++;
		}
        v = ug->g->arc[z].v^1; w = (ug->g->arc[z].ul>>32)^1;
        nv = asg_arc_n(ug->g, v); av = asg_arc_a(ug->g, v);
        for (k = 0; k < nv; ++k) {
			if (av[k].del) continue;
			// fprintf(stderr, "found <%lu> -> <%u>\n", av[k].ul>>32, av[k].v);
			if (av[k].v == w) break;
		}
            
        if (k == nv) {
			ug->g->arc[z].del = 1, ++n_asymm;
			// fprintf(stderr, "# lack of <%u> -> <%u>, should be <%u> -> <%u>\n\n", w^1, v^1, v, w);
		}
    }

	if(n_asymm || n_disconnect) {
		asg_cleanup(ug->g); 
		fprintf(stderr, "[M::%s] # asymm edges: %u, # disconnect edges: %u\n", 
		__func__, n_asymm, n_disconnect);
		// exit(1);
	}
}

void append_inexact_edges(ma_ug_t *ug, const ug_opt_t *uopt, asg_t *rg)
{
	uint32_t *idx = NULL, n_read = R_INF.total_reads, z, v, k, qn, tn, tu, ut_v, ut_w; 
	ma_utg_t *u = NULL; ma_hit_t_alloc *src = uopt->sources, *s = NULL;
	int32_t r; asg_arc_t t, *p = NULL;
	int64_t min_ovlp = uopt->min_ovlp, max_hang = uopt->max_hang, occ = 0;

	MALLOC(idx, n_read); memset(idx, -1, n_read*sizeof(*(idx)));
	for (z = 0; z < ug->u.n; z++) {
		u = &(ug->u.a[z]);
		if(u->circ) continue;
		idx[u->start>>1] = idx[u->end>>1] = z;
	}

	for (z = 0; z < ug->u.n; z++) {
		u = &(ug->u.a[z]);
		if(u->circ) continue;

		v = u->end^1; s = &(src[v>>1]); ut_v = (z<<1);
		for (k = 0; k < s->length; k++) {
			if(s->buffer[k].el) continue;///we just need inexact edges
			qn = Get_qn(s->buffer[k]); tn = Get_tn(s->buffer[k]); tu = idx[tn]; ut_w = (uint32_t)-1;
			if(tu == (uint32_t)-1 || ug->g->seq[tu].del) continue;
			if((Get_qe(s->buffer[k]) - Get_qs(s->buffer[k])) < min_ovlp) continue;
            if((Get_te(s->buffer[k]) - Get_ts(s->buffer[k])) < min_ovlp) continue;
			r = ma_hit2arc(&(s->buffer[k]), rg->seq[qn].len, rg->seq[tn].len, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
			if(r < 0 || (t.ul>>32) != v) continue;
			if(t.v == ug->u.a[tu].start) ut_w = tu<<1;
			if(t.v == ug->u.a[tu].end) ut_w = (tu<<1)+1;
			if(ut_w==(uint32_t)-1) continue;
			p = asg_arc_pushp(ug->g); memset(p, 0, sizeof(*p));
			*p = t; p->ul = ut_v; p->ul <<= 32; p->ul += ((uint32_t)(t.ul)); p->v = ut_w;
			occ++;
			// if((p->v>>1)>=ug->g->n_seq || (p->ul>>33)>=ug->g->n_seq) {
			// 	fprintf(stderr, "+ug->g->n_seq:%u, (p->ul>>33):%u, (p->v>>1):%u\n", 
			// 	(uint32_t)ug->g->n_seq, (uint32_t)(p->ul>>33), (uint32_t)(p->v>>1));
			// }
			// assert((p->v>>1)<ug->g->n_seq && (p->ul>>33)<ug->g->n_seq);
		}

		v = u->start^1; s = &(src[v>>1]); ut_v = (z<<1) + 1;
		for (k = 0; k < s->length; k++) {
			if(s->buffer[k].el) continue;///we just need inexact edges
			qn = Get_qn(s->buffer[k]); tn = Get_tn(s->buffer[k]); tu = idx[tn]; ut_w = (uint32_t)-1;
			if(tu == (uint32_t)-1 || ug->g->seq[tu].del) continue;
			if((Get_qe(s->buffer[k]) - Get_qs(s->buffer[k])) < min_ovlp) continue;
            if((Get_te(s->buffer[k]) - Get_ts(s->buffer[k])) < min_ovlp) continue;
			r = ma_hit2arc(&(s->buffer[k]), rg->seq[qn].len, rg->seq[tn].len, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
			if(r < 0 || (t.ul>>32) != v) continue;
			if(t.v == ug->u.a[tu].start) ut_w = tu<<1;
			if(t.v == ug->u.a[tu].end) ut_w = (tu<<1)+1;
			if(ut_w==(uint32_t)-1) continue;
			p = asg_arc_pushp(ug->g); memset(p, 0, sizeof(*p));
			*p = t; p->ul = ut_v; p->ul <<= 32; p->ul += ((uint32_t)(t.ul)); p->v = ut_w;
			occ++;
			// if((p->v>>1)>=ug->g->n_seq || (p->ul>>33)>=ug->g->n_seq) {
			// 	fprintf(stderr, "+ug->g->n_seq:%u, (p->ul>>33):%u, (p->v>>1):%u\n", 
			// 	(uint32_t)ug->g->n_seq, (uint32_t)(p->ul>>33), (uint32_t)(p->v>>1));
			// }
			// assert((p->v>>1)<ug->g->n_seq && (p->ul>>33)<ug->g->n_seq);
		}
	}
	if(occ) {
		free(ug->g->idx);
        ug->g->idx = 0;
        ug->g->is_srt = 0;
		asg_cleanup(ug->g); 
	}
	
	free(idx);
	///for debug
	debug_append_inexact_edges(ug, uopt);
	fprintf(stderr, "[M::%s] # inserted inexact edges: %ld\n", __func__, occ);
}

typedef struct {
	ucov_t *cr;
	ma_hit_t_alloc* src;
	int64_t min_ovlp; 
	int64_t max_hang;
	uint64_t is_el;
	uint64_t is_del;
	uint64_t is_src_cc;
	asg_t *rg;
	ma_ug_t *ug;
} r_contain_aux;

static void update_gen_r_contain(void *data, long i, int tid) // callback for kt_for()
{
	r_contain_aux *s = (r_contain_aux *)data;
	ma_hit_t_alloc *src = s->src; ma_hit_t *t; int32_t r; asg_arc_t x;
	uint64_t *a = s->cr->interval.a + s->cr->idx[i], a_n = s->cr->idx[i+1] - s->cr->idx[i], k, dp, l, z, qn, tn;
	uint64_t is_el = s->is_el, is_del = s->is_del, min_ovlp = s->min_ovlp, max_hang = s->max_hang, qs, qe, cs, ce, sum;
	int64_t ii; asg_t *rg = s->rg; uint64_t *b, b_n, ti;
	// if(a_n == 0 || rg->seq[i].del) return;
	if(s->is_src_cc) {
		for (z = 0; z < src[i].length; z++) {
			t = &(src[i].buffer[z]); t->cc = 0;
			if(a_n == 0 || rg->seq[i].del) continue;
			qn = Get_qn((*t)); tn = Get_tn((*t));
			if(qn > tn) continue;
			if(is_el && (!(t->el))) continue;
			if(is_del && (!(t->del))) continue;
			if((Get_qe((*t)) - Get_qs((*t))) < min_ovlp) continue;
			if((Get_te((*t)) - Get_ts((*t))) < min_ovlp) continue;
			if(rg->seq[tn].del) continue;
			b = s->cr->interval.a + s->cr->idx[tn]; b_n = s->cr->idx[tn+1] - s->cr->idx[tn];
			if(b_n == 0) continue;
			r = ma_hit2arc(t, Get_READ_LENGTH(R_INF, qn), Get_READ_LENGTH(R_INF, tn), max_hang, asm_opt.max_hang_rate, min_ovlp, &x);
			if(r < 0) continue;	
			qs = Get_qs((*t)); qe = Get_qe((*t));
			for (k = 0; k < a_n; k += 2) {
				cs = a[k]>>33; ce = a[k+1]>>33; assert(rg->seq[(uint32_t)(a[k])].del);
				if(qs<=cs+128 && qe+128>=ce) { ///128 is the offset for indel
					for (ti = 0; ti < b_n; ti+=2) {
						if((uint32_t)(b[ti]) == (uint32_t)(a[k])) {
							sum = t->cc; sum += (ce - cs);
							if(sum > 0x3fffffffU) sum = 0x3fffffffU;
							t->cc = sum;
							break;
						}
					}
				}
			}
		}
	} else {
		radix_sort_gfa64(a, a + a_n);
		for (k = 0, dp = 0; k < a_n; ++k) {
			///if a[j] is qe
			if ((a[k]>>32)&1) --dp;
			else ++dp;
			l = a[k]>>33; l <<= 32; l += dp;
			a[k] = l;
		}
		for (z = 0; z < src[i].length; z++) {
			t = &(src[i].buffer[z]); 
			qn = Get_qn((*t)); tn = Get_tn((*t));
			if(qn > tn) continue;
			if(t->cc == 0) continue;
			ii = get_specific_overlap(&(src[tn]), tn, qn);
			src[tn].buffer[ii].cc = t->cc;
		}
	}
}

static void update_ug_uo_t(void *data, long i, int tid)
{
	r_contain_aux *sl = (r_contain_aux *)data; int32_t r;
	ma_hit_t_alloc *src = sl->src, *x; uint32_t k, qn, tn, uv, uw, v, w;
	asg_arc_t *e = &(sl->ug->g->arc[i]), t; 
	uv = e->ul>>32; uw = e->v; e->ou = 0;
	if(sl->ug->u.a[uv>>1].circ || sl->ug->u.a[uw>>1].circ) return;
	v = ((uv&1)?(sl->ug->u.a[uv>>1].start^1):(sl->ug->u.a[uv>>1].end^1));
	w = ((uw&1)?(sl->ug->u.a[uw>>1].end):(sl->ug->u.a[uw>>1].start));
	x = &(src[v>>1]);

	for (k = 0; k < x->length; k++) {
		qn = Get_qn(x->buffer[k]);
		tn = Get_tn(x->buffer[k]);
		if(qn == (v>>1) && tn == (w>>1)) {
			r = ma_hit2arc(&(x->buffer[k]), sl->rg->seq[v>>1].len, sl->rg->seq[w>>1].len,
            							sl->max_hang, asm_opt.max_hang_rate, sl->min_ovlp, &t);
			if(r < 0) continue;
            if((t.ul>>32)!=v || t.v!=w) continue;
			e->ou = (x->buffer[k].cc>OU_MASK?OU_MASK:x->buffer[k].cc);
			break;
		}
	}
	assert(k < x->length);
}

ucov_t *gen_r_contain(ma_ug_t *ug, asg_t *rg, ma_hit_t_alloc* src, uint64_t n_read, int64_t min_ovlp, int64_t max_hang, uint64_t n_thread, uint64_t is_el, uint64_t is_del)
{
	ucov_t *cr = NULL; uint64_t i, z, qn, tn, qs, qe; 
	ma_hit_t *t = NULL; int32_t r; asg_arc_t x;
	CALLOC(cr, 1); MALLOC(cr->idx, n_read+1); kv_init(cr->interval);
	for (i = 0; i < n_read; i++) {
		cr->idx[i] = cr->interval.n;
		if(rg->seq[i].del) continue;
		for (z = 0; z < src[i].length; z++) {
			t = &(src[i].buffer[z]); t->cc = 0;
			if(is_el && (!(t->el))) continue;
			if(is_del && (!(t->del))) continue;
			if((Get_qe((*t)) - Get_qs((*t))) < min_ovlp) continue;
            if((Get_te((*t)) - Get_ts((*t))) < min_ovlp) continue;
			qn = Get_qn((*t)); tn = Get_tn((*t));
			if(!rg->seq[tn].del) continue;
			r = ma_hit2arc(t, Get_READ_LENGTH(R_INF, qn), Get_READ_LENGTH(R_INF, tn), max_hang, asm_opt.max_hang_rate, min_ovlp, &x);
            if(r == MA_HT_TCONT) { ///tn is contained
				qs = Get_qs((*t)); qe = Get_qe((*t));
				kv_push(uint64_t, cr->interval, ((qs<<1)<<32)|tn);
				kv_push(uint64_t, cr->interval, (((qe<<1)|1)<<32)|tn);
			}
		}
	}
	cr->idx[i] = cr->interval.n;
	// fprintf(stderr, "+++[M::%s]n_read:%lu\n", __func__, n_read);
	r_contain_aux aux;
	aux.cr = cr; aux.src = src; aux.min_ovlp = min_ovlp; aux.rg = rg; aux.ug = ug;
	aux.max_hang = max_hang; aux.is_el = 0/**is_el**/; aux.is_del = 0/**is_del**/;
	aux.is_src_cc = 1;
	kt_for(n_thread, update_gen_r_contain, &aux, n_read);///note: here we should set is_el = is_del = 0

	aux.is_src_cc = 0;
	kt_for(n_thread, update_gen_r_contain, &aux, n_read);

	if(ug) kt_for(n_thread, update_ug_uo_t, &aux, ug->g->n_arc);

	return cr;
}

ucov_t *gen_cov_track(ma_ug_t *ug, asg_t *rg, ul_contain *ct, ma_hit_t_alloc* src, int64_t min_ovlp, int64_t max_hang, uint64_t is_el, uint64_t is_del)
{
	uint64_t i, k;
	ucov_t *cc = NULL; CALLOC(cc, 1); MALLOC(cc->idx, ug->u.n+1); kv_init(cc->interval);
	for (i = k = 0; i < ug->u.n; i++) {
		k += ug->u.a[i].len;
		push_coverage_track(cc, ct, i, &(ug->u.a[i]), rg, src, min_ovlp, max_hang, is_el, is_del);
	}
	fprintf(stderr, "[M::%s::] # bases: %lu\n", __func__, k);
	return cc;
}

ul_idx_t *dedup_HiFis(const ug_opt_t *uopt, uint64_t is_el, uint64_t is_del)
{
	uint64_t i, k, qn, tn, n_read = R_INF.total_reads, cc_num = 0;
	int32_t r; asg_arc_t t, *p = NULL; 
	uint8_t *rset = NULL; CALLOC(rset, n_read<<1);
	asg_t *rg = asg_init();
	ma_hit_t_alloc* src = uopt->sources;
	int64_t min_ovlp = uopt->min_ovlp;
	int64_t max_hang = uopt->max_hang;
	int64_t gap_fuzz = uopt->gap_fuzz;

	rg->m_seq = rg->n_seq = n_read; MALLOC(rg->seq, rg->m_seq);
	for (i = 0; i < n_read; ++i) {
		rg->seq[i].len = Get_READ_LENGTH(R_INF, i);
		rg->seq[i].del = rg->seq[i].c = 0;
	}
	
	for (i = 0; i < n_read; i++) {
		if(rg->seq[i].del) continue;
		for (k = 0; k < src[i].length; k++) {
			if(is_el && (!src[i].buffer[k].el)) continue;
			if(is_del && (!src[i].buffer[k].del)) continue;
			qn = Get_qn(src[i].buffer[k]); tn = Get_tn(src[i].buffer[k]);
			if(rg->seq[qn].del || rg->seq[tn].del) continue;
			if((Get_qe(src[i].buffer[k]) - Get_qs(src[i].buffer[k])) < min_ovlp) continue;
			if((Get_te(src[i].buffer[k]) - Get_ts(src[i].buffer[k])) < min_ovlp) continue;
			r = ma_hit2arc(&(src[i].buffer[k]), rg->seq[qn].len, rg->seq[tn].len, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
			if (r == MA_HT_QCONT/** && check_if_fully_contain(qn, tn, src[i].buffer[k].rev, rset, rg, src, min_ovlp, max_hang, gap_fuzz)**/) {
				rg->seq[qn].del = 1;
			} else if(r == MA_HT_TCONT/** && check_if_fully_contain(tn, qn, src[i].buffer[k].rev, rset, rg, src, min_ovlp, max_hang, gap_fuzz)**/) {
				rg->seq[tn].del = 1;
			}
			if(rg->seq[i].del) break;
		}
	}
	

	for (i = 0; i < n_read; i++) {
		if(rg->seq[i].del) {cc_num++; continue;}
		for (k = 0; k < src[i].length; k++) {
			if(is_el && (!src[i].buffer[k].el)) continue;
			if(is_del && (!src[i].buffer[k].del)) continue;
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
	append_inexact_edges(ug, uopt, rg);
	
	ul_idx_t *uu = NULL; CALLOC(uu, 1);
	uu->ug = ug; 
	
	uu->ct = ul_contain_gen(ug, rg, src, min_ovlp, max_hang, is_el, is_del);
	uu->cc = gen_cov_track(ug, rg, uu->ct, src, min_ovlp, max_hang, is_el, is_del);
	uu->cr = gen_r_contain(ug, rg, src, n_read, min_ovlp, max_hang, asm_opt.thread_num, is_el, is_del);
	
	// uu->ov = compress_dedup_HiFis(ug, src);
	asg_destroy(rg); free(rset);
	// uu->nug = cvert_t_gen(uopt);

	fprintf(stderr, "[M::%s::] # unitigs: %lu, # edges: %lu, # cc_num: %lu\n", __func__, (uint64_t)ug->u.n, (uint64_t)ug->g->n_arc, cc_num);
	// print_dedup_HiFis_seq(ug);
	return uu;
}

ul_idx_t *gen_ul_idx(const ug_opt_t *uopt, ma_ug_t *ug, asg_t *sg)
{
	ma_hit_t_alloc* src = uopt->sources;
	int64_t min_ovlp = uopt->min_ovlp;
	int64_t max_hang = uopt->max_hang;
	
	ul_idx_t *uu = NULL; CALLOC(uu, 1); uu->ug = ug; 
	uu->ct = ul_contain_gen(ug, sg, src, min_ovlp, max_hang, 0, 1);
	uu->cc = gen_cov_track(ug, sg, uu->ct, src, min_ovlp, max_hang, 0, 1);
	uu->cr = gen_r_contain(ug, sg, src, R_INF.total_reads, min_ovlp, max_hang, asm_opt.thread_num, 0, 1);
	return uu;
}



utg_rid_t *gen_r_ug_idx(ma_ug_t *ug, asg_t *rg)
{
	uint64_t i, k, l, m, rid, a_n; utg_rid_dt *a; ma_utg_t *u = NULL;
	utg_rid_t *cc = NULL; CALLOC(cc, 1); CALLOC(cc->idx, rg->n_seq+1); kv_init(cc->p); cc->rg = rg;
	for (i = 0; i < ug->u.n; i++) {
		u = &(ug->u.a[i]);
		for (k = 0; k < u->n; k++) cc->idx[u->a[k]>>33]++;
	}
	for (k = l = 0; k <= rg->n_seq; k++) {
        m = cc->idx[k];
        cc->idx[k] = l;
        l += m;
    }
	cc->p.n = cc->p.m = l; CALLOC(cc->p.a, cc->p.n);
	for (i = 0; i < ug->u.n; i++) {
		u = &(ug->u.a[i]);
		for (k = l = 0; k < u->n; k++) {
			rid = u->a[k]>>33;
			a = cc->p.a + cc->idx[rid];
            a_n = cc->idx[rid+1] - cc->idx[rid];
			if(a_n) {
                if(a[a_n-1].off == a_n-1) {
					a[a_n-1].u = (i<<1)|((u->a[k]>>32)&1);
					a[a_n-1].pos = l; a[a_n-1].off = l;					
				}
                else {
					a[a[a_n-1].off].u = (i<<1)|((u->a[k]>>32)&1);
					a[a[a_n-1].off].pos = l; a[a[a_n-1].off].off = l;
					a[a_n-1].off++;
				}
            }
			l += (uint32_t)u->a[k];
		}
	}
	return cc;
}

ul_idx_t *gen_ul_idx_t(const ug_opt_t *uopt, asg_t *sg, uint64_t is_el, uint64_t is_del)
{
	uint64_t n_read = R_INF.total_reads;
	ma_hit_t_alloc* src = uopt->sources;
	int64_t min_ovlp = uopt->min_ovlp;
	int64_t max_hang = uopt->max_hang;
	// int64_t gap_fuzz = uopt->gap_fuzz;
	ul_idx_t *uu = NULL; CALLOC(uu, 1);	uu->ug = ma_ug_gen(sg);
	uu->ct = ul_contain_gen(uu->ug, sg, src, min_ovlp, max_hang, is_el, is_del);
    uu->cc = gen_cov_track(uu->ug, sg, uu->ct, src, min_ovlp, max_hang, is_el, is_del);
	uu->cr = gen_r_contain(uu->ug, sg, src, n_read, min_ovlp, max_hang, asm_opt.thread_num, is_el, is_del);
	uu->r_ug = gen_r_ug_idx(uu->ug, sg);
	return uu;
}


void destroy_ul_idx_t(ul_idx_t *uu)
{
	if(!uu) return;
	if(uu->cc) {
		if(uu->cc) {
			free(uu->cc->idx);
			free(uu->cc->interval.a);
			free(uu->cc);
		}

		if(uu->cr) {
			free(uu->cr->idx);
			free(uu->cr->interval.a);
			free(uu->cr);
		}

		if(uu->ct) {
			free(uu->ct->idx.a);
			free(uu->ct->rids.a);
			free(uu->ct->is_c.a);
			free(uu->ct);
		}

		if(uu->r_ug) {
			free(uu->r_ug->idx);
			free(uu->r_ug->p.a);
			free(uu->r_ug);
		}
		// if(uu->ov) {
		// 	free(uu->ov->a);
		// 	free(uu->ov);
		// }
	}
	ma_ug_destroy(uu->ug);
	// if(uu->nug) {
	// 	free(uu->nug->idx);
	// 	ma_ug_destroy(uu->nug->ug);
	// 	free(uu->nug);
	// }
	free(uu);
}


void gen_UL_ovlps(uldat_t *sl, int32_t cutoff)
{
	ul_idx_t *uu = dedup_HiFis(sl->uopt, 1, 0);
	int exist = (asm_opt.load_index_from_disk? uidx_load(&ha_flt_tab, &ha_idx, asm_opt.output_file_name, NULL) : 0);
    if(exist == 0) uidx_l_build(uu->ug, (mg_idxopt_t *)sl->opt, cutoff);
	if(exist == 0) uidx_write(ha_flt_tab, ha_idx, asm_opt.output_file_name, NULL);
	sl->ha_flt_tab = ha_flt_tab; sl->ha_idx = (ha_pt_t *)ha_idx; sl->uu = uu;		
	ul_v_call(sl, asm_opt.ar);
	destroy_ul_idx_t(uu); ha_ft_destroy(ha_flt_tab); ha_pt_destroy(ha_idx);
	sl->ha_flt_tab = NULL; sl->ha_idx = NULL; sl->uu = NULL;
}


void gen_UL_reovlps(uldat_t *sl, ma_ug_t *ug, asg_t *sg, char* gfa_name, int32_t cutoff)
{
	ul_idx_t *uu = gen_ul_idx(sl->uopt, ug, sg);
	int exist = (asm_opt.load_index_from_disk? uidx_load(&ha_flt_tab, &ha_idx, gfa_name, ug) : 0);
    if(exist == 0) uidx_l_build(uu->ug, (mg_idxopt_t *)sl->opt, cutoff);
	if(exist == 0) uidx_write(ha_flt_tab, ha_idx, gfa_name, ug);
	sl->ha_flt_tab = ha_flt_tab; sl->ha_idx = (ha_pt_t *)ha_idx; sl->uu = uu;	

	init_ucr_file_t(sl, gfa_name, 1);
	ul_v_recall(sl, asm_opt.ar);
	destory_ucr_file_t(sl);
	///do not free ug
	uu->ug = NULL; destroy_ul_idx_t(uu); ha_ft_destroy(ha_flt_tab); ha_pt_destroy(ha_idx);
	sl->ha_flt_tab = NULL; sl->ha_idx = NULL; sl->uu = NULL;
	
	// exit(1);
}

void init_uldat_t(uldat_t *sl, void *ha_flt_tab, void *ha_idx, mg_idxopt_t *opt, uint64_t chunk_size, uint64_t n_thread, const ug_opt_t *uopt, ul_idx_t *uu)
{
	memset(sl, 0, sizeof(uldat_t));
    sl->ha_flt_tab = ha_flt_tab;
    sl->ha_idx = (ha_pt_t *)ha_idx;
    sl->opt = opt;   
    sl->chunk_size = chunk_size;
    sl->n_thread = n_thread;
    sl->uu = uu;
    sl->uopt = uopt;
}

int32_t write_all_ul_t(all_ul_t *x, char* file_name, ma_ug_t *ug)
{
	char* gfa_name = NULL; MALLOC(gfa_name, strlen(file_name)+50);
	sprintf(gfa_name, "%s.ul.ovlp.bin", file_name);
	FILE* fp = fopen(gfa_name, "w"); free(gfa_name);
	if (!fp) return 0;
	uint64_t k; ul_vec_t *p = NULL;

	if(ug) write_dbug(ug, fp);

	fwrite(&x->nid.n, sizeof(x->nid.n), 1, fp);
	for (k = 0; k < x->nid.n; k++) {
		fwrite(&x->nid.a[k].n, sizeof(x->nid.a[k].n), 1, fp);
		fwrite(x->nid.a[k].a, sizeof((*(x->nid.a[k].a))), x->nid.a[k].n, fp);
	}

	fwrite(&x->ridx.idx.n, sizeof(x->ridx.idx.n), 1, fp);
	fwrite(x->ridx.idx.a, sizeof((*(x->ridx.idx.a))), x->ridx.idx.n, fp);

	fwrite(&x->ridx.occ.n, sizeof(x->ridx.occ.n), 1, fp);
	fwrite(x->ridx.occ.a, sizeof((*(x->ridx.occ.a))), x->ridx.occ.n, fp);

	fwrite(&x->n, sizeof(x->n), 1, fp);
	for (k = 0; k < x->n; k++) {
		p = &(x->a[k]);
		fwrite(&p->dd, sizeof(p->dd), 1, fp);
		fwrite(&p->rlen, sizeof(p->rlen), 1, fp);

		fwrite(&p->r_base.n, sizeof(p->r_base.n), 1, fp);
		fwrite(p->r_base.a, sizeof((*(p->r_base.a))), p->r_base.n, fp);

		fwrite(&p->bb.n, sizeof(p->bb.n), 1, fp);
		fwrite(p->bb.a, sizeof((*(p->bb.a))), p->bb.n, fp);

		fwrite(&p->N_site.n, sizeof(p->N_site.n), 1, fp);
		fwrite(p->N_site.a, sizeof((*(p->N_site.a))), p->N_site.n, fp);
	}

	fprintf(stderr, "[M::%s] Index has been written.\n", __func__);
	fclose(fp);
	return 1;
}


int32_t load_all_ul_t(all_ul_t *x, char* file_name, All_reads *hR, ma_ug_t *ug)
{
	char* gfa_name = NULL; MALLOC(gfa_name, strlen(file_name)+50);
	sprintf(gfa_name, "%s.ul.ovlp.bin", file_name);
	FILE* fp = fopen(gfa_name, "r"); free(gfa_name);
	if (!fp) return 0;
	if(ug && (!test_dbug(ug, fp))) {
		fprintf(stderr, "[M::%s] Renew UL Index\n", __func__);
		fclose(fp);
		return 0;
    }

	memset(x, 0, sizeof(*x)); x->hR = hR; init_aux_table();
	uint64_t k; ul_vec_t *p = NULL;
	
	fread(&x->nid.n, sizeof(x->nid.n), 1, fp); x->nid.m = x->nid.n; MALLOC(x->nid.a, x->nid.n); 
	for (k = 0; k < x->nid.n; k++) {
		fread(&x->nid.a[k].n, sizeof(x->nid.a[k].n), 1, fp); MALLOC(x->nid.a[k].a, x->nid.a[k].n);
		fread(x->nid.a[k].a, sizeof((*(x->nid.a[k].a))), x->nid.a[k].n, fp);
	}

	fread(&x->ridx.idx.n, sizeof(x->ridx.idx.n), 1, fp); x->ridx.idx.m = x->ridx.idx.n; MALLOC(x->ridx.idx.a, x->ridx.idx.n);
	fread(x->ridx.idx.a, sizeof((*(x->ridx.idx.a))), x->ridx.idx.n, fp);

	fread(&x->ridx.occ.n, sizeof(x->ridx.occ.n), 1, fp); x->ridx.occ.m = x->ridx.occ.n; MALLOC(x->ridx.occ.a, x->ridx.occ.n);
	fread(x->ridx.occ.a, sizeof((*(x->ridx.occ.a))), x->ridx.occ.n, fp);

	fread(&x->n, sizeof(x->n), 1, fp); x->m = x->n; MALLOC(x->a, x->n);
	for (k = 0; k < x->n; k++) {
		p = &(x->a[k]);
		fread(&p->dd, sizeof(p->dd), 1, fp);
		fread(&p->rlen, sizeof(p->rlen), 1, fp);

		fread(&p->r_base.n, sizeof(p->r_base.n), 1, fp); p->r_base.m = p->r_base.n; MALLOC(p->r_base.a, p->r_base.n);
		fread(p->r_base.a, sizeof((*(p->r_base.a))), p->r_base.n, fp);

		fread(&p->bb.n, sizeof(p->bb.n), 1, fp); p->bb.m = p->bb.n; MALLOC(p->bb.a, p->bb.n);
		fread(p->bb.a, sizeof((*(p->bb.a))), p->bb.n, fp);

		fread(&p->N_site.n, sizeof(p->N_site.n), 1, fp); p->N_site.m = p->N_site.n; MALLOC(p->N_site.a, p->N_site.n);
		fread(p->N_site.a, sizeof((*(p->N_site.a))), p->N_site.n, fp);
	}

	fprintf(stderr, "[M::%s] Index has been loaded.\n", __func__);
	fclose(fp);
	return 1;
}

void ul_load(const ug_opt_t *uopt)
{
	fprintf(stderr, "[M::%s::] ==> UL\n", __func__);
	mg_idxopt_t opt; uldat_t sl;
	int32_t cutoff;
	init_aux_table(); ha_opt_update_cov(&asm_opt, asm_opt.hom_cov);
	cutoff = asm_opt.max_n_chain;
	init_mg_opt(&opt, !(asm_opt.flag&HA_F_NO_HPC), 19, 10, cutoff, asm_opt.max_n_chain, asm_opt.ul_error_rate, asm_opt.ul_error_rate, asm_opt.ul_error_rate_low, asm_opt.ul_ec_round);
	init_uldat_t(&sl, NULL, NULL, &opt, CHUNK_SIZE, asm_opt.thread_num, uopt, NULL);

	if(!load_all_ul_t(&UL_INF, asm_opt.output_file_name, &R_INF, NULL)) {
		gen_UL_ovlps(&sl, cutoff);
		write_all_ul_t(&UL_INF, asm_opt.output_file_name, NULL);
	}

	// print_all_ul_t_stat(&UL_INF);
	// fprintf(stderr, "**1**\n");
	kt_for(sl.n_thread, update_ovlp_src, &sl, R_INF.total_reads);
	// fprintf(stderr, "**2**\n");
	kt_for(sl.n_thread, update_ovlp_src_bl, &sl, R_INF.total_reads);
	// fprintf(stderr, "**3**\n");
	
	print_ovlp_src_bl_stat(&UL_INF, sl.uopt);
	// print_ul_ovlps(&UL_INF, 0); print_ul_ovlps(&UL_INF, 1);

	// destory_all_ul_t(&UL_INF); 
}


uint64_t ul_refine_alignment(const ug_opt_t *uopt, asg_t *sg)
{
	fprintf(stderr, "[M::%s::] ==> UL refinement...\n", __func__);
	mg_idxopt_t opt; uldat_t sl; int32_t cutoff;
	init_aux_table(); ha_opt_update_cov(&asm_opt, asm_opt.hom_cov);
	cutoff = asm_opt.max_n_chain;
	init_mg_opt(&opt, !(asm_opt.flag&HA_F_NO_HPC), 19, 10, cutoff, asm_opt.max_n_chain, asm_opt.ul_error_rate, asm_opt.ul_error_rate, asm_opt.ul_error_rate_low, asm_opt.ul_ec_round);
	ul_idx_t *uu = gen_ul_idx_t(uopt, sg, 0, 0);///record contained reads; is_el = is_del = 0	
	init_uldat_t(&sl, NULL, NULL, &opt, CHUNK_SIZE, asm_opt.thread_num, uopt, uu); sl.rg = sg;
	if(work_ul_gchains(&sl)) {
		free(UL_INF.ridx.idx.a); free(UL_INF.ridx.occ.a); memset(&(UL_INF.ridx), 0, sizeof(UL_INF.ridx));
		gen_ul_vec_rid_t(&UL_INF);
		kt_for(sl.n_thread, update_ovlp_src, &sl, R_INF.total_reads);
		kt_for(sl.n_thread, update_ovlp_src_bl, &sl, R_INF.total_reads);
		destroy_ul_idx_t(uu);
		return 1;
	} else {
		destroy_ul_idx_t(uu);
		return 0;
	}
}

uint32_t dd_ug(asg_t *sg, ma_ug_t *ug, ma_sub_t* coverage_cut, ma_hit_t_alloc* sources, R_to_U* ruIndex, const char* output_file_name)
{
    fprintf(stderr, "Writing raw unitig GFA to disk... \n");
    char* gfa_name = (char*)malloc(strlen(output_file_name)+25);
    sprintf(gfa_name, "%s.r_utg.noseq.gfa", output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    ma_ug_print_simple(ug, sg, coverage_cut, sources, ruIndex, "utg", output_file);
    fclose(output_file);

    free(gfa_name);
    exit(1);
}



ma_ug_t *ul_realignment(const ug_opt_t *uopt, asg_t *sg)
{
	fprintf(stderr, "[M::%s::] ==> UL\n", __func__);
	mg_idxopt_t opt; uldat_t sl;
	int32_t cutoff;
	char* gfa_name = NULL; MALLOC(gfa_name, strlen(asm_opt.output_file_name)+50);
	sprintf(gfa_name, "%s.re", asm_opt.output_file_name);

	init_aux_table(); ha_opt_update_cov(&asm_opt, asm_opt.hom_cov);
	cutoff = REA_ALIGN_CUTOFF;
	init_mg_opt(&opt, !(asm_opt.flag&HA_F_NO_HPC), 19, 10, cutoff, asm_opt.max_n_chain, asm_opt.ul_error_rate, asm_opt.ul_error_rate, asm_opt.ul_error_rate_low, asm_opt.ul_ec_round);
	init_uldat_t(&sl, NULL, NULL, &opt, CHUNK_SIZE, asm_opt.thread_num, uopt, NULL);
	ma_ug_t *ug = gen_polished_ug(uopt, sg);
	// dd_ug(sg, ug, uopt->coverage_cut, uopt->sources, uopt->ruIndex, "UL.sa");
	// debug_sl_compress_base_disk_0(&sl, asm_opt.ar);

	if(!load_all_ul_t(&UL_INF, gfa_name, &R_INF, ug)) {
		gen_UL_reovlps(&sl, ug, sg, gfa_name, cutoff);
		write_all_ul_t(&UL_INF, gfa_name, ug);
	}

	// print_all_ul_t_stat(&UL_INF);
	// kt_for(sl.n_thread, update_ovlp_src, &sl, R_INF.total_reads);
	// kt_for(sl.n_thread, update_ovlp_src_bl, &sl, R_INF.total_reads);
	
	// print_ovlp_src_bl_stat(&UL_INF, sl.uopt);
	// print_ul_ovlps(&UL_INF, 0); print_ul_ovlps(&UL_INF, 1);

	// destory_all_ul_t(&UL_INF); 
	free(gfa_name);
	return ug;
}