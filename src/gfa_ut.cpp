#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <zlib.h>
#include <math.h>
#include "kdq.h"
#include "kthread.h"
#include "gfa_ut.h" 
#include "CommandLines.h"
#include "Correct.h"
#include "inter.h"
#include "Overlaps.h"
#include "hic.h"
#include "Purge_Dups.h"

#define generic_key(x) (x)
KRADIX_SORT_INIT(srt64, uint64_t, generic_key, 8)
#define OU_NOISY 2
#define ASG_ET_MERGEABLE 0
#define ASG_ET_TIP       1
#define ASG_ET_MULTI_OUT 2
#define ASG_ET_MULTI_NEI 3
#define UL_TRAV_HERATE 0.2
#define UL_TRAV_FT_RATE 0.8

KDQ_INIT(uint64_t)

typedef struct { size_t n, m; char *a; } asgc8_v;

typedef struct {
    asg64_v cnt;
    asg64_v idx_a; 
    uint64_t idx_n;
    asg_t *ext;
    uint64_t a_n;
} asg_ext_t;

typedef struct {
    uint32_t v, uid, off;
} usg_arc_mm_t;

typedef struct {
    size_t n, m;
    usg_arc_mm_t *a;
} usg_arc_mm_warp;

typedef struct {
    uint64_t ul;
    uint32_t v;
    uint32_t ol:31, del:1;
    uint32_t ou;
    uint64_t idx;
} usg_arc_t;

typedef struct {
    size_t n, m;
    usg_arc_t *a; 
} usg_arc_warp;

typedef struct {
    uint32_t mm, occ;
    uint32_t len;
    usg_arc_warp arc[2];
    usg_arc_mm_warp arc_mm[2];
    uint8_t del, telo;
} usg_seq_t;

#define usg_arc_key(p) ((p).v)
KRADIX_SORT_INIT(usg_arc_srt, usg_arc_t, usg_arc_key, member_size(usg_arc_t, v))

#define usg_arc_mm_key(p) ((p).v)
KRADIX_SORT_INIT(usg_arc_mm_srt, usg_arc_mm_t, usg_arc_mm_key, member_size(usg_arc_mm_t, v))

typedef struct {
    uint32_t *a;
    size_t n, m;
} mmap_t;

typedef struct {
    usg_seq_t *a;
    size_t n, m;
    kvec_t(mmap_t) mp;
} usg_t;

#define usg_arc_a(g, v) ((g)->a[(v)>>1].arc[(v)&1].a)
#define usg_arc_n(g, v) ((g)->a[(v)>>1].arc[(v)&1].n)

typedef struct{
    int64_t tipsLen; 
    float tip_drop_ratio; 
    int64_t stops_threshold; 
    float chimeric_rate;
    float drop_ratio; 

    bub_label_t* b_mask_t;
    int64_t clean_round;
    double min_ovlp_drop_ratio;
    double max_ovlp_drop_ratio;
    double hom_check_drop_rate;
    double min_path_drop_ratio;
    double max_path_drop_ratio;
    int64_t max_tip, max_tip_hifi;
    uint32_t is_trio;
}ulg_opt_t;

typedef struct {
    uint32_t hid;
    uint32_t qs, qe, ts, te;
    uint32_t qs_k, qe_k, ts_k, te_k;
    uint8_t is_rev:6, is_del:1, is_ct:1;
} ul2ul_t;

typedef struct {
    ul2ul_t *a;
    size_t n, m;
    uint32_t id:31, is_del:1;
    uint32_t cn;
    // uint8_t is_consist;
} ul2ul_item_t;

typedef struct {
    uint32_t v, s, e, n;
} uinfo_srt_t;

typedef struct {
    size_t n, m;
    uinfo_srt_t *a;
} uinfo_srt_warp_t;

typedef struct {
    uint32_t *uc, *hc, *raw_uc;
    uinfo_srt_warp_t *iug_a;
    uint32_t *iug_idx;
    uint64_t *iug_b;
} ul_cov_t;


typedef struct {
    asg_t *bg;
    uint32_t *w_n, *a_n;
} ul_bg_t;

// typedef struct {
//     uint32_t v, n, wn;
// } ul_tra_t;

// typedef struct {
//     kvec_t(ul_tra_t) arc;
//     kvec_t(uint32_t) idx;
// } ul_tra_idx_t;

// #define iug_tra_arc_n(z, v) ((z)->idx.a[(v)+1]-(z)->idx.a[(v)])
// #define iug_tra_arc_a(z, v) ((z)->arc.a + (z)->idx.a[(v)])

typedef struct {
    ul2ul_item_t *a;
    size_t n, m;
    uint64_t uln, gn, tot;
    uint32_t *item_idx;
    asg_t *i_g; ma_ug_t *i_ug;
    ma_ug_t *hybrid_ug;
    usg_t *h_usg;
    ul_cov_t cc;
    ul_bg_t bg;
    asg64_v *iug_tra;
    uinfo_srt_warp_t *iug_seq; uint64_t iug_cov_thre;
    uint8_t *telo;
} ul2ul_idx_t;

#define ul2ul_srt_key(p) ((p).hid)
KRADIX_SORT_INIT(ul2ul_srt, ul2ul_t, ul2ul_srt_key, member_size(ul2ul_t, hid))

typedef struct {
    kvec_t(uint64_t) ref;
    kvec_t(uint64_t) pat;
    kvec_t(uint64_t) pat_cor;
    kvec_t(uint8_t) g_flt;
    kvec_t(uint8_t) m_dir;
    kvec_t(int64_t) m_score;
    uint64_t n, m;
} path_dp_t;

typedef struct {
    uint32_t p; // the optimal parent vertex
    uint32_t d; // the shortest distance from the initial vertex
    uint64_t c; // max count of positive reads
    uint64_t m; // max count of negative reads
    // uint32_t np; // max count of non-positive reads
    uint32_t nc; // max count of reads, no matter positive or negative
    uint32_t r:31, s:1; // r: the number of remaining incoming arc; s: state
    //s: state, s=0, this edge has not been visited, otherwise, s=1
} uinfo_t;



// #define uinfo_srt_t_c_key(p) ((p).se)
// KRADIX_SORT_INIT(uinfo_srt_t_c, uinfo_srt_t, uinfo_srt_t_c_key, member_size(uinfo_srt_t, se))

typedef struct {
    ///all information for each node
    kvec_t(uinfo_t) a;
    // kvec_t(uint32_t) u;
    kvec_t(uint32_t) S; // set of vertices without parents, nodes with all incoming edges visited
    kvec_t(uint32_t) T; // set of tips
    kvec_t(uint32_t) b; // visited vertices
    kvec_t(uint32_t) e; // visited edges/arcs
    // kvec_t(uinfo_srt_t) srt;
    kvec_t(uint8_t) us;
    path_dp_t dp;
} ubuf_t;

typedef struct{
    uint32_t bid, beg, occ;
    uint32_t n_path, path_idx, path_occ;
}ul_sub_path_t;


typedef struct{
    ma_ug_t *buf_ug;
    kvec_t(uint64_t) buf;
}ul_path_t;

typedef struct{
    kvec_t(uint64_t) idx;
    kvec_t(uint64_t) srt;
    uint64_t ul_n;
}ul_path_srt_t;

typedef struct{
    size_t n, m;
    uint64_t *a;
    uint32_t cn:31, is_cir:1;
}ul_str_t;

typedef struct{
    kvec_t(uint64_t) idx;
    kvec_t(uint64_t) occ;
    kvec_t(ul_str_t) str;
}ul_str_idx_t;

typedef struct{
    uint32_t v, pi, ai;
    uint32_t k:31, is_gc:1;
    int32_t dis;
} integer_seq_t;

typedef struct{
    size_t n, m;
    integer_seq_t *a;
} kv_integer_seq_t;

typedef struct{
    uint32_t tk, vq, sc;
    uint64_t tn_rev_qk;
} integer_aln_t;

#define integer_aln_t_vqk_key(x) ((x).tn_rev_qk)
KRADIX_SORT_INIT(integer_aln_t_srt, integer_aln_t, integer_aln_t_vqk_key, member_size(integer_aln_t, tn_rev_qk))

typedef struct {
    size_t n, m;
    integer_aln_t *a;
} integer_aln_vec_t;

typedef struct{
    uint32_t s, e, v;
    uint64_t sc;
    uint32_t q_sidx, q_eidx;
    uint32_t t_sidx, t_eidx;
} ul_chain_t;



typedef struct{
    uint64_t qidx_occ;
    uint32_t tidx_occ;
    uint32_t chain_id:31, is_rev:1;
} ul_snp_t;

#define ul_snp_t_srt_key(x) ((x).qidx_occ)
KRADIX_SORT_INIT(ul_snp_t_srt, ul_snp_t, ul_snp_t_srt_key, member_size(ul_snp_t, qidx_occ))


typedef struct {
    uint32_t occ, nid;
} poa_nid_t;

typedef struct {
    uint64_t ul;
    uint32_t v;
} poa_arc_t;

#define poa_arc_key(a) ((a).ul)
KRADIX_SORT_INIT(poa_arc_srt, poa_arc_t, poa_arc_key, member_size(poa_arc_t, ul))

typedef struct {
    kvec_t(uint32_t) ind;
    kvec_t(uint32_t) stack;
    kvec_t(uint32_t) res;
    kvec_t(uint32_t) res2nid;
    kvec_t(uint64_t) aln;
} topo_srt_t;

typedef struct {
    kvec_t(int64_t) sc;
    kvec_t(uint8_t) dir;
    kvec_t(uint64_t) prefix;
    uint64_t n, m;
} poa_dp_t;

#define poa_dp_idx(dp, x, y) ((dp).m*(x)+(y))
#define e_pdp 0
#define ue_pdp 1
#define lstr_dp 2
#define lg_dp 3

typedef struct{
    uint64_t pge, ule;
    uint32_t ulid;
} emap_t;

#define emap_t_srt_key(x) ((x).pge)
KRADIX_SORT_INIT(emap_t_srt, emap_t, emap_t_srt_key, member_size(emap_t, pge))

typedef struct {
    kvec_t(poa_nid_t) seq;
    kvec_t(poa_arc_t) arc;
    kvec_t(uint64_t) idx;
    uint32_t update_seq;
    uint32_t update_arc;
    topo_srt_t srt_b;
    // poa_dp_t dp;
    kvec_t(emap_t) e_idx;
    ubuf_t bb;
} poa_g_t;

typedef struct {
    kv_integer_seq_t q;
    kv_integer_seq_t t;
    integer_aln_vec_t b;
    kvec_t(int64_t) f;
    kvec_t(int64_t) p;
    kvec_t(uint64_t) o;
    kvec_t(uint64_t) u;
    kvec_t(uint32_t) vis;
    // kvec_t(uint64_t) srt;
    // kvec_t(uint64_t) v;
    // kvec_t(uint64_t) u;
    // kvec_t(uint64_t) d;
    kvec_t(ul_chain_t) sc;
    kvec_t(ul_snp_t) snp;
    poa_g_t pg;
    kvec_t(uint64_t) res_dump;
    uint64_t n_correct, n_circle;
}integer_t;

typedef struct {
    // ul_resolve_t *u;
    integer_t *buf;
    uint64_t n_thread;
}integer_ml_t;

typedef struct {
    asg_t *g;
    uint64_t n[2], tot;
} cul_g_t;

typedef struct{
    ug_opt_t *uopt;
    ma_ug_t *init_ug;
    ma_ug_t *l0_ug;
    ma_ug_t *l1_ug;
    asg_t *sg;
    bubble_type *bub;
    all_ul_t *idx;
    ul_path_t path;
    uint8_t *r_het;
    ubuf_t buf;
    ul_str_idx_t pstr;
    integer_ml_t str_b;
    cul_g_t *cg;
    ul2ul_idx_t uovl;
    // ul_path_srt_t psrt;
}ul_resolve_t;

typedef struct{
    ug_opt_t *uopt;
    asg_t *sg;
    ma_ug_t *ug;
    buf_t b;
    uint32_t *idx, *bid;
    uint32_t idx_n, bid_n;
    uint8_t is_ou;
    uint8_t is_trio;
    int64_t max_ext;
    uint64_t tlen;
    double len_rat;
    double ou_rat;
    int64_t min_ou;
}ug_clean_t;

void deep_graph_clean(ug_opt_t *uopt, asg_t *sg, uint8_t is_ou, uint8_t is_trio, int64_t max_ext, 
double min_ovlp_drop_ratio, double max_ovlp_drop_ratio, double ou_rat, int64_t min_ou, int64_t clean_round, int64_t long_tip);

void init_integer_ml_t(integer_ml_t *x, ul_resolve_t *u, uint64_t n_thread)
{
    memset(x, 0, sizeof((*x))); 
    x->n_thread = n_thread; ///x->u = u;
    CALLOC(x->buf, n_thread);
}

int32_t if_sup_chimeric(ma_hit_t_alloc* src, uint64_t rLen, asg64_v *b, int if_exact);

void print_edge(asg_arc_t *t, const char *cmd)
{
    uint32_t v = t->ul>>32, w = t->v;
    fprintf(stderr, "%s: v->%u(%c)[%u], w->%u(%c)[%u], el->%u, del->%u\n", cmd, v>>1, "+-"[v&1], v, w>>1, "+-"[w&1], w, t->el, t->del);
}

void stats_chimeric(asg_t *g, ma_hit_t_alloc* src, asg64_v *in)
{
    asg64_v tx = {0,0,0}, *b = NULL;
    uint32_t v, s[2] = {0};
    if(in) b = in;
    else b = &tx;
    b->n = 0;

    for (v = 0; v < g->n_seq; ++v) {
        if (g->seq[v].del) continue;
        s[if_sup_chimeric(&(src[v]), g->seq[v].len, b, 1)]++;
    }

    fprintf(stderr, "[M::%s::] ==> # non-chimeric:%u, # chimeric:%u\n", __func__, s[0], s[1]);
    if(!in) free(tx.a);
}

static void stats_sysm_worker(void *_data, long eid, int tid)
{
    asg_t *g = (asg_t*)_data;
    asg_arc_t *p = &(g->arc[eid]);
    if(p->del) return;
    uint32_t k, v = p->v^1, w = (p->ul>>32)^1, nv; asg_arc_t *av;
    av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
    for (k = 0; k < nv; k++) {
        if (av[k].del || av[k].v!=w) continue;
        break;
    }
    assert(k < nv);

    v = p->ul>>32; w = p->v;
    av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
    for (k = 0; k < nv; k++) {
        if (av[k].del || av[k].v!=w) continue;
        assert((uint32_t)eid == av-g->arc+k);
    }
}

void stats_sysm(asg_t *g) {
    kt_for(asm_opt.thread_num, stats_sysm_worker, g, g->n_arc);
    fprintf(stderr, "[M::%s::]",  __func__);
}

uint32_t get_arcs(asg_t *g, uint32_t v, uint32_t* idx, uint32_t idx_n)
{
    uint32_t i, kv = 0, an = asg_arc_n(g, v), beg = g->idx[v]>>32;
    for (i = 0, kv = 0; i < an; i++) {
        if(g->arc[beg+i].del) continue;
        if(idx && kv<idx_n) idx[kv] = beg+i;
        kv++;
    }
    return kv;
}

#define flex_arcs0(res, fg, id) ((res) = ((!((id)&(((uint32_t)(0x80000000)))))?(&((fg).g->arc[(id)])):(&((fg).a[(id)-(((uint32_t)(0x80000000)))]))));

uint32_t get_flex_arcs(flex_asg_t *fg, uint32_t v, uint32_t* idx, uint32_t idx_n)
{
    asg_t *g = fg->g; uint32_t i, kv = 0;
    uint32_t an = asg_arc_n(g, v), beg = g->idx[v]>>32;
    for (i = 0, kv = 0; i < an; i++) {
        if(g->arc[beg+i].del) continue;
        if(idx && kv<idx_n) idx[kv] = beg+i;
        kv++;
    }
    for(i = fg->idx[v]; i != ((uint32_t)-1); i = fg->pi.a[i]) {
        if(fg->a[i].del) continue;
        if(idx && kv<idx_n) idx[kv] = (i)|((uint32_t)(0x80000000));
        kv++;
    }
    return kv;
}

uint32_t follow_limit_path(asg_t *g, uint32_t s, uint32_t *e, uint32_t *occ, asg64_v *b, uint32_t lim)
{
    
    uint32_t v = s, w = 0;
    uint32_t kv, kw;
    (*occ) = 0;


    while (1) {
        (*occ)++;
        kv = get_arcs(g, v, &w, 1);        
        (*e) = v;

        ///if(b) kv_push(uint32_t, b->b, v>>1);
        if(b) kv_push(uint64_t, *b, v);

        if(kv == 0) return END_TIPS;
        if(kv == 2) return TWO_OUTPUT;
        if(kv > 2) return MUL_OUTPUT;
        if((*occ) > lim) return LONG_TIPS;
        w = g->arc[w].v;
        ///up to here, kv=1
        ///kw must >= 1
        kw = get_arcs(g, w^1, NULL, 0);
        v = w;

        if(kw == 2) return TWO_INPUT;
        if(kw > 2) return MUL_INPUT;
        if(v == s) return LOOP;
    }

    return LONG_TIPS;   
}

static inline int asg_end(const asg_t *g, uint32_t v, uint64_t *lw, uint32_t *ou)
{
	///v^1 is the another direction of v
	uint32_t w, nv, nw, nw0, nv0 = asg_arc_n(g, v^1);
	int i, i0 = -1;
	asg_arc_t *aw, *av = asg_arc_a(g, v^1);

	///if this arc has not been deleted
	for (i = nv = 0; i < (int)nv0; ++i)
		if (!av[i].del) i0 = i, ++nv;

	///end without any out-degree
	if (nv == 0) return ASG_ET_TIP; // tip
	if (nv > 1) return ASG_ET_MULTI_OUT; // multiple outgoing arcs
	///until here, nv == 1
	if (lw) *lw = av[i0].ul<<32 | av[i0].v;
	if (ou) *ou = av[i0].ou;
	w = av[i0].v ^ 1;
	nw0 = asg_arc_n(g, w);
	aw = asg_arc_a(g, w);
	for (i = nw = 0; i < (int)nw0; ++i)
		if (!aw[i].del) ++nw;

	if (nw != 1) return ASG_ET_MULTI_NEI;
	return ASG_ET_MERGEABLE;
}

uint32_t asg_arc_cut_tips(asg_t *g, uint32_t max_ext, asg64_v *in, uint32_t is_ou, R_to_U *ru, telo_end_t *te)
{
    asg64_v tx = {0,0,0}, *b = NULL;
    uint32_t n_vtx = g->n_seq<<1, v, w, i, k, cnt = 0, nv, kv, pb, ou, mm_ou, rr, is_u, is_telo;
    asg_arc_t *av = NULL; uint64_t lw;
	if(in) b = in;
    else b = &tx;
    b->n = 0;
    for (v = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;
        if(te && te->hh[v>>1]) continue;

        av = asg_arc_a(g, v^1); nv = asg_arc_n(g, v^1);
        for (i = kv = 0; i < nv; i++) {
            if (av[i].del) continue;
            kv++; break;
        }
        if(kv) continue;
        
        kv = 1; mm_ou = (uint32_t)-1; ou = 0; is_telo = 0;
        if(te && te->hh[v>>1]) is_telo = 1;
        for (i = 0, w = v; i < max_ext; i++) {
            if(asg_end(g, w^1, &lw, is_ou?&ou:NULL)!=0) break;
            w = (uint32_t)lw; kv++; mm_ou = MIN(mm_ou, ou);
            if(te && te->hh[w>>1]) is_telo = 1;
        }

		if(mm_ou == (uint32_t)-1) mm_ou = 0;
		kv += mm_ou; i += mm_ou;
        if((i < max_ext/** + (!!is_ou)**/) && (!is_telo)) kv_push(uint64_t, *b, (((uint64_t)kv)<<32)|v);  
    }

    radix_sort_srt64(b->a, b->a + b->n);

    for (k = 0; k < b->n; k++) {
        v = (uint32_t)(b->a[k]);

        if (g->seq[v>>1].del) continue;
        if(te && te->hh[v>>1]) continue;

        av = asg_arc_a(g, v^1); nv = asg_arc_n(g, v^1);
        for (i = kv = 0; i < nv; i++) {
            if (av[i].del) continue;
            kv++; break;
        }

        if(kv) continue;
        
        pb = b->n; kv_push(uint64_t, *b, v); mm_ou = (uint32_t)-1; ou = 0; is_telo = 0;
        if(te && te->hh[v>>1]) is_telo = 1;
        for (i = 0, w = v; i < max_ext; i++) {
            if(asg_end(g, w^1, &lw, is_ou?&ou:NULL)!=0) break;
            w = (uint32_t)lw; kv_push(uint64_t, *b, lw); mm_ou = MIN(mm_ou, ou);
            if(te && te->hh[w>>1]) is_telo = 1;
        }
		if(mm_ou == (uint32_t)-1) mm_ou = 0;
		i += mm_ou;

        if((i < max_ext/** + (!!is_ou)**/) && (!is_telo)) {
            for (i = pb; i < b->n; i++) asg_seq_del(g, ((uint32_t)b->a[i])>>1);
            cnt++;
        }
        b->n = pb; 
    }

    if(ru && is_ou) {
        for (v = b->n = 0; v < n_vtx; ++v) {
            if (g->seq[v>>1].del) continue;
            av = asg_arc_a(g, v^1); nv = asg_arc_n(g, v^1);
            for (i = kv = 0; i < nv; i++) {
                if (av[i].del) continue;
                kv++; break;
            }
            if(kv) continue;


            get_R_to_U(ru, v>>1, &rr, &is_u);
            if(rr == (uint32_t)-1 || is_u == 1) continue;
            if(te && te->hh[v>>1]) continue;
            kv_push(uint64_t, *b, v); 

            for (i = 0, w = v; i < max_ext; i++) {
                if(asg_end(g, w^1, &lw, NULL)!=0) break;
                w = (uint32_t)lw; 

                get_R_to_U(ru, w>>1, &rr, &is_u);
                if(rr == (uint32_t)-1 || is_u == 1) break;
                if(te && te->hh[w>>1]) break;

                kv_push(uint64_t, *b, lw);
            }

            for (i = 0; i < b->n; i++) {
                asg_seq_del(g, ((uint32_t)b->a[i])>>1);
            }
            if(b->n) cnt++;
        }
    }
    

    /**
    for (v = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;

        av = asg_arc_a(g, v^1); nv = asg_arc_n(g, v^1);
        for (i = kv = 0; i < nv; i++) {
            if (av[i].del) continue;
            kv++;
        }

        if(kv) continue;
        pb = b->n; kv_push(uint64_t, *b, v);
        for (i = 0, w = v; i < max_ext; i++) {
            if(asg_is_utg_end(g, w^1, &lw)!=0) break;
            w = (uint32_t)lw; kv_push(uint64_t, *b, lw);
        }

        if(i < max_ext) {
            for (i = pb; i < b->n; i++) asg_seq_del(g, ((uint32_t)b->a[i])>>1);
            cnt++;
        }
        b->n = pb; 
    }
    **/
    // stats_sysm(g);
    if(!in) free(tx.a);
    if (cnt > 0) asg_cleanup(g);

    return cnt;
}

static void update_sg_contain(void *data, long i, int tid)
{
    sset_aux *sl = (sset_aux *)data; ma_hit_t *h, *z; asg_arc_t t;
	ma_hit_t_alloc *src = sl->src; asg_t *g = sl->g; int64_t r, idx;
    R_to_U *ridx = sl->ruIndex; uint32_t k, rr, qn, tn, is_u;
    g->seq_vis[i] = 0; 
    if(!(g->seq[i].del)) return;
    get_R_to_U(ridx, i, &rr, &is_u);
    if(rr == (uint32_t)-1 || is_u == 1) return;
    ma_hit_t_alloc *x = &(src[i]);
    for (k = 0; k < x->length; k++) {
        h = &(x->buffer[k]);
		qn = Get_qn((*h)); tn = Get_tn((*h));
        if(h->bl < sl->ul_occ) continue;
        if(g->seq[tn].del) {
            get_R_to_U(ridx, tn, &rr, &is_u);
            if(rr == (uint32_t)-1 || is_u == 1) continue;
        }
        r = ma_hit2arc(h, g->seq[qn].len, g->seq[tn].len, sl->max_hang, asm_opt.max_hang_rate, sl->min_ovlp, &t);
        if(r < 0) continue;
        idx = get_specific_overlap(&(src[tn]), tn, qn);
        z = &(src[tn].buffer[idx]);
        assert(z->bl == h->bl);
		r = ma_hit2arc(z, g->seq[tn].len, g->seq[qn].len, sl->max_hang, asm_opt.max_hang_rate, sl->min_ovlp, &t);
        if(r < 0) continue;
        h->del = 0; if(!(g->seq[tn].del)) z->del = 0;
        g->seq_vis[i] = 1;
    }
}

void recover_contain_g(asg_t *g, ma_hit_t_alloc *src, R_to_U* ruIndex, int64_t max_hang, int64_t min_ovlp, int64_t ul_occ)
{
    sset_aux s; s.g = g; s.src = src; s.ruIndex = ruIndex;
    s.max_hang = max_hang; s.min_ovlp = min_ovlp; s.ul_occ = ul_occ;
    kt_for(asm_opt.thread_num, update_sg_contain, &s, g->n_seq);
    uint32_t k;
    for (k = 0; k < g->n_seq; k++) {
        if(g->seq_vis[k]) g->seq[k].del = 0;
    }
    memset(g->seq_vis, 0, (sizeof(*(g->seq_vis))*(g->n_seq<<1)));
}

// static void update_norm_arc(void *data, long i, int tid)
// {
//     sset_aux *sl = (sset_aux *)data; ma_hit_t *h, *z; asg_arc_t t;
// 	ma_hit_t_alloc *src = sl->src; int64_t r, idx; 
//     uint32_t k, rr, qn, tn, is_u;
//     ma_hit_t_alloc *x = &(src[i]);
//     for (k = 0; k < x->length; k++) {
//         h = &(x->buffer[k]);
// 		qn = Get_qn((*h)); tn = Get_tn((*h));
//         if(h->bl < sl->ul_occ) continue;
//         if(g->seq[tn].del) {
//             get_R_to_U(ridx, tn, &rr, &is_u);
//             if(rr == (uint32_t)-1 || is_u == 1) continue;
//         }
//         r = ma_hit2arc(h, g->seq[qn].len, g->seq[tn].len, sl->max_hang, asm_opt.max_hang_rate, sl->min_ovlp, &t);
//         if(r < 0) continue;
//         idx = get_specific_overlap(&(src[tn]), tn, qn);
//         z = &(src[tn].buffer[idx]);
//         assert(z->bl == h->bl);
// 		r = ma_hit2arc(z, g->seq[tn].len, g->seq[qn].len, sl->max_hang, asm_opt.max_hang_rate, sl->min_ovlp, &t);
//         if(r < 0) continue;
//         h->del = 0; if(!(g->seq[tn].del)) z->del = 0;
//         g->seq_vis[i] = 1;
//     }
// }

// void normalize_ma_hit_t_mul(ma_hit_t_alloc *src, uint32_t n_src)
// {
//     sset_aux s; s.src = src;
//     kt_for(asm_opt.thread_num, update_norm_arc, &s, n_src);
// }

static void normalize_gou0(void *data, long i, int tid)
{
    sset_aux *sl = (sset_aux *)data;
    asg_t *g = sl->g;
    asg_arc_t *e = &(g->arc[i]);
    if(e->v > (e->ul>>32)) return;
    uint32_t k, v = e->v^1, w = (e->ul>>32)^1, ou;
    uint32_t nv = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v);
    for (k = 0; k < nv; ++k) {
        if (av[k].v == w) {
            ou = MAX(av[k].ou, e->ou);
            av[k].ou = e->ou = ou;
            break;
        }
    }
}

void normalize_gou(asg_t *g)
{
    sset_aux s; s.g = g; 
    kt_for(asm_opt.thread_num, normalize_gou0, &s, g->n_arc);
}

static void update_sg_uo_t(void *data, long i, int tid)
{
	sset_aux *sl = (sset_aux *)data;
	ma_hit_t_alloc *src = sl->src; asg_t *g = sl->g;
	asg_arc_t *e = &(g->arc[i]); uint32_t k, qn, tn;
    ma_hit_t_alloc *x = &(src[e->ul>>33]);
    e->ou = 0;
    if(e->del) return;
	for (k = 0; k < x->length; k++) {
		qn = Get_qn(x->buffer[k]);
		tn = Get_tn(x->buffer[k]);
		if(qn == (e->ul>>33) && tn == (e->v>>1)) {
			e->ou = (x->buffer[k].bl>OU_MASK?OU_MASK:x->buffer[k].bl);
			break;
		}
	}
	assert(k < x->length);
}

void update_sg_uo(asg_t *g, ma_hit_t_alloc *src)
{
	sset_aux s; s.g = g; s.src = src;
	kt_for(asm_opt.thread_num, update_sg_uo_t, &s, g->n_arc);
    uint32_t k, z, nv, occ_a = 0, occ_n = 0; asg_arc_t *av = NULL;
    for (k = 0; k < g->n_seq; k++) {
        if(g->seq[k].del) continue;
        occ_n++;

        av = asg_arc_a(g, (k<<1)); nv = asg_arc_n(g, (k<<1));
        for (z = 0; z < nv; z++) {
            if(av[z].del || av[z].ou == 0) continue;
            break;
        }
        if(z < nv) {
            occ_a++;
            continue;
        }

        av = asg_arc_a(g, ((k<<1)+1)); nv = asg_arc_n(g, ((k<<1)+1));
        for (z = 0; z < nv; z++) {
            if(av[z].del || av[z].ou == 0) continue;
            break;
        }
        if(z < nv) {
            occ_a++;
        }
    }

    fprintf(stderr, "[M::%s::] ==> # gfa reads:%u, # covered gfa reads:%u\n", __func__, occ_n, occ_a);

    // asg_arc_t *e; uint32_t v, w;
    // for (k = 0; k < g->n_arc; k++) {
    //     e = &(g->arc[k]);
    //     v = e->v^1; w = (e->ul>>32)^1;
    //     av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
    //     for (z = 0; z < nv; z++) {
    //         if(av[z].v == w) break;
    //     }
    //     if(z >= nv || av[z].ou != e->ou) fprintf(stderr, "[M::%s::asymmetry]\n", __func__);
    // }
}

int32_t if_sup_chimeric(ma_hit_t_alloc* src, uint64_t rLen, asg64_v *b, int if_exact)
{
    uint32_t k, qs, qe, l[2], r[2], st, bn;
    int32_t dp, op;
    l[0] = r[0] = rLen; l[1] = r[1] = 0;
    for (k = 0; k < src->length; k++){
        if(src->buffer[k].del) continue;
        if(if_exact && !(src->buffer[k].el)) continue;

        qs = Get_qs(src->buffer[k]); qe = Get_qe(src->buffer[k]);

        ///overlaps from left side
        if(qs == 0){
            if(qs < l[0]) l[0] = qs;
            if(qe > l[1]) l[1] = qe;
        }

        ///overlaps from right side
        if(qe == rLen){
            if(qs < r[0]) r[0] = qs;
            if(qe > r[1]) r[1] = qe;
        }

        ///note: if (qs == 0 && qe == rLen)
        ///this overlap would be added to both b_left and b_right
        ///that is what we want
    }
    if (l[1] > r[0]) return 0;
    if (l[1] <= l[0] || r[1] <= r[0]) return 1;

    bn = b->n;
    if(l[1] > l[0]) {
        kv_push(uint64_t, *b, (l[0]<<1)); kv_push(uint64_t, *b, (l[1]<<1)|1);
    }
    if(r[1] > r[0]) {
        kv_push(uint64_t, *b, (r[0]<<1)); kv_push(uint64_t, *b, (r[1]<<1)|1);
    }

    ///check contained overlaps
    for (k = 0; k < src->length; k++) {
        if(src->buffer[k].del) continue;
        if(if_exact && !(src->buffer[k].el)) continue;

        qs = Get_qs(src->buffer[k]); qe = Get_qe(src->buffer[k]);
        if(qs == 0 || qe == rLen) continue;

        kv_push(uint64_t, *b, (qs<<1)); kv_push(uint64_t, *b, (qe<<1)|1);
    }
    radix_sort_srt64(b->a + bn, b->a + b->n);
    l[0] = r[0] = rLen; l[1] = r[1] = 0;

    for (k = bn, dp = st = 0; k < b->n; k++) {
        op = dp;
        ///if a[j] is qe
        if (b->a[k]&1) --dp;
        else ++dp;

        if(op < 1 && dp >= 1) {
            st = b->a[k]>>1;
        } else if(op >= 1 && dp < 1) {
            if(st == 0) l[0] = st, l[1] = b->a[k]>>1;
            if((b->a[k]>>1) == rLen) r[0] = st, r[1] = b->a[k]>>1;
        }
    }
    
    b->n = bn;
    if (l[1] > r[0]) return 0;

    return 1;
}

///remove single node
void asg_arc_cut_chimeric(asg_t *g, ma_hit_t_alloc* src, asg64_v *in, uint32_t ou_thres, telo_end_t *te)
{
	asg64_v tx = {0,0,0}, *b = NULL;
	uint32_t v, w, ei[2] = {0}, k, i, n_vtx = g->n_seq<<1;
    uint32_t nw, el_n, cnt = 0; asg_arc_t *aw;
	if(in) b = in;
    else b = &tx;
	b->n = 0;

	for (v = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;
        if (te && te->hh[v>>1]) continue;
		if(g->seq_vis[v] == 0) {
            if((get_arcs(g, v, &(ei[0]), 1)!=1) || (get_arcs(g, v^1, &(ei[1]), 1)!=1)) continue;
            assert((g->arc[ei[0]].ul>>32) == v && (g->arc[ei[1]].ul>>32) == (v^1));
            if((get_arcs(g, g->arc[ei[0]].v^1, NULL, 0)<2) || (get_arcs(g, g->arc[ei[1]].v^1, NULL, 0)<2)) continue;
            if(g->arc[ei[0]].el) continue;       
            if(ou_thres!=(uint32_t)-1&&g->arc[ei[0]].ou>=ou_thres&&g->arc[ei[1]].ou>=ou_thres) continue;///UL
            if(!if_sup_chimeric(&(src[v>>1]), g->seq[v>>1].len, b, 1)) continue;///HiFi
            kv_push(uint64_t, *b, (((uint64_t)(g->arc[ei[0]].ol))<<32)|((uint64_t)(ei[0])));
        }
	}

    radix_sort_srt64(b->a, b->a + b->n);
    ///here all edges are inexact matches 
    for (k = 0; k < b->n; k++) {
        if(g->arc[(uint32_t)b->a[k]].del) continue;
        v = g->arc[(uint32_t)b->a[k]].ul>>32; w = g->arc[(uint32_t)b->a[k]].v^1;
        if(g->seq[v>>1].del || g->seq[w>>1].del) continue;
        aw = asg_arc_a(g, w); nw = asg_arc_n(g, w);
        if((get_arcs(g, v, &(ei[0]), 1)!=1) || (get_arcs(g, v^1, &(ei[1]), 1)!=1)) continue;
        if((get_arcs(g, g->arc[ei[0]].v^1, NULL, 0)<2) || (get_arcs(g, g->arc[ei[1]].v^1, NULL, 0)<2)) continue;

        for (i = el_n = 0; i < nw; i++) {
            if ((aw[i].del) || (aw[i].v==(v^1)) || (!aw[i].el)) continue;
            el_n++; break;
        }

        if(!el_n) continue;
        if(te && te->hh[v>>1]) continue;
        asg_seq_del(g, v>>1);
        cnt++;
    }
    // stats_sysm(g);
	if(!in) free(tx.a);
    if (cnt > 0) asg_cleanup(g);
}

void asg_arc_cut_inexact(asg_t *g, ma_hit_t_alloc* src, asg64_v *in, int32_t max_ext, uint32_t is_ou, uint32_t is_trio, uint32_t min_diff, float ou_rat/**, asg64_v *dbg**/)
{
    asg64_v tx = {0,0,0}, *b = NULL;
    uint32_t v, w, i, k, n_vtx = g->n_seq<<1;
    asg_arc_t *av, *aw, *ve, *vmax, *we; uint32_t nv, nw, kv, kw, ol_max, ou_max, to_del, cnt = 0, mm_ol, mm_ou;
    uint32_t trioF = (uint32_t)-1, ntrioF = (uint32_t)-1;
    if(in) b = in;
    else b = &tx;
	b->n = 0;

    for (v = 0; v < n_vtx; ++v) {
        if(g->seq[v>>1].del) continue;
        if(g->seq_vis[v] == 0) {
            av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
            if (nv < 2) continue;

            for (i = kv = 0; i < nv; ++i) {
                if(av[i].del) continue;
                kv++;
            }
            if(kv < 2) continue;

            for (i = 0; i < nv; ++i) {
                if(av[i].del || av[i].el) continue;
                kv_push(uint64_t, *b, (uint64_t)((((uint64_t)av[i].ol)<<32)|((uint64_t)(av-g->arc+i))));   
            }
        }
	}

    radix_sort_srt64(b->a, b->a + b->n);
	for (k = 0; k < b->n; k++) {
        if(g->arc[(uint32_t)b->a[k]].del) continue;
        assert((!g->arc[(uint32_t)b->a[k]].el));

        v = g->arc[(uint32_t)b->a[k]].ul>>32; w = g->arc[(uint32_t)b->a[k]].v^1;
        if(g->seq[v>>1].del || g->seq[w>>1].del) continue;
        nv = asg_arc_n(g, v); nw = asg_arc_n(g, w);
        av = asg_arc_a(g, v); aw = asg_arc_a(g, w);
        // if(((v>>1) == 50356 && (w>>1) == 1276292)||((v>>1) == 1276292 && (w>>1) == 50356)) {
        //     fprintf(stderr, "[0]v->%u, w->%u, nv->%u, nw->%u\n", v, w, nv, nw);
        // }
        if(nv<=1 && nw <= 1) continue;
        if(is_trio) {
            if(get_arcs(g, v, NULL, 0)<=1 && get_arcs(g, w, NULL, 0)<=1) continue;///speedup
            trioF = get_tip_trio_infor(g, v^1);
            ntrioF = (trioF==FATHER? MOTHER : (trioF==MOTHER? FATHER : (uint32_t)-1));
        }
        ve = &(g->arc[(uint32_t)b->a[k]]);
		for (i = 0; i < nw; ++i) {
			if (aw[i].v == (v^1)) {
                we = &(aw[i]);
                break;
            }
        }
        ///mm_ol and mm_ou are used to make edge with long indel more easy to be cutted
        mm_ol = MIN(ve->ol, we->ol); mm_ou = MIN(ve->ou, we->ou);
        for (i = kv = ol_max = ou_max = 0, vmax = NULL; i < nv; ++i) {
            if(av[i].del) continue;
            kv++;
            if(is_trio && get_tip_trio_infor(g, av[i].v) == ntrioF) continue;
            if(ol_max < av[i].ol) ol_max = av[i].ol, vmax = &(av[i]);
            if(ou_max < av[i].ou) ou_max = av[i].ou;
        }
        // if(((v>>1) == 50356 && (w>>1) == 1276292)||((v>>1) == 1276292 && (w>>1) == 50356)) {
        //     fprintf(stderr, "[1]v->%u, w->%u, kv->%u, ve->ol->%u, ol_max->%u\n", v, w, kv, ve->ol, ol_max);
        // }
        if (kv < 1) continue;
        if (kv >= 2) {
            if (mm_ol >= ol_max) continue;
            if (is_ou && mm_ou > ou_max*ou_rat) continue;
            if ((mm_ol + min_diff) > ol_max) continue;
        }
        
        for (i = kw = ol_max = ou_max = 0; i < nw; ++i) {
            if(aw[i].del) continue;
            kw++;
            if(is_trio && get_tip_trio_infor(g, aw[i].v) == ntrioF) continue;
            if(ol_max < aw[i].ol) ol_max = aw[i].ol;
            if(ou_max < aw[i].ou) ou_max = aw[i].ou;
        }
        // if(((v>>1) == 50356 && (w>>1) == 1276292)||((v>>1) == 1276292 && (w>>1) == 50356)) {
        //     fprintf(stderr, "[1]v->%u, w->%u, kw->%u, we->ol->%u, ol_max->%u\n", v, w, kw, we->ol, ol_max);
        // }
        if (kw < 1) continue;
        if (kw >= 2) {
            if (mm_ol >= ol_max) continue;
            if (is_ou && mm_ou > ou_max*ou_rat) continue;
            if ((mm_ol + min_diff) > ol_max) continue;
        }
        if (kv <= 1 && kw <= 1) continue;

        to_del = 0;
        ///if there is an inexact edge between two good reads
        if(src[v>>1].is_fully_corrected == 1 && src[w>>1].is_fully_corrected == 1) {
            if (kv > 1 && kw > 1) {
				to_del = 1;
            } else if (kw == 1) {
                if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
            } else if (kv == 1) {
                if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
            }
        }
        ///TODO: should check if the edge wmax also works
        if(src[v>>1].is_fully_corrected == 1 && src[w>>1].is_fully_corrected == 0) {
            if(vmax && vmax->v != ve->v && vmax->el == 1 && src[vmax->v>>1].is_fully_corrected == 1) {
                if (kv > 1 && kw > 1) {
                    to_del = 1;
                } else if (kw == 1) {
                    if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
                } else if (kv == 1) {
                    if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
                }
            }
        }

        // if(((v>>1) == 50356 && (w>>1) == 1276292)||((v>>1) == 1276292 && (w>>1) == 50356)) {
        //     fprintf(stderr, "[3]v->%u, w->%u, to_del->%u, el->%u, src[v>>1].is_fully_corrected->%u, src[w>>1].is_fully_corrected->%u\n", 
        //             v, w, to_del, g->arc[(uint32_t)b->a[k]].el, src[v>>1].is_fully_corrected, src[w>>1].is_fully_corrected);
        //     fprintf(stderr, "[4]v->%u, w->%u, to_del->%u, vmax->v->%u, vmax->el->%u, src[vmax->v>>1].is_fully_corrected->%u\n", 
        //             v, w, to_del, vmax->v, vmax->el, src[vmax->v>>1].is_fully_corrected);
        // }

        if (to_del) {
            ve->del = we->del = 1, ++cnt;
            /**
            if(dbg) {
                kv_push(uint64_t, *dbg, ve - g->arc); 
                kv_push(uint64_t, *dbg, we - g->arc); 
                // if(((v>>1) == 50356 && (w>>1) == 1276292)||((v>>1) == 1276292 && (w>>1) == 50356)) {
                //     fprintf(stderr, "[5]v->%u, w->%u, (ve - g->arc)->%u, (we - g->arc)->%u\n", 
                //     v, w, (uint32_t)(ve - g->arc), (uint32_t)(we - g->arc));
                // }
            }
            **/
        }
        
    }
    // stats_sysm(g);
    if(!in) free(tx.a);
    if (cnt > 0) asg_cleanup(g);
}

void asg_arc_cut_inexact_debug(asg_t *g, ma_hit_t_alloc* src, asg64_v *in, int32_t max_ext, uint32_t is_ou, uint32_t is_trio, asg64_v *dbg)
{
    asg64_v tx = {0,0,0}, *b = NULL;
    uint32_t v, w, i, k, nv, nw, kv, kw, iv, iw, n_vtx = g->n_seq<<1, to_del, cnt = 0, ov_max = 0, ow_max = 0, ov_max_i = 0;
    asg_arc_t *av = NULL, *aw = NULL, *a = NULL;
    if(in) b = in;
    else b = &tx;
    b->n = 0;

    for (v = 0; v < n_vtx; ++v)  {
        if(g->seq[v>>1].del) continue;
        if(g->seq_vis[v] == 0) {
            av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
            if (nv < 2) continue;
            for (i = kv = 0; i < nv; ++i) {
                if(av[i].del) continue;
                kv++;
            }
            if(kv < 2) continue;

            for (i = 0; i < nv; ++i) {
                if(av[i].del || av[i].el) continue;
                kv_push(uint64_t, *b, (uint64_t)((uint64_t)av[i].ol << 32 | (av - g->arc + i)));   
            }
        }
	}
    radix_sort_srt64(b->a, b->a + b->n);

    for (k = 0; k < b->n; k++) {
        a = &g->arc[(uint32_t)b->a[k]];
        if(a->del) continue;
        v = (a->ul)>>32, w = a->v^1; to_del = 0;
        nv = asg_arc_n(g, v), nw = asg_arc_n(g, w);
        if(((v>>1) == 50356 && (w>>1) == 1276292)||((v>>1) == 1276292 && (w>>1) == 50356)) {
            fprintf(stderr, "[0]v->%u, w->%u, nv->%u, nw->%u\n", v, w, nv, nw);
        }
        if (nv == 1 && nw == 1) continue;
        av = asg_arc_a(g, v), aw = asg_arc_a(g, w);
        ov_max = ow_max = ov_max_i = 0;

        for (i = 0, kv = 0; i < nv; ++i) {
			if (av[i].del) continue;
			if (ov_max < av[i].ol) {
                ov_max = av[i].ol;
                ov_max_i = i;
            }
			++kv;
		}
        if(((v>>1) == 50356 && (w>>1) == 1276292)||((v>>1) == 1276292 && (w>>1) == 50356)) {
            fprintf(stderr, "[1]v->%u, w->%u, kv->%u, a->ol->%u, ov_max->%u\n", v, w, kv, a->ol, ov_max);
        }
		if (kv >= 2 && a->ol == ov_max) continue;

        for (i = 0, kw = 0; i < nw; ++i) {
			if (aw[i].del) continue;
			if (ow_max < aw[i].ol) {
                ow_max = aw[i].ol;
            }
			++kw;
		}
        if(((v>>1) == 50356 && (w>>1) == 1276292)||((v>>1) == 1276292 && (w>>1) == 50356)) {
            fprintf(stderr, "[2]v->%u, w->%u, kw->%u, a->ol->%u, ow_max->%u\n", v, w, kw, a->ol, ow_max);
        }
		if (kw >= 2 && a->ol == ow_max) continue;

        if (kv <= 1 && kw <= 1) continue;

        ///to see which one is the current edge (from v and w)
		for (iv = 0; iv < nv; ++iv)
			if (av[iv].v == (w^1)) break;
		for (iw = 0; iw < nw; ++iw)
			if (aw[iw].v == (v^1)) break;
        ///if one edge has been deleted, it should be deleted in both direction
        if (av[iv].del && aw[iw].del) continue;

        ///if this edge is an inexact edge
        if(a->el == 0 && src[v>>1].is_fully_corrected == 1 && src[w>>1].is_fully_corrected == 1) {
            if (kv > 1 && kw > 1) {
				to_del = 1;
            } else if (kw == 1) {
                if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
            } else if (kv == 1) {
                if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
            }
        }

        if(a->el == 0 && src[v>>1].is_fully_corrected == 1 && src[w>>1].is_fully_corrected == 0) {
            if(av[ov_max_i].el == 1 && src[av[ov_max_i].v>>1].is_fully_corrected == 1) {
                if (kv > 1 && kw > 1) {
                    to_del = 1;
                } else if (kw == 1) {
                    if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
                } else if (kv == 1) {
                    if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
                }
            }
        }

        if(((v>>1) == 50356 && (w>>1) == 1276292)||((v>>1) == 1276292 && (w>>1) == 50356)) {
            fprintf(stderr, "[3]v->%u, w->%u, to_del->%u, el->%u, src[v>>1].is_fully_corrected->%u, src[w>>1].is_fully_corrected->%u\n", 
                    v, w, to_del, a->el, src[v>>1].is_fully_corrected, src[w>>1].is_fully_corrected);
            fprintf(stderr, "[4]v->%u, w->%u, to_del->%u, av[ov_max_i].v->%u, av[ov_max_i].el->%u, src[av[ov_max_i].v>>1].is_fully_corrected->%u\n", 
                    v, w, to_del, av[ov_max_i].v, av[ov_max_i].el, src[av[ov_max_i].v>>1].is_fully_corrected);
        }

        if (to_del) {
            av[iv].del = aw[iw].del = 1, ++cnt;
            if(dbg) {
                kv_push(uint64_t, *dbg, (av - g->arc + iv)); 
                kv_push(uint64_t, *dbg, (aw - g->arc + iw)); 
                if(((v>>1) == 50356 && (w>>1) == 1276292)||((v>>1) == 1276292 && (w>>1) == 50356)) {
                    fprintf(stderr, "[5]v->%u, w->%u, (av-g->arc+iv)->%u, (aw-g->arc+iw)->%u\n", 
                    v, w, (uint32_t)(av - g->arc + iv), (uint32_t)(aw - g->arc + iw));
                }
            }
        }
    }

    if(!in) free(tx.a);
    if(dbg) {
        for (i = 0; i < dbg->n; i++) g->arc[dbg->a[i]].del = 0;
    }
    // if (cnt > 0) asg_cleanup(g);
}

uint32_t trans_path_check(uint32_t a, uint32_t b, asg_t *g, ma_hit_t_alloc *rev, R_to_U* rI, 
uint32_t minLen, asg64_v *t)
{
    if(a == b) return -1;
    uint32_t tn = t->n, m, e, l[2] = {0};
    uint64_t *x[2];
    m = follow_limit_path(g, a, &e, &(l[0]), t, (uint32_t)-1);
    if(m == LOOP || l[0] <= minLen) {
        t->n = tn; return -1;
    }

    m = follow_limit_path(g, b, &e, &(l[1]), t, (uint32_t)-1);
    if(m == LOOP || l[1] <= minLen) {
        t->n = tn; return -1;
    }

    x[0] = t->a + tn; x[1] = t->a + tn + l[0];
    if(l[0] > l[1]) {
        x[0] = t->a + tn + l[0]; x[1] = t->a + tn;
        m = l[0]; l[0] = l[1]; l[1] = m;
    }
    assert(l[0]+l[1]+tn==t->n);

    uint32_t i, k, qi, ti, isU; double max_count = 0, min_count = 0;
    for (i = 0; i < l[1]; i++) g->seq_vis[x[1][i]>>1] = 1;
    for (i = 0; i < l[0]; i++) {
        qi = x[0][i]>>1;
        for (k = 0; k < rev[qi].length; k++) {
            ti = Get_tn(rev[qi].buffer[k]);
            if(g->seq[ti].del == 1) {
                get_R_to_U(rI, ti, &ti, &isU);
                if(ti == (uint32_t)-1 || isU == 1 || g->seq[ti].del == 1) continue;
            }
            min_count++; max_count += g->seq_vis[ti];
        }
    }

    for (i = 0; i < l[1]; i++) g->seq_vis[x[1][i]>>1] = 0;
    t->n = tn;

    if(min_count == 0) return -1;
    if(max_count == 0) return 0;
    if((max_count/min_count)>0.3) return 1;
    return 0;
}

void asg_arc_cut_length(asg_t *g, asg64_v *in, int32_t max_ext, float len_rat, float ou_rat, uint32_t is_ou, uint32_t is_trio, 
uint32_t is_topo, uint32_t min_diff, uint32_t min_ou, ma_hit_t_alloc *rev, R_to_U* rI, uint32_t *max_drop_len)
{
    asg64_v tx = {0,0,0}, *b = NULL;
    uint32_t i, k, v, w, n_vtx = g->n_seq<<1, nv, nw, kv, kw, trioF = (uint32_t)-1, ntrioF = (uint32_t)-1, ol_max, ou_max, to_del, cnt = 0, mm_ol, mm_ou;
    asg_arc_t *av, *aw, *ve, *we, *vl_max, *wl_max;

    if(in) b = in;
    else b = &tx;
	b->n = 0;

    for (v = 0; v < n_vtx; ++v) {
        // if((v>>1)==17078) fprintf(stderr, "[M::%s::] v:%u, del:%u, seq_vis:%u\n", __func__, v, g->seq[v>>1].del, g->seq_vis[v]);
        if (g->seq[v>>1].del) continue;
        if(g->seq_vis[v] == 0) {
            av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
            if (nv < 2) continue;

            for (i = kv = 0; i < nv; ++i) {
                if(av[i].del) continue;
                kv++;
            }
            if(kv < 2) continue;

            for (i = 0; i < nv; ++i) {
                if(av[i].del) continue;
                // if((av[i].ul>>33)==287) {
                //     fprintf(stderr, "++++++%.*s(%c)\t%.*s(%c)\tol:%u\tou:%u\n", 
                //     (int32_t)Get_NAME_LENGTH(R_INF, (av[i].ul>>33)), Get_NAME(R_INF, (av[i].ul>>33)), "+-"[(av[i].ul>>32)&1], 
                //     (int32_t)Get_NAME_LENGTH(R_INF, (av[i].v>>1)), Get_NAME(R_INF, (av[i].v>>1)), "+-"[av[i].v&1], av[i].ol, av[i].ou);
                // }
                if(max_drop_len && av[i].ol >= (*max_drop_len)) continue;
                kv_push(uint64_t, *b, (((uint64_t)av[i].ol)<<32) | ((uint64_t)(av-g->arc+i)));   
            }
        }
	}

    if(rev && rI) memset(g->seq_vis, 0, g->n_seq*2*sizeof(uint8_t));
    radix_sort_srt64(b->a, b->a + b->n);
    for (k = 0; k < b->n; k++) {
        if(g->arc[(uint32_t)b->a[k]].del) continue;

        v = g->arc[(uint32_t)b->a[k]].ul>>32; w = g->arc[(uint32_t)b->a[k]].v^1;
        if(g->seq[v>>1].del || g->seq[w>>1].del) continue;
        nv = asg_arc_n(g, v); nw = asg_arc_n(g, w);
        av = asg_arc_a(g, v); aw = asg_arc_a(g, w);
        if(nv<=1 && nw <= 1) continue;

        if(is_trio) {
            if(get_arcs(g, v, NULL, 0)<=1 && get_arcs(g, w, NULL, 0)<=1) continue;///speedup
            trioF = get_tip_trio_infor(g, v^1);
            ntrioF = (trioF==FATHER? MOTHER : (trioF==MOTHER? FATHER : (uint32_t)-1));
        }

        ve = &(g->arc[(uint32_t)b->a[k]]);
        for (i = 0; i < nw; ++i) {
            if (aw[i].v == (v^1)) {
                we = &(aw[i]);
                break;
            }
        }
        ///mm_ol and mm_ou are used to make edge with long indel more easy to be cutted
        mm_ol = MIN(ve->ol, we->ol); mm_ou = MIN(ve->ou, we->ou);

        for (i = kv = ol_max = ou_max = 0, /**ve =**/ vl_max = NULL; i < nv; ++i) {
            if(av[i].del) continue;
            // if(av[i].v == (w^1)) ve = &(av[i]);
            // if((av[i].ul>>33)==287) {
            //     fprintf(stderr, "++++++%.*s(%c)\t%.*s(%c)\tol:%u\tou:%u\n", 
            //     (int32_t)Get_NAME_LENGTH(R_INF, (av[i].ul>>33)), Get_NAME(R_INF, (av[i].ul>>33)), "+-"[(av[i].ul>>32)&1], 
            //     (int32_t)Get_NAME_LENGTH(R_INF, (av[i].v>>1)), Get_NAME(R_INF, (av[i].v>>1)), "+-"[av[i].v&1], av[i].ol, av[i].ou);
            // }
            kv++;
            if(is_trio && get_tip_trio_infor(g, av[i].v) == ntrioF) continue;
            if(ol_max < av[i].ol) ol_max = av[i].ol, vl_max = &(av[i]);
            if(ou_max < av[i].ou) ou_max = av[i].ou;
        }
        if (kv < 1) continue;
        if (kv >= 2) {
            if (mm_ol > ol_max*len_rat) continue;
            if (is_ou && mm_ou > ou_max*ou_rat && mm_ou > min_ou) continue;
            if ((mm_ol + min_diff) > ol_max) continue;
        }
        

        for (i = kw = ol_max = ou_max = 0, wl_max = NULL; i < nw; ++i) {
            if(aw[i].del) continue;
            kw++;
            if(is_trio && get_tip_trio_infor(g, aw[i].v) == ntrioF) continue;
            if(ol_max < aw[i].ol) ol_max = aw[i].ol, wl_max = &(aw[i]);
            if(ou_max < aw[i].ou) ou_max = aw[i].ou;
        }
        if (kw < 1) continue;
        if (kw >= 2) {
            if (mm_ol > ol_max*len_rat) continue;
            if (is_ou && mm_ou > ou_max*ou_rat && mm_ou > min_ou) continue;
            if ((mm_ol + min_diff) > ol_max) continue;
        }

        if (kv <= 1 && kw <= 1) continue;

        to_del = 0;
        if(is_topo) {
            if (kv > 1 && kw > 1) {
                to_del = 1;
            } else if (kw == 1) {
                if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
            } else if (kv == 1) {
                if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
            }
        }

        if(rev && rI) {
            if((to_del == 0) && vl_max && (ve->v!=vl_max->v) && (trans_path_check(ve->v, vl_max->v, g, rev, rI, max_ext, b)==0)) {
                to_del = 1;
            }
            if((to_del == 0) && wl_max && (we->v!=wl_max->v) && (trans_path_check(we->v, wl_max->v, g, rev, rI, max_ext, b)==0)) {
                to_del = 1;
            }
            if(vl_max && wl_max) assert(ve->v!=vl_max->v||we->v!=wl_max->v);
        }
        

        if (to_del) {
            ve->del = we->del = 1, ++cnt;
        }
    }
    // stats_sysm(g);
    if(!in) free(tx.a);
    if (cnt > 0) asg_cleanup(g);
}


void asg_arc_cut_length_adv(asg_t *g, asg64_v *in, int32_t max_ext, float len_rat, float ou_rat, uint32_t is_ou, uint32_t is_trio, 
uint32_t is_topo, uint32_t min_diff, ma_hit_t_alloc *rev, R_to_U* rI, uint32_t *max_drop_len)
{
    asg64_v tx = {0,0,0}, *b = NULL;
    uint32_t i, k, v, w, n_vtx = g->n_seq<<1, nv, nw, kv, kw, trioF = (uint32_t)-1, ntrioF = (uint32_t)-1, ol_max, ou_max, to_del, cnt = 0, mm_ol, mm_ou;
    asg_arc_t *av, *aw, *ve, *we, *vl_max, *wl_max;

    if(in) b = in;
    else b = &tx;
	b->n = 0;

    for (v = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;
        if(g->seq_vis[v] == 0) {
            av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
            if (nv < 2) continue;

            for (i = kv = 0; i < nv; ++i) {
                if(av[i].del) continue;
                kv++;
            }
            if(kv < 2) continue;

            for (i = 0; i < nv; ++i) {
                if(av[i].del) continue;
                if(max_drop_len && av[i].ol >= (*max_drop_len)) continue;
                kv_push(uint64_t, *b, (((uint64_t)av[i].ol)<<32) | ((uint64_t)(av-g->arc+i)));   
            }
        }
	}

    if(rev && rI) memset(g->seq_vis, 0, g->n_seq*2*sizeof(uint8_t));
    radix_sort_srt64(b->a, b->a + b->n);
    for (k = 0; k < b->n; k++) {
        if(g->arc[(uint32_t)b->a[k]].del) continue;

        v = g->arc[(uint32_t)b->a[k]].ul>>32; w = g->arc[(uint32_t)b->a[k]].v^1;
        if(g->seq[v>>1].del || g->seq[w>>1].del) continue;
        nv = asg_arc_n(g, v); nw = asg_arc_n(g, w);
        av = asg_arc_a(g, v); aw = asg_arc_a(g, w);
        if(nv<=1 && nw <= 1) continue;

        if(is_trio) {
            if(get_arcs(g, v, NULL, 0)<=1 && get_arcs(g, w, NULL, 0)<=1) continue;///speedup
            trioF = get_tip_trio_infor(g, v^1);
            ntrioF = (trioF==FATHER? MOTHER : (trioF==MOTHER? FATHER : (uint32_t)-1));
        }

        ve = &(g->arc[(uint32_t)b->a[k]]);
        for (i = 0; i < nw; ++i) {
            if (aw[i].v == (v^1)) {
                we = &(aw[i]);
                break;
            }
        }
        ///mm_ol and mm_ou are used to make edge with long indel more easy to be cutted
        mm_ol = MIN(ve->ol, we->ol); mm_ou = MIN(ve->ou, we->ou);

        for (i = kv = ol_max = ou_max = 0, /**ve =**/ vl_max = NULL; i < nv; ++i) {
            if(av[i].del) continue;
            kv++;
            if(is_trio && get_tip_trio_infor(g, av[i].v) == ntrioF) continue;
            if(ol_max < av[i].ol) ol_max = av[i].ol, vl_max = &(av[i]);
            if(ou_max < av[i].ou) ou_max = av[i].ou;
        }
        if (kv < 1) continue;
        if (kv >= 2) {
            if (mm_ol > ol_max*len_rat) continue;
            if (is_ou && mm_ou > ou_max*ou_rat) continue;
            if ((mm_ol + min_diff) > ol_max) continue;
        }
        

        for (i = kw = ol_max = ou_max = 0, wl_max = NULL; i < nw; ++i) {
            if(aw[i].del) continue;
            kw++;
            if(is_trio && get_tip_trio_infor(g, aw[i].v) == ntrioF) continue;
            if(ol_max < aw[i].ol) ol_max = aw[i].ol, wl_max = &(aw[i]);
            if(ou_max < aw[i].ou) ou_max = aw[i].ou;
        }
        if (kw < 1) continue;
        if (kw >= 2) {
            if (mm_ol > ol_max*len_rat) continue;
            if (is_ou && mm_ou > ou_max*ou_rat) continue;
            if ((mm_ol + min_diff) > ol_max) continue;
        }

        if (kv <= 1 && kw <= 1) continue;

        to_del = 0;
        if(is_topo) {
            if (kv > 1 && kw > 1) {
                to_del = 1;
            } else if (kw == 1) {
                if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
            } else if (kv == 1) {
                if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
            }
        }

        if(rev && rI) {
            if((to_del == 0) && vl_max && (ve->v!=vl_max->v) && (trans_path_check(ve->v, vl_max->v, g, rev, rI, max_ext, b)==0)) {
                to_del = 1;
            }
            if((to_del == 0) && wl_max && (we->v!=wl_max->v) && (trans_path_check(we->v, wl_max->v, g, rev, rI, max_ext, b)==0)) {
                to_del = 1;
            }
            if(vl_max && wl_max) assert(ve->v!=vl_max->v||we->v!=wl_max->v);
        }
        

        if (to_del) {
            ve->del = we->del = 1, ++cnt;
        }
    }
    // stats_sysm(g);
    if(!in) free(tx.a);
    if (cnt > 0) asg_cleanup(g);
}

uint8_t get_ug_tip_trio_infor(ma_ug_t *ug, uint32_t begNode)
{
    uint32_t v = begNode, w, k, s;
    uint32_t kv;
    uint32_t eLen = 0, uLen = 0;
    uint32_t father_occ = 0, mother_occ = 0, ambigious_occ = 0;
    ma_utg_t *u;

    while (1) {
        kv = get_real_length(ug->g, v, NULL);
        u = &(ug->u.a[v>>1]); eLen += u->n; 

        for (k = 0; k < u->n; k++) {
            s = u->a[k]>>33;
            if(R_INF.trio_flag[s]==FATHER) father_occ++;
            else if(R_INF.trio_flag[s]==MOTHER) mother_occ++;
            else if((R_INF.trio_flag[s]==AMBIGU) || (R_INF.trio_flag[s]==DROP)) ambigious_occ++;
        }
        if(kv!=1) break;
        ///kv must be 1 here
        kv = get_real_length(ug->g, v, &w);
        if(get_real_length(ug->g, w^1, NULL)!=1) break;
        v = w;
        if(v == begNode) break;
    }

    uLen = eLen;///how many nodes
    eLen = father_occ + mother_occ;///haplotype-sepcific nodes
    if(eLen == 0) return AMBIGU;
    if(father_occ >= mother_occ) {
        if((father_occ > TRIO_THRES*eLen) && (father_occ >= DOUBLE_CHECK_THRES*uLen)) return FATHER;
    } else {
        if((mother_occ > TRIO_THRES*eLen) && (mother_occ >= DOUBLE_CHECK_THRES*uLen)) return MOTHER;
    }
    return AMBIGU;
}


asg_arc_t *iter_flex_asg(flex_asg_t *fg, flex_asg_e_retrive_t *rr, uint32_t v)
{
    asg_arc_t *av; uint32_t nv; asg_arc_t *z;
    av = asg_arc_a(fg->g, v); nv = asg_arc_n(fg->g, v);
    while(rr->i[0] < nv) {
        return &(av[rr->i[0]++]);
    }

    if(rr->i[0] >= nv && rr->i[0] != ((uint32_t)-1)) {
        rr->i[0] = ((uint32_t)-1); rr->i[1] = fg->idx[v];
    }
    while (rr->i[1] != ((uint32_t)-1)) {
        z = &(fg->a[rr->i[1]]); rr->i[1] = fg->pi.a[rr->i[1]];
        return z;
    }
    return NULL;
}

uint32_t detect_tip2(flex_asg_t *fg, uint32_t id, asg64_v *st, float ou_rat)
{
    uint32_t v, w, kw, ou_max, mm_ou; asg_arc_t *sv, *sw;
    flex_asg_e_retrive_t rv, rw;
    v = id<<1; rv.i[0] = 0; rv.i[1] = (uint32_t)-1;
    while(1) {
        sv = iter_flex_asg(fg, &rv, v);
        if(!sv) break;
        if(sv->del) continue;
        w = sv->v^1;
        if(fg->g->seq_vis[w]&128) continue;
        if(fg->g->seq_vis[w>>1]&2) continue;
        
        rw.i[0] = 0; rw.i[1] = (uint32_t)-1; kw = 0; ou_max = mm_ou = 0;
        while(1) {
            sw = iter_flex_asg(fg, &rw, w);
            if(!sw) break;
            if(sw->del) continue;
            if(ou_max < sw->ou) ou_max = sw->ou;
            if(fg->g->seq_vis[sw->v>>1]&2) {
                if(mm_ou < sw->ou) mm_ou = sw->ou;
                continue;
            }
            kw++;
        }
        if(kw < 1) return 1;
        if(ou_rat >= 0 && mm_ou > ou_max*ou_rat) return 1;
        fg->g->seq_vis[w] |= 128; kv_push(uint64_t, *st, w);
    }

    v = (id<<1)+1; rv.i[0] = 0; rv.i[1] = (uint32_t)-1;
    while(1) {
        sv = iter_flex_asg(fg, &rv, v);
        if(!sv) break;
        if(sv->del) continue;
        w = sv->v^1;
        if(fg->g->seq_vis[w]&128) continue;
        if(fg->g->seq_vis[w>>1]&2) continue;
        
        rw.i[0] = 0; rw.i[1] = (uint32_t)-1; kw = 0; ou_max = mm_ou = 0;
        while(1) {
            sw = iter_flex_asg(fg, &rw, w);
            if(!sw) break;
            if(sw->del) continue;
            if(ou_max < sw->ou) ou_max = sw->ou;
            if(fg->g->seq_vis[sw->v>>1]&2) {
                if(mm_ou < sw->ou) mm_ou = sw->ou;
                continue;
            }
            kw++;
        }
        if(kw < 1) return 1;
        if(ou_rat >= 0 && mm_ou > ou_max*ou_rat) return 1;
        fg->g->seq_vis[w] |= 128; kv_push(uint64_t, *st, w);
    }

    // av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
    // for (i = 0; i < nv; ++i) {
    //     if(av[i].del) continue;
    //     w = av[i].v^1;
    //     if(g->seq_vis[w]&128) continue;
    //     if((g->seq_vis[w>>1]&3)==2) continue;
    //     aw = asg_arc_a(g, w); nw = asg_arc_n(g, w);
    //     for (k = kw = 0; k < nw && kw < 1; k++) {
    //         if(aw[k].del) continue;
    //         if((g->seq_vis[aw[k].v>>1]&3)==2) continue;
    //         kw++;
    //     }
    //     if(kw < 1) return 1;
    //     g->seq_vis[w] |= 128; kv_push(uint64_t, *st, w);
    // }


    // v = (id<<1)+1;
    // av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
    // for (i = 0; i < nv; ++i) {
    //     if(av[i].del) continue;
    //     w = av[i].v^1;
    //     if(g->seq_vis[w]&128) continue;
    //     if((g->seq_vis[w>>1]&3)==2) continue;
    //     aw = asg_arc_a(g, w); nw = asg_arc_n(g, w);
    //     for (k = kw = 0; k < nw && kw < 1; k++) {
    //         if(aw[k].del) continue;
    //         if((g->seq_vis[aw[k].v>>1]&3)==2) continue;
    //         kw++;
    //     }
    //     if(kw < 1) return 1;
    //     g->seq_vis[w] |= 128; kv_push(uint64_t, *st, w);
    // }

    return 0;
}

uint32_t trans_check(flex_asg_t *fg, asg_arc_t *z, asg64_v *b)
{
    flex_asg_e_retrive_t rv, rw; asg_arc_t *sv, *sw; 
    uint32_t v = z->ul>>32, w;

    rv.i[0] = 0; rv.i[1] = (uint32_t)-1;
    while(1) {
        sv = iter_flex_asg(fg, &rv, v);
        if(!sv) break;
        if(sv->del) continue;
        if(sv->v == z->v) return 0;
        w = sv->v;
        // if((z->ul>>33) == 8340 && (z->v>>1) == 8352) {
        //     fprintf(stderr, "+[M::%s] v>>1::%u(%c), w>>1::%u(%c), sv->v::%u\n", 
        //     __func__, (z->ul>>33), "+-"[(z->ul>>32)&1], (z->v>>1), "+-"[z->v&1], sv->v);
        // }

        rw.i[0] = 0; rw.i[1] = (uint32_t)-1;
        while (1) {
            sw = iter_flex_asg(fg, &rw, w);
            if(!sw) break;
            if(sw->del) continue;
            if(sw->v == z->v) return 0;
        }
    }


    uint32_t bn = b->n; uint64_t vl, d, L = asg_arc_len((*z)) + fg->gap_fuzz;
    kv_push(uint64_t, *b, v);
    while (b->n > bn) {
        vl = kv_pop(*b); v = (uint32_t)vl; vl >>= 32;
        rv.i[0] = 0; rv.i[1] = (uint32_t)-1;
         while(1) {
            sv = iter_flex_asg(fg, &rv, v);
            if(!sv) break;
            if(sv->del) continue;
            d = vl + asg_arc_len((*sv));
            if(d > L) continue;
            if(sv->v == z->v) {
                b->n = bn;
                return 0;
            }
            d <<= 32; d |= sv->v;
            kv_push(uint64_t, *b, d); 
         }
    }

    b->n = bn;
    return 1;
}

void push_flex_asg_t(flex_asg_t *fg, asg_arc_t *z)
{
    asg_arc_t *av; uint32_t nv, k, v = z->ul>>32;
    av = asg_arc_a(fg->g, v); nv = asg_arc_n(fg->g, v);
    for (k = 0; k < nv; k++) {
        if(av[k].del) {
            av[k] = *z; 
            if(!(fg->need_srt)) {
                if(k > 0 && av[k].ul < av[k-1].ul) fg->need_srt = 1;
                if(k+1 < nv && av[k].ul > av[k+1].ul) fg->need_srt = 1;
            }
            return;
        }
    }

    k = fg->n; fg->need_srt = 1;
    kv_push(asg_arc_t, *fg, *z); 
    kv_push(uint32_t, fg->pi, fg->idx[v]); fg->idx[v] = k;
}

void append_notrans_e(flex_asg_t *fg, uint64_t *a, uint64_t a_n, asg64_v *b)
{
    ma_hit_t_alloc* src = fg->src; int32_t idx, r; 
    uint32_t i, k, v, w; ma_hit_t_alloc *z; asg_arc_t t0, t1; 
    for (i = 0; i < a_n; i++) {
        v = a[i]; 
        z = &(src[v>>1]);
        for (k = i + 1; k < a_n; k++) {
            w = a[k]^1;
            idx = get_specific_overlap(z, v>>1, w>>1);
            if(idx < 0 || z->buffer[idx].del) continue;
            r = ma_hit2arc(&(z->buffer[idx]), Get_READ_LENGTH(R_INF, v>>1), Get_READ_LENGTH(R_INF, w>>1),
						fg->max_hang, fg->max_hang_rate, fg->min_ovlp, &t0);
			if(r < 0) continue;
            if((t0.ul>>32) != v || t0.v != w) continue;
            t0.ou = ((z->buffer[idx].bl>OU_MASK)?OU_MASK:z->buffer[idx].bl);
            // if((v>>1) == 8340 && (w>>1) == 8352) {
            //     fprintf(stderr, "[M::%s] v>>1::%u(%c), w>>1::%u(%c), t0.ou::%u\n", 
            //     __func__, (v>>1), "+-"[v&1], (w>>1), "+-"[w&1], t0.ou);
            // }


            idx = get_specific_overlap(&(src[w>>1]), w>>1, v>>1);
            if(idx < 0 || src[w>>1].buffer[idx].del) continue;
            r = ma_hit2arc(&(src[w>>1].buffer[idx]), Get_READ_LENGTH(R_INF, w>>1), Get_READ_LENGTH(R_INF, v>>1),
						fg->max_hang, fg->max_hang_rate, fg->min_ovlp, &t1);
            if(r < 0) continue;
            if((t1.ul>>32) != (w^1) || t1.v != (v^1)) continue;
            t1.ou = ((src[w>>1].buffer[idx].bl>OU_MASK)?OU_MASK:src[w>>1].buffer[idx].bl);
            // if((v>>1) == 8340 && (w>>1) == 8352) {
            //     fprintf(stderr, "[M::%s] v>>1::%u(%c), w>>1::%u(%c), t1.ou::%u\n", 
            //     __func__, (v>>1), "+-"[v&1], (w>>1), "+-"[w&1], t1.ou);
            // }

            if(trans_check(fg, &t0, b) && trans_check(fg, &t1, b)) {
                push_flex_asg_t(fg, &t0); push_flex_asg_t(fg, &t1);
            }
        }
    }
}

uint32_t append_trans_check(flex_asg_t *fg, uint64_t *a, uint64_t a_n)
{
    ma_hit_t_alloc* src = fg->src; int32_t idx, r; 
    uint32_t i, k, v, w; ma_hit_t_alloc *z; asg_arc_t t0, t1; 
    for (i = 0; i < a_n; i++) {
        v = a[i]; 
        z = &(src[v>>1]);
        for (k = i + 1; k < a_n; k++) {
            w = a[k]^1;
            idx = get_specific_overlap(z, v>>1, w>>1);
            if(idx < 0 || z->buffer[idx].del) return 0;
            r = ma_hit2arc(&(z->buffer[idx]), Get_READ_LENGTH(R_INF, v>>1), Get_READ_LENGTH(R_INF, w>>1),
                        fg->max_hang, fg->max_hang_rate, fg->min_ovlp, &t0);
            if(r < 0) return 0;
            if((t0.ul>>32) != v || t0.v != w) return 0;
            t0.ou = ((z->buffer[idx].bl>OU_MASK)?OU_MASK:z->buffer[idx].bl);
            // if((v>>1) == 8340 && (w>>1) == 8352) {
            //     fprintf(stderr, "[M::%s] v>>1::%u(%c), w>>1::%u(%c), t0.ou::%u\n", 
            //     __func__, (v>>1), "+-"[v&1], (w>>1), "+-"[w&1], t0.ou);
            // }


            idx = get_specific_overlap(&(src[w>>1]), w>>1, v>>1);
            if(idx < 0 || src[w>>1].buffer[idx].del) return 0;
            r = ma_hit2arc(&(src[w>>1].buffer[idx]), Get_READ_LENGTH(R_INF, w>>1), Get_READ_LENGTH(R_INF, v>>1),
                        fg->max_hang, fg->max_hang_rate, fg->min_ovlp, &t1);
            if(r < 0) return 0;
            if((t1.ul>>32) != (w^1) || t1.v != (v^1)) return 0;
            t1.ou = ((src[w>>1].buffer[idx].bl>OU_MASK)?OU_MASK:src[w>>1].buffer[idx].bl);
        }
    }
    return 1;
}

uint32_t iter_contain_g(R_to_U* rI, flex_asg_t *fg, uint32_t v0, asg64_v *b, asg64_v *st, float ou_rat, uint32_t only_trans_nn)
{
    uint32_t v, w = (uint32_t)-1, x, i, kv, kw, ulen, cnt = 0, st_n, m, is_purge; 
    asg_arc_t *s; flex_asg_e_retrive_t rr;
    b->n = st->n = ulen = 0; kv_push(uint64_t, *st, v0);
    while (st->n) {
        v = kv_pop(*st); kv = kw = (uint32_t)-1; ulen = 0;
        if(fg->g->seq_vis[v>>1]) continue;
        while ((!(fg->g->seq_vis[v>>1])) && (is_contain_r((*rI), (v>>1)))) {///push a unitig
            kv = get_flex_arcs(fg, v, &w, 1); ulen++;
            kv_push(uint64_t, *b, v); 
            if(kv == 1) {
                flex_arcs0(s, (*fg), w);
                w = s->v;
                kw = get_flex_arcs(fg, w^1, NULL, 0);
                if(kw == 1) v = w;
                else break;
            } else {
                break;
            }
        }
        // if((v0>>1) == 12321 || (v0>>1) == 12334) {
        //     fprintf(stderr, "0[M::%s] v0>>1::%u(%c), b->n::%u, ulen::%u\n", 
        //     __func__, (v0>>1), "+-"[v0&1], (uint32_t)b->n, ulen);
        // }   
        if((!ulen) || (fg->g->seq_vis[v>>1]) || (!(is_contain_r((*rI), (v>>1))))) {
            for (i = b->n - ulen; i < b->n; i++) {
                fg->g->seq_vis[b->a[i]>>1] = 1; b->a[i] |= ((uint64_t)0x100000000);
            }
            continue;
        }
        // fprintf(stderr, "1[M::%s] v0>>1::%u(%c), b->n::%u, ulen::%u\n", 
        //     __func__, (v0>>1), "+-"[v0&1], (uint32_t)b->n, ulen);
        
        for (i = 0; i < b->n; i++) {
            if(b->a[i]&(0x100000000)) continue;
            fg->g->seq_vis[b->a[i]>>1] = 2;
        }
        // fprintf(stderr, "2[M::%s] v0>>1::%u(%c), b->n::%u, ulen::%u\n", 
        //     __func__, (v0>>1), "+-"[v0&1], (uint32_t)b->n, ulen);
        for (i = 0, st_n = st->n; i < b->n; i++) {
            if(b->a[i]&(0x100000000)) continue;
            if(detect_tip2(fg, b->a[i]>>1, st, ou_rat)) break;
        }
        if(i >= b->n) is_purge = 1;
        else is_purge = 0;
        // fprintf(stderr, "3[M::%s] v0>>1::%u(%c), b->n::%u, ulen::%u\n", 
        //     __func__, (v0>>1), "+-"[v0&1], (uint32_t)b->n, ulen);
        for (i = m = st_n; i < st->n; i++) {
            if(fg->g->seq_vis[st->a[i]]&128) fg->g->seq_vis[st->a[i]] -= 128;
            ///append edges to all nodes, instead of non-contained only
            // if(is_contain_r((*rI), (st->a[i]>>1))) continue;
            st->a[m++] = st->a[i];
        }
        st->n = m;
        // if((v0>>1) == 12321 || (v0>>1) == 12334) {
        //     fprintf(stderr, "4[M::%s] v0>>1::%u(%c), b->n::%u, i::%u\n", 
        //         __func__, (v0>>1), "+-"[v0&1], (uint32_t)b->n, i);
        // }
        if(is_purge && only_trans_nn) {
            is_purge = append_trans_check(fg, st->a + st_n, st->n-st_n);
        }
        
        if(is_purge) {
            for (i = 0; i < b->n; i++) {
                x = ((uint32_t)b->a[i])>>1;
                // if((v0>>1) == 12321 || (v0>>1) == 12334) {
                //     fprintf(stderr, "del::[M::%s] x::%u, seq_vis::%u\n", __func__, x, fg->g->seq_vis[x]);
                // }
                if(fg->g->seq_vis[x]&2) {
                    asg_seq_del(fg->g, x); cnt++;
                }
                fg->g->seq_vis[x] = 0;
            }
            append_notrans_e(fg, st->a + st_n, st->n-st_n, b);
            st->n = st_n;
            return cnt;
        }
        for (i = 0; i < b->n; i++) {
            if(b->a[i]&(0x100000000)) continue;
            fg->g->seq_vis[b->a[i]>>1] = 1;
        }
        st->n = st_n;

        rr.i[0] = 0; rr.i[1] = (uint32_t)-1;
        while(1) {
            s = iter_flex_asg(fg, &rr, v);
            if(!s) break;
            if(s->del) continue;
            if(fg->g->seq_vis[s->v>>1]) continue;
            if(!is_contain_r((*rI), (s->v>>1))) continue;
            kv_push(uint64_t, *st, s->v); 
        }
    }
    for (i = 0; i < b->n; i++) fg->g->seq_vis[((uint32_t)b->a[i])>>1] = 0;
    return cnt;
}

void flex_asg_t_cleanup(flex_asg_t *fg)
{
    asg_arc_t *p; uint32_t i;
    if(fg->n) {
        for (i = 0; i < fg->n; i++) {
            if(fg->a[i].del) continue;
            p = asg_arc_pushp(fg->g);
            *p = fg->a[i];
        }
        free(fg->g->idx);
        fg->g->idx = 0;
        fg->g->is_srt = 0;
    } else if(fg->need_srt){
        fg->g->is_srt = 0;
    }
    asg_cleanup(fg->g);
    // asg_symm(fg->g);
    // asg_arc_del_trans_ul(fg->g, fg->gap_fuzz);
}

void asg_arc_cut_contain(flex_asg_t *fg, asg64_v *in, asg64_v *in0, R_to_U* rI, float ou_rat, uint32_t only_trans_nn)
{
    // fprintf(stderr, "+[M::%s]\n", __func__);
    asg64_v tx = {0,0,0}, tx0 = {0,0,0}, *b = NULL, *b0 = NULL; 
    uint32_t v, w = (uint32_t)-1, n_vtx = fg->g->n_seq<<1, cnt = 0; asg_arc_t *s;
    b = in?in:&tx; b0 = in0?in0:&tx0; b->n = b0->n = 0; 
    fg->pi.n = fg->n = 0; memset(fg->idx, -1, (fg->g->n_seq<<1)*sizeof(*(fg->idx)));
    fg->need_srt = 0;

    memset(fg->g->seq_vis, 0, sizeof(*(fg->g->seq_vis))*n_vtx);
    for (v = 0; v < n_vtx; ++v) {
        // if((v>>1) == 6236) {
        //     fprintf(stderr, "[M::%s] v>>1::%u(%c), del::%u, contain::%u, get_arcs(v)::%u, get_arcs(v^1)::%u\n", 
        //     __func__, (v>>1), "+-"[v&1], g->seq[v>>1].del, is_contain_r((*rI), (v>>1)), 
        //     get_arcs(g, v, NULL, 0), get_arcs(g, v^1, NULL, 0));
        // }  
        // if((v>>1) == 749651) {
        //     fprintf(stderr, "+[M::%s] v>>1::%u(%c), del::%u, contain::%u, fg->n::%u, fg->need_srt::%u\n", 
        //     __func__, (v>>1), "+-"[v&1], fg->g->seq[v>>1].del, is_contain_r((*rI), (v>>1)), (uint32_t)fg->n, (uint32_t)fg->need_srt);
        // }
        if (fg->g->seq[v>>1].del) continue;
        if(!is_contain_r((*rI), (v>>1))) continue;
        if(get_flex_arcs(fg, v^1, &w, 1) == 1) {
            flex_arcs0(s, (*fg), w);
            if(get_flex_arcs(fg, s->v^1, NULL, 0) == 1) continue;
        }
        // if((v>>1) == 12321 || (v>>1) == 12334) {
        //     fprintf(stderr, "-[M::%s] v>>1::%u(%c), del::%u, contain::%u, fg->n::%u, fg->need_srt::%u\n", 
        //     __func__, (v>>1), "+-"[v&1], fg->g->seq[v>>1].del, is_contain_r((*rI), (v>>1)), (uint32_t)fg->n, (uint32_t)fg->need_srt);
        // }
        
        cnt += iter_contain_g(rI, fg, v, b, b0, ou_rat, only_trans_nn);
        // if((v>>1) == 12321 || (v>>1) == 12334) {
        //     fprintf(stderr, "*[M::%s] v>>1::%u(%c), del::%u, contain::%u, fg->n::%u, fg->need_srt::%u\n", 
        //     __func__, (v>>1), "+-"[v&1], fg->g->seq[v>>1].del, is_contain_r((*rI), (v>>1)), (uint32_t)fg->n, (uint32_t)fg->need_srt);
        // }
    }
    // stats_sysm(g);
    if(!in) free(tx.a); if(!in0) free(tx0.a);
    if(cnt > 0) flex_asg_t_cleanup(fg);
    // fprintf(stderr, "-[M::%s]\n", __func__);
}

void label_contain_dup(asg_t *g, R_to_U* rI, uint32_t v0, asg64_v *b, asg64_v *dump)
{
    asg_arc_t *av; uint32_t nv, v, i;
    if(!is_contain_r((*rI), (v0>>1))) return;
    b->n = 0; kv_push(uint64_t, *b, v0);
    while (b->n) {
        v = kv_pop(*b);
        if(g->seq_vis[v]&1) continue;
        kv_push(uint64_t, *dump, v);
        g->seq_vis[v] |= 1;
        av = asg_arc_a(g, v); 
        nv = asg_arc_n(g, v);
        for (i = 0; i < nv; ++i) {
            if(av[i].del || (g->seq_vis[av[i].v]&1) || (!is_contain_r((*rI), (av[i].v>>1)))) continue;
            kv_push(uint64_t, *b, av[i].v);
        }
    }
}

/**
void asg_arc_contain_trans_del(asg_t *g, asg64_v *in, asg64_v *in0, R_to_U* rI, float ou_rat)
{
    uint64_t n_vtx = g->n_seq<<1, i, k, v, w, nv, kv; asg_arc_t *av;
    memset(g->seq_vis, 0, sizeof((*g->seq_vis))*n_vtx);
    for (v = 0; v < n_vtx; ++v) {
        av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
        for (i = kv = 0; i < nv; i++) {
            if(av[i].del) continue;
            if(is_contain_r((*rI), (av[i].v>>1))) kv++;
            g->seq_vis[av[i].v] = 1;
        }
        if(kv <= 0) {
            for (i = 0; i < nv; i++) {
                if(av[i].del) continue;
                g->seq_vis[av[i].v] = 0;
            }
            continue;
        }
        for (i = 0; i < nv; i++) {
            if(av[i].del) continue;
            if(!(is_contain_r((*rI), (av[i].v>>1)))) continue;

        }

        

        if ((g->seq[v>>1].del) || (g->seq_vis[v]&1)) continue;
        if(!is_contain_r((*rI), (v>>1))) continue;
        // if(get_arcs(g, v, &w, 1) == 1) {
        //     w = g->arc[w].v;
        //     if(get_arcs(g, w^1, NULL, 0) == 1) continue;
        // }
        in0->n = 0;
        label_contain_dup(g, rI, v, in, in0);
        label_contain_dup(g, rI, v^1, in, in0);
    }
}
**/

uint32_t if_false_bub_links(uint32_t v, asg_t *g, buf_t *x, asg64_v *b, uint32_t bs, int32_t check_dist)
{
    uint32_t i, mm = 1;
    if (g->seq[v>>1].del) return 0;
    for (i = bs; i < b->n; i++) {
        g->arc[b->a[i]].del = 1;
        asg_arc_del(g, g->arc[b->a[i]].v^1, (g->arc[b->a[i]].ul>>32)^1, 1);
    }
    if (asg_arc_n(g, v) < 2 || get_arcs(g, v, NULL, 0) < 2) mm = 0;
    
    if(mm) {
        mm = 0;
        if(asg_bub_pop1_primary_trio(g, NULL, v, check_dist, x, (uint32_t)-1, (uint32_t)-1, 0, 
                                                                    NULL, NULL, NULL, 0, 0, NULL)) {

            for (i = bs; i < b->n; i++) {
                g->arc[b->a[i]].del = 0;
                asg_arc_del(g, g->arc[b->a[i]].v^1, (g->arc[b->a[i]].ul>>32)^1, 0);
            }

            asg_arc_t *av = asg_arc_a(g, v); uint32_t nv = asg_arc_n(g, v);
            for (i = 0, b->n = bs; i < nv; i++) {
                if (av[i].del) continue;
                av[i].del = 1; asg_arc_del(g, av[i].v^1, (av[i].ul>>32)^1, 1);
                kv_push(uint64_t, *b, ((uint64_t)(av-g->arc+i)));     
            }

            if(asg_bub_pop1_primary_trio(g, NULL, x->S.a[0]^1, check_dist, x, (uint32_t)-1, (uint32_t)-1, 0, NULL, NULL, NULL, 0, 0, NULL)) {
                mm = 1;
            }
        }
    }
    
    for (i = bs; i < b->n; i++) {
        g->arc[b->a[i]].del = 0;
        asg_arc_del(g, g->arc[b->a[i]].v^1, (g->arc[b->a[i]].ul>>32)^1, 0);
    }

    return mm;
}

void asg_arc_cut_bub_links(asg_t *g, asg64_v *in, float len_rat, float sec_len_rat, float ou_rat, uint32_t is_ou, uint64_t check_dist, ma_hit_t_alloc *rev, R_to_U* rI, int32_t max_ext)
{
    asg64_v tx = {0,0,0}, *b = NULL;
    uint32_t v, w, t, k, i, n_vtx = g->n_seq<<1, nv, nw, kv, kw, kol, bn, me, mu, cnt = 0, sec_check;
    asg_arc_t *av, *aw, *ref;
    buf_t x; memset(&x, 0, sizeof(x)); x.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));

    if(in) b = in;
    else b = &tx;
	b->n = 0;

    for (v = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;
        if(g->seq_vis[v] == 0) {
            av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
            if (nv < 2) continue;

            for (i = kv = kol = 0; i < nv; ++i) {
                if(av[i].del) continue;
                kv++; kol += av[i].ol; 
            }
            if(kv < 2) continue;//must have at least one exact and one inexact

            kv_push(uint64_t, *b, ((((uint64_t)(kol))<<32) | v));  
        }
    }

    if(rev && rI) memset(g->seq_vis, 0, g->n_seq*2*sizeof(uint8_t));
    radix_sort_srt64(b->a, b->a + b->n); bn = b->n;
    for (k = 0; k < bn; k++) {
        v = (uint32_t)b->a[k];
        if (g->seq[v>>1].del) continue;
        nv = asg_arc_n(g, v); av = asg_arc_a(g, v);
        if (nv < 2 || get_arcs(g, v, NULL, 0) < 2) continue;

        for (i = 0, b->n = bn, sec_check = 0; i < nv; i++) {
            if (av[i].del) continue;

            w = av[i].v^1; nw = asg_arc_n(g, w); aw = asg_arc_a(g, w);
            if(nw < 2) break;

            for (t = kw = 0, me = mu = (uint32_t)-1; t < nw; t++) {
                if(aw[t].del) continue;
                kw++; 
                if(aw[t].v == (v^1)) continue;//note: me is the shortest edge except aw[t], so here is continue
                if(aw[t].ol < me) me = aw[t].ol;
                if(aw[t].ou < mu) mu = aw[t].ou;
                kv_push(uint64_t, *b, ((uint64_t)(aw-g->arc+t)));                
            }
            if(kw < 2) break;

            if(av[i].ol > me*len_rat && av[i].ol > me*sec_len_rat) break;
            if(av[i].ol > me*len_rat) sec_check++;
            if(is_ou && av[i].ou > mu*ou_rat) break;
        }

        if(i < nv) continue;

        if(sec_check) {
            for (i = 0, ref = NULL; i < nv; i++) {//forward
                if (av[i].del) continue;
                if(!ref) {
                    ref = &(av[i]);
                } else {
                    if(trans_path_check(ref->v, av[i].v, g, rev, rI, max_ext, b)!=1) break;
                }
            }

            if(i < nv) {
                if (b->n < bn + 2) continue;///less than two edges
                for (i = bn, ref = NULL; i < b->n; i++) {
                    if(g->arc[b->a[i]].del) continue;
                    if(get_arcs(g, g->arc[b->a[i]].ul>>32, NULL, 0)!=2) break;
                    if(!ref) {
                        ref = &(g->arc[b->a[i]]);
                    } else {
                        if(trans_path_check(ref->v, g->arc[b->a[i]].v, g, rev, rI, max_ext, b)!=1) break;
                    }
                }
                if(i < b->n) continue;
            }
        }

        if(if_false_bub_links(v, g, &x, b, bn, check_dist)) {
            for (i = 0; i < nv; ++i) {
                if (av[i].del) continue;
                av[i].del = 1; asg_arc_del(g, av[i].v^1, (av[i].ul>>32)^1, 1);
            }
            cnt++;
        }
    }

    // stats_sysm(g);
    if(!in) free(tx.a);
    free(x.a); free(x.S.a); free(x.T.a); free(x.b.a); free(x.e.a);
    if(cnt > 0) asg_cleanup(g);
}

void asg_arc_cut_complex_bub_links(asg_t *g, asg64_v *in, float len_rat, float ou_rat, uint32_t is_ou, bub_label_t *b_mask_t)
{
    asg64_v tx = {0,0,0}, *b = NULL;
    uint32_t v, w, t, k, i, n_vtx = g->n_seq<<1, nv, nw, kv, kw, kol, me, mu, cnt = 0, bn;
    asg_arc_t *av, *aw;

    if(in) b = in;
    else b = &tx;
    b->n = 0;

    for (v = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;
        if(g->seq_vis[v] == 0) {
            av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
            if (nv < 2) continue;

            for (i = kv = kol = 0; i < nv; ++i) {
                if(av[i].del) continue;
                kv++; kol += av[i].ol; 
            }
            if(kv < 2) continue;//must have at least one exact and one inexact

            kv_push(uint64_t, *b, ((((uint64_t)(kol))<<32) | v));  
        }
    }

    radix_sort_srt64(b->a, b->a + b->n); bn = b->n;
    for (k = 0; k < bn; k++) {
        v = (uint32_t)b->a[k];
        if (g->seq[v>>1].del) continue;
        nv = asg_arc_n(g, v); av = asg_arc_a(g, v);
        if (nv < 2 || get_arcs(g, v, NULL, 0) < 2) continue;

        for (i = 0; i < nv; i++) {
            if (av[i].del) continue;

            w = av[i].v^1; nw = asg_arc_n(g, w); aw = asg_arc_a(g, w);
            if(nw < 2) break;

            for (t = kw = 0, me = mu = (uint32_t)-1; t < nw; t++) {
                if(aw[t].del) continue;
                kw++; 
                if(aw[t].v == (v^1)) continue;//note: me is the shortest edge except aw[t], so here is continue
                if(aw[t].ol < me) me = aw[t].ol;
                if(aw[t].ou < mu) mu = aw[t].ou;
            }
            if(kw < 2) break;

            if(av[i].ol > me*len_rat) break;
            if(is_ou && av[i].ou > mu*ou_rat) break;
        }

        if(i < nv) continue;

        for (i = 0; i < nv; ++i) {
            if (av[i].del) continue;
            av[i].del = 1; asg_arc_del(g, av[i].v^1, (av[i].ul>>32)^1, 1);
            kv_push(uint64_t, *b, ((uint64_t)(av-g->arc+i)));
        }
        // b->a[cnt++] = v;
    }

    if(b->n > bn) {
        asg_arc_identify_simple_bubbles_multi(g, b_mask_t, 0);
        for (k = bn, cnt = 0; k < b->n; k++) {
            // if(g->arc[b->a[k]].del) continue;
            v = g->arc[b->a[k]].ul>>32; w = g->arc[b->a[k]].v; 
            if(g->seq_vis[v] || g->seq_vis[v^1] || g->seq_vis[w] || g->seq_vis[w^1]) {
                cnt++; continue;
            }
            g->arc[b->a[k]].del = 0; asg_arc_del(g, g->arc[b->a[k]].v^1, (g->arc[b->a[k]].ul>>32)^1, 0);
        }
    }
    // stats_sysm(g);
    if(!in) free(tx.a);
    if(cnt > 0) asg_cleanup(g);
}

#define LIM_LEN 100

uint32_t asg_cut_semi_circ(asg_t *g, uint32_t lim_len, uint32_t is_clean)
{
    uint32_t v, t, k, e, ss, i, n_vtx = g->n_seq<<1, nv, kv, nw, cnt = 0;
    asg_arc_t *av, *aw;

    for (v = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;

        av = asg_arc_a(g, v^1); nv = asg_arc_n(g, v^1);
        if(nv <= 1) continue;
        for (i = kv = 0; i < nv; ++i) {
            if(av[i].del) continue;
            kv++; if(kv > 1) break;
        }
        if(kv <= 1) continue;

        av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
        if(nv < 1) continue;
        for (i = kv = 0; i < nv; ++i) {
            if(av[i].del) continue;
            kv++; if(kv > 1) break;
        }
        if(kv != 1) continue;

        for (i = 0; i < nv; ++i) {
            if(av[i].del) continue;
            t = follow_limit_path(g, v, &e, &ss, NULL, lim_len);
            if(ss > lim_len || t == LONG_TIPS || t == LOOP || t == END_TIPS) break;//as kv == 1
            aw = asg_arc_a(g, v^1); nw = asg_arc_n(g, v^1);
            for (k = 0; k < nw; k++) {
                if (aw[k].del) continue;
                if (aw[k].v == (e^1)) {
                    aw[k].del = 1;
                    asg_arc_del(g, aw[k].v^1, (aw[k].ul>>32)^1, 1);
                    cnt++;
                }
            }
            break; //as kv == 1
        }
    }

    if(cnt > 0 && is_clean) asg_cleanup(g);
    return cnt;
}

uint32_t asg_cut_chimeric_bub(asg_t *g, ma_hit_t_alloc* src, asg64_v *in, uint32_t normal_len, uint32_t is_clean, telo_end_t *te)
{
    asg64_v tx = {0,0,0}, *b = NULL;
    uint32_t v, w, nw, k, n_vtx = g->n_seq<<1, ei[2] = {0}, e, ss, cnt = 0;
    asg_arc_t *aw;
    if(in) b = in;
    else b = &tx;
    b->n = 0;
    // fprintf(stderr, "[M::%s]\n", __func__);
    for (v = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;
        if(te && te->hh[v>>1]) continue;
        ///note: ei[0] and ei[1] are the edge idx
        if((get_arcs(g, v, &(ei[0]), 1)!=1) || (get_arcs(g, v^1, &(ei[1]), 1)!=1)) continue;
        assert((g->arc[ei[0]].ul>>32) == v && (g->arc[ei[1]].ul>>32) == (v^1));
        if((get_arcs(g, g->arc[ei[0]].v^1, NULL, 0)!=2) || (get_arcs(g, g->arc[ei[1]].v^1, NULL, 0)!=2)) continue;
        if(!if_sup_chimeric(&(src[v>>1]), g->seq[v>>1].len, b, 1)) continue;
        w = g->arc[ei[0]].v^1;
        aw = asg_arc_a(g, w); nw = asg_arc_n(g, w);
        for (k = 0; k < nw; k++) {
            if (aw[k].del) continue;
            if (aw[k].v == (v^1)) {
                ss = k; continue;
            }
            break;
        }
        // assert(aw[ss].v == (v^1));//this assert does not work, just ignore
        if(follow_limit_path(g, aw[k].v, &e, &ss, NULL, (uint32_t)-1) != TWO_INPUT) continue;
        if(ss > normal_len) {
            w = e;
            aw = asg_arc_a(g, w); nw = asg_arc_n(g, w);
            for (k = ss = 0; k < nw; k++) {
                if (aw[k].del) continue;
                ss++; e = aw[k].v;
                if(ss > 1) break;
            }
            if(ss == 1 && e == g->arc[ei[1]].v) asg_seq_del(g, v>>1), cnt++;
        }
    }

    if(!in) free(tx.a);
    if (is_clean && cnt > 0) asg_cleanup(g);
    return cnt;
}

void asg_iterative_semi_circ(asg_t *g, ma_hit_t_alloc* src, asg64_v *in, uint32_t normal_len, uint32_t pop_chimer, telo_end_t *te)
{
    uint64_t occ = 0, s = 1;
    while (s) {
        s = asg_cut_semi_circ(g, LIM_LEN, 0);
        if(pop_chimer) s = s + asg_cut_chimeric_bub(g, src, in, normal_len, 0, te);
        occ += s;
    }

    // stats_sysm(g);
    if(occ) asg_cleanup(g);
}

uint32_t asg_cut_large_indel(asg_t *g, asg64_v *in, int32_t max_ext, float ou_rat, uint32_t is_ou, uint32_t min_diff)
{
    asg64_v tx = {0,0,0}, *b = NULL;
    uint32_t v, w, n_vtx = g->n_seq<<1, i, k, kv, kw, nv, nw, ou_max, ol_max, to_del, cnt = 0;
    asg_arc_t *av, *aw, *ve, *we;
    if(in) b = in;
    else b = &tx;
    b->n = 0;

    for (v = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;
        if(g->seq_vis[v] == 0) {
            av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
            if (nv < 2) continue;

            for (i = kv = 0; i < nv; ++i) {
                if(av[i].del) continue;
                kv++;
            }
            if(kv < 2) continue;

            for (i = 0; i < nv; ++i) {
                if(av[i].del || av[i].no_l_indel) continue;
                ///means there is a large indel at this edge
                kv_push(uint64_t, *b, (((uint64_t)av[i].ol)<<32) | ((uint64_t)(av-g->arc+i)));   
            }
        }
    }

    radix_sort_srt64(b->a, b->a + b->n);
    for (k = 0; k < b->n; k++) {
        if(g->arc[(uint32_t)b->a[k]].del) continue;

        v = g->arc[(uint32_t)b->a[k]].ul>>32; w = g->arc[(uint32_t)b->a[k]].v^1;
        if(g->seq[v>>1].del || g->seq[w>>1].del) continue;
        nv = asg_arc_n(g, v); nw = asg_arc_n(g, w);
        av = asg_arc_a(g, v); aw = asg_arc_a(g, w);
        if(nv<=1 && nw <= 1) continue;

        for (i = kv = ou_max = ol_max = 0, ve = NULL; i < nv; ++i) {
            if(av[i].del) continue;
            if(av[i].v == (w^1)) ve = &(av[i]);
            kv++;
            if(ou_max < av[i].ou) ou_max = av[i].ou;
            if(ol_max < av[i].ol) ol_max = av[i].ol;
        }
        if (kv < 1) continue;
        if (kv >= 2) {
            if (is_ou && ve->ou > ou_max*ou_rat) continue;
            if ((ve->ol + min_diff) > ol_max) continue;
        }
        

        for (i = kw = ou_max = ol_max = 0, we = NULL; i < nw; ++i) {
            if(aw[i].del) continue;
            if(aw[i].v == (v^1)) we = &(aw[i]);
            kw++;
            if(ou_max < aw[i].ou) ou_max = aw[i].ou;
            if(ol_max < aw[i].ol) ol_max = aw[i].ol;
        }
        if (kw < 1) continue;
        if (kw >= 2) {
            if (is_ou && we->ou > ou_max*ou_rat) continue;
            if ((we->ol + min_diff) > ol_max) continue;
        }

        if (kv <= 1 && kw <= 1) continue;

        to_del = 0;
        if (kv > 1 && kw > 1) {
            to_del = 1;
        } else if (kw == 1) {
            if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
        } else if (kv == 1) {
            if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
        }

        if (to_del) {
            ve->del = we->del = 1, ++cnt;
        }
    }
    // stats_sysm(g);
    if(!in) free(tx.a);
    if (cnt > 0) asg_cleanup(g);
    return cnt;
}

void debug_edges(asg64_v *dbg, uint32_t *l, uint32_t l_n) {
    uint32_t k, k_n, i, m;
    for (i = k = 0; i < l_n; i++) {
        fprintf(stderr, "# gid-%u: %u\n", i, l[i]);
        for (k_n = k + l[i]; k < k_n; k++) {
            dbg->a[k] <<= 32; dbg->a[k] += i;
        }
    }
    fprintf(stderr, "# dbg->n: %u\n", (uint32_t)dbg->n);
    
    radix_sort_srt64(dbg->a, dbg->a + dbg->n);
    for (k = 1, i = 0; k <= dbg->n; k++) {
        if(k == dbg->n || (dbg->a[k]>>32) != (dbg->a[i]>>32)) {
            if(k - i < l_n) {
                
                for (m = i; m < k; m++) {
                    fprintf(stderr, "eid->%lu, gid->%u\n", dbg->a[m]>>32, (uint32_t)dbg->a[m]);
                }
            }
            i = k;
        }
    }
}

void print_node(asg_t *sg, uint32_t src)
{
    asg_arc_t *av; uint32_t nv, v, i;
    v = src<<1;
    av = asg_arc_a(sg, v); nv = asg_arc_n(sg, v);
    fprintf(stderr, "\n%.*s(%c)\tnv:%u\n", 
                (int32_t)Get_NAME_LENGTH(R_INF, v>>1), Get_NAME(R_INF, v>>1), "+-"[v&1], nv);
    for (i = 0; i < nv; i++) {
        fprintf(stderr, "++++++%.*s(%c)<id:%lu>\t%.*s(%c)<id:%u>\tol:%u\tou:%u\tdel:%u\n", 
                (int32_t)Get_NAME_LENGTH(R_INF, (av[i].ul>>33)), Get_NAME(R_INF, (av[i].ul>>33)), "+-"[(av[i].ul>>32)&1], av[i].ul>>33, 
                (int32_t)Get_NAME_LENGTH(R_INF, (av[i].v>>1)), Get_NAME(R_INF, (av[i].v>>1)), "+-"[av[i].v&1], av[i].v>>1, av[i].ol, av[i].ou, av[i].del);
    }

    v = (src<<1)+1;
    av = asg_arc_a(sg, v); nv = asg_arc_n(sg, v);
    fprintf(stderr, "\n%.*s(%c)\tnv:%u\n", 
                (int32_t)Get_NAME_LENGTH(R_INF, v>>1), Get_NAME(R_INF, v>>1), "+-"[v&1], nv);
    for (i = 0; i < nv; i++) {
        fprintf(stderr, "------%.*s(%c)<id:%lu>\t%.*s(%c)<id:%u>\tol:%u\tou:%u\tdel:%u\n", 
                (int32_t)Get_NAME_LENGTH(R_INF, (av[i].ul>>33)), Get_NAME(R_INF, (av[i].ul>>33)), "+-"[(av[i].ul>>32)&1], av[i].ul>>33, 
                (int32_t)Get_NAME_LENGTH(R_INF, (av[i].v>>1)), Get_NAME(R_INF, (av[i].v>>1)), "+-"[av[i].v&1], av[i].v>>1, av[i].ol, av[i].ou, av[i].del);
    }
}

void print_vw_edge(asg_t *sg, uint32_t vid, uint32_t wid, const char *cmd)
{
    asg_arc_t *av; uint32_t nv, i, sid, eid;
    if(vid >= sg->n_seq || wid >= sg->n_seq) return;

    sid = vid; eid = wid;
    av = asg_arc_a(sg, (sid<<1)); nv = asg_arc_n(sg, (sid<<1));
    for (i = 0; i < nv; i++) {
        if((av[i].v>>1) == eid) {
            fprintf(stderr, "[%s]\t%.*s(%c)<id:%lu>\t%.*s(%c)<id:%u>\tol:%u\tou:%u\tdel:%u\n", cmd, 
                (int32_t)Get_NAME_LENGTH(R_INF, (av[i].ul>>33)), Get_NAME(R_INF, (av[i].ul>>33)), "+-"[(av[i].ul>>32)&1], av[i].ul>>33, 
                (int32_t)Get_NAME_LENGTH(R_INF, (av[i].v>>1)), Get_NAME(R_INF, (av[i].v>>1)), "+-"[av[i].v&1], av[i].v>>1, av[i].ol, av[i].ou, av[i].del);
            break;
        }
    }
    av = asg_arc_a(sg, ((sid<<1)+1)); nv = asg_arc_n(sg, ((sid<<1)+1));
    for (i = 0; i < nv; i++) {
        if((av[i].v>>1) == eid) {
            fprintf(stderr, "[%s]\t%.*s(%c)<id:%lu>\t%.*s(%c)<id:%u>\tol:%u\tou:%u\tdel:%u\n", cmd, 
                (int32_t)Get_NAME_LENGTH(R_INF, (av[i].ul>>33)), Get_NAME(R_INF, (av[i].ul>>33)), "+-"[(av[i].ul>>32)&1], av[i].ul>>33, 
                (int32_t)Get_NAME_LENGTH(R_INF, (av[i].v>>1)), Get_NAME(R_INF, (av[i].v>>1)), "+-"[av[i].v&1], av[i].v>>1, av[i].ol, av[i].ou, av[i].del);
            break;
        }
    }

    sid = wid; eid = vid;
    av = asg_arc_a(sg, (sid<<1)); nv = asg_arc_n(sg, (sid<<1));
    for (i = 0; i < nv; i++) {
        if((av[i].v>>1) == eid) {
            fprintf(stderr, "[%s]\t%.*s(%c)<id:%lu>\t%.*s(%c)<id:%u>\tol:%u\tou:%u\tdel:%u\n", cmd, 
                (int32_t)Get_NAME_LENGTH(R_INF, (av[i].ul>>33)), Get_NAME(R_INF, (av[i].ul>>33)), "+-"[(av[i].ul>>32)&1], av[i].ul>>33, 
                (int32_t)Get_NAME_LENGTH(R_INF, (av[i].v>>1)), Get_NAME(R_INF, (av[i].v>>1)), "+-"[av[i].v&1], av[i].v>>1, av[i].ol, av[i].ou, av[i].del);
            break;
        }
    }
    av = asg_arc_a(sg, ((sid<<1)+1)); nv = asg_arc_n(sg, ((sid<<1)+1));
    for (i = 0; i < nv; i++) {
        if((av[i].v>>1) == eid) {
            fprintf(stderr, "[%s]\t%.*s(%c)<id:%lu>\t%.*s(%c)<id:%u>\tol:%u\tou:%u\tdel:%u\n", cmd, 
                (int32_t)Get_NAME_LENGTH(R_INF, (av[i].ul>>33)), Get_NAME(R_INF, (av[i].ul>>33)), "+-"[(av[i].ul>>32)&1], av[i].ul>>33, 
                (int32_t)Get_NAME_LENGTH(R_INF, (av[i].v>>1)), Get_NAME(R_INF, (av[i].v>>1)), "+-"[av[i].v&1], av[i].v>>1, av[i].ol, av[i].ou, av[i].del);
            break;
        }
    }
    // if(i >= nv) fprintf(stderr, "[%s]\tno edges\n", cmd);
}

void prt_spec_edge(asg_t *rg, ma_hit_t_alloc *src, uint32_t tot_rid, uint32_t vid, uint32_t wid, ug_opt_t *uopt, const char *cmd)
{
    if(vid >= tot_rid) return;
    uint32_t k, qn, tn; int32_t r; asg_arc_t p; ma_hit_t_alloc *s = &(src[vid]); 
    for (k = 0; k < s->length; k++) {
        qn = Get_qn(s->buffer[k]); tn = Get_tn(s->buffer[k]);
        if(qn == vid && tn == wid) {
            r = ma_hit2arc(&(s->buffer[k]), rg->seq[qn].len, rg->seq[tn].len, uopt->max_hang, asm_opt.max_hang_rate, uopt->min_ovlp, &p);
            fprintf(stderr, "[M::%s::]%s\tqn::%u(%c)\tq::[%u,\t%u)\t%c\ttn::%u(%c)\tt::[%u,\t%u)\n", 
                    __func__, cmd, qn, (r>=0)?("+-"[(p.ul>>32)&1]):('*'), Get_qs(s->buffer[k]), Get_qe(s->buffer[k]), 
                    "+-"[s->buffer[k].rev], 
                    tn, (r>=0)?("+-"[p.v&1]):('*'), Get_ts(s->buffer[k]), Get_te(s->buffer[k]));
        }
    }
}

int32_t gen_spec_edge(asg_t *rg, ug_opt_t *uopt, uint32_t v, uint32_t w, asg_arc_t *t)
{
    uint32_t k, qn, tn; int32_t r; ma_hit_t_alloc *s = &(uopt->sources[v>>1]); asg_arc_t p;
    for (k = 0; k < s->length; k++) {
        qn = Get_qn(s->buffer[k]); tn = Get_tn(s->buffer[k]);
        if(tn != (w>>1)) continue;
        r = ma_hit2arc(&(s->buffer[k]), rg->seq[qn].len, rg->seq[tn].len, uopt->max_hang, asm_opt.max_hang_rate, uopt->min_ovlp, &p);
        if(r < 0) continue;
        if((p.ul>>32) != v || p.v != w) continue;
        *t = p; t->ou = 0;
        return 1;
    }
    return -1;
}

void filter_sg_by_ug(asg_t *rg, ma_ug_t *ug, ug_opt_t *uopt)
{
    uint32_t i, m, v, w, nv, n_vx, vx, wx; int32_t r;
    asg_arc_t *av = NULL; ma_utg_t *u = NULL; asg_arc_t *p, t;
    n_vx = rg->n_seq; rg->n_arc = 0;
    for (v = 0; v < n_vx; v++) rg->seq[v].del = (!!1);

    for (i = 0; i < ug->g->n_seq; ++i) {
        ug->g->seq[i].c = PRIMARY_LABLE;;
        if(ug->g->seq[i].del) continue;
        u = &(ug->u.a[i]);
        for (m = 0; m < u->n; m++) rg->seq[u->a[m]>>33].del = (!!0);
        for (m = 0; (m + 1) < u->n; m++) {
            v = u->a[m]>>32; w = u->a[m+1]>>32;
            r = gen_spec_edge(rg, uopt, v, w, &t);
            assert(r >= 0); p = asg_arc_pushp(rg); *p = t;

            r = gen_spec_edge(rg, uopt, w^1, v^1, &t);
            assert(r >= 0); p = asg_arc_pushp(rg); *p = t;
        }


        v = i<<1; nv = asg_arc_n(ug->g, v); av = asg_arc_a(ug->g, v);
        for (m = 0; m < nv; m++) {
            if(av[m].del) continue;
            w = av[m].v;
            vx = (v&1?((ug->u.a[v>>1].a[0]>>32)^1):(ug->u.a[v>>1].a[ug->u.a[v>>1].n-1]>>32));
            wx = (w&1?((ug->u.a[w>>1].a[ug->u.a[w>>1].n-1]>>32)^1):(ug->u.a[w>>1].a[0]>>32));

            r = gen_spec_edge(rg, uopt, vx, wx, &t);
            assert(r >= 0); p = asg_arc_pushp(rg); *p = t;

            // r = gen_spec_edge(rg, uopt, wx^1, vx^1, &t);
            // assert(r >= 0); p = asg_arc_pushp(rg); *p = t;
        }

        v = (i<<1)+1; nv = asg_arc_n(ug->g, v); av = asg_arc_a(ug->g, v);
        for (m = 0; m < nv; m++) {
            if(av[m].del) continue;
            w = av[m].v;
            vx = (v&1?((ug->u.a[v>>1].a[0]>>32)^1):(ug->u.a[v>>1].a[ug->u.a[v>>1].n-1]>>32));
            wx = (w&1?((ug->u.a[w>>1].a[ug->u.a[w>>1].n-1]>>32)^1):(ug->u.a[w>>1].a[0]>>32));

            r = gen_spec_edge(rg, uopt, vx, wx, &t);
            assert(r >= 0); p = asg_arc_pushp(rg); *p = t;

            // r = gen_spec_edge(rg, uopt, wx^1, vx^1, &t);
            // assert(r >= 0); p = asg_arc_pushp(rg); *p = t;
        }

        if(u->circ) {
            v = w = i<<1;
            vx = (v&1?((ug->u.a[v>>1].a[0]>>32)^1):(ug->u.a[v>>1].a[ug->u.a[v>>1].n-1]>>32));
            wx = (w&1?((ug->u.a[w>>1].a[ug->u.a[w>>1].n-1]>>32)^1):(ug->u.a[w>>1].a[0]>>32));
            r = gen_spec_edge(rg, uopt, vx, wx, &t);
            assert(r >= 0); p = asg_arc_pushp(rg); *p = t;

            v = w = (i<<1)^1;
            vx = (v&1?((ug->u.a[v>>1].a[0]>>32)^1):(ug->u.a[v>>1].a[ug->u.a[v>>1].n-1]>>32));
            wx = (w&1?((ug->u.a[w>>1].a[ug->u.a[w>>1].n-1]>>32)^1):(ug->u.a[w>>1].a[0]>>32));
            r = gen_spec_edge(rg, uopt, vx, wx, &t);
            assert(r >= 0); p = asg_arc_pushp(rg); *p = t;
        }
    }

    free(rg->idx);
    rg->idx = 0;
    rg->is_srt = 0;
    asg_cleanup(rg);
    asg_symm(rg);


    /*******************************for debug************************************/
    // ma_ug_t *dbg = ma_ug_gen(rg);
    // for (i = 0; i < dbg->g->n_seq; ++i) dbg->g->seq[i].c = PRIMARY_LABLE;  
    // cmp_untig_graph(dbg, ug);
    /*******************************for debug************************************/
}

void prt_specfic_sge(asg_t *g, uint32_t src, uint32_t dst, const char* cmd)
{
    uint32_t k, v, w, nv; asg_arc_t *av;
    fprintf(stderr, "[M::%s::%s] src::%.*s(id::%u), dst::%.*s(id::%u)\n", __func__, cmd, 
    (int)Get_NAME_LENGTH(R_INF, src), Get_NAME(R_INF, src), src, 
    (int)Get_NAME_LENGTH(R_INF, dst), Get_NAME(R_INF, dst), dst);
    v = src<<1; 
    nv = asg_arc_n(g, v); av = asg_arc_a(g, v);
    for (k = 0; k < nv; ++k) {
        if ((av[k].v>>1) == dst) {
            w = av[k].v;
            fprintf(stderr, "[M::%s::]\t%.*s(%c)\t%.*s(%c)\tou::%u\n", __func__,
					(int)Get_NAME_LENGTH(R_INF, (v>>1)), Get_NAME(R_INF, (v>>1)), "+-"[v&1],
					(int)Get_NAME_LENGTH(R_INF, (w>>1)), Get_NAME(R_INF, (w>>1)), "+-"[w&1], av[k].ou);
        }
    }

    v = (src<<1)+1; 
    nv = asg_arc_n(g, v); av = asg_arc_a(g, v);
    for (k = 0; k < nv; ++k) {
        if ((av[k].v>>1) == dst) {
            w = av[k].v;
            fprintf(stderr, "[M::%s::]\t%.*s(%c)\t%.*s(%c)\tou::%u\n", __func__,
					(int)Get_NAME_LENGTH(R_INF, (v>>1)), Get_NAME(R_INF, (v>>1)), "+-"[v&1],
					(int)Get_NAME_LENGTH(R_INF, (w>>1)), Get_NAME(R_INF, (w>>1)), "+-"[w&1], av[k].ou);
        }
    }
}

flex_asg_t *init_flex_asg_t(asg_t *g, ma_hit_t_alloc* src, int64_t min_ovlp, int64_t max_hang, int64_t max_hang_rate, int64_t gap_fuzz)
{
    flex_asg_t *z; CALLOC(z, 1);
    z->g = g; z->src = src; z->min_ovlp = min_ovlp; z->gap_fuzz = gap_fuzz;
    z->max_hang = max_hang; z->max_hang_rate = max_hang_rate; 
    MALLOC(z->idx, (z->g->n_seq<<1)); memset(z->idx, -1, (z->g->n_seq<<1)*sizeof(*(z->idx)));
    return z;
}

void des_flex_asg_t(flex_asg_t *z)
{
    free(z->idx); free(z->pi.a); free(z->a);
}


void print_raw_u2rgfa_seq(all_ul_t *aln, R_to_U* rI, uint32_t is_detail)
{
	uint64_t id, a_n, k, z; uc_block_t *a = NULL; 
	kvec_t(uint8_t) ff; kv_init(ff);
	for (id = 0; id < aln->n; id++) {
		a = aln->a[id].bb.a; a_n = aln->a[id].bb.n;
		if(a_n == 0) continue;
		fprintf(stderr,"\n%.*s\tid::%lu\trlen::%u", (int32_t)aln->nid.a[id].n, aln->nid.a[id].a, id, aln->a[id].rlen);
		kv_resize(uint8_t, ff, a_n); memset(ff.a, 0, a_n*sizeof((*(ff.a))));
		if(is_detail) {
			fprintf(stderr, "\n");
			for (k = 0; k < a_n; k++) {
				if(ff.a[k]) continue;
				for (z = k; z != (uint32_t)-1; z = a[z].aidx) {
					fprintf(stderr, "%.*s\t%c\tq::[%u, %u)\tt::[%u, %u)\tid::%u\ttl::%lu\tc::%u\n", 
					(int)Get_NAME_LENGTH(R_INF, a[z].hid), Get_NAME(R_INF, a[z].hid), "+-"[a[z].rev],
					a[z].qs, a[z].qe, a[z].ts, a[z].te, a[z].hid, Get_READ_LENGTH(R_INF, a[z].hid),
                    rI?is_contain_r((*rI), a[z].hid):0);
					assert(ff.a[z] == 0);
					ff.a[z] = 1;
				}
				fprintf(stderr, "************\n");
			}
		} else {
			fprintf(stderr, "\t");
			for (k = 0; k < a_n; k++) {
				if(ff.a[k]) continue;
				for (z = k; z != (uint32_t)-1; z = a[z].aidx) {
					fprintf(stderr, "%.*s\t", 
					(int)Get_NAME_LENGTH(R_INF, a[z].hid), Get_NAME(R_INF, a[z].hid));
					assert(ff.a[z] == 0);
					ff.a[z] = 1;
				}
				fprintf(stderr, "\n");
			}
		}
	}
	kv_destroy(ff);
}


void post_rescue(ug_opt_t *uopt, asg_t *sg, ma_hit_t_alloc *src, ma_hit_t_alloc *rev, R_to_U* rI, bub_label_t *b_mask_t, long long no_trio_recover)
{
    rescue_contained_reads_aggressive(NULL, sg, src, uopt->coverage_cut, rI, uopt->max_hang, uopt->min_ovlp, 10, 1, 0, NULL, NULL, b_mask_t);
    rescue_missing_overlaps_aggressive(NULL, sg, src, uopt->coverage_cut, rI, uopt->max_hang, uopt->min_ovlp, 1, 0, NULL, b_mask_t);
    rescue_missing_overlaps_backward(NULL, sg, src, uopt->coverage_cut, rI, uopt->max_hang, uopt->min_ovlp, 10, 1, 0, b_mask_t);
    // rescue_wrong_overlaps_to_unitigs(NULL, sg, sources, reverse_sources, coverage_cut, ruIndex, 
    // max_hang_length, mini_overlap_length, bubble_dist, NULL);
    // rescue_no_coverage_aggressive(sg, sources, reverse_sources, &coverage_cut, ruIndex, max_hang_length, 
    // mini_overlap_length, bubble_dist, 10);
    set_hom_global_coverage(&asm_opt, sg, uopt->coverage_cut, src, rev, rI, uopt->max_hang, uopt->min_ovlp);
    rescue_bubble_by_chain(sg, uopt->coverage_cut, src, rev, (asm_opt.max_short_tip*2), 0.15, 3, rI, 0.05, 0.9, uopt->max_hang, 
    uopt->min_ovlp, 10, uopt->gap_fuzz, b_mask_t, no_trio_recover);
}

void ul_clean_gfa(ug_opt_t *uopt, asg_t *sg, ma_hit_t_alloc *src, ma_hit_t_alloc *rev, R_to_U* rI, int64_t clean_round, double min_ovlp_drop_ratio, double max_ovlp_drop_ratio, 
double ou_drop_rate, int64_t max_tip, int64_t gap_fuzz, bub_label_t *b_mask_t, int32_t is_ou, int32_t is_trio, uint32_t ou_thres, char *o_file)
{
    #define HARD_OU_DROP 0.75
    #define HARD_OL_DROP 0.6
    #define HARD_OL_SEC_DROP 0.85
    #define HARD_ORTHOLOGY_DROP 0.4
    double step = 
        (clean_round==1?max_ovlp_drop_ratio:((max_ovlp_drop_ratio-min_ovlp_drop_ratio)/(clean_round-1)));
    double drop = min_ovlp_drop_ratio;
    int64_t i; asg64_v bu = {0,0,0}, ba = {0,0,0}; uint32_t l_drop = 2000; flex_asg_t *fg = NULL; uint32_t min_diff = 0, step_diff = 2000;
    if(is_ou) fg = init_flex_asg_t(sg, uopt->sources, uopt->min_ovlp, uopt->max_hang, asm_opt.max_hang_rate, gap_fuzz);
    // if(is_ou) update_sg_uo(sg, src);///do not do it here
    // print_debug_gfa(sg, NULL, uopt->coverage_cut, "UL.dirty.debug", uopt->sources, uopt->ruIndex, uopt->max_hang, uopt->min_ovlp, 1, 0, 0);
    // exit(1);
    // print_raw_u2rgfa_seq(&UL_INF, rI, 1);
    // exit(1);
    // fprintf(stderr, "%.*s\tid::%u\tis_c::%u\n", 
    //                 (int)Get_NAME_LENGTH(R_INF, 10785), Get_NAME(R_INF, 10785), 10785, is_contain_r((*rI), 10785));
    // fprintf(stderr, "%.*s\tid::%u\tis_c::%u\n", 
    //             (int)Get_NAME_LENGTH(R_INF, 10790), Get_NAME(R_INF, 10790), 10790, is_contain_r((*rI), 10790));
    // fprintf(stderr, "%.*s\tid::%u\tis_c::%u\n", 
    //         (int)Get_NAME_LENGTH(R_INF, 10805), Get_NAME(R_INF, 10805), 10805, is_contain_r((*rI), 10805));
    // fprintf(stderr, "%.*s\tid::%u\tis_c::%u\n", 
    //         (int)Get_NAME_LENGTH(R_INF, 10809), Get_NAME(R_INF, 10809), 10809, is_contain_r((*rI), 10809));
    // fprintf(stderr, "%.*s\tid::%u\tis_c::%u\n", 
    //         (int)Get_NAME_LENGTH(R_INF, 10819), Get_NAME(R_INF, 10819), 10819, is_contain_r((*rI), 10819));
    // debug_info_of_specfic_node("m64011_190830_220126/47516220/ccs", sg, rI, "beg-0");
    // debug_info_of_specfic_node("m64012_190920_173625/163644465/ccs", sg, rI, "beg-0");

	asg_arc_cut_tips(sg, max_tip, &bu, is_ou, is_ou?rI:NULL, uopt->te);///p_telo
    // fprintf(stderr, "[M::%s] count_edges_v_w(sg, 49778, 49847)->%ld\n", __func__, count_edges_v_w(sg, 49778, 49847));
    // if(is_ou) dedup_contain_g(uopt, sg);
    for (i = 0; i < clean_round; i++, drop += step) {
        if(drop > max_ovlp_drop_ratio) drop = max_ovlp_drop_ratio;
        if(is_ou) {
            if(drop <= 0.500001) min_diff = step_diff>>1;
            else min_diff = step_diff;
        }
        // fprintf(stderr, "(0):i->%ld, drop->%f\n", i, drop);
        // prt_specfic_sge(sg, 10531, 10519, "--0--");

        // print_vw_edge(sg, 34156, 34090, "0");
        // stats_chimeric(sg, src, &bu);
        if(!is_ou) asg_iterative_semi_circ(sg, src, &bu, max_tip, 1, uopt->te);///p_telo
        asg_arc_identify_simple_bubbles_multi(sg, b_mask_t, 1);
        asg_arc_cut_chimeric(sg, src, &bu, is_ou?ou_thres:(uint32_t)-1, uopt->te);///p_telo

        asg_arc_cut_tips(sg, max_tip, &bu, is_ou, is_ou?rI:NULL, uopt->te);
        // prt_specfic_sge(sg, 10531, 10519, "--1--");
        asg_arc_identify_simple_bubbles_multi(sg, b_mask_t, 0);
        asg_arc_cut_inexact(sg, src, &bu, max_tip, is_ou, is_trio, min_diff, ou_drop_rate/**, NULL**//**&dbg**/);
        // debug_edges(&dbg, d, 2);
        asg_arc_cut_tips(sg, max_tip, &bu, is_ou, is_ou?rI:NULL, uopt->te);
        // prt_specfic_sge(sg, 10531, 10519, "--2--");

        asg_arc_identify_simple_bubbles_multi(sg, b_mask_t, 1);
        asg_arc_cut_length(sg, &bu, max_tip, drop, ou_drop_rate, is_ou, is_trio, 1, min_diff, 1, NULL, NULL, NULL);
        asg_arc_cut_tips(sg, max_tip, &bu, is_ou, is_ou?rI:NULL, uopt->te);

        // prt_specfic_sge(sg, 10531, 10519, "--3--");
        // if(is_ou) asg_arc_cut_contain(fg, &bu, &ba, rI, ((i+1)<clean_round)?ou_drop_rate:-1);
        if(is_ou) {
            asg_arc_cut_contain(fg, &bu, &ba, rI, ou_drop_rate, 0);
            // dedup_contain_g(uopt, sg);
        }

        asg_arc_identify_simple_bubbles_multi(sg, b_mask_t, 1);
        asg_arc_cut_bub_links(sg, &bu, HARD_OL_DROP, HARD_OL_SEC_DROP, HARD_OU_DROP, is_ou, asm_opt.large_pop_bubble_size, rev, rI, max_tip);
        // prt_specfic_sge(sg, 10531, 10519, "--4--");

        asg_arc_identify_simple_bubbles_multi(sg, b_mask_t, 1);
        asg_arc_cut_complex_bub_links(sg, &bu, HARD_OL_DROP, HARD_OU_DROP, is_ou, b_mask_t);
        asg_arc_cut_tips(sg, max_tip, &bu, is_ou, is_ou?rI:NULL, uopt->te);
        // prt_specfic_sge(sg, 10531, 10519, "--5--");

        // if(i == 3) {
        //     print_debug_gfa(sg, NULL, uopt->coverage_cut, "UL.dirty3.debug", uopt->sources, uopt->ruIndex, uopt->max_hang, uopt->min_ovlp, 1, 0, 0);
        //     // exit(1);
        // }

        if(is_ou) {
            if(ul_refine_alignment(uopt, sg)) update_sg_uo(sg, src);
            if(clean_contain_g(uopt, sg, 1)) update_sg_uo(sg, src);
        }
    }

    if(is_ou) min_diff = step_diff;
    // print_debug_gfa(sg, NULL, uopt->coverage_cut, "UL.dirty4.debug", uopt->sources, uopt->ruIndex, uopt->max_hang, uopt->min_ovlp, 1, 0, 0);
    // debug_info_of_specfic_node("m64012_190921_234837/111673711/ccs", sg, rI, "end");
    // debug_info_of_specfic_node("m64011_190830_220126/95028102/ccs", sg, rI, "end");
    if(is_ou) {
        asg_arc_cut_contain(fg, &bu, &ba, rI, ou_drop_rate, 0);
        asg_arc_cut_contain(fg, &bu, &ba, rI, -1, 1);
        // dedup_contain_g(uopt, sg);
        if(clean_contain_g(uopt, sg, 1)) update_sg_uo(sg, src);
    }
    if(!is_ou) asg_iterative_semi_circ(sg, src, &bu, max_tip, 1, uopt->te);

    // print_debug_gfa(sg, NULL, uopt->coverage_cut, "UL.dirty5.debug", uopt->sources, uopt->ruIndex, uopt->max_hang, uopt->min_ovlp, 1, 0, 0);

    asg_arc_identify_simple_bubbles_multi(sg, b_mask_t, 0);
    asg_cut_large_indel(sg, &bu, max_tip, HARD_OU_DROP, is_ou, min_diff);///shoule we ignore ou here?
    asg_arc_cut_tips(sg, max_tip, &bu, is_ou, is_ou?rI:NULL, uopt->te);

    // print_debug_gfa(sg, NULL, uopt->coverage_cut, "UL.dirty6.debug", uopt->sources, uopt->ruIndex, uopt->max_hang, uopt->min_ovlp, 1, 0, 0);

    if(!is_ou) {
        ///asg_arc_del_triangular_directly might be unnecessary
        asg_arc_identify_simple_bubbles_multi(sg, b_mask_t, 0);
        asg_arc_cut_length(sg, &bu, max_tip, HARD_ORTHOLOGY_DROP/**min_ovlp_drop_ratio**/, ou_drop_rate, is_ou, 0/**is_trio**/, is_ou?1:0, min_diff, 1, rev, rI, NULL);
        asg_arc_cut_tips(sg, max_tip, &bu, is_ou, is_ou?rI:NULL, uopt->te);

        // print_debug_gfa(sg, NULL, uopt->coverage_cut, "UL.dirty7.debug", uopt->sources, uopt->ruIndex, uopt->max_hang, uopt->min_ovlp, 1, 0, 0);

        asg_arc_identify_simple_bubbles_multi(sg, b_mask_t, 0);
        asg_arc_cut_length(sg, &bu, max_tip, min_ovlp_drop_ratio, ou_drop_rate, is_ou, 0/**is_trio**/, is_ou?1:0, min_diff, 1, rev, rI, &l_drop);
        asg_arc_cut_tips(sg, max_tip, &bu, is_ou, is_ou?rI:NULL, uopt->te);
    } else {
        min_diff = step_diff; l_drop = 6000;
        asg_arc_identify_simple_bubbles_multi(sg, b_mask_t, 0);
        asg_arc_cut_length(sg, &bu, max_tip, 0.3, 0.9, is_ou, is_trio, 1, min_diff, 8, NULL, NULL, &l_drop);
        asg_arc_cut_tips(sg, max_tip, &bu, is_ou, is_ou?rI:NULL, uopt->te);
    }

    // print_debug_gfa(sg, NULL, uopt->coverage_cut, "UL.dirty8.debug", uopt->sources, uopt->ruIndex, uopt->max_hang, uopt->min_ovlp, 1, 0, 0);

    if(!is_ou) asg_cut_semi_circ(sg, LIM_LEN, 1);
    /**
    rescue_contained_reads_aggressive(NULL, sg, src, uopt->coverage_cut, rI, uopt->max_hang, uopt->min_ovlp, 10, 1, 0, NULL, NULL, b_mask_t);
    rescue_missing_overlaps_aggressive(NULL, sg, src, uopt->coverage_cut, rI, uopt->max_hang, uopt->min_ovlp, 1, 0, NULL, b_mask_t);
    rescue_missing_overlaps_backward(NULL, sg, src, uopt->coverage_cut, rI, uopt->max_hang, uopt->min_ovlp, 10, 1, 0, b_mask_t);
    // rescue_wrong_overlaps_to_unitigs(NULL, sg, sources, reverse_sources, coverage_cut, ruIndex, 
    // max_hang_length, mini_overlap_length, bubble_dist, NULL);
    // rescue_no_coverage_aggressive(sg, sources, reverse_sources, &coverage_cut, ruIndex, max_hang_length, 
    // mini_overlap_length, bubble_dist, 10);
    set_hom_global_coverage(&asm_opt, sg, uopt->coverage_cut, src, rev, rI, uopt->max_hang, uopt->min_ovlp);
    rescue_bubble_by_chain(sg, uopt->coverage_cut, src, rev, (asm_opt.max_short_tip*2), 0.15, 3, rI, 0.05, 0.9, uopt->max_hang, uopt->min_ovlp, 10, uopt->gap_fuzz, b_mask_t);
    **/
    post_rescue(uopt, sg, src, rev, rI, b_mask_t, is_ou);
    
    ug_ext_gfa(uopt, sg, ug_ext_len);

    // if(is_ou) dedup_contain_g(uopt, sg);
    // exit(1)

    output_unitig_graph(sg, uopt->coverage_cut, o_file, src, rI, uopt->max_hang, uopt->min_ovlp);
    // exit(1);
    // flat_bubbles(sg, ruIndex->is_het); free(ruIndex->is_het); ruIndex->is_het = NULL;
    flat_soma_v(sg, src, rI);

    ///note: although the above functions will not change the UL part, but it will change the read graph
    ///so it is necessary to run update_sg_uo
    if(is_ou) {
        update_sg_uo(sg, src);
    }
    if(is_ou) {
        des_flex_asg_t(fg); free(fg);
    }
    // print_debug_gfa(sg, NULL, uopt->coverage_cut, "UL.debug", uopt->sources, uopt->ruIndex, uopt->max_hang, uopt->min_ovlp, 1, 0, 0);
    // exit(1);
    // print_node(sg, 17078); //print_node(sg, 8311); print_node(sg, 8294);
    
    free(bu.a); free(ba.a);
}

int32_t gen_ext_tip(uint32_t v, asg_t *sg, ug_opt_t *uopt, uint8_t *ff, uint32_t min_ovlp, uint64_t *res, uint64_t *path_len)
{
    uint32_t k, qn, tn, cc, ccu; int32_t r; asg_arc_t p, pmax;
    ma_hit_t_alloc *s = &(uopt->sources[v>>1]); memset(&pmax, 0, sizeof(pmax)); pmax.ol = 0;
    for (k = 0; k < s->length; k++) {
        qn = Get_qn(s->buffer[k]); tn = Get_tn(s->buffer[k]);
        r = ma_hit2arc(&(s->buffer[k]), sg->seq[qn].len, sg->seq[tn].len, uopt->max_hang, asm_opt.max_hang_rate, uopt->min_ovlp, &p);
        if(r < 0) continue;
        if((p.ul>>32) != v) continue;
        if(p.ol < min_ovlp) continue;
        if(p.ol > pmax.ol) pmax = p;
    }
    if(pmax.ol == 0) return 0;

    tn = pmax.v>>1;
    if((!(sg->seq[tn].del)) || (ff[tn])) return -1;///not a tip
    get_R_to_U(uopt->ruIndex, tn, &cc, &ccu);
    if(ccu == 1) return -1;
    if((cc != (uint32_t)-1) && 
                    ((!(sg->seq[cc].del)) || (ff[cc]))) {
        return -1;///contained in an existing read 
    }

    if((*path_len) >= sg->seq[v>>1].len) (*path_len) -= sg->seq[v>>1].len;
    else (*path_len) = 0;
    (*path_len) += (uint32_t)(pmax.ul) + sg->seq[pmax.v>>1].len;
    
    uint32_t v0 = v; (*res) = pmax.v;
    v = (pmax.v^1); s = &(uopt->sources[v>>1]); pmax.ol = 0;
    for (k = 0; k < s->length; k++) {
        qn = Get_qn(s->buffer[k]); tn = Get_tn(s->buffer[k]);
        r = ma_hit2arc(&(s->buffer[k]), sg->seq[qn].len, sg->seq[tn].len, uopt->max_hang, asm_opt.max_hang_rate, uopt->min_ovlp, &p);
        if(r < 0) continue;
        if((p.ul>>32) != v) continue;
        if(p.ol < min_ovlp) continue;
        if(p.ol > pmax.ol) pmax = p;
    }
    if((pmax.v^1) == v0) {
        ff[(*res)>>1] = 1; return 1;
    }
    (*res) = (uint64_t)-1; return -1;///not the longest overlap
}

void ug_ext_gfa(ug_opt_t *uopt, asg_t *sg, uint32_t max_len)
{
    asg_arc_t *av; asg64_v res, idx; kv_init(res); kv_init(idx);
    uint32_t v, nv, k, z, nvtx = sg->n_seq<<1, tip_n = uopt->tipsLen + 1, s, e;
    uint8_t *ff; CALLOC(ff, sg->n_seq); int32_t rr; uint64_t v0, w, bn, plen;
    asg_arc_t t, *p;
    for (v = 0; v < nvtx; v++) {
        if(sg->seq[v>>1].del) continue;
        av = asg_arc_a(sg, v); nv = asg_arc_n(sg, v);
        for (k = 0; k < nv && av[k].del; k++);
        if(k < nv) continue;
        for (z = 0, v0 = v, rr = 0, bn = res.n, plen = sg->seq[v0>>1].len; z < tip_n || plen < max_len; z++) {
            rr = gen_ext_tip(v0, sg, uopt, ff, 2000, &w, &plen);
            if(rr <= 0) break;
            kv_push(uint64_t, res, (v0<<32)|w); v0 = w;
        }
        if(rr < 0 || (z >= tip_n && plen >= max_len)) {
            for (z = bn; z < res.n; z++) {
                ff[((uint32_t)res.a[z])>>1] = 0;
            }
            res.n = bn; 
        } 
        if(res.n > bn) {
            bn <<= 32; bn |= (uint64_t)res.n; 
            kv_push(uint64_t, idx, bn);
        }
    }

    if(idx.n > 0) {
        for (k = 0; k < idx.n; k++) {
            s = idx.a[k]>>32; e = (uint32_t)idx.a[k];
            for (z = s; z < e; z++) {
                v = res.a[z]>>32; w = (uint32_t)res.a[z];
                assert((!sg->seq[v>>1].del) && (sg->seq[w>>1].del));
                sg->seq[w>>1].del = 0;
                rr = gen_spec_edge(sg, uopt, v, w, &t); assert(rr >= 0);
                p = asg_arc_pushp(sg); *p = t;
                rr = gen_spec_edge(sg, uopt, w^1, v^1, &t); assert(rr >= 0);
                p = asg_arc_pushp(sg); *p = t;
            }
        }
        
        free(sg->idx); sg->idx = 0; sg->is_srt = 0; asg_cleanup(sg);
        fprintf(stderr, "[M::%s::] # tips::%u\n", __func__, (uint32_t)idx.n);
    }

    free(ff); kv_destroy(res); kv_destroy(idx);
}


bubble_type *gen_bubble_chain(asg_t *sg, ma_ug_t *ug, ug_opt_t *uopt, uint8_t **ir_het, uint8_t avoid_het)
{
    kvec_asg_arc_t_warp new_rtg_edges;
    kv_init(new_rtg_edges.a);
    hap_cov_t *cov = NULL; 
    bubble_type *bub = NULL; 

    asg_t *copy_sg = copy_read_graph(sg);
    ma_ug_t *copy_ug = copy_untig_graph(ug);   
    // fprintf(stderr, "0[M::%s]\n", __func__);
    adjust_utg_by_primary(&copy_ug, copy_sg, TRIO_THRES, uopt->sources, uopt->reverse_sources, uopt->coverage_cut, 
    uopt->tipsLen, uopt->tip_drop_ratio, uopt->stops_threshold, uopt->ruIndex, uopt->chimeric_rate, uopt->drop_ratio,
    uopt->max_hang, uopt->min_ovlp, &new_rtg_edges, &cov, uopt->b_mask_t, 0, 0);
    ma_ug_destroy(copy_ug); copy_ug = NULL;
    asg_destroy(copy_sg); copy_sg = NULL; 
    // fprintf(stderr, "1[M::%s]\n", __func__);
    CALLOC(bub, 1); (*ir_het) = cov->t_ch->ir_het; 
    cov->t_ch->ir_het = NULL; cov->is_r_het = NULL;
    if(!avoid_het) {
        identify_bubbles(ug, bub, (*ir_het), NULL);
    } else {
        uint8_t *pr_het; CALLOC(pr_het, sg->n_seq);
        identify_bubbles_recal_poy(sg, ug, bub, pr_het, uopt->sources, uopt->ruIndex, NULL);
        free(pr_het);
    }
    
    // fprintf(stderr, "2[M::%s]\n", __func__);
    kv_destroy(new_rtg_edges.a); destory_hap_cov_t(&cov);
    if(asm_opt.purge_level_primary == 0) {///all nodes are het
        uint32_t k; 
        for (k = 0; k < ug->g->n_seq; k++) {
            if(IF_HOM(k, *bub)) bub->index[k] = bub->f_bub+1;
        }
    }
    // fprintf(stderr, "3[M::%s]\n", __func__);
    return bub;
}

void clear_path_dp_t(path_dp_t *x, asg_t *g)
{
    uint32_t n_vx = g->n_seq<<1;
    x->ref.n = x->pat.n = x->pat_cor.n = 0; 
    kv_resize(uint8_t, x->g_flt, n_vx); x->g_flt.n = n_vx;
    memset(x->g_flt.a, 0, sizeof(*(x->g_flt.a))*x->g_flt.n);
}

void clear_ubuf_t(ubuf_t *x, asg_t *g, all_ul_t *ul_idx, int32_t up_dp)
{
    uint32_t n_vx = g->n_seq<<1;
    x->a.n = x->S.n = x->T.n = x->b.n = x->e.n = 0;
    kv_resize(uinfo_t, x->a, n_vx); x->a.n = n_vx; memset(x->a.a, 0, sizeof(*(x->a.a))*x->a.n);
    if(up_dp) clear_path_dp_t(&(x->dp), g);
}

uint64_t get_ul_read_weight(all_ul_t *ul, uint32_t *prg, uinfo_t *g_idx, ul_vec_t *p, uint32_t ii, uint32_t v, uint32_t w)
{
    assert((!p->bb.a[ii].base)&&(p->bb.a[ii].hid == (v>>1))&&(p->bb.a[ii].el)&&(p->bb.a[ii].pchain));
    if(p->bb.a[ii].aidx == (uint32_t)-1) return 0; ///not connected
    uc_block_t *li = NULL, *lk = NULL; uint32_t li_v, lk_v;
    li = &(p->bb.a[ii]); li_v = (((uint32_t)(li->hid))<<1)|((uint32_t)(li->rev)); //li_v^=1;
    if(li_v != v) return 0;

    for (lk = &(p->bb.a[li->aidx]), w^=1; lk; ) {
        lk_v = (((uint32_t)(lk->hid))<<1)|((uint32_t)(lk->rev)); //lk_v^=1;
        while (lk_v != (w^1)) {
            if(g_idx[w].p==(uint32_t)-1) break;
            w = g_idx[w].p;
        }

        if(lk_v == (w^1)) {

        }

        lk = ((lk->aidx==(uint32_t)-1)?NULL:&(p->bb.a[lk->aidx]));
        
    }
    return 1;
}

#define arc_first(g, v) ((g)->arc[(g)->idx[(v)]>>32])
#define arc_cnt(g, v) ((uint32_t)(g)->idx[(v)])

uint64_t ulg_len_check(asg_t *g, uint32_t v, uc_block_t *p)
{
    // if(!(((v&1) && (p->ts==0)) || (((v&1)==0) && (p->te==g->seq[v>>1].len)))) {
    //     fprintf(stderr, "\n[M::%s::] v>>1:%u, v&1:%u, ts:%u, te:%u, tlen:%u, qs:%u, qe:%u, pidx:%u, aidx:%u, flag:%u\n", 
    //     __func__, v>>1, v&1, p->ts, p->te, g->seq[v>>1].len, p->qs, p->qe, p->pidx, p->aidx, flag);
    // }
    // assert(((v&1) && (p->ts==0)) || (((v&1)==0) && (p->te==g->seq[v>>1].len)));
    if(!(((v&1) && (p->ts==0)) || (((v&1)==0) && (p->te==g->seq[v>>1].len)))) return 0;
    // if((v&1) && (p->ts!=0)) return 0;
    // if(((v&1)==0) && (p->te!=g->seq[v>>1].len)) return 0;
    if(p->ts==0 && p->te==g->seq[v>>1].len) return 1;
    int32_t ol = arc_first(g, v).ol;///max len
    int32_t tl = p->te - p->ts;
    tl -= ol;
    if(tl > 20000 || tl > (g->seq[v>>1].len*0.2)) return 1;
    return 0;
}   

void update_l_coord(asg_t *g, uint64_t o_s, uint32_t v, uc_block_t *p, uint64_t *s, uint64_t *e)
{
    // if(!(((v&1) && (p->ts==0)) || (((v&1)==0) && (p->te==g->seq[v>>1].len)))) {
    //     fprintf(stderr, "\n[M::%s::] v>>1:%u, v&1:%u, ts:%u, te:%u, tlen:%u\n", 
    //     __func__, v>>1, v&1, p->ts, p->te, g->seq[v>>1].len);
    // }
    // assert(((v&1) && (p->ts==0)) || (((v&1)==0) && (p->te==g->seq[v>>1].len)));
    // assert((v&1) && (p->ts==0));
    // assert(((v&1)==0) && (p->te==g->seq[v>>1].len));
    if(((v&1)==0)) {
        (*s) = o_s + p->ts; (*e) = o_s + p->te;
    } else {
        (*s) = o_s + g->seq[v>>1].len - p->te;
        (*e) = o_s + g->seq[v>>1].len - p->ts;
    }
}

int64_t gen_ul_pat_seq(uc_block_t *a, path_dp_t *b, all_ul_t *ul, uint64_t idx, asg_t *g, uint32_t v, uint32_t is_backward, uint64_t *rl)
{
    uint32_t m = 0, pv, ai, pi; uint64_t *t, l = 0, tt, s, e; b->pat.n = 0; b->pat_cor.n = 0;
    if(is_backward) {
        for (ai = idx, l = 0; ai != (uint32_t)-1; ai = a[ai].pidx) {
            pv = (((uint32_t)(a[ai].hid))<<1)|((uint32_t)(a[ai].rev)); pv^=1;
            if((a[ai].pidx != (uint32_t)-1) && (!ulg_len_check(g, pv, &(a[ai])))) break;
            if((b->pat.n > 0) && (!ulg_len_check(g, pv^1, &(a[ai])))) break;
            if(pv == v) m++;
            kv_pushp(uint64_t, b->pat, &t); 
            update_l_coord(g, l, pv, &(a[ai]), &s, &e);
            kv_push(uint64_t, b->pat_cor, ((s<<32)|e));
            (*t) = pv; (*t) |= (l<<32);

            if(a[ai].pidx != (uint32_t)-1) l += a[ai].pdis;
            else l += g->seq[pv>>1].len;
        }
    } else {
        for (ai = idx, pi = 0; ai != (uint32_t)-1; ai = a[ai].aidx) {
            pv = (((uint32_t)(a[ai].hid))<<1)|((uint32_t)(a[ai].rev));
            if((a[ai].aidx != (uint32_t)-1) && (!ulg_len_check(g, pv, &(a[ai])))) break;
            if((pi > 0) && (!ulg_len_check(g, pv^1, &(a[ai])))) break;
            l = ai; pi++;
        }

        ai = l;
        kv_resize(uint64_t, b->pat, pi); b->pat.n = pi; 
        kv_resize(uint64_t, b->pat_cor, pi); b->pat_cor.n = pi; 
        for (l = 0; ai != (uint32_t)-1 && ai >= idx; ai = a[ai].pidx) {
            pv = (((uint32_t)(a[ai].hid))<<1)|((uint32_t)(a[ai].rev));
            pi--; if(pv == v) m++;
            b->pat.a[pi] = pv; b->pat.a[pi] |= (l<<32); 
            update_l_coord(g, l, pv, &(a[ai]), &s, &e);
            b->pat_cor.a[pi] = ((s<<32)|e);
            if(a[ai].pidx != (uint32_t)-1 && a[ai].pidx >= idx) l += a[ai].pdis;
            else l += g->seq[pv>>1].len;
        }
        assert(pi == 0);
        
        for (ai = 0; ai < b->pat.n; ai++) {
            tt = (b->pat.a[ai]>>32) + g->seq[((uint32_t)b->pat.a[ai])>>1].len;
            // if(l < tt) {
            //     fprintf(stderr, "\n[M::%s::] ai::%u, b->pat.n::%u, l::%lu, tt::%lu, beg::%lu\n", 
            //     __func__, ai, (uint32_t)b->pat.n, l, tt, (b->pat.a[ai]>>32));
            // }
            assert(l >= tt);
            tt = l - tt;
            b->pat.a[ai] = (uint32_t)b->pat.a[ai];
            b->pat.a[ai] += (tt<<32);
        }
    }
    (*rl) = l;

    assert(m > 0);
    return m;
}

uint32_t get_arch_len(asg_t *g, uint32_t v, uint32_t w)
{

    uint32_t i, an; asg_arc_t *av; 
    av = asg_arc_a(g, v); an = asg_arc_n(g, v);     
    for (i = 0; i < an; i++) {
        if(av[i].del) continue;
        if(av[i].v == w) return ((uint32_t)av[i].ul);
    }
    return (uint32_t)-1;
}

void gen_ref_pat_seq(asg_t *g, ubuf_t *b, uint64_t rul, double diff_rate, uint32_t v, uint32_t w/**, uint32_t is_debug**/)
{
    int64_t ref_n = b->dp.ref.n, k; uint64_t x, l; uint32_t p, arc_l;
    if(ref_n >= 2 && ((uint32_t)b->dp.ref.a[0]) == v && ((uint32_t)b->dp.ref.a[1]) == w) {
        // if(is_debug) fprintf(stderr, "+[M::%s::] ref_n::%ld, rul::%lu\n", __func__, ref_n, rul);
        if((ref_n) > 0 && ((b->dp.ref.a[ref_n-1]>>32) >= (rul*(1.0+diff_rate)))) {
            for (k = ref_n-1; k >= 0; k--) {
                if((b->dp.ref.a[k]>>32) < (rul*(1.0+diff_rate))) break;
                b->dp.g_flt.a[(uint32_t)b->dp.ref.a[k]] = 0;
            }
            ref_n = k + 1; b->dp.ref.n = ref_n;
        } else {
            p = (uint32_t)b->dp.ref.a[ref_n-1]; p ^= 1; p = b->a.a[p].p;
            for (l = b->dp.ref.a[ref_n-1]>>32; p != (uint32_t)-1; p = b->a.a[p].p) {
                x = p^1; arc_l = get_arch_len(g, ((uint32_t)b->dp.ref.a[b->dp.ref.n-1]), x); 
                assert(arc_l != (uint32_t)-1); l += arc_l; 
                if(l >= (rul*(1.0+diff_rate))) break;
                b->dp.g_flt.a[x] = 1; x += (l<<32); kv_push(uint64_t, b->dp.ref, x); 
            }
            ref_n = b->dp.ref.n;
        }
    } else {
        // if(is_debug) fprintf(stderr, "-[M::%s::] ref_n::%ld, rul::%lu\n", __func__, ref_n, rul);
        b->dp.ref.n = 0;
        for (k = 0; k < ref_n; k++) {
            b->dp.g_flt.a[(uint32_t)b->dp.ref.a[k]] = 0;
        }

        x = v; kv_push(uint64_t, b->dp.ref, x); b->dp.g_flt.a[x] = 1;
        for (p = w^1, l = 0; p != (uint32_t)-1; p = b->a.a[p].p) {
            x = p^1; arc_l = get_arch_len(g, ((uint32_t)b->dp.ref.a[b->dp.ref.n-1]), x); 
            assert(arc_l != (uint32_t)-1); l += arc_l; 
            if(l >= (rul*(1.0+diff_rate))) break;
            b->dp.g_flt.a[x] = 1; x += (l<<32); kv_push(uint64_t, b->dp.ref, x); 
        }
        ref_n = b->dp.ref.n;
    }
}

uint32_t quick_check(uc_block_t *a, asg_t *g, path_dp_t *b, double diff_rate, double filter_rate)
{
    uint64_t ref_l = 0, pat_l = 0, mm, min, max, l, k, s, e, match_s, match_e, unmatch_s, unmatch_e, match_l, unmatch_l;
    if(b->ref.n == 0 || b->pat.n == 0) return 0;
    assert(b->g_flt.a[(uint32_t)b->pat.a[0]]);
    for (k = 1; k < b->pat.n; k++) {///k = 0, must be matched
        if(b->g_flt.a[(uint32_t)b->pat.a[k]]) break;
    }
    if(k >= b->pat.n) return 0;

    ref_l = (b->ref.a[b->ref.n-1]>>32) + g->seq[((uint32_t)b->ref.a[b->ref.n-1])>>1].len;
    pat_l = (b->pat.a[b->pat.n-1]>>32) + g->seq[((uint32_t)b->pat.a[b->pat.n-1])>>1].len;
    mm = MIN(ref_l, pat_l); min = mm * (1.0 - diff_rate); max = mm * (1.0 + diff_rate);

    match_s = match_e = unmatch_s = unmatch_e = (uint64_t)-1; match_l = unmatch_l = 0;
    for (k = 0; k < b->pat.n; k++) {
        l = (b->pat.a[k]>>32) + g->seq[((uint32_t)b->pat.a[k])>>1].len; 
        s = b->pat_cor.a[k]>>32; e = (uint32_t)b->pat_cor.a[k];
        if(l <= min) {
            if(unmatch_e == (uint64_t)-1 || s >= unmatch_e) {
                if(unmatch_e != (uint64_t)-1) unmatch_l += unmatch_e - unmatch_s;
                unmatch_s = s; unmatch_e = e;
            } else {
                if(e > unmatch_e) unmatch_e = e;
            }
        }

        if(l <= max) {
            if(b->g_flt.a[(uint32_t)b->pat.a[k]]) {
                if(match_e == (uint64_t)-1 || s >= match_e) {
                    if(match_e != (uint64_t)-1) match_l += match_e - match_s;
                    match_s = s; match_e = e;
                } else {
                    if(e > match_e) match_e = e;
                }
            }
        }
        if(l > max) break;
    }
    
    if(unmatch_e != (uint64_t)-1) unmatch_l += unmatch_e - unmatch_s;
    if(match_e != (uint64_t)-1) match_l += match_e - match_s;
    // fprintf(stderr, "[M::%s::] min::%lu, max::%lu, match_l::%lu, unmatch_l::%lu\n", __func__, 
    //         min, max, match_l, unmatch_l);
    if(match_l >= unmatch_l*filter_rate) return 1;
    return 0;
}

#define ul_dp_idx(dp, x, y) ((dp).m*(x)+(y))
#define e_mdp 0
#define ue_mdp 1
#define lpat_dp 2
#define lref_dp 3
//0->match; 1->mismatch; 2->up (longer pat); 3->left (longer ref)
void init_ul_dp(asg_t *g, path_dp_t *dp)
{
    kv_resize(uint8_t, dp->m_dir, (dp->pat.n+1)*(dp->ref.n+1));
    kv_resize(int64_t, dp->m_score, (dp->pat.n+1)*(dp->ref.n+1));
    dp->n = dp->pat.n+1; dp->m = dp->ref.n+1;
    /**
    uint64_t k, s, e, sp, ep, sc;
    dp->m_dir.a[0] = dp->m_score.a[0] = 0; ///[0, 0]
    for (k = 1, sp = ep = (uint64_t)-1; k < dp->n; k++) {///pat; ul read
        s = dp->pat.a[k-1]>>32; e = (dp->pat.a[k-1]>>32) + g->seq[(uint32_t)dp->pat.a[k-1]].len; sc = 0;
        if(ep == (uint64_t)-1 || s >= ep) {
            sc = e - s; sp = s; ep = e;
        } else {
            if(e > ep) sc = e - ep;
        }
        if(((uint32_t)dp->pat.a[k-1]) == ((uint32_t)dp->ref.a[0])) {
            dp->m_score.a[ul_dp_idx(*dp, k, 0)] = dp->m_dir.a[ul_dp_idx(*dp, k, 0)] = 0;
        } else {
            dp->m_score.a[ul_dp_idx(*dp, k, 0)] = dp->m_score.a[ul_dp_idx(*dp, k-1, 0)] + sc;
            dp->m_dir.a[ul_dp_idx(*dp, k, 0)] = lpat_dp;
        }
    }

    for (k = 1, sp = ep = (uint64_t)-1; k < dp->m; k++) {///ref; graph
        s = dp->ref.a[k-1]>>32; e = (dp->ref.a[k-1]>>32) + g->seq[(uint32_t)dp->ref.a[k-1]].len; sc = 0;
        if(ep == (uint64_t)-1 || s >= ep) {
            sc = e - s; sp = s; ep = e;
        } else {
            if(e > ep) sc = e - ep;
        }
        dp->m_score.a[ul_dp_idx(*dp, 0, k)] = dp->m_score.a[ul_dp_idx(*dp, 0, k-1)] + sc;
        dp->m_dir.a[ul_dp_idx(*dp, 0, k)] = lref_dp;
    }
    **/
   uint64_t k; int64_t min_weight = -1*((int64_t)(0xffffffff)); dp->m_dir.a[0] = dp->m_score.a[0] = 0; ///[0, 0]

   for (k = 1; k < dp->n; k++) {///pat; ul read
        if(((uint32_t)dp->pat.a[k-1]) == ((uint32_t)dp->ref.a[0])) {
            dp->m_score.a[ul_dp_idx(*dp, k, 0)] = dp->m_dir.a[ul_dp_idx(*dp, k, 0)] = 0;
        } else {
            dp->m_score.a[ul_dp_idx(*dp, k, 0)] = min_weight;
                        // dp->m_score.a[ul_dp_idx(*dp, k-1, 0)] - g->seq[((uint32_t)dp->pat.a[k-1])>>1].len;
            dp->m_dir.a[ul_dp_idx(*dp, k, 0)] = lpat_dp;
        }
    }

    for (k = 1; k < dp->m; k++) {///ref; graph
        dp->m_score.a[ul_dp_idx(*dp, 0, k)] = min_weight;
                // dp->m_score.a[ul_dp_idx(*dp, 0, k-1)] - g->seq[((uint32_t)dp->ref.a[k-1])>>1].len;
        dp->m_dir.a[ul_dp_idx(*dp, 0, k)] = lref_dp;
    }
}

uint64_t node_check(uint32_t ul_v, uint32_t g_v, uint32_t ul_v_len, uint32_t g_v_len, double diff_len, uint32_t ul_weight)
{
    if(ul_v != g_v) return 0;
    uint32_t x = ((ul_v_len >= g_v_len)? (ul_v_len - g_v_len): (g_v_len - ul_v_len));
    if(x > g_v_len*diff_len) return 0;
    if(g_v_len == 0) {
        return ul_weight;
    } else {
        return (((double)(g_v_len-x))/((double)g_v_len))*ul_weight;
    }
}

void print_dp_matrix(path_dp_t *dp)
{
    uint64_t i, j;
    for (i = 0; i < dp->n; i++) {//pat
        for (j = 0; j < dp->m; j++) {//mat
            fprintf(stderr, "%ld<%u>,\t", dp->m_score.a[ul_dp_idx(*dp, i, j)], dp->m_dir.a[ul_dp_idx(*dp, i, j)]);
        }
        fprintf(stderr, "\n");
    }
    
}

int64_t ul_dp0(bubble_type* bub, asg_t *g, path_dp_t *dp, double pass_thres/**, uint64_t is_debug**/)
{
    if(dp->pat.n < 2 || dp->ref.n < 2) return -1;
    uint64_t i, j, d, w; int64_t sc0, sc1, sc2, sc, sc_i, sc_j;
    init_ul_dp(g, dp);
    // if(is_debug) {
    //     print_dp_matrix(dp);
    // }
    for (i = 1; i < dp->n; i++) {//pat
        for (j = 1; j < dp->m; j++) {//mat
            w = node_check(((uint32_t)dp->pat.a[i-1]), ((uint32_t)dp->ref.a[j-1]), 
                        dp->pat.a[i-1]>>32, dp->ref.a[j-1]>>32, 0.04, g->seq[((uint32_t)dp->pat.a[i-1])>>1].len);

            if(w) {
                // if(is_debug) fprintf(stderr, "\n[M::%s::] i::%lu, j::%lu, w::%lu\n", __func__, i, j, w);
                dp->m_score.a[ul_dp_idx(*dp, i, j)] = dp->m_score.a[ul_dp_idx(*dp, i-1, j-1)] + w;
                dp->m_dir.a[ul_dp_idx(*dp, i, j)] = e_mdp;
            } else {
                sc0 = dp->m_score.a[ul_dp_idx(*dp, i-1, j-1)] - 
                        (int64_t)(MIN(g->seq[((uint32_t)dp->pat.a[i-1])>>1].len, g->seq[((uint32_t)dp->ref.a[j-1])>>1].len));
                sc1 = dp->m_score.a[ul_dp_idx(*dp, i-1, j)] - (int64_t)(g->seq[((uint32_t)dp->ref.a[j-1])>>1].len);
                sc2 = dp->m_score.a[ul_dp_idx(*dp, i, j-1)] - (int64_t)(g->seq[((uint32_t)dp->pat.a[i-1])>>1].len);

                d = ue_mdp; sc = sc0;
                if(sc < sc1) {
                    sc = sc1; d = lref_dp;
                }
                if(sc < sc2) {
                    sc = sc2; d = lpat_dp;
                }
                dp->m_score.a[ul_dp_idx(*dp, i, j)] = sc;
                dp->m_dir.a[ul_dp_idx(*dp, i, j)] = d;
            }
        }
    }

    sc = 0; sc_i = sc_j = -1;
    for (j = 1, i = dp->n - 1; j < dp->m; j++) {
        if(sc_i < 0 || sc < dp->m_score.a[ul_dp_idx(*dp, i, j)]) {
            sc_i = i; sc_j = j; sc = dp->m_score.a[ul_dp_idx(*dp, i, j)];
        }
    }
    for (i = 1, j = dp->m - 1; i < dp->n; i++) {
        if(sc_i < 0 || sc < dp->m_score.a[ul_dp_idx(*dp, i, j)]) {
            sc_i = i; sc_j = j; sc = dp->m_score.a[ul_dp_idx(*dp, i, j)];
        }
    }

    // if(is_debug) {
    //     fprintf(stderr, "[M::%s::] sc_i::%ld, sc_j::%ld, sc::%ld\n", __func__, sc_i, sc_j, sc);
    //     // print_dp_matrix(dp);
    // }

    if (sc_i < 0 || sc_i == 1 || sc_j == 1) return 0;
    int64_t match_all = 0, unmatch_all = 0, het_match_all = 0, het_unmatch_all = 0, het_occ = 0;
    //pat_s = sc_i - 1; ref_s = ref_e = sc_j - 1;
    for(i = sc_i, j = sc_j; i > 0 && j > 0;) {
        d = dp->m_dir.a[ul_dp_idx(*dp, i, j)];
        //pat_s = i - 1; ref_s = j - 1;
        // if(is_debug) fprintf(stderr, "[M::%s::] i::%lu, j::%lu, dir::%lu\n", __func__, i, j, d);
        // if(pat_s == 0) continue;
        if(d == e_mdp) {
            // if(i != 1 || j != 1) {//skip (i == 1 && j == 1)
            match_all += g->seq[((uint32_t)dp->pat.a[i-1])>>1].len;
            if(!IF_HOM((((uint32_t)dp->pat.a[i-1])>>1), *bub)) {
                het_match_all += g->seq[((uint32_t)dp->pat.a[i-1])>>1].len;
                het_occ++;
            }
            i--; j--; 
            // if(is_debug) {
            //     fprintf(stderr, "[M::%s::] PAT::utg%.6dl(%c::len->%lu), REF::utg%.6dl(%c::len->%lu)\n", __func__, 
            //         ((((uint32_t)dp->pat.a[i])>>1))+1, "+-"[((uint32_t)dp->pat.a[i])&1], dp->pat.a[i]>>32,
            //         ((((uint32_t)dp->ref.a[j])>>1))+1, "+-"[((uint32_t)dp->ref.a[j])&1], dp->ref.a[j]>>32);
            // }
        }
        if(d == ue_mdp) {
            unmatch_all += g->seq[((uint32_t)dp->pat.a[i-1])>>1].len;
            if(!IF_HOM((((uint32_t)dp->pat.a[i-1])>>1), *bub)) {
                het_unmatch_all += g->seq[((uint32_t)dp->pat.a[i-1])>>1].len;
            }
            i--; j--;
        }
        if(d == lpat_dp) {
            j--;
        }
        if(d == lref_dp) {
            unmatch_all += g->seq[((uint32_t)dp->pat.a[i-1])>>1].len;
            if(!IF_HOM((((uint32_t)dp->pat.a[i-1])>>1), *bub)) {
                het_unmatch_all += g->seq[((uint32_t)dp->pat.a[i-1])>>1].len;
            }
            i--;
        }
    }
    // if(is_debug) {
    //     fprintf(stderr, "[M::%s::] het_occ::%ld, match_all::%ld, unmatch_all::%ld, het_match_all::%ld, het_unmatch_all::%ld\n", __func__, 
    //     het_occ, match_all, unmatch_all, het_match_all, het_unmatch_all);
    // }

    assert((j == 0) && (((uint32_t)dp->pat.a[i]) == ((uint32_t)dp->ref.a[j])));
    ///skip the beg node
    match_all -= g->seq[((uint32_t)dp->pat.a[i])>>1].len;
    if(!IF_HOM((((uint32_t)dp->pat.a[i])>>1), *bub)) {
        het_match_all -= g->seq[((uint32_t)dp->pat.a[i])>>1].len; het_occ--;
    }


    if(het_occ < 1 || match_all == 0 || het_match_all == 0) return -1;
    if(match_all <= (match_all + unmatch_all)*pass_thres) return -1;
    if(het_match_all <= (het_match_all + het_unmatch_all)*pass_thres) return -1;
    return (((double)het_match_all)/((double)(het_match_all + het_unmatch_all)))*
                            (((uint32_t)dp->pat_cor.a[i]) - (dp->pat_cor.a[i]>>32));
}
uint64_t ul_dp(bubble_type* bub, all_ul_t *ul, ubuf_t *b, uint64_t *a, int64_t a_n, asg_t *g, uint32_t v, uint32_t w, uint32_t is_backward)
{
    int64_t i, m, a_i = -1; uc_block_t *p; uint32_t gv; uint64_t url; int64_t we;
    for (i = m = 0; i < a_n; i++) {
        p = &(ul->a[a[i]>>32].bb.a[(uint32_t)(a[i])]);
        assert((!p->base)&&(p->hid == (v>>1))&&(p->el)/**&&(p->pchain)**/);
        if(!p->pchain) continue;
        gv = (((uint32_t)(p->hid))<<1)|((uint32_t)(p->rev));
        if(is_backward) {
            gv ^= 1;
            if(p->pidx == (uint32_t)-1) continue;
        } else {
            if(p->aidx == (uint32_t)-1) continue;
        }
        if(gv != v) continue;
        // fprintf(stderr, "[M::%s::] i:%ld, a_n:%lu, v:%u, w:%u, ulid:%lu, uidx:%u, is_backward:%u\n", 
        // __func__, i, a_n, v, w, a[i]>>32, (uint32_t)(a[i]), is_backward);
        // if(v == 6 && w == 5 && (a[i]>>32) == 95) {
        // int64_t z;
        //     for (z = 0; z < ul->a[a[i]>>32].bb.n; z++) {
        //         fprintf(stderr, "[M::%s::vid->%u] qs:%u, qe:%u, qlen:%u, ts:%u, te:%u, tlen:%u\n", __func__, 
        //         (((uint32_t)(ul->a[a[i]>>32].bb.a[z].hid))<<1)|((uint32_t)(ul->a[a[i]>>32].bb.a[z].rev)), 
        //         ul->a[a[i]>>32].bb.a[z].qs, ul->a[a[i]>>32].bb.a[z].qe, ul->a[a[i]>>32].rlen, 
        //         ul->a[a[i]>>32].bb.a[z].ts, ul->a[a[i]>>32].bb.a[z].te, g->seq[ul->a[a[i]>>32].bb.a[z].hid].len);
        //     }
        // }
        if(ulg_len_check(g, v, p) == 0) continue;///too short
        m++; a_i = i;
    }
    
    if(m == 0) return 0;
    // uint32_t is_debug = (v == 231);
    // if(is_debug) {
    //     fprintf(stderr, "\n+[M::%s::] m->%ld, ulid->%lu\n", __func__, m, a[a_i]>>32);
    // }

    p = &(ul->a[a[a_i]>>32].bb.a[(uint32_t)(a[a_i])]);//longest one
    m = gen_ul_pat_seq(ul->a[a[a_i]>>32].bb.a, &(b->dp), ul, (uint32_t)(a[a_i]), g, v, is_backward, &url);
    assert(b->dp.pat.n > 0);

    // if(is_debug) {
    //     for (i = 0; i < (int64_t)b->dp.pat.n; i++) {
    //         fprintf(stderr, "pat[M::%s::] utg%.6dl(%c), d::%lu\n", __func__, 
    //         (((uint32_t)b->dp.pat.a[i])>>1)+1, "+-"[((uint32_t)b->dp.pat.a[i])&1], b->dp.pat.a[i]>>32);
    //     }
    // }
    if(b->dp.pat.n == 1) return 0;
    gen_ref_pat_seq(g, b, url, UL_TRAV_HERATE, v, w/**, is_debug**/);

    // if(is_debug) {
    //     for (i = 0; i < (int64_t)b->dp.ref.n; i++) {
    //         fprintf(stderr, "ref[M::%s::] utg%.6dl(%c), d::%lu\n", __func__, 
    //         (((uint32_t)b->dp.ref.a[i])>>1)+1, "+-"[((uint32_t)b->dp.ref.a[i])&1], b->dp.ref.a[i]>>32);
    //     }
    // }

    if(!quick_check(ul->a[a[a_i]>>32].bb.a, g, &(b->dp), UL_TRAV_HERATE, UL_TRAV_FT_RATE)) return 0;
    we = ul_dp0(bub, g, &(b->dp), 0.9/**, is_debug**/);
    if(we < 0) return 0;
    return we;
}

uint64_t get_eul_weight(bubble_type* bub, uint32_t v, uint32_t w, uinfo_t *g_idx, ubuf_t *b, all_ul_t *ul, asg_t *g)
{
    uint64_t *a, to = 0; int64_t a_n, l, k;
    a = ul->ridx.occ.a + ul->ridx.idx.a[v>>1];
	a_n = ul->ridx.idx.a[(v>>1)+1] - ul->ridx.idx.a[v>>1]; 
    b->dp.pat.n = b->dp.pat_cor.n = b->dp.ref.n = 0;
    for (l = 0, k = 1; k <= a_n; k++) {
        if((k == a_n) || ((a[k]>>32) != (a[l]>>32))) {
            // if((v == 106 && w == 103) || (v == 24 && w == 23)) {
            //     fprintf(stderr, "+[M::%s::] l->%ld, k->%ld, a_n->%ld\n", __func__, l, k, a_n);
            // }
            to += ul_dp(bub, ul, b, a + l, k - l, g, v, w, 1);
            // if((v == 106 && w == 103) || (v == 24 && w == 23)) { 
            //     fprintf(stderr, "-[M::%s::] l->%ld, k->%ld, a_n->%ld\n", __func__, l, k, a_n);
            // }
            to += ul_dp(bub, ul, b, a + l, k - l, g, v, w, 0);
            // if((v == 106 && w == 103) || (v == 24 && w == 23)) {
            //     fprintf(stderr, "*[M::%s::] l->%ld, k->%ld, a_n->%ld\n", __func__, l, k, a_n);
            // }
            l = k;
        }
    }

    return to;
}

void extract_paths(ubuf_t *b, ma_utg_v *gu, ul_path_t *res, uint32_t src, uint32_t dest)
{
    uint32_t i, v; uinfo_t *t; uint64_t m = 0, mi = (uint64_t)-1, avn; ///kv_resize(uinfo_srt_t, b->srt, b->b.n); b->srt.n = b->b.n;
    for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
        t = &b->a.a[b->b.a[i]]; ///memset(t, 0, sizeof(*(t)));
        //b->srt.a[i].c = ((uint64_t)-1) - t->c; b->srt.a[i].i = i;
        if(m < t->c) {
            m = t->c; mi = i; ///b->b.a[i];
        }
    }
    //radix_sort_uinfo_srt_t_c(b->srt.a, b->srt.a + b->srt.n);
    // fprintf(stderr, "\n[M::%s::]->beg\n", __func__);
    if(mi != (uint64_t)-1) {
        ma_utg_t *p; kv_pushp(ma_utg_t, res->buf_ug->u, &p); 
        p->len = p->circ = p->n = 0; p->start = p->end = UINT32_MAX; p->a = NULL; p->s = NULL;
        for(v = b->b.a[mi], avn = 0; v != (uint32_t)-1; v = b->a.a[v].p) {
            if(b->a.a[v].c == 0 && avn == 0) break;
            b->us.a[v>>1] = 1; 
            if(v != src && v != dest) kv_push(uint64_t, (*p), v); 
            // fprintf(stderr, "[M::%s::] utg%.6d%c(%c)\tC::%lu\n", __func__, (v>>1)+1, "lc"[gu->a[v>>1].circ], "+-"[v&1], b->a.a[v].c);
            avn = b->a.a[v].c;
        }

        if(p->n == 0) res->buf_ug->u.n = 0;
    }
}

uint32_t hc_simple_traversal(bubble_type* bub, asg_t *g, ma_utg_v *gu, ubuf_t *b, all_ul_t *ul, uint32_t src, uint32_t dest, ul_path_t *res)
{
    uint32_t v, nv, i, w, n_pending = 0, is_update = 0; uint64_t l, d, c, cc, nc, c_nc; asg_arc_t *av; uinfo_t *t;
    if (g->seq[src>>1].del || g->seq[dest>>1].del) return 0;
    clear_ubuf_t(b, g, ul, 1);  b->a.a[src].p = (uint32_t)-1;
    kv_push(uint32_t, b->S, src);
    while (b->S.n > 0) {
        v = kv_pop(b->S); d = b->a.a[v].d; c = b->a.a[v].c; nc = b->a.a[v].nc;
        nv = asg_arc_n(g, v); av = asg_arc_a(g, v);
        for (i = 0; i < nv; ++i) {
            if (av[i].del || g->seq[av[i].v>>1].del) continue;
            w = av[i].v; l = (uint32_t)av[i].ul; t = &b->a.a[w];
            kv_push(uint32_t, b->e, ((g->idx[v]>>32)+i)); ///push the edge
            // fprintf(stderr, "+[M::%s::] utg%.6d%c(%c:%u) -> utg%.6d%c(%c:%u)\n", __func__, 
            // (v>>1)+1, "lc"[gu->a[v>>1].circ], "+-"[v&1], v, 
            // (w>>1)+1, "lc"[gu->a[w>>1].circ], "+-"[w&1], w);
            cc = c + get_eul_weight(bub, w^1, v^1, b->a.a, b, ul, g);
            c_nc = nc + (gu?gu->a[w>>1].n:1);
            // fprintf(stderr, "-[M::%s::] utg%.6d%c(%c:%u) -> utg%.6d%c(%c:%u)\n", __func__, 
            // (v>>1)+1, "lc"[gu->a[v>>1].circ], "+-"[v&1], v, 
            // (w>>1)+1, "lc"[gu->a[w>>1].circ], "+-"[w&1], w);

            if (t->s == 0) {///a new node
                kv_push(uint32_t, b->b, w); // save it for revert
                t->p = v; t->s = 1; t->d = d + l;
                t->r = get_arcs(g, w^1, NULL, 0);
                t->nc = c_nc; t->c = cc;
                ++n_pending; 
            } else {
                is_update = 0;
                if(b->us.a[t->p>>1] == b->us.a[v>>1]) {
                    if(cc > t->c) is_update = 1;
                    if(cc == t->c && c_nc > t->nc) is_update = 1;
                } else {
                    if(b->us.a[t->p>>1]) is_update = 1;
                    else is_update = 0;
                }

                if(is_update) {
                    t->p = v; t->s = 1; t->d = d + l; t->nc = c_nc; t->c = cc;
                }
            }

            if (--(t->r) == 0) {
                if(get_arcs(g, w, NULL, 0) > 0) kv_push(uint32_t, b->S, w);
                --n_pending;
                if(w == dest && n_pending == 0) goto pp_end;
            }
        }
    }
    pp_end:
    extract_paths(b, gu, res, src, dest);
    for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
        t = &b->a.a[b->b.a[i]]; memset(t, 0, sizeof(*(t)));
    }
    return 1;
}


// void clean_path_g(ul_path_t *p, asg_t *ref)
// {
//     if(!(p->pg)) CALLOC(p->pg, 1); 
//     free(p->pg->idx); p->pg->idx = 0; p->pg->is_srt = 0; p->pg->is_symm = 0;
//     REALLOC(p->pg->seq, ref->n_seq); p->pg->n_seq = ref->n_seq; p->pg->n_arc = 0;
//     memcpy(p->pg->seq, ref->seq, p->pg->n_seq*sizeof(*(p->pg->seq)));
//     asg_cleanup(p->pg);
// }

void build_sub_graph(ul_resolve_t *uidx, ul_path_t *res)
{

}

uint64_t get_chain_ul_cov(ul_resolve_t *uidx, uint64_t *a, uint64_t a_n, uint64_t bub_id, uint64_t inner_beg)
{
    if(a_n <= 0) return 0;
    uint64_t w = 0, bid, rev; uint32_t rt, beg, sink; 
    uidx->path.buf.n = 0;
    // kv_pushp(ul_sub_path_t, uidx->path, &p); memset(p, 0, sizeof((*p)));
    // p->bid = bub_id; p->beg = inner_beg; p->occ = a_n;

    bid = a[0]>>33; rev = (a[0]>>32)&1; 
    get_bubbles(uidx->bub, bid, rev==0?&rt:NULL, rev==1?&rt:NULL, NULL, NULL, NULL); beg = rt;

    bid = a[a_n-1]>>33; rev = (a[a_n-1]>>32)&1; 
    get_bubbles(uidx->bub, bid, rev==1?&rt:NULL, rev==0?&rt:NULL, NULL, NULL, NULL); sink = rt^1;








    // fprintf(stderr, "\n[M::%s::] # bubbles->%lu\n", __func__, a_n);
    // uint64_t k
    // for (k = 0; k < a_n; k++) {
    //     bid = a[k]>>33; rev = (a[k]>>32)&1;
    //     get_bubbles(uidx->bub, bid, rev==0?&rt:NULL, rev==1?&rt:NULL, NULL, NULL, NULL);
    //     // w += get_bub_ul_cov(a[k]>>33, ug, idx);
    //     fprintf(stderr, "[M::%s::] utg%.6d%c(%c)\n", __func__, (rt>>1)+1, "lc"[uidx->l1_ug->u.a[rt>>1].circ], "+-"[rt&1]);
    // }
    // if (a_n > 0) {
    //     bid = a[k-1]>>33; rev = (a[k-1]>>32)&1;
    //     get_bubbles(uidx->bub, bid, rev==1?&rt:NULL, rev==0?&rt:NULL, NULL, NULL, NULL); rt ^= 1;
    //     fprintf(stderr, "[M::%s::] utg%.6d%c(%c)\n", __func__, (rt>>1)+1, "lc"[uidx->l1_ug->u.a[rt>>1].circ], "+-"[rt&1]);
    // }
    // clean_path_g(&(uidx->path), uidx->l1_ug->g);
    ma_utg_t *p; 
    kv_resize(uint8_t, uidx->buf.us, uidx->l1_ug->g->n_seq); uidx->buf.us.n = uidx->l1_ug->g->n_seq; memset(uidx->buf.us.a, 0, uidx->buf.us.n*sizeof(*(uidx->buf.us.a)));
    hc_simple_traversal(uidx->bub, uidx->l1_ug->g, &(uidx->l1_ug->u), &(uidx->buf), uidx->idx, beg, sink, &(uidx->path));
    hc_simple_traversal(uidx->bub, uidx->l1_ug->g, &(uidx->l1_ug->u), &(uidx->buf), uidx->idx, beg, sink, &(uidx->path));
    
    kv_pushp(ma_utg_t, uidx->path.buf_ug->u, &p); p->len = p->circ = p->n = 0; 
    p->start = p->end = UINT32_MAX; p->a = NULL; p->s = NULL; kv_push(uint64_t, *p, beg);

    kv_pushp(ma_utg_t, uidx->path.buf_ug->u, &p); p->len = p->circ = p->n = 0; 
    p->start = p->end = UINT32_MAX; p->a = NULL; p->s = NULL; kv_push(uint64_t, *p, sink);

    return w;
}

void phrase_exact_chains(ul_resolve_t *uidx, uint64_t *a, int64_t a_n, asg_t *bg, uint64_t bub_id, uint64_t inner_beg)
{
    int64_t k, l;
    for (k = l = 0; k < a_n; k++) {
        if(k+1 >= a_n) continue;
        if((arc_first(bg, a[k]>>32)).el) continue;
        get_chain_ul_cov(uidx, a+l, k+1-l, bub_id, inner_beg + l);
        l = k + 1;
    }
    if(l < a_n) get_chain_ul_cov(uidx, a+l, a_n-l, bub_id, inner_beg + l);
}
void resolve_dip_bub_chains(ul_resolve_t *uidx)
{
    bubble_type *bub = uidx->bub; 
    uint32_t i; int32_t k, l, un; ma_utg_t *u; 
    for (i = 0; i < bub->b_ug->u.n; i++) {
        u = &(bub->b_ug->u.a[i]); un = u->n;
        if(un <= 1) continue;///single bubble
        for (l = -1, k = 0; k <= un; k++) {
            if((k == un) || ((u->a[k]>>33) >= bub->f_bub)) {///not a complete bubble
                // if(k < un) {
                //     fprintf(stderr, "++++++[M::%s::] l->%d, k->%d, bid->%lu, f_bub->%lu\n", 
                //     __func__, l, k, (u->a[k]>>33), bub->f_bub);
                // }
                if(k - l > 1) {///at least a complete bubble
                    // fprintf(stderr, "\n[M::%s::] l->%d, k->%d\n", __func__, l, k);
                    phrase_exact_chains(uidx, u->a+l+1, k-l-1, bub->b_g, i, l+1);
                }
                l = k;
            }
        }

    }
}

/**
static void fill_ul_path_srt_t(void *data, long i, int tid) // callback for kt_for()
{
    all_ul_t *idx = ((ul_resolve_t *)data)->idx; 
    ul_path_srt_t *idx_srt = &(((ul_resolve_t *)data)->psrt);
	uc_block_t *a = NULL; uc_block_t *p; int64_t k, a_n, b_n, occ; uint64_t *b = NULL;

	a = idx->a[i].bb.a; a_n = idx->a[i].bb.n;
    if(idx_srt->ul_n == 0) {
        for (k = occ = 0; k < a_n; k++) {
            p = &(a[k]);
            if(p->base || (!p->el) || (!p->pchain)) continue;
            occ++;
        }
        if(occ == 1) occ = 0;///UL read is too short
        idx_srt->idx.a[i] = occ;
    } else {
        b_n = idx_srt->idx.a[i]>>33; b = idx_srt->srt.a + (uint32_t)idx_srt->idx.a[i];
        if(b_n == 0) return;//UL covers 0 or 1 node; not useful
        for (k = occ = 0; k < a_n; k++) {
            p = &(a[k]);
            if(p->base || (!p->el) || (!p->pchain)) continue;
            b[occ] = p->hid; b[occ] <<= 1; b[occ] |= p->rev; b[occ] <<= 32; b[occ] += (uint64_t)k;
            occ++;
        }
        assert(occ == b_n);
        radix_sort_srt64(b, b + b_n);
        for (k = 1; k < b_n; k++) {
            if((b[k]>>32) == (b[k-1]>>32)) break;
        }

        if(k < b_n) idx_srt->idx.a[i] |= (uint64_t)(0x100000000);
    }
}


void init_ul_path_srt_t(ul_resolve_t *p, all_ul_t *idx, ul_path_srt_t *idx_srt)
{
    idx_srt->ul_n = 0;
    idx_srt->idx.n = idx_srt->idx.m = idx->n; CALLOC(idx_srt->idx.a, idx_srt->idx.n);

    kt_for(asm_opt.thread_num, fill_ul_path_srt_t, p, idx->n);
    uint64_t l, k;
    for (k = l = 0; k < idx_srt->idx.n; k++) {
        idx_srt->idx.a[k] <<= 33; idx_srt->idx.a[k] += l; l += (idx_srt->idx.a[k]>>33);
    }
    MALLOC(idx_srt->srt.a, l); idx_srt->srt.n = idx_srt->srt.m = l;

    idx_srt->ul_n = idx->n;
    kt_for(asm_opt.thread_num, fill_ul_path_srt_t, p, idx->n);
}
**/

uint64_t ug_occ_w(uint64_t is, uint64_t ie, ma_utg_t *u)
{
    if(is == 0 && ie == u->len) return u->n;
    uint64_t l, i, us, ue, occ;
    for (i = l = occ = 0; i < u->n; i++) {
		us = l; ue = l + Get_READ_LENGTH(R_INF, (u->a[i]>>33));
        // if(is == 15390 && ie == 31730) {
        //     fprintf(stderr, "[M::%s::i->%lu] is->%lu, ie->%lu, us->%lu, ue->%lu, u->len->%u\n", 
        //     __func__, i, is, ie, us, ue, u->len);
        // }
		if(is <= us && ie >= ue) occ++;
		if(us >= ie) break;
		l += (uint32_t)u->a[i];
	}
    return occ;
}

// static void gen_ul_str_idx_t(void *data, long i, int tid) // callback for kt_for()
// {
//     all_ul_t *idx = ((ul_resolve_t *)data)->idx; 
//     ma_ug_t *ug = ((ul_resolve_t *)data)->l1_ug;
//     ul_str_t *str = &(((ul_resolve_t *)data)->pstr.str.a[i]);
//     // uint64_t *str_idx = ((ul_resolve_t *)data)->pstr.idx.a;
// 	   uc_block_t *a = NULL; uc_block_t *xk; uint64_t t;
//     uint64_t k, a_n; 

//     str->n = str->m = 0; str->a = NULL; 
// 	   a = idx->a[i].bb.a; a_n = idx->a[i].bb.n;
//     kv_resize(uint64_t, *str, idx->a[i].bb.n);
//     for (k = 0; k < a_n; k++) {
//         xk = &(a[k]);
//         if(xk->base || (!xk->el) || (!xk->pchain)) continue;
//         t = k; t <<= 32; t += (xk->hid<<1); t += xk->rev; 
//         kv_push(uint64_t, *str, t);
//     }
//     str->cn = str->n; str->is_cir = 0;
// }

void init_ul_str_idx_t(ul_resolve_t *p)
{
    uint64_t k, l, i, m, *a, a_n, x_n; uc_block_t *x_a; 
    all_ul_t *idx = p->idx; ul_str_idx_t *str = &(p->pstr); ul_str_t *z; uint32_t v0, v1;
    MALLOC(str->str.a, idx->n); str->str.n = str->str.m = idx->n;
    CALLOC(str->idx.a, p->l1_ug->u.n+1); str->idx.n = str->idx.m = p->l1_ug->u.n+1;
    for (i = 0; i < str->str.n; i++) {
        x_a = idx->a[i].bb.a; x_n = idx->a[i].bb.n; 
        memset(&(str->str.a[i]), 0, sizeof(str->str.a[i]));
        for (k = 0; k < x_n; k++) {
            if(x_a[k].base || (!x_a[k].el) || (!x_a[k].pchain)) continue;
            l = k; l <<= 32; l += (x_a[k].hid<<1); l += x_a[k].rev; str->idx.a[x_a[k].hid]++;
            kv_push(uint64_t, str->str.a[i], l);///idx|hid|rev
        }
        str->str.a[i].cn = str->str.a[i].n; str->str.a[i].is_cir = 0;
    }
    
    // kt_for(asm_opt.thread_num, gen_ul_str_idx_t, p, str->str.n);
    for (k = l = 0; k < str->idx.n; k++) {
		m = str->idx.a[k];
		str->idx.a[k] = l;
		l += m;
	}

    MALLOC(str->occ.a, l); str->occ.n = str->occ.m = l; 
	for (k = 0; k < p->l1_ug->u.n; k++) {
		a = str->occ.a + str->idx.a[k];
		a_n = str->idx.a[k+1] - str->idx.a[k];
		if(a_n) a[a_n-1] = 0;
	}

	for (k = 0; k < str->str.n; k++) {
        z = &(str->str.a[k]);
		for (i = 0; i < z->cn; i++) {
			a = str->occ.a + str->idx.a[((uint32_t)z->a[i])>>1];
			a_n = str->idx.a[(((uint32_t)z->a[i])>>1)+1] - str->idx.a[((uint32_t)z->a[i])>>1];
			if(a_n) {
				if(a[a_n-1] == a_n-1) {
                    m = a_n-1; a[a_n-1] = (k<<32)|i; 
                } else {
                    m = a[a_n-1]; a[a[a_n-1]++] = (k<<32)|i;
                }

                v0 = (uint32_t)z->a[(uint32_t)a[m]];
                while (m > 0 && z->is_cir == 0) {
                    m--;
                    if((a[m]>>32) != k) break;
                    v1 = (uint32_t)z->a[(uint32_t)a[m]];
                    if(v0 == v1) z->is_cir = 1;
                }
			}
		}
	}
}

ul_resolve_t *init_ul_resolve_t(asg_t *sg, ma_ug_t *init_ug, bubble_type* bub, all_ul_t *idx, ug_opt_t *uopt, uint8_t *r_het)
{
    ul_resolve_t *p = NULL; CALLOC(p, 1);
    p->sg = sg; p->init_ug = init_ug; p->bub = bub; p->idx = idx; p->r_het = r_het; p->uopt = uopt;
    p->l1_ug = copy_untig_graph(p->init_ug);  
    init_integer_ml_t(&p->str_b, p, asm_opt.thread_num);
    // init_ul_path_srt_t(p);
    init_ul_str_idx_t(p);
    return p;
}

void print_bubble_gfa(FILE *fp, bubble_type *bub, const char* utg_pre, const char* bub_pre, const char* chain_pre)
{
    uint32_t i, k, m, *a, n, beg, sink, x; ma_utg_t *p; uint64_t occ;
    ma_ug_t *b_ug = bub->b_ug; char name[32], bname[32]; uint8_t *f; CALLOC(f, bub->ug->u.n);
    for (i = 0; i < b_ug->u.n; i++) {
        p = &b_ug->u.a[i];
        if(p->n == 0) continue;
        for (k = occ = 0; k < p->n; k++){
            x = p->a[k]>>33;
            get_bubbles(bub, x, &beg, &sink, &a, &n, NULL);

            for (m = 0; m < n; m++) {
                occ += bub->ug->u.a[a[m]>>1].n; f[a[m]>>1] = 1;
            }
            if(beg != (uint32_t)-1 && f[beg>>1] == 0) {
                occ += bub->ug->u.a[beg>>1].n; f[beg>>1] = 1;
            }
            if(sink != (uint32_t)-1 && f[sink>>1] == 0) {
                occ += bub->ug->u.a[sink>>1].n; f[sink>>1] = 1;
            }
        }

        sprintf(name, "%s%.6d%c", chain_pre, i + 1, "lc"[p->circ]);
        fprintf(fp, "S\t%s\t*\tLN:i:%lu\n", name, occ);
        for (k = 0; k < p->n; k++) {
            x = p->a[k]>>33;
            sprintf(bname, "%s%.6d", bub_pre, x + 1);
            fprintf(fp, "B\t%s\t%c\tcid:i:%s\tsm:%c\n", bname, "+-"[(p->a[k]>>32)&1], name, "01"[x<bub->f_bub]);

            get_bubbles(bub, x, &beg, &sink, &a, &n, NULL);
            if(beg != (uint32_t)-1) {
                fprintf(fp, "U\t%s%.6d%c\t%c\tcid:i:%s\tbid:b:%s\thom:%c\n", 
                utg_pre, (beg>>1)+1, "lc"[bub->ug->u.a[(beg>>1)].circ], "+-"[beg&1], name, bname, "10"[IF_HOM((beg>>1), *bub)]);
            }

            if(sink != (uint32_t)-1) {
                fprintf(fp, "U\t%s%.6d%c\t%c\tcid:i:%s\tbid:s:%s\thom:%c\n", 
                utg_pre, (sink>>1)+1, "lc"[bub->ug->u.a[(sink>>1)].circ], "+-"[sink&1], name, bname, "10"[IF_HOM((sink>>1), *bub)]);
            }
            for (m = 0; m < n; m++) {
                occ += bub->ug->u.a[a[m]>>1].n;
                fprintf(fp, "U\t%s%.6d%c\t%c\tcid:i:%s\tbid:c:%s\thom:%c\n", 
                utg_pre, (a[m]>>1)+1, "lc"[bub->ug->u.a[(a[m]>>1)].circ], "+-"[a[m]&1], name, bname, "10"[IF_HOM((a[m]>>1), *bub)]);
            }
        }
    }

    asg_arc_t* au = NULL;
    uint32_t nu, u, v, j;
    for (i = 0; i < b_ug->u.n; ++i) {
        if(b_ug->u.a[i].m == 0) continue;
        if(b_ug->u.a[i].circ)
        {
            fprintf(fp, "L\t%s%.6dc\t+\t%s%.6dc\t+\t%dM\tL1:i:%d\n", 
            chain_pre, i+1, chain_pre, i+1, 0, 0);
            fprintf(fp, "L\t%s%.6dc\t-\t%s%.6dc\t-\t%dM\tL1:i:%d\n", 
            chain_pre, i+1, chain_pre, i+1, 0, 0);
        } 
        u = i<<1;
        au = asg_arc_a(b_ug->g, u);
        nu = asg_arc_n(b_ug->g, u);
        for (j = 0; j < nu; j++)
        {
            if(au[j].del) continue;
            v = au[j].v;
            fprintf(fp, "L\t%s%.6d%c\t%c\t%s%.6d%c\t%c\t%dM\tL1:i:%d\n", 
            chain_pre, (u>>1)+1, "lc"[b_ug->u.a[u>>1].circ], "+-"[u&1],
            chain_pre, (v>>1)+1, "lc"[b_ug->u.a[v>>1].circ], "+-"[v&1], 0, 0);
        }


        u = (i<<1) + 1;
        au = asg_arc_a(b_ug->g, u);
        nu = asg_arc_n(b_ug->g, u);
        for (j = 0; j < nu; j++)
        {
            if(au[j].del) continue;
            v = au[j].v;
            fprintf(fp, "L\t%s%.6d%c\t%c\t%s%.6d%c\t%c\t%dM\tL1:i:%d\n", 
            chain_pre, (u>>1)+1, "lc"[b_ug->u.a[u>>1].circ], "+-"[u&1],
            chain_pre, (v>>1)+1, "lc"[b_ug->u.a[v>>1].circ], "+-"[v&1], 0, 0);
        }
    }
    
    for (i = 0; i < bub->ug->u.n; i++) {
        if(f[i]) continue;
        fprintf(fp, "U\t%s%.6d%c\t+\tcid:i:*\tbid:c:*\thom:%c\n", 
                utg_pre, i+1, "lc"[bub->ug->u.a[i].circ], "10"[IF_HOM(i, *bub)]);
    }
    free(f);
}

void print_uls_ovlp(FILE *fp, all_ul_t *uls, const char* utg_pre, ma_ug_t *ug)
{
    ul_vec_t *p = NULL; nid_t *z = NULL; uc_block_t *m = NULL; uint64_t k, i; uint32_t a, la, occ;
    kvec_t(uint8_t) f; kv_init(f);
    for (k = 0; k < uls->n; k++) {
        z = &(uls->nid.a[k]);
        p = &(uls->a[k]);
        for (i = 0; i < p->bb.n; i++) {
            m = &(p->bb.a[i]);
            if(m->base || (!m->el)) continue;
            if(m->pchain) break;
        }
        if(i >= p->bb.n) continue;
        kv_resize(uint8_t, f, p->bb.n); memset(f.a, 0, sizeof(*(f.a))*p->bb.n);
        fprintf(fp, ">\t%.*s\n", (int32_t)z->n, z->a);
        for (i = 0; i < p->bb.n; i++) {
            if(f.a[i]) continue;
            m = &(p->bb.a[i]);
            if(m->base || (!m->el) || (!m->pchain)) continue;
            for (a = i, la = i, occ = 0; a != (uint32_t)-1; a = p->bb.a[a].aidx) {
                la = a; f.a[a] = 1; occ++;
            }
            fprintf(fp, "C\tLEN:%u\tS:%u\tE:%u\tOCC:%u\n", p->rlen, p->bb.a[i].qs, p->bb.a[la].qe, occ);
            for (a = i; a != (uint32_t)-1; a = p->bb.a[a].aidx) {
                m = &(p->bb.a[a]);
                fprintf(fp, "%s%.6d%c[%c::", utg_pre, m->hid+1, "lc"[ug->u.a[(m->hid>>1)].circ], "+-"[m->rev]);
                if(m->pdis == (uint32_t)-1) fprintf(fp, "*]\t");
                else fprintf(fp, "%u]\t", m->pdis);
            }
            fprintf(fp, "\n");
        }
    }

    free(f.a);
}

void print_debug_ul(const char* o_n, ma_ug_t *ug, asg_t *rg, const ma_sub_t *cov, ma_hit_t_alloc* src, R_to_U* ridx, 
bubble_type *bub, all_ul_t *uls)
{
    char* gfa_name = (char*)malloc(strlen(o_n)+50); FILE *fn;

    if(ug && rg && cov && src && ridx) {
        uint32_t k, m, v, w, an; asg_arc_t *p, *av;
        for (k = 0; k < ug->g->n_arc; k++) {
            p = &(ug->g->arc[k]); v = p->v^1; w = (p->ul>>32)^1;
            av = asg_arc_a(ug->g, v); an = asg_arc_n(ug->g, v);
            for (m = 0; m < an; m++) {
                if(av[m].del == 0 && av[m].v == w) break;
            }
            assert(m < an && av[m].ou == p->ou);
        }
        sprintf(gfa_name, "%s.r_utg.noseq.gfa", o_n); fn = fopen(gfa_name, "w");
        ma_ug_print_simple(ug, rg, cov, src, ridx, "utg", fn);
        fclose(fn);
    }

    if(bub) {
        sprintf(gfa_name, "%s.bub.noseq.gfa", o_n); fn = fopen(gfa_name, "w");
        print_bubble_gfa(fn, bub, "utg", "btg", "ctg");
        fclose(fn);
    }

    if(uls && ug) {
        sprintf(gfa_name, "%s.uls.ovlp", o_n); fn = fopen(gfa_name, "w");
        print_uls_ovlp(fn, uls, "utg", ug);
        fclose(fn);
    }

    free(gfa_name);
    fprintf(stderr, "[M::%s::] done\n", __func__);
    exit(1);
}

int64_t normlize_gdis(ma_ug_t *ug, uc_block_t *i, uc_block_t *k, int64_t is_i2k_forward)
{
    int64_t i_len, k_len;
    ///i > k
    // if(is_i2k_forward == 0){
    //     if(!i->rev) {
    //         i_len = i->qe + (ug->g->seq[i->hid].len - i->te);
    //     } else {
    //         i_len = i->qe + i->ts;
    //     }

    //     if(!k->rev) {
    //         k_len = k->qe + (ug->g->seq[k->hid].len - k->te);
    //     } else {
    //         k_len = k->qe + k->ts;
    //     }
    // } else {
    //     if(!i->rev) {
    //         i_len = (int64_t)i->qs - (int64_t)i->ts;
    //     } else {
    //         i_len = (int64_t)i->qs - (int64_t)(ug->g->seq[i->hid].len - i->te);
    //     }
    //     if(i_len < 0) i_len = 0;

    //     if(!k->rev) {
    //         k_len = (int64_t)k->qs - (int64_t)k->ts;
    //     } else {
    //         k_len = (int64_t)k->qs - (int64_t)(ug->g->seq[k->hid].len - k->te);
    //     }
    //     if(k_len < 0) k_len = 0;
    // }

    if(!i->rev) {
        i_len = i->qe + (ug->g->seq[i->hid].len - i->te);
    } else {
        i_len = i->qe + i->ts;
    }
    

    if(!k->rev) {
        k_len = k->qe + (ug->g->seq[k->hid].len - k->te);
    } else {
        k_len = k->qe + k->ts;
    }

    if(is_i2k_forward) {
        i_len -= ug->g->seq[i->hid].len;
        k_len -= ug->g->seq[k->hid].len;
    }
    

    if(i_len >= k_len) return i_len - k_len;
    // assert(i_len >= k_len);
    return 0;
}

int64_t normlize_gdis_exact(ma_ug_t *ug, uc_block_t *a, uint32_t i, uint32_t k, int64_t is_i2k_forward)
{
    assert(i > k);
    uint32_t li, lk, pk, bi = i; int64_t l;
    for (li = i, l = 0; i != (uint32_t)-1 && i >= k; i = a[i].pidx) {
        li = i; 
        if(a[i].pidx != (uint32_t)-1 && i > k) l += a[i].pdis;
    }
    
    if(is_i2k_forward) {
        l += (int64_t)(ug->g->seq[a[li].hid].len);
        l -= (int64_t)(ug->g->seq[a[bi].hid].len);
    }
    i = li;
    if(i == k) return l;
    assert(i > k); pk = k;
    for (lk = k; k != (uint32_t)-1 && k <= i; k = a[k].aidx) lk = k;
    for (k = lk; k != (uint32_t)-1 && k != pk; k = a[k].pidx) l += a[k].pdis;
    // assert(k == pk);
    if(is_i2k_forward) {
        l += (int64_t)(ug->g->seq[a[pk].hid].len);
        l -= (int64_t)(ug->g->seq[a[lk].hid].len);
    }
    if(l < 0) l = 0;
    k = lk;
    assert(i > k);
    return l + normlize_gdis(ug, &(a[i]), &(a[k]), is_i2k_forward);
}

uint32_t dis_check_integer_aln_t(all_ul_t *ul_idx, ul_str_idx_t *str_idx, ma_ug_t *ug, integer_aln_t *li, integer_aln_t *lk, 
uint32_t qid, uint32_t tid, uint32_t is_rev, double diff_rate, int64_t hard_thres)
{
    assert(((uint32_t)li->tn_rev_qk) >= ((uint32_t)lk->tn_rev_qk));
    if(((uint32_t)li->tn_rev_qk) == ((uint32_t)lk->tn_rev_qk)) return 0;
    if(li->tk <= lk->tk) return 0;
    ul_str_t *q = &(str_idx->str.a[qid]), *t = &(str_idx->str.a[tid]);
    uc_block_t *iq, *it, *kq, *kt; int64_t qlen, tlen, mm, dd, i_qk, k_qk, i_tk, k_tk;
    i_qk = ((uint32_t)li->tn_rev_qk); k_qk = ((uint32_t)lk->tn_rev_qk);
    if(!is_rev) {
        i_tk = li->tk; k_tk = lk->tk;
    } else {
        i_tk = t->cn - lk->tk - 1; k_tk = t->cn - li->tk - 1;
    }
    iq = &(ul_idx->a[qid].bb.a[(q->a[i_qk]>>32)]);
    it = &(ul_idx->a[tid].bb.a[(t->a[i_tk]>>32)]);

    kq = &(ul_idx->a[qid].bb.a[(q->a[k_qk]>>32)]); 
    kt = &(ul_idx->a[tid].bb.a[(t->a[k_tk]>>32)]);
    qlen = normlize_gdis(ug, iq, kq, 0); tlen = normlize_gdis(ug, it, kt, is_rev);
    if(qlen < tlen) {
        mm = qlen; dd = tlen - qlen;
    } else {
        mm = tlen; dd = qlen - tlen;
    }

    // if(tid == 3074) {
    //     fprintf(stderr, "\n+++[M::%s::] qlen::%ld, tlen::%ld, iq->hid::%u, kq->hid::%u, (q->a[i_qk]>>32)::%lu, (q->a[k_qk]>>32)::%lu, (t->a[i_tk]>>32)::%lu, (t->a[k_tk]>>32)::%lu\n", 
    //     __func__, qlen, tlen, iq->hid, kq->hid, 
    //     (q->a[i_qk]>>32), (q->a[k_qk]>>32), (t->a[i_tk]>>32), (t->a[k_tk]>>32));
    // }

    if(dd <= (mm*diff_rate) || dd < hard_thres) return 1;
    
    qlen = normlize_gdis_exact(ug, ul_idx->a[qid].bb.a, (q->a[i_qk]>>32), (q->a[k_qk]>>32), 0);
    tlen = normlize_gdis_exact(ug, ul_idx->a[tid].bb.a, (t->a[i_tk]>>32), (t->a[k_tk]>>32), is_rev);
    if(qlen < tlen) {
        mm = qlen; dd = tlen - qlen;
    } else {
        mm = tlen; dd = qlen - tlen;
    }

    if(dd <= (mm*diff_rate) || dd < hard_thres) return 1;

    // if(tid == 3074) {
    //     fprintf(stderr, "---[M::%s::] qlen::%ld, tlen::%ld, iq->hid::%u, kq->hid::%u\n", 
    //     __func__, qlen, tlen, iq->hid, kq->hid);
    // }
    

    return 0;
}

uint32_t integer_chain(uint32_t qid, integer_aln_t *a, int64_t a_n, int64_t offset, integer_t *buf, 
ma_ug_t *ug, ul_str_idx_t *str_idx, all_ul_t *ul_idx, ul_chain_t *res)
{
    res->v = res->s = res->e = (uint32_t)-1; res->sc = (uint64_t)-1; 
    res->q_sidx = res->q_eidx = res->t_sidx = res->t_eidx = (uint32_t)-1;
    if(a_n <= 0) return 0;
    ///already sorted by qe
    int64_t i, k, max_f, max_k, sc, csc, *p, *f, tf, ti/**, is_circle**/;
    integer_aln_t *li, *lk; uint32_t tid = a[0].tn_rev_qk>>33; uint32_t is_rev = (a[0].tn_rev_qk>>32)&1;
    for (i = 1, sc = 0; i < a_n; ++i) {
        sc += a[i].sc;
        if(a[i].tk <= a[i-1].tk) break;///== means there is a circle
        if(((uint32_t)a[i].tn_rev_qk) <= ((uint32_t)a[i-1].tn_rev_qk)) break;
        if(str_idx && ul_idx && (!dis_check_integer_aln_t(ul_idx, str_idx, 
                                    ug, &(a[i]), &(a[i-1]), qid, tid, is_rev, 0.08, 2000))) {
            break;
        }
    }
    // if(tid == 269 || tid == 276 || tid == 277 || tid == 278) fprintf(stderr, "tid->%u, i->%ld, an->%ld\n", tid, i, a_n);
    if(i >= a_n) {//the whole chain is co-inear
        sc += a[0].sc; 
        res->v = a[0].tn_rev_qk>>32; res->s = offset; res->e = offset + a_n; res->sc = sc; 
        return 1;
    }
    // is_circle = 1;
    // if(str_idx->str.a[qid].is_cir == 0 && str_idx->str.a[tid].is_cir == 0) is_circle = 0;

    buf->p.n = buf->f.n = 0;
    kv_resize(int64_t, buf->p, (uint64_t)a_n); p = buf->p.a;
    kv_resize(int64_t, buf->f, (uint64_t)a_n); f = buf->f.a;

    tf = ti = -1;
    for (i = 0; i < a_n; ++i) {
        li = &(a[i]); csc = a[i].sc;
        max_f = csc; max_k = -1;
        for (k = i-1; k >= 0; --k) { 
            lk = &(a[k]);
            ///qk of lk and li might be equal
            if(lk->tk >= li->tk || ((uint32_t)lk->tn_rev_qk) >= ((uint32_t)li->tn_rev_qk)) continue;
            // if(is_circle && (!dis_check_integer_aln_t(ul_idx, str_idx, ug, li, lk, qid, tid, is_rev, 0.08))) continue;
            if(str_idx && ul_idx && (!dis_check_integer_aln_t(ul_idx, str_idx, 
                                                    ug, li, lk, qid, tid, is_rev, 0.08, 2000))) {
                continue;
            }
            sc = csc + f[k];
            if(sc > max_f) {
                max_f = sc; max_k = k;
            }
        }
        f[i] = max_f; p[i] = max_k;
        if(tf < max_f) {
            tf = max_f; ti = i;
        }
    }

    if(ti < 0) return 0;
    
    for (i = ti, k = 0; i >= 0; i = p[i]) f[k++] = i;
    assert(k > 0);
    for (i = sc = 0, k--; k >= 0; k--, i++) {
        a[i] = a[f[k]]; sc += a[i].sc;
        // if(tid == 269) {
        //     fprintf(stderr, "[%ld] qk->%u, tk->%u\n", i, (uint32_t)a[i].tn_rev_qk, a[i].tk);
        // }
    }
    res->v = a[0].tn_rev_qk>>32; res->s = offset; res->e = offset + i; res->sc = sc;
    return 1;
}

/**
void calculate_boundary_integer_length(all_ul_t *ul_idx, ul_str_t *str, int64_t qid, int64_t tid, int64_t qk, int64_t tk, 
int64_t is_rev, int64_t is_prefix, int64_t is_suffix, int64_t *r_qoff, int64_t *r_toff)
{
    (*r_qoff) = (*r_toff) = -1;
    int64_t qlen = ul_idx->a[qid].rlen, tlen = ul_idx->a[tid].rlen, qoff, toff, q_ext, t_ext;
    ul_str_t *qstr = &(str[qid]), *tstr = &(str[tid]); 
    if(is_rev) tk = tstr->cn - tk;
    uc_block_t *q_b = &(ul_idx->a[qid].bb.a[qstr->a[qk]>>32]);
    uc_block_t *t_b = &(ul_idx->a[tid].bb.a[tstr->a[tk]>>32]);
    if(is_prefix) {
        qoff = q_b->qs; toff = (is_rev?(tlen-t_b->qe):(t_b->qs));
    }
    if(is_suffix) {
        qoff = qlen - q_b->qe; toff = (is_rev?(t_b->qs):(tlen-t_b->qe));
    }

    if(qoff <= toff) {
        t_ext = qoff; q_ext = -1;
    } else {
        q_ext = toff; t_ext = -1;
    }

    if(q_ext == -1) {
        if(is_prefix) (*r_qoff) = 0;
        if(is_suffix) (*r_qoff) = (int64_t)(qstr->cn);
    }

}

void push_ul_snp_t(int64_t chain_id, ul_str_t *str, integer_t *buf, integer_aln_t *a, int64_t a_n, int64_t qid, int64_t tid, int64_t is_rev)
{
    if(a_n <= 0) return;
    int64_t k, qk, tk, p_qk, p_tk; ul_snp_t *p;
    ul_str_t *qstr = &(str[qid]), *tstr = &(str[tid]);
    for (k = 0, p_qk = 0, p_tk = 0; k < a_n; k++) {
        qk = (uint32_t)a[k].tn_rev_qk; tk = a[k].tk; 
        if(qk - p_qk > 0 || tk - p_tk > 0) {
            if(p_qk == 0 && p_tk == 0) {///first window
                ///qk == 0 || tk == 0 means we already reach the end
                if(qk > 0 && tk > 0) {
                    ;
                }
            } else {
                kv_pushp(ul_snp_t, buf->snp, &p);
                p->chain_id = chain_id; p->is_rev = is_rev;
                p->qidx_occ = p_qk; p->qidx_occ += (qk - p_qk);            
                p->tidx_occ = p_tk; p->tidx_occ += (tk - p_tk);
            }
        }
        p_qk = qk + 1; p_tk = tk + 1;
    }
    
    qk = qstr->cn; tk = tstr->cn;
    if(qk - p_qk > 0 || tk - p_tk > 0) {
        kv_pushp(ul_snp_t, buf->snp, &p);
        p->chain_id = chain_id; p->is_rev = is_rev;
        p->qidx_occ = p_qk; p->qidx_occ += (qk - p_qk);            
        p->tidx_occ = p_tk; p->tidx_occ += (tk - p_tk);
    }
}

void integer_variant_call(integer_t *buf, ul_snp_t *a, int64_t a_n)
{
    kv_resize(int64_t, buf->f, (uint64_t)a_n);
    int64_t *f = buf->f.a; int64_t k, m; ul_snp_t *p;
    memset(f, -1, sizeof((*f))*a_n);
    for (k = 0; k < a_n; k++) {
        if(f[k] != (uint64_t)-1) continue;
        for (m = k+1, p = &(a[k]); m < a_n; m++) {
            
        }
    }
    
    
}

void integer_phase(ul_str_t *str, integer_t *buf, ul_chain_t *idx, int64_t idx_n, integer_aln_t *aln, int64_t qid)
{
    int64_t k; integer_aln_t *a; buf->snp.n = 0;
    for (k = 0; k < idx_n; k++) {
        push_ul_snp_t(k, str, buf, aln + idx[k].s, idx[k].e - idx[k].s, qid, idx[k].v>>1, idx[k].v&1);
    }   

    int64_t z, snp_n = buf->snp.n;
    radix_sort_ul_snp_t_srt(buf->snp.a, buf->snp.a + buf->snp.n);
    for (z = 0, k = 1; k <= snp_n; k++) {
        if(k == snp_n || buf->snp.a[z].qidx_occ != buf->snp.a[k].qidx_occ) {
            integer_variant_call(buf, buf->snp.a + z, k - z);
            z = k;
        }
    }
}
**/

int64_t append_connective(integer_aln_t *aln, ul_chain_t *idx, int64_t str_i, uint64_t occ_thres, uint64_t *res)
{
    int64_t k, kl = idx->e - idx->s, str_k = -1, fp = 0; integer_aln_t *a = aln + idx->s; 
    assert(idx->sc <= (uint64_t)kl);
    for (k = idx->sc; k < kl; k++) {
        str_k = (uint32_t)a[k].tn_rev_qk;
        if(str_k >= str_i) break;
    }
    idx->sc = k;
    if(k >= kl || str_k != str_i) return 0;

    for (k -= 1; k >= 0; k--) {
        str_k = (uint32_t)a[k].tn_rev_qk;
        if(res[str_k] < occ_thres) {
            res[str_k]++;
            if(res[str_k] == occ_thres) fp++;
        }
    }
    return fp;
}

int64_t connective_conform(integer_aln_t *aln, ul_chain_t *idx_a, int64_t idx_n, int64_t str_i0, int64_t occ_thres, uint64_t is_cov_check)
{
    if(str_i0 == 0) return 1;
    int64_t z; ul_chain_t *x; int64_t k, kl, str_k, match, exact; integer_aln_t *a;
    for (z = match = exact = 0; z < idx_n; z++) {
        x = &(idx_a[z]);
        kl = x->e - x->s; str_k = -1; a = aln + x->s; 
        assert(x->sc <= (uint64_t)kl);
        for (k = x->sc; k < kl; k++) {
            str_k = (uint32_t)a[k].tn_rev_qk;
            if(str_k >= str_i0) break;
        }
        x->sc = k;
        if(k >= kl || str_k != str_i0) continue;
        k--;
        if(k >= 0) {
            str_k = (uint32_t)a[k].tn_rev_qk;
            if(str_k+1 == str_i0) {
                match++;
                if(a[k].tk+1 == a[k+1].tk) exact++;
            }
        }
    }

    if(is_cov_check) {
        if(match < occ_thres) return 0;
    } else {
        assert(match >= occ_thres);
    }

    if(exact == match) return 1;
    if(exact > (match*0.51) && exact > (match/2)) return 1;
    return 0;
}

///occ_thres does not consider reference read itself; so the real coverage is (occ_thres+1)
int64_t integer_chain_dp(bubble_type *bub, integer_t *buf, ul_str_t *str, integer_aln_t *aln, ul_chain_t *idx, int64_t idx_n, int64_t qid, uint32_t is_hom, 
uint64_t occ_thres, uint64_t *corrected)
{
    ul_str_t *qstr = &(str[qid]); int64_t k, q_n = qstr->cn, *f, *p, z, max_f, tf, tk, max_p, done_z, sc, csc, n_skip;
    kv_resize(int64_t, buf->f, qstr->cn); kv_resize(int64_t, buf->p, qstr->cn); 
    kv_resize(uint64_t, buf->o, qstr->cn); uint64_t *o;
    f = buf->f.a; p = buf->p.a; o = buf->o.a; if(corrected) (*corrected) = 0;

    // radix_sort_ul_chain_t_srt(idx, idx + idx_n);
    for (k = 0; k < idx_n; k++) idx[k].sc = 0;

    for (k = 0, tf = tk = -1, n_skip = 0; k < q_n; k++) {
        csc = buf->u.a[k]; 
        if((!is_hom) && (IF_HOM((((uint32_t)qstr->a[k])>>1), (*bub)))) {
            csc = -1; n_skip++;
        }
        max_p = -1; max_f = csc;
        if(k > 0 && max_f >= 0) {
            done_z = 0;
            if(is_hom) {//check all nodes
                memset(o, 0, sizeof((*o))*k); 
            } else {///mask hom nodes
                for (z = 0; z < k; z++) {
                    o[z] = 0;
                    if(f[z] == -1) {
                        o[z] = occ_thres; ///if node is hom, ignore it
                        done_z++;
                    }
                }
            }
            
            for (z = 0; z < idx_n && done_z < k; z++) {
                done_z += append_connective(aln, &(idx[z]), k, occ_thres, o);
            }
            assert(done_z <= k);
            for (z = k - 1; z >= 0; z--) {
                if(o[z] < occ_thres) continue;
                if(f[z] == -1) continue;///masked hom nodes
                sc = csc + f[z];
                if(sc > max_f) {
                    max_f = sc; max_p = z;
                }
            }
        }
        f[k] = max_f; p[k] = max_p;
        if(tf < max_f && max_f >= 0) {
            tf = max_f; tk = k;
        }

        // fprintf(stderr, "[M::%s::k->%ld] f[k]->%ld, p[k]->%ld, csc->%ld\n", __func__, k, f[k], p[k], csc);
    }
    if(tk < 0) return 0;
    for (k = tk, done_z = 0; k >= 0; k = p[k]) done_z++;
    ///might be fully corrected; but some of hom nodes have been skipped
    if((corrected) && ((n_skip+done_z) == q_n) ) {
        for (k = 0; k < idx_n; k++) idx[k].sc = 0;
        for (k = 1, (*corrected) = 0; k < q_n; k++) {
            if(!connective_conform(aln, idx, idx_n, k, occ_thres, (f[k] >= 0 && f[k-1] >= 0)?0:1)) break;
        }
        if(k >= q_n) (*corrected) = 1;
    }

    for (k = tk, sc = done_z; k >= 0; k = p[k]) o[--done_z] = k;
    return sc;
}

void print_integer_seq(ma_ug_t *ug, ul_str_t *str, int64_t id, int64_t is_header)
{
    uint64_t i;
    if(is_header) {
        fprintf(stderr,"[M::%s::tid->%ld] occ::%u\n", __func__, id, str[id].cn);
    }
    for (i = 0; i < str[id].cn; i++) {
        fprintf(stderr, "utg%.6d%c(%c)\t", (((uint32_t)str[id].a[i])>>1)+1, 
        "lc"[ug->u.a[(((uint32_t)str[id].a[i])>>1)].circ], "+-"[(((uint32_t)str[id].a[i])&1)]);
    }
    fprintf(stderr,"\n");
}


#define poa_arc_n(g, v) ((uint32_t)(g)->idx.a[(v)])
#define poa_arc_a(g, v) (&(g)->arc.a[(g)->idx.a[(v)]>>32])

void print_integer_g(poa_g_t *g, ma_ug_t *ug, uint64_t is_gfa)
{
    uint64_t i, v, w;
    if(!is_gfa) {
        for (i = 0; i < g->seq.n; i++) {
            fprintf(stderr, "Node::[M::%s::i->%lu] utg%.6d%c(%c)\n", __func__, i, (int32_t)(g->seq.a[i].nid>>1)+1, 
            "lc"[ug->u.a[g->seq.a[i].nid>>1].circ], "+-"[g->seq.a[i].nid&1]);
        }
        for (i = 0; i < g->arc.n; i++) {
            v = g->arc.a[i].ul>>32; w = g->arc.a[i].v;
            if(v&1) continue;
            v = g->seq.a[v>>1].nid; w = g->seq.a[w>>1].nid;
            fprintf(stderr, "Arch::[M::%s::w->%u] utg%.6d%c(%c)<nid::%lu> -> utg%.6d%c(%c)<nid::%u>\n", __func__, (uint32_t)g->arc.a[i].ul, 
            (int32_t)(v>>1)+1, "lc"[ug->u.a[v>>1].circ], "+-"[v&1], g->arc.a[i].ul>>33, 
            (int32_t)(w>>1)+1, "lc"[ug->u.a[w>>1].circ], "+-"[w&1], g->arc.a[i].v>>1);
        }
        for (i = 0; i < g->seq.n; i++) {
            fprintf(stderr, "Srt::[M::%s::i->%lu] utg%.6d%c(%c)<nid::%u>\n", __func__, i, 
            (int32_t)(g->seq.a[g->srt_b.res.a[i]].nid>>1)+1, 
            "lc"[ug->u.a[g->seq.a[g->srt_b.res.a[i]].nid>>1].circ], "+-"[g->seq.a[g->srt_b.res.a[i]].nid&1], g->srt_b.res.a[i]);
        }
    } else {
        char name[32];
        for (i = 0; i < g->seq.n; i++) {
            sprintf(name, "dtg%.6d", (int)i);
            fprintf(stderr, "S\t%s\t*\tLN:i:%u\trd:i:utg%.6d%c(%c)\n", name, g->srt_b.res2nid.a[i],
            (int32_t)(g->seq.a[i].nid>>1)+1, "lc"[ug->u.a[g->seq.a[i].nid>>1].circ], "+-"[g->seq.a[i].nid&1]);
        }
        uint32_t u, nu, j; poa_arc_t *au;
        for (i = 0; i < g->seq.n; ++i) {
            u = i<<1;
            au = poa_arc_a(g, u); nu = poa_arc_n(g, u);
            for (j = 0; j < nu; j++) {
                v = au[j].v;
                fprintf(stderr, "L\tdtg%.6d\t%c\tdtg%.6d\t%c\t%dM\n", 
                (int)(u>>1), "+-"[u&1], (int)(v>>1), "+-"[v&1], asg_arc_len(au[j]));
            }

            u = (i<<1) + 1;
            au = poa_arc_a(g, u); nu = poa_arc_n(g, u);
            for (j = 0; j < nu; j++) {
                v = au[j].v;
                fprintf(stderr, "L\tdtg%.6d\t%c\tdtg%.6d\t%c\t%dM\n", 
                (int)(u>>1), "+-"[u&1], (int)(v>>1), "+-"[v&1], asg_arc_len(au[j]));
            }
        }
    }
}

void print_cns_seq(ma_ug_t *ug, ul_str_t *str, uint64_t *cns_seq, uint64_t cns_occ)
{
    fprintf(stderr,"[M::%s::]\t", __func__);
    uint64_t k, ck;
    for (k = ck = 0; k < str->cn; k++) {
        for (; ck < cns_occ; ck++) {
            if(cns_seq[ck] >= k) break;
        }
        if(ck < cns_occ && cns_seq[ck] == k) {
            fprintf(stderr, "utg%.6d%c(%c)\t", (((uint32_t)str->a[k])>>1)+1, 
                    "lc"[ug->u.a[(((uint32_t)str->a[k])>>1)].circ], "+-"[(((uint32_t)str->a[k])&1)]);
        } else {
            fprintf(stderr, "*\t");
        }
    }
    fprintf(stderr,"\n");
}

void print_res_seq(poa_g_t *pg, ma_ug_t *ug, uint32_t *cns_seq, uint32_t cns_occ)
{
    fprintf(stderr,"[M::%s::]\t", __func__);
    uint64_t k;
    for (k = 0; k < cns_occ; k++) {
        fprintf(stderr, "utg%.6d%c(%c)\t", (pg->seq.a[(cns_seq[k]>>1)].nid>>1)+1, 
                "lc"[ug->u.a[(pg->seq.a[(cns_seq[k]>>1)].nid>>1)].circ], 
                "+-"[(pg->seq.a[(cns_seq[k]>>1)].nid&1)]);
        
    }
    fprintf(stderr,"\n");
}

int64_t utg_cover_read_occ_by_qs(ma_ug_t *ug, int64_t oqs, int64_t oqe, uc_block_t *x)
{
    assert(oqs >= (int64_t)x->qs && oqs <= (int64_t)x->qe);
    assert(oqe >= (int64_t)x->qs && oqe <= (int64_t)x->qe);
    assert(oqe > oqs);
    int64_t s_off, e_off, k;
    s_off = get_offset_adjust(oqs-x->qs, x->qe-x->qs, x->te-x->ts);
	e_off = get_offset_adjust(x->qe-oqe, x->qe-x->qs, x->te-x->ts);
    if(x->rev) {
        k = s_off; s_off = e_off; e_off = k;
    }

    return ug_occ_w(x->ts + s_off, x->te - e_off, &(ug->u.a[x->hid]));
} 

int64_t estimate_ul_len(ma_ug_t *ug, ul_vec_t *raw_ov, ul_str_t *str, int64_t idx, int64_t ext_len, uint64_t *cov_buf)
{
    if(ext_len == 0) return idx;
    int64_t aim_s, aim_e, k = idx, qs, qe, ovlp_s, ovlp_e, cov_occ, str_n = str->cn;
    if(ext_len < 0) {
        aim_s = ((int64_t)(raw_ov->bb.a[str->a[idx]>>32].qs)) + ext_len; 
        if(aim_s < 0) aim_s = 0; 
        aim_e = raw_ov->bb.a[str->a[idx]>>32].qe;
        for (k = idx - 1; k >= 0; k--) {
            qs = raw_ov->bb.a[str->a[k]>>32].qs; qe = raw_ov->bb.a[str->a[k]>>32].qe;
            // if(qe <= aim_s) break;
            ovlp_s = MAX(aim_s, qs); ovlp_e = MIN(aim_e, qe);
            if(ovlp_e <= ovlp_s) {
                k++;
                break;
            }
            cov_occ = utg_cover_read_occ_by_qs(ug, ovlp_s, ovlp_e, &(raw_ov->bb.a[str->a[k]>>32]));
            if(cov_occ == 0) {
                k++;
                break;
            }
            cov_buf[k] = cov_occ;
        }
        if(k < 0) k++;
    } else {
        aim_s = raw_ov->bb.a[str->a[idx]>>32].qs;
        aim_e = ((int64_t)(raw_ov->bb.a[str->a[idx]>>32].qe)) + ext_len; 
        if(aim_e > (int64_t)(raw_ov->rlen)) aim_e = raw_ov->rlen;
        for (k = idx + 1; k < str_n; k++) {
            qs = raw_ov->bb.a[str->a[k]>>32].qs; qe = raw_ov->bb.a[str->a[k]>>32].qe;
            // if(qs >= aim_e)
            ovlp_s = MAX(aim_s, qs); ovlp_e = MIN(aim_e, qe);
            if(ovlp_e <= ovlp_s) {
                k--;
                break;
            }
            cov_occ = utg_cover_read_occ_by_qs(ug, ovlp_s, ovlp_e, &(raw_ov->bb.a[str->a[k]>>32]));
            if(cov_occ == 0) {
                k--;
                break;
            }
            cov_buf[k] = cov_occ;
        }
        if(k >= str_n) k--;
    }
    return k;
}

void print_integer_ovlps(ma_ug_t *ug, ul_str_t *str, integer_aln_t *aln, int64_t aln_occ, ul_chain_t *idx, int64_t idx_n, int64_t qid, int64_t consenus_occ)
{
    fprintf(stderr, "\n[M::%s::qid->%ld] qstr->cn::%u, aln_occ::%ld, idx_n::%ld, consenus_occ::%ld\n", 
                                __func__, qid, str[qid].cn, aln_occ, idx_n, consenus_occ);
    // print_integer_seq(ug, str, qid, 0);
    int64_t k, z, z_n, qk, tk, is_rev, tid; //uint64_t z;
    for (k = 0; k < idx_n; k++) {
        fprintf(stderr, "\n[M::%s::tid->%lu] rev->%lu, aln_n->%u\n", __func__, aln[idx[k].s].tn_rev_qk>>33, (aln[idx[k].s].tn_rev_qk>>32)&1, idx[k].e - idx[k].s);
        print_integer_seq(ug, str, aln[idx[k].s].tn_rev_qk>>33, 0);
        z = idx[k].s; z_n = idx[k].e; tid = aln[idx[k].s].tn_rev_qk>>33;
        for (; z < z_n; z++) {
            qk = (uint32_t)aln[z].tn_rev_qk; is_rev = ((aln[z].tn_rev_qk>>32)&1);
            tk = ((is_rev == 0)? (aln[z].tk):(str[tid].cn - aln[z].tk - 1));
            fprintf(stderr, "[qk::%ld]utg%.6d%c(%c) <---> [tk::%ld]utg%.6d%c(%c)\n", 
            qk, (((uint32_t)str[qid].a[qk])>>1)+1, "lc"[ug->u.a[(((uint32_t)str[qid].a[qk])>>1)].circ], "+-"[(((uint32_t)str[qid].a[qk])&1)], 
            tk, (((uint32_t)str[tid].a[tk])>>1)+1, "lc"[ug->u.a[(((uint32_t)str[tid].a[tk])>>1)].circ], "+-"[(((uint32_t)str[tid].a[tk])&1)]);
        }
        
        // for (z = idx[k].s; z < idx[k].e; z++) {
        //     aln[z].tn_rev_qk
        // }
    }
}

void integer_align_extention(ma_ug_t *ug, all_ul_t *ul_idx, int64_t qid, ul_str_t *q_str, int64_t tid, ul_str_t *t_str, int64_t is_rev,
uint64_t *q_cov_buf, uint64_t *t_cov_buf, integer_aln_t *aln_pair, int64_t is_backward, int64_t *r_q_end, int64_t *r_t_end)
{
    int64_t qk, tk, qoff, toff, qlen, tlen, ext, q_end, t_end; uc_block_t *q_b, *t_b;
    qk = (uint32_t)(aln_pair->tn_rev_qk); tk = ((is_rev == 0)? (aln_pair->tk):(t_str->cn - aln_pair->tk - 1));
    q_b = &(ul_idx->a[qid].bb.a[q_str->a[qk]>>32]); t_b = &(ul_idx->a[tid].bb.a[t_str->a[tk]>>32]);
    qlen = ul_idx->a[qid].rlen; tlen = ul_idx->a[tid].rlen;

    if(is_backward) {
        qoff = q_b->qs; toff = (is_rev?(tlen-t_b->qe):(t_b->qs));
    } else {
        qoff = qlen - q_b->qe; toff = (is_rev?(t_b->qs):(tlen-t_b->qe));
    }

    if (qoff <= toff) ext = qoff;
    else ext = toff;
    
    if(is_backward) {
        q_end = estimate_ul_len(ug, &(ul_idx->a[qid]), q_str, qk, -ext, q_cov_buf);
        t_end = estimate_ul_len(ug, &(ul_idx->a[tid]), t_str, tk, ((is_rev)?(ext):(-ext)), t_cov_buf);
    } else {
        q_end = estimate_ul_len(ug, &(ul_idx->a[qid]), q_str, qk, ext, q_cov_buf);
        t_end = estimate_ul_len(ug, &(ul_idx->a[tid]), t_str, tk, ((is_rev)?(-ext):(ext)), t_cov_buf);
    }

    (*r_q_end) = q_end; (*r_t_end) = t_end;
}

int64_t gap_chain_check(ma_ug_t *ug, integer_aln_t *i0, integer_aln_t *i1, int64_t is_rev, ul_str_t *q_str, ul_str_t *t_str, 
uc_block_t *q_block, uc_block_t *t_block, int64_t qid, int64_t tid, double diff_rate, int64_t hard_thres)
{
    int64_t i0_qk, i0_tk, i1_qk, i1_tk, q_near, t_near, lq, lt, min, max, dif; uint32_t e, k;
    i0_qk = (uint32_t)(i0->tn_rev_qk); i1_qk = (uint32_t)(i1->tn_rev_qk); 
    if(is_rev == 0) {
        i0_tk = i0->tk; i1_tk = i1->tk;
    } else {
        i1_tk = t_str->cn - i0->tk - 1; i0_tk = t_str->cn - i1->tk - 1;
    }
    assert(i0_qk < i1_qk); assert(i0_tk < i1_tk);
    q_near = t_near = 0;
    if((i0_qk + 1) == i1_qk) q_near = 1;
    if((i0_tk + 1) == i1_tk) t_near = 1;
    if((q_near + t_near) != 1) return 0;

    i0_qk = q_str->a[i0_qk]>>32; i1_qk = q_str->a[i1_qk]>>32; 
    i0_tk = t_str->a[i0_tk]>>32; i1_tk = t_str->a[i1_tk]>>32; 

    k = i1_qk; e = i0_qk; lq = 0;
    while ((k != (uint32_t)-1) && (k != e)) {
        if(q_block[k].pidx != (uint32_t)-1) lq += q_block[k].pdis;
        k = q_block[k].pidx;
    }
    if(k != e) return 0;
    
    k = i1_tk; e = i0_tk; lt = 0;
    while ((k != (uint32_t)-1) && (k != e)) {
        if(t_block[k].pidx != (uint32_t)-1) lt += t_block[k].pdis;
        k = t_block[k].pidx;
    }
    if(k != e) return 0;
    if(is_rev) {
        lt += (int64_t)ug->g->seq[t_block[i0_tk].hid].len;
        lt -= (int64_t)ug->g->seq[t_block[i1_tk].hid].len;
    }
    if(lt < 0) lt = 0;
    if(lq <= lt) {
        min = lq; max = lt;
    } else {
        min = lt; max = lq;
    }

    dif = max - min;
    if(dif >= (min*diff_rate)) {
        // fprintf(stderr, "+[M::%s] qid::%ld, tid::%ld, i0_qk::%ld, i1_qk::%ld, lq::%ld, lt::%ld\n", 
        // __func__, qid, tid, i0_qk, i1_qk, lq, lt);
        lq = normlize_gdis(ug, &(q_block[i1_qk]), &(q_block[i0_qk]), 0);
        lt = normlize_gdis(ug, &(t_block[i1_tk]), &(t_block[i0_tk]), is_rev);
        // fprintf(stderr, "-[M::%s] qid::%ld, tid::%ld, i0_qk::%ld, i1_qk::%ld, lq::%ld, lt::%ld\n", 
        // __func__, qid, tid, i0_qk, i1_qk, lq, lt);
        // fprintf(stderr, "[M::%s] q_block[i0_qk].qs::%u, q_block[i0_qk].qe::%u, q_block[i1_qk].qs::%u, q_block[i1_qk].qe::%u\n", 
        // __func__, q_block[i0_qk].qs, q_block[i0_qk].qe, q_block[i1_qk].qs, q_block[i1_qk].qe);

        // fprintf(stderr, "[M::%s] t_block[i0_tk].qs::%u, t_block[i0_tk].qe::%u, t_block[i1_tk].qs::%u, t_block[i1_tk].qe::%u\n", 
        // __func__, t_block[i0_tk].qs, t_block[i0_tk].qe, t_block[i1_tk].qs, t_block[i1_tk].qe);
        if(lq <= lt) {
            min = lq; max = lt;
        } else {
            min = lt; max = lq;
        }

        dif = max - min;
        if(dif >= (min*diff_rate) || dif < hard_thres) return 0;
    }
    return 1;
}

int64_t refine_integer_ovlps(all_ul_t *ul_idx, bubble_type *bub, ma_ug_t *ug, ul_str_t *str, integer_aln_t *aln, ul_chain_t *idx, int64_t qid, integer_t *buf, 
uint64_t *cns, uint64_t cns_occ)
{
    if(idx->e<=idx->s) return 0;
    int64_t qk, tk, is_rev, tid, q_end, t_end, z; 
    ul_str_t *q_str, *t_str; uc_block_t *t_b; integer_aln_t *x, *y; 
    tid = aln[idx->s].tn_rev_qk>>33; is_rev = ((aln[idx->s].tn_rev_qk>>32)&1);
    q_str = &(str[qid]); t_str = &(str[tid]); 
    kv_resize(uint64_t, buf->u, buf->u.n + q_str->cn + t_str->cn);
    uint64_t *q_cov_buf = buf->u.a + buf->u.n, *t_cov_buf = buf->u.a + buf->u.n + q_str->cn, k, rg_occ, cn_k, mm;
    memset(q_cov_buf, -1, sizeof((*q_cov_buf))*q_str->cn); memset(t_cov_buf, -1, sizeof((*t_cov_buf))*t_str->cn);
    // if(tid != 392) return 0;
    //beg
    x = &(aln[idx->s]);
    qk = (uint32_t)(x->tn_rev_qk); tk = ((is_rev == 0)? (x->tk):(t_str->cn - x->tk - 1));

    ///direction
    if((qk > 0) && (x->tk > 0)) {
        integer_align_extention(ug, ul_idx, qid, q_str, tid, t_str, is_rev, q_cov_buf, t_cov_buf, x, 1, &q_end, &t_end);
        idx->q_sidx = q_end; idx->t_sidx = ((is_rev == 0)? (t_end):(t_str->cn - t_end - 1));
    } else {
        idx->q_sidx = (uint32_t)(x->tn_rev_qk); idx->t_sidx = x->tk;
    }
    // if(!((idx->q_sidx <= ((uint32_t)(x->tn_rev_qk))) && (idx->t_sidx <= x->tk))) {
    //     fprintf(stderr, "[M::%s] qid::%ld, tid::%ld\n", __func__, qid, tid);
    // }
    assert((idx->q_sidx <= ((uint32_t)(x->tn_rev_qk)))); assert(idx->t_sidx <= x->tk);
    // assert((is_rev && idx->t_sidx >= x->tk) || (is_rev == 0 && idx->t_sidx <= x->tk));

    //end
    x = &(aln[idx->e-1]);
    qk = (uint32_t)(x->tn_rev_qk); tk = ((is_rev == 0)? (x->tk):(t_str->cn - x->tk - 1));
    // if(tid == 316) {
    //     fprintf(stderr, "[end-M::%s::tid->%ld] qk::%ld, tk::%ld, x->tk::%u, t_str->cn::%u\n", 
    //     __func__, tid, qk, tk, x->tk, t_str->cn);
    // }
    ///direction
    if((((uint32_t)qk + 1) < q_str->cn) && ((x->tk + 1) < t_str->cn)) {
        integer_align_extention(ug, ul_idx, qid, q_str, tid, t_str, is_rev, q_cov_buf, t_cov_buf, x, 0, &q_end, &t_end);
        idx->q_eidx = q_end + 1; idx->t_eidx = ((is_rev == 0)? (t_end + 1):(t_str->cn - t_end));
    } else {
        idx->q_eidx = (uint32_t)(x->tn_rev_qk)+1; idx->t_eidx = x->tk+1;
    }

    assert(idx->q_eidx > ((uint32_t)(x->tn_rev_qk))); assert(idx->t_eidx > x->tk);
    // assert((is_rev && idx->t_eidx < x->tk) || (is_rev == 0 && idx->t_eidx > x->tk));

    assert(idx->q_eidx > idx->q_sidx); assert(idx->t_eidx > idx->t_sidx);
    // assert((is_rev && idx->t_eidx < idx->t_sidx) || (is_rev == 0 && idx->t_eidx > idx->t_sidx));
    // if(tid == 284) {
    //     fprintf(stderr, "[M::%s::tid->%ld] idx->q_sidx::%u, idx->q_eidx::%u, idx->t_sidx::%u, idx->t_eidx::%u\n", 
    //     __func__, tid, idx->q_sidx, idx->q_eidx, idx->t_sidx, idx->t_eidx);
    // }

    ///mid
    for (k = idx->s; k < idx->e; k++) {///go through all alignment pairs
        // if(qid == 3113 && tid == 3075) {
        //     fprintf(stderr, "+[M::%s] qid::%ld, tid::%ld, idx->s::%u, idx->e::%u, k::%lu, qk::%u, tk::%u, (idx->v>>1)::%u\n", 
        //     __func__, qid, tid, idx->s, idx->e, k, (uint32_t)(aln[k].tn_rev_qk), aln[k].tk, idx->v>>1);
        // }   
        
        x = &(aln[k]); y = ((k > idx->s)? (&(aln[k-1])):(NULL));
        qk = (uint32_t)(x->tn_rev_qk); 
        q_cov_buf[qk] = buf->u.a[qk]; 

        tk = ((is_rev == 0)? (x->tk):(t_str->cn - x->tk - 1));
        t_b = &(ul_idx->a[tid].bb.a[t_str->a[tk]>>32]);
        assert(((t_b->hid<<1)+t_b->rev)==((uint32_t)t_str->a[tk]));
        t_cov_buf[tk] = ug_occ_w(t_b->ts, t_b->te, &(ug->u.a[t_b->hid]));
    
        q_cov_buf[qk] += ((uint64_t)(0x8000000000000000));
        t_cov_buf[tk] += ((uint64_t)(0x8000000000000000));
        
        if(!y) continue;
        mm = 0;
        if(gap_chain_check(ug, y, x, is_rev, q_str, t_str, ul_idx->a[qid].bb.a, ul_idx->a[tid].bb.a, qid, tid, 0.08, 2000)) {
            mm = ((uint64_t)(0x8000000000000000));
        }
        for (z = ((uint32_t)(y->tn_rev_qk)) + 1; z < qk; z++) {
            q_cov_buf[z] = buf->u.a[z] + mm; 
        }

        z = ((is_rev == 0)? (y->tk):(t_str->cn - x->tk - 1)) + 1;
        tk = ((is_rev == 0)? (x->tk):(t_str->cn - y->tk - 1));
        assert(z <= tk);
        for (; z < tk; z++) {
            t_b = &(ul_idx->a[tid].bb.a[t_str->a[z]>>32]);
            assert(((t_b->hid<<1)+t_b->rev)==((uint32_t)t_str->a[z]));
            t_cov_buf[z] = ug_occ_w(t_b->ts, t_b->te, &(ug->u.a[t_b->hid])) + mm;
        }
    }

    if(cns) {
        int64_t cns_cov_occ = 0, cns_het_occ = 0, match_cns_occ = 0, match_cns_het_occ = 0, hm;
        rg_occ = idx->q_eidx;
        for (k = idx->q_sidx, cn_k = 0; k < rg_occ; k++) {
            assert((q_cov_buf[k]!=(uint64_t)-1) && (q_cov_buf[k] > 0));
            hm = ((q_cov_buf[k]<<1)>>1);
            for (; cn_k < cns_occ; cn_k++) {
                if(cns[cn_k] == k) {
                    cns_cov_occ += hm; 
                    if(q_cov_buf[k]&((uint64_t)(0x8000000000000000))) match_cns_occ += hm; 
                    if(!IF_HOM((((uint32_t)q_str->a[k])>>1), (*bub))) {
                        cns_het_occ += hm;
                        if(q_cov_buf[k]&((uint64_t)(0x8000000000000000))) match_cns_het_occ += hm; 
                    }
                } else if(cns[cn_k] > k) {
                    break;
                }
            }
        }

        if(cns_cov_occ <= 0) return 0;
        if((cns_het_occ > 0) && (match_cns_het_occ <= (cns_het_occ*0.5))) return 0;
        if(match_cns_occ <= (cns_cov_occ*0.5)) return 0;

        k = idx->t_sidx; rg_occ = idx->t_eidx; 
        if(is_rev) {
            k = t_str->cn - idx->t_eidx; rg_occ = t_str->cn - idx->t_sidx;
        }
        cns_cov_occ = match_cns_occ = 0; 
        assert(k < rg_occ);
        for (; k < rg_occ; k++) {
            assert((t_cov_buf[k]!=(uint64_t)-1) && (t_cov_buf[k] > 0));
            hm = ((t_cov_buf[k]<<1)>>1);
            cns_cov_occ += hm; 
            if(t_cov_buf[k]&((uint64_t)(0x8000000000000000))) match_cns_occ += hm; 
        }

        if(cns_cov_occ <= 0 || match_cns_occ <= 0) return 0;
        if(match_cns_occ <= (cns_cov_occ*0.5)) return 0;
    }
    ///not useful for correction
    if((idx->q_sidx + 1 == idx->q_eidx) && (idx->t_sidx + 1 == idx->t_eidx)) {///must matched with at least one cns node
        return 0;
    }

    // /**if(tid == 3074)**/ {
    //     print_integer_ovlps(ug, str, buf->b.a, buf->b.n, idx, 1, qid, buf->o.n);
    //     fprintf(stderr, "[M::%s::tid->%ld] idx->q_sidx::%u, idx->q_eidx::%u, idx->t_sidx::%u, idx->t_eidx::%u\n******************************************************\n", 
    //         __func__, tid, idx->q_sidx, idx->q_eidx, idx->t_sidx, idx->t_eidx);
    // }

    return 1;
}


void uc_block_qse_cutoff(uint32_t c_ts, uint32_t c_te, uc_block_t *x, uint32_t *nqs, uint32_t *nqe)
{
    uint32_t l, r;
    assert(c_ts >= x->ts && c_te <= x->te);
    if(!(x->rev)) {
        l = c_ts - x->ts; r = x->te - c_te;
    } else {
        r = c_ts - x->ts; l = x->te - c_te;
    }

    (*nqs) = x->qs + get_offset_adjust(l, x->te - x->ts, x->qe - x->qs);
    (*nqe) = x->qe - get_offset_adjust(r, x->te - x->ts, x->qe - x->qs);
    assert((*nqs) >= x->qs && (*nqe) <= x->qe);
}

int64_t update_exact_ul_ovlps(all_ul_t *ul_idx, ma_ug_t *ug, ul_str_t *str, integer_aln_t *aln, ul_chain_t *idx, int64_t qid, 
integer_t *buf, ul2ul_t *res)
{
    if(idx->e<=idx->s) return 0;
    int64_t qk, tk, is_rev, tid, q_end, t_end; 
    ul_str_t *q_str, *t_str; integer_aln_t *x; 
    tid = aln[idx->s].tn_rev_qk>>33; is_rev = ((aln[idx->s].tn_rev_qk>>32)&1);
    q_str = &(str[qid]); t_str = &(str[tid]); 
    uint64_t os, oe; uint32_t v, w, q_ns, q_ne, t_ns, t_ne, acc, chained, qbs, qbe, tbs, tbe;
    // if(qid == 95 || qid == 36) {
    //     fprintf(stderr,"[M::%s::qid->%ld::tid->%ld] qs_idx::%u, qe_idx::%u, ts_idx::%u, te_idx::%u\n", __func__, qid, tid, 
    //     (uint32_t)(aln[idx->s].tn_rev_qk), (uint32_t)(aln[idx->e-1].tn_rev_qk), 
    //     aln[idx->s].tk, aln[idx->e-1].tk);
    // }
    //beg
    x = &(aln[idx->s]);
    qk = (uint32_t)(x->tn_rev_qk); tk = ((is_rev == 0)? (x->tk):(t_str->cn - x->tk - 1));
    if((qk > 0) && (x->tk > 0)) return 0;
    idx->q_sidx = (uint32_t)(x->tn_rev_qk); idx->t_sidx = x->tk;

    //end
    x = &(aln[idx->e-1]);
    qk = (uint32_t)(x->tn_rev_qk); tk = ((is_rev == 0)? (x->tk):(t_str->cn - x->tk - 1));
    if((((uint32_t)qk + 1) < q_str->cn) && ((x->tk + 1) < t_str->cn)) return 0;
    idx->q_eidx = (uint32_t)(x->tn_rev_qk)+1; idx->t_eidx = x->tk+1;

    if((idx->q_eidx - idx->q_sidx) != (idx->t_eidx - idx->t_sidx)) return 0;
    if(idx->q_eidx <= idx->q_sidx) return 0;

    qk = idx->q_sidx; q_end = idx->q_eidx;
    tk = idx->t_sidx; t_end = idx->t_eidx; 
    for (; qk < q_end && tk < t_end; qk++, tk++){
        v = ((uint32_t)q_str->a[qk]); 
        if(!is_rev) w = ((uint32_t)t_str->a[tk]);
        else w = ((uint32_t)t_str->a[t_str->cn-tk-1])^1;
        if(v != w) break;
    }
    if(qk!=q_end || tk!=t_end) return 0;
    
    uc_block_t *qi, *ti;
    qk = idx->q_sidx; q_end = idx->q_eidx;
    tk = idx->t_sidx; t_end = idx->t_eidx; 

    buf->p.n = buf->f.n = buf->o.n = buf->u.n = 0;
    kv_resize(int64_t, buf->p, idx->q_eidx - idx->q_sidx); int64_t *p = buf->p.a;
    kv_resize(int64_t, buf->f, idx->q_eidx - idx->q_sidx); int64_t *f = buf->f.a;
    kv_resize(uint64_t, buf->o, idx->q_eidx - idx->q_sidx); uint64_t *qc = buf->o.a;
    kv_resize(uint64_t, buf->u, idx->q_eidx - idx->q_sidx); uint64_t *tc = buf->u.a;
    int64_t max_f, max_z, csc, ps, sc, tf, tz, z;
    for (acc = 0, chained = 1, tf = tz = -1; qk < q_end && tk < t_end; qk++, tk++) {
        qi = &(ul_idx->a[qid].bb.a[q_str->a[qk]>>32]); 
        if(!is_rev) ti = &(ul_idx->a[tid].bb.a[t_str->a[tk]>>32]);
        else ti = &(ul_idx->a[tid].bb.a[t_str->a[t_str->cn-tk-1]>>32]);
        assert(qi->hid == ti->hid);
        os = MAX(qi->ts, ti->ts); oe = MIN(qi->te, ti->te);
        if(oe <= os) continue;
        uc_block_qse_cutoff(os, oe, qi, &q_ns, &q_ne);
        uc_block_qse_cutoff(os, oe, ti, &t_ns, &t_ne);
        qc[acc] = q_ns; qc[acc] <<= 32; qc[acc] |= q_ne;
        tc[acc] = t_ns; tc[acc] <<= 32; tc[acc] |= t_ne;

        z = acc; z -= 1; csc = q_ne - q_ns; if(csc <= 0) csc = 1;
        max_f = csc; max_z = -1; ps = 0;
        if(chained && z >= 0) {
            qbs = qc[z]>>32; qbe = (uint32_t)qc[z]; 
            tbs = tc[z]>>32; tbe = (uint32_t)tc[z];
            if(!is_rev) {
                if(qbs <= q_ns && qbe <= q_ne && tbs <= t_ns && tbe <= t_ne) ps = 1;
            } else {
                if(qbs <= q_ns && qbe <= q_ne && tbs >= t_ns && tbe >= t_ne) ps = 1;
            }
            if(ps) {
                sc = csc + f[z];
                if(sc > max_f) {
                    max_f = sc; max_z = z;
                }
            } else {
                chained = 0;
            }
        }

        if((!ps) || (!chained)) {
            for (; z >= 0; z--) {
                qbs = qc[z]>>32; qbe = (uint32_t)qc[z]; 
                tbs = tc[z]>>32; tbe = (uint32_t)tc[z];
                ps = 0;
                if(!is_rev) {
                    if(qbs <= q_ns && qbe <= q_ne && tbs <= t_ns && tbe <= t_ne) ps = 1;
                } else {
                    if(qbs <= q_ns && qbe <= q_ne && tbs >= t_ns && tbe >= t_ne) ps = 1;
                }
                if(!ps) continue;
                sc = csc + f[z];
                if(sc > max_f) {
                    max_f = sc; max_z = z;
                }
            }
        }

        f[acc] = max_f; p[acc] = max_z;
        if(tf < max_f) {
            tf = max_f; tz = acc;
        }
        acc++;
    }
    if(acc <= 0) return 0;
    assert(tz >= 0);
    res->hid = res->qs = res->qe = res->ts = res->te = (uint32_t)-1;
    res->hid = tid; res->is_rev = !!is_rev; res->is_del = 0; res->is_ct = 0;
    res->qs_k = idx->q_sidx; res->qe_k = idx->q_eidx;
    if(!is_rev) {
        res->ts_k = idx->t_sidx; res->te_k = idx->t_eidx;
    } else {
        res->ts_k = t_str->cn - idx->t_eidx; res->te_k = t_str->cn - idx->t_sidx;
    }
    

    res->qe = (uint32_t)qc[tz];
    if(!is_rev) res->te = (uint32_t)tc[tz];
    else res->ts = tc[tz]>>32;
    
    for (z = tz; z >= 0; z = p[z]) {
        res->qs = qc[z]>>32; 
        if(!is_rev) res->ts = tc[z]>>32;
        else res->te = (uint32_t)tc[z];
    }
    assert(res->qs <= res->qe && res->ts <= res->te);


    int64_t qs = 0, qe = 0, rs = 0, re = 0, qtail = 0, rtail = 0; 
    qs = res->qs; qe = res->qe; rs = res->ts; re = res->te;
    if(is_rev) {
        rs = ul_idx->a[tid].rlen - res->te; re = ul_idx->a[tid].rlen - res->ts;
    }
    if(qs <= rs) {
        rs -= qs; qs = 0;
    } else {
        qs -= rs; rs = 0;
    }

    qtail = ul_idx->a[qid].rlen - qe; rtail = ul_idx->a[tid].rlen - re;
    if(qtail <= rtail) {
        qe = ul_idx->a[qid].rlen; re += qtail;
    }
    else {
        re = ul_idx->a[tid].rlen; qe += rtail; 
    }
    res->qs = qs; res->qe = qe;
    res->ts = rs; res->te = re;
    if(is_rev) {
        res->ts = ul_idx->a[tid].rlen - re; 
        res->te = ul_idx->a[tid].rlen - rs;
    }
    return 1;
}



void reset_poa_g_t(poa_g_t *g)
{
    g->seq.n = g->arc.n = g->idx.n = 0; g->update_arc = g->update_seq = g->e_idx.n = 0;
}

void clean_poa_g_t(poa_g_t *g)
{
    uint64_t set_n = g->seq.n<<1;
    assert(g->arc.n >= g->update_arc); assert(g->seq.n >= g->update_seq);
    if(g->seq.n > g->update_seq) {
        kv_resize(uint64_t, g->idx, (g->seq.n<<1)); set_n = g->update_seq<<1;
        memset(g->idx.a + (g->update_seq<<1), 0, sizeof((*g->idx.a))*((g->seq.n<<1)-(g->update_seq<<1)));
        g->update_seq = g->seq.n; g->idx.n = (g->seq.n<<1);
    }

    if(g->arc.n > g->update_arc) {
        radix_sort_poa_arc_srt(g->arc.a, g->arc.a + g->arc.n);
        memset(g->idx.a, 0, sizeof((*g->idx.a))*set_n);
        int64_t k, last, n = g->arc.n;
        for (k = 1, last = 0; k <= n; ++k){
            if (k == n || g->arc.a[k-1].ul>>32 != g->arc.a[k].ul>>32) {
                g->idx.a[g->arc.a[k-1].ul>>32] = (((uint64_t)last<<32)) | (k - last), last = k;
            }
        }
        g->update_arc = g->arc.n;
    }
}

void append_unmatch_integer_seq(poa_g_t *g, ma_ug_t *ug, uc_block_t *raw, ul_str_t *str, int64_t s, int64_t e, int64_t is_rev, int64_t str_id)
{
    int64_t k; poa_nid_t *nn; poa_arc_t *ae; uint32_t v; uc_block_t *z; emap_t *em;
    for (k = s; k < e; k++) {
        if(is_rev) {
            v = (((uint32_t)str->a[str->cn-k-1])^1);
            z = &(raw[str->a[str->cn-k-1]>>32]);
        }
        else {
            v = ((uint32_t)str->a[k]); z = &(raw[str->a[k]>>32]);
        }
        kv_pushp(poa_nid_t, g->seq, &nn);
        nn->nid = v; nn->occ = ug_occ_w(z->ts, z->te, &(ug->u.a[z->hid])); 
        
        // ug->u.a[v>>1].n;
        if(k > s) {
            kv_pushp(poa_arc_t, g->arc, &ae);
            ae->ul = g->seq.n-1; ae->ul <<= 33; ae->ul += ((uint64_t)(0x100000000)); ae->ul += 1;
            ae->v = g->seq.n-2; ae->v <<= 1; ae->v += 1;

            kv_pushp(poa_arc_t, g->arc, &ae);
            ae->ul = g->seq.n-2; ae->ul <<= 33; ae->ul += 1;
            ae->v = g->seq.n-1; ae->v <<= 1;

            kv_pushp(emap_t, g->e_idx, &em);
            em->pge = g->seq.n-2; em->pge <<= 32; em->pge += g->seq.n-1;
            em->ule = (is_rev?(str->cn-k):(k-1)); em->ule <<= 32; em->ule += (is_rev?(str->cn-k-1):(k));
            em->ulid = str_id;
            // assert((g->seq.a[em->pge>>32].nid^(is_rev?1:0)) == ((uint32_t)str->a[em->ule>>32]));
            // assert((g->seq.a[(uint32_t)em->pge].nid^(is_rev?1:0)) == ((uint32_t)str->a[(uint32_t)em->ule]));
        }
    }
}

void update_poa_nid_occ(ma_ug_t *ug, uc_block_t *raw, poa_g_t *g, int64_t gidx, uint64_t *str, int64_t str_idx, int64_t is_rev)
{
    uint32_t g_v, str_v, new_occ; uc_block_t *z; 
    ///update 
    g_v = g->seq.a[gidx].nid;
    z = &(raw[str[str_idx]>>32]); str_v = ((uint32_t)str[str_idx]); if(is_rev) str_v ^= 1;
    assert(g_v == str_v); assert(g->seq.a[gidx].occ <= ug->u.a[g->seq.a[gidx].nid>>1].n);
    if(g->seq.a[gidx].occ != ug->u.a[g->seq.a[gidx].nid>>1].n) {
        new_occ = ug_occ_w(z->ts, z->te, &(ug->u.a[z->hid]));
        if(new_occ > g->seq.a[gidx].occ) g->seq.a[gidx].occ = new_occ;
    }
}

void insert_poa_nodes_0(ma_ug_t *ug, uc_block_t *raw, poa_g_t *g, uint64_t *str, int64_t str_occ, uint64_t is_rev, int64_t str_id, int64_t str_off)
{
    int64_t k; poa_nid_t *nn; poa_arc_t *ae; uint32_t v; uc_block_t *z; emap_t *em;
    for (k = 0; k < str_occ; k++) {
        if(is_rev == 0) {
            v = ((uint32_t)str[k]); z = &(raw[str[k]>>32]);
        } else {
            v = ((uint32_t)str[str_occ-k-1])^1; z = &(raw[str[str_occ-k-1]>>32]);
        }
        
        kv_pushp(poa_nid_t, g->seq, &nn);
        nn->nid = v; nn->occ = ug_occ_w(z->ts, z->te, &(ug->u.a[z->hid]));

        if(k > 0) {
            kv_pushp(poa_arc_t, g->arc, &ae);
            ae->ul = g->seq.n-1; ae->ul <<= 33; ae->ul += ((uint64_t)(0x100000000)); ae->ul += 1;
            ae->v = g->seq.n-2; ae->v <<= 1; ae->v += 1;

            kv_pushp(poa_arc_t, g->arc, &ae);
            ae->ul = g->seq.n-2; ae->ul <<= 33; ae->ul += 1;
            ae->v = g->seq.n-1; ae->v <<= 1;

            kv_pushp(emap_t, g->e_idx, &em);
            em->pge = g->seq.n-2; em->pge <<= 32; em->pge += g->seq.n-1;
            em->ule = str_off + (is_rev?(str_occ-k):(k-1)); em->ule <<= 32; 
            em->ule += str_off + (is_rev?(str_occ-k-1):(k));
            em->ulid = str_id;

            // if(!((g->seq.a[em->pge>>32].nid^(is_rev?1:0)) == ((uint32_t)debug_str->a[em->ule>>32]))) {
            //     fprintf(stderr, "[M::%s::k->%ld::str_occ->%ld] is_rev->%lu, pg_v->%u, str_v->%u, str_off->%ld, str_k->%lu, str_cn->%u, address_diff->%u\n", 
            //             __func__, k, str_occ, is_rev, (g->seq.a[em->pge>>32].nid^(is_rev?1:0)), 
            //             ((uint32_t)debug_str->a[em->ule>>32]), str_off, em->ule>>32, debug_str->cn, (uint32_t)(str - debug_str->a));
            //     uint64_t debug_k;
            //     for (debug_k = 0; debug_k < debug_str->cn; debug_k++) {
            //         fprintf(stderr, "[M::%s::debug_k->%lu] str_v->%u\n", 
            //             __func__, debug_k, ((uint32_t)debug_str->a[debug_k]));
            //     }
                
            // }
            // assert((g->seq.a[em->pge>>32].nid^(is_rev?1:0)) == ((uint32_t)debug_str->a[em->ule>>32]));
            // assert((g->seq.a[(uint32_t)em->pge].nid^(is_rev?1:0)) == ((uint32_t)debug_str->a[(uint32_t)em->ule]));
        }
    }
}

void push_poa_arch_0(poa_g_t *g, uint32_t src, uint32_t des, int64_t str_id, int64_t str_src, int64_t str_des)
{
    poa_arc_t *ae; emap_t *em;
    kv_pushp(poa_arc_t, g->arc, &ae);
    ae->ul = des; ae->ul <<= 33; ae->ul += ((uint64_t)(0x100000000)); ae->ul += 1;
    ae->v = src; ae->v <<= 1; ae->v += 1;

    kv_pushp(poa_arc_t, g->arc, &ae);
    ae->ul = src; ae->ul <<= 33; ae->ul += 1;
    ae->v = des; ae->v <<= 1;

    kv_pushp(emap_t, g->e_idx, &em);
    em->pge = src; em->pge <<= 32; em->pge += des;
    em->ule = str_src; em->ule <<= 32; em->ule += str_des;
    em->ulid = str_id;

    // assert((g->seq.a[em->pge>>32].nid^(debug_is_rev?1:0)) == ((uint32_t)debug_str->a[em->ule>>32]));
    // assert((g->seq.a[(uint32_t)em->pge].nid^(debug_is_rev?1:0)) == ((uint32_t)debug_str->a[(uint32_t)em->ule]));
}

void update_poa_arch_0(poa_g_t *g, uint32_t src, uint32_t des, int64_t str_id, int64_t str_src, int64_t str_des)
{
    uint32_t k, v, w, a_n; poa_arc_t *a; emap_t *em;
    v = src<<1; w = des<<1;
    a_n = poa_arc_n(g, v); a = poa_arc_a(g, v);
    for (k = 0; k < a_n; k++) {
        if(a[k].v == w) break;
    }
    if(k >= a_n) {
        push_poa_arch_0(g, src, des, str_id, str_src, str_des);
    } else {
        a[k].ul++;
        v = des<<1; v^=1; 
        w = src<<1; w^=1;
        a_n = poa_arc_n(g, v); a = poa_arc_a(g, v);
        for (k = 0; k < a_n; k++) {
            if(a[k].v == w) break;
        }
        assert(k < a_n);
        a[k].ul++;

        kv_pushp(emap_t, g->e_idx, &em);
        em->pge = src; em->pge <<= 32; em->pge += des;
        em->ule = str_src; em->ule <<= 32; em->ule += str_des;
        em->ulid = str_id;

        // assert((g->seq.a[em->pge>>32].nid^(debug_is_rev?1:0)) == ((uint32_t)debug_str->a[em->ule>>32]));
        // assert((g->seq.a[(uint32_t)em->pge].nid^(debug_is_rev?1:0)) == ((uint32_t)debug_str->a[(uint32_t)em->ule]));
    }
}

void append_integer_seq_frag(ma_ug_t *ug, uc_block_t *raw, poa_g_t *g, int64_t g_beg, int64_t g_end, uint64_t *str, int64_t str_occ, uint64_t is_rev, int64_t str_id, int64_t str_off)
{
    // if(str_occ <= 0) return;
    uint32_t nid; 
    // if((g_beg >= 0 || g_end >= 0) && str_occ < 2) return;
    if(g_beg < 0 && g_end < 0) {//add new nodes
        insert_poa_nodes_0(ug, raw, g, str, str_occ, is_rev, str_id, str_off);
        return;
    }

    if(g_beg < 0 && g_end >= 0) {///add nodes to the left end
        assert(str_occ >= 1);
        update_poa_nid_occ(ug, raw, g, g_end, str, (is_rev?(0):(str_occ-1)), is_rev);
        if(str_occ < 2) return;
        insert_poa_nodes_0(ug, raw, g, (is_rev?(str+1):(str)), str_occ-1, is_rev, str_id, str_off+(is_rev?(1):(0)));
        
        push_poa_arch_0(g, g->seq.n-1, g_end, str_id, str_off + (is_rev?(1):(str_occ-2)), 
                                                                str_off + (is_rev?(0):(str_occ-1)));
        return;
    }

    if(g_beg >= 0 && g_end < 0) {///add nodes to the right end
        assert(str_occ >= 1);
        update_poa_nid_occ(ug, raw, g, g_beg, str, (is_rev?(str_occ-1):(0)), is_rev);
        if(str_occ < 2) return;
        nid = g->seq.n;///backup
        insert_poa_nodes_0(ug, raw, g, (is_rev?(str):(str+1)), str_occ-1, is_rev, str_id, str_off+(is_rev?(0):(1)));
        push_poa_arch_0(g, g_beg, nid, str_id, str_off + (is_rev?(str_occ-1):(0)), 
                                                        str_off + (is_rev?(str_occ-2):(1)));
        return;
    }

    if(g_beg >= 0 && g_end >= 0) {///add nodes to the middle
        assert(str_occ >= 2);
        update_poa_nid_occ(ug, raw, g, g_beg, str, (is_rev?(str_occ-1):(0)), is_rev);
        update_poa_nid_occ(ug, raw, g, g_end, str, (is_rev?(0):(str_occ-1)), is_rev);
        if(str_occ > 2) {///insert new nodes
            nid = g->seq.n;///backup
            insert_poa_nodes_0(ug, raw, g, str+1, str_occ-2, is_rev, str_id, str_off+1);
            push_poa_arch_0(g, g_beg, nid, str_id, str_off + (is_rev?(str_occ-1):(0)),
                                                            str_off + (is_rev?(str_occ-2):(1))); 
            push_poa_arch_0(g, g->seq.n-1, g_end, str_id, str_off + (is_rev?(1):(str_occ-2)),
                                                            str_off + (is_rev?(0):(str_occ-1)));
        } else {
            ///add an edge between g_beg and g_end
            update_poa_arch_0(g, g_beg, g_end, str_id, str_off + (is_rev?(str_occ-1):(0)), 
                                                        str_off + (is_rev?(0):(str_occ-1)));
        }
    }
}

void append_aligned_integer_seq(ma_ug_t *ug, uc_block_t *raw, poa_g_t *g, int64_t g_occ, uint64_t *str, int64_t str_occ, uint64_t is_rev,
uint32_t *match_g, uint32_t *match_str, int64_t match_occ)
{
    if(str_occ <= 0) return;
    int64_t k, p_str, p_g;
    
    if(is_rev == 0) {
        for (k = match_occ-1, p_str = 0, p_g = -1; k >= 0; k--) {
            fprintf(stderr, "+[M::%s::] match_str[%ld]::%u, match_occ::%ld, p_g::%ld\n", 
            __func__, k, match_str[k], match_occ, p_g);
            assert(((int64_t)match_str[k]) >= p_str); 
            // append_integer_seq_frag(ug, raw, g, p_g, match_g[k], str+p_str, match_str[k]+1-p_str, is_rev);
            p_str = match_str[k]; p_g = match_g[k]; 
        }
        // append_integer_seq_frag(ug, raw, g, p_g, -1, str+p_str, str_occ-p_str, is_rev);
    } else {
        for (k = match_occ-1, p_str = str_occ-1, p_g = -1; k >= 0; k--) {
            fprintf(stderr, "-[M::%s::] match_str[%ld]::%u, match_occ::%ld\n", 
            __func__, k, match_str[k], match_occ);
            assert(((int64_t)match_str[k]) <= p_str); 
            // append_integer_seq_frag(ug, raw, g, p_g, match_g[k], str+match_str[k], p_str+1-match_str[k], is_rev);
            p_str = match_str[k]; p_g = match_g[k];  
        }
        // append_integer_seq_frag(ug, raw, g, p_g, -1, str, p_str+1, is_rev);
    }
}


#define poa_str_idx(i, occ, is_rev) (((is_rev))?((occ)-(i)-1):(i))
void append_aligned_integer_seq_by_aln_pair(ma_ug_t *ug, uc_block_t *raw, poa_g_t *g, int64_t g_occ, uint64_t *str, int64_t str_occ, uint64_t is_rev, integer_aln_t *a, int64_t a_n, int64_t str_id, int64_t str_off)
{
    if(str_occ <= 0 || a_n <= 0) return;
    int64_t k, p_str, p_g; uint64_t *qstr_a, qstr_n, qoff; uint32_t *gidx = g->srt_b.res.a;

    for (k = 0, p_g = -1, p_str = (is_rev?(str_occ-1):(0)); k < a_n; k++) {
        if(!is_rev) {
            qstr_a = str+p_str; qstr_n = ((uint32_t)a[k].tn_rev_qk)+1-p_str; qoff = p_str + str_off;
        } else {
            qstr_a = str+poa_str_idx(((uint32_t)a[k].tn_rev_qk), str_occ, is_rev);
            qstr_n = p_str+1-(poa_str_idx(((uint32_t)a[k].tn_rev_qk), str_occ, is_rev));
            qoff = poa_str_idx(((uint32_t)a[k].tn_rev_qk), str_occ, is_rev) + str_off;
        }
        // if(str_id == 47072) {
        //     fprintf(stderr, "+[M::%s::k->%ld] p_g::%ld, c_g::%u, p_str::%ld, c_str::%ld, str_off::%ld, qoff::%lu, qstr_n::%lu\n", 
        //         __func__, k, p_g, gidx[a[k].tk], p_str, poa_str_idx(((uint32_t)a[k].tn_rev_qk), str_occ, is_rev), 
        //         str_off, qoff, qstr_n);
        // }

        append_integer_seq_frag(ug, raw, g, p_g, gidx[a[k].tk], qstr_a, qstr_n, is_rev, str_id, qoff);

        p_str = (uint32_t)a[k].tn_rev_qk;
        if(is_rev) p_str = poa_str_idx(p_str, str_occ, is_rev);
        p_g = gidx[a[k].tk];
    }
    if(!is_rev) {
        qstr_a = str+p_str; qstr_n = str_occ-p_str; qoff = p_str + str_off;
    } else {
        qstr_a = str; qstr_n = p_str+1; qoff = str_off;
    }
    append_integer_seq_frag(ug, raw, g, p_g, -1, qstr_a, qstr_n, is_rev, str_id, qoff);
}


uint32_t topo_srt_gen(poa_g_t *g, uint32_t debug_qid, ma_ug_t *ug)
{
    uint32_t k, v, w, a_n; poa_arc_t *a;
    kv_resize(uint32_t, g->srt_b.ind, g->seq.n);
    kv_resize(uint32_t, g->srt_b.stack, g->seq.n);
    kv_resize(uint32_t, g->srt_b.res, g->seq.n);
    kv_resize(uint32_t, g->srt_b.res2nid, g->seq.n);
    // kv_resize(uint64_t, g->srt_b.aln, g->seq.n); g->srt_b.aln.n = 0;
    g->srt_b.ind.n = g->srt_b.stack.n = g->srt_b.res.n = g->srt_b.res2nid.n = 0;

    for (k = 0; k < g->seq.n; k++) {
        v = (k<<1) + 1; 
        g->srt_b.ind.a[k] = poa_arc_n(g, v);
        if(g->srt_b.ind.a[k] == 0) kv_push(uint32_t, g->srt_b.stack, k);
    }
    
    while (g->srt_b.stack.n > 0) {
        v = g->srt_b.stack.a[--g->srt_b.stack.n]; kv_push(uint32_t, g->srt_b.res, v); v <<= 1;
        a_n = poa_arc_n(g, v); a = poa_arc_a(g, v);
        for (k = 0; k < a_n; k++) {
            w = a[k].v>>1;
            g->srt_b.ind.a[w]--;
            if(g->srt_b.ind.a[w] == 0) kv_push(uint32_t, g->srt_b.stack, w);
        }
    }
    // if(!(g->srt_b.res.n == g->seq.n)) {
    //     fprintf(stderr, "[M::%s::] debug_qid::%u, g->seq.n::%u, g->srt_b.res.n::%u\n", 
    //     __func__, debug_qid, (uint32_t)g->seq.n, (uint32_t)g->srt_b.res.n);
    //     print_integer_g(g, ug, 1);
    //     for (k = 0; k < g->seq.n; k++) {
    //         if(g->srt_b.ind.a[k] > 0) fprintf(stderr, "circle->nid::%u\n", k);
    //     }
    // }
    if(g->srt_b.res.n != g->seq.n) return 0;///there is a circle
    assert(g->srt_b.res.n == g->seq.n);
    for (k = 0; k < g->seq.n; k++) {
        g->srt_b.res2nid.a[g->srt_b.res.a[k]] = k;
        // g->srt_b.aln.a[k] = g->seq.a[g->srt_b.res.a[k]].nid;
        // g->srt_b.aln.a[k] <<= 32; g->srt_b.aln.a[k] += k;
    }
    // radix_sort_srt64(g->srt_b.aln.a, g->srt_b.aln.a + g->srt_b.aln.n);
    return 1;
}

void init_poa_dp(ma_ug_t *ug, poa_dp_t *dp, poa_g_t *g, uint64_t g_occ, uint64_t *str, uint64_t str_occ, uint64_t is_rev, uc_block_t *raw, integer_t *buf)
{
    kv_resize(uint8_t, dp->dir, (str_occ+1)*(g_occ+1));
    kv_resize(int64_t, dp->sc, (str_occ+1)*(g_occ+1));
    kv_resize(uint64_t, dp->prefix, (str_occ+1)*(g_occ+1));
    dp->n = str_occ+1; dp->m = g_occ+1;
    kv_resize(uint64_t, buf->u, str_occ);
    int64_t *sc = dp->sc.a, bsc, ss; uint8_t *dir = dp->dir.a; uint64_t k, l, *str_w = buf->u.a, *prefix = dp->prefix.a, m; 
    uc_block_t *z; uint32_t *g_idx = g->srt_b.res.a, *n2gidx = g->srt_b.res2nid.a, v, w, a_n, bsc_i; poa_arc_t *a; 
    for (k = 0; k < str_occ; k++) {
        z = &(raw[str[k]>>32]); str_w[k] = ug_occ_w(z->ts, z->te, &(ug->u.a[z->hid]));
    }
    
    sc[poa_dp_idx(*dp, 0, 0)] = dir[poa_dp_idx(*dp, 0, 0)] = 0;
    for (k = 1, l = 0; k < dp->n; k++) {///pat; new ul read
        l += str_w[poa_str_idx(k-1, str_occ, is_rev)];
        sc[poa_dp_idx(*dp, k, 0)] = l;
        dir[poa_dp_idx(*dp, k, 0)] = lstr_dp;
        prefix[poa_dp_idx(*dp, k, 0)] = ((k - 1)<<32);
    }
    
    for (k = 1; k < dp->m; k++) {///ref; graph
        v = (g_idx[k-1]<<1) + 1;
        a_n = poa_arc_n(g, v); a = poa_arc_a(g, v);
        for (m = 0, bsc = bsc_i = 0; m < a_n; m++) {
            w = n2gidx[a[m].v>>1]; assert(w+1 < k);
            ss = sc[poa_dp_idx(*dp, 0, w+1)];
            if(ss < bsc || bsc_i == 0) {
                bsc = ss; bsc_i = w + 1;
            }
        }
        sc[poa_dp_idx(*dp, 0, k)] = bsc + g->seq.a[v>>1].occ;
        dir[poa_dp_idx(*dp, 0, k)] = lref_dp;
        prefix[poa_dp_idx(*dp, 0, k)] = bsc_i;
    }
}

uint32_t update_poa_dp(poa_g_t *g, uint32_t debug_qid, ma_ug_t *debug_ug)
{
    uint32_t is_srt = 0/**, is_up_aln = 0**/, k, is_circle = 0;
    if(g->seq.n > g->update_seq || g->arc.n > g->update_arc) {
        if(g->seq.n > g->update_seq) {
            kv_resize(uint32_t, g->srt_b.res, g->seq.n);
            kv_resize(uint32_t, g->srt_b.res2nid, g->seq.n);
            // kv_resize(uint64_t, g->srt_b.aln, g->seq.n); g->srt_b.aln.n = g->seq.n;
            g->srt_b.res.n = g->srt_b.res2nid.n = g->seq.n;
            for (k = g->update_seq; k < g->seq.n; k++) {
                g->srt_b.res.a[k] = g->srt_b.res2nid.a[k] = k;
                // g->srt_b.aln.a[k] = (((uint64_t)(g->seq.a[k].nid))<<32)+k;
                // if(k > 0 && is_up_aln == 0) {
                //     if((g->srt_b.aln.a[k]>>32) < (g->srt_b.aln.a[k-1]>>32)) is_up_aln = 1;
                // }
            }
        }
        if(g->arc.n > g->update_arc) { ///check if it is necessary to resort
            for (k = g->update_arc; k < g->arc.n; k++) {
                if((g->arc.a[k].ul>>32)&1) continue;
                if(g->srt_b.res2nid.a[g->arc.a[k].ul>>33] >= g->srt_b.res2nid.a[g->arc.a[k].v>>1]) break;
            }
            if(k < g->arc.n) is_srt = 1;
        }
        // fprintf(stderr, "[M::%s::] is_srt::%u\n", __func__, is_srt);
        clean_poa_g_t(g); 
        if(is_srt) {
            is_circle = 1 - topo_srt_gen(g, debug_qid, debug_ug);
        } 
        // else if(is_up_aln) {
        //     radix_sort_srt64(g->srt_b.aln.a, g->srt_b.aln.a + g->srt_b.aln.n);
        // }
    }

    return is_circle;
}

/**
void poa_dp(poa_g_t *g, ma_ug_t *ug, uc_block_t *raw, ul_str_t *str, int64_t s, int64_t e, int64_t is_rev, integer_t *buf)
{
    if(e <= s) return;

    

    uint32_t *g_idx = g->srt_b.res.a, *n2gidx = g->srt_b.res2nid.a, gnid, v, a_n, pat_v, g_v;
    uint64_t *pat = (is_rev?(str->a + str->cn - e):(str->a + s)), pat_n = e - s, pp;
    init_poa_dp(ug, &(g->dp), g, g->seq.n, pat, pat_n, is_rev, raw, buf);
    int64_t *sc = g->dp.sc.a, min_w, c_w; uint8_t *dir = g->dp.dir.a; 
    uint64_t *str_w = buf->u.a, *prefix = g->dp.prefix.a;
    poa_arc_t *a; poa_dp_t *dp = &(g->dp); 
    int64_t min_i, min_k, min_d, i, k, n = g->dp.n - 1, m = g->dp.m - 1, match_sc, z, pidx, g_k, pat_k;
    fprintf(stderr, "[M::%s::] ts::%ld, te::%ld, is_rev::%ld, m::%ld(g->seq.n::%u), n::%ld(pat_n::%lu)\n", 
    __func__, s, e, is_rev, m, (uint32_t)g->seq.n, n, pat_n);
    for (i = 0; i < m; i++) {///graph
        gnid = g_idx[i]; v = (gnid<<1) + 1; g_v = g->seq.a[gnid].nid;
        a_n = poa_arc_n(g, v); a = poa_arc_a(g, v);
        for (k = 0; k < n; k++) {//read 
            pat_v = (uint32_t)pat[poa_str_idx(k, pat_n, is_rev)]; if(is_rev) pat_v ^= 1;
            if(g_v == pat_v) {
                match_sc = str_w[poa_str_idx(k, pat_n, is_rev)]; match_sc *= -1;
            } else {
                match_sc = str_w[poa_str_idx(k, pat_n, is_rev)];
                if(match_sc < g->seq.a[gnid].occ) match_sc = g->seq.a[gnid].occ;
            }
            
            ///longer graph/shorter read
            min_w = sc[poa_dp_idx(*dp, k, i+1)] + g->seq.a[gnid].occ; 
            min_i = i + 1; min_k = k; min_d = lg_dp;

            for (z = 0; z < a_n; z++) {
                pidx = n2gidx[a[z].v>>1]; assert(pidx < i);
                ///match
                c_w = sc[poa_dp_idx(*dp, k, pidx+1)] + match_sc;
                if(c_w < min_w) {
                    min_w = c_w; min_i = pidx+1; min_k = k; 
                    if(match_sc <= 0) min_d = e_pdp;
                    else min_d = ue_pdp;
                }
                ///longer read/shorter graph
                c_w = sc[poa_dp_idx(*dp, k+1, pidx+1)] + str_w[poa_str_idx(k, pat_n, is_rev)];
                if(c_w < min_w) {
                    min_w = c_w; min_i = pidx+1; min_k = k + 1; min_d = lstr_dp;
                }
            }

            if(a_n == 0) {///no prefix
                pidx = -1;
                ///match
                c_w = sc[poa_dp_idx(*dp, k, pidx+1)] + match_sc;
                if(c_w < min_w) {
                    min_w = c_w; min_i = pidx+1; min_k = k; 
                    if(match_sc <= 0) min_d = e_pdp;
                    else min_d = ue_pdp;
                }
                ///longer read/shorter graph
                c_w = sc[poa_dp_idx(*dp, k+1, pidx+1)] + str_w[poa_str_idx(k, pat_n, is_rev)];
                if(c_w < min_w) {
                    min_w = c_w; min_i = pidx+1; min_k = k + 1; min_d = lstr_dp;
                }
            }

            sc[poa_dp_idx(*dp, k+1, i+1)] = min_w; dir[poa_dp_idx(*dp, k+1, i+1)] = min_d;
            prefix[poa_dp_idx(*dp, k+1, i+1)] = (((uint64_t)min_k)<<32)|((uint64_t)min_i);
            fprintf(stderr, "[M::%s::] i::%ld(graph->utg%.6d%c, min_i->%ld), k::%ld(str->utg%.6d%c, min_k->%ld), match_sc::%ld, min_w::%ld, min_d::%ld\n", 
                                            __func__, i+1, (int32_t)(g_v>>1)+1, "lc"[ug->u.a[g_v>>1].circ], min_i,
                                            k+1, (int32_t)(pat_v>>1)+1, "lc"[ug->u.a[pat_v>>1].circ], min_k, match_sc, min_w, min_d);
        }        
    }

    ///backtrack; global alignment
    min_i = -1; min_k = n; min_w = 0;
    for (i = 0; i < m; i++) {///go through graph
        // gnid = g_idx[i]; v = (gnid<<1);
        // if(poa_arc_n(g, v) > 0) continue;
        c_w = sc[poa_dp_idx(*dp, min_k, i+1)];
        if(min_i < 0 || c_w < min_w) {
            min_w = c_w; min_i = i+1;
        }
    }
    assert(min_i > 0);

    g->srt_b.ind.n = g->srt_b.stack.n = 0;
    while (min_i > 0 || min_k > 0) {
        pat_k = poa_str_idx((min_k-1), pat_n, is_rev);///read
        g_k = g_idx[min_i-1];///graph
        fprintf(stderr, "******[M::%s::] min_i::%ld(gid->%ld, m->%ld), min_k::%ld(str_id->%ld, n->%ld), dir::%u, sc::%ld\n", 
                                            __func__, min_i, g_k, m, min_k, pat_k, n,  
                                            dir[poa_dp_idx(*dp, min_k, min_i)], sc[poa_dp_idx(*dp, min_k, min_i)]);
        if(dir[poa_dp_idx(*dp, min_k, min_i)] == e_pdp) {
            assert(g->seq.a[g_k].nid == (((uint32_t)pat[pat_k])^(is_rev?1:0)));
            kv_push(uint32_t, g->srt_b.ind, pat_k); ///read
            kv_push(uint32_t, g->srt_b.stack, g_k); ///graph
        }
        pp = prefix[poa_dp_idx(*dp, min_k, min_i)];
        assert((int64_t)(pp>>32)<=min_k); assert((int64_t)((uint32_t)pp)<=min_i);
        min_k = pp >> 32; min_i = (uint32_t)pp;
    }
    if(!(g->srt_b.ind.n > 0 && g->srt_b.ind.n == g->srt_b.stack.n)) {
        fprintf(stderr, "[M::%s::] g->srt_b.ind.n::%u, g->srt_b.stack.n::%u, ts::%ld, te::%ld\n", 
                                            __func__, (uint32_t)g->srt_b.ind.n, (uint32_t)g->srt_b.stack.n, s, e);
    }
    assert(g->srt_b.ind.n > 0 && g->srt_b.ind.n == g->srt_b.stack.n);
    append_aligned_integer_seq(ug, raw, g, g->seq.n, pat, pat_n, is_rev, g->srt_b.stack.a, g->srt_b.ind.a, g->srt_b.ind.n);
    update_poa_dp(g);
}


void poa_cns_dp(poa_g_t *g, all_ul_t *ul_idx, ma_ug_t *ug, ul_str_t *str, integer_aln_t *aln, 
ul_chain_t *idx, int64_t idx_n, int64_t qid, integer_t *buf)
{
    int64_t k, tid, is_rev;
    reset_poa_g_t(g);

    append_unmatch_integer_seq(g, ug, ul_idx->a[qid].bb.a, &(str[qid]), 0, str[qid].cn, 0, qid, 0);
    clean_poa_g_t(g); topo_srt_gen(g);

    for (k = 0; k < idx_n; k++) {
        tid = aln[idx[k].s].tn_rev_qk>>33; is_rev = ((aln[idx[k].s].tn_rev_qk>>32)&1);
        fprintf(stderr, "\n[M::%s::] k::%ld, tid::%ld, is_rev::%ld\n", __func__, k, tid, is_rev);
        print_integer_seq(ug, str, tid, 1);
        poa_dp(g, ug, ul_idx->a[tid].bb.a, &(str[tid]), idx[k].t_sidx, idx[k].t_eidx, is_rev, buf);
    }

    gen_cns_by_poa(g);
}
**/

int64_t suffix_gorder_check(poa_g_t *g, integer_t *buf, uint64_t gk_0, uint64_t gk_1, int64_t update_vis, uint64_t set_flag)
{
    uint32_t *g_idx = g->srt_b.res.a, *n2gidx = g->srt_b.res2nid.a, a_n, v, init_n, k; 
    poa_arc_t *a;
    if(buf->vis.n != g->seq.n) {
        kv_resize(uint32_t, buf->vis, g->seq.n); buf->vis.n = g->seq.n;
        memset(buf->vis.a, -1, sizeof((*buf->vis.a))*buf->vis.n);
    }

    if(update_vis) {
        v = g_idx[gk_0]<<1; init_n = buf->vis.n;
        kv_push(uint32_t, buf->vis, v); 
        while (buf->vis.n > init_n) {
            v = buf->vis.a[--buf->vis.n];
            buf->vis.a[n2gidx[v>>1]] = set_flag;
            a_n = poa_arc_n(g, v); a = poa_arc_a(g, v);
            for (k = 0; k < a_n; k++) {
                if(buf->vis.a[n2gidx[a[k].v>>1]] == set_flag) continue;
                kv_push(uint32_t, buf->vis, a[k].v);
            }
        }
        assert(buf->vis.n == init_n);
    }

    if(buf->vis.a[gk_1] == set_flag) return 1;
    return 0;
}

int64_t integer_g_chain(poa_g_t *g, ma_ug_t *ug, integer_aln_t *a, int64_t a_n, integer_t *buf, ul_chain_t *res)
{
    res->v = res->s = res->e = (uint32_t)-1; res->sc = (uint64_t)-1; 
    res->q_sidx = res->q_eidx = res->t_sidx = res->t_eidx = (uint32_t)-1;
    if(a_n <= 0) return 0;
    int64_t i, k, *p, *f, tf, ti, csc, sc, max_f, max_k, vis_i, pas; integer_aln_t *li, *lk; 
    buf->vis.n = 0; vis_i = -1;
    for (i = 1; i < a_n; ++i) {//already sorted by qk
        if(((uint32_t)a[i].tn_rev_qk) <= ((uint32_t)a[i-1].tn_rev_qk)) break; ///== means there is a circle
        if(a[i].tk == a[i-1].tk) break;
        if(a[i].tk < a[i-1].tk) {
            pas = suffix_gorder_check(g, buf, a[i].tk, a[i-1].tk, vis_i==i?0:1, i); vis_i = i;
            if(pas) break;
        }
    }
    // fprintf(stderr, "[M::%s::] i::%ld, a_n::%ld\n", __func__, i, a_n);
    
    if(i >= a_n) {
        res->s = 0; res->e = a_n;
        return 1;
    }

    buf->p.n = buf->f.n = 0; 
    kv_resize(int64_t, buf->p, (uint64_t)a_n); p = buf->p.a;
    kv_resize(int64_t, buf->f, (uint64_t)a_n); f = buf->f.a;
    
    tf = ti = -1; buf->vis.n = 0; vis_i = -1;
    for (i = 0; i < a_n; ++i) {
        li = &(a[i]); csc = li->sc;
        max_f = csc; max_k = -1;
        
        for (k = i-1; k >= 0; --k) {
            lk = &(a[k]);
            ///qk of lk and li might be equal
            if(((uint32_t)lk->tn_rev_qk) >= ((uint32_t)li->tn_rev_qk)) continue;
            if(lk->tk == li->tk) continue;
            if(lk->tk > li->tk) {
                pas = suffix_gorder_check(g, buf, a[i].tk, a[k].tk, vis_i==i?0:1, i);
                // if(i == 17 && a[i].tk == 6 && k == 16 && a[k].tk == 60) {
                //     fprintf(stderr, "******i::%ld, a[i].tk::%u, k::%ld, a[k].tk::%u, vis_i::%ld, pas::%ld\n", 
                //         i, a[i].tk, k, a[k].tk, vis_i, pas);
                // }
                vis_i = i; if(pas) continue;
            }
            sc = csc + f[k];
            if(sc > max_f) {
                max_f = sc; max_k = k;
            }
        }
        f[i] = max_f; p[i] = max_k;
        if(tf < max_f) {
            tf = max_f; ti = i;
        }
        // fprintf(stderr, "[M::%s::i->%ld] qk::%u, tk::%u, max_k::%ld, max_f::%ld, ti::%ld\n", 
        // __func__, i, (uint32_t)li->tn_rev_qk, li->tk, max_k, max_f, ti);
    }

    if(ti < 0) return 0;
    for (i = ti, k = 0; i >= 0; i = p[i]) f[k++] = i;
    assert(k > 0);
    for (i = 0, k--; k >= 0; k--) {
        a[i] = a[f[k]]; i++;
    }
    res->s = 0; res->e = i;
    return 1;
}


void poa_chain_0(poa_g_t *g, ma_ug_t *ug, uc_block_t *raw, ul_str_t *str, int64_t s, int64_t e, int64_t is_rev, integer_t *buf, int64_t str_id, int64_t debug_qid, uint32_t *is_circle)
{
    (*is_circle) = 0;
    if(e <= s) return;
    uint32_t *g_idx = g->srt_b.res.a; uc_block_t *z; integer_aln_t *b; ul_chain_t rr; int64_t i, k, n;
    uint64_t *pat = (is_rev?(str->a + str->cn - e):(str->a + s)), pp; int64_t pat_n = e - s;
    // fprintf(stderr, "[M::%s::] ts::%ld, te::%ld, is_rev::%ld, g->seq.n::%u, pat_n::%lu\n", 
    // __func__, s, e, is_rev, (uint32_t)g->seq.n, pat_n);

    g->srt_b.aln.n = 0; n = g->seq.n;
    for (k = 0; k < n; k++) {///graph
        pp = (((uint64_t)(g->seq.a[g_idx[k]].nid))<<32); pp += ((uint64_t)(k)); pp += ((uint64_t)(0x80000000));
        kv_push(uint64_t, g->srt_b.aln, pp);
    }
    for (k = 0; k < pat_n; k++) {
        pp = (uint32_t)pat[poa_str_idx(k, pat_n, is_rev)]; if(is_rev) pp ^= 1; pp <<= 32; pp += ((uint64_t)(k));
        kv_push(uint64_t, g->srt_b.aln, pp);
    }

    radix_sort_srt64(g->srt_b.aln.a, g->srt_b.aln.a + g->srt_b.aln.n); n = g->srt_b.aln.n;
    for (k = 0, buf->b.n = 0; k < n; k++) {
        if(g->srt_b.aln.a[k]&((uint64_t)(0x80000000))) continue;///skip nodes in the graph
        for (i = k+1; (i < n) && ((g->srt_b.aln.a[k]>>32) == (g->srt_b.aln.a[i]>>32)); i++) {
            if((g->srt_b.aln.a[i]&((uint64_t)(0x80000000))) == 0) continue;///skip nodes in the read
            ///a[k] is read (q); a[i] is graph (t)
            kv_pushp(integer_aln_t, buf->b, &b);
            b->vq = g->srt_b.aln.a[k]>>32;
            b->tn_rev_qk = (uint32_t)g->srt_b.aln.a[k];
            b->tk = ((g->srt_b.aln.a[i]<<33)>>33);
            z = &(raw[pat[poa_str_idx(b->tn_rev_qk, pat_n, is_rev)]>>32]);
            b->sc = ug_occ_w(z->ts, z->te, &(ug->u.a[z->hid]));
        }
    }

    n = buf->b.n;
    radix_sort_integer_aln_t_srt(buf->b.a, buf->b.a + buf->b.n); ///sorted by qk
    // fprintf(stderr, "[M::%s::] # align pairs::%ld\n", __func__, n);

    // for (i = 0; i < n; i++) {///sort score
    //     b = &(buf->b.a[i]); z = &(raw[pat[poa_str_idx(b->tn_rev_qk, pat_n, is_rev)]>>32]);
    //     pp = ug_occ_w(z->ts, z->te, &(ug->u.a[z->hid]));
    //     b->tn_rev_qk += (pp<<32);
    // }
    

    i = integer_g_chain(g, ug, buf->b.a, buf->b.n, buf, &rr); assert(i);
    if(!i) return; buf->b.n = rr.e; assert(buf->b.n);
    append_aligned_integer_seq_by_aln_pair(ug, raw, g, g->seq.n, pat, pat_n, is_rev, buf->b.a, buf->b.n, str_id, (is_rev?(str->cn - e):(s)));
    (*is_circle) = update_poa_dp(g, debug_qid, ug);
}

void gen_cns_by_poa(poa_g_t *g)
{
    uint32_t n_vx = g->seq.n<<1, i, k, v, nv, w; uint64_t c, cc, n_pending = 0;
    ubuf_t *b = &(g->bb); poa_arc_t *av; uinfo_t *t;
    b->a.n = b->S.n = b->T.n = b->b.n = b->e.n = 0;
    kv_resize(uinfo_t, b->a, n_vx); b->a.n = n_vx; memset(b->a.a, 0, sizeof(*(b->a.a))*b->a.n);
    for (k = 0; k < g->seq.n; k++) {
        v = (k<<1) + 1;
        if(poa_arc_n(g, v)) continue;
        v ^= 1; kv_push(uint32_t, b->S, v); b->a.a[v].p = (uint32_t)-1;
    }
    assert(b->S.n);
    while (b->S.n > 0) {
        v = kv_pop(b->S); c = b->a.a[v].c;
        nv = poa_arc_n(g, v); av = poa_arc_a(g, v);
        for (i = 0; i < nv; ++i) {
            w = av[i].v; t = &b->a.a[w];
            kv_push(uint32_t, b->e, ((g->idx.a[v]>>32)+i)); ///push the edge
            cc = c + (uint32_t)av[i].ul;
            if (t->s == 0) {///a new node
                kv_push(uint32_t, b->b, w); // save it for revert
                t->p = v; t->s = 1; //t->d = d + l;
                t->r = poa_arc_n(g, w^1); t->c = cc; ///t->nc = c_nc; 
                ++n_pending; 
            } else {
                if(cc > t->c) {
                    t->p = v; t->c = cc; //t->s = 1; t->d = d + l; t->nc = c_nc; 
                }
            }

            if (--(t->r) == 0) {
                if(poa_arc_n(g, w) > 0) kv_push(uint32_t, b->S, w);
                --n_pending;
                // if(w == dest && n_pending == 0) goto pp_end;
            }
        }
    }
    assert(!n_pending);
    uint64_t m = 0, mi = (uint64_t)-1;
    for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
        t = &b->a.a[b->b.a[i]]; ///memset(t, 0, sizeof(*(t)));
        //b->srt.a[i].c = ((uint64_t)-1) - t->c; b->srt.a[i].i = i;
        if(m < t->c) {
            m = t->c; mi = i; ///b->b.a[i];
        }
    }

    if(mi != (uint64_t)-1) {
        g->srt_b.res.n = 0;
        for(v = b->b.a[mi]; v != (uint32_t)-1; v = b->a.a[v].p) {
            kv_push(uint32_t, g->srt_b.res, v);
        }
        m = g->srt_b.res.n>>1;
        for (i = 0; i < m; i++) {
            v = g->srt_b.res.a[i]; 
            g->srt_b.res.a[i] = g->srt_b.res.a[g->srt_b.res.n-i-1];
            g->srt_b.res.a[g->srt_b.res.n-i-1] = v;
        }
    }

    // for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
    //     t = &b->a.a[b->b.a[i]]; memset(t, 0, sizeof(*(t)));
    // }
}

void poa_cns_chain(poa_g_t *g, all_ul_t *ul_idx, ma_ug_t *ug, ul_str_t *str, ul_chain_t *idx, int64_t idx_n, int64_t qid, integer_t *buf, uint32_t *is_circle)
{
    (*is_circle) = 0;
    int64_t k, tid, is_rev; 
    reset_poa_g_t(g); 

    append_unmatch_integer_seq(g, ug, ul_idx->a[qid].bb.a, &(str[qid]), 0, str[qid].cn, 0, qid);
    clean_poa_g_t(g); (*is_circle) = 1 - topo_srt_gen(g, qid, ug);
    if((*is_circle)) return;

    for (k = 0; k < idx_n; k++) {
        tid = idx[k].v>>1; is_rev = idx[k].v&1;
        // fprintf(stderr, "\n[M::%s::] k::%ld, tid::%ld, is_rev::%ld\n", __func__, k, tid, is_rev);
        // print_integer_seq(ug, str, tid, 1);
        poa_chain_0(g, ug, ul_idx->a[tid].bb.a, &(str[tid]), idx[k].t_sidx, idx[k].t_eidx, is_rev, buf, tid, qid, is_circle);
        if((*is_circle)) return;
        // print_integer_g(g, ug, 1);
    }

    gen_cns_by_poa(g);
}

uint64_t cal_forward_dis(asg_t *g, uc_block_t *a, uint32_t s, uint32_t e)
{
    uint32_t i, v, w, nv, z; int64_t l; asg_arc_t *av;
    //fprintf(stderr, "\n[M::%s::] s::%u, e::%u\n", __func__, s, e);
    ///TODO: a[i].aidx might be < i; if s == e, then i <= e might be wrong, cannot pass the assert(li == e);
    for(i = s, l = 0, v = w = (uint32_t)-1; i != (uint32_t)-1 && i != e; i = a[i].aidx) {
        w = (((uint32_t)(a[i].hid))<<1)|((uint32_t)(a[i].rev));
        if(v != (uint32_t)-1) {
            av = asg_arc_a(g, v); nv = asg_arc_n(g, v); 
            for (z = 0; z < nv; z++) {
                if(av[z].del) continue;
                if(av[z].v == w) break;
            }
            if(z < nv) {//found
                l += (uint32_t)av[z].ul;
            } else {
                l += a[i].pdis + g->seq[v>>1].len - g->seq[w>>1].len;
            }
        }
        v = w; 
        // fprintf(stderr, "[M::%s::] i::%u, a[i].aidx::%u, a[i].qs::%u, a[i].qe::%u\n", __func__, i, a[i].aidx, a[i].qs, a[i].qe);
    }


    assert(i == e);
    w = (((uint32_t)(a[i].hid))<<1)|((uint32_t)(a[i].rev));
    if(v != (uint32_t)-1) {
        av = asg_arc_a(g, v); nv = asg_arc_n(g, v); 
        for (z = 0; z < nv; z++) {
            if(av[z].del) continue;
            if(av[z].v == w) break;
        }
        if(z < nv) {//found
            l += (uint32_t)av[z].ul;
        } else {
            l += a[i].pdis + g->seq[v>>1].len - g->seq[w>>1].len;
        }
    }
    v = w; 


    if(l < 0) l = 0;
    return l;
}

uint64_t cal_integer_match_dis(ma_ug_t *ug, uc_block_t *a, int64_t k_0, int64_t k_1, int64_t is_rev, uint32_t *is_g_connect)
{
    uint32_t i, k;
    assert((is_rev && k_1 < k_0) || ((!is_rev) && k_1 > k_0)); (*is_g_connect) = 0;
    if(is_rev) {
        i = k_0; k = k_1;
    } else {
        i = k_1; k = k_0;
    }
    // fprintf(stderr, "k_0::%ld, k_1::%ld\n", k_0, k_1);
    uint32_t li, lk, pk, bi = i; int64_t l;
    for (li = i, l = 0; i != (uint32_t)-1 && i >= k; i = a[i].pidx) {
        li = i; if(a[i].pidx != (uint32_t)-1 && i > k) l += a[i].pdis;
    }
    // fprintf(stderr, "+bi::%u, li::%u, i::%u, l::%ld\n", bi, li, i, l);
    if(is_rev) l = cal_forward_dis(ug->g, a, li, bi);
    // fprintf(stderr, "++bi::%u, li::%u, i::%u, l::%ld\n", bi, li, i, l);
    if(li == k) {///direct path
        (*is_g_connect) = 1;    
        return l;
    }
    i = li;
    assert(i > k); pk = k;
    for (lk = k; k != (uint32_t)-1 && k <= i; k = a[k].aidx) lk = k;
    // fprintf(stderr, "-pk::%u, lk::%u, k::%u\n", pk, lk, k);
    if(!is_rev) {
        for (k = lk; k != (uint32_t)-1 && k != pk; k = a[k].pidx) l += a[k].pdis;
    } else {
        l += cal_forward_dis(ug->g, a, pk, lk);
    }
    // fprintf(stderr, "--pk::%u, lk::%u, k::%u, l::%ld\n", pk, lk, k, l);
    if(l < 0) l = 0;
    k = lk;
    assert(i > k);
    return l + normlize_gdis(ug, &(a[i]), &(a[k]), is_rev);
}

uint64_t cal_integer_most_dis(uint32_t qid, uint64_t *a, uint64_t a_n, double cluster_rate)
{
    if(a_n <= 0) return (uint64_t)-1;
    // fprintf(stderr, "\n[M::%s::] a_n::%lu\n", __func__, a_n);
    uint64_t k, l, i, m, r_an = a_n, max_m, max_i, cc, cd, nd; int64_t z;
    for (k = 1, l = m = max_m = 0, max_i = (uint64_t)-1; k <= a_n; k++) {
        a[k-1] <<= 1; a[k-1] >>= 1;
        // fprintf(stderr, "[k->%lu] d::%lu, rev::%lu\n", k-1, a[k-1]>>1, a[k-1]&1);
        if((k == a_n) || (((a[k]&((uint64_t)(0x7fffffffffffffff)))>>1) 
                                        != ((a[l]&((uint64_t)(0x7fffffffffffffff)))>>1))) {
            for (i = l; i < k; i++) {
                if(!(a[i]&1)) break;
            }
            a[m] = a[l]; a[m] >>= 1; a[m] <<= 1; 
            ///if all integer sequence are mapped reversely 
            if(i >= k) a[m] += 1;
            a[m] |= ((uint64_t)(k-l))<<32;
            if((k-l) > max_m) {
                max_m = k - l; max_i = m;
            } else if((k-l) == max_m && i < k) {///i < k means this distance is supported by forward sequences
                max_m = k - l; max_i = m;
            }
            l = k; m++;
        }
    }
    a_n = m; assert(a_n > 0);
    // fprintf(stderr, "[M::%s::] m::%lu\n", __func__, m);
    if(max_m > 0 && max_m > (r_an>>1)) {
        // fprintf(stderr, "[M::%s::] dis::%u\n", __func__, ((uint32_t)a[max_i])>>1);
        return ((uint32_t)a[max_i])>>1;
    }

    for (k = 0, max_m = 0, max_i = (uint64_t)-1; k < a_n; k++) {
        cd = (((uint32_t)a[k])>>1); z = k;
        // fprintf(stderr, "[mk->%lu] cd::%lu\n", k, cd);
        for (cc = cd, z--; z >= 0; z--) {
            nd = (((uint32_t)a[z])>>1); assert(nd < cd);
            if(((cd-nd) > (nd*cluster_rate)) && ((cd-nd) > 512)) break;
            cc += (a[z]>>32);
        }
        for (i = k+1; i < a_n; i++) {
            nd = (((uint32_t)a[i])>>1); 
            // if(!(nd > cd)) {
            //     fprintf(stderr, "[M::%s::] qid::%u, a_n::%lu, i::%lu, k::%lu, cd::%lu, nd::%lu\n", 
            //     __func__, qid, a_n, i, k, cd, nd);
            // }
            assert(nd > cd);
            if(((nd-cd) > (nd*cluster_rate)) && ((nd-cd) > 512)) break;
            cc += (a[i]>>32);
        }

        if(cc > max_m) {
            max_m = cc; max_i = k;
        } else if(cc == max_m && (!(a[k]&1))) {///means this distance is supported by forward sequences
            max_m = cc; max_i = k;
        }
    }

    k = max_i;
    cd = (((uint32_t)a[k])>>1); z = k;
    for (max_m = cd, max_i = k, z--; z >= 0; z--) {
        nd = (((uint32_t)a[z])>>1); assert(nd < cd);
        if(((cd-nd) > (nd*cluster_rate)) && ((cd-nd) > 512)) break;
        cc = (a[z]>>32);
        if(cc > max_m) {
            max_m = cc; max_i = z;
        } else if(cc == max_m && (!(a[z]&1))) {///means this distance is supported by forward sequences
            max_m = cc; max_i = z;
        }
    }
    for (i = k+1; i < a_n; i++) {
        nd = (((uint32_t)a[i])>>1); assert(nd > cd);
        if(((nd-cd) > (nd*cluster_rate)) && ((nd-cd) > 512)) break;
        cc = (a[i]>>32);
        if(cc > max_m) {
            max_m = cc; max_i = i;
        } else if(cc == max_m && (!(a[i]&1))) {///means this distance is supported by forward sequences
            max_m = cc; max_i = i;
        }
    }

    // fprintf(stderr, "[M::%s::] dis::%u\n", __func__, ((uint32_t)a[max_i])>>1);
    return ((uint32_t)a[max_i])>>1;
}

uint32_t poa_g_arc_w(poa_g_t *pg, uint32_t v, uint32_t w)
{
    poa_arc_t *av; uint32_t an, k;
    an = poa_arc_n(pg, v); av = poa_arc_a(pg, v);
    for (k = 0; k < an; k++) {
        if(av[k].v == w) return (uint32_t)av[k].ul;
    }
    return 0;
}

void dump_cns_res(poa_g_t *pg, ma_ug_t *ug, uint32_t *cns_seq, uint32_t cns_occ, all_ul_t *ul_idx, ul_str_t *str, uint32_t qid, integer_t *buf,
uint64_t *arc_idx, uint64_t arc_idx_n)
{   
    uint64_t t, dd, k, i; uint32_t e_s, e_e, is_rev, is_g_connect, con_occ; emap_t *g_arc;
    if(cns_occ > 0) {
        t = qid; t <<= 32; t |= ((uint64_t)(0xffffffff));
        kv_push(uint64_t, buf->res_dump, t);

        t = pg->seq.a[cns_seq[0]>>1].nid; 
        // t = pg->seq.a[cns_seq[0]>>1].nid; t |= ((uint64_t)(0xffffffff00000000));
        t |= ((uint64_t)(ug->g->seq[pg->seq.a[cns_seq[0]>>1].nid>>1].len))<<32;
        kv_push(uint64_t, buf->res_dump, t);
        buf->n_correct++;
    }
    for (k = 0; k < arc_idx_n; k++) {
        // csn_v = pg->seq.a[cns_seq[k]>>1].nid; cns_w = pg->seq.a[cns_seq[k+1]>>1].nid;
        e_s = arc_idx[k]>>32; e_e = (uint32_t)arc_idx[k]; 
        // fprintf(stderr, "[M::%s::] k::%lu, arc_idx_n::%lu, e_s::%u, e_e::%u\n", __func__, k, arc_idx_n, e_s, e_e);
        assert(poa_g_arc_w(pg, cns_seq[k], cns_seq[k+1]) == (e_e - e_s)); assert(e_e > e_s); 
        buf->o.n = 0; kv_resize(uint64_t, buf->o, e_e - e_s); con_occ = 0;
        for (i = e_s; i < e_e; i++) {
            g_arc = &(pg->e_idx.a[i]);
            assert((g_arc->pge>>32) == (cns_seq[k]>>1) && ((uint32_t)g_arc->pge) == (cns_seq[k+1]>>1));
            is_rev = ((g_arc->ule>>32) > ((uint32_t)g_arc->ule)?1:0);
            assert((pg->seq.a[g_arc->pge>>32].nid^(is_rev?1:0)) == ((uint32_t)str[g_arc->ulid].a[g_arc->ule>>32]));
            assert((pg->seq.a[(uint32_t)g_arc->pge].nid^(is_rev?1:0)) == ((uint32_t)str[g_arc->ulid].a[(uint32_t)g_arc->ule]));
            // fprintf(stderr, "+++i::%lu, target_ulid::%u, str_sidx::%u, str_eidx::%u, str_cn::%u, ul_idx->a[g_arc->ulid].bb.n::%u\n", 
            // i - e_s, g_arc->ulid, (uint32_t)(g_arc->ule>>32), (uint32_t)g_arc->ule, str[g_arc->ulid].cn, (uint32_t)ul_idx->a[g_arc->ulid].bb.n);
            dd = cal_integer_match_dis(ug, ul_idx->a[g_arc->ulid].bb.a, str[g_arc->ulid].a[g_arc->ule>>32]>>32, 
                                                                    str[g_arc->ulid].a[(uint32_t)g_arc->ule]>>32, is_rev, &is_g_connect);
            // fprintf(stderr, "---i::%lu, dd::%lu, is_g_connect::%u\n", i - e_s, dd, is_g_connect);
            dd <<= 1; if(is_rev) dd += 1; 
            if(is_g_connect) {
                con_occ++; buf->o.a[buf->o.n] = dd;
            } else {
                buf->o.a[buf->o.n] = dd; buf->o.a[buf->o.n] |= ((uint64_t)(0x8000000000000000));
            } 
            buf->o.n++;
        }

        radix_sort_srt64(buf->o.a, buf->o.a + buf->o.n);
        if(con_occ > 0) buf->o.n = con_occ;
        dd = cal_integer_most_dis(qid, buf->o.a, buf->o.n, 0.08);


        t = pg->seq.a[cns_seq[k+1]>>1].nid; t |= ((uint64_t)(dd<<32)); 
        if(con_occ > 0) t |= ((uint64_t)(0x8000000000000000));
        kv_push(uint64_t, buf->res_dump, t);
    }

}

void update_raw_integer_seq(poa_g_t *pg, ma_ug_t *ug, uint32_t *cns_seq, uint32_t cns_occ, all_ul_t *ul_idx, ul_str_t *str, uint32_t qid, integer_t *buf, ul_chain_t *idx_a, uint64_t idx_n)
{
    if(cns_occ <= 0) return;
    ///no need to check edges if there is just one node -> no dege
    if(cns_occ > 1) radix_sort_emap_t_srt(pg->e_idx.a, pg->e_idx.a + pg->e_idx.n); 
    kv_resize(uint64_t, buf->u, cns_occ); buf->u.n = cns_occ - 1; memset(buf->u.a, 0, sizeof(*(buf->u.a)*buf->u.n));
    uint64_t *arc_idx = buf->u.a, arc_idx_n = buf->u.n, x; uint64_t k, l, i, n = pg->e_idx.n, fe = cns_occ - 1;
    for (k = 1, l = 0; k <= n && fe > 0; k++) {
        if(k == n || pg->e_idx.a[k].pge != pg->e_idx.a[l].pge) {
            if(k > l) {
                for (i = 0; i < arc_idx_n; i++) {
                    x = (cns_seq[i]>>1); x <<= 32; x += (cns_seq[i+1]>>1);
                    if(x == pg->e_idx.a[l].pge) {
                        arc_idx[i] = l; arc_idx[i] <<= 32; arc_idx[i] += k; fe--;
                        break;
                    }
                }
            }
            l = k;
        }
    }
    assert(fe == 0);
    // for (k = 0; k < pg->e_idx.n; k++) {
    //     assert((pg->seq.a[pg->e_idx.a[k].pge>>32].nid>>1) == 
    //                             (((uint32_t)str[pg->e_idx.a[k].ulid].a[pg->e_idx.a[k].ule>>32])>>1));
    //     assert((pg->seq.a[(uint32_t)pg->e_idx.a[k].pge].nid>>1) == 
    //                             (((uint32_t)str[pg->e_idx.a[k].ulid].a[(uint32_t)pg->e_idx.a[k].ule])>>1));
    //     assert(((pg->seq.a[pg->e_idx.a[k].pge>>32].nid&1)^(((uint32_t)str[pg->e_idx.a[k].ulid].a[pg->e_idx.a[k].ule>>32])&1)) ==
    //         ((pg->seq.a[(uint32_t)pg->e_idx.a[k].pge].nid&1)^(((uint32_t)str[pg->e_idx.a[k].ulid].a[(uint32_t)pg->e_idx.a[k].ule])&1)));
    // }

    uint32_t e_s, e_e; emap_t *g_arc; uint64_t l_clip = 0, r_clip = 0; int64_t e_occ;
    if(arc_idx_n > 0) {
        ///clip unreliable left/right end
        for (l_clip = 0; l_clip < arc_idx_n; l_clip++) {
            k = l_clip;
            e_s = arc_idx[k]>>32; e_e = (uint32_t)arc_idx[k]; 
            assert(poa_g_arc_w(pg, cns_seq[k], cns_seq[k+1]) == (e_e - e_s)); 
            e_occ = ((int64_t)e_e) - ((int64_t)e_s);
            if(e_occ > 1) break;///more than one read supporting this edge
            for (i = e_s; i < e_e; i++) {
                g_arc = &(pg->e_idx.a[i]);
                if(g_arc->ulid == qid) break;
            }
            if(i < e_e) break;///the query read itself supports this edge
        }
        
        for (r_clip = 0; r_clip < arc_idx_n; r_clip++) {
            k = arc_idx_n - r_clip - 1;
            e_s = arc_idx[k]>>32; e_e = (uint32_t)arc_idx[k]; 
            assert(poa_g_arc_w(pg, cns_seq[k], cns_seq[k+1]) == (e_e - e_s)); 
            e_occ = ((int64_t)e_e) - ((int64_t)e_s);
            if(e_occ > 1) break;///more than one read supporting this edge
            for (i = e_s; i < e_e; i++) {
                g_arc = &(pg->e_idx.a[i]);
                if(g_arc->ulid == qid) break;
            }
            if(i < e_e) break;///the query read itself supports this edge
        }
    }
    if(l_clip+r_clip >= cns_occ) return;
    dump_cns_res(pg, ug, cns_seq+l_clip, cns_occ-l_clip-r_clip, ul_idx, str, qid, buf, arc_idx+l_clip, arc_idx_n-l_clip-r_clip);
}

void integer_candidate(ul_resolve_t *uidx, integer_t *buf, uint32_t qid, uint32_t is_hom)
{
    uint64_t k, z, m_het, m_het_occ, ref_occ, b_n, m; uint32_t vk, vz, is_circle = 0; integer_aln_t *p; ul_chain_t sc;
    ul_str_idx_t *str_idx = &(uidx->pstr); ma_ug_t *ug = uidx->l1_ug; 
    uint64_t *hid_a, hid_n; uc_block_t *xi;
    ul_str_t *str = &(str_idx->str.a[qid]);
    if(str->cn < 2) return; ///directly filter out too short UL
    kv_resize(uint64_t, buf->u, str->cn); buf->u.n = str->cn;
    for (k = m_het = m_het_occ = ref_occ = 0; k < str->cn; k++) {
        xi = &(uidx->idx->a[qid].bb.a[str->a[k]>>32]);
        assert(((xi->hid<<1)+xi->rev)==((uint32_t)str->a[k]));
        buf->u.a[k] = ug_occ_w(xi->ts, xi->te, &(ug->u.a[xi->hid]));
        assert(buf->u.a[k] > 0);
        if((!is_hom) && (!IF_HOM((((uint32_t)str->a[k])>>1), (*uidx->bub)))) {
            m_het++; m_het_occ += buf->u.a[k];///m_het_occ: how many het HiFi reads
        }
        ref_occ += buf->u.a[k];
        // fprintf(stderr, "[M::%s::k->%lu] buf->u.a[k]->%lu, ts->%u, te->%u, pchain->%u\n", __func__, k, buf->u.a[k], xi->ts, xi->te, xi->pchain);
    }
    // if((!is_hom) && (m_het < 2) && (m_het > 0)) return;///if all matched unitigs are hom, is ok
    if(m_het == 0 || m_het_occ == 0) is_hom = 1;
    // print_ul_alignment(ug, &UL_INF, 27512, "inner-0");
    for (k = 0, buf->b.n = 0; k < str->cn; k++) {
        vk = (uint32_t)str->a[k];
        hid_a = str_idx->occ.a + str_idx->idx.a[vk>>1];
        hid_n = str_idx->idx.a[(vk>>1)+1] - str_idx->idx.a[vk>>1];
        for (z = 0; z < hid_n; z++) {
            if((hid_a[z]>>32) == qid) continue;
            if(str_idx->str.a[hid_a[z]>>32].cn < 2) continue;
            vz = (uint32_t)(str_idx->str.a[hid_a[z]>>32].a[(uint32_t)hid_a[z]]);
            assert((vk>>1) == (vz>>1));
            kv_pushp(integer_aln_t, buf->b, &p);
            p->vq = vk; p->tk = (uint32_t)hid_a[z]; 
            if((vk^vz)&1) p->tk = str_idx->str.a[hid_a[z]>>32].cn - p->tk - 1;///rev
            p->tn_rev_qk = (hid_a[z]>>32); p->tn_rev_qk <<= 1; p->tn_rev_qk |= ((vk^vz)&1);
            p->tn_rev_qk <<= 32; p->tn_rev_qk += k; 
            ///set score of this pair
            p->sc = buf->u.a[k];
            xi = &(uidx->idx->a[hid_a[z]>>32].bb.a[(str_idx->str.a[hid_a[z]>>32].a[(uint32_t)hid_a[z]])>>32]);
            assert(((xi->hid<<1)+xi->rev)==vz); m = ug_occ_w(xi->ts, xi->te, &(ug->u.a[xi->hid]));
            if(p->sc > m) p->sc = m;
        }
    }
    // print_ul_alignment(ug, &UL_INF, 27512, "inner-1");
    radix_sort_integer_aln_t_srt(buf->b.a, buf->b.a + buf->b.n); 
    b_n = buf->b.n; buf->sc.n = 0;
    for (k = 1, z = 0; k <= b_n; k++) {
        if(k == b_n || (buf->b.a[z].tn_rev_qk>>32) != (buf->b.a[k].tn_rev_qk>>32)) {
            ///get the chain for <qn, tn>
            if(integer_chain(qid, buf->b.a + z, k - z, z, buf, ug, str_idx, uidx->idx, &sc) && sc.v != (uint32_t)-1) {
                if((buf->sc.n > 0) && ((buf->sc.a[buf->sc.n-1].v>>1) == (sc.v>>1))) {
                    if(buf->sc.a[buf->sc.n-1].sc < sc.sc) {
                        buf->sc.a[buf->sc.n-1] = sc;
                    }
                } else {
                    kv_push(ul_chain_t, buf->sc, sc);
                }
            }
            z = k;
        }
    }
    // print_ul_alignment(ug, &UL_INF, 27512, "inner-2");
    uint64_t *o, o_n, cns_het, cns_het_occ, ref_cns_occ, corrected = 0;
    o_n = integer_chain_dp(uidx->bub, buf, str_idx->str.a, buf->b.a, buf->sc.a, buf->sc.n, qid, is_hom, 2, &corrected);
    assert(o_n <= str->cn); 
    // if(qid == 2062 || qid == 2093) fprintf(stderr,"[M::%s::] o_n::%lu, str->cn::%u, corrected::%lu\n", __func__, o_n, str->cn, corrected);
    if(o_n <= 0) return; 
    if(corrected) return;
    // print_ul_alignment(ug, &UL_INF, 27512, "inner-3");
    
    for (k = cns_het = cns_het_occ = ref_cns_occ = 0, o = buf->o.a; k < o_n; k++) {
        if((!is_hom) && (!IF_HOM((((uint32_t)str->a[o[k]])>>1), (*uidx->bub)))) {
            cns_het++; cns_het_occ += buf->u.a[o[k]];
        }
        ref_cns_occ += buf->u.a[o[k]];     
    }
    buf->o.n = o_n;
    ///1. if the ref read only has hom unitigs, is fine
    ///2. otherwise need to have consenus het untigs
    if(!is_hom) {
        if((cns_het <= 0) || (cns_het_occ <= 0) || (cns_het_occ <= (m_het_occ*0.25))) return;
    } else {
        if(ref_cns_occ <= (ref_occ*0.25)) return;
    }
    // print_ul_alignment(ug, &UL_INF, 27512, "inner-4");
    // fprintf(stderr, "\n");
    // print_integer_seq(ug, str_idx->str.a, qid, 1);
    // print_aln_seq(ug, uidx->idx, qid, 1);
    // if(qid == 2062 || qid == 2093) print_cns_seq(ug, str, o, o_n);
    // print_ul_alignment(ug, &UL_INF, 27512, "inner-5");

    for (k = m = 0; k < buf->sc.n; k++) {
        // fprintf(stderr, "[M::%s::k->%lu] m::%lu\n", __func__, k, m);
        // if(k == 72) {
        //     print_integer_ovlps(ug, str_idx->str.a, buf->b.a, buf->b.n, buf->sc.a+k, 1, qid, o_n);
        // }
        if(refine_integer_ovlps(uidx->idx, uidx->bub, ug, str_idx->str.a, buf->b.a, 
                                                            &(buf->sc.a[k]), qid, buf, o, o_n)) {                                                
            buf->sc.a[m++] = buf->sc.a[k];
        } 
        // else {
        //     print_integer_ovlps(ug, str_idx->str.a, buf->b.a, buf->b.n, buf->sc.a+k, 1, qid, o_n);
        //     fprintf(stderr, "[M::%s::] idx->q_sidx::%u, idx->q_eidx::%u, idx->t_sidx::%u, idx->t_eidx::%u\n******************************************************\n", 
        //         __func__, buf->sc.a[k].q_sidx, buf->sc.a[k].q_eidx, buf->sc.a[k].t_sidx, buf->sc.a[k].t_eidx);
        // }
    }
    buf->sc.n = m;
    if(m <= 0) return;
    // print_ul_alignment(ug, &UL_INF, 27512, "inner-6");
    // if(m != str->cn) print_integer_ovlps(uidx->l1_ug, str_idx->str.a, buf->b.a, buf->b.n, buf->sc.a, buf->sc.n, qid, m);
    poa_cns_chain(&(buf->pg), uidx->idx, ug, str_idx->str.a, buf->sc.a, buf->sc.n, qid, buf, &is_circle);
    if(is_circle) {
        buf->n_circle++;
        return;
    }
    // print_ul_alignment(ug, &UL_INF, 27512, "inner-7");
    // print_res_seq(&(buf->pg), ug, buf->pg.srt_b.res.a, buf->pg.srt_b.res.n);
    // radix_sort_ul_chain_t_srt(buf->sc.a, buf->sc.a + buf->sc.n);
    
    // integer_phase(str_idx->str.a, buf, buf->sc.a, buf->sc.n, buf->b.a, qid);

    // radix_sort_ul_chain_t_srt(buf->sc.a, buf->sc.a + buf->sc.n);
    o = NULL; o_n = 0;
    update_raw_integer_seq(&(buf->pg), ug, buf->pg.srt_b.res.a, buf->pg.srt_b.res.n, uidx->idx, str_idx->str.a, qid, buf, buf->sc.a, buf->sc.n);
}

ul2ul_item_t *get_ul_ovlp(ul2ul_idx_t *z, uint64_t id, uint64_t is_ul)
{
    uint64_t x = (is_ul?(id):(id+z->uln));
    if(z->item_idx[x] == (uint32_t)-1) return NULL;
    return &(z->a[z->item_idx[x]]);
}

ul2ul_t *get_ul_o(ul2ul_idx_t *idx, uint32_t qid, uint32_t is_q_ul, uint32_t tid, uint32_t is_t_ul)
{
    ul2ul_item_t *z = get_ul_ovlp(idx, qid, is_q_ul); uint64_t k;
    if(!is_t_ul) tid += idx->uln;
    for (k = 0; k < z->n; k++) {
        if(z->a[k].hid == tid) return &(z->a[k]);
    }
    return NULL;
}

void integer_gen_ovlp(ul_resolve_t *uidx, integer_t *buf, uint32_t qid, ul2ul_item_t *o, uint32_t ug_offset)
{
    ul_str_idx_t *str_idx = &(uidx->pstr); ma_ug_t *ug = uidx->l1_ug; ul2ul_t res;
    ul_str_t *str = &(str_idx->str.a[qid]); integer_aln_t *p; ul_chain_t sc;
    uint64_t k, z, *hid_a, hid_n, sck, b_n, m/**, ovn = 0**/; uint32_t vk, vz; uc_block_t *xi;
    o->n = o->cn = 0; o->id = qid; o->is_del = 0; ///o->is_consist = 0;
    // if(qid == 95 || qid == 36) print_integer_seq(ug, str_idx->str.a, qid, 1);
    if(str->cn < 2) return;
    for (k = 0, buf->b.n = 0; k < str->cn; k++) {
        xi = &(uidx->idx->a[qid].bb.a[str->a[k]>>32]);
        assert(((xi->hid<<1)+xi->rev)==((uint32_t)str->a[k]));
        sck = ug_occ_w(xi->ts, xi->te, &(ug->u.a[xi->hid]));

        vk = (uint32_t)str->a[k];
        hid_a = str_idx->occ.a + str_idx->idx.a[vk>>1];
        hid_n = str_idx->idx.a[(vk>>1)+1] - str_idx->idx.a[vk>>1];
        for (z = 0; z < hid_n; z++) {
            // if(qid == 95 || qid == 36) {
            //     fprintf(stderr,"[M::%s::qid->%u::tid->%lu] k->%lu, tid_occ->%u\n", 
            //                 __func__, qid, (hid_a[z]>>32), k, str_idx->str.a[hid_a[z]>>32].cn);
            // }
            if((hid_a[z]>>32) == qid) continue;
            if(str_idx->str.a[hid_a[z]>>32].cn < 2) continue;
            vz = (uint32_t)(str_idx->str.a[hid_a[z]>>32].a[(uint32_t)hid_a[z]]);
            assert((vk>>1) == (vz>>1));
            kv_pushp(integer_aln_t, buf->b, &p);
            p->vq = vk; p->tk = (uint32_t)hid_a[z]; 
            if((vk^vz)&1) p->tk = str_idx->str.a[hid_a[z]>>32].cn - p->tk - 1;///rev
            p->tn_rev_qk = (hid_a[z]>>32); p->tn_rev_qk <<= 1; p->tn_rev_qk |= ((vk^vz)&1);
            p->tn_rev_qk <<= 32; p->tn_rev_qk += k; 
            ///set score of this pair
            p->sc = sck;
            xi = &(uidx->idx->a[hid_a[z]>>32].bb.a[(str_idx->str.a[hid_a[z]>>32].a[(uint32_t)hid_a[z]])>>32]);
            assert(((xi->hid<<1)+xi->rev)==vz); m = ug_occ_w(xi->ts, xi->te, &(ug->u.a[xi->hid]));
            if(p->sc > m) p->sc = m;
        }
    }

    radix_sort_integer_aln_t_srt(buf->b.a, buf->b.a + buf->b.n); 

    // b_n = buf->b.n;
    // for (k = 1, z = ovn = 0; k <= b_n; k++) {
    //     if(k == b_n || (buf->b.a[z].tn_rev_qk>>33) != (buf->b.a[k].tn_rev_qk>>33)) {
    //         ovn++; z = k;
    //     }
    // }

    b_n = buf->b.n; buf->sc.n = 0;
    for (k = 1, z = 0; k <= b_n; k++) {
        if(k == b_n || (buf->b.a[z].tn_rev_qk>>32) != (buf->b.a[k].tn_rev_qk>>32)) {
            // if(qid == 95 || qid == 36) {
            //     fprintf(stderr,"[M::%s::qid->%u::tid->%lu]\n", __func__, qid, buf->b.a[z].tn_rev_qk>>33);
            // }
            if(integer_chain(qid, buf->b.a + z, k - z, z, buf, ug, str_idx, uidx->idx, &sc) && sc.v != (uint32_t)-1) {
                if((buf->sc.n > 0) && ((buf->sc.a[buf->sc.n-1].v>>1) == (sc.v>>1))) {
                    if(buf->sc.a[buf->sc.n-1].sc < sc.sc) {
                        buf->sc.a[buf->sc.n-1] = sc;
                    }
                } else {
                    kv_push(ul_chain_t, buf->sc, sc);
                }
            }
            z = k;
        }
    }

    for (k = o->n = 0; k < buf->sc.n; k++) {
        if(update_exact_ul_ovlps(uidx->idx, ug, str_idx->str.a, buf->b.a, &(buf->sc.a[k]), qid, buf, &res)) {                                                
            kv_push(ul2ul_t, *o, res);
        }
    }
    // if(str->cn > 0) {
    //     for (k = 0, buf->b.n = 0; k < str->cn; k++) {
    //         xi = &(uidx->idx->a[qid].bb.a[str->a[k]>>32]);

    //         // if(xi->ts)
    //     }
    //     // uint32_t v, nv; asg_arc_t *av; asg_t *g = ug->g; ul2ul_t *r;
    //     // v = ((uint32_t)str->a[0])^1;
    //     // nv = asg_arc_n(g, v); av = asg_arc_a(g, v);
    //     // for (k = 0; k < nv; k++) {
    //     //     if(av[k].del) continue;
    //     //     kv_pushp(ul2ul_t, *o, &r);
    //     //     r->hid = (av[k].v>>1) + ug_offset; 
    //     //     r->is_rev = (av[k].v^v)&1;
    //     //     r->qs = 0; r->qe = av[k].ol;

    //     // }
    //     // v = ((uint32_t)str->a[str->cn-1]);
    //     // nv = asg_arc_n(g, v); av = asg_arc_a(g, v);
    // }
    o->cn = o->n; 
    // if(ovn == o->cn) o->is_consist = 1;
    if(o->n <= 0) return;
}

#define ulg_id(ul2, i) (((i)>=(ul2).uln)?((i)-(ul2).uln):(i))
#define ulg_type(ul2, i) (((i)>=(ul2).uln)?(0):(1))
#define ulg_len(uidx, i) (ulg_type((uidx).uovl,(i))?((uidx).idx->a[ulg_id((uidx).uovl,(i))].rlen):((uidx).l1_ug->u.a[ulg_id((uidx).uovl,(i))].len))
#define ulg_occ(uidx, i) (ulg_type((uidx).uovl,(i))?((uidx).pstr.str.a[ulg_id((uidx).uovl,(i))].cn):(1))

void integer_normalize_ovlp(ul_resolve_t *uidx, uint32_t qid, ul2ul_item_t *o, ul2ul_idx_t *ul2)
{
    ul2ul_t *z; uint64_t k; ul2ul_item_t *t;
    for (k = 0; k < o->n; k++) {
        assert(ul2->item_idx[o->a[k].hid] != (uint32_t)-1);
        z = get_ul_o(ul2, ulg_id(*ul2, o->a[k].hid), ulg_type(*ul2, o->a[k].hid), 
                                                            ulg_id(*ul2, qid), ulg_type(*ul2, qid));
        // assert(z && z->hid == qid);
        // if(o->a[k].is_rev == z->is_rev && (*z).qs == o->a[k].ts && (*z).qe == o->a[k].te && 
        //     (*z).ts == o->a[k].qs && (*z).te == o->a[k].qe) {
        //     fprintf(stderr, "good::[M::%s] fn::%u(%c), fqs::%u, fqe::%u, fts::%u, fte::%u, rn::%u(%c), rqs::%u, rqe::%u, rts::%u, rte::%u\n", __func__, 
        //         o->a[k].hid, "+-"[o->a[k].is_rev], o->a[k].qs, o->a[k].qe, o->a[k].ts, o->a[k].te, 
        //         z->hid, "+-"[z->is_rev], z->qs, z->qe, z->ts, z->te);
        // } else {
        //     fprintf(stderr, "bad::[M::%s] fn::%u(%c), fqs::%u, fqe::%u, fts::%u, fte::%u, rn::%u(%c), rqs::%u, rqe::%u, rts::%u, rte::%u\n", __func__, 
        //         o->a[k].hid, "+-"[o->a[k].is_rev], o->a[k].qs, o->a[k].qe, o->a[k].ts, o->a[k].te, 
        //         z->hid, "+-"[z->is_rev], z->qs, z->qe, z->ts, z->te);
        // }
        // continue;
        if(z && o->a[k].hid < qid) continue;
        if(!z) {
            t = get_ul_ovlp(ul2, ulg_id(*ul2, o->a[k].hid), ulg_type(*ul2, o->a[k].hid));
            assert(t);
            kv_pushp(ul2ul_t, *t, &z); 
            (*z).hid = qid; 
            (*z).qs = o->a[k].ts; (*z).qe = o->a[k].te; 
            (*z).ts = o->a[k].qs; (*z).te = o->a[k].qe;
            (*z).qs_k = o->a[k].ts_k; (*z).qe_k = o->a[k].te_k; 
            (*z).ts_k = o->a[k].qs_k; (*z).te_k = o->a[k].qe_k;
            (*z).is_rev = o->a[k].is_rev; (*z).is_del = o->a[k].is_del;
            (*z).is_ct = o->a[k].is_ct;
            if(!((*z).is_del)) {
                t->cn++;
            }
        } else {
            if((o->a[k].qe - o->a[k].qs) >= (z->te - z->ts)) {
                (*z).qs = o->a[k].ts; (*z).qe = o->a[k].te;
                (*z).ts = o->a[k].qs; (*z).te = o->a[k].qe;
                (*z).qs_k = o->a[k].ts_k; (*z).qe_k = o->a[k].te_k; 
                (*z).ts_k = o->a[k].qs_k; (*z).te_k = o->a[k].qe_k;
                (*z).is_rev = o->a[k].is_rev; (*z).is_del = o->a[k].is_del;
                (*z).is_ct = o->a[k].is_ct;
            } else {
                o->a[k].qs = (*z).ts; o->a[k].qe = (*z).te;
                o->a[k].ts = (*z).qs; o->a[k].te = (*z).qe;
                o->a[k].qs_k = (*z).ts_k; o->a[k].qe_k = (*z).te_k;
                o->a[k].ts_k = (*z).qs_k; o->a[k].te_k = (*z).qe_k;
                o->a[k].is_rev = (*z).is_rev; o->a[k].is_del = (*z).is_del;
                o->a[k].is_ct = (*z).is_ct;
            }
        }
    }
}


void integer_normalize_ovlp_purge(ul_resolve_t *uidx, uint32_t qid, ul2ul_item_t *o, ul2ul_idx_t *ul2)
{
    if(o->is_del) return;
    ul2ul_t *z; uint64_t k, is_del;
    for (k = is_del = 0; k < o->n; k++) {
        if(o->a[k].is_del) continue;
        z = get_ul_o(ul2, o->a[k].hid, 1, qid, 1);
        if((!z) || (z->is_del)) {
            o->a[k].is_del = 1; is_del++;
        } 
    }

    if(is_del) {
        for (k = o->cn = 0; k < o->n; k++) {
            if(o->a[k].is_del) o->a[k].hid |= ((uint32_t)(0x80000000));
            else o->cn++;
        }

        radix_sort_ul2ul_srt(o->a, o->a + o->n);

        for (k = 0; k < o->n; k++) {
            if(o->a[k].hid&((uint32_t)(0x80000000))) {
                o->a[k].hid -= ((uint32_t)(0x80000000));
            }
        }
    }

    if(o->cn == 0) o->is_del = 1;
}

static inline int64_t integer_hit2arc_idx_contain(const ul2ul_t *z, int64_t qn, int64_t tn, uint64_t *dir)
{
    int64_t tn5, tn3, extn5, extn3, qsn = z->qs_k, qen = qn - ((int64_t)z->qe_k);
    if (z->is_rev) tn5 = tn - z->te_k, tn3 = z->ts_k;
    else tn5 = z->ts_k, tn3 = tn - z->te_k;
    (*dir) = (uint64_t)-1;

    extn5 = ((qsn<tn5)?qsn:tn5);
    extn3 = ((qen<tn3)?qen:tn3);
    if(extn5 > 0 || extn3 > 0) return MA_HT_INT;///overhang
    if (qsn < tn5 && qen < tn3) { // query contained in target
        return MA_HT_QCONT; 
    } else if (qsn > tn5 && qen > tn3) { // target contained in query
        return MA_HT_TCONT;  
    } else if(qsn == tn5 && qen == tn3) {
        (*dir) = (uint64_t)-1;
    } else if (qsn > tn5) { ///query-to-target overlap
        (*dir) = 0;
    } else if(qsn < tn5) { ///target-to-query overlaps
        (*dir) = 1;
    } else if(qen > tn3) { ///target-to-query overlaps
        (*dir) = 1;
    } else if(qen < tn3) { ///query-to-target overlap
        (*dir) = 0;
    }
    return MA_HT_DOUBLE;
}

static inline int integer_hit2arc(const ul2ul_t *z, int64_t ql, int64_t tl, int64_t qocc, int64_t tocc, uint64_t qid, uint64_t tid, 
int64_t min_ovlp, asg_arc_t *p)
{
    int64_t tl5, tl3, ext5, ext3, qs = z->qs, rf;
    uint64_t u, v, l, rr; // u: query end; v: target end; l: length from u to v

    ///if query and target are in different strand
    if (z->is_rev) tl5 = tl - z->te, tl3 = z->ts; // tl5: 5'-end overhang (on the query strand); tl3: similar
    else tl5 = z->ts, tl3 = tl - z->te;

    ///ext5 and ext3 is the hang on left side and right side, respectively
    ext5 = qs < tl5? qs : tl5;
    ext3 = (((ql - ((int64_t)z->qe)) < tl3)? (ql - ((int64_t)z->qe)) : tl3);

    if(ext5 > 0 || ext3 > 0) return MA_HT_INT;///overhang
    if (qs <= tl5 && (ql - (int64_t)z->qe) <= tl3) { // query contained in target
        return MA_HT_QCONT; 
    } else if (qs >= tl5 && (ql - (int64_t)z->qe) >= tl3) { // target contained in query
        return MA_HT_TCONT;  
    } else if (qs > tl5) { ///u = 0 means query-to-target overlap, l is the length of node in string graph (not the overlap length)
        u = 0, v = !!(z->is_rev), l = qs - tl5; 
    } else { ///u = 1 means target-to-query overlaps, l is the length of node in string graph (not the overlap length)
        u = 1, v = !(z->is_rev), l = (ql - z->qe) - tl3; 
    }
    if ((int64_t)z->qe - qs + ext5 + ext3 < min_ovlp || (int64_t)z->te - (int64_t)z->ts + ext5 + ext3 < min_ovlp) {
        return MA_HT_SHORT_OVLP; // short overlap
    }
    rf = integer_hit2arc_idx_contain(z, qocc, tocc, &rr);
    if(rf != MA_HT_DOUBLE) return rf;
    if(rr != (uint64_t)-1 && rr != u) {
        // fprintf(stderr, "[M::%s::] z->is_rev::%u, z->qs::%u, z->qe::%u, z->ts::%u, z->te::%u, z->qs_k::%u, z->qe_k::%u, z->ts_k::%u, z->te_k::%u, ql::%ld, tl::%ld, qocc::%ld, tocc::%ld\n", __func__, 
        //     z->is_rev, z->qs, z->qe, z->ts, z->te, z->qs_k, z->qe_k, z->ts_k, z->te_k, ql, tl, qocc, tocc);
        return MA_HT_INT;///overhang
    }
    ///u = 0 / 1 means query-to-target / target-to-query overlaps, 
    ///l is the length of node in string graph (not the overlap length between two reads)
    u |= qid<<1, v |= tid<<1;
    /**
    p->ul: |____________31__________|__________1___________|______________32_____________|
                        qn            direction of overlap       length of this node (not overlap length)
                                        (in the view of query)
    p->v : |___________31___________|__________1___________|
                        tn             reverse direction of overlap
                                      (in the view of target)
    p->ol: overlap length
    **/
    if(p) {
        p->ul = (uint64_t)u<<32 | l, p->v = v, p->ol = ql - l, p->del = 0;
        ///l is the length of node in string graph (not the overlap length)

        p->strong = 1; p->el = 1; p->no_l_indel = 1;
    }
    return l;
}


void integer_append_ug_ovlp(ul_resolve_t *uidx, uint32_t qid, ul2ul_item_t *o, ul2ul_idx_t *ul2)
{
    // if(o->is_del) return;
    uint64_t k; uc_block_t *xi; ul2ul_t *z; ul2ul_item_t *t; int32_t r;
    ma_ug_t *ug = uidx->l1_ug; ul_str_t *str = &(uidx->pstr.str.a[qid]);
    for (k = 0; k < str->cn; k++) {
        xi = &(uidx->idx->a[qid].bb.a[str->a[k]>>32]);
        assert(((xi->hid<<1)+xi->rev)==((uint32_t)str->a[k]));

        ///ul side
        kv_pushp(ul2ul_t, *o, &z);
        z->hid = xi->hid + ul2->uln; z->is_rev = xi->rev; z->is_del = o->is_del; z->is_ct = 0;
        z->qs = xi->qs; z->qe = xi->qe; z->ts = xi->ts; z->te = xi->te; 
        z->qs_k = k; z->qe_k = k + 1; z->ts_k = 0; z->te_k = 1;

        r = integer_hit2arc(z, uidx->idx->a[qid].rlen, ug->u.a[xi->hid].len, uidx->pstr.str.a[qid].cn, 
        1, qid, z->hid, 0, NULL);
        if(r == MA_HT_INT) {
            o->n--; continue;
        }
        ///ug side
        t = get_ul_ovlp(ul2, xi->hid, 0);
        assert(t);
        kv_pushp(ul2ul_t, *t, &z); 
        z->hid = qid; z->is_rev = xi->rev; z->is_del = o->is_del; z->is_ct = 0;
        z->qs = xi->ts; z->qe = xi->te; z->ts = xi->qs; z->te = xi->qe; 
        z->qs_k = 0; z->qe_k = 1; z->ts_k = k; z->te_k = k + 1;
    }
}

void integer_node_del(ul2ul_idx_t *ul2, uint64_t id, uint64_t is_ct)
{
    ul2ul_item_t *o = get_ul_ovlp(ul2, ulg_id(*ul2, id), ulg_type(*ul2, id));
    if(o) {
        uint64_t k; ul2ul_t *z;
        for (k = 0; k < o->cn; k++) {
            if(is_ct) o->a[k].is_ct = 1;
            else o->a[k].is_del = 1;
            z = get_ul_o(ul2, ulg_id(*ul2, o->a[k].hid), ulg_type(*ul2, o->a[k].hid), ulg_id(*ul2, id), ulg_type(*ul2, id));
            if(is_ct) z->is_ct = 1;
            else z->is_del = 1;
        }
        o->is_del = 1;
    }
}

void integer_containment_purge(ul_resolve_t *uidx, uint32_t qid, ul2ul_item_t *q, ul2ul_idx_t *ul2, uint32_t keep_raw_utg)
{
    if(q->is_del) return;
    uint64_t k; int32_t r; ul2ul_item_t *t; ul2ul_t *z; assert(qid == q->id);
    for (k = 0; k < q->cn; k++) {
        if(q->a[k].is_del || q->a[k].is_ct) continue;
        t = get_ul_ovlp(ul2, ulg_id(*ul2, q->a[k].hid), ulg_type(*ul2, q->a[k].hid));
        assert(t);
        if(t->is_del) continue;
        r = integer_hit2arc(&(q->a[k]), ulg_len(*uidx, qid), ulg_len(*uidx, q->a[k].hid), 
        ulg_occ(*uidx, qid), ulg_occ(*uidx, q->a[k].hid), qid, q->a[k].hid, 0, NULL);
        // assert(r != MA_HT_INT);
        if (r == MA_HT_QCONT) {
            q->a[k].is_ct = 1;
            z = get_ul_o(ul2, ulg_id(*ul2, q->a[k].hid), ulg_type(*ul2, q->a[k].hid), ulg_id(*ul2, qid), ulg_type(*ul2, qid));
            assert(z && (!z->is_del) && (!z->is_ct)); z->is_ct = 1;
            if((!keep_raw_utg) || ulg_type(*ul2, qid)) integer_node_del(ul2, qid, 1);
        } else if (r == MA_HT_TCONT) {
            q->a[k].is_ct = 1;
            z = get_ul_o(ul2, ulg_id(*ul2, q->a[k].hid), ulg_type(*ul2, q->a[k].hid), ulg_id(*ul2, qid), ulg_type(*ul2, qid));
            assert(z && (!z->is_del) && (!z->is_ct)); z->is_ct = 1;
            if((!keep_raw_utg) || ulg_type(*ul2, q->a[k].hid)) integer_node_del(ul2, q->a[k].hid, 1);
        }
    }

    // if(ulg_type(*ul2, k)) {
    //     for (k = 0; k < q->cn; k++) {
    //         if((!(q->a[k].is_del)) && (!(q->a[k].is_ct))) break;
    //     }
    //     if(k >= q->cn) q->is_del = 1;
    // }
}


static void worker_integer_correction(void *data, long i, int tid) // callback for kt_for()
{
    ul_resolve_t *uidx = (ul_resolve_t *)data;
    integer_ml_t *sl = &(uidx->str_b); 
    integer_t *buf = &(sl->buf[tid]);
    // uc_block_t *uls; uint64_t uls_n; uint32_t v;
    // uint64_t *srt_a, srt_n, is_circle = ((uidx->psrt.idx.a[i]&((uint64_t)(0x100000000)))?1:0);
    // srt_a = uidx->psrt.srt.a + (uint32_t)uidx->psrt.idx.a[i]; srt_n = uidx->psrt.idx.a[i]>>33;
    // if(srt_n == 0) return;
    // integer_candidate(uidx, srt_a, srt_n, is_circle, buf);
    integer_candidate(uidx, buf, i, (asm_opt.purge_level_primary == 0?1:0));
}

static void worker_integer_postprecess(void *data, long i, int tid) // callback for kt_for()
{
    ul_resolve_t *uidx = (ul_resolve_t *)data;
    integer_ml_t *sl = &(uidx->str_b); 
    integer_t *buf = &(sl->buf[tid]);
    ul2ul_item_t *it = get_ul_ovlp(&(uidx->uovl), i, 1);
    if(!it) return;
    assert(uidx->pstr.str.a[i].cn > 1);
    // uc_block_t *uls; uint64_t uls_n; uint32_t v;
    // uint64_t *srt_a, srt_n, is_circle = ((uidx->psrt.idx.a[i]&((uint64_t)(0x100000000)))?1:0);
    // srt_a = uidx->psrt.srt.a + (uint32_t)uidx->psrt.idx.a[i]; srt_n = uidx->psrt.idx.a[i]>>33;
    // if(srt_n == 0) return;
    // integer_candidate(uidx, srt_a, srt_n, is_circle, buf);
    integer_gen_ovlp(uidx, buf, i, it, uidx->uovl.uln);
}


void gen_integer_normalize(ul_resolve_t *uidx) 
{
    uint64_t k; ul2ul_idx_t *u2o = &(uidx->uovl); ul2ul_item_t *it;
    for (k = 0; k < u2o->uln; k++) {
        it = get_ul_ovlp(u2o, k, 1);
        if(!it) continue;
        assert(uidx->pstr.str.a[k].cn > 1);
        if(it->is_del) continue;
        integer_normalize_ovlp(uidx, k, it, u2o);
    }
}

uint64_t dd_path_connect(asg_t *g, ul_str_t *str)
{
    if(str->cn < 2) return 1;///actually should return 1, doesn't matter
    uint64_t i, v, w, nv, k; asg_arc_t *av;
    v = ((uint32_t)str->a[0]);
    for (i = 1; i < str->cn; i++) {
        w = ((uint32_t)str->a[i]);

        nv = asg_arc_n(g, v); av = asg_arc_a(g, v);
        for (k = 0; k < nv; k++) {
            if(av[k].del) continue;
            if(av[k].v == w) break;
        }
        if(k >= nv) return 0;

        nv = asg_arc_n(g, w^1); av = asg_arc_a(g, w^1);
        for (k = 0; k < nv; k++) {
            if(av[k].del) continue;
            if(av[k].v == (v^1)) break;
        }
        if(k >= nv) return 0;

        v = w;
    }
    return 1;
}

void clip_integer_chimeric(ul_resolve_t *uidx, uint32_t qid, ul2ul_item_t *o, ul2ul_idx_t *ul2, integer_t *buf, int64_t min_dp)
{
    uint64_t k, is_srt = 0, is_del = 0/**, rid, zs, ze**/; assert(o->id == qid);
    ul_str_t *str = &(uidx->pstr.str.a[qid]); ///int64_t z, z_n;
    for (k = 0; k < o->cn; k++) {
        if(o->a[k].is_del) break;
        if(k > 0 && o->a[k].hid < o->a[k-1].hid) break;
    }
    if(k < o->cn) {
        is_srt = 1;
    } else {
        for (; k < o->n; k++) {
            if(!(o->a[k].is_del)) break;
            if(k > o->cn && o->a[k].hid < o->a[k-1].hid) break;
        }
        if(k < o->n) is_srt = 1;
    }

    if(is_srt) {
        for (k = o->cn = 0; k < o->n; k++) {
            if(o->a[k].is_del) o->a[k].hid |= ((uint32_t)(0x80000000));
            else o->cn++;
        }

        radix_sort_ul2ul_srt(o->a, o->a + o->n);

        for (k = 0; k < o->n; k++) {
            if(o->a[k].hid&((uint32_t)(0x80000000))) {
                o->a[k].hid -= ((uint32_t)(0x80000000));
            }
        }
    }

    for (k = 0, buf->u.n = 0; k < o->cn; k++) {
        // kv_push(uint64_t, buf->u, (o->a[k].qs<<1));
        kv_push(uint64_t, buf->u, (o->a[k].qs_k<<1));
        // kv_push(uint64_t, buf->u, (o->a[k].qe<<1)|1);
        kv_push(uint64_t, buf->u, (o->a[k].qe_k<<1)|1);
    }
    radix_sort_srt64(buf->u.a, buf->u.a + buf->u.n);

    int64_t dp, old_dp; uint64_t start, end, b_n = buf->u.n;
    for (k = 0, dp = 0, start = 0; k < b_n; ++k) {
        old_dp = dp;
        ///if a[j] is qe
        if (buf->u.a[k]&1) --dp;
        else ++dp;

        if (old_dp < min_dp && dp >= min_dp) {///old_dp < dp, b.a[j] is qs
            start = buf->u.a[k]>>1;
        } else if (old_dp >= min_dp && dp < min_dp) {///old_dp > min_dp, b.a[j] is qe
            end = buf->u.a[k]>>1;
            kv_push(uint64_t, buf->u, ((start<<32)|(end)));
        }
    }

    is_del = 0;
    if(buf->u.n == b_n) {
        is_del = 1;
    } else {
        uint32_t is_left = 0, is_right = 0, is_middle = 0;
        for (k = b_n; k < buf->u.n; ++k) {
            start = buf->u.a[k]>>32; end = (uint32_t)buf->u.a[k];
            if(start == 0) is_left = 1;
            else is_middle = 1;
            if(end == str->cn/**uidx->idx->a[qid].rlen**/) is_right = 1;
            else is_middle = 1;
        }
        if(is_left == 0 && is_right == 0) is_del = 1;
        if(is_left && is_right && is_middle) is_del = 1;
    }

    // if(is_del && str->cn > 1 && o->is_consist && dd_path_connect(uidx->l1_ug->g, str)) {
    //     is_del = 0;
    // }

    // fprintf(stderr, "[M::%s] qid::%u, o->is_consist::%u, str->cn::%u, is_connect::%lu\n", 
    //     __func__, qid, o->is_consist, str->cn, dd_path_connect(uidx->l1_ug->g, str));

    /**
    if(!is_del) {
        for (k = 0; k < str->cn; k++) {
            rid = (((uint32_t)str->a[k])>>1);
            if(!IF_HOM(rid, *(uidx->bub))) break;
        }

        if(k < str->cn) {///at least a het node covered
            for (k = 0, buf->u.n = 0; k < o->cn; k++) {
                z = o->a[k].qs_k; z_n = o->a[k].qe_k;
                for (; z < z_n; z++) {
                    rid = (((uint32_t)str->a[z])>>1);
                    if(!IF_HOM(rid, *(uidx->bub))) break;
                }
                if(z >= z_n) continue;
                zs = z;

                for (z = z_n - 1; z >= (int64_t)zs; z--) {
                    rid = (((uint32_t)str->a[z])>>1);
                    if(!IF_HOM(rid, *(uidx->bub))) break;
                }
                ze = z + 1;
                assert(zs < ze);

                kv_push(uint64_t, buf->u, (zs<<1));
                kv_push(uint64_t, buf->u, (ze<<1)|1);
                
            }
            radix_sort_srt64(buf->u.a, buf->u.a + buf->u.n);

            b_n = buf->u.n;
            for (k = 0, dp = 0, start = 0; k < b_n; ++k) {
                old_dp = dp;
                ///if a[j] is qe
                if (buf->u.a[k]&1) --dp;
                else ++dp;

                if (old_dp < min_dp && dp >= min_dp) {///old_dp < dp, b.a[j] is qs
                    start = buf->u.a[k]>>1;
                } else if (old_dp >= min_dp && dp < min_dp) {///old_dp > min_dp, b.a[j] is qe
                    end = buf->u.a[k]>>1;
                    kv_push(uint64_t, buf->u, ((start<<32)|(end)));
                }
            }

            is_del = 0;
            if(buf->u.n == b_n) {
                is_del = 1;
            } else {
                uint32_t is_left = 0, is_right = 0, is_middle = 0;
                for (k = b_n; k < buf->u.n; ++k) {
                    start = buf->u.a[k]>>32; end = (uint32_t)buf->u.a[k];
                    if(start == 0) is_left = 1;
                    else is_middle = 1;
                    if(end == str->cn) is_right = 1;
                    else is_middle = 1;
                }
                if(is_left == 0 && is_right == 0) is_del = 1;
                if(is_left && is_right && is_middle) is_del = 1;
            }
        }
    }
    **/

    if(is_del) {
        o->is_del = 1;
        for (k = 0; k < o->cn; k++) o->a[k].is_del = 1;
        if(o->cn && o->n > o->cn) radix_sort_ul2ul_srt(o->a, o->a + o->n);
        o->cn = 0;
    }
}

void clean_srt_integer(ul_resolve_t *uidx, uint32_t qid, ul2ul_item_t *o, ul2ul_idx_t *ul2)
{
    uint64_t k, l, i, m; ul2ul_t *p;
    radix_sort_ul2ul_srt(o->a, o->a + o->n);
    for(k = 1, l = m = o->cn = 0; k <= o->n; k++) {
        if(k == o->n || o->a[k].hid != o->a[l].hid) {
            for (i = l, p = &(o->a[l]); i < k; i++) {
                if(((p->is_del) && (!o->a[i].is_del)) || 
                            ((o->a[i].qe - o->a[i].qs + o->a[i].te - o->a[i].ts) > (p->qe - p->qs + p->te - p->ts))) {
                    p = &(o->a[i]);
                }
            }
            if(!(p->is_del)) o->cn++;
            o->a[m++] = *p;
            l = k;
        }
    }
    o->n = m;
    if(o->n == o->cn) return;

    for (k = o->cn = 0; k < o->n; k++) {
        if(o->a[k].is_del) o->a[k].hid |= ((uint32_t)(0x80000000));
        else o->cn++;
    }
    radix_sort_ul2ul_srt(o->a, o->a + o->n);

    for (k = 0; k < o->n; k++) {
        if(o->a[k].hid&((uint32_t)(0x80000000))) {
            o->a[k].hid -= ((uint32_t)(0x80000000));
        }
    }
}


static void worker_detect_chimeric(void *data, long i, int tid) // callback for kt_for()
{
    ul_resolve_t *uidx = (ul_resolve_t *)data;
    integer_ml_t *sl = &(uidx->str_b); 
    integer_t *buf = &(sl->buf[tid]);
    ul2ul_item_t *it = get_ul_ovlp(&(uidx->uovl), i, 1);
    if(!it) return;
    assert(uidx->pstr.str.a[i].cn > 1);
    if(it->is_del) return;
    clip_integer_chimeric(uidx, i, it, &(uidx->uovl), buf, 1);
}


void chimeric_integer_deal(ul_resolve_t *uidx)
{
    uint64_t k; ul2ul_idx_t *u2o = &(uidx->uovl); ul2ul_item_t *it;
    kt_for(uidx->str_b.n_thread, worker_detect_chimeric, uidx, uidx->idx->n);
    for (k = 0; k < u2o->uln; k++) {
        it = get_ul_ovlp(u2o, k, 1);
        if(!it) continue;
        assert(uidx->pstr.str.a[k].cn > 1);
        integer_normalize_ovlp_purge(uidx, k, it, u2o);
    }
}




static void worker_integert_clean(void *data, long i, int tid) // callback for kt_for()
{
    ul_resolve_t *uidx = (ul_resolve_t *)data;
    ul2ul_item_t *it = get_ul_ovlp(&(uidx->uovl), ulg_id(uidx->uovl, (uint32_t)i), 
                                                            ulg_type(uidx->uovl, (uint32_t)i));
    if(!it) return;
    if((uint32_t)i >= uidx->uovl.uln) it->id = i;
    clean_srt_integer(uidx, i, it, &(uidx->uovl));
}

static void worker_integert_debug_sym(void *data, long i, int tid) // callback for kt_for()
{
    ul_resolve_t *uidx = (ul_resolve_t *)data;
    ul2ul_item_t *o = get_ul_ovlp(&(uidx->uovl), ulg_id(uidx->uovl, (uint32_t)i), 
                                                            ulg_type(uidx->uovl, (uint32_t)i));
    if(!o) return;
    ul2ul_idx_t *ul2 = &(uidx->uovl); 
    // if(o->is_del) {
    //     assert(o->cn == 0);
    // }

    ul2ul_t *z; uint64_t k, qid = i, ct_n = 0, del_n = 0; assert(o->id == qid);
    for (k = 0; k < o->n; k++) {
        if(o->a[k].is_del) del_n++;
        else if(o->a[k].is_ct) ct_n++;
        assert(ul2->item_idx[o->a[k].hid] != (uint32_t)-1);
        z = get_ul_o(ul2, ulg_id(*ul2, o->a[k].hid), ulg_type(*ul2, o->a[k].hid), 
                                            ulg_id(*ul2, qid), ulg_type(*ul2, qid));
        assert(z && z->hid == qid && z->is_del == o->a[k].is_del && z->is_ct == o->a[k].is_ct);

        // if(!(o->a[k].is_rev == z->is_rev && (*z).qs == o->a[k].ts && (*z).qe == o->a[k].te && 
        //                                                 (*z).ts == o->a[k].qs && (*z).te == o->a[k].qe)){
        //         fprintf(stderr, "uln::%u[M::%s] fn::%u(%c), fqs::%u, fqe::%u, fts::%u, fte::%u, fdel::%u, rn::%u(%c), rqs::%u, rqe::%u, rts::%u, rte::%u, rdel::%u\n", ul2->uln, __func__, 
        //         o->a[k].hid, "+-"[o->a[k].is_rev], o->a[k].qs, o->a[k].qe, o->a[k].ts, o->a[k].te, o->a[k].is_del,
        //         z->hid, "+-"[z->is_rev], z->qs, z->qe, z->ts, z->te, z->is_del);
        // }
        // if(o->a[k].is_rev == z->is_rev && (*z).qs == o->a[k].ts && (*z).qe == o->a[k].te && 
        //     (*z).ts == o->a[k].qs && (*z).te == o->a[k].qe) {
        //     fprintf(stderr, "good::[M::%s] fn::%u(%c), fqs::%u, fqe::%u, fts::%u, fte::%u, rn::%u(%c), rqs::%u, rqe::%u, rts::%u, rte::%u\n", __func__, 
        //         o->a[k].hid, "+-"[o->a[k].is_rev], o->a[k].qs, o->a[k].qe, o->a[k].ts, o->a[k].te, 
        //         z->hid, "+-"[z->is_rev], z->qs, z->qe, z->ts, z->te);
        // } else {
        //     fprintf(stderr, "bad::[M::%s] fn::%u(%c), fqs::%u, fqe::%u, fts::%u, fte::%u, rn::%u(%c), rqs::%u, rqe::%u, rts::%u, rte::%u\n", __func__, 
        //         o->a[k].hid, "+-"[o->a[k].is_rev], o->a[k].qs, o->a[k].qe, o->a[k].ts, o->a[k].te, 
        //         z->hid, "+-"[z->is_rev], z->qs, z->qe, z->ts, z->te);
        // }
        assert(o->a[k].is_rev == z->is_rev && (*z).qs == o->a[k].ts && (*z).qe == o->a[k].te && 
                                                        (*z).ts == o->a[k].qs && (*z).te == o->a[k].qe);
        if(k<o->cn) assert(!o->a[k].is_del);
        else assert(o->a[k].is_del);
    }

    if(o->is_del) assert((ct_n + del_n) == o->n);
}

void print_primary_ul_chain(ul_str_t *p_str, ul_vec_t *ul, ma_ug_t *ug)
{
    uc_block_t *z; uint64_t k, dd = 0;
    if(p_str) {
        for (k = 0; k < p_str->cn; k++) {
            z = &(ul->bb.a[p_str->a[k]>>32]);
            fprintf(stderr, "[M::%s::pk->%lu] utg%.6d%c(%c::len->%u), qs::%u, qe::%u, ts::%u, te::%u\n", __func__, k,
                    z->hid+1, "lc"[ug->u.a[z->hid].circ], "+-"[z->rev], ug->u.a[z->hid].len, 
                    z->qs, z->qe, z->ts, z->te);
            dd++;
        }
    } else {
        for (k = 0; k < ul->bb.n; k++) {
            z = &(ul->bb.a[k]);
            if(!z->pchain) continue;
            fprintf(stderr, "[M::%s::rk->%lu] utg%.6d%c(%c::len->%u), qs::%u, qe::%u, ts::%u, te::%u\n", __func__, k,
                    z->hid+1, "lc"[ug->u.a[z->hid].circ], "+-"[z->rev], ug->u.a[z->hid].len, 
                    z->qs, z->qe, z->ts, z->te);
            dd++;
        }
    }

    if(dd) fprintf(stderr, "*******************************************\n");
}

void push_integer_seq_exact(ul_vec_t *res, ma_ug_t *ug, uint64_t *seq, uint64_t seq_n, uint64_t *off, uint32_t tid)
{
    if(seq_n <= 0) return;
    uint64_t k, v, ql, ul; uc_block_t *x;
    // fprintf(stderr, "[M::%s::] old_len::%u, new_len::%u\n", __func__, res->rlen, (uint32_t)off[seq_n-1]);
    res->rlen = (uint32_t)off[seq_n-1]; res->bb.n = 0; kv_resize(uc_block_t, res->bb, seq_n);
    for (k = 0; k < seq_n; k++) {
        v = (uint32_t)seq[k];
        kv_pushp(uc_block_t, res->bb, &x);
        x->hid = v>>1; x->rev = !!(v&1); x->pchain = 1; x->el = 1; x->base = 0;
        x->qs = off[k]>>32; x->qe = (uint32_t)off[k]; 
        ql = x->qe - x->qs; ul = ug->g->seq[x->hid].len;
        if((ul == ql) || (k > 0 && k + 1 < seq_n)) {
            x->ts = 0; x->te = ul;
        } else {
            // if(!(k == 0 || k + 1 == seq_n)){
            //     fprintf(stderr, "[M::%s::] k::%lu, ql::%lu, ul::%lu\n", __func__, k, ql, ul);
            // }
            // assert(k == 0 || k + 1 == seq_n);
            if(k == 0) {
                if(x->rev) {
                    x->ts = 0; x->te = MIN(ql, ul);                 
                } else {
                    x->ts = ul - MIN(ql, ul); x->te = ul;
                }
            }

            if(k + 1 == seq_n) {
                if(!x->rev) {
                    x->ts = 0; x->te = MIN(ql, ul);
                } else {
                    x->ts = ul - MIN(ql, ul); x->te = ul;
                }
            }
        }
        x->pidx = x->aidx = x->pdis = (uint32_t)-1;
        if(k > 0 && (seq[k]&((uint64_t)(0x8000000000000000)))) {
            x->pidx = k - 1; x->pdis = (seq[k]<<1)>>33;
            res->bb.a[x->pidx].aidx = k;
        }
    }
    assert(res->rlen == x->qe);
    // res->rlen = ((uint32_t)-1) - tid;
    // print_primary_ul_chain(NULL, res, ug);
}

uint32_t double_check_gconnect(asg_t *g, uint32_t v, uint32_t w, uint32_t v2w_d)
{
    asg_arc_t *av = asg_arc_a(g, v); uint32_t nv = asg_arc_n(g, v), k;
    for (k = 0; k < nv; k++) {
        if(av[k].del) continue;
        if(av[k].v == w) break;
    }
    if(k < nv) return 1;
    return 0;
}


void update_integer_seq(ul_resolve_t *uidx, integer_t *buf, uint32_t id, uint64_t *seq, uint64_t seq_n, uint32_t tid)
{
    // if((((uint32_t)-1) - uidx->idx->a[id].rlen) < asm_opt.thread_num) {
    //     fprintf(stderr, "id->%u, c_tid->%u, l_tid->%u\n", id, tid, (((uint32_t)-1) - uidx->idx->a[id].rlen));
    //     exit(0);
    // }
    // assert(uidx->idx->a[id].rlen!=(uint32_t)-1);
    // fprintf(stderr, "dd->%u\n", uidx->idx->a[id].dd);
    // if(id != 1487) return;
    if(seq_n <= 0) return;
    buf->n_correct++;
    ul_str_t *str = &(uidx->pstr.str.a[id]); uc_block_t *z; ul_chain_t sc, msc;
    uint64_t k, i, pp; integer_aln_t *b; ma_ug_t *ug = uidx->l1_ug;
    assert((seq[0]>>32) == (ug->g->seq[((uint32_t)seq[0])>>1].len));
    // fprintf(stderr, "\n[M::%s::id->%u] seq_n::%lu\n", __func__, id, seq_n);
    // print_primary_ul_chain(str, &(uidx->idx->a[id]), ug);
    if(seq_n == 1) {
        pp = ug->g->seq[((uint32_t)seq[0])>>1].len;
        push_integer_seq_exact(&(uidx->idx->a[id]), ug, seq, seq_n, &pp, tid);
        return;
    }
    
    buf->u.n = 0;
    for (k = 0; k < str->cn; k++) {///old seq
        pp = ((uint32_t)str->a[k])>>1; pp <<= 33; pp += ((uint64_t)(0x100000000));
        pp += (k<<1); pp += ((uint32_t)str->a[k])&1;
        kv_push(uint64_t, buf->u, pp);
        // fprintf(stderr, "ok::%lu, ov::%u\n", k, ((uint32_t)str->a[k]));
    }
    for (k = 0; k < seq_n; k++) {///new seq
        pp = ((uint32_t)seq[k])>>1; pp <<= 33; 
        pp += (k<<1); pp += ((uint32_t)seq[k])&1;
        kv_push(uint64_t, buf->u, pp);
        // fprintf(stderr, "nk::%lu, nv::%u\n", k, ((uint32_t)seq[k]));

        // fprintf(stderr, ">>>>>>[M::%s::c_k->%lu] utg%.6d%c(%c::len->%u)\n", __func__, k,
        //             (((uint32_t)seq[k])>>1)+1, "lc"[ug->u.a[(((uint32_t)seq[k])>>1)].circ], 
        //             "+-"[((uint32_t)seq[k])&1], ug->u.a[(((uint32_t)seq[k])>>1)].len);
    }

    radix_sort_srt64(buf->u.a, buf->u.a + buf->u.n);
    for (k = 0, buf->b.n = 0; k < buf->u.n; k++) {
        if(buf->u.a[k]&((uint64_t)(0x100000000))) continue;///skip nodes in the old seq
        for (i = k+1; (i < buf->u.n) && ((buf->u.a[k]>>33) == (buf->u.a[i]>>33)); i++) {
            if((buf->u.a[i]&((uint64_t)(0x100000000))) == 0) continue;///skip nodes in the new seq
            if((buf->u.a[k]&1) != (buf->u.a[i]&1)) continue;///ignore reverse alignment;
            ///a[k] is the new seq (q); a[i] is old seq (t)
            kv_pushp(integer_aln_t, buf->b, &b);
            b->vq = ((buf->u.a[k]>>33)<<1) + (buf->u.a[k]&1);
            b->tn_rev_qk = ((uint32_t)buf->u.a[k])>>1;
            b->tk = ((uint32_t)buf->u.a[i])>>1; 
            z = &(uidx->idx->a[id].bb.a[str->a[b->tk]>>32]);
            assert(((z->hid<<1)+z->rev)==(((buf->u.a[i]>>33)<<1) + (buf->u.a[i]&1)));
            if((buf->u.a[k]&1) != (buf->u.a[i]&1)) {
                b->tk = str->cn - b->tk - 1;
                b->tn_rev_qk |= ((uint64_t)(0x100000000));
            }            
            b->sc = ug_occ_w(z->ts, z->te, &(ug->u.a[z->hid]));
            assert(b->sc > 0);
        }
    }

    msc.v = msc.s = msc.e = (uint32_t)-1; msc.sc = 0;
    if(buf->b.n > 0) {
        radix_sort_integer_aln_t_srt(buf->b.a, buf->b.a + buf->b.n); ///sorted by rev|qk
        for (k = 1, i = 0; k <= buf->b.n; k++) {
            // fprintf(stderr, "[M::%s::z->%lu] rev::%u, qk::%u, tk::%u\n", 
            // __func__, k-1, (uint32_t)(!!(buf->b.a[k-1].tn_rev_qk>>32)), (uint32_t)buf->b.a[k-1].tn_rev_qk, buf->b.a[k-1].tk);
            if(k == buf->b.n || (buf->b.a[i].tn_rev_qk>>32) != (buf->b.a[k].tn_rev_qk>>32)) {
                if(integer_chain(id, buf->b.a + i, k - i, i, buf, ug, NULL, NULL, &sc) && sc.v != (uint32_t)-1) {
                    if(msc.s == (uint32_t)-1 || msc.sc < sc.sc) msc = sc;
                }
                i = k;
            }
        }
    }
    // fprintf(stderr, "[M::%s::id->%u] seq_n::%lu, align_n::%u, buf->b.n::%u\n", 
    // __func__, id, seq_n, msc.e - msc.s, (uint32_t)buf->b.n);

    
    uint64_t l, v, pd, uls, ule, p_ls, p_le, is_rev, qk, tk, bc; integer_aln_t *x; uc_block_t *tb;
    buf->u.n = 0; kv_resize(uint64_t, buf->u, seq_n); 
    // assert((seq[0]>>32) == (ug->g->seq[((uint32_t)seq[0])>>1].len));
    for (k = l = 0, p_ls = p_le = (uint64_t)-1; k < seq_n; k++) {
        v = (uint32_t)seq[k]; pd = (seq[k]<<1)>>33;
        ule = l + pd; uls = ((ule >= ug->g->seq[v>>1].len)?(ule - ug->g->seq[v>>1].len):(0));
        bc = 0;
        if(k > 0){
            bc = (!!(seq[k]&((uint64_t)(0x8000000000000000))));///if connected in the graph
            if(bc) bc = double_check_gconnect(ug->g, ((uint32_t)seq[k])^1, ((uint32_t)seq[k-1])^1, pd);
            // fprintf(stderr, "bc->%lu\n", bc);
        }
        if(p_ls != (uint64_t)-1) {
            ///uls should uls>=p_ls && uls<=p_le
            if(uls < p_ls) uls = p_ls;
            if(bc && uls > p_le) uls = p_le;
            ///ule should ule > p_le
            if(ule <= p_le) ule = p_le + 1;
        }
        buf->u.a[k] = uls; buf->u.a[k] <<= 32; buf->u.a[k] |= ule;
        p_ls = uls; p_le = ule;
        l = ule;
        // fprintf(stderr, "[init_k->%lu] uls::%lu, ule::%lu, pd::%lu\n", k, uls, ule, pd);
    }
    buf->u.n = seq_n;

    int64_t beg_nl, end_nl; ma_utg_t *u;
    v = (uint32_t)seq[0]; u = &(ug->u.a[v>>1]); beg_nl = u->len;
    v = (uint32_t)seq[seq_n-1]; u = &(ug->u.a[v>>1]); end_nl = u->len;
    // fprintf(stderr, "+[M::%s] beg_nl::%lu, end_nl::%lu\n", __func__, beg_nl, end_nl);

    ///q -> new seq; t -> old seq
    if(msc.s != (uint32_t)-1 && msc.e > msc.s) {
        int64_t qoff, toff, beg_cut = 0, end_cut = 0; uint32_t r_end;
        is_rev = ((buf->b.a[msc.s].tn_rev_qk>>32)&1);

        ///beg
        x = &(buf->b.a[msc.s]); 
        qk = (uint32_t)(x->tn_rev_qk); tk = ((is_rev == 0)? (x->tk):(str->cn - x->tk - 1));
        tb = &(uidx->idx->a[id].bb.a[str->a[tk]>>32]);
        if((!!tb->rev) == (!!is_rev)) {
            if(tb->te == ug->u.a[tb->hid].len) r_end = 1;
            else r_end = 0;
        } else {
            if(tb->ts == 0) r_end = 1;
            else r_end = 0;
        }
        if(r_end) {///start from right end
            qoff = (uint32_t)buf->u.a[qk]; toff = (is_rev?(uidx->idx->a[id].rlen-tb->qs):(tb->qe));
        } else {///start from left end
            qoff = buf->u.a[qk]>>32; toff = (is_rev?(uidx->idx->a[id].rlen-tb->qe):(tb->qs));
        }
        
        if(qk == 0) toff = tb->te - tb->ts;
        if(toff < qoff) beg_cut = qoff - toff;
        // fprintf(stderr, "[M::%s::beg::is_rev->%lu] qoff::%ld, toff::%ld, r_end::%u\n", __func__, is_rev, qoff, toff, r_end);
        

        ///end
        x = &(buf->b.a[msc.e-1]);
        qk = (uint32_t)(x->tn_rev_qk); tk = ((is_rev == 0)? (x->tk):(str->cn - x->tk - 1));
        tb = &(uidx->idx->a[id].bb.a[str->a[tk]>>32]);
        if((!!tb->rev) == (!!is_rev)) {
            if(tb->ts == 0) r_end = 0;
            else r_end = 1;
        } else {
            if(tb->te == ug->u.a[tb->hid].len) r_end = 0;
            else r_end = 1;
        }
        if(!r_end) {
            qoff = l - (buf->u.a[qk]>>32); toff = (is_rev?(tb->qe):(uidx->idx->a[id].rlen-tb->qs));
        } else {
            qoff = l - (uint32_t)(buf->u.a[qk]); toff = (is_rev?(tb->qs):(uidx->idx->a[id].rlen-tb->qe));
        }
        if(qk+1==seq_n) toff = tb->te - tb->ts;
        if(toff < qoff) end_cut = qoff - toff;
        // fprintf(stderr, "[M::%s::end::is_rev->%lu] qoff::%ld, toff::%ld, r_end::%u\n", __func__, is_rev, qoff, toff, r_end);

        if(beg_cut > 0) {
            v = (uint32_t)seq[0]; u = &(ug->u.a[v>>1]);
            if(beg_cut < u->len) beg_cut = u->len - beg_cut;
            else beg_cut = 0;
            if(v&1) { 
                beg_nl = (int64_t)(Get_READ_LENGTH(R_INF, (u->a[0]>>33)));
            } else {
                beg_nl = (int64_t)(Get_READ_LENGTH(R_INF, (u->a[u->n-1]>>33)));
            }
            if(beg_cut > beg_nl) beg_nl = beg_cut;
            if(beg_nl == u->len) beg_cut = 0;
            else beg_cut = 1;
        }

        
        if(end_cut > 0) {
            v = (uint32_t)seq[seq_n-1]; u = &(ug->u.a[v>>1]);
            if(end_cut < u->len) end_cut = u->len - end_cut;
            else end_cut = 0;
            if(v&1) { 
                end_nl = Get_READ_LENGTH(R_INF, (u->a[u->n-1]>>33));
            } else {
                end_nl = Get_READ_LENGTH(R_INF, (u->a[0]>>33));
            }
            if(end_cut > end_nl) end_nl = end_cut;
            if(end_nl == u->len) end_cut = 0;
            else end_cut = 1;
        }
        // fprintf(stderr, "++[M::%s] beg_nl::%lu, end_nl::%lu\n", __func__, beg_nl, end_nl);
        if(beg_cut || end_cut) {
            k = 0; uls = 0; ule = beg_nl; 
            buf->u.a[k] = uls; buf->u.a[k] <<= 32; buf->u.a[k] |= ule;
            l = ule; p_ls = uls; p_le = ule;
            for (k += 1; k + 1 < seq_n; k++) {
                v = (uint32_t)seq[k]; pd = (seq[k]<<1)>>33;
                ule = l + pd; uls = ((ule >= ug->g->seq[v>>1].len)?(ule - ug->g->seq[v>>1].len):(0));
                bc = 0;
                if(k > 0){
                    bc = (!!(seq[k]&((uint64_t)(0x8000000000000000))));///if connected in the graph
                    if(bc) bc = double_check_gconnect(ug->g, ((uint32_t)seq[k])^1, ((uint32_t)seq[k-1])^1, pd);
                    // fprintf(stderr, "bc->%lu\n", bc);
                }
                if(p_ls != (uint64_t)-1) {
                    ///uls should uls>=p_ls && uls<=p_le
                    if(uls < p_ls) uls = p_ls;
                    if(bc && uls > p_le) uls = p_le;
                    ///ule should ule > p_le
                    if(ule <= p_le) ule = p_le + 1;
                }
                buf->u.a[k] = uls; buf->u.a[k] <<= 32; buf->u.a[k] |= ule;
                p_ls = uls; p_le = ule;
                l = ule;
            }

            assert(k < seq_n);
            v = (uint32_t)seq[k]; pd = (seq[k]<<1)>>33;
            ule = l + pd; uls = ((ule >= ug->g->seq[v>>1].len)?(ule - ug->g->seq[v>>1].len):(0));
            ule = uls + end_nl;
            bc = 0;
            if(k > 0){
                bc = (!!(seq[k]&((uint64_t)(0x8000000000000000))));///if connected in the graph
                if(bc) bc = double_check_gconnect(ug->g, ((uint32_t)seq[k])^1, ((uint32_t)seq[k-1])^1, pd);
                // fprintf(stderr, "bc->%lu\n", bc);
            }
            if(p_ls != (uint64_t)-1) {
                ///uls should uls>=p_ls && uls<=p_le
                if(uls < p_ls) uls = p_ls;
                if(bc && uls > p_le) uls = p_le;
                ///ule should ule > p_le
                if(ule <= p_le) ule = p_le + 1;
            }
            if(ule - uls < (uint64_t)end_nl) ule = uls + end_nl;
            buf->u.a[k] = uls; buf->u.a[k] <<= 32; buf->u.a[k] |= ule;
            p_ls = uls; p_le = ule;
            l = ule;
        }
    }
    
    push_integer_seq_exact(&(uidx->idx->a[id]), ug, seq, seq_n, buf->u.a, tid);
}

static void worker_integer_update(void *data, long i, int tid) // callback for kt_for()
{
    ul_resolve_t *uidx = (ul_resolve_t *)data;
    integer_ml_t *sl = &(uidx->str_b); 
    integer_t *buf = &(sl->buf[i]);///normally should be buf = &(sl->buf[tid])
    uint64_t k, l;
    for (k = 0, l = (uint64_t)-1; k <= buf->res_dump.n; k++) {
        if(k == buf->res_dump.n || (buf->res_dump.a[k]&((uint64_t)(0xffffffff))) == ((uint64_t)(0xffffffff))) {
            if(l != (uint64_t)-1 && k - l > 1) {
                update_integer_seq(uidx, buf, buf->res_dump.a[l]>>32, buf->res_dump.a + l + 1, k - l - 1, tid);
            }
            l = k;
        }
    }
}


void integer_correction(ul_resolve_t *uidx)
{
    integer_ml_t sl; 
    init_integer_ml_t(&sl, uidx, asm_opt.thread_num);
}

uint64_t clean_ul_re_correct_buf(ul_resolve_t *uidx, uint64_t clean_dump, uint64_t *tot_circle)
{
    uint64_t k, occ, n_circle;
    for (k = occ = n_circle = 0; k < uidx->str_b.n_thread; k++) {
        occ += uidx->str_b.buf[k].n_correct; 
        n_circle += uidx->str_b.buf[k].n_circle;
        uidx->str_b.buf[k].n_correct = 0;
        uidx->str_b.buf[k].n_circle = 0;
        if(clean_dump) uidx->str_b.buf[k].res_dump.n = 0;
        uidx->str_b.buf[k].q.n = 0;
        uidx->str_b.buf[k].t.n = 0;
        uidx->str_b.buf[k].b.n = 0;
        uidx->str_b.buf[k].f.n = 0;
        uidx->str_b.buf[k].p.n = 0;
        uidx->str_b.buf[k].o.n = 0;
        uidx->str_b.buf[k].u.n = 0;
        uidx->str_b.buf[k].vis.n = 0;
        uidx->str_b.buf[k].sc.n = 0;
        uidx->str_b.buf[k].snp.n = 0;
    }
    (*tot_circle) = n_circle;
    return occ;
}


void rebuid_idx(ul_resolve_t *uidx)
{
    uint32_t k;
    free(uidx->idx->ridx.idx.a); free(uidx->idx->ridx.occ.a); 
    memset(&(uidx->idx->ridx), 0, sizeof((uidx->idx->ridx)));
    filter_ul_ug(uidx->l1_ug);
	gen_ul_vec_rid_t(uidx->idx, NULL, uidx->l1_ug);
    update_ug_arch_ul_mul(uidx->l1_ug);

    free(uidx->pstr.idx.a); free(uidx->pstr.occ.a); 
    for (k = 0; k < uidx->pstr.str.n; k++) {
        free(uidx->pstr.str.a[k].a);
    }
    free(uidx->pstr.str.a); 
    memset(&(uidx->pstr), 0, sizeof(uidx->pstr));
    init_ul_str_idx_t(uidx);
}

void shrink_1b(ma_ug_t *ug, uc_block_t *z, uc_block_t *lim, uint32_t is_forward)
{
    if(z->ts != 0 || z->te != ug->g->seq[z->hid].len) return;
    uc_block_t bc = *z;
    if(z->ts + 1 < z->te && z->qs + 1 < z->qe) {
        uint32_t off = get_offset_adjust(1, z->te-z->ts, z->qe-z->qs);
        if(is_forward) {
            if((!z->rev)) {
                z->ts += 1; z->qs += off;
            } else {
                z->te -= 1; z->qs += off;
            }
            if((lim) && (!((z->qs <= lim->qs) && (z->qe <= lim->qe)))) *z = bc;
        } else {
            if((!z->rev)) {
                z->te -= 1; z->qe -= off;
            } else {
                z->ts += 1; z->qe -= off;
            }
            if((lim) && (!((z->qs >= lim->qs) && (z->qe >= lim->qe)))) *z = bc;
        }
        if((ugl_cover_check(bc.ts, bc.te, &(ug->u.a[bc.hid]))) && 
                                            (!ugl_cover_check(z->ts, z->te, &(ug->u.a[z->hid])))) {
            *z = bc;
        }
    }
    
}

void renew_ul_vec_t(ul_vec_t *x, ma_ug_t *ug)
{
    if(x->bb.n <= 0) return;
    int64_t l, k, z, bn = x->bb.n, dd; int64_t qs = x->bb.a[0].qs, qe;
    for (l = 0, k = 1, qs = 0; k <= bn; k++) {
        if(k == bn || x->bb.a[k].pidx == (uint32_t)-1) {
            if(k - l > 0) {
                for (z = k - 1, dd = 0; z >= l; z--) {
                    if(z > l) {
                        assert(x->bb.a[z].pidx == z-1);
                        dd += x->bb.a[z].pdis;
                    } else {
                        assert(x->bb.a[z].pidx == (uint32_t)-1);
                        dd += ug->g->seq[x->bb.a[z].hid].len;
                    }
                    dd -= (int64_t)(ug->g->seq[x->bb.a[z].hid].len - (x->bb.a[z].te - x->bb.a[z].ts));
                }
                if(dd < 0) dd = 0;
                qe = qs + dd; 
                if(k < bn) {
                    dd = qe + (int64_t)(x->bb.a[k].qs) - (int64_t)(x->bb.a[k-1].qe);
                    if(dd < 0) dd = 0;
                }
                for (z = k - 1; z >= l; z--) {
                    qs = qe - (int64_t)(x->bb.a[z].te - x->bb.a[z].ts);
                    if(qs < 0) qs = 0;
                    x->bb.a[z].qs = qs; x->bb.a[z].qe = qe; 
                    if(z > l) {
                        qe += (int64_t)(ug->g->seq[x->bb.a[z].hid].len - (x->bb.a[z].te - x->bb.a[z].ts));
                        qe -= (int64_t)(x->bb.a[z].pdis);
                        if(qe < 0) qe = 0;
                    }
                }

                if(k < bn) qs = dd;
            }
            l = k;
        }
    }
    x->rlen = x->bb.a[bn-1].qe;
}

void shrink_ul0(all_ul_t *uls, ul_str_t *str, uint64_t id, integer_t *buf, ma_ug_t *ug)
{
    uint32_t k, c_k, p_k, cv, pv, bl, i; uc_block_t *xi, *yi; buf->u.n = 0; nid_t *np = NULL;
    asg_arc_t *av; uint32_t nv, s, e, m, d, mm, is_conn; uint64_t *z; ul_vec_t *x;
    if(str->cn < 2) return;
    for (k = 0, bl = 0, c_k = p_k = pv = (uint32_t)-1; k < str->cn; k++) {
        xi = &(uls->a[id].bb.a[str->a[k]>>32]); is_conn = 0;
        c_k = k; cv = ((uint32_t)str->a[k])^1;
        assert(((xi->hid<<1)+xi->rev)==((uint32_t)str->a[k]));
        if(p_k != (uint32_t)-1) {
            if(xi->pidx == (str->a[p_k]>>32)) {
                av = asg_arc_a(ug->g, cv); nv = asg_arc_n(ug->g, cv);
                for (i = 0; i < nv; i++) {
                    if(av[i].del) continue;
                    if(av[i].v == pv) break;
                }
                if(i < nv) {
                    s = (uint32_t)av[i].ul;
                    if(s == xi->pdis) {
                        is_conn = 1;
                    } else {
                        d = (s>=xi->pdis?s-xi->pdis:xi->pdis-s);
                        mm = MAX(s, xi->pdis);
                        if((d <= (mm*0.08)) || (d <= 512)) is_conn = 1;
                    }
                }
            } else if(xi->pidx == (uint32_t)-1) {
                is_conn = 1;
            }
        }
        if(is_conn) {
            is_conn = 0; assert(k);
            yi = &(uls->a[id].bb.a[str->a[k-1]>>32]);
            if((xi->qs >= yi->qs) && (xi->qe >= yi->qe)) is_conn = 1;
        }
        if(is_conn) {
            bl++;
        } else {
            if(bl > 0) {
                kv_pushp(uint64_t, buf->u, &z);
                (*z) = k - 1 - bl; (*z) <<= 32; (*z) += bl + 1;
            }
            bl = 0;
        }
        p_k = c_k; pv = cv;
    }

    if(bl > 0) {
        kv_pushp(uint64_t, buf->u, &z);
        (*z) = k - 1 - bl; (*z) <<= 32; (*z) += bl + 1;
    }
    // fprintf(stderr, "[M::%s::] buf->u.n::%u, str->cn::%u\n", __func__, buf->u.n, str->cn);
    if(buf->u.n <= 0) {
        str->cn = str->n = 0;
    } else if(buf->u.n > 0) {
        uint32_t pidx, aidx, pdis;
        for (i = 0; i + 1 < buf->u.n; i++) {
            s = (buf->u.a[i]>>32); e = s + ((uint32_t)buf->u.a[i]); 
            ///new id
            kv_pushp(nid_t, uls->nid, &np);
            np->n = uls->nid.a[id].n; MALLOC(np->a, np->n+1);
            memcpy(np->a, uls->nid.a[id].a, np->n); np->a[np->n] = '\0';
            ///new ovlps
            kv_pushp(ul_vec_t, *uls, &x); memset(x, 0, sizeof(*x));
            MALLOC(x->bb.a, e - s); x->bb.n = x->bb.m = e - s;

            for (k = s, m = 0; k < e; k++, m++) {
                pidx = uls->a[id].bb.a[str->a[k]>>32].pidx; 
                aidx = uls->a[id].bb.a[str->a[k]>>32].aidx; 
                pdis = uls->a[id].bb.a[str->a[k]>>32].pdis;
                x->bb.a[m] = uls->a[id].bb.a[str->a[k]>>32];
                if(k == s) {
                    x->bb.a[m].pidx = x->bb.a[m].pdis = (uint32_t)-1;
                } else if(pidx != (uint32_t)-1) {
                    x->bb.a[m].pidx = m - 1; x->bb.a[m].pdis = pdis;
                }
                if(k + 1 == e) {
                    x->bb.a[m].aidx = (uint32_t)-1;
                } else if(aidx != (uint32_t)-1) {
                    x->bb.a[m].aidx = m + 1;
                }
            }
            x->bb.n = m;
            assert(x->bb.n > 1);

            shrink_1b(ug, &(x->bb.a[0]), ((x->bb.n>=2)?&(x->bb.a[1]):(NULL)), 1);
            shrink_1b(ug, &(x->bb.a[x->bb.n-1]), ((x->bb.n>=2)?&(x->bb.a[x->bb.n-2]):(NULL)), 0);
            
            d = x->bb.a[0].qs;
            for (k = 0; k < x->bb.n; k++) {
                x->bb.a[k].qs -= d; x->bb.a[k].qe -= d;
            }
            x->rlen = x->bb.a[x->bb.n-1].qe;
            renew_ul_vec_t(x, ug);
        }

        ///the last one; update in-place
        s = (buf->u.a[i]>>32); e = s + ((uint32_t)buf->u.a[i]); x = &(uls->a[id]);
        for (k = s, m = 0; k < e; k++, m++) {
            pidx = x->bb.a[str->a[k]>>32].pidx; 
            aidx = x->bb.a[str->a[k]>>32].aidx; 
            pdis = x->bb.a[str->a[k]>>32].pdis;
            x->bb.a[m] = x->bb.a[str->a[k]>>32];
            if(k == s) {
                x->bb.a[m].pidx = x->bb.a[m].pdis = (uint32_t)-1;
            } else if(pidx != (uint32_t)-1) {
                x->bb.a[m].pidx = m - 1; x->bb.a[m].pdis = pdis;
            }
            if(k + 1 == e) {
                x->bb.a[m].aidx = (uint32_t)-1;
            } else if(aidx != (uint32_t)-1) {
                x->bb.a[m].aidx = m + 1;
            }
        }
        x->bb.n = m;
        assert(x->bb.n > 1);
        shrink_1b(ug, &(x->bb.a[0]), ((x->bb.n>=2)?&(x->bb.a[1]):(NULL)), 1);
        shrink_1b(ug, &(x->bb.a[x->bb.n-1]), ((x->bb.n>=2)?&(x->bb.a[x->bb.n-2]):(NULL)), 0);
        d = x->bb.a[0].qs;
        for (k = 0; k < x->bb.n; k++) {
            x->bb.a[k].qs -= d; x->bb.a[k].qe -= d;
        }
        x->rlen = x->bb.a[x->bb.n-1].qe;
        renew_ul_vec_t(x, ug);
    }
}

void shrink_uls(ul_resolve_t *uidx)
{
    all_ul_t *uls = uidx->idx; uint64_t k, ul_n = uls->n; 
    for (k = 0; k < ul_n; k++) {
        // fprintf(stderr, "+[M::%s::k->%lu] ul_n::%lu, uls->n::%u\n", __func__, k, ul_n, uls->n);
        shrink_ul0(uls, &(uidx->pstr.str.a[k]), k, &(uidx->str_b.buf[0]), uidx->l1_ug);
        // fprintf(stderr, "-[M::%s::k->%lu] ul_n::%lu, uls->n::%u\n", __func__, k, ul_n, uls->n);
    }
    rebuid_idx(uidx);
}

void print_ul_seq(ul_resolve_t *uidx, uint64_t id)
{
    uint64_t i; ul_str_t *str = &(uidx->pstr.str.a[id]); ma_ug_t *ug = uidx->init_ug;
    fprintf(stderr, "%.*s\tid::%lu\t", (int32_t)uidx->idx->nid.a[id].n, 
                                                        uidx->idx->nid.a[id].a, id);
    for (i = 0; i < str->cn; i++) {
        fprintf(stderr, "utg%.6d%c(%c)\t", (((uint32_t)str->a[i])>>1)+1, 
        "lc"[ug->u.a[(((uint32_t)str->a[i])>>1)].circ], "+-"[(((uint32_t)str->a[i])&1)]);
    }
    fprintf(stderr,"\n");
}

void ul_re_correct(ul_resolve_t *uidx, uint64_t n_r)
{   
    // print_ul_seq(uidx, 2062); print_ul_seq(uidx, 2093);
    uint64_t k, occ, n_circle; ///uidx->str_b.n_thread = 1;
    for (k = 0; k < n_r; k++) {
        kt_for(uidx->str_b.n_thread, worker_integer_correction, uidx, uidx->idx->n);
        occ = clean_ul_re_correct_buf(uidx, 0, &n_circle);
        fprintf(stderr, "+[M::%s::round->%lu] # corrected UL reads::%lu, # circle UL reads::%lu\n", 
                                                                                __func__, k, occ, n_circle);
        kt_for(uidx->str_b.n_thread, worker_integer_update, uidx, uidx->str_b.n_thread);
        occ = clean_ul_re_correct_buf(uidx, 1, &n_circle);
        fprintf(stderr, "-[M::%s::round->%lu] # corrected UL reads::%lu, # circle UL reads::%lu\n", 
                                                                                __func__, k, occ, n_circle);
        rebuid_idx(uidx);
        // print_ul_seq(uidx, 2062); print_ul_seq(uidx, 2093);
        // exit(1);
    }

    shrink_uls(uidx);
}

void append_utg_es(ul_resolve_t *uidx)
{
    uint64_t k; ul2ul_idx_t *u2o = &(uidx->uovl); ul2ul_item_t *it;
    for (k = 0; k < u2o->uln; k++) {
        it = get_ul_ovlp(u2o, k, 1);
        if(!it) continue;
        assert(uidx->pstr.str.a[k].cn > 1);
        integer_append_ug_ovlp(uidx, k, it, u2o);
    }

    kt_for(uidx->str_b.n_thread, worker_integert_clean, uidx, u2o->tot);///all ul + ug

    for (k = 0; k < u2o->tot; k++) {
        it = get_ul_ovlp(u2o, ulg_id(*u2o, k), ulg_type(*u2o, k));
        if(!it) continue;
        integer_normalize_ovlp(uidx, k, it, u2o);
    }

    kt_for(uidx->str_b.n_thread, worker_integert_clean, uidx, u2o->tot);///all ul + ug
}


void remove_integert_containment(ul_resolve_t *uidx, uint32_t keep_raw_utg)
{
    ul2ul_idx_t *u2o = &(uidx->uovl); ul2ul_item_t *o;
    uint64_t k, kn = keep_raw_utg?u2o->uln:u2o->tot; 
    for (k = 0; k < kn; k++) {
        o = get_ul_ovlp(u2o, ulg_id(*u2o, k), ulg_type(*u2o, k));
        if(!o) continue;
        integer_containment_purge(uidx, k, o, u2o, keep_raw_utg);
    }

    kt_for(uidx->str_b.n_thread, worker_integert_clean, uidx, u2o->tot);///all ul + ug
}

void print_integert_ovlp_stat(ul2ul_idx_t *ul2)
{
    uint64_t i, k, occ_r = 0, occ_o = 0; ul2ul_item_t *o;
    for (i = 0; i < ul2->uln; i++) {
        o = get_ul_ovlp(ul2, ulg_id(*ul2, i), ulg_type(*ul2, i));
        if((!o) || (o->is_del)) continue;
        occ_r++;
        for (k = 0; k < o->n; k++) {
            if(o->a[k].is_del || o->a[k].is_ct || (!ulg_type(*ul2, o->a[k].hid))) continue;
            occ_o++;
        }
    }
    fprintf(stderr, "[M::%s::] # UL reads::%lu, # UL ovlps::%lu\n", __func__, occ_r, occ_o);
}

asg_t *integer_sg_gen(ul_resolve_t *uidx, uint64_t min_ovlp)
{
    ul2ul_idx_t *ul2 = &(uidx->uovl); ma_ug_t *ug = uidx->l1_ug; asg_t *raw_g = ug->g; ///uc_block_t *xi; 
    uint64_t i, k, is_del, v, w, nv, z; ul2ul_item_t *o, *ow; int32_t r; asg_arc_t t, *p; asg_arc_t *av;
    // ul_str_t *str;
    asg_t *g = asg_init();
    for (i = 0; i < ul2->tot; i++) {
        is_del = 0;
        o = get_ul_ovlp(ul2, ulg_id(*ul2, i), ulg_type(*ul2, i));
        if((!o) || (o->is_del)) is_del = 1;
        asg_seq_set(g, i, ulg_len(*uidx, i), is_del);
        g->seq[i].c = 0;
    }
    // CALLOC(g->seq_vis, g->n_seq*2);

    for (i = 0; i < ul2->tot; ++i) {
        o = get_ul_ovlp(ul2, ulg_id(*ul2, i), ulg_type(*ul2, i));
        if((!o) || (o->is_del)) continue;
        for (k = 0; k < o->n; k++) {
            if(o->a[k].is_del || o->a[k].is_ct) continue;
            r = integer_hit2arc(&(o->a[k]), ulg_len(*uidx, i), ulg_len(*uidx, o->a[k].hid), 
            ulg_occ(*uidx, i), ulg_occ(*uidx, o->a[k].hid), i, o->a[k].hid, min_ovlp, &t);
            if (r >= 0) {
                p = asg_arc_pushp(g);
                *p = t;
            }
        }
        if(i >= ul2->uln) {///is a node of ug
            v = (i-ul2->uln)<<1;
            nv = asg_arc_n(raw_g, v); av = asg_arc_a(raw_g, v);
            for (z = 0; z < nv; z++) {
                if(av[z].del) continue;
                w = av[z].v + (ul2->uln<<1);
                ow = get_ul_ovlp(ul2, ulg_id(*ul2, (w>>1)), ulg_type(*ul2, (w>>1)));
                if((!ow) || (ow->is_del)) continue;
                p = asg_arc_pushp(g);
                *p = av[z]; p->ul += (((uint64_t)ul2->uln)<<33); p->v += (ul2->uln<<1);
            }

            v = ((i-ul2->uln)<<1)+1;
            nv = asg_arc_n(raw_g, v); av = asg_arc_a(raw_g, v);
            for (z = 0; z < nv; z++) {
                if(av[z].del) continue;
                w = av[z].v + (ul2->uln<<1);
                ow = get_ul_ovlp(ul2, ulg_id(*ul2, (w>>1)), ulg_type(*ul2, (w>>1)));
                if((!ow) || (ow->is_del)) continue;
                p = asg_arc_pushp(g);
                *p = av[z]; p->ul += (((uint64_t)ul2->uln)<<33); p->v += (ul2->uln<<1);
            }
        } 
        // else { ///is a node of ug
        //     str = &(uidx->pstr.str.a[i]);
        //     if(str->cn > 0) {
        //         xi = &(uidx->idx->a[i].bb.a[str->a[0]>>32]); v = ((uint32_t)str->a[0])^1;
        //         ow = get_ul_ovlp(ul2, ulg_id(*ul2, (v>>1)), ulg_type(*ul2, (v>>1)));
        //         if((!ow) || (ow->is_del)) {
        //             nv = asg_arc_n(raw_g, v); av = asg_arc_a(raw_g, v);
        //         }

        //         xi = &(uidx->idx->a[i].bb.a[str->a[str->cn-1]>>32]); v = ((uint32_t)str->a[str->cn-1]);
        //         ow = get_ul_ovlp(ul2, ulg_id(*ul2, (v>>1)), ulg_type(*ul2, (v>>1)));
        //         if((!ow) || (ow->is_del)) {

        //         }
        //     }
        // }
    }
    asg_cleanup(g);
    g->r_seq = g->n_seq;
    // fprintf(stderr, "[M::%s::] # ig nodes::%u, # ig archs::%u\n", __func__, g->n_seq, g->n_arc);
    return g;
}

inline uint64_t get_ul_occ(ul_resolve_t *uidx, uint64_t id)
{
    return uidx->uovl.cc.uc[id];
}

inline void get_iug_u_raw_occ(ul_resolve_t *uidx, uint32_t id, uint32_t *ul_occ, uint32_t *raw_ug_occ)
{
    if(ul_occ) *ul_occ = uidx->uovl.cc.raw_uc[id];
    if(raw_ug_occ) *raw_ug_occ = uidx->uovl.i_ug->u.a[id].n - uidx->uovl.cc.raw_uc[id];
}

inline asg_arc_t* get_specfic_edge(asg_t *g, uint32_t v, uint32_t w)
{
    asg_arc_t *av = asg_arc_a(g, v); uint32_t nv = asg_arc_n(g, v), k;
    for (k = 0; k < nv; k++) {
        if(av[k].del) continue;
        if(av[k].v == w) break;
    }

    if(k < nv) return (&av[k]);
    return NULL;    
}


void gen_ul_seq(ul_resolve_t *uidx, uint64_t iug_id, asgc8_v * res)
{
    ul2ul_idx_t *idx = &(uidx->uovl); ma_ug_t *raw = uidx->l1_ug; uint64_t k, m, s, e, ol, rev, Ns;
    uinfo_srt_warp_t *seq = &(idx->cc.iug_a[iug_id]); ma_utg_t *ru; asg_arc_t *z; 
    for (k = res->n = 0; k < seq->n; k++) {
        s = seq->a[k].s; e = seq->a[k].e; rev = seq->a[k].v&1; 
        ru = &(raw->u.a[seq->a[k].v>>1]); Ns = 0;
        if(k + 1 < seq->n) {
            z = get_specfic_edge(raw->g, seq->a[k].v, seq->a[k+1].v);
            if(z) {
                ol = z->ol;
                if(!rev) {
                    // assert(seq->a[k].e == ru->len);
                    e = (ru->len > ol)?(ru->len-ol):(0);
                } else {
                    // assert(seq->a[k].s == 0);
                    s = ol;
                }
            } else {
                Ns = 50;
            }
        }
        if(s < e) {
            kv_resize(char, (*res), res->n + e - s);
            retrieve_u_seq(NULL, res->a + res->n, ru, rev, (rev)?(ru->len-e):(s), e - s, NULL);
            res->n += e - s;
        }

        if(Ns) {
            kv_resize(char, (*res), res->n + Ns);
            for (m = 0; m < Ns; m++) res->a[res->n++] = 'N';
        }        
    }
    
    kv_push(char, *res, '\0');    
    // fprintf(stderr, "[M::%s::iug_id->%lu] # u->len::%u, # res->n::%u, # strlen(res->a)::%u\n", 
    // __func__, iug_id, idx->i_ug->u.a[iug_id].len, (uint32_t)res->n, (uint32_t)strlen(res->a));
}


void ma_integer_ug_print0(const ma_ug_t *ug, ul_resolve_t *uidx, int print_seq, const char* prefix, FILE *fp, uint32_t is_seq)
{
    uint32_t i, j, l, x; ma_utg_t *p, *s; ul2ul_idx_t *idx = &(uidx->uovl);
    char name[32]; uinfo_srt_warp_t *seq; asgc8_v t; kv_init(t); //uint64_t tot = 0;
    for (i = 0; i < ug->u.n; ++i) { // the Segment lines in GFA
        p = &ug->u.a[i];
        if(p->m == 0) continue;
        if(ug->g && ug->g->seq[i].del) continue;
        sprintf(name, "%s%.6d%c", prefix, i + 1, "lc"[p->circ]);
        if(is_seq) {
            gen_ul_seq(uidx, i, &t);
            fprintf(fp, "S\t%s\t%s\tLN:i:%d\trd:i:%lu\n", name, t.a, p->len, get_ul_occ(uidx, i));
        } else {
            fprintf(fp, "S\t%s\t*\tLN:i:%d\trd:i:%lu\n", name, p->len, get_ul_occ(uidx, i));
        }
        // tot += p->len;
        
        for (j = l = 0; j < p->n; j++) {
            if(p->a[j] != (uint64_t)-1) {
                x = p->a[j]>>33;
                if(ulg_type(uidx->uovl, x)) {//read
                    x = ulg_id(uidx->uovl, x);
                    fprintf(fp, "A\t%s\t%d\t%c\t%.*s\t%d\t%d\tid:i:%d\tHG:A:*\n", name, l, "+-"[p->a[j]>>32&1],
                    (int32_t)uidx->idx->nid.a[x].n, uidx->idx->nid.a[x].a, 0, uidx->idx->a[x].rlen, x);
                } else { ///node
                    x = ulg_id(uidx->uovl, x); s = &(uidx->init_ug->u.a[x]);
                    fprintf(fp, "A\t%s\t%d\t%c\tutg%.6d%c\t%d\t%d\tid:i:%d\tHG:A:*\n", name, l, "+-"[p->a[j]>>32&1],
                    x + 1, "lc"[s->circ], 0, s->len, x);
                }
            }
            else
            {
                fprintf(fp, "A\t%s\t%d\t*\t*\t*\t*\tid:i:*\tHG:A:*\n", name, l);
            }
            l += (uint32_t)p->a[j];
        }

        seq = &(idx->cc.iug_a[i]);
        for (j = 0; j < seq->n; j++) {
            x = seq->a[j].v>>1; s = &(uidx->init_ug->u.a[x]);
            fprintf(fp, "U\t%s\t%c\tutg%.6d%c\t%d\t%d\tHG:A:*\n", 
            name, "+-"[seq->a[j].v&1], x + 1, "lc"[s->circ], seq->a[j].s, seq->a[j].e);
        }
        
    }

    if(ug->g)
    {
        asg_arc_t* au = NULL;
        uint32_t nu, u, v;
        for (i = 0; i < ug->u.n; ++i) {
            if(ug->u.a[i].m == 0) continue;
            if(ug->u.a[i].circ)
            {
                fprintf(fp, "L\t%s%.6dc\t+\t%s%.6dc\t+\t%dM\tL1:i:%d\n", 
                prefix, i+1, prefix, i+1, 0, ug->u.a[i].len);
                fprintf(fp, "L\t%s%.6dc\t-\t%s%.6dc\t-\t%dM\tL1:i:%d\n", 
                prefix, i+1, prefix, i+1, 0, ug->u.a[i].len);
            } 
            u = i<<1;
            au = asg_arc_a(ug->g, u);
            nu = asg_arc_n(ug->g, u);
            for (j = 0; j < nu; j++)
            {
                if(au[j].del) continue;
                v = au[j].v;
                fprintf(fp, "L\t%s%.6d%c\t%c\t%s%.6d%c\t%c\t%dM\tL1:i:%d\tL2:i:%u\n", 
                prefix, (u>>1)+1, "lc"[ug->u.a[u>>1].circ], "+-"[u&1],
                prefix, (v>>1)+1, "lc"[ug->u.a[v>>1].circ], "+-"[v&1], au[j].ol, asg_arc_len(au[j]), au[j].ou);
            }


            u = (i<<1) + 1;
            au = asg_arc_a(ug->g, u);
            nu = asg_arc_n(ug->g, u);
            for (j = 0; j < nu; j++)
            {
                if(au[j].del) continue;
                v = au[j].v;
                fprintf(fp, "L\t%s%.6d%c\t%c\t%s%.6d%c\t%c\t%dM\tL1:i:%d\tL2:i:%u\n", 
                prefix, (u>>1)+1, "lc"[ug->u.a[u>>1].circ], "+-"[u&1],
                prefix, (v>>1)+1, "lc"[ug->u.a[v>>1].circ], "+-"[v&1], au[j].ol, asg_arc_len(au[j]), au[j].ou);
            }
        }
    }
    kv_destroy(t);
    // fprintf(stderr, "[M::%s::] tot::%lu\n", __func__, tot);
}


void gen_u2g_seq(ul_resolve_t *uidx, uint64_t iug_id, asgc8_v * res)
{
    ul2ul_idx_t *idx = &(uidx->uovl); ma_ug_t *raw = uidx->l1_ug; uint64_t k, m, s, e, ol, rev, Ns;
    uinfo_srt_warp_t *seq = &(idx->cc.iug_a[iug_id]); ma_utg_t *ru; asg_arc_t *z; 
    for (k = res->n = 0; k < seq->n; k++) {
        s = seq->a[k].s; e = seq->a[k].e; rev = seq->a[k].v&1; 
        ru = &(raw->u.a[seq->a[k].v>>1]); Ns = 0;
        if(k + 1 < seq->n) {
            z = get_specfic_edge(raw->g, seq->a[k].v, seq->a[k+1].v);
            if(z) {
                ol = z->ol;
                if(!rev) {
                    // assert(seq->a[k].e == ru->len);
                    e = (ru->len > ol)?(ru->len-ol):(0);
                } else {
                    // assert(seq->a[k].s == 0);
                    s = ol;
                }
            } else {
                Ns = 50;
            }
        }
        if(s < e) {
            kv_resize(char, (*res), res->n + e - s);
            retrieve_u_seq(NULL, res->a + res->n, ru, rev, (rev)?(ru->len-e):(s), e - s, NULL);
            res->n += e - s;
        }

        if(Ns) {
            kv_resize(char, (*res), res->n + Ns);
            for (m = 0; m < Ns; m++) res->a[res->n++] = 'N';
        }        
    }
    
    kv_push(char, *res, '\0');    
    // fprintf(stderr, "[M::%s::iug_id->%lu] # u->len::%u, # res->n::%u, # strlen(res->a)::%u\n", 
    // __func__, iug_id, idx->i_ug->u.a[iug_id].len, (uint32_t)res->n, (uint32_t)strlen(res->a));
}


void output_integer_graph(ul_resolve_t *uidx, ma_ug_t *iug, const char *nn, uint32_t is_seq)
{
    char* gfa_name = NULL; MALLOC(gfa_name, strlen(nn)+50);
	sprintf(gfa_name, "%s.integer.noseq.gfa", nn);
    FILE* fp = fopen(gfa_name, "w"); free(gfa_name);
	if (!fp) return;
    ma_integer_ug_print0(iug, uidx, 0, "itg", fp, is_seq);
    fclose(fp);
}


void print_raw_uls_seq(ul_resolve_t *uidx, const char *nn)
{
    char* gfa_name = NULL; MALLOC(gfa_name, strlen(nn)+70);
    sprintf(gfa_name, "%s.raw.integer.seq.log", nn);
    FILE* fp = fopen(gfa_name, "w"); free(gfa_name);
    if (!fp) return;
    ma_ug_t *ug = uidx->init_ug; all_ul_t *aln = uidx->idx;
    uint64_t id; uc_block_t *a = NULL; int64_t k, a_n;
    for (id = 0; id < aln->n; id++) {
        a = aln->a[id].bb.a; a_n = aln->a[id].bb.n; k = 0;
        if(a_n == 0) continue;
        fprintf(fp,"%.*s\tid::%lu\t", (int32_t)aln->nid.a[id].n, aln->nid.a[id].a, id);
        // for (k = 0; k < a_n && ug_occ_w(a[k].ts, a[k].te, &(ug->u.a[a[k].hid])) == 0; k++);
        for (; k < a_n; k++) {
            // if(ug_occ_w(a[k].ts, a[k].te, &(ug->u.a[a[k].hid])) == 0) break;
            fprintf(fp, "utg%.6d%c(%c)(n::%lu)\t", a[k].hid + 1, "lc"[ug->u.a[a[k].hid].circ], "+-"[a[k].rev], ug_occ_w(a[k].ts, a[k].te, &(ug->u.a[a[k].hid])));
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}

void print_raw_uls_aln(ul_resolve_t *uidx, const char *nn)
{
    char* gfa_name = NULL; MALLOC(gfa_name, strlen(nn)+70);
    sprintf(gfa_name, "%s.raw.aln.log", nn);
    FILE* fp = fopen(gfa_name, "w"); free(gfa_name);
    if (!fp) return;
    ma_ug_t *ug = uidx->init_ug; all_ul_t *aln = uidx->idx;
    uint64_t id; uc_block_t *a = NULL; int64_t k, a_n;
    for (id = 0; id < aln->n; id++) {
        a = aln->a[id].bb.a; a_n = aln->a[id].bb.n; k = 0;
        if(a_n == 0) continue;
        fprintf(fp,"%.*s\tid::%lu\tdd::%u\t", (int32_t)aln->nid.a[id].n, aln->nid.a[id].a, id, aln->a[id].dd);
        // for (k = 0; k < a_n && ug_occ_w(a[k].ts, a[k].te, &(ug->u.a[a[k].hid])) == 0; k++);
        for (; k < a_n; k++) {
            // if(ug_occ_w(a[k].ts, a[k].te, &(ug->u.a[a[k].hid])) == 0) break;
            fprintf(fp, "utg%.6d%c(%c)q::[%u, %u)\t", a[k].hid + 1, "lc"[ug->u.a[a[k].hid].circ], "+-"[a[k].rev], a[k].qs, a[k].qe);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}


void print_uls_seq(ul_resolve_t *uidx, const char *nn)
{
    char* gfa_name = NULL; MALLOC(gfa_name, strlen(nn)+50);
    sprintf(gfa_name, "%s.integer.seq.log", nn);
    FILE* fp = fopen(gfa_name, "w"); free(gfa_name);
    if (!fp) return;
    uint64_t k, i; all_ul_t *uls = uidx->idx; ul2ul_item_t *o;
    ul_str_idx_t *str_idx = &(uidx->pstr); ma_ug_t *ug = uidx->init_ug;
    for (k = 0; k < uls->n; k++) {
        if(str_idx->str.a[k].cn < 2) continue;
        o = get_ul_ovlp(&(uidx->uovl), ulg_id(uidx->uovl, k), ulg_type(uidx->uovl, k));
        if(!o) continue;
        fprintf(fp,"%.*s\tid::%lu\tdel::%u\t", (int32_t)uidx->idx->nid.a[k].n, uidx->idx->nid.a[k].a, k, o->is_del);
        for (i = 0; i < str_idx->str.a[k].cn; i++) {
            fprintf(fp, "utg%.6d%c(%c)\t", (((uint32_t)str_idx->str.a[k].a[i])>>1)+1, 
            "lc"[ug->u.a[(((uint32_t)str_idx->str.a[k].a[i])>>1)].circ], "+-"[(((uint32_t)str_idx->str.a[k].a[i])&1)]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}


void print_uls_ovs(ul_resolve_t *uidx, const char *nn)
{
    char* gfa_name = NULL; MALLOC(gfa_name, strlen(nn)+50);
    sprintf(gfa_name, "%s.integer.ovlp.log", nn);
    FILE* fp = fopen(gfa_name, "w"); free(gfa_name);
    if (!fp) return;
    uint64_t i, k, qv, tv; ma_ug_t *ug = uidx->l1_ug; ul2ul_idx_t *ul2 = &(uidx->uovl); ul2ul_item_t *o;
    for (i = 0; i < ul2->uln; ++i) {
        o = get_ul_ovlp(ul2, ulg_id(*ul2, i), ulg_type(*ul2, i));
        if(!o) continue;
        for (k = 0; k < o->n; k++) {
            qv = i; tv = o->a[k].hid;
            if(ulg_type(uidx->uovl, qv)) {//ul read
                qv = ulg_id(uidx->uovl, qv);
                fprintf(fp, "%.*s(id::%lu)\t%u\t%u(i::%u)\t%u(i::%u)\t%c\t", (int32_t)uidx->idx->nid.a[qv].n, uidx->idx->nid.a[qv].a, 
                    i, uidx->idx->a[qv].rlen, o->a[k].qs, o->a[k].qs_k, o->a[k].qe, o->a[k].qe_k, "+-"[o->a[k].is_rev]);
            } else {///ug node
                qv = ulg_id(uidx->uovl, qv);
                fprintf(fp, "utg%.6d%c(id::%lu)\t%u\t%u(i::%u)\t%u(i::%u)\t%c\t", (int32_t)qv + 1, "lc"[ug->u.a[qv].circ], 
                    i, ug->u.a[qv].len, o->a[k].qs, o->a[k].qs_k, o->a[k].qe, o->a[k].qe_k, "+-"[o->a[k].is_rev]);
            }

            if(ulg_type(uidx->uovl, tv)) {//ul read
                tv = ulg_id(uidx->uovl, tv);
                fprintf(fp, "%.*s(id::%u)\t%u\t%u(i::%u)\t%u(i::%u)\t", (int32_t)uidx->idx->nid.a[tv].n, uidx->idx->nid.a[tv].a, 
                    o->a[k].hid, uidx->idx->a[tv].rlen, o->a[k].ts, o->a[k].ts_k, o->a[k].te, o->a[k].te_k);

            } else {
                tv = ulg_id(uidx->uovl, tv);
                fprintf(fp, "utg%.6d%c(id::%u)\t%u\t%u(i::%u)\t%u(i::%u)\t\t", (int32_t)tv + 1, "lc"[ug->u.a[tv].circ], 
                o->a[k].hid, ug->u.a[tv].len, o->a[k].ts, o->a[k].ts_k, o->a[k].te, o->a[k].te_k);
            }

            fprintf(fp, "del::%u\tct::%u\n", o->a[k].is_del, o->a[k].is_ct);
        }
    }
    fclose(fp);
}



uint64_t gen_srt_cov_interval(uint64_t *b, uint64_t b_n)
{
    uint64_t i, m, start, end; int64_t dp, old_dp;
    radix_sort_srt64(b, b + b_n);
    for (i = m = 0, dp = 0, start = 0; i < b_n; ++i) {
        old_dp = dp;
        if (b[i]&1) --dp;
        else ++dp;

        if (old_dp < 1 && dp >= 1) {///old_dp < dp, b.a[j] is qs
            start = b[i]>>1;
        } else if (old_dp >= 1 && dp < 1){
            end = b[i]>>1;
            b[m] = start; b[m] <<= 32; b[m] += end; m++;
        } 
    }

    return m;
}

inline uint64_t get_remove_hifi_occ_back(ul_resolve_t *uidx, uint64_t thres, uint64_t *v_a, uint64_t v_n, asg64_v *buf)
{
    ul2ul_idx_t *idx = &(uidx->uovl);
    ma_ug_t *iug = idx->i_ug, *raw = uidx->l1_ug;
    uint64_t k, l, z, m, p, *raw_a, raw_n, raw_id, iug_id, iug_off, del_n, keep_ns[2], dup_raw, n_mask, pn, fn, full_del, is, ie; 
    uinfo_srt_warp_t *iu; uint64_t *s_a, s_n, ms, me; uinfo_srt_t *ps;
    for (k = del_n = fn = 0, buf->n = dup_raw = 0; k < v_n; k++) {
        iu = &(idx->cc.iug_a[v_a[k]]); ///integer unitigs
        for (z = 0; z < iu->n; z++) {
            raw_id = (iu->a[z].v>>1);
            raw_a = idx->cc.iug_b + idx->cc.iug_idx[raw_id];
            raw_n = idx->cc.iug_idx[raw_id+1] - idx->cc.iug_idx[raw_id];
            assert(raw_n > 0);
            for (m = keep_ns[0] = keep_ns[1] = 0, pn = buf->n, n_mask = 0; m < raw_n; m++) {
                iug_id = raw_a[m]>>32; iug_off = (uint32_t)raw_a[m];
                if(iug->g->seq[iug_id].del) continue;
                if(raw_a[m]&((uint64_t)(0x8000000000000000))) {
                    n_mask++;
                    continue;
                }
                if(iug_id == v_a[k] && iug_off == z) {///query itself
                    raw_a[m] |= ((uint64_t)(0x8000000000000000));
                    p = (uint32_t)-1; p <<= 32; p += idx->cc.iug_idx[raw_id] + m; 
                    kv_push(uint64_t, *buf, p);
                } else {
                    if(!(idx->cc.iug_a[iug_id].a[iug_off].v&1)) {
                        keep_ns[0] = MAX(idx->cc.iug_a[iug_id].a[iug_off].n, keep_ns[0]);
                    } else {
                        keep_ns[1] = MAX(idx->cc.iug_a[iug_id].a[iug_off].n, keep_ns[1]);
                    }                    
                }
            }
            assert(buf->n == pn + 1);
            if(keep_ns[0] + keep_ns[1] >= raw->u.a[raw_id].n) {
                keep_ns[0] = keep_ns[1] = raw->u.a[raw_id].n;
            }

            is = keep_ns[0]; ie = raw->u.a[raw_id].n - keep_ns[1];///[is, ie) -> uncovered coordinates
            if(ie > is) {
                if(!(iu->a[z].v&1)) {
                    ms = 0; me = iu->a[z].n;
                } else {
                    ms = raw->u.a[raw_id].n - iu->a[z].n; me = raw->u.a[raw_id].n;
                }
                if(MAX(is, ms) < MIN(ie, me)) {
                    del_n += MIN(ie, me) - MAX(is, ms);
                    buf->a[pn] <<= 32; buf->a[pn] >>= 32; buf->a[pn] |= (raw_id<<32); fn++;
                    if(n_mask) dup_raw = 1;
                }
            }
        }
    }

    if(del_n > thres && dup_raw) {
        radix_sort_srt64(buf->a, buf->a + buf->n); del_n = 0;
        for (l = 0, k = 1; k <= fn; k++) {
            if(k == fn || (buf->a[k]>>32) != (buf->a[l]>>32)) {
                raw_id = (buf->a[l]>>32); pn = buf->n;
                raw_a = idx->cc.iug_b + idx->cc.iug_idx[raw_id];
                raw_n = idx->cc.iug_idx[raw_id+1] - idx->cc.iug_idx[raw_id];
                assert(raw_n > 0);
                for (m = full_del = 0; m < raw_n; m++) {
                    iug_id = raw_a[m]>>32; iug_off = (uint32_t)raw_a[m];
                    if(iug->g->seq[iug_id].del) continue;
                    if(raw_a[m]&((uint64_t)(0x8000000000000000))) {
                        kv_push(uint64_t, *buf, (idx->cc.iug_a[iug_id].a[iug_off].s<<1));
                        kv_push(uint64_t, *buf, (idx->cc.iug_a[iug_id].a[iug_off].e<<1)+1);
                        if(idx->cc.iug_a[iug_id].a[iug_off].s == 0 && 
                                idx->cc.iug_a[iug_id].a[iug_off].e == raw->u.a[raw_id].len) {
                            full_del = 1; break;
                        }
                    }
                }
                
                s_a = buf->a + buf->n; s_n = buf->n - pn;
                assert(s_n > 0);
                
                if(s_n > 2 && full_del == 0) {
                    s_n = gen_srt_cov_interval(s_a, s_n);
                } else if(full_del) {///an interval has already cover the whole unitig
                    s_a[0] = raw->u.a[raw_id].len; s_n = 1;
                } else {///sn == 2;
                    s_a[0] >>= 1; s_a[0] <<= 32; s_a[0] += (s_a[1]>>1); s_n = 1;
                }


                for (m = full_del = 0; m < raw_n && full_del < s_n; m++) {
                    iug_id = raw_a[m]>>32; iug_off = (uint32_t)raw_a[m];
                    if(iug->g->seq[iug_id].del) continue;
                    if(raw_a[m]&((uint64_t)(0x8000000000000000))) continue;
                    ps = &(idx->cc.iug_a[iug_id].a[iug_off]);
                    for (z = 0; z < s_n; z++) {
                        if(s_a[z] == ((uint64_t)-1)) continue;
                        is = s_a[z]>>32; ie = (uint32_t)s_a[z];
                        ms = MAX(is, ps->s); me = MIN(ie, ps->e);
                        if(ms >= me) continue;
                        ///[ms, me) is unlikely to be contained in the [is, ie)
                        assert(ms == is || me == ie);
                        if(ms == is) is = me;
                        else if(me == ie) ie = ms;
                        if(is >= ie) {
                            s_a[z] = ((uint64_t)-1); full_del++;
                        } else {
                            s_a[z] = is; s_a[z] <<= 32; s_a[z] += ie;
                        } 
                    }
                }

                for (z = 0; z < s_n; z++) {
                    if(s_a[z] == ((uint64_t)-1)) continue;
                    is = s_a[z]>>32; ie = (uint32_t)s_a[z];
                    del_n += ug_occ_w(is, ie, &(raw->u.a[raw_id]));
                }
                if(del_n > thres) break;
                l = k; buf->n = pn;
            }
        }
    }

    for (k = 0; k < buf->n; k++) {
        assert(idx->cc.iug_b[(uint32_t)buf->a[k]]&((uint64_t)(0x8000000000000000)));
        idx->cc.iug_b[(uint32_t)buf->a[k]] <<= 1; idx->cc.iug_b[(uint32_t)buf->a[k]] >>= 1;
    }
    

    return ((del_n > thres)?0:1);
}


inline uint64_t get_remove_hifi_occ(ul_resolve_t *uidx, uint64_t thres, uint64_t *v_a, uint64_t v_n, asg64_v *buf, uint32_t *del_occ)
{
    ul2ul_idx_t *idx = &(uidx->uovl);
    ma_ug_t *iug = idx->i_ug, *raw = uidx->l1_ug;
    uint64_t k, l, z, m, p, *raw_a, raw_n, raw_id, iug_id, iug_off, del_n, keep_ns[2], del_ns[2], dup_raw, n_mask, pn, fn, is, ie; 
    uinfo_srt_warp_t *iu; uint64_t ms, me; 
    for (k = del_n = fn = 0, buf->n = dup_raw = 0; k < v_n; k++) {
        iu = &(idx->cc.iug_a[v_a[k]>>1]); ///integer unitigs
        for (z = 0; z < iu->n; z++) {
            raw_id = (iu->a[z].v>>1);
            raw_a = idx->cc.iug_b + idx->cc.iug_idx[raw_id];
            raw_n = idx->cc.iug_idx[raw_id+1] - idx->cc.iug_idx[raw_id];
            assert(raw_n > 0);
            for (m = keep_ns[0] = keep_ns[1] = 0, pn = buf->n, n_mask = 0; m < raw_n; m++) {
                iug_id = (raw_a[m]<<1)>>33; iug_off = (uint32_t)raw_a[m];
                // if(iug_id >= iug->g->n_seq){
                //     fprintf(stderr, ">>>[M::%s::] v_a[%lu]>>1::%lu, raw_id::%lu, raw_n::%lu, m::%lu, raw_a[m]::%lu, uidx->uovl.cc.iug_b[455]::%lu, iug_id::%lu, iug_off::%lu, iug->g->n_seq::%u\n", __func__, 
                //         k, v_a[k]>>1, raw_id, raw_n, m, raw_a[m], uidx->uovl.cc.iug_b[455], iug_id, iug_off, (uint32_t)iug->g->n_seq);
                // }
                if(iug->g->seq[iug_id].del) continue;
                if(raw_a[m]&((uint64_t)(0x8000000000000000))) {
                    n_mask++;
                    continue;
                }
                if(iug_id == (v_a[k]>>1) && iug_off == z) {///query itself
                    raw_a[m] |= ((uint64_t)(0x8000000000000000));
                    p = (uint32_t)-1; p <<= 32; p += idx->cc.iug_idx[raw_id] + m; 
                    kv_push(uint64_t, *buf, p);
                } else {
                    if(!(idx->cc.iug_a[iug_id].a[iug_off].v&1)) {
                        keep_ns[0] = MAX(idx->cc.iug_a[iug_id].a[iug_off].n, keep_ns[0]);
                    } else {
                        keep_ns[1] = MAX(idx->cc.iug_a[iug_id].a[iug_off].n, keep_ns[1]);
                    }                    
                }
            }
            assert(buf->n == pn + 1);
            if(keep_ns[0] + keep_ns[1] >= raw->u.a[raw_id].n) {
                keep_ns[0] = keep_ns[1] = raw->u.a[raw_id].n;
            }

            is = keep_ns[0]; ie = raw->u.a[raw_id].n - keep_ns[1];///[is, ie) -> uncovered coordinates
            if(ie > is) {
                if(!(iu->a[z].v&1)) {
                    ms = 0; me = iu->a[z].n;
                } else {
                    ms = raw->u.a[raw_id].n - iu->a[z].n; me = raw->u.a[raw_id].n;
                }
                if(MAX(is, ms) < MIN(ie, me)) {
                    del_n += MIN(ie, me) - MAX(is, ms);
                    buf->a[pn] <<= 32; buf->a[pn] >>= 32; buf->a[pn] |= (raw_id<<32); fn++;
                    if(n_mask) dup_raw = 1;
                }
            }
        }
    }

    if(del_n > thres && dup_raw) {
        radix_sort_srt64(buf->a, buf->a + buf->n); del_n = 0;
        for (l = 0, k = 1; k <= fn; k++) {
            if(k == fn || (buf->a[k]>>32) != (buf->a[l]>>32)) {
                raw_id = (buf->a[l]>>32); 
                raw_a = idx->cc.iug_b + idx->cc.iug_idx[raw_id];
                raw_n = idx->cc.iug_idx[raw_id+1] - idx->cc.iug_idx[raw_id];
                assert(raw_n > 0);
                for (m = del_ns[0] = del_ns[1] = keep_ns[0] = keep_ns[1] = 0; m < raw_n; m++) {
                    iug_id = (raw_a[m]<<1)>>33; iug_off = (uint32_t)raw_a[m];
                    if(iug->g->seq[iug_id].del) continue;
                    if(raw_a[m]&((uint64_t)(0x8000000000000000))) {
                        if(!(idx->cc.iug_a[iug_id].a[iug_off].v&1)) {
                            del_ns[0] = MAX(idx->cc.iug_a[iug_id].a[iug_off].n, del_ns[0]);
                        } else {
                            del_ns[1] = MAX(idx->cc.iug_a[iug_id].a[iug_off].n, del_ns[1]);
                        } 
                    } else {
                        if(!(idx->cc.iug_a[iug_id].a[iug_off].v&1)) {
                            keep_ns[0] = MAX(idx->cc.iug_a[iug_id].a[iug_off].n, keep_ns[0]);
                        } else {
                            keep_ns[1] = MAX(idx->cc.iug_a[iug_id].a[iug_off].n, keep_ns[1]);
                        }   
                    }
                }
                assert(del_ns[0] + del_ns[1] > 0);
                if(keep_ns[0] + keep_ns[1] >= raw->u.a[raw_id].n) {
                    keep_ns[0] = keep_ns[1] = raw->u.a[raw_id].n;
                }
                

                is = keep_ns[0]; ie = raw->u.a[raw_id].n - keep_ns[1];///[is, ie) -> uncovered coordinates
                if(ie > is) {
                    if(del_ns[0] + del_ns[1] >= raw->u.a[raw_id].n) {
                        del_n += ie - is;
                    } else {
                        ms = 0; me = del_ns[0];
                        if(ms < me && MAX(is, ms) < MIN(ie, me)) {
                            del_n += MIN(ie, me) - MAX(is, ms);
                        }

                        ms = raw->u.a[raw_id].n - del_ns[1]; me = raw->u.a[raw_id].n;
                        if(ms < me && MAX(is, ms) < MIN(ie, me)) {
                            del_n += MIN(ie, me) - MAX(is, ms);
                        }
                    }
                }
                if(del_n > thres) break;
                l = k; 
            }
        }
    }

    for (k = 0; k < buf->n; k++) {
        assert(idx->cc.iug_b[(uint32_t)buf->a[k]]&((uint64_t)(0x8000000000000000)));
        idx->cc.iug_b[(uint32_t)buf->a[k]] <<= 1; idx->cc.iug_b[(uint32_t)buf->a[k]] >>= 1;
    }
    // fprintf(stderr, "*[M::%s::] del_n::%lu\n", __func__, del_n);
    if(del_occ) *del_occ = del_n;
    return ((del_n > thres)?0:1);
}
#define bg_correct (1)
#define bg_wrong (2)
#define bg_ambiguous (3)
#define bg_unavailable ((uint32_t)-1)
inline uint32_t get_bg_flag(ul_resolve_t *uidx, uint32_t v, uint32_t w)
{
    asg_t *g = uidx->uovl.bg.bg; asg_arc_t *av; uint32_t k, nv; 
    if(g->seq[v>>1].del || g->seq[w>>1].del) return bg_unavailable;

    av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
    for (k = 0; k < nv; k++) {
        if(av[k].v == w) break;
    }
    if(k >= nv) return bg_unavailable;
    if(av[k].ou == 3) return bg_correct;
    if(av[k].ou == 2) return bg_wrong;
    return bg_ambiguous;
} 

inline void get_bridges(ul_resolve_t *uidx, uint64_t *v_a, uint64_t v_n, uint32_t *w_occ, uint32_t *am_occ)
{
    ul_bg_t *bg = &(uidx->uovl.bg); uint64_t k, w, am, l; uint32_t uv, uw, bv, bw; 
    for (k = w = am = 0, uv = uw = (uint32_t)-1; k < v_n; k++) {
        w += bg->w_n[v_a[k]>>1]; am += bg->a_n[v_a[k]>>1]; uw = v_a[k];
        if(uv != (uint32_t)-1) {
            bv = uidx->uovl.cc.iug_a[uv>>1].a[((uv&1)?(0):(uidx->uovl.cc.iug_a[uv>>1].n-1))].v; if(uv&1) bv ^= 1;
            bw = uidx->uovl.cc.iug_a[uw>>1].a[((uw&1)?(uidx->uovl.cc.iug_a[uw>>1].n-1):(0))].v; if(uw&1) bw ^= 1;
            if((ulg_type(uidx->uovl, (bv>>1))) && (ulg_type(uidx->uovl, (bw>>1)))) {
                l = get_bg_flag(uidx, bv, bw);
                if(l == bg_wrong) w++;
                if(l == bg_ambiguous) am++;
            }
        }
        uv = uw;
    }
    if(w_occ) *w_occ = w;
    if(am_occ) *am_occ = am;
}

static inline void ulg_seq_del(ma_ug_t *ug, uint32_t s)
{
	uint32_t k; asg_t *g = ug->g;
	g->seq[s].del = 1;
	for (k = 0; k < 2; ++k) {
		uint32_t i, v = s<<1 | k;
		uint32_t nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		for (i = 0; i < nv; ++i) {
			av[i].del = 1;
			asg_arc_del(g, av[i].v^1, v^1, 1);
		}
	}
    free(ug->u.a[s].a); memset(&(ug->u.a[s]), 0, sizeof(ug->u.a[s]));
}


uint32_t ulg_arc_cut_tips(ul_resolve_t *uidx, ma_ug_t *ug, uint32_t max_ext, uint32_t max_ext_hifi, uint32_t is_double_check, asg64_v *in, asg64_v *ib)
{
    asg64_v tx = {0,0,0}, tb = {0,0,0}, *b = NULL, *ub = NULL; asg_t *g = ug->g;
    uint32_t n_vtx = g->n_seq<<1, v, w, i, k, cnt = 0, nv, kv, pb, w_occ, a_occ, ul_occ, del_occ, is_telo;
    asg_arc_t *av = NULL; uint64_t lw;
    b = (in?(in):(&tx)); ub = (ib?(ib):(&tb));
    for (v = 0, b->n = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;

        av = asg_arc_a(g, v^1); nv = asg_arc_n(g, v^1);
        for (i = kv = 0; i < nv; i++) {
            if (av[i].del) continue;
            kv++; break;
        }

        if(kv) continue;
        // kv = ug->u.a[v>>1].n;/// get_ul_occ(uidx, v>>1);
        get_iug_u_raw_occ(uidx, v>>1, &ul_occ, NULL); kv = ul_occ;
        is_telo = 0; if(uidx->uovl.telo && uidx->uovl.telo[v>>1]) is_telo = 1;
        for (i = 0, w = v; i < max_ext; i++) {
            if(asg_end(g, w^1, &lw, NULL)!=0) break;
            w = (uint32_t)lw; 
            // kv += ug->u.a[w>>1].n; ///get_ul_occ(uidx, w>>1);
            get_iug_u_raw_occ(uidx, w>>1, &ul_occ, NULL); kv += ul_occ;
            if(uidx->uovl.telo && uidx->uovl.telo[w>>1]) is_telo = 1;
        }

        if((kv <= max_ext) && (!is_telo)) kv_push(uint64_t, *b, (((uint64_t)kv)<<32)|v);  
    }

    radix_sort_srt64(b->a, b->a + b->n);

    for (k = 0; k < b->n; k++) {
        v = (uint32_t)(b->a[k]);
        if (g->seq[v>>1].del) continue;

        av = asg_arc_a(g, v^1); nv = asg_arc_n(g, v^1);
        for (i = kv = 0; i < nv; i++) {
            if (av[i].del) continue;
            kv++; break;
        }

        if(kv) continue;

        pb = b->n; kv_push(uint64_t, *b, v); 
        // kv = ug->u.a[v>>1].n;/// get_ul_occ(uidx, v>>1);
        get_iug_u_raw_occ(uidx, v>>1, &ul_occ, NULL); kv = ul_occ;
        is_telo = 0; if(uidx->uovl.telo && uidx->uovl.telo[v>>1]) is_telo = 1;
        for (i = 0, w = v; i < max_ext; i++) {
            if(asg_end(g, w^1, &lw, NULL)!=0) break;
            w = (uint32_t)lw; kv_push(uint64_t, *b, w); 
            // kv += ug->u.a[w>>1].n; ///get_ul_occ(uidx, w>>1);
            get_iug_u_raw_occ(uidx, w>>1, &ul_occ, NULL); kv += ul_occ;
            if(uidx->uovl.telo && uidx->uovl.telo[w>>1]) is_telo = 1;
        }

        

        if((!is_telo) && (kv <= max_ext) && (get_remove_hifi_occ(uidx, max_ext_hifi, b->a + pb, b->n - pb, ub, &del_occ))) {
            // fprintf(stderr, "*[M::%s::] k::%u, v>>1::%u, v&1::%u, kv::%u, max_ext::%u, max_ext_hifi::%u\n\n", 
            // __func__, k, v>>1, v&1, kv, max_ext, max_ext_hifi);
            if(is_double_check) get_bridges(uidx, b->a + pb, b->n - pb, &w_occ, &a_occ);
            if(is_double_check == 0 || w_occ > 0 || a_occ > 0 || kv == 0 || del_occ == 0) {
                for (i = pb; i < b->n; i++) ulg_seq_del(ug, (b->a[i]>>1));
                cnt += b->n - pb;
            }
        }
        b->n = pb; 
    }

    // stats_sysm(g);
    if(!in) free(tx.a); if(!ib) free(tb.a);
    if (cnt > 0) asg_cleanup(g);

    return cnt;
}

uint32_t is_het_ulg_edge(ul_resolve_t *uidx, uint32_t uv, uint32_t uw)
{
    ma_ug_t *iug = uidx->uovl.i_ug; bubble_type *bub = uidx->bub; uint32_t qid, tid;
    qid = iug->u.a[uv>>1].a[((uv&1)?(0):(iug->u.a[uv>>1].n-1))]>>33;
    tid = iug->u.a[uw>>1].a[((uw&1)?(iug->u.a[uw>>1].n-1):(0))]>>33;

    if(!ulg_type(uidx->uovl, qid)) return (!IF_HOM(ulg_id(uidx->uovl, qid), *bub));
    if(!ulg_type(uidx->uovl, tid)) return (!IF_HOM(ulg_id(uidx->uovl, tid), *bub));

    if(uidx->uovl.item_idx[qid] == (uint32_t)-1) return (uint32_t)-1;
    ul2ul_item_t *o = &(uidx->uovl.a[uidx->uovl.item_idx[qid]]); 
    ul2ul_t *z; uint64_t k, n_hom, n_het; ul_str_t *str;
    for (k = 0, z = NULL; k < o->cn; k++) {
        if(o->a[k].hid != tid) continue;
        z = &(o->a[k]);
        break;
    }
    if(!z) return (uint32_t)-1;
    str = &(uidx->pstr.str.a[ulg_id(uidx->uovl, qid)]);

    for (k = z->qs_k, n_hom = n_het = 0; k < z->qe_k; k++) {
        if(IF_HOM((((uint32_t)str->a[k])>>1), *bub)) {
            n_hom++;
        } else {
            n_het++; break;
        }
    }
    
    if(n_het) return 1;
    if(n_hom) return 0;
    return (uint32_t)-1;
}


int32_t usg_topocut_aux(ul_resolve_t *uidx, ma_ug_t *ug, uint32_t v, int32_t max_ext, int32_t max_ext_hifi, asg64_v *b, asg64_v *ub)
{
    int32_t n_ext; asg_arc_t *av; uint32_t w = v, nv, i, kv, ul_occ, pn = b->n;
    for (n_ext = 0; n_ext < max_ext; v = w) {
        av = asg_arc_a(ug->g, v^1); nv = asg_arc_n(ug->g, v^1);
        for (i = kv = 0; i < nv && kv <= 1; i++) {
            if (av[i].del) continue;
            kv++; 
        }
        if(kv!=1) break;
        get_iug_u_raw_occ(uidx, v>>1, &ul_occ, NULL); n_ext += ul_occ; kv_push(uint64_t, *b, v);

        av = asg_arc_a(ug->g, v); nv = asg_arc_n(ug->g, v);
        for (i = kv = 0; i < nv && kv <= 1; i++) {
            if (av[i].del) continue;
            kv++; w = av[i].v;
        }
        if(kv!=1) break;
    }

    if(n_ext < max_ext) {
        if(get_remove_hifi_occ(uidx, max_ext_hifi, b->a + pn, b->n - pn, ub, NULL)) {
            b->n = pn;
            return 1;
        }
    }

    b->n = pn;
    return 0;
}

uint32_t ulg_arc_cut_length(ul_resolve_t *uidx, ma_ug_t *ug, int32_t max_ext, uint32_t max_ext_hifi, 
float len_rat, uint32_t is_trio, uint32_t topo_level, uint32_t hom_check, uint32_t *max_drop_len, asg64_v *in, asg64_v *ib)
{
    asg64_v tx = {0,0,0}, tb = {0,0,0}, *b = NULL, *ub = NULL; asg_t *g = ug->g;
    uint32_t v, w, i, k, kv, kw, nv, nw, cnt = 0, n_vtx = g->n_seq<<1, n_het, n_hom, ff, ol_max, mm_ol, to_del;
    asg_arc_t *av, *aw, *ve, *we, *vl_max, *wl_max;

    b = (in?(in):(&tx)); ub = (ib?(ib):(&tb));
    for (v = b->n = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;
        if(g->seq_vis[v] == 0) {
            av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
            if (nv < 2) continue;

            for (i = kv = 0; i < nv && kv < 2; ++i) {
                if(av[i].del) continue;
                kv++;
            }
            if(kv < 2) continue;
            if(hom_check) {
                for (i = n_het = n_hom = 0; i < nv; ++i) {
                    if(av[i].del) continue;
                    ff = is_het_ulg_edge(uidx, av[i].ul>>32, av[i].v); assert(ff != (uint32_t)-1);
                    if(ff) n_het++;
                    else n_hom++;
                    if(n_het > 0 && n_hom > 0) break;
                }
                if(n_het == 0 || n_hom == 0) continue;
            }

            for (i = 0; i < nv; ++i) {
                if(av[i].del) continue;
                kv_push(uint64_t, *b, (((uint64_t)av[i].ol)<<32) | ((uint64_t)(av-g->arc+i)));
            }
        }
    }

    radix_sort_srt64(b->a, b->a + b->n);
    for (k = 0; k < b->n; k++) {
        if(g->arc[(uint32_t)b->a[k]].del) continue;
        v = g->arc[(uint32_t)b->a[k]].ul>>32; w = g->arc[(uint32_t)b->a[k]].v^1;
        if(g->seq[v>>1].del || g->seq[w>>1].del) continue;
        nv = asg_arc_n(g, v); nw = asg_arc_n(g, w);
        av = asg_arc_a(g, v); aw = asg_arc_a(g, w);
        if(nv<=1 && nw <= 1) continue;
        if(hom_check) {
            ff = is_het_ulg_edge(uidx, g->arc[(uint32_t)b->a[k]].ul>>32, g->arc[(uint32_t)b->a[k]].v); 
            assert(ff != (uint32_t)-1);
            if(ff) continue;
        }

        ve = &(g->arc[(uint32_t)b->a[k]]);
        for (i = 0; i < nw; ++i) {
            if (aw[i].v == (v^1)) {
                we = &(aw[i]);
                break;
            }
        }
        mm_ol = MIN(ve->ol, we->ol); ///ve and we are hom edges

        for (i = kv = ol_max = 0, vl_max = NULL; i < nv; ++i) {
            if(av[i].del) continue;
            kv++;
            if(hom_check) { 
                ff = is_het_ulg_edge(uidx, av[i].ul>>32, av[i].v); assert(ff != (uint32_t)-1);
                if(!ff) continue;///vl_max must be a het edge
            }
            if(ol_max < av[i].ol) ol_max = av[i].ol, vl_max = &(av[i]);
        }
        if (kv < 1 || (!vl_max)) continue;
        if (kv >= 2) {
            if (mm_ol > ol_max*len_rat) continue;
        }

        for (i = kw = ol_max = 0, wl_max = NULL; i < nw; ++i) {
            if(aw[i].del) continue;
            kw++;
            if(hom_check) { 
                ff = is_het_ulg_edge(uidx, aw[i].ul>>32, aw[i].v); assert(ff != (uint32_t)-1);
                if(!ff) continue;///wl_max must be a het edge
            }
            if(ol_max < aw[i].ol) ol_max = aw[i].ol, wl_max = &(aw[i]);
        }
        if (kw < 1 || (!wl_max)) continue;
        if (kw >= 2) {
            if (mm_ol > ol_max*len_rat) continue;
        }

        if (kv <= 1 && kw <= 1) continue;

        to_del = 0;
        if(topo_level == 0) {
            to_del = 1;
        } else if(topo_level == 2) {
            if (kv > 1 && kw > 1) to_del = 1;
        } else {
            if (kv > 1 && kw > 1) {
                to_del = 1;
            } else if (kw == 1) {
                if (usg_topocut_aux(uidx, ug, w^1, max_ext, max_ext_hifi, b, ub)) to_del = 1;                    
            } else if (kv == 1) {
                if (usg_topocut_aux(uidx, ug, v^1, max_ext, max_ext_hifi, b, ub)) to_del = 1;   
            }
        }

        if (to_del) {
            ve->del = we->del = 1, ++cnt;
        }
    }


    if(!in) free(tx.a); if(!ib) free(tb.a);
    if (cnt > 0) asg_cleanup(g);

    return cnt;
}

uint32_t is_ul_edge(ul_resolve_t *uidx, uint32_t uv, uint32_t uw)
{
    ma_ug_t *iug = uidx->uovl.i_ug; uint32_t qid, tid;
    qid = iug->u.a[uv>>1].a[((uv&1)?(0):(iug->u.a[uv>>1].n-1))]>>33;
    tid = iug->u.a[uw>>1].a[((uw&1)?(iug->u.a[uw>>1].n-1):(0))]>>33;
    if((qid != tid) && (!ulg_type(uidx->uovl, qid)) && (!ulg_type(uidx->uovl, tid))) return 0;
    return 1;
}

uint32_t ulg_arc_cut_occ(ul_resolve_t *uidx, ma_ug_t *ug, int32_t max_ext, uint32_t max_ext_hifi, 
uint32_t is_trio, uint32_t topo_level, asg64_v *in, asg64_v *ib)
{
    asg_t *g = ug->g; asg64_v tx = {0,0,0}, tb = {0,0,0}, *b = NULL, *ub = NULL;
    uint32_t v, w, i, z, kv, kw, nv, nw, cnt = 0, n_vtx = g->n_seq<<1, to_del, n_ul, n_ug;
    asg_arc_t *av, *aw; b = (in?(in):(&tx)); ub = (ib?(ib):(&tb));

    for (v = b->n = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;
        av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
        if (nv < 2) continue;
        
        for (i = n_ul = n_ug = kv = 0; i < nv; ++i) {
            if(av[i].del) continue;
            if(is_ul_edge(uidx, av[i].ul>>32, av[i].v)) n_ul++;
            else n_ug++;
            kv++;
        }
        if(n_ul == 0 || n_ug == 0 || kv < 2) continue;

        for (i = 0; i < nv; ++i) {
            if(av[i].del) continue;
            if(is_ul_edge(uidx, av[i].ul>>32, av[i].v)) continue;
            w = av[i].v^1; if(g->seq[w>>1].del) continue;
            kw = get_arcs(g, w, NULL, 0);
            if (kv <= 1 && kw <= 1) continue;

            to_del = 0;
            if(topo_level == 0) {
                to_del = 1;
            } else if(topo_level == 2) {
                if (kv > 1 && kw > 1) to_del = 1;
            } else {
                if (kv > 1 && kw > 1) {
                    to_del = 1;
                } else if (kw == 1) {
                    if (usg_topocut_aux(uidx, ug, w^1, max_ext, max_ext_hifi, b, ub)) to_del = 1;                    
                } else if (kv == 1) {
                    if (usg_topocut_aux(uidx, ug, v^1, max_ext, max_ext_hifi, b, ub)) to_del = 1;   
                }
            }

            if (to_del) {
                ++cnt; av[i].del = 1;
                aw = asg_arc_a(g, w); nw = asg_arc_n(g, w);
                for (z = 0; z < nw; ++z) {
                    if (aw[z].v == (v^1)) {
                        aw[z].del = 1;
                        break;
                    }
                }
                assert(z < nw);
            }
        }
    }

    if(!in) free(tx.a); if(!ib) free(tb.a);
    if (cnt > 0) asg_cleanup(g);

    return cnt;
}


uint32_t get_ul_path_info(ul_resolve_t *uidx, ma_ug_t *ug, uint32_t s, uint32_t *e, uint32_t *occ, 
uint32_t *ul_cnt, uint32_t *bridge_w, uint32_t *bridge_am, asg64_v *b)
{
    uint32_t v = s, w = (uint32_t)-1, kv, kw, uv = (uint32_t)-1, uw = (uint32_t)-1, bv, bw, l; 
    ul_bg_t *bg = &(uidx->uovl.bg); 
    if(occ) (*occ) = 0; if(ul_cnt) (*ul_cnt) = 0; if(bridge_w) (*bridge_w) = 0; if(bridge_am) (*bridge_am) = 0;

    while (1) {
        if(occ) (*occ)++;
        kv = get_arcs(ug->g, v, &w, 1);        
        if(e) (*e) = v;

        ///if(b) kv_push(uint32_t, b->b, v>>1);
        if(b) kv_push(uint64_t, *b, v);
        if(ul_cnt) (*ul_cnt) += get_ul_occ(uidx, v>>1);
        if(bridge_w || bridge_am) {
            uw = v; 
            if(bridge_w) (*bridge_w) += bg->w_n[v>>1];
            if(bridge_am) (*bridge_am) += bg->a_n[v>>1];
            if(uv != (uint32_t)-1) {
                bv = uidx->uovl.cc.iug_a[uv>>1].a[((uv&1)?(0):(uidx->uovl.cc.iug_a[uv>>1].n-1))].v; if(uv&1) bv ^= 1;
                bw = uidx->uovl.cc.iug_a[uw>>1].a[((uw&1)?(uidx->uovl.cc.iug_a[uw>>1].n-1):(0))].v; if(uw&1) bw ^= 1;
                if((ulg_type(uidx->uovl, (bv>>1))) && (ulg_type(uidx->uovl, (bw>>1)))) {
                    l = get_bg_flag(uidx, bv, bw);
                    if(l == bg_wrong && bridge_w) (*bridge_w)++;
                    if(l == bg_ambiguous && bridge_am) (*bridge_am)++;
                }
            }
            uv = uw;
        }

        if(kv == 0) return END_TIPS;
        if(kv == 2) return TWO_OUTPUT;
        if(kv > 2) return MUL_OUTPUT;
        w = ug->g->arc[w].v;
        ///up to here, kv=1
        ///kw must >= 1
        kw = get_arcs(ug->g, w^1, NULL, 0);
        v = w;

        if(kw == 2) return TWO_INPUT;
        if(kw > 2) return MUL_INPUT;
        if(v == s) return LOOP;
    }

    return LONG_TIPS;  
}

#define ul_path_w(ul_cnt, bridge_w, bridge_am) (((bridge_w)+(bridge_w))==0?((uint32_t)-1):((ul_cnt)/((bridge_w)+(bridge_w))))


uint32_t usg_bridge_topocut_aux(ul_resolve_t *uidx, ma_ug_t *ug, uint32_t v, uint32_t max_ext, uint32_t max_ext_hifi, uint32_t topo_level, uint32_t double_check, asg64_v *b, asg64_v *ub)
{
    uint32_t k, z, bn = b->n, w, bridge_w = 1, bridge_am = 1, raw_ul, ul, del_occ, is_del = 0, kk;
    asg_arc_t *av; uint32_t nv;
    get_ul_path_info(uidx, ug, v, &w, NULL, NULL, (double_check?(&bridge_w):(NULL)), 
                                                                    (double_check?(&bridge_am):(NULL)), b);
    for (k = bn, raw_ul = 0; k < b->n; k++) {
        get_iug_u_raw_occ(uidx, b->a[k]>>1, &ul, NULL); raw_ul += ul;
    }
    
    if(raw_ul <= max_ext && get_remove_hifi_occ(uidx, max_ext_hifi, b->a + bn, b->n - bn, ub, &del_occ)) {
        b->n = bn;
        if(topo_level == 0) {
            is_del = 1;
        } else if(double_check == 0 || bridge_w > 0 || bridge_am > 0 || raw_ul == 0 || del_occ == 0) {
            av = asg_arc_a(ug->g, v^1); nv = asg_arc_n(ug->g, v^1);
            for (z = 0; z < nv; z++) {
                if(av[z].del) continue;
                kk = get_arcs(ug->g, av[z].v^1, NULL, 0);
                if(topo_level == 2 && kk > 1) continue;
                ///==1, might be a tip
                if (topo_level == 1 && usg_topocut_aux(uidx, ug, av[z].v, max_ext, max_ext_hifi, b, ub)) continue;
                break;
            }

            if(z >= nv) {
                av = asg_arc_a(ug->g, w); nv = asg_arc_n(ug->g, w);
                for (z = 0; z < nv; z++) {
                    if(av[z].del) continue;
                    kk = get_arcs(ug->g, av[z].v^1, NULL, 0);
                    if(topo_level == 2 && kk > 1) continue;
                    ///==1, might be a tip
                    if (topo_level == 1 && usg_topocut_aux(uidx, ug, av[z].v, max_ext, max_ext_hifi, b, ub)) continue;
                    break;
                }

                if(z >= nv) is_del = 1;
            }
        }
    }
    b->n = bn;

    return is_del;
}

uint32_t ulg_arc_cut_bridge(ul_resolve_t *uidx, ma_ug_t *ug, int32_t max_ext, uint32_t max_ext_hifi, 
float len_rat, uint32_t is_trio, uint32_t topo_level, uint32_t *max_drop_len, asg64_v *in, asg64_v *ib)
{
    asg64_v tx = {0,0,0}, tb = {0,0,0}, *b = NULL, *ub = NULL; asg_t *g = ug->g;
    uint32_t v, w, i, k, kv, nv, cnt = 0, n_vtx = g->n_seq<<1, ff, ol_max, ul_max, mm_ol;
    asg_arc_t *av, *ve, *vl_max; uint32_t ul_cnt, bridge_w, bridge_am, pb;
    // fprintf(stderr, "+++[M::%s::] idx->cc.iug_b[455]::%lu\n", __func__, uidx->uovl.cc.iug_b[455]);
    b = (in?(in):(&tx)); ub = (ib?(ib):(&tb));
    for (v = b->n = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;
        if(g->seq_vis[v] == 0) {
            av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
            if (nv < 2) continue;

            for (i = kv = 0; i < nv && kv < 2; ++i) {
                if(av[i].del) continue;
                kv++;
            }
            if(kv < 2) continue;

            for (i = 0; i < nv; ++i) {
                if(av[i].del) continue;
                get_ul_path_info(uidx, ug, av[i].v, NULL, NULL, &ul_cnt, &bridge_w, &bridge_am, NULL);
                ff = ul_path_w(ul_cnt, bridge_w, bridge_am);
                // if((v>>1) == 5) {
                //     fprintf(stderr, "+[M::%s::] v>>1::%u, v&1::%u, w>>1::%u, w&1::%u, ul_cnt::%u, bridge_w::%u, bridge_am::%u, ff::%u\n", __func__, 
                //     v>>1, v&1, av[i].v>>1, av[i].v&1, ul_cnt, bridge_w, bridge_am, ff);
                // }
                kv_push(uint64_t, *b, (((uint64_t)ff)<<32) | ((uint64_t)(av-g->arc+i)));
            }
        }
    }

    // fprintf(stderr, "---[M::%s::] idx->cc.iug_b[455]::%lu\n", __func__, uidx->uovl.cc.iug_b[455]);
    radix_sort_srt64(b->a, b->a + b->n);
    for (k = 0; k < b->n; k++) {
        if(g->arc[(uint32_t)b->a[k]].del) continue;
        v = g->arc[(uint32_t)b->a[k]].ul>>32; w = g->arc[(uint32_t)b->a[k]].v^1;
        if(g->seq[v>>1].del || g->seq[w>>1].del) continue;
        nv = asg_arc_n(g, v); av = asg_arc_a(g, v); 
        if(nv<=1 && asg_arc_n(g, w) <= 1) continue;
        if(get_arcs(ug->g, w, NULL, 0) != 1) continue;

        ve = &(g->arc[(uint32_t)b->a[k]]);
        get_ul_path_info(uidx, ug, ve->v, NULL, NULL, &ul_cnt, &bridge_w, &bridge_am, NULL);
        ff = ul_path_w(ul_cnt, bridge_w, bridge_am);
        if(ff == (uint32_t)-1) continue;
        mm_ol = ff; 

        for (i = kv = ol_max = ul_max = 0, vl_max = NULL; i < nv; ++i) {
            if(av[i].del) continue;
            kv++;
            get_ul_path_info(uidx, ug, av[i].v, NULL, NULL, &ul_cnt, &bridge_w, &bridge_am, NULL);
            ff = ul_path_w(ul_cnt, bridge_w, bridge_am);
            if((ol_max < ff) || (ol_max == ff && ul_max < ul_cnt)) {
                ol_max = ff; ul_max = ul_cnt; vl_max = &(av[i]);
            } 
        }

        // if((v>>1) == 5) {
        //     fprintf(stderr, "***[M::%s::] v>>1::%u, v&1::%u, w>>1::%u, w&1::%u, ul_cnt::%u, bridge_w::%u, bridge_am::%u, ff::%u, mm_ol::%u, ol_max::%u, kv::%u\n", __func__, 
        //     v>>1, v&1, w>>1, w&1, ul_cnt, bridge_w, bridge_am, ff, mm_ol, ol_max, kv);
        // }
        
        if (kv <= 1 || (!vl_max)) continue;
        if (kv >= 2) {
            if (mm_ol > ol_max*len_rat) continue;
        }

        // if((v>>1) == 5) {
        //     fprintf(stderr, "###[M::%s::] v>>1::%u, v&1::%u, w>>1::%u, w&1::%u, ul_cnt::%u, bridge_w::%u, bridge_am::%u, ff::%u, mm_ol::%u, ol_max::%u\n", __func__, 
        //     v>>1, v&1, w>>1, w&1, ul_cnt, bridge_w, bridge_am, ff, mm_ol, ol_max);
        // }

        
        if(usg_bridge_topocut_aux(uidx, ug, ve->v, max_ext, max_ext_hifi, topo_level, 1, b, ub)) {
            pb = b->n; 
            get_ul_path_info(uidx, ug, ve->v, NULL, NULL, NULL, NULL, NULL, b);
            for (i = pb; i < b->n; i++) ulg_seq_del(ug, (b->a[i]>>1));
            cnt += b->n - pb; b->n = pb;
        }
    }


    if(!in) free(tx.a); if(!ib) free(tb.a);
    if (cnt > 0) asg_cleanup(g);
    return cnt;
}



ul2ul_t* get_ul_spec_ovlp(ul2ul_idx_t *z, uint64_t qid, uint64_t tid)
{
    // uint64_t x = (is_qul?(qid):(qid+z->uln));
    if(z->item_idx[qid] == (uint32_t)-1) return NULL;
    ul2ul_item_t *o = &(z->a[z->item_idx[qid]]); uint64_t k;
    for (k = 0; k < o->cn; k++) {
        if(o->a[k].hid != tid) continue;
        return &(o->a[k]);
    }
    return NULL;
}

uint64_t gen_ug_integer_seq_on_fly(ul_resolve_t *uidx, uint64_t *u_a, uint64_t u_n, asg64_v *res)
{
    ul2ul_idx_t *idx = &(uidx->uovl); ma_ug_t *raw = uidx->l1_ug;
    uint64_t k, rn = res->n, nn, u_rev, uv, uw, bv, bw, is_bv_ul, is_bw_ul, x; int64_t z, zn;
    ma_ug_t *iug = idx->i_ug; uinfo_srt_warp_t *seq; ul2ul_t *m;
    for (k = 0; k < u_n; k++) {
        seq = &(idx->cc.iug_a[u_a[k]>>1]); u_rev = u_a[k]&1;
        // if(u_n == 6 && u_a[0] == 37683 && u_a[5] == 45894) {
        //     fprintf(stderr, "\n[M::%s::] k::%lu, u_a[k]>>1::%lu, u_a[k]&1::%lu\n", __func__, k, u_a[k]>>1, u_a[k]&1);
        //     for (z = 0, zn = seq->n; z < zn; z++) {
        //         fprintf(stderr, "[M::%s::] v>>1::%u, v&1::%u\n", __func__, seq->a[z].v>>1, seq->a[z].v&1);
        //     }
        // }
        if(k > 0) {
            uv = u_a[k-1]; uw = u_a[k];
            bv = iug->u.a[uv>>1].a[((uv&1)?(0):(iug->u.a[uv>>1].n-1))]>>32; if(uv&1) bv ^= 1; //pre
            bw = iug->u.a[uw>>1].a[((uw&1)?(iug->u.a[uw>>1].n-1):(0))]>>32; if(uw&1) bw ^= 1; //current
            is_bv_ul = ulg_type((*idx), (bv>>1)); is_bw_ul = ulg_type((*idx), (bw>>1)); 
            if(is_bw_ul) {
                m = get_ul_spec_ovlp(idx, bv>>1, bw>>1);
                assert(m && (!m->is_del));
                nn = res->n + seq->n - (m->te_k - m->ts_k); 
                kv_resize(uint64_t, *res, nn); 
                nn = m->te_k - m->ts_k; ///skipped length
                zn = seq->n; 
                // fprintf(stderr, "[M::%s::] bv>>1::%lu, bv&1::%lu, bw>>1::%lu, bw&1::%lu, m->qs_k::%u, m->qe_k::%u, m->ts_k::%u, m->te_k::%u\n", __func__, 
                // bv>>1, bv&1, bw>>1, bw&1, m->qs_k, m->qe_k, m->ts_k, m->te_k);
                if(!u_rev) {
                    ///debug
                    for (z = 0; z < ((int64_t)nn); z++) {
                        // if(res->a[res->n-nn+z] != seq->a[z].v) {
                        //     fprintf(stderr, "****[M::%s::] k::%lu, z::%ld, nn::%lu, res->n::%u, res->v>>1::%lu, res->v&1::%lu, seq->a[z].v>>1::%u, seq->a[z].v&1::%u\n", 
                        //     __func__, k, z, nn, (uint32_t)res->n, res->a[res->n-nn+z]>>1, res->a[res->n-nn+z]&1, seq->a[z].v>>1, seq->a[z].v&1);
                        // }
                        assert(res->a[res->n-nn+z] == seq->a[z].v);
                    }

                    for (z = nn; z < zn; z++) {
                        res->a[res->n++] = seq->a[z].v;
                    }
                } else {
                    ///debug
                    for (z = zn - 1; z >= zn - ((int64_t)nn); z--) {
                        // fprintf(stderr, "[M::%s::z->%ld::idx->%ld] nn::%lu, zn::%ld, res->n::%ld, res>>1::%lu, res&1::%lu, seq>>1::%u, seq&1::%u\n", 
                        // __func__, z, (int64_t)(res->n-nn+(zn-1-z)), nn, zn, (int64_t)(res->n), res->a[res->n-nn+(zn-1-z)]>>1, res->a[res->n-nn+(zn-1-z)]&1, seq->a[z].v>>1, seq->a[z].v&1);
                        assert(res->a[res->n-nn+(zn-1-z)] == (seq->a[z].v^1));
                    }

                    for (z = zn - 1 - ((int64_t)nn); z >= 0; z--) {
                        res->a[res->n++] = seq->a[z].v^1;
                    }
                }
            } else {
                x = ulg_id(uidx->uovl, (bw>>1)); x <<= 1; x += (bw&1);
                if(!is_bv_ul) {///if previous ul is also a raw utg node
                    assert(get_specfic_edge(raw->g, res->a[res->n-1], x));
                } else {
                    // assert(((res->a[res->n-1].v) == x));
                    if(((res->a[res->n-1]) == x)) {
                        res->n--;
                    } else {
                        assert(get_specfic_edge(raw->g, res->a[res->n-1], x));
                    }
                }
                nn = res->n + seq->n; zn = seq->n;
                kv_resize(uint64_t, *res, nn);
                if(!u_rev) {
                    for (z = 0; z < zn; z++) {
                        res->a[res->n++] = seq->a[z].v;
                    }
                } else {
                    for (z = zn - 1; z >= 0; z--) {
                        res->a[res->n++] = seq->a[z].v^1;
                    }
                }
            }
        } else {
            nn = res->n + seq->n; zn = seq->n;
            kv_resize(uint64_t, *res, nn);
            if(!u_rev) {
                for (z = 0; z < zn; z++) {
                    // fprintf(stderr, ">+<[M::%s::res->n->%ld] a>>1::%u, a&1::%u\n", 
                    //     __func__, (int64_t)res->n, seq->a[z].v>>1, seq->a[z].v&1);
                    res->a[res->n++] = seq->a[z].v;
                }
            } else {
                for (z = zn - 1; z >= 0; z--) {
                    // fprintf(stderr, ">-<[M::%s::res->n->%ld] a>>1::%u, a&1::%u\n", 
                    //     __func__, (int64_t)res->n, seq->a[z].v>>1, (seq->a[z].v^1)&1);
                    res->a[res->n++] = seq->a[z].v^1;
                }
            }
        }
    }
    
    return res->n - rn;
}

uint32_t get_integer_seq_ovlps(ul_resolve_t *uidx, uint64_t *p_a, int64_t p_n, int64_t match_bound, uint64_t skip_hom, asg64_v *res, uint64_t *r_w)
{
    ul_str_idx_t *str_idx = &(uidx->pstr); uint32_t v = p_a[match_bound]; uc_block_t *xi;
    uint64_t *hid_a, hid_n, z, vz, ps, pe, ww[2], sw[2], occ; ul_str_t *str; int64_t s_n, s, p;
    bubble_type *bub = uidx->bub;

    hid_a = str_idx->occ.a + str_idx->idx.a[v>>1];
    hid_n = str_idx->idx.a[(v>>1)+1] - str_idx->idx.a[v>>1];
    for (z = ww[0] = ww[1] = (*r_w) = occ = 0; z < hid_n; z++) {
        str = &(str_idx->str.a[hid_a[z]>>32]); s_n = str->cn;
        if(s_n < 2) continue;
        vz = (uint32_t)(str->a[(uint32_t)hid_a[z]]);
        assert((v>>1) == (vz>>1)); ps = pe = (uint64_t)-1;
        sw[0] = sw[1] = 0;
        if(v == vz) {
            s = ((uint32_t)hid_a[z]) + 1; p = match_bound + 1; 
            for (; (s < s_n) && (p < p_n) && ((uint32_t)(str->a[s]) == p_a[p]); s++, p++) {
                xi = &(uidx->idx->a[hid_a[z]>>32].bb.a[str->a[s]>>32]);
                assert(((xi->hid<<1)+xi->rev)==((uint32_t)str->a[s]));
                if(skip_hom && IF_HOM(xi->hid, *bub)) continue;
                sw[1] += ug_occ_w(xi->ts, xi->te, &(uidx->l1_ug->u.a[xi->hid]));
            }
            if(s < s_n && p < p_n) continue;
            if(p <= match_bound + 1) continue;///bridging the two sides of the breakpoint
            if(sw[1] == 0) continue;
            pe = p;

            s = (uint32_t)hid_a[z]; p = match_bound; 
            for (; (s >= 0) && (p >= 0) && ((uint32_t)(str->a[s]) == p_a[p]); s--, p--) {
                xi = &(uidx->idx->a[hid_a[z]>>32].bb.a[str->a[s]>>32]);
                assert(((xi->hid<<1)+xi->rev)==((uint32_t)str->a[s]));
                if(skip_hom && IF_HOM(xi->hid, *bub)) continue;
                sw[0] += ug_occ_w(xi->ts, xi->te, &(uidx->l1_ug->u.a[xi->hid]));
            }
            assert(p < match_bound);
            if(s >= 0 && p >= 0) continue;
            if(sw[0] == 0) continue;
            ps = p + 1;
        } else {
            s = ((int32_t)((uint32_t)hid_a[z]))-1; p = match_bound + 1; 
            for (; (s >= 0) && (p < p_n) && ((uint32_t)(str->a[s]) == (p_a[p]^1)); s--, p++) {
                xi = &(uidx->idx->a[hid_a[z]>>32].bb.a[str->a[s]>>32]);
                assert(((xi->hid<<1)+xi->rev)==((uint32_t)str->a[s]));
                if(skip_hom && IF_HOM(xi->hid, *bub)) continue;
                sw[1] += ug_occ_w(xi->ts, xi->te, &(uidx->l1_ug->u.a[xi->hid]));
            }
            if(s >= 0 && p < p_n) continue;
            if(p <= match_bound + 1) continue;///bridging the two sides of the breakpoint
            if(sw[1] == 0) continue;
            pe = p;

            s = (uint32_t)hid_a[z]; p = match_bound; 
            for (; (s < s_n) && (p >= 0) && ((uint32_t)(str->a[s]) == (p_a[p]^1)); s++, p--) {
                xi = &(uidx->idx->a[hid_a[z]>>32].bb.a[str->a[s]>>32]);
                assert(((xi->hid<<1)+xi->rev)==((uint32_t)str->a[s]));
                if(skip_hom && IF_HOM(xi->hid, *bub)) continue;
                sw[0] += ug_occ_w(xi->ts, xi->te, &(uidx->l1_ug->u.a[xi->hid]));
            }
            assert(p < match_bound);
            if(s < s_n && p >= 0) continue;
            if(sw[0] == 0) continue;
            ps = p + 1;
        }
        if(res) kv_push(uint64_t, *res, ((ps<<32)|(pe)));
        occ++;
        ww[0] += sw[0]; ww[1] += sw[1]; 
        // if(match_bound == 15 && p_n == 36) {
        //     fprintf(stderr, "[M::%s::] z::%lu, ulid::%lu, ps::%lu, pe::%lu\n", __func__, z, hid_a[z]>>32, ps, pe);
        // }
    }

    (*r_w) = MIN(ww[0], ww[1]);
    return occ;
}

uint32_t usg_misjoin_topocut_aux(ul_resolve_t *uidx, ma_ug_t *ug, uint32_t v, uint32_t ref_v, uint32_t beg_v, 
asg64_v *b, asg64_v *ub)
{
    #define cut_rate 0.49999
    uint32_t k, bn = b->n, is_del = 0, kk, ks, ke;
    uint64_t n_ref, n_v, *a_ref, *a_v, n_min, n_ref_ov, n_v_ov, *a_ref_ov, *a_v_ov, l_occ, r_occ, ref_occ, v_occ, ww[2]; 
    ub->n = 0; kv_push(uint64_t, *ub, beg_v);
    get_ul_path_info(uidx, ug, ref_v, NULL, NULL, NULL, NULL, NULL, ub); 
    // fprintf(stderr, "+++[M::%s::] start..., v>>1::%u, v&1::%u, ref_v>>1::%u, ref_v&1::%u, beg_v>>1::%u, beg_v&1::%u\n", __func__, v>>1, v&1, ref_v>>1, ref_v&1, beg_v>>1, beg_v&1);
    n_ref = gen_ug_integer_seq_on_fly(uidx, ub->a, ub->n, b);
    // fprintf(stderr, "+++[M::%s::] done...\n", __func__);

    ub->n = 0; kv_push(uint64_t, *ub, beg_v);
    get_ul_path_info(uidx, ug, v, NULL, NULL, NULL, NULL, NULL, ub);
    // fprintf(stderr, "---[M::%s::] start...\n", __func__);
    n_v = gen_ug_integer_seq_on_fly(uidx, ub->a, ub->n, b);
    // fprintf(stderr, "---[M::%s::] done...\n", __func__);
    
    a_ref = b->a + bn; a_v = b->a + bn + n_ref; n_min = MIN(n_ref, n_v); ub->n = 0;
    for (k = 0; k < n_min && a_ref[k] == a_v[k]; k++); 

    // if((beg_v>>1) == 5958) {
    //     fprintf(stderr, "[M::%s::] beg_v>>1::%u, beg_v&1::%u, v>>1::%u, v&1::%u, n_v::%lu, ref_v>>1::%u, ref_v&1::%u, n_ref::%lu, # prefix::%u\n", 
    //     __func__, beg_v>>1, beg_v&1, v>>1, v&1, n_v, ref_v>>1, ref_v&1, n_ref, k);
    //     for (kk = 0; kk < n_v; kk++) {
    //         fprintf(stderr, "+[M::%s::] a_v[%u]>>1::%lu, a_v[%u]&1::%lu\n", 
    //                                                 __func__, kk, a_v[kk]>>1, kk, a_v[kk]&1);
    //     }

    //     for (kk = 0; kk < n_ref; kk++) {
    //         fprintf(stderr, "+[M::%s::] a_ref[%u]>>1::%lu, a_ref[%u]&1::%lu\n", 
    //                                                 __func__, kk, a_ref[kk]>>1, kk, a_ref[kk]&1);
    //     }
    // }


    if(k > 0) {
        n_ref_ov = get_integer_seq_ovlps(uidx, a_ref, n_ref, k - 1, 1, ub, &ww[0]);
        n_v_ov = get_integer_seq_ovlps(uidx, a_v, n_v, k - 1, 1, ub, &ww[1]);
        if(n_v_ov <= (n_ref_ov*cut_rate)) {
            is_del = 1;
        } else {
            a_ref_ov = ub->a; a_v_ov = ub->a + n_ref_ov; kk = k;
            
            for (k = ref_occ = 0; k < n_ref_ov; k++) {
                ks = a_ref_ov[k]>>32; ke = (uint32_t)a_ref_ov[k];
                l_occ = kk - ks; r_occ = ke - kk;
                ref_occ += MIN(l_occ, r_occ);
            }

            for (k = v_occ = 0; k < n_v_ov; k++) {
                ks = a_v_ov[k]>>32; ke = (uint32_t)a_v_ov[k];
                l_occ = kk - ks; r_occ = ke - kk;
                v_occ += MIN(l_occ, r_occ);
            }

            if(v_occ <= (ref_occ*cut_rate)) is_del = 1;
        }
    }
    b->n = bn;
    return is_del;
}


uint32_t ulg_arc_cut_misjoin(ul_resolve_t *uidx, ma_ug_t *ug, int32_t max_ext, uint32_t max_ext_hifi, 
float len_rat, uint32_t is_trio, uint32_t topo_level, uint32_t *max_drop_len, asg64_v *in, asg64_v *ib)
{
    asg64_v tx = {0,0,0}, tb = {0,0,0}, *b = NULL, *ub = NULL; asg_t *g = ug->g;
    uint32_t v, w, i, k, kv, nv, kw, nw, cnt = 0, n_vtx = g->n_seq<<1, mm_ul, ul_max, to_del;
    asg_arc_t *av, *aw, *ve, *we; uint32_t ul_cnt, pb;

    b = (in?(in):(&tx)); ub = (ib?(ib):(&tb));
    for (v = b->n = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;

        if(g->seq_vis[v] == 0) {
            av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
            if (nv < 2) continue;

            for (i = kv = 0, ul_max = -1; i < nv; ++i) {
                if(av[i].del) continue; kv++;
                get_ul_path_info(uidx, ug, av[i].v, NULL, NULL, &ul_cnt, NULL, NULL, NULL);
                if(ul_cnt > ul_max) ul_max = ul_cnt;
            }
            // if((v>>1) == 5958) {
            //     fprintf(stderr, "[M::%s::] v>>1::%u, v&1::%u, kv::%u, ul_max::%u\n", 
            //     __func__, v>>1, v&1, kv, ul_max);
            // }
            if(kv < 2 || ul_max == 0) continue; //note: ul_cnt/ul_max might be 0

            for (i = 0; i < nv; ++i) {
                if(av[i].del) continue;
                get_ul_path_info(uidx, ug, av[i].v, NULL, NULL, &ul_cnt, NULL, NULL, NULL);
                // if((v>>1) == 5958) {
                //     fprintf(stderr, "+++[M::%s::] v>>1::%u, v&1::%u, av[i].v>>1::%u, av[i].v&1::%u, ul_cnt::%u\n", 
                //     __func__, v>>1, v&1, av[i].v>>1, av[i].v&1, ul_cnt);
                // }
                if(ul_cnt > 0) ul_cnt = ul_max/ul_cnt;
                else ul_cnt = ul_max<<2;
                ul_cnt = ((uint32_t)-1) - ul_cnt;
                kv_push(uint64_t, *b, (((uint64_t)ul_cnt)<<32) | ((uint64_t)(av-g->arc+i)));
            }
        }
    }

    radix_sort_srt64(b->a, b->a + b->n);
    for (k = 0; k < b->n; k++) {
        if(g->arc[(uint32_t)b->a[k]].del) continue;
        v = g->arc[(uint32_t)b->a[k]].ul>>32; w = g->arc[(uint32_t)b->a[k]].v^1;
        if(g->seq[v>>1].del || g->seq[w>>1].del) continue;
        nv = asg_arc_n(g, v); av = asg_arc_a(g, v); 
        if(nv<=1 && asg_arc_n(g, w) <= 1) continue;

        ve = &(g->arc[(uint32_t)b->a[k]]); pb = b->n; ub->n = 0;
        get_ul_path_info(uidx, ug, ve->v, NULL, NULL, &mm_ul, NULL, NULL, NULL);
        // if((v>>1) == 5958 && (w>>1) == 11105) {
        //     fprintf(stderr, ">>>[M::%s::] v>>1::%u, v&1::%u, w>>1::%u, w&1::%u, mm_ul::%u\n", 
        //     __func__, v>>1, v&1, w>>1, w&1, mm_ul);
        // }
        if(mm_ul == 0) continue;///no UL read support this path

        for (i = kv = 0; i < nv; ++i) {
            if(av[i].del) continue; 
            kv++;
            if(av[i].v == ve->v) continue;
            get_ul_path_info(uidx, ug, av[i].v, NULL, NULL, &ul_cnt, NULL, NULL, NULL);
            // if((v>>1) == 5958 && (w>>1) == 11105) {
            //     fprintf(stderr, "*[M::%s::i->%u] v>>1::%u, v&1::%u, w>>1::%u, w&1::%u, av[i].v>>1::%u, av[i].v&1::%u, ul_cnt::%u\n", 
            //     __func__, i, v>>1, v&1, w>>1, w&1, av[i].v>>1, av[i].v&1, ul_cnt);
            // }
            if (mm_ul <= ul_cnt*len_rat) {
                ul_cnt = ((uint32_t)-1) - ul_cnt;
                kv_push(uint64_t, *b, (((uint64_t)ul_cnt)<<32)|((uint64_t)(i)));
            }
        }
        if(b->n == pb) continue;
        to_del = 0; assert(kv >= 2); 

        aw = asg_arc_a(g, w); nw = asg_arc_n(g, w);
        for (i = kw = 0, we = NULL; i < nw; ++i) {
            if (aw[i].del) continue;
            if (aw[i].v == (v^1)) we = &(aw[i]);
            kw++;
        }

        if((kv > 1 && kw > 1) || (usg_bridge_topocut_aux(uidx, ug, ve->v, max_ext, max_ext_hifi, topo_level, 0, b, ub))) {
            radix_sort_srt64(b->a + pb, b->a + b->n);
            for (i = pb; i < b->n; i++) {
                // if((v>>1) == 5958 && (w>>1) == 11105) {
                //     fprintf(stderr, "#[M::%s::srt_i->%u] v>>1::%u, v&1::%u, w>>1::%u, w&1::%u, av[(uint32_t)b->a[i]].v>>1::%u, av[(uint32_t)b->a[i]].v&1::%u\n", 
                //     __func__, (uint32_t)b->a[i], v>>1, v&1, w>>1, w&1, av[(uint32_t)b->a[i]].v>>1, av[(uint32_t)b->a[i]].v&1);
                // }
                if(usg_misjoin_topocut_aux(uidx, ug, ve->v, av[(uint32_t)b->a[i]].v, v, b, ub)) {
                    to_del = 1;
                    break;
                }
            }
        }   

        b->n = pb;
        if(to_del) {
            if(kv > 1 && kw > 1) {
                ve->del = we->del = 1; ++cnt;
            } else {
                pb = b->n; 
                get_ul_path_info(uidx, ug, ve->v, NULL, NULL, NULL, NULL, NULL, b);
                for (i = pb; i < b->n; i++) ulg_seq_del(ug, (b->a[i]>>1));
                cnt += b->n - pb; b->n = pb; 
            }
        }
    }


    if(!in) free(tx.a); if(!ib) free(tb.a);
    if (cnt > 0) asg_cleanup(g);
    return cnt;
}

uint32_t is_het_bridge(ul_resolve_t *uidx, uint64_t *p_a, int64_t p_n, int64_t match_bound)
{
    bubble_type *bub = uidx->bub; int64_t z;
    for (z = match_bound; z >= 0 && IF_HOM((p_a[z]>>1), *bub); z--) {
        // if(is_debug){
        //     fprintf(stderr, "+[M::%s::] match_bound::%ld, p_n::%ld, (p_a[%ld]>>1)::%lu, IF_HOM::%u\n", 
        //     __func__, match_bound, p_n, z, (p_a[z]>>1), IF_HOM((p_a[z]>>1), *bub));
        // }
    } 
    if(z < 0) return 0;
    for (z = match_bound + 1; z < p_n && IF_HOM((p_a[z]>>1), *bub); z++) {
        // if(is_debug){
        //     fprintf(stderr, "-[M::%s::] match_bound::%ld, p_n::%ld, (p_a[%ld]>>1)::%lu, IF_HOM::%u\n", 
        //     __func__, match_bound, p_n, z, (p_a[z]>>1), IF_HOM((p_a[z]>>1), *bub));
        // }
    }
    // if(is_debug) {
    //     fprintf(stderr, "*[M::%s::] z::%ld, p_n::%ld\n", __func__, z, p_n);
    // }
    if(z >= p_n) return 0;
    return 1;
}

uint32_t get_ul_arc_supports(ul_resolve_t *uidx, asg_arc_t *ve, asg64_v *b_int, asg64_v *b_raw, uint64_t skip_hom, uint64_t *retrun_w_v, uint64_t *retrun_w_r)
{
    ma_ug_t *iug = uidx->uovl.i_ug; asg_t *g = iug->g; uint32_t v, k, z, zn, nv, n_pre, l_v, l_r, skip_hom_local; 
    asg_arc_t *av; v = ve->ul>>32; uint32_t b_int_s = b_int->n, b_raw_s = b_raw->n, is_collapse = 0, v_occ, r_occ;
    (*retrun_w_v) = (*retrun_w_r) = (uint64_t)-1; uint64_t *raw_v, *raw_r, w_v, w_r, min_w_v, min_w_r;

    get_ul_path_info(uidx, iug, v^1, NULL, NULL, NULL, NULL, NULL, b_int); nv = b_int->n-b_int_s;
    for (k = 0, n_pre = b_int->n; k < (nv>>1); k++) {
        v = b_int->a[k+b_int_s]; 
        b_int->a[k+b_int_s] = b_int->a[b_int_s+nv-k-1]^1; 
        b_int->a[b_int_s+nv-k-1] = v^1;
    }
    if(nv&1) b_int->a[k+b_int_s] ^= 1;
    v = ve->ul>>32;

    get_ul_path_info(uidx, iug, ve->v, NULL, NULL, &v_occ, NULL, NULL, b_int);
    // if(((ve->ul>>32) == 24093) && (ve->v == 61472)) {
    //     for (z = 0; z < b_int->n-b_int_s; z++) {
    //         fprintf(stderr, "[M::%s::] integ_v[%u]>>1:%lu, integ_v[%u]&1:%lu\n", __func__, 
    //         z, b_int->a[b_int_s+z]>>1, z, b_int->a[b_int_s+z]&1);
    //     }
    // }
    l_v = gen_ug_integer_seq_on_fly(uidx, b_int->a+b_int_s, b_int->n-b_int_s, b_raw);
    
    av = asg_arc_a(g, v); nv = asg_arc_n(g, v); 
    for (k = 0, min_w_v = min_w_r = (uint64_t)-1; k < nv; k++) {
        if(av[k].del || av[k].v == ve->v) continue;///skip ve->v

        b_int->n = n_pre; b_raw->n = b_raw_s + l_v;
        get_ul_path_info(uidx, iug, av[k].v, NULL, NULL, &r_occ, NULL, NULL, b_int);
        // if(((ve->ul>>32) == 24093) && (ve->v == 61472)) {
        //     for (z = 0; z < b_int->n-b_int_s; z++) {
        //         fprintf(stderr, "[M::%s::k->%u] integ_r[%u]>>1:%lu, integ_r[%u]&1:%lu\n", __func__, 
        //         k, z, b_int->a[b_int_s+z]>>1, z, b_int->a[b_int_s+z]&1);
        //     }
        // }

        l_r = gen_ug_integer_seq_on_fly(uidx, b_int->a+b_int_s, b_int->n-b_int_s, b_raw);
        ///raw_v -> ve->v; raw_r -> av[k].v
        raw_v = b_raw->a + b_raw_s; raw_r = b_raw->a + b_raw_s + l_v; zn = MIN(l_v, l_r);
        // if((v>>1) == 77 && (ve->v>>1) == 78) {
        //     fprintf(stderr, "[M::%s::] l_v::%u\n", __func__, l_v);
        //     for (z = 0; z < l_v; z++) {
        //         fprintf(stderr, "raw_v[%u]>>1::utg%.6dl, raw_v[%u]&1::%lu\n", 
        //         z, (int32_t)(raw_v[z]>>1)+1, z , raw_v[z]&1);
        //     }
        //     fprintf(stderr, "[M::%s::] l_r::%u\n", __func__, l_r);
        //     for (z = 0; z < l_r; z++) {
        //         fprintf(stderr, "raw_r[%u]>>1::utg%.6dl, raw_r[%u]&1::%lu\n", 
        //         z, (int32_t)(raw_r[z]>>1)+1, z, raw_r[z]&1);
        //     }
        // }
        for (z = 0; z < zn && raw_v[z] == raw_r[z]; z++); ///z: first raw unitig that is different between two paths
        skip_hom_local = skip_hom;
        if(skip_hom_local && z < l_v) {
            skip_hom_local = is_het_bridge(uidx, raw_v, l_v, z - 1);///, (v>>1) == 191 && (ve->v>>1) == 450);
        }
        if(skip_hom_local && z < l_r) {
            skip_hom_local = is_het_bridge(uidx, raw_r, l_r, z - 1);///, (v>>1) == 191 && (ve->v>>1) == 450);
        }

        assert(z > 0); w_v = w_r = (uint64_t)-1;
        if(z < l_v) get_integer_seq_ovlps(uidx, raw_v, l_v, z - 1, skip_hom_local, NULL, &w_v);
        if(z < l_r) get_integer_seq_ovlps(uidx, raw_r, l_r, z - 1, skip_hom_local, NULL, &w_r);
        if(w_v == (uint64_t)-1) w_v = 0;
        if(w_r == (uint64_t)-1) w_r = 0;
        // if((v>>1) == 409 && (ve->v>>1) == 407) {
        //     fprintf(stderr, "[M::%s::v>>1::%u] l_v::%u, l_r::%u, z::%u, w_v::%lu, w_r::%lu, skip_hom_local::%u\n", 
        //     __func__, av[k].v>>1, l_v, l_r, z, w_v, w_r, skip_hom_local);
        // }
        ///z == zn: -> prefer collapse
        if((min_w_v == (uint64_t)-1) || (z == zn) || (min_w_v > w_v) || (min_w_v == w_v && min_w_r < w_r)) {
            min_w_v = w_v; min_w_r = w_r;
            if(z == zn && v_occ < r_occ) is_collapse = 1;
        }
    }
    if(is_collapse) min_w_v = 0;

    (*retrun_w_v) = min_w_v; (*retrun_w_r) = min_w_r;
    return is_collapse;
}

static void worker_update_ul_arc_supports(void *data, long i, int tid) // callback for kt_for()
{
    ul_resolve_t *uidx = (ul_resolve_t *)data; integer_t *buf = &(uidx->str_b.buf[tid]);
    uint64_t *x = &(uidx->uovl.iug_tra->a[i]); asg_arc_t *ve = &(uidx->uovl.i_ug->g->arc[*x]);
    asg64_v b_v, b_r; uint64_t w_v, w_r;
    b_v.a = buf->u.a; b_v.n = buf->u.n; b_v.m = buf->u.m; 
    b_r.a = buf->o.a; b_r.n = buf->o.n; b_r.m = buf->o.m; 

    b_v.n = b_r.n = 0;
    get_ul_arc_supports(uidx, ve, &b_v, &b_r, 1, &w_v, &w_r);

    buf->u.a = b_v.a; buf->u.n = b_v.n; buf->u.m = b_v.m; 
    buf->o.a = b_r.a; buf->o.n = b_r.n; buf->o.m = b_r.m; 

    (*x) |= (w_v<<32);
}


uint32_t check_ul_contain_arc_supports(ul_resolve_t *uidx, asg_arc_t *ve, asg64_v *b_int, asg64_v *b_raw, uint64_t skip_hom, uint64_t *retrun_w_v, uint64_t *retrun_w_r)
{
    ma_ug_t *iug = uidx->uovl.i_ug; asg_t *g = iug->g; uint32_t v, k, z, zn, nv, n_pre, l_v, l_r, skip_hom_local; 
    asg_arc_t *av; v = ve->ul>>32; uint32_t b_int_s = b_int->n, b_raw_s = b_raw->n, is_collapse = 0, v_occ, r_occ;
    (*retrun_w_v) = (*retrun_w_r) = (uint64_t)-1; uint64_t *raw_v, *raw_r, w_v, w_r, min_w_v, min_w_r;

    get_ul_path_info(uidx, iug, v^1, NULL, NULL, NULL, NULL, NULL, b_int); nv = b_int->n-b_int_s;
    for (k = 0, n_pre = b_int->n; k < (nv>>1); k++) {
        v = b_int->a[k+b_int_s]; 
        b_int->a[k+b_int_s] = b_int->a[b_int_s+nv-k-1]^1; 
        b_int->a[b_int_s+nv-k-1] = v^1;
    }
    if(nv&1) b_int->a[k+b_int_s] ^= 1;
    v = ve->ul>>32;

    get_ul_path_info(uidx, iug, ve->v, NULL, NULL, &v_occ, NULL, NULL, b_int);
    l_v = gen_ug_integer_seq_on_fly(uidx, b_int->a+b_int_s, b_int->n-b_int_s, b_raw);
    
    av = asg_arc_a(g, v); nv = asg_arc_n(g, v); 
    for (k = 0, min_w_v = min_w_r = (uint64_t)-1; k < nv; k++) {
        if(av[k].del || av[k].v == ve->v) continue;

        b_int->n = n_pre; b_raw->n = b_raw_s + l_v;
        get_ul_path_info(uidx, iug, av[k].v, NULL, NULL, &r_occ, NULL, NULL, b_int);
        l_r = gen_ug_integer_seq_on_fly(uidx, b_int->a+b_int_s, b_int->n-b_int_s, b_raw);

        raw_v = b_raw->a + b_raw_s; raw_r = b_raw->a + b_raw_s + l_v; zn = MIN(l_v, l_r);
        for (z = 0; z < zn && raw_v[z] == raw_r[z]; z++); 
        skip_hom_local = skip_hom;
        if(skip_hom_local && z < l_v) {
            skip_hom_local = is_het_bridge(uidx, raw_v, l_v, z - 1);///, (v>>1) == 191 && (ve->v>>1) == 450);
        }
        if(skip_hom_local && z < l_r) {
            skip_hom_local = is_het_bridge(uidx, raw_r, l_r, z - 1);///, (v>>1) == 191 && (ve->v>>1) == 450);
        }

        assert(z > 0); w_v = w_r = (uint64_t)-1;
        if(z < l_v) get_integer_seq_ovlps(uidx, raw_v, l_v, z - 1, skip_hom_local, NULL, &w_v);
        if(z < l_r) get_integer_seq_ovlps(uidx, raw_r, l_r, z - 1, skip_hom_local, NULL, &w_r);
        if(w_v == (uint64_t)-1) w_v = 0;
        if(w_r == (uint64_t)-1) w_r = 0;
        // if((v>>1) == 409 && (ve->v>>1) == 407) {
        //     fprintf(stderr, "[M::%s::v>>1::%u] l_v::%u, l_r::%u, z::%u, w_v::%lu, w_r::%lu, v_occ::%u, r_occ::%u, skip_hom_local::%u\n", 
        //     __func__, av[k].v>>1, l_v, l_r, z, w_v, w_r, v_occ, r_occ, skip_hom_local);
        // }
        ///z == zn: -> prefer collapse
        if((min_w_v == (uint64_t)-1) || (z == zn) || (min_w_v > w_v) || (min_w_v == w_v && min_w_r < w_r)) {
            min_w_v = w_v; min_w_r = w_r;
            if(z == zn && v_occ < r_occ) {
                is_collapse = 1;
            }
        }
    }

    (*retrun_w_v) = min_w_v; (*retrun_w_r) = min_w_r;
    return is_collapse;
}

uint32_t check_ulg_to_del(ul_resolve_t *uidx, ma_ug_t *ug, uint32_t v, uint32_t w, uint32_t kv, uint32_t kw, 
uint32_t max_ext, uint32_t max_ext_hifi, uint32_t topo_level, uint32_t collapse, asg64_v *b, asg64_v *ub)
{
    uint32_t to_del = 0;
    if(collapse) topo_level = 3; 
    if(topo_level == 0) {
        to_del = 1;
    } else if(topo_level == 2) {
        if (kv > 1 && kw > 1) to_del = 1;
    } else {
        if (kv > 1 && kw > 1) {
            to_del = 1;
        } else if (kw == 1) {
            if (usg_topocut_aux(uidx, ug, w^1, max_ext, max_ext_hifi, b, ub)) to_del = 1;                    
        } else if (kv == 1) {
            if (usg_topocut_aux(uidx, ug, v^1, max_ext, max_ext_hifi, b, ub)) to_del = 1;   
        }
    }

    return to_del;
}

uint32_t ulg_arc_cut_supports(ul_resolve_t *uidx, ma_ug_t *ug, int32_t max_ext, uint32_t max_ext_hifi, 
float len_rat, uint32_t is_trio, uint32_t topo_level, uint32_t skip_hom, uint32_t *max_drop_len, uint32_t collapse_check,
asg64_v *in, asg64_v *ib)
{
    // fprintf(stderr, "\n[M::%s::] max_ext::%d, max_ext_hifi::%d, len_rat::%f, is_trio::%u, topo_level::%u, skip_hom::%u, collapse_check::%u\n", 
    //         __func__, max_ext, max_ext_hifi, len_rat, is_trio, topo_level, skip_hom, collapse_check);

    asg64_v tx = {0,0,0}, tb = {0,0,0}, *b = NULL, *ub = NULL; asg_t *g = ug->g;
    uint32_t v, w, i, k, kv, nv, kw, nw, cnt = 0, n_vtx = g->n_seq<<1, to_del, collapse;
    asg_arc_t *av, *aw, *ve, *we; uint64_t w_q, w_t, pb;

    b = (in?(in):(&tx)); ub = (ib?(ib):(&tb));
    for (v = b->n = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;

        if(g->seq_vis[v] == 0) {
            av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
            if (nv < 2) continue;

            for (i = kv = 0; i < nv && kv < 2; ++i) {
                if(av[i].del) continue; kv++;
            }
            if(kv < 2) continue;

            for (i = 0; i < nv; ++i) {
                if(av[i].del) continue;
                kv_push(uint64_t, *b, ((uint64_t)(av-g->arc+i)));
            }
        }
    }

    // fprintf(stderr, "\n#[M::%s::] Starting...\n", __func__);
    uidx->uovl.iug_tra = b;
    kt_for(uidx->str_b.n_thread, worker_update_ul_arc_supports, uidx, b->n);///all ul + ug
    uidx->uovl.iug_tra = NULL;
    // fprintf(stderr, "#[M::%s::] Done\n", __func__);
    // fprintf(stderr, "#[M::%s::] collapse_check::%u, len_rat::%f\n", __func__, collapse_check, len_rat);
    
    radix_sort_srt64(b->a, b->a + b->n);
    for (k = 0; k < b->n; k++) {
        if(g->arc[(uint32_t)b->a[k]].del) continue;
        v = g->arc[(uint32_t)b->a[k]].ul>>32; w = g->arc[(uint32_t)b->a[k]].v^1;
        if(g->seq[v>>1].del || g->seq[w>>1].del) continue;
        nv = asg_arc_n(g, v); av = asg_arc_a(g, v); 
        nw = asg_arc_n(g, w); aw = asg_arc_a(g, w); 
        if(nv <= 1 && nw <= 1) continue;
        ve = &(g->arc[(uint32_t)b->a[k]]);
        for (i = kv = 0; i < nv; ++i) {
            if(av[i].del) continue; 
            kv++;
        }
        for (i = kw = 0, we = NULL; i < nw; ++i) {
            if (aw[i].del) continue;
            if (aw[i].v == (v^1)) we = &(aw[i]);
            kw++;
        }
        if(kv <= 1 && kw <= 1) continue;
        collapse = 0;

        if(collapse_check) {
            if(kv > 1) {
                pb = b->n; ub->n = 0;
                collapse = check_ul_contain_arc_supports(uidx, ve, b, ub, skip_hom, &w_q, &w_t);
                b->n = pb; ub->n = 0;
            }

            if(collapse == 0 && kw > 1) {
                pb = b->n; ub->n = 0;
                collapse = check_ul_contain_arc_supports(uidx, we, b, ub, skip_hom, &w_q, &w_t);
                b->n = pb; ub->n = 0;
            }
        }

        // if(((v>>1) == 6788 && (w>>1) == 17213) || ((w>>1) == 6788 && (v>>1) == 17213)) {
        //     fprintf(stderr, "#[M::%s::] v>>1::%u, v&1::%u, kv::%u, w>>1::%u, w&1::%u, kw::%u, collapse::%u\n", 
        //     __func__, v>>1, v&1, kv, w>>1, w&1, kw, collapse);
        // }

        if(collapse == 0) {
            if(kv > 1) {
                pb = b->n; ub->n = 0;
                get_ul_arc_supports(uidx, ve, b, ub, skip_hom, &w_q, &w_t);
                b->n = pb; ub->n = 0;
                // if((v>>1) == 409 && (w>>1) == 407) {
                //     fprintf(stderr, "+[M::%s::] v>>1::%u, v&1::%u, kv::%u, w>>1::%u, w&1::%u, kw::%u, w_q::%lu, w_t::%lu\n", 
                //     __func__, v>>1, v&1, kv, w>>1, w&1, kw, w_q, w_t);
                // }
                // if(((v>>1) == 6788 && (w>>1) == 17213) || ((w>>1) == 6788 && (v>>1) == 17213)) {
                //     fprintf(stderr, "+[M::%s::] v>>1::%u, v&1::%u, w_q::%lu, w>>1::%u, w&1::%u, w_t::%lu\n", 
                //     __func__, v>>1, v&1, w_q, w>>1, w&1, w_t);
                // }
                if(w_q == (uint64_t)-1) continue;
                if(w_q > w_t*len_rat) continue;
            }

            if(kw > 1) {
                pb = b->n; ub->n = 0;
                get_ul_arc_supports(uidx, we, b, ub, skip_hom, &w_q, &w_t);
                b->n = pb; ub->n = 0;
                // if((v>>1) == 409 && (w>>1) == 407) {
                //     fprintf(stderr, "-[M::%s::] v>>1::%u, v&1::%u, kv::%u, w>>1::%u, w&1::%u, kw::%u, w_q::%lu, w_t::%lu\n", 
                //     __func__, v>>1, v&1, kv, w>>1, w&1, kw, w_q, w_t);
                // }
                // if(((v>>1) == 6788 && (w>>1) == 17213) || ((w>>1) == 6788 && (v>>1) == 17213)) {
                //     fprintf(stderr, "-[M::%s::] v>>1::%u, v&1::%u, w_q::%lu, w>>1::%u, w&1::%u, w_t::%lu\n", 
                //     __func__, v>>1, v&1, w_q, w>>1, w&1, w_t);
                // }
                if(w_q == (uint64_t)-1) continue;
                if(w_q > w_t*len_rat) continue;
            }
        }

        to_del = check_ulg_to_del(uidx, ug, v, w, kv, kw, max_ext, max_ext_hifi, topo_level, collapse, b, ub);
        
        // if(((v>>1) == 6788 && (w>>1) == 17213) || ((w>>1) == 6788 && (v>>1) == 17213)) {
        //     fprintf(stderr, "#[M::%s::] v>>1::%u, v&1::%u, kv::%u, w>>1::%u, w&1::%u, kw::%u, to_del::%u\n", 
        //     __func__, v>>1, v&1, kv, w>>1, w&1, kw, to_del);
        // }
        if (to_del) {
            ve->del = we->del = 1, ++cnt;
        }
    }

    if(!in) free(tx.a); if(!ib) free(tb.a);
    if (cnt > 0) asg_cleanup(g);
    return cnt;
}

uint64_t ulg_bub_pop_cut_aux(ul_resolve_t *uidx, ma_ug_t *ug, uint32_t v0, buf_t *x, uint32_t max_ext, uint32_t max_ext_hifi, asg64_v *b, asg64_v *ub)
{
    uint32_t i, v, u, bn = b->n, n_ext, ul_occ, r = 0; binfo_t *t;
    v = x->S.a[0]; 
    do {
		u = x->a[v].p; // u->v
        x->a[v].d = (uint32_t)-1;
		v = u;
	} while (v != v0);

    for (i = n_ext = 0; i < x->b.n; ++i) { // clear the states of visited vertices
		t = &x->a[x->b.a[i]];
        if(t->d == (uint32_t)-1) continue;
        get_iug_u_raw_occ(uidx, x->b.a[i]>>1, &ul_occ, NULL); 
        n_ext += ul_occ; kv_push(uint64_t, *b, x->b.a[i]);
	}

    if(n_ext <= max_ext) {
        if(get_remove_hifi_occ(uidx, max_ext_hifi, b->a + bn, b->n - bn, ub, NULL)) r = 1;
    }
    b->n = bn;
    return r;
}

uint32_t ulg_bub_pop_backtrack(ma_ug_t *ug, uint32_t v0, buf_t *b)
{
	uint32_t i, v, u, cnt = b->e.n; asg_t *g = ug->g;

    ///b->S.a[0] is the sink of this bubble
	for (i = 0; i < b->b.n; ++i) g->seq[b->b.a[i]>>1].del = 1;

    ///remove all edges (self/reverse for each edge) in this bubble
	for (i = 0; i < b->e.n; ++i) {
        g->arc[b->e.a[i]].del = 1;
        asg_arc_del(g, g->arc[b->e.a[i]].v^1, (g->arc[b->e.a[i]].ul>>32)^1, 1);
	}

    ///v is the sink of this bubble
	v = b->S.a[0];
	do {
		u = b->a[v].p; // u->v
		g->seq[v>>1].del = 0;
		asg_arc_del(g, u, v, 0);
		asg_arc_del(g, v^1, u^1, 0);
        cnt--;
		v = u;
	} while (v != v0);

    for (i = 0; i < b->b.n; ++i) {
        if(!g->seq[b->b.a[i]>>1].del) continue;
        ulg_seq_del(ug, (b->b.a[i]>>1));
    }

    return cnt;
}

uint64_t ulg_bub_pop1(ul_resolve_t *uidx, ma_ug_t *ug, uint32_t v0, uint64_t max_dist, buf_t *x, 
uint32_t max_ext, uint32_t max_ext_hifi, uint32_t check_bubble_only, uint32_t is_pop, uint32_t skip_hom, 
asg64_v *b, asg64_v *ub, uint32_t *r_w_c, uint32_t *r_w_m)
{
    asg_t *g = ug->g; uint32_t pb = b->n; uint64_t w_q, w_t, wc, wm, ww; (*r_w_c) = (*r_w_m) = (uint32_t)-1;
    uint32_t v, w, i, kv, nv, kw, cnt = 0, fail_b = 0, n_tips = 0, tip_end = (uint32_t)-1;
    uint32_t l, d, c, m, n_pending = 0, z, to_replace; asg_arc_t *av, *ve, *we; binfo_t *t;
    if(g->seq[v0>>1].del || get_arcs(g, v0, NULL, 0) < 2) return 0; // already deleted
    // fprintf(stderr, "sbsbsbsbsbsbsb[M::%s::(v0>>1)->%u::(v0&1)->%u] check_bubble_only::%u, is_pop::%u\n", 
    //             __func__, v0>>1, v0&1, check_bubble_only, is_pop);

    x->S.n = x->T.n = x->b.n = x->e.n = 0;
    x->a[v0].c = x->a[v0].d = x->a[v0].m = x->a[v0].nc = x->a[v0].np = 0;
    kv_push(uint32_t, x->S, v0);

    do {
        v = kv_pop(x->S); d = x->a[v].d; c = x->a[v].c; m = x->a[v].m;
        nv = asg_arc_n(g, v); av = asg_arc_a(g, v); kv = get_arcs(g, v, NULL, 0);
        for (i = 0; i < nv; ++i) {
            if (av[i].del) continue;
            w = av[i].v; t = &(x->a[w]); l = ((v == v0)?(0):((uint32_t)av[i].ul));
            if ((w>>1) == (v0>>1)) {
                fail_b = 1;
                break;
            }
            kv_push(uint32_t, x->e, (g->idx[v]>>32) + i); ///for backtracking
            if (d + l > max_dist) {
                fail_b = 1;
                break;
            }
            kw = get_arcs(g, w^1, NULL, 0);
            wc = 0; wm = get_ul_occ(uidx, w>>1);
            if(!check_bubble_only) {
                if(kv > 1) {
                    pb = b->n; ub->n = 0; ve = &(av[i]);
                    // fprintf(stderr, "[M::%s::] ve->v>>1::%lu, ve->v&1::%lu, ve->w>>1::%u, ve->w&1::%u\n", 
                    // __func__, ve->ul>>33, (ve->ul>>32)&1, ve->v>>1, ve->v&1);
                    get_ul_arc_supports(uidx, ve, b, ub, skip_hom, &w_q, &w_t);
                    b->n = pb; ub->n = 0;
                    if(w_q < w_t) {
                        if(w_q == 0) {
                            ww = 10;
                        } else {
                            if(w_t != 0) ww = ((w_t - w_q)*10)/w_t;
                            else ww = 0;
                        } 
                        if(ww > wc) wc = ww;
                    }
                }

                if(kw > 1) {
                    pb = b->n; ub->n = 0; we = get_specfic_edge(g, w^1, v^1); assert(we);
                    // fprintf(stderr, "[M::%s::] we->v>>1::%lu, we->v&1::%lu, we->w>>1::%u, we->w&1::%u\n", 
                    // __func__, we->ul>>33, (we->ul>>32)&1, we->v>>1, we->v&1);
                    get_ul_arc_supports(uidx, we, b, ub, skip_hom, &w_q, &w_t);
                    b->n = pb; ub->n = 0;
                    if(w_q < w_t) {
                        if(w_q == 0) {
                            ww = 10;
                        } else {
                            if(w_t != 0) ww = ((w_t - w_q)*10)/w_t;
                            else ww = 0;
                        } 
                        if(ww > wc) wc = ww;
                    }
                }
            }

            if (t->s == 0) {
                kv_push(uint32_t, x->b, w); 
                t->p = v, t->s = 1, t->d = d + l;
                t->c = c + wc; t->m = m + wm;
                t->r = kw;
				++n_pending;
            } else {
                to_replace = 0;
                if((c + wc) < t->c) {
                    to_replace = 1;
                } else if(((c + wc) == t->c) && (m + wm > t->m)) {
                    to_replace = 1;
                } else if(((c + wc) == t->c) && (m + wm == t->m) && (d + l > t->d)) {
                    to_replace = 1;
                }
                if(to_replace) {
                    t->p = v; t->c = c + wc; t->m = m + wm;
                }
                if (d + l < t->d) t->d = d + l; // update dist
            }

            if (--(t->r) == 0) {
                z = get_arcs(g, w, NULL, 0);
                if(z > 0) {
                    kv_push(uint32_t, x->S, w);
                }
                else {
                    ///at most one tip
                    if(n_tips != 0) {
                        fail_b = 1;
                        break;
                    }
                    n_tips++; tip_end = w;
                }
				--n_pending;
            }
        }
        if(fail_b) break;
        if(n_tips == 1) {
            if(tip_end != (uint32_t)-1 && n_pending == 0 && x->S.n == 0) {
                kv_push(uint32_t, x->S, tip_end);
                break;
            }
            fail_b = 1;
            break;
        }

        if (i < nv || x->S.n == 0) {
            fail_b = 1;
            break;
        }
    } while (x->S.n > 1 || n_pending);

    if(!fail_b) {//there is a bubble
        cnt = 1; (*r_w_c) = x->a[x->S.a[0]].c; (*r_w_m) = x->a[x->S.a[0]].m; 
        if(!check_bubble_only) {
            if(is_pop) {
                cnt = ulg_bub_pop_backtrack(ug, v0, x);
            } else {
                cnt = ulg_bub_pop_cut_aux(uidx, ug, v0, x, max_ext, max_ext_hifi, b, ub);
            }
        }
    }
    for (i = 0; i < x->b.n; ++i) { // clear the states of visited vertices
        // if(v0 == 119) {
        //     fprintf(stderr, "-[M::%s::(v0>>1)->%u::(v0&1)->%u] v[%u]>>1::%u, v[%u]&1::%u, cnt::%u, is_pop::%u\n", 
        //         __func__, v0>>1, v0&1, i, x->b.a[i]>>1, i, x->b.a[i]&1, cnt, is_pop);
        // }
		t = &x->a[x->b.a[i]];
		t->s = t->c = t->d = t->m = t->nc = t->np = 0;
	}
    if(!cnt) (*r_w_c) = (*r_w_m) = (uint32_t)-1;
    return cnt;
}

uint64_t ulg_pop_bubble(ul_resolve_t *uidx, ma_ug_t *ug, uint64_t* i_max_dist, uint32_t max_ext, uint32_t max_ext_hifi, uint32_t skip_hom, asg64_v *in, asg64_v *ib)
{
    // fprintf(stderr, "[M::%s::] Starting...\n", __func__);
    asg_t *g = ug->g; 
    asg64_v tx = {0,0,0}, tb = {0,0,0}, *ob = NULL, *ub = NULL; 
    uint32_t v, w, n_vtx = g->n_seq<<1, n_arc, nv, i, wc[2], wm[2], mm_c, mm_m, mm_v;
    uint64_t n_pop = 0, max_dist;
    asg_arc_t *av = NULL; if (!g->is_symm) asg_symm(g);
    buf_t b; memset(&b, 0, sizeof(buf_t));
    b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    if(i_max_dist) max_dist = (*i_max_dist);
    else max_dist = get_bub_pop_max_dist_advance(g, &b);
    ob = (in?(in):(&tx)); ub = (ib?(ib):(&tb)); ob->n = ub->n = 0;

    if(max_dist > 0) {
        for (v = 0; v < n_vtx; ++v) {
            nv = asg_arc_n(g, v); av = asg_arc_a(g, v);
            if (nv < 2 || g->seq[v>>1].del) continue;
            for (i = n_arc = 0; i < nv; ++i) {
                if (!av[i].del) ++n_arc;
            }
            if (n_arc < 2) continue;
            ///find a bubble
            ob->n = ub->n = 0;
            if(ulg_bub_pop1(uidx, ug, v, max_dist, &b, max_ext, max_ext_hifi, 1, 0, skip_hom, ob, ub, &(wc[0]), &(wm[0]))) {
                w = b.S.a[0]^1; mm_c = mm_m = mm_v = (uint32_t)-1;
                
                ob->n = ub->n = 0;
                if(ulg_bub_pop1(uidx, ug, v, max_dist, &b, max_ext, max_ext_hifi, 0, 0, skip_hom, 
                                                                                ob, ub, &(wc[0]), &(wm[0]))) {
                    mm_c = wc[0]; mm_m = wm[0]; mm_v = v;
                }

                ob->n = ub->n = 0;
                if(ulg_bub_pop1(uidx, ug, w, max_dist, &b, max_ext, max_ext_hifi, 0, 0, skip_hom, 
                                                                                ob, ub, &(wc[1]), &(wm[1]))) {
                    if((wc[1] < mm_c) && (wc[1] == mm_c && wm[1] > mm_m)) {
                        mm_c = wc[1]; mm_m = wm[1]; mm_v = w;
                    }                                                               
                }

                if(mm_v != (uint32_t)-1) {
                    ob->n = ub->n = 0;
                    w = ulg_bub_pop1(uidx, ug, mm_v, max_dist, &b, max_ext, max_ext_hifi, 0, 1, skip_hom, ob, ub, &(wc[0]), &(wm[0]));
                    assert(w);
                    n_pop += w;
                }
            }
        }
    }

    free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
    if(n_pop) asg_cleanup(g);
    if(!in) free(tx.a); if(!ib) free(tb.a);
    // fprintf(stderr, "[M::%s::] Done...\n", __func__);
    return n_pop;
}

/**
void infer_reliable_regions(ul_resolve_t *uidx, asg64_v *b)
{
    ul2ul_idx_t *idx = &(uidx->uovl); ma_ug_t *iug = idx->i_ug; ma_ug_t *raw = uidx->l1_ug;
    uint64_t k, z, v, l, nv, uv, uw, bv, bw, is_bv_ul, is_bw_ul; ma_utg_t *iu; uinfo_srt_warp_t *seq; 
    uint8_t *raw_idx; CALLOC(raw_idx, raw->u.n<<1); 
    uint64_t *iu_idx; CALLOC(iu_idx, iug->u.n); uint32_t *iu_a;
    asg_t *g = asg_init(); asg_arc_t *av;
    

    for (k = l = 0; k < iug->u.n; k++) {
        seq = &(idx->cc.iug_a[k]);
        iu_idx[k] = l; iu_idx[k] <<= 32; iu_idx[k] += (l + seq->n);
        l += seq->n;
    }

    MALLOC(iu_a, l);
    for (k = b->n = 0; k < iug->u.n; k++) {
        seq = &(idx->cc.iug_a[k]);
        v = k<<1; uv = v;
        bv = iug->u.a[uv>>1].a[((uv&1)?(0):(iug->u.a[uv>>1].n-1))]>>32; if(uv&1) bv ^= 1; 
        is_bv_ul = ulg_type((*idx), (bv>>1)); av = asg_arc_a(iug->g, v); nv = asg_arc_n(iug->g, v); 
        for (z = 0; z < nv; z++) {
            if(av[z].del || (av[z].v>>1) >= k) continue;
            uw = av[z].v; bw = iug->u.a[uw>>1].a[((uw&1)?(iug->u.a[uw>>1].n-1):(0))]>>32; 
            if(uw&1) bw ^= 1; is_bw_ul = ulg_type((*idx), (bw>>1));
        }
        



        v = (k<<1)+1; uv = v;
        bv = iug->u.a[uv>>1].a[((uv&1)?(0):(iug->u.a[uv>>1].n-1))]>>32; if(uv&1) bv ^= 1;
        is_bv_ul = ulg_type((*idx), (bv>>1)); av = asg_arc_a(iug->g, v); nv = asg_arc_n(iug->g, v); 
        for (z = 0; z < nv; z++) {
            if(av[z].del || (av[z].v>>1) >= k) continue;
            uw = av[z].v; bw = iug->u.a[uw>>1].a[((uw&1)?(iug->u.a[uw>>1].n-1):(0))]>>32; 
            if(uw&1) bw ^= 1; is_bw_ul = ulg_type((*idx), (bw>>1));
        }
    }
    
}

void fill_u2g(ul_resolve_t *uidx, asg64_v *b, asg64_v *ub)
{
    renew_ul2_utg(uidx);
    infer_reliable_regions(uidx);
}
**/

uint64_t get_ug_integer_seq_occ(ul_resolve_t *uidx, uint64_t *u_a, uint64_t u_n, asg64_v *b)
{
    uint32_t bn = b->n, k, occ; ma_ug_t *raw = uidx->l1_ug;
    gen_ug_integer_seq_on_fly(uidx, u_a, u_n, b);
    for (k = bn, occ = 0; k < b->n; k++) occ += raw->u.a[b->a[k]>>1].n;
    b->n = bn;
    return occ;
}

uint32_t ul_occ_check(ul_resolve_t *uidx, ma_ug_t *ug, uint32_t qocc_ul, uint32_t qocc_hifi, uint32_t tv, float occ_rate, asg64_v *b, asg64_v *ub)
{
    uint32_t bn = b->n, z, tocc_ul, tocc_hifi, ul;
    get_ul_path_info(uidx, ug, tv, NULL, NULL, NULL, NULL, NULL, b);
    for (z = bn, tocc_ul = 0; z < b->n; z++) {
        get_iug_u_raw_occ(uidx, b->a[z]>>1, &ul, NULL); tocc_ul += ul;
    }
    tocc_hifi = get_ug_integer_seq_occ(uidx, b->a + bn, b->n - bn, ub);
    b->n = bn;
    if((qocc_ul <= (tocc_ul*occ_rate)) && (qocc_hifi <= (tocc_hifi*occ_rate))) return 1;
    return 0;
}

uint32_t idx_check_rate(uint64_t *a, uint64_t a_n, uint64_t tot, uint64_t winlen, float match_rate)
{
    if(a_n == 0) return 0;
    if(tot < winlen) winlen = tot;
    uint64_t k, mm = a[0]>>32, ks, s, e; int64_t p;
    s = 0; e = winlen;
    for (k = ks = 0; k < a_n; k++) {
        assert(a[k] != (uint64_t)-1);
        // if(mm != (a[k]>>32)) {
        //     fprintf(stderr, "***[M::%s::] mm::%lu, (a[%lu]>>32)::%lu\n", 
        //     __func__, mm, k, (a[k]>>32));
        // }
        assert(mm == (a[k]>>32));
        if(((uint32_t)a[k]) >= s && ((uint32_t)a[k]) < e) continue;
        break;
    }
    p = k;
    if(p > 0 && p >= (int64_t)(winlen*match_rate)) return 1;

    for (; k < a_n; k++) {
        assert(a[k] != (uint64_t)-1);
        assert(mm == (a[k]>>32));
        e = ((uint32_t)a[k]) + 1; p++;
        for (; ks < k; ks++) {
            if(e - ((uint32_t)a[ks]) <= winlen) break;
            p--;
        }
        assert(p >= 0);
        if(p > 0 && p >= (int64_t)(winlen*match_rate)) return 1;
    }
    return 0;
}

uint32_t ul_homo_path_check(ul_resolve_t *uidx, ma_ug_t *ug, uint32_t v, uint32_t w, uint32_t raw_ug_occ, float match_rate, asg64_v *b, asg64_v *rb)
{
    // if(((v>>1) == 8158 && (w>>1) == 22318) || ((v>>1) == 16783 && (w>>1) == 28464)) {
    //     fprintf(stderr, "\n++++++[M::%s::] v>>1::%u, v&1::%u, w>>1::%u, w&1::%u\n", 
    //     __func__, v>>1, v&1, w>>1, w&1);
    // }
    
    bubble_type *bub = uidx->bub; 
    uint64_t k, l, z, bn = b->n, vn = 0, wn = 0, *va, *wa, rbn = rb->n, *rva, *rwa, rvn, rwn, rid, x;
    get_ul_path_info(uidx, ug, v, NULL, NULL, NULL, NULL, NULL, b); vn = b->n - bn;
    get_ul_path_info(uidx, ug, w, NULL, NULL, NULL, NULL, NULL, b); wn = b->n - bn - vn;
    va = b->a + bn; wa = b->a + bn + vn;
    // if(((v>>1) == 8158 && (w>>1) == 22318) || ((v>>1) == 16783 && (w>>1) == 28464)) {
    //     fprintf(stderr, "[M::%s::] vn::%lu, wn::%lu\n", __func__, vn, wn);
    // }

    gen_ug_integer_seq_on_fly(uidx, va, vn, rb); rvn = rb->n - rbn;
    gen_ug_integer_seq_on_fly(uidx, wa, wn, rb); rwn = rb->n - rbn - rvn;
    rva = rb->a + rbn; rwa = rb->a + rbn + rvn;
    // if(((v>>1) == 8158 && (w>>1) == 22318) || ((v>>1) == 16783 && (w>>1) == 28464)) {
    //     fprintf(stderr, "[M::%s::] rvn::%lu, rwn::%lu, raw_ug_occ::%u\n", __func__, rvn, rwn, raw_ug_occ);
    // }

    b->n = bn;
    for (k = 0; k < rvn; k++) {
        rid = rva[k]>>1;
        if(IF_BUB(rid, *bub)) {
            x = bub->index[rid]; 
            x |= ((uint64_t)(0x80000000)); x <<= 32; ///bubble id
        } else {
            x = rid; ///node id
            x <<= 32;
        }
        x |= (k<<1); kv_push(uint64_t, *b, x);
    }

    for (k = 0; k < rwn; k++) {
        rid = rwa[k]>>1;
        if(IF_BUB(rid, *bub)) {
            x = bub->index[rid]; 
            x |= ((uint64_t)(0x80000000)); x <<= 32; ///bubble id
        } else {
            x = rid; ///node id
            x <<= 32;
        }
        x |= 1; x |= (k<<1); kv_push(uint64_t, *b, x);
    }
    radix_sort_srt64(b->a + bn, b->a + b->n);

    uint64_t o[2], c[2];
    for (l = bn, k = bn + 1, o[0] = o[1] = 0; k <= b->n; k++) {
        if (k == b->n || (b->a[k]>>32) != (b->a[l]>>32)) {
            if((k - l > 1)) {
                for (z = l, c[0] = c[1] = 0; z < k; z++) {
                    c[b->a[z]&1]++;
                    if(c[0] > 0 && c[1] > 0) break;
                }
                if(c[0] > 0 && c[1] > 0) {
                    for (z = l; z < k; z++) {
                        o[b->a[z]&1]++;

                        x = b->a[z]; 
                        x <<= 32; x >>= 32; x >>= 1; 
                        if(b->a[z]&1) x|= ((uint64_t)(0x100000000));
                        b->a[z] = x;
                    }
                } else {
                    for (z = l; z < k; z++) b->a[z] = (uint64_t)-1;
                }
            } else {
                for (z = l; z < k; z++) b->a[z] = (uint64_t)-1;
            }
            l = k;
        }
    }
    // fprintf(stderr, "\n[M::%s::] rvn::%lu, rwn::%lu, bn::%lu, b->n::%lu, o[0]::%lu, o[1]::%lu\n", __func__, 
    // rvn, rwn, bn, (uint64_t)b->n, o[0], o[1]);
    radix_sort_srt64(b->a + bn, b->a + b->n); b->n = bn + o[0] + o[1]; 
    b->n = bn; rb->n = rbn;

    if(idx_check_rate(b->a + bn, o[0], rvn, raw_ug_occ, match_rate)) return 1;
    if(idx_check_rate(b->a + bn + o[0], o[1], rwn, raw_ug_occ, match_rate)) return 1;

    // if(o[0] >= (rvn*match_rate)) return 1;
    // if(o[1] >= (rwn*match_rate)) return 1;

    return 0;
}

///small_occ_rate = 0.15; len_rat = 1.5
uint32_t ulg_arc_cut_z(ul_resolve_t *uidx, ma_ug_t *ug, uint32_t max_ext, uint32_t max_ext_hifi, 
float len_rat, float small_occ_rate, uint32_t raw_ug_occ, float raw_match_rate, uint32_t is_trio, 
uint32_t skip_hom, uint32_t *max_drop_len, asg64_v *in, asg64_v *ib)
{
    asg64_v tx = {0,0,0}, tb = {0,0,0}, *b = NULL, *ub = NULL; asg_t *g = ug->g; 
    uint32_t v, w, wt, z, i, k, kv, nv, kw, kwt, nw, cnt = 0, n_vtx = g->n_seq<<1, ul, vp[2], wp[2];
    asg_arc_t *av, *aw, *ve, *we, *wte; uint64_t w_q, w_t, pb, raw_ul, raw_hifi;

    b = (in?(in):(&tx)); ub = (ib?(ib):(&tb));
    for (v = b->n = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;
        av = asg_arc_a(g, v); nv = asg_arc_n(g, v);
        if (nv < 2) continue;
        for (i = kv = 0; i < nv && kv <= 2; ++i) {
            if(av[i].del) continue; kv++;
        }
        if(kv != 2) continue;
        for (i = 0; i < nv; ++i) {
            if(av[i].del) continue; 
            kw = get_arcs(ug->g, av[i].v^1, NULL, 0);
            // if((v>>1) == 16783) {
            //     fprintf(stderr, "sss[M::%s::] v>>1::%u, v&1::%u, av[i].v>>1::%u, av[i].v&1::%u, kw::%u\n", 
            //     __func__, v>>1, v&1, av[i].v>>1, av[i].v&1, kw);
            // }
            if(kw == 2) {
                kv_push(uint64_t, *b, ((uint64_t)(av-g->arc+i)));
            } else if(kw == 1) {
                ub->n = 0;
                if(get_ul_path_info(uidx, ug, av[i].v, NULL, NULL, NULL, NULL, NULL, ub)==TWO_INPUT) {
                    for (z = raw_ul = 0; z < ub->n; z++) {
                        get_iug_u_raw_occ(uidx, ub->a[z]>>1, &ul, NULL); raw_ul += ul;
                    }
                    kv_push(uint64_t, *b, ((raw_ul<<32)|((uint64_t)(av-g->arc+i))));
                }
            }
        }
    }

    radix_sort_srt64(b->a, b->a + b->n);
    for (k = 0; k < b->n; k++) {
        if(g->arc[(uint32_t)b->a[k]].del) continue;
        v = g->arc[(uint32_t)b->a[k]].ul>>32; w = g->arc[(uint32_t)b->a[k]].v^1;
        if(g->seq[v>>1].del || g->seq[w>>1].del) continue;
        nv = asg_arc_n(g, v); av = asg_arc_a(g, v); 
        nw = asg_arc_n(g, w); aw = asg_arc_a(g, w); 
        if(nv <= 1 && nw <= 1) continue;

        vp[0] = v^1; vp[1] = (uint32_t)-1;
        ve = &(g->arc[(uint32_t)b->a[k]]);
        for (i = kv = 0; i < nv; ++i) {
            if(av[i].del) continue; 
            if(av[i].v != (w^1)) vp[1] = av[i].v;
            kv++;
        }

        wp[0] = w^1; wp[1] = (uint32_t)-1;
        for (i = kw = 0, we = NULL; i < nw; ++i) {
            if (aw[i].del) continue;
            if (aw[i].v == (v^1)) we = &(aw[i]);
            else wp[1] = aw[i].v;
            kw++;
        }
        if(kv <= 1 && kw <= 1) continue;
        if(kv != 2 || kw > 2) continue;
        // if((v>>1) == 16783) {
        //     fprintf(stderr, "bbb[M::%s::] v>>1::%u, v&1::%u, w>>1::%u, w&1::%u, kv::%u, kw::%u\n", 
        //     __func__, v>>1, v&1, w>>1, w&1, kv, kw);
        // }
        raw_ul = 0; wt = w; wte = we; kwt = kw; pb = b->n; ub->n = 0;
        if(kw == 1) {
            if(get_ul_path_info(uidx, ug, w^1, &wt, NULL, NULL, NULL, NULL, b)==TWO_INPUT) {
                for (z = pb, raw_ul = 0; z < b->n; z++) {
                    get_iug_u_raw_occ(uidx, b->a[z]>>1, &ul, NULL); raw_ul += ul;
                }
                ul = wt; get_arcs(ug->g, wt, &wt, 1); wt = ug->g->arc[wt].v^1;
                wp[0] = wt^1; wp[1] = (uint32_t)-1;
                nw = asg_arc_n(g, wt); aw = asg_arc_a(g, wt); 
                for (i = kwt = 0; i < nw; ++i) {
                    if (aw[i].del) continue;
                    if (aw[i].v == (ul^1)) wte = &(aw[i]);
                    else wp[1] = aw[i].v;
                    kwt++;
                }
                // if((v>>1) == 16783) {
                //     fprintf(stderr, "bbb[M::%s::] v>>1::%u, v&1::%u, w>>1::%u, w&1::%u, wt>>1::%u, wt&1::%u, kwt::%u\n", 
                //     __func__, v>>1, v&1, w>>1, w&1, wt>>1, wt&1, kwt);
                // }
                nw = asg_arc_n(g, w); aw = asg_arc_a(g, w); 
            } else {
                b->n = pb;
                continue;
            }
        }
        // if((v>>1) == 16783) {
        //     fprintf(stderr, "bbb[M::%s::] v>>1::%u, v&1::%u, raw_ul::%lu, max_ext::%u, ***1***\n", 
        //     __func__, v>>1, v&1, raw_ul, max_ext);
        // }
        assert(kv == 2 && kwt == 2);
        if(raw_ul <= max_ext && get_remove_hifi_occ(uidx, max_ext_hifi, b->a + pb, b->n - pb, ub, NULL)) {
            raw_hifi = get_ug_integer_seq_occ(uidx, b->a + pb, b->n - pb, ub);
            // if((v>>1) == 16783) {
            //     fprintf(stderr, "bbb[M::%s::] v>>1::%u, v&1::%u, ***2***\n", 
            //     __func__, v>>1, v&1);
            // }

            b->n = pb;
            if(raw_ul > 0) {
                if(!ul_occ_check(uidx, ug, raw_ul, raw_hifi, vp[0], small_occ_rate, b, ub)) continue;
                if(!ul_occ_check(uidx, ug, raw_ul, raw_hifi, vp[1], small_occ_rate, b, ub)) continue;
                if(!ul_occ_check(uidx, ug, raw_ul, raw_hifi, wp[0], small_occ_rate, b, ub)) continue;
                if(!ul_occ_check(uidx, ug, raw_ul, raw_hifi, wp[1], small_occ_rate, b, ub)) continue;
            }
            
            // if((v>>1) == 16783) {
            //     fprintf(stderr, "bbb[M::%s::] v>>1::%u, v&1::%u, ***3***\n", 
            //     __func__, v>>1, v&1);
            // }
            pb = b->n; ub->n = 0;
            get_ul_arc_supports(uidx, ve, b, ub, skip_hom, &w_q, &w_t);
            b->n = pb; ub->n = 0;
            if((w_q == (uint64_t)-1) || (w_q > w_t*len_rat)) continue;

            // if((v>>1) == 16783) {
            //     fprintf(stderr, "bbb[M::%s::] v>>1::%u, v&1::%u, ***4***\n", 
            //     __func__, v>>1, v&1);
            // }

            pb = b->n; ub->n = 0;
            get_ul_arc_supports(uidx, wte, b, ub, skip_hom, &w_q, &w_t);
            b->n = pb; ub->n = 0;
            if((w_q == (uint64_t)-1) || (w_q > w_t*len_rat)) continue;

            // if((v>>1) == 16783) {
            //     fprintf(stderr, "bbb[M::%s::] v>>1::%u, v&1::%u, ***5***\n", 
            //     __func__, v>>1, v&1);
            // }

            if(!ul_homo_path_check(uidx, ug, vp[0], wp[1], raw_ug_occ, raw_match_rate, b, ub)) continue;


            // if((v>>1) == 16783) {
            //     fprintf(stderr, "bbb[M::%s::] v>>1::%u, v&1::%u, ***6***\n", 
            //     __func__, v>>1, v&1);
            // }
            if(!ul_homo_path_check(uidx, ug, wp[0], vp[1], raw_ug_occ, raw_match_rate, b, ub)) continue;

            // if((v>>1) == 16783) {
            //     fprintf(stderr, "bbb[M::%s::] v>>1::%u, v&1::%u, ***7***\n", 
            //     __func__, v>>1, v&1);
            // }

            if(kw == 1) {
                get_ul_path_info(uidx, ug, w^1, &wt, NULL, NULL, NULL, NULL, b);
                for (z = pb; z < b->n; z++) ulg_seq_del(ug, (b->a[z]>>1));
            }
            ve->del = we->del = 1; cnt++;
        }

        b->n = pb;
    }



    if(!in) free(tx.a); if(!ib) free(tb.a);
    if (cnt > 0) asg_cleanup(g);
    // fprintf(stderr, "[M::%s::] cnt::%u\n", __func__, cnt);
    return cnt;
}


usg_seq_t *push_usg_t_node(usg_t *ng, uint64_t id)
{
    if(id >= ng->m) {
        uint64_t m = ng->m;
        kv_resize(usg_seq_t, *ng, id + 1); 
        memset(ng->a + m, 0, (ng->m - m)*(sizeof((*ng->a))));
    }
    if(id >= ng->n) ng->n = id + 1;

    return ng->a + id;
}

static void worker_update_ul_arc_drop(void *data, long i, int tid) // callback for kt_for()
{
    ul_resolve_t *uidx = (ul_resolve_t *)data; integer_t *buf = &(uidx->str_b.buf[tid]);
    ul_str_idx_t *str_idx = &(uidx->pstr); uinfo_srt_warp_t *seq = uidx->uovl.iug_seq;
    uint32_t v = seq->a[i].v, z, vz; uint64_t *hid_a, hid_n; ul_str_t *str;
    int64_t s_n, s, p, p_n = seq->n; uint64_t cutoff = uidx->uovl.iug_cov_thre, occ;

    hid_a = str_idx->occ.a + str_idx->idx.a[v>>1];
    hid_n = str_idx->idx.a[(v>>1)+1] - str_idx->idx.a[v>>1];
    for (z = occ = 0; z < hid_n; z++) {
        str = &(str_idx->str.a[hid_a[z]>>32]); s_n = str->cn;
        if(s_n < 2) continue;
        vz = (uint32_t)(str->a[(uint32_t)hid_a[z]]);
        assert((v>>1) == (vz>>1)); 

        if(v == vz) {
            s = ((uint32_t)hid_a[z]) + 1; p = i + 1; 
            if((s < s_n) && (p < p_n) && ((uint32_t)(str->a[s]) == seq->a[p].v)) {
                occ++;
            }
        } else {
            s = ((int32_t)((uint32_t)hid_a[z]))-1; p = i + 1; 
            if((s >= 0) && (p < p_n) && ((uint32_t)(str->a[s]) == (seq->a[p].v^1))) {
                occ++;
            }
        }
        if(occ >= cutoff) break;
    }

    if(occ < cutoff) kv_push(uint64_t, buf->res_dump, i);
}

inline usg_arc_t* get_usg_arc(usg_t *g, uint32_t v, uint32_t w)
{
    usg_arc_t *av = usg_arc_a(g, v); uint32_t nv = usg_arc_n(g, v), k;
    for (k = 0; k < nv; k++) {
        if(av[k].v == w) break;
    }

    if(k < nv) return (&av[k]);
    return NULL;    
}

void pushp_usg_arc_mm(usg_t *g, uint32_t v, uint32_t w, uint32_t uid, uint32_t off)
{
    usg_arc_mm_t *pm;
    kv_pushp(usg_arc_mm_t, g->a[v>>1].arc_mm[v&1], &pm);
    pm->v = w; pm->uid = uid; pm->off = off;
}

static inline void usg_seq_del(usg_t *g, uint32_t s)
{
	uint32_t i, nv, v; usg_arc_t *av, *p;
	g->a[s].del = 1;
    // fprintf(stderr, "\n#[M::%s::] s::%u\n",  __func__, s);
    v = s<<1; av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
    for (i = 0; i < nv; ++i) {
        av[i].del = 1;
        p = get_usg_arc(g, av[i].v^1, v^1); 
        // if(!p) {
        //     fprintf(stderr, "**+**[M::%s::] v>>1::%u, v&1::%u, av[i].v>>1::%u, av[i].v&1::%u\n", 
        //                                                 __func__, v>>1, v&1, av[i].v>>1, av[i].v&1);
        // }
        p->del = 1;
    }

	v = (s<<1)+1; av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
    for (i = 0; i < nv; ++i) {
        av[i].del = 1;
        p = get_usg_arc(g, av[i].v^1, v^1); 
        // if(!p) {
        //     fprintf(stderr, "**-**[M::%s::] v>>1::%u, v&1::%u, av[i].v>>1::%u, av[i].v&1::%u\n", 
        //                                                 __func__, v>>1, v&1, av[i].v>>1, av[i].v&1);
        // }
        p->del = 1;
    }

    // g->a[s].arc[0].n = g->a[s].arc[1].n = 0;
}

static void worker_clean_usg(void *data, long i, int tid) // callback for kt_for()
{
    usg_t *g = (usg_t *)data; usg_seq_t *z = g->a + i; uint32_t k, m, l, srt; 
    usg_arc_warp *x; usg_arc_mm_warp *y;
    if(z->del) z->arc[0].n = z->arc[1].n = 0;

    x = &(z->arc[0]); y = &(z->arc_mm[0]);
    for (k = m = srt = 0; k < x->n; k++) {
        if(x->a[k].del) continue;
        x->a[m] = x->a[k]; x->a[m].idx = 0;
        if(m > 0 && x->a[m].v < x->a[m-1].v) srt = 1;
        m++;
    }
    x->n = m; if(srt) radix_sort_usg_arc_srt(x->a, x->a + x->n);

    radix_sort_usg_arc_mm_srt(y->a, y->a + y->n); k = l = m = 0;
    while (k < y->n && l < x->n) {
        if(y->a[k].v < x->a[l].v) {
            k++;
        } else if(y->a[k].v > x->a[l].v) {
            l++;
        } else {
            y->a[m].v = y->a[k].v;
            if(m > 0 && y->a[m].v == y->a[m-1].v) {
                x->a[l].idx++;
            } else {
                x->a[l].idx = m; x->a[l].idx <<= 32; x->a[l].idx++;
            }
            m++; k++;
        }
    }
    y->n = m;    



    x = &(z->arc[1]); y = &(z->arc_mm[1]);
    for (k = m = srt = 0; k < x->n; k++) {
        if(x->a[k].del) continue;
        x->a[m] = x->a[k]; x->a[m].idx = 0;
        if(m > 0 && x->a[m].v < x->a[m-1].v) srt = 1;
        m++;
    }
    x->n = m; if(srt) radix_sort_usg_arc_srt(x->a, x->a + x->n);
    
    radix_sort_usg_arc_mm_srt(y->a, y->a + y->n); k = l = m = 0;
    while (k < y->n && l < x->n) {
        if(y->a[k].v < x->a[l].v) {
            k++;
        } else if(y->a[k].v > x->a[l].v) {
            l++;
        } else {
            y->a[m].v = y->a[k].v;
            if(m > 0 && y->a[m].v == y->a[m-1].v) {
                x->a[l].idx++;
            } else {
                x->a[l].idx = m; x->a[l].idx <<= 32; x->a[l].idx++;
            }
            m++; k++;
        }
    }
    y->n = m;    

}

void usg_cleanup(usg_t *g)
{
    kt_for(asm_opt.thread_num, worker_clean_usg, g, g->n);
}

static inline int usg_end(const usg_t *g, uint32_t v, uint64_t *lw)
{
    ///v^1 is the another direction of v
    uint32_t w, nv, nw, nw0, nv0 = usg_arc_n(g, v^1);
    int i, i0 = -1;
    usg_arc_t *aw, *av = usg_arc_a(g, v^1);

    ///if this arc has not been deleted
    for (i = nv = 0; i < (int)nv0; ++i)
        if (!av[i].del) i0 = i, ++nv;

    ///end without any out-degree
    if (nv == 0) return ASG_ET_TIP; // tip
    if (nv > 1) return ASG_ET_MULTI_OUT; // multiple outgoing arcs
    ///until here, nv == 1
    if (lw) *lw = ((uint64_t)(v^1))<<32 | av[i0].v;
    w = av[i0].v^1;
    nw0 = usg_arc_n(g, w); aw = usg_arc_a(g, w);
    for (i = nw = 0; i < (int)nw0; ++i)
        if (!aw[i].del) ++nw;

    if (nw != 1) return ASG_ET_MULTI_NEI;
    return ASG_ET_MERGEABLE;
}

uint32_t usg_real_tip(usg_t *g, uint32_t v0, double rate)
{
    usg_arc_t *av, *aw; uint32_t nv, nw, i, k, ov, ow;
    av = usg_arc_a(g, v0); 
    nv = usg_arc_n(g, v0);
    for (i = 0; i < nv; i++) {
        if(av[i].del) continue;
        aw = usg_arc_a(g, (av[i].v^1)); 
        nw = usg_arc_n(g, (av[i].v^1)); 
        for (k = ov = ow = 0; k < nw; k++) {
            if(aw[k].del) continue;
            if((aw[k].v>>1) == (v0>>1)) {
                if(aw[k].ol > ov) ov = aw[k].ol;
            } else {
                if(aw[k].ol > ow) ow = aw[k].ol;
            }
        }
        if(ow == 0 || ov >= ow) return 0;
        if(ov > (ow*rate)) return 0;
    }
    return 1;
}


uint32_t usg_arc_cut_tips(usg_t *g, uint32_t max_ext, uint32_t ignore_ul, asg64_v *in)
{
    asg64_v tx = {0,0,0}, *b = NULL;
    uint32_t n_vtx = g->n<<1, v, w, i, k, cnt = 0, nv, kv, pb, ff, is_telo;
    usg_arc_t *av = NULL, *p; uint64_t lw;
    if(in) b = in;
    else b = &tx;

    for (v = b->n = 0; v < n_vtx; ++v) {
        if (g->a[v>>1].del) continue;

        av = usg_arc_a(g, v^1); nv = usg_arc_n(g, v^1);
        for (i = kv = 0; i < nv; i++) {
            if (av[i].del) continue;
            kv++; break;
        }

        if(kv) continue;
        
        is_telo = 0; if(g->a[v>>1].telo) is_telo = 1;
        for (i = 0, w = v, kv = g->a[v>>1].occ; i < max_ext; i++) {
            if(usg_end(g, w^1, &lw)!=0) break;
            w = (uint32_t)lw; kv += g->a[w>>1].occ;
            if(g->a[w>>1].telo) is_telo = 1;
        }
        if((kv <= max_ext) && (!is_telo)) kv_push(uint64_t, *b, (((uint64_t)kv)<<32)|v);  
    }

    radix_sort_srt64(b->a, b->a + b->n);
    for (k = 0; k < b->n; k++) {
        v = (uint32_t)(b->a[k]);
        if (g->a[v>>1].del) continue;

        av = usg_arc_a(g, v^1); nv = usg_arc_n(g, v^1);
        for (i = kv = 0; i < nv; i++) {
            if (av[i].del) continue;
            kv++; break;
        }

        if(kv) continue;
        pb = b->n; kv_push(uint64_t, *b, v); 
        is_telo = 0; if(g->a[v>>1].telo) is_telo = 1;
        for (i = 0, w = v, kv = g->a[v>>1].occ; i < max_ext; i++) {
            if(usg_end(g, w^1, &lw)!=0) break;
            w = (uint32_t)lw; kv += g->a[w>>1].occ; kv_push(uint64_t, *b, lw);
            if(g->a[w>>1].telo) is_telo = 1;
        }
        // if((v>>1) == 308) {
        //     fprintf(stderr, "[M::%s::] v>>1::%u, v&1::%u, kv::%u\n", 
        //     __func__, v>>1, v&1, kv);
        // }


        if((kv <= max_ext) && (!is_telo)) {
            ff = 0;
            if(!ignore_ul) {///consider UL
                for (i = pb; i + 1 < b->n; i++) {
                    p = get_usg_arc(g, ((uint32_t)b->a[i]), ((uint32_t)b->a[i+1])); assert(p);
                    if(p->ou > 1) {//ignore ou == 1
                        ff = 1;
                        break;
                    }
                }

                if(ff == 0 && i < b->n) {
                    av = usg_arc_a(g, ((uint32_t)b->a[i])); nv = usg_arc_n(g, ((uint32_t)b->a[i]));
                    for (i = 0; i < nv; i++) {
                        if(av[i].del) continue;
                        if(av[i].ou > 1) {//ignore ou == 1
                            ff = 1;
                            break;
                        }
                    }
                }
            }

            // if((v>>1) == 139131) {
            //     fprintf(stderr, "+[M::%s::] v>>1::%u, v&1::%u, kv::%u, ff::%u, pb::%u, b->n::%u, w>>1::%u, w&1::%u\n", __func__, v>>1, v&1, kv,
            //     ff, pb, (uint32_t)b->n, ((uint32_t)b->a[b->n-1])>>1, ((uint32_t)b->a[b->n-1])&1);
            // }

            if((ff == 0) && (pb < b->n) && (!usg_real_tip(g, ((uint32_t)b->a[b->n-1]), 0.75))) {
                ff = 1;
            }

            // if((v>>1) == 139131) {
            //     fprintf(stderr, "-[M::%s::] v>>1::%u, v&1::%u, kv::%u, ff::%u, pb::%u, b->n::%u, w>>1::%u, w&1::%u\n", __func__, v>>1, v&1, kv,
            //     ff, pb, (uint32_t)b->n, ((uint32_t)b->a[b->n-1])>>1, ((uint32_t)b->a[b->n-1])&1);
            // }

            if(ff == 0) {
                for (i = pb; i < b->n; i++) usg_seq_del(g, ((uint32_t)b->a[i])>>1);
                cnt++;
            }

            // if((v>>1) == 139131) {
            //     fprintf(stderr, ">[M::%s::] v>>1::%u, v&1::%u, kv::%u, ff::%u, pb::%u, b->n::%u, w>>1::%u, w&1::%u, del::%u\n", __func__, v>>1, v&1, kv,
            //     ff, pb, (uint32_t)b->n, ((uint32_t)b->a[b->n-1])>>1, ((uint32_t)b->a[b->n-1])&1, 
            //     g->a[((uint32_t)b->a[b->n-1])>>1].del);
            // }
        }
        b->n = pb; 
    }

    if(!in) free(tx.a);
    if(cnt > 0) usg_cleanup(g);

    return cnt;
}

///check if v has only one branch

int32_t usg_tip_detect(usg_t *g, uint32_t v0, int max_ext, uint8_t *f, asg64_v *z_a, asg64_v *z_b)
{
    uint64_t a_n = z_a->n, b_n = z_b->n, v = v0, kv, nv, i; 
    int32_t n_ext = 0; usg_arc_t *av;

    av = usg_arc_a(g, v^1); nv = usg_arc_n(g, v^1);
    for (i = kv = 0; i < nv && kv <= 1; i++) {
        if (av[i].del || f[av[i].v>>1]) continue;
        kv++; 
    }
    if(kv > 1) return 0;
    n_ext += g->a[v>>1].occ; f[v>>1] = 1; kv_push(uint64_t, *z_b, v>>1);

    av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
    for (i = 0; i < nv; i++) {
        if (av[i].del || f[av[i].v>>1]) continue;
        kv_push(uint64_t, *z_a, av[i].v);
    }   

    while (z_a->n > a_n && n_ext < max_ext) {
        v = z_a->a[--z_a->n]; if(f[v>>1]) continue;

        av = usg_arc_a(g, v^1); nv = usg_arc_n(g, v^1);
        for (i = kv = 0; i < nv && kv < 1; i++) {
            if (av[i].del || f[av[i].v>>1]) continue;
            kv++; 
        }
        if(kv > 0) continue;
        n_ext += g->a[v>>1].occ; f[v>>1] = 1; kv_push(uint64_t, *z_b, v>>1);

        av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
        for (i = 0; i < nv; i++) {
            if (av[i].del || f[av[i].v>>1]) continue;
            kv_push(uint64_t, *z_a, av[i].v);
        }        
    }

    for (i = b_n; i < z_b->n; i++) f[z_b->a[i]] = 0;
    z_b->n = b_n; z_a->n = a_n;
    return n_ext;
}

int usg_naive_topocut_aux(usg_t *g, uint32_t v, int max_ext, uint8_t *f, asg64_v *b0, asg64_v *b1)
{
    int32_t n_ext; usg_arc_t *av; uint32_t w = v, v0 = v, nv, i, kv, tip; 
    for (n_ext = tip = 0; n_ext < max_ext; v = w) {
        tip = 0;
        av = usg_arc_a(g, v^1); nv = usg_arc_n(g, v^1);
        for (i = kv = 0; i < nv && kv <= 1; i++) {
            if (av[i].del) continue;
            kv++; 
        }
        if(kv!=1) break;
        n_ext += g->a[v>>1].occ;

        av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
        for (i = kv = 0; i < nv && kv <= 1; i++) {
            if (av[i].del) continue;
            kv++; w = av[i].v;
        }
        if(kv!=1) {
            if(kv > 1) tip = 1;
            break;
        }
    }

    if(n_ext < max_ext && tip) {
        n_ext = usg_tip_detect(g, v0, max_ext, f, b0, b1);
    }
	return n_ext;
}


int usg_naive_topocut_aux_sec(usg_t *g, uint32_t v0, int max_ext)
{
    int32_t n_ext = 0; usg_arc_t *av; uint32_t w = (uint32_t)-1, v = v0, nv, i, kv, tip[2] = {0}; 
    v = v0;
    while (1) {
        av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
        for (i = kv = 0; i < nv; i++) {
            if (av[i].del) continue;
            kv++; w = av[i].v;
        }
        n_ext += g->a[v>>1].occ;
        // if((v0>>1) == 306) {
        //     fprintf(stderr, "[M::%s::] v>>1::%u, v&1::%u, n_ext::%d\n", __func__, v>>1, v&1, n_ext);
        // }
        if(kv != 1) {
            if(kv == 0) tip[0] = 1;
            break;
        }
        v = w;

        av = usg_arc_a(g, v^1); nv = usg_arc_n(g, v^1);
        for (i = kv = 0; i < nv; i++) {
            if (av[i].del) continue;
            kv++; 
        }
        if(kv != 1) break;
        if(v == v0) return 0;///circle, it is ok to remove it
    }
    
    v = v0^1;
    while (1) {
        av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
        for (i = kv = 0; i < nv; i++) {
            if (av[i].del) continue;
            kv++; w = av[i].v;
        }
        n_ext += g->a[v>>1].occ;
        // if((v0>>1) == 306) {
        //     fprintf(stderr, "[M::%s::] v>>1::%u, v&1::%u, n_ext::%d\n", __func__, v>>1, v&1, n_ext);
        // }
        if(kv != 1) {
            if(kv == 0) tip[1] = 1;
            break;
        }
        v = w;

        av = usg_arc_a(g, v^1); nv = usg_arc_n(g, v^1);
        for (i = kv = 0; i < nv; i++) {
            if (av[i].del) continue;
            kv++; 
        }
        if(kv != 1) break;
        if(v == (v0^1)) return 0;///circle, it is ok to remove it
    }

    n_ext -= g->a[v0>>1].occ;
    // if((v0>>1) == 306) {
    //     fprintf(stderr, "[M::%s::] n_ext::%d, tip[0]::%u, tip[1]::%u\n", __func__, n_ext, tip[0], tip[1]);
    // }
    if((tip[0] || tip[1]) && (n_ext >= max_ext)) return n_ext;

    return 0;
}

void usg_arc_cut_length(usg_t *g, asg64_v *in_0, asg64_v *in_1, int32_t max_ext, float len_rat, uint32_t is_trio, 
uint32_t is_topo, uint32_t *max_drop_len)
{
    // if(len_rat > 0.7) {
    //     fprintf(stderr, "+[M::%s::] max_ext::%d, len_rat::%f\n", __func__, max_ext, len_rat);
    // }
    asg64_v tx = {0,0,0}, tz = {0,0,0}, *b = NULL, *ub = NULL;
    uint32_t i, k, v, w, n_vtx = g->n<<1, nv, nw, kv, kw, /**trioF = (uint32_t)-1, ntrioF = (uint32_t)-1,**/ ol_max, ou_max, to_del, cnt = 0, mm_ol;
    usg_arc_t *av, *aw, *ve, *we; uint64_t x, kocc[2], ou; uint8_t *f; CALLOC(f, g->n);
    b = ((in_0)?(in_0):(&tx)); ub = ((in_1)?(in_1):(&tz));
    
    for (v = 0, b->n = ub->n = 0; v < n_vtx; ++v) {
        if (g->a[v>>1].del) continue;
        av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
        if (nv < 2) continue;
        for (i = kv = kocc[0] = kocc[1] = 0; i < nv; ++i) {
            if(av[i].del) continue;
            kv++; 
            if((av[i].ou>>1) > 0) { ///if av[i].ou == 1, ignore it
                if(kocc[1] < (av[i].ou>>1)) kocc[1] = (av[i].ou>>1);
            } else if(av[i].ou == 0) { 
                kocc[0]++;
            }            
        }
        if(kv < 2 || kocc[0] == 0 || kocc[1] == 0) continue;
        ou = kocc[1];
        for (i = 0; i < nv; ++i) {
            if(av[i].del || av[i].ou) continue;
            if(max_drop_len && av[i].ol >= (*max_drop_len)) continue;
            x = (((uint64_t)av[i].ol)*10)/ou; x <<= 32;
            kv_push(uint64_t, *b, ((x)|((uint64_t)(ub->n))));
            kv_push(uint64_t, *ub, ((((uint64_t)(v))<<32)|((uint64_t)(i))));   
        }
    }
    // if(len_rat > 0.7) {
    //     fprintf(stderr, "[M::%s::] max_ext::%d, len_rat::%f, b->n::%u, ub->n::%u\n", 
    //     __func__, max_ext, len_rat, (uint32_t)b->n, (uint32_t)ub->n);
    // }

    radix_sort_srt64(b->a, b->a + b->n);
    for (k = 0; k < b->n; k++) {
        v = ub->a[(uint32_t)b->a[k]]>>32; 
        ve = &(usg_arc_a(g, v)[(uint32_t)(ub->a[(uint32_t)b->a[k]])]);
        w = ve->v^1;
        if(ve->del || g->a[v>>1].del || g->a[w>>1].del || ve->ou) continue;
        nv = usg_arc_n(g, v); nw = usg_arc_n(g, w);
        av = usg_arc_a(g, v); aw = usg_arc_a(g, w);
        if(nv<=1 && nw <= 1) continue;    

        // if(is_trio) {
        //     if(get_arcs(g, v, NULL, 0)<=1 && get_arcs(g, w, NULL, 0)<=1) continue;///speedup
        //     trioF = get_tip_trio_infor(g, v^1);
        //     ntrioF = (trioF==FATHER? MOTHER : (trioF==MOTHER? FATHER : (uint32_t)-1));
        // }    
        for (i = 0; i < nw; ++i) {
            if (aw[i].v == (v^1)) {
                we = &(aw[i]);
                break;
            }
        }
        mm_ol = MIN(ve->ol, we->ol); kocc[0] = kocc[1] = 0;

        for (i = kv = ol_max = ou_max = 0; i < nv; ++i) {
            if(av[i].del) continue;
            kv++; 
            if(av[i].ou != 1) kocc[!!(av[i].ou)]++; ///if av[i].ou == 1, ignore it
            // if(is_trio && get_tip_trio_infor(g, av[i].v) == ntrioF) continue;
            if(ol_max < av[i].ol) ol_max = av[i].ol;
        }
        if (kv < 1 || kocc[0] < 1 || kocc[1] < 1) continue;
        if (kv >= 2) {
            if (mm_ol > ol_max*len_rat) continue;
        }

        
        for (i = kw = ol_max = ou_max = 0; i < nw; ++i) {
            if(aw[i].del) continue;
            kw++;
            // if(is_trio && get_tip_trio_infor(g, aw[i].v) == ntrioF) continue;
            if(ol_max < aw[i].ol) ol_max = aw[i].ol;
        }
        if (kw < 1) continue;
        if (kw >= 2) {
            if (mm_ol > ol_max*len_rat) continue;
        }

        if (kv <= 1 && kw <= 1) continue;
        // if(len_rat > 0.7) {
        //     fprintf(stderr, "0[M::%s::] v>>1::%u(%c), w>>1::%u(%c), kv::%u, kw::%u\n", 
        //                             __func__, v>>1, "+-"[v&1], w>>1, "+-"[w&1], kv, kw);
        // }

        to_del = 1;
        if(is_topo) {
            to_del = 0;
            if (kv > 1 && kw > 1) {
                to_del = 1;
            } else if (kw == 1) {
                if (usg_naive_topocut_aux(g, w^1, max_ext, f, b, ub) < max_ext) to_del = 1;
            } else if (kv == 1) {
                if (usg_naive_topocut_aux(g, v^1, max_ext, f, b, ub) < max_ext) to_del = 1;
            }
        } 
        // if(len_rat > 0.7) {
        //     fprintf(stderr, "1[M::%s::] v>>1::%u(%c), w>>1::%u(%c), kv::%u, kw::%u, to_del::%u\n", 
        //                             __func__, v>>1, "+-"[v&1], w>>1, "+-"[w&1], kv, kw, to_del);
        // }

        if (to_del) {
            ve->del = we->del = 1; 
            if((usg_naive_topocut_aux_sec(g, v, max_ext) < max_ext) && 
                        (usg_naive_topocut_aux_sec(g, w, max_ext) < max_ext)) {
                // if((((v>>1) == 308) && ((w>>1) == 311)) || (((w>>1) == 308) && ((v>>1) == 311))) {
                //     fprintf(stderr, "[M::%s::] v>>1::%u, v&1::%u, kv::%u, w>>1::%u, w&1::%u, kw::%u, max_ext::%d\n", 
                //     __func__, v>>1, v&1, kv, w>>1, w&1, kw, max_ext);
                // }
                // if((((v>>1) == 306) && ((w>>1) == 310)) || (((w>>1) == 306) && ((v>>1) == 310))) {
                //     fprintf(stderr, "[M::%s::] v>>1::%u, v&1::%u, kv::%u, w>>1::%u, w&1::%u, kw::%u, max_ext::%d\n", 
                //     __func__, v>>1, v&1, kv, w>>1, w&1, kw, max_ext);
                // }
                ++cnt;
            } else {
                ve->del = we->del = 0; 
            }
            
        }
        // if(len_rat > 0.7) {
        //     fprintf(stderr, "2[M::%s::] v>>1::%u(%c), w>>1::%u(%c), kv::%u, kw::%u, to_del::%u\n", 
        //                             __func__, v>>1, "+-"[v&1], w>>1, "+-"[w&1], kv, kw, to_del);
        // }
    }

    if(in_0) free(tx.a); if(in_1) free(tz.a);
    if (cnt > 0) usg_cleanup(g);
    free(f);
    // fprintf(stderr, "-[M::%s::] max_ext::%d, len_rat::%f\n", __func__, max_ext, len_rat);
}

void usg_arc_cut_srt_length(usg_t *g, asg64_v *in_0, asg64_v *in_1, int32_t max_ext, float len_rat, uint32_t is_trio, 
uint32_t is_topo, uint32_t *max_drop_len, uint8_t *ff)
{
    asg64_v tx = {0,0,0}, tz = {0,0,0}, *b = NULL, *ub = NULL;
    uint32_t i, k, v, w, n_vtx = g->n<<1, nv, nw, kv, kw, /**trioF = (uint32_t)-1, ntrioF = (uint32_t)-1,**/ ol_max, ou_max, to_del, cnt = 0, mm_ol;
    usg_arc_t *av, *aw, *ve, *we; uint64_t x, kocc[2]; uint8_t *f; CALLOC(f, g->n);
    b = ((in_0)?(in_0):(&tx)); ub = ((in_1)?(in_1):(&tz));
    
    for (v = 0, b->n = ub->n = 0; v < n_vtx; ++v) {
        if (g->a[v>>1].del || ff[v]) continue;
        av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
        if (nv < 2) continue;
        for (i = kv = kocc[0] = kocc[1] = 0; i < nv; ++i) {
            if(av[i].del) continue;
            kv++; 
            if((av[i].ou>>1) > 0) { ///if av[i].ou == 1, ignore it
                if(kocc[1] < (av[i].ou>>1)) kocc[1] = (av[i].ou>>1);
            } else if(av[i].ou == 0) { 
                kocc[0]++;
            }            
        }
        if(kv < 2 || kocc[0] == 0 || kocc[1] == 0) continue;
        for (i = 0; i < nv; ++i) {
            if(av[i].del || av[i].ou) continue;
            if(max_drop_len && av[i].ol >= (*max_drop_len)) continue;
            x = av[i].ol; x <<= 32;
            kv_push(uint64_t, *b, ((x)|((uint64_t)(ub->n))));
            kv_push(uint64_t, *ub, ((((uint64_t)(v))<<32)|((uint64_t)(i))));   
        }
    }

    radix_sort_srt64(b->a, b->a + b->n);
    for (k = 0; k < b->n; k++) {
        v = ub->a[(uint32_t)b->a[k]]>>32; 
        ve = &(usg_arc_a(g, v)[(uint32_t)(ub->a[(uint32_t)b->a[k]])]);
        w = ve->v^1;
        if(ve->del || g->a[v>>1].del || g->a[w>>1].del || ve->ou) continue;
        nv = usg_arc_n(g, v); nw = usg_arc_n(g, w);
        av = usg_arc_a(g, v); aw = usg_arc_a(g, w);
        if(nv<=1 && nw <= 1) continue;    

        // if(is_trio) {
        //     if(get_arcs(g, v, NULL, 0)<=1 && get_arcs(g, w, NULL, 0)<=1) continue;///speedup
        //     trioF = get_tip_trio_infor(g, v^1);
        //     ntrioF = (trioF==FATHER? MOTHER : (trioF==MOTHER? FATHER : (uint32_t)-1));
        // }    
        for (i = 0; i < nw; ++i) {
            if (aw[i].v == (v^1)) {
                we = &(aw[i]);
                break;
            }
        }
        mm_ol = MIN(ve->ol, we->ol); kocc[0] = kocc[1] = 0;

        for (i = kv = ol_max = ou_max = 0; i < nv; ++i) {
            if(av[i].del) continue;
            kv++; 
            if(av[i].ou != 1) kocc[!!(av[i].ou)]++; ///if av[i].ou == 1, ignore it
            // if(is_trio && get_tip_trio_infor(g, av[i].v) == ntrioF) continue;
            if(ol_max < av[i].ol) ol_max = av[i].ol;
        }
        if (kv < 1 || kocc[0] < 1 || kocc[1] < 1) continue;
        if (kv >= 2) {
            if (mm_ol > ol_max*len_rat) continue;
        }

        
        for (i = kw = ol_max = ou_max = 0; i < nw; ++i) {
            if(aw[i].del) continue;
            kw++;
            // if(is_trio && get_tip_trio_infor(g, aw[i].v) == ntrioF) continue;
            if(ol_max < aw[i].ol) ol_max = aw[i].ol;
        }
        if (kw < 1) continue;
        if (kw >= 2) {
            if (mm_ol > ol_max*len_rat) continue;
        }

        if (kv <= 1 && kw <= 1) continue;
        // if(len_rat > 0.7) {
        //     fprintf(stderr, "0[M::%s::] v>>1::%u(%c), w>>1::%u(%c), kv::%u, kw::%u\n", 
        //                             __func__, v>>1, "+-"[v&1], w>>1, "+-"[w&1], kv, kw);
        // }

        to_del = 1;
        if(is_topo) {
            to_del = 0;
            if (kv > 1 && kw > 1) {
                to_del = 1;
            } else if (kw == 1) {
                if (usg_naive_topocut_aux(g, w^1, max_ext, f, b, ub) < max_ext) to_del = 1;
            } else if (kv == 1) {
                if (usg_naive_topocut_aux(g, v^1, max_ext, f, b, ub) < max_ext) to_del = 1;
            }
        } 
        // if(len_rat > 0.7) {
        //     fprintf(stderr, "1[M::%s::] v>>1::%u(%c), w>>1::%u(%c), kv::%u, kw::%u, to_del::%u\n", 
        //                             __func__, v>>1, "+-"[v&1], w>>1, "+-"[w&1], kv, kw, to_del);
        // }

        if (to_del) {
            ve->del = we->del = 1; 
            if((usg_naive_topocut_aux_sec(g, v, max_ext) < max_ext) && 
                        (usg_naive_topocut_aux_sec(g, w, max_ext) < max_ext)) {
                // if((((v>>1) == 308) && ((w>>1) == 311)) || (((w>>1) == 308) && ((v>>1) == 311))) {
                //     fprintf(stderr, "[M::%s::] v>>1::%u, v&1::%u, kv::%u, w>>1::%u, w&1::%u, kw::%u, max_ext::%d\n", 
                //     __func__, v>>1, v&1, kv, w>>1, w&1, kw, max_ext);
                // }
                // if((((v>>1) == 306) && ((w>>1) == 310)) || (((w>>1) == 306) && ((v>>1) == 310))) {
                //     fprintf(stderr, "[M::%s::] v>>1::%u, v&1::%u, kv::%u, w>>1::%u, w&1::%u, kw::%u, max_ext::%d\n", 
                //     __func__, v>>1, v&1, kv, w>>1, w&1, kw, max_ext);
                // }
                ++cnt;
            } else {
                ve->del = we->del = 0; 
            }
            
        }
        // if(len_rat > 0.7) {
        //     fprintf(stderr, "2[M::%s::] v>>1::%u(%c), w>>1::%u(%c), kv::%u, kw::%u, to_del::%u\n", 
        //                             __func__, v>>1, "+-"[v&1], w>>1, "+-"[w&1], kv, kw, to_del);
        // }
    }

    if(in_0) free(tx.a); if(in_1) free(tz.a);
    if (cnt > 0) usg_cleanup(g);
    free(f);
    // fprintf(stderr, "-[M::%s::] max_ext::%d, len_rat::%f\n", __func__, max_ext, len_rat);
}


inline int undel_arcs(usg_t *g, uint32_t v, uint32_t* v_s)
{
    uint32_t i, nv = usg_arc_n(g, v), kv; 
    usg_arc_t *av = usg_arc_a(g, v);
    for (i = kv = 0; i < nv; i++) {
        if(av[i].del) continue;
        if(v_s) v_s[kv] = av[i].v;
        kv++; 
    }
    return kv;
}

inline uint32_t get_usg_unitig(usg_t *g, uint32_t begNode, uint32_t* endNode, 
uint64_t* nodeLen, uint64_t* baseLen, uint64_t *occ, asg64_v* b)
{
    uint32_t v = begNode, w, k; usg_arc_t *av;
    uint32_t nv, kv, return_flag;
    if(endNode) (*endNode) = (uint32_t)-1;
    if(nodeLen) (*nodeLen) = 0;
    if(baseLen) (*baseLen) = 0;
    if(occ) (*occ) = 0;

    while (1) {
        kv = undel_arcs(g, v, NULL);
        if(endNode) (*endNode) = v; 
        if(nodeLen) (*nodeLen) += g->a[v>>1].occ;
		
        if(b) kv_push(uint64_t, *b, v);
        if(occ) (*occ)++;
		///means reach the end of a unitig
		if(kv!=1 && baseLen) (*baseLen) += g->a[v>>1].len;
		if(kv==0) {
			return_flag = END_TIPS; break;
		} 
		if(kv>1) {
			return_flag = MUL_OUTPUT; break;
		}
        ///kv must be 1 here
        kv = undel_arcs(g, v, &w);
        ///means reach the end of a unitig
        if(undel_arcs(g, w^1, NULL)!=1) {
            if(baseLen) (*baseLen) += g->a[v>>1].len;
            return_flag = MUL_INPUT; break;
        } else if(baseLen) {
            av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
            for (k = 0; k < nv; k++) {
                if(av[k].del) continue;
                ///here is just one undeleted edge
                (*baseLen) += asg_arc_len(av[k]);
                break;
            }
        }
        v = w;
        if(v == begNode){
			return_flag = LOOP; break;
		} 
    }

	return return_flag;
}

uint64_t dfs_max_bub(usg_t *g, buf_t *b, uint32_t x, asg64_v *nb, uint32_t *p_bub)
{
    uint64_t len = 0, baseLen, uLen; uint32_t c_v, e_v, nv, convex, v, i, kv_0, kv_1, flag_0 = 0, flag_1 = 0, op;
    usg_arc_t *av = NULL; (*p_bub) = 0; b->S.n = 0;
    if(b->a[x>>1].s || g->a[x>>1].del) return 0;
    kv_push(uint32_t, b->S, x);
    // fprintf(stderr, "\n[M::%s::] g->n::%u, x::%u\n", __func__, (uint32_t)g->n, x);
    while (b->S.n > 0) {
        c_v = b->S.a[--b->S.n];
        // fprintf(stderr, "[M::%s::] b->S.n::%u, c_v::%u\n", __func__, (uint32_t)b->S.n, c_v);
        // if(c_v >= g->n) {
        //     fprintf(stderr, "+++++[M::%s::] g->n::%u, c_v::%u\n", __func__, (uint32_t)g->n, c_v);
        // }
        if(b->a[c_v>>1].s) continue;

        nb->n = 0; op = get_usg_unitig(g, c_v, &convex, NULL, &baseLen, NULL, nb);
        uLen = baseLen;
        for(i = 0; i < nb->n; i++) b->a[nb->a[i]>>1].s = 1;
        if(op == LOOP) return 0;

        e_v = convex^1;
        nb->n = 0; op = get_usg_unitig(g, e_v, &convex, NULL, &baseLen, NULL, nb);

        uLen = MAX(uLen, baseLen); len += uLen;


        v = c_v^1; nv = usg_arc_n(g, v); av = usg_arc_a(g, v);
        for (i = kv_0 = 0; i < nv; i++) {
            if(av[i].del) continue;
            kv_0++;
            if(b->a[av[i].v>>1].s) continue;
            kv_push(uint32_t, b->S, av[i].v);
        }

        v = e_v^1; nv = usg_arc_n(g, v); av = usg_arc_a(g, v);
        for (i = kv_1 = 0; i < nv; i++) {
            if(av[i].del) continue;
            kv_1++;
            if(b->a[av[i].v>>1].s) continue;
            kv_push(uint32_t, b->S, av[i].v);
        }

        if(kv_0 > 0 && kv_1 > 0) flag_0++;
        if(kv_0 > 1) flag_1++;
        if(kv_1 > 1) flag_1++;
    }

    if(flag_0 > 0 && flag_1 > 1) (*p_bub) = 1;
    return len;
}

uint64_t usg_max_bub(usg_t *g, buf_t *b, asg64_v *nb)
{
    usg_arc_t *av, *aw; uint64_t cLen = 0, mLen = 0;
    uint32_t n_vtx = g->n<<1, k, v, w, kv, nv, kw, nw, p_bub;
    for (v = 0; v < n_vtx; ++v) {
        if(b->a[v>>1].s || g->a[v>>1].del) continue;

        av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
        for (k = kv = 0; k < nv && kv <= 1; k++) {
            if(av[k].del) continue;
            w = av[k].v^1; kv++;
        }

        if(kv == 1) {
            aw = usg_arc_a(g, w); nw = usg_arc_n(g, w);
            for (k = kw = 0; k < nw && kw <= 1; k++) {
                if(aw[k].del) continue;
                kw++;
            }
            if(kw == 1) continue;
        }

        cLen = dfs_max_bub(g, b, v^1, nb, &p_bub);
        if(p_bub == 0) continue;///no bubble
        if(cLen > mLen) mLen = cLen;
    }

    for (k = 0; k < g->n; ++k) b->a[k].s = 0;

    b->S.n = b->b.n = 0;
    return mLen;
}

uint64_t usg_bub_pop1(usg_t *g, uint32_t v0, uint64_t max_dist, buf_t *x)
{
    uint32_t v, w, i, nv, kw, cnt = 0, fail_b = 0, n_tips = 0, tip_end = (uint32_t)-1;
    uint32_t l, d, c, n_pending = 0, z, to_replace, wc; usg_arc_t *av; binfo_t *t;
    if(g->a[v0>>1].del || undel_arcs(g, v0, NULL) < 2) return 0; // already deleted

    x->S.n = x->T.n = x->b.n = x->e.n = 0;
    x->a[v0].c = x->a[v0].d = x->a[v0].m = x->a[v0].nc = x->a[v0].np = 0;
    kv_push(uint32_t, x->S, v0);

    do {
        v = kv_pop(x->S); d = x->a[v].d; c = x->a[v].c; 
        nv = usg_arc_n(g, v); av = usg_arc_a(g, v); 
        for (i = 0; i < nv; ++i) {
            if (av[i].del) continue;
            w = av[i].v; t = &(x->a[w]); l = ((v == v0)?(0):((uint32_t)av[i].ul));
            if ((w>>1) == (v0>>1)) {
                fail_b = 1;
                break;
            }
            // kv_push(uint32_t, x->e, (g->idx[v]>>32) + i); ///for backtracking
            if (d + l > max_dist) {
                fail_b = 1;
                break;
            }
            kw = undel_arcs(g, w^1, NULL); wc = g->a[w>>1].occ;
            if (t->s == 0) {
                kv_push(uint32_t, x->b, w); 
                t->p = v, t->s = 1, t->d = d + l;
                t->c = c + wc; 
                t->r = kw;
				++n_pending;
            } else {
                to_replace = 0;
                if((c + wc) < t->c) {
                    to_replace = 1;
                } else if(((c + wc) == t->c) && (d + l > t->d)) {
                    to_replace = 1;
                }
                if(to_replace) {
                    t->p = v; t->c = c + wc; 
                }
                if (d + l < t->d) t->d = d + l; // update dist
            }

            if (--(t->r) == 0) {
                z = undel_arcs(g, w, NULL);
                if(z > 0) {
                    kv_push(uint32_t, x->S, w);
                }
                else {
                    ///at most one tip
                    if(n_tips != 0) {
                        fail_b = 1;
                        break;
                    }
                    n_tips++; tip_end = w;
                }
				--n_pending;
            }
        }
        if(fail_b) break;
        if(n_tips == 1) {
            if(tip_end != (uint32_t)-1 && n_pending == 0 && x->S.n == 0) {
                kv_push(uint32_t, x->S, tip_end);
                break;
            }
            fail_b = 1;
            break;
        }

        if (i < nv || x->S.n == 0) {
            fail_b = 1;
            break;
        }
    } while (x->S.n > 1 || n_pending);

    if(!fail_b) {//there is a bubble
        cnt = 1; 
    }
    for (i = 0; i < x->b.n; ++i) { // clear the states of visited vertices
		t = &x->a[x->b.a[i]];
		t->s = t->c = t->d = t->m = t->nc = t->np = 0;
	}
    return cnt;
}

uint32_t get_usg_arc_mm(usg_t *g, usg_arc_t *z, usg_arc_mm_t **res)
{
    uint32_t v = z->ul>>32; (*res) = NULL;
    if(((uint32_t)z->idx) == 0) return 0;

    // fprintf(stderr, "[M::%s::] g->n::%u, v>>1::%u, v&1::%u, idx::%u, idx_n::%u, arc_mm.n::%u\n", __func__, 
    // (uint32_t)g->n, v>>1, v&1, (uint32_t)(z->idx>>32), (uint32_t)(z->idx), (uint32_t)(g->a[v>>1].arc_mm[v&1].n));

    (*res) = g->a[v>>1].arc_mm[v&1].a + (z->idx>>32);
    return ((uint32_t)z->idx);
}

uint32_t usg_arc_mm_consist(usg_t *g, usg_arc_t *v, usg_arc_t *w, uint32_t *inconsist, asg64_v *b)
{
    usg_arc_mm_t *v_a = NULL, *w_a = NULL; uint32_t v_n, w_n, v_k, w_k, occ = 0, n_occ = 0, cov; uint64_t l = 0;
    get_usg_unitig(g, (v->ul>>32)^1, &cov, NULL, NULL, &l, NULL); assert(cov == (w->ul>>32));
    
    v_n = get_usg_arc_mm(g, v, &v_a); w_n = get_usg_arc_mm(g, w, &w_a);
    // if((v->ul>>33) == 257 || (v->ul>>33) == 256) {
    //     fprintf(stderr, "****[M::%s::] v>>1::%u, v&1::%u, v->des::%u, v_n::%u, v_ou::%u, w>>1::%u, w&1::%u, w->des::%u, w_n::%u, w_ou::%u\n", 
    //     __func__, (uint32_t)(v->ul>>33), (uint32_t)(v->ul>>32)&1, v->v>>1, v_n, v->ou, cov>>1, cov&1, w->v>>1, w_n, w->ou);
    // }

    for (v_k = 0; v_k < v_n; v_k++) {
        for (w_k = 0; w_k < w_n; w_k++) {
            // if((v->ul>>33) == 257 || (v->ul>>33) == 256) {
            //     fprintf(stderr, "[M::%s::] v_a->uid::%u, v_a->off::%u, w_a->uid::%u, w_a->off::%u\n", 
            //     __func__, v_a[v_k].uid, v_a[v_k].off, w_a[w_k].uid, w_a[w_k].off);
            // }
            if(v_a[v_k].uid != w_a[w_k].uid) continue;
            if(v_a[v_k].off + l == w_a[w_k].off) {
                if(b) kv_push(uint64_t, *b, (((uint64_t)v_a[v_k].uid)<<32)|((uint64_t)v_a[v_k].off));
                occ++; 
            } else if(v_a[v_k].off == w_a[w_k].off + l) {
                if(b) kv_push(uint64_t, *b, (((uint64_t)v_a[v_k].uid)<<32)|((uint64_t)w_a[w_k].off));
                occ++;
            } else {
                n_occ++;
            }
        }
    }

    if(inconsist) (*inconsist) = n_occ;
    return occ;
}

uint32_t is_junction_circle(usg_t *g, uint32_t v, uint32_t w)
{
    uint32_t k, nv, cov, op, cov_w; usg_arc_t *av;
    av = usg_arc_a(g, v); nv = usg_arc_n(g, v);

    for (k = 0; k < nv; k++){
        if(av[k].del) continue;
        if(av[k].v == (w^1)) return 1;///circle
        op = get_usg_unitig(g, av[k].v, &cov, NULL, NULL, NULL, NULL);
        if(op == LOOP || cov == (w^1)) return 1;///circle
        if(op == MUL_INPUT) {
            undel_arcs(g, cov, &cov_w);
            if(cov_w == (w^1)) return 1;///circle
        }
    }

    return 0;
}

uint32_t get_junction_w(usg_t *g, uint32_t v0, uint32_t v1, uint32_t no_inconsist, asg64_v *res)
{
    uint32_t k, i, v[2], nv[2], ou[2], i0, i1, k0, k1, occ, ff, is_mul, nn[2], *c, ww, inconsist; 
    usg_arc_t *av[2]; uint64_t m = 0;

    if(is_junction_circle(g, v0, v1)) return (uint32_t)-1;
    // if((v0>>1) == 257 || (v1>>1) == 256) {
    //     fprintf(stderr, "[M::%s::] v0>>1::%u, v0&1::%u, v0pid::%u, v1>>1::%u, v1&1::%u, v1pid::%u\n", 
    //     __func__, v0>>1, v0&1, g->a[v0>>1].mm, v1>>1, v1&1, g->a[v1>>1].mm);
    // }

    v[0] = v0; v[1] = v1;
    nv[0] = usg_arc_n(g, v[0]); nv[1] = usg_arc_n(g, v[1]); 
    av[0] = usg_arc_a(g, v[0]); av[1] = usg_arc_a(g, v[1]); 

    i = 0;
    for (k = 0, ou[i] = 0; k < nv[i]; k++) {
        if(!av[i][k].ou) continue;
        if(av[i][k].del) continue;
        ou[i]++;
    }
    // if((v0>>1) == 257 || (v1>>1) == 256) {
    //     fprintf(stderr, "[M::%s::] ou[0]::%u\n", __func__, ou[0]);
    // }
    if(!ou[i]) return (uint32_t)-1;
    
    i = 1;
    for (k = 0, ou[i] = 0; k < nv[i]; k++) {
        if(!av[i][k].ou) continue;
        if(av[i][k].del) continue;
        ou[i]++;
    }
    // if((v0>>1) == 257 || (v1>>1) == 256) {
    //     fprintf(stderr, "[M::%s::] ou[1]::%u\n", __func__, ou[1]);
    // }
    if(!ou[i]) return (uint32_t)-1;

    //this is not ture; some UL may not be able to go through nid
    // if(ou[0] != ou[1]) return (uint32_t)-1;
    nn[0] = nn[1] = ww = 0;

    i0 = 0; i1 = 1; is_mul = 0; c = &(nn[0]);
    for (k0 = 0; k0 < nv[i0]; k0++) {
        if(!av[i0][k0].ou) continue;
        if(av[i0][k0].del) continue;
        for (k1 = 0, ff = 0; k1 < nv[i1] && ff <= 1; k1++) {
            if(!av[i1][k1].ou) continue;
            if(av[i1][k1].del) continue;
            occ = usg_arc_mm_consist(g, &(av[i0][k0]), &(av[i1][k1]), &inconsist, NULL);
            if(occ == 0) continue;
            if(no_inconsist && inconsist > 0) {
                ff = 2; break;
            }
            ff++;
        }

        if(ff > 1) {
            is_mul = 1; break;
        } else if(ff == 1) {
            (*c) += 1;
        }
    }
    // if((v0>>1) == 257 || (v1>>1) == 256) {
    //     fprintf(stderr, "[M::%s::] c[0]::%u\n", __func__, *c);
    // }
    if((*c) == 0 || is_mul) return (uint32_t)-1;


    i0 = 1; i1 = 0; is_mul = 0; c = &(nn[1]);
    for (k0 = 0; k0 < nv[i0]; k0++) {
        if(!av[i0][k0].ou) continue;
        if(av[i0][k0].del) continue;
        for (k1 = 0, ff = 0; k1 < nv[i1] && ff <= 1; k1++) {
            if(!av[i1][k1].ou) continue;
            if(av[i1][k1].del) continue;
            occ = usg_arc_mm_consist(g, &(av[i0][k0]), &(av[i1][k1]), &inconsist, NULL);
            if(occ == 0) continue;
            if(no_inconsist && inconsist > 0) {
                ff = 2; break;
            }
            ff++; ww += MIN((av[i0][k0].ou>>1), (av[i1][k1].ou>>1)); 
            m = k0; m <<= 32; m |= k1;
        }

        if(ff > 1) {
            is_mul = 1; break;
        } else if(ff == 1) {
            (*c) += 1;
            if(res) {
                kv_push(uint64_t, *res, (((uint64_t)v[i0])<<32)|(((uint64_t)(m>>32))));
                kv_push(uint64_t, *res, (((uint64_t)v[i1])<<32)|(((uint64_t)((uint32_t)m))));
            }
        }
    }
    // if((v0>>1) == 257 || (v1>>1) == 256) {
    //     fprintf(stderr, "[M::%s::] c[1]::%u\n", __func__, *c);
    // }
    if((*c) == 0 || is_mul) return (uint32_t)-1;
    assert(nn[0] == nn[1]);

    return ww;
}

void remap_gen_arcs(usg_t *g, uint32_t ov, uint32_t ow, uint32_t nv, uint32_t nw)
{
    usg_arc_t *op = NULL, *np = NULL; usg_arc_warp *sv; uint32_t k; usg_arc_mm_t *v_a, *pm; uint32_t v_n;
    sv = &(g->a[nv>>1].arc[nv&1]); kv_pushp(usg_arc_t, *sv, &np);
    op = get_usg_arc(g, ov, ow); assert(op);///must get op ater since op might be changed by kv_pushp

    np->ul = (uint32_t)op->ul; np->ul |= (((uint64_t)nv)<<32); np->v = nw;
    np->ol = op->ol; np->del = 0; np->ou = op->ou; 
    np->idx = g->a[nv>>1].arc_mm[nv&1].n; np->idx <<= 32; np->idx |= (uint32_t)op->idx;

    v_n = (uint32_t)op->idx;
    for (k = 0; k < v_n; k++) {
        kv_pushp(usg_arc_mm_t, g->a[nv>>1].arc_mm[nv&1], &pm);
        v_a = g->a[op->ul>>33].arc_mm[(op->ul>>32)&1].a + (op->idx>>32);//va might be changed
        pm->v = nw; pm->uid = v_a[k].uid; pm->off = v_a[k].off;
    }
}

void update_dual_junction(usg_t *g, uint32_t v, uint32_t v_id, uint32_t w, uint32_t w_id, asg64_v *buf)
{
    usg_arc_t *av, *aw, *z; uint32_t k, nv, nw, kv, kw, nid, pnid; usg_seq_t *s; 
    av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
    for (k = kv = 0; k < nv; k++) {
        if(av[k].del) continue;
        kv++;
    }

    aw = usg_arc_a(g, w); nw = usg_arc_n(g, w);
    for (k = kw = 0; k < nw; k++) {
        if(aw[k].del) continue;
        kw++;
    }
    assert(kv && kw); 
    assert((!av[v_id].del) && (!aw[w_id].del));
    assert(av[v_id].ou && aw[w_id].ou);
    if(kv == 1 && kw == 1) return;///no need to 
    
    buf->n = 0; pnid = g->n;
    get_usg_unitig(g, v^1, NULL, NULL, NULL, NULL, buf); 
    // if(!(buf->n > 0 && buf->a[buf->n-1] == w)) {
    //     fprintf(stderr, "[M::%s::] buf->n::%u, v>>1::%u, v&1::%u, kv::%u, w>>1::%u, w&1::%u, kw::%u\n", __func__, 
    //             (uint32_t)buf->n, v>>1, v&1, kv, w>>1, w&1, kw);
    // }
    assert(buf->n > 0 && buf->a[buf->n-1] == w);

    for (k = 0; k < buf->n; k++) {
        nid = buf->a[k]>>1; 
        kv_push(uint32_t, g->mp.a[g->a[nid].mm], g->n); s = push_usg_t_node(g, g->n); 
        s->mm = g->a[nid].mm; s->occ = g->a[nid].occ; s->len = g->a[nid].len; s->telo = g->a[nid].telo; s->del = 0;
        s->arc[0].n = s->arc[1].n = 0; s->arc_mm[0].n = s->arc_mm[1].n = 0;
    }
    
    for (k = 0; k + 1 < buf->n; k++) {
        remap_gen_arcs(g, buf->a[k], buf->a[k+1], (((pnid+k)<<1)|(buf->a[k]&1)), (((pnid+k+1)<<1)|(buf->a[k+1]&1)));
        remap_gen_arcs(g, buf->a[k+1]^1, buf->a[k]^1, (((pnid+k+1)<<1)|(buf->a[k+1]&1))^1, (((pnid+k)<<1)|(buf->a[k]&1))^1);
    }

    //av might be changed
    av = usg_arc_a(g, v); 
    remap_gen_arcs(g, av[v_id].ul>>32, av[v_id].v, (pnid<<1)|((av[v_id].ul>>32)&1), av[v_id].v);
    av = usg_arc_a(g, v); 
    remap_gen_arcs(g, av[v_id].v^1, (av[v_id].ul>>32)^1, av[v_id].v^1, ((pnid<<1)|((av[v_id].ul>>32)&1))^1);

    pnid = pnid + buf->n - 1;
    //aw might be changed
    aw = usg_arc_a(g, w); 
    remap_gen_arcs(g, aw[w_id].ul>>32, aw[w_id].v, (pnid<<1)|((aw[w_id].ul>>32)&1), aw[w_id].v);
    aw = usg_arc_a(g, w);
    remap_gen_arcs(g, aw[w_id].v^1, (aw[w_id].ul>>32)^1, aw[w_id].v^1, ((pnid<<1)|((aw[w_id].ul>>32)&1))^1);

    ///drop edges from the current node
    av[v_id].del = 1; z = get_usg_arc(g, av[v_id].v^1, (av[v_id].ul>>32)^1); z->del = 1;
    aw[w_id].del = 1; z = get_usg_arc(g, aw[w_id].v^1, (aw[w_id].ul>>32)^1); z->del = 1;
}

uint32_t u2g_n_hybrid_thread(usg_t *ng, uint32_t no_inconsist, asg64_v *in, asg64_v *buf)
{
    if(in->n < 1) return 0;
    uint32_t k, i, v, w, in_n = in->n, mm, ov[2], kv[2], cnt = 0; uint64_t *a, a_n;
    for (k = 0; k < in->n; k++) {
        v = (uint32_t)in->a[k]; in_n = in->n;
        if(get_usg_unitig(ng, v^1, &w, NULL, NULL, NULL, NULL) == LOOP) continue;
        mm = get_junction_w(ng, v, w, no_inconsist, in);
        if(mm == (uint32_t)-1) {
            in->n = in_n; continue;
        }
        if(in->n == in_n) continue;

        ov[0] = ov[1] = 0; a = in->a + in_n; a_n = in->n - in_n;
        for (i = 0; i < a_n; i++) {
            if((a[i]>>32) == v) ov[0]++;
            if((a[i]>>32) == w) ov[1]++;
        }

        kv[0] = undel_arcs(ng, v, NULL); kv[1] = undel_arcs(ng, w, NULL);
        assert(ov[0] == ov[1] && ov[0] <= kv[0] && ov[1] <= kv[1]); 
        if((ov[0] == kv[0] && ov[1] < kv[1]) || (ov[0] < kv[0] && ov[1] == kv[1])) {
            in->n = in_n; continue;
        }
        
        for (i = 0; i < a_n; i += 2) {
            update_dual_junction(ng, a[i]>>32, (uint32_t)a[i], a[i+1]>>32, (uint32_t)a[i+1], buf);
        }
        cnt++; in->n = in_n;
    }

    return cnt;
}


void u2g_hybrid_extend(usg_t *ng, uint64_t* i_max_dist, asg64_v *in, asg64_v *ib)
{
    uint32_t v, w, n_vtx = ng->n<<1, n_arc, nv, i, mm;
    uint64_t n_pop = 0, max_dist; usg_arc_t *av = NULL; 
    asg64_v tx = {0,0,0}, tb = {0,0,0}, *ob = NULL, *ub = NULL; 
    ob = (in?(in):(&tx)); ub = (ib?(ib):(&tb)); ob->n = ub->n = 0;
    buf_t b; memset(&b, 0, sizeof(buf_t)); CALLOC(b.a, n_vtx);
    if(i_max_dist) max_dist = (*i_max_dist);
    else max_dist = usg_max_bub(ng, &b, ob);
    uint8_t* bs_flag = NULL; CALLOC(bs_flag, n_vtx);

    if(max_dist > 0) {
        for(v = 0; v < n_vtx; ++v) {
            if(bs_flag[v] != 0) continue;
            nv = usg_arc_n(ng, v); av = usg_arc_a(ng, v);
            if(nv < 2 || ng->a[v>>1].del) continue;
            for (i = n_arc = 0; i < nv; ++i) {
                if (!av[i].del) ++n_arc;
            }
            if (n_arc < 2) continue;
            if(usg_bub_pop1(ng, v, max_dist, &b)) {
                //beg is v, end is b.S.a[0]; note b.b include end, does not include beg
                for (i = 0; i < b.b.n; i++) {
                    if(b.b.a[i]==v || b.b.a[i]==b.S.a[0]) continue;
                    bs_flag[b.b.a[i]] = bs_flag[b.b.a[i]^1] = 1;
                }
                bs_flag[v] = 2; bs_flag[b.S.a[0]^1] = 3;
            }
        }

        for(v = ob->n = 0; v < n_vtx; ++v) {
            if(bs_flag[v] <= 1) continue;
            if(get_usg_unitig(ng, v^1, &w, NULL, NULL, NULL, NULL) != LOOP && bs_flag[w] > 1) {
                mm = get_junction_w(ng, v, w, 0, NULL);
                if(mm == (uint32_t)-1) continue;
                mm = ((uint32_t)-1) - mm;
                kv_push(uint64_t, *ob, ((((uint64_t)mm)<<32)|((uint64_t)v)));
                bs_flag[v] = bs_flag[w] = 0;
            }
        }

        radix_sort_srt64(ob->a, ob->a + ob->n); 
        n_pop += u2g_n_hybrid_thread(ng, 0, ob, ub);
    }

    free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
    if(n_pop) usg_cleanup(ng);
    if(!in) free(tx.a); if(!ib) free(tb.a);
}



// uint32_t usg_path_pop_1(ul_resolve_t *uidx, usg_t *g, uint32_t v0, buf_t *x, asg64_v *arc_b)
// {
//     uint32_t v, w, i, nv, kv, kw, cnt = 0, fail_b = 0, n_tips = 0, tip_end = (uint32_t)-1;
//     uint32_t l, d, c, n_pending = 0, z, to_replace, wc, i_id, i_off, i_rev; 
//     usg_arc_t *av; binfo_t *t; usg_arc_mm_t *arc_a; uint32_t arc_n;
//     if(g->a[v0>>1].del) return 0; // already deleted
//     if(undel_arcs(g, v0, NULL) != 1) return 0;
//     undel_arcs(g, v0, &w); w ^= 1;
//     if(undel_arcs(g, w, NULL) < 2) return 0;

//     x->S.n = x->T.n = x->b.n = x->e.n = 0;
//     x->a[v0].c = x->a[v0].d = x->a[v0].m = x->a[v0].nc = x->a[v0].np = 0;
//     kv_push(uint32_t, x->S, v0); arc_b->n = 0; i_id = i_off = i_rev = (uint32_t)-1;

//     do {
//         v = kv_pop(x->S); d = x->a[v].d; c = x->a[v].c; 
//         nv = usg_arc_n(g, v); av = usg_arc_a(g, v); kv = undel_arcs(g, v, NULL);
//         for (i = 0; i < nv; ++i) {
//             if (av[i].del) continue;
//             w = av[i].v; t = &(x->a[w]); l = ((v == v0)?(0):((uint32_t)av[i].ul));
//             if ((w>>1) == (v0>>1)) {
//                 fail_b = 1;
//                 break;
//             }
//             kw = undel_arcs(g, w^1, NULL); wc = g->a[w>>1].occ;
//             kv_push(uint64_t, *arc_b, (((uint64_t)v)<<32)|((uint64_t)i)); ///for backtracking
//             arc_n = get_usg_arc_mm(g, &(av[i]), &arc_a);
//             if(v == v0 && arc_n == 0) {
//                 fail_b = 1;
//                 break;
//             }
//             if(arc_n > 0) {

//             } else if(kv == 1) {

//             }





//             kw = undel_arcs(g, w^1, NULL); wc = g->a[w>>1].occ;
//             if (t->s == 0) {
//                 kv_push(uint32_t, x->b, w); 
//                 t->p = v, t->s = 1, t->d = d + l;
//                 t->c = c + wc; 
//                 t->r = kw;
//                 ++n_pending;
//             } else {
//                 to_replace = 0;
//                 if((c + wc) < t->c) {
//                     to_replace = 1;
//                 } else if(((c + wc) == t->c) && (d + l > t->d)) {
//                     to_replace = 1;
//                 }
//                 if(to_replace) {
//                     t->p = v; t->c = c + wc; 
//                 }
//                 if (d + l < t->d) t->d = d + l; // update dist
//             }

//             if (--(t->r) == 0) {
//                 z = undel_arcs(g, w, NULL);
//                 if(z > 0) {
//                     kv_push(uint32_t, x->S, w);
//                 }
//                 else {
//                     ///at most one tip
//                     if(n_tips != 0) {
//                         fail_b = 1;
//                         break;
//                     }
//                     n_tips++; tip_end = w;
//                 }
//                 --n_pending;
//             }
//         }
//         if(fail_b) break;
//         if(n_tips == 1) {
//             if(tip_end != (uint32_t)-1 && n_pending == 0 && x->S.n == 0) {
//                 kv_push(uint32_t, x->S, tip_end);
//                 break;
//             }
//             fail_b = 1;
//             break;
//         }

//         if (i < nv || x->S.n == 0) {
//             fail_b = 1;
//             break;
//         }
//     } while (x->S.n > 1 || n_pending);
// }


///***debug-hybrid***
uint32_t check_hybrid_connect(usg_t *ng, uint32_t i_uid, uint32_t v, uint32_t vidx, uint32_t w, uint32_t widx)
{
    usg_arc_mm_t *z_a = NULL; uint32_t z_n, k;
    usg_arc_t *z = get_usg_arc(ng, v, w);
    if(!z) return 0;
    z_n = get_usg_arc_mm(ng, z, &z_a);
    if(!z_n) return 0;
    for (k = 0; k < z_n; k++) {
        if(z_a[k].uid != i_uid || z_a[k].off != vidx) continue;
        return 1;
    }
    return 0;
}

uint32_t usg_arc_occ(usg_t *ng, uint32_t v)
{
    usg_arc_t *p; uint32_t nv, kv, i;
    p = usg_arc_a(ng, v); nv = usg_arc_n(ng, v);
    for (i = kv = 0; i < nv; i++) {
        if(p[i].del) continue;
        kv++;
    }
    return kv;
}

void integer_realign_g(ul_resolve_t *uidx, usg_t *ng, uinfo_srt_warp_t *seq, uint32_t seq_id, integer_t *buf)
{
    if(seq->n < 2) return;
    // if(seq_id != 1137) return;
    uint64_t *srt, *track, k, z, m, t, i, l, seq_n, j, sc, csc, mm_sc, mm_idx, n_v, n_u, n_v0, *p, lin, pn;
    uint32_t vi, vj; mmap_t *zm, *zt; 
    ///ng->map: the nodes in the new graph that are mapped to the initial HiFi graph 
    for (i = seq_n = 0; i < seq->n; ++i) seq_n += ng->mp.a[seq->a[i].v>>1].n;
    // fprintf(stderr, "[M::%s::] seq_id::%u, seq_n::%u, seq->n::%u\n", 
    // __func__, seq_id, (uint32_t)seq_n, (uint32_t)seq->n);

    buf->o.n = buf->u.n = 0; 
    kv_resize(uint64_t, buf->o, (seq_n<<1)); 
    srt = buf->o.a; track = buf->o.a + seq_n;

    for (z = l = 0; z < seq->n; ++z) {
        csc = seq->a[z].n; 
        zm = &(ng->mp.a[seq->a[z].v>>1]);///current 
        zt = ((z>0)?(&(ng->mp.a[seq->a[z-1].v>>1])):(NULL));///prefix
        for (m = 0; m < zm->n; m++) {
            i = l + m; mm_sc = csc; mm_idx = ((uint64_t)0x7FFFFFFF); 
            vi = (zm->a[m]<<1)|(seq->a[z].v&1);
            if(zt && zt->n > 0) {
                for (t = 0; t < zt->n; t++) {
                    j = l + t - zt->n;
                    vj = (zt->a[t]<<1)|(seq->a[z-1].v&1);
                    if(check_hybrid_connect(ng, seq_id, vj, z-1, vi, z)) {///seq_id:: integer contig id
                        sc = csc + (((uint32_t)-1) - (track[j]>>32));
                        if(sc > mm_sc) {
                            mm_sc = sc; mm_idx = j;
                        }
                    }
                }
            }
            mm_sc = ((uint32_t)-1) - mm_sc;
            track[i] = mm_sc; track[i] <<= 32; track[i] |= mm_idx;
            srt[i] = track[i]>>32; srt[i] <<= 32; srt[i] |= i;
        }
        l += zm->n;
    }
    assert(l == seq_n);

    radix_sort_srt64(srt, srt + seq_n); 
    kv_resize(uint64_t, buf->res_dump, buf->res_dump.n + seq_n); 
    l = buf->res_dump.n; m = buf->u.n;
    for (k = n_v = n_u = 0; k < seq_n; k++) {
        n_v0 = n_v; i = (uint32_t)srt[k];
        for (; i != ((uint64_t)0x7FFFFFFF) && (!(track[i]&((uint64_t)0x80000000)));) {
            kv_push(uint64_t, buf->res_dump, i); n_v++; 
            track[i] |= ((uint64_t)0x80000000);
            i = track[i]&((uint64_t)0x7FFFFFFF); 
        }
        if(n_v - n_v0 <= 1) {///not useful to resolve anything if the UL covers less than 1 nodes 
            buf->res_dump.n -= (n_v-n_v0); n_v = n_v0;
            continue;
        }
        kv_pushp(uint64_t, buf->u, &p); n_u++; 
        (*p) = n_v-n_v0; (*p) = ((uint32_t)-1) - (*p); (*p) <<= 32; (*p) |= ((uint64_t)(l+n_v0));
    }

    for (z = i = 0; z < seq->n; ++z) {
        zm = &(ng->mp.a[seq->a[z].v>>1]);
        for (k = 0; k < zm->n; k++, i++) {
            srt[i] = z; srt[i] <<= 32; srt[i] |= ((zm->a[k]<<1)|(seq->a[z].v&1));
        }
    }
    assert(i == seq_n);
    ///srt[]: (idx in seq)|(node id)
    uint64_t *r, nt;
    for (k = n_u = nt = 0; k < buf->u.n; k++) {
        n_v0 = (uint32_t)buf->u.a[k]; r = buf->res_dump.a + n_v0; 
        n_v = (((uint32_t)-1) - (buf->u.a[k]>>32));
        // if(n_v < 2) continue;
        // fprintf(stderr, "+[M::%s::k->%lu] buf->u.n::%u, n_v0::%lu, n_v::%lu, nt::%lu\n", 
        // __func__, k, (uint32_t)buf->u.n, n_v0, n_v, nt);
        assert(n_v >= 2);
        buf->u.a[n_u] = nt<<32;
        for (i = 0; i < n_v; i++, nt++) {
            // fprintf(stderr, ">[M::%s::] i::%lu, nt::%lu, n_v-i-1::%lu, r[n_v-i-1]::%lu\n", 
            //         __func__, i, nt, n_v-i-1, r[n_v-i-1]);
            track[nt] = srt[r[n_v-i-1]];
        }
        buf->u.a[n_u] |= nt; n_u++;
        // fprintf(stderr, "-[M::%s::k->%lu] buf->u.n::%u, n_v0::%lu, n_v::%lu, nt::%lu\n", 
        // __func__, k, (uint32_t)buf->u.n, n_v0, n_v, nt);
    }
    buf->u.n = n_u;

    uint64_t ceq_s, ceq_e, peq_s, peq_e, min_e, max_s;
    for (k = n_u = 0; k < buf->u.n; k++) {
        n_v0 = buf->u.a[k]>>32; n_v = (uint32_t)buf->u.a[k];///[n_v0, n_v)
        ceq_s = track[n_v0]>>32; ceq_e = (track[n_v-1]>>32) + 1;
        assert((ceq_e-ceq_s) == (n_v-n_v0));
        for (z = 0; z < n_u; z++) {
            n_v0 = buf->u.a[z]>>32; n_v = (uint32_t)buf->u.a[z];
            peq_s = track[n_v0]>>32; peq_e = (track[n_v-1]>>32) + 1;
            assert((peq_e-peq_s) == (n_v-n_v0));
            max_s = MAX(peq_s, ceq_s); min_e = MIN(peq_e, ceq_e);
            if(min_e > max_s) break;
        }
        if(z < n_u) continue;
        buf->u.a[n_u++] = buf->u.a[k];
    }
    buf->u.n = n_u; buf->res_dump.n = l;
    for (k = n_u = 0; k < buf->u.n; k++) {
        // t = (uint32_t)-1; t <<= 32; 
        t = seq_id; t <<= 32; t |= ((uint64_t)0x8000000000000000);
        t |= (((uint32_t)buf->u.a[k])-(buf->u.a[k]>>32)); 
        pn = buf->res_dump.n; kv_push(uint64_t, buf->res_dump, t);
        n_v0 = buf->u.a[k]>>32; n_v = (uint32_t)buf->u.a[k];
        for (z = n_v0, lin = 1; z < n_v; z++) {
            if(lin && (z > n_v0)) {
                if((usg_arc_occ(ng, ((uint32_t)track[z])^1) > 1) || 
                            (usg_arc_occ(ng, ((uint32_t)track[z-1])) > 1)) {
                    lin = 0;
                }
            }
            kv_push(uint64_t, buf->res_dump, (uint32_t)track[z]);
        }

        if(lin == 1) buf->res_dump.n = pn;
    }
}

static void worker_integer_realign_g(void *data, long i, int tid) // callback for kt_for()
{
    ul_resolve_t *uidx = (ul_resolve_t *)data; 
    integer_realign_g(uidx, uidx->uovl.h_usg, &(uidx->uovl.cc.iug_a[i]), i, &(uidx->str_b.buf[tid]));
}

uint64_t get_ug_occ_v(uint32_t i_ug_occ)
{
    uint64_t v = (uint64_t)-1;
    if(!(i_ug_occ&((uint32_t)0x80000000))) return v;
    if(!(i_ug_occ&3)) return v;
    if((i_ug_occ&1)&&(i_ug_occ&2)) {
        v = (((i_ug_occ<<1)>>3)<<1);
        v <<= 32; v |= ((((i_ug_occ<<1)>>3)<<1)|1);
        return v;
    }
    if(i_ug_occ&1) {
        v <<= 32; v |= (((i_ug_occ<<1)>>3)<<1);
        return v;
    }
    if(i_ug_occ&2) {
        v <<= 32; v |= ((((i_ug_occ<<1)>>3)<<1)|1);
        return v;
    }
    return v;
}

///this function might be wrong
uint32_t usg_unique_arcs_cluster(asg64_v *b64, uint64_t a_n, uint64_t *idx, uint64_t *integ_seq)
{
    uint64_t bn = b64->n, k, v; 
    kv_resize(uint64_t, *b64, b64->n + a_n); b64->n += a_n; 
    memset(b64->a + bn, -1, sizeof((*(b64->a)))*a_n);
    uint64_t *cidx = b64->a + bn, s, e, i, zs, ze, z;

    ///b64.a[0, a_n]:: all resolvable paths with unique beg && end
    ///b64.a[a_n, bn]:: (raw unitig/non-unqiue node id)|(resolvable path id)
    ///idx:: the idx for b64.a[a_n, bn]
    for (i = 0; i < a_n; i++) {///available interval with beg/end with unique arcs
        s = b64->a[i]>>32; e = (uint32_t)(b64->a[i]); assert(e > s);
        if(cidx[i] != (uint64_t)-1) continue;
        for (k = s + 1; k < e; k++) {///note: here is [s, e]
            v = integ_seq[k];///v is the raw unitig id
            zs = (idx[v]<<1)>>33; ze = (uint32_t)idx[v];
            for (z = zs; z < ze; z++) {
                // if((b64->a[z]>>32)!=v) {
                //     fprintf(stderr, "[M::%s::] v::%lu, b64->a[z]::%lu, zs::%lu, ze::%lu, a_n::%lu\n", 
                //                                     __func__, v, b64->a[z]>>32, zs, ze, a_n);
                // }
                ///(b64->a[z]>>32):: raw unitig id
                ///(uint32_t)b64->a[z]:: available interval id
                assert((b64->a[z]>>32)==v);
                cidx[((uint32_t)b64->a[z])] = (i<<32)|((uint32_t)b64->a[z]);
            }

            v = integ_seq[k]^1;
            zs = (idx[v]<<1)>>33; ze = (uint32_t)idx[v];
            for (z = zs; z < ze; z++) {
                assert((b64->a[z]>>32)==v);
                cidx[((uint32_t)b64->a[z])] = (i<<32)|((uint32_t)b64->a[z]);
            }
        }
        v = integ_seq[s];
        zs = (idx[v]<<1)>>33; ze = (uint32_t)idx[v];
        for (z = zs; z < ze; z++) {
            assert((b64->a[z]>>32)==v);
            cidx[((uint32_t)b64->a[z])] = (i<<32)|((uint32_t)b64->a[z]);
        }

        v = integ_seq[e]^1;
        zs = (idx[v]<<1)>>33; ze = (uint32_t)idx[v];
        for (z = zs; z < ze; z++) {
            assert((b64->a[z]>>32)==v);
            cidx[((uint32_t)b64->a[z])] = (i<<32)|((uint32_t)b64->a[z]);
        }
    }

    radix_sort_srt64(cidx, cidx + a_n); b64->n = a_n;
    for (i = 0, k = 1, v = 0; k <= a_n; k++) {
        if(k == a_n || (cidx[i]>>32) != (cidx[k]>>32)) {
            for (z = i; z < k; z++) {
                assert(cidx[z] != (uint64_t)-1);
                ///cluest integer seqs-> (cluster id)|(available interval id) 
                b64->a[b64->n++] = (v<<32)|((uint32_t)cidx[z]);
            }
            i = k; v++;
        }
    }

    assert(b64->n == (a_n<<1));
    return v;///how many cluster
}

void iter_unique_arcs(asg64_v *buf, asg64_v *b64, uint64_t a_n, uint64_t *idx, uint64_t *integ_seq, uint64_t *cidx, uint64_t i0)
{
    uint64_t s, e, x, k, v, zs, ze, z;
    buf->n = 0;
    kv_push(uint64_t, *buf, i0);
    while(buf->n) {
        x = kv_pop(*buf);
        if(cidx[x] != (uint64_t)-1) continue;
        cidx[x] = (i0<<32)|(x);
        s = b64->a[x]>>32; e = (uint32_t)(b64->a[x]); assert(e > s);
        for (k = s + 1; k < e; k++) {///note: here is [s, e]
            v = integ_seq[k];///v is the raw unitig id
            zs = (idx[v]<<1)>>33; ze = (uint32_t)idx[v];
            for (z = zs; z < ze; z++) {
                ///(b64->a[z]>>32):: raw unitig id
                ///(uint32_t)b64->a[z]:: available interval id
                assert((b64->a[z]>>32)==v);
                if(cidx[((uint32_t)b64->a[z])] != (uint64_t)-1) {
                    assert((cidx[((uint32_t)b64->a[z])]>>32)==i0);
                    continue;
                }
                // cidx[((uint32_t)b64->a[z])] = (i0<<32)|((uint32_t)b64->a[z]);
                kv_push(uint64_t, *buf, ((uint32_t)b64->a[z]));
            }

            v = integ_seq[k]^1;
            zs = (idx[v]<<1)>>33; ze = (uint32_t)idx[v];
            for (z = zs; z < ze; z++) {
                ///(b64->a[z]>>32):: raw unitig id
                ///(uint32_t)b64->a[z]:: available interval id
                assert((b64->a[z]>>32)==v);
                if(cidx[((uint32_t)b64->a[z])] != (uint64_t)-1) {
                    assert((cidx[((uint32_t)b64->a[z])]>>32)==i0);
                    continue;
                }
                // cidx[((uint32_t)b64->a[z])] = (i0<<32)|((uint32_t)b64->a[z]);
                kv_push(uint64_t, *buf, ((uint32_t)b64->a[z]));
            }
        }

        v = integ_seq[s];
        zs = (idx[v]<<1)>>33; ze = (uint32_t)idx[v];
        for (z = zs; z < ze; z++) {
            assert((b64->a[z]>>32)==v);
            if(cidx[((uint32_t)b64->a[z])] != (uint64_t)-1) {
                assert((cidx[((uint32_t)b64->a[z])]>>32)==i0);
                continue;
            }
            // cidx[((uint32_t)b64->a[z])] = (i0<<32)|((uint32_t)b64->a[z]);
            kv_push(uint64_t, *buf, ((uint32_t)b64->a[z]));
        }

        v = integ_seq[e]^1;
        zs = (idx[v]<<1)>>33; ze = (uint32_t)idx[v];
        for (z = zs; z < ze; z++) {
            assert((b64->a[z]>>32)==v);
            if(cidx[((uint32_t)b64->a[z])] != (uint64_t)-1) {
                assert((cidx[((uint32_t)b64->a[z])]>>32)==i0);
                continue;
            }
            // cidx[((uint32_t)b64->a[z])] = (i0<<32)|((uint32_t)b64->a[z]);
            kv_push(uint64_t, *buf, ((uint32_t)b64->a[z]));
        }
    }
}

uint32_t usg_unique_arcs_cluster_adv(asg64_v *b64, uint64_t a_n, uint64_t *idx, uint64_t *integ_seq, asg64_v *buf)
{
    uint64_t bn = b64->n, k, z, v; buf->n = 0;
    kv_resize(uint64_t, *b64, b64->n + a_n); b64->n += a_n; 
    memset(b64->a + bn, -1, sizeof((*(b64->a)))*a_n);
    uint64_t *cidx = b64->a + bn, i;

    ///b64.a[0, a_n]:: all resolvable paths with unique beg && end
    ///b64.a[a_n, bn]:: (raw unitig/non-unqiue node id)|(resolvable path id)
    ///idx:: the idx for b64.a[a_n, bn]
    for (i = 0; i < a_n; i++) {///available interval with beg/end with unique arcs
        if(cidx[i] != (uint64_t)-1) continue;
        iter_unique_arcs(buf, b64, a_n, idx, integ_seq, cidx, i);
    }

    radix_sort_srt64(cidx, cidx + a_n); b64->n = a_n;
    for (i = 0, k = 1, v = 0; k <= a_n; k++) {
        if(k == a_n || (cidx[i]>>32) != (cidx[k]>>32)) {
            for (z = i; z < k; z++) {
                assert(cidx[z] != (uint64_t)-1);
                ///cluest integer seqs-> (cluster id)|(available interval id) 
                b64->a[b64->n++] = (v<<32)|((uint32_t)cidx[z]);
            }
            i = k; v++;
        }
    }

    buf->n = 0;
    assert(b64->n == (a_n<<1));
    return v;///how many cluster
}

uint32_t ava_pass_unique_bridge(uint64_t *idx, uint64_t *integer_seq, uint64_t s, uint64_t e)
{
    uint64_t k;
    for (k = s + 1; k < e; k++) {///note: here is [s, e]
        if((idx[integer_seq[k]]&((uint64_t)0x8000000000000000)) || 
                        (idx[integer_seq[k]^1]&((uint64_t)0x8000000000000000))) {
            return 0;
        }
    }
    if((idx[integer_seq[s]]&((uint64_t)0x8000000000000000)) || 
                        (idx[integer_seq[e]^1]&((uint64_t)0x8000000000000000))) {
        return 0;
    }
    return 1;
}

uint32_t ava_pass_unique_bridge_tips(usg_t *g, asg64_v *b64, uint64_t g_s, uint64_t g_e, uint64_t *integer_seq, uint8_t *f, uint64_t max_ext, double max_ext_rate, uint64_t ext_up)
{
    uint64_t i, k, z, s, e, v, nv, bn = b64->n, kv, kv_t, n_ext = 0, nkeep = 0, ncut = 0, is_telo; usg_arc_t *av;
    for (i = g_s; i < g_e; i++) {///available intervals within the same cluster
        s = b64->a[((uint32_t)b64->a[i])]>>32; e = ((uint32_t)b64->a[((uint32_t)b64->a[i])]); assert(s < e);
        for (k = s + 1; k < e; k++) {///note: here is [s, e]
            f[integer_seq[k]] = f[integer_seq[k]^1] = 1; 
            nkeep += g->a[integer_seq[k]>>1].occ;
        }
        f[integer_seq[s]] = f[integer_seq[e]^1] = 1;
    }
    nkeep = nkeep*max_ext_rate;

    ///collect nodes within raw unitig graph that are linked by the clusters but not in the cluster
    for (i = g_s; i < g_e; i++) {
        s = b64->a[((uint32_t)b64->a[i])]>>32; e = ((uint32_t)b64->a[((uint32_t)b64->a[i])]); assert(s < e);
        for (k = s + 1; k < e; k++) {///note: here is [s, e]
            v = integer_seq[k]; av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
            for (z = 0; z < nv; z++) {
                if(av[z].del || f[av[z].v^1]) continue;
                kv_push(uint64_t, *b64, av[z].v);
            }

            v = integer_seq[k]^1; av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
            for (z = 0; z < nv; z++) {
                if(av[z].del || f[av[z].v^1]) continue;
                kv_push(uint64_t, *b64, av[z].v);
            }
        }

        v = integer_seq[s]; av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
        for (z = 0; z < nv; z++) {
            if(av[z].del || f[av[z].v^1]) continue;
            kv_push(uint64_t, *b64, av[z].v);
        }

        v = integer_seq[e]^1; av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
        for (z = 0; z < nv; z++) {
            if(av[z].del || f[av[z].v^1]) continue;
            kv_push(uint64_t, *b64, av[z].v);
        }
    }

    ncut = MAX(nkeep, max_ext); is_telo = 0;
    if(ncut > ext_up) ncut = ext_up;
    if(ncut < max_ext) ncut = max_ext;
    while ((b64->n > bn) && (n_ext < ncut) && (!is_telo)) {
        v = b64->a[--b64->n]; if(f[v]) continue;
        av = usg_arc_a(g, v^1); nv = usg_arc_n(g, v^1);
        for (i = kv = 0; i < nv && kv < 1; i++) {
            if (av[i].del || f[av[i].v^1]) continue;
            kv++; 
        }
        if(kv > 0) continue;
        n_ext += g->a[v>>1].occ; f[v] = f[v^1] = 1;

        av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
        for (i = kv_t = 0; i < nv; i++) {
            if(av[i].del) continue;
            kv_t++;
            if (f[av[i].v^1]) continue;
            kv_push(uint64_t, *b64, av[i].v);
        }
        if((!kv_t) && (g->a[v>>1].telo)) is_telo = 1;     
    }

    ///reset f[] to 0
    for (i = g_s; i < g_e; i++) { 
        s = b64->a[((uint32_t)b64->a[i])]>>32; e = ((uint32_t)b64->a[((uint32_t)b64->a[i])]); assert(s < e);
        for (k = s + 1; k < e; k++) {///note: here is [s, e]
            kv_push(uint64_t, *b64, integer_seq[k]); 
            kv_push(uint64_t, *b64, (integer_seq[k]^1));
            f[integer_seq[k]] = f[integer_seq[k]^1] = 0; 
        }
        kv_push(uint64_t, *b64, integer_seq[s]); 
        kv_push(uint64_t, *b64, (integer_seq[e]^1));
        f[integer_seq[s]] = f[integer_seq[e]^1] = 0;
    }

    while (b64->n > bn) {
        v = b64->a[--b64->n]; 
        f[v] = f[v^1] = 0;

        av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
        for (i = 0; i < nv; i++) {
            if (av[i].del) continue;
            if(f[av[i].v] || f[av[i].v^1]) {
                kv_push(uint64_t, *b64, av[i].v);
            }
        }

        av = usg_arc_a(g, v^1); nv = usg_arc_n(g, v^1);
        for (i = 0; i < nv; i++) {
            if (av[i].del) continue;
            if(f[av[i].v] || f[av[i].v^1]) {
                kv_push(uint64_t, *b64, av[i].v);
            }
        }
    }
    if((n_ext < ncut) && (n_ext > max_ext - 1)) n_ext = max_ext - 1;
    if(is_telo) n_ext = max_ext + 1;
    return n_ext;
}

void debug_sysm_usg_t(usg_t *ng, const char *cmd) 
{
    uint32_t v, i; usg_arc_t *p, *q; uint32_t nv;
    for (v = 0; v < (ng->n<<1); v++) {
        p = usg_arc_a(ng, v); nv = usg_arc_n(ng, v);
        for (i = 0; i < nv; i++) {
            q = get_usg_arc(ng, p[i].v^1, v^1);
            if((!q) || (p[i].del != q->del)) {
                fprintf(stderr, "[M::%s::%s] utg%.6dl(%c)%u -> utg%.6dl(%c)%u, p->del::%u, q->del::%u\n", 
                __func__, cmd, 
                (int32_t)(v>>1)+1, "+-"[v&1], v, 
                (int32_t)(p[i].v>>1)+1, "+-"[p[i].v&1], p[i].v, 
                p[i].del, q?q->del:1);
            }
        }
    }
}

void update_usg_t_threading_0(ul_resolve_t *uidx, usg_t *ng, uint64_t *a, uint64_t a_n, uint32_t *occ, asg64_v *b)
{
    uint64_t k, i, *p, nid, nnid, v, w; b->n = 0; usg_seq_t *s; 

    // fprintf(stderr, "[M::%s::] a_n::%lu\n", __func__, a_n);
    // for (k = 0; k < a_n; k++) {
    //     fprintf(stderr, "utg%.6dl(%c)(hom::%u)\t", (int32_t)(a[k]>>1)+1, "+-"[a[k]&1], (IF_HOM((a[k]>>1), (*(uidx->bub)))));
    // }
    // fprintf(stderr, "\n");

    assert(a_n > 1); 
    assert(occ[a[0]] == 1); 
    kv_pushp(uint64_t, *b, &p); (*p) = a[0]; (*p) <<= 32; (*p) |= a[0]; occ[a[0]]--; 
    
    for (k = 1; k + 1 < a_n; k++) {
        assert(occ[a[k]] == occ[a[k]^1]); assert(occ[a[k]] > 0);
        if(occ[a[k]] > 1) {///copy a node
            nid = a[k]>>1; nnid = ng->n;
            kv_push(uint32_t, ng->mp.a[ng->a[nid].mm], nnid);
            s = push_usg_t_node(ng, nnid); 
            s->mm = ng->a[nid].mm; s->occ = ng->a[nid].occ; s->len = ng->a[nid].len; s->telo = ng->a[nid].telo; s->del = 0;
            s->arc[0].n = s->arc[1].n = 0; s->arc_mm[0].n = s->arc_mm[1].n = 0;
            nnid <<= 1; nnid |= a[k]&1;
            kv_pushp(uint64_t, *b, &p); (*p) = a[k]^1; (*p) <<= 32; (*p) |= nnid^1;
            kv_pushp(uint64_t, *b, &p); (*p) = a[k]; (*p) <<= 32; (*p) |= nnid;
        } else {
            kv_pushp(uint64_t, *b, &p); (*p) = a[k]^1; (*p) <<= 32; (*p) |= a[k]^1;
            kv_pushp(uint64_t, *b, &p); (*p) = a[k]; (*p) <<= 32; (*p) |= a[k];
        }
        occ[a[k]]--; occ[a[k]^1]--;
    }
    // if(occ[a[a_n-1]^1] != 1) {
    //     fprintf(stderr, "[M::%s] utg%.6dl(%c)(occ::%u)\n", __func__,
    //             (int32_t)((a[a_n-1]^1)>>1)+1, "+-"[(a[a_n-1]^1)&1], occ[a[a_n-1]^1]);
    // }
    assert(occ[a[a_n-1]^1] == 1); 
    kv_pushp(uint64_t, *b, &p); (*p) = a[a_n-1]^1; (*p) <<= 32; (*p) |= a[a_n-1]^1; occ[a[a_n-1]^1]--; 
    
    for (k = 0; k < b->n; k += 2) {
        v = b->a[k]; w = b->a[k+1];
        if(((v>>32) == ((uint32_t)v)) && ((w>>32) == ((uint32_t)w))) continue;//no need arc
        remap_gen_arcs(ng, v>>32, (w>>32)^1, ((uint32_t)v), ((uint32_t)w)^1);///v>>32: old id; ((uint32_t)v): new id
        remap_gen_arcs(ng, w>>32, (v>>32)^1, ((uint32_t)w), ((uint32_t)v)^1);
    }

    // usg_arc_t *z = get_usg_arc(ng, 220, 157075), *q = get_usg_arc(ng, 157074, 221);
    // fprintf(stderr, "\n+[M::%s::a_n->%lu] p->del::%u, q->del::%u\n", 
    //             __func__, a_n, z?z->del:1, q?q->del:1);

    usg_arc_t *av, *aw, *z; uint32_t kv, kw, nv, nw, v0, v1, w0, w1;
    for (k = 0; k < b->n; k += 2) {
        v0 = b->a[k]>>32; v1 = (uint32_t)b->a[k];
        w0 = b->a[k+1]>>32; w1 = (uint32_t)b->a[k+1];

        av = usg_arc_a(ng, v1); nv = usg_arc_n(ng, v1);
        for (i = kv = 0; i < nv; i++) {
            if(av[i].v == (w1^1)) {
                av[i].del = 0; kv++;
            } else if(av[i].del == 0) {
                av[i].del = 1;
                z = get_usg_arc(ng, av[i].v^1, (av[i].ul>>32)^1); 
                // if(!(z && (!z->del))) {
                //     fprintf(stderr, "[M::%s::] z::%u, z->del::%u\n", __func__, z?1:0, z?z->del:1);
                //     fprintf(stderr, "[M::%s::arc] utg%.6dl(%c)->utg%.6dl(%c)\n", 
                //                                     __func__, (int32_t)(av[i].ul>>33)+1, "+-"[(av[i].ul>>32)&1], 
                //                                               (int32_t)(av[i].v>>1)+1, "+-"[(av[i].v)&1]);
                //     fprintf(stderr, "[M::%s::new] utg%.6dl(%c)->utg%.6dl(%c)\n", 
                //                                     __func__, (int32_t)(v1>>1)+1, "+-"[v1&1], 
                //                                               (int32_t)(w1>>1)+1, "+-"[w1&1]);
                //     fprintf(stderr, "[M::%s::old] utg%.6dl(%c)->utg%.6dl(%c)\n", 
                //                                     __func__, (int32_t)(v0>>1)+1, "+-"[v0&1], 
                //                                               (int32_t)(w0>>1)+1, "+-"[w0&1]);
                // }
                assert(z && (!z->del)); z->del = 1;
            }
        }
        assert(kv == 1);
        
        aw = usg_arc_a(ng, w1); nw = usg_arc_n(ng, w1);
        for (i = kw = 0; i < nw; i++) {
            if(aw[i].v == (v1^1)) {
                aw[i].del = 0; kw++;
            } else if(aw[i].del == 0) {
                aw[i].del = 1;
                z = get_usg_arc(ng, aw[i].v^1, (aw[i].ul>>32)^1); 
                assert(z && (!z->del)); z->del = 1;
            }
        }
        assert(kw == 1);

        if(v0 == v1 && w0 == w1) continue;

        av = usg_arc_a(ng, v0); nv = usg_arc_n(ng, v0);
        for (i = 0; i < nv; i++) {
            if(av[i].v == (w0^1)) {
                av[i].del = 1; break;
            }
        }
        assert(i < nv);

        aw = usg_arc_a(ng, w0); nw = usg_arc_n(ng, w0);
        for (i = 0; i < nw; i++) {
            if(aw[i].v == (v0^1)) {
                aw[i].del = 1; break;
            }
        }
        assert(i < nw);
    }

    // z = get_usg_arc(ng, 220, 157075), q = get_usg_arc(ng, 157074, 221);
    // fprintf(stderr, "-[M::%s::a_n->%lu] p->del::%u, q->del::%u\n", 
    //             __func__, a_n, z?z->del:1, q?q->del:1);
    
}

void update_usg_t_threading(ul_resolve_t *uidx, usg_t *ng, uint64_t *arcs, uint64_t *arcs_g, uint64_t arcs_gn, uint64_t *integ_seq, uint32_t *occ, asg64_v *b)
{
    uint64_t i, k, s, e, nvtx = ng->n<<1; memset(occ, 0, sizeof((*occ))*nvtx);
    for (i = 0; i < arcs_gn; i++) {///set cluster
        // prt = 0;
        s = arcs[((uint32_t)arcs_g[i])]>>32; e = ((uint32_t)arcs[((uint32_t)arcs_g[i])]); assert(s < e);
        for (k = s + 1; k < e; k++) {///note: here is [s, e]
            occ[integ_seq[k]]++; occ[integ_seq[k]^1]++; 
            // if((integ_seq[k]>>1) == 2736) prt = 1;
        }
        occ[integ_seq[s]]++; occ[integ_seq[e]^1]++;
        // if((integ_seq[s]>>1) == 2736) prt = 1;
        // if(((integ_seq[e]^1)>>1) == 2736) prt = 1;
        // if(occ[integ_seq[s]] > 1 || occ[integ_seq[e]^1] > 1) {
        //     fprintf(stderr, "s::utg%.6dl(%c)(occ::%u), e::utg%.6dl(%c)(occ::%u)\n", 
        //         (int32_t)(integ_seq[s]>>1)+1, "+-"[integ_seq[s]&1], occ[integ_seq[s]],
        //         (int32_t)((integ_seq[e]^1)>>1)+1, "+-"[((integ_seq[e]^1)&1)], occ[integ_seq[e]^1]);
        // }
        // if(prt) {
        //     fprintf(stderr, "[M::%s::] a_n::%lu\n", __func__, e + 1 - s);
        //     for (k = s; k <= e; k++) {
        //         fprintf(stderr, "utg%.6dl(%c)(hom::%u)\t", 
        //             (int32_t)(integ_seq[k]>>1)+1, "+-"[integ_seq[k]&1], 
        //             (IF_HOM((integ_seq[k]>>1), (*(uidx->bub)))));
        //     }
        //     fprintf(stderr, "\n");
        // }
    }

    for (i = 0; i < nvtx; i++) {
        if(!occ[i]) occ[i] = (uint32_t)-1;        
    }
    for (i = 0; i < arcs_gn; i++) {///arcs_g[]>>32 is the group id; (uint32_t)arcs_g[]
        s = arcs[((uint32_t)arcs_g[i])]>>32; e = ((uint32_t)arcs[((uint32_t)arcs_g[i])]); 
        update_usg_t_threading_0(uidx, ng, integ_seq + s, e + 1 - s, occ, b);
    }
}

ma_ug_t *ma_ug_hybrid_gen(usg_t *g)
{
    int32_t *mark; uint32_t i, v, n_vtx = g->n<<1;
    uint32_t w, x, l, start, end, len; ma_utg_t *p;
    kdq_t(uint64_t) *q; ///is a queue
    ma_ug_t *ug;

    ug = (ma_ug_t*)calloc(1, sizeof(ma_ug_t));
    ug->g = asg_init();
    ///each node has two directions
    mark = (int32_t*)calloc(n_vtx, 4);

    q = kdq_init(uint64_t);
    for (v = 0; v < n_vtx; ++v) {
        // fprintf(stderr, "+[M::%s::] v::%u, n_vtx::%u\n", __func__, v, n_vtx);
        if (g->a[v>>1].del || mark[v]) continue;
        if (usg_arc_n(g, v) == 0 && usg_arc_n(g, (v^1)) != 0) continue;
        // fprintf(stderr, "-[M::%s::] v::%u, n_vtx::%u\n", __func__, v, n_vtx);
        mark[v] = 1; q->count = 0, start = v, end = v^1, len = 0;
        // forward
        w = v;
        while (1) {
            /**
             * w----->x
             * w<-----x
             * that means the only suffix of w is x, and the only prefix of x is w
             **/
            if (usg_arc_n(g, w) != 1) break;
            x = usg_arc_a(g, w)[0].v; // w->x
            if (usg_arc_n(g, x^1) != 1) break;

            /**
             * another direction of w would be marked as used (since w has been used)
            **/
            mark[x] = mark[w^1] = 1;
            ///l is the edge length, instead of overlap length
            ///note: edge length is different with overlap length
            l = asg_arc_len(usg_arc_a(g, w)[0]);
            kdq_push(uint64_t, q, (uint64_t)w<<32 | l);
            end = x^1, len += l;
            w = x;
            if (x == v) break;
        }
        if (start != (end^1) || kdq_size(q) == 0) { // linear unitig
            ///length of seq, instead of edge
            l = g->a[end>>1].len;
            kdq_push(uint64_t, q, (uint64_t)(end^1)<<32 | l);
            len += l;
        } else { // circular unitig
            start = end = UINT32_MAX;
            goto add_usg_unitig; // then it is not necessary to do the backward
        }

        // backward
        x = v;
        while (1) { // similar to forward but not the same
            if (usg_arc_n(g, x^1) != 1) break;
            w = usg_arc_a(g, x^1)[0].v ^ 1; // w->x
            if (usg_arc_n(g, w) != 1) break;
            mark[x] = mark[w^1] = 1;
            l = asg_arc_len(usg_arc_a(g, w)[0]);
            ///w is the seq id + direction, l is the length of edge
            ///push element to the front of a queue
            kdq_unshift(uint64_t, q, (uint64_t)w<<32 | l);

            // fprintf(stderr, "uId: %u, >%.*s (%u)\n", 
            // ug->u.n, (int)Get_NAME_LENGTH((R_INF), w>>1), Get_NAME((R_INF), w>>1), w>>1);

            start = w, len += l;
            x = w;
        }

        add_usg_unitig:
        if (start != UINT32_MAX) mark[start] = mark[end] = 1;
        kv_pushp(ma_utg_t, ug->u, &p);
        p->s = 0, p->start = start, p->end = end, p->len = len, p->n = kdq_size(q), p->circ = (start == UINT32_MAX);
        p->m = p->n;
        kv_roundup32(p->m);
        p->a = (uint64_t*)malloc(8 * p->m);
        //all elements are saved here
        for (i = 0; i < kdq_size(q); ++i)
            p->a[i] = kdq_at(q, i);
        // fprintf(stderr, "*[M::%s::] v::%u, n_vtx::%u\n", __func__, v, n_vtx);
    }
    kdq_destroy(uint64_t, q);
    // fprintf(stderr, "-1-[M::%s::] **************\n", __func__);
    // add arcs between unitigs; reusing mark for a different purpose
    //ug saves all unitigs
    for (v = 0; v < n_vtx; ++v) mark[v] = -1;
    // fprintf(stderr, "-2-[M::%s::] **************\n", __func__);

    //mark all start nodes and end nodes of all unitigs
    for (i = 0; i < ug->u.n; ++i) {
        if (ug->u.a[i].circ) continue;
        mark[ug->u.a[i].start] = i<<1 | 0;
        mark[ug->u.a[i].end] = i<<1 | 1;
    }
    // fprintf(stderr, "-3-[M::%s::] **************\n", __func__);
    //scan all edges
    usg_arc_t *av; uint32_t nv;
    for (v = 0; v < n_vtx; v++) {
        av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
        for (i = 0; i < nv; i++) {
            if(av[i].del) continue;
            if (mark[av[i].ul>>32^1] >= 0 && mark[av[i].v] >= 0) {
                uint32_t u = mark[(av[i].ul>>32)^1]^1;
                int l = ug->u.a[u>>1].len - av[i].ol;
                if (l < 0) l = 1;
                asg_arc_t *q = asg_arc_pushp(ug->g);
                q->ol = av[i].ol, q->del = 0;
                q->ul = (uint64_t)u<<32 | l;
                q->v = mark[av[i].v]; q->ou = 0;
            }
        }
    }
    // fprintf(stderr, "-4-[M::%s::] **************\n", __func__);
    for (i = 0; i < ug->u.n; ++i)
        asg_seq_set(ug->g, i, ug->u.a[i].len, 0);
    // fprintf(stderr, "-5-[M::%s::] **************\n", __func__);
    asg_cleanup(ug->g);
    free(mark);
    return ug;
}


// ma_ug_t *gen_unique_g(ul_resolve_t *uidx, usg_t *ng, uint32_t *hg_occ, uint64_t *integ_seq_idx, uint64_t integ_seq_n, uint64_t *integ_seq_a, uint32_t max_ext)
// {
//     ma_ug_t *un_g = ma_ug_hybrid_gen(ng); ma_utg_t *u; uint64_t v, w, pi, ei, *pz, bn, a_n, ua_n;
//     uint32_t i, k, n_vtx = un_g->g->n_seq<<1, s, e, vw[2], rev, sid, arc_id, z, zs, ze, mm; 
//     uint8_t *arc_del; MALLOC(arc_del, n_vtx); asg_arc_t *p; asg64_v b64; kv_init(b64);
//     asg64_v b_za, b_zb; kv_init(b_za); kv_init(b_zb); uint8_t *ff; CALLOC(ff, ng->n<<1);
//     for (i = 0; i < un_g->u.n; i++) {
//         un_g->g->seq[i].c = PRIMARY_LABLE;
//         u = &(un_g->u.a[i]);
//         arc_del[i<<1] = arc_del[(i<<1)+1] = 1;
//         if(u->n > 0) {
//             if(hg_occ[u->a[u->n-1]>>32]==1) arc_del[i<<1] = 0;
//             if(hg_occ[(u->a[0]>>32)^1]==1) arc_del[(i<<1)+1] = 0;
//         }
//         if(arc_del[i<<1] && arc_del[(i<<1)+1]) un_g->g->seq[i].del = 1;
//     }
//     free(un_g->g->idx); un_g->g->idx = 0; un_g->g->is_srt = 0; un_g->g->n_arc = 0;///release all edges

//     for (i = 0; i < un_g->u.n; i++) {
//         if(un_g->g->seq[i].del) continue;
//         if(arc_del[i<<1] && arc_del[(i<<1)+1]) continue;
//         u = &(un_g->u.a[i]);
//         if(!arc_del[i<<1]) {
//             assert(hg_occ[u->a[u->n-1]>>32] == 1);
//             hg_occ[u->a[u->n-1]>>32] = ((uint32_t)0x80000000);
//             hg_occ[u->a[u->n-1]>>32] |= (i<<1);
//         }

//         if(!arc_del[(i<<1)+1]) {
//             assert(hg_occ[(u->a[0]>>32)^1] == 1);
//             hg_occ[(u->a[0]>>32)^1] = ((uint32_t)0x80000000);
//             hg_occ[(u->a[0]>>32)^1] |= ((i<<1)|1);
//         }
//     }

//     for (i = b64.n = ua_n = a_n = 0; i < integ_seq_n; i++) {
//         s = integ_seq_idx[i]; e = s + ((uint32_t)integ_seq_idx[i]); 
//         assert(e > s + 1);
//         for (k = s, bn = b64.n, pi = (uint64_t)-1; k < e; k++) {
//             v = integ_seq_a[k];
//             if((!(hg_occ[v]&((uint32_t)0x80000000)))&&(!(hg_occ[v^1]&((uint32_t)0x80000000)))) {
//                 continue;///it must be a unique node
//             }
//             if((pi != (uint64_t)-1) && (hg_occ[v^1]&((uint32_t)0x80000000))) {
//                 kv_pushp(uint64_t, b64, &pz); *pz = pi; (*pz) <<= 32; (*pz) |= k; a_n++;
//             }
//             if(hg_occ[v]&((uint32_t)0x80000000)) pi = k;
//         }

//         for (k = bn, pi = s; k < b64.n; k++) {///note here is [s, e]
//             ei = b64.a[k]>>32;
//             if(pi < ei) {
//                 kv_pushp(uint64_t, b64, &pz); ua_n++;
//                 (*pz) = pi; (*pz) <<= 32; (*pz) |= ei; (*pz) |= ((uint64_t)0x8000000000000000);
//             }
//             pi = (uint32_t)b64.a[k];
//         }
//         if(pi < e) {///note here is [s, e]
//             ei = e - 1;
//             if(pi < ei) {
//                 kv_pushp(uint64_t, b64, &pz); ua_n++;
//                 (*pz) = pi; (*pz) <<= 32; (*pz) |= ei; (*pz) |= ((uint64_t)0x8000000000000000);
//             }
//         }
//     }
//     assert(a_n + ua_n == b64.n);

//     uint64_t *i_idx, b64_n, iid, n_clus; CALLOC(i_idx, ng->n<<1);
//     radix_sort_srt64(b64.a, b64.a + b64.n);
//     for (i = ua_n; i < b64.n; i++) {///unavailable intervals; mask all unavailable nodes
//         assert(b64.a[i]&((uint64_t)0x8000000000000000));
//         s = (b64.a[i]<<1)>>33; e = (uint32_t)b64.a[i]; assert(e > s);//[s, e]
//         for (k = s + 1; k < e; k++) {///note: here is [s, e]
//             i_idx[integ_seq_a[k]] |= ((uint64_t)0x8000000000000000);
//             i_idx[integ_seq_a[k]^1] |= ((uint64_t)0x8000000000000000);
//         }
//         i_idx[integ_seq_a[s]] |= ((uint64_t)0x8000000000000000);
//         i_idx[integ_seq_a[e]^1] |= ((uint64_t)0x8000000000000000);
//     }
//     ///unavailable intervals are useless
//     b64.n = a_n; ua_n = 0;
//     for (i = 0; i < a_n; i++) { ///available intervals
//         s = b64.a[i]>>32; e = (uint32_t)b64.a[i]; assert(e > s);///[s, e]
//         assert(!(b64.a[i]&((uint64_t)0x8000000000000000)));
//         for (k = s + 1; k < e; k++) {///note: here is [s, e]
//             kv_pushp(uint64_t, b64, &pz); i_idx[integ_seq_a[k]]++;
//             (*pz) = integ_seq_a[k]; (*pz) <<= 32; (*pz) |= i;

//             kv_pushp(uint64_t, b64, &pz); i_idx[integ_seq_a[k]^1]++;
//             (*pz) = integ_seq_a[k]^1; (*pz) <<= 32; (*pz) |= i;
//         }
//         kv_pushp(uint64_t, b64, &pz); i_idx[integ_seq_a[s]]++;
//         (*pz) = integ_seq_a[s]; (*pz) <<= 32; (*pz) |= i;

//         kv_pushp(uint64_t, b64, &pz); i_idx[integ_seq_a[e]^1]++;
//         (*pz) = integ_seq_a[e]^1; (*pz) <<= 32; (*pz) |= i;
//     }

//     ///index
//     radix_sort_srt64(b64.a + a_n, b64.a + b64.n);
//     for (k = a_n + 1, i = a_n; k <= b64.n; k++) {
//         if(k == b64.n || (b64.a[k]>>32) != (b64.a[i]>>32)) {
//             v = b64.a[i]>>32;
//             i_idx[v] |= (((uint64_t)i)<<32)|((uint64_t)k);
//             i = k;
//         }
//     }

//     n_clus = usg_unique_arcs_cluster(&b64, a_n, i_idx, integ_seq_a);
//     assert(b64.n == (a_n<<1));
//     for (k = a_n + 1, i = a_n, mm = a_n; k <= b64.n; k++) {
//         if(k == b64.n || (b64.a[k]>>32) != (b64.a[i]>>32)) {
//             for (z = i; z < k; z++) {
//                 s = b64.a[((uint32_t)b64.a[z])]>>32; e = ((uint32_t)b64.a[((uint32_t)b64.a[z])]); assert(e > s);
//                 if(!ava_pass_unique_bridge(i_idx, integ_seq_a, s, e)) break;
//             }
//             if(z >= k) {///all arcs in this cluster is fine
//                 for (z = i; z < k; z++) b64.a[mm++] = b64.a[z];
//             }
//             i = k; n_clus--;
//         }
//     }
//     assert(n_clus == 0);
//     b64.n = mm; n_clus = 0;
//     for (k = a_n + 1, i = a_n, mm = a_n; k <= b64.n; k++) {
//         if(k == b64.n || (b64.a[k]>>32) != (b64.a[i]>>32)) {
//             if(ava_pass_unique_bridge_tips(ng, &b64, i, k, integ_seq_a, ff, max_ext) < max_ext) {///no long tip
//                 for (z = i; z < k; z++) b64.a[mm++] = b64.a[z];
//                 n_clus++;
//             }
//             i = k; 
//         }
//     }
//     b64.n = mm;

//     if(n_clus > 0) {
//         update_usg_t_threading(uidx, ng, b64.a, b64.a + a_n, b64.n - a_n, integ_seq_a, hg_occ, NULL);
//     }

//     free(arc_del); free(i_idx); kv_destroy(b64);
// }









void prt_thread_info(uint64_t *interval, uint64_t interval_n, uint64_t *cluster, uint64_t cluster_n, 
const char *nn)
{
    char* gfa_name = NULL; MALLOC(gfa_name, strlen(nn)+70);
    sprintf(gfa_name, "%s.thread_info.log", nn);
    FILE* fp = fopen(gfa_name, "w"); free(gfa_name);
    if (!fp) return;
    uint64_t k;
    if(interval) {
        for (k = 0; k < interval_n; k++) fprintf(fp,"it_val::[%lu, %u)\n", interval[k]>>32, (uint32_t)interval[k]);
    }
    if(cluster) {
        for (k = 0; k < cluster_n; k++) fprintf(fp,"cluster::%lu\tit_id::%u)\n", cluster[k]>>32, (uint32_t)cluster[k]);
    }
    fclose(fp);
}

void prt_intg_info(uint64_t *int_idx, uint64_t int_idx_n, uint64_t *int_a, const char *nn)
{
    char* gfa_name = NULL; MALLOC(gfa_name, strlen(nn)+70);
    sprintf(gfa_name, "%s.intg_info.log", nn);
    FILE* fp = fopen(gfa_name, "w"); free(gfa_name);
    if (!fp) return;
    uint64_t k, i, s, e;
    for (i = 0; i < int_idx_n; i++) {///scan all integer contigs
        s = int_idx[i]>>32; e = s + ((uint32_t)int_idx[i]); 
        assert(e > s + 1);//the length is at least 2
        fprintf(fp,"idx::[%lu, %lu)\n", s, e);
        for (k = s; k < e; k++) {
            fprintf(fp,"%lu\n", int_a[k]);
        }
    }
    fclose(fp);
}


void prt_usg_t(ul_resolve_t *uidx, usg_t *ng, const char *cmd);


#define occ_m(x) ((x)&((uint32_t)0x7fffffff))
#define c_unqiue_m(v, occ) (occ_m((occ)[(v)]) == 1 && occ_m((occ)[(v)^1]) <= 1)
#define reli_pass(x) (!((x)&((uint64_t)0x8000000000000000)))

uint64_t get_ext_tip(usg_t *g, uint64_t *seq, uint64_t s, uint64_t e, asg64_v *buf, asg64_v *set, uint8_t *f, uint64_t max_ext)
{
    uint64_t bn = buf->n, sn = set->n, v, k, nv, z, n_ext = 0, i, kv; usg_arc_t *av;
    for (k = s + 1; k < e; k++) {///note: here is [s, e]
        v = seq[k]; av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
        for (z = 0; z < nv; z++) {
            if(av[z].del || f[av[z].v^1]) continue;
            kv_push(uint64_t, *buf, av[z].v);
        }

        v = seq[k]^1; av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
        for (z = 0; z < nv; z++) {
            if(av[z].del || f[av[z].v^1]) continue;
            kv_push(uint64_t, *buf, av[z].v);
        }
    }

    v = seq[s]; av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
    for (z = 0; z < nv; z++) {
        if(av[z].del || f[av[z].v^1]) continue;
        kv_push(uint64_t, *buf, av[z].v);
    }

    v = seq[e]^1; av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
    for (z = 0; z < nv; z++) {
        if(av[z].del || f[av[z].v^1]) continue;
        kv_push(uint64_t, *buf, av[z].v);
    }


    while (buf->n > bn && n_ext < max_ext) {
        v = buf->a[--buf->n]; if(f[v]) continue;
        av = usg_arc_a(g, v^1); nv = usg_arc_n(g, v^1);
        for (i = kv = 0; i < nv && kv < 1; i++) {
            if (av[i].del || f[av[i].v^1]) continue;
            kv++; 
        }
        if(kv > 0) continue;
        n_ext += g->a[v>>1].occ; 
        if(!(f[v])) {
            f[v] = 1; kv_push(uint64_t, *set, v);
        }
        if(!(f[v^1])) {
            f[v^1] = 1; kv_push(uint64_t, *set, (v^1));
        }

        av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
        for (i = 0; i < nv; i++) {
            if (av[i].del || f[av[i].v^1]) continue;
            kv_push(uint64_t, *buf, av[i].v);
        }    
    }
    buf->n = bn;
    for (k = sn; k < set->n; k++) f[set->a[k]] = 0;
    return n_ext;
}

void select_unqiue_path_on_fly(usg_t *g, uint64_t s, uint64_t e, asg64_v *b64, asg64_v *ub64,
uint64_t *integer_seq, uint32_t *ng_occ, uint64_t *is_reli, uint8_t *f, uint64_t max_ext, asg64_v *res)
{
    uint64_t k, sk = (uint64_t)-1, ek = (uint64_t)-1, *p = NULL; ub64->n = 0;
    assert(occ_m(ng_occ[integer_seq[s]]) == 1);
    assert(occ_m(ng_occ[integer_seq[e]^1]) == 1);
    
    if(reli_pass(is_reli[integer_seq[s]])) {
        sk = s; 
    } else {
        for (k = s + 1; k < e; k++) {
            if((reli_pass(is_reli[integer_seq[k]])) && (occ_m(ng_occ[integer_seq[k]]) == 1)) {
                sk = k;
                break;
            }
        }
    }
    if((sk == ((uint64_t)-1)) || (sk >= e)) return;

    for (k = sk + 1; k <= e; k++) {
        if(sk != ((uint64_t)-1)) {
            if(reli_pass(is_reli[(integer_seq[k]^1)])) {
                if(occ_m((ng_occ[integer_seq[k]]^1)) == 1) {
                    ek = k; ///check and extend [sk, ek]
                    if((ek > sk) && (get_ext_tip(g, integer_seq, sk, ek, b64, ub64, f, max_ext) < max_ext)) {
                        p = ((res->n > 0)? (&(res->a[res->n-1])):NULL);
                        if(p && (((uint32_t)(*p)) == sk)) {
                            (*p) >>= 32; (*p) <<= 32; (*p) |= ek;
                        } else {
                            kv_pushp(uint64_t, *res, &p); *p = sk; (*p) <<= 32; (*p) |= ek;
                        }
                    }
                    sk = ek = ((uint64_t)-1);
                }
            } else {
                sk = ek = ((uint64_t)-1);
            }
        }

        ///case2: sk != (uint64_t)-1 && ek == (uint64_t)-1 -> wait for an available ek
        ///case3: sk == (uint64_t)-1 && ek == (uint64_t)-1 -> not available
        if(k < e) {
            if((!(reli_pass(is_reli[integer_seq[k]])))) {
                sk = ek = (uint64_t)-1; continue;
            }
            if((sk == (uint64_t)-1) && (occ_m(ng_occ[integer_seq[k]]) == 1)) {
                sk = k;
            }
        }
    }
}

void select_unqiue_path_on_fly_raw(usg_t *g, uint64_t s, uint64_t e, asg64_v *b64, asg64_v *ub64,
uint64_t *integer_seq, uint32_t *ng_occ, uint64_t *is_reli, uint8_t *f, uint64_t max_ext, asg64_v *res)
{
    uint64_t k, sk = (uint64_t)-1, ek = (uint64_t)-1, *p = NULL; ub64->n = 0;
    assert(occ_m(ng_occ[integer_seq[s]]) == 1);
    assert(occ_m(ng_occ[integer_seq[e]^1]) == 1);
    sk = s; ek = e;
    if(get_ext_tip(g, integer_seq, sk, ek, b64, ub64, f, max_ext) < max_ext) {
        kv_pushp(uint64_t, *res, &p); *p = sk; (*p) <<= 32; (*p) |= ek;
        // fprintf(stderr, "+ext_k::[%lu, %lu]\ts::%lu\te::%lu\n", sk, ek, s, e);
        // fprintf(stderr, "occ[sk]::%u\tocc[ek]::%u\tocc[s]::%u\tocc[e]::%u\n", 
        // occ_m(ng_occ[integer_seq[sk]]), occ_m(ng_occ[integer_seq[ek]^1]), 
        // occ_m(ng_occ[integer_seq[s]]), occ_m(ng_occ[integer_seq[e]^1]));
    }
    return;
    


    for (k = sk + 1; k <= e; k++) {
        if(sk != ((uint64_t)-1)) {
            if(reli_pass(is_reli[(integer_seq[k]^1)])) {
                if(occ_m((ng_occ[integer_seq[k]]^1)) == 1) {
                    ek = k; ///check and extend [sk, ek]
                    if((ek > sk) && (get_ext_tip(g, integer_seq, sk, ek, b64, ub64, f, max_ext) < max_ext)) {
                        p = ((res->n > 0)? (&(res->a[res->n-1])):NULL);
                        if(p && (((uint32_t)(*p)) == sk)) {
                            (*p) >>= 32; (*p) <<= 32; (*p) |= ek;
                        } else {
                            kv_pushp(uint64_t, *res, &p); *p = sk; (*p) <<= 32; (*p) |= ek;
                        }
                    }
                    sk = ek = ((uint64_t)-1);
                }
            } else {
                sk = ek = ((uint64_t)-1);
            }
        }

        ///case2: sk != (uint64_t)-1 && ek == (uint64_t)-1 -> wait for an available ek
        ///case3: sk == (uint64_t)-1 && ek == (uint64_t)-1 -> not available
        if(k < e) {
            if((!(reli_pass(is_reli[integer_seq[k]])))) {
                sk = ek = (uint64_t)-1; continue;
            }
            if((sk == (uint64_t)-1) && (occ_m(ng_occ[integer_seq[k]]) == 1)) {
                sk = k;
            }
        }
    }
}

void update_thread_path(usg_t *g, asg64_v *b64, asg64_v *ub64, uint64_t g_s, uint64_t g_e, uint64_t *integer_seq, 
uint32_t *ng_occ, uint64_t *reliable_idx, uint8_t *f, uint64_t max_ext, asg64_v *res)
{
    uint64_t i, k, s, e; 
    for (i = g_s; i < g_e; i++) {///available intervals within the same cluster
        s = b64->a[((uint32_t)b64->a[i])]>>32; e = ((uint32_t)b64->a[((uint32_t)b64->a[i])]); assert(s < e);
        for (k = s + 1; k < e; k++) {///note: here is [s, e]
            f[integer_seq[k]] = f[integer_seq[k]^1] = 1;
        }
        f[integer_seq[s]] = f[integer_seq[e]^1] = 1;
    }

    for (i = g_s; i < g_e; i++) {///available intervals within the same cluster
        s = b64->a[((uint32_t)b64->a[i])]>>32; e = ((uint32_t)b64->a[((uint32_t)b64->a[i])]); assert(s < e);
        // fprintf(stderr, "[M::%s::] i::%lu, chain_id::%u\n", __func__, i, ((uint32_t)b64->a[i]));
        // select_unqiue_path_on_fly(g, s, e, b64, ub64, integer_seq, ng_occ, reliable_idx, f, max_ext, res);
        select_unqiue_path_on_fly_raw(g, s, e, b64, ub64, integer_seq, ng_occ, reliable_idx, f, max_ext, res);
    }

    for (i = g_s; i < g_e; i++) {///available intervals within the same cluster
        s = b64->a[((uint32_t)b64->a[i])]>>32; e = ((uint32_t)b64->a[((uint32_t)b64->a[i])]); assert(s < e);
        for (k = s + 1; k < e; k++) {///note: here is [s, e]
            f[integer_seq[k]] = f[integer_seq[k]^1] = 0;
        }
        f[integer_seq[s]] = f[integer_seq[e]^1] = 0;
    }
}


typedef struct {
    ul_resolve_t *uidx;
    uint64_t *integer_seq;
    uint64_t *gidx;
    uint64_t *interval_idx;
    uint64_t gidx_n;
    uint8_t *f;
    usg_t *ng;
    uint64_t *int_idx;
    uint64_t int_idx_n; 
    uint64_t *int_a; 
    uint64_t int_an;
    uint64_t *ridx_a; 
    uint64_t *ridx;
    // uint8_t is_double_check;
} unique_bridge_check_t;

uint64_t get_arc_support(ul_resolve_t *uidx, uint64_t v, uint64_t w)
{
    uint64_t *hid_a, hid_n; ul_str_idx_t *str_idx = &(uidx->pstr);
    uint64_t z, vz, occ = 0; ul_str_t *str; int64_t s, s_n;
    hid_a = str_idx->occ.a + str_idx->idx.a[v>>1];
    hid_n = str_idx->idx.a[(v>>1)+1] - str_idx->idx.a[v>>1];

    for (z = occ = 0; z < hid_n; z++) {
        str = &(str_idx->str.a[hid_a[z]>>32]); s_n = str->cn;
        if(s_n < 2) continue;
        vz = (uint32_t)(str->a[(uint32_t)hid_a[z]]);
        assert((v>>1) == (vz>>1)); 

        if(v == vz) {
            for(s = ((uint32_t)hid_a[z]) + 1; s < s_n; s++) {
                if(((uint32_t)(str->a[s])) == w) {
                    occ++; break;
                }
            }
        } else {
            for(s = ((int32_t)((uint32_t)hid_a[z]))-1; s >= 0; s--) {
                if(((uint32_t)(str->a[s])) == (w^1)) {
                    occ++; break;
                }
            }
        }
    }
    return occ;
}
#define unique_bridge_occ 2
#define unique_bridge_rate 0.499999

// static void worker_unique_bridge_check(void *data, long i, int tid) // callback for kt_for()
// {
//     unique_bridge_check_t *uaux = (unique_bridge_check_t *)data;
//     uint64_t *integer_seq = uaux->integer_seq, s, e, v, w, nse, self_k = i, k, kv;
//     uint64_t *gidx = uaux->gidx, *interval = uaux->interval_idx;
//     uint8_t *f = uaux->f;

//     s = interval[((uint32_t)gidx[i])]>>32; v = integer_seq[s];
//     e = ((uint32_t)interval[((uint32_t)gidx[i])]); w = integer_seq[e]^1;
//     assert(s < e);

//     nse = get_arc_support(uaux->uidx, v, w^1); 
//     if(nse < unique_bridge_occ) {
//         f[v] = f[w] = 1; return;
//     }

//     for (k = 0; k < uaux->gidx_n; k++) {
//         if(k == self_k) continue;
//         s = interval[((uint32_t)gidx[k])]>>32;
//         e = ((uint32_t)interval[((uint32_t)gidx[k])]); 
//         kv = get_arc_support(uaux->uidx, v, integer_seq[s]^1);
//         if((kv > 0) && (kv >= (nse*unique_bridge_rate))) {
//             f[v] = f[w] = 1; return;
//         }

//         kv = get_arc_support(uaux->uidx, v, integer_seq[e]);
//         if((kv > 0) && (kv >= (nse*unique_bridge_rate))) {
//             f[v] = f[w] = 1; return;
//         }

//         kv = get_arc_support(uaux->uidx, w, integer_seq[s]^1);
//         if((kv > 0) && (kv >= (nse*unique_bridge_rate))) {
//             f[v] = f[w] = 1; return;
//         }

//         kv = get_arc_support(uaux->uidx, w, integer_seq[e]);
//         if((kv > 0) && (kv >= (nse*unique_bridge_rate))) {
//             f[v] = f[w] = 1; return;
//         }
//     }
// }

uint64_t get_arc_support_chain(ul_resolve_t *uidx, uint64_t *a, uint64_t a_n, usg_t *ng)
{
    if(a_n <= 0) return 0;
    uint64_t *hid_a, hid_n; ul_str_idx_t *str_idx = &(uidx->pstr);
    uint64_t z, vz, occ = 0, v, w, ai; ul_str_t *str; int64_t s, s_n;
    v = (ng->a[a[0]>>1].mm<<1)|(a[0]&1);
    hid_a = str_idx->occ.a + str_idx->idx.a[v>>1];
    hid_n = str_idx->idx.a[(v>>1)+1] - str_idx->idx.a[v>>1];

    for (z = occ = 0; z < hid_n; z++) {
        str = &(str_idx->str.a[hid_a[z]>>32]); s_n = str->cn;
        if(s_n < 2) continue;
        vz = (uint32_t)(str->a[(uint32_t)hid_a[z]]);
        assert((v>>1) == (vz>>1)); 

        if(v == vz) {
            for(s = ((uint32_t)hid_a[z]) + 1, ai = 1; s < s_n && ai < a_n; s++, ai++) {
                w = (ng->a[a[ai]>>1].mm<<1)|(a[ai]&1);
                if(((uint32_t)(str->a[s])) != w) break;
            }
            if(ai >= a_n) occ++;
        } else {
            for(s = ((int32_t)((uint32_t)hid_a[z]))-1, ai = 1; s >= 0 && ai < a_n; s--, ai++) {
                w = (ng->a[a[ai]>>1].mm<<1)|(a[ai]&1);
                if(((uint32_t)(str->a[s])) != (w^1)) break;
            }
            if(ai >= a_n) occ++;
        }
    }
    return occ;
}

// uint64_t is_consist_ul(ul_resolve_t *uidx, uint64_t *a, uint64_t a_n, usg_t *ng)
// {
//     if(a_n <= 0) return 0;
//     uint64_t k, v, w, z, *hid_a, hid_n, ulid, i; ul2ul_item_t *it;;
//     ul_str_idx_t *str_idx = &(uidx->pstr);
//     for (k = 0; k < a_n; k++) {
//         v = (ng->a[a[k]>>1].mm<<1)|(a[k]&1);
//         hid_a = str_idx->occ.a + str_idx->idx.a[v>>1];
//         hid_n = str_idx->idx.a[(v>>1)+1] - str_idx->idx.a[v>>1];
//         for (z = 0; z < hid_n; z++) {
//             ulid = hid_a[z]>>32;
//             if(ulid < uidx->uovl.uln) {
//                 it = get_ul_ovlp(&(uidx->uovl), ulid, 1);
//                 if(!it) continue;
//                 assert(uidx->pstr.str.a[ulid].cn > 1);
//                 if(it->is_consist == 0) return 0;
//             }
//         }
//     }

//     uint64_t nv; asg_arc_t *av;
//     v = (ng->a[a[0]>>1].mm<<1)|(a[0]&1);
//     for (i = 1; i < a_n; i++) {
//         w = (ng->a[a[i]>>1].mm<<1)|(a[i]&1);

//         nv = asg_arc_n(uidx->l1_ug->g, v); 
//         av = asg_arc_a(uidx->l1_ug->g, v);
//         for (k = 0; k < nv; k++) {
//             if(av[k].del) continue;
//             if(av[k].v == w) break;
//         }
//         if(k >= nv) return 0;

//         nv = asg_arc_n(uidx->l1_ug->g, w^1); 
//         av = asg_arc_a(uidx->l1_ug->g, w^1);
//         for (k = 0; k < nv; k++) {
//             if(av[k].del) continue;
//             if(av[k].v == (v^1)) break;
//         }
//         if(k >= nv) return 0;

//         v = w;
//     }

//     return 1;
// }

static void worker_unique_bridge_check_s(void *data, long i, int tid) // callback for kt_for()
{
    unique_bridge_check_t *uaux = (unique_bridge_check_t *)data;
    uint64_t *integer_seq = uaux->integer_seq, s, e, vk, wk, nse, k;
    uint64_t *gidx = uaux->gidx, *interval = uaux->interval_idx;
    bubble_type *bub = uaux->uidx->bub; usg_t *ng = uaux->ng;

    s = interval[((uint32_t)gidx[i])]>>32; 
    e = ((uint32_t)interval[((uint32_t)gidx[i])]); 
    assert(s < e);
    for (k = s, vk = wk = (uint32_t)-1; k <= e; k++) {
        if(k > s) {
            // nse = get_arc_support(uaux->uidx, (ng->a[integer_seq[k-1]>>1].mm<<1)|(integer_seq[k-1]&1), 
            // (ng->a[integer_seq[k]>>1].mm<<1)|(integer_seq[k]&1));
            nse = get_arc_support_chain(uaux->uidx, integer_seq+k-1, 2, ng);
            if(nse < unique_bridge_occ) {
                // if((!(uaux->is_double_check)) || (!is_consist_ul(uaux->uidx, integer_seq+k-1, 2, ng))) {
                    gidx[i] |= ((uint64_t)0x8000000000000000); return;
                // }
            }
        }
        if(IF_HOM((integer_seq[k]>>1), *bub)) continue;
        // w = (ng->a[integer_seq[k]>>1].mm<<1)|(integer_seq[k]&1);
        wk = k;
        if(vk != (uint32_t)-1) {
            // nse = get_arc_support(uaux->uidx, v, w);
            nse = get_arc_support_chain(uaux->uidx, integer_seq+vk, wk+1-vk, ng);
            if(nse < unique_bridge_occ) {
                // if((!(uaux->is_double_check)) || (!is_consist_ul(uaux->uidx, integer_seq+vk, wk+1-vk, ng))) {
                    gidx[i] |= ((uint64_t)0x8000000000000000); return;
                // }
            }
        }
        vk = wk;
    }
}

uint32_t ava_pass_unique_bridge_cov(ul_resolve_t *uidx, usg_t *g, asg64_v *b64, uint64_t g_s, uint64_t g_e, uint64_t *integer_seq, uint64_t tip_l)
{
    if((tip_l == 0) && (g_e - g_s == 1)) return 1; ///if only one path, go through in anyway 
    uint64_t i, is_del = 0; unique_bridge_check_t uaux;

    uaux.uidx = uidx; uaux.integer_seq = integer_seq; 
    uaux.gidx = b64->a + g_s; uaux.gidx_n = g_e - g_s;
    uaux.interval_idx = b64->a; uaux.ng = g;
    ///if tip_l == 0, it is more likely to be right, so give more chance by double checking
    // if(!tip_l) uaux.is_double_check = 1;
    // else uaux.is_double_check = 0;

    kt_for(uidx->str_b.n_thread, worker_unique_bridge_check_s, (&uaux), g_e-g_s);///seq->n > 1

    for (i = g_s; i < g_e; i++) {///available intervals within the same cluster
        if(b64->a[i]&((uint64_t)0x8000000000000000)) {
            b64->a[i] -= ((uint64_t)0x8000000000000000); is_del = 1;
        }
    }
    return (!is_del);
}



uint64_t old_path_ext(ul_resolve_t *uidx, usg_t *ng, asg64_v *b64, asg64_v *ub64, uint64_t a_n, uint64_t n_clus, uint64_t *i_idx, uint64_t *int_a, 
uint8_t *ff, uint32_t *ng_occ, uint32_t max_ext)
{
    uint64_t k, i, z, /**s, e,**/ mm, tip_l;
    /**
    for (k = a_n + 1, i = a_n, mm = a_n; k <= b64->n; k++) {
        if(k == b64->n || (b64->a[k]>>32) != (b64->a[i]>>32)) {
            for (z = i; z < k; z++) {///all intger seqs within the same cluster
                s = b64->a[((uint32_t)b64->a[z])]>>32; 
                e = ((uint32_t)b64->a[((uint32_t)b64->a[z])]); assert(e > s);
                ///[s, e]:: available interval
                if(!ava_pass_unique_bridge(i_idx, int_a, s, e)) break;
            }
            if(z >= k) {///all arcs in this cluster is fine -> each of arch is reliable
                for (z = i; z < k; z++) b64->a[mm++] = b64->a[z];
            }
            i = k; n_clus--;
        }
    }
    assert(n_clus == 0);
    b64->n = mm; 
    **/
    // prt_thread_info(b64.a, a_n, NULL, 0, "tt5");
    // fprintf(stderr, "**4**[M::%s::] a_n::%lu, n_clus::%lu\n", __func__, a_n, n_clus);
    n_clus = 0;
    for (k = a_n + 1, i = a_n, mm = a_n; k <= b64->n; k++) {
        if(k == b64->n || (b64->a[k]>>32) != (b64->a[i]>>32)) {
            tip_l = ava_pass_unique_bridge_tips(ng, b64, i, k, int_a, ff, max_ext, 0.03, 16);
            if(tip_l < max_ext) {///no long tip
                if(/**(!tip_l) || (**/ava_pass_unique_bridge_cov(uidx, ng, b64, i, k, int_a, tip_l)) {
                    for (z = i; z < k; z++) {
                        b64->a[mm++] = b64->a[z];
                    }
                    n_clus++;
                }
            }
            i = k; 
        }
    }
    // prt_thread_info(b64.a, a_n, NULL, 0, "tt6");
    b64->n = mm;
    // fprintf(stderr, "**5**[M::%s::] a_n::%lu, n_clus::%lu, ng->n::%lu\n", __func__, a_n, n_clus, ng->n);
    // prt_thread_info(b64.a, a_n, b64.a + a_n, b64.n - a_n, "thred");
    if(n_clus > 0) {
        update_usg_t_threading(uidx, ng, b64->a, b64->a + a_n, b64->n - a_n, int_a, ng_occ, ub64);
    }
    // fprintf(stderr, "**6**[M::%s::] a_n::%lu, n_clus::%lu, ng->n::%lu\n", __func__, a_n, n_clus, ng->n);
    return n_clus;
}



void new_path_ext(ul_resolve_t *uidx, usg_t *ng, asg64_v *b64, asg64_v *ub64, uint64_t a_n, uint64_t n_clus, uint64_t *i_idx, uint64_t *int_a, 
uint8_t *ff, uint32_t *ng_occ, uint32_t max_ext)
{
    asg64_v res; kv_init(res); uint64_t k, i, s, e, nvtx = ng->n<<1;
    for (k = a_n + 1, i = a_n, res.n = 0; k <= b64->n; k++) {
        if(k == b64->n || (b64->a[k]>>32) != (b64->a[i]>>32)) {
            update_thread_path(ng, b64, ub64, i, k, int_a, ng_occ, i_idx, ff, max_ext, &res);
            i = k; n_clus--;
        }
    }
    assert(n_clus == 0);
    

    if(res.n > 0) {
        memset(ng_occ, 0, sizeof((*ng_occ))*nvtx);
        for (i = 0; i < res.n; i++) {
            s = res.a[i]>>32; e = ((uint32_t)res.a[i]); assert(s < e);
            for (k = s + 1; k < e; k++) {///note: here is [s, e]
                ng_occ[int_a[k]]++; ng_occ[int_a[k]^1]++; 
            }
            ng_occ[int_a[s]]++; ng_occ[int_a[e]^1]++;
        }
        for (i = 0; i < nvtx; i++) {
            if(!ng_occ[i]) ng_occ[i] = (uint32_t)-1;        
        }

        for (i = 0, ub64->n = 0; i < res.n; i++) {
            s = res.a[i]>>32; e = ((uint32_t)res.a[i]); 
            update_usg_t_threading_0(uidx, ng, int_a + s, e + 1 - s, ng_occ, ub64);
        }
    }
    kv_destroy(res);
}

uint64_t gen_unique_g_adv_old(ul_resolve_t *uidx, usg_t *ng, uint64_t *int_idx, uint64_t int_idx_n, uint64_t *int_a, uint32_t max_ext)
{
    uint64_t pi, ei, *pz, bn, a_n, ua_n, *r_a, r_n; uint8_t *ff; CALLOC(ff, ng->n<<1);
    uint32_t i, k, n_vtx = ng->n<<1, s, e, zs, ze, v, b64_n; 
    uint32_t *ng_occ; CALLOC(ng_occ, n_vtx); 
    asg64_v b64, ub64; kv_init(b64); kv_init(ub64); 
    
    for (k = 0; k < int_idx_n; k++) {///scan all integer contigs
        r_a = int_a + (int_idx[k]>>32); r_n = (uint32_t)int_idx[k];
        assert(r_n >= 2);
        for (i = 1; i + 1 < r_n; i++) {
            ng_occ[r_a[i]]++; ng_occ[r_a[i]^1]++;
        }
        ng_occ[r_a[0]]++; ng_occ[r_a[r_n-1]^1]++;
    }

    ma_ug_t *un_g = ma_ug_hybrid_gen(ng); ma_utg_t *u;
    // print_debug_gfa(uidx->sg, un_g, uidx->uopt->coverage_cut, "iig0", uidx->uopt->sources, 
    // uidx->uopt->ruIndex, uidx->uopt->max_hang, uidx->uopt->min_ovlp, 0, 0, 0);
    int32_t ui, un;
    for (k = 0; k < un_g->u.n; k++) {///all unitigs of raw utg 
        u = &(un_g->u.a[k]); 
        zs = ze = (uint32_t)-1; un = u->n;
        for (ui = 0; ui < un; ui++) {
            if(occ_m(ng_occ[(u->a[ui]>>32)^1]) == 1) {
                zs = ui; break;
            }
        }

        for (ui = ((int32_t)un)-1; ui >= 0; ui--) {
            if(occ_m(ng_occ[u->a[ui]>>32]) == 1) {
                ze = ui; break;
            }
        }
        if(zs != (uint32_t)-1 && ze != (uint32_t)-1 && zs > ze) continue;
        if(zs != (uint32_t)-1) ng_occ[(u->a[zs]>>32)^1] |= ((uint32_t)0x80000000);
        if(ze != (uint32_t)-1) ng_occ[(u->a[ze]>>32)] |= ((uint32_t)0x80000000);
    }
    ma_ug_destroy(un_g);
    
    // fprintf(stderr, ">>>>>>[M::%s::] int_idx[0]::%lu, int_idx[1]::%lu\n", __func__, int_idx[0], int_idx[1]);
    // prt_intg_info(int_idx, int_idx_n, int_a, "intg");
    for (i = b64.n = ua_n = a_n = 0; i < int_idx_n; i++) {///scan all integer contigs
        s = int_idx[i]>>32; e = s + ((uint32_t)int_idx[i]); 
        assert(e > s + 1);//the length is at least 2
        for (k = s, bn = b64.n, pi = (uint64_t)-1; k < e; k++) {
            v = int_a[k];
            if((!(ng_occ[v]&((uint32_t)0x80000000)))&&(!(ng_occ[v^1]&((uint32_t)0x80000000)))) {
                continue;///it must be a unique node
            }
            if((pi != (uint64_t)-1) && (ng_occ[v^1]&((uint32_t)0x80000000))) {
                // fprintf(stderr, "+++[M::%s::i->%u::s->%u::e->%u] pi::%lu, k::%u\n", __func__, i, s, e, pi, k);
                kv_pushp(uint64_t, b64, &pz); *pz = pi; (*pz) <<= 32; (*pz) |= k; a_n++;
            }
            if(ng_occ[v]&((uint32_t)0x80000000)) pi = k;
        }

        b64_n = b64.n;
        for (k = bn, pi = s; k < b64_n; k++) {
            ei = b64.a[k]>>32;
            if(pi < ei) {
                kv_pushp(uint64_t, b64, &pz); ua_n++;
                (*pz) = pi; (*pz) <<= 32; (*pz) |= ei; (*pz) |= ((uint64_t)0x8000000000000000);
            }
            pi = (uint32_t)b64.a[k];
        }

        ei = e - 1;
        if(pi < ei) {
            kv_pushp(uint64_t, b64, &pz); ua_n++;
            (*pz) = pi; (*pz) <<= 32; (*pz) |= ei; (*pz) |= ((uint64_t)0x8000000000000000);
        }
    }
    assert(a_n + ua_n == b64.n);
    // prt_thread_info(b64.a, a_n+ua_n, NULL, 0, "tt_minus");
    // fprintf(stderr, "**0**[M::%s::] a_n::%lu, ua_n::%lu\n", __func__, a_n, ua_n);

    uint64_t *i_idx, n_clus; CALLOC(i_idx, ng->n<<1);
    radix_sort_srt64(b64.a, b64.a + b64.n);///keeps the coordinates within int_a[]

    // prt_thread_info(b64.a, a_n, NULL, 0, "tt0");
    /*********debugging*********/
    // for (i = 0; i < a_n; i++) { ///available intervals
    //     s = b64.a[i]>>32; e = (uint32_t)b64.a[i];
    //     fprintf(stderr, "[M::%s::ava->%u] s::%u, e::%u\n", __func__, i, s, e);
    //     for (k = s; k <= e; k++) {
    //         fprintf(stderr, "utg%.6dl(%c)\n", ((int32_t)(int_a[k]>>1))+1, "+-"[int_a[k]&1]);
    //     }
    // }
    // for (i = a_n; i < b64.n; i++) {///unavailable intervals
    //     s = (b64.a[i]<<1)>>33; e = (uint32_t)b64.a[i];
    //     fprintf(stderr, "[M::%s::uava->%u] s::%u, e::%u\n", __func__, i - (uint32_t)a_n, s, e);
    //     for (k = s; k <= e; k++) {
    //         fprintf(stderr, "utg%.6dl(%c)\n", ((int32_t)(int_a[k]>>1))+1, "+-"[int_a[k]&1]);
    //     }
    // }
    /*********debugging*********/






    for (i = a_n; i < b64.n; i++) {///unavailable intervals; mask all unavailable nodes
        assert(b64.a[i]&((uint64_t)0x8000000000000000)); ua_n--;
        s = (b64.a[i]<<1)>>33; e = (uint32_t)b64.a[i]; assert(e > s);//[s, e]
        for (k = s + 1; k < e; k++) {///note: here is [s, e]
            i_idx[int_a[k]] |= ((uint64_t)0x8000000000000000);
            i_idx[int_a[k]^1] |= ((uint64_t)0x8000000000000000);
        }
        i_idx[int_a[s]] |= ((uint64_t)0x8000000000000000);
        i_idx[int_a[e]^1] |= ((uint64_t)0x8000000000000000);
    }
    // prt_thread_info(b64.a, a_n, NULL, 0, "tt1");
    ///unavailable intervals are useless
    b64.n = a_n; assert(ua_n == 0);
    for (i = 0; i < a_n; i++) { ///available intervals
        s = b64.a[i]>>32; e = (uint32_t)b64.a[i]; assert(e > s);///[s, e] -> coordinates within int_a[]
        if(i == 3201 || i == 3202) {
            fprintf(stderr, "[M::%s::] s::%u, e::%u, int_a[s]::%lu(occ::%u), int_a[e]^1::%lu(occ::%u)\n", 
            __func__, s, e, int_a[s], occ_m(ng_occ[int_a[s]]),
            int_a[e]^1, occ_m(ng_occ[int_a[e]^1]));
        }
        
        assert(!(b64.a[i]&((uint64_t)0x8000000000000000)));
        for (k = s + 1; k < e; k++) {///note: here is [s, e]; s && e are unique, but [s+1, e-1] are not unique
            kv_pushp(uint64_t, b64, &pz); //i_idx[int_a[k]]++;
            (*pz) = int_a[k]; (*pz) <<= 32; (*pz) |= i;///(raw unitig/non-unqiue node id)|(integer contig id)

            kv_pushp(uint64_t, b64, &pz); //i_idx[int_a[k]^1]++;
            (*pz) = int_a[k]^1; (*pz) <<= 32; (*pz) |= i;
        }
        kv_pushp(uint64_t, b64, &pz); //i_idx[int_a[s]]++;
        (*pz) = int_a[s]; (*pz) <<= 32; (*pz) |= i;

        kv_pushp(uint64_t, b64, &pz); //i_idx[int_a[e]^1]++;
        (*pz) = int_a[e]^1; (*pz) <<= 32; (*pz) |= i;
    }
    // prt_thread_info(b64.a, a_n, NULL, 0, "tt2");
    // fprintf(stderr, "**1**[M::%s::] a_n::%lu, ua_n::%lu\n", __func__, a_n, ua_n);
    ///index
    radix_sort_srt64(b64.a + a_n, b64.a + b64.n);///(raw unitig node id)|(integer contig id)
    for (k = a_n + 1, i = a_n; k <= b64.n; k++) {
        if(k == b64.n || (b64.a[k]>>32) != (b64.a[i]>>32)) {
            i_idx[b64.a[i]>>32] |= (((uint64_t)i)<<32)|((uint64_t)k);///b64.a[i]>>32 appear once (unique ends)/multipe times
            i = k;
        }
    }
    // prt_thread_info(b64.a, a_n, NULL, 0, "tt3");
    // fprintf(stderr, "**2**[M::%s::] a_n::%lu, ua_n::%lu\n", __func__, a_n, ua_n);
    ///b64.a[0, a_n]:: all resolvable paths with unique beg && end
    ///b64.a[a_n, b64.n]:: (raw unitig/non-unqiue node id)|(resolvable path id)
    // n_clus = usg_unique_arcs_cluster(&b64, a_n, i_idx, int_a);///this function might be wrong
    n_clus = usg_unique_arcs_cluster_adv(&b64, a_n, i_idx, int_a, &ub64);
    fprintf(stderr, "**3**[M::%s::] a_n::%lu, n_clus::%lu\n", __func__, a_n, n_clus);
    // prt_thread_info(b64.a, a_n, NULL, 0, "tt4");
    assert(b64.n == (a_n<<1));
    old_path_ext(uidx, ng, &b64, &ub64, a_n, n_clus, i_idx, int_a, ff, ng_occ, max_ext);
    // new_path_ext(uidx, ng, &b64, &ub64, a_n, n_clus, i_idx, int_a, ff, ng_occ, max_ext);

    free(ng_occ); free(i_idx); free(ff); kv_destroy(b64); kv_destroy(ub64); 
    return b64.n - a_n;
}

uint64_t gen_unique_g_adv(ul_resolve_t *uidx, usg_t *ng, uint64_t *int_idx, uint64_t int_idx_n, uint64_t *int_a, uint32_t max_ext)
{
    uint64_t pi, ei, *pz, bn, a_n, ua_n, *r_a, r_n; uint8_t *ff; CALLOC(ff, ng->n<<1);
    uint32_t i, k, n_vtx = ng->n<<1, s, e, zs, ze, v, b64_n; 
    uint32_t *ng_occ; CALLOC(ng_occ, n_vtx); 
    asg64_v b64, ub64; kv_init(b64); kv_init(ub64); 
    uint64_t *i_idx, n_clus; CALLOC(i_idx, ng->n<<1);
    
    for (k = 0; k < int_idx_n; k++) {///scan all integer contigs
        r_a = int_a + (int_idx[k]>>32); r_n = (uint32_t)int_idx[k];
        assert(r_n >= 2);
        for (i = 1; i + 1 < r_n; i++) {
            ng_occ[r_a[i]]++; ng_occ[r_a[i]^1]++;
        }
        ng_occ[r_a[0]]++; ng_occ[r_a[r_n-1]^1]++;
    }

    ma_ug_t *un_g = ma_ug_hybrid_gen(ng); ma_utg_t *u;
    // print_debug_gfa(uidx->sg, un_g, uidx->uopt->coverage_cut, "iig0", uidx->uopt->sources, 
    // uidx->uopt->ruIndex, uidx->uopt->max_hang, uidx->uopt->min_ovlp, 0, 0, 0);
    int32_t ui, un;
    for (k = 0; k < un_g->u.n; k++) {///all unitigs of raw utg 
        u = &(un_g->u.a[k]); 
        zs = ze = (uint32_t)-1; un = u->n;
        for (ui = 0; ui < un; ui++) {
            if(occ_m(ng_occ[(u->a[ui]>>32)^1]) == 1) {
                zs = ui; break;
            }
        }

        for (ui = ((int32_t)un)-1; ui >= 0; ui--) {
            if(occ_m(ng_occ[u->a[ui]>>32]) == 1) {
                ze = ui; break;
            }
        }
        if(zs != (uint32_t)-1 && ze != (uint32_t)-1 && zs > ze) continue;
        if(zs != (uint32_t)-1) ng_occ[(u->a[zs]>>32)^1] |= ((uint32_t)0x80000000);
        if(ze != (uint32_t)-1) ng_occ[(u->a[ze]>>32)] |= ((uint32_t)0x80000000);
    }
    ma_ug_destroy(un_g);
    
    // fprintf(stderr, ">>>>>>[M::%s::] int_idx[0]::%lu, int_idx[1]::%lu\n", __func__, int_idx[0], int_idx[1]);
    // prt_intg_info(int_idx, int_idx_n, int_a, "intg");
    for (i = b64.n = ua_n = a_n = 0; i < int_idx_n; i++) {///scan all integer contigs
        s = int_idx[i]>>32; e = s + ((uint32_t)int_idx[i]); 
        assert(e > s + 1);//the length is at least 2
        for (k = s, bn = b64.n, pi = (uint64_t)-1; k < e; k++) {
            v = int_a[k]; 
            if((!(ng_occ[v]&((uint32_t)0x80000000)))&&(!(ng_occ[v^1]&((uint32_t)0x80000000)))) {
                continue;///it must be a unique node
            }
            if((pi != (uint64_t)-1) && (ng_occ[v^1]&((uint32_t)0x80000000))) {
                // fprintf(stderr, "+++[M::%s::i->%u::s->%u::e->%u] pi::%lu, k::%u\n", __func__, i, s, e, pi, k);
                pz = (b64.n>0)? &(b64.a[b64.n-1]):(NULL);
                if((!pz) || (((*pz)>>32) != pi)) {//keep the shortest pi<->k 
                    kv_pushp(uint64_t, b64, &pz); 
                    *pz = pi; (*pz) <<= 32; (*pz) |= k;
                    a_n++;
                }
            }
            if(ng_occ[v]&((uint32_t)0x80000000)) pi = k;//keep the shortest pi<->k
        }

        b64_n = b64.n;
        for (k = bn, pi = s; k < b64_n; k++) {
            ei = b64.a[k]>>32;
            if(pi < ei) {
                kv_pushp(uint64_t, b64, &pz); ua_n++;
                (*pz) = pi; (*pz) <<= 32; (*pz) |= ei; (*pz) |= ((uint64_t)0x8000000000000000);
            }
            pi = (uint32_t)b64.a[k];
        }

        ei = e - 1;
        if(pi < ei) {
            kv_pushp(uint64_t, b64, &pz); ua_n++;
            (*pz) = pi; (*pz) <<= 32; (*pz) |= ei; (*pz) |= ((uint64_t)0x8000000000000000);
        }
    }
    assert(a_n + ua_n == b64.n);
    // prt_thread_info(b64.a, a_n+ua_n, NULL, 0, "tt_minus");
    // fprintf(stderr, "**0**[M::%s::] a_n::%lu, ua_n::%lu\n", __func__, a_n, ua_n);
    radix_sort_srt64(b64.a, b64.a + b64.n);///keeps the coordinates within int_a[]

    for (i = a_n; i < b64.n; i++) {///unavailable intervals; mask all unavailable nodes
        assert(b64.a[i]&((uint64_t)0x8000000000000000)); ua_n--;
        s = (b64.a[i]<<1)>>33; e = (uint32_t)b64.a[i]; assert(e > s);//[s, e]
        for (k = s + 1; k < e; k++) {///note: here is [s, e]
            i_idx[int_a[k]] |= ((uint64_t)0x8000000000000000);
            i_idx[int_a[k]^1] |= ((uint64_t)0x8000000000000000);
        }
        i_idx[int_a[s]] |= ((uint64_t)0x8000000000000000);
        i_idx[int_a[e]^1] |= ((uint64_t)0x8000000000000000);
    }
    // prt_thread_info(b64.a, a_n, NULL, 0, "tt1");
    ///unavailable intervals are useless
    b64.n = a_n; assert(ua_n == 0);
    for (i = 0; i < a_n; i++) { ///available intervals
        s = b64.a[i]>>32; e = (uint32_t)b64.a[i]; assert(e > s);///[s, e] -> coordinates within int_a[]
        assert(!(b64.a[i]&((uint64_t)0x8000000000000000)));

        for (k = s + 1; k < e; k++) {///note: here is [s, e]; s && e are unique, but [s+1, e-1] are not unique
            kv_pushp(uint64_t, b64, &pz); //i_idx[int_a[k]]++;
            (*pz) = int_a[k]; (*pz) <<= 32; (*pz) |= i;///(raw unitig/non-unqiue node id)|(integer contig id)

            kv_pushp(uint64_t, b64, &pz); //i_idx[int_a[k]^1]++;
            (*pz) = int_a[k]^1; (*pz) <<= 32; (*pz) |= i;
        }
        kv_pushp(uint64_t, b64, &pz); //i_idx[int_a[s]]++;
        (*pz) = int_a[s]; (*pz) <<= 32; (*pz) |= i;

        kv_pushp(uint64_t, b64, &pz); //i_idx[int_a[e]^1]++;
        (*pz) = int_a[e]^1; (*pz) <<= 32; (*pz) |= i;
    }
    // prt_thread_info(b64.a, a_n, NULL, 0, "tt2");
    // fprintf(stderr, "**1**[M::%s::] a_n::%lu, ua_n::%lu\n", __func__, a_n, ua_n);
    ///index
    radix_sort_srt64(b64.a + a_n, b64.a + b64.n);///(raw unitig node id)|(integer contig id)
    for (k = a_n + 1, i = a_n; k <= b64.n; k++) {
        if(k == b64.n || (b64.a[k]>>32) != (b64.a[i]>>32)) {
            i_idx[b64.a[i]>>32] |= (((uint64_t)i)<<32)|((uint64_t)k);///b64.a[i]>>32 appear once (unique ends)/multipe times
            i = k;
        }
    }
    // prt_thread_info(b64.a, a_n, NULL, 0, "tt3");
    // fprintf(stderr, "**2**[M::%s::] a_n::%lu, ua_n::%lu\n", __func__, a_n, ua_n);
    ///b64.a[0, a_n]:: all resolvable paths with unique beg && end
    ///b64.a[a_n, b64.n]:: (raw unitig/non-unqiue node id)|(resolvable path id)
    // n_clus = usg_unique_arcs_cluster(&b64, a_n, i_idx, int_a);///this function might be wrong
    n_clus = usg_unique_arcs_cluster_adv(&b64, a_n, i_idx, int_a, &ub64);
    // fprintf(stderr, "**3**[M::%s::] a_n::%lu, n_clus::%lu\n", __func__, a_n, n_clus);
    // prt_thread_info(b64.a, a_n, NULL, 0, "tt4");
    assert(b64.n == (a_n<<1));
    old_path_ext(uidx, ng, &b64, &ub64, a_n, n_clus, i_idx, int_a, ff, ng_occ, max_ext);
    // new_path_ext(uidx, ng, &b64, &ub64, a_n, n_clus, i_idx, int_a, ff, ng_occ, max_ext);//wrong

    free(ng_occ); free(i_idx); free(ff); kv_destroy(b64); kv_destroy(ub64); 
    return b64.n - a_n;
}

void u2g_hybrid_detan(ul_resolve_t *uidx, usg_t *ng, uint32_t max_ext, asg64_v *in, asg64_v *ib)
{
    uint64_t k, i, x, m, *tmp, sn; asg64_v tx = {0,0,0}, tb = {0,0,0}, *ob = NULL, *ub = NULL; 
    ob = (in?(in):(&tx)); ub = (ib?(ib):(&tb)); ob->n = ub->n = 0;

    for (k = 0; k < uidx->str_b.n_thread; k++) {
        uidx->str_b.buf[k].res_dump.n = uidx->str_b.buf[k].u.n = uidx->str_b.buf[k].o.n = 0;
    }
    kt_for(uidx->str_b.n_thread, worker_integer_realign_g, uidx, uidx->uovl.i_ug->u.n);

    for (k = ob->n = ub->n = m = 0; k < uidx->str_b.n_thread; k++) {
        for (i = 0; i < uidx->str_b.buf[k].res_dump.n; i++) {
            x = uidx->str_b.buf[k].res_dump.a[i];
            kv_push(uint64_t, *ob, x);//aln details
            if(x&((uint64_t)0x8000000000000000)) {
                x -= ((uint64_t)0x8000000000000000); x >>= 32; x <<= 32;//seq_id
                x |= ob->n;//offset
                kv_push(uint64_t, *ub, x);///idx:: seq_id|offset_in_ob
            } else {
                m++;
            }
        }
    }

    radix_sort_srt64(ub->a, ub->a + ub->n);
    kv_resize(uint64_t, *ub, ub->n+m); tmp = ub->a + ub->n; m = 0;
    for (k = 0; k < ub->n; k++) {
        sn = ((uint32_t)(ob->a[((uint32_t)ub->a[k])-1]));
        memcpy(tmp + m, ob->a + ((uint32_t)ub->a[k]), sn*sizeof((*tmp))); 
        ub->a[k] = m; ub->a[k] <<= 32; ub->a[k] |= sn;//offset_in_ob|occ
        m += sn;
    }
    assert(m <= ob->n);
    memcpy(ob->a, tmp, m*sizeof((*tmp))); ob->n = m;


    // for (k = 0, p = NULL; k < ub->n; k++) {
    //     p = &(ob->a[((uint32_t)ub->a[k])-1]);
    //     x = ub->a[k]<<32;//offset
    //     x |= ((uint32_t)(*p));///occ
    //     ub->a[k] = x;
    //     (*p) >>= 32; (*p) <<= 32; (*p) |= k; 
    // }
    // for (k = m = 0; k < ob->n; k++) {
    //     if(ob->a[k]&((uint64_t)0x8000000000000000)) {
    //         x = m; x <<= 32; x |= ((uint32_t)ub->a[(uint32_t)ob->a[k]]);
    //         ub->a[(uint32_t)ob->a[k]] = x;
    //     } else {
    //         ob->a[m++] = ob->a[k];
    //     }
    // }
    // ob->n = m;

    // for (k = ob->n = ub->n = 0; k < uidx->str_b.n_thread; k++) {
    //     for (i = 0; i < uidx->str_b.buf[k].res_dump.n; i++) {
    //         if((uidx->str_b.buf[k].res_dump.a[i]>>32)==((uint32_t)-1)) {
    //             x = ob->n; x <<= 32; x |= ((uint32_t)uidx->str_b.buf[k].res_dump.a[i]);
    //             kv_push(uint64_t, *ub, x);///idx:: offset_in_ob|occ
    //         } else {
    //             kv_push(uint64_t, *ob, uidx->str_b.buf[k].res_dump.a[i]);//aln details
    //         }
    //     }        
    // }

    /*********debugging*********/
    fprintf(stderr, "[M::%s::] # iug::%u, # gchain::%u\n", __func__, (uint32_t)uidx->uovl.i_ug->u.n, (uint32_t)ub->n);
    // for (k = 0; k < ub->n; k++) {
    //     // fprintf(stderr, "[M::%s::k->%lu] # iug::%u, # chain::%u, off chain::%u\n", 
    //     // __func__, k, (uint32_t)uidx->uovl.cc.iug_a[k].n, (uint32_t)ub->a[k], (uint32_t)(ub->a[k]>>32));
    //     for (i = 0; i < (uint32_t)ub->a[k]; i++) {
    //         // fprintf(stderr, "utg%.6dl(%c) <------> utg%.6dl(%c)\n", 
    //         // ((int32_t)(uidx->uovl.cc.iug_a[k].a[i].v>>1))+1, "+-"[uidx->uovl.cc.iug_a[k].a[i].v&1],
    //         // ((int32_t)(ob->a[(ub->a[k]>>32)+i]>>1))+1, "+-"[ob->a[(ub->a[k]>>32)+i]&1]);
    //         assert(uidx->uovl.cc.iug_a[k].a[i].v == ob->a[(ub->a[k]>>32)+i]);
    //     }
    // }
    /*********debugging*********/

    ///ub->idx; ob->nodes
    // u_ug = gen_unique_g(uidx, ng, ng_occ, ub->a, ub->n, ob->a);//ma_ug_hybrid_gen(ng);
    ///debug
    debug_sysm_usg_t(ng, __func__);
    // prt_usg_t(uidx, ng, "ng4");
    if(gen_unique_g_adv(uidx, ng, ub->a, ub->n, ob->a, max_ext)) {
        // usg_arc_t *z = get_usg_arc(ng, 2, 576), *q = get_usg_arc(ng, 577, 3);
        // fprintf(stderr, "xxxx0xxx[M::%s::] p->del::%u, q->del::%u\n", 
        //             __func__, z?z->del:1, q?q->del:1);
        ///debug
        debug_sysm_usg_t(ng, __func__);
        // z = get_usg_arc(ng, 2, 576); q = get_usg_arc(ng, 577, 3);
        // fprintf(stderr, "xxxx1xxx[M::%s::] p->del::%u, q->del::%u\n", 
        //             __func__, z?z->del:1, q?q->del:1);
        // fprintf(stderr, "+[M::%s::] ng->n::%u\n", __func__, (uint32_t)ng->n);
        usg_cleanup(ng);
        ///debug
        debug_sysm_usg_t(ng, __func__);
        // fprintf(stderr, "-[M::%s::] ng->n::%u\n", __func__, (uint32_t)ng->n);
    }
    // prt_usg_t(uidx, ng, "ng_dbg");
    if(!in) free(tx.a); if(!ib) free(tb.a);
}

void u2g_hybrid_aln(ul_resolve_t *uidx, usg_t *ng, asg64_v *ob, asg64_v *ub)
{
    uint64_t k, m, i, x, *tmp, sn;
    ob->n = ub->n = 0;
    for (k = 0; k < uidx->str_b.n_thread; k++) {
        uidx->str_b.buf[k].res_dump.n = uidx->str_b.buf[k].u.n = uidx->str_b.buf[k].o.n = 0;
    }
    kt_for(uidx->str_b.n_thread, worker_integer_realign_g, uidx, uidx->uovl.i_ug->u.n);

    for (k = ob->n = ub->n = m = 0; k < uidx->str_b.n_thread; k++) {
        for (i = 0; i < uidx->str_b.buf[k].res_dump.n; i++) {
            x = uidx->str_b.buf[k].res_dump.a[i];
            kv_push(uint64_t, *ob, x);//aln details
            if(x&((uint64_t)0x8000000000000000)) {
                x -= ((uint64_t)0x8000000000000000); x >>= 32; x <<= 32;//seq_id
                x |= ob->n;//offset
                kv_push(uint64_t, *ub, x);///idx:: seq_id|offset_in_ob
            } else {
                m++;
            }
        }
    }

    radix_sort_srt64(ub->a, ub->a + ub->n);
    kv_resize(uint64_t, *ub, ub->n+m); tmp = ub->a + ub->n; m = 0;
    for (k = 0; k < ub->n; k++) {
        sn = ((uint32_t)(ob->a[((uint32_t)ub->a[k])-1]));
        memcpy(tmp + m, ob->a + ((uint32_t)ub->a[k]), sn*sizeof((*tmp))); 
        ub->a[k] = m; ub->a[k] <<= 32; ub->a[k] |= sn;//offset_in_ob|occ
        m += sn;
    }
    assert(m <= ob->n);
    memcpy(ob->a, tmp, m*sizeof((*tmp))); ob->n = m;
}

uint32_t ug_ext(ul_resolve_t *uidx, usg_t *ng, uint64_t *int_idx, uint64_t int_idx_n, uint64_t *int_a, uint32_t max_ext, 
uint8_t *ff, uint32_t *ng_occ, uint64_t *i_idx, asg64_v *b64, asg64_v *ub64)
{
    uint32_t n_vtx = ng->n<<1, k, i; uint64_t pi, ei, *pz, v, b64_n;
    uint64_t *r_a, r_n, zs, ze, s, e, a_n, ua_n, bn, n_clus;
    memset(ff, 0, sizeof((*ff))*n_vtx);
    memset(ng_occ, 0, sizeof((*ng_occ))*n_vtx); 
    memset(i_idx, 0, sizeof((*i_idx))*n_vtx); 
    b64->n = ub64->n = 0;

    for (k = 0; k < int_idx_n; k++) {///scan all integer contigs
        r_a = int_a + (int_idx[k]>>32); r_n = (uint32_t)int_idx[k];
        assert(r_n >= 2);
        for (i = 1; i + 1 < r_n; i++) {
            ng_occ[r_a[i]]++; ng_occ[r_a[i]^1]++;
        }
        ng_occ[r_a[0]]++; ng_occ[r_a[r_n-1]^1]++;
    }

    ma_ug_t *un_g = ma_ug_hybrid_gen(ng); 
    int32_t ui, un; ma_utg_t *u = NULL;
    for (k = 0; k < un_g->u.n; k++) {///all unitigs of raw utg 
        u = &(un_g->u.a[k]); 
        zs = ze = (uint32_t)-1; un = u->n;
        for (ui = 0; ui < un; ui++) {
            if(occ_m(ng_occ[(u->a[ui]>>32)^1]) == 1) {
                zs = ui; break;
            }
        }

        for (ui = ((int32_t)un)-1; ui >= 0; ui--) {
            if(occ_m(ng_occ[u->a[ui]>>32]) == 1) {
                ze = ui; break;
            }
        }
        if(zs != (uint32_t)-1 && ze != (uint32_t)-1 && zs > ze) continue;
        if(zs != (uint32_t)-1) ng_occ[(u->a[zs]>>32)^1] |= ((uint32_t)0x80000000);
        if(ze != (uint32_t)-1) ng_occ[(u->a[ze]>>32)] |= ((uint32_t)0x80000000);
    }
    ma_ug_destroy(un_g);

        // fprintf(stderr, ">>>>>>[M::%s::] int_idx[0]::%lu, int_idx[1]::%lu\n", __func__, int_idx[0], int_idx[1]);
    // prt_intg_info(int_idx, int_idx_n, int_a, "intg");
    for (i = b64->n = ua_n = a_n = 0; i < int_idx_n; i++) {///scan all integer contigs
        s = int_idx[i]>>32; e = s + ((uint32_t)int_idx[i]); 
        assert(e > s + 1);//the length is at least 2
        for (k = s, bn = b64->n, pi = (uint64_t)-1; k < e; k++) {
            v = int_a[k]; 
            if((!(ng_occ[v]&((uint32_t)0x80000000)))&&(!(ng_occ[v^1]&((uint32_t)0x80000000)))) {
                continue;///it must be a unique node
            }
            if((pi != (uint64_t)-1) && (ng_occ[v^1]&((uint32_t)0x80000000))) {
                pz = (b64->n>0)? &(b64->a[b64->n-1]):(NULL);
                if((!pz) || (((*pz)>>32) != pi)) {//keep the shortest pi<->k 
                    kv_pushp(uint64_t, *b64, &pz); 
                    *pz = pi; (*pz) <<= 32; (*pz) |= k;
                    a_n++;
                }
            }
            if(ng_occ[v]&((uint32_t)0x80000000)) pi = k;//keep the shortest pi<->k
        }

        b64_n = b64->n;
        for (k = bn, pi = s; k < b64_n; k++) {
            ei = b64->a[k]>>32;
            if(pi < ei) {
                kv_pushp(uint64_t, *b64, &pz); ua_n++;
                (*pz) = pi; (*pz) <<= 32; (*pz) |= ei; (*pz) |= ((uint64_t)0x8000000000000000);
            }
            pi = (uint32_t)b64->a[k];
        }

        ei = e - 1;
        if(pi < ei) {
            kv_pushp(uint64_t, *b64, &pz); ua_n++;
            (*pz) = pi; (*pz) <<= 32; (*pz) |= ei; (*pz) |= ((uint64_t)0x8000000000000000);
        }
    }
    assert(a_n + ua_n == b64->n);
    // prt_thread_info(b64.a, a_n+ua_n, NULL, 0, "tt_minus");
    // fprintf(stderr, "**0**[M::%s::] a_n::%lu, ua_n::%lu\n", __func__, a_n, ua_n);
    radix_sort_srt64(b64->a, b64->a + b64->n);///keeps the coordinates within int_a[]


    for (i = a_n; i < b64->n; i++) {///unavailable intervals; mask all unavailable nodes
        assert(b64->a[i]&((uint64_t)0x8000000000000000)); ua_n--;
        s = (b64->a[i]<<1)>>33; e = (uint32_t)b64->a[i]; assert(e > s);//[s, e]
        for (k = s + 1; k < e; k++) {///note: here is [s, e]
            i_idx[int_a[k]] |= ((uint64_t)0x8000000000000000);
            i_idx[int_a[k]^1] |= ((uint64_t)0x8000000000000000);
        }
        i_idx[int_a[s]] |= ((uint64_t)0x8000000000000000);
        i_idx[int_a[e]^1] |= ((uint64_t)0x8000000000000000);
    }
    // prt_thread_info(b64.a, a_n, NULL, 0, "tt1");
    ///unavailable intervals are useless
    b64->n = a_n; assert(ua_n == 0);
    for (i = 0; i < a_n; i++) { ///available intervals
        s = b64->a[i]>>32; e = (uint32_t)b64->a[i]; assert(e > s);///[s, e] -> coordinates within int_a[]
        assert(!(b64->a[i]&((uint64_t)0x8000000000000000)));

        for (k = s + 1; k < e; k++) {///note: here is [s, e]; s && e are unique, but [s+1, e-1] are not unique
            kv_pushp(uint64_t, *b64, &pz); //i_idx[int_a[k]]++;
            (*pz) = int_a[k]; (*pz) <<= 32; (*pz) |= i;///(raw unitig/non-unqiue node id)|(integer contig id)

            kv_pushp(uint64_t, *b64, &pz); //i_idx[int_a[k]^1]++;
            (*pz) = int_a[k]^1; (*pz) <<= 32; (*pz) |= i;
        }
        kv_pushp(uint64_t, *b64, &pz); //i_idx[int_a[s]]++;
        (*pz) = int_a[s]; (*pz) <<= 32; (*pz) |= i;

        kv_pushp(uint64_t, *b64, &pz); //i_idx[int_a[e]^1]++;
        (*pz) = int_a[e]^1; (*pz) <<= 32; (*pz) |= i;
    }
    // prt_thread_info(b64.a, a_n, NULL, 0, "tt2");
    // fprintf(stderr, "**1**[M::%s::] a_n::%lu, ua_n::%lu\n", __func__, a_n, ua_n);
    ///index
    radix_sort_srt64(b64->a + a_n, b64->a + b64->n);///(raw unitig node id)|(integer contig id)
    for (k = a_n + 1, i = a_n; k <= b64->n; k++) {
        if(k == b64->n || (b64->a[k]>>32) != (b64->a[i]>>32)) {
            i_idx[b64->a[i]>>32] |= (((uint64_t)i)<<32)|((uint64_t)k);///b64.a[i]>>32 appear once (unique ends)/multipe times
            i = k;
        }
    }
    // prt_thread_info(b64.a, a_n, NULL, 0, "tt3");
    // fprintf(stderr, "**2**[M::%s::] a_n::%lu, ua_n::%lu\n", __func__, a_n, ua_n);
    ///b64.a[0, a_n]:: all resolvable paths with unique beg && end
    ///b64.a[a_n, b64.n]:: (raw unitig/non-unqiue node id)|(resolvable path id)
    // n_clus = usg_unique_arcs_cluster(&b64, a_n, i_idx, int_a);///this function might be wrong
    n_clus = usg_unique_arcs_cluster_adv(b64, a_n, i_idx, int_a, ub64);
    // fprintf(stderr, "**3**[M::%s::] a_n::%lu, n_clus::%lu\n", __func__, a_n, n_clus);
    // prt_thread_info(b64.a, a_n, NULL, 0, "tt4");
    assert(b64->n == (a_n<<1));
    old_path_ext(uidx, ng, b64, ub64, a_n, n_clus, i_idx, int_a, ff, ng_occ, max_ext);
    // new_path_ext(uidx, ng, &b64, &ub64, a_n, n_clus, i_idx, int_a, ff, ng_occ, max_ext);//wrong
    return b64->n - a_n;
}

void merge_hybrid_utg_content(ma_utg_t* cc, ma_ug_t* raw, asg_t* rg, usg_t *ng, kvec_asg_arc_t_warp* edge);
void renew_usg_t_bub(ul_resolve_t *uidx, usg_t *ng, uint32_t *id_map, uint8_t *ff, uint32_t rocc_cut)
{
    ma_ug_t *ug = ma_ug_hybrid_gen(ng); 
    uint32_t i, k, v, p[2]; ma_utg_t *u; p[0] = 1; p[1] = 2;
    memset(id_map, -1, sizeof((*id_map))*ng->n);
    memset(ff, 0, sizeof((*ff))*(ng->n<<1));
    kvec_asg_arc_t_warp e; kv_init(e.a); e.i = 0;
    // for (i = dn = 0; i < ng->n; i++) {
    //     if(ng->a[i].del) continue;
    //     dn++;
    // }
    // fprintf(stderr, "ng->n::%u, ug->u.n::%u, dn::%u\n", (uint32_t)ng->n, (uint32_t)ug->u.n, dn);
    // dn = 0
    for (i = 0; i < ug->u.n; i++) {
        ug->g->seq[i].c = PRIMARY_LABLE;
        u = &(ug->u.a[i]);
        if(u->m == 0) continue;
        for (k = 0; k < u->n; k++) {
            v = u->a[k]>>32; 
            id_map[v>>1] = i<<2; 
            if(k == 0) id_map[v>>1] |= p[(v^1)&1];
            if(k + 1 == u->n) id_map[v>>1] |= p[(v)&1];
            // dn++;
        }
        // fprintf(stderr, "+[M::%s] i::%u\n", __func__, i);
        merge_hybrid_utg_content(u, uidx->l1_ug, uidx->sg, ng, &e);
        // fprintf(stderr, "-[M::%s] i::%u\n", __func__, i);
        ug->g->seq[i].len = u->len;
    }
    kv_destroy(e.a);

    destory_bubbles(uidx->bub); free(uidx->bub); CALLOC(uidx->bub, 1);
    // fprintf(stderr, "[M::%s] homozygous read coverage threshold: %d\n", __func__, asm_opt.hom_global_coverage_set?
    //                         asm_opt.hom_global_coverage:(int)(((double)asm_opt.hom_global_coverage)/((double)HOM_PEAK_RATE)));
    // identify_bubbles(ug, uidx->bub, uidx->r_het, NULL);
    if(asm_opt.polyploidy <= 2) {
        identify_bubbles_recal(uidx->sg, ug, uidx->bub, uidx->r_het, uidx->uopt->sources, uidx->uopt->ruIndex, NULL);
    } else {
        identify_bubbles_recal_poy(uidx->sg, ug, uidx->bub, uidx->r_het, uidx->uopt->sources, uidx->uopt->ruIndex, NULL);
    }
    
    // fprintf(stderr, "0[M::%s::] f[51]::%u\n", __func__, ff[51]);
    for (i = 0; i < ng->n; i++) {
        if(id_map[i] != (uint32_t)-1) {
            k = id_map[i]>>2; 
            // fprintf(stderr, "k::%u, ng->n::%u, ug->u.n::%u, dn::%u\n", k, (uint32_t)ng->n, (uint32_t)ug->u.n, dn);
            // if(i == 7075 || i == 28174 || i == 77111 || i == 3826 || i == 12150 || i == 58312 || i == 59190 || i == 72134) {
            //     fprintf(stderr, "[M::%s::k->%u] utg%.6dl(%c), ug->u.a[k].n::%u, rocc_cut::%u, is_hom::%u\n", __func__, k, 
            //         i+1, "+-"[0], (uint32_t)ug->u.a[k].n, rocc_cut, IF_HOM(k, (*(uidx->bub))));
            // }
            if((ug->u.a[k].n >= rocc_cut) && 
                            ((!(IF_HOM(k, (*(uidx->bub))))) || (asm_opt.purge_level_primary == 0))) {
                if(id_map[i]&p[0]) {
                    ff[i<<1] = 2;
                    // fprintf(stderr, "[M::%s::] utg%.6dl(%c), ug->u.a[k].n::%u, rocc_cut::%u\n", __func__, 
                    // i+1, "+-"[0], (uint32_t)ug->u.a[k].n, rocc_cut);
                }
                if(id_map[i]&p[1]) {
                    ff[(i<<1)+1] = 2;
                    // fprintf(stderr, "[M::%s::] utg%.6dl(%c), ug->u.a[k].n::%u, rocc_cut::%u\n", __func__, 
                    // i+1, "+-"[1], (uint32_t)ug->u.a[k].n, rocc_cut);
                }
            }
            id_map[i] = uidx->bub->index[k];
        }
    }
    // fprintf(stderr, "1[M::%s::] f[51]::%u\n", __func__, ff[51]);
    free(uidx->bub->index); MALLOC(uidx->bub->index, ng->n); 
    memcpy(uidx->bub->index, id_map, sizeof((*id_map))*ng->n);
    if(asm_opt.purge_level_primary == 0) {///all nodes are het
        for (k = 0; k < ug->g->n_seq; k++) {
            if(IF_HOM(k, *(uidx->bub))) uidx->bub->index[k] = uidx->bub->f_bub+1;
        }
    }
    ma_ug_destroy(ug);
}

void gen_hybrid_aln_idx(usg_t *ng, uint64_t *int_idx, uint64_t int_idx_n, uint64_t *int_a, uint64_t int_an, asg64_v *b64, uint64_t *ridx)
{
    uint64_t k, i, l, m, v, s, e, *a, a_n; 
    kv_resize(uint64_t, *b64, int_an); b64->n = int_an; 
    memset(ridx, 0, sizeof((*ridx))*ng->n);
    for (k = 0; k < int_idx_n; k++) {///scan all integer contigs
        s = int_idx[k]>>32; e = s + ((uint32_t)int_idx[k]); assert(e > s + 1);//the length is at least 2
        for (i = s; i < e; i++) ridx[int_a[i]>>1]++;
    }

    for (k = l = 0; k < ng->n; k++) {
        m = ridx[k];
        ridx[k] = l; ridx[k] <<= 32; ridx[k] |= m;
        l += m;
    }

    for (k = 0; k < ng->n; k++) {
        a = b64->a + (ridx[k]>>32); a_n = (uint32_t)ridx[k];
        if(a_n) a[a_n-1] = 0;
    }

    for (k = 0; k < int_idx_n; k++) {
        s = int_idx[k]>>32; e = s + ((uint32_t)int_idx[k]); assert(e > s + 1);//the length is at least 2
        for (i = s; i < e; i++) {
            v = int_a[i]>>1;
            a = b64->a + (ridx[v]>>32); 
            a_n = (uint32_t)ridx[v];
            if(a_n) {
                if(a[a_n-1] == a_n-1) {
                    a[a_n-1] = (k<<32)|i; 
                } else {
                    a[a[a_n-1]++] = (k<<32)|i;
                }
            }
        }
    }
}

void gen_sub_integer_path(uint64_t *int_a, int64_t s, int64_t e, int64_t it, uint64_t v, asg64_v *res)
{
    int64_t k; res->n = 0; assert((int_a[it]>>1) == (v>>1));
    if(int_a[it] == v) {///forward
        kv_resize(uint64_t, *res, (uint64_t)(e - it));
        for (k = it; k < e; k++) res->a[res->n++] = int_a[k];
    } else {//reverse
        kv_resize(uint64_t, *res, (uint64_t)(it - s));
        for (k = it; k >= s; k--) res->a[res->n++] = int_a[k]^1;
    }
}

void prt_sub_integer_path(asg64_v *res, uint64_t it, uint64_t w0, uint64_t w1, uint64_t z_n, uint64_t z)
{
    uint64_t k;
    fprintf(stderr, "(it::%lu) res->n::%u, w0::%lu, w1::%lu, z_n::%lu, z::%lu\n", it, (uint32_t)res->n, 
    w0, w1, z_n, z);
    for (k = 0; k < res->n; k++) {
        fprintf(stderr, "(it::%lu) utg%.6dl,", it, (int32_t)(res->a[k]>>1)+1);
    }
    fprintf(stderr, "\n");
}

uint64_t is_best_path(ul_resolve_t *uidx, usg_t *ng, uint64_t *int_idx, uint64_t *int_a, 
uint64_t s, uint64_t e, uint64_t it, uint64_t v, uint64_t *ridx_a, uint64_t *ridx, 
asg64_v *b0, asg64_v *b1, double cutoff)
{
    uint64_t *arc_a, arc_n, k, is, ie, z, zn, w0, w1, min_w0, min_w1, alt_n = 0, is_contain = 1;
    gen_sub_integer_path(int_a, s, e, it, v, b0);

    arc_a = ridx_a + (ridx[v>>1]>>32); 
    arc_n = (uint32_t)ridx[v>>1];


    // uint64_t is_debug = 0;
    // if((((v>>1) == 2736) && (v&1))/** && (b0->n > 1 && (b0->a[b0->n-1]>>1) == 34944)**/) {
    //     // is_debug = 1;
    //     fprintf(stderr, "[M::%s::] utg%.6dl(%c)\n", __func__, 
    //                 (int32_t)(v>>1)+1, "+-"[v&1]);
    //     prt_sub_integer_path(b0, it, (uint64_t)-1, (uint64_t)-1, (uint64_t)-1, (uint64_t)-1);
    // }


    for (k = 0, min_w0 = min_w1 = (uint64_t)-1, is_contain = 1; k < arc_n; k++) {
        if(it == ((uint32_t)arc_a[k])) continue;
        is = int_idx[arc_a[k]>>32]>>32; alt_n++;
        ie = is + ((uint32_t)(int_idx[arc_a[k]>>32])); 
        gen_sub_integer_path(int_a, is, ie, ((uint32_t)arc_a[k]), v, b1);
        zn = MIN(b0->n, b1->n);
        for (z = 0; z < zn && b0->a[z] == b1->a[z]; z++); ///z: first raw unitig that is different between two paths
        assert(z > 0); w0 = w1 = (uint64_t)-1;
        if(z < b0->n) get_integer_seq_ovlps(uidx, b0->a, b0->n, z - 1, 0, NULL, &w0);
        else return 0;///b0 is contained
        if(z < b1->n) get_integer_seq_ovlps(uidx, b1->a, b1->n, z - 1, 0, NULL, &w1);
        else continue;///b1 is contained
        // if(((v>>1) == 2736) && (v&1)) {
        //     prt_sub_integer_path(b0, it, w0, w1, zn, z);
        //     prt_sub_integer_path(b1, it, w0, w1, zn, z);
        // }
        if(w0 == (uint64_t)-1) w0 = 0;
        if(w1 == (uint64_t)-1) w1 = 0;
        if((w0 <= w1) || (w1 > (w0*cutoff))) return 0;
        // if(b0->n == zn) return 0;///b0 is shorter
        // if(b1->n == zn) continue;///b1 is shorter
        if((min_w0 == (uint64_t)-1) || (z == zn) || (min_w0 > w0) || (min_w0 == w0 && min_w1 < w1)) {
            min_w0 = w0; min_w1 = w1;
        }
        is_contain = 0;
    }

    // if(((v>>1) == 2736) && (v&1)) {
    //     fprintf(stderr, "[M::%s::] min_w0::%lu, min_w1::%lu\n\n", __func__, min_w0, min_w1);
    // }
    // fprintf(stderr, "[M::%s::] utg%.6dl(%c)->utg%.6dl(%c)\n", __func__, 
    //                 (int32_t)(v>>1)+1, "+-"[v&1], (int32_t)(int_a[k]>>1)+1, "+-"[int_a[k]&1]);
    if(alt_n == 0) return 1;
    if(is_contain) return 1;
    
    if(min_w0 == (uint64_t)-1) min_w0 = 0;
    if(min_w1 == (uint64_t)-1) min_w1 = 0;
    if((min_w0 > min_w1) && (min_w1 <= (min_w0*cutoff))) return 1;
    return 0;
}


void get_best_path(ul_resolve_t *uidx, usg_t *ng, uint64_t *int_idx, uint64_t *int_a, 
uint8_t *f, uint64_t s, uint64_t e, asg64_v *b0, asg64_v *b1, uint64_t *ridx_a, 
uint64_t *ridx, asg64_v *res)
{
    uint64_t pi, v, k, res_n = res->n, *pz;
    b0->n = b1->n = 0; 
    for (k = s, pi = (uint64_t)-1; k < e; k++) {
        v = int_a[k]; 
        if((!f[v])&&(!f[v^1])) continue;
        if((pi != (uint64_t)-1) && (f[v^1])) {
            // if(((v>>1) == 2736)) {
            //     fprintf(stderr, "+[M::%s::ii[%lu, %lu)] utg%.6dl(%c), f[v^1]::%u, v^1::%lu, putg%.6dl(%c), f[pv]::%u, pv::%lu\n", 
            //     __func__, s, e, (int32_t)(int_a[k]>>1)+1, "+-"[int_a[k]&1], f[v^1], v^1, 
            //     (int32_t)(int_a[pi]>>1)+1, "+-"[int_a[pi]&1], f[int_a[pi]], int_a[pi]);
            // }
            if(is_best_path(uidx, ng, int_idx, int_a, s, e, k, v^1, ridx_a, ridx, b0, b1, 0.51)) {
                // if(((v>>1) == 2736)) {
                //     fprintf(stderr, "-[M::%s::ii[%lu, %lu)] utg%.6dl(%c), f[v^1]::%u, v^1::%lu, putg%.6dl(%c), f[pv]::%u, pv::%lu\n", 
                //     __func__, s, e, (int32_t)(int_a[k]>>1)+1, "+-"[int_a[k]&1], f[v^1], v^1, 
                //     (int32_t)(int_a[pi]>>1)+1, "+-"[int_a[pi]&1], f[int_a[pi]], int_a[pi]);
                // }
                // fprintf(stderr, "[M::%s::] utg%.6dl(%c)->utg%.6dl(%c)\n", __func__, 
                //     (int32_t)(int_a[pi]>>1)+1, "+-"[int_a[pi]&1], (int32_t)(int_a[k]>>1)+1, "+-"[int_a[k]&1]);
                pz = (res->n > res_n)? &(res->a[res->n-1]):(NULL);
                if((!pz) || (((*pz)>>32) != pi)) {//keep the shortest pi<->k 
                    kv_pushp(uint64_t, *res, &pz); 
                    *pz = pi; (*pz) <<= 32; (*pz) |= k;
                }
            } else {
                pi = (uint64_t)-1;
            }            
        }
        if(f[v]) {
            // fprintf(stderr, "-[M::%s::] utg%.6dl(%c), f[v]::%u, v::%lu\n", 
            // __func__, (int32_t)(int_a[k]>>1)+1, "+-"[int_a[k]&1], f[v], v);
            if(is_best_path(uidx, ng, int_idx, int_a, s, e, k, v, ridx_a, ridx, b0, b1, 0.51)) {
                pi = k;//keep the shortest pi<->k
            } else {
                pi = (uint64_t)-1;
            }
        }
    }
}


uint64_t gen_sub_integer_path_ff(uint64_t *int_a, int64_t s, int64_t e, int64_t it, uint64_t v, uint8_t *f, asg64_v *res)
{
    int64_t k; res->n = 0; assert((int_a[it]>>1) == (v>>1));
    if(int_a[it] == v) {///forward
        kv_resize(uint64_t, *res, (uint64_t)(e - it));
        for (k = it; k < e; k++) {
            res->a[res->n++] = int_a[k];
            if(res->n > 1 && f[res->a[res->n-1]^1]) return res->n;
        }
    } else {//reverse
        kv_resize(uint64_t, *res, (uint64_t)(it - s));
        for (k = it; k >= s; k--) {
            res->a[res->n++] = int_a[k]^1;
            if(res->n > 1 && f[res->a[res->n-1]^1]) return res->n;
        }
    }
    res->n = 0; return 0;
}

uint64_t is_best_pair(ul_resolve_t *uidx, usg_t *ng, uint64_t *int_idx, uint64_t *int_a, 
uint64_t s, uint64_t e, uint64_t it, uint64_t v, uint64_t *ridx_a, uint64_t *ridx, 
asg64_v *b0, asg64_v *b1, uint8_t *f, double cutoff)
{
    uint64_t *arc_a, arc_n, k, is, ie, w0 = 0, w1 = 0, sup_cut = 3; 
    if(!gen_sub_integer_path_ff(int_a, s, e, it, v, f, b0)) return 0;
    assert(b0->n >= 2); assert(f[b0->a[0]] && f[b0->a[b0->n-1]^1]);
    w0 = get_arc_support_chain(uidx, b0->a, b0->n, ng);
    if(w0 < sup_cut) return 0;

    arc_a = ridx_a + (ridx[v>>1]>>32); 
    arc_n = (uint32_t)ridx[v>>1];
    for (k = 0; k < arc_n; k++) {
        if(it == ((uint32_t)arc_a[k])) continue;
        is = int_idx[arc_a[k]>>32]>>32; 
        ie = is + ((uint32_t)(int_idx[arc_a[k]>>32])); 
        if(!gen_sub_integer_path_ff(int_a, is, ie, ((uint32_t)arc_a[k]), v, f, b1)) continue;
        assert(b1->n >= 2); assert(f[b1->a[0]] && f[b1->a[b1->n-1]^1]);
        if((b0->n == b1->n) && (memcmp(b0->a, b1->a, sizeof((*(b0->a)))*b0->n) == 0)) {
            if(it < ((uint32_t)arc_a[k])) continue;
            else return 0;///only keep 1 equal interval
        }
        w1 = get_arc_support_chain(uidx, b1->a, b1->n, ng);
        // fprintf(stderr, "[M::%s::] w0::%lu, w1::%lu\n", __func__, w0, w1);
        if((w0 <= w1) || (w1 > (w0*cutoff))) return 0;
    }
    return 1;
}

void get_best_pair(ul_resolve_t *uidx, usg_t *ng, uint64_t *int_idx, uint64_t *int_a, 
uint8_t *f, uint64_t s, uint64_t e, asg64_v *b0, asg64_v *b1, uint64_t *ridx_a, 
uint64_t *ridx, asg64_v *res)
{
    uint64_t pi, v, k, *pz;
    b0->n = b1->n = 0; 
    for (k = s, pi = (uint64_t)-1; k < e; k++) {
        v = int_a[k]; 
        if((f[v^1])) {
            if(pi != (uint64_t)-1) {
                if(is_best_pair(uidx, ng, int_idx, int_a, s, e, k, v^1, ridx_a, ridx, b0, b1, f, 0.51)) {
                    kv_pushp(uint64_t, *res, &pz); 
                    *pz = pi; (*pz) <<= 32; (*pz) |= k;
                } 
            }
            pi = (uint64_t)-1;    
        }
        if(f[v]) {
            pi = (uint64_t)-1;
            if(is_best_pair(uidx, ng, int_idx, int_a, s, e, k, v, ridx_a, ridx, b0, b1, f, 0.51)) {
                pi = k;//keep the shortest pi<->k
            } 
        }
    }
}

static void worker_ul_aln_path(void *data, long i, int tid) // callback for kt_for()
{
    unique_bridge_check_t *u_aux = (unique_bridge_check_t*)data;
    ul_resolve_t *uidx = u_aux->uidx; 
    integer_t *buf = &(uidx->str_b.buf[tid]);
    uint64_t s, e;
    // uint64_t *x = &(uidx->uovl.iug_tra->a[i]); 
    // asg_arc_t *ve = &(uidx->uovl.i_ug->g->arc[*x]);

    asg64_v b_v, b_r, res; 
    b_v.a = buf->u.a; b_v.n = buf->u.n; b_v.m = buf->u.m; 
    b_r.a = buf->o.a; b_r.n = buf->o.n; b_r.m = buf->o.m; 
    res.a = buf->res_dump.a; res.n = buf->res_dump.n; res.m = buf->res_dump.m;

    b_v.n = b_r.n = 0; 
    s = u_aux->int_idx[i]>>32; ///the i-th integer contig/path 
    e = s + ((uint32_t)(u_aux->int_idx[i])); 
    assert(e > s + 1);//the length is at least 2

    get_best_path(uidx, u_aux->ng, u_aux->int_idx, u_aux->int_a, u_aux->f, s, e, &b_v, &b_r, u_aux->ridx_a, u_aux->ridx, &res);
    // get_ul_arc_supports(uidx, ve, &b_v, &b_r, 1, &w_v, &w_r);

    buf->u.a = b_v.a; buf->u.n = b_v.n; buf->u.m = b_v.m; 
    buf->o.a = b_r.a; buf->o.n = b_r.n; buf->o.m = b_r.m; 
    buf->res_dump.a = res.a; buf->res_dump.n = res.n; buf->res_dump.m = res.m;
}

static void worker_ul_aln_pair(void *data, long i, int tid) // callback for kt_for()
{
    unique_bridge_check_t *u_aux = (unique_bridge_check_t*)data;
    ul_resolve_t *uidx = u_aux->uidx; 
    integer_t *buf = &(uidx->str_b.buf[tid]);
    uint64_t s, e;
    // uint64_t *x = &(uidx->uovl.iug_tra->a[i]); 
    // asg_arc_t *ve = &(uidx->uovl.i_ug->g->arc[*x]);

    asg64_v b_v, b_r, res; 
    b_v.a = buf->u.a; b_v.n = buf->u.n; b_v.m = buf->u.m; 
    b_r.a = buf->o.a; b_r.n = buf->o.n; b_r.m = buf->o.m; 
    res.a = buf->res_dump.a; res.n = buf->res_dump.n; res.m = buf->res_dump.m;

    b_v.n = b_r.n = 0; 
    s = u_aux->int_idx[i]>>32; ///the i-th integer contig/path 
    e = s + ((uint32_t)(u_aux->int_idx[i])); 
    assert(e > s + 1);//the length is at least 2

    get_best_pair(uidx, u_aux->ng, u_aux->int_idx, u_aux->int_a, u_aux->f, s, e, &b_v, &b_r, u_aux->ridx_a, u_aux->ridx, &res);
    // get_ul_arc_supports(uidx, ve, &b_v, &b_r, 1, &w_v, &w_r);

    buf->u.a = b_v.a; buf->u.n = b_v.n; buf->u.m = b_v.m; 
    buf->o.a = b_r.a; buf->o.n = b_r.n; buf->o.m = b_r.m; 
    buf->res_dump.a = res.a; buf->res_dump.n = res.n; buf->res_dump.m = res.m;
}

uint32_t ug_ext_0(ul_resolve_t *uidx, usg_t *ng, uint64_t *int_idx, uint64_t int_idx_n, uint64_t *int_a, uint32_t max_ext, 
uint8_t *ff, uint32_t *ng_occ, uint64_t *i_idx, asg64_v *b64, asg64_v *ub64, uint64_t a_n)
{
    uint32_t n_vtx = ng->n<<1, k, i; uint64_t *pz;
    uint64_t s, e, n_clus;
    memset(ff, 0, sizeof((*ff))*n_vtx);
    // memset(ng_occ, 0, sizeof((*ng_occ))*n_vtx); 
    memset(i_idx, 0, sizeof((*i_idx))*n_vtx); 
    ub64->n = 0;

    for (i = 0; i < a_n; i++) { ///available intervals
        s = b64->a[i]>>32; e = (uint32_t)b64->a[i]; assert(e > s);///[s, e] -> coordinates within int_a[]
        assert(!(b64->a[i]&((uint64_t)0x8000000000000000)));
        // fprintf(stderr, "\n[M::%s::] occ::%lu\n", __func__, e - s);
        // for (k = s; k <= e; k++) {
        //     // fprintf(stderr, "utg%.6dl(%c),", (int32_t)(int_a[k]>>1)+1, "+-"[int_a[k]&1]);
        //     fprintf(stderr, "utg%.6dl,", (int32_t)(int_a[k]>>1)+1);
        // }
        // fprintf(stderr, "\n");

        for (k = s + 1; k < e; k++) {///note: here is [s, e]; s && e are unique, but [s+1, e-1] are not unique
            kv_pushp(uint64_t, *b64, &pz); //i_idx[int_a[k]]++;
            (*pz) = int_a[k]; (*pz) <<= 32; (*pz) |= i;///(raw unitig/non-unqiue node id)|(integer contig id)

            kv_pushp(uint64_t, *b64, &pz); //i_idx[int_a[k]^1]++;
            (*pz) = int_a[k]^1; (*pz) <<= 32; (*pz) |= i;
        }
        kv_pushp(uint64_t, *b64, &pz); //i_idx[int_a[s]]++;
        (*pz) = int_a[s]; (*pz) <<= 32; (*pz) |= i;

        kv_pushp(uint64_t, *b64, &pz); //i_idx[int_a[e]^1]++;
        (*pz) = int_a[e]^1; (*pz) <<= 32; (*pz) |= i;
    }
    // prt_thread_info(b64.a, a_n, NULL, 0, "tt2");
    // fprintf(stderr, "**1**[M::%s::] a_n::%lu, ua_n::%lu\n", __func__, a_n, ua_n);
    ///index
    radix_sort_srt64(b64->a + a_n, b64->a + b64->n);///(raw unitig node id)|(integer contig id)
    for (k = a_n + 1, i = a_n; k <= b64->n; k++) {
        if(k == b64->n || (b64->a[k]>>32) != (b64->a[i]>>32)) {
            i_idx[b64->a[i]>>32] |= (((uint64_t)i)<<32)|((uint64_t)k);///b64.a[i]>>32 appear once (unique ends)/multipe times
            i = k;
        }
    }

    n_clus = usg_unique_arcs_cluster_adv(b64, a_n, i_idx, int_a, ub64);
    // fprintf(stderr, "**3**[M::%s::] a_n::%lu, n_clus::%lu\n", __func__, a_n, n_clus);
    // prt_thread_info(b64.a, a_n, NULL, 0, "tt4");
    assert(b64->n == (a_n<<1));
    
    // new_path_ext(uidx, ng, &b64, &ub64, a_n, n_clus, i_idx, int_a, ff, ng_occ, max_ext);//wrong
    return old_path_ext(uidx, ng, b64, ub64, a_n, n_clus, i_idx, int_a, ff, ng_occ, max_ext);
}

void debug_prt_renew_aln(usg_t *ng, uint64_t *int_idx, uint64_t int_idx_n, uint64_t *int_a, uint64_t int_an, asg64_v *b64, uint64_t *ridx)
{
    uint64_t k, s, e, i;
    for (i = 0; i < int_idx_n; i++) {///scan all integer contigs
        s = int_idx[i]>>32; e = s + ((uint32_t)int_idx[i]); assert(e > s + 1);//the length is at least 2
        fprintf(stderr, "\n[M::%s::] occ::%lu\n", __func__, e - s);
        for (k = s; k <= e; k++) {
            // fprintf(stderr, "utg%.6dl(%c),", (int32_t)(int_a[k]>>1)+1, "+-"[int_a[k]&1]);
            fprintf(stderr, "utg%.6dl,", (int32_t)(int_a[k]>>1)+1);
        }
        fprintf(stderr, "\n");
    }
}

uint32_t ug_ext_free(asg64_v *ob, asg64_v *ub, ul_resolve_t *uidx, usg_t *ng, uint32_t max_ext, 
uint8_t **ff, uint32_t **ng_occ, uint64_t **i_idx, asg64_v *b64, asg64_v *ub64, uint32_t rocc_cut, 
uint32_t thread_path)
{
    // fprintf(stderr, "[M::%s::] rocc_cut::%u\n", __func__, rocc_cut);
    uint32_t k, n_vtx = ng->n<<1, a_n; unique_bridge_check_t u_aux;
    REALLOC((*ff), n_vtx); REALLOC((*ng_occ), n_vtx); REALLOC((*i_idx), n_vtx);
    renew_usg_t_bub(uidx, ng, *ng_occ, *ff, rocc_cut);
    u2g_hybrid_aln(uidx, ng, ob, ub);
    gen_hybrid_aln_idx(ng, ub->a, ub->n, ob->a, ob->n, b64, *i_idx);
    // if(rocc_cut == 10) 
    // {
    //     fprintf(stderr, "\n[M::%s::] rocc_cut::%u\n", __func__, rocc_cut);
    //     debug_prt_renew_aln(ng, ub->a, ub->n, ob->a, ob->n, b64, *i_idx);
    // }

    u_aux.f = *ff; u_aux.ng = ng; u_aux.int_idx = ub->a; u_aux.int_idx_n = ub->n; u_aux.uidx = uidx;
    u_aux.int_a = ob->a; u_aux.int_an = ob->n; u_aux.ridx_a = b64->a; u_aux.ridx = *i_idx;
    for (k = 0; k < uidx->str_b.n_thread; k++) {
        uidx->str_b.buf[k].res_dump.n = uidx->str_b.buf[k].u.n = uidx->str_b.buf[k].o.n = 0;
    }

    if(thread_path) {
        kt_for(uidx->str_b.n_thread, worker_ul_aln_path, &u_aux, u_aux.int_idx_n);
    } else {
        kt_for(uidx->str_b.n_thread, worker_ul_aln_pair, &u_aux, u_aux.int_idx_n);
    }

    for (k = b64->n = a_n = 0; k < uidx->str_b.n_thread; k++) {
        a_n += uidx->str_b.buf[k].res_dump.n;
        kv_resize(uint64_t, *b64, a_n);
        memcpy(b64->a+b64->n, uidx->str_b.buf[k].res_dump.a, uidx->str_b.buf[k].res_dump.n*(sizeof(*(b64->a))));
        b64->n = a_n;
    }
    radix_sort_srt64(b64->a, b64->a + b64->n);///keeps the coordinates within int_a[]
    // fprintf(stderr, "[M::%s::] a_n::%u, int_idx_n::%u, int_n::%u\n", 
    //                 __func__, a_n, (uint32_t)ub->n, (uint32_t)ob->n);
    a_n = ug_ext_0(uidx, ng, ub->a, ub->n, ob->a, max_ext, *ff, *ng_occ, *i_idx, b64, ub64, a_n);
    // if(a_n) usg_cleanup(ng);
    return a_n;
}

uint32_t ug_ext_strict(asg64_v *ob, asg64_v *ub, ul_resolve_t *uidx, usg_t *ng, uint32_t max_ext, 
uint8_t **ff, uint32_t **ng_occ, uint64_t **i_idx, asg64_v *b64, asg64_v *ub64)
{
    // fprintf(stderr, "[M::%s::]\n", __func__);
    uint32_t k, i, n_vtx = ng->n<<1, a_n;
    REALLOC((*ff), n_vtx); REALLOC((*ng_occ), n_vtx); REALLOC((*i_idx), n_vtx);
    u2g_hybrid_aln(uidx, ng, ob, ub); b64->n = ub64->n = 0; 
    uint64_t *int_idx = ub->a, int_idx_n = ub->n, *int_a = ob->a; 
    uint64_t *r_a, r_n, zs, ze, s, e, pi, v, *pz;

    memset((*ng_occ), 0, sizeof((*(*ng_occ)))*n_vtx);  
    for (k = 0; k < int_idx_n; k++) {///scan all integer contigs
        r_a = int_a + (int_idx[k]>>32); r_n = (uint32_t)int_idx[k];
        assert(r_n >= 2);
        for (i = 1; i + 1 < r_n; i++) {
            (*ng_occ)[r_a[i]]++; (*ng_occ)[r_a[i]^1]++;
        }
        (*ng_occ)[r_a[0]]++; (*ng_occ)[r_a[r_n-1]^1]++;
    }

    ma_ug_t *un_g = ma_ug_hybrid_gen(ng); 
    int32_t ui, un; ma_utg_t *u = NULL;
    for (k = 0; k < un_g->u.n; k++) {///all unitigs of raw utg 
        u = &(un_g->u.a[k]); 
        zs = ze = (uint32_t)-1; un = u->n;
        for (ui = 0; ui < un; ui++) {
            if(occ_m((*ng_occ)[(u->a[ui]>>32)^1]) == 1) {
                zs = ui; break;
            }
        }

        for (ui = ((int32_t)un)-1; ui >= 0; ui--) {
            if(occ_m((*ng_occ)[u->a[ui]>>32]) == 1) {
                ze = ui; break;
            }
        }
        if(zs != (uint32_t)-1 && ze != (uint32_t)-1 && zs > ze) continue;
        if(zs != (uint32_t)-1) (*ng_occ)[(u->a[zs]>>32)^1] |= ((uint32_t)0x80000000);
        if(ze != (uint32_t)-1) (*ng_occ)[(u->a[ze]>>32)] |= ((uint32_t)0x80000000);
    }
    ma_ug_destroy(un_g);

        // fprintf(stderr, ">>>>>>[M::%s::] int_idx[0]::%lu, int_idx[1]::%lu\n", __func__, int_idx[0], int_idx[1]);
    // prt_intg_info(int_idx, int_idx_n, int_a, "intg");
    for (i = b64->n = a_n = 0; i < int_idx_n; i++) {///scan all integer contigs
        s = int_idx[i]>>32; e = s + ((uint32_t)int_idx[i]); 
        assert(e > s + 1);//the length is at least 2
        for (k = s, pi = (uint64_t)-1; k < e; k++) {
            v = int_a[k]; 
            if((!((*ng_occ)[v]&((uint32_t)0x80000000)))&&(!((*ng_occ)[v^1]&((uint32_t)0x80000000)))) {
                continue;///it must be a unique node
            }
            if((pi != (uint64_t)-1) && ((*ng_occ)[v^1]&((uint32_t)0x80000000))) {
                pz = (b64->n>0)? &(b64->a[b64->n-1]):(NULL);
                if((!pz) || (((*pz)>>32) != pi)) {//keep the shortest pi<->k 
                    kv_pushp(uint64_t, *b64, &pz); 
                    *pz = pi; (*pz) <<= 32; (*pz) |= k;
                    a_n++;
                }
            }
            if((*ng_occ)[v]&((uint32_t)0x80000000)) pi = k;//keep the shortest pi<->k
        }
    }
    // prt_thread_info(b64.a, a_n+ua_n, NULL, 0, "tt_minus");
    // fprintf(stderr, "**0**[M::%s::] a_n::%lu, ua_n::%lu\n", __func__, a_n, ua_n);
    radix_sort_srt64(b64->a, b64->a + b64->n);///keeps the coordinates within int_a[]

    a_n = ug_ext_0(uidx, ng, ub->a, ub->n, ob->a, max_ext, *ff, *ng_occ, *i_idx, b64, ub64, a_n);
    // if(a_n) usg_cleanup(ng);
    return a_n;
}

typedef struct {
    uint32_t nid, ulid;
    uint32_t raw_sid, raw_eid;
    uint32_t raw_sof, raw_eof;
} usc_t;

typedef struct {
    usc_t *a;
    size_t n, m;
} usc_vec_t;

typedef struct {
    // ul_resolve_t *uidx;
    usc_t *a;
    uint64_t a_n;
    ma_ug_t *ug;
} scaf_mul_t;



#define B4Lg(x) (((x)>>2)+(((x)&3)?1:0))
uint32_t load_scaf_base(all_ul_t *x, char* file_name, const char *bin_file)
{
    char *gfa_name = (char*)malloc(strlen(file_name)+50);
	sprintf(gfa_name, "%s.%s.uidx.ucr.bin", file_name, bin_file);
    // fprintf(stderr, "[M::%s] open %s...\n", __func__, gfa_name);
    FILE *fp = fopen(gfa_name, "r"); free(gfa_name);
    if (!fp) return 0;
    // fprintf(stderr, "[M::%s] open sucess\n", __func__);
    uint64_t rid; uint32_t len; ul_vec_t ss, *z; memset(&ss, 0, sizeof(ss));
    while(1) {
        fread(&rid, sizeof(rid), 1, fp);
        if(feof(fp)) break;
        z = &(x->a[rid]);
        fread(&len, sizeof(len), 1, fp); assert(z->rlen == len);
        fread(&(ss.N_site.n), sizeof(ss.N_site.n), 1, fp);
        kv_resize(uint32_t, ss.N_site, ss.N_site.n);
        fread(ss.N_site.a, sizeof((*(ss.N_site.a))), ss.N_site.n, fp);
        ss.r_base.n = B4Lg(len); kv_resize(uint8_t, ss.r_base, ss.r_base.n); 
        fread(ss.r_base.a, sizeof((*(ss.r_base.a))), ss.r_base.n, fp);
        // if(rid == 37238) {
        //     fprintf(stderr, "[M::%s::37238] z->dd::%u\n", __func__, z->dd);
        // }
        // fprintf(stderr, "[M::%s::rid->%lu] z->dd::%u\n", __func__, rid, z->dd);
        if(z->dd != 4) continue;

        
        kv_resize(uint32_t, z->N_site, ss.N_site.n); z->N_site.n = ss.N_site.n;
        memcpy(z->N_site.a, ss.N_site.a, sizeof((*(z->N_site.a)))*z->N_site.n);

        kv_resize(uint8_t, z->r_base, ss.r_base.n); z->r_base.n = ss.r_base.n;
        memcpy(z->r_base.a, ss.r_base.a, sizeof((*(z->r_base.a)))*z->r_base.n);
    }
    // load_compress_base_disk(fp, &rid, des.a, &ulen, &(sl->ucr_s->u));

    fclose(fp); free(ss.N_site.a); free(ss.r_base.a);
    return 1;
}   


uint64_t reset_scaf_node_uinfo_srt_t(usc_t *z, int64_t min_arc_len, int64_t scaf_len, int64_t *nlen)
{
    (*nlen) = scaf_len + (min_arc_len<<1);
    if(z->ulid == ((uint32_t)-1)) return 0;///a scaffold node
    if(z->raw_eof + min_arc_len <= z->raw_sof) {///has an overlap longer than min_arc_len; 
        return 1;///a scaffold node; no need this node, could directly use existing nodes
    }
    int64_t s, e;
    s = z->raw_sof; e = z->raw_eof; 
    (*nlen) = (min_arc_len<<1) + (e - s);  ///(e - s)>=(-min_arc_len)
    return 2;
}


void get_end_hifi(char *des, uint32_t v, All_reads* rdb, uint32_t len, uint32_t is_beg)
{
    uint32_t rlen = Get_READ_LENGTH((*rdb), (v>>1));
    assert(rlen >= len);
    recover_UC_Read_sub_region(des, ((is_beg)?(rlen-len):(0)), len, v&1, rdb, v>>1);
}

void get_ul_subregion(all_ul_t *x, char *des, uint32_t id, uint32_t rev, int64_t ssp, int64_t sep, uint32_t reset_Ns)
{
    ul_vec_t *z = &(x->a[id]); 
    int64_t sl, slr, offset, begLen, tailLen, a_n, src_i, des_i, i;
    sl = sep - ssp;
	offset = ssp&3; begLen = 4-offset;
	if(begLen > sl) begLen = sl;
	tailLen = (sl-begLen)&3;
	a_n = (sl - begLen - tailLen)>>2;

    // fprintf(stderr, "[M::%s] z->r_base.n::%u, begLen::%ld, tailLen::%ld, a_n::%ld\n", 
    // __func__, (uint32_t)z->r_base.n, begLen, tailLen, a_n);
    src_i = ssp; i = 0; des_i = 0;
	if(begLen > 0) {
        // fprintf(stderr, "[M::%s] des_i::%ld, src_i>>2::%ld\n", __func__, des_i, src_i>>2);
		memcpy(des+des_i, bit_t_seq_table[z->r_base.a[src_i>>2]]+offset, begLen);
		des_i += begLen; src_i += begLen;
	}

	for (i = 0; i < a_n; i++) {
        // fprintf(stderr, "[M::%s] des_i::%ld, src_i>>2::%ld\n", __func__, des_i, src_i>>2);
		memcpy(des+des_i, bit_t_seq_table[z->r_base.a[src_i>>2]], 4);
		des_i += 4; src_i += 4;
	}

	if(tailLen > 0) {
        // fprintf(stderr, "[M::%s] des_i::%ld, src_i>>2::%ld\n", __func__, des_i, src_i>>2);
		memcpy(des+des_i, bit_t_seq_table[z->r_base.a[src_i>>2]], tailLen);
		des_i += tailLen; src_i += tailLen;
	}

    uint64_t k, sk = ssp, ek = sep;
    for (k = 0; k < z->N_site.n; k++) {
        if(z->N_site.a[k] >= sk && z->N_site.a[k] < ek){
            des[z->N_site.a[k]-sk] = 'N';
        } else if(z->N_site.a[k] >= ek) {
            break;
        }
    }
    
    if(reset_Ns) {
        for (i = 0; i < sl; i++) {
            if (seq_nt4_table[(uint8_t)des[i]] >= 4) des[i] = 'A';
        }
    }

    if(rev) {
        char t; slr = sl>>1;
        for (i = 0; i < slr; i++) {
            t = des[sl-i-1];
            des[sl-i-1] = RC_CHAR(des[i]);
            des[i] = RC_CHAR(t);
        }
        if(sl&1) des[i] = RC_CHAR(des[i]);
    }
}

int64_t push_scaf_bases(uint8_t *des, uint64_t** N_site, All_reads* rdb, all_ul_t *x, usc_t *z, int64_t min_arc_len, int64_t scaf_len, int64_t nlen, UC_Read *tu, 
uint64_t *rmap, ma_ug_t *raw_g)
{   
    int64_t nlen0, ff, Nocc, i; char *da = NULL;
    ff = reset_scaf_node_uinfo_srt_t(z, min_arc_len, scaf_len, &nlen0); assert(nlen0 == nlen);
    if(ff == 1) {///has an overlap longer than min_arc_len between z->raw_sid and z->raw_eid
        memset(des, 0, sizeof((*(des)))*(nlen/4+1));
        return ff;
    }
    uint32_t uv, uw, sv, sw; int64_t s, e, ul;
    uv = z->raw_sid; uw = z->raw_eid;

    sv = (uv&1?((raw_g->u.a[uv>>1].a[0]>>32)^1):(raw_g->u.a[uv>>1].a[raw_g->u.a[uv>>1].n-1]>>32));
    sw = (uw&1?((raw_g->u.a[uw>>1].a[raw_g->u.a[uw>>1].n-1]>>32)^1):(raw_g->u.a[uw>>1].a[0]>>32));
    resize_UC_Read(tu, nlen); da = tu->seq;
    assert((Get_READ_LENGTH((*rdb), (sv>>1))) > ((uint32_t)min_arc_len));
    assert((Get_READ_LENGTH((*rdb), (sw>>1))) > ((uint32_t)min_arc_len));
    ///need introduce some UL bases; 
    ///or the overlap length between sv and sw is shorter than min_arc_len
    if(ff == 2) {
        s = z->raw_sof; e = z->raw_eof; ul = e - s;
        assert(ul>=(-min_arc_len));
        // assert((Get_READ_LENGTH((*rdb), (sv>>1))) > ul);
        // assert((Get_READ_LENGTH((*rdb), (sw>>1))) > ul);
        get_end_hifi(da, sv, rdb, min_arc_len, 1);
        get_end_hifi(da+min_arc_len+ul, sw, rdb, min_arc_len, 0);    
        if(ul > 0) {///need UL
            // fprintf(stderr, "[M::%s] ulid::%u(%c), s::%ld, e::%ld, nlen::%ld, min_arc_len::%ld\n", __func__, 
            // z->ulid>>1, "+-"[z->ulid&1], s, e, nlen, min_arc_len);
            get_ul_subregion(x, da+min_arc_len, z->ulid>>1, z->ulid&1, s, e, 1);
        }
    } else {///no UL found; ff = 0
        get_end_hifi(da, sv, rdb, min_arc_len, 1);
        memset(da+min_arc_len, 'A', scaf_len);
        get_end_hifi(da+min_arc_len+scaf_len, sw, rdb, min_arc_len, 0);        
    }
    for (i = Nocc = 0; i < nlen; i++) {
        if(seq_nt6_table[(uint8_t)da[i]] >= 4) Nocc++;
    }
    ha_compress_base(des, da, nlen, N_site, Nocc);
    return ff;
}

void realloc_rdb_adv(All_reads* rdb, all_ul_t *x, ma_sub_t **cov, R_to_U *ruI, uint64_t *rmap, 
uint64_t rid_n, uint64_t scaf_len, char *scaf_id, asg_t *ng, ug_opt_t *uopt, ma_ug_t *raw_g, usc_t *a)
{
    uint64_t i, rid_n0 = rdb->total_reads, tname, cname; usc_t *z; 
    uint64_t scaf_id_len = strlen(scaf_id); char *des, *src; int64_t nlen, ff;
    UC_Read g_read; init_UC_Read(&g_read);
    rdb->total_reads = rid_n;
    tname = rdb->name_index[rid_n0];
    fprintf(stderr, "+[M::%s] rid_n0::%lu, rid_n::%lu\n", __func__, rid_n0, rid_n);
    ///for read bases
    REALLOC(rdb->N_site, rdb->total_reads);
    REALLOC(rdb->read_length, rdb->total_reads);
    REALLOC(rdb->read_size, rdb->total_reads);
    REALLOC(rdb->read_sperate, rdb->total_reads);
    REALLOC(rdb->trio_flag, rdb->total_reads);
    REALLOC(rdb->name_index, rdb->total_reads+1);///total_reads+1
    REALLOC((*cov), rdb->total_reads);
    if(uopt->te) {
        assert(uopt->te->n == rid_n0);
        uopt->te->n = rdb->total_reads;
        REALLOC(uopt->te->hh, rdb->total_reads);
    }
    for (i = rid_n0; i < rdb->total_reads; i++) {
        // fprintf(stderr, "[M::%s] i::%lu\n", __func__, i);
        rdb->N_site[i] = NULL;
        rdb->trio_flag[i] = AMBIGU;
        (*cov)[i].c = (*cov)[i].del = 0; 
        // rdb->read_length[i] = scaf_len;
        // rdb->read_size[i] = scaf_len;
        // (*cov)[i].s = 0; (*cov)[i].e = scaf_len;
        // cname = scaf_id_len;
        if(rmap[i] != ((uint64_t)-1)) {///not a scaffold node
            if(rdb->N_site[rmap[i]] != NULL) {
                MALLOC(rdb->N_site[i], rdb->N_site[rmap[i]][0]+1);
                memcpy(rdb->N_site[i], rdb->N_site[rmap[i]], 
                    sizeof((*(rdb->N_site[i])))*(rdb->N_site[rmap[i]][0]+1));
            }
            rdb->read_length[i] = rdb->read_length[rmap[i]];
            rdb->read_size[i] = rdb->read_length[rmap[i]];
            (*cov)[i].s = (*cov)[rmap[i]].s; 
            (*cov)[i].e = (*cov)[rmap[i]].e;
            rdb->trio_flag[i] = rdb->trio_flag[rmap[i]];
            (*cov)[i] = (*cov)[rmap[i]];
            cname = Get_NAME_LENGTH((*rdb), (rmap[i]));
            if(uopt->te) uopt->te->hh[i] = uopt->te->hh[rmap[i]];
        } else {
            z = &(a[(i-ng->r_seq)]); assert(z->nid == i); cname = scaf_id_len;
            if(reset_scaf_node_uinfo_srt_t(z, uopt->min_ovlp, scaf_len, &nlen) == 2) {
                cname = UL_INF.nid.a[z->ulid>>1].n;///ul name length
            }
            rdb->read_length[i] = nlen;
            rdb->read_size[i] = nlen;
            (*cov)[i].s = 0; (*cov)[i].e = nlen;
            ng->seq[i].len = nlen;///update length for the scaffold node
            if(uopt->te) uopt->te->hh[i] = 0;
        } 
        rdb->name_index[i] = tname; tname += cname;
        rdb->total_reads_bases += rdb->read_length[i];

        MALLOC(rdb->read_sperate[i], (rdb->read_length[i]/4+1));
        if(rmap[i] != ((uint64_t)-1)) {///not a scaffold node
            memcpy(rdb->read_sperate[i], rdb->read_sperate[rmap[i]], 
                sizeof((*(rdb->read_sperate[i])))*(rdb->read_length[i]/4+1));
        } else {
            ///set to A
            z = &(a[(i-ng->r_seq)]); assert(z->nid == i);
            ff = push_scaf_bases(rdb->read_sperate[i], &(rdb->N_site[i]), rdb, x, z, uopt->min_ovlp, scaf_len, rdb->read_length[i], &g_read, rmap, raw_g);
            if(ff == 1) ng->seq[i].del = 1;///could directly reuse HiFi; no need new node
            // memset(rdb->read_sperate[i], 0, sizeof((*(rdb->read_sperate[i])))*(rdb->read_length[i]/4+1));
        }
        ng->seq[i].len = rdb->read_length[i];
    }

    rdb->index_size = rdb->total_reads;
    rdb->name_index[i] = tname; 
    rdb->total_name_length = tname; 
    rdb->name_index_size = rdb->total_reads+1;
    REALLOC(rdb->name, tname);
    for (i = rid_n0; i < rdb->total_reads; i++) {///only need to update names for new reads
        des = Get_NAME((*rdb), i); 
        if(rmap[i] != ((uint64_t)-1)) {
            src = Get_NAME((*rdb), rmap[i]); cname = Get_NAME_LENGTH((*rdb), (rmap[i]));
        } else {
            z = &(a[(i-ng->r_seq)]); assert(z->nid == i);
            if(reset_scaf_node_uinfo_srt_t(z, uopt->min_ovlp, scaf_len, &nlen) == 2) {
                cname = UL_INF.nid.a[z->ulid>>1].n;///ul name length
                src = UL_INF.nid.a[z->ulid>>1].a;
            } else {
                src = scaf_id; cname = scaf_id_len;
            }
            // assert(nlen == (int64_t)cname);
        }
        memcpy(des, src, sizeof((*(des)))*cname);
    }

    REALLOC(rdb->paf, rdb->total_reads); 
    memset(rdb->paf+rid_n0, 0, (rdb->total_reads-rid_n0)*sizeof((*rdb->paf)));
    REALLOC(rdb->reverse_paf, rdb->total_reads);
    memset(rdb->reverse_paf+rid_n0, 0, (rdb->total_reads-rid_n0)*sizeof((*rdb->reverse_paf)));

    ruI->len = rdb->total_reads;
    REALLOC(ruI->index, ruI->len); 
    memset(ruI->index, -1, sizeof((*(ruI->index)))*(ruI->len));

    reset_bub_label_t(uopt->b_mask_t, ng, 0, 0);
    uopt->coverage_cut = (*cov);
    uopt->reverse_sources = rdb->reverse_paf;
    uopt->sources = rdb->paf;
    destory_UC_Read(&g_read);
    fprintf(stderr, "-[M::%s] rid_n0::%lu, rid_n::%lu\n", __func__, rid_n0, rid_n);
}


int64_t reload_uovl(all_ul_t *x, char* file_name)
{
    char* gfa_name = NULL; MALLOC(gfa_name, strlen(file_name)+100);
    sprintf(gfa_name, "%s.%s.ul.ovlp.bin", file_name, "re");
    FILE* fp = fopen(gfa_name, "r"); free(gfa_name);
    if (!fp) return 0;
    uint64_t k; size_t kn; ul_vec_t *p = NULL;

    ///skip test ug
    size_t tt, i; uint32_t t; ma_utg_t ua;
    fread(&tt, sizeof(tt), 1, fp);
    for (i = 0; i < tt; i++) {
        fread(&t, sizeof(t), 1, fp);
        fread(&t, sizeof(t), 1, fp);
        fread(&(ua.start), sizeof(ua.start), 1, fp);
        fread(&(ua.end), sizeof(ua.end), 1, fp);
        fread(&(ua.n), sizeof(ua.n), 1, fp);
        fseek(fp, sizeof(uint64_t)*ua.n, SEEK_CUR);
        // fread(ua.a, sizeof(uint64_t), ua.n, fp);
    }



    fseek(fp, sizeof(x->nid.n), SEEK_CUR);
	for (k = 0; k < x->nid.n; k++) {
        fseek(fp, sizeof(x->nid.a[k].n), SEEK_CUR); 
        fseek(fp, sizeof((*(x->nid.a[k].a)))*x->nid.a[k].n, SEEK_CUR);
	}

    // free(x->ridx.idx.a); x->ridx.idx.a = NULL;
    fread(&kn, sizeof(kn), 1, fp); 
    fseek(fp, sizeof((*(x->ridx.idx.a)))*kn, SEEK_CUR);

    // free(x->ridx.occ.a); x->ridx.occ.a = NULL;
    fread(&kn, sizeof(kn), 1, fp); 
    fseek(fp, sizeof((*(x->ridx.occ.a)))*kn, SEEK_CUR);

    fread(&x->n, sizeof(x->n), 1, fp); assert(x->nid.n == x->n);
    if(x->n > x->m) kv_resize(ul_vec_t, *x, x->n);    
    fprintf(stderr, "[M::%s] x->n::%u, x->nid.n::%u\n", __func__, (uint32_t)x->n, (uint32_t)x->nid.n);
    for (k = 0; k < x->n; k++) {
        p = &(x->a[k]);
        fread(&p->dd, sizeof(p->dd), 1, fp);
        fread(&p->rlen, sizeof(p->rlen), 1, fp);

        fread(&p->r_base.n, sizeof(p->r_base.n), 1, fp); 
        if(p->r_base.n > p->r_base.m) kv_resize(uint8_t, p->r_base, p->r_base.n);        
        fread(p->r_base.a, sizeof((*(p->r_base.a))), p->r_base.n, fp);

        fread(&p->bb.n, sizeof(p->bb.n), 1, fp); 
        if(p->bb.n > p->bb.m) kv_resize(uc_block_t, p->bb, p->bb.n);
        fread(p->bb.a, sizeof((*(p->bb.a))), p->bb.n, fp);

        fread(&p->N_site.n, sizeof(p->N_site.n), 1, fp); 
        if(p->N_site.n > p->N_site.m) kv_resize(uint32_t, p->N_site, p->N_site.n);        
        fread(p->N_site.a, sizeof((*(p->N_site.a))), p->N_site.n, fp);
    }
    fclose(fp);
    return 1;
}

static void gen_scaffold_id(void *data, long i, int tid) // callback for kt_for()
{
    scaf_mul_t *ss = (scaf_mul_t *)data;
    ma_ug_t *ug = ss->ug;
    usc_t *z = &(ss->a[i]);
    uint32_t v = z->raw_sid, w = z->raw_eid, sv, sw;
    uint64_t *a, a_n, k; int64_t ql, mmql, max_ql; 
    uc_block_t *p, *n, *m0, *m1;
    sv = (v&1?((ug->u.a[v>>1].a[0]>>32)^1):(ug->u.a[v>>1].a[ug->u.a[v>>1].n-1]>>32));
    sw = (w&1?((ug->u.a[w>>1].a[ug->u.a[w>>1].n-1]>>32)^1):(ug->u.a[w>>1].a[0]>>32));
    // if((((v>>1) == 15016) && ((w>>1) == 2778)) || (((v>>1) == 2778) && ((w>>1) == 15016))) {
    //     fprintf(stderr, "[M::%s::]\t%.*s(%c)(sv>>1::%u)(utg%.6ul(%c))\t%.*s(%c)(sw>>1::%u)(utg%.6ul(%c))\n", __func__, 
    //         (int)Get_NAME_LENGTH(R_INF, (sv>>1)), Get_NAME(R_INF, (sv>>1)), "+-"[sv&1], sv>>1, (v>>1)+1, "+-"[v&1],
    //         (int)Get_NAME_LENGTH(R_INF, (sw>>1)), Get_NAME(R_INF, (sw>>1)), "+-"[sw&1], sw>>1, (w>>1)+1, "+-"[w&1]);
    // }
    sv = Get_READ_LENGTH(R_INF, (sv>>1)); sw = Get_READ_LENGTH(R_INF, (sw>>1));
    max_ql = MIN(sv, sw); max_ql = -max_ql;
    // if((((v>>1) == 235) && ((w>>1) == 153)) || (((v>>1) == 153) && ((w>>1) == 235))) {
    //     fprintf(stderr, "[M::%s] sv::%u, sw::%u, max_ql::%ld, v>>1::%u(%c), w>>1::%u(%c)\n", 
    //     __func__, sv, sw, max_ql, v>>1, "+-"[v&1], w>>1, "+-"[w&1]);
    // }
    // if((((v>>1) == 235) && ((w>>1) == 153)) || (((v>>1) == 153) && ((w>>1) == 235))) {
    //     print_ul_alignment(ug, &UL_INF, 5625, "+++");
    // }
    uint32_t uv, uw, rev = 0, id = (uint32_t)-1;
    p = n = m0 = m1 = NULL;
    a = UL_INF.ridx.occ.a + UL_INF.ridx.idx.a[v>>1];
	a_n = UL_INF.ridx.idx.a[(v>>1)+1] - UL_INF.ridx.idx.a[v>>1]; 
	for (k = 0, mmql = INT32_MAX; k < a_n; k++) {
        p = &(UL_INF.a[a[k]>>32].bb.a[(uint32_t)(a[k])]); assert(p->hid == (v>>1));
        if(UL_INF.a[a[k]>>32].dd != 3) continue;///dd = 3 is saved
		if(p->base || (!p->el) || (!p->pchain)) continue;
        uv = (((uint32_t)(p->hid))<<1)|((uint32_t)(p->rev));

        // if((uv == v) && (p->aidx != (uint32_t)-1)) {
        if((uv == v) && (p->aidx == (uint32_t)-1) && (((uint32_t)(a[k]))+1 < UL_INF.a[a[k]>>32].bb.n)) {
			n = &(UL_INF.a[a[k]>>32].bb.a[((uint32_t)(a[k]))+1]);
            if(n->base || (!n->el) || (!n->pchain) || (n->pidx != (uint32_t)-1)) continue;
			uw = (((uint32_t)(n->hid))<<1)|((uint32_t)(n->rev));
			if(uw == w) {
                // if((((v>>1) == 15016) && ((w>>1) == 2778)) || (((v>>1) == 2778) && ((w>>1) == 15016))) {
                //     fprintf(stderr, "+[M::%s] ul_id::%lu, qs::%u, qe::%u\n", __func__, a[k]>>32, n->qs, n->qe);
                //     // print_ul_alignment(ug, &UL_INF, a[k]>>32, "+++");
                // }
                ql = ((int64_t)n->qs) - ((int64_t)p->qe);
                if(ql > max_ql && n->qe > p->qe && n->qs > p->qs) {
                    if(ql < mmql) {
                        mmql = ql; rev = 0; m0 = p; m1 = n; id = a[k]>>32;
                    }
                }
            }
		}

		// if(((uv^1) == v) && (p->pidx != (uint32_t)-1)) {
        if(((uv^1) == v) && (p->pidx == (uint32_t)-1) && (((uint32_t)(a[k])) > 0)) {
			n = &(UL_INF.a[a[k]>>32].bb.a[((uint32_t)(a[k]))-1]);
            if(n->base || (!n->el) || (!n->pchain) || (n->aidx != (uint32_t)-1)) continue;
			uw = (((uint32_t)(n->hid))<<1)|((uint32_t)(n->rev)); uw ^= 1;
			if(uw == w) {
                // if((((v>>1) == 15016) && ((w>>1) == 2778)) || (((v>>1) == 2778) && ((w>>1) == 15016))) {
                //     fprintf(stderr, "-[M::%s] ul_id::%lu, qs::%u, qe::%u\n", __func__, a[k]>>32, n->qs, n->qe);
                //     // print_ul_alignment(ug, &UL_INF, a[k]>>32, "+++");
                // }
                ql = ((int64_t)p->qs) - ((int64_t)n->qe);
                if(ql > max_ql && n->qe < p->qe && n->qs < p->qs) {
                    if(ql < mmql) {
                        mmql = ql; rev = 1; m0 = n; m1 = p; id = a[k]>>32;
                    }
                }
            }
		}
    }

    if(m0 && m1 && id != ((uint32_t)-1)) {
        z->ulid = (id<<1)|rev;
        z->raw_sof = m0->qe;
        z->raw_eof = m1->qs;
        // if((((v>>1) == 15016) && ((w>>1) == 2778)) || (((v>>1) == 2778) && ((w>>1) == 15016))) {
        //     fprintf(stderr, "[M::%s] id::%u, v>>1::%u(%c), w>>1::%u(%c), raw_sof::%u, raw_eof::%u\n", 
        //     __func__, id, v>>1, "+-"[v&1], w>>1, "+-"[w&1], z->raw_sof, z->raw_eof);
        // }
    }
}

void reload_uu(all_ul_t *x, ma_ug_t *raw_g, char* file_name, const char *bin_file)
{
    char* gfa_name = NULL; MALLOC(gfa_name, strlen(asm_opt.output_file_name)+50);
	sprintf(gfa_name, "%s.%s", asm_opt.output_file_name, bin_file);
    clear_all_ul_t(x);
    load_all_ul_t(x, gfa_name, &R_INF, raw_g);
    // filter_ul_ug(raw_g);//no filter since we would like to use all alignments
    gen_ul_vec_rid_t(x, NULL, raw_g);
    free(gfa_name);
}

void fill_scaffolds(usc_t *a, uint64_t a_n, ma_ug_t *raw_g, int64_t min_arc_len, const char *bin_file)
{   
    scaf_mul_t ss; uint32_t k, occ = 0;
    reload_uu(&UL_INF, raw_g, asm_opt.output_file_name, bin_file);
    // reload_uovl(&UL_INF, asm_opt.output_file_name);
    // filter_ul_ug(raw_g);
    // free(UL_INF.ridx.idx.a); UL_INF.ridx.idx.n = UL_INF.ridx.idx.m = 0;
    // free(UL_INF.ridx.occ.a); UL_INF.ridx.occ.n = UL_INF.ridx.occ.m = 0;
    // gen_ul_vec_rid_t(&UL_INF, NULL, raw_g);
    ss.a = a; ss.a_n = a_n; ss.ug = raw_g;
    kt_for(asm_opt.thread_num, gen_scaffold_id, &ss, a_n);
    for (k = 0; k < a_n; k++) {
        if(a[k].ulid == ((uint32_t)-1)) continue;
        // if(a[k].raw_eof + min_arc_len <= a[k].raw_sof) {///has an overlap longer than min_arc_len; 
        if(a[k].raw_eof <= a[k].raw_sof) {
            continue;///no need ul, could directly use hifi
        }
        // if(a[k].raw_eof <= a[k].raw_sof) continue;///no need UL reads
        // fprintf(stderr, "[M::%s] ulid::%u(%c), dd::%u\n", 
        // __func__, a[k].ulid>>1, "+-"[a[k].ulid&1], UL_INF.a[a[k].ulid>>1].dd);
        UL_INF.a[a[k].ulid>>1].dd = 4;///load
        occ++;
    }

    // fprintf(stderr, "[M::%s] occ::%u\n", __func__, occ);
    if(occ) load_scaf_base(&UL_INF, asm_opt.output_file_name, bin_file);
}

void update_paf(ma_hit_t_alloc *src, ma_hit_t_alloc *r_src, uint64_t *rmap, uint64_t rid_n, uint64_t pre_gn, asg_t *ng);
void push_scaff_node(ma_hit_t_alloc *src, uint64_t v, uint64_t w, uint64_t ol, asg_t *ng);

ma_ug_t *convert_usg_t(usg_t *ng, ma_ug_t *ref)
{
    ma_ug_t *ug = NULL; uint32_t k, nv, z; usg_arc_t *av; asg_arc_t *p;
    ma_utg_t *des, *src;
    CALLOC(ug, 1); ug->g = asg_init(); 
    ug->u.n = ug->u.m = ng->n; CALLOC(ug->u.a, ug->u.n);
    for (k = 0; k < ng->n; k++) {
        asg_seq_set(ug->g, k, ng->a[k].len, ng->a[k].del);
        ug->g->seq[k].c = ref->g->seq[ng->a[k].mm].c;

        av = usg_arc_a(ng, (k<<1)); nv = usg_arc_n(ng, (k<<1));
        for (z = 0; z < nv; z++) {
            if(av[z].del) continue;
            p = asg_arc_pushp(ug->g); memset(p, 0, sizeof((*p)));
            p->ul = av[z].ul; p->v = av[z].v; p->ol = av[z].ol; p->del = av[z].del;
        }

        av = usg_arc_a(ng, (k<<1)+1); nv = usg_arc_n(ng, (k<<1)+1);
        for (z = 0; z < nv; z++) {
            if(av[z].del) continue;
            p = asg_arc_pushp(ug->g); memset(p, 0, sizeof((*p)));
            p->ul = av[z].ul; p->v = av[z].v; p->ol = av[z].ol; p->del = av[z].del;
        }

        des = &(ug->u.a[k]); src = &(ref->u.a[ng->a[k].mm]);
        (*des) = (*src); des->a = NULL; des->s = NULL;
        if(src->a) {
            MALLOC(des->a, des->m); 
            memcpy(des->a, src->a, src->n*sizeof((*(src->a))));
        }

        if(src->s) {
            MALLOC(des->s, des->len);
            memcpy(des->s, src->s, src->len*sizeof((*(src->s))));
        }
    }

    asg_cleanup(ug->g);
    return ug;
}

int32_t gen_spec_rc_edge(asg_t *rg, ug_opt_t *uopt, uint32_t v, uint32_t w, asg_arc_t *t, ma_hit_t **te)
{
    uint32_t k, qn, tn; int32_t r; ma_hit_t_alloc *s = &(uopt->reverse_sources[v>>1]); asg_arc_t p;
    for (k = 0; k < s->length; k++) {
        qn = Get_qn(s->buffer[k]); tn = Get_tn(s->buffer[k]);
        if(tn != (w>>1)) continue;
        r = ma_hit2arc(&(s->buffer[k]), rg->seq[qn].len, rg->seq[tn].len, uopt->max_hang, asm_opt.max_hang_rate, uopt->min_ovlp, &p);
        if(r < 0) continue;
        if((p.ul>>32) != v || p.v != w) continue;
        *t = p; t->ou = 0; if(te) (*te) = &(s->buffer[k]);
        return 1;
    }
    return -1;
}

void push_direct_scaff_node(usc_t *psa, asg_t *ng, ug_opt_t *uopt, uint64_t vx, uint64_t wx)
{
    uint64_t uol, dif; int32_t r; asg_arc_t t, *p; ma_hit_t *rc0, *rc1;
    uol = psa->raw_sof - psa->raw_eof;
    assert(uol >= (uint64_t)uopt->min_ovlp);
    r = gen_spec_rc_edge(ng, uopt, vx, wx, &t, &rc0);
    if(r >= 0) {
        if(uol > t.ol) dif = uol - t.ol;
        else dif = t.ol - uol;
        if((dif > (uol*0.05)) || (dif > (t.ol*0.05))) {
            if(dif > 16) r = -1;
        }
    }
    if(r >= 0) {
        r = gen_spec_rc_edge(ng, uopt, wx^1, vx^1, &t, &rc1);
        if(uol > t.ol) dif = uol - t.ol;
        else dif = t.ol - uol;
        if((dif > (uol*0.05)) || (dif > (t.ol*0.05))) {
            if(dif > 16) r = -1;
        }
    }

    if(r >= 0) {
        r = gen_spec_rc_edge(ng, uopt, vx, wx, &t, &rc0); assert(r >= 0);
        p = asg_arc_pushp(ng); *p = t;

        r = gen_spec_rc_edge(ng, uopt, wx^1, vx^1, &t, &rc1); assert(r >= 0);
        p = asg_arc_pushp(ng); *p = t;

        add_ma_hit_t_alloc(&(R_INF.paf[Get_qn((*rc0))]), rc0);
        add_ma_hit_t_alloc(&(R_INF.paf[Get_qn((*rc1))]), rc1);
    } else {
        push_scaff_node(R_INF.paf, vx, wx, psa->raw_sof-psa->raw_eof, ng);

        r = gen_spec_edge(ng, uopt, vx, wx, &t);
        assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
        r = gen_spec_edge(ng, uopt, wx^1, vx^1, &t);
        assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
    }
}

void print_debug_scaffold_nodes(ug_opt_t *uopt, usc_t *a, uint32_t a_n, asg_t *ng)
{
    uint32_t i, k, z, qn, tn, rl, hrl; usc_t *psa = NULL; ma_hit_t_alloc *s = NULL;
    UC_Read r0, r1; init_UC_Read(&r0); init_UC_Read(&r1);
    fprintf(stderr, "[M::%s] a_n::%u\n", __func__, a_n);
    for (i = 0; i < a_n; i++) {
        psa = &(a[i]);
        s = &(uopt->sources[psa->nid]);
        fprintf(stderr, "[M::%s::%u]%.*s(id::%u)\tlen::%lu\tulid::%u\trev::%u\tdel::%u\n", 
        __func__, i, (int32_t)Get_NAME_LENGTH(R_INF, psa->nid), 
        Get_NAME(R_INF, psa->nid), psa->nid, Get_READ_LENGTH(R_INF, psa->nid), psa->ulid>>1, psa->ulid&1, 
        ng->seq[psa->nid].del);

        recover_UC_Read(&r0, &R_INF, psa->nid);
        fprintf(stderr, "%.*s\n", (int32_t)Get_READ_LENGTH(R_INF, psa->nid), r0.seq);

        for (k = 0; k < s->length; k++) {
            qn = Get_qn(s->buffer[k]); tn = Get_tn(s->buffer[k]); assert(qn == psa->nid);
            fprintf(stderr, "q::[%u, %u)\t%c\tt::[%u, %u)\t%.*s(id::%u)\tlen::%lu\n", 
            Get_qs(s->buffer[k]), Get_qe(s->buffer[k]), "+-"[s->buffer[k].rev], 
            Get_ts(s->buffer[k]), Get_te(s->buffer[k]),
            (int32_t)Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn), tn, Get_READ_LENGTH(R_INF, tn));
            if(psa->raw_eof >= psa->raw_sof) {
                resize_UC_Read(&r0, Get_qe(s->buffer[k])-Get_qs(s->buffer[k]));
                recover_UC_Read_sub_region(r0.seq, Get_qs(s->buffer[k]), Get_qe(s->buffer[k])-Get_qs(s->buffer[k]), 0, &R_INF, qn);
                
                resize_UC_Read(&r1, Get_te(s->buffer[k])-Get_ts(s->buffer[k]));
                recover_UC_Read_sub_region(r1.seq, Get_ts(s->buffer[k]), Get_te(s->buffer[k])-Get_ts(s->buffer[k]), 0, &R_INF, tn);
                
                assert(Get_te(s->buffer[k])-Get_ts(s->buffer[k]) == Get_qe(s->buffer[k])-Get_qs(s->buffer[k]));
                rl = Get_te(s->buffer[k])-Get_ts(s->buffer[k]);
                if(s->buffer[k].rev) {
                    char t; hrl = rl>>1;
                    for (z = 0; z < hrl; z++) {
                        t = r1.seq[rl-z-1];
                        r1.seq[rl-z-1] = RC_CHAR(r1.seq[z]);
                        r1.seq[z] = RC_CHAR(t);
                    }
                    if(rl&1) r1.seq[z] = RC_CHAR(r1.seq[z]);
                }
                assert(memcmp(r0.seq, r1.seq, rl) == 0);
            }
        }
    }
    destory_UC_Read(&r0); destory_UC_Read(&r1);
}

asg_t *renew_ng(usg_t *eg, ma_ug_t *rug, asg_t *sg, ug_opt_t *uopt, ma_sub_t **cov, R_to_U *ruI, uint64_t scaffold_len, const char *bin_file)
{
    init_aux_table();
    ma_utg_t *u; uint64_t i, v, w, m, h, z, raw_v, raw_w, nocc, nv, vx, wx; int32_t r;
    asg_arc_t t, *p; usg_arc_t *av = NULL; usc_t *psa; ma_ug_t *ug1 = NULL;
    asg_ext_t ext; memset(&ext, 0, sizeof(ext)); ext.ext = asg_init(); asg_t *ng = ext.ext;
    usc_vec_t sa; memset(&sa, 0, sizeof(sa)); 
    ext.a_n = sg->n_seq; 
    ext.cnt.n = ext.cnt.m = ext.a_n;
    CALLOC(ext.cnt.a, ext.cnt.n);///count

    ext.idx_a.n = ext.idx_a.m = sg->n_seq; 
    MALLOC(ext.idx_a.a, ext.idx_a.n); 
    memset(ext.idx_a.a, -1, sizeof(*(ext.idx_a.a))*ext.idx_a.n);//map
    ug1 = convert_usg_t(eg, rug);

    // fprintf(stderr, "-0-[M::%s]\n", __func__);
    // print_vw_edge(sg, 681356, 73351, __func__);
    // prt_spec_edge(sg, R_INF.paf, R_INF.total_reads, 681356, 73351, uopt, "src");
    // // prt_spec_edge(sg, R_INF.reverse_paf, R_INF.total_reads, 681356, 73351, uopt, "rsrc");

    for (i = 0; i < ug1->u.n; ++i) {
        if(ug1->g->seq[i].del) continue;
        u = &(ug1->u.a[i]);        
        for (m = 0; m < u->n; m++) {
            v = u->a[m]>>32;
            ext.cnt.a[v>>1]++;
            if(ext.cnt.a[v>>1] == 1) {
                h = v>>1; ext.idx_a.a[h] = v>>1;
            } else {
                h = ext.idx_a.n; 
                kv_push(uint64_t, ext.idx_a, (v>>1));
            }

            z = (h<<1)|(v&1); z <<= 32; z |= ((uint32_t)u->a[m]);
            u->a[m] = z;
        }
    }

    for (i = 0; i < sg->n_seq; ++i) {
        asg_seq_set(ng, i, sg->seq[i].len, ext.idx_a.a[i]==((uint64_t)-1)?1:0);
        ng->seq[i].c = 0;
        if(ext.idx_a.a[i]!=((uint64_t)-1)) assert(ext.idx_a.a[i] == i);
        ext.idx_a.a[i] = i;///set for the deleted read
    }
    for (; i < ext.idx_a.n; i++) {
        asg_seq_set(ng, i, sg->seq[ext.idx_a.a[i]].len, 0);
        ng->seq[i].c = 0; assert(!(ng->seq[ext.idx_a.a[i]].del));
    }
    ng->r_seq = ng->n_seq;
    assert(ng->n_seq == ext.idx_a.n);


    for (i = 0, nocc = ng->r_seq, sa.n = 0; i < eg->n; ++i) {
        if(ug1->g->seq[i].del) continue;
        ///there shouldn't any scaffolding within the nodes
        v = i<<1; nv = usg_arc_n(eg, v); av = usg_arc_a(eg, v);
        for (m = 0; m < nv; m++) {
            if(av[m].del || eg->a[av[m].v>>1].del) continue;
            w = av[m].v;
            vx = (v&1?((ug1->u.a[v>>1].a[0]>>32)^1):(ug1->u.a[v>>1].a[ug1->u.a[v>>1].n-1]>>32));
            wx = (w&1?((ug1->u.a[w>>1].a[ug1->u.a[w>>1].n-1]>>32)^1):(ug1->u.a[w>>1].a[0]>>32));
            raw_v = (ext.idx_a.a[vx>>1]<<1)|(vx&1);
            raw_w = (ext.idx_a.a[wx>>1]<<1)|(wx&1);
            assert((raw_v>>1) < sg->n_seq);
            assert((raw_w>>1) < sg->n_seq);
            if(vx > wx) continue;
            if((vx == wx) && (vx&1)) continue;
            if(gen_spec_edge(sg, uopt, raw_v, raw_w, &t) < 0) {
                asg_seq_set(ng, nocc, 0, 0);///this is a scaffold node
                ng->seq[nocc].c = 0;
                kv_push(uint64_t, ext.idx_a, ((uint64_t)-1));

                kv_pushp(usc_t, sa, &psa); 
                psa->ulid = (uint32_t)-1; psa->nid = nocc;
                psa->raw_sid = (eg->a[v>>1].mm<<1)|(v&1);///raw unitig id
                psa->raw_eid = (eg->a[w>>1].mm<<1)|(w&1);///raw unitig id
                psa->raw_sof = psa->raw_eof = (uint32_t)-1;
                nocc++;///a scaffold node
            }
        }

        v = (i<<1)+1; nv = usg_arc_n(eg, v); av = usg_arc_a(eg, v);
        for (m = 0; m < nv; m++) {
            if(av[m].del || eg->a[av[m].v>>1].del) continue;
            w = av[m].v;
            vx = (v&1?((ug1->u.a[v>>1].a[0]>>32)^1):(ug1->u.a[v>>1].a[ug1->u.a[v>>1].n-1]>>32));
            wx = (w&1?((ug1->u.a[w>>1].a[ug1->u.a[w>>1].n-1]>>32)^1):(ug1->u.a[w>>1].a[0]>>32));
            raw_v = (ext.idx_a.a[vx>>1]<<1)|(vx&1);
            raw_w = (ext.idx_a.a[wx>>1]<<1)|(wx&1);
            assert((raw_v>>1) < sg->n_seq);
            assert((raw_w>>1) < sg->n_seq);
            if(vx > wx) continue;
            if((vx == wx) && (vx&1)) continue;
            if(gen_spec_edge(sg, uopt, raw_v, raw_w, &t) < 0) {
                asg_seq_set(ng, nocc, 0, 0);///this is a scaffold node
                ng->seq[nocc].c = 0;
                kv_push(uint64_t, ext.idx_a, ((uint64_t)-1));

                kv_pushp(usc_t, sa, &psa); 
                psa->ulid = (uint32_t)-1; psa->nid = nocc;
                psa->raw_sid = (eg->a[v>>1].mm<<1)|(v&1);///raw unitig id
                psa->raw_eid = (eg->a[w>>1].mm<<1)|(w&1);///raw unitig id
                psa->raw_sof = psa->raw_eof = (uint32_t)-1;
                nocc++;///a scaffold node
            }
        }
    }

    assert(ng->n_seq == ext.idx_a.n); 
    CALLOC(ng->seq_vis, (ng->n_seq<<1));
    ///# scaffolding nodes
    if(sa.n > 0) fill_scaffolds(sa.a, sa.n, rug, uopt->min_ovlp, bin_file);
    realloc_rdb_adv(&(R_INF), &UL_INF, cov, ruI, ext.idx_a.a, ext.idx_a.n, scaffold_len, (char *)"scaf", ng, uopt, rug, sa.a);     
    update_paf(R_INF.paf, R_INF.reverse_paf, ext.idx_a.a, ext.idx_a.n, sg->n_seq, ng);

    for (i = 0, nocc = ng->r_seq; i < eg->n; ++i) {
        if(ug1->g->seq[i].del) continue;
        u = &(ug1->u.a[i]);
        for (m = 1; m < u->n; m++) {
            // fprintf(stderr, "[M::%s] i::%lu, m::%lu\n", __func__, i, m);
            v = u->a[m-1]>>32; w = u->a[m]>>32;
            r = gen_spec_edge(ng, uopt, v, w, &t);
            ///there shouldn't any scaffolding within the nodes
            assert(r >= 0);
            p = asg_arc_pushp(ng); *p = t;
            r = gen_spec_edge(ng, uopt, w^1, v^1, &t);
            assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
        }

        v = i<<1; nv = usg_arc_n(eg, v); av = usg_arc_a(eg, v);
        for (m = 0; m < nv; m++) {
            if(av[m].del || eg->a[av[m].v>>1].del) continue;
            w = av[m].v;
            vx = (v&1?((ug1->u.a[v>>1].a[0]>>32)^1):(ug1->u.a[v>>1].a[ug1->u.a[v>>1].n-1]>>32));
            wx = (w&1?((ug1->u.a[w>>1].a[ug1->u.a[w>>1].n-1]>>32)^1):(ug1->u.a[w>>1].a[0]>>32));
            if(vx > wx) continue;
            if((vx == wx) && (vx&1)) continue;///it is possible

            r = gen_spec_edge(ng, uopt, vx, wx, &t);
            if(r >= 0) {
                p = asg_arc_pushp(ng); *p = t;
                r = gen_spec_edge(ng, uopt, wx^1, vx^1, &t);
                // if(!(r >= 0)) {
                //     fprintf(stderr, "[M::%s]\tvx>>1::%lu(%c)\twx>>1::%lu(%c)\n", 
                //     __func__, vx>>1, "+-"[vx&1], wx>>1, "+-"[wx&1]);
                //     fprintf(stderr, "[M::%s]\t%.*s(%c)->%.*s(%c)\n", 
                //     __func__, (int)Get_NAME_LENGTH(R_INF, (vx>>1)), Get_NAME(R_INF, (vx>>1)), "+-"[vx&1],
                //     (int)Get_NAME_LENGTH(R_INF, (wx>>1)), Get_NAME(R_INF, (wx>>1)), "+-"[wx&1]);
                // }
                assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
            } else {
                psa = &(sa.a[nocc-ng->r_seq]); z = nocc<<1; 
                assert(psa->nid == nocc);
                assert(psa->raw_sid == ((eg->a[v>>1].mm<<1)|(v&1)));
                assert(psa->raw_eid == ((eg->a[w>>1].mm<<1)|(w&1)));
                if(ng->seq[nocc].del) {///direct link vx and wx
                    push_direct_scaff_node(psa, ng, uopt, vx, wx);
                } else {
                    push_scaff_node(R_INF.paf, vx, z, uopt->min_ovlp, ng);
                    push_scaff_node(R_INF.paf, z, wx, uopt->min_ovlp, ng);

                    r = gen_spec_edge(ng, uopt, vx, z, &t);
                    assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
                    r = gen_spec_edge(ng, uopt, z^1, vx^1, &t);
                    assert(r >= 0); p = asg_arc_pushp(ng); *p = t;

                    r = gen_spec_edge(ng, uopt, z, wx, &t);
                    assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
                    r = gen_spec_edge(ng, uopt, wx^1, z^1, &t);
                    assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
                }
                nocc++;
            }
        }

        v = (i<<1)+1; nv = usg_arc_n(eg, v); av = usg_arc_a(eg, v);
        for (m = 0; m < nv; m++) {
            if(av[m].del || eg->a[av[m].v>>1].del) continue;
            w = av[m].v;
            vx = (v&1?((ug1->u.a[v>>1].a[0]>>32)^1):(ug1->u.a[v>>1].a[ug1->u.a[v>>1].n-1]>>32));
            wx = (w&1?((ug1->u.a[w>>1].a[ug1->u.a[w>>1].n-1]>>32)^1):(ug1->u.a[w>>1].a[0]>>32));
            if(vx > wx) continue;
            if((vx == wx) && (vx&1)) continue;///it is possible

            r = gen_spec_edge(ng, uopt, vx, wx, &t);
            if(r >= 0) {
                p = asg_arc_pushp(ng); *p = t;
                r = gen_spec_edge(ng, uopt, wx^1, vx^1, &t);
                assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
            } else {
                psa = &(sa.a[nocc-ng->r_seq]); z = nocc<<1; 
                assert(psa->nid == nocc);
                assert(psa->raw_sid == ((eg->a[v>>1].mm<<1)|(v&1)));
                assert(psa->raw_eid == ((eg->a[w>>1].mm<<1)|(w&1)));
                if(ng->seq[nocc].del) {///direct link vx and wx
                    push_direct_scaff_node(psa, ng, uopt, vx, wx);
                } else {
                    push_scaff_node(R_INF.paf, vx, z, uopt->min_ovlp, ng);
                    push_scaff_node(R_INF.paf, z, wx, uopt->min_ovlp, ng);

                    r = gen_spec_edge(ng, uopt, vx, z, &t);
                    assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
                    r = gen_spec_edge(ng, uopt, z^1, vx^1, &t);
                    assert(r >= 0); p = asg_arc_pushp(ng); *p = t;

                    r = gen_spec_edge(ng, uopt, z, wx, &t);
                    assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
                    r = gen_spec_edge(ng, uopt, wx^1, z^1, &t);
                    assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
                }
                nocc++;
            }
        }
    }

    asg_cleanup(ng); ng->r_seq = ng->n_seq;
    free(ext.cnt.a); free(ext.idx_a.a); free(sa.a); ma_ug_destroy(ug1);

    fprintf(stderr, "[M::%s] nocc::%lu, ng->n_seq::%u, sg->n_seq::%u\n", 
    __func__, nocc, (uint32_t)ng->n_seq, (uint32_t)sg->n_seq);
    assert(nocc == ng->n_seq);
    return ng;
}

void u2g_hybrid_detan_iter(ul_resolve_t *uidx, usg_t *ng, uint32_t max_ext, uint32_t clean_round, asg64_v *in, asg64_v *ib)
{
    uint32_t n_vtx = ng->n<<1, ncut = 0, k;
    uint8_t *ff; CALLOC(ff, n_vtx);
    uint32_t *ng_occ; CALLOC(ng_occ, n_vtx); 
    uint64_t *i_idx; CALLOC(i_idx, n_vtx);
    asg64_v b64, ub64; kv_init(b64); kv_init(ub64); 
    asg64_v tx = {0,0,0}, tb = {0,0,0}, *ob = NULL, *ub = NULL;
    ob = (in?(in):(&tx)); ub = (ib?(ib):(&tb)); ob->n = ub->n = 0;
    // fprintf(stderr, "\n[M::%s::] asm_opt.is_low_het_ul::%u, max_ext::%u\n", 
    //     __func__, asm_opt.is_low_het_ul, max_ext);
    // prt_usg_t(uidx, ng, "ng0");
    for (k = 0; k < clean_round; k++) {
        ncut += ug_ext_strict(ob, ub, uidx, ng, max_ext, &ff, &ng_occ, &i_idx, &b64, &ub64);
        // if(asm_opt.is_low_het_ul) break;
        ncut += ug_ext_free(ob, ub, uidx, ng, max_ext, &ff, &ng_occ, &i_idx, &b64, &ub64, 48, 1);
        // prt_usg_t(uidx, ng, "ng_python");
        ncut += ug_ext_free(ob, ub, uidx, ng, max_ext, &ff, &ng_occ, &i_idx, &b64, &ub64, 16, 1);

        ncut += ug_ext_free(ob, ub, uidx, ng, max_ext, &ff, &ng_occ, &i_idx, &b64, &ub64, 0, 0);
        ///renew bubble for ug_ext_strict
        n_vtx = ng->n<<1; 
        REALLOC(ff, n_vtx); REALLOC(ng_occ, n_vtx); REALLOC(i_idx, n_vtx);
        renew_usg_t_bub(uidx, ng, ng_occ, ff, 0);
    }
    // // ug_ext_strict(ob, ub, uidx, ng, max_ext, &ff, &ng_occ, &i_idx, &b64, &ub64);
    // ncut += ug_ext_free(ob, ub, uidx, ng, max_ext, &ff, &ng_occ, &i_idx, &b64, &ub64, 50);
    // prt_usg_t(uidx, ng, "ng.db");
    // ncut += ug_ext_free(ob, ub, uidx, ng, max_ext, &ff, &ng_occ, &i_idx, &b64, &ub64, 10);
    // ncut += ug_ext_strict(ob, ub, uidx, ng, max_ext, &ff, &ng_occ, &i_idx, &b64, &ub64);
    // prt_usg_t(uidx, ng, "ng.db");
    if(ncut) {
        usg_cleanup(ng); 
        // if((max_ext>>1) > 0) usg_arc_cut_tips(ng, (max_ext>>1), 1, ub);
        usg_arc_cut_tips(ng, max_ext, 1, ub);

        // if(ng->n > 139131) {
        //     fprintf(stderr, ">[M::%s::] ng->a[139131].del::%u\n", __func__, 
        //     ng->a[139131].del);
        // }
    }
    // prt_usg_t(uidx, ng, "ng2");
    // u2g_hybrid_aln(uidx, ng, ob, ub);
    // if(ug_ext(uidx, ng, ub->a, ub->n, ob->a, max_ext, ff, ng_occ, i_idx, &b64, &ub64)) {
    //     // debug_sysm_usg_t(ng, __func__);
    //     usg_cleanup(ng);
    //     // debug_sysm_usg_t(ng, __func__);
    // }
    
    // prt_usg_t(uidx, ng, "ng0");
    // if(ug_ext_free(ob, ub, uidx, ng, max_ext, &ff, &ng_occ, &i_idx, &b64, &ub64, 50)) {
    //     usg_cleanup(ng);
    // }
    // prt_usg_t(uidx, ng, "ng1");
    
    if(!in) free(tx.a); if(!ib) free(tb.a);
    free(ng_occ); free(i_idx); free(ff); kv_destroy(b64); kv_destroy(ub64); 
}


/**

void u2g_hybrid_detan(ul_resolve_t *uidx, usg_t *ng, uint32_t max_ext, asg64_v *in, asg64_v *ib)
{
    uint64_t k, i, x, m, *tmp, sn; asg64_v tx = {0,0,0}, tb = {0,0,0}, *ob = NULL, *ub = NULL; 
    ob = (in?(in):(&tx)); ub = (ib?(ib):(&tb)); ob->n = ub->n = 0;

    for (k = 0; k < uidx->str_b.n_thread; k++) {
        uidx->str_b.buf[k].res_dump.n = uidx->str_b.buf[k].u.n = uidx->str_b.buf[k].o.n = 0;
    }
    kt_for(uidx->str_b.n_thread, worker_integer_realign_g, uidx, uidx->uovl.i_ug->u.n);

    for (k = ob->n = ub->n = m = 0; k < uidx->str_b.n_thread; k++) {
        for (i = 0; i < uidx->str_b.buf[k].res_dump.n; i++) {
            x = uidx->str_b.buf[k].res_dump.a[i];
            kv_push(uint64_t, *ob, x);//aln details
            if(x&((uint64_t)0x8000000000000000)) {
                x -= ((uint64_t)0x8000000000000000); x >>= 32; x <<= 32;//seq_id
                x |= ob->n;//offset
                kv_push(uint64_t, *ub, x);///idx:: seq_id|offset_in_ob
            } else {
                m++;
            }
        }
    }

    radix_sort_srt64(ub->a, ub->a + ub->n);
    kv_resize(uint64_t, *ub, ub->n+m); tmp = ub->a + ub->n; m = 0;
    for (k = 0; k < ub->n; k++) {
        sn = ((uint32_t)(ob->a[((uint32_t)ub->a[k])-1]));
        memcpy(tmp + m, ob->a + ((uint32_t)ub->a[k]), sn*sizeof((*tmp))); 
        ub->a[k] = m; ub->a[k] <<= 32; ub->a[k] |= sn;//offset_in_ob|occ
        m += sn;
    }
    assert(m <= ob->n);
    memcpy(ob->a, tmp, m*sizeof((*tmp))); ob->n = m;

    debug_sysm_usg_t(ng, __func__);
    // prt_usg_t(uidx, ng, "ng4");
    if(gen_unique_g_adv(uidx, ng, ub->a, ub->n, ob->a, max_ext)) {
        // usg_arc_t *z = get_usg_arc(ng, 2, 576), *q = get_usg_arc(ng, 577, 3);
        // fprintf(stderr, "xxxx0xxx[M::%s::] p->del::%u, q->del::%u\n", 
        //             __func__, z?z->del:1, q?q->del:1);
        ///debug
        debug_sysm_usg_t(ng, __func__);
        // z = get_usg_arc(ng, 2, 576); q = get_usg_arc(ng, 577, 3);
        // fprintf(stderr, "xxxx1xxx[M::%s::] p->del::%u, q->del::%u\n", 
        //             __func__, z?z->del:1, q?q->del:1);
        // fprintf(stderr, "+[M::%s::] ng->n::%u\n", __func__, (uint32_t)ng->n);
        usg_cleanup(ng);
        ///debug
        debug_sysm_usg_t(ng, __func__);
        // fprintf(stderr, "-[M::%s::] ng->n::%u\n", __func__, (uint32_t)ng->n);
    }
    // prt_usg_t(uidx, ng, "ng_dbg");
    if(!in) free(tx.a); if(!ib) free(tb.a);
}
**/

ma_ug_t *gen_debug_hybrid_ug(ul_resolve_t *uidx, usg_t *ng)
{
    ma_ug_t *ug = NULL; uint32_t k, nv, z; usg_arc_t *av; asg_arc_t *p;
    CALLOC(ug, 1); ug->g = asg_init(); 
    ug->u.n = ug->u.m = ng->n; CALLOC(ug->u.a, ug->u.n);
    for (k = 0; k < ng->n; k++) {
        asg_seq_set(ug->g, k, ng->a[k].len, ng->a[k].del);
        ug->g->seq[k].c = 0;

        av = usg_arc_a(ng, (k<<1)); nv = usg_arc_n(ng, (k<<1));
        for (z = 0; z < nv; z++) {
            if(av[z].del) continue;
            p = asg_arc_pushp(ug->g); memset(p, 0, sizeof((*p)));
            p->ul = av[z].ul; p->v = av[z].v; p->ol = av[z].ol; p->del = av[z].del;
        }

        av = usg_arc_a(ng, (k<<1)+1); nv = usg_arc_n(ng, (k<<1)+1);
        for (z = 0; z < nv; z++) {
            if(av[z].del) continue;
            p = asg_arc_pushp(ug->g); memset(p, 0, sizeof((*p)));
            p->ul = av[z].ul; p->v = av[z].v; p->ol = av[z].ol; p->del = av[z].del;
        }

        ug->u.a[k].len = ug->g->seq[k].len; 
        ug->u.a[k].n = ug->u.a[k].m = 1; CALLOC(ug->u.a[k].a, 1); 
        ug->u.a[k].a[0] = (((uint64_t)k)<<33)|((uint64_t)(ug->u.a[k].len));
    }
    asg_cleanup(ug->g);
    

    uint32_t i; ma_utg_t *u; kvec_asg_arc_t_warp e; kv_init(e.a); e.i = 0;
    for (i = 0; i < ug->u.n; i++) {
        ug->g->seq[i].c = PRIMARY_LABLE;
        u = &(ug->u.a[i]);
        if(u->m == 0) continue;
        merge_hybrid_utg_content(u, uidx->l1_ug, uidx->sg, ng, &e);
        ug->g->seq[i].len = u->len;
    }
    kv_destroy(e.a);
    return ug;
}


void prt_usg_t(ul_resolve_t *uidx, usg_t *ng, const char *cmd)
{
    ma_ug_t *ug = gen_debug_hybrid_ug(uidx, ng);
    print_debug_gfa(uidx->sg, ug, uidx->uopt->coverage_cut, cmd, 
    uidx->uopt->sources, uidx->uopt->ruIndex, uidx->uopt->max_hang, uidx->uopt->min_ovlp, 0, 0, 0);
    ma_ug_destroy(ug);
    // exit(1);
}

inline uint32_t usg_arc_occ(usg_t *g, uint32_t v, uint32_t *res)
{
    uint32_t i, kv, nv = usg_arc_n(g, v);; 
    usg_arc_t *av = usg_arc_a(g, v);
    for (i = kv = 0; i < nv; i++) {
        if(av[i].del) continue;
        if(res) res[kv] = av[i].v;
        kv++;
    }
    return kv;
}

inline uint32_t gen_usg_tig(usg_t *g, uint32_t sid, uint32_t *eid, int64_t *baseLen, buf_t* b)
{
    uint32_t v = sid, w, k, kv, nv, return_flag;
    usg_arc_t *av;
    (*baseLen) = 0; (*eid) = (uint32_t)-1;
    while (1) {
        kv = usg_arc_occ(g, v, NULL);
        (*eid) = v;
        if(b) kv_push(uint32_t, b->b, v);
        ///means reach the end of a unitig
        if(kv!=1) (*baseLen) += g->a[v>>1].len;
        if(kv==0) {
            return_flag = END_TIPS;
            break;
        } else if(kv>1) {
            return_flag = MUL_OUTPUT;
            break;
        }

        ///kv must be 1 here
        kv = usg_arc_occ(g, v, &w);
        if(usg_arc_occ(g, w^1, NULL)!=1) {
            (*baseLen) += g->a[v>>1].len;
            return_flag = MUL_INPUT;
            break;
        } else {
            av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
            for (k = 0; k < nv; k++) {
                if(av[k].del) continue;
                ///here is just one undeleted edge
                (*baseLen) += asg_arc_len(av[k]);
                break;
            }
        }

        v = w;
        if(v == sid) {
            return_flag = LOOP;
            break;
        } 
    }
    return return_flag;
}

uint64_t dfs_usg_t_dis(usg_t *g, buf_t *b, uint32_t x, uint32_t *p_bub)
{
    uint64_t len = 0; usg_arc_t *av = NULL; int64_t baseLen, uLen;
    uint32_t c_v, e_v, nv, convex, v, i, kv_0, kv_1, flag_0 = 0, flag_1 = 0, op;
    
    (*p_bub) = 0;
    if(b->a[x>>1].s || g->a[x>>1].del) return 0;
    b->S.n = 0;
    kv_push(uint32_t, b->S, x);

    while (b->S.n > 0) {
        b->S.n--;
        c_v = b->S.a[b->S.n];
        if(b->a[c_v>>1].s) continue;

        b->b.n = 0; //uint32_t gen_usg_tig(usg_t *g, uint32_t sid, uint32_t *eid, int64_t *baseLen, buf_t* b)
        op = gen_usg_tig(g, c_v, &convex, &baseLen, b);
        uLen = baseLen;
        for(i = 0; i < b->b.n; i++) b->a[b->b.a[i]>>1].s = 1;
        if(op == LOOP) return 0;


        e_v = convex^1;
        b->b.n = 0;
        op = gen_usg_tig(g, e_v, &convex, &baseLen, b);
        uLen = MAX(uLen, baseLen);

        len += uLen;

        
        v = c_v^1;
        nv = usg_arc_n(g, v); av = usg_arc_a(g, v);
        for (i = kv_0 = 0; i < nv; i++) {
            if(av[i].del) continue;
            kv_0++;
            if(b->a[av[i].v>>1].s) continue;
            kv_push(uint32_t, b->S, av[i].v);
        }

        v = e_v^1;
        nv = usg_arc_n(g, v); av = usg_arc_a(g, v);
        for (i = kv_1 = 0; i < nv; i++) {
            if(av[i].del) continue;
            kv_1++;
            if(b->a[av[i].v>>1].s) continue;
            kv_push(uint32_t, b->S, av[i].v);
        }

        if(kv_0 > 0 && kv_1 > 0) flag_0++;
        if(kv_0 > 1) flag_1++;
        if(kv_1 > 1) flag_1++;
    }

    if(flag_0 > 0 && flag_1 > 1) (*p_bub) = 1;
    return len;
}

uint64_t usg_bub_dis(usg_t *g, buf_t *b)
{
    usg_arc_t *av = NULL; uint64_t cLen = 0, mLen = 0;
    uint32_t n_vtx = g->n<<1, k, v, w, kv, nv, p_bub;

    for (v = 0; v < n_vtx; ++v) {
        if(b->a[v>>1].s || g->a[v>>1].del) continue;
        av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
        for (k = kv = 0; k < nv && kv <= 1; k++) {
            if(av[k].del) continue;
            w = av[k].v^1; kv++;
        }
        if(kv == 1) {
            av = usg_arc_a(g, w); nv = usg_arc_n(g, w);
            for (k = kv = 0; k < nv && kv <= 1; k++) {
                if(av[k].del) continue;
                w = av[k].v^1; kv++;
            }
            if(kv == 1) continue;
        }

        cLen = dfs_usg_t_dis(g, b, v^1, &p_bub);
        if(p_bub == 0) continue;///no bubble
        if(cLen > mLen) mLen = cLen;
    }

    for (k = 0; k < g->n; ++k) b->a[k].s = 0;
    b->S.n = b->b.n = 0;
    return mLen;
}

uint64_t usg_bub_identify(usg_t *g, uint32_t v0, uint64_t max_dist, uint64_t max_occ, buf_t *b)
{
    uint32_t i, n_pending = 0, is_first = 1, n_tips, tip_end;
	uint64_t n_pop = 0; uint32_t v, nv, d, x; usg_arc_t *av;
	if (g->a[v0>>1].del) return 0; // already deleted
    if(usg_arc_occ(g, v0, NULL) < 2) return 0;
	///S saves nodes with all incoming edges visited
	b->S.n = b->T.n = b->b.n = b->e.n = 0;
	///for each node, b->a saves all related information
	b->a[v0].c = b->a[v0].d = b->a[v0].m = b->a[v0].nc = b->a[v0].np = 0;
	///b->S is the nodes with all incoming edges visited
	kv_push(uint32_t, b->S, v0);
    n_tips = 0; tip_end = (uint32_t)-1;
    
    do {
		///v is a node that all incoming edges have been visited
		///d is the distance from v0 to v
		v = kv_pop(b->S); d = b->a[v].d;
		nv = usg_arc_n(g, v); av = usg_arc_a(g, v);
		///why we have this assert?
		///assert(nv > 0);
		///all out-edges of v
		for (i = 0; i < nv; ++i) { // loop through v's neighbors
			/**
			p->ul: |____________31__________|__________1___________|______________32_____________|
								qn            direction of overlap       length of this node (not overlap length)
												(in the view of query)
			p->v : |___________31___________|__________1___________|
								tn             reverse direction of overlap
											(in the view of target)
			p->ol: overlap length
			**/
            ///if this edge has been deleted
			if (av[i].del) continue;
            
			uint32_t w = av[i].v, l = (uint32_t)av[i].ul; // v->w with length l
			binfo_t *t = &b->a[w];
			///that means there is a circle, directly terminate the whole bubble poping
			///if (w == v0) goto pop_reset;
            if ((w>>1) == (v0>>1)) goto usg_clean_reset;
            /****************************may have bugs********************************/
            ///important when poping at long untig graph
            if(is_first) l = 0;
            /****************************may have bugs********************************/

			///push the edge
            ///high 32-bit of g->idx[v] is the start point of v's edges
            //so here is the point of this specfic edge
			// kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);
			///find a too far path? directly terminate the whole bubble poping
			if (d + l > max_dist) break; // too far
            if (b->b.n > max_occ) break; // too far

            ///if this node
			if (t->s == 0) { // this vertex has never been visited
				kv_push(uint32_t, b->b, w); // save it for revert
				///t->p is the parent node of 
				///t->s = 1 means w has been visited
				///d is len(v0->v), l is len(v->w), so t->d is len(v0->w)
				t->p = v, t->s = 1, t->d = d + l;
				///incoming edges of w
                t->r = usg_arc_occ(g, w^1, NULL);
				++n_pending;
			} else { // visited before
				if (d + l < t->d) t->d = d + l; // update dist
			}
			///assert(t->r > 0);
			//if all incoming edges of w have visited
			//push it to b->S
			if (--(t->r) == 0) {
                x = usg_arc_occ(g, w, NULL);
                /****************************may have bugs for bubble********************************/
                if(x > 0) {
                    kv_push(uint32_t, b->S, w);
                }
                else {
                    ///at most one tip
                    if(n_tips != 0) goto usg_clean_reset;
                    n_tips++; tip_end = w;
                }
                /****************************may have bugs for bubble********************************/
				--n_pending;
			}
		}
        is_first = 0;
        //if found a tip
        /****************************may have bugs for bubble********************************/
        if(n_tips == 1) {
            if(tip_end != (uint32_t)-1 && n_pending == 0 && b->S.n == 0) {
                kv_push(uint32_t, b->S, tip_end);
                break;
            } else {
                goto usg_clean_reset;
            }
        }
        /****************************may have bugs for bubble********************************/
		///if i < nv, that means (d + l > max_dist)
		if (i < nv || b->S.n == 0) goto usg_clean_reset;
	} while (b->S.n > 1 || n_pending);

    n_pop = 1;
    usg_clean_reset:

    for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
		binfo_t *t = &b->a[b->b.a[i]];
		t->s = t->c = t->d = t->m = t->nc = t->np = 0;
	}
	return n_pop;
}

void extracr_clean_arc(usg_t *g, uint32_t v, asg64_v *b, asg64_v *ub)
{
    usg_arc_t *av; uint32_t i, kv, nv; uint64_t x, kocc[2];
    av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
    if (nv < 2) return;
    for (i = kv = kocc[0] = kocc[1] = 0; i < nv; ++i) {
        if(av[i].del) continue;
        kv++; 
        if((av[i].ou>>1) > 0) { ///if av[i].ou == 1, ignore it
            if(kocc[1] < (av[i].ou>>1)) kocc[1] = (av[i].ou>>1);
        } else if(av[i].ou == 0) { 
            kocc[0]++;
        }            
    }
    if(kv < 2 || kocc[0] == 0 || kocc[1] == 0) return;
    for (i = 0; i < nv; ++i) {
        if(av[i].del || av[i].ou) continue;
        // if(max_drop_len && av[i].ol >= (*max_drop_len)) continue;
        x = av[i].ol; x <<= 32;
        kv_push(uint64_t, *b, ((x)|((uint64_t)(ub->n))));
        kv_push(uint64_t, *ub, ((((uint64_t)(v))<<32)|((uint64_t)(i))));   
    }
}


int usg_tip_del(usg_t *g, uint32_t v, asg64_v *z_a, asg64_v *z_b)
{
    uint64_t a_n = z_a->n, b_n = z_b->n, kv, nv, i; 
    int32_t n_ext = 0; usg_arc_t *av;

    av = usg_arc_a(g, v^1); nv = usg_arc_n(g, v^1);
    for (i = kv = 0; i < nv && kv <= 1; i++) {
        if (av[i].del || g->a[av[i].v>>1].del) continue;
        kv++; 
    }
    if(kv > 1) return 0;
    n_ext += g->a[v>>1].occ; g->a[v>>1].del = 1; kv_push(uint64_t, *z_b, v>>1);

    av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
    for (i = 0; i < nv; i++) {
        if (av[i].del || g->a[av[i].v>>1].del) continue;
        kv_push(uint64_t, *z_a, av[i].v);
    }   

    while (z_a->n > a_n) {
        v = z_a->a[--z_a->n]; if(g->a[v>>1].del) continue;

        av = usg_arc_a(g, v^1); nv = usg_arc_n(g, v^1);
        for (i = kv = 0; i < nv && kv < 1; i++) {
            if (av[i].del || g->a[av[i].v>>1].del) continue;
            kv++; 
        }
        if(kv > 0) continue;
        n_ext += g->a[v>>1].occ; g->a[v>>1].del = 1; kv_push(uint64_t, *z_b, v>>1);

        av = usg_arc_a(g, v); nv = usg_arc_n(g, v);
        for (i = 0; i < nv; i++) {
            if (av[i].del || g->a[av[i].v>>1].del) continue;
            kv_push(uint64_t, *z_a, av[i].v);
        }        
    }

    for (i = b_n; i < z_b->n; i++) usg_seq_del(g, z_b->a[i]);
    z_b->n = b_n; z_a->n = a_n;
    return n_ext;
}

uint64_t usg_bub_clean0(usg_t *g, uint32_t *a, uint32_t a_n, uint32_t sid, uint32_t eid, 
uint32_t is_topo, int32_t max_ext, float len_rat, uint8_t *f, asg64_v *in_0, asg64_v *in_1)
{
    uint32_t i, k, v, w, nv, nw, kv, kw, ol_max, ou_max, to_del, cnt = 0, tip_v, mm_ol; 
    asg64_v tx = {0,0,0}, tz = {0,0,0}, *b = NULL, *ub = NULL; int32_t n_tip, r_tip;
    uint64_t kocc[2]; usg_arc_t *av, *aw, *ve, *we;
    b = ((in_0)?(in_0):(&tx)); ub = ((in_1)?(in_1):(&tz));

    for (k = b->n = ub->n = 0; k < a_n; ++k) {
        extracr_clean_arc(g, a[k], b, ub); 
        extracr_clean_arc(g, a[k]^1, b, ub);
    }
    extracr_clean_arc(g, sid, b, ub); 
    extracr_clean_arc(g, eid, b, ub);
    

    radix_sort_srt64(b->a, b->a + b->n);
    for (k = 0; k < b->n; k++) {
        v = ub->a[(uint32_t)b->a[k]]>>32; 
        ve = &(usg_arc_a(g, v)[(uint32_t)(ub->a[(uint32_t)b->a[k]])]);
        w = ve->v^1;
        if(ve->del || g->a[v>>1].del || g->a[w>>1].del || ve->ou) continue;
        nv = usg_arc_n(g, v); nw = usg_arc_n(g, w);
        av = usg_arc_a(g, v); aw = usg_arc_a(g, w);
        if(nv<=1 && nw <= 1) continue;    

        // if(is_trio) {
        //     if(get_arcs(g, v, NULL, 0)<=1 && get_arcs(g, w, NULL, 0)<=1) continue;///speedup
        //     trioF = get_tip_trio_infor(g, v^1);
        //     ntrioF = (trioF==FATHER? MOTHER : (trioF==MOTHER? FATHER : (uint32_t)-1));
        // }    
        for (i = 0; i < nw; ++i) {
            if (aw[i].v == (v^1)) {
                we = &(aw[i]);
                break;
            }
        }
        mm_ol = MIN(ve->ol, we->ol); kocc[0] = kocc[1] = 0;

        for (i = kv = ol_max = ou_max = 0; i < nv; ++i) {
            if(av[i].del) continue;
            kv++; 
            if(av[i].ou != 1) kocc[!!(av[i].ou)]++; ///if av[i].ou == 1, ignore it
            // if(is_trio && get_tip_trio_infor(g, av[i].v) == ntrioF) continue;
            if(ol_max < av[i].ol) ol_max = av[i].ol;
        }
        if (kv < 1 || kocc[0] < 1 || kocc[1] < 1) continue;
        if (kv >= 2) {
            if (mm_ol > ol_max*len_rat) continue;
        }

        
        for (i = kw = ol_max = ou_max = 0; i < nw; ++i) {
            if(aw[i].del) continue;
            kw++;
            // if(is_trio && get_tip_trio_infor(g, aw[i].v) == ntrioF) continue;
            if(ol_max < aw[i].ol) ol_max = aw[i].ol;
        }
        if (kw < 1) continue;
        if (kw >= 2) {
            if (mm_ol > ol_max*len_rat) continue;
        }

        if (kv <= 1 && kw <= 1) continue;

        to_del = 1; n_tip = 0; tip_v = (uint32_t)-1;
        if(is_topo) {
            to_del = 0;
            if (kv > 1 && kw > 1) {
                to_del = 1;
            } else if (kw == 1) {
                n_tip = usg_naive_topocut_aux(g, w^1, max_ext, f, b, ub);
                if (n_tip < max_ext) to_del = 1; tip_v = w^1;
            } else if (kv == 1) {
                n_tip = usg_naive_topocut_aux(g, v^1, max_ext, f, b, ub);
                if (n_tip < max_ext) to_del = 1; tip_v = v^1;
            }
        } 

        if (to_del) {
            // ve->del = we->del = 1; 
            if(n_tip > 0 && tip_v != ((uint32_t)-1)) {
                r_tip = usg_tip_del(g, tip_v, b, ub);
                assert(n_tip == r_tip);
            }
            cnt++;
        }
    }

    if(in_0) free(tx.a); if(in_1) free(tz.a);
    return cnt;
}

uint64_t usg_bub_clean(usg_t *g, buf_t *b, asg64_v *in_0, asg64_v *in_1, int32_t max_ext, float len_rat, 
uint32_t is_topo, uint8_t *bs, uint8_t *f)
{
    uint32_t v, m, n_vtx = g->n<<1, n_arc, nv, i, n_cut = 0;
    uint64_t max_dist; usg_arc_t *av = NULL;
    // for (i = 0; i < g->n; ++i) b->a[i].s = 0;
    max_dist = usg_bub_dis(g, b); 
    memset(bs, 0, sizeof((*bs))*n_vtx);

    if(max_dist > 0) {
        for (v = 0; v < n_vtx; ++v) {
            if(bs[v]) continue;
            nv = usg_arc_n(g, v); av = usg_arc_a(g, v);
            if (nv < 2 || g->a[v>>1].del) continue;
            for (i = n_arc = 0; i < nv; ++i) // asg_bub_pop1() may delete some edges/arcs
                if (!av[i].del) ++n_arc;
            if (n_arc < 2) continue;
            if(usg_bub_identify(g, v, max_dist, (uint64_t)-1, b)) {
                //beg is v, end is b.S.a[0]
                //note b.b include end, does not include beg
                for (i = 0; i < b->b.n; i++) {
                    if(b->b.a[i]==v || b->b.a[i]==b->S.a[0]) continue;
                    bs[b->b.a[i]] = bs[b->b.a[i]^1] = 1;
                }
                bs[v] = 2; bs[b->S.a[0]^1] = 3;
            }
        }

        //traverse all node with two directions 
        for (v = 0; v < n_vtx; ++v) {
            if(bs[v] != 2) continue;
            nv = usg_arc_n(g, v); av = usg_arc_a(g, v);
            for (i = n_arc = 0; i < nv; ++i) // asg_bub_pop1() may delete some edges/arcs
                if (!av[i].del) ++n_arc;
            if (n_arc < 2) continue;
            if(usg_bub_identify(g, v, max_dist, (uint64_t)-1, b)) {
                for (i = m = 0; i < b->b.n; i++) {
                    if(b->b.a[i]==v || b->b.a[i]==b->S.a[0]) continue;
                    bs[b->b.a[i]] = bs[b->b.a[i]^1] = 1;
                    b->b.a[m++] = b->b.a[i];
                }
                b->b.n = m; bs[v] = 2; bs[b->S.a[0]^1] = 3;
                n_cut += usg_bub_clean0(g, b->b.a, b->b.n, v, b->S.a[0]^1, is_topo, max_ext, len_rat, f, in_0, in_1);
            }
        }
    }
    if (n_cut > 0) usg_cleanup(g);
    return n_cut;
}

void u2g_hybrid_clean(ul_resolve_t *uidx, ulg_opt_t *ulopt, usg_t *ng, asg64_v *b, asg64_v *ub)
{
    int64_t i, ss, mm_tip = ulopt->max_tip_hifi;///ulopt->max_tip; 
    double step = (ulopt->clean_round==1?ulopt->max_ovlp_drop_ratio:
                        ((ulopt->max_ovlp_drop_ratio-ulopt->min_ovlp_drop_ratio)/(ulopt->clean_round-1)));
    double drop = ulopt->min_ovlp_drop_ratio; ///CALLOC(iug->g->seq_vis, iug->g->n_seq*2);
    buf_t bb; memset(&bb, 0, sizeof(buf_t)); CALLOC(bb.a, ng->n<<1);
    uint8_t *bs, *f; CALLOC(bs, ng->n<<1); CALLOC(f, ng->n);
    // fprintf(stderr, "\n[M::%s::] Starting hybrid clean, mm_tip::%ld\n", __func__, mm_tip);
    // prt_usg_t(uidx, ng, "ng_h0");
    usg_arc_cut_tips(ng, mm_tip, 0, b);///p_telo
    // prt_usg_t(uidx, ng, "ng_h1");
    // char sb[1000]; 

    for (ss = 1; ss <= 1/**6**/; ss++) {
        mm_tip = ulopt->max_tip_hifi*ss;
        // fprintf(stderr, "\n[M::%s::] ss::%ld, mm_tip::%ld, ulopt->clean_round::%ld\n", 
        // __func__, ss, mm_tip, ulopt->clean_round);
        for (i = 0, drop = ulopt->min_ovlp_drop_ratio; i < ulopt->clean_round; i++, drop += step) {
            if(drop > ulopt->max_ovlp_drop_ratio) drop = ulopt->max_ovlp_drop_ratio;
            // fprintf(stderr, "-0-[M::%s::] i::%ld, drop::%f\n", __func__, i, drop);
            // sprintf(sb, "ng_ss::%ld_i::%ld_drop::%f_a", ss, i, drop);
            // prt_usg_t(uidx, ng, sb);
            // usg_arc_cut_length(ng, b, ub, mm_tip>>1, drop, ulopt->is_trio, 1, NULL);
            usg_bub_clean(ng, &bb, b, ub, mm_tip>>1, drop, 1, bs, f);
            usg_arc_cut_srt_length(ng, b, ub, mm_tip>>1, drop, ulopt->is_trio, 1, NULL, bs);
            // fprintf(stderr, "-1-[M::%s::] i::%ld, drop::%f\n", __func__, i, drop);
            // sprintf(sb, "ng_ss::%ld_i::%ld_drop::%f_b", ss, i, drop);
            // prt_usg_t(uidx, ng, sb);
            usg_arc_cut_tips(ng, mm_tip, 0, b);
            // fprintf(stderr, "-2-[M::%s::] i::%ld, drop::%f\n", __func__, i, drop);
            // sprintf(sb, "ng_ss::%ld_i::%ld_drop::%f_c", ss, i, drop);
            // prt_usg_t(uidx, ng, sb);
            // usg_arc_cut_length(ng, b, ub, mm_tip, drop, ulopt->is_trio, 1, NULL);
            usg_bub_clean(ng, &bb, b, ub, mm_tip, drop, 1, bs, f);
            usg_arc_cut_srt_length(ng, b, ub, mm_tip, drop, ulopt->is_trio, 1, NULL, bs);
            // fprintf(stderr, "-3-[M::%s::] i::%ld, drop::%f\n", __func__, i, drop);
            // sprintf(sb, "ng_ss::%ld_i::%ld_drop::%f_d", ss, i, drop);
            // prt_usg_t(uidx, ng, sb);
            usg_arc_cut_tips(ng, mm_tip, 1, b);
            // sprintf(sb, "ng_ss::%ld_i::%ld_drop::%f_e", ss, i, drop);
            // prt_usg_t(uidx, ng, sb);
            // fprintf(stderr, "-4-[M::%s::] i::%ld, drop::%f\n", __func__, i, drop);
        }

        drop = 1;
        // fprintf(stderr, "-0-[M::%s::] i::%ld, drop::%f\n", __func__, i, drop);
        // sprintf(sb, "ng_ss::%ld_i::%ld_drop::%f_a", ss, i, drop);
        // prt_usg_t(uidx, ng, sb);
        // usg_arc_cut_length(ng, b, ub, mm_tip>>1, drop, ulopt->is_trio, 1, NULL);
        usg_bub_clean(ng, &bb, b, ub, mm_tip>>1, drop, 1, bs, f);
        // fprintf(stderr, "-1bub-[M::%s::] i::%ld, drop::%f\n", __func__, i, drop);
        // sprintf(sb, "ng_ss::%ld_i::%ld_drop::%f_b::bub", ss, i, drop);
        // prt_usg_t(uidx, ng, sb);
        usg_arc_cut_srt_length(ng, b, ub, mm_tip>>1, drop, ulopt->is_trio, 1, NULL, bs);
        // fprintf(stderr, "-1-[M::%s::] i::%ld, drop::%f\n", __func__, i, drop);
        // sprintf(sb, "ng_ss::%ld_i::%ld_drop::%f_b", ss, i, drop);
        // prt_usg_t(uidx, ng, sb);
        usg_arc_cut_tips(ng, mm_tip, 0, b);
        // fprintf(stderr, "-2-[M::%s::] i::%ld, drop::%f\n", __func__, i, drop);
        // sprintf(sb, "ng_ss::%ld_i::%ld_drop::%f_c", ss, i, drop);
        // prt_usg_t(uidx, ng, sb);
        // usg_arc_cut_length(ng, b, ub, mm_tip, drop, ulopt->is_trio, 1, NULL);
        usg_bub_clean(ng, &bb, b, ub, mm_tip, drop, 1, bs, f);
        usg_arc_cut_srt_length(ng, b, ub, mm_tip, drop, ulopt->is_trio, 1, NULL, bs);
        // fprintf(stderr, "-3-[M::%s::] i::%ld, drop::%f\n", __func__, i, drop);
        // sprintf(sb, "ng_ss::%ld_i::%ld_drop::%f_d", ss, i, drop);
        // prt_usg_t(uidx, ng, sb);
        usg_arc_cut_tips(ng, mm_tip, 1, b);
        // sprintf(sb, "ng_ss::%ld_i::%ld_drop::%f_e", ss, i, drop);
        // prt_usg_t(uidx, ng, sb);
        // fprintf(stderr, "-4-[M::%s::] i::%ld, drop::%f\n", __func__, i, drop);
    }
    free(bb.a); free(bb.S.a); free(bb.T.a); free(bb.b.a); free(bb.e.a); free(bs); free(f);
    // prt_usg_t(uidx, ng, "ng_h2");
    ///debug
    debug_sysm_usg_t(ng, __func__);

    /******for debug******/
    // prt_usg_t(uidx, ng, "ng_dbg");
    /******for debug******/

    // u2g_hybrid_extend(ng, NULL, b, ub);
    // u2g_hybrid_detan(uidx, ng, mm_tip, b, ub);
    u2g_hybrid_detan_iter(uidx, ng, mm_tip, ulopt->clean_round, b, ub);
}

void merge_hybrid_utg_content(ma_utg_t* cc, ma_ug_t* raw, asg_t* rg, usg_t *ng, kvec_asg_arc_t_warp* edge)
{
    if(cc->m == 0) return;
    uint32_t i, j, index, uId, ori, uv, uw, bv, bw;
    uint64_t tot, z; asg_arc_t *p; 
    ma_utg_t* q = NULL;
    for (i = index = 0; i < cc->n; i++) {
        z = (ng->a[cc->a[i]>>33].mm<<1)|((cc->a[i]>>32)&1); z <<= 32; z += ((uint32_t)cc->a[i]);
        index += ng->a[cc->a[i]>>33].occ; cc->a[i] = z;
    }

    uint64_t *buffer, *aim = NULL; MALLOC(buffer, index);
    for (i = index = edge->a.n = 0; i < cc->n; i++) {
        uId = cc->a[i]>>33; ori = cc->a[i]>>32&1;
        q = &(raw->u.a[uId]); aim = buffer + index;
        if(ori == 1) {
            for (j = 0; j < q->n; j++) {
                aim[q->n - j - 1] = (q->a[j])^(uint64_t)(0x100000000);
            }
        } else {
            for (j = 0; j < q->n; j++) {
                aim[j] = q->a[j];
            }
        }
        index += q->n;
        if(i > 0) {
            uv = cc->a[i-1]>>32; uw = cc->a[i]>>32;
            p = get_specfic_edge(raw->g, uv, uw);
            if(!p) {
                uv = cc->a[i-1]>>32; uw = cc->a[i]>>32;
                bv = raw->u.a[uv>>1].a[((uv&1)?(0):(raw->u.a[uv>>1].n-1))]>>32; if(uv&1) bv ^= 1;
                bw = raw->u.a[uw>>1].a[((uw&1)?(raw->u.a[uw>>1].n-1):(0))]>>32; if(uw&1) bw ^= 1;
                kv_pushp(asg_arc_t, edge->a, &p); memset(p, 0, sizeof((*p)));
                p->ul = bv; p->ul <<= 32; p->ul += rg->seq[bv>>1].len; p->v = bw;

                uv = (cc->a[i]>>32)^1; uw = (cc->a[i-1]>>32)^1;
                bv = raw->u.a[uv>>1].a[((uv&1)?(0):(raw->u.a[uv>>1].n-1))]>>32; if(uv&1) bv ^= 1;
                bw = raw->u.a[uw>>1].a[((uw&1)?(raw->u.a[uw>>1].n-1):(0))]>>32; if(uw&1) bw ^= 1;
                kv_pushp(asg_arc_t, edge->a, &p); memset(p, 0, sizeof((*p)));
                p->ul = bv; p->ul <<= 32; p->ul += rg->seq[bv>>1].len; p->v = bw;
            }
        }
    }

    if(index == 0) return;

    fill_unitig(buffer, index, rg, edge, (cc->n == 1 && raw->u.a[cc->a[0]>>33].circ), &tot);

    ///important. must be here
    if(cc->n == 1 && raw->u.a[cc->a[0]>>33].circ) cc->circ = 1;
    
    free(cc->a);
    cc->a = buffer; cc->n = cc->m = index; cc->len = tot;
    if(!cc->circ) {
        cc->start = cc->a[0]>>32;
        cc->end = (cc->a[cc->n-1]>>32)^1;
    }
    else {
        cc->start = cc->end = UINT32_MAX;
    }
}

ma_ug_t *gen_hybrid_ug(ul_resolve_t *uidx, usg_t *ng)
{
    // fprintf(stderr, "[M::%s::] ng->n::%u\n", __func__, (uint32_t)ng->n);
    ma_ug_t *ug = ma_ug_hybrid_gen(ng); 
    // fprintf(stderr, "[M::%s::] ug->g->n_seq::%u\n", __func__, (uint32_t)ug->g->n_seq);
    uint32_t i; ma_utg_t *u; kvec_asg_arc_t_warp e; kv_init(e.a); e.i = 0;
    for (i = 0; i < ug->u.n; i++) {
        ug->g->seq[i].c = PRIMARY_LABLE;
        u = &(ug->u.a[i]);
        if(u->m == 0) continue;
        // fprintf(stderr, "+[M::%s::] i::%u\n", __func__, i);
        merge_hybrid_utg_content(u, uidx->l1_ug, uidx->sg, ng, &e);
        // fprintf(stderr, "-[M::%s::] i::%u\n", __func__, i);
        ug->g->seq[i].len = u->len;
    }
    kv_destroy(e.a);
    return ug;
}

uint32_t gen_utg_telo(ma_utg_t *u, telo_end_t *te)
{
    uint32_t i;
    for (i = 0; i < u->n; ++i) {
        if(te->hh[u->a[i]>>33]) return 1;
    }
    return 0;
}

void renew_ul2_utg(ul_resolve_t *uidx);

void u2g_threading(ul_resolve_t *uidx, ulg_opt_t *ulopt, uint64_t cov_cutoff, uint64_t is_bridg, asg64_v *b, asg64_v *ub)
{
    renew_ul2_utg(uidx);
    ul2ul_idx_t *idx = &(uidx->uovl); ma_ug_t *iug = idx->i_ug; ma_ug_t *raw = uidx->l1_ug;
    uint64_t k, z, i, t_s, t_e; uinfo_srt_warp_t *seq; 
    usg_t *ng; CALLOC(ng, 1); usg_seq_t *s; usg_arc_warp *sv; usg_arc_t *p; 
    asg_arc_t *av; uint32_t nv, v, w, vx, wx; int64_t tt, tl, tm; 

    ng->mp.n = ng->mp.m = raw->g->n_seq; MALLOC(ng->mp.a, ng->mp.n);
    for (k = 0; k < raw->g->n_seq; k++) {
        s = push_usg_t_node(ng, k); 
        s->mm = k; s->arc[0].n = s->arc[1].n = 0; s->occ = raw->u.a[k].n;
        s->arc_mm[0].n = s->arc_mm[1].n;
        s->del = raw->g->seq[k].del; s->len = raw->g->seq[k].len;
        s->telo = 0; 
        if(uidx->uopt->te) {
            s->telo = gen_utg_telo(&(raw->u.a[k]), uidx->uopt->te);
        }

        av = asg_arc_a(raw->g, (k<<1)); nv = asg_arc_n(raw->g, (k<<1)); sv = &(s->arc[0]);
        for (z = 0; z < nv; z++) {
            if(av[z].del) continue;
            kv_pushp(usg_arc_t, *sv, &p); 
            p->del = 0; p->ou = 0; p->v = av[z].v; p->ol = av[z].ol; p->ul = av[z].ul; p->idx = 0;
        }

        av = asg_arc_a(raw->g, ((k<<1)+1)); nv = asg_arc_n(raw->g, ((k<<1)+1)); sv = &(s->arc[1]);
        for (z = 0; z < nv; z++) {
            if(av[z].del) continue;
            kv_pushp(usg_arc_t, *sv, &p); 
            p->del = 0; p->ou = 0; p->v = av[z].v; p->ol = av[z].ol; p->ul = av[z].ul; p->idx = 0;
        }
        ///map: ng id -> raw id
        ng->mp.a[k].n = ng->mp.a[k].m = 1; MALLOC(ng->mp.a[k].a, 1); ng->mp.a[k].a[0] = k;
    }

    for (k = b->n = 0; k < iug->u.n; k++) {
        seq = &(idx->cc.iug_a[k]);
        if(seq->n <= 1) continue;

        uidx->uovl.iug_seq = seq; uidx->uovl.iug_cov_thre = cov_cutoff;
        for (z = 0; z < uidx->str_b.n_thread; z++) {
            uidx->str_b.buf[z].res_dump.n = 0;
        }
        kt_for(uidx->str_b.n_thread, worker_update_ul_arc_drop, uidx, seq->n-1);///seq->n > 1
        for (z = ub->n = 0; z < uidx->str_b.n_thread; z++) {
            for (i = 0; i < uidx->str_b.buf[z].res_dump.n; i++) {
                kv_push(uint64_t, *ub, uidx->str_b.buf[z].res_dump.a[i]);
            }
        }

        for (z = 0; z < ub->n; z++) {///all unreliable arcs
            i = ub->a[z]; tm = 1; p = get_usg_arc(ng, seq->a[i].v, seq->a[i+1].v);
            if(p) {
                if(p->ou < (uint64_t)tm) p->ou = tm;///since the initial ou is 0, p->ou = 1
                pushp_usg_arc_mm(ng, seq->a[i].v, seq->a[i+1].v, k, i); 

                p = get_usg_arc(ng, seq->a[i+1].v^1, seq->a[i].v^1);
                if(p->ou < (uint64_t)tm) p->ou = tm;///since the initial ou is 0, p->ou = 1
                pushp_usg_arc_mm(ng, seq->a[i+1].v^1, seq->a[i].v^1, k, i); 
            } 
            ///give up unreliable arcs if they are not adjacent
        }

        radix_sort_srt64(ub->a, ub->a + ub->n); 
        if(ub->n == 0 || ub->a[ub->n-1] < seq->n-1) kv_push(uint64_t, *ub, seq->n-1);
        for (z = t_s = t_e = 0; z < ub->n; z++) {
            t_e = ub->a[z];
            if(t_s < t_e) {
                for (i = t_s, tt = seq->a[t_e].n; i < t_e; i++) {
                    tt += seq->a[i].n;///how many HiFi reads are covered
                }
                for (i = t_s, tl = 0; i < t_e; i++) {
                    tl += seq->a[i].n; 
                    tm = tt - tl; if(tm > tl) tm = tl; assert(tm > 0); tm <<= 1; tm += 1;
                    p = get_usg_arc(ng, seq->a[i].v, seq->a[i+1].v);
                    if(p) {
                        if(p->ou < (uint64_t)tm) p->ou = tm;
                        p = get_usg_arc(ng, seq->a[i+1].v^1, seq->a[i].v^1);
                        if(p->ou < (uint64_t)tm) p->ou = tm;
                    } else if(is_bridg) {
                        v = seq->a[i].v; w = seq->a[i+1].v; p = NULL;
                        vx = (v&1?((raw->u.a[v>>1].a[0]>>32)^1):(raw->u.a[v>>1].a[raw->u.a[v>>1].n-1]>>32));
                        wx = (w&1?((raw->u.a[w>>1].a[raw->u.a[w>>1].n-1]>>32)^1):(raw->u.a[w>>1].a[0]>>32));
                        ///(v>>1) != (w>>1) or v == w, need to handle (v>>1) == (w>>1) && (v&1) != (w&1) later
                        if(((v^w) != 1) && (uidx->sg->seq[vx>>1].len > uidx->uopt->min_ovlp) && 
                                (uidx->sg->seq[wx>>1].len > uidx->uopt->min_ovlp)) {
                            kv_pushp(usg_arc_t, (ng->a[v>>1].arc[v&1]), &p); 
                            p->del = 0; p->ou = tm; p->v = w; p->ol = 0; p->idx = 0;
                            p->ul = (((uint64_t)v)<<32)|(raw->g->seq[v>>1].len);

                            v = seq->a[i+1].v^1; w = seq->a[i].v^1; 
                            kv_pushp(usg_arc_t, (ng->a[v>>1].arc[v&1]), &p); 
                            p->del = 0; p->ou = tm; p->v = w; p->ol = 0; p->idx = 0;
                            p->ul = (((uint64_t)v)<<32)|(raw->g->seq[v>>1].len);
                        }
                    }
                    if(p) {
                        pushp_usg_arc_mm(ng, seq->a[i].v, seq->a[i+1].v, k, i); 
                        pushp_usg_arc_mm(ng, seq->a[i+1].v^1, seq->a[i].v^1, k, i); 
                    }
                }
            }
            t_s = ub->a[z] + 1;
        }
    }

    // for (v = 0; v < (ng->n<<1); v++) {
    //     p = usg_arc_a(ng, v); nv = usg_arc_n(ng, v);
    //     for (i = 0; i < nv; i++) {
    //         assert(get_usg_arc(ng, p[i].v^1, v^1));
    //     }
    // }

    usg_cleanup(ng);

    ///debug
    debug_sysm_usg_t(ng, __func__);

    idx->h_usg = ng;

    /******for debug******/
    // prt_usg_t(uidx, ng, "ng0");
    /******for debug******/



    u2g_hybrid_clean(uidx, ulopt, ng, b, ub);
    // idx->hybrid_ug = gen_hybrid_ug(uidx, ng);
    
    
    // idx->hybrid_ug = gen_debug_hybrid_ug(uidx, ng);
    // fprintf(stderr, "-[M::%s::] idx->hybrid_ug->g->n_seq::%u\n", __func__, (uint32_t)idx->hybrid_ug->g->n_seq);
}

void u2g_clean(ul_resolve_t *uidx, ulg_opt_t *ulopt, uint32_t keep_raw_utg, uint64_t is_bridg)
{
    ul2ul_idx_t *idx = &(uidx->uovl); asg64_v bu = {0,0,0}, uu = {0,0,0}; int64_t max_tip_hifi0 = ulopt->max_tip_hifi;
    ma_ug_t *iug = idx->i_ug; int64_t i, mm_tip = ulopt->max_tip; uint64_t cnt = 1, topo_level, ss = 0;
    double step = (ulopt->clean_round==1?ulopt->max_path_drop_ratio:
                        ((ulopt->max_path_drop_ratio-ulopt->min_path_drop_ratio)/(ulopt->clean_round-1)));
    double drop = ulopt->min_path_drop_ratio; CALLOC(iug->g->seq_vis, iug->g->n_seq*2);


    // char sb[1000]; 
    // output_integer_graph(uidx, iug, "ig_h0", 0);

    // fprintf(stderr, "\n[M::%s::] max_path_drop_ratio::%f, min_path_drop_ratio::%f, max_tip_hifi::%ld, max_tip::%ld\n", 
    // __func__, ulopt->max_path_drop_ratio, ulopt->min_path_drop_ratio, 
    // ulopt->max_tip_hifi, ulopt->max_tip);
    for (ss = 0; ss < 2; ss++) {
        for (i = 0, drop = ulopt->min_path_drop_ratio; i < ulopt->clean_round; i++, drop += step) {
            if(drop > ulopt->max_path_drop_ratio) drop = ulopt->max_path_drop_ratio;
            // fprintf(stderr, "[M::%s::] Starting round-%ld, drop::%f, max_tip_hifi::%ld\n", 
            // __func__, i, drop, ulopt->max_tip_hifi);
            cnt = 1; topo_level = 2; mm_tip = ulopt->max_tip;
            while (cnt) {
                cnt = 0;
                asg_arc_identify_simple_bubbles_multi(iug->g, ulopt->b_mask_t, 0);
                cnt += ulg_arc_cut_supports(uidx, iug, mm_tip, ulopt->max_tip_hifi, drop, ulopt->is_trio, topo_level, 1, NULL, keep_raw_utg, &bu, &uu);
                cnt += ulg_arc_cut_tips(uidx, iug, mm_tip, ulopt->max_tip_hifi, 1, &bu, &uu);///p_telo
            }
            
            cnt = 1; topo_level = 1; mm_tip = ulopt->max_tip;
            while (cnt) {
                cnt = 0;
                asg_arc_identify_simple_bubbles_multi(iug->g, ulopt->b_mask_t, 0);
                cnt += ulg_arc_cut_supports(uidx, iug, mm_tip, ulopt->max_tip_hifi, drop, ulopt->is_trio, topo_level, 1, NULL, keep_raw_utg, &bu, &uu);
                cnt += ulg_arc_cut_tips(uidx, iug, mm_tip, ulopt->max_tip_hifi, 1, &bu, &uu);
            }

            cnt = 1; topo_level = 1; mm_tip = ((int64_t)0x7fffffff);
            while (cnt) {
                cnt = 0;
                asg_arc_identify_simple_bubbles_multi(iug->g, ulopt->b_mask_t, 0);
                cnt += ulg_arc_cut_supports(uidx, iug, mm_tip, ulopt->max_tip_hifi, drop, ulopt->is_trio, topo_level, 1, NULL, keep_raw_utg, &bu, &uu);
                cnt += ulg_arc_cut_tips(uidx, iug, mm_tip, ulopt->max_tip_hifi, 0, &bu, &uu);
            }
            // fprintf(stderr, "[M::%s::] Done round-%ld, drop::%f\n", __func__, i, drop);
            // sprintf(sb, "ig_ss_%lu_i_%ld", ss, i);
            // output_integer_graph(uidx, iug, sb, 0);
        }

        ulg_pop_bubble(uidx, iug, NULL, ((int64_t)0x7fffffff), ulopt->max_tip_hifi, 1, &bu, &uu);

        while(ulg_arc_cut_z(uidx, iug, ((int64_t)0x7fffffff), ulopt->max_tip_hifi, 1.5, 0.15, 100, 0.8, ulopt->is_trio, 1, NULL, &bu, &uu));///non_p_telo

        // sprintf(sb, "ig_ss_%lu_i_%ld", ss, i);
        // output_integer_graph(uidx, iug, sb, 0);

        ulopt->max_tip_hifi <<= 1;
    }
    ulopt->max_tip_hifi = max_tip_hifi0;

    // output_integer_graph(uidx, iug, "ig_h1", 0);

    u2g_threading(uidx, ulopt, 3, is_bridg, &bu, &uu);


    // fill_u2g(uidx, &bu, &uu);

    // while (cnt) {
        // for (i = cnt = 0, mm_tip = ulopt->max_tip, drop = ulopt->min_ovlp_drop_ratio; i < ulopt->clean_round; i++, drop += step) {
        //     if(drop > ulopt->max_ovlp_drop_ratio) drop = ulopt->max_ovlp_drop_ratio;
            // cnt += ulg_arc_cut_occ(uidx, iug, mm_tip, ulopt->max_tip_hifi, ulopt->is_trio, 2, &bu, &uu);
            // asg_arc_identify_simple_bubbles_multi(iug->g, ulopt->b_mask_t, 0);
            // cnt += ulg_arc_cut_length(uidx, iug, mm_tip, ulopt->max_tip_hifi, drop/**MIN(drop, ulopt->hom_check_drop_rate)**/, ulopt->is_trio, 2, 1, NULL, &bu, &uu);
            // cnt += ulg_arc_cut_tips(uidx, iug, mm_tip, ulopt->max_tip_hifi, &bu, &uu);

            // asg_arc_identify_simple_bubbles_multi(iug->g, ulopt->b_mask_t, 0);
            // cnt += ulg_arc_cut_bridge(uidx, iug, mm_tip, ulopt->max_tip_hifi, 0.5, ulopt->is_trio, 2, NULL, &bu, &uu);

            // asg_arc_identify_simple_bubbles_multi(iug->g, ulopt->b_mask_t, 0);
            // cnt += ulg_arc_cut_misjoin(uidx, iug, mm_tip, ulopt->max_tip_hifi, 0.5, ulopt->is_trio, 2, NULL, &bu, &uu);
        // }

        // for (i = cnt = 0, mm_tip = ((int64_t)0x7fffffff), drop = ulopt->min_ovlp_drop_ratio; i < ulopt->clean_round; i++, drop += step) {
        //     if(drop > ulopt->max_ovlp_drop_ratio) drop = ulopt->max_ovlp_drop_ratio;
            // cnt += ulg_arc_cut_occ(uidx, iug, mm_tip, ulopt->max_tip_hifi, ulopt->is_trio, 2, &bu, &uu);
            // asg_arc_identify_simple_bubbles_multi(iug->g, ulopt->b_mask_t, 0);
            // cnt += ulg_arc_cut_length(uidx, iug, mm_tip, ulopt->max_tip_hifi, drop/**MIN(drop, ulopt->hom_check_drop_rate)**/, ulopt->is_trio, 2, 1, NULL, &bu, &uu);
            // cnt += ulg_arc_cut_tips(uidx, iug, mm_tip, ulopt->max_tip_hifi, &bu, &uu);

            // asg_arc_identify_simple_bubbles_multi(iug->g, ulopt->b_mask_t, 0);
            // cnt += ulg_arc_cut_bridge(uidx, iug, mm_tip, ulopt->max_tip_hifi, 0.5, ulopt->is_trio, 2, NULL, &bu, &uu);

            // asg_arc_identify_simple_bubbles_multi(iug->g, ulopt->b_mask_t, 0);
            // cnt += ulg_arc_cut_misjoin(uidx, iug, mm_tip, ulopt->max_tip_hifi, 0.5, ulopt->is_trio, 2, NULL, &bu, &uu);
        // }
    // }

    free(bu.a); free(uu.a); 
}

void clc_contain(ul_resolve_t *uidx, uint64_t id, uint64_t is_ul, integer_t *buf)
{
    ma_ug_t *raw = uidx->l1_ug; ma_utg_t *ru; ug_opt_t *uopt = uidx->uopt;
    uint64_t k, m, qn, tn; ma_hit_t_alloc *x; asg_arc_t e; ul2ul_item_t *o;
    int64_t min_ovlp = uopt->min_ovlp, max_hang = uopt->max_hang, r;

    if(is_ul) {
        o = get_ul_ovlp(&(uidx->uovl), id, 1);
        assert(o && (!o->is_del));
        for (k = 0; k < o->cn; k++) {
            if((o->a[k].is_del) || (!ulg_type(uidx->uovl, o->a[k].hid))) continue;
            r = integer_hit2arc(&(o->a[k]), ulg_len(*uidx, id), ulg_len(*uidx, o->a[k].hid), 
                                        ulg_occ(*uidx, id), ulg_occ(*uidx, o->a[k].hid), id, o->a[k].hid, 0, NULL);
            if(r != MA_HT_TCONT) continue;
            kv_push(uint64_t, buf->o, o->a[k].hid);
        }
    } else {
        ru = &(raw->u.a[id]); 
        for (k = 0; k < ru->n; k++) {
            x = &(uopt->sources[ru->a[k]>>33]);
            for (m = 0; m < x->length; m++) {
                qn = Get_qn(x->buffer[m]); tn = Get_tn(x->buffer[m]);
                r = ma_hit2arc(&(x->buffer[m]), Get_READ_LENGTH(R_INF, qn), Get_READ_LENGTH(R_INF, tn), 
                                                    max_hang, asm_opt.max_hang_rate, min_ovlp, &e);
                if(r != MA_HT_TCONT) continue;
                kv_push(uint64_t, buf->u, tn);
            }
        }
    }
}

static void worker_renew_u2g_cov(void *data, long i, int tid) // callback for kt_for()
{
    ul_resolve_t *uidx = (ul_resolve_t *)data; 
    integer_t *buf = &(uidx->str_b.buf[tid]);
    ma_utg_t *iu = &(uidx->uovl.i_ug->u.a[i]);
    uint64_t z, ri, k, ul_n, ug_n, *a, a_n; ma_ug_t *raw = uidx->l1_ug;
    ul_cov_t *idx = &(uidx->uovl.cc); 

    idx->uc[i] = idx->raw_uc[i] = idx->hc[i] = ul_n = ug_n = 0; buf->u.n = buf->o.n = 0;
    for (z = 0; z < iu->n; z++) {
        ri = iu->a[z]>>33;
        if(ulg_type(uidx->uovl, ri)) {///ul
            clc_contain(uidx, ulg_id(uidx->uovl, ri), 1, buf); ul_n++;
        } else {///ug node
            clc_contain(uidx, ulg_id(uidx->uovl, ri), 0, buf); 
            ug_n += raw->u.a[ulg_id(uidx->uovl, ri)].n;///HiFi reads occ
        }
    }
    idx->raw_uc[i] = ul_n;

    ///buf->o: UL contained reads
    ///buf->u: HiFi reads
    if(buf->o.n > 0) {
        a = buf->o.a; a_n = buf->o.n;
        radix_sort_srt64(a, a + a_n);
        for (z = 0, k = 1; k <= a_n; k++) {
            if(k == a_n || a[z] != a[k]) {
                ul_n++; z = k;
            }
        }
    }
    if(buf->u.n > 0) {
        a = buf->u.a; a_n = buf->u.n;
        radix_sort_srt64(a, a + a_n);
        for (z = 0, k = 1; k <= a_n; k++) {
            if(k == a_n || a[z] != a[k]) {
                ug_n++; z = k;
            }
        }
    }
    idx->uc[i] = ul_n; idx->hc[i] = ug_n;
}

ul2ul_t* get_ulg_spec_ovlp(ul_resolve_t *uidx, uint64_t uv, uint64_t uw)
{
    ma_ug_t *iug = uidx->uovl.i_ug; uint32_t qid, tid;
    qid = iug->u.a[uv>>1].a[((uv&1)?(0):(iug->u.a[uv>>1].n-1))]>>33;
    tid = iug->u.a[uw>>1].a[((uw&1)?(iug->u.a[uw>>1].n-1):(0))]>>33;

    if(uidx->uovl.item_idx[qid] == (uint32_t)-1) return NULL;
    ul2ul_item_t *o = &(uidx->uovl.a[uidx->uovl.item_idx[qid]]); uint64_t k;
    for (k = 0; k < o->cn; k++) {
        if(o->a[k].hid != tid) continue;
        return &(o->a[k]);
    }
    return NULL;
}





void gen_raw_ug_seq(ul_resolve_t *uidx, ul_str_t *str, ma_utg_t *u, ma_ug_t *raw, uinfo_srt_warp_t *res, uint64_t iug_id)
{
    uint64_t k, cd, nd, rev, x, os, oe; uinfo_srt_t *p; ma_utg_t *ru;
    uc_block_t *xi; ul2ul_t *z = NULL; ul_str_t *c_str; int64_t m, cs, ce;
    for (k = 0, res->n = 0; k < u->n; k++) {
        cd = u->a[k]>>33; rev = ((u->a[k]>>32)&1); ru = NULL; c_str = NULL;
        if(ulg_type(uidx->uovl, cd)) c_str = &(str[ulg_id(uidx->uovl, cd)]); ///ul
        else ru = &(raw->u.a[ulg_id(uidx->uovl, cd)]); //ug

        if(c_str) {
            if(k + 1 < u->n) {
                nd = u->a[k+1]>>33;
                z = get_ul_spec_ovlp(&(uidx->uovl), cd, nd);
                assert(z && (!z->is_del));
                if(!rev) {
                    cs = 0; ce = z->qs_k + 1;
                } else {
                    cs = z->qe_k - 1; ce = c_str->cn;
                }
            } else {
                cs = 0; ce = c_str->cn;
            }

            if(!rev) {
                xi = &(uidx->idx->a[cd].bb.a[c_str->a[cs]>>32]);
                x = (uint32_t)c_str->a[cs];
            } else {
                xi = &(uidx->idx->a[cd].bb.a[c_str->a[ce-1]>>32]);
                x = (uint32_t)c_str->a[ce-1]; x ^= 1;
            }

            if(res->n) {
                assert((res->a[res->n-1].v) == x);
                os = MAX(xi->ts, res->a[res->n-1].s); oe = MIN(xi->te, res->a[res->n-1].e);
                assert(oe > os);
                os = MIN(xi->ts, res->a[res->n-1].s); oe = MAX(xi->te, res->a[res->n-1].e);
                res->a[res->n-1].s = os; res->a[res->n-1].e = oe;
                if(!rev) cs++;
                else ce--;
            }

            if(!rev) {
                for (m = cs; m < ce; m++) {
                    xi = &(uidx->idx->a[cd].bb.a[c_str->a[m]>>32]);
                    x = (uint32_t)c_str->a[m];
                    kv_pushp(uinfo_srt_t, *res, &p);
                    p->v = x; p->s = xi->ts; p->e = xi->te;
                }
            } else {
                for (m = ce-1; m >= cs; m--) {
                    xi = &(uidx->idx->a[cd].bb.a[c_str->a[m]>>32]);
                    x = ((uint32_t)c_str->a[m])^1;
                    kv_pushp(uinfo_srt_t, *res, &p);
                    p->v = x; p->s = xi->ts; p->e = xi->te;
                }
            }
        } else if(ru) {
            x = ulg_id(uidx->uovl, cd); x <<= 1; if(rev) x^=1;
            if(res->n) {
                if((!ulg_type(uidx->uovl, (u->a[k-1]>>33)))) {///in previous ul is also a raw utg node
                    assert(get_specfic_edge(raw->g, res->a[res->n-1].v, x));
                    kv_pushp(uinfo_srt_t, *res, &p);
                    p->v = x; p->s = 0; p->e = ru->len;
                } else {
                    // assert(((res->a[res->n-1].v) == x));
                    if(((res->a[res->n-1].v) == x)) {
                        res->a[res->n-1].s = 0; res->a[res->n-1].e = ru->len;
                    } else {
                        assert(get_specfic_edge(raw->g, res->a[res->n-1].v, x));
                        kv_pushp(uinfo_srt_t, *res, &p);
                        p->v = x; p->s = 0; p->e = ru->len;
                    }
                }
            } else {
                kv_pushp(uinfo_srt_t, *res, &p);
                p->v = x; p->s = 0; p->e = ru->len;
            }
        }
    }


    for (k = 0; k < res->n; k++) {
        res->a[k].n = ug_occ_w(res->a[k].s, res->a[k].e, &(raw->u.a[res->a[k].v>>1]));
    }
}


void renew_u2g_cov(ul_resolve_t *uidx, telo_end_t *te)
{
    ul2ul_idx_t *idx = &(uidx->uovl); uinfo_srt_warp_t *x;
    ma_ug_t *i_ug = idx->i_ug, *raw = uidx->l1_ug; uint64_t k, z, iug_occ, m, l, *a, a_n; 
    free(idx->cc.uc); free(idx->cc.hc); free(idx->cc.raw_uc); 
    free(idx->cc.iug_idx); free(idx->cc.iug_b); free(idx->cc.iug_a); 
    if(te) {
        free(idx->telo); idx->telo = NULL;
    }
    memset(&(idx->cc), 0, sizeof(idx->cc));

    MALLOC(idx->cc.uc, i_ug->u.n); ///number of ul (contained+non-contained)
    MALLOC(idx->cc.raw_uc, i_ug->u.n); ///number of ul (non-contained)
    MALLOC(idx->cc.hc, i_ug->u.n); ///number of HiFi reads
    if(te) MALLOC(idx->telo, i_ug->u.n); ///is telo

    kt_for(uidx->str_b.n_thread, worker_renew_u2g_cov, uidx, i_ug->u.n);

    CALLOC(idx->cc.iug_a, i_ug->u.n); ///integer sequence of each integer unitig
    CALLOC(idx->cc.iug_idx, raw->u.n+1); ///idx for integer sequences
    // fprintf(stderr, "malloc::[M::%s::] i_ug->u.n:%u\n", __func__, (uint32_t)i_ug->u.n);
    for (k = iug_occ = 0; k < i_ug->u.n; k++) {
        // fprintf(stderr, "-[M::%s::] k::%lu, i_ug->u.n:%u\n", __func__, k, (uint32_t)i_ug->u.n);
        gen_raw_ug_seq(uidx, uidx->pstr.str.a, &(i_ug->u.a[k]), raw, &(idx->cc.iug_a[k]), k);
        iug_occ += idx->cc.iug_a[k].n; x = &(idx->cc.iug_a[k]);
        for (z = 0; z < x->n; z++) idx->cc.iug_idx[x->a[z].v>>1]++;

        if(te) {
            idx->telo[k] = 0;
            for (z = 0; z < x->n; z++) {
                if(gen_utg_telo(&(raw->u.a[x->a[z].v>>1]), te)) break;
            }
            if(z < x->n) idx->telo[k] = 1;
        }
    }

    MALLOC(idx->cc.iug_b, iug_occ); ///idx for integer sequences
    for (k = l = 0; k < raw->u.n+1; k++) {
        m = idx->cc.iug_idx[k];
        idx->cc.iug_idx[k] = l;
        l += m;
        if(m > 0) idx->cc.iug_b[l-1] = 0;
    }
    assert(l == iug_occ);

    for (k = 0; k < i_ug->u.n; k++) {
        x = &(idx->cc.iug_a[k]);
        for (z = 0; z < x->n; z++) {
            a = idx->cc.iug_b + idx->cc.iug_idx[x->a[z].v>>1];
            a_n = idx->cc.iug_idx[(x->a[z].v>>1)+1] - idx->cc.iug_idx[x->a[z].v>>1];
            if(a_n) {
                if(a[a_n-1] == a_n-1) a[a_n-1] = (k<<32)|z;///id of unitig | offset within the unitig
                else a[a[a_n-1]++] = (k<<32)|z;
            }        
        }
    }
}

static void worker_renew_integer_bridge(void *data, long i, int tid) // callback for kt_for()
{
    ul_resolve_t *uidx = (ul_resolve_t *)data; bubble_type *bub = uidx->bub;
    asg_t *bg = uidx->uovl.bg.bg; ul_str_idx_t *str_idx = &(uidx->pstr);
    asg_arc_t *p = &(bg->arc[i]); uint64_t *hid_a, hid_n, z; ul_str_t *str;
    uint32_t v = p->ul>>32, w = p->v, vz, wz; int64_t s, s_n, occ = 0;

    hid_a = str_idx->occ.a + str_idx->idx.a[v>>1];
    hid_n = str_idx->idx.a[(v>>1)+1] - str_idx->idx.a[v>>1];
    for (z = 0; z < hid_n; z++) {
        str = &(str_idx->str.a[hid_a[z]>>32]); s_n = str->cn;
        if(s_n < 2) continue;
        vz = (uint32_t)(str->a[(uint32_t)hid_a[z]]);
        assert((v>>1) == (vz>>1));
        s = (uint32_t)hid_a[z]; s -= 1; vz = ((uint32_t)(str->a[(uint32_t)hid_a[z]]))^1;
        for (; s >= 0; s--) {
            wz = (uint32_t)(str->a[s]); wz ^= 1;
            if(IF_HOM((wz>>1), *bub)) continue;
            if(vz == v && wz == w) occ++;
            else if(vz == (v^1) && wz == (w^1)) occ++;
            break;
        }

        s = (uint32_t)hid_a[z]; s += 1; vz = (uint32_t)(str->a[(uint32_t)hid_a[z]]);
        for (; s < s_n; s++) {
            wz = (uint32_t)(str->a[s]);
            if(IF_HOM((wz>>1), *bub)) continue;
            if(vz == v && wz == w) occ++;
            else if(vz == (v^1) && wz == (w^1)) occ++;
            break;
        }
    }
    p->ol = occ; p->ul >>= 32; p->ul <<= 32; p->ul += (((uint32_t)-1) - p->ol);
}

void renew_u2g_bg(ul_resolve_t *uidx)
{
    ul2ul_idx_t *idx = &(uidx->uovl); ul_bg_t *bg = &(idx->bg);
    ma_ug_t *iug = idx->i_ug, *raw = uidx->l1_ug; bubble_type *bub = uidx->bub;
    uint64_t k, z, l, *raw_a, raw_n, raw_id, iug_id, iug_off, v, w, nv, nw, n_vtx; 
    uinfo_srt_warp_t *x; asg64_v buf = {0,0,0}; int64_t s, s_n; asg_arc_t *p, *av, *aw;
    asg_destroy(bg->bg); free(bg->w_n); free(bg->a_n); memset(bg, 0, sizeof((*bg)));

    bg->bg = asg_init();
    bg->bg->n_seq = 0; bg->bg->m_seq = raw->g->n_seq; MALLOC(bg->bg->seq, bg->bg->m_seq);
    for (k = 0; k < raw->g->n_seq; k++) {
        raw_id = k; 
        if(IF_HOM(raw_id, *bub)) {///updated-line
            // asg_seq_set(bg->bg, k, raw->g->seq[k].len, 1);
            // bg->bg->seq[k].c = 0;
            continue;
        }
        raw_a = idx->cc.iug_b + idx->cc.iug_idx[raw_id];
        raw_n = idx->cc.iug_idx[raw_id+1] - idx->cc.iug_idx[raw_id];
        for (z = 0, buf.n = 0; z < raw_n; z++) {
            iug_id = raw_a[z]>>32; iug_off = (uint32_t)raw_a[z];
            x = &(idx->cc.iug_a[iug_id]); s_n = x->n;
            assert(raw_id == (x->a[iug_off].v>>1));
            if(iug->g->seq[iug_id].del) continue;
            for (s = ((int64_t)iug_off) - 1, v = x->a[iug_off].v^1; s >= 0; s--) {
                if(IF_HOM((x->a[s].v>>1), *bub)) continue;
                kv_push(uint64_t, buf, ((v<<32)|((uint64_t)(x->a[s].v^1))));
                break;
            }
            for (s = ((int64_t)iug_off) + 1, v = x->a[iug_off].v; s < s_n; s++) {
                if(IF_HOM((x->a[s].v>>1), *bub)) continue;
                kv_push(uint64_t, buf, ((v<<32)|((uint64_t)(x->a[s].v))));
                break;
            }
        }

        asg_seq_set(bg->bg, k, raw->g->seq[k].len, ((buf.n>0)?0:1));
        bg->bg->seq[k].c = 0;

        radix_sort_srt64(buf.a, buf.a + buf.n);
        for (l = 0, z = 1; z <= buf.n; z++) {
            if(z == buf.n || buf.a[l] != buf.a[z]) {
                p = asg_arc_pushp(bg->bg); memset(p, 0, sizeof(*p));
                p->v = (uint32_t)buf.a[l]; p->ol = z - l;
                p->ul = buf.a[l]>>32; p->ul <<= 32; 
                p->ul += (((uint32_t)-1) - p->ol);
                l = z;
            }
        }
    }

    kt_for(uidx->str_b.n_thread, worker_renew_integer_bridge, uidx, bg->bg->n_arc);
    asg_cleanup(bg->bg); bg->bg->r_seq = bg->bg->n_seq;

    /***********debug***********/
    for (k = 0; k < bg->bg->n_arc; k++) {
        p = &(bg->bg->arc[k]); v = p->v^1; w = (p->ul>>32)^1;
        av = asg_arc_a(bg->bg, v); nv = asg_arc_n(bg->bg, v);
        for (z = 0; z < nv; z++) {
            if(av[z].v == w) break;
        }
        assert(z < nv && p->ol == av[z].ol);
    }
    /***********debug***********/

    uint64_t v_occ[2], w_occ[2]; double ss = 0.500001;
    for (v = 0, n_vtx = bg->bg->n_seq<<1; v < n_vtx; v++) {
        if(bg->bg->seq[v>>1].del) continue;
        av = asg_arc_a(bg->bg, v); nv = asg_arc_n(bg->bg, v);
        if(!nv) continue; w = av[0].v;
        v_occ[0] = v_occ[1] = w_occ[0] = w_occ[1] = (uint32_t)-1; 

        av = asg_arc_a(bg->bg, v); nv = asg_arc_n(bg->bg, v);
        for (k = 0; k < nv && k < 2; k++) {
            v_occ[k] = av[k].ol; v_occ[k] <<= 32; v_occ[k] += av[k].v;
        }
    
        aw = asg_arc_a(bg->bg, (w^1)); nw = asg_arc_n(bg->bg, (w^1));
        for (k = 0; k < nw && k < 2; k++) {
            w_occ[k] = aw[k].ol; w_occ[k] <<= 32; w_occ[k] += aw[k].v;
        }

        if((((uint32_t)v_occ[0]) == w) && (((uint32_t)w_occ[0]) == (v^1))) {
            if(((((v_occ[0]>>32)*ss) >= (v_occ[1]>>32)) && (((w_occ[0]>>32)*ss) >= (w_occ[1]>>32))) || 
                                ((((v_occ[0]>>32)+(w_occ[0]>>32))*ss) >= ((v_occ[1]>>32)+(w_occ[1]>>32)))) {
                for (k = 0; k < nv; k++) av[k].ou = 2;///wrong 
                for (k = 0; k < nw; k++) aw[k].ou = 2;///wrong 
                av[0].ou = aw[0].ou = 3;//correct
            } else {
                for (k = 0; k < nv; k++) av[k].ou = 1;///ambg
                for (k = 0; k < nw; k++) aw[k].ou = 1;///ambg
            }
        }
    }

    for (k = 0; k < bg->bg->n_arc; k++) {
        p = &(bg->bg->arc[k]); v = p->v^1; w = (p->ul>>32)^1;
        av = asg_arc_a(bg->bg, v); nv = asg_arc_n(bg->bg, v);
        for (z = 0; z < nv; z++) {
            if(av[z].v == w) break;
        }
        assert(z < nv && p->ol == av[z].ol);
        l = MAX(p->ou, av[z].ou); if(l == 0) l = 1;
        p->ou = av[z].ou = l;
    }
    

    MALLOC(bg->w_n, iug->u.n); MALLOC(bg->a_n, iug->u.n);
    for (k = 0; k < iug->u.n; k++) {
        x = &(idx->cc.iug_a[k]); bg->w_n[k] = bg->a_n[k] = 0;
        for (z = nv = 0; z < x->n; z++) {
            if(IF_HOM((x->a[z].v>>1), *bub)) continue;
            v_occ[nv&1] = x->a[z].v; nv++;
            if(nv < 2) continue;
            l = get_bg_flag(uidx, v_occ[(nv-2)&1], v_occ[(nv-1)&1]);
            assert(l != bg_unavailable);
            if(l == bg_wrong) bg->w_n[k]++;
            if(l == bg_ambiguous) bg->a_n[k]++;
        }
        // fprintf(stderr, "-[M::%s::] k::%lu, x->n::%u, w_n[k]::%u, a_n[k]::%u\n", 
        //                 __func__, k, (uint32_t)x->n, bg->w_n[k], bg->a_n[k]);
    }
    
    kv_destroy(buf);
}

/**
static void worker_update_ul_tra_idx(void *data, long i, int tid) // callback for kt_for()
{
    ul_resolve_t *uidx = (ul_resolve_t *)data; integer_t *buf = &(uidx->str_b.buf[tid]);
    ma_ug_t *iug = uidx->uovl.i_ug; ul_tra_idx_t *iug_tra = &(uidx->uovl.iug_tra);
    uint32_t v = i, n_tra = iug_tra_arc_n(iug_tra, v), nv, k; 
    if(n_tra == 0) return;
    ul_tra_t *a_tra = iug_tra_arc_a(iug_tra, v); asg_arc_t *av;

    av = asg_arc_a(iug->g, v); nv = asg_arc_n(iug->g, v); kv_resize(uint64_t, buf->u, 2);
    for (k = 0; k < nv; k++) {
        buf->u.n = 0; buf->u.a[buf->u.n++] = v; buf->u.a[buf->u.n++] = av[k].v;
    }
}

void update_ul_tra_idx_t(ul_resolve_t *uidx)
{
    ma_ug_t *iug = uidx->uovl.i_ug; ul_tra_idx_t *iug_tra = &(uidx->uovl.iug_tra);
    uint64_t v, n_vtx = iug->g->n_seq<<1, l, m; 
    iug_tra->arc.n = iug_tra->idx.n = 0;
    kv_resize(uint32_t, iug_tra->idx, n_vtx + 1); iug_tra->idx.n = n_vtx + 1;
    memset(iug_tra->idx.a, 0, iug_tra->idx.n *sizeof(*(iug_tra->idx.a)));
    for (v = l = 0; v < n_vtx; v++) {
        m = asg_arc_n(iug->g, v);
        if(m < 2) m = 0;
        iug_tra->idx.a[v] = l;
        l += m;
    }
    iug_tra->idx.a[v] = l;
    kv_resize(ul_tra_t, iug_tra->arc, l); iug_tra->arc.n = l;
    
    kt_for(uidx->str_b.n_thread, worker_update_ul_tra_idx, uidx, n_vtx);///all ul + ug
}
**/
void renew_ul2_utg(ul_resolve_t *uidx)
{
    ul2ul_idx_t *z = &(uidx->uovl);
    if(uidx->uovl.cc.iug_a) {
        uint32_t k;
        for (k = 0; k < z->i_ug->u.n; k++) free(uidx->uovl.cc.iug_a[k].a);
    }
    renew_utg(&(z->i_ug), z->i_g, NULL); 
    free(z->i_ug->g->seq_vis); CALLOC(z->i_ug->g->seq_vis, z->i_ug->g->n_seq*2);
    renew_u2g_cov(uidx, uidx->uopt->te); 
    renew_u2g_bg(uidx);
}


ul2ul_idx_t *gen_ul2ul(ul_resolve_t *uidx, ug_opt_t *uopt, ulg_opt_t *ulopt, uint32_t keep_raw_utg, uint64_t is_bridg)
{
    uint64_t k, m;
    ma_ug_t *ug = uidx->l1_ug; all_ul_t *uls = uidx->idx; 
    ul2ul_idx_t *z = &(uidx->uovl); ul_str_idx_t *str_idx = &(uidx->pstr);
    z->uln = uls->n; z->gn = ug->g->n_seq; 
    z->tot = uls->n + ug->g->n_seq; MALLOC(z->item_idx, z->tot);
    for (k = m = 0; k < z->tot; k++) {
        z->item_idx[k] = (uint32_t)-1;
        if(k < z->uln) {///is a ul read
            if(str_idx->str.a[k].cn > 1) {
                z->item_idx[k] = m; m++;
            }
        } else {//is a node in graph
            z->item_idx[k] = m; m++;
        }
    }
    CALLOC(z->a, m); z->n = z->m = m;

    kt_for(uidx->str_b.n_thread, worker_integer_postprecess, uidx, uls->n);
    gen_integer_normalize(uidx);

    chimeric_integer_deal(uidx);
    append_utg_es(uidx);


    kt_for(uidx->str_b.n_thread, worker_integert_debug_sym, uidx, z->tot);///all ul + ug
    print_integert_ovlp_stat(z);
    remove_integert_containment(uidx, keep_raw_utg);

    kt_for(uidx->str_b.n_thread, worker_integert_debug_sym, uidx, z->tot);///all ul + ug
    print_integert_ovlp_stat(z);
    // print_uls_seq(uidx, asm_opt.output_file_name);
    // print_uls_ovs(uidx, asm_opt.output_file_name);
    
    z->i_g = integer_sg_gen(uidx, uopt->min_ovlp);
    asg_arc_del_trans(z->i_g, uopt->gap_fuzz);

    z->i_ug = ma_ug_gen(z->i_g);
    // CALLOC(z->i_ug->g->seq_vis, z->i_ug->g->n_seq*2);
    renew_u2g_cov(uidx, uidx->uopt->te); 
    renew_u2g_bg(uidx);
    
    // output_integer_graph(uidx, z->i_ug, asm_opt.output_file_name);
    u2g_clean(uidx, ulopt, keep_raw_utg, is_bridg);
    // renew_ul2_utg(uidx);

    // output_integer_graph(uidx, z->i_ug, asm_opt.output_file_name, 0);
    return z;
}

uint64_t str_occ_w(ul_str_t *str, ul_vec_t *raw, ma_ug_t *ug)
{
    uint64_t occ, k; uc_block_t *z;
    for (k = occ = 0; k < str->cn; k++) {
        z = &(raw->bb.a[str->a[k]>>32]);
        assert(((z->hid<<1)+z->rev)==((uint32_t)str->a[k]));
        occ += ug_occ_w(z->ts, z->te, &(ug->u.a[z->hid]));
    }
    return occ;
}

void gen_cul_g_t(ul_resolve_t *uidx)
{
    ma_ug_t *ug = uidx->l1_ug; uint64_t k; ul_str_idx_t *str_idx = &(uidx->pstr);
    CALLOC(uidx->cg, 1);
    uidx->cg->n[0] = uidx->idx->n; uidx->cg->n[1] = ug->g->n_seq; 
    uidx->cg->tot = uidx->cg->n[0] + uidx->cg->n[1];
    uidx->cg->g = asg_init();
    for (k = 0; k < uidx->cg->tot; k++) {
        if(k < uidx->cg->n[0]) {
            asg_seq_set(uidx->cg->g, k, str_occ_w(&(str_idx->str.a[k]), &(uidx->idx->a[k]), ug), 
                                                                        ((str_idx->str.a[k].cn>1)?0:1));
        } else {
            asg_seq_set(uidx->cg->g, k, ug->u.a[k-uidx->cg->n[0]].n, 0);
        }
    }
}

void init_ulg_opt_t(ulg_opt_t *z, ug_opt_t *uopt, int64_t clean_round, 
double min_path_drop_ratio, double max_path_drop_ratio,
double min_ovlp_drop_ratio, double max_ovlp_drop_ratio, double hom_check_drop_rate, 
int64_t max_tip, int64_t max_tip_hifi, bub_label_t *b_mask_t, uint32_t is_trio)
{
    z->tipsLen = uopt->tipsLen;
    z->tip_drop_ratio = uopt->tip_drop_ratio;
    z->stops_threshold = uopt->stops_threshold;
    z->chimeric_rate = uopt->chimeric_rate;
    z->drop_ratio = uopt->drop_ratio;


    z->b_mask_t = b_mask_t;
    z->clean_round = clean_round;
    z->min_path_drop_ratio = min_path_drop_ratio;
    z->max_path_drop_ratio = max_path_drop_ratio;
    z->min_ovlp_drop_ratio = min_ovlp_drop_ratio;
    z->max_ovlp_drop_ratio = max_ovlp_drop_ratio;
    z->hom_check_drop_rate = hom_check_drop_rate;
    z->max_tip = max_tip;
    z->max_tip_hifi = max_tip_hifi;
    z->is_trio = is_trio;
}

ma_ug_t* output_trio_unitig_graph_ul(ug_opt_t *uopt, ul_resolve_t *uidx, char* ou, uint8_t flag)
{
    char* gfa_name; MALLOC(gfa_name, strlen(ou)+100);
    sprintf(gfa_name, "%s.%s.p_ctg.gfa", ou, (flag==FATHER?"hap1":"hap2"));
    FILE* output_file = fopen(gfa_name, "w");

    ma_ug_t *ug = copy_untig_graph(uidx->uovl.hybrid_ug); 
    kvec_asg_arc_t_warp ne; kv_init(ne.a);

    adjust_utg_by_trio(&ug, uidx->sg, flag, TRIO_THRES, uopt->sources, uopt->reverse_sources, 
    uopt->coverage_cut, uopt->tipsLen, uopt->tip_drop_ratio, uopt->stops_threshold, uopt->ruIndex, 
    uopt->chimeric_rate, uopt->drop_ratio, uopt->max_hang, uopt->min_ovlp, uopt->gap_fuzz, &ne, uopt->b_mask_t);    

    // if(asm_opt.b_low_cov > 0) {
    //     break_ug_contig(&ug, uidx->sg, &R_INF, uopt->coverage_cut, uopt->sources, uopt->ruIndex, &ne, 
    //     uopt->max_hang, uopt->min_ovlp, &asm_opt.b_low_cov, NULL, asm_opt.m_rate);
    // }

    // if(asm_opt.b_high_cov > 0)
    // {
    //     break_ug_contig(&ug, uidx->sg, &R_INF, uopt->coverage_cut, uopt->sources, uopt->ruIndex, &ne, 
    //     uopt->max_hang, uopt->min_ovlp, NULL, &asm_opt.b_high_cov, asm_opt.m_rate);
    // }

    fprintf(stderr, "Writing %s to disk... \n", gfa_name);
    ma_ug_seq(ug, uidx->sg, uopt->coverage_cut, uopt->sources, &ne, uopt->max_hang, uopt->min_ovlp, 0, 1);
    
    ma_ug_print(ug, uidx->sg, uopt->coverage_cut, uopt->sources, uopt->ruIndex, (flag==FATHER?"h1tg":"h2tg"), output_file);
    fclose(output_file);

    sprintf(gfa_name, "%s.%s.p_ctg.noseq.gfa", ou, (flag==FATHER?"hap1":"hap2"));
    output_file = fopen(gfa_name, "w");
    ma_ug_print_simple(ug, uidx->sg, uopt->coverage_cut, uopt->sources, uopt->ruIndex, (flag==FATHER?"h1tg":"h2tg"), output_file);
    fclose(output_file);
    // if(asm_opt.bed_inconsist_rate != 0)
    // {
    //     sprintf(gfa_name, "%s.%s.p_ctg.lowQ.bed", output_file_name, f_prefix?f_prefix:(flag==FATHER?"hap1":"hap2"));
    //     output_file = fopen(gfa_name, "w");
    //     ma_ug_print_bed(ug, sg, &R_INF, coverage_cut, sources, &new_rtg_edges, 
    //     max_hang, min_ovlp, asm_opt.bed_inconsist_rate, (flag==FATHER?"h1tg":"h2tg"), output_file, NULL);
    //     fclose(output_file);
    // }

    free(gfa_name);
    ma_ug_destroy(ug);
    kv_destroy(ne.a);
    return NULL;
}



void realloc_rdb(All_reads* rdb, ma_sub_t **cov, R_to_U *ruI, uint64_t *rmap, uint64_t rid_n, uint64_t scaf_len, char *scaf_id, 
asg_t *ng, ug_opt_t *uopt)
{
    uint64_t i, rid_n0 = rdb->total_reads, tname, cname;
    uint64_t scaf_id_len = strlen(scaf_id); char *des, *src;
    rdb->total_reads = rid_n;
    tname = rdb->name_index[rid_n0];
    fprintf(stderr, "+[M::%s] rid_n0::%lu, rid_n::%lu\n", __func__, rid_n0, rid_n);
    ///for read bases
    REALLOC(rdb->N_site, rdb->total_reads);
    REALLOC(rdb->read_length, rdb->total_reads);
    REALLOC(rdb->read_size, rdb->total_reads);
    REALLOC(rdb->read_sperate, rdb->total_reads);
    REALLOC(rdb->trio_flag, rdb->total_reads);
    REALLOC(rdb->name_index, rdb->total_reads+1);///total_reads+1
    REALLOC((*cov), rdb->total_reads);
	for (i = rid_n0; i < rdb->total_reads; i++) {
        // fprintf(stderr, "[M::%s] i::%lu\n", __func__, i);
        rdb->N_site[i] = NULL;
        rdb->read_length[i] = scaf_len;
        rdb->read_size[i] = scaf_len;
        rdb->trio_flag[i] = AMBIGU;
        (*cov)[i].c = (*cov)[i].del = 0; 
        (*cov)[i].s = 0; (*cov)[i].e = scaf_len;
        cname = scaf_id_len;
        if(rmap[i] != ((uint64_t)-1)) {///not a scaffold node
            if(rdb->N_site[rmap[i]] != NULL) {
                MALLOC(rdb->N_site[i], rdb->N_site[rmap[i]][0]+1);
                memcpy(rdb->N_site[i], rdb->N_site[rmap[i]], 
                    sizeof((*(rdb->N_site[i])))*(rdb->N_site[rmap[i]][0]+1));
            }
            rdb->read_length[i] = rdb->read_length[rmap[i]];
            rdb->read_size[i] = rdb->read_length[rmap[i]];
            rdb->trio_flag[i] = rdb->trio_flag[rmap[i]];
            (*cov)[i] = (*cov)[rmap[i]];
            cname = Get_NAME_LENGTH((*rdb), (rmap[i]));
        } 
        rdb->name_index[i] = tname; tname += cname;
        rdb->total_reads_bases += rdb->read_length[i];

        MALLOC(rdb->read_sperate[i], (rdb->read_length[i]/4+1));
        if(rmap[i] != ((uint64_t)-1)) {///not a scaffold node
            memcpy(rdb->read_sperate[i], rdb->read_sperate[rmap[i]], 
                sizeof((*(rdb->read_sperate[i])))*(rdb->read_length[i]/4+1));
        } else {
            ///set to A
            memset(rdb->read_sperate[i], 0, sizeof((*(rdb->read_sperate[i])))*(rdb->read_length[i]/4+1));
        }
	}

    rdb->index_size = rdb->total_reads;
    rdb->name_index[i] = tname; 
    rdb->total_name_length = tname; 
    rdb->name_index_size = rdb->total_reads+1;
    REALLOC(rdb->name, tname);
    for (i = rid_n0; i < rdb->total_reads; i++) {
        des = Get_NAME((*rdb), i); src = scaf_id; cname = scaf_id_len;
        if(rmap[i] != ((uint64_t)-1)) {
            src = Get_NAME((*rdb), rmap[i]); cname = Get_NAME_LENGTH((*rdb), (rmap[i]));
        }
        memcpy(des, src, sizeof((*(des)))*cname);
    }

    REALLOC(rdb->paf, rdb->total_reads); 
    memset(rdb->paf+rid_n0, 0, (rdb->total_reads-rid_n0)*sizeof((*rdb->paf)));
    REALLOC(rdb->reverse_paf, rdb->total_reads);
    memset(rdb->reverse_paf+rid_n0, 0, (rdb->total_reads-rid_n0)*sizeof((*rdb->reverse_paf)));

    ruI->len = rdb->total_reads;
    REALLOC(ruI->index, ruI->len); 
    memset(ruI->index, -1, sizeof((*(ruI->index)))*(ruI->len));

    reset_bub_label_t(uopt->b_mask_t, ng, 0, 0);
    uopt->coverage_cut = (*cov);
    uopt->reverse_sources = rdb->reverse_paf;
    uopt->sources = rdb->paf;
    fprintf(stderr, "-[M::%s] rid_n0::%lu, rid_n::%lu\n", __func__, rid_n0, rid_n);
}

inline void update_qtn(ma_hit_t *z, uint64_t qn, uint64_t tn)
{
    z->qns <<= 32; z->qns >>= 32; z->qns |= (qn<<32); z->tn = tn;
}

inline uint64_t dup_paf_check(ma_hit_t_alloc *x, ma_hit_t *p)
{
    int64_t i; ma_hit_t *z;
    for (i = 0; i < x->length; i++) {
        z = &(x->buffer[i]);
        if((z->qns == p->qns) && (z->tn == p->tn) && (z->qe == p->qe) && (z->ts == p->ts) && (z->te == p->te) && 
            (z->cc == p->cc) && (z->ml == p->ml) && (z->rev == p->rev) && (z->bl == p->bl) && (z->del == p->del) && 
            (z->el == p->el) && (z->no_l_indel == p->no_l_indel)) {
            return 0;
        }    
    }
    return 1;
}

void renew_paf0(ma_hit_t_alloc *paf, uint64_t *a0, uint64_t a0n, uint64_t *a1, uint64_t a1n, asg64_v *buf)
{
    if(a0n <= 1 && a1n <= 1) return;
    uint64_t qn = ((uint32_t)a0[0]), tn = ((uint32_t)a1[0]), i, k, sf = 0;
    ma_hit_t e01, e10; uint64_t qi, ti, *qa, *ta, qlen, tlen; buf->n = 0; ma_hit_t_alloc *qo, *to;
    if(qn == tn) sf = 1;
    qo = &(paf[qn]); qi = qo->length;
    for (i = qlen = 0, qa = NULL; i < qi; i++) {
        if(qo->buffer[i].tn == tn) kv_push(uint64_t, *buf, i);
    }
    qlen = buf->n; 

    to = &(paf[tn]); ti = to->length;
    for (i = tlen = 0, ta = NULL; i < ti; i++) {
        if(to->buffer[i].tn == qn) kv_push(uint64_t, *buf, i);
    }
    tlen = buf->n - qlen;
    assert(qlen && tlen);
    qa = buf->a; ta = buf->a + qlen;

    for (qi = 0; qi < qlen; qi++) {
        e01 = qo->buffer[qa[qi]];
        for (ti = 0; ti < tlen; ti++) {
            e10 = to->buffer[ta[ti]];

            for (i = 0; i < a0n; i++) {
                qn = ((uint32_t)a0[i]);
                for (k = 0; k < a1n; k++) {
                    if(i == 0 && k == 0) continue;
                    tn = ((uint32_t)a1[k]);

                    update_qtn(&e01, qn, tn); 
                    if((!sf) || (dup_paf_check(&(paf[qn]), &e01))) {
                        add_ma_hit_t_alloc(&(paf[qn]), &e01);
                    }

                    update_qtn(&e10, tn, qn); 
                    if((!sf) || (dup_paf_check(&(paf[tn]), &e10))) {
                        add_ma_hit_t_alloc(&(paf[tn]), &e10);
                    }
                }
            }
        }
    }
    
    /**
    idx = get_specific_overlap(&(paf[qn]), qn, tn);
    if(idx < 0) return;
    e01 = paf[qn].buffer[idx];
    idx = get_specific_overlap(&(paf[tn]), tn, qn);
    e10 = paf[tn].buffer[idx];

    for (i = 0; i < a0n; i++) {
        qn = ((uint32_t)a0[i]);
        for (k = 0; k < a1n; k++) {
            if(i == 0 && k == 0) continue;
            tn = ((uint32_t)a1[k]);
            if((qn == 5619628 && tn == 5619629) || (tn == 5619628 && qn == 5619629)) {
                fprintf(stderr, "[M::%s]\tqn::%lu(qg::%u)\ttn::%lu(tg::%u)\n", __func__, 
                qn, ((uint32_t)a0[0]), tn, ((uint32_t)a1[0]));
            }
            update_qtn(&e01, qn, tn); add_ma_hit_t_alloc(&(paf[qn]), &e01);
            update_qtn(&e10, tn, qn); add_ma_hit_t_alloc(&(paf[tn]), &e10);
        }
    }
    **/
}

void renew_paf1(ma_hit_t_alloc *paf, uint64_t *a, uint64_t an, uint64_t len)
{
    uint64_t i, k, qn, tn; ma_hit_t arc;
    arc.qns = 0; arc.qe = len;
    arc.tn = 0; arc.ts = 0; arc.te = len;
    arc.rev = arc.el = arc.ml = arc.no_l_indel = arc.bl = 0;

    for (i = 0; i < an; i++) {
        qn = ((uint32_t)a[i]);
        for (k = 0; k < an; k++) {
            tn = ((uint32_t)a[k]);
            update_qtn(&arc, qn, tn); 
            add_ma_hit_t_alloc(&(paf[qn]), &arc);
        }
    }
}

void update_paf(ma_hit_t_alloc *src, ma_hit_t_alloc *r_src, uint64_t *rmap, uint64_t rid_n, uint64_t pre_gn, asg_t *ng)
{
    asg64_v clus, bb; uint64_t k, l, dn = 0, j, qn, tn; 
    uint64_t *idx; ma_hit_t_alloc *z; 
    uint64_t *a0, *a1, a0n, a1n;
    CALLOC(idx, pre_gn); kv_init(bb);
    kv_init(clus); kv_resize(uint64_t, clus, rid_n);
    for (k = 0; k < rid_n; k++) {
        if(rmap[k] == ((uint64_t)-1)) continue;///scaffold
        clus.a[clus.n++] = (rmap[k]<<32)|k;
    }
    radix_sort_srt64(clus.a, clus.a+clus.n);
    for (k = 1, l = 0; k <= clus.n; k++) {
       if(k == clus.n || (clus.a[k]>>32) != (clus.a[l]>>32)) {
            idx[(clus.a[l]>>32)] = (l<<32)|k;
            l = k; dn++;
       }
    }
    assert(dn == pre_gn);
    fprintf(stderr, "+[M::%s] dn::%lu\n", __func__, dn);

    for (k = 0; k < dn; k++) {
        a0 = clus.a + (idx[k]>>32); a0n = ((uint32_t)idx[k]) - (idx[k]>>32);
        assert(a0n > 0 && ((uint32_t)a0[0]) == k);
        // if(k == 4265018) {
        //     fprintf(stderr, "\n[M::%s]\tgroup::%lu\n", __func__, k);
        //     for (j = 0; j < a0n; j++) {
        //         fprintf(stderr, "[M::%s]\tnid::%uu\n", __func__, (uint32_t)a0[j]);
        //     }
        // }
        z = &(src[k]);
        for (j = 0; j < z->length; j++) {
            qn = Get_qn(z->buffer[j]);
            tn = Get_tn(z->buffer[j]);
            if(tn >= dn) continue;
            if(qn > tn) continue;
            a1 = clus.a + (idx[tn]>>32); 
            a1n = ((uint32_t)idx[tn]) - (idx[tn]>>32);
            // fprintf(stderr, "+[M::%s] qn::%lu, tn::%lu, a0n::%lu, a1n::%lu\n", __func__, qn, tn, a0n, a1n);
            assert(a1n > 0 && ((uint32_t)a1[0]) == tn);
            renew_paf0(src, a0, a0n, a1, a1n, &bb);
        }

        z = &(r_src[k]);
        for (j = 0; j < z->length; j++) {
            qn = Get_qn(z->buffer[j]);
            tn = Get_tn(z->buffer[j]);
            if(tn >= dn) continue;
            if(qn > tn) continue;
            a1 = clus.a + (idx[tn]>>32); 
            a1n = ((uint32_t)idx[tn]) - (idx[tn]>>32);
            // fprintf(stderr, "-[M::%s] qn::%lu, tn::%lu, a0n::%lu, a1n::%lu\n", __func__, qn, tn, a0n, a1n);
            assert(a1n > 0 && ((uint32_t)a1[0]) == tn);
            renew_paf0(r_src, a0, a0n, a1, a1n, &bb);
        }



        // for (l = k; l < dn; l++) {
        //     a1 = clus.a + (idx[l]>>32); a1n = ((uint32_t)idx[l]) - (idx[l]>>32);
        //     fprintf(stderr, "[M::%s] k::%lu, l::%lu, a0n::%lu, a1n::%lu\n", __func__, k, l, a0n, a1n);
        //     assert(a1n > 0 && ((uint32_t)a1[0]) == l);
        //     renew_paf0(src, a0, a0n, a1, a1n);
        //     renew_paf0(r_src, a0, a0n, a1, a1n);
        // }
        if(a0n > 1) renew_paf1(r_src, a0, a0n, ng->seq[k].len);
    }
    free(idx); free(clus.a); free(bb.a);
}

void push_scaff_node(ma_hit_t_alloc *src, uint64_t v, uint64_t w, uint64_t ol, asg_t *ng)
{
    uint64_t vl = ng->seq[v>>1].len, wl = ng->seq[w>>1].len; 
    ma_hit_t arc, arc1;
    arc.qns = (v>>1)<<32; 
    if(!(v&1)) {
        arc.qns += vl - ol; arc.qe = vl;
    } else {
        arc.qns += 0; arc.qe = ol;
    }
    arc.tn = w>>1;
    if(!(w&1)) {
        arc.ts = 0; arc.te = ol;
    } else {
        arc.ts = wl - ol; arc.te = wl;
    }

    arc.rev = (v^w)&1; arc.el = 0;
    arc.ml = arc.no_l_indel = arc.bl = 0;

    add_ma_hit_t_alloc(&(src[v>>1]), &arc);
    set_reverse_overlap(&arc1, &arc);
    add_ma_hit_t_alloc(&(src[w>>1]), &arc1);
}

asg_t *gen_ng(ma_ug_t *ug, asg_t *sg, ug_opt_t *uopt, ma_sub_t **cov, R_to_U *ruI, 
uint64_t scaffold_len)
{
    ma_utg_t *u; uint64_t i, v, w, m, h, z, raw_v, raw_w, nocc, nv, vx, wx; int32_t r;
    asg_arc_t t, *p; asg_arc_t *av = NULL; uint64_t slen = scaffold_len + (uopt->min_ovlp*2);
    asg_ext_t ext; memset(&ext, 0, sizeof(ext)); ext.ext = asg_init(); asg_t *ng = ext.ext;
    ext.a_n = sg->n_seq; 
    ext.cnt.n = ext.cnt.m = ext.a_n;
    CALLOC(ext.cnt.a, ext.cnt.n);///count

    ext.idx_a.n = ext.idx_a.m = sg->n_seq; 
    MALLOC(ext.idx_a.a, ext.idx_a.n); 
    memset(ext.idx_a.a, -1, sizeof(*(ext.idx_a.a))*ext.idx_a.n);//map
    fprintf(stderr, "\n+[M::%s] 0\n", __func__);

    for (i = 0; i < ug->g->n_seq; ++i) {
        u = &(ug->u.a[i]);
        for (m = 0; m < u->n; m++) {
            v = u->a[m]>>32;
            ext.cnt.a[v>>1]++;
            if(ext.cnt.a[v>>1] == 1) {
                h = v>>1; ext.idx_a.a[h] = v>>1;
            } else {
                h = ext.idx_a.n; 
                kv_push(uint64_t, ext.idx_a, (v>>1));
            }

            z = (h<<1)|(v&1); z <<= 32; z |= ((uint32_t)u->a[m]);
            u->a[m] = z;
        }
    }

    for (i = 0; i < sg->n_seq; ++i) {
        asg_seq_set(ng, i, sg->seq[i].len, ext.idx_a.a[i]==((uint64_t)-1)?1:0);
        ng->seq[i].c = 0;
        if(ext.idx_a.a[i]!=((uint64_t)-1)) assert(ext.idx_a.a[i] == i);
        ext.idx_a.a[i] = i;///set for delted read
    }
    for (; i < ext.idx_a.n; i++) {
        asg_seq_set(ng, i, sg->seq[ext.idx_a.a[i]].len, 0);
        ng->seq[i].c = 0; assert(!(ng->seq[ext.idx_a.a[i]].del));
    }
    ng->r_seq = ng->n_seq;
    // fprintf(stderr, "+[M::%s] ng->n_seq::%u, ext.idx_a.n::%u\n", 
    // __func__, (uint32_t)ng->n_seq, (uint32_t)ext.idx_a.n);
    assert(ng->n_seq == ext.idx_a.n); 
    fprintf(stderr, "\n+[M::%s] 1\n", __func__);

    for (i = 0, nocc = ng->r_seq; i < ug->g->n_seq; ++i) {
        u = &(ug->u.a[i]);
        for (m = 1; m < u->n; m++) {
            v = u->a[m-1]>>32; w = u->a[m]>>32;
            raw_v = (ext.idx_a.a[v>>1]<<1)|(v&1);
            raw_w = (ext.idx_a.a[w>>1]<<1)|(w&1);
            assert((raw_v>>1) < sg->n_seq);
            assert((raw_w>>1) < sg->n_seq);
            if(gen_spec_edge(sg, uopt, raw_v, raw_w, &t) < 0) {
                asg_seq_set(ng, nocc, slen, 0);
                ng->seq[nocc].c = 0;
                kv_push(uint64_t, ext.idx_a, ((uint64_t)-1));
                nocc++;///a scaffold node
            }
        }

        v = i<<1; nv = asg_arc_n(ug->g, v); av = asg_arc_a(ug->g, v);
        for (m = 0; m < nv; m++) {
            if(av[m].del) continue;
            w = av[m].v;
            vx = (v&1?((ug->u.a[v>>1].a[0]>>32)^1):(ug->u.a[v>>1].a[ug->u.a[v>>1].n-1]>>32));
            wx = (w&1?((ug->u.a[w>>1].a[ug->u.a[w>>1].n-1]>>32)^1):(ug->u.a[w>>1].a[0]>>32));
            raw_v = (ext.idx_a.a[vx>>1]<<1)|(vx&1);
            raw_w = (ext.idx_a.a[wx>>1]<<1)|(wx&1);
            assert((raw_v>>1) < sg->n_seq);
            assert((raw_w>>1) < sg->n_seq);
            if(vx > wx) continue;
            if((vx == wx) && (vx&1)) continue;
            if(gen_spec_edge(sg, uopt, raw_v, raw_w, &t) < 0) {
                asg_seq_set(ng, nocc, slen, 0);
                ng->seq[nocc].c = 0;
                kv_push(uint64_t, ext.idx_a, ((uint64_t)-1));
                nocc++;///a scaffold node
            }
        }

        v = (i<<1)+1; nv = asg_arc_n(ug->g, v); av = asg_arc_a(ug->g, v);
        for (m = 0; m < nv; m++) {
            if(av[m].del) continue;
            w = av[m].v;
            vx = (v&1?((ug->u.a[v>>1].a[0]>>32)^1):(ug->u.a[v>>1].a[ug->u.a[v>>1].n-1]>>32));
            wx = (w&1?((ug->u.a[w>>1].a[ug->u.a[w>>1].n-1]>>32)^1):(ug->u.a[w>>1].a[0]>>32));
            raw_v = (ext.idx_a.a[vx>>1]<<1)|(vx&1);
            raw_w = (ext.idx_a.a[wx>>1]<<1)|(wx&1);
            assert((raw_v>>1) < sg->n_seq);
            assert((raw_w>>1) < sg->n_seq);
            if(vx > wx) continue;
            if((vx == wx) && (vx&1)) continue;
            if(gen_spec_edge(sg, uopt, raw_v, raw_w, &t) < 0) {
                asg_seq_set(ng, nocc, slen, 0);
                ng->seq[nocc].c = 0;
                kv_push(uint64_t, ext.idx_a, ((uint64_t)-1));
                nocc++;///a scaffold node
            }
        }
    }
    fprintf(stderr, "\n+[M::%s] 2\n", __func__);
    // fprintf(stderr, "+[M::%s] ng->n_seq::%u, ext.idx_a.n::%u\n", 
    // __func__, (uint32_t)ng->n_seq, (uint32_t)ext.idx_a.n);
    assert(ng->n_seq == ext.idx_a.n); 
    CALLOC(ng->seq_vis, (ng->n_seq<<1));
    realloc_rdb(&(R_INF), cov, ruI, ext.idx_a.a, ext.idx_a.n, slen, (char *)"scaf", ng, uopt);

    fprintf(stderr, "\n+[M::%s] 3\n", __func__);
    update_paf(R_INF.paf, R_INF.reverse_paf, ext.idx_a.a, ext.idx_a.n, sg->n_seq, ng);
    fprintf(stderr, "\n+[M::%s] 4\n", __func__);

    for (i = 0, nocc = ng->r_seq; i < ug->g->n_seq; ++i) {
        u = &(ug->u.a[i]);
        for (m = 1; m < u->n; m++) {
            // fprintf(stderr, "[M::%s] i::%lu, m::%lu\n", __func__, i, m);
            v = u->a[m-1]>>32; w = u->a[m]>>32;
            r = gen_spec_edge(ng, uopt, v, w, &t);
            if(r >= 0) {
                p = asg_arc_pushp(ng); *p = t;
                r = gen_spec_edge(ng, uopt, w^1, v^1, &t);
                assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
            } else {
                z = nocc<<1; nocc++;
                push_scaff_node(R_INF.paf, v, z, uopt->min_ovlp, ng);
                push_scaff_node(R_INF.paf, z, w, uopt->min_ovlp, ng);

                r = gen_spec_edge(ng, uopt, v, z, &t);
                assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
                r = gen_spec_edge(ng, uopt, z^1, v^1, &t);
                assert(r >= 0); p = asg_arc_pushp(ng); *p = t;

                r = gen_spec_edge(ng, uopt, z, w, &t);
                assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
                r = gen_spec_edge(ng, uopt, w^1, z^1, &t);
                assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
            }
        }

        v = i<<1; nv = asg_arc_n(ug->g, v); av = asg_arc_a(ug->g, v);
        for (m = 0; m < nv; m++) {
            if(av[m].del) continue;
            w = av[m].v;
            vx = (v&1?((ug->u.a[v>>1].a[0]>>32)^1):(ug->u.a[v>>1].a[ug->u.a[v>>1].n-1]>>32));
            wx = (w&1?((ug->u.a[w>>1].a[ug->u.a[w>>1].n-1]>>32)^1):(ug->u.a[w>>1].a[0]>>32));
            if(vx > wx) continue;
            if((vx == wx) && (vx&1)) continue;///it is possible

            r = gen_spec_edge(ng, uopt, vx, wx, &t);
            if(r >= 0) {
                p = asg_arc_pushp(ng); *p = t;
                r = gen_spec_edge(ng, uopt, wx^1, vx^1, &t);
                assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
            } else {
                z = nocc<<1; nocc++;
                push_scaff_node(R_INF.paf, vx, z, uopt->min_ovlp, ng);
                push_scaff_node(R_INF.paf, z, wx, uopt->min_ovlp, ng);

                r = gen_spec_edge(ng, uopt, vx, z, &t);
                assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
                r = gen_spec_edge(ng, uopt, z^1, vx^1, &t);
                assert(r >= 0); p = asg_arc_pushp(ng); *p = t;

                r = gen_spec_edge(ng, uopt, z, wx, &t);
                assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
                r = gen_spec_edge(ng, uopt, wx^1, z^1, &t);
                assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
            }
        }

        v = (i<<1)+1; nv = asg_arc_n(ug->g, v); av = asg_arc_a(ug->g, v);
        for (m = 0; m < nv; m++) {
            if(av[m].del) continue;
            w = av[m].v;
            vx = (v&1?((ug->u.a[v>>1].a[0]>>32)^1):(ug->u.a[v>>1].a[ug->u.a[v>>1].n-1]>>32));
            wx = (w&1?((ug->u.a[w>>1].a[ug->u.a[w>>1].n-1]>>32)^1):(ug->u.a[w>>1].a[0]>>32));
            if(vx > wx) continue;
            if((vx == wx) && (vx&1)) continue;///it is possible

            r = gen_spec_edge(ng, uopt, vx, wx, &t);
            if(r >= 0) {
                p = asg_arc_pushp(ng); *p = t;
                r = gen_spec_edge(ng, uopt, wx^1, vx^1, &t);
                assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
            } else {
                z = nocc<<1; nocc++;
                push_scaff_node(R_INF.paf, vx, z, uopt->min_ovlp, ng);
                push_scaff_node(R_INF.paf, z, wx, uopt->min_ovlp, ng);

                r = gen_spec_edge(ng, uopt, vx, z, &t);
                assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
                r = gen_spec_edge(ng, uopt, z^1, vx^1, &t);
                assert(r >= 0); p = asg_arc_pushp(ng); *p = t;

                r = gen_spec_edge(ng, uopt, z, wx, &t);
                assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
                r = gen_spec_edge(ng, uopt, wx^1, z^1, &t);
                assert(r >= 0); p = asg_arc_pushp(ng); *p = t;
            }
        }
    }
    fprintf(stderr, "\n+[M::%s] 5\n", __func__);
    asg_cleanup(ng); ng->r_seq = ng->n_seq;
    free(ext.cnt.a); free(ext.idx_a.a); 
    fprintf(stderr, "[M::%s] nocc::%lu, ng->n_seq::%u, sg->n_seq::%u\n", 
    __func__, nocc, (uint32_t)ng->n_seq, (uint32_t)sg->n_seq);
    assert(nocc == ng->n_seq);
    return ng;
}

void gen_ul_trio_graph(ug_opt_t *uopt, ul_resolve_t *uidx, char *o_file)
{
    output_trio_unitig_graph_ul(uopt, uidx, o_file, FATHER);
    output_trio_unitig_graph_ul(uopt, uidx, o_file, MOTHER);
}

void destroy_integer_t(integer_t *z)
{
    free(z->f.a); free(z->p.a); free(z->o.a); free(z->u.a);
    free(z->vis.a); free(z->sc.a); free(z->snp.a); free(z->res_dump.a);
    free(z->q.a); free(z->t.a); free(z->b.a); 
    
    free(z->pg.seq.a); free(z->pg.arc.a); free(z->pg.idx.a); free(z->pg.e_idx.a);
    free(z->pg.srt_b.ind.a);
    free(z->pg.srt_b.stack.a);
    free(z->pg.srt_b.res.a);
    free(z->pg.srt_b.res2nid.a);
    free(z->pg.srt_b.aln.a);

    free(z->pg.bb.a.a);
    free(z->pg.bb.S.a);
    free(z->pg.bb.T.a);
    free(z->pg.bb.b.a);
    free(z->pg.bb.e.a);
    free(z->pg.bb.us.a);

    free(z->pg.bb.dp.ref.a);
    free(z->pg.bb.dp.pat.a);
    free(z->pg.bb.dp.pat_cor.a);
    free(z->pg.bb.dp.g_flt.a);
    free(z->pg.bb.dp.m_dir.a);
    free(z->pg.bb.dp.m_score.a);
}

void usg_t_destroy(usg_t *g) {
    uint32_t k;
    for (k = 0; k < g->mp.n; k++) {
        free(g->mp.a[k].a);
    }
    free(g->mp.a);

    for (k = 0; k < g->n; k++) {
        free(g->a[k].arc[0].a); free(g->a[k].arc[1].a);
        free(g->a[k].arc_mm[0].a); free(g->a[k].arc_mm[1].a);
    }
    free(g->a);
    
    free(g);
}


void destroy_ul_resolve_t(ul_resolve_t *uidx)
{
    ma_ug_destroy(uidx->l1_ug);
    uint32_t k;
    for (k = 0; k < uidx->str_b.n_thread; k++) {
        destroy_integer_t(&(uidx->str_b.buf[k]));
    }
    free(uidx->str_b.buf);

    for (k = 0; k < uidx->pstr.str.n; k++) {
        free(uidx->pstr.str.a[k].a);
    }
    free(uidx->pstr.str.a); 
    free(uidx->pstr.occ.a); 
    free(uidx->pstr.idx.a);

    
    free(uidx->uovl.telo);
    for (k = 0; k < uidx->uovl.n; k++) {
        free(uidx->uovl.a[k].a);
    }
    free(uidx->uovl.a);
    free(uidx->uovl.item_idx);
    asg_destroy(uidx->uovl.i_g);
    

    free(uidx->uovl.cc.uc); 
    free(uidx->uovl.cc.hc); 
    free(uidx->uovl.cc.raw_uc); 
    if(uidx->uovl.cc.iug_a) {
        for (k = 0; k < uidx->uovl.i_ug->u.n; k++) free(uidx->uovl.cc.iug_a[k].a);
    }
    free(uidx->uovl.cc.iug_a); 
    free(uidx->uovl.cc.iug_idx); 
    free(uidx->uovl.cc.iug_b);

    asg_destroy(uidx->uovl.bg.bg); 
    free(uidx->uovl.bg.w_n); 
    free(uidx->uovl.bg.a_n);
    ma_ug_destroy(uidx->uovl.i_ug);
    usg_t_destroy(uidx->uovl.h_usg);

    free(uidx);
}

static void clear_ma_hit_t_alloc(void *data, long i, int tid)
{
    ma_hit_t_alloc *src = (ma_hit_t_alloc *)data;
    ma_hit_t_alloc *z = &(src[i]); uint32_t k;
    for (k = 0; k < z->length; k++) z->buffer[k].del = 0;
}

static void reset_ma_sub_t(void *data, long i, int tid)
{
    sset_aux *s = (sset_aux *)data;
    if(s->g && s->g->seq[i].del) s->cov[i].del = 1;
    else s->cov[i].del = 0;
}

void renew_R_to_U(asg_t *ng, ma_hit_t_alloc* src, ma_hit_t_alloc* r_src, int64_t n_read, ma_sub_t *coverage_cut, 
R_to_U* ruIndex, int max_hang, int min_ovlp)
{
    sset_aux s; memset(&s, 0, sizeof(s)); 
    kt_for(asm_opt.thread_num, clear_ma_hit_t_alloc, src, n_read);
    kt_for(asm_opt.thread_num, clear_ma_hit_t_alloc, r_src, n_read);

    s.g = NULL; s.cov = coverage_cut;
    kt_for(asm_opt.thread_num, reset_ma_sub_t, &s, n_read);

    ma_hit_contained_advance(src, n_read, coverage_cut, ruIndex, max_hang, min_ovlp);

    // s.g = ng; s.cov = coverage_cut;
    // kt_for(asm_opt.thread_num, reset_ma_sub_t, &s, n_read);
}

void ul_realignment_gfa(ug_opt_t *uopt, asg_t *sg, int64_t clean_round, double min_ovlp_drop_ratio, 
double max_ovlp_drop_ratio, int64_t max_tip, int64_t max_ul_tip, bub_label_t *b_mask_t, uint32_t is_trio, char *o_file,
ul_renew_t *ropt, const char *bin_file, uint64_t free_uld, uint64_t is_bridg, uint64_t deep_clean)
{
    uint64_t i, bn = 0, idn = 0; uint32_t *bl = NULL; uint8_t *r_het = NULL; bubble_type *bub = NULL; ulg_opt_t uu;
    for (i = 0; i < sg->n_seq; ++i) {
        if(sg->seq[i].del) continue;
        sg->seq[i].c = PRIMARY_LABLE;
    }
    // fprintf(stderr, "0[M::%s]\n", __func__);
    hic_clean(sg);
    if(deep_clean) {
        // print_debug_gfa(sg, NULL, uopt->coverage_cut, "bclean", uopt->sources, uopt->ruIndex, uopt->max_hang, uopt->min_ovlp, 0, 0, 0);
        deep_graph_clean(uopt, sg, 1, is_trio, asm_opt.max_short_tip, asm_opt.min_drop_rate, 
        MIN(asm_opt.max_drop_rate, 0.7), 0.75, 1, asm_opt.clean_round, 10/**20**/);
        // print_debug_gfa(sg, NULL, uopt->coverage_cut, "aclean", uopt->sources, uopt->ruIndex, uopt->max_hang, uopt->min_ovlp, 0, 0, 0);
    }
    
    // fprintf(stderr, "1[M::%s]\n", __func__);
    ma_ug_t *init_ug = ul_realignment(uopt, sg, 0, bin_file);
    // fprintf(stderr, "2[M::%s]\n", __func__);
    // exit(1);

    // char* gfa_name = NULL; MALLOC(gfa_name, strlen(o_file)+strlen(bin_file)+50);
    // sprintf(gfa_name, "%s.%s", o_file, bin_file);
    // print_debug_gfa(sg, init_ug, uopt->coverage_cut, gfa_name, uopt->sources, uopt->ruIndex, uopt->max_hang, uopt->min_ovlp, 0, 0, 0);
    // print_debug_gfa(sg, init_ug, uopt->coverage_cut, gfa_name, uopt->sources, uopt->ruIndex, uopt->max_hang, uopt->min_ovlp, 0, 0, 1);
    // free(gfa_name);
    
    
    filter_sg_by_ug(sg, init_ug, uopt);
    // fprintf(stderr, "-0-[M::%s]\tUL_INF.a[25].rlen::%u\n", __func__, UL_INF.a[25].rlen);
    // print_debug_gfa(sg, init_ug, uopt->coverage_cut, "UL.debug", uopt->sources, uopt->ruIndex, uopt->max_hang, uopt->min_ovlp, 0, 0, 0);
    // print_ul_alignment(init_ug, &UL_INF, 47072, "after-0");
    bub = gen_bubble_chain(sg, init_ug, uopt, &r_het, ((asm_opt.polyploidy>2)?1:0));
    // fprintf(stderr, "4[M::%s]\n", __func__);
    // print_ul_alignment(init_ug, &UL_INF, 47072, "after-1");
    ul_resolve_t *uidx = init_ul_resolve_t(sg, init_ug, bub, &UL_INF, uopt, r_het);
    if(!free_uld) {///backup
        bn = UL_INF.n; idn = UL_INF.nid.n; MALLOC(bl, bn);
        for (i = 0; i < bn; i++) bl[i] = UL_INF.a[i].rlen;
    }
    // fprintf(stderr, "5[M::%s]\n", __func__);
    // print_ul_alignment(init_ug, &UL_INF, 47072, "after-2");
    // exit(1);
    // if(free_uld) {
        // print_raw_uls_seq(uidx, asm_opt.output_file_name);
        // print_raw_uls_aln(uidx, asm_opt.output_file_name);
    // }
    
    ul_re_correct(uidx, asm_opt.integer_correct_round/**3**/); 
    init_ulg_opt_t(&uu, uopt, clean_round, asm_opt.min_path_drop_rate, asm_opt.max_path_drop_rate, min_ovlp_drop_ratio, max_ovlp_drop_ratio, 0.55, max_tip, max_ul_tip, b_mask_t, is_trio);
    // print_debug_gfa(sg, init_ug, uopt->coverage_cut, "UL.debug0", uopt->sources, uopt->ruIndex, uopt->max_hang, uopt->min_ovlp, 0, 0, 1);
    /**ul2ul_idx_t *u2o = **/gen_ul2ul(uidx, uopt, &uu, 0, is_bridg);
    // print_ul_alignment(init_ug, &UL_INF, 47072, "after-3");

    // print_debug_ul("UL.debug", init_ug, sg, uopt->coverage_cut, uopt->sources, uopt->ruIndex, bub, &UL_INF);

    // resolve_dip_bub_chains(uidx);

    // free(r_het); destory_bubbles(bub); free(bub);
    // if(free_uld) {
        // uidx->uovl.hybrid_ug = gen_hybrid_ug(uidx, uidx->uovl.h_usg);
    //     // print_debug_gfa(sg, uidx->uovl.hybrid_ug, uopt->coverage_cut, "hybrid_ug", uopt->sources, uopt->ruIndex, uopt->max_hang, uopt->min_ovlp, 0, 0, 1);
        // print_debug_gfa(sg, uidx->uovl.hybrid_ug, uopt->coverage_cut, "hybrid_ug", uopt->sources, uopt->ruIndex, uopt->max_hang, uopt->min_ovlp, 0, 0, 0);
    //     // print_debug_gfa(sg, init_ug, uopt->coverage_cut, bin_file, uopt->sources, uopt->ruIndex, uopt->max_hang, uopt->min_ovlp, 0, 0, 0);
    // }
    // if(is_trio) gen_ul_trio_graph(uopt, uidx, o_file);    
    // exit(0);
    // return uidx->uovl.hybrid_ug;

    asg_t *ng = renew_ng(uidx->uovl.h_usg, init_ug, sg, uopt, ropt->cov, ropt->ruIndex, 16, bin_file);
    destory_bubbles(uidx->bub); free(uidx->bub); ma_ug_destroy(init_ug); free(r_het);
    destroy_ul_resolve_t(uidx); 
    if(free_uld) {
        destory_all_ul_t(&UL_INF); memset((&UL_INF), 0, sizeof(UL_INF));
    } else {
        for (i = 0; i < bn; i++) UL_INF.a[i].rlen = bl[i];
        for (i = bn; i < UL_INF.n; i++) {
            free(UL_INF.a[i].N_site.a); free(UL_INF.a[i].r_base.a); 
            free(UL_INF.a[i].bb.a); memset(&(UL_INF.a[i]), 0, sizeof(UL_INF.a[i]));
        }
        UL_INF.n = bn;

        for (i = idn; i < UL_INF.nid.n; i++) {
            free(UL_INF.nid.a[i].a); memset(&(UL_INF.nid.a[i]), 0, sizeof(UL_INF.nid.a[i])); 
        }
        UL_INF.nid.n = idn; free(bl);
    }

    asg_destroy((*(ropt->sg))); (*(ropt->sg)) = ng; 
    (*(ropt->src)) = R_INF.paf; (*(ropt->r_src)) = R_INF.reverse_paf;
    (*(ropt->n_read)) = R_INF.total_reads; (*(ropt->readLen)) = R_INF.read_length;
    renew_R_to_U(ng, (*(ropt->src)), (*(ropt->r_src)), (*(ropt->n_read)), (*(ropt->cov)), ropt->ruIndex, ropt->max_hang, ropt->mini_ovlp);
    post_rescue(uopt, (*(ropt->sg)), (*(ropt->src)), (*(ropt->r_src)), ropt->ruIndex, ropt->b_mask_t, 0);
    // print_raw_uls_aln(uidx, asm_opt.output_file_name);
    // exit(0);
}

ug_clean_t *init_ug_clean_t(ug_opt_t *uopt, asg_t *sg, uint8_t is_ou, uint8_t is_trio, int64_t max_ext, double len_rat, double ou_rat, int64_t min_ou)
{
    uint64_t i;
    ug_clean_t *sl; CALLOC(sl, 1);
    ma_ug_t *ug = ma_ug_gen(sg);
    for (i = 0; i < ug->g->n_seq; ++i) {
        if(ug->g->seq[i].del) continue;
        ug->g->seq[i].c = PRIMARY_LABLE;
    }
    MALLOC(sl->idx, (ug->g->n_seq<<1)); memset(sl->idx, -1, sizeof((*(sl->idx)))*(ug->g->n_seq<<1));
    MALLOC(sl->bid, (ug->g->n_seq<<1)); memset(sl->bid, -1, sizeof((*(sl->bid)))*(ug->g->n_seq<<1));
    sl->uopt = uopt; sl->sg = sg; sl->ug = ug; sl->idx_n = sl->bid_n = ug->g->n_seq<<1;
    memset(&(sl->b), 0, sizeof(sl->b)); CALLOC(sl->b.a, (ug->g->n_seq<<1));
    sl->is_ou = is_ou;
    sl->is_trio = is_trio;
    sl->max_ext = max_ext;
    sl->len_rat = len_rat;
    sl->ou_rat = ou_rat;
    sl->min_ou = min_ou;
    if(is_ou) update_ug_ou(sl->ug, sl->sg);
    return sl;
}

void destroy_ug_clean_t(ug_clean_t *sl)
{
    free(sl->idx); free(sl->bid); ma_ug_destroy(sl->ug);
    free(sl->b.a); free(sl->b.S.a); free(sl->b.T.a); free(sl->b.b.a); free(sl->b.e.a); 
}

void update_ug_clean_t(ug_clean_t *sl)
{
    uint32_t v, k, i, n_vtx = sl->ug->g->n_seq<<1;
    sl->tlen = get_bub_pop_max_dist_advance(sl->ug->g, &(sl->b));
    for (k = 0; k < sl->ug->g->n_seq; ++k) {
        sl->idx[k<<1] = sl->idx[(k<<1)+1] = sl->bid[k<<1] = sl->bid[(k<<1)+1] = (uint32_t)-1;
        if(sl->ug->g->seq[k].del) continue;
        sl->ug->g->seq[k].c = PRIMARY_LABLE;
    }

    for (v = 0; v < n_vtx; ++v) {
        if(sl->ug->g->seq[v>>1].del) continue;
        if(asg_arc_n(sl->ug->g, v) < 2) continue;
        if(sl->bid[v] != ((uint32_t)-1)) continue;
        if(asg_bub_pop1_primary_trio(sl->ug->g, NULL, v, sl->tlen, &(sl->b), (uint32_t)-1, (uint32_t)-1, 0, NULL, NULL, NULL, 0, 0, NULL)) {
            //beg is v, end is b.S.a[0]
            //note b.b include end, does not include beg
            for (i = 0; i < sl->b.b.n; i++) {
                if(sl->b.b.a[i]==v || sl->b.b.a[i]==sl->b.S.a[0]) continue;
                sl->bid[sl->b.b.a[i]] = sl->bid[sl->b.b.a[i]^1] = 1;
            }
            sl->bid[v] = 2; sl->bid[sl->b.S.a[0]^1] = 3;
        }
    }

    for (v = 0; v < n_vtx; ++v) {
        if(sl->bid[v] !=2) continue;
        if(asg_bub_pop1_primary_trio(sl->ug->g, NULL, v, sl->tlen, &(sl->b), (uint32_t)-1, (uint32_t)-1, 0, NULL, NULL, NULL, 0, 0, NULL)) {
            //note b.b include end, does not include beg
            for (i = 0; i < sl->b.b.n; i++) {
                if(sl->b.b.a[i]==v || sl->b.b.a[i]==sl->b.S.a[0]) continue;
                sl->idx[sl->b.b.a[i]] = v;
                sl->idx[sl->b.b.a[i]^1] = sl->b.S.a[0]^1;
            }
            sl->idx[v] = v; sl->idx[sl->b.S.a[0]^1] = sl->b.S.a[0]^1;
        }
    }
    memcpy(sl->bid, sl->idx, sizeof((*(sl->bid)))*n_vtx);
    // memset(sl->idx, -1, sizeof((*(sl->idx)))*n_vtx);
}

uint32_t get_best_het_node(ma_ug_t *ug, uint32_t v, uint32_t *bid, uint8_t is_ou, uint8_t is_trio, double len_rat, double ou_rat, uint32_t min_ou)
{
    uint32_t w = (uint32_t)-1, k, nv, kv, ol_max, ol_k;
    asg_arc_t *av; uint32_t trioF = (uint32_t)-1, ntrioF = (uint32_t)-1;
    if(ug->g->seq[v>>1].del) return (uint32_t)-1;
    av = asg_arc_a(ug->g, v); nv = asg_arc_n(ug->g, v);
    for (k = kv = 0; k < nv; k++) {
        if(av[k].del) continue;
        ///could not connect to the beg/sink node
        if((av[k].v>>1)==(bid[v]>>1) || (av[k].v>>1)==(bid[v^1]>>1)) return (uint32_t)-1;
        kv++; w = av[k].v;
    }
    if(kv < 2) return w;///kv == 0, return (uint32_t)-1; kv == 1, return node; 
    // if(v == 351) fprintf(stderr, "-0-[M::%s]\tutg%.6ul(%c)\tv::%u\tkv::%u\n", __func__, (v>>1)+1, "+-"[v&1], v, kv);

    if(is_trio) {
        trioF = get_ug_tip_trio_infor(ug, v^1);
        ntrioF = (trioF==FATHER? MOTHER : (trioF==MOTHER? FATHER : (uint32_t)-1));
    }

    ol_max = 0; ol_k = (uint32_t)-1;
    for (k = 0; k < nv; ++k) {
        if(av[k].del) continue;
        // if(is_trio && get_ug_tip_trio_infor(ug, av[k].v) == ntrioF) continue;
        if(ol_max < av[k].ol) ol_max = av[k].ol, ol_k = k;
    }
    // if(v == 351) fprintf(stderr, "-1-[M::%s]\tutg%.6ul(%c)\tv::%u\tol_k::%u\n", __func__, (v>>1)+1, "+-"[v&1], v, ol_k);
    if(ol_k == (uint32_t)-1) return (uint32_t)-1;
    for (k = 0; k < nv; ++k) {
        if(av[k].del || k == ol_k) continue;
        // if(is_trio && get_ug_tip_trio_infor(ug, av[k].v) == ntrioF) continue;
        if(av[k].ol > av[ol_k].ol*len_rat) return (uint32_t)-1;
        if((is_ou) && (av[k].ou > min_ou) && ((av[k].ou) > (av[ol_k].ou*ou_rat))) return (uint32_t)-1;
    }
    if(is_trio && get_ug_tip_trio_infor(ug, av[ol_k].v) == ntrioF) return (uint32_t)-1;
    // if(v == 351) fprintf(stderr, "-2-[M::%s]\tutg%.6ul(%c)\tv::%u\tav[ol_k].v::%u\n", __func__, (v>>1)+1, "+-"[v&1], v, av[ol_k].v);
    return av[ol_k].v;
}

static void cal_bub_best(void *data, long i, int tid)
{
    ug_clean_t *sl = (ug_clean_t *)data;
    ma_ug_t *ug = sl->ug; uint32_t v, w;
    sl->idx[i<<1] = sl->idx[(i<<1)+1] = (uint32_t)-1;
    if(ug->g->seq[i].del) return;
    if((sl->bid[i<<1] == (uint32_t)-1) || (sl->bid[(i<<1)+1] == (uint32_t)-1)) return;///not within a bubble
    
    v = i<<1;
    if(sl->bid[v] != v) {///not the beg/sink node
        w = get_best_het_node(ug, v, sl->bid, sl->is_ou, sl->is_trio, sl->len_rat, sl->ou_rat, sl->min_ou);
        // if(v == 352) fprintf(stderr, "-v-[M::%s]\tutg%.6ul(%c)\tv::%u\tw::%u\n", __func__, (v>>1)+1, "+-"[v&1], v, w);
        if((w != (uint32_t)-1) && ((v^1) == get_best_het_node(ug, w^1, sl->bid, sl->is_ou, sl->is_trio, sl->len_rat, sl->ou_rat, sl->min_ou))) {
            sl->idx[v] = w;
        }
    }

    v = (i<<1)+1;
    if(sl->bid[v] != v) {///not the beg/sink node
        w = get_best_het_node(ug, v, sl->bid, sl->is_ou, sl->is_trio, sl->len_rat, sl->ou_rat, sl->min_ou);
        if((w != (uint32_t)-1) && ((v^1) == get_best_het_node(ug, w^1, sl->bid, sl->is_ou, sl->is_trio, sl->len_rat, sl->ou_rat, sl->min_ou))) {
            sl->idx[v] = w;
        }
    }
}

uint32_t usg_topocut_aux_unambi1(asg_t *g, uint32_t v)
{
	asg_arc_t *av = asg_arc_a(g, v);
	uint32_t i, nv = asg_arc_n(g, v);
	uint32_t k = nv, kv;
	for (i = 0, kv = 0; i < nv; ++i)
		if (!av[i].del) ++kv, k = i;
	if (kv != 1) return (uint32_t)-1;
	return av[k].v;
}

int usg_topocut_aux_del(ma_ug_t *ug, uint32_t v, int max_ext, asg64_v *b)
{
	int32_t n_ext; 
	for (n_ext = 0; n_ext < max_ext && v != (uint32_t)-1; ) {
		if (usg_topocut_aux_unambi1(ug->g, v^1) == (uint32_t)-1) break;
        if(b) kv_push(uint64_t, *b, v); n_ext += ug->u.a[v>>1].n;
		v = usg_topocut_aux_unambi1(ug->g, v);
	}
	return n_ext;
}

#define is_best_arc(sl, v, w) (((sl).idx[(v)]==(w))&&((sl).idx[((w)^1)]==((v)^1)))

inline void asg_arc_del_by_ug(asg_t *sg, ma_ug_t *ug, uint32_t uv, uint32_t uw, uint32_t del)
{
    uint32_t rv, rw, i, nv; asg_arc_t *av;
    if(uv&1) rv = ug->u.a[uv>>1].start^1;
    else rv = ug->u.a[uv>>1].end^1;

    if(uw&1) rw = ug->u.a[uw>>1].end;
    else rw = ug->u.a[uw>>1].start;

    av = asg_arc_a(sg, rv);
    nv = asg_arc_n(sg, rv);
    for (i = 0; i < nv; i++) {
        if(av[i].v == rw) break;
    }
    assert(i < nv);
    av[i].del = del;
}

uint32_t cal_utg_occ(ma_ug_t *ug, uint32_t begNode)
{
    uint32_t v = begNode, w, kv, occ = 0; ma_utg_t *u;

    while (1) {
        kv = get_real_length(ug->g, v, NULL);
        u = &(ug->u.a[v>>1]); occ += u->n; 
        if(kv!=1) break;
        ///kv must be 1 here
        kv = get_real_length(ug->g, v, &w);
        if(get_real_length(ug->g, w^1, NULL)!=1) break;
        v = w;
        if(v == begNode) break;
    }
    return occ;
}

void cal_bub_best_by_len(ug_clean_t *sl, asg64_v *in, uint32_t max_ext, uint32_t is_trio, uint32_t is_ou, double len_rat, double ou_rat, uint32_t min_ou, 
uint32_t min_node)
{
    // fprintf(stderr, "[M::%s]\tStart\n", __func__);
    ma_ug_t *ug = sl->ug; asg_t *g = sl->ug->g; uint32_t ol_max, ou_max, lnid;
    uint32_t v, w, n_vtx = (g->n_seq<<1), nv, nw, i, k, z, kv, kw, bb, to_del, bn, tip;
    asg64_v tx = {0,0,0}, *b = NULL; asg_arc_t *av, *aw, *ve, *we; ma_utg_t *u;
    uint32_t trioF = (uint32_t)-1, ntrioF = (uint32_t)-1, mm_ol, mm_ou, cnt = 0, del_v, del_w;

    if(in) b = in;
    else b = &tx;
    b->n = 0;

    for (v = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;
        av = asg_arc_a(g, v); 
        nv = asg_arc_n(g, v);
        if (nv < 2) continue;

        for (i = kv = bb = 0; i < nv; ++i) {
            if(av[i].del) continue;
            if(is_best_arc((*sl), v, av[i].v)) bb++;
            kv++;
        }
        if((kv < 2) || (!bb) || (kv<=bb)) continue;///it is impossible that kv <= bv

        for (i = 0; i < nv; ++i) {
            if(av[i].del) continue;
            if(is_best_arc((*sl), v, av[i].v)) continue;
            kv_push(uint64_t, *b, (((uint64_t)av[i].ol)<<32) | ((uint64_t)(av-g->arc+i)));   
        }
    }

    radix_sort_srt64(b->a, b->a + b->n);
    for (k = 0; k < b->n; k++) {
        if(g->arc[(uint32_t)b->a[k]].del) continue;

        v = g->arc[(uint32_t)b->a[k]].ul>>32; w = g->arc[(uint32_t)b->a[k]].v^1;
        if(g->seq[v>>1].del || g->seq[w>>1].del) continue;
        nv = asg_arc_n(g, v); nw = asg_arc_n(g, w);
        av = asg_arc_a(g, v); aw = asg_arc_a(g, w);
        if(nv<=1 && nw <= 1) continue;

        if(is_trio) {
            if(get_arcs(g, v, NULL, 0)<=1 && get_arcs(g, w, NULL, 0)<=1) continue;///speedup
            trioF = get_ug_tip_trio_infor(ug, v^1);
            ntrioF = (trioF==FATHER? MOTHER : (trioF==MOTHER? FATHER : (uint32_t)-1));
        }

        ve = &(g->arc[(uint32_t)b->a[k]]);
        for (i = 0; i < nw; ++i) {
            if (aw[i].v == (v^1)) {
                we = &(aw[i]);
                break;
            }
        }
        ///mm_ol and mm_ou are used to make edge with long indel more easy to be cutted
        mm_ol = MIN(ve->ol, we->ol); mm_ou = MIN(ve->ou, we->ou); lnid = 0;

        for (i = kv = ol_max = ou_max = 0; i < nv; ++i) {
            if(av[i].del) continue;
            kv++; 
            if(is_trio && get_ug_tip_trio_infor(ug, av[i].v) == ntrioF) continue;
            if(ol_max < av[i].ol) ol_max = av[i].ol; 
            if(ou_max < av[i].ou) ou_max = av[i].ou;
            if((!lnid) && ((is_best_arc((*sl), v, av[i].v)))) {
                if((cal_utg_occ(ug, v^1) >= min_node) || (cal_utg_occ(ug, av[i].v) >= min_node)) lnid = 1;
            }
        }
        if (kv < 1) continue;
        if (kv >= 2) {
            if (mm_ol > ol_max*len_rat) continue;
            if ((is_ou) && (mm_ou > min_ou) && (mm_ou > (ou_max*ou_rat))) continue;
        }
        

        for (i = kw = ol_max = ou_max = 0; i < nw; ++i) {
            if(aw[i].del) continue;
            kw++; 
            if(is_trio && get_ug_tip_trio_infor(ug, aw[i].v) == ntrioF) continue;
            if(ol_max < aw[i].ol) ol_max = aw[i].ol; 
            if(ou_max < aw[i].ou) ou_max = aw[i].ou;
            if((!lnid) && ((is_best_arc((*sl), w, aw[i].v)))) {
                if((cal_utg_occ(ug, w^1) >= min_node) || (cal_utg_occ(ug, aw[i].v) >= min_node)) lnid = 1;
            }
        }
        if (kw < 1) continue;
        if (kw >= 2) {
            if (mm_ol > ol_max*len_rat) continue;
            if ((is_ou) && (mm_ou > min_ou) && (mm_ou > (ou_max*ou_rat))) continue;
        }

        if (kv <= 1 && kw <= 1) continue;
        if(!lnid) continue;

        to_del = 0; del_v = del_w = (uint32_t)-1; tip = 0;
        if (kv > 1 && kw > 1) {
            to_del = 1;
        } else if (kw == 1) {
            // tip = asg_topocut_aux(g, w^1, max_ext);
            tip = usg_topocut_aux_del(ug, w^1, max_ext, NULL);
            if (tip < max_ext) {
                to_del = 1; del_w = w^1;
            }
        } else if (kv == 1) {
            // tip = asg_topocut_aux(g, v^1, max_ext);
            tip = usg_topocut_aux_del(ug, v^1, max_ext, NULL);
            if (tip < max_ext) {
                to_del = 1; del_v = v^1;
            }
        }

        
        if (to_del) {
            bn = b->n;
            if(del_v != ((uint32_t)-1)) usg_topocut_aux_del(ug, del_v, max_ext, b);
            if(del_w != ((uint32_t)-1)) usg_topocut_aux_del(ug, del_w, max_ext, b);
            for (i = bn; i < b->n; i++) {
                u = &(ug->u.a[b->a[i]>>1]);
                if(u->m == 0) continue;
                for (z = 0; z < u->n; z++) asg_seq_del(sl->sg, u->a[z]>>33);
                asg_seq_del(ug->g, (b->a[i]>>1));
                if(u->m) {
                    u->m = u->n = 0; free(u->a); u->a = NULL;
                }
            }
            // assert(tip == (b->n-bn));
            b->n = bn;

            ve->del = we->del = 1, ++cnt;
            asg_arc_del_by_ug(sl->sg, ug, ve->ul>>32, ve->v, 1);
            asg_arc_del_by_ug(sl->sg, ug, we->ul>>32, we->v, 1);
            
        }
    }

    if(!in) free(tx.a);
    if (cnt > 0) {
        asg_cleanup(g); asg_cleanup(sl->sg);
    }
    // fprintf(stderr, "[M::%s]\tEnd\n", __func__);
}

/**
void cal_bub_best_by_topo(ug_clean_t *sl, asg64_v *in, uint32_t max_ext, uint32_t is_trio, uint32_t is_ou, double len_rat, double ou_rat, uint32_t min_ou, uint32_t long_tip)
{
    // fprintf(stderr, "[M::%s]\tStart\n", __func__);
    ma_ug_t *ug = sl->ug; asg_t *g = sl->ug->g; uint32_t ol_max, ou_max;
    uint32_t v, w, n_vtx = (g->n_seq<<1), nv, nw, i, k, z, kv, kw, to_del, bn, tip, avi, awi;
    asg64_v tx = {0,0,0}, *b = NULL; asg_arc_t *av, *aw, *ve, *we; ma_utg_t *u;
    uint32_t trioF = (uint32_t)-1, mm_ol, mm_ou, cnt = 0, del_v, del_w;

    if(in) b = in;
    else b = &tx;
    b->n = 0;

    for (v = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;
        if((sl->bid[v] == (uint32_t)-1) || (sl->bid[v^1] == (uint32_t)-1)) continue;///not within a bubble

        av = asg_arc_a(g, v); 
        nv = asg_arc_n(g, v);
        if (nv < 2) continue;

        for (i = kv = 0; i < nv; ++i) {
            if(av[i].del) continue;
            ///could not connect to the beg/sink node
            if((av[k].v>>1)==(sl->bid[v]>>1) || (av[k].v>>1)==(sl->bid[v^1]>>1)) break;
            kv++;
        }
        if(kv < 2 || i < nv) continue;///it is impossible that kv <= bv

        for (i = 0; i < nv; ++i) {
            if(av[i].del) continue;
            kv_push(uint64_t, *b, (((uint64_t)av[i].ol)<<32) | ((uint64_t)(av-g->arc+i)));   
        }
    }

    radix_sort_srt64(b->a, b->a + b->n);
    for (k = 0; k < b->n; k++) {
        if(g->arc[(uint32_t)b->a[k]].del) continue;

        v = g->arc[(uint32_t)b->a[k]].ul>>32; w = g->arc[(uint32_t)b->a[k]].v^1;
        if(g->seq[v>>1].del || g->seq[w>>1].del) continue;
        nv = asg_arc_n(g, v); nw = asg_arc_n(g, w);
        av = asg_arc_a(g, v); aw = asg_arc_a(g, w);
        if(nv < 2 || nw < 2) continue;
        kv = get_arcs(g, v, NULL, 0); kw = get_arcs(g, w, NULL, 0);
        if(kv < 2 || kw < 2) continue;

        if(is_trio) {
            trioF = get_ug_tip_trio_infor(ug, v^1);
            if(trioF == FATHER || trioF == MOTHER) {
                if(get_ug_tip_trio_infor(ug, w^1) == trioF) continue;
            }
        }

        avi = awi = (uint32_t)-1;

        avi = ((uint32_t)b->a[k]) - ((uint64_t)(av-g->arc));
        ve = &(g->arc[(uint32_t)b->a[k]]);

        for (i = 0; i < nw; ++i) {
            if (aw[i].v == (v^1)) {
                we = &(aw[i]); awi = i;
                break;
            }
        }

        for (i = kv = ol_max = ou_max = 0; i < nv; ++i) {
            if(av[i].del) continue;
        }

        ///mm_ol and mm_ou are used to make edge with long indel more easy to be cutted
        mm_ol = MIN(ve->ol, we->ol); mm_ou = MIN(ve->ou, we->ou);

        for (i = kv = ol_max = ou_max = 0; i < nv; ++i) {
            if(av[i].del) continue;
            kv++; 
            if(is_trio && get_ug_tip_trio_infor(ug, av[i].v) == ntrioF) continue;
            if(ol_max < av[i].ol) ol_max = av[i].ol; 
            if(ou_max < av[i].ou) ou_max = av[i].ou;
        }
        if (kv < 1) continue;
        if (kv >= 2) {
            if (mm_ol > ol_max*len_rat) continue;
            if ((is_ou) && (mm_ou > min_ou) && (mm_ou > (ou_max*ou_rat))) continue;
        }
        

        for (i = kw = ol_max = ou_max = 0; i < nw; ++i) {
            if(aw[i].del) continue;
            kw++; 
            if(is_trio && get_ug_tip_trio_infor(ug, aw[i].v) == ntrioF) continue;
            if(ol_max < aw[i].ol) ol_max = aw[i].ol; 
            if(ou_max < aw[i].ou) ou_max = aw[i].ou;
        }
        if (kw < 1) continue;
        if (kw >= 2) {
            if (mm_ol > ol_max*len_rat) continue;
            if ((is_ou) && (mm_ou > min_ou) && (mm_ou > (ou_max*ou_rat))) continue;
        }

        if (kv <= 1 && kw <= 1) continue;

        to_del = 0; del_v = del_w = (uint32_t)-1; tip = 0;
        if (kv > 1 && kw > 1) {
            to_del = 1;
        } else if (kw == 1) {
            tip = asg_topocut_aux(g, w^1, max_ext);
            if (tip < max_ext) {
                to_del = 1; del_w = w^1;
            }
        } else if (kv == 1) {
            tip = asg_topocut_aux(g, v^1, max_ext);
            if (tip < max_ext) {
                to_del = 1; del_v = v^1;
            }
        }

        if (to_del) {
            bn = b->n;
            if(del_v != ((uint32_t)-1)) usg_topocut_aux_del(ug, del_v, max_ext, b);
            if(del_w != ((uint32_t)-1)) usg_topocut_aux_del(ug, del_w, max_ext, b);
            for (i = bn; i < b->n; i++) {
                u = &(ug->u.a[b->a[i]>>1]);
                if(u->m == 0) continue;
                for (z = 0; z < u->n; z++) asg_seq_del(sl->sg, u->a[z]>>33);
                asg_seq_del(ug->g, (b->a[i]>>1));
                if(u->m) {
                    u->m = u->n = 0; free(u->a); u->a = NULL;
                }
            }
            // assert(tip == (b->n-bn));
            b->n = bn;

            ve->del = we->del = 1, ++cnt;
            asg_arc_del_by_ug(sl->sg, ug, ve->ul>>32, ve->v, 1);
            asg_arc_del_by_ug(sl->sg, ug, we->ul>>32, we->v, 1);
            
        }
    }

    if(!in) free(tx.a);
    if (cnt > 0) {
        asg_cleanup(g); asg_cleanup(sl->sg);
    }
    // fprintf(stderr, "[M::%s]\tEnd\n", __func__);
}
**/



void deep_graph_clean(ug_opt_t *uopt, asg_t *sg, uint8_t is_ou, uint8_t is_trio, int64_t max_ext, 
double min_ovlp_drop_ratio, double max_ovlp_drop_ratio, double ou_rat, int64_t min_ou, int64_t clean_round, int64_t long_tip)
{
    ug_clean_t *sl; asg64_v b; double step; int64_t i;
    kv_init(b); hic_clean_adv(sg, uopt);
    step = (clean_round==1?max_ovlp_drop_ratio:((max_ovlp_drop_ratio-min_ovlp_drop_ratio)/(clean_round-1)));
    sl = init_ug_clean_t(uopt, sg, is_ou, is_trio, max_ext, min_ovlp_drop_ratio, ou_rat, min_ou);

    // print_debug_gfa(sg, NULL, uopt->coverage_cut, "debug", uopt->sources, uopt->ruIndex, uopt->max_hang, uopt->min_ovlp, 0, 0, 0);
    for (i = 0, sl->len_rat = min_ovlp_drop_ratio; i < clean_round; i++, sl->len_rat += step) {
        if(sl->len_rat > max_ovlp_drop_ratio) sl->len_rat = max_ovlp_drop_ratio;
        update_ug_clean_t(sl);
        kt_for(asm_opt.thread_num, cal_bub_best, sl, sl->ug->g->n_seq);
        cal_bub_best_by_len(sl, &b, max_ext, is_trio, sl->is_ou, sl->len_rat, sl->ou_rat, min_ou, long_tip);
    }

    // if(is_ou) {
    //     uint64_t k; sl->is_ou = 0;
    //     for (k = 0; k < sl->ug->g->n_arc; k++) sl->ug->g->arc[k].ou = 0;
    //     max_ovlp_drop_ratio = 0.6;
    //     if(max_ovlp_drop_ratio > min_ovlp_drop_ratio) {
    //         step = (clean_round==1?max_ovlp_drop_ratio:((max_ovlp_drop_ratio-min_ovlp_drop_ratio)/(clean_round-1)));
    //         for (i = 0, sl->len_rat = min_ovlp_drop_ratio; i < clean_round; i++, sl->len_rat += step) {
    //             if(sl->len_rat > max_ovlp_drop_ratio) sl->len_rat = max_ovlp_drop_ratio;
    //             update_ug_clean_t(sl);
    //             kt_for(asm_opt.thread_num, cal_bub_best, sl, sl->ug->g->n_seq);
    //             cal_bub_best_by_len(sl, &b, max_ext, is_trio, sl->is_ou, sl->len_rat, sl->ou_rat, min_ou);
    //         }
    //     }
    // }

    // update_ug_clean_t(sl);
    // cal_bub_best_by_topo(sl, &b, max_ext, is_trio, sl->is_ou, sl->len_rat, sl->ou_rat, min_ou, long_tip);

    kv_destroy(b); destroy_ug_clean_t(sl); free(sl);
    hic_clean_adv(sg, uopt);
}
