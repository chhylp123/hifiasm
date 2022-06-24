#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <zlib.h>
#include <math.h>
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

typedef struct {
	asg_t *g; 
	ma_hit_t_alloc *src;
} sset_aux;


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

typedef struct {
    uint64_t c; // max count of positive reads
    uint32_t i;
} uinfo_srt_t;

#define uinfo_srt_t_c_key(p) ((p).c)
KRADIX_SORT_INIT(uinfo_srt_t_c, uinfo_srt_t, uinfo_srt_t_c_key, member_size(uinfo_srt_t, c))

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
    // ul_path_srt_t psrt;
}ul_resolve_t;

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

uint32_t asg_arc_cut_tips(asg_t *g, uint32_t max_ext, asg64_v *in, uint32_t is_ou)
{
    asg64_v tx = {0,0,0}, *b = NULL;
    uint32_t n_vtx = g->n_seq<<1, v, w, i, k, cnt = 0, nv, kv, pb, ou, mm_ou;
    asg_arc_t *av = NULL; uint64_t lw;
	if(in) b = in;
    else b = &tx;
    b->n = 0;
    for (v = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;

        av = asg_arc_a(g, v^1); nv = asg_arc_n(g, v^1);
        for (i = kv = 0; i < nv; i++) {
            if (av[i].del) continue;
            kv++; break;
        }

        if(kv) continue;
        kv = 1; mm_ou = (uint32_t)-1; ou = 0;
        for (i = 0, w = v; i < max_ext; i++) {
            if(asg_end(g, w^1, &lw, is_ou?&ou:NULL)!=0) break;
            w = (uint32_t)lw; kv++; mm_ou = MIN(mm_ou, ou);
        }
		if(mm_ou == (uint32_t)-1) mm_ou = 0;
		kv += mm_ou; i += mm_ou;
        if(i < max_ext/** + (!!is_ou)**/) kv_push(uint64_t, *b, (((uint64_t)kv)<<32)|v);  
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
        pb = b->n; kv_push(uint64_t, *b, v); mm_ou = (uint32_t)-1; ou = 0;
        for (i = 0, w = v; i < max_ext; i++) {
            if(asg_end(g, w^1, &lw, is_ou?&ou:NULL)!=0) break;
            w = (uint32_t)lw; kv_push(uint64_t, *b, lw); mm_ou = MIN(mm_ou, ou);
        }
		if(mm_ou == (uint32_t)-1) mm_ou = 0;
		i += mm_ou;

        if(i < max_ext/** + (!!is_ou)**/) {
            for (i = pb; i < b->n; i++) asg_seq_del(g, ((uint32_t)b->a[i])>>1);
            cnt++;
        }
        b->n = pb; 
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
void asg_arc_cut_chimeric(asg_t *g, ma_hit_t_alloc* src, asg64_v *in, uint32_t ou_thres)
{
	asg64_v tx = {0,0,0}, *b = NULL;
	uint32_t v, w, ei[2] = {0}, k, i, n_vtx = g->n_seq<<1;
    uint32_t nw, el_n, cnt = 0; asg_arc_t *aw;
	if(in) b = in;
    else b = &tx;
	b->n = 0;

	for (v = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;
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
        asg_seq_del(g, v>>1);
        cnt++;
    }
    // stats_sysm(g);
	if(!in) free(tx.a);
    if (cnt > 0) asg_cleanup(g);
}

void asg_arc_cut_inexact(asg_t *g, ma_hit_t_alloc* src, asg64_v *in, int32_t max_ext, uint32_t is_ou, uint32_t is_trio/**, asg64_v *dbg**/)
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
            if (is_ou && mm_ou >= ou_max) continue;
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
            if (is_ou && mm_ou >= ou_max) continue;
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
uint32_t is_topo, ma_hit_t_alloc *rev, R_to_U* rI, uint32_t *max_drop_len)
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
            if (is_ou && mm_ou > ou_max*ou_rat) continue;
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

uint32_t asg_cut_chimeric_bub(asg_t *g, ma_hit_t_alloc* src, asg64_v *in, uint32_t normal_len, uint32_t is_clean)
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

void asg_iterative_semi_circ(asg_t *g, ma_hit_t_alloc* src, asg64_v *in, uint32_t normal_len, uint32_t pop_chimer)
{
    uint64_t occ = 0, s = 1;
    while (s) {
        s = asg_cut_semi_circ(g, LIM_LEN, 0);
        if(pop_chimer) s = s + asg_cut_chimeric_bub(g, src, in, normal_len, 0);
        occ += s;
    }

    // stats_sysm(g);
    if(occ) asg_cleanup(g);
}

uint32_t asg_cut_large_indel(asg_t *g, asg64_v *in, int32_t max_ext, float ou_rat, uint32_t is_ou)
{
    asg64_v tx = {0,0,0}, *b = NULL;
    uint32_t v, w, n_vtx = g->n_seq<<1, i, k, kv, kw, nv, nw, ou_max, to_del, cnt = 0;
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

        for (i = kv = ou_max = 0, ve = NULL; i < nv; ++i) {
            if(av[i].del) continue;
            if(av[i].v == (w^1)) ve = &(av[i]);
            kv++;
            if(ou_max < av[i].ou) ou_max = av[i].ou;
        }
        if (kv < 1) continue;
        if (kv >= 2) {
            if (is_ou && ve->ou > ou_max*ou_rat) continue;
        }
        

        for (i = kw = ou_max = 0, we = NULL; i < nw; ++i) {
            if(aw[i].del) continue;
            if(aw[i].v == (v^1)) we = &(aw[i]);
            kw++;
            if(ou_max < aw[i].ou) ou_max = aw[i].ou;
        }
        if (kw < 1) continue;
        if (kw >= 2) {
            if (is_ou && we->ou > ou_max*ou_rat) continue;
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

void print_vw_edge(asg_t *sg, uint32_t v, uint32_t w, const char *cmd)
{
    asg_arc_t *av; uint32_t nv, i;
    av = asg_arc_a(sg, v); nv = asg_arc_n(sg, v);
    for (i = 0; i < nv; i++) {
        if(av[i].v == w) {
            fprintf(stderr, "[%s]\t%.*s(%c)<id:%lu>\t%.*s(%c)<id:%u>\tol:%u\tou:%u\tdel:%u\n", cmd, 
                (int32_t)Get_NAME_LENGTH(R_INF, (av[i].ul>>33)), Get_NAME(R_INF, (av[i].ul>>33)), "+-"[(av[i].ul>>32)&1], av[i].ul>>33, 
                (int32_t)Get_NAME_LENGTH(R_INF, (av[i].v>>1)), Get_NAME(R_INF, (av[i].v>>1)), "+-"[av[i].v&1], av[i].v>>1, av[i].ol, av[i].ou, av[i].del);
            break;
        }
    }
    if(i >= nv) fprintf(stderr, "[%s]\tno edges\n", cmd);
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
    uint32_t i, m, v, w, nv, n_vx = ug->g->n_seq<<1, vx, wx; int32_t r;
    asg_arc_t *av = NULL; ma_utg_t *u = NULL; asg_arc_t *p, t;
    n_vx = rg->n_seq; rg->n_arc = 0;
    for (v = 0; v < n_vx; v++) rg->seq[v].del = (!!1);

    for (i = 0; i < ug->g->n_seq; ++i) {
        ug->g->seq[i].c = PRIMARY_LABLE;;
        if(ug->g->seq[i].del) continue;
        u = &(ug->u.a[i]);
        for (m = 0; m < u->n; m++) rg->seq[u->a[m]>>33].del = (!!0);
        for (m = 0; (m + 1) < u->n; m++){
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
    }

    free(rg->idx);
    rg->idx = 0;
    rg->is_srt = 0;
    asg_cleanup(rg);


    /*******************************for debug************************************/
    // ma_ug_t *dbg = ma_ug_gen(rg);
    // for (i = 0; i < dbg->g->n_seq; ++i) dbg->g->seq[i].c = PRIMARY_LABLE;  
    // cmp_untig_graph(dbg, ug);
    /*******************************for debug************************************/
}

void ul_clean_gfa(ug_opt_t *uopt, asg_t *sg, ma_hit_t_alloc *src, ma_hit_t_alloc *rev, R_to_U* rI, int64_t clean_round, double min_ovlp_drop_ratio, double max_ovlp_drop_ratio, 
double ou_drop_rate, int64_t max_tip, bub_label_t *b_mask_t, int32_t is_ou, int32_t is_trio, uint32_t ou_thres, char *o_file)
{
    #define HARD_OU_DROP 0.75
    #define HARD_OL_DROP 0.6
    #define HARD_OL_SEC_DROP 0.85
    #define HARD_ORTHOLOGY_DROP 0.4
    double step = 
        (clean_round==1?max_ovlp_drop_ratio:((max_ovlp_drop_ratio-min_ovlp_drop_ratio)/(clean_round-1)));
    double drop = min_ovlp_drop_ratio;
    int64_t i; asg64_v bu = {0,0,0}; uint32_t l_drop = 2000;
	if(is_ou) update_sg_uo(sg, src);

    // debug_info_of_specfic_node("m64012_190921_234837/111673711/ccs", sg, rI, "beg");
    // debug_info_of_specfic_node("m64011_190830_220126/95028102/ccs", sg, rI, "beg");

	asg_arc_cut_tips(sg, max_tip, &bu, is_ou);
    // fprintf(stderr, "[M::%s] count_edges_v_w(sg, 49778, 49847)->%ld\n", __func__, count_edges_v_w(sg, 49778, 49847));
    
    for (i = 0; i < clean_round; i++, drop += step) {
        if(drop > max_ovlp_drop_ratio) drop = max_ovlp_drop_ratio;
        // fprintf(stderr, "(0):i->%ld, drop->%f\n", i, drop);
        // print_vw_edge(sg, 34156, 34090, "0");
        // stats_chimeric(sg, src, &bu);
        if(!is_ou) asg_iterative_semi_circ(sg, src, &bu, max_tip, 1);

        asg_arc_identify_simple_bubbles_multi(sg, b_mask_t, 1);
        asg_arc_cut_chimeric(sg, src, &bu, is_ou?ou_thres:(uint32_t)-1);
        asg_arc_cut_tips(sg, max_tip, &bu, is_ou);

        asg_arc_identify_simple_bubbles_multi(sg, b_mask_t, 0);
        asg_arc_cut_inexact(sg, src, &bu, max_tip, is_ou, is_trio/**, NULL**//**&dbg**/);
        // debug_edges(&dbg, d, 2);
        asg_arc_cut_tips(sg, max_tip, &bu, is_ou);

        asg_arc_identify_simple_bubbles_multi(sg, b_mask_t, 1);
        asg_arc_cut_length(sg, &bu, max_tip, drop, ou_drop_rate, is_ou, is_trio, 1, NULL, NULL, NULL);
        asg_arc_cut_tips(sg, max_tip, &bu, is_ou);

        asg_arc_identify_simple_bubbles_multi(sg, b_mask_t, 1);
        asg_arc_cut_bub_links(sg, &bu, HARD_OL_DROP, HARD_OL_SEC_DROP, HARD_OU_DROP, is_ou, asm_opt.large_pop_bubble_size, rev, rI, max_tip);
        
        asg_arc_identify_simple_bubbles_multi(sg, b_mask_t, 1);
        asg_arc_cut_complex_bub_links(sg, &bu, HARD_OL_DROP, HARD_OU_DROP, is_ou, b_mask_t);
        asg_arc_cut_tips(sg, max_tip, &bu, is_ou);

        if(is_ou) {
            if(ul_refine_alignment(uopt, sg)) update_sg_uo(sg, src);
        }
    }
    // debug_info_of_specfic_node("m64012_190921_234837/111673711/ccs", sg, rI, "end");
    // debug_info_of_specfic_node("m64011_190830_220126/95028102/ccs", sg, rI, "end");

    if(!is_ou) asg_iterative_semi_circ(sg, src, &bu, max_tip, 1);

    asg_arc_identify_simple_bubbles_multi(sg, b_mask_t, 0);
    asg_cut_large_indel(sg, &bu, max_tip, HARD_OU_DROP, is_ou);///shoule we ignore ou here?
    asg_arc_cut_tips(sg, max_tip, &bu, is_ou);

    ///asg_arc_del_triangular_directly might be unnecessary
    asg_arc_identify_simple_bubbles_multi(sg, b_mask_t, 0);
    asg_arc_cut_length(sg, &bu, max_tip, HARD_ORTHOLOGY_DROP/**min_ovlp_drop_ratio**/, ou_drop_rate, is_ou, 0/**is_trio**/, 0, rev, rI, NULL);
    asg_arc_cut_tips(sg, max_tip, &bu, is_ou);

    asg_arc_identify_simple_bubbles_multi(sg, b_mask_t, 0);
    asg_arc_cut_length(sg, &bu, max_tip, min_ovlp_drop_ratio, ou_drop_rate, is_ou, 0/**is_trio**/, 0, rev, rI, &l_drop);
    asg_arc_cut_tips(sg, max_tip, &bu, is_ou);

    if(!is_ou) asg_cut_semi_circ(sg, LIM_LEN, 1);


    rescue_contained_reads_aggressive(NULL, sg, src, uopt->coverage_cut, rI, uopt->max_hang, uopt->min_ovlp, 10, 1, 0, NULL, NULL, b_mask_t);
    rescue_missing_overlaps_aggressive(NULL, sg, src, uopt->coverage_cut, rI, uopt->max_hang, uopt->min_ovlp, 1, 0, NULL, b_mask_t);
    rescue_missing_overlaps_backward(NULL, sg, src, uopt->coverage_cut, rI, uopt->max_hang, uopt->min_ovlp, 10, 1, 0, b_mask_t);
    // rescue_wrong_overlaps_to_unitigs(NULL, sg, sources, reverse_sources, coverage_cut, ruIndex, 
    // max_hang_length, mini_overlap_length, bubble_dist, NULL);
    // rescue_no_coverage_aggressive(sg, sources, reverse_sources, &coverage_cut, ruIndex, max_hang_length, 
    // mini_overlap_length, bubble_dist, 10);
    set_hom_global_coverage(&asm_opt, sg, uopt->coverage_cut, src, rev, rI, uopt->max_hang, uopt->min_ovlp);
    rescue_bubble_by_chain(sg, uopt->coverage_cut, src, rev, (asm_opt.max_short_tip*2), 0.15, 3, rI, 0.05, 0.9, uopt->max_hang, uopt->min_ovlp, 10, uopt->gap_fuzz, b_mask_t);
    output_unitig_graph(sg, uopt->coverage_cut, o_file, src, rI, uopt->max_hang, uopt->min_ovlp);
    // flat_bubbles(sg, ruIndex->is_het); free(ruIndex->is_het); ruIndex->is_het = NULL;
    flat_soma_v(sg, src, rI);

    ///note: although the above functions will not change the UL part, but it will change the read graph
    ///so it is necessary to run update_sg_uo
    if(is_ou) {
        update_sg_uo(sg, src);
    }
    
    // print_node(sg, 17078); //print_node(sg, 8311); print_node(sg, 8294);
    
    free(bu.a); 
}


bubble_type *gen_bubble_chain(asg_t *sg, ma_ug_t *ug, ug_opt_t *uopt, uint8_t **ir_het)
{
    kvec_asg_arc_t_warp new_rtg_edges;
    kv_init(new_rtg_edges.a);
    hap_cov_t *cov = NULL; 
    bubble_type *bub = NULL; 

    asg_t *copy_sg = copy_read_graph(sg);
    ma_ug_t *copy_ug = copy_untig_graph(ug);   

    adjust_utg_by_primary(&copy_ug, copy_sg, TRIO_THRES, uopt->sources, uopt->reverse_sources, uopt->coverage_cut, 
    uopt->tipsLen, uopt->tip_drop_ratio, uopt->stops_threshold, uopt->ruIndex, uopt->chimeric_rate, uopt->drop_ratio,
    uopt->max_hang, uopt->min_ovlp, &new_rtg_edges, &cov, uopt->b_mask_t, 0, 0);
    ma_ug_destroy(copy_ug); copy_ug = NULL;
    asg_destroy(copy_sg); copy_sg = NULL; 

    CALLOC(bub, 1); (*ir_het) = cov->t_ch->ir_het; 
    cov->t_ch->ir_het = NULL; cov->is_r_het = NULL;
    identify_bubbles(ug, bub, (*ir_het), NULL);

    kv_destroy(new_rtg_edges.a); destory_hap_cov_t(&cov);
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
            kv_push(uint64_t, str->str.a[i], l);
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

ul_resolve_t *init_ul_resolve_t(asg_t *sg, ma_ug_t *init_ug, bubble_type* bub, all_ul_t *idx, uint8_t *r_het)
{
    ul_resolve_t *p = NULL; CALLOC(p, 1);
    p->sg = sg; p->init_ug = init_ug; p->bub = bub; p->idx = idx; p->r_het = r_het;
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
    if(i >= a_n) {
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
        if(str_id == 47072) {
            fprintf(stderr, "+[M::%s::k->%ld] p_g::%ld, c_g::%u, p_str::%ld, c_str::%ld, str_off::%ld, qoff::%lu, qstr_n::%lu\n", 
                __func__, k, p_g, gidx[a[k].tk], p_str, poa_str_idx(((uint32_t)a[k].tn_rev_qk), str_occ, is_rev), 
                str_off, qoff, qstr_n);
        }

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
    fprintf(stderr, "[M::%s::] i::%ld, a_n::%ld\n", __func__, i, a_n);
    
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
        fprintf(stderr, "[M::%s::i->%ld] qk::%u, tk::%u, max_k::%ld, max_f::%ld, ti::%ld\n", 
        __func__, i, (uint32_t)li->tn_rev_qk, li->tk, max_k, max_f, ti);
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
        fprintf(stderr, "\n[M::%s::] k::%ld, tid::%ld, is_rev::%ld\n", __func__, k, tid, is_rev);
        print_integer_seq(ug, str, tid, 1);
        poa_chain_0(g, ug, ul_idx->a[tid].bb.a, &(str[tid]), idx[k].t_sidx, idx[k].t_eidx, is_rev, buf, tid, qid, is_circle);
        if((*is_circle)) return;
        print_integer_g(g, ug, 1);
    }

    gen_cns_by_poa(g);
}

uint64_t cal_forward_dis(asg_t *g, uc_block_t *a, uint32_t s, uint32_t e)
{
    uint32_t i, li, v, w, nv, z; int64_t l; asg_arc_t *av;
    fprintf(stderr, "\n[M::%s::] s::%u, e::%u\n", __func__, s, e);
    ///TODO: a[i].aidx might be < i; if s == e, then i <= e might be wrong, cannot pass the assert(li == e);
    for(i = s, l = 0, v = w = (uint32_t)-1, li = s; i != (uint32_t)-1 && i <= /**!=**/ e; i = a[i].aidx) {
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
        v = w; li = i;
        fprintf(stderr, "[M::%s::] i::%u, li::%u, a[i].aidx::%u, a[i].qs::%u, a[i].qe::%u\n", __func__, i, li, a[i].aidx, 
        a[i].qs, a[i].qe);
    }
    assert(li == e);
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
    fprintf(stderr, "k_0::%ld, k_1::%ld\n", k_0, k_1);
    uint32_t li, lk, pk, bi = i; int64_t l;
    for (li = i, l = 0; i != (uint32_t)-1 && i >= k; i = a[i].pidx) {
        li = i; if(a[i].pidx != (uint32_t)-1 && i > k) l += a[i].pdis;
    }
    fprintf(stderr, "+bi::%u, li::%u, i::%u, l::%ld\n", bi, li, i, l);
    if(is_rev) l = cal_forward_dis(ug->g, a, li, bi);
    fprintf(stderr, "++bi::%u, li::%u, i::%u, l::%ld\n", bi, li, i, l);
    if(li == k) {///direct path
        (*is_g_connect) = 1;    
        return l;
    }
    i = li;
    assert(i > k); pk = k;
    for (lk = k; k != (uint32_t)-1 && k <= i; k = a[k].aidx) lk = k;
    fprintf(stderr, "-pk::%u, lk::%u, k::%u\n", pk, lk, k);
    if(!is_rev) {
        for (k = lk; k != (uint32_t)-1 && k != pk; k = a[k].pidx) l += a[k].pdis;
    } else {
        l += cal_forward_dis(ug->g, a, pk, lk);
    }
    fprintf(stderr, "--pk::%u, lk::%u, k::%u, l::%ld\n", pk, lk, k, l);
    if(l < 0) l = 0;
    k = lk;
    assert(i > k);
    return l + normlize_gdis(ug, &(a[i]), &(a[k]), is_rev);
}

uint64_t cal_integer_most_dis(uint64_t *a, uint64_t a_n, double cluster_rate)
{
    if(a_n <= 0) return (uint64_t)-1;
    // fprintf(stderr, "\n[M::%s::] a_n::%lu\n", __func__, a_n);
    uint64_t k, l, i, m, r_an = a_n, max_m, max_i, cc, cd, nd; int64_t z;
    for (k = 1, l = m = max_m = 0, max_i = (uint64_t)-1; k <= a_n; k++) {
        a[k-1] <<= 1; a[k-1] >>= 1;
        // fprintf(stderr, "[k->%lu] d::%lu, rev::%lu\n", k-1, a[k-1]>>1, a[k-1]&1);
        if(k == a_n || (a[k]>>1) != (a[l]>>1)) {
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
            nd = (((uint32_t)a[i])>>1); assert(nd > cd);
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

    uint32_t e_s, e_e, is_rev, is_g_connect, con_occ; emap_t *g_arc; uint64_t t, dd;
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
        fprintf(stderr, "[M::%s::] k::%lu, arc_idx_n::%lu, e_s::%u, e_e::%u\n", __func__, k, arc_idx_n, e_s, e_e);
        assert(poa_g_arc_w(pg, cns_seq[k], cns_seq[k+1]) == (e_e - e_s)); assert(e_e > e_s); 
        buf->o.n = 0; kv_resize(uint64_t, buf->o, e_e - e_s); con_occ = 0;
        for (i = e_s; i < e_e; i++) {
            g_arc = &(pg->e_idx.a[i]);
            assert((g_arc->pge>>32) == (cns_seq[k]>>1) && ((uint32_t)g_arc->pge) == (cns_seq[k+1]>>1));
            is_rev = ((g_arc->ule>>32) > ((uint32_t)g_arc->ule)?1:0);
            assert((pg->seq.a[g_arc->pge>>32].nid^(is_rev?1:0)) == ((uint32_t)str[g_arc->ulid].a[g_arc->ule>>32]));
            assert((pg->seq.a[(uint32_t)g_arc->pge].nid^(is_rev?1:0)) == ((uint32_t)str[g_arc->ulid].a[(uint32_t)g_arc->ule]));
            fprintf(stderr, "+++i::%lu, target_ulid::%u, str_sidx::%u, str_eidx::%u, str_cn::%u, ul_idx->a[g_arc->ulid].bb.n::%u\n", 
            i - e_s, g_arc->ulid, (uint32_t)(g_arc->ule>>32), (uint32_t)g_arc->ule, str[g_arc->ulid].cn, (uint32_t)ul_idx->a[g_arc->ulid].bb.n);
            dd = cal_integer_match_dis(ug, ul_idx->a[g_arc->ulid].bb.a, str[g_arc->ulid].a[g_arc->ule>>32]>>32, 
                                                                    str[g_arc->ulid].a[(uint32_t)g_arc->ule]>>32, is_rev, &is_g_connect);
            fprintf(stderr, "---i::%lu, dd::%lu, is_g_connect::%u\n", i - e_s, dd, is_g_connect);
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
        dd = cal_integer_most_dis(buf->o.a, buf->o.n, 0.08);


        t = pg->seq.a[cns_seq[k+1]>>1].nid; t |= ((uint64_t)(dd<<32)); 
        if(con_occ > 0) t |= ((uint64_t)(0x8000000000000000));
        kv_push(uint64_t, buf->res_dump, t);
    }
}

void integer_candidate(ul_resolve_t *uidx, integer_t *buf, uint32_t qid, uint32_t is_hom)
{
    // if(qid != 3165) return;
    // if(qid != 17165) return;
    // if(qid != 24100) return;///circle
    if(qid != 27512) return;
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
            m_het++; m_het_occ += buf->u.a[k];
        }
        ref_occ += buf->u.a[k];
        // fprintf(stderr, "[M::%s::k->%lu] buf->u.a[k]->%lu, ts->%u, te->%u, pchain->%u\n", __func__, k, buf->u.a[k], xi->ts, xi->te, xi->pchain);
    }
    // if((!is_hom) && (m_het < 2) && (m_het > 0)) return;///if all matched unitigs are hom, is ok
    if(m_het == 0 || m_het_occ == 0) is_hom = 1;
    print_ul_alignment(ug, &UL_INF, 27512, "inner-0");
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
    print_ul_alignment(ug, &UL_INF, 27512, "inner-1");
    radix_sort_integer_aln_t_srt(buf->b.a, buf->b.a + buf->b.n); 
    b_n = buf->b.n; buf->sc.n = 0;
    for (k = 1, z = 0; k < b_n; k++) {
        if(k == b_n || (buf->b.a[z].tn_rev_qk>>32) != (buf->b.a[k].tn_rev_qk>>32)) {
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
    print_ul_alignment(ug, &UL_INF, 27512, "inner-2");
    uint64_t *o, o_n, cns_het, cns_het_occ, ref_cns_occ, corrected = 0;
    o_n = integer_chain_dp(uidx->bub, buf, str_idx->str.a, buf->b.a, buf->sc.a, buf->sc.n, qid, is_hom, 2, &corrected);
    assert(o_n <= str->cn); 
    if(o_n <= 0) return; 
    if(corrected) return;
    print_ul_alignment(ug, &UL_INF, 27512, "inner-3");
    
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
    print_ul_alignment(ug, &UL_INF, 27512, "inner-4");
    fprintf(stderr, "\n");
    print_integer_seq(ug, str_idx->str.a, qid, 1);
    // print_aln_seq(ug, uidx->idx, qid, 1);
    print_cns_seq(ug, str, o, o_n);
    print_ul_alignment(ug, &UL_INF, 27512, "inner-5");

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
    print_ul_alignment(ug, &UL_INF, 27512, "inner-6");
    // if(m != str->cn) print_integer_ovlps(uidx->l1_ug, str_idx->str.a, buf->b.a, buf->b.n, buf->sc.a, buf->sc.n, qid, m);
    poa_cns_chain(&(buf->pg), uidx->idx, ug, str_idx->str.a, buf->sc.a, buf->sc.n, qid, buf, &is_circle);
    if(is_circle) {
        buf->n_circle++;
        return;
    }
    print_ul_alignment(ug, &UL_INF, 27512, "inner-7");
    print_res_seq(&(buf->pg), ug, buf->pg.srt_b.res.a, buf->pg.srt_b.res.n);
    // radix_sort_ul_chain_t_srt(buf->sc.a, buf->sc.a + buf->sc.n);
    
    // integer_phase(str_idx->str.a, buf, buf->sc.a, buf->sc.n, buf->b.a, qid);

    // radix_sort_ul_chain_t_srt(buf->sc.a, buf->sc.a + buf->sc.n);
    o = NULL; o_n = 0;
    update_raw_integer_seq(&(buf->pg), ug, buf->pg.srt_b.res.a, buf->pg.srt_b.res.n, uidx->idx, str_idx->str.a, qid, buf, buf->sc.a, buf->sc.n);
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
    // uc_block_t *uls; uint64_t uls_n; uint32_t v;
    // uint64_t *srt_a, srt_n, is_circle = ((uidx->psrt.idx.a[i]&((uint64_t)(0x100000000)))?1:0);
    // srt_a = uidx->psrt.srt.a + (uint32_t)uidx->psrt.idx.a[i]; srt_n = uidx->psrt.idx.a[i]>>33;
    // if(srt_n == 0) return;
    // integer_candidate(uidx, srt_a, srt_n, is_circle, buf);
    integer_candidate(uidx, buf, i, (asm_opt.purge_level_primary == 0?1:0));
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
    free(uidx->idx->ridx.idx.a); free(uidx->idx->ridx.occ.a); 
    memset(&(uidx->idx->ridx), 0, sizeof((uidx->idx->ridx)));
    filter_ul_ug(uidx->l1_ug);
	gen_ul_vec_rid_t(uidx->idx, NULL, uidx->l1_ug);
    update_ug_arch_ul_mul(uidx->l1_ug);

    free(uidx->pstr.idx.a); free(uidx->pstr.occ.a); free(uidx->pstr.str.a); 
    memset(&(uidx->pstr), 0, sizeof(uidx->pstr));
    init_ul_str_idx_t(uidx);
}

void ul_re_correct(ul_resolve_t *uidx, uint64_t n_r)
{   
    uint64_t k, occ, n_circle; uidx->str_b.n_thread = 1;
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
    }

    // kt_for(uidx->str_b.n_thread, worker_integer_postprecess, uidx, uidx->idx->n);
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

void ul_realignment_gfa(ug_opt_t *uopt, asg_t *sg)
{
    uint64_t i; uint8_t *r_het = NULL; bubble_type *bub = NULL;
    for (i = 0; i < sg->n_seq; ++i) {
        if(sg->seq[i].del) continue;
        sg->seq[i].c = PRIMARY_LABLE;
    }
    hic_clean(sg);
    ma_ug_t *init_ug = ul_realignment(uopt, sg);
    filter_sg_by_ug(sg, init_ug, uopt);
    print_ul_alignment(init_ug, &UL_INF, 47072, "after-0");
    bub = gen_bubble_chain(sg, init_ug, uopt, &r_het);
    print_ul_alignment(init_ug, &UL_INF, 47072, "after-1");
    ul_resolve_t *uidx = init_ul_resolve_t(sg, init_ug, bub, &UL_INF, r_het);
    print_ul_alignment(init_ug, &UL_INF, 47072, "after-2");
    ul_re_correct(uidx, 3);
    print_ul_alignment(init_ug, &UL_INF, 47072, "after-3");

    // print_debug_ul("UL.debug", init_ug, sg, uopt->coverage_cut, uopt->sources, uopt->ruIndex, bub, &UL_INF);

    // resolve_dip_bub_chains(uidx);

    // free(r_het); destory_bubbles(bub); free(bub);
}