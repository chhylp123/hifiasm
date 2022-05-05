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

#define generic_key(x) (x)
KRADIX_SORT_INIT(srt64, uint64_t, generic_key, 8)
#define OU_NOISY 2
#define ASG_ET_MERGEABLE 0
#define ASG_ET_TIP       1
#define ASG_ET_MULTI_OUT 2
#define ASG_ET_MULTI_NEI 3

typedef struct {
	asg_t *g; 
	ma_hit_t_alloc *src;
} sset_aux;

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

void fill_containment_by_ul(asg_t *g, ma_hit_t_alloc *src)
{

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

void ul_realignment_gfa(ug_opt_t *uopt, asg_t *sg)
{
    uint64_t i;
    for (i = 0; i < sg->n_seq; ++i) {
        if(sg->seq[i].del) continue;
        sg->seq[i].c = PRIMARY_LABLE;
    }
    hic_clean(sg);
    ul_realignment(uopt, sg);
    // if(ul_refine_alignment(uopt, sg)) update_sg_uo(sg, src);
}