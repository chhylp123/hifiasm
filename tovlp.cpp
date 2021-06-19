#define __STDC_LIMIT_MACROS
#include "float.h"
#include "horder.h"
#include <math.h>
#include "hic.h"
#include "htab.h"
#include "assert.h"
#include "Overlaps.h"
#include "Hash_Table.h"
#include "Correct.h"
#include "Purge_Dups.h"
#include "rcut.h"
#include "khashl.h"
#include "kthread.h"
#include "ksort.h"
#include "kseq.h" // FASTA/Q parser
#include "kdq.h"
#include "tovlp.h"
#include "hic.h"
KSEQ_INIT(gzFile, gzread)
KDQ_INIT(uint64_t)


typedef struct {
	kvec_t_u32_warp tt;
    kvec_t_u32_warp stack;
    uint8_t *vis;
    pdq pq_p;
    pdq pq_a;
    buf_t b;
    uint32_t min_v, offset1, offset2;
    long long min_d;
    int is_b;
} clean_t;

typedef struct {
	clean_t *a;
    size_t n, m;
    uint64_t max_dist;
    uint8_t *bs_flag, fp, fa;
    asg_t *g;
    ma_ug_t *ug;
} clean_mul_t;

clean_mul_t *init_clean_mul_t(asg_t *g, ma_ug_t *ug, uint64_t n_threads, uint8_t fp, uint8_t fa)
{
    uint32_t i;
    clean_t *kt = NULL;
    clean_mul_t *p = NULL; CALLOC(p, 1);
    p->g = g; p->ug = ug;
    CALLOC(p->bs_flag, p->g->n_seq<<1);
    p->n = p->m = n_threads;
    CALLOC(p->a, p->n);
    for (i = 0, kt = NULL; i < p->n; i++)
    {
        kt = &(p->a[i]);
        kv_init(kt->tt.a); kv_init(kt->stack.a);
        CALLOC(kt->vis, p->g->n_seq<<1);
        init_pdq(&kt->pq_p, p->g->n_seq<<1);
        init_pdq(&kt->pq_a, p->g->n_seq<<1);
        memset(&kt->b, 0, sizeof(kt->b));
        kt->b.a = (binfo_t*)calloc(p->g->n_seq<<1, sizeof(binfo_t));
    }
    p->max_dist = get_bub_pop_max_dist_advance(g, &kt->b);
    p->fp = fp; p->fa = fa;
    return p;
}

void destroy_clean_mul_t(clean_mul_t **p)
{
    uint32_t i;
    clean_t *kt = NULL;
    free((*p)->bs_flag);
    for (i = 0; i < (*p)->n; i++)
    {
        kt = &((*p)->a[i]);
        kv_destroy(kt->tt.a);
        kv_destroy(kt->stack.a);
        free(kt->vis);
        destory_pdq(&kt->pq_p);
        destory_pdq(&kt->pq_a);
        free(kt->b.a); free(kt->b.S.a); free(kt->b.T.a); free(kt->b.b.a); free(kt->b.e.a);
    }
    free((*p)->a);
    free((*p));
}

utg_trans_t *init_utg_trans_t(ma_ug_t *ug, ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, R_to_U* ruIndex, asg_t *read_g, int max_hang, int min_ovlp)
{
    utg_trans_t *o = NULL; 
    CALLOC(o, 1);
    o->reverse_sources = reverse_sources;
    o->coverage_cut = coverage_cut;
    o->ruIndex = ruIndex;
    o->read_g = read_g;
    o->rn = read_g->n_seq;
    memset(&(o->b0), 0, sizeof(o->b0));
    memset(&(o->b1), 0, sizeof(o->b1));
    MALLOC(o->pos_idx, o->rn); memset(o->pos_idx, -1, o->rn*sizeof(uint64_t));
    o->cug = copy_untig_graph(ug);
    o->max_hang = max_hang;
    o->min_ovlp = min_ovlp;
    return o;
}

void destroy_utg_trans_t(utg_trans_t **o)
{
    ;
}

uint64_t get_utg_chain_length(ma_ug_t *ug, uint32_t *a, uint32_t a_n)
{
    uint64_t i, k, l, qn, v, w, nv, offset;
    asg_arc_t *av  = NULL;
    for (i = 0, offset = 0; i < a_n; i++)
    {
        qn = a[i]>>1;
        l = ug->g->seq[qn].len;
        if(i + 1 < a_n)
        {
            v = a[i]; w = a[i+1];
            av = asg_arc_a(ug->g, v);
            nv = asg_arc_n(ug->g, v);
            for (k = 0; k < nv; k++)
            {
                if(av[k].del) continue;
                if(av[k].v == w) 
                {
                    l = asg_arc_len(av[k]);
                    break;
                }
            }
            if(k >= nv) fprintf(stderr, "ERROR-mc\n");
        }
        offset += l;
    }
    return offset;
}

void refine_utg_thit_t(utg_thit_t *q, kv_ca_buf_t* cb)
{
    ///already know [qScur, qEcur), [qSpre, qEpre)
    uint32_t s, e, i, si, ei;
    s = q->qs; e = q->qe;///[s, e)
    for (i = 0, si = ei = cb->n; i < cb->n; i++)
    {
        if(cb->a[i].c_x_p > s && si == cb->n)
        {
            si = i;
        } 

        if(cb->a[i].c_x_p > (e-1) && ei == cb->n)
        {
            ei = i;
        }

        if(si != cb->n && ei != cb->n) break;        
    }

    if(si == 0 || ei == 0) fprintf(stderr, "ERROR-si-ei-0\n");
    if(si >= cb->n || ei >= cb->n) fprintf(stderr, "ERROR-si-ei-1\n");
    si--; ei--;

    if(s < cb->a[si].c_x_p || ((si + 1) < cb->n && s >= cb->a[si + 1].c_x_p))
    {
        fprintf(stderr, "ERROR3\n");
    }
    if(e < cb->a[ei].c_x_p || ((ei + 1) < cb->n && e > cb->a[ei + 1].c_x_p))
    {
        fprintf(stderr, "ERROR4\n");
    }
    ///si and ei must be less than (cb->n-1)    
    q->ts = cb->a[si].c_y_p + 
        get_offset_adjust(s-cb->a[si].c_x_p, cb->a[si+1].c_x_p-cb->a[si].c_x_p, cb->a[si+1].c_y_p-cb->a[si].c_y_p);
 
    q->te = cb->a[ei].c_y_p + 
        get_offset_adjust(e-cb->a[ei].c_x_p, cb->a[ei+1].c_x_p-cb->a[ei].c_x_p, cb->a[ei+1].c_y_p-cb->a[ei].c_y_p);
    
    ///might be equal
    if(q->ts > q->te) fprintf(stderr, "ERROR5\n");
    // if(q->tScur >= q->tEcur)
    // {
    //     fprintf(stderr, "\n###q->tScur: %u, s: %u, si: %u, q->tEcur: %u, e: %u, ei: %u\n", 
    //     q->tScur, s, si, q->tEcur, e, ei);

    //     fprintf(stderr, "###cb->a[si].c_x_p: %u, cb->a[si].c_y_p: %u, cb->a[si+1].c_x_p: %u, cb->a[si+1].c_y_p: %u\n", 
    //     cb->a[si].c_x_p, cb->a[si].c_y_p, cb->a[si+1].c_x_p, cb->a[si+1].c_y_p);

    //     fprintf(stderr, "###cb->a[ei].c_x_p: %u, cb->a[ei].c_y_p: %u, cb->a[ei+1].c_x_p: %u, cb->a[ei+1].c_y_p: %u\n", 
    //     cb->a[ei].c_x_p, cb->a[ei].c_y_p, cb->a[ei+1].c_x_p, cb->a[ei+1].c_y_p);

    //     fprintf(stderr, "ERROR5\n");
    // }
}

///[ts, te)
void extract_sub_overlaps_utg_thit(uint32_t i_ts, uint32_t i_te, uint32_t i_tus, uint32_t i_tue,
uint32_t tn, kv_utg_thit_t_t* ktb, uint32_t bn)
{
    uint32_t i, ovlp, found, beg, end, offS, offE;
    utg_thit_t *q = NULL, x;
    for (i = found = 0; i < bn; i++)
    {
        q = &(ktb->a[i]);///for q, already know [qs, qe), [qus, que), [ts, te)

        ovlp = ((MIN(i_te, q->te) > MAX(i_ts, q->ts))?
                                        MIN(i_te, q->te) - MAX(i_ts, q->ts):0);
        if(found == 1 && ovlp == 0) break;
        if(ovlp > 0) found = 1;
        if(ovlp == 0) continue;

        
        beg = MAX(i_ts, q->ts); end = MIN(i_te, q->te);
        offS = beg - q->ts; offE = q->te - end;
        x.ts = q->ts + offS; 
        x.te = q->te - offE;
        //x.qs = q->qs + offS; 
        x.qs = q->qs + get_offset_adjust(offS, q->te-q->ts, q->qe-q->qs);        
        ///x.qe = q->qe - offE;
        x.qe = q->qe - get_offset_adjust(offE, q->te-q->ts, q->qe-q->qs);  

        x.qn = q->qn;
        offS = beg - q->ts; offE = q->te - end;
        if((x.qn&1) == 0)
        {
            // x.qus = q->qus + offS; 
            x.qus = q->qus + get_offset_adjust(offS, q->te-q->ts, q->que-q->qus);  
            // x.que = q->que - offE;
            x.que = q->que - get_offset_adjust(offE, q->te-q->ts, q->que-q->qus);
        }
        else
        {
            // x.qus = q->qus + offE; 
            x.qus = q->qus + get_offset_adjust(offE, q->te-q->ts, q->que-q->qus);
            // x.que = q->que - offS;
            x.que = q->que - get_offset_adjust(offS, q->te-q->ts, q->que-q->qus);
        }

        x.tn = tn; 
        offS = beg - i_ts; offE = i_te - end;
        if((x.tn&1) == 0)
        {
            // x.tus = i_tus + offS; 
            x.tus = i_tus +  get_offset_adjust(offS, i_te-i_ts, i_tue-i_tus); 
            // x.tue = i_tue - offE;
            x.tue = i_tue - get_offset_adjust(offE, i_te-i_ts, i_tue-i_tus); 
        }
        else
        {
            // x.tus = i_tus + offE; 
            x.tus = i_tus + get_offset_adjust(offE, i_te-i_ts, i_tue-i_tus); 
            // x.tue = i_tue - offS;
            x.tue = i_tue - get_offset_adjust(offS, i_te-i_ts, i_tue-i_tus);
        }
        
        kv_push(utg_thit_t, *ktb, x);

        // if(x.tus >= x.tue || x.qus >= x.que)
        // {
        //     fprintf(stderr, "\n*********x.qn: %u, x.tn: %u\n", x.qn, x.tn);
        //     fprintf(stderr, "x.qus: %u, x.que: %u, x.tus: %u, x.tue: %u\n", 
        //     x.qus, x.que, x.tus, x.tue);
        //     fprintf(stderr, "q->qs: %u, q->qe: %u, q->qus: %u, q->que: %u\n", 
        //     q->qs, q->qe, q->qus, q->que);
        //     fprintf(stderr, "q->ts: %u, q->te: %u, q->tus: %u, q->tue: %u\n", 
        //     q->ts, q->te, q->tus, q->tue);
        //     fprintf(stderr, "i_ts: %u, i_te: %u, i_tus: %u, i_tue: %u\n", 
        //     i_ts, i_te, i_tus, i_tue);
        // }
    }
}

uint32_t get_utg_trans_hit(utg_trans_hit_idx *t, utg_thit_t *hit)
{
    uint32_t r_beg, r_end, ovlp;
    uint32_t k, l, qn, v, w, nv;
    uint32_t *a = t->a, a_n = t->an;
    asg_arc_t *av  = NULL;
    hit->qs = hit->qe = hit->qn = hit->qus = hit->que = (uint32_t)-1;
    hit->ts = hit->te = hit->tn = hit->tus = hit->tue = (uint32_t)-1;
    
    while (t->ui < a_n)
    {
        qn = a[t->ui];
        l = t->ug->g->seq[qn>>1].len;
        if(t->ui + 1 < a_n)
        {
            v = a[t->ui]; w = a[t->ui+1];
            av = asg_arc_a(t->ug->g, v);
            nv = asg_arc_n(t->ug->g, v);
            for (k = 0; k < nv; k++)
            {
                if(av[k].del) continue;
                if(av[k].v == w) 
                {
                    l = asg_arc_len(av[k]);
                    break;
                }
            }
            if(k >= nv) fprintf(stderr, "ERROR-mc-1\n");
        }

        r_beg = t->len; r_end = t->len + t->ug->g->seq[qn>>1].len;
        t->len += l; t->ui++;

        ovlp = ((MIN(t->cEnd, r_end) > MAX(t->cBeg, r_beg))? (MIN(t->cEnd, r_end) - MAX(t->cBeg, r_beg)) : 0);
        if(ovlp == 0 && r_beg >= t->cEnd) return 0;
        if(ovlp == 0) continue;

        hit->qn = qn; hit->qs = r_beg; hit->qe = r_end; hit->qus = 0; hit->que = t->ug->g->seq[qn>>1].len;
        if(hit->qs < t->cBeg)
        {
            hit->qus += (t->cBeg - hit->qs);
            hit->qs = t->cBeg;
        } 
        
        if(hit->qe > t->cEnd)
        {
            hit->que -= (hit->qe - t->cEnd);
            hit->qe = t->cEnd;
        } 
        
        return 1;
    }
    return 0;
}

void reset_utg_trans_hit_idx(utg_trans_hit_idx *t, uint32_t* i_x_a, uint32_t i_x_n, ma_ug_t *i_ug, 
utg_trans_t *i_o, uint32_t i_cBeg, uint32_t i_cEnd)
{
    t->a = i_x_a;
    t->an = i_x_n;
    t->ug = i_ug;
    t->o = i_o;
    t->cBeg = i_cBeg;
    t->cEnd = i_cEnd;
    t->ui = t->len = 0;
}

void chain_trans_d(utg_trans_t *o, 
uint32_t *pri_a, uint32_t pri_n, uint32_t pri_beg, uint64_t *i_pri_len, 
uint32_t *aux_a, uint32_t aux_n, uint32_t aux_beg, uint64_t *i_aux_len,
ma_ug_t *ug, uint32_t flag, double score, const char* cmd)
{
    uint32_t i, len, bn;
    uint64_t pri_len, aux_len;
    kvec_asg_arc_t_offset* u_buffer = &(o->u_buffer);
    kvec_t_i32_warp* tailIndex = &(o->tailIndex);
    asg_arc_t_offset *tt = NULL;
    ca_buf_t *tx = NULL;

    if(i_pri_len) pri_len = (*i_pri_len);
    else pri_len = get_utg_chain_length(ug, pri_a, pri_n);
    
    if(i_aux_len) aux_len = (*i_aux_len);
    else aux_len = get_utg_chain_length(ug, aux_a, aux_n);
        
    tt = NULL; o->c_buf.n = 0; o->k_t_b.n = 0;
    if(tailIndex->a.n > 0) tt = &(u_buffer->a.a[tailIndex->a.a[0]]);
    if(!tt || (pri_beg < (tt->Off>>32) && aux_beg < ((uint32_t)tt->Off)))
    {
        kv_pushp(ca_buf_t, o->c_buf, &tx);
        tx->c_x_p = pri_beg; 
        tx->c_y_p = aux_beg;

        for (i = 0; i < tailIndex->a.n; i++)
        {
            tt = &(u_buffer->a.a[tailIndex->a.a[i]]);
            kv_pushp(ca_buf_t, o->c_buf, &tx);

            tx->c_x_p = tt->Off>>32;
            tx->c_y_p = (uint32_t)tt->Off;
        }
    }
    else if(tailIndex->a.n == 1)//1 ele in chain
    {
        kv_pushp(ca_buf_t, o->c_buf, &tx);
        tx->c_x_p = pri_beg; tx->c_y_p = aux_beg;
    }
    else if(tailIndex->a.n > 0)
    {
        uint32_t cx, cy, ax, ay, found = 0;
        for (i = 0; i < tailIndex->a.n; i++)
        {
            cx = cy = ax = ay = (uint32_t)-1;

            cx = u_buffer->a.a[tailIndex->a.a[i]].Off>>32;
            cy = (uint32_t)u_buffer->a.a[tailIndex->a.a[i]].Off;
            if((i + 1) < tailIndex->a.n)
            {
                ax = u_buffer->a.a[tailIndex->a.a[i+1]].Off>>32;
                ay = (uint32_t)u_buffer->a.a[tailIndex->a.a[i+1]].Off;
            }

            kv_pushp(ca_buf_t, o->c_buf, &tx);
            tx->c_x_p = cx; tx->c_y_p = cy;

            if(found) continue;
            if(pri_beg > cx && pri_beg < ax && aux_beg > cy && aux_beg < ay)
            {
                kv_pushp(ca_buf_t, o->c_buf, &tx);
                tx->c_x_p = pri_beg; tx->c_y_p = aux_beg;
                found = 1;
            } 
        }
    }

    
    // fprintf(stderr, "\ncmd-%s\n", cmd);
    // fprintf(stderr, "pri_beg=%u, pri_len=%lu\n", pri_beg, pri_len);
    // fprintf(stderr, "aux_beg=%u, aux_len=%lu\n", aux_beg, aux_len);

    // print_buf_t(ug, pri, "pri");
    // print_buf_t(ug, aux, "aux");


    tx = &(o->c_buf.a[o->c_buf.n-1]);
    len = MIN(pri_len - tx->c_x_p, aux_len - tx->c_y_p);
    if(len > 0)///insert boundary
    {
        kv_pushp(ca_buf_t, o->c_buf, &tx);
        tx->c_x_p = o->c_buf.a[o->c_buf.n-2].c_x_p + len;
        tx->c_y_p = o->c_buf.a[o->c_buf.n-2].c_y_p + len;     
    }

    tx = &(o->c_buf.a[0]);///insert boundary
    if(tx->c_x_p != 0 && tx->c_y_p != 0)///already at boundary
    {
        len = MIN(tx->c_x_p, tx->c_y_p);///offset
        kv_pushp(ca_buf_t, o->c_buf, &tx);
        for (i = 0; (i + 1)< o->c_buf.n; i++)
        {
            o->c_buf.a[o->c_buf.n - i - 1] = o->c_buf.a[o->c_buf.n - i - 2];
        }
        o->c_buf.a[0].c_x_p = o->c_buf.a[1].c_x_p - len;
        o->c_buf.a[0].c_y_p = o->c_buf.a[1].c_y_p - len;
    }

    ///chain is [s, e)
    if(o->c_buf.a[0].c_x_p != 0 && o->c_buf.a[0].c_y_p != 0) fprintf(stderr, "ERROR1\n");
    if(o->c_buf.a[o->c_buf.n-1].c_x_p!= pri_len && 
                                            o->c_buf.a[o->c_buf.n-1].c_y_p!= aux_len)
    {
        fprintf(stderr, "ERROR2\n");
    }

    utg_trans_hit_idx iter;
    utg_thit_t hit, *kh = NULL;
    ////////prx
    reset_utg_trans_hit_idx(&iter, pri_a, pri_n, ug, o, o->c_buf.a[0].c_x_p, o->c_buf.a[o->c_buf.n-1].c_x_p);    
    while(get_utg_trans_hit(&iter, &hit))//get [qs, qe), [qus, que)
    {
        refine_utg_thit_t(&hit, &(o->c_buf)); ///get [ts, te)
        kv_push(utg_thit_t, o->k_t_b, hit);
    }
    bn = o->k_t_b.n;
    
    ////////aux
    reset_utg_trans_hit_idx(&iter, aux_a, aux_n, ug, o, o->c_buf.a[0].c_y_p, o->c_buf.a[o->c_buf.n-1].c_y_p);    
    while(get_utg_trans_hit(&iter, &hit))
    {
        extract_sub_overlaps_utg_thit(hit.qs, hit.qe, hit.qus, hit.que, hit.qn, &(o->k_t_b), bn);
    }

    if(o->k_t_b.n - bn == 0) fprintf(stderr, "ERROR7\n");
    
    // fprintf(stderr, "\n******o->k_t_b.n: %u, bn: %u\n", (uint32_t)o->k_t_b.n, bn);
    // for (i = 0; i < pri_n; i++)
    // {
    //     fprintf(stderr, "p-utg%.6ul\n", (pri_a[i]>>1) + 1);
    // }
    // for (i = 0; i < aux_n; i++)
    // {
    //     fprintf(stderr, "a-utg%.6ul\n", (aux_a[i]>>1) + 1);
    // }

    u_trans_t *kt = NULL;
    double x_score, y_score;
    for (i = bn; i < o->k_t_b.n; i++)
    {
        kh = &(o->k_t_b.a[i]);
        if(kh->que <= kh->qus) continue;
        if(kh->tue <= kh->tus) continue;
        kv_pushp(u_trans_t, o->k_trans, &kt);
        kt->f = flag; kt->rev = ((kh->qn ^ kh->tn) & 1); kt->del = 0; 
        kt->qn = kh->qn>>1; kt->qs = kh->qus; kt->qe = kh->que; ///kt->qo = kh->qn&1;
        kt->tn = kh->tn>>1; kt->ts = kh->tus; kt->te = kh->tue; ///kt->to = kh->tn&1;
        if(score < 0)
        {
            kt->nw = (MIN((kt->qe - kt->qs), (kt->te - kt->ts)))*CHAIN_MATCH;
        } 
        else
        {
            x_score = ((double)(kt->qe-kt->qs)/(double)(o->c_buf.a[o->c_buf.n-1].c_x_p-o->c_buf.a[0].c_x_p))*score;
            y_score = ((double)(kt->te-kt->ts)/(double)(o->c_buf.a[o->c_buf.n-1].c_y_p-o->c_buf.a[0].c_y_p))*score;
            kt->nw = MIN(x_score, y_score);
        }
    }
}


void chain_bubble(buf_t *pri, uint64_t pri_len, buf_t* aux, uint64_t aux_len, 
uint32_t beg, uint32_t sink, ma_ug_t *ug, utg_trans_t *o)
{
    if(pri->b.n == 0 || aux->b.n == 0) return;
    asg_arc_t *av = NULL;
    uint32_t pri_v, aux_v, nv, i, priBeg, priEnd, auxBeg, auxEnd;

    priBeg = priEnd = auxBeg = auxEnd = (uint32_t)-1;
    pri_v = pri->b.a[0]; aux_v = aux->b.a[0];
    av = asg_arc_a(ug->g, beg);
    nv = asg_arc_n(ug->g, beg);
    for (i = 0; i < nv; ++i)
    {
        if(av[i].del) continue;
        if(av[i].v == pri_v) priBeg = av[i].ol;
        if(av[i].v == aux_v) auxBeg = av[i].ol;
    }

    pri_v = pri->b.a[pri->b.n-1]^1; aux_v = aux->b.a[aux->b.n-1]^1;
    av = asg_arc_a(ug->g, sink);
    nv = asg_arc_n(ug->g, sink);
    for (i = 0; i < nv; ++i)
    {
        if(av[i].del) continue;
        if(av[i].v == pri_v) priEnd = ((pri_len > av[i].ol)? (pri_len - av[i].ol - 1) : 0);
        if(av[i].v == aux_v) auxEnd = ((aux_len > av[i].ol)? (aux_len - av[i].ol - 1) : 0);
    }

    if(priBeg == (uint32_t)-1 || priEnd == (uint32_t)-1 || auxBeg == (uint32_t)-1 || auxEnd == (uint32_t)-1)
    {
        fprintf(stderr, "ERROR-s_bubble\n");
    }

    o->u_buffer.a.n = o->tailIndex.a.n = 0;

    kv_resize(asg_arc_t_offset, o->u_buffer.a, 1);
    o->u_buffer.a.n = 1;
    o->u_buffer.a.a[0].Off = priEnd;
    o->u_buffer.a.a[0].Off <<= 32;
    o->u_buffer.a.a[0].Off |= auxEnd;

    kv_resize(int32_t, o->tailIndex.a, 1);
    o->tailIndex.a.n = 1;
    o->tailIndex.a.a[0] = 0;

    chain_trans_d(o, pri->b.a, pri->b.n, priBeg, &pri_len, aux->b.a, aux->b.n, auxBeg, &aux_len, ug, RC_0, -1024, __func__);
}

void chain_c_bubble(uint32_t query, buf_t *target, buf_t *idx, ma_ug_t *ug, utg_trans_t *o)
{
    if(target->b.n == 0) return;
    uint32_t qs, qe, ts, te, i, v, ovlp;
    qs = idx->a[query].d; qe = qs + ug->g->seq[query>>1].len;
    uint64_t qlen = ug->g->seq[query>>1].len, tlen;
    o->u_buffer.a.n = o->tailIndex.a.n = 0;

    for (i = 0; i < target->b.n; ++i)
    {
        v = target->b.a[i];
        if(v < query) continue; //avoid dup
        ts = idx->a[v].d; te = ts + ug->g->seq[v>>1].len; tlen = ug->g->seq[v>>1].len;
        ovlp = ((MIN(qe, te) > MAX(qs, ts))? (MIN(qe, te) - MAX(qs, ts)) : 0);
        if(ovlp == 0) continue;
        chain_trans_d(o, &query, 1, MAX(qs, ts) - qs, &qlen, 
                                                    &v, 1, MAX(qs, ts) - ts, &tlen, ug, RC_0, -1024, __func__);
    }
    // fprintf(stderr, "\nocc: %u\n", (uint32_t)(cov->t_ch->k_trans.n - i_n));
    // for (i = i_n; i < cov->t_ch->k_trans.n; i++)
    // {
    //     fprintf(stderr, "s-utg%.6ul\t%u\t%u\td-utg%.6ul\t%u\t%u\trev(%u)\n", 
    //     cov->t_ch->k_trans.a[i].qn+1, cov->t_ch->k_trans.a[i].qs, cov->t_ch->k_trans.a[i].qe, 
    //     cov->t_ch->k_trans.a[i].tn+1, cov->t_ch->k_trans.a[i].ts, cov->t_ch->k_trans.a[i].te, 
    //     cov->t_ch->k_trans.a[i].rev);
    // }
}

void tpSort(asg_t *g, utg_trans_t *o, uint32_t beg, uint32_t sink)
{
    buf_t *b = &(o->b0);
    uint64_t *visited = o->pos_idx;
    uint32_t v = beg, nv, kv, i;
    asg_arc_t *av = NULL;

    b->b.n = 0; o->topo_res.n = 0;
    kv_push(uint32_t, b->b, v);
    while (b->b.n > 0)
    {
        ///b->b.n--;
        v = b->b.a[b->b.n-1];
        if(visited[v>>1] == (uint64_t)-1)
        {
            visited[v>>1] = 0;
        } 
        
        nv = asg_arc_n(g, v);
        av = asg_arc_a(g, v);
        for (i = kv = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            if((av[i].v>>1) == (beg>>1) || (av[i].v>>1) == (sink>>1)) continue;
            if(visited[av[i].v>>1] != (uint64_t)-1) continue;
            kv_push(uint32_t, b->b, av[i].v);
            kv++;
        }

        if(kv != 0) continue;
        b->b.n--;
        if(visited[v>>1] != 1)
        {
            kv_push(uint32_t, o->topo_res, v);
            visited[v>>1] = 1;
        }
    }
    for (i = 0; i < o->topo_res.n; ++i) 
    {
        visited[o->topo_res.a[i]>>1] = (uint64_t)-1;
    }

    o->topo_res.n--;//remove beg
    for (i = 0; i < (o->topo_res.n>>1); ++i) 
    {
        v = o->topo_res.a[i];
        o->topo_res.a[i] = o->topo_res.a[o->topo_res.n - i - 1];
        o->topo_res.a[o->topo_res.n - i - 1] = v;
    }

    /*******************************for debug************************************/
    // for (i = 0; i < cov->n; ++i) 
    // {
    //     if(visited[i] != (uint64_t)-1) fprintf(stderr, "ERROR-2\n");
    // }
    // debug_topo_sorting(g, cov, beg, sink);
    /*******************************for debug************************************/
}

void dfs_trans_chain(asg_t *g, utg_trans_t *o, uint32_t v, uint32_t beg, uint32_t sink)
{
    buf_t *b = &(o->b0);
    b->b.n = 0;
    if(v == beg || v == sink) return;
    uint64_t *flag = o->pos_idx;
    asg_arc_t *acur = NULL;
    uint32_t cur, ncur, i;
    v = v << 1;
    kv_push(uint32_t, b->b, v);
    while (b->b.n > 0)
    {
        b->b.n--;
        cur = b->b.a[b->b.n];
        if(flag[cur>>1] == 0 && (cur>>1) != (v>>1)) continue;
        flag[cur>>1] = 0;

        ncur = asg_arc_n(g, cur);
        acur = asg_arc_a(g, cur);
        for (i = 0; i < ncur; i++)
        {
            if(acur[i].del) continue;
            if((acur[i].v>>1) == beg || (acur[i].v>>1) == sink) continue;
            if(flag[acur[i].v>>1] == 0) continue;
            kv_push(uint32_t, b->b, acur[i].v);
        }
    }

    v = v + 1;
    kv_push(uint32_t, b->b, v);
    while (b->b.n > 0)
    {
        b->b.n--;
        cur = b->b.a[b->b.n];
        if(flag[cur>>1] == 0 && (cur>>1) != (v>>1)) continue;
        flag[cur>>1] = 0;

        ncur = asg_arc_n(g, cur);
        acur = asg_arc_a(g, cur);
        for (i = 0; i < ncur; i++)
        {
            if(acur[i].del) continue;
            if((acur[i].v>>1) == beg || (acur[i].v>>1) == sink) continue;
            if(flag[acur[i].v>>1] == 0) continue;
            kv_push(uint32_t, b->b, acur[i].v);
        }
    }

    b->b.n = 0;
    for (i = 0; i < o->topo_res.n; ++i) //has been sorted
    {
        if(flag[o->topo_res.a[i]>>1] == 0)
        {
            flag[o->topo_res.a[i]>>1] = (uint64_t)-1;  
        } 
        else
        {
            kv_push(uint32_t, b->b, o->topo_res.a[i]);
        }
    }


    /*******************************for debug************************************/
    // for (i = 0; i < cov->n; ++i) 
    // {
    //     if(flag[i] != (uint64_t)-1) fprintf(stderr, "ERROR-0\n");
    // }
    /*******************************for debug************************************/
}

void asg_bub_collect_ovlp(ma_ug_t *ug, uint32_t v0, buf_t *b, utg_trans_t *o)
{
    uint32_t i, uId;
    ///b->S.a[0] is the sink of this bubble
    if(get_real_length(ug->g, v0, NULL) == 2 && get_real_length(ug->g, b->S.a[0]^1, NULL) == 2)
    {
        long long tmp, max_stop_nodeLen, max_stop_baseLen, bch_occ[2], len[2];
        uint32_t bch[2], convex[2];
        get_real_length(ug->g, v0, bch);

        ///in rare cases, one side of a bubble might be empty
        if((bch[0]>>1)!=(b->S.a[0]>>1) && (bch[1]>>1)!=(b->S.a[0]>>1))
        {
            get_unitig(ug->g, NULL, bch[0], &convex[0], &bch_occ[0], &tmp, 
                                                &max_stop_nodeLen, &max_stop_baseLen, 1, NULL);
            get_unitig(ug->g, NULL, bch[1], &convex[1], &bch_occ[1], &tmp, 
                                                &max_stop_nodeLen, &max_stop_baseLen, 1, NULL);
            if(((bch_occ[0] + bch_occ[1] + 1) == (uint32_t)b->b.n) && 
                    get_real_length(ug->g, convex[0], NULL) == 1 && get_real_length(ug->g, convex[1], NULL) == 1)
            {
                get_real_length(ug->g, convex[0], &convex[0]);
                get_real_length(ug->g, convex[1], &convex[1]);
                if(convex[0] == b->S.a[0] && convex[1] == b->S.a[0])
                {
                    o->b0.b.n = 0;
                    get_unitig(ug->g, NULL, bch[0], &convex[0], &bch_occ[0], &len[0], &max_stop_nodeLen, &max_stop_baseLen, 1, &(o->b0));
                    o->b1.b.n = 0;
                    get_unitig(ug->g, NULL, bch[1], &convex[1], &bch_occ[1], &len[1], &max_stop_nodeLen, &max_stop_baseLen, 1, &(o->b1));
                    
                    chain_bubble(&(o->b0), len[0], &(o->b1), len[1], v0, b->S.a[0]^1, ug, o);
                    return;
                }
            }
        }
    }

    tpSort(ug->g, o, v0, b->S.a[0]);
    ///if(o->topo_res.n != b->b.n - 1) fprintf(stderr, "ERROR-4\n");
    if(o->topo_res.n == 0) return;
    for (i = 0; i < o->topo_res.n; ++i)
    {
        uId = o->topo_res.a[i]>>1;
        if(uId == (b->S.a[0]>>1)) continue;

        dfs_trans_chain(ug->g, o, uId, v0>>1, b->S.a[0]>>1);

        if(o->b0.b.n == 0) continue;
        chain_c_bubble(o->topo_res.a[i], &(o->b0), b, ug, o);
    }
}


void collect_trans_ovlp(const char* cmd, buf_t* pri, uint64_t pri_offset, buf_t* aux, uint64_t aux_offset,
ma_ug_t *ug, utg_trans_t *o)
{
    uint32_t thre_pri;
    uint64_t len_aux;
    if(pri->b.n == 0 || aux->b.n == 0) return;

    len_aux = set_utg_offset(aux->b.a, aux->b.n, ug, o->read_g, o->pos_idx, 0, 0);
    chain_trans_ovlp(NULL, o, ug, o->read_g, pri, len_aux, &thre_pri);

    
    if(thre_pri > 0)
    {
        chain_trans_d(o, pri->b.a, pri->b.n, pri_offset, NULL, aux->b.a, aux->b.n, aux_offset, &len_aux, 
        ug, RC_1, -1024, __func__);
        // fprintf(stderr, "\n%s, thre_pri: %u, len_aux: %lu\n", cmd, thre_pri, len_aux);
        // print_buf_t(ug, pri, "pri");
        // print_buf_t(ug, aux, "aux");
    }
    set_utg_offset(aux->b.a, aux->b.n, ug, o->read_g, o->pos_idx, 1, 0);
}


uint32_t dfs_set(asg_t *g, uint32_t v0, uint32_t fbv, kvec_t_u32_warp *stack, kvec_t_u32_warp *tmp, uint8_t *vis, uint8_t flag)
{
    uint32_t cur, ncur, i, occ = 0;;
    asg_arc_t *acur = NULL;
    stack->a.n = 0;
    kv_push(uint32_t, stack->a, v0);
    // vis[v0] = 0;

    while (stack->a.n > 0)
    {
        stack->a.n--;
        cur = stack->a.a[stack->a.n];
        if(vis[cur] && (!(vis[cur]&flag)))
        {
            vis[cur] |= flag;
            kv_push(uint32_t, tmp->a, cur);
            occ++;
        }
        if(vis[cur]) continue;
        else kv_push(uint32_t, tmp->a, cur);

        vis[cur] |= flag;
        ncur = asg_arc_n(g, cur);
        acur = asg_arc_a(g, cur);
        for (i = 0; i < ncur; i++)
        {
            if(acur[i].del) continue;
            if((acur[i].v>>1) == fbv) continue;
            if(vis[acur[i].v] && (!(vis[acur[i].v]&flag)))
            {
                vis[acur[i].v] |= flag;
                kv_push(uint32_t, tmp->a, acur[i].v);
                occ++;
                continue;
            } 
            
            kv_push(uint32_t, stack->a, acur[i].v);
        }
    }

    return occ;
}


uint32_t dfs_reach(asg_t *g, uint32_t src, uint32_t dest, kvec_t_u32_warp *stack, kvec_t_u32_warp *tmp, uint8_t *vis)
{
    uint32_t cur, ncur, i, occ = 0;;
    asg_arc_t *acur = NULL;
    stack->a.n = 0; tmp->a.n = 0;
    kv_push(uint32_t, stack->a, src);

    while (stack->a.n > 0)
    {
        stack->a.n--;
        cur = stack->a.a[stack->a.n];
        if(vis[cur]) continue;
        else kv_push(uint32_t, tmp->a, cur);
        vis[cur] = 1;
        if(cur == dest)
        {
            occ = 1;
            break;
        }


        ncur = asg_arc_n(g, cur);
        acur = asg_arc_a(g, cur);
        for (i = 0; i < ncur; i++)
        {
            if(acur[i].del) continue;
            if(vis[acur[i].v]) continue;
            kv_push(uint32_t, stack->a, acur[i].v);
            if(acur[i].v == dest)
            {
                occ = 1;
                break;
            }
        }
        if(occ) break;
    }

    for (i = 0; i < tmp->a.n; i++)
    {
        vis[tmp->a.a[i]] = 0;
    }
    
    stack->a.n = 0; tmp->a.n = 0;
    return occ;
}

int is_hap_ovlp(ma_ug_t *ug, asg_t *read_sg, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex, 
 uint32_t *a, uint32_t a_n, uint32_t *b, uint32_t b_n)
{
    uint32_t i, k, j, qn, tn, is_Unitig, uId, found = 0;
    ma_utg_t *u = NULL;
    for (i = 0; i < b_n; i++)
    {
        u = &(ug->u.a[b[i]]);
        for (k = 0; k < u->n; k++)
        {
            qn = (u->a[k]>>33);
            set_R_to_U(ruIndex, qn, b[i], 1, NULL);
        }
    }


    for (i = 0; i < a_n; i++)
    {
        u = &(ug->u.a[a[i]]);
        for (k = 0; k < u->n; k++)
        {
            qn = (u->a[k]>>33);
            for (j = 0; j < reverse_sources[qn].length; j++)
            {
                tn = Get_tn(reverse_sources[qn].buffer[j]);
                if(read_sg->seq[tn].del == 1)
                {
                    get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
                    if(tn == (uint32_t)-1 || is_Unitig == 1 || read_sg->seq[tn].del == 1) continue;
                }

                get_R_to_U(ruIndex, tn, &uId, &is_Unitig);
                if(uId!=(uint32_t)-1 && is_Unitig == 1)
                {
                    found = 1;
                    break;
                }
            }
            if(found) break;
        }
        if(found) break;
    }




    for (i = 0; i < b_n; i++)
    {
        u = &(ug->u.a[b[i]]);
        for (k = 0; k < u->n; k++)
        {
            qn = (u->a[k]>>33);
            ruIndex->index[qn] = (uint32_t)-1;
        }
    }

    return found;
}

uint32_t get_dir_v(uint32_t x, buf_t *b)
{
    uint32_t i;
    for (i = 0; i < b->b.n; i++)
    {
        if((b->b.a[i]>>1) == x) return b->b.a[i];
    }

    return (uint32_t)-1;
}

int is_pop_unitig(uint32_t v, asg_t *g, ma_ug_t *ug, kvec_t_u32_warp *tt, kvec_t_u32_warp *stack, uint8_t *vis,
pdq *pq_p, pdq *pq_a, uint8_t fp, uint8_t fa, long long *d)
{
    asg_arc_t *as = NULL;
    uint32_t s, sv, return_flag, convex, ns, i, k, nc, p_n, found, *a_a = NULL, a_n;
    long long ll, tmp, max_stop_nodeLen, max_stop_baseLen;
    if(g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) return 0;
    if(asg_arc_n(g, v) == 0 || get_real_length(g, v, NULL) != 1) return 0;
    get_real_length(g, v, &s);
    if(get_real_length(g, s^1, NULL) < 2) return 0;
    return_flag = get_unitig(g, ug, v^1, &convex, &ll, &tmp, &max_stop_nodeLen,&max_stop_baseLen, 1, NULL);
    if(return_flag == LOOP) return 0;

    as = asg_arc_a(g, convex); ns = asg_arc_n(g, convex);
    for (i = nc = 0; i < ns; i++)
    {
        if(as[i].del) continue;
        if(get_real_length(g, as[i].v^1, NULL) < 2) break;
        nc++;
    }
    if(nc == 0 || i < ns) return 0;

    tt->a.n = 0; s^=1; sv = v^1;
    dfs_set(g, sv, s>>1, stack, tt, vis, fp);
    p_n = tt->a.n;
    as = asg_arc_a(g, s); ns = asg_arc_n(g, s); found = 0;
    for (i = 0; i < ns; i++)
    {
        if(as[i].del || as[i].v == sv) continue;
        nc = dfs_set(g, as[i].v, s>>1, stack, tt, vis, fa);
        if(nc && check_trans_relation_by_path(sv, as[i].v, pq_p, NULL, NULL, pq_a, NULL, NULL, g, 
        vis, fp+fa, nc, 0.45, d))
        {
            found = 1;
        }
        a_a = tt->a.a + p_n; a_n = tt->a.n - p_n;
        for (k = 0; k < a_n; k++)
        {
            if(vis[a_a[k]]&fa) vis[a_a[k]] -= fa;
        }
        tt->a.n = p_n;
        if(found) break;
    }
    a_a = tt->a.a; a_n = p_n;
    for (k = 0; k < a_n; k++) vis[a_a[k]] = 0;
    return found;
}

static void unitig_iden_worker(void *_data, long eid, int tid)
{
    long long d;
    clean_mul_t *buf = (clean_mul_t*)_data;
    clean_t *b = &(buf->a[tid]);
    if(is_pop_unitig(eid, buf->g, buf->ug, &(b->tt), &(b->stack), b->vis, &(b->pq_p), &(b->pq_a), buf->fp, buf->fa, &d))
    {
        if(b->min_v == (uint32_t)-1 || (b->min_d > d) || 
        (b->min_d == d && buf->g->seq[b->min_v>>1].len > buf->g->seq[eid>>1].len) ||
        (b->min_d == d && buf->g->seq[b->min_v>>1].len == buf->g->seq[eid>>1].len && b->min_v > eid))
        {
            b->min_v = eid; b->min_d = d;
        }
    }
}

int is_pop_bub(uint32_t v, asg_t *g, ma_ug_t *ug, buf_t *b, uint64_t max_dist, uint8_t *bs_flag)
{
    bs_flag[v] = 4;
    asg_arc_t *av = NULL;
    uint32_t nv, i, n_arc;
    int occ = 0;
    nv = asg_arc_n(g, v);
    av = asg_arc_a(g, v);
    ///some node could be deleted
    if (nv < 2 || g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) return occ;
    ///some edges could be deleted
    for (i = n_arc = 0; i < nv; ++i) // asg_bub_pop1() may delete some edges/arcs
        if (!av[i].del) ++n_arc;
    if (n_arc < 2) return occ;
    if(asg_bub_pop1_primary_trio(ug->g, NULL, v, max_dist, b, (uint32_t)-1, (uint32_t)-1, 0, NULL, NULL, NULL, 0, 0, NULL))
    {
        //beg is v, end is b.S.a[0]
        //note b.b include end, does not include beg
        bs_flag[v] = 0;
        occ = 1;
    }
    return occ;
}

static void bub_iden_worker(void *_data, long eid, int tid)
{
    clean_mul_t *buf = (clean_mul_t*)_data;
    clean_t *b = &(buf->a[tid]);
    if(is_pop_bub(eid, buf->g, buf->ug, &(b->b), buf->max_dist, buf->bs_flag))
    {
        b->is_b = 1;
    }
}

int asg_pop_bubble_primary_trio(asg_t *g, ma_ug_t *ug, uint8_t *bs_flag, buf_t *b,
uint64_t max_dist, uint32_t positive_flag, uint32_t negative_flag, utg_trans_t *o)
{
	uint32_t v, n_vtx = g->n_seq * 2, n_arc, nv, i;
	uint64_t n_pop = 0;
    asg_arc_t *av = NULL;

    if(max_dist > 0)
    {
        for (v = 0; v < n_vtx; ++v) 
        {
            if(bs_flag[v] != 0) continue;
            nv = asg_arc_n(g, v);
            av = asg_arc_a(g, v);
            ///some node could be deleted
            if (nv < 2 || g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
            ///some edges could be deleted
            for (i = n_arc = 0; i < nv; ++i) // asg_bub_pop1() may delete some edges/arcs
                if (!av[i].del) ++n_arc;
            if (n_arc < 2) continue;
            if(asg_bub_pop1_primary_trio(ug->g, NULL, v, max_dist, b, (uint32_t)-1, (uint32_t)-1, 0, NULL, NULL, NULL, 0, 0, NULL))
            {
                //beg is v, end is b.S.a[0]
                //note b.b include end, does not include beg
                for (i = 0; i < b->b.n; i++)
                {
                    if(b->b.a[i]==v || b->b.a[i]==b->S.a[0]) continue;
                    bs_flag[b->b.a[i]] = bs_flag[b->b.a[i]^1] = 1;
                }
                bs_flag[v] = 2; bs_flag[b->S.a[0]^1] = 3;
            }
        }

        //traverse all node with two directions 
        for (v = 0; v < n_vtx; ++v) {
            if(bs_flag[v] !=2) continue;
            nv = asg_arc_n(g, v);
            av = asg_arc_a(g, v);
            ///some node could be deleted
            if (nv < 2 || g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
            ///some edges could be deleted
            for (i = n_arc = 0; i < nv; ++i) // asg_bub_pop1() may delete some edges/arcs
                if (!av[i].del) ++n_arc;
            if (n_arc > 1)
                n_pop += asg_bub_pop1_primary_trio(ug->g, ug, v, max_dist, b, positive_flag, negative_flag, 1, NULL, NULL, NULL, 0, 0, o);
        }
        
        if(VERBOSE >= 1)
        {
            fprintf(stderr, "[M::%s] popped %lu bubbles\n", __func__, (unsigned long)n_pop);
        }
    }
    
    if (n_pop) asg_cleanup(g);
    return n_pop;
}

int get_min_dec(clean_mul_t *cl, uint32_t positive_flag, uint32_t negative_flag, utg_trans_t *o, uint32_t *v)
{
    uint32_t i, n_pop;
    clean_t *p = NULL;
    (*v) = (uint32_t)-1;
    for (i = 0; i < cl->n; i++) cl->a[i].is_b = 0, cl->a[i].min_v = (uint32_t)-1;    
    kt_for(cl->n, bub_iden_worker, cl, cl->g->n_seq<<1);
    for (i = 0; i < cl->n; i++)
    {
        if(cl->a[i].is_b) ///pop bubble
        {
            n_pop = asg_pop_bubble_primary_trio(cl->g, cl->ug, cl->bs_flag, &(cl->a[i].b),
            cl->max_dist, positive_flag, negative_flag, o);
            if(n_pop ==0) fprintf(stderr, "ERROR-n_pop\n");
            break;
        } 
    }
    kt_for(cl->n, unitig_iden_worker, cl, cl->g->n_seq<<1);
    for (i = 0; i < cl->n; i++)
    {
        if(cl->a[i].min_v == (uint32_t)-1) continue;
        // if(b->min_v == (uint32_t)-1 || (b->min_d > d) || (b->min_d == d && b->min_v < eid))
        if(!p || p->min_d > cl->a[i].min_d || 
        (p->min_d == cl->a[i].min_d && cl->g->seq[p->min_v>>1].len > cl->g->seq[cl->a[i].min_v>>1].len) ||
        (p->min_d == cl->a[i].min_d && cl->g->seq[p->min_v>>1].len == cl->g->seq[cl->a[i].min_v>>1].len && p->min_v > cl->a[i].min_v))
        {
            p = &(cl->a[i]);
        }
    }
    if(p) (*v) = p->min_v;
    return p?1:0;
}

int asg_arc_decompress_mul(asg_t *g, ma_ug_t *ug, asg_t *read_sg, uint32_t positive_flag, uint32_t negative_flag,
ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex, utg_trans_t *o)
{
    double startTime = Get_T();
    asg_arc_t *as = NULL, *pm = NULL;
    u_trans_t *kt = NULL;
    uint32_t v, s, sv, nc, ns, n_vtx = g->n_seq * 2, n_reduced = 0, convex, pi, qn, tn;
    uint32_t return_flag, k, m, i, p_n, a_n, *a_a = NULL;
    uint8_t fp = 1, fa = 2, found;
    long long ll, tmp, max_stop_nodeLen, max_stop_baseLen;
    kvec_t_u32_warp tt; kv_init(tt.a);
    kvec_t_u32_warp stack; kv_init(stack.a);
    uint8_t *vis = NULL; CALLOC(vis, n_vtx);
    buf_t b;
    memset(&b, 0, sizeof(buf_t));
    pdq pq_p, pq_a;
    init_pdq(&pq_p, g->n_seq<<1);
    init_pdq(&pq_a, g->n_seq<<1);
    uint32_t *path_p = NULL, *path_q = NULL; 
    CALLOC(path_p, g->n_seq<<1);
    CALLOC(path_q, g->n_seq<<1);
    clean_mul_t *cl = init_clean_mul_t(g, ug, asm_opt.thread_num, fp, fa);
    
    while (get_min_dec(cl, positive_flag, negative_flag, o, &v))
    {
        if(g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        if(asg_arc_n(g, v) == 0 || get_real_length(g, v, NULL) != 1) continue;
        get_real_length(g, v, &s);
        if(get_real_length(g, s^1, NULL) < 2) continue;
        return_flag = get_unitig(g, ug, v^1, &convex, &ll, &tmp, &max_stop_nodeLen,&max_stop_baseLen, 1, NULL);
        if(return_flag == LOOP) continue;

        as = asg_arc_a(g, convex); ns = asg_arc_n(g, convex);
        for (i = nc = 0; i < ns; i++)
        {
            if(as[i].del) continue;
            if(get_real_length(g, as[i].v^1, NULL) < 2) break;
            nc++;
        }
        if(nc == 0 || i < ns) continue;

        tt.a.n = 0; s^=1; sv = v^1;
        dfs_set(g, sv, s>>1, &stack, &tt, vis, fp);
        p_n = tt.a.n;
        as = asg_arc_a(g, s); ns = asg_arc_n(g, s); found = 0;
        for (i = 0, pm = NULL; i < ns; i++)
        {
            if(as[i].del || as[i].v != sv) continue;
            pm = &(as[i]);
            break;
        }
 
        for (i = 0; i < ns; i++)
        {
            if(as[i].del || as[i].v == sv) continue;
            nc = dfs_set(g, as[i].v, s>>1, &stack, &tt, vis, fa);
            if(nc && check_trans_relation_by_path(sv, as[i].v, &pq_p, path_p, &(o->b0), &pq_a, path_q, &(o->b1), g, 
            vis, fp+fa, nc, 0.45, NULL))
            {
                found = 1;
            }
            a_a = tt.a.a + p_n; a_n = tt.a.n - p_n;
            for (k = 0; k < a_n; k++)
            {
                if(vis[a_a[k]]&fa) vis[a_a[k]] -= fa;
            }
            tt.a.n = p_n;
            if(found) break;
        }
        a_a = tt.a.a; a_n = p_n;
        for (k = 0; k < a_n; k++) vis[a_a[k]] = 0;

        if(found == 0) fprintf(stderr, "ERROR-found\n");
        if(found)
        {
            b.b.n = 0;
            get_unitig(g, ug, sv, &convex, &ll, &tmp, &max_stop_nodeLen, &max_stop_baseLen, 1, &b);
            for (k = 0; k < b.b.n; k++)
            {
                g->seq[b.b.a[k]>>1].c = ALTER_LABLE;
            }

            for (k = 0; k < b.b.n; k++)
            {
                asg_seq_drop(g, b.b.a[k]>>1);
            }

            // uint32_t ui;
            // fprintf(stderr, "\n");
            // for (ui = 0; ui < o->b0.b.n; ui++)
            // {
            //     fprintf(stderr, "p-utg%.6ul\n", (o->b0.b.a[ui]>>1)+1);
            // }
            // for (ui = 0; ui < o->b1.b.n; ui++)
            // {
            //     fprintf(stderr, "a-utg%.6ul\n", (o->b1.b.a[ui]>>1)+1);
            // }


            if(o->b0.b.n < b.b.n)  fprintf(stderr, "ERROR-ll\n"); 
            o->b0.b.n = b.b.n; pi = o->k_trans.n;
            collect_trans_ovlp(__func__, &(o->b0), pm->ol, &(o->b1), as[i].ol, ug, o);

            for (k = m = pi; k < o->k_trans.n; k++)
            {
                kt = &(o->k_trans.a[k]);
                if(is_hap_ovlp(ug, read_sg, reverse_sources, ruIndex, &(kt->qn), 1, &(kt->tn), 1) || 
                      is_hap_ovlp(ug, read_sg, reverse_sources, ruIndex, &(kt->tn), 1, &(kt->qn), 1))
                {
                    qn = get_dir_v(kt->qn, &(o->b0));
                    tn = get_dir_v(kt->tn, &(o->b1));

                    if(dfs_reach(o->cug->g, qn, tn, &stack, &tt, vis) || dfs_reach(o->cug->g, tn, qn, &stack, &tt, vis))
                    {
                        // fprintf(stderr, "<<<<<<DEL: q-utg%.6ul (qs: %u, qe: %u), t-utg%.6ul (ts: %u, te: %u)\n", 
                        // kt->qn + 1, kt->qs, kt->qe, kt->tn + 1, kt->ts, kt->te);
                        continue;
                    }

                    // fprintf(stderr, "q-utg%.6ul (qs: %u, qe: %u), t-utg%.6ul (ts: %u, te: %u)\n", 
                    // kt->qn + 1, kt->qs, kt->qe, kt->tn + 1, kt->ts, kt->te);
                    o->k_trans.a[m] = *kt;
                    m++;
                }
                // else
                // {
                //     fprintf(stderr, "******delete: q-utg%.6ul (qs: %u, qe: %u), t-utg%.6ul (ts: %u, te: %u)\n", 
                //     kt->qn + 1, kt->qs, kt->qe, kt->tn + 1, kt->ts, kt->te);
                // }      
            }
            o->k_trans.n = m;
        }
    }
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
    }
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    kv_destroy(tt.a);
    kv_destroy(stack.a);
    free(vis);
    destory_pdq(&pq_p);
    destory_pdq(&pq_a);
    free(path_p);
    free(path_q);
    destroy_clean_mul_t(&cl);
    return n_reduced;
}


int asg_arc_decompress(asg_t *g, ma_ug_t *ug, asg_t *read_sg, ma_hit_t_alloc* reverse_sources,
R_to_U* ruIndex, utg_trans_t *o)
{
    double startTime = Get_T();
    asg_arc_t *as = NULL, *pm = NULL;
    u_trans_t *kt = NULL;
    uint32_t v, s, sv, nc, ns, n_vtx = g->n_seq * 2, n_reduced = 0, convex, pi, qn, tn;
    uint32_t return_flag, k, m, i, p_n, a_n, *a_a = NULL;
    uint8_t fp = 1, fa = 2, found;
    long long ll, tmp, max_stop_nodeLen, max_stop_baseLen;
    kvec_t_u32_warp tt; kv_init(tt.a);
    kvec_t_u32_warp stack; kv_init(stack.a);
    uint8_t *vis = NULL; CALLOC(vis, n_vtx);
    buf_t b;
    memset(&b, 0, sizeof(buf_t));
    pdq pq_p, pq_a;
    init_pdq(&pq_p, g->n_seq<<1);
    init_pdq(&pq_a, g->n_seq<<1);
    uint32_t *path_p = NULL, *path_q = NULL; 
    CALLOC(path_p, g->n_seq<<1);
    CALLOC(path_q, g->n_seq<<1);
    
    for (v = 0; v < n_vtx; v++)
    {
        if(g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        if(asg_arc_n(g, v) == 0 || get_real_length(g, v, NULL) != 1) continue;
        get_real_length(g, v, &s);
        if(get_real_length(g, s^1, NULL) < 2) continue;
        return_flag = get_unitig(g, ug, v^1, &convex, &ll, &tmp, &max_stop_nodeLen,&max_stop_baseLen, 1, NULL);
        if(return_flag == LOOP) continue;

        as = asg_arc_a(g, convex); ns = asg_arc_n(g, convex);
        for (i = nc = 0; i < ns; i++)
        {
            if(as[i].del) continue;
            if(get_real_length(g, as[i].v^1, NULL) < 2) break;
            nc++;
        }
        if(nc == 0 || i < ns) continue;

        tt.a.n = 0; s^=1; sv = v^1;
        dfs_set(g, sv, s>>1, &stack, &tt, vis, fp);
        p_n = tt.a.n;
        as = asg_arc_a(g, s); ns = asg_arc_n(g, s); found = 0;
        for (i = 0, pm = NULL; i < ns; i++)
        {
            if(as[i].del || as[i].v != sv) continue;
            pm = &(as[i]);
            break;
        }
 
        for (i = 0; i < ns; i++)
        {
            if(as[i].del || as[i].v == sv) continue;
            nc = dfs_set(g, as[i].v, s>>1, &stack, &tt, vis, fa);
            if(nc && check_trans_relation_by_path(sv, as[i].v, &pq_p, path_p, &(o->b0), &pq_a, path_q, &(o->b1), g, 
            vis, fp+fa, nc, 0.45, NULL))
            {
                found = 1;
            }
            a_a = tt.a.a + p_n; a_n = tt.a.n - p_n;
            for (k = 0; k < a_n; k++)
            {
                if(vis[a_a[k]]&fa) vis[a_a[k]] -= fa;
            }
            tt.a.n = p_n;
            if(found) break;
        }
        a_a = tt.a.a; a_n = p_n;
        for (k = 0; k < a_n; k++) vis[a_a[k]] = 0;

        if(found)
        {
            b.b.n = 0;
            get_unitig(g, ug, sv, &convex, &ll, &tmp, &max_stop_nodeLen, &max_stop_baseLen, 1, &b);
            for (k = 0; k < b.b.n; k++)
            {
                g->seq[b.b.a[k]>>1].c = ALTER_LABLE;
            }

            for (k = 0; k < b.b.n; k++)
            {
                asg_seq_drop(g, b.b.a[k]>>1);
            }

            // uint32_t ui;
            // fprintf(stderr, "\n");
            // for (ui = 0; ui < o->b0.b.n; ui++)
            // {
            //     fprintf(stderr, "p-utg%.6ul\n", (o->b0.b.a[ui]>>1)+1);
            // }
            // for (ui = 0; ui < o->b1.b.n; ui++)
            // {
            //     fprintf(stderr, "a-utg%.6ul\n", (o->b1.b.a[ui]>>1)+1);
            // }


            if(o->b0.b.n < b.b.n)  fprintf(stderr, "ERROR-ll\n"); 
            o->b0.b.n = b.b.n; pi = o->k_trans.n;
            collect_trans_ovlp(__func__, &(o->b0), pm->ol, &(o->b1), as[i].ol, ug, o);

            for (k = m = pi; k < o->k_trans.n; k++)
            {
                kt = &(o->k_trans.a[k]);
                if(is_hap_ovlp(ug, read_sg, reverse_sources, ruIndex, &(kt->qn), 1, &(kt->tn), 1) || 
                      is_hap_ovlp(ug, read_sg, reverse_sources, ruIndex, &(kt->tn), 1, &(kt->qn), 1))
                {
                    qn = get_dir_v(kt->qn, &(o->b0));
                    tn = get_dir_v(kt->tn, &(o->b1));

                    if(dfs_reach(o->cug->g, qn, tn, &stack, &tt, vis) || dfs_reach(o->cug->g, tn, qn, &stack, &tt, vis))
                    {
                        // fprintf(stderr, "<<<<<<DEL: q-utg%.6ul (qs: %u, qe: %u), t-utg%.6ul (ts: %u, te: %u)\n", 
                        // kt->qn + 1, kt->qs, kt->qe, kt->tn + 1, kt->ts, kt->te);
                        continue;
                    }

                    // fprintf(stderr, "q-utg%.6ul (qs: %u, qe: %u), t-utg%.6ul (ts: %u, te: %u)\n", 
                    // kt->qn + 1, kt->qs, kt->qe, kt->tn + 1, kt->ts, kt->te);
                    o->k_trans.a[m] = *kt;
                    m++;
                }
                // else
                // {
                //     fprintf(stderr, "******delete: q-utg%.6ul (qs: %u, qe: %u), t-utg%.6ul (ts: %u, te: %u)\n", 
                //     kt->qn + 1, kt->qs, kt->qe, kt->tn + 1, kt->ts, kt->te);
                // }      
            }
            o->k_trans.n = m;
        }
    }

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    kv_destroy(tt.a);
    kv_destroy(stack.a);
    free(vis);
    destory_pdq(&pq_p);
    destory_pdq(&pq_a);
    free(path_p);
    free(path_q);
    return n_reduced;
}

typedef struct {
	uint64_t *idx;
    kvec_t(uint64_t) pos;
} mz_ds_t;

typedef struct {
	uint64_t x, y;
} pt128_t;

typedef struct {
	uint64_t x;
	uint64_t rid:32, span:32;
    uint64_t pos:63, rev:1;
} pt_mz1_t;

typedef struct {
	///cnt1: how many unique minimizers
	///cnt2: how many non-unique minimizers
	uint32_t cnt2, cnt1;
	uint32_t m[2];
	int8_t s;
} pt_uinfo_t;

typedef struct {
	uint32_t n_seq; // number of segments; same as gfa_t::n_seg
	pt_uinfo_t *info; // of size n_seg
    kv_u_trans_t *ma;
} pt_match_t;

#define mz_key(z) ((z).x)
KRADIX_SORT_INIT(mz, pt_mz1_t, mz_key, 8)

#define pt128x_key(z) ((z).x)
KRADIX_SORT_INIT(pt128x, pt128_t, pt128x_key, 8)

#define generic_key(x) (x)
KRADIX_SORT_INIT(tb64, uint64_t, generic_key, 8)

typedef struct { uint32_t n, m; pt128_t *a; } pt128_v;
typedef struct { uint32_t n, m; pt_mz1_t *a; } pt_mz1_v;

static inline int mzcmp(const pt_mz1_t *a, const pt_mz1_t *b)
{
	return (a->x > b->x) - (a->x < b->x);
}

void pt_sketch(const char *str, int len, int w, int k, uint32_t rid, int is_hpc, pt_mz1_v *p)
{
    static const pt_mz1_t dummy = { UINT64_MAX, (1<<28) - 1, 0, 0 };
	uint64_t shift1 = k - 1, mask = (1ULL<<k) - 1, kmer[4] = {0,0,0,0};
	int i, j, l, buf_pos, min_pos, kmer_span = 0;
	pt_mz1_t buf[256], min = dummy;
	tiny_queue_t tq;
    assert(len > 0 && rid < (uint64_t)1<<32 && (w > 0 && w < 256) && (k > 0 && k <= 63));

    memset(buf, 0xff, w * sizeof(pt_mz1_t));
    memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(pt_mz1_t, *p, p->n + len/w);

    for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		pt_mz1_t info = dummy;
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
				if (tq.count > k) kmer_span -= tq_shift(&tq);
			} else kmer_span = l + 1 < k? l + 1 : k; 

			kmer[0] = (kmer[0] << 1 | (c&1))  & mask;                  // forward k-mer
			kmer[1] = (kmer[1] << 1 | (c>>1)) & mask;
			kmer[2] = kmer[2] >> 1 | (uint64_t)(1 - (c&1))  << shift1; // reverse k-mer
			kmer[3] = kmer[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift1;
			if (kmer[1] == kmer[3]) continue; // skip "symmetric k-mers" as we don't know its strand
			z = kmer[1] < kmer[3]? 0 : 1; // strand
			++l;
			if (l >= k && kmer_span < 256) {
				uint64_t y;
				y = yak_hash64_64(kmer[z<<1|0]) + yak_hash64_64(kmer[z<<1|1]);
				info.x = y, info.rid = rid, info.pos = i, info.rev = z, info.span = kmer_span; // initially pt_mz1_t::rid keeps the k-mer count
			}
		} else l = 0, tq.count = tq.front = 0, kmer_span = 0;

		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (mzcmp(&min, &buf[j]) == 0 && buf[j].pos != min.pos) kv_push(pt_mz1_t, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (mzcmp(&min, &buf[j]) == 0 && buf[j].pos != min.pos) kv_push(pt_mz1_t, *p, buf[j]);
		}
		///three cases: 1. 
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k && min.x != UINT64_MAX) kv_push(pt_mz1_t, *p, min);
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1 && min.x != UINT64_MAX) kv_push(pt_mz1_t, *p, min);
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (mzcmp(&min, &buf[j]) >= 0) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (mzcmp(&min, &buf[j]) >= 0) min = buf[j], min_pos = j;

			if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (mzcmp(&min, &buf[j]) == 0 && min.pos != buf[j].pos) kv_push(pt_mz1_t, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (mzcmp(&min, &buf[j]) == 0 && min.pos != buf[j].pos) kv_push(pt_mz1_t, *p, buf[j]);
			}
		}
		if (++buf_pos == w) buf_pos = 0;
	}
	if (min.x != UINT64_MAX)
		kv_push(pt_mz1_t, *p, min);
}

pt_mz1_v *pt_collect_minimizers(ma_ug_t *ug, asg_t *read_g, ma_sub_t *coverage_cut, ma_hit_t_alloc* sources, 
kvec_asg_arc_t_warp* edge, int max_hang, int min_ovlp)
{
    uint32_t i;
    pt_mz1_v *mz = NULL; CALLOC(mz, 1);
    
    for (i = 0; i < ug->u.n; i++)
    {
        if(!ug->u.a[i].len) continue;
        pt_sketch(ug->u.a[i].s, ug->u.a[i].len, asm_opt.mz_win, asm_opt.k_mer_length, i, 0, mz);
    }

    radix_sort_mz(mz->a, mz->a + mz->n);
    return mz;
}

pt128_v *pt_collect_anchors(ma_ug_t *ug, pt_mz1_v *mz, mz_ds_t *mz_idx, uint32_t max_occ)
{
	uint32_t st, j;
    pt128_v *pa = NULL; CALLOC(pa, 1);
    pt128_t *p = NULL;
    mz_idx->pos.n = 0;
	///all minimizers
	for (j = 1, st = 0; j <= mz->n; ++j) {
		if (j == mz->n || mz->a[j].x != mz->a[st].x) {
			uint32_t k, l;
			// if (j - st == 1) ++info[mz->a[st].rid].cnt1; ///of size n_seg
			///max_occ is the frquency threshold of minimizer
			///if (j - st) == 1, means minimizer only occurs in one read, it is not useful
			if (j - st == 1 || j - st > max_occ) goto end_anchor;
			for (k = st; k < j; ++k) {
				// ++info[mz->a[k].rid].cnt2;
                kv_push(uint64_t, mz_idx->pos, (uint64_t)(mz->a[k].rid)<<32|(uint64_t)(mz->a[k].pos));
				for (l = k + 1; l < j; ++l) {
					///k is current minimizer
					uint32_t span, rev = (mz->a[k].rev != mz->a[l].rev);
					int32_t lk = ug->u.a[mz->a[k].rid].len, ll = ug->u.a[mz->a[l].rid].len;
                    kv_pushp(pt128_t, *pa, &p);
					span = mz->a[l].span;
                    p->x = (uint64_t)mz->a[k].rid << 33 | mz->a[l].rid << 1 | rev;
					p->y = (uint64_t)mz->a[k].pos << 32 | (rev? ll - (mz->a[l].pos + 1 - span) - 1 : mz->a[l].pos);
					kv_pushp(pt128_t, *pa, &p);
                    span = mz->a[k].span;
					p->x = (uint64_t)mz->a[l].rid << 33 | mz->a[k].rid << 1 | rev;
					p->y = (uint64_t)mz->a[l].pos << 32 | (rev? lk - (mz->a[k].pos + 1 - span) - 1 : mz->a[k].pos);
				}
			}
            end_anchor:	st = j;
		}
	}

	radix_sort_pt128x(pa->a, pa->a + pa->n);

    radix_sort_tb64(mz_idx->pos.a, mz_idx->pos.a + mz_idx->pos.n);
    CALLOC(mz_idx->idx, ug->u.n);
    for (st = 0, j = 1; j <= mz_idx->pos.n; ++j)
    {
        if (j == mz_idx->pos.n || (mz_idx->pos.a[j]>>32) != (mz_idx->pos.a[st]>>32))
        {
            mz_idx->idx[mz_idx->pos.a[st]>>32] = (uint64_t)st << 32 | (j - st), st = j;
        }
    }
	
    return pa;
}

int32_t pt_lis_64(int32_t n, const uint64_t *a_idx, int32_t *b, int32_t *M)
{
	int32_t i, k, L = 0, *P = b;
	// MALLOC(M, n+1);
	for (i = 0; i < n; ++i) {
		int32_t lo = 1, hi = L, newL;
		while (lo <= hi) {
			int32_t mid = (lo + hi + 1) >> 1;
			if ((uint32_t)a_idx[M[mid]] < (uint32_t)a_idx[i]) lo = mid + 1;
			else hi = mid - 1;
		}
		newL = lo, P[i] = M[newL - 1], M[newL] = i;
		if (newL > L) L = newL;
	}
	k = M[L];
	memcpy(M, P, n * sizeof(int32_t));
	for (i = L - 1; i >= 0; --i) b[i] = k, k = M[k];
	// free(M);
	return L;
}

uint32_t debug_lis_64(const uint64_t *a, int32_t *b, uint64_t n)
{
    uint32_t i;
    if(n <= 1) return 1;
    for (i = 0; i+1 < n; i++)
    {
        if(((a[b[i]]>>32) > (a[b[i+1]]>>32)) || (((uint32_t)a[b[i]]) > ((uint32_t)a[b[i+1]])))
        {
            i = (uint32_t)-1;
            break;
        }
    } 
    if(i == (uint32_t)-1)
    {
        fprintf(stderr, "\nERROR-chain\n");   
        for (i = 0; i < n; i++)
        {
            fprintf(stderr, "x-%lu, y-%lu\n", (a[b[i]]>>32), (uint64_t)((uint32_t)a[b[i]]));
        }
    }
    return 1;
}

void update_mz_ovlp(uint32_t* n_x_beg, uint32_t* n_x_end, int64_t xLen, 
uint32_t* n_y_beg, uint32_t* n_y_end, int64_t yLen, uint32_t rev)
{
    int64_t x_beg = (*n_x_beg), x_end = (*n_x_end);
    int64_t y_beg = (*n_y_beg), y_end = (*n_y_end);
    if(x_beg <= y_beg)
    {
        y_beg = y_beg - x_beg;
        x_beg = 0;
    }
    else
    {
        x_beg = x_beg - y_beg;
        y_beg = 0;
    }

    long long x_right_length = xLen - x_end - 1;
    long long y_right_length = yLen - y_end - 1;


    if(x_right_length <= y_right_length)
    {
        x_end = xLen - 1;
        y_end = y_end + x_right_length;        
    }
    else
    {
        x_end = x_end + y_right_length;
        y_end = yLen - 1;
    }

    if(rev == 0)
    {
        (*n_y_beg) = y_beg;
        (*n_y_end) = y_end + 1;
    }
    else
    {
        (*n_y_beg) = yLen - y_end - 1;
        (*n_y_end) = yLen - y_beg - 1 + 1;
    }

    (*n_x_beg) = x_beg;
    (*n_x_end) = x_end + 1;
}

int64_t get_insert_pos(uint64_t *a, int64_t n, int64_t target)
{
    int64_t left, right, ans, mid;
    left = 0; right = n - 1; ans = n;
    while (left <= right) {
        mid = ((right - left) >> 1) + left;
        if (target <= (uint32_t)a[mid]) {
            ans = mid;
            right = mid - 1;
        } else {
            left = mid + 1;
        }
    }
    return ans;
}

uint32_t get_mz_occ(mz_ds_t *mz_idx, uint64_t uid, uint64_t s, uint64_t e)
{
    uint64_t *a = mz_idx->pos.a + (mz_idx->idx[uid]>>32), n = (uint32_t)(mz_idx->idx[uid]);
    if(n == 0) return 0;
    int64_t sid, eid;
    sid = get_insert_pos(a, n, s);
    eid = get_insert_pos(a, n, e);
    return eid + 1 - sid;
}

kv_u_trans_t *pt_cal_sim(pt128_v *pa, mz_ds_t *mz_idx, ma_ug_t *ug, uint32_t min_cnt, double min_sim)
{
	int64_t st, i, j;
    kvec_t(uint64_t) a; kv_init(a);
    kvec_t(int32_t) b; kv_init(b);
    kvec_t(int32_t) M; kv_init(M);
	kv_u_trans_t *ma = NULL; CALLOC(ma, 1);
    u_trans_t m;
	///x = (uint64_t)mz[l].rid << 33 | mz[k].rid << 1 | rev;
	for (st = 0, i = 1; i <= pa->n; ++i) {
		if (i == pa->n || pa->a[i].x != pa->a[st].x) {///minimizers between a pair of unitigs
            if((pa->a[st].x>>33) == (((uint32_t)pa->a[st].x)>>1)) goto end_chain;
			uint32_t nn[2], nn_min;
            a.n = 0; memset(&m, 0, sizeof(m));
			if (i - st < min_cnt) goto end_chain;
            //(uint64_t)mz[l].pos << 32 | (rev? lk - (mz[k].pos + 1 - span) - 1 : mz[k].pos);
			for (j = st; j < i; ++j) kv_push(uint64_t, a, pa->a[j].y);
			radix_sort_tb64(a.a, a.a + a.n);///sort by query pos + target pos
			// for (l = 0; l < a.n; ++l) a.a[l] = (uint32_t)a.a[l];///only need target pos to do LIS
            kv_resize(int32_t, b, a.n); kv_resize(int32_t, M, a.n+1);
			m.occ = pt_lis_64(a.n, a.a, b.a, M.a);
            /*******************************for debug************************************/
            // debug_lis_64(a.a, b.a, m.occ);
            /*******************************for debug************************************/
			if (m.occ == 0 || m.occ < min_cnt) goto end_chain;//chain occ
            m.qn = pa->a[st].x >> 33;///query id
            m.tn = ((uint32_t)pa->a[st].x) >> 1;///target id
			if (m.qn == m.tn) goto end_chain;
            m.qs = a.a[b.a[0]]>>32; m.qe = a.a[b.a[m.occ-1]]>>32;
            m.ts = (uint32_t)(a.a[b.a[0]]); m.te = (uint32_t)(a.a[b.a[m.occ-1]]);
            m.rev = (pa->a[st].x>>32&1) ^ (pa->a[st].x&1);
            update_mz_ovlp(&(m.qs), &(m.qe), ug->u.a[m.qn].len, &(m.ts), &(m.te), ug->u.a[m.tn].len, m.rev);

            nn[0] = get_mz_occ(mz_idx, m.qn, m.qs, m.qe-1);
            nn[1] = get_mz_occ(mz_idx, m.tn, m.ts, m.te-1);
			nn_min = MIN(nn[0], nn[1]);
            nn_min = MAX(nn_min, m.occ);
            // m.sim = pow(2.0 * m.m / (nn[0] + nn[1]), 1.0 / k);
            if(m.occ >= nn_min*min_sim){
                m.nw = (double)(m.occ) - (double)(nn_min-m.occ)*0.2;
                if(m.nw > 0) kv_push(u_trans_t, *ma, m);
			}
            end_chain:	st = i;
		}
	}
	kv_destroy(b); kv_destroy(a); kv_destroy(M);
	return ma;
}

void clean_mz_ovlp(kv_u_trans_t *ta, ma_ug_t *ug)
{
    u_trans_t *a = NULL;
    asg_arc_t *as = NULL;
    uint32_t k, i, n, v, ns;
    pdq pq; init_pdq(&pq, ug->g->n_seq<<1);
    uint8_t *vis = NULL; CALLOC(vis, ug->g->n_seq);
    kvec_t_u32_warp p; kv_init(p.a);
    for (k = 0; k < ta->idx.n; k++)
    {
        p.a.n = 0;
        a = u_trans_a(*ta, k);
        n = u_trans_n(*ta, k);
        if(n == 0) continue;

        v = k<<1;
        as = asg_arc_a(ug->g, v); ns = asg_arc_n(ug->g, v);
        if(ns > 0)
        {
            for (i = 0; i < ns; i++)
            {
                if(as[i].del) continue;
                kv_push(uint32_t, p.a, as[i].v);
            }
            set_utg_by_dis(v, &pq, ug->g, &p, ug->g->seq[v>>1].len);
        }


        v = (k<<1) + 1;
        as = asg_arc_a(ug->g, v); ns = asg_arc_n(ug->g, v);
        if(ns > 0)
        {
            for (i = 0; i < ns; i++)
            {
                if(as[i].del) continue;
                kv_push(uint32_t, p.a, as[i].v);
            }
            set_utg_by_dis(v, &pq, ug->g, &p, ug->g->seq[v>>1].len);
        }

        for (i = 0; i < p.a.n; i++) vis[p.a.a[i]>>1] = 1; 

        for (i = 0; i < n; i++)
        {
            if(vis[a[i].tn]) a[i].del = 1;
        }

        for (i = 0; i < p.a.n; i++) vis[p.a.a[i]>>1] = 0; 
    }
    destory_pdq(&pq);
    kv_destroy(p.a);
    free(vis);

    for (k = n = 0; k < ta->n; ++k)
    {
        if(ta->a[k].del) continue;
        ta->a[n] = ta->a[k];
        n++;
    }
    ta->n = n;
    kt_u_trans_t_idx(ta, ug->g->n_seq);
    kt_u_trans_t_simple_symm(ta, ug->g->n_seq, 0);
}

int cmp_u_trans_nw(const void * a, const void * b)
{
    if((*(u_trans_t*)a).nw == (*(u_trans_t*)b).nw) return 0;
    return (*(u_trans_t*)a).nw < (*(u_trans_t*)b).nw ? 1 : -1;
}

/**
void flat_mz_ovlp(kv_u_trans_t *ta, ma_ug_t *ug, asg_t *read_g, ma_sub_t* coverage_cut, 
ma_hit_t_alloc* sources, R_to_U* ruIndex, uint32_t min_cnt, uint32_t cov_thres, 
double cov_pass_rate, double nw_pass_rate)
{
    u_trans_t *a = NULL, *p = NULL;
    uint8_t *vis = NULL; CALLOC(vis, ug->g->n_seq);
    uint32_t *cov = NULL; CALLOC(cov, ug->g->n_seq);
    kvec_t(uint32_t) tc; kv_init(tc);
    uint32_t k, i, j, n, v, ns, qs, qe, occ, hapN, pass;
    double w;
    for (i = 0; i < ug->u.n; i++)
    {
        cov[i] = get_utg_cov(ug, i, read_g, coverage_cut, sources, ruIndex, vis);
    }

    for (k = 0; k < ta->idx.n; k++)
    {
        a = u_trans_a(*ta, k);
        n = u_trans_n(*ta, k);
        if(n == 0) continue;
        qsort(a, n, sizeof(u_trans_t), cmp_u_trans_nw);
        kv_resize(uint32_t, tc, ug->u.a[k].len);
        tc.n = ug->u.a[k].len;
        memset(tc.a, 0, tc.n*sizeof(uint32_t));
        for (i = 0, p = NULL; i < n; i++)
        {
            qs = a[i].qs; qe = a[i].qe; w = a[i].nw; occ = a[i].occ; pass = 0;
            for (j = qs; j < qe; j++)
            {
                if(tc.a[j] + cov[k] >= cov_thres) pass++;
            }
            if(pass > cov_pass_rate*(qe - qs))
            {
                if(p && nw_pass_rate*) a[i].del = 1;
            }
            else
            {
                for (j = qs; j < qe; j++) tc.a[j] += cov[k];
            }
        }
    }
    free(vis); free(cov); kv_destroy(tc);
}
**/


void print_u_trans(kv_u_trans_t *ta)
{
    uint32_t i;
    u_trans_t *p = NULL;
    for (i = 0; i < ta->n; i++)
    {
        p = &(ta->a[i]);
        fprintf(stderr, "q-utg%.6ul\tqs(%u)\tqe(%u)\tt-utg%.6ul\tts(%u)\tte(%u)\trev(%u)\tw(%f)\tf(%u)\n", 
        p->qn+1, p->qs, p->qe, p->tn+1, p->ts, p->te, p->rev, p->nw, p->f);
    }
    fprintf(stderr, "[M::%s::] \n", __func__);
}

kv_u_trans_t *pt_pdist(ma_ug_t *ug, asg_t *read_g, ma_sub_t *coverage_cut, ma_hit_t_alloc* sources, 
kvec_asg_arc_t_warp* edge, int max_hang, int min_ovlp, uint32_t min_chain_cnt)
{
    pt_mz1_v *mz = NULL;
    mz_ds_t mz_idx; memset(&mz_idx, 0, sizeof(mz_idx));
    pt128_v *an = NULL;
    kv_u_trans_t *ma = NULL; CALLOC(ma, 1);
    mz = pt_collect_minimizers(ug, read_g, coverage_cut, sources, edge, max_hang, min_ovlp);
    an = pt_collect_anchors(ug, mz, &mz_idx, asm_opt.polyploidy*10);
    kv_destroy(*mz); free(mz);

    ma = pt_cal_sim(an, &mz_idx, ug, min_chain_cnt, asm_opt.purge_simi_thres);
    kv_destroy(*an); free(an); 
    kv_destroy(mz_idx.pos); free(mz_idx.idx);
    
    kt_u_trans_t_idx(ma, ug->g->n_seq);
    kt_u_trans_t_simple_symm(ma, ug->g->n_seq, 1);
    clean_mz_ovlp(ma, ug);
    // print_u_trans(ma);
    // exit(1);
    return ma;
}