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
typedef struct {///[cBeg, cEnd)
	uint32_t ui, len, cBeg, cEnd;
    uint32_t *a, an;
    ma_ug_t *ug;
    utg_trans_t *o;
} utg_trans_hit_idx;

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
        kt->qn = kh->qn>>1; kt->qs = kh->qus; kt->qe = kh->que; kt->qo = kh->qn&1;
        kt->tn = kh->tn>>1; kt->ts = kh->tus; kt->te = kh->tue; kt->to = kh->tn&1;
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