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
KSEQ_INIT(gzFile, gzread)
KDQ_INIT(uint64_t)
#define pe_hit_an1_idx_key(x) ((x).s<<1)
KRADIX_SORT_INIT(pe_hit_idx_hn1, pe_hit, pe_hit_an1_idx_key, member_size(pe_hit, s))
#define pe_hit_an2_idx_key(x) ((x).e<<1)
KRADIX_SORT_INIT(pe_hit_idx_hn2, pe_hit, pe_hit_an2_idx_key, member_size(pe_hit, e))
#define generic_key(x) (x)
KRADIX_SORT_INIT(ho64, uint64_t, generic_key, 8)
#define osg_arc_key(a) ((a).u)
KRADIX_SORT_INIT(osg, osg_arc_t, osg_arc_key, member_size(osg_arc_t, u))

#define OVL(s_0, e_0, s_1, e_1) ((MIN((e_0), (e_1)) > MAX((s_0), (s_1)))? MIN((e_0), (e_1)) - MAX((s_0), (s_1)):0)
#define BREAK_THRES 5000000
#define BREAK_CUTOFF 0.1
#define BREAK_BOUNDARY 0.015


typedef struct {
	uint64_t ruid;
    uint64_t off;
} hit_aux_t;

typedef struct {
	hit_aux_t *a;
    size_t n, m;
    kvec_t(uint64_t) idx;
} u_hits_t;

typedef struct {
	uint64_t e, d;
    double w;
} hw_aux_t;

typedef struct {
	hw_aux_t *a;
    size_t n, m;
} h_w_t;

#define hw_e_key(x) ((x).e)
KRADIX_SORT_INIT(hw_e, hw_aux_t, hw_e_key, member_size(hw_aux_t, e))

#define hw_d_key(x) ((x).d)
KRADIX_SORT_INIT(hw_d, hw_aux_t, hw_d_key, member_size(hw_aux_t, d))

#define hw_ew_key(x) ((uint32_t)((x).e))
KRADIX_SORT_INIT(hw_ew, hw_aux_t, hw_ew_key, member_size(hw_aux_t, e))

#define hw_dw_key(x) ((uint32_t)((x).d))
KRADIX_SORT_INIT(hw_dw, hw_aux_t, hw_dw_key, member_size(hw_aux_t, d))

typedef struct {
    kvec_t(uint64_t) pos;
	uint64_t *a;
    size_t n, m;
} dens_idx_t;

typedef struct {
	uint64_t s, e, dp;
} h_cov_t;

typedef struct {
	h_cov_t *a;
    size_t n, m;
} h_covs;

typedef struct {
    uint32_t *a;
    size_t n, m;
}lay_t;

typedef struct {
    lay_t *a;
    size_t n, m;
}sc_lay_t;

typedef struct {
	uint64_t uid, sid;
    uint64_t iid:63, ori:1;
} sc_id_t;

typedef struct {
    sc_id_t *a;
    size_t n, m;
    sc_lay_t *sl;
    osg_t *sg;
    uint32_t n_thread;
} sc_mul;

#define h_cov_s_key(x) ((x).s)
KRADIX_SORT_INIT(h_cov_s, h_cov_t, h_cov_s_key, member_size(h_cov_t, s))
#define h_cov_e_key(x) ((x).e)
KRADIX_SORT_INIT(h_cov_e, h_cov_t, h_cov_e_key, member_size(h_cov_t, e))
#define h_cov_dp_key(x) ((x).dp)
KRADIX_SORT_INIT(h_cov_dp, h_cov_t, h_cov_dp_key, member_size(h_cov_t, dp))
#define hit_aux_ruid_key(x) ((x).ruid)
KRADIX_SORT_INIT(hit_aux_ruid, hit_aux_t, hit_aux_ruid_key, member_size(hit_aux_t, ruid))

typedef struct {
    kv_u_trans_t *ref;
    trans_chain* idx;
} trans_col_t;

typedef struct {
    uint32_t Spre, Epre, Scur, Ecur, uCur, uPre;///[qSp, qEp) && [qSn, qEn]
} u_hit_t;

#define u_hit_t_key(x) ((x).uPre)
KRADIX_SORT_INIT(u_hit, u_hit_t, u_hit_t_key, member_size(u_hit_t, uPre))

typedef struct {
    ma_ug_t *ug; 
    kvec_pe_hit hits;
    kv_u_trans_t k_trans;
    h_covs *b_points;
} debug_phasing_t;

debug_phasing_t *init_debug_phasing(ma_ug_t *ug, kvec_pe_hit *hits, kv_u_trans_t *k_trans, h_covs *b_points)
{
    debug_phasing_t *p = NULL; CALLOC(p, 1);
    p->ug = copy_untig_graph(ug);
    p->b_points = b_points;

    p->hits.pos_mode = hits->pos_mode; p->hits.uID_bits = hits->pos_mode;
    p->hits.a.n = p->hits.a.m = hits->a.n = hits->a.m;
    MALLOC(p->hits.a.a, p->hits.a.n); 
    memcpy(p->hits.a.a, hits->a.a, p->hits.a.n*sizeof(pe_hit));
    
    p->k_trans.n = p->k_trans.m = k_trans->n;
    MALLOC(p->k_trans.a, p->k_trans.n); 
    memcpy(p->k_trans.a, k_trans->a, p->k_trans.n*sizeof(u_trans_t));

    p->k_trans.idx.n = p->k_trans.idx.m = k_trans->idx.n;
    MALLOC(p->k_trans.idx.a, p->k_trans.idx.n); 
    memcpy(p->k_trans.idx.a, k_trans->idx.a, p->k_trans.idx.n*sizeof(uint64_t));

    return p;
}

void destory_debug_phasing_t(debug_phasing_t **x)
{
    ma_ug_destroy((*x)->ug);
    free((*x)->hits.a.a); free((*x)->hits.idx.a); free((*x)->hits.occ.a);
    free((*x)->k_trans.a); free((*x)->k_trans.idx.a);
    free(*x);
}

uint64_t new_node(uint64_t v, h_covs *join, uint64_t *idx, uint64_t is_ul)
{
    if(v&1) return (((uint32_t)join->a[(idx[v>>1]>>32)+(is_ul?0:(((uint32_t)idx[v>>1])-1))].dp)<<1)+1;
    else return (((uint32_t)join->a[(idx[v>>1]>>32)+(is_ul?(((uint32_t)idx[v>>1])-1):0)].dp)<<1);
}

uint32_t iter_rid(buf_t *b, uint64_t *ui, uint64_t *ri, ma_ug_t *ug)
{
    ma_utg_t *u = NULL;
    while ((*ui) < b->b.n)
    {
        u =  &(ug->u.a[b->b.a[(*ui)]>>1]);
        while ((*ri) < u->n) return u->a[(*ri)++]>>32;
        (*ui)++; (*ri) = 0;
    }

    return (uint32_t)-1;
}

void debug_debug_phasing_t(debug_phasing_t *x, ma_ug_t *cug, kvec_pe_hit *chits, kv_u_trans_t *ck_trans, 
h_covs *join, uint64_t *idx)
{
    uint64_t i, k, m, pv, cv, pui, pri, cui, cri, e;
    uint32_t p_cvx, c_cvx, pr, cr, occ;
    long long p_nodeLen, p_baseLen, t1, t2;
    long long c_nodeLen, c_baseLen;
    asg_arc_t *ap, *ac;
    u_trans_t *pk, *ck, *p;
    buf_t pb, cb;
    memset(&pb, 0, sizeof(buf_t));
    memset(&cb, 0, sizeof(buf_t));

    uint32_t *cnt = NULL; CALLOC(cnt, cug->g->n_seq<<1);
    for (cv = 0; cv < (uint64_t)(cug->g->n_seq<<1); ++cv) {
        ///out-nodes of v
        ac = asg_arc_a(cug->g, cv);
        ///if v just have one out-node, there is no muti-edge
        if (asg_arc_n(cug->g, cv) < 2) continue;
        for (i = 0; i < asg_arc_n(cug->g, cv); ++i) ++cnt[ac[i].v];
        for (i = 0; i < asg_arc_n(cug->g, cv); ++i)
            if (--cnt[ac[i].v] != 0) fprintf(stderr, "ERROR-9\n");
    }
    free(cnt);

    for (e = 0; e < cug->g->n_arc; ++e) {
        uint32_t v = cug->g->arc[e].v^1, u = cug->g->arc[e].ul>>32^1;
        asg_arc_t *av = asg_arc_a(cug->g, v);
        for (i = 0; i < asg_arc_n(cug->g, v); ++i)
            if (av[i].v == u) break;
        if (i == asg_arc_n(cug->g, v)) fprintf(stderr, "ERROR-10\n");
    }


    for (i = 0; i < x->ug->g->n_seq; i++)
    {
        pv = (i<<1);
        cv = new_node(pv, join, idx, 1);
        if(asg_arc_n(x->ug->g, pv) != asg_arc_n(cug->g, cv)) fprintf(stderr, "ERROR-1\n");
        ap = asg_arc_a(x->ug->g, pv); ac = asg_arc_a(cug->g, cv);
        for (k = 0; k < asg_arc_n(x->ug->g, pv); k++)
        {
            for (m = 0; m < asg_arc_n(cug->g, cv); m++)
            {
                if(new_node(ap[k].v, join, idx, 0) == ac[m].v) break;
            }
            if(m >= asg_arc_n(cug->g, cv))
            {
                fprintf(stderr, "\n+ERROR-2\n");
                fprintf(stderr, "+putg%.6lul (%lu), cutg%.6lul (%lu)\n", (pv>>1) + 1, pv&1, (cv>>1) + 1, cv&1);
                fprintf(stderr, "+p-occ: %u, c-occ: %u\n", asg_arc_n(x->ug->g, pv), asg_arc_n(cug->g, cv));
                fprintf(stderr, "+ap[k]-utg%.6ul (%u), new-utg%.6lul (%lu)\n", (ap[k].v>>1)+1, ap[k].v&1, 
                                (new_node(ap[k].v, join, idx, 0)>>1)+1, new_node(ap[k].v, join, idx, 0)&1);
                for (m = 0; m < asg_arc_n(cug->g, cv); m++)
                {
                    fprintf(stderr, "+ac[%lu]-utg%.6ul (%u)\n", m, (ac[m].v>>1) + 1, ac[m].v&1);
                }
            }
        }
        


        pv = (i<<1) + 1;
        cv = new_node(pv, join, idx, 1);
        if(asg_arc_n(x->ug->g, pv) != asg_arc_n(cug->g, cv)) fprintf(stderr, "ERROR-1\n");
        ap = asg_arc_a(x->ug->g, pv); ac = asg_arc_a(cug->g, cv);
        for (k = 0; k < asg_arc_n(x->ug->g, pv); k++)
        {
            for (m = 0; m < asg_arc_n(cug->g, cv); m++)
            {
                if(new_node(ap[k].v, join, idx, 0) == ac[m].v) break;
            }
            if(m >= asg_arc_n(cug->g, cv))
            {
                fprintf(stderr, "\n-ERROR-2\n");
                fprintf(stderr, "-putg%.6lul (%lu), cutg%.6lul (%lu)\n", (pv>>1) + 1, pv&1, (cv>>1) + 1, cv&1);
                fprintf(stderr, "-p-occ: %u, c-occ: %u\n", asg_arc_n(x->ug->g, pv), asg_arc_n(cug->g, cv));
                fprintf(stderr, "-ap[k]-utg%.6ul (%u), new-utg%.6lul (%lu)\n", (ap[k].v>>1)+1, ap[k].v&1, 
                                (new_node(ap[k].v, join, idx, 0)>>1)+1, new_node(ap[k].v, join, idx, 0)&1);
                for (m = 0; m < asg_arc_n(cug->g, cv); m++)
                {
                    fprintf(stderr, "-ac[%lu]-utg%.6ul (%u)\n", m, (ac[m].v>>1) + 1, ac[m].v&1);
                }
            }
        }

        pb.b.n = cb.b.n = 0;
        if(get_unitig(x->ug->g, NULL, i<<1, &p_cvx, &p_nodeLen, &p_baseLen, &t1, &t2, 1, &pb) !=
           get_unitig(cug->g, NULL, new_node((i<<1)+1, join, idx, 1)^1, &c_cvx, &c_nodeLen, &c_baseLen, &t1, &t2, 1, &cb))
        {
            fprintf(stderr, "ERROR-3\n");
        }
        if(new_node(p_cvx, join, idx, 1) != c_cvx) fprintf(stderr, "ERROR-4\n");
        if(p_baseLen != c_baseLen)
        {
            fprintf(stderr, "ERROR-6\n");
            fprintf(stderr, "-putg%.6lul, pb.b.n: %u, p_nodeLen: %lld, p_baseLen: %lld, cb.b.n: %u, c_nodeLen: %lld, c_baseLen: %lld\n", 
            i+1, (uint32_t)pb.b.n, p_nodeLen, p_baseLen, (uint32_t)cb.b.n, c_nodeLen, c_baseLen);
        } 
        
        

        pui = pri = cui = cri = 0;
        while(1)
        {
            pr = iter_rid(&pb, &pui, &pri, x->ug);
            cr = iter_rid(&cb, &cui, &cri, cug);
            if(pr != cr) fprintf(stderr, "ERROR-7\n");
            if(pr == (uint32_t)-1 || cr == (uint32_t)-1) break;
        }
    }

    for (i = 0; i < x->k_trans.n; i++)
    {
        pk = &(x->k_trans.a[i]);
        if(((uint32_t)idx[pk->qn]) <= 1 && ((uint32_t)idx[pk->tn]) <= 1)
        {
            get_u_trans_spec(ck_trans, ((uint32_t)join->a[idx[pk->qn]>>32].dp), 
            ((uint32_t)join->a[idx[pk->tn]>>32].dp), &ck, &occ);
            if(occ != 1 || !ck) fprintf(stderr, "ERROR-8\n");

            if(pk->qs != ck->qs || pk->qe != ck->qe || pk->ts != ck->ts || pk->te != ck->te ||
               pk->f != ck->f || pk->rev != ck->rev || pk->del != ck->del)
            {
                fprintf(stderr, "ERROR-9\n");
            }
        }
        else
        {
            p = pk;
            fprintf(stderr, "\n+q-utg%.6ul\tqs(%u)\tqe(%u)\tt-utg%.6ul\tts(%u)\tte(%u)\trev(%u)\tw(%f)\tf(%u)\n", 
            p->qn+1, p->qs, p->qe, p->tn+1, p->ts, p->te, p->rev, p->nw, p->f);

            uint64_t qi, ti, qid = p->qn, tid = p->tn;
            for (qi = 0; qi < (uint32_t)idx[qid]; qi++)
            {                
                for (ti = 0; ti < (uint32_t)idx[tid]; ti++)
                {
                    get_u_trans_spec(ck_trans, (uint32_t)(join->a[(idx[qid]>>32)+qi].dp),
                    (uint32_t)(join->a[(idx[tid]>>32)+ti].dp), &ck, &occ);
                    // fprintf(stderr, "s-utg%.6ul\td-utg%.6ul\tocc:%u\n", 
                    // (uint32_t)(join->a[(idx[qid]>>32)+qi].dp)+1,
                    // (uint32_t)(join->a[(idx[tid]>>32)+ti].dp)+1,
                    // occ);
                    for (k = 0; k < occ; k++)
                    {
                        p = &(ck[k]);
                        fprintf(stderr, "-q-utg%.6ul\tqs(%u)\tqe(%u)\tt-utg%.6ul\tts(%u)\tte(%u)\trev(%u)\tw(%f)\tf(%u)\n", 
                        p->qn+1, p->qs, p->qe, p->tn+1, p->ts, p->te, p->rev, p->nw, p->f);
                    }
                }
                
            }
            
            
            // get_u_trans_spec(ck_trans, ((uint32_t)join->a[idx[pk->qn]>>32].dp), 
            // ((uint32_t)join->a[idx[pk->tn]>>32].dp), &ck, &occ);
            // for (k = 0; k < occ; k++)
            // {
            //     p = &(ck[k]);
            //     fprintf(stderr, "-q-utg%.6ul\tqs(%u)\tqe(%u)\tt-utg%.6ul\tts(%u)\tte(%u)\trev(%u)\tw(%f)\tf(%u)\n", 
            //     p->qn+1, p->qs, p->qe, p->tn+1, p->ts, p->te, p->rev, p->nw, p->f);
            // }
        }
    }

    free(pb.b.a); free(cb.b.a);
}

void print_N50(ma_ug_t* ug)
{
    kvec_t(uint64_t) b; kv_init(b);
    uint64_t i, s, len;
    for (i = s = 0; i < ug->u.n; ++i) 
    {
        kv_push(uint64_t, b, ug->u.a[i].len);
        s += ug->u.a[i].len;
    }
    len = s;
    
    radix_sort_ho64(b.a, b.a+b.n);
    fprintf(stderr, "[M::%s::] Genome Size: %lu, # Contigs: %u, Largest Contig: %lu\n", 
    __func__, len, (uint32_t)ug->u.n, b.a[b.n-1]);
    i = b.n; s = 0;
    while (i > 0)
    {
        i--;
        s += b.a[i];
        if(s >= (len>>1))
        {
            fprintf(stderr, "[M::%s::] N50: %lu\n", __func__, b.a[i]);
            break;
        }
    }
    
    kv_destroy(b);
}

void print_N50_layout(ma_ug_t* ug, sc_lay_t* sl)
{
    kvec_t(uint64_t) b; kv_init(b);
    lay_t *p = NULL;
    uint64_t i, k, s, ulen, len, occ = 0;
    for (i = s = 0; i < sl->n; ++i) 
    {
        p = &(sl->a[i]);
        for (k = ulen = 0; k < p->n; k+=2)
        {
            ulen += ug->u.a[p->a[k]>>1].len;
        }
        occ += p->n;
        
        kv_push(uint64_t, b, ulen);
        s += ulen;
    }
    len = s;
    
    radix_sort_ho64(b.a, b.a+b.n);
    fprintf(stderr, "[M::%s::] Scaffold Size: %lu, # Scaffolds: %u (occ-%lu), Largest Scaffold: %lu\n", 
    __func__, len, (uint32_t)sl->n, occ, b.a[b.n-1]);
    i = b.n; s = 0;
    while (i > 0)
    {
        i--;
        s += b.a[i];
        if(s >= (len>>1))
        {
            fprintf(stderr, "[M::%s::] N50: %lu\n", __func__, b.a[i]);
            break;
        }
    }
    
    kv_destroy(b);
}

trans_col_t *init_trans_col(ma_ug_t *ug, uint64_t r_num, kv_u_trans_t *ref)
{
    trans_col_t *p = NULL; 
    CALLOC(p, 1);
    p->ref = ref;
    p->idx = init_trans_chain(ug, r_num);
    return p;
}

void destory_trans_col(trans_col_t **p)
{
    destory_trans_chain(&((*p)->idx));
    free(*p);
}

void resolve_hit(uint64_t x, uint32_t rLen, uint64_t uID_bits, uint64_t pos_mode, uint64_t *uid, uint64_t *beg, uint64_t *end)
{
    if(uid) (*uid) = ((x<<1)>>(64 - uID_bits)); 
    uint32_t rev = (x>>63);
    long long ref_p = x & pos_mode;
    long long p_beg, p_end;

    if(rev)
    {
        p_end = ref_p;
        p_beg = p_end + 1 - rLen;
    }
    else
    {
        p_beg = ref_p; 
        p_end = p_beg + rLen - 1;
    }
    if(p_beg < 0) p_beg = 0;
    if(p_end < 0) p_end = 0;
    if(beg) (*beg) = p_beg; 
    if(end) (*end) = p_end + 1;
}

void idx_hits(kvec_pe_hit* hits, uint64_t n)
{
    uint64_t k, l;
    kv_resize(uint64_t, hits->idx, n);
    hits->idx.n = n;
    memset(hits->idx.a, 0, hits->idx.n*sizeof(uint64_t));

    radix_sort_pe_hit_idx_hn1(hits->a.a, hits->a.a + hits->a.n);
    for (k = 1, l = 0; k <= hits->a.n; ++k) 
    {   
        if (k == hits->a.n || (get_hit_suid(*hits, k) != get_hit_suid(*hits, l))) 
        {
            if (k - l > 1) radix_sort_pe_hit_idx_hn2(hits->a.a + l, hits->a.a + k);
            
            hits->idx.a[get_hit_suid(*hits, l)] = (uint64_t)l << 32 | (k - l);            
            l = k;
        }
    }
}

kvec_pe_hit *get_r_hits_for_trio(kvec_pe_hit *u_hits, asg_t* r_g, ma_ug_t* ug, bubble_type* bub, uint64_t uID_bits, uint64_t pos_mode)
{
    kvec_pe_hit *r_hits = NULL;
    CALLOC(r_hits, 1);
    uint64_t k, l, i, r_i, offset, rid, rev, rBeg, rEnd, ubits, p_mode, upos, rpos, update, ubeg, uend, suid, euid;
    ma_utg_t *u = NULL;
    memset(r_hits, 0, sizeof(*r_hits));
    r_hits->uID_bits = uID_bits;
    r_hits->pos_mode = pos_mode;
    //reset for reads
    for (ubits=1; (uint64_t)(1<<ubits)<(uint64_t)r_g->n_seq; ubits++);
    p_mode = ((uint64_t)-1) >> (ubits + 1);

    u_hits->uID_bits = uID_bits; u_hits->pos_mode = pos_mode;
    for (i = r_i = 0; i < u_hits->a.n; i++)
    {
        suid = get_hit_suid(*u_hits, i);
        euid = get_hit_euid(*u_hits, i);
        if(IF_HOM(suid, *bub)) continue;
        if(IF_HOM(euid, *bub)) continue;
        if(suid == euid) continue;
        kv_push(pe_hit, r_hits->a, u_hits->a.a[i]);

        resolve_hit(r_hits->a.a[r_i].s, r_hits->a.a[r_i].len>>32, r_hits->uID_bits, 
                                                                r_hits->pos_mode, NULL, &ubeg, &uend);
        upos = (ubeg+uend-1)>>1;
        r_hits->a.a[r_i].s -= get_hit_spos(*r_hits, r_i);
        r_hits->a.a[r_i].s += upos;


        resolve_hit(r_hits->a.a[r_i].e, (uint32_t)r_hits->a.a[r_i].len, r_hits->uID_bits, 
                                                                r_hits->pos_mode, NULL, &ubeg, &uend);
        upos = (ubeg+uend-1)>>1;
        r_hits->a.a[r_i].e -= get_hit_epos(*r_hits, r_i);///pos at unitig
        r_hits->a.a[r_i].e += upos;
        
        r_hits->a.a[r_i].id = (suid<<32)|euid;
        r_i++;
    }

    radix_sort_pe_hit_idx_hn1(r_hits->a.a, r_hits->a.a + r_hits->a.n);
    for (k = 1, l = 0; k <= r_hits->a.n; ++k) 
    {   
        if (k == r_hits->a.n || get_hit_suid(*r_hits, k) != get_hit_suid(*r_hits, l))//same suid
        {
            ///already sort by spos
            u = &(ug->u.a[get_hit_suid(*r_hits, l)]);
            update = 0;
            for (i = offset = 0, r_i = l; i < u->n; i++)
            {
                rid = u->a[i]>>33; 
                rBeg = offset;
                rEnd = rBeg + r_g->seq[rid].len - 1;
                for (; r_i < k; r_i++)
                {
                    upos = get_hit_spos(*r_hits, r_i);///pos at unitig
                    
                    if(upos > rEnd) break;
                    if(upos >= rBeg && upos <= rEnd)
                    {
                        rpos = (((u->a[i]>>32)&1)? rEnd - upos : upos - rBeg);///pos at read
                        rev = ((u->a[i]>>32)&1) ^ (r_hits->a.a[r_i].s>>63);
                        r_hits->a.a[r_i].s = (rev<<63) | ((rid << (64-ubits))>>1) | (rpos & p_mode);
                        
                        update++;
                    }
                }
                offset += (uint32_t)u->a[i];
            }

            if(r_i != k || update != k - l) fprintf(stderr, "ERROR-r_i\n");
            l = k;
        }
    }

    radix_sort_pe_hit_idx_hn2(r_hits->a.a, r_hits->a.a + r_hits->a.n);
    for (k = 1, l = 0; k <= r_hits->a.n; ++k) 
    {
        if (k == r_hits->a.n || get_hit_euid(*r_hits, k) != get_hit_euid(*r_hits, l))//same euid
        {
            ///already sort by epos
            u = &(ug->u.a[get_hit_euid(*r_hits, l)]);
            update = 0;
            for (i = offset = 0, r_i = l; i < u->n; i++)
            {
                rid = u->a[i]>>33; 
                rBeg = offset;
                rEnd = rBeg + r_g->seq[rid].len - 1;
                for (; r_i < k; r_i++)
                {
                    upos = get_hit_epos(*r_hits, r_i);///pos at unitig

                    if(upos > rEnd) break;
                    if(upos >= rBeg && upos <= rEnd)
                    {
                        rpos = (((u->a[i]>>32)&1)? rEnd - upos : upos - rBeg);///pos at read
                        rev = ((u->a[i]>>32)&1) ^ (r_hits->a.a[r_i].e>>63);
                        r_hits->a.a[r_i].e = (rev<<63) | ((rid << (64-ubits))>>1) | (rpos & p_mode);
                        
                        update++;
                    }
                }
                offset += (uint32_t)u->a[i];
            }

            if(r_i != k || update != k - l) fprintf(stderr, "ERROR-r_i\n");
            l = k;
        }
    }


    r_hits->uID_bits = ubits;
    r_hits->pos_mode = p_mode;
    return r_hits;
}

void get_r_hits(kvec_pe_hit *u_hits, kvec_pe_hit *r_hits, asg_t* r_g, ma_ug_t* ug, bubble_type* bub, uint64_t uID_bits, uint64_t pos_mode)
{
    uint64_t k, l, i, r_i, offset, rid, rev, rBeg, rEnd, ubits, p_mode, upos, rpos, update;
    ma_utg_t *u = NULL;
    memset(r_hits, 0, sizeof(*r_hits));
    r_hits->uID_bits = uID_bits;
    r_hits->pos_mode = pos_mode;
    //reset for reads
    for (ubits=1; (uint64_t)(1<<ubits)<(uint64_t)r_g->n_seq; ubits++);
    p_mode = ((uint64_t)-1) >> (ubits + 1);

    kv_malloc(r_hits->a, u_hits->a.n); r_hits->a.n = r_hits->a.m = u_hits->a.n;
    memcpy(r_hits->a.a, u_hits->a.a, r_hits->a.n*sizeof(pe_hit));
    radix_sort_pe_hit_idx_hn1(r_hits->a.a, r_hits->a.a + r_hits->a.n);
    for (k = 1, l = 0; k <= r_hits->a.n; ++k) 
    {   
        if (k == r_hits->a.n || get_hit_suid(*r_hits, k) != get_hit_suid(*r_hits, l))//same suid
        {
            ///already sort by spos
            u = &(ug->u.a[get_hit_suid(*r_hits, l)]);
            update = 0;
            for (i = offset = 0, r_i = l; i < u->n; i++)
            {
                rid = u->a[i]>>33; 
                rBeg = offset;
                rEnd = rBeg + r_g->seq[rid].len - 1;
                for (; r_i < k; r_i++)
                {
                    upos = get_hit_spos(*r_hits, r_i);///pos at unitig
                    
                    if(upos > rEnd) break;
                    if(upos >= rBeg && upos <= rEnd)
                    {
                        if(bub)
                        {
                            r_hits->a.a[r_i].id = (uint32_t)r_hits->a.a[r_i].id;
                            if(!IF_HOM(get_hit_suid(*r_hits, r_i), *bub))
                            {
                                r_hits->a.a[r_i].id += ((uint64_t)(1)<<32);
                            }
                        }

                        rpos = (((u->a[i]>>32)&1)? rEnd - upos : upos - rBeg);///pos at read
                        rev = ((u->a[i]>>32)&1) ^ (r_hits->a.a[r_i].s>>63);
                        r_hits->a.a[r_i].s = (rev<<63) | ((rid << (64-ubits))>>1) | (rpos & p_mode);
                        
                        update++;
                    }
                }
                offset += (uint32_t)u->a[i];
            }

            if(r_i != k || update != k - l) fprintf(stderr, "ERROR-r_i\n");
            l = k;
        }
    }

    radix_sort_pe_hit_idx_hn2(r_hits->a.a, r_hits->a.a + r_hits->a.n);
    for (k = 1, l = 0; k <= r_hits->a.n; ++k) 
    {
        if (k == r_hits->a.n || get_hit_euid(*r_hits, k) != get_hit_euid(*r_hits, l))//same euid
        {
            ///already sort by epos
            u = &(ug->u.a[get_hit_euid(*r_hits, l)]);
            update = 0;
            for (i = offset = 0, r_i = l; i < u->n; i++)
            {
                rid = u->a[i]>>33; 
                rBeg = offset;
                rEnd = rBeg + r_g->seq[rid].len - 1;
                for (; r_i < k; r_i++)
                {
                    upos = get_hit_epos(*r_hits, r_i);///pos at unitig

                    if(upos > rEnd) break;
                    if(upos >= rBeg && upos <= rEnd)
                    {
                        if(bub)
                        {
                            r_hits->a.a[r_i].id >>= 32;
                            r_hits->a.a[r_i].id <<= 32;
                            if(!IF_HOM(get_hit_euid(*r_hits, r_i), *bub))
                            {
                                r_hits->a.a[r_i].id += 1;
                            }
                        }

                        rpos = (((u->a[i]>>32)&1)? rEnd - upos : upos - rBeg);///pos at read
                        rev = ((u->a[i]>>32)&1) ^ (r_hits->a.a[r_i].e>>63);
                        r_hits->a.a[r_i].e = (rev<<63) | ((rid << (64-ubits))>>1) | (rpos & p_mode);
                        
                        update++;
                    }
                }
                offset += (uint32_t)u->a[i];
            }

            if(r_i != k || update != k - l) fprintf(stderr, "ERROR-r_i\n");
            l = k;
        }
    }

    r_hits->uID_bits = ubits;
    r_hits->pos_mode = p_mode;
    idx_hits(r_hits, r_g->n_seq);
}

uint64_t get_corresp_usite(uint64_t rid, uint64_t rpos, uint64_t rev, uint64_t rlen, u_hits_t *x, uint64_t ubits, uint64_t p_mode, kvec_t_u64_warp *buf)
{
    hit_aux_t *a = NULL;
    uint64_t a_n, i, new_uid, new_pos, new_rev;
    
    a = x->a + (x->idx.a[rid]>>32);
    a_n = (uint32_t)x->idx.a[rid];
    for (i = 0; i < a_n; i++)
    {
        new_uid = (uint32_t)a[i].ruid;
        new_pos = (((a[i].ruid>>32)&1)? a[i].off + rlen - 1 - rpos : a[i].off + rpos);
        new_rev = ((a[i].ruid>>32)&1)^rev;
        kv_push(uint64_t, buf->a, (new_rev<<63) | ((new_uid << (64-ubits))>>1) | (new_pos & p_mode));
    }

    return a_n;
}

void update_u_hits(kvec_pe_hit *u_hits, kvec_pe_hit *r_hits, ma_ug_t* ug, asg_t* r_g)
{
    u_hits_t x; memset(&x, 0, sizeof(x));
    hit_aux_t *p = NULL;
    ma_utg_t *u = NULL;
    pe_hit *t = NULL;
    uint64_t v, i, l, k, offset, occ_1, occ_2, *a_1, *a_2, i_1, i_2;
    
    for (v = 0; v < ug->u.n; v++)
    {
        u = &(ug->u.a[v]);
        for (i = offset = 0; i < u->n; i++)
        {
            if(u->a[i] != (uint64_t)-1)
            {
                kv_pushp(hit_aux_t, x, &p);
                p->ruid = u->a[i]>>32; 
                p->ruid <<= 32;
                p->ruid |= v;
                p->off = offset;
                offset += (uint32_t)u->a[i];
            }
            else
            {
                offset += GAP_LEN;
            }
        }
    }

    radix_sort_hit_aux_ruid(x.a, x.a + x.n);
    x.idx.n = x.idx.m = (x.n?(x.a[x.n-1].ruid>>33)+1:0);///how many unitigs
    CALLOC(x.idx.a, x.idx.n);
    for (k = 1, l = 0; k <= x.n; ++k) 
    {   
        if (k == x.n || (x.a[k].ruid>>33) != (x.a[l].ruid>>33))//same rid
        {
            x.idx.a[x.a[l].ruid>>33] = (uint64_t)l << 32 | (k - l);
            l = k;
        }
    }

    u_hits->a.n = u_hits->idx.n = u_hits->occ.n = 0;
    for (u_hits->uID_bits=1; (uint64_t)(1<<u_hits->uID_bits)<(uint64_t)ug->u.n; u_hits->uID_bits++);
    u_hits->pos_mode = ((uint64_t)-1) >> (u_hits->uID_bits + 1);
    kvec_t_u64_warp buf; kv_init(buf.a);

    for (i = 0; i < r_hits->a.n; i++)
    {
        buf.a.n = 0;
        occ_1 = get_corresp_usite(get_hit_suid(*r_hits, i), get_hit_spos(*r_hits, i), 
        r_hits->a.a[i].s>>63, r_g->seq[get_hit_suid(*r_hits, i)].len, &x, u_hits->uID_bits, 
        u_hits->pos_mode, &buf);
        occ_2 = get_corresp_usite(get_hit_euid(*r_hits, i), get_hit_epos(*r_hits, i),
        r_hits->a.a[i].e>>63, r_g->seq[get_hit_euid(*r_hits, i)].len, &x, u_hits->uID_bits, 
        u_hits->pos_mode, &buf);
        if(occ_1 == 0 || occ_2 == 0) continue;
        a_1 = buf.a.a; a_2 = buf.a.a + occ_1;
        for (i_1 = 0; i_1 < occ_1; i_1++)
        {
            for (i_2 = 0; i_2 < occ_2; i_2++)
            {
                kv_pushp(pe_hit, u_hits->a, &t);
                t->id = ((occ_1 == 1) && (occ_2 == 1));
                t->len = r_hits->a.a[i].len;
                t->s = a_1[i_1];
                t->e = a_2[i_2];
            }
        }
    }
    
    free(x.a); free(x.idx.a); kv_destroy(buf.a);
    idx_hits(u_hits, ug->u.n);
}

ma_ug_t* get_trio_unitig_graph(asg_t *sg, uint8_t flag, ug_opt_t *opt)
{
    kvec_asg_arc_t_warp new_rtg_edges;
    kv_init(new_rtg_edges.a);

    ma_ug_t *ug = NULL;
    ug = ma_ug_gen(sg);

    adjust_utg_by_trio(&ug, sg, flag, TRIO_THRES, opt->sources, opt->reverse_sources, 
    opt->coverage_cut, opt->tipsLen, opt->tip_drop_ratio, opt->stops_threshold, 
    opt->ruIndex, opt->chimeric_rate, opt->drop_ratio, opt->max_hang, opt->min_ovlp,
    &new_rtg_edges, opt->b_mask_t);

    kv_destroy(new_rtg_edges.a);
    return ug;
}

static inline void asg_arc_unique_del(asg_t *g, uint32_t v, uint32_t w, int del)
{
	uint32_t i, nv = asg_arc_n(g, v);
	asg_arc_t *av = asg_arc_a(g, v);
	for (i = 0; i < nv; ++i)
    {
        if (av[i].v == w)
        {
            av[i].del = !!del;
            break;
        }
    }
}
void horder_clean_sg_by_utg(asg_t *sg, ma_ug_t *ug)
{
    uint32_t i, v, n_vx, w, k, m, nv, vx, wx;
    asg_arc_t *av = NULL;
    ma_utg_t *u = NULL; 
    
    n_vx = sg->n_seq<<1;
    for (v = 0; v < n_vx; v++)
    {
        nv = asg_arc_n(sg, v);
        av = asg_arc_a(sg, v);
        for (m = 0; m < nv; m++) av[m].del = (!!1);
        sg->seq[v>>1].del = (!!1);
    }
    
    for (i = 0; i < ug->g->n_seq; ++i) 
    {
        if(ug->g->seq[i].del) continue;
        u = &(ug->u.a[i]);
        for (k = 0; (k + 1) < u->n; k++)
        {
            v = u->a[k]>>32; w = u->a[k+1]>>32;

            asg_arc_unique_del(sg, v, w, 0);
            asg_arc_unique_del(sg, w^1, v^1, 0);
        }
        for (k = 0; k < u->n; k++)
        {
            sg->seq[u->a[k]>>33].del = (!!0);
        }

        v = i<<1;
        nv = asg_arc_n(ug->g, v); av = asg_arc_a(ug->g, v);
        for (k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            w = av[k].v;

            vx = (v&1?((ug->u.a[v>>1].a[0]>>32)^1):(ug->u.a[v>>1].a[ug->u.a[v>>1].n-1]>>32));
            wx = (w&1?((ug->u.a[w>>1].a[ug->u.a[w>>1].n-1]>>32)^1):(ug->u.a[w>>1].a[0]>>32));
            asg_arc_unique_del(sg, vx, wx, 0); asg_arc_unique_del(sg, wx^1, vx^1, 0);
        }

        v = (i<<1)+1;
        nv = asg_arc_n(ug->g, v); av = asg_arc_a(ug->g, v);
        for (k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            w = av[k].v;

            vx = (v&1?((ug->u.a[v>>1].a[0]>>32)^1):(ug->u.a[v>>1].a[ug->u.a[v>>1].n-1]>>32));
            wx = (w&1?((ug->u.a[w>>1].a[ug->u.a[w>>1].n-1]>>32)^1):(ug->u.a[w>>1].a[0]>>32));
            asg_arc_unique_del(sg, vx, wx, 0); asg_arc_unique_del(sg, wx^1, vx^1, 0);
        }         
    }
    asg_cleanup(sg);

    /*******************************for debug************************************/
    // ma_ug_t *dbg = ma_ug_gen(sg);
    // print_N50(dbg);
    // print_N50(ug);

    // uint8_t *end = NULL; CALLOC(end, sg->n_seq<<1);
    // for (i = 0; i < dbg->g->n_seq; ++i) 
    // {
    //     u = &(dbg->u.a[i]);
    //     if(u->n == 0) continue;
    //     end[(u->a[0]>>32)^1] = 1;
    //     end[u->a[u->n-1]>>32] = 2;
    // }

    // for (i = 0; i < ug->g->n_seq; ++i) 
    // {
    //     u = &(ug->u.a[i]);
    //     if(u->n == 0) continue;
    //     for (k = 1; (k + 1) < u->n; k++)
    //     {
    //         if(end[(u->a[k]>>32)])
    //         {
    //             fprintf(stderr, "(1) node-%lu, v-%lu, w-%lu, sg(v).n: %u\n", 
    //             u->a[k]>>33, (u->a[k]>>32), (u->a[k+1]>>32), asg_arc_n(sg, (u->a[k]>>32)));
    //         }

    //         if(end[(u->a[k]>>32)^1])
    //         {
    //             fprintf(stderr, "(2) node-%lu, v-%lu, w-%lu, sg(v).n: %u\n", 
    //             u->a[k]>>33, (u->a[k]>>32)^1, (u->a[k-1]>>32)^1, asg_arc_n(sg, (u->a[k]>>32)^1));
    //         }
    //     }
    // }
    // free(end);
    // ma_ug_destroy(dbg);
    /*******************************for debug************************************/
}

uint64_t get_hic_cov_interval(uint64_t *b, uint64_t b_n, int64_t min_dp, int64_t *boundS, int64_t *boundE,
h_covs *res)
{
    if(res) res->n = 0;
    if(min_dp == 0 || b_n == 0) return (uint64_t)-1;
    uint64_t i, len = 0;
    int64_t dp, old_dp, start = 0, bs = b[0]>>1, be = b[b_n-1]>>1, olen;
    h_cov_t *p = NULL;
    if(boundS) bs = (*boundS);
    if(boundE) be = (*boundE);
    for (i = 0, dp = 0, start = 0; i < b_n; ++i) 
    {
        old_dp = dp;
        ///if a[j] is qe
        if (b[i]&1) --dp;
        else ++dp;

        if (old_dp < min_dp && dp >= min_dp) ///old_dp < dp, b.a[j] is qs
        { 
            ///case 2, a[j] is qs
            start = b[i]>>1;
        } 
        else if (old_dp >= min_dp && dp < min_dp) ///old_dp > min_dp, b.a[j] is qe
        {
            olen = OVL(start, (int64_t)(b[i]>>1), bs, be);
            if(olen == 0) continue;
            if(res)
            {
                kv_pushp(h_cov_t, *res, &p);
                p->s = MAX(start, bs); 
                p->e = MIN((int64_t)(b[i]>>1), be);
                p->dp = old_dp;
            }
            len += olen;
        }
    }
    return len;
}

void get_hic_breakpoint(uint64_t *b, uint64_t b_n, int64_t cutoff, h_covs *res, 
int64_t cov_s_pos, int64_t cov_e_pos, uint64_t *s, uint64_t *e)
{
    uint64_t i;
    (*s) = (*e) = (uint64_t)-1;
    res->n = 0;
    get_hic_cov_interval(b, b_n, cutoff, &cov_s_pos, &cov_e_pos, res);
    if(res->n == 0) return;
    int64_t max = -1, max_cur = 0;
    int64_t max_s_idx, max_e_idx, cur_s_idx;
    for (i = 0; i < res->n; i++)//all intervals have cov >= cutoff
    {
        if(i > 0 && (res->a[i].s - res->a[i-1].e) > 0)//cov < cutoff
        {
            max_cur += (res->a[i].s - res->a[i-1].e);//at least positive
            if(max < max_cur)
            {
                max_s_idx = (max < 0? res->a[i-1].e:cur_s_idx);
                max_e_idx = res->a[i].s;
                max = max_cur;
                cur_s_idx = max_s_idx;
            }
        }

        //cov >= cutoff
        max_cur -= (res->a[i].e - res->a[i].s);
        if(max_cur < 0)
        {
            max_cur = 0;
            cur_s_idx = res->a[i].e;
        }
    }

    if(max > 0)
    {
        (*s) = max_s_idx;
        (*e) = max_e_idx;
    }
}

void get_consensus_break(h_covs *res, h_covs *tmp)
{
    uint64_t i, k, n, m, max_cut = 0;
    h_cov_t *p = NULL;
    tmp->n = 0;
    if(res->n == 0) return;
    n = res->n;
    for (i = 0; i < n; i++)
    {
        res->a[i].s <<= 1;
        if(max_cut < res->a[i].dp) max_cut = res->a[i].dp;
        kv_pushp(h_cov_t, *res, &p);
        *p = res->a[i];
        p->s = (res->a[i].e<<1)|1;
    }

    radix_sort_h_cov_s(res->a, res->a+res->n);
    int64_t dp, old_dp, start = 0, max_dp;
    p = NULL; max_dp = -1; tmp->n = 0;
    for (i = 0, dp = 0, start = 0; i < res->n; ++i) 
    {
        old_dp = dp;
        ///if a[j] is qe
        if (res->a[i].s&1) --dp;
        else ++dp;

        if (old_dp < dp) ///old_dp < dp, b.a[j] is qs
        { 
            ///case 2, a[j] is qs
            start = res->a[i].s>>1;
        } 
        else if (old_dp > dp) ///old_dp > min_dp, b.a[j] is qe
        {
            if(max_dp < old_dp)
            {
                max_dp = old_dp;
                tmp->n = 0;
                kv_pushp(h_cov_t, *tmp, &p);
                p->s = start; p->e = res->a[i].s>>1; p->dp = max_cut;
            }
            else if(max_dp == old_dp)
            {
                kv_pushp(h_cov_t, *tmp, &p);
                p->s = start; p->e = res->a[i].s>>1; p->dp = max_cut;
            }
        }
    }
    if(tmp->n == 0) fprintf(stderr, "ERROR-break-0\n");
    // if(tmp->n == 1) return;

    for (i = m = 0; i < res->n; ++i)
    {
        if(res->a[i].s&1) continue;
        res->a[m] = res->a[i];
        res->a[m].s >>= 1;
        m++;
    }
    res->n = m;
    if(res->n != n) fprintf(stderr, "ERROR-break-1\n");
    for (i = 0; i < res->n; ++i)
    {
        for (k = 0; k < tmp->n; k++)
        {
            if(OVL(res->a[i].s, res->a[i].e, tmp->a[k].s, tmp->a[k].e) == 0) continue;
            tmp->a[k].dp = MIN(tmp->a[k].dp, res->a[i].dp);
            if(max_cut > tmp->a[k].dp) max_cut = tmp->a[k].dp;
        }
    }

    for (k = m = 0; k < tmp->n; k++)
    {
        if(max_cut != tmp->a[k].dp) continue;
        tmp->a[m] = tmp->a[k];
        m++;
    }
    tmp->n = m;
    if(tmp->n == 0) fprintf(stderr, "ERROR-break-2\n");
}

int64_t update_r_break(uint64_t rs, uint64_t re, h_covs *hits)
{
    uint64_t i, hs, he;
    int64_t dp, old_dp, max_dp = 0;
    for (i = 0, dp = 0; i < hits->n; ++i) 
    {
        if(hits->a[i].s&1)
        {
            hs = hits->a[i].e;
            he = hits->a[i].s>>1;
        }
        else
        {
            hs = hits->a[i].s>>1;
            he = hits->a[i].e;
        }
        if(hs <= rs && he >= re)
        {
            old_dp = dp;
            ///if a[j] is qe
            if (hits->a[i].s&1) --dp;
            else ++dp;
            ///hits->a[i].s is qe
            if (old_dp > dp && max_dp < old_dp)
            {
                 max_dp = old_dp;
            }
        }
    }

    return max_dp;
}

void get_read_breaks(ma_utg_t *u, asg_t* r_g, h_covs *cov, h_covs *hit_tmp, 
kvec_pe_hit *hits, uint64_t sidx, uint64_t eidx, uint64_t ulen, uint64_t *idx, uint64_t *rdp)
{
    (*idx) = (*rdp) = (uint64_t)-1;
    uint64_t i, k, offset, beg, end, n = cov->n, min_ovlp, o, mi, dp, min_dp;
    uint64_t p0s, p0e, p1s, p1e, span_s, span_e;
    h_cov_t *p = NULL, *a = NULL;
    for (i = offset = 0; i < u->n; i++)
    {
        end = offset + ((u->a[i] != (uint64_t)-1?r_g->seq[u->a[i]>>33].len:GAP_LEN));
        offset += (u->a[i] != (uint64_t)-1? (uint32_t)u->a[i]:GAP_LEN);
        beg = offset;

        if(u->a[i] == (uint64_t)-1) beg -= GAP_LEN;
        if(end <= beg && i + 1 < u->n && u->a[i+1] == (uint64_t)-1)
        {
            end = beg + GAP_LEN;
        }
         
        for (k = 0; k < n; k++)
        {
            if(OVL(beg, end, cov->a[k].s, cov->a[k].e) == 0) continue;
            kv_pushp(h_cov_t, *cov, &p);
            p->s = beg; p->e = end; p->dp = i; 
            break;
        }
    }
    a = cov->a + n;
    n = cov->n - n;
    if(n == 0) fprintf(stderr, "ERROR-r-break\n");
    hit_tmp->n = 0; 
    for (i = sidx; i < eidx; i++)///keep all hic hits that contain interval we want
    {
        if(get_hit_suid(*hits, i) != get_hit_euid(*hits, i)) continue;
        p0s = get_hit_spos(*hits, i);
        p0e = get_hit_spos_e(*hits, i);
        p1s = get_hit_epos(*hits, i);
        p1e = get_hit_epos_e(*hits, i);

        span_s = MIN(MIN(p0s, p0e), MIN(p1s, p1e));
        span_s = MIN(span_s, ulen-1);
        span_e = MAX(MAX(p0s, p0e), MAX(p1s, p1e));
        span_e = MIN(span_e, ulen-1) + 1;

        //if(span_e - span_s <= ulen*BREAK_CUTOFF)//need it or not?
        {
            for (k = 0; k < n; k++)
            {
                if(span_s <= a[k].s && span_e >= a[k].e) break;
            }
            if(k >= n) continue;
            kv_pushp(h_cov_t, *hit_tmp, &p);
            p->s = (span_s<<1); p->e = span_e;
            kv_pushp(h_cov_t, *hit_tmp, &p);
            p->s = ((span_e<<1)|1); p->e = span_s;
        }
    }

    radix_sort_h_cov_s(hit_tmp->a, hit_tmp->a+hit_tmp->n);
    for (k = 0, min_dp = (uint64_t)-1; k < n; k++)
    {
        //a[k].s, a[k].e
        dp = update_r_break(a[k].s, a[k].e, hit_tmp);
        a[k].dp = (uint32_t)a[k].dp;
        a[k].dp += (dp << 32);
        if(dp < min_dp) min_dp = dp;
    }

    min_ovlp = mi = (uint64_t)-1;
    for (k = 0; k < n; k++)
    {
        if((a[k].dp>>32) != min_dp) continue;
        i = (uint32_t)a[k].dp;
        o = 0;
        if(u->a[i] != (uint64_t)-1)
        {
            o = r_g->seq[u->a[i]>>33].len - (uint32_t)u->a[i];
        }
        if(o < min_ovlp)
        {
            min_ovlp = o;
            mi = i;
        }
    }
    if(mi != (uint64_t)-1) (*idx) = mi, (*rdp) = min_dp;
}

void debug_sub_cov(kvec_pe_hit *hits, uint64_t sidx, uint64_t eidx, uint64_t ulen, ma_utg_t *u, asg_t* r_g,
uint64_t rid, uint64_t i_cnt)
{
    uint64_t i, offset, beg, end, rs, re;
    uint64_t p0s, p0e, p1s, p1e, span_s, span_e, cnt = 0;
    rs = re = (uint64_t)-1;
    for (i = offset = 0; i < u->n; i++)
    {
        end = offset + ((u->a[i] != (uint64_t)-1?r_g->seq[u->a[i]>>33].len:GAP_LEN));
        offset += (u->a[i] != (uint64_t)-1? (uint32_t)u->a[i]:GAP_LEN);
        beg = offset;

        if(u->a[i] == (uint64_t)-1) beg -= GAP_LEN;
        if(end <= beg && i + 1 < u->n && u->a[i+1] == (uint64_t)-1)
        {
            end = beg + GAP_LEN;
        }
        if(rid == i)
        {
            rs = beg;
            re = end;
        }
    }

    for (i = sidx; i < eidx; i++)///keep all hic hits that contain interval we want
    {
        if(get_hit_suid(*hits, i) != get_hit_euid(*hits, i)) continue;
        p0s = get_hit_spos(*hits, i);
        p0e = get_hit_spos_e(*hits, i);
        p1s = get_hit_epos(*hits, i);
        p1e = get_hit_epos_e(*hits, i);

        span_s = MIN(MIN(p0s, p0e), MIN(p1s, p1e));
        span_s = MIN(span_s, ulen-1);
        span_e = MAX(MAX(p0s, p0e), MAX(p1s, p1e));
        span_e = MIN(span_e, ulen-1) + 1;

        //if(span_e - span_s <= ulen*BREAK_CUTOFF)//need it or not?
        {
            if(span_s <= rs && span_e >= re) cnt++;
        }
    }

    // if(cnt != i_cnt) fprintf(stderr, "cnt-%lu, i_cnt-%lu\n", cnt, i_cnt);
    fprintf(stderr, "******cnt-%lu, i_cnt-%lu\n", cnt, i_cnt);
}

int append_sub_utg(ma_ug_t *ug, asg_t *rg, uint64_t uid, uint64_t sidx, uint64_t eidx, uint64_t *rsidx, uint64_t *reidx)
{
    if(eidx <= sidx) return 0;
    uint64_t i, offset;
    ma_utg_t *u = &(ug->u.a[uid]), *p = NULL;
    if(u->a[sidx] == (uint64_t)-1) sidx++;
    if(u->a[eidx-1] == (uint64_t)-1) eidx--;
    if(eidx <= sidx) return 0 ;
    kv_pushp(ma_utg_t, ug->u, &p);
    memset(p, 0, sizeof(*p));
    if(rsidx) (*rsidx) = sidx;
    if(reidx) (*reidx) = eidx;
    p->m = p->n = eidx - sidx;
    MALLOC(p->a, p->m);
    memcpy(p->a, u->a + sidx, p->n*sizeof(uint64_t));
    p->start = p->a[0]>>32;
    p->end = (p->a[p->n-1]>>32)^1;
    p->a[p->n-1] >>= 32; p->a[p->n-1] <<= 32; 
    p->a[p->n-1] += rg->seq[(p->a[p->n-1]>>33)].len;

    p->circ = 0;
    for (i = offset = 0; i < p->n; i++)
    {
        offset += (p->a[i] != (uint64_t)-1? (uint32_t)p->a[i]:GAP_LEN);
    }
    p->len = offset;
    return 1;
}

uint64_t get_utg_len(ma_ug_t *ug)
{   
    uint64_t i, s;
    for (i = s = 0; i < ug->u.n; ++i) 
    {
        if(!ug->u.a[i].a) continue;
        s += ug->u.a[i].len;
    }

    return s;
}
void break_utg_horder(horder_t *h, h_covs *b_points)
{
    if(b_points->n == 0) return;
    kvec_t(uint64_t) join; kv_init(join);
    ma_ug_t *ug = h->ug;
    uint64_t k, l, i, idx, m, pidx, de_u, u_n, oug_n = ug->u.n, dug_n = 0, puid, nuid[2], ps, pe;
    radix_sort_h_cov_s(b_points->a, b_points->a+b_points->n);

    for (k = 1, l = 0; k <= b_points->n; ++k) 
    {
        if (k == b_points->n || b_points->a[k].s != b_points->a[l].s)
        {
            de_u = 0;
            radix_sort_h_cov_e(b_points->a+l, b_points->a+k);
            u_n = ug->u.a[b_points->a[l].s].n;
            for (i = l, pidx = 0; i < k; i++)
            {
                idx = b_points->a[i].e + 1;
                if(idx > pidx && idx - pidx < u_n)
                {
                    if(append_sub_utg(h->ug, h->r_g, b_points->a[l].s, pidx, idx, NULL, NULL))
                    {
                        de_u++;
                        kv_push(uint64_t, join, (b_points->a[l].s<<32)|(ug->u.n-1));
                    }
                }
                pidx = idx;
            }

            idx = u_n;
            if(idx > pidx && idx - pidx < u_n)
            {
                if(append_sub_utg(h->ug, h->r_g, b_points->a[l].s, pidx, idx, NULL, NULL))
                {
                    de_u++;
                    kv_push(uint64_t, join, (b_points->a[l].s<<32)|(ug->u.n-1));
                }
            }

            if(de_u)
            {
                free(ug->u.a[b_points->a[l].s].a); free(ug->u.a[b_points->a[l].s].s);
                memset(&(ug->u.a[b_points->a[l].s]), 0, sizeof(ug->u.a[b_points->a[l].s]));
                dug_n++;
            }

            l = k;
        }
    }

    // fprintf(stderr, "oug_n-%lu, dug_n-%lu\n", oug_n, dug_n);
    // for (i = 0; i < join.n; i++)
    // {
    //     fprintf(stderr, "+puid-%lu, nuid-%u\n", join.a[i]>>32, (uint32_t)join.a[i]);
    // }

    for (i = 0; i < join.n; i++)
    {
        join.a[i] -= dug_n;
    }

    for (i = m = 0; i < ug->u.n; i++)
    {
        if(!ug->u.a[i].a) continue;
        if(i < oug_n) kv_push(uint64_t, join, (i<<32)|(m));
        ug->u.a[m] = ug->u.a[i];
        m++;
    }

    if(m < ug->u.n)
    {
        for (i = m; i < ug->u.n; i++)
        {
            memset(&(ug->u.a[i]), 0, sizeof(ug->u.a[i]));
        }
        ug->u.n = m;
    }

    // for (i = 0; i < join.n; i++)
    // {
    //     fprintf(stderr, "-puid-%lu, nuid-%u\n", join.a[i]>>32, (uint32_t)join.a[i]);
    // }

    oug_n = h->avoid.n;
    radix_sort_ho64(join.a, join.a+join.n);
    for (k = 1, l = 0; k <= join.n; ++k) 
    {
        if (k == join.n || ((join.a[k]>>32) != (join.a[l]>>32))) 
        {
            puid = (join.a[l]>>32);
            nuid[0] = (uint32_t)join.a[l]; 
            nuid[0] <<= 1;

            nuid[1] = (uint32_t)join.a[k-1];
            nuid[1] <<= 1; nuid[1] += 1;

            for (i = 0; i < oug_n; i++)
            {
                ps = h->avoid.a[i]>>32;
                pe = (uint32_t)h->avoid.a[i];

                if((ps>>1) == puid) ps = nuid[ps&1];
                if((pe>>1) == puid) pe = nuid[pe&1];

                h->avoid.a[i] = (ps<<32)|pe;
            }

            if(k - l > 1)
            {
                for (i = l; i + 1 < k; i++)
                {
                    nuid[0] = (uint32_t)join.a[i]; 
                    nuid[0] <<=1; nuid[0] += 1;

                    nuid[1] = (uint32_t)join.a[i+1]; 
                    nuid[1] <<=1;
                    kv_push(uint64_t, h->avoid, (nuid[0]<<32)|(nuid[1]));
                }
            }

            l = k;
        }
    }
    
    radix_sort_ho64(h->avoid.a, h->avoid.a + h->avoid.n);

    // for (i = 0; i < h->avoid.n; i++)
    // {
    //     fprintf(stderr, "break-s-%lu (dir: %lu), break-e-%u (dir: %u)\n", 
    //     h->avoid.a[i]>>33, (h->avoid.a[i]>>32)&1, 
    //     ((uint32_t)h->avoid.a[i])>>1, ((uint32_t)h->avoid.a[i])&1);
    // }

    kv_destroy(join);
}


void update_nus(ma_utg_t *u, ma_ug_t *ug, h_cov_t *a, uint64_t a_n)
{
    uint64_t i, k, offset;
    for (i = offset = 0; i < u->n; i++)
    {   
        for (k = 0; k < a_n; k++)
        {
            if(a[k].s == i)
            {
                a[k].s = offset; a[k].e += offset;
                if(k + 1 == a_n) return;
            } 
            
        }
        offset += (u->a[i] != (uint64_t)-1? (uint32_t)u->a[i]:GAP_LEN);
    }
}

void break_contig(horder_t *h, uint64_t cutoff_s, uint64_t cutoff_e)
{
    uint64_t k, l, i, p0s, p0e, p1s, p1e, ulen, cov_hic, cov_utg, cov_ava, span_s, span_e, cutoff, bs, be, dp;
    kvec_t(uint64_t) b; kv_init(b);
    h_covs cov_buf; kv_init(cov_buf);
    h_covs res; kv_init(res);
    h_covs b_points; kv_init(b_points);
    h_cov_t *p = NULL;
    ma_ug_t *ug = h->ug;
    kvec_pe_hit *hits = &(h->u_hits);
    b_points.n = 0;
    for (k = 1, l = 0; k <= hits->a.n; ++k) 
    {   
        if (k == hits->a.n || (get_hit_suid(*hits, k) != get_hit_suid(*hits, l))) 
        {
            ulen = ug->u.a[get_hit_suid(*hits, l)].len;
            b.n = 0; cov_hic = cov_utg = 0;
            if(ulen >= BREAK_THRES)
            {
                for (i = l; i < k; i++)
                {
                    if(get_hit_suid(*hits, i) != get_hit_euid(*hits, i)) continue;
                    p0s = get_hit_spos(*hits, i);
                    p0e = get_hit_spos_e(*hits, i);
                    p1s = get_hit_epos(*hits, i);
                    p1e = get_hit_epos_e(*hits, i);

                    span_s = MIN(MIN(p0s, p0e), MIN(p1s, p1e));
                    span_s = MIN(span_s, ulen-1);
                    span_e = MAX(MAX(p0s, p0e), MAX(p1s, p1e));
                    span_e = MIN(span_e, ulen-1) + 1;
                    //if(span_e - span_s <= ulen*BREAK_CUTOFF)//need it or not?
                    {
                        kv_push(uint64_t, b, (span_s<<1));
                        kv_push(uint64_t, b, (span_e<<1)|1);
                        cov_hic += (span_e - span_s);
                    }
                }

                radix_sort_ho64(b.a, b.a+b.n);
                cov_utg = get_hic_cov_interval(b.a, b.n, 1, NULL, NULL, NULL);
                cov_ava = (cov_utg? cov_hic/cov_utg:0);
                ///if cov_ava == 0, do nothing or break?

                /*******************************for debug************************************/
                // fprintf(stderr, "\n[M::%s::] utg%.6lul, ulen: %lu, # hic hits: %lu, map cov: %lu, utg cov: %lu, average: %lu\n", 
                // __func__, get_hit_suid(*hits, l)+1, ulen, (uint64_t)(b.n>>1), cov_hic, cov_utg, cov_ava);
                /*******************************for debug************************************/

                res.n = 0;
                for (i = cutoff_s; i <= cutoff_e; i++)
                {
                    if(i == 0) continue;
                    cutoff = cov_ava/i;
                    if(cutoff == 0) continue;
                    get_hic_breakpoint(b.a, b.n, cutoff, &cov_buf, ulen*BREAK_BOUNDARY, ulen - ulen*BREAK_BOUNDARY, &bs, &be);
                    if(bs != (uint64_t)-1 && be != (uint64_t)-1)
                    {
                        kv_pushp(h_cov_t, res, &p);
                        p->s = bs; p->e = be; p->dp = cutoff;
                        /*******************************for debug************************************/
                        // fprintf(stderr, "cutoff: %lu, bs: %lu, be: %lu\n", cutoff, bs, be);
                        /*******************************for debug************************************/
                    }
                }

                if(res.n > 0)
                {
                    get_consensus_break(&res, &cov_buf);
                    /*******************************for debug************************************/
                    // for (i = 0; i < cov_buf.n; i++)
                    // {
                    //     fprintf(stderr, "consensus_break-s: %lu, e: %lu\n", cov_buf.a[i].s, cov_buf.a[i].e);
                    // }
                    /*******************************for debug************************************/
                    get_read_breaks(&(ug->u.a[get_hit_suid(*hits, l)]), h->r_g, &cov_buf, 
                    &res, hits, l, k, ulen, &bs, &dp);
                    if(bs == (uint64_t)-1) fprintf(stderr, "ERROR-read\n");
                    
                    kv_pushp(h_cov_t, b_points, &p);
                    p->s = get_hit_suid(*hits, l); p->e = bs; p->dp = dp;
                    /*******************************for debug************************************/
                    // fprintf(stderr, "consensus_break-rid: %lu, cov: %lu\n", bs, dp);
                    // debug_sub_cov(hits, l, k, ulen, &(ug->u.a[get_hit_suid(*hits, l)]), h->r_g, bs, dp);
                    /*******************************for debug************************************/
                }
                
            }  
            l = k;
        }
    }

    break_utg_horder(h, &b_points);

    kv_destroy(b);
    kv_destroy(cov_buf);
    kv_destroy(res);
    kv_destroy(b_points);
}

asg_arc_t *get_r_edge(asg_t *rg, uint64_t v, uint64_t w, uint64_t *id)
{
    asg_arc_t *av = asg_arc_a(rg, v);
    uint64_t nv = asg_arc_n(rg, v), k;
    if(id) (*id) = (uint64_t)-1;
    for (k = 0; k < nv; k++)
    {
        if(av[k].del) continue;
        if(av[k].v == w) 
        {
            if(id) (*id) = (rg->idx[v]>>32) + k;
            return &(av[k]);
        }
    }
    return NULL;
}

void update_unitig_ends(asg_t *g, asg_arc_t *arc, uint64_t puid, uint64_t nuid_0, uint64_t nuid_1)
{
    uint64_t nv, k, v, idx, ridx;
    

    v = puid<<1;
    idx = g->idx[v]>>32; nv = asg_arc_n(g, v);
    for (k = 0; k < nv; k++)
    {
        arc[idx+k].ul = ((uint32_t)arc[idx+k].ul) + (nuid_0<<32);
        get_r_edge(g, g->arc[idx+k].v^1, v^1, &ridx);
        arc[ridx].v = nuid_0^1;
    }


    v = (puid<<1)+1;
    idx = g->idx[v]>>32; nv = asg_arc_n(g, v);
    for (k = 0; k < nv; k++)
    {
        arc[idx+k].ul = ((uint32_t)arc[idx+k].ul) + (nuid_1<<32);
        get_r_edge(g, g->arc[idx+k].v^1, v^1, &ridx);
        arc[ridx].v = nuid_1^1;
    }
}



int get_switch_ovlp_hits(uint64_t *i, uint64_t pid, uint64_t ls, uint64_t le, u_trans_hit_t *hit, 
h_covs *join, uint64_t *idx)
{
    uint64_t os, oe;
    hit->qSpre = hit->qEpre = hit->qScur = hit->qEcur = hit->qn = (uint32_t)-1;
    hit->tSpre = hit->tEpre = hit->tScur = hit->tEcur = hit->tn = (uint32_t)-1;

    while ((*i) < ((uint32_t)idx[pid]))
    {
        os = join->a[(idx[pid]>>32)+(*i)].s;
        oe = join->a[(idx[pid]>>32)+(*i)].e;
        if(OVL(os, oe, ls, le) == 0)
        {
            (*i)++;
            continue;
        }
        
        hit->qn = (uint32_t)(join->a[(idx[pid]>>32)+(*i)].dp);
        hit->qScur = MAX(os, ls); hit->qEcur = MIN(oe, le);
        hit->qSpre = hit->qScur - os; hit->qEpre = hit->qEcur - os;
        (*i)++;
        return 1;
    }

    return 0;
}

///[ts, te)
void extract_switch_sub(uint32_t i_tScur, uint32_t i_tEcur, uint32_t i_tSpre, uint32_t i_tEpre,
uint32_t tn, kv_u_trans_hit_t* ktb, uint32_t bn, uint32_t rev)
{
    uint32_t i, ovlp, found, beg, end, offS, offE;
    u_trans_hit_t *q = NULL, x;
    for (i = found = 0; i < bn; i++)
    {
        q = &(ktb->a[i]);///for q, already know [qScur, qEcur), [qSpre, qEpre), [tScur, tEcur)

        ovlp = ((MIN(i_tEcur, q->tEcur) > MAX(i_tScur, q->tScur))?
                                        MIN(i_tEcur, q->tEcur) - MAX(i_tScur, q->tScur):0);
        if(found == 1 && ovlp == 0) break;
        if(ovlp > 0) found = 1;
        if(ovlp == 0) continue;

        
        beg = MAX(i_tScur, q->tScur); end = MIN(i_tEcur, q->tEcur);
        offS = beg - q->tScur; offE = q->tEcur - end;
        x.tScur = q->tScur + offS; //beg
        x.tEcur = q->tEcur - offE; //end
        
        
        x.qn = q->qn;
        offS = beg - q->tScur; offE = q->tEcur - end;
        if(rev == 0)
        {
            // x.qSpre = q->qSpre + offS; 
            x.qSpre = q->qSpre + get_offset_adjust(offS, q->tEcur-q->tScur, q->qEpre-q->qSpre);  
            // x.qEpre = q->qEpre - offE;
            x.qEpre = q->qEpre - get_offset_adjust(offE, q->tEcur-q->tScur, q->qEpre-q->qSpre);
        }
        else
        {
            // x.qSpre = q->qSpre + offE; 
            x.qSpre = q->qSpre + get_offset_adjust(offE, q->tEcur-q->tScur, q->qEpre-q->qSpre);
            // x.qEpre = q->qEpre - offS;
            x.qEpre = q->qEpre - get_offset_adjust(offS, q->tEcur-q->tScur, q->qEpre-q->qSpre);
        }

        x.tn = tn; 
        offS = beg - i_tScur; offE = i_tEcur - end;
        // x.tSpre = i_tSpre + offS; 
        x.tSpre = i_tSpre +  get_offset_adjust(offS, i_tEcur-i_tScur, i_tEpre-i_tSpre); 
        // x.tEpre = i_tEpre - offE;
        x.tEpre = i_tEpre - get_offset_adjust(offE, i_tEcur-i_tScur, i_tEpre-i_tSpre); 
        
        kv_push(u_trans_hit_t, *ktb, x);

        // if(x.tSpre >= x.tEpre || x.qSpre >= x.qEpre)
        // {
        //     fprintf(stderr, "\n*********x.qn: %u, x.tn: %u\n", x.qn, x.tn);
        //     fprintf(stderr, "x.qSpre: %u, x.qEpre: %u, x.tSpre: %u, x.tEpre: %u\n", 
        //     x.qSpre, x.qEpre, x.tSpre, x.tEpre);
        //     fprintf(stderr, "q->qScur: %u, q->qEcur: %u, q->qSpre: %u, q->qEpre: %u\n", 
        //     q->qScur, q->qEcur, q->qSpre, q->qEpre);
        //     fprintf(stderr, "q->tScur: %u, q->tEcur: %u, q->tSpre: %u, q->tEpre: %u\n", 
        //     q->tScur, q->tEcur, q->tSpre, q->tEpre);
        //     fprintf(stderr, "i_tScur: %u, i_tEcur: %u, i_tSpre: %u, i_tEpre: %u\n", 
        //     i_tScur, i_tEcur, i_tSpre, i_tEpre);
        // }
    }
}


void extract_novlp(u_trans_t *e, h_covs *join, uint64_t *idx, kv_u_trans_hit_t *kv, kv_u_trans_t *n_trans)
{
    uint64_t i, bn; 
    u_trans_hit_t hit, *kh = NULL;
    u_trans_t *kt = NULL;
    kv->n = 0;

    i = 0;
    while (get_switch_ovlp_hits(&i, e->qn, e->qs, e->qe, &hit, join, idx)) //get [qScur, qEcur), [qSpre, qEpre)
    {
        if(e->rev == 0)
        {
            hit.tScur = e->ts + get_offset_adjust(hit.qScur-e->qs, e->qe-e->qs, e->te-e->ts);  
            hit.tEcur = e->te - get_offset_adjust(e->qe-hit.qEcur, e->qe-e->qs, e->te-e->ts); 
        }
        else
        {
            hit.tScur = e->ts + get_offset_adjust(e->qe-hit.qEcur, e->qe-e->qs, e->te-e->ts); 
            hit.tEcur = e->te - get_offset_adjust(hit.qScur-e->qs, e->qe-e->qs, e->te-e->ts);  
        }
        
        kv_push(u_trans_hit_t, *kv, hit);
    }
    bn = kv->n;

    i = 0;
    while (get_switch_ovlp_hits(&i, e->tn, e->ts, e->te, &hit, join, idx))
    {
        extract_switch_sub(hit.qScur, hit.qEcur, hit.qSpre, hit.qEpre, hit.qn, kv, bn, e->rev);
    }
    
    double x_score, y_score;
    for (i = bn; i < kv->n; i++)
    {
        kh = &(kv->a[i]);
        if(kh->qEpre <= kh->qSpre) continue;
        if(kh->tEpre <= kh->tSpre) continue;
        kv_pushp(u_trans_t, *n_trans, &kt);

        kt->f = e->f; kt->rev = e->rev; kt->del = 0; 
        kt->qn = kh->qn; kt->qs = kh->qSpre; kt->qe = kh->qEpre;
        kt->tn = kh->tn; kt->ts = kh->tSpre; kt->te = kh->tEpre;

        x_score = ((double)(kt->qe-kt->qs)/(double)(e->qe-e->qs))*e->nw;
        y_score = ((double)(kt->te-kt->ts)/(double)(e->te-e->ts))*e->nw;
        kt->nw = MIN(x_score, y_score);
        kt->occ = 0;

        // if(kv->n - bn > 1)
        // {
        //     fprintf(stderr, "-kt-utg%.6ul\tqs(%u)\tqe(%u)\tt-utg%.6ul\tts(%u)\tte(%u)\trev(%u)\tw(%f)\tf(%u)\n", 
        //         kt->qn+1, kt->qs, kt->qe, kt->tn+1, kt->ts, kt->te, kt->rev, kt->nw, kt->f);
        // }
        
    }
}

int get_new_offset(kvec_pe_hit *hits, uint64_t id, uint64_t *index, h_covs *join, uint64_t new_uID_bits)
{
    uint64_t suid, sbeg, send, k, uid, os, oe, ls, le, occ, nuid, noff;
    uint64_t euid, ebeg, eend;
    resolve_hit(hits->a.a[id].s, hits->a.a[id].len>>32, hits->uID_bits, 
                                                                hits->pos_mode, &suid, &sbeg, &send);
    uid = suid; ls = sbeg; le = send; occ = 0; nuid = noff = (uint64_t)-1; 
    for (k = 0; k < ((uint32_t)index[uid]); k++)
    {
        os = join->a[(index[uid]>>32)+k].s;
        oe = join->a[(index[uid]>>32)+k].e;
        if(ls >= os && le <= oe)
        {
            occ++;
            nuid = (uint32_t)join->a[(index[uid]>>32)+k].dp;
            noff = os;
        }
    }
    if(occ != 1) return 0;
    noff = get_hit_spos(*hits, id) - noff;
    nuid <<= (64 - new_uID_bits - 1);
    hits->a.a[id].s >>= 63; hits->a.a[id].s <<= 63; 
    hits->a.a[id].s += nuid + noff;

    
    resolve_hit(hits->a.a[id].e, (uint32_t)hits->a.a[id].len, hits->uID_bits, 
                                                                hits->pos_mode, &euid, &ebeg, &eend);
    uid = euid; ls = ebeg; le = eend; occ = 0; nuid = noff = (uint64_t)-1;
    for (k = 0; k < ((uint32_t)index[uid]); k++)
    {
        os = join->a[(index[uid]>>32)+k].s;
        oe = join->a[(index[uid]>>32)+k].e;
        if(ls >= os && le <= oe)
        {
            occ++;
            nuid = (uint32_t)join->a[(index[uid]>>32)+k].dp;
            noff = os;
        }
    }
    if(occ != 1) return 0;
    noff = get_hit_epos(*hits, id) - noff;
    nuid <<= (64 - new_uID_bits - 1);
    hits->a.a[id].e >>= 63; hits->a.a[id].e <<= 63; 
    hits->a.a[id].e += nuid + noff;
    return 1;
}

void verbose_misjoin(uint64_t *a, uint64_t a_n)
{
    radix_sort_ho64(a, a + a_n);
    uint64_t k, l, i, len, s, C[2], N[2];
    C[0] = C[1] = N[0] = N[1] = 0;
    for (k = 1, l = 0; k <= a_n; ++k) 
    {   
        if (k == a_n || (a[k]>>63) != (a[l]>>63))
        {
            for (i = l, s = 0; i < k; i++) s += ((a[i]<<1)>>1);
            len = s;

            i = k; s = 0;
            while (i > l)
            {
                i--;
                s += ((a[i]<<1)>>1);
                if(s >= (len>>1))
                {
                    N[a[l]>>63] = ((a[i]<<1)>>1);
                    C[a[l]>>63] = k - l;
                    // fprintf(stderr, "[M::%s::] N50: %lu\n", __func__, b.a[i]);
                    break;
                }
            }
            
            l = k;
        }
    }

    fprintf(stderr, "[M::stat] # misjoined unitigs: %lu (N50: %lu); # corrected unitigs: %lu (N50: %lu)\n", 
    C[0], N[0], C[1], N[1]);
}

void break_phasing_utg(ma_ug_t *ug, asg_t *rg, kvec_pe_hit *hits, kv_u_trans_t *k_trans, h_covs *b_points)
{
    if(b_points->n == 0) return;
    kv_u_trans_hit_t kv; kv_init(kv);
    asg_arc_t *rt = NULL, *ut = NULL;
    asg_t *utg = NULL;
    h_cov_t *p = NULL;
    h_covs join; kv_init(join);
    uint64_t k, l, i, idx, m, pidx, de_u, u_n, oug_n = ug->u.n, dug_n = 0, puid, pn, rsi, prid, nrid, x, y, *index = NULL;
    kv_u_trans_t *n_trans = NULL; CALLOC(n_trans, 1);
    kvec_t(uint64_t) dbg_N50; kv_init(dbg_N50);
    /*******************************for debug************************************/
    // debug_phasing_t *dbp = init_debug_phasing(ug, hits, k_trans, b_points);
    /*******************************for debug************************************/
    radix_sort_h_cov_s(b_points->a, b_points->a+b_points->n);
    for (k = 1, l = 0; k <= b_points->n; ++k) 
    {
        if (k == b_points->n || b_points->a[k].s != b_points->a[l].s)
        {
            de_u = 0; pn = join.n;
            radix_sort_h_cov_e(b_points->a+l, b_points->a+k);
            u_n = ug->u.a[b_points->a[l].s].n;
            for (i = l, pidx = 0; i < k; i++)
            {
                idx = b_points->a[i].e + 1;
                if(idx > pidx && idx - pidx < u_n)
                {
                    // fprintf(stderr, "sa-utg%.6lul, pidx: %lu, idx: %lu\n", b_points->a[l].s+1, pidx, idx);
                    if(append_sub_utg(ug, rg, b_points->a[l].s, pidx, idx, &rsi, NULL))
                    {
                        de_u++;
                        kv_pushp(h_cov_t, join, &p);
                        p->dp = (b_points->a[l].s<<32)|(ug->u.n-1);
                        p->s = rsi; p->e = ug->u.a[ug->u.n-1].len;
                        kv_push(uint64_t, dbg_N50, (uint64_t)(ug->u.a[ug->u.n-1].len)+((uint64_t)(1)<<63));
                    }
                }
                pidx = idx;
            }

            idx = u_n;
            if(idx > pidx && idx - pidx < u_n)
            {
                // fprintf(stderr, "sa-utg%.6lul, pidx: %lu, idx: %lu\n", b_points->a[l].s+1, pidx, idx);
                if(append_sub_utg(ug, rg, b_points->a[l].s, pidx, idx, &rsi, NULL))
                {
                    de_u++;
                    kv_pushp(h_cov_t, join, &p);
                    p->dp = (b_points->a[l].s<<32)|(ug->u.n-1);
                    p->s = rsi; p->e = ug->u.a[ug->u.n-1].len;
                    kv_push(uint64_t, dbg_N50, (uint64_t)(ug->u.a[ug->u.n-1].len)+((uint64_t)(1)<<63));
                }
            }

            if(de_u)
            {
                kv_push(uint64_t, dbg_N50, ug->u.a[b_points->a[l].s].len);
                update_nus(&(ug->u.a[b_points->a[l].s]), ug, join.a + pn, join.n - pn);
                free(ug->u.a[b_points->a[l].s].a); free(ug->u.a[b_points->a[l].s].s);
                memset(&(ug->u.a[b_points->a[l].s]), 0, sizeof(ug->u.a[b_points->a[l].s]));
                dug_n++;
            }

            l = k;
        }
    }
    // fprintf(stderr, "oug_n-%lu, dug_n-%lu\n", oug_n, dug_n);
    // for (i = 0; i < join.n; i++)
    // {
    //     fprintf(stderr, "+p-utg%.6lul, n-utg%.6ul\n", (join.a[i].dp>>32)+1, ((uint32_t)join.a[i].dp)+1);
    // }

    for (i = 0; i < join.n; i++)
    {
        join.a[i].dp -= dug_n;
    }
    for (i = m = 0; i < ug->u.n; i++)
    {
        if(!ug->u.a[i].a) continue;
        if(i < oug_n)
        {
            kv_pushp(h_cov_t, join, &p);
            p->dp = (i<<32)|(m);
            p->s = 0; p->e = ug->u.a[i].len;
        }

        ug->u.a[m] = ug->u.a[i];
        m++;
    }
    if(m < ug->u.n)
    {
        for (i = m; i < ug->u.n; i++)
        {
            memset(&(ug->u.a[i]), 0, sizeof(ug->u.a[i]));
        }
        ug->u.n = m;
    }
    radix_sort_h_cov_dp(join.a, join.a+join.n);
    CALLOC(index, ug->u.n);
    // for (i = 0; i < join.n; i++)
    // {
    //     fprintf(stderr, "+p-utg%.6lul (len: %u), n-utg%.6ul, s-%lu, e-%lu\n", 
    //             (join.a[i].dp>>32)+1, ug->g->seq[(join.a[i].dp>>32)].len, ((uint32_t)join.a[i].dp)+1, join.a[i].s, join.a[i].e);
    // }
    utg = asg_init();
    utg->n_arc = utg->m_arc = ug->g->n_arc;
    MALLOC(utg->arc, utg->n_arc);
    memcpy(utg->arc, ug->g->arc, utg->n_arc*sizeof(asg_arc_t));
    for (k = 1, l = 0; k <= join.n; ++k) 
    {
        if (k == join.n || ((join.a[k].dp>>32) != (join.a[l].dp>>32))) 
        {
            puid = (join.a[l].dp>>32);
            index[puid] = (l<<32) | (k-l);

            for (i = l; i < k; i++)
            {
                asg_seq_set(utg, (uint32_t)join.a[i].dp, ug->u.a[(uint32_t)join.a[i].dp].len, 0);
                utg->seq[(uint32_t)join.a[i].dp].c = ug->g->seq[puid].c;

                if(i + 1 >= k) continue;
                prid = ug->u.a[(uint32_t)join.a[i].dp].a[ug->u.a[(uint32_t)join.a[i].dp].n-1]>>32;
                nrid = ug->u.a[(uint32_t)join.a[i+1].dp].a[0]>>32;

                rt = get_r_edge(rg, prid, nrid, NULL);
                ut = asg_arc_pushp(utg); 
                x = (uint32_t)join.a[i].dp; x <<= 1;
                y = (uint32_t)join.a[i+1].dp; y <<= 1;
			    ut->ol = rt->ol, ut->del = 0;
			    ut->ul = (uint64_t)x<<32 | (ug->u.a[x>>1].len - ut->ol);
			    ut->v = y;

                rt = get_r_edge(rg, nrid^1, prid^1, NULL);
                ut = asg_arc_pushp(utg); 
                x = (uint32_t)join.a[i+1].dp; x <<= 1; x++;
                y = (uint32_t)join.a[i].dp; y <<= 1; y++;
                ut->ol = rt->ol, ut->del = 0;
                ut->ul = (uint64_t)x<<32 | (ug->u.a[x>>1].len - ut->ol);
			    ut->v = y;
            }
            update_unitig_ends(ug->g, utg->arc, puid, ((uint32_t)join.a[k-1].dp)<<1, (((uint32_t)join.a[l].dp)<<1)+1);
            l = k;
        }
    }
    asg_cleanup(utg);
    asg_destroy(ug->g);
    ug->g = utg;

    for (k = 0; k < k_trans->n; k++)
    {
        extract_novlp(&(k_trans->a[k]), &join, index, &kv, n_trans);
    }
    free(k_trans->a); free(k_trans->idx.a); memset(k_trans, 0, sizeof(*k_trans));
    k_trans->a = n_trans->a; k_trans->n = n_trans->n; k_trans->m = n_trans->m; 
    free(n_trans);
    clean_u_trans_t_idx(k_trans, ug, rg);
    uint64_t uID_bits, pos_mode;
    for (uID_bits=1; (uint64_t)(1<<uID_bits)<(uint64_t)ug->u.n; uID_bits++);
    pos_mode = ((uint64_t)-1)>>(uID_bits+1);
    
    for (k = m = 0; k < hits->a.n; k++)
    {
        if(get_new_offset(hits, k, index, &join, uID_bits))
        {
            hits->a.a[m] = hits->a.a[k];
            m++;
        }
    }
    // fprintf(stderr, "hits->a.n: %u, m: %lu\n", (uint32_t)hits->a.n, m);
    hits->a.n = m;  

    dedup_hits(hits, 0); 
    hits->uID_bits = uID_bits; hits->pos_mode = pos_mode;
    free(hits->idx.a); hits->idx.a = NULL; hits->idx.n = hits->idx.m = 0;
    free(hits->occ.a); hits->occ.a = NULL; hits->occ.n = hits->occ.m = 0;

    /*******************************for debug************************************/
    // debug_debug_phasing_t(dbp, ug, hits, k_trans, &join, index);
    // destory_debug_phasing_t(&dbp);
    /*******************************for debug************************************/
    verbose_misjoin(dbg_N50.a, dbg_N50.n);
    kv_destroy(join); kv_destroy(kv); free(index); kv_destroy(dbg_N50);
}

///min_ulen = BREAK_THRES
///boundaryRate = BREAK_BOUNDARY
void update_switch_unitig(ma_ug_t *ug, asg_t *rg, kvec_pe_hit *hits, kv_u_trans_t *k_trans, uint64_t cutoff_s, uint64_t cutoff_e,
uint64_t min_ulen, double boundaryRate)
{
    uint64_t k, l, i, p0s, p0e, p1s, p1e, ulen, cov_hic, cov_utg, cov_ava, span_s, span_e, cutoff, bs, be, dp;
    kvec_t(uint64_t) b; kv_init(b);
    h_covs cov_buf; kv_init(cov_buf);
    h_covs res; kv_init(res);
    h_covs b_points; kv_init(b_points);
    h_cov_t *p = NULL;
    radix_sort_pe_hit_idx_hn1(hits->a.a, hits->a.a + hits->a.n);
    b_points.n = 0;
    for (k = 1, l = 0; k <= hits->a.n; ++k) 
    {   
        if (k == hits->a.n || (get_hit_suid(*hits, k) != get_hit_suid(*hits, l))) 
        {
            ulen = ug->u.a[get_hit_suid(*hits, l)].len;
            b.n = 0; cov_hic = cov_utg = 0;
            if(ulen >= min_ulen)
            {
                for (i = l; i < k; i++)
                {
                    if(get_hit_suid(*hits, i) != get_hit_euid(*hits, i)) continue;
                    p0s = get_hit_spos(*hits, i);
                    p0e = get_hit_spos_e(*hits, i);
                    p1s = get_hit_epos(*hits, i);
                    p1e = get_hit_epos_e(*hits, i);

                    span_s = MIN(MIN(p0s, p0e), MIN(p1s, p1e));
                    span_s = MIN(span_s, ulen-1);
                    span_e = MAX(MAX(p0s, p0e), MAX(p1s, p1e));
                    span_e = MIN(span_e, ulen-1) + 1;
                    //if(span_e - span_s <= ulen*BREAK_CUTOFF)//need it or not?
                    {
                        kv_push(uint64_t, b, (span_s<<1));
                        kv_push(uint64_t, b, (span_e<<1)|1);
                        cov_hic += (span_e - span_s);
                    }
                }

                radix_sort_ho64(b.a, b.a+b.n);
                cov_utg = get_hic_cov_interval(b.a, b.n, 1, NULL, NULL, NULL);
                cov_ava = (cov_utg? cov_hic/cov_utg:0);
                ///if cov_ava == 0, do nothing or break?

                /*******************************for debug************************************/
                // fprintf(stderr, "\n[M::%s::] utg%.6lul, ulen: %lu, # hic hits: %lu, map cov: %lu, utg cov: %lu, average: %lu\n", 
                // __func__, get_hit_suid(*hits, l)+1, ulen, (uint64_t)(b.n>>1), cov_hic, cov_utg, cov_ava);
                /*******************************for debug************************************/

                res.n = 0;
                for (i = cutoff_s; i <= cutoff_e; i++)
                {
                    if(i == 0) continue;
                    cutoff = cov_ava/i;
                    if(cutoff == 0) continue;
                    get_hic_breakpoint(b.a, b.n, cutoff, &cov_buf, ulen*boundaryRate, ulen - ulen*boundaryRate, &bs, &be);
                    if(bs != (uint64_t)-1 && be != (uint64_t)-1)
                    {
                        kv_pushp(h_cov_t, res, &p);
                        p->s = bs; p->e = be; p->dp = cutoff;
                        /*******************************for debug************************************/
                        // fprintf(stderr, "cutoff: %lu, bs: %lu, be: %lu\n", cutoff, bs, be);
                        /*******************************for debug************************************/
                    }
                }

                if(res.n > 0)
                {
                    get_consensus_break(&res, &cov_buf);
                    /*******************************for debug************************************/
                    // for (i = 0; i < cov_buf.n; i++)
                    // {
                    //     fprintf(stderr, "consensus_break-s: %lu, e: %lu\n", cov_buf.a[i].s, cov_buf.a[i].e);
                    // }
                    /*******************************for debug************************************/
                    get_read_breaks(&(ug->u.a[get_hit_suid(*hits, l)]), rg, &cov_buf, 
                    &res, hits, l, k, ulen, &bs, &dp);
                    if(bs == (uint64_t)-1) fprintf(stderr, "ERROR-read\n");
                    
                    if(bs > 0 && bs < ug->u.a[get_hit_suid(*hits, l)].n)
                    {
                        kv_pushp(h_cov_t, b_points, &p);
                        p->s = get_hit_suid(*hits, l); p->e = bs; p->dp = dp;
                        /*******************************for debug************************************/
                        // fprintf(stderr, "\n[M::%s::] utg%.6lul, consensus_break-rid: %lu, cov: %lu, un: %u\n", 
                        // __func__, get_hit_suid(*hits, l)+1, bs, dp, ug->u.a[get_hit_suid(*hits, l)].n);
                        // debug_sub_cov(hits, l, k, ulen, &(ug->u.a[get_hit_suid(*hits, l)]), h->r_g, bs, dp);
                        /*******************************for debug************************************/
                    }
                    
                }
                
            }  
            l = k;
        }
    }
    // break_utg_horder(h, &b_points);
    break_phasing_utg(ug, rg, hits, k_trans, &b_points);

    kv_destroy(b);
    kv_destroy(cov_buf);
    kv_destroy(res);
    kv_destroy(b_points);
}

void get_Ns(ma_utg_t *u, h_covs *Ns)
{
    uint64_t i, offset;
    h_cov_t *p = NULL;
    Ns->n = 0;
    for (i = offset = 0; i < u->n; i++)
    {
        if(u->a[i] == (uint64_t)-1)
        {
            kv_pushp(h_cov_t, *Ns, &p);
            p->s = offset;
            p->e = offset + GAP_LEN;
            p->dp = i;
        }
        offset += (u->a[i] != (uint64_t)-1? (uint32_t)u->a[i]:GAP_LEN);
    }
}

void debug_sub_cov(kvec_pe_hit *hits, uint64_t sidx, uint64_t eidx, uint64_t ulen, 
uint64_t rs, uint64_t re, uint64_t limit_s, uint64_t limit_e, int unique_only)
{
    uint64_t i, p0s, p0e, p1s, p1e, span_s, span_e, cnt = 0, cnt_no_lim = 0;


    for (i = sidx; i < eidx; i++)///keep all hic hits that contain interval we want
    {
        if(get_hit_suid(*hits, i) != get_hit_euid(*hits, i)) continue;
        if(unique_only && hits->a.a[i].id == 0) continue;
        p0s = get_hit_spos(*hits, i);
        p0e = get_hit_spos_e(*hits, i);
        p1s = get_hit_epos(*hits, i);
        p1e = get_hit_epos_e(*hits, i);

        span_s = MIN(MIN(p0s, p0e), MIN(p1s, p1e));
        span_s = MIN(span_s, ulen-1);
        span_e = MAX(MAX(p0s, p0e), MAX(p1s, p1e));
        span_e = MIN(span_e, ulen-1) + 1;

        //if(span_e - span_s <= ulen*BREAK_CUTOFF)//need it or not?
        if(span_s <= rs && span_e >= re)
        {
            cnt_no_lim++;
            if(span_s >= limit_s && span_e <= limit_e) cnt++;
        } 
    }

    fprintf(stderr, "******cnt-%lu, cnt_no_lim-%lu, rs-%lu, re-%lu, limit_s-%lu, limit_e-%lu\n", cnt, cnt_no_lim, rs, re, limit_s, limit_e);
}

uint64_t get_sub_cov(kvec_pe_hit *hits, uint64_t sidx, uint64_t eidx, uint64_t ulen, 
uint64_t rs, uint64_t re, uint64_t limit_s, uint64_t limit_e, int unique_only)
{
    uint64_t i, p0s, p0e, p1s, p1e, span_s, span_e, cnt = 0;


    for (i = sidx; i < eidx; i++)///keep all hic hits that contain interval we want
    {
        if(get_hit_suid(*hits, i) != get_hit_euid(*hits, i)) continue;
        if(unique_only && hits->a.a[i].id == 0) continue;
        p0s = get_hit_spos(*hits, i);
        p0e = get_hit_spos_e(*hits, i);
        p1s = get_hit_epos(*hits, i);
        p1e = get_hit_epos_e(*hits, i);

        span_s = MIN(MIN(p0s, p0e), MIN(p1s, p1e));
        span_s = MIN(span_s, ulen-1);
        span_e = MAX(MAX(p0s, p0e), MAX(p1s, p1e));
        span_e = MIN(span_e, ulen-1) + 1;

        //if(span_e - span_s <= ulen*BREAK_CUTOFF)//need it or not?
        if(span_s <= rs && span_e >= re)
        {
            if(span_s >= limit_s && span_e <= limit_e) cnt++;
        } 
    }

    return cnt;
}

int detect_lowNs(kvec_pe_hit *hit, uint64_t sHit, uint64_t eHit, kvec_t_u64_warp *b, 
h_cov_t *Np, uint64_t len, uint64_t cutoff_s, uint64_t cutoff_e, uint64_t force_cutoff, uint64_t force_cutoff_cov,
h_covs *res, h_covs *cov_buf, h_covs *b_points, uint64_t local_bound, int unique_only)
{
    uint64_t cov_hic, cov_utg, cov_ava, i, p0s, p0e, p1s, p1e, span_s, span_e, cutoff, bs, be, occ = 0;
    uint64_t sPos, ePos, min_cutoff;
    h_cov_t *p = NULL;
    b->a.n = 0; cov_hic = cov_utg = 0;
    sPos = (Np->s>=local_bound? Np->s-local_bound:0);
    ePos = (Np->e+local_bound<=len? Np->e+local_bound:len);
    for (i = sHit; i < eHit; i++)
    {
        if(get_hit_suid(*hit, i) != get_hit_euid(*hit, i)) continue;
        
        p0s = get_hit_spos(*hit, i);
        p0e = get_hit_spos_e(*hit, i);
        p1s = get_hit_epos(*hit, i);
        p1e = get_hit_epos_e(*hit, i);

        span_s = MIN(MIN(p0s, p0e), MIN(p1s, p1e));
        span_s = MIN(span_s, len-1);
        span_e = MAX(MAX(p0s, p0e), MAX(p1s, p1e));
        span_e = MIN(span_e, len-1) + 1;
        //if(span_e - span_s <= ulen*BREAK_CUTOFF)//need it or not?
        {
            if(span_s >= sPos && span_e <= ePos)
            {
                occ++;
                if(unique_only && hit->a.a[i].id == 0) continue;
                kv_push(uint64_t, b->a, (span_s<<1));
                kv_push(uint64_t, b->a, (span_e<<1)|1);
                cov_hic += (span_e - span_s);
            }
        }
    }

    

    radix_sort_ho64(b->a.a, b->a.a+b->a.n);
    cov_utg = get_hic_cov_interval(b->a.a, b->a.n, 1, NULL, NULL, NULL);
    cov_ava = (cov_utg? cov_hic/cov_utg:0);
    if(force_cutoff != (uint64_t)-1 || force_cutoff_cov != (uint64_t)-1)
    {
        min_cutoff = get_sub_cov(hit, sHit, eHit, len, Np->s, Np->e, sPos, ePos, unique_only);
        if((force_cutoff != (uint64_t)-1 && min_cutoff <= (cov_ava/force_cutoff)) || 
                (force_cutoff_cov != (uint64_t)-1 && min_cutoff <= force_cutoff_cov))
        {
            kv_pushp(h_cov_t, *b_points, &p);
            p->s = get_hit_suid(*hit, sHit); p->e = Np->dp; p->dp = 0;
            return 1;
        }
    }
    
    ///if cov_ava == 0, do nothing or break?
    /*******************************for debug************************************/
    // fprintf(stderr, "\n[M::%s::] utg%.6lul, ulen: %lu, # hic hits: %lu, map cov: %lu, utg cov: %lu, average: %lu\n", 
    // __func__, get_hit_suid(*hit, sHit)+1, len, (uint64_t)(b->a.n>>1), cov_hic, cov_utg, cov_ava);
    /*******************************for debug************************************/

    // fprintf(stderr, "[M::%s::] sPos: %lu, ePos: %lu, # hits: %lu, # non-unique hits: %lu, cov_hic: %lu, cov_utg: %lu, cov_ava: %lu\n", 
    //                     __func__, sPos, ePos, (uint64_t)(b->a.n>>1), occ, cov_hic, cov_utg, cov_ava);
    // debug_sub_cov(hit, sHit, eHit, len, Np->s, Np->e, sPos, ePos, unique_only);
    res->n = 0;
    for (i = cutoff_s; i <= cutoff_e; i++)
    {
        if(i == 0) continue;
        cutoff = cov_ava/i;
        if(cutoff == 0) continue;
        get_hic_breakpoint(b->a.a, b->a.n, cutoff, cov_buf, sPos + (ePos-sPos)*BREAK_BOUNDARY, 
        ePos - (ePos-sPos)*BREAK_BOUNDARY, &bs, &be);
        
        if(bs != (uint64_t)-1 && be != (uint64_t)-1)
        {
            kv_pushp(h_cov_t, *res, &p);
            p->s = bs; p->e = be; p->dp = cutoff;
            /*******************************for debug************************************/
            // fprintf(stderr, "cutoff: %lu, bs: %lu, be: %lu\n", cutoff, bs, be);
            /*******************************for debug************************************/
        }
    }

    if(res->n > 0)
    {
        get_consensus_break(res, cov_buf);
        /*******************************for debug************************************/
        // for (i = 0; i < cov_buf->n; i++)
        // {
        //     fprintf(stderr, "consensus_break-s: %lu, e: %lu\n", cov_buf->a[i].s, cov_buf->a[i].e);
        // }
        /*******************************for debug************************************/
        for (i = 0, min_cutoff = (uint64_t)-1; i < cov_buf->n; i++)
        {
            if(cov_buf->a[i].s<=Np->s && cov_buf->a[i].e>=Np->e)
            {
                break;
            }
            min_cutoff = MIN(min_cutoff, cov_buf->a[i].dp);
        }
        if(i < cov_buf->n || 
                get_sub_cov(hit, sHit, eHit, len, Np->s, Np->e, sPos, ePos, unique_only) <= min_cutoff)
        {
            kv_pushp(h_cov_t, *b_points, &p);
            p->s = get_hit_suid(*hit, sHit); p->e = Np->dp; p->dp = 0;
            return 1;
            /*******************************for debug************************************/
            // fprintf(stderr, "consensus_break-rid: %lu\n", p->e);
            /*******************************for debug************************************/
        }
    }
    return 0;
}

uint64_t break_scaffold(horder_t *h, uint64_t cutoff_s, uint64_t cutoff_e, uint64_t force_cutoff, uint64_t force_cutoff_cov,
uint64_t local_bound, int unique_only, h_covs *r_b_points)
{
    uint64_t k, l, i, ulen;
    kvec_t_u64_warp b; kv_init(b.a);
    h_covs cov_buf; kv_init(cov_buf);
    h_covs res; kv_init(res);
    h_covs *b_points = NULL; 
    if(r_b_points) b_points = r_b_points;
    else CALLOC(b_points, 1);
    b_points->n = 0;
    h_covs Ns; kv_init(Ns);
    ma_ug_t *ug = h->ug;
    kvec_pe_hit *hits = &(h->u_hits);
    for (k = 1, l = 0; k <= hits->a.n; ++k) 
    {   
        if (k == hits->a.n || (get_hit_suid(*hits, k) != get_hit_suid(*hits, l))) 
        {
            ulen = ug->u.a[get_hit_suid(*hits, l)].len;
            Ns.n = 0;
            if(ulen >= BREAK_THRES)
            {
                get_Ns(&(ug->u.a[get_hit_suid(*hits, l)]), &Ns);
                // fprintf(stderr, "\n[M::%s::] utg%.6lul, ulen: %lu, # Ns: %lu\n", 
                //         __func__, get_hit_suid(*hits, l)+1, ulen, (uint64_t)(Ns.n));
                if(Ns.n)
                {
                    for (i = 0; i < Ns.n; i++)
                    {
                        if(detect_lowNs(hits, l, k, &b, &(Ns.a[i]), ulen, cutoff_s, cutoff_e, force_cutoff,
                        force_cutoff_cov, &res, &cov_buf, b_points, local_bound, unique_only) && r_b_points){
                            b_points->a[b_points->n-1].dp = i;
                        }
                    }
                }
            }  
            l = k;
        }
    }

    if(!r_b_points){
        break_utg_horder(h, b_points);
        kv_destroy(*b_points);
    }
    
    kv_destroy(b.a);
    kv_destroy(cov_buf);
    kv_destroy(res);
    kv_destroy(Ns);

    return b_points->n;
}

uint64_t break_scaffold_mean(horder_t *h, uint64_t cutoff_s, uint64_t cutoff_e, uint64_t force_cutoff, uint64_t force_cutoff_cov,
uint64_t local_bound, int unique_only, h_covs *r_b_points)
{
    uint64_t k, l, i, ulen;
    kvec_t_u64_warp b; kv_init(b.a);
    h_covs cov_buf; kv_init(cov_buf);
    h_covs res; kv_init(res);
    h_covs *b_points = NULL; 
    if(r_b_points) b_points = r_b_points;
    else CALLOC(b_points, 1);
    b_points->n = 0;
    h_covs Ns; kv_init(Ns);
    ma_ug_t *ug = h->ug;
    kvec_pe_hit *hits = &(h->u_hits);
    for (k = 1, l = 0; k <= hits->a.n; ++k) 
    {   
        if (k == hits->a.n || (get_hit_suid(*hits, k) != get_hit_suid(*hits, l))) 
        {
            ulen = ug->u.a[get_hit_suid(*hits, l)].len;
            Ns.n = 0;
            // if(ulen >= BREAK_THRES)
            {
                get_Ns(&(ug->u.a[get_hit_suid(*hits, l)]), &Ns);
                if(Ns.n)
                {
                    for (i = 0; i < Ns.n; i++)
                    {
                        if(detect_lowNs(hits, l, k, &b, &(Ns.a[i]), ulen, cutoff_s, cutoff_e, force_cutoff,
                        force_cutoff_cov, &res, &cov_buf, b_points, local_bound, unique_only) && r_b_points){
                            b_points->a[b_points->n-1].dp = i;
                        }
                    }
                }
            }  
            l = k;
        }
    }

    if(!r_b_points){
        break_utg_horder(h, b_points);
        kv_destroy(*b_points);
    }
    
    kv_destroy(b.a);
    kv_destroy(cov_buf);
    kv_destroy(res);
    kv_destroy(Ns);

    return b_points->n;
}


void generate_haplotypes(horder_t *h, ug_opt_t *opt)
{
    uint64_t i;
    ma_ug_t *ug_1 = NULL, *ug_2 = NULL;

    ug_1 = get_trio_unitig_graph(h->r_g, FATHER, opt);
    ug_2 = get_trio_unitig_graph(h->r_g, MOTHER, opt);
    // kv_push(uint64_t, h->occ, ug_1->u.n);
    // kv_push(uint64_t, h->occ, ug_2->u.n);

    h->ug = ug_1; 
    
    ///update unitigs
    ma_utg_t *pu = NULL;
    for (i = 0; i < ug_2->u.n; i++)
    {
        kv_pushp(ma_utg_t, h->ug->u, &pu);
        *pu = ug_2->u.a[i];
        ug_2->u.a[i].a = NULL;
        ug_2->u.a[i].s = NULL;
    }

    ug_1 = NULL;
    ma_ug_destroy(ug_2);
    asg_destroy(h->ug->g);
    h->ug->g = NULL;

    // MALLOC(h->hf.a, h->ug->u.n); h->hf.n = h->hf.m = h->ug->u.n;
    // for (i = 0; i < h->ug->u.n; i++) h->hf.a[i] = (i < h->occ.a[0]? 0 : 1);
}

osg_t *osg_init(void)
{
	return (osg_t*)calloc(1, sizeof(osg_t));
}

void osg_destroy(osg_t *g)
{
	if (g == 0) return;
	free(g->seq); free(g->idx); free(g->arc); 
    free(g);
}

void osg_seq_set(osg_t *g, int sid, int del)
{
	///just malloc size
	if (sid >= (int)g->m_seq) {
		g->m_seq = sid + 1;
		kv_roundup32(g->m_seq);
		g->seq = (osg_seq_t*)realloc(g->seq, g->m_seq * sizeof(osg_seq_t));
	}

	if (sid >= g->n_seq) g->n_seq = sid + 1;
    g->seq[sid].del = !!del;
}

static inline osg_arc_t *osg_arc_pushp(osg_t *g)
{
	if (g->n_arc == g->m_arc) {
		g->m_arc = g->m_arc? g->m_arc<<1 : 16;
		g->arc = (osg_arc_t*)realloc(g->arc, g->m_arc * sizeof(osg_arc_t));
	}
	return &g->arc[g->n_arc++];
}

void osg_arc_rm(osg_t *g)
{
	uint32_t e, n;
	///just clean arc requiring: 1. arc it self must be available 2. both the query and target are available
	for (e = n = 0; e < g->n_arc; ++e) {
		//u and v is the read id
		uint32_t u = g->arc[e].u, v = g->arc[e].v;
		if (!g->arc[e].del && !g->seq[u>>1].del && !g->seq[v>>1].del)
			g->arc[n++] = g->arc[e];
	}
	if (n < g->n_arc) { // arc index is out of sync
		if (g->idx) free(g->idx);
		g->idx = 0;
	}
	g->n_arc = n;
}

uint64_t *osg_arc_index_core(size_t max_seq, size_t n, const osg_arc_t *a)
{
	size_t i, last;
	uint64_t *idx;
	idx = (uint64_t*)calloc(max_seq * 2, 8);

	for (i = 1, last = 0; i <= n; ++i)
		if (i == n || a[i-1].u != a[i].u)
			idx[a[i-1].u] = (uint64_t)last<<32 | (i - last), last = i;
	return idx;
}

void osg_arc_index(osg_t *g)
{
	if (g->idx) free(g->idx);
	g->idx = osg_arc_index_core(g->n_seq, g->n_arc, g->arc);
}

void osg_cleanup(osg_t *g)
{
	osg_arc_rm(g);
	if (!g->is_srt) {
        radix_sort_osg(g->arc, g->arc + g->n_arc);
		g->is_srt = 1;
	}

	if (g->idx == 0) osg_arc_index(g);
}

double get_max_weight(uint32_t u, uint32_t v, osg_t *g)
{
    double max = 0;
    uint32_t i, nv;
	osg_arc_t *av = NULL;

    nv = asg_arc_n(g, u);
    av = asg_arc_a(g, u);
    for (i = 0; i < nv; i++)
    {
        if(av[i].v == v) continue;
        max = MAX(max, av[i].w);
    }
    
    nv = asg_arc_n(g, v);
    av = asg_arc_a(g, v);
    for (i = 0; i < nv; i++)
    {
        if(av[i].v == u) continue;
        max = MAX(max, av[i].w);
    }

    return max;
}

dens_idx_t *build_interval_idx(kvec_pe_hit *hits, ma_ug_t *ug)
{
    uint64_t i, k, l, p0s, p0e, p1s, p1e, suid, euid, slen, elen, span_s, span_e;
    dens_idx_t *idx = NULL;
    CALLOC(idx, 1);

    for (i = 0; i < hits->a.n; i++)
    {
        if(!hits->a.a[i].id) continue;
        suid = get_hit_suid(*hits, i);
        euid = get_hit_euid(*hits, i);
        slen = ug->u.a[suid].len;
        elen = ug->u.a[euid].len;

        p0s = get_hit_spos(*hits, i);
        p0e = get_hit_spos_e(*hits, i);
        span_s = MIN(p0s, p0e);
        span_s = MIN(span_s, slen-1);
        span_e = MAX(p0s, p0e);
        span_e =  MIN(span_e, slen-1);
        span_s = ((span_s+span_e)>>1);
        kv_push(uint64_t, idx->pos, (suid<<32)|span_s);


        p1s = get_hit_epos(*hits, i);
        p1e = get_hit_epos_e(*hits, i);
        span_s = MIN(p1s, p1e);
        span_s = MIN(span_s, elen-1);
        span_e = MAX(p1s, p1e);
        span_e =  MIN(span_e, elen-1);
        span_s = ((span_s+span_e)>>1);
        kv_push(uint64_t, idx->pos, (euid<<32)|span_s);
    }
    radix_sort_ho64(idx->pos.a, idx->pos.a + idx->pos.n);
    idx->n = idx->m = ug->u.n;
    CALLOC(idx->a, idx->n);
    for (k = 1, l = 0; k <= idx->pos.n; ++k) 
    {
        if (k == idx->pos.n || ((idx->pos.a[k]>>32) != (idx->pos.a[l]>>32))) 
        {
            idx->a[idx->pos.a[l]>>32] = (uint64_t)l << 32 | (k - l);
            l = k;
        }
    }
    return idx;
} 

uint64_t get_vw_hits_num(h_w_t *e, uint64_t x, uint64_t y)
{
    uint64_t i, occ;
    hw_aux_t *ep = NULL;
    for (i = occ = 0; i < e->n; i++)
    {
        ep = &(e->a[i]);
        if(((ep->e>>33) == x && (((uint32_t)(ep->e))>>1) == y) || 
           ((ep->e>>33) == y && (((uint32_t)(ep->e))>>1) == x))
        {
            occ++;
        }
    }
    return occ;
}

void update_h_w(h_w_t *e, dens_idx_t *idx, double *max_div)
{
    uint64_t k, l, ii, pi, pos, ori, uid, *id = NULL, idn;
    if(max_div) (*max_div) = 0;
    radix_sort_hw_e(e->a, e->a+e->n);
    for (k = 1, l = 0; k <= e->n; ++k) 
    {   
        if (k == e->n || (e->a[k].e>>32) != (e->a[l].e>>32)) ///same uid
        {
            if(k - l > 1) radix_sort_hw_d(e->a+l, e->a+k);
            ori = (e->a[l].e>>32)&1;
            uid = e->a[l].e>>33;
            id = idx->pos.a + (idx->a[uid]>>32);
            idn = (uint32_t)(idx->a[uid]);
            ii = 0;
            for (pi = l; pi < k; pi++)
            {
                pos = e->a[pi].d>>32;///
                while (ii < idn)
                {
                    if(((uint32_t)id[ii]) == pos)
                    {
                        if(ori)
                        {
                            break;
                        }
                        else
                        {
                            while (ii < idn && (((uint32_t)id[ii]) == pos))
                            {
                                ii++;
                            }
                            ii--;
                            break;
                        } 
                    }
                    ii++;
                }
                if(ii >= idn) fprintf(stderr, "ERROR-1\n");
                e->a[pi].w += (ori? idn-ii: ii+1);
                if(max_div) (*max_div) = MAX((*max_div), e->a[pi].w);
            }
            l = k;
        }
    }

    radix_sort_hw_ew(e->a, e->a+e->n);
    for (k = 1, l = 0; k <= e->n; ++k) 
    {   
        if (k == e->n || ((uint32_t)(e->a[k].e)) != ((uint32_t)(e->a[l].e))) ///same uid
        {
            if(k - l > 1) radix_sort_hw_dw(e->a+l, e->a+k);
            ori = e->a[l].e&1;
            uid = (((uint32_t)e->a[l].e)>>1);
            id = idx->pos.a + (idx->a[uid]>>32);
            idn = (uint32_t)(idx->a[uid]);
            ii = 0;
            for (pi = l; pi < k; pi++)
            {
                pos = (uint32_t)(e->a[pi].d);///
                while (ii < idn)
                {
                    if(((uint32_t)id[ii]) == pos)
                    {
                        if(ori)
                        {
                            break;
                        }
                        else
                        {
                            while (ii < idn && (((uint32_t)id[ii]) == pos))
                            {
                                ii++;
                            }
                            ii--;
                            break;
                        } 
                    }
                    ii++;
                }
                if(ii >= idn) fprintf(stderr, "ERROR-2\n");
                e->a[pi].w += (ori? idn-ii: ii+1);
                if(max_div) (*max_div) = MAX((*max_div), e->a[pi].w);
            }
            l = k;
        }
    }
    // radix_sort_hw_e(e->a, e->a+e->n);
    if(max_div) (*max_div) *= 2;///different with slsa2
}

void print_specfic_hic_hits(kvec_pe_hit *hits, uint64_t v, uint64_t w, ma_ug_t *ug)
{
    uint64_t i, suid, euid, slen, elen, p0s, p0e, p1s, p1e, span_s, span_e, sd, ed;
    for (i = 0; i < hits->a.n; i++)
    {
        if(!hits->a.a[i].id) continue;
        suid = get_hit_suid(*hits, i);
        euid = get_hit_euid(*hits, i);
        if((suid == v && euid == w) || (suid == w && euid == v))
        {
            slen = ug->u.a[suid].len;
            elen = ug->u.a[euid].len;

            p0s = get_hit_spos(*hits, i);
            p0e = get_hit_spos_e(*hits, i);
            span_s = MIN(p0s, p0e);
            span_s = MIN(span_s, slen-1);
            span_e = MAX(p0s, p0e);
            span_e =  MIN(span_e, slen-1);
            span_s = ((span_s+span_e)>>1);
            sd = span_s;


            p1s = get_hit_epos(*hits, i);
            p1e = get_hit_epos_e(*hits, i);
            span_s = MIN(p1s, p1e);
            span_s = MIN(span_s, elen-1);
            span_e = MAX(p1s, p1e);
            span_e =  MIN(span_e, elen-1);
            span_s = ((span_s+span_e)>>1);
            ed = span_s;

            fprintf(stderr, "u-stg%.6lul\t%lu\tv-stg%.6lul\t%lu\n", suid+1, sd, euid + 1, ed);
        }
    }
}

kv_u_trans_t *get_update_trans_idx(horder_t *h, trans_col_t *t_idx)
{
    uint64_t i, k, l;
    uint32_t uid/**, q_n, t_n**/;
    ma_ug_t *ug = h->ug;
    asg_t *rg = h->r_g;
    kv_u_trans_t *idx = NULL; CALLOC(idx, 1);
    u_trans_hit_idx iter;
    u_trans_hit_t hit;
    kvec_t(u_hit_t) b; kv_init(b);
    uint64_t *b_idx = NULL; CALLOC(b_idx, t_idx->ref->idx.n);
    /**u_hit_t *q, *t;**/
    u_hit_t *p = NULL;
    for (i = 0; i < ug->u.n; i++)
    {
        uid = i<<1;
        reset_u_trans_hit_idx(&iter, &uid, 1, ug, rg, t_idx->idx, 0, ug->u.a[i].len); 
        while(get_u_trans_hit(&iter, &hit))///get [qScur, qEcur), [qSpre, qEpre)
        {
            kv_pushp(u_hit_t, b, &p);
            p->Scur = hit.qScur;
            p->Ecur = hit.qEcur;
            p->Spre = hit.qSpre;
            p->Epre = hit.qEpre;
            p->uPre = hit.qn;
            p->uCur = uid;
        }
    }

    radix_sort_u_hit(b.a, b.a + b.n);
    for (k = 1, l = 0; k <= b.n; ++k) 
    {   
        if (k == b.n || (b.a[k].uPre>>1) != (b.a[l].uPre>>1))
        {            
            b_idx[b.a[l].uPre>>1] = (uint64_t)l << 32 | (k - l);
            l = k;
        }
    }

    /**
    for (i = 0; i < t_idx->ref->idx.n; i++)
    {
        q_n = (uint32_t)b_idx[i];
        q = b.a + (b_idx[i]>>32);
        if(!q_n) continue;
    }
    **/
    
    // for (i = 0; i < b.n; i++)
    // {
        
    // }
    
    
    free(b_idx);
    kv_destroy(b);
    return idx;
}

void update_scg(horder_t *h, trans_col_t *t_idx)
{
    uint64_t i, k, l, p0s, p0e, p1s, p1e, span_s, span_e, suid, euid, v, w, slen, elen, sd, ed;
    uint64_t t_hits = 0, a_hits = 0;
    double max_div, we;
    h_w_t e; kv_init(e);
    dens_idx_t *idx = NULL;
    hw_aux_t *ep = NULL;
    ma_ug_t *ug = h->ug;
    kvec_pe_hit *hits = &(h->u_hits);
    osg_arc_t *p = NULL;
    osg_destroy(h->sg.g); 
    h->sg.g = osg_init();
    for (i = 0; i < ug->u.n; i++)
    {
        osg_seq_set(h->sg.g, i, 0);
        h->sg.g->seq[i].mw[0] = h->sg.g->seq[i].mw[1] = 0;
        h->sg.g->seq[i].ez[0] = ug->u.a[i].len>>1;
        h->sg.g->seq[i].ez[1] = ug->u.a[i].len - (ug->u.a[i].len>>1);
    }
    idx = build_interval_idx(hits, ug);///idx is used to get density
    
    for (i = 0, e.n = 0; i < hits->a.n; i++)
    {
        if(!hits->a.a[i].id) continue;//hom hits
        suid = get_hit_suid(*hits, i);
        euid = get_hit_euid(*hits, i);
        if(suid == euid) continue;
        slen = ug->u.a[suid].len;
        elen = ug->u.a[euid].len;

        p0s = get_hit_spos(*hits, i);
        p0e = get_hit_spos_e(*hits, i);
        span_s = MIN(p0s, p0e);
        span_s = MIN(span_s, slen-1);
        span_e = MAX(p0s, p0e);
        span_e =  MIN(span_e, slen-1);
        span_s = ((span_s+span_e)>>1);
        sd = span_s;
        v = suid << 1;
        if(span_s > (slen>>1)) v++;


        p1s = get_hit_epos(*hits, i);
        p1e = get_hit_epos_e(*hits, i);
        span_s = MIN(p1s, p1e);
        span_s = MIN(span_s, elen-1);
        span_e = MAX(p1s, p1e);
        span_e =  MIN(span_e, elen-1);
        span_s = ((span_s+span_e)>>1);
        ed = span_s;
        w = euid << 1;
        if(span_s > (elen>>1)) w++;
        t_hits++;
        
        kv_pushp(hw_aux_t, e, &ep);
        ep->w = 0;
        ep->e = (v<<32)|w; 
        ep->d = (sd<<32)|ed;

        if(v > w)
        {
            ep->e = (w<<32)|v;
            ep->d = (ed<<32)|sd;
        } 
        
        // div = h->sg.g->seq[v>>1].ez[v&1] + h->sg.g->seq[w>>1].ez[w&1];
        // max_div = MAX(max_div, div);
        a_hits++;
    }

    // fprintf(stderr, "sa-0-sa: occ-%lu\n", get_vw_hits_num(&e, 12, 690));
    // print_specfic_hic_hits(hits, 12, 690, ug);

    
    
    update_h_w(&e, idx, &max_div);

    // fprintf(stderr, "sa-1-sa: occ-%lu\n", get_vw_hits_num(&e, 676, 738));

    radix_sort_hw_e(e.a, e.a+e.n);
    for (k = 1, l = 0; k <= e.n; ++k) 
    {   
        if (k == e.n || e.a[k].e != e.a[l].e) 
        {
            for (i = 0; i < h->avoid.n; i++)
            {
                v = e.a[l].e;
                w = e.a[l].e<<32; w |= (e.a[l].e>>32);
                if(h->avoid.a[i] == v || h->avoid.a[i] == w)
                {
                    break;
                }
            }

            if(i >= h->avoid.n)
            {
                for (i = l, we = 0; i < k; i++)
                {
                    if(e.a[i].w == 0) fprintf(stderr, "ERROR-3\n");
                    we += (max_div/e.a[i].w);
                }
                
                p = osg_arc_pushp(h->sg.g);
                p->u = p->v = p->occ = p->del = p->w = p->nw = 0;
                p->u = e.a[l].e>>32; p->v = (uint32_t)e.a[l].e;
                p->occ = k - l; p->w = we;
                // if(div != 0) p->w = (double)(k - l)*(max_div/div);
                p = osg_arc_pushp(h->sg.g);
                p->u = p->v = p->occ = p->del = p->w = p->nw = 0;
                p->u = (uint32_t)e.a[l].e; p->v = e.a[l].e>>32;
                p->occ = k - l; p->w = we;
                

                h->sg.g->seq[e.a[l].e>>33].mw[(e.a[l].e>>32)&1] 
                                = MAX(h->sg.g->seq[e.a[l].e>>33].mw[(e.a[l].e>>32)&1], p->w);
                h->sg.g->seq[((uint32_t)e.a[l].e)>>1].mw[e.a[l].e&1]
                                = MAX(h->sg.g->seq[((uint32_t)e.a[l].e)>>1].mw[e.a[l].e&1], p->w);
            }

            l = k;
        }
    }

    // fprintf(stderr, "sa-2-sa: occ-%lu\n", get_vw_hits_num(&e, 676, 738));
    osg_cleanup(h->sg.g);
    double bestAlt;
    uint64_t eg_edges = 0;
    for (i = 0; i < h->sg.g->n_arc; i++)///all p->w should be >= 2
    {
        p = &(h->sg.g->arc[i]);
        bestAlt = MAX(h->sg.g->seq[p->u>>1].mw[p->u&1], h->sg.g->seq[p->v>>1].mw[p->v&1]);
        if(p->w >= bestAlt*0.95)///acutally should be p->w == bestAlt
        {
            bestAlt = get_max_weight(p->u, p->v, h->sg.g);
        }

        if(bestAlt == 0) bestAlt = 1;
        p->nw = p->w/bestAlt;
        if(p->nw > 1) eg_edges++;
    }

    fprintf(stderr, "[M::%s::] # Nodes: %u, # Edges: %u, # Best Edges: %lu, t_hits: %lu, a_hits: %lu\n", 
    __func__, h->sg.g->n_seq, h->sg.g->n_arc, eg_edges, t_hits, a_hits);

    /*******************************for debug************************************/
    /**
    for (i = 0; i < h->sg.g->n_arc; i++)
    {
        p = &(h->sg.g->arc[i]);
        fprintf(stderr, "u-stg%.6ul(%c)(div:%f)\tv-stg%.6ul(%c)(div:%f)\tocc:%u\tw:%f\tnw:%f\n", 
        (p->u>>1)+1, "+-"[p->u&1], h->sg.g->seq[p->u>>1].ez[p->u&1],
        (p->v>>1)+1, "+-"[p->v&1], h->sg.g->seq[p->v>>1].ez[p->v&1], p->occ, p->w, p->nw);
    }
    fprintf(stderr, "sbsbsbsb\n\n\n\n\n\n");
    **/
    
    // uint32_t u, nv, f;
    // osg_arc_t *av = NULL;
    // for (k = 0; k < h->sg.g->n_arc; k++)
    // {
    //     p = &(h->sg.g->arc[k]);
    //     u = p->u; v = p->v;
    //     f = 0;

    //     nv = asg_arc_n(h->sg.g, u);
    //     av = asg_arc_a(h->sg.g, u);
    //     for (i = 0; i < nv; i++)
    //     {
    //         if(av[i].v == v) continue;
    //         if(av[i].w > p->w) f = 1;
    //     }

    //     nv = asg_arc_n(h->sg.g, v);
    //     av = asg_arc_a(h->sg.g, v);
    //     for (i = 0; i < nv; i++)
    //     {
    //         if(av[i].v == u) continue;
    //         if(av[i].w > p->w) f = 1;
    //     }

    //     if(p->nw > 1 && f == 1) fprintf(stderr, "ERROR1\n");
    //     if(p->nw <= 1 && f == 0)
    //     {
    //         fprintf(stderr, "\nERROR2, nw-%f, w-%f, u-%u, v-%u\n", p->nw, p->w, p->u, p->v);
    //          nv = asg_arc_n(h->sg.g, u);
    //         av = asg_arc_a(h->sg.g, u);
    //         for (i = 0; i < nv; i++)
    //         {
    //             if(av[i].v == v) continue;
    //             fprintf(stderr, "+u-%u, v-%u, w-%f\n", av[i].u, av[i].v, av[i].w);
    //         }

    //         nv = asg_arc_n(h->sg.g, v);
    //         av = asg_arc_a(h->sg.g, v);
    //         for (i = 0; i < nv; i++)
    //         {
    //             if(av[i].v == u) continue;
    //             fprintf(stderr, "-u-%u, v-%u, w-%f\n", av[i].u, av[i].v, av[i].w);
    //         }
    //     } 
        
    // }
    /*******************************for debug************************************/
    kv_destroy(e); 
    free(idx->pos.a);
    free(idx->a);
    free(idx);
}

int cmp_arc_nw(const void * a, const void * b)
{
    if((*(osg_arc_t*)a).nw == (*(osg_arc_t*)b).nw) return 0;
    return (*(osg_arc_t*)a).nw < (*(osg_arc_t*)b).nw ? 1 : -1;
}

#define arc_first(g, v) ((g)->arc[(g)->idx[(v)]>>32])
void get_backbone_layout(horder_t *h, sc_lay_t *sl, osg_t *lg, uint8_t *vis)
{
    uint32_t k, v, nc = 0, c = 0;
    lay_t *p = NULL;
    osg_arc_t *t = NULL;
    sl->n = 0;
    ///in lg, there might be single-path paths or cycles
    memset(vis, 0, sizeof(uint8_t)*(lg->n_seq<<1));
    for (k = 0; k < lg->n_seq; k++)
    {
        if((asg_arc_n(lg, k<<1))^(asg_arc_n(lg, (k<<1)+1)))
        {
            v = (asg_arc_n(lg, k<<1)?(k<<1):((k<<1)+1));
            if(vis[k<<1] || vis[(k<<1)+1]) continue;
            kv_pushp(lay_t, *sl, &p);
            kv_init(*p);
            kv_push(uint32_t, *p, v^1);
            kv_push(uint32_t, *p, v);
            vis[v] = vis[v^1] = 1;

            while (asg_arc_n(lg, v))
            {
                v = (arc_first(lg, v).v)^1;
                kv_push(uint32_t, *p, v^1);
                kv_push(uint32_t, *p, v);
                vis[v] = vis[v^1] = 1;
            }
        }
    }
    
    nc = sl->n;
    for (k = 0; k < lg->n_seq; k++)
    {
        if(vis[k<<1] || vis[(k<<1)+1]) continue;
        if(asg_arc_n(lg, k<<1) && asg_arc_n(lg, (k<<1)+1))//circle
        {
            v = k<<1; t = NULL;
            while (asg_arc_n(lg, v))
            {
                if((!t) || (t->nw > arc_first(lg, v).nw) 
                            || (t->nw == arc_first(lg, v).nw && t->w > arc_first(lg, v).w))
                {
                    t = &(arc_first(lg, v));
                } 
                v = (arc_first(lg, v).v)^1;
                if(v == (k<<1)) break;
            }

            v = t->v^1;
            kv_pushp(lay_t, *sl, &p);
            kv_init(*p);
            kv_push(uint32_t, *p, v^1);
            kv_push(uint32_t, *p, v);
            vis[v] = vis[v^1] = 1;

            while (1)
            {
                v = (arc_first(lg, v).v)^1;
                if(vis[v]) break;
                kv_push(uint32_t, *p, v^1);
                kv_push(uint32_t, *p, v);
                vis[v] = vis[v^1] = 1;
            }
        }
    }

    c = sl->n - nc;
    fprintf(stderr, "[M::%s::] # Scaffolds: %u, # non-circles: %u, # circles: %u\n", 
    __func__, (uint32_t)sl->n, nc, c);
    /*******************************for debug************************************/
    // for (k = 0; k < sl->n; k++)
    // {
    //     p = &(sl->a[k]);
    //     if(k >= nc)
    //     {
    //         fprintf(stderr, "%s:\t", k < nc?"non-circle":"circle");
    //         for (i = 0; i < p->n; i+=2)
    //         {
    //             if((p->a[i]>>1) != (p->a[i+1]>>1)) fprintf(stderr, "ERROR-S\n");
    //             fprintf(stderr, "utg%.6ul[%u%u](%u)#", p->a[i]>>1, p->a[i]&1, p->a[i+1]&1, h->ug->u.a[p->a[i]>>1].len);
    //         }
    //         fprintf(stderr, "\n");
    //     }

    //     for (i = 1; i < p->n; i+=2)
    //     {
    //         if(k < nc)
    //         {
    //             if(i < p->n - 1)
    //             {
    //                 if(asg_arc_n(lg, p->a[i])!=1) fprintf(stderr, "ERROR-A\n");
    //                 if(arc_first(lg, p->a[i]).v!=p->a[i+1]) fprintf(stderr, "ERROR-B\n");
    //             }
                
    //             if(i == p->n - 1)
    //             {
    //                 if(asg_arc_n(lg, p->a[i])!=0) fprintf(stderr, "ERROR-A-0\n");
    //             }
    //         }


    //         if(k >= nc)
    //         {
    //             if(i < p->n - 1)
    //             {
    //                 if(asg_arc_n(lg, p->a[i])!=1) fprintf(stderr, "ERROR-A\n");
    //                 if(arc_first(lg, p->a[i]).v!=p->a[i+1]) fprintf(stderr, "ERROR-B\n");
    //                 fprintf(stderr, "i-%u, nw-%f\n", i, arc_first(lg, p->a[i]).nw);
    //             }

    //             if(i == p->n - 1)
    //             {
    //                 if(asg_arc_n(lg, p->a[i])!=1) fprintf(stderr, "ERROR-A\n");
    //                 if(arc_first(lg, p->a[i]).v!=p->a[0]) fprintf(stderr, "ERROR-B-0\n");
    //                 fprintf(stderr, "i-%u, nw-%f\n", i, arc_first(lg, p->a[i]).nw);
    //             }
    //         }
    //     }
        
    // }
    

    // for (k = 0; k < lg->n_seq; k++)
    // {
    //     if(!asg_arc_n(lg, k<<1) && !asg_arc_n(lg, (k<<1)+1))
    //     {
    //         if(vis[k<<1] || vis[(k<<1)+1]) fprintf(stderr, "ERROR-bone\n");
    //     }
    // }
    /*******************************for debug************************************/
}

/**
static void worker_for_insert(void *data, long i, int tid) // callback for kt_for()
{
    sc_id_t *s = &((*(sc_mul*)(data)).a[i]);
    osg_t *sg = (*(sc_mul*)(data)).sg;
    sc_lay_t *sl = (*(sc_mul*)(data)).sl;
    lay_t *p = NULL;
    uint32_t k, i, uid = s->uid;
    for (k = 0; k < sl->n; k++)
    {
        p = &(sl->a[k]);
        for (i = 0; i < p->n; i += 2)
        {
        }
    }
}


void refine_layout(horder_t *h, sc_lay_t *sl, uint8_t *vis)
{
    uint32_t i;
    sc_id_t *p = NULL;
    sc_mul st; kv_init(st);
    st.sl = sl; st.sg = h->sg.g; st.n_thread = asm_opt.thread_num;
    for (i = 0; i < h->sg.g->n_seq; i++)
    {
        if(vis[i<<1]) continue;
        kv_pushp(sc_id_t, st, &p);
        p->uid = i;
        p->iid = p->ori = p->sid = 0;
    }
    
    while (st.n)
    {
        kt_for(st.n_thread, worker_for_insert, &st, st.n);
    }
    sc_id_t

    kv_destroy(st);
}
**/
void refine_layout_back(horder_t *h, sc_lay_t *sl, uint8_t *vis)
{
    uint32_t i, k, m, v, nv, max_k;
    osg_arc_t *av = NULL;
    lay_t *p = NULL;
    uint8_t *sgv = NULL; MALLOC(sgv, sl->n);
    double *w = NULL; MALLOC(w, sl->n);
    uint32_t *idx = NULL; MALLOC(idx, h->sg.g->n_seq);
    memset(idx, -1, sizeof(uint32_t)*h->sg.g->n_seq);
    kvec_t(uint64_t) p_refine; kv_init(p_refine);

    for (k = 0; k < sl->n; k++)
    {
        p = &(sl->a[k]);
        for (m = 0; m < p->n; m++)
        {
            idx[p->a[m]>>1] = k;
        }
    }

    for (i = 0, p_refine.n = 0; i < h->sg.g->n_seq; i++)
    {
        if(vis[i<<1]) continue;
        if(!asg_arc_n(h->sg.g, i<<1)&&!asg_arc_n(h->sg.g, (i<<1)+1)) continue;
        for (k = 0; k < sl->n; k++) w[k] = 0, sgv[k] = 0;

        v = i<<1;
        nv = asg_arc_n(h->sg.g, v);
        av = asg_arc_a(h->sg.g, v);
        for (k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            w[idx[av[k].v>>1]] += av[k].nw;
            sgv[idx[av[k].v>>1]] = 1;
        }
        

        v = (i<<1) + 1;
        nv = asg_arc_n(h->sg.g, v);
        av = asg_arc_a(h->sg.g, v);
        for (k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            w[idx[av[k].v>>1]] += av[k].nw;
            sgv[idx[av[k].v>>1]] = 1;
        }

        for (k = 0, max_k = (uint32_t)-1; k < sl->n; k++)
        {
            if(!sgv[k]) continue;
            if(max_k == (uint32_t)-1 || w[max_k] < w[k]) max_k = k;
        }     

        if(max_k != (uint32_t)-1)
        {
            vis[i<<1] = vis[(i<<1) + 1] = 1;
            kv_push(uint64_t, p_refine, (((uint64_t)(max_k))<<32)|((uint64_t)(i)));
        }    
    }
    free(w); free(idx); free(sgv);

    for (i = 0; i < p_refine.n; i++)
    {
        p = &(sl->a[p_refine.a[i]>>32]);
        // uid = (uint32_t)p_refine.a[i];
    }

    kv_destroy(p_refine);
}

uint32_t get_max_anchor(horder_t *h, sc_lay_t *sl, uint8_t *vis, double *w, uint8_t *sgv, uint32_t *idx,
uint32_t *max_utg, uint32_t *max_sc)
{
    (*max_utg) = (*max_sc) = (uint32_t)-1;
    double max_utg_w = -1;
    uint32_t i, k, v, nv, max_k;
    osg_arc_t *av = NULL;
    for (i = 0; i < h->sg.g->n_seq; i++)
    {
        if(vis[i<<1]) continue;
        if(!asg_arc_n(h->sg.g, i<<1)&&!asg_arc_n(h->sg.g, (i<<1)+1)) continue;
        for (k = 0; k < sl->n; k++) w[k] = 0, sgv[k] = 0;

        v = i<<1;
        nv = asg_arc_n(h->sg.g, v);
        av = asg_arc_a(h->sg.g, v);
        for (k = 0; k < nv; k++)
        {
            if(av[k].del || idx[av[k].v>>1] == (uint32_t)-1) continue;
            w[idx[av[k].v>>1]] += av[k].nw;
            sgv[idx[av[k].v>>1]] = 1;
        }
        

        v = (i<<1) + 1;
        nv = asg_arc_n(h->sg.g, v);
        av = asg_arc_a(h->sg.g, v);
        for (k = 0; k < nv; k++)
        {
            if(av[k].del || idx[av[k].v>>1] == (uint32_t)-1) continue;
            w[idx[av[k].v>>1]] += av[k].nw;
            sgv[idx[av[k].v>>1]] = 1;
        }

        for (k = 0, max_k = (uint32_t)-1; k < sl->n; k++)
        {
            if(!sgv[k]) continue;
            if(max_k == (uint32_t)-1 || w[max_k] < w[k]) max_k = k;
        }     

        if(max_k != (uint32_t)-1)
        {
            if((*max_utg) == (uint32_t)-1 || max_utg_w < w[max_k])
            {
                max_utg_w = w[max_k];
                (*max_utg) = i;
                (*max_sc) = max_k;
            }
        }    
    }

    return (*max_utg) == (uint32_t)-1?0:1;
}

osg_arc_t *get_osg_arc(osg_t *g, uint32_t u, uint32_t v)
{
    osg_arc_t *au = asg_arc_a(g, u);
    uint32_t i, nu = asg_arc_n(g, u);
    for (i = 0; i < nu; i++)
    {
        if(au[i].del) continue;
        if(au[i].v == v) return &(au[i]);
    }
    
    return NULL;
}

void insert_sc(osg_t *g, lay_t *p, uint32_t uid)
{
    uint32_t i, b, e, vb, ve, max_i = (uint32_t)-1, is_found, ori, max_ori = (uint32_t)-1;
    osg_arc_t *bE = NULL, *eE = NULL;
    double w[2], s_w, max_w = -1;
    for (i = 1; i+1 < p->n; i++)///middle points
    {
        b = p->a[i]; e = p->a[i+1]; w[0] = w[1] = 0; 
        s_w = 0; ori = 0; is_found = 0;

        vb = (uid<<1); ve = (uid<<1)+1;
        bE = get_osg_arc(g, vb, b);
        eE = get_osg_arc(g, ve, e);
        if(bE || eE) ///different with slsa2
        {
            if(bE) w[0] += bE->nw;
            if(eE) w[0] += eE->nw;
            is_found++;
        }


        vb = (uid<<1)+1; ve = (uid<<1);
        bE = get_osg_arc(g, vb, b);
        eE = get_osg_arc(g, ve, e);
        if(bE || eE) ///different with slsa2
        {
            if(bE) w[1] += bE->nw;
            if(eE) w[1] += eE->nw;
            is_found++;
        }

        if(is_found > 0)
        {
            s_w = MAX(w[0], w[1]);
            ori = ((w[0] >= w[1])? 0 : 1);
            
            if(max_i == (uint32_t)-1 || max_w < s_w)
            {
                max_w = s_w;
                max_i = i;
                max_ori = ori;
            }
        }
    }


    ///beg point, ///different with slsa2
    i = 0;
    b = (uint32_t)-1; e = p->a[0]; w[0] = w[1] = 0; 
    s_w = 0; ori = 0; is_found = 0;

    vb = (uint32_t)-1; ve = (uid<<1)+1;
    bE = NULL;
    eE = get_osg_arc(g, ve, e);
    if(bE || eE) ///different with slsa2
    {
        if(bE) w[0] += bE->nw;
        if(eE) w[0] += eE->nw;
        is_found++;
    }

    vb = (uint32_t)-1; ve = (uid<<1);
    bE = NULL;
    eE = get_osg_arc(g, ve, e);
    if(bE || eE) ///different with slsa2
    {
        if(bE) w[1] += bE->nw;
        if(eE) w[1] += eE->nw;
        is_found++;
    }

    if(is_found > 0)
    {
        s_w = MAX(w[0], w[1]);
        ori = ((w[0] >= w[1])? 0 : 1);
        
        if(max_i == (uint32_t)-1 || max_w < s_w)
        {
            max_w = s_w;
            max_i = i;
            max_ori = ori;
        }
    }

    ///end point, ///different with slsa2
    i = p->n - 1;
    b = p->a[p->n - 1]; e = (uint32_t)-1; w[0] = w[1] = 0; 
    s_w = 0; ori = 0; is_found = 0;

    vb = (uid<<1); ve = (uint32_t)-1;
    bE = get_osg_arc(g, vb, b);
    eE = NULL;
    if(bE || eE) ///different with slsa2
    {
        if(bE) w[0] += bE->nw;
        if(eE) w[0] += eE->nw;
        is_found++;
    }

    vb = (uid<<1)+1; ve = (uint32_t)-1;
    bE = get_osg_arc(g, vb, b);
    eE = NULL;
    if(bE || eE) ///different with slsa2
    {
        if(bE) w[1] += bE->nw;
        if(eE) w[1] += eE->nw;
        is_found++;
    }

    if(is_found > 0)
    {
        s_w = MAX(w[0], w[1]);
        ori = ((w[0] >= w[1])? 0 : 1);
        
        if(max_i == (uint32_t)-1 || max_w < s_w)
        {
            max_w = s_w;
            max_i = i;
            max_ori = ori;
        }
    }

    if(max_i != 0 && max_i != (uint32_t)-1) max_i++;

    i = p->n;
    kv_resize(uint32_t, *p, p->n+2);
    p->n += 2; 
    while (i > max_i)
    {
        i--;
        p->a[i+2] = p->a[i];
    }
    
    p->a[max_i] = (uid<<1) + max_ori;
    p->a[max_i+1] = (uid<<1) + 1 - max_ori;
}


uint32_t get_vis_occ(uint8_t *vis, uint32_t n)
{
    uint32_t i, occ;
    n <<= 1;
    for (i = occ = 0; i < n; i++)
    {
        if(vis[i]) occ++;
    }
    return occ;
}

uint32_t get_sl_occ(sc_lay_t *sl)
{
    uint32_t i, occ;
    for (i = occ = 0; i < sl->n; i++)
    {
        occ += sl->a[i].n;
    }
    return occ;
}

void refine_layout(horder_t *h, sc_lay_t *sl, uint8_t *vis)
{
    uint32_t k/**m, max_utg, max_sc**/;
    lay_t *p = NULL;
    uint8_t *sgv = NULL; MALLOC(sgv, sl->n);
    double *w = NULL; MALLOC(w, sl->n);
    
    /**
    uint32_t *idx = NULL; MALLOC(idx, h->sg.g->n_seq);
    memset(idx, -1, sizeof(uint32_t)*h->sg.g->n_seq);

    for (k = 0; k < sl->n; k++)
    {
        p = &(sl->a[k]);
        for (m = 0; m < p->n; m++)
        {
            idx[p->a[m]>>1] = k;
        }
    }
    while (get_max_anchor(h, sl, vis, w, sgv, idx, &max_utg, &max_sc))
    {
        
        insert_sc(h->sg.g, &(sl->a[max_sc]), max_utg);

        vis[max_utg<<1] = vis[(max_utg<<1)+1] = 1;
        idx[max_utg] = max_sc;
    }
    free(idx); 
    **/

    for (k = 0; k < h->sg.g->n_seq; k++)
    {
        if(vis[k<<1]) continue;
        kv_pushp(lay_t, *sl, &p);
        kv_init(*p);
        kv_push(uint32_t, *p, (k<<1));
        kv_push(uint32_t, *p, (k<<1)+1);
        vis[(k<<1)] = vis[(k<<1)+1] = 1;
    }

    free(w); free(sgv);
}



void generate_scaffold(ma_utg_t *su, lay_t *ly, ma_ug_t *pug, asg_t *rg)
{
    ma_utg_t *uu = NULL;
    uint32_t i, k, r_i, uid, ori, nv, is_circle = 0;
    uint64_t v, w, l, totalLen;
    asg_arc_t *av  = NULL;
    memset(su, 0, sizeof(*su));
    for (i = 0; i < ly->n; i += 2)
    {
        ori = ly->a[i]&1;
        uid = ly->a[i]>>1;
        uu = &(pug->u.a[uid]);
        is_circle = uu->circ;
        for (r_i = 0; r_i < uu->n; r_i++)
        {   
            v = (ori?uu->a[uu->n - r_i - 1]:uu->a[r_i]);
            if(v != (uint64_t)-1 && ori) v ^= (uint64_t)(0x100000000);
            kv_push(uint64_t, *su, v);
        }
        if(i < ly->n - 2) kv_push(uint64_t, *su, (uint64_t)-1);
    }
    if(ly->n != 2) is_circle = 0;

    for (i = 0, totalLen = 0; i < su->n-1; i++)
    {
        if(su->a[i] == (uint64_t)-1)
        {
            totalLen += GAP_LEN;
            continue;
        }
        v = su->a[i]>>32;
        if(su->a[i+1] == (uint64_t)-1)
        {
            l = rg->seq[v>>1].len;
        }
        else
        {
            w = su->a[i+1]>>32; 
            av = asg_arc_a(rg, v);
            nv = asg_arc_n(rg, v);

            l = 0;
            for (k = 0; k < nv; k++)
            {
                if(av[k].del) continue;
                if(av[k].v == w) 
                {
                    l = asg_arc_len(av[k]);
                    break;
                }
            }
            if(k == nv) fprintf(stderr, "ERROR-scf-0, v-%lu, w-%lu\n", v, w);
        }
        
        su->a[i] = v; su->a[i] = su->a[i]<<32; su->a[i] = su->a[i] | (uint64_t)(l);
        totalLen += l;
    }
    if(i < su->n)
    {
        if(su->a[i] == (uint64_t)-1)
        {
            totalLen += GAP_LEN;
        }
        else
        {
            if(is_circle && su->a[0] != (uint64_t)-1)
            {
                v = su->a[i]>>32;
                w = su->a[0]>>32; 
                av = asg_arc_a(rg, v);
                nv = asg_arc_n(rg, v);

                l = 0;
                for (k = 0; k < nv; k++)
                {
                    if(av[k].del) continue;
                    if(av[k].v == w) 
                    {
                        l = asg_arc_len(av[k]);
                        break;
                    }
                }
                if(k == nv) fprintf(stderr, "ERROR-scf-1, v-%lu, w-%lu\n", v, w);

                su->a[i] = v; su->a[i] = su->a[i]<<32; su->a[i] = su->a[i] | (uint64_t)(l);
                totalLen += l;
            }
            else
            {
                v = su->a[i]>>32;
                l = rg->seq[v>>1].len;
                su->a[i] = v;
                su->a[i] = su->a[i]<<32;
                su->a[i] = su->a[i] | (uint64_t)(l);
                totalLen += l;
            }
        }
    }

    su->circ = is_circle;
    su->len = totalLen;
    if(!su->circ)
    {
        su->start = su->a[0]>>32;
        su->end = (su->a[su->n-1]>>32)^1;
    }
    else
    {
        su->start = su->end = UINT32_MAX;
    }
}

uint64_t get_nuid(sc_lay_t *sl, uint64_t *p)
{
    uint64_t k, ouid[2];
    for (k = 0; k < sl->n; k++)
    {
        ouid[0] = sl->a[k].a[0]; 
        ouid[1] = sl->a[k].a[sl->a[k].n - 1];
        if((*p) == ouid[0])
        {
            (*p) = (k<<1);
            return 1;
        } 

        if((*p) == ouid[1])
        {
            (*p) = (k<<1)+1;
            return 1;
        }
    }

    return 0;
}

void update_avoids(horder_t *h, sc_lay_t *sl)
{
    uint64_t i, m, ps, pe;
    for (i = m = 0; i < h->avoid.n; i++)
    {
        ps = h->avoid.a[i]>>32; 
        pe = (uint32_t)h->avoid.a[i];
        if(get_nuid(sl, &ps) && get_nuid(sl, &pe))
        {
            h->avoid.a[m] = (ps<<32)|pe;
            m++;
        }
    }
    h->avoid.n = m;
}

void update_ug_by_layout(horder_t *h, sc_lay_t *sl, ma_ug_t* i_ug)
{
    uint32_t i;
    lay_t *p = NULL;
    ma_utg_t *pu = NULL;
    ma_ug_t *sug = NULL;
    kvec_t_u64_warp n_avoids; kv_init(n_avoids.a);
    sug = (ma_ug_t*)calloc(1, sizeof(ma_ug_t));
    for (i = 0; i < sl->n; i++)
    {
        p = &(sl->a[i]);
        kv_pushp(ma_utg_t, sug->u, &pu);
        generate_scaffold(pu, p, i_ug?i_ug:h->ug, h->r_g);
    }
    ma_ug_destroy(h->ug);
    h->ug = sug;
    kv_destroy(n_avoids.a);
    update_avoids(h, sl);
}

void get_long_switch_scaffolds(horder_t *h, sc_lay_t *sl, osg_t *lg)
{
    fprintf(stderr, "\n[M::%s::]\n", __func__);
    uint32_t i, k, r_i, ori, sw[3], sw_inner[3];
    uint64_t v, len;
    kvec_t(uint64_t) idx; kv_init(idx);
    lay_t *p = NULL;
    ma_utg_t *u = NULL;
    for (i = 0; i < sl->n; i++)
    {
        p = &(sl->a[i]);
        sw[0] = sw[1] = sw[2] = 0;
        for (k = 0; k < p->n; k += 2)
        {
            u = &(h->ug->u.a[p->a[k]>>1]);
            ori = p->a[k]&1;
            for (r_i = 0; r_i < u->n; r_i++)
            {
                v = (ori?u->a[u->n - r_i - 1]:u->a[r_i]);
                if(v != (uint64_t)-1)
                {
                    v >>= 33;
                    sw[R_INF.trio_flag[v]]++;
                }
            }
        }

        v = MIN(sw[FATHER], sw[MOTHER]);
        v = ((uint32_t)-1) - v;
        v <<= 32; v |= i;
        kv_push(uint64_t, idx, v);
    }

    radix_sort_ho64(idx.a, idx.a+idx.n);
    for (i = 0; i < sl->n; i++)
    {
        p = &(sl->a[(uint32_t)(idx.a[i])]);
        fprintf(stderr, "\nscaf-%u-th, occ-%u\n", (uint32_t)(idx.a[i]), (uint32_t)(p->n>>1));
        sw[0] = sw[1] = sw[2] = len = 0;
        for (k = 0; k < p->n; k += 2)
        {
            sw_inner[0] = sw_inner[1] = sw_inner[2] = 0;
            u = &(h->ug->u.a[p->a[k]>>1]);
            ori = p->a[k]&1;
            for (r_i = 0; r_i < u->n; r_i++)
            {
                v = (ori?u->a[u->n - r_i - 1]:u->a[r_i]);
                if(v != (uint64_t)-1)
                {
                    v >>= 33;
                    if(R_INF.trio_flag[v] == FATHER || R_INF.trio_flag[v] == MOTHER)
                    {
                        sw[R_INF.trio_flag[v]]++;
                        sw_inner[R_INF.trio_flag[v]]++;
                    }
                }
            }
            fprintf(stderr, "utg%.6ul (ori: %u), u->len-%u, sw_in[FATHER]-%u, sw_in[MOTHER]-%u\n", 
            (p->a[k]>>1)+1, p->a[k]&1, u->len, sw_inner[FATHER], sw_inner[MOTHER]);
            len += u->len + ((k + 2)< p->n? GAP_LEN:0);
        }
        fprintf(stderr, "sw[FATHER]-%u, sw[MOTHER]-%u, len-%lu\n", sw[FATHER], sw[MOTHER], len);
    }
    kv_destroy(idx);
}

void destory_sc_lay_t(sc_lay_t *sl)
{
    uint32_t i;
    for (i = 0; i < sl->n; i++) kv_destroy(sl->a[i]);
    kv_destroy(*sl);
}

void layout_scg(horder_t *h, double nw_thres, uint32_t occ_thres, sc_lay_t *r_sl)
{
    uint32_t k;
    osg_arc_t *p = NULL, *lp = NULL;
    uint8_t *vis = NULL; CALLOC(vis, h->sg.g->n_seq<<1);
    sc_lay_t sl; kv_init(sl);
    osg_t *lg = osg_init();
    qsort(h->sg.g->arc, h->sg.g->n_arc, sizeof(osg_arc_t), cmp_arc_nw);
    for (k = 0; k < h->sg.g->n_arc; k++)
    {
        p = &(h->sg.g->arc[k]);
        if(vis[p->u] || vis[p->v]) continue;
        ///different with slsa2
        if(p->nw <= nw_thres || p->occ <= occ_thres) continue;
        vis[p->u] = vis[p->v] = 1;
        lp = osg_arc_pushp(lg);
        (*lp) = (*p);
        lp = osg_arc_pushp(lg);
        (*lp) = (*p);
        lp->u = p->v;
        lp->v = p->u;
    }
    for (k = 0; k < h->sg.g->n_seq; k++)
    {
        osg_seq_set(lg, k, 0);
    }
    osg_cleanup(lg);
    radix_sort_osg(h->sg.g->arc, h->sg.g->arc + h->sg.g->n_arc);

    get_backbone_layout(h, &sl, lg, vis);

    // get_long_switch_scaffolds(h, &sl, lg);
    
    refine_layout(h, &sl, vis);
    
    // print_N50_layout(h->ug, &sl);

    update_ug_by_layout(h, &sl, NULL);

    print_N50(h->ug);

    if(r_sl){
        r_sl->a = sl.a; sl.a = NULL;
        r_sl->m = sl.m; sl.m = 0;
        r_sl->n = sl.n; sl.n = 0;

    }
    destory_sc_lay_t(&sl);
    osg_destroy(lg);
    free(vis);
}

void renew_scaffold(horder_t *h)
{
    double index_time = yak_realtime();
    while (1)
    {
        update_u_hits(&(h->u_hits), &(h->r_hits), h->ug, h->r_g);
        if(!break_scaffold(h, 5, 15, (uint64_t)-1, (uint64_t)-1, 2500000, 1, NULL)) break;
        print_N50(h->ug);
    }
    fprintf(stderr, "[M::%s::%.3f] \n", __func__, yak_realtime()-index_time);
}

void print_scaffold(ma_ug_t *ug, asg_t *sg, ma_sub_t* coverage_cut, char* output_file_name, 
ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, 
long long tipsLen, float tip_drop_ratio, long long stops_threshold, R_to_U* ruIndex, 
float chimeric_rate, float drop_ratio, int max_hang, int min_ovlp)
{
    char* gfa_name = (char*)malloc(strlen(output_file_name)+100);
    sprintf(gfa_name, "%s.%s.p_ctg.gfa", output_file_name, "stg");
    fprintf(stderr, "Writing %s to disk... \n", gfa_name);
    FILE* output_file = NULL;
    output_file = fopen(gfa_name, "w");


    kvec_asg_arc_t_warp new_rtg_edges;
    kv_init(new_rtg_edges.a);
    ma_ug_seq_scaffold(ug, sg, coverage_cut, sources, &new_rtg_edges, max_hang, min_ovlp, 0, 1);
    ma_ug_print(ug, sg, coverage_cut, sources, ruIndex, "stg", output_file);
    fclose(output_file);

    sprintf(gfa_name, "%s.%s.p_ctg.noseq.gfa", output_file_name, "stg");
    output_file = fopen(gfa_name, "w");
    ma_ug_print_simple(ug, sg, coverage_cut, sources, ruIndex, "stg", output_file);
    fclose(output_file);

    free(gfa_name);
    kv_destroy(new_rtg_edges.a);
}

void scaffold_hap(horder_t *h, ug_opt_t *opt, trans_col_t *t_idx, uint32_t round, char *output_file_name, uint8_t flag)
{
    uint32_t i;
    kv_destroy(h->u_hits.a);
    kv_destroy(h->u_hits.idx);
    kv_destroy(h->u_hits.occ);
    memset(&(h->u_hits), 0, sizeof(h->u_hits));
    kv_destroy(h->avoid);
    h->avoid.m = h->avoid.n = 0;
    h->avoid.a = NULL;
    osg_destroy(h->sg.g);
    h->sg.g = NULL;
    ma_ug_destroy(h->ug);
    h->ug = NULL;

    h->ug = get_trio_unitig_graph(h->r_g, flag, opt);
    asg_destroy(h->ug->g);
    h->ug->g = NULL;
    // kv_push(uint64_t, h->occ, h->ug->u.n);
    // CALLOC(h->hf.a, h->ug->u.n); h->hf.n = h->hf.m = h->ug->u.n;

    print_N50(h->ug);

    for (i = 0; i < round; i++)
    {
        update_u_hits(&(h->u_hits), &(h->r_hits), h->ug, h->r_g);
        update_scg(h, t_idx);
        layout_scg(h, 1.001, 19, NULL);
        renew_scaffold(h);
    }
    
    char* gfa_name = (char*)malloc(strlen(output_file_name)+100);
    sprintf(gfa_name, "%s.%s", output_file_name, (flag==FATHER?"hap1":"hap2"));

    print_scaffold(h->ug, h->r_g, opt->coverage_cut, gfa_name, 
    opt->sources, opt->reverse_sources, opt->tipsLen, opt->tip_drop_ratio, 
    opt->stops_threshold, opt->ruIndex, opt->chimeric_rate, opt->drop_ratio, 
    opt->max_hang, opt->min_ovlp);

    free(gfa_name);
}

void scaffold_ug(horder_t *h, ma_ug_t *ug, ug_opt_t *opt, uint32_t round, char *output_file_name, uint8_t flag)
{
    uint32_t i;
    kv_destroy(h->u_hits.a);
    kv_destroy(h->u_hits.idx);
    kv_destroy(h->u_hits.occ);
    memset(&(h->u_hits), 0, sizeof(h->u_hits));
    kv_destroy(h->avoid);
    h->avoid.m = h->avoid.n = 0;
    h->avoid.a = NULL;
    osg_destroy(h->sg.g);
    h->sg.g = NULL;

    h->ug = ug;
    asg_destroy(h->ug->g);
    h->ug->g = NULL;
    // kv_push(uint64_t, h->occ, h->ug->u.n);
    // CALLOC(h->hf.a, h->ug->u.n); h->hf.n = h->hf.m = h->ug->u.n;

    print_N50(h->ug);
    fprintf(stderr, "[M::%s::]***0***\n", __func__);
    for (i = 0; i < round; i++)
    {
        fprintf(stderr, "[M::%s::]*i->%u*\n", __func__, i);
        update_u_hits(&(h->u_hits), &(h->r_hits), h->ug, h->r_g);
        fprintf(stderr, "[M::%s::]**i->%u**\n", __func__, i);
        update_scg(h, NULL);
        fprintf(stderr, "[M::%s::]***i->%u***\n", __func__, i);
        layout_scg(h, 1.001, 19, NULL);
        fprintf(stderr, "[M::%s::]****i->%u****\n", __func__, i);
        renew_scaffold(h);
        fprintf(stderr, "[M::%s::]*****i->%u*****\n", __func__, i);
    }
    
    char* gfa_name = (char*)malloc(strlen(output_file_name)+100);
    sprintf(gfa_name, "%s.%s", output_file_name, (flag==FATHER?"hap1":"hap2"));

    print_scaffold(h->ug, h->r_g, opt->coverage_cut, gfa_name, 
    opt->sources, opt->reverse_sources, opt->tipsLen, opt->tip_drop_ratio, 
    opt->stops_threshold, opt->ruIndex, opt->chimeric_rate, opt->drop_ratio, 
    opt->max_hang, opt->min_ovlp);

    free(gfa_name);
}

void output_hic_rtg(ma_ug_t *ug, asg_t *rg, ug_opt_t *opt, char* output_file_name)
{
    char* gfa_name = (char*)malloc(strlen(output_file_name)+50);
    sprintf(gfa_name, "%s.all.noseq.gfa", output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    ma_ug_print_simple(ug, rg, opt->coverage_cut, opt->sources, opt->ruIndex, "utg", output_file);
    fclose(output_file);
    free(gfa_name);
}

horder_t *init_horder_t(kvec_pe_hit *i_hits, uint64_t i_hits_uid_bits, uint64_t i_hits_pos_mode, 
asg_t *i_rg, ma_ug_t* i_ug, bubble_type* bub, kv_u_trans_t *ref, ug_opt_t *opt, uint32_t round)
{
    uint32_t i;
    trans_col_t *t_idx = NULL;
    horder_t *h = NULL; CALLOC(h, 1);
    get_r_hits(i_hits, &(h->r_hits), i_rg, i_ug, bub, i_hits_uid_bits, i_hits_pos_mode);
    h->r_g = copy_read_graph(i_rg);
    horder_clean_sg_by_utg(h->r_g, i_ug);
    t_idx = init_trans_col(i_ug, h->r_g->n_seq, ref);
    // output_hic_rtg(i_ug, h->r_g, opt, asm_opt.output_file_name);

    reduce_hamming_error(h->r_g, opt->sources, opt->coverage_cut, opt->max_hang, opt->min_ovlp, opt->gap_fuzz);
    /**
    scaffold_hap(h, t_idx, opt, round, asm_opt.output_file_name, FATHER);
    scaffold_hap(h, t_idx, opt, round, asm_opt.output_file_name, MOTHER);
    **/
    
    
    generate_haplotypes(h, opt);
    print_N50(h->ug);
    
    // update_u_hits(&(h->u_hits), &(h->r_hits), h->ug, h->r_g);
    // break_contig(h, 10, 20);  
    for (i = 0; i < round; i++)
    {
        update_u_hits(&(h->u_hits), &(h->r_hits), h->ug, h->r_g);
        update_scg(h, t_idx);
        layout_scg(h, 1.001, 19, NULL);
        renew_scaffold(h);
    }
    

    print_scaffold(h->ug, h->r_g, opt->coverage_cut, asm_opt.output_file_name, 
    opt->sources, opt->reverse_sources, opt->tipsLen, opt->tip_drop_ratio, 
    opt->stops_threshold, opt->ruIndex, opt->chimeric_rate, opt->drop_ratio, 
    opt->max_hang, opt->min_ovlp);
    
    destory_trans_col(&t_idx);
    exit(1);
    return h;
}

void cpy_u_hits(kvec_pe_hit *u_hits, kvec_pe_hit *i_hits, uint32_t u_n)
{
    uint64_t i;
    memset(u_hits, 0, sizeof(kvec_pe_hit));
    u_hits->pos_mode = i_hits->pos_mode;
    u_hits->uID_bits = i_hits->uID_bits;
    u_hits->a.n = u_hits->a.m = i_hits->a.n;
    MALLOC(u_hits->a.a, u_hits->a.n);
    memcpy(u_hits->a.a, i_hits->a.a, u_hits->a.n*sizeof(pe_hit));
    for (i = 0; i < u_hits->a.n; i++) u_hits->a.a[i].id = 1;
    idx_hits(u_hits, u_n);
}

void update_sc_lay(sc_lay_t *sl, h_covs *b)
{
    uint64_t k, l, i, pidx, cidx;
    lay_t *s = NULL;
    lay_t *p = NULL;
    for (k = 1, l = 0; k <= b->n; ++k) 
    {   
        if (k == b->n || (b->a[k].s != b->a[l].s)) 
        {
            s = &(sl->a[b->a[l].s]);
            for (i = l, pidx = 0; i < k; i++){
                cidx = b->a[i].dp;
                kv_pushp(lay_t, *sl, &p);
                kv_init(*p);
                p->n = p->m = (cidx - pidx + 1)<<1;
                MALLOC(p->a, p->n);
                memcpy(p->a, s->a + (pidx<<1), sizeof(*(p->a))*p->n);
                pidx = cidx + 1;
            }

            if(pidx >= (s->n>>1)) fprintf(stderr, "ERROR-update\n");
            cidx = (s->n>>1)-1;
            kv_pushp(lay_t, *sl, &p);
            kv_init(*p);
            p->n = p->m = (cidx - pidx + 1)<<1;
            MALLOC(p->a, p->n);
            memcpy(p->a, s->a + (pidx<<1), sizeof(*(p->a))*p->n);
            free(s->a); s->n = s->m = 0;
            l = k;
        }
    }

    for (i = k = 0; i < sl->n; i++){
        if(!sl->a[i].a) continue;
        sl->a[k] = sl->a[i];
        sl->a[i].a = NULL; sl->a[i].n = sl->a[i].m = 0;
        if(sl->a[k].n == 2){
            sl->a[k].a[0] >>= 1; sl->a[k].a[0] <<= 1;
            sl->a[k].a[1] >>= 1; sl->a[k].a[1] <<= 1; sl->a[k].a[1]++;
        }
        k++;
    }
    sl->n = k;
}

void renew_scaffold_utg(horder_t *h, sc_lay_t *sl, ma_ug_t* i_ug)
{
    double index_time = yak_realtime();
    h_covs b; kv_init(b);
    while (1)
    {
        update_u_hits(&(h->u_hits), &(h->r_hits), h->ug, h->r_g);
        if(!break_scaffold(h, 5, 15, 15, 25, 2500000, 1, &b)) break;
        update_sc_lay(sl, &b);
        update_ug_by_layout(h, sl, i_ug);
        print_N50(h->ug);
    }
    fprintf(stderr, "[M::%s::%.3f] \n", __func__, yak_realtime()-index_time);
    kv_destroy(b);
}

spg_t *scf_g(sc_lay_t *sl, ma_ug_t* ug)
{
    spg_t *scg = NULL; CALLOC(scg, 1); scg->ug = ug;    
    uint64_t i, k;
    lay_t *p = NULL;
    for (i = 0; i < sl->n; i++){
        p = &(sl->a[i]);
        kv_push(uint64_t, scg->idx, (uint64_t)scg->dst.n << 32 | (p->n>>1));
        for (k = 0; k < p->n; k+=2) kv_push(uint32_t, scg->dst, p->a[k]);        
    }
    return scg;
}

spg_t *horder_utg(kvec_pe_hit *i_hits, uint64_t i_hits_uid_bits, uint64_t i_hits_pos_mode, 
asg_t *i_rg, ma_ug_t* i_ug, bubble_type* bub, ug_opt_t *opt)
{
    horder_t *h = NULL; CALLOC(h, 1);
    sc_lay_t sl; kv_init(sl);
    get_r_hits(i_hits, &(h->r_hits), i_rg, i_ug, bub, i_hits_uid_bits, i_hits_pos_mode);
    h->r_g = copy_read_graph(i_rg);
    horder_clean_sg_by_utg(h->r_g, i_ug);
    h->ug = copy_untig_graph(i_ug); asg_destroy(h->ug->g); h->ug->g = NULL;
    cpy_u_hits(&(h->u_hits), i_hits, h->ug->u.n);

    update_scg(h, NULL);
    layout_scg(h, ((double)1)/((double)0.75), 19, &sl);
    renew_scaffold_utg(h, &sl, i_ug);
    
    spg_t *scg = scf_g(&sl, i_ug);
    destory_sc_lay_t(&sl);
    destory_horder_t(&h);
    return scg;
}


void ha_aware_order(kvec_pe_hit *r_hits, asg_t *rg, ma_ug_t *ug_fa, ma_ug_t *ug_mo, kv_u_trans_t *ref, 
ug_opt_t *opt, uint32_t round)
{
    horder_t *h = NULL; CALLOC(h, 1);
    h->r_hits = *r_hits;
    h->r_g = rg;
    scaffold_ug(h, ug_fa, opt, round, asm_opt.output_file_name, FATHER);
    scaffold_ug(h, ug_mo, opt, round, asm_opt.output_file_name, MOTHER);
}

void destory_horder_t(horder_t **h)
{
    kv_destroy((*h)->r_hits.a);
    kv_destroy((*h)->r_hits.idx);
    kv_destroy((*h)->r_hits.occ);

    kv_destroy((*h)->u_hits.a);
    kv_destroy((*h)->u_hits.idx);
    kv_destroy((*h)->u_hits.occ);

    kv_destroy((*h)->avoid);

    osg_destroy((*h)->sg.g);

    ma_ug_destroy((*h)->ug);
    asg_destroy((*h)->r_g);
    free((*h));
}

kvec_pe_hit *get_r_hits_order(kvec_pe_hit *uhits, uint64_t hits_uid_bits, uint64_t hits_pos_mode, 
asg_t *rg, ma_ug_t* ug, bubble_type* bub)
{
    kvec_pe_hit *r_hits = NULL; CALLOC(r_hits, 1);
    get_r_hits(uhits, r_hits, rg, ug, bub, hits_uid_bits, hits_pos_mode);
    return r_hits;
}