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
KSEQ_INIT(gzFile, gzread)
KDQ_INIT(uint64_t)
#define pe_hit_an1_idx_key(x) ((x).s<<1)
KRADIX_SORT_INIT(pe_hit_idx_hn1, pe_hit, pe_hit_an1_idx_key, member_size(pe_hit, s))
#define pe_hit_an2_idx_key(x) ((x).e<<1)
KRADIX_SORT_INIT(pe_hit_idx_hn2, pe_hit, pe_hit_an2_idx_key, member_size(pe_hit, e))
#define generic_key(x) (x)
KRADIX_SORT_INIT(ho64, uint64_t, generic_key, 8)
#define get_hit_srev(x, k) ((x).a.a[(k)].s>>63)
#define get_hit_slen(x, k) ((x).a.a[(k)].len>>32)
#define get_hit_suid(x, k) (((x).a.a[(k)].s<<1)>>(64 - (x).uID_bits))
#define get_hit_spos(x, k) ((x).a.a[(k)].s & (x).pos_mode)
#define get_hit_spos_e(x, k) (get_hit_srev((x),(k))?\
        ((get_hit_spos((x),(k))+1>=get_hit_slen((x),(k)))?\
            (get_hit_spos((x),(k))+1-get_hit_slen((x),(k))):0)\
                    :(get_hit_spos((x),(k))+get_hit_slen((x),(k))-1))

#define get_hit_erev(x, k) ((x).a.a[(k)].e>>63)
#define get_hit_elen(x, k) ((uint32_t)((x).a.a[(k)].len))
#define get_hit_euid(x, k) (((x).a.a[(k)].e<<1)>>(64 - (x).uID_bits))
#define get_hit_epos(x, k) ((x).a.a[(k)].e & (x).pos_mode)
#define get_hit_epos_e(x, k) (get_hit_erev((x),(k))?\
        ((get_hit_epos((x),(k))+1>=get_hit_elen((x),(k)))?\
            (get_hit_epos((x),(k))+1-get_hit_elen((x),(k))):0)\
                    :(get_hit_epos((x),(k))+get_hit_elen((x),(k))-1))
#define OVL(s_0, e_0, s_1, e_1) ((MIN((e_0), (e_1)) > MAX((s_0), (s_1)))? MIN((e_0), (e_1)) - MAX((s_0), (s_1)):0)
#define BREAK_THRES 1000000
#define BREAK_CUTOFF 0.1
#define BREAK_BOUNDARY 0.01
#define GAP_LEN 100

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
	uint64_t s, e, dp;
} h_cov_t;

typedef struct {
	h_cov_t *a;
    size_t n, m;
} h_covs;

#define h_cov_s_key(x) ((x).s)
KRADIX_SORT_INIT(h_cov_s, h_cov_t, h_cov_s_key, member_size(h_cov_t, s))
#define hit_aux_ruid_key(x) ((x).ruid)
KRADIX_SORT_INIT(hit_aux_ruid, hit_aux_t, hit_aux_ruid_key, member_size(hit_aux_t, ruid))

void print_N50(ma_ug_t* ug)
{
    kvec_t(uint64_t) b; kv_init(b);
    uint64_t i, s, len;
    for (i = s = 0; i < ug->g->n_seq; ++i) 
    {
        kv_push(uint64_t, b, ug->u.a[i].len);
        s += ug->u.a[i].len;
    }
    len = s;
    fprintf(stderr, "[M::%s::] Genome Size: %lu, # Contigs: %u\n", __func__, len, ug->g->n_seq);

    radix_sort_ho64(b.a, b.a+b.n);
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
    x.idx.n = x.idx.m = (x.n?(x.a[x.n-1].ruid>>33)+1:0);
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
    for (u_hits->uID_bits=1; (uint64_t)(1<<u_hits->uID_bits)<(uint64_t)ug->g->n_seq; u_hits->uID_bits++);
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
                t->id = r_hits->a.a[i].id;
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
    if(tmp->n == 1) return;

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

void break_utg_horder(horder_t *h, h_covs *b_points)
{
    uint64_t i, rid, uid;
    for (i = 0; i < b_points->n; i++)
    {
        rid = b_points->a[i].s;
        uid = b_points->a[i].e;
    }
    
}

void break_contig_init(horder_t *h, uint64_t cutoff_s, uint64_t cutoff_e)
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

                fprintf(stderr, "\n[M::%s::] utg%.6lul, ulen: %lu, # hic hits: %lu, map cov: %lu, utg cov: %lu, average: %lu\n", 
                __func__, get_hit_suid(*hits, l)+1, ulen, (uint64_t)(b.n>>1), cov_hic, cov_utg, cov_ava);


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
                        fprintf(stderr, "cutoff: %lu, bs: %lu, be: %lu\n", cutoff, bs, be);
                    }
                }

                if(res.n > 0)
                {
                    get_consensus_break(&res, &cov_buf);
                    for (i = 0; i < cov_buf.n; i++)
                    {
                        fprintf(stderr, "consensus_break-s: %lu, e: %lu\n", cov_buf.a[i].s, cov_buf.a[i].e);
                    }
                    get_read_breaks(&(ug->u.a[get_hit_suid(*hits, l)]), h->r_g, &cov_buf, 
                    &res, hits, l, k, ulen, &bs, &dp);
                    if(bs == (uint64_t)-1) fprintf(stderr, "ERROR-read\n");
                    
                    kv_pushp(h_cov_t, b_points, &p);
                    p->s = bs; p->e = get_hit_suid(*hits, l); p->dp = dp;
                    fprintf(stderr, "consensus_break-rid: %lu, cov: %lu\n", bs, dp);
                    // debug_sub_cov(hits, l, k, ulen, &(ug->u.a[get_hit_suid(*hits, l)]), h->r_g, bs, dp);
                }
                
            }  
            l = k;
        }
    }

    

    kv_destroy(b);
    kv_destroy(cov_buf);
    kv_destroy(res);
    kv_destroy(b_points);
}


void generate_haplotypes(horder_t *h, ug_opt_t *opt)
{
    uint64_t i, off;
    ma_ug_t *ug_1 = NULL, *ug_2 = NULL;
    // asg_cleanup(h->r_g);
    // asg_arc_del_trans(h->r_g, asm_opt.gap_fuzz);///must

    ug_1 = get_trio_unitig_graph(h->r_g, FATHER, opt);
    ug_2 = get_trio_unitig_graph(h->r_g, MOTHER, opt);
    h->ug = ug_1; 

    off = ug_1->g->n_seq;
    ///update graph
    for (i = 0; i < ug_2->g->n_seq; i++)
    {
        asg_seq_set(h->ug->g, off + i, ug_2->g->seq[i].len, ug_2->g->seq[i].del);
    }
    
    asg_arc_t *p = NULL;
    for (i = 0; i < ug_2->g->n_arc; i++)
    {
        p = asg_arc_pushp(h->ug->g);
        *p = ug_2->g->arc[i];
        p->v += (off<<1);
        p->ul += (off<<33);
    }
    free(h->ug->g->idx);
    h->ug->g->idx = 0;
    h->ug->g->is_srt = 0;
    asg_cleanup(h->ug->g);
    


    ///update unitigs
    ma_utg_t *pu = NULL;
    for (i = 0; i < ug_2->u.n; i++)
    {
        kv_pushp(ma_utg_t, h->ug->u, &pu);
        *pu = ug_2->u.a[i];
        ug_2->u.a[i].a = NULL;
        ug_2->u.a[i].s = NULL;
    }

    for (i = 0; i < h->ug->g->n_seq; i++)
    {
        h->ug->g->seq[i].c = (i < off? FATHER:MOTHER);
    }
    /*******************************for debug************************************/
    // if(h->ug->g->n_seq != h->ug->u.n)
    // {
    //     fprintf(stderr, "ERROR-non-equal-length\n");
    // }
    // fprintf(stderr, "h->ug->u.n-%u, off-%lu, ug_2->u.n-%u\n", (uint32_t)h->ug->u.n, off, (uint32_t)ug_2->u.n);
    // for (i = 0; i < h->ug->g->n_seq; i++)
    // {
    //     if(h->ug->g->seq[i].len != h->ug->u.a[i].len)
    //     {
    //         fprintf(stderr, "****ERROR-1-%lu: g->seq[i].len-%u, u.a[i].len-%u\n", 
    //         i, (uint32_t)h->ug->g->seq[i].len, (uint32_t)h->ug->u.a[i].len);
    //     }
    //     pu = i < off? (&ug_1->u.a[i]) : (&ug_2->u.a[i - off]);
    //     if(h->ug->g->seq[i].len != pu->len)
    //     {
    //         fprintf(stderr, "ERROR-2\n");
    //     } 
    // }
    
    /*******************************for debug************************************/
    ug_1 = NULL;
    ma_ug_destroy(ug_2);
    asg_destroy(h->ug->g);
    h->ug->g = NULL;
}

horder_t *init_horder_t(kvec_pe_hit *i_hits, uint64_t i_hits_uid_bits, uint64_t i_hits_pos_mode, 
asg_t *i_rg, ma_ug_t* i_ug, bubble_type* bub, ug_opt_t *opt)
{
    horder_t *h = NULL; CALLOC(h, 1);
    get_r_hits(i_hits, &(h->r_hits), i_rg, i_ug, bub, i_hits_uid_bits, i_hits_pos_mode);
    h->r_g = copy_read_graph(i_rg);
    horder_clean_sg_by_utg(h->r_g, i_ug);
    generate_haplotypes(h, opt);
    print_N50(h->ug);
    update_u_hits(&(h->u_hits), &(h->r_hits), h->ug, h->r_g);  
    break_contig_init(h, 10, 20);  

    return h;
}

void destory_horder_t(horder_t **h)
{
    kv_destroy((*h)->r_hits.a);
    kv_destroy((*h)->r_hits.idx);
    kv_destroy((*h)->r_hits.occ);

    kv_destroy((*h)->u_hits.a);
    kv_destroy((*h)->u_hits.idx);
    kv_destroy((*h)->u_hits.occ);

    kv_destroy((*h)->hp);

    ma_ug_destroy((*h)->ug);
    asg_destroy((*h)->r_g);
    free((*h));
}