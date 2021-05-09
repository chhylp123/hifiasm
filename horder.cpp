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
#define get_hit_suid(x, k) (((x).a.a[(k)].s<<1)>>(64 - (x).uID_bits))
#define get_hit_spos(x, k) ((x).a.a[(k)].s & (x).pos_mode)
#define get_hit_euid(x, k) (((x).a.a[(k)].e<<1)>>(64 - (x).uID_bits))
#define get_hit_epos(x, k) ((x).a.a[(k)].e & (x).pos_mode)

typedef struct {
	uint64_t ruid;
    uint64_t off;
} hit_aux_t;

typedef struct {
	hit_aux_t *a;
    size_t n, m;
    kvec_t(uint64_t) idx;
} u_hits_t;


#define hit_aux_ruid_key(x) ((x).ruid)
KRADIX_SORT_INIT(hit_aux_ruid, hit_aux_t, hit_aux_ruid_key, member_size(hit_aux_t, ruid))


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

void update_u_hits(kvec_pe_hit *u_hits, kvec_pe_hit *r_hits, ma_ug_t* ug, asg_t* r_g)
{
    u_hits_t x; memset(&x, 0, sizeof(x));
    hit_aux_t *p = NULL;
    ma_utg_t *u = NULL;
    pe_hit *t = NULL;
    uint64_t v, i, l, k, offset, occ_1, occ_2, *a_1, *a_2, i_1, i_2;
    
    for (v = 0; v < ug->g->n_seq; v++)
    {
        u = &(ug->u.a[v]);
        for (i = offset = 0; i < u->n; i++)
        {
            kv_pushp(hit_aux_t, x, &p);
            p->ruid = u->a[i]>>32; 
            p->ruid <<= 32;
            p->ruid |= v;
            p->off = offset;
            offset += (uint32_t)u->a[i];
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

void generate_haplotypes(horder_t *h, ug_opt_t *opt)
{
    uint64_t i, off;
    ma_ug_t *ug_1 = NULL, *ug_2 = NULL;
    asg_cleanup(h->r_g);
    asg_arc_del_trans(h->r_g, asm_opt.gap_fuzz);///must

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
    if(h->ug->g->n_seq != h->ug->u.n)
    {
        fprintf(stderr, "ERROR-non-equal-length\n");
    }
    fprintf(stderr, "h->ug->u.n-%u, off-%lu, ug_2->u.n-%u\n", (uint32_t)h->ug->u.n, off, (uint32_t)ug_2->u.n);
    for (i = 0; i < h->ug->g->n_seq; i++)
    {
        if(h->ug->g->seq[i].len != h->ug->u.a[i].len)
        {
            fprintf(stderr, "****ERROR-1-%lu: g->seq[i].len-%u, u.a[i].len-%u\n", 
            i, (uint32_t)h->ug->g->seq[i].len, (uint32_t)h->ug->u.a[i].len);
        }
        pu = i < off? (&ug_1->u.a[i]) : (&ug_2->u.a[i - off]);
        if(h->ug->g->seq[i].len != pu->len)
        {
            fprintf(stderr, "ERROR-2\n");
        } 
    }
    
    /*******************************for debug************************************/
    ug_1 = NULL;
    ma_ug_destroy(ug_2);
}

horder_t *init_horder_t(kvec_pe_hit *i_hits, uint64_t i_hits_uid_bits, uint64_t i_hits_pos_mode, 
asg_t *i_rg, ma_ug_t* i_ug, bubble_type* bub, ug_opt_t *opt)
{
    horder_t *h = NULL; CALLOC(h, 1);
    get_r_hits(i_hits, &(h->r_hits), i_rg, i_ug, bub, i_hits_uid_bits, i_hits_pos_mode);
    h->r_g = copy_read_graph(i_rg);
    generate_haplotypes(h, opt);
    update_u_hits(&(h->u_hits), &(h->r_hits), h->ug, h->r_g);    

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

    ma_ug_destroy((*h)->ug);
    asg_destroy((*h)->r_g);
    free((*h));
}