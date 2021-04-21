#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <stdlib.h>
///#include "partig.h"
#include "rcut.h"
#include "Purge_Dups.h"
#include "Correct.h"
#include "ksort.h"
#include "kthread.h"
#include "hic.h"

#define mc_edge_key(e) ((e).x)
KRADIX_SORT_INIT(mce, mc_edge_t, mc_edge_key, member_size(mc_edge_t, x))
#define mb_edge_key(e) ((e).x)
KRADIX_SORT_INIT(mbe, mb_edge_t, mb_edge_key, member_size(mb_edge_t, x))
#define mc_generic_key(x) (x)
KRADIX_SORT_INIT(mc64, uint64_t, mc_generic_key, 8)
KRADIX_SORT_INIT(mc32, uint32_t, mc_generic_key, 4)

#define pt_a(x, id) ((x).ma.a + ((x).idx.a[(id)]>>32))
#define pt_n(x, id) ((uint32_t)((x).idx.a[(id)]))
#define ma_x(z) (((z).x>>32))
#define ma_y(z) (((uint32_t)((z).x)))

uint8_t bit_filed[8] = {1, 2, 4, 8, 16, 32, 64, 128};
#define is_bit_set(id, a) ((a)[(id)>>3]&bit_filed[(id)&7])

typedef struct {
	int32_t max_iter;
    int32_t n_perturb;
    double f_perturb;
    uint64_t seed;
} mc_opt_t;


typedef struct {
	t_w_t z[2];
} mc_pairsc_t;

typedef struct {
	uint64_t x; // RNG
	uint32_t cc_off, cc_size;
	kvec_t(uint64_t) cc_edge;
	uint32_t *cc_node;
	uint32_t *bfs, *bfs_mark;
	mc_pairsc_t *z, *z_opt;///keep scores to nodes(1) and nodes(-1)
	int8_t *s, *s_opt;
	uint8_t *f;
} mc_svaux_t;

typedef struct{
	uint64_t chain_id, bid, uid;
}mc_bp_iter;

typedef struct{
    uint32_t chain_id;
    uint32_t f_bid;
    uint32_t f_uid;
    uint32_t l_bid; 
    uint32_t l_uid;
    uint32_t id;
    t_w_t w;
}mc_bp_res;

typedef struct{
	size_t n, m;
	uint8_t *a;
}bits_p;

typedef struct{
    bubble_type* b_b;
	uint64_t *idx, idx_n, occ, n_thread;
	mc_bp_res *res;
	uint8_t *lock;
	mc_svaux_t *b_aux;
	mc_match_t *ma;
	bits_p *vis;
}mc_bp_t;

typedef struct {
	uint64_t x; // RNG
	uint32_t cc_off, cc_size;
	kvec_t(uint32_t) cc_node;
	kvec_t(uint32_t) bfs;
	///share
	uint32_t *bfs_mark;
	mc_pairsc_t *z, *z_opt;///keep scores to nodes(1) and nodes(-1)
	int8_t *s, *s_opt;
} mc_svaux_t_s;

typedef struct {
	kvec_t(mc_svaux_t_s) bs;
	kvec_t(uint64_t) cc_edge;
	uint32_t *bfs_mark;
	mc_pairsc_t *z, *z_opt;///keep scores to nodes(1) and nodes(-1)
	int8_t *s, *s_opt;
	uint32_t n_t, n_g;
} mc_svaux_t_mul;

void mc_opt_init(mc_opt_t *opt, int32_t n_perturb, double f_perturb, uint64_t seed)
{
	memset(opt, 0, sizeof(mc_opt_t));
	// opt->n_perturb = 50000;
	opt->n_perturb = n_perturb;
	// opt->f_perturb = 0.1;
	opt->f_perturb = f_perturb;
	opt->max_iter = 1000;
	// opt->seed = 11;
	opt->seed = seed;
}

void mc_merge_dup(mc_g_t *mg) // MUST BE sorted
{
	uint32_t i, j, k, st;
	w_t w;
	for (st = 0, i = 1, k = 0; i <= mg->e->ma.n; ++i) {
		if (i == mg->e->ma.n || mg->e->ma.a[i].x != mg->e->ma.a[st].x) {
			if (i - st > 1) {
				for (j = st, w = 0; j < i; ++j) {
					w += mg->e->ma.a[j].w;
				}
				mg->e->ma.a[k] = mg->e->ma.a[st];
				mg->e->ma.a[k++].w = w;
			} else mg->e->ma.a[k++] = mg->e->ma.a[st];
			st = i;
		}
	}
	mg->e->ma.n = k;
}

void mb_merge_dup(mb_g_t *mbg) // MUST BE sorted
{
	uint32_t i, j, k, st;
	t_w_t w[4];
	for (st = 0, i = 1, k = 0; i <= mbg->e->ma.n; ++i) {
		if (i == mbg->e->ma.n || mbg->e->ma.a[i].x != mbg->e->ma.a[st].x) {
			if (i - st > 1) {
				w[0] = w[1] = w[2] = w[3] = 0;
				for (j = st; j < i; ++j) {
					w[0] += mbg->e->ma.a[j].w[0];
					w[1] += mbg->e->ma.a[j].w[1];
					w[2] += mbg->e->ma.a[j].w[2];
					w[3] += mbg->e->ma.a[j].w[3];
				}
				mbg->e->ma.a[k] = mbg->e->ma.a[st];
				mbg->e->ma.a[k].w[0] = w[0];
				mbg->e->ma.a[k].w[1] = w[1];
				mbg->e->ma.a[k].w[2] = w[2];
				mbg->e->ma.a[k].w[3] = w[3];
				k++;
			} else mbg->e->ma.a[k++] = mbg->e->ma.a[st];
			st = i;
		}
	}
	mbg->e->ma.n = k;
}

static void mc_edges_idx(mc_match_t *ma)
{
	uint32_t st, i;
	kv_resize(uint64_t, ma->idx, ma->n_seq);
	ma->idx.n = ma->n_seq;
	memset(ma->idx.a, 0, ma->idx.n*sizeof(uint64_t));
	for (st = 0, i = 1; i <= ma->ma.n; ++i)
		if (i == ma->ma.n || (ma->ma.a[i].x>>32) != (ma->ma.a[st].x>>32))
			ma->idx.a[ma->ma.a[st].x>>32] = (uint64_t)st << 32 | (i - st), st = i;
}

static void mb_edges_idx(mb_match_t *ma)
{
	uint32_t st, i;
	kv_resize(uint64_t, ma->idx, ma->n_seq);
	ma->idx.n = ma->n_seq;
	memset(ma->idx.a, 0, ma->idx.n*sizeof(uint64_t));
	for (st = 0, i = 1; i <= ma->ma.n; ++i)
		if (i == ma->ma.n || (ma->ma.a[i].x>>32) != (ma->ma.a[st].x>>32))
			ma->idx.a[ma->ma.a[st].x>>32] = (uint64_t)st << 32 | (i - st), st = i;
}


mc_g_t *init_mc_g_t(ma_ug_t *ug, asg_t *read_g, int8_t *s, uint32_t renew_s)
{
	mc_g_t *p = NULL; CALLOC(p, 1);
	p->ug = ug;
	p->rg = read_g;
	kv_init(p->s);
	if(s)
	{
		p->s.a = s;
		p->s.n = ug->g->n_seq;
		p->s.m = 0;
		if(renew_s) memset(p->s.a, 0, p->s.n);
	}
	else
	{
		CALLOC(p->s.a, ug->g->n_seq);
		p->s.n = p->s.m = ug->g->n_seq;
	} 
	return p;
}

static mc_edge_t *get_mc_edge(const mc_match_t *ma, uint32_t sid1, uint32_t sid2)
{
	mc_edge_t *o = pt_a(*ma, sid1);
	uint32_t n = pt_n(*ma, sid1), k;
	for (k = 0; k < n; ++k)
		if (((uint32_t)o[k].x) == sid2)
			return &(o[k]);
	return NULL;
}

static mb_edge_t *get_mb_edge(const mb_match_t *ma, uint32_t sid1, uint32_t sid2)
{
	mb_edge_t *o = pt_a(*ma, sid1);
	uint32_t n = pt_n(*ma, sid1), k;
	for (k = 0; k < n; ++k)
		if (((uint32_t)o[k].x) == sid2)
			return &(o[k]);
	return NULL;
}

uint32_t mb_edges_symm(mb_match_t *ma);

uint32_t debug_mb_edges_symm(mb_match_t *ma)
{
	uint32_t i, n = 0;
	mb_edge_t *t = NULL, *m = NULL;

	for (i = 0; i < ma->ma.n; ++i) {
		m = &ma->ma.a[i];
		if (ma_x(*m) == ma_y(*m))
		{
			fprintf(stderr, "ERROR-0\n"); 
			continue;
		}
		t = get_mb_edge(ma, ma_y(*m), ma_x(*m));
		if(!t)
		{
			fprintf(stderr, "ERROR-1\n"); 
			continue;
		} 

		if(m->w[0] != t->w[0]) fprintf(stderr, "ERROR-2\n"); 
		if(m->w[3] != t->w[3]) fprintf(stderr, "ERROR-3\n"); 
		if(m->w[1] != t->w[2]) fprintf(stderr, "ERROR-4\n"); 
		if(m->w[2] != t->w[1]) fprintf(stderr, "ERROR-5\n"); 
	}

	return n;
}

inline void decode_mb_node(mb_g_t *mbg, uint32_t id, uint32_t **a0, uint32_t *n0, mc_node_t *s0, 
uint32_t **a1, uint32_t *n1, mc_node_t *s1)
{
	if(a0) (*a0) = mbg->u->bid.a + mbg->u->u.a[id].a[0];
	if(n0) (*n0) = mbg->u->u.a[id].occ[0];
	if(s0) (*s0) = mbg->u->u.a[id].s[0];

	if(a1) (*a1) = mbg->u->bid.a + mbg->u->u.a[id].a[1];
	if(n1) (*n1) = mbg->u->u.a[id].occ[1];
	if(s1) (*s1) = mbg->u->u.a[id].s[1];
}

mb_g_t *init_mb_g_t(mc_g_t *mg, mb_nodes_t* u, uint32_t is_sys)
{
	if(u == NULL || mg == NULL) return NULL;
	mc_edge_t *o = NULL;
	uint32_t i, k, m, n, a_n[2], *a[2], qn, tn, qb, tb;
	mb_edge_t *t = NULL;
	mb_g_t *p = NULL; CALLOC(p, 1);
	p->u = u;
	p->e = NULL; CALLOC(p->e, 1);
	p->e->n_seq = u->u.n;
	kv_init(p->e->ma); kv_init(p->e->idx);
	for (i = 0; i < p->u->u.n; i++)///each block
	{
		qb = i;
		p->u->u.a[i].s[0] = p->u->u.a[i].s[1] = 0;
		decode_mb_node(p, i, &(a[0]), &(a_n[0]), NULL, &(a[1]), &(a_n[1]), NULL);
		for (k = 0; k < a_n[0]; k++)
		{
			qn = a[0][k];
			o = pt_a(*(mg->e), qn);
    		n = pt_n(*(mg->e), qn);
			for (m = 0; m < n; m++)
			{
				tn = ma_y(o[m]);
				tb = p->u->idx.a[tn]>>1;
				if(tb == qb)
				{
					continue;
				} 
				
				kv_pushp(mb_edge_t, p->e->ma, &t);
				t->x = (uint64_t)qb << 32 | tb;
				t->w[0] = t->w[1] = t->w[2] = t->w[3] = 0;
				t->w[p->u->idx.a[tn]&1] = o[m].w;
			}
		}

		for (k = 0; k < a_n[1]; k++)
		{
			qn = a[1][k];
			o = pt_a(*(mg->e), qn);
    		n = pt_n(*(mg->e), qn);
			for (m = 0; m < n; m++)
			{
				tn = ma_y(o[m]);
				tb = p->u->idx.a[tn]>>1;
				if(tb == qb)
				{
					continue;
				} 

				kv_pushp(mb_edge_t, p->e->ma, &t);
				t->x = (uint64_t)qb << 32 | tb;
				t->w[0] = t->w[1] = t->w[2] = t->w[3] = 0;
				t->w[(p->u->idx.a[tn]&1)+2] = o[m].w;
			}
		}
	}

	radix_sort_mbe(p->e->ma.a, p->e->ma.a + p->e->ma.n);
	mb_merge_dup(p); // MUST BE sorted
	mb_edges_idx(p->e);
	/*******************************for debug************************************/
	debug_mb_edges_symm(p->e);
	/*******************************for debug************************************/
	if(is_sys) mb_edges_symm(p->e);

	return p;
}

void destory_mc_g_t(mc_g_t **p)
{
	if(!p || !(*p)) return;
	if((*p)->s.m == 0) (*p)->s.a = NULL;
	kv_destroy((*p)->s);
	if((*p)->e)
	{
		kv_destroy((*p)->e->idx);
		kv_destroy((*p)->e->ma);
		free((*p)->e->cc);
		free((*p)->e);
	}
	free((*p));
}

void destory_mb_g_t(mb_g_t **p)
{
	if(!p || !(*p)) return;
	kv_destroy((*p)->e->idx);
	kv_destroy((*p)->e->ma);
	free((*p)->e->cc);
	free((*p)->e);
	free((*p));
}


static void ks_shuffle_uint32_t(size_t n, uint32_t a[], uint64_t *x)
{
	size_t i, j;
	for (i = n; i > 1; --i) {
		uint32_t tmp;
		j = (size_t)(kr_drand_r(x) * i);
		tmp = a[j]; a[j] = a[i-1]; a[i-1] = tmp;
	}
}



static void normalize_mc_edge(mc_edge_t *a, mc_edge_t *b)
{
	if(a->w >= b->w)
	{
		b->x = (uint32_t)a->x;
		b->x <<= 32;
		b->x |= (a->x>>32);
		b->w = a->w;
	}
	else
	{
		a->x = (uint32_t)b->x;
		a->x <<= 32;
		a->x |= (b->x>>32);
		a->w = b->w;
	}
}

uint32_t mc_edges_symm(mc_match_t *ma)
{
	uint8_t *del = NULL;
	uint32_t i, k, n = 0;
	mc_edge_t *t = NULL, *m = NULL;
	CALLOC(del, ma->ma.n);

	for (i = 0; i < ma->ma.n; ++i) {
		m = &ma->ma.a[i];
		if (ma_x(*m) == ma_y(*m))
		{
			del[i] = 1, ++n;///self overlap
			continue;
		}
		t = get_mc_edge(ma, ma_y(*m), ma_x(*m));
		if(!t)
		{
			del[i] = 1, ++n;///self overlap
			continue;
		} 
		normalize_mc_edge(m, t);
	}

	if (n > 0) {
		for (i = k = 0; i < ma->ma.n; ++i)
			if (!del[i]) ma->ma.a[k++] = ma->ma.a[i];
		ma->ma.n = k;
		mc_edges_idx(ma);
	}

	free(del);
	return n;
}


static void normalize_mb_edge(mb_edge_t *a, mb_edge_t *b)
{
	if(a->w >= b->w)
	{
		b->x = (uint32_t)a->x;
		b->x <<= 32;
		b->x |= (a->x>>32);
		b->w[0] = a->w[0];
		b->w[3] = a->w[3];
		b->w[1] = a->w[2];
		b->w[2] = a->w[1];
	}
	else
	{
		a->x = (uint32_t)b->x;
		a->x <<= 32;
		a->x |= (b->x>>32);
		a->w[0] = b->w[0];
		a->w[3] = b->w[3];
		a->w[1] = b->w[2];
		a->w[2] = b->w[1];
	}
}

uint32_t mb_edges_symm(mb_match_t *ma)
{
	uint8_t *del = NULL;
	uint32_t i, k, n = 0;
	mb_edge_t *t = NULL, *m = NULL;
	CALLOC(del, ma->ma.n);

	for (i = 0; i < ma->ma.n; ++i) {
		m = &ma->ma.a[i];
		if (ma_x(*m) == ma_y(*m))
		{
			del[i] = 1, ++n;///self overlap
			continue;
		}
		t = get_mb_edge(ma, ma_y(*m), ma_x(*m));
		if(!t)
		{
			del[i] = 1, ++n;///self overlap
			continue;
		} 
		normalize_mb_edge(m, t);
	}

	if (n > 0) {
		for (i = k = 0; i < ma->ma.n; ++i)
			if (!del[i]) ma->ma.a[k++] = ma->ma.a[i];
		ma->ma.n = k;
		mb_edges_idx(ma);
	}

	free(del);
	return n;
}

void debug_mc_interval_t(mc_interval_t *p, uint32_t p_n, uint32_t *p_idx, ma_ug_t *ug, asg_t *rg, 
trans_chain* t_ch)
{
	fprintf(stderr, "0----------[M::%s]----------\n", __func__);
	uint32_t i, offset, v, sid, eid, spos, epos, p_status, p_uid, occ;
    ma_utg_t *u = NULL;
	mc_interval_t *a = NULL, *t = NULL;

	for (v = 0; v < ug->u.n; v++)
    {
		a = p + p_idx[v]; 
		occ = p_idx[v+1] - p_idx[v];
        for (i = 0; i < occ; i++)
        {
            if(a[i].uID != v) fprintf(stderr, "ERROR-s\n");
        }
    }

    for (v = 0, p_status = (uint32_t)-1, p_uid = (uint32_t)-1; v < p_n; v++)
    {
        t = &(p[v]);
        sid = t->nS;
        eid = t->nE;
        spos = t->bS;
        epos = t->bE;
        if(p_uid == t->uID && p_status == t->hs)
        {
            fprintf(stderr, "ERROR-a\n");
        } 
        p_status = t->hs;
        p_uid = t->uID;
        u = &(ug->u.a[t->uID]);
        for (i = offset = 0; i < u->n; i++)
        {
            if(i == sid)
            {
                if(spos != offset)
                {
                    fprintf(stderr, "ERROR-b\n");
                }
            }

            if(i == eid)
            {
                if(epos != (offset+rg->seq[u->a[i]>>33].len - 1))
                {
                    fprintf(stderr, "ERROR-c, real end: %u\n", 
                                    (uint32_t)(offset+rg->seq[u->a[i]>>33].len - 1));
                }
            }

            offset += (uint32_t)u->a[i];
            if(i >= sid && i <= eid)
            {
                if(t_ch->is_r_het[u->a[i]>>33] != t->hs)
                {
                    fprintf(stderr, "ERROR-d: is_r_het: %u, h_status: %u\n", t_ch->is_r_het[u->a[i]>>33], t->hs);
                }
            }
        }
    }
}

void update_mc_edges(mc_g_t *mg, hap_overlaps_list* ha, kv_u_trans_t *ta, trans_chain* t_ch, double f_rate, uint32_t is_sys)
{
	uint32_t v, i, k, qn, tn, qs, qe, ts, te, occ, as, ae, l, offset, l_pos;
	uint64_t hetLen, homLen, oLen;
	mc_interval_t *a = NULL;
	mc_edge_t *ma = NULL;
	asg_t* nsg = mg->ug->g;
	ma_utg_t *u = NULL;
	mc_interval_t *t = NULL;

	kvec_t(mc_interval_t) p; kv_init(p);
	kvec_t(uint32_t) p_idx; kv_init(p_idx);
	if(t_ch)
	{
		kv_push(uint32_t, p_idx, 0);
		for (v = 0; v < nsg->n_seq; v++)
		{
			u = &(mg->ug->u.a[v]);
			for (k = 1, l = 0, offset = 0, l_pos = 0; k <= u->n; ++k) 
			{   
				if (k == u->n || t_ch->is_r_het[u->a[k]>>33] != t_ch->is_r_het[u->a[l]>>33])
				{
					kv_pushp(mc_interval_t, p, &t);
					t->uID = v;
					t->hs = t_ch->is_r_het[u->a[l]>>33];

					t->bS = l_pos;
					t->bE = offset + mg->rg->seq[u->a[k-1]>>33].len - 1;

					t->nS = l;
					t->nE = k - 1;
					l = k;
					l_pos = offset + (uint32_t)u->a[k-1];
				}
				offset += (uint32_t)u->a[k-1];
			}
			kv_push(uint32_t, p_idx, p.n);
		}

		///debug_mc_interval_t(p.a, p.n, p_idx.a, mg->ug, mg->rg, t_ch);
	}

	if(!mg->e)
	{
		CALLOC(mg->e, 1);
		mg->e->n_seq = mg->ug->g->n_seq;
		kv_init(mg->e->idx); kv_init(mg->e->ma); 
	}

	if(ha)
	{
		for (v = 0; v < ha->num; v++)
		{
			for (i = 0; i < ha->x[v].a.n; i++)
			{
				if(ha->x[v].a.a[i].score <= 0) continue;
				if(ha->x[v].a.a[i].xUid == ha->x[v].a.a[i].yUid) continue;
				if(p.n > 0 && p_idx.n > 0)
				{
					/*****************qn*****************/
					qn = ha->x[v].a.a[i].xUid;
					qs = ha->x[v].a.a[i].x_beg_pos;
					qe = ha->x[v].a.a[i].x_end_pos - 1;

					a = p.a + p_idx.a[qn]; 
					occ = p_idx.a[qn+1] - p_idx.a[qn];
					for (k = 0, hetLen = 0, homLen = 0; k < occ; k++)
					{
						as = a[k].bS;
						ae = a[k].bE;
						oLen = ((MIN(qe, ae) >= MAX(qs, as))? MIN(qe, ae) - MAX(qs, as) + 1 : 0);
						if(homLen + hetLen > 0 && oLen == 0) break;
						if(oLen == 0) continue;
						if(a[k].hs == N_HET)
						{
							homLen += oLen;
						} 
						else if(asm_opt.polyploidy <= 2 && (a[k].hs&P_HET))///if(asm_opt.polyploidy <= 2 && (a[k].hs&S_HET))
						{
							homLen += oLen;
						}
						else
						{
							hetLen += oLen;
						}                
					}

					if(hetLen <= ((hetLen + homLen)*f_rate)) continue;
					/*****************qn*****************/


					/*****************tn*****************/
					tn = ha->x[v].a.a[i].yUid;
					ts = ha->x[v].a.a[i].y_beg_pos;
					te = ha->x[v].a.a[i].y_end_pos - 1;

					a = p.a + p_idx.a[tn]; 
					occ = p_idx.a[tn+1] - p_idx.a[tn];
					for (k = 0, hetLen = 0, homLen = 0; k < occ; k++)
					{
						as = a[k].bS;
						ae = a[k].bE;
						oLen = ((MIN(te, ae) >= MAX(ts, as))? MIN(te, ae) - MAX(ts, as) + 1 : 0);
						if(homLen + hetLen > 0 && oLen == 0) break;
						if(oLen == 0) continue;
						if(a[k].hs == N_HET)
						{
							homLen += oLen;
						} 
						else if(asm_opt.polyploidy <= 2 && (a[k].hs&P_HET))///if(asm_opt.polyploidy <= 2 && (a[k].hs&S_HET))
						{
							homLen += oLen;
						}
						else
						{
							hetLen += oLen;
						}                 
					}

					if(hetLen <= ((hetLen + homLen)*f_rate)) continue;
					/*****************tn*****************/
				}
				kv_pushp(mc_edge_t, mg->e->ma, &ma);
				ma->x = (uint64_t)ha->x[v].a.a[i].xUid << 32 | ha->x[v].a.a[i].yUid;
				ma->w = ha->x[v].a.a[i].score;
			}
		}
	}

	if(ta)
	{
		for (i = 0; i < ta->n; ++i)
		{
			if(ta->a[i].del) continue;
			if(p.n > 0 && p_idx.n > 0)
			{
				/*****************qn*****************/
				qn = ta->a[i].qn;
				qs = ta->a[i].qs;
				qe = ta->a[i].qe - 1;

				a = p.a + p_idx.a[qn]; 
				occ = p_idx.a[qn+1] - p_idx.a[qn];
				for (k = 0, hetLen = 0, homLen = 0; k < occ; k++)
				{
					as = a[k].bS;
					ae = a[k].bE;
					oLen = ((MIN(qe, ae) >= MAX(qs, as))? MIN(qe, ae) - MAX(qs, as) + 1 : 0);
					if(homLen + hetLen > 0 && oLen == 0) break;
					if(oLen == 0) continue;
					if(a[k].hs == N_HET)
					{
						homLen += oLen;
					} 
					else if(asm_opt.polyploidy <= 2 && (a[k].hs&P_HET))///if(asm_opt.polyploidy <= 2 && (a[k].hs&S_HET))
					{
						homLen += oLen;
					}
					else
					{
						hetLen += oLen;
					}                
				}

				if(hetLen <= ((hetLen + homLen)*f_rate)) continue;
				/*****************qn*****************/

				/*****************tn*****************/
				tn = ta->a[i].tn;
				ts = ta->a[i].ts;
				te = ta->a[i].te - 1;

				a = p.a + p_idx.a[tn]; 
				occ = p_idx.a[tn+1] - p_idx.a[tn];
				for (k = 0, hetLen = 0, homLen = 0; k < occ; k++)
				{
					as = a[k].bS;
					ae = a[k].bE;
					oLen = ((MIN(te, ae) >= MAX(ts, as))? MIN(te, ae) - MAX(ts, as) + 1 : 0);
					if(homLen + hetLen > 0 && oLen == 0) break;
					if(oLen == 0) continue;
					if(a[k].hs == N_HET)
					{
						homLen += oLen;
					} 
					else if(asm_opt.polyploidy <= 2 && (a[k].hs&P_HET))///if(asm_opt.polyploidy <= 2 && (a[k].hs&S_HET))
					{
						homLen += oLen;
					}
					else
					{
						hetLen += oLen;
					}                 
				}

				if(hetLen <= ((hetLen + homLen)*f_rate)) continue;
				/*****************tn*****************/
			}
			kv_pushp(mc_edge_t, mg->e->ma, &ma);
			ma->x = (uint64_t)ta->a[i].qn << 32 | ta->a[i].tn;
			ma->w = w_cast(ta->a[i].nw);
		}
	}
	
	radix_sort_mce(mg->e->ma.a, mg->e->ma.a + mg->e->ma.n);
	mc_merge_dup(mg);
	mc_edges_idx(mg->e);
	if(is_sys) mc_edges_symm(mg->e);
	kv_destroy(p); kv_destroy(p_idx);
}

void debug_mc_g_t(mc_g_t *mg)
{
    fprintf(stderr, "0----------[M::%s]----------\n", __func__);
	mc_edge_t *o = NULL, *s = NULL;
	uint32_t i, k, n, cnt;
	for (i = 0; i < mg->e->n_seq; ++i) 
	{
		o = pt_a(*(mg->e), i); n = pt_n(*(mg->e), i);
		for (k = 0; k < n; ++k)
		{
			if(ma_x(o[k]) != i) fprintf(stderr, "ERROR-g\n");
			s = get_mc_edge(mg->e, ma_y(o[k]), ma_x(o[k]));
			if(!s) fprintf(stderr, "ERROR-e\n");
			if(s)
			{
				if(!(ma_x(*s) == ma_y(o[k]) && ma_y(*s) == ma_x(o[k]) && s->w == o[k].w))
				{
					fprintf(stderr, "ERROR-f\n");
				}
			}
		}

		for (k = cnt = 0; k < mg->e->ma.n; ++k)
		{
			if(ma_x(mg->e->ma.a[k]) == i) cnt++;
		}

		if(cnt != n) fprintf(stderr, "ERROR-h\n");
	}
}

uint64_t *mc_g_cc_core(mc_match_t *ma)
{
	uint32_t i, x, y, *flag;
	uint64_t *group;
	mc_edge_t *o = NULL;
	kvec_t(uint32_t) stack; kv_init(stack);

	MALLOC(flag, ma->n_seq);
	for (i = 0; i < ma->n_seq; ++i)
		flag[i] = (uint32_t)-1;

	// connected componets
	for (i = 0; i < ma->n_seq; ++i) {
		if (flag[i] != (uint32_t)-1) continue;
		stack.n = 0;
		kv_push(uint32_t, stack, i);
		while (stack.n > 0) {
			uint32_t k, j, n;
			stack.n--;
			k = stack.a[stack.n];
			flag[k] = i;///group id
			// n = (uint32_t)ma->idx[k];
			// s = ma->idx[k] >> 32;
			o = pt_a(*ma, k); 
			n = pt_n(*ma, k);
			for (j = 0; j < n; ++j) {
				uint32_t t = ma_y(o[j]);
				if (flag[t] != (uint32_t)-1) continue;
				// if (ns == ms) PT_EXPAND(stack, ms);
				// stack[ns++] = t;
				kv_push(uint32_t, stack, t);
			}
		}
	}
	kv_destroy(stack);

	// precalculate the size of each group
	CALLOC(group, ma->n_seq);
	for (i = 0; i < ma->n_seq; ++i)
		group[i] = (uint64_t)flag[i] << 32 | i;
	radix_sort_mc64(group, group + ma->n_seq);
	for (i = 1, x = y = 0; i <= ma->n_seq; ++i) {
		if (i == ma->n_seq || group[i]>>32 != group[x]>>32) {
			uint32_t j;
			for (j = x; j < i; ++j)
				group[j] = (uint64_t)y << 32 | (uint32_t)group[j];///(group id)|first element in this group
			++y, x = i;
		}
	}
	free(flag);
	return group;
}

void mc_g_cc(mc_match_t *ma)
{
	ma->cc = mc_g_cc_core(ma);
}
mc_bp_t *mc_bp_t_init(mc_match_t *ma, mc_svaux_t *b_aux, bubble_type* bub, uint64_t n_thread)
{
	uint32_t i, k, n;
	mc_bp_t *bp = NULL; 
	ma_utg_t *u = NULL;
	CALLOC(bp, 1);
	bp->b_aux = b_aux;
	bp->ma = ma;
	bp->b_b = bub;
	bp->n_thread = n_thread;
	CALLOC(bp->lock, bub->ug->g->n_seq);
	CALLOC(bp->res, n_thread);
	CALLOC(bp->vis, n_thread);
	for (i = 0; i < n_thread; i++)
	{
		bp->vis[i].m = bp->vis[i].n = bub->ug->g->n_seq;
		CALLOC(bp->vis[i].a, bp->vis[i].n);
	}
	

	MALLOC(bp->idx, bub->chain_weight.n+1);
	bp->idx_n = bp->occ = 0;
	for (i = 0; i < bub->chain_weight.n; i++)
	{
		bp->idx[i] = bp->occ;
		bp->idx_n++;
		if(bub->chain_weight.a[i].del) continue;
		u = &(bub->b_ug->u.a[bub->chain_weight.a[i].id]);///list of bubbles
		for (k = 0; k < u->n; k++)
		{
			get_bubbles(bub, u->a[k]>>33, NULL, NULL, NULL, &n, NULL);
			bp->occ += n;
		}
	}
	bp->idx[i] = bp->occ;
	fprintf(stderr, "# nodes in chains: %lu, # chains: %lu\n", bp->occ, bp->idx_n);
	return bp;
}

void destroy_mc_bp_t(mc_bp_t **bp)
{
	uint32_t i;
	for (i = 0; i < (*bp)->n_thread; i++)
	{
		free((*bp)->vis[i].a);
	}
	free((*bp)->idx);
	free((*bp)->res);
	free((*bp)->lock);
	free((*bp));
}

mc_svaux_t *mc_svaux_init(const mc_g_t *mg, uint64_t x)
{
	uint32_t st, i, max_cc = 0;
	mc_match_t *ma = mg->e;
	mc_svaux_t *b;
	CALLOC(b, 1);
	b->x = x;
	for (st = 0, i = 1; i <= ma->n_seq; ++i)
		if (i == ma->n_seq || ma->cc[st]>>32 != ma->cc[i]>>32)
			max_cc = max_cc > i - st? max_cc : i - st, st = i;
	kv_init(b->cc_edge);
	MALLOC(b->cc_node, max_cc);
	///CALLOC(b->s, ma->n_seq);
	b->s = mg->s.a;
	CALLOC(b->s_opt, ma->n_seq);
	MALLOC(b->bfs, ma->n_seq);

	MALLOC(b->bfs_mark, ma->n_seq); 
	memset(b->bfs_mark, -1, ma->n_seq*sizeof(uint32_t));
	
	CALLOC(b->z, ma->n_seq);
	CALLOC(b->z_opt, ma->n_seq);
	CALLOC(b->f, ma->n_seq);
	return b;
}
void mc_svaux_destroy(mc_svaux_t *b)
{
	b->s = NULL;
	kv_destroy(b->cc_edge); free(b->cc_node);
	free(b->s); free(b->s_opt);
	free(b->z); free(b->z_opt);
	free(b->bfs); free(b->bfs_mark);
	free(b->f);
	free(b);
}


mc_svaux_t_mul *init_mc_svaux_t_mul(const mc_g_t *mg, uint64_t n_threads)
{
	uint32_t st, i;
	mc_match_t *ma = mg->e;
	mc_svaux_t_mul *b; CALLOC(b, 1);
	for (st = 0, i = 1, b->n_t = 0; i <= ma->n_seq; ++i)
	{
		if (i == ma->n_seq || ma->cc[st]>>32 != ma->cc[i]>>32)
		{
			b->n_t++;
		}
	}
	b->n_g = b->n_t;
	if(n_threads < b->n_t) b->n_t = n_threads;
	
	kv_init(b->cc_edge);
	MALLOC(b->bfs_mark, ma->n_seq);
	memset(b->bfs_mark, -1, ma->n_seq*sizeof(uint32_t));
	CALLOC(b->s_opt, ma->n_seq);
	CALLOC(b->z, ma->n_seq);
	CALLOC(b->z_opt, ma->n_seq);
	b->s = mg->s.a;

	kv_init(b->bs); CALLOC(b->bs.a, b->n_t); 
	b->bs.n = b->bs.m = b->n_t;
	for (i = 0; i < b->bs.n; i++)
	{
		kv_init(b->bs.a[i].cc_node);
		kv_init(b->bs.a[i].bfs);
		b->bs.a[i].bfs_mark = b->bfs_mark;
		b->bs.a[i].z = b->z;
		b->bs.a[i].z_opt = b->z_opt;
		b->bs.a[i].s = b->s;
		b->bs.a[i].s_opt = b->s_opt;
	}
	return b;
}


void destroy_mc_svaux_t_mul(mc_svaux_t_mul *b)
{
	uint32_t i;
	kv_destroy(b->cc_edge);
	free(b->bfs_mark);
	free(b->s_opt);
	free(b->z);
	free(b->z_opt);
	b->s = NULL;
	for (i = 0; i < b->bs.n; i++)
	{
		kv_destroy(b->bs.a[i].cc_node);
		kv_destroy(b->bs.a[i].bfs);
	}
	kv_destroy(b->bs);
	free(b);
}

uint32_t mc_best(const mc_match_t *ma, mc_svaux_t *b)
{
	uint32_t i, max_i = (uint32_t)-1;
	t_w_t w, max_w;
	for (i = 0, max_w = -1; i < b->cc_size; ++i) {
		uint32_t k = (uint32_t)ma->cc[b->cc_off + i];///uid
		if(b->f[k] || b->s[k] == 0) continue;
		///z += -((t_w_t)(b->s[k])) * (b->z[k].z[0] - b->z[k].z[1]);
		///-((t_w_t)(b->s[k])) * (b->z[k].z[0] - b->z[k].z[1]) current
		///((t_w_t)(b->s[k])) * (b->z[k].z[0] - b->z[k].z[1]) flipped
		w = ((t_w_t)(b->s[k])) * (b->z[k].z[0] - b->z[k].z[1]) * 2;
		if(w <= 0) continue;
		if(w > max_w)
		{
			w = max_w; max_i = k;
		}
	}

	return max_i;
}


t_w_t mc_score(const mc_match_t *ma, mc_svaux_t *b)
{
	uint32_t i;
	t_w_t z = 0;
	for (i = 0; i < b->cc_size; ++i) {
		uint32_t k = (uint32_t)ma->cc[b->cc_off + i];///uid
		z += -((t_w_t)(b->s[k])) * (b->z[k].z[0] - b->z[k].z[1]);
	}
	return z;
}

t_w_t mc_score_all(const mc_match_t *ma, mc_svaux_t *b)
{
	uint32_t k;
	t_w_t z = 0;
	for (k = 0; k < ma->n_seq; ++k) 
	{
		z += -((t_w_t)(b->s[k])) * (b->z[k].z[0] - b->z[k].z[1]);
	}
	return z;
}

void mc_reset_z(const mc_match_t *ma, mc_svaux_t *b)
{
	uint32_t i;
	for (i = 0; i < b->cc_size; ++i) {
		uint32_t k = (uint32_t)ma->cc[b->cc_off + i];///uid
		uint32_t o = ma->idx.a[k] >> 32;
		uint32_t j, n = (uint32_t)ma->idx.a[k];
		b->z[k].z[0] = b->z[k].z[1] = 0;
		for (j = 0; j < n; ++j) {
			const mc_edge_t *e = &ma->ma.a[o + j];
			uint32_t t = ma_y(*e);
			if (b->s[t] > 0) b->z[k].z[0] += e->w;
			else if (b->s[t] < 0) b->z[k].z[1] += e->w;
		}
	}
}

void mc_reset_z_debug(const mc_match_t *ma, mc_svaux_t *b)
{
	uint32_t i;
	t_w_t z[2];
	for (i = 0; i < b->cc_size; ++i) {
		uint32_t k = (uint32_t)ma->cc[b->cc_off + i];///uid
		uint32_t o = ma->idx.a[k] >> 32;
		uint32_t j, n = (uint32_t)ma->idx.a[k];
		z[0] = b->z[k].z[0]; z[1] = b->z[k].z[1];
		b->z[k].z[0] = b->z[k].z[1] = 0;
		for (j = 0; j < n; ++j) {
			const mc_edge_t *e = &ma->ma.a[o + j];
			uint32_t t = ma_y(*e);
			if (b->s[t] > 0) b->z[k].z[0] += e->w;
			else if (b->s[t] < 0) b->z[k].z[1] += e->w;
		}
		if(z[0] != b->z[k].z[0]) fprintf(stderr, "ERROR1\n");
		if(z[1] != b->z[k].z[1]) fprintf(stderr, "ERROR2\n");
	}
}


t_w_t mc_init_spin(const mc_match_t *ma, mc_svaux_t *b)
{
	uint32_t i;
	b->cc_edge.n = 0;
	if(b->cc_size <= 2)//mannually clear
	{
		for (i = 0; i < b->cc_size; ++i) {///how many nodes
			b->s[(uint32_t)ma->cc[b->cc_off + i]] = 0;
		}
	}

	for (i = 0; i < b->cc_size; ++i) {///how many nodes
		uint32_t k = (uint32_t)ma->cc[b->cc_off + i];///node id
		b->cc_node[i] = k;
		if(b->s[k] == 0) break;
	}
	if(i >= b->cc_size) goto passed;

	for (i = 0; i < b->cc_size; ++i) {///how many nodes
		uint32_t k = (uint32_t)ma->cc[b->cc_off + i];///node id
		uint32_t o = ma->idx.a[k] >> 32;///cc group id
		uint32_t n = (uint32_t)ma->idx.a[k], j;
		b->cc_node[i] = k;
		for (j = 0; j < n; ++j) {
			w_t w = ma->ma.a[o + j].w;
			w = w > 0? w : -w;
			kv_push(uint64_t, b->cc_edge, (uint64_t)((uint32_t)-1 - ((uint32_t)w)) << 32 | (o + j));
		}
	}
	radix_sort_mc64(b->cc_edge.a, b->cc_edge.a + b->cc_edge.n);
	for (i = 0; i < b->cc_edge.n; ++i) { // from the strongest edge to the weakest
		const mc_edge_t *e = &ma->ma.a[(uint32_t)b->cc_edge.a[i]];
		uint32_t n1 = ma_x(*e), n2 = ma_y(*e);
		if (b->s[n1] == 0 && b->s[n2] == 0) {
			b->x = kr_splitmix64(b->x);
			b->s[n1] = b->x&1? 1 : -1;
			b->s[n2] = e->w > 0? -b->s[n1] : b->s[n1];
		}/****************************may have bugs********************************/
		else if(b->s[n1] == 0)
		{
			b->s[n1] = e->w > 0? -b->s[n2] : b->s[n2];
		}
		else if(b->s[n2] == 0)
		{
			b->s[n2] = e->w > 0? -b->s[n1] : b->s[n1];
		}
		/****************************may have bugs********************************/
	}

	passed:
	mc_reset_z(ma, b);
	return mc_score(ma, b);
}
///k is uid
static void mc_set_spin(const mc_match_t *ma, mc_svaux_t *b, uint32_t k, int8_t s)
{
	uint32_t o, j, n;
	int8_t s0 = b->s[k];
	if (s0 == s) return;
	o = ma->idx.a[k] >> 32;
	n = (uint32_t)ma->idx.a[k];
	for (j = 0; j < n; ++j) {
		const mc_edge_t *e = &ma->ma.a[o + j];
		uint32_t t = ma_y(*e);///1->z[0]; (-1)->z[1];
		if (s0 != 0) b->z[t].z[(s0 < 0)] -= e->w; 
		if (s  != 0) b->z[t].z[(s  < 0)] += e->w;
	}
	b->s[k] = s;
}

void mc_best_flip(const mc_match_t *ma, mc_svaux_t *b)
{
	uint32_t idx;
	for (idx = 0; idx < b->cc_size; ++idx) {
		b->f[(uint32_t)ma->cc[b->cc_off + idx]] = 0;
        ///uint32_t k = (uint32_t)ma->cc[b->cc_off + idx];///uid
	}
	while (1)
	{
		idx = mc_best(ma, b);
		if(idx == (uint32_t)-1) break;
		mc_set_spin(ma, b, idx, -b->s[idx]);
		b->f[idx] = 1;
	}
}

static t_w_t mc_optimize_local(const mc_opt_t *opt, const mc_match_t *ma, mc_svaux_t *b, uint32_t *n_iter)
{
	uint32_t i, n_flip = 0;
	int32_t n_iter_local = 0;
	while (n_iter_local < opt->max_iter) {
		++(*n_iter);
		ks_shuffle_uint32_t(b->cc_size, b->cc_node, &b->x);
		for (i = n_flip = 0; i < b->cc_size; ++i) {
			uint32_t k = b->cc_node[i];///uid
			int8_t s;
			if (b->z[k].z[0] == b->z[k].z[1]) continue;
			s = b->z[k].z[0] > b->z[k].z[1]? -1 : 1;
			if (b->s[k] != s) {
				mc_set_spin(ma, b, k, s);///no need to change the score of k itself
				++n_flip;
			}
		}
		++n_iter_local;
		if (n_flip == 0) break;
	}

	// if(n_flip != 0) mc_best_flip(ma, b);
	return mc_score(ma, b);
}

static void mc_perturb(const mc_opt_t *opt, const mc_match_t *ma, mc_svaux_t *b)
{
	uint32_t i;
	for (i = 0; i < b->cc_size; ++i) {
		uint32_t k = (uint32_t)ma->cc[b->cc_off + i];///node id
		double y;
		y = kr_drand_r(&b->x);
		if (y < opt->f_perturb)
			mc_set_spin(ma, b, k, -b->s[k]);
	}
}

static uint32_t mc_bfs(const mc_match_t *ma, mc_svaux_t *b, uint32_t k0, uint32_t bfs_round, uint32_t max_size)
{
	uint32_t i, n_bfs = 0, st, en, r;
	b->bfs[n_bfs++] = k0, b->bfs_mark[k0] = k0;
	st = 0, en = n_bfs;
	for (r = 0; r < bfs_round; ++r) {
		for (i = st; i < en; ++i) {
			uint32_t k = b->bfs[i];
			uint32_t o = ma->idx.a[k] >> 32;
			uint32_t n = (uint32_t)ma->idx.a[k], j;
			for (j = 0; j < n; ++j) {
				uint32_t t = (uint32_t)ma->ma.a[o + j].x;
				if (b->bfs_mark[t] != k0)
					b->bfs[n_bfs++] = t, b->bfs_mark[t] = k0;
			}
		}
		st = en, en = n_bfs;
		if (max_size > 0 && n_bfs > max_size) break;
	}
	return n_bfs;
}

///bfs_round is 3
static void mc_perturb_node(const mc_opt_t *opt, const mc_match_t *ma, mc_svaux_t *b, int32_t bfs_round)
{
	uint32_t i, k, n_bfs = 0;
	k = (uint32_t)(kr_drand_r(&b->x) * b->cc_size + .499);
	k = (uint32_t)ma->cc[b->cc_off + k];///node id
	n_bfs = mc_bfs(ma, b, k, bfs_round, (int32_t)(b->cc_size * opt->f_perturb));
	for (i = 0; i < n_bfs; ++i)
		mc_set_spin(ma, b, b->bfs[i], -b->s[b->bfs[i]]);
}

void clean_mc_bp_res(mc_bp_res *res)
{
	res->chain_id = (uint32_t)-1;
	res->f_bid = res->f_uid = res->l_bid = res->l_uid = (uint32_t)-1;
	res->id = (uint32_t)-1; res->w = -1;
}

void reset_mc_bp_iter(mc_bp_t* bp, mc_bp_iter *x, uint64_t id)
{
	ma_utg_t *u = NULL;
	uint32_t i, n, occ;
	for (i = 0; i < bp->idx_n; i++)
	{
		if(id >= bp->idx[i]) break;
	}

	id -= bp->idx[i];
	x->chain_id = bp->b_b->chain_weight.a[i].id;
	u = &(bp->b_b->b_ug->u.a[x->chain_id]);
	for (i = occ = 0; i < u->n; i++)
	{
		get_bubbles(bp->b_b, u->a[i]>>33, NULL, NULL, NULL, &n, NULL);
		occ += n;
		if(id < occ)
		{
			x->bid = i;
			x->uid = id - (occ -n);
			break;
		}
	}
}
inline uint32_t next_uid(mc_bp_iter *iter, bubble_type* bub, uint32_t *c_bid, uint32_t *c_uid)
{
	ma_utg_t *u = &(bub->b_ug->u.a[iter->chain_id]);
	uint32_t *a, n, uid;
	while (1)
	{
		if(iter->bid >= u->n) break;
		get_bubbles(bub, u->a[iter->bid]>>33, NULL, NULL, &a, &n, NULL);
		while (1)
		{
			if(iter->uid >= n) break;
			uid = a[iter->uid]>>1;
			if(c_bid) (*c_bid) = iter->bid;
			if(c_uid) (*c_uid) = iter->uid;
			iter->uid++;
			return uid;
		}
		iter->bid++, iter->uid = 0;
	}
	return (uint32_t)-1;
}

t_w_t incre_weight(mc_svaux_t *b_aux, mc_match_t *ma, uint8_t* vis, uint32_t uid)
{
	mc_edge_t *o = NULL;
	uint32_t n, i, t;
	t_w_t w = ((t_w_t)(b_aux->s[uid])) * (b_aux->z[uid].z[0] - b_aux->z[uid].z[1]) * 2;
	t_w_t w_off = 0;
	o = pt_a(*ma, uid);
	n = pt_n(*ma, uid);
	for (i = 0; i < n; ++i) 
	{
		t = ma_y(o[i]);
		if(vis[t] == 0) continue;
        if(t == uid) continue;
		w_off += (b_aux->s[uid]*b_aux->s[t]*o[i].w);
	}
	return w - (w_off*4);//2 for self; 4 for both directions
}

void select_min_bp(bubble_type* bub, mc_match_t *ma, mc_svaux_t *b_aux, mc_bp_t* bp,
uint8_t *lock, bits_p *vis, uint64_t id, mc_bp_res* r)
{
	uint32_t uid, val = 0, max_bid, max_uid, c_bid, c_uid, f_bid, f_uid;
	t_w_t w = 0, max_w = -1;
	mc_bp_iter i;
	reset_mc_bp_iter(bp, &i, id);
	memset(vis->a, 0, vis->n);
	max_bid = max_uid = (uint32_t)-1;
	f_bid = i.bid; f_uid = i.uid;

	while (1)
	{
		uid = next_uid(&i, bub, &c_bid, &c_uid);
		if(uid == (uint32_t)-1) break;
		if(vis->a[uid]) continue; ///already flip uid
		///update w
		w += incre_weight(b_aux, ma, vis->a, uid);
		vis->a[uid] = 1;
		if(lock[uid] == 0) val = 1;
        if(val == 0) continue;
		///update max_w
		if(max_w < w)
		{
			max_w = w;
			max_bid = c_bid;
			max_uid = c_uid;
		}
	}

	if(max_w <= 0 || max_bid == (uint32_t)-1 || max_uid == (uint32_t)-1) return;
	if((max_w > r->w) || (max_w == r->w && id < r->id))
	{
		r->w = max_w;
		r->id = id;
		r->chain_id = i.chain_id;
		r->l_bid = max_bid;
		r->l_uid = max_uid;
		r->f_bid = f_bid;
		r->f_uid = f_uid;
	}
	///if((res->min_w > i_b->weight) || (res->min_w == i_b->weight && id < res->min_idx))
}

static void worker_for_min_bp(void *data, long i, int tid) // callback for kt_for()
{
    mc_bp_t* bp = (mc_bp_t *)data;
	select_min_bp(bp->b_b, bp->ma, bp->b_aux, bp, bp->lock, &(bp->vis[tid]), i, &(bp->res[tid]));
}

uint32_t best_bp(mc_bp_t *bp, mc_bp_res *res)
{
	uint32_t i;
	clean_mc_bp_res(res);
	for (i = 0; i < bp->n_thread; i++)
	{
		clean_mc_bp_res(&(bp->res[i]));
	} 
	kt_for(bp->n_thread, worker_for_min_bp, bp, bp->occ);
	
	for (i = 0; i < bp->n_thread; i++)
	{
		if(bp->res[i].chain_id == (uint32_t)-1) continue;
		if(bp->res[i].w <= 0) continue;
		if((bp->res[i].w > res->w) || (bp->res[i].w == res->w && bp->res[i].id < res->id))
		{
			(*res) = bp->res[i];
		}
	}
	if(res->chain_id != (uint32_t)-1) return 1;
	return 0;
}

void mc_set_bp_spin(bubble_type* bub, mc_match_t *ma, mc_svaux_t *b_aux, mc_bp_t* bp,
uint8_t *lock, bits_p *vis, mc_bp_res *res)
{
	uint32_t uid, val = 0, c_bid, c_uid;
	mc_bp_iter i;
	i.chain_id = res->chain_id;
	i.bid = res->f_bid;
	i.uid = res->f_uid;
	memset(vis->a, 0, vis->n);
	while (1)
    {
        uid = next_uid(&i, bub, &c_bid, &c_uid);
        if(uid == (uint32_t)-1) break;
		if(lock[uid] == 0)
		{
			val = 1;
			break;
		}
        if(c_bid == res->l_bid && c_uid == res->l_uid) break;
    }

	if(val == 0)
	{
		fprintf(stderr, "ERROR-1\n");
		return;
	}

	i.bid = res->f_bid;
	i.uid = res->f_uid;
	while (1)
	{
		uid = next_uid(&i, bub, &c_bid, &c_uid);
        if(uid == (uint32_t)-1) break;
        if(vis->a[uid] == 1) continue;
        lock[uid] = 1;
        vis->a[uid] = 1;
		mc_set_spin(ma, b_aux, uid, -b_aux->s[uid]);
        if(c_bid == res->l_bid && c_uid == res->l_uid) break;
	}
}

double mc_solve_bp_cc(mc_bp_t *bp)
{
	mc_bp_res res;
	memset(bp->lock, 0, bp->b_b->ug->g->n_seq);
	while (best_bp(bp, &res))
	{
		mc_set_bp_spin(bp->b_b, bp->ma, bp->b_aux, bp, bp->lock, &(bp->vis[0]), &res);
	}

	return mc_score_all(bp->ma, bp->b_aux);
}

void mc_reset_z_all(const mc_match_t *ma, mc_svaux_t *b)
{
	// t_w_t z[2];
	uint32_t k;
	for (k = 0; k < ma->n_seq; ++k) 
	{
		uint32_t o = ma->idx.a[k] >> 32;
        uint32_t j, n = (uint32_t)ma->idx.a[k];
		// z[0] = b->z[k].z[0]; z[1] = b->z[k].z[1];
        b->z[k].z[0] = b->z[k].z[1] = 0;
        for (j = 0; j < n; ++j) {
            const mc_edge_t *e = &ma->ma.a[o + j];
            uint32_t t = ma_y(*e);
            if (b->s[t] > 0) b->z[k].z[0] += e->w;
            else if (b->s[t] < 0) b->z[k].z[1] += e->w;
        }
		// if(z[0] != b->z[k].z[0]) fprintf(stderr, "ERROR1-all\n");
        // if(z[1] != b->z[k].z[1]) fprintf(stderr, "ERROR2-all\n");
	}
}

void mc_solve_bp(mc_bp_t *bp)
{
	double index_time = yak_realtime();
	uint32_t r = 1;
	double sc_opt, sc;
	mc_reset_z_all(bp->ma, bp->b_aux);
	sc_opt = mc_score_all(bp->ma, bp->b_aux);

	while (1)
	{
		sc = mc_solve_bp_cc(bp);
		fprintf(stderr, "[M::%s::# round: %u] sc_opt: %f, sc: %f\n", __func__, r, sc_opt, sc);
		if(sc <= sc_opt) break;
		sc_opt = sc;
		r++;
	}
	fprintf(stderr, "[M::%s::%.3f] ==> round %u\n", __func__, yak_realtime()-index_time, r);
}

void print_sc(const mc_opt_t *opt, const mc_match_t *ma, mc_svaux_t *b, t_w_t sc_opt, uint32_t n_iter)
{
	t_w_t w = mc_score(ma, b);
	if(w != sc_opt) fprintf(stderr, "ERROR\n");
	fprintf(stderr, "# iter: %u, sc_opt: %f, sc-local: %f, sc-global: %f\n", 
					n_iter, sc_opt, w, mc_score_all(ma, b));
}

uint32_t mc_solve_cc(const mc_opt_t *opt, const mc_g_t *mg, mb_g_t *mbg, mc_svaux_t *b, uint32_t cc_off, uint32_t cc_size)
{
	uint32_t j, k, n_iter = 0;
	t_w_t sc_opt = -(1<<30), sc;///problem-w
	b->cc_off = cc_off, b->cc_size = cc_size;
	if (b->cc_size < 2) return 0;
	// print_sc(opt, mg->e, b, sc_opt, (uint32_t)-1);
	sc_opt = mc_init_spin(mg->e, b);
	if (b->cc_size == 2) return 0;
	for (j = 0; j < b->cc_size; ++j) {///backup s and z in s_opt and z_opt
		b->s_opt[b->cc_node[j]] = b->s[b->cc_node[j]]; ///hap status of each unitig
		b->z_opt[b->cc_node[j]] = b->z[b->cc_node[j]]; ///z[0]: positive weight; z[1]: positive weight
	}
	// print_sc(opt, mg->e, b, sc_opt, n_iter);
	sc = mc_optimize_local(opt, mg->e, b, &n_iter);
	if (sc > sc_opt) {
		for (j = 0; j < b->cc_size; ++j) {
			b->s_opt[b->cc_node[j]] = b->s[b->cc_node[j]];
			b->z_opt[b->cc_node[j]] = b->z[b->cc_node[j]];
		}
		sc_opt = sc;
	} else {
		for (j = 0; j < b->cc_size; ++j) {
			b->s[b->cc_node[j]] = b->s_opt[b->cc_node[j]];
			b->z[b->cc_node[j]] = b->z_opt[b->cc_node[j]];
		}
	}
	// mc_reset_z_debug(mg->e, b);
	// print_sc(opt, mg->e, b, sc_opt, n_iter);
	// fprintf(stderr, "\ncc_size: %u, cc_off: %u\n", b->cc_size, b->cc_off);
	for (k = 0; k < (uint32_t)opt->n_perturb; ++k) {
		if (k&1) mc_perturb(opt, mg->e, b);
		else mc_perturb_node(opt, mg->e, b, 3);
		sc = mc_optimize_local(opt, mg->e, b, &n_iter);
		// fprintf(stderr, "(%u) sc_opt: %f, sc: %f\n", k, sc_opt, sc);
		if (sc > sc_opt) {
			for (j = 0; j < b->cc_size; ++j) {
				b->s_opt[b->cc_node[j]] = b->s[b->cc_node[j]];
				b->z_opt[b->cc_node[j]] = b->z[b->cc_node[j]];
			}
			sc_opt = sc;
		} else {
			for (j = 0; j < b->cc_size; ++j) {
				b->s[b->cc_node[j]] = b->s_opt[b->cc_node[j]];
				b->z[b->cc_node[j]] = b->z_opt[b->cc_node[j]];
			}
		}
		// print_sc(opt, mg->e, b, sc_opt, n_iter);
	}
	for (j = 0; j < b->cc_size; ++j)
		b->s[b->cc_node[j]] = b->s_opt[b->cc_node[j]];
	return n_iter;
}



uint32_t mb_solve_cc(const mc_opt_t *opt, const mc_g_t *mg, mb_g_t *mbg, mc_svaux_t *b, uint32_t cc_off, uint32_t cc_size)
{
	uint32_t j, k, n_iter = 0;
	t_w_t sc_opt = -(1<<30), sc;///problem-w
	b->cc_off = cc_off, b->cc_size = cc_size;
	if (b->cc_size < 2) return 0;
	// print_sc(opt, mg->e, b, sc_opt, (uint32_t)-1);
	sc_opt = mc_init_spin(mg->e, b);
	if (b->cc_size == 2) return 0;
	for (j = 0; j < b->cc_size; ++j) {///backup s and z in s_opt and z_opt
		b->s_opt[b->cc_node[j]] = b->s[b->cc_node[j]]; ///hap status of each unitig
		b->z_opt[b->cc_node[j]] = b->z[b->cc_node[j]]; ///z[0]: positive weight; z[1]: positive weight
	}
	// print_sc(opt, mg->e, b, sc_opt, n_iter);
	sc = mc_optimize_local(opt, mg->e, b, &n_iter);
	if (sc > sc_opt) {
		for (j = 0; j < b->cc_size; ++j) {
			b->s_opt[b->cc_node[j]] = b->s[b->cc_node[j]];
			b->z_opt[b->cc_node[j]] = b->z[b->cc_node[j]];
		}
		sc_opt = sc;
	} else {
		for (j = 0; j < b->cc_size; ++j) {
			b->s[b->cc_node[j]] = b->s_opt[b->cc_node[j]];
			b->z[b->cc_node[j]] = b->z_opt[b->cc_node[j]];
		}
	}
	// mc_reset_z_debug(mg->e, b);
	// print_sc(opt, mg->e, b, sc_opt, n_iter);
	// fprintf(stderr, "\ncc_size: %u, cc_off: %u\n", b->cc_size, b->cc_off);
	for (k = 0; k < (uint32_t)opt->n_perturb; ++k) {
		if (k&1) mc_perturb(opt, mg->e, b);
		else mc_perturb_node(opt, mg->e, b, 3);
		sc = mc_optimize_local(opt, mg->e, b, &n_iter);
		// fprintf(stderr, "(%u) sc_opt: %f, sc: %f\n", k, sc_opt, sc);
		if (sc > sc_opt) {
			for (j = 0; j < b->cc_size; ++j) {
				b->s_opt[b->cc_node[j]] = b->s[b->cc_node[j]];
				b->z_opt[b->cc_node[j]] = b->z[b->cc_node[j]];
			}
			sc_opt = sc;
		} else {
			for (j = 0; j < b->cc_size; ++j) {
				b->s[b->cc_node[j]] = b->s_opt[b->cc_node[j]];
				b->z[b->cc_node[j]] = b->z_opt[b->cc_node[j]];
			}
		}
		// print_sc(opt, mg->e, b, sc_opt, n_iter);
	}
	for (j = 0; j < b->cc_size; ++j)
		b->s[b->cc_node[j]] = b->s_opt[b->cc_node[j]];
	return n_iter;
}

void reset_mb_g_t_z(mb_g_t *mbg)
{
	uint32_t k;
	mb_match_t *ma = mbg->e;
	for (k = 0; k < ma->n_seq; ++k) 
	{
		uint32_t o = ma->idx.a[k] >> 32;
        uint32_t j, n = (uint32_t)ma->idx.a[k];
		mbg->u->u.a[k].z[0] = mbg->u->u.a[k].z[1] = 0;
		mbg->u->u.a[k].z[2] = mbg->u->u.a[k].z[3] = 0;
		for (j = 0; j < n; ++j) {
            const mb_edge_t *e = &ma->ma.a[o + j];
            uint32_t t = ma_y(*e);
			///mbg->u->u.a[t].s[t]

			///a[0]->b[0]
			if(mbg->u->u.a[t].s[0] > 0) mbg->u->u.a[k].z[0] += e->w[0];
			else if(mbg->u->u.a[t].s[0] < 0) mbg->u->u.a[k].z[1] += e->w[0];

			///a[0]->b[1]
			if(mbg->u->u.a[t].s[1] > 0) mbg->u->u.a[k].z[0] += e->w[1];
			else if(mbg->u->u.a[t].s[1] < 0) mbg->u->u.a[k].z[1] += e->w[1];

			///a[1]->b[0]
			if(mbg->u->u.a[t].s[0] > 0) mbg->u->u.a[k].z[2] += e->w[2];
			else if(mbg->u->u.a[t].s[0] < 0) mbg->u->u.a[k].z[3] += e->w[2];

			///a[1]->b[1]
			if(mbg->u->u.a[t].s[1] > 0) mbg->u->u.a[k].z[2] += e->w[3];
			else if(mbg->u->u.a[t].s[1] < 0) mbg->u->u.a[k].z[3] += e->w[3];
        }
	}
}

void mc_init_spin_all(const mc_opt_t *opt, mc_g_t *mg, mb_g_t *mbg, mc_svaux_t *b)
{
	uint32_t st, i, k, m, a_n[2], *a[2], qn;
	for (st = 0, i = 1; i <= mg->e->n_seq; ++i) {
		if (i == mg->e->n_seq || mg->e->cc[st]>>32 != mg->e->cc[i]>>32) {
			b->cc_off = st, b->cc_size = i - st;
			if (b->cc_size >= 2)
			{
				mc_init_spin(mg->e, b);
			}
			st = i;
		}
	}

	if(!mbg) return;

	kvec_t(uint32_t) s; kv_init(s);
	uint32_t x[2], y[2], x_i, y_i, ws, we;
	uint8_t *vis = NULL; CALLOC(vis, mg->e->n_seq);
	t_w_t w = 0, max_w = 0;
	for (i = 0; i < mbg->u->u.n; i++)///each block
	{
		decode_mb_node(mbg, i, &(a[0]), &(a_n[0]), NULL, &(a[1]), &(a_n[1]), NULL);
		if(a_n[0] == 0 || a_n[1] == 0) continue;
		s.n = 0; x[0] = x[1] = y[0] = y[1] = 0;
		for (k = 0; k < a_n[0]; k++)
        {
            qn = a[0][k];
			if(b->s[qn] == 0) continue;

			x[b->s[qn] < 0]++;
			// qn <<= 1;
			if(b->s[qn] < 0) qn += ((uint32_t)1<<31);			
			kv_push(uint32_t, s, qn);
		}

		for (k = 0; k < a_n[1]; k++)
        {
            qn = a[1][k];
			if(b->s[qn] == 0) continue;

			y[b->s[qn] < 0]++;
			// qn <<= 1; qn++;
			if(b->s[qn] > 0) qn += ((uint32_t)1<<31);
			kv_push(uint32_t, s, qn);
		}

		if(x[0] + x[1] == 0) continue;
		if(y[0] + y[1] == 0) continue;
		radix_sort_mc32(s.a, s.a + s.n);

		x_i = y_i = (uint32_t)-1;
		if(x[0] == 0) x_i = 0;
		if(x[1] == 0) x_i = 1;

		if(y[0] == 0) y_i = 0;
		if(y[1] == 0) y_i = 1;

		if(x_i != (uint32_t)-1 && y_i != (uint32_t)-1 && x_i != y_i) continue;

		ws = we = (uint32_t)-1;
		for (st = 0, k = 1; k <= s.n; k++)
		{
			if(k == s.n || (s.a[k]>>31) != (s.a[st]>>31))
			{
				w = 0;
				for (m = st; m < k; m++) vis[(s.a[m]<<1)>>1] = 1;
				for (m = st; m < k; m++) w += incre_weight(b, mg->e, vis, (s.a[m]<<1)>>1);
				for (m = st; m < k; m++) vis[(s.a[m]<<1)>>1] = 0;
				if(ws == (uint32_t)-1 || max_w < w) max_w = w, ws = st, we = k;
				st = k;
			}
		}

		/*******************************for debug************************************/
		t_w_t sc_opt = mc_score_all(mg->e, b);
		/*******************************for debug************************************/
		for (m = ws; m < we; m++) mc_set_spin(mg->e, b, (s.a[m]<<1)>>1, -b->s[(s.a[m]<<1)>>1]);
		/*******************************for debug************************************/
		t_w_t sc_cur = mc_score_all(mg->e, b);
		fprintf(stderr, "sc_opt: %f, sc_cur: %f, max_w: %f\n", sc_opt, sc_cur, max_w);
		/*******************************for debug************************************/
	}

	/*******************************for debug************************************/
	for (i = 0; i < mbg->u->u.n; i++)///each block
	{
		decode_mb_node(mbg, i, &(a[0]), &(a_n[0]), NULL, &(a[1]), &(a_n[1]), NULL);
		if(a_n[0] == 0 || a_n[1] == 0) continue;
		x[0] = x[1] = y[0] = y[1] = 0;
		for (k = 0; k < a_n[0]; k++)
        {
            qn = a[0][k];
			if(b->s[qn] == 0) continue;
			x[b->s[qn] < 0]++;
		}

		for (k = 0; k < a_n[1]; k++)
        {
            qn = a[1][k];
			if(b->s[qn] == 0) continue;
			y[b->s[qn] < 0]++;
		}
		if(x[0] + x[1] == 0) continue;
		if(y[0] + y[1] == 0) continue;


		x_i = y_i = (uint32_t)-1;
		if(x[0] == 0) x_i = 0;
		if(x[1] == 0) x_i = 1;

		if(y[0] == 0) y_i = 0;
		if(y[1] == 0) y_i = 1;

		if(x_i != (uint32_t)-1 && y_i != (uint32_t)-1 && x_i != y_i) continue;

		fprintf(stderr, "ERROR-1\n");
	}
	/*******************************for debug************************************/
	kv_destroy(s); free(vis);

	for (i = 0; i < mbg->u->u.n; i++)///each block
	{
		decode_mb_node(mbg, i, &(a[0]), &(a_n[0]), NULL, &(a[1]), &(a_n[1]), NULL);
		if(a_n[0] == 0 || a_n[1] == 0) continue;
		x[0] = x[1] = y[0] = y[1] = 0;
		mbg->u->u.a[i].s[0] = mbg->u->u.a[i].s[1] = 0;
		mbg->u->u.a[i].z[0] = mbg->u->u.a[i].z[1] = 0;
		mbg->u->u.a[i].z[2] = mbg->u->u.a[i].z[3] = 0;
		for (k = 0; k < a_n[0]; k++)
        {
            qn = a[0][k];
			if(b->s[qn] == 0) continue;
			mbg->u->u.a[i].s[0] = b->s[qn];
			break;
		}

		for (k = 0; k < a_n[1]; k++)
        {
            qn = a[1][k];
			if(b->s[qn] == 0) continue;
			mbg->u->u.a[i].s[1] = b->s[qn];
			break;
		}
	}

	reset_mb_g_t_z(mbg);
}


uint64_t *mb_g_cc_core(mb_match_t *ma)
{
    uint32_t i, x, y, *flag;
    uint64_t *group;
    mb_edge_t *o = NULL;
    kvec_t(uint32_t) stack; kv_init(stack);

    MALLOC(flag, ma->n_seq);
    for (i = 0; i < ma->n_seq; ++i)
        flag[i] = (uint32_t)-1;

    // connected componets
    for (i = 0; i < ma->n_seq; ++i) {
        if (flag[i] != (uint32_t)-1) continue;
        stack.n = 0;
        kv_push(uint32_t, stack, i);
        while (stack.n > 0) {
            uint32_t k, j, n;
            stack.n--;
            k = stack.a[stack.n];
            flag[k] = i;///group id
            // n = (uint32_t)ma->idx[k];
            // s = ma->idx[k] >> 32;
            o = pt_a(*ma, k); 
            n = pt_n(*ma, k);
            for (j = 0; j < n; ++j) {
                uint32_t t = ma_y(o[j]);
                if (flag[t] != (uint32_t)-1) continue;
                // if (ns == ms) PT_EXPAND(stack, ms);
                // stack[ns++] = t;
                kv_push(uint32_t, stack, t);
            }
        }
    }
    kv_destroy(stack);

    // precalculate the size of each group
    CALLOC(group, ma->n_seq);
    for (i = 0; i < ma->n_seq; ++i)
        group[i] = (uint64_t)flag[i] << 32 | i;
    radix_sort_mc64(group, group + ma->n_seq);
    for (i = 1, x = y = 0; i <= ma->n_seq; ++i) {
        if (i == ma->n_seq || group[i]>>32 != group[x]>>32) {
            uint32_t j;
            for (j = x; j < i; ++j)
                group[j] = (uint64_t)y << 32 | (uint32_t)group[j];///(group id)|first element in this group
            ++y, x = i;
        }
    }
    free(flag);
    return group;
}

void mb_g_cc(mb_g_t *mbg)
{
	mbg->e->cc = mb_g_cc_core(mbg->e);
}

void mc_solve_core(const mc_opt_t *opt, mc_g_t *mg, bubble_type* bub, mb_g_t *mbg)
{
	double index_time = yak_realtime();
	uint32_t st, i;
	mc_svaux_t *b;
	mc_bp_t *bp = NULL;
	mc_g_cc(mg->e);
	b = mc_svaux_init(mg, opt->seed);
	if(mbg) mb_g_cc(mbg);
	if(bub) bp = mc_bp_t_init(mg->e, b, bub, asm_opt.thread_num);
	/*******************************for debug************************************/
	if(mbg || bp) mc_init_spin_all(opt, mg, mbg, b);
	if(bp) mc_solve_bp(bp);
	/*******************************for debug************************************/
	// fprintf(stderr, "\n\n\n\n\n*************beg-[M::%s::score->%f] ==> Partition\n", __func__, mc_score_all(mg->e, b));
	for (st = 0, i = 1; i <= mg->e->n_seq; ++i) {
		if (i == mg->e->n_seq || mg->e->cc[st]>>32 != mg->e->cc[i]>>32) {
			mc_solve_cc(opt, mg, mbg, b, st, i - st);
			st = i;
		}
	}
	// fprintf(stderr, "##############end-[M::%s::score->%f] ==> Partition\n", __func__, mc_score_all(mg->e, b));
	if(bp) mc_solve_bp(bp);	
	///mc_write_info(g, b);
	mc_svaux_destroy(b);
	if(bp) destroy_mc_bp_t(&bp);
	fprintf(stderr, "[M::%s::%.3f] ==> Partition\n", __func__, yak_realtime()-index_time);
}

void set_p_flag(mc_g_t *mg, uint32_t uID, uint8_t* trio_flag, trans_chain* t_ch, int8_t s)
{
	uint32_t i;
	ma_utg_t *u = &(mg->ug->u.a[uID]);
	for (i = 0; i < u->n; i++)
	{
		trio_flag[u->a[i]>>33] |= SET_TRIO;
		if(t_ch->is_r_het[u->a[i]>>33] == N_HET) continue;
		if(s == 0)
		{
			if(t_ch->is_r_het[u->a[i]>>33]&P_HET)//special case
			{
				trio_flag[u->a[i]>>33] |= FATHER;
			}
			continue;
		} 
		trio_flag[u->a[i]>>33] |= (s > 0? FATHER:MOTHER);
	}
}

void filter_ovlp_by_mc(ma_ug_t *ug, asg_t *read_g, uint32_t uID, hap_overlaps_list* ha, mc_match_t *ma, 
int8_t *s)
{
	mc_edge_t *o = pt_a(*ma, uID);
	uint32_t n = pt_n(*ma, uID), k, qn, tn;
	hap_overlaps *p = NULL;
	int index;

	for (k = 0; k < n; ++k)
	{
		qn = ma_x(o[k]); tn = ma_y(o[k]); p = NULL;
		if((s[qn]*s[tn])!=-1) continue;
		index = get_specific_hap_overlap(&(ha->x[qn]), qn, tn);
		if(index != -1 && ha->x[qn].a.a[index].score == (long long)o[k].w)
		{
			p = &(ha->x[qn].a.a[index]);
		}
		else
		{
			index = get_specific_hap_overlap(&(ha->x[tn]), tn, qn);
			if(index != -1 && ha->x[tn].a.a[index].score == (long long)o[k].w)
			{
				p = &(ha->x[tn].a.a[index]);
			}
		}
		if(!p) fprintf(stderr, "ERROR\n");
		p->status = FLIP;
	}
}

void clean_ovlp_by_mc(mc_g_t *mg, hap_overlaps_list* ha)
{
	uint32_t v, i, k, qn, tn, types[4];
    types[X2Y] = Y2X; types[Y2X] = X2Y; types[XCY] = YCX; types[YCX] = XCY;
	int index;
	hap_overlaps *x = NULL, *y = NULL;
	for (i = 0; i < mg->e->n_seq; ++i) 
	{
		filter_ovlp_by_mc(mg->ug, mg->rg, i, ha, mg->e, mg->s.a);
	}

	for (v = 0; v < ha->num; v++)
	{
		for (i = 0; i < ha->x[v].a.n; i++)
		{
			qn = ha->x[v].a.a[i].xUid;
			tn = ha->x[v].a.a[i].yUid;
			x = &(ha->x[v].a.a[i]);
			if(x->status != FLIP) continue;
			index = get_specific_hap_overlap(&(ha->x[tn]), tn, qn);
			if(index != -1)
			{
				y = &(ha->x[tn].a.a[index]);
				set_reverse_hap_overlap(y, x, types);
				y->status = FLIP;
			}
			if(index == -1) fprintf(stderr, "ERROR\n");
		}
	}

	for (v = 0; v < ha->num; v++)
    {
        for (i = k = 0; i < ha->x[v].a.n; i++)
        {
            if(ha->x[v].a.a[i].status != FLIP) continue;

            ha->x[v].a.a[k] = ha->x[v].a.a[i];
            k++;
        }
        ha->x[v].a.n = k;
    }
}


void p_nodes(mc_g_t *mg, trans_chain* t_ch, uint8_t* trio_flag)
{
	uint32_t i;
	for (i = 0; i < mg->e->n_seq; ++i) 
	{
		set_p_flag(mg, i, trio_flag, t_ch, mg->s.a[i]);
	}
}

void mc_solve(hap_overlaps_list* ovlp, trans_chain* t_ch, kv_u_trans_t *ta, ma_ug_t *ug, asg_t *read_g, double f_rate, uint8_t* trio_flag, uint32_t renew_s, int8_t *s, uint32_t is_sys, bubble_type* bub, mb_nodes_t* u)
{
	mc_opt_t opt;
	mc_opt_init(&opt, asm_opt.n_perturb, asm_opt.f_perturb, asm_opt.seed);
	mc_g_t *mg = init_mc_g_t(ug, read_g, s, renew_s);
	update_mc_edges(mg, ovlp, ta, t_ch, f_rate, is_sys);
	mb_g_t *mbg = init_mb_g_t(mg, u, is_sys);
	///debug_mc_g_t(mg);
	mc_solve_core(&opt, mg, bub, mbg);

	if((asm_opt.flag & HA_F_PARTITION) && t_ch)
	{
		p_nodes(mg, t_ch, trio_flag);
	}

	if(ovlp) clean_ovlp_by_mc(mg, ovlp);

	destory_mc_g_t(&mg);
	destory_mb_g_t(&mbg);
}