#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <stdlib.h>
#include "rcut.h"
#include "Purge_Dups.h"
#include "Correct.h"
#include "ksort.h"
#include "kthread.h"
#include "hic.h"

#define VERBOSE_CUT 0

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

#define mcb_pat(x, id) (&((x).m.a[(id)-1]))
#define mcp_de(x, id, m) ((x).z[((id)<<(x).hapN)+(m)])

uint8_t bit_filed[8] = {1, 2, 4, 8, 16, 32, 64, 128};
#define is_bit_set(id, a) ((a)[(id)>>3]&bit_filed[(id)&7])

typedef struct {
	int32_t max_iter;
    int32_t n_perturb;
	int32_t n_b_perturb;
	int32_t n_s_perturb;
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


typedef struct {
	uint64_t x; // RNG
	uint32_t cc_off, cc_size;
	kvec_t(uint64_t) cc_edge;
	uint32_t *cc_node;
	uint32_t *bfs, *bfs_mark;
	mb_node_t *u, *u_opt;///keep status
	uint8_t *f;
} mb_svaux_t;

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
	opt->n_perturb = n_perturb;
	opt->f_perturb = f_perturb;
	opt->max_iter = 1000;
	opt->seed = seed;
	opt->n_s_perturb = n_perturb;
	opt->n_b_perturb = n_perturb*0.5;
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
			fprintf(stderr, "ERROR-0-::%s\n", __func__);
			continue;
		}
		t = get_mb_edge(ma, ma_y(*m), ma_x(*m));
		if(!t)
		{
			fprintf(stderr, "ERROR-1-::%s\n", __func__);
			continue;
		} 

		if(m->w[0] != t->w[0])
		{
			fprintf(stderr, "\nERROR-2-::%s\n", __func__);
			fprintf(stderr, "ma_x(*m): %lu, ma_y(*m): %u\n", ma_x(*m), ma_y(*m));
			// fprintf(stderr, "m->w[0]: %ld, m->w[1]: %ld, m->w[2]: %ld, m->w[3]: %ld\n", 
			// 										m->w[0], m->w[1], m->w[2], m->w[3]);
			// fprintf(stderr, "t->w[0]: %ld, t->w[1]: %ld, t->w[2]: %ld, t->w[3]: %ld\n", 
			// 										t->w[0], t->w[1], t->w[2], t->w[3]);
			fprintf(stderr, "m->w[0]: %f, m->w[1]: %f, m->w[2]: %f, m->w[3]: %f\n", 
													m->w[0], m->w[1], m->w[2], m->w[3]);
			fprintf(stderr, "t->w[0]: %f, t->w[1]: %f, t->w[2]: %f, t->w[3]: %f\n", 
													t->w[0], t->w[1], t->w[2], t->w[3]);
		} 
		
		if(m->w[3] != t->w[3])
		{
			fprintf(stderr, "\nERROR-3-::%s\n", __func__);
			fprintf(stderr, "ma_x(*m): %lu, ma_y(*m): %u\n", ma_x(*m), ma_y(*m));
			// fprintf(stderr, "m->w[0]: %ld, m->w[1]: %ld, m->w[2]: %ld, m->w[3]: %ld\n", 
			// 										m->w[0], m->w[1], m->w[2], m->w[3]);
			// fprintf(stderr, "t->w[0]: %ld, t->w[1]: %ld, t->w[2]: %ld, t->w[3]: %ld\n", 
			// 										t->w[0], t->w[1], t->w[2], t->w[3]);
			fprintf(stderr, "m->w[0]: %f, m->w[1]: %f, m->w[2]: %f, m->w[3]: %f\n", 
													m->w[0], m->w[1], m->w[2], m->w[3]);
			fprintf(stderr, "t->w[0]: %f, t->w[1]: %f, t->w[2]: %f, t->w[3]: %f\n", 
													t->w[0], t->w[1], t->w[2], t->w[3]);
		} 
		
		if(m->w[1] != t->w[2])
		{
			fprintf(stderr, "\nERROR-4-::%s\n", __func__);
			fprintf(stderr, "ma_x(*m): %lu, ma_y(*m): %u\n", ma_x(*m), ma_y(*m));
			// fprintf(stderr, "m->w[0]: %ld, m->w[1]: %ld, m->w[2]: %ld, m->w[3]: %ld\n", 
			// 										m->w[0], m->w[1], m->w[2], m->w[3]);
			// fprintf(stderr, "t->w[0]: %ld, t->w[1]: %ld, t->w[2]: %ld, t->w[3]: %ld\n", 
			// 										t->w[0], t->w[1], t->w[2], t->w[3]);
			fprintf(stderr, "m->w[0]: %f, m->w[1]: %f, m->w[2]: %f, m->w[3]: %f\n", 
													m->w[0], m->w[1], m->w[2], m->w[3]);
			fprintf(stderr, "t->w[0]: %f, t->w[1]: %f, t->w[2]: %f, t->w[3]: %f\n", 
													t->w[0], t->w[1], t->w[2], t->w[3]);
		} 
		
		
		
		if(m->w[2] != t->w[1])
		{
			fprintf(stderr, "\nERROR-5-::%s\n", __func__);
			fprintf(stderr, "ma_x(*m): %lu, ma_y(*m): %u\n", ma_x(*m), ma_y(*m));
			// fprintf(stderr, "m->w[0]: %ld, m->w[1]: %ld, m->w[2]: %ld, m->w[3]: %ld\n", 
			// 										m->w[0], m->w[1], m->w[2], m->w[3]);
			// fprintf(stderr, "t->w[0]: %ld, t->w[1]: %ld, t->w[2]: %ld, t->w[3]: %ld\n", 
			// 										t->w[0], t->w[1], t->w[2], t->w[3]);
			fprintf(stderr, "m->w[0]: %f, m->w[1]: %f, m->w[2]: %f, m->w[3]: %f\n", 
													m->w[0], m->w[1], m->w[2], m->w[3]);
			fprintf(stderr, "t->w[0]: %f, t->w[1]: %f, t->w[2]: %f, t->w[3]: %f\n", 
													t->w[0], t->w[1], t->w[2], t->w[3]);
		} 
		
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


uint64_t *mb_nodes_core(kv_u_trans_t *ref, mc_match_t* ma, uint32_t occ, uint32_t *fg, kvec_t_u32_warp *sk)
{
	uint32_t i, x, y;
	uint64_t *group;
	u_trans_t *o = NULL;

	for (i = 0; i < occ; ++i)
		fg[i] = (uint32_t)-1;

	// connected componets
	for (i = 0; i < occ; ++i) {
		if (fg[i] != (uint32_t)-1) continue;
		if (pt_n(*ma, i) == 0) 
		{
			fg[i] = i;///group id
			continue;
		}
		sk->a.n = 0;
		kv_push(uint32_t, sk->a, i);
		while (sk->a.n > 0) {
			uint32_t k, j, st, n, t;
			sk->a.n--;
			k = sk->a.a[sk->a.n];
			fg[k] = i;///group id
			o = u_trans_a(*ref, k); 
			n = u_trans_n(*ref, k);
			for (st = 0, j = 1; j <= n; ++j)
			{
				if(j == n || o[j].tn != o[st].tn)
				{
					t = o[st].tn;
					if ((pt_n(*ma, t) != 0) && (fg[t] == (uint32_t)-1))
					{
						kv_push(uint32_t, sk->a, t);
					}
					st = j;
				}
			}
		}
	}

	// precalculate the size of each group
	CALLOC(group, occ);
	for (i = 0; i < occ; ++i)
		group[i] = (uint64_t)fg[i] << 32 | i;
	radix_sort_mc64(group, group + occ);
	for (i = 1, x = y = 0; i <= occ; ++i) {
		if (i == occ || group[i]>>32 != group[x]>>32) {
			uint32_t j;
			for (j = x; j < i; ++j)
				group[j] = (uint64_t)y << 32 | (uint32_t)group[j];///(group id)|first element in this group
			++y, x = i;
		}
	}
	return group;
}

void assgin_mb_node(mb_nodes_t *x, kv_u_trans_t *ref, mc_match_t* ma, uint32_t *fg, kvec_t_u32_warp *sk, uint64_t *cc, uint32_t cc_off, uint32_t cc_size)
{
	uint32_t i, v, n, st, j, t, pass, uid;
	u_trans_t *o = NULL;
	mb_node_t *p = NULL;
	for (i = 0; i < cc_size; ++i)
	{
		v = (uint32_t)cc[cc_off + i];///node id
		fg[v] = (uint32_t)-1;
	}

	sk->a.n = 0;
	v = (uint32_t)cc[cc_off];
	kv_push(uint32_t, sk->a, v);
	fg[v] = 0;
	pass = 1;
	while (sk->a.n > 0)
	{
		sk->a.n--;
		v = sk->a.a[sk->a.n];
		if(pt_n(*ma, v) == 0) continue;
		o = u_trans_a(*ref, v); 
		n = u_trans_n(*ref, v);
		for (st = 0, j = 1; j <= n; ++j)
		{
			if(j == n || o[j].tn != o[st].tn)
			{
				t = o[st].tn;
				if (pt_n(*ma, t) == 0)
				{
					st = j;
					continue;
				} 
				
				if (fg[t] == (uint32_t)-1)///uncolor
				{
					fg[t] = 1 - fg[v];
					kv_push(uint32_t, sk->a, t);
				}
				else if(fg[t] == fg[v]) ///color
				{
					pass = 0;
					break;
				} 
				st = j;
			}
		}
		if(pass == 0) break;
	}
	
	if(pass)
	{
		uid = x->u.n;
		kv_pushp(mb_node_t, x->u, &p);
		p->z[0] = p->z[1] = p->z[2] = p->z[3] = 0;
		p->s[0] = p->s[1] = 0;
		p->occ[0] = p->occ[1] = 0;

		for (i = 0, p->a[0] = x->bid.n; i < cc_size; ++i)
		{
			v = (uint32_t)cc[cc_off + i];///node id
			if(fg[v] == 0)
			{
				kv_push(uint32_t, x->bid, v);
				p->occ[0]++;
				x->idx.a[v] = (uid<<1);
			}
		}

		for (i = 0, p->a[1] = x->bid.n; i < cc_size; ++i)
		{
			v = (uint32_t)cc[cc_off + i];///node id
			if(fg[v] == 1)
			{
				kv_push(uint32_t, x->bid, v);
				p->occ[1]++;
				x->idx.a[v] = (uid<<1) + 1;
			}
		}
	}
	else
	{
		for (i = 0; i < cc_size; ++i)
		{
			v = (uint32_t)cc[cc_off + i];///node id
			uid = x->u.n;
			kv_pushp(mb_node_t, x->u, &p);
			p->z[0] = p->z[1] = p->z[2] = p->z[3] = 0;
			p->s[0] = p->s[1] = 0;
			p->occ[0] = p->occ[1] = 0;

			p->a[0] = x->bid.n;
			kv_push(uint32_t, x->bid, v);
			p->occ[0]++;
			x->idx.a[v] = (uid<<1);

			p->a[1] = x->bid.n;
		}
	}
}

void debug_mb_nodes(mb_nodes_t *x, kv_u_trans_t *ref, mc_match_t* ma)
{
	uint32_t i, k, bid, ori, *a = NULL, a_n, v;
	for (i = 0; i < x->idx.n; i++)
	{
		if(x->idx.a[i] == (uint32_t)-1)
		{
			fprintf(stderr, "ERROR-0-::%s\n", __func__);
			continue;
		}
		bid = x->idx.a[i]>>1; ori = x->idx.a[i] & 1;
		a = x->bid.a + x->u.a[bid].a[ori];
		a_n = x->u.a[bid].occ[ori];
		for (k = 0; k < a_n; k++)
		{
			if(a[k] == i) break;
		}
		if(k >= a_n) fprintf(stderr, "ERROR-1-::%s\n", __func__);
	}
	
	uint32_t *tt[2], t_n[2], o_n, m, t, found;
	int8_t *vis = NULL; CALLOC(vis, x->idx.n);
	u_trans_t *o = NULL;
	for (i = 0; i < x->u.n; i++)///each block
	{
		tt[0] = x->bid.a + x->u.a[i].a[0];
		t_n[0] = x->u.a[i].occ[0];
		tt[1] = x->bid.a + x->u.a[i].a[1];
		t_n[1] = x->u.a[i].occ[1];

		if(t_n[0] == 1 && t_n[1] == 0) continue;

		for (k = 0; k < t_n[0]; k++)
		{
			vis[tt[0][k]] = 1;
		}
		for (k = 0; k < t_n[1]; k++)
		{
			vis[tt[1][k]] = -1;
		}





		for (k = 0; k < t_n[0]; k++)
		{
			v = tt[0][k];
			o = u_trans_a(*ref, v); 
        	o_n = u_trans_n(*ref, v);
			found = 0;
			for (m = 0; m < o_n; m++)
			{
				t = o[m].tn;
				if (pt_n(*ma, t) == 0) continue;
				if (vis[t] != -1) fprintf(stderr, "ERROR-2-::%s\n", __func__);
				else found = 1;
			}
			if(found == 0) fprintf(stderr, "ERROR-2-*::%s\n", __func__);
		}

		for (k = 0; k < t_n[1]; k++)
		{
			v = tt[1][k];
			o = u_trans_a(*ref, v); 
        	o_n = u_trans_n(*ref, v);
			found = 0;
			for (m = 0; m < o_n; m++)
			{
				t = o[m].tn;
				if (pt_n(*ma, t) == 0) continue;
				if (vis[t] != 1) fprintf(stderr, "ERROR-3-::%s\n", __func__);
				else found = 1;
			}
			if(found == 0) fprintf(stderr, "ERROR-3-*::%s\n", __func__);
		}


		for (k = 0; k < t_n[0]; k++)
		{
			vis[tt[0][k]] = 0;
		}
		for (k = 0; k < t_n[1]; k++)
		{
			vis[tt[1][k]] = 0;
		}
	}

	free(vis);
}

mb_nodes_t *update_mb_nodes_t(kv_u_trans_t *ref, mc_match_t* ma, uint32_t occ)
{
	uint32_t i, st;
	uint64_t *cc = NULL;
	mb_nodes_t *x = NULL; CALLOC(x, 1);
	x->bid.n = x->u.n = x->idx.n = 0;
	kv_resize(uint32_t, x->idx, occ);
	x->idx.n = occ;
	memset(x->idx.a, -1, sizeof(uint32_t)*x->idx.n);
	kvec_t_u32_warp stack; kv_init(stack.a);
	uint32_t *flag = NULL; MALLOC(flag, occ);

	cc = mb_nodes_core(ref, ma, occ, flag, &stack);

	for (st = 0, i = 1; i <= occ; ++i) {
        if (i == occ || cc[st]>>32 != cc[i]>>32) {
			assgin_mb_node(x, ref, ma, flag, &stack, cc, st, i - st);
            st = i;
        }
	}
	
	free(cc); free(flag);
	kv_destroy(stack.a);

	/*******************************for debug************************************/
	// debug_mb_nodes(x, ref, ma);
	/*******************************for debug************************************/
	return x;
}

mb_g_t *init_mb_g_t(mc_match_t* e, kv_u_trans_t *ref, uint32_t is_sys)
{
	mc_edge_t *o = NULL;
	uint32_t i, k, m, n, a_n[2], *a[2], qn, tn, qb, tb;
	mb_edge_t *t = NULL;
	mb_g_t *p = NULL; CALLOC(p, 1);
	p->u = update_mb_nodes_t(ref, e, e->n_seq);
	p->e = NULL; CALLOC(p->e, 1);
	p->e->n_seq = p->u->u.n;
	kv_init(p->e->ma); kv_init(p->e->idx);
	for (i = 0; i < p->u->u.n; i++)///each block
	{
		qb = i;
		p->u->u.a[i].s[0] = p->u->u.a[i].s[1] = 0;
		decode_mb_node(p, i, &(a[0]), &(a_n[0]), NULL, &(a[1]), &(a_n[1]), NULL);
		for (k = 0; k < a_n[0]; k++)
		{
			qn = a[0][k];
			o = pt_a(*e, qn);
    		n = pt_n(*e, qn);
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
			o = pt_a(*e, qn);
    		n = pt_n(*e, qn);
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
	// debug_mb_edges_symm(p->e);
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

	kv_destroy((*p)->u->bid);
	kv_destroy((*p)->u->idx);
	kv_destroy((*p)->u->u);
	free((*p)->u);

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
                if(t_ch->ir_het[u->a[i]>>33] != t->hs)
                {
                    fprintf(stderr, "ERROR-d: is_r_het: %u, h_status: %u\n", t_ch->ir_het[u->a[i]>>33], t->hs);
                }
            }
        }
    }
}

double get_w_scale(kv_u_trans_t *ta)
{
	uint32_t i, max_w_i, min_w_i;
	double max_w = 0, min_w = 0, w;
	max_w_i = min_w_i = (uint32_t)-1;
	for (i = 0; i < ta->n; ++i)
	{
		if(ta->a[i].del) continue;
		if(ta->a[i].nw == 0) continue;
		w = (ta->a[i].nw >= 0? ta->a[i].nw:-ta->a[i].nw);
		if(max_w_i == (uint32_t)-1 || max_w < w) max_w = w, max_w_i = i;
		if(min_w_i == (uint32_t)-1 || min_w > w) min_w = w, min_w_i = i;
	}

	if(max_w_i == (uint32_t)-1 || min_w_i == (uint32_t)-1) return 1;
	if(min_w > 1.1) return 1;
	double sc_max = (double)(1<<30), sc_min = 1.1;

	return MIN(sc_max/max_w, sc_min/min_w) + sc_min;
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
				if (k == u->n || t_ch->ir_het[u->a[k]>>33] != t_ch->ir_het[u->a[l]>>33])
				{
					kv_pushp(mc_interval_t, p, &t);
					t->uID = v;
					t->hs = t_ch->ir_het[u->a[l]>>33];

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
		// double sc = get_w_scale(ta);
		// fprintf(stderr, "sc: %f\n", sc);
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
			// ma->w = w_cast((ta->a[i].nw*sc));
			ma->w = w_cast((ta->a[i].nw));
		}
	}

	for (i = k = 0; i < mg->e->ma.n; i++)
    {
		if(mg->e->ma.a[i].w == 0) continue;
		mg->e->ma.a[k] = mg->e->ma.a[i];
        k++;
    }
    mg->e->ma.n = k;
	
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


mb_svaux_t *mb_svaux_init(const mb_g_t *mg, uint64_t x)
{
	uint32_t st, i, max_cc = 0;
	mb_match_t *ma = mg->e;
	mb_svaux_t *b;
	CALLOC(b, 1);
	b->x = x;
	for (st = 0, i = 1; i <= ma->n_seq; ++i)
		if (i == ma->n_seq || ma->cc[st]>>32 != ma->cc[i]>>32)
			max_cc = max_cc > i - st? max_cc : i - st, st = i;
	kv_init(b->cc_edge);
	MALLOC(b->cc_node, max_cc);

	b->u = mg->u->u.a;
	CALLOC(b->u_opt, ma->n_seq);

	MALLOC(b->bfs, ma->n_seq);
	MALLOC(b->bfs_mark, ma->n_seq); 
	memset(b->bfs_mark, -1, ma->n_seq*sizeof(uint32_t));
	
	CALLOC(b->f, ma->n_seq);
	return b;
}

void mb_svaux_destroy(mb_svaux_t *b)
{
	b->u = NULL;
	kv_destroy(b->cc_edge); free(b->cc_node);
	free(b->u); free(b->u_opt);
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

t_w_t mb_score(mb_g_t *mbg, mb_svaux_t *b)
{
	uint32_t i;
	t_w_t z = 0;
	mb_match_t *ma = mbg->e;
	for (i = 0; i < b->cc_size; ++i) {
		uint32_t k = (uint32_t)ma->cc[b->cc_off + i];///uid
		///a[0]
		z += -(t_w_t)(mbg->u->u.a[k].s[0]) * (mbg->u->u.a[k].z[0] - mbg->u->u.a[k].z[1]);
		///a[1]
		z += -(t_w_t)(mbg->u->u.a[k].s[1]) * (mbg->u->u.a[k].z[2] - mbg->u->u.a[k].z[3]);
	}
	return z;
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

void debug_mb_z(mb_g_t *mbg, uint32_t k)
{
	t_w_t z[4];
	z[0] = z[1] = z[2] = z[3] = 0;
	mb_match_t *ma = mbg->e;
	uint32_t o = ma->idx.a[k] >> 32;
    uint32_t j, n = (uint32_t)ma->idx.a[k];

	for (j = 0; j < n; ++j) {
		const mb_edge_t *e = &ma->ma.a[o + j];
		uint32_t t = ma_y(*e);

		///a[0]->b[0]
		if(mbg->u->u.a[t].s[0] > 0) z[0] += e->w[0];
		else if(mbg->u->u.a[t].s[0] < 0) z[1] += e->w[0];

		///a[0]->b[1]
		if(mbg->u->u.a[t].s[1] > 0) z[0] += e->w[1];
		else if(mbg->u->u.a[t].s[1] < 0) z[1] += e->w[1];

		///a[1]->b[0]
		if(mbg->u->u.a[t].s[0] > 0) z[2] += e->w[2];
		else if(mbg->u->u.a[t].s[0] < 0) z[3] += e->w[2];

		///a[1]->b[1]
		if(mbg->u->u.a[t].s[1] > 0) z[2] += e->w[3];
		else if(mbg->u->u.a[t].s[1] < 0) z[3] += e->w[3];
	}

	if(mbg->u->u.a[k].z[0] != z[0]) fprintf(stderr, "ERROR-0-::%s, z[0]: %f, n-z[0]: %f\n", __func__, z[0], mbg->u->u.a[k].z[0]);
	if(mbg->u->u.a[k].z[1] != z[1]) fprintf(stderr, "ERROR-1-::%s, z[1]: %f, n-z[1]: %f\n", __func__, z[1], mbg->u->u.a[k].z[1]);
	if(mbg->u->u.a[k].z[2] != z[2]) fprintf(stderr, "ERROR-2-::%s, z[2]: %f, n-z[2]: %f\n", __func__, z[2], mbg->u->u.a[k].z[2]);
	if(mbg->u->u.a[k].z[3] != z[3]) fprintf(stderr, "ERROR-3-::%s, z[3]: %f, n-z[3]: %f\n", __func__, z[3], mbg->u->u.a[k].z[3]);
}

///k is uid
static void mb_flip_spin(mb_g_t *mbg, mb_svaux_t *b, uint32_t k)
{
	if(mbg->u->u.a[k].s[0] == 0 && mbg->u->u.a[k].s[1] == 0) return;
	mb_match_t *ma = mbg->e;
	uint32_t o, j, n;
	o = ma->idx.a[k] >> 32;
	n = (uint32_t)ma->idx.a[k];
	mbg->u->u.a[k].s[0] = -mbg->u->u.a[k].s[0];
	mbg->u->u.a[k].s[1] = -mbg->u->u.a[k].s[1];
	for (j = 0; j < n; ++j) {
		const mb_edge_t *e = &ma->ma.a[o + j];
		uint32_t t = ma_y(*e);///1->z[0]; (-1)->z[1];
		
		if(mbg->u->u.a[k].s[0] != 0)
		{
			///note: e is from k to t
			///t[0] ---> k[0], so update z[0 + x] and e[0]
			mbg->u->u.a[t].z[mbg->u->u.a[k].s[0] < 0] += e->w[0];
			mbg->u->u.a[t].z[mbg->u->u.a[k].s[0] > 0] -= e->w[0];

			///t[1] ---> k[0], so update z[2 + x] and e[1]
			mbg->u->u.a[t].z[(mbg->u->u.a[k].s[0] < 0) + 2] += e->w[1];
			mbg->u->u.a[t].z[(mbg->u->u.a[k].s[0] > 0) + 2] -= e->w[1];
		}

		if(mbg->u->u.a[k].s[1] != 0)
		{
			///note: e is from k to t
			///t[0] ---> k[1], so update z[0 + x] and e[2]
			mbg->u->u.a[t].z[mbg->u->u.a[k].s[1] < 0] += e->w[2];
			mbg->u->u.a[t].z[mbg->u->u.a[k].s[1] > 0] -= e->w[2];

			///t[1] ---> k[1], so update z[2 + x] and e[3]
			mbg->u->u.a[t].z[(mbg->u->u.a[k].s[1] < 0) + 2] += e->w[3];
			mbg->u->u.a[t].z[(mbg->u->u.a[k].s[1] > 0) + 2] -= e->w[3];
		}

		/*******************************for debug************************************/
		// debug_mb_z(mbg, t);
		/*******************************for debug************************************/
	}
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

static t_w_t mb_optimize_local(const mc_opt_t *opt, mb_g_t *mbg, mb_svaux_t *b, uint32_t *n_iter)
{
	uint32_t i, n_flip = 0;
	int32_t n_iter_local = 0;
	mb_node_t *u = NULL;
	t_w_t z = 0;
	while (n_iter_local < opt->max_iter) {
		++(*n_iter);
		ks_shuffle_uint32_t(b->cc_size, b->cc_node, &b->x);
		for (i = n_flip = 0; i < b->cc_size; ++i) {
			uint32_t k = b->cc_node[i];///uid
			u = &(mbg->u->u.a[k]);
			if(u->z[0] == u->z[1] && u->z[2] == u->z[3]) continue;
			z = 0;
			///a[0]
			z += -(t_w_t)(u->s[0]) * (u->z[0] - u->z[1]);
			///a[1]
			z += -(t_w_t)(u->s[1]) * (u->z[2] - u->z[3]);
			if(z >= 0) continue;

			mb_flip_spin(mbg, b, k);///no need to change the score of k itself
			++n_flip;
		}
		++n_iter_local;
		if (n_flip == 0) break;
	}

	// if(n_flip != 0) mc_best_flip(ma, b);
	return mb_score(mbg, b);
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

static void mb_perturb(const mc_opt_t *opt, mb_g_t *mbg, mb_svaux_t *b)
{
	uint32_t i;
	for (i = 0; i < b->cc_size; ++i) {
		uint32_t k = (uint32_t)mbg->e->cc[b->cc_off + i];///node id
		double y;
		y = kr_drand_r(&b->x);
		if (y < opt->f_perturb)
			mb_flip_spin(mbg, b, k);
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
	if(k >= b->cc_size) k = b->cc_size - 1;
	k = (uint32_t)ma->cc[b->cc_off + k];///node id
	n_bfs = mc_bfs(ma, b, k, bfs_round, (int32_t)(b->cc_size * opt->f_perturb));
	for (i = 0; i < n_bfs; ++i)
		mc_set_spin(ma, b, b->bfs[i], -b->s[b->bfs[i]]);
}


static uint32_t mb_bfs(const mb_match_t *ma, mb_svaux_t *b, uint32_t k0, uint32_t bfs_round, uint32_t max_size)
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
static void mb_perturb_node(const mc_opt_t *opt, mb_g_t *mbg, mb_svaux_t *b, int32_t bfs_round)
{
	uint32_t i, k, n_bfs = 0;
	k = (uint32_t)(kr_drand_r(&b->x) * b->cc_size + .499);
	if(k >= b->cc_size) k = b->cc_size - 1;
	k = (uint32_t)mbg->e->cc[b->cc_off + k];///node id
	n_bfs = mb_bfs(mbg->e, b, k, bfs_round, (int32_t)(b->cc_size * opt->f_perturb));
	for (i = 0; i < n_bfs; ++i)
		mb_flip_spin(mbg, b, b->bfs[i]);
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

t_w_t mc_solve_bp_cc(mc_bp_t *bp)
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
	t_w_t sc_opt, sc;
	mc_reset_z_all(bp->ma, bp->b_aux);
	sc_opt = mc_score_all(bp->ma, bp->b_aux);

	while (1)
	{
		sc = mc_solve_bp_cc(bp);
		// fprintf(stderr, "[M::%s::# round: %u] sc_opt: %ld, sc: %ld\n", __func__, r, sc_opt, sc);
		fprintf(stderr, "[M::%s::# round: %u] sc_opt: %f, sc: %f\n", __func__, r, sc_opt, sc);
		if(sc <= sc_opt) break;
		sc_opt = sc;
		r++;
	}
	fprintf(stderr, "[M::%s::%.3f] ==> round %u\n", __func__, yak_realtime()-index_time, r);
}

t_w_t mc_score_all_advance(const mc_match_t *ma, int8_t *s)
{
	uint32_t k;
	t_w_t z[2], zt = 0;
	for (k = 0; k < ma->n_seq; ++k) 
	{
		uint32_t o = ma->idx.a[k] >> 32;
        uint32_t j, n = (uint32_t)ma->idx.a[k];
		z[0] = z[1] = 0;
		for (j = 0; j < n; ++j) {
            const mc_edge_t *e = &ma->ma.a[o + j];
            uint32_t t = ma_y(*e);
            if (s[t] > 0) z[0] += e->w;
            else if (s[t] < 0) z[1] += e->w;
        }
		zt += -((t_w_t)(s[k])) * (z[0] - z[1]);
	}
	return zt;
}

t_w_t mb_score_all_advance(const mc_match_t *ma, mb_g_t *mbg)
{
	uint32_t k;
	t_w_t z[2], zt = 0;
	for (k = 0; k < ma->n_seq; ++k) 
	{
		uint32_t o = ma->idx.a[k] >> 32;
        uint32_t j, n = (uint32_t)ma->idx.a[k];
		z[0] = z[1] = 0;
		for (j = 0; j < n; ++j) {
            const mc_edge_t *e = &ma->ma.a[o + j];
            uint32_t t = ma_y(*e);
            if (mbg->u->u.a[mbg->u->idx.a[t]>>1].s[mbg->u->idx.a[t]&1] > 0) z[0] += e->w;
            else if (mbg->u->u.a[mbg->u->idx.a[t]>>1].s[mbg->u->idx.a[t]&1] < 0) z[1] += e->w;
        }
		zt += -((t_w_t)(mbg->u->u.a[mbg->u->idx.a[k]>>1].s[mbg->u->idx.a[k]&1])) * (z[0] - z[1]);
	}
	return zt;
}

void print_sc(const mc_opt_t *opt, const mc_g_t *mg, mc_svaux_t *b, t_w_t sc_opt, uint32_t n_iter)
{
	t_w_t w = mc_score(mg->e, b);
	// if(w != sc_opt) fprintf(stderr, "ERROR\n");
	fprintf(stderr, "# iter: %u, sc_opt: %f, sc-local: %f, sc-global: %f\n", n_iter, sc_opt, w, mc_score_all_advance(mg->e, mg->s.a));
}

void print_mc_node(const mc_match_t *ma, mc_svaux_t *b, uint32_t id)
{
	fprintf(stderr, "[M::%s::utg%.6ul-hap%u]\n", __func__, id, b->s[id]>0?1:b->s[id]<0?2:0);
	w_t w[128], z[4];
	uint32_t o, n, i, hn = 4;
	int8_t s;
	for (i = 0; i < hn; i++) w[i] = 0;
	o = ma->idx.a[id] >> 32;
    n = (uint32_t)ma->idx.a[id];
	for (i = 0; i < n; ++i)
	{
		s = b->s[ma_y(ma->ma.a[o + i])];
		w[s>0?1:s<0?2:0] += ma->ma.a[o + i].w;
	}	
	z[0] = z[3] = 0; z[1] = b->z[id].z[0]; z[2] = b->z[id].z[1];
	for (i = 0; i < hn; i++) fprintf(stderr, "w[%u]-%f, z[%u]-%f\n", i, w[i], i, z[i]);
}

uint32_t mc_solve_cc(const mc_opt_t *opt, const mc_g_t *mg, mc_svaux_t *b, uint32_t cc_off, uint32_t cc_size)
{
	uint32_t j, k, n_iter = 0, flush = opt->max_iter * 50;
	t_w_t sc_opt = -(1<<30), sc;///problem-w
	b->cc_off = cc_off, b->cc_size = cc_size;
	if (b->cc_size < 2) return 0;
	sc_opt = mc_init_spin(mg->e, b);
	// print_sc(opt, mg, b, sc_opt, n_iter);
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

    // print_mc_node(mg->e, b, 36880);
	// print_sc(opt, mg, b, sc_opt, n_iter);
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
			// print_sc(opt, mg, b, sc_opt, n_iter);
		} else {
			for (j = 0; j < b->cc_size; ++j) {
				b->s[b->cc_node[j]] = b->s_opt[b->cc_node[j]];
				b->z[b->cc_node[j]] = b->z_opt[b->cc_node[j]];
			}
		}

		if((n_iter%flush) == 0)
		{
			mc_reset_z(mg->e, b);
			sc = mc_score(mg->e, b);

			for (j = 0; j < b->cc_size; ++j) {
				b->z_opt[b->cc_node[j]] = b->z[b->cc_node[j]];
			}
			sc_opt = sc;
		}
	}

	for (j = 0; j < b->cc_size; ++j)
	{
		b->s[b->cc_node[j]] = b->s_opt[b->cc_node[j]];
		b->z[b->cc_node[j]] = b->z_opt[b->cc_node[j]];
	}
		
	return n_iter;
}


void reset_mb_g_t_z(mb_g_t *mbg);
uint32_t mb_solve_cc(const mc_opt_t *opt, mb_g_t *mbg, mb_svaux_t *b, uint32_t cc_off, uint32_t cc_size)
{
	uint32_t j, k, n_iter = 0, flush = opt->max_iter * 50;
	t_w_t sc_opt = -(1<<30), sc;///problem-w
	b->cc_off = cc_off, b->cc_size = cc_size;
	if (b->cc_size <= 2) return 0;

	for (j = 0; j < b->cc_size; ++j) {///how many nodes
		k = (uint32_t)mbg->e->cc[b->cc_off + j];///node id
		b->cc_node[j] = k;
	}

	sc_opt = mb_score(mbg, b);
	for (j = 0; j < b->cc_size; ++j) {///backup s and z in s_opt and z_opt
		memcpy(b->u_opt[b->cc_node[j]].s, b->u[b->cc_node[j]].s, 2*sizeof(mc_node_t));
		memcpy(b->u_opt[b->cc_node[j]].z, b->u[b->cc_node[j]].z, 4*sizeof(t_w_t));
	}

	// print_sc(opt, mg->e, b, sc_opt, n_iter);
	sc = mb_optimize_local(opt, mbg, b, &n_iter);
	if (sc > sc_opt) {
		for (j = 0; j < b->cc_size; ++j) {
			memcpy(b->u_opt[b->cc_node[j]].s, b->u[b->cc_node[j]].s, 2*sizeof(mc_node_t));
			memcpy(b->u_opt[b->cc_node[j]].z, b->u[b->cc_node[j]].z, 4*sizeof(t_w_t));
		}
		sc_opt = sc;
	} else {
		for (j = 0; j < b->cc_size; ++j) {
			memcpy(b->u[b->cc_node[j]].s, b->u_opt[b->cc_node[j]].s, 2*sizeof(mc_node_t));
			memcpy(b->u[b->cc_node[j]].z, b->u_opt[b->cc_node[j]].z, 4*sizeof(t_w_t));
		}
	}
	// mc_reset_z_debug(mg->e, b);
	// print_sc(opt, mg->e, b, sc_opt, n_iter);
	// fprintf(stderr, "\ncc_size: %u, cc_off: %u\n", b->cc_size, b->cc_off);
	for (k = 0; k < (uint32_t)opt->n_perturb; ++k) {
		if (k&1) mb_perturb(opt, mbg, b);
		else mb_perturb_node(opt, mbg, b, 3);
		sc = mb_optimize_local(opt, mbg, b, &n_iter);
		// fprintf(stderr, "(%u) sc_opt: %f, sc: %f\n", k, sc_opt, sc);
		if (sc > sc_opt) {
			for (j = 0; j < b->cc_size; ++j) {
				memcpy(b->u_opt[b->cc_node[j]].s, b->u[b->cc_node[j]].s, 2*sizeof(mc_node_t));
				memcpy(b->u_opt[b->cc_node[j]].z, b->u[b->cc_node[j]].z, 4*sizeof(t_w_t));
			}
			sc_opt = sc;
		} else {
			for (j = 0; j < b->cc_size; ++j) {
				memcpy(b->u[b->cc_node[j]].s, b->u_opt[b->cc_node[j]].s, 2*sizeof(mc_node_t));
				memcpy(b->u[b->cc_node[j]].z, b->u_opt[b->cc_node[j]].z, 4*sizeof(t_w_t));
			}
		}

		if((n_iter%flush) == 0)
		{
			reset_mb_g_t_z(mbg);
			sc = mb_score(mbg, b);

			for (j = 0; j < b->cc_size; ++j) {
				memcpy(b->u_opt[b->cc_node[j]].z, b->u[b->cc_node[j]].z, 4*sizeof(t_w_t));
			}
			sc_opt = sc;
		}

		// print_sc(opt, mg->e, b, sc_opt, n_iter);
	}
	for (j = 0; j < b->cc_size; ++j)
	{
		memcpy(b->u[b->cc_node[j]].s, b->u_opt[b->cc_node[j]].s, 2*sizeof(mc_node_t));
		memcpy(b->u[b->cc_node[j]].z, b->u_opt[b->cc_node[j]].z, 4*sizeof(t_w_t));
	}
		
	return n_iter;
}

void reset_mb_g_t_z(mb_g_t *mbg)
{
	uint32_t k;
	mb_match_t *ma = mbg->e;
	for (k = 0; k < ma->n_seq; ++k) ///each block
	{
		uint32_t o = ma->idx.a[k] >> 32;
        uint32_t j, n = (uint32_t)ma->idx.a[k];
		mbg->u->u.a[k].z[0] = mbg->u->u.a[k].z[1] = 0;
		mbg->u->u.a[k].z[2] = mbg->u->u.a[k].z[3] = 0;
		for (j = 0; j < n; ++j) {
            const mb_edge_t *e = &ma->ma.a[o + j];
            uint32_t t = ma_y(*e);

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

void debug_mbg(mb_g_t *mbg, mc_g_t *mg, mc_svaux_t *b)
{
	uint32_t i, k, a_n[2], *a[2], qn, found;
	int8_t s[2];
	for (i = 0; i < mbg->u->u.n; i++)///each block
	{
		fprintf(stderr, "i: %u, mbg->u->u.n: %u\n", i, (uint32_t)mbg->u->u.n);
		decode_mb_node(mbg, i, &(a[0]), &(a_n[0]), &(s[0]), &(a[1]), &(a_n[1]), &(s[1]));

		found = 0;
		for (k = 0; k < a_n[0]; k++)
        {
            qn = a[0][k];
			if(b->s[qn] == s[0]) found = 1;
			if(b->s[qn] == 0) continue;
			if(b->s[qn] != s[0]) fprintf(stderr, "ERROR-0-::%s\n", __func__);
		}
		if(found == 0) fprintf(stderr, "ERROR-0-*::%s\n", __func__);

		found = 0;
		for (k = 0; k < a_n[1]; k++)
        {
            qn = a[1][k];
			if(b->s[qn] == s[1]) found = 1;
			if(b->s[qn] == 0) continue;
			if(b->s[qn] != s[1]) fprintf(stderr, "ERROR-1-::%s\n", __func__);
		}
		if(found == 0 && a_n[1] > 0) fprintf(stderr, "ERROR-1-*::%s\n", __func__);
	}

	// t_w_t mcw = mc_score_all(mg->e, b);
	// t_w_t mbw = 0;
	// mb_match_t *ma = mbg->e;
	// for (k = 0; k < ma->n_seq; ++k) {
	// 	///a[0]
	// 	mbw += -(t_w_t)(mbg->u->u.a[k].s[0]) * (mbg->u->u.a[k].z[0] - mbg->u->u.a[k].z[1]);
	// 	///a[1]
	// 	mbw += -(t_w_t)(mbg->u->u.a[k].s[1]) * (mbg->u->u.a[k].z[2] - mbg->u->u.a[k].z[3]);
	// }
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
	// memcpy(b->s_opt, b->s, sizeof(int8_t)*mg->e->n_seq);

	///adjust by block
	kvec_t(uint32_t) s; kv_init(s);
	uint32_t ws, we, n_flip;
	uint8_t *vis = NULL; CALLOC(vis, mg->e->n_seq);
	t_w_t w = 0, max_w = 0;
	for (i = 0; i < mbg->u->u.n; i++)///each block
	{
		decode_mb_node(mbg, i, &(a[0]), &(a_n[0]), NULL, &(a[1]), &(a_n[1]), NULL);
		///if(a_n[0] == 0 || a_n[1] == 0) continue;
		if(a_n[0] + a_n[1] <= 1) continue;
		s.n = 0;
		for (k = 0; k < a_n[0]; k++)
        {
            qn = a[0][k];
			if(b->s[qn] == 0) continue;
			// qn <<= 1;
			if(b->s[qn] < 0) qn += ((uint32_t)1<<31);			
			kv_push(uint32_t, s, qn);
		}

		for (k = 0; k < a_n[1]; k++)
        {
            qn = a[1][k];
			if(b->s[qn] == 0) continue;
			// qn <<= 1; qn++;
			if(b->s[qn] > 0) qn += ((uint32_t)1<<31);
			kv_push(uint32_t, s, qn);
		}

		if(s.n == 0) continue;
		radix_sort_mc32(s.a, s.a + s.n);

		ws = we = (uint32_t)-1; n_flip = 0;
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
				n_flip++;
			}
		}

		if(n_flip <= 1) continue;

		/*******************************for debug************************************/
		// t_w_t sc_opt = mc_score_all(mg->e, b);
		/*******************************for debug************************************/
		for (m = ws; m < we; m++) mc_set_spin(mg->e, b, (s.a[m]<<1)>>1, -b->s[(s.a[m]<<1)>>1]);
		/*******************************for debug************************************/
		// t_w_t sc_cur = mc_score_all(mg->e, b);
		// if(sc_cur - sc_opt != max_w*2)
		// {
		// 	fprintf(stderr, "ERROR-0-::%s\n", __func__);
		// 	// fprintf(stderr, "sc_opt: %ld, sc_cur: %ld, max_w: %ld\n", sc_opt, sc_cur, max_w);
		// 	fprintf(stderr, "sc_opt: %f, sc_cur: %f, max_w: %f\n", sc_opt, sc_cur, max_w);
		// }
		/*******************************for debug************************************/
	}
	kv_destroy(s); free(vis);

	for (i = 0; i < mbg->u->u.n; i++)///each block
	{
		decode_mb_node(mbg, i, &(a[0]), &(a_n[0]), NULL, &(a[1]), &(a_n[1]), NULL);
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

	/*******************************for debug************************************/
	// debug_mbg(mbg, mg, b);
	/*******************************for debug************************************/
	// memcpy(b->s, b->s_opt, sizeof(int8_t)*mg->e->n_seq);
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

void mc_set_by_mbg(mc_g_t *mg, mb_g_t *mbg)
{
	// t_w_t z0 = mc_score_all_advance(mg->e, mg->s.a);
	// t_w_t z1 = mb_score_all_advance(mg->e, mbg);
	// if(z0 >= z1) return;
	uint32_t i, k, qn, *a[2], a_n[2];
	int8_t s[2];
	for (i = 0; i < mbg->u->u.n; i++)
	{
		decode_mb_node(mbg, i, &(a[0]), &(a_n[0]), &(s[0]), &(a[1]), &(a_n[1]), &(s[1]));
		for (k = 0; k < a_n[0]; k++)
		{
			qn = a[0][k];
			if(mg->s.a[qn] == 0) continue;
			mg->s.a[qn] = s[0];
		}

		for (k = 0; k < a_n[1]; k++)
		{
			qn = a[1][k];
			if(mg->s.a[qn] == 0) continue;
			mg->s.a[qn] = s[1];
		}
	}
}


void debug_mb_solve_core(mb_g_t *mbg)
{
	uint32_t i;
	for (i = 0; i < mbg->u->u.n; i++)
	{
		debug_mb_z(mbg, i);
	}
}

void print_mb_g_blcok(mb_g_t *mbg)
{
	uint32_t i, k, a_n[2], *a[2];
	for (i = 0; i < mbg->u->u.n; i++)///each block
	{
		decode_mb_node(mbg, i, &(a[0]), &(a_n[0]), NULL, &(a[1]), &(a_n[1]), NULL);
		{
			if(a_n[0] + a_n[1] <= 1) continue;
			fprintf(stderr, "\nB0-%s\n", a_n[0]>1?"mul":"single");
			for (k = 0; k < a_n[0]; k++)
        	{
				fprintf(stderr,"s-utg%.6ul\n", a[0][k]+1);
			}

			fprintf(stderr, "B1-%s\n", a_n[0]>1?"mul":"single");
			for (k = 0; k < a_n[1]; k++)
			{
				fprintf(stderr,"d-utg%.6ul\n", a[1][k]+1);
			}
		}
	}
}



void mb_solve_core(mc_opt_t *opt, mc_g_t *mg, kv_u_trans_t *ref, uint32_t is_sys)
{
	if(!ref) return;
	double index_time = yak_realtime();
	uint32_t st, i;
	mb_g_t *mbg = init_mb_g_t(mg->e, ref, is_sys);
	mb_svaux_t *bb;
	/**************************init**************************/
	if(VERBOSE_CUT)
	{
		fprintf(stderr, "\n\n\n\n\n*************beg-[M::%s::score->%f] ==> Partition\n", __func__, mc_score_all_advance(mg->e, mg->s.a));
	} 
	
	mc_svaux_t *b;
	mc_g_cc(mg->e);
	b = mc_svaux_init(mg, opt->seed);
	mc_init_spin_all(opt, mg, mbg, b);
	mc_svaux_destroy(b);
	free(mg->e->cc); 
	mg->e->cc = NULL;
	/**************************init**************************/
	mb_g_cc(mbg);
	bb = mb_svaux_init(mbg, opt->seed);
	
	if(VERBOSE_CUT)
	{
		fprintf(stderr, "*********before-[M::%s::mc_score->%f] ==> Partition\n", __func__, mc_score_all_advance(mg->e, mg->s.a));
		fprintf(stderr, "*********before-[M::%s::mb_score->%f] ==> Partition\n", __func__, mb_score_all_advance(mg->e, mbg));
		/*******************************for debug************************************/
		// print_mb_g_blcok(mbg);
		/*******************************for debug************************************/
	}
	
	
	opt->n_perturb = opt->n_b_perturb;
	for (st = 0, i = 1; i <= mbg->e->n_seq; ++i) {
		if (i == mbg->e->n_seq || mbg->e->cc[st]>>32 != mbg->e->cc[i]>>32) {
			mb_solve_cc(opt, mbg, bb, st, i - st);
			st = i;
		}
	}
	opt->n_perturb = opt->n_s_perturb - opt->n_b_perturb;
	
	mc_set_by_mbg(mg, mbg);
	if(VERBOSE_CUT)
	{
		/*******************************for debug************************************/
		// debug_mb_solve_core(mbg);
		/*******************************for debug************************************/
		fprintf(stderr, "##############end-[M::%s::score->%f] ==> Partition\n", __func__, mc_score_all_advance(mg->e, mg->s.a));
	}
	destory_mb_g_t(&mbg);
	mb_svaux_destroy(bb);
	fprintf(stderr, "[M::%s::%.3f] ==> Partition\n", __func__, yak_realtime()-index_time);
}

void mc_solve_core(const mc_opt_t *opt, mc_g_t *mg, bubble_type* bub)
{
	double index_time = yak_realtime();
	uint32_t st, i;
	mc_svaux_t *b;
	mc_bp_t *bp = NULL;
	mc_g_cc(mg->e);
	b = mc_svaux_init(mg, opt->seed);
	if(bub) bp = mc_bp_t_init(mg->e, b, bub, asm_opt.thread_num);
	/*******************************for debug************************************/
	if(bp) mc_init_spin_all(opt, mg, NULL, b);
	if(bp) mc_solve_bp(bp);
	/*******************************for debug************************************/
	if(VERBOSE_CUT)
	{
		fprintf(stderr, "\n\n\n\n\n*************beg-[M::%s::score->%f] ==> Partition\n", __func__, mc_score_all_advance(mg->e, mg->s.a));
	}
	
	for (st = 0, i = 1; i <= mg->e->n_seq; ++i) {
		if (i == mg->e->n_seq || mg->e->cc[st]>>32 != mg->e->cc[i]>>32) {
			mc_solve_cc(opt, mg, b, st, i - st);
			st = i;
		}
	}

	if(VERBOSE_CUT)
	{
		fprintf(stderr, "##############end-[---M::%s::score->%f] ==> Partition\n", __func__, mc_score_all(mg->e, b));
	}
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
		if(t_ch->ir_het[u->a[i]>>33] == N_HET) continue;
		if(s == 0)
		{
			if(t_ch->ir_het[u->a[i]>>33]&P_HET)//special case
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


void filter_ta_by_mc(ma_ug_t *ug, asg_t *read_g, uint32_t uID, kv_u_trans_t* ta, mc_match_t *ma, int8_t *s)
{
	mc_edge_t *o = pt_a(*ma, uID);
	uint32_t n = pt_n(*ma, uID), k, qn, tn;
	u_trans_t *p = NULL;
	for (k = 0; k < n; ++k)
	{
		qn = ma_x(o[k]); tn = ma_y(o[k]); p = NULL;
		if((s[qn]*s[tn])!=-1) continue;
		get_u_trans_spec(ta, qn, tn, &p, NULL);
		if(p && p->nw == o[k].w) {
			p->del = 0; 
		}
		else {
			get_u_trans_spec(ta, tn, qn, &p, NULL);
			if(p && p->nw == o[k].w) p->del = 0;
		}
		if(!p) fprintf(stderr, "ERROR-ta-p\n");
	}
}


void clean_ta_by_mc(mc_g_t *mg, kv_u_trans_t *ta)
{
	uint32_t v, i;
	u_trans_t *p = NULL;
	for (i = 0; i < ta->n; i++) ta->a[i].del = 1;
	for (i = 0; i < mg->e->n_seq; ++i) filter_ta_by_mc(mg->ug, mg->rg, i, ta, mg->e, mg->s.a);
	for (i = v = 0; i < ta->n; i++) {
		if(ta->a[i].del) continue;
		ta->a[v++] = ta->a[i];
	}
	ta->n = v;
	for (i = 0; i < v; i++) {
		kv_pushp(u_trans_t, *ta, &p);
		(*p) = ta->a[i];
		p->qn = ta->a[i].tn; p->qs = ta->a[i].ts; p->qe = ta->a[i].te;
        p->tn = ta->a[i].qn; p->ts = ta->a[i].qs; p->te = ta->a[i].qe;
	}
	kt_u_trans_t_idx(ta, mg->ug->g->n_seq);
}

void p_nodes(mc_g_t *mg, trans_chain* t_ch, uint8_t* trio_flag)
{
	uint32_t i;
	for (i = 0; i < mg->e->n_seq; ++i) 
	{
		set_p_flag(mg, i, trio_flag, t_ch, mg->s.a[i]);
	}
}

void write_mc_g_t(mc_opt_t *opt, mc_g_t *mg, const char *name)
{
	FILE* fp = fopen(name, "w");

    fwrite(opt, sizeof(mc_opt_t), 1, fp);
	fwrite(&(mg->s.n), sizeof(mg->s.n), 1, fp);
	fwrite(mg->s.a, sizeof(mc_node_t), mg->s.n, fp);
	fwrite(&(mg->e->n_seq), sizeof(mg->e->n_seq), 1, fp);
	fwrite(&(mg->e->ma.n), sizeof(mg->e->ma.n), 1, fp);
	fwrite(mg->e->ma.a, sizeof(mc_edge_t), mg->e->ma.n, fp);
	fwrite(&(mg->e->idx.n), sizeof(mg->e->idx.n), 1, fp);
	fwrite(mg->e->idx.a, sizeof(uint64_t), mg->e->idx.n, fp);

    fclose(fp);
}

mc_g_t* load_mc_g_t(mc_opt_t *opt, const char *name)
{
	FILE* fp = NULL; 
    fp = fopen(name, "r"); 
    if(!fp) return NULL;

	uint64_t flag = 0;
	mc_g_t *mg = NULL; CALLOC(mg, 1);
	kv_init(mg->s); CALLOC(mg->e, 1);

    flag += fread(opt, sizeof(mc_opt_t), 1, fp);

	flag += fread(&(mg->s.n), sizeof(mg->s.n), 1, fp);
	mg->s.m = mg->s.n; MALLOC(mg->s.a, mg->s.n);
	flag += fread(mg->s.a, sizeof(mc_node_t), mg->s.n, fp);

	flag += fread(&(mg->e->n_seq), sizeof(mg->e->n_seq), 1, fp);

	flag += fread(&(mg->e->ma.n), sizeof(mg->e->ma.n), 1, fp);
	mg->e->ma.m = mg->e->ma.n; MALLOC(mg->e->ma.a, mg->e->ma.n);
	flag += fread(mg->e->ma.a, sizeof(mc_edge_t), mg->e->ma.n, fp);

	flag += fread(&(mg->e->idx.n), sizeof(mg->e->idx.n), 1, fp);
	mg->e->idx.m = mg->e->idx.n; MALLOC(mg->e->idx.a, mg->e->idx.n);
	flag += fread(mg->e->idx.a, sizeof(uint64_t), mg->e->idx.n, fp);

    fclose(fp);
	return mg;
}

void debug_mc_g_t(const char* name)
{
	mc_opt_t opt;
	mc_g_t *mg = load_mc_g_t(&opt, name);
	mc_solve_core(&opt, mg, NULL);
	destory_mc_g_t(&mg);
	exit(1);
}

void print_hap_s(int8_t *s, uint32_t sn)
{
	fprintf(stderr, "\n[M::%s]\n", __func__);
	uint32_t i;
	for (i = 0; i < sn; i++)
	{
		fprintf(stderr,"utg%.6ul\t", i + 1);
		if(s[i] > 0) fprintf(stderr, "h%u\n", 1);
		else if(s[i] < 0) fprintf(stderr, "h%u\n", 2);
		else fprintf(stderr,"\n");
	}
}

void mc_solve(hap_overlaps_list* ovlp, trans_chain* t_ch, kv_u_trans_t *ta, ma_ug_t *ug, asg_t *read_g, double f_rate, uint8_t* trio_flag, uint32_t renew_s, int8_t *s, uint32_t is_sys, bubble_type* bub, kv_u_trans_t *ref, int clean_ov)
{
	mc_opt_t opt;
	mc_opt_init(&opt, asm_opt.n_perturb, asm_opt.f_perturb, asm_opt.seed);
	mc_g_t *mg = init_mc_g_t(ug, read_g, s, renew_s);
	update_mc_edges(mg, ovlp, ta, t_ch, f_rate, is_sys);

	// fprintf(stderr, "[M::%s:: # edges: %u]\n", __func__, (uint32_t)mg->e->ma.n);
	
	mb_solve_core(&opt, mg, ref, is_sys);
	///debug_mc_g_t(mg);
	// if(renew_s == 0) write_mc_g_t(&opt, mg, MC_NAME);
	mc_solve_core(&opt, mg, bub);

	if((asm_opt.flag & HA_F_PARTITION) && t_ch)
	{
		p_nodes(mg, t_ch, trio_flag);
	}

	if(clean_ov){
		if(ovlp) clean_ovlp_by_mc(mg, ovlp);
		if(ta) clean_ta_by_mc(mg, ta);
	}
	// print_hap_s(s, ug->u.n);

	destory_mc_g_t(&mg);	
}



void comp(int m, int N, int M, mcb_t *p, int *c)
{
    if (m == M + 1)
    {
		int i;
		mcg_node_t x = 0;
		for (i = 0; i < M; i++) x |= ((mcg_node_t)1<<(c[i+1]-1));
		kv_push(mcg_node_t, *p, x);
    }
    else
    {
        for (c[m] = c[m - 1] + 1; c[m] <= N - M + m; c[m]++) 
        {
            comp(m + 1, N, M, p, c);
        }
    }
}


void get_mcb(uint32_t n, uint32_t m, mcb_t *p, int *c)
{
	memset(c, 0, sizeof(int)*n+1);
	p->n = 0;
	comp(1, n, m, p, c);
}

mc_gg_t *init_mc_gg_t(uint32_t un, kv_gg_status *s, uint16_t hapN)
{
	uint32_t i; 
	int *c = NULL; CALLOC(c, hapN+1);
	mc_gg_t *p = NULL; CALLOC(p, 1);
	// p->ug = ug; p->rg = read_g; 
	p->un = un; p->s = s; p->hN = hapN;
	CALLOC(p->m.a, p->hN); p->m.n = p->m.m = hapN;
	for (i = 0; i < p->hN; i++) get_mcb(hapN, i+1, &(p->m.a[i]), c);
	p->mask = (1<<hapN); p->mask--;
	free(c);
	return p;
}

void destory_mc_gg_t(mc_gg_t **p)
{
	uint32_t i; 
	if(!p || !(*p)) return;
	for (i = 0; i < (*p)->m.m; i++)
	{
		free((*p)->m.a[i].a);
	}
	free((*p)->m.a);
	
	if((*p)->e)
	{
		kv_destroy((*p)->e->idx);
		kv_destroy((*p)->e->ma);
		free((*p)->e->cc);
		free((*p)->e);
	}
	free((*p));
}

kv_gg_status *init_mc_gg_status(ma_ug_t *ug, asg_t *read_g, ma_sub_t* coverage_cut, 
ma_hit_t_alloc* sources, R_to_U* ruIndex, uint64_t t_cov, uint16_t hapN)
{
	fprintf(stderr, "t_cov-%lu\n", t_cov);
	uint64_t *covs = NULL, i, k, k_i, c = t_cov/hapN, c_min, c_max; 
	uint32_t len[2];
	uint8_t *vis = NULL; CALLOC(vis, read_g->n_seq);
	kv_gg_status *p = NULL; 
	CALLOC(covs, hapN);
	for (i = 0; i < hapN; i++)
	{
		c_min = ((i+1)*c) - (0.6*c);
		c_max = ((i+1)*c) + (0.6*c);
		if(i == 0) c_min = 0;
		if(i == hapN) c_max = (uint32_t)-1;
		covs[i] = (c_max<<32)|c_min; 
	}
	CALLOC(p, 1);
	p->n = p->m = ug->u.n;
	CALLOC(p->a, p->n);
	for (i = 0; i < p->n; i++)
	{
		c = get_utg_cov(ug, i, read_g, coverage_cut, sources, ruIndex, vis);
		p->a[i].h[0] = p->a[i].h[1] = (uint16_t)-1; p->a[i].s = 0; k_i = 0;
		p->a[i].hw[0] = 1; p->a[i].hw[1] = 0; p->a[i].hc = 0;
		for (k = 0; k < hapN; k++)
		{
			c_min = (uint32_t)covs[k]; c_max = covs[k]>>32;
			if(c < c_min || c >= c_max) continue;
			p->a[i].h[k_i++] = k + 1;
		}
		if(p->a[i].h[1] == (uint16_t)-1) continue;
		
		len[0] = (c >= (p->a[i].h[0]*(t_cov/hapN))? 
							c - (p->a[i].h[0]*(t_cov/hapN)) : (p->a[i].h[0]*(t_cov/hapN)) - c);
		len[1] = (c >= (p->a[i].h[1]*(t_cov/hapN))? 
							c - (p->a[i].h[1]*(t_cov/hapN)) : (p->a[i].h[1]*(t_cov/hapN)) - c);
		if(len[0] > len[1])
		{
			k_i = p->a[i].h[0];
			p->a[i].h[0] = p->a[i].h[1];
			p->a[i].h[1] = k_i;

			k_i = len[0];
			len[0] = len[1];
			len[1] = k_i;
		}
		
		p->a[i].hw[0] = (double)len[1]/(double)(len[0]+len[1]);
		p->a[i].hw[1] = (double)len[0]/(double)(len[0]+len[1]);
	}
	free(covs); free(vis);
	return p;
}

void update_mc_edges_general(mc_gg_t *mg, kv_u_trans_t *ta, uint16_t hapN)
{
	uint32_t i, k;
	mc_edge_t *ma = NULL;
	CALLOC(mg->e, 1);
	mg->e->n_seq = mg->un;
	kv_init(mg->e->idx); kv_init(mg->e->ma); 

	for (i = 0; i < ta->n; ++i)
	{
		if(ta->a[i].del) continue;
		if(mg->s->a[ta->a[i].qn].h[0] >= hapN && mg->s->a[ta->a[i].qn].h[1] >= hapN) continue;
		if(mg->s->a[ta->a[i].tn].h[0] >= hapN && mg->s->a[ta->a[i].tn].h[1] >= hapN) continue;
		kv_pushp(mc_edge_t, mg->e->ma, &ma);
		ma->x = (uint64_t)ta->a[i].qn << 32 | ta->a[i].tn;
		ma->w = w_cast((ta->a[i].nw));
	}

	for (i = k = 0; i < mg->e->ma.n; i++)
    {
        if(mg->e->ma.a[i].w == 0) continue;
        mg->e->ma.a[k] = mg->e->ma.a[i];
        k++;
    }
    mg->e->ma.n = k;
    
    radix_sort_mce(mg->e->ma.a, mg->e->ma.a + mg->e->ma.n);
    mc_edges_idx(mg->e);
    mc_edges_symm(mg->e);
}


typedef struct {
	t_w_t *z;
	uint32_t hapN;
} mc_poy_t;

typedef struct {
	uint64_t x; // RNG
	uint32_t cc_off, cc_size;
	kvec_t(uint64_t) cc_edge;
	uint32_t *cc_node;
	uint32_t *bfs, *bfs_mark;
	mc_poy_t *z, *z_opt;///keep scores to nodes(1) and nodes(-1)
	mc_gg_status *s, *s_opt;
	mcb_t *m;
	mcg_node_t mask;
	uint32_t hapN;
} mcgg_svaux_t;

mc_poy_t *init_mc_poy_t(uint32_t un, uint32_t hapN)
{
	uint32_t i;
	mc_poy_t *z = NULL; CALLOC(z, 1);
	MALLOC(z->z, ((uint32_t)un<<hapN));
	for (i = 0; i < ((uint32_t)un<<hapN); i++) z->z[i] = 0;
	z->hapN = hapN;
	return z;
}

void destroy_mc_poy_t(mc_poy_t **p)
{
	if(!p || !(*p)) return;
	free((*p)->z); free(*p);
}

mcgg_svaux_t *mcgg_svaux_init(const mc_gg_t *mg, uint64_t x, uint32_t hapN)
{
    uint32_t st, i, max_cc = 0;
    mc_match_t *ma = mg->e;
    mcgg_svaux_t *b;
    CALLOC(b, 1);
    b->x = x;
    for (st = 0, i = 1; i <= ma->n_seq; ++i)
        if (i == ma->n_seq || ma->cc[st]>>32 != ma->cc[i]>>32)
            max_cc = max_cc > i - st? max_cc : i - st, st = i;
    kv_init(b->cc_edge);
    MALLOC(b->cc_node, max_cc);
    b->s = mg->s->a;
    CALLOC(b->s_opt, ma->n_seq);
    MALLOC(b->bfs, ma->n_seq);

    MALLOC(b->bfs_mark, ma->n_seq); 
    memset(b->bfs_mark, -1, ma->n_seq*sizeof(uint32_t));
    
	b->z = init_mc_poy_t(ma->n_seq, hapN);
	b->z_opt = init_mc_poy_t(ma->n_seq, hapN);

	b->m = mg->m.a;
	b->mask = ((mcg_node_t)1)<<hapN; b->mask--;
	b->hapN = hapN;
    return b;
}

void mcgg_svaux_destroy(mcgg_svaux_t *b)
{
	b->s = NULL;
    kv_destroy(b->cc_edge); free(b->cc_node);
    free(b->s); free(b->s_opt);
    destroy_mc_poy_t(&(b->z)); 
	destroy_mc_poy_t(&(b->z_opt));
    free(b->bfs); free(b->bfs_mark);
    free(b);
}

static inline mcg_node_t kr_drand_node(uint64_t id, mcgg_svaux_t *b, uint16_t *hc)
{
	uint16_t hn = b->s[id].h[0]-1;
	if(hc) (*hc) = 0;
	b->x = kr_splitmix64(b->x);
	if(b->s[id].h[1] != (uint16_t)-1)
	{
		union { uint64_t i; double d; } u;
		u.i = 0x3FFULL << 52 | (b->x) >> 12;
		if((u.d - 1.0) > b->s[id].hw[0])
		{
			hn = b->s[id].h[1]-1;
			if(hc) (*hc) = 1;
		}
	} 
    return b->m[hn].a[(b->x)%b->m[hn].n];
}

static inline mcg_node_t kr_drand_node_ref(uint64_t id, mcgg_svaux_t *b, uint16_t *hc, mcg_node_t ref, uint64_t rev)
{
	uint64_t i, k, m, mm, mn, mi, cn, hn = b->s[id].h[0];
	mcg_node_t t;
	if(hc) (*hc) = 0;
	b->x = kr_splitmix64(b->x);
	if(b->s[id].h[1] != (uint16_t)-1)
	{
		union { uint64_t i; double d; } u;
		u.i = 0x3FFULL << 52 | (b->x) >> 12;
		if((u.d - 1.0) > b->s[id].hw[0])
		{
			hn = b->s[id].h[1];
			if(hc) (*hc) = 1;
		}
	} 
	if(rev) ref ^= (mcg_node_t)-1;
	ref &= b->mask;

	t = ref; cn = 0;
	while (t)
	{	
		cn += t&1;
		t >>= 1;
	} 
	if(cn == hn) return ref;

	if(cn>=hn) mm=cn, m=1, mn=cn-hn;///1->0
	else mm=b->hapN-cn, m=0, mn=hn-cn;///0->1

	for (i = 0; i < mn; i++)
	{
		mi = (b->x%mm);
		for (k = 0; k < b->hapN; k++)
		{
			if(((ref>>k)&1)!=m) continue;
			if(mi == 0)
			{
				ref ^= ((mcg_node_t)1<<k);
				break;
			}
			mi--;
		}
		mm--;
	}
    
	return ref;
}


void debug_hapM(mc_gg_status *s, const char* cmd)
{
	uint32_t i, cn, hn = s->h[s->hc];
	mcg_node_t ref = s->s;
	for (i = cn = 0; i < 32; i++) cn += ((ref>>i)&1);
	if (cn != hn) fprintf(stderr, "%s-ERROR-cn, cn-%u, hn-%u\n", cmd, cn, hn);
	else fprintf(stderr, "%s-pass-cn, cn-%u, hn-%u\n", cmd, cn, hn);
}


void mcgg_reset_z(const mc_match_t *ma, mcgg_svaux_t *b)
{
	uint32_t i;
	for (i = 0; i < b->cc_size; ++i) {
		uint32_t k = (uint32_t)ma->cc[b->cc_off + i];///uid
		uint32_t o = ma->idx.a[k] >> 32;
		uint32_t j, n = (1<<b->z->hapN);
		for (j = 0; j < n; j++) mcp_de(*(b->z), k, j) = 0;
		n = (uint32_t)ma->idx.a[k];
		for (j = 0; j < n; ++j) {
			const mc_edge_t *e = &ma->ma.a[o + j];
			uint32_t t = ma_y(*e);
			mcp_de(*(b->z), k, b->s[t].s) += e->w;
		}
	}
}

t_w_t mcgg_score(const mc_match_t *ma, mcgg_svaux_t *b)
{
	uint32_t i, j, n = (1<<b->z->hapN);
	t_w_t z = 0;
	for (i = 0; i < b->cc_size; ++i) {
		uint32_t k = (uint32_t)ma->cc[b->cc_off + i];///uid
		for (j = 0; j < n; j++)
		{
			z += ((b->s[k].s&((mcg_node_t)j))?-mcp_de(*(b->z), k, j):mcp_de(*(b->z), k, j));
		}		
	}
	return z;
}

t_w_t mcgg_init_spin(const mc_match_t *ma, mcgg_svaux_t *b)
{
	uint32_t i;
	b->cc_edge.n = 0;
	for (i = 0; i < b->cc_size; ++i) {///how many nodes
		uint32_t k = (uint32_t)ma->cc[b->cc_off + i];///node id
		b->cc_node[i] = k;
		if(b->s[k].s == 0) break;
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
		if (b->s[n1].s == 0 && b->s[n2].s == 0) {
			b->s[n1].s = kr_drand_node(n1, b, &(b->s[n1].hc));
			// debug_hapM(&(b->s[n1]), "s0");
			b->s[n2].s = kr_drand_node_ref(n2, b, &(b->s[n2].hc), b->s[n1].s, e->w>0?1:0);
			// debug_hapM(&(b->s[n2]), "s1");
		}
		else if(b->s[n1].s == 0)
		{
			b->s[n1].s = kr_drand_node_ref(n1, b, &(b->s[n1].hc), b->s[n2].s, e->w>0?1:0);	
			// debug_hapM(&(b->s[n1]), "s2");	
		}
		else if(b->s[n2].s == 0)
		{
			b->s[n2].s = kr_drand_node_ref(n2, b, &(b->s[n2].hc), b->s[n1].s, e->w>0?1:0);
			// debug_hapM(&(b->s[n2]), "s3");
		}
	}

	passed:
	mcgg_reset_z(ma, b);
	return mcgg_score(ma, b);
}

static mcg_node_t get_max_m(uint64_t id, mcgg_svaux_t *b)
{
	uint32_t hn = b->s[id].h[0]-1, k, j, n = (1<<b->z->hapN);
	mcg_node_t m, *p = NULL;
	t_w_t z, max_z = -(1<<30);
	for (k = 0; k < b->m[hn].n; k++)
	{
		m = b->m[hn].a[k]; z = 0;
		for (j = 0; j < n; j++)
		{
			z += ((m&((mcg_node_t)j))?-mcp_de(*(b->z), id, j):mcp_de(*(b->z), id, j));
		}
		if(!p || max_z < z || (max_z == z && m == b->s[id].s)) p = &(b->m[hn].a[k]), max_z = z;
	}

	if(b->s[id].h[1] != (uint16_t)-1)
	{
		hn = b->s[id].h[1]-1;
		for (k = 0; k < b->m[hn].n; k++)
		{
			m = b->m[hn].a[k]; z = 0;
			for (j = 0; j < n; j++)
			{
				z += ((m&((mcg_node_t)j))?-mcp_de(*(b->z), id, j):mcp_de(*(b->z), id, j));
			}
			if(!p || max_z < z || (max_z == z && m == b->s[id].s)) p = &(b->m[hn].a[k]), max_z = z;
		}
	}
	return (*p);
}

///k is uid
static void mcgg_set_spin(const mc_match_t *ma, mcgg_svaux_t *b, uint32_t k, mcg_node_t s, const char* cmd)
{
	uint32_t o, j, n;
	mcg_node_t s0 = b->s[k].s;
	/*******************************for debug************************************/
	// mcg_node_t t = s;
	// o = 0;
	// while (t) {
	// 	o += (t&1); t>>=1;
	// }
	// if(o != b->s[k].h[0] && o != b->s[k].h[1]) fprintf(stderr, "cmd-%s, ERROR-mcgg-1\n", cmd);
	// if(s0 == s) fprintf(stderr, "cmd-%s, ERROR-mcgg-2\n", cmd);
	/*******************************for debug************************************/
	if (s0 == s) return;
	o = ma->idx.a[k] >> 32;
	n = (uint32_t)ma->idx.a[k];
	for (j = 0; j < n; ++j) {
		const mc_edge_t *e = &ma->ma.a[o + j];
		uint32_t t = ma_y(*e);///1->z[0]; (-1)->z[1];
		mcp_de(*(b->z), t, s0) -= e->w; 
		mcp_de(*(b->z), t, s) += e->w; 
	}
	b->s[k].s = s;
}

static t_w_t mcgg_optimize_local(const mc_opt_t *opt, const mc_match_t *ma, mcgg_svaux_t *b, uint32_t *n_iter)
{
	uint32_t i, n_flip = 0;
	int32_t n_iter_local = 0;
	mcg_node_t ms;
	while (n_iter_local < opt->max_iter) {
		++(*n_iter);
		ks_shuffle_uint32_t(b->cc_size, b->cc_node, &b->x);
		for (i = n_flip = 0; i < b->cc_size; ++i) {
			uint32_t k = b->cc_node[i];///uid
			ms = get_max_m(k, b);
			if(ms != b->s[k].s)
			{
				mcgg_set_spin(ma, b, k, ms, __func__);///no need to change the score of k itself
				// debug_hapM(&(b->s[k]), "s4");
				++n_flip;
			}
		}
		++n_iter_local;
		if (n_flip == 0) break;
	}

	return mcgg_score(ma, b);
}

void inline back_status(mcgg_svaux_t *b, uint32_t id, uint32_t to_opt)
{
	if(to_opt)
	{
		b->s_opt[id] = b->s[id];
		memcpy(b->z_opt->z+(id<<b->z->hapN), b->z->z+(id<<b->z->hapN), (1<<b->z->hapN)*sizeof(t_w_t));
	}
	else
	{
		b->s[id] = b->s_opt[id];
		memcpy(b->z->z+(id<<b->z->hapN), b->z_opt->z+(id<<b->z->hapN), (1<<b->z->hapN)*sizeof(t_w_t));
	}
	
}
/**
static inline mcg_node_t kr_drand_node_new(uint64_t id, mcgg_svaux_t *b)
{
	if(id == 10070) fprintf(stderr, "id-%lu, h[0]-%u, h[1]-%u\n", id, b->s[id].h[0], b->s[id].h[1]);
	uint32_t hn = b->s[id].h[0]-1, is_old = 1, k;
	if((b->s[id].h[1] != (uint16_t)-1) && (kr_drand_r(&b->x) > b->s[id].hw[0]))
	{
		hn = b->s[id].h[1]-1; is_old = 0;
	}
	if(id == 10070) fprintf(stderr, "id-%lu, hn-%u, b->m[hn].n-%u, is_old-%u\n", id, hn, b->m[hn].n, is_old);
	b->x = kr_splitmix64(b->x);
	k = b->x%(b->m[hn].n-is_old);
	if(id == 10070) fprintf(stderr, "id-%lu, hn-%u, k-%u\n", id, hn, k);
	if(b->m[hn].a[k] == b->s[id].s) k = b->m[hn].n-1;
	return b->m[hn].a[k];
}
**/

static inline mcg_node_t kr_drand_node_new(uint64_t id, mcgg_svaux_t *b)
{
	uint32_t hn = b->s[id].h[0]-1, k;
	if((b->s[id].h[1] != (uint16_t)-1) && (kr_drand_r(&b->x) > b->s[id].hw[0]))
	{
		hn = b->s[id].h[1]-1;
	}
	b->x = kr_splitmix64(b->x);
	k = b->x%(b->m[hn].n);
	if(b->m[hn].a[k] == b->s[id].s && b->m[hn].n > 1)
	{
		k = b->x%(b->m[hn].n-1);
		if(b->m[hn].a[k] == b->s[id].s) k = b->m[hn].n-1;
	} 
	return b->m[hn].a[k];
}

static void mcgg_perturb(const mc_opt_t *opt, const mc_match_t *ma, mcgg_svaux_t *b)
{
    uint32_t i;
    for (i = 0; i < b->cc_size; ++i) {
        uint32_t k = (uint32_t)ma->cc[b->cc_off + i];///node id
        double y;
        y = kr_drand_r(&b->x);
        if (y < opt->f_perturb)
			mcgg_set_spin(ma, b, k, kr_drand_node_new(k, b), __func__);
    }
}

static uint32_t mcgg_bfs(const mc_match_t *ma, mcgg_svaux_t *b, uint32_t k0, uint32_t bfs_round, uint32_t max_size)
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
static void mcgg_perturb_node(const mc_opt_t *opt, const mc_match_t *ma, mcgg_svaux_t *b, int32_t bfs_round)
{
	uint32_t i, k, n_bfs = 0;
	k = (uint32_t)(kr_drand_r(&b->x) * b->cc_size + .499);
	if(k >= b->cc_size) k = b->cc_size - 1;
	k = (uint32_t)ma->cc[b->cc_off + k];///node id
	n_bfs = mcgg_bfs(ma, b, k, bfs_round, (int32_t)(b->cc_size * opt->f_perturb));
	for (i = 0; i < n_bfs; ++i)
		mcgg_set_spin(ma, b, b->bfs[i], kr_drand_node_new(b->bfs[i], b), __func__);
}

void print_mcgg_node(const mc_match_t *ma, mcgg_svaux_t *b, uint32_t id)
{
	fprintf(stderr, "[M::%s::utg%.6ul-hap%u]\n", __func__, id, b->s[id].s);
	w_t w[128];
	uint32_t o, n, i, hn = (1<<b->z->hapN);
	for (i = 0; i < hn; i++) w[i] = 0;
	o = ma->idx.a[id] >> 32;
    n = (uint32_t)ma->idx.a[id];
	for (i = 0; i < n; ++i) w[b->s[ma_y(ma->ma.a[o + i])].s] += ma->ma.a[o + i].w;
	for (i = 0; i < hn; i++) fprintf(stderr, "w[%u]-%f, z[%u]-%f\n", i, w[i], i, mcp_de(*(b->z), id, i));
}

uint32_t mcgg_solve_cc(const mc_opt_t *opt, const mc_gg_t *mg, mcgg_svaux_t *b, uint32_t cc_off, uint32_t cc_size)
{
	// double t0, t1, tt0, tt1;
	uint32_t j, k, n_iter = 0, flush = opt->max_iter * 50;
    t_w_t sc_opt = -(1<<30), sc;///problem-w
    b->cc_off = cc_off, b->cc_size = cc_size;
    if (b->cc_size < 2) return 0;
	sc_opt = mcgg_init_spin(mg->e, b);
	if (b->cc_size == 2) return 0;
	for (j = 0; j < b->cc_size; ++j) back_status(b, b->cc_node[j], 1);

	sc = mcgg_optimize_local(opt, mg->e, b, &n_iter);
	if (sc > sc_opt)
	{
		for (j = 0; j < b->cc_size; ++j) back_status(b, b->cc_node[j], 1);
		sc_opt = sc;
	}
	else
	{
		for (j = 0; j < b->cc_size; ++j) back_status(b, b->cc_node[j], 0);
	}
	fprintf(stderr, "\nBeg-[M::%s::score->%f]\n", __func__, mcgg_score(mg->e, b));
	// print_mcgg_node(mg->e, b, 3838);
	// print_mcgg_node(mg->e, b, 36880);

	// tt0 = tt1 = 0;
	for (k = 0; k < (uint32_t)opt->n_perturb; ++k) {
		// t0 = yak_realtime();
        if (k&1) mcgg_perturb(opt, mg->e, b);
        else mcgg_perturb_node(opt, mg->e, b, 3);
		// tt0 += yak_realtime()-t0;

		// t1 = yak_realtime();
        sc = mcgg_optimize_local(opt, mg->e, b, &n_iter);
		// tt1 += yak_realtime()-t1;
		// fprintf(stderr, "++(%u) sc_after: %f\n", k, mcgg_score(mg->e, b));

        if (sc > sc_opt) {
            for (j = 0; j < b->cc_size; ++j) back_status(b, b->cc_node[j], 1);
			sc_opt = sc;
        } else {
            for (j = 0; j < b->cc_size; ++j) back_status(b, b->cc_node[j], 0);
        }

		
        if((n_iter%flush) == 0)
        {
            mcgg_reset_z(mg->e, b);
            sc = mcgg_score(mg->e, b);
            for (j = 0; j < b->cc_size; ++j) back_status(b, b->cc_node[j], 1);
			sc_opt = sc;
        }
		// fprintf(stderr, "\n");
		// print_mcgg_node(mg->e, b, 3838);
		// print_mcgg_node(mg->e, b, 36880);
		// if((k&31)==0)fprintf(stderr, "+++(%u) sc: %f, sc_opt: %f, tt0: %.3f, tt1: %.3f\n", k, sc, sc_opt, tt0, tt1);
    }

	for (j = 0; j < b->cc_size; ++j) back_status(b, b->cc_node[j], 0);
	fprintf(stderr, "End-[M::%s::score->%f]\n", __func__, mcgg_score(mg->e, b));
	return n_iter;
}

void mc_solve_core_genral(const mc_opt_t *opt, mc_gg_t *mg, uint32_t hapN)
{
	double index_time = yak_realtime();
	uint32_t st, i;
	mcgg_svaux_t *b;
	mc_g_cc(mg->e);
	b = mcgg_svaux_init(mg, opt->seed, hapN);
	// if(VERBOSE_CUT)
	// {
	// 	fprintf(stderr, "\n\n\n\n\n*************beg-[M::%s::score->%f] ==> Partition\n", __func__, mc_score_all_advance(mg->e, mg->s.a));
	// }
	
	for (st = 0, i = 1; i <= mg->e->n_seq; ++i) {
		if (i == mg->e->n_seq || mg->e->cc[st]>>32 != mg->e->cc[i]>>32) {
			mcgg_solve_cc(opt, mg, b, st, i - st);
			st = i;
		}
	}

	// if(VERBOSE_CUT)
	// {
	// 	fprintf(stderr, "##############end-[---M::%s::score->%f] ==> Partition\n", __func__, mc_score_all(mg->e, b));
	// }
	///mc_write_info(g, b);
	mcgg_svaux_destroy(b);
	fprintf(stderr, "[M::%s::%.3f] ==> Partition\n", __func__, yak_realtime()-index_time);
}

void print_mcb(mc_gg_t *mg)
{
	uint32_t i, k, m;
	mcg_node_t t;
	mcb_t *p;
	for (i = 0; i < mg->m.n; i++)
	{
		p = mcb_pat(*mg, i + 1);
		fprintf(stderr, "# haplotypes: %u, # combination: %u\n", i+1, (uint32_t)p->n);
		for (k = 0; k < p->n; k++)
		{
			t = p->a[k];
			for (m = 0; m < 32; m++)
			{
				if((t>>m)&1) fprintf(stderr, "%u\t", m);
			}
			fprintf(stderr, "\n");
		}
	}
}

void print_hap_p(kv_gg_status *s)
{
	fprintf(stderr, "\n[M::%s]\n", __func__);
	uint32_t i, h;
	mcg_node_t m;
	for (i = 0; i < s->n; i++)
	{
		fprintf(stderr,"utg%.6ul\t", i + 1);
		m = s->a[i].s; h = 0;
		while (m) {
			h++;
			if(m&1) fprintf(stderr, "h%u\t", h);
			m>>=1;
		}
		fprintf(stderr,"\n");
	}
}

void write_mc_gg_dump(kv_u_trans_t *ta, uint32_t un, kv_gg_status *s, uint16_t hapN, const char* fn)
{
	fprintf(stderr, "\n[M::%s]\n", __func__);
	char *buf = (char*)calloc(strlen(fn) + 50, 1);
    sprintf(buf, "%s.hic.dbg.dump.bin", fn);
    FILE* fp = fopen(buf, "w");

	fwrite(&ta->n, sizeof(ta->n), 1, fp);
    fwrite(ta->a, sizeof(u_trans_t), ta->n, fp);

    fwrite(&ta->idx.n, sizeof(ta->idx.n), 1, fp);
    fwrite(ta->idx.a, sizeof(uint64_t), ta->idx.n, fp);

	fwrite(&un, sizeof(un), 1, fp);

	fwrite(&(s->n), sizeof(s->n), 1, fp);
    fwrite(s->a, sizeof(mc_gg_status), s->n, fp);

	fwrite(&hapN, sizeof(hapN), 1, fp);

	fclose(fp);
    free(buf);
}

uint32_t load_mc_gg_dump(kv_u_trans_t **rta, uint32_t *un, kv_gg_status **rs, uint16_t *hapN, const char* fn)
{
	fprintf(stderr, "\n[M::%s]\n", __func__);
	kv_u_trans_t *ta = NULL;
	kv_gg_status *s = NULL;
	uint64_t flag = 0;
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
    sprintf(buf, "%s.hic.dbg.dump.bin", fn);

    FILE* fp = NULL; 
    fp = fopen(buf, "r"); 
    if(!fp)
	{
		free(buf);
		return 0;
	} 
	CALLOC(ta, 1);
	flag += fread(&ta->n, sizeof(ta->n), 1, fp);
	ta->m = ta->n; MALLOC(ta->a, ta->n);
	flag += fread(ta->a, sizeof(u_trans_t), ta->n, fp);

	flag += fread(&ta->idx.n, sizeof(ta->idx.n), 1, fp);
	ta->idx.m = ta->idx.n; MALLOC(ta->idx.a, ta->idx.n);
    flag += fread(ta->idx.a, sizeof(uint64_t), ta->idx.n, fp);

	flag += fread(un, sizeof(*un), 1, fp);

	CALLOC(s, 1);
	flag += fread(&(s->n), sizeof(s->n), 1, fp);
	s->m = s->n; MALLOC(s->a, s->n);
    flag += fread(s->a, sizeof(mc_gg_status), s->n, fp);

	flag += fread(hapN, sizeof(*hapN), 1, fp);
	*rta = ta; *rs = s;
	fclose(fp);
    free(buf);
	return 1;
}

mc_g_t* to_mc_g_t(kv_u_trans_t *ta, kv_gg_status *s, uint32_t un)
{
	fprintf(stderr, "[M::%s]\n", __func__);
	uint32_t i, k;
	mc_edge_t *ma = NULL;
	mc_g_t *mg = NULL; CALLOC(mg, 1);
    

	CALLOC(mg->e, 1);
	mg->e->n_seq = un;
	kv_init(mg->e->idx); kv_init(mg->e->ma); 
	// double sc = get_w_scale(ta);
	// fprintf(stderr, "sc: %f\n", sc);
	for (i = 0; i < ta->n; ++i)
	{
		if(ta->a[i].del) continue;
		kv_pushp(mc_edge_t, mg->e->ma, &ma);
		ma->x = (uint64_t)ta->a[i].qn << 32 | ta->a[i].tn;
		// ma->w = w_cast((ta->a[i].nw*sc));
		ma->w = w_cast((ta->a[i].nw));
	}
	

	for (i = k = 0; i < mg->e->ma.n; i++)
    {
		if(mg->e->ma.a[i].w == 0) continue;
		mg->e->ma.a[k] = mg->e->ma.a[i];
        k++;
    }
    mg->e->ma.n = k;
	
	radix_sort_mce(mg->e->ma.a, mg->e->ma.a + mg->e->ma.n);
	mc_merge_dup(mg);
	mc_edges_idx(mg->e);
	mc_edges_symm(mg->e);

	kv_init(mg->s); 
	mg->s.m = mg->s.n = s->n; CALLOC(mg->s.a, mg->s.n);
	for (i = 0; i < mg->s.n; ++i)
	{
		if(s->a[i].s != 1 && s->a[i].s != 2) continue;
		mg->s.a[i] = (s->a[i].s == 1? 1:-1);
	}
	return mg;
}

void debug_mc_gg_t(const char* fn, uint32_t update_ta, uint32_t convert_mc_g_t)
{
	kv_u_trans_t *ta = NULL;
	kv_gg_status *s = NULL;
	uint32_t un;
	uint16_t hapN;
	if(load_mc_gg_dump(&ta, &un, &s, &hapN, fn))
	{
		if(convert_mc_g_t)
		{
			mc_opt_t opt;
			mc_opt_init(&opt, asm_opt.n_perturb, asm_opt.f_perturb, asm_opt.seed);
			mc_g_t *mg = NULL;
			mg = to_mc_g_t(ta, s, un);
			char *o_file = get_outfile_name(asm_opt.output_file_name);
			trans_chain* t_ch = load_hc_trans(o_file);
			mb_solve_core(&opt, mg, &(t_ch->k_trans), 1);
			mc_solve_core(&opt, mg, NULL);
		}
		else
		{
			mc_solve_general(ta, un, s, hapN, update_ta, 0);
		}
	}
	exit(1);
}

void clean_solve_general_ovlp(kv_u_trans_t *ta, uint32_t un, kv_gg_status *s)
{
	uint32_t i;
	for (i = 0; i < ta->n; i++) ta->a[i].del = !!(s->a[ta->a[i].qn].s&s->a[ta->a[i].tn].s);
	kt_u_trans_t_simple_symm(ta, un, 0);
}

void mc_solve_general(kv_u_trans_t *ta, uint32_t un, kv_gg_status *s, uint16_t hapN, uint16_t update_ta, uint16_t write_dump)
{
	print_hap_p(s);
	exit(1);
	mc_opt_t opt;
	mc_opt_init(&opt, asm_opt.n_perturb, asm_opt.f_perturb, asm_opt.seed);
	mc_gg_t *mg = init_mc_gg_t(un, s, hapN);
	// print_mcb(mg);
	
	update_mc_edges_general(mg, ta, hapN);
	mc_solve_core_genral(&opt, mg, hapN);
	if(update_ta) clean_solve_general_ovlp(ta, un, s);
	// print_hap_p(s);
	
	destory_mc_gg_t(&mg);
	if(write_dump) write_mc_gg_dump(ta, un, s, hapN, MC_NAME);
}