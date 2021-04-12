#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <stdlib.h>
///#include "partig.h"
#include "rcut.h"
#include "Purge_Dups.h"
#include "Correct.h"
#include "ksort.h"

#define mc_edge_key(e) ((e).x)
KRADIX_SORT_INIT(mce, mc_edge_t, mc_edge_key, member_size(mc_edge_t, x))
#define mc_generic_key(x) (x)
KRADIX_SORT_INIT(mc64, uint64_t, mc_generic_key, 8)

#define pt_a(x, id) ((x).ma.a + ((x).idx.a[(id)]>>32))
#define pt_n(x, id) ((uint32_t)((x).idx.a[(id)]))
#define ma_x(z) (((z).x>>32))
#define ma_y(z) (((uint32_t)((z).x)))


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
	///uint32_t n_cc_edge, m_cc_edge;
	///uint64_t *cc_edge;
	kvec_t(uint64_t) cc_edge;
	uint32_t n_sub;
	uint32_t *cc_node;
	uint32_t *bfs, *bfs_mark;
	mc_pairsc_t *z, *z_opt;///keep scores to nodes(1) and nodes(-1)
	int8_t *s, *s_opt;
} mc_svaux_t;

void mc_opt_init(mc_opt_t *opt)
{
	memset(opt, 0, sizeof(mc_opt_t));
	opt->n_perturb = 5000;
	opt->f_perturb = 0.1;
	opt->max_iter = 1000;
	opt->seed = 11;
}

mc_g_t *init_mc_g_t(ma_ug_t *ug, asg_t *read_g)
{
	mc_g_t *p = NULL; CALLOC(p, 1);
	p->ug = ug;
	p->rg = read_g;
	kv_init(p->s);
	CALLOC(p->s.a, ug->g->n_seq);
	p->s.n = p->s.m = ug->g->n_seq;
	return p;
}

void destory_mc_g_t(mc_g_t **p)
{
	if(!p || !(*p)) return;
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


static inline uint64_t kr_splitmix64(uint64_t x)
{
	uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
	return z ^ (z >> 31);
}

static inline double kr_drand_r(uint64_t *x)
{
    union { uint64_t i; double d; } u;
	*x = kr_splitmix64(*x);
    u.i = 0x3FFULL << 52 | (*x) >> 12;
    return u.d - 1.0;
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

static mc_edge_t *get_mc_edge(const mc_match_t *ma, uint32_t sid1, uint32_t sid2)
{
	mc_edge_t *o = pt_a(*ma, sid1);
	uint32_t n = pt_n(*ma, sid1), k;
	for (k = 0; k < n; ++k)
		if (((uint32_t)o[k].x) == sid2)
			return &(o[k]);
	return NULL;
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

void update_mc_edges(mc_g_t *mg, hap_overlaps_list* ha, kv_u_trans_t *ta, trans_chain* t_ch, double f_rate)
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

		debug_mc_interval_t(p.a, p.n, p_idx.a, mg->ug, mg->rg, t_ch);
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
				kv_pushp(mc_edge_t, mg->e->ma, &ma);
                ma->x = (uint64_t)ta->a[i].qn << 32 | ta->a[i].tn;
                ma->w = ta->a[i].nw;
			}
		}
	}
	
	radix_sort_mce(mg->e->ma.a, mg->e->ma.a + mg->e->ma.n);
	mc_merge_dup(mg);
	mc_edges_idx(mg->e);
	mc_edges_symm(mg->e);
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
	CALLOC(b->z, ma->n_seq);
	CALLOC(b->z_opt, ma->n_seq);
	return b;
}
void mc_svaux_destroy(mc_svaux_t *b)
{
	b->s = NULL;
	kv_destroy(b->cc_edge); free(b->cc_node);
	free(b->s); free(b->s_opt);
	free(b->z); free(b->z_opt);
	free(b->bfs); free(b->bfs_mark);
	free(b);
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

t_w_t mc_init_spin(const mc_opt_t *opt, const mc_match_t *ma, mc_svaux_t *b)
{
	uint32_t i;
	b->cc_edge.n = 0;
	for (i = 0; i < b->cc_size; ++i) {///how many nodes
		uint32_t k = (uint32_t)ma->cc[b->cc_off + i];///node id
		uint32_t o = ma->idx.a[k] >> 32;///cc group id
		uint32_t n = (uint32_t)ma->idx.a[k], j;
		b->cc_node[i] = k;
		for (j = 0; j < n; ++j) {
			w_t w = ma->ma.a[o + j].w;
			w = w > 0? w : -w;
			kv_push(uint64_t, b->cc_edge, (uint64_t)((uint32_t)-1 - ((uint32_t)w)) << 32 | (o + j));
			///b->cc_edge[b->n_cc_edge++] = (uint64_t)((uint32_t)-1 - w) << 32 | (o + j);
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

static t_w_t mc_optimize_local(const mc_opt_t *opt, const mc_match_t *ma, mc_svaux_t *b, uint32_t *n_iter)
{
	int32_t n_iter_local = 0;
	while (n_iter_local < opt->max_iter) {
		uint32_t i, n_flip = 0;
		++(*n_iter);
		ks_shuffle_uint32_t(b->cc_size, b->cc_node, &b->x);
		for (i = 0; i < b->cc_size; ++i) {
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

uint32_t mc_solve_cc(const mc_opt_t *opt, const mc_g_t *mg, mc_svaux_t *b, uint32_t cc_off, uint32_t cc_size)
{
	uint32_t j, k, n_iter = 0;
	t_w_t /**sc_ori,**/ sc_opt = -(1<<30), sc;///problem-w

	b->cc_off = cc_off, b->cc_size = cc_size;
	if (b->cc_size < 2) return 0;

	// first guess
	/**sc_ori =**/ mc_init_spin(opt, mg->e, b);
	if (b->cc_size == 2) return 0;

	// optimize
	sc_opt = mc_optimize_local(opt, mg->e, b, &n_iter);
	for (j = 0; j < b->cc_size; ++j) {///backup s and z in s_opt and z_opt
		b->s_opt[b->cc_node[j]] = b->s[b->cc_node[j]]; ///hap status of each unitig
		b->z_opt[b->cc_node[j]] = b->z[b->cc_node[j]]; ///z[0]: positive weight; z[1]: positive weight
	}
	for (k = 0; k < (uint32_t)opt->n_perturb; ++k) {
		if (k&1) mc_perturb(opt, mg->e, b);
		else mc_perturb_node(opt, mg->e, b, 3);
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
	}
	for (j = 0; j < b->cc_size; ++j)
		b->s[b->cc_node[j]] = b->s_opt[b->cc_node[j]];
	return n_iter;
}

void mc_solve_core(const mc_opt_t *opt, mc_g_t *mg)
{
	double index_time = yak_realtime();
	uint32_t st, i;
	mc_svaux_t *b;
	mc_g_cc(mg->e);
	b = mc_svaux_init(mg, opt->seed);
	for (st = 0, i = 1; i <= mg->e->n_seq; ++i) {
		if (i == mg->e->n_seq || mg->e->cc[st]>>32 != mg->e->cc[i]>>32) {
			mc_solve_cc(opt, mg, b, st, i - st);
			st = i;
		}
	}
	///mc_write_info(g, b);
	mc_svaux_destroy(b);
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

void mc_solve(hap_overlaps_list* ovlp, trans_chain* t_ch, kv_u_trans_t *ta, ma_ug_t *ug, asg_t *read_g, double f_rate, uint8_t* trio_flag)
{
	mc_opt_t opt;
	mc_opt_init(&opt);
	mc_g_t *mg = init_mc_g_t(ug, read_g);
	update_mc_edges(mg, ovlp, ta, t_ch, f_rate);
	debug_mc_g_t(mg);
	mc_solve_core(&opt, mg);

	if((asm_opt.flag & HA_F_PARTITION) && t_ch)
	{
		p_nodes(mg, t_ch, trio_flag);
	}

	if(ovlp) clean_ovlp_by_mc(mg, ovlp);


	destory_mc_g_t(&mg);
}