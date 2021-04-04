#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <stdlib.h>
#include "partig.h"
#include "Purge_Dups.h"
#include "Correct.h"
#include "ksort.h"

#define generic_key(x) (x)
KRADIX_SORT_INIT(gfa64, uint64_t, generic_key, 8)

#define pt_a(x, id) ((x).ma.a + ((x).idx.a[(id)]>>32))
#define pt_n(x, id) ((uint32_t)((x).idx.a[(id)]))

typedef struct {
	int32_t topn;
	int32_t n_perturb;
	uint64_t seed;
	double f_perturb;
} pt_svopt_t;

typedef struct {
	///uint32_t m, n, *shuffled;
	uint32_t *shuffled;
	uint32_t off, size; // offset in pt_match_t::cc; size of the component
	uint64_t *buf;
	kvec_t(uint64_t) a;
	int8_t *s, *s_tmp;///s is the status (haplotype) of each unitig: for backup
} solve_aux_t;

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

void pt_svopt_init(pt_svopt_t *opt)
{
	memset(opt, 0, sizeof(pt_svopt_t));
	opt->seed = 11;
	opt->topn = 1<<30;
	opt->n_perturb = 1000;
	opt->f_perturb = 0.1;
}

static void pt_pdist_idx(pt_match_t *ma)
{
	uint32_t st, i;
	kv_resize(uint64_t, ma->idx, ma->n_seq);
	ma->idx.n = ma->n_seq;
	memset(ma->idx.a, 0, ma->idx.n*sizeof(uint64_t));
	for (st = 0, i = 1; i <= ma->ma.n; ++i)
		if (i == ma->ma.n || ma->ma.a[i].sid[0] != ma->ma.a[st].sid[0])
			ma->idx.a[ma->ma.a[st].sid[0]] = (uint64_t)st << 32 | (i - st), st = i;
}

static pt_match1_t *pt_pdist(const pt_match_t *ma, uint32_t sid1, uint32_t sid2)
{
	pt_match1_t *o = pt_a(*ma, sid1);
	uint32_t n = pt_n(*ma, sid1), k;
	for (k = 0; k < n; ++k)
		if (o[k].sid[1] == sid2)
			return &(o[k]);
	return NULL;
}

static void normalize_pdist(pt_match1_t *a, pt_match1_t *b)
{
	if(a->w >= b->w)
	{
		b->sid[0] = a->sid[1];
		b->sid[1] = a->sid[0];
		b->w = a->w;
	}
	else
	{
		a->sid[0] = b->sid[1];
		a->sid[1] = b->sid[0];
		a->w = b->w;
	}
}

uint32_t pt_pdist_symm(pt_match_t *ma)
{
	uint8_t *del = NULL;
	uint32_t i, k, n = 0;
	pt_match1_t *t = NULL, *m = NULL;
	CALLOC(del, ma->ma.n);

	for (i = 0; i < ma->ma.n; ++i) {
		m = &ma->ma.a[i];
		if (m->sid[0] == m->sid[1])
		{
			del[i] = 1, ++n;///self overlap
			continue;
		}
		t = pt_pdist(ma, m->sid[1], m->sid[0]);
		if(!t)
		{
			del[i] = 1, ++n;///self overlap
			continue;
		} 
		normalize_pdist(m, t);
	}

	if (n > 0) {
		for (i = k = 0; i < ma->ma.n; ++i)
			if (!del[i]) ma->ma.a[k++] = ma->ma.a[i];
		ma->ma.n = k;
		pt_pdist_idx(ma);
	}

	free(del);
	return n;
}


static void pt_pdist_flt(pt_match_t *ma, uint32_t min_cnt, double drop_thres)
{
	uint32_t i, k, j, n, max, o;
	uint8_t *mark = NULL;
	CALLOC(mark, ma->ma.n);
	for (i = 0; i < ma->n_seq; ++i) {
		o = ma->idx.a[i] >> 32;
		n = (uint32_t)ma->idx.a[i];
		max = 0;
		if (n == 0) continue;
		for (j = o; j < o + n; ++j)
			max = max > ma->ma.a[j].w? max : ma->ma.a[j].w;
		for (j = o; j < o + n; ++j)
			if (ma->ma.a[j].w >= (max*drop_thres) || ma->ma.a[j].w + min_cnt >= max)
				mark[j] = 1;
	}
	for (i = 0; i < ma->ma.n; ++i)
	{
		if (mark[i] == 0) continue;
		o = ma->idx.a[ma->ma.a[i].sid[1]]>>32;
		n = (uint32_t)ma->idx.a[ma->ma.a[i].sid[1]];
		for (j = o; j < o + n; ++j)
		{
			if (ma->ma.a[j].sid[1] == ma->ma.a[i].sid[0]) mark[j] = 1;
		}
	}
		
	for (i = k = 0; i < ma->ma.n; ++i)
		if (mark[i]) ma->ma.a[k++] = ma->ma.a[i];
	ma->ma.n = k;
	free(mark);
	pt_pdist_idx(ma);
	pt_pdist_symm(ma);
}

pt_match_t *init_pt_match_t(hap_overlaps_list* ha, pt_g_t *x, double f_rate)
{
	pt_match_t *p = NULL; CALLOC(p, 1); p->n_seq = x->ug->g->n_seq;
	kv_init(p->idx); kv_init(p->ma);
	uint32_t v, i, k, qn, tn, qs, qe, ts, te, occ, as, ae;
	uint64_t hetLen, homLen, oLen;
	pt_node_t *a = NULL;
	pt_match1_t *ma = NULL;

	for (v = 0; v < ha->num; v++)
    {
        for (i = 0; i < ha->x[v].a.n; i++)
        {
			if(ha->x[v].a.a[i].score <= 0) continue;
			if(ha->x[v].a.a[i].xUid == ha->x[v].a.a[i].yUid) continue;
            /*****************qn*****************/
            qn = ha->x[v].a.a[i].xUid;
            qs = ha->x[v].a.a[i].x_beg_pos;
            qe = ha->x[v].a.a[i].x_end_pos - 1;

			a = x->p.a + x->p_idx.a[qn]; 
			occ = x->p_idx.a[qn+1] - x->p_idx.a[qn];
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

			a = x->p.a + x->p_idx.a[tn]; 
			occ = x->p_idx.a[tn+1] - x->p_idx.a[tn];
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

			kv_pushp(pt_match1_t, p->ma, &ma);
			ma->sid[0] = ha->x[v].a.a[i].xUid;
			ma->sid[1] = ha->x[v].a.a[i].yUid;
			ma->w = ha->x[v].a.a[i].score;
        }
    }

	pt_pdist_idx(p);
	pt_pdist_symm(p);
	pt_pdist_flt(p, 5, 0.5);
	return p;
}


void debug_pt_g_t(pt_g_t *pg)
{
    fprintf(stderr, "0----------[M::%s]----------\n", __func__);
    uint32_t i, offset, v, sid, eid, spos, epos, p_status, p_uid, occ;
    ma_utg_t *u = NULL;
    pt_node_t *a = NULL, *t = NULL;

    for (v = 0; v < pg->ug->u.n; v++)
    {
		a = pg->p.a + pg->p_idx.a[v]; 
		occ = pg->p_idx.a[v+1] - pg->p_idx.a[v];
        for (i = 0; i < occ; i++)
        {
            if(a[i].uID != v) fprintf(stderr, "ERROR-s\n");
        }
    }

    for (v = 0, p_status = (uint32_t)-1, p_uid = (uint32_t)-1; v < pg->p.n; v++)
    {
        t = &(pg->p.a[v]);
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
        u = &(pg->ug->u.a[t->uID]);
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
                if(epos != (offset+pg->rg->seq[u->a[i]>>33].len - 1))
                {
                    fprintf(stderr, "ERROR-c, real end: %u\n", 
                                    (uint32_t)(offset+pg->rg->seq[u->a[i]>>33].len - 1));
                }
            }

            offset += (uint32_t)u->a[i];
            if(i >= sid && i <= eid)
            {
                if(pg->t_ch->is_r_het[u->a[i]>>33] != t->hs)
                {
                    fprintf(stderr, "ERROR-d: is_r_het: %u, h_status: %u\n", pg->t_ch->is_r_het[u->a[i]>>33], t->hs);
                }
            }
        }
    }


	pt_match1_t *o = NULL, *s = NULL;
	uint32_t k, n, cnt;
	for (i = 0; i < pg->e->n_seq; ++i) 
	{
		o = pt_a(*(pg->e), i); n = pt_n(*(pg->e), i);
		for (k = 0; k < n; ++k)
		{
			if(o[k].sid[0] != i) fprintf(stderr, "ERROR-g\n");
			s = pt_pdist(pg->e, o[k].sid[1], o[k].sid[0]);
			if(!s) fprintf(stderr, "ERROR-e\n");
			if(s)
			{
				if(!(s->sid[0] == o[k].sid[1] && s->sid[1] == o[k].sid[0] && s->w == o[k].w))
				{
					fprintf(stderr, "ERROR-f\n");
				}
			}
		}

		for (k = cnt = 0; k < pg->e->ma.n; ++k)
		{
			if(pg->e->ma.a[k].sid[0] == i) cnt++;
		}

		if(cnt != n) fprintf(stderr, "ERROR-h\n");
	}
}

pt_g_t *init_pt_g_t(hap_overlaps_list* ovlp, trans_chain* t_ch, ma_ug_t *ug, asg_t *read_g, double f_rate)
{	
	uint32_t v, l, k, offset, l_pos;
	pt_g_t *p = NULL; CALLOC(p, 1);
	asg_t* nsg = ug->g;
	ma_utg_t *u = NULL;
	pt_node_t *t = NULL;
	p->ug = ug;
	p->rg = read_g;
	p->t_ch = t_ch;
	kv_init(p->info); p->info.n = p->info.m = p->ug->g->n_seq; CALLOC(p->info.a, p->info.n);
	kv_init(p->p);
	kv_init(p->p_idx); kv_push(uint32_t, p->p_idx, 0);

	for (v = 0; v < nsg->n_seq; v++)
    {
        if(nsg->seq[v].del || nsg->seq[v].c == ALTER_LABLE) continue;

        u = &(ug->u.a[v]);
        for (k = 1, l = 0, offset = 0, l_pos = 0; k <= u->n; ++k) 
        {   
            if (k == u->n || t_ch->is_r_het[u->a[k]>>33] != t_ch->is_r_het[u->a[l]>>33])
            {
                kv_pushp(pt_node_t, p->p, &t);
				t->uID = v;
				t->hs = t_ch->is_r_het[u->a[l]>>33];

				t->bS = l_pos;
				t->bE = offset + read_g->seq[u->a[k-1]>>33].len - 1;

				t->nS = l;
				t->nE = k - 1;
                l = k;
                l_pos = offset + (uint32_t)u->a[k-1];
            }
            offset += (uint32_t)u->a[k-1];
        }
		kv_push(uint32_t, p->p_idx, p->p.n);
    }

	p->e = init_pt_match_t(ovlp, p, f_rate);
	return p;
}

void destory_pt_g_t(pt_g_t **p)
{
	if(!p || !(*p)) return;
	kv_destroy((*p)->p);
	kv_destroy((*p)->info);
	kv_destroy((*p)->p_idx);
	kv_destroy((*p)->e->idx);
	kv_destroy((*p)->e->ma);
	free((*p)->e->cc);
	free((*p)->e);
	free((*p));
}

uint64_t *pt_cc_core(const pt_match_t *ma)
{
	uint32_t i, x, y, *flag;
	uint64_t *group;
	pt_match1_t *o = NULL;
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
				uint32_t t = o[j].sid[1];
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
	radix_sort_gfa64(group, group + ma->n_seq);
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

void pt_cc(pt_match_t *ma)
{
	ma->cc = pt_cc_core(ma);
}

///inspect top INT edges
static int64_t pt_score(const pt_match_t *ma, uint32_t topn, solve_aux_t *aux)
{
	uint32_t i;
	int64_t z = 0;
	for (i = 0; i < aux->size; ++i) {///aux->size: how many unitigs in this group
		uint32_t k = (uint32_t)ma->cc[aux->off + i];
		uint32_t o = ma->idx.a[k] >> 32;
		uint32_t n = (uint32_t)ma->idx.a[k], j;
		for (j = 0; j < n; ++j)
			aux->buf[j] = (uint64_t)((uint32_t)-1 - ma->ma.a[o + j].w) << 32 | (o + j);
		radix_sort_gfa64(aux->buf, aux->buf + n);
		for (j = 0; j < n && j < topn; ++j) {
			const pt_match1_t *m = &ma->ma.a[(uint32_t)aux->buf[j]];
			z += -(int64_t)m->w * aux->s[m->sid[0]] * aux->s[m->sid[1]];
		}
	}
	return z;
}

static int64_t pt_solve1_init_phase(const pt_match_t *ma, int32_t topn, uint64_t *x, solve_aux_t *aux)
{
	uint32_t i;
	aux->a.n = 0;
	for (i = 0; i < aux->size; ++i) {
		uint32_t k = (uint32_t)ma->cc[aux->off + i];///unitig id
		uint32_t o = ma->idx.a[k] >> 32;///group id
		uint32_t n = (uint32_t)ma->idx.a[k], j;
		aux->shuffled[i] = k;///init shuffled status
		for (j = 0; j < n; ++j) {
			const pt_match1_t *m = &ma->ma.a[o + j];
			///first is weight, second part is edge id
			kv_push(uint64_t, aux->a, (uint64_t)((uint32_t)-1 - m->w) << 32 | (o + j));
		}
	}
	radix_sort_gfa64(aux->a.a, aux->a.a + aux->a.n);///sort all edges in this group
	///randomly assign each unitig
	for (i = 0; i < aux->a.n; ++i) { // from the strongest edge to the weakest
		const pt_match1_t *m = &ma->ma.a[(uint32_t)aux->a.a[i]];
		///sid[0]: query id, sid[1]: target id
		///the initial results of aux->s is 0
		if (aux->s[m->sid[0]] == 0 && aux->s[m->sid[1]] == 0) {
			*x = kr_splitmix64(*x);// random number
			aux->s[m->sid[0]] = *x&1? 1 : -1;
			aux->s[m->sid[1]] = -aux->s[m->sid[0]];
		} else if (aux->s[m->sid[0]] == 0) {
			aux->s[m->sid[0]] = -aux->s[m->sid[1]];
		} else if (aux->s[m->sid[1]] == 0) {
			aux->s[m->sid[1]] = -aux->s[m->sid[0]];
		}
	}
	return pt_score(ma, topn, aux);
}

///size is how many unitigs in this group
static void ks_shuffle_uint32_t(size_t n, uint32_t a[], uint64_t *x)
{
	size_t i, j;
	for (i = n; i > 1; --i) {
		uint32_t tmp;
		j = (size_t)(kr_drand_r(x) * i);///semms 
		tmp = a[j]; a[j] = a[i-1]; a[i-1] = tmp;
	}
}

static void pt_solve1_perturb(const pt_svopt_t *opt, const pt_match_t *ma, uint64_t *x, solve_aux_t *aux)
{
	uint32_t i;
	double y;
	for (i = 0; i < aux->size; ++i) {
		uint32_t k = (uint32_t)ma->cc[aux->off + i];
		y = kr_drand_r(x);
		if (y < opt->f_perturb)
			aux->s[k] = -aux->s[k];
	}
}

static int64_t pt_solve1_optimize(const pt_match_t *ma, uint32_t topn, uint64_t *x, solve_aux_t *aux, uint32_t *n_iter)
{
	uint32_t i;
	while (1) {
		uint32_t n_flip = 0;
		++(*n_iter);
		ks_shuffle_uint32_t(aux->size, aux->shuffled, x);
		for (i = 0; i < aux->size; ++i) {
			uint32_t k = aux->shuffled[i];
			uint32_t o = ma->idx.a[k] >> 32;
			uint32_t n = (uint32_t)ma->idx.a[k], j;
			uint64_t z[2];
			int8_t s;
			for (j = 0; j < n; ++j) {
				const pt_match1_t *m = &ma->ma.a[o + j];
				///assert(m->sid[0] == k);
				aux->buf[j] = (uint64_t)((uint32_t)-1 - m->w) << 32 | (o + j);
			}
			radix_sort_gfa64(aux->buf, aux->buf + n);///still sort by edge weight
			for (j = 0, z[0] = z[1] = 0; j < n && j < topn; ++j) {
				const pt_match1_t *m = &ma->ma.a[(uint32_t)aux->buf[j]];
				if (aux->s[m->sid[1]] > 0) z[0] += m->w;
				else if (aux->s[m->sid[1]] < 0) z[1] += m->w;
			}
			if (z[0] == z[1]) continue;
			s = z[0] > z[1]? -1 : 1;
			if (aux->s[k] != s)
				aux->s[k] = s, ++n_flip;
		}
		if (n_flip == 0) break;
	}
	return pt_score(ma, topn, aux);
}

uint32_t pt_solve1(const pt_svopt_t *opt, const pt_match_t *ma, uint64_t *x, solve_aux_t *aux)
{
	uint32_t j, k, n_iter = 0;
	int64_t sc_ori, sc_opt = -(1<<30), sc;
	if (aux->size < 2) return 0;///how many unitigs

	// first guess
	///randomly assign haplotype status, and get a score
	sc_ori = pt_solve1_init_phase(ma, opt->topn, x, aux);
	if (aux->size == 2) return 0;

	// optimize
	sc_opt = pt_solve1_optimize(ma, opt->topn, x, aux, &n_iter);
	for (j = 0; j < aux->size; ++j)
		aux->s_tmp[aux->shuffled[j]] = aux->s[aux->shuffled[j]];
	for (k = 0; k < (uint32_t)opt->n_perturb; ++k) {
		pt_solve1_perturb(opt, ma, x, aux);
		sc = pt_solve1_optimize(ma, opt->topn, x, aux, &n_iter);
		if (sc > sc_opt) {
			for (j = 0; j < aux->size; ++j)
				aux->s_tmp[aux->shuffled[j]] = aux->s[aux->shuffled[j]];
			sc_opt = sc;
		} else {
			for (j = 0; j < aux->size; ++j)
				aux->s[aux->shuffled[j]] = aux->s_tmp[aux->shuffled[j]];
		}
	}
	for (j = 0; j < aux->size; ++j)
		aux->s[aux->shuffled[j]] = aux->s_tmp[aux->shuffled[j]];
	fprintf(stderr, "[%s] group:%d, size:%d, #edges:%u, #iter:%d, sc_ori:%ld, sc_opt:%ld\n", __func__,
			(uint32_t)(ma->cc[aux->off]>>32), (uint32_t)(aux->size), (uint32_t)(aux->a.n), n_iter, (long)sc_ori, (long)sc_opt);
	return n_iter;
}

int8_t *pt_solve_core(const pt_svopt_t *opt, const pt_match_t *ma)
{
	int8_t *s;
	uint32_t st, i, max = 0;
	uint64_t x = opt->seed;
	solve_aux_t *aux;
	CALLOC(aux, 1); CALLOC(aux->s, ma->n_seq); CALLOC(aux->s_tmp, ma->n_seq);
	kv_init(aux->a);
	for (i = 0; i < ma->n_seq; ++i) {///count how many links for each unitig
		uint32_t n = pt_n(*ma, i);
		max = max > n? max : n;
	}
	MALLOC(aux->buf, max);
	MALLOC(aux->shuffled, ma->n_seq); // FIXME: this is over-allocation for convenience
	for (st = 0, i = 1; i <= ma->n_seq; ++i) {
		if (i == ma->n_seq || ma->cc[st]>>32 != ma->cc[i]>>32) {///at same group
			if (i - st >= 2) {///all unitigs in the same group
				aux->off = st, aux->size = i - st;
				pt_solve1(opt, ma, &x, aux);
			}
			st = i;
		}
	}
	s = aux->s;
	kv_destroy(aux->a); free(aux->buf); free(aux->shuffled); free(aux->s_tmp);
	free(aux);
	return s;
}

void set_trio_flag(ma_ug_t *ug, asg_t *read_g, uint32_t uID, uint8_t* trio_flag, trans_chain* t_ch, 
													hap_overlaps_list* ha, pt_match_t *ma, int8_t s)
{
	uint32_t i;
	ma_utg_t *u = &(ug->u.a[uID]);
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

	// pt_match1_t *o = pt_a(*ma, uID);
	// uint32_t n = pt_n(*ma, uID), k, qn, tn, qs, qe, r_beg, r_end, offset, oLen, found;
	// int index;
	// for (k = 0; k < n; ++k)
	// {
	// 	qn = o[k].sid[0]; tn = o[k].sid[1]; qs = qe = (uint32_t)-1;
	// 	index = get_specific_hap_overlap(&(ha->x[qn]), qn, tn);
	// 	if(index != -1 && ha->x[qn].a.a[index].score == (long long)o[k].w)
	// 	{
	// 		qs = ha->x[qn].a.a[index].x_beg_pos;
    //         qe = ha->x[qn].a.a[index].x_end_pos - 1;
	// 	}
	// 	else
	// 	{
	// 		index = get_specific_hap_overlap(&(ha->x[tn]), tn, qn);
	// 		if(index != -1 && ha->x[tn].a.a[index].score == (long long)o[k].w)
	// 		{
	// 			qs = ha->x[qn].a.a[index].y_beg_pos;
    //         	qe = ha->x[qn].a.a[index].y_end_pos - 1;
	// 		}
	// 	}
		
	// 	if(qs == (uint32_t)-1 || qe == (uint32_t)-1) fprintf(stderr, "ERROR\n");
	// 	for (i = 0, offset = 0, found = 0; i < u->n; i++)
	// 	{
	// 		r_beg = offset; r_end = offset + (long long)(read_g->seq[u->a[i]>>33].len) - 1;
    //     	offset += (uint32_t)u->a[i];
	// 		oLen = ((MIN(qe, r_end) >= MAX(qs, r_beg))? MIN(qe, r_end) - MAX(qs, r_beg) + 1 : 0);
	// 		if(oLen > 0) found = 1;
	// 		if(found == 1 && oLen == 0) break;

	// 		if(oLen > 0 && t_ch->is_r_het[u->a[i]>>33] != N_HET)
	// 		{
	// 			trio_flag[u->a[i]>>33] |= (s > 0? FATHER:MOTHER);
	// 		} 
			
	// 	}
	// }
}

void filter_ovlp(ma_ug_t *ug, asg_t *read_g, uint32_t uID, hap_overlaps_list* ha, pt_match_t *ma, 
int8_t *s)
{
	pt_match1_t *o = pt_a(*ma, uID);
	uint32_t n = pt_n(*ma, uID), k, qn, tn;
	hap_overlaps *p = NULL;
	int index;

	for (k = 0; k < n; ++k)
	{
		qn = o[k].sid[0]; tn = o[k].sid[1]; p = NULL;
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


void clean_ovlp(ma_ug_t *ug, asg_t *read_g, hap_overlaps_list* ha, pt_g_t *pg, int8_t* s)
{
	uint32_t v, i, k, qn, tn, types[4];
    types[X2Y] = Y2X; types[Y2X] = X2Y; types[XCY] = YCX; types[YCX] = XCY;
	int index;
	hap_overlaps *x = NULL, *y = NULL;
	for (i = 0; i < pg->e->n_seq; ++i) 
	{
		filter_ovlp(ug, read_g, i, ha, pg->e, s);
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


	// for (v = 0; v < ha->num; v++)
    // {
    //     for (i = k = 0; i < ha->x[v].a.n; i++)
    //     {
    //         if(ha->x[v].a.a[i].status != FLIP)
	// 		{
	// 			if(s[ha->x[v].a.a[i].xUid]*s[ha->x[v].a.a[i].yUid] == -1)
	// 			{
	// 				fprintf(stderr, "\ns[0]=%d, s[1]=%d\n", s[ha->x[v].a.a[i].xUid], s[ha->x[v].a.a[i].yUid]);
	// 				print_hap_paf(ug, &(ha->x[v].a.a[i]));
	// 			}
				
	// 			continue;
	// 		} 
			
			
    //         ha->x[v].a.a[k] = ha->x[v].a.a[i];
    //         k++;
    //     }
    //     ha->x[v].a.n = k;
    // }

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

void pt_solve(hap_overlaps_list* ovlp, trans_chain* t_ch, ma_ug_t *ug, asg_t *read_g, double f_rate, uint8_t* trio_flag)
{
	pt_svopt_t opt;
	int8_t *s = NULL;
	uint64_t *buf = NULL, i;
	pt_svopt_init(&opt);
	pt_g_t *pg = init_pt_g_t(ovlp, t_ch, ug, read_g, f_rate);
	///debug_pt_g_t(pg);

	pt_cc(pg->e);
	s = pt_solve_core(&opt, pg->e);

	if(asm_opt.flag & HA_F_PARTITION)
	{
		MALLOC(buf, pg->e->ma.n); // FIXME: this is over-allocation for convenience
		for (i = 0; i < pg->e->n_seq; ++i) {
			uint64_t z[2];
			uint32_t o = pg->e->idx.a[i] >> 32;
			uint32_t n = (uint32_t)pg->e->idx.a[i], j;

			set_trio_flag(ug, read_g, i, trio_flag, t_ch, ovlp, pg->e, s[i]);

			pg->info.a[i].s = s[i];
			for (j = 0; j < n; ++j) {
				const pt_match1_t *m = &pg->e->ma.a[o + j];
				buf[j] = (uint64_t)((uint32_t)-1 - m->w) << 32 | (o + j);
			}
			radix_sort_gfa64(buf, buf + n);
			for (j = 0, z[0] = z[1] = 0; j < n; ++j) {
				const pt_match1_t *m = &pg->e->ma.a[(uint32_t)buf[j]];
				if (s[m->sid[1]] > 0) z[0] += m->w;
				else if (s[m->sid[1]] < 0) z[1] += m->w;
			}
			pg->info.a[i].m[0] = z[0], pg->info.a[i].m[1] = z[1];
		}
	}
	

	clean_ovlp(ug, read_g, ovlp, pg, s);

	free(buf);
	free(s);
	destory_pt_g_t(&pg);
}