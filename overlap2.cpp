#include <assert.h>
#include "utils.h"
#include "CommandLines.h"
#include "Overlaps.h"

/*******************************
 * Dropping strong containment *
 *******************************/

static ma_hit_t *get_specific_overlap_with_del(ma_hit_t_alloc *sources, const ma_sub_t *coverage_cut, uint32_t qn, uint32_t tn)
{
	if (coverage_cut[qn].del || coverage_cut[tn].del) return NULL;
	ma_hit_t_alloc *x = &sources[qn];
	uint32_t i;
	for (i = 0; i < x->length; i++) {
		if (x->buffer[i].del) continue;
		if (coverage_cut[Get_qn(x->buffer[i])].del) continue;
		if (coverage_cut[Get_tn(x->buffer[i])].del) continue;
		if (Get_tn(x->buffer[i]) == tn && Get_qn(x->buffer[i]) == qn)
			return &x->buffer[i];
	}
	return NULL;
}

void delete_single_edge(ma_hit_t_alloc *sources, const ma_sub_t *coverage_cut, uint32_t qn, uint32_t tn)
{
	ma_hit_t *tmp = get_specific_overlap_with_del(sources, coverage_cut, qn, tn);
	if (tmp != NULL) tmp->del = 1;
}

void delete_all_edges(ma_hit_t_alloc *sources, ma_sub_t *coverage_cut, uint32_t qn)
{
	ma_hit_t_alloc* x = &sources[qn];
	uint32_t i;
	for (i = 0; i < x->length; i++) {
		x->buffer[i].del = 1;
		delete_single_edge(sources, coverage_cut, Get_tn(x->buffer[i]), Get_qn(x->buffer[i]));
	}
	coverage_cut[qn].del = 1;
}

void ma_hit_contained_advance(ma_hit_t_alloc *sources, long long n_read, ma_sub_t *coverage_cut, R_to_U *ruIndex, int max_hang, int min_ovlp)
{
	int32_t r;
	long long i, j, n_strong_contain = 0, n_weak_contain = 0;
	asg_arc_t t;
	ma_hit_t *h = NULL;
	ma_sub_t *sq = NULL;
	ma_sub_t *st = NULL;

	for (i = 0; i < n_read; ++i) {
		if (coverage_cut[i].del) continue;
		for (j = 0; j < (long long)sources[i].length; j++) {
			h = &sources[i].buffer[j];
			//check the corresponding two reads 
			sq = &coverage_cut[Get_qn(*h)];
			st = &coverage_cut[Get_tn(*h)];
			/****************************may have trio bugs********************************/
			if (sq->del || st->del) continue;
			if (h->del) continue;
			/****************************may have trio bugs********************************/
			r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
			//assert(r != MA_HT_INT && r != MA_HT_SHORT_OVLP);
			if (r == MA_HT_QCONT) {
				if (h->ml || 1) {
					h->del = 1;
					++n_strong_contain;
					delete_single_edge(sources, coverage_cut, Get_tn(*h), Get_qn(*h));
					delete_all_edges(sources, coverage_cut, Get_qn(*h));
					set_R_to_U(ruIndex, Get_qn(*h), Get_tn(*h), 0);

					// if (delete_all_edges_carefully(sources, coverage_cut, max_hang, min_ovlp, Get_qn(*h)) == 0)
					//     set_R_to_U(ruIndex, Get_qn(*h), Get_tn(*h), 0);
					// sq->del = 1;
					// set_R_to_U(ruIndex, Get_qn(*h), Get_tn(*h), 0);
				} else {
					sq->weak_contain = 1;
					++n_weak_contain;
				}
			} else if (r == MA_HT_TCONT) {
				if (h->ml || 1) {
					h->del = 1;
					++n_strong_contain;
					delete_single_edge(sources, coverage_cut, Get_tn(*h), Get_qn(*h));
					delete_all_edges(sources, coverage_cut, Get_tn(*h));
					set_R_to_U(ruIndex, Get_tn(*h), Get_qn(*h), 0);

					// if (delete_all_edges_carefully(sources, coverage_cut, max_hang, min_ovlp, Get_tn(*h)) == 0)
					//     set_R_to_U(ruIndex, Get_tn(*h), Get_qn(*h), 0);
					// st->del = 1;
					// set_R_to_U(ruIndex, Get_tn(*h), Get_qn(*h), 0);
				} else {
					st->weak_contain = 1;
					++n_weak_contain;
				}
			}
		}
	}

	transfor_R_to_U(ruIndex);

	for (i = 0; i < n_read; ++i) {
		int m = 0;
		for (j = 0; j < (long long)sources[i].length; j++) {
			ma_hit_t *h = &(sources[i].buffer[j]);
			if (h->del) continue;
			/// both the qn and tn have not been deleted
			if (coverage_cut[Get_qn(*h)].del != 1 && coverage_cut[Get_tn(*h)].del != 1)
				h->del = 0, ++m;
			else h->del = 1;
		}
		/// sources[i].length == 0 means all overlapped reads with read i are the contained reads
		if (m == 0) coverage_cut[i].del = 1;
	}

	fprintf(stderr, "[M::%s] %lld strong containments; %lld weak containments\n", __func__,
			n_strong_contain, n_weak_contain);
}

/************************************
 * Graph construction and reduction *
 ************************************/

static inline void asg_con_push(asg_t *g, uint32_t lower, uint32_t upper, int rev)
{
	if (g->n_con == g->m_con) {
		g->m_con = g->m_con? g->m_con<<1 : 16;
		REALLOC(g->contain, g->m_con);
	}
	g->contain[g->n_con++] = (uint64_t)lower << 32 | upper << 1 | (!!rev);
}

void asg_con_sort(asg_t *g)
{
	if (g->n_con > 1) radix_sort_ha64(g->contain, g->contain + g->n_con);
}

void asg_con_index(asg_t *g)
{
	uint32_t i, k;
	if (g->n_con == 0 || g->contain == 0) return;
	if (g->con_idx) free(g->con_idx);
	CALLOC(g->con_idx, g->n_seq);
	for (k = 0, i = 1; i < g->n_con; ++i)
		if (g->contain[k] != g->contain[i])
			g->contain[k++] = g->contain[i];
	g->n_con = k;
	for (i = 1, k = 0; i <= g->n_con; ++i)
		if (i == g->n_con || g->contain[i-1]>>32 != g->contain[i]>>32)
			g->con_idx[g->contain[i-1]>>32] = (uint64_t)k << 32 | (i - k), k = i;
}

asg_t *ma_sg_gen(const ma_hit_t_alloc* sources, long long n_read, const ma_sub_t *coverage_cut, int max_hang, int min_ovlp)
{
	size_t i, j;
	asg_t *g;
	g = asg_init();

	// add seq to graph, seq just save the length of each read
	for (i = 0; i < (uint64_t)n_read; ++i) {
		///if a read has been deleted, should we still add them?
		asg_seq_set(g, i, coverage_cut[i].e - coverage_cut[i].s, coverage_cut[i].del);
		g->seq[i].c = coverage_cut[i].c;
	}

	g->seq_vis = (uint8_t*)calloc(g->n_seq*2, sizeof(uint8_t));

	for (i = 0; i < (uint64_t)n_read; ++i) {
		for (j = 0; j < sources[i].length; ++j) {
			int r, ql, tl;
			asg_arc_t t, *p;
			const ma_hit_t *h = &sources[i].buffer[j];
			uint32_t qn, tn;
			if (h->del) continue;
			qn = Get_qn(*h);
			tn = Get_tn(*h);
			ql = coverage_cut[qn].e - coverage_cut[qn].s;
			tl = coverage_cut[tn].e - coverage_cut[tn].s;
			r = ma_hit2arc(h, ql, tl, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
			assert(r >= 0);
			if (r >= 0) {
				p = asg_arc_pushp(g);
				*p = t;
			} else if (r == MA_HT_QCONT) {
				assert(h->ml == 0);
				asg_con_push(g, h->qns>>32, h->tn, h->rev);
			} else if (r == MA_HT_TCONT) {
				assert(h->ml == 0);
				asg_con_push(g, h->tn, h->qns>>32, h->rev);
			}
		}
	}
	asg_cleanup(g);
	g->r_seq = g->n_seq;
	return g;
}

typedef struct {
	uint32_t len;
	uint8_t mark; // can only be 0, 1 or 2
} trinfo_t;

// transitive reduction; see Myers, 2005
int asg_arc_del_trans(asg_t *g, int fuzz)
{
	trinfo_t *info;
	///n_vtx = number of seq * 2; the reason is that each read has two direction (query->target, target->query)
	uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0;
	///at first, all nodes should be set to vacant
	CALLOC(info, n_vtx);

	/**v is the id+direction of a node, 
	 * the high 31-bit is the id, 
	 * and the lowest 1-bit is the direction
	 * (0 means query-to-target, 1 means target-to-query)**/
	for (v = 0; v < n_vtx; ++v) {
		///nv is the number of overlaps with v(qn+direction)
		uint32_t L, i, nv = asg_arc_n(g, v);
		///av is the array of v
		asg_arc_t *av = asg_arc_a(g, v);
		///that means in this direction, read v is not overlapped with any other reads
		if (nv == 0) continue; // no hits

		// if the read itself has been removed
		if (g->seq[v>>1].del) {
			for (i = 0; i < nv; ++i) av[i].del = 1, ++n_reduced;
			continue;
		}

    /**
	********************************query-to-target overlap****************************
	case 1: u = 0, rev = 0                           in the view of target: direction is 1 
	query: CCCCCCCCTAATTAAAAT                        target: TAATTAAAATGGGGGG (use ex-target as query)
	               ||||||||||         <--->                  ||||||||||
	       target: TAATTAAAATGGGGGG           query: CCCCCCCCTAATTAAAAT (use ex-query as target)

	case 2: u = 0, rev = 1                           in the view of target: direction is 0
	query: CCCCCCCCTAATTAAAAT					     target: CCCCCCATTTTAATTA  (use ex-target as query)
                   ||||||||||        <--->                         ||||||||||
	       target: TAATTAAAATGGGGGG                         query: ATTTTAATTAGGGGGGGG  (use ex-query as target)
	********************************query-to-target overlap****************************

	********************************target-to-query overlap****************************
	case 3: u = 1, rev = 0                           in the view of target: direction is 0
			 query: AAATAATATCCCCCCGCG                target: GGGCCGGCAAATAATAT (use ex-target as query)
					|||||||||          <--->                          |||||||||
	target: GGGCCGGCAAATAATAT                                  query: AAATAATATCCCCCCGCG (use ex-query as target)

	case 4: u = 1, rev = 1                          in the view of target: direction is 1
	         query: AAATAATATCCCCCCGCG                        target: ATATTATTTGCCGGCCC (use ex-target as query)
                    |||||||||             <--->                       |||||||||
	target: GGGCCGGCAAATAATAT                          query: CGCGGGGGATATTATTT (use ex-query as target)
	********************************target-to-query overlap****************************

	p->ul: |____________31__________|__________1___________|______________32_____________|
	                    qns            direction of overlap       length of this node (not overlap length)
						                (in the view of query)
	p->v : |___________31___________|__________1___________|
				        tns             reverse direction of overlap
						              (in the view of target)
	p->ol: overlap length
    **/

		// all outnode of v should be set to "not reduce"
		for (i = 0; i < nv; ++i) {
			uint32_t w = av[i].v;
			info[w].mark = g->seq[w>>1].del? 2 : 1;
			info[w].len = asg_arc_len(av[i]);
		}

		// length of node (not overlap length)
		// av[nv-1] is longest out-dege
		/**
		 * v---------------
		 *   w1---------------
		 *      w2--------------
		 *         w3--------------
		 *            w4--------------
		 *               w5-------------
		 * for v, the longest out-edge is v->w5
		 **/
		L = asg_arc_len(av[nv-1]) + fuzz;

		for (i = 0; i < nv; ++i) {
			uint32_t w = av[i].v;
			uint32_t j, nw = asg_arc_n(g, w);
			asg_arc_t *aw = asg_arc_a(g, w);
			if (info[w].mark != 1) continue;
			for (j = 0; j < nw; ++j) {
				uint32_t x, sum = asg_arc_len(aw[j]) + asg_arc_len(av[i]);
				if (sum > L) break;
				x = aw[j].v;
				if (info[x].mark == 1 && sum < info[x].len + fuzz && sum + fuzz > info[x].len)
					info[x].mark = 2;
			}
		}
		#if 0
		for (i = 0; i < nv; ++i) {
			uint32_t w = av[i].v;
			uint32_t j, nw = asg_arc_n(g, w);
			asg_arc_t *aw = asg_arc_a(g, w);
			for (j = 0; j < nw && (j == 0 || asg_arc_len(aw[j]) < fuzz); ++j)
				if (info[aw[j].v].mark) info[aw[j].v].mark = 2;
		}
		#endif
		// remove edges
		for (i = 0; i < nv; ++i) {
			if (info[av[i].v].mark == 2) av[i].del = 1, ++n_reduced;
			info[av[i].v].mark = 0;
		}
	}
	free(info);

	if (n_reduced) {
		asg_cleanup(g);
		asg_symm(g);
	}
	return n_reduced;
}
