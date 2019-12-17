#include <stdio.h>
#include <stdlib.h>
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include "Overlaps.h"
#include "ksort.h"
#include "Process_Read.h"
#include "CommandLines.h"

KDQ_INIT(uint64_t)

#define ma_hit_key_tn(a) ((a).tn)
KRADIX_SORT_INIT(hit_tn, ma_hit_t, ma_hit_key_tn, member_size(ma_hit_t, tn))

#define ma_hit_key_qns(a) ((a).qns)
KRADIX_SORT_INIT(hit_qns, ma_hit_t, ma_hit_key_qns, member_size(ma_hit_t, qns))

#define asg_arc_key(a) ((a).ul)
KRADIX_SORT_INIT(asg, asg_arc_t, asg_arc_key, 8)

#define generic_key(x) (x)
KRADIX_SORT_INIT(arch64, uint64_t, generic_key, 8)

KSORT_INIT_GENERIC(uint32_t)

///actually min_thres = MAX_SHORT_TIPS + 1 there are MAX_SHORT_TIPS reads
long long min_thres = MAX_SHORT_TIPS + 1;

void ma_hit_sort_tn(ma_hit_t *a, long long n)
{
	radix_sort_hit_tn(a, a + n);
}

void ma_hit_sort_qns(ma_hit_t *a, long long n)
{
	radix_sort_hit_qns(a, a + n);
}




asg_t *asg_init(void)
{
	return (asg_t*)calloc(1, sizeof(asg_t));
}

void asg_destroy(asg_t *g)
{
	if (g == 0) return;
	free(g->seq); free(g->idx); free(g->arc); free(g->seq_vis); free(g);
}

void asg_arc_sort(asg_t *g)
{
	radix_sort_asg(g->arc, g->arc + g->n_arc);
}

void ma_ug_destroy(ma_ug_t *ug)
{
	uint32_t i;
	if (ug == 0) return;
	for (i = 0; i < ug->u.n; ++i) {
		free(ug->u.a[i].a);
		free(ug->u.a[i].s);
	}
	free(ug->u.a);
	asg_destroy(ug->g);
	free(ug);
}

uint64_t *asg_arc_index_core(size_t max_seq, size_t n, const asg_arc_t *a)
{
	size_t i, last;
	uint64_t *idx;
	idx = (uint64_t*)calloc(max_seq * 2, 8);


    /**
		 * ul: |____________31__________|__________1___________|______________32_____________|
	                       qns            direction of overlap       length of this node (not overlap length)
	**/
    ///so if we use high 32-bit, we store the index of each qn with two direction
	for (i = 1, last = 0; i <= n; ++i)
		if (i == n || a[i-1].ul>>32 != a[i].ul>>32)
			idx[a[i-1].ul>>32] = (uint64_t)last<<32 | (i - last), last = i;

    
	return idx;
}

void asg_arc_index(asg_t *g)
{
	if (g->idx) free(g->idx);
	g->idx = asg_arc_index_core(g->n_seq, g->n_arc, g->arc);
}

void asg_seq_set(asg_t *g, int sid, int len, int del)
{
	///just malloc size
	if (sid >= g->m_seq) {
		g->m_seq = sid + 1;
		kv_roundup32(g->m_seq);
		g->seq = (asg_seq_t*)realloc(g->seq, g->m_seq * sizeof(asg_seq_t));
	}


	if (sid >= g->n_seq) g->n_seq = sid + 1;
	
    g->seq[sid].del = !!del;
    g->seq[sid].len = len;
    // if(g->seq[sid].del)
    // {
    //     g->seq[sid].len = 0;
    // }
}


// hard remove arcs marked as "del"
void asg_arc_rm(asg_t *g)
{
	/**
	p->ul: |____________31__________|__________1___________|______________32_____________|
	                    qns            direction of overlap       length of this node (not overlap length)
	p->v : |___________31___________|__________1___________|
				        tns              relative strand between query and target
	p->ol: overlap length
	**/
	uint32_t e, n;
	///just clean arc requiring: 1. arc it self must be available 2. both the query and target are available
	for (e = n = 0; e < g->n_arc; ++e) {
		//u and v is the read id
		uint32_t u = g->arc[e].ul>>32, v = g->arc[e].v;
		if (!g->arc[e].del && !g->seq[u>>1].del && !g->seq[v>>1].del)
			g->arc[n++] = g->arc[e];
	}
	if (n < g->n_arc) { // arc index is out of sync
		if (g->idx) free(g->idx);
		g->idx = 0;
	}
	g->n_arc = n;
}

void asg_cleanup(asg_t *g)
{
    ///remove overlaps, instead of reads
	asg_arc_rm(g);
	if (!g->is_srt) {
		/**
		 * sort by ul, that is, sort by qns + direction
		 * ul: |____________31__________|__________1___________|______________32_____________|
	                       qns            direction of overlap       length of this node (not overlap length)
		**/
		asg_arc_sort(g);
		g->is_srt = 1;
	}
	///index the overlaps in graph with query id
	if (g->idx == 0) asg_arc_index(g);
}


// delete multi-arcs
/**
 * remove edges like:   v has two out-edges to w
**/
int asg_arc_del_multi(asg_t *g)
{
	//the number of nodes are number of read times 2
	uint32_t *cnt, n_vtx = g->n_seq * 2, n_multi = 0, v;
	cnt = (uint32_t*)calloc(n_vtx, 4);
	for (v = 0; v < n_vtx; ++v) {
		///out-nodes of v
		asg_arc_t *av = asg_arc_a(g, v);
		int32_t i, nv = asg_arc_n(g, v);
		///if v just have one out-node, there is no muti-edge
		if (nv < 2) continue;
		for (i = nv - 1; i >= 0; --i) ++cnt[av[i].v];
		for (i = nv - 1; i >= 0; --i)
			if (--cnt[av[i].v] != 0)
				av[i].del = 1, ++n_multi;
	}
	free(cnt);
	if (n_multi) asg_cleanup(g);
	fprintf(stderr, "[M::%s] removed %d multi-arcs\n", __func__, n_multi);
	return n_multi;
}

// remove asymmetric arcs: u->v is present, but v'->u' not
int asg_arc_del_asymm(asg_t *g)
{
	uint32_t e, n_asymm = 0;
	///g->n_arc is the number of overlaps 
	for (e = 0; e < g->n_arc; ++e) {
		uint32_t v = g->arc[e].v^1, u = g->arc[e].ul>>32^1;
		uint32_t i, nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		for (i = 0; i < nv; ++i)
			if (av[i].v == u) break;
		if (i == nv) g->arc[e].del = 1, ++n_asymm;
	}
	if (n_asymm) asg_cleanup(g);
	fprintf(stderr, "[M::%s] removed %d asymmetric arcs\n", __func__, n_asymm);
	return n_asymm;
}


void asg_symm(asg_t *g)
{
	asg_arc_del_multi(g);
	asg_arc_del_asymm(g);
	g->is_symm = 1;
}

void init_ma_hit_t_alloc(ma_hit_t_alloc* x)
{
    x->size = 0;
    x->buffer = NULL;
    x->length = 0;
}

void clear_ma_hit_t_alloc(ma_hit_t_alloc* x)
{
    x->length = 0;
}

void resize_ma_hit_t_alloc(ma_hit_t_alloc* x, uint64_t size)
{
    if(size > x->size)
    {
        x->size = size;
        x->buffer = (ma_hit_t*)realloc(x->buffer, x->size*sizeof(ma_hit_t));
    }
}

void destory_ma_hit_t_alloc(ma_hit_t_alloc* x)
{
    free(x->buffer);
}


void add_ma_hit_t_alloc(ma_hit_t_alloc* x, ma_hit_t* element)
{
    if(x->length + 1 > x->size)
    {
        x->size = (x->length + 1) * 2;
        x->buffer = (ma_hit_t*)realloc(x->buffer, x->size*sizeof(ma_hit_t));
    }

    x->buffer[x->length] = (*element);
    x->length++;
}



void init_Assembly_Graph(Assembly_Graph* x)
{
    init_ma_hit_t_alloc(&(x->overlaps));
}

void destory_Assembly_Graph(Assembly_Graph* x)
{
    destory_ma_hit_t_alloc(&(x->overlaps));
}

long long get_specific_overlap(ma_hit_t_alloc* x, uint32_t qn, uint32_t tn)
{
    long long i;
    for (i = 0; i < x->length; i++)
    {
        if(x->buffer[i].tn == tn 
        &&
          ((uint32_t)(x->buffer[i].qns>>32)) == qn)
        {
            return i;
        }        
    }

    return -1;
}



void collect_ma_hit_t(ma_hit_t_alloc* dest, ma_hit_t_alloc* sources, long long num_sources)
{
    long long bi_overlaps = 0;
    long long si_overlaps = 0;
    long long i, j, index;
    uint32_t qn, tn;
    ma_hit_t new_element;
    for (i = 0; i < num_sources; i++)
    {
        resize_ma_hit_t_alloc(dest, dest->length + sources[i].length);
        for (j = 0; j < sources[i].length; j++)
        {
            qn = sources[i].buffer[j].qns>>32;
            tn = sources[i].buffer[j].tn;

            index = get_specific_overlap(&(sources[tn]), tn, qn);


            if(index != -1)
            {
                
                // fprintf(stderr, "\n+qn: %d, tn: %d, qs: %d, qe: %d, ts: %d, te: %d\n", qn, tn,
                // (uint32_t)(sources[i].buffer[j].qns), sources[i].buffer[j].qe,
                // sources[i].buffer[j].ts, sources[i].buffer[j].te);

                // fprintf(stderr, "-qn: %d, tn: %d, qs: %d, qe: %d, ts: %d, te: %d\n\n", 
                // sources[tn].buffer[index].qns>>32, sources[tn].buffer[index].tn,
                // (uint32_t)(sources[tn].buffer[index].qns), sources[tn].buffer[index].qe,
                // sources[tn].buffer[index].ts, sources[tn].buffer[index].te);
                
                if(qn <= tn)
                {

                }
                

                bi_overlaps++;
            }
            else
            {
                si_overlaps++;
            }
        }
    }


    fprintf(stderr, "bi_overlaps: %d, si_overlaps: %d\n", bi_overlaps, si_overlaps);
    
}


inline void set_reverse_overlap(ma_hit_t* dest, ma_hit_t* source)
{
    dest->qns = Get_tn(*source);
    dest->qns = dest->qns << 32;
    dest->qns = dest->qns | Get_ts(*source);
    dest->qe = Get_te(*source);


    dest->tn = Get_qn(*source);
    dest->ts = Get_qs(*source);
    dest->te = Get_qe(*source);

    dest->rev = source->rev;
    dest->el = source->el;


    /****************************may have bugs********************************/
    if(dest->ml == 0 || source->ml == 0)
    {
        dest->ml = source->ml = 0;
    }
    else
    {
        dest->ml = source->ml = 1;
    }


    if(dest->no_l_indel == 0 || source->no_l_indel == 0)
    {
        dest->no_l_indel = source->no_l_indel = 0;
    }
    else
    {
        dest->no_l_indel = source->no_l_indel = 1;
    }
    
    /****************************may have bugs********************************/
    dest->bl = Get_qe(*dest) - Get_qs(*dest);
}

void normalize_ma_hit_t(ma_hit_t_alloc* sources, long long num_sources)
{
    long long bi_overlaps = 0;
    long long si_overlaps = 0;
    long long i, j, index;
    uint32_t qn, tn;
    ma_hit_t new_element;
    long long qLen_0, qLen_1;
    for (i = 0; i < num_sources; i++)
    {
        for (j = 0; j < sources[i].length; j++)
        {
            qn = Get_qn(sources[i].buffer[j]);
            tn = Get_tn(sources[i].buffer[j]);
            
            sources[i].buffer[j].bl = Get_qe(sources[i].buffer[j]) - Get_qs(sources[i].buffer[j]);

            index = get_specific_overlap(&(sources[tn]), tn, qn);


            if(index != -1)
            {
                qLen_0 = Get_qe(sources[i].buffer[j]) - Get_qs(sources[i].buffer[j]);
                qLen_1 = Get_qe(sources[tn].buffer[index]) - Get_qs(sources[tn].buffer[index]);

                if(qLen_0 == qLen_1)
                {
                    ///qn must be not equal to tn
                    ///make sources[qn] = sources[tn] if qn > tn
                    if(qn < tn)
                    {   
                        set_reverse_overlap(&(sources[tn].buffer[index]), &(sources[i].buffer[j]));
                    }
                }
                else if(qLen_0 > qLen_1)
                {
                    
                    set_reverse_overlap(&(sources[tn].buffer[index]), &(sources[i].buffer[j]));
                }
                

                bi_overlaps++;
            }
            else
            {
                ///must have this line
                new_element.ml = 1;
                new_element.no_l_indel = 1;
                set_reverse_overlap(&new_element, &(sources[i].buffer[j]));
                add_ma_hit_t_alloc(&(sources[tn]), &new_element);
                si_overlaps++;
            }
        }
    }


    /**
    for (i = 0; i < num_sources; i++)
    {
        ma_hit_sort_qns(sources[i].buffer, sources[i].length);
    }
    **/    
}


void normalize_ma_hit_t_single_side(ma_hit_t_alloc* sources, long long num_sources)
{
    double startTime = Get_T();

    long long bi_overlaps = 0;
    long long si_overlaps = 0;
    long long i, j, index;
    uint32_t qn, tn;
    ma_hit_t new_element;
    long long qLen_0, qLen_1, m;
    for (i = 0; i < num_sources; i++)
    {
        m = 0;
        for (j = 0; j < sources[i].length; j++)
        {
            qn = Get_qn(sources[i].buffer[j]);
            tn = Get_tn(sources[i].buffer[j]);

            sources[i].buffer[j].bl = Get_qe(sources[i].buffer[j]) - Get_qs(sources[i].buffer[j]);

            index = get_specific_overlap(&(sources[tn]), tn, qn);


            if(index != -1)
            {
                qLen_0 = Get_qe(sources[i].buffer[j]) - Get_qs(sources[i].buffer[j]);
                qLen_1 = Get_qe(sources[tn].buffer[index]) - Get_qs(sources[tn].buffer[index]);

                if(qLen_0 == qLen_1)
                {
                    ///qn must be not equal to tn
                    ///make sources[qn] = sources[tn] if qn > tn
                    if(qn < tn)
                    {
                        set_reverse_overlap(&(sources[tn].buffer[index]), &(sources[i].buffer[j]));
                    }
                }
                else if(qLen_0 > qLen_1)
                {
                    set_reverse_overlap(&(sources[tn].buffer[index]), &(sources[i].buffer[j]));
                }
                
                sources[i].buffer[m] = sources[i].buffer[j];
                m++;
                bi_overlaps++;
            }
        }

        sources[i].length = m;
    }


    fprintf(stderr, "[M::%s] takes %0.2fs\n\n", __func__, Get_T()-startTime);
}




void ma_hit_contained(ma_hit_t_alloc* sources, long long n_read, ma_sub_t *coverage_cut, int max_hang, int min_ovlp)
{
    double startTime = Get_T();
	int32_t r;
	size_t i, j, m;
	asg_arc_t t;
    for (i = 0; i < n_read; ++i) 
    {
        for (j = 0; j < sources[i].length; j++)
        {
            ma_hit_t *h = &(sources[i].buffer[j]);
            //check the corresponding two reads 
		    ma_sub_t *sq = &(coverage_cut[Get_qn(*h)]);
            ma_sub_t *st = &(coverage_cut[Get_tn(*h)]);
            r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, MAX_HANG_PRE, min_ovlp, &t);
            ///r could not be MA_HT_SHORT_OVLP or MA_HT_INT
            if (r == MA_HT_QCONT) 
            {
                sq->del = 1;
            }
		    else if (r == MA_HT_TCONT) 
            {
                st->del = 1;
            }
        }
    }



    for (i = 0; i < n_read; ++i) 
    {
        m = 0;
        for (j = 0; j < sources[i].length; j++)
        {
            ma_hit_t *h = &(sources[i].buffer[j]);
            ///both the qn and tn have not been deleted
            if(coverage_cut[Get_qn(*h)].del != 1 && coverage_cut[Get_tn(*h)].del != 1)
            {
                sources[i].buffer[m] = *h;
                m++;
            }
        }
        sources[i].length = m;
        ///may have bugs here 
        ///if sources[i].length == 0, that means all overlapped reads with read i are the contained reads
        if(sources[i].length == 0)
        {
            coverage_cut[i].del = 1;
        }
    }
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

}



void ma_hit_flt(ma_hit_t_alloc* sources, long long n_read, ma_sub_t *coverage_cut, int max_hang, int min_ovlp)
{
    double startTime = Get_T();
	size_t i, j, m;
	asg_arc_t t;
	uint64_t tot_dp = 0, tot_len = 0;

    for (i = 0; i < n_read; ++i) 
    {
        m = 0;
        for (j = 0; j < sources[i].length; j++)
        {
            ma_hit_t *h = &(sources[i].buffer[j]);
            //check the corresponding two reads 
		    const ma_sub_t *sq = &(coverage_cut[Get_qn(*h)]);
            const ma_sub_t *st = &(coverage_cut[Get_tn(*h)]);
            int r;
		    if (sq->del || st->del) continue;
            ///[sq->s, sq->e) and [st->s, st->e) are the high coverage region in query and target
            ///here just exculde the overhang?
            ///in miniasm the 5-th option is 0.5, instead of 0.8
            /**note!!! h->qn and h->qs have been normalized by sq->s
             * h->ts and h->tn have been normalized by sq->e
             **/  
            ///here the max_hang = 1000, MAX_HANG_PRE = 0.8, min_ovlp = 500
            ///for me, there should not have any overhang..so r cannot be equal to MA_HT_INT
            r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, MAX_HANG_PRE, min_ovlp, &t);


            ///for me, there should not have any overhang..so r cannot be equal to MA_HT_INT
            ///and I think if we use same min_ovlp in all functions, r also cannot be MA_HT_SHORT_OVLP
            ///so it does not matter we have ma_hit2arc or not
		    if (r >= 0 || r == MA_HT_QCONT || r == MA_HT_TCONT)
            {
                sources[i].buffer[m] = *h;
                m++;
            }/**
            else
            {
                fprintf(stderr, "shit\n");
            }
            **/
            
        }
        sources[i].length = m;
        if(sources[i].length == 0)
        {
            (coverage_cut)[i].del = 1;
        }
    }
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
}



///a is the overlap vector, n is the length of overlap vector
///min_dp is used for coverage droping
///select reads with coverage >= min_dp
void ma_hit_sub_back(int min_dp, ma_hit_t_alloc* sources, long long n_read, uint64_t* readLen, 
long long mini_overlap_length, ma_sub_t** coverage_cut)
{
    (*coverage_cut) = (ma_sub_t*)malloc(sizeof(ma_sub_t)*n_read);

	size_t i, j, n_remained = 0;
	kvec_t(uint32_t) b = {0,0,0};

	///all overlaps in vector a has been sorted by qns
	///so for overlaps of one reads, it must be contiguous
	for (i = 0; i < n_read; ++i) 
    {
        // fprintf(stderr, "i: %d, n_read: %d\n", i, n_read);

        kv_resize(uint32_t, b, sources[i].length);
        b.n = 0;
        for (j = 0; j < sources[i].length; j++)
        {
            uint32_t qs, qe;
            qs = Get_qs(sources[i].buffer[j]);
            qe = Get_qe(sources[i].buffer[j]);
            kv_push(uint32_t, b, qs<<1);
			kv_push(uint32_t, b, qe<<1|1);
        }

        ///we can identify the qs and qe by the 0-th bit
		ks_introsort_uint32_t(b.n, b.a);
        ma_sub_t max, max2;
        max.s = max.e = max.del = max2.s = max2.e = max2.del = 0;
        int dp, start;
        ///max is the longest subregion, max2 is the second longest subregion
        for (j = 0, dp = 0; j < b.n; ++j) 
        {
            int old_dp = dp;
            ///if a[j] is qe
            if (b.a[j]&1) 
            {
                --dp;
            }
            else
            {
                ++dp;
            } 
            
            ///min_dp is the coverage drop threshold
            ///there are two cases: 1. old_dp = dp + 1 (b.a[j] is qe); 2. old_dp = dp - 1 (b.a[j] is qs);
            ///if one read has multiple separate sub-regions with coverage > 3, does miniasm only select the longest one?
            if (old_dp < min_dp && dp >= min_dp) ///old_dp < dp, b.a[j] is qs
            { ///case 2, a[j] is qs
                start = b.a[j]>>1;
            } 
            else if (old_dp >= min_dp && dp < min_dp) ///old_dp > min_dp, b.a[j] is qe
            {
                int len = (b.a[j]>>1) - start;
                if (len > max.e - max.s) 
                {
                    max2 = max; 
                    max.s = start;
                    max.e = b.a[j]>>1;
                }
                else if (len > max2.e - max2.s) 
                {
                    max2.s = start; 
                    max2.e = b.a[j]>>1;
                }
            }
        }



        if (max.e - max.s > 0) 
        {
            (*coverage_cut)[i].s = max.s;
            (*coverage_cut)[i].e = max.e;
            (*coverage_cut)[i].del = 0;
            ++n_remained;
		} 
        else 
        {
            (*coverage_cut)[i].del = 1;
        }
	}




    ma_hit_t* p;
    ma_sub_t* rq;
    ma_sub_t* rt;
    long long m = 0;
    for (i = 0; i < n_read; ++i) 
    {
        m = 0;
        for (j = 0; j < sources[i].length; j++)
        {
            ///this is a overlap
            p = &(sources[i].buffer[j]);

            rq = &((*coverage_cut)[Get_qn(*p)]);
            rt = &((*coverage_cut)[Get_tn(*p)]);
            ///if any of target read and the query read has no enough coverage
		    if (rq->del || rt->del) continue;
            int qs, qe, ts, te;




            ///target and query in different strand
            if (p->rev) 
            {
                qs = p->te < rt->e? Get_qs(*p): Get_qs(*p) + (p->te - rt->e);
                qe = p->ts > rt->s? p->qe : p->qe - (rt->s - p->ts);
                ts = p->qe < rq->e? p->ts : p->ts + (p->qe - rq->e);
                te = Get_qs(*p) > rq->s? p->te : p->te - (rq->s - Get_qs(*p));
            }
            else ///target and query in same strand
            { 
                ///note: ts is the targe start in this overlap, 
                ///while rt->s is the high coverage start in the whole target (not only in this overlap)
                ///so this line is to normalize the qs in quey to high coverage region
                qs = p->ts > rt->s? Get_qs(*p): Get_qs(*p) + (rt->s - p->ts); //(rt->s - p->ts) is the offset        
                qe = p->te < rt->e? p->qe : p->qe - (p->te - rt->e);//(p->te - rt->e) is the offset
                ts = Get_qs(*p) > rq->s? p->ts : p->ts + (rq->s - Get_qs(*p));//(rq->s - Get_qs(*p) is the offset
                te = p->qe < rq->e? p->te : p->te - (p->qe - rq->e);//(p->qe - rq->e) is the offset
		    }

            

            //cut by self coverage
            qs = (qs > rq->s? qs : rq->s) - rq->s;
            qe = (qe < rq->e? qe : rq->e) - rq->s;
            ts = (ts > rt->s? ts : rt->s) - rt->s;
            te = (te < rt->e? te : rt->e) - rt->s;

            if (qe - qs >= mini_overlap_length && te - ts >= mini_overlap_length) 
            {
                ///p->qns = p->qns>>32<<32 | qs;
                p->qns = p->qns>>32;
                p->qns = p->qns << 32;
                p->qns = p->qns | qs;

                p->qe = qe;
                p->ts = ts;
                p->te = te;
                sources[i].buffer[m] = *p;
                m++;
		    }
        }
        sources[i].length = m;
    }







    
    


	free(b.a);
    ///free((*coverage_cut));
}


///a is the overlap vector, n is the length of overlap vector
///min_dp is used for coverage droping
///select reads with coverage >= min_dp
void ma_hit_sub(int min_dp, ma_hit_t_alloc* sources, long long n_read, uint64_t* readLen, 
long long mini_overlap_length, ma_sub_t** coverage_cut)
{
    double startTime = Get_T();

    (*coverage_cut) = (ma_sub_t*)malloc(sizeof(ma_sub_t)*n_read);

	size_t i, j, n_remained = 0;
	kvec_t(uint32_t) b = {0,0,0};

	///all overlaps in vector a has been sorted by qns
	///so for overlaps of one reads, it must be contiguous
	for (i = 0; i < n_read; ++i) 
    {
        if(min_dp <= 1)
        {
            (*coverage_cut)[i].s = 0;
            (*coverage_cut)[i].e = readLen[i];
            (*coverage_cut)[i].del = 0;
            ++n_remained;
            continue;
        }


        kv_resize(uint32_t, b, sources[i].length);
        b.n = 0;
        for (j = 0; j < sources[i].length; j++)
        {
            uint32_t qs, qe;
            qs = Get_qs(sources[i].buffer[j]);
            qe = Get_qe(sources[i].buffer[j]);
            kv_push(uint32_t, b, qs<<1);
			kv_push(uint32_t, b, qe<<1|1);
        }

        ///we can identify the qs and qe by the 0-th bit
		ks_introsort_uint32_t(b.n, b.a);
        ma_sub_t max, max2;
        max.s = max.e = max.del = max2.s = max2.e = max2.del = 0;
        int dp, start;
        ///max is the longest subregion, max2 is the second longest subregion
        for (j = 0, dp = 0; j < b.n; ++j) 
        {
            int old_dp = dp;
            ///if a[j] is qe
            if (b.a[j]&1) 
            {
                --dp;
            }
            else
            {
                ++dp;
            } 
            
            /**
            min_dp is the coverage drop threshold
            there are two cases: 
                1. old_dp = dp + 1 (b.a[j] is qe); 2. old_dp = dp - 1 (b.a[j] is qs);
            if one read has multiple separate sub-regions with coverage >= min_dp, 
            does miniasm only select the longest one?
            **/
            if (old_dp < min_dp && dp >= min_dp) ///old_dp < dp, b.a[j] is qs
            { 
                ///case 2, a[j] is qs
                start = b.a[j]>>1;
            } 
            else if (old_dp >= min_dp && dp < min_dp) ///old_dp > min_dp, b.a[j] is qe
            {
                int len = (b.a[j]>>1) - start;
                if (len > max.e - max.s) 
                {
                    max2 = max; 
                    max.s = start;
                    max.e = b.a[j]>>1;
                }
                else if (len > max2.e - max2.s) 
                {
                    max2.s = start; 
                    max2.e = b.a[j]>>1;
                }
            }
        }


        ///max.e - max.s is the 
        if (max.e - max.s > 0) 
        {
            (*coverage_cut)[i].s = max.s;
            (*coverage_cut)[i].e = max.e;
            (*coverage_cut)[i].del = 0;
            ++n_remained;
		} 
        else 
        {
            (*coverage_cut)[i].s = (*coverage_cut)[i].e = 0;

            (*coverage_cut)[i].del = 1;
        }
	}

	free(b.a);

    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
}




void ma_hit_chimeric(int min_dp, ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, 
long long n_read, uint64_t* readLen, ma_sub_t* coverage_cut)
{
    double startTime = Get_T();

	int i, j, k, n_remove = 0;
	kvec_t(uint32_t) b = {0,0,0};

	for (i = 0; i < n_read; ++i) 
    {

        kv_resize(uint32_t, b, readLen[i]);
        memset(b.a, 0, sizeof(uint32_t)*readLen[i]);

        
        
        for (j = 0; j < sources[i].length; j++)
        {
            uint32_t qs, qe;
            qs = Get_qs(sources[i].buffer[j]);
            qe = Get_qe(sources[i].buffer[j]);


            ///if(qe - qs < 1000) continue;
            for (k = qs; k < qe; k++)
            {
                b.a[k]++;
            }
        }

    
        // for (j = 0; j < reverse_sources[i].length; j++)
        // {
        //     uint32_t qs, qe;
        //     qs = Get_qs(reverse_sources[i].buffer[j]);
        //     qe = Get_qe(reverse_sources[i].buffer[j]);

        //     for (k = qs; k < qe; k++)
        //     {
        //         b.a[k]++;
        //     }
        // }
    
        
        int left, right;
        left = -1;
        right = readLen[i];
        for (k = 0; k < readLen[i]; k++)
        {
            if(b.a[k] < min_dp)
            {
                left = k;
                break;
            }
        }
        
        for (k = readLen[i] - 1; k >= 0; k--)
        {
            if(b.a[k] < min_dp)
            {
                right = k;
                break;
            }
        }

        // if(i == 4616942 || i == 4024299 || i == 6135193)
        // {
        //     fprintf(stderr, "i: %d, left: %d, right: %d, rLen: %d, min_dp: %d\n", 
        //     i, left, right, readLen[i], min_dp);
        //     for (k = 0; k < readLen[i]; k++)
        //     {
        //         fprintf(stderr, "a[%d]: %d\n", k, b.a[k]);
        //     }
        // }
    
        if( (left <= right) && (left > 0) && (right < readLen[i] - 1))
        {
            coverage_cut[i].c = 1;
            n_remove++;
            /****************************may have bugs********************************/
            coverage_cut[i].del = 1;
            sources[i].length = 0;
            /****************************may have bugs********************************/
        }
        else
        {
            coverage_cut[i].c = 0;
        }
	}

	free(b.a);

    fprintf(stderr, "[M::%s] takes %0.2f s, n_remove: %d\n\n", __func__, Get_T()-startTime, n_remove);
}






void ma_hit_cut(int min_dp, ma_hit_t_alloc* sources, long long n_read, uint64_t* readLen, 
long long mini_overlap_length, ma_sub_t** coverage_cut)
{
    double startTime = Get_T();
	size_t i, j;
    ma_hit_t* p;
    ma_sub_t* rq;
    ma_sub_t* rt;
    long long m = 0;
    for (i = 0; i < n_read; ++i) 
    {
        m = 0;
        for (j = 0; j < sources[i].length; j++)
        {
            ///this is a overlap
            p = &(sources[i].buffer[j]);

            rq = &((*coverage_cut)[Get_qn(*p)]);
            rt = &((*coverage_cut)[Get_tn(*p)]);
            ///if any of target read and the query read has no enough coverage
		    if (rq->del || rt->del) continue;
            int qs, qe, ts, te;




            ///target and query in different strand
            if (p->rev) 
            {
                /**
                here is an example in different strand:
                                      
                                       (te) (rt->e)           (rt->s)    (ts)
                                        |      |                 |        |
                target ----------------------------------------------------
                                        ------------------------------------------- query
                                        qs                                qe
                **/
                qs = p->te < rt->e? Get_qs(*p): Get_qs(*p) + (p->te - rt->e);
                qe = p->ts > rt->s? p->qe : p->qe - (rt->s - p->ts);
                ts = p->qe < rq->e? p->ts : p->ts + (p->qe - rq->e);
                te = Get_qs(*p) > rq->s? p->te : p->te - (rq->s - Get_qs(*p));
            }
            else ///target and query in same strand
            { 
                /**
                note: ts is the targe start in this overlap, 
                while rt->s is the high coverage start in the whole target (not only in this overlap)
                so this line is to normalize the qs in quey to high coverage region
                **/
                //(rt->s - p->ts) is the offset        
                qs = p->ts > rt->s? Get_qs(*p): Get_qs(*p) + (rt->s - p->ts); 
                //(p->te - rt->e) is the offset
                qe = p->te < rt->e? p->qe : p->qe - (p->te - rt->e);
                //(rq->s - Get_qs(*p) is the offset
                ts = Get_qs(*p) > rq->s? p->ts : p->ts + (rq->s - Get_qs(*p));
                //(p->qe - rq->e) is the offset
                te = p->qe < rq->e? p->te : p->te - (p->qe - rq->e);
		    }

            

            //cut by self coverage
            //and normalize the qs, qe, ts, te by rq->s and rt->e
            qs = (qs > rq->s? qs : rq->s) - rq->s;
            qe = (qe < rq->e? qe : rq->e) - rq->s;
            ts = (ts > rt->s? ts : rt->s) - rt->s;
            te = (te < rt->e? te : rt->e) - rt->s;

            if (qe - qs >= mini_overlap_length && te - ts >= mini_overlap_length) 
            {
                ///p->qns = p->qns>>32<<32 | qs;
                p->qns = p->qns>>32;
                p->qns = p->qns << 32;
                p->qns = p->qns | qs;

                p->qe = qe;
                p->ts = ts;
                p->te = te;
                sources[i].buffer[m] = *p;
                ///fprintf(stderr, "p->del: %d\n", p->del);
                m++;
		    }
        }
        sources[i].length = m;
        if(sources[i].length == 0)
        {
            (*coverage_cut)[i].del = 1;
        }
    }

    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
}



void debug_normalize_ma_hit_t(ma_hit_t_alloc* sources, long long num_sources)
{
    long long i, j, index;
    uint32_t qn, tn;
    long long total_overlaps = 0;
    long long total_reads = 0;

    for (i = 0; i < num_sources; i++)
    {
        for (j = 0; j < sources[i].length; j++)
        {
            if(Get_qn(sources[i].buffer[j]) != i)
            {
                fprintf(stderr, "1 error 2\n");
            }
            total_overlaps += sources[i].length;
            
        }

        if(sources[i].length != 0)
        {
            total_reads++;
        }
    }

    for (i = 0; i < num_sources; i++)
    {
        for (j = 0; j < sources[i].length; j++)
        {
            qn = sources[i].buffer[j].qns>>32;
            tn = sources[i].buffer[j].tn;

            long long k;
            for (k = 0; k < sources[i].length; k++)
            {
                //here can be improved, since ma_hit_t_alloc has been sorted by tn
                if(sources[i].buffer[k].tn == tn 
                &&
                ((uint32_t)(sources[i].buffer[k].qns>>32)) == qn)
                {
                    if(k != j)
                    {
                        fprintf(stderr, "2 ERROR\n");
                    }
                }
            }
        }
    }










    for (i = 0; i < num_sources; i++)
    {
        for (j = 0; j < sources[i].length; j++)
        {
            

            
            if(sources[i].buffer[j].bl != Get_qe(sources[i].buffer[j]) - Get_qs(sources[i].buffer[j]))
            {
                fprintf(stderr, "3 error2, bl: %d, qs: %d, qe: %d\n",
                sources[i].buffer[j].bl, Get_qs(sources[i].buffer[j]),
                Get_qe(sources[i].buffer[j]));
            }
        

            qn = sources[i].buffer[j].qns>>32;
            tn = sources[i].buffer[j].tn;

            index = get_specific_overlap(&(sources[tn]), tn, qn);


            if(index == -1)
            {
                fprintf(stderr, "4 error3\n");
            }
            else
            {

                if(sources[i].buffer[j].rev != sources[tn].buffer[index].rev)
                {
                    fprintf(stderr, "5 hahaha\n");
                }

                if(sources[i].buffer[j].el != sources[tn].buffer[index].el)
                {
                    fprintf(stderr, "el hahaha\n");
                }

               if(Get_qn(sources[i].buffer[j]) != Get_tn(sources[tn].buffer[index]))
               {
                   fprintf(stderr, "6 error4\n");
               }

               if(Get_tn(sources[i].buffer[j]) != Get_qn(sources[tn].buffer[index]))
               {
                   fprintf(stderr, "7 error5\n");
               }

               if(Get_ts(sources[i].buffer[j]) != Get_qs(sources[tn].buffer[index]))
               {
                   fprintf(stderr, "error6, Get_ts(%d, %d)=%d, Get_qs(%d, %d)=%d, rev: %d\n", 
                   i, j, Get_ts(sources[i].buffer[j]), 
                   tn, index, Get_qs(sources[tn].buffer[index]), sources[i].buffer[j].rev);
               }

               if(Get_te(sources[i].buffer[j]) != Get_qe(sources[tn].buffer[index]))
               {
                   fprintf(stderr, "error7, Get_te(%d, %d)=%d, Get_qe(%d, %d)=%d, rev: %d\n", 
                   i, j, Get_te(sources[i].buffer[j]), 
                   tn, index, Get_qe(sources[tn].buffer[index]), 
                   sources[i].buffer[j].rev);
               }

               

                if(sources[i].buffer[j].ml != sources[tn].buffer[index].ml)
                {
                    fprintf(stderr, "9 hahaha\n");
                }

            }
            
        }
    }

    fprintf(stderr, "total_reads:%d, total_overlaps: %d\n", total_reads, total_overlaps);
}



void debug_cut_ma_hit_t(ma_hit_t_alloc* sources, long long num_sources, ma_sub_t* coverage_cut)
{
    long long i, j, index;
    uint32_t qn, tn;
    long long total_overlaps = 0;
    long long total_reads = 0;

    for (i = 0; i < num_sources; i++)
    {
        for (j = 0; j < sources[i].length; j++)
        {
            if(Get_qn(sources[i].buffer[j]) != i)
            {
                fprintf(stderr, "error 2\n");
            }
            total_overlaps += sources[i].length;
            
        }

        if(sources[i].length != 0)
        {
            total_reads++;
        }

        if((sources[i].length == 0 && coverage_cut[i].del == 1)
            ||
            (sources[i].length != 0 && coverage_cut[i].del == 0))
        {
            ;
        }
        else
        {
            fprintf(stderr, "i: %d, sources[i].length: %d, coverage_cut[i].del: %d\n", i, sources[i].length, coverage_cut[i].del);
        }
        
    }

    for (i = 0; i < num_sources; i++)
    {
        for (j = 0; j < sources[i].length; j++)
        {
            qn = sources[i].buffer[j].qns>>32;
            tn = sources[i].buffer[j].tn;

            long long k;
            for (k = 0; k < sources[i].length; k++)
            {
                //here can be improved, since ma_hit_t_alloc has been sorted by tn
                if(sources[i].buffer[k].tn == tn 
                &&
                ((uint32_t)(sources[i].buffer[k].qns>>32)) == qn)
                {
                    if(k != j)
                    {
                        fprintf(stderr, "ERROR\n");
                    }
                }
            }
        }
    }










    for (i = 0; i < num_sources; i++)
    {
        for (j = 0; j < sources[i].length; j++)
        {
            /**
            if(sources[i].buffer[j].bl != sources[i].buffer[j].ml)
            {
                fprintf(stderr, "error1\n");
            }
            **/

            qn = sources[i].buffer[j].qns>>32;
            tn = sources[i].buffer[j].tn;

            index = get_specific_overlap(&(sources[tn]), tn, qn);


            if(index == -1)
            {
                fprintf(stderr, "error3\n");
            }
            else
            {
                if(sources[i].buffer[j].rev != sources[tn].buffer[index].rev)
                {
                    fprintf(stderr, "hahaha\n");
                }

               if(Get_qn(sources[i].buffer[j]) != Get_tn(sources[tn].buffer[index]))
               {
                   fprintf(stderr, "error4\n");
               }

               if(Get_tn(sources[i].buffer[j]) != Get_qn(sources[tn].buffer[index]))
               {
                   fprintf(stderr, "error5\n");
               }

               if(Get_ts(sources[i].buffer[j]) != Get_qs(sources[tn].buffer[index]))
               {
                   fprintf(stderr, "\nerror6, Get_ts(%d, %d)=%d, Get_qs(%d, %d)=%d, rev: %d\n", 
                   i, j, Get_ts(sources[i].buffer[j]), 
                   tn, index, Get_qs(sources[tn].buffer[index]), 
                   sources[i].buffer[j].rev);
               }

               if(Get_te(sources[i].buffer[j]) != Get_qe(sources[tn].buffer[index]))
               {
                   fprintf(stderr, "\nerror7, Get_te(%d, %d)=%d, Get_qe(%d, %d)=%d, rev: %d\n", 
                   i, j, Get_te(sources[i].buffer[j]), 
                   tn, index, Get_qe(sources[tn].buffer[index]), 
                   sources[i].buffer[j].rev);
               }

            }
            
        }
    }

    fprintf(stderr, "total_reads:%d, total_overlaps: %d\n", total_reads, total_overlaps);
}

/**********************************
 * Filter short potential unitigs *
 **********************************/
#define ASG_ET_MERGEABLE 0
#define ASG_ET_TIP       1
#define ASG_ET_MULTI_OUT 2
#define ASG_ET_MULTI_NEI 3

static inline int asg_is_utg_end(const asg_t *g, uint32_t v, uint64_t *lw)
{
    
	/**
	..............................
	.  w1---------------         .
	.    w2--------------        .
	.       w3--------------     .--->asg_arc_a(g, v^1)
	.          w4-------------   .
	.             w5------------ .
	..............................
						v---------------
						..............................
						.  w1---------------         .
						.    w2--------------        .
						.       w3--------------     .--->asg_arc_a(g, v)
						.          w4-------------   .
						.             w5------------ .
						..............................
    !!!!!note here the graph has already been cleaned by transitive reduction, so idealy:

        .........................
        .    w1---------------  .--->asg_arc_a(g, v^1)
        .........................
                            v---------------
                                    .........................
                                    .    w5---------------  .--->asg_arc_a(g, v)
                                    .........................

	**/
	///v^1 is the another direction of v
	uint32_t w, nv, nw, nw0, nv0 = asg_arc_n(g, v^1);
	int i, i0 = -1;
	asg_arc_t *aw, *av = asg_arc_a(g, v^1);

	///if this arc has not been deleted
	for (i = nv = 0; i < nv0; ++i)
		if (!av[i].del) i0 = i, ++nv;

	///see the example below
	if (nv == 0) return ASG_ET_TIP; // tip

    /**
        since the graph has already been cleaned by transitive reduction,
        w1 and w2 should not be overlapped with each other
        that mean v has mutiple in-edges, and each of them is not overlapped with others
        .........................
        .    w2---------------  .--->asg_arc_a(g, v^1)
        .    w1---------------  .
        .........................
                            v---------------

	**/
	if (nv > 1) return ASG_ET_MULTI_OUT; // multiple outgoing arcs



	
    /**
     * ///until here, nv == 1
        note the graph has already been cleaned by transitive reduction,
        .........................
        .    w1---------------  .--->asg_arc_a(g, v^1)
        .........................
                            v---------------
	**/
	/**
	p->ul: |____________31__________|__________1___________|______________32_____________|
	                    qn            direction of overlap       length of this node (not overlap length)
						                 (based on query)
	p->v : |___________31___________|__________1___________|
				        tn             reverse direction of overlap
						                  (based on target)
	p->ol: overlap length
	**/
	///until here, nv == 1
	if (lw) *lw = av[i0].ul<<32 | av[i0].v;
	/**
	p->ul: |____________31__________|__________1___________|______________32_____________|
	                    qns            direction of overlap       length of this node (not overlap length)
						                 (based on query)
	p->v : |___________31___________|__________1___________|
				        tns             reverse direction of overlap
						                  (based on target)
	p->ol: overlap length
	**/
	w = av[i0].v ^ 1;
	nw0 = asg_arc_n(g, w);
	aw = asg_arc_a(g, w);
	for (i = nw = 0; i < nw0; ++i)
		if (!aw[i].del) ++nw;


    /**
        note nw is at least 1, since we have v 
        nw > 1 means
        .........................
        . av[i0].v^1----------  .--->asg_arc_a(g, v^1)
        .........................
                        v---------------
                        w---------------
                        z---------------
        asg_arc_a(av[i0].v^1) is the (v, w, z), and v, w, z are not overlapped with each others
	**/
	if (nw != 1) return ASG_ET_MULTI_NEI;

	/**
     *  nw == 1 means
        note the graph has already been cleaned by transitive reduction,
        .........................
        .    w1---------------  .--->asg_arc_a(g, v^1)
        .........................
                            v---------------
                                    .........................
                                    .    w5---------------  .--->asg_arc_a(g, v)
                                    .........................

	**/
	return ASG_ET_MERGEABLE;
}


int asg_extend(const asg_t *g, uint32_t v, int max_ext, asg64_v *a)
{
	int ret;
	uint64_t lw;
	a->n = 0;
	kv_push(uint64_t, *a, v);
	do {
        /**
		 note that here the graph has been cleanned by transitive reduction
		 the following first line is to find the prefix of v^1:
		         (v^1)--->()---->()---->()----->....
         that is the suffix of v:
             ...>(v)
        **/
		ret = asg_is_utg_end(g, v^1, &lw);
        /**
        #define ASG_ET_MERGEABLE 0
        #define ASG_ET_TIP       1
        #define ASG_ET_MULTI_OUT 2
        #define ASG_ET_MULTI_NEI 3
        **/
		if (ret != 0) break;
		kv_push(uint64_t, *a, lw);
		/**
		 ret == 0 means: 
		   v^1 and is the only prefix of (uint32_t)lw, 
		   and (uint32_t)lw is the only prefix of v^1
		**/
		
		v = (uint32_t)lw;
	} while (--max_ext > 0);
	return ret;
}


static inline int asg_is_single_edge(const asg_t *g, uint32_t v, uint32_t start_node)
{
    
	/**
	..............................
	.  w1---------------         .
	.    w2--------------        .
	.       w3--------------     .--->asg_arc_a(g, v^1)
	.          w4-------------   .
	.             w5------------ .
	..............................
						v---------------
						..............................
						.  w1---------------         .
						.    w2--------------        .
						.       w3--------------     .--->asg_arc_a(g, v)
						.          w4-------------   .
						.             w5------------ .
						..............................
    !!!!!note here the graph has already been cleaned by transitive reduction, so idealy:

        .........................
        .    w1---------------  .--->asg_arc_a(g, v^1)
        .........................
                            v---------------
                                    .........................
                                    .    w5---------------  .--->asg_arc_a(g, v)
                                    .........................

	**/
	///v^1 is the another direction of v
	uint32_t w, nv, nw, nw0, nv0 = asg_arc_n(g, v^1);
	int i, i0 = -1;
	asg_arc_t *av = asg_arc_a(g, v^1);

    int flag = 0;
	///if this arc has not been deleted
	for (i = nv = 0; i < nv0; ++i)
    {
        ///if (!av[i].del)
        {
            i0 = i;
            ++nv;
            if(av[i].v>>1 == start_node)
            {
                flag = 1;
            }
        } 
    }
		
    if(flag == 0)
    {
        fprintf(stderr, "****ERROR\n");
    }

    return nv;
}


asg_t *ma_sg_gen(const ma_hit_t_alloc* sources, long long n_read, const ma_sub_t *coverage_cut, 
int max_hang, int min_ovlp)
{
    double startTime = Get_T();
	size_t i, j;
	asg_t *g;
	///just calloc
	g = asg_init();

	///add seq to graph, seq just save the length of each read
	for (i = 0; i < n_read; ++i) 
    {
        ///if a read has been deleted, should we still add them?
		asg_seq_set(g, i, coverage_cut[i].e - coverage_cut[i].s, coverage_cut[i].del);
        g->seq[i].c = coverage_cut[i].c;
	}

    g->seq_vis = (uint8_t*)calloc(g->n_seq*2, sizeof(uint8_t));

    for (i = 0; i < n_read; ++i) 
    {
        for (j = 0; j < sources[i].length; j++)
        {
            int r;
		    asg_arc_t t, *p;
		    const ma_hit_t *h = &(sources[i].buffer[j]);
            //high coverage region [sub[qn].e, sub[qn].s) in query
            int ql = coverage_cut[Get_qn(*h)].e - coverage_cut[Get_qn(*h)].s;
            //high coverage region [sub[qn].e, sub[qn].s) in target
            int tl = coverage_cut[Get_tn(*h)].e - coverage_cut[Get_tn(*h)].s;
		    r = ma_hit2arc(h, ql, tl, max_hang, MAX_HANG_PRE, min_ovlp, &t);
            /**
            #define MA_HT_INT       (-1)
            #define MA_HT_QCONT      (-2)
            #define MA_HT_TCONT      (-3)
            #define MA_HT_SHORT_OVLP (-4)
            the short overlaps and the overlaps with contain reads have already been removed
            here we should have overhang
            so r should always >= 0
            **/
            if (r >= 0) 
            {
                ///push node?
                p = asg_arc_pushp(g);
                *p = t;
		    }
            else
            {
                fprintf(stderr, "error\n");
            }
        }
    }

    asg_cleanup(g);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

	return g;
}


// pop bubbles from vertex v0; the graph MJUST BE symmetric: if u->v present, v'->u' must be present as well
//note!!!!!!!! here we don't exculde the deleted edges
static uint64_t asg_bub_finder_with_del(asg_t *g, uint32_t v0, int max_dist, buf_t *b,
uint32_t cut_in_node)
{
	uint32_t i, n_pending = 0;
	uint64_t n_pop = 0;
	///if this node has been deleted
	if (g->seq[v0>>1].del) return 0; // already deleted
	///asg_arc_n(n0)
	if ((uint32_t)g->idx[v0] < 2) return 0; // no bubbles
	///S saves nodes with all incoming edges visited
	b->S.n = b->T.n = b->b.n = b->e.n = 0;
	///for each node, b->a saves all related information
	b->a[v0].c = b->a[v0].d = 0;
	///b->S is the nodes with all incoming edges visited
	kv_push(uint32_t, b->S, v0);

	do {
		///v is a node that all incoming edges have been visited
		///d is the distance from v0 to v
		uint32_t v = kv_pop(b->S), d = b->a[v].d, c = b->a[v].c;
		uint32_t nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		///why we have this assert?
        /****************************may have bugs********************************/
		///assert(nv > 0);
        /****************************may have bugs********************************/

		///all out-edges of v
		for (i = 0; i < nv; ++i) { // loop through v's neighbors
			/**
			p->ul: |____________31__________|__________1___________|______________32_____________|
								qn            direction of overlap       length of this node (not overlap length)
												(in the view of query)
			p->v : |___________31___________|__________1___________|
								tn             reverse direction of overlap
											(in the view of target)
			p->ol: overlap length
			**/
			uint32_t w = av[i].v, l = (uint32_t)av[i].ul; // v->w with length l
			binfo_t *t = &b->a[w];
			///that means there is a circle, directly terminate the whole bubble poping
			if (w == v0)
            {
                //fprintf(stderr, "n_pop error1\n");
                goto pop_reset;
            } 
            
			///if this edge has been deleted
            /****************************may have bugs********************************/
			///if (av[i].del) continue;
            /****************************may have bugs********************************/

			///push the edge
			kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);

			///find a too far path? directly terminate the whole bubble poping
			if (d + l > max_dist)
            {
                //fprintf(stderr, "n_pop error2\n");
                break; // too far
            } 
            


			if (t->s == 0) { // this vertex has never been visited
				kv_push(uint32_t, b->b, w); // save it for revert
				///t->p means the in-node of w is v
				///t->s = 1 means w has been visited
				///d is len(v0->v), l is len(v->w), so t->d is len(v0->w)
				t->p = v, t->s = 1, t->d = d + l;
				///incoming edges of w
				t->r = count_out_with_del(g, w^1);
                if((w>>1) == cut_in_node)
                {
                    t->r--;
                }
				++n_pending;
			} else { // visited before
				///c seems the max weight of node
				if (c + 1 > t->c || (c + 1 == t->c && d + l > t->d)) t->p = v;
				if (c + 1 > t->c) t->c = c + 1;
				///update len(v0->w)
				if (d + l < t->d) t->d = d + l; // update dist
			}
            /****************************may have bugs********************************/
			///assert(t->r > 0);
            /****************************may have bugs********************************/
			//if all incoming edges of w have visited
			//push it to b->S
			if (--(t->r) == 0) {
				uint32_t x = asg_arc_n(g, w);
				if (x) kv_push(uint32_t, b->S, w);
				else kv_push(uint32_t, b->T, w); // a tip
				--n_pending;
			}
		}
		///if i < nv, that means (d + l > max_dist)
		if (i < nv || b->S.n == 0)
        {
            ///fprintf(stderr, "n_pop error3\n");
            goto pop_reset;
        } 
        
	} while (b->S.n > 1 || n_pending);
	///asg_bub_backtrack(g, v0, b);
	///n_pop = 1 | (uint64_t)b->T.n<<32;
    n_pop = 1;
pop_reset:
	for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
		binfo_t *t = &b->a[b->b.a[i]];
		t->s = t->c = t->d = 0;
	}
    ///fprintf(stderr, "n_pop: %d\n", n_pop);
	return n_pop;
}


int asg_arc_del_triangular(asg_t *g, long long max_dist)
{
    ///the reason is that each read has two direction (query->target, target->query)
	uint32_t v, w, n_vtx = g->n_seq * 2, n_reduced = 0;
    buf_t b;
	if (!g->is_symm) asg_symm(g);
	memset(&b, 0, sizeof(buf_t));
	///set information for each node
	b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    int flag0, flag1, node;
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if (g->seq[v>>1].del)
        {
            continue;
        } 

        if(nv != 2)
        {
            continue;
        }

        /**********************test first node************************/
        flag0 = asg_is_single_edge(g, av[0].v, v>>1);
        flag1 = asg_is_single_edge(g, av[1].v, v>>1);
        if(flag0 == flag1)
        {
            continue;
        }
        if(flag0 < 1 || flag0 > 2)
        {
            continue;
        }
        if(flag1 < 1 || flag1 > 2)
        {
            continue;
        }
        if(flag0 == 2)
        {
            node = 0;
        }
        else if(flag1 == 2)
        {
            node = 1;
        }
       /**********************test first node************************/



        /**********************test second node************************/
        w = av[node].v^1;
        asg_arc_t *aw = asg_arc_a(g, w);
        uint32_t nw = asg_arc_n(g, w);
        flag0 = asg_is_single_edge(g, aw[0].v, w>>1);
        flag1 = asg_is_single_edge(g, aw[1].v, w>>1);
        if(flag0 == flag1)
        {
            continue;
        }
        if(flag0 < 1 || flag0 > 2)
        {
            continue;
        }
        if(flag1 < 1 || flag1 > 2)
        {
            continue;
        }


        if(flag0 == 2 && (v>>1) != (aw[0].v>>1))
        {
            fprintf(stderr, "error 0\n");
        }


        if(flag1 == 2 && (v>>1) != (aw[1].v>>1))
        {
            fprintf(stderr, "error 1\n");
        }
        
        

        /**********************test second node************************/
        
        ///if not a bubble
        if(asg_bub_finder_with_del(g, v, max_dist, &b, (w>>1)) == 0)
        {
            continue;
        }
        
        
        av[node].del = 1;
        ///remove the reverse direction
		asg_arc_del(g, av[node].v^1, av[node].ul>>32^1, 1);
        n_reduced++;
        
    }

    free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);

    if (n_reduced) {
		asg_cleanup(g);
		asg_symm(g);
	}

    fprintf(stderr, "[M::%s] removed %d short overlaps\n", __func__, n_reduced);
}



// pop bubbles from vertex v0; the graph MJUST BE symmetric: if u->v present, v'->u' must be present as well
//note!!!!!!!! here we don't exculde the deleted edges
static uint64_t asg_bub_finder_with_del_advance(asg_t *g, uint32_t v0, int max_dist, buf_t *b)
{
	uint32_t i, n_pending = 0;
	uint64_t n_pop = 0;
	///if this node has been deleted
	if (g->seq[v0>>1].del) return 0; // already deleted
	///asg_arc_n(n0)
	if ((uint32_t)g->idx[v0] < 2) return 0; // no bubbles
	///S saves nodes with all incoming edges visited
	b->S.n = b->T.n = b->b.n = b->e.n = 0;
	///for each node, b->a saves all related information
	b->a[v0].c = b->a[v0].d = 0;
	///b->S is the nodes with all incoming edges visited
	kv_push(uint32_t, b->S, v0);

	do {
		///v is a node that all incoming edges have been visited
		///d is the distance from v0 to v
		uint32_t v = kv_pop(b->S), d = b->a[v].d, c = b->a[v].c;
		uint32_t nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		///why we have this assert?
        /****************************may have bugs********************************/
		///assert(nv > 0);
        /****************************may have bugs********************************/

		///all out-edges of v
		for (i = 0; i < nv; ++i) { // loop through v's neighbors
			/**
			p->ul: |____________31__________|__________1___________|______________32_____________|
								qn            direction of overlap       length of this node (not overlap length)
												(in the view of query)
			p->v : |___________31___________|__________1___________|
								tn             reverse direction of overlap
											(in the view of target)
			p->ol: overlap length
			**/
			uint32_t w = av[i].v, l = (uint32_t)av[i].ul; // v->w with length l
			binfo_t *t = &b->a[w];
			///that means there is a circle, directly terminate the whole bubble poping
			if (w == v0)
            {
                //fprintf(stderr, "n_pop error1\n");
                goto pop_reset;
            } 
            
			///if this edge has been deleted
            /****************************may have bugs********************************/
			///if (av[i].del) continue;
            /****************************may have bugs********************************/

			///push the edge
			kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);

			///find a too far path? directly terminate the whole bubble poping
			if (d + l > max_dist)
            {
                //fprintf(stderr, "n_pop error2\n");
                break; // too far
            } 
            


			if (t->s == 0) { // this vertex has never been visited
				kv_push(uint32_t, b->b, w); // save it for revert
				///t->p means the in-node of w is v
				///t->s = 1 means w has been visited
				///d is len(v0->v), l is len(v->w), so t->d is len(v0->w)
				t->p = v, t->s = 1, t->d = d + l;
				///incoming edges of w
				t->r = count_out_with_del(g, w^1);
				++n_pending;
			} else { // visited before
				///c seems the max weight of node
				if (c + 1 > t->c || (c + 1 == t->c && d + l > t->d)) t->p = v;
				if (c + 1 > t->c) t->c = c + 1;
				///update len(v0->w)
				if (d + l < t->d) t->d = d + l; // update dist
			}
            /****************************may have bugs********************************/
			///assert(t->r > 0);
            /****************************may have bugs********************************/
			//if all incoming edges of w have visited
			//push it to b->S
			if (--(t->r) == 0) {
				uint32_t x = asg_arc_n(g, w);
				if (x) kv_push(uint32_t, b->S, w);
				else kv_push(uint32_t, b->T, w); // a tip
				--n_pending;
			}
		}
		///if i < nv, that means (d + l > max_dist)
		if (i < nv || b->S.n == 0)
        {
            ///fprintf(stderr, "n_pop error3\n");
            goto pop_reset;
        } 
        
	} while (b->S.n > 1 || n_pending);
	///asg_bub_backtrack(g, v0, b);
	///n_pop = 1 | (uint64_t)b->T.n<<32;
    n_pop = 1;
pop_reset:
	for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
		binfo_t *t = &b->a[b->b.a[i]];
		t->s = t->c = t->d = 0;
	}
    ///fprintf(stderr, "n_pop: %d\n", n_pop);
	return n_pop;
}


// pop bubbles from vertex v0; the graph MJUST BE symmetric: if u->v present, v'->u' must be present as well
//note!!!!!!!! here we don't exculde the deleted edges
static uint64_t asg_bub_finder_without_del_advance(asg_t *g, uint32_t v0, int max_dist, 
buf_t *b)
{
	uint32_t i, n_pending = 0;
	uint64_t n_pop = 0;
	///if this node has been deleted
	if (g->seq[v0>>1].del) return 0; // already deleted
	///asg_arc_n(n0)
	if ((uint32_t)g->idx[v0] < 2) return 0; // no bubbles
    if(count_out_without_del(g, v0) < 2) return 0; // no bubbles


	///S saves nodes with all incoming edges visited
	b->S.n = b->T.n = b->b.n = b->e.n = 0;
	///for each node, b->a saves all related information
	b->a[v0].c = b->a[v0].d = 0;
	///b->S is the nodes with all incoming edges visited
	kv_push(uint32_t, b->S, v0);

	do {
		///v is a node that all incoming edges have been visited
		///d is the distance from v0 to v
		uint32_t v = kv_pop(b->S), d = b->a[v].d, c = b->a[v].c;
		uint32_t nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		///why we have this assert?
        /****************************may have bugs********************************/
		///assert(nv > 0);
        /****************************may have bugs********************************/

		///all out-edges of v
		for (i = 0; i < nv; ++i) { // loop through v's neighbors
			/**
			p->ul: |____________31__________|__________1___________|______________32_____________|
								qn            direction of overlap       length of this node (not overlap length)
												(in the view of query)
			p->v : |___________31___________|__________1___________|
								tn             reverse direction of overlap
											(in the view of target)
			p->ol: overlap length
			**/
			uint32_t w = av[i].v, l = (uint32_t)av[i].ul; // v->w with length l
			binfo_t *t = &b->a[w];
			///that means there is a circle, directly terminate the whole bubble poping
			if (w == v0)
            {
                //fprintf(stderr, "n_pop error1\n");
                goto pop_reset;
            } 
            
			///if this edge has been deleted
            /****************************may have bugs********************************/
			if (av[i].del) continue;
            /****************************may have bugs********************************/

			///push the edge
			kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);

			///find a too far path? directly terminate the whole bubble poping
			if (d + l > max_dist)
            {
                //fprintf(stderr, "n_pop error2\n");
                break; // too far
            } 
            


			if (t->s == 0) { // this vertex has never been visited
				kv_push(uint32_t, b->b, w); // save it for revert
				///t->p means the in-node of w is v
				///t->s = 1 means w has been visited
				///d is len(v0->v), l is len(v->w), so t->d is len(v0->w)
				t->p = v, t->s = 1, t->d = d + l;
				///incoming edges of w
				t->r = count_out_without_del(g, w^1);
				++n_pending;
			} else { // visited before
				///c seems the max weight of node
				if (c + 1 > t->c || (c + 1 == t->c && d + l > t->d)) t->p = v;
				if (c + 1 > t->c) t->c = c + 1;
				///update len(v0->w)
				if (d + l < t->d) t->d = d + l; // update dist
			}
            /****************************may have bugs********************************/
			///assert(t->r > 0);
            /****************************may have bugs********************************/
			//if all incoming edges of w have visited
			//push it to b->S
			if (--(t->r) == 0) {
				uint32_t x = asg_arc_n(g, w);
				if (x) kv_push(uint32_t, b->S, w);
				else kv_push(uint32_t, b->T, w); // a tip
				--n_pending;
			}
		}
		///if i < nv, that means (d + l > max_dist)
		if (i < nv || b->S.n == 0)
        {
            ///fprintf(stderr, "n_pop error3\n");
            goto pop_reset;
        } 
        
	} while (b->S.n > 1 || n_pending);
	///asg_bub_backtrack(g, v0, b);
	///n_pop = 1 | (uint64_t)b->T.n<<32;
    n_pop = 1;
pop_reset:
	for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
		binfo_t *t = &b->a[b->b.a[i]];
		t->s = t->c = t->d = 0;
	}
    ///fprintf(stderr, "n_pop: %d\n", n_pop);
	return n_pop;
}


// pop bubbles from vertex v0; the graph MJUST BE symmetric: if u->v present, v'->u' must be present as well
//note!!!!!!!! here we don't exculde the deleted edges
static uint64_t asg_bub_end_finder_with_del_advance(asg_t *g, uint32_t* v_Ns, uint32_t occ, 
int max_dist, buf_t *b, uint32_t exculde_init, uint32_t exclude_node, uint32_t* sink)
{
	uint32_t i, j, n_pending = 0;
	uint64_t n_pop = 0;
    ///S saves nodes with all incoming edges visited
    b->S.n = b->T.n = b->b.n = b->e.n = 0;
    for (j = 0; j < occ; j++)
    {
        ///if this node has been deleted
        if (g->seq[v_Ns[j]>>1].del) return 0; // already deleted
        ///for each node, b->a saves all related information
        b->a[v_Ns[j]].c = b->a[v_Ns[j]].d = 0;
        ///b->S is the nodes with all incoming edges visited
        kv_push(uint32_t, b->S, (v_Ns[j]<<1)|exculde_init);
    }
    

	
	do {
		///v is a node that all incoming edges have been visited
		///d is the distance from v0 to v
		uint32_t v = kv_pop(b->S), f = v & (uint32_t)1;
        v = v >> 1;
        uint32_t d = b->a[v].d, c = b->a[v].c;


		uint32_t nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		
		///all out-edges of v
		for (i = 0; i < nv; ++i) { // loop through v's neighbors
			/**
			p->ul: |____________31__________|__________1___________|______________32_____________|
								qn            direction of overlap       length of this node (not overlap length)
												(in the view of query)
			p->v : |___________31___________|__________1___________|
								tn             reverse direction of overlap
											(in the view of target)
			p->ol: overlap length
			**/
			uint32_t w = av[i].v, l = (uint32_t)av[i].ul; // v->w with length l
			binfo_t *t = &b->a[w];

            for (j = 0; j < occ; j++)
            {
                if(w == v_Ns[j]) goto pop_reset;
            }
            


            if(f && (exclude_node) == (w)) continue;

            ///if (av[i].del) continue;

			///push the edge
			kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);

			///find a too far path? directly terminate the whole bubble poping
			if (d + l > max_dist)
            {
                break; // too far
            } 
            


			if (t->s == 0) { // this vertex has never been visited
				kv_push(uint32_t, b->b, w); // save it for revert
				///t->p means the in-node of w is v
				///t->s = 1 means w has been visited
				///d is len(v0->v), l is len(v->w), so t->d is len(v0->w)
				t->p = v, t->s = 1, t->d = d + l;
				///incoming edges of w
				t->r = count_out_with_del(g, w^1);
                ///t->r = count_out_without_del(g, w^1);
				++n_pending;
			} else { // visited before
				///c seems the max weight of node
				if (c + 1 > t->c || (c + 1 == t->c && d + l > t->d)) t->p = v;
				if (c + 1 > t->c) t->c = c + 1;
				///update len(v0->w)
				if (d + l < t->d) t->d = d + l; // update dist
			}
            /****************************may have bugs********************************/
			///assert(t->r > 0);
            /****************************may have bugs********************************/
			//if all incoming edges of w have visited
			//push it to b->S
			if (--(t->r) == 0) {
				uint32_t x = asg_arc_n(g, w);
				//if (x) kv_push(uint32_t, b->S, w);
                if (x) kv_push(uint32_t, b->S, w<<1);
				else kv_push(uint32_t, b->T, w); // a tip
				--n_pending;
			}
		}
		///if i < nv, that means (d + l > max_dist)
		if (i < nv || b->S.n == 0)
        {
            goto pop_reset;
        } 
        
	} while (b->S.n > 1 || n_pending);

    (*sink) = b->S.a[0]>>1;
	///asg_bub_backtrack(g, v0, b);
	///n_pop = 1 | (uint64_t)b->T.n<<32;
    n_pop = 1;
pop_reset:
	for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
		binfo_t *t = &b->a[b->b.a[i]];
		t->s = t->c = t->d = 0;
	}
	return n_pop;
}

// pop bubbles from vertex v0; the graph MJUST BE symmetric: if u->v present, v'->u' must be present as well
//note!!!!!!!! here we don't exculde the deleted edges
static uint64_t asg_bub_end_finder_with_del_advance_debug(asg_t *g, uint32_t* v_Ns, uint32_t occ, 
int max_dist, buf_t *b, uint32_t exculde_init, uint32_t exclude_node)
{
	uint32_t i, j, n_pending = 0;
	uint64_t n_pop = 0;
    ///S saves nodes with all incoming edges visited
    b->S.n = b->T.n = b->b.n = b->e.n = 0;
    for (j = 0; j < occ; j++)
    {
        ///if this node has been deleted
        if (g->seq[v_Ns[j]>>1].del) return 0; // already deleted
        ///for each node, b->a saves all related information
        b->a[v_Ns[j]].c = b->a[v_Ns[j]].d = 0;
        ///b->S is the nodes with all incoming edges visited
        kv_push(uint32_t, b->S, (v_Ns[j]<<1)|exculde_init);

        fprintf(stderr, "init: %.*s\n", Get_NAME_LENGTH(R_INF, v_Ns[j]>>1), 
        Get_NAME(R_INF, v_Ns[j]>>1));
    }
    

	
	do {
		///v is a node that all incoming edges have been visited
		///d is the distance from v0 to v
		uint32_t v = kv_pop(b->S), f = v & (uint32_t)1;
        v = v >> 1;
        uint32_t d = b->a[v].d, c = b->a[v].c;


		uint32_t nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);

        fprintf(stderr, "v: %.*s, nv: %d, f: %d\n", 
        Get_NAME_LENGTH(R_INF, v>>1), Get_NAME(R_INF, v>>1), nv, f);
		
		///all out-edges of v
		for (i = 0; i < nv; ++i) { // loop through v's neighbors
			/**
			p->ul: |____________31__________|__________1___________|______________32_____________|
								qn            direction of overlap       length of this node (not overlap length)
												(in the view of query)
			p->v : |___________31___________|__________1___________|
								tn             reverse direction of overlap
											(in the view of target)
			p->ol: overlap length
			**/
			uint32_t w = av[i].v, l = (uint32_t)av[i].ul; // v->w with length l
			binfo_t *t = &b->a[w];

            fprintf(stderr, "w: %.*s\n", Get_NAME_LENGTH(R_INF, w>>1), 
            Get_NAME(R_INF, w>>1));

            if(f && (exclude_node) == (w)) 
            {
                fprintf(stderr, "exclude_node: %.*s\n", Get_NAME_LENGTH(R_INF, exclude_node>>1), 
            Get_NAME(R_INF, exclude_node>>1));
                continue;
            }

            ///if (av[i].del) continue;

			///push the edge
			kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);

			///find a too far path? directly terminate the whole bubble poping
			if (d + l > max_dist)
            {
                break; // too far
            } 
            


			if (t->s == 0) { // this vertex has never been visited
				kv_push(uint32_t, b->b, w); // save it for revert
				///t->p means the in-node of w is v
				///t->s = 1 means w has been visited
				///d is len(v0->v), l is len(v->w), so t->d is len(v0->w)
				t->p = v, t->s = 1, t->d = d + l;
				///incoming edges of w
				t->r = count_out_with_del(g, w^1);
                ///t->r = count_out_without_del(g, w^1);
				++n_pending;
			} else { // visited before
				///c seems the max weight of node
				if (c + 1 > t->c || (c + 1 == t->c && d + l > t->d)) t->p = v;
				if (c + 1 > t->c) t->c = c + 1;
				///update len(v0->w)
				if (d + l < t->d) t->d = d + l; // update dist
			}
            /****************************may have bugs********************************/
			///assert(t->r > 0);
            /****************************may have bugs********************************/
			//if all incoming edges of w have visited
			//push it to b->S
			if (--(t->r) == 0) {
				uint32_t x = asg_arc_n(g, w);
				//if (x) kv_push(uint32_t, b->S, w);
                if (x) kv_push(uint32_t, b->S, w<<1);
				else kv_push(uint32_t, b->T, w); // a tip
				--n_pending;
			}
		}
		///if i < nv, that means (d + l > max_dist)
		if (i < nv || b->S.n == 0)
        {
            goto pop_reset;
        } 
        
	} while (b->S.n > 1 || n_pending);
	///asg_bub_backtrack(g, v0, b);
	///n_pop = 1 | (uint64_t)b->T.n<<32;
    n_pop = 1;
pop_reset:
	for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
		binfo_t *t = &b->a[b->b.a[i]];
		t->s = t->c = t->d = 0;
	}
	return n_pop;
}



int if_node_exist(uint32_t* nodes, uint32_t length, uint32_t query)
{
    uint32_t i;
    for (i = 0; i < length; ++i)
    {
        if((nodes[i]>>1) == query)
        {
            return 1;
        }
    }

    return 0;
}


int test_triangular(asg_t *g, uint32_t* nodes, uint32_t length, 
uint32_t startNode, uint32_t endNode)
{
    
    uint32_t i, v, w;
    int flag0, flag1, node;
    int n_reduced = 0;
    for (i = 0; i < length; ++i)
    {
        v = nodes[i];

        if((v>>1) == (startNode>>1) || (v>>1) == (endNode>>1))
        {
            continue;
        }

        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if(nv != 2)
        {
            continue;
        }
        /**********************test first node************************/
        flag0 = asg_is_single_edge(g, av[0].v, v>>1);
        flag1 = asg_is_single_edge(g, av[1].v, v>>1);
        if(flag0 == flag1)
        {
            continue;
        }
        if(flag0 < 1 || flag0 > 2)
        {
            continue;
        }
        if(flag1 < 1 || flag1 > 2)
        {
            continue;
        }
        if(flag0 == 2)
        {
            node = 0;
        }
        else if(flag1 == 2)
        {
            node = 1;
        }
       /**********************test first node************************/



        /**********************test second node************************/
        w = av[node].v^1;
        asg_arc_t *aw = asg_arc_a(g, w);
        uint32_t nw = asg_arc_n(g, w);
        flag0 = asg_is_single_edge(g, aw[0].v, w>>1);
        flag1 = asg_is_single_edge(g, aw[1].v, w>>1);
        if(flag0 == flag1)
        {
            continue;
        }
        if(flag0 < 1 || flag0 > 2)
        {
            continue;
        }
        if(flag1 < 1 || flag1 > 2)
        {
            continue;
        }


        if(flag0 == 2 && (v>>1) != (aw[0].v>>1))
        {
            fprintf(stderr, "error 0\n");
        }


        if(flag1 == 2 && (v>>1) != (aw[1].v>>1))
        {
            fprintf(stderr, "error 1\n");
        }
        /**********************test second node************************/

        if(if_node_exist(nodes, length, (w>>1)) && ((w>>1) != (endNode>>1)))
        {
            av[node].del = 1;
            ///remove the reverse direction
            asg_arc_del(g, av[node].v^1, av[node].ul>>32^1, 1);
            n_reduced++;

            ////fprintf(stderr, "v>>1: %u, w>>1: %u\n", v>>1, w>>1);
        }

    }

    return n_reduced;
}

long long single_edge_length(asg_t *g, uint32_t begNode, uint32_t endNode, long long edgeLen)
{
    
    uint32_t v = begNode;
    uint32_t nv = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v);
    long long rLen = 0;

    while (rLen < edgeLen && nv == 1)
    {
        rLen++;
        if((av[0].v>>1) == endNode)
        {
            return rLen;
        }

        if(asg_is_single_edge(g, av[0].v, v>>1) != 1)
        {
            return -1;
        }

        v = av[0].v;
        nv = asg_arc_n(g, v);
        av = asg_arc_a(g, v);      
    }

    return -1;   

}

uint32_t detect_single_path(asg_t *g, uint32_t begNode, uint32_t* endNode, long long* Len, buf_t* b)
{
    
    uint32_t v = begNode;
    uint32_t nv, rnv;
    asg_arc_t *av;
    (*Len) = 0;


    while (1)
    {
        (*Len)++;
        nv = asg_arc_n(g, v);
        av = asg_arc_a(g, v);   
        (*endNode) = v;

        if(b) kv_push(uint32_t, b->b, v>>1);

        if(nv == 0)
        {
            return END_TIPS;
        }

        if(nv == 2)
        {
            return TWO_OUTPUT;
        }

        if(nv > 2)
        {
            return MUL_OUTPUT;
        }

        ///up to here, nv=1
        ///rnv must >= 1
        rnv = asg_is_single_edge(g, av[0].v, v>>1);
        v = av[0].v;
        (*endNode) = v;
        if(rnv == 2)
        {
            (*Len)++;
            if(b) kv_push(uint32_t, b->b, v>>1);
            return TWO_INPUT;
        }

        if(rnv > 2)
        {
            (*Len)++;
            if(b) kv_push(uint32_t, b->b, v>>1);
            return MUL_INPUT;
        }   


        if((v>>1) == (begNode>>1))
        {
            return LOOP;
        } 
    }

    return LONG_TIPS;   
}

int detect_bubble_end(asg_t *g, uint32_t begNode1, uint32_t begNode2, uint32_t* endNode, 
long long* minLen, buf_t* b)
{
    uint32_t e1, e2;
    long long l1, l2;

    if(detect_single_path(g, begNode1, &e1, &l1, b) == TWO_INPUT 
       && 
       detect_single_path(g, begNode2, &e2, &l2, b) == TWO_INPUT)
    {
        if(e1 == e2)
        {
            (*endNode) = e1;
            (*minLen) = (l1 <= l2)? l1: l2;
            return 1;
        }
    }

    return 0;
}


int detect_simple_bubble(asg_t *g, uint32_t begNode, uint32_t* endNode, long long* minLen, buf_t* b)
{
    uint32_t e1, e2;
    long long l1, l2;

    if(asg_arc_n(g, begNode) != 2)
    {
        return 0;
    }

    if(asg_is_single_edge(g, asg_arc_a(g, begNode)[0].v, begNode>>1)!=1
      || 
       asg_is_single_edge(g, asg_arc_a(g, begNode)[1].v, begNode>>1)!=1)
    {
        return 0;
    }

    if(b) kv_push(uint32_t, b->b, begNode>>1);

    // if(begNode == 17939188)
    // {
    //     fprintf(stderr, "#begNode: %d, nv: %d\n", begNode, asg_arc_n(g, begNode));
    //     fprintf(stderr, "#asg_arc_a(g, begNode)[0].v: %d\n", asg_arc_a(g, begNode)[0].v);
    //     fprintf(stderr, "#asg_arc_a(g, begNode)[1].v: %d\n", asg_arc_a(g, begNode)[1].v);
    //     fflush(stderr);
    // }


    if(detect_single_path(g, asg_arc_a(g, begNode)[0].v, &e1, &l1, b) == TWO_INPUT 
       && 
       detect_single_path(g, asg_arc_a(g, begNode)[1].v, &e2, &l2, b) == TWO_INPUT)
    {
        if(e1 == e2)
        {
            (*endNode) = e1;
            (*minLen) = (l1 <= l2)? l1: l2;
            (*minLen)++;
            return 1;
        }
    }

    return 0;
}


uint32_t detect_single_path_with_single_bubbles(asg_t *g, uint32_t begNode, uint32_t* endNode, 
long long* Len, buf_t* b, uint32_t max_ext)
{
    
    uint32_t v = begNode;
    uint32_t nv, rnv;
    asg_arc_t *av;
    long long bLen;
    long long pre_b_n;
    (*Len) = 0;


    // if(begNode == 18046652)
    // {
    //     fprintf(stderr, "inner begNode: %d, vn: %d\n", begNode, asg_arc_n(g, begNode));
    //     fflush(stderr);
    // }

    // if(begNode == 18014006)
    // {
    //     fprintf(stderr, "inner begNode: %d, vn: %d\n", begNode, asg_arc_n(g, begNode));
    //     fflush(stderr);
    // }
    
    
    while (1)
    {
        nv = asg_arc_n(g, v);
        av = asg_arc_a(g, v);   
        (*endNode) = v;
        (*Len)++;

        // if(begNode == 18046652)
        // {
        //     fprintf(stderr, "(*Len): %d, v: %d, begNode: %d\n", (*Len), v, begNode);
        //     fflush(stderr);
        // }

        // if(begNode == 18014006)
        // {
        //     fprintf(stderr, "(*Len): %d, v: %d, begNode: %d\n", (*Len), v, begNode);
        //     fflush(stderr);
        // }

        if((*Len) > max_ext)
        {
            return LONG_TIPS_UNDER_MAX_EXT;
        }

        
        if(b) kv_push(uint32_t, b->b, v>>1);
        /**
        if(b && b->b.n > 1000000)
        {
            fprintf(stderr, "begNode>>1: %u, v>>1: %u, b->b.n: %u\n", begNode>>1, v>>1, b->b.n);
        }
        **/
 

        if(nv == 0)
        {
            return END_TIPS;
        }

        if(nv == 2)
        {

            // if(begNode == 18014006)
            // {
            //     fprintf(stderr, "bubble (*Len): %d, v: %d, begNode: %d\n", (*Len), v, begNode);
            //     fflush(stderr);
            // }

            if(b) pre_b_n = b->b.n;
            if(!detect_simple_bubble(g, v, &v, &bLen, b))
            {
                if(b) b->b.n = pre_b_n;
                return TWO_OUTPUT;
            }


            // if(begNode == 18014006)
            // {
            //     fprintf(stderr, "bubble (*Len): %d, v: %d, begNode: %d\n", (*Len), v, begNode);
            //     fflush(stderr);
            // }



            (*Len) = (*Len) + bLen - 2;
            continue;
        }

        if(nv > 2)
        {
            return MUL_OUTPUT;
        }

        ///up to here, nv=1
        ///rnv must >= 1
        rnv = asg_is_single_edge(g, av[0].v, v>>1);
        v = av[0].v;
        (*endNode) = v;
        if(rnv == 2)
        {
            if(b) kv_push(uint32_t, b->b, v>>1);
            (*Len)++;
            return TWO_INPUT;
        }

        if(rnv > 2)
        {
            if(b) kv_push(uint32_t, b->b, v>>1);
            (*Len)++;
            return MUL_INPUT;
        } 

        if((v>>1) == (begNode>>1))
        {
            return LOOP;
        }  
    }

    return LONG_TIPS;   
}

int detect_bubble_end_with_bubbles(asg_t *g, uint32_t begNode1, uint32_t begNode2, 
uint32_t* endNode, long long* minLen, buf_t* b)
{
    uint32_t e1, e2;
    long long l1, l2;

    if(detect_single_path_with_single_bubbles(g, begNode1, &e1, &l1, b, (uint32_t)-1) == TWO_INPUT 
       && 
       detect_single_path_with_single_bubbles(g, begNode2, &e2, &l2, b, (uint32_t)-1) == TWO_INPUT)
    {
        if(e1 == e2)
        {
            (*endNode) = e1;
            (*minLen) = (l1 <= l2)? l1: l2;
            return 1;
        }
    }

    return 0;
}


int detect_mul_bubble_end_with_bubbles(asg_t *g, uint32_t* begs, uint32_t occ, 
uint32_t* endNode, long long* minLen, buf_t* b)
{
    uint32_t e, flag, e_s;
    long long l, i, l_s;

    if(occ < 1) return 0;

    flag = detect_single_path_with_single_bubbles(g, begs[0], &e, &l, b, (uint32_t)-1);

    if(flag == TWO_INPUT || flag == MUL_INPUT)
    {
        e_s = e;
        l_s = l;
    }
    else
    {
        return 0;
    }
    


    for (i = 1; i < occ; i++)
    {
        flag = detect_single_path_with_single_bubbles(g, begs[i], &e, &l, b, (uint32_t)-1);
        if(flag == TWO_INPUT || flag == MUL_INPUT)
        {
            if(e != e_s) return 0;
            if(l < l_s) l_s = l;
        }
        else
        {
            return 0;
        }
    }

    if(asg_arc_n(g, e_s^1) == occ)
    {
        (*endNode) = e_s;
        (*minLen) = l_s;
        return 1;
    }
    
    return 0;
}


int detect_bubble_with_bubbles(asg_t *g, uint32_t begNode, uint32_t* endNode, long long* minLen, 
buf_t* b, uint32_t max_ext)
{
    uint32_t e1, e2;
    long long l1, l2;

    if(asg_arc_n(g, begNode) != 2)
    {
        return 0;
    }

    if(asg_is_single_edge(g, asg_arc_a(g, begNode)[0].v, begNode>>1)!=1
      || 
       asg_is_single_edge(g, asg_arc_a(g, begNode)[1].v, begNode>>1)!=1)
    {
        return 0;
    }

    if(b) kv_push(uint32_t, b->b, begNode>>1);


    if(detect_single_path_with_single_bubbles(g, asg_arc_a(g, begNode)[0].v, &e1, &l1, b, max_ext) == TWO_INPUT)
    {
        b->b.n--;
        if(detect_single_path_with_single_bubbles(g, asg_arc_a(g, begNode)[1].v, &e2, &l2, b, max_ext) == TWO_INPUT)
        {
            b->b.n--;
            if(e1 == e2)
            {
                (*endNode) = e1;
                (*minLen) = (l1 <= l2)? l1: l2;
                (*minLen)++;
                return 1;
            }
        }
    }

    return 0;
}

int test_triangular_exact(asg_t *g, uint32_t* nodes, uint32_t length, 
uint32_t startNode, uint32_t endNode, int max_dist, buf_t* bub)
{
    
    uint32_t i, v, w;
    ///int flag0, flag1, node;
    int n_reduced = 0, todel;
    long long NodeLen_first[3];
    long long NodeLen_second[3];

    uint32_t Ns_first[3];
    uint32_t Ns_second[3];
    
    
    for (i = 0; i < length; ++i)
    {
        v = nodes[i];

        if((v>>1) == (startNode>>1) || (v>>1) == (endNode>>1))
        {
            continue;
        }

        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if(nv != 2)
        {
            continue;
        }
        if(av[0].v == av[1].v)
        {
            continue;
        }
        /**********************test first node************************/
        NodeLen_first[0] = NodeLen_first[1] = NodeLen_first[2] = -1;
        if(asg_is_single_edge(g, av[0].v, v>>1) <= 2 && asg_is_single_edge(g, av[1].v, v>>1) <= 2)
        {
            NodeLen_first[asg_is_single_edge(g, av[0].v, v>>1)] = 0;
            NodeLen_first[asg_is_single_edge(g, av[1].v, v>>1)] = 1;
        }
        ///one node has one out-edge, another node has two out-edges
        if(NodeLen_first[1] == -1 || NodeLen_first[2] == -1)
        {
            continue;
        }
       /**********************test first node************************/

        ///if the potiential edge has already been removed 
        if(av[NodeLen_first[2]].del == 1)
        {
            continue;
        }

        /**********************test second node************************/
        w = av[NodeLen_first[2]].v^1;
        asg_arc_t *aw = asg_arc_a(g, w);
        uint32_t nw = asg_arc_n(g, w);
        if(nw != 2)
        {
            fprintf(stderr, "error\n");
        }
        NodeLen_second[0] = NodeLen_second[1] = NodeLen_second[2] = -1;
        if(asg_is_single_edge(g, aw[0].v, w>>1) <= 2 && asg_is_single_edge(g, aw[1].v, w>>1) <= 2)
        {
            NodeLen_second[asg_is_single_edge(g, aw[0].v, w>>1)] = 0;
            NodeLen_second[asg_is_single_edge(g, aw[1].v, w>>1)] = 1;
        }
        ///one node has one out-edge, another node has two out-edges
        if(NodeLen_second[1] == -1 || NodeLen_second[2] == -1)
        {
            continue;
        }

        /**********************test second node************************/

        if(if_node_exist(nodes, length, (w>>1)) && ((w>>1) != (endNode>>1)))
        {
            
                //if(av[NodeLen_first[1]].el == 1 && av[NodeLen_first[2]].el == 0 &&
                //aw[NodeLen_second[1]].el == 1 && aw[NodeLen_second[2]].el == 0)
               ///if(av[NodeLen_first[2]].el == 0 && aw[NodeLen_second[2]].el == 0)
                //    if(
                //    (av[NodeLen_first[1]].strong == 1 && av[NodeLen_first[2]].strong == 0 &&
                //    aw[NodeLen_second[1]].strong == 1 && aw[NodeLen_second[2]].strong == 0) 
                //    ||
                //    (av[NodeLen_first[1]].el == 1 && av[NodeLen_first[2]].el == 0 &&
                //    aw[NodeLen_second[1]].el == 1 && aw[NodeLen_second[2]].el == 0))
               uint32_t convex1, convex2, f1, f2;
               long long l1, l2;
               todel = 0;
               f1 = detect_bubble_end_with_bubbles(g, av[0].v, av[1].v, &convex1, &l1, NULL);
               f2 = detect_bubble_end_with_bubbles(g, aw[0].v, aw[1].v, &convex2, &l2, NULL);
               if(f1 && f2)
               {
                    if(((convex1>>1) == (startNode>>1) || (convex1>>1) == (endNode>>1)) && 
                     ((convex2>>1) == (startNode>>1) || (convex2>>1) == (endNode>>1)))
                     {
                         if(l1 <= min_thres || l2 <= min_thres)
                         {
                             continue;
                         } 

                         todel = 1;
                     } 
               }
               else if(f1)
               {
                  if(((convex1>>1) == (startNode>>1) || (convex1>>1) == (endNode>>1)))
                   {
                       if(l1 <= min_thres)
                       {
                           continue;
                       }

                       todel = 1;
                   }   
               }
               else if(f2)
               {
                   if(((convex2>>1) == (startNode>>1) || (convex2>>1) == (endNode>>1)))
                   {
                       if(l2 <= min_thres)
                       {
                           continue;
                       }

                       todel = 1;
                   }   
               }

            
               if(todel == 0)
               {
                   if(!f1)
                   {
                       Ns_first[0] = av[0].v; Ns_first[1] = av[1].v;
                       f1 = asg_bub_end_finder_with_del_advance(g, Ns_first, 2, max_dist, 
                       bub, 0, (u_int32_t)-1, &convex1);
                       l1 = min_thres + 10;
                   }

                   if(!f2)
                   {
                       Ns_second[0] = aw[0].v; Ns_second[1] = aw[1].v;
                       f2 = asg_bub_end_finder_with_del_advance(g, Ns_second, 2, max_dist, 
                       bub, 0, (u_int32_t)-1, &convex2);
                       l2 = min_thres + 10;
                   }




                    
                    if(f1 && f2)
                    {
                        if(((convex1>>1) == (startNode>>1) || (convex1>>1) == (endNode>>1)) && 
                        ((convex2>>1) == (startNode>>1) || (convex2>>1) == (endNode>>1)))
                        {
                            if(l1 <= min_thres || l2 <= min_thres)
                            {
                                continue;
                            } 

                            todel = 1;
                        } 
                    }
                    else if(f1)
                    {
                        if(((convex1>>1) == (startNode>>1) || (convex1>>1) == (endNode>>1)))
                        {
                            if(l1 <= min_thres)
                            {
                                continue;
                            }

                            todel = 1;
                        }   
                    }
                    else if(f2)
                    {
                        if(((convex2>>1) == (startNode>>1) || (convex2>>1) == (endNode>>1)))
                        {
                            if(l2 <= min_thres)
                            {
                                continue;
                            }

                            todel = 1;
                        }   
                    }
                    
               }
               
               
               if(todel)
               {
                   av[NodeLen_first[2]].del = 1;
                    ///remove the reverse direction
                    asg_arc_del(g, av[NodeLen_first[2]].v^1, av[NodeLen_first[2]].ul>>32^1, 1);
                    n_reduced++;
                    ///fprintf(stderr, "***rm v>>1: %u, w>>1: %u\n", v>>1, w>>1);
               }
 
        }

    }

    return n_reduced;
}





int find_single_link(asg_t *g, uint32_t link_beg, int linkLen, uint32_t* link_end)
{
    uint32_t v, w;
    int i;
    i = 0;
    v = link_beg^1;
    uint32_t nv, nw;
    asg_arc_t *av;
    int edgeLen = 0;
    int flag = -1;

    nv = asg_arc_n(g, v);
    av = asg_arc_a(g, v);
    if(nv != 1)
    {
        return 0;
    }
    v = av[0].v;

    while (edgeLen < linkLen)
    {
        nv = asg_arc_n(g, v);
        av = asg_arc_a(g, v);
        if(nv != 1)
        {
            return 0;
        }

        w = v^1;
        nw = asg_arc_n(g, w);
        if(nw == 2)
        {
            (*link_end) = w;
            return 1;
        }
        else if(nw > 2)
        {
            return 0;
        }
        
        v = av[0].v;
        edgeLen++;
    }


    return 0;    
}


int if_edge_exist(asg_arc_t* edges, uint32_t length, uint32_t query)
{
    uint32_t i;
    for (i = 0; i < length; ++i)
    {
        if((edges[i].v>>1) == query)
        {
            return 1;
        }
    }

    return 0;
}

int test_quadangular_with_addition_node(asg_t *g, uint32_t* nodes, uint32_t length, 
uint32_t addition_node_length)
{
    
    uint32_t i, v, w;
    int flag, occ_v_0, occ_v_1, occ_w_0, occ_w_1;
    int n_reduced = 0;
    uint32_t v_out2_node, w_out2_node;
    uint32_t cut_edge_v, cut_edge_w;
    for (i = 0; i < length; ++i)
    {
        v = nodes[i];
        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if(nv != 1)
        {
            continue;
        }
        /**********************test first node************************/
        flag = asg_is_single_edge(g, av[0].v, v>>1);
        if(flag != 2)
        {
            continue;
        }
       /**********************test first node************************/

        if(!find_single_link(g, v, addition_node_length, &w))
        {
            continue;
        }
        v = av[0].v^1;
        ///up to now, v and w is the node what we want
        asg_arc_t *aw = asg_arc_a(g, w);
        uint32_t nw = asg_arc_n(g, w);
        nv = asg_arc_n(g, v);
        av = asg_arc_a(g, v);
        if(nv!=2 || nw != 2)
        {
            fprintf(stderr, "error\n");
        }

        if(!if_node_exist(nodes, length, (v>>1)))
        {
            continue;
        }
        if(!if_node_exist(nodes, length, (w>>1)))
        {
            continue;
        }
        /**********************for v************************/
        occ_v_0 = asg_is_single_edge(g, av[0].v, v>>1);
        occ_v_1 = asg_is_single_edge(g, av[1].v, v>>1);
        if(occ_v_0 == occ_v_1)
        {
            continue;
        }
        if(occ_v_0 < 1 || occ_v_0 > 2)
        {
            continue;
        }
        if(occ_v_1 < 1 || occ_v_1 > 2)
        {
            continue;
        }

        if(occ_v_0 == 2)
        {
            v_out2_node = av[0].v^1;
            cut_edge_v = 0;
        }
        else
        {
            v_out2_node = av[1].v^1;
            cut_edge_v = 1;
        }
        if(!if_node_exist(nodes, length, (v_out2_node>>1)))
        {
            continue;
        }
        /**********************for v************************/

        /**********************for w************************/
        occ_w_0 = asg_is_single_edge(g, aw[0].v, w>>1);
        occ_w_1 = asg_is_single_edge(g, aw[1].v, w>>1);
        if(occ_w_0 == occ_w_1)
        {
            continue;
        }
        if(occ_w_0 < 1 || occ_w_0 > 2)
        {
            continue;
        }
        if(occ_w_1 < 1 || occ_w_1 > 2)
        {
            continue;
        }

        if(occ_w_0 == 2)
        {
            w_out2_node = aw[0].v^1;
            cut_edge_w = 0;
        }
        else
        {
            w_out2_node = aw[1].v^1;
            cut_edge_w = 1;
        }
        if(!if_node_exist(nodes, length, (w_out2_node>>1)))
        {
            continue;
        }
        /**********************for w************************/

        
        if(!if_edge_exist(asg_arc_a(g, w_out2_node), asg_arc_n(g, w_out2_node), (v_out2_node>>1)))
        {
            continue;
        }

        if(!if_edge_exist(asg_arc_a(g, v_out2_node), asg_arc_n(g, v_out2_node), (w_out2_node>>1)))
        {
            continue;
        }


        av[cut_edge_v].del = 1;
        ///remove the reverse direction
        asg_arc_del(g, av[cut_edge_v].v^1, av[cut_edge_v].ul>>32^1, 1);

        aw[cut_edge_w].del = 1;
        ///remove the reverse direction
        asg_arc_del(g, aw[cut_edge_w].v^1, aw[cut_edge_w].ul>>32^1, 1);


        n_reduced++;
    }

    return n_reduced;
}

int test_triangular_addition_exact(asg_t *g, uint32_t* nodes, uint32_t length, 
uint32_t startNode, uint32_t endNode, int max_dist, buf_t* bub)
{
    
    uint32_t i, j, v, w;
    ///int flag0, flag1, node;
    int n_reduced = 0, todel;
    uint32_t Nodes1[2];
    uint32_t Nodes2[2];

    uint32_t Ns_first[2];
    uint32_t Ns_second[2];

    /**
    if(startNode == 203)
    {
        uint32_t nv;
        asg_arc_t *av;

        fprintf(stderr, "start: %d, end:%d, length: %d\n", startNode>>1, endNode>>1, length);

        v = startNode;
        nv = asg_arc_n(g, v);
        av = asg_arc_a(g, v);
        fprintf(stderr, "******start: %d, nv: %d, dir: %d******\n", v>>1, nv, v&1);
        for (j = 0; j < nv; j++)
        {
            fprintf(stderr, "j: %d, w: %d, dir: %d\n", j, av[j].v>>1, av[j].v&1);
        }
        fprintf(stderr, "*********************************************\n");

        for (i = 0; i < length; ++i)
        {
            v = nodes[i];
            
            nv = asg_arc_n(g, v);
            av = asg_arc_a(g, v);
            fprintf(stderr, "******v: %d, nv: %d, dir: %d******\n", v>>1, nv, v&1);
            for (j = 0; j < nv; j++)
            {
                fprintf(stderr, "j: %d, w: %d, dir: %d\n", j, av[j].v>>1, av[j].v&1);
            }
            fprintf(stderr, "*********************************************\n");
        }


        v = endNode;
        nv = asg_arc_n(g, v);
        av = asg_arc_a(g, v);
        fprintf(stderr, "******end: %d, nv: %d, dir: %d******\n", v>>1, nv, v&1);
        for (j = 0; j < nv; j++)
        {
            fprintf(stderr, "j: %d, w: %d, dir: %d\n", j, av[j].v>>1, av[j].v&1);
        }
        fprintf(stderr, "*********************************************\n");
    }
    **/
    
    
    for (i = 0; i < length; ++i)
    {
        v = nodes[i];
        /**
        if(startNode == 203)
        {
            fprintf(stderr, "0 v: %d, i: %d\n", v, i);
        }
        **/

        if((v>>1) == (startNode>>1) || (v>>1) == (endNode>>1))
        {
            continue;
        }

        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if(nv != 1)
        {
            continue;
        }
        if(asg_is_single_edge(g, av[0].v, v>>1) != 2)
        {
            continue;
        }


        w = v^1;
        asg_arc_t *aw = asg_arc_a(g, w);
        uint32_t nw = asg_arc_n(g, w);
        if(nw != 1)
        {
            continue;
        }
        if(asg_is_single_edge(g, aw[0].v, w>>1) != 2)
        {
            continue;
        }

        if((av[0].v>>1) == (aw[0].v>>1))
        {
            continue;
        }



        Nodes1[0] = av[0].v^1;
        Nodes2[0] = aw[0].v^1;

        
        for(j = 0; j < 2; j++)
        {
            if((asg_arc_a(g, Nodes1[0])[j].v>>1)!= (v>>1))
            {
                Nodes1[1] = asg_arc_a(g, Nodes1[0])[j].v^1;
            }
        }

        for(j = 0; j < 2; j++)
        {
            if((asg_arc_a(g, Nodes2[0])[j].v>>1)!= (v>>1))
            {
                Nodes2[1] = asg_arc_a(g, Nodes2[0])[j].v^1;
            }
        }
        
        if(asg_arc_n(g, Nodes1[1]) != 1 || asg_arc_n(g, Nodes2[1]) != 1)
        {
            continue;
        }

        if((Nodes1[1]>>1) == (Nodes2[1]>>1))
        {
            continue;
        }



        if(asg_arc_a(g, Nodes1[1])[0].el == 0 || asg_arc_a(g, Nodes2[1])[0].el == 0)
        {
            continue;
        }

        

        /**
        if(startNode == 203)
        {
            fprintf(stderr, "1 v>>1: %d, i: %d\n", v>>1, i);
        }
        **/


        uint32_t convex1, convex2, f1, f2;
        long long l1, l2;
        todel = 0;

        if(Nodes1[0]^1 == startNode^1 || Nodes1[0]^1 == endNode)
        {
            continue;
        }
        if(Nodes2[1]^1 == startNode^1 || Nodes2[1]^1 == endNode)
        {
            continue;
        }

        f1 = detect_bubble_end_with_bubbles(g, Nodes1[0]^1, Nodes2[1]^1, &convex1, &l1, NULL);


        

        /**
        if(startNode == 203)
        {
            fprintf(stderr, "2 v: %d, i: %d\n", v, i);
        }
        **/

        if(Nodes2[0]^1 == startNode^1 || Nodes2[0]^1 == endNode)
        {
            continue;
        }
        if(Nodes1[1]^1 == startNode^1 || Nodes1[1]^1 == endNode)
        {
            continue;
        }

        f2 = detect_bubble_end_with_bubbles(g, Nodes2[0]^1, Nodes1[1]^1, &convex2, &l2, NULL);


        

        /**
        if(startNode == 203)
        {
            fprintf(stderr, "3 v: %d, i: %d\n", v, i);
        }
        **/

        if(f1 && f2)
        {
            if(((convex1>>1) == (startNode>>1) || (convex1>>1) == (endNode>>1)) && 
                ((convex2>>1) == (startNode>>1) || (convex2>>1) == (endNode>>1)))
            {
                if(l1 <= min_thres || l2 <= min_thres)
                {
                    continue;
                }

                todel = 1;
            }
        }
        else if(f1)
        {
           if(((convex1>>1) == (startNode>>1) || (convex1>>1) == (endNode>>1)))
            {
                if(l1 <= min_thres)
                {
                    continue;
                }

                todel = 1;
            }
        }
        else if(f2)
        {
            if(((convex2>>1) == (startNode>>1) || (convex2>>1) == (endNode>>1)))
            {
                if(l2 <= min_thres)
                {
                    continue;
                }

                todel = 1;
            }
        }

        /**
        if(startNode == 203)
        {
            fprintf(stderr, "4 v: %d, i: %d\n", v, i);
        }
        **/


        if(todel == 0)
        {
            if(!f1)
            {
                Ns_first[0] = Nodes1[0]^1; Ns_first[1] = Nodes2[1]^1;
                f1 = asg_bub_end_finder_with_del_advance(g, Ns_first, 2, max_dist, 
                bub, 0, (u_int32_t)-1, &convex1);
                l1 = min_thres + 10;
            }

            if(!f2)
            {
                Ns_second[0] = Nodes2[0]^1; Ns_second[1] = Nodes1[1]^1;
                f2 = asg_bub_end_finder_with_del_advance(g, Ns_second, 2, max_dist, 
                bub, 0, (u_int32_t)-1, &convex2);
                l2 = min_thres + 10;
            }

            if(f1 && f2)
            {
                if(((convex1>>1) == (startNode>>1) || (convex1>>1) == (endNode>>1)) && 
                ((convex2>>1) == (startNode>>1) || (convex2>>1) == (endNode>>1)))
                {
                    if(l1 <= min_thres || l2 <= min_thres)
                    {
                        continue;
                    } 

                    todel = 1;
                } 
            }
            else if(f1)
            {
                if(((convex1>>1) == (startNode>>1) || (convex1>>1) == (endNode>>1)))
                {
                    if(l1 <= min_thres)
                    {
                        continue;
                    }

                    todel = 1;
                }   
            }
            else if(f2)
            {
                if(((convex2>>1) == (startNode>>1) || (convex2>>1) == (endNode>>1)))
                {
                    if(l2 <= min_thres)
                    {
                        continue;
                    }

                    todel = 1;
                }   
            }
            
        }
























        
        if(todel)
        {
            if(av[0].el == 0 || aw[0].el == 0)
            {
                av[0].del = 1;
                asg_arc_del(g, av[0].v^1, av[0].ul>>32^1, 1);

                aw[0].del = 1;
                asg_arc_del(g, aw[0].v^1, aw[0].ul>>32^1, 1);
                
                n_reduced++;
            }
        }
    }

    return n_reduced;
}

int asg_arc_del_triangular_advance(asg_t *g, long long max_dist)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
	uint32_t v, w, n_vtx = g->n_seq * 2, n_reduced = 0, n_reduced_a = 0;


    if (!g->is_symm) asg_symm(g);


    buf_t b;
	memset(&b, 0, sizeof(buf_t));
	b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));


    buf_t bub;
    memset(&bub, 0, sizeof(buf_t));
    bub.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));



    int flag0, flag1, node;
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if (g->seq[v>>1].del)
        {
            continue;
        } 

        if(nv < 2)
        {
            continue;
        }

        // if(v == 17929187)
        // {
        //     fprintf(stderr, "*********v: %d\n", v);
        //     fflush(stderr);
        // }

        

        ///if this is a bubble
        if(asg_bub_finder_with_del_advance(g, v, max_dist, &b) == 1)
        {
            n_reduced += test_triangular_exact(g, b.b.a, b.b.n, v, b.S.a[0], max_dist, &bub);

            // if(v == 17929187)
            // {
            //     fprintf(stderr, "??????v: %d, b.S.a[0]: %d\n", v, b.S.a[0]);
            //     fflush(stderr);
            // }

            n_reduced_a += test_triangular_addition_exact(g, b.b.a, b.b.n, v, b.S.a[0],max_dist, &bub);
        }

        // if(v == 17929187)
        // {
        //     fprintf(stderr, "##########v: %d\n", v);
        //     fflush(stderr);
        // }
    
    }

    free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
    free(bub.a); free(bub.S.a); free(bub.T.a); free(bub.b.a); free(bub.e.a);

    if (n_reduced + n_reduced_a) {
		asg_cleanup(g);
		asg_symm(g);
	}

    fprintf(stderr, "[M::%s] removed %d/%d triangular/triangular_a overlaps\n", 
    __func__, n_reduced, n_reduced_a);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

    return n_reduced + n_reduced_a;
}




int asg_arc_del_triangular_advance_debug(asg_t *g, long long max_dist)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
	uint32_t v, w, n_vtx = g->n_seq * 2, n_reduced = 0, n_reduced_a = 0;


    if (!g->is_symm) asg_symm(g);


    buf_t b;
	memset(&b, 0, sizeof(buf_t));
	b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));


    buf_t bub;
    memset(&bub, 0, sizeof(buf_t));
    bub.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));

    fprintf(stderr, "n_vtx: %d\n", n_vtx);

    int flag0, flag1, node;
    for (v = 0; v < n_vtx; ++v) 
    {

        fprintf(stderr, "0 v: %d\n", v);
        fflush(stderr);

        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if (g->seq[v>>1].del)
        {
            continue;
        } 

        if(nv < 2)
        {
            continue;
        }

        // if(v == 17929187)
        // {
        //     fprintf(stderr, "*********v: %d\n", v);
        //     fflush(stderr);
        // }

        fprintf(stderr, "1 v: %d\n", v);
        fflush(stderr);

        ///if this is a bubble
        if(asg_bub_finder_with_del_advance(g, v, max_dist, &b) == 1)
        {
            fprintf(stderr, "2 v: %d\n", v);
            fflush(stderr);

            n_reduced += test_triangular_exact(g, b.b.a, b.b.n, v, b.S.a[0], max_dist, &bub);

            fprintf(stderr, "3 v: %d\n", v);
            fflush(stderr);

            if(v == 203)
            {
                fprintf(stderr, "??????v: %d, b.S.a[0]: %d, b.b.n: %d\n", v, b.S.a[0], b.b.n);
                fflush(stderr);
            }

            n_reduced_a += test_triangular_addition_exact(g, b.b.a, b.b.n, v, b.S.a[0],max_dist, &bub);
        
            fprintf(stderr, "5 v: %d\n", v);
            fflush(stderr);
        
        }

        fprintf(stderr, "5 v: %d\n", v);
        fflush(stderr);

        // if(v == 17929187)
        // {
        //     fprintf(stderr, "##########v: %d\n", v);
        //     fflush(stderr);
        // }
    
    }

    free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
    free(bub.a); free(bub.S.a); free(bub.T.a); free(bub.b.a); free(bub.e.a);

    if (n_reduced + n_reduced_a) {
		asg_cleanup(g);
		asg_symm(g);
	}

    fprintf(stderr, "[M::%s] removed %d/%d triangular/triangular_a overlaps\n", 
    __func__, n_reduced, n_reduced_a);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

    return n_reduced + n_reduced_a;
}


int check_if_cross(asg_t *g, uint32_t v)
{
    uint32_t N_list[5];
    if (g->seq[v>>1].del) return 0;
    uint32_t nv = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v);
    if(nv != 2) return 0;
    if(asg_is_single_edge(g, av[0].v, v>>1) != 2 || asg_is_single_edge(g, av[1].v, v>>1) != 2)
    {
        return 0;
    }
    if(av[0].v == av[1].v)
    {
        return 0;
    }
    N_list[0] = v;
    N_list[1] = av[0].v^1;
    N_list[2] = av[1].v^1;
    if(asg_arc_n(g, N_list[0]) != 2 || 
       asg_arc_n(g, N_list[1]) != 2 ||
       asg_arc_n(g, N_list[2]) != 2 )
    {
        return 0;
    }

    if(asg_arc_a(g, N_list[1])[0].v == asg_arc_a(g, N_list[1])[1].v)
    {
        return 0;
    }

    if(asg_arc_a(g, N_list[2])[0].v == asg_arc_a(g, N_list[2])[1].v)
    {
        return 0;
    }

    if(asg_arc_a(g, N_list[1])[0].v == (N_list[0]^1))
    {
        N_list[3] = asg_arc_a(g, N_list[1])[1].v^1;
    }
    else if(asg_arc_a(g, N_list[1])[1].v == (N_list[0]^1))
    {
        N_list[3] = asg_arc_a(g, N_list[1])[0].v^1;
    }

    if(asg_arc_a(g, N_list[2])[0].v == (N_list[0]^1))
    {
        N_list[4] = asg_arc_a(g, N_list[2])[1].v^1;
    }
    else if(asg_arc_a(g, N_list[2])[1].v == (N_list[0]^1))
    {
        N_list[4] = asg_arc_a(g, N_list[2])[0].v^1;
    }

    if(N_list[3] != N_list[4])
    {
        return 0;
    }

    if(asg_arc_n(g, N_list[0]) != 2 || 
        asg_arc_n(g, N_list[1]) != 2 ||
        asg_arc_n(g, N_list[2]) != 2 ||
        asg_arc_n(g, N_list[3]) != 2)
    {
        return 0;
    }


    uint32_t convex1, convex2, f1, f2;
    long long l1, l2;
    l1 = l2 = 0;
    int todel = 0;
    
    f1 = detect_bubble_end_with_bubbles(g, N_list[0]^1, N_list[3]^1, &convex1, &l1, NULL);
    f2 = detect_bubble_end_with_bubbles(g, N_list[1]^1, N_list[2]^1, &convex2, &l2, NULL);

    if(f1 && f2)
    {
        if(l1 > min_thres && l2 > min_thres)
        {
            todel = 1;
        }     
    }
    else if(f1)
    {
        if(l1 > min_thres)
        {
            todel = 1;
        }        
    }
    else if(f2)
    {
        if(l2 > min_thres)
        {
            todel = 1;
        }  
    }
    

    return todel;
}

int asg_arc_identify_simple_bubbles(asg_t *g)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
	uint32_t v, w, n_vtx = g->n_seq * 2, n_reduced = 0;
    buf_t b;
    memset(&b, 0, sizeof(buf_t));
    memset(g->seq_vis, 0, g->n_seq*2*sizeof(uint8_t));
    long long l, i;
    for (v = 0; v < n_vtx; ++v) 
    {
        if (g->seq[v>>1].del) continue;

        b.b.n = 0;

        if(g->seq_vis[v] != 1)
        {
            ///if(detect_bubble_with_bubbles(g, v, &w, &l, &b, (uint32_t)-1))
            if(detect_bubble_with_bubbles(g, v, &w, &l, &b, SMALL_BUBBLE_SIZE))
            {
                for (i = 0; i < b.b.n; i++)
                {
                    if(b.b.a[i] != (v>>1) && b.b.a[i] != (w>>1))
                    {
                        g->seq_vis[b.b.a[i]<<1] = 1;
                        g->seq_vis[(b.b.a[i]<<1) + 1] = 1;
                    }
                }
                g->seq_vis[v] = 1;
                g->seq_vis[w^1] = 1;

                
                // if(asg_arc_n(g, v) != 2 || asg_arc_n(g, w^1) != 2)
                // {
                //     fprintf(stderr, "error\n");
                // }
                // if(v>>1 != b.b.a[0] || w>>1 != b.b.a[b.b.n-1])
                // {
                //     fprintf(stderr, "sbsbsbs\n");
                // }
                
            }
        }
        

        if(check_if_cross(g, v))
        {
            g->seq_vis[v] = 2;
        }
    }
    free(b.b.a);

    long long nodes, bub_nodes, cross_nodes;
    bub_nodes = nodes = cross_nodes = 0;
    for (v = 0; v < n_vtx; ++v) 
    {
        if (g->seq[v>>1].del) continue;
        nodes++;
        if(g->seq_vis[v] == 1) bub_nodes++;
        if(g->seq_vis[v] == 2) cross_nodes++;
    }


    fprintf(stderr, "[M::%s]  nodes:%d, bub_nodes: %d, cross_nodes: %d\n", __func__, 
    nodes, bub_nodes, cross_nodes);
    

    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

    return n_reduced;
}


typedef struct {
	int threadID;
    int thread_num;
    int check_cross;
    asg_t *g;
} para_for_simple_bub;

void* asg_arc_identify_simple_bubbles_pthread(void* arg)
{
    int thr_ID = ((para_for_simple_bub*)arg)->threadID;
    int thr_num = ((para_for_simple_bub*)arg)->thread_num;
    asg_t *g = ((para_for_simple_bub*)arg)->g;
    int check_cross = ((para_for_simple_bub*)arg)->check_cross;
    ///the reason is that each read has two direction (query->target, target->query)    

	uint32_t v, w, n_vtx = g->n_seq * 2, n_reduced = 0;
    buf_t b;
    memset(&b, 0, sizeof(buf_t));
    long long l, i;
    ///for (v = 0; v < n_vtx; ++v)
    for (v = thr_ID; v < n_vtx; v = v + thr_num) 
    {
        if (g->seq[v>>1].del) continue;

        b.b.n = 0;

        if(g->seq_vis[v] != 1)
        {
            ///if(detect_bubble_with_bubbles(g, v, &w, &l, &b, (uint32_t)-1))
            if(detect_bubble_with_bubbles(g, v, &w, &l, &b, SMALL_BUBBLE_SIZE))
            {
                for (i = 0; i < b.b.n; i++)
                {
                    if(b.b.a[i] != (v>>1) && b.b.a[i] != (w>>1))
                    {
                        g->seq_vis[b.b.a[i]<<1] = 1;
                        g->seq_vis[(b.b.a[i]<<1) + 1] = 1;
                    }
                }
                g->seq_vis[v] = 1;
                g->seq_vis[w^1] = 1;
            }
        }


        ///if(v%10000 == 0)
        //if(v >= 18020000)
        // if(v == 18021291)
        // {
        //     fprintf(stderr, "1 v: %d, thr_ID: %d\n", v, thr_ID);
        //     fflush(stderr);
        // }
        

        if(check_cross == 1 && check_if_cross(g, v))
        {
            g->seq_vis[v] = 2;
        }

        ///if(v%10000 == 0)
        //if(v >= 18020000)
        // if(v == 18021291)
        // {
        //     fprintf(stderr, "2 v: %d, thr_ID: %d\n", v, thr_ID);
        //     fflush(stderr);
        // }
    }
    free(b.b.a);

    free(arg);

    // fprintf(stderr, "thr_ID: %d end\n", thr_ID);
    // fflush(stderr);
}

int asg_arc_identify_simple_bubbles_multi(asg_t *g, int check_cross)
{
    double startTime = Get_T();
    memset(g->seq_vis, 0, g->n_seq*2*sizeof(uint8_t));

    pthread_t *_r_threads;

	_r_threads = (pthread_t *)malloc(sizeof(pthread_t)*thread_num);

    int i = 0;

	for (i = 0; i < thread_num; i++)
	{
        para_for_simple_bub* arg = (para_for_simple_bub*)malloc(sizeof(*arg));
        arg->g = g;
        arg->thread_num = thread_num;
        arg->threadID = i;
        arg->check_cross = check_cross;

        pthread_create(_r_threads + i, NULL, asg_arc_identify_simple_bubbles_pthread, (void*)arg);
	}
    

    for (i = 0; i<thread_num; i++)
		pthread_join(_r_threads[i], NULL);


    free(_r_threads);

    uint32_t v, n_vtx = g->n_seq * 2;
    long long nodes, bub_nodes, cross_nodes;
    bub_nodes = nodes = cross_nodes = 0;
    for (v = 0; v < n_vtx; ++v) 
    {
        if (g->seq[v>>1].del) continue;
        nodes++;
        if(g->seq_vis[v] == 1) bub_nodes++;
        if(g->seq_vis[v] == 2) cross_nodes++;
    }


    fprintf(stderr, "[M::%s]  nodes:%d, bub_nodes: %d, cross_nodes: %d\n", __func__, 
    nodes, bub_nodes, cross_nodes);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

    return bub_nodes+cross_nodes;
}



int check_small_bubble(asg_t *g, uint32_t begNode, uint32_t v, uint32_t w, 
long long* vLen, long long* wLen, uint32_t* endNode)
{
    uint32_t nv = asg_arc_n(g, v);
    uint32_t nw = asg_arc_n(g, w);

    asg_arc_t *av = asg_arc_a(g, v);
    asg_arc_t *aw = asg_arc_a(g, w);
    if(nv != 1 || nw != 1)
    {
        return 0;
    }

    ///first node
    ///nv must be 1
    if(asg_is_single_edge(g, av[0].v, v>>1) == 2)
    {
        uint32_t vv;
        vv = av[0].v^1;

        if(
            asg_is_single_edge(g, asg_arc_a(g, vv)[0].v, vv>>1) == 1
            &&
            asg_is_single_edge(g, asg_arc_a(g, vv)[1].v, vv>>1) == 1
            )
        {
            ///walk along first path
            long long pLen1;
            pLen1 = single_edge_length(g, asg_arc_a(g, vv)[0].v, begNode>>1, 1000);

            ///walk along first path
            long long pLen2;
            pLen2 = single_edge_length(g, asg_arc_a(g, vv)[1].v, begNode>>1, 1000);



            if(pLen1 >= 0 && pLen2 >= 0)
            {
                if(((asg_arc_a(g, vv)[0].v) == (v^1)) && pLen1 == 1)
                {
                    (*vLen) = pLen1;
                    (*wLen) = pLen2;
                }
                else if(((asg_arc_a(g, vv)[1].v) == (v^1)) && pLen2 == 1)
                {
                    (*vLen) = pLen2;
                    (*wLen) = pLen1;
                }
                else
                {
                    fprintf(stderr, "error\n");
                }
                ///(*endNode) = vv>>1;
                (*endNode) = vv;
                return 1;
            }

            /**
            if(pLen1 != 1 && pLen2 != 1)
            {
                fprintf(stderr, "v: pLen1: %d, pLen2: %d\n", pLen1, pLen2);
                pLen1 = single_edge_length(g, v^1, begNode>>1);
                fprintf(stderr, "pLen_v: %d\n", pLen1);
                fprintf(stderr, "vv>>1: %u, v>>1: %u, begNode>>1: %u\n", 
                vv>>1, v>>1, begNode>>1);
            }

            if((asg_arc_a(g, vv)[0].v) != (v^1) && (asg_arc_a(g, vv)[1].v) != (v^1))
            {
                fprintf(stderr, "vv: %u\n", vv);
            }
            **/

        }
    }


    ///second node
    ///nw must be 1
    if(asg_is_single_edge(g, aw[0].v, w>>1) == 2)
    {
        uint32_t ww;
        ww = aw[0].v^1;

        if(
            asg_is_single_edge(g, asg_arc_a(g, ww)[0].v, ww>>1) == 1
            &&
            asg_is_single_edge(g, asg_arc_a(g, ww)[1].v, ww>>1) == 1
            )
        {
            ///walk along first path
            long long pLen1;
            pLen1 = single_edge_length(g, asg_arc_a(g, ww)[0].v, begNode>>1, 1000);

            ///walk along first path
            long long pLen2;
            pLen2 = single_edge_length(g, asg_arc_a(g, ww)[1].v, begNode>>1, 1000);

            if(pLen1 >= 0 && pLen2 >= 0)
            {
                if(((asg_arc_a(g, ww)[0].v) == (w^1)) && pLen1 == 1)
                {
                    (*wLen) = pLen1;
                    (*vLen) = pLen2;
                }
                else if(((asg_arc_a(g, ww)[1].v) == (w^1)) && pLen2 == 1)
                {
                    (*wLen) = pLen2;
                    (*vLen) = pLen1;
                }
                else
                {
                    fprintf(stderr, "error\n");
                }

                //(*endNode) = ww>>1;
                (*endNode) = ww;

                return 1;
            }

            /**
            if(pLen1 != 1 && pLen2 != 1)
            {
                fprintf(stderr, "w: pLen1: %d, pLen2: %d\n", pLen1, pLen2);
                pLen1 = single_edge_length(g, w^1, begNode>>1);
                fprintf(stderr, "pLen_w: %d\n", pLen1);
            }

            if((asg_arc_a(g, ww)[0].v) != (w^1) && (asg_arc_a(g, ww)[1].v) != (w^1))
            {
                fprintf(stderr, "ww: %u\n", ww);
            }
            **/
       }
    }

    return 0;

}


int test_single_node_bubble(asg_t *g, uint32_t* nodes, uint32_t length, 
uint32_t startNode, uint32_t endNode)
{
    
    uint32_t i, v, w;
    uint32_t vEnd;
    int flag0, flag1, node;
    int n_reduced = 0;
    long long Len[2], longLen;
    long long longLen_thres = 4;
    for (i = 0; i < length; ++i)
    {
        v = nodes[i];

        if((v>>1) == (startNode>>1) || (v>>1) == (endNode>>1))
        {
            continue;
        }

        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if(nv != 2)
        {
            continue;
        }

        
        flag0 = asg_is_single_edge(g, av[0].v, v>>1);
        flag1 = asg_is_single_edge(g, av[1].v, v>>1);
        if(flag0 != 1 || flag1 != 1)
        {
            continue;
        }

        if(check_small_bubble(g, v, av[0].v, av[1].v, &(Len[0]), &(Len[1]), &vEnd))
        {

            if(if_node_exist(nodes, length, vEnd>>1) && ((vEnd>>1) != (endNode>>1)))
           {
               
               if(Len[0] == 1 && Len[1] != 1)
               {
                   w = av[0].v;
                   longLen = Len[1];
                    
                  /****************************may have bugs********************************/
                  if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0)
                   {
                    //    fprintf(stderr, "w>>1: %u, beg: %u, end: %u\n", 
                    //    w>>1, startNode>>1, endNode>>1);
                       asg_seq_del(g, w>>1);
                       n_reduced++;
                   }///up to here w is exactly overlapped in both directions
                   else if(longLen >= longLen_thres)
                   {
                       if(av[0].el == 1 && av[1].el == 1
                          &&
                          asg_arc_a(g, vEnd)[0].el == 1 && asg_arc_a(g, vEnd)[1].el == 1)
                       {
                           asg_seq_del(g, w>>1);
                           n_reduced++;
                       }
                   }
                   
                  /****************************may have bugs********************************/
                   
               }
               else if(Len[0] != 1 && Len[1] == 1)
               {
                   w = av[1].v;
                   longLen = Len[0];
                   
                  /****************************may have bugs********************************/
                  if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0)
                   {
                    //    fprintf(stderr, "w>>1: %u, beg: %u, end: %u\n", 
                    //    w>>1, startNode>>1, endNode>>1);
                       asg_seq_del(g, w>>1);
                       n_reduced++;
                   }///up to here w is exactly overlapped in both directions
                   else if(longLen >= longLen_thres)
                   {
                       if(av[0].el == 1 && av[1].el == 1
                          &&
                          asg_arc_a(g, vEnd)[0].el == 1 && asg_arc_a(g, vEnd)[1].el == 1)
                       {
                           asg_seq_del(g, w>>1);
                           n_reduced++;
                       }
                   }
                  /****************************may have bugs********************************/
               }
               else if(Len[0] == 1 && Len[1] == 1)
               {
                   w = av[0].v;
                   flag0 = asg_arc_a(g, w)[0].el + asg_arc_a(g, w^1)[0].el;
                   w = av[1].v;
                   flag1 = asg_arc_a(g, w)[0].el + asg_arc_a(g, w^1)[0].el;
                    ///>=2 means this is an exact overlap
                   if(flag0 < 2 && flag1 >= 2)
                   {
                       w = av[0].v;
                       
                       /****************************may have bugs********************************/
                        if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0)
                        {
                            // fprintf(stderr, "w>>1: %u, beg: %u, end: %u\n", 
                            // w>>1, startNode>>1, endNode>>1);
                            asg_seq_del(g, w>>1);
                            n_reduced++;
                        }
                        /****************************may have bugs********************************/
                   }

                   if(flag0 >= 2 && flag1 < 2)
                   {
                       w = av[1].v;
                      
                       /****************************may have bugs********************************/
                        if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0)
                        {
                            // fprintf(stderr, "w>>1: %u, beg: %u, end: %u\n", 
                            // w>>1, startNode>>1, endNode>>1);
                            asg_seq_del(g, w>>1);
                            n_reduced++;
                        }
                        /****************************may have bugs********************************/
                   }
               }
               else
               {
                   fprintf(stderr, "error\n");
               }
               
           }
        }

    }

    return n_reduced;
}


int test_single_node_bubble_directly(asg_t *g, uint32_t v, long long longLen_thres, ma_hit_t_alloc* sources)
{
    
    uint32_t w, vEnd;
    int flag0, flag1, node;
    int n_reduced = 0;
    long long Len[2], longLen;
    ///long long longLen_thres = 4;


    

    uint32_t nv = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v);
    if(nv != 2)
    {
        return 0;
    }

    
    flag0 = asg_is_single_edge(g, av[0].v, v>>1);
    flag1 = asg_is_single_edge(g, av[1].v, v>>1);
    if(flag0 != 1 || flag1 != 1)
    {
        return 0;
    }

    if(check_small_bubble(g, v, av[0].v, av[1].v, &(Len[0]), &(Len[1]), &vEnd))
    {
        if(Len[0] == 1 && Len[1] != 1)
        {
            w = av[0].v;
            longLen = Len[1];
            
            /****************************may have bugs********************************/
            ///if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0)
            if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0 || sources[w>>1].is_abnormal == 1)
            {
                asg_seq_del(g, w>>1);
                n_reduced++;
            }///up to here w is exactly overlapped in both directions
            else if(longLen >= longLen_thres)
            {
                if(av[0].el == 1 && av[1].el == 1
                    &&
                    asg_arc_a(g, vEnd)[0].el == 1 && asg_arc_a(g, vEnd)[1].el == 1)
                {
                    asg_seq_del(g, w>>1);
                    n_reduced++;
                }
            }
            
            /****************************may have bugs********************************/
            
        }
        else if(Len[0] != 1 && Len[1] == 1)
        {
            w = av[1].v;
            longLen = Len[0];
            
            /****************************may have bugs********************************/
            ///if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0)
            if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0 || sources[w>>1].is_abnormal == 1)
            {
                asg_seq_del(g, w>>1);
                n_reduced++;
            }///up to here w is exactly overlapped in both directions
            else if(longLen >= longLen_thres)
            {
                if(av[0].el == 1 && av[1].el == 1
                    &&
                    asg_arc_a(g, vEnd)[0].el == 1 && asg_arc_a(g, vEnd)[1].el == 1)
                {
                    asg_seq_del(g, w>>1);
                    n_reduced++;
                }
            }
            /****************************may have bugs********************************/
        }
        else if(Len[0] == 1 && Len[1] == 1)
        {
            flag0 = sources[av[0].v>>1].is_abnormal;
            flag1 = sources[av[1].v>>1].is_abnormal;

            if(flag0 == 1 && flag1 == 0)
            {
                asg_seq_del(g, av[0].v>>1);
                n_reduced++;
            }
            else if(flag0 == 0 && flag1 == 1)
            {
                asg_seq_del(g, av[1].v>>1);
                n_reduced++;
            }
            else
            {
                w = av[0].v;
                flag0 = asg_arc_a(g, w)[0].el + asg_arc_a(g, w^1)[0].el;
                w = av[1].v;
                flag1 = asg_arc_a(g, w)[0].el + asg_arc_a(g, w^1)[0].el;
                ///>=2 means this is an exact overlap
                if(flag0 < 2 && flag1 >= 2)
                {
                    w = av[0].v;
                    
                    /****************************may have bugs********************************/
                    if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0)
                    {
                        asg_seq_del(g, w>>1);
                        n_reduced++;
                    }
                    /****************************may have bugs********************************/
                }

                if(flag0 >= 2 && flag1 < 2)
                {
                    w = av[1].v;
                    
                    /****************************may have bugs********************************/
                    if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0)
                    {
                        asg_seq_del(g, w>>1);
                        n_reduced++;
                    }
                    /****************************may have bugs********************************/
                }
            }
        }
        
    }
    return n_reduced;
}


int asg_arc_del_single_node_bubble(asg_t *g, long long max_dist)
{
    ///the reason is that each read has two direction (query->target, target->query)
	uint32_t v, w, n_vtx = g->n_seq * 2, n_reduced = 0;
    buf_t b;
	if (!g->is_symm) asg_symm(g);
	memset(&b, 0, sizeof(buf_t));
	///set information for each node
	b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    int flag0, flag1, node;
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if (g->seq[v>>1].del)
        {
            continue;
        } 

        if(nv < 2)
        {
            continue;
        }

        ///if this is a bubble
        if(asg_bub_finder_with_del_advance(g, v, max_dist, &b) == 1)
        {
            n_reduced += test_single_node_bubble(g, b.b.a, b.b.n, v, b.S.a[0]);
        }
    
    }

    free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);

    if (n_reduced) {
		asg_cleanup(g);
		asg_symm(g);
	}

    fprintf(stderr, "[M::%s] removed %d short bubbles\n\n", __func__, n_reduced);

    return n_reduced;
}

int asg_arc_del_single_node_directly(asg_t *g, long long longLen_thres, ma_hit_t_alloc* sources)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
	uint32_t v, w, n_vtx = g->n_seq * 2, n_reduced = 0;
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if (g->seq[v>>1].del)
        {
            continue;
        } 

        if(nv != 2)
        {
            continue;
        }

        n_reduced += test_single_node_bubble_directly(g, v, longLen_thres, sources);
    }


    if (n_reduced) {
		asg_cleanup(g);
		asg_symm(g);
	}

    fprintf(stderr, "[M::%s] removed %d small bubbles\n", __func__, n_reduced);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

    return n_reduced;
}


int asg_arc_del_self_circle_contig(asg_t *g)
{
    double startTime = Get_T();
    uint32_t v;
    uint32_t n_vtx = g->n_seq * 2, n_reduced = 0;
    uint32_t vEnd;
    int flag0, flag1, node;
    long long Len[3];
    for (v = 0; v < n_vtx; ++v) 
    {
        if (g->seq[v>>1].del)
        {
            continue;
        } 

        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if(nv != 2)
        {
            continue;
        }
        if(av[0].v == av[1].v)
        {
            continue;
        }
        Len[0] = Len[1] = Len[2] = -1;
        if(asg_is_single_edge(g, av[0].v, v>>1) <= 2 
        && asg_is_single_edge(g, av[1].v, v>>1) <= 2)
        {
            Len[asg_is_single_edge(g, av[0].v, v>>1)] = 0;
            Len[asg_is_single_edge(g, av[1].v, v>>1)] = 1;
        }

        if(Len[1] == -1 || Len[2] == -1)
        {
            continue;
        }


        if(asg_arc_n(g, av[Len[2]].v) == 1 && 
        single_edge_length(g, av[Len[2]].v, v>>1, 100)!=-1)
        {
            av[Len[2]].del = 1;
            ///remove the reverse direction
            asg_arc_del(g, av[Len[2]].v^1, av[Len[2]].ul>>32^1, 1);
            n_reduced++;
        }
        
    }

    if (n_reduced) {
		asg_cleanup(g);
		asg_symm(g);
	}

    fprintf(stderr, "[M::%s] removed %d self-circle contig\n", __func__, n_reduced);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

    return n_reduced;
}


int test_cross(asg_t *g, uint32_t* nodes, uint32_t length, 
uint32_t startNode, uint32_t endNode)
{
    uint32_t a1, a2;
    uint32_t N_list[5];
    uint32_t i, v;
    int flag0, flag1;
    int n_reduced = 0;
    for (i = 0; i < length; ++i)
    {
        v = nodes[i];

        if((v>>1) == (startNode>>1) || (v>>1) == (endNode>>1))
        {
            continue;
        }

        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if(nv != 2)
        {
            continue;
        }
        if(av[0].v == av[1].v)
        {
            continue;
        }
        flag0 = asg_is_single_edge(g, av[0].v, v>>1);
        flag1 = asg_is_single_edge(g, av[1].v, v>>1);

        if(flag0 != 2 || flag1 != 2)
        {
            continue;
        }
        

        N_list[0] = v;
        N_list[1] = av[0].v^1;
        N_list[2] = av[1].v^1;

        if(asg_arc_n(g, N_list[0]) != 2 || 
           asg_arc_n(g, N_list[1]) != 2 ||
           asg_arc_n(g, N_list[2]) != 2 )
        {
            continue;
        }

        if(asg_arc_a(g, N_list[1])[0].v == asg_arc_a(g, N_list[1])[1].v)
        {
            continue;
        }

        if(asg_arc_a(g, N_list[2])[0].v == asg_arc_a(g, N_list[2])[1].v)
        {
            continue;
        }


        if(asg_arc_a(g, N_list[1])[0].v == (N_list[0]^1))
        {
            N_list[3] = asg_arc_a(g, N_list[1])[1].v^1;
        }
        else if(asg_arc_a(g, N_list[1])[1].v == (N_list[0]^1))
        {
            N_list[3] = asg_arc_a(g, N_list[1])[0].v^1;
        }
        else
        {
            fprintf(stderr, "ERROR\n");
        }
        
        if(asg_arc_a(g, N_list[2])[0].v == (N_list[0]^1))
        {
            N_list[4] = asg_arc_a(g, N_list[2])[1].v^1;
        }
        else if(asg_arc_a(g, N_list[2])[1].v == (N_list[0]^1))
        {
            N_list[4] = asg_arc_a(g, N_list[2])[0].v^1;
        }
        else
        {
            fprintf(stderr, "ERROR\n");
        }

        if(N_list[3] != N_list[4])
        {
            continue;
        }

        if(asg_arc_n(g, N_list[0]) != 2 || 
           asg_arc_n(g, N_list[1]) != 2 ||
           asg_arc_n(g, N_list[2]) != 2 ||
           asg_arc_n(g, N_list[3]) != 2)
        {
            continue;
        }
        /**
        N_list[3]        N_list[0]

        N_list[2]        N_list[1]
        **/
       if(asg_arc_a(g, N_list[0])[0].el == asg_arc_a(g, N_list[0])[1].el)
       {
           continue;
       }

       if(asg_arc_a(g, N_list[0])[0].el == 1)
       {
           //a1 = asg_arc_a(g, N_list[0])[0].v >> 1;
           a1 = 0;
       }
       else
       {
           ///a1 = asg_arc_a(g, N_list[0])[1].v >> 1;
           a1 = 1;
       }





       
       if(asg_arc_a(g, N_list[3])[0].el == asg_arc_a(g, N_list[3])[1].el)
       {
           continue;
       }

       if(asg_arc_a(g, N_list[3])[0].el == 1)
       {
           //a2 = asg_arc_a(g, N_list[3])[0].v >> 1;
           a2 = 0;
       }
       else
       {
           //a2 = asg_arc_a(g, N_list[3])[1].v >> 1;
           a2 = 1;
       }

        if(
            (asg_arc_a(g, N_list[0])[a1].v >> 1)
            !=
            (asg_arc_a(g, N_list[3])[a2].v >> 1)
           )
        {
            if(((N_list[0]>>1) != (endNode>>1)) &&
               ((N_list[1]>>1) != (endNode>>1)) && 
               ((N_list[2]>>1) != (endNode>>1)) &&
               ((N_list[3]>>1) != (endNode>>1)))
            {
                asg_arc_a(g, N_list[0])[a1].del = 1;
                asg_arc_del(g, asg_arc_a(g, N_list[0])[a1].v^1, 
                asg_arc_a(g, N_list[0])[a1].ul>>32^1, 1);


                asg_arc_a(g, N_list[3])[a2].del = 1;
                asg_arc_del(g, asg_arc_a(g, N_list[3])[a2].v^1, 
                asg_arc_a(g, N_list[3])[a2].ul>>32^1, 1);
                /**
                fprintf(stderr, "(N_list[0]>>1): %u, (N_list[1]>>1): %u, (N_list[2]>>1): %u, (N_list[3]>>1): %u\n", 
                (N_list[0]>>1), (N_list[1]>>1), (N_list[2]>>1), (N_list[3]>>1));
                fprintf(stderr, "a1: %u\n", asg_arc_a(g, N_list[0])[a1].v>>1);
                fprintf(stderr, "a2: %u\n", asg_arc_a(g, N_list[3])[a2].v>>1);
                **/

                n_reduced++;
            }
        }
    }

    return n_reduced;
}

int asg_arc_del_cross_bubble(asg_t *g, long long max_dist)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
	uint32_t v, w, n_vtx = g->n_seq * 2, n_reduced = 0;
    buf_t b;
	if (!g->is_symm) asg_symm(g);
	memset(&b, 0, sizeof(buf_t));
	///set information for each node
	b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    int flag0, flag1, node;
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if (g->seq[v>>1].del)
        {
            continue;
        } 

        if(nv < 2)
        {
            continue;
        }

        ///if this is a bubble
        if(asg_bub_finder_with_del_advance(g, v, max_dist, &b) == 1)
        {
            n_reduced += test_cross(g, b.b.a, b.b.n, v, b.S.a[0]);
        }
    
    }

    free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);

    if (n_reduced) {        
		asg_cleanup(g);
		asg_symm(g);
	}

    fprintf(stderr, "[M::%s] removed %d cross\n", __func__, n_reduced);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

    return n_reduced;
}



// transitive reduction; see Myers, 2005
int asg_arc_del_trans(asg_t *g, int fuzz)
{
    double startTime = Get_T();

	uint8_t *mark;
	///n_vtx = number of seq * 2
	///the reason is that each read has two direction (query->target, target->query)
	uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0;
	///at first, all nodes should be set to vacant
	mark = (uint8_t*)calloc(n_vtx, 1);

	/**v is the id+direction of a node, 
     * the high 32-bit is the id, 
     * and the lowest 1-bit is the direction
     * (0 means query-to-target, 1 means target-to-query)**/
	for (v = 0; v < n_vtx; ++v) {
		///nv is the number of overlaps with v(qn+direction)
		uint32_t L, i, nv = asg_arc_n(g, v);
		///av is the array of v
		asg_arc_t *av = asg_arc_a(g, v);
        ///that means in this direction, read v is not overlapped with any other reads
		if (nv == 0) continue; // no hits
		
        ///if the read itself has been removed
		if (g->seq[v>>1].del) 
        {
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


		//all outnode of v should be set to "not reduce"
		for (i = 0; i < nv; ++i) mark[av[i].v] = 1;

		///length of node (not overlap length)
		///av[nv-1] is longest out-dege
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
			//w is an out-node of v
			uint32_t w = av[i].v;
			
			uint32_t j, nw = asg_arc_n(g, w);
			asg_arc_t *aw = asg_arc_a(g, w);
			///if w has already been reduced
			if (mark[av[i].v] != 1) continue;

			for (j = 0; j < nw && asg_arc_len(aw[j]) + asg_arc_len(av[i]) <= L; ++j)
				if (mark[aw[j].v]) mark[aw[j].v] = 2;
		}
		#if 0
		for (i = 0; i < nv; ++i) {
			uint32_t w = av[i].v;
			uint32_t j, nw = asg_arc_n(g, w);
			asg_arc_t *aw = asg_arc_a(g, w);
			for (j = 0; j < nw && (j == 0 || asg_arc_len(aw[j]) < fuzz); ++j)
				if (mark[aw[j].v]) mark[aw[j].v] = 2;
		}
		#endif
        //remove edges
		for (i = 0; i < nv; ++i) {
			if (mark[av[i].v] == 2) av[i].del = 1, ++n_reduced;
			mark[av[i].v] = 0;
		}
	}
	free(mark);
	fprintf(stderr, "[M::%s] transitively reduced %d arcs\n", __func__, n_reduced);
	if (n_reduced) {
		asg_cleanup(g);
		asg_symm(g);
	}


    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
	return n_reduced;
}





///max_ext is 4
int asg_cut_tip(asg_t *g, int max_ext)
{
    double startTime = Get_T();

	asg64_v a = {0,0,0};
	uint32_t n_vtx = g->n_seq * 2, v, i, cnt = 0;
    
	for (v = 0; v < n_vtx; ++v) {
		//if this seq has been deleted
		if (g->seq[v>>1].del) continue;
		///check if the another direction of v has no overlaps
		///if the self direction of v has no overlaps, we don't have the overlaps of them
		///here is check if the reverse direction of v 
		/**
		 the following first line is to find (means v is a node has no prefix):
		        (v)--->()---->()---->()----->....
         another case is:
            ......()---->()---->()----->()------>(v)
         this case can be found by (v^1), so we don't need to process this case here
		**/
		if (asg_is_utg_end(g, v, 0) != ASG_ET_TIP) continue; // not a tip
        /**
         the following second line is: 
		 (v)--->()---->()---->()----->()
		        |--------max_ext-------|
        **/
       ///that means here is a long tip, which is longer than max_ext
		if (asg_extend(g, v, max_ext, &a) == ASG_ET_MERGEABLE) continue; // not a short unitig

		/**
		 * so combining the last two lines, they are designed to reomve(n(0), n(1), n(2)):
		 *                  		    ----->n(4)
		 *                  		   |
		 * n(0)--->n(1)---->n(2)---->n(3)
		 *                             |
		 *                              ----->n(5)
		**/
		for (i = 0; i < a.n; ++i)
			asg_seq_del(g, (uint32_t)a.a[i]>>1);
		++cnt;
	}
	free(a.a);
	if (cnt > 0) asg_cleanup(g);
	fprintf(stderr, "[M::%s] cut %d tips\n", __func__, cnt);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);


	return cnt;
}


///max_ext is 4
int debug_asg_cut_tip(asg_t *g, int max_ext)
{
	asg64_v a = {0,0,0};
	uint32_t n_vtx = g->n_seq * 2, v, i, cnt = 0;
    
	for (v = 0; v < n_vtx; ++v) {
		//if this seq has been deleted
		if (g->seq[v>>1].del) continue;
		///check if the another direction of v has no overlaps
		///if the self direction of v has no overlaps, we don't have the overlaps of them
		///here is check if the reverse direction of v 
		/**
		 the following first line is to find (means v is a node has no prefix):
		        (v)--->()---->()---->()----->....
         another case is:
            ......()---->()---->()----->()------>(v)
         this case can be found by (v^1), so we don't need to process this case here
		**/
		if (asg_is_utg_end(g, v, 0) != ASG_ET_TIP) continue; // not a tip
        /**
         the following second line is: 
		 (v)--->()---->()---->()----->()
		        |--------max_ext-------|
        **/
       ///that means here is a long tip, which is longer than max_ext
		if (asg_extend(g, v, max_ext, &a) == ASG_ET_MERGEABLE) continue; // not a short unitig

		/**
		 * so combining the last two lines, they are designed to reomve(n(0), n(1), n(2)):
		 *                  		    ----->n(4)
		 *                  		   |
		 * n(0)--->n(1)---->n(2)---->n(3)
		 *                             |
		 *                              ----->n(5)
		**/
		for (i = 0; i < a.n; ++i)
        {
            asg_seq_del(g, (uint32_t)a.a[i]>>1);
            fprintf(stderr, "removed node: %u\n", (uint32_t)a.a[i]>>1);
        }
			
		++cnt;
	}
	free(a.a);
	if (cnt > 0) asg_cleanup(g);
	fprintf(stderr, "[M::%s] cut %d tips\n", __func__, cnt);
	return cnt;
}


// delete short arcs
///for best graph?
int asg_arc_del_short(asg_t *g, float drop_ratio)
{
	uint32_t v, n_vtx = g->n_seq * 2, n_short = 0;
	for (v = 0; v < n_vtx; ++v) {
		asg_arc_t *av = asg_arc_a(g, v);
		uint32_t i, thres, nv = asg_arc_n(g, v);
		///if there is just one overlap, do nothing
		if (nv < 2) continue;
		//av[0] has the most overlap length
		///remove short overlaps
		thres = (uint32_t)(av[0].ol * drop_ratio + .499);
		///av has been sorted by overlap length
		for (i = nv - 1; i >= 1 && av[i].ol < thres; --i);
		for (i = i + 1; i < nv; ++i)
			av[i].del = 1, ++n_short;
	}
	if (n_short) {
		asg_cleanup(g);
		asg_symm(g);
	}
	fprintf(stderr, "[M::%s] removed %d short overlaps\n", __func__, n_short);
	return n_short;
}


inline int check_weak_ma_hit(ma_hit_t_alloc* aim_paf, ma_hit_t_alloc* reverse_paf_list, 
long long weakID, uint32_t w_qs, uint32_t w_qe)
{
    long long i = 0;
    long long strongID, index;
    for (i = 0; i < aim_paf->length; i++)
    {
        ///if this is a strong overlap
        if (
        aim_paf->buffer[i].ml == 1
        && 
        Get_qs(aim_paf->buffer[i]) <= w_qs
        && 
        Get_qe(aim_paf->buffer[i]) >= w_qe)
        {
            strongID = Get_tn(aim_paf->buffer[i]);
            index = get_specific_overlap(&(reverse_paf_list[strongID]), strongID, weakID);
            if(index != -1)
            {
                return 0;
            }
        }


    }

    return 1;
}


inline int check_weak_ma_hit_reverse(ma_hit_t_alloc* r_paf, ma_hit_t_alloc* r_paf_source, 
long long weakID)
{
    long long i = 0;
    long long strongID, index;
    ///all overlaps coming from another haplotye are strong
    for (i = 0; i < r_paf->length; i++)
    {
        strongID = Get_tn(r_paf->buffer[i]);
        index = get_specific_overlap
            (&(r_paf_source[strongID]), strongID, weakID);
        ///must be a strong overlap
        if(index != -1 && r_paf_source[strongID].buffer[index].ml == 1)
        {
            return 0;
        }
    }

    return 1;
}


inline int check_weak_ma_hit_debug(ma_hit_t_alloc* aim_paf, ma_hit_t_alloc* reverse_paf_list, 
long long weakID)
{
    long long i = 0;
    long long strongID, index;
    for (i = 0; i < aim_paf->length; i++)
    {
        ///if this is a strong overlap
        if (aim_paf->buffer[i].ml == 1)
        {
            strongID = Get_tn(aim_paf->buffer[i]);
            index = get_specific_overlap(&(reverse_paf_list[strongID]), strongID, weakID);
            if(index != -1)
            {
                return strongID;
            }
        }


    }

    return 0;
}

// delete short arcs
///for best graph?
int asg_arc_del_short_diploid(asg_t *g, float drop_ratio, ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources)
{
    float second_drop_ratio = 0.3;
	uint32_t v, n_vtx = g->n_seq * 2, n_short = 0;
    uint32_t last_e, flag;
	for (v = 0; v < n_vtx; ++v) 
    {
		asg_arc_t *av = asg_arc_a(g, v);
		uint32_t i, thres, nv = asg_arc_n(g, v);
		///if there is just one overlap, do nothing
		if (nv < 2) continue;
		//av[0] has the most overlap length
		///remove short overlaps
		thres = (uint32_t)(av[0].ol * drop_ratio + .499);
		///av has been sorted by overlap length
		for (i = nv - 1; i >= 1 && av[i].ol < thres; --i);
        last_e = i + 1;

		for (i = i + 1; i < nv; ++i)
			av[i].del = 1, ++n_short;
        
        
        if(nv >= 2 && av[1].del == 1)
        {
            thres = (uint32_t)(av[0].ol * second_drop_ratio + .499);
            if(av[1].ol >= thres)
            {
                ///second longest
                av[1].del = 0;
                --n_short;
                last_e++;
            }
        }

       /**
       if(last_e > 1)
       {
           flag = 0;
           ///at least one exact edge
           for (i = 0; i < last_e; i++)
           {
               if(av[i].el)
               {
                   flag = 1;
                   break;
               }
           }

           if(flag)
           {
               for (i = 0; i < last_e; i++)
               {
                   //drop inexact overlaps
                   if(av[i].el == 0)
                   {
                       av[i].del = 1;
                       ++n_short;
                   }
               }
           }
       }
       **/
        
	}
	///if (n_short) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}
	fprintf(stderr, "[M::%s] removed %d short overlaps\n", __func__, n_short);
	return n_short;
}


// delete short arcs
///for best graph?
int asg_arc_del_short_diploid_unclean(asg_t *g, float drop_ratio, ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources)
{
    double startTime = Get_T();

	uint32_t v, n_vtx = g->n_seq * 2, n_short = 0;
    uint32_t last_e, flag;
	for (v = 0; v < n_vtx; ++v) 
    {
        if (g->seq[v>>1].del) continue;

		asg_arc_t *av = asg_arc_a(g, v);
		uint32_t i, thres, nv = asg_arc_n(g, v);
		///if there is just one overlap, do nothing
		if (nv < 2) continue;
		//av[0] has the most overlap length
		///remove short overlaps
		thres = (uint32_t)(av[0].ol * drop_ratio + .499);
		///av has been sorted by overlap length
		for (i = nv - 1; i >= 1 && av[i].ol < thres; --i);
        last_e = i + 1;

		for (i = i + 1; i < nv; ++i)
			av[i].del = 1, ++n_short;
        
        
        if(nv >= 2 && av[1].del == 1)
        {
            ///second longest
            av[1].del = 0;
            --n_short;
            last_e++;
        }

       /**
       if(last_e > 1)
       {
           flag = 0;
           ///at least one exact edge
           for (i = 0; i < last_e; i++)
           {
               if(av[i].el)
               {
                   flag = 1;
                   break;
               }
           }

           if(flag)
           {
               for (i = 0; i < last_e; i++)
               {
                   //drop inexact overlaps
                   if(av[i].el == 0)
                   {
                       av[i].del = 1;
                       ++n_short;
                   }
               }
           }
       }
       **/
        
	}
	///if (n_short) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}
	fprintf(stderr, "[M::%s] removed %d short overlaps\n", __func__, n_short);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

	return n_short;
}




inline int get_real_length(asg_t *g, uint32_t v, uint32_t* v_s)
{
    uint32_t i, kv = 0;
    for (i = 0, kv = 0; i < asg_arc_n(g, v); i++)
    {
        if(!asg_arc_a(g, v)[i].del)
        {
            if(v_s) v_s[kv] = asg_arc_a(g, v)[i].v;
            kv++;
        }
    }

    return kv;
}

uint32_t detect_single_path_with_dels(asg_t *g, uint32_t begNode, uint32_t* endNode, long long* Len, buf_t* b)
{
    
    uint32_t v = begNode, w;
    uint32_t kv, kw;
    (*Len) = 0;


    while (1)
    {
        (*Len)++;
        kv = get_real_length(g, v, NULL);
        (*endNode) = v;

        if(b) kv_push(uint32_t, b->b, v>>1);

        if(kv == 0)
        {
            return END_TIPS;
        }

        if(kv == 2)
        {
            return TWO_OUTPUT;
        }

        if(kv > 2)
        {
            return MUL_OUTPUT;
        }

        ///up to here, kv=1
        ///kw must >= 1
        get_real_length(g, v, &w);
        kw = get_real_length(g, w^1, NULL);
        v = w;
        (*endNode) = v;


        if(kw == 2)
        {
            (*Len)++;
            if(b) kv_push(uint32_t, b->b, v>>1);
            return TWO_INPUT;
        }

        if(kw > 2)
        {
            (*Len)++;
            if(b) kv_push(uint32_t, b->b, v>>1);
            return MUL_INPUT;
        }   


        if((v>>1) == (begNode>>1))
        {
            return LOOP;
        } 
    }

    return LONG_TIPS;   
}



long long check_if_diploid(uint32_t v1, uint32_t v2, asg_t *g, 
ma_hit_t_alloc* reverse_sources, long long min_edge_length)
{
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    uint32_t convex1, convex2, f1, f2;
    long long l1, l2;

    b_0.b.n = 0;
    b_1.b.n = 0;
    ///uint32_t flag1 = detect_single_path(g, v1, &convex1, &l1, &b_0);
    uint32_t flag1 = detect_single_path_with_dels(g, v1, &convex1, &l1, &b_0);
    
    ///uint32_t flag2 = detect_single_path(g, v2, &convex2, &l2, &b_1);
    uint32_t flag2 = detect_single_path_with_dels(g, v2, &convex2, &l2, &b_1);
    
    if(flag1 == LOOP || flag2 == LOOP)
    {
        return -1;
    }

    if(flag1 != END_TIPS && flag1 != LONG_TIPS)
    {
        l1--;
        b_0.b.n--;
    }

    if(flag2 != END_TIPS && flag2 != LONG_TIPS)
    {
        l2--;
        b_1.b.n--;
    }


    if(l1 <= min_edge_length || l2 <= min_edge_length)
    {
        return -1;
    }

    buf_t* b_min;
    buf_t* b_max;
    if(l1<=l2)
    {
        b_min = &b_0;
        b_max = &b_1;
    }
    else
    {
        b_min = &b_1;
        b_max = &b_0;
    }

    long long i, j, k;
    double max_count = 0;
    double min_count = 0;
    uint32_t qn, tn;
    for (i = 0; i < b_min->b.n; i++)
    {
        qn = b_min->b.a[i];
        for (j = 0; j < reverse_sources[qn].length; j++)
        {
            tn = Get_tn(reverse_sources[qn].buffer[j]);
            if(g->seq[tn].del == 1) continue;
            min_count++;
            for (k = 0; k < b_max->b.n; k++)
            {
                if(b_max->b.a[k]==tn)
                {
                    max_count++;
                    break;
                }
            }
        }
    }


    free(b_0.b.a);
    free(b_1.b.a);

    if(min_count == 0) return -1;
    if(max_count == 0) return 0;
    if(max_count/min_count>0.3) return 1;
    return 0;

}


long long check_if_diploid_aggressive(uint32_t v1, uint32_t v2, asg_t *g, 
ma_hit_t_alloc* reverse_sources, long long min_edge_length)
{
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    uint32_t convex1, convex2, f1, f2;
    long long l1, l2;

    b_0.b.n = 0;
    b_1.b.n = 0;
    ///uint32_t flag1 = detect_single_path(g, v1, &convex1, &l1, &b_0);
    uint32_t flag1 = detect_single_path_with_dels(g, v1, &convex1, &l1, &b_0);
    
    ///uint32_t flag2 = detect_single_path(g, v2, &convex2, &l2, &b_1);
    uint32_t flag2 = detect_single_path_with_dels(g, v2, &convex2, &l2, &b_1);
    
    if(flag1 == LOOP || flag2 == LOOP)
    {
        return -1;
    }

    if(flag1 != END_TIPS && flag1 != LONG_TIPS)
    {
        l1--;
        b_0.b.n--;
    }

    if(flag2 != END_TIPS && flag2 != LONG_TIPS)
    {
        l2--;
        b_1.b.n--;
    }


    if(l1 <= min_edge_length || l2 <= min_edge_length)
    {
        return -1;
    }

    buf_t* b_min;
    buf_t* b_max;
    if(l1<=l2)
    {
        b_min = &b_0;
        b_max = &b_1;
    }
    else
    {
        b_min = &b_1;
        b_max = &b_0;
    }

    long long i, j, k;
    double max_count = 0;
    double min_count = 0;
    uint32_t qn, tn;
    for (i = 0; i < b_min->b.n; i++)
    {
        qn = b_min->b.a[i];
        for (j = 0; j < reverse_sources[qn].length; j++)
        {
            tn = Get_tn(reverse_sources[qn].buffer[j]);
            if(g->seq[tn].del == 1) continue;
            min_count++;
            for (k = 0; k < b_max->b.n; k++)
            {
                if(b_max->b.a[k]==tn)
                {
                    max_count++;
                    break;
                }
            }
        }
    }



    free(b_0.b.a);
    free(b_1.b.a);

    if(min_count == 0) return -1;
    if(max_count == 0) return 0;

    return 1;
    /**
    if(max_count/min_count>0.3) return 1;
    return 0;
    **/

}

long long check_if_diploid_debug(uint32_t v1, uint32_t v2, asg_t *g, 
ma_hit_t_alloc* reverse_sources, long long min_edge_length)
{
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    uint32_t convex1, convex2, f1, f2;
    long long l1, l2;

    b_0.b.n = 0;
    b_1.b.n = 0;
    uint32_t flag1 = detect_single_path(g, v1, &convex1, &l1, &b_0);
    uint32_t flag2 = detect_single_path(g, v2, &convex2, &l2, &b_1);
    
    if(flag1 == LOOP || flag2 == LOOP)
    {
        return -1;
    }

    if(flag1 != END_TIPS && flag1 != LONG_TIPS)
    {
        l1--;
        b_0.b.n--;
    }

    if(flag2 != END_TIPS && flag2 != LONG_TIPS)
    {
        l2--;
        b_1.b.n--;
    }


    if(l1 <= min_edge_length || l2 <= min_edge_length)
    {
        return 0;
    }

    buf_t* b_min;
    buf_t* b_max;
    if(l1<=l2)
    {
        b_min = &b_0;
        b_max = &b_1;
    }
    else
    {
        b_min = &b_1;
        b_max = &b_0;
    }


    fprintf(stderr, "b_0.n: %d\n", b_0.b.n);
    fprintf(stderr, "b_1.n: %d\n", b_1.b.n);

    fprintf(stderr, "b_min.n: %d\n", b_min->b.n);
    fprintf(stderr, "b_max.n: %d\n", b_max->b.n);

    long long i, j, k;
    double max_count = 0;
    double min_count = 0;
    uint32_t qn, tn;
    for (i = 0; i < b_min->b.n; i++)
    {
        qn = b_min->b.a[i];
        for (j = 0; j < reverse_sources[qn].length; j++)
        {
            tn = Get_tn(reverse_sources[qn].buffer[j]);
            if(g->seq[tn].del == 1) continue;
            min_count++;
            for (k = 0; k < b_max->b.n; k++)
            {
                if(b_max->b.a[k]==tn)
                {
                    max_count++;
                    break;
                }
            }
        }
    }


    free(b_0.b.a);
    free(b_1.b.a);

    if(min_count == 0) return -1;
    if(max_count == 0) return 0;
    if(max_count/min_count>0.3) return 1;
    return 0;

}

int asg_arc_del_too_short_overlaps(asg_t *g, long long dropLen, float drop_ratio, ma_hit_t_alloc* reverse_sources, long long min_edge_length)
{
    double startTime = Get_T();

	uint32_t v, v_max, v_maxLen, n_vtx = g->n_seq * 2, n_short = 0;
    long long drop_ratio_Len;
    /**
	for (v = 0; v < n_vtx; ++v) 
    {
        if (g->seq[v>>1].del) continue;
        if (g->seq_vis[v] != 0) continue;

		asg_arc_t *av = asg_arc_a(g, v);
		uint32_t i, thres, nv = asg_arc_n(g, v);
		///if there is just one overlap, do nothing
		if (nv < 2) continue;
		//av[0] has the most overlap length
		///remove short overlaps
        if(av[0].ol < dropLen) continue;

        drop_ratio_Len = av[0].ol * drop_ratio;
        if(dropLen < drop_ratio_Len)
        {
            drop_ratio_Len = dropLen;
        }

		for (i = nv - 1; i >= 1 && av[i].ol < drop_ratio_Len; --i);

		// for (i = i + 1; i < nv; ++i)
		// 	av[i].del = 1, ++n_short;
        for (i = i + 1; i < nv; ++i)
        {
            if(check_if_diploid(av[0].v, av[i].v, g, reverse_sources, min_edge_length) != 1)
            {
                av[i].del = 1, ++n_short;
            }
        }
	}
    **/


    for (v = 0; v < n_vtx; ++v) 
    {
        if (g->seq[v>>1].del) continue;
        if (g->seq_vis[v] != 0) continue;

        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (nv < 2 || g->seq[v>>1].del) continue;



        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;
        v_max = (uint32_t)-1;

        for (i = 0, n_arc = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                if(v_max == (uint32_t)-1)
                {
                    v_max = av[i].v;
                    v_maxLen = av[i].ol;
                    if(v_maxLen < dropLen) break;
                    drop_ratio_Len = v_maxLen * drop_ratio;
                    if(dropLen < drop_ratio_Len)
                    {
                        drop_ratio_Len = dropLen;
                    }
                }
                else if(av[i].ol < drop_ratio_Len && 
                check_if_diploid(v_max, av[i].v, g, reverse_sources, min_edge_length) != 1)
                {
                    av[i].ol = 1;
                    asg_arc_del(g, av[i].v^1, av[i].ul>>32^1, 1);
                    ++n_short;
                }
            }
        }
    }


    asg_cleanup(g);
    asg_symm(g);
	
	fprintf(stderr, "[M::%s] removed %d short overlaps\n", __func__, n_short);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

	return n_short;
}



int asg_arc_del_short_diploid_unclean_exact(asg_t *g, float drop_ratio, ma_hit_t_alloc* sources)
{
	uint32_t v, n_vtx = g->n_seq * 2, n_short = 0;
    uint32_t last_e, flag;
	for (v = 0; v < n_vtx; ++v) 
    {
        if (g->seq[v>>1].del) continue;

		asg_arc_t *av = asg_arc_a(g, v);
		uint32_t i, nv = asg_arc_n(g, v);
		///if there is just one overlap, do nothing
		if (nv < 2) continue;
		///keep the longest one
        for (i = 1; i < nv; i++)
        {
            ///if it is an inexact overlap
            if(av[i].el == 0 && 
              sources[v>>1].is_fully_corrected == 1&& 
              sources[(av[i].v>>1)].is_fully_corrected == 1)
            {
                av[i].del = 1;
                ++n_short;
            }
        }

	}

	if (n_short) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}
	fprintf(stderr, "[M::%s] removed %d inexact overlaps\n", __func__, n_short);
	return n_short;
}

long long single_edge(asg_t *g, uint32_t begNode, long long edgeLen)
{
    
    uint32_t v = begNode;
    uint32_t nv = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v);
    long long rLen = 0;

    while (rLen < edgeLen && nv == 1)
    {
        rLen++;
        

        if(asg_is_single_edge(g, av[0].v, v>>1) != 1)
        {
            return -1;
        }

        if(rLen == edgeLen)
        {
            return rLen;
        }

        v = av[0].v;
        nv = asg_arc_n(g, v);
        av = asg_arc_a(g, v);      
    }

    return -1;   

}


// delete short arcs
///for best graph?
int asg_arc_del_short_diploid_based_on_length_back(asg_t *g, float drop_ratio, ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources)
{
	uint32_t v, n_vtx = g->n_seq * 2, n_short = 0;
    uint32_t last_e, flag;
	for (v = 0; v < n_vtx; ++v) 
    {
		asg_arc_t *av = asg_arc_a(g, v);
		uint32_t i, thres, nv = asg_arc_n(g, v);
		///if there is just one overlap, do nothing
		if (nv < 2) continue;
		//av[0] has the most overlap length
		///remove short overlaps
		thres = (uint32_t)(av[0].ol * drop_ratio + .499);
		///av has been sorted by overlap length
		for (i = nv - 1; i >= 1 && av[i].ol < thres; --i);
        last_e = i + 1;

		for (i = i + 1; i < nv; ++i)
			av[i].del = 1, ++n_short;
        
        
        if(nv >= 2 && av[1].del == 1)
        {
            if(single_edge(g, av[1].v, 4) != -1)
            {
                ///second longest
                av[1].del = 0;
                --n_short;
                last_e++;
            }
        }

        
	}
	///if (n_short) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}
	fprintf(stderr, "[M::%s] removed %d short overlaps\n", __func__, n_short);
	return n_short;
}


///check if v has only one branch
static uint32_t asg_check_unambi1(asg_t *g, uint32_t v)
{
	asg_arc_t *av = asg_arc_a(g, v);
	uint32_t i, nv = asg_arc_n(g, v);
	uint32_t k = nv, kv;
	for (i = 0, kv = 0; i < nv; ++i)
		if (!av[i].del) ++kv, k = i;
	if (kv != 1) return (uint32_t)-1;
	return av[k].v;
}
///to see if it is a long tip
static int asg_topocut_aux(asg_t *g, uint32_t v, int max_ext)
{
	int32_t n_ext;
	for (n_ext = 1; n_ext < max_ext && v != (uint32_t)-1; ++n_ext) {
		if (asg_check_unambi1(g, v^1) == (uint32_t)-1) {
			--n_ext;
			break;
		}
		v = asg_check_unambi1(g, v);
	}
    /**
    if(v == (uint32_t)-1 && n_ext < max_ext) //return max_ext + 1;
    {
        fprintf(stderr, "v: %llu, n_ext: %llu, max_ext: %llu \n", v, n_ext, max_ext);
    }
    **/
    ///if(v == (uint32_t)-1) return max_ext + 1;
	return n_ext;
}

// delete short arcs
///for best graph?
int asg_arc_del_short_diploid_by_length(asg_t *g, float drop_ratio, int max_ext, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen)
{
    double startTime = Get_T();
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));

	uint32_t v, n_vtx = g->n_seq * 2;
    long long n_cut = 0;
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(g->seq_vis[v] == 0)
        {
            asg_arc_t *av = asg_arc_a(g, v);
            uint32_t nv = asg_arc_n(g, v);
            if (nv < 2) continue;
            long long i;
            for (i = 0; i < nv; ++i)
            {
                kv_push(uint64_t, b, (uint64_t)(av[i].ol << 32 | (av - g->arc + i)));   
            }
        }
	}

    fprintf(stderr, "[M::%s] %lld unsorted pending overlaps\n", __func__, b.n);

    radix_sort_arch64(b.a, b.a + b.n);

    fprintf(stderr, "[M::%s] %lld sorted pending overlaps\n", __func__, b.n);

    long long k;
    for (k = 0; k < b.n; k++)
    {
        
        asg_arc_t *a = &g->arc[(uint32_t)b.a[k]];
		///v is self id, w is the id of another end
		uint32_t i, iv, iw, v = (a->ul)>>32, w = a->v^1, to_del = 0;
		uint32_t nv = asg_arc_n(g, v), nw = asg_arc_n(g, w), kv, kw;
		uint32_t ov_max = 0, ov_max_i, ow_max = 0, ow_max_i;
		asg_arc_t *av, *aw;
		///nv must be >= 2
		if (nv == 1 && nw == 1) continue;
        av = asg_arc_a(g, v);
		aw = asg_arc_a(g, w);


        ///calculate the longest edge for v and w
		for (i = 0, kv = 0; i < nv; ++i) {
			if (av[i].del) continue;
			if (ov_max < av[i].ol) ov_max = av[i].ol, ov_max_i = i; 
			++kv;
		}
		if (kv >= 2 && a->ol > ov_max * drop_ratio) continue;


		for (i = 0, kw = 0; i < nw; ++i) {
			if (aw[i].del) continue;
			if (ow_max < aw[i].ol) ow_max = aw[i].ol, ow_max_i = i;
			++kw;
		}
		if (kw >= 2 && a->ol > ow_max * drop_ratio) continue;
		///if (kv == 1 && kw == 1) continue;
        if (kv <= 1 && kw <= 1) continue;


        ///to see which one is the current edge (from v and w)
		for (iv = 0; iv < nv; ++iv)
			if (av[iv].v == (w^1)) break;
		for (iw = 0; iw < nw; ++iw)
			if (aw[iw].v == (v^1)) break;
        ///if one edge has been deleted, it should be deleted in both direction
        if (av[iv].del && aw[iw].del) continue;

        ///kv and kw is the avialiable 
		if (kv > 1 && kw > 1) {
			if (a->ol < ov_max * drop_ratio && a->ol < ow_max * drop_ratio)
				to_del = 1;
            // if(to_del == 1)
            // {
            //     if(check_if_diploid(av[ov_max_i].v, w^1, g, reverse_sources, miniedgeLen) == 1 
            //        || 
            //        check_if_diploid(aw[ow_max_i].v, v^1, g, reverse_sources, miniedgeLen) == 1)
            //     {
            //         to_del = 0;
            //     }
            // }
        
		} else if (kw == 1) {
			if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
            ///kv > 1
            // if(to_del == 1)
            // {
            //     if(check_if_diploid(av[ov_max_i].v, w^1, g, reverse_sources, miniedgeLen) == 1)
            //     {
            //         to_del = 0;
            //     }
            // }
		} else if (kv == 1) {
			if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
            ///kw > 1
            // if(to_del == 1)
            // {
            //     if(check_if_diploid(aw[ow_max_i].v, v^1, g, reverse_sources, miniedgeLen) == 1)
            //     {
            //         to_del = 0;
            //     }
            // }
		}
		if (to_del)
			av[iv].del = aw[iw].del = 1, ++n_cut;
        
    }
    
    free(b.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}
	fprintf(stderr, "[M::%s] removed %d short overlaps\n", __func__, n_cut);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
	return n_cut;
}




// delete short arcs
///for best graph?
int asg_arc_del_short_false_link_back(asg_t *g, float drop_ratio, int max_ext)
{
    double startTime = Get_T();
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));

	uint32_t v, n_vtx = g->n_seq * 2;
    long long n_cut = 0;
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(g->seq_vis[v] == 0)
        {
            asg_arc_t *av = asg_arc_a(g, v);
            uint32_t nv = asg_arc_n(g, v);
            if (nv != 2) continue;
            if(asg_arc_n(g, v^1)!=1) continue;

            
            long long i;
            for (i = 0; i < nv; ++i)
            {
                kv_push(uint64_t, b, (uint64_t)(av[i].ol << 32 | (av - g->arc + i)));   
            }
        }
	}

    fprintf(stderr, "[M::%s] %lld unsorted pending overlaps\n", __func__, b.n);

    radix_sort_arch64(b.a, b.a + b.n);

    fprintf(stderr, "[M::%s] %lld sorted pending overlaps\n", __func__, b.n);

    uint32_t v_s[2];
    uint32_t w_s[4];

    long long k;
    for (k = 0; k < b.n; k++)
    {
        asg_arc_t *a = &g->arc[(uint32_t)b.a[k]];
        if(a->del) continue;
		///v is self id, w is the id of another end
		uint32_t i, iv, v = (a->ul)>>32, to_del = 0;
		uint32_t nv = asg_arc_n(g, v), kv;
		asg_arc_t *av, *aw;
		if (nv < 2) continue;
        av = asg_arc_a(g, v);


        kv = get_real_length(g, v, NULL);
        if (kv != 2) continue;
        
        if(get_real_length(g, v^1, NULL)!=1) continue;


        get_real_length(g, v^1, v_s);
        if(get_real_length(g, v_s[0]^1, NULL) < 2) continue;




        //check the length
        for (i = 0, kv = 0; i < nv; ++i) {
			if (av[i].del) continue;
            v_s[kv] = av[i].ol;
			++kv;
		}

        uint32_t s_max = 0;
        for (i = 0; i < asg_arc_n(g, v^1); i++)
        {
            if(asg_arc_a(g, v^1)[i].del) continue;
            s_max = asg_arc_a(g, v^1)[i].ol;
            break;
        }

        uint32_t ov_max, ov_min;
        if(v_s[0] >= v_s[1])
        {
            ov_max = v_s[0];
            ov_min = v_s[1];
        }
        else
        {
            ov_max = v_s[1];
            ov_min = v_s[0];
        }

        if(ov_min < ov_max * drop_ratio)
        {
            continue;
        }

        if(ov_max > s_max * 0.5)
        {
            continue;
        }
        //check the length
        



        get_real_length(g, v, v_s);
        v_s[0] = v_s[0]^1;
        v_s[1] = v_s[1]^1;

        if(v_s[0] == v_s[1]) continue;

        if(get_real_length(g, v_s[0], NULL)!=2) continue;
        if(get_real_length(g, v_s[1], NULL)!=2) continue;
        

        get_real_length(g, v_s[0], w_s);
        get_real_length(g, v_s[1], w_s + 2);

        for (i = 0; i < 2; i++)
        {
            if((w_s[i]>>1) != (v>>1))
            {
                w_s[0] = w_s[i];
            }
        }

        for (i = 2; i < 4; i++)
        {
            if((w_s[i]>>1) != (v>>1))
            {
                w_s[1] = w_s[i];
            }
        }

        if(w_s[0] == w_s[1])
        { 
            to_del = 1;
        }

        uint32_t convex1, f1;
        long long l1;
        if(to_del == 0)
        {
            f1 = detect_bubble_end_with_bubbles(g, w_s[0], w_s[1], &convex1, &l1, NULL);
            if(f1)
            {
                to_del = 1;
            }
        }

        if(to_del == 0)
        {
            v_s[0] = v_s[0]^1;
            v_s[1] = v_s[1]^1;

            if(get_real_length(g, v_s[0], NULL)!=1) continue;
            if(get_real_length(g, v_s[1], NULL)!=1) continue;

            if(v_s[0] == v_s[1])
            {
                to_del = 1;
            }

            f1 = detect_bubble_end_with_bubbles(g, v_s[0], v_s[1], &convex1, &l1, NULL);
            if(f1)
            {
                to_del = 1;
            }
        }


		if (to_del)
        {
            for (i = 0; i < nv; ++i) 
            {
			    if (av[i].del) continue;
                ++n_cut;
                av[i].del = 1;
                asg_arc_del(g, av[i].v^1, av[i].ul>>32^1, 1);
		    }

        }
    }
    
    free(b.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}
	fprintf(stderr, "[M::%s] removed %d false overlaps\n", __func__, n_cut);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
	return n_cut;
}


int asg_arc_del_short_false_link(asg_t *g, float drop_ratio, float o_drop_ratio, int max_dist, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen)
{
    double startTime = Get_T();
    
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));


    kvec_t(uint32_t) b_f;
    memset(&b_f, 0, sizeof(b_f));

    kvec_t(uint32_t) b_r;
    memset(&b_r, 0, sizeof(b_r));


	uint32_t v, w, n_vtx = g->n_seq * 2, n_cut = 0;
    uint32_t sink;

    buf_t bub;
    if (!g->is_symm) asg_symm(g);
    memset(&bub, 0, sizeof(buf_t));
    bub.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(g->seq_vis[v] == 0)
        {
            asg_arc_t *av = asg_arc_a(g, v);
            uint32_t nv = asg_arc_n(g, v);
            if(nv == 1 && asg_arc_n(g, v^1) == 1) continue;

            uint64_t t_ol = 0;
            long long i;
            for (i = 0; i < nv; ++i)
            {
                t_ol += av[i].ol;
            }
            kv_push(uint64_t, b, (uint64_t)(t_ol << 32 | v));  
        }
	}

    fprintf(stderr, "[M::%s] %lld unsorted pending overlaps\n", __func__, b.n);

    radix_sort_arch64(b.a, b.a + b.n);

    fprintf(stderr, "[M::%s] %lld sorted pending overlaps\n", __func__, b.n);

    

    uint32_t min_edge;
    
    
    long long k, t;
    for (k = 0; k < b.n; k++)
    {
        ///v is the node
        v = (uint32_t)b.a[k];
        if (g->seq[v>>1].del) continue;
        uint32_t nv = asg_arc_n(g, v), nw, to_del_l, to_del_r;
        if (nv < 2) continue;
        uint32_t kv = get_real_length(g, v, NULL), kw;
        if (kv < 2) continue;
        uint32_t i;
        asg_arc_t *av = asg_arc_a(g, v), *aw;

        b_f.n = 0;
        b_r.n = 0;
        to_del_l = 0;
        for (i = 0; i < nv; i++)
        {
            if (av[i].del) continue;
            
            w = av[i].v^1;
            nw = asg_arc_n(g, w);
            if(nw < 2) break;
            kw = get_real_length(g, w, NULL);
            if(kw < 2) break;

            kv_push(uint32_t, b_f, av[i].v);
            kv_push(uint32_t, b_r, w);

            aw = asg_arc_a(g, w);
            min_edge = (u_int32_t)-1;
            for (t = 0; t < nw; t++)
            {
                if(aw[t].del) continue;
                if((aw[t].v>>1) == (v>>1)) continue;
                if(aw[t].ol < min_edge) min_edge = aw[t].ol;
                ///kv_push(uint32_t, b_r, aw[t].v);
            }

            if(av[i].ol < min_edge * drop_ratio) to_del_l++;
        }












        /****************************may have bugs********************************/
        if(to_del_l != kv)
        {
            b_f.n = 0;
            b_r.n = 0;
            to_del_l = 0;

            for (i = 0; i < nv; i++)
            {
                if (av[i].del) continue;

                w = av[i].v^1;
                nw = asg_arc_n(g, w);
                if(nw < 2) break;
                kw = get_real_length(g, w, NULL);
                if(kw < 2) break;

                kv_push(uint32_t, b_f, av[i].v);
                kv_push(uint32_t, b_r, w);

                aw = asg_arc_a(g, w);
                min_edge = (u_int32_t)-1;
                for (t = 0; t < nw; t++)
                {
                    if(aw[t].del) continue;
                    if((aw[t].v>>1) == (v>>1)) continue;
                    if(aw[t].ol < min_edge) min_edge = aw[t].ol;
                }

                if(av[i].ol < min_edge * o_drop_ratio) to_del_l++;
            }

            if(to_del_l == kv)
            {
                ///forward
                to_del_l = 1;
                for (i = 1; i < b_f.n; i++)
                {
                    if(check_if_diploid(b_f.a[0], b_f.a[i], g, reverse_sources, miniedgeLen) == 1)
                    {
                        to_del_l++;
                    }
                }

                ///backward
                if(to_del_l != kv && b_r.n >= 2)
                {
                    to_del_l = 0;
                    uint32_t w0, w1;


                    w = b_r.a[0];
                    kw = get_real_length(g, w, NULL);
                    if(kw != 2) goto terminal;
                    aw = asg_arc_a(g, w);
                    nw = asg_arc_n(g, w);
                    for (t = 0; t < nw; t++)
                    {
                        if(aw[t].del) continue;
                        if((aw[t].v>>1) == (v>>1)) continue;
                        w0 = aw[t].v;
                    }
                    to_del_l = 1;


                    for (i = 1; i < b_r.n; i++)
                    {
                        w = b_r.a[i];
                        kw = get_real_length(g, w, NULL);
                        if(kw != 2) goto terminal;
                        aw = asg_arc_a(g, w);
                        nw = asg_arc_n(g, w);
                        for (t = 0; t < nw; t++)
                        {
                            if(aw[t].del) continue;
                            if((aw[t].v>>1) == (v>>1)) continue;
                            w1 = aw[t].v;
                        }

                        if(check_if_diploid(w0, w1, g, reverse_sources, miniedgeLen) == 1)
                        {
                            to_del_l++;
                        }
                    }

                }
            }
        }

        terminal:
        /****************************may have bugs********************************/


        if(to_del_l != kv) continue;





        uint32_t convex1;
        long long l1;




        ////forward bubble
        to_del_l = 0;
        for (i = 0; i < b_f.n; i++)
        {
            if(b_f.a[i] == b_f.a[0])
            {
                to_del_l = 1;
            } 
            else
            {
                to_del_l = 0;
                break;
            }
        }
        //check the length
        if(to_del_l == 0 && asg_bub_end_finder_with_del_advance(g, 
        b_f.a, b_f.n, max_dist, &bub, 0, (u_int32_t)-1, &sink)==1)
        {
            to_del_l = 1;
        }
        if(to_del_l == 0 && detect_mul_bubble_end_with_bubbles(g, b_f.a, b_f.n, &convex1, &l1, NULL))
        {
            to_del_l = 1;
        }

        ///if(v>>1 == 4581428) fprintf(stderr, "to_del_l: %d, b_f.n: %d\n", to_del_l, b_f.n);

        ///if(to_del_l == 0) continue;







        ////backward bubble
        to_del_r = 0;
        for (i = 0; i < b_r.n; i++)
        {
            if(b_r.a[i] == b_r.a[0])
            {
                to_del_r = 1;
            } 
            else
            {
                to_del_r = 0;
                break;
            }
        }
        if(to_del_r == 0 && asg_bub_end_finder_with_del_advance
        (g, b_r.a, b_r.n, max_dist, &bub, 1, v^1, &sink)==1)
        {
            to_del_r = 1;
        }
        if(to_del_r == 0 && detect_mul_bubble_end_with_bubbles(g, b_r.a, b_r.n, &convex1, &l1, NULL))
        {
            to_del_r = 1;
        }

        
        // if(v>>1 == 4581428)
        // {
        //     fprintf(stderr, "to_del_l: %d, b_f.n: %d\n", to_del_l, b_f.n);
        //     asg_bub_end_finder_with_del_advance_debug(g, b_r.a, b_r.n, max_dist, &bub, 1, v^1);
        // }
        
        // if(v>>1 == 7318036)
        // {
        //     fprintf(stderr, "to_del_l: %d, to_del_r: %d, b_f.n: %d\n", to_del_l, to_del_r, b_f.n);
        // } 



		if (to_del_l && to_del_r)
        {
            for (i = 0; i < nv; ++i) 
            {
			    if (av[i].del) continue;
                ///fprintf(stderr, "%.*s\n", Get_NAME_LENGTH((R_INF), v>>1), Get_NAME((R_INF), v>>1));
                ++n_cut;
                av[i].del = 1;
                asg_arc_del(g, av[i].v^1, av[i].ul>>32^1, 1);
		    }

        }
    }
    
    free(b.a); free(b_f.a); free(b_r.a);
    free(bub.a); free(bub.S.a); free(bub.T.a); free(bub.b.a); free(bub.e.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}
	fprintf(stderr, "[M::%s] removed %d false overlaps\n", __func__, n_cut);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
	return n_cut;
}


int asg_arc_del_short_false_link_advance(asg_t *g, float drop_ratio, float o_drop_ratio, int max_dist, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen)
{
    double startTime = Get_T();
    
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));


    kvec_t(uint32_t) b_f;
    memset(&b_f, 0, sizeof(b_f));

    kvec_t(uint32_t) b_r;
    memset(&b_r, 0, sizeof(b_r));


	uint32_t v, w, n_vtx = g->n_seq * 2, n_cut = 0;
    uint32_t sink;

    buf_t bub;
    if (!g->is_symm) asg_symm(g);
    memset(&bub, 0, sizeof(buf_t));
    bub.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(g->seq_vis[v] == 0)
        {
            asg_arc_t *av = asg_arc_a(g, v);
            uint32_t nv = asg_arc_n(g, v);
            if(nv == 1 && asg_arc_n(g, v^1) == 1) continue;

            uint64_t t_ol = 0;
            long long i;
            for (i = 0; i < nv; ++i)
            {
                t_ol += av[i].ol;
            }
            kv_push(uint64_t, b, (uint64_t)(t_ol << 32 | v));  
        }
	}

    fprintf(stderr, "[M::%s] %lld unsorted pending overlaps\n", __func__, b.n);

    radix_sort_arch64(b.a, b.a + b.n);

    fprintf(stderr, "[M::%s] %lld sorted pending overlaps\n", __func__, b.n);

    

    uint32_t min_edge;
    
    
    long long k, t;
    for (k = 0; k < b.n; k++)
    {
        ///v is the node
        v = (uint32_t)b.a[k];
        if (g->seq[v>>1].del) continue;
        uint32_t nv = asg_arc_n(g, v), nw, to_del_l, to_del_r;
        if (nv < 2) continue;
        uint32_t kv = get_real_length(g, v, NULL), kw;
        if (kv < 2) continue;
        uint32_t i;
        asg_arc_t *av = asg_arc_a(g, v), *aw;

        b_f.n = 0;
        b_r.n = 0;
        to_del_l = 0;
        for (i = 0; i < nv; i++)
        {
            if (av[i].del) continue;
            
            w = av[i].v^1;
            nw = asg_arc_n(g, w);
            if(nw < 2) break;
            kw = get_real_length(g, w, NULL);
            if(kw < 2) break;

            kv_push(uint32_t, b_f, av[i].v);
            kv_push(uint32_t, b_r, w);

            aw = asg_arc_a(g, w);
            min_edge = (u_int32_t)-1;
            for (t = 0; t < nw; t++)
            {
                if(aw[t].del) continue;
                if((aw[t].v>>1) == (v>>1)) continue;
                if(aw[t].ol < min_edge) min_edge = aw[t].ol;
                ///kv_push(uint32_t, b_r, aw[t].v);
            }

            if(av[i].ol < min_edge * drop_ratio) to_del_l++;
        }












        /****************************may have bugs********************************/
        if(to_del_l != kv)
        {
            b_f.n = 0;
            b_r.n = 0;
            to_del_l = 0;

            for (i = 0; i < nv; i++)
            {
                if (av[i].del) continue;

                w = av[i].v^1;
                nw = asg_arc_n(g, w);
                if(nw < 2) break;
                kw = get_real_length(g, w, NULL);
                if(kw < 2) break;

                kv_push(uint32_t, b_f, av[i].v);
                kv_push(uint32_t, b_r, w);

                aw = asg_arc_a(g, w);
                min_edge = (u_int32_t)-1;
                for (t = 0; t < nw; t++)
                {
                    if(aw[t].del) continue;
                    if((aw[t].v>>1) == (v>>1)) continue;
                    if(aw[t].ol < min_edge) min_edge = aw[t].ol;
                }

                if(av[i].ol < min_edge * o_drop_ratio) to_del_l++;
            }

            if(to_del_l == kv)
            {
                ///forward
                to_del_l = 1;
                for (i = 1; i < b_f.n; i++)
                {
                    if(check_if_diploid(b_f.a[0], b_f.a[i], g, reverse_sources, miniedgeLen) == 1)
                    {
                        to_del_l++;
                    }
                }

                ///backward
                if(to_del_l != kv && b_r.n >= 2)
                {
                    to_del_l = 0;
                    uint32_t w0, w1;


                    w = b_r.a[0];
                    kw = get_real_length(g, w, NULL);
                    if(kw != 2) goto terminal;
                    aw = asg_arc_a(g, w);
                    nw = asg_arc_n(g, w);
                    for (t = 0; t < nw; t++)
                    {
                        if(aw[t].del) continue;
                        if((aw[t].v>>1) == (v>>1)) continue;
                        w0 = aw[t].v;
                    }
                    to_del_l = 1;


                    for (i = 1; i < b_r.n; i++)
                    {
                        w = b_r.a[i];
                        kw = get_real_length(g, w, NULL);
                        if(kw != 2) goto terminal;
                        aw = asg_arc_a(g, w);
                        nw = asg_arc_n(g, w);
                        for (t = 0; t < nw; t++)
                        {
                            if(aw[t].del) continue;
                            if((aw[t].v>>1) == (v>>1)) continue;
                            w1 = aw[t].v;
                        }

                        if(check_if_diploid(w0, w1, g, reverse_sources, miniedgeLen) == 1)
                        {
                            to_del_l++;
                        }
                    }

                }
            }
        }

        terminal:
        /****************************may have bugs********************************/


        if(to_del_l != kv) continue;





        uint32_t convex1;
        long long l1;




        ////forward bubble
        to_del_l = 0;
        for (i = 0; i < b_f.n; i++)
        {
            if(b_f.a[i] == b_f.a[0])
            {
                to_del_l = 1;
            } 
            else
            {
                to_del_l = 0;
                break;
            }
        }
        //check the length
        if(to_del_l == 0 && asg_bub_end_finder_with_del_advance(g, 
        b_f.a, b_f.n, max_dist, &bub, 0, (u_int32_t)-1, &sink)==1)
        {
            to_del_l = 1;
        }
        if(to_del_l == 0 && detect_mul_bubble_end_with_bubbles(g, b_f.a, b_f.n, &convex1, &l1, NULL))
        {
            to_del_l = 1;
        }

        ///if(v>>1 == 4581428) fprintf(stderr, "to_del_l: %d, b_f.n: %d\n", to_del_l, b_f.n);

        ///if(to_del_l == 0) continue;







        ////backward bubble
        to_del_r = 0;
        for (i = 0; i < b_r.n; i++)
        {
            if(b_r.a[i] == b_r.a[0])
            {
                to_del_r = 1;
            } 
            else
            {
                to_del_r = 0;
                break;
            }
        }
        if(to_del_r == 0 && asg_bub_end_finder_with_del_advance
        (g, b_r.a, b_r.n, max_dist, &bub, 1, v^1, &sink)==1)
        {
            to_del_r = 1;
        }
        if(to_del_r == 0 && detect_mul_bubble_end_with_bubbles(g, b_r.a, b_r.n, &convex1, &l1, NULL))
        {
            to_del_r = 1;
        }

        
        // if(v>>1 == 4581428)
        // {
        //     fprintf(stderr, "to_del_l: %d, b_f.n: %d\n", to_del_l, b_f.n);
        //     asg_bub_end_finder_with_del_advance_debug(g, b_r.a, b_r.n, max_dist, &bub, 1, v^1);
        // }
        
        



		if (to_del_l && to_del_r)
        {
            fprintf(stderr, "%.*s\n", Get_NAME_LENGTH((R_INF), v>>1), Get_NAME((R_INF), v>>1));
            for (i = 0; i < nv; ++i) 
            {
			    if (av[i].del) continue;
                
                ++n_cut;
                av[i].del = 1;
                asg_arc_del(g, av[i].v^1, av[i].ul>>32^1, 1);
		    }

        }
    }
    
    free(b.a); free(b_f.a); free(b_r.a);
    free(bub.a); free(bub.S.a); free(bub.T.a); free(bub.b.a); free(bub.e.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}
	fprintf(stderr, "[M::%s] removed %d false overlaps\n", __func__, n_cut);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
	return n_cut;
}



int asg_arc_del_tri_link(asg_t *g, int max_dist)
{
    double startTime = Get_T();
    
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));


    kvec_t(uint32_t) b_f;
    memset(&b_f, 0, sizeof(b_f));

    kvec_t(uint32_t) b_r;
    memset(&b_r, 0, sizeof(b_r));


	uint32_t v, w, n_vtx = g->n_seq * 2, n_cut = 0;
    uint32_t sink;

    buf_t bub;
    if (!g->is_symm) asg_symm(g);
    memset(&bub, 0, sizeof(buf_t));
    bub.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    
    uint32_t Ns_first[3];
    uint32_t Ns_second[3];
    long long NodeLen_first[3];
    long long NodeLen_second[3];
    for (v = 0; v < n_vtx; ++v) 
    {
        if(g->seq_vis[v] == 0)
        {
            asg_arc_t *av = asg_arc_a(g, v);
            uint32_t nv = asg_arc_n(g, v);

            if(nv != 2) continue;
            if(av[0].v == av[1].v) continue;

            /**********************test first node************************/
            NodeLen_first[0] = NodeLen_first[1] = NodeLen_first[2] = -1;
            if(asg_is_single_edge(g, av[0].v, v>>1) <= 2 && asg_is_single_edge(g, av[1].v, v>>1) <= 2)
            {
                NodeLen_first[asg_is_single_edge(g, av[0].v, v>>1)] = 0;
                NodeLen_first[asg_is_single_edge(g, av[1].v, v>>1)] = 1;
            }
            ///one node has one out-edge, another node has two out-edges
            if(NodeLen_first[1] == -1 || NodeLen_first[2] == -1)
            {
                continue;
            }
            /**********************test first node************************/



            /**********************test second node************************/
            w = av[NodeLen_first[2]].v^1;
            asg_arc_t *aw = asg_arc_a(g, w);
            uint32_t nw = asg_arc_n(g, w);
            if(nw != 2)
            {
                fprintf(stderr, "error\n");
            }
            NodeLen_second[0] = NodeLen_second[1] = NodeLen_second[2] = -1;
            if(asg_is_single_edge(g, aw[0].v, w>>1) <= 2 && asg_is_single_edge(g, aw[1].v, w>>1) <= 2)
            {
                NodeLen_second[asg_is_single_edge(g, aw[0].v, w>>1)] = 0;
                NodeLen_second[asg_is_single_edge(g, aw[1].v, w>>1)] = 1;
            }
            ///one node has one out-edge, another node has two out-edges
            if(NodeLen_second[1] == -1 || NodeLen_second[2] == -1)
            {
                continue;
            }
            /**********************test second node************************/
            
            uint64_t t_ol = av[NodeLen_first[2]].ol;
            kv_push(uint64_t, b, (uint64_t)(t_ol << 32 | v));  
        }
	}

    fprintf(stderr, "[M::%s] %lld unsorted pending overlaps\n", __func__, b.n);

    radix_sort_arch64(b.a, b.a + b.n);

    fprintf(stderr, "[M::%s] %lld sorted pending overlaps\n", __func__, b.n);

    

    uint32_t min_edge;
    
    
    long long k, t;
    for (k = 0; k < b.n; k++)
    {
        ///v is the node
        v = (uint32_t)b.a[k];
        if (g->seq[v>>1].del) continue;
        uint32_t nv = asg_arc_n(g, v), nw, to_del;
        uint32_t kv = get_real_length(g, v, NULL), kw;


        ///at the begining, the nv of all nodes must be == 2;
        ///here kv == 2, that means all edges are kept
        ///so we can use normal method to delete edges
        if (nv != 2) continue;
        if (kv != 2) continue;
        


        uint32_t i;
        asg_arc_t *av = asg_arc_a(g, v), *aw;
        if(av[0].v == av[1].v)
        {
            continue;
        }

        /**********************test first node************************/
        NodeLen_first[0] = NodeLen_first[1] = NodeLen_first[2] = -1;
        if(asg_is_single_edge(g, av[0].v, v>>1) <= 2 && asg_is_single_edge(g, av[1].v, v>>1) <= 2)
        {
            NodeLen_first[asg_is_single_edge(g, av[0].v, v>>1)] = 0;
            NodeLen_first[asg_is_single_edge(g, av[1].v, v>>1)] = 1;
        }
        ///one node has one out-edge, another node has two out-edges
        if(NodeLen_first[1] == -1 || NodeLen_first[2] == -1)
        {
            continue;
        }
       /**********************test first node************************/

       /**********************test second node************************/
        w = av[NodeLen_first[2]].v^1;
        aw = asg_arc_a(g, w);
        nw = asg_arc_n(g, w);
        kw = get_real_length(g, w, NULL);

        ///at the begining, the nw of all nodes must be == 2;
        ///here kw == 2, that means all edges are kept
        ///so we can use normal method to delete edges
        if(nw != 2) continue;
        if(kw != 2) continue;
        
        
        NodeLen_second[0] = NodeLen_second[1] = NodeLen_second[2] = -1;
        if(asg_is_single_edge(g, aw[0].v, w>>1) <= 2 && asg_is_single_edge(g, aw[1].v, w>>1) <= 2)
        {
            NodeLen_second[asg_is_single_edge(g, aw[0].v, w>>1)] = 0;
            NodeLen_second[asg_is_single_edge(g, aw[1].v, w>>1)] = 1;
        }
        ///one node has one out-edge, another node has two out-edges
        if(NodeLen_second[1] == -1 || NodeLen_second[2] == -1)
        {
            continue;
        }

        /**********************test second node************************/






        uint32_t convex1, convex2, f1, f2;
        long long l1, l2;
        to_del = 0;
        f1 = detect_bubble_end_with_bubbles(g, av[0].v, av[1].v, &convex1, &l1, NULL);
        f2 = detect_bubble_end_with_bubbles(g, aw[0].v, aw[1].v, &convex2, &l2, NULL);
        if(f1 && f2)
        {
            if(l1 <= min_thres || l2 <= min_thres)
            {
                continue;
            } 
            to_del = 1;        
        }
        else if(f1)
        {
            if(l1 <= min_thres)
            {
                continue;
            }
            to_del = 1;   
        }
        else if(f2)
        {
            if(l2 <= min_thres)
            {
                continue;
            }
            to_del = 1;   
        }

    
        if(to_del == 0)
        {
            if(!f1)
            {
                Ns_first[0] = av[0].v; Ns_first[1] = av[1].v;
                f1 = asg_bub_end_finder_with_del_advance(g, Ns_first, 2, max_dist, 
                &bub, 0, (u_int32_t)-1, &convex1);
                l1 = min_thres + 10;
            }

            if(!f2)
            {
                Ns_second[0] = aw[0].v; Ns_second[1] = aw[1].v;
                f2 = asg_bub_end_finder_with_del_advance(g, Ns_second, 2, max_dist, 
                &bub, 0, (u_int32_t)-1, &convex2);
                l2 = min_thres + 10;
            }




            
            if(f1 && f2)
            {
                if(l1 <= min_thres || l2 <= min_thres)
                {
                    continue;
                } 
                to_del = 1;
            }
            else if(f1)
            {

                if(l1 <= min_thres)
                {
                    continue;
                }
                to_del = 1;   
            }
            else if(f2)
            {

                if(l2 <= min_thres)
                {
                    continue;
                }

                to_del = 1;  
            }
            
        }


		if (to_del)
        {
            ++n_cut;
            av[NodeLen_first[2]].del = 1;            
            asg_arc_del(g, av[NodeLen_first[2]].v^1, av[NodeLen_first[2]].ul>>32^1, 1);
        }
    }
    
    free(b.a); free(b_f.a); free(b_r.a);
    free(bub.a); free(bub.S.a); free(bub.T.a); free(bub.b.a); free(bub.e.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}
	fprintf(stderr, "[M::%s] removed %d false overlaps\n", __func__, n_cut);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
	return n_cut;
}

int asg_arc_del_complex_false_link(asg_t *g, float drop_ratio, float o_drop_ratio, int max_dist, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen)
{
    double startTime = Get_T();
    
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));


    kvec_t(uint32_t) b_f;
    memset(&b_f, 0, sizeof(b_f));

    kvec_t(uint32_t) b_r;
    memset(&b_r, 0, sizeof(b_r));


	uint32_t v, w, n_vtx = g->n_seq * 2, n_cut = 0;


    buf_t bub;
    if (!g->is_symm) asg_symm(g);
    memset(&bub, 0, sizeof(buf_t));
    bub.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(g->seq_vis[v] == 0)
        {
            asg_arc_t *av = asg_arc_a(g, v);
            uint32_t nv = asg_arc_n(g, v);

            if(nv == 1 && asg_arc_n(g, v^1) == 1) continue;

            uint64_t t_ol = 0;
            long long i;
            for (i = 0; i < nv; ++i)
            {
                t_ol += av[i].ol;
            }
            kv_push(uint64_t, b, (uint64_t)(t_ol << 32 | v));  
        }
	}

    fprintf(stderr, "[M::%s] %lld unsorted pending overlaps\n", __func__, b.n);

    radix_sort_arch64(b.a, b.a + b.n);

    fprintf(stderr, "[M::%s] %lld sorted pending overlaps\n", __func__, b.n);

    

    uint32_t min_edge;
    
    
    long long k, t;
    for (k = 0; k < b.n; k++)
    {
        ///v is the node
        v = (uint32_t)b.a[k];
        if (g->seq[v>>1].del) continue;
        uint32_t nv = asg_arc_n(g, v), nw, to_del;
        if (nv < 2) continue;
        uint32_t kv = get_real_length(g, v, NULL), kw;
        if (kv < 2) continue;
        uint32_t i;
        asg_arc_t *av = asg_arc_a(g, v), *aw;

        b_f.n = 0;
        b_r.n = 0;
        to_del = 0;
        for (i = 0; i < nv; i++)
        {
            if (av[i].del) continue;
            
            w = av[i].v^1;
            nw = asg_arc_n(g, w);
            if(nw < 2) break;
            kw = get_real_length(g, w, NULL);
            if(kw < 2) break;

            kv_push(uint32_t, b_f, av[i].v);
            kv_push(uint32_t, b_r, w);

            aw = asg_arc_a(g, w);
            min_edge = (u_int32_t)-1;
            for (t = 0; t < nw; t++)
            {
                if(aw[t].del) continue;
                if((aw[t].v>>1) == (v>>1)) continue;
                if(aw[t].ol < min_edge) min_edge = aw[t].ol;
            }

            if(av[i].ol < min_edge * drop_ratio) to_del++;
        }

        if(to_del != kv) continue;

        for (i = 0; i < nv; ++i) 
        {
            if (av[i].del) continue;
            ++n_cut;
            av[i].del = 1;
            asg_arc_del(g, av[i].v^1, av[i].ul>>32^1, 1);
        }
    }


    if(n_cut > 0)
    {
        for (v = 0; v < n_vtx; ++v) 
        {
            uint32_t nv = asg_arc_n(g, v);
            asg_arc_t *av = asg_arc_a(g, v);
            if (g->seq[v>>1].del)
            {
                continue;
            } 

            if(nv < 2)
            {
                continue;
            }


            if(asg_bub_finder_without_del_advance(g, v, max_dist, &bub) == 1)
            {
                uint32_t i;
                g->seq_vis[v] = 3;
                g->seq_vis[v^1] = 3;
                for (i = 0; i < bub.b.n; i++)
                {
                    g->seq_vis[bub.b.a[i]] = 3;
                    g->seq_vis[bub.b.a[i]^1] = 3;
                }
            }
        }


        for (v = 0; v < n_vtx; ++v) 
        {
            if(g->seq_vis[v] == 3) continue;

            uint32_t nv = asg_arc_n(g, v);
            asg_arc_t *av = asg_arc_a(g, v);
            uint32_t i;
            for (i = 0; i < nv; ++i) 
            {
                if (av[i].del && g->seq_vis[av[i].v] != 3)
                {
                    av[i].del = 0;
                    asg_arc_del(g, av[i].v^1, av[i].ul>>32^1, 0);
                } 
            }
        }
    }
    


    
    free(b.a); free(b_f.a); free(b_r.a);
    free(bub.a); free(bub.S.a); free(bub.T.a); free(bub.b.a); free(bub.e.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}
	fprintf(stderr, "[M::%s] removed %d false overlaps\n", __func__, n_cut);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
	return n_cut;
}



int asg_arc_del_short_diploid_by_exact(asg_t *g, int max_ext, ma_hit_t_alloc* sources)
{
    double startTime = Get_T();
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));

	uint32_t v, n_vtx = g->n_seq * 2;
    long long n_cut = 0;
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(g->seq_vis[v] == 0)
        {
            asg_arc_t *av = asg_arc_a(g, v);
            uint32_t nv = asg_arc_n(g, v);
            if (nv < 2) continue;
            long long i;
            for (i = 0; i < nv; ++i)
            {
                kv_push(uint64_t, b, (uint64_t)(av[i].ol << 32 | (av - g->arc + i)));   
            }
        }
	}

    fprintf(stderr, "[M::%s] %lld unsorted pending overlaps\n", __func__, b.n);

    radix_sort_arch64(b.a, b.a + b.n);

    fprintf(stderr, "[M::%s] %lld sorted pending overlaps\n", __func__, b.n);

    long long k;
    for (k = 0; k < b.n; k++)
    {
        
        asg_arc_t *a = &g->arc[(uint32_t)b.a[k]];
		///v is self id, w is the id of another end
		uint32_t i, iv, iw, v = (a->ul)>>32, w = a->v^1, to_del = 0;
		uint32_t nv = asg_arc_n(g, v), nw = asg_arc_n(g, w), kv, kw;
		uint32_t ov_max = 0, ow_max = 0, ov_max_i, ow_max_i;
		asg_arc_t *av, *aw;
		///nv must be >= 2
		if (nv == 1 && nw == 1) continue;
        av = asg_arc_a(g, v);
		aw = asg_arc_a(g, w);


        ///calculate the longest edge for v and w
		for (i = 0, kv = 0; i < nv; ++i) {
			if (av[i].del) continue;
			if (ov_max < av[i].ol) 
            {
                ov_max = av[i].ol;
                ov_max_i = i;
            }
			++kv;
		}
		if (kv >= 2 && a->ol == ov_max) continue;


		for (i = 0, kw = 0; i < nw; ++i) {
			if (aw[i].del) continue;
			if (ow_max < aw[i].ol) 
            {
                ow_max = aw[i].ol;
                ow_max_i = i;
            }
			++kw;
		}
		if (kw >= 2 && a->ol == ow_max) continue;

		///if (kv == 1 && kw == 1) continue;
        if (kv <= 1 && kw <= 1) continue;


        ///to see which one is the current edge (from v and w)
		for (iv = 0; iv < nv; ++iv)
			if (av[iv].v == (w^1)) break;
		for (iw = 0; iw < nw; ++iw)
			if (aw[iw].v == (v^1)) break;
        ///if one edge has been deleted, it should be deleted in both direction
        if (av[iv].del && aw[iw].del) continue;



        if(a->el == 0 && 
        sources[v>>1].is_fully_corrected == 1 && 
        sources[w>>1].is_fully_corrected == 1)
        {
            if (kv > 1 && kw > 1) {
				to_del = 1;
            } else if (kw == 1) {
                if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
            } else if (kv == 1) {
                if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
            }
        }

        if(a->el == 0 && 
        sources[v>>1].is_fully_corrected == 1 && 
        sources[w>>1].is_fully_corrected == 0)
        {
            /****************************may have bugs********************************/
            ///if(av[ov_max_i].el == 1 && sources[av[ov_max_i].v>>1].is_fully_corrected)
            /****************************may have bugs********************************/
            if(av[ov_max_i].el == 1 && sources[av[ov_max_i].v>>1].is_fully_corrected == 1)
            {
                if (kv > 1 && kw > 1) {
                    to_del = 1;
                } else if (kw == 1) {
                    if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
                } else if (kv == 1) {
                    if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
                }
            }
        }


        
		if (to_del)
			av[iv].del = aw[iw].del = 1, ++n_cut;
        
    }
    
    free(b.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}
	fprintf(stderr, "[M::%s] removed %d inexact overlaps\n", __func__, n_cut);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
	return n_cut;
}



int asg_arc_del_short_diploi_by_suspect_edge(asg_t *g, int max_ext, ma_hit_t_alloc* sources)
{
    double startTime = Get_T();
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));

	uint32_t v, n_vtx = g->n_seq * 2;
    long long n_cut = 0;
    
    for (v = 0; v < n_vtx; ++v) 
    {
        ///if(g->seq_vis[v] == 0)
        {
            asg_arc_t *av = asg_arc_a(g, v);
            uint32_t nv = asg_arc_n(g, v);
            if (nv < 2) continue;

            long long i;
            for (i = 0; i < nv; ++i)
            {
                ///means there is a large indel at this edge
                if(av[i].no_l_indel == 0)
                {
                    kv_push(uint64_t, b, (uint64_t)(av[i].ol << 32 | (av - g->arc + i)));   
                }
            }
        }
	}

    fprintf(stderr, "[M::%s] %lld unsorted pending overlaps\n", __func__, b.n);

    radix_sort_arch64(b.a, b.a + b.n);

    fprintf(stderr, "[M::%s] %lld sorted pending overlaps\n", __func__, b.n);

    long long k;
    for (k = 0; k < b.n; k++)
    {
        
        asg_arc_t *a = &g->arc[(uint32_t)b.a[k]];
		///v is self id, w is the id of another end
		uint32_t i, iv, iw, v = (a->ul)>>32, w = a->v^1, to_del = 0;
		uint32_t nv = asg_arc_n(g, v), nw = asg_arc_n(g, w), kv, kw;
		asg_arc_t *av, *aw;
		///nv must be >= 2
		if (nv == 1 && nw == 1) continue;
        av = asg_arc_a(g, v);
		aw = asg_arc_a(g, w);

        ///calculate the longest edge for v and w
		for (i = 0, kv = 0; i < nv; ++i) 
        {
			if (av[i].del) continue;
			++kv;
		}

		for (i = 0, kw = 0; i < nw; ++i) 
        {
			if (aw[i].del) continue;
			++kw;
		}

		if (kv == 1 && kw == 1) continue;

        ///to see which one is the current edge (from v and w)
		for (iv = 0; iv < nv; ++iv)
			if (av[iv].v == (w^1)) break;
		for (iw = 0; iw < nw; ++iw)
			if (aw[iw].v == (v^1)) break;

        ///if one edge has been deleted, it should be deleted in both direction
        if (av[iv].del && aw[iw].del) continue;

        
        if (kv > 1 && kw > 1) {
			to_del = 1;
        } else if (kw == 1) {
            if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
        } else if (kv == 1) {
            if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
        }
        
		if (to_del)
			av[iv].del = aw[iw].del = 1, ++n_cut;
    }
    
    free(b.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}
	fprintf(stderr, "[M::%s] removed %d suspect overlaps\n", __func__, n_cut);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
	return n_cut;
}

int asg_arc_del_false_node(asg_t *g, int max_ext)
{
    double startTime = Get_T();
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));

	uint32_t v, n_vtx = g->n_seq * 2;
    long long n_cut = 0;
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(g->seq_vis[v] == 0)
        {
            if(asg_arc_n(g, v)!=1 || asg_arc_n(g, v^1)!=1)
            {
                continue;
            }


            if(asg_is_single_edge(g, asg_arc_a(g, v)[0].v, v>>1) < 2)
            {
                continue;
            }

            if(asg_is_single_edge(g, asg_arc_a(g, v^1)[0].v, (v^1)>>1) < 2)
            {
                continue;
            }

            if(asg_arc_a(g, v)[0].el == 1)
            {
                continue;
            }

            asg_arc_t *av = asg_arc_a(g, v);
            kv_push(uint64_t, b, (uint64_t)(av[0].ol << 32 | (av - g->arc)));
        }
	}

    radix_sort_arch64(b.a, b.a + b.n);

    long long k;
    ///here all edges are inexact matches 
    for (k = 0; k < b.n; k++)
    {
        
        asg_arc_t *a = &g->arc[(uint32_t)b.a[k]];
		///v is self id, w is the id of another end
		uint32_t i, iv, iw, v = (a->ul)>>32, w = a->v^1, to_del = 0;
		uint32_t nv = asg_arc_n(g, v), nw = asg_arc_n(g, w), kv, kw;
		asg_arc_t *av, *aw;
        av = asg_arc_a(g, v);
		aw = asg_arc_a(g, w);


        /**
        uint32_t en;
        long long pathLen;
        detect_single_path(g, w^1, &en, &pathLen, NULL);
        ///<=2 means there is just one single read from w
        if(pathLen <= 2)
        {
            continue;
        }
        **/



        ///calculate the longest edge for v and w
		for (i = 0, kv = 0; i < nv; ++i) {
			if (av[i].del) continue;
			++kv;
		}

		for (i = 0, kw = 0; i < nw; ++i) {
			if (aw[i].del) continue;
			++kw;
		}

		if (kv < 1 || kw < 2) continue;

        ///to see which one is the current edge (from v and w)
		for (iv = 0; iv < nv; ++iv)
			if (av[iv].v == (w^1)) break;
		for (iw = 0; iw < nw; ++iw)
			if (aw[iw].v == (v^1)) break;
        ///if one edge has been deleted, it should be deleted in both direction
        if (av[iv].del && aw[iw].del) continue;

        

        uint32_t el_edges = 0;
        ///there should be at least two available edges in aw
        for (i = 0; i < nw; i++)
        {
            if (aw[i].del) continue;

            if(i != iw && aw[i].el == 1)
            {
                el_edges++;
            }
        }

        if(el_edges > 0 && av[iv].el == 0)
        {
            asg_seq_del(g, v>>1);
            ++n_cut;
        }
    }
    
    free(b.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}
	fprintf(stderr, "[M::%s] removed %d single nodes\n", __func__, n_cut);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
	return n_cut;
}



#define arc_cnt(g, v) ((uint32_t)(g)->idx[(v)])
#define arc_first(g, v) ((g)->arc[(g)->idx[(v)]>>32])

ma_ug_t *ma_ug_gen(asg_t *g)
{
	int32_t *mark;
	uint32_t i, v, n_vtx = g->n_seq * 2;
	///is a queue
	kdq_t(uint64_t) *q;
	ma_ug_t *ug;

	ug = (ma_ug_t*)calloc(1, sizeof(ma_ug_t));
	ug->g = asg_init();
	mark = (int32_t*)calloc(n_vtx, 4);

	q = kdq_init(uint64_t);
	for (v = 0; v < n_vtx; ++v) {
		uint32_t w, x, l, start, end, len;
		ma_utg_t *p;
        ///what's the usage of mark array
        ///mark array is used to select another direction of node
		if (g->seq[v>>1].del || arc_cnt(g, v) == 0 || mark[v]) continue;
		mark[v] = 1;
		q->count = 0, start = v, end = v^1, len = 0;
		// forward
		w = v;
		while (1) {
			/**
			 * w----->x
			 * w<-----x
			 * that means the only suffix of w is x, and the only prefix of x is w
			 **/
			if (arc_cnt(g, w) != 1) break;
			x = arc_first(g, w).v; // w->x
			if (arc_cnt(g, x^1) != 1) break;

			/**
			 * another direction of w would be marked as used (since w has been used)
			**/
			mark[x] = mark[w^1] = 1;
			///l is the edge length, instead of overlap length
            ///note: edge length is different with overlap length
			l = asg_arc_len(arc_first(g, w));
			kdq_push(uint64_t, q, (uint64_t)w<<32 | l);
			end = x^1, len += l;
			w = x;
			if (x == v) break;
		}
		if (start != (end^1) || kdq_size(q) == 0) { // linear unitig
			///length of seq, instead of edge
			l = g->seq[end>>1].len;
			kdq_push(uint64_t, q, (uint64_t)(end^1)<<32 | l);
			len += l;
		} else { // circular unitig
			start = end = UINT32_MAX;
			goto add_unitig; // then it is not necessary to do the backward
		}
		// backward
		x = v;
		while (1) { // similar to forward but not the same
			if (arc_cnt(g, x^1) != 1) break;
			w = arc_first(g, x^1).v ^ 1; // w->x
			if (arc_cnt(g, w) != 1) break;
			mark[x] = mark[w^1] = 1;
			l = asg_arc_len(arc_first(g, w));
			///w is the seq id + direction, l is the length of edge
			///push element to the front of a queue
			kdq_unshift(uint64_t, q, (uint64_t)w<<32 | l);
			start = w, len += l;
			x = w;
		}
add_unitig:
		if (start != UINT32_MAX) mark[start] = mark[end] = 1;
		kv_pushp(ma_utg_t, ug->u, &p);
		p->s = 0, p->start = start, p->end = end, p->len = len, p->n = kdq_size(q), p->circ = (start == UINT32_MAX);
		p->m = p->n;
		kv_roundup32(p->m);
		p->a = (uint64_t*)malloc(8 * p->m);
		//all elements are saved here
		for (i = 0; i < kdq_size(q); ++i)
			p->a[i] = kdq_at(q, i);
	}
	kdq_destroy(uint64_t, q);

	// add arcs between unitigs; reusing mark for a different purpose
	//ug saves all unitigs
	for (v = 0; v < n_vtx; ++v) mark[v] = -1;


	for (i = 0; i < ug->u.n; ++i) {
		if (ug->u.a[i].circ) continue;
		mark[ug->u.a[i].start] = i<<1 | 0;
		mark[ug->u.a[i].end] = i<<1 | 1;
	}

	//scan all edges
	for (i = 0; i < g->n_arc; ++i) {
		asg_arc_t *p = &g->arc[i];
		if (p->del) continue;
		/**
	p->ul: |____________31__________|__________1___________|______________32_____________|
	                    qns            direction of overlap       length of this node (not overlap length)
						                 (based on query)
	p->v : |___________31___________|__________1___________|
				        tns             reverse direction of overlap
						                  (based on target)
	p->ol: overlap length
	**/
		///to connect two unitigs, we need to connect the end of unitig x to the start of unitig y
		///so we need to ^1 to get the reverse direction of (x's end)?
		if (mark[p->ul>>32^1] >= 0 && mark[p->v] >= 0) {
			asg_arc_t *q;
			uint32_t u = mark[p->ul>>32^1]^1;
			int l = ug->u.a[u>>1].len - p->ol;
			if (l < 0) l = 1;
			q = asg_arc_pushp(ug->g);
			q->ol = p->ol, q->del = 0;
			q->ul = (uint64_t)u<<32 | l;
			q->v = mark[p->v];
		}
	}
	for (i = 0; i < ug->u.n; ++i)
		asg_seq_set(ug->g, i, ug->u.a[i].len, 0);
	asg_cleanup(ug->g);
	free(mark);
	return ug;
}

static char comp_tab[] = { // complement base
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

// generate unitig sequences
int ma_ug_seq(ma_ug_t *g, All_reads *RNF, const ma_sub_t *coverage_cut,
const long long n_read)
{
    UC_Read g_read;
    init_UC_Read(&g_read);
	utg_intv_t *tmp;
	uint32_t i, j;

    

	tmp = (utg_intv_t*)calloc(n_read, sizeof(utg_intv_t));
	///number of unitigs
	for (i = 0; i < g->u.n; ++i) {
		ma_utg_t *u = &g->u.a[i];
		uint32_t l = 0;
		u->s = (char*)calloc(1, u->len + 1);
		memset(u->s, 'N', u->len);
		for (j = 0; j < u->n; ++j) {
			utg_intv_t *t = &tmp[u->a[j]>>33];
			///assert(t->len == 0);
			t->utg = i, t->ori = u->a[j]>>32&1;
			t->start = l, t->len = (uint32_t)u->a[j];
			l += t->len;
		}
	}

    

    int32_t id;
    for (id = 0; id < n_read; id++)
    {
		utg_intv_t *t;
		ma_utg_t *u;
		if (id < 0 || tmp[id].len == 0) continue;

		t = &tmp[id];
		u = &g->u.a[t->utg];
        recover_UC_Read(&g_read, RNF, id);
		
        memmove(g_read.seq, g_read.seq + coverage_cut[id].s, coverage_cut[id].e - coverage_cut[id].s);
        g_read.length = coverage_cut[id].e - coverage_cut[id].s;

		if (!t->ori) { // forward strand
			for (i = 0; i < t->len; ++i)
				u->s[t->start + i] = g_read.seq[i];
		} else {
			for (i = 0; i < t->len; ++i) {
				int c = (uint8_t)g_read.seq[g_read.length - 1 - i];
				u->s[t->start + i] = c >= 128? 'N' : comp_tab[c];
			}
		}
    }
    

	free(tmp);

    destory_UC_Read(&g_read);
	return 0;
}

void ma_ug_print(const ma_ug_t *ug, All_reads *RNF, const ma_sub_t *coverage_cut, FILE *fp)
{
	uint32_t i, j, l;
	char name[32];
	for (i = 0; i < ug->u.n; ++i) { // the Segment lines in GFA
		ma_utg_t *p = &ug->u.a[i];
		sprintf(name, "utg%.6d%c", i + 1, "lc"[p->circ]);
		fprintf(fp, "S\t%s\t%s\tLN:i:%d\n", name, p->s? p->s : "*", p->len);
        

    
		for (j = l = 0; j < p->n; l += (uint32_t)p->a[j++]) {
			uint32_t x = p->a[j]>>33;
            fprintf(fp, "a\t%s\t%d\t%.*s(%u):%d-%d\t%c\t%d\n", name, l, 
            Get_NAME_LENGTH((*RNF), x), Get_NAME((*RNF), x), x,
            coverage_cut[x].s + 1, coverage_cut[x].e, "+-"[p->a[j]>>32&1], (uint32_t)p->a[j]);
            
			// if (sub) fprintf(fp, "a\t%s\t%d\t%s:%d-%d\t%c\t%d\n", name, l, d->seq[x].name, sub[x].s + 1, sub[x].e, "+-"[p->a[j]>>32&1], (uint32_t)p->a[j]);
			// else fprintf(fp, "a\t%s\t%d\t%s\t%c\t%d\n", name, l, d->seq[x].name, "+-"[p->a[j]>>32&1], (uint32_t)p->a[j]);
            
        }
        
        
	}
	for (i = 0; i < ug->g->n_arc; ++i) { // the Link lines in GFA
		uint32_t u = ug->g->arc[i].ul>>32, v = ug->g->arc[i].v;
		fprintf(fp, "L\tutg%.6d%c\t%c\tutg%.6d%c\t%c\t%dM\tSD:i:%d\n", (u>>1)+1, "lc"[ug->u.a[u>>1].circ], "+-"[u&1],
				(v>>1)+1, "lc"[ug->u.a[v>>1].circ], "+-"[v&1], ug->g->arc[i].ol, asg_arc_len(ug->g->arc[i]));
	}


    /**
	for (i = 0; i < ug->u.n; ++i) { // summary of unitigs
		uint32_t cnt[2];
		ma_utg_t *u = &ug->u.a[i];
		if (u->start == UINT32_MAX) {
			fprintf(fp, "x\tutg%.6dc\t%d\t%d\n", i + 1, u->len, u->n);
		} else 
        {
			for (j = 0; j < 2; ++j) cnt[j] = asg_arc_n(ug->g, i<<1|j);

            fprintf(fp, "x\tutg%.6dl\t%d\t%d\t%d\t%d\t%.*s:%d-%d\t%c\t%.*s:%d-%d\t%c\n", 
                        i + 1, u->len, u->n, cnt[1], cnt[0],
						///d->seq[u->start>>1].name, 
                        Get_NAME_LENGTH((*RNF), u->start>>1), Get_NAME((*RNF), u->start>>1),
                        coverage_cut[u->start>>1].s + 1, coverage_cut[u->start>>1].e, 
                        "+-"[u->start&1],
						///d->seq[u->end>>1].name, 
                        Get_NAME_LENGTH((*RNF), u->end>>1), Get_NAME((*RNF), u->end>>1),
                        coverage_cut[u->end>>1].s + 1, coverage_cut[u->end>>1].e, 
                        "+-"[u->end&1]);
            
			// if (sub)
			// 	fprintf(fp, "x\tutg%.6dl\t%d\t%d\t%d\t%d\t%s:%d-%d\t%c\t%s:%d-%d\t%c\n", i + 1, u->len, u->n, cnt[1], cnt[0],
			// 			d->seq[u->start>>1].name, sub[u->start>>1].s + 1, sub[u->start>>1].e, "+-"[u->start&1],
			// 			d->seq[u->end>>1].name, sub[u->end>>1].s + 1, sub[u->end>>1].e, "+-"[u->end&1]);
			// else
			// 	fprintf(fp, "x\tutg%.6dl\t%d\t%d\t%d\t%d\t%s\t%c\t%s\t%c\n", i + 1, u->len, u->n, cnt[1], cnt[0],
			// 			d->seq[u->start>>1].name, "+-"[u->start&1], d->seq[u->end>>1].name, "+-"[u->end&1]);
        
        
        }
	}
    **/
    
}

int asg_cut_internal(asg_t *g, int max_ext)
{
	asg64_v a = {0,0,0};
	uint32_t n_vtx = g->n_seq * 2, v, i, cnt = 0;
	for (v = 0; v < n_vtx; ++v) {
		if (g->seq[v>>1].del) continue;
		if (asg_is_utg_end(g, v, 0) != ASG_ET_MULTI_NEI) continue;
		if (asg_extend(g, v, max_ext, &a) != ASG_ET_MULTI_NEI) continue;
		/**
		 * so combining the last two lines, they are designed to reomve(n(1), n(2))?
							   ----->              <-------
							  |                		      |
							 n(0)--->n(1)---->n(2)---->n(3)
							  |                           |
							   ------>             <-------       
		**/
		for (i = 0; i < a.n; ++i)
			asg_seq_del(g, (uint32_t)a.a[i]>>1);
		++cnt;
	}
	free(a.a);
	if (cnt > 0) asg_cleanup(g);
	fprintf(stderr, "[M::%s] cut %d internal sequences\n", __func__, cnt);
	return cnt;
}






void clean_weak_ma_hit_t(ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, long long num_sources)
{
    double startTime = Get_T();
    long long bi_overlaps = 0;
    long long si_overlaps = 0;
    long long i, j, index;
    uint32_t qn, tn;
    ma_hit_t new_element;
    long long qLen_0, qLen_1;

    // if(memcmp("m64016_190918_162737/92668450/ccs", Get_NAME(R_INF, i), 
    //     Get_NAME_LENGTH(R_INF, i)) == 0)
    // {
    //     debug_info_of_specfic_read("m64016_190918_162737/92668450/ccs", 
    //     sources, reverse_sources, -1, "clean_weak_ma_hit_t");

    //     debug_info_of_specfic_read("m64016_190918_162737/53545052/ccs", 
    //     sources, reverse_sources, -1, "clean_weak_ma_hit_t");
    // }

    for (i = 0; i < num_sources; i++)
    {
        for (j = 0; j < sources[i].length; j++)
        {
            qn = Get_qn(sources[i].buffer[j]);
            tn = Get_tn(sources[i].buffer[j]);

            //if this is a weak overlap
            if(sources[i].buffer[j].ml == 0)
            {   
                if(
                !check_weak_ma_hit(&(sources[qn]), reverse_sources, tn, 
                Get_qs(sources[i].buffer[j]), Get_qe(sources[i].buffer[j]))
                /**
                ||
                !check_weak_ma_hit_reverse(&(reverse_sources[qn]), sources, tn)**/)
                {
                    sources[i].buffer[j].bl = 0;
                    index = get_specific_overlap(&(sources[tn]), tn, qn);
                    sources[tn].buffer[index].bl = 0;
                }
            }
        }
    }



    long long m = 0;
    long long pre_overlaps, current_overlaps, exact_overlaps;
    exact_overlaps = pre_overlaps = current_overlaps = 0;
    for (i = 0; i < num_sources; i++)
    {
        m = 0;
        
        for (j = 0; j < sources[i].length; j++)
        {
            if(sources[i].buffer[j].bl != 0)
            {
                sources[i].buffer[m] = sources[i].buffer[j];
                if(sources[i].buffer[m].el)
                {
                    exact_overlaps++;
                }
                m++;
            }
        }

        pre_overlaps += sources[i].length;
        sources[i].length = m;
        current_overlaps += sources[i].length;
    }


    // if(memcmp("m64016_190918_162737/92668450/ccs", Get_NAME(R_INF, i), 
    //     Get_NAME_LENGTH(R_INF, i)) == 0)
    // {
    //     debug_info_of_specfic_read("m64016_190918_162737/92668450/ccs", 
    //     sources, reverse_sources, -1, "clean_weak_ma_hit_t");

    //     debug_info_of_specfic_read("m64016_190918_162737/53545052/ccs", 
    //     sources, reverse_sources, -1, "clean_weak_ma_hit_t");
    // }


    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
}


void build_string_graph(int min_dp, ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, long long n_read, uint64_t* readLen, 
long long mini_overlap_length, long long max_hang_length,
long long clean_round, float min_ovlp_drop_ratio, float max_ovlp_drop_ratio,
float final_ovlp_drop_ratio, char* output_file_name, long long bubble_dist)
{
    ma_sub_t* coverage_cut;
    normalize_ma_hit_t(sources, n_read);
    ///normalize_ma_hit_t_single_side(sources, n_read);
    ///debug_normalize_ma_hit_t(sources, n_read);
    clean_weak_ma_hit_t(sources, reverse_sources, n_read);
    ///debug_normalize_ma_hit_t(sources, n_read);

    ma_hit_sub(min_dp, sources, n_read, readLen, mini_overlap_length, &coverage_cut);
    ma_hit_cut(min_dp, sources, n_read, readLen, mini_overlap_length, &coverage_cut);
    ///it seems we do not need ma_hit_flt
    ma_hit_flt(sources, n_read, coverage_cut, max_hang_length, mini_overlap_length);
    ///debug_cut_ma_hit_t(sources, n_read, coverage_cut);
    ma_hit_contained(sources, n_read, coverage_cut, max_hang_length, mini_overlap_length);
    ///debug_cut_ma_hit_t(sources, n_read, coverage_cut);
    asg_t *sg = NULL;
    sg = ma_sg_gen(sources, n_read, coverage_cut, max_hang_length, mini_overlap_length);
    asg_arc_del_trans(sg, GAP_FUZZ);
    asg_cut_tip(sg, MAX_SHORT_TIPS);

    asg_arc_del_triangular_advance(sg, bubble_dist);
    

    if (asg_arc_del_short_diploid(sg, final_ovlp_drop_ratio, sources, reverse_sources) != 0) 
    {
        ///asg_cut_tip(sg, MAX_SHORT_TIPS);
        
	}

    asg_arc_del_triangular_advance(sg, bubble_dist);

    while(asg_cut_tip(sg, MAX_SHORT_TIPS)!=0 && asg_arc_del_triangular_advance(sg, bubble_dist) != 0)
    {
        ;
    }

    



    ma_ug_t *ug = NULL;
    ug = ma_ug_gen(sg);
    ma_ug_seq(ug, &R_INF, coverage_cut, n_read);


    fprintf(stdout, "Writing GFA to disk ...... \n");
    char* gfa_name = (char*)malloc(strlen(output_file_name)+5);
    sprintf(gfa_name, "%s.gfa", output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    ma_ug_print(ug, &R_INF, coverage_cut, output_file);


    asg_destroy(sg);
	ma_ug_destroy(ug);
    free(coverage_cut);

    free(gfa_name);
    fclose(output_file);
}



void debug_info_of_specfic_node(char* name, asg_t *g, char* command)
{
    fprintf(stderr, "\n\n\n");
    uint32_t v, n_vtx = g->n_seq * 2;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(memcmp(name, Get_NAME(R_INF, (v>>1)), Get_NAME_LENGTH(R_INF, (v>>1))) == 0)
        {
            fprintf(stderr, "\nafter %s\n****************graph ref_read: %.*s, dir: %u****************\n", 
            command, Get_NAME_LENGTH(R_INF, (v>>1)), Get_NAME(R_INF, (v>>1)), v&1);
            if(g->seq[v>>1].del)
            {
                fprintf(stderr, "read has already been deleted.\n");
                continue;
            }

            asg_arc_t *av = asg_arc_a(g, v);
            uint32_t i, nv = asg_arc_n(g, v);
            for (i = 0; i < nv; ++i)
            {
                fprintf(stderr, "target: %.*s, el: %u, strong: %u, ol: %u, del: %u\n", 
                Get_NAME_LENGTH(R_INF, (av[i].v>>1)), 
                Get_NAME(R_INF, (av[i].v>>1)),
                av[i].el, av[i].strong, av[i].ol, av[i].del);
            }
        }
    }
}

void debug_info_of_specfic_read(char* name, ma_hit_t_alloc* sources, 
ma_hit_t_alloc* reverse_sources, int id, char* command)
{
    long long i, j, Len;
    uint32_t qn, tn;

    if(id == -1)
    {
        i = 0;
        Len = R_INF.total_reads;
    }
    else
    {
        i = id;
        Len = id + 1;
    }
    

    for (; i < Len; i++)
    {
        if(memcmp(name, Get_NAME(R_INF, i), Get_NAME_LENGTH(R_INF, i)) == 0)
        {
            fprintf(stderr, "\n\n\nafter %s\n", command);

            fprintf(stderr, "****************ma_hit_t ref_read: %.*s****************\n", 
            Get_NAME_LENGTH(R_INF, i), Get_NAME(R_INF, i));


            fprintf(stderr, "sources Len: %d, is_fully_corrected: %d\n", 
            sources[i].length, sources[i].is_fully_corrected);

            for (j = 0; j < sources[i].length; j++)
            {
                qn = Get_qn(sources[i].buffer[j]);
                tn = Get_tn(sources[i].buffer[j]);
                fprintf(stderr, "target: %.*s, qs: %d, qe: %d, ts: %d, te: %d, ml: %d, rev: %d, el: %d\n", 
                Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn),
                Get_qs(sources[i].buffer[j]),
                Get_qe(sources[i].buffer[j]),
                Get_ts(sources[i].buffer[j]),
                Get_te(sources[i].buffer[j]),
                sources[i].buffer[j].ml,
                sources[i].buffer[j].rev,
                sources[i].buffer[j].el);

                

                /**
                //if this is a weak overlap
                if(sources[i].buffer[j].ml == 0)
                {
                    if(!check_weak_ma_hit(&(sources[qn]), reverse_sources, tn,
                    Get_qs(sources[i].buffer[j]), Get_qe(sources[i].buffer[j])))
                    {
                        int c_id = check_weak_ma_hit_debug(&(sources[qn]), reverse_sources, tn);
                        fprintf(stderr, "*************************conflict with %.*s\n", 
                        Get_NAME_LENGTH(R_INF, c_id), Get_NAME(R_INF, c_id));
                    }
                }
                **/
            }


            
            fprintf(stderr, "######reverse_query_read Len: %d\n", reverse_sources[i].length);



            for (j = 0; j < reverse_sources[i].length; j++)
            {
                qn = Get_qn(reverse_sources[i].buffer[j]);
                tn = Get_tn(reverse_sources[i].buffer[j]);
                fprintf(stderr, "target: %.*s, qs: %u, qe: %u, ts: %u, te: %u\n", 
                Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn),
                Get_qs(reverse_sources[i].buffer[j]),
                Get_qe(reverse_sources[i].buffer[j]),
                Get_ts(reverse_sources[i].buffer[j]),
                Get_te(reverse_sources[i].buffer[j]));
            }
            
        }
    }

    fflush(stderr);

    
}

void ma_sg_print(const asg_t *g, const All_reads *RNF, const ma_sub_t *sub, FILE *fp)
{
	uint32_t i;
    for (i = 0; i < g->n_seq; ++i) 
    {
        if(!g->seq[i].del)
        {
             fprintf(fp, 
             "S\t%.*s\t*\tLN:i:%d\n",
             Get_NAME_LENGTH((*RNF), i), 
             Get_NAME((*RNF), i),
             g->seq[i].len);
        }
    }

	for (i = 0; i < g->n_arc; ++i) {
		const asg_arc_t *p = &g->arc[i];
		if (sub) {
			const ma_sub_t *sq = &sub[p->ul>>33], *st = &sub[p->v>>1];
            /**
			fprintf(fp, "L\t%s:%d-%d\t%c\t%s:%d-%d\t%c\t%d:\tL1:i:%d\n", 
            d->seq[p->ul>>33].name, sq->s + 1, sq->e, "+-"[p->ul>>32&1],
					d->seq[p->v>>1].name, st->s + 1, st->e, "+-"[p->v&1], p->ol, (uint32_t)p->ul);
            **/
           /**
           fprintf(stderr, "Get_NAME_LENGTH((*RNF), p->ul>>33): %u, p->ul>>33: %u\n", 
           Get_NAME_LENGTH((*RNF), (p->ul>>33)),
           (p->ul>>33));

           fprintf(stderr, "Get_NAME_LENGTH((*RNF), p->v>>1): %u, p->v>>1: %u\n", 
           Get_NAME_LENGTH((*RNF), (p->v>>1)),
           (p->v>>1));
           **/

            fprintf(fp, 
            "L\t%.*s:%d-%d\t%c\t%.*s:%d-%d\t%c\t%d:\tL1:i:%d\n", 
            Get_NAME_LENGTH((*RNF), p->ul>>33), 
            Get_NAME((*RNF), p->ul>>33),
            sq->s + 1, sq->e, "+-"[p->ul>>32&1],
            Get_NAME_LENGTH((*RNF), p->v>>1), 
            Get_NAME((*RNF), p->v>>1),
            st->s + 1, st->e, "+-"[p->v&1], p->ol, (uint32_t)p->ul);
            

		} 
        else 
        {
            /**
			fprintf(fp, "L\t%s\t%c\t%s\t%c\t%d:\tL1:i:%d\n", 
            d->seq[p->ul>>33].name, "+-"[p->ul>>32&1],
					d->seq[p->v>>1].name, "+-"[p->v&1], p->ol, (uint32_t)p->ul);
            **/
           fprintf(fp, "L\t%.*s\t%c\t%.*s\t%c\t%d:\tL1:i:%d\n", 
            Get_NAME_LENGTH((*RNF), p->ul>>33), 
            Get_NAME((*RNF), p->ul>>33),
            "+-"[p->ul>>32&1],
            Get_NAME_LENGTH((*RNF), p->v>>1), 
            Get_NAME((*RNF), p->v>>1),
            "+-"[p->v&1], p->ol, (uint32_t)p->ul);
		}
	}
}


void ma_ug_print_simple(const ma_ug_t *ug, All_reads *RNF, const ma_sub_t *coverage_cut, FILE *fp)
{
	uint32_t i, j, l;
	char name[32];
	for (i = 0; i < ug->u.n; ++i) { // the Segment lines in GFA
		ma_utg_t *p = &ug->u.a[i];
		sprintf(name, "utg%.6d%c", i + 1, "lc"[p->circ]);
		fprintf(fp, "S\t%s\t%s\tLN:i:%d\n", name, "*", p->len);        
	}
	for (i = 0; i < ug->g->n_arc; ++i) { // the Link lines in GFA
		uint32_t u = ug->g->arc[i].ul>>32, v = ug->g->arc[i].v;
		fprintf(fp, "L\tutg%.6d%c\t%c\tutg%.6d%c\t%c\t%dM\tSD:i:%d\n", (u>>1)+1, "lc"[ug->u.a[u>>1].circ], "+-"[u&1],
				(v>>1)+1, "lc"[ug->u.a[v>>1].circ], "+-"[v&1], ug->g->arc[i].ol, asg_arc_len(ug->g->arc[i]));
	}    
}


int asg_arc_cut_long_tip(asg_t *g, float drop_ratio)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, v_max, w, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag;
    long long ll, v_maxLen;

    buf_t b;
    memset(&b, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (nv < 2 || g->seq[v>>1].del) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;

        v_maxLen = -1;

        for (i = 0, n_arc = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                flag = detect_single_path_with_dels(g, av[i].v, &convex, &ll, NULL);
                if(v_maxLen < ll)
                {
                    v_maxLen = ll;
                }
            }
        }

        for (i = 0, n_arc = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                b.b.n = 0;
                if(detect_single_path_with_dels(g, av[i].v, &convex, &ll, &b) == END_TIPS)
                {
                    if(v_maxLen*drop_ratio > ll)
                    {
                        n_reduced++;
                        long long k;
                        for (k = 0; k < b.b.n; k++)
                        {
                            asg_seq_del(g, b.b.a[k]);
                        }
                    }
                }
            }
        }
    }


    asg_cleanup(g);
    asg_symm(g);
    

    fprintf(stderr, "[M::%s] removed %d long tips\n", 
    __func__, n_reduced);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

    return n_reduced;
}

uint32_t detect_single_path_with_dels_contigLen(asg_t *g, uint32_t begNode, uint32_t* endNode, long long* baseLen, buf_t* b)
{
    
    uint32_t v = begNode, w;
    uint32_t kv, kw, k;
    (*baseLen) = 0;


    while (1)
    {
        ///(*Len)++;
        kv = get_real_length(g, v, NULL);
        (*endNode) = v;

        if(b) kv_push(uint32_t, b->b, v>>1);

        if(kv == 0)
        {
            (*baseLen) += g->seq[v>>1].len;
            return END_TIPS;
        }

        if(kv == 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            return TWO_OUTPUT;
        }

        if(kv > 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            return MUL_OUTPUT;
        }

        ///kv must be 1
        for (k = 0; k < asg_arc_n(g, v); k++)
        {
            if(!asg_arc_a(g, v)[k].del)
            {
                w = asg_arc_a(g, v)[k].v;
                (*baseLen) += ((uint32_t)(asg_arc_a(g, v)[k].ul));
                break;
            }
        }

        ///up to here, kv=1
        ///kw must >= 1
        kw = get_real_length(g, w^1, NULL);
        v = w;
        (*endNode) = v;


        if(kw == 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            if(b) kv_push(uint32_t, b->b, v>>1);
            return TWO_INPUT;
        }

        if(kw > 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            if(b) kv_push(uint32_t, b->b, v>>1);
            return MUL_INPUT;
        }   


        if((v>>1) == (begNode>>1))
        {
            return LOOP;
        } 
    }

    return LONG_TIPS;   
}


int asg_arc_cut_long_equal_tips_only_tips(asg_t *g, ma_hit_t_alloc* reverse_sources, long long miniedgeLen)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, v_max, w, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag;
    long long ll, base_maxLen, base_maxLen_i;

    buf_t b;
    memset(&b, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (nv < 2 || g->seq[v>>1].del) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;

        base_maxLen = -1;
        base_maxLen_i = -1;

        for (i = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                flag = detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, NULL);
                if(flag != END_TIPS)
                {
                    base_maxLen = -1;
                    base_maxLen_i = -1;
                    break;
                }
                if(base_maxLen < ll)
                {
                    base_maxLen = ll;
                    base_maxLen_i = i;
                }
            }
        }

        ///all outedges are tips
        if(base_maxLen != -1)
        {
            for (i = 0; i < nv; i++)
            {
                if(i == base_maxLen_i) continue;
                if (!av[i].del)
                {
                    ///check_if_diploid_aggressive
                    if(check_if_diploid(av[base_maxLen_i].v, av[i].v, g, 
                    reverse_sources, miniedgeLen)==1)
                    {
                        b.b.n = 0;
                        if(detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b) 
                        == END_TIPS)
                        {
                            n_reduced++;
                            long long k;
                            for (k = 0; k < b.b.n; k++)
                            {
                                asg_seq_del(g, b.b.a[k]);
                            }
                        }
                    }
                }
            }
        }
    }


    asg_cleanup(g);
    asg_symm(g);
    

    fprintf(stderr, "[M::%s] removed %d long tips\n", 
    __func__, n_reduced);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

    return n_reduced;
}

int asg_arc_cut_long_equal_tips(asg_t *g, ma_hit_t_alloc* reverse_sources, long long miniedgeLen)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, v_max, w, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag;
    long long ll, base_maxLen, base_maxLen_i;

    buf_t b;
    memset(&b, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v), n_tips;
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (nv < 2 || g->seq[v>>1].del) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;

        base_maxLen = -1;
        base_maxLen_i = -1;
        n_tips = 0;

        for (i = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                flag = detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, NULL);
                
                if(base_maxLen < ll)
                {
                    base_maxLen = ll;
                    base_maxLen_i = i;
                }

                if(flag == END_TIPS)
                {
                    n_tips++;
                }
            }
        }

        ///at least one tip
        if(n_tips > 0)
        {
            for (i = 0; i < nv; i++)
            {
                if(i == base_maxLen_i) continue;
                if (!av[i].del)
                {
                    b.b.n = 0;
                    if(detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b) 
                    != END_TIPS)
                    {
                        continue;
                    }
                    //we can only cut tips
                    if(check_if_diploid(av[base_maxLen_i].v, av[i].v, g, 
                    reverse_sources, miniedgeLen)==1)
                    {
                        n_reduced++;
                        long long k;
                        for (k = 0; k < b.b.n; k++)
                        {
                            asg_seq_del(g, b.b.a[k]);
                        }
                    }
                }
            }
        }
    }


    asg_cleanup(g);
    asg_symm(g);
    

    fprintf(stderr, "[M::%s] removed %d long tips\n", 
    __func__, n_reduced);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

    return n_reduced;
}

void output_unitig_graph(asg_t *sg, ma_sub_t* coverage_cut, char* output_file_name, long long n_read)
{
    ma_ug_t *ug = NULL;
    ug = ma_ug_gen(sg);
    ma_ug_seq(ug, &R_INF, coverage_cut, n_read);

    fprintf(stdout, "Writing unitig GFA to disk ...... \n");
    char* gfa_name = (char*)malloc(strlen(output_file_name)+25);
    sprintf(gfa_name, "%s.gfa", output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    ma_ug_print(ug, &R_INF, coverage_cut, output_file);
    fclose(output_file);
    
    



    sprintf(gfa_name, "%s.simple.gfa", output_file_name);
    output_file = fopen(gfa_name, "w");
    ma_ug_print_simple(ug, &R_INF, coverage_cut, output_file);
    fclose(output_file);

    free(gfa_name);
    ma_ug_destroy(ug);
}


void output_read_graph(asg_t *sg, ma_sub_t* coverage_cut, char* output_file_name, long long n_read)
{
    fprintf(stdout, "Writing read GFA to disk ...... \n");
    char* gfa_name = (char*)malloc(strlen(output_file_name)+25);
    sprintf(gfa_name, "%s.read.gfa", output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    ma_sg_print(sg, &R_INF, coverage_cut, output_file);
    free(gfa_name);
    fclose(output_file);
}


void read_ma(ma_hit_t* x, FILE* fp)
{
    fread(&(x->qns), sizeof(x->qns), 1, fp);
    fread(&(x->qe), sizeof(x->qe), 1, fp);
    fread(&(x->tn), sizeof(x->tn), 1, fp);
    fread(&(x->ts), sizeof(x->ts), 1, fp);
    fread(&(x->te), sizeof(x->te), 1, fp);
    fread(&(x->el), sizeof(x->el), 1, fp);
    fread(&(x->no_l_indel), sizeof(x->no_l_indel), 1, fp);

    uint32_t t;
    fread(&(t), sizeof(t), 1, fp);
    x->ml = t;

    fread(&(t), sizeof(t), 1, fp);
    x->rev = t;

    
    fread(&(t), sizeof(t), 1, fp);
    x->bl = t;

    fread(&(t), sizeof(t), 1, fp);
    x->del = t;
}

int load_ma_hit_ts(ma_hit_t_alloc** x, char* read_file_name)
{
    fprintf(stdout, "Loading ma_hit_ts to disk ...... \n");
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "r");
    if(!fp)
    {
        return 0;
    }


    long long n_read;
    long long i, k;
    fread(&n_read, sizeof(n_read), 1, fp);
    (*x) = (ma_hit_t_alloc*)malloc(sizeof(ma_hit_t_alloc)*n_read);


    for (i = 0; i < n_read; i++)
    {        
        fread(&((*x)[i].is_fully_corrected), sizeof((*x)[i].is_fully_corrected), 1, fp);
        fread(&((*x)[i].is_abnormal), sizeof((*x)[i].is_abnormal), 1, fp);
        fread(&((*x)[i].length), sizeof((*x)[i].length), 1, fp);

        (*x)[i].buffer = (ma_hit_t*)malloc(sizeof(ma_hit_t)*(*x)[i].length);

        for (k = 0; k < (*x)[i].length; k++)
        {
            read_ma(&((*x)[i].buffer[k]), fp);
        }  
    }

    free(index_name);    
    fclose(fp);
    fprintf(stdout, "ma_hit_ts has been read.\n");
}




void write_ma(ma_hit_t* x, FILE* fp)
{
    fwrite(&(x->qns), sizeof(x->qns), 1, fp);
    fwrite(&(x->qe), sizeof(x->qe), 1, fp);
    fwrite(&(x->tn), sizeof(x->tn), 1, fp);
    fwrite(&(x->ts), sizeof(x->ts), 1, fp);
    fwrite(&(x->te), sizeof(x->te), 1, fp);
    fwrite(&(x->el), sizeof(x->el), 1, fp);
    fwrite(&(x->no_l_indel), sizeof(x->no_l_indel), 1, fp);

    uint32_t t = x->ml;
    fwrite(&(t), sizeof(t), 1, fp);
    t = x->rev;
    fwrite(&(t), sizeof(t), 1, fp);

    t = x->bl;
    fwrite(&(t), sizeof(t), 1, fp);
    t =x->del;
    fwrite(&(t), sizeof(t), 1, fp);
}


void write_ma_hit_ts(ma_hit_t_alloc* x, long long n_read, char* read_file_name)
{
    fprintf(stdout, "Writing ma_hit_ts to disk ...... \n");
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "w");
    long long i, k;
    fwrite(&n_read, sizeof(n_read), 1, fp);


    for (i = 0; i < n_read; i++)
    {
        fwrite(&(x[i].is_fully_corrected), sizeof(x[i].is_fully_corrected), 1, fp);
        fwrite(&(x[i].is_abnormal), sizeof(x[i].is_abnormal), 1, fp);
        fwrite(&(x[i].length), sizeof(x[i].length), 1, fp);
        for (k = 0; k < x[i].length; k++)
        {
            write_ma(x[i].buffer + k, fp);
        }  
    }


    free(index_name);
    fflush(fp);    
    fclose(fp);
    fprintf(stdout, "ma_hit_ts has been written.\n");
}

void write_all_data_to_disk(ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, 
All_reads *RNF, char* output_file_name)
{
    char* gfa_name = (char*)malloc(strlen(output_file_name)+25);
    sprintf(gfa_name, "%s.gfa.aux", output_file_name);
    write_All_reads(RNF, gfa_name);
    
    sprintf(gfa_name, "%s.gfa.aux.source", output_file_name);
    write_ma_hit_ts(sources, RNF->total_reads, gfa_name);

    sprintf(gfa_name, "%s.gfa.aux.reverse", output_file_name);
    write_ma_hit_ts(reverse_sources, RNF->total_reads, gfa_name);

    free(gfa_name);
}

int load_all_data_from_disk(ma_hit_t_alloc **sources, ma_hit_t_alloc **reverse_sources, 
char* output_file_name)
{
    char* gfa_name = (char*)malloc(strlen(output_file_name)+25);
    sprintf(gfa_name, "%s.gfa.aux", output_file_name);
    if(!load_All_reads(&R_INF, gfa_name))
    {
        return 0;
    }

    sprintf(gfa_name, "%s.gfa.aux.source", output_file_name);
    if(!load_ma_hit_ts(sources, gfa_name))
    {
        return 0;
    }


    sprintf(gfa_name, "%s.gfa.aux.reverse", output_file_name);
    if(!load_ma_hit_ts(reverse_sources, gfa_name))
    {
        return 0;
    }

    return 1;
}


// in a resolved bubble, mark unused vertices and arcs as "reduced"
static void asg_bub_backtrack(asg_t *g, uint32_t v0, buf_t *b)
{
	uint32_t i, v;
	///assert(b->S.n == 1);
	///first remove all nodes in this bubble
	for (i = 0; i < b->b.n; ++i)
		g->seq[b->b.a[i]>>1].del = 1;

	///second remove all edges (self/reverse for each edge) in this bubble
	for (i = 0; i < b->e.n; ++i) {
		asg_arc_t *a = &g->arc[b->e.a[i]];
		///remove this edge self
		a->del = 1;
		///remove the reverse direction
		asg_arc_del(g, a->v^1, a->ul>>32^1, 1);
	}
	///v is the sink of this bubble
	v = b->S.a[0];
	///recover node
	do {
		uint32_t u = b->a[v].p; // u->v
		g->seq[v>>1].del = 0;
		asg_arc_del(g, u, v, 0);
		asg_arc_del(g, v^1, u^1, 0);
		v = u;
	} while (v != v0);
}


// count the number of outgoing arcs, excluding reduced arcs
static inline int count_out(const asg_t *g, uint32_t v)
{
	uint32_t i, n, nv = asg_arc_n(g, v);
	const asg_arc_t *av = asg_arc_a(g, v);
	for (i = n = 0; i < nv; ++i)
		if (!av[i].del) ++n;
	return n;
}

// pop bubbles from vertex v0; the graph MJUST BE symmetric: if u->v present, v'->u' must be present as well
static uint64_t asg_bub_pop1(asg_t *g, uint32_t v0, int max_dist, buf_t *b)
{
	uint32_t i, n_pending = 0;
	uint64_t n_pop = 0;
	///if this node has been deleted
	if (g->seq[v0>>1].del) return 0; // already deleted
	///asg_arc_n(n0)
	if ((uint32_t)g->idx[v0] < 2) return 0; // no bubbles
	///S saves nodes with all incoming edges visited
	b->S.n = b->T.n = b->b.n = b->e.n = 0;
	///for each node, b->a saves all related information
	b->a[v0].c = b->a[v0].d = 0;
	///b->S is the nodes with all incoming edges visited
	kv_push(uint32_t, b->S, v0);

	do {
		///v is a node that all incoming edges have been visited
		///d is the distance from v0 to v
		uint32_t v = kv_pop(b->S), d = b->a[v].d, c = b->a[v].c;
		uint32_t nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		///why we have this assert?
		///assert(nv > 0);

		///all out-edges of v
		for (i = 0; i < nv; ++i) { // loop through v's neighbors
			/**
			p->ul: |____________31__________|__________1___________|______________32_____________|
								qn            direction of overlap       length of this node (not overlap length)
												(in the view of query)
			p->v : |___________31___________|__________1___________|
								tn             reverse direction of overlap
											(in the view of target)
			p->ol: overlap length
			**/
			uint32_t w = av[i].v, l = (uint32_t)av[i].ul; // v->w with length l
			binfo_t *t = &b->a[w];
			///that means there is a circle, directly terminate the whole bubble poping
			///if (w == v0) goto pop_reset;
            if ((w>>1) == (v0>>1)) goto pop_reset;
			///if this edge has been deleted
			if (av[i].del) continue;

			///push the edge
            ///high 32-bit of g->idx[v] is the start point of v's edges
            //so here is the point of this specfic edge
			kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);

			///find a too far path? directly terminate the whole bubble poping
			if (d + l > max_dist) break; // too far

            ///if this node
			if (t->s == 0) { // this vertex has never been visited
				kv_push(uint32_t, b->b, w); // save it for revert
				///t->p is the parent node of 
				///t->s = 1 means w has been visited
				///d is len(v0->v), l is len(v->w), so t->d is len(v0->w)
				t->p = v, t->s = 1, t->d = d + l;
				///incoming edges of w
				t->r = count_out(g, w^1);
				++n_pending;
			} else { // visited before
				///c is the weight (is very likely the number of node in this edge) of the parent node
				///select the longest edge (longest meams most reads/longest edge)
                if (c + 1 > t->c || (c + 1 == t->c && d + l > t->d)) t->p = v;
				if (c + 1 > t->c) t->c = c + 1;
				///update len(v0->w)
                ///node: t->d is not the length from this node's parent
                ///it is the shortest edge
				if (d + l < t->d) t->d = d + l; // update dist
			}
			///assert(t->r > 0);
			//if all incoming edges of w have visited
			//push it to b->S
			if (--(t->r) == 0) {
				uint32_t x = asg_arc_n(g, w);
				if (x) kv_push(uint32_t, b->S, w);
				///else kv_push(uint32_t, b->T, w); // a tip
                else goto pop_reset;
				--n_pending;
			}
		}
		///if i < nv, that means (d + l > max_dist)
		if (i < nv || b->S.n == 0) goto pop_reset;
	} while (b->S.n > 1 || n_pending);
    ///fprintf(stderr, "\nbeg>>1: %u, sink>>1: %u, num nodes: %d\n", v0>>1, (b->S.a[0])>>1, b->b.n);
	asg_bub_backtrack(g, v0, b);
	///n_pop = 1 | (uint64_t)b->T.n<<32;
    n_pop = 1;
    
    // for (i = 0; i < b->b.n; ++i)
    // {
    //     fprintf(stderr, "b->b.a[%d]: %d\n", i, b->b.a[i]>>1);
    // }
    
pop_reset:
	for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
		binfo_t *t = &b->a[b->b.a[i]];
		t->s = t->c = t->d = 0;
	}
	return n_pop;
}

// pop bubbles
int asg_pop_bubble(asg_t *g, int max_dist)
{
	uint32_t v, n_vtx = g->n_seq * 2;
	uint64_t n_pop = 0;
	buf_t b;
	if (!g->is_symm) asg_symm(g);
	memset(&b, 0, sizeof(buf_t));
	///set information for each node
	b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
	//traverse all node with two directions 
	for (v = 0; v < n_vtx; ++v) {
		uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		///some node could be deleted
		if (nv < 2 || g->seq[v>>1].del) continue;
		///some edges could be deleted
		for (i = 0; i < nv; ++i) // asg_bub_pop1() may delete some edges/arcs
			if (!av[i].del) ++n_arc;
		if (n_arc > 1)
			n_pop += asg_bub_pop1(g, v, max_dist, &b);
	}
	free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
	if (n_pop) asg_cleanup(g);
	///fprintf(stderr, "[M::%s] popped %d bubbles and trimmed %d tips\n", __func__, (uint32_t)n_pop, (uint32_t)(n_pop>>32));
	fprintf(stderr, "[M::%s] popped %d bubbles\n", __func__, n_pop);
    return n_pop;
}










int test_triangular_directly(asg_t *g, uint32_t v, 
long long min_edge_length, ma_hit_t_alloc* reverse_sources)
{
    
    uint32_t w;
    int todel;
    long long NodeLen_first[3];
    long long NodeLen_second[3];

    uint32_t Ns_first[3];
    uint32_t Ns_second[3];
    
    uint32_t nv = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v);
    if(av[0].v == av[1].v)
    {
        return 0;
    }
    /**********************test first node************************/
    NodeLen_first[0] = NodeLen_first[1] = NodeLen_first[2] = -1;
    if(asg_is_single_edge(g, av[0].v, v>>1) <= 2 && asg_is_single_edge(g, av[1].v, v>>1) <= 2)
    {
        NodeLen_first[asg_is_single_edge(g, av[0].v, v>>1)] = 0;
        NodeLen_first[asg_is_single_edge(g, av[1].v, v>>1)] = 1;
    }
    ///one node has one out-edge, another node has two out-edges
    if(NodeLen_first[1] == -1 || NodeLen_first[2] == -1)
    {
        return 0;
    }
    /**********************test first node************************/
    ///if the potiential edge has already been removed 
    if(av[NodeLen_first[2]].del == 1)
    {
        return 0;
    }

    /**********************test second node************************/
    w = av[NodeLen_first[2]].v^1;
    asg_arc_t *aw = asg_arc_a(g, w);
    uint32_t nw = asg_arc_n(g, w);
    if(nw != 2)
    {
        fprintf(stderr, "error\n");
    }
    NodeLen_second[0] = NodeLen_second[1] = NodeLen_second[2] = -1;
    if(asg_is_single_edge(g, aw[0].v, w>>1) <= 2 && asg_is_single_edge(g, aw[1].v, w>>1) <= 2)
    {
        NodeLen_second[asg_is_single_edge(g, aw[0].v, w>>1)] = 0;
        NodeLen_second[asg_is_single_edge(g, aw[1].v, w>>1)] = 1;
    }
    ///one node has one out-edge, another node has two out-edges
    if(NodeLen_second[1] == -1 || NodeLen_second[2] == -1)
    {
        return 0;
    }



    todel = 0;
    // if(check_if_diploid(av[0].v, av[1].v, g, reverse_sources, min_edge_length) && 
    // check_if_diploid(aw[0].v, aw[1].v, g, reverse_sources, min_edge_length))
    if(check_if_diploid(av[0].v, av[1].v, g, reverse_sources, min_edge_length) == 1||
    check_if_diploid(aw[0].v, aw[1].v, g, reverse_sources, min_edge_length) == 1)
    {
        todel = 1;
    }
    
    
    
    if(todel)
    {
        ///fprintf(stderr, "v: %u\n", v>>1);
        av[NodeLen_first[2]].del = 1;
        ///remove the reverse direction
        asg_arc_del(g, av[NodeLen_first[2]].v^1, av[NodeLen_first[2]].ul>>32^1, 1);
    }

    return todel;
}



int asg_arc_del_triangular_directly(asg_t *g, long long min_edge_length, ma_hit_t_alloc* reverse_sources)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, w, n_vtx = g->n_seq * 2, n_reduced = 0;


    int flag0, flag1, node;
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if (g->seq[v>>1].del)
        {
            continue;
        } 

        if(nv < 2)
        {
            continue;
        }


        n_reduced += test_triangular_directly(g, v, min_edge_length, reverse_sources);
    }

    

    if (n_reduced) {
        asg_cleanup(g);
        asg_symm(g);
    }

    fprintf(stderr, "[M::%s] removed %d triangular overlaps\n", 
    __func__, n_reduced);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

    return n_reduced;
}




int asg_arc_del_orthology(asg_t *g, ma_hit_t_alloc* reverse_sources, float drop_ratio, long long miniedgeLen)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, w, n_vtx = g->n_seq * 2, n_reduced = 0;
    uint32_t idx[2];

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (nv < 2 || g->seq[v>>1].del) continue;
        ///some edges could be deleted
        for (i = 0; i < nv; ++i) // asg_bub_pop1() may delete some edges/arcs
            if (!av[i].del) ++n_arc;
        if (n_arc != 2) continue;

        for (i = 0, n_arc = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                idx[n_arc] = i;
                n_arc++;
            }
        }

        if(check_if_diploid(av[idx[0]].v, av[idx[1]].v, g, reverse_sources, miniedgeLen) == 0)
        {
            float max = av[idx[0]].ol;
            float min = av[idx[1]].ol;

            if(min < drop_ratio * max)
            {
                av[idx[1]].del = 1;
                asg_arc_del(g, av[idx[1]].v^1, av[idx[1]].ul>>32^1, 1);
                n_reduced++;
            }
        }
    }


    
    

    if (n_reduced) {
        asg_cleanup(g);
        asg_symm(g);
    }

    fprintf(stderr, "[M::%s] removed %d different hap overlaps\n", 
    __func__, n_reduced);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

    return n_reduced;
}




int asg_arc_del_orthology_multiple_way(asg_t *g, ma_hit_t_alloc* reverse_sources, float drop_ratio, long long miniedgeLen)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, v_max, v_maxLen, w, n_vtx = g->n_seq * 2, n_reduced = 0;

    for (v = 0; v < n_vtx; ++v) 
    {
        if (g->seq_vis[v] != 0) continue;
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (nv < 2 || g->seq[v>>1].del) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;
        v_max = (uint32_t)-1;

        for (i = 0, n_arc = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                if(v_max == (uint32_t)-1)
                {
                    v_max = av[i].v;
                    v_maxLen = av[i].ol;
                }
                else if(check_if_diploid(v_max, av[i].v, g, reverse_sources, miniedgeLen) == 0)
                {
                    if(av[i].ol < drop_ratio * v_maxLen)
                    {
                        ///fprintf(stderr, "v: %u, v_max: %u, av[%d].v: %u\n", v>>1, v_max>>1, i, av[i].v>>1);
                        
                        av[i].ol = 1;
                        asg_arc_del(g, av[i].v^1, av[i].ul>>32^1, 1);
                        n_reduced++;
                    }
                }
            }
        }
    }

    if (n_reduced) {
        asg_cleanup(g);
        asg_symm(g);
    }

    fprintf(stderr, "[M::%s] removed %d different hap overlaps\n", 
    __func__, n_reduced);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

    return n_reduced;
}



int asg_arc_del_chimeric_read(asg_t *g, long long miniedgeLen)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, v_max, v_maxLen, w, n_vtx = g->n_seq * 2, n_reduced = 0, n_arc;

    for (v = 0; v < n_vtx; ++v) 
    {
        ///if (g->seq_vis[v] != 0) continue;
        if (g->seq[v>>1].del) continue;
        if (g->seq[v>>1].c == 0) continue;
        ///fprintf(stderr, "v>>1: %d\n", v>>1);
        if((get_real_length(g, v, NULL) == 0) || (get_real_length(g, v^1, NULL) == 0))
        {
            continue;
        }


        uint32_t convex1, convex2, flag1, flag2;
        long long l1, l2, ll;
        flag1 = detect_single_path_with_dels(g, v, &convex1, &l1, NULL);
        if(flag1 == END_TIPS || flag1 == LONG_TIPS) continue;
        flag2 = detect_single_path_with_dels(g, v^1, &convex2, &l2, NULL);
        if(flag2 == END_TIPS || flag2 == LONG_TIPS) continue;
        ll = l1 + l2 - 1;
        if(ll <= miniedgeLen)
        {
            ///fprintf(stderr, "***v>>1: %d\n", v>>1);
            asg_seq_del(g, v>>1);
            n_reduced++;
        }
    }

    if (n_reduced) {
        asg_cleanup(g);
        asg_symm(g);
    }

    fprintf(stderr, "[M::%s] removed %d chimeric reads\n", 
    __func__, n_reduced);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

    return n_reduced;
}


uint32_t detect_single_path_with_dels_by_length
(asg_t *g, uint32_t begNode, uint32_t* endNode, long long* Len, buf_t* b, long long maxLen)
{
    
    uint32_t v = begNode, w;
    uint32_t kv, kw;
    (*Len) = 0;


    while (1)
    {
        (*Len)++;
        kv = get_real_length(g, v, NULL);
        (*endNode) = v;

        if(b) kv_push(uint32_t, b->b, v>>1);

        if(kv == 0)
        {
            return END_TIPS;
        }

        if(kv == 2)
        {
            return TWO_OUTPUT;
        }

        if(kv > 2)
        {
            return MUL_OUTPUT;
        }

        if((*Len) > maxLen)
        {
            return LONG_TIPS;
        }

        ///up to here, kv=1
        ///kw must >= 1
        get_real_length(g, v, &w);
        kw = get_real_length(g, w^1, NULL);
        v = w;
        (*endNode) = v;


        if(kw == 2)
        {
            (*Len)++;
            if(b) kv_push(uint32_t, b->b, v>>1);
            return TWO_INPUT;
        }

        if(kw > 2)
        {
            (*Len)++;
            if(b) kv_push(uint32_t, b->b, v>>1);
            return MUL_INPUT;
        }   


        if((v>>1) == (begNode>>1))
        {
            return LOOP;
        } 
    }

    return LONG_TIPS;   
}



long long asg_arc_del_self_circle_untig(asg_t *g, long long circleLen)
{
    uint32_t v, v_max, w, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag;
    long long ll;
    asg_arc_t *aw;
    uint32_t nw, k;
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (g->seq[v>>1].del) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc != 1) continue;

        for (i = 0; i < nv; i++)
        {
            ///actually there is just one un-del edge
            if (!av[i].del)
            {    
                flag = detect_single_path_with_dels_by_length(g, v, &convex, &ll, NULL, circleLen);
                if(ll > circleLen || flag == LONG_TIPS)
                {
                    break;
                }
                if(flag == LOOP)
                {
                    break;
                }
                if(flag != END_TIPS && flag != LONG_TIPS)
                {
                    w = v^1;
                    n_arc = get_real_length(g, w, NULL);
                    if(n_arc == 0)
                    {
                        break;
                    }
                    aw = asg_arc_a(g, w);
                    nw = asg_arc_n(g, w);
                    for (k = 0; k < nw; k++)
                    {
                        if ((!aw[k].del) && (aw[k].v == (convex^1)))
                        {
                            // fprintf(stderr, "v: %u, w: %u, aw[%d].v: %u, convex: %u, ll: %d\n", 
                            // v>>1, w>>1, k, aw[k].v>>1, convex>>1, ll);

                            // fprintf(stderr, "*aw[k].v: %u, *convex: %u\n\n", 
                            // aw[k].v, convex);

                            aw[k].del = 1;
                            asg_arc_del(g, aw[k].v^1, aw[k].ul>>32^1, 1);
                            n_reduced++;
                        }
                    }
                }
            }
        }
    }

    if (n_reduced) {
        asg_cleanup(g);
        asg_symm(g);
    }

    fprintf(stderr, "[M::%s] removed %d self-circles\n", 
    __func__, n_reduced);

    return n_reduced;
}


void output_unitig_graph_without_small_bubbles(asg_t *sg, ma_sub_t* coverage_cut, 
char* output_file_name, long long n_read, long long bubble_dist, long long tipsLen)
{
    asg_cut_tip(sg, tipsLen);
    asg_pop_bubble(sg, bubble_dist);
    asg_cut_tip(sg, tipsLen);
    
    ma_ug_t *ug = NULL;
    ug = ma_ug_gen(sg);
    ma_ug_seq(ug, &R_INF, coverage_cut, n_read);

    fprintf(stdout, "Writing unitig GFA to disk ...... \n");
    char* gfa_name = (char*)malloc(strlen(output_file_name)+35);
    sprintf(gfa_name, "%s.no_s_bub.gfa", output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    ma_ug_print(ug, &R_INF, coverage_cut, output_file);
    fclose(output_file);
    
    
    sprintf(gfa_name, "%s.simple.no_s_bub.gfa", output_file_name);
    output_file = fopen(gfa_name, "w");
    ma_ug_print_simple(ug, &R_INF, coverage_cut, output_file);
    fclose(output_file);

    free(gfa_name);
    ma_ug_destroy(ug);
}


void output_contig_graph(asg_t *sg, ma_sub_t* coverage_cut, char* output_file_name, long long n_read, long long bubble_dist, long long tipsLen, float tip_drop_ratio, long long circleLen,
ma_hit_t_alloc* reverse_sources, long long miniedgeLen)
{
    asg_cut_tip(sg, tipsLen);
    // asg_pop_bubble(sg, bubble_dist);
    // asg_arc_del_self_circle_untig(sg, circleLen);
    long long n_ac = 1;
    long long pre_cons = sg->n_seq + sg->n_arc;
    long long cur_cons = 0;
    ///while(n_ac > 0)
    while(pre_cons != cur_cons)
    {
        pre_cons = sg->n_seq + sg->n_arc;
        n_ac = 0;
        n_ac += asg_pop_bubble(sg, bubble_dist);
        n_ac += asg_arc_del_self_circle_untig(sg, circleLen);
        n_ac += asg_arc_cut_long_tip(sg, tip_drop_ratio);
        n_ac += asg_arc_cut_long_equal_tips(sg, reverse_sources, 2);

        cur_cons = sg->n_seq + sg->n_arc;
    }

    asg_arc_identify_simple_bubbles_multi(sg, 1);
    asg_arc_del_short_false_link(sg, 0.6, 0.85, bubble_dist, reverse_sources, MAX_SHORT_TIPS);

    ///asg_arc_del_self_circle_untig(sg, circleLen);
    
    ma_ug_t *ug = NULL;
    ug = ma_ug_gen(sg);
    ma_ug_seq(ug, &R_INF, coverage_cut, n_read);

    fprintf(stdout, "Writing unitig GFA to disk ...... \n");
    char* gfa_name = (char*)malloc(strlen(output_file_name)+35);
    sprintf(gfa_name, "%s.contig.gfa", output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    ma_ug_print(ug, &R_INF, coverage_cut, output_file);
    fclose(output_file);
    
    



    sprintf(gfa_name, "%s.simple.contig.gfa", output_file_name);
    output_file = fopen(gfa_name, "w");
    ma_ug_print_simple(ug, &R_INF, coverage_cut, output_file);
    fclose(output_file);

    free(gfa_name);
    ma_ug_destroy(ug);
}


int output_tips(asg_t *g, const All_reads *RNF)
{
	uint32_t v, w, n_vtx = g->n_seq * 2, n_reduced = 0;
    for (v = 0; v < n_vtx; ++v)
    {
        if (g->seq[v>>1].del) continue;

        if(asg_arc_n(g, v) == 0)
        {
            fprintf(stderr, "%.*s\n",
             Get_NAME_LENGTH((*RNF), v>>1), 
             Get_NAME((*RNF), v>>1));
        }
    }
}

void build_string_graph_without_clean(
int min_dp, ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, 
long long n_read, uint64_t* readLen, 
long long mini_overlap_length, long long max_hang_length, long long clean_round, 
float min_ovlp_drop_ratio, float max_ovlp_drop_ratio, float corase_ovlp_drop_ratio, 
char* output_file_name, long long bubble_dist, int read_graph, int write)
{

    if (write_index_to_disk && write)
    {
        write_all_data_to_disk(sources, reverse_sources, 
        &R_INF, output_file_name);
    }


    // debug_info_of_specfic_read("m64016_190918_162737/179635219/ccs", 
    // sources, reverse_sources, -1, "init");

    // debug_info_of_specfic_read("m64016_190918_162737/130811282/ccs", 
    // sources, reverse_sources, -1, "init");
    

    
    
    ma_sub_t* coverage_cut;
    ///normalize_ma_hit_t(sources, n_read);
    normalize_ma_hit_t_single_side(sources, n_read);


    
    







    ///debug_normalize_ma_hit_t(sources, n_read);
    clean_weak_ma_hit_t(sources, reverse_sources, n_read);
    ///debug_normalize_ma_hit_t(sources, n_read);

    
    
    
    // debug_info_of_specfic_read("m64013_190322_203854/82051959/ccs", 
    // sources, reverse_sources, -1, "clean");

    // debug_info_of_specfic_read("m64013_190322_203854/74385680/ccs", 
    // sources, reverse_sources, -1, "clean");

    // debug_info_of_specfic_read("m64011_190329_072846/80545633/ccs", 
    // sources, reverse_sources, -1, "clean");


    
    

    ma_hit_sub(min_dp, sources, n_read, readLen, mini_overlap_length, &coverage_cut);
    ////这个会断开
    ma_hit_chimeric(1, sources, reverse_sources, n_read, readLen, coverage_cut);

    ma_hit_cut(min_dp, sources, n_read, readLen, mini_overlap_length, &coverage_cut);
    ///it seems we do not need ma_hit_flt
    ma_hit_flt(sources, n_read, coverage_cut, max_hang_length, mini_overlap_length);
    ///debug_cut_ma_hit_t(sources, n_read, coverage_cut);
    ma_hit_contained(sources, n_read, coverage_cut, max_hang_length, mini_overlap_length);
    ///debug_cut_ma_hit_t(sources, n_read, coverage_cut);

    

    // debug_info_of_specfic_read("m64016_190918_162737/179635219/ccs", 
    // sources, reverse_sources, -1, "contain");

    // debug_info_of_specfic_read("m64016_190918_162737/130811282/ccs", 
    // sources, reverse_sources, -1, "contain");
    
    // debug_info_of_specfic_read("m64016_190918_162737/72220752/ccs", 
    // sources, reverse_sources, -1, "contain");
    


    asg_t *sg = NULL;
    sg = ma_sg_gen(sources, n_read, coverage_cut, max_hang_length, mini_overlap_length);

    
    // debug_info_of_specfic_node("m64016_190918_162737/72220752/ccs", sg, "sg_gen");
    

    asg_arc_del_trans(sg, GAP_FUZZ);

    // debug_info_of_specfic_node("m64016_190918_162737/72220752/ccs", sg, "del_trans");

    char* unlean_name = (char*)malloc(strlen(output_file_name)+25);
    sprintf(unlean_name, "%s.unclean", output_file_name);
    output_read_graph(sg, coverage_cut, unlean_name, n_read);
    free(unlean_name);

    

    asg_cut_tip(sg, MAX_SHORT_TIPS);


    ///goto out;
    
    // debug_info_of_specfic_node("m64016_190918_162737/72220752/ccs", sg, "cut_tip");
    




    ///asg_arc_del_short_diploid_unclean(sg, corase_ovlp_drop_ratio, sources, reverse_sources);


    
    // debug_info_of_specfic_node("m64016_190918_162737/72220752/ccs", sg, "cut_corase");
    



    // asg_arc_del_single_node_bubble(sg, bubble_dist);
    // asg_cut_tip(sg, MAX_SHORT_TIPS);
    ///asg_cut_tip(sg, MAX_SHORT_TIPS);

    ///clean_round = 0;
    // fprintf(stderr, "\n\nWill perform %d round of clean...**********\n", 
    //         clean_round);

    
    if(clean_round > 0)
    {
        double cut_step;
        if(clean_round == 1)
        {
            cut_step = max_ovlp_drop_ratio;
        }
        else
        {
            cut_step = (max_ovlp_drop_ratio - min_ovlp_drop_ratio) / (clean_round - 1);
        }
        double drop_ratio = min_ovlp_drop_ratio;
        int i = 0;
        for (i = 0; i < clean_round; i++, drop_ratio += cut_step)
        {
            if(drop_ratio > max_ovlp_drop_ratio)
            {
                drop_ratio = max_ovlp_drop_ratio;
            }

            fprintf(stderr, "\n\n**********%d-th round drop: drop_ratio = %f**********\n", 
            i, drop_ratio);


            while(1)
            {
                int tri_flag = 0;
                tri_flag += asg_arc_del_self_circle_contig(sg);
                ///asg_arc_del_single_node_bubble(sg, bubble_dist);
                tri_flag += asg_arc_del_single_node_directly(sg, MAX_SHORT_TIPS, sources);
                tri_flag += asg_arc_del_triangular_advance(sg, bubble_dist);
                tri_flag += asg_arc_del_cross_bubble(sg, bubble_dist);
                ///asg_arc_del_single_node_bubble(sg, bubble_dist);
                tri_flag += asg_arc_del_single_node_directly(sg, MAX_SHORT_TIPS, sources);

                if(tri_flag == 0)
                {
                    break;
                }
            }


            ///asg_arc_del_orthology(sg, reverse_sources, drop_ratio, MAX_SHORT_TIPS);
            // asg_arc_del_orthology_multiple_way(sg, reverse_sources, drop_ratio, MAX_SHORT_TIPS);
            // asg_cut_tip(sg, MAX_SHORT_TIPS);
            
            /****************************may have bugs********************************/
            //asg_arc_identify_simple_bubbles(sg);
            asg_arc_identify_simple_bubbles_multi(sg, 1);
            //reomve edge between two chromesomes
            asg_arc_del_false_node(sg, MAX_SHORT_TIPS);
            asg_cut_tip(sg, MAX_SHORT_TIPS);
            /****************************may have bugs********************************/

            /****************************may have bugs********************************/
            ///asg_arc_identify_simple_bubbles_multi(sg, 1);
            asg_arc_identify_simple_bubbles_multi(sg, 0);
            ///asg_arc_del_short_diploid_unclean_exact(sg, drop_ratio, sources);
            asg_arc_del_short_diploid_by_exact(sg, MAX_SHORT_TIPS, sources);
            asg_cut_tip(sg, MAX_SHORT_TIPS);
            /****************************may have bugs********************************/


            //asg_arc_identify_simple_bubbles(sg);
            asg_arc_identify_simple_bubbles_multi(sg, 1);
            asg_arc_del_short_diploid_by_length(sg, drop_ratio, MAX_SHORT_TIPS, reverse_sources, MAX_SHORT_TIPS);
            asg_cut_tip(sg, MAX_SHORT_TIPS);

            asg_arc_identify_simple_bubbles_multi(sg, 1);
            asg_arc_del_short_false_link(sg, 0.6, 0.85, bubble_dist, reverse_sources, MAX_SHORT_TIPS);

            asg_arc_identify_simple_bubbles_multi(sg, 1);
            asg_arc_del_complex_false_link(sg, 0.6, 0.85, bubble_dist, reverse_sources, MAX_SHORT_TIPS);

            asg_cut_tip(sg, MAX_SHORT_TIPS);
        }
    }


    fprintf(stderr, "\n\n**********final clean**********\n");

    ///debug_info_of_specfic_node("m64016_190918_162737/72220752/ccs", sg, "before final clean");

    while(1)
    {
        int tri_flag = 0;
        tri_flag += asg_arc_del_self_circle_contig(sg);
        // fprintf(stderr, "tri_flag: %d\n", tri_flag);
        // fflush(stderr);
        ///asg_arc_del_single_node_bubble(sg, bubble_dist);
        tri_flag += asg_arc_del_single_node_directly(sg, MAX_SHORT_TIPS, sources);
        // fprintf(stderr, "tri_flag: %d\n", tri_flag);
        // fflush(stderr);
        tri_flag += asg_arc_del_triangular_advance(sg, bubble_dist);
        ///tri_flag += asg_arc_del_triangular_advance_debug(sg, bubble_dist);
        
        // fprintf(stderr, "tri_flag: %d\n", tri_flag);
        // fflush(stderr);
        tri_flag += asg_arc_del_cross_bubble(sg, bubble_dist);
        // fprintf(stderr, "tri_flag: %d\n", tri_flag);
        // fflush(stderr);
        ///asg_arc_del_single_node_bubble(sg, bubble_dist);
        tri_flag += asg_arc_del_single_node_directly(sg, MAX_SHORT_TIPS, sources);
        // fprintf(stderr, "tri_flag: %d\n", tri_flag);
        // fflush(stderr);

        if(tri_flag == 0)
        {
            break;
        }
    }

    

    asg_arc_del_short_diploi_by_suspect_edge(sg, MAX_SHORT_TIPS, sources);
    asg_cut_tip(sg, MAX_SHORT_TIPS);


    asg_arc_del_triangular_directly(sg, MAX_SHORT_TIPS, reverse_sources);





    ///asg_arc_identify_simple_bubbles_multi(sg, 0);
    // asg_arc_del_chimeric_read(sg, MAX_SHORT_TIPS*2);
    // asg_cut_tip(sg, MAX_SHORT_TIPS);
    



    
    asg_arc_identify_simple_bubbles_multi(sg, 0);
    asg_arc_del_orthology_multiple_way(sg, reverse_sources, 0.4, MAX_SHORT_TIPS);
    asg_cut_tip(sg, MAX_SHORT_TIPS);

    



    asg_arc_identify_simple_bubbles_multi(sg, 0);
    asg_arc_del_too_short_overlaps(sg, 2000, min_ovlp_drop_ratio, reverse_sources, MAX_SHORT_TIPS);
    asg_cut_tip(sg, MAX_SHORT_TIPS);


    /**
    asg_arc_identify_simple_bubbles_multi(sg, 1);
    asg_arc_del_short_false_link_advance(sg, 0.6, 0.85, bubble_dist, reverse_sources, MAX_SHORT_TIPS);
    **/
    
    



    ///asg_arc_del_triangular_advance_debug(sg, bubble_dist);

    /**
    fprintf(stderr, "\n\n**********final aggressive clean**********\n");
    while(1)
    {
        int tri_flag = 0;
        asg_arc_identify_simple_bubbles_multi(sg, 0);
        tri_flag = asg_arc_del_tri_link(sg, bubble_dist);
        if(tri_flag == 0)
        {
            break;
        }
    }
    **/



    


    /****************************may have bugs********************************/
    /**
    long long c_tips = 1;
    int i = 0;
    while (c_tips && i < clean_round)
    {
        asg_arc_identify_simple_bubbles_multi(sg, 0);
        c_tips = asg_arc_del_short_false_link(sg, 0.7, bubble_dist);

        asg_arc_identify_simple_bubbles_multi(sg, 0);
        c_tips += asg_arc_del_complex_false_link(sg, 0.7, bubble_dist);

        if(c_tips) asg_cut_tip(sg, MAX_SHORT_TIPS);
        i++;
    }
    **/
    
    /****************************may have bugs********************************/

    /**
    memset(sg->seq_vis, 0, sg->n_seq*2*sizeof(uint8_t));
    asg_arc_del_short_diploid_by_exact(sg, MAX_SHORT_TIPS, sources);
    asg_cut_tip(sg, MAX_SHORT_TIPS);
    **/
    
    // debug_info_of_specfic_node("m64016_190918_162737/141297762/ccs", sg);
    
    out:
    ///output_tips(sg, &R_INF);


    output_unitig_graph(sg, coverage_cut, output_file_name, n_read);
    output_read_graph(sg, coverage_cut, output_file_name, n_read);

    /****************************may have bugs********************************/
    output_unitig_graph_without_small_bubbles(sg, coverage_cut, output_file_name, n_read, 
    100000, MAX_SHORT_TIPS);
    /****************************may have bugs********************************/

    output_contig_graph(sg, coverage_cut, output_file_name, n_read, 10000000, MAX_SHORT_TIPS, 0.1, 20,
    reverse_sources, MAX_SHORT_TIPS);
    ///output_contig_graph(sg, coverage_cut, output_file_name, n_read, 10000000);

    
    asg_destroy(sg);
    free(coverage_cut);
}
