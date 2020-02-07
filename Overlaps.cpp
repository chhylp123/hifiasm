#include <stdio.h>
#include <stdlib.h>
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include "Overlaps.h"
#include "ksort.h"
#include "Process_Read.h"
#include "CommandLines.h"
#include "Hash_Table.h"
#include "Correct.h"

KDQ_INIT(uint64_t)
KDQ_INIT(uint32_t)

#define ma_hit_key_tn(a) ((a).tn)
KRADIX_SORT_INIT(hit_tn, ma_hit_t, ma_hit_key_tn, member_size(ma_hit_t, tn))

#define ma_hit_key_qns(a) ((a).qns)
KRADIX_SORT_INIT(hit_qns, ma_hit_t, ma_hit_key_qns, member_size(ma_hit_t, qns))

#define asg_arc_key(a) ((a).ul)
KRADIX_SORT_INIT(asg, asg_arc_t, asg_arc_key, 8)

#define generic_key(x) (x)
KRADIX_SORT_INIT(arch64, uint64_t, generic_key, 8)


///#define Hap_Align_key(a) ((((uint64_t)((a).is_color))<<63)|(((uint64_t)((a).t_id))<<32)|((uint64_t)((a).q_pos)))
///#define Hap_Align_key(a) ((((uint64_t)((a).is_color))<<32)|(((uint64_t)((a).t_id))<<33)|((uint64_t)((a).q_pos)))
#define Hap_Align_key(a) ((((uint64_t)((a).t_id))<<33)|((uint64_t)((a).q_pos)))
KRADIX_SORT_INIT(Hap_Align_sort, Hap_Align, Hap_Align_key, 8)


KSORT_INIT_GENERIC(uint32_t)

///this value has been updated at the first line of build_string_graph_without_clean
long long min_thres;

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


void add_overlaps(ma_hit_t_alloc* source_paf, ma_hit_t_alloc* dest_paf, uint64_t* source_index, long long listLen)
{
    long long i;
    ma_hit_t* tmp;
    for (i = 0; i < listLen; i++)
    {
        tmp = &(source_paf->buffer[(uint32_t)(source_index[i])]);
        add_ma_hit_t_alloc(dest_paf, tmp);
    }
}


void remove_overlaps(ma_hit_t_alloc* source_paf, uint64_t* source_index, long long listLen)
{
    long long i, m;
    for (i = 0; i < listLen; i++)
    {
        source_paf->buffer[(uint32_t)(source_index[i])].qns = (uint64_t)(-1);
    }

    m = 0;
    for (i = 0; i < source_paf->length; i++)
    {
        if(source_paf->buffer[i].qns != (uint64_t)(-1))
        {
            source_paf->buffer[m] = source_paf->buffer[i];
            m++;
        }
    }
    source_paf->length = m;
}


void add_overlaps_from_different_sources(ma_hit_t_alloc* source_paf_list, ma_hit_t_alloc* dest_paf, 
uint64_t* source_index, long long listLen)
{
    long long i;
    ma_hit_t ele;
    ma_hit_t* tmp;
    uint32_t source_n, source_i;
    for (i = 0; i < listLen; i++)
    {
        source_n = source_index[i] >> 32;
        source_i = (uint32_t)(source_index[i]);
        tmp = &(source_paf_list[source_n].buffer[source_i]);

        ele.del = 0;
        ele.rev = tmp->rev;
        ele.qns = Get_tn((*tmp));
        ele.qns = ele.qns << 32;
        ele.qns = ele.qns | (uint64_t)(Get_ts((*tmp)));
        ele.qe = Get_te((*tmp));

        ele.tn = Get_qn((*tmp));
        ele.ts = Get_qs((*tmp));
        ele.te = Get_qe((*tmp));

        ele.bl = R_INF.read_length[ele.tn];
        ele.ml = tmp->ml;
        ele.el = tmp->el;
        ele.no_l_indel = tmp->no_l_indel;
        
        add_ma_hit_t_alloc(dest_paf, &ele);
    }
}

void print_revise_edges(ma_hit_t_alloc* source_paf, uint64_t* source_index, long long listLen)
{
    long long i;
    ma_hit_t* tmp;
    for (i = 0; i < listLen; i++)
    {
        tmp = &(source_paf->buffer[(uint32_t)(source_index[i])]);

        fprintf(stderr, "%.*s(%d) ---(+)--> %.*s(%d), Len: %d\n", 
        (int)Get_NAME_LENGTH(R_INF, Get_qn((*tmp))), Get_NAME(R_INF, Get_qn((*tmp))), Get_qn((*tmp)),
        (int)Get_NAME_LENGTH(R_INF, Get_tn((*tmp))), Get_NAME(R_INF, Get_tn((*tmp))), Get_tn((*tmp)),
        Get_qe((*tmp)) - Get_qs((*tmp)));
    }
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
	if (sid >= (int)g->m_seq) {
		g->m_seq = sid + 1;
		kv_roundup32(g->m_seq);
		g->seq = (asg_seq_t*)realloc(g->seq, g->m_seq * sizeof(asg_seq_t));
	}

	if (sid >= g->n_seq) g->n_seq = sid + 1;
	
    g->seq[sid].del = !!del;
    g->seq[sid].len = len;
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
    ///remove edges with del, and free idx
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

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d multi-arcs\n", __func__, n_multi);
    }
	
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
    if(VERBOSE >= 1)
    {
	    fprintf(stderr, "[M::%s] removed %d asymmetric arcs\n", __func__, n_asymm);
    }
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
    /**
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
    **/
    dest->ml = source->ml;
    dest->no_l_indel = source->no_l_indel;
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
                new_element.del = 0;
                new_element.ml = 1;
                new_element.no_l_indel = 1;
                set_reverse_overlap(&new_element, &(sources[i].buffer[j]));
                add_ma_hit_t_alloc(&(sources[tn]), &new_element);
                si_overlaps++;
            }
        }
    } 
}


void normalize_ma_hit_t_single_side(ma_hit_t_alloc* sources, long long num_sources)
{
    double startTime = Get_T();

    long long bi_overlaps = 0;
    long long i, j, index;
    uint32_t qn, tn;
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

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2fs\n\n", __func__, Get_T()-startTime);
    }
}




void ma_hit_contained(ma_hit_t_alloc* sources, long long n_read, ma_sub_t *coverage_cut, 
R_to_U* ruIndex, int max_hang, int min_ovlp)
{
    double startTime = Get_T();
	int32_t r;
	long long i, j, m;
	asg_arc_t t;
    for (i = 0; i < n_read; ++i) 
    {
        for (j = 0; j < (long long)sources[i].length; j++)
        {
            ma_hit_t *h = &(sources[i].buffer[j]);
            //check the corresponding two reads 
		    ma_sub_t *sq = &(coverage_cut[Get_qn(*h)]);
            ma_sub_t *st = &(coverage_cut[Get_tn(*h)]);
            r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
            ///r could not be MA_HT_SHORT_OVLP or MA_HT_INT
            if (r == MA_HT_QCONT) 
            {
                if(st->del == 1) continue;
                sq->del = 1;
                set_R_to_U(ruIndex, Get_qn(*h), Get_tn(*h), 0);
            }
		    else if (r == MA_HT_TCONT) 
            {
                if(sq->del == 1) continue;
                st->del = 1;
                set_R_to_U(ruIndex, Get_tn(*h), Get_qn(*h), 0);
            }
        }
    }

    transfor_R_to_U(ruIndex);

    for (i = 0; i < n_read; ++i) 
    {
        m = 0;
        for (j = 0; j < (long long)sources[i].length; j++)
        {
            ma_hit_t *h = &(sources[i].buffer[j]);
            ///both the qn and tn have not been deleted
            if(coverage_cut[Get_qn(*h)].del != 1 && coverage_cut[Get_tn(*h)].del != 1)
            {
                h->del = 0;
                m++;
            }
            else
            {
                h->del = 1;
            }
        }

        ///may have bugs here 
        ///if sources[i].length == 0, that means all overlapped reads with read i are the contained reads
        if(m == 0)
        {
            coverage_cut[i].del = 1;
        }
    }

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

}



void ma_hit_flt(ma_hit_t_alloc* sources, long long n_read, ma_sub_t *coverage_cut, int max_hang, int min_ovlp)
{
    double startTime = Get_T();
	long long i, j, m;
	asg_arc_t t;

    for (i = 0; i < n_read; ++i) 
    {
        m = 0;
        for (j = 0; j < (long long)sources[i].length; j++)
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
            ///here the max_hang = 1000, asm_opt.max_hang_rate = 0.8, min_ovlp = 50
            ///for me, there should not have any overhang..so r cannot be equal to MA_HT_INT
            ///sq->e - sq->s = the length of query; st->e - st->s = the length od target
            r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);


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

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
}





///a is the overlap vector, n is the length of overlap vector
///min_dp is used for coverage droping
///select reads with coverage >= min_dp
void ma_hit_sub(int min_dp, ma_hit_t_alloc* sources, long long n_read, uint64_t* readLen, 
long long mini_overlap_length, ma_sub_t** coverage_cut)
{
    double startTime = Get_T();

    (*coverage_cut) = (ma_sub_t*)malloc(sizeof(ma_sub_t)*n_read);

	uint64_t i, j, n_remained = 0;
	kvec_t(uint32_t) b = {0,0,0};

	///all overlaps in vector a has been sorted by qns
	///so for overlaps of one reads, it must be contiguous
	for (i = 0; i < (uint64_t)n_read; ++i) 
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
        int dp, start = 0;
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
                if (len > (int)(max.e - max.s)) 
                {
                    max2 = max; 
                    max.s = start;
                    max.e = b.a[j]>>1;
                }
                else if (len > int(max2.e - max2.s)) 
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
    if(VERBOSE >= 1)
    {   
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
}




int boundary_verify(uint32_t x_interval_s, uint32_t x_interval_e, ma_hit_t* map, 
char* x_buffer, char* y_buffer, All_reads* R_INF)
{
    uint32_t xs, ys, dir, x_id, y_id, x_interval_Len, y_interval_Len, y_interval_s, y_interval_e;
    dir = (*map).rev;
    xs = Get_qs((*map));
    x_id = Get_qn((*map));
    y_id = Get_tn((*map));
    long long yLen = Get_READ_LENGTH((*R_INF), y_id);

    if(dir == 1)
    {
        ys = yLen - (Get_te((*map)) - 1) - 1;         
    }
    else
    {
        ys = Get_ts((*map));
    }
    ///[x_interval_s, x_interval_e)
    x_interval_Len = x_interval_e - x_interval_s;

    ///[y_interval_s, y_interval_e]
    y_interval_s = (x_interval_s - xs) + ys;
    if(y_interval_s >= yLen)
    {
        return 0;
    }
    y_interval_e = y_interval_s + x_interval_Len - 1;
    if(y_interval_e >= yLen)
    {
        y_interval_e = yLen - 1;
    }

    if(y_interval_e < y_interval_s)
    {
        return 0;
    }
    
    y_interval_Len = y_interval_e - y_interval_s + 1;


    if(y_interval_Len <= WINDOW)
    {
        return verify_single_window(y_interval_s, y_interval_e, ys, xs, y_id, x_id, 
                                                        dir, y_buffer, x_buffer, R_INF);
    }
    else
    {
        
        if(verify_single_window(y_interval_s, y_interval_s + WINDOW - 1, ys, xs, y_id, x_id, 
                                                        dir, y_buffer, x_buffer, R_INF) == 0)
        {
            return 0;
        }

        if(verify_single_window(y_interval_e - WINDOW + 1, y_interval_e, ys, xs, y_id, x_id, 
                                                        dir, y_buffer, x_buffer, R_INF) == 0)
        {
            return 0;
        }

        return 1;
    }
}

void collect_sides(ma_hit_t_alloc* paf, uint64_t rLen, ma_sub_t* max_left, ma_sub_t* max_right)
{
    long long j;
    uint32_t qs, qe;
    for (j = 0; j < paf->length; j++)
    {
        qs = Get_qs(paf->buffer[j]);
        qe = Get_qe(paf->buffer[j]);


        ///overlaps from left side
        if(qs == 0)
        {
            if(qs < max_left->s) max_left->s = qs;
            if(qe > max_left->e) max_left->e = qe;
        }

        ///overlaps from right side
        if(qe == rLen)
        {
            if(qs < max_right->s) max_right->s = qs;
            if(qe > max_right->e) max_right->e = qe;
        }

        ///note: if (qs == 0 && qe == rLen)
        ///this overlap would be added to both b_left and b_right
        ///that is what we want
    }
}

void collect_contain(ma_hit_t_alloc* paf1, ma_hit_t_alloc* paf2, uint64_t rLen, 
ma_sub_t* max_left, ma_sub_t* max_right, float overlap_rate)
{
    long long j, new_left_e, new_right_s;
    new_left_e = max_left->e;
    new_right_s = max_right->s;
    uint32_t qs, qe;
    ma_hit_t_alloc* paf;

    if(paf1 != NULL)
    {
        paf = paf1;
        for (j = 0; j < paf->length; j++)
        {
            qs = Get_qs(paf->buffer[j]);
            qe = Get_qe(paf->buffer[j]);
            ///check contained overlaps
            if(qs != 0 && qe != rLen)
            {
                ///[qs, qe), [max_left.s, max_left.e)
                if(qs < max_left->e && qe > max_left->e && max_left->e - qs > (overlap_rate * (qe -qs)))
                {
                    ///if(qe > max_left->e) max_left->e = qe;
                    if(qe > max_left->e && qe > new_left_e) new_left_e = qe;
                }

                ///[qs, qe), [max_right.s, max_right.e)
                if(qs < max_right->s && qe > max_right->s && qe - max_right->s > (overlap_rate * (qe -qs)))
                {
                    ///if(qs < max_right->s) max_right->s = qs;
                    if(qs < max_right->s && qs < new_right_s) new_right_s = qs;
                }
            }
        }
    }

    if(paf2 != NULL)
    {
        paf = paf2;
        for (j = 0; j < paf->length; j++)
        {
            qs = Get_qs(paf->buffer[j]);
            qe = Get_qe(paf->buffer[j]);
            ///check contained overlaps
            if(qs != 0 && qe != rLen)
            {
                ///[qs, qe), [max_left.s, max_left.e)
                if(qs < max_left->e && qe > max_left->e && max_left->e - qs > (overlap_rate * (qe -qs)))
                {
                    ///if(qe > max_left->e) max_left->e = qe;
                    if(qe > max_left->e && qe > new_left_e) new_left_e = qe;
                }

                ///[qs, qe), [max_right.s, max_right.e)
                if(qs < max_right->s && qe > max_right->s && qe - max_right->s > (overlap_rate * (qe -qs)))
                {
                    ///if(qs < max_right->s) max_right->s = qs;
                    if(qs < max_right->s && qs < new_right_s) new_right_s = qs;
                }
            }
        }
    }
    

    max_left->e = new_left_e;
    max_right->s = new_right_s;
}



int intersection_check(ma_hit_t_alloc* paf, uint64_t rLen, uint32_t interval_s, uint32_t interval_e)
{
    long long j, cov = 0;
    uint32_t qs, qe;

    for (j = 0; j < paf->length; j++)
    {
        qs = Get_qs(paf->buffer[j]);
        qe = Get_qe(paf->buffer[j]);
        ///[interval_s, interval_e) must be at least contained at one of the [qs, qe)
        if(qs<=interval_s && qe>=interval_e)
        {
            cov++;
        }
    }

    return cov;
}


int intersection_check_by_base(ma_hit_t_alloc* paf, uint64_t rLen, uint32_t interval_s, uint32_t interval_e,
char* bq, char* bt)
{
    long long j;
    uint32_t qs, qe;

    for (j = 0; j < paf->length; j++)
    {
        qs = Get_qs(paf->buffer[j]);
        qe = Get_qe(paf->buffer[j]);
        ///[interval_s, interval_e) must be at least contained at one of the [qs, qe)
        if(qs<=interval_s && qe>=interval_e)
        {
            if(boundary_verify(interval_s, interval_e, &(paf->buffer[j]), bq, bt, &R_INF) == 0)
            {
                return 1;
            }
        }
    }

    return 0;
}

void print_overlaps(ma_hit_t_alloc* paf, long long rLen, long long interval_s, long long interval_e)
{
    long long j;

    fprintf(stderr, "left: \n");
    for (j = 0; j < paf->length; j++)
    {
        if(Get_qs(paf->buffer[j]) == 0)
        {
            fprintf(stderr, "?????? interval_s: %lld, interval_e: %lld, qn: %u, tn: %u, j: %lld, qs: %u, qe: %u, ts: %u, te: %u, dir: %u\n",
            interval_s, interval_e,
            Get_qn(paf->buffer[j]), Get_tn(paf->buffer[j]),
            j, Get_qs(paf->buffer[j]), Get_qe(paf->buffer[j]),
            Get_ts(paf->buffer[j]), Get_te(paf->buffer[j]),
            paf->buffer[j].rev);

            fprintf(stderr, "%.*s\n", (int)Get_NAME_LENGTH(R_INF, Get_tn(paf->buffer[j])),
                Get_NAME(R_INF, Get_tn(paf->buffer[j])));
        }   
    }

    fprintf(stderr, "right: \n");
    for (j = 0; j < paf->length; j++)
    {
        if(Get_qe(paf->buffer[j]) == rLen)
        {
            fprintf(stderr, "?????? interval_s: %lld, interval_e: %lld, qn: %u, tn: %u, j: %lld, qs: %u, qe: %u, ts: %u, te: %u, dir: %u\n",
            interval_s, interval_e,
            Get_qn(paf->buffer[j]), Get_tn(paf->buffer[j]),
            j, Get_qs(paf->buffer[j]), Get_qe(paf->buffer[j]),
            Get_ts(paf->buffer[j]), Get_te(paf->buffer[j]),
            paf->buffer[j].rev);
            
            fprintf(stderr, "%.*s\n", (int)Get_NAME_LENGTH(R_INF, Get_tn(paf->buffer[j])),
                Get_NAME(R_INF, Get_tn(paf->buffer[j])));
        }
    }


    fprintf(stderr, "middle: \n");
    for (j = 0; j < paf->length; j++)
    {
        if(Get_qs(paf->buffer[j]) != 0 && Get_qe(paf->buffer[j]) != rLen)
        {
            fprintf(stderr, "?????? interval_s: %lld, interval_e: %lld, qn: %u, tn: %u, j: %lld, qs: %u, qe: %u, ts: %u, te: %u, dir: %u\n",
            interval_s, interval_e,
            Get_qn(paf->buffer[j]), Get_tn(paf->buffer[j]),
            j, Get_qs(paf->buffer[j]), Get_qe(paf->buffer[j]),
            Get_ts(paf->buffer[j]), Get_te(paf->buffer[j]),
            paf->buffer[j].rev);
            
            fprintf(stderr, "%.*s\n", (int)Get_NAME_LENGTH(R_INF, Get_tn(paf->buffer[j])),
                Get_NAME(R_INF, Get_tn(paf->buffer[j])));
        }
    }

}


void detect_chimeric_reads(ma_hit_t_alloc* paf, ma_hit_t_alloc* rev_paf_x, 
long long n_read, uint64_t* readLen, ma_sub_t* coverage_cut, float shift_rate)
{
    double startTime = Get_T();
    init_aux_table();
    long long i, rLen, /**cov,**/ n_simple_remove = 0, n_complex_remove = 0, n_complex_remove_real = 0;
    uint32_t interval_s, interval_e;
    ma_sub_t max_left, max_right;
    kvec_t(char) b_q = {0,0,0};
    kvec_t(char) b_t = {0,0,0};
    for (i = 0; i < n_read; ++i) 
    {
        coverage_cut[i].c = PRIMARY_LABLE;
        rLen = readLen[i];


        max_left.s = max_right.s = rLen;
        max_left.e = max_right.e = 0;

        // if(i == 6937245)
        // {
        //     fprintf(stderr, "\n\npaf: \n");
        //     print_overlaps(&paf[i], rLen, interval_s, interval_e);
        //     fprintf(stderr, "\n\nrev_paf: \n");
        //     print_overlaps(&rev_paf[i], rLen, interval_s, interval_e);
        // }


        collect_sides(&(paf[i]), rLen, &max_left, &max_right);
        ///collect_sides(&(rev_paf[i]), rLen, &max_left, &max_right);
        ///that means this read is an end node
        if(max_left.s == rLen || max_right.s == rLen)
        {
            continue;
        }


        // if(i == 6937245)
        // {
        //     fprintf(stderr, "max_left.s: %d, max_left.e: %d, max_right.s: %d, max_right.e: %d\n", 
        //     max_left.s, max_left.e, max_right.s, max_right.e);
        // }


        collect_contain(&(paf[i]), NULL, rLen, &max_left, &max_right, 0.1);
        ///collect_contain(&(paf[i]), &(rev_paf[i]), rLen, &max_left, &max_right, 0.1);

        // if(i == 6937245)
        // {
        //     fprintf(stderr, "max_left.s: %d, max_left.e: %d, max_right.s: %d, max_right.e: %d\n", 
        //     max_left.s, max_left.e, max_right.s, max_right.e);
        // }


        ////shift_rate should be (FINAL_OVERLAP_ERROR_RATE*2)
        ///this read is a normal read
        if(max_left.e > max_right.s && 
           (max_left.e - max_right.s >= rLen * shift_rate))
        {
            continue;
        }

        ///simple chimeric reads
        if(max_left.e <= max_right.s)
        {
            coverage_cut[i].del = 1;
            paf[i].length = 0;
            n_simple_remove++;
            continue;
        }

        // if(i == 6937245)
        // {
        //     fprintf(stderr, "max_left.s: %d, max_left.e: %d, max_right.s: %d, max_right.e: %d\n", 
        //     max_left.s, max_left.e, max_right.s, max_right.e);
        // }

        ///now max_left.e >  max_right.s && max_left.e - max_right.s is small enough
        //[interval_s, interval_e)
        interval_s = max_right.s;
        interval_e = max_left.e;

        /**
        cov = 0;
        cov += intersection_check(&(paf[i]), rLen, interval_s, interval_e);
        cov += intersection_check(&(rev_paf[i]), rLen, interval_s, interval_e);
        if(interval_e - interval_s < WINDOW && cov <= 2)
        {
            coverage_cut[i].del = 1;
            paf[i].length = 0;
            n_complex_remove_real++;
        }
        else**/
        {
            kv_resize(char, b_q, WINDOW*4+20);
            kv_resize(char, b_t, WINDOW*4+20);
            if(intersection_check_by_base(&(paf[i]), rLen, interval_s, interval_e, b_q.a, b_t.a)
              /**||
              intersection_check_by_base(&(rev_paf[i]), rLen, interval_s, interval_e, b_q.a, b_t.a)**/)
            {
                coverage_cut[i].del = 1;
                paf[i].length = 0;
                n_complex_remove_real++;
            }
        }

        n_complex_remove++;
    }

    free(b_q.a);
    free(b_t.a);

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2f s, n_simple_remove: %lld, n_complex_remove: %lld/%lld\n\n", 
        __func__, Get_T()-startTime, n_simple_remove, n_complex_remove_real, n_complex_remove);
    }
}





void ma_hit_cut(ma_hit_t_alloc* sources, long long n_read, uint64_t* readLen, 
long long mini_overlap_length, ma_sub_t** coverage_cut)
{
    double startTime = Get_T();
	size_t i, j;
    ma_hit_t* p;
    ma_sub_t* rq;
    ma_sub_t* rt;
    long long m = 0;
    for (i = 0; i < (uint64_t)n_read; ++i) 
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
            qs = ((uint32_t)qs > rq->s? qs : rq->s) - rq->s;
            qe = ((uint32_t)qe < rq->e? qe : rq->e) - rq->s;
            ts = ((uint32_t)ts > rt->s? ts : rt->s) - rt->s;
            te = ((uint32_t)te < rt->e? te : rt->e) - rt->s;

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

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    
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
	for (i = nv = 0; i < (int)nv0; ++i)
		if (!av[i].del) i0 = i, ++nv;

	///end without any out-degree
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
	for (i = nw = 0; i < (int)nw0; ++i)
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
	uint32_t nv, nv0 = asg_arc_n(g, v^1);
	int i;
	asg_arc_t *av = asg_arc_a(g, v^1);

    int flag = 0;
	///if this arc has not been deleted
	for (i = nv = 0; i < (int)nv0; ++i)
    {
        ///if (!av[i].del)
        {
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
	for (i = 0; i < (uint64_t)n_read; ++i) 
    {
        ///if a read has been deleted, should we still add them?
		asg_seq_set(g, i, coverage_cut[i].e - coverage_cut[i].s, coverage_cut[i].del);
        g->seq[i].c = coverage_cut[i].c;
	}

    g->seq_vis = (uint8_t*)calloc(g->n_seq*2, sizeof(uint8_t));

    for (i = 0; i < (uint64_t)n_read; ++i) 
    {
        for (j = 0; j < sources[i].length; j++)
        {
            int r;
		    asg_arc_t t, *p;
		    const ma_hit_t *h = &(sources[i].buffer[j]);
            if(h->del) continue;

            //high coverage region [sub[qn].e, sub[qn].s) in query
            int ql = coverage_cut[Get_qn(*h)].e - coverage_cut[Get_qn(*h)].s;
            //high coverage region [sub[qn].e, sub[qn].s) in target
            int tl = coverage_cut[Get_tn(*h)].e - coverage_cut[Get_tn(*h)].s;
		    r = ma_hit2arc(h, ql, tl, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
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
    g->r_seq = g->n_seq;

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

	return g;
}


// pop bubbles from vertex v0; the graph MJUST BE symmetric: if u->v present, v'->u' must be present as well
//note!!!!!!!! here we don't exculde the deleted edges
uint64_t asg_bub_finder_with_del(asg_t *g, uint32_t v0, int max_dist, buf_t *b,
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
			if (d + l > (uint32_t)max_dist)
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
			if (d + l > (uint32_t)max_dist)
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
			if (d + l > (uint32_t)max_dist)
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
			if (d + l > (uint32_t)max_dist)
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
    long long pre_b_n = 0;
    (*Len) = 0;

    
    
    while (1)
    {
        nv = asg_arc_n(g, v);
        av = asg_arc_a(g, v);   
        (*endNode) = v;
        (*Len)++;


        if((*Len) > max_ext)
        {
            return LONG_TIPS_UNDER_MAX_EXT;
        }

        
        if(b) kv_push(uint32_t, b->b, v>>1);

 

        if(nv == 0)
        {
            return END_TIPS;
        }

        if(nv == 2)
        {


            if(b) pre_b_n = b->b.n;
            if(!detect_simple_bubble(g, v, &v, &bLen, b))
            {
                if(b) b->b.n = pre_b_n;
                return TWO_OUTPUT;
            }

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

               uint32_t convex1 = 0, convex2 = 0, f1, f2;
               long long l1 = 0, l2 = 0;
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
               }
 
        }

    }

    return n_reduced;
}





int find_single_link(asg_t *g, uint32_t link_beg, int linkLen, uint32_t* link_end)
{
    uint32_t v, w;
    v = link_beg^1;
    uint32_t nv, nw;
    asg_arc_t *av;
    int edgeLen = 0;

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
    uint32_t Nodes1[2]={0};
    uint32_t Nodes2[2]={0};

    uint32_t Ns_first[2]={0};
    uint32_t Ns_second[2]={0};


    
    for (i = 0; i < length; ++i)
    {
        v = nodes[i];

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

        uint32_t convex1, convex2, f1, f2;
        long long l1, l2;
        todel = 0;

        if((Nodes1[0]^1) == (startNode^1) || (Nodes1[0]^1) == endNode)
        {
            continue;
        }
        if((Nodes2[1]^1) == (startNode^1) || (Nodes2[1]^1) == endNode)
        {
            continue;
        }

        f1 = detect_bubble_end_with_bubbles(g, Nodes1[0]^1, Nodes2[1]^1, &convex1, &l1, NULL);


        if((Nodes2[0]^1) == (startNode^1) || (Nodes2[0]^1) == endNode)
        {
            continue;
        }
        if((Nodes1[1]^1) == (startNode^1) || (Nodes1[1]^1) == endNode)
        {
            continue;
        }

        f2 = detect_bubble_end_with_bubbles(g, Nodes2[0]^1, Nodes1[1]^1, &convex2, &l2, NULL);


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
	uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, n_reduced_a = 0;


    if (!g->is_symm) asg_symm(g);


    buf_t b;
	memset(&b, 0, sizeof(buf_t));
	b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));


    buf_t bub;
    memset(&bub, 0, sizeof(buf_t));
    bub.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t nv = asg_arc_n(g, v);
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
            n_reduced += test_triangular_exact(g, b.b.a, b.b.n, v, b.S.a[0], max_dist, &bub);
            n_reduced_a += test_triangular_addition_exact(g, b.b.a, b.b.n, v, b.S.a[0],max_dist, &bub);
        }

    }

    free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
    free(bub.a); free(bub.S.a); free(bub.T.a); free(bub.b.a); free(bub.e.a);

    if (n_reduced + n_reduced_a) {
		asg_cleanup(g);
		asg_symm(g);
	}

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d/%d triangular/triangular_a overlaps\n", 
        __func__, n_reduced, n_reduced_a);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

    return n_reduced + n_reduced_a;
}






int check_if_cross(asg_t *g, uint32_t v)
{
    uint32_t N_list[5] = {0};
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

	uint32_t v, w, n_vtx = g->n_seq * 2;
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
                for (i = 0; i < (long long)b.b.n; i++)
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

        if(check_cross == 1 && check_if_cross(g, v))
        {
            g->seq_vis[v] = 2;
        }
    }
    free(b.b.a);

    free(arg);

    return NULL;
}

int asg_arc_identify_simple_bubbles_multi(asg_t *g, int check_cross)
{
    double startTime = Get_T();
    memset(g->seq_vis, 0, g->n_seq*2*sizeof(uint8_t));

    pthread_t *_r_threads;

	_r_threads = (pthread_t *)malloc(sizeof(pthread_t)*asm_opt.thread_num);

    int i = 0;

	for (i = 0; i < asm_opt.thread_num; i++)
	{
        para_for_simple_bub* arg = (para_for_simple_bub*)malloc(sizeof(*arg));
        arg->g = g;
        arg->thread_num = asm_opt.thread_num;
        arg->threadID = i;
        arg->check_cross = check_cross;

        pthread_create(_r_threads + i, NULL, asg_arc_identify_simple_bubbles_pthread, (void*)arg);
	}
    

    for (i = 0; i<asm_opt.thread_num; i++)
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


    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
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
       }
    }

    return 0;

}


int test_single_node_bubble(asg_t *g, uint32_t* nodes, uint32_t length, 
uint32_t startNode, uint32_t endNode)
{
    
    uint32_t i, v, w;
    uint32_t vEnd;
    int flag0, flag1;
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
    int flag0, flag1;
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
	uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0;
    buf_t b;
	if (!g->is_symm) asg_symm(g);
	memset(&b, 0, sizeof(buf_t));
	///set information for each node
	b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t nv = asg_arc_n(g, v);
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
	uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0;
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t nv = asg_arc_n(g, v);
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

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d small bubbles\n", __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

    return n_reduced;
}



int test_cross(asg_t *g, uint32_t* nodes, uint32_t length, 
uint32_t startNode, uint32_t endNode)
{
    uint32_t a1, a2;
    uint32_t N_list[5] = {0};
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
	uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0;
    buf_t b;
	if (!g->is_symm) asg_symm(g);
	memset(&b, 0, sizeof(buf_t));
	///set information for each node
	b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t nv = asg_arc_n(g, v);
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

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d cross\n", __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

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

    if(VERBOSE >= 1)
    {
	    fprintf(stderr, "[M::%s] transitively reduced %d arcs\n", __func__, n_reduced);
    }

	if (n_reduced) {
		asg_cleanup(g);
		asg_symm(g);
	}

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
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
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] cut %d tips\n", __func__, cnt);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
	return cnt;
}


///max_ext is 4
int asg_cut_tip_primary(asg_t *g, ma_ug_t *ug, int max_ext)
{
    double startTime = Get_T();

	asg64_v a = {0,0,0};
	uint32_t n_vtx = g->n_seq * 2, v, i, cnt = 0, tipEvaluateLen;
    
	for (v = 0; v < n_vtx; ++v) {
		//if this seq has been deleted
		if (g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
		///asg_arc_n(v^1) == 0
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
        if(ug!=NULL)
        {
            tipEvaluateLen = 0;
            for (i = 0; i < a.n; ++i)
            {
                tipEvaluateLen += EvaluateLen(ug->u, ((uint32_t)a.a[i]>>1));
            }
            if(tipEvaluateLen > (uint32_t)max_ext) continue;
        }
    
        for (i = 0; i < a.n; ++i)
        {
            g->seq[((uint32_t)a.a[i]>>1)].c = ALTER_LABLE;
        }

        for (i = 0; i < a.n; ++i)
        {
            asg_seq_drop(g, (uint32_t)a.a[i]>>1);
        }

		++cnt;
	}
	free(a.a);
	if (cnt > 0) asg_cleanup(g);
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] cut %d tips\n", __func__, cnt);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
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
int asg_arc_del_short_diploid_unclean(asg_t *g, float drop_ratio, ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources)
{
    double startTime = Get_T();

	uint32_t v, n_vtx = g->n_seq * 2, n_short = 0;
    uint32_t last_e;
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
		for (i = nv - 1; i >= 1 && av[i].ol < thres; --i) {}
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




/*****************************read graph*****************************/

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

uint32_t detect_single_path_with_dels_n_stops(asg_t *g, uint32_t begNode, uint32_t* endNode, 
long long* Len, long long* max_stop_Len, buf_t* b, uint32_t stops_threshold)
{
    
    uint32_t v = begNode, w;
    uint32_t kv, kw, n_stops = 0, flag = LONG_TIPS;
    (*Len) = 0;
    (*max_stop_Len) = 0;
    long long preLen = 0, currentLen;

    while (1)
    {
        (*Len)++;
        kv = get_real_length(g, v, NULL);
        (*endNode) = v;

        if(b) kv_push(uint32_t, b->b, v>>1);

        if(kv == 0)
        {
            flag = END_TIPS;
            break;
        }

        if(kv == 2)
        {
            flag = TWO_OUTPUT;
            break;
        }

        if(kv > 2)
        {
            flag = MUL_OUTPUT;
            break;
        }

        ///up to here, kv=1
        ///kw must >= 1
        get_real_length(g, v, &w);
        kw = get_real_length(g, w^1, NULL);
        v = w;
        (*endNode) = v;

        if(kw >= 2)
        {
            n_stops++;
            currentLen = (*Len) - preLen;
            preLen = (*Len);
            if(currentLen > (*max_stop_Len))
            {
                (*max_stop_Len) = currentLen;
            }
        }
        if(kw >= 2 && n_stops >= stops_threshold)
        {
            (*Len)++;
            if(b) kv_push(uint32_t, b->b, v>>1);
            if(kw == 2) flag = TWO_INPUT;
            if(kw > 2) flag = MUL_INPUT;
            break;
        }
        
        if((v>>1) == (begNode>>1))
        {
            flag = LOOP;
            break;
        } 
    }


    currentLen = (*Len) - preLen;
    preLen = (*Len);
    if(currentLen > (*max_stop_Len))
    {
        (*max_stop_Len) = currentLen;
    }
    return flag;   
}

uint32_t detect_single_path_with_dels_contigLen(asg_t *g, uint32_t begNode, uint32_t* endNode, long long* baseLen, buf_t* b)
{
    
    uint32_t v = begNode, w = 0;
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

uint32_t detect_single_path_with_dels_contigLen_complex(asg_t *g, uint32_t begNode, uint32_t* endNode, 
long long* baseLen, long long* max_stop_base_Len, buf_t* b, uint32_t stops_threshold)
{
    
    uint32_t v = begNode, w = 0;
    uint32_t kv, kw, k, n_stops = 0, flag = LONG_TIPS;
    (*baseLen) = 0;
    (*max_stop_base_Len) = 0;
    long long preBaseLen = 0, currentBaseLen;


    while (1)
    {
        ///(*Len)++;
        kv = get_real_length(g, v, NULL);
        (*endNode) = v;

        if(b) kv_push(uint32_t, b->b, v>>1);

        if(kv == 0)
        {
            (*baseLen) += g->seq[v>>1].len;
            flag = END_TIPS;
            break;
        }

        if(kv == 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            flag = TWO_OUTPUT;
            break;
        }

        if(kv > 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            flag = MUL_OUTPUT;
            break;
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

        if(kw >= 2)
        {
            n_stops++;
            currentBaseLen = (*baseLen) - preBaseLen;
            preBaseLen = (*baseLen);
            if(currentBaseLen > (*max_stop_base_Len))
            {
                (*max_stop_base_Len) = currentBaseLen;
            }
        }

        if(kw >= 2 && n_stops >= stops_threshold)
        {
            (*baseLen) += g->seq[v>>1].len;
            if(b) kv_push(uint32_t, b->b, v>>1);
            if(kw == 2) flag = TWO_INPUT;
            if(kw > 2) flag = MUL_INPUT;
            break;
        }


        if((v>>1) == (begNode>>1))
        {
            flag = LOOP;
            break;
        } 
    }

    currentBaseLen = (*baseLen) - preBaseLen;
    preBaseLen = (*baseLen);
    if(currentBaseLen > (*max_stop_base_Len))
    {
        (*max_stop_base_Len) = currentBaseLen;
    }

    return flag;   
}

/*****************************read graph*****************************/




/*****************************untig graph*****************************/

uint32_t untig_detect_single_path_with_dels(asg_t *g, ma_ug_t *ug, uint32_t begNode, 
uint32_t* endNode, long long* Len, buf_t* b)
{
    
    uint32_t v = begNode, w;
    uint32_t kv, kw;
    (*Len) = 0;


    while (1)
    {
        (*Len) = (*Len) + EvaluateLen((*ug).u, v>>1);
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

        if(kw == 2)
        {
            return TWO_INPUT;
        }

        if(kw > 2)
        {
            return MUL_INPUT;
        }   

        v = w;

        if((v>>1) == (begNode>>1))
        {
            return LOOP;
        } 
    }

    return LONG_TIPS;   
}

uint32_t untig_detect_single_path_with_dels_n_stops(asg_t *g, ma_ug_t *ug, uint32_t begNode, 
uint32_t* endNode, long long* Len, long long* max_stop_Len, buf_t* b, uint32_t stops_threshold)
{
    
    uint32_t v = begNode, w;
    uint32_t kv, kw, n_stops = 0, flag = LONG_TIPS;
    (*Len) = 0;
    (*max_stop_Len) = 0;
    long long preLen = 0, currentLen;

    while (1)
    {
        (*Len) = (*Len) + EvaluateLen((*ug).u, v>>1);
        kv = get_real_length(g, v, NULL);
        (*endNode) = v;

        if(b) kv_push(uint32_t, b->b, v>>1);

        if(kv == 0)
        {
            flag = END_TIPS;
            break;
        }

        if(kv == 2)
        {
            flag = TWO_OUTPUT;
            break;
        }

        if(kv > 2)
        {
            flag = MUL_OUTPUT;
            break;
        }

        ///up to here, kv=1
        ///kw must >= 1
        get_real_length(g, v, &w);
        kw = get_real_length(g, w^1, NULL);

        ///just calculate the max_stop_Len
        if(kw >= 2)
        {
            n_stops++;
            currentLen = (*Len) - preLen;
            preLen = (*Len);
            if(currentLen > (*max_stop_Len))
            {
                (*max_stop_Len) = currentLen;
            }
        }
        if(kw >= 2 && n_stops >= stops_threshold)
        {
            if(kw == 2) flag = TWO_INPUT;
            if(kw > 2) flag = MUL_INPUT;
            break;
        }
        
        v = w; 

        if((v>>1) == (begNode>>1))
        {
            flag = LOOP;
            break;
        }
    }


    currentLen = (*Len) - preLen;
    preLen = (*Len);
    if(currentLen > (*max_stop_Len))
    {
        (*max_stop_Len) = currentLen;
    }
    return flag;   
}

uint32_t untig_detect_single_path_with_dels_contigLen(asg_t *g, uint32_t begNode, uint32_t* endNode, 
long long* baseLen, buf_t* b)
{
    
    uint32_t v = begNode, w = 0;
    uint32_t kv, kw, k;
    (*baseLen) = 0;


    while (1)
    {
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


        ///up to here, kv=1
        ///kw must >= 1
        get_real_length(g, v, &w);
        kw = get_real_length(g, w^1, NULL);
        
        if(kw == 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            return TWO_INPUT;
        }

        if(kw > 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            return MUL_INPUT;
        }


        for (k = 0; k < asg_arc_n(g, v); k++)
        {
            if(!asg_arc_a(g, v)[k].del)
            {
                w = asg_arc_a(g, v)[k].v;
                (*baseLen) += ((uint32_t)(asg_arc_a(g, v)[k].ul));
                break;
            }
        }


        v = w; 

        if((v>>1) == (begNode>>1))
        {
            return LOOP;
        }
    }

    return LONG_TIPS;   
}

uint32_t untig_detect_single_path_with_dels_contigLen_complex(asg_t *g, uint32_t begNode, uint32_t* endNode, 
long long* baseLen, long long* max_stop_base_Len, buf_t* b, uint32_t stops_threshold)
{
    
    uint32_t v = begNode, w = 0;
    uint32_t kv, kw, k, n_stops = 0, flag = LONG_TIPS;
    (*baseLen) = 0;
    (*max_stop_base_Len) = 0;
    long long preBaseLen = 0, currentBaseLen;


    while (1)
    {
        kv = get_real_length(g, v, NULL);
        (*endNode) = v;

        if(b) kv_push(uint32_t, b->b, v>>1);

        if(kv == 0)
        {
            (*baseLen) += g->seq[v>>1].len;
            flag = END_TIPS;
            break;
        }

        if(kv == 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            flag = TWO_OUTPUT;
            break;
        }

        if(kv > 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            flag = MUL_OUTPUT;
            break;
        }

        

        ///up to here, kv=1
        ///kw must >= 1
        get_real_length(g, v, &w);
        kw = get_real_length(g, w^1, NULL);

        if(kw == 1)
        {
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
        }

        
        ///just calculate the max_stop_Len
        if(kw >= 2)
        {
            n_stops++;
            if(n_stops >= stops_threshold)
            {
                (*baseLen) += g->seq[v>>1].len;
            }
            else
            {
                for (k = 0; k < asg_arc_n(g, v); k++)
                {
                    if(!asg_arc_a(g, v)[k].del)
                    {
                        w = asg_arc_a(g, v)[k].v;
                        (*baseLen) += ((uint32_t)(asg_arc_a(g, v)[k].ul));
                        break;
                    }
                }
            }
            
            currentBaseLen = (*baseLen) - preBaseLen;
            preBaseLen = (*baseLen);
            if(currentBaseLen > (*max_stop_base_Len))
            {
                (*max_stop_base_Len) = currentBaseLen;
            }
        }

        if(kw >= 2 && n_stops >= stops_threshold)
        {
            if(kw == 2) flag = TWO_INPUT;
            if(kw > 2) flag = MUL_INPUT;
            break;
        }
        
        
        v = w; 

        if((v>>1) == (begNode>>1))
        {
            flag = LOOP;
            break;
        }
    }

    currentBaseLen = (*baseLen) - preBaseLen;
    preBaseLen = (*baseLen);
    if(currentBaseLen > (*max_stop_base_Len))
    {
        (*max_stop_base_Len) = currentBaseLen;
    }

    return flag;   
}

/*****************************untig graph*****************************/







long long check_if_diploid(uint32_t v1, uint32_t v2, asg_t *g, 
ma_hit_t_alloc* reverse_sources, long long min_edge_length, R_to_U* ruIndex)
{
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    uint32_t convex1, convex2;
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
    uint32_t qn, tn, is_Unitig;
    for (i = 0; i < (long long)b_min->b.n; i++)
    {
        qn = b_min->b.a[i];
        for (j = 0; j < (long long)reverse_sources[qn].length; j++)
        {
            tn = Get_tn(reverse_sources[qn].buffer[j]);
            /****************************may have bugs********************************/
            ///if(g->seq[tn].del == 1) continue;
            if(g->seq[tn].del == 1)
            {
                get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
                if(tn == (uint32_t)-1 || is_Unitig == 1 || g->seq[tn].del == 1) continue;
            }
            /****************************may have bugs********************************/
            min_count++;
            for (k = 0; k < (long long)b_max->b.n; k++)
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


long long check_if_diploid_primary_complex(uint32_t v1, uint32_t v2, asg_t *g, 
ma_hit_t_alloc* reverse_sources, long long min_edge_length, uint32_t stops_threshold,
int if_drop, R_to_U* ruIndex)
{
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    uint32_t convex1, convex2;
    long long l1, l2, max_stop_Len;

    b_0.b.n = 0;
    b_1.b.n = 0;
    uint32_t flag1 = detect_single_path_with_dels_n_stops(g, v1, &convex1, &l1, 
    &max_stop_Len, &b_0, stops_threshold);
    
    uint32_t flag2 = detect_single_path_with_dels_n_stops(g, v2, &convex2, &l2, 
    &max_stop_Len, &b_1, stops_threshold);
    
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


    long long i, j, k;
    i = b_0.b.n; i--;
    j = b_1.b.n; j--;
    while (i>=0 && j>=0)
    {
        if(b_0.b.a[i] == b_1.b.a[j])
        {
            i--;
            j--;
        }
        else
        {
            break;
        }
        
    }
    b_0.b.n = i+1; l1 = i+1;
    b_1.b.n = j+1; l2 = j+1;


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

    
    double max_count = 0;
    double min_count = 0;
    uint32_t qn, tn, is_Unitig;
    for (i = 0; i < (long long)b_min->b.n; i++)
    {
        qn = b_min->b.a[i];
        for (j = 0; j < (long long)reverse_sources[qn].length; j++)
        {
            tn = Get_tn(reverse_sources[qn].buffer[j]);
            /****************************may have bugs********************************/
            ///if(g->seq[tn].del == 1 || (if_drop == 1 && g->seq[tn].c == ALTER_LABLE)) continue;
            if(g->seq[tn].del == 1)
            {
                get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
                if(tn == (uint32_t)-1 || is_Unitig == 1 || g->seq[tn].del == 1) continue;
            }
            /****************************may have bugs********************************/
            min_count++;
            for (k = 0; k < (long long)b_max->b.n; k++)
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


uint32_t get_long_tip_length_stops(asg_t *sg, ma_utg_v* u, uint32_t begNode, uint32_t* endNode, 
buf_t* b, uint32_t stops_threshold)
{
    uint32_t v = begNode, w, n_stops = 0;
    uint32_t kv, kw;
    uint32_t eLen = 0;
    (*endNode) = (uint32_t)-1;
    if(u->a[v>>1].circ) return eLen;
    while (1)
    {
        kv = get_real_length(sg, v, NULL);
        (*endNode) = v;
        eLen += EvaluateLen((*u), v>>1);
        if(b) kv_push(uint32_t, b->b, v);
        if(kv!=1) return eLen;
        ///kw must be 1 here
        kw = get_real_length(sg, v, &w);
        kw = get_real_length(sg, w^1, NULL);
        ///if(get_real_length(sg, w^1, NULL)!=1) return eLen;
        if(kw >= 2)
        {
            n_stops++;
        }
        if(kw >= 2 && n_stops >= stops_threshold)
        {
            return eLen;
        }
        v = w;
        if(v == begNode) return eLen;
    }
}



uint32_t get_long_tip_length(asg_t *sg, ma_utg_v* u, uint32_t begNode, uint32_t* endNode, buf_t* b)
{
    uint32_t v = begNode, w;
    uint32_t kv;
    uint32_t eLen = 0;
    (*endNode) = (uint32_t)-1;
    if(u->a[v>>1].circ) return eLen;
    while (1)
    {
        kv = get_real_length(sg, v, NULL);
        (*endNode) = v;
        eLen += EvaluateLen((*u), v>>1);
        if(b) kv_push(uint32_t, b->b, v);
        if(kv!=1) return eLen;
        ///kv must be 1 here
        kv = get_real_length(sg, v, &w);
        if(get_real_length(sg, w^1, NULL)!=1) return eLen;
        v = w;
        if(v == begNode)
        {
            u->a[begNode>>1].circ = 1;
            return eLen;
        } 
    }
}

uint32_t if_long_tip_length(asg_t *sg, ma_utg_v* u, uint32_t begNode, uint32_t* untigLen, 
long long minLongUntig, long long maxShortUntig, float ShortUntigRate, long long mainLen)
{
    uint32_t Len, endNode;
    if(untigLen == NULL)
    {
        Len = get_long_tip_length(sg, u, begNode, &endNode, NULL);
    }
    else
    {
        Len = (*untigLen);
    }
    
    
    if(Len == 0) return 0;
    if(Len < (ShortUntigRate*mainLen)) return 0;
    if(Len >= maxShortUntig) return 1;
    if(Len >= minLongUntig && Len >= (ShortUntigRate*mainLen)) return 1;
    return 0;
}


long long check_if_diploid_untigs(asg_t *nsg, asg_t *read_sg, uint32_t v_0, uint32_t v_1, 
ma_utg_v* ut_v,  ma_hit_t_alloc* reverse_sources, long long min_edge_length, 
float hap_rate, buf_t* b_0, buf_t* b_1, int if_drop, R_to_U* ruIndex)
{
    uint32_t vEnd;
    b_0->b.n = b_1->b.n = 0;
    if(get_long_tip_length(nsg, ut_v, v_0, &vEnd, b_0) == 0) return -1;
    if(get_long_tip_length(nsg, ut_v, v_1, &vEnd, b_1) == 0) return -1;
    uint32_t l_0 = 0, l_1 = 0, i;
    for (i = 0; i < b_0->b.n; i++)
    {
        l_0 += ut_v->a[(b_0->b.a[i])>>1].n;
    }

    for (i = 0; i < b_1->b.n; i++)
    {
        l_1 += ut_v->a[(b_1->b.a[i])>>1].n;
    }
    
    if((long long)(l_0) <= min_edge_length || (long long)(l_1) <= min_edge_length)
    {
        return -1;
    }

    buf_t* b_max;
    buf_t* b_min;
    
    if(l_0 <= l_1)
    {
        b_min = b_0;
        b_max = b_1;
    }
    else
    {
        b_min = b_1;
        b_max = b_0;
    }


    

    uint32_t max_count = 0;
    uint32_t min_count = 0;
    ma_utg_t* node_a;
    uint32_t i_a_1, i_a_2, untigID_a, qn, tn, j;
    ma_utg_t* node_b;
    uint32_t i_b_1, i_b_2, untigID_b;
    uint32_t is_Unitig;
    for (i_a_1 = 0; i_a_1 < b_min->b.n; i_a_1++)
    {
        untigID_a = b_min->b.a[i_a_1]>>1;
        node_a = &(ut_v->a[untigID_a]);
        for (i_a_2 = 0; i_a_2 < node_a->n; i_a_2++)
        {
            qn = node_a->a[i_a_2]>>33;

            for (j = 0; j < reverse_sources[qn].length; j++)
            {
                tn = Get_tn(reverse_sources[qn].buffer[j]);
                /****************************may have bugs********************************/
                ///here is bug
                // if(read_sg->seq[tn].del == 1 || (if_drop == 1 && read_sg->seq[tn].c == ALTER_LABLE)) 
                // {
                //     continue;
                // }
                if(read_sg->seq[tn].del == 1)
                {
                    get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
                    if(tn == (uint32_t)-1 || is_Unitig == 1 || read_sg->seq[tn].del == 1) continue;
                }
                /****************************may have bugs********************************/
                min_count++;



                for (i_b_1 = 0; i_b_1 < b_max->b.n; i_b_1++)
                {
                    untigID_b = b_max->b.a[i_b_1]>>1;
                    node_b = &(ut_v->a[untigID_b]);
                    for (i_b_2 = 0; i_b_2 < node_b->n; i_b_2++)
                    {
                        if(tn == (node_b->a[i_b_2]>>33))
                        {
                            max_count++;
                            goto end_reverse;
                        }
                    }
                }
                end_reverse:;
            }
        }
    }
    


    if(min_count == 0) return -1;
    if(max_count == 0) return 0;
    if(max_count > min_count*hap_rate) return 1;
    return 0;
}

long long check_if_diploid_untigs_complex(asg_t *nsg, asg_t *read_sg, uint32_t v_0, uint32_t v_1, 
ma_utg_v* ut_v,  ma_hit_t_alloc* reverse_sources, long long min_edge_length, long long stops_threshold,
float hap_rate, buf_t* b_0, buf_t* b_1, int if_drop, R_to_U* ruIndex)
{
    uint32_t vEnd;
    b_0->b.n = b_1->b.n = 0;

    
    if(get_long_tip_length_stops(nsg, ut_v, v_0, &vEnd, b_0, stops_threshold) == 0) return -1;
    if(get_long_tip_length_stops(nsg, ut_v, v_1, &vEnd, b_1, stops_threshold) == 0) return -1;


    uint32_t l_0 = 0, l_1 = 0;
    long long i, j;
    i = b_0->b.n; i--;
    j = b_1->b.n; j--;
    while (i>=0 && j >=0)
    {
        if(b_0->b.a[i] == b_1->b.a[j])
        {
            i--;
            j--;
        }
        else
        {
            break;
        }
    }
    b_0->b.n = i+1;
    b_1->b.n = j+1;


    l_0 = 0; l_1 = 0;
    for (i = 0; i < (long long)b_0->b.n; i++)
    {
        l_0 += ut_v->a[(b_0->b.a[i])>>1].n;
    }

    for (i = 0; i < (long long)b_1->b.n; i++)
    {
        l_1 += ut_v->a[(b_1->b.a[i])>>1].n;
    }
    
    if((long long)(l_0) <= min_edge_length || (long long)(l_1) <= min_edge_length)
    {
        return -1;
    }

    buf_t* b_max;
    buf_t* b_min;
    
    if(l_0 <= l_1)
    {
        b_min = b_0;
        b_max = b_1;
    }
    else
    {
        b_min = b_1;
        b_max = b_0;
    }


    

    uint32_t max_count = 0;
    uint32_t min_count = 0;
    ma_utg_t* node_a;
    uint32_t i_a_1, i_a_2, untigID_a, qn, tn;
    ma_utg_t* node_b;
    uint32_t i_b_1, i_b_2, untigID_b;
    uint32_t is_Unitig;
    for (i_a_1 = 0; i_a_1 < b_min->b.n; i_a_1++)
    {
        untigID_a = b_min->b.a[i_a_1]>>1;
        node_a = &(ut_v->a[untigID_a]);
        for (i_a_2 = 0; i_a_2 < node_a->n; i_a_2++)
        {
            qn = node_a->a[i_a_2]>>33;

            for (j = 0; j < (long long)reverse_sources[qn].length; j++)
            {
                tn = Get_tn(reverse_sources[qn].buffer[j]);
                ///here is bug
                /****************************may have bugs********************************/
                // if(read_sg->seq[tn].del == 1 || (if_drop == 1 && read_sg->seq[tn].c == ALTER_LABLE)) 
                // {
                //     continue;
                // }
                if(read_sg->seq[tn].del == 1)
                {
                    get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
                    if(tn == (uint32_t)-1 || is_Unitig == 1 || read_sg->seq[tn].del == 1) continue;
                }
                /****************************may have bugs********************************/
                min_count++;



                for (i_b_1 = 0; i_b_1 < b_max->b.n; i_b_1++)
                {
                    untigID_b = b_max->b.a[i_b_1]>>1;
                    node_b = &(ut_v->a[untigID_b]);
                    for (i_b_2 = 0; i_b_2 < node_b->n; i_b_2++)
                    {
                        if(tn == (node_b->a[i_b_2]>>33))
                        {
                            max_count++;
                            goto end_reverse;
                        }
                    }
                }
                end_reverse:;
            }
        }
    }
    


    if(min_count == 0) return -1;
    if(max_count == 0) return 0;
    if(max_count > min_count*hap_rate) return 1;
    return 0;
}



long long check_if_diploid_aggressive(uint32_t v1, uint32_t v2, asg_t *g, 
ma_hit_t_alloc* reverse_sources, long long min_edge_length)
{
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    uint32_t convex1, convex2;
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
    for (i = 0; i < (long long)b_min->b.n; i++)
    {
        qn = b_min->b.a[i];
        for (j = 0; j < (long long)reverse_sources[qn].length; j++)
        {
            tn = Get_tn(reverse_sources[qn].buffer[j]);
            if(g->seq[tn].del == 1) continue;
            min_count++;
            for (k = 0; k < (long long)b_max->b.n; k++)
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

int asg_arc_del_too_short_overlaps(asg_t *g, long long dropLen, float drop_ratio, 
ma_hit_t_alloc* reverse_sources, long long min_edge_length, R_to_U* ruIndex)
{
    double startTime = Get_T();

	uint32_t v, v_max, v_maxLen, n_vtx = g->n_seq * 2, n_short = 0;
    long long drop_ratio_Len = 0;
    
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
                check_if_diploid(v_max, av[i].v, g, reverse_sources, min_edge_length, ruIndex) != 1)
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
	
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d short overlaps\n", __func__, n_short);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

	return n_short;
}



int asg_arc_del_short_diploid_unclean_exact(asg_t *g, float drop_ratio, ma_hit_t_alloc* sources)
{
	uint32_t v, n_vtx = g->n_seq * 2, n_short = 0;
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
    
	return n_ext;
}

// delete short arcs
///for best graph?
int asg_arc_del_short_diploid_by_length(asg_t *g, float drop_ratio, int max_ext, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen, uint32_t stops_threshold,
int if_skip_bubble, int if_drop, int if_check_hap, R_to_U* ruIndex)
{
    double startTime = Get_T();
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));

	uint32_t v, n_vtx = g->n_seq * 2;
    long long n_cut = 0;
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(if_skip_bubble && g->seq_vis[v] != 0) continue;
        if(if_drop && g->seq[v>>1].c == ALTER_LABLE) continue;
        asg_arc_t *av = asg_arc_a(g, v);
        uint32_t nv = asg_arc_n(g, v);
        if (nv < 2) continue;
        uint64_t i;
        for (i = 0; i < nv; ++i)
        {
            kv_push(uint64_t, b, (uint64_t)((uint64_t)av[i].ol << 32 | (av - g->arc + i)));   
        }
	}

    radix_sort_arch64(b.a, b.a + b.n);

    uint64_t k;
    for (k = 0; k < b.n; k++)
    {
        
        asg_arc_t *a = &g->arc[(uint32_t)b.a[k]];
		///v is self id, w is the id of another end
		uint32_t i, iv, iw, v = (a->ul)>>32, w = a->v^1, to_del = 0;
		uint32_t nv = asg_arc_n(g, v), nw = asg_arc_n(g, w), kv, kw;
		uint32_t ov_max = 0, ow_max = 0, ov_max_i = 0, ow_max_i = 0;
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
            if(if_check_hap == 1 && to_del == 1)
            {
                
                if(check_if_diploid_primary_complex(av[ov_max_i].v, w^1, g, reverse_sources, 
                miniedgeLen, stops_threshold, if_drop, ruIndex) == 1 
                   || 
                   check_if_diploid_primary_complex(aw[ow_max_i].v, v^1, g, reverse_sources, 
                miniedgeLen, stops_threshold, if_drop, ruIndex) == 1)
                {
                    to_del = 0;
                }
            }
        
		} else if (kw == 1) {
			if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
            ///kv > 1
            if(if_check_hap == 1 && to_del == 1)
            {
                if(check_if_diploid_primary_complex(av[ov_max_i].v, w^1, g, reverse_sources, 
                miniedgeLen, stops_threshold, if_drop, ruIndex) == 1)
                {
                    to_del = 0;
                }
            }
		} else if (kv == 1) {
			if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
            ///kw > 1
            if(if_check_hap == 1 && to_del == 1)
            {
                if(check_if_diploid_primary_complex(aw[ow_max_i].v, v^1, g, reverse_sources, 
                miniedgeLen, stops_threshold, if_drop, ruIndex) == 1)
                {
                    to_del = 0;
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

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %lld short overlaps\n", __func__, n_cut);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
	return n_cut;
}


int asg_arc_del_short_false_link(asg_t *g, float drop_ratio, float o_drop_ratio, int max_dist, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen, R_to_U* ruIndex)
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

    radix_sort_arch64(b.a, b.a + b.n);

    uint32_t min_edge;
    
    
    uint64_t k, t;
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
                    if(check_if_diploid(b_f.a[0], b_f.a[i], g, reverse_sources, miniedgeLen, ruIndex) == 1)
                    {
                        to_del_l++;
                    }
                }

                ///backward
                if(to_del_l != kv && b_r.n >= 2)
                {
                    to_del_l = 0;
                    uint32_t w0 = 0, w1 = 0;


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

                        if(check_if_diploid(w0, w1, g, reverse_sources, miniedgeLen, ruIndex) == 1)
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

        
    



		if (to_del_l && to_del_r)
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
    
    free(b.a); free(b_f.a); free(b_r.a);
    free(bub.a); free(bub.S.a); free(bub.T.a); free(bub.b.a); free(bub.e.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %u false overlaps\n", __func__, n_cut);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
	return n_cut;
}


int asg_arc_del_short_false_link_primary(asg_t *g, float drop_ratio, float o_drop_ratio, int max_dist, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen, R_to_U* ruIndex)
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
        if(g->seq[v>>1].c == ALTER_LABLE) continue;

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

    radix_sort_arch64(b.a, b.a + b.n);

    uint32_t min_edge;
    
    
    uint64_t k, t;
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
                    /****************************may have bugs********************************/
                    if(check_if_diploid(b_f.a[0], b_f.a[i], g, reverse_sources, miniedgeLen, ruIndex) == 1)
                    {/****************************may have bugs********************************/
                        to_del_l++;
                    }
                }

                ///backward
                if(to_del_l != kv && b_r.n >= 2)
                {
                    to_del_l = 0;
                    uint32_t w0 = 0, w1 = 0;


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
                        /****************************may have bugs********************************/
                        if(check_if_diploid(w0, w1, g, reverse_sources, miniedgeLen, ruIndex) == 1)
                        {/****************************may have bugs********************************/
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

        
    



		if (to_del_l && to_del_r)
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
    
    free(b.a); free(b_f.a); free(b_r.a);
    free(bub.a); free(bub.S.a); free(bub.T.a); free(bub.b.a); free(bub.e.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %u false overlaps\n", __func__, n_cut);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
	return n_cut;
}

int asg_arc_del_short_false_link_advance(asg_t *g, float drop_ratio, float o_drop_ratio, int max_dist, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen, R_to_U* ruIndex)
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

    radix_sort_arch64(b.a, b.a + b.n);

    uint32_t min_edge;
    
    uint64_t k, t;
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
                    if(check_if_diploid(b_f.a[0], b_f.a[i], g, reverse_sources, miniedgeLen, ruIndex) == 1)
                    {
                        to_del_l++;
                    }
                }

                ///backward
                if(to_del_l != kv && b_r.n >= 2)
                {
                    to_del_l = 0;
                    uint32_t w0 = 0, w1 = 0;


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

                        if(check_if_diploid(w0, w1, g, reverse_sources, miniedgeLen, ruIndex) == 1)
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

  


		if (to_del_l && to_del_r)
        {
            ///fprintf(stderr, "%.*s\n", Get_NAME_LENGTH((R_INF), v>>1), Get_NAME((R_INF), v>>1));
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


    radix_sort_arch64(b.a, b.a + b.n);
    
    
    uint64_t k;
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
            if((l1 <= min_thres) || (l2 <= min_thres))
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

        ///here both f1 and f2 must be false
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
            
            ///since both f1 and f2 must be false before 'if(to_del == 0)'
            ///here l1 = min_thres + 10, l2 = min_thres + 10
            if(f1 && f2)
            {
                // if(l1 <= min_thres || l2 <= min_thres) // TODO: check this block: it is always false
                // {
                //     continue;
                // } 
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


    radix_sort_arch64(b.a, b.a + b.n);

    uint32_t min_edge;
    
    
    uint64_t k, t;
    for (k = 0; k < b.n; k++)
    {
        ///v is the node
        v = (uint32_t)b.a[k];
        if (g->seq[v>>1].del) continue;
        uint32_t nv = asg_arc_n(g, v), nw, to_del;
        if (nv < 2) continue;
        uint32_t kv = get_real_length(g, v, NULL), kw;
        if (kv < 2) continue;///V must have two out-nodes
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

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d false overlaps\n", __func__, n_cut);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
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
                kv_push(uint64_t, b, (uint64_t)((uint64_t)av[i].ol << 32 | (av - g->arc + i)));   
            }
        }
	}

    radix_sort_arch64(b.a, b.a + b.n);

    uint64_t k;
    for (k = 0; k < b.n; k++)
    {
        
        asg_arc_t *a = &g->arc[(uint32_t)b.a[k]];
		///v is self id, w is the id of another end
		uint32_t i, iv, iw, v = (a->ul)>>32, w = a->v^1, to_del = 0;
		uint32_t nv = asg_arc_n(g, v), nw = asg_arc_n(g, w), kv, kw;
		uint32_t ov_max = 0, ow_max = 0, ov_max_i = 0;
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

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %lld inexact overlaps\n", __func__, n_cut);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
	return n_cut;
}



int asg_arc_del_short_diploi_by_suspect_edge(asg_t *g, int max_ext)
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
                    kv_push(uint64_t, b, (uint64_t)((uint64_t)av[i].ol << 32 | (av - g->arc + i)));   
                }
            }
        }
	}


    radix_sort_arch64(b.a, b.a + b.n);


    uint64_t k;
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
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %lld suspect overlaps\n", __func__, n_cut);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
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
            kv_push(uint64_t, b, (uint64_t)((uint64_t)av[0].ol << 32 | (av - g->arc)));
        }
	}

    radix_sort_arch64(b.a, b.a + b.n);

    uint64_t k;
    ///here all edges are inexact matches 
    for (k = 0; k < b.n; k++)
    {
        
        asg_arc_t *a = &g->arc[(uint32_t)b.a[k]];
		///v is self id, w is the id of another end
		uint32_t i, iv, iw, v = (a->ul)>>32, w = a->v^1;
		uint32_t nv = asg_arc_n(g, v), nw = asg_arc_n(g, w), kv, kw;
		asg_arc_t *av, *aw;
        av = asg_arc_a(g, v);
		aw = asg_arc_a(g, w);


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

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %lld single nodes\n", __func__, n_cut);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
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
    ///each node has two directions
	mark = (int32_t*)calloc(n_vtx, 4);

	q = kdq_init(uint64_t);
	for (v = 0; v < n_vtx; ++v) {
		uint32_t w, x, l, start, end, len;
		ma_utg_t *p;
        ///what's the usage of mark array
        ///mark array is used to mark if this node has already been included in a contig
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

    //mark all start nodes and end nodes of all unitigs
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
        ///>=0 means this node is a start/end node of an unitig
        ///means this node is a intersaction node
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





ma_ug_t *ma_ug_gen_primary(asg_t *g, uint8_t flag)
{
	int32_t *mark;
	uint32_t i, v, n_vtx = g->n_seq * 2;
	///is a queue
	kdq_t(uint64_t) *q;
	ma_ug_t *ug;

	ug = (ma_ug_t*)calloc(1, sizeof(ma_ug_t));
	ug->g = asg_init();
    ///each node has two directions
	mark = (int32_t*)calloc(n_vtx, 4);

    ///for each untig, all node have the same direction
    ///and all node except the last one just have one edge
    ///the last one may have multiple edges
	q = kdq_init(uint64_t);
	for (v = 0; v < n_vtx; ++v) {
		uint32_t w, x, l, start, end, len;
		ma_utg_t *p;
        ///what's the usage of mark array
        ///mark array is used to mark if this node has already been included in a contig
        /****************************may have hap bugs********************************/
		///if (g->seq[v>>1].del || arc_cnt(g, v) == 0 || mark[v] || g->seq[v>>1].c != flag) continue;
        ///if (g->seq[v>>1].del || arc_cnt(g, v) == 0 || mark[v] || (g->seq[v>>1].c & flag) == 0) continue;
        if (g->seq[v>>1].del || arc_cnt(g, v) == 0 || mark[v]) continue;
        if(flag == PRIMARY_LABLE && g->seq[v>>1].c == ALTER_LABLE) continue;
        if(flag == ALTER_LABLE && g->seq[v>>1].c != ALTER_LABLE) continue;
        /****************************may have hap bugs********************************/
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
        ///kdq_size(q) == 0 means there is just one read
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

    //mark all start nodes and end nodes of all unitigs
    ///note ug->u.a[i].start == ug->u.a[i].a[0]
    ///but ug->u.a[i].end = ug->u.a[i].a[ug->u.a[i].n-1]^1
    ///ug->u.a[i].start has the same direction as other non-end node
    ///while ug->u.a[i].end is the only one with reverse direction
	for (i = 0; i < ug->u.n; ++i) {
		if (ug->u.a[i].circ) continue;
        ///i is the untig id
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
        ///>=0 means this node is a start/end node of an unitig
        ///means this node is a intersaction node
        /**for one untig,
          start                end
         ----->              <------
         so if we want to find the edge between two untigs, their start and end look like: 

         start                end
         ----->              <------
         ---------------------------

                                    start                end
                                    ----->              <------
                                    ---------------------------
        p->ul>>32^1 is the end of the first untig,
        p->v is the start of the second untig
        **/
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
int ma_ug_seq_back(ma_ug_t *g, All_reads *RNF, const ma_sub_t *coverage_cut,
const long long n_read)
{
    UC_Read g_read;
    init_UC_Read(&g_read);
	utg_intv_t *tmp;
	uint32_t i, j;

    
    ///why we need n_read here? it is just beacuse one read can only be included in one untig
    ///but it is not true
	tmp = (utg_intv_t*)calloc(n_read, sizeof(utg_intv_t));
	///number of unitigs
	for (i = 0; i < g->u.n; ++i) {
		ma_utg_t *u = &g->u.a[i];
		uint32_t l = 0;
		u->s = (char*)calloc(1, u->len + 1);
		memset(u->s, 'N', u->len);
		for (j = 0; j < u->n; ++j) {
            ///u->a[j]>>33 is the readID
			utg_intv_t *t = &tmp[u->a[j]>>33];
			///assert(t->len == 0);
			t->utg = i, t->ori = u->a[j]>>32&1;
            ///l is the start pos of this read at its corresponding untig
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


// generate unitig sequences
int ma_ug_seq(ma_ug_t *g, All_reads *RNF, const ma_sub_t *coverage_cut,
const long long n_read)
{
    UC_Read g_read;
    init_UC_Read(&g_read);
	///utg_intv_t *tmp;
	uint32_t i, j, k;
    uint32_t rId, /**uId,**/ori, start, eLen, readLen;
    char* readS = NULL;

    
    ///why we need n_read here? it is just beacuse one read can only be included in one untig
    ///but it is not true
	///tmp = (utg_intv_t*)calloc(n_read, sizeof(utg_intv_t));
	///number of unitigs
	for (i = 0; i < g->u.n; ++i) {
		ma_utg_t *u = &g->u.a[i];
        if(u->m == 0) continue;
		uint32_t l = 0;
		u->s = (char*)calloc(1, u->len + 1);
		memset(u->s, 'N', u->len);
		for (j = 0; j < u->n; ++j) {
            rId = u->a[j]>>33;
            ///uId = i;
            ori = u->a[j]>>32&1;
            start = l;
            eLen = (uint32_t)u->a[j];
			l += eLen;

            if(eLen == 0) continue;
            recover_UC_Read(&g_read, RNF, rId);
            readS = g_read.seq + coverage_cut[rId].s;
            readLen = coverage_cut[rId].e - coverage_cut[rId].s;
            if (!ori) // forward strand
            {
                for (k = 0; k < eLen; k++)
                {
                    u->s[start + k] = readS[k];
                }
            }
            else
            {
                for (k = 0; k < eLen; k++)
                {
                    uint8_t c = (uint8_t)readS[readLen - 1 - k];
                    u->s[start + k] = c >= 128? 'N' : comp_tab[c];
                }
            }
            

		}
	}

    destory_UC_Read(&g_read);
	return 0;
}


void ma_ug_print2(const ma_ug_t *ug, All_reads *RNF, const ma_sub_t *coverage_cut, int print_seq, FILE *fp)
{
	uint32_t i, j, l;
	char name[32];
	for (i = 0; i < ug->u.n; ++i) { // the Segment lines in GFA
		ma_utg_t *p = &ug->u.a[i];
        if(p->m == 0) continue;
		sprintf(name, "utg%.6d%c", i + 1, "lc"[p->circ]);
		if (print_seq) fprintf(fp, "S\t%s\t%s\tLN:i:%d\n", name, p->s? p->s : "*", p->len);
		else fprintf(fp, "S\t%s\t*\tLN:i:%d\n", name, p->len);
        
		for (j = l = 0; j < p->n; l += (uint32_t)p->a[j++]) {
			uint32_t x = p->a[j]>>33;
			fprintf(fp, "A\t%s\t%d\t%c\t%.*s\t%d\t%d\tid:i:%d\n", name, l, "+-"[p->a[j]>>32&1], (int)Get_NAME_LENGTH((*RNF), x),
					Get_NAME((*RNF), x), coverage_cut[x].s, coverage_cut[x].e, x);
        }
	}
	for (i = 0; i < ug->g->n_arc; ++i) { // the Link lines in GFA
		uint32_t u = ug->g->arc[i].ul>>32, v = ug->g->arc[i].v;
		fprintf(fp, "L\tutg%.6d%c\t%c\tutg%.6d%c\t%c\t%dM\tL1:i:%d\n", (u>>1)+1, "lc"[ug->u.a[u>>1].circ], "+-"[u&1],
				(v>>1)+1, "lc"[ug->u.a[v>>1].circ], "+-"[v&1], ug->g->arc[i].ol, asg_arc_len(ug->g->arc[i]));
	}    
}

void ma_ug_print(const ma_ug_t *ug, All_reads *RNF, const ma_sub_t *coverage_cut, FILE *fp)
{
	ma_ug_print2(ug, RNF, coverage_cut, 1, fp);
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
    long long i, j, index;
    uint32_t qn, tn;

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

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
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
            command, (int)Get_NAME_LENGTH(R_INF, (v>>1)), Get_NAME(R_INF, (v>>1)), v&1);
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
                (int)Get_NAME_LENGTH(R_INF, (av[i].v>>1)), 
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
    uint32_t tn;

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

            fprintf(stderr, "****************ma_hit_t (%lld)ref_read: %.*s****************\n", 
            i, (int)Get_NAME_LENGTH(R_INF, i), Get_NAME(R_INF, i));


            fprintf(stderr, "sources Len: %d, is_fully_corrected: %d\n", 
            sources[i].length, sources[i].is_fully_corrected);

            for (j = 0; j < sources[i].length; j++)
            {
                tn = Get_tn(sources[i].buffer[j]);
                fprintf(stderr, "target: %.*s, qs: %u, qe: %u, ts: %u, te: %u, ml: %u, rev: %u, el: %u\n", 
                (int)Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn),
                Get_qs(sources[i].buffer[j]),
                Get_qe(sources[i].buffer[j]),
                Get_ts(sources[i].buffer[j]),
                Get_te(sources[i].buffer[j]),
                sources[i].buffer[j].ml,
                sources[i].buffer[j].rev,
                sources[i].buffer[j].el);
            }


            
            fprintf(stderr, "######reverse_query_read Len: %d\n", reverse_sources[i].length);



            for (j = 0; j < reverse_sources[i].length; j++)
            {
                tn = Get_tn(reverse_sources[i].buffer[j]);
                fprintf(stderr, "target: %.*s, qs: %u, qe: %u, ts: %u, te: %u\n", 
                (int)Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn),
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
             "S\t%.*s\t*\tLN:i:%u\n",
             (int)Get_NAME_LENGTH((*RNF), i), 
             Get_NAME((*RNF), i),
             g->seq[i].len);
        }
    }

	for (i = 0; i < g->n_arc; ++i) {
		const asg_arc_t *p = &g->arc[i];
		if (sub) {
			const ma_sub_t *sq = &sub[p->ul>>33], *st = &sub[p->v>>1];

            fprintf(fp, 
            "L\t%.*s:%d-%d\t%c\t%.*s:%d-%d\t%c\t%d:\tL1:i:%u\n", 
            (int)Get_NAME_LENGTH((*RNF), p->ul>>33), 
            Get_NAME((*RNF), p->ul>>33),
            sq->s + 1, sq->e, "+-"[p->ul>>32&1],
            (int)Get_NAME_LENGTH((*RNF), p->v>>1), 
            Get_NAME((*RNF), p->v>>1),
            st->s + 1, st->e, "+-"[p->v&1], p->ol, (uint32_t)p->ul);
            

		} 
        else 
        {
           fprintf(fp, "L\t%.*s\t%c\t%.*s\t%c\t%d:\tL1:i:%u\n", 
            (int)Get_NAME_LENGTH((*RNF), p->ul>>33), 
            Get_NAME((*RNF), p->ul>>33),
            "+-"[p->ul>>32&1],
            (int)Get_NAME_LENGTH((*RNF), p->v>>1), 
            Get_NAME((*RNF), p->v>>1),
            "+-"[p->v&1], p->ol, (uint32_t)p->ul);
		}
	}
}


void ma_ug_print_simple(const ma_ug_t *ug, All_reads *RNF, const ma_sub_t *coverage_cut, FILE *fp)
{
	ma_ug_print2(ug, RNF, coverage_cut, 0, fp);
}




int asg_arc_cut_long_tip_primary(asg_t *g, ma_ug_t *ug, float drop_ratio)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, v_maxLen_i = (uint32_t)-1;
    long long ll, v_maxLen;

    buf_t b;
    memset(&b, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v), flag;
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (nv < 2 || g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;

        v_maxLen = -1;
        v_maxLen_i = (uint32_t)-1;

        for (i = 0, n_arc = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                if(ug == NULL)
                {
                    detect_single_path_with_dels(g, av[i].v, &convex, &ll, NULL);
                }
                else
                {
                    untig_detect_single_path_with_dels(g, ug, av[i].v, &convex, &ll, NULL);
                }
                
                
                if(v_maxLen < ll)
                {
                    v_maxLen = ll;
                    v_maxLen_i = i;
                }
            }
        }

        for (i = 0, n_arc = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                if(v_maxLen_i == i) continue;

                b.b.n = 0;

                if(ug == NULL)
                {
                    flag = detect_single_path_with_dels(g, av[i].v, &convex, &ll, &b);
                }
                else
                {
                    flag = untig_detect_single_path_with_dels(g, ug, av[i].v, &convex, &ll, &b);
                }
                
                if(flag == END_TIPS)
                {
                    if(v_maxLen*drop_ratio > ll)
                    {
                        n_reduced++;
                        uint64_t k;

                        for (k = 0; k < b.b.n; k++)
                        {
                            g->seq[b.b.a[k]].c = ALTER_LABLE;
                        }

                        for (k = 0; k < b.b.n; k++)
                        {
                            asg_seq_drop(g, b.b.a[k]);
                        }
                    }
                }
            }
        }
    }


    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

    return n_reduced;
}


int asg_arc_cut_long_tip_primary_complex(asg_t *g, float drop_ratio, uint32_t stops_threshold)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag;
    long long ll, max_stopLen;

    buf_t b;
    memset(&b, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i;
        ///some node could be deleted
        if (g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        ///tip
        if (get_real_length(g, v, NULL) != 0) continue;
        if(get_real_length(g, v^1, NULL) != 1) continue;      
        flag = detect_single_path_with_dels(g, v^1, &convex, &ll, NULL);
        if(flag != TWO_INPUT && flag != MUL_INPUT) continue;
        convex = convex^1;ll--;
        uint32_t n_convex = asg_arc_n(g, convex), convexLen = ll;
        asg_arc_t *a_convex = asg_arc_a(g, convex);


        for (i = 0; i < n_convex; i++)
        {
            if (!a_convex[i].del)
            {
                ///if stops_threshold = 1, 
                ///detect_single_path_with_dels_n_stops() is detect_single_path_with_dels()
                detect_single_path_with_dels_n_stops(g, a_convex[i].v, &convex, &ll, &max_stopLen,
                NULL, stops_threshold);

                if(convex == v) continue;

                if(ll*drop_ratio > convexLen && max_stopLen*2>ll)
                {

                    b.b.n = 0; 
                    flag = detect_single_path_with_dels(g, v^1, &convex, &ll, &b);
                    if(b.b.n < 2) break;
                    b.b.n--;


                    n_reduced++;
                    uint64_t k;

                    for (k = 0; k < b.b.n; k++)
                    {
                        g->seq[b.b.a[k]].c = ALTER_LABLE;
                    }

                    for (k = 0; k < b.b.n; k++)
                    {
                        asg_seq_drop(g, b.b.a[k]);
                    }
                    break;
                }
            }
        }
    }


    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

    return n_reduced;
}


int untig_asg_arc_cut_long_tip_primary_complex(ma_ug_t *ug, float drop_ratio, 
uint32_t stops_threshold)
{
    double startTime = Get_T();
    asg_t *g = ug->g;
    
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag;
    long long ll, max_stopLen;

    buf_t buf;
    memset(&buf, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i;
        ///some node could be deleted
        if (g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        ///tip
        if (get_real_length(g, v, NULL) != 0) continue;
        if(get_real_length(g, v^1, NULL) != 1) continue;
        /****************************may have bugs********************************/
        buf.b.n = 0;       
        flag = untig_detect_single_path_with_dels(g, ug, v^1, &convex, &ll, &buf);
        if(flag != TWO_INPUT && flag != MUL_INPUT) continue;
        get_real_length(g, convex, &convex);
        /****************************may have bugs********************************/
        convex = convex^1;
        ///note the convexLen here
        uint32_t n_convex = asg_arc_n(g, convex), convexLen = ll;
        asg_arc_t *a_convex = asg_arc_a(g, convex);


        for (i = 0; i < n_convex; i++)
        {
            if (!a_convex[i].del)
            {
                ///if stops_threshold = 1, 
                ///untig_detect_single_path_with_dels_n_stops() is untig_detect_single_path_with_dels()
                untig_detect_single_path_with_dels_n_stops(g, ug, a_convex[i].v, &convex, &ll, 
                &max_stopLen, NULL, stops_threshold);
                if(convex == v) continue;

                if(ll*drop_ratio > convexLen && max_stopLen*2>ll)
                {
                    n_reduced++;
                    uint64_t k;

                    for (k = 0; k < buf.b.n; k++)
                    {
                        g->seq[buf.b.a[k]].c = ALTER_LABLE;
                    }

                    for (k = 0; k < buf.b.n; k++)
                    {
                        asg_seq_drop(g, buf.b.a[k]);
                    }
                    break;
                }
            }
        }
    }


    asg_cleanup(g);
    asg_symm(g);
    free(buf.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

    return n_reduced;
}



int untig_asg_arc_cut_long_equal_tips_assembly(ma_ug_t *ug, asg_t *read_sg, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen, R_to_U* ruIndex)
{
    asg_t *g = ug->g;
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag, is_hap;
    long long ll, base_maxLen, base_maxLen_i;
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    buf_t b;
    memset(&b, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v), n_tips;
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (nv < 2 || g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;

        base_maxLen = -1;
        base_maxLen_i = -1;
        n_tips = 0;
        is_hap = 0;

        ///there must be more than 1 out-edges
        for (i = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                flag = untig_detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, NULL);
                
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
                    if(untig_detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b) 
                    != END_TIPS)
                    {
                        continue;
                    }
                    //we can only cut tips
                    if(check_if_diploid_untigs(g, read_sg, av[base_maxLen_i].v, av[i].v, &(ug->u),
                    reverse_sources, miniedgeLen, 0.3, &b_0, &b_1, 1, ruIndex) == 1)
                    {
                        n_reduced++;
                        uint64_t k;

                        for (k = 0; k < b.b.n; k++)
                        {
                            g->seq[b.b.a[k]].c = ALTER_LABLE;
                        }
                        for (k = 0; k < b.b.n; k++)
                        {
                            asg_seq_drop(g, b.b.a[k]);
                        }

                        is_hap++;
                    }
                }
            }
        }

        if(is_hap > 0)
        {
            i = base_maxLen_i;
            b.b.n = 0;
            untig_detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b);
            uint64_t k;
            for (k = 0; k < b.b.n; k++)
            {
                g->seq[b.b.a[k]].c = HAP_LABLE;
            }
        }
    }


    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    free(b_0.b.a);
    free(b_1.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    

    return n_reduced;
}

int asg_arc_cut_long_equal_tips_assembly(asg_t *g, ma_hit_t_alloc* reverse_sources, 
long long miniedgeLen, R_to_U* ruIndex)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag, is_hap;
    long long ll, base_maxLen, base_maxLen_i;

    buf_t b;
    memset(&b, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v), n_tips;
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (nv < 2 || g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;

        base_maxLen = -1;
        base_maxLen_i = -1;
        n_tips = 0;
        is_hap = 0;

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
                    /****************************may have bugs********************************/
                    if(check_if_diploid(av[base_maxLen_i].v, av[i].v, g, 
                    reverse_sources, miniedgeLen, ruIndex)==1)
                    {/****************************may have bugs********************************/
                        n_reduced++;
                        uint64_t k;

                        for (k = 0; k < b.b.n; k++)
                        {
                            g->seq[b.b.a[k]].c = ALTER_LABLE;
                        }
                        for (k = 0; k < b.b.n; k++)
                        {
                            asg_seq_drop(g, b.b.a[k]);
                        }
                        is_hap++;
                    }
                }
            }
        }


        if(is_hap > 0)
        {
            i = base_maxLen_i;
            b.b.n = 0;
            detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b);

            uint64_t k;
            for (k = 0; k < b.b.n; k++)
            {
                g->seq[b.b.a[k]].c = HAP_LABLE;
            }
        }


    }


    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    

    return n_reduced;
}


int asg_arc_simple_large_bubbles(asg_t *g, ma_hit_t_alloc* reverse_sources, long long miniedgeLen, 
R_to_U* ruIndex)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag, is_hap;
    long long ll, base_maxLen, base_maxLen_i, all_covex;

    buf_t b;
    memset(&b, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (nv < 2 || g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;

        base_maxLen = -1;
        base_maxLen_i = -1;
        all_covex = -1;
        is_hap = 0;

        for (i = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                flag = detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, NULL);
                if(flag != TWO_INPUT && flag != MUL_INPUT)
                {
                    break;
                }

                if(all_covex != -1 && (uint32_t)all_covex != convex)
                {
                    break;
                }

                if(all_covex == -1)
                {
                    all_covex = convex;
                }
                
                if(base_maxLen < ll)
                {
                    base_maxLen = ll;
                    base_maxLen_i = i;
                }

            }
        }


        if(i == nv)
        {
            for (i = 0; i < nv; i++)
            {
                if(i == base_maxLen_i) continue;
                if (!av[i].del)
                {
                    b.b.n = 0;
                    detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b);
                    if(b.b.n < 2) continue;
                    b.b.n--;

                    //we can only cut tips
                    /****************************may have bugs********************************/
                    if(check_if_diploid(av[base_maxLen_i].v, av[i].v, g, 
                    reverse_sources, miniedgeLen, ruIndex)==1)
                    {/****************************may have bugs********************************/
                        n_reduced++;
                        uint64_t k;

                        for (k = 0; k < b.b.n; k++)
                        {
                            g->seq[b.b.a[k]].c = ALTER_LABLE;
                        }
                        for (k = 0; k < b.b.n; k++)
                        {
                            asg_seq_drop(g, b.b.a[k]);
                        }

                        is_hap++;
                    }
                }
            }

            if(is_hap > 0)
            {
                i = base_maxLen_i;
                b.b.n = 0;
                detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b);
                uint64_t k;
                for (k = 0; k < b.b.n; k++)
                {
                    g->seq[b.b.a[k]].c = HAP_LABLE;
                }
            }
        }

    }


    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    

    return n_reduced;
}

int untig_asg_arc_simple_large_bubbles(ma_ug_t *ug, asg_t *read_sg, ma_hit_t_alloc* reverse_sources, 
long long miniedgeLen, R_to_U* ruIndex)
{
    asg_t *g = ug->g;
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag, is_hap;
    long long ll, base_maxLen, base_maxLen_i, all_covex;

    buf_t b;
    memset(&b, 0, sizeof(buf_t));
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (nv < 2 || g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;

        base_maxLen = -1;
        base_maxLen_i = -1;
        all_covex = -1;
        is_hap = 0;

        for (i = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                flag = untig_detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, NULL);
                if(flag != TWO_INPUT && flag != MUL_INPUT)
                {
                    break;
                }

                get_real_length(g, convex, &convex);

                if(all_covex != -1 && (uint32_t)all_covex != convex)
                {
                    break;
                }

                if(all_covex == -1)
                {
                    all_covex = convex;
                }
                
                if(base_maxLen < ll)
                {
                    base_maxLen = ll;
                    base_maxLen_i = i;
                }

            }
        }


        if(i == nv)
        {
            for (i = 0; i < nv; i++)
            {
                if(i == base_maxLen_i) continue;
                if (!av[i].del)
                {
                    b.b.n = 0;
                    untig_detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b);
                    

                    //we can only cut tips
                    if(check_if_diploid_untigs(g, read_sg, av[base_maxLen_i].v, av[i].v, &(ug->u),
                    reverse_sources, miniedgeLen, 0.3, &b_0, &b_1, 1, ruIndex) == 1)
                    {
                        n_reduced++;
                        uint64_t k;

                        for (k = 0; k < b.b.n; k++)
                        {
                            g->seq[b.b.a[k]].c = ALTER_LABLE;
                        }
                        for (k = 0; k < b.b.n; k++)
                        {
                            asg_seq_drop(g, b.b.a[k]);
                        }

                        is_hap++;
                    }
                }
            }

            if(is_hap > 0)
            {
                i = base_maxLen_i;
                b.b.n = 0;
                untig_detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b);
                uint64_t k;
                for (k = 0; k < b.b.n; k++)
                {
                    g->seq[b.b.a[k]].c = HAP_LABLE;
                }
            }
        }

    }


    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    free(b_0.b.a);
    free(b_1.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    

    return n_reduced;
}




int asg_arc_cut_long_equal_tips_assembly_complex(asg_t *g, ma_hit_t_alloc* reverse_sources, 
long long miniedgeLen, uint32_t stops_threshold, R_to_U* ruIndex)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag;
    long long ll, max_stopLen;

    buf_t b;
    memset(&b, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {

        uint32_t i;
        ///some node could be deleted
        if (g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        ///tip
        if (get_real_length(g, v, NULL) != 0) continue;
        if(get_real_length(g, v^1, NULL) != 1) continue;
        flag = detect_single_path_with_dels_contigLen(g, v^1, &convex, &ll, NULL);
        if(flag != TWO_INPUT && flag != MUL_INPUT) continue;
        convex = convex^1;
        ///uint32_t n_convex = asg_arc_n(g, convex), convexLen = ll;
        uint32_t n_convex = asg_arc_n(g, convex), convexLen = (uint32_t)-1, convex_i = (uint32_t)-1;
        asg_arc_t *a_convex = asg_arc_a(g, convex);
        for (i = 0; i < n_convex; i++)
        {
            if (!a_convex[i].del)
            {
                detect_single_path_with_dels_contigLen(g, a_convex[i].v, &convex, &ll, NULL);
                if(convex == v)
                {
                    convexLen = ll;
                    convex_i = i;
                    break;
                }
            }
        }

        
        for (i = 0; i < n_convex; i++)
        {
            if (!a_convex[i].del)
            {
                if(i == convex_i) continue;
                
                detect_single_path_with_dels_contigLen_complex(g, a_convex[i].v, &convex, &ll, 
                &max_stopLen, NULL, stops_threshold);
                ///threshold = 0.8
                if(ll > convexLen  && max_stopLen*1.25>ll)
                {
                    ///keep all nodes of this tip
                    b.b.n = 0; 
                    detect_single_path_with_dels_contigLen(g, v^1, &convex, &ll, &b);
                    if(b.b.n < 2) break;
                    b.b.n--;
                    ///keep all nodes of this tip

                     //we can only cut tips
                    if(check_if_diploid_primary_complex(v^1, a_convex[i].v, g, 
                    reverse_sources, miniedgeLen, stops_threshold, 1, ruIndex)==1)
                    {
                        n_reduced++;
                        uint64_t k;

                        for (k = 0; k < b.b.n; k++)
                        {
                            g->seq[b.b.a[k]].c = ALTER_LABLE;
                        }
                        for (k = 0; k < b.b.n; k++)
                        {
                            asg_seq_drop(g, b.b.a[k]);
                        }


                        ///lable the primary one
                        b.b.n = 0; 
                        detect_single_path_with_dels_contigLen_complex(g, a_convex[i].v, &convex, 
                        &ll, &max_stopLen, &b, stops_threshold);
                        for (k = 0; k < b.b.n; k++)
                        {
                            g->seq[b.b.a[k]].c = HAP_LABLE;
                        }

                        break;
                    }

                }
                
            }
        }    
    }


    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    

    return n_reduced;
}



int untig_asg_arc_cut_long_equal_tips_assembly_complex(ma_ug_t *ug, asg_t *read_sg, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen, uint32_t stops_threshold, R_to_U* ruIndex)
{
    asg_t *g = ug->g;
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag;
    long long ll, max_stopLen;

    buf_t buf, b_0, b_1;
    memset(&buf, 0, sizeof(buf_t));
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));


    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i;
        ///some node could be deleted
        if (g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        ///tip
        if (get_real_length(g, v, NULL) != 0) continue;
        if(get_real_length(g, v^1, NULL) != 1) continue;
        /****************************may have bugs********************************/
        buf.b.n = 0; 
        flag = untig_detect_single_path_with_dels_contigLen(g, v^1, &convex, &ll, &buf);
        if(flag != TWO_INPUT && flag != MUL_INPUT) continue;
        get_real_length(g, convex, &convex);
        /****************************may have bugs********************************/
        convex = convex^1;
        uint32_t n_convex = asg_arc_n(g, convex), convexLen = ll;
        asg_arc_t *a_convex = asg_arc_a(g, convex);

        for (i = 0; i < n_convex; i++)
        {
            if (!a_convex[i].del)
            {
                untig_detect_single_path_with_dels_contigLen_complex(g, a_convex[i].v, &convex, &ll, 
                &max_stopLen, NULL, stops_threshold);

                if(convex == v) continue;

                ///threshold = 0.8
                if(ll > convexLen  && max_stopLen*1.25>ll)
                {
                    
                     //we can only cut tips
                    if(check_if_diploid_untigs_complex(g, read_sg, v^1, a_convex[i].v, 
                    &(ug->u), reverse_sources, miniedgeLen, stops_threshold,
                    0.3, &b_0, &b_1, 1, ruIndex)==1)
                    {
                        n_reduced++;
                        uint64_t k;
                        for (k = 0; k < buf.b.n; k++)
                        {
                            g->seq[buf.b.a[k]].c = ALTER_LABLE;
                        }
                        for (k = 0; k < buf.b.n; k++)
                        {
                            asg_seq_drop(g, buf.b.a[k]);
                        }




                        ///lable the primary one
                        b_0.b.n = 0; 
                        untig_detect_single_path_with_dels_contigLen_complex(g, a_convex[i].v, 
                        &convex, &ll, &max_stopLen, &b_0, stops_threshold);
                        for (k = 0; k < b_0.b.n; k++)
                        {
                            g->seq[b_0.b.a[k]].c = HAP_LABLE;
                        }
                        
                        break;
                    }

                }
                
            }
        }    
    }


    asg_cleanup(g);
    asg_symm(g);
    free(buf.b.a);
    free(b_0.b.a);
    free(b_1.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    

    return n_reduced;
}



int untig_asg_arc_cut_chimeric(ma_ug_t *ug, asg_t *read_sg, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen, uint32_t stops_threshold, float drop_rate,
R_to_U* ruIndex)
{
    asg_t *g = ug->g;
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, w1, w2, wv, nw, n_vtx = g->n_seq * 2, n_reduced = 0, convex, convex_T;
    asg_arc_t *aw;
    long long ll, max_stopLen;

    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));


    for (v = 0; v < n_vtx; ++v) 
    {

        uint32_t i;
        ///some node could be deleted
        if (g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        if(get_real_length(g, v, NULL) != 1 || get_real_length(g, v^1, NULL) != 1 ) continue;

        get_real_length(g, v, &w1);
        if(get_real_length(g, w1^1, NULL)<=1) continue;

        get_real_length(g, v^1, &w2);
        if(get_real_length(g, w2^1, NULL)<=1) continue;


        untig_detect_single_path_with_dels_n_stops(g, ug, w1, &convex, &ll, &max_stopLen, NULL, 
        stops_threshold);
        if(ll*drop_rate < EvaluateLen((*ug).u, v>>1)) continue;
        if(ll*0.4 > max_stopLen) continue;


        untig_detect_single_path_with_dels_n_stops(g, ug, w2, &convex, &ll, &max_stopLen, NULL, 
        stops_threshold);
        if(ll*drop_rate < EvaluateLen((*ug).u, v>>1)) continue;
        if(ll*0.4 > max_stopLen) continue;

        stops_threshold++;


        w1 = w1^1;wv=v^1;aw = asg_arc_a(g, w1); nw = asg_arc_n(g, w1);convex_T = (uint32_t)-1;
        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;

            untig_detect_single_path_with_dels_n_stops(g, ug, aw[i].v, &convex, &ll, 
            &max_stopLen, NULL, stops_threshold);
            if(ll*drop_rate < EvaluateLen((*ug).u, v>>1)) break;
            if(ll*0.4 > max_stopLen) continue;

            if((aw[i].v>>1)==(v>>1))
            {
                if(convex_T != (uint32_t)-1) break;
                convex_T = convex;
            }
        }
        if(i!=nw) continue;

        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;
            if((aw[i].v>>1) == (v>>1)) continue;

            untig_detect_single_path_with_dels_n_stops(g, ug, aw[i].v, &convex, &ll, 
            &max_stopLen, NULL, stops_threshold);
            if(convex == convex_T) break;

            if(check_if_diploid_untigs_complex(g, read_sg, wv, aw[i].v, &(ug->u), reverse_sources, 
            miniedgeLen, stops_threshold, 0.3, &b_0, &b_1, 1, ruIndex)==1)
            {
                break;
            }
        }
        if(i!=nw) continue;






        w2 = w2^1;wv=v;aw = asg_arc_a(g, w2); nw = asg_arc_n(g, w2);convex_T = (uint32_t)-1;
        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;

            untig_detect_single_path_with_dels_n_stops(g, ug, aw[i].v, &convex, &ll, 
            &max_stopLen, NULL, stops_threshold);
            if(ll*drop_rate < EvaluateLen((*ug).u, v>>1)) break;
            if(ll*0.4 > max_stopLen) continue;

            if((aw[i].v>>1) == (v>>1))
            {
                if(convex_T != (uint32_t)-1) break;
                convex_T = convex;
            }
        }
        if(i!=nw) continue;

        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;
            if((aw[i].v>>1) == (v>>1)) continue;

            untig_detect_single_path_with_dels_n_stops(g, ug, aw[i].v, &convex, &ll, 
            &max_stopLen, NULL, stops_threshold);
            if(convex == convex_T) break;

            if(check_if_diploid_untigs_complex(g, read_sg, wv, aw[i].v, &(ug->u), reverse_sources, 
            miniedgeLen, stops_threshold, 0.3, &b_0, &b_1, 1, ruIndex)==1)
            {
                break;
            }
        }
        if(i!=nw) continue;






        n_reduced++;g->seq[v>>1].c = ALTER_LABLE;asg_seq_drop(g, v>>1);
        ///fprintf(stderr, "v>>1: %u\n", v>>1);
    }


    asg_cleanup(g);
    free(b_0.b.a);
    free(b_1.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    

    return n_reduced;
}




void output_unitig_graph(asg_t *sg, ma_sub_t* coverage_cut, char* output_file_name, long long n_read)
{
    ma_ug_t *ug = NULL;
    ug = ma_ug_gen(sg);
    ma_ug_seq(ug, &R_INF, coverage_cut, n_read);

    fprintf(stderr, "Writing raw unitig GFA to disk... \n");
    char* gfa_name = (char*)malloc(strlen(output_file_name)+25);
    sprintf(gfa_name, "%s.r_utg.gfa", output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    ma_ug_print(ug, &R_INF, coverage_cut, output_file);
    fclose(output_file);
    sprintf(gfa_name, "%s.r_utg.noseq.gfa", output_file_name);
    output_file = fopen(gfa_name, "w");
    ma_ug_print_simple(ug, &R_INF, coverage_cut, output_file);
    fclose(output_file);

    free(gfa_name);
    ma_ug_destroy(ug);
}


void output_read_graph(asg_t *sg, ma_sub_t* coverage_cut, char* output_file_name, long long n_read)
{
    fprintf(stderr, "Writing read GFA to disk... \n");
    char* gfa_name = (char*)malloc(strlen(output_file_name)+25);
    sprintf(gfa_name, "%s.read.gfa", output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    ma_sg_print(sg, &R_INF, coverage_cut, output_file);
    free(gfa_name);
    fclose(output_file);
}


void read_ma(ma_hit_t* x, FILE* fp)
{
    int f_flag;
    f_flag = fread(&(x->qns), sizeof(x->qns), 1, fp);
    f_flag += fread(&(x->qe), sizeof(x->qe), 1, fp);
    f_flag += fread(&(x->tn), sizeof(x->tn), 1, fp);
    f_flag += fread(&(x->ts), sizeof(x->ts), 1, fp);
    f_flag += fread(&(x->te), sizeof(x->te), 1, fp);
    f_flag += fread(&(x->el), sizeof(x->el), 1, fp);
    f_flag += fread(&(x->no_l_indel), sizeof(x->no_l_indel), 1, fp);

    uint32_t t;
    f_flag += fread(&(t), sizeof(t), 1, fp);
    x->ml = t;

    f_flag += fread(&(t), sizeof(t), 1, fp);
    x->rev = t;
    
    f_flag += fread(&(t), sizeof(t), 1, fp);
    x->bl = t;

    f_flag += fread(&(t), sizeof(t), 1, fp);
    x->del = t;
}

int load_ma_hit_ts(ma_hit_t_alloc** x, char* read_file_name)
{
    fprintf(stderr, "Loading ma_hit_ts from disk... \n");
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "r");
    if(!fp)
    {
        return 0;
    }


    long long n_read;
    long long i, k;
    int f_flag;
    f_flag = fread(&n_read, sizeof(n_read), 1, fp);
    (*x) = (ma_hit_t_alloc*)malloc(sizeof(ma_hit_t_alloc)*n_read);


    for (i = 0; i < n_read; i++)
    {        
        f_flag += fread(&((*x)[i].is_fully_corrected), sizeof((*x)[i].is_fully_corrected), 1, fp);
        f_flag += fread(&((*x)[i].is_abnormal), sizeof((*x)[i].is_abnormal), 1, fp);
        f_flag += fread(&((*x)[i].length), sizeof((*x)[i].length), 1, fp);
        (*x)[i].size = (*x)[i].length;

        (*x)[i].buffer = (ma_hit_t*)malloc(sizeof(ma_hit_t)*(*x)[i].length);

        for (k = 0; k < (*x)[i].length; k++)
        {
            read_ma(&((*x)[i].buffer[k]), fp);
        }  
    }

    free(index_name);    
    fclose(fp);
    fprintf(stderr, "ma_hit_ts has been read.\n");

    return 1;
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
    fprintf(stderr, "Writing ma_hit_ts to disk... \n");
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
    fprintf(stderr, "ma_hit_ts has been written.\n");
}

void write_all_data_to_disk(ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, 
All_reads *RNF, char* output_file_name)
{   
    char* gfa_name = (char*)malloc(strlen(output_file_name)+25);
    sprintf(gfa_name, "%s.ovlp", output_file_name);
    write_All_reads(RNF, gfa_name);
    
    sprintf(gfa_name, "%s.ovlp.source", output_file_name);
    write_ma_hit_ts(sources, RNF->total_reads, gfa_name);

    sprintf(gfa_name, "%s.ovlp.reverse", output_file_name);
    write_ma_hit_ts(reverse_sources, RNF->total_reads, gfa_name);

    free(gfa_name);
}

int load_all_data_from_disk(ma_hit_t_alloc **sources, ma_hit_t_alloc **reverse_sources, 
char* output_file_name)
{
    char* gfa_name = (char*)malloc(strlen(output_file_name)+25);
    sprintf(gfa_name, "%s.ovlp", output_file_name);
    if(!load_All_reads(&R_INF, gfa_name))
    {
        return 0;
    }

    sprintf(gfa_name, "%s.ovlp.source", output_file_name);
    if(!load_ma_hit_ts(sources, gfa_name))
    {
        return 0;
    }


    sprintf(gfa_name, "%s.ovlp.reverse", output_file_name);
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
			if (d + l > (uint32_t)max_dist) break; // too far

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
	if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] popped %lu bubbles\n", __func__, (unsigned long)n_pop);
    }
    return n_pop;
}


// in a resolved bubble, mark unused vertices and arcs as "reduced"
static void asg_bub_backtrack_primary(asg_t *g, uint32_t v0, buf_t *b)
{
	uint32_t i, v, qn, tn;

	///assert(b->S.n == 1);
	///first remove all nodes in this bubble
	for (i = 0; i < b->b.n; ++i)
    {
        g->seq[b->b.a[i]>>1].c = ALTER_LABLE;
    }

    ///v is the sink of this bubble
	v = b->S.a[0];
	///recover node
	do {
		uint32_t u = b->a[v].p; // u->v
        /****************************may have hap bugs********************************/
		////g->seq[v>>1].c = PRIMARY_LABLE;
        g->seq[v>>1].c = HAP_LABLE;
        /****************************may have hap bugs********************************/
		v = u;
	} while (v != v0);


    ///remove all edges (self/reverse for each edge) in this bubble
	for (i = 0; i < b->e.n; ++i) {
		asg_arc_t *a = &g->arc[b->e.a[i]];
        qn = a->ul>>33;
        tn = a->v>>1;
        if(g->seq[qn].c == ALTER_LABLE && g->seq[tn].c == ALTER_LABLE)
        {
            continue;
        }
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

// pop bubbles from vertex v0; the graph MJUST BE symmetric: if u->v present, v'->u' must be present as well
static uint64_t asg_bub_pop1_primary(asg_t *g, uint32_t v0, int max_dist, buf_t *b)
{
	uint32_t i, n_pending = 0, is_first = 1;
	uint64_t n_pop = 0;
	///if this node has been deleted
	if (g->seq[v0>>1].del || g->seq[v0>>1].c == ALTER_LABLE) return 0; // already deleted
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
            /****************************may have bugs********************************/
            ///important when poping at long untig graph
            if(is_first) l = 0;
            /****************************may have bugs********************************/

			///if this edge has been deleted
			if (av[i].del) continue;

			///push the edge
            ///high 32-bit of g->idx[v] is the start point of v's edges
            //so here is the point of this specfic edge
			kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);
			///find a too far path? directly terminate the whole bubble poping
			if (d + l > (uint32_t)max_dist) break; // too far

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
        is_first = 0;
		///if i < nv, that means (d + l > max_dist)
		if (i < nv || b->S.n == 0) goto pop_reset;
	} while (b->S.n > 1 || n_pending);
	asg_bub_backtrack_primary(g, v0, b);
    n_pop = 1;
pop_reset:
	for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
		binfo_t *t = &b->a[b->b.a[i]];
		t->s = t->c = t->d = 0;
	}
	return n_pop;
}


// pop bubbles
int asg_pop_bubble_primary(asg_t *g, int max_dist)
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
		if (nv < 2 || g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
		///some edges could be deleted
		for (i = 0; i < nv; ++i) // asg_bub_pop1() may delete some edges/arcs
			if (!av[i].del) ++n_arc;
		if (n_arc > 1)
			n_pop += asg_bub_pop1_primary(g, v, max_dist, &b);
	}
	free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
	if (n_pop) asg_cleanup(g);

	if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] popped %lu bubbles\n", __func__, (unsigned long)n_pop);
    }
    return n_pop;
}





int test_triangular_directly(asg_t *g, uint32_t v, 
long long min_edge_length, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex)
{
    
    uint32_t w;
    int todel;
    long long NodeLen_first[3];
    long long NodeLen_second[3];
    
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
    if(check_if_diploid(av[0].v, av[1].v, g, reverse_sources, min_edge_length, ruIndex) == 1||
    check_if_diploid(aw[0].v, aw[1].v, g, reverse_sources, min_edge_length, ruIndex) == 1)
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



int asg_arc_del_triangular_directly(asg_t *g, long long min_edge_length, 
ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0;

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t nv = asg_arc_n(g, v);
        if (g->seq[v>>1].del)
        {
            continue;
        } 

        if(nv < 2)
        {
            continue;
        }


        n_reduced += test_triangular_directly(g, v, min_edge_length, reverse_sources, ruIndex);
    }

    

    if (n_reduced) {
        asg_cleanup(g);
        asg_symm(g);
    }

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d triangular overlaps\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    return n_reduced;
}




int asg_arc_del_orthology(asg_t *g, ma_hit_t_alloc* reverse_sources, float drop_ratio, 
long long miniedgeLen, R_to_U* ruIndex)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0;
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

        if(check_if_diploid(av[idx[0]].v, av[idx[1]].v, g, reverse_sources, miniedgeLen, ruIndex) == 0)
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




int asg_arc_del_orthology_multiple_way(asg_t *g, ma_hit_t_alloc* reverse_sources, float drop_ratio, 
long long miniedgeLen, R_to_U* ruIndex)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, v_max, v_maxLen, n_vtx = g->n_seq * 2, n_reduced = 0;

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
        v_maxLen = 0;

        for (i = 0, n_arc = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                if(v_max == (uint32_t)-1)
                {
                    v_max = av[i].v;
                    v_maxLen = av[i].ol;
                }
                else if(check_if_diploid(v_max, av[i].v, g, reverse_sources, miniedgeLen, ruIndex) == 0)
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

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d different hap overlaps\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

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

        ///if(b) kv_push(uint32_t, b->b, v>>1);
        if(b) kv_push(uint32_t, b->b, v);

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
            ///if(b) kv_push(uint32_t, b->b, v>>1);
            if(b) kv_push(uint32_t, b->b, v);
            return TWO_INPUT;
        }

        if(kw > 2)
        {
            (*Len)++;
            ///if(b) kv_push(uint32_t, b->b, v>>1);
            if(b) kv_push(uint32_t, b->b, v);
            return MUL_INPUT;
        }   


        if((v>>1) == (begNode>>1))
        {
            return LOOP;
        } 
    }

    return LONG_TIPS;   
}



long long asg_arc_del_self_circle_untig(asg_t *g, long long circleLen, int is_drop)
{
    uint32_t v, w, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag;
    long long ll;
    asg_arc_t *aw;
    uint32_t nw, k;
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (g->seq[v>>1].del) continue;
        if(is_drop && g->seq[v>>1].c == ALTER_LABLE) continue;

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

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d self-circles\n", 
        __func__, n_reduced);
    }

    return n_reduced;
}


long long get_untig_coverage(ma_hit_t_alloc* sources, ma_sub_t* coverage_cut, uint32_t* b, uint64_t n)
{
    uint64_t k, j;
    uint32_t v;
    ma_hit_t *h;
    long long R_bases = 0, C_bases = 0;
    for (k = 0; k < n; ++k) 
    {
        v = b[k]>>1;
        R_bases += coverage_cut[v].e - coverage_cut[v].s;
        for (j = 0; j < (uint64_t)(sources[v].length); j++)
        {
            h = &(sources[v].buffer[j]);
            C_bases += Get_qe((*h)) - Get_qs((*h));
        }
    }
    return C_bases/R_bases;
}

/**
void copy_untig(asg_t *g, long long times, uint32_t* b, long long n, C_graph* cg)
{
    if(times <= 1) return;

    times--;
    long long i, j;
    uint64_t tmp;
    for (i = 0; i < times; i++)
    {
        for (j = 0; j < n; j++)
        {
            tmp = 
            kv_push(uint64_t, cg->Node, );
            asg_add_auxiliary_seq_set(g, (b[j]>>1), 0);
        }
    }
}
**/

void init_C_graph(C_graph* g, uint32_t n_seq)
{
    kv_init(g->Nodes);
    kv_init(g->Edges);
    g->pre_n_seq = n_seq;
    g->seqID = n_seq;
}

void destory_C_graph(C_graph* g)
{
    kv_destroy(g->Nodes);
    kv_destroy(g->Edges);
}


long long asg_arc_del_simple_circle_untig(ma_hit_t_alloc* sources, ma_sub_t* coverage_cut, asg_t *g, long long circleLen, int is_drop)
{
    uint32_t v, w, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag;
    long long ll/**, coverage**/;
    asg_arc_t *aw;
    uint32_t nw, k;
    buf_t b;
    memset(&b, 0, sizeof(buf_t));

    C_graph cg;
    init_C_graph(&cg, g->n_seq);

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (g->seq[v>>1].del) continue;
        if(is_drop && g->seq[v>>1].c == ALTER_LABLE) continue;
        if(get_real_length(g, v^1, NULL)<=1) continue;

        n_arc = get_real_length(g, v, NULL);
        if (n_arc != 1) continue;

        for (i = 0; i < nv; i++)
        {
            ///actually there is just one un-del edge
            if (!av[i].del)
            {
                b.b.n = 0;    
                flag = detect_single_path_with_dels_by_length(g, v, &convex, &ll, &b, circleLen);
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
                    if(v == convex && b.b.n > 1)
                    {
                        convex = b.b.a[b.b.n - 2];
                        b.b.n--;
                    }

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
                            // coverage = get_untig_coverage(sources, coverage_cut, b.b.a, b.b.n);
                            // copy_untig(g, (coverage/asm_opt.coverage), b.b.a, b.b.n, &cg);

                            
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

    free(b.b.a);
    destory_C_graph(&cg);

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d self-circles\n", 
        __func__, n_reduced);
    }

    return n_reduced;
}


void output_unitig_graph_without_small_bubbles_primary(asg_t *sg, ma_sub_t* coverage_cut, 
char* output_file_name, long long n_read, long long bubble_dist, long long tipsLen)
{
    asg_cut_tip_primary(sg, NULL, tipsLen);
    asg_pop_bubble_primary(sg, bubble_dist);
    asg_cut_tip_primary(sg, NULL, tipsLen);
    
    ma_ug_t *ug = NULL;
    ug = ma_ug_gen_primary(sg, PRIMARY_LABLE);
    ma_ug_seq(ug, &R_INF, coverage_cut, n_read);

    fprintf(stderr, "Writing processed unitig GFA to disk... \n");
    char* gfa_name = (char*)malloc(strlen(output_file_name)+35);
    sprintf(gfa_name, "%s.p_utg.gfa", output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    ma_ug_print(ug, &R_INF, coverage_cut, output_file);
    fclose(output_file);
    
    sprintf(gfa_name, "%s.p_utg.noseq.gfa", output_file_name);
    output_file = fopen(gfa_name, "w");
    ma_ug_print_simple(ug, &R_INF, coverage_cut, output_file);
    fclose(output_file);

    free(gfa_name);
    ma_ug_destroy(ug);
}



long long get_graph_statistic(asg_t *g)
{
    long long num_arc = 0;
	uint32_t n_vtx = g->n_seq * 2, v;
    
	for (v = 0; v < n_vtx; ++v) 
    {
		if (g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        num_arc += asg_arc_n(g, v);
    }

    return num_arc;
}

uint32_t* build_unitig_index(asg_t *sg, ma_ug_t *ug)
{
    if(sg == NULL && ug == NULL) return NULL;
    uint32_t* index = (uint32_t*)calloc(sg->n_seq, sizeof(uint32_t));
    uint32_t i, j;
    for (i = 0; i < ug->u.n; ++i) 
    {
        ma_utg_t *u = &(ug->u.a[i]);
        for (j = 0; j < u->n; j++)
        {
            index[(u->a[j]>>33)] = i;
        }
    }
    return index;
}




int explore_graph(asg_t *nsg, uint32_t vBeg, float single_threshold,
float l_untig_rate_threshold, long long minLongUntig, long long maxShortUntig, uint32_t ignore_d,
buf_t* bb, uint32_t** r, size_t* rm, size_t* rn, uint8_t* visit, ma_utg_v* ut_v, uint32_t* r_ID)
{
    (*r_ID) = (uint32_t)-1;
    uint32_t nv;
    asg_arc_t *av;

    kvec_t(uint32_t) u_vecs;
    kv_init(u_vecs);
    if(r && rm && rn) kv_reuse(u_vecs, 0, (*rm), (*r));
    
    memset(visit, 0, nsg->n_seq);
    kdq_t(uint32_t) *buf;
    buf = kdq_init(uint32_t);

    uint32_t vEnd, threshold, num_reads = 0, i, k, v, end = (uint32_t)-1, in = 0, vELen, tmp;
    bb->b.n = 0;
    threshold = get_long_tip_length(nsg, ut_v, vBeg, &vEnd, bb);
    for (i = 0; i < bb->b.n; i++)
    {
        Set_vis(visit, bb->b.a[i], ignore_d);
        Set_vis(visit, bb->b.a[i]^1, ignore_d);
    }    
    v = vBeg;


    kdq_push(uint32_t, buf, v);
    while (kdq_size(buf) != 0)
    {
        in++;
        v = *(kdq_pop(uint32_t, buf));

        ///in == 1 means the start node, it is useless
        if(in != 1)
        {
            ///get current untig length
            bb->b.n = 0;
            vELen = get_long_tip_length(nsg, ut_v, v, &vEnd, bb);

            ///first long unitig except the start node
            if(if_long_tip_length(nsg, ut_v, v, &vELen, 
                    minLongUntig, maxShortUntig, l_untig_rate_threshold, threshold)==1)
            {
                if(end == (uint32_t)-1)
                {
                    end = v;
                    continue;
                }
                else
                {
                    in = (uint32_t)-1;
                    break;
                }
            }

            if(vELen > (threshold * single_threshold))
            {
                in = (uint32_t)-1;
                break;
            }

            num_reads = num_reads + vELen;
            if(num_reads > threshold)
            {
                in = (uint32_t)-1;
                break;
            }

            if(r && rm && rn)
            {
                for (i = 0; i < bb->b.n; i++)
                {
                    kv_push(uint32_t, u_vecs, bb->b.a[i]);
                }   
            } 
        }

        tmp = v^1;
        v = vEnd;
        nv = asg_arc_n(nsg, v);
        av = asg_arc_a(nsg, v);
        for(k = 0; k < nv; k++)
        {
            if(av[k].del) continue;

            if(av[k].v == vBeg)
            {
                in = (uint32_t)-1;
                goto termi;
            }
            
            if(Get_vis(visit,av[k].v,ignore_d)==0)
            {
                kdq_push(uint32_t, buf, av[k].v);


                bb->b.n = 0;
                get_long_tip_length(nsg, ut_v, av[k].v, &vEnd, bb);
                for (i = 0; i < bb->b.n; i++)
                {
                    Set_vis(visit, bb->b.a[i], ignore_d);
                }
            }
        }

        ///for start node, we just need one direction
        if(in != 1 && ignore_d)
        {
            v = tmp;
            nv = asg_arc_n(nsg, v);
            av = asg_arc_a(nsg, v);
            for(k = 0; k < nv; k++)
            {
                if(av[k].del) continue;

                if(av[k].v == vBeg)
                {
                    in = (uint32_t)-1;
                    goto termi;
                }

                if(Get_vis(visit,av[k].v,ignore_d)==0)
                {
                    kdq_push(uint32_t, buf, av[k].v);


                    bb->b.n = 0;
                    get_long_tip_length(nsg, ut_v, av[k].v, &vEnd, bb);
                    for (i = 0; i < bb->b.n; i++)
                    {
                        Set_vis(visit, bb->b.a[i], ignore_d);
                    }
                }
            }
        }
    }

    termi:
    kdq_destroy(uint32_t, buf);
    if(r && rm && rn)
    {
        (*rn) = u_vecs.n;
        (*rm) = u_vecs.m;
        (*r) = u_vecs.a;
    }
    

    (*r_ID) = end;
    ///in == 1 means the end subgraph is the long untig itself
    ///so here is no tangles
    if(in == (uint32_t)-1 || in == 1)
    {
        return 0;
    }
    else
    {
        return 1;
    }
    
}


int explore_graph_back(asg_t *nsg, uint32_t start, uint32_t threshold, float single_threshold,
float l_untig_rate_threshold, long long minLongUntig, long long maxShortUntig, uint32_t ignore_d,
uint32_t** r, size_t* rm, size_t* rn, uint8_t* visit, ma_utg_v* ut_v, uint32_t* r_ID)
{
    (*r_ID) = (uint32_t)-1;
    uint32_t nv;
    asg_arc_t *av;

    kvec_t(uint32_t) u_vecs;
    kv_init(u_vecs);
    if(r && rm && rn) kv_reuse(u_vecs, 0, (*rm), (*r));
    
    memset(visit, 0, nsg->n_seq);
    kdq_t(uint32_t) *buf;
    buf = kdq_init(uint32_t);

    uint32_t num_reads = 0;
    uint32_t v = start, k;
    uint32_t end = (uint32_t)-1;
    uint32_t in = 0;
    Set_vis(visit, v, ignore_d);
    Set_vis(visit, v^1, ignore_d);

    if(nsg->seq[v>>1].del)
    {
        in = (uint32_t)-1;
        goto termi;
    }

    kdq_push(uint32_t, buf, v);
    while (kdq_size(buf) != 0)
    {
        in++;
        v = *(kdq_pop(uint32_t, buf));

        ///in == 1 means the start node, it is useless
        if(in != 1)
        {
            ///first long unitig except the start node
            if(check_long_tip(*ut_v, v>>1, minLongUntig, 
                maxShortUntig, l_untig_rate_threshold, threshold))
            {
                if(end == (uint32_t)-1)
                {
                    end = v;
                    continue;
                }
                else
                {
                    in = (uint32_t)-1;
                    break;
                }
            }

            if(EvaluateLen(*ut_v, v>>1) > (threshold * single_threshold))
            {
                in = (uint32_t)-1;
                break;
            }

            num_reads = num_reads + EvaluateLen(*ut_v, v>>1);
            if(num_reads > threshold)
            {
                in = (uint32_t)-1;
                break;
            }

            if(r && rm && rn) kv_push(uint32_t, u_vecs, v);
        }

        
        
        nv = asg_arc_n(nsg, v);
        av = asg_arc_a(nsg, v);
        for(k = 0; k < nv; k++)
        {
            if(av[k].del) continue;

            if(av[k].v == start)
            {
                in = (uint32_t)-1;
                goto termi;
            }
            
            if(Get_vis(visit,av[k].v,ignore_d)==0)
            {
                Set_vis(visit,av[k].v,ignore_d);
                kdq_push(uint32_t, buf, av[k].v);
            }
        }

        ///for start node, we just need one direction
        if(in != 1 && ignore_d)
        {
            v = v^1;
            nv = asg_arc_n(nsg, v);
            av = asg_arc_a(nsg, v);
            for(k = 0; k < nv; k++)
            {
                if(av[k].del) continue;

                if(av[k].v == start)
                {
                    in = (uint32_t)-1;
                    goto termi;
                }

                if(Get_vis(visit,av[k].v,ignore_d)==0)
                {
                    Set_vis(visit,av[k].v,ignore_d);
                    kdq_push(uint32_t, buf, av[k].v);
                }
            }
        }
    }

    termi:
    kdq_destroy(uint32_t, buf);
    if(r && rm && rn)
    {
        (*rn) = u_vecs.n;
        (*rm) = u_vecs.m;
        (*r) = u_vecs.a;
    }
    

    (*r_ID) = end;

    if(in == (uint32_t)-1)
    {
        return 0;
    }
    else
    {
        return 1;
    }
    
}



void output_tangles(uint32_t startID, uint32_t endId, uint32_t* a, uint32_t n, const char* lable)
{
    kvec_t(uint32_t) u_vecs;
    u_vecs.a = a, u_vecs.n = n;

    fprintf(stderr, "\n%sstartID: %u, dir: %u\n", lable, startID>>1, startID&1);
    fprintf(stderr, "%sendID: %u, dir: %u\n", lable, endId>>1, endId&1);
    uint32_t ijk;
    for (ijk = 0; ijk < u_vecs.n; ijk++)
    {
        fprintf(stderr, "%stangleID: %u, dir: %u\n", lable, u_vecs.a[ijk]>>1, u_vecs.a[ijk]&1);
    }
}


int get_arc(asg_t *g, uint32_t src, uint32_t dest, asg_arc_t* result)
{
    uint32_t i;
    if(g->seq[src>>1].del) return 0;

    uint32_t nv = asg_arc_n(g, src);
    asg_arc_t *av = asg_arc_a(g, src);
    for (i = 0; i < nv; i++)
    {
        if(av[i].del) continue;
        if(av[i].v == dest)
        {
            (*result) = av[i];
            break;
        } 
    }

    if(i != nv) return 1;
    return 0;
}

uint32_t insert_index(asg_t *g, uint32_t v)
{
    uint32_t v_tx = g->n_seq * 2;
    if(v >= v_tx) return (uint32_t)-1;
    if(asg_arc_n(g, v)!=0) return (g->idx[v]>>32) + asg_arc_n(g, v);
    ///now v itself does not have any edge
    while (v < v_tx && asg_arc_n(g, v) == 0){v++;}
    ///means there are no edge at the whole graph
    if(v>=v_tx) return g->n_arc;
    return (g->idx[v]>>32);
}

asg_arc_t* insert_index_p(asg_t *g, long long index, uint32_t m_distance)
{
    ///each edge has two direction
    if (g->n_arc + m_distance > g->m_arc) 
    {
        ///g->m_arc = g->n_arc + (m_distance<<1);
        g->m_arc = (g->n_arc + m_distance)<<1;
		g->arc = (asg_arc_t*)realloc(g->arc, g->m_arc * sizeof(asg_arc_t));
	}
    long long i = g->n_arc; i--;
    for (; i >= index; i--)
    {
        g->arc[i+m_distance] = g->arc[i];
    }

    g->n_arc = g->n_arc + m_distance;

    return &(g->arc[index]);
}

void exchange_arcs(asg_arc_t* x, asg_arc_t* y)
{
    asg_arc_t k;
    k = (*x);(*x) = (*y);(*y) = k;
}

void insert_arc(asg_t *g, long long index, uint32_t m_distance, uint32_t src, uint32_t srcLen, uint32_t dest, 
uint32_t oLen, uint8_t strong, uint8_t el, uint8_t no_l_indel)
{
    asg_arc_t* p;
    uint32_t l;
    long long i;
    p = insert_index_p(g, index, m_distance);
    p->del = !!(0);p->el = el; p->no_l_indel = no_l_indel; p->strong = strong; p->ol = oLen;
    p->v = dest;
    p->ul = src; p->ul = p->ul << 32; l = srcLen - oLen; p->ul = p->ul | l;
    for (i = index - 1; i >= 0 && p->ul <= g->arc[i].ul; i--)
    {
        if((!g->arc[i].del)&&(g->arc[i].v==p->v) && ((g->arc[i].ul>>32)==(p->ul>>32)))
        {
            p->del = !!(1);
            break;
        }
        exchange_arcs(p, &(g->arc[i]));p = &(g->arc[i]);
    }
    

    if(asg_arc_n(g, src) == 0)
    {
        g->idx[src] = index; g->idx[src] = g->idx[src]<<32; g->idx[src] = g->idx[src] | 1;
    }
    else
    {
        g->idx[src] = g->idx[src] + 1;
    }
      
    uint32_t v, v_tx = g->n_seq*2;
    uint64_t add = 1; add = add << 32;
    for (v = src+1; v < v_tx; v++)
    {
        if(asg_arc_n(g, v) == 0) continue;
        g->idx[v] += add;
    }
}

int asg_append_edges_to_srt(asg_t *g, uint32_t src, 
uint32_t srcLen, uint32_t dest, uint32_t oLen,
uint8_t strong, uint8_t el, uint8_t no_l_indel)
{
    uint32_t v_tx = g->n_seq * 2;
    /**
    uint32_t d_i = 0, current_i, next_i, src_n, s_index, e_index;
    while(d_i < v_tx)
    {
        current_i = d_i;
        d_i++;
        src_n = asg_arc_n(g, current_i);
        if(src_n == 0) continue;
        s_index = (g->idx[current_i]>>32);
        e_index = s_index + src_n;

        while (d_i < v_tx && asg_arc_n(g, d_i) == 0){d_i++;}
        if(d_i>=v_tx) break;
        next_i = d_i;


        if(current_i != v_tx)
        {
            if(e_index != (g->idx[next_i]>>32))
            {
                fprintf(stderr, "v_tx: %u, current_i: %u, node: %u, s_index: %u, e_index: %u, g->idx[next_i]>>32: %u\n",
                v_tx, current_i, current_i>>1, s_index, e_index, g->idx[next_i]>>32);
            }
        }
        else 
        {
            if(e_index != g->n_arc)
            {
                fprintf(stderr, "ERROR2\n");
            }
        }
    }

    for (d_i = 0; d_i < v_tx; d_i++)
    {
        uint32_t nv = asg_arc_n(g, d_i), s_i;
        asg_arc_t *av = asg_arc_a(g, d_i);
        asg_arc_t forward, backward;
        for (s_i = 0; s_i < nv; s_i++)
        {
            if(av[s_i].del) continue;
            forward = av[s_i];
            if(get_arc(g, forward.ul>>32, forward.v, &forward) == 0)
            {
                fprintf(stderr, "ERROR1\n");
            }
            if(forward.del != av[s_i].del || forward.el != av[s_i].el || 
               forward.no_l_indel != av[s_i].no_l_indel ||forward.ol != av[s_i].ol || 
               forward.strong != av[s_i].strong || forward.ul != av[s_i].ul || 
               forward.v != av[s_i].v)
            {
                fprintf(stderr, "ERROR2\n");
            }

            if(get_arc(g, forward.v^1, (forward.ul>>32)^1, &backward) == 0)
            {
                fprintf(stderr, "ERROR3\n");
            }

            if(forward.ol != backward.ol)
            {
                fprintf(stderr, "forward.ol: %u, backward.ol: %u\n", 
                forward.ol, backward.ol);
            }

        }
    }
    **/


    if (src >= v_tx || dest >= v_tx)
    {
        return 0;
    }

    ///we need to link src--->dest and (dest^1)----->(src^1)
    uint32_t src_i = insert_index(g, src);
    insert_arc(g, src_i, 1, src, srcLen, dest, oLen, strong, el, no_l_indel);


    /**
    if (src >= v_tx || (src^1) >= v_tx || dest >= v_tx || (dest^1) >= v_tx)
    {
        return 0;
    }

    ///we need to link src--->dest and (dest^1)----->(src^1)
    uint32_t src_i = insert_index(g, src);
    insert_arc(g, src_i, 1, src, srcLen, dest, oLen);
    uint32_t dest_i = insert_index(g, (dest^1));
    insert_arc(g, dest_i, 1, dest^1, destLen, src^1, oLen);
    **/
    
    
    ///we need to link src--->dest and (dest^1)----->(src^1)
    /**
    if(src > (dest^1))
    {
        uint32_t src_i = insert_index(g, src);
        insert_arc(g, src_i, src, srcLen, dest, oLen);
        uint32_t dest_i = insert_index(g, (dest^1));
        insert_arc(g, dest_i, dest^1, destLen, src^1, oLen);
    }
    else
    {
        uint32_t dest_i = insert_index(g, (dest^1));
        insert_arc(g, dest_i, dest^1, destLen, src^1, oLen);
        uint32_t src_i = insert_index(g, src);
        insert_arc(g, src_i, src, srcLen, dest, oLen);
    }
    **/

   return 1;
}

void append_ma_utg_t(ma_utg_t* v_x, ma_utg_t* v_y)
{
    if(v_x->m < (v_x->n + v_y->n))
    {
        v_x->m = v_x->n + v_y->n;
        v_x->a = (uint64_t*)realloc(v_x->a, v_x->m * sizeof(uint64_t));
    }
    memcpy(v_x->a+v_x->n, v_y->a, v_y->n*sizeof(uint64_t));
    v_x->n = v_x->n + v_y->n;
    free(v_y->a);
    v_y->m=v_y->n=0;v_y->a=NULL;
}


#define UNROLL 0
#define CONVEX 1
void merge_nodes(ma_ug_t *ug, uint32_t startID, uint32_t endId, uint32_t* a, uint32_t n, 
uint32_t type)
{
    ma_utg_t *v_x = NULL, *v_y = NULL;
    asg_t* nsg = ug->g;
    kvec_t(uint32_t) u_vecs;
    uint32_t i = 0, maxEvaluateLen = 0, maxBaseLen = 0, v, totalEvaluateLen = 0;
    u_vecs.a = a, u_vecs.n = n;
    if(u_vecs.n > 0) 
    {
        v_x = &(ug->u.a[u_vecs.a[0]>>1]);
        asg_seq_del(nsg, u_vecs.a[0]>>1);
        maxEvaluateLen = EvaluateLen(ug->u, u_vecs.a[0]>>1);
        maxBaseLen = v_x->len;
        totalEvaluateLen += EvaluateLen(ug->u, u_vecs.a[0]>>1);
    }

    
    for (i = 1; i < u_vecs.n; i++)
    {
        v_y = &(ug->u.a[u_vecs.a[i]>>1]);
        if(maxEvaluateLen <= EvaluateLen(ug->u, u_vecs.a[i]>>1))
        {
            maxEvaluateLen = EvaluateLen(ug->u, u_vecs.a[i]>>1);
            maxBaseLen = v_y->len;
        }
        totalEvaluateLen += EvaluateLen(ug->u, u_vecs.a[i]>>1);

        
        append_ma_utg_t(v_x, v_y);
        asg_seq_del(nsg, u_vecs.a[i]>>1);
    }


    if(u_vecs.n > 0)
    {
        i = 0;
        nsg->seq[u_vecs.a[0]>>1].del = !!(0);


        EvaluateLen(ug->u, u_vecs.a[0]>>1) = maxEvaluateLen;
        totalEvaluateLen = totalEvaluateLen / 2;
        if(totalEvaluateLen > EvaluateLen(ug->u, u_vecs.a[0]>>1))
        {
            EvaluateLen(ug->u, u_vecs.a[0]>>1) = totalEvaluateLen;
        }
        IsMerge(ug->u, u_vecs.a[0]>>1)++;


        v_x->len = maxBaseLen;
        /****************************may have bugs********************************/
        nsg->seq[u_vecs.a[0]>>1].len = maxBaseLen;
        /****************************may have bugs********************************/
        v = u_vecs.a[0];

        if(type == UNROLL)
        {
            if(startID != (uint32_t)-1)
            {
                asg_append_edges_to_srt(nsg, startID, ug->u.a[startID>>1].len, v, 0, 0, 0 ,0);
                asg_append_edges_to_srt(nsg, v^1, ug->u.a[v>>1].len, startID^1, 0, 0, 0, 0);
            }
            
            if(endId != (uint32_t)-1)
            {
                asg_append_edges_to_srt(nsg, v, ug->u.a[v>>1].len, endId, 0, 0, 0 ,0);
                asg_append_edges_to_srt(nsg, endId^1, ug->u.a[endId>>1].len, v^1, 0, 0, 0 ,0);
            }
        }

        if(type == CONVEX)
        {
            if(startID != (uint32_t)-1)
            {
                asg_append_edges_to_srt(nsg, v, ug->u.a[v>>1].len, startID^1, 0, 0, 0, 0);
                asg_append_edges_to_srt(nsg, startID, ug->u.a[startID>>1].len, v^1, 0, 0, 0 ,0);
            }

            if(endId != (uint32_t)-1)
            {
                asg_append_edges_to_srt(nsg, v, ug->u.a[v>>1].len, endId, 0, 0, 0 ,0);
                asg_append_edges_to_srt(nsg, endId^1, ug->u.a[endId>>1].len, v^1, 0, 0, 0 ,0);
            }
        }
        
    }
}

///return 1 is what we want
///as for return value: 0: do nothing, 1: unroll, 2: convex
#define CONVEX_M 1
#define UNROLL_M 2
#define UNROLL_E 3
inline uint32_t walk_through(asg_t *sg, ma_ug_t *ug, ma_hit_t_alloc* reverse_sources, long long minLongUntig, 
long long maxShortUntig, float l_untig_rate, float max_node_threshold, buf_t* b_0, buf_t* b_1,
kvec_t_u32_warp* u_vecs, uint8_t* visit, uint32_t v, uint32_t* r_beg, uint32_t* r_end, 
uint32_t* r_next_uID, R_to_U* ruIndex)
{
    (*r_beg) = (*r_end) = (uint32_t)-1;
    asg_t* nsg = ug->g;
    uint32_t i, beg, end, primaryLen, returnFlag;

    /*****************************simple checking**********************************/
    beg = v;
    if(nsg->seq[beg>>1].del || asg_arc_n(nsg, beg) <= 0 || get_real_length(nsg, beg, NULL)<=0)
    {
        return 0;
    }
    if(get_real_length(nsg, beg^1, NULL) == 1)
    {
        get_real_length(nsg, beg^1, &end);
        if(get_real_length(nsg, end^1, NULL) == 1)
        {
            return 0;
        }
    }
    ///if the contig here is too small
    primaryLen = get_long_tip_length(nsg, &(ug->u), beg, &end, NULL);
    if(primaryLen == 0 || primaryLen < minLongUntig || ug->u.a[beg>>1].circ)
    {
        return 0;
    }
    if(get_real_length(nsg, end, NULL) <= 0)
    {
        return 0;
    }
    /*****************************simple checking**********************************/
    (*r_beg) = beg; (*r_end) = end;
    
    /*****************************adjacent checking**********************************/
    uint32_t nw = asg_arc_n(nsg, end), n_arc = 0, next_uID, next_uID_verify;
    asg_arc_t *aw = asg_arc_a(nsg, end);
    for (i = 0; i < nw; i++)
    {
        if(!aw[i].del)
        {
            n_arc++;
            ////we don't want any long tip here
            if(if_long_tip_length(nsg, &(ug->u), aw[i].v, NULL, 
                minLongUntig, maxShortUntig, l_untig_rate, primaryLen)==1)
            {
                n_arc = 0;
                break;
            }
        }
    }
    if(n_arc == 0)
    {
        return 0;
    }
    /*****************************adjacent checking**********************************/

    ///here all out-nodes of v are small untigs
    if(explore_graph(nsg, beg, max_node_threshold, l_untig_rate, 
    minLongUntig, maxShortUntig, 1, b_0, &(u_vecs->a.a), &(u_vecs->a.m), 
    &(u_vecs->a.n), visit, &(ug->u), &next_uID) == 1)
    {
        (*r_next_uID) = next_uID;
        if(next_uID != (uint32_t)-1 && u_vecs->a.n == 1 
                        && IsMerge(ug->u, u_vecs->a.a[0]>>1) > 0)
        {
            return 0;
        }
        ///if next_uID == (uint32_t)-1, that means we found an end subgraph
        if(next_uID != (uint32_t)-1)
        {
            ///need to avoid missassembly
            ///merge all tangle as two types: 1. convex; 2, unroll them as a line
            ///should do somthing here, but we just skip it for covenience
            returnFlag = explore_graph(nsg, beg, max_node_threshold, l_untig_rate, 
            minLongUntig, maxShortUntig, 0, b_0, NULL, NULL, NULL, visit, &(ug->u), 
            &next_uID_verify);


            if((returnFlag == 1 && next_uID_verify != next_uID) || (returnFlag == 0))
            {
                //output_tangles(beg, next_uID, u_vecs->a.a, u_vecs->a.n, (char*)("###"));
                returnFlag = 0;
            }
            else
            {
                returnFlag = 1;
            }
            



            if(returnFlag == 1 && check_if_diploid_untigs(nsg, sg, beg, next_uID, &(ug->u),
                reverse_sources, minLongUntig-1, 0.3, b_0, b_1, 0, ruIndex) == 1)
            {
                ///output_tangles(beg, next_uID, u_vecs->a.a, u_vecs->a.n, (char*)("???"));
                returnFlag = 0;
            }

            if(returnFlag == 0)
            {
                merge_nodes(ug, end, next_uID, u_vecs->a.a, u_vecs->a.n, CONVEX);
                return CONVEX_M;
            }
        }

        ///actually it is not possible here
        if(u_vecs->a.n <= 0 || next_uID>>1 == beg>>1)
        {
            return 0;
        }
        ///output_tangles(beg, next_uID, u_vecs->a.a, u_vecs->a.n, (char*)(""));
        merge_nodes(ug, end, next_uID, u_vecs->a.a, u_vecs->a.n, UNROLL);
        if(next_uID != (uint32_t)-1)
        {
            (*r_next_uID) = next_uID;
            return UNROLL_M;
        }
        else ///end tangle
        {
            if(u_vecs->a.n > 0)
            {
                (*r_next_uID) = u_vecs->a.a[0];
            } 
            return UNROLL_E;
        }
    }

    return 0;
}


void adjust_asg_by_ug(ma_ug_t *ug, asg_t *read_g)
{
    uint32_t v, n_vtx, i = 0, qn;
    asg_t* nsg = ug->g;
    n_vtx = nsg->n_seq;
    ma_utg_t* node = NULL;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].c == ALTER_LABLE)
        {
            node = &(ug->u.a[v]);
            for (i = 0; i < node->n; i++)
            {
                qn = node->a[i]>>33;
                read_g->seq[qn].c = ALTER_LABLE;
            }
        }
    }


    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].c == ALTER_LABLE)
        {
            node = &(ug->u.a[v]);
            for (i = 0; i < node->n; i++)
            {
                qn = node->a[i]>>33;
                ///read_g->seq[qn].c = ALTER_LABLE;
                asg_seq_drop(read_g, qn);
            }
        }
    }

    asg_cleanup(read_g);
	asg_symm(read_g);
}

void lable_hap_asg_by_ug(ma_ug_t *ug, asg_t *read_g)
{
    uint32_t v, n_vtx, i = 0, qn;
    asg_t* nsg = ug->g;
    n_vtx = nsg->n_seq;
    ma_utg_t* node = NULL;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].c == HAP_LABLE)
        {
            node = &(ug->u.a[v]);
            for (i = 0; i < node->n; i++)
            {
                qn = node->a[i]>>33;
                read_g->seq[qn].c = HAP_LABLE;
            }
        }
    }

}

///qn is the read Id
void query_reverse_sources(asg_t *read_g, ma_hit_t_alloc* reverse_sources,
R_to_U* ruIndex, uint32_t qn, uint32_t self_offset, kvec_t_u64_warp* u_vecs)
{
    uint32_t i, rId, is_Unitig, cId;
    uint64_t mode;
    if(read_g->seq[qn].c == HAP_LABLE)
    {
        cId = (uint32_t)-1;
        mode = self_offset; mode = mode << 33; mode = mode|cId; mode = mode | (uint64_t)(0x100000000);
        kv_push(uint64_t, u_vecs->a, mode);
        return;
    }


    for (i = 0; i < reverse_sources[qn].length; i++)
    {
        rId = Get_tn(reverse_sources[qn].buffer[i]);
        ///there are three cases: 
        ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
        ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
        ///3. read has bee deleted, get_R_to_U() return the id of read that contains it 
        if(read_g->seq[rId].del == 1)
        {
            ///get the id of read that contains it 
            get_R_to_U(ruIndex, rId, &rId, &is_Unitig);
            if(rId == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[rId].del == 1) continue;
        }

        ///here rId is the id of the read coming from the different haplotype
        ///cId is the id of the corresponding contig (note here is the contig, instead of untig)
        get_R_to_U(ruIndex, rId, &cId, &is_Unitig);
        if(is_Unitig == 0) continue;
        ///if the read is at alternative contigs, cId might be (uint32_t)-1
        ///if(cId == (uint32_t)-1)
        mode = self_offset; mode = mode << 33; mode = mode|(uint64_t)(cId);
        kv_push(uint64_t, u_vecs->a, mode);
    }
    
}



uint32_t get_rId_from_contig_by_offset(ma_ug_t *ug, rIdContig* array, uint32_t offset)
{
    uint32_t uId, rId;
    ma_utg_t* reads;
    for (;array->untigI < array->b_0->b.n; array->untigI++)
    {
        uId = array->b_0->b.a[array->untigI]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        for (;array->readI < reads->n; array->readI++, array->offset++)
        {
            if(array->offset == offset)
            {
                rId = reads->a[array->readI]>>33;
                return rId;
            }
        }
        array->readI = 0;
    }


    array->offset = 0;
    for (array->untigI = 0;array->untigI < array->b_0->b.n; array->untigI++)
    {
        uId = array->b_0->b.a[array->untigI]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        for (array->readI = 0;array->readI < reads->n; array->readI++, array->offset++)
        {
            if(array->offset == offset)
            {
                rId = reads->a[array->readI]>>33;
                return rId;
            }
        }
    }

    array->untigI = array->readI = array->offset = 0;
    return (uint32_t)-1;
}

inline uint32_t get_contig_len(ma_ug_t *ug, buf_t* b_0)
{
    uint32_t uId, untigI, Len = 0;
    ma_utg_t* reads;
    
    for (untigI = 0;untigI < b_0->b.n; untigI++)
    {
        uId = b_0->b.a[untigI]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        Len = Len + reads->n;
    }

    return Len;
}

uint32_t get_offset_from_contig_by_rId(ma_ug_t *ug, rIdContig* array, uint32_t query_rId)
{
    uint32_t uId, rId;
    ma_utg_t* reads;
    for (;array->untigI < array->b_0->b.n; array->untigI++)
    {
        uId = array->b_0->b.a[array->untigI]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        for (;array->readI < reads->n; array->readI++, array->offset++)
        {
            rId = reads->a[array->readI]>>33;
            if(rId == query_rId) return array->offset;
        }
        array->readI = 0;
    }


    array->offset = 0;
    for (array->untigI = 0;array->untigI < array->b_0->b.n; array->untigI++)
    {
        uId = array->b_0->b.a[array->untigI]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        for (array->readI = 0;array->readI < reads->n; array->readI++, array->offset++)
        {
            rId = reads->a[array->readI]>>33;
            if(rId == query_rId) return array->offset;
        }
    }

    array->untigI = array->readI = array->offset = 0;
    return (uint32_t)-1;
}

///tn is the cId, tn_off is the order of rId with the same cId
uint32_t get_reverseId(asg_t *read_g, ma_hit_t_alloc* reverse_sources,
R_to_U* ruIndex, uint32_t qn, uint32_t tn, uint32_t tn_off)
{
    uint32_t i, rId, is_Unitig, cId, cId_off = 0;
    for (i = 0; i < reverse_sources[qn].length; i++)
    {
        rId = Get_tn(reverse_sources[qn].buffer[i]);
        ///there are three cases: 
        ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
        ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
        ///3. read has bee deleted, get_R_to_U() return the id of read that contains it 
        if(read_g->seq[rId].del == 1)
        {
            ///get the id of read that contains it 
            get_R_to_U(ruIndex, rId, &rId, &is_Unitig);
            if(rId == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[rId].del == 1) continue;
        }

        ///here rId is the id of the read coming from the different haplotype
        ///cId is the id of the corresponding contig (note here is the contig, instead of untig)
        get_R_to_U(ruIndex, rId, &cId, &is_Unitig);
        if(is_Unitig == 0) continue;
        ///if the read is at alternative contigs, cId might be (uint32_t)-1
        ///if(cId == (uint32_t)-1)
        if(cId == tn)
        {
            if(cId_off == tn_off) return rId;
            cId_off++;
        }
    }

    return (uint32_t)-1;
}

uint32_t get_contig_overlap_interval(uint32_t is_reverse, 
uint32_t q_beg, uint32_t q_end, uint32_t qLen, uint32_t t_beg, uint32_t t_end, uint32_t tLen, 
uint32_t* r_q_beg, uint32_t* r_q_end, uint32_t* r_t_beg, uint32_t* r_t_end)
{
    uint32_t k;
    (*r_q_beg) = (*r_q_end) = (*r_t_beg) = (*r_t_end) = (uint32_t)-1;
    if(q_end >= qLen || t_end >= tLen) return 0;
    if(q_beg >= qLen || t_beg >= tLen) return 0;

    if(is_reverse)
    {
        ///k = q_beg; q_beg = q_end; q_end = k;
        ///q_beg = qLen - q_beg - 1; q_end = qLen - q_end - 1; k = q_beg; q_beg = q_end; q_end = k;
        ///k = t_beg; t_beg = t_end; t_end = k;
        t_beg = tLen - t_beg - 1; t_end = tLen - t_end - 1; k = t_beg; t_beg = t_end; t_end = k;
    }




    if(q_beg <= t_beg)
    {
        t_beg = t_beg - q_beg; q_beg = 0;
    }
    else
    {
        q_beg = q_beg - t_beg; t_beg = 0;
    }

    uint32_t q_right_length = qLen - q_end - 1;
    uint32_t t_right_length = tLen - t_end - 1;
    if(q_right_length <= t_right_length)
    {
        q_end = qLen - 1; t_end = t_end + q_right_length;
    }
    else
    {
        t_end = tLen - 1; q_end = q_end + t_right_length;
    }

    if(is_reverse)
    {
        ///q_beg = qLen - q_beg - 1; q_end = qLen - q_end - 1; k = q_beg; q_beg = q_end; q_end = k;
        t_beg = tLen - t_beg - 1; t_end = tLen - t_end - 1; k = t_beg; t_beg = t_end; t_end = k;
    }

    (*r_q_beg) = q_beg; (*r_q_end) = q_end;
    (*r_t_beg) = t_beg; (*r_t_end) = t_end;
    return 1;
}


int get_haplotype_rate(ma_ug_t *ug, asg_t *read_g, ma_hit_t_alloc* reverse_sources,
R_to_U* ruIndex, buf_t* b_0, uint32_t beg, uint32_t end, uint32_t query_cId, float match_rate)
{
    uint32_t i, j, k, num, match_num, offset, uId, qn, rId, is_Unitig, cId;
    ma_utg_t* reads;

    offset = match_num = num = 0;
    for (i = 0; i < b_0->b.n; i++)
    {
        uId = b_0->b.a[i]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        for (j = 0; j < reads->n; j++, offset++)
        {
            qn = reads->a[j]>>33;
            for (k = 0; k < reverse_sources[qn].length; k++)
            {
                rId = Get_tn(reverse_sources[qn].buffer[k]);
                ///there are three cases: 
                ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
                ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
                ///3. read has bee deleted, get_R_to_U() return the id of read that contains it 
                if(read_g->seq[rId].del == 1)
                {
                    ///get the id of read that contains it 
                    get_R_to_U(ruIndex, rId, &rId, &is_Unitig);
                    if(rId == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[rId].del == 1) continue;
                }

                ///here rId is the id of the read coming from the different haplotype
                ///cId is the id of the corresponding contig (note here is the contig, instead of untig)
                get_R_to_U(ruIndex, rId, &cId, &is_Unitig);
                if(is_Unitig == 0) continue;
                ///if the read is at alternative contigs, cId might be (uint32_t)-1
                if(offset>=beg && offset <=end)
                {
                    num++;
                    if(cId == query_cId) match_num++;
                }
            }
        }
    } 

    if(match_num >= num*match_rate) return 1;
    return 0; 
}

#define fully_cover_rate 0.9
#define extraord_rate 0.2
#define hap_seed 20
#define GetCid(x) ((uint64_t)(0x7fffffff)&(uint64_t)(x))
#define GetOff(x) ((uint64_t)(x)>>33)
#define IfAlter(x) (GetCid((x))==(0x7fffffff))
#define IfColor(x) (((uint64_t)(0x100000000)&(uint64_t)(x))!=0)
#define IfVisit(x) (((uint64_t)(0x80000000)&(uint64_t)(x))!=0)
#define SetVisit(x) ((x) = ((uint64_t)(0x80000000)|(uint64_t)(x)))
inline int get_useful_contig_advance(ma_ug_t *ug, asg_t *read_g, ma_hit_t_alloc* reverse_sources,
R_to_U* ruIndex, kvec_t_u64_warp* r_vecs, buf_t* b_0, uint64_t* contigBeg, asg_t *bi_g, 
uint32_t currentId, float Hap_rate, uint32_t long_hap_overlap, float long_hap_overlap_rate)
{
    uint32_t i = 0, j, k, sLen, is_found, rId, beg, end, is_build;
    uint32_t q_beg, q_end, t_beg, t_end, qLen, tLen;
    uint32_t interval_q_beg, interval_q_end, interval_t_beg, interval_t_end;
    uint64_t anchor_offset, offset;
    uint64_t anchor_cId, cId;

    buf_t b_tid;
    memset(&b_tid, 0, sizeof(buf_t));

    kvec_t(Hap_Align) seed;
    kv_init(seed);

    kvec_t(Hap_Align) alignment;
    kv_init(alignment);

    Hap_Align x;
    rIdContig iterator, iterator_tid;
    iterator.b_0 = b_0;
    iterator.offset = iterator.readI = iterator.untigI = 0;
    if(r_vecs->a.n == 0) return -1;

    qLen = get_contig_len(ug, b_0);


    /**********************for debug****************************/
    ///if(r_vecs->a.n > 20 && r_vecs->a.n < 50)
    // {
    //     i = 0;
    //     fprintf(stderr,"*******\n");
    //     for (i = 0; i < r_vecs->a.n; i++)
    //     {
    //         offset = GetOff(r_vecs->a.a[i]); 
    //         cId = GetCid(r_vecs->a.a[i]);
    //         fprintf(stderr, "((%u) off: %u, cId: %d, HAP: %u)\n", i, offset, (int)cId, 
    //         IfColor(r_vecs->a.a[i]));
    //     }
    // }
    /**********************for debug****************************/




    i = 0;
    alignment.n = 0;
    while (i < r_vecs->a.n)
    {
        ///there are three cases:
        /**
         * 1) u_vecs->a.a[i] has haplotype infor at the primary contigs
         * 2) u_vecs->a.a[i] has haplotype infor at the alternative contigs, 
         *    IfAlter(u_vecs->a.a[i])==1
         * 3) u_vecs->a.a[i] itself is labled with HAP_LABLE, 
         *    IfColor(u_vecs->a.a[i])==1
         * 4) if u_vecs->a.a[i] has been labled
         *    IfVisit(u_vecs->a.a[i])==1
         * 4) u_vecs->a.a[i] does not have any haplotype infor. In this case, 
         *    self_offsetLen is not consecutive 
        **/

        if(IfAlter(r_vecs->a.a[i]) || IfColor(r_vecs->a.a[i]) || IfVisit(r_vecs->a.a[i]))
        {
            i++;
            continue;
        }
        

        anchor_offset = GetOff(r_vecs->a.a[i]); anchor_cId = GetCid(r_vecs->a.a[i]); 
        sLen = 0; is_found = 0; seed.n = 0;
        for (j = i; j < r_vecs->a.n; j++)
        {
            ///offset is the self offset at query
            ///cId is the target Id (contig)
            offset = GetOff(r_vecs->a.a[j]); 
            cId = GetCid(r_vecs->a.a[j]);

            ///means we go ahead one step at query
            if(offset != anchor_offset)
            {
                ///sLen must be >=1
                if(is_found == 0) 
                {
                    break;
                }
                is_found = 0;
                sLen++;
                anchor_offset = offset;
            } 

            if(IfVisit(r_vecs->a.a[j]) || IfColor(r_vecs->a.a[j])) continue;

            if(cId == anchor_cId)
            {
                x.q_pos = GetOff(r_vecs->a.a[j]);
                x.t_id = GetCid(r_vecs->a.a[j]);
                x.is_color = 1;
                x.t_pos = (uint32_t)-1;
                SetVisit(r_vecs->a.a[j]);
                is_found++;
                kv_push(Hap_Align, seed, x);
            } 
        }

        if(j == r_vecs->a.n && is_found > 0)
        {
            sLen++;
        }


        // for (j = 1; j < seed.n; j++)
        // {
        //     if(seed.a[j].t_id != seed.a[j-1].t_id)
        //     {
        //         fprintf(stderr, "hehe\n");
        //         fprintf(stderr, "sLen: %u, seed.n: %u\n", sLen, seed.n);
        //     }

        //     if(seed.a[j].q_pos < seed.a[j-1].q_pos)
        //     {
        //         fprintf(stderr, "haha\n");
        //         fprintf(stderr, "sLen: %u, seed.n: %u\n", sLen, seed.n);
        //     }
        // }

        ///don't need to know how long of this seed
        if(sLen >= hap_seed && seed.n > 0)
        {
            /**********************for debug****************************/
            ///if(r_vecs->a.n > 20 && r_vecs->a.n < 50)
            ///fprintf(stderr, "\nsLen: %u, seed.n: %u, anchor_cId: %u\n", sLen, seed.n, anchor_cId);
            /**********************for debug****************************/


            uint32_t inner_off = 0, pre_q_pos = (uint32_t)-1, RrId;
            iterator.offset = iterator.readI = iterator.untigI = 0;

            ///for one vector, all t_id should be same; that is why we recover tid here
            beg = (uint32_t)(contigBeg[seed.a[0].t_id]);
            b_tid.b.n = 0;ug->u.a[beg>>1].circ = 0;
            get_long_tip_length(ug->g, &(ug->u), beg, &end, &b_tid);
            iterator_tid.b_0 = &b_tid;
            iterator_tid.offset = iterator_tid.readI = iterator_tid.untigI = 0;
            ///we have already got qLen at the begining
            tLen = get_contig_len(ug, &b_tid);

            for (k = 0; k < seed.n; k++)
            {

                rId = get_rId_from_contig_by_offset(ug, &iterator, seed.a[k].q_pos);
                ///actually we shouldn't have this case
                if(rId == (uint32_t)-1) continue;
                if(pre_q_pos != seed.a[k].q_pos)
                {
                    inner_off = 0;
                    pre_q_pos = seed.a[k].q_pos;
                }
                else
                {
                    inner_off++;
                }
                
                RrId = get_reverseId(read_g, reverse_sources, ruIndex, rId, seed.a[k].t_id, 
                inner_off);
                ///if(RrId == (uint32_t)-1) fprintf(stderr, "ERROR1\n");
                ///for one vector, all t_id should be same
                seed.a[k].t_pos = get_offset_from_contig_by_rId(ug, &iterator_tid, RrId);
                // if(seed.a[k].t_pos == (uint32_t)-1) fprintf(stderr, "ERROR2\n");
                // uint32_t debug_uID, debug_is_Unitig;
                // get_R_to_U(ruIndex, RrId, &debug_uID, &debug_is_Unitig);
                // if(seed.a[k].t_id != debug_uID) fprintf(stderr, "ERROR1\n");
                

                /**********************for debug****************************/
                ///if(r_vecs->a.n > 20 && r_vecs->a.n < 50)
                // fprintf(stderr, "%u, q_pos: %u, t_id: %u, t_pos: %u, inner_off: %u, rId: %u, RrId: %u\n", 
                // k, seed.a[k].q_pos, seed.a[k].t_id, seed.a[k].t_pos, inner_off, rId, RrId);
                /**********************for debug****************************/
            }

            q_beg = t_beg = (uint32_t)-1; q_end = t_end = 0;
            for (k = 0; k < seed.n; k++)
            {
                if(seed.a[k].t_pos < t_beg) t_beg = seed.a[k].t_pos;
                if(seed.a[k].t_pos > t_end) t_end = seed.a[k].t_pos;

                if(seed.a[k].q_pos < q_beg) q_beg = seed.a[k].q_pos;
                if(seed.a[k].q_pos > q_end) q_end = seed.a[k].q_pos;
            }

            

            if(q_beg<=q_end && t_beg<=t_end && tLen > 0 && qLen > 0)
            {
                ///exclude extraordinary points
                if((DIFF((q_end+1-q_beg), (t_end+1-t_beg))) > (q_end+1-q_beg)*extraord_rate)
                {
                    /**********************for debug****************************/
                    // fprintf(stderr, "\nsLen: %u, seed.n: %u\n", sLen, seed.n);
                    // for (k = 0; k < seed.n; k++)
                    // {
                    //     fprintf(stderr, "%u, q_pos: %u, t_id: %u, t_pos: %u\n", 
                    //     k, seed.a[k].q_pos, seed.a[k].t_id, seed.a[k].t_pos);
                    // }
                    // fprintf(stderr, "-q_beg: %u, q_end: %u, qLen: %u, t_beg: %u, t_end: %u, tLen: %u, anchor_cId: %u\n", 
                    // q_beg, q_end, qLen, t_beg, t_end, tLen, anchor_cId);
                    /**********************for debug****************************/

                    uint32_t leftLen = 0, rightLen = 0, m = 0, median = 0, diff = (q_end+1-q_beg)*(1+extraord_rate)*0.5;
                    for (k = 0; k < seed.n; k++)
                    {
                        median += seed.a[k].t_pos;
                    }
                    median = median / seed.n;

                    leftLen = median - diff;
                    if(median < diff) leftLen = 0;

                    rightLen = median + diff;
                    if(rightLen >= tLen) rightLen = tLen - 1;

                    t_beg = (uint32_t)-1; t_end = 0; m = 0;
                    for (k = 0; k < seed.n; k++)
                    {
                        if(seed.a[k].t_pos >= leftLen && seed.a[k].t_pos <= rightLen)
                        {
                            if(seed.a[k].t_pos < t_beg) t_beg = seed.a[k].t_pos;
                            if(seed.a[k].t_pos > t_end) t_end = seed.a[k].t_pos;
                            seed.a[m] = seed.a[k];
                            m++;
                        }
                    }
                    seed.n = m;

                    if(t_beg>t_end) goto direct_skip;
                    /**********************for debug****************************/
                    // fprintf(stderr, "+q_beg: %u, q_end: %u, qLen: %u, t_beg: %u, t_end: %u, tLen: %u\n", 
                    // q_beg, q_end, qLen, t_beg, t_end, tLen);
                    /**********************for debug****************************/
                }
                ///here we know q_beg, q_end, qLen
                ///and          t_beg, t_end, tLen
                ///there might be two directions: 
                ///a) query and target at the same direction
                ///b) query and target at different direction
                get_contig_overlap_interval(0, q_beg, q_end, qLen, 
                t_beg, t_end, tLen, &interval_q_beg, &interval_q_end,
                &interval_t_beg, &interval_t_end);

                uint32_t flag_forward = get_haplotype_rate(ug, read_g, reverse_sources, ruIndex, b_0, 
                interval_q_beg, interval_q_end, anchor_cId, Hap_rate);
                if(flag_forward)
                {
                    x.t_id = anchor_cId;
                    x.q_pos = interval_q_beg;
                    x.t_pos = interval_q_end;
                    ///x.is_color = 0;
                    x.is_color = tLen;
                    kv_push(Hap_Align, alignment, x);
                }
                

                /**********************for debug****************************/
                ////if(r_vecs->a.n > 20 && r_vecs->a.n < 50)
                // if(flag_forward == 1)
                // {
                //     fprintf(stderr, "\nanchor_cId: %u\n", anchor_cId);
                //     fprintf(stderr, "(0) q_beg: %u, q_end: %u, qLen: %u, t_beg: %u, t_end: %u, tLen: %u\n", 
                //     q_beg, q_end, qLen, t_beg, t_end, tLen);

                //     fprintf(stderr, "(0) interval_q_beg: %u, interval_q_end: %u, interval_t_beg: %u, interval_t_end: %u, flag_forward: %u\n", 
                //     interval_q_beg, interval_q_end, interval_t_beg, interval_t_end, flag_forward);   
                // }
                /**********************for debug****************************/


                get_contig_overlap_interval(1, q_beg, q_end, qLen, 
                t_beg, t_end, tLen, &interval_q_beg, &interval_q_end,
                &interval_t_beg, &interval_t_end);

                uint32_t flag_backward = get_haplotype_rate(ug, read_g, reverse_sources, ruIndex, b_0, 
                interval_q_beg, interval_q_end, anchor_cId, Hap_rate);
                if(flag_backward)
                {
                    x.t_id = anchor_cId;
                    x.q_pos = interval_q_beg;
                    x.t_pos = interval_q_end;
                    ///x.is_color = 1;
                    x.is_color = tLen;
                    kv_push(Hap_Align, alignment, x);
                }

                /**********************for debug****************************/
                ////if(r_vecs->a.n > 20 && r_vecs->a.n < 50)
                // if(flag_backward == 1)
                // {
                //     fprintf(stderr, "\nanchor_cId: %u\n", anchor_cId);
                //     fprintf(stderr, "(1) q_beg: %u, q_end: %u, qLen: %u, t_beg: %u, t_end: %u, tLen: %u\n", 
                //     q_beg, q_end, qLen, t_beg, t_end, tLen);

                //     fprintf(stderr, "(1) interval_q_beg: %u, interval_q_end: %u, interval_t_beg: %u, interval_t_end: %u, flag_backward: %u\n", 
                //     interval_q_beg, interval_q_end, interval_t_beg, interval_t_end, flag_backward);   
                // }
                /**********************for debug****************************/
                
            
            }
            
        }

        direct_skip:
        i++;
    }

    ///must sort here
    radix_sort_Hap_Align_sort(alignment.a, alignment.a+alignment.n);
    /**********************for debug****************************/
    // fprintf(stderr, "\nalignment.n: %u\n", (uint32_t)alignment.n);
    // for (i = 0; i < alignment.n; i++)
    // {
    //     fprintf(stderr, "c_id: %u, beg: %u, end: %u, tLen: %u\n", 
    //      alignment.a[i].t_id, alignment.a[i].q_pos, 
    //     alignment.a[i].t_pos, alignment.a[i].is_color);
    // }
    /**********************for debug****************************/

    uint32_t m = 0;
    sLen = 0; anchor_cId = tLen = q_beg = q_end = (uint32_t)-1;
    for (i = 0; i < alignment.n; i++)
    {
        if(anchor_cId != alignment.a[i].t_id)
        {
            if(anchor_cId != (uint32_t)-1)
            {
                alignment.a[m].t_id = anchor_cId;
                alignment.a[m].q_pos = q_beg;
                alignment.a[m].t_pos = q_end;
                alignment.a[m].is_color = tLen;
                m++;
            }
            anchor_cId = alignment.a[i].t_id;
            sLen = 0;
        }
        if(alignment.a[i].t_pos < alignment.a[i].q_pos) continue;
        if(alignment.a[i].t_pos - alignment.a[i].q_pos + 1 > sLen)
        {
            sLen = alignment.a[i].t_pos - alignment.a[i].q_pos + 1;
            q_beg = alignment.a[i].q_pos;
            q_end = alignment.a[i].t_pos;
            tLen = alignment.a[i].is_color;
        }
    }
    if(anchor_cId != (uint32_t)-1)
    {
        alignment.a[m].t_id = anchor_cId;
        alignment.a[m].q_pos = q_beg;
        alignment.a[m].t_pos = q_end;
        alignment.a[m].is_color = tLen;
        m++;
    }
    alignment.n = m;
    /**********************for debug****************************/
    // fprintf(stderr, "***alignment.m: %u\n", (uint32_t)alignment.n);
    // for (i = 0; i < alignment.n; i++)
    // {
    //     fprintf(stderr, "c_id: %u, beg: %u, end: %u, tLen: %u\n", 
    //      alignment.a[i].t_id, alignment.a[i].q_pos, 
    //     alignment.a[i].t_pos, alignment.a[i].is_color);
    // }
    /**********************for debug****************************/
    asg_arc_t *e;
    for (i = 0; i < alignment.n; i++)
    {
        is_build = 0;
        if(alignment.a[i].t_pos < alignment.a[i].q_pos) continue;
        sLen = alignment.a[i].t_pos - alignment.a[i].q_pos + 1;
        tLen = alignment.a[i].is_color;
        ///if overlap is short, must fully cover one of the read
        if(sLen < long_hap_overlap)
        {
            if(sLen >= (fully_cover_rate*tLen) || sLen >= (fully_cover_rate*qLen))
            {
                is_build = 1;
            }
        }
        else///if is long, can be partly cover one of the read
        {
            if(sLen >= (long_hap_overlap_rate*tLen) || sLen >= (long_hap_overlap_rate*qLen))
            {
                is_build = 1;
            }
        }
        if(is_build)
        {
            e = asg_arc_pushp(bi_g);
            e->del = 0;
            e->ol = sLen;
            e->ul = currentId; e->ul = e->ul << 32; e->ul = e->ul | (uint64_t)(qLen);
            e->v = alignment.a[i].t_id;
        }
    }
    /**
    is_build = 0;
    ///uint32_t long_hap_overlap, float long_hap_overlap_rate
    ///if overlap is short, must fully cover one of the read
    if(sLen < long_hap_overlap)
    {
        if(sLen >= (fully_cover_rate*tLen) || sLen >= (fully_cover_rate*qLen))
        {
            is_build = 1;
        }
    }
    else///if is long, can be partly cover one of the read
    {
        if(sLen >= (long_hap_overlap_rate*tLen) || sLen >= (long_hap_overlap_rate*qLen))
        {
            is_build = 1;
        }
    }

    if(is_build)
    {
        asg_arc_t *e;
        e = asg_arc_pushp(bi_g);
        e->del = 0;
        e->ol = flag;
        e->ul = i; e->ul = e->ul << 32; e->ul = e->ul | (uint64_t)(u_vecs.a.n);
        e->v = cId;
    }
    **/
    

    /**
    radix_sort_Hap_Align_sort(alignment.a, alignment.a+alignment.n);
    fprintf(stderr, "alignment.n: %u, sLen: %u, q_beg: %u, q_end: %u\n", 
    (uint32_t)alignment.n, sLen, q_beg, q_end);
    for (i = 0; i < alignment.n; i++)
    {
        fprintf(stderr, "dir: %u, c_id: %u, beg: %u, end: %u\n", 
        alignment.a[i].is_color, alignment.a[i].t_id, alignment.a[i].q_pos, 
        alignment.a[i].t_pos);
    }
    **/

    


    /**
    uint32_t self_offset, uId;
    ma_utg_t* reads;
    iterator.offset = iterator.readI = iterator.untigI = 0;
    for (j = 0, self_offset = 0; j < b_0->b.n; j++)
    {
        uId = b_0->b.a[j]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        ///scan all reads
        ///self_offset will skip fake(merge) nodes, but not skip HAP_LABLE
        for (k = 0; k < reads->n; k++, self_offset++)
        {
            rId = reads->a[k]>>33;
            if(get_rId_from_contig_by_offset(ug, &iterator, self_offset)!=rId)
            {
                fprintf(stderr, "ERROR1\n");
            }
        }
    }

    iterator.offset = iterator.readI = iterator.untigI = 0;
    for (long long debug_i = self_offset; debug_i >= 0; debug_i--)
    {
        for (j = 0, self_offset = 0; j < b_0->b.n; j++)
        {
            uId = b_0->b.a[j]>>1;
            ///if(IsMerge(ug->u, uId)>0) continue;
            reads = &(ug->u.a[uId]);
            ///scan all reads
            ///self_offset will skip fake(merge) nodes, but not skip HAP_LABLE
            for (k = 0; k < reads->n; k++, self_offset++)
            {
                if(debug_i != self_offset) continue;
                rId = reads->a[k]>>33;
                if(get_rId_from_contig_by_offset(ug, &iterator, self_offset)!=rId)
                {
                    fprintf(stderr, "ERROR2: self_offset: %u\n", self_offset);
                }
            }
        }
    }
    


    


    iterator.offset = iterator.readI = iterator.untigI = 0;
    iterator.offset = iterator.readI = iterator.untigI = 0;
    for (j = 0, self_offset = 0; j < b_0->b.n; j++)
    {
        uId = b_0->b.a[j]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        ///scan all reads
        ///self_offset will skip fake(merge) nodes, but not skip HAP_LABLE
        for (k = 0; k < reads->n; k++, self_offset++)
        {
            rId = reads->a[k]>>33;
            if(get_offset_from_contig_by_rId(ug, &iterator, rId)!=self_offset)
            {
                fprintf(stderr, "ERROR3\n");
            }
        }
    }

    iterator.offset = iterator.readI = iterator.untigI = 0;
    for (long long debug_i = self_offset; debug_i >= 0; debug_i--)
    {
        for (j = 0, self_offset = 0; j < b_0->b.n; j++)
        {
            uId = b_0->b.a[j]>>1;
            ///if(IsMerge(ug->u, uId)>0) continue;
            reads = &(ug->u.a[uId]);
            ///scan all reads
            ///self_offset will skip fake(merge) nodes, but not skip HAP_LABLE
            for (k = 0; k < reads->n; k++, self_offset++)
            {
                if(debug_i != self_offset) continue;
                rId = reads->a[k]>>33;
                if(get_offset_from_contig_by_rId(ug, &iterator, rId)!=self_offset)
                {
                    fprintf(stderr, "ERROR33: self_offset: %u\n", self_offset);
                }
            }
        }
    }
    **/


    kv_destroy(alignment);
    kv_destroy(seed);
    free(b_tid.b.a);
    return 1;
}

inline int get_useful_contig_advance_back(ma_ug_t *ug, asg_t *read_g, ma_hit_t_alloc* reverse_sources,
R_to_U* ruIndex, kvec_t_u64_warp* u_vecs, buf_t* b_0, uint64_t* contigBeg)
{
    uint32_t i = 0, j, k, sLen, is_found, rId, beg, end;
    uint32_t q_beg, q_end, t_beg, t_end, qLen, tLen;
    uint32_t interval_q_beg, interval_q_end, interval_t_beg, interval_t_end;
    uint64_t anchor_offset, offset;
    uint64_t anchor_cId, cId;

    buf_t b_tid;
    memset(&b_tid, 0, sizeof(buf_t));

    kvec_t(Hap_Align) seed;
    kv_init(seed);
    Hap_Align x;
    rIdContig iterator, iterator_tid;
    iterator.b_0 = b_0;
    iterator.offset = iterator.readI = iterator.untigI = 0;
    if(u_vecs->a.n == 0) return -1;

    qLen = get_contig_len(ug, b_0);
    
    /**********************for debug****************************/
    // if(u_vecs->a.n > 20 && u_vecs->a.n < 50)
    // {
    //     i = 0;
    //     fprintf(stderr,"*******\n");
    //     for (i = 0; i < u_vecs->a.n; i++)
    //     {
    //         offset = u_vecs->a.a[i]>>33; 
    //         cId = (uint32_t)(u_vecs->a.a[i]);
    //         fprintf(stderr, "((%u) off: %u, cId: %d, HAP: %u)\n", i, offset, (int)cId, 
    //         ((u_vecs->a.a[i]&(uint64_t)(0x100000000))!=0));
    //     }
    // }
    /**********************for debug****************************/


    i = 0;
    while (i < u_vecs->a.n)
    {
        ///there are three cases:
        /**
         * 1) u_vecs->a.a[i] has haplotype infor at the primary contigs
         * 2) u_vecs->a.a[i] has haplotype infor at the alternative contigs, 
         *    (uint32_t)u_vecs->a.a[i]==(uint32_t)-1
         * 3) u_vecs->a.a[i] itself is labled with HAP_LABLE, 
         *    u_vecs->a.a[i]&(uint64_t)(0x100000000)>0 && (uint32_t)u_vecs->a.a[i]==(uint32_t)-1
         * 4) u_vecs->a.a[i] does not have any haplotype infor. In this case, 
         *    self_offsetLen is not consecutive 
         * 
        **/
        if((u_vecs->a.a[i]&(uint64_t)(0x100000000)) || ((uint32_t)u_vecs->a.a[i]==((uint32_t)-1)))
        {
            i++;
            continue;
        }
        

        anchor_offset = (u_vecs->a.a[i]>>33); anchor_cId = (uint32_t)(u_vecs->a.a[i]); 
        sLen = 0; is_found = 0; seed.n = 0;
        for (j = i; j < u_vecs->a.n; j++)
        {
            ///offset is the self offset at query
            ///cId is the target Id (contig)
            offset = u_vecs->a.a[j]>>33; 
            cId = (uint32_t)(u_vecs->a.a[j]);

            ///means we go ahead one step at query
            if(offset != anchor_offset)
            {
                ///sLen must be >=1
                if(is_found == 0) 
                {
                    break;
                }
                is_found = 0;
                sLen++;
                anchor_offset = offset;
            } 

            if(cId == anchor_cId)
            {
                x.q_pos = u_vecs->a.a[j]>>33;
                x.t_id = (uint32_t)(u_vecs->a.a[j]);
                x.is_color = 1;
                x.t_pos = (uint32_t)-1;
                u_vecs->a.a[j] = u_vecs->a.a[j]|0xffffffff;
                is_found++;
                // if(seed.n > 0 && seed.a[seed.n-1].t_id == x.t_id && 
                // seed.a[seed.n-1].q_pos == x.q_pos)
                // {
                //     seed.a[seed.n-1].is_color++;
                //     continue;
                // }
                kv_push(Hap_Align, seed, x);
            } 
        }

        if(j == u_vecs->a.n && is_found > 0)
        {
            sLen++;
        }


        // for (j = 1; j < seed.n; j++)
        // {
        //     if(seed.a[j].t_id != seed.a[j-1].t_id)
        //     {
        //         fprintf(stderr, "hehe\n");
        //         fprintf(stderr, "sLen: %u, seed.n: %u\n", sLen, seed.n);
        //     }

        //     if(seed.a[j].q_pos < seed.a[j-1].q_pos)
        //     {
        //         fprintf(stderr, "haha\n");
        //         fprintf(stderr, "sLen: %u, seed.n: %u\n", sLen, seed.n);
        //     }
        // }

        ///don't need to know how long of this seed
        if(sLen >= hap_seed && seed.n > 0)
        {
            /**********************for debug****************************/
            ///if(u_vecs->a.n > 20 && u_vecs->a.n < 50)
            ///fprintf(stderr, "\nsLen: %u, seed.n: %u, anchor_cId: %u\n", sLen, seed.n, anchor_cId);
            /**********************for debug****************************/
            
            uint32_t inner_off = 0, pre_q_pos = (uint32_t)-1, RrId;
            iterator.offset = iterator.readI = iterator.untigI = 0;

            ///for one vector, all t_id should be same; that is why we recover tid here
            beg = (uint32_t)(contigBeg[seed.a[0].t_id]);
            b_tid.b.n = 0;ug->u.a[beg>>1].circ = 0;
            get_long_tip_length(ug->g, &(ug->u), beg, &end, &b_tid);
            iterator_tid.b_0 = &b_tid;
            iterator_tid.offset = iterator_tid.readI = iterator_tid.untigI = 0;
            ///we have already got qLen at the begining
            tLen = get_contig_len(ug, &b_tid);

            for (k = 0; k < seed.n; k++)
            {

                rId = get_rId_from_contig_by_offset(ug, &iterator, seed.a[k].q_pos);
                ///actually we shouldn't have this case
                if(rId == (uint32_t)-1) continue;
                if(pre_q_pos != seed.a[k].q_pos)
                {
                    inner_off = 0;
                    pre_q_pos = seed.a[k].q_pos;
                }
                else
                {
                    inner_off++;
                }
                
                RrId = get_reverseId(read_g, reverse_sources, ruIndex, rId, seed.a[k].t_id, 
                inner_off);
                ///if(RrId == (uint32_t)-1) fprintf(stderr, "ERROR1\n");
                ///for one vector, all t_id should be same
                seed.a[k].t_pos = get_offset_from_contig_by_rId(ug, &iterator_tid, RrId);
                // if(seed.a[k].t_pos == (uint32_t)-1) fprintf(stderr, "ERROR2\n");
                // uint32_t debug_uID, debug_is_Unitig;
                // get_R_to_U(ruIndex, RrId, &debug_uID, &debug_is_Unitig);
                // if(seed.a[k].t_id != debug_uID) fprintf(stderr, "ERROR1\n");
                
                /**********************for debug****************************/
                ///if(u_vecs->a.n > 20 && u_vecs->a.n < 50)
                // fprintf(stderr, "%u, q_pos: %u, t_id: %u, t_pos: %u, inner_off: %u, rId: %u, RrId: %u\n", 
                // k, seed.a[k].q_pos, seed.a[k].t_id, seed.a[k].t_pos, inner_off, rId, RrId);
                /**********************for debug****************************/
            }

            q_beg = t_beg = (uint32_t)-1; q_end = t_end = 0;
            for (k = 0; k < seed.n; k++)
            {
                if(seed.a[k].t_pos < t_beg) t_beg = seed.a[k].t_pos;
                if(seed.a[k].t_pos > t_end) t_end = seed.a[k].t_pos;

                if(seed.a[k].q_pos < q_beg) q_beg = seed.a[k].q_pos;
                if(seed.a[k].q_pos > q_end) q_end = seed.a[k].q_pos;
            }

            

            if(q_beg<=q_end && t_beg<=t_end)
            {
                ///here we know q_beg, q_end, qLen
                ///and          t_beg, t_end, tLen
                ///there might be two directions: 
                ///a) query and target at the same direction
                ///b) query and target at different direction
                get_contig_overlap_interval(0, q_beg, q_end, qLen, 
                t_beg, t_end, tLen, &interval_q_beg, &interval_q_end,
                &interval_t_beg, &interval_t_end);


                /**********************for debug****************************/
                ///if(u_vecs->a.n > 20 && u_vecs->a.n < 50)
                {
                    fprintf(stderr, "q_beg: %u, q_end: %u, qLen: %u, t_beg: %u, t_end: %u, tLen: %u\n", 
                    q_beg, q_end, qLen, t_beg, t_end, tLen);

                    fprintf(stderr, "interval_q_beg: %u, interval_q_end: %u, interval_t_beg: %u, interval_t_end: %u\n", 
                    interval_q_beg, interval_q_end, interval_t_beg, interval_t_end);   
                }
                /**********************for debug****************************/
            
            }
            
        }
        i++;
    }






    /**
    uint32_t self_offset, uId;
    ma_utg_t* reads;
    iterator.offset = iterator.readI = iterator.untigI = 0;
    for (j = 0, self_offset = 0; j < b_0->b.n; j++)
    {
        uId = b_0->b.a[j]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        ///scan all reads
        ///self_offset will skip fake(merge) nodes, but not skip HAP_LABLE
        for (k = 0; k < reads->n; k++, self_offset++)
        {
            rId = reads->a[k]>>33;
            if(get_rId_from_contig_by_offset(ug, &iterator, self_offset)!=rId)
            {
                fprintf(stderr, "ERROR1\n");
            }
        }
    }

    iterator.offset = iterator.readI = iterator.untigI = 0;
    for (long long debug_i = self_offset; debug_i >= 0; debug_i--)
    {
        for (j = 0, self_offset = 0; j < b_0->b.n; j++)
        {
            uId = b_0->b.a[j]>>1;
            ///if(IsMerge(ug->u, uId)>0) continue;
            reads = &(ug->u.a[uId]);
            ///scan all reads
            ///self_offset will skip fake(merge) nodes, but not skip HAP_LABLE
            for (k = 0; k < reads->n; k++, self_offset++)
            {
                if(debug_i != self_offset) continue;
                rId = reads->a[k]>>33;
                if(get_rId_from_contig_by_offset(ug, &iterator, self_offset)!=rId)
                {
                    fprintf(stderr, "ERROR2: self_offset: %u\n", self_offset);
                }
            }
        }
    }
    


    


    iterator.offset = iterator.readI = iterator.untigI = 0;
    iterator.offset = iterator.readI = iterator.untigI = 0;
    for (j = 0, self_offset = 0; j < b_0->b.n; j++)
    {
        uId = b_0->b.a[j]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        ///scan all reads
        ///self_offset will skip fake(merge) nodes, but not skip HAP_LABLE
        for (k = 0; k < reads->n; k++, self_offset++)
        {
            rId = reads->a[k]>>33;
            if(get_offset_from_contig_by_rId(ug, &iterator, rId)!=self_offset)
            {
                fprintf(stderr, "ERROR3\n");
            }
        }
    }

    iterator.offset = iterator.readI = iterator.untigI = 0;
    for (long long debug_i = self_offset; debug_i >= 0; debug_i--)
    {
        for (j = 0, self_offset = 0; j < b_0->b.n; j++)
        {
            uId = b_0->b.a[j]>>1;
            ///if(IsMerge(ug->u, uId)>0) continue;
            reads = &(ug->u.a[uId]);
            ///scan all reads
            ///self_offset will skip fake(merge) nodes, but not skip HAP_LABLE
            for (k = 0; k < reads->n; k++, self_offset++)
            {
                if(debug_i != self_offset) continue;
                rId = reads->a[k]>>33;
                if(get_offset_from_contig_by_rId(ug, &iterator, rId)!=self_offset)
                {
                    fprintf(stderr, "ERROR33: self_offset: %u\n", self_offset);
                }
            }
        }
    }
    **/



    kv_destroy(seed);
    free(b_tid.b.a);
    return 1;
}


#define contig_seed 20
inline int get_useful_contig(kvec_t_u64_warp* u_vecs, float density, uint32_t miniLen, uint32_t* r_cId)
{
    (*r_cId) = (uint32_t)-1;
    uint32_t cId;
    uint32_t intervalLen = 0, realLen = 0, useful_index = (uint32_t)-1;
    uint32_t i = u_vecs->i, end;
    float T_density = density/2;
    
    if(i >= u_vecs->a.n) return -1;
    if((uint32_t)(u_vecs->a.a[i]) == (uint32_t)-1)
    {
        u_vecs->i++;
        return 0;
    } 

    end = i + contig_seed;
    if(end > u_vecs->a.n) end = u_vecs->a.n;
    cId = (uint32_t)(u_vecs->a.a[i]);
    for (; i < end; i++)
    {
        if(cId == (uint32_t)(u_vecs->a.a[i])) realLen++;
        intervalLen++;
        if(realLen >= intervalLen*density) useful_index = i;
    }

    if(realLen < intervalLen*T_density)
    {
        u_vecs->i++;
        return 0;
    } 
    ///here we found a useful seed
    ///it seems we don't need to scan backward
    for (; i < u_vecs->a.n; i++)
    {
        if(cId == (uint32_t)(u_vecs->a.a[i])) realLen++;
        intervalLen++;
        if(realLen < intervalLen*T_density) break;
        if(realLen >= intervalLen*density) useful_index = i;
    }

    if(useful_index == (uint32_t)-1 || (useful_index - u_vecs->i) < miniLen)
    {
        u_vecs->i++;
        return 0;
    }

    ///don't need to set backward
    end = (uint32_t)-1;
    for (i = u_vecs->i; i < u_vecs->a.n; i++)
    {
        if(cId == (uint32_t)(u_vecs->a.a[i]))
        {
            ///u_vecs->a.a[i] = (uint32_t)-1;
            u_vecs->a.a[i] = u_vecs->a.a[i]|(uint64_t)(0xffffffff);
        } 
    }
    (*r_cId) = cId;
    u_vecs->i++;
    return (useful_index - u_vecs->i);
}





#define UNVISIT (uint32_t)(0x7fffffff)
#define RED 0
#define BLACK 1
#define ISO 2
#define LABLE 3
void bi_paration(asg_t *bi_g, uint64_t* array, uint32_t bi_graph_Len)
{
    
    kvec_t(uint64_t) a;
    kv_init(a);
    a.a = array;
    a.n = a.m = bi_g->n_seq;
    kdq_t(uint32_t) *buf;
    buf = kdq_init(uint32_t);
    uint32_t i, len, v, w, k, nv;
    asg_arc_t *av;
    for (i = 0; i < bi_g->n_seq; i++)
    {
        bi_g->seq[i].c = 0;
    }

    for (i = 0; i < bi_g->n_seq; i++)
    {
        /****************************may have bugs********************************/
        ///how many reads contained in this contig
        len = (a.a[i]>>33);
        /****************************may have bugs********************************/
        ///if the contig is too small, or the contig has already been visited
        ///if(len < bi_graph_Len || bi_g->seq[i].len != UNVISIT) continue;
        if(len < bi_graph_Len || bi_g->seq[i].c == 1) continue;

        ///set the color of this node
        bi_g->seq[i].len = RED; bi_g->seq[i].c = 1;
        kdq_push(uint32_t, buf, i);
        while (kdq_size(buf) != 0)
        {
            v = *(kdq_pop(uint32_t, buf)); bi_g->seq[v].c = 1;
            if(bi_g->seq[v].len == ISO) continue;
            nv = asg_arc_n(bi_g, v);
            av = asg_arc_a(bi_g, v);


            ///get all out-nodes of v
            for(k = 0; k < nv; k++)
            {
                w = av[k].v;
                /****************************may have bugs********************************/
                ///if(bi_g->seq[w].len == bi_g->seq[v].len)
                len = (a.a[w]>>33);
                ///only check large contig
                if(len >= bi_graph_Len && bi_g->seq[w].len == bi_g->seq[v].len)
                { /****************************may have bugs********************************/
                    break;
                }
            }
            ///means v is conflict 
            if(k != nv)
            {
                bi_g->seq[v].len = ISO; bi_g->seq[v].c = 1;
                continue;
            }



            ///if v is not conflict with others
            for(k = 0; k < nv; k++)
            {
                w = av[k].v;
                ///if the out-node has not been visited
                if(bi_g->seq[w].len == UNVISIT)
                {
                    bi_g->seq[w].len = 1 - bi_g->seq[v].len;
                    bi_g->seq[w].c = 1;
                    /****************************may have bugs********************************/
                    len = (a.a[w]>>33);
                    /****************************may have bugs********************************/
                    if(len < bi_graph_Len) continue;
                    kdq_push(uint32_t, buf, w);
                }
                ///here bi_g->seq[w].len might be ISO or another color
                ///don't need to do anything here
            }
        }
    }


    //secondary checking
    // uint32_t a_color;
    // kdq_size(buf) = 0;
    // for (i = 0; i < bi_g->n_seq; i++)
    // {
    //     /****************************may have bugs********************************/
    //     len = (a.a[i]>>33);
    //     /****************************may have bugs********************************/
    //     ///just check end contig
    //     if(bi_g->seq[i].len != UNVISIT || (a.a[i] & (uint64_t)(0x100000000)) == 0) continue;
        
    //     v = i;
    //     nv = asg_arc_n(bi_g, v);
    //     av = asg_arc_a(bi_g, v);
    //     a_color = UNVISIT;

    //     for(k = 0; k < nv; k++)
    //     {
    //         w = av[k].v;
    //         ///check all out-nodes that has already been colored
    //         if(bi_g->seq[w].len != UNVISIT)
    //         {
    //             ///if this is the first colored node
    //             if(a_color == UNVISIT)
    //             {
    //                 a_color = bi_g->seq[w].len;
    //             }///if this is not
    //             else if(a_color != bi_g->seq[w].len)
    //             {
    //                 a_color = ISO;
    //                 bi_g->seq[v].len = ISO;
    //                 break;
    //             }
    //         }
    //     }

    //     if(a_color == ISO) continue;
    //     ///no colored out-node
    //     if(a_color == UNVISIT)
    //     {
    //         bi_g->seq[v].len = RED;
    //     }
    //     else
    //     {
    //         bi_g->seq[v].len = 1 - a_color;
    //     }



        

    //     ///check if all reachable nodes are end-contig
    //     kdq_push(uint32_t, buf, i);
    //     while (kdq_size(buf) != 0)
    //     {
    //         v = *(kdq_pop(uint32_t, buf));
    //         ///find a non-end contig
    //         if((a.a[v] & (uint64_t)(0x100000000)) == 0)
    //         {
    //             bi_g->seq[i].len = UNVISIT;
    //             kdq_size(buf) = 0;
    //             break;
    //         }

    //         nv = asg_arc_n(bi_g, v);
    //         av = asg_arc_a(bi_g, v);
    //         for(k = 0; k < nv; k++)
    //         {
    //             w = av[k].v;
    //             ///if the out-node has not been visited
    //             if(bi_g->seq[w].len == UNVISIT)
    //             {
    //                 bi_g->seq[v].len = LABLE;
    //                 kdq_push(uint32_t, buf, w);
    //             }
    //         }
    //     }

    //     for (v = 0; v < bi_g->n_seq; v++)
    //     {
    //         if(bi_g->seq[v].len == LABLE) bi_g->seq[v].len = UNVISIT;
    //     }

    //     if(bi_g->seq[i].len == UNVISIT)
    //     {
    //         continue;
    //     }






    //     ///now all reachable nodes of i are end-contigs
    //     kdq_push(uint32_t, buf, i);
    //     while (kdq_size(buf) != 0)
    //     {
    //         ///each v here is the uncolored node in the first round
    //         v = *(kdq_pop(uint32_t, buf));
    //         if(bi_g->seq[v].len == ISO) continue;
    //         nv = asg_arc_n(bi_g, v);
    //         av = asg_arc_a(bi_g, v);


    //         ///get all out-nodes of v
    //         for(k = 0; k < nv; k++)
    //         {
    //             w = av[k].v;
    //             if(bi_g->seq[w].len == bi_g->seq[v].len)
    //             {
    //                 break;
    //             }
    //         }
    //         ///means v is conflict 
    //         if(k != nv)
    //         {
    //             bi_g->seq[v].len = ISO;
    //             continue;
    //         }



    //         ///if v is not conflict with others
    //         for(k = 0; k < nv; k++)
    //         {
    //             w = av[k].v;
    //             ///if the out-node has not been visited
    //             if(bi_g->seq[w].len == UNVISIT)
    //             {   
    //                 bi_g->seq[w].len = 1 - bi_g->seq[v].len;
    //                 kdq_push(uint32_t, buf, w);
    //             }
    //         }
    //     }
    // }

    kdq_destroy(uint32_t, buf);
}

void further_clean_untig_graph(ma_ug_t *ug, asg_t *read_g, ma_hit_t_alloc* reverse_sources,
R_to_U* ruIndex, buf_t* b_0, uint8_t* visit, float density, uint32_t miniLen, uint32_t bi_graph_Len,
uint32_t long_hap_overlap, float long_hap_overlap_rate, float lable_match_rate)
{
    asg_t *bi_g = NULL;
    bi_g = asg_init();
    kvec_t(uint64_t) a;
    kv_init(a);
    kvec_t_u64_warp u_vecs;
    kv_init(u_vecs.a);
    asg_t* nsg = ug->g;
    uint32_t v, n_vtx = nsg->n_seq * 2, beg = 0, end, uId, cId, rId, i, j, k, self_offset, self_label_offset;
    uint64_t uInfor = 0;
    ma_utg_t* reads;
    ///int flag;
    memset(visit, 0, nsg->n_seq);
    /******************************set ruIndex*********************************/
    for (v = 0; v < n_vtx; ++v) 
    {
        beg = v;
        if(nsg->seq[v>>1].c == ALTER_LABLE || nsg->seq[beg>>1].del || visit[beg>>1])
        {
            continue;
        }

        if(get_real_length(nsg, beg, NULL)<=0 && get_real_length(nsg, beg^1, NULL)>0)
        {
            continue;
        }


        if(get_real_length(nsg, beg^1, NULL) == 1)
        {
            get_real_length(nsg, beg^1, &end);
            if(get_real_length(nsg, end^1, NULL) == 1)
            {
                continue;
            }
        }

        b_0->b.n = 0;ug->u.a[beg>>1].circ = 0;
        get_long_tip_length(nsg, &(ug->u), beg, &end, b_0);

        uInfor = 0;
        //scan all untigs
        for (i = 0; i < b_0->b.n; i++)
        {
            uId = b_0->b.a[i]>>1;
            visit[uId] = 1;
            reads = &(ug->u.a[uId]);
            uInfor += reads->n;
        }
        /****************************may have bugs********************************/
        ///uInfor = uInfor << 32; uInfor = uInfor | (uint64_t)beg;
        uInfor = uInfor << 33; uInfor = uInfor | (uint64_t)beg;
        /****************************may have bugs********************************/
        kv_push(uint64_t, a, uInfor);
    }


    radix_sort_arch64(a.a, a.a + a.n);
    for (i = 0; i < (a.n>>1); ++i) 
    {
        uInfor = a.a[i];
        a.a[i] = a.a[a.n - i - 1];
        a.a[a.n - i - 1] = uInfor;
    }


    for (v = 0; v < a.n; v++)
    {

        ///all untig ID of this contig
        beg = (uint32_t)(a.a[v]);
        b_0->b.n = 0;ug->u.a[beg>>1].circ = 0;
        get_long_tip_length(nsg, &(ug->u), beg, &end, b_0);
        cId = v;
        if(get_real_length(nsg, beg^1, NULL) == 0 && get_real_length(nsg, end, NULL) == 0)
        {
            a.a[v] = a.a[v] | (uint64_t)(0x100000000);
        }
        
        for (i = 0; i < b_0->b.n; i++)
        {
            uId = b_0->b.a[i]>>1;
            reads = &(ug->u.a[uId]);
            for (j = 0; j < reads->n; j++)
            {
                rId = reads->a[j]>>33;
                set_R_to_U(ruIndex, rId, cId, 1);
            }
        }
        
    }
    /******************************set ruIndex*********************************/

    
    ///scan all contigs
    for (i = 0; i < a.n; i++)
    {
        ///all untig ID of this contig
        beg = (uint32_t)(a.a[i]);
        b_0->b.n = 0;ug->u.a[beg>>1].circ = 0;
        get_long_tip_length(nsg, &(ug->u), beg, &end, b_0);
        
        
        u_vecs.a.n = 0;
        ///scan all untigs
        for (j = 0, self_offset = 0, self_label_offset = 0; j < b_0->b.n; j++)
        {
            uId = b_0->b.a[j]>>1;
            ///if(IsMerge(ug->u, uId)>0) continue;
            reads = &(ug->u.a[uId]);
            ///scan all reads
            ///self_offset will skip fake(merge) nodes, but not skip HAP_LABLE
            for (k = 0; k < reads->n; k++, self_offset++)
            {
                rId = reads->a[k]>>33;
                ///if(read_g->seq[rId].c == HAP_LABLE) continue;
                query_reverse_sources(read_g, reverse_sources, ruIndex, rId, self_offset, &u_vecs);

                if(read_g->seq[rId].c == HAP_LABLE) self_label_offset++;
            }
        }

        if(self_offset >= long_hap_overlap && self_label_offset >= (self_offset*lable_match_rate))
        {
            asg_seq_set(bi_g, i, RED, 0);
        }
        else
        {
            asg_seq_set(bi_g, i, UNVISIT, 0);
        }

        ///fprintf(stderr, "self_label_offset: %u, self_offset: %u\n", self_label_offset, self_offset);
        
        get_useful_contig_advance(ug, read_g, reverse_sources, ruIndex, &u_vecs, b_0, a.a, 
        bi_g, i, density, long_hap_overlap, long_hap_overlap_rate);
    }

    asg_cleanup(bi_g);

    ///fprintf(stderr, "***********n_seq: %u, n_arc: %u\n", bi_g->n_seq, bi_g->n_arc);
    for (v = 0; v < bi_g->n_seq; v++)
    {
        uint32_t nv = asg_arc_n(bi_g, v), nw, w;
        asg_arc_t *av = asg_arc_a(bi_g, v), *aw;
        for (i = 0; i < nv; ++i)
        {
            w = av[i].v;
            if(w == v)
            {
                av[i].del = 1;
                continue;
            }
            nw = asg_arc_n(bi_g, w);
            aw = asg_arc_a(bi_g, w);
            for (j = 0; j < nw; j++)
            {
                if(aw[j].v == v) break;
            }

            if(j == nw) av[i].del = 1;
        }
    }

    asg_cleanup(bi_g);
    ///fprintf(stderr, "***********n_seq: %u, n_arc: %u\n", bi_g->n_seq, bi_g->n_arc);

    bi_paration(bi_g, a.a, bi_graph_Len);

    for (v = 0; v < bi_g->n_seq; v++)
    {
        beg = (uint32_t)(a.a[v]);
        b_0->b.n = 0;ug->u.a[beg>>1].circ = 0;
        get_long_tip_length(nsg, &(ug->u), beg, &end, b_0);

        if(bi_g->seq[v].len == BLACK)
        {
            for (i = 0; i < b_0->b.n; i++)
            {
                nsg->seq[(b_0->b.a[i])>>1].c = ALTER_LABLE;
            }

            for (i = 0; i < b_0->b.n; i++)
            {
                asg_seq_drop(nsg, ((b_0->b.a[i])>>1));
            }
        }
        







        // fprintf(stderr, "cId: %u, beg>>1: %u, end>>1: %u, b_0->b.n: %u, Len: %u, type: %u\n", 
        // v, beg>>1, end>>1, (uint32_t)b_0->b.n, (uint32_t)(a.a[v]>>33), bi_g->seq[v].len);
        // uint32_t nv = asg_arc_n(bi_g, v), w;
        // asg_arc_t *av = asg_arc_a(bi_g, v);
        // for (i = 0; i < nv; ++i)
        // {
        //     w = av[i].v;
        //     fprintf(stderr, "w: %u\n", w);
        // }
    }

    
    /**
    for (i = 1; i < a.n; i++)
    {
        if((a.a[i]>>33) > (a.a[i-1]>>33))
        {
            fprintf(stderr, "hehe\n");
        }
    }

    
    for (i = 0; i < a.n; i++)
    {
        uint32_t k, get_cId, tLen = 0, is_Unitig;
        beg = (uint32_t)(a.a[i]);
        b_0->b.n = 0;ug->u.a[beg>>1].circ = 0;
        get_long_tip_length(nsg, &(ug->u), beg, &end, b_0);
        for (j = 0; j < b_0->b.n; j++)
        {
            uId = b_0->b.a[j]>>1;
            reads = &(ug->u.a[uId]);
            tLen += reads->n;
            for (k = 0; k < reads->n; k++)
            {
                rId = reads->a[k]>>33;
                get_R_to_U(ruIndex, rId, &get_cId, &is_Unitig);
                if(is_Unitig != 1 || get_cId != i)
                {
                    fprintf(stderr, "###is_Unitig: %u, get_cId: %u, i: %u\n", 
                    is_Unitig, get_cId, i);
                }
            }
        }

        uint32_t qLen = 0, m;
        for (m = 0; m < ruIndex->len; m++)
        {
            get_R_to_U(ruIndex, m, &get_cId, &is_Unitig);
            if(get_cId == (uint32_t)-1 || is_Unitig != 1) continue;
            if(get_cId == i)
            {
                qLen = 0;
                for (j = 0; j < b_0->b.n; j++)
                {
                    uId = b_0->b.a[j]>>1;
                    reads = &(ug->u.a[uId]);
                    tLen += reads->n;
                    for (k = 0; k < reads->n; k++)
                    {
                        rId = reads->a[k]>>33;
                        if(m == rId)
                        {
                            qLen = 1;
                            goto found;
                        }
                    }
                }
                found:;
                if(qLen == 0)
                {
                    fprintf(stderr, "***is_Unitig: %u, get_cId: %u, m: %u\n", 
                    is_Unitig, get_cId, m);
                }
            }
        }
    }
    **/

    asg_cleanup(nsg);
    free(a.a);
    kv_destroy(u_vecs.a);
    asg_destroy(bi_g);
}

void clean_untig_graph(ma_ug_t *ug, asg_t *read_g, ma_hit_t_alloc* reverse_sources,
long long bubble_dist, long long tipsLen, float tip_drop_ratio, long long stops_threshold, 
R_to_U* ruIndex, buf_t* b_0, uint8_t* visit, float density, uint32_t miniHapLen, 
uint32_t miniBiGraph, float chimeric_rate, int is_final_clean)
{
    asg_t *g = ug->g;
    asg_cut_tip_primary(g, ug, tipsLen);
    long long pre_cons = get_graph_statistic(g);
    long long cur_cons = 0;
    while(pre_cons != cur_cons)
    {
        pre_cons = get_graph_statistic(g);
        ///need consider tangles
        asg_pop_bubble_primary(g, bubble_dist);
        ///need consider tangles
        asg_arc_cut_long_tip_primary(g, ug, tip_drop_ratio);
        ///need consider tangles
        ///note we need both the read graph and the untig graph
        untig_asg_arc_cut_long_equal_tips_assembly(ug, read_g, reverse_sources, 2, ruIndex);
        untig_asg_arc_cut_long_tip_primary_complex(ug, tip_drop_ratio, stops_threshold);
        untig_asg_arc_cut_long_equal_tips_assembly_complex(ug, read_g, reverse_sources, 2, 
        stops_threshold, ruIndex);
        if(is_final_clean)
        {
            untig_asg_arc_cut_chimeric(ug, read_g, reverse_sources, 2, stops_threshold, chimeric_rate,
            ruIndex);
        }
        cur_cons = get_graph_statistic(g);
    }

    asg_cut_tip_primary(g, ug, tipsLen);
    untig_asg_arc_simple_large_bubbles(ug, read_g, reverse_sources, 2, ruIndex);

    if(is_final_clean)
    {
        lable_hap_asg_by_ug(ug, read_g);
        further_clean_untig_graph(ug, read_g, reverse_sources, ruIndex, b_0, visit, density, miniHapLen,
        miniBiGraph, 200, 0.4, 0.5);
        adjust_asg_by_ug(ug, read_g);
    }
    
}


void clean_untig_graph_bubbles(ma_ug_t *ug, asg_t *read_g, ma_hit_t_alloc* reverse_sources,
long long bubble_dist, long long tipsLen, float tip_drop_ratio, long long stops_threshold, 
R_to_U* ruIndex, buf_t* b_0, uint8_t* visit, float density, uint32_t miniHapLen, 
uint32_t miniBiGraph, float chimeric_rate, int is_final_clean)
{
    asg_t *g = ug->g;
    asg_cut_tip_primary(g, ug, tipsLen);
    long long pre_cons = get_graph_statistic(g);
    long long cur_cons = 0;
    while(pre_cons != cur_cons)
    {
        pre_cons = get_graph_statistic(g);
        ///need consider tangles
        asg_pop_bubble_primary(g, bubble_dist);
        cur_cons = get_graph_statistic(g);
    }

    asg_cut_tip_primary(g, ug, tipsLen);
    untig_asg_arc_simple_large_bubbles(ug, read_g, reverse_sources, 2, ruIndex);

    if(is_final_clean)
    {
        lable_hap_asg_by_ug(ug, read_g);
        adjust_asg_by_ug(ug, read_g);
    }
    
}



void resolve_simple_case(ma_ug_t *ug, asg_t* nsg)
{
    uint32_t v, n_vtx = nsg->n_seq * 2, nw, nv, w1, w2, i;
    asg_arc_t *av, *aw;
    ma_utg_t *v_x = NULL, *v_y = NULL;

    for (v = 0; v < n_vtx; ++v) 
    {
        if (nsg->seq[v>>1].del || nsg->seq[v>>1].c == ALTER_LABLE) continue;
        if(asg_arc_n(nsg, v) < 1 || asg_arc_n(nsg, v^1) < 1) continue;
        if(get_real_length(nsg, v, NULL) != 1 || get_real_length(nsg, v^1, NULL) != 1) continue;
        get_real_length(nsg, v, &w1); get_real_length(nsg, v^1, &w2);

        ///for simple circle
        if(w1 == (w2^1))
        {
            ///fprintf(stderr, "* v>>1: %u\n", v>>1);
            EvaluateLen(ug->u, w1>>1) = EvaluateLen(ug->u, w1>>1) + EvaluateLen(ug->u, v>>1);
            ///nsg->seq[w1>>1].len = nsg->seq[w1>>1].len + nsg->seq[v>>1].len;

            v_x = &(ug->u.a[w1>>1]);
            v_y = &(ug->u.a[v>>1]);
            append_ma_utg_t(v_x, v_y);
            asg_seq_del(nsg, v>>1);
        }
        else if(w1 == w2)
        {
            if(get_real_length(nsg, w1^1, NULL) != 2) continue;
            if(get_real_length(nsg, w1, NULL) != 1 && get_real_length(nsg, w1, NULL) != 2) continue;

            ///fprintf(stderr, "# v>>1: %u\n", v>>1);
            EvaluateLen(ug->u, w1>>1) = EvaluateLen(ug->u, w1>>1) + EvaluateLen(ug->u, v>>1);
            ///nsg->seq[w1>>1].len = nsg->seq[w1>>1].len + nsg->seq[v>>1].len;
            
            v_x = &(ug->u.a[w1>>1]);
            v_y = &(ug->u.a[v>>1]);

            asg_seq_del(nsg, v>>1);
            if(get_real_length(nsg, w1, NULL) == 1) continue;

            aw = asg_arc_a(nsg, w1);
            nw = asg_arc_n(nsg, w1);
            av = asg_arc_a(nsg, w1^1);
            for (i = 0; i < nw; i++)
            {
                if(!aw[i].del)
                {
                    av[0].del = 0;aw[i].del = 1;
                    av[0].el = aw[i].el;
                    av[0].no_l_indel = aw[i].no_l_indel;
                    av[0].ol = aw[i].ol;
                    av[0].strong = aw[i].strong;
                    av[0].v = aw[i].v;
                    av[0].ul = (aw[i].ul)^(0x100000000);


                    av = asg_arc_a(nsg, aw[i].v^1);
                    nv = asg_arc_n(nsg, aw[i].v^1);
                    uint32_t k = 0;
                    for (k = 0; k < nv; k++)
                    {
                        if(av[k].v == (aw[i].ul>>32^1))
                        {
                            av[k].v = av[k].v^1;
                            break;
                        }
                    }
                    
                    break;
                }
            }

        }
    }
}

void print_node(asg_t* g, ma_ug_t *ug)
{
    int input_iv;
    uint32_t iv, v;
    asg_arc_t *av, *aw;
    uint32_t nv, nw, i, k, w;
    while (1)
    {
        fprintf(stderr, "\n\ninput v: ");
        if(scanf("%d", &input_iv) == 0) break;
        if(input_iv == -1) break;
        iv = input_iv;
        iv = iv << 1;

        v = iv;
        av = asg_arc_a(g, v);
        nv = asg_arc_n(g, v);
        fprintf(stderr, "0**************v>>1: %u, v&1: %u, c: %u, del: %u, len: %u**************\n", 
        v>>1, v&1, (uint32_t)g->seq[v>>1].c, (uint32_t)g->seq[v>>1].del, g->seq[v>>1].len);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            w = av[i].v;
            fprintf(stderr, "(%u) w>>1: %u, w&1: %u, ol: %u, eLen: %u\n", 
            i, w>>1, w&1, av[i].ol, (uint32_t)av[i].ul);
            

            aw = asg_arc_a(g, w^1);
            nw = asg_arc_n(g, w^1);
            for (k = 0; k < nw; k++)
            {
                if(aw[k].del) continue;
                if(aw[k].v == (v^1)) break;
            }

            fprintf(stderr, "self>>1: %u, self&1: %u, out>>1: %u, out^1: %u, is_sym: %u\n", 
            (uint32_t)(av[i].ul>>33), (uint32_t)((av[i].ul>>32)&1), av[i].v>>1, av[i].v&1, 
            (uint32_t)(k != nw));

        }






        v = iv^1;
        av = asg_arc_a(g, v);
        nv = asg_arc_n(g, v);
        fprintf(stderr, "1**************v>>1: %u, v&1: %u, c: %u, del: %u, len: %u**************\n", 
        v>>1, v&1, (uint32_t)g->seq[v>>1].c, (uint32_t)g->seq[v>>1].del, (uint32_t)g->seq[v>>1].len);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            w = av[i].v;
            fprintf(stderr, "w>>1: %u, w&1: %u, ol: %u, eLen: %u\n", 
            w>>1, w&1, (uint32_t)av[i].ol, (uint32_t)av[i].ul);
            

            aw = asg_arc_a(g, w^1);
            nw = asg_arc_n(g, w^1);
            for (k = 0; k < nw; k++)
            {
                if(aw[k].del) continue;
                if(aw[k].v == (v^1)) break;
            }

            fprintf(stderr, "self>>1: %u, self&1: %u, out>>1: %u, out^1: %u, is_sym: %u\n", 
            (uint32_t)(av[i].ul>>33), (uint32_t)((av[i].ul>>32)&1), av[i].v>>1, av[i].v&1, 
            (uint32_t)(k != nw));

        }

        if(ug != NULL)
        {
            ma_utg_t *ma_v = &(ug->u.a[v>>1]);
            fprintf(stderr, "\n###\nma_v->n: %u, ma_v->len: %u\n", ma_v->n, ma_v->len);
            fprintf(stderr, "start>>1: %u, start&1: %u, end>>1: %u, end&1: %u\n", 
            ma_v->start>>1, ma_v->start&1, ma_v->end>>1, ma_v->end&1);
            for (i = 0; i < ma_v->n; i++)
            {
                v = ma_v->a[i]>>32;
                fprintf(stderr, "i: %u, v>>1: %u, v&1: %u, len: %u\n", 
                i, v>>1, v&1, (uint32_t)ma_v->a[i]);
            }
            
        }
    }
}

void label_tangles(asg_t *sg, ma_hit_t_alloc* reverse_sources, long long minLongUntig, 
long long maxShortUntig, float l_untig_rate, float max_node_threshold, long long bubble_dist,
long long tipsLen, float tip_drop_ratio, long long stops_threshold, R_to_U* ruIndex, int just_bubble)
{
    double startTime = Get_T();
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));
    uint32_t v, sv, n_vtx, beg, end, next_uID = (uint32_t)-1;
    kvec_t_u32_warp u_vecs;
    kv_init(u_vecs.a);
    ma_ug_t *ug = NULL;
    ug = ma_ug_gen_primary(sg, PRIMARY_LABLE);

    ///for each untig, all node have the same direction
    ///and all node except the last one just have one edge
    ///the last one may have multiple edges
    ///for the useful untig, the signal is (u->n >= LongUntigThreshold && !u->circ)
    asg_t* nsg = ug->g;
    n_vtx = nsg->n_seq;
    for (v = 0; v < n_vtx; ++v) 
    {
        nsg->seq[v].c = PRIMARY_LABLE;
        EvaluateLen(ug->u, v) = ug->u.a[v].n;
        IsMerge(ug->u, v) = 0;
    }

    resolve_simple_case(ug, nsg);
    asg_cleanup(nsg);
    asg_symm(nsg);

    ///print_node(nsg);

    if(just_bubble)
    {
        clean_untig_graph_bubbles(ug, sg, reverse_sources, bubble_dist, tipsLen, 
        tip_drop_ratio, stops_threshold, ruIndex, &b_0, NULL, 0.8, 20, 200, 0.05, 0);
    }
    else
    {
        clean_untig_graph(ug, sg, reverse_sources, bubble_dist, tipsLen, 
        tip_drop_ratio, stops_threshold, ruIndex, &b_0, NULL, 0.8, 20, 200, 0.05, 0);
    }
    
    

    uint8_t* visit = NULL;
    visit = (uint8_t*)malloc(sizeof(uint8_t) * nsg->n_seq);
    uint32_t n_reduce, flag;
    n_vtx = nsg->n_seq * 2;
    while (1)
    {
        n_reduce = 0;
        for (v = 0; v < n_vtx; ++v) 
        {
            //as for return value: 0: do nothing, 1: unroll, 2: convex
            //we just need 1
            sv = v;
            flag = 0;
            while (1)
            {
                flag = walk_through(sg, ug, reverse_sources, minLongUntig, 
                    maxShortUntig, l_untig_rate, max_node_threshold, &b_0, &b_1, 
                    &u_vecs, visit, sv, &beg, &end, &next_uID, ruIndex);
                n_reduce += flag;
                if(flag != UNROLL_M)
                {
                    break;
                }
            }
        }
        if(n_reduce == 0) break;
    }

    asg_cleanup(nsg);
	asg_symm(nsg);

    if(just_bubble)
    {
        clean_untig_graph_bubbles(ug, sg, reverse_sources, bubble_dist, tipsLen, 
        tip_drop_ratio, stops_threshold, ruIndex, &b_0, visit, 0.8, 20, 200, 0.05, 1);
    }
    else
    {
        clean_untig_graph(ug, sg, reverse_sources, bubble_dist, tipsLen, 
        tip_drop_ratio, stops_threshold, ruIndex, &b_0, visit, 0.8, 20, 200, 0.05, 1);
    }
    




    kv_destroy(u_vecs.a);
    ma_ug_destroy(ug);


    free(visit);
    free(b_0.b.a);
    free(b_1.b.a);
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2f s\n", __func__, Get_T()-startTime);
    }
}


uint32_t copy_ug_node(ma_ug_t *ug, asg_t* nsg, uint32_t v)
{
    ma_utg_t *p;
    ma_utg_t *o = &(ug->u.a[v]);
    uint32_t n_node = nsg->n_seq;
    asg_seq_set(nsg, n_node, nsg->seq[v].len, 0);
    kv_pushp(ma_utg_t, ug->u, &p);

    p->s = 0, p->start = o->start, p->end = o->end, p->len = o->len;
    p->n = o->n, p->circ = o->circ, p->m = o->m;
    p->a = (uint64_t*)malloc(8 * p->m);
    memcpy(p->a, o->a, (8*p->m));

    return n_node;
}


///v and w here have directions
uint32_t collect_ma_utg_ts(ma_ug_t *ug, uint32_t v, uint32_t w, ma_utg_t* result)
{
    uint32_t des_dir = w&1;
    long long i;
    uint64_t* aim = NULL;
    ma_utg_t *ma_w = &(ug->u.a[w>>1]);

    if(v == (uint32_t)-1)
    {
        result->start = ma_w->start;
        result->circ = ma_w->circ;
        result->end = ma_w->end;
        result->len = ma_w->len;
        result->n = 0;
    }

    
    if(result->m < (result->n + ma_w->n))
    {
        result->m = (result->n + ma_w->n);
        result->a = (uint64_t*)realloc(result->a, result->m * sizeof(uint64_t));
    }

    aim = result->a + result->n;

    if(des_dir == 1)
    {
        for (i = 0; i < (long long)ma_w->n; i++)
        {
            aim[ma_w->n - i - 1] = (ma_w->a[i])^(uint64_t)(0x100000000);
        }
    }
    else
    {
        for (i = 0; i < (long long)ma_w->n; i++)
        {
            aim[i] = ma_w->a[i];
        }
    }

    result->n = result->n + ma_w->n;
    return 1;
}


void debug_utg_graph(ma_ug_t *ug, asg_t* read_g, int test_tangle)
{
    asg_t* nsg = ug->g;
    uint32_t n_vtx = nsg->n_seq, i, j, k, l, totalLen, v, nv, nw, w, untig_v, rid_v;
    asg_arc_t *aw = NULL, *av = NULL;
    for (i = 0; i < n_vtx; i++)
    {
        if(ug->g->seq[i].del) continue;

        totalLen = 0;
        ma_utg_t* result = &(ug->u.a[i]);
        if(result->n == 0) continue;
        v = (uint64_t)(result->a[0])>>32;
        if(result->start != UINT32_MAX && result->start != v) fprintf(stderr, "hehe\n");
        v = (uint64_t)(result->a[result->n-1])>>32;
        if(result->end != UINT32_MAX && result->end != (v^1)) fprintf(stderr, "haha\n");
        for (j = 0; j < result->n-1; j++)
        {
            v = (uint64_t)(result->a[j])>>32;
            w = (uint64_t)(result->a[j+1])>>32; 

            av = asg_arc_a(read_g, v);
            nv = asg_arc_n(read_g, v);
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
            if(k == nv) fprintf(stderr ,"******error, j: %u, k: %u, nv: %u\n", j, k, nv);
            if(l != (uint32_t)(result->a[j])) fprintf(stderr ,"ERROR Length\n");
            totalLen = totalLen + l;
        }


        if(result->start != UINT32_MAX && result->end != UINT32_MAX && j < result->n)
        {
            v = (uint64_t)(result->a[j])>>32;
            l = read_g->seq[v>>1].len;
            if(l != (uint32_t)(result->a[j])) fprintf(stderr ,"*** ERROR Length\n");
            totalLen = totalLen + l;
        }

        if(result->start != UINT32_MAX && result->end != UINT32_MAX && totalLen != result->len)
        {
            fprintf(stderr ,"ERROR Total Length\n");
        } 
                


        if(ug->u.a[i].start == UINT32_MAX && ug->u.a[i].end == UINT32_MAX) continue;

        v = i<<1; v=v^1;
        av = asg_arc_a(ug->g, v);
        nv = asg_arc_n(ug->g, v);

        w = (ug->u.a[v>>1].start^1);
        aw = asg_arc_a(read_g, w);
        nw = asg_arc_n(read_g, w);


        if(get_real_length(ug->g, v, NULL) != get_real_length(read_g, w, NULL))
        {
            fprintf(stderr, "#########ERROR: i: %u, nv: %u, nw: %u\n", i, nv, nw);
        }

        for (j = 0; j < nv; j++)
        {
            if(av[j].del) continue;
            untig_v = av[j].v;
            if(untig_v&1) rid_v = ug->u.a[untig_v>>1].end;
            else rid_v = ug->u.a[untig_v>>1].start;
            
            for (k = 0; k < nw; k++)
            {
                if(aw[k].del) continue;
                if(aw[k].v == rid_v) break;
            }

            if(k == nw) fprintf(stderr, "#########ERROR: i: %u\n", i);
            if((k != nw) && (av[j].ol != aw[k].ol))
            {
                fprintf(stderr, "#########????????ERROR\n");
                fprintf(stderr, "nv: %u, nw: %u\n", nv, nw);
                fprintf(stderr, "av[%u].ol: %u, aw[%u].ol: %u, untig_v>>1: %u, untig_v&1: %u\n", 
                j, av[j].ol, k, aw[k].ol, untig_v>>1, untig_v&1);
            }
        }






        v = v^1;
        av = asg_arc_a(ug->g, v);
        nv = asg_arc_n(ug->g, v);

        w = (ug->u.a[v>>1].end^1);
        aw = asg_arc_a(read_g, w);
        nw = asg_arc_n(read_g, w);
        if(get_real_length(ug->g, v, NULL) != get_real_length(read_g, w, NULL))
        {
            fprintf(stderr, "*******ERROR: i: %u, nv: %u, nw: %u\n", i, nv, nw);
        }
        for (j = 0; j < nv; j++)
        {
            if(av[j].del) continue;
            untig_v = av[j].v;
            if(untig_v&1) rid_v = ug->u.a[untig_v>>1].end;
            else rid_v = ug->u.a[untig_v>>1].start;
            
            for (k = 0; k < nw; k++)
            {
                if(aw[k].del) continue;
                if(aw[k].v == rid_v) break;
            }

            if(k == nw) fprintf(stderr, "***********ERROR: i: %u\n", i);
            if((k != nw) && (av[j].ol != aw[k].ol))
            {
                fprintf(stderr, "***********????????ERROR\n");
                fprintf(stderr, "nv: %u, nw: %u\n", nv, nw);
                fprintf(stderr, "av[%u].ol: %u, aw[%u].ol: %u, untig_v>>1: %u, untig_v&1: %u\n", 
                j, av[j].ol, k, aw[k].ol, untig_v>>1, untig_v&1);
            }
        }

    }


    if(test_tangle != 1) return;
    fprintf(stderr, "test_tangle: %u\n", test_tangle);
    n_vtx = nsg->n_seq * 2;
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t w1, w2, beg, end, rnw;
        if(nsg->seq[v>>1].del) continue;
        if(asg_arc_n(nsg, v) < 1 || asg_arc_n(nsg, v^1) < 1) continue;
        if(get_real_length(nsg, v, NULL) != 1 || get_real_length(nsg, v^1, NULL) != 1) continue;
        get_real_length(nsg, v, &w1); get_real_length(nsg, v^1, &w2);

        ///for simple circle
        if(w1 == (w2^1))
        {
            beg = end = 0;
            aw = asg_arc_a(nsg, w1^1);
            nw = asg_arc_n(nsg, w1^1);
            for (i = 0, rnw = 0; i < nw; i++)
            {
                if(aw[i].del) continue;
                rnw++;
                if(aw[i].v == (v^1)) continue;
                beg = aw[i].v;
            }
            if(rnw != 2) continue;

            aw = asg_arc_a(nsg, w2^1);
            nw = asg_arc_n(nsg, w2^1);
            for (i = 0, rnw = 0; i < nw; i++)
            {
                if(aw[i].del) continue;
                rnw++;
                if(aw[i].v == v) continue;
                end = aw[i].v;
            }
            if(rnw != 2) continue;

            if(get_real_length(nsg, beg^1, NULL)!=1) continue;
            if(get_real_length(nsg, end^1, NULL)!=1) continue;
            if((beg>>1) == (end>>1)) continue;
            fprintf(stderr, "\n************\n");
        }
        else if(w1 == w2)
        {
            if(get_real_length(nsg, w1^1, NULL) != 2) continue;
            if(get_real_length(nsg, w1, NULL) != 2) continue;
            end = beg = (uint32_t)-1;

            aw = asg_arc_a(nsg, w1);
            nw = asg_arc_n(nsg, w1);
            for (i = 0, rnw = 0; i < nw && rnw < 2; i++)
            {
                if(aw[i].del) continue;
                if(rnw == 0) beg = aw[i].v;
                if(rnw == 1) end = aw[i].v;
                rnw++;
            }
            if((beg>>1) == (end>>1)) continue;
            if(get_real_length(nsg, beg^1, NULL)!=1) continue;
            if(get_real_length(nsg, end^1, NULL)!=1) continue;
            fprintf(stderr, "\n#################\n");
        }
    }
}

///just merge, don't delete anything
void merge_ug_nodes(ma_ug_t *ug, asg_t* read_g, kvec_t_u64_warp* array)
{
    if(array->a.n == 0) return;
    uint32_t i, k, v, w;
    ma_utg_t result;
    memset(&result, 0, sizeof(ma_utg_t));
    v = w = (uint32_t)-1;
    for (i = 0; i < array->a.n; i++)
    {
        w = array->a.a[i];
        collect_ma_utg_ts(ug, v, w, &result);
        v = w;
    }

    

    if(result.n == 0) return;
    asg_arc_t *av = NULL;
    uint32_t nv, l;
    result.len = 0;
    for (i = 0; i < result.n - 1; i++)
    {
        
        v = (uint64_t)(result.a[i])>>32;
        w = (uint64_t)(result.a[i + 1])>>32; 
        av = asg_arc_a(read_g, v);
        nv = asg_arc_n(read_g, v);
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
        if(k == nv) fprintf(stderr ,"******error, i: %u, k: %u, nv: %u\n", i, k, nv);
        result.a[i] = v; result.a[i] = result.a[i]<<32; result.a[i] = result.a[i] | (uint64_t)(l);
        result.len += l;
    }
    

    if(i < result.n)
    {
        v = (uint64_t)(result.a[i])>>32;
        l = read_g->seq[v>>1].len;
        result.a[i] = v;
        result.a[i] = result.a[i]<<32;
        result.a[i] = result.a[i] | (uint64_t)(l);
        result.len += l;
    }
    //has already set result.a, result.len, result.n, result.m
    result.circ = 0;
    result.start = result.a[0]>>32;
    result.end = (result.a[result.n-1]>>32)^1;

    
    uint32_t beg_uid = array->a.a[0];
    uint32_t end_uid = array->a.a[array->a.n-1];
    uint32_t realLen = array->a.n;
    uint32_t new_uid = array->a.a[0]>>1;
    uint64_t kmp;





    /*******************************just for debug**********************************/
    /**
    if(beg_uid&1)
    {
        if(result.start != ug->u.a[beg_uid>>1].end)
        {
            fprintf(stderr, "ERROR\n");
        }
    }
    else
    {
        if(result.start != ug->u.a[beg_uid>>1].start)
        {
            fprintf(stderr, "ERROR\n");
        }
    }

    if(end_uid&1)
    {
        if(result.end != ug->u.a[end_uid>>1].start)
        {
            fprintf(stderr, "ERROR\n");
        }
    }
    else
    {
        if(result.end != ug->u.a[end_uid>>1].end)
        {
            fprintf(stderr, "ERROR\n");
        }
    }
    **/
    /*******************************just for debug**********************************/

    asg_arc_t *aw = NULL;
    uint32_t nw = 0;
    ///corresponding to direction 1 of new node
    v = beg_uid^1;
    av = asg_arc_a(ug->g, v);
    nv = asg_arc_n(ug->g, v);
    ///fprintf(stderr, "beg_uid_v>>1: %u, beg_uid_v&1: %u, nv: %u\n", v>>1, v&1, nv);
    for (k = 0; k < nv; k++)
    {
        if(av[k].del) continue;
        kmp = new_uid<<1; kmp = kmp^1; kmp = kmp << 32;

        ///if((av[k].v>>1) == (end_uid>>1)) continue;
        w = av[k].v^1; aw = asg_arc_a(ug->g, w); nw = asg_arc_n(ug->g, w);
        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;
            if(aw[i].v == (v^1)) break;
        }
        if(i == nw) fprintf(stderr, "ERROR\n");
        kmp = kmp | (uint64_t)(aw[i].ol);


        ///kmp = kmp | (uint64_t)(((ug->g)->idx[v]>>32) + k);/**kmp = kmp | av[k].v;**/
        ///here kmp is ul
        kv_push(uint64_t, array->a, kmp);
        kmp = av[k].ol; kmp = kmp<<32; kmp = kmp|(uint64_t)(av[k].v);
        ///here kmp is ol + v
        kv_push(uint64_t, array->a, kmp);
        ///fprintf(stderr, "*av[%u].v>>1: %u, v&1: %u, ol: %u\n", k, av[k].v>>1, av[k].v&1, av[k].ol);
    }

    ///corresponding to direction 0 of new node
    v = end_uid;
    av = asg_arc_a(ug->g, v);
    nv = asg_arc_n(ug->g, v);
    ///fprintf(stderr, "end_uid_v>>1: %u, end_uid_v&1: %u, nv: %u\n", v>>1, v&1, nv);
    for (k = 0; k < nv; k++)
    {
        if(av[k].del) continue;
        kmp = new_uid<<1; kmp = kmp << 32;

        w = av[k].v^1; 
        aw = asg_arc_a(ug->g, w); 
        nw = asg_arc_n(ug->g, w);
        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;
            if(aw[i].v == (v^1)) break;
        }
        if(i == nw) fprintf(stderr, "ERROR\n");
        kmp = kmp | (uint64_t)(aw[i].ol);


        ///here kmp is ul
        kv_push(uint64_t, array->a, kmp);
        kmp = av[k].ol; kmp = kmp<<32; kmp = kmp|(uint64_t)(av[k].v);
        ///here kmp is ol + v
        kv_push(uint64_t, array->a, kmp);
        ///fprintf(stderr, "#av[%u].v>>1: %u, v&1: %u, ol: %u\n", k, av[k].v>>1, av[k].v&1, av[k].ol);
    }

    
    ma_utg_t* tmp;
    for (i = 0; i < realLen; i++)
    {
        w = array->a.a[i];
        tmp = &(ug->u.a[w>>1]);
        if(tmp->m != 0)
        {
            tmp->circ = tmp->end = tmp->len = tmp->m = tmp->n = tmp->start = 0;
            free(tmp->a);
            tmp->a = NULL;
        }
        asg_seq_del(ug->g, w>>1);
    }
    
    ug->u.a[beg_uid>>1] = result;
    ug->g->seq[beg_uid>>1].del = 0;

    uint32_t oLen = 0;
    for (; i < array->a.n; i += 2)
    {
        v = array->a.a[i]>>32;
        w = (uint32_t)array->a.a[i+1];
        /****************************may have bugs********************************/
        ///may have bug here, if there is an edge between beg_uid and end_uid
        ///if(((w>>1) == (beg_uid>>1)) || ((w>>1) == (end_uid>>1))) continue;
        if(((w>>1) == (beg_uid>>1)) || ((w>>1) == (end_uid>>1))) w = v;
        /****************************may have bugs********************************/

        oLen = array->a.a[i+1]>>32;
        asg_append_edges_to_srt(ug->g, v, ug->u.a[v>>1].len, w, oLen, 0, 0, 0);
        oLen = (uint32_t)array->a.a[i];
        asg_append_edges_to_srt(ug->g, w^1, ug->u.a[w>>1].len, v^1, oLen, 0, 0, 0);
    }
}

void unroll_simple_case(ma_ug_t *ug, asg_t* read_g)
{
    asg_t* nsg = ug->g;
    uint32_t v, n_vtx = nsg->n_seq * 2, rnw, nw, w1, w2, beg, end, i;
    asg_arc_t *aw;
    kvec_t_u64_warp u_vecs;
    kv_init(u_vecs.a);

    for (v = 0; v < n_vtx; ++v) 
    {
        ///if (nsg->seq[v>>1].del || nsg->seq[v>>1].c == ALTER_LABLE) continue;
        if (nsg->seq[v>>1].del) continue;
        if(asg_arc_n(nsg, v) < 1 || asg_arc_n(nsg, v^1) < 1) continue;
        if(get_real_length(nsg, v, NULL) != 1 || get_real_length(nsg, v^1, NULL) != 1) continue;
        get_real_length(nsg, v, &w1); get_real_length(nsg, v^1, &w2);
        if((v>>1) == (w1>>1) || (v>>1) == (w2>>1)) continue;

        ///for simple circle
        if(w1 == (w2^1))
        {
            beg = end = 0;
            aw = asg_arc_a(nsg, w1^1);
            nw = asg_arc_n(nsg, w1^1);
            for (i = 0, rnw = 0; i < nw; i++)
            {
                if(aw[i].del) continue;
                rnw++;
                if(aw[i].v == (v^1)) continue;
                beg = aw[i].v;
            }
            if(rnw != 2) continue;

            aw = asg_arc_a(nsg, w2^1);
            nw = asg_arc_n(nsg, w2^1);
            for (i = 0, rnw = 0; i < nw; i++)
            {
                if(aw[i].del) continue;
                rnw++;
                if(aw[i].v == v) continue;
                end = aw[i].v;
            }
            if(rnw != 2) continue;

            if(get_real_length(nsg, beg^1, NULL)!=1) continue;
            if(get_real_length(nsg, end^1, NULL)!=1) continue;
            if((beg>>1) == (end>>1)) continue;
            ///fprintf(stderr, "\n***v>>1: %u\n", v>>1);
            u_vecs.a.n = 0;
            kv_push(uint64_t, u_vecs.a, beg^1);
            kv_push(uint64_t, u_vecs.a, w1);
            kv_push(uint64_t, u_vecs.a, v);
            kv_push(uint64_t, u_vecs.a, w1);
            kv_push(uint64_t, u_vecs.a, end);
            merge_ug_nodes(ug, read_g, &u_vecs);
        }
        else if(w1 == w2)
        {
            if(get_real_length(nsg, w1^1, NULL) != 2) continue;
            if(get_real_length(nsg, w1, NULL) != 2) continue;
            end = beg = (uint32_t)-1;

            aw = asg_arc_a(nsg, w1);
            nw = asg_arc_n(nsg, w1);
            for (i = 0, rnw = 0; i < nw && rnw < 2; i++)
            {
                if(aw[i].del) continue;
                if(rnw == 0) beg = aw[i].v;
                if(rnw == 1) end = aw[i].v;
                rnw++;
            }
            if((beg>>1) == (end>>1)) continue;
            if(get_real_length(nsg, beg^1, NULL)!=1) continue;
            if(get_real_length(nsg, end^1, NULL)!=1) continue;
            ///fprintf(stderr, "\n###v>>1: %u\n", v>>1);
            u_vecs.a.n = 0;
            kv_push(uint64_t, u_vecs.a, beg^1);
            kv_push(uint64_t, u_vecs.a, w1^1);
            kv_push(uint64_t, u_vecs.a, v);
            kv_push(uint64_t, u_vecs.a, w1);
            kv_push(uint64_t, u_vecs.a, end);
            merge_ug_nodes(ug, read_g, &u_vecs);
        }
    }

    kv_destroy(u_vecs.a);
}

void adjust_utg(asg_t *sg, ma_ug_t *ug)
{
    double startTime = Get_T();

    asg_t* nsg = ug->g;
    unroll_simple_case(ug, sg);
    asg_cleanup(nsg);
    asg_symm(nsg);

    ///debug_utg_graph(ug, sg, 1);

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2f s\n", __func__, Get_T()-startTime);
    }
}








asg_t *copy_graph(asg_t* src, int round)
{
    asg_t *rg;
	///just calloc
	rg = asg_init();
    uint64_t i, k;
    for (i = 0; i < src->n_seq; ++i) 
    {
        ///if a read has been deleted, should we still add them?
		asg_seq_set(rg, i, src->seq[i].len, src->seq[i].del);
        rg->seq[i].c = src->seq[i].c;
	}
    asg_cleanup(rg);

    
    uint32_t v, n_vtx = src->n_seq * 2, totalL = 0;;


    for (k = 0; k < (uint32_t)round; k++)
    {
        for (i = k; i < src->n_arc; i = i + round)
        {
            if(totalL%50000==0) fprintf(stderr, "totalL: %u, n_arc: %u\n", totalL, src->n_arc);
            totalL++;
            if(src->arc[i].del) continue;
            asg_append_edges_to_srt(rg, src->arc[i].ul>>32, 
            (uint32_t)(src->arc[i].ul) + src->arc[i].ol,
            src->arc[i].v, src->arc[i].ol, src->arc[i].strong, 
            src->arc[i].el, src->arc[i].no_l_indel);
        }
    }
    
    
    totalL = 0;
    for (v = 0; v < n_vtx; ++v)
    {
        if (src->seq[v>>1].del) continue;
        uint32_t nv = asg_arc_n(src, v);
        asg_arc_t *av = asg_arc_a(src, v);
        for (i = 0; i < nv; i++)
        {
            if(totalL%50000==0) fprintf(stderr, "-totalL: %u, n_arc: %u\n", totalL, src->n_arc);
            totalL++;
            if(av[i].del) continue;
            asg_append_edges_to_srt(rg, av[i].ul>>32, (uint32_t)(av[i].ul) + av[i].ol,
            av[i].v, av[i].ol, av[i].strong, av[i].el, av[i].no_l_indel);
        }
    }
    

    

    /**
    for (v = 0; v < n_vtx; ++v)
    {
        if(src->idx[v] != rg->idx[v])
        {
            fprintf(stderr, "*****v: %u, src_i: %u, srcLen: %u, rg_i: %u, rgLen: %u\n", 
            v, src->idx[v]>>32, (uint32_t)src->idx[v], 
            rg->idx[v]>>32, (uint32_t)rg->idx[v]); 
        }
    }
    **/


    for (v = 0; v < n_vtx; ++v)
    {
        
        if(src->seq[v>>1].del != rg->seq[v>>1].del)
        {
            fprintf(stderr, "ERROR1\n");
        }

        uint32_t nv = asg_arc_n(src, v);
        asg_arc_t *av = asg_arc_a(src, v);
        ///if(nv != rnv)
        if(get_real_length(src, v, NULL) != get_real_length(rg, v, NULL))
        {
            fprintf(stderr, "ERROR2\n");
        }

        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            asg_arc_t forward = av[i], backward;
            memset(&backward, 0, sizeof(backward));
            if(get_arc(rg, (forward.ul>>32), forward.v, &backward)!=1)
            {
                fprintf(stderr, "ERROR3\n");
            }

            if(forward.del != backward.del || forward.el != backward.el || 
               forward.no_l_indel != backward.no_l_indel || forward.ol != backward.ol ||
               forward.strong != backward.strong || forward.ul != backward.ul || forward.v != backward.v)
            {
                fprintf(stderr, "ERROR4\n");
                fprintf(stderr, "forward, del: %u, el: %u, no_l_indel: %u, ol: %u, strong: %u, ul: %u, v: %u\n",
                (uint32_t)forward.del, (uint32_t)forward.el, (uint32_t)forward.no_l_indel, 
                (uint32_t)forward.ol, (uint32_t)forward.strong,
                (uint32_t)forward.ul, (uint32_t)forward.v);
                fprintf(stderr, "backward, del: %u, el: %u, no_l_indel: %u, ol: %u, strong: %u, ul: %u, v: %u\n",
                (uint32_t)backward.del, (uint32_t)backward.el, (uint32_t)backward.no_l_indel, 
                (uint32_t)backward.ol, (uint32_t)backward.strong,
                (uint32_t)backward.ul, (uint32_t)backward.v);
            }
        }

        
        nv = asg_arc_n(rg, v);
        av = asg_arc_a(rg, v);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            asg_arc_t forward = av[i], backward;
            memset(&backward, 0, sizeof(backward));
            if(get_arc(src, (forward.ul>>32), forward.v, &backward)!=1)
            {
                fprintf(stderr, "ERROR5\n");
            }

            if(forward.del != backward.del || forward.el != backward.el || 
               forward.no_l_indel != backward.no_l_indel || forward.ol != backward.ol ||
               forward.strong != backward.strong || forward.ul != backward.ul || forward.v != backward.v)
            {
                fprintf(stderr, "ERROR6\n");
            }
        }
        
        
        ///get_arc(g, forward.v^1, (forward.ul>>32)^1, &backward) 
        ///fprintf(stderr, "+v: %u\n", v);
        
    }

    
    if(src->seq_vis)
    {
        rg->seq_vis = (uint8_t*)malloc(src->n_seq*2*sizeof(uint8_t));
        memcpy(rg->seq_vis, src->seq_vis, src->n_seq*2*sizeof(uint8_t));
    }
    

    fprintf(stderr, "end_copy\n");
    return rg;
}

int load_coverage_cut(ma_sub_t** coverage_cut, char* read_file_name)
{
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "r");
    if(!fp)
    {
        return 0;
    }
    int f_flag = 0;
    uint64_t n_read;
    f_flag += fread(&n_read, sizeof(n_read), 1, fp);
    (*coverage_cut) = (ma_sub_t*)malloc(sizeof(ma_sub_t)*n_read);

    uint64_t i = 0, tmp;
    for (i = 0; i < n_read; i++)
    {
        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*coverage_cut)[i].c = tmp;
        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*coverage_cut)[i].del = tmp;
        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*coverage_cut)[i].e = tmp;
        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*coverage_cut)[i].s = tmp;
    }
    free(index_name);    
    fflush(fp);
    fclose(fp);
    return 1;
}


int write_coverage_cut(ma_sub_t* coverage_cut, char* read_file_name, uint64_t n_read)
{
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "w");
    fwrite(&n_read, sizeof(n_read), 1, fp);
    uint64_t i = 0, tmp;
    for (i = 0; i < n_read; i++)
    {
        tmp = coverage_cut[i].c;
        fwrite(&tmp, sizeof(tmp), 1, fp);
        tmp = coverage_cut[i].del;
        fwrite(&tmp, sizeof(tmp), 1, fp);
        tmp = coverage_cut[i].e;
        fwrite(&tmp, sizeof(tmp), 1, fp);
        tmp = coverage_cut[i].s;
        fwrite(&tmp, sizeof(tmp), 1, fp);
    }
    free(index_name);    
    fflush(fp);
    fclose(fp);

    return 1;
}

int write_ruIndex(R_to_U* ruIndex, char* read_file_name)
{
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "w");
    fwrite(&ruIndex->len, sizeof(ruIndex->len), 1, fp);
    fwrite(ruIndex->index, sizeof(ruIndex->index[0]), ruIndex->len, fp);

    free(index_name);    
    fflush(fp);
    fclose(fp);

    return 1;
}

int load_ruIndex(R_to_U* ruIndex, char* read_file_name)
{
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "r");
    if(!fp)
    {
        return 0;
    }
    int f_flag = 0;
    f_flag += fread(&(ruIndex)->len, sizeof((ruIndex)->len), 1, fp);
    (ruIndex)->index = (uint32_t*)malloc(sizeof(uint32_t)*(ruIndex)->len);
    f_flag += fread((ruIndex)->index, sizeof((ruIndex)->index[0]), (ruIndex)->len, fp);

    free(index_name);    
    fflush(fp);
    fclose(fp);

    return 1;
}

int write_asg_t(asg_t *sg, char* read_file_name)
{
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "w");
    uint32_t tmp, i;
    
    tmp = sg->n_arc;
    fwrite(&tmp, sizeof(tmp), 1, fp);
    // tmp = sg->m_arc;
    fwrite(&tmp, sizeof(tmp), 1, fp);

    tmp = sg->is_srt;
    fwrite(&tmp, sizeof(tmp), 1, fp);


    tmp = sg->n_seq;
    fwrite(&tmp, sizeof(tmp), 1, fp);
    // tmp = sg->m_seq;
    fwrite(&tmp, sizeof(tmp), 1, fp);
    
    tmp = sg->is_symm;
    fwrite(&tmp, sizeof(tmp), 1, fp);
    tmp = sg->r_seq;
    fwrite(&tmp, sizeof(tmp), 1, fp);


    uint32_t Len;

    Len = sg->n_seq*2;
    fwrite(sg->seq_vis, sizeof(sg->seq_vis[0]), Len, fp);

    Len = sg->n_seq*2;
    fwrite(sg->idx, sizeof(sg->idx[0]), Len, fp);


    for (i = 0; i < sg->n_arc; i++)
    {
        tmp = sg->arc[i].del;
        fwrite(&tmp, sizeof(tmp), 1, fp);
        tmp = sg->arc[i].el;
        fwrite(&tmp, sizeof(tmp), 1, fp);
        tmp = sg->arc[i].no_l_indel;
        fwrite(&tmp, sizeof(tmp), 1, fp);
        tmp = sg->arc[i].ol;
        fwrite(&tmp, sizeof(tmp), 1, fp);
        tmp = sg->arc[i].strong;
        fwrite(&tmp, sizeof(tmp), 1, fp);

        uint64_t tmp_64;
        tmp_64 = sg->arc[i].ul;
        fwrite(&tmp_64, sizeof(tmp_64), 1, fp);

        tmp = sg->arc[i].v;
        fwrite(&tmp, sizeof(tmp), 1, fp);
    }


    for (i = 0; i < sg->n_seq; i++)
    {
        tmp = sg->seq[i].c;
        fwrite(&tmp, sizeof(tmp), 1, fp);
        tmp = sg->seq[i].del;
        fwrite(&tmp, sizeof(tmp), 1, fp);
        tmp = sg->seq[i].len;
        fwrite(&tmp, sizeof(tmp), 1, fp);
    }
    
    free(index_name);    
    fflush(fp);
    fclose(fp);

    return 1;
}


int load_asg_t(asg_t **sg, char* read_file_name)
{
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "r");
    if(!fp)
    {
        return 0;
    }
    uint32_t tmp, i;

    (*sg) = (asg_t*)calloc(1, sizeof(asg_t));
    int f_flag = 0;
    f_flag += fread(&tmp, sizeof(tmp), 1, fp);
    (*sg)->n_arc = tmp;
    f_flag += fread(&tmp, sizeof(tmp), 1, fp);
    (*sg)->m_arc = tmp;
    f_flag += fread(&tmp, sizeof(tmp), 1, fp);
    (*sg)->is_srt = tmp;
    f_flag += fread(&tmp, sizeof(tmp), 1, fp);
    (*sg)->n_seq = tmp;
    f_flag += fread(&tmp, sizeof(tmp), 1, fp);
    (*sg)->m_seq = tmp;
    f_flag += fread(&tmp, sizeof(tmp), 1, fp);
    (*sg)->is_symm = tmp;
    f_flag += fread(&tmp, sizeof(tmp), 1, fp);
    (*sg)->r_seq = tmp;

    uint32_t Len;

    Len = (*sg)->n_seq*2;
    (*sg)->seq_vis = (uint8_t*)malloc(sizeof(uint8_t)*Len);
    f_flag += fread((*sg)->seq_vis, sizeof((*sg)->seq_vis[0]), Len, fp);



    Len = (*sg)->n_seq*2;
    (*sg)->idx = (uint64_t*)malloc(sizeof(uint64_t)*Len);
    f_flag += fread((*sg)->idx, sizeof((*sg)->idx[0]), Len, fp);
    
    



    (*sg)->arc = (asg_arc_t*)malloc(sizeof(asg_arc_t)*(*sg)->m_arc);
    (*sg)->seq = (asg_seq_t*)malloc(sizeof(asg_seq_t)*(*sg)->m_seq);
    

    for (i = 0; i < (*sg)->n_arc; i++)
    {
        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*sg)->arc[i].del = tmp;

        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*sg)->arc[i].el = tmp;

        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*sg)->arc[i].no_l_indel = tmp;

        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*sg)->arc[i].ol = tmp;

        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*sg)->arc[i].strong = tmp;

        uint64_t tmp_64;
        f_flag += fread(&tmp_64, sizeof(tmp_64), 1, fp);
        (*sg)->arc[i].ul = tmp_64;

        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*sg)->arc[i].v = tmp;
    }


    for (i = 0; i < (*sg)->n_seq; i++)
    {
        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*sg)->seq[i].c = tmp;
        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*sg)->seq[i].del = tmp;
        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*sg)->seq[i].len = tmp;
    }
    
    free(index_name);    
    fflush(fp);
    fclose(fp);

    return 1;
}


int write_debug_graph(asg_t *sg, ma_hit_t_alloc* sources, ma_sub_t* coverage_cut, 
char* output_file_name, long long n_read, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex)
{

    char* gfa_name = (char*)malloc(strlen(output_file_name)+55);

    ////write_All_reads(&R_INF, gfa_name);
    sprintf(gfa_name, "%s.all.debug.source", output_file_name);
    write_ma_hit_ts(sources, R_INF.total_reads, gfa_name);
    sprintf(gfa_name, "%s.all.debug.reverse", output_file_name);
    write_ma_hit_ts(reverse_sources, R_INF.total_reads, gfa_name);
    sprintf(gfa_name, "%s.all.debug.coverage_cut", output_file_name);
    write_coverage_cut(coverage_cut, gfa_name, n_read);
    sprintf(gfa_name, "%s.all.debug.ruIndex", output_file_name);
    write_ruIndex(ruIndex, gfa_name);
    sprintf(gfa_name, "%s.all.debug.asg_t", output_file_name);
    write_asg_t(sg, gfa_name);
    free(gfa_name);
    return 1;
}


int load_debug_graph(asg_t** sg, ma_hit_t_alloc** sources, ma_sub_t** coverage_cut, 
char* output_file_name, ma_hit_t_alloc** reverse_sources, R_to_U* ruIndex)
{
    FILE* fp = NULL;
    char* gfa_name = (char*)malloc(strlen(output_file_name)+55);
    sprintf(gfa_name, "%s.all.debug.source.bin", output_file_name);
    fp = fopen(gfa_name, "r"); if(!fp) return 0;
    sprintf(gfa_name, "%s.all.debug.reverse.bin", output_file_name);
    fp = fopen(gfa_name, "r"); if(!fp) return 0;
    sprintf(gfa_name, "%s.all.debug.coverage_cut.bin", output_file_name);
    fp = fopen(gfa_name, "r"); if(!fp) return 0;
    sprintf(gfa_name, "%s.all.debug.ruIndex.bin", output_file_name);
    fp = fopen(gfa_name, "r"); if(!fp) return 0;
    sprintf(gfa_name, "%s.all.debug.asg_t.bin", output_file_name);
    fp = fopen(gfa_name, "r"); if(!fp) return 0;

    if((*sources)!=NULL)
    {
        destory_ma_hit_t_alloc((*sources));
    }

    if((*reverse_sources)!=NULL)
    {
        destory_ma_hit_t_alloc((*reverse_sources));
    }

    if((*coverage_cut)!=NULL)
    {
        free((*coverage_cut));
    }

    if((ruIndex)!=NULL)
    {
        destory_R_to_U((ruIndex));
    }

    if((*sg)!=NULL)
    {
        asg_destroy(*sg);
    }


    

    sprintf(gfa_name, "%s.all.debug.source", output_file_name);
    if(!load_ma_hit_ts(sources, gfa_name))
    {
        return 0;
    }

    sprintf(gfa_name, "%s.all.debug.reverse", output_file_name);
    if(!load_ma_hit_ts(reverse_sources, gfa_name))
    {
        return 0;
    }

    sprintf(gfa_name, "%s.all.debug.coverage_cut", output_file_name);
    if(!load_coverage_cut(coverage_cut, gfa_name))
    {
        return 0;
    }

    sprintf(gfa_name, "%s.all.debug.ruIndex", output_file_name);
    if(!load_ruIndex(ruIndex, gfa_name))
    {
        return 0;
    }

    sprintf(gfa_name, "%s.all.debug.asg_t", output_file_name);
    if(!load_asg_t(sg, gfa_name))
    {
        return 0;
    }

    return 1;
}


void output_contig_graph_primary(asg_t *sg, ma_hit_t_alloc* sources, ma_sub_t* coverage_cut, 
char* output_file_name, long long n_read, long long bubble_dist, long long tipsLen, 
float tip_drop_ratio, 
long long circleLen, long long stops_threshold, float drop_ratio,
ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex)
{
    asg_cut_tip_primary(sg, NULL, tipsLen);
    long long pre_cons = get_graph_statistic(sg);
    long long cur_cons = 0;
    while(pre_cons != cur_cons)
    {
        pre_cons = get_graph_statistic(sg);
        ///need consider tangles
        asg_pop_bubble_primary(sg, bubble_dist);
        asg_arc_del_simple_circle_untig(sources, coverage_cut, sg, circleLen, 1);
        ///need consider tangles
        asg_arc_cut_long_tip_primary(sg, NULL, tip_drop_ratio);
        ///need consider tangles
        asg_arc_cut_long_equal_tips_assembly(sg, reverse_sources, 2, ruIndex);
        asg_arc_cut_long_tip_primary_complex(sg, tip_drop_ratio, stops_threshold);
        asg_arc_cut_long_equal_tips_assembly_complex(sg, reverse_sources, 2, stops_threshold, ruIndex);
        cur_cons = get_graph_statistic(sg);
    }

    
    asg_arc_simple_large_bubbles(sg, reverse_sources, 2, ruIndex);
    asg_arc_identify_simple_bubbles_multi(sg, 1);
    ///we don't need a special function here since it just removes edges instead of nodes
    asg_arc_del_short_false_link_primary(sg, 0.6, 0.85, bubble_dist, reverse_sources, 
    asm_opt.max_short_tip, ruIndex);
    

    
    asg_arc_del_short_diploid_by_length(sg, drop_ratio, asm_opt.max_short_tip, reverse_sources, 
        asm_opt.max_short_tip, stops_threshold, 0, 1, 1, ruIndex);
    asg_cut_tip_primary(sg, NULL, tipsLen);
    ///second round
    pre_cons = get_graph_statistic(sg);
    cur_cons = 0;
    while(pre_cons != cur_cons)
    {
        pre_cons = get_graph_statistic(sg);
        asg_pop_bubble_primary(sg, bubble_dist);
        asg_arc_del_simple_circle_untig(sources, coverage_cut, sg, circleLen, 1);
        asg_arc_cut_long_tip_primary(sg, NULL, tip_drop_ratio);
        asg_arc_cut_long_equal_tips_assembly(sg, reverse_sources, 2, ruIndex);
        asg_arc_cut_long_tip_primary_complex(sg, tip_drop_ratio, stops_threshold);
        asg_arc_cut_long_equal_tips_assembly_complex(sg, reverse_sources, 2, stops_threshold, ruIndex);
        ////here is the difference
        asg_arc_del_short_diploid_by_length(sg, drop_ratio, asm_opt.max_short_tip, reverse_sources, 
        asm_opt.max_short_tip, stops_threshold, 0, 1, 1, ruIndex);
        asg_cut_tip_primary(sg, NULL, tipsLen);
        cur_cons = get_graph_statistic(sg);
    }
    ///unroll_tangles(sg, reverse_sources, 1, 2);
    label_tangles(sg, reverse_sources, 20, 100, 0.05, 0.2, bubble_dist, tipsLen, tip_drop_ratio,
    stops_threshold, ruIndex, 0);


    ma_ug_t *ug = NULL;
    ug = ma_ug_gen_primary(sg, PRIMARY_LABLE);
    adjust_utg(sg, ug);
    ma_ug_seq(ug, &R_INF, coverage_cut, n_read);

    
    fprintf(stderr, "Writing primary contig GFA to disk... \n");
    char* gfa_name = (char*)malloc(strlen(output_file_name)+35);
    sprintf(gfa_name, "%s.p_ctg.gfa", output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    ma_ug_print(ug, &R_INF, coverage_cut, output_file);
    fclose(output_file);

    sprintf(gfa_name, "%s.p_ctg.noseq.gfa", output_file_name);
    output_file = fopen(gfa_name, "w");
    ma_ug_print_simple(ug, &R_INF, coverage_cut, output_file);
    fclose(output_file);

    free(gfa_name);
    ma_ug_destroy(ug);
}



void output_contig_graph_primary_bubble_poping(asg_t *sg, ma_hit_t_alloc* sources, ma_sub_t* coverage_cut, 
char* output_file_name, long long n_read, long long bubble_dist, long long tipsLen, 
float tip_drop_ratio, 
long long circleLen, long long stops_threshold, float drop_ratio,
ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex)
{
    asg_cut_tip_primary(sg, NULL, tipsLen);
    long long pre_cons = get_graph_statistic(sg);
    long long cur_cons = 0;
    while(pre_cons != cur_cons)
    {
        pre_cons = get_graph_statistic(sg);
        ///need consider tangles
        asg_pop_bubble_primary(sg, bubble_dist);
        asg_arc_del_simple_circle_untig(sources, coverage_cut, sg, circleLen, 1);
        cur_cons = get_graph_statistic(sg);
    }

    
    asg_arc_simple_large_bubbles(sg, reverse_sources, 2, ruIndex);
    asg_arc_identify_simple_bubbles_multi(sg, 1);
    ///we don't need a special function here since it just removes edges instead of nodes
    asg_arc_del_short_false_link_primary(sg, 0.6, 0.85, bubble_dist, reverse_sources, 
    asm_opt.max_short_tip, ruIndex);
    

    
    asg_arc_del_short_diploid_by_length(sg, drop_ratio, asm_opt.max_short_tip, reverse_sources, 
        asm_opt.max_short_tip, stops_threshold, 0, 1, 1, ruIndex);
    asg_cut_tip_primary(sg, NULL, tipsLen);
    ///second round
    pre_cons = get_graph_statistic(sg);
    cur_cons = 0;
    while(pre_cons != cur_cons)
    {
        pre_cons = get_graph_statistic(sg);
        asg_pop_bubble_primary(sg, bubble_dist);
        asg_arc_del_simple_circle_untig(sources, coverage_cut, sg, circleLen, 1);
        ////here is the difference
        asg_arc_del_short_diploid_by_length(sg, drop_ratio, asm_opt.max_short_tip, reverse_sources, 
        asm_opt.max_short_tip, stops_threshold, 0, 1, 1, ruIndex);
        asg_cut_tip_primary(sg, NULL, tipsLen);
        cur_cons = get_graph_statistic(sg);
    }
    label_tangles(sg, reverse_sources, 20, 100, 0.05, 0.2, bubble_dist, tipsLen, tip_drop_ratio,
    stops_threshold, ruIndex, 1);
    
    
    ma_ug_t *ug = NULL;
    ug = ma_ug_gen_primary(sg, PRIMARY_LABLE);
    ma_ug_seq(ug, &R_INF, coverage_cut, n_read);

    
    fprintf(stderr, "Writing primary contig GFA to disk... \n");
    char* gfa_name = (char*)malloc(strlen(output_file_name)+35);
    sprintf(gfa_name, "%s.p_ctg.gfa", output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    ma_ug_print(ug, &R_INF, coverage_cut, output_file);
    fclose(output_file);

    sprintf(gfa_name, "%s.p_ctg.noseq.gfa", output_file_name);
    output_file = fopen(gfa_name, "w");
    ma_ug_print_simple(ug, &R_INF, coverage_cut, output_file);
    fclose(output_file);

    free(gfa_name);
    ma_ug_destroy(ug);
}


void output_contig_graph_alternative(asg_t *sg, ma_sub_t* coverage_cut, char* output_file_name, long long n_read)
{
    ma_ug_t *ug = NULL;
    ug = ma_ug_gen_primary(sg, ALTER_LABLE);
    ma_ug_seq(ug, &R_INF, coverage_cut, n_read);

    fprintf(stderr, "Writing alternate contig GFA to disk... \n");
    char* gfa_name = (char*)malloc(strlen(output_file_name)+35);
    sprintf(gfa_name, "%s.a_ctg.gfa", output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    ma_ug_print(ug, &R_INF, coverage_cut, output_file);
    fclose(output_file);
    
    sprintf(gfa_name, "%s.a_ctg.noseq.gfa", output_file_name);
    output_file = fopen(gfa_name, "w");
    ma_ug_print_simple(ug, &R_INF, coverage_cut, output_file);
    fclose(output_file);

    free(gfa_name);
    ma_ug_destroy(ug);
}

int output_tips(asg_t *g, const All_reads *RNF)
{
	uint32_t v, n_vtx = g->n_seq * 2;
    for (v = 0; v < n_vtx; ++v)
    {
        if (g->seq[v>>1].del) continue;

        if(asg_arc_n(g, v) == 0)
        {
            fprintf(stderr, "%.*s\n",
             (int)Get_NAME_LENGTH((*RNF), v>>1), 
             Get_NAME((*RNF), v>>1));
        }
    }

    return 1;
}

void collect_abnormal_edges(ma_hit_t_alloc* paf, ma_hit_t_alloc* rev_paf, long long readNum)
{
    double startTime = Get_T();
    long long T_edges, T_Single_Dir_Edges_0, T_Single_Dir_Edges_1, T_Conflict_Equal_Edges, T_Conflict_Strong_Edges;
    T_edges = T_Single_Dir_Edges_0 = T_Single_Dir_Edges_1 = T_Conflict_Equal_Edges = T_Conflict_Strong_Edges = 0;
    long long T_Single_Dir_Edges_1_1000 = 0;
    long long related_reads = 0;
    long long related_overlaps = 0;
    long long i, j;
    uint32_t qn, tn;
    int is_equal_f, is_strong_f; 
    int is_equal_b, is_strong_b, is_exist_b; 

    kvec_t(uint64_t) edge_vector;
    kv_init(edge_vector);
    
    for (i = 0; i < readNum; i++)
    {
        for (j = 0; j < (long long)paf[i].length; j++)
        {
            qn = Get_qn(paf[i].buffer[j]);
            tn = Get_tn(paf[i].buffer[j]);
            T_edges++;

            is_equal_f = paf[i].buffer[j].el;
            is_strong_f = paf[i].buffer[j].ml;
            
            is_exist_b = get_specific_overlap(&(paf[tn]), tn, qn);
            if(is_exist_b == -1)
            {
                is_exist_b = get_specific_overlap(&(rev_paf[tn]), tn, qn);
                if(is_exist_b != -1)
                {
                    T_Single_Dir_Edges_0++;
                    kv_push(uint64_t, edge_vector, qn);
                    kv_push(uint64_t, edge_vector, tn);
                    ///related_overlaps += paf[qn].length + rev_paf[qn].length + paf[tn].length + rev_paf[tn].length;
                    // fprintf(stderr, "%.*s(%d) ---(+)--> %.*s(%d), Len: %d\n", 
                    // Get_NAME_LENGTH(R_INF, qn), Get_NAME(R_INF, qn), qn,
                    // Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn), tn,
                    // Get_qe(paf[i].buffer[j]) - Get_qs(paf[i].buffer[j]));

                    // fprintf(stderr, "%.*s(%d) ---(-)--> %.*s(%d), Len: %d\n\n", 
                    // Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn), tn,
                    // Get_NAME_LENGTH(R_INF, qn), Get_NAME(R_INF, qn), qn,
                    // Get_qe(rev_paf[tn].buffer[is_exist_b]) - Get_qs(rev_paf[tn].buffer[is_exist_b]));
                } 
                else
                {
                    T_Single_Dir_Edges_1++;
                    if(Get_qe(paf[i].buffer[j]) - Get_qs(paf[i].buffer[j]) >= 1000)
                    {
                        T_Single_Dir_Edges_1_1000++;
                        // fprintf(stderr, "%.*s(%d) ---(%d)--> %.*s(%d), Len: %d\n\n", 
                        // Get_NAME_LENGTH(R_INF, qn), Get_NAME(R_INF, qn), qn,
                        // is_strong_f,
                        // Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn), tn,
                        // Get_qe(paf[i].buffer[j]) - Get_qs(paf[i].buffer[j]));
                    }
                }

                related_reads = related_reads + 2;
            }
            else
            {
                is_equal_b = paf[tn].buffer[is_exist_b].el;
                is_strong_b = paf[tn].buffer[is_exist_b].ml;

                if(is_equal_f != is_equal_b)
                {
                    T_Conflict_Equal_Edges++;

                    // fprintf(stderr, "%.*s(%d) ---(%d)--> %.*s(%d), Len: %d\n", 
                    // Get_NAME_LENGTH(R_INF, qn), Get_NAME(R_INF, qn), qn,
                    // is_equal_f,
                    // Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn), tn,
                    // Get_qe(paf[i].buffer[j]) - Get_qs(paf[i].buffer[j]));

                    // fprintf(stderr, "%.*s(%d) ---(%d)--> %.*s(%d), Len: %d\n\n", 
                    // Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn), tn,
                    // is_equal_b,
                    // Get_NAME_LENGTH(R_INF, qn), Get_NAME(R_INF, qn), qn,
                    // Get_qe(paf[tn].buffer[is_exist_b]) - Get_qs(paf[tn].buffer[is_exist_b]));
                } 
                
                if(is_strong_f != is_strong_b)
                {
                    T_Conflict_Strong_Edges++;
                    // fprintf(stderr, "%.*s(%d) ---(%d)--> %.*s(%d), Len: %d\n", 
                    // Get_NAME_LENGTH(R_INF, qn), Get_NAME(R_INF, qn), qn,
                    // is_strong_f,
                    // Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn), tn,
                    // Get_qe(paf[i].buffer[j]) - Get_qs(paf[i].buffer[j]));

                    // fprintf(stderr, "%.*s(%d) ---(%d)--> %.*s(%d), Len: %d\n\n", 
                    // Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn), tn,
                    // is_strong_b,
                    // Get_NAME_LENGTH(R_INF, qn), Get_NAME(R_INF, qn), qn,
                    // Get_qe(paf[tn].buffer[is_exist_b]) - Get_qs(paf[tn].buffer[is_exist_b]));
                } 

                if(is_equal_f != is_equal_b || is_strong_f != is_strong_b)
                {
                    related_reads++;
                }
            }
        }
    }

    radix_sort_arch64(edge_vector.a, edge_vector.a + edge_vector.n);
    uint64_t pre = (uint64_t)-1;
    long long mn = 0;
    for (i = 0; i < (long long)edge_vector.n; i++)
    {
        if(pre != edge_vector.a[i])
        {
            mn++;
            pre = edge_vector.a[i];
            related_overlaps += paf[pre].length + rev_paf[pre].length;
        }

        if(i>0 && edge_vector.a[i] < edge_vector.a[i-1]) fprintf(stderr, "hehe\n");
    }
    
    

    fprintf(stderr, "****************statistic for abnormal overlaps****************\n");
    fprintf(stderr, "overlaps #: %lld\n", T_edges);
	fprintf(stderr, "one direction overlaps (different phasing)#: %lld\n", T_Single_Dir_Edges_0);
    fprintf(stderr, "one direction overlaps (missing)#: %lld\n", T_Single_Dir_Edges_1);
    fprintf(stderr, "one direction overlaps (missing) >= 1000#: %lld\n", T_Single_Dir_Edges_1_1000);
    fprintf(stderr, "conflict strong/weak overlaps #: %lld\n", T_Conflict_Strong_Edges);
    fprintf(stderr, "conflict exact/inexact overlaps #: %lld\n", T_Conflict_Equal_Edges);
    fprintf(stderr, "related_reads #: %lld/%lld\n", related_reads, mn);
    fprintf(stderr, "related_overlaps #: %lld\n", related_overlaps);
    fprintf(stderr, "****************statistic for abnormal overlaps****************\n");

    fprintf(stderr, "[M::%s] took %0.2fs\n\n", __func__, Get_T()-startTime);

    kv_destroy(edge_vector);
}

///if we don't have this function, we just simply remove all one-direction edges
///by utilizing this function, some one-direction edges can be recovered as two-direction edges  
void try_rescue_overlaps(ma_hit_t_alloc* paf, ma_hit_t_alloc* rev_paf, long long readNum,
long long rescue_threshold)
{
    double startTime = Get_T();
    long long i, j, revises = 0;
    uint32_t qn, tn, qs, qe;

    kvec_t(uint64_t) edge_vector;
    kv_init(edge_vector);
    kvec_t(uint64_t) edge_vector_index;
    kv_init(edge_vector_index);
    kvec_t(uint32_t) b;
    kv_init(b);
    uint64_t flag;
    int index;
    
    for (i = 0; i < readNum; i++)
    {
        edge_vector.n = 0;
        edge_vector_index.n = 0;
        for (j = 0; j < (long long)rev_paf[i].length; j++)
        {
            qn = Get_qn(rev_paf[i].buffer[j]);
            tn = Get_tn(rev_paf[i].buffer[j]);
            index = get_specific_overlap(&(paf[tn]), tn, qn);
            if(index != -1)
            {
                flag = tn;
                flag = flag << 32;
                flag = flag | (uint64_t)(index);
                kv_push(uint64_t, edge_vector, flag);
                kv_push(uint64_t, edge_vector_index, j);
            }
        }

        ///based on qn, all edges at edge_vector/edge_vector_index come from different haplotype
        ///but at another direction, all these edges come from the same haplotype
        //here we want to recover these edges
        if((long long)edge_vector_index.n >= rescue_threshold)
        {
            kv_resize(uint32_t, b, edge_vector_index.n);
            b.n = 0;
            for (j = 0; j < (long long)edge_vector_index.n; j++)
            {
                qs = Get_qs(rev_paf[i].buffer[edge_vector_index.a[j]]);
                qe = Get_qe(rev_paf[i].buffer[edge_vector_index.a[j]]);
                kv_push(uint32_t, b, qs<<1);
                kv_push(uint32_t, b, qe<<1|1);
            }

            ks_introsort_uint32_t(b.n, b.a);
            int dp = 0, start = 0, max_dp = 0;
            ma_sub_t max_interval;
            max_interval.s = max_interval.e = 0;
            for (j = 0, dp = 0; j < (long long)b.n; ++j) 
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
                there are two cases: 
                    1. old_dp = dp + 1 (b.a[j] is qe); 2. old_dp = dp - 1 (b.a[j] is qs);
                **/
                ///if (old_dp < min_dp && dp >= min_dp) ///old_dp < dp, b.a[j] is qs
                if(old_dp < dp) ///b.a[j] is qs
                { 
                    ///case 2, a[j] is qs
                    //here should use dp >= max_dp, instead of dp > max_dp
                    if(dp >= max_dp)
                    {
                        start = b.a[j]>>1;
                        max_dp = dp;
                    }
                } 
                ///else if (old_dp >= min_dp && dp < min_dp) ///old_dp > min_dp, b.a[j] is qe
                else if (old_dp > dp) ///old_dp > min_dp, b.a[j] is qe
                {
                    if(old_dp == max_dp)
                    {
                        max_interval.s = start;
                        max_interval.e = b.a[j]>>1;
                    }
                }
                // else
                // {
                //     fprintf(stderr, "error\n");
                // }
                
            }

            if(max_dp>= rescue_threshold)
            {
                long long m = 0;
                for (j = 0; j < (long long)edge_vector_index.n; j++)
                {
                    qs = Get_qs(rev_paf[i].buffer[edge_vector_index.a[j]]);
                    qe = Get_qe(rev_paf[i].buffer[edge_vector_index.a[j]]);
                    if(qs <= max_interval.s && qe >= max_interval.e)
                    {
                        edge_vector_index.a[m] = edge_vector_index.a[j];
                        edge_vector.a[m] = edge_vector.a[j];
                        m++;
                    }
                }
                edge_vector_index.n = m;
                edge_vector.n = m;

                ///the read itself do not have these overlaps, but all related reads have
                ///we need to remove all overlaps from rev_paf[i], and then add all overlaps to paf[i]
                // fprintf(stderr,"\nadd following %d edges...\n", edge_vector.n);
                // print_revise_edges(&(rev_paf[i]), edge_vector_index.a, edge_vector_index.n);
                remove_overlaps(&(rev_paf[i]), edge_vector_index.a, edge_vector_index.n);
                add_overlaps_from_different_sources(paf, &(paf[i]), edge_vector.a, edge_vector.n);
                revises = revises + edge_vector.n;
            }
        }
    }


    kv_destroy(edge_vector);
    kv_destroy(edge_vector_index);
    kv_destroy(b);
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] took %0.2fs, rescue edges #: %lld\n\n", 
                                    __func__, Get_T()-startTime, revises);
    }
    
}


long long get_coverage(ma_hit_t_alloc* sources, ma_sub_t* coverage_cut, uint64_t n_read)
{
    uint64_t i, j;
    ma_hit_t *h;
    long long R_bases = 0, C_bases = 0;
    for (i = 0; i < n_read; ++i) 
    {
        R_bases += coverage_cut[i].e - coverage_cut[i].s;
        for (j = 0; j < (uint64_t)(sources[i].length); j++)
        {
            h = &(sources[i].buffer[j]);
            C_bases += Get_qe((*h)) - Get_qs((*h));
        }
    }

    return C_bases/R_bases;
}


void pre_clean(ma_hit_t_alloc* sources, ma_sub_t* coverage_cut, asg_t *sg, long long bubble_dist)
{
    while(1)
    {
        int tri_flag = 0;
        tri_flag += asg_arc_del_simple_circle_untig(sources, coverage_cut, sg, 100, 0);
        ///asg_arc_del_single_node_bubble(sg, bubble_dist);
        tri_flag += asg_arc_del_single_node_directly(sg, asm_opt.max_short_tip, sources);
        tri_flag += asg_arc_del_triangular_advance(sg, bubble_dist);
        tri_flag += asg_arc_del_cross_bubble(sg, bubble_dist);
        ///asg_arc_del_single_node_bubble(sg, bubble_dist);
        tri_flag += asg_arc_del_single_node_directly(sg, asm_opt.max_short_tip, sources);
        if(tri_flag == 0)
        {
            break;
        }
    }
}


void init_R_to_U(R_to_U* x, uint64_t len)
{
	x->len = len;
	x->index = (uint32_t*)malloc(sizeof(uint32_t)*(x->len));
	memset(x->index, -1, sizeof(uint32_t)*(x->len));
}

void destory_R_to_U(R_to_U* x)
{
	free(x->index);
}

void set_R_to_U(R_to_U* x, uint32_t rID, uint32_t uID, uint32_t is_Unitig)
{
    // if(rID >= x->len)
	// {
    //     fprintf(stderr, "s: rID: %u, len: %lu\n", rID, x->len);
    //     fflush(stderr);
    // }
	if(rID >= x->len)
	{
		x->index = (uint32_t*)realloc(x->index, (rID + 1)*sizeof(uint32_t));
        memset(x->index + x->len, -1, sizeof(uint32_t)*((rID + 1) - x->len));
        x->len = rID + 1;
	}

	x->index[rID] = uID & (uint32_t)(0x7fffffff);
	x->index[rID] = x->index[rID] | (uint32_t)(is_Unitig<<31);
}


void get_R_to_U(R_to_U* x, uint32_t rID, uint32_t* uID, uint32_t* is_Unitig)
{
    // if(rID >= x->len)
	// {
    //     fprintf(stderr, "r: rID: %u, len: %lu\n", rID, x->len);
    //     fflush(stderr);
    // }
	if(rID >= x->len || (x->index[rID] == (uint32_t)(-1)))
	{
		(*uID) = (uint32_t)-1;
		(*is_Unitig) = (uint32_t)-1;
		return;
	}
    
    (*uID) = x->index[rID] & (uint32_t)(0x7fffffff);
	(*is_Unitig) = (x->index[rID]>>31);
}

void transfor_R_to_U(R_to_U* x)
{
    uint64_t i = 0;
    uint32_t rID, uID, is_Unitig;
    for (i = 0; i < x->len; i++)
    {
        rID = i;
        get_R_to_U(x, rID, &uID, &is_Unitig);
        if(uID == (uint32_t)-1) continue;
        if(is_Unitig == 1) continue;


        ///here i/rID is contained in uID
        rID = uID;
        while (1)
        {
            get_R_to_U(x, rID, &uID, &is_Unitig);
            if(uID == (uint32_t)-1) break;
            if(is_Unitig == 1) break;
            rID = uID;
        }

        set_R_to_U(x, i, rID, 0);
    }


    /**
    for (i = 0; i < x->len; i++)
    {
        rID = i;
        get_R_to_U(x, rID, &uID, &is_Unitig);
        if(uID == (uint32_t)-1) continue;
        if(is_Unitig == 1) continue;
        ///here i/rID is contained in uID
        rID = uID;
        get_R_to_U(x, rID, &uID, &is_Unitig);
        if(uID != (uint32_t)-1)
        {
            fprintf(stderr, "ERROR\n");
        }
    }
    **/
    
}

void build_string_graph_without_clean(
int min_dp, ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, 
long long n_read, uint64_t* readLen, long long mini_overlap_length, 
long long max_hang_length, long long clean_round, long long gap_fuzz,
float min_ovlp_drop_ratio, float max_ovlp_drop_ratio, char* output_file_name, 
long long bubble_dist, int read_graph, int write)
{
    R_to_U ruIndex;
    asg_t *sg = NULL;
    ma_sub_t* coverage_cut = NULL;
    
    ///actually min_thres = asm_opt.max_short_tip + 1 there are asm_opt.max_short_tip reads
    min_thres = asm_opt.max_short_tip + 1;

    // if(load_debug_graph(&sg, &sources, &coverage_cut, output_file_name, &reverse_sources, &ruIndex))
    // {
    //     fprintf(stderr, "debug gfa has been loaded\n");
    //     goto debug_gfa;
    // }

    if (asm_opt.write_index_to_disk && write)
    {
        write_all_data_to_disk(sources, reverse_sources, 
        &R_INF, output_file_name);
    }
    
    init_R_to_U(&ruIndex, n_read);
    
    try_rescue_overlaps(sources, reverse_sources, n_read, 4);        
    
    normalize_ma_hit_t_single_side(sources, n_read);
    clean_weak_ma_hit_t(sources, reverse_sources, n_read);

    // debug_info_of_specfic_read("m64011_190329_072846/80545633/ccs", 
    // sources, reverse_sources, -1, "clean");

    ma_hit_sub(min_dp, sources, n_read, readLen, mini_overlap_length, &coverage_cut);
    detect_chimeric_reads(sources, reverse_sources, n_read, readLen, coverage_cut, 
    FINAL_OVERLAP_ERROR_RATE*2);
    ma_hit_cut(sources, n_read, readLen, mini_overlap_length, &coverage_cut);
    ///it seems we do not need ma_hit_flt
    ma_hit_flt(sources, n_read, coverage_cut, max_hang_length, mini_overlap_length);
    ma_hit_contained(sources, n_read, coverage_cut, &ruIndex, max_hang_length, mini_overlap_length);
    
    sg = ma_sg_gen(sources, n_read, coverage_cut, max_hang_length, mini_overlap_length);
    asg_arc_del_trans(sg, gap_fuzz);
    asm_opt.coverage = get_coverage(sources, coverage_cut, n_read);


    if(VERBOSE >= 1)
    {
        char* unlean_name = (char*)malloc(strlen(output_file_name)+25);
        sprintf(unlean_name, "%s.unclean", output_file_name);
        output_read_graph(sg, coverage_cut, unlean_name, n_read);
        free(unlean_name);
    }


    
    asg_cut_tip(sg, asm_opt.max_short_tip);
    
    // debug_info_of_specfic_node("m64016_190918_162737/72220752/ccs", sg, "cut_tip");
    ///asg_arc_del_short_diploid_unclean(sg, corase_ovlp_drop_ratio, sources, reverse_sources);

    // asg_arc_del_single_node_bubble(sg, bubble_dist);
    // asg_cut_tip(sg, asm_opt.max_short_tip);
    ///asg_cut_tip(sg, asm_opt.max_short_tip);
    
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

            if(VERBOSE >= 1)
            {
                fprintf(stderr, "\n\n**********%d-th round drop: drop_ratio = %f**********\n", 
                i, drop_ratio);
            }
            

            pre_clean(sources, coverage_cut, sg, bubble_dist);


            ///asg_arc_del_orthology(sg, reverse_sources, drop_ratio, asm_opt.max_short_tip);
            // asg_arc_del_orthology_multiple_way(sg, reverse_sources, drop_ratio, asm_opt.max_short_tip);
            // asg_cut_tip(sg, asm_opt.max_short_tip);
            
            /****************************may have bugs********************************/
            asg_arc_identify_simple_bubbles_multi(sg, 1);
            //reomve edge between two chromesomes
            asg_arc_del_false_node(sg, asm_opt.max_short_tip);
            asg_cut_tip(sg, asm_opt.max_short_tip);
            /****************************may have bugs********************************/

            /****************************may have bugs********************************/
            ///asg_arc_identify_simple_bubbles_multi(sg, 1);
            asg_arc_identify_simple_bubbles_multi(sg, 0);
            ///asg_arc_del_short_diploid_unclean_exact(sg, drop_ratio, sources);
            asg_arc_del_short_diploid_by_exact(sg, asm_opt.max_short_tip, sources);
            asg_cut_tip(sg, asm_opt.max_short_tip);
            /****************************may have bugs********************************/


            asg_arc_identify_simple_bubbles_multi(sg, 1);
            asg_arc_del_short_diploid_by_length(sg, drop_ratio, asm_opt.max_short_tip, reverse_sources, 
            asm_opt.max_short_tip, 1, 1, 0, 0, &ruIndex);
            asg_cut_tip(sg, asm_opt.max_short_tip);

            asg_arc_identify_simple_bubbles_multi(sg, 1);
            asg_arc_del_short_false_link(sg, 0.6, 0.85, bubble_dist, reverse_sources, 
            asm_opt.max_short_tip, &ruIndex);

            asg_arc_identify_simple_bubbles_multi(sg, 1);
            asg_arc_del_complex_false_link(sg, 0.6, 0.85, bubble_dist, reverse_sources, asm_opt.max_short_tip);

            asg_cut_tip(sg, asm_opt.max_short_tip);
        }
    }

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "\n\n**********final clean**********\n");
    }


    pre_clean(sources, coverage_cut, sg, bubble_dist);


    asg_arc_del_short_diploi_by_suspect_edge(sg, asm_opt.max_short_tip);
    asg_cut_tip(sg, asm_opt.max_short_tip);
    asg_arc_del_triangular_directly(sg, asm_opt.max_short_tip, reverse_sources, &ruIndex);


    asg_arc_identify_simple_bubbles_multi(sg, 0);
    asg_arc_del_orthology_multiple_way(sg, reverse_sources, 0.4, asm_opt.max_short_tip, &ruIndex);
    asg_cut_tip(sg, asm_opt.max_short_tip);

    



    asg_arc_identify_simple_bubbles_multi(sg, 0);
    asg_arc_del_too_short_overlaps(sg, 2000, min_ovlp_drop_ratio, reverse_sources, 
    asm_opt.max_short_tip, &ruIndex);
    asg_cut_tip(sg, asm_opt.max_short_tip);

    asg_arc_del_simple_circle_untig(sources, coverage_cut, sg, 100, 0);

    /**
    asg_arc_identify_simple_bubbles_multi(sg, 1);
    asg_arc_del_short_false_link_advance(sg, 0.6, 0.85, bubble_dist, reverse_sources, asm_opt.max_short_tip);
    **/
    
    




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

        if(c_tips) asg_cut_tip(sg, asm_opt.max_short_tip);
        i++;
    }
    **/
    
    /****************************may have bugs********************************/

    /**
    memset(sg->seq_vis, 0, sg->n_seq*2*sizeof(uint8_t));
    asg_arc_del_short_diploid_by_exact(sg, asm_opt.max_short_tip, sources);
    asg_cut_tip(sg, asm_opt.max_short_tip);
    **/
    
    
    ///out:
    ///output_tips(sg, &R_INF);

    /***********************debug************************/
    // asg_t* debug_g = NULL;
    // debug_g = copy_graph(sg, 10);
    // asg_destroy(sg);
    // sg = debug_g;
    /***********************debug************************/






    output_unitig_graph(sg, coverage_cut, output_file_name, n_read);

    if(VERBOSE >= 1)
    {
        output_read_graph(sg, coverage_cut, output_file_name, n_read);
    }

    output_unitig_graph_without_small_bubbles_primary(sg, coverage_cut, output_file_name, n_read, 
    asm_opt.small_pop_bubble_size, asm_opt.max_short_tip);


    // write_debug_graph(sg, sources, coverage_cut, output_file_name, n_read, reverse_sources, &ruIndex);
    // debug_gfa:


    output_contig_graph_primary(sg, sources, coverage_cut, output_file_name, n_read, bubble_dist, 
    (asm_opt.max_short_tip*2), 0.15, 20, 3, 0.9, reverse_sources, &ruIndex);
    // output_contig_graph_primary_bubble_poping(sg, sources, coverage_cut, output_file_name, n_read, bubble_dist, 
    // (asm_opt.max_short_tip*2), 0.15, 20, 3, 0.9, reverse_sources, &ruIndex);

    output_contig_graph_alternative(sg, coverage_cut, output_file_name, n_read);

    asg_destroy(sg);
    free(coverage_cut);
    destory_R_to_U(&ruIndex);
}
