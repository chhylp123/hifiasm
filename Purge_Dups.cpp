#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <stdlib.h>
#include "ksort.h"
#include "Purge_Dups.h"
#include "Overlaps.h"
#include "Correct.h"
#include "kthread.h"
#include "kdq.h"
#include "hic.h"
#include "rcut.h"
#include "tovlp.h"

KDQ_INIT(uint64_t)
KSORT_INIT_GENERIC(uint64_t)

uint8_t debug_enable = 0;

typedef struct {
    uint64_t weight;
    uint32_t x_beg_pos;
    uint32_t x_end_pos;
    uint32_t y_beg_pos;
    uint32_t y_end_pos;
    uint32_t index_beg;
    uint32_t index_end;
    long long score;
    uint8_t rev;
    asg_arc_t t;
}hap_candidates;

typedef struct {
    kvec_t(hap_candidates) a;
    uint64_t i;
}kvec_hap_candidates;


typedef struct {
    uint64_t* vote_counting;
    uint8_t* visit;
    kvec_t_u64_warp u_vecs;
    kvec_asg_arc_t_offset u_buffer;
    kvec_t_i32_warp u_buffer_tailIndex;
    kvec_t_i32_warp u_buffer_prevIndex;
    kvec_t_u8_warp u_buffer_flag;
    kvec_t_i32_warp u_buffer_beg;
    kvec_hap_candidates u_can;
}hap_alignment_struct;


typedef struct {
    hap_alignment_struct* buf;
    uint32_t num_threads;
    uint8_t *hh;
    ma_ug_t *ug;
    asg_t *read_g;
    ma_hit_t_alloc* sources;
    ma_hit_t_alloc* reverse_sources;
    R_to_U* ruIndex;
    ma_sub_t *coverage_cut; 
    uint64_t* position_index;
    float Hap_rate;
    int max_hang;
    int min_ovlp;
    float chain_rate;
    hap_overlaps_list* all_ovlp;
    long long cov_threshold;
    hap_cov_t *cov;
}hap_alignment_struct_pip;

typedef struct {
    uint32_t baseBeg, baseEnd;
    uint32_t nodeBeg, nodeEnd;
    uint32_t h_lev_idx;
    uint32_t h_status, c_ug_id;
}p_node_t;

typedef struct {
    uint32_t beg;
    uint32_t occ;
}p_g_in_t;

typedef struct {
    ma_ug_t *ug;
    kvec_t(p_node_t) pg_het_node;
    asg_t *pg_het;
    asg_t *pg_h_lev;
    kvec_t(p_g_in_t) pg_h_lev_idx;
}p_g_t;

void print_peak_line(int c, int x, int exceed, int64_t cnt)
{
	int j;
	if (c >= 0) fprintf(stderr, "[M::%s] %5d: ", __func__, c);
	else fprintf(stderr, "[M::%s] %5s: ", __func__, "rest");
	for (j = 0; j < x; ++j) fputc('*', stderr);
	if (exceed) fputc('>', stderr);
	fprintf(stderr, " %lld\n", (long long)cnt);
}

void print_peak(long long* cov_buf, long long cov_buf_length, long long max_i)
{
    long long i;
    const long long hist_max = 100;
    // print histogram
	for (i = 0; i < cov_buf_length; ++i) 
    {
		long long x, exceed = 0;
		x = (int)((double)hist_max * cov_buf[i] / cov_buf[max_i] + .499);
		if (x > hist_max) exceed = 1, x = hist_max; // may happen if cnt[2] is higher
		if (i > max_i && x == 0) break;
		print_peak_line(i, x, exceed, cov_buf[i]);
	}
	{
		long long x, exceed = 0;
		long long rest = 0;
		for (; i < cov_buf_length; ++i) rest += cov_buf[i];
		x = (int)((double)hist_max * rest / cov_buf[max_i] + .499);
		if (x > hist_max) exceed = 1, x = hist_max;
		print_peak_line(-1, x, exceed, rest);
	}
}


void get_read_peak(asg_t *read_g, long long* cov_buf, long long cov_buf_length, long long* topo_peak_cov, 
long long* hom_peak, long long* het_peak, long long* k_mer_only, long long* coverage_only, long long g_size)
{
    long long i, start, err_i, max_i, max2_i, max3_i, topo_peak_i, max, max2, max3, topo_peak, min;

    i = start = err_i = max_i = max2_i = max3_i = topo_peak_i = -1;
    max = max2 = max3 = topo_peak = min = -1;

    ///cov_buf[0] is usually very large
    for (i = 1; i < cov_buf_length; ++i)
    {
        if(cov_buf[i] > cov_buf[i-1]) break;
    }
    err_i = i - 1;
    // find the global highest peak
	max_i = err_i + 1, max = cov_buf[max_i];
    for (i = max_i; i < cov_buf_length; ++i)
    {
        if (cov_buf[i] > max)
        {
            max = cov_buf[i];
            max_i = i;
        }
    }

    ///print_peak(cov_buf, cov_buf_length, max_i);

    // look for smaller peak on the low end
	max2 = -1; max2_i = -1;
	for (i = max_i - 1; i > err_i; --i) 
    {
		///at first, it should be a peak
		if (cov_buf[i] >= cov_buf[i-1] && cov_buf[i] >= cov_buf[i+1]) 
        {
			if (cov_buf[i] > max2)
            {
                max2 = cov_buf[i];
                max2_i = i;
            } 
		}
	}

    ///fprintf(stderr, "***max2: %lld, max2_i: %lld\n", max2, max2_i);

    if (max2_i != -1 && max2_i > err_i && max2_i < max_i) 
    {
		for (i = max2_i + 1, min = max; i < max_i; ++i)
        {
            if (cov_buf[i] < min) min = cov_buf[i];
        }
			
		///if the second peak is not significant
		if(max2 < max * 0.05 || min > max2 * 0.95) max2 = max2_i = -1;
	}
    if(max2 < max*0.0075) max2 = max2_i = -1;
		

    // look for smaller peak on the high end
	max3 = -1; max3_i = -1;
    // we'd better use i < cov_buf_length - 1, since cov_buf[cov_buf_length-1] may have problem
	for (i = max_i + 1; i < cov_buf_length - 1; ++i) 
    {
		//at first, it should be a peak
		if (cov_buf[i] >= cov_buf[i-1] && cov_buf[i] >= cov_buf[i+1]) 
        {
			if (cov_buf[i] > max3)
            {
                max3 = cov_buf[i], max3_i = i;
            } 
		}
	}

    ///fprintf(stderr, "***max3: %lld, max3_i: %lld\n", max3, max3_i);
    
    //if found a peak
	if (max3 != -1 && max3_i > max_i) 
    {
		for (i = max_i + 1, min = max; i < max3_i; ++i)
        {
            if (cov_buf[i] < min) min = cov_buf[i];
        }
			
		if (max3 < max * 0.05 || min > max3 * 0.95 || max3_i > max_i * 3) max3 = max3_i = -1;
	}
    if (max3 < max*0.0075) max3 = max3_i = -1;

    

    (*hom_peak) = (*het_peak) = -1;
    if (topo_peak_cov && (*topo_peak_cov) < cov_buf_length)
    {
        topo_peak_i = (*topo_peak_cov);
        topo_peak = cov_buf[topo_peak_i];
        if (topo_peak <= max * 0.05) topo_peak_i = topo_peak = -1;
    }

    if(asm_opt.purge_level_primary == 0) 
    {
        (*hom_peak) = max_i;
        return;
    }
    
    
    long long k_mer_het, k_mer_hom, coverage_het, coverage_hom, alter_peak;
    k_mer_het = k_mer_hom = coverage_het = coverage_hom = alter_peak = -1;

    alter_peak = topo_peak_i;
    k_mer_het = asm_opt.het_cov;
    k_mer_hom = asm_opt.hom_cov;
    if(max3_i > 0)
    {
        coverage_het = max_i;
        coverage_hom = max3_i;
    }
    else
    {
        coverage_het = max2_i;
        coverage_hom = max_i;
    }

    if(g_size > 0) {
        long long n_bs, m_peak_hom = -1;
        int p_ht = -1;
        for (i = n_bs = 0; i < read_g->n_seq; i++) n_bs += read_g->seq[i].len;
        m_peak_hom = n_bs/g_size;
        if(m_peak_hom > 0) {
            p_ht = -1;
            coverage_hom = adj_m_peak_hom(m_peak_hom, max_i, max2_i, max3_i, &p_ht);
            coverage_het = p_ht;
        }   
    }

    if(k_mer_het != -1)
    {
        (*het_peak) = k_mer_het;
        (*hom_peak) = k_mer_hom;
        return;
    }
    else if(coverage_het != -1)
    {
        (*het_peak) = coverage_het;
        (*hom_peak) = coverage_hom;
        return;
    }
    else if(k_mer_hom > coverage_hom*1.5)
    {
        (*het_peak) = coverage_hom;
        (*hom_peak) = k_mer_hom;
        return;
    }
    else if(alter_peak != -1)
    {
        ///if peak is het, coverage peak is more reliable 
        if(coverage_hom >= alter_peak*0.8 && coverage_hom <= alter_peak*1.2) 
        {
            (*het_peak) = coverage_hom;
            return;
        }///if peak is homo, k-mer peak is more reliable 
        else if(k_mer_hom >= alter_peak*0.8*2 && k_mer_hom <= alter_peak*1.2*2)
        {
            (*hom_peak) = k_mer_hom;
            return;
        }
    }

    (*k_mer_only) = k_mer_hom;
    (*coverage_only) = coverage_hom;
    

    // fprintf(stderr, "max: %lld, max_i: %lld\n", max, max_i);
    // fprintf(stderr, "max2: %lld, max2_i: %lld\n", max2, max2_i);
    // fprintf(stderr, "max3: %lld, max3_i: %lld\n", max3, max3_i);
    // fprintf(stderr, "[M::%s] Heterozygous k-mer peak: %d\n", __func__, asm_opt.het_cov);
    // fprintf(stderr, "[M::%s] Homozygous k-mer peak: %d\n", __func__, asm_opt.hom_cov);
    // fprintf(stderr, "[M::%s] Heterozygous coverage peak: %lld\n", __func__, (*het_peak));
    // fprintf(stderr, "[M::%s] Homozygous coverage peak: %lld\n", __func__, (*hom_peak));
    // fprintf(stderr, "[M::%s] Alter coverage peak: %lld\n", __func__, topo_peak_i);
}




long long get_alter_peak(ma_ug_t *ug, asg_t *read_g, R_to_U* ruIndex, uint64_t* position_index, 
ma_hit_t_alloc* sources, ma_sub_t* coverage_cut, long long cov_buf_length)
{
    
    ma_utg_t* u = NULL;
    asg_t* nsg = ug->g;
    uint64_t v, j, k, qn, n_vtx = nsg->n_seq, primary_bases = 0, alter_bases = 0;
    uint32_t tn, is_Unitig;
    long long* cov_buf = NULL;
    ma_hit_t *h;
    cov_buf = (long long*)calloc(cov_buf_length, sizeof(long long));
    long long R_bases = 0, C_bases_primary = 0, C_bases_alter = 0, C_bases = 0;
    memset(position_index, -1, sizeof(uint64_t)*read_g->n_seq);


    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].del) continue;
        if(nsg->seq[v].c == ALTER_LABLE) continue;
        u = &(ug->u.a[v]);
        if(u->m == 0) continue;
        for (k = 0; k < u->n; k++)
        {
            qn = u->a[k]>>33;
            position_index[qn] = 0;
            R_bases = coverage_cut[qn].e - coverage_cut[qn].s;
            primary_bases += R_bases;
        }
    }

    for (qn = 0; qn < read_g->n_seq; qn++)
    {
        if(position_index[qn] == 0) continue;
        if(read_g->seq[qn].del) continue;

        C_bases = C_bases_primary = C_bases_alter = 0;
        R_bases = coverage_cut[qn].e - coverage_cut[qn].s;
        alter_bases += R_bases;
        for (j = 0; j < (uint64_t)(sources[qn].length); j++)
        {
            h = &(sources[qn].buffer[j]);
            if(h->el != 1) continue;
            tn = Get_tn((*h));

            if(read_g->seq[tn].del == 1)
            {
                ///get the id of read that contains it 
                get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
                if(tn == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[tn].del == 1) continue;
            }

            if(position_index[tn] == 0)
            {
                C_bases_primary += Get_qe((*h)) - Get_qs((*h));
            }
            else
            {
                C_bases_alter += Get_qe((*h)) - Get_qs((*h));
            }
        }

        C_bases = C_bases_primary + C_bases_alter;
        if(C_bases_alter < C_bases * ALTER_COV_THRES) continue;

        C_bases = C_bases/R_bases;
        if(C_bases < 0 || C_bases >= cov_buf_length) continue;
        cov_buf[C_bases]++;
    }
    
    long long max_i = -1, max = -1;
    for (j = 0; (long long)j < cov_buf_length; ++j)
    {
        if (cov_buf[j] > max)
        {
            max = cov_buf[j];
            max_i = j;
        }
    }

    if(alter_bases < primary_bases * REAL_ALTER_THRES) max_i = max = -1;

    free(cov_buf);
    memset(position_index, -1, sizeof(uint64_t)*read_g->n_seq);

    return max_i;
}

long long get_read_coverage_thres(ma_ug_t *ug, asg_t *read_g, R_to_U* ruIndex, uint64_t* position_index, 
ma_hit_t_alloc* sources, ma_sub_t* coverage_cut, uint64_t n_read, long long cov_buf_length,
long long* k_mer_only, long long* coverage_only)
{
    uint64_t i, j;
    long long* cov_buf = NULL;
    ma_hit_t *h;
    cov_buf = (long long*)calloc(cov_buf_length, sizeof(long long));
    long long R_bases = 0, C_bases = 0;
    for (i = 0; i < n_read; ++i) 
    {
        C_bases = 0;
        R_bases = coverage_cut[i].e - coverage_cut[i].s;
        for (j = 0; j < (uint64_t)(sources[i].length); j++)
        {
            h = &(sources[i].buffer[j]);
            if(h->el != 1) continue;
            C_bases += Get_qe((*h)) - Get_qs((*h));
        }
        C_bases = C_bases/R_bases;
        if(C_bases < 0 || C_bases >= cov_buf_length) continue;
        cov_buf[C_bases]++;
    }

    long long alter_peak = -1, hom_peak = -1, het_peak = -1;
    if(position_index)
    {
        alter_peak = get_alter_peak(ug, read_g, ruIndex, position_index, sources, coverage_cut, 
        cov_buf_length);
    }
    
    get_read_peak(read_g, cov_buf, cov_buf_length, alter_peak == -1? NULL: &alter_peak, &hom_peak, &het_peak,
    k_mer_only, coverage_only, asm_opt.hg_size);

    free(cov_buf);

    if(hom_peak != -1) return hom_peak*HOM_PEAK_RATE;
    if(het_peak != -1) return het_peak*HET_PEAK_RATE;
    return -1;
}










void init_hap_alignment_struct(hap_alignment_struct* x, uint32_t size)
{
    x->vote_counting = (uint64_t*)malloc(sizeof(uint64_t)*size);
    memset(x->vote_counting, 0, sizeof(uint64_t)*size);

    x->visit = (uint8_t*)malloc(sizeof(uint8_t)*size);
    memset(x->visit, 0, size);

    kv_init(x->u_vecs.a);
    kv_init(x->u_buffer.a);
    kv_init(x->u_buffer_tailIndex.a);
    kv_init(x->u_buffer_prevIndex.a);
    kv_init(x->u_buffer_beg.a);
    kv_init(x->u_buffer_flag.a);
    kv_init(x->u_can.a);
}

void destory_hap_alignment_struct(hap_alignment_struct* x)
{
    free(x->vote_counting);
    free(x->visit);
    kv_destroy(x->u_vecs.a);
    kv_destroy(x->u_buffer.a);
    kv_destroy(x->u_buffer_tailIndex.a);
    kv_destroy(x->u_buffer_prevIndex.a);
    kv_destroy(x->u_buffer_beg.a);
    kv_destroy(x->u_buffer_flag.a);
    kv_destroy(x->u_can.a);
}

uint8_t *init_pip_hh(asg_t *rg, ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, long long sc)
{
    uint8_t *c = NULL; CALLOC(c, rg->n_seq);
    ma_hit_t *h = NULL;
    uint32_t i, k;
    long long R_Base, C_Base;
    for (i = 0; i < rg->n_seq; i++) {
        if(sc <= 0){
            c[i] = 1;
            continue;
        }
        R_Base = (coverage_cut[i].e - coverage_cut[i].s);
        for (k = 0, C_Base = 0; k < reverse_sources[i].length; k++){
            h = &(reverse_sources[i].buffer[k]);
            C_Base += (Get_qe((*h)) - Get_qs((*h)));
        }
        C_Base = (R_Base!=0?(C_Base/R_Base):0);
        C_Base /= sc;
        c[i] = 1;
        if(C_Base < REV_W) c[i] = REV_W -  C_Base;
    }
    return c;
}

void init_hap_alignment_struct_pip(hap_alignment_struct_pip* x, uint32_t num_threads, uint32_t n_seq,
ma_ug_t *ug, asg_t *read_g,  ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex, ma_sub_t *coverage_cut, 
uint64_t* position_index, float Hap_rate, int max_hang, int min_ovlp, float chain_rate, hap_overlaps_list* all_ovlp, hap_cov_t *cov)
{
    uint32_t i;
    x->num_threads = num_threads;
    x->buf = (hap_alignment_struct*)malloc(sizeof(hap_alignment_struct)*x->num_threads);
    for (i = 0; i < x->num_threads; i++)
    {
        init_hap_alignment_struct(&(x->buf[i]), n_seq);
    }

    x->ug = ug;
    x->read_g = read_g;
    x->sources = sources;
    x->reverse_sources = reverse_sources;
    x->ruIndex = ruIndex;
    x->coverage_cut = coverage_cut;
    x->position_index = position_index;
    x->Hap_rate = Hap_rate;
    x->max_hang = max_hang;
    x->min_ovlp = min_ovlp;
    x->chain_rate = chain_rate;
    x->all_ovlp = all_ovlp;
    x->cov = cov;
    long long sc = -1;
    if(asm_opt.hom_global_coverage_set) {
        sc = asm_opt.hom_global_coverage*1.75;
    }
    else {
        if(asm_opt.hom_global_coverage > 0){
            sc = ((int)(((double)asm_opt.hom_global_coverage)/((double)HOM_PEAK_RATE)))*1.75;
        }
    }
    if(sc <= 0) sc = -1;
    x->hh = init_pip_hh(read_g, reverse_sources, coverage_cut, sc);
}


void destory_hap_alignment_struct_pip(hap_alignment_struct_pip* x)
{
    uint32_t i;
    for (i = 0; i < x->num_threads; i++)
    {
        destory_hap_alignment_struct(&(x->buf[i]));
    }

    free(x->buf);
    free(x->hh);
}

void init_hap_overlaps_list(hap_overlaps_list* x, uint32_t num)
{
    uint32_t i = 0;
    x->num = num;
    x->x = (kvec_hap_overlaps*)malloc(sizeof(kvec_hap_overlaps)*x->num);
    for (i = 0; i < x->num; i++)
    {
        kv_init(x->x[i].a);
    }
}

void enable_debug_mode(uint32_t mode)
{
    debug_enable = mode;
}

void destory_hap_overlaps_list(hap_overlaps_list* x)
{
    uint32_t i = 0;
    for (i = 0; i < x->num; i++)
    {
        kv_destroy(x->x[i].a);
    }
    free(x->x);
}



inline void clean_visit_flag(uint8_t* visit, asg_t *read_g, R_to_U* ruIndex, uint32_t contigNum, 
ma_hit_t_alloc* x)
{
    uint32_t k, rId, is_Unitig, Hap_cId;

    if(x->length*2 > contigNum)
    {
        memset(visit, 0, contigNum);
    }
    else
    {
        for (k = 0; k < x->length; k++)
        {
            rId = Get_tn(x->buffer[k]);
        
            if(read_g->seq[rId].del == 1)
            {
                ///get the id of read that contains it 
                get_R_to_U(ruIndex, rId, &rId, &is_Unitig);
                if(rId == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[rId].del == 1) continue;
            }

            ///there are two cases: 
            ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
            ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
            get_R_to_U(ruIndex, rId, &Hap_cId, &is_Unitig);
            if(is_Unitig == 0 || Hap_cId == (uint32_t)-1) continue;
            ///here rId is the id of the read coming from the different haplotype
            ///Hap_cId is the id of the corresponding contig (note here is the contig, instead of untig)
            visit[Hap_cId] = 0;
        }
    }
}

uint32_t prefilter(uint32_t x_pos, uint32_t y_pos, uint32_t xLen, uint32_t yLen, uint32_t dir, 
float Hap_rate, uint32_t seedOcc)
{
    uint32_t max_count = 0, min_count = 0;
    uint32_t /**xLeftBeg, **/xLeftLen, yLeftBeg, yLeftLen;
    uint32_t xRightBeg, xRightLen, yRightBeg, yRightLen;
    if(dir == 0)
    {
        /**xLeftBeg = 0;**/ xLeftLen = x_pos; xRightBeg = x_pos; xRightLen = xLen - xRightBeg;
        yLeftBeg = 0; yLeftLen = y_pos; yRightBeg = y_pos; yRightLen = yLen - yRightBeg;
    }
    else
    {
        /**xLeftBeg = 0;**/ xLeftLen = x_pos; xRightBeg = x_pos; xRightLen = xLen - xRightBeg;

        yLeftBeg = y_pos + 1; yLeftLen = yLen - yLeftBeg;
        yRightBeg = 0; yRightLen = y_pos + 1;
    }



    max_count = seedOcc;
    min_count = MIN(xLeftLen, yLeftLen) + MIN(xRightLen, yRightLen);
    if(min_count == 0) return NON_PLOID;
    if(max_count <= min_count*Hap_rate) return NON_PLOID;
    return PLOID;
}



inline uint64_t get_xy_pos(asg_t *read_g, asg_arc_t* t, uint32_t v_in_unitig, uint32_t w_in_unitig, 
uint32_t xUnitigLen, uint32_t yUnitigLen, uint64_t* position_index, uint8_t* rev)
{
    uint32_t x_pos, y_pos, x_dir = 0, y_dir = 0;
    uint64_t tmp;
    x_pos = y_pos = (uint32_t)-1;
    if((t->ul>>32)==v_in_unitig)///end pos
    {
        x_pos = (position_index[v_in_unitig>>1]>>32) + read_g->seq[v_in_unitig>>1].len - 1;
        x_dir = 0;
    }
    else if((t->ul>>32)==(v_in_unitig^1))///start pos
    {
        x_pos = (position_index[v_in_unitig>>1]>>32);
        x_dir = 1;
    }
    else
    {
        fprintf(stderr, "ERROR\n");
    }
    
    if(t->v == w_in_unitig)
    {
        y_pos = (position_index[w_in_unitig>>1]>>32) + t->ol - 1;
        y_dir = 0;
    }
    else if(t->v == (w_in_unitig^1))
    {
        y_pos = (position_index[w_in_unitig>>1]>>32) + read_g->seq[w_in_unitig>>1].len - t->ol;
        y_dir = 1;
    }
    else
    {
        fprintf(stderr, "ERROR\n");
    }

    (*rev) = x_dir^y_dir;
    if((*rev))
    {
        if(yUnitigLen <= y_pos)
        {
            y_pos = (uint32_t)-1;
        }
        else
        {
            y_pos = yUnitigLen - y_pos - 1;
        }
    } 

    if(x_pos>=xUnitigLen) x_pos = (uint32_t)-1;
    if(y_pos>=yUnitigLen) y_pos = (uint32_t)-1;
    
    tmp = x_pos; tmp = tmp << 32; tmp = tmp | y_pos;
    return tmp;
}


void print_debug_unitig(ma_utg_t *xReads, uint64_t* position_index, const char* infor)
{
    uint32_t k;
    fprintf(stderr, "\n%s: n = %u\n", infor, xReads->n);
    for (k = 0; k < xReads->n; k++)
    {
        fprintf(stderr, "(%u)v: %u, len: %u, index: %u, pos: %u\n", 
        k, (uint32_t)(xReads->a[k]>>32), (uint32_t)xReads->a[k], (uint32_t)(position_index[xReads->a[k]>>33]),
        (uint32_t)(position_index[xReads->a[k]>>33]>>32));
    }
}

void deduplicate_edge(kvec_asg_arc_t_offset* u_buffer)
{
    if(u_buffer->a.n == 0) return;
    long long i = u_buffer->a.n - 1, k, i_off;
    uint32_t v = u_buffer->a.a[i].x.ul>>33, m;

    for (; i >= 0; i--) {
        if((u_buffer->a.a[i].x.ul>>33) != v) break;
    }

    i = i + 1;
    for (m = i; i < (long long)u_buffer->a.n; i++)
    {
        if(u_buffer->a.a[i].x.del) continue;

        i_off = Cal_Off(u_buffer->a.a[i].Off);

        for (k = i + 1; k < (long long)u_buffer->a.n; k++)
        {
            if(u_buffer->a.a[k].x.del) continue;
            if(u_buffer->a.a[i].x.el != u_buffer->a.a[k].x.el) continue;
            if(i_off != Cal_Off(u_buffer->a.a[k].Off)) continue;
            u_buffer->a.a[k].x.del = 1;
            u_buffer->a.a[i].weight += u_buffer->a.a[k].weight;
        }
        u_buffer->a.a[m] = u_buffer->a.a[i];
        m++;
    }

    u_buffer->a.n = m;
    ///fprintf(stderr, "u_buffer->a.n: %u, i: %lld\n", u_buffer->a.n, i);
}

int cmp_hap_alignment(const void * a, const void * b)
{
    if((*(asg_arc_t_offset*)a).x.el > (*(asg_arc_t_offset*)b).x.el) return 1;
    if((*(asg_arc_t_offset*)a).x.el < (*(asg_arc_t_offset*)b).x.el) return -1;

    long long aOff = Cal_Off((*(asg_arc_t_offset*)a).Off);
    long long bOff = Cal_Off((*(asg_arc_t_offset*)b).Off);

    if(aOff > bOff) return 1;
    if(aOff < bOff) return -1;

    if(((*(asg_arc_t_offset*)a).Off>>32) > ((*(asg_arc_t_offset*)b).Off>>32)) return 1;
    if(((*(asg_arc_t_offset*)a).Off>>32) < ((*(asg_arc_t_offset*)b).Off>>32)) return -1;

    if((uint32_t)((*(asg_arc_t_offset*)a).Off) > (uint32_t)((*(asg_arc_t_offset*)b).Off)) return 1;
    if((uint32_t)((*(asg_arc_t_offset*)a).Off) < (uint32_t)((*(asg_arc_t_offset*)b).Off)) return -1;

    if((*(asg_arc_t_offset*)a).weight < (*(asg_arc_t_offset*)b).weight) return 1;
    if((*(asg_arc_t_offset*)a).weight > (*(asg_arc_t_offset*)b).weight) return -1;

    return 0;
}


int cmp_hap_alignment_chaining(const void * a, const void * b)
{
    if((*(asg_arc_t_offset*)a).x.el > (*(asg_arc_t_offset*)b).x.el) return 1;
    if((*(asg_arc_t_offset*)a).x.el < (*(asg_arc_t_offset*)b).x.el) return -1;

    if(((*(asg_arc_t_offset*)a).Off>>32) > ((*(asg_arc_t_offset*)b).Off>>32)) return 1;
    if(((*(asg_arc_t_offset*)a).Off>>32) < ((*(asg_arc_t_offset*)b).Off>>32)) return -1;

    if((uint32_t)((*(asg_arc_t_offset*)a).Off) > (uint32_t)((*(asg_arc_t_offset*)b).Off)) return 1;
    if((uint32_t)((*(asg_arc_t_offset*)a).Off) < (uint32_t)((*(asg_arc_t_offset*)b).Off)) return -1;

    return 0;
}

int cmp_hap_candidates(const void * a, const void * b)
{
    if((*(hap_candidates*)a).weight < (*(hap_candidates*)b).weight) return 1;
    if((*(hap_candidates*)a).weight > (*(hap_candidates*)b).weight) return -1;

    if((*(hap_candidates*)a).index_beg > (*(hap_candidates*)b).index_beg) return 1;
    if((*(hap_candidates*)a).index_beg < (*(hap_candidates*)b).index_beg) return -1;
    
    return 0;
}
inline long long get_hap_overlapLen(long long x_beg, long long x_end, long long xLen, 
long long y_beg, long long y_end, long long yLen, long long* n_x_beg, long long* n_x_end,
long long* n_y_beg, long long* n_y_end)
{
    if(x_beg <= y_beg)
    {
        y_beg = y_beg - x_beg;
        x_beg = 0;
    }
    else
    {
        x_beg = x_beg - y_beg;
        y_beg = 0;
    }

    long long x_right_length = xLen - x_end - 1;
    long long y_right_length = yLen - y_end - 1;


    if(x_right_length <= y_right_length)
    {
        x_end = xLen - 1;
        y_end = y_end + x_right_length;        
    }
    else
    {
        x_end = x_end + y_right_length;
        y_end = yLen - 1;
    }

    if(n_x_beg) (*n_x_beg) = x_beg;
    if(n_x_end) (*n_x_end) = x_end;
    if(n_y_beg) (*n_y_beg) = y_beg;
    if(n_y_end) (*n_y_end) = y_end;

    return x_end - x_beg + 1;
}



uint32_t classify_hap_overlap(long long xBeg, long long xEnd, long long xLen,
long long yBeg, long long yEnd, long long yLen, long long* r_xBeg, long long* r_xEnd, 
long long* r_yBeg, long long* r_yEnd)
{
    long long n_x_beg, n_x_end, n_y_beg, n_y_end;
    get_hap_overlapLen(xBeg, xEnd, xLen, yBeg, yEnd, yLen, &n_x_beg, &n_x_end, &n_y_beg, &n_y_end);
    if(r_xBeg) (*r_xBeg) = n_x_beg;
    if(r_xEnd) (*r_xEnd) = n_x_end;
    if(r_yBeg) (*r_yBeg) = n_y_beg;
    if(r_yEnd) (*r_yEnd) = n_y_end;
    if(n_x_beg == 0 && n_x_end == xLen - 1) return YCX;
    if(n_y_beg == 0 && n_y_end == yLen - 1) return XCY;
    if(n_y_beg == 0 && n_x_end == xLen - 1) return X2Y;
    if(n_x_beg == 0 && n_y_end == yLen - 1) return Y2X;
    return XCY;
}


uint64_t get_pair_hap_coverage(uint64_t* readIDs, uint32_t Len, ma_hit_t_alloc* sources, ma_sub_t* coverage_cut)
{
    uint32_t m, n, qn;
    ma_hit_t *h;
    uint64_t R_bases = 0, C_bases = 0;

    for (m = 0; m < Len; m++)
    {
        qn = readIDs[m]>>33;
        R_bases += coverage_cut[qn].e - coverage_cut[qn].s;
        for (n = 0; n < (uint64_t)(sources[qn].length); n++)
        {
            h = &(sources[qn].buffer[n]);
            C_bases += Get_qe((*h)) - Get_qs((*h));
        }
    }

    return C_bases/R_bases;
}


uint64_t get_pair_purge_coverage(ma_utg_t *xReads, long long xPosBeg, long long xPosEnd, 
ma_utg_t *yReads, long long yPosBeg, long long yPosEnd, uint32_t rev, asg_t *read_g, hap_cov_t *cov)
{
    long long offset, r_beg, r_end, i_beg, i_end, ovlp, IdxBeg, IdxEnd;
    uint64_t i, rId, uCov, uLen;
    ma_utg_t *x = NULL;
    uCov = uLen = 0;
    if(rev) 
    {
        yPosBeg = yReads->len - yPosBeg - 1;
        yPosEnd = yReads->len - yPosEnd - 1;
        offset = yPosBeg; yPosBeg = yPosEnd; yPosEnd = offset;
    }
    

    IdxBeg = IdxEnd = -1; 
    x = xReads; i_beg = xPosBeg; i_end = xPosEnd;
    for (i = 0, offset = 0; i < x->n; i++)
    {
        rId = x->a[i]>>33;
        r_beg = offset; r_end = offset + (long long)(read_g->seq[rId].len) - 1;
        offset += (uint32_t)x->a[i];
        ovlp = (long long)(MIN(r_end, i_end)) - (long long)(MAX(r_beg, i_beg)) + 1;
        if(ovlp <= 0 || ovlp < read_g->seq[rId].len * 0.8)
        {
            if(IdxBeg != -1 && IdxEnd != -1) break;
            continue;
        }

        if(IdxBeg == -1) IdxBeg = i;
        IdxEnd = i;
    }
    if(IdxBeg != -1 && IdxEnd != -1)
    {
        for (i = IdxBeg; (long long)i <= IdxEnd; i++)
        {
            rId = x->a[i]>>33;
            uCov += cov->cov[rId];
            uLen += cov->read_g->seq[rId].len;
        }
    }    



    IdxBeg = IdxEnd = -1; 
    x = yReads; i_beg = yPosBeg; i_end = yPosEnd;
    for (i = 0, offset = 0; i < x->n; i++)
    {
        rId = x->a[i]>>33;
        r_beg = offset; r_end = offset + (long long)(read_g->seq[rId].len) - 1;
        offset += (uint32_t)x->a[i];
        ovlp = (long long)(MIN(r_end, i_end)) - (long long)(MAX(r_beg, i_beg)) + 1;
        if(ovlp <= 0 || ovlp < read_g->seq[rId].len * 0.8)
        {
            if(IdxBeg != -1 && IdxEnd != -1) break;
            continue;
        }

        if(IdxBeg == -1) IdxBeg = i;
        IdxEnd = i;
    }
    if(IdxBeg != -1 && IdxEnd != -1)
    {
        for (i = IdxBeg; (long long)i <= IdxEnd; i++)
        {
            rId = x->a[i]>>33;
            uCov += cov->cov[rId];
            uLen += cov->read_g->seq[rId].len;
        }
    }

    return (uLen == 0? 0 : uCov / uLen);
}


void get_pair_hap_similarity_by_base(ma_utg_t *xReads, asg_t *read_g, uint32_t target_uId, 
ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex, long long xBegPos, long long xEndPos, 
double* Match, double* Total)
{
    uint32_t i, j, qn, tn, is_Unitig, uId, min_count = 0, max_count = 0;
    long long offset, r_beg, r_end, ovlp;
    
    for (i = 0, offset = 0; i < xReads->n; i++)
    {
        qn = xReads->a[i]>>33;
        r_beg = offset; r_end = offset + (long long)(read_g->seq[qn].len) - 1;
        offset += (uint32_t)xReads->a[i];
        
        ovlp = (long long)(MIN(r_end, xEndPos)) - (long long)(MAX(r_beg, xBegPos)) + 1;
        if(ovlp <= 0) continue;

        if(reverse_sources[qn].length > 0) min_count++;
        if(reverse_sources[qn].length == 0) continue;
        for (j = 0; j < reverse_sources[qn].length; j++)
        {
            tn = Get_tn(reverse_sources[qn].buffer[j]);
            if(read_g->seq[tn].del == 1)
            {
                get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
                if(tn == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[tn].del == 1) continue;
            }


            get_R_to_U(ruIndex, tn, &uId, &is_Unitig);
            if(uId!=(uint32_t)-1 && is_Unitig == 1 && uId == target_uId)
            {
                max_count++;
                break;
            }
        }
    }

    (*Match) = max_count;
    (*Total) = min_count;
}

void get_pair_hap_similarity(uint64_t* readIDs, uint32_t Len, uint32_t target_uId, 
ma_hit_t_alloc* reverse_sources, asg_t *read_g, R_to_U* ruIndex, double* Match, double* Total)
{
    #define CUTOFF_THRES 1000
    uint32_t i, j, qn, tn, is_Unitig, uId, min_count = 0, max_count = 0, cutoff = 0;;
    for (i = 0; i < Len; i++)
    {
        if(cutoff > CUTOFF_THRES && cutoff > (Len>>1))
        {
            max_count = 0;
            min_count = Len;
            break;
        }
        qn = readIDs[i]>>33;
        if(reverse_sources[qn].length > 0) min_count++;
        if(reverse_sources[qn].length == 0) continue;
        for (j = 0; j < reverse_sources[qn].length; j++)
        {
            tn = Get_tn(reverse_sources[qn].buffer[j]);
            if(read_g->seq[tn].del == 1)
            {
                get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
                if(tn == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[tn].del == 1) continue;
            }


            get_R_to_U(ruIndex, tn, &uId, &is_Unitig);
            if(uId!=(uint32_t)-1 && is_Unitig == 1 && uId == target_uId)
            {
                max_count++;
                break;
            }
        }

        //means no match
        if(j == reverse_sources[qn].length)
        {
            cutoff++;
        }
        else
        {
            cutoff = 0;
        }
    }

    (*Match) = max_count;
    (*Total) = min_count;
}
/**
void get_pair_hap_similarity_deduplicate(uint64_t* readIDs, uint32_t Len, uint32_t target_uId, 
ma_hit_t_alloc* reverse_sources, asg_t *read_g, R_to_U* ruIndex, double* Match, double* Total)
{
    get_pair_hap_similarity(readIDs, Len, target_uId, reverse_sources, read_g, ruIndex, Match, Total);
    return;

    #define CUTOFF_THRES 100
    uint32_t i, j, qn, tn, is_Unitig, uId, min_count = 0, max_count = 0, cutoff = 0, is_found;
    for (i = 0; i < Len; i++)
    {
        if(cutoff > CUTOFF_THRES)
        {
            max_count = 0;
            min_count = Len;
            break;
        }
        qn = readIDs[i]>>33;
        is_found = 0;

        for (j = 0; j < reverse_sources[qn].length; j++)
        {
            tn = Get_tn(reverse_sources[qn].buffer[j]);
            if(read_g->seq[tn].del == 1)
            {
                get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
                if(tn == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[tn].del == 1) continue;
            }


            get_R_to_U(ruIndex, tn, &uId, &is_Unitig);
            if(uId!=(uint32_t)-1 && is_Unitig == 1 && uId == target_uId)
            {
                max_count++;
            }
            min_count++;
            is_found = 1;
        }

        //means there is a match
        if(is_found)
        {
            cutoff = 0;
        }
        else
        {
            cutoff++;
        }
    }

    (*Match) = max_count;
    (*Total) = min_count;
}
**/

inline void check_hap_match(uint32_t qn, uint32_t targetBeg, uint32_t targetEnd, uint32_t targetID, 
uint64_t* position_index, ma_hit_t_alloc* reverse_sources, asg_t *read_g, R_to_U* ruIndex, uint32_t* is_found, uint32_t* is_match)
{
    uint32_t j, tn, uId, is_Unitig, offset;
    (*is_found) = (*is_match) = 0;
    if(reverse_sources[qn].length > 0) (*is_found) = 1;
    for (j = 0; j < reverse_sources[qn].length; j++)
    {
        tn = Get_tn(reverse_sources[qn].buffer[j]);
        if(read_g->seq[tn].del == 1)
        {
            get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
            if(tn == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[tn].del == 1) continue;
        }


        get_R_to_U(ruIndex, tn, &uId, &is_Unitig);
        if(uId!=(uint32_t)-1 && is_Unitig == 1 && uId == targetID)
        {
            offset = (uint32_t)(position_index[tn]);
            if(offset >= targetBeg && offset <= targetEnd)
            {
                (*is_match) = 1;
                break;
            }                    
        }
    }

}

/**
inline void check_hap_match_deduplicate(uint32_t qn, uint32_t targetBeg, uint32_t targetEnd, uint32_t targetID, 
uint64_t* position_index, ma_hit_t_alloc* reverse_sources, asg_t *read_g, R_to_U* ruIndex, uint32_t* is_found, uint32_t* is_match)
{
    check_hap_match(qn, targetBeg, targetEnd, targetID, position_index, reverse_sources, read_g, 
    ruIndex, is_found, is_match);
    return;


    uint32_t j, tn, uId, is_Unitig, offset;
    (*is_found) = (*is_match) = 0;

    for (j = 0; j < reverse_sources[qn].length; j++)
    {
        tn = Get_tn(reverse_sources[qn].buffer[j]);
        if(read_g->seq[tn].del == 1)
        {
            get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
            if(tn == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[tn].del == 1) continue;
        }


        get_R_to_U(ruIndex, tn, &uId, &is_Unitig);
        if(uId!=(uint32_t)-1 && is_Unitig == 1 && uId == targetID)
        {
            offset = (uint32_t)(position_index[tn]);
            if(offset >= targetBeg && offset <= targetEnd)
            {
                (*is_match)++;
            }                    
        }
        (*is_found)++;
    }

}
**/

void determin_hap_alignment_boundary_single_side(uint64_t* readIDs, long long queryLen, long long targetBeg, 
long long targetEnd, long long targetID, long long eMatch, long long eTotal, long long dir, 
float H_rate, int is_local, uint64_t* position_index, ma_hit_t_alloc* reverse_sources, asg_t *read_g, 
R_to_U* ruIndex, uint32_t* n_matchLen, uint32_t* n_max_count, uint32_t* n_min_count)
{
    if(queryLen == 0)
    {
        (*n_matchLen) = (*n_min_count) = (*n_max_count) = 0;
        return;
    }
    long long i, maxId, min_count = eTotal, max_count = eMatch, matchLen = 0;
    long long rLen, score = 0, max_score = 0;
    uint32_t is_found, is_match;
    if(dir == 0)
    {
        for (i = 0, maxId = 0; i < queryLen; i++)
        {
            
            check_hap_match(readIDs[i]>>33, targetBeg, targetEnd, targetID, position_index, reverse_sources, 
            read_g, ruIndex, &is_found, &is_match);

            min_count += is_found;
            max_count += is_match;
            if(max_count > min_count*H_rate) maxId = i;

            if(is_local && is_found)
            {
                rLen = read_g->seq[readIDs[i]>>33].len;
                score += (is_match? rLen : (rLen*(-1)));
                if(score >= max_score) max_score = score, maxId = i;
            }
        }

        for (i = maxId; i >= 0; i--)
        {
            check_hap_match(readIDs[i]>>33, targetBeg, targetEnd, targetID, position_index, reverse_sources, 
            read_g, ruIndex, &is_found, &is_match);

            ///if(is_found > 0 && is_match > 0 && is_match > is_found*Hap_rate)
            if(is_found == 1 && is_match == 1)
            {
                break;
            }
            min_count -= is_found;
            max_count -= is_match;
        }
        
        matchLen = i+1;
    }
    else
    {
        for (i = queryLen - 1, maxId = queryLen - 1; i >= 0; i--)
        {
            check_hap_match(readIDs[i]>>33, targetBeg, targetEnd, targetID, position_index, reverse_sources, 
            read_g, ruIndex, &is_found, &is_match);

            min_count += is_found;
            max_count += is_match;
            if(max_count > min_count*H_rate) maxId = i;

            if(is_local && is_found)
            {
                rLen = read_g->seq[readIDs[i]>>33].len;
                score += (is_match? rLen : (rLen*(-1)));
                if(score >= max_score) max_score = score, maxId = i;
            }
        }

        for (i = maxId; i < queryLen; i++)
        {
            check_hap_match(readIDs[i]>>33, targetBeg, targetEnd, targetID, position_index, reverse_sources, 
            read_g, ruIndex, &is_found, &is_match);

            ///if(is_found > 0 && is_match > 0 && is_match > is_found*Hap_rate)
            if(is_found == 1 && is_match == 1)
            {
                break;
            }
            min_count -= is_found;
            max_count -= is_match;
        }

        matchLen = queryLen - i;
    }

    ///need to check if min_count == 0
    if(min_count == 0)
    {
        (*n_matchLen) = (*n_min_count) = (*n_max_count) = 0;
        return;
    }

    (*n_matchLen) = matchLen;
    (*n_min_count) = min_count;
    (*n_max_count) = max_count;
}

inline void modify_target_interval(long long beg, long long end, long long len,
long long* target_beg, long long* target_end)
{
    #define TARGET_SGIFT 3
    beg -= TARGET_SGIFT;
    end += TARGET_SGIFT;
    if(beg < 0) beg = 0;
    if(end >= len) end = len - 1;
    (*target_beg) = beg;
    (*target_end) = end;
}


void bi_direction_hap_alignment_extention(ma_utg_t* xReads, uint32_t xLeftBeg, uint32_t xLeftLen,
uint32_t xRightBeg, uint32_t xRightLen, uint32_t targetUid, uint32_t target_beg, uint32_t target_end,
float Hap_rate, int is_local, uint64_t* position_index, ma_hit_t_alloc* reverse_sources, asg_t *read_g, R_to_U* ruIndex,
uint32_t rev, long long* x_interval_beg, long long* x_interval_end)
{
    if(rev)
    {
        uint32_t k;
        k = xLeftBeg; xLeftBeg = xRightBeg; xRightBeg = k;
        k = xLeftLen; xLeftLen = xRightLen; xRightLen = k;
    }
    uint32_t n_matchLenLeft, x_max_countLeft, x_min_countLeft;
    uint32_t n_matchLenRight, x_max_countRight, x_min_countRight;
    n_matchLenLeft = x_max_countLeft = x_min_countLeft = 0;
    determin_hap_alignment_boundary_single_side(xReads->a+xLeftBeg, xLeftLen,
    target_beg, target_end, targetUid, x_max_countLeft, x_min_countLeft, 1, Hap_rate, is_local,
    position_index, reverse_sources, read_g, ruIndex, &n_matchLenLeft, &x_max_countLeft, 
    &x_min_countLeft);

    n_matchLenRight = x_max_countRight = x_min_countRight = 0;
    determin_hap_alignment_boundary_single_side(xReads->a+xRightBeg, xRightLen,
    target_beg, target_end, targetUid, x_max_countRight, x_min_countRight, 0, Hap_rate, is_local,
    position_index, reverse_sources, read_g, ruIndex, &n_matchLenRight, &x_max_countRight,
    &x_min_countRight);

    if(x_max_countLeft >= x_max_countRight)
    {
        determin_hap_alignment_boundary_single_side(xReads->a+xRightBeg, xRightLen,
        target_beg, target_end, targetUid, x_max_countLeft, x_min_countLeft, 0, Hap_rate, is_local,
        position_index, reverse_sources, read_g, ruIndex, &n_matchLenRight, &x_max_countRight,
        &x_min_countRight);
    }
    else
    {
        determin_hap_alignment_boundary_single_side(xReads->a+xLeftBeg, xLeftLen,
        target_beg, target_end, targetUid, x_max_countRight, x_min_countRight, 1, Hap_rate, is_local,
        position_index, reverse_sources, read_g, ruIndex, &n_matchLenLeft, &x_max_countLeft, 
        &x_min_countLeft);
    }

    (*x_interval_beg) = xLeftBeg + xLeftLen; (*x_interval_beg) -= n_matchLenLeft;
    (*x_interval_end) = xRightBeg + n_matchLenRight; (*x_interval_end) -= 1;
}

void get_hap_alignment_boundary(ma_utg_t* xReads, ma_utg_t* yReads, uint32_t type, 
uint32_t xLeftMatch, uint32_t xLeftTotal, uint32_t yLeftMatch, uint32_t yLeftTotal, 
uint32_t xRightMatch, uint32_t xRightTotal, uint32_t yRightMatch, uint32_t yRightTotal,
uint32_t xLeftBeg, uint32_t xLeftLen, uint32_t yLeftBeg, uint32_t yLeftLen,
uint32_t xRightBeg, uint32_t xRightLen, uint32_t yRightBeg, uint32_t yRightLen,
uint32_t xUid, uint32_t yUid, float Hap_rate, int is_local, uint64_t* position_index,
ma_hit_t_alloc* reverse_sources, asg_t *read_g, R_to_U* ruIndex, uint32_t rev,
long long* r_x_interval_beg, long long* r_x_interval_end, 
long long* r_y_interval_beg, long long* r_y_interval_end)
{
    
    uint32_t x_max_count, x_min_count, y_max_count, y_min_count, n_matchLen;
    long long x_interval_beg, x_interval_end, y_interval_beg, y_interval_end;
    long long target_beg, target_end;
    x_max_count = x_min_count = y_max_count = y_min_count = 0;
    if(type == X2Y)
    {
        /********************x*********************/
        x_max_count = xRightMatch;
        x_min_count = xRightTotal;

        modify_target_interval(yLeftBeg, yLeftBeg+yLeftLen-1, yReads->n, &target_beg, &target_end);
        determin_hap_alignment_boundary_single_side(xReads->a+xLeftBeg, xLeftLen, 
        /**yLeftBeg, yLeftBeg+yLeftLen-1,**/ target_beg, target_end, yUid, 
        x_max_count, x_min_count, 1, Hap_rate, is_local, position_index, reverse_sources, 
        read_g, ruIndex, &n_matchLen, &x_max_count, &x_min_count);

        x_interval_beg = xLeftBeg + xLeftLen; x_interval_beg -= n_matchLen;
        x_interval_end = xRightBeg + xRightLen; x_interval_end -= 1;
        /********************x*********************/

        /********************y*********************/
        y_max_count = yLeftMatch;
        y_min_count = yLeftTotal;
        modify_target_interval(xRightBeg, xRightBeg+xRightLen-1, xReads->n, &target_beg, &target_end);
        determin_hap_alignment_boundary_single_side(yReads->a+yRightBeg, yRightLen,
        /**xRightBeg, xRightBeg+xRightLen-1,**/ target_beg, target_end, xUid, 
        y_max_count, y_min_count, rev, Hap_rate, is_local, position_index, reverse_sources, 
        read_g, ruIndex, &n_matchLen, &y_max_count, &y_min_count);
        if(rev == 0)
        {
            y_interval_beg = yLeftBeg;
            y_interval_end = yRightBeg + n_matchLen; y_interval_end -= 1;
        }
        else
        {
            y_interval_beg = yRightBeg + yRightLen; y_interval_beg -= n_matchLen;
            y_interval_end = yLeftBeg + yLeftLen; y_interval_end -= 1;
        }
        
        /********************y*********************/
    }
    else if(type == Y2X)
    {
        /********************x*********************/
        x_max_count = xLeftMatch;
        x_min_count = xLeftTotal;
        modify_target_interval(yRightBeg, yRightBeg+yRightLen-1, yReads->n, &target_beg, &target_end);
        determin_hap_alignment_boundary_single_side(xReads->a+xRightBeg, xRightLen,
        /**yRightBeg, yRightBeg+yRightLen-1,**/ target_beg, target_end, yUid, 
        x_max_count, x_min_count, 0, Hap_rate, is_local, position_index, reverse_sources, 
        read_g, ruIndex, &n_matchLen, &x_max_count, &x_min_count);

        x_interval_beg = xLeftBeg;
        x_interval_end = xRightBeg + n_matchLen; x_interval_end -= 1;
        /********************x*********************/

        /********************y*********************/
        y_max_count = yRightMatch;
        y_min_count = yRightTotal;
        modify_target_interval(xLeftBeg, xLeftBeg+xLeftLen-1, xReads->n, &target_beg, &target_end);
        determin_hap_alignment_boundary_single_side(yReads->a+yLeftBeg, yLeftLen,
        /**xLeftBeg, xLeftBeg+xLeftLen-1,**/ target_beg, target_end, xUid, 
        y_max_count, y_min_count, 1-rev, Hap_rate, is_local, position_index, reverse_sources, 
        read_g, ruIndex, &n_matchLen, &y_max_count, &y_min_count);
        if(rev == 0)
        {
            y_interval_beg = yLeftBeg + yLeftLen; y_interval_beg -= n_matchLen;
            y_interval_end = yRightBeg + yRightLen; y_interval_end -= 1;
        }
        else
        {
            y_interval_beg = yRightBeg;
            y_interval_end = yLeftBeg + n_matchLen; y_interval_end -= 1;
        }
        /********************y*********************/
    }
    else if(type == XCY)
    {
        /********************x*********************/   
        bi_direction_hap_alignment_extention(xReads, xLeftBeg, xLeftLen, xRightBeg, xRightLen, 
        yUid, 0, yReads->n - 1, Hap_rate, is_local, position_index, reverse_sources, read_g, ruIndex, 0,
        &x_interval_beg, &x_interval_end);
        /********************x*********************/

        /********************y*********************/
        y_interval_beg = 0;
        y_interval_end = yReads->n; y_interval_end -= 1;
        /********************y*********************/
    }
    else if(type == YCX)
    {
        /********************x*********************/  
        x_interval_beg = 0;
        x_interval_end = xReads->n; x_interval_end -= 1;
        /********************x*********************/  

        /********************y*********************/
        bi_direction_hap_alignment_extention(yReads, yLeftBeg, yLeftLen, yRightBeg, yRightLen,
        xUid, 0, xReads->n - 1, Hap_rate, is_local, position_index, reverse_sources, read_g, ruIndex, rev,
        &y_interval_beg, &y_interval_end);
        /********************y*********************/
    } else abort();

    (*r_x_interval_beg) = x_interval_beg;
    (*r_x_interval_end) = x_interval_end;
    (*r_y_interval_beg) = y_interval_beg;
    (*r_y_interval_end) = y_interval_end;    
}

uint32_t vote_overlap_type(kvec_asg_arc_t_offset* u_buffer, hap_candidates* hap_can, 
uint64_t* position_index, ma_utg_t* xReads, ma_utg_t* yReads)
{
    uint32_t i, xBasePos, yBasePos;
    asg_arc_t_offset* arch = NULL;
    uint32_t flag[4];
    flag[X2Y] = flag[Y2X] = flag[XCY] = flag[YCX] = 0;


    for (i = hap_can->index_beg; i <= hap_can->index_end; i++)
    {
        arch = &(u_buffer->a.a[i]);
        xBasePos = (uint32_t)(arch->Off>>32);
        yBasePos = (uint32_t)(arch->Off);
        flag[classify_hap_overlap(xBasePos, xBasePos, xReads->len, yBasePos, yBasePos, yReads->len,
        NULL, NULL, NULL, NULL)]++;
    }
    
    uint32_t max_flag_i = 0;
    for (i = 0; i < 4; i++)
    {
        if(i == max_flag_i) continue;
        if(flag[i] > flag[max_flag_i])
        {
            max_flag_i = i;
        }
    }
    
    return max_flag_i;
}


void get_base_boundary(R_to_U* ruIndex, ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, 
asg_t *read_g, uint64_t* position_index, int max_hang, int min_ovlp, ma_utg_t *xReads, ma_utg_t *yReads,
uint32_t xUid, uint32_t yUid, long long xBegIndex, long long xEndIndex, long long yBegIndex, long long yEndIndex, 
uint32_t dir, uint32_t rev, uint32_t* x_off, uint32_t* y_off)
{
    long long k, j, offset;
    ma_hit_t_alloc *xR = NULL;
    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL, *st = NULL;
    int32_t r;
    asg_arc_t t;
    uint32_t rId, Hap_uId, is_Unitig, v, w, v_dir, w_dir, is_found = 0, oLen = 0;
    uint64_t tmp;
    (*x_off) = (*y_off) = (uint32_t)-1;
    if(dir == 1)
    {
        for (k = xEndIndex; k >= xBegIndex; k--)
        {
            xR = &(reverse_sources[xReads->a[k]>>33]);
            is_found = 0; oLen = 0;
            for (j = 0; j < xR->length; j++)
            {
                h = &(xR->buffer[j]);
                sq = &(coverage_cut[Get_qn(*h)]);
                st = &(coverage_cut[Get_tn(*h)]);
                if(st->del || read_g->seq[Get_tn(*h)].del) continue;

                r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, &t);
                ///if it is a contained overlap, skip
                if(r < 0) continue;

                rId = t.v>>1;
                if(read_g->seq[rId].del == 1) continue;
                ///there are two cases: 
                ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
                ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
                get_R_to_U(ruIndex, rId, &Hap_uId, &is_Unitig);
                if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
                if(Hap_uId != yUid) continue;

                v = xReads->a[k]>>32;
                get_R_to_U(ruIndex, v>>1, &Hap_uId, &is_Unitig);
                if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
                if(Hap_uId != xUid) continue;
                if((uint32_t)(position_index[v>>1]) != k) continue;

                w = (yReads->a[(uint32_t)(position_index[rId])])>>32;

                v_dir = ((t.ul>>32)==v)?1:0;
                w_dir = (t.v == w)?1:0;
                if(rev == 0 && v_dir != w_dir) continue;
                if(rev == 1 && v_dir == w_dir) continue;

                /****************************may have bugs********************************/
                offset = (uint32_t)(position_index[rId]);
                if(offset < yBegIndex || offset > yEndIndex) continue;
                /****************************may have bugs********************************/

                tmp = get_xy_pos(read_g, &t, v, w, xReads->len, yReads->len, position_index, &(t.el));
                if(((tmp>>32) == (uint32_t)-1) || (((uint32_t)tmp) == (uint32_t)-1)) continue;

                ///if(is_found == 0 || ((uint32_t)(tmp>>32) > (*x_off) && ((uint32_t)tmp) > (*y_off)))
                if(is_found == 0 || t.ol > oLen)
                {
                    (*x_off) = tmp>>32;
                    (*y_off) = (uint32_t)tmp;
                    oLen = t.ol;
                }

                is_found = 1; 
            }
            if(is_found) return;
        }   
    }
    else
    {
        for (k = xBegIndex; k <= xEndIndex; k++)
        {
            xR = &(reverse_sources[xReads->a[k]>>33]);
            is_found = 0; oLen = 0;
            for (j = 0; j < xR->length; j++)
            {
                h = &(xR->buffer[j]);
                sq = &(coverage_cut[Get_qn(*h)]);
                st = &(coverage_cut[Get_tn(*h)]);
                if(st->del || read_g->seq[Get_tn(*h)].del) continue;

                r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, &t);
                ///if it is a contained overlap, skip
                if(r < 0) continue;

                rId = t.v>>1;
                if(read_g->seq[rId].del == 1) continue;
                ///there are two cases: 
                ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
                ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
                get_R_to_U(ruIndex, rId, &Hap_uId, &is_Unitig);
                if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
                if(Hap_uId != yUid) continue;

                v = xReads->a[k]>>32;
                get_R_to_U(ruIndex, v>>1, &Hap_uId, &is_Unitig);
                if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
                if(Hap_uId != xUid) continue;
                if((uint32_t)(position_index[v>>1]) != k) continue;

                w = (yReads->a[(uint32_t)(position_index[rId])])>>32;

                v_dir = ((t.ul>>32)==v)?1:0;
                w_dir = (t.v == w)?1:0;
                if(rev == 0 && v_dir != w_dir) continue;
                if(rev == 1 && v_dir == w_dir) continue;
            
                /****************************may have bugs********************************/
                offset = (uint32_t)(position_index[rId]);
                if(offset < yBegIndex || offset > yEndIndex) continue;
                /****************************may have bugs********************************/

                tmp = get_xy_pos(read_g, &t, v, w, xReads->len, yReads->len, position_index, &(t.el));
                if(((tmp>>32) == (uint32_t)-1) || (((uint32_t)tmp) == (uint32_t)-1)) continue;

                ///if(is_found == 0 || ((uint32_t)(tmp>>32) < (*x_off) && ((uint32_t)tmp) < (*y_off)))
                if(is_found == 0 || t.ol > oLen)
                {
                    (*x_off) = tmp>>32;
                    (*y_off) = (uint32_t)tmp;
                    oLen = t.ol;
                }  

                is_found = 1;         
            }
            if(is_found) return;
        }
    }


    (*x_off) = (*y_off) = (uint32_t)-1;
    
}

void print_asg_arc_t_offset(asg_arc_t_offset* x, long long n, const char* info)
{
    fprintf(stderr,"\n\n(%s)n: %lld\n", info, n);
    long long i, x_off, y_off;
    for (i = 0; i < n; i++)
    {
        x_off = (long long)(x[i].Off>>32);
        y_off = (long long)((uint32_t)x[i].Off);
        fprintf(stderr, "i: %lld, x_off: %lld, y_off: %lld, weight: %lu, rev: %u, ol: %u\n", 
        i, x_off, y_off, (unsigned long)x[i].weight, x[i].x.el, x[i].x.ol);
    }
}


// Binary search 
inline int GetCeilIndex(asg_arc_t_offset* arr, kvec_t_i32_warp* T, int l, int r, uint32_t key) 
{ 
    while (r - l > 1) { 
        int m = l + (r - l) / 2; 
        if (Get_yOff(arr[T->a.a[m]].Off) >= key) 
            r = m; 
        else
            l = m; 
    } 
  
    return r; 
} 


void quick_LIS(asg_arc_t_offset* x, uint32_t n, kvec_t_i32_warp* tailIndex, kvec_t_i32_warp* prevIndex)
{
    tailIndex->a.n = prevIndex->a.n = 0;
    if(n == 0) return;

    kv_resize(int32_t, tailIndex->a, n);
    kv_resize(int32_t, prevIndex->a, n);

    long long len = 1, i, pos, m; ///the length of chain must be >=1
    tailIndex->a.a[0] = 0;
    prevIndex->a.a[0] = -1;

    ///x has already sorted by x_pos
    for(i = 1; i < (long long)n; i++) 
    { 
        if(Get_yOff(x[i].Off) < Get_yOff(x[tailIndex->a.a[0]].Off))  
        { 
            // new smallest value 
            tailIndex->a.a[0] = i; ///doesn't matter too much
        } 
        else if(Get_yOff(x[i].Off) > Get_yOff(x[tailIndex->a.a[len - 1]].Off)) 
        { 
            // arr[i] wants to extend largest subsequence 
            prevIndex->a.a[i] = tailIndex->a.a[len - 1]; 
            tailIndex->a.a[len++] = i; 
        } 
        else 
        { 
            // arr[i] wants to be a potential condidate of 
            // future subsequence 
            // It will replace ceil value in tailIndices 
            pos = GetCeilIndex(x, tailIndex, -1, len - 1, Get_yOff(x[i].Off)); 
            prevIndex->a.a[i] = pos > 0? tailIndex->a.a[pos - 1] : -1; 
            tailIndex->a.a[pos] = i; 
        } 
    }


    for (m = 0, i = tailIndex->a.a[len - 1]; m < len; i = prevIndex->a.a[i], m++)
    {
        tailIndex->a.a[len-m-1] = i;
    }

    tailIndex->a.n = len;
}

uint64_t get_xy_pos_by_pos(asg_t *read_g, asg_arc_t* t, uint32_t v_in_unitig, uint32_t w_in_unitig, 
uint32_t v_in_pos, uint32_t w_in_pos, uint32_t xUnitigLen, uint32_t yUnitigLen, uint8_t* rev)
{
    uint32_t x_pos, y_pos, x_dir = 0, y_dir = 0;
    uint64_t tmp;
    x_pos = y_pos = (uint32_t)-1;
    if((t->ul>>32)==v_in_unitig)///end pos
    {
        x_pos = v_in_pos + read_g->seq[v_in_unitig>>1].len - 1;
        x_dir = 0;
    }
    else if((t->ul>>32)==(v_in_unitig^1))///start pos
    {
        x_pos = v_in_pos;
        x_dir = 1;
    }
    else
    {
        fprintf(stderr, "ERROR\n");
    }
    
    if(t->v == w_in_unitig)
    {
        y_pos = w_in_pos + t->ol - 1;
        y_dir = 0;
    }
    else if(t->v == (w_in_unitig^1))
    {
        y_pos = w_in_pos + read_g->seq[w_in_unitig>>1].len - t->ol;
        y_dir = 1;
    }
    else
    {
        fprintf(stderr, "ERROR\n");
    }

    (*rev) = x_dir^y_dir;
    if((*rev))
    {
        if(yUnitigLen <= y_pos)
        {
            y_pos = (uint32_t)-1;
        }
        else
        {
            y_pos = yUnitigLen - y_pos - 1;
        }
    } 

    if(x_pos>=xUnitigLen) x_pos = (uint32_t)-1;
    if(y_pos>=yUnitigLen) y_pos = (uint32_t)-1;
    
    tmp = x_pos; tmp = tmp << 32; tmp = tmp | y_pos;
    return tmp;
}

void chain_trans_ovlp(hap_cov_t *cover, utg_trans_t *o, ma_ug_t *ug, asg_t *read_sg, buf_t* xReads, uint32_t targetBaseLen, uint32_t* xEnd)
{
    ma_hit_t_alloc* reverse_sources = (o? o->reverse_sources:cover->reverse_sources);
    ma_sub_t *coverage_cut = (o? o->coverage_cut:cover->coverage_cut);
    int max_hang = (o? o->max_hang:cover->max_hang);
    int min_ovlp = (o? o->min_ovlp:cover->min_ovlp);
    kvec_asg_arc_t_offset* u_buffer = (o? &(o->u_buffer):&(cover->u_buffer));
    kvec_t_i32_warp* tailIndex = (o? &(o->tailIndex):&(cover->tailIndex));
    kvec_t_i32_warp* prevIndex = (o? &(o->prevIndex):&(cover->prevIndex));
    uint64_t *pos_idx = (o? o->pos_idx:cover->pos_idx);
    ma_hit_t_alloc *xR = NULL;
    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL, *st = NULL;
    int32_t r;
    asg_arc_t t;
    uint32_t rId, v, w;
    uint64_t tmp;
    asg_arc_t_offset t_offset;
    u_buffer->a.n = 0;
    ///(*xEnd) = (uint32_t)-1;
    (*xEnd) = 0;
    uint32_t u_i, r_i, k, j, len, p_v, *a = xReads->b.a, uid, ori, l, m, aOcc, nv, xOcc = (uint32_t)-1;
    ma_utg_t* u = NULL;
    asg_arc_t *av  = NULL;


    for (u_i = r_i = len = aOcc = 0, xOcc = (uint32_t)-1, p_v = (uint32_t)-1; u_i < xReads->b.n; u_i++)
    {
        uid = a[u_i] >> 1;
        ori = a[u_i] & 1;
        u = &(ug->u.a[uid]);
        if(u->n == 0) continue;

        for (r_i = 0; r_i < u->n; r_i++, aOcc++)
        {
            l = 0;
            v = (ori == 1?((uint64_t)((u->a[u->n - r_i - 1])^(uint64_t)(0x100000000)))>>32:((uint64_t)(u->a[r_i]))>>32);

            if(p_v != (uint32_t)-1)
            {
                av = asg_arc_a(read_sg, p_v);
                nv = asg_arc_n(read_sg, p_v);
                
                for (k = 0; k < nv; k++)
                {
                    if(av[k].del) continue;
                    if(av[k].v == v) 
                    {
                        l = asg_arc_len(av[k]);
                        break;
                    }
                }
                if(k == nv) fprintf(stderr, "ERROR\n");
            }

            p_v = v; len += l;
            if(len >= targetBaseLen)
            {
                xOcc = aOcc;
                break;
            }
        }

        if(xOcc != (uint32_t)-1) break;
    }
    if(xOcc == (uint32_t)-1) xOcc = aOcc;
    if(xOcc == 0) xOcc = 1;



    for (u_i = r_i = len = aOcc = 0, p_v = (uint32_t)-1; u_i < xReads->b.n; u_i++)
    {
        uid = a[u_i] >> 1;
        ori = a[u_i] & 1;
        u = &(ug->u.a[uid]);
        if(u->n == 0) continue;

        for (r_i = 0; r_i < u->n; r_i++, aOcc++)
        {
            if(aOcc >= xOcc) break;
            l = 0;
            v = (ori == 1?((uint64_t)((u->a[u->n - r_i - 1])^(uint64_t)(0x100000000)))>>32:((uint64_t)(u->a[r_i]))>>32);

            if(p_v != (uint32_t)-1)
            {
                av = asg_arc_a(read_sg, p_v);
                nv = asg_arc_n(read_sg, p_v);
                
                for (k = 0; k < nv; k++)
                {
                    if(av[k].del) continue;
                    if(av[k].v == v) 
                    {
                        l = asg_arc_len(av[k]);
                        break;
                    }
                }
                if(k == nv) fprintf(stderr, "ERROR\n");
            }

            p_v = v; len += l;

            xR = &(reverse_sources[v>>1]);
            for (j = 0; j < xR->length; j++)
            {
                h = &(xR->buffer[j]);
                sq = &(coverage_cut[Get_qn(*h)]);
                st = &(coverage_cut[Get_tn(*h)]);
                if(st->del || read_sg->seq[Get_tn(*h)].del) continue;
                r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                            asm_opt.max_hang_rate, min_ovlp, &t);
                ///if it is a contained overlap, skip
                if(r < 0) continue;

                rId = t.v>>1;
                if(read_sg->seq[rId].del == 1) continue;
                if(pos_idx[rId] == (uint64_t)-1) continue;
                w = (uint32_t)(pos_idx[rId]);
                if(rId != (w>>1)) continue;

                tmp = get_xy_pos_by_pos(read_sg, &t, v, w, len, pos_idx[w>>1]>>32,
                                                    (uint32_t)-1, targetBaseLen, &(t.el));
                if(((tmp>>32) == (uint32_t)-1) || (((uint32_t)tmp) == (uint32_t)-1)) continue;
                if(t.el) continue; ///must

                t_offset.Off = tmp;
                t_offset.x = t;
                t_offset.weight = 1;
                kv_push(asg_arc_t_offset, u_buffer->a, t_offset); 
            }
        }

        if(aOcc >= xOcc) break;
    }

    if(u_buffer->a.n == 0) return;

    qsort(u_buffer->a.a, u_buffer->a.n, sizeof(asg_arc_t_offset), cmp_hap_alignment_chaining);

    ///print_asg_arc_t_offset(u_buffer->a.a, u_buffer->a.n, "before");


    for (k = 1, l = 0, m = 0; k <= u_buffer->a.n; ++k)
    {
        if (k == u_buffer->a.n || u_buffer->a.a[k].x.el != u_buffer->a.a[l].x.el || 
                                            u_buffer->a.a[k].Off != u_buffer->a.a[l].Off) 
        {
            u_buffer->a.a[m] = u_buffer->a.a[l];
            for (l += 1; l < k; l++)
            {
                u_buffer->a.a[m].weight += u_buffer->a.a[l].weight;
                if(u_buffer->a.a[l].x.ol > u_buffer->a.a[m].x.ol)
                {
                    u_buffer->a.a[m].x = u_buffer->a.a[l].x;
                }
            }
            l = k;
            m++;
        }
    } 
    u_buffer->a.n = m;

    ///print_asg_arc_t_offset(u_buffer->a.a, u_buffer->a.n, "after");
    quick_LIS(u_buffer->a.a, u_buffer->a.n, tailIndex, prevIndex);
    if(tailIndex->a.n == 0) return;


    uint32_t xLen_thres = (uint32_t)-1;
    asg_arc_t_offset* best = &(u_buffer->a.a[tailIndex->a.a[tailIndex->a.n-1]]);
    for (u_i = r_i = len = aOcc = 0, p_v = (uint32_t)-1; u_i < xReads->b.n; u_i++)
    {
        uid = a[u_i] >> 1;
        ori = a[u_i] & 1;
        u = &(ug->u.a[uid]);
        if(u->n == 0) continue;

        for (r_i = 0; r_i < u->n; r_i++, aOcc++)
        {
            l = 0;
            v = (ori == 1?((uint64_t)((u->a[u->n - r_i - 1])^(uint64_t)(0x100000000)))>>32:((uint64_t)(u->a[r_i]))>>32);

            if(p_v != (uint32_t)-1)
            {
                av = asg_arc_a(read_sg, p_v);
                nv = asg_arc_n(read_sg, p_v);
                
                for (k = 0; k < nv; k++)
                {
                    if(av[k].del) continue;
                    if(av[k].v == v) 
                    {
                        l = asg_arc_len(av[k]);
                        break;
                    }
                }
                if(k == nv) fprintf(stderr, "ERROR\n");
            }

            p_v = v; len += l;

            if((v>>1) == (best->x.ul>>33) && xLen_thres == (uint32_t)-1)
            {
                ///cov->pos_idx[v>>1] = len; 
                xR = &(reverse_sources[v>>1]);
                for (j = 0; j < xR->length; j++)
                {
                    h = &(xR->buffer[j]);
                    sq = &(coverage_cut[Get_qn(*h)]);
                    st = &(coverage_cut[Get_tn(*h)]);
                    if(st->del || read_sg->seq[Get_tn(*h)].del) continue;
                    r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, &t);
                    ///if it is a contained overlap, skip
                    if(r < 0) continue;

                    rId = t.v>>1;
                    if(read_sg->seq[rId].del == 1) continue;
                    if(pos_idx[rId] == (uint64_t)-1) continue;
                    w = (uint32_t)(pos_idx[rId]);
                    if(rId != (w>>1)) continue;

                    tmp = get_xy_pos_by_pos(read_sg, &t, v, w, len, pos_idx[w>>1]>>32,
                                                        (uint32_t)-1, targetBaseLen, &(t.el));
                    if(((tmp>>32) == (uint32_t)-1) || (((uint32_t)tmp) == (uint32_t)-1)) continue;
                    if(t.el) continue; ///must

                    t_offset.Off = tmp;
                    t_offset.x = t;
                    t_offset.weight = 1;
                    if(t_offset.Off == best->Off && t_offset.x.v == best->x.v && t_offset.x.ul == best->x.ul)
                    {
                        xLen_thres = Get_xOff(best->Off) + targetBaseLen - Get_yOff(best->Off);
                    }
                }
            }

            if(len >= xLen_thres)
            {
                (*xEnd) = aOcc;
                return;
            }
        }
    }

    (*xEnd) = aOcc;
}


void get_base_boundary_advance_back(R_to_U* ruIndex, ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, 
asg_t *read_g, uint64_t* position_index, int max_hang, int min_ovlp, ma_utg_t *xReads, ma_utg_t *yReads,
uint32_t xUid, uint32_t yUid, long long xBegIndex, long long xEndIndex, long long yBegIndex, long long yEndIndex, 
uint32_t rev, kvec_asg_arc_t_offset* u_buffer, kvec_t_i32_warp* tailIndex, kvec_t_i32_warp* prevIndex,
uint32_t* xBeg, uint32_t* xEnd, uint32_t* yBeg, uint32_t* yEnd)
{
    long long k, j, offset, m;
    ma_hit_t_alloc *xR = NULL;
    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL, *st = NULL;
    int32_t r;
    asg_arc_t t;
    uint32_t rId, Hap_uId, is_Unitig, v, w, v_dir, w_dir;
    uint64_t tmp;
    asg_arc_t_offset t_offset;
    u_buffer->a.n = 0;
    (*xBeg) = (*xEnd) = (*yBeg) = (*yEnd) = (uint32_t)-1;
    for (k = xBegIndex; k <= xEndIndex; k++)
    {
        xR = &(reverse_sources[xReads->a[k]>>33]);
        for (j = 0; j < xR->length; j++)
        {
            h = &(xR->buffer[j]);
            sq = &(coverage_cut[Get_qn(*h)]);
            st = &(coverage_cut[Get_tn(*h)]);
            if(st->del || read_g->seq[Get_tn(*h)].del) continue;

            r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                            asm_opt.max_hang_rate, min_ovlp, &t);
            ///if it is a contained overlap, skip
            if(r < 0) continue;

            rId = t.v>>1;
            if(read_g->seq[rId].del == 1) continue;
            ///there are two cases: 
            ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
            ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
            get_R_to_U(ruIndex, rId, &Hap_uId, &is_Unitig);
            if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
            if(Hap_uId != yUid) continue;

            v = xReads->a[k]>>32;
            get_R_to_U(ruIndex, v>>1, &Hap_uId, &is_Unitig);
            if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
            if(Hap_uId != xUid) continue;
            if((uint32_t)(position_index[v>>1]) != k) continue;

            w = (yReads->a[(uint32_t)(position_index[rId])])>>32;

            v_dir = ((t.ul>>32)==v)?1:0;
            w_dir = (t.v == w)?1:0;
            if(rev == 0 && v_dir != w_dir) continue;
            if(rev == 1 && v_dir == w_dir) continue;
        
            /****************************may have bugs********************************/
            offset = (uint32_t)(position_index[rId]);
            if(offset < yBegIndex || offset > yEndIndex) continue;
            /****************************may have bugs********************************/

            tmp = get_xy_pos(read_g, &t, v, w, xReads->len, yReads->len, position_index, &(t.el));
            if(((tmp>>32) == (uint32_t)-1) || (((uint32_t)tmp) == (uint32_t)-1)) continue;

            t_offset.Off = tmp;
            t_offset.x = t;
            t_offset.weight = 1;
            kv_push(asg_arc_t_offset, u_buffer->a, t_offset);      
        }
    }
    if(u_buffer->a.n == 0) return;

    qsort(u_buffer->a.a, u_buffer->a.n, sizeof(asg_arc_t_offset), cmp_hap_alignment_chaining);

    ///print_asg_arc_t_offset(u_buffer->a.a, u_buffer->a.n, "before");

    for (k = 1, m = 1; k < (long long)u_buffer->a.n; k++)
    {
        if(u_buffer->a.a[m-1].Off == u_buffer->a.a[k].Off)
        {
            u_buffer->a.a[m-1].weight += u_buffer->a.a[k].weight;
            if(u_buffer->a.a[k].x.ol > u_buffer->a.a[m-1].x.ol)
            {
                u_buffer->a.a[m-1].x = u_buffer->a.a[k].x;
            }
            continue;
        }
        u_buffer->a.a[m] = u_buffer->a.a[k];
        m++;
    }
    u_buffer->a.n = m;

    ///print_asg_arc_t_offset(u_buffer->a.a, u_buffer->a.n, "after");
    quick_LIS(u_buffer->a.a, u_buffer->a.n, tailIndex, prevIndex);

    if(tailIndex->a.n == 0) return;

    (*xBeg) = Get_xOff(u_buffer->a.a[tailIndex->a.a[0]].Off);
    (*yBeg) = Get_yOff(u_buffer->a.a[tailIndex->a.a[0]].Off);
    (*xEnd) = Get_xOff(u_buffer->a.a[tailIndex->a.a[tailIndex->a.n-1]].Off);
    (*yEnd) = Get_yOff(u_buffer->a.a[tailIndex->a.a[tailIndex->a.n-1]].Off);
}

uint32_t determine_hap_overlap_type_advance_back(hap_candidates* hap_can, ma_utg_t *xReads, ma_utg_t *yReads,
R_to_U* ruIndex, ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, asg_t *read_g, uint64_t* position_index, 
int max_hang, int min_ovlp, uint32_t xUid, uint32_t yUid, kvec_asg_arc_t_offset* u_buffer, kvec_t_i32_warp* tailIndex, 
kvec_t_i32_warp* prevIndex, long long* r_x_pos_beg, long long* r_x_pos_end, long long* r_y_pos_beg, long long* r_y_pos_end)
{
    uint32_t x_pos_beg, y_pos_beg, x_pos_end, y_pos_end;
    /*************************x***************************/
    get_base_boundary_advance_back(ruIndex, reverse_sources, coverage_cut, read_g, position_index, 
    max_hang, min_ovlp, xReads, yReads, xUid, yUid, Get_x_beg(*hap_can), Get_x_end(*hap_can), 
    Get_y_beg(*hap_can), Get_y_end(*hap_can), Get_rev(*hap_can), u_buffer, tailIndex, prevIndex, 
    &x_pos_beg, &x_pos_end, &y_pos_beg, &y_pos_end);
    /*************************x***************************/

    if(x_pos_beg == (uint32_t)-1 || y_pos_beg == (uint32_t)-1 
                    || x_pos_end == (uint32_t)-1 || y_pos_end == (uint32_t)-1)
    {
        return (uint32_t)-1;
    }
    if(x_pos_beg > x_pos_end || y_pos_beg > y_pos_end) return (uint32_t)-1;
    
    /**
    #define X2Y 0
    #define Y2X 1
    #define XCY 2
    #define YCX 3
    **/
    return classify_hap_overlap(x_pos_beg, x_pos_end, xReads->len, y_pos_beg, y_pos_end, yReads->len, 
    r_x_pos_beg, r_x_pos_end, r_y_pos_beg, r_y_pos_end);
}

void get_idx_by_base(ma_utg_t *x, asg_t *read_g, long long beg_base, long long end_base, 
                                            long long* beg_idx, long long* end_idx)
{
    long long offset, r_beg, r_end;
    uint64_t i, rId;
    (*beg_idx) = (*end_idx) = -1;
    for (i = 0, offset = 0; i < x->n; i++)
    {
        rId = x->a[i]>>33;
        r_beg = offset; r_end = offset + (long long)(read_g->seq[rId].len) - 1;
        offset += (uint32_t)x->a[i];
        if(beg_base > r_end || r_beg > end_base)
        {
            if((*beg_idx) != -1 && (*end_idx) != -1) break;
            continue;
        } 
        if((*beg_idx) == -1) (*beg_idx) = i;
        (*end_idx) = i;
    }
}

int get_base_boundary_chain(R_to_U* ruIndex, ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, 
asg_t *read_g, uint64_t* position_index, int max_hang, int min_ovlp, ma_utg_t *xReads, ma_utg_t *yReads,
uint32_t xUid, uint32_t yUid, long long xBegIndex, long long xEndIndex, long long yBegIndex, long long yEndIndex, 
uint32_t rev, kvec_asg_arc_t_offset* u_buffer, kvec_t_i32_warp* tailIndex, kvec_t_i32_warp* prevIndex)
{
    long long k, j, l, offset, m;
    ma_hit_t_alloc *xR = NULL;
    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL, *st = NULL;
    int32_t r;
    asg_arc_t t;
    uint32_t rId, Hap_uId, is_Unitig, v, w, v_dir, w_dir;
    uint64_t tmp;
    asg_arc_t_offset t_offset;
    u_buffer->a.n = 0;
    for (k = xBegIndex; k <= xEndIndex; k++)
    {
        xR = &(reverse_sources[xReads->a[k]>>33]);
        for (j = 0; j < xR->length; j++)
        {
            h = &(xR->buffer[j]);
            sq = &(coverage_cut[Get_qn(*h)]);
            st = &(coverage_cut[Get_tn(*h)]);
            if(st->del || read_g->seq[Get_tn(*h)].del) continue;

            r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                            asm_opt.max_hang_rate, min_ovlp, &t);
            ///if it is a contained overlap, skip
            if(r < 0) continue;

            rId = t.v>>1;
            if(read_g->seq[rId].del == 1) continue;
            ///there are two cases: 
            ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
            ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
            get_R_to_U(ruIndex, rId, &Hap_uId, &is_Unitig);
            if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
            if(Hap_uId != yUid) continue;

            v = xReads->a[k]>>32;
            get_R_to_U(ruIndex, v>>1, &Hap_uId, &is_Unitig);
            if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
            if(Hap_uId != xUid) continue;
            if((uint32_t)(position_index[v>>1]) != k) continue;

            w = (yReads->a[(uint32_t)(position_index[rId])])>>32;

            v_dir = ((t.ul>>32)==v)?1:0;
            w_dir = (t.v == w)?1:0;
            if(rev == 0 && v_dir != w_dir) continue;
            if(rev == 1 && v_dir == w_dir) continue;
        
            /****************************may have bugs********************************/
            offset = (uint32_t)(position_index[rId]);
            if(offset < yBegIndex || offset > yEndIndex) continue;
            /****************************may have bugs********************************/

            tmp = get_xy_pos(read_g, &t, v, w, xReads->len, yReads->len, position_index, &(t.el));
            if(((tmp>>32) == (uint32_t)-1) || (((uint32_t)tmp) == (uint32_t)-1)) continue;

            t_offset.Off = tmp;
            t_offset.x = t;
            t_offset.weight = 1;
            kv_push(asg_arc_t_offset, u_buffer->a, t_offset);      
        }
    }
    if(u_buffer->a.n == 0) return 0;

    qsort(u_buffer->a.a, u_buffer->a.n, sizeof(asg_arc_t_offset), cmp_hap_alignment_chaining);

    ///print_asg_arc_t_offset(u_buffer->a.a, u_buffer->a.n, "before");
    for (k = 1, l = 0, m = 0; k <= (long long)u_buffer->a.n; ++k)
    {
        if (k == (long long)u_buffer->a.n || u_buffer->a.a[k].Off != u_buffer->a.a[l].Off) 
        {
            u_buffer->a.a[m] = u_buffer->a.a[l];
            for (l += 1; l < k; l++)
            {
                u_buffer->a.a[m].weight += u_buffer->a.a[l].weight;
                if(u_buffer->a.a[l].x.ol > u_buffer->a.a[m].x.ol)
                {
                    u_buffer->a.a[m].x = u_buffer->a.a[l].x;
                }
            }
            l = k;
            m++;
        }
    } 
    u_buffer->a.n = m;

    ///print_asg_arc_t_offset(u_buffer->a.a, u_buffer->a.n, "after");
    quick_LIS(u_buffer->a.a, u_buffer->a.n, tailIndex, prevIndex);

    if(tailIndex->a.n == 0) return 0;
    return 1;
}


void get_base_boundary_advance(R_to_U* ruIndex, ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, 
asg_t *read_g, uint64_t* position_index, int max_hang, int min_ovlp, ma_utg_t *xReads, ma_utg_t *yReads,
uint32_t xUid, uint32_t yUid, long long xBegIndex, long long xEndIndex, long long yBegIndex, long long yEndIndex, 
uint32_t rev, kvec_asg_arc_t_offset* u_buffer, kvec_t_i32_warp* tailIndex, kvec_t_i32_warp* prevIndex,
uint32_t* xBeg, uint32_t* xEnd, uint32_t* yBeg, uint32_t* yEnd)
{
    long long offset;
    long long new_xBeg, new_yBeg, new_xEnd, new_yEnd;
    long long new_xIdxBeg, new_yIdxBeg, new_xIdxEnd, new_yIdxEnd;

    (*xBeg) = (*xEnd) = (*yBeg) = (*yEnd) = (uint32_t)-1;
    if(!get_base_boundary_chain(ruIndex, reverse_sources, coverage_cut, read_g, position_index, 
    max_hang, min_ovlp, xReads, yReads, xUid, yUid, xBegIndex, xEndIndex, yBegIndex, yEndIndex, 
    rev, u_buffer, tailIndex, prevIndex))
    {
        return;
    }
    ///base
    new_xBeg = Get_xOff(u_buffer->a.a[tailIndex->a.a[0]].Off);
    new_yBeg = Get_yOff(u_buffer->a.a[tailIndex->a.a[0]].Off);
    new_xEnd = Get_xOff(u_buffer->a.a[tailIndex->a.a[tailIndex->a.n-1]].Off);
    new_yEnd = Get_yOff(u_buffer->a.a[tailIndex->a.a[tailIndex->a.n-1]].Off);

    if(new_xBeg > new_xEnd || new_yBeg > new_yEnd) return;
    
    classify_hap_overlap(new_xBeg, new_xEnd, xReads->len, new_yBeg, new_yEnd, yReads->len, 
                                                        &new_xBeg, &new_xEnd, &new_yBeg, &new_yEnd);

    if(rev) 
    {
        new_yBeg = yReads->len - new_yBeg - 1;
        new_yEnd = yReads->len - new_yEnd - 1;
        offset = new_yBeg; new_yBeg = new_yEnd; new_yEnd = offset;
    }
    ///idx
    get_idx_by_base(xReads, read_g, new_xBeg, new_xEnd, &new_xIdxBeg, &new_xIdxEnd);
    get_idx_by_base(yReads, read_g, new_yBeg, new_yEnd, &new_yIdxBeg, &new_yIdxEnd);
    if(new_xIdxBeg == -1 || new_xIdxEnd == -1 || new_yIdxBeg == -1 || new_yIdxEnd == -1) return;
    
    if(!get_base_boundary_chain(ruIndex, reverse_sources, coverage_cut, read_g, position_index, 
    max_hang, min_ovlp, xReads, yReads, xUid, yUid, new_xIdxBeg, new_xIdxEnd, new_yIdxBeg, 
    new_yIdxEnd, rev, u_buffer, tailIndex, prevIndex))
    {
        return;
    }

    (*xBeg) = Get_xOff(u_buffer->a.a[tailIndex->a.a[0]].Off);
    (*yBeg) = Get_yOff(u_buffer->a.a[tailIndex->a.a[0]].Off);
    (*xEnd) = Get_xOff(u_buffer->a.a[tailIndex->a.a[tailIndex->a.n-1]].Off);
    (*yEnd) = Get_yOff(u_buffer->a.a[tailIndex->a.a[tailIndex->a.n-1]].Off);
}

#define generic_key(x) (x)
KRADIX_SORT_INIT(i32, int32_t, generic_key, sizeof(int32_t))
KRADIX_SORT_INIT(ru32, uint32_t, generic_key, sizeof(uint32_t))

long long get_chain_score(ma_utg_t *xReads, asg_t *read_g, kvec_asg_arc_t_offset* u_buffer, kvec_t_i32_warp* tailIndex, kvec_t_i32_warp* idx, 
ma_hit_t_alloc* reverse_sources, long long xBegPos, long long xEndPos)
{
    long long offset, r_beg, r_end, inp_beg, inp_end, hap_beg, hap_end, inp_match, hap_match, ovlp;
    uint64_t i, k, rId;
    idx->a.n = 0;
    
    for (i = k = 0; i < tailIndex->a.n; i++)
    {
        rId = u_buffer->a.a[tailIndex->a.a[i]].x.ul>>33;
        
        
        for (; k < xReads->n; k++)
        {
            if(rId == (xReads->a[k]>>33)) break;
        }

        if(k >= xReads->n)
        {
            for (k = 0; k < xReads->n; k++)
            {
                if(rId == (xReads->a[k]>>33)) break;
            }
        }

        if(k < xReads->n) kv_push(int32_t, idx->a, k);
        else
        {
            fprintf(stderr, "\nERROR-get_chain_score: tailIndex->a.n: %lu, xReads->n: %lu\n", (uint64_t)tailIndex->a.n, (uint64_t)xReads->n);
        } 
    }
    

    radix_sort_i32(idx->a.a, idx->a.a + idx->a.n);
    inp_beg = -1; inp_end = -2;
    hap_beg = -1; hap_end = -2;
    for (i = k = 0, offset = 0, inp_match = hap_match = 0; i < xReads->n; i++)
    {
        rId = xReads->a[i]>>33;
        r_beg = offset; r_end = offset + (long long)(read_g->seq[rId].len) - 1;
        offset += (uint32_t)xReads->a[i];

        if(reverse_sources[rId].length > 0)
        {
            if(r_beg <= hap_end)
            {
                hap_end = MAX(hap_end, r_end);
            } 
            else
            {
                ///match += (hap_end - hap_beg + 1);
                ovlp = (long long)(MIN(hap_end, xEndPos)) - (long long)(MAX(hap_beg, xBegPos)) + 1;
                hap_match += (ovlp >= 0? ovlp : 0);
                hap_beg = r_beg; hap_end = r_end;
            }
        }

        for (; k < idx->a.n; k++)
        {
            if(i <= (uint64_t)idx->a.a[k]) break;
        }

        if(k >= idx->a.n) continue;
        
        if(i == (uint64_t)idx->a.a[k])
        {
            if(r_beg <= inp_end)
            {
                inp_end = MAX(inp_end, r_end);
            } 
            else
            {
                ovlp = (long long)(MIN(inp_end, xEndPos)) - (long long)(MAX(inp_beg, xBegPos)) + 1;
                inp_match += (ovlp >= 0? ovlp : 0);
                inp_beg = r_beg; inp_end = r_end;
            }
        }
    }
    
    ovlp = (long long)(MIN(inp_end, xEndPos)) - (long long)(MAX(inp_beg, xBegPos)) + 1;
    inp_match += (ovlp >= 0? ovlp : 0);
    
    ovlp = (long long)(MIN(hap_end, xEndPos)) - (long long)(MAX(hap_beg, xBegPos)) + 1;
    hap_match += (ovlp >= 0? ovlp : 0);

    // if(inp_match > (xEndPos - xBegPos + 1)) fprintf(stderr, "ERROR1\n");
    // if(hap_match > (xEndPos - xBegPos + 1)) fprintf(stderr, "ERRO2\n");
    // if(inp_match > hap_match) fprintf(stderr, "ERROR3\n");
    // fprintf(stderr, "tailIndex->a.n: %u, xReads->n: %u, total_match: %lld, hap_match: %lld, inp_match: %lld\n", 
    //                     tailIndex->a.n, xReads->n, (xEndPos - xBegPos + 1), hap_match, inp_match);

    return ((double)(inp_match)*CHAIN_MATCH) - ((double)(hap_match-inp_match)*CHAIN_UNMATCH);
}


uint32_t determine_hap_overlap_type_advance(hap_candidates* hap_can, ma_utg_t *xReads, ma_utg_t *yReads,
R_to_U* ruIndex, ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, asg_t *read_g, uint64_t* position_index, 
int max_hang, int min_ovlp, uint32_t xUid, uint32_t yUid, kvec_asg_arc_t_offset* u_buffer, kvec_t_i32_warp* tailIndex, 
kvec_t_i32_warp* prevIndex, long long* r_x_pos_beg, long long* r_x_pos_end, long long* r_y_pos_beg, long long* r_y_pos_end)
{
    uint32_t x_pos_beg, y_pos_beg, x_pos_end, y_pos_end;
    /*************************x***************************/
    get_base_boundary_advance(ruIndex, reverse_sources, coverage_cut, read_g, position_index, 
    max_hang, min_ovlp, xReads, yReads, xUid, yUid, Get_x_beg(*hap_can), Get_x_end(*hap_can), 
    Get_y_beg(*hap_can), Get_y_end(*hap_can), Get_rev(*hap_can), u_buffer, tailIndex, prevIndex, 
    &x_pos_beg, &x_pos_end, &y_pos_beg, &y_pos_end);
    /*************************x***************************/
    if(x_pos_beg == (uint32_t)-1 || y_pos_beg == (uint32_t)-1 
                    || x_pos_end == (uint32_t)-1 || y_pos_end == (uint32_t)-1)
    {
        return (uint32_t)-1;
    }
    if(x_pos_beg > x_pos_end || y_pos_beg > y_pos_end) return (uint32_t)-1;

    /**
    #define X2Y 0
    #define Y2X 1
    #define XCY 2
    #define YCX 3
    **/
    hap_can->index_end = classify_hap_overlap(x_pos_beg, x_pos_end, xReads->len, y_pos_beg, y_pos_end, yReads->len, 
    r_x_pos_beg, r_x_pos_end, r_y_pos_beg, r_y_pos_end);
    hap_can->x_beg_pos = MIN((uint32_t)(position_index[u_buffer->a.a[tailIndex->a.a[0]].x.ul>>33]), 
                                (uint32_t)(position_index[u_buffer->a.a[tailIndex->a.a[tailIndex->a.n-1]].x.ul>>33]));
    hap_can->x_end_pos = MAX((uint32_t)(position_index[u_buffer->a.a[tailIndex->a.a[0]].x.ul>>33]), 
                                (uint32_t)(position_index[u_buffer->a.a[tailIndex->a.a[tailIndex->a.n-1]].x.ul>>33]));
    hap_can->y_beg_pos = MIN((uint32_t)(position_index[u_buffer->a.a[tailIndex->a.a[0]].x.v>>1]), 
                                (uint32_t)(position_index[u_buffer->a.a[tailIndex->a.a[tailIndex->a.n-1]].x.v>>1]));
    hap_can->y_end_pos = MAX((uint32_t)(position_index[u_buffer->a.a[tailIndex->a.a[0]].x.v>>1]), 
                                (uint32_t)(position_index[u_buffer->a.a[tailIndex->a.a[tailIndex->a.n-1]].x.v>>1]));
    double xLeftMatch, xLeftTotal;
    get_pair_hap_similarity(xReads->a + hap_can->x_beg_pos, hap_can->x_end_pos + 1 - hap_can->x_beg_pos, 
        yUid, reverse_sources, read_g, ruIndex, &xLeftMatch, &xLeftTotal);
    if(xLeftMatch == 0 || xLeftTotal == 0) return (uint32_t)-1;
    hap_can->weight = xLeftMatch;
    hap_can->index_beg = xLeftTotal;
    hap_can->score = get_chain_score(xReads, read_g, u_buffer, tailIndex, prevIndex, reverse_sources, 
                                                                        (*r_x_pos_beg), (*r_x_pos_end));
    if(hap_can->score <= 0) return (uint32_t)-1;
    return hap_can->index_end;
}



uint32_t determine_hap_overlap_type(hap_candidates* hap_can, ma_utg_t *xReads, ma_utg_t *yReads,
R_to_U* ruIndex, ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, asg_t *read_g, 
uint64_t* position_index, int max_hang, int min_ovlp, uint32_t xUid, uint32_t yUid,
long long* r_x_pos_beg, long long* r_x_pos_end, long long* r_y_pos_beg, long long* r_y_pos_end)
{
    uint32_t x_pos_beg, y_pos_beg, x_pos_end, y_pos_end;
    /*************************x***************************/
    get_base_boundary(ruIndex, reverse_sources, coverage_cut, read_g, position_index, max_hang, 
    min_ovlp, xReads, yReads, xUid, yUid, Get_x_beg(*hap_can), Get_x_end(*hap_can), 
    Get_y_beg(*hap_can), Get_y_end(*hap_can), 0, Get_rev(*hap_can), &x_pos_beg, &y_pos_beg);

    get_base_boundary(ruIndex, reverse_sources, coverage_cut, read_g, position_index, max_hang, 
    min_ovlp, xReads, yReads, xUid, yUid, Get_x_beg(*hap_can), Get_x_end(*hap_can), 
    Get_y_beg(*hap_can), Get_y_end(*hap_can), 1, Get_rev(*hap_can), &x_pos_end, &y_pos_end);
    /*************************x***************************/

    if(x_pos_beg == (uint32_t)-1 || y_pos_beg == (uint32_t)-1 
                    || x_pos_end == (uint32_t)-1 || y_pos_end == (uint32_t)-1)
    {
        return (uint32_t)-1;
    }
    if(x_pos_beg > x_pos_end || y_pos_beg > y_pos_end) return (uint32_t)-1;
    
    /**
    #define X2Y 0
    #define Y2X 1
    #define XCY 2
    #define YCX 3
    **/
    return classify_hap_overlap(x_pos_beg, x_pos_end, xReads->len, y_pos_beg, y_pos_end, yReads->len, 
    r_x_pos_beg, r_x_pos_end, r_y_pos_beg, r_y_pos_end);
}

void adjust_hap_overlaps_score(ma_utg_t* xReads, float *sim, long long *score, 
long long xUid, long long yUid, long long xBeg, long long xEnd)
{
    uint64_t all, found;
    if(count_unique_k_mers(xReads->s + xBeg, xEnd+1-xBeg, xUid, yUid, &all, &found))
    {
        double k_w = 1;
        if(sim) (*sim) = MAX((*sim), (all == 0?0:(((double)found)/((double)all))));
        if(all) k_w += ((double)(found)/(double)(all));
        if(score) (*score) = ((*score)*k_w)/2;
    }
}

uint32_t calculate_pair_hap_similarity_advance(hap_candidates* hap_can, 
uint64_t* position_index, uint32_t xUid, uint32_t yUid, ma_utg_t* xReads, ma_utg_t* yReads,
ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, asg_t *read_g, R_to_U* ruIndex, ma_sub_t *coverage_cut, 
float Hap_rate, int is_local, int max_hang, int min_ovlp, uint64_t cov_threshold, kvec_asg_arc_t_offset* u_buffer, 
kvec_t_i32_warp* tailIndex, kvec_t_i32_warp* prevIndex, hap_cov_t *cov, long long* r_x_pos_beg, long long* r_x_pos_end, 
long long* r_y_pos_beg, long long* r_y_pos_end, float *sim)
{
    uint32_t max_count = 0, min_count = 0, flag;
    uint32_t xLen = xReads->n, xIndex;
    uint32_t yLen = yReads->n, yIndex;
    uint32_t xLeftBeg, xLeftLen, yLeftBeg, yLeftLen;
    uint32_t xRightBeg, xRightLen, yRightBeg, yRightLen;
    double xLeftMatch = 0, xLeftTotal = 0, yLeftMatch = 0, yLeftTotal = 0;
    double xRightMatch = 0, xRightTotal = 0, yRightMatch = 0, yRightTotal = 0;
    asg_arc_t* arch = NULL;
    

    arch = &(hap_can->t);
    xIndex = (uint32_t)(position_index[arch->ul>>33]);
    yIndex = (uint32_t)(position_index[arch->v>>1]);
    
    if(hap_can->rev == 0)
    {
        xLeftBeg = 0; xLeftLen = xIndex; xRightBeg = xIndex; xRightLen = xLen - xRightBeg;
        yLeftBeg = 0; yLeftLen = yIndex; yRightBeg = yIndex; yRightLen = yLen - yRightBeg;
    }
    else
    {
        xLeftBeg = 0; xLeftLen = xIndex; xRightBeg = xIndex; xRightLen = xLen - xRightBeg;

        yLeftBeg = yIndex + 1; yLeftLen = yLen - yLeftBeg;
        yRightBeg = 0; yRightLen = yIndex + 1;
    }

    flag = Get_type(*hap_can);


    if(flag == XCY)
    {
        get_pair_hap_similarity(yReads->a, yLen, xUid, reverse_sources, read_g, ruIndex, 
        &yLeftMatch, &yLeftTotal);
        max_count = yLeftMatch;
        min_count = yLeftTotal;
    }
    else if(flag == YCX)
    {
        get_pair_hap_similarity(xReads->a, xLen, yUid, reverse_sources, read_g, ruIndex, 
        &xLeftMatch, &xLeftTotal);
        max_count = xLeftMatch;
        min_count = xLeftTotal;
    }
    else if(flag == X2Y)
    {
        get_pair_hap_similarity(yReads->a+yLeftBeg, yLeftLen, xUid, reverse_sources, read_g, ruIndex,
        &yLeftMatch, &yLeftTotal);
        get_pair_hap_similarity(xReads->a+xRightBeg, xRightLen, yUid, reverse_sources, read_g, ruIndex, 
        &xRightMatch, &xRightTotal);
        max_count = yLeftMatch + xRightMatch;
        min_count = yLeftTotal + xRightTotal;
    }
    else if(flag == Y2X)
    {
        get_pair_hap_similarity(xReads->a+xLeftBeg, xLeftLen, yUid, reverse_sources, read_g, ruIndex, 
        &xLeftMatch, &xLeftTotal);
        get_pair_hap_similarity(yReads->a+yRightBeg, yRightLen, xUid, reverse_sources, read_g, ruIndex,
        &yRightMatch, &yRightTotal);
        max_count = xLeftMatch + yRightMatch;
        min_count = xLeftTotal + yRightTotal;
    } else abort();

    hap_can->weight = hap_can->index_beg = 0;
    if(min_count == 0) return NON_PLOID;
    if((max_count > min_count*Hap_rate) || is_local)
    {
        long long r_x_interval_beg, r_x_interval_end, r_y_interval_beg, r_y_interval_end;
        uint64_t ploid_coverage = 0;

        ///for containment, don't need to do anything
        get_hap_alignment_boundary(xReads, yReads, flag, xLeftMatch, xLeftTotal, 
        yLeftMatch, yLeftTotal, xRightMatch, xRightTotal, yRightMatch, yRightTotal,
        xLeftBeg, xLeftLen, yLeftBeg, yLeftLen, xRightBeg, xRightLen, yRightBeg, yRightLen,
        xUid, yUid, Hap_rate, is_local, position_index, reverse_sources, read_g, ruIndex, 
        hap_can->rev, &r_x_interval_beg, &r_x_interval_end, &r_y_interval_beg, &r_y_interval_end);

        if(r_x_interval_beg < 0 || r_x_interval_end < 0 || r_y_interval_beg < 0 || r_y_interval_end < 0)
        {
            return NON_PLOID;
        }

        get_pair_hap_similarity(xReads->a + r_x_interval_beg, r_x_interval_end + 1 - r_x_interval_beg, 
        yUid, reverse_sources, read_g, ruIndex, &xLeftMatch, &xLeftTotal);
        if(xLeftMatch == 0 || xLeftTotal == 0 || (is_local == 0 && xLeftMatch <= xLeftTotal*Hap_rate))
        {
            return NON_PLOID;
        }
                
        
        hap_can->weight = xLeftMatch;
        hap_can->index_beg = xLeftTotal;
        hap_can->index_end = flag;
        hap_can->x_beg_pos = r_x_interval_beg;
        hap_can->x_end_pos = r_x_interval_end;
        hap_can->y_beg_pos = r_y_interval_beg;
        hap_can->y_end_pos = r_y_interval_end;
        hap_can->index_end = determine_hap_overlap_type_advance(hap_can, xReads, yReads,
        ruIndex, reverse_sources, coverage_cut, read_g, position_index, max_hang, min_ovlp, 
        xUid, yUid, u_buffer, tailIndex, prevIndex, r_x_pos_beg, r_x_pos_end, r_y_pos_beg, 
        r_y_pos_end);
        if(hap_can->index_end == XCY && yReads->len > (xReads->len*2)) return NON_PLOID;
        if(hap_can->index_end == YCX && xReads->len > (yReads->len*2)) return NON_PLOID;
        if(hap_can->index_end == (uint32_t)-1) return NON_PLOID;

        ploid_coverage = get_pair_purge_coverage(xReads, *r_x_pos_beg, *r_x_pos_end, 
        yReads, *r_y_pos_beg, *r_y_pos_end, hap_can->rev, read_g, cov);

        if(cov_threshold > 0 && ploid_coverage >= cov_threshold) return NON_PLOID;

        get_pair_hap_similarity_by_base(xReads, read_g, yUid, reverse_sources, ruIndex, 
        *r_x_pos_beg, *r_x_pos_end, &xLeftMatch, &xLeftTotal);
        (*sim) = (xLeftTotal== 0? 0:((double)xLeftMatch)/((double)xLeftTotal));

        // adjust_hap_overlaps_score(xReads, sim, &(hap_can->score), xUid, yUid, (*r_x_pos_beg), (*r_x_pos_end));

        if(xLeftMatch == 0 || xLeftTotal == 0 || (*sim) <= Hap_rate)
        {
            return NON_PLOID;
        }

        return PLOID;
    } 
	return NON_PLOID;
}


void print_hap_paf(ma_ug_t *ug, hap_overlaps* ovlp)
{
    fprintf(stderr, "utg%.6d%c\t%u(%u)\t%u(%u)\t%u(%u)\t%c\tutg%.6d%c\t%u(%u)\t%u(%u)\t%u(%u)\t%u\t%u\t%lld(%u)\n", 
    ovlp->xUid+1, "lc"[ug->u.a[ovlp->xUid].circ], ug->u.a[ovlp->xUid].len, ug->u.a[ovlp->xUid].n,
    ovlp->x_beg_pos, ovlp->x_beg_id, ovlp->x_end_pos, ovlp->x_end_id, "+-"[ovlp->rev], 
    ovlp->yUid+1, "lc"[ug->u.a[ovlp->yUid].circ], ug->u.a[ovlp->yUid].len, ug->u.a[ovlp->yUid].n,
    ovlp->y_beg_pos, ovlp->y_beg_id, ovlp->y_end_pos, ovlp->y_end_id, ovlp->type, ovlp->weight, 
    ovlp->score, ovlp->status);
}

inline long long get_max_index(asg_arc_t_offset* x, int32_t* Scores, uint8_t* Flag, long long n,
long long x_readLen, long long y_readLen)
{
    long long i = 0, max_result = -1, max_i = -1, min_xLen = x_readLen * 2 + 2, x_off, y_off, tmp_xLen;
    for (i = 0; i < n; i++)
    {
        if(Flag[i] != 0) continue;
        x_off = (long long)(x[i].Off>>32);
        y_off = (long long)((uint32_t)x[i].Off);
        if(Scores[i] > max_result)
        {
            max_result = Scores[i];
            max_i = i;
            min_xLen = get_hap_overlapLen(x_off, x_off, x_readLen, y_off, y_off, y_readLen, 
            NULL, NULL, NULL, NULL);
        }
        else if(Scores[i] == max_result)
        {
            tmp_xLen = get_hap_overlapLen(x_off, x_off, x_readLen, y_off, y_off, y_readLen, 
            NULL, NULL, NULL, NULL);

            if(tmp_xLen < min_xLen)
            {
                max_result = Scores[i];
                max_i = i;
                min_xLen = tmp_xLen;
            }
        }
    }

    return max_i;
}


inline void get_chain_details(int32_t* Pres, int32_t* Results, uint8_t* Flag, long long max_i,
long long* chainLen, long long* dup)
{
    long long i = max_i;
    (*chainLen) = 0;
    (*dup) = 0;

    while (i >= 0)
    {
        if(Flag[i] == 1) (*dup)++;
        Results[(*chainLen)] = i;
        i = Pres[i];
        (*chainLen)++;
    }
}

inline void push_hap_can(asg_arc_t_offset* x, kvec_hap_candidates* u_can, int32_t* Results, 
long long chainLen, long long x_readLen, long long y_readLen)
{
    if(chainLen <= 0) return;
    hap_candidates hap_can;
    hap_can.rev = x[Results[0]].x.el;
    hap_can.x_beg_pos = hap_can.x_end_pos = (uint32_t)(x[Results[0]].Off>>32);
    hap_can.y_beg_pos = hap_can.y_end_pos = (uint32_t)(x[Results[0]].Off);
    hap_can.weight = 0;
    long long i = 0;
    uint64_t totalWeigth = 0;
    ///fprintf(stderr, "^^^chainLen: %lld\n", chainLen);
    for (i = 0; i < chainLen; i++)
    {
        ///fprintf(stderr, "i: %lld, Results[i]: %d\n", i, Results[i]);
        hap_can.x_beg_pos = (uint32_t)(x[Results[i]].Off>>32);
        hap_can.y_beg_pos = (uint32_t)(x[Results[i]].Off);
        hap_can.weight += x[Results[i]].weight;
    }

    for (i = 0; i < chainLen; i++)
    {
        totalWeigth += x[Results[i]].weight;
        if(totalWeigth >= (hap_can.weight/2)) break;
    }

    if(i >= chainLen) i = chainLen-1;
    ///Get_total(hap_can) = Results[i];
    hap_can.t = x[Results[i]].x;
    Get_type(hap_can) = classify_hap_overlap(hap_can.x_beg_pos, hap_can.x_end_pos, 
    x_readLen, hap_can.y_beg_pos, hap_can.y_end_pos, y_readLen, NULL, NULL, NULL, NULL);
    kv_push(hap_candidates, u_can->a, hap_can);
    if(hap_can.x_beg_pos > hap_can.x_end_pos || hap_can.y_beg_pos > hap_can.y_end_pos)
    {
        fprintf(stderr, "ERROR\n");
    }
}



void print_chain_data(int32_t* Scores, int32_t* Pres, int32_t* Begs, long long n)
{
    fprintf(stderr,"*****\nn_chain: %lld\n", n);
    long long i;
    for (i = 0; i < n; i++)
    {
        fprintf(stderr, "i: %lld, Scores: %d, Pres: %d, Begs: %d\n", 
        i, Scores[i], Pres[i], Begs[i]);
    }
}



void hap_chaining(asg_arc_t_offset* x, uint32_t n, kvec_t_i32_warp* score_vc, kvec_t_i32_warp* prevIndex_vec, 
kvec_t_i32_warp* begIndex_vec, kvec_t_u8_warp* flag_vec, float band_width_threshold, long long max_skip, 
long long x_readLen, long long y_readLen, kvec_hap_candidates* u_can)
{
    #define DUP_OVLP_RATE 0.75
    score_vc->a.n = prevIndex_vec->a.n = begIndex_vec->a.n = flag_vec->a.n = 0;
    if(n == 0) return;
    kv_resize(int32_t, score_vc->a, n);
    kv_resize(int32_t, prevIndex_vec->a, n);
    kv_resize(int32_t, begIndex_vec->a, n);
    kv_resize(uint8_t, flag_vec->a, n);
    
    int32_t* Scores = score_vc->a.a;
    int32_t* Pres = prevIndex_vec->a.a;
    int32_t* Begs = begIndex_vec->a.a;
    uint8_t* Flag = flag_vec->a.a;
    long long i, j, n_max_skip, x_off, y_off, max_beg, max_j = -1, max_score, score;
    long long distance_x, distance_y, total_distance_x, total_distance_y, distance_gap;
    float gap_rate, band_width_penalty = 1 / band_width_threshold;
    long long max_result, max_i, min_xLen, tmp_xLen, chainLen = 0, dup = 0;
    max_result = max_i = -1; min_xLen = x_readLen * 2 + 2;
    for (i = 0; i < n; i++)
    {
		n_max_skip = 0;

        x_off = (long long)(x[i].Off>>32);
        y_off = (long long)((uint32_t)x[i].Off);
        max_j = -1;
        max_score = x[i].weight;
        max_beg = i; //i itself
        ///may have a pre-cut condition for j
        for (j = i - 1; j >= 0; --j) 
        {
            distance_x = x_off - (long long)(x[j].Off>>32);
            distance_y = y_off - (long long)((uint32_t)x[j].Off);
            ///x has been sorted by x_off
            if(distance_x <= 0 || distance_y <= 0) continue;

            total_distance_x = x_off - (long long)(x[Begs[j]].Off>>32);
            total_distance_y = y_off - (long long)((uint32_t)x[Begs[j]].Off);

            distance_gap = total_distance_x - total_distance_y;
            if(distance_gap < 0) distance_gap = -distance_gap;
            if(distance_gap > band_width_threshold * total_distance_x)
            {
                continue;
            }

            score = x[i].weight;
            gap_rate = (float)((float)(distance_gap)/(float)(total_distance_x));
            score -= (long long)(score * gap_rate * band_width_penalty);
            score += Scores[j];

            ///find a new max score
			if (score > max_score) {
				max_score = score;
				max_j = j;
				max_beg = Begs[j];
				n_max_skip = 0;
			} 
            else 
            {
				if (++n_max_skip > max_skip) break;
			}
        }

        Scores[i] = max_score;
        Pres[i] = max_j;
        Begs[i] = max_beg;

        if(Scores[i] > max_result)
        {
            max_result = Scores[i];
            max_i = i;
            min_xLen = get_hap_overlapLen(x_off, x_off, x_readLen, y_off, y_off, y_readLen, 
            NULL, NULL, NULL, NULL);
        }
        else if(Scores[i] == max_result)
        {
            tmp_xLen = get_hap_overlapLen(x_off, x_off, x_readLen, y_off, y_off, y_readLen, 
            NULL, NULL, NULL, NULL);

            if(tmp_xLen < min_xLen)
            {
                max_result = Scores[i];
                max_i = i;
                min_xLen = tmp_xLen;
            }
        }
        Flag[i] = 0;
    }

    // print_asg_arc_t_offset(x, n);
    // print_chain_data(Scores, Pres, Begs, n);
    while (max_i != -1)
    {
        get_chain_details(Pres, Begs, Flag, max_i, &chainLen, &dup);
        if(chainLen == 0) break;
        if(dup > chainLen*DUP_OVLP_RATE)
        {
            for (i = 0; i < chainLen; i++)
            {
                if(Flag[Begs[i]] == 1) continue;
                Flag[Begs[i]] = 2;
            }
            
        }
        else
        {
            push_hap_can(x, u_can, Begs, chainLen, x_readLen, y_readLen);

            for (i = 0; i < chainLen; i++)
            {
                Flag[Begs[i]] = 1;
            }
        }

        max_i = get_max_index(x, Scores, Flag, n, x_readLen, y_readLen);
    }
}

void get_candidate_hap_alignment(kvec_hap_candidates* u_can, kvec_asg_arc_t_offset* u_buffer,
kvec_t_i32_warp* score_vc, kvec_t_i32_warp* prevIndex_vec, kvec_t_i32_warp* begIndex_vec,
kvec_t_u8_warp* flag_vec, float band_width_threshold, long long max_skip, long long x_readLen, 
long long y_readLen)
{
    u_can->a.n = 0;
    if(u_buffer->a.n == 0) return;
    uint32_t i = 0, /**anchor_i = 0, **/m = 1, break_point = (uint32_t)-1, is_merge;

    qsort(u_buffer->a.a, u_buffer->a.n, sizeof(asg_arc_t_offset), cmp_hap_alignment_chaining);

    ///print_asg_arc_t_offset(u_buffer->a.a, u_buffer->a.n);

    for (i = 1; i < u_buffer->a.n; i++)
    {
        is_merge = 0;
        if(u_buffer->a.a[m-1].x.el == u_buffer->a.a[i].x.el)
        {
            if(u_buffer->a.a[m-1].Off == u_buffer->a.a[i].Off) is_merge = 1;
            ///I think we don't need the following merging
            // if(is_merge == 0 && (Get_xOff(u_buffer->a.a[m-1].Off)==Get_xOff(u_buffer->a.a[i].Off)))
            // {
            //     if((Get_yOff(u_buffer->a.a[i].Off)-(Get_yOff(u_buffer->a.a[m-1].Off))) == (i-anchor_i))///not sure why, does it use for tolerate indels in overlaps? 
            //     {
            //         is_merge = 1;
            //     }
            // }

            if(is_merge)
            {
                u_buffer->a.a[m-1].weight += u_buffer->a.a[i].weight;
                continue;
            }
        } 
        u_buffer->a.a[m] = u_buffer->a.a[i];
        // anchor_i = i;
        if(u_buffer->a.a[m].x.el != u_buffer->a.a[m-1].x.el) break_point = m;
        m++;
    }
    u_buffer->a.n = m;
    if(break_point > u_buffer->a.n) break_point = u_buffer->a.n;

    ///print_asg_arc_t_offset(u_buffer->a.a, u_buffer->a.n);

    hap_chaining(u_buffer->a.a, break_point, score_vc, prevIndex_vec, begIndex_vec, flag_vec,
    band_width_threshold, max_skip, x_readLen, y_readLen, u_can);
    
    hap_chaining(u_buffer->a.a + break_point, u_buffer->a.n - break_point, score_vc, prevIndex_vec, 
    begIndex_vec, flag_vec, band_width_threshold, max_skip, x_readLen, y_readLen, u_can);
}

int filter_secondary_chain(long long max_score, long long cur_score, double rate)
{
    if(cur_score >= max_score) return 1;
    long long diff = max_score - cur_score;
    if(max_score < 0) max_score *= -1;
    if(diff >= max_score*(1-rate)) return 0;
    return 1;
}


void filter_secondary_ovlp(kvec_hap_overlaps *x, kvec_t_u64_warp *a, float sim_flt, float ovlp_flt)
{
    if(sim_flt == 0 || ovlp_flt == 0 || x->a.n == 0) return;
    #define f_ovlp(s_0, e_0, s_1, e_1) ((MIN((e_0), (e_1)) > MAX((s_0), (s_1)))? MIN((e_0), (e_1)) - MAX((s_0), (s_1)):0)
    uint32_t i, m, k;
    uint64_t t, ovlp;
    hap_overlaps *p = NULL;
    a->a.n = 0;
    for (i = 0; i < x->a.n; i++)
    {
        if(x->a.a[i].s < sim_flt) continue;
        t = x->a.a[i].x_beg_pos; t<<=32; t |= x->a.a[i].x_end_pos;
        kv_push(uint64_t, a->a, t);
    }

    if(a->a.n == 0) return;
    ks_introsort_uint64_t(a->a.n, a->a.a);

    for (i = m = 1; i < a->a.n; ++i) 
    {
        t = a->a.a[m-1];
        ovlp = f_ovlp(t>>32, (uint32_t)t, a->a.a[i]>>32, (uint32_t)a->a.a[i]);
        if(ovlp == 0)
        {
            a->a.a[m] = a->a.a[i];
            m++;
        } 
        else
        {
            t = MIN(a->a.a[m-1]>>32, a->a.a[i]>>32);
            t<<=32;
            t |= MAX((uint32_t)a->a.a[m-1], (uint32_t)a->a.a[i]);
            a->a.a[m-1] = t;
        }
    }
    a->a.n = m;

    for (i = m = 0; i < x->a.n; i++)
    {
        p = &(x->a.a[i]);
        if(p->s < sim_flt)
        {
            for (k = ovlp = 0; k < a->a.n; k++)
            {
                ovlp += f_ovlp(p->x_beg_pos, p->x_end_pos, a->a.a[k]>>32, (uint32_t)a->a.a[k]);
                if(ovlp >= ovlp_flt*(p->x_end_pos-p->x_beg_pos)) break;
            }
            if(k < a->a.n) continue;
            if(ovlp >= ovlp_flt*(p->x_end_pos-p->x_beg_pos)) continue;
        } 
        x->a.a[m] = x->a.a[i];
        m++;
    }
    x->a.n = m;
}

static void hap_alignment_advance_worker(void *_data, long eid, int tid)
{
    hap_alignment_struct_pip* hap_buf = (hap_alignment_struct_pip*)_data;
    ma_ug_t *ug = hap_buf->ug;
    asg_t *read_g = hap_buf->read_g;
    ma_hit_t_alloc* sources = hap_buf->sources;
    ma_hit_t_alloc* reverse_sources = hap_buf->reverse_sources;
    R_to_U* ruIndex = hap_buf->ruIndex;
    ma_sub_t *coverage_cut = hap_buf->coverage_cut;
    uint64_t* position_index = hap_buf->position_index;
    float Hap_rate = hap_buf->Hap_rate/**MIN(hap_buf->Hap_rate, 0.2)**/, sim;
    int max_hang = hap_buf->max_hang;
    int min_ovlp = hap_buf->min_ovlp;
    float chain_rate = hap_buf->chain_rate;
    hap_overlaps_list* all_ovlp = hap_buf->all_ovlp;
    uint32_t Input_uId = eid;
    uint64_t* vote_counting = hap_buf->buf[tid].vote_counting;
    uint8_t* visit = hap_buf->buf[tid].visit;
    kvec_t_u64_warp* u_vecs = &(hap_buf->buf[tid].u_vecs);
    kvec_asg_arc_t_offset* u_buffer = &(hap_buf->buf[tid].u_buffer);
    kvec_hap_candidates* u_can = &(hap_buf->buf[tid].u_can);
    kvec_t_i32_warp* score_vc = &(hap_buf->buf[tid].u_buffer_tailIndex);
    kvec_t_i32_warp* prevIndex_vec = &(hap_buf->buf[tid].u_buffer_prevIndex);
    kvec_t_i32_warp* begIndex_vec = &(hap_buf->buf[tid].u_buffer_beg);
    kvec_t_u8_warp* flag_vec = &(hap_buf->buf[tid].u_buffer_flag);
    uint64_t cov_threshold = hap_buf->cov_threshold;
    hap_cov_t *cov = hap_buf->cov;
    uint8_t *hh = hap_buf->hh, hhc;
    if(hap_buf->cov_threshold < 0) cov_threshold = (uint64_t)-1;
    ma_utg_t *xReads = NULL, *yReads = NULL;
    ma_hit_t_alloc *xR = NULL;
    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL, *st = NULL;
    asg_t* nsg = ug->g;
    uint32_t i, j, v, rId, k, is_Unitig, Hap_uId, xUid, yUid, seedOcc;
    uint64_t tmp, max_weight, m;
    long long r_x_pos_beg, r_x_pos_end, r_y_pos_beg, r_y_pos_end, max_score;
    int32_t r;
    asg_arc_t t;
    asg_arc_t_offset t_offset;
    hap_overlaps hap_align;
    hap_overlaps *hap_align_x = NULL;
    xUid = Input_uId;
    if(nsg->seq[xUid].del || nsg->seq[xUid].c == ALTER_LABLE) return;
    memset(vote_counting, 0, sizeof(uint64_t)*nsg->n_seq);
    memset(visit, 0, nsg->n_seq);
    u_vecs->a.n = 0;
    u_can->a.n = 0;

    xReads = &(ug->u.a[xUid]);
    for (i = 0; i < xReads->n; i++)
    {
        xR = &(reverse_sources[xReads->a[i]>>33]);
        
        for (k = 0; k < xR->length; k++)
        {
            rId = Get_tn(xR->buffer[k]);
        
            if(read_g->seq[rId].del == 1)
            {
                ///get the id of read that contains it 
                get_R_to_U(ruIndex, rId, &rId, &is_Unitig);
                if(rId == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[rId].del == 1) continue;
            }

            ///there are two cases: 
            ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
            ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
            get_R_to_U(ruIndex, rId, &Hap_uId, &is_Unitig);
            if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
            ///here rId is the id of the read coming from the different haplotype
            ///Hap_cId is the id of the corresponding contig (note here is the contig, instead of untig)
            if(visit[Hap_uId]!=0) continue; ///one read only has one vote for one hap unitig
            visit[Hap_uId] = 1;
            if(vote_counting[Hap_uId] < UINT64_MAX) vote_counting[Hap_uId]++;
        }

        clean_visit_flag(visit, read_g, ruIndex, nsg->n_seq, xR);
    }



    u_vecs->a.n = 0;
    for (i = 0; i < nsg->n_seq; i++)
    {
        if(i == xUid) continue;
        if(vote_counting[i] == 0) continue;
        tmp = vote_counting[i]; tmp = tmp << 32; tmp = tmp | (uint64_t)i;
        kv_push(uint64_t, u_vecs->a, tmp);
    }

    if(u_vecs->a.n == 0) return;
    sort_kvec_t_u64_warp(u_vecs, 1);


    ///scan each candidate unitig
    for (i = 0; i < u_vecs->a.n; i++)
    {
        yUid = (uint32_t)u_vecs->a.a[i];
        seedOcc = u_vecs->a.a[i]>>32;
        xReads = &(ug->u.a[xUid]);
        yReads = &(ug->u.a[yUid]);
        u_buffer->a.n = 0;

        for (k = 0; k < xReads->n; k++)
        {
            xR = &(reverse_sources[xReads->a[k]>>33]);
            hhc = hh[xReads->a[k]>>33];
            for (j = 0; j < xR->length; j++)
            {
                h = &(xR->buffer[j]);
                sq = &(coverage_cut[Get_qn(*h)]);
                st = &(coverage_cut[Get_tn(*h)]);
                if(st->del || read_g->seq[Get_tn(*h)].del) continue;

                r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, &t);
                ///if it is a contained overlap, skip
                if(r < 0) continue;

                rId = t.v>>1;
                if(read_g->seq[rId].del == 1) continue;
                ///there are two cases: 
                ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
                ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
                get_R_to_U(ruIndex, rId, &Hap_uId, &is_Unitig);
                if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
                if(Hap_uId != yUid) continue;

                v = xReads->a[k]>>32;
                get_R_to_U(ruIndex, v>>1, &Hap_uId, &is_Unitig);
                if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
                if(Hap_uId != xUid) continue;
                if((uint32_t)(position_index[v>>1]) != k) continue;

                if(asm_opt.purge_level_primary <= 2 &&
                    (prefilter((uint32_t)(position_index[v>>1]), (uint32_t)(position_index[rId]), 
                    xReads->n, yReads->n, 0, Hap_rate, seedOcc)==NON_PLOID) && 
                   (prefilter((uint32_t)(position_index[v>>1]), (uint32_t)(position_index[rId]), 
                    xReads->n, yReads->n, 1, Hap_rate, seedOcc)==NON_PLOID))
                {
                    continue;
                }
            
                t_offset.Off = get_xy_pos(read_g, &t, v, (yReads->a[(uint32_t)(position_index[rId])])>>32, 
                xReads->len, yReads->len, position_index, &(t.el));
                if(((t_offset.Off>>32) == (uint32_t)-1) || (((uint32_t)t_offset.Off) == (uint32_t)-1)) continue;

                t_offset.x = t;
                t_offset.weight = hhc;

                kv_push(asg_arc_t_offset, u_buffer->a, t_offset);
            }

            deduplicate_edge(u_buffer);
        }

        if(u_buffer->a.n == 0) continue;

        
        get_candidate_hap_alignment(u_can, u_buffer, score_vc, prevIndex_vec, begIndex_vec,
        flag_vec, chain_rate, 50, xReads->len, yReads->len);

        if(u_can->a.n == 0) continue;
        
        qsort(u_can->a.a, u_can->a.n, sizeof(hap_candidates), cmp_hap_candidates);

        memset(&hap_align, 0, sizeof(hap_overlaps));    
        m = all_ovlp->x[xUid].a.n;     
        max_weight = 0; max_score = 0;
        for (k = 0; k < u_can->a.n; k++)
        {
            if(u_can->a.a[k].weight < max_weight*0.33) continue;
            if(calculate_pair_hap_similarity_advance(&(u_can->a.a[k]), position_index, xUid, yUid, 
            xReads, yReads, sources, reverse_sources, read_g, ruIndex, coverage_cut, Hap_rate,
            (asm_opt.purge_level_primary<=2? 0:1), max_hang, min_ovlp, cov_threshold, u_buffer, 
            score_vc, prevIndex_vec, cov, &r_x_pos_beg, &r_x_pos_end, &r_y_pos_beg, &r_y_pos_end, &sim)!=PLOID)
            {
                continue;
            }
            ///max_weight == 0 means the first matched chain
            if(max_weight == 0 || max_score < u_can->a.a[k].score) max_score = u_can->a.a[k].score;
            if(max_weight < u_can->a.a[k].weight) max_weight = u_can->a.a[k].weight;
            ///if one is positive and another one is negative, it is wrong
            if(!filter_secondary_chain(max_score, u_can->a.a[k].score, CHAIN_FILTER_RATE)) continue;
                            
            hap_align.rev = Get_rev(u_can->a.a[k]);
            hap_align.type = Get_type(u_can->a.a[k]);
            hap_align.x_beg_id = Get_x_beg(u_can->a.a[k]);
            hap_align.x_end_id = Get_x_end(u_can->a.a[k]) + 1;
            hap_align.y_beg_id = Get_y_beg(u_can->a.a[k]);
            hap_align.y_end_id = Get_y_end(u_can->a.a[k]) + 1;
            hap_align.weight = Get_match(u_can->a.a[k]);
            hap_align.score = u_can->a.a[k].score;
            hap_align.x_beg_pos = r_x_pos_beg;
            hap_align.x_end_pos = r_x_pos_end + 1;
            if(hap_align.rev == 0)
            {
                hap_align.y_beg_pos = r_y_pos_beg;
                hap_align.y_end_pos = r_y_pos_end + 1;
            }
            else
            {
                hap_align.y_beg_pos = yReads->len - r_y_pos_end - 1;
                hap_align.y_end_pos = yReads->len - r_y_pos_beg - 1 + 1;
            }
            hap_align.xUid = xUid;
            hap_align.yUid = yUid;
            hap_align.status = SELF_EXIST;
            hap_align.s = sim;
            kv_push(hap_overlaps, all_ovlp->x[hap_align.xUid].a, hap_align);
        }
        /**
        ///chains with same xUid && yUid
        for (k = m; k < all_ovlp->x[xUid].a.n; k++)
        {
            if(!filter_secondary_chain(max_score, 
                        all_ovlp->x[xUid].a.a[k].score, CHAIN_FILTER_RATE))
            {
                continue;
            } 
            all_ovlp->x[xUid].a.a[m] = all_ovlp->x[xUid].a.a[k];
            m++;
        }
        all_ovlp->x[xUid].a.n = m;
        **/

        hap_align_x = NULL;
        for (k = m; k < all_ovlp->x[xUid].a.n; k++)
        {
            if(all_ovlp->x[xUid].a.a[k].score != max_score) continue;
            if(hap_align_x == NULL || all_ovlp->x[xUid].a.a[k].weight > hap_align_x->weight)
            {
                hap_align_x = &(all_ovlp->x[xUid].a.a[k]);
            }
            else if(all_ovlp->x[xUid].a.a[k].weight == hap_align_x->weight)
            {
                if((all_ovlp->x[xUid].a.a[k].x_end_pos 
                            - all_ovlp->x[xUid].a.a[k].x_beg_pos) < 
                    (hap_align_x->x_end_pos - hap_align_x->x_beg_pos))
                {
                    hap_align_x = &(all_ovlp->x[xUid].a.a[k]);
                }
            }
        }
        if(hap_align_x)
        {
            all_ovlp->x[xUid].a.a[m] = (*hap_align_x);
            all_ovlp->x[xUid].a.n = m + 1;
        }
    }
    // filter_secondary_ovlp(&all_ovlp->x[xUid], u_vecs, hap_buf->Hap_rate, 0.7);
}

int get_specific_hap_overlap(kvec_hap_overlaps* x, uint32_t qn, uint32_t tn)
{
    uint32_t i;
    for (i = 0; i < x->a.n; i++)
    {
        if(x->a.a[i].xUid == qn && x->a.a[i].yUid == tn)
        {
            return i;
        }        
    }

    return -1;
}

void set_reverse_hap_overlap(hap_overlaps* dest, hap_overlaps* source, uint32_t* types)
{
    dest->status = REVE_EXIST;
    dest->rev = source->rev;
    dest->type = types[source->type];
    dest->weight = source->weight;
    dest->xUid = source->yUid;
    dest->yUid = source->xUid;
    dest->x_beg_pos = source->y_beg_pos;
    dest->x_end_pos = source->y_end_pos;
    dest->y_beg_pos = source->x_beg_pos;
    dest->y_end_pos = source->x_end_pos;
    dest->x_beg_id = source->y_beg_id;
    dest->x_end_id = source->y_end_id;
    dest->y_beg_id = source->x_beg_id;
    dest->y_end_id = source->x_end_id;
    dest->score = source->score;
}

/**
#define X2Y 0
#define Y2X 1
#define XCY 2
#define YCX 3
**/
void normalize_hap_overlaps(hap_overlaps_list* all_ovlp, hap_overlaps_list* back_all_ovlp)
{
    hap_overlaps *x = NULL, *y = NULL;
    uint32_t v, i, uId, qn, tn;
    uint32_t types[4];
    types[X2Y] = Y2X; types[Y2X] = X2Y; types[XCY] = YCX; types[YCX] = XCY;
    int index;
    for (v = 0; v < all_ovlp->num; v++)
    {
        uId = v;
        for (i = 0; i < all_ovlp->x[uId].a.n; i++)
        {
            qn = all_ovlp->x[uId].a.a[i].xUid;
            tn = all_ovlp->x[uId].a.a[i].yUid;
            x = &(all_ovlp->x[uId].a.a[i]);
            index = get_specific_hap_overlap(&(all_ovlp->x[tn]), tn, qn);
            if(index != -1)
            {
                y = &(all_ovlp->x[tn].a.a[index]);
                if(x->rev == y->rev && types[x->type]==y->type) continue;
                if(x->weight >= y->weight)
                {
                    kv_push(hap_overlaps, back_all_ovlp->x[tn].a, (*y));
                    set_reverse_hap_overlap(y, x, types);
                }
                else
                {
                    kv_push(hap_overlaps, back_all_ovlp->x[qn].a, (*x));
                    set_reverse_hap_overlap(x, y, types);
                }
            }
            else
            {
                kv_pushp(hap_overlaps, all_ovlp->x[tn].a, &y);
                set_reverse_hap_overlap(y, x, types);
            }
        }
    }
}

inline uint64_t calculate_bi_weight(hap_overlaps *x, ma_ug_t *ug, asg_t *read_g, ma_hit_t_alloc* reverse_sources, 
R_to_U* ruIndex)
{
    double xMatch, xTotal;
    uint64_t weight = x->weight;
    ma_utg_t *yReads = &(ug->u.a[x->yUid]);
    get_pair_hap_similarity(yReads->a + x->y_beg_id, x->y_end_id - x->y_beg_id, 
    x->xUid, reverse_sources, read_g, ruIndex, &xMatch, &xTotal);
    weight += xMatch;
    return weight;
}

void normalize_hap_overlaps_advance(hap_overlaps_list* all_ovlp, hap_overlaps_list* back_all_ovlp,
ma_ug_t *ug, asg_t *read_g, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex)
{
    hap_overlaps *x = NULL, *y = NULL;
    uint32_t v, i, uId, qn, tn;
    uint32_t types[4];
    
    types[X2Y] = Y2X; types[Y2X] = X2Y; types[XCY] = YCX; types[YCX] = XCY;
    int index;
    for (v = 0; v < all_ovlp->num; v++)
    {
        uId = v;
        for (i = 0; i < all_ovlp->x[uId].a.n; i++)
        {
            qn = all_ovlp->x[uId].a.a[i].xUid;
            tn = all_ovlp->x[uId].a.a[i].yUid;
            x = &(all_ovlp->x[uId].a.a[i]);
            index = get_specific_hap_overlap(&(all_ovlp->x[tn]), tn, qn);
            if(index != -1)
            {
                y = &(all_ovlp->x[tn].a.a[index]);
                if(x->rev == y->rev && types[x->type]==y->type) continue;
                ///if(x->weight >= y->weight)
                // if((calculate_bi_weight(x, ug, read_g, reverse_sources, ruIndex)) >= 
                //    (calculate_bi_weight(y, ug, read_g, reverse_sources, ruIndex)))
                if(x->score >= y->score)
                {
                    kv_push(hap_overlaps, back_all_ovlp->x[tn].a, (*y));
                    set_reverse_hap_overlap(y, x, types);
                }
                else
                {
                    kv_push(hap_overlaps, back_all_ovlp->x[qn].a, (*x));
                    set_reverse_hap_overlap(x, y, types);
                }
            }
            else
            {
                kv_pushp(hap_overlaps, all_ovlp->x[tn].a, &y);
                set_reverse_hap_overlap(y, x, types);
            }
        }
    }
}

void get_p_nodes(p_g_t *pg, p_node_t **x, uint32_t *x_occ, uint32_t id)
{
    if(x) (*x) = pg->pg_het_node.a + pg->pg_h_lev_idx.a[id].beg;
    if(x_occ) (*x_occ) = pg->pg_h_lev_idx.a[id].occ;
}

void normalize_hap_overlaps_advance_by_p_g_t(hap_overlaps_list* all_ovlp, hap_overlaps_list* back_all_ovlp,
ma_ug_t *ug, asg_t *read_g, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex, p_g_t *pg, hap_cov_t *cov, 
double filter_rate)
{
    hap_overlaps *x = NULL, *y = NULL;
    uint32_t v, i, uId, qn, tn;
    uint32_t types[4];
    

    types[X2Y] = Y2X; types[Y2X] = X2Y; types[XCY] = YCX; types[YCX] = XCY;
    int index;
    uint32_t k, qs, qe, ts, te, occ, as, ae, ovlp, hetLen, homLen;
    p_node_t *a = NULL;


    for (v = 0; v < all_ovlp->num; v++)
    {
        uId = v;
        for (i = 0; i < all_ovlp->x[uId].a.n; i++)
        {
            /*****************qn*****************/
            qn = all_ovlp->x[uId].a.a[i].xUid;
            qs = all_ovlp->x[uId].a.a[i].x_beg_pos;
            qe = all_ovlp->x[uId].a.a[i].x_end_pos - 1;
            get_p_nodes(pg, &a, &occ, qn);
            for (k = 0, hetLen = 0, homLen = 0; k < occ; k++)
            {
                as = a[k].baseBeg;
                ae = a[k].baseEnd;
                ovlp = ((MIN(qe, ae) >= MAX(qs, as))? MIN(qe, ae) - MAX(qs, as) + 1 : 0);
                if(homLen + hetLen > 0 && ovlp == 0) break;
                if(ovlp == 0) continue;
                if(a[k].h_status == N_HET)
                {
                    homLen += ovlp;
                } 
                else if(asm_opt.polyploidy <= 2 && (a[k].h_status&P_HET))///if(asm_opt.polyploidy <= 2 && (a[k].h_status&S_HET))
                {
                    homLen += ovlp;
                }
                else
                {
                    hetLen += ovlp;
                }                
            }

            if(hetLen <= ((hetLen + homLen) * filter_rate))
            {
                all_ovlp->x[uId].a.a[i].status = DELETE;
                continue;
            }
            /*****************qn*****************/


            /*****************tn*****************/
            tn = all_ovlp->x[uId].a.a[i].yUid;
            ts = all_ovlp->x[uId].a.a[i].y_beg_pos;
            te = all_ovlp->x[uId].a.a[i].y_end_pos - 1;
            get_p_nodes(pg, &a, &occ, tn);
            for (k = 0, hetLen = 0, homLen = 0; k < occ; k++)
            {
                as = a[k].baseBeg;
                ae = a[k].baseEnd;
                ovlp = ((MIN(te, ae) >= MAX(ts, as))? MIN(te, ae) - MAX(ts, as) + 1 : 0);
                if(homLen + hetLen > 0 && ovlp == 0) break;
                if(ovlp == 0) continue;
                if(a[k].h_status == N_HET)
                {
                    homLen += ovlp;
                } 
                else if(asm_opt.polyploidy <= 2 && (a[k].h_status&P_HET))///if(asm_opt.polyploidy <= 2 && (a[k].h_status&S_HET))
                {
                    homLen += ovlp;
                }
                else
                {
                    hetLen += ovlp;
                }                 
            }

            if(hetLen <= ((hetLen + homLen) * filter_rate))
            {
                all_ovlp->x[uId].a.a[i].status = DELETE;
                continue;
            }
            /*****************tn*****************/
        }
    }

    for (v = 0; v < all_ovlp->num; v++)
    {
        uId = v;
        k = 0;
        for (i = 0; i < all_ovlp->x[uId].a.n; i++)
        {
            if(all_ovlp->x[uId].a.a[i].status == DELETE) continue;
            all_ovlp->x[uId].a.a[k] = all_ovlp->x[uId].a.a[i];
            k++;
        }
        all_ovlp->x[uId].a.n = k;
    }


    for (v = 0; v < all_ovlp->num; v++)
    {
        uId = v;
        for (i = 0; i < all_ovlp->x[uId].a.n; i++)
        {
            qn = all_ovlp->x[uId].a.a[i].xUid;
            tn = all_ovlp->x[uId].a.a[i].yUid;
            x = &(all_ovlp->x[uId].a.a[i]);
            index = get_specific_hap_overlap(&(all_ovlp->x[tn]), tn, qn);
            if(index != -1)
            {
                y = &(all_ovlp->x[tn].a.a[index]);
                if(x->rev == y->rev && types[x->type]==y->type) continue;
                if(x->score >= y->score)
                {
                    kv_push(hap_overlaps, back_all_ovlp->x[tn].a, (*y));
                    set_reverse_hap_overlap(y, x, types);
                }
                else
                {
                    kv_push(hap_overlaps, back_all_ovlp->x[qn].a, (*x));
                    set_reverse_hap_overlap(x, y, types);
                }
            }
            else
            {
                kv_pushp(hap_overlaps, all_ovlp->x[tn].a, &y);
                set_reverse_hap_overlap(y, x, types);
            }
        }
    }
}

void filter_hap_overlaps_by_length(hap_overlaps_list* all_ovlp, uint32_t minLen)
{
    if(minLen == 0) return;
    hap_overlaps *x = NULL;
    uint32_t v, i, m, uId;    

    for (v = 0; v < all_ovlp->num; v++)
    {
        uId = v;
        m = 0;
        for (i = 0; i < all_ovlp->x[uId].a.n; i++)
        {
            x = &(all_ovlp->x[uId].a.a[i]);
            if(x->x_end_id - x->x_beg_id < minLen) continue;
            all_ovlp->x[uId].a.a[m] = (*x);
            m++;
        }
        all_ovlp->x[uId].a.n = m;
    }
}

void debug_hap_overlaps(hap_overlaps_list* all_ovlp, hap_overlaps_list* back_all_ovlp)
{
    hap_overlaps *x = NULL, *y = NULL;
    uint32_t v, i, uId, qn, tn;
    uint32_t types[4];
    types[X2Y] = Y2X; types[Y2X] = X2Y; types[XCY] = YCX; types[YCX] = XCY;
    int index;
    for (v = 0; v < all_ovlp->num; v++)
    {
        uId = v;
        for (i = 0; i < all_ovlp->x[uId].a.n; i++)
        {
            qn = all_ovlp->x[uId].a.a[i].xUid;
            tn = all_ovlp->x[uId].a.a[i].yUid;
            x = &(all_ovlp->x[uId].a.a[i]);
            index = get_specific_hap_overlap(&(all_ovlp->x[tn]), tn, qn);
            if(index == -1)
            {
                fprintf(stderr, "ERROR 0\n");
                continue;
            }

            y = &(all_ovlp->x[tn].a.a[index]);
            if(x->rev != y->rev || types[x->type] != y->type)
            {
                fprintf(stderr, "ERROR 1\n");
                continue;
            }

            if(x->status == REVE_EXIST && y->status != SELF_EXIST)
            {
                fprintf(stderr, "ERROR 2\n");
                continue;
            }

            if(x->status == REVE_EXIST)
            {
                if(x->weight != y->weight) fprintf(stderr, "ERROR 3\n"); 
                if(x->xUid != y->yUid) fprintf(stderr, "ERROR 4\n"); 
                if(x->yUid != y->xUid) fprintf(stderr, "ERROR 5\n"); 
                if(x->x_beg_pos != y->y_beg_pos) fprintf(stderr, "ERROR 6\n"); 
                if(x->x_end_pos != y->y_end_pos) fprintf(stderr, "ERROR 7\n"); 
                if(x->y_beg_pos != y->x_beg_pos) fprintf(stderr, "ERROR 8\n");
                if(x->y_end_pos != y->x_end_pos) fprintf(stderr, "ERROR 9\n");
                if(x->x_beg_id != y->y_beg_id) fprintf(stderr, "ERROR 10\n");
                if(x->x_end_id != y->y_end_id) fprintf(stderr, "ERROR 11\n");
                if(x->y_beg_id != y->x_beg_id) fprintf(stderr, "ERROR 12\n");
                if(x->y_beg_id != y->x_beg_id) fprintf(stderr, "ERROR 13\n");
                if(x->y_end_id != y->x_end_id) fprintf(stderr, "ERROR 14\n");

                index = get_specific_hap_overlap(&(back_all_ovlp->x[qn]), qn, tn);
                if(index != -1)
                {
                    if(back_all_ovlp->x[qn].a.a[index].weight > x->weight) fprintf(stderr, "ERROR 15\n");
                }
            }
           
        }
    }
}

void print_purge_gfa(ma_ug_t *ug, asg_t *purge_g)
{
    uint32_t v, i, n_vtx = purge_g->n_seq * 2;
    for (v = 0; v < n_vtx; v++)
    {
        if(v%2==0) fprintf(stderr, "\n");
        if(purge_g->seq[v>>1].del)
        {
            fprintf(stderr, "(D) v>>1: %u, v&1: %u, utg%.6d%c\n", v>>1, v&1, (v>>1)+1,
            "lc"[ug->u.a[v>>1].circ]);
            continue;
        } 

        fprintf(stderr, "(E) v>>1: %u, v&1: %u, utg%.6dl%c\n", v>>1, v&1, (v>>1)+1,
        "lc"[ug->u.a[v>>1].circ]);
        
        uint32_t nv = asg_arc_n(purge_g, v);
		asg_arc_t *av = asg_arc_a(purge_g, v);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            fprintf(stderr, "av[i].ul: %u (utg%.6d%c, dir: %u, len: %u), av[i].v: %u (utg%.6d%c, dir: %u, len: %u), ol: %u\n", 
                (uint32_t)(av[i].ul>>33), (uint32_t)(av[i].ul>>33)+1, "lc"[ug->u.a[av[i].ul>>33].circ], (uint32_t)(av[i].ul>>32)&1, ug->u.a[av[i].ul>>33].len, 
                av[i].v>>1, (av[i].v>>1)+1, "lc"[ug->u.a[av[i].v>>1].circ], av[i].v&1, ug->u.a[av[i].v>>1].len, av[i].ol);
        }

    }
    
}

long long decode_score(uint32_t h_bits, uint32_t l_bits)
{
    uint64_t x; 
    x = h_bits; x <<= 32; x += l_bits;
    long long score = ((uint64_t)((uint64_t)x<<1)>>1);
    if((x>>63) == 0) score *= -1;
    return score;
}

void encode_score(long long i_s, uint32_t *h_bits, uint32_t *l_bits)
{
    uint64_t score = (i_s >= 0? (i_s) : (i_s*(-1)));
    if(i_s >= 0) score += (((uint64_t)1)<<63);
    (*l_bits) = (uint32_t)score; (*h_bits)= (score>>32);
}

uint64_t asg_bub_pop1_purge_graph(asg_t *g, uint32_t v0, int max_dist, buf_t *b)
{   
	uint32_t i, n_pending = 0, n_tips, tip_end;
	uint64_t n_pop = 0;
	///if this node has been deleted
	if (g->seq[v0>>1].del || g->seq[v0>>1].c == ALTER_LABLE) return 0; // already deleted
	///if ((uint32_t)g->idx[v0] < 2) return 0; // no bubbles
    if(get_real_length(g, v0, NULL)<2) return 0;
	///S saves nodes with all incoming edges visited
	b->S.n = b->T.n = b->b.n = b->e.n = 0;
	///for each node, b->a saves all related information
	b->a[v0].c = b->a[v0].d = b->a[v0].m = b->a[v0].nc = b->a[v0].np = 0;
	///b->S is the nodes with all incoming edges visited
	kv_push(uint32_t, b->S, v0);
    n_tips = 0;
    tip_end = (uint32_t)-1;

	do {
		///v is a node that all incoming edges have been visited
		///d is the distance from v0 to v
		uint32_t v = kv_pop(b->S), d = b->a[v].d;
		uint32_t nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
        long long t_s = decode_score(b->a[v].c, b->a[v].m), c_s;
		///why we have this assert?
		///assert(nv > 0);
		///all out-edges of v
		for (i = 0; i < nv; ++i) { // loop through v's neighbors
            ///if this edge has been deleted
			if (av[i].del) continue;
			uint32_t w = av[i].v; // v->w with length l
			binfo_t *t = &b->a[w];
			///that means there is a circle, directly terminate the whole bubble poping
			///if (w == v0) goto pop_reset;
            if ((w>>1) == (v0>>1)) goto pop_reset;
            c_s = decode_score((uint32_t)av[i].ul, av[i].ol);
			///push the edge
            ///high 32-bit of g->idx[v] is the start point of v's edges
            //so here is the point of this specfic edge
			kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);
			///find a too far path? directly terminate the whole bubble poping
			if (d + 1 > (uint32_t)max_dist) break; // too far

            ///if this node
			if (t->s == 0) { // this vertex has never been visited
				kv_push(uint32_t, b->b, w); // save it for revert
				///t->p is the parent node of 
				///t->s = 1 means w has been visited
				///d is len(v0->v), l is len(v->w), so t->d is len(v0->w)
				t->p = v, t->s = 1, t->d = d + 1;
                encode_score(t_s + c_s, &(t->c), &(t->m));
				///incoming edges of w
				///t->r = count_out(g, w^1);
                t->r = get_real_length(g, w^1, NULL);
				++n_pending;
			} else { // visited before
                if((t_s + c_s)> decode_score(t->c, t->m))
                {
                    t->p = v;
                    encode_score(t_s + c_s, &(t->c), &(t->m));
                }
                ///it is the shortest edge
				if (d + 1 < t->d) t->d = d + 1; // update dist
			}
			///assert(t->r > 0);
			//if all incoming edges of w have visited
			//push it to b->S
			if (--(t->r) == 0) {
                uint32_t x = get_real_length(g, w, NULL);
                /****************************may have bugs for bubble********************************/
                if(x > 0)
                {
                    kv_push(uint32_t, b->S, w);
                }
                else
                {
                    ///at most one tip
                    if(n_tips != 0) goto pop_reset;
                    n_tips++;
                    tip_end = w;
                }
                /****************************may have bugs for bubble********************************/
				--n_pending;
			}
		}
        //if found a tip
        /****************************may have bugs for bubble********************************/
        if(n_tips == 1)
        {
            if(tip_end != (uint32_t)-1 && n_pending == 0 && b->S.n == 0)
            {
                kv_push(uint32_t, b->S, tip_end);
                break;
            }
            else
            {
                goto pop_reset;
            }
        }
        /****************************may have bugs for bubble********************************/
		///if i < nv, that means (d + l > max_dist)
		if (i < nv || b->S.n == 0) goto pop_reset;
	} while (b->S.n > 1 || n_pending);

	asg_bub_backtrack_primary(g, v0, b);

    n_pop = 1;
pop_reset:
	for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
		binfo_t *t = &b->a[b->b.a[i]];
		t->s = t->c = t->d = t->m = t->nc = t->np = 0;
	}
	return n_pop;
}



// pop bubbles
int asg_pop_bubble_purge_graph(asg_t *purge_g)
{
	uint32_t v, n_vtx = purge_g->n_seq * 2;
	uint64_t n_pop = 0;
	buf_t b;
	if (!purge_g->is_symm) asg_symm(purge_g);
	memset(&b, 0, sizeof(buf_t));
	///set information for each node
	b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
	//traverse all node with two directions 
	for (v = 0; v < n_vtx; ++v) {
		uint32_t i, n_arc = 0, nv = asg_arc_n(purge_g, v);
		asg_arc_t *av = asg_arc_a(purge_g, v);
		///some node could be deleted
		if (nv < 2 || purge_g->seq[v>>1].del || purge_g->seq[v>>1].c == ALTER_LABLE) continue;
		///some edges could be deleted
		for (i = 0; i < nv; ++i) // asg_bub_pop1() may delete some edges/arcs
			if (!av[i].del) ++n_arc;
		if (n_arc > 1)
			n_pop += asg_bub_pop1_purge_graph(purge_g, v, purge_g->n_seq, &b);            
	}
	free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
	if (n_pop) asg_cleanup(purge_g);

	if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] popped %lu bubbles\n", __func__, (unsigned long)n_pop);
    }
    return n_pop;
}

int get_hap_arch(hap_overlaps* hap, uint32_t qLen, uint32_t tLen, int max_hang, float max_hang_rate, 
int min_ovlp, asg_arc_t* t)
{
    int r;
    ma_hit_t h;
    h.qns = hap->xUid;
    h.qns = h.qns << 32;
    h.qns = h.qns | hap->x_beg_pos;
    h.qe = hap->x_end_pos;
    h.tn = hap->yUid;
    h.ts = hap->y_beg_pos;
    h.te = hap->y_end_pos;
    h.rev = hap->rev;
    h.del = 0;
    h.bl = h.el = h.ml = h.no_l_indel = 0;

    r = ma_hit2arc(&h, qLen, tLen, max_hang, max_hang_rate, min_ovlp, t);
    if(r < 0) return r;
    uint64_t score = (hap->score >= 0? (hap->score) : (hap->score*(-1)));
    if(hap->score >= 0) score += (((uint64_t)1)<<63);
    t->ol = (uint32_t)score;
    t->ul >>= 32; t->ul <<= 32; t->ul |= (score>>32);
    return r;
}

typedef struct {
    uint64_t eid;
    uint64_t score;
}e_score;

typedef struct {
    size_t n, m;
    e_score* a;
}e_score_warp;

#define e_score_key(a) ((a).score)
KRADIX_SORT_INIT(e_score, e_score, e_score_key, member_size(e_score, score))

int purge_g_arc_del_short_diploid_by_score(asg_t *g, float drop_ratio)
{
    e_score_warp b;
    kv_init(b);
    e_score *p = NULL;

	uint32_t v, n_vtx = g->n_seq * 2;
    long long n_cut = 0;
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(g->seq[v>>1].c == ALTER_LABLE || g->seq[v>>1].del) continue;
        asg_arc_t *av = asg_arc_a(g, v);
        uint32_t nv = asg_arc_n(g, v);
        if (nv < 2) continue;
        uint64_t i;
        for (i = 0; i < nv; ++i)
        {
            kv_pushp(e_score, b, &p);
            p->eid = av - g->arc + i;
            p->score = (uint32_t)av[i].ul; 
            p->score <<= 32;
            p->score |= av[i].ol;
        }
	}

    radix_sort_e_score(b.a, b.a + b.n);

    uint64_t k;
    for (k = 0; k < b.n; k++)
    {
        asg_arc_t *a = &g->arc[b.a[k].eid];
		///v is self id, w is the id of another end
		uint32_t i, v = (a->ul)>>32;
		uint32_t nv = asg_arc_n(g, v), kv;
		long long ovlp_max = 0, ovlp;
		asg_arc_t *av = NULL;
		///nv must be >= 2
		if (nv <= 1) continue;
        av = asg_arc_a(g, v);

        ///calculate the longest edge for v and w
		for (i = 0, kv = 0; i < nv; ++i) {
			if (av[i].del) continue;
            ovlp = decode_score((uint32_t)av[i].ul, av[i].ol);
            if (kv == 0 || ovlp_max < ovlp) ovlp_max = ovlp;
			++kv;
		}

        if (kv <= 1) continue;
        ovlp = decode_score((uint32_t)a->ul, a->ol);
		if (kv >= 2)
        {   
            if(ovlp >= 0 && ovlp_max >= 0 && ovlp > ovlp_max * drop_ratio) continue;
        }
         
        a->del = 1;
        asg_arc_del(g, a->v^1, av->ul>>32^1, 1);
        ++n_cut; 
    }

    kv_destroy(b); 
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}

    
	return n_cut;
}


void clean_purge_graph(asg_t *purge_g, float drop_ratio, uint32_t is_force_break)
{
    uint64_t operation = 1;
    while (operation > 0)
    {
        operation = 0;
        operation += asg_pop_bubble_purge_graph(purge_g);
        operation += purge_g_arc_del_short_diploid_by_score(purge_g, drop_ratio);
    }
    
    if(is_force_break) purge_g_arc_del_short_diploid_by_score(purge_g, 1);
}


void get_node_boundary_advance(R_to_U* ruIndex, ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, 
asg_t *read_g, uint64_t* position_index, int max_hang, int min_ovlp, ma_utg_t *xReads, ma_utg_t *yReads,
uint32_t xUid, uint32_t yUid, long long xBegIndex, long long xEndIndex, long long yBegIndex, 
long long yEndIndex, uint32_t dir, uint32_t rev, kvec_asg_arc_t_offset* u_buffer, 
kvec_t_i32_warp* tailIndex, kvec_t_i32_warp* prevIndex, asg_arc_t* reture_t_f, asg_arc_t* reture_t_r)
{
    long long k, j, offset, m;
    ma_hit_t_alloc *xR = NULL;
    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL, *st = NULL;
    int r, index;
    asg_arc_t t_f, t_r;
    uint32_t rId, Hap_uId, is_Unitig, v, w, v_dir, w_dir;
    uint64_t tmp;
    asg_arc_t_offset t_offset;
    reture_t_f->del = reture_t_r->del = 1;
    u_buffer->a.n = 0;


    for (k = xBegIndex; k <= xEndIndex; k++)
    {
        xR = &(reverse_sources[xReads->a[k]>>33]);
        for (j = 0; j < xR->length; j++)
        {
            h = &(xR->buffer[j]);
            sq = &(coverage_cut[Get_qn(*h)]);
            st = &(coverage_cut[Get_tn(*h)]);
            if(st->del || read_g->seq[Get_tn(*h)].del) continue;

            r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                            asm_opt.max_hang_rate, min_ovlp, &t_f);
            ///if it is a contained overlap, skip
            if(r < 0) continue;

            rId = t_f.v>>1;
            if(read_g->seq[rId].del == 1) continue;
            ///there are two cases: 
            ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
            ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
            get_R_to_U(ruIndex, rId, &Hap_uId, &is_Unitig);
            if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
            if(Hap_uId != yUid) continue;

            v = xReads->a[k]>>32;
            get_R_to_U(ruIndex, v>>1, &Hap_uId, &is_Unitig);
            if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
            if(Hap_uId != xUid) continue;
            if((uint32_t)(position_index[v>>1]) != k) continue;

            w = (yReads->a[(uint32_t)(position_index[rId])])>>32;

            v_dir = ((t_f.ul>>32)==v)?1:0;
            w_dir = (t_f.v == w)?1:0;
            if(rev == 0 && v_dir != w_dir) continue;
            if(rev == 1 && v_dir == w_dir) continue;
            if(dir == v_dir) continue;

            /****************************may have bugs********************************/
            offset = (uint32_t)(position_index[rId]);
            if(offset < yBegIndex || offset > yEndIndex) continue;
            /****************************may have bugs********************************/
        
            /************************get reverse edge*************************/
            index = get_specific_overlap(&(reverse_sources[Get_tn(*h)]), Get_tn(*h), Get_qn(*h));
            if(index == -1) continue;
            h = &(reverse_sources[Get_tn(*h)].buffer[index]);
            sq = &(coverage_cut[Get_qn(*h)]);
            st = &(coverage_cut[Get_tn(*h)]);
            if(st->del || read_g->seq[Get_tn(*h)].del) continue;
            r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                            asm_opt.max_hang_rate, min_ovlp, &t_r);
            if(r < 0) continue;
            /************************get reverse edge*************************/

            tmp = get_xy_pos(read_g, &t_f, v, w, xReads->len, yReads->len, position_index, &(t_f.el));
            if(((tmp>>32) == (uint32_t)-1) || (((uint32_t)tmp) == (uint32_t)-1)) continue;

            t_offset.Off = tmp;
            t_offset.x = t_f;
            t_offset.weight = 1;
            kv_push(asg_arc_t_offset, u_buffer->a, t_offset);

        }
    }

    if(u_buffer->a.n == 0) return;

    qsort(u_buffer->a.a, u_buffer->a.n, sizeof(asg_arc_t_offset), cmp_hap_alignment_chaining);

    for (k = 1, m = 1; k < (long long)u_buffer->a.n; k++)
    {
        if(u_buffer->a.a[m-1].Off == u_buffer->a.a[k].Off)
        {
            u_buffer->a.a[m-1].weight += u_buffer->a.a[k].weight;
            if(u_buffer->a.a[k].x.ol > u_buffer->a.a[m-1].x.ol)
            {
                u_buffer->a.a[m-1].x = u_buffer->a.a[k].x;
            }
            continue;
        }
        u_buffer->a.a[m] = u_buffer->a.a[k];
        m++;
    }
    u_buffer->a.n = m;

    quick_LIS(u_buffer->a.a, u_buffer->a.n, tailIndex, prevIndex);

    if(tailIndex->a.n == 0) return;

    if(dir == 0)
    {
        for (k = 0; k < (long long)tailIndex->a.n; k++)
        {
            v = u_buffer->a.a[tailIndex->a.a[k]].x.v>>1;
            w = u_buffer->a.a[tailIndex->a.a[k]].x.ul>>33;
            index = get_specific_overlap(&(reverse_sources[v]), v, w);
            if(index == -1) continue;
            h = &(reverse_sources[v].buffer[index]);
            sq = &(coverage_cut[Get_qn(*h)]);
            st = &(coverage_cut[Get_tn(*h)]);
            if(st->del || read_g->seq[Get_tn(*h)].del) continue;
            r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                            asm_opt.max_hang_rate, min_ovlp, &t_r);
            if(r < 0) continue;

            (*reture_t_f) = u_buffer->a.a[tailIndex->a.a[k]].x;
            (*reture_t_r) = t_r;
            return;
        }
    }
    else
    {
        for (k = tailIndex->a.n-1; k >= 0; k--)
        {
            v = u_buffer->a.a[tailIndex->a.a[k]].x.v>>1;
            w = u_buffer->a.a[tailIndex->a.a[k]].x.ul>>33;
            index = get_specific_overlap(&(reverse_sources[v]), v, w);
            if(index == -1) continue;
            h = &(reverse_sources[v].buffer[index]);
            sq = &(coverage_cut[Get_qn(*h)]);
            st = &(coverage_cut[Get_tn(*h)]);
            if(st->del || read_g->seq[Get_tn(*h)].del) continue;
            r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                            asm_opt.max_hang_rate, min_ovlp, &t_r);
            if(r < 0) continue;

            (*reture_t_f) = u_buffer->a.a[tailIndex->a.a[k]].x;
            (*reture_t_r) = t_r;
            return;
        }
    }
}


void fill_unitig(uint64_t* buffer, uint32_t bufferLen, asg_t* read_g, kvec_asg_arc_t_warp* edge,
uint32_t is_circle, uint64_t* rLen)
{
    uint32_t i, k, totalLen, v, w, nv, l;
    asg_arc_t *av  = NULL;
    (*rLen) = totalLen = 0;
    for (i = 0; i < bufferLen - 1; i++)
    {
        v = (uint64_t)(buffer[i])>>32;
        w = (uint64_t)(buffer[i + 1])>>32; 
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

        if(k == nv)
        {
            for (k = 0; k < edge->a.n; k++)
            {
                if(edge->a.a[k].del) continue;
                if((edge->a.a[k].ul>>32) == v && edge->a.a[k].v == w)
                {
                    l = asg_arc_len(edge->a.a[k]);
                    break;
                }
            }

            if(k == edge->a.n)
            {
                fprintf(stderr, "####ERROR1-fill: i: %u, v>>1: %u, v&1: %u, w>>1: %u, w&1: %u\n",
                i, v>>1, v&1, w>>1, w&1);
            }
        } 


        
        buffer[i] = v; buffer[i] = buffer[i]<<32; buffer[i] = buffer[i] | (uint64_t)(l);
        totalLen += l;
    }
    
    if(i < bufferLen)
    {
        if(is_circle)
        {
            v = (uint64_t)(buffer[i])>>32;
            w = (uint64_t)(buffer[0])>>32; 
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

            if(k == nv)
            {
                for (k = 0; k < edge->a.n; k++)
                {
                    if(edge->a.a[k].del) continue;
                    if((edge->a.a[k].ul>>32) == v && edge->a.a[k].v == w)
                    {
                        l = asg_arc_len(edge->a.a[k]);
                        break;
                    }
                }

                if(k == edge->a.n)
                {
                    fprintf(stderr, "####ERROR2-fill: i: %u, v>>1: %u, v&1: %u, w>>1: %u, w&1: %u\n",
                    i, v>>1, v&1, w>>1, w&1);
                }
            } 

            buffer[i] = v; buffer[i] = buffer[i]<<32; buffer[i] = buffer[i] | (uint64_t)(l);
            totalLen += l;
        }
        else
        {
            v = (uint64_t)(buffer[i])>>32;
            l = read_g->seq[v>>1].len;
            buffer[i] = v;
            buffer[i] = buffer[i]<<32;
            buffer[i] = buffer[i] | (uint64_t)(l);
            totalLen += l;
        }
    }

    (*rLen) = totalLen;

}

void collect_trans_purge_cov(hap_cov_t *cov, ma_ug_t *ug, hap_overlaps* x, uint32_t is_keep_X)
{
    if(ug->u.a[x->xUid].n == 0 || ug->u.a[x->yUid].n == 0) return;
    uint64_t *pri = NULL, pri_n, *aux = NULL, aux_n, i, rId, uCov = 0, uLen = 0;

    if(is_keep_X)
    {
        pri = ug->u.a[x->xUid].a + x->x_beg_id;
        pri_n = x->x_end_id - x->x_beg_id;

        aux = ug->u.a[x->yUid].a + x->y_beg_id;
        aux_n = x->y_end_id - x->y_beg_id;
    }
    else
    {
        pri = ug->u.a[x->yUid].a + x->y_beg_id;
        pri_n = x->y_end_id - x->y_beg_id;

        aux = ug->u.a[x->xUid].a + x->x_beg_id;
        aux_n = x->x_end_id - x->x_beg_id;
    }


    uCov = uLen = 0;
    for (i = 0; i < aux_n; i++)
    {
        rId = aux[i]>>33;
        uCov += cov->cov[rId];
    }

    for (i = 0; i < pri_n; i++)
    {
        rId = pri[i]>>33;
        uLen += cov->read_g->seq[rId].len;
    }
    
    uCov = (uLen == 0? 0 : uCov / uLen);


    for (i = 0; i < pri_n; i++)
    {
        rId = pri[i]>>33;
        cov->cov[rId] += (uCov * cov->read_g->seq[rId].len);
    }
}


void collect_trans_purge_joint_cov(hap_cov_t *cov, ma_ug_t *ug, hap_overlaps* x)
{
    if(ug->u.a[x->xUid].n == 0 || ug->u.a[x->yUid].n == 0) return;
    uint64_t *a[2], a_n[2], uCov[2], uLen[2], uDepth[2], i, rId;
    
    a[0] = ug->u.a[x->xUid].a + x->x_beg_id;
    a_n[0] = x->x_end_id - x->x_beg_id;
    uCov[0] = uLen[0] = 0;
    for (i = 0; i < a_n[0]; i++)
    {
        rId = a[0][i]>>33;
        uCov[0] += cov->cov[rId];
        uLen[0] += cov->read_g->seq[rId].len;
    }
    
    a[1] = ug->u.a[x->yUid].a + x->y_beg_id;
    a_n[1] = x->y_end_id - x->y_beg_id;
    uCov[1] = uLen[1] = 0;
    for (i = 0; i < a_n[1]; i++)
    {
        rId = a[1][i]>>33;
        uCov[1] += cov->cov[rId];
        uLen[1] += cov->read_g->seq[rId].len;
    }
    
    uDepth[0] = (uLen[0] == 0? 0 : uCov[1] / uLen[0]);
    uDepth[1] = (uLen[1] == 0? 0 : uCov[0] / uLen[1]);

    for (i = 0; i < a_n[0]; i++)
    {
        rId = a[0][i]>>33;
        cov->cov[rId] += (uDepth[0] * cov->read_g->seq[rId].len);
    }

    for (i = 0; i < a_n[1]; i++)
    {
        rId = a[1][i]>>33;
        cov->cov[rId] += (uDepth[1] * cov->read_g->seq[rId].len);
    }
}



void purge_merge(asg_t *purge_g, ma_ug_t *ug, hap_overlaps_list* all_ovlp, buf_t* b_0,
R_to_U* ruIndex, ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, asg_t *read_g, 
uint64_t* position_index, kvec_asg_arc_t_offset* u_buffer, kvec_t_i32_warp* tailIndex, 
kvec_t_i32_warp* prevIndex, int max_hang, int min_ovlp, kvec_asg_arc_t_warp* edge, uint8_t* visit,
hap_cov_t *cov)
{
    uint32_t i, nv, k, v, w, x_beg_index, x_end_index, y_beg_index, y_end_index, cut_beg, cut_end, begIndex, endIndex, keepUid;
    hap_overlaps *x = NULL/**, *y = NULL**/;
    ma_utg_t *xReads = NULL, *yReads = NULL;
    asg_arc_t t_forward, t_backward;
    asg_arc_t *av = NULL;
    kvec_t(uint64_t) buffer;
    uint64_t totalLen;
    int index = 0;
    i = 0;
    while (i < b_0->b.n)
    {
        cut_beg = 0; cut_end = (uint32_t)-1;
        kv_init(buffer);
        /********************for the first node********************/
        v = b_0->b.a[i];
        keepUid = v>>1;
        xReads = &(ug->u.a[v>>1]);
        if(v&1)
        {
            for (k = 0; k < xReads->n; k++)
            {
                kv_push(uint64_t, buffer, (xReads->a[xReads->n - k - 1])^(uint64_t)(0x100000000));
            }
        }
        else
        {
            for (k = 0; k < xReads->n; k++)
            {
                kv_push(uint64_t, buffer, xReads->a[k]);
            }
        }
        cut_beg = 0; cut_end = xReads->n - 1;
        i++;
        /********************for the first node********************/

        
        for (; i < b_0->b.n; i++)
        {
            ///x = y = NULL;
            x = NULL;

            v = b_0->b.a[i-1];
            w = b_0->b.a[i];
            

            index = get_specific_hap_overlap(&(all_ovlp->x[v>>1]), v>>1, w>>1);
            x = &(all_ovlp->x[v>>1].a.a[index]);


        
            xReads = &(ug->u.a[v>>1]);
            yReads = &(ug->u.a[w>>1]);

            begIndex = x->x_beg_id;
            if(cut_beg > begIndex) begIndex = cut_beg;

            endIndex = x->x_end_id-1;
            if(cut_end < endIndex) endIndex = cut_end;

            get_node_boundary_advance(ruIndex, reverse_sources, coverage_cut, read_g, position_index, max_hang, 
            min_ovlp, xReads, yReads, v>>1, w>>1, begIndex, endIndex, x->y_beg_id, x->y_end_id-1, v&1, 
            x->rev, u_buffer, tailIndex, prevIndex, &t_forward, &t_backward);

            if(t_forward.del || t_backward.del) break;
            
            kv_push(asg_arc_t, edge->a, t_forward);
            kv_push(asg_arc_t, edge->a, t_backward);

            x_beg_index = 0; x_end_index = xReads->n - 1;
            y_beg_index = 0; y_end_index = yReads->n - 1;
            
            if((v&1) == 0)
            {
                x_end_index = (uint32_t)position_index[t_forward.ul>>33];
                buffer.n = buffer.n - (cut_end - x_end_index);
            }
            else
            {
                x_beg_index = (uint32_t)position_index[t_forward.ul>>33];
                buffer.n = buffer.n - (x_beg_index - cut_beg);
            }

            if((w&1) == 1)
            {
                y_end_index = (uint32_t)position_index[t_forward.v>>1];
            }
            else
            {
                y_beg_index = (uint32_t)position_index[t_forward.v>>1];
            }

            cut_beg = y_beg_index;
            cut_end = y_end_index;

            if((w&1) == 1)
            {
                for (k = y_end_index; k >= y_beg_index; k--)
                {
                    kv_push(uint64_t, buffer, (yReads->a[k])^(uint64_t)(0x100000000));
                    if(k==0) break;
                }
            }
            else
            {
                for (k = y_beg_index; k <= y_end_index; k++)
                {
                    kv_push(uint64_t, buffer, yReads->a[k]);
                }
            }

            purge_g->seq[w>>1].c = ALTER_LABLE;
            if(cov) collect_trans_purge_joint_cov(cov, ug, x);

            // if(buffer.n > 1)
            // {
            //     for (k = 0; k < buffer.n - 1; k++)
            //     {
            //         if((buffer.a[k]>>32) == 854769 && (buffer.a[k+1]>>32) == 64486)
            //         {
            //             fprintf(stderr, "+++++++v: %u, w: %u, xReads->n: %u, yReads->n: %u\n", 
            //                         v, w, xReads->n, yReads->n);
            //             fprintf(stderr, "x->rev: %u, x->x_beg_id: %u, x->x_end_id: %u, x->y_beg_id: %u, x->y_end_id: %u\n", 
            //                         x->rev, x->x_beg_id, x->x_end_id, x->y_beg_id, x->y_end_id);
            //             fprintf(stderr, "t_forward.ul>>32: %u, t_forward.v: %u, y_beg_index: %u, y_end_index: %u\n", 
            //                         t_forward.ul>>32, t_forward.v, y_beg_index, y_end_index);
            //             fprintf(stderr, "type: %u, x->x_beg_pos: %u, x->x_end_pos: %u, xReads->len: %u\n", 
            //                         x->type, x->x_beg_pos, x->x_end_pos, xReads->len);
            //             fprintf(stderr, "x->y_beg_pos: %u, x->y_end_pos: %u, yReads->len: %u\n", 
            //                         x->y_beg_pos, x->y_end_pos, yReads->len);
            //         }
            //     }
            // }
        }

        // fprintf(stderr, "+keepUid: %u, i: %u, b_0->b.n: %u, buffer.n: %u\n", 
        //                             keepUid, i, (uint32_t)b_0->b.n, (uint32_t)buffer.n);
        fill_unitig(buffer.a, buffer.n, read_g, edge, 0, &totalLen);
        ///fprintf(stderr, "-keepUid: %u\n", keepUid);

        xReads = &(ug->u.a[keepUid]);
        free(xReads->a);
        xReads->a = buffer.a;
        xReads->n = buffer.n;
        xReads->m = buffer.m;
        xReads->len = totalLen;
        xReads->circ = 0;
        if(xReads->start != (xReads->a[0]>>32))
        {
            xReads->start = xReads->a[0]>>32;
            v = (keepUid<<1)+1;
            av = asg_arc_a(ug->g, v);
            nv = asg_arc_n(ug->g, v);
            for (k = 0; k < nv; k++)
            {
                if(av[k].del) continue;
                asg_arc_del(ug->g, av[k].ul>>32, av[k].v, 1);
                asg_arc_del(ug->g, av[k].v^1, av[k].ul>>32^1, 1);
            }
        }
        
        if(xReads->end != ((xReads->a[xReads->n-1]>>32)^1))
        {
            xReads->end = ((xReads->a[xReads->n-1]>>32)^1);
            v = (keepUid<<1);
            av = asg_arc_a(ug->g, v);
            nv = asg_arc_n(ug->g, v);
            for (k = 0; k < nv; k++)
            {
                if(av[k].del) continue;
                asg_arc_del(ug->g, av[k].ul>>32, av[k].v, 1);
                asg_arc_del(ug->g, av[k].v^1, av[k].ul>>32^1, 1);
            }
        }
    }

    for (i = 0; i < b_0->b.n; i++)
    {
        v = b_0->b.a[i];
        visit[v>>1] = 1;
        if(purge_g->seq[v>1].c != ALTER_LABLE) continue;
        asg_seq_drop(purge_g, v>1);
    }
}

void print_het_ovlp(p_g_t *pg, ma_ug_t *ug, hap_overlaps_list* ha, double filter_rate)
{
    uint32_t v, i, k, n_vtx = pg->pg_h_lev->n_seq * 2, nv, qn, qs, qe, tn, ts, te, as, ae, occ, ovlp, hetLen, homLen;
    asg_arc_t *av = NULL;
    hap_overlaps *x = NULL;
    p_node_t *a = NULL;
    int index;
    for (v = 0; v < n_vtx; v++)
    {
        av = asg_arc_a(pg->pg_h_lev, v);
        nv = asg_arc_n(pg->pg_h_lev, v);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            index = get_specific_hap_overlap(&(ha->x[av[i].ul>>33]), av[i].ul>>33, av[i].v>>1);
            if(index == -1) fprintf(stderr, "ERROR\n");
            x = &(ha->x[av[i].ul>>33].a.a[index]);

            qn = x->xUid;
            qs = x->x_beg_pos;
            qe = x->x_end_pos - 1;
            get_p_nodes(pg, &a, &occ, qn);
            for (k = 0, hetLen = 0, homLen = 0; k < occ; k++)
            {
                as = a[k].baseBeg;
                ae = a[k].baseEnd;
                ovlp = ((MIN(qe, ae) >= MAX(qs, as))? MIN(qe, ae) - MAX(qs, as) + 1 : 0);
                if(homLen + hetLen > 0 && ovlp == 0) break;
                if(ovlp == 0) continue;
                if(a[k].h_status == N_HET)
                {
                    homLen += ovlp;
                } 
                else if(asm_opt.polyploidy <= 2 && (a[k].h_status&P_HET))
                {
                    homLen += ovlp;
                }
                else
                {
                    hetLen += ovlp;
                }  
            }

            if(hetLen <= ((hetLen + homLen) * filter_rate))
            {
                ///all_ovlp->x[uId].a.a[i].status = DELETE;
                fprintf(stderr, "********XY********\n");
                print_hap_paf(ug, x);
            }


            tn = x->yUid;
            ts = x->y_beg_pos;
            te = x->y_end_pos - 1;
            get_p_nodes(pg, &a, &occ, tn);
            for (k = 0, hetLen = 0, homLen = 0; k < occ; k++)
            {
                as = a[k].baseBeg;
                ae = a[k].baseEnd;
                ovlp = ((MIN(te, ae) >= MAX(ts, as))? MIN(te, ae) - MAX(ts, as) + 1 : 0);
                if(homLen + hetLen > 0 && ovlp == 0) break;
                if(ovlp == 0) continue;
                if(a[k].h_status == N_HET)
                {
                    homLen += ovlp;
                } 
                else if(asm_opt.polyploidy <= 2 && (a[k].h_status&P_HET))
                {
                    homLen += ovlp;
                }
                else
                {
                    hetLen += ovlp;
                }                 
            }

            if(hetLen <= ((hetLen + homLen) * filter_rate))
            {
                ///all_ovlp->x[uId].a.a[i].status = DELETE;
                fprintf(stderr, "********YX********\n");
                print_hap_paf(ug, x);
            }
        }
    }

    for (v = 0; v < ha->num; v++)
    {
        for (i = 0; i < ha->x[v].a.n; i++)
        {
            if(ha->x[v].a.a[i].status == DELETE)
            {
                x = &(ha->x[v].a.a[i]);

                qn = x->xUid;
                qs = x->x_beg_pos;
                qe = x->x_end_pos - 1;
                get_p_nodes(pg, &a, &occ, qn);
                for (k = 0, hetLen = 0, homLen = 0; k < occ; k++)
                {
                    as = a[k].baseBeg;
                    ae = a[k].baseEnd;
                    ovlp = ((MIN(qe, ae) >= MAX(qs, as))? MIN(qe, ae) - MAX(qs, as) + 1 : 0);
                    if(homLen + hetLen > 0 && ovlp == 0) break;
                    if(ovlp == 0) continue;
                    if(a[k].h_status == N_HET)
                    {
                        homLen += ovlp;
                    } 
                    else if(asm_opt.polyploidy <= 2 && (a[k].h_status&P_HET))
                    {
                        homLen += ovlp;
                    }
                    else
                    {
                        hetLen += ovlp;
                    }  
                }

                if(hetLen <= ((hetLen + homLen) * filter_rate))
                {
                    fprintf(stderr, "********C(X)********hetLen-%u, homLen-%u\n", hetLen, homLen);
                    print_hap_paf(ug, x);
                }


                tn = x->yUid;
                ts = x->y_beg_pos;
                te = x->y_end_pos - 1;
                get_p_nodes(pg, &a, &occ, tn);
                for (k = 0, hetLen = 0, homLen = 0; k < occ; k++)
                {
                    as = a[k].baseBeg;
                    ae = a[k].baseEnd;
                    ovlp = ((MIN(te, ae) >= MAX(ts, as))? MIN(te, ae) - MAX(ts, as) + 1 : 0);
                    if(homLen + hetLen > 0 && ovlp == 0) break;
                    if(ovlp == 0) continue;
                    if(a[k].h_status == N_HET)
                    {
                        homLen += ovlp;
                    } 
                    else if(asm_opt.polyploidy <= 2 && (a[k].h_status&P_HET))
                    {
                        homLen += ovlp;
                    }
                    else
                    {
                        hetLen += ovlp;
                    }                 
                }

                if(hetLen <= ((hetLen + homLen) * filter_rate))
                {
                    fprintf(stderr, "********C(Y)********hetLen-%u, homLen-%u\n", hetLen, homLen);
                    print_hap_paf(ug, x);
                }
            }
        }
    }
}

void link_unitigs(asg_t *purge_g, ma_ug_t *ug, hap_overlaps_list* all_ovlp,
R_to_U* ruIndex, ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, asg_t *read_g, 
uint64_t* position_index, kvec_asg_arc_t_offset* u_buffer, kvec_t_i32_warp* tailIndex, 
kvec_t_i32_warp* prevIndex, int max_hang, int min_ovlp, kvec_asg_arc_t_warp* edge, uint8_t* visit,
hap_cov_t *cov)
{
    uint32_t v, n_vtx = purge_g->n_seq * 2, beg, end;
    long long nodeLen, baseLen, max_stop_nodeLen, max_stop_baseLen;
    buf_t b_0;
    memset(&b_0, 0, sizeof(buf_t));
    memset(visit, 0, purge_g->n_seq);
    for (v = 0; v < n_vtx; ++v) 
    {
        if(purge_g->seq[v>>1].c == ALTER_LABLE || purge_g->seq[v>>1].del || visit[v>>1]) continue;
        if(get_real_length(purge_g, v, NULL) != 1) continue;
        if(get_real_length(purge_g, v^1, NULL) != 0) continue;

        beg = v;
        b_0.b.n = 0;
        if(get_unitig(purge_g, NULL, beg, &end, &nodeLen, &baseLen, &max_stop_nodeLen, 
                    &max_stop_baseLen, 1, &b_0) == LOOP)
        {
            continue;
        }

        ///if(cov->link) collect_reverse_unitigs_purge(&b_0, cov->link, ug, all_ovlp);
        purge_merge(purge_g, ug, all_ovlp, &b_0, ruIndex, reverse_sources, coverage_cut, 
        read_g, position_index, u_buffer, tailIndex, prevIndex,max_hang, min_ovlp, edge, visit, cov);
    }
    free(b_0.b.a);
}

void print_all_purge_ovlp(ma_ug_t *ug, hap_overlaps_list* all_ovlp, const char* cmd)
{
    fprintf(stderr, "\n%s--->ug->u.n: %u\n", cmd, (uint32_t)ug->u.n);
    uint32_t v, uId, i;
    for (v = 0; v < all_ovlp->num; v++)
    {
        uId = v;
        for (i = 0; i < all_ovlp->x[uId].a.n; i++)
        {
            print_hap_paf(ug, &(all_ovlp->x[uId].a.a[i]));
        }
    }

}

inline int get_available_cnt(asg_t *g, uint32_t v, uint8_t* del, asg_arc_t* v_s)
{
    //v has direction
    if(del && del[v>>1]) return 0;
    uint32_t i, kv = 0;
    asg_arc_t *av = asg_arc_a(g, v);
    uint32_t nv = asg_arc_n(g, v);

    for (i = 0, kv = 0; i < nv; i++)
    {
        if(!av[i].del)
        {
            if(del && del[av[i].v>>1]) continue;
            if(v_s) v_s[kv] = av[i];
            kv++;
        }
    }

    return kv;
}

long long get_specific_contig_length(asg_t *g, uint8_t *del)
{
    asg_cleanup(g);
    uint32_t v, n_vtx = g->n_seq * 2, q_occ;
    uint8_t *mark = NULL;
    ///is a queue
	//kdq_t(uint64_t) *q;
    ///each node has two directions
    //q = kdq_init(uint64_t);


    mark = (uint8_t*)calloc(n_vtx, 1);
    
    long long totalLen = 0;
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t w, x, l, start, end, len;
        asg_arc_t arc;
        if (g->seq[v>>1].del || mark[v]) continue;
        if (get_available_cnt(g, v, del, NULL) == 0 && get_available_cnt(g, (v^1), del, NULL) != 0) continue;
        if (del[v>>1]) continue;

        mark[v] = 1;
		//q->count = 0, start = v, end = v^1, len = 0;
        q_occ =0, start = v, end = v^1, len = 0;
		// forward
		w = v;


        while (1) 
        {
			/**
			 * w----->x
			 * w<-----x
			 * that means the only suffix of w is x, and the only prefix of x is w
			 **/
			if (get_available_cnt(g, w, del, NULL) != 1) break;
            get_available_cnt(g, w, del, &arc);
			x = arc.v; // w->x
			if (get_available_cnt(g, x^1, del, NULL) != 1) break;

			/**
			 * another direction of w would be marked as used (since w has been used)
			**/
			mark[x] = mark[w^1] = 1;
			///l is the edge length, instead of overlap length
            ///note: edge length is different with overlap length
			///l = asg_arc_len(arc_first(g, w));
            get_available_cnt(g, w, del, &arc);
            l = ((uint32_t)((arc).ul));
			//kdq_push(uint64_t, q, (uint64_t)w<<32 | l);
            q_occ++;
			end = x^1, len += l;
			w = x;
			if (x == v) break;
		}


        //if (start != (end^1) || kdq_size(q) == 0) { // linear unitig
        if (start != (end^1) || q_occ == 0) { // linear unitig
			///length of seq, instead of edge
			l = g->seq[end>>1].len;
			//kdq_push(uint64_t, q, (uint64_t)(end^1)<<32 | l);
            q_occ++;
			len += l;
		} else { // circular unitig
			start = end = UINT32_MAX;
			goto add_unitig; // then it is not necessary to do the backward
		}

        // backward
		x = v;
		while (1) { // similar to forward but not the same
			if (get_available_cnt(g, x^1, del, NULL) != 1) break;
            get_available_cnt(g, x^1, del, &arc);
            w = arc.v ^ 1;
			if (get_available_cnt(g, w, del, NULL) != 1) break;
			mark[x] = mark[w^1] = 1;
			///l = asg_arc_len(arc_first(g, w));
            get_available_cnt(g, w, del, &arc);
            l = ((uint32_t)((arc).ul));
			///w is the seq id + direction, l is the length of edge
			///push element to the front of a queue
			//kdq_unshift(uint64_t, q, (uint64_t)w<<32 | l);
            q_occ++;

			start = w, len += l;
			x = w;
		}


        add_unitig:
        if (start != UINT32_MAX) mark[start] = mark[end] = 1;
        totalLen += len;
    }
    //kdq_destroy(uint64_t, q);
    free(mark);
    return totalLen;
}


void get_contig_length(ma_ug_t *ug, asg_t *g, uint64_t* primaryLen, uint64_t* alterLen)
{
    uint8_t *del = (uint8_t *)malloc(sizeof(uint8_t)*g->n_seq);
    uint32_t v, k;
    ma_utg_t* u = NULL;
    memset(del, 1, g->n_seq);
    (*primaryLen) = (*alterLen) = 0;

    for (v = 0; v < ug->g->n_seq; ++v) 
    {
        if(ug->g->seq[v].del) continue;
        if(ug->g->seq[v].c == ALTER_LABLE) continue;
        u = &(ug->u.a[v]);
        if(u->m == 0) continue;
        for (k = 0; k < u->n; k++)
        {
            del[u->a[k]>>33] = 0;
        }
    }
    (*primaryLen) = get_specific_contig_length(g, del);


    for (v = 0; v < g->n_seq; ++v) 
    {
        del[v] = 1 - del[v];
    }

    (*alterLen) = get_specific_contig_length(g, del);

    free(del);
}


int if_ploid_sample(ma_ug_t *ug, asg_t *read_g, R_to_U* ruIndex, 
ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, ma_sub_t* coverage_cut, 
hap_alignment_struct_pip* hap_buf, hap_overlaps_list* all_ovlp, hap_overlaps_list* back_all_ovlp, 
uint32_t minLen, double purge_threshold)
{
    asg_t* nsg = ug->g;
    uint64_t v, k, total_bases = 0, alter_bases = 0, primary_bases = 0, purge_bases = 0;
    kt_for(asm_opt.thread_num, hap_alignment_advance_worker, hap_buf, nsg->n_seq);
    
    filter_hap_overlaps_by_length(all_ovlp, minLen);
    normalize_hap_overlaps_advance(all_ovlp, back_all_ovlp, ug, read_g, reverse_sources, ruIndex);

    get_contig_length(ug, read_g, &primary_bases, &alter_bases);
    total_bases = primary_bases + alter_bases;
    // fprintf(stderr, "primary_bases: %lu\n", primary_bases);
    // fprintf(stderr, "alter_bases: %lu\n", alter_bases);
    // fprintf(stderr, "total_bases: %lu\n", total_bases);


    for (v = 0; v < all_ovlp->num; v++)
    {
        for (k = 0; k < all_ovlp->x[v].a.n; k++)
        {
            purge_bases += all_ovlp->x[v].a.a[k].x_end_pos - all_ovlp->x[v].a.a[k].x_beg_pos;
        }
    }
    purge_bases = purge_bases/2;
    ///fprintf(stderr, "purge_bases: %lu\n", purge_bases);
    alter_bases = alter_bases + purge_bases;
    ///fprintf(stderr, "new alter_bases: %lu\n", alter_bases);


    for (v = 0; v < all_ovlp->num; v++)
    {
        all_ovlp->x[v].a.n = 0;
    }

    for (v = 0; v < back_all_ovlp->num; v++)
    {
        back_all_ovlp->x[v].a.n = 0;
    }

    if(alter_bases > total_bases * purge_threshold) return 1;
    return 0;
}

int cmp_chain_score(const void * a, const void * b)
{
    if((*(hap_overlaps*)a).score < (*(hap_overlaps*)b).score) return 1;
    if((*(hap_overlaps*)a).score > (*(hap_overlaps*)b).score) return -1;
    
    return 0;
}
long long get_ovlp_len(long long a_beg, long long a_end, long long b_beg, long long b_end)
{
    long long ovlp = (long long)(MIN(a_end, b_end)) - (long long)(MAX(a_beg, b_beg)) + 1;
    return ovlp <= 0? 0 : ovlp;
}
void sort_hap_chain(hap_overlaps_list* all_ovlp)
{
    hap_overlaps *x = NULL, *p = NULL;
    uint32_t v, i, k, uId;    
    long long ovlp, xLen, pLen;
    kvec_t(hap_overlaps) pri; kv_init(pri);
    kvec_t(hap_overlaps) alt; kv_init(alt);

    for (v = 0; v < all_ovlp->num; v++)
    {
        uId = v;
        qsort(all_ovlp->x[uId].a.a, all_ovlp->x[uId].a.n, sizeof(hap_overlaps), cmp_chain_score);
        pri.n = alt.n = 0;        
        for (i = 0; i < all_ovlp->x[uId].a.n; i++)
        {
            x = &(all_ovlp->x[uId].a.a[i]);
            xLen = x->x_end_pos - x->x_beg_pos;
            for (k = 0; k < pri.n; k++)
            {
                p = &(pri.a[k]);
                pLen = p->x_end_pos - p->x_beg_pos;
                ovlp = get_ovlp_len(x->x_beg_pos, x->x_end_pos-1, p->x_beg_pos, p->x_end_pos-1);
                if(ovlp == 0) continue;
                if(ovlp >= (MIN(xLen, pLen))*0.5) break;
            }

            if(k < pri.n)
            {
                x->xUid = k;
                kv_push(hap_overlaps, alt, *x);
            }
            else
            {
                kv_push(hap_overlaps, pri, *x);
            } 
        }

    }

    kv_destroy(pri); kv_destroy(alt);
}

void remove_contained_haplotig(hap_overlaps_list* all_ovlp, ma_ug_t *ug, asg_t* nsg, asg_t *purge_g, hap_cov_t *cov)
{
    uint32_t v, i, uId, xUid;
    hap_overlaps *p = NULL;
    for (v = 0; v < all_ovlp->num; v++)
    {
        uId = v; p = NULL;
        for (i = 0; i < all_ovlp->x[uId].a.n; i++)
        {
            if(p == NULL || p->score < all_ovlp->x[uId].a.a[i].score)
            {
                p = &(all_ovlp->x[uId].a.a[i]);
            }
        }

        for (i = 0; i < all_ovlp->x[uId].a.n; i++)
        {
            if(all_ovlp->x[uId].a.a[i].type == YCX)
            {
                if(!filter_secondary_chain(p->score, all_ovlp->x[uId].a.a[i].score, 0.95))
                {
                    continue;
                } 
                
                xUid = all_ovlp->x[uId].a.a[i].xUid;

                nsg->seq[xUid].c = ALTER_LABLE;
                purge_g->seq[xUid].c = ALTER_LABLE;
                purge_g->seq[xUid].del = 1;

                all_ovlp->x[uId].a.a[i].status = DELETE;
                ///if(cov->link) collect_reverse_unitig_pair(cov->link, ug, &(all_ovlp->x[uId].a.a[i]));
                collect_trans_purge_cov(cov, ug, &(all_ovlp->x[uId].a.a[i]), 0);
            }

            ///print_hap_paf(ug, &(all_ovlp.x[uId].a.a[i]));
        }
    }

    // for (v = 0; v < all_ovlp.num; v++)
    // {
    //     uId = v;
    //     for (i = 0; i < all_ovlp.x[uId].a.n; i++)
    //     {
    //         if(all_ovlp.x[uId].a.a[i].type == YCX)
    //         {
    //             nsg->seq[all_ovlp.x[uId].a.a[i].xUid].c = ALTER_LABLE;
    //             purge_g->seq[all_ovlp.x[uId].a.a[i].xUid].c = ALTER_LABLE;
    //             purge_g->seq[all_ovlp.x[uId].a.a[i].xUid].del = 1;
    //             all_ovlp.x[uId].a.a[i].status = DELETE;
    //             if(link) collect_reverse_unitig_pair(link, ug, &(all_ovlp.x[uId].a.a[i]));
    //             collect_trans_purge_cov(cov, ug, &(all_ovlp.x[uId].a.a[i]), 0);
    //         }

    //         if(all_ovlp.x[uId].a.a[i].type == XCY)
    //         {
    //             nsg->seq[all_ovlp.x[uId].a.a[i].yUid].c = ALTER_LABLE;
    //             purge_g->seq[all_ovlp.x[uId].a.a[i].yUid].c = ALTER_LABLE;
    //             purge_g->seq[all_ovlp.x[uId].a.a[i].yUid].del = 1;
    //             all_ovlp.x[uId].a.a[i].status = DELETE;
    //             if(link) collect_reverse_unitig_pair(link, ug, &(all_ovlp.x[uId].a.a[i]));
    //             collect_trans_purge_cov(cov, ug, &(all_ovlp.x[uId].a.a[i]), 1);
    //         }
    //         ///print_hap_paf(ug, &(all_ovlp.x[uId].a.a[i]));
    //     }
    // }
}

void debug_p_g_t(p_g_t* pg, hap_cov_t *cov, asg_t *read_g)
{
    fprintf(stderr, "----------[M::%s]----------\n", __func__);
    uint32_t i, offset, v, sid, eid, spos, epos, p_status, p_uid, occ;
    p_node_t *t = NULL;
    ma_utg_t *u = NULL;
    p_node_t *a = NULL;

    for (v = 0; v < pg->ug->u.n; v++)
    {
        ///fprintf(stderr, "\nu->n: %u, uid: %u\n", (uint32_t)(pg->ug->u.a[v].n), v);
        get_p_nodes(pg, &a, &occ, v);
        for (i = 0; i < occ; i++)
        {
            if(a[i].c_ug_id != v) fprintf(stderr, "sbsbsbsbsb\n");
            ///fprintf(stderr, "sid: %u, eid: %u\n", a[i].nodeBeg, a[i].nodeEnd);
        }
    }
    

    for (v = 0, p_status = (uint32_t)-1, p_uid = (uint32_t)-1; v < pg->pg_het_node.n; v++)
    {
        t = &(pg->pg_het_node.a[v]);
        sid = t->nodeBeg;
        eid = t->nodeEnd;
        spos = t->baseBeg;
        epos = t->baseEnd;
        // fprintf(stderr, "sid: %u, eid: %u, spos: %u, epos: %u, t->b_ug_id: %u\n", 
        //                                                 sid, eid, spos, epos, (uint32_t)t->b_ug_id);
        // fprintf(stderr, "pg->pg_het_node.n: %u\n", (uint32_t)pg->pg_het_node.n);
        if(p_uid == t->c_ug_id && p_status == t->h_status)
        {
            fprintf(stderr, "ERROR-(-1)\n");
        } 
        p_status = t->h_status;
        p_uid = t->c_ug_id;

        u = &(pg->ug->u.a[t->c_ug_id]);
        ///fprintf(stderr, "u->n: %u, sid: %u, eid: %u\n", (uint32_t)u->n, sid, eid);
        for (i = offset = 0; i < u->n; i++)
        {
            if(i == sid)
            {
                if(spos != offset)
                {
                    fprintf(stderr, "ERROR-1\n");
                }
            }

            if(i == eid)
            {
                if(epos != (offset+read_g->seq[u->a[i]>>33].len - 1))
                {
                    fprintf(stderr, "ERROR-2, real end: %u\n", 
                                    (uint32_t)(offset+read_g->seq[u->a[i]>>33].len - 1));
                }
            }
            offset += (uint32_t)u->a[i];
            if(i >= sid && i <= eid)
            {
                if(cov->t_ch->ir_het[u->a[i]>>33] != t->h_status)
                {
                    fprintf(stderr, "ERROR-(-3): is_r_het: %u, h_status: %u\n", cov->t_ch->ir_het[u->a[i]>>33], t->h_status);
                }
            }
        }
    }
}


void print_p_g_t_interval(p_g_t* pg, hap_cov_t *cov)
{
    fprintf(stderr, "----------[M::%s]----------\n", __func__);
    uint32_t i, v, sid, eid;
    p_node_t *t = NULL;
    ma_utg_t *u = NULL;

    for (v = 0; v < pg->pg_het_node.n; v++)
    {
        t = &(pg->pg_het_node.a[v]);
        sid = t->nodeBeg;
        eid = t->nodeEnd;

        u = &(pg->ug->u.a[t->c_ug_id]);
        fprintf(stderr, "\nu->n=%u, sid=%u, eid=%u, h_status=%u\n", 
                                                (uint32_t)u->n, sid, eid, t->h_status);
        for (i = sid; i <= eid; i++)
        {
            fprintf(stderr, "id:i:%u------>utg%.6ul\n", 
                        (uint32_t)(u->a[i]>>33), (get_origin_uid(u->a[i]>>32, cov->t_ch, NULL, NULL)>>1)+1);
        }
    }
    fprintf(stderr, "----------[M::%s]----------\n", __func__);
}


p_g_t *init_p_g_t(ma_ug_t *ug, hap_cov_t *cov, asg_t *read_g)
{
    uint32_t v, uId, k, l, offset, l_pos, g_beg_idx, occ/**, ovlp, tLen, zLen**/;
    p_g_t *pg = NULL; CALLOC(pg, 1);
    pg->ug = ug;
    asg_t* nsg = pg->ug->g;
    ma_utg_t *u = NULL;
    p_node_t *t = NULL/**, *z = NULL**/;
    ///asg_arc_t *e = NULL;
    p_g_in_t *x = NULL;
    ///pg->pg_het = asg_init();
    pg->pg_h_lev = asg_init();
    kv_init(pg->pg_het_node);
    kv_init(pg->pg_h_lev_idx);

    for (v = 0; v < nsg->n_seq; v++)
    {
        uId = v;
        if(nsg->seq[uId].del || nsg->seq[uId].c == ALTER_LABLE)
        {
            asg_seq_set(pg->pg_h_lev, uId, 0, 1);
            pg->pg_h_lev->seq[uId].c = ALTER_LABLE;
            continue;
        } 
        
        asg_seq_set(pg->pg_h_lev, uId, ug->u.a[uId].len, 0);
        pg->pg_h_lev->seq[uId].c = PRIMARY_LABLE;
    }

    // if(asm_opt.polyploidy <= 2)
    // {
    //     for (v = 0; v < cov->t_ch->r_num; v++)
    //     {
    //         if(cov->t_ch->is_r_het[v]&P_HET) cov->t_ch->is_r_het[v] |= S_HET;
    //     }
    // } 

    for (v = 0; v < nsg->n_seq; v++)
    {
        uId = v;
        if(nsg->seq[uId].del || nsg->seq[uId].c == ALTER_LABLE) continue;

        u = &(ug->u.a[uId]);
        g_beg_idx = pg->pg_het_node.n;
        ///fprintf(stderr, "\n+v: %u, pg->pg_het_node.n: %u\n", v, (uint32_t)pg->pg_het_node.n);
        for (k = 1, l = 0, offset = 0, l_pos = 0; k <= u->n; ++k) 
        {   
            ///if (k == u->n || (!!cov->t_ch->is_r_het[u->a[k]>>33]) != (!!cov->t_ch->is_r_het[u->a[l]>>33]))
            if (k == u->n || cov->t_ch->ir_het[u->a[k]>>33] != cov->t_ch->ir_het[u->a[l]>>33])
            {
                kv_pushp(p_node_t, pg->pg_het_node, &t);
                t->c_ug_id = uId;
                t->h_status = cov->t_ch->ir_het[u->a[l]>>33];
                t->baseBeg = l_pos;
                t->baseEnd = offset + read_g->seq[u->a[k-1]>>33].len - 1;
                t->nodeBeg = l;
                t->nodeEnd = k - 1;
                ///if(t->b_ug_id == (uint32_t)-1) fprintf(stderr, "xxxx\n");
                ///asg_seq_set(pg->pg_het, pg->pg_het_node.n-1, t->baseEnd+1-t->baseBeg, 0);
                ///fprintf(stderr, "l: %u, k: %u, u->n: %u, t->h_status: %u\n", l, k, u->n, t->h_status);
                l = k;
                l_pos = offset + (uint32_t)u->a[k-1];
            }
            offset += (uint32_t)u->a[k-1];
        }

        occ = pg->pg_het_node.n - g_beg_idx;
        kv_pushp(p_g_in_t, pg->pg_h_lev_idx, &x);
        x->beg = g_beg_idx; x->occ = occ;
        ///fprintf(stderr, "-v: %u, pg->pg_het_node.n: %u\n", v, (uint32_t)pg->pg_het_node.n);
        
        // if(occ > 1)
        // {
        //     for (k = g_beg_idx; (k + 1) < pg->pg_het_node.n; ++k) 
        //     {
        //         t = &(pg->pg_het_node.a[k]); tLen = t->baseEnd + 1 - t->baseBeg;
        //         z = &(pg->pg_het_node.a[k+1]); zLen = z->baseEnd + 1 - z->baseBeg;

        //         ovlp = ((MIN(t->baseEnd, z->baseEnd) >= MAX(t->baseBeg, z->baseBeg))? 
        //                         MIN(t->baseEnd, z->baseEnd) - MAX(t->baseBeg, z->baseBeg) + 1 : 0);

        //         e = asg_arc_pushp(pg->pg_het);
        //         e->ol = ovlp;
        //         e->ul = (k<<1); e->ul <<= 32; e->ul += (tLen - ovlp);
        //         e->v = ((k+1)<<1); e->del = 0; e->el = e->no_l_indel = e->strong = 1;

        //         e = asg_arc_pushp(pg->pg_het);
        //         e->ol = ovlp;
        //         e->ul = ((k+1)<<1)+1; e->ul <<= 32; e->ul += (zLen - ovlp);
        //         e->v = (k<<1)+1; e->del = 0; e->el = e->no_l_indel = e->strong = 1;
        //     }
        // }
    }
    ///asg_cleanup(pg->pg_het);
    ///debug_p_g_t(pg, cov, read_g);
    ///print_p_g_t_interval(pg, cov);

    return pg;
}

void destory_p_g_t(p_g_t **pg)
{
    if(pg && (*pg))
    {
        kv_destroy((*pg)->pg_het_node);
        kv_destroy((*pg)->pg_h_lev_idx);
        asg_destroy((*pg)->pg_het);
        asg_destroy((*pg)->pg_h_lev);
        free((*pg)); 
        (*pg) = NULL;
    }
}

void chain_origin_trans_uid_by_purge(hap_overlaps *x, ma_ug_t *ug, hap_cov_t *cov, uint64_t* position_index)
{
    uint32_t pri_uid, aux_uid, r_x, r_y;
    hap_candidates hap_for, hap_rev, *hap = NULL;
    long long x_pos_beg, x_pos_end, y_pos_beg, y_pos_end;

    Get_rev(hap_for) = x->rev;
    Get_x_beg(hap_for) = x->x_beg_id; Get_x_end(hap_for) = x->x_end_id - 1;
    Get_y_beg(hap_for) = x->y_beg_id; Get_y_end(hap_for) = x->y_end_id - 1;
    r_x = determine_hap_overlap_type_advance(&hap_for, &(ug->u.a[x->xUid]), &(ug->u.a[x->yUid]),
    cov->ruIndex, cov->reverse_sources, cov->coverage_cut, cov->read_g, position_index, 
    cov->max_hang, cov->min_ovlp, x->xUid, x->yUid, &(cov->u_buffer), &(cov->tailIndex), 
    &(cov->prevIndex), &x_pos_beg, &x_pos_end, &y_pos_beg, &y_pos_end);

    // if(r_x != (uint32_t)-1)
    // {
    //     adjust_hap_overlaps_score(&(ug->u.a[x->xUid]), NULL, &(hap_for.score), 
    //     x->xUid, x->yUid, x_pos_beg, x_pos_end);
    // }   
    

    Get_rev(hap_rev) = x->rev;
    Get_x_beg(hap_rev) = x->y_beg_id; Get_x_end(hap_rev) = x->y_end_id - 1;
    Get_y_beg(hap_rev) = x->x_beg_id; Get_y_end(hap_rev) = x->x_end_id - 1;
    r_y = determine_hap_overlap_type_advance(&hap_rev, &(ug->u.a[x->yUid]), &(ug->u.a[x->xUid]),
    cov->ruIndex, cov->reverse_sources, cov->coverage_cut, cov->read_g, position_index, 
    cov->max_hang, cov->min_ovlp, x->yUid, x->xUid, &(cov->u_buffer), &(cov->tailIndex), 
    &(cov->prevIndex), &y_pos_beg, &y_pos_end, &x_pos_beg, &x_pos_end);

    // if(r_y != (uint32_t)-1)
    // {
    //     adjust_hap_overlaps_score(&(ug->u.a[x->yUid]), NULL, &(hap_rev.score), 
    //     x->yUid, x->xUid, y_pos_beg, y_pos_end);
    // }

    if(r_x == (uint32_t)-1 && r_y == (uint32_t)-1)
    {
        fprintf(stderr, "ERROR-purge\n");
        return;
    }

    Get_rev(hap_for) = x->rev;
    Get_x_beg(hap_for) = x->x_beg_id; Get_x_end(hap_for) = x->x_end_id - 1;
    Get_y_beg(hap_for) = x->y_beg_id; Get_y_end(hap_for) = x->y_end_id - 1;

    Get_rev(hap_rev) = x->rev;
    Get_x_beg(hap_rev) = x->y_beg_id; Get_x_end(hap_rev) = x->y_end_id - 1;
    Get_y_beg(hap_rev) = x->x_beg_id; Get_y_end(hap_rev) = x->x_end_id - 1;


    if(r_x != (uint32_t)-1 && r_y == (uint32_t)-1)
    {
        pri_uid = x->xUid; aux_uid = x->yUid; hap = &hap_for;
    }
    else if(r_x == (uint32_t)-1 && r_y != (uint32_t)-1)
    {
        aux_uid = x->xUid; pri_uid = x->yUid; hap = &hap_rev;
    }
    else
    {
        if(hap_for.score >= hap_rev.score)
        {
            pri_uid = x->xUid; aux_uid = x->yUid; hap = &hap_for;
        }
        else
        {
            aux_uid = x->xUid; pri_uid = x->yUid; hap = &hap_rev;
        }
    }

    determine_hap_overlap_type_advance(hap, &(ug->u.a[pri_uid]), &(ug->u.a[aux_uid]),
    cov->ruIndex, cov->reverse_sources, cov->coverage_cut, cov->read_g, position_index, 
    cov->max_hang, cov->min_ovlp, pri_uid, aux_uid, &(cov->u_buffer), &(cov->tailIndex), 
    &(cov->prevIndex), &x_pos_beg, &x_pos_end, &y_pos_beg, &y_pos_end);

    // adjust_hap_overlaps_score(&(ug->u.a[pri_uid]), NULL, &(hap->score), 
    // pri_uid, aux_uid, x_pos_beg, x_pos_end);

    uint64_t pri_len = ug->u.a[pri_uid].len, aux_len = ug->u.a[aux_uid].len;
    pri_uid <<= 1; aux_uid <<= 1; aux_uid += hap->rev;

    // uint32_t i_n = cov->t_ch->k_trans.n, i;

    chain_origin_trans_uid_by_distance(cov, cov->read_g, &pri_uid, 1, x_pos_beg, &pri_len, 
                                                    &aux_uid, 1, y_pos_beg, &aux_len, ug, RC_2, hap->score, __func__);
    
    // fprintf(stderr, "\nocc: %u\n", (uint32_t)(cov->t_ch->k_trans.n - i_n));
    // fprintf(stderr, "#s-utg%.6ul\t%u\t%u\td-utg%.6ul\t%u\t%u\trev(%u)\n", 
    // x->xUid+1, x->x_beg_pos, x->x_end_pos, x->yUid+1, x->y_beg_pos, x->y_end_pos, x->rev);
    // for (i = i_n; i < cov->t_ch->k_trans.n; i++)
    // {
    //     fprintf(stderr, "s-utg%.6ul\t%u\t%u\td-utg%.6ul\t%u\t%u\trev(%u)\n", 
    //     cov->t_ch->k_trans.a[i].qn+1, cov->t_ch->k_trans.a[i].qs, cov->t_ch->k_trans.a[i].qe, 
    //     cov->t_ch->k_trans.a[i].tn+1, cov->t_ch->k_trans.a[i].ts, cov->t_ch->k_trans.a[i].te, 
    //     cov->t_ch->k_trans.a[i].rev);
    // }
}

void collect_purge_trans_cov(ma_ug_t *ug, hap_overlaps_list* ha, hap_cov_t *cov, uint64_t* position_index)
{   
    uint32_t v, i;
    hap_overlaps *x = NULL;
    for (v = 0; v < ha->num; v++)
    {
        for (i = 0; i < ha->x[v].a.n; i++)
        {
            x = &(ha->x[v].a.a[i]);
            if(x->yUid < x->xUid) continue;
            chain_origin_trans_uid_by_purge(x, ug, cov, position_index);
        }
    }
}

/**
typedef struct {
    uint32_t qn, qs, qe;
    uint32_t tn, ts, te;
    uint32_t oid;
    uint8_t rev;
}scg_hits;

typedef struct {
    scg_hits *a;
    size_t n,m;
    kvec_t(uint64_t) idx;
}scg_hits_v;

#define scg_key_qtn(a) ((((uint64_t)(a).qn)<<32)|((uint64_t)(a).tn))
KRADIX_SORT_INIT(scg_qtn, scg_hits, scg_key_qtn, 8)
#define scg_key_qts(a) ((((uint64_t)(a).qs)<<32)|((uint64_t)(a).ts))
KRADIX_SORT_INIT(scg_qts, scg_hits, scg_key_qts, 8)
#define scg_key_qte(a) ((((uint64_t)(a).qe)<<32)|((uint64_t)(a).te))
KRADIX_SORT_INIT(scg_qte, scg_hits, scg_key_qte, 8)
#define scg_key_rev(a) ((a).rev)
KRADIX_SORT_INIT(scg_rev, scg_hits, scg_key_rev, member_size(scg_hits, rev))

inline void rev_scg_hits(scg_hits *p, spg_t *scg)
{
    if(p->rev){
        uint32_t t;
        p->ts = scg->ug->u.a[p->tn].len - p->ts - 1;
        p->te = scg->ug->u.a[p->tn].len - (p->te - 1) - 1;
        t = p->ts; p->ts = p->te; p->te = t; p->te++;
    }
}

#define arc_first(g, v) ((g)->arc[(g)->idx[(v)]>>32])
scg_hits_v *get_scg_hits_v(scg_hits_v *vp, spg_t *scg)
{   
    scg_hits_v *hh = NULL; CALLOC(hh, 1);
    ma_utg_v *u = &(scg->ug->u);
    uint32_t i, mn, *ma = NULL;
    uint64_t offset, *idx = NULL;

    for (i = 0; i < scg->idx.n; i++) {
        mn = (uint32_t)scg->idx.a[i];
        ma = scg->dst.a + (scg->idx.a[i]>>32);
    }
    



    return hh;
}

void refine_scg(spg_t *scg, ma_ug_t *lug, hap_overlaps_list *ha, hap_cov_t *cov, uint64_t* position_index)
{
    scg_hits_v vp; kv_init(vp);
    uint32_t v, i, k, st, c[2];
    hap_overlaps *x = NULL;
    u_trans_t *z = NULL;
    scg_hits *p = NULL;

    for (v = 0; v < cov->t_ch->k_trans.n; v++){
        z = &(cov->t_ch->k_trans.a[v]);
        if(z->del) continue;
        kv_pushp(scg_hits, vp, &p);
        p->rev = z->rev; p->oid = (uint32_t)-1;
        p->qn = z->qn; p->qs = z->qs; p->qe = z->qe;
        p->tn = z->tn; p->ts = z->ts; p->te = z->te;
        // rev_scg_hits(p, scg);
        kv_pushp(scg_hits, vp, &p);
        p->rev = z->rev; p->oid = (uint32_t)-1;
        p->qn = z->tn; p->qs = z->ts; p->qe = z->te;
        p->tn = z->qn; p->ts = z->qs; p->te = z->qe;
        // rev_scg_hits(p, scg);
    }
    
    for (v = 0; v < ha->num; v++){
        for (i = 0; i < ha->x[v].a.n; i++){
            x = &(ha->x[v].a.a[i]);
            st = cov->t_ch->k_trans.n;
            chain_origin_trans_uid_by_purge(x, lug, cov, position_index);
            for (k = st; k < cov->t_ch->k_trans.n; k++){
                z = &(cov->t_ch->k_trans.a[k]);
                if(z->del) continue;
                kv_pushp(scg_hits, vp, &p);
                p->rev = z->rev; p->oid = v;
                p->qn = z->qn; p->qs = z->qs; p->qe = z->qe;
                p->tn = z->tn; p->ts = z->ts; p->te = z->te;
                // rev_scg_hits(p, scg);
                kv_pushp(scg_hits, vp, &p);
                p->rev = z->rev; p->oid = v;
                p->qn = z->tn; p->qs = z->ts; p->qe = z->te;
                p->tn = z->qn; p->ts = z->qs; p->te = z->qe;
                // rev_scg_hits(p, scg);
            }
            cov->t_ch->k_trans.n = st;
        }
    }

    ///two scg_hits might be totally equal; must remove first
    radix_sort_scg_qtn(vp.a, vp.a + vp.n);
    for (st = 0, i = 1; i <= vp.n; ++i){
        if (i == vp.n || vp.a[i].qn != vp.a[st].qn || vp.a[i].tn != vp.a[st].tn){
           if(i - st > 1) radix_sort_scg_rev(vp.a+st, vp.a+i);
           for (v = st, c[0] = c[1] = 0; v < i; v++) c[vp.a[v].rev]++;
           if(c[0]>1) radix_sort_scg_qts(vp.a+st, vp.a+st+c[0]);
           if(c[1]>1) radix_sort_scg_qts(vp.a+st+c[0], vp.a+st+c[0]+c[1]);
           st = i;
        }
    }
    for (st = 0, i = 1; i <= vp.n; ++i){
        if (i == vp.n || vp.a[i].rev != vp.a[st].rev || 
                vp.a[i].qn != vp.a[st].qn || vp.a[i].tn != vp.a[st].tn || 
                vp.a[i].qs != vp.a[st].qs || vp.a[i].ts != vp.a[st].ts)
        {
           if(i - st > 1) radix_sort_scg_qte(vp.a+st, vp.a+i);
           st = i;
        }
    }
    for (st = 0, i = 1, k = 0; i <= vp.n; ++i){
        if (i == vp.n || vp.a[i].rev != vp.a[st].rev ||
                vp.a[i].qn != vp.a[st].qn || vp.a[i].tn != vp.a[st].tn || 
                vp.a[i].qs != vp.a[st].qs || vp.a[i].ts != vp.a[st].ts ||
                vp.a[i].qe != vp.a[st].qe || vp.a[i].te != vp.a[st].te)
        {
           vp.a[k] = vp.a[st];
           k++;
           st = i;
        }
    }
    ///build idx
    vp.n = k; kv_resize(uint64_t, vp.idx, scg->ug->u.n); vp.idx.n = scg->ug->u.n; 
    memset(vp.idx.a, 0, vp.idx.n*sizeof(uint64_t));
    for (st = 0, i = 1; i <= vp.n; ++i)
    {
        if (i == vp.n || vp.a[i].qn != vp.a[st].qn)
        {
            vp.idx.a[vp.a[st].qn] = (uint64_t)st << 32 | (i - st);
            st = i;
        }
    }


    kv_destroy(vp); kv_destroy(vp.idx);
}
**/
uint32_t seed_uid(ma_utg_t *vu, uint64_t* ps_idx, R_to_U* ruIndex, ma_ug_t *rug)
{
    int64_t v_i, v, w, w_i, wb, we, vb, ve, k;
    uint32_t uid, is_u;
    ma_utg_t *wu = NULL;
    for (v_i = 0; v_i < vu->n; v_i++) {
        v = vu->a[v_i]>>32;
        get_R_to_U(ruIndex, v>>1, &uid, &is_u);
        if(is_u == 0 || uid == (uint32_t)-1 || ps_idx[v>>1] == (uint64_t)-1) continue;
        w_i = (uint32_t)ps_idx[v>>1];
        wu = &(rug->u.a[uid]);
        w = wu->a[w_i]>>32;
        if((v>>1)!=(w>>1)) continue;
        vb = 0; ve = vu->n; ///[vb, ve)
        if(v == w){ ///[wb, we)
            wb = w_i - v_i;
            we = wb + vu->n;
            if(wb < 0 || we > wu->n) continue;
            for (k = 0; k < vu->n; k++){
                if((vu->a[k+vb]>>32) != (wu->a[k+wb]>>32)) break;
            }
            if(k >= vu->n) return uid;
        } else {
            wb = w_i + 1 - (ve - v_i);
            we = wb + vu->n;
            if(wb < 0 || we > wu->n) continue;
            for (k = 0; k < vu->n; k++){
                if((vu->a[k+vb]>>32) != ((wu->a[we-k-1]>>32)^1)) break;
            }
            if(k >= vu->n) return uid;
        }
    }
    return (uint32_t)-1;
}

void filter_ovlp_vecs(hap_overlaps_list* ha, uint32_t *a, uint32_t a_n)
{
    uint32_t st, k, i, m, v, w;
    int idx;
    radix_sort_ru32(a, a + a_n);
    for (st = 0, m = 0, k = 1; k <= a_n; k++){
        if(k == a_n || a[k] != a[st]){
            a[m++] = a[st];
            st = k;
        }
    }
    a_n = m;
    if(a_n < 2) return;
    for (k = 0; k < a_n; k++){
        v = a[k];
        for (i = k+1; i < a_n; i++) {
            w = a[i];
            idx = get_specific_hap_overlap(&(ha->x[v]), v, w);
            if(idx != -1) ha->x[v].a.a[idx].status = DELETE;
            idx = get_specific_hap_overlap(&(ha->x[w]), w, v);
            if(idx != -1) ha->x[w].a.a[idx].status = DELETE;
        }
    }
}

void filter_ovlp_scg(hap_overlaps_list* ha, uint64_t* ps_idx, R_to_U* ruIndex, ma_ug_t *rug, spg_t *scg)
{
    uint32_t i, k, v, *ma = NULL, mn, luid;
    ma_utg_v *pp = &(scg->ug->u);
    kvec_t(uint32_t) vv; kv_init(vv);
    for (i = 0; i < scg->idx.n; i++){
        ma = scg->dst.a + (scg->idx.a[i]>>32);
        mn = (uint64_t)scg->idx.a[i];
        if(mn < 2) continue;
        for (k = 0, vv.n = 0; k < mn; k++){
            luid = seed_uid(&(pp->a[ma[k]>>1]), ps_idx, ruIndex, rug);
            if(luid == (uint32_t)-1) {
                fprintf(stderr, "ERROR-scg\n");
                continue;
            }
            kv_push(uint32_t, vv, luid);
        }
        if(vv.n < 2) continue;
        filter_ovlp_vecs(ha, vv.a, vv.n);
    }

    for (v = 0; v < ha->num; v++){
        for (i = 0, k = 0; i < ha->x[v].a.n; i++){
            if(ha->x[v].a.a[i].status == DELETE) continue;
            ha->x[v].a.a[k++] = ha->x[v].a.a[i];
        }
        ha->x[v].a.n = k;
    }

    kv_destroy(vv);
}

void purge_dups(ma_ug_t *ug, asg_t *read_g, ma_sub_t* coverage_cut, ma_hit_t_alloc* sources, 
ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex, kvec_asg_arc_t_warp* edge, float density, 
uint32_t purege_minLen, int max_hang, int min_ovlp, float drop_ratio, uint32_t just_contain, 
uint32_t just_coverage, hap_cov_t *cov, uint32_t collect_p_trans, uint32_t collect_p_trans_f)
{
    p_g_t *pg = NULL;
    asg_t* nsg = ug->g;
    uint32_t v, rId, uId, i, offset;
    ma_utg_t* reads = NULL;
    uint64_t* position_index = NULL;
    if(cov) position_index = cov->pos_idx;
    else position_index = (uint64_t*)malloc(sizeof(uint64_t)*read_g->n_seq);
    memset(position_index, -1, sizeof(uint64_t)*read_g->n_seq);
    
    hap_overlaps_list all_ovlp;
    init_hap_overlaps_list(&all_ovlp, nsg->n_seq);
    hap_overlaps_list back_all_ovlp;
    init_hap_overlaps_list(&back_all_ovlp, nsg->n_seq);
    asg_arc_t t, *p = NULL;
    int r;
    hap_alignment_struct_pip hap_buf;
    long long k_mer_only, coverage_only;

    if(asm_opt.hom_global_coverage != -1)
    {
        hap_buf.cov_threshold = (asm_opt.hom_global_coverage_set?
                (((double)asm_opt.hom_global_coverage)*((double)HOM_PEAK_RATE)):(asm_opt.hom_global_coverage));
    }
    else
    {
        hap_buf.cov_threshold = get_read_coverage_thres(ug, read_g, ruIndex, position_index, 
        sources, coverage_cut, read_g->n_seq, COV_COUNT, &k_mer_only, &coverage_only);
    }
    
    
    for (v = 0; v < nsg->n_seq; v++)
    {
        uId = v;
        reads = &(ug->u.a[uId]);
        for (i = 0, offset = 0; i < reads->n; i++)
        {
            rId = reads->a[i]>>33;
            set_R_to_U(ruIndex, rId, uId, 1, &(read_g->seq[rId].c));

            position_index[rId] = offset;
            position_index[rId] = position_index[rId] << 32;
            position_index[rId] = position_index[rId] | (uint64_t)i;

            offset += (uint32_t)reads->a[i];
        }
    }

    // if(just_coverage == 0)
    // {
    //     ma_ug_seq(ug, read_g, coverage_cut, sources, edge, max_hang, min_ovlp, 0, 0);
    // }
    // init_ug_idx(ug, asm_opt.k_mer_length, asm_opt.polyploidy, 2, !just_coverage);
    
    init_hap_alignment_struct_pip(&hap_buf, asm_opt.thread_num, nsg->n_seq, ug, read_g,  
    sources, reverse_sources, ruIndex, coverage_cut, position_index, density, max_hang, min_ovlp, 
    0.1, &all_ovlp, cov);
    
    if(hap_buf.cov_threshold < 0)
    {
        if(if_ploid_sample(ug, read_g, ruIndex, sources, reverse_sources, coverage_cut, 
        &hap_buf, &all_ovlp, &back_all_ovlp, purege_minLen, 0.333))
        {
            ///if peak is het, coverage peak is more reliable 
            hap_buf.cov_threshold = coverage_only * HET_PEAK_RATE;
        }
        else
        {
            ///if peak is homo, k-mer peak is more reliable 
            hap_buf.cov_threshold = k_mer_only * HOM_PEAK_RATE;
        }
    }
    if(asm_opt.hom_global_coverage == -1) asm_opt.hom_global_coverage = hap_buf.cov_threshold;
    if(asm_opt.pur_global_coverage != -1) hap_buf.cov_threshold = asm_opt.pur_global_coverage;
    fprintf(stderr, "[M::%s] homozygous read coverage threshold: %d\n", __func__, asm_opt.hom_global_coverage_set?
                            asm_opt.hom_global_coverage:(int)(((double)asm_opt.hom_global_coverage)/((double)HOM_PEAK_RATE)));
    fprintf(stderr, "[M::%s] purge duplication coverage threshold: %lld\n", __func__, hap_buf.cov_threshold);
    if(just_coverage) goto end_coverage;

    kt_for(asm_opt.thread_num, hap_alignment_advance_worker, &hap_buf, nsg->n_seq);
    
    ///if(debug_enable) print_all_purge_ovlp(ug, &all_ovlp);
    filter_hap_overlaps_by_length(&all_ovlp, purege_minLen);

    // normalize_hap_overlaps_advance(&all_ovlp, &back_all_ovlp, ug, read_g, reverse_sources, ruIndex);
    pg = init_p_g_t(ug, cov, read_g);
    normalize_hap_overlaps_advance_by_p_g_t(&all_ovlp, &back_all_ovlp, ug, read_g, reverse_sources, ruIndex, pg, cov, 0.8);

    if(collect_p_trans && collect_p_trans_f == 0)
    {
        collect_purge_trans_cov(ug, &all_ovlp, cov, position_index);
    } 
    
    if(asm_opt.polyploidy <= 2)
    {
        mc_solve(&all_ovlp, cov->t_ch, NULL, ug, read_g, 0.8, R_INF.trio_flag, 1, NULL, 1, NULL, NULL, 1);
    } 
    
    if(collect_p_trans && collect_p_trans_f == 1)
    {
        collect_purge_trans_cov(ug, &all_ovlp, cov, position_index);
    } 

    ///normalize_hap_overlaps_advance(&all_ovlp, &back_all_ovlp, ug, read_g, reverse_sources, ruIndex);
    ///debug_hap_overlaps(&all_ovlp, &back_all_ovlp);

    remove_contained_haplotig(&all_ovlp, ug, nsg, pg->pg_h_lev, cov);
    
    if(just_contain == 0)
    {
        for (v = 0; v < all_ovlp.num; v++)
        {
            uId = v;
            if(pg->pg_h_lev->seq[uId].del || pg->pg_h_lev->seq[uId].c == ALTER_LABLE) continue;
            for (i = 0; i < all_ovlp.x[uId].a.n; i++)
            {
                if(all_ovlp.x[uId].a.a[i].status == DELETE) continue;
                if(pg->pg_h_lev->seq[all_ovlp.x[uId].a.a[i].xUid].c == ALTER_LABLE||
                   pg->pg_h_lev->seq[all_ovlp.x[uId].a.a[i].xUid].del||
                   pg->pg_h_lev->seq[all_ovlp.x[uId].a.a[i].yUid].c == ALTER_LABLE||
                   pg->pg_h_lev->seq[all_ovlp.x[uId].a.a[i].yUid].del)
                {
                    continue;
                }

                
                ///print_hap_paf(ug, &(all_ovlp.x[uId].a.a[i]));
                
                r = get_hap_arch(&(all_ovlp.x[uId].a.a[i]), ug->u.a[all_ovlp.x[uId].a.a[i].xUid].len, 
                ug->u.a[all_ovlp.x[uId].a.a[i].yUid].len, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
                
                if(r < 0) continue;
                p = asg_arc_pushp(pg->pg_h_lev);
                *p = t;
            }
        }

        asg_cleanup(pg->pg_h_lev);
		asg_symm(pg->pg_h_lev);
        ///may need to do transitive reduction
        clean_purge_graph(pg->pg_h_lev, drop_ratio, 1);

        // if(debug_enable) print_purge_gfa(ug, purge_g);
        // if(debug_enable) print_all_purge_ovlp(ug, &all_ovlp);
        /*******************************for debug************************************/
        // print_het_ovlp(pg, ug, &all_ovlp, 0.8);
        /*******************************for debug************************************/

        link_unitigs(pg->pg_h_lev, ug, &all_ovlp, ruIndex, reverse_sources, coverage_cut, read_g, position_index, 
        &(hap_buf.buf[0].u_buffer), &(hap_buf.buf[0].u_buffer_tailIndex), &(hap_buf.buf[0].u_buffer_prevIndex),
        max_hang, min_ovlp, edge, hap_buf.buf[0].visit, cov);
    }

    for (v = 0; v < all_ovlp.num; v++)
    {
        uId = v;
        if(pg->pg_h_lev->seq[uId].c == ALTER_LABLE)
        {
            ug->g->seq[uId].c = ALTER_LABLE;
        }
    }

    end_coverage:
    uint32_t is_Unitig;
    for (v = 0; v < ruIndex->len; v++)
    {
        get_R_to_U(ruIndex, v, &uId, &is_Unitig);
        if(is_Unitig == 1) ruIndex->index[v] = (uint32_t)-1;
    }
    asg_cleanup(nsg);
    destory_hap_overlaps_list(&all_ovlp);
    destory_hap_overlaps_list(&back_all_ovlp);
    if(cov) memset(position_index, -1, sizeof(uint64_t)*read_g->n_seq);
    else free(position_index);
    destory_hap_alignment_struct_pip(&hap_buf);
    destory_p_g_t(&pg);
    // if(just_coverage == 0)
    // {
    //     des_ug_idx();
    //     for (i = 0; i < ug->u.n; i++)
    //     {
    //         free(ug->u.a[i].s);
    //         ug->u.a[i].s = NULL;
    //     }
    // }
}
