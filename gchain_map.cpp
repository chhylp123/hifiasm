#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <zlib.h>
#include <math.h>
#include "kseq.h" // FASTA/Q parser
#include "kavl.h"
#include "khash.h"
#include "kalloc.h"
#include "kthread.h"
#include "inter.h"
#include "Overlaps.h"
#include "CommandLines.h"
#include "htab.h"
#include "Hash_Table.h"
#include "Correct.h"
#include "Process_Read.h"
#include "Assembly.h"
#include "gchain_map.h"
KSEQ_INIT(gzFile, gzread)
void ul_map_lchain(ha_abufl_t *ab, uint32_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k, const ul_idx_t *uref, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, 
								 int max_n_chain, int apend_be, kvec_t_u8_warp* k_flag, overlap_region* f_cigar, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t is_accurate, uint32_t gen_off, double mcopy_rate, uint32_t chain_cutoff, uint32_t mcopy_khit_cut);

typedef struct { // global data structure for kt_pipeline()
    const void *ha_flt_tab;
    const ha_pt_t *ha_idx;
    const mg_idxopt_t *opt;
    const ma_ug_t *ug;
	const asg_t *rg;
	const ug_opt_t *uopt;
	const ul_idx_t *uu;
	ucr_file_t *ucr_s;
	kseq_t *ks;
    int64_t chunk_size;
    uint64_t n_thread;
    uint64_t total_base;
    uint64_t total_pair;
	uint64_t num_bases, num_corrected_bases, num_recorrected_bases;
	uint64_t remap, mini_cut;
} gmap_t;

typedef struct { // data structure for each step in kt_pipeline()
    const mg_idxopt_t *opt;
    const void *ha_flt_tab;
    const ha_pt_t *ha_idx;
    const ma_ug_t *ug;
	const asg_t *rg;
	const ug_opt_t *uopt;
	const ul_idx_t *uu;
	int n, m, sum_len;
	uint64_t *len, id;
	char **seq;
    // ha_mzl_v *mzs;///useless
    // st_mt_t *sps;///useless
	// mg_gchains_t **gcs;///useless
    mg_tbuf_t **buf;///useless
	ha_ovec_buf_t **hab;
	kv_ul_ov_t *res;
	// glchain_t *ll;
	// gdpchain_t *gdp;
	// glchain_t *sec_ll;
	uint64_t num_bases, num_corrected_bases, num_recorrected_bases, mini_cut;
	int64_t n_thread;
} sstep_t;

/**
static void worker_ul_map(void *data, long i, int tid) // callback for kt_for()
{
    sstep_t *s = (sstep_t*)data;
    ha_ovec_buf_t *b = s->hab[tid];
    kv_ul_ov_t *res = (s->res?(&(s->res[tid])):(NULL));
    mg_tbuf_t *buf = (s->buf?s->buf[tid]:NULL);
    int64_t winLen = MIN((((double)THRESHOLD_MAX_SIZE)/s->opt->diff_ec_ul), WINDOW);
    int fully_cov, abnormal; 
    assert(UL_INF.a[s->id+i].rlen == s->len[i]);
    // if(s->id+i!=43) return;

    ul_map_lchain(b->abl, i, s->seq[i], s->len[i], s->opt->w, s->opt->k, s->uu, &b->olist, &b->olist_hp, &b->clist, s->opt->bw_thres, 
        s->opt->max_n_chain, 1, NULL, &(b->tmp_region), NULL, &(b->sp), s->mini_cut, 0);
    
    clear_Cigar_record(&b->cigar1);
    clear_Round2_alignment(&b->round2);

    b->self_read.seq = s->seq[i]; b->self_read.length = s->len[i]; b->self_read.size = 0;
    lchain_align(&b->olist, s->uu, &b->self_read, &b->correct, &b->ovlp_read, &b->POA_Graph, &b->DAGCon,
            &b->cigar1, &b->hap, &b->round2, &b->r_buf, &(b->tmp_region.w_list), 0, 1, &fully_cov, &abnormal, s->opt->diff_ec_ul, winLen, NULL);

    // gl_chain_refine(&b->olist, &b->correct, &b->hap, bl, s->uu, s->opt->diff_ec_ul, winLen, s->len[i], km);
    gl_chain_refine_advance_combine(s->buf[tid], &(UL_INF.a[s->id+i]), &b->olist, &b->correct, &b->hap, &(s->sps[tid]), bl, &(s->gdp[tid]), s->uu, s->opt->diff_ec_ul, winLen, s->len[i], s->uopt, s->id+i, tid, NULL);
    
    memset(&b->self_read, 0, sizeof(b->self_read));
    if(UL_INF.a[s->id+i].dd) {
        free(s->seq[i]); s->seq[i] = NULL; b->num_correct_base++;
    }
    s->hab[tid]->num_read_base++;
}

static void *worker_gmap_work_ovec_pip(void *data, int step, void *in) // callback for kt_pipeline()
{
    gmap_t *p = (gmap_t*)data;
    if (step == 0) { // step 1: read a block of sequences
        int32_t ret; uint64_t l; sstep_t *s; CALLOC(s, 1);
        s->ha_flt_tab = p->ha_flt_tab; s->ha_idx = p->ha_idx; s->id = p->total_pair; 
        s->opt = p->opt; s->uu = p->uu; s->uopt = p->uopt; s->rg = p->rg; s->mini_cut = p->mini_cut;///need set
        while ((ret = kseq_read(p->ks)) >= 0) {
            if (p->ks->seq.l < (uint64_t)p->opt->k) continue;
            if (s->n == s->m) {
                s->m = s->m < 16? 16 : s->m + (s->n>>1);
                REALLOC(s->len, s->m);
                REALLOC(s->seq, s->m);
            }
            if(!(p->remap)) {
				append_ul_t(&UL_INF, NULL, p->ks->name.s, p->ks->name.l, NULL, 0, NULL, 0, P_CHAIN_COV, s->uopt, 0);        
			}
			l = p->ks->seq.l;
            MALLOC(s->seq[s->n], l);
            s->sum_len += l;
            memcpy(s->seq[s->n], p->ks->seq.s, l);
            s->len[s->n++] = l;
            if (s->sum_len >= p->chunk_size) break;            
        }
        p->total_pair += s->n;
        if (s->sum_len == 0) free(s);
        else return s;
    }
    else if (step == 1) { // step 2: alignment
        sstep_t *s = (sstep_t*)in; uint64_t i;
        CALLOC(s->hab, p->n_thread); 
		CALLOC(s->buf, p->n_thread); 
		if(!(p->remap)) CALLOC(s->res, p->n_thread);//for results
        for (i = 0; i < p->n_thread; ++i) {
			s->hab[i] = ha_ovec_init(0, 0, 1); s->buf[i] = mg_tbuf_init();
		}
        // kt_for(p->n_thread, worker_for_ul_scall_alignment, s, s->n);
       
        for (i = 0; i < p->n_thread; ++i) {
            s->num_bases += s->hab[i]->num_read_base;
            s->num_corrected_bases += s->hab[i]->num_correct_base;
            s->num_recorrected_bases += s->hab[i]->num_recorrect_base;
			ha_ovec_destroy(s->hab[i]); mg_tbuf_destroy(s->buf[i]); 
        }
        free(s->hab); free(s->buf); 
        return s;
    }
    else if (step == 2) { // step 3: dump
        sstep_t *s = (sstep_t*)in;
        uint64_t i, rid, sn = s->n;
        p->num_bases += s->num_bases;
        p->num_corrected_bases += s->num_corrected_bases;
        p->num_recorrected_bases += s->num_recorrected_bases;
		if(!(p->remap)) {
			for (i = 0; i < p->n_thread; ++i) {
				push_uc_block_t(s->uopt, &(s->res[i]), s->seq, s->len, s->id);
				kv_destroy(s->res[i]);
			}
			free(s->res);
        
			for (i = 0; i < sn; ++i) {
				rid = s->id + i; 
				if((UL_INF.n <= rid) || (UL_INF.n > rid && UL_INF.a[rid].rlen != s->len[i])) {///reads without alignment
					append_ul_t(&UL_INF, &rid, NULL, 0, s->seq[i], s->len[i], NULL, 0, P_CHAIN_COV, s->uopt, 0);
				} 
				free(s->seq[i]); 
			}
		} else {
			for (i = 0; i < sn; ++i) {
				rid = s->id + i;
				if(UL_INF.a[rid].dd == 0 && p->ucr_s && p->ucr_s->flag == 1) {
					assert(s->seq[i]);
					///for debug interval
					write_compress_base_disk(p->ucr_s->fp, rid, s->seq[i], s->len[i], &(p->ucr_s->u));
				}
				free(s->seq[i]);
			}
		}
		free(s->len); free(s->seq); free(s);
    }
    return 0;
}
**/

// int gmap_work_ovec(gmap_t* sl, const enzyme *fn)
// {
//     double index_time = yak_realtime();
//     int i;

// 	init_all_ul_t(&UL_INF, &R_INF);
//     for (i = 0; i < fn->n; i++){
//         gzFile fp;
//         if ((fp = gzopen(fn->a[i], "r")) == 0) return 0;
//         sl->ks = kseq_init(fp);
//         kt_pipeline(3, worker_gmap_work_ovec_pip, sl, 3);
//         kseq_destroy(sl->ks);
//         gzclose(fp);
//     }
// 	sl->hits.total_base = sl->total_base;
// 	sl->hits.total_pair = sl->total_pair;
//     fprintf(stderr, "[M::%s::%.3f] ==> Qualification\n", __func__, yak_realtime()-index_time);
// 	fprintf(stderr, "[M::%s::] ==> # reads: %lu, # bases: %lu\n", __func__, UL_INF.n, sl->total_base);
// 	fprintf(stderr, "[M::%s::] ==> # bases: %lu; # corrected bases: %lu; # recorrected bases: %lu\n", 
// 	__func__, sl->num_bases, sl->num_corrected_bases, sl->num_recorrected_bases);
// 	gen_ul_vec_rid_t(&UL_INF, &R_INF, NULL);
//     return 1;
// }