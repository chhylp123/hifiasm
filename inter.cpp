#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <zlib.h>
#include "kseq.h" // FASTA/Q parser
#include "kthread.h"
#include "inter.h"
#include "Overlaps.h"
#include "CommandLines.h"
#include "htab.h"
#include "Hash_Table.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
	int w, k, bw, max_gap, is_HPC, hap_n;
} mg_idxopt_t;

typedef struct {
    ha_abufl_t *abl;
    st_mt_t sp;
    Candidates_list clist;
    overlap_region_alloc olist;
} ma_ov_buf_t;

typedef struct { // global data structure for kt_pipeline()
    const void *ha_flt_tab;
    const ha_pt_t *ha_idx;
    const mg_idxopt_t *opt;
	kseq_t *ks;
    int64_t chunk_size;
    uint64_t n_thread;
    uint64_t total_base;
    uint64_t total_pair;
} uldat_t;

typedef struct { // data structure for each step in kt_pipeline()
    const void *ha_flt_tab;
    const ha_pt_t *ha_idx;
	int n, m, sum_len;
	uint64_t *len, id;
	char **seq;
    ma_ov_buf_t **mo;
} utepdat_t;

ma_ov_buf_t **init_ma_ov_buf_t_arr(uint64_t n)
{
    uint64_t i;
    ma_ov_buf_t **p = NULL; CALLOC(p, 1); CALLOC(*p, n);
    for (i = 0; i < n; i++){
        kv_init((*p)[i].sp); (*p)[i].abl = ha_abufl_init();
        init_Candidates_list(&((*p)[i].clist));
	    init_overlap_region_alloc(&((*p)[i].olist));
    }
    return p;
}

void destory_ma_ov_buf_t_arr(ma_ov_buf_t ***p, uint64_t n)
{
    uint64_t i;
    for (i = 0; i < n; i++){
        kv_destroy((**p)[i].sp); ha_abufl_destroy((**p)[i].abl);
        destory_Candidates_list(&((**p)[i].clist));
	    destory_overlap_region_alloc(&((**p)[i].olist));
    }
    free((**p)); free((*p));
}


void init_mg_opt(mg_idxopt_t *opt, int is_HPC, int k, int w, int hap_n)
{
    opt->k = k; 
    opt->w = w;
    opt->hap_n = hap_n;
    opt->is_HPC = is_HPC;
    opt->bw = 2000;
    opt->max_gap = 5000;
}

void uidx_build(ma_ug_t *ug, mg_idxopt_t *opt)
{
    int flag = asm_opt.flag;
    asm_opt.flag |= HA_F_NO_HPC;
    ha_flt_tab = ha_ft_ug_gen(&asm_opt, &(ug->u), opt->is_HPC, opt->k, opt->w, 1, opt->hap_n*10);
    ha_idx = ha_pt_ug_gen(&asm_opt, ha_flt_tab, &(ug->u), opt->is_HPC, opt->k, opt->w, 1);
    asm_opt.flag = flag;
}

void uidx_destory()
{
    ha_ft_destroy(ha_flt_tab); 
    ha_pt_destroy(ha_idx);
}

static void worker_for_ul_alignment(void *data, long i, int tid) // callback for kt_for()
{
    utepdat_t *s = (utepdat_t*)data;
    ma_ov_buf_t *b = s->mo[tid];
    uint64_t len = s->len[i];
    char *r = s->seq[i];


}

static void *worker_ul_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
    uldat_t *p = (uldat_t*)data;
    ///uint64_t total_base = 0, total_pair = 0;
    if (step == 0) { // step 1: read a block of sequences
        int ret;
        uint64_t l;
		utepdat_t *s;
		CALLOC(s, 1);
        s->ha_flt_tab = p->ha_flt_tab; s->ha_idx = p->ha_idx; s->id = p->total_pair; 
        while ((ret = kseq_read(p->ks)) >= 0) 
        {
            if (p->ks->seq.l < (uint64_t)p->opt->k) continue;
            if (s->n == s->m) {
                s->m = s->m < 16? 16 : s->m + (s->n>>1);
                REALLOC(s->len, s->m);
                REALLOC(s->seq, s->m);
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
        utepdat_t *s = (utepdat_t*)in;
        s->mo = init_ma_ov_buf_t_arr(p->n_thread);
        /**
        CALLOC(s->pos, s->n);
        **/
        int i;
        kt_for(p->n_thread, worker_for_ul_alignment, s, s->n);
        for (i = 0; i < s->n; ++i) {
            free(s->seq[i]);
            p->total_base += s->len[i];
        }
        free(s->seq); free(s->len);
        destory_ma_ov_buf_t_arr(&(s->mo), p->n_thread);
		return s;
    }
    else if (step == 2) { // step 3: dump
        utepdat_t *s = (utepdat_t*)in;
        /**
        int i;
        for (i = 0; i < s->n; ++i) {
            // if(s->pos[i].a == NULL) continue;
            // kv_push(pe_hit_hap, p->hits, s->pos[i]);
            if(s->pos[i].s == (uint64_t)-1) continue;
            kv_push(pe_hit, p->hits.a, s->pos[i]);
        }
        free(s->pos);
        **/
        free(s);
    }
    return 0;
}

int alignment_ul_pipeline(uldat_t* sl, const enzyme *fn)
{
    double index_time = yak_realtime();
    int i;
    for (i = 0; i < fn->n; i++){
        gzFile fp;
        if ((fp = gzopen(fn->a[i], "r")) == 0) return 0;
        sl->ks = kseq_init(fp);
        kt_pipeline(3, worker_ul_pipeline, sl, 3);
        kseq_destroy(sl->ks);
        gzclose(fp);
    }
    fprintf(stderr, "[M::%s::%.3f] ==> Qualification\n", __func__, yak_realtime()-index_time);
    return 1;
}

int ul_align(mg_idxopt_t *opt, const enzyme *fn, void *ha_flt_tab, ha_pt_t *ha_idx)
{
    uldat_t sl; memset(&sl, 0, sizeof(sl));
    sl.ha_flt_tab = ha_flt_tab;
    sl.ha_idx = ha_idx;
    sl.opt = opt;   
    sl.chunk_size = 20000000;
    sl.n_thread = asm_opt.thread_num;
    alignment_ul_pipeline(&sl, fn);
    return 1;
}

void ul_resolve(ma_ug_t *ug, int hap_n)
{
    mg_idxopt_t opt;
    init_mg_opt(&opt, 0, 19, 10, hap_n);
    uidx_build(ug, &opt);
    ul_align(&opt, asm_opt.ar, ha_flt_tab, ha_idx);




    uidx_destory();
}