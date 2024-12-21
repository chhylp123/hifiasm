#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "Correct.h"
#include "Process_Read.h"
#include "ecovlp.h"
#include "kthread.h"
#include "htab.h"
#define HA_KMER_GOOD_RATIO 0.333
#define E_KHIT 31
#define CNS_DEL_E (0x7fffffffu)
#define del_cns_arc(z, arc_i) ((z).arc.a[(arc_i)].v == CNS_DEL_E)
#define CNS_DEL_V (0x1fffffffu)
#define del_cns_nn(z, nn_i) ((z).a[(nn_i)].sc == CNS_DEL_V)
#define REFRESH_N 128 
#define COV_W 3072

KDQ_INIT(uint32_t)

typedef struct {
	uint32_t v:31, f:1;
	uint32_t sc;
} cns_arc;
typedef struct {size_t n, m, nou; cns_arc *a; } cns_arc_v;

typedef struct {
	// uint16_t c:2, t:2, f:1, sc:3;
	uint32_t c:2, f:1, sc:29;
	cns_arc_v arc;
}cns_t;

typedef struct {
	size_t n, m; 
	cns_t *a;
	uint32_t si, ei, off, bn, bb0, bb1, cns_g_wl;
	kdq_t(uint32_t) *q;
}cns_gfa;

typedef struct {
	// chaining and overlapping related buffers
	UC_Read self_read, ovlp_read;
	Candidates_list clist;
	overlap_region_alloc olist;
	ha_abuf_t *ab;
	// int64_t num_read_base, num_correct_base, num_recorrect_base;
	uint64_t cnt[6], rr;
	haplotype_evdience_alloc hap;
	bit_extz_t exz;
	kv_ul_ov_t pidx;
	asg64_v v64;
	asg32_v v32;
	asg16_v v16;
	asg8_v v8q, v8t;

    kvec_t_u8_warp k_flag;
	st_mt_t sp;
	cns_gfa cns;	
} ec_ovec_buf_t0;

typedef struct {
	ec_ovec_buf_t0 *a;
	uint32_t n, rev;
    uint8_t *cr;
} ec_ovec_buf_t;

typedef struct {
	uint32_t n_thread, n_a, chunk_size, cn;
    FILE *fp;
} cal_ec_r_dbg_t;

typedef struct {
	ma_hit_t *a;
    size_t n, m;
    asg16_v ec;
} r_dbg_step_res_t;

typedef struct { // data structure for each step in kt_pipeline()
    ec_ovec_buf_t *buf;
    r_dbg_step_res_t *res;
    uint32_t si, ei;
} cal_ec_r_dbg_step_t;

ec_ovec_buf_t* gen_ec_ovec_buf_t(uint32_t n);
void destroy_ec_ovec_buf_t(ec_ovec_buf_t *p);


#define generic_key(x) (x)
KRADIX_SORT_INIT(ec16, uint16_t, generic_key, 2)
KRADIX_SORT_INIT(ec32, uint32_t, generic_key, 4)
KRADIX_SORT_INIT(ec64, uint64_t, generic_key, 8)

#define kdq_clear(q) ((q)->count = (q)->front = 0)

typedef struct {size_t n, m; asg16_v *a; uint8_t *f; } cc_v;
cc_v scc = {0, 0, NULL, NULL};
cc_v scb = {0, 0, NULL, NULL};
cc_v sca = {0, 0, NULL, NULL};

typedef struct {size_t n, m; char *a; UC_Read z; asg8_v q;} sl_v;


void h_ec_lchain(ha_abuf_t *ab, uint32_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k, All_reads *rref, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, 
								 int max_n_chain, int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t is_accurate, uint32_t gen_off, int64_t mcopy_num, double mcopy_rate, uint32_t chain_cutoff, uint32_t mcopy_khit_cut, uint64_t ocv_w);
void h_ec_lchain_amz(ha_abuf_t *ab, uint32_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k, All_reads *rref, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, 
								 int max_n_chain, int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t is_accurate, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t chain_cutoff, uint32_t mcopy_khit_cut, uint64_t ocv_w);
void h_ec_lchain_re_gen(ha_abuf_t *ab, uint32_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k, ha_pt_t *ha_idx, All_reads *rref, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, 
								 int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, 
								 int64_t max_skip, int64_t max_iter, int64_t max_dis, int64_t quick_check, double chn_pen_gap, double chn_pen_skip, UC_Read *tu, asg64_v *oidx, asg16_v *scc);
void h_ec_lchain_re_gen3(ha_abuf_t *ab, uint32_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k, ha_pt_t *ha_idx, All_reads *rref, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, 
								 int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, 
								 int64_t max_skip, int64_t max_iter, int64_t max_dis, int64_t quick_check, double chn_pen_gap, double chn_pen_skip, UC_Read *tu, asg64_v *oidx, asg16_v *scc);
uint64_t get_mz1(const char *str, int len, int w, int k, uint32_t rid, int is_hpc, ha_abuf_t *ab, const void *hf, ha_pt_t *ha_idx, int sample_dist, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, ha_pt_t *pt, int min_freq, int32_t dp_min_len, float dp_e, st_mt_t *mt, int32_t ws, int32_t is_unique, void *km, uint64_t beg_i);
void get_pi_ec_chain(ha_abuf_t *ab, uint64_t rid, uint64_t rl, uint32_t tid, char* ts, uint64_t tl, uint64_t mz_w, uint64_t mz_k, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, 
								int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, /**uint32_t is_accurate,**/ uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, 
								int64_t max_skip, int64_t max_iter, int64_t max_dis, int64_t quick_check, double chn_pen_gap, double chn_pen_skip);
void set_lchain_dp_op(uint32_t is_accurate, uint32_t mz_k, int64_t *max_skip, int64_t *max_iter, int64_t *max_dis, double *chn_pen_gap, double *chn_pen_skip, int64_t *quick_check);
void h_ec_lchain_re_gen_srt(ha_abuf_t *ab, ha_pt_t *ha_idx, overlap_region_alloc *olst, Candidates_list *cl);
uint64_t h_ec_lchain_re_gen_qry(ha_abuf_t *ab, uint64_t *k, uint64_t *l, uint64_t *i, uint64_t *idx_a, uint64_t idx_n, uint64_t *tid, uint64_t *trev);
uint64_t h_ec_lchain_re_chn(ha_abuf_t *ab, uint64_t si, uint64_t ei, uint32_t rid, char* rs, uint64_t rl, uint64_t tid, char* ts, uint64_t tl, uint64_t trev, uint64_t mz_w, uint64_t mz_k, overlap_region_alloc *olst, Candidates_list *cl, double bw_thres, 
								 int apend_be, uint64_t max_cnt, uint64_t min_cnt, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, int64_t max_skip, int64_t max_iter, int64_t max_dis, int64_t quick_check, double chn_pen_gap, double chn_pen_skip, tiny_queue_t *tq, asg16_v *scc, int64_t *n, int64_t *zn);
overlap_region* h_ec_lchain_fast(ha_abuf_t *ab, uint32_t rid, UC_Read *qu, UC_Read *tu, uint64_t mz_w, uint64_t mz_k, All_reads *rref, overlap_region_alloc *ol, Candidates_list *cl, bit_extz_t *exz, asg16_v *buf, asg64_v *srt_i, double bw_thres, 
								 int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t is_accurate, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, ma_hit_t_alloc *in0, ma_hit_t_alloc *in1, double sh);
void h_ec_lchain_fast_new(ha_abuf_t *ab, uint32_t rid, UC_Read *qu, UC_Read *tu, All_reads *rref, overlap_region_alloc *ol, Candidates_list *cl, bit_extz_t *exz, asg16_v *buf, asg64_v *srt_i, ma_hit_t_alloc *in0, ma_hit_t_alloc *in1, double sh);

ec_ovec_buf_t* gen_ec_ovec_buf_t(uint32_t n)
{
    uint32_t k; ec_ovec_buf_t0 *z = NULL;
    ec_ovec_buf_t *p = NULL; CALLOC(p, 1);
    p->n = n; CALLOC(p->a, p->n);
    for (k = 0; k < p->n; k++) {
        z = &(p->a[k]);
        init_UC_Read(&z->self_read);
	    init_UC_Read(&z->ovlp_read);
	    init_Candidates_list(&z->clist);
	    init_overlap_region_alloc(&z->olist);

        // init_fake_cigar(&(z->tmp.f_cigar));
        // memset(&(z->tmp.w_list), 0, sizeof(z->tmp.w_list));
        // CALLOC(z->tmp.w_list.a, 1); z->tmp.w_list.n = z->tmp.w_list.m = 1;

        // kv_init(z->b_buf.a);
        // kv_init(z->r_buf.a);
        kv_init(z->k_flag.a);
        kv_init(z->sp);
        kv_init(z->pidx);
	    kv_init(z->v64);
        kv_init(z->v32);
        kv_init(z->v16);
        kv_init(z->v8q);
        kv_init(z->v8t);
        init_bit_extz_t(&(z->exz), 31);

        z->ab = ha_abuf_init();
        
        InitHaplotypeEvdience(&z->hap);
        z->cns.q = kdq_init(uint32_t);
    }
    
    return p;
}

void destroy_cns_gfa(cns_gfa *p)
{
    size_t k;
    for (k = 0; k < p->m; k++) {
        kv_destroy(p->a[k].arc);
    }
    free(p->a); kdq_destroy(uint32_t, p->q);
}

void destroy_ec_ovec_buf_t(ec_ovec_buf_t *p)
{
    uint32_t k; ec_ovec_buf_t0 *z = NULL;
    for (k = 0; k < p->n; k++) {
        z = &(p->a[k]); z->rr = 0;
        destory_UC_Read(&z->self_read);
        destory_UC_Read(&z->ovlp_read);
        destory_Candidates_list(&z->clist);
	    destory_overlap_region_alloc(&z->olist);

        // destory_fake_cigar(&(z->tmp.f_cigar));
        // free(z->tmp.w_list.a); free(z->tmp.w_list.c.a);

        // kv_destroy(z->r_buf.a);
        kv_destroy(z->k_flag.a);
        kv_destroy(z->sp);
        kv_destroy(z->pidx);
	    kv_destroy(z->v64);
        kv_destroy(z->v32);
        kv_destroy(z->v16);
        kv_destroy(z->v8q);
        kv_destroy(z->v8t);
        destroy_bit_extz_t(&(z->exz));

        ha_abuf_destroy(z->ab);
        
        destoryHaplotypeEvdience(&z->hap);
        destroy_cns_gfa(&(z->cns));

    }
    free(p->a); free(p->cr); free(p);

    // fprintf(stderr, "[M::%s-chains] #->%lld\n", __func__, asm_opt.num_bases);
    // fprintf(stderr, "[M::%s-passed-chains-0] #->%lld\n", __func__, asm_opt.num_corrected_bases);
    // fprintf(stderr, "[M::%s-cis-chains-1] #->%lld\n", __func__, asm_opt.num_recorrected_bases);
}

inline void refresh_ec_ovec_buf_t0(ec_ovec_buf_t0 *z, uint64_t n)
{
    z->rr++;
    if((z->rr%n) == 0) {
        free(z->self_read.seq); memset(&(z->self_read), 0, sizeof(z->self_read));
        free(z->ovlp_read.seq); memset(&(z->ovlp_read), 0, sizeof(z->ovlp_read));

        destory_Candidates_list(&z->clist); memset(&(z->clist), 0, sizeof(z->clist));
        destory_overlap_region_alloc(&z->olist); memset(&(z->olist), 0, sizeof(z->olist)); init_overlap_region_alloc(&z->olist);

        kv_destroy(z->k_flag.a); kv_init(z->k_flag.a);
        kv_destroy(z->sp); kv_init(z->sp);
        kv_destroy(z->pidx); kv_init(z->pidx);
        kv_destroy(z->v64); kv_init(z->v64);
        kv_destroy(z->v32); kv_init(z->v32);
        kv_destroy(z->v16); kv_init(z->v16);
        kv_destroy(z->v8q); kv_init(z->v8q);
        kv_destroy(z->v8t); kv_init(z->v8t);

        destroy_bit_extz_t(&(z->exz)); init_bit_extz_t(&(z->exz), 31);

        ha_abuf_destroy(z->ab); z->ab = ha_abuf_init();

        destoryHaplotypeEvdience(&z->hap); memset(&(z->hap), 0, sizeof(z->hap)); InitHaplotypeEvdience(&z->hap);

        destroy_cns_gfa(&(z->cns)); memset(&(z->cns), 0, sizeof(z->cns)); z->cns.q = kdq_init(uint32_t);

        // z->rr = 1;
    }
}


void prt_chain(overlap_region_alloc *o)
{
    uint64_t k;
    for (k = 0; k < o->length; k++) {
        fprintf(stderr, "[M::%s]\t#%u\tlen::%lu\t%u\t%u\t%c\t#%u\tlen::%lu\t%u\t%u\tsc::%d\taln::%u\terr::%u\n", __func__, o->list[k].x_id, Get_READ_LENGTH(R_INF, o->list[k].x_id), o->list[k].x_pos_s, o->list[k].x_pos_e+1, "+-"[o->list[k].y_pos_strand], 
                                                                                        o->list[k].y_id, Get_READ_LENGTH(R_INF, o->list[k].y_id), o->list[k].y_pos_s, o->list[k].y_pos_e+1, o->list[k].shared_seed, o->list[k].align_length, o->list[k].non_homopolymer_errors);
    }    
}

overlap_region *fetch_aux_ovlp(overlap_region_alloc* ol) /// exactly same to gen_aux_ovlp
{
	if (ol->length + 1 >= ol->size) {
		uint64_t sl = ol->size;
        ol->size = ol->length + 1;
        kroundup64(ol->size);
        REALLOC(ol->list, ol->size);
        /// need to set new space to be 0
        memset(ol->list + sl, 0, sizeof(overlap_region)*(ol->size - sl));
	}
    ///debug for memory
    // if(ol->length + 1 >= ol->size) {
    //     fprintf(stderr, "[M::%s] length::%lu, size::%lu\n", __func__, ol->length, ol->size);
    // }
	return &(ol->list[ol->length+1]);
}

typedef struct {
	ul_ov_t *c_idx;
    asg64_v *idx;
    int64_t i, i0, srt_n, rr, ru;
    uint64_t mms, mme;
} cc_idx_t;


///[s, e)
int64_t extract_sub_cigar_mm(overlap_region *z, int64_t s, int64_t e, ul_ov_t *p, uint64_t *ct)
{
    int64_t wk = ovlp_cur_wid(*p), xk = ovlp_cur_xoff(*p), yk = ovlp_cur_yoff(*p), ck = ovlp_cur_coff(*p), os, oe, t;
    bit_extz_t ez; int64_t bd = ovlp_bd(*p), s0, e0; 
    s0 = ((int64_t)(z->w_list.a[wk].x_start)) + bd; 
    e0 = ((int64_t)(z->w_list.a[wk].x_end)) + 1 - bd;
    if(s < s0) s = s0; if(e > e0) e = e0;///exclude boundary
    if(s >= e) return -1;
    os = MAX(s, s0); oe = MIN(e, e0);
    if(oe <= os) return -1;

    set_bit_extz_t(ez, (*z), wk);
    if(!ez.cigar.n) return -1;
    int64_t cn = ez.cigar.n, op; int64_t ws, we, ovlp; 
    if((ck < 0) || (ck > cn)) {//(*ck) == cn is allowed
        ck = 0; xk = ez.ts; yk = ez.ps;
    }
    
    while (ck > 0 && xk >= s) {///x -> t; y -> p; first insertion and then match/mismatch
        --ck;
        op = ez.cigar.a[ck]>>14;
        if(op!=2) xk -= (ez.cigar.a[ck]&(0x3fff));
        if(op!=3) yk -= (ez.cigar.a[ck]&(0x3fff)); 
    }

     //some cigar will span s or e
    while (ck < cn && xk < e) {//[s, e)
        ws = xk; 
        op = ez.cigar.a[ck]>>14;
        ///op == 3: -> x; op == 2: -> y;
        if(op!=2) xk += (ez.cigar.a[ck]&(0x3fff));
        if(op!=3) yk += (ez.cigar.a[ck]&(0x3fff)); 
        ck++; we = xk;

        os = MAX(s, ws); oe = MIN(e, we);
        ovlp = ((oe>os)? (oe-os):0);
        if(op != 2) {
            if(!ovlp) continue;
        } else {///ws == we
            if(ws < s || ws >= e) continue;
        }
        

        if(op == 0) {
            for (t = os + 1; t < oe; t++) {
                ct[(t-s)<<1]++; ct[(t-s)<<1] += ((uint64_t)(0x100000000));
                ct[((t-s)<<1)+1]++; ct[((t-s)<<1)+1] += ((uint64_t)(0x100000000));
            }

            t = os;
            if(t < oe) {
                ct[(t-s)<<1]++; ct[(t-s)<<1] += ((uint64_t)(0x100000000));
                if(os > ws) {
                    ct[((t-s)<<1)+1]++; ct[((t-s)<<1)+1] += ((uint64_t)(0x100000000));
                }
            }
        } else if(op!=2) {
            for (t = os + 1; t < oe; t++) {
                ct[(t-s)<<1]++; 
                ct[((t-s)<<1)+1]++;
            }

            t = os;
            if(t < oe) {
                ct[(t-s)<<1]++; 
                if(os > ws) {
                    ct[((t-s)<<1)+1]++;
                }
            }
        } else {
            ct[((ws-s)<<1)+1]++; ///ct[((ws-s)<<1)+1] += ((uint64_t)(0x100000000));
        }
    }

    ovlp_cur_xoff(*p) = xk; ovlp_cur_yoff(*p) = yk; ovlp_cur_coff(*p) = ck; ovlp_cur_ylen(*p) = 0;
    
    return 1;
}

#define simp_vote_len 6

///[s, e)
uint32_t extract_sub_cigar_ii(overlap_region *z, int64_t ql, All_reads *rref, int64_t s, int64_t e, int64_t iws, int64_t iwe, UC_Read* tu, ul_ov_t *p)
{
    int64_t wk = ovlp_cur_wid(*p), xk = ovlp_cur_xoff(*p), yk = ovlp_cur_yoff(*p), ck = ovlp_cur_coff(*p), os, oe, ol;
    bit_extz_t ez; int64_t bd = ovlp_bd(*p), s0, e0, ii[2], it[2]; uint32_t res = (uint32_t)-1;
    s0 = ((int64_t)(z->w_list.a[wk].x_start)) + bd; 
    e0 = ((int64_t)(z->w_list.a[wk].x_end)) + 1 - bd;
    if(s < s0) s = s0; if(e > e0) e = e0;///exclude boundary
    if(s > e) return -1;///it is possible s == e
    os = MAX(s, s0); oe = MIN(e, e0);
    if(oe < os) return -1;///it is possible os == oe
    // fprintf(stderr, "[M::%s] s0::%ld, e0::%ld, iws::%ld, iwe::%ld\n", __func__, s0, e0, iws, iwe);
    ///make sure that this alignment block could cover the whole [iws, iwe) -> s0 < iws && e0 > iwe
    // if((s0 >= iws) || (e0 <= iwe)) return -1;///!(s0 < iws && e0 > iwe) -> only consider the alignment that could cover the whole [s, e)
    if(!(((s0 < iws) || (s0 == 0)) && ((e0 > iwe) || (e0 == ql)))) return -1;///!(s0 < iws && e0 > iwe) -> only consider the alignment that could cover the whole [s, e)

    set_bit_extz_t(ez, (*z), wk);
    if(!ez.cigar.n) return -1;
    int64_t cn = ez.cigar.n; uint16_t op; int64_t ws, we, wts, wte, ovlp, cc = 0, cci; 
    if((ck < 0) || (ck > cn)) {//(*ck) == cn is allowed
        ck = 0; xk = ez.ts; yk = ez.ps;
    }
    
    while (ck > 0 && xk >= s) {///x -> t; y -> p; first insertion and then match/mismatch
        --ck;
        op = ez.cigar.a[ck]>>14;
        if(op!=2) xk -= (ez.cigar.a[ck]&(0x3fff));
        if(op!=3) yk -= (ez.cigar.a[ck]&(0x3fff)); 
    }

    // char cm[4]; cm[0] = 'M'; cm[1] = 'S'; cm[2] = 'I'; cm[3] = 'D';  
    //some cigar will span s or e
    ii[0] = ii[1] = it[0] = it[1] = -1; res = cc = 0;
    while (ck < cn && xk < e) {//[s, e)
        ws = xk; wts = yk;
        op = ez.cigar.a[ck]>>14; ol = (ez.cigar.a[ck]&(0x3fff));
        ///op == 3: -> x; op == 2: -> y;
        if(op!=2) xk += ol;
        if(op!=3) yk += ol; 
        ck++; we = xk; wte = yk;

        // if(s == 10480) {
        //     fprintf(stderr, "[%ld, %ld)\t%c\n", ws, we, cm[op]);
        // }

        os = MAX(s, ws); oe = MIN(e, we);
        ovlp = ((oe>os)? (oe-os):0);

        if(s == e) {///insertion in comparsion with the reference
            if(op != 0 || ws >= s || we <= e || e != iwe || s != iws) continue;///must be a match
        } else {
            if(op != 2) {
                if(!ovlp) continue;
            } else {///ws == we
                if(ws < s || ws >= e) continue;
            }
        }
        
        
        if(ii[0] == -1) {
            ii[0] = os;
            if(op < 2) {
                it[0] = os - ws + wts;
            } else {///op == 2: more y; p == 3: more x
                it[0] = wts;
            } 
        }

        ii[1] = oe;
        if(op < 2) {
            it[1] = oe - ws + wts;
        } else {///op == 2: more y; p == 3: more x
            it[1] = wte;
        } 




        if(op != 2) ol = oe-os;
        cc += ol;
        // if(s == 11851 && e == 11853) {
        // if(!ol) fprintf(stderr, "%ld%c", ol, cm[op]);
        // }
        if(cc <= simp_vote_len) {
            for (cci = 0; cci < ol; cci++) {
                res <<= 2; res |= op;
            }
        }
    }

    while (ck < cn && xk <= e) {//[s, e)
        ws = xk; wts = yk;
        op = ez.cigar.a[ck]>>14; ol = (ez.cigar.a[ck]&(0x3fff));
        if(op != 2) break;
        yk += (ez.cigar.a[ck]&(0x3fff)); 
        ck++; we = xk; wte = yk;
        if(ws >= s && ws <= e) {
            
            if(ii[0] == -1) {
                ii[0] = ws; it[0] = wts;
            }
            ii[1] = we; it[1] = wte;
            
            cc += ol;
            // if(s == 11851 && e == 11853) {
            //     fprintf(stderr, "%ld%c", ol, cm[op]);
            // }
            if(cc <= simp_vote_len) {
                for (cci = 0; cci < ol; cci++) {
                    res <<= 2; res |= op;
                }
            }
        }
    }
    // if(s == 11851 && e == 11853) {
    //     fprintf(stderr, "\tx::[%ld, %ld)\ty::[%ld, %ld)\tcc::%ld\n", ii[0], ii[1], it[0], it[1], cc);
    // }
    if((cc <= simp_vote_len) 
            && (ii[1] >= ii[0]) && (ii[1] - ii[0] <= simp_vote_len) 
            && (it[1] >= it[0]) && (it[1] - it[0] <= simp_vote_len)) {
        // ii[0] = ii[0] - s; ii[1] = e - ii[1];
        if((ii[0] == iws) && (ii[1] == iwe)) {
            op = cc; op <<= 12; res |= op;

            char *ystr = NULL; res <<= 16; cc = it[1] - it[0]; op = 0;
            if(cc > 0) {
                UC_Read_resize(*tu, (it[1] - it[0])); ystr = tu->seq;
                recover_UC_Read_sub_region(ystr, it[0], (it[1] - it[0]), z->y_pos_strand, rref, z->y_id);
                
                for (cci = 0; cci < cc; cci++) {
                    op <<= 2; op |= seq_nt6_table[(uint32_t)(ystr[cci])];
                }
            } 
            res |= op;

            op = it[1] - it[0]; op <<= 12; res |= op;
        } else {
            res = (uint32_t)-1;
        }
    } else {
        res = (uint32_t)-1;
    }
    

    ovlp_cur_xoff(*p) = xk; ovlp_cur_yoff(*p) = yk; ovlp_cur_coff(*p) = ck; ovlp_cur_ylen(*p) = 0;
    
    return res;
}

typedef struct {
	All_reads *rref;
    UC_Read *tu;    
    uint64_t s, e, n0, n1, id, rev;
} rr_seq_t;

inline void insert_cns_arc(cns_gfa *cns, uint32_t src, uint32_t des, uint32_t is_ou, uint32_t plus0, uint32_t rid)
{
    if(src >= cns->n) {
        fprintf(stderr, "[M::%s] rid::%u, src::%u, des::%u, (*cns).n::%u\n", __func__, rid, src, des, (uint32_t)(*cns).n);
        exit(1);
    }
    cns_arc *p, t; kv_pushp(cns_arc, (*cns).a[src].arc, &p);
    p->f = 0; p->sc = plus0; p->v = des;
    if(is_ou) {
        (*cns).a[src].arc.nou++;
        if((*cns).a[src].arc.nou < (*cns).a[src].arc.n) {
            t = (*cns).a[src].arc.a[(*cns).a[src].arc.nou-1]; 
            (*cns).a[src].arc.a[(*cns).a[src].arc.nou-1] = *p;
            *p = t;
        }
    } 
}

inline uint32_t insert_cns_node(cns_gfa *cns)
{
    cns_t *p; uint32_t m0; 
    if (((*cns)).n == ((*cns)).m) { 
        m0 = ((*cns)).m;
        ((*cns)).m = ((*cns)).m? ((*cns)).m<<1 : 2; 
        ((*cns)).a = (cns_t*)realloc(((*cns)).a, sizeof(cns_t) * ((*cns)).m); 
        if(((*cns)).m > m0) {
            memset(((*cns)).a + m0, 0, sizeof(cns_t)*(((*cns)).m-m0));
        }
    } 
    *(&p) = &((*cns)).a[((*cns)).n++];
    p->arc.n = p->arc.nou = 0;
    p->c = p->f = p->sc = 0;
    return ((*cns)).n - 1;
}

inline uint32_t add_cns_arc(cns_gfa *cns, uint32_t src, uint32_t des, uint32_t is_ou, uint32_t plus)
{
    uint32_t k, s, e;
    if(is_ou) {
        s = 0; e = (*cns).a[src].arc.nou;
    } else {
        s = (*cns).a[src].arc.nou; e = (*cns).a[src].arc.n;
    }

    for (k = s; k < e; k++) {
        if((*cns).a[src].arc.a[k].v == des) {
            (*cns).a[src].arc.a[k].sc += plus;
            break;
        }
    }
    
    return ((k < e)?(1):(0));
}

inline void prt_cns_arc(cns_gfa *cns, uint32_t src, const char* cmd)
{
    uint32_t k;
    fprintf(stderr, "\n%s\t[M::%s] src::%u, sc::%u, c::%u\n", cmd, __func__, src, (*cns).a[src].sc, (*cns).a[src].c);
    for (k = 0; k < (*cns).a[src].arc.n; k++) {
        fprintf(stderr, "%s\t[M::%s] des::%u, sc::%u, is_ou::%u\n", cmd, __func__, (*cns).a[src].arc.a[k].v, (*cns).a[src].arc.a[k].sc, k<(*cns).a[src].arc.nou?1:0);
    }
}

inline uint32_t get_cns_arc_bp(cns_gfa *cns, uint32_t src, uint32_t bp, uint32_t is_ou, uint32_t av_bp)
{
    uint32_t k, s, e;
    if(is_ou) {
        s = 0; e = (*cns).a[src].arc.nou;
    } else {
        s = (*cns).a[src].arc.nou; e = (*cns).a[src].arc.n;
    }

    for (k = s; k < e; k++) {
        if((*cns).a[src].arc.a[k].v == 0 || (*cns).a[src].arc.a[k].v == 1) continue;
        if(av_bp && (*cns).a[src].arc.a[k].v >= (*cns).bb0 && (*cns).a[src].arc.a[k].v < (*cns).bb1) continue;///no backbone
        if((*cns).a[(*cns).a[src].arc.a[k].v].c == bp) {
            return k;
        }
    }
    
    return ((uint32_t)-1);
}

inline uint32_t add_cns_arc_bp(cns_gfa *cns, uint32_t src, uint32_t bp, uint32_t plus0, uint32_t rid, uint32_t av_bp)
{
    uint32_t rr, des;
    rr = get_cns_arc_bp(cns, src, bp, 1, av_bp);
    if(rr != ((uint32_t)-1)) {///find an existing node
        des = (*cns).a[src].arc.a[rr].v;
        (*cns).a[des].sc++;
        (*cns).a[src].arc.a[rr].sc += plus0;

        rr = add_cns_arc(cns, des, src, 0, plus0); 
        // if(rr == 0) {
        //     fprintf(stderr, "[M::%s] src::%u -> des::%u\n", __func__, src, des);
        //     prt_cns_arc(cns, src);
        //     prt_cns_arc(cns, des);
        // }
        assert(rr);

        return des;
    } else {///create a new node
        des = insert_cns_node(cns);
        (*cns).a[des].sc++; (*cns).a[des].c = bp;
        insert_cns_arc(cns, src, des, 1, plus0, rid); 
        insert_cns_arc(cns, des, src, 0, plus0, rid);
    }

    return des;
}

void init_cns_g(cns_gfa *cns, char *s, uint64_t sl, uint32_t rid)
{
    uint32_t m0 = cns->m, m1 = sl + 2, k; cns_t *p;
    if ((*cns).m < (m1)) { ///equal to kv_resize()
        (*cns).m = (m1); 
        (--((*cns).m), ((*cns).m)|=((*cns).m)>>1, ((*cns).m)|=((*cns).m)>>2, ((*cns).m)|=((*cns).m)>>4, ((*cns).m)|=((*cns).m)>>8, ((*cns).m)|=((*cns).m)>>16, ++((*cns).m)); 
        (*cns).a = (cns_t*)realloc((*cns).a, sizeof(cns_t) * (*cns).m); 
        if((*cns).m > m0) {
            memset((*cns).a + m0, 0, sizeof(cns_t)*((*cns).m-m0));
        }
    }
    (*cns).n = 0; (*cns).si = 0; (*cns).ei = 1; (*cns).off = 2;

    p = &((*cns).a[(*cns).n++]); p->arc.nou = p->arc.n = p->c = p->f = p->sc = 0; ///beg
    p = &((*cns).a[(*cns).n++]); p->arc.nou = p->arc.n = p->c = p->f = p->sc = 0; ///end
    (*cns).bb0 = (*cns).n;

    for (k = 0; k < sl; k++) {
        p = &((*cns).a[(*cns).n++]); p->arc.nou = p->arc.n = p->f = 0;
        p->c = seq_nt6_table[(uint32_t)(s[k])]; p->sc = 1;
        
        if(k + 1 < sl) insert_cns_arc(cns, k + (*cns).off, k + 1 + (*cns).off, 1, 1, rid);

        if(k > 0) insert_cns_arc(cns, k + (*cns).off, k - 1 + (*cns).off, 0, 1, rid);
    }
    
    if(sl) {
        insert_cns_arc(cns, (*cns).si, 0 + (*cns).off, 1, 1, rid); insert_cns_arc(cns, 0 + (*cns).off, (*cns).si, 0, 1, rid);
        insert_cns_arc(cns, sl - 1 + (*cns).off, (*cns).ei, 1, 1, rid); insert_cns_arc(cns, (*cns).ei, sl - 1 + (*cns).off, 0, 1, rid);
    } else {
        insert_cns_arc(cns, (*cns).si, (*cns).ei, 1, 1, rid); 
        insert_cns_arc(cns, (*cns).ei, (*cns).si, 0, 1, rid);
    }

    // prt_cns_arc(cns, 0, __func__);
    // prt_cns_arc(cns, 1, __func__);

    (*cns).bn = (*cns).n; (*cns).bb1 = (*cns).n;
}

///[s, e)
uint32_t push_cns_c0(cns_gfa *cns, uint64_t s0, uint64_t s, uint64_t e, uint32_t plus0, uint32_t rid)
{
    if(s > e) return s0;///it is possible that s == e
    uint32_t rr, k, re; 
    
    // rr = add_cns_arc(cns, s0, s, 1, plus0); assert(rr);
    // rr = add_cns_arc(cns, s, s0, 0, plus0); assert(rr);
    if(!add_cns_arc(cns, s0, s, 1, plus0)) {
        insert_cns_arc(cns, s0, s, 1, plus0, rid);
        insert_cns_arc(cns, s, s0, 0, plus0, rid);
    } else {
        rr = add_cns_arc(cns, s, s0, 0, plus0); assert(rr);
    }
    (*cns).a[s].sc++; re = s;

    for (k = s + 1; k < e; k++) {
        rr = add_cns_arc(cns, k-1, k, 1, 1); 
        // if(!rr) {
        //     fprintf(stderr, "[M::%s] s0::%u, s::%u\n", __func__, k-1, k);
        //     prt_cns_arc(cns, k-1, __func__); prt_cns_arc(cns, k, __func__);
        // }
        assert(rr);



        rr = add_cns_arc(cns, k, k-1, 0, 1); assert(rr);
        (*cns).a[k].sc++; re = k;
    }
    
    return re;
}

uint32_t trace_cns_bp(cns_gfa *cns, uint64_t s0, char *tstr, uint64_t tl, asg32_v* b32, uint32_t plus0, uint32_t *rn, uint64_t max_trace, uint32_t av_bp)
{
    (*rn) = s0;
    if(tl <= 0) return 0;
    // fprintf(stderr, "\n[M::%s] tl::%lu\n", __func__, tl);
    uint32_t k, i, s, e, m, bp, nm, bi, bn0, src, des, ff = 0; b32->n = 0;

    kv_push(uint32_t, (*b32), s0); kv_push(uint32_t, (*b32), ((uint32_t)-1)); nm = 2;
    // if(s0 == 2863) {
    //     fprintf(stderr, "***0***[M::%s] s::%lu\tb32->n::%u\n", __func__, s0, (uint32_t)b32->n);
    // }

    for (i = 0; (i < tl) && (!ff); i++) {
        bp = seq_nt6_table[(uint32_t)(tstr[i])]; bn0 = b32->n;
        for (bi = bn0 - nm; bi < bn0; bi += 2) {
            m = b32->a[bi]; s = 0; e = (*cns).a[m].arc.nou;
            for (k = s; k < e; k++) {
                if((*cns).a[m].arc.a[k].v == 0 || (*cns).a[m].arc.a[k].v == 1) continue;
                if(av_bp && (*cns).a[m].arc.a[k].v >= (*cns).bb0 && (*cns).a[m].arc.a[k].v < (*cns).bb1) continue;///no backbone
                if((*cns).a[(*cns).a[m].arc.a[k].v].c == bp) {
                    kv_push(uint32_t, (*b32), (*cns).a[m].arc.a[k].v);
                    kv_push(uint32_t, (*b32), bi);

                    // if(s0 == 2863) {
                    //     fprintf(stderr, "***1***[M::%s] s::%u\tb32->n::%u\n", __func__, (*cns).a[m].arc.a[k].v, (uint32_t)b32->n);
                    // }
                    // if((i + 1) == tl) break;///quick end
                    // if(b32->n > max_trace) break;///redue the size of b32
                    if(((i + 1) == tl) || (b32->n > max_trace)) {
                        ff = 1; break;
                    }
                }
            }
            if(ff) break;
        }

        if(b32->n <= bn0) {///no node
            break;
        } else {
            nm = b32->n - bn0;
        } 
    }
    // fprintf(stderr, "[M::%s] b32->n::%u, nm::%u\n", __func__, (uint32_t)b32->n, nm);

    // if(s0 == 2863) {
    //     fprintf(stderr, "[M::%s] i::%u\tnm::%u\tb32->n::%u\n", __func__, i, nm, (uint32_t)b32->n);
    // }
    if(i > 0 && nm > 0) {
        (*rn) = b32->a[b32->n - nm];
        for (bi = b32->n - nm; b32->a[bi + 1] != ((uint32_t)-1); bi = b32->a[bi + 1]) {
            // fprintf(stderr, "[M::%s] bi::%u, p_bi::%u\n", __func__, bi, b32->a[bi + 1]);
            des = b32->a[bi]; src = b32->a[b32->a[bi + 1]]; bp = ((src!=s0)?(1):(plus0));
            m = add_cns_arc(cns, src, des, 1, bp); assert(m);
            m = add_cns_arc(cns, des, src, 0, bp); assert(m);
            (*cns).a[des].sc++;

            // if(s0 == 2863) {
            //     fprintf(stderr, "***2***[M::%s] src::%u\tdes::%u\n", __func__, src, des);
            // }
        }
    } else {
        i = 0;
    }

    return i;
} 

///[s, e)
uint32_t push_cns_c1(cns_gfa *cns, uint64_t s0, char *tstr, uint64_t tl, uint32_t plus0, asg32_v* b32, uint64_t max_trace, uint32_t rid)
{
    if(tl <= 0) return s0;
    uint32_t rr = plus0, k, re = s0; 
    k = trace_cns_bp(cns, s0, tstr, tl, b32, plus0, &re, max_trace, 1);
    if(k > 0) rr = 1;

    // if(s0 == 2863) {
    //     fprintf(stderr, "[M::%s] (%.*s)\ts0::%lu\tk::%u\ttl::%lu\tre::%u\n", __func__, tstr?((int)(tl)):0, tstr, s0, k, tl, re);
    // }
    //     fprintf(stderr, "[M::%s] s0::%u, s::%u\n", __func__, k-1, k);
        //     prt_cns_arc(cns, k-1, __func__); prt_cns_arc(cns, k, __func__);

    for (; k < tl; k++) {///the weight of (s0 -> tstr[0]) might be 0
        re = add_cns_arc_bp(cns, re, seq_nt6_table[(uint32_t)(tstr[k])], rr, rid, 1); rr = 1;
    }
    
    return re;
}

uint64_t append_cns_g(cns_gfa *cns, char *tstr, uint64_t tl, uint64_t qs, uint64_t qe, uint64_t cp, uint64_t cl, uint64_t pe, asg32_v* b32, uint64_t max_trace, uint32_t rid, int64_t insert_pos)
{
    // fprintf(stderr, ">q::[%lu, %lu)\n", qs, qe);
    uint64_t s0 = pe, plus0 = 1, ns = qs + cns->off, ne = qe + cns->off;
    if(pe == ((uint64_t)-1)) {
        if(qs > 0) {
            s0 = qs - 1 + cns->off;///just before node in backbone
        } else {
            s0 = 0;//beg
        }
        // plus0 = 0;
    }

    // fprintf(stderr, "+n_nodes::%u, tl::%lu, qs::%lu, qe::%lu\n", (uint32_t)cns->n, tl, qs, qe);

    if(cp == 0) {
        if((cl == 0) && (cp == 0) && (qs == qe) && (((int64_t)qe) == insert_pos)) {
            s0 = 0;//beg
            ns = ne = 1;//end
            plus0 = 1;
        }
        // if(cp == 0) {
        //     fprintf(stderr, "cp::%lu, cl::%lu, qs::%lu, qe::%lu, s0::%lu, ns::%lu, ne::%lu, plus0::%lu\n", cp, cl, qs, qe, s0, ns, ne, plus0);
        // }
        return push_cns_c0(cns, s0, ns, ne, plus0, rid);
    } else if(cp == 1 || cp == 2) { ///cp == 2: more y -> insertion
        return push_cns_c1(cns, s0, tstr, tl, plus0, b32, max_trace, rid);
    } else { ///more x -> do nothing
        return s0;
    }
}

char *get_sub_seq(rr_seq_t *ssq, uint64_t s, uint64_t e)
{
    if(s >= e) return NULL;

    if(s >= ssq->s && e <= ssq->e) {
        return ssq->tu->seq + s - ssq->s;
    }

    uint64_t l = e - s;
    if(ssq->s >= ssq->e) {
        if(l < ssq->n0) l = ssq->n0;
    } else {
        if(l < ssq->n1) l = ssq->n1;
    }

    ssq->s = s; ssq->e = s + l;
    if(ssq->e > Get_READ_LENGTH((*(ssq->rref)), ssq->id)) {
        ssq->e = Get_READ_LENGTH((*(ssq->rref)), ssq->id);
        l = ssq->e - ssq->s;
    }
    UC_Read_resize((*(ssq->tu)), ((int64_t)l)); 
    recover_UC_Read_sub_region(ssq->tu->seq, ssq->s, l, ssq->rev, ssq->rref, ssq->id);
    return ssq->tu->seq;
}


///[s, e)
uint32_t extract_sub_cigar_cns(overlap_region *z, int64_t s, int64_t e, int64_t iws, int64_t iwe, int64_t s_end, rr_seq_t *ssq, ul_ov_t *p, cns_gfa *cns, asg32_v* b32, uint64_t max_trace, uint32_t rid)
{
    // if(s == 10539 && e == 10760) {
        // fprintf(stderr, "\n>>>>>>iw::[%ld, %ld)\tw::[%ld, %ld)\tox::[%u, %u)<<<<<<\n", iws, iwe, s, e, z->x_pos_s, z->x_pos_e + 1);
    // }
    
    int64_t wk = ovlp_cur_wid(*p), xk = ovlp_cur_xoff(*p), yk = ovlp_cur_yoff(*p), ck = ovlp_cur_coff(*p), os, oe, ots, ote, ol, insert_pos = ((iws == iwe)? (0): (-1));
    bit_extz_t ez; int64_t bd = ovlp_bd(*p), s0, e0, ii[2], it[2]; uint64_t pe = (uint64_t)-1;
    s0 = ((int64_t)(z->w_list.a[wk].x_start)) + bd; 
    e0 = ((int64_t)(z->w_list.a[wk].x_end)) + 1 - bd;
    if(s < s0) s = s0; if(e > e0) e = e0;///exclude boundary
    if(s > e) return -1;///it is possible s == e
    os = MAX(s, s0); oe = MIN(e, e0);
    if(oe < os) return -1;///it is possible os == oe

    set_bit_extz_t(ez, (*z), wk);
    if(!ez.cigar.n) return -1;
    int64_t cn = ez.cigar.n; uint16_t op; int64_t ws, we, wts, wte, ovlp; char *tstr;
    if((ck < 0) || (ck > cn)) {//(*ck) == cn is allowed
        ck = 0; xk = ez.ts; yk = ez.ps;
    }

    while (ck > 0 && xk >= s) {///x -> t; y -> p; first insertion and then match/mismatch
        --ck;
        op = ez.cigar.a[ck]>>14;
        if(op!=2) xk -= (ez.cigar.a[ck]&(0x3fff));
        if(op!=3) yk -= (ez.cigar.a[ck]&(0x3fff)); 
    }

    if(s_end == 0 && s == iws) s_end = 0;
    else s_end = 1;

    // if(s_end == 0 || s != iws) {///do not conside the insertion before s
    //     while (ck < cn && xk < s) {
    //     }
    // }

    // if(s == 10539 && e == 10760) {
        // fprintf(stderr, "ck::%ld, cn::%ld, xk::%ld, yk::%ld\n", ck, cn, xk, yk);
    // }

    // char cm[4]; cm[0] = 'M'; cm[1] = 'S'; cm[2] = 'I'; cm[3] = 'D';  
    //some cigar will span s or e
    ii[0] = ii[1] = it[0] = it[1] = -1; 
    ssq->s = ssq->e = 0; ssq->n0 = e - s; 
    ssq->id = z->y_id; ssq->rev = z->y_pos_strand;
    if(ssq->n0 == 0) {ssq->n0 = ssq->n1;} 

    while (ck < cn && xk < e) {//[s, e)
        ws = xk; wts = yk;
        op = ez.cigar.a[ck]>>14; ol = (ez.cigar.a[ck]&(0x3fff));
        
        for (ck++; (ck < cn) && (op == (ez.cigar.a[ck]>>14)); ck++) {
            ol += (ez.cigar.a[ck]&(0x3fff));
        }
        ///op == 3: -> x; op == 2: -> y;
        if(op!=2) xk += ol;
        if(op!=3) yk += ol; 
        we = xk; wte = yk;

        // fprintf(stderr, "ck::%ld, cn::%ld, op::%u, ol::%ld\n", ck, cn, op, ol);

        os = MAX(s, ws); oe = MIN(e, we);
        ovlp = ((oe>os)? (oe-os):0);

        if(s == e) {///insertion in comparsion with the reference
            if(op != 0 || ws >= s || we <= e || e != iwe || s != iws) continue;///must be a match
        } else {
            if(op != 2) {
                if(!ovlp) continue;
            } else {///ws == we
                if(ws < s || ws >= e) continue;
            }
        }
        
        if((s_end == 0) && (op == 2) && (ws == s)) continue;///skip the insertion just before s
        
        

        if(op < 2) {
            ots = os - ws + wts; ote = oe - ws + wts;
        } else {///op == 2: more y; p == 3: more x
            ots = wts; ote = wte;
        } 



        if(ii[0] == -1) {
            ii[0] = os; it[0] = ots;
            // if(op < 2) {
            //     it[0] = os - ws + wts;
            // } else {///op == 2: more y; p == 3: more x
            //     it[0] = wts;
            // } 
        }

        ii[1] = oe; it[1] = ote;
        // if(op < 2) {
        //     it[1] = oe - ws + wts;
        // } else {///op == 2: more y; p == 3: more x
        //     it[1] = wte;
        // } 




        if(op != 2) ol = oe-os;

        tstr = NULL;
        if(op != 0) tstr = get_sub_seq(ssq, ots, ote);

        // if(s == 10539 && e == 10760) {
            // fprintf(stderr, "+0-%ld%c(%.*s)\tpe::%lu", ol, cm[op], tstr?((int)(ote - ots)):0, tstr, pe);
        // }
        // fprintf(stderr, ">q::[%ld, %ld)\n", ws, we);
        // fprintf(stderr, "%ld%c(%.*s)", ol, cm[op], tstr?((int)(ote - ots)):0, tstr);

        pe = append_cns_g(cns, tstr, ote - ots, os - iws, oe - iws, op, ol, pe, b32, max_trace, rid, insert_pos);

        // if(s == 10539 && e == 10760) {
            // fprintf(stderr, "+1-pe::%lu\n", pe);
        // }
    }

    while (ck < cn && xk <= e) {//[s, e)
        ws = xk; wts = yk;
        op = ez.cigar.a[ck]>>14; ol = (ez.cigar.a[ck]&(0x3fff));
        if(op != 2) break;

        for (ck++; (ck < cn) && (op == (ez.cigar.a[ck]>>14)); ck++) {
            ol += (ez.cigar.a[ck]&(0x3fff));
        }
        yk += ol;//yk += (ez.cigar.a[ck]&(0x3fff)); 
        we = xk; wte = yk;


        if(ws >= s && ws <= e) {
            ots = wts; ote = wte;
            if(ii[0] == -1) {
                ii[0] = ws; it[0] = ots;
            }
            ii[1] = we; it[1] = ote;
            

            tstr = NULL;
            if(op != 0) tstr = get_sub_seq(ssq, ots, ote);
            // if(s == 10539 && e == 10760) {
                // fprintf(stderr, "-0-%ld%c(%.*s)\tpe::%lu", ol, cm[op], tstr?((int)(ote - ots)):0, tstr, pe);
            // }
            // fprintf(stderr, "%ld%c(%.*s)", ol, cm[op], tstr?((int)(ote - ots)):0, tstr);

            pe = append_cns_g(cns, tstr, ote - ots, ws - iws, we - iws, op, ol, pe, b32, max_trace, rid, insert_pos);

            // if(s == 10539 && e == 10760) {
                // fprintf(stderr, "-1-pe::%lu\n", pe);
            // }
        }
    }

    // if(s == 10539 && e == 10760) {
        // fprintf(stderr, "\tx::[%ld, %ld)\ty::[%ld, %ld)\tiw::[%ld, %ld)\n", ii[0], ii[1], it[0], it[1], iws, iwe);
    // }

    // prt_cns_arc(cns, 0, __func__);
    // prt_cns_arc(cns, 1, __func__);
    if(ii[1] == -1) return -1;///it is possible when s == e and the cigar here is not a match

    uint64_t ae = 1;
    if((ii[1] == iwe)) {
        ae = 1;///end node
    } else {
        ae = ii[1] + cns->off - iws;
    }

    // if(s == 10539 && e == 10760) {
        // fprintf(stderr, "pe::%lu, ae::%lu, n_nodes::%u, ii[1]::%ld\n", pe, ae, (uint32_t)cns->n, ii[1]);
    // }

    if(pe == ((uint64_t)-1)) pe = 0;///start node

    ///if iws == iwe and the cigar is a match, pe will be equal to ae
    if(pe != ae) {
        if(!add_cns_arc(cns, pe, ae, 1, /**ae==1?1:0**/1)) {
            insert_cns_arc(cns, pe, ae, 1, /**ae==1?1:0**/1, rid); 
            insert_cns_arc(cns, ae, pe, 0, /**ae==1?1:0**/1, rid);
        } else {
            add_cns_arc(cns, ae, pe, 0, /**ae==1?1:0**/1);
        }
    }    


    // prt_cns_arc(cns, 0, __func__);
    // prt_cns_arc(cns, 1, __func__);


    // if(s == 10539 && e == 10760) {
    //     fprintf(stderr, "-end-pe::%lu, ae::%lu\n", pe, ae);
    // }

    ovlp_cur_xoff(*p) = xk; ovlp_cur_yoff(*p) = yk; ovlp_cur_coff(*p) = ck; ovlp_cur_ylen(*p) = 0;
    
    return 1;
}


uint64_t iter_cc_idx_t(overlap_region* ol, cc_idx_t *z, int64_t s, int64_t e, uint64_t is_reduce, uint64_t is_insert, uint64_t **ra)
{
    int64_t rm_n, q[2], os, oe; ul_ov_t *cp; uint64_t m; *ra = NULL;

    // if(s == 15816 && e == 15819) {
    //     fprintf(stderr, "[M::%s] is_reduce::%lu\n", __func__, is_reduce);
    // }
    if(z->ru == 0) {
        if(is_reduce) {
            for (m = rm_n = z->srt_n; m < z->idx->n; m++) {
                cp = &(z->c_idx[z->idx->a[m]]);
                // if(s == 15816 && e == 15819) {
                //     fprintf(stderr, "-0-[M::%s] ii::%lu, ii0::%ld\n", __func__, z->idx->a[m], z->i0);
                // }
                
                q[0] = ol[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_start+ovlp_bd(*cp);
                q[1] = ol[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_end+1-ovlp_bd(*cp);
                os = MAX(q[0], s); oe = MIN(q[1], e);
                if((oe > os) || ((is_insert) && (s == e) && (s >= q[0]) && (s <= q[1]))) {
                    z->idx->a[rm_n++] = z->idx->a[m];
                }
            }
            z->idx->n = rm_n;
        }

        for (; z->i < z->srt_n; ++z->i) {
            cp = &(z->c_idx[(uint32_t)z->idx->a[z->i]]);
            // if(s == 15816 && e == 15819) {
            //     fprintf(stderr, "-1-[M::%s] ii::%u, ii0::%ld\n", __func__, (uint32_t)z->idx->a[z->i], z->i0);
            // }
            q[0] = ol[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_start+ovlp_bd(*cp);
            q[1] = ol[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_end+1-ovlp_bd(*cp);
            if(q[0] > e) break;
            if((!is_insert) && (q[0] >= e)) break;
            os = MAX(q[0], s); oe = MIN(q[1], e);
            if((oe > os) || ((is_insert) && (s == e) && (s >= q[0]) && (s <= q[1]))) {
                kv_push(uint64_t, *(z->idx), ((uint32_t)z->idx->a[z->i]));
            }
        }
    } else {
        z->ru = 0;
    }

    (*ra) = z->idx->a + z->srt_n; 
    return z->idx->n - z->srt_n;
}

void debug_inter0(overlap_region* ol, ul_ov_t *c_idx, uint64_t *idx, int64_t idx_n, uint64_t *res, int64_t res_n, int64_t s, int64_t e, uint64_t is_insert, uint64_t is_hard_check, const char *cmd)
{
    ul_ov_t *cp; int64_t q[2], a_n = 0, i, k = 0, os, oe; 
    for (i = 0; i < idx_n; i++) {
        cp = &(c_idx[(uint32_t)idx[i]]);
        q[0] = ol[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_start+ovlp_bd(*cp);
        q[1] = ol[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_end+1-ovlp_bd(*cp);

        // fprintf(stderr, "%s[M::%s] tid::%u\t%.*s\twid::%u\tq::[%u, %u)\terr::%d\toerr::%u\n", cmd, __func__, ol[ovlp_id(*cp)].y_id, (int)Get_NAME_LENGTH(R_INF, ol[ovlp_id(*cp)].y_id), Get_NAME(R_INF, ol[ovlp_id(*cp)].y_id),
        // ovlp_cur_wid(*cp), ol[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_start, ol[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_end+1, ol[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].error, ol[ovlp_id(*cp)].non_homopolymer_errors);

        os = MAX(q[0], s); oe = MIN(q[1], e);
        if((oe > os) || ((is_insert) && (s == e) && (s >= q[0]) && (s <= q[1]))) {
            a_n++; 
            // if(!(((uint32_t)idx[i]) == res[k])) {
            //     fprintf(stderr, "[M::%s] a_n::%ld\tres_n::%ld\ts::%ld\te::%ld\ti::%ld\tk::%ld\n", __func__, a_n, res_n, s, e, i, k);
            // }
            if(is_hard_check) {
                assert(((uint32_t)idx[i]) == res[k++]);
            } else {
                for (; (k < res_n) && (((uint32_t)idx[i]) != res[k]); k++);
                assert(k < res_n);
            }
        }
    }
    // if(a_n != res_n) {
    //     fprintf(stderr, "[M::%s] a_n::%ld\tres_n::%ld\ts::%ld\te::%ld\tidx_n::%ld\n", __func__, a_n, res_n, s, e, idx_n);
    // }
    if(is_hard_check) {
        assert(a_n == res_n);
    } else {
        assert(a_n <= res_n);
    }
}

void prt_cigar0(uint64_t in, int64_t len)
{
    int64_t k; uint64_t mp;
    char cm[4]; cm[0] = 'M'; cm[1] = 'S'; cm[2] = 'I'; cm[3] = 'D';  
    for (k = 0; k < len; k++) {
        mp = len - 1 - k; mp <<= 1;
        fprintf(stderr, "%c", cm[(in >> mp)&3]);
    }
    fprintf(stderr, "\n");
}

void prt_bp0(uint64_t in, int64_t len)
{
    int64_t k; uint64_t mp;
    char cm[4]; cm[0] = 'A'; cm[1] = 'C'; cm[2] = 'G'; cm[3] = 'T';  
    for (k = 0; k < len; k++) {
        mp = len - 1 - k; mp <<= 1;
        fprintf(stderr, "%c", cm[(in >> mp)&3]);
    }
    fprintf(stderr, "\n");
}

uint64_t cns_gen0(overlap_region* ol, All_reads *rref, uint64_t s, uint64_t e, uint64_t ql, UC_Read* tu, cc_idx_t *idx, uint64_t occ_tot, double occ_max, asg32_v* b32, uint32_t *rc)
{
    if(e > s + simp_vote_len) return 0;///too long

    uint64_t *id_a = NULL, id_n, an = 0, oc[2]; b32->n = 0; uint32_t m, *a = NULL;
    id_n = iter_cc_idx_t(ol, idx, s, e, idx->rr, ((s==e)?1:0), &id_a);
    // debug_inter0(ol, idx->c_idx, idx->idx->a + idx->i0, idx->srt_n - idx->i0, id_a, id_n, s, e, ((s==e)?1:0), 0, "-1-");
    uint64_t k, l, q[2], os, oe; ul_ov_t *p; overlap_region *z; idx->rr = 0;
    // fprintf(stderr, "[M::%s] [%lu, %lu) id_n::%lu\n", __func__, s, e, id_n);
    for (k = 0; k < id_n; k++) {
        p = &(idx->c_idx[id_a[k]]); z = &(ol[ovlp_id(*p)]); 
        q[0] = z->w_list.a[ovlp_cur_wid(*p)].x_start+ovlp_bd(*p);
        q[1] = z->w_list.a[ovlp_cur_wid(*p)].x_end+1-ovlp_bd(*p);

        // if(s == 11851 && e == 11853) {
        //     fprintf(stderr, "[M::%s] tid::%u\t%.*s\twid::%u\tq::[%u, %u)\terr::%d\toerr::%u\n", __func__, ol[ovlp_id(*p)].y_id, (int)Get_NAME_LENGTH(R_INF, ol[ovlp_id(*p)].y_id), Get_NAME(R_INF, ol[ovlp_id(*p)].y_id),
        //         ovlp_cur_wid(*p), ol[ovlp_id(*p)].w_list.a[ovlp_cur_wid(*p)].x_start, ol[ovlp_id(*p)].w_list.a[ovlp_cur_wid(*p)].x_end+1, ol[ovlp_id(*p)].w_list.a[ovlp_cur_wid(*p)].error, ol[ovlp_id(*p)].non_homopolymer_errors);
        // }
        
        if(q[1] <= e) idx->rr = 1;
        os = MAX(q[0], s); oe = MIN(q[1], e);
        // if((oe > os) || ((s == e) && (s >= q[0]) && (s <= q[1]))) {
        if((oe > os) || ((s == e) && (s > q[0]) && (s < q[1]))) {
        // if(oe >= os) {
            ///[-4-][-12-][-4-][-12-]
            ///[cigar_len][cigar][base_len][base]
            m = extract_sub_cigar_ii(z, ql, rref, os, oe, s, e, tu, p); an++;
            if(m != ((uint32_t)-1)) {///no gap in both sides
                kv_push(uint32_t, *b32, m);
            }
        }
    }

    oc[0] = b32->n; oc[1] = an + 1; //+1 for the reference read
    // if(s == 11851 && e == 11853) {
    //     fprintf(stderr, "-0-[M::%s] oc[0]::%lu, oc[1]::%lu\n", __func__, oc[0], oc[1]);
    // }
    if(((oc[0] > (oc[1]*occ_max)) && (oc[0] > (oc[1]-oc[0])) && (oc[1] >= occ_tot) && (oc[0] > 1))) {
        radix_sort_ec32(b32->a, b32->a+b32->n); an = 0;
        for (k = 1, l = 0; k <= b32->n; ++k) {
            if (k == b32->n || b32->a[k] != b32->a[l]) {
                if(k - l > an) {
                    an = k - l; a = b32->a + l;
                }
                l = k;
            }
        }      
        oc[0] = an;
        // fprintf(stderr, "-1-[M::%s] oc[0]::%lu, oc[1]::%lu\n", __func__, oc[0], oc[1]);
        if(((oc[0] > (oc[1]*occ_max)) && (oc[0] > (oc[1]-oc[0])) && (oc[1] >= occ_tot) && (oc[0] > 1))) {
            (*rc) = a[0];
            // prt_cigar0((a[0]<<4)>>20, a[0]>>28);
            // prt_bp0((a[0]<<20)>>20, (a[0]<<16)>>28);
            return 1;
        }
    }

    idx->ru = 1;
    return 0;
}

inline void gen_mm_cns_arc(cns_gfa *cns, uint32_t src, uint32_t des, uint32_t sc, uint32_t f)
{
    cns_t *av = &((*cns).a[src]), *aw;
    uint32_t vk, wk;
    for (vk = 0; vk < av->arc.nou; vk++) {///out-edge of src
        if((av->arc.a[vk].v != des) || (del_cns_arc((*av), vk))) continue;
        // av->arc.a[vk].f = 1; ///not sure if we should set these edges as visited
        av->arc.a[vk].f = f;
        av->arc.a[vk].sc += sc;

        aw = &((*cns).a[des]);
        for (wk = aw->arc.nou; wk < aw->arc.n; wk++) {///in-edge of des
            if((aw->arc.a[wk].v != src) || (del_cns_arc((*aw), wk))) continue;
            // aw->arc.a[wk].f = 1; ///not sure if we should set these edges as visited
            aw->arc.a[wk].f = f;
            aw->arc.a[wk].sc += sc; 
            break;
        }

        assert(wk < aw->arc.n);
        return;
    }

    cns_arc *p, t;
    ///src -> des
    kv_pushp(cns_arc, (*cns).a[src].arc, &p);
    p->sc = sc; p->v = des;
    // p->f = 1; ///not sure if we should set these edges as visited
    p->f = f; 
    ///ou-edge
    (*cns).a[src].arc.nou++;
    if((*cns).a[src].arc.nou < (*cns).a[src].arc.n) {
        t = (*cns).a[src].arc.a[(*cns).a[src].arc.nou-1]; 
        (*cns).a[src].arc.a[(*cns).a[src].arc.nou-1] = *p;
        *p = t;
    }

    ///src <- des; in-edge
    kv_pushp(cns_arc, (*cns).a[des].arc, &p);
    p->sc = sc; p->v = src;
    // p->f = 1; ///not sure if we should set these edges as visited
    p->f = f; 
}

void del_cns_g_nn(cns_gfa *cns, uint32_t v)
{
    uint32_t w, vk, wk; cns_t *av = &((*cns).a[v]), *aw = NULL;

    for (vk = 0; vk < av->arc.nou; vk++) {///out-edge of src
        if(del_cns_arc((*av), vk)) continue;
        w = av->arc.a[vk].v; av->arc.a[vk].v = CNS_DEL_E;

        aw = &((*cns).a[w]);
        for (wk = aw->arc.nou; wk < aw->arc.n; wk++) {///in-edge of des
            if((aw->arc.a[wk].v != v) || (del_cns_arc((*aw), wk))) continue;
            aw->arc.a[wk].v = CNS_DEL_E; break;
        }
        assert(wk < aw->arc.n);
    }

    for (vk = av->arc.nou; vk < av->arc.n; vk++) {///in-edge of src
        if(del_cns_arc((*av), vk)) continue;
        w = av->arc.a[vk].v; av->arc.a[vk].v = CNS_DEL_E;

        aw = &((*cns).a[w]);
        for (wk = 0; wk < aw->arc.nou; wk++) {///out-edge of des
            if((aw->arc.a[wk].v != v) || (del_cns_arc((*aw), wk))) continue;
            aw->arc.a[wk].v = CNS_DEL_E; break;
        }
        assert(wk < aw->arc.n);
    }


    cns->a[v].arc.n = cns->a[v].arc.nou = 0;
    cns->a[v].c = cns->a[v].f = 0; cns->a[v].sc = CNS_DEL_V;
}

void merge_cns_g_in(cns_gfa *cns, uint32_t v0, asg32_v* b32)
{
    cns_t *av, *aw; uint32_t v, bp, vk, wk, wka, w, wn, nn, mn, wh, mn_k[2];

    b32->n = 0;
    kv_push(uint32_t, *b32, v0);
    while (b32->n) {
        v = b32->a[--b32->n];
        if(del_cns_nn((*cns), v)) continue;

        av = &((*cns).a[v]);
        for (bp = 0; bp < 4; bp++) {
            //nn: number of node; wh: weight
            nn = wh = 0; mn = mn_k[0] = mn_k[1] = wka = (uint32_t)-1;
            for (vk = av->arc.nou; vk < av->arc.n; vk++) {///in-edge of v
                if(del_cns_arc((*av), vk)) continue;
                w = av->arc.a[vk].v; aw = &((*cns).a[w]);
                if(aw->c != bp) continue;
                if(w == cns->si || w == cns->ei) continue;

                for (wk = wn = 0; wk < aw->arc.nou; wk++) {///out-edge of w
                    if(del_cns_arc((*aw), wk)) continue;
                    wn++; wka = wk; if(wn > 1) break;
                }

                if(wn != 1) continue;
                
                assert(aw->arc.a[wka].v == v);

                ///deal with out-edge of w
                if(nn == 0) {
                    mn = w; mn_k[0] = vk; mn_k[1] = wka;
                    wh = av->arc.a[vk].sc;
                    ///not sure if we should set these edges as visited
                    // av->arc.a[vk].f = 1; aw->arc.a[wka].f = 1;
                } else {
                    wh += aw->arc.a[wka].sc;
                }

                ///deal with in-edge of w
                ///all edges to w, should be move to mn
                if(nn > 0) {///not sure if we should set these edges as visited; affect when nn == 0
                    for (wk = aw->arc.nou; wk < aw->arc.n; wk++) {
                        if(del_cns_arc((*aw), wk)) continue;
                        ///previously, aw->arc.a[wk].v -> w
                        ///currently, aw->arc.a[wk].v -> mn
                        /// if(nn == 0), then mn = w
                        gen_mm_cns_arc(cns, aw->arc.a[wk].v, mn, aw->arc.a[wk].sc/**(nn?(aw->arc.a[wk].sc):(0))**/, aw->arc.a[wk].f);///not sure if we should set these edges as visited
                    }
                }

                ///mn != w
                if(nn > 0) del_cns_g_nn(cns, w);

                nn++;
            }

            if(nn) {
                aw = &((*cns).a[mn]);
                av->arc.a[mn_k[0]].sc = aw->arc.a[mn_k[1]].sc = wh;
                // merge_cns_g_in(cns_gfa *cns, uint32_t v, asg32_v* b32)
                kv_push(uint32_t, *b32, mn);
            }
        }
    }
}

void merge_cns_g_ou(cns_gfa *cns, uint32_t v0, asg32_v* b32)
{
    cns_t *av, *aw; uint32_t v, bp, vk, wk, wka, w, wn, nn, mn, wh, mn_k[2];

    b32->n = 0;
    kv_push(uint32_t, *b32, v0);
    while (b32->n) {
        v = b32->a[--b32->n];
        if(del_cns_nn((*cns), v)) continue;
    
        av = &((*cns).a[v]);
        for (bp = 0; bp < 4; bp++) {
            //nn: number of node; wh: weight
            nn = wh = 0; mn = mn_k[0] = mn_k[1] = wka = (uint32_t)-1;
            for (vk = 0; vk < av->arc.nou; vk++) {///ou-edge of v
                if(del_cns_arc((*av), vk)) continue;
                w = av->arc.a[vk].v; aw = &((*cns).a[w]);
                if(aw->c != bp) continue;
                if(w == cns->si || w == cns->ei) continue;

                for (wk = aw->arc.nou, wn = 0; wk < aw->arc.n; wk++) {///in-edge of w
                    if(del_cns_arc((*aw), wk)) continue;
                    wn++; wka = wk; if(wn > 1) break;
                }

                if(wn != 1) continue;
                
                assert(aw->arc.a[wka].v == v);


                ///deal with in-edge of w
                if(nn == 0) {
                    mn = w; mn_k[0] = vk; mn_k[1] = wka;
                    wh = av->arc.a[vk].sc;
                    ///not sure if we should set these edges as visited
                    // av->arc.a[vk].f = 1; aw->arc.a[wka].f = 1;
                } else {
                    wh += aw->arc.a[wka].sc;
                }


                ///deal with ou-edge of w
                ///all edges from w, should be move to mn
                if(nn > 0) {///not sure if we should set these edges as visited; affect when nn == 0
                    for (wk = 0; wk < aw->arc.nou; wk++) {
                        if(del_cns_arc((*aw), wk)) continue;
                        ///previously, w -> aw->arc.a[wk].v
                        ///currently, mn -> aw->arc.a[wk].v 
                        /// if(nn == 0), then mn = w
                        gen_mm_cns_arc(cns, mn, aw->arc.a[wk].v, aw->arc.a[wk].sc/**(nn?(aw->arc.a[wk].sc):(0))**/, aw->arc.a[wk].f);///not sure if we should set these edges as visited
                    }
                }

                ///mn != w
                if(nn > 0) del_cns_g_nn(cns, w);

                nn++;
            }

            if(nn) {
                aw = &((*cns).a[mn]);
                // fprintf(stderr, "\n[M::%s] nn::%u, mn::%u\n", __func__, nn, mn);
                // fprintf(stderr, "[M::%s] vi::%u, vn::%u\n", __func__, mn_k[0], (uint32_t)av->arc.n);
                // fprintf(stderr, "[M::%s] wi::%u, wn::%u\n", __func__, mn_k[1], (uint32_t)aw->arc.n);

                av->arc.a[mn_k[0]].sc = aw->arc.a[mn_k[1]].sc = wh;
                // merge_cns_g_in(cns_gfa *cns, uint32_t v, asg32_v* b32)
                kv_push(uint32_t, *b32, mn);
            }
        }
    }
}

void refine_cns_g(cns_gfa *cns, asg32_v *b32)
{
    uint32_t v, w, vk, wk, *p = NULL; cns_t *av = NULL, *aw = NULL;
    kdq_clear(cns->q); 
    kdq_push(uint32_t, cns->q, cns->si); ///in-degree == 0

    while (1) {
        p = kdq_shift(uint32_t, cns->q);
        if(!p) break; v = *p;

        if(del_cns_nn((*cns), v)) continue;
        ///merge in
        merge_cns_g_in(cns, v, b32);
        ///merge out
        merge_cns_g_ou(cns, v, b32);

        av = &((*cns).a[v]);
        ///set arcs
        for (vk = 0; vk < av->arc.nou; vk++) {///out-edge of v
            if(del_cns_arc((*av), vk)) continue;
            if(av->arc.a[vk].f == 0) continue;

            av->arc.a[vk].f = 1; w = av->arc.a[vk].v; aw = &((*cns).a[w]);

            for (wk = aw->arc.nou; wk < aw->arc.n; wk++) {///in-edge of w
                if((aw->arc.a[wk].v != v) || (del_cns_arc((*aw), wk))) continue;
                aw->arc.a[wk].f = 1; break;
            }
        }
        ///set node
        av->f = 1;


        for (vk = 0; vk < av->arc.nou; vk++) {///out-edge of v
            if(del_cns_arc((*av), vk)) continue;
            w = av->arc.a[vk].v;

            aw = &((*cns).a[w]);
            for (wk = aw->arc.nou; wk < aw->arc.n; wk++) {///in-edge of w
                if((del_cns_arc((*aw), wk))) continue;
                if(aw->arc.a[wk].f) continue;///test arcs
                if(cns->a[aw->arc.a[wk].v].f) continue;///test node
                break;
            }
            if(wk >= aw->arc.n) {
                kdq_push(uint32_t, cns->q, w); ///in-degree == 0
            }
        }
    }
}

void gseq_cns_g(cns_gfa *cns, asg32_v *b32, uint32_t bl)
{
    b32->n = 0; kv_resize(uint32_t, *b32, cns->n);
    uint32_t v, vk, w, mme, mmn, mmk, mmw, *ii = b32->a, *p; cns_t *av;
    uint32_t bs = cns->off, be = bl + cns->off, sw;
    for (v = 0; v < cns->n; v++) {
        ii[v] = 0;///score
        if(del_cns_nn((*cns), v)) continue;
        cns->a[v].sc = 0; ///in-degree or prefix
        cns->a[v].f = 0;

        av = &((*cns).a[v]);
        for (vk = av->arc.nou; vk < av->arc.n; vk++) {///in-edge of v
            if(del_cns_arc((*av), vk)) continue;
            cns->a[v].sc++;///in-degree
        }
    }

    kdq_clear(cns->q); 
    kdq_push(uint32_t, cns->q, cns->si); ///in-degree == 0
    assert((cns->a[cns->si].sc == 0) && (!del_cns_nn((*cns), cns->si)));

    while (1) {
        p = kdq_shift(uint32_t, cns->q);
        if(!p) break; v = *p;

        if(del_cns_nn((*cns), v)) continue;
        assert(cns->a[v].sc == 0); ///in-degree == 0

        av = &((*cns).a[v]);
        for (vk = av->arc.nou, mme = mmn = mmw = 0, mmk = (uint32_t)-1; vk < av->arc.n; vk++) {///in-edge of v
            if(del_cns_arc((*av), vk)) continue;
            w = av->arc.a[vk].v;
            assert((*cns).a[w].f);
            sw = ((w >= bs && w < be)?1:0); ///backbone

            if((mmk == ((uint32_t)-1)) || (av->arc.a[vk].sc > mme) || 
                (((av->arc.a[vk].sc == mme) && (ii[w] > mmn))) ||
                (((av->arc.a[vk].sc == mme) && (ii[w] == mmn) && (sw == 1) && (mmw == 0)))) {
                mmk = vk; mme = av->arc.a[vk].sc; mmn = ii[w]; mmw = sw;
            }
        }

        ii[v] = mme + mmn; cns->a[v].f = 1;
        if(mmk != ((uint32_t)-1)) cns->a[v].sc = av->arc.a[mmk].v; 
        else cns->a[v].sc = v;

        for (vk = 0; vk < av->arc.nou; vk++) {///out-edge of v
            if(del_cns_arc((*av), vk)) continue;
            w = av->arc.a[vk].v;
            assert(cns->a[w].f == 0);
            assert(cns->a[w].sc);
            cns->a[w].sc--;
            if(cns->a[w].sc == 0) {
                kdq_push(uint32_t, cns->q, w); ///in-degree == 0
            }
        }

        // fprintf(stderr, "[M::%s] sc[%u]::%u\n", __func__, v, ii[v]);
    }

    for (v = cns->a[cns->ei].sc, b32->n = 0; v != cns->si; v = cns->a[v].sc) {
        // fprintf(stderr, "[M::%s] v::%u\n", __func__, v);
        kv_push(uint32_t, *b32, v);
    }
    assert(v == cns->si);
    
    mmn = b32->n; mmn >>= 1;
    for (vk = 0; vk < mmn; vk++) {
        v = b32->a[vk]; b32->a[vk] = b32->a[b32->n-vk-1]; b32->a[b32->n-vk-1] = v;
    }
}

uint64_t push_correct1(window_list *idx, window_list_alloc *res, cns_gfa *cns, asg32_v *rc, uint32_t bl)
{
    // fprintf(stderr, "[M::%s]\t", __func__);
    uint64_t nec = 0; uint32_t k, l, i, ff, sl, sk, bs = cns->off, be = bl + cns->off, bend = cns->off;///[bs, be)
    if(rc->n) {///it is possible that rc->n == 0, which means there is a deletion
        for (k = 1, l = 0; k <= rc->n; ++k) {
            ff = 0; sl = sk = 0;
            if(k == rc->n) {
                if(l < rc->n) sl = ((rc->a[l] >= bs && rc->a[l] < be)?1:0); 
                ff = 1;
            } else {
                sl = ((rc->a[l] >= bs && rc->a[l] < be)?1:0); 
                sk = ((rc->a[k] >= bs && rc->a[k] < be)?1:0);
                if(sl != sk) {
                    ff = 1;
                } else if((sl == 1) && ((rc->a[k] - rc->a[l]) != (k - l))) {
                    ff = 1;
                }
            }

            if(!ff) continue;

            if(sl) { ///match
                if(rc->a[l] > bend) {///deltetion [bend, rc->a[l])
                    push_trace_bp(((asg16_v *)(&(res->c))), 3, (uint16_t)-1, rc->a[l] - bend, ((idx->clen>0)?1:0)); 
                    idx->clen = res->c.n - idx->cidx; nec += rc->a[l] - bend;
                    // fprintf(stderr, "%uD", rc->a[l] - bend);
                }

                ///push match
                push_trace_bp(((asg16_v *)(&(res->c))), 0, (uint16_t)-1, rc->a[k-1] + 1 - rc->a[l], ((idx->clen>0)?1:0)); 
                idx->clen = res->c.n - idx->cidx;

                // fprintf(stderr, "%uM", rc->a[k-1] + 1 - rc->a[l]);

                bend = rc->a[k-1] + 1;
            } else { ///unmatch
                for (i = l; i < k; i++) {
                    push_trace_bp(((asg16_v *)(&(res->c))), 2, cns->a[rc->a[i]].c, 1, ((idx->clen>0)?1:0)); 
                    idx->clen = res->c.n - idx->cidx; nec++;
                    // fprintf(stderr, "I");
                }
            }
            l = k;
        }
    } 
    
    ///push remaining deletion
    if(be > bend) {
        push_trace_bp(((asg16_v *)(&(res->c))), 3, (uint16_t)-1, be - bend, ((idx->clen>0)?1:0)); 
        idx->clen = res->c.n - idx->cidx; nec += be - bend;
        // fprintf(stderr, "%uD", be - bend);
    }
    // fprintf(stderr, "\n");
    return nec;
}

uint64_t push_correct1_fhc_indel_exz(asg16_v *sc, int64_t sc0, window_list *idx, cns_gfa *cns, char *ostr, UC_Read* tu, bit_extz_t *exz, uint64_t gbeg, uint64_t c0, int64_t cl0, int64_t ok0)
{
    int64_t ck = sc->n, k, ok = 0, nk = 0, cn, cn0, nl, ol, diff, diff0, ml, ml0, e0 = 0; uint32_t on, f = 0, nec = 0; uint16_t bq, bt, op; 
    assert(c0 == 3);

    ///debug
    // char cm[4]; cm[0] = 'M'; cm[1] = 'S'; cm[2] = 'I'; cm[3] = 'D';

    // if(c0 != 2) ok += cl0;
    // if(c0 != 3) nk += cl0;
    ok += cl0; e0 += cl0;

    // fprintf(stderr, "\n%lu%c", cl0, cm[c0]);

    for (ck--; ck >= sc0/**0**/; ck--) {
        op = sc->a[ck]>>14;
        if(!op) break;

        if((op == 2) || (op == 3)) {
            on = sc->a[ck]&(0xfff);
        } else if(op == 1) {
            on = sc->a[ck]&(0x3ff);
        } else {
            on = sc->a[ck]&(0x3fff);
        }
        if(op != 2) ok += on;
        if(op != 3) nk += on;
        if(op != 0) e0 += on;
        if(c0 != op) f = 1;

        // fprintf(stderr, "%u%c(%c)", on, cm[op], "ACGT"[((sc->a[ck]>>12)&3)]);
    }
    cn0 = ck + 1; cn = sc->n;
    // fprintf(stderr, "\n");

    // fprintf(stderr, "+[M::%s] cn0::%ld, cn::%ld, ok::%ld, sc->n::%u, ok0::%ld, cl0::%ld\n", __func__, cn0, cn, ok, (uint32_t)sc->n, ok0, cl0);
    // f = 0;

    if((!f) || (!ok) || (!nk)) {
        for (k = 0, ck = ok0 + gbeg; k < cl0; k++, ck++) {
            push_trace_bp_f(sc, c0, cns->a[ck].c, (uint16_t)-1,  1, ((idx->clen>0)?1:0));
            idx->clen = sc->n - idx->cidx; nec++;
            // fprintf(stderr, "%c\n", "ACGT"[cns->a[ck].c]);
        }
    } else {
        char *oseq = NULL, *nseq = NULL; int64_t wo[2], wn[2];
        UC_Read_resize((*tu), nk); nseq = tu->seq; wo[0] = wo[1] = wn[0] = wn[1] = 0; 
        // if(c0 != 2) ok0 += cl0;
        ok0 += cl0;

        ok0 -= ok; 
        // if(!(ok0 >= 0)) {
        //     fprintf(stderr, "+[M::%s] rid::%u, cn0::%ld, cn::%ld, ok0::%ld, ok::%ld, sc->n::%u\n", __func__, rid, cn0, cn, ok0, ok, (uint32_t)sc->n);
        // }
        assert(ok0 >= 0);
        
        oseq = ostr + ok0;
        ol = ok; nl = nk; 
        ck = cn0; ok = nk = 0;
        while (ck < cn) {
            wo[0] = ok; wn[0] = nk;
            ck = pop_trace_bp_f(sc, ck, &op, &bq, &bt, &on);
            if(op != 2) ok += on;
            if(op != 3) nk += on;
            wo[1] = ok; wn[1] = nk;

            if(op == 0) {
                memcpy(nseq + wn[0], oseq + wo[0], (wo[1]-wo[0])*sizeof((*nseq)));
            } else if(op == 1 || op == 2) {
                for (k = wn[0]; k < wn[1]; k++) nseq[k] = s_H[bt];
            }
            // fprintf(stderr, "[M::%s] ck::%ld, wo::[%ld, %ld), wn::[%ld, %ld)\n", __func__, ck, wo[0], wo[1], wn[0], wn[1]);
        }
        // fprintf(stderr, "[M::%s] qstr::%.*s*\n", __func__, (int32_t)ol, oseq);
        // fprintf(stderr, "[M::%s] tstr::%.*s*\n", __func__, (int32_t)nl, nseq);
        
        assert(wo[1] + cl0 == ol); 
        assert(wn[1] == nl); ///since c0 must be 3
        // memcpy(nseq + wn[1], oseq + wo[1], cl0*sizeof((*nseq)));
        if(nl == ol && nl == 1) {
            sc->n = cn0;
            if(oseq[0] == nseq[0]) {
                push_trace_bp_f(sc, 0, (uint16_t)-1, (uint16_t)-1, 1, ((idx->clen>0)?1:0)); 
                idx->clen = sc->n - idx->cidx; 
            } else {
                push_trace_bp_f(sc, 1, seq_nt6_table[(uint32_t)(oseq[0])], seq_nt6_table[(uint32_t)(nseq[0])], 1, ((idx->clen>0)?1:0)); 
                idx->clen = sc->n - idx->cidx; nec++;
            }
        } else {
            ml = MAX(ol, nl); f = 0;
            
            diff = 31;
            if(diff > ml) diff = ml; 
            diff0 = diff; clear_align(*exz);
            cal_exz_global(nseq, nl, oseq, ol, diff, exz);
            if(is_align(*exz)) f = 1; 

            if(!f) {
                diff = 63;
                if(diff > ml) diff = ml; 
                if(diff > diff0) {
                    diff0 = diff; clear_align(*exz);
                    cal_exz_global(nseq, nl, oseq, ol, diff, exz);
                    if(is_align(*exz)) f = 1; 
                }
            }

            // fprintf(stderr, "[M::%s] f::%u, exz->err::%d, e0::%ld\n", __func__, f, exz->err, e0);

            if(f && exz->err < e0) {
                sc->n = cn0;
                cn = exz->cigar.n; ok = nk = 0;
                for (ck = 0; ck < cn;) {
                    wo[0] = ok; wn[0] = nk;
                    ck = pop_trace(&(exz->cigar), ck, &op, &on);
                    if(op!=2) ok += on;
                    if(op!=3) nk += on; 
                    wo[1] = ok; wn[1] = nk;

                    // fprintf(stderr, "%u%c(", on, cm[op]);

                    if(op == 0) {
                        push_trace_bp_f(sc, op, (uint16_t)-1, (uint16_t)-1, on, ((idx->clen>0)?1:0)); 
                        idx->clen = sc->n - idx->cidx; 
                    } else if(op == 1) {
                        for (k = 0; k < on; k++) {
                            push_trace_bp_f(sc, op, seq_nt6_table[(uint32_t)(oseq[wo[0]+k])], seq_nt6_table[(uint32_t)(nseq[wn[0]+k])], 1, ((idx->clen>0)?1:0)); 
                            idx->clen = sc->n - idx->cidx; nec++;
                            // fprintf(stderr, "<%c|%c>)", oseq[wo[0]+k], nseq[wn[0]+k]);
                        }
                    } else if(op == 2) {
                        for (k = 0; k < on; k++) {
                            push_trace_bp_f(sc, op, (uint16_t)-1, seq_nt6_table[(uint32_t)(nseq[wn[0]+k])], 1, ((idx->clen>0)?1:0)); 
                            idx->clen = sc->n - idx->cidx; nec++;
                            // fprintf(stderr, "<|%c>)", nseq[wn[0]+k]);
                        }
                    } else if(op == 3) {
                        for (k = 0; k < on; k++) {
                            push_trace_bp_f(sc, op, seq_nt6_table[(uint32_t)(oseq[wo[0]+k])], (uint16_t)-1, 1, ((idx->clen>0)?1:0)); 
                            idx->clen = sc->n - idx->cidx; nec++;
                            // fprintf(stderr, "<%c|>)", oseq[wo[0]+k]);
                        }
                    }
                    // fprintf(stderr, ")");
                }
                // fprintf(stderr, "\n");
            } else {
                if(ml < e0) {
                    sc->n = cn0;
                    ml0 = MIN(ol, nl); op = 1;
                    for (k = 0; k < ml0; k++) {
                        push_trace_bp_f(sc, op, seq_nt6_table[(uint32_t)(oseq[k])], seq_nt6_table[(uint32_t)(nseq[k])], 1, ((idx->clen>0)?1:0)); 
                        idx->clen = sc->n - idx->cidx; nec++;
                    }
                    
                    if(ol > ml0) {///op = 3
                        for (k = ml0, op = 3; k < ol; k++) {
                            push_trace_bp_f(sc, op, seq_nt6_table[(uint32_t)(oseq[k])], (uint16_t)-1, 1, ((idx->clen>0)?1:0)); 
                            idx->clen = sc->n - idx->cidx; nec++;
                        }
                    } else if(nl > ml0) {///op = 2
                        for (k = ml0, op = 2; k < nl; k++) {
                            push_trace_bp_f(sc, op, (uint16_t)-1, seq_nt6_table[(uint32_t)(nseq[k])], 1, ((idx->clen>0)?1:0)); 
                            idx->clen = sc->n - idx->cidx; nec++;
                        }
                    }
                } else {
                    for (k = 0, ck = ok0 + gbeg; k < cl0; k++, ck++) {
                        push_trace_bp_f(sc, c0, cns->a[ck].c, (uint16_t)-1,  1, ((idx->clen>0)?1:0));
                        idx->clen = sc->n - idx->cidx; nec++;
                    }
                }
            }
        }
    }

    // fprintf(stderr, "-[M::%s] cn0::%ld, cn::%ld, ok::%ld, sc->n::%u\n", __func__, cn0, cn, ok, (uint32_t)sc->n);
    return nec;
}

uint64_t push_correct1_fhc(window_list *idx, window_list_alloc *res, cns_gfa *cns, char* qstr, UC_Read* tu, bit_extz_t *exz, asg32_v *rc, uint32_t bl, uint32_t rid)
{
    // fprintf(stderr, "[M::%s]\trc->n::%u\tbl::%u\n", __func__, (uint32_t)rc->n, bl);
    uint64_t nec = 0; uint32_t k, l, i, ff, sl, sk, bs = cns->off, be = bl + cns->off, bend = cns->off, is_i = 0, sc0 = res->c.n;///[bs, be)
    if(rc->n) {///it is possible that rc->n == 0, which means there is a deletion
        for (k = 1, l = 0; k <= rc->n; ++k) {
            ff = 0; sl = sk = 0;
            if(k == rc->n) {
                if(l < rc->n) sl = ((rc->a[l] >= bs && rc->a[l] < be)?1:0); 
                ff = 1;
            } else {
                sl = ((rc->a[l] >= bs && rc->a[l] < be)?1:0); 
                sk = ((rc->a[k] >= bs && rc->a[k] < be)?1:0);
                if(sl != sk) {
                    ff = 1;
                } else if((sl == 1) && ((rc->a[k] - rc->a[l]) != (k - l))) {
                    ff = 1;
                }
            }

            if(!ff) continue;

            if(sl) { ///match
                if(rc->a[l] > bend) {///deltetion [bend, rc->a[l])
                    if(is_i && exz) {
                        nec += push_correct1_fhc_indel_exz(((asg16_v *)(&(res->c))), sc0, idx, cns, qstr, tu, exz, cns->off, 3, rc->a[l] - bend, bend-cns->off);
                    } else {
                        for (i = bend; i < rc->a[l]; i++) {
                            push_trace_bp_f(((asg16_v *)(&(res->c))), 3, cns->a[i].c, (uint16_t)-1,  1, ((idx->clen>0)?1:0));
                            idx->clen = res->c.n - idx->cidx; nec++;
                        }
                    }
                }

                ///push match
                push_trace_bp_f(((asg16_v *)(&(res->c))), 0, (uint16_t)-1, (uint16_t)-1, rc->a[k-1] + 1 - rc->a[l], ((idx->clen>0)?1:0)); 
                idx->clen = res->c.n - idx->cidx;

                bend = rc->a[k-1] + 1; is_i = 0;
            } else { ///unmatch
                for (i = l; i < k; i++) {
                    push_trace_bp_f(((asg16_v *)(&(res->c))), 2, (uint16_t)-1, cns->a[rc->a[i]].c, 1, ((idx->clen>0)?1:0)); 
                    idx->clen = res->c.n - idx->cidx; nec++;
                }
                is_i = 1;
            }
            l = k;
        }
    } 
    
    ///push remaining deletion
    if(be > bend) {
        if(is_i && exz) {
            nec += push_correct1_fhc_indel_exz(((asg16_v *)(&(res->c))), sc0, idx, cns, qstr, tu, exz, cns->off, 3, be - bend, bend-cns->off);
        } else {
            for (i = bend; i < be; i++) {
                push_trace_bp_f(((asg16_v *)(&(res->c))), 3, cns->a[i].c, (uint16_t)-1, 1, ((idx->clen>0)?1:0));
                idx->clen = res->c.n - idx->cidx; nec++;
            }
        }
    }
    // fprintf(stderr, "\n");
    return nec;
}

///no e_end since e always covers end; s may not have end
uint64_t cns_gen_full0(overlap_region* ol, All_reads *rref, uint64_t s, uint64_t e, uint64_t s_end, char* qstr, UC_Read* tu, bit_extz_t *exz, cc_idx_t *idx, uint64_t occ_tot, double occ_max, asg32_v* b32, cns_gfa *cns, uint64_t max_trace, window_list *ridx, window_list_alloc *res, uint32_t rid)
{
    init_cns_g(cns, qstr + s, e - s, rid);
    // if(e -s > 100) {
        // fprintf(stderr, "[M::%s]::[%lu, %lu), cns->n::%u, qstr::%.*s\n", __func__, s, e, (uint32_t)cns->n, (int)(e-s), qstr+s);
    // }
    // fprintf(stderr, "[M::%s]::[%lu, %lu)\n", __func__, s, e);
    
    uint64_t *id_a = NULL, id_n, nec = 0; b32->n = 0; 
    rr_seq_t ssq; ssq.rref = rref; ssq.tu = tu; ssq.s = ssq.e = 0; ssq.n0 = ssq.n1 = 32; ssq.id = ssq.rev = 0;
    id_n = iter_cc_idx_t(ol, idx, s, e, idx->rr, ((s==e)?1:0), &id_a);
    // debug_inter0(ol, idx->c_idx, idx->idx->a + idx->i0, idx->srt_n - idx->i0, id_a, id_n, s, e, ((s==e)?1:0), 0, "-1-");
    uint64_t k, q[2], os, oe; ul_ov_t *p; overlap_region *z; idx->rr = 0;
    // fprintf(stderr, "[M::%s] [%lu, %lu) id_n::%lu\n", __func__, s, e, id_n);
    for (k = 0; k < id_n; k++) {
        p = &(idx->c_idx[id_a[k]]); z = &(ol[ovlp_id(*p)]); 
        q[0] = z->w_list.a[ovlp_cur_wid(*p)].x_start+ovlp_bd(*p);
        q[1] = z->w_list.a[ovlp_cur_wid(*p)].x_end+1-ovlp_bd(*p);

        // fprintf(stderr, "[M::%s] tid::%u\t%.*s\twid::%u\tq::[%u, %u)\terr::%d\toerr::%u\n", __func__, ol[ovlp_id(*p)].y_id, (int)Get_NAME_LENGTH(R_INF, ol[ovlp_id(*p)].y_id), Get_NAME(R_INF, ol[ovlp_id(*p)].y_id),
        // ovlp_cur_wid(*p), ol[ovlp_id(*p)].w_list.a[ovlp_cur_wid(*p)].x_start, ol[ovlp_id(*p)].w_list.a[ovlp_cur_wid(*p)].x_end+1, ol[ovlp_id(*p)].w_list.a[ovlp_cur_wid(*p)].error, ol[ovlp_id(*p)].non_homopolymer_errors);

        if(q[1] <= e) idx->rr = 1;
        os = MAX(q[0], s); oe = MIN(q[1], e);
        // if((oe > os) || ((s == e) && (s >= q[0]) && (s <= q[1]))) {
        if((oe > os) || ((s == e) && (s > q[0]) && (s < q[1]))) {
        // if(oe >= os) {
            ///[-4-][-12-][-4-][-12-]
            ///[cigar_len][cigar][base_len][base]
            extract_sub_cigar_cns(z, os, oe, s, e, s_end, &ssq, p, cns, b32, max_trace, rid); 
            // if(m != ((uint32_t)-1)) {///no gap in both sides
            //     kv_push(uint32_t, *b32, m);
            // }
        }
    }

    // fprintf(stderr, "-2-[M::%s] cns->n::%u\n", __func__, (uint32_t)cns->n);

    // return;

    refine_cns_g(cns, b32);

    // fprintf(stderr, "-3-[M::%s] cns->n::%u\n", __func__, (uint32_t)cns->n);

    gseq_cns_g(cns, b32, e - s);

    // fprintf(stderr, "-4-[M::%s] cns->n::%u\n", __func__, (uint32_t)cns->n);

    // nec += push_correct1(ridx, res, cns, b32, e - s);
    nec += push_correct1_fhc(ridx, res, cns, qstr + s, tu, exz, b32, e - s, rid);

    // fprintf(stderr, "-5-[M::%s] cns->n::%u\n", __func__, (uint32_t)cns->n);
    return nec;
}

uint64_t cns_gen_full(overlap_region* ol, All_reads *rref, uint64_t s0, uint64_t e0, uint64_t wl, char* qstr, UC_Read* tu, bit_extz_t *exz, cc_idx_t *idx, uint64_t occ_tot, double occ_max, asg32_v* b32, cns_gfa *cns, uint64_t max_trace, window_list *ridx, window_list_alloc *res, uint32_t rid)
{
    uint64_t nec = 0;
    if(e0 - s0 <= wl) {
        nec += cns_gen_full0(ol, rref, s0, e0, 1, qstr, tu, exz, idx, occ_tot, occ_max, b32, cns, max_trace, ridx, res, rid);
    } else {
        uint64_t s, e;
        s = s0; e = s0 + wl; e = ((e<=e0)?e:e0);
        for (; s < e0; ) {
            // rn = iter_cc_idx_t(ol->list, &ii_a, s, e, rr, 0, &ra);
            // debug_inter0(ol->list, ii_a.c_idx, ii_a.idx->a + ii_a.i0, ii_a.srt_n - ii_a.i0, ra, rn, s, e, 0, 1, "-0-");
            // rr = wcns_vote(ol->list, rref, qu->seq, ql, tu, ra, rn, s, e, ii_a.c_idx, &ii_b, occ_tot, occ_exact, aux_o, b32, cns, rid);
            nec += cns_gen_full0(ol, rref, s, e, (s==s0)?1:0, qstr, tu, exz, idx, occ_tot, occ_max, b32, cns, max_trace, ridx, res, rid);
            s += wl; e += wl; e = ((e<=e0)?e:e0);
        }
    }
    return nec;
}

uint64_t push_correct0(window_list *idx, window_list_alloc *res, uint32_t len0, uint32_t rc)
{
    uint64_t nec = 0;
    // fprintf(stderr, "[M::%s] NULL(idx)::%u\n", __func__, idx?0:1);
    if(len0 != ((uint32_t)-1)) {
        push_trace_bp(((asg16_v *)(&(res->c))), 0, (uint16_t)-1, len0, ((idx->clen>0)?1:0)); 
        idx->clen = res->c.n - idx->cidx;
    } else if(rc != ((uint32_t)-1)) {
        // fprintf(stderr, "[M::%s]\t", __func__);
        uint32_t cc = (rc<<4)>>20, cn = rc>>28, ck = 0, cs, cp;
        uint32_t bc = (rc<<20)>>20, bn = (rc<<16)>>28, bk = 0, bs, bp;
        ///debug
        // prt_cigar0(cc, cn);
        // prt_bp0(bc, bn);

        for (ck = 0; ck < cn; ck++) {
            cs = (cn-1-ck)<<1; cp = (cc>>cs)&3;

            bp = (uint32_t)-1;
            if(cp != 3) {///bp == 3: more x
                bs = (bn-1-bk)<<1; bp = (bc>>bs)&3; 
                bk++;
            }
            // fprintf(stderr, "[M::%s] cp::%u, bp::%u\n", __func__, cp, bp);

            push_trace_bp(((asg16_v *)(&(res->c))), cp, bp, 1, ((idx->clen>0)?1:0)); 
            idx->clen = res->c.n - idx->cidx;
            // fprintf(stderr, "%c", cm[(in >> mp)&3]);

            // fprintf(stderr, "%c", "MSID"[cp]);
            if(cp != 0) nec++;
        }

        // fprintf(stderr, "\n");
    }

    return nec;
}

uint64_t push_correct0_fhc(window_list *idx, window_list_alloc *res, uint32_t len0, uint32_t rc, char *qstr)
{
    uint64_t nec = 0;
    // fprintf(stderr, "[M::%s] NULL(idx)::%u\n", __func__, idx?0:1);
    if(len0 != ((uint32_t)-1)) {
        push_trace_bp_f(((asg16_v *)(&(res->c))), 0, (uint16_t)-1, (uint16_t)-1, len0, ((idx->clen>0)?1:0)); 
        idx->clen = res->c.n - idx->cidx;
    } else if(rc != ((uint32_t)-1)) {
        // fprintf(stderr, "[M::%s]\t", __func__);
        uint32_t cc = (rc<<4)>>20, cn = rc>>28, ck = 0, cs, cp;
        uint32_t bc = (rc<<20)>>20, bn = (rc<<16)>>28, btk = 0, bqk = 0, bs, bqp, btp;
        ///debug
        // prt_cigar0(cc, cn);
        // prt_bp0(bc, bn);

        for (ck = 0; ck < cn; ck++) {
            cs = (cn-1-ck)<<1; cp = (cc>>cs)&3;

            bqp = btp = (uint32_t)-1;
            if(cp != 3) {///bp == 3: more x
                bs = (bn-1-btk)<<1; btp = (bc>>bs)&3; btk++;
            }
            if(cp != 2) {///bp == 2: more y
                bqp = seq_nt6_table[(uint32_t)qstr[bqk]]; bqk++;
            }
            // fprintf(stderr, "[M::%s] cp::%u, btp::%u, bqp::%u\n", __func__, cp, btp, bqp);

            push_trace_bp_f(((asg16_v *)(&(res->c))), cp, bqp, btp, 1, ((idx->clen>0)?1:0)); 
            idx->clen = res->c.n - idx->cidx;
            // fprintf(stderr, "%c", cm[(in >> mp)&3]);

            // fprintf(stderr, "%c", "MSID"[cp]);
            if(cp != 0) nec++;
        }

        // fprintf(stderr, "\n");
    }

    return nec;
}

void output_cns_g(cns_gfa *cns, uint64_t s, uint64_t e)
{
    // char *gfa_id = NULL, p; MALLOC(gfa_id, Get_NAME_LENGTH(R_INF, qid) + 128); 
    // sprintf(gfa_id, "%.*s.%lu_%lu.cns.gfa", (int)Get_NAME_LENGTH(R_INF, qid), Get_NAME(R_INF, qid), s, e);
    char *gfa_id = NULL, p, f; MALLOC(gfa_id, 128); 
    sprintf(gfa_id, "ec.%lu_%lu.cns.gfa", s, e);

    FILE *fp = fopen(gfa_id, "w"); uint64_t k, z; char cm[4]; cm[0] = 'A'; cm[1] = 'C'; cm[2] = 'G'; cm[3] = 'T';  

    sprintf(gfa_id, "s_0_%c", cm[cns->a[0].c]);
    fprintf(fp, "S\t%s\t*\tLN:i:%u\trd:i:%c\n", gfa_id, 0/**cns->a[0].sc**/, cm[cns->a[0].c]);

    sprintf(gfa_id, "e_1_%c", cm[cns->a[1].c]);
    fprintf(fp, "S\t%s\t*\tLN:i:%u\trd:i:%c\n", gfa_id, 0/**cns->a[1].sc**/, cm[cns->a[1].c]);

    for (k = 2; k < cns->n; k++) {
        if(del_cns_nn((*cns), k)) continue;
        sprintf(gfa_id, "%c_%lu_%c", ((k<cns->bn)?'b':'n'), k, cm[cns->a[k].c]);
        fprintf(fp, "S\t%s\t*\tLN:i:%d\trd:i:%c\n", gfa_id, 0/**cns->a[k].sc**/, cm[cns->a[k].c]);
    }
    
    for (k = 0; k < cns->n; k++) {
        if(del_cns_nn((*cns), k)) continue;

        if(k == 0) {
            p = 's';
        } else if(k == 1) {
            p = 'e';
        } else if (k<cns->bn) {
            p = 'b';
        } else {
            p = 'n';
        }
        
        sprintf(gfa_id, "%c_%lu_%c", p, k, cm[cns->a[k].c]);
        for (z = 0; z < cns->a[k].arc.n; z++) {
            if(del_cns_arc((cns->a[k]), z)) continue;

            if(cns->a[k].arc.a[z].v == 0) {
                p = 's';
            } else if(cns->a[k].arc.a[z].v == 1) {
                p = 'e';
            } else if (cns->a[k].arc.a[z].v < cns->bn) {
                p = 'b';
            } else {
                p = 'n';
            }

            f = ((z < cns->a[k].arc.nou)?('+'):('-'));
            fprintf(fp, "L\t%s\t%c\t%c_%u_%c\t%c\t0M\tL1:i:%u\n", 
                    gfa_id, f, 
                    p, cns->a[k].arc.a[z].v, cm[cns->a[cns->a[k].arc.a[z].v].c], f, 
                    cns->a[k].arc.a[z].sc);
        }
    }

    fclose(fp); free(gfa_id);
}

uint32_t cal_cigar_xlen(overlap_region *in)
{
    assert(in->w_list.n == 1);
    uint32_t cn = in->w_list.a[0].clen; uint16_t *a = in->w_list.c.a + in->w_list.a[0].cidx;
    uint32_t ci, len, xk, yk; uint16_t c, b; asg16_v scc; scc.n = scc.m = cn; scc.a = a;

    ci = 0; xk = yk = 0;
    while (ci < cn) {
        ci = pop_trace_bp(&scc, ci, &c, &b, &len);
        if(c != 2) xk += len;
        if(c != 3) yk += len;
    }

    return xk;
}

uint32_t cal_cigar_xlen_fhc(overlap_region *in)
{
    assert(in->w_list.n == 1);
    uint32_t cn = in->w_list.a[0].clen; uint16_t *a = in->w_list.c.a + in->w_list.a[0].cidx;
    uint32_t ci, len, xk, yk; uint16_t c, bq, bt; asg16_v scc; scc.n = scc.m = cn; scc.a = a;

    ci = 0; xk = yk = 0;
    while (ci < cn) {
        ci = pop_trace_bp_f(&scc, ci, &c, &bq, &bt, &len);
        if(c != 2) xk += len;
        if(c != 3) yk += len;
    }

    return xk;
}

uint64_t push_cns_anchor(overlap_region* ol, All_reads *rref, uint64_t s, uint64_t e, char* qstr, uint64_t ql, UC_Read* tu, bit_extz_t *exz, cc_idx_t *idx, overlap_region *aux_o, uint64_t is_tail, uint64_t occ_tot, double occ_max, asg32_v* b32, cns_gfa *cns, uint32_t rid)
{
    if((!is_tail) && (s >= e)) return 0;//if s >= e && is_tail = 1, gen the cns of the last a few bases -> s == e == ql
    // fprintf(stderr, "\n****************[M::%s::M] [%lu, %lu)****************\n", __func__, s, e);
    window_list *p = NULL; uint64_t e0 = 0, nec = 0; uint32_t rc;
    if(aux_o->w_list.n > 0) {
        p = &(aux_o->w_list.a[aux_o->w_list.n-1]);
        e0 = p->x_end+1;
        ///make sure e > s
    }
    assert(s >= e0);
    if((s == e) && (is_tail == 1) && (s == e0)) return 0;///in this case, s == e == ql

    if(aux_o->w_list.n == 0) {
        kv_pushp(window_list, aux_o->w_list, &p);
        p->x_start = -1; p->x_end = -1;
        p->clen = p->cidx = 0;
    }

    if(((!is_tail) && (s > 0)) || ((is_tail) && (s > e0))) {
        // fprintf(stderr, ">>>>>>[M::%s::M] [%lu, %lu), ec::%u\n", __func__, e0, s, cal_cigar_xlen_fhc(aux_o));
        // fprintf(stderr, ">>>>>>[M::%s::M] [%lu, %lu)\n", __func__, e0, s);
        if (cns_gen0(ol, rref, e0, s, ql, tu, idx, occ_tot, occ_max, b32, &rc)) {
            if(p->x_start == -1 || p->x_end == -1) {///hasn't neem set
                p->x_start = e0; p->x_end = s-1;
            }
            // nec += push_correct0(p, &(aux_o->w_list), (uint32_t)-1, rc);
            nec += push_correct0_fhc(p, &(aux_o->w_list), (uint32_t)-1, rc, qstr + e0);
            // fprintf(stderr, "-0-[M::%s::M]\n", __func__);
        } else {
            nec += cns_gen_full(ol, rref, e0, s, cns->cns_g_wl, qstr, tu, exz, idx, occ_tot, occ_max, b32, cns, ql, p, &(aux_o->w_list), rid);
            // output_cns_g(cns, e0, s); exit(1);
            // fprintf(stderr, "-1-[M::%s::M]\n", __func__);
        }
        p->x_end = s-1;

        // if((uint32_t)p->x_end + 1 != cal_cigar_xlen_fhc(aux_o)) {
        //     fprintf(stderr, "sb1::%u\n", cal_cigar_xlen_fhc(aux_o));
        // }
    }

    
    // fprintf(stderr, "------[M::%s::M] [%lu, %lu), ec::%u\n", __func__, s, e, cal_cigar_xlen_fhc(aux_o));
    // fprintf(stderr, "------[M::%s::M] [%lu, %lu)\n", __func__, s, e);
    if(p->x_start == -1 || p->x_end == -1) {///hasn't neem set
        p->x_start = s; p->x_end = e-1;
    }
    // nec += push_correct0(p, &(aux_o->w_list), e-s, (uint32_t)-1);    
    nec += push_correct0_fhc(p, &(aux_o->w_list), e-s, (uint32_t)-1, NULL);    
    p->x_end = e-1;
    // if((uint32_t)p->x_end + 1 != cal_cigar_xlen_fhc(aux_o)) {
    //     fprintf(stderr, "ta2::%u\n", cal_cigar_xlen_fhc(aux_o));
    // }
    return nec;
}

void prt_correct0_dbg(overlap_region *in)
{
    int64_t wn = in->w_list.n, k; uint32_t ci, len; asg16_v ff; uint16_t c, b; 
    char cm[4], cc[4]; 
    cm[0] = 'M'; cm[1] = 'S'; cm[2] = 'I'; cm[3] = 'D';
    cc[0] = 'A'; cc[1] = 'C'; cc[2] = 'G'; cc[3] = 'T';

    for (k = 0; k < wn; k++) {
        fprintf(stderr, "\n[M::%s] w[%ld] [%d, %d) clen::%u\n", __func__, k, in->w_list.a[k].x_start, in->w_list.a[k].x_end + 1, in->w_list.a[k].clen);

        ci = 0; ff.n = ff.m = in->w_list.a[k].clen; ff.a = in->w_list.c.a + in->w_list.a[k].cidx;
        while (ci < ff.n) {
            ci = pop_trace_bp(&ff, ci, &c, &b, &len);
            fprintf(stderr, "|%u%c(%c)", len, cm[c], ((c==1)||(c==2))?(cc[b]):('*')); 
        }

        fprintf(stderr, "|\n");
    }
}

uint64_t wcns_vote(overlap_region* ol, All_reads *rref, char* qstr, uint64_t ql, UC_Read* tu, bit_extz_t *exz, uint64_t *id_a, uint64_t id_n, uint64_t s, uint64_t e, ul_ov_t *c_idx, cc_idx_t *occ, uint64_t occ_tot, double occ_exact, overlap_region *aux_o, asg32_v* b32, cns_gfa *cns, uint32_t rid, uint64_t *nec)
{
    uint64_t k, q[2], rr = 0, os, oe, wl, oc[2], fI; ul_ov_t *p, *gp; overlap_region *z; 
    // uint64_t *ct = occ->idx->a;///occ->idx->a[0, wl<<1)
    for (k = 0; k < id_n; k++) {
        p = &(c_idx[id_a[k]]); z = &(ol[ovlp_id(*p)]); 
        q[0] = z->w_list.a[ovlp_cur_wid(*p)].x_start+ovlp_bd(*p);
        q[1] = z->w_list.a[ovlp_cur_wid(*p)].x_end+1-ovlp_bd(*p);
        if(q[1] <= e) rr = 1;
        os = MAX(q[0], s); oe = MIN(q[1], e);
        if(oe > os) {
            ///prepare for CNS
            gp = &(occ->c_idx[id_a[k]]);
            ovlp_cur_xoff(*gp) = ovlp_cur_xoff(*p); ovlp_cur_yoff(*gp) = ovlp_cur_yoff(*p); ovlp_cur_coff(*gp) = ovlp_cur_coff(*p); ovlp_cur_ylen(*gp) = ovlp_cur_ylen(*p);
            // assert(ovlp_cur_wid(*p) == ovlp_cur_wid(*gp));
            // assert(ovlp_id(*p) == ovlp_id(*gp));
            // fprintf(stderr, "[M::%s] tid::%u\t%.*s\twid::%u\tq::[%u, %u)\tos::%lu\toe::%lu\n", __func__, ol[ovlp_id(*p)].y_id, (int)Get_NAME_LENGTH(R_INF, ol[ovlp_id(*p)].y_id), Get_NAME(R_INF, ol[ovlp_id(*p)].y_id),
            // ovlp_cur_wid(*p), ol[ovlp_id(*p)].w_list.a[ovlp_cur_wid(*p)].x_start, ol[ovlp_id(*p)].w_list.a[ovlp_cur_wid(*p)].x_end+1, os, oe);
            extract_sub_cigar_mm(z, os, oe, p, occ->idx->a + os - s);
        }
    }

    wl = e - s;
    os = occ->mms; oe = occ->mme;

    // fprintf(stderr, "[M::%s] s::%lu\te::%lu\n", __func__, s, e);
    
    for (k = 0; k < wl; k++) {
        //+1 for the reference read
        oc[0] = (occ->idx->a[(k<<1)]>>32) + 1;
        oc[1] = ((uint32_t)occ->idx->a[(k<<1)]) + 1;
        // fprintf(stderr, "-0-p::%lu\toc[0]::%lu\toc[1]::%lu\tgoc[0]::%lu\tgoc[1]::%u\n", s + k, oc[0], oc[1], (ct[(k<<1)+1]>>32) + 1, ((uint32_t)ct[(k<<1)+1]) + 1);
        // if(oc[1] < occ_tot || oc[0] <= 1) {
        //     ct[(k<<1)] = ct[(k<<1)+1] = 0;
        //     continue;
        // }
        // fprintf(stderr, "-1-p::%lu\toc[0]::%lu\toc[1]::%lu\n", s + k, oc[0], oc[1]);

        ///a) pass coverage check; b) no enough coverage
        if(((oc[0] > (oc[1]*occ_exact)) && (oc[0] > (oc[1]-oc[0])) && (oc[1] >= occ_tot) && (oc[0] > 1)) ||
            (oc[1] < occ_tot)) {
            ///note: there might be insertions at q[k-1, k], insead if q[k, k+1]
            fI = 1;
            ///make sure there is no insertion
            //+1 for the reference read
            oc[0] = (occ->idx->a[(k<<1)+1]>>32) + 1;
            oc[1] = ((uint32_t)occ->idx->a[(k<<1)+1]) + 1;
            ///a) pass coverage check; b) no enough coverage
            if((((oc[0] > (oc[1]*occ_exact)) && (oc[0] > (oc[1]-oc[0])) && (oc[1] >= occ_tot) && (oc[0] > 1))) || 
                (oc[1] < occ_tot)) {
                fI = 0;
            }

            if(fI) {
                // fprintf(stderr, "-1-p::%lu\toc[0]::%lu\toc[1]::%u\tgoc[0]::%lu\tgoc[1]::%u\n", s + k, (occ->idx->a[(k<<1)]>>32) + 1, ((uint32_t)occ->idx->a[(k<<1)]) + 1, (occ->idx->a[(k<<1)+1]>>32) + 1, ((uint32_t)occ->idx->a[(k<<1)+1]) + 1);
                if(oe > os && os != ((uint64_t)-1)) {///push previous intervals
                    (*nec) += push_cns_anchor(ol, rref, os, oe, qstr, ql, tu, exz, occ, aux_o, 0, occ_tot, occ_exact, b32, cns, rid);
                }
                os = oe = (uint64_t)-1;
            }

            //+1 for the reference read
            oc[0] = (occ->idx->a[(k<<1)]>>32) + 1;
            oc[1] = ((uint32_t)occ->idx->a[(k<<1)]) + 1;
            if((s+k) == oe) {
                oe++;
            } else {
                if(oe > os && os != ((uint64_t)-1)) {///push previous intervals
                    (*nec) += push_cns_anchor(ol, rref, os, oe, qstr, ql, tu, exz, occ, aux_o, 0, occ_tot, occ_exact, b32, cns, rid);
                }
                os = s+k; oe = s+k+1;
            }
        } else {
            // fprintf(stderr, "-2-p::%lu\toc[0]::%lu\toc[1]::%u\tgoc[0]::%lu\tgoc[1]::%u\n", s + k, (occ->idx->a[(k<<1)]>>32) + 1, ((uint32_t)occ->idx->a[(k<<1)]) + 1, (occ->idx->a[(k<<1)+1]>>32) + 1, ((uint32_t)occ->idx->a[(k<<1)+1]) + 1);
            if(oe > os && os != ((uint64_t)-1)) {///push previous intervals
                (*nec) += push_cns_anchor(ol, rref, os, oe, qstr, ql, tu, exz, occ, aux_o, 0, occ_tot, occ_exact, b32, cns, rid);
            }
            os = oe = (uint64_t)-1;
        }
        occ->idx->a[(k<<1)] = occ->idx->a[(k<<1)+1] = 0;
    }

    occ->mms = occ->mme = (uint64_t)-1;
    if(oe > os && os != ((uint64_t)-1)) {
        occ->mms = os; occ->mme = oe;
    }
    return rr;
}



void print_debug_ovlp_cigar(overlap_region_alloc* ol, asg64_v* idx, kv_ul_ov_t *c_idx)
{
    uint64_t k, ci; uint32_t cl; ul_ov_t *cp; bit_extz_t ez; uint16_t c; char cm[4]; 
    cm[0] = 'M'; cm[1] = 'S'; cm[2] = 'I'; cm[3] = 'D';  
    for (k = 0; k < idx->n; k++) {
        cp = &(c_idx->a[(uint32_t)idx->a[k]]);
        fprintf(stderr, "**********[M::%s] tid::%u\t%.*s\twid::%u\tq::[%u, %u)\terr::%d\toerr::%u**********\n", __func__, ol->list[ovlp_id(*cp)].y_id, (int)Get_NAME_LENGTH(R_INF, ol->list[ovlp_id(*cp)].y_id), Get_NAME(R_INF, ol->list[ovlp_id(*cp)].y_id),
        ovlp_cur_wid(*cp), ol->list[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_start, ol->list[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_end+1, ol->list[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].error, ol->list[ovlp_id(*cp)].non_homopolymer_errors);
        set_bit_extz_t(ez, ol->list[ovlp_id(*cp)], ovlp_cur_wid(*cp)); ci = 0;
        while (ci < ez.cigar.n) {
            ci = pop_trace(&(ez.cigar), ci, &c, &cl);
            fprintf(stderr, "%u%c", cl, cm[c]);
        }
        fprintf(stderr, "\n");
    }    
}

uint64_t wcns_gen(overlap_region_alloc* ol, All_reads *rref, UC_Read* qu, UC_Read* tu, bit_extz_t *exz, kv_ul_ov_t *c_idx, asg64_v* idx, asg64_v* buf, int64_t bd, uint64_t wl, int64_t ql, uint64_t occ_tot, double occ_exact, overlap_region *aux_o, asg32_v* b32, cns_gfa *cns, uint64_t cns_g_wl, uint32_t rid)
{
    int64_t on = ol->length, k, i, zwn, q[2]; cns->cns_g_wl = cns_g_wl;
    uint64_t m, *ra, rn, nec = 0, n_id, l_nid, p[2], li; overlap_region *z; ul_ov_t *cp; 
    bit_extz_t ez; uint64_t ci; uint32_t cl; uint16_t c;

    for (k = idx->n = c_idx->n = 0; k < on; k++) {
        z = &(ol->list[k]); zwn = z->w_list.n; z->without_large_indel = l_nid = 0;
        if((!zwn) || (z->is_match != 1)) continue;
        for (i = 0, li = (uint64_t)-1; i < zwn; i++) {
            if(is_ualn_win(z->w_list.a[i])) {
                n_id = z->w_list.a[i].x_end + 1 - z->w_list.a[i].x_start; 
                if(n_id >= 6) l_nid = 1;
                continue;
            }
            q[0] = z->w_list.a[i].x_start; q[1] = z->w_list.a[i].x_end;
            q[0] += bd; q[1] -= bd;
            if(q[1] >= q[0]) {
                m = ((uint64_t)q[0]); m <<= 32; 
                m += c_idx->n; kv_push(uint64_t, *idx, m);

                kv_pushp(ul_ov_t, *c_idx, &cp);
                ovlp_id(*cp) = k; ///ovlp id
                // ovlp_min_wid(*cp) = i; ///beg id of windows
                // ovlp_max_wid(*cp) = i; ///end id of windows
                ovlp_cur_wid(*cp) = i; ///cur id of windows
                ovlp_cur_xoff(*cp) = z->w_list.a[i].x_start; ///cur xpos
                ovlp_cur_yoff(*cp) = z->w_list.a[i].y_start; ///cur xpos
                ovlp_cur_ylen(*cp) = 0;
                ovlp_cur_coff(*cp) = 0; ///cur cigar off in cur window
                ovlp_bd(*cp) = bd;
            }

            if(l_nid == 0) {
                if(i == 0) {
                    p[0] = z->w_list.a[i].x_start; p[1] = z->x_pos_s;
                    n_id = ((p[0] >= p[1])? (p[0] - p[1]): (p[1] - p[0]));
                    if(n_id >= 6) l_nid = 1;
                } 
                
                if(li != (uint64_t)-1) {
                    p[0] = z->w_list.a[i].x_start; p[1] = z->w_list.a[li].x_end + 1;
                    n_id = ((p[0] >= p[1])? (p[0] - p[1]): (p[1] - p[0]));
                    if(n_id >= 6) l_nid = 1;

                    p[0] = z->w_list.a[i].y_start; p[1] = z->w_list.a[li].y_end + 1;
                    n_id = ((p[0] >= p[1])? (p[0] - p[1]): (p[1] - p[0]));
                    if(n_id >= 6) l_nid = 1;
                }

                if(i + 1 == zwn) {
                    p[0] = z->w_list.a[i].x_end; p[1] = z->x_pos_e;
                    n_id = ((p[0] >= p[1])? (p[0] - p[1]): (p[1] - p[0]));
                    if(n_id >= 6) l_nid = 1;
                } 

                if(l_nid == 0) {
                    set_bit_extz_t(ez, (*z), i); ci = 0;
                    while (ci < ez.cigar.n && l_nid == 0) {
                        ci = pop_trace(&(ez.cigar), ci, &c, &cl);
                        if(c >= 2 && cl >= 6) l_nid = 1;
                    }
                }
            }

            li = i;
        }
        z->without_large_indel = (l_nid?0:1);
    }

    int64_t srt_n = idx->n, s, e, t, rr; i = 0;
    radix_sort_ec64(idx->a, idx->a+idx->n);
    for (k = 1, i = 0; k < srt_n; k++) {
        if (k == srt_n || (idx->a[k]>>32) != (idx->a[i]>>32)) {
            if(k - i > 1) {
                for (t = i; t < k; t++) {
                    cp = &(c_idx->a[(uint32_t)idx->a[t]]);
                    // s = ol->list[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_start+ovlp_bd(*cp);
                    // assert(s == (int64_t)(idx->a[i]>>32));
                    m = ol->list[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_end+1-ovlp_bd(*cp);
                    m <<= 32; m += ((uint32_t)idx->a[t]); idx->a[t] = m;
                    // fprintf(stderr, "[M::%s] s::%ld\tsi::%lu\n", __func__, s, (idx->a[i]>>32));
                }
                radix_sort_ec64(idx->a + i, idx->a + k);
            }
            i = k;
        }
    }

    // print_debug_ovlp_cigar(ol, idx, c_idx);

    ///second index
    kv_resize(ul_ov_t, *c_idx, (c_idx->n<<1));
    ul_ov_t *idx_a = NULL, *idx_b = NULL;
    idx_a = c_idx->a; idx_b = c_idx->a;
    memcpy(idx_b, idx_a, c_idx->n * (sizeof((*(idx_a)))));

    kv_resize(uint64_t, *buf, ((wl<<1) + idx->n)); buf->n = ((wl<<1) + idx->n);
    memcpy(buf->a + (wl<<1), idx->a, idx->n * (sizeof((*(idx->a)))));
    memset(buf->a, 0, (wl<<1)*(sizeof((*(idx->a)))));

    cc_idx_t ii_a, ii_b; memset(&ii_a, 0, sizeof(ii_a)); memset(&ii_b, 0, sizeof(ii_b));

    ii_a.c_idx = idx_a; ii_a.idx = idx; ii_a.i = ii_a.i0 = 0; ii_a.srt_n = ii_a.idx->n; ii_a.mms = ii_a.mme = (uint64_t)-1;
    ii_b.c_idx = idx_b; ii_b.idx = buf; ii_b.i = ii_b.i0 = (wl<<1); ii_b.srt_n = ii_b.idx->n; ii_b.mms = ii_b.mme = (uint64_t)-1;

    s = 0; e = wl; e = ((e<=ql)?e:ql); rr = 0; 
    aux_o->w_list.n = aux_o->w_list.c.n = 0; ///for cigar
    for (; s < ql; ) {
        rn = iter_cc_idx_t(ol->list, &ii_a, s, e, rr, 0, &ra);
        // debug_inter0(ol->list, ii_a.c_idx, ii_a.idx->a + ii_a.i0, ii_a.srt_n - ii_a.i0, ra, rn, s, e, 0, 1, "-0-");
        rr = wcns_vote(ol->list, rref, qu->seq, ql, tu, exz, ra, rn, s, e, ii_a.c_idx, &ii_b, occ_tot, occ_exact, aux_o, b32, cns, rid, &nec);
        s += wl; e += wl; e = ((e<=ql)?e:ql);
    }

    if(ii_b.mme > ii_b.mms && ii_b.mms != (uint64_t)-1) {
        nec += push_cns_anchor(ol->list, rref, ii_b.mms, ii_b.mme, qu->seq, ql, tu, exz, &ii_b, aux_o, 0, occ_tot, occ_exact, b32, cns, rid);
    }

    nec += push_cns_anchor(ol->list, rref, ql, ql, qu->seq, ql, tu, exz, &ii_b, aux_o, 1, occ_tot, occ_exact, b32, cns, rid);



    ///for debug
    // zwn = aux_o->w_list.n;
    // fprintf(stderr, "\n******[M::%s]****** wn::%ld, ql::%ld\n", __func__, zwn, ql);
    // prt_correct0_dbg(aux_o);
    // for (i = 0; i < zwn; i++) {
    //     fprintf(stderr, "[M::%s] (%ld)[%d, %d)\n", __func__, i, aux_o->w_list.a[i].x_start, aux_o->w_list.a[i].x_end+1);
    // }
    return nec;
}

void push_nec_re(overlap_region *in, asg16_v *ou)
{
    assert(in->w_list.n == 1);
    uint32_t n1 = in->w_list.a[0].clen;
    if(n1 > ou->m) {
        REALLOC(ou->a, n1); ou->m = n1;
    }
    ou->n = n1;
    memcpy(ou->a, in->w_list.c.a + in->w_list.a[0].cidx, n1*sizeof((*(in->w_list.c.a))));
}

uint32_t extract_max_exact_sub(asg16_v *in, int64_t xs0, int64_t xe0, int64_t ys0, int64_t ye0, int64_t *exk, int64_t *eyk, int64_t *eck, uint64_t *rxs, uint64_t *rxe, uint64_t *rys, uint64_t *rye)
{
    int64_t xk = *exk, yk = *eyk, ck = *eck, cn = in->n, ol, wx[2], wy[2], os, oe; uint16_t op, bq, bt; uint32_t cl, ovlp;
    *rxs = *rxe = *rys = *rye = 0;

    if((ck < 0) || (ck > cn)) {//(*ck) == cn is allowed
        ck = 0; xk = 0; yk = 0;
    }

    while (ck > 0 && xk >= xs0) {///x -> t; y -> p; first insertion and then match/mismatch
        --ck;
        op = in->a[ck]>>14;
        // ol = (((op == 1) || (op == 2))?(in->a[ck]&(0xfff)):(in->a[ck]&(0x3fff)));
        if((op == 2) || (op == 3)) {
            ol = in->a[ck]&(0xfff);
        } else if(op == 1) {
            ol = in->a[ck]&(0x3ff);
        } else {
            ol = in->a[ck]&(0x3fff);
        }
        if(op != 2) xk -= ol;
        if(op != 3) yk -= ol;
    }

    while (ck < cn && xk < xe0) {
        wx[0] = xk; wy[0] = yk; 
        // ck = pop_trace_bp(in, ck, &op, &b, &cl);
        ck = pop_trace_bp_f(in, ck, &op, &bq, &bt, &cl);
        if(op != 2) xk += cl;
        if(op != 3) yk += cl;
        wx[1] = xk; wy[1] = yk;  
        if(op == 0) {
            os = MAX(xs0, wx[0]); oe = MIN(xe0, wx[1]);
            ovlp = ((oe>os)? (oe-os):0);
            if((ovlp > 0) && (ovlp > (*rxe) - (*rxs))) {
                (*rxs) = wy[0] + os - wx[0]; 
                (*rxe) = wy[0] + oe - wx[0]; 

                (*rys) = ys0 + os - xs0; 
                (*rye) = ys0 + oe - xs0; 
            }
        }
    }

    *exk = xk; *eyk = yk; *eck = ck;
    if((*rxe) > (*rxs)) return 1;
    return 0; 
}

void debug_extract_max_exact_sub(uint32_t qid, UC_Read* qu, UC_Read* tu)
{
    uint32_t ci = 0, len, xk, yk, wx[2], wy[2], k; uint16_t c, b; 

    recover_UC_Read(tu, &R_INF, qid);

    ci = 0; yk = 0;
    while (ci < scc.a[qid].n) {
        ci = pop_trace_bp(&scc.a[qid], ci, &c, &b, &len);
        if(c != 3) yk += len;
    }
    resize_UC_Read(qu, yk);

    ci = 0; xk = yk = 0; 
    while (ci < scc.a[qid].n) {
        wx[0] = xk; wy[0] = yk;
        ci = pop_trace_bp(&scc.a[qid], ci, &c, &b, &len);
        if(c != 2) xk += len;
        if(c != 3) yk += len;
        wx[1] = xk; wy[1] = yk;
        if(c == 0) {
            // memcpy(p->a + wy[0], p->z.seq + wx[0], (wx[1]-wx[0])*sizeof((*(p->a))));
            for (; wx[0] < wx[1]; wx[0]++, wy[0]++) {
                qu->seq[wy[0]] = tu->seq[wx[0]];
            }
        } else if(c == 1 || c == 2) {
            for (k = wy[0]; k < wy[1]; k++) {
                qu->seq[k] = s_H[b]; 
            }
        }
        // if(i == 700) fprintf(stderr, "|%u%c(%c)(x::%u)(y::%u)", len, cm[c], ((c==1)||(c==2))?(cc[b]):('*'), wx[1], wy[1]); // s_H
    }
}

uint64_t extract_max_exact(overlap_region *z, asg16_v *ec, /**UC_Read *qu, UC_Read *tu,**/ uint32_t *rxs, uint32_t *rxe, uint32_t *rys, uint32_t *rye)
{
    // if(z->x_id == 75 && z->y_id == 59) {
    //     fprintf(stderr, "[M::%s]\tz->y_id::%u\n", __func__, z->y_id);
    //     if(z->y_pos_strand) {
    //         recover_UC_Read_RC(tu, &R_INF, z->y_id);
    //     } else {
    //         recover_UC_Read(tu, &R_INF, z->y_id); ///b->z.length
    //     }
    // }

    *rxs = *rxe = *rys = *rye = (uint32_t)-1;
    uint64_t k, rx[2], ry[2], xk, yk, wx[2], wy[2], mx[2], my[2]; 
    uint32_t cl, ck; uint16_t c; asg16_v ct; int64_t exk, eyk, eck;
    exk = eyk = eck = 0; mx[0] = mx[1] = my[0] = my[1] = 0;
    for (k = 0; k < z->w_list.n; k++) {
        ct.a = z->w_list.c.a + z->w_list.a[k].cidx;
        ct.n = ct.m = z->w_list.a[k].clen; ck = 0; 
        xk = z->w_list.a[k].x_start; yk = z->w_list.a[k].y_start;
        while (ck < ct.n) {
            wx[0] = xk; wy[0] = yk;
            ck = pop_trace(&ct, ck, &c, &cl);
            if(c != 2) xk += cl;
            if(c != 3) yk += cl;
            wx[1] = xk; wy[1] = yk;
            // if(c == 0) {
            //     if(memcmp(qref + wx[0], tu->seq + wy[0], wx[1] - wx[0])) {
            //         fprintf(stderr, "-0-[M::%s]\teq::[%lu,\t%lu)\tet::[%lu,\t%lu)\t%c\n", __func__, wx[0], wx[1], wy[0], wy[1], "+-"[z->y_pos_strand]);
            //         // exit(1);
            //     }
            // }
            if(wx[1] <= wx[0]) continue;
            if((wx[1] - wx[0]) <= (mx[1] - mx[0])) continue;
            if((c == 0) && (extract_max_exact_sub(ec, wx[0], wx[1], wy[0], wy[1], &exk, &eyk, &eck, &(rx[0]), &(rx[1]), &(ry[0]), &(ry[1])))) {
                // if(z->x_id == 75 && z->y_id == 59) {
                //     if(memcmp(qu->seq + rx[0], tu->seq + ry[0], rx[1] - rx[0])) {
                //         fprintf(stderr, "-1-[M::%s]\tzq::[%lu,\t%lu)\tzt::[%lu,\t%lu)\teq::[%lu,\t%lu)\tet::[%lu,\t%lu)\n", 
                //         __func__, wx[0], wx[1], wy[0], wy[1], rx[0], rx[1], ry[0], ry[1]);
                //         exit(1);
                //     } 
                //     else {
                //         fprintf(stderr, "-2-[M::%s]\tzq::[%lu,\t%lu)\tzt::[%lu,\t%lu)\teq::[%lu,\t%lu)\tet::[%lu,\t%lu)\n", 
                //         __func__, wx[0], wx[1], wy[0], wy[1], rx[0], rx[1], ry[0], ry[1]);

                //         fprintf(stderr, "[M::%s] qstr::%.*s\n", __func__, ((int)(rx[1] - rx[0])), qu->seq + rx[0]);
                //         fprintf(stderr, "[M::%s] tstr::%.*s\n", __func__, ((int)(ry[1] - ry[0])), tu->seq + ry[0]);
                //     }
                // }
                if((rx[1] - rx[0]) > (mx[1] - mx[0])) {
                    mx[0] = rx[0]; mx[1] = rx[1];
                    my[0] = ry[0]; my[1] = ry[1];
                }
            }
        }
    }

    if(mx[1] > mx[0]) {
        *rxs = mx[0]; *rxe = mx[1]; 
        *rys = my[0]; *rye = my[1]; 
        return 1;
    }
    
    return 0;
}

void push_ne_ovlp(ma_hit_t_alloc* paf, overlap_region_alloc* ov, uint32_t flag, All_reads* R_INF, asg16_v *ec/**, uint64_t qid, UC_Read *qu, UC_Read *tu**/)
{
    // if(qu && tu) {
    //     debug_extract_max_exact_sub(qid, qu, tu);
    // }
    uint64_t k, n; ma_hit_t *z; uint32_t rxs, rxe, rys, rye;
    for (k = n = 0; k < ov->length; k++) {
        if(ov->list[k].is_match == flag) n++;
    }

    if(n > paf->size) {
        paf->size = n;
        REALLOC(paf->buffer, paf->size);
    }

    for (k = paf->length = 0; k < ov->length; k++) {
        if(ov->list[k].is_match == flag) {
            // fprintf(stderr, "@%s\tSN:%.*s(id::%u)\terr::%u\n", flag==1?"SQ":"RQ", (int32_t)Get_NAME_LENGTH((*R_INF), ov->list[k].y_id), Get_NAME((*R_INF), ov->list[k].y_id), ov->list[k].y_id, ov->list[k].non_homopolymer_errors);

            z = &(paf->buffer[paf->length++]);

            z->qns = ov->list[k].x_id;
            z->qns = z->qns << 32;
            z->tn = ov->list[k].y_id;

            z->qns = z->qns | (uint64_t)(ov->list[k].x_pos_s);
            z->qe = ov->list[k].x_pos_e + 1;
            z->ts = ov->list[k].y_pos_s;
            z->te = ov->list[k].y_pos_e + 1;

            ///for overlap_list, the x_strand of all overlaps are 0, so the tmp.rev is the same as the y_strand
            z->rev = ov->list[k].y_pos_strand;

            z->bl = Get_READ_LENGTH((*R_INF), ov->list[k].y_id);
            z->ml = ov->list[k].strong;
            z->no_l_indel = ov->list[k].without_large_indel;

            if(ec) {
                extract_max_exact(&ov->list[k], ec, /**qu, tu,**/ &rxs, &rxe, &rys, &rye);
                z->el = 0;
                // fprintf(stderr, "[M::%s]\tq::[%u,\t%u)\tt::[%u,\t%u)\teq::[%u,\t%u)\tet::[%u,\t%u)\n", __func__, ov->list[k].x_pos_s, ov->list[k].x_pos_e + 1, ov->list[k].y_pos_s, ov->list[k].y_pos_e + 1, rxs, rxe, rys, rye);
                if(rxe > rxs) {
                    z->qns = ov->list[k].x_id;
                    z->qns = z->qns << 32;
                    z->qns = z->qns | (uint64_t)(rxs);
                    z->qe = rxe;
                    z->ts = rys;
                    z->te = rye;

                    z->el = 1;
                }
            }
        }
    }
}

void push_ff_ovlp(ma_hit_t_alloc* paf, overlap_region_alloc* ov, uint32_t flag, All_reads* R_INF, uint64_t *cnt)
{
    // if(qu && tu) {
    //     debug_extract_max_exact_sub(qid, qu, tu);
    // }
    uint64_t k, n; ma_hit_t *z; 
    for (k = n = 0; k < ov->length; k++) {
        if(ov->list[k].is_match == flag) n++;
    }

    if(n > paf->size) {
        paf->size = n;
        REALLOC(paf->buffer, paf->size);
    }

    for (k = paf->length = 0; k < ov->length; k++) {
        if(ov->list[k].is_match == flag) {
            z = &(paf->buffer[paf->length++]);

            z->qns = ov->list[k].x_id;
            z->qns = z->qns << 32;
            z->tn = ov->list[k].y_id;

            z->qns = z->qns | (uint64_t)(ov->list[k].x_pos_s);
            z->qe = ov->list[k].x_pos_e + 1;
            z->ts = ov->list[k].y_pos_s;
            z->te = ov->list[k].y_pos_e + 1;

            ///for overlap_list, the x_strand of all overlaps are 0, so the tmp.rev is the same as the y_strand
            z->rev = ov->list[k].y_pos_strand;

            z->bl = Get_READ_LENGTH((*R_INF), ov->list[k].y_id);
            z->ml = ov->list[k].strong;
            z->no_l_indel = ov->list[k].without_large_indel;
            z->el = ov->list[k].shared_seed;

            if(z->rev) {
                z->ts = z->bl - ov->list[k].y_pos_e - 1;
                z->te = z->bl - ov->list[k].y_pos_s;
            }

            if(flag == 1) {
                if(z->ml == 1) cnt[2]++;
                if(z->ml == 0) cnt[3]++;
                if(z->el == 1) cnt[4]++;
                if(z->no_l_indel) cnt[5]++;
            }
            z->del = 0;
        }
    }

    if(flag == 1) cnt[0] += paf->length;
    if(flag == 2) cnt[1] += paf->length;
}

void debug_mm_exact_cigar(overlap_region_alloc* ol, uint32_t qid, UC_Read *qu, UC_Read *tu)
{
    int64_t on = ol->length, k, i, zwn, xk, yk; uint32_t cl, ck; ///bit_extz_t ez; 
    overlap_region *z; asg16_v ct; uint64_t wx[2], wy[2]; uint16_t c;
    recover_UC_Read(qu, &R_INF, qid);

    for (i = 0; i < on; i++) {
        z = &(ol->list[i]); zwn = z->w_list.n; 
        if(z->y_pos_strand) {
            recover_UC_Read_RC(tu, &R_INF, z->y_id);
        } else {
            recover_UC_Read(tu, &R_INF, z->y_id); ///b->z.length
        }

        // fprintf(stderr, "[M::%s] x_id::%u\ty_id::%u\tx::[%u, %u)\ty::[%u, %u)\tzwn::%ld\n", 
        //         __func__, z->x_id, z->y_id, z->x_pos_s, z->x_pos_e + 1, z->y_pos_s, z->y_pos_e + 1, zwn);
        // fprintf(stderr, "qstr(%lld)::%.*s\n", qu->length, (int32_t)(qu->length), qu->seq);
        // fprintf(stderr, "tstr(%lld)::%.*s\n", tu->length, (int32_t)(tu->length), tu->seq);

        for (k = 0; k < zwn; k++) {
            // set_bit_extz_t(ez, (*z), k);
            // if(!cigar_check(tu->seq, qu->seq, &ez)) {
            //     fprintf(stderr, "\n[M::%s] x_id::%u, y_id::%u, x::[%u, %u), y::[%u, %u)\n", __func__, z->x_id, z->y_id, z->x_pos_s, z->x_pos_e + 1, z->y_pos_s, z->y_pos_e + 1);
            //     exit(1);
            // }
            // continue;
            ct.a = z->w_list.c.a + z->w_list.a[k].cidx;
            ct.n = ct.m = z->w_list.a[k].clen; ck = 0; 
            xk = z->w_list.a[k].x_start; yk = z->w_list.a[k].y_start;
            while (ck < ct.n) {
                wx[0] = xk; wy[0] = yk;
                ck = pop_trace(&ct, ck, &c, &cl);
                if(c != 2) xk += cl;
                if(c != 3) yk += cl;
                wx[1] = xk; wy[1] = yk;
                if(c == 0) {
                    if(memcmp(qu->seq + wx[0], tu->seq + wy[0], wx[1] - wx[0])) {
                        fprintf(stderr, "\n-0-[M::%s]\teq::[%lu,\t%lu)\tet::[%lu,\t%lu)\t%c\n", __func__, wx[0], wx[1], wy[0], wy[1], "+-"[z->y_pos_strand]);
                        fprintf(stderr, "qstr(%u)::%.*s\n", z->x_id, (int32_t)(wx[1] - wx[0]), qu->seq + wx[0]);
                        fprintf(stderr, "tstr(%u)::%.*s\n", z->y_id, (int32_t)(wy[1] - wy[0]), tu->seq + wy[0]);
                        // exit(1);
                    } 
                    // else {
                    //     fprintf(stderr, "\n-1-[M::%s]\teq::[%lu,\t%lu)\tet::[%lu,\t%lu)\t%c\n", __func__, wx[0], wx[1], wy[0], wy[1], "+-"[z->y_pos_strand]);
                    //     fprintf(stderr, "qstr(%u)::%.*s\n", z->x_id, (int32_t)(wx[1] - wx[0]), qu->seq + wx[0]);
                    //     fprintf(stderr, "tstr(%u)::%.*s\n", z->y_id, (int32_t)(wy[1] - wy[0]), tu->seq + wy[0]);
                    // }
                }
            }
        }
        
    }
}

void check_well_cal(asg16_v *sc, asg64_v *idx, uint8_t *f_ec, uint8_t *abnormal, int64_t len, int64_t min_dp, ma_hit_t_alloc *in)
{
    uint64_t k, s, e; int64_t dp, old_dp, st = 0, ed; ma_hit_t *z; 

    (*f_ec) = 1; (*abnormal) = 0; idx->n = 0;
    for (k = 0; k < in->length; k++) {
        z = &(in->buffer[k]);
        s = ((uint32_t)(z->qns)); e = z->qe;
        kv_push(uint64_t, (*idx), (s<<1));
        kv_push(uint64_t, (*idx), (e<<1)|1);
    }

    radix_sort_ec64(idx->a, idx->a + idx->n);
    for (k = 0, dp = 0, st = ed = 0; k < idx->n; ++k) {
        old_dp = dp;
        ///if a[j] is qe
        if (idx->a[k]&1) --dp;
        else ++dp;

        ed = idx->a[k]>>1;
        if(ed > st) {
            if(old_dp < min_dp) (*f_ec) = 0;
            if(old_dp == 0) {
                if(st > 0 && ed < len) {
                    (*abnormal) = 1;
                }else if((*abnormal)==0){
                    (*abnormal) = 2;
                }
            }
        }

        st = ed;
    }


    ed = len; old_dp = dp;
    if(ed > st) {
        if(old_dp < min_dp) (*f_ec) = 0;
        if(old_dp == 0) {
            if(st > 0 && ed < len) {
                (*abnormal) = 1;
            }else if((*abnormal) == 0){
                (*abnormal) = 2;
            }
        }
    }

    if((*f_ec)) {
        for (k = 0; (k < sc->n) && ((sc->a[k]>>14) == 0); k++);
        if(k < sc->n) (*f_ec) = 0;
    }    
}

inline uint64_t exact_ec_check(char *qstr, uint64_t ql, char *tstr, uint64_t tl, int64_t qs, int64_t qe, int64_t ts, int64_t te)
{
    if(qe - qs != te - ts) return 0;
    if(memcmp(qstr + qs, tstr + ts, qe - qs) == 0) return 1;
    return 0;
}

void gen_hc_r_alin_ea(overlap_region_alloc* ol, Candidates_list *cl, All_reads *rref, UC_Read* qu, UC_Read* tu, bit_extz_t *exz, overlap_region *aux_o, double e_rate, int64_t wl, int64_t rid, int64_t khit, int64_t move_gap, asg16_v *buf, asg64_v *srt, ma_hit_t_alloc *in, uint8_t chem_drop)
{
    if(ol->length <= 0) return;

    // uint64_t k, l, i, s, m, mm_k, *ei, en, *oi, on, tid, trev, nec; int64_t sc, mm_sc, plus, minus; overlap_region *z, t; ma_hit_t *p;
    uint64_t k, i, m, *ei, en, *oi, on, tid, trev, nec; overlap_region *z; ma_hit_t *p;
    srt->n = 0;
    for (k = 0; k < in->length; k++) {
        if(in->buffer[k].el) {
            m = in->buffer[k].tn; m <<= 1; m |= in->buffer[k].rev; 
            m <<= 32; m |= k; kv_push(uint64_t, (*srt), m);
        }
    }

    if(!(srt->n)) {
        gen_hc_r_alin(ol, cl, rref, qu, tu, exz, aux_o, e_rate, wl, rid, khit, move_gap, buf, chem_drop);
    } else {
        ///debug for memory
        // snprintf(NULL, 0, "dwn::%u\tdcn::%u", (uint32_t)aux_o->w_list.n, (uint32_t)aux_o->w_list.c.n);

        kv_resize(uint64_t, *srt, (srt->n + ol->length));
        ei = srt->a; en = srt->n; oi = srt->a + srt->n; on = ol->length;
        for (k = 0; k < on; k++) {
            z = &(ol->list[k]); z->is_match = z->strong = z->without_large_indel = 0;
            oi[k] = z->y_id; oi[k] <<= 1; oi[k] |= z->y_pos_strand;
            oi[k] <<= 32; oi[k] |= k;
        }

        radix_sort_ec64(ei, ei + en); radix_sort_ec64(oi, oi + on);
        for (k = i = nec = 0; k < on; k++) {
            z = &(ol->list[(uint32_t)oi[k]]); tid = z->y_id; trev = z->y_pos_strand;
            for (; (i < en) && ((ei[i]>>32) < ((tid<<1)|trev)); i++);
            if((i < en) && ((ei[i]>>32) == ((tid<<1)|trev))) {
                p = &(in->buffer[(uint32_t)ei[i]]);
                if((z->x_pos_s == ((uint32_t)p->qns)) && (z->x_pos_e + 1 == p->qe) && 
                            (z->y_pos_s == p->ts) && (z->y_pos_e + 1 == p->te)) {
                    resize_UC_Read(tu, p->te - p->ts); recover_UC_Read_sub_region(tu->seq, p->ts, p->te - p->ts, trev, rref, tid);
                    if(exact_ec_check(qu->seq, qu->length, tu->seq, p->te - p->ts, ((uint32_t)p->qns), p->qe, 0, p->te - p->ts)) {
                        z->is_match = 1; z->shared_seed = z->non_homopolymer_errors;///for index
                        z->non_homopolymer_errors = 0; z->strong = z->without_large_indel = 0;
                        set_exact_exz(exz, z->x_pos_s, z->x_pos_e + 1, z->y_pos_s, z->y_pos_e + 1); push_alnw(z, exz);
                        nec++;
                    }
                }
            }
        }
        ///debug for memory
        // snprintf(NULL, 0, "dwn::%u\tdcn::%u", (uint32_t)aux_o->w_list.n, (uint32_t)aux_o->w_list.c.n);

        if(on > nec) gen_hc_r_alin_nec(ol, cl, rref, qu, tu, exz, aux_o, e_rate, wl, rid, khit, move_gap, buf, chem_drop);

        ///debug for memory
        // snprintf(NULL, 0, "dwn::%u\tdcn::%u", (uint32_t)aux_o->w_list.n, (uint32_t)aux_o->w_list.c.n);
    }

    /**
    if(ol->length > 1) {///for duplicated chains
        overlap_region_sort_y_id(ol->list, ol->length);
        for (k = 1, l = m = 0; k <= ol->length; k++) {
            if(k == ol->length || ol->list[k].y_id != ol->list[l].y_id) {
                // fprintf(stderr, "\n[M::%s::tid->%u] n->%lu\n", __func__, ol->list[l].y_id, k - l);
                mm_k = l;
                if(k - l > 1) {
                    for (s = l, mm_sc = INT32_MIN, mm_k = ((uint64_t)-1); s < k; s++) {
                        z = &(ol->list[s]);
                        plus = z->x_pos_e + 1 - z->x_pos_s; minus = (z->non_homopolymer_errors) * 12;
                        sc = plus - minus;
                        if((sc > mm_sc) || ((sc == mm_sc) && ((ol->list[mm_k].x_pos_e+1-ol->list[mm_k].x_pos_s) < (z->x_pos_e+1-z->x_pos_s)))) {
                            mm_sc = sc; mm_k = s;
                        }
                        // fprintf(stderr, "[M::%s::%c] q::[%u, %u), t::[%u, %u), sc::%ld, err::%u, s::%lu\n", __func__, "+-"[z->y_pos_strand], z->x_pos_s, z->x_pos_e + 1, z->y_pos_s, z->y_pos_e + 1, sc, z->non_homopolymer_errors, s);
                    }
                }
                // fprintf(stderr, "[M::%s::tid->%u] mm_k::%lu\n", __func__, ol->list[l].y_id, mm_k);
                if(mm_k != ((uint64_t)-1)) {
                    if(mm_k != m) {
                        t = ol->list[mm_k];
                        ol->list[mm_k] = ol->list[m];
                        ol->list[m] = t;
                    }
                    m++;
                }
                l = k;
            }
        }
        ol->length = m;
    }
    **/
}

void prt_ovlp_sam_0(char *cm, FILE *fp, char *ref_id, int32_t ref_id_n, char *qry_id, int32_t qry_id_n, char *qry_seq, uint64_t qry_seq_n, uint64_t rs, uint64_t re, uint64_t qs, uint64_t qe, uint64_t flag, uint64_t err, bit_extz_t *ez)
{
    uint64_t ci = 0; uint16_t c; uint32_t cl;
    fprintf(fp, "%.*s\t%lu\t%.*s\t%lu\t60\t", qry_id_n, qry_id, flag, ref_id_n, ref_id, rs + 1);

    if(qs) fprintf(fp, "%luS", qs);
    while (ci < ez->cigar.n) {
        ci = pop_trace(&(ez->cigar), ci, &c, &cl);
        fprintf(fp, "%u%c", cl, cm[c]);
    }
    if(qry_seq_n > qe) fprintf(fp, "%luS", qry_seq_n - qe);
    fprintf(fp, "\t*\t0\t0\t%.*s\t", (int32_t)qry_seq_n, qry_seq);
    for (ci = 0; ci < qry_seq_n; ci++) fprintf(fp, "~");
    fprintf(fp, "\tNM:i:%lu\n", err);
}


void prt_ovlp_sam(overlap_region_alloc* ol, UC_Read* tu, char *ref_seq, int32_t ref_seq_n)
{
    int64_t on = ol->length, k, i, zwn; overlap_region *z; bit_extz_t ez;
    char *qry = NULL, *ref = Get_NAME(R_INF, ol->list[0].x_id); 
    uint64_t qry_n = 0, ref_n = Get_NAME_LENGTH(R_INF, ol->list[0].x_id), qid, rev;
    char cm[4]; cm[0] = 'M'; cm[1] = 'S'; cm[2] = 'I'; cm[3] = 'D';  
    FILE *fp = fopen("aln.sam", "w");
    fprintf(fp, "@HD\tVN:1.6\tSO:unknown\n");
    fprintf(fp, "@SQ\tSN:%.*s\tLN:%lu\n", (int32_t)ref_n, ref, Get_READ_LENGTH(R_INF, ol->list[0].x_id));
    for (k = 0; k < on; k++) {
        z = &(ol->list[k]); zwn = z->w_list.n; 
        if(!zwn) continue; 
        qid = ol->list[k].y_id; 
        if(z->y_pos_strand) {
            recover_UC_Read_RC(tu, &R_INF, qid); rev = 16;
        } else {
            recover_UC_Read(tu, &R_INF, qid); rev = 0;
        }
        for (i = 0; i < zwn; i++) {
            if(is_ualn_win(z->w_list.a[i])) continue;
            qry_n = Get_NAME_LENGTH(R_INF, qid); 
            qry = Get_NAME(R_INF, qid);
            set_bit_extz_t(ez, (*z), i);
            
            prt_ovlp_sam_0(cm, fp, ref, ref_n, qry, qry_n, tu->seq, tu->length, z->w_list.a[i].x_start, z->w_list.a[i].x_end + 1, z->w_list.a[i].y_start, z->w_list.a[i].y_end + 1, rev, z->w_list.a[i].error, &ez);
        }
    }
    fclose(fp);

    fp = fopen("ref.fa", "w");
    fprintf(fp, ">%.*s\n", (int32_t)ref_n, ref);
    fprintf(fp, "%.*s\n", ref_seq_n, ref_seq);
    fclose(fp);
}

void stderr_phase_ovlp(overlap_region_alloc* ol)
{
    int64_t on = ol->length, k; overlap_region *z;
    if(!on) return;
    uint64_t qry_n = 0, rid, ref_n, qid;
    rid = ol->list[0].x_id; ref_n = Get_NAME_LENGTH(R_INF, rid);
    
    for (k = 0; k < on; k++) {
        z = &(ol->list[k]); qid = ol->list[k].y_id; 
        qry_n = Get_NAME_LENGTH(R_INF, qid); 

        fprintf(stderr, "%.*s(qid::%lu)\tql::%lu\tq::[%u,\t%u)\t%c\t%.*s(tid::%lu)\ttl::%lu\tt::[%u,\t%u)\ttrans::%u\n", 
        (int32_t)Get_NAME_LENGTH(R_INF, rid), Get_NAME(R_INF, rid), rid, ref_n, z->x_pos_s, z->x_pos_e + 1, "+-"[z->y_pos_strand], 
        (int32_t)Get_NAME_LENGTH(R_INF, qid), Get_NAME(R_INF, qid), qid, qry_n, z->y_pos_s, z->y_pos_e + 1, ((z->is_match==1)?(0):(1)));
    }
}

void dedup_chains(overlap_region_alloc* ol)
{
    uint64_t k, l, s, m, mm_k, mm_m, sf; int64_t sc, mm_sc, plus, minus; overlap_region *z, t; 
    if(ol->length > 1) {///for duplicated chains
        overlap_region_sort_y_id(ol->list, ol->length);
        for (k = 1, l = m = 0; k <= ol->length; k++) {
            if(k == ol->length || ol->list[k].y_id != ol->list[l].y_id) {
                // fprintf(stderr, "\n[M::%s::tid->%u] n->%lu\n", __func__, ol->list[l].y_id, k - l);
                mm_k = l;
                if(k - l > 1) {
                    for (s = l, mm_sc = INT32_MIN, mm_k = ((uint64_t)-1), mm_m = 3; s < k; s++) {
                        z = &(ol->list[s]);
                        plus = z->x_pos_e + 1 - z->x_pos_s; minus = (z->non_homopolymer_errors) * 12;
                        sc = plus - minus; 
                        
                        sf = 0;
                        if(z->is_match < mm_m) {
                            sf = 1;
                        } else if(z->is_match == mm_m) {
                            if(sc > mm_sc) {
                                sf = 1;
                            } else if((sc == mm_sc) && ((z->x_pos_e+1-z->x_pos_s) > (ol->list[mm_k].x_pos_e+1-ol->list[mm_k].x_pos_s))) {
                                sf = 1;
                            }
                        }
                        if(sf) {
                            mm_sc = sc; mm_k = s; mm_m = z->is_match;
                        }
                        
                        // fprintf(stderr, "[M::%s::%c] q::[%u, %u), t::[%u, %u), sc::%ld, err::%u, mm::%u, s::%lu\n", __func__, "+-"[z->y_pos_strand], z->x_pos_s, z->x_pos_e + 1, z->y_pos_s, z->y_pos_e + 1, sc, z->non_homopolymer_errors, z->is_match, s);
                    }
                }
                // fprintf(stderr, "[M::%s::tid->%u] mm_k::%lu\n", __func__, ol->list[l].y_id, mm_k);
                if(mm_k != ((uint64_t)-1)) {
                    if(mm_k != m) {
                        t = ol->list[mm_k];
                        ol->list[mm_k] = ol->list[m];
                        ol->list[m] = t;
                    }
                    m++;
                }
                l = k;
            }
        }
        ol->length = m;
    }
}

void debug_retrive_bqual(asg8_v *vq, asg8_v *vt, uint64_t id, uint64_t rn)
{
    uint64_t k, n, z[2], s, e, rev;
    retrive_bqual(vq, NULL, id, -1, -1, 0, sc_bn); n = vq->n;
    retrive_bqual(vt, NULL, id, -1, -1, 1, sc_bn);
    assert(vq->n == vt->n);
    for (k = 0; k < vq->n && vq->a[k] == vt->a[vt->n - k - 1]; k++);
    // if((k == vq->n)) {
    //     fprintf(stderr, "[M::%s] id::%lu, k::%lu, n::%lu\n", __func__, id, k, ((uint64_t)vq->n));
    // }
    assert(k == vq->n);

    // if(id == 0) {
    //     s = 21519; e = 22332; rev = 0;
    //     retrive_bqual(vt, NULL, id, s, e, rev, sc_bn);
    //     if(memcmp(vq->a + s, vt->a, e - s)) {
    //         fprintf(stderr, "+[M::%s] id::%lu, t::[%lu, %lu), rev::%lu\n", __func__, id, s, e, rev);
    //         exit(1);
    //     }

    // }
    // return;

    for (k = 0, rev = 0; k < rn; k++) {
        z[0] = rand()%(n + 1); 
        z[1] = rand()%(n + 1); 
        if(z[0] == z[1]) continue;
        s = MIN(z[0], z[1]); e = MAX(z[0], z[1]);
        retrive_bqual(vt, NULL, id, s, e, rev, sc_bn);
        if(memcmp(vq->a + s, vt->a, e - s)) {
            fprintf(stderr, "[M::%s] id::%lu, t::[%lu, %lu), rev::%lu\n", __func__, id, s, e, rev);
            exit(1);
        }
    }
    

    retrive_bqual(vq, NULL, id, -1, -1, 1, sc_bn);
    for (k = 0, rev = 1; k < rn; k++) {
        z[0] = rand()%(n + 1); 
        z[1] = rand()%(n + 1); 
        if(z[0] == z[1]) continue;
        s = MIN(z[0], z[1]); e = MAX(z[0], z[1]);
        retrive_bqual(vt, NULL, id, s, e, rev, sc_bn);
        if(memcmp(vq->a + s, vt->a, e - s)) {
            fprintf(stderr, "[M::%s] id::%lu, t::[%lu, %lu), rev::%lu\n", __func__, id, s, e, rev);
            exit(1);
        }
    }
}

uint32_t is_uncorrected_read(overlap_region_alloc* ov, asg64_v *idx, int64_t len, int64_t min_len)
{
    uint64_t k, s, e; int64_t dp, old_dp, st = 0, ed;
    for (k = idx->n = 0; k < ov->length; k++) {
        s = ov->list[k].x_pos_s; e = ov->list[k].x_pos_e + 1;
        kv_push(uint64_t, (*idx), (s<<1));
        kv_push(uint64_t, (*idx), (e<<1)|1);
    }

    radix_sort_ec64(idx->a, idx->a + idx->n);
    for (k = 0, dp = 0, st = ed = 0; k < idx->n; ++k) {
        old_dp = dp;
        ///if a[j] is qe
        if (idx->a[k]&1) --dp;
        else ++dp;

        ed = idx->a[k]>>1;
        if(ed > st) {
            if(old_dp == 0) {
                if((ed - st) >= min_len) return 1;
            }
        }
        st = ed;
    }


    ed = len; old_dp = dp;
    if(ed > st) {
        if(old_dp == 0) {
            if((ed - st) >= min_len) return 1;
            if((ed - st) >= len) return 1;
        }
    }
    
    return 0;
}

uint32_t is_chemical_r_qual(overlap_region_alloc *ov, asg64_v *idx, int64_t len, int64_t cov, int64_t frank_len, asg8_v *qv, uint64_t rid)
{
    uint64_t k, s, e; int64_t dp, old_dp, st = 0, ed, s0, s1, e0, e1, rr, qk;
    if((frank_len) > (len *0.01)) frank_len = len *0.01;
    for (k = idx->n = 0; k < ov->length; k++) {
        s = ov->list[k].x_pos_s; e = ov->list[k].x_pos_e + 1;
        kv_push(uint64_t, (*idx), (s<<1));
        kv_push(uint64_t, (*idx), (e<<1)|1);
        // fprintf(stderr, "[M::%s]\trid::%lu\ts::%lu\te::%lu\n", __func__, rid, s, e);
    }

    radix_sort_ec64(idx->a, idx->a + idx->n); s0 = s1 = e0 = e1 = rr = -1;
    for (k = 0, dp = 0, st = ed = 0; k < idx->n; ++k) {
        old_dp = dp;
        ///if a[j] is qe
        if (idx->a[k]&1) --dp;
        else ++dp;

        ed = idx->a[k]>>1;
        if(ed > st) {
            if(old_dp <= cov) {
                rr = 1;
            } else {
                if(s0 < 0) {
                    s0 = st; s1 = ed;
                }
                e0 = st; e1 = ed;
            }
        }
        st = ed;
    }


    ed = len; old_dp = dp;
    if(ed > st) {
        if(old_dp <= cov) {
            rr = 1;
        } else {
            if(s0 < 0) {
                s0 = st; s1 = ed;
            }
            e0 = st; e1 = ed;
        }
    }


    if((s0 != e0) && (s1 != e1) && (s0 <= frank_len) && ((len - e1) <= frank_len) && (rr > 0)) {
        for (k = 0, dp = 0, st = ed = 0; k < idx->n; ++k) {
            old_dp = dp;
            ///if a[j] is qe
            if (idx->a[k]&1) --dp;
            else ++dp;

            ed = idx->a[k]>>1;
            if(ed > st) {
                if((old_dp <= cov) && (st >= s0) && (ed <= e1)) {
                    retrive_bqual(qv, NULL, rid, -1, -1, 0, sc_bn); 
                    fprintf(stderr, "[M::%s]\tlf::[%ld,%ld)\trt::[%ld,%ld)\tmd::[%ld,%ld)\tcov::%ld\n", __func__, s0, s1, e0, e1, st, ed, old_dp);
                    for (qk = st; qk < ed; qk++) fprintf(stderr, "%u", qv->a[qk]);
                    fprintf(stderr, "\n");
                    return 1;
                }
            }
            st = ed;
        }


        ed = len; old_dp = dp;
        if(ed > st) {
            if((old_dp <= cov) && (st >= s0) && (ed <= e1)) {
                retrive_bqual(qv, NULL, rid, -1, -1, 0, sc_bn);
                fprintf(stderr, "[M::%s]\tlf::[%ld,%ld)\trt::[%ld,%ld)\tmd::[%ld,%ld)\tcov::%ld\n", __func__, s0, s1, e0, e1, st, ed, old_dp);
                for (qk = st; qk < ed; qk++) fprintf(stderr, "%u", qv->a[qk]);
                fprintf(stderr, "\n");
                return 1;
            }
        }
    }


    // for (k = 0, dp = 0, st = ed = 0; k < idx->n; ++k) {
    //     old_dp = dp;
    //     ///if a[j] is qe
    //     if (idx->a[k]&1) --dp;
    //     else ++dp;

    //     ed = idx->a[k]>>1;
    //     if(ed > st) {
    //         if((old_dp <= cov)) {
    //             retrive_bqual(qv, NULL, rid, -1, -1, 0, sc_bn); 
    //             fprintf(stderr, "[M::%s]\tlf::[%ld,%ld)\trt::[%ld,%ld)\tmd::[%ld,%ld)\tcov::%ld\n", __func__, s0, s1, e0, e1, st, ed, old_dp);
    //             for (qk = st; qk < ed; qk++) fprintf(stderr, "%u", qv->a[qk]);
    //             fprintf(stderr, "\n");
    //             return 1;
    //         }
    //     }
    //     st = ed;
    // }


    // ed = len; old_dp = dp;
    // if(ed > st) {
    //     if((old_dp <= cov)) {
    //         retrive_bqual(qv, NULL, rid, -1, -1, 0, sc_bn);
    //         fprintf(stderr, "[M::%s]\tlf::[%ld,%ld)\trt::[%ld,%ld)\tmd::[%ld,%ld)\tcov::%ld\n", __func__, s0, s1, e0, e1, st, ed, old_dp);
    //         for (qk = st; qk < ed; qk++) fprintf(stderr, "%u", qv->a[qk]);
    //         fprintf(stderr, "\n");
    //         return 1;
    //     }
    // }
    
    return 0;
}


static void worker_hap_ec(void *data, long i, int tid)
{
	ec_ovec_buf_t0 *b = &(((ec_ovec_buf_t*)data)->a[tid]);
    uint32_t high_occ = asm_opt.hom_cov * (2.0 - HA_KMER_GOOD_RATIO);
    uint32_t low_occ = asm_opt.hom_cov * HA_KMER_GOOD_RATIO;
    overlap_region *aux_o = NULL; asg64_v buf0; uint32_t qlen = 0;

    /**
    if((i != 1129685) && (i != 1137865) && (i != 1137917) && (i != 1140647) && (i != 1144740) && (i != 1148936) && (i != 1149134) && (i != 1151224) && (i != 1151386) && (i != 1152960) && (i != 1154846) && (i != 1154881) && (i != 1155112) && 
        (i != 1156823) && (i != 1157099) && (i != 1157393) && (i != 1158300) && (i != 1158368) && (i != 1160411) && (i != 1160659) && (i != 1161458) && (i != 1163595) && (i != 1164084) && (i != 1164230) && (i != 1165050) && (i != 1168249) && 
        (i != 1168304) && (i != 1168514) && (i != 1170447) && (i != 1171377) && (i != 1171387) && (i != 1172376) && (i != 1173566) && (i != 1174275) && (i != 1174434) && (i != 1174511) && (i != 1174860) && (i != 1175306) && (i != 1175314) && 
        (i != 1177101) && (i != 1178196) && (i != 1178470) && (i != 1179327) && (i != 1179626) && (i != 1180357) && (i != 1181347) && (i != 1181422) && (i != 1181725) && (i != 1183135) && (i != 1183734) && (i != 1185569) && (i != 1185604) && 
        (i != 1186192) && (i != 1188441) && (i != 1188487) && (i != 1189865) && (i != 1189943) && (i != 1192819) && (i != 1193133) && (i != 1196908) && (i != 1197590) && (i != 1200549) && (i != 1200757) && (i != 1205500)) return;

    fprintf(stderr, "%ld\t+++\n", i);
    **/
    // if(i < 1100000 || i > 1400000) return;
    // if(i % 100000 == 0) fprintf(stderr, "-a-[M::%s-beg] rid->%ld\n", __func__, i);
    // if (memcmp("c42804f3-0e13-43a0-8a71-b91b40accf9a", Get_NAME((R_INF), i), Get_NAME_LENGTH((R_INF),i)) == 0) {
    // if (memcmp("b2e68ecf-381a-439c-b676-c1e6831d6acf", Get_NAME((R_INF), i), Get_NAME_LENGTH((R_INF),i)) == 0) {
    // if (memcmp("64b2c27d-86b8-451e-9330-6ba62be2ffcc", Get_NAME((R_INF), i), Get_NAME_LENGTH((R_INF),i)) == 0) {
    //     fprintf(stderr, "-a-[M::%s-beg] rid->%ld\n", __func__, i);
    // } else {
    //     return;
    // }

    // if(i != 3028559) return;
    // if(i != 306) return;
    // if(i != 1124) return;
    // if(i != 700) return;
    // if(i != 2243244) return;
    // if(i != 15139) return;

    // debug_retrive_bqual(D, &b->v8t, i, 256); return;

    recover_UC_Read(&b->self_read, &R_INF, i); qlen = b->self_read.length; 

    h_ec_lchain(b->ab, i, b->self_read.seq, b->self_read.length, asm_opt.mz_win, asm_opt.k_mer_length, &R_INF, &b->olist, &b->clist, ((asm_opt.is_ont)?(0.05):(0.02)), asm_opt.max_n_chain, 1, NULL, NULL, &(b->sp), &high_occ, &low_occ, 1, 1, 3, 0.7, 2, 32, COV_W);///ONT high error

    // b->num_read_base += b->olist.length;
    b->cnt[0] += b->self_read.length;

    aux_o = fetch_aux_ovlp(&b->olist);///must be here

    // stderr_phase_ovlp(&b->olist);

    ///debug for memory
    // snprintf(NULL, 0, "dwn::%u\tdcn::%u", (uint32_t)aux_o->w_list.n, (uint32_t)aux_o->w_list.c.n);

    gen_hc_r_alin_ea(&b->olist, &b->clist, &R_INF, &b->self_read, &b->ovlp_read, &b->exz, aux_o, asm_opt.max_ov_diff_ec, (asm_opt.is_ont)?(WINDOW_OHC):(WINDOW_HC), i, E_KHIT/**asm_opt.k_mer_length**/, 1, &b->v16, &b->v64, &(R_INF.paf[i]), asm_opt.is_ont);

    // prt_ovlp_sam(&b->olist, &b->ovlp_read, b->self_read.seq, b->self_read.length);


    // fprintf(stderr, "\n[M::%s] rid::%ld\t%.*s\tlen::%lld\tocc::%lu\n", __func__, i, (int)Get_NAME_LENGTH(R_INF, i), 
    //             Get_NAME(R_INF, i), b->self_read.length, b->olist.length);

    // fprintf(stderr, "[M::%s] rid::%ld\n", __func__, i);
    // debug_mm_exact_cigar(&b->olist, i, &b->self_read, &b->ovlp_read);

    // b->num_correct_base += b->olist.length;

    copy_asg_arr(buf0, b->sp); 
    rphase_hc(&b->olist, &R_INF, &b->hap, &b->self_read, &b->ovlp_read, &b->pidx, &b->v64, &buf0, 0, WINDOW_MAX_SIZE, b->self_read.length, 1/**, 0**/, i, (asm_opt.is_ont)?HPC_PL:0, asm_opt.is_ont, ((asm_opt.is_ont)?&(b->clist.chainDP):NULL), ((asm_opt.is_sc)?&(b->v8q):NULL), ((asm_opt.is_sc)?&(b->v8t):NULL));
    copy_asg_arr(b->sp, buf0);

    // stderr_phase_ovlp(&b->olist);

    dedup_chains(&b->olist);

    copy_asg_arr(buf0, b->sp);
    b->cnt[1] += wcns_gen(&b->olist, &R_INF, &b->self_read, &b->ovlp_read, &b->exz, &b->pidx, &b->v64, &buf0, 0, 512, b->self_read.length, 3, 0.500001, aux_o, &b->v32, &b->cns, 256, i);
    copy_asg_arr(b->sp, buf0);

    push_nec_re(aux_o, &(scc.a[i]));
    push_nec_re(aux_o, &(scb.a[i]));

    // if((asm_opt.is_ont) && is_chemical_r_qual(&b->olist, &b->v64, qlen, 1, 16, &(b->v8q), i)/**(is_uncorrected_read(&b->olist, &b->v64, qlen, 1600))**/) {
    //     // b->olist.length = 0;
    //     fprintf(stderr, "[M::%s] rid::%ld\t%.*s\n\n", __func__, i, (int)Get_NAME_LENGTH(R_INF, i), Get_NAME(R_INF, i));
    // }

    push_ne_ovlp(&(R_INF.paf[i]), &b->olist, 1, &R_INF, &(scc.a[i])/**, i, &b->self_read, &b->ovlp_read**/);
    push_ne_ovlp(&(R_INF.reverse_paf[i]), &b->olist, 2, &R_INF, NULL/**, i, NULL, NULL**/);


    check_well_cal(&(scc.a[i]), &b->v64, &(R_INF.paf[i].is_fully_corrected), &(R_INF.paf[i].is_abnormal), qlen, (MIN_COVERAGE_THRESHOLD*2), &(R_INF.paf[i]));
    R_INF.trio_flag[i] = AMBIGU;

    // uint32_t k;
    // for (k = 0; k < b->olist.length; k++) {
    //     if(b->olist.list[k].is_match == 1) b->num_recorrect_base++;
    // }

    // exit(1);
    
    

    // prt_chain(&b->olist);

    // ul_map_lchain(b->abl, (uint32_t)-1, s->seq[i], s->len[i], s->opt->w, s->opt->k, s->uu, &b->olist, &b->clist, s->opt->bw_thres, 
    //         s->opt->max_n_chain, 1, NULL, &(b->tmp_region), NULL, &(b->sp), &high_occ, NULL, 0, 1, 0.2/**0.75**/, 2, 3);

    /**
	int fully_cov, abnormal;
    // if(i != 12578) return;
    // fprintf(stderr, "[M::%s-beg] rid->%ld\n", __func__, i);
    // if (memcmp("7897e875-76e5-42c8-bc37-94b370c4cc8d", Get_NAME((R_INF), i), Get_NAME_LENGTH((R_INF),i)) == 0) {
    //     fprintf(stderr, "[M::%s-beg] rid->%ld\n", __func__, i);
    // } else {
    //     return;
    // }

    ha_get_candidates_interface(b->ab, i, &b->self_read, &b->olist, &b->olist_hp, &b->clist, 
    0.02, asm_opt.max_n_chain, 1, NULL, &b->r_buf, &(R_INF.paf[i]), &(R_INF.reverse_paf[i]), &(b->tmp_region), NULL, &(b->sp));

	clear_Cigar_record(&b->cigar1);
	clear_Round2_alignment(&b->round2);

	correct_overlap(&b->olist, &R_INF, &b->self_read, &b->correct, &b->ovlp_read, &b->POA_Graph, &b->DAGCon,
			&b->cigar1, &b->hap, &b->round2, &b->r_buf, &(b->tmp_region.w_list), 0, 1, &fully_cov, &abnormal);

	b->num_read_base += b->self_read.length;
	b->num_correct_base += b->correct.corrected_base;
	b->num_recorrect_base += b->round2.dumy.corrected_base;

	push_cigar(R_INF.cigars, i, &b->cigar1);
	push_cigar(R_INF.second_round_cigar, i, &b->round2.cigar);

	R_INF.paf[i].is_fully_corrected = 0;
	if (fully_cov) {
		if (get_cigar_errors(&b->cigar1) == 0 && get_cigar_errors(&b->round2.cigar) == 0)
			R_INF.paf[i].is_fully_corrected = 1;
	}
	R_INF.paf[i].is_abnormal = abnormal;

    R_INF.trio_flag[i] = AMBIGU;
    
    ///need to be fixed in r305
    // if(ha_idx_hp == NULL)
    // {
    //     R_INF.trio_flag[i] += collect_hp_regions(&b->olist, &R_INF, &(b->k_flag), RESEED_HP_RATE, Get_READ_LENGTH(R_INF, i), NULL);
    // }

    if (R_INF.trio_flag[i] != AMBIGU || b->save_ov) {
		int is_rev = (asm_opt.number_of_round % 2 == 0);
		push_overlaps(&(R_INF.paf[i]), &b->olist, 1, &R_INF, is_rev);
		push_overlaps(&(R_INF.reverse_paf[i]), &b->olist, 2, &R_INF, is_rev);
	}

    if(het_cnt) het_cnt[i] = get_het_cnt(&b->hap);
    // fprintf(stderr, "[M::%s-end] rid->%ld\n", __func__, i);
    **/
    // exit(1);
    refresh_ec_ovec_buf_t0(b, REFRESH_N);

    /**
    fprintf(stderr, "%ld\t---\n", i);
    **/
}

static void worker_hap_ec_dbg_paf(void *data, long i, int tid)
{
	ec_ovec_buf_t0 *b = &(((cal_ec_r_dbg_step_t*)data)->buf->a[tid]);
    r_dbg_step_res_t *rr = &(((cal_ec_r_dbg_step_t*)data)->res[tid]);
    uint32_t high_occ = asm_opt.hom_cov * (2.0 - HA_KMER_GOOD_RATIO);
    uint32_t low_occ = asm_opt.hom_cov * HA_KMER_GOOD_RATIO;
    overlap_region *aux_o = NULL; i += ((cal_ec_r_dbg_step_t*)data)->si;

    // debug_retrive_bqual(D, &b->v8t, i, 256); return;

    recover_UC_Read(&b->self_read, &R_INF, i); 

    h_ec_lchain(b->ab, i, b->self_read.seq, b->self_read.length, asm_opt.mz_win, asm_opt.k_mer_length, &R_INF, &b->olist, &b->clist, ((asm_opt.is_ont)?(0.05):(0.02)), asm_opt.max_n_chain, 1, NULL, NULL, &(b->sp), &high_occ, &low_occ, 1, 1, 3, 0.7, 2, 32, COV_W);///ONT high error

    aux_o = fetch_aux_ovlp(&b->olist);///must be here

    // stderr_phase_ovlp(&b->olist);
    gen_hc_r_alin_ea(&b->olist, &b->clist, &R_INF, &b->self_read, &b->ovlp_read, &b->exz, aux_o, asm_opt.max_ov_diff_ec, (asm_opt.is_ont)?(WINDOW_OHC):(WINDOW_HC), i, E_KHIT, 1, &b->v16, &b->v64, &(R_INF.paf[i]), 0);

    uint32_t k, m, tl; overlap_region *z; bit_extz_t ez; ma_hit_t *t; 
    for (k = 0; k < b->olist.length; k++) {
        z = &(b->olist.list[k]);
        if(!(z->w_list.n)) continue;

        tl = Get_READ_LENGTH((R_INF), z->y_id);
        for (m = 0; m < z->w_list.n; m++) {
            if(is_ualn_win(z->w_list.a[m])) continue;
            set_bit_extz_t(ez, (*z), m);
            kv_pushp(ma_hit_t, *rr, &t);

            t->qns = z->x_id; t->qns = t->qns << 32;
            t->tn = z->y_id;

            t->qns = t->qns | (uint64_t)(z->w_list.a[m].x_start);
            t->qe = z->w_list.a[m].x_end + 1;

            t->ts = z->w_list.a[m].y_start;
            t->te = z->w_list.a[m].y_end + 1;

            t->rev = z->y_pos_strand;
            
            t->bl = rr->ec.n;
            kv_resize(uint16_t, rr->ec, rr->ec.n + ez.cigar.n);
            memcpy(rr->ec.a + rr->ec.n, ez.cigar.a, ez.cigar.n * sizeof((*(rr->ec.a))));
            rr->ec.n += ez.cigar.n;
            t->cc = rr->ec.n - t->bl;

            if(t->rev) {
                t->ts = tl - z->w_list.a[m].y_end - 1;
                t->te = tl - z->w_list.a[m].y_start;
            }
        }
    }
}

uint32_t adjust_exact_match(asg16_v *in, int64_t xs0, int64_t xe0, int64_t ys0, int64_t ye0, uint64_t *rxs, uint64_t *rxe, uint64_t *rys, uint64_t *rye, uint32_t rev)
{
    *rxs = *rxe = *rys = *rye = 0; 
    if((xe0 <= xs0) || (ye0 <= ys0)) return 0;
    int64_t xk, yk, ck, cn = in->n, wx[2], wy[2], os, oe; uint16_t op, bq, bt; uint32_t cl; uint64_t ovlp;
    xk = yk = 0;

    // fprintf(stderr, "[M::%s]\tx0::[%ld,%ld)\ty0::[%ld,%ld)\trev::%u\n", __func__, xs0, xe0, ys0, ye0, rev);
    
    if(!rev) {
        ck = 0;
        while (ck < cn && xk < xe0) {
            wx[0] = xk; wy[0] = yk; 
            // ck = pop_trace_bp(in, ck, &op, &b, &cl);
            ck = pop_trace_bp_f(in, ck, &op, &bq, &bt, &cl);
            if(op != 2) xk += cl;
            if(op != 3) yk += cl;
            wx[1] = xk; wy[1] = yk;  
            // fprintf(stderr, "[M::%s]\told::[%ld,%ld]\tnew::[%ld,%ld]\t%u%c\n", __func__, wx[0], wx[1], wy[0], wy[1], cl, "MSID"[op]);
            if(op == 0) {
                os = MAX(xs0, wx[0]); oe = MIN(xe0, wx[1]);
                ovlp = ((oe>os)? (oe-os):0);
                if((ovlp > 0) && (ovlp > (*rxe) - (*rxs))) {
                    // fprintf(stderr, "[M::%s]\to::[%ld,%ld)\n", __func__, os, oe);
                    (*rxs) = wy[0] + os - wx[0]; 
                    (*rxe) = wy[0] + oe - wx[0]; 

                    (*rys) = ys0 + os - xs0; 
                    (*rye) = ys0 + oe - xs0; 
                }
            }
            // if((op == 0) && (wx[0] <= s) && (wx[1] >= e)) {
            //     (*rs) = wy[0] + s - wx[0]; 
            //     (*re) = wy[0] + e - wx[0]; 
            // }
        }
    } else {
        ck = cn - 1;
        while (ck >= 0 && xk < xe0) {
            wx[0] = xk; wy[0] = yk; 
            // ck = pop_trace_bp_rev(in, ck, &op, &b, &cl);
            ck = pop_trace_bp_rev_f(in, ck, &op, &bq, &bt, &cl);
            if(op != 2) xk += cl;
            if(op != 3) yk += cl;
            wx[1] = xk; wy[1] = yk;  
            if(op == 0) {
                os = MAX(xs0, wx[0]); oe = MIN(xe0, wx[1]);
                ovlp = ((oe>os)? (oe-os):0);
                if((ovlp > 0) && (ovlp > (*rxe) - (*rxs))) {
                    (*rxs) = wy[0] + os - wx[0]; 
                    (*rxe) = wy[0] + oe - wx[0]; 

                    (*rys) = ys0 + os - xs0; 
                    (*rye) = ys0 + oe - xs0; 
                }
            }
            // if((op == 0) && (wx[0] <= s) && (wx[1] >= e)) {
            //     (*rs) = wy[0] + s - wx[0]; 
            //     (*re) = wy[0] + e - wx[0]; 
            // }
        }
    }
    

    return (*rxe) - (*rxs);
}

uint32_t quick_exact_match(ma_hit_t *z, All_reads *rref, UC_Read* qu, UC_Read* tu, cc_v *sc)
{
    uint64_t rts, rte, rqs, rqe, f = 0; int64_t ql, tl, qr, tr, qs, qe, ts, te;

    // fprintf(stderr, "-0-[M::%s]\tf::%lu\n", __func__, f);
    if(adjust_exact_match(&(sc->a[z->tn]), z->ts, z->te, ((uint32_t)(z->qns)), z->qe, &rts, &rte, &rqs, &rqe, z->rev)) {
        z->ts = rts; z->te = rte; f = 1;
        z->qns >>= 32; z->qns <<= 32; z->qns |= ((uint64_t)(rqs)); z->qe = rqe;

        ///debug
        // qs = rqs; qe = rqe; ts = rts; te = rte; 
        // resize_UC_Read(tu, te - ts); 
        // recover_UC_Read_sub_region(tu->seq, ts, te - ts, z->rev, rref, z->tn);
        
        // fprintf(stderr, "[M::%s]\trq::[%lu,%lu)\trt::[%lu,%lu)\n", __func__, rqs, rqe, rts, rte);
        // fprintf(stderr, "-0-[M::%s] qstr::%.*s\n", __func__, ((int)(qe - qs)), qu->seq + qs);
        // fprintf(stderr, "-0-[M::%s] tstr::%.*s\n", __func__, ((int)(te - ts)), tu->seq);

        // if(memcmp(qu->seq + qs, tu->seq, qe - qs) == 0) {
        //     fprintf(stderr, "-0-[M::%s]\tsb\n", __func__);
        // } else {
        //     fprintf(stderr, "-1-[M::%s]\tsa\n", __func__);
        // }
    }
    // fprintf(stderr, "-1-[M::%s]\tf::%lu\n", __func__, f);
    ql = qu->length; tl = Get_READ_LENGTH((*rref), z->tn);
    qs = ((uint32_t)(z->qns)); qe = z->qe; ts = z->ts; te = z->te;
    if(qs >= ql) qs = ql; if(qe > ql) qe = ql; if(qe <= qs) f = 0;
    // fprintf(stderr, "-2-[M::%s]\tf::%lu\n", __func__, f);
    if(ts >= tl) ts = tl; if(te > tl) te = tl; if(te <= ts) f = 0;
    // fprintf(stderr, "-3-[M::%s]\tf::%lu\n", __func__, f);
    if((qe - qs) != (te - ts)) f = 0;
    // fprintf(stderr, "-4-[M::%s]\tf::%lu\n", __func__, f);

    if(qs <= ts) {
        ts -= qs; qs = 0;
    } else {
        qs -= ts; ts = 0;
    }


    qr = ql - qe; tr = tl - te;
    if(qr <= tr) {
        qe = ql; te += qr;        
    } else {
        te = tl; qe += tr; 
    }

    // fprintf(stderr, "-5-[M::%s]\tzq::[%ld,\t%ld)\tzt::[%ld,\t%ld)\teq::[%u,\t%u)\tet::[%u,\t%u)\tql::%ld\ttl::%ld\tf::%lu\n", 
    // __func__, qs, qe, ts, te, ((uint32_t)(z->qns)), z->qe, z->ts, z->te, ql, tl, f);

    z->qns >>= 32; z->qns <<= 32; z->qns |= ((uint64_t)(qs)); z->qe = qe;
    z->ts = ts; z->te = te;

    if((f) && ((te - ts) == (qe - qs)) && (qe > qs)) {
        resize_UC_Read(tu, te - ts); 
        recover_UC_Read_sub_region(tu->seq, ts, te - ts, z->rev, rref, z->tn);

        // fprintf(stderr, "[M::%s] qstr::%.*s\n", __func__, ((int)(qe - qs)), qu->seq + qs);
        // fprintf(stderr, "[M::%s] tstr::%.*s\n", __func__, ((int)(te - ts)), tu->seq);

        if(memcmp(qu->seq + qs, tu->seq, qe - qs) == 0) return 1;
    }

    return 0;
}

uint64_t cal_cov_re(asg64_v *idx, int64_t *k, uint64_t s, uint64_t e)
{
    uint64_t *a = idx->a, ws, we, cn, os, oe, ovlp, tot = 0; int64_t n = idx->n;
    if(n <= 0) return 0;
    if((*k) >= n) (*k) = 0;
    while (((*k) > 0) && (((uint32_t)a[*k]) > s)) (*k) -= 2;
    if((*k) < 0) (*k) = 0;

    while (((*k) < n) && ((a[*k]>>32) < e)) {
        ws = (a[*k]>>32); we = ((uint32_t)a[*k]); cn = a[(*k) + 1];
        os = MAX(s, ws); oe = MIN(e, we);
        ovlp = ((oe>os)? (oe-os):0);
        tot += (ovlp*cn);
        (*k) += 2;
    }

    return tot;
}

uint64_t gen_hap_dc_cov(asg64_v *be, asg64_v *ba, ma_hit_t_alloc *paf, All_reads *rref, uint64_t wl, int64_t occ_exact, double occ_exact_rate, UC_Read* qu, UC_Read* tu, cc_v *sc, uint64_t rid)
{
    ma_hit_t *z; be->n = ba->n = 0; uint64_t k, s, e, vn, m; asg64_v *v = NULL;
    int64_t dp, old_dp, st = 0, ed, ql = qu->length, qs, qe, ff;
    for (k = 0; k < paf->length; k++) {
        z = &(paf->buffer[k]);

        // if(rid == 16) {
        //     fprintf(stderr, "[M::%s]\tqn::%u::%.*s\ttn::%u::%.*s\t%c\tq::[%u,%u)\tq::[%u,%u)\tel::%u\n", __func__, 
        //         (uint32_t)(z->qns>>32), (int)Get_NAME_LENGTH(R_INF, (uint32_t)(z->qns>>32)), Get_NAME(R_INF, (uint32_t)(z->qns>>32)), z->tn, (int)Get_NAME_LENGTH(R_INF, z->tn), Get_NAME(R_INF, z->tn), "+-"[z->rev],
        //         (uint32_t)z->qns, z->qe, z->ts, z->te, z->el);
        // }
        
        if((z->el) && (quick_exact_match(z, rref, qu, tu, sc))) {
            s = ((uint32_t)(z->qns)); e = z->qe;
            kv_push(uint64_t, (*be), (s<<1));
            kv_push(uint64_t, (*be), (e<<1)|1);
            z->el = 1;
            // fprintf(stderr, "-gmm-[M::%s]\tqn::%u::%.*s\ttn::%u::%.*s\t%c\tq::[%u,%u)\tq::[%u,%u)\n", __func__, 
            //     (uint32_t)(z->qns>>32), (int)Get_NAME_LENGTH(R_INF, (uint32_t)(z->qns>>32)), Get_NAME(R_INF, (uint32_t)(z->qns>>32)), z->tn, (int)Get_NAME_LENGTH(R_INF, z->tn), Get_NAME(R_INF, z->tn), "+-"[z->rev],
            //     (uint32_t)z->qns, z->qe, z->ts, z->te);
        } else {
            s = ((uint32_t)(z->qns)); e = z->qe;
            kv_push(uint64_t, (*ba), (s<<1));
            kv_push(uint64_t, (*ba), (e<<1)|1);
            z->el = 0;
            // fprintf(stderr, "-gum-[M::%s]\tqn::%u::%.*s\ttn::%u::%.*s\t%c\tq::[%u,%u)\tq::[%u,%u)\n", __func__, 
            //     (uint32_t)(z->qns>>32), (int)Get_NAME_LENGTH(R_INF, (uint32_t)(z->qns>>32)), Get_NAME(R_INF, (uint32_t)(z->qns>>32)), z->tn, (int)Get_NAME_LENGTH(R_INF, z->tn), Get_NAME(R_INF, z->tn), "+-"[z->rev],
            //     (uint32_t)z->qns, z->qe, z->ts, z->te);
        }
        /**
        s = ((uint32_t)(z->qns)); e = z->qe;
        // sc->a[z->tn]
        if(z->el) {
            // if((z->qns>>32) == 75 && z->tn == 59) {
            if(quick_exact_match(z, rref, qu, tu, sc)) {
                s = ((uint32_t)(z->qns)); e = z->qe;
                // fprintf(stderr, "-mm-[M::%s]\tqn::%u::%.*s\ttn::%u::%.*s\t%c\n", __func__, 
                // (uint32_t)(z->qns>>32), (int)Get_NAME_LENGTH(R_INF, (uint32_t)(z->qns>>32)), Get_NAME(R_INF, (uint32_t)(z->qns>>32)), z->tn, (int)Get_NAME_LENGTH(R_INF, z->tn), Get_NAME(R_INF, z->tn), "+-"[z->rev]);
            } else {
                s = ((uint32_t)(z->qns)); e = z->qe;
                // fprintf(stderr, "-um-[M::%s]\tqn::%u::%.*s\ttn::%u::%.*s\t%c\n", __func__, 
                // (uint32_t)(z->qns>>32), (int)Get_NAME_LENGTH(R_INF, (uint32_t)(z->qns>>32)), Get_NAME(R_INF, (uint32_t)(z->qns>>32)), z->tn, (int)Get_NAME_LENGTH(R_INF, z->tn), Get_NAME(R_INF, z->tn), "+-"[z->rev]);
            }
            // }
        } 
        **/
    }


    v = be; vn = v->n;
    radix_sort_ec64(v->a, v->a+v->n);
    for (k = 0, dp = 0, st = ed = 0; k < vn; ++k) {
        old_dp = dp;
        ///if a[j] is qe
        if (v->a[k]&1) --dp;
        else ++dp;

        ed = v->a[k]>>1;
        if(ed > st) {
            m = st; m <<= 32; m |= ((uint64_t)ed);
            kv_push(uint64_t, (*v), m); 
            kv_push(uint64_t, (*v), old_dp); 
            if((old_dp + 1) < occ_exact) return 0;///+1 for self
        }
        st = ed;
    }
    ed = ql; old_dp = dp;
    if(ed > st) {
        m = st; m <<= 32; m |= ((uint64_t)ed);
        kv_push(uint64_t, (*v), m); 
        kv_push(uint64_t, (*v), old_dp); 
        if((old_dp + 1) < occ_exact) return 0;///+1 for self
    }

    for (k = vn, m = 0; k < v->n; k++) {
        v->a[m++] = v->a[k];
    }
    v->n = m;


    v = ba; vn = v->n;
    radix_sort_ec64(v->a, v->a+v->n);
    for (k = 0, dp = 0, st = ed = 0; k < vn; ++k) {
        old_dp = dp;
        ///if a[j] is qe
        if (v->a[k]&1) --dp;
        else ++dp;

        ed = v->a[k]>>1;
        if(ed > st) {
            m = st; m <<= 32; m |= ((uint64_t)ed);
            kv_push(uint64_t, (*v), m); 
            kv_push(uint64_t, (*v), old_dp); 
        }
        st = ed;
    }
    ed = ql; old_dp = dp;
    if(ed > st) {
        m = st; m <<= 32; m |= ((uint64_t)ed);
        kv_push(uint64_t, (*v), m); 
        kv_push(uint64_t, (*v), old_dp); 
    }

    for (k = vn, m = 0; k < v->n; k++) {
        v->a[m++] = v->a[k];
    }
    v->n = m;
    






    ///debug
    // v = be; 
    // for (k = 0; k < v->n; k += 2) {
    //     fprintf(stderr, "-be-[M::%s]\tq::[%lu,%u)\tocc::%lu\n", __func__, v->a[k]>>32, (uint32_t)v->a[k], v->a[k+1]);
    // }

    // v = ba; 
    // for (k = 0; k < v->n; k += 2) {
    //     fprintf(stderr, "-ba-[M::%s]\tq::[%lu,%u)\tocc::%lu\n", __func__, v->a[k]>>32, (uint32_t)v->a[k], v->a[k+1]);
    // }









    int64_t ke = 0, ka = 0, cc[2];
    qs = 0; qe = wl; qe = ((qe<=ql)?qe:ql);
    for (; qs < ql; ) {
        cc[0] = cal_cov_re(be, &ke, qs, qe); 
        cc[1] = cal_cov_re(ba, &ka, qs, qe);

        // fprintf(stderr, "[M::%s]\tq::[%ld,%ld)\tcc[0]::%ld\tcc[1]::%ld\n", __func__, qs, qe, cc[0], cc[1]);

        ff = 0;
        if((cc[0]) && (cc[0] > cc[1])) {
            cc[0] += qe - qs;///+(qe - qs) for self
            cc[1] += cc[0];
            if(cc[0] > (cc[1]*occ_exact_rate)) {
                ff = 1;
            }
        }

        if(ff == 0) return 0;
        // fprintf(stderr, "[M::%s-beg]\tq::[%ld,%ld)\tcc[0]::%ld\tcc[1]::%ld\n", __func__, qs, qe, cc[0], cc[1]);
        qs += wl; qe += wl; qe = ((qe<=ql)?qe:ql);
    }


    return 1;
}

static void worker_hap_dc_ec(void *data, long i, int tid)
{
    ec_ovec_buf_t0 *b = &(((ec_ovec_buf_t*)data)->a[tid]);
    // fprintf(stderr, "-0-[M::%s-beg] rid->%ld\n", __func__, i);
    // if (memcmp("m64012_190921_234837/139067658/ccs", Get_NAME((R_INF), i), Get_NAME_LENGTH((R_INF),i)) == 0) {
    //     fprintf(stderr, "-0-[M::%s-beg] rid->%ld\n", __func__, i);
    // } else if (memcmp("m64012_190921_234837/28968323/ccs", Get_NAME((R_INF), i), Get_NAME_LENGTH((R_INF),i)) == 0) {
    //     fprintf(stderr, "-1-[M::%s-beg] rid->%ld\n", __func__, i);
    // } else {
    //     return;
    // }
    // if(i != 2851) return;

    // if(scb.a[i].m < scc.a[i].n) {
    //     scb.a[i].m = scc.a[i].n;
    //     REALLOC(scb.a[i].a, scb.a[i].m);
    // }
    // scb.a[i].n = scc.a[i].n;
    // memcpy(scb.a[i].a, scc.a[i].a, scc.a[i].n*sizeof((*(scb.a[i].a))));

    scc.f[i] = 0;

    if(!(R_INF.paf[i].length)) return;
    // if(scc.f[i]) return;
    asg64_v buf0;

    recover_UC_Read(&b->self_read, &R_INF, i);

    copy_asg_arr(buf0, b->sp);
    if(gen_hap_dc_cov(&(b->v64), &buf0, &(R_INF.paf[i]), &R_INF, WINDOW_HC_FAST, 4, 0.7, &b->self_read, &b->ovlp_read, &scc, i)) {
        scc.f[i] = 1; b->cnt[0]++;
        // fprintf(stderr, "-mm-[M::%s]\tqn::%u::%.*s\n", __func__, (uint32_t)(i), (int)Get_NAME_LENGTH(R_INF, i), Get_NAME((R_INF), i));
    } else {
        scc.f[i] = 0; b->cnt[1]++;
        // fprintf(stderr, "-um-[M::%s]\tqn::%u::%.*s\n", __func__, (uint32_t)(i), (int)Get_NAME_LENGTH(R_INF, i), Get_NAME((R_INF), i));
    }
    copy_asg_arr(b->sp, buf0);

    refresh_ec_ovec_buf_t0(b, REFRESH_N);
}


void flip_paf_rc(uint64_t rid, ma_hit_t_alloc *paf, All_reads *rref)
{
    ma_hit_t *z; uint64_t k, m; int64_t ql = Get_READ_LENGTH((*rref), rid), tl, qs, qe, ts, te;
    for (k = m = 0; k < paf->length; k++) {
        z = &(paf->buffer[k]); tl = Get_READ_LENGTH((*rref), z->tn);
        qs = (uint32_t)z->qns; if(qs < 0) qs = 0; if(qs > ql) qs = ql;
        qe = z->qe; if(qe < 0) qe = 0; if(qe > ql) qe = ql;
        ts = z->ts; if(ts < 0) ts = 0; if(ts > tl) ts = tl;
        te = z->te; if(te < 0) te = 0; if(te > tl) te = tl;
        if(qe > qs && te > ts) {
            z->qns >>= 32; z->qns <<= 32;  
            z->qns |= ((uint64_t)(ql - qe));
            z->qe = ql - qs;
            z->ts = tl - te;
            z->te = tl - ts;
            paf->buffer[m++] = *z;
        }
    }
    paf->length = m;
}

static void worker_hap_post_rev(void *data, long i, int tid)
{
    ec_ovec_buf_t0 *b = &(((ec_ovec_buf_t*)data)->a[tid]);
    uint64_t k, l, kl, nn; char *a, c;
    // fprintf(stderr, "-0-[M::%s-beg] rid->%ld\n", __func__, i);
    // if (memcmp("m64012_190921_234837/139067658/ccs", Get_NAME((R_INF), i), Get_NAME_LENGTH((R_INF),i)) == 0) {
    //     fprintf(stderr, "-0-[M::%s-beg] rid->%ld\n", __func__, i);
    // } else if (memcmp("m64012_190921_234837/28968323/ccs", Get_NAME((R_INF), i), Get_NAME_LENGTH((R_INF),i)) == 0) {
    //     fprintf(stderr, "-1-[M::%s-beg] rid->%ld\n", __func__, i);
    // } else {
    //     return;
    // }
    // if(i != 2851) return;

    // if(scb.a[i].m < scc.a[i].n) {
    //     scb.a[i].m = scc.a[i].n;
    //     REALLOC(scb.a[i].a, scb.a[i].m);
    // }
    // scb.a[i].n = scc.a[i].n;
    // memcpy(scb.a[i].a, scc.a[i].a, scc.a[i].n*sizeof((*(scb.a[i].a))));

    flip_paf_rc(i, &(R_INF.paf[i]), &R_INF);
    flip_paf_rc(i, &(R_INF.reverse_paf[i]), &R_INF);

    recover_UC_Read(&b->self_read, &R_INF, i); 
    l = b->self_read.length; kl = l>>1; a = b->self_read.seq;
    for (k = nn = 0; k < kl; k++) {
        c = a[l-k-1]; a[l-k-1] = RC_CHAR(a[k]); a[k] = RC_CHAR(c);
        if(a[k] == 'N') nn++;
        if(a[l-k-1] == 'N') nn++;
    }
    if(l&1) {
        a[k] = RC_CHAR(a[k]); 
        if(a[k] == 'N') nn++;
    }

    ha_compress_base(Get_READ(R_INF, i), a, l, &R_INF.N_site[i], nn);

    if(asm_opt.is_sc) {
        retrive_bqual(&(b->v8q), NULL, i, -1, -1, 0, sc_bn);
        for (k = 0; k < l; k++) a[l - k - 1] = b->v8q.a[k];
        ha_compress_qual_bit(Get_QUAL(R_INF, i), a, l, sc_bn);
    }
}

static void worker_hap_dc_ec_gen(void *data, long i, int tid)
{
    
    ec_ovec_buf_t0 *b = &(((ec_ovec_buf_t*)data)->a[tid]);
    uint32_t high_occ = asm_opt.hom_cov * (2.0 - HA_KMER_GOOD_RATIO);
    uint32_t low_occ = asm_opt.hom_cov * HA_KMER_GOOD_RATIO;

    recover_UC_Read(&b->self_read, &R_INF, i); 

    // overlap_region_sort_y_id(b->olist.list, b->olist.length);
	// ma_hit_sort_tn(R_INF.paf[i].buffer, R_INF.paf[i].length);
	// ma_hit_sort_tn(R_INF.reverse_paf[i].buffer, R_INF.reverse_paf[i].length);

    // R_INF.paf[i].is_fully_corrected = is_well_cal(&b->v64, &(R_INF.paf[i]), &(R_INF.reverse_paf[i]), b->self_read.length, 4);

    // R_INF.paf[i].is_abnormal = abnormal;
    // R_INF.trio_flag[i] = AMBIGU;

    h_ec_lchain_fast(b->ab, i, &b->self_read, &b->ovlp_read, asm_opt.mz_win, asm_opt.k_mer_length, &R_INF, &b->olist, &b->clist, &b->exz, &b->v16, &b->v64, 0.02, 1, NULL, NULL, &(b->sp), &high_occ, &low_occ, 1, 1, 0, 2, UINT32_MAX, &(R_INF.paf[i]), &(R_INF.reverse_paf[i]), 0.866666);

    push_ff_ovlp(&(R_INF.paf[i]), &b->olist, 1, &R_INF, b->cnt);
    push_ff_ovlp(&(R_INF.reverse_paf[i]), &b->olist, 2, &R_INF, b->cnt);

    /**
    copy_asg_arr(buf0, b->sp);
    if(gen_hap_dc_cov(&(b->v64), &buf0, &(R_INF.paf[i]), &R_INF, WINDOW_HC_FAST, 4, 0.7, &b->self_read, &b->ovlp_read, &scc, i)) {
        scc.f[i] = 1; b->num_read_base++;
        // fprintf(stderr, "-mm-[M::%s]\tqn::%u::%.*s\n", __func__, (uint32_t)(i), (int)Get_NAME_LENGTH(R_INF, i), Get_NAME((R_INF), i));
    } else {
        scc.f[i] = 0; b->num_correct_base++;
        // fprintf(stderr, "-um-[M::%s]\tqn::%u::%.*s\n", __func__, (uint32_t)(i), (int)Get_NAME_LENGTH(R_INF, i), Get_NAME((R_INF), i));
    }
    copy_asg_arr(b->sp, buf0);
    **/
   refresh_ec_ovec_buf_t0(b, REFRESH_N);
}

static void worker_hap_dc_ec_gen_new_idx(void *data, long i, int tid)
{
    
    ec_ovec_buf_t0 *b = &(((ec_ovec_buf_t*)data)->a[tid]);
    uint32_t high_occ = asm_opt.hom_cov * (2.0 - HA_KMER_GOOD_RATIO);
    uint32_t low_occ = asm_opt.hom_cov * HA_KMER_GOOD_RATIO; uint32_t qlen = 0;

    recover_UC_Read(&b->self_read, &R_INF, i); qlen = b->self_read.length; 

    h_ec_lchain(b->ab, i, b->self_read.seq, b->self_read.length, asm_opt.mz_win, asm_opt.k_mer_length, &R_INF, &b->olist, &b->clist, /**0.02**/0.001, asm_opt.max_n_chain, 1, NULL, NULL, &(b->sp), &high_occ, &low_occ, 1, 1, 3, 0.7, 2, 32, COV_W);

    overlap_region_sort_y_id(b->olist.list, b->olist.length);

    // overlap_region_sort_y_id(b->olist.list, b->olist.length);
	// ma_hit_sort_tn(R_INF.paf[i].buffer, R_INF.paf[i].length);
	// ma_hit_sort_tn(R_INF.reverse_paf[i].buffer, R_INF.reverse_paf[i].length);

    // R_INF.paf[i].is_fully_corrected = is_well_cal(&b->v64, &(R_INF.paf[i]), &(R_INF.reverse_paf[i]), b->self_read.length, 4);

    // R_INF.paf[i].is_abnormal = abnormal;
    // R_INF.trio_flag[i] = AMBIGU;

    h_ec_lchain_fast_new(b->ab, i, &b->self_read, &b->ovlp_read, &R_INF, &b->olist, &b->clist, &b->exz, &b->v16, &b->v64, &(R_INF.paf[i]), &(R_INF.reverse_paf[i]), 0.866666);

    if((asm_opt.is_ont) && (is_uncorrected_read(&b->olist, &b->v64, qlen, 1600))) {
        b->olist.length = 0;
        // fprintf(stderr, "[M::%s] rid::%ld\t%.*s\n", __func__, i, (int)Get_NAME_LENGTH(R_INF, i), Get_NAME(R_INF, i));
    }

    push_ff_ovlp(&(R_INF.paf[i]), &b->olist, 1, &R_INF, b->cnt);
    push_ff_ovlp(&(R_INF.reverse_paf[i]), &b->olist, 2, &R_INF, b->cnt);

    /**
    copy_asg_arr(buf0, b->sp);
    if(gen_hap_dc_cov(&(b->v64), &buf0, &(R_INF.paf[i]), &R_INF, WINDOW_HC_FAST, 4, 0.7, &b->self_read, &b->ovlp_read, &scc, i)) {
        scc.f[i] = 1; b->num_read_base++;
        // fprintf(stderr, "-mm-[M::%s]\tqn::%u::%.*s\n", __func__, (uint32_t)(i), (int)Get_NAME_LENGTH(R_INF, i), Get_NAME((R_INF), i));
    } else {
        scc.f[i] = 0; b->num_correct_base++;
        // fprintf(stderr, "-um-[M::%s]\tqn::%u::%.*s\n", __func__, (uint32_t)(i), (int)Get_NAME_LENGTH(R_INF, i), Get_NAME((R_INF), i));
    }
    copy_asg_arr(b->sp, buf0);
    **/
   refresh_ec_ovec_buf_t0(b, REFRESH_N);
}

uint32_t is_chemical_r(ma_hit_t_alloc *ov, asg64_v *idx, int64_t len, int64_t cov, int64_t frank_len)
{
    uint64_t k, s, e; int64_t dp, old_dp, st = 0, ed, s0, s1, e0, e1, rr;
    if((frank_len) > (len *0.01)) frank_len = len *0.01;
    for (k = idx->n = 0; k < ov->length; k++) {
        s = (uint32_t)ov->buffer[k].qns; e = ov->buffer[k].qe;
        kv_push(uint64_t, (*idx), (s<<1));
        kv_push(uint64_t, (*idx), (e<<1)|1);
    }

    radix_sort_ec64(idx->a, idx->a + idx->n); s0 = s1 = e0 = e1 = rr = -1;
    for (k = 0, dp = 0, st = ed = 0; k < idx->n; ++k) {
        old_dp = dp;
        ///if a[j] is qe
        if (idx->a[k]&1) --dp;
        else ++dp;

        ed = idx->a[k]>>1;
        if(ed > st) {
            if(old_dp <= cov) {
                rr = 1;
            } else {
                if(s0 < 0) {
                    s0 = st; s1 = ed;
                }
                e0 = st; e1 = ed;
            }
        }
        st = ed;
    }


    ed = len; old_dp = dp;
    if(ed > st) {
        if(old_dp <= cov) {
            rr = 1;
        } else {
            if(s0 < 0) {
                s0 = st; s1 = ed;
            }
            e0 = st; e1 = ed;
        }
    }


    if((s0 != e0) && (s1 != e1) && (s0 <= frank_len) && ((len - e1) <= frank_len) && (rr > 0)) {
        for (k = 0, dp = 0, st = ed = 0; k < idx->n; ++k) {
            old_dp = dp;
            ///if a[j] is qe
            if (idx->a[k]&1) --dp;
            else ++dp;

            ed = idx->a[k]>>1;
            if(ed > st) {
                if((old_dp <= cov) && (st >= s0) && (ed <= e1)) {
                    // if((ov->buffer[0].qns>>32) == 3364) fprintf(stderr, "[M::%s]\tlf::[%ld,%ld)\trt::[%ld,%ld)\tmd::[%ld,%ld)\tcov::%ld\n", __func__, s0, s1, e0, e1, st, ed, old_dp);
                    return 1;
                }
            }
            st = ed;
        }


        ed = len; old_dp = dp;
        if(ed > st) {
            if((old_dp <= cov) && (st >= s0) && (ed <= e1)) {
                // if((ov->buffer[0].qns>>32) == 3364) fprintf(stderr, "[M::%s]\tlf::[%ld,%ld)\trt::[%ld,%ld)\tmd::[%ld,%ld)\tcov::%ld\n", __func__, s0, s1, e0, e1, st, ed, old_dp);
                return 1;
            }
        }
    }
    
    return 0;
}


uint32_t is_chemical_r_adv(ma_hit_t_alloc *ov, asg64_v *idx, int64_t len, int64_t cov, int64_t cut_len, double dup_rate, uint64_t is_del)
{
    uint64_t k, s, e; int64_t dp, old_dp, st = 0, ed, s0, e0, rr, lt;
    for (k = idx->n = 0; k < ov->length; k++) {
        if(is_del && ov->buffer[k].del) continue;
        s0 = (uint32_t)ov->buffer[k].qns; e0 = ov->buffer[k].qe;
        if(s0 > 0) s0 += cut_len;
        if(e0 < len) e0 -= cut_len;
        if(e0 <= s0) continue;
        s = s0; e = e0;

        lt = Get_READ_LENGTH((R_INF), ov->buffer[k].tn);
        rr = (lt >= len)?(lt - len):(len - lt);
        if((rr <= (len*dup_rate)) && (rr <= (lt*dup_rate)) && (ov->buffer[k].rev)) {
            dp = (ov->buffer[k].qe) - ((uint32_t)ov->buffer[k].qns); dp = len - dp;
            old_dp = ov->buffer[k].te - ov->buffer[k].ts; old_dp = lt - old_dp;
            if((dp <= (len*dup_rate)) && (old_dp <= (lt*dup_rate))) continue;
        } 

        kv_push(uint64_t, (*idx), (s<<1));
        kv_push(uint64_t, (*idx), (e<<1)|1);
    }

    radix_sort_ec64(idx->a, idx->a + idx->n); s0 = e0 = rr = -1;
    for (k = 0, dp = 0, st = ed = 0; k < idx->n; ++k) {
        old_dp = dp;
        ///if a[j] is qe
        if (idx->a[k]&1) --dp;
        else ++dp;

        ed = idx->a[k]>>1;
        if(ed > st) {
            // if(ov->length && ((ov->buffer[0].qns>>32) == 5045637)) {
            //     fprintf(stderr, "[M::%s]\tmd::[%ld,%ld)\tcov::%ld\tlen::%ld\tid::%lu\n", __func__, st, ed, old_dp, len, ov->buffer[0].qns>>32);
            // }
            if(old_dp <= cov) {
                // if(ov->length && (ov->buffer[0].qns>>32) == 22344) fprintf(stderr, "[M::%s]\tmd::[%ld,%ld)\tcov::%ld\tlen::%ld\n", __func__, st, ed, old_dp, len);
                return 1;
            } 
        }
        st = ed;
    }


    ed = len; old_dp = dp;
    if(ed > st) {
        // if(ov->length && ((ov->buffer[0].qns>>32) == 5045637)) {
        //     fprintf(stderr, "[M::%s]\tmd::[%ld,%ld)\tcov::%ld\tlen::%ld\tid::%lu\n", __func__, st, ed, old_dp, len, ov->buffer[0].qns>>32);
        // }
        if(old_dp <= cov) {
            // if(ov->length && (ov->buffer[0].qns>>32) == 22344) fprintf(stderr, "[M::%s]\tmd::[%ld,%ld)\tcov::%ld\tlen::%ld\n", __func__, st, ed, old_dp, len);
            return 1;
        }
    }
    
    return 0;
}

int64_t cal_chemical_r_adv(ma_hit_t_alloc *ov, asg64_v *idx, int64_t len, int64_t cut_len, double dup_rate, uint64_t is_del)
{
    uint64_t k, s, e; int64_t dp, old_dp, st = 0, ed, s0, e0, rr, lt, min_cov;
    for (k = idx->n = 0; k < ov->length; k++) {
        if(is_del && ov->buffer[k].del) continue;
        s0 = (uint32_t)ov->buffer[k].qns; e0 = ov->buffer[k].qe;
        if(s0 > 0) s0 += cut_len;
        if(e0 < len) e0 -= cut_len;
        if(e0 <= s0) continue;
        s = s0; e = e0;

        lt = Get_READ_LENGTH((R_INF), ov->buffer[k].tn);
        rr = (lt >= len)?(lt - len):(len - lt);
        if((rr <= (len*dup_rate)) && (rr <= (lt*dup_rate)) && (ov->buffer[k].rev)) {
            dp = (ov->buffer[k].qe) - ((uint32_t)ov->buffer[k].qns); dp = len - dp;
            old_dp = ov->buffer[k].te - ov->buffer[k].ts; old_dp = lt - old_dp;
            if((dp <= (len*dup_rate)) && (old_dp <= (lt*dup_rate))) continue;
        } 

        kv_push(uint64_t, (*idx), (s<<1));
        kv_push(uint64_t, (*idx), (e<<1)|1);
    }

    radix_sort_ec64(idx->a, idx->a + idx->n); s0 = e0 = rr = -1; min_cov = INT64_MAX;
    for (k = 0, dp = 0, st = ed = 0; k < idx->n; ++k) {
        old_dp = dp;
        ///if a[j] is qe
        if (idx->a[k]&1) --dp;
        else ++dp;

        ed = idx->a[k]>>1;
        if(ed > st) {
            // if(ov->length && ((ov->buffer[0].qns>>32) == 5045637)) {
            //     fprintf(stderr, "[M::%s]\tmd::[%ld,%ld)\tcov::%ld\tlen::%ld\tid::%lu\n", __func__, st, ed, old_dp, len, ov->buffer[0].qns>>32);
            // }
            if(old_dp <= min_cov) {
                // if(ov->length && (ov->buffer[0].qns>>32) == 22344) fprintf(stderr, "[M::%s]\tmd::[%ld,%ld)\tcov::%ld\tlen::%ld\n", __func__, st, ed, old_dp, len);
                min_cov = old_dp;
            } 
        }
        st = ed;
    }


    ed = len; old_dp = dp;
    if(ed > st) {
        // if(ov->length && ((ov->buffer[0].qns>>32) == 5045637)) {
        //     fprintf(stderr, "[M::%s]\tmd::[%ld,%ld)\tcov::%ld\tlen::%ld\tid::%lu\n", __func__, st, ed, old_dp, len, ov->buffer[0].qns>>32);
        // }
        if(old_dp <= min_cov) {
            // if(ov->length && (ov->buffer[0].qns>>32) == 22344) fprintf(stderr, "[M::%s]\tmd::[%ld,%ld)\tcov::%ld\tlen::%ld\n", __func__, st, ed, old_dp, len);
            min_cov = old_dp;
        }
    }
    
    return min_cov;
}

void prt_dbg_rid_paf(ma_hit_t_alloc *ov, UC_Read *ra, asg8_v *qa)
{
    if(!(ov->length)) return;
    uint64_t k, qn = (ov->buffer[0].qns>>32), qn_n, i, m; char *nn = NULL; FILE *fp = NULL; ma_hit_t *h = NULL;
    qn_n = Get_NAME_LENGTH((R_INF), qn) + 64; MALLOC(nn, qn_n);

    sprintf(nn, "%.*s.qry.fq", (int)Get_NAME_LENGTH(R_INF, qn), Get_NAME((R_INF), qn)); fp = fopen(nn, "w");
    for (k = 0; k < ov->length; k++) {
        i = ov->buffer[k].tn;
        recover_UC_Read(ra, &R_INF, i);
        fprintf(fp, "@%.*s\n", (int32_t)Get_NAME_LENGTH(R_INF, i), Get_NAME(R_INF, i));
        fprintf(fp, "%.*s\n", (int32_t)ra->length, ra->seq);        
        fprintf(fp, "+\n");
        retrive_bqual(qa, NULL, i, -1, -1, 0, sc_bn);
        for (m = 0; m < qa->n; m++) fprintf(fp, "%c", (char)(sc_tb[qa->a[m]] + 33 - 1));
        fprintf(fp, "\n");
    }
    fclose(fp);

    sprintf(nn, "%.*s.ref.fq", (int)Get_NAME_LENGTH(R_INF, qn), Get_NAME((R_INF), qn)); fp = fopen(nn, "w");
    i = qn;
    recover_UC_Read(ra, &R_INF, i);
    fprintf(fp, "@%.*s\n", (int32_t)Get_NAME_LENGTH(R_INF, i), Get_NAME(R_INF, i));
    fprintf(fp, "%.*s\n", (int32_t)ra->length, ra->seq);        
    fprintf(fp, "+\n");
    retrive_bqual(qa, NULL, i, -1, -1, 0, sc_bn);
    for (m = 0; m < qa->n; m++) fprintf(fp, "%c", (char)(sc_tb[qa->a[m]] + 33 - 1));
    fprintf(fp, "\n");
    fclose(fp);

    sprintf(nn, "%.*s.ref.fa", (int)Get_NAME_LENGTH(R_INF, qn), Get_NAME((R_INF), qn)); fp = fopen(nn, "w");
    i = qn;
    recover_UC_Read(ra, &R_INF, i);
    fprintf(fp, ">%.*s\n", (int32_t)Get_NAME_LENGTH(R_INF, i), Get_NAME(R_INF, i));
    fprintf(fp, "%.*s\n", (int32_t)ra->length, ra->seq);        
    // fprintf(fp, "+\n");
    // retrive_bqual(qa, NULL, i, -1, -1, 0, sc_bn);
    // for (m = 0; m < qa->n; m++) fprintf(fp, "%c", (char)(sc_tb[qa->a[m]] + 33 - 1));
    // fprintf(fp, "\n");
    fclose(fp);

    sprintf(nn, "%.*s.ov.paf", (int)Get_NAME_LENGTH(R_INF, qn), Get_NAME((R_INF), qn)); fp = fopen(nn, "w");
    for (k = 0; k < ov->length; k++) {
        h = &(ov->buffer[k]);
        fprintf(fp, "%.*s(qn::%u)\t%u\t%u\t%u\t%c\t%.*s(tn::%u)\t%u\t%u\t%u\t%u\t%u\t255\n", (int)Get_NAME_LENGTH(R_INF, Get_qn(*h)), Get_NAME((R_INF), Get_qn(*h)), Get_qn(*h), (uint32_t)Get_READ_LENGTH(R_INF, Get_qn(*h)), Get_qs(*h), Get_qe(*h), "+-"[h->rev], 
        (int)Get_NAME_LENGTH(R_INF, Get_tn(*h)), Get_NAME((R_INF), Get_tn(*h)), Get_tn(*h), (uint32_t)Get_READ_LENGTH(R_INF, Get_tn(*h)), Get_ts(*h), Get_te(*h), h->ml, h->bl);
    }
    fclose(fp);

    free(nn);
}

static void worker_hap_dc_ec_chemical_r(void *data, long i, int tid)
{
    ec_ovec_buf_t0 *b = &(((ec_ovec_buf_t*)data)->a[tid]);
    ma_hit_t_alloc *paf = &(R_INF.paf[i]); uint64_t k, m;

    // if (memcmp("3ed80bc4-1169-4948-a9ff-9c2463b7f7a2", Get_NAME((R_INF), i), Get_NAME_LENGTH((R_INF),i)) == 0) {
    //     fprintf(stderr, "-a-[M::%s-beg] rid->%ld, b->rr->%lu\n", __func__, i, b->rr);
    // } 
    if(b->cnt[1] == 0) {
        // if(i == 6204620) prt_dbg_rid_paf(&(R_INF.paf[i]), &(b->self_read), &(b->v8q));
        // if(is_chemical_r(&(R_INF.paf[i]), &b->v64, Get_READ_LENGTH((R_INF), i), 3, 16)) {
        if(is_chemical_r_adv(&(R_INF.paf[i]), &b->v64, Get_READ_LENGTH((R_INF), i), asm_opt.chemical_cov, asm_opt.chemical_flank, 0.02, 0)) {
            // fprintf(stderr, "-um-[M::%s]\tqn::%u::%.*s\n\n", __func__, (uint32_t)(i), (int)Get_NAME_LENGTH(R_INF, i), Get_NAME((R_INF), i));
            R_INF.paf[i].length = 0; b->cnt[0]++;
        } 
    } else if(b->cnt[1] == 1) {
        for (k = 0; k < paf->length; k++) {
            if(R_INF.paf[paf->buffer[k].tn].length == 0) {
                paf->buffer[k].tn = (uint32_t)-1; b->cnt[0]++;
            }
        }
    } else {
        for (k = m = 0; k < paf->length; k++) {
            if(paf->buffer[k].tn == ((uint32_t)-1)) continue;
            paf->buffer[m++] = paf->buffer[k];
        }
        paf->length = m;
    }

   refresh_ec_ovec_buf_t0(b, REFRESH_N);
}

static void worker_hap_dc_ec_chemical_arc(void *data, long i, int tid)
{
    ec_ovec_buf_t0 *b = &(((ec_ovec_buf_t*)data)->a[tid]);
    ma_hit_t_alloc *paf = &(R_INF.paf[i]), *rev; uint64_t k, z;

    if(b->cnt[1] == 0) {
        if(is_chemical_r_adv(&(R_INF.paf[i]), &b->v64, Get_READ_LENGTH((R_INF), i), asm_opt.chemical_cov, asm_opt.chemical_flank, 0.02, 1)) {
            // fprintf(stderr, "-um-[M::%s]\tqn::%u::%.*s\n\n", __func__, (uint32_t)(i), (int)Get_NAME_LENGTH(R_INF, i), Get_NAME((R_INF), i));   
            for (k = 0; k < paf->length; k++) paf->buffer[k].del = 1; b->cnt[0]++;
        }
    } else if(b->cnt[1] == 1) {
        for (k = 0; k < paf->length; k++) {
            if((Get_qn(paf->buffer[k])) > (Get_tn(paf->buffer[k]))) continue;
            rev = &(R_INF.paf[paf->buffer[k].tn]);
            for (z = 0; z < rev->length; z++) {
                if((rev->buffer[z].tn == (Get_qn(paf->buffer[k])))) {
                    if(paf->buffer[k].del != rev->buffer[z].del) {
                        paf->buffer[k].del = rev->buffer[z].del = 1;
                    }
                }
            }
        }
    }

   refresh_ec_ovec_buf_t0(b, REFRESH_N);
}

static void worker_hap_dc_ec_chemical_arc_mark(void *data, long i, int tid)
{
    ec_ovec_buf_t0 *b = &(((ec_ovec_buf_t*)data)->a[tid]);
    ma_hit_t_alloc *paf = &(R_INF.paf[i]), *rev; uint64_t k, z; int64_t cov, msk_cut = asm_opt.chemical_cov;
    uint8_t *msk = ((ec_ovec_buf_t*)data)->cr;

    if(b->cnt[1] == 0) {
        msk[i] = (uint8_t)-1;
        cov = cal_chemical_r_adv(&(R_INF.paf[i]), &b->v64, Get_READ_LENGTH((R_INF), i), asm_opt.chemical_flank, 0.02, 1);
        if(cov <= msk_cut) msk[i] = cov;
        if(cov <= msk_cut/**FORCE_CUT**/) {
            // fprintf(stderr, "-um-[M::%s]\tqn::%u::%.*s\n\n", __func__, (uint32_t)(i), (int)Get_NAME_LENGTH(R_INF, i), Get_NAME((R_INF), i));   
            for (k = 0; k < paf->length; k++) paf->buffer[k].del = 1; b->cnt[0]++;
        }
    } else if(b->cnt[1] == 1) {
        for (k = 0; k < paf->length; k++) {
            if((Get_qn(paf->buffer[k])) > (Get_tn(paf->buffer[k]))) continue;
            rev = &(R_INF.paf[paf->buffer[k].tn]);
            for (z = 0; z < rev->length; z++) {
                if((rev->buffer[z].tn == (Get_qn(paf->buffer[k])))) {
                    if((paf->buffer[k].del != rev->buffer[z].del) || (msk[Get_qn(paf->buffer[k])] <= msk_cut/**FORCE_CUT**/) || (msk[Get_tn(paf->buffer[k])] <= msk_cut/**FORCE_CUT**/)) {
                        paf->buffer[k].del = rev->buffer[z].del = 1;
                    }
                }
            }
        }
    }

   refresh_ec_ovec_buf_t0(b, REFRESH_N);
}

void gen_ovlst_paf(ma_hit_t_alloc *in_e, ma_hit_t_alloc *in_r, asg64_v *ou)
{
    uint32_t n = 0, k; 

    for (k = 0; k < in_e->length; k++) {
        if(!(in_e->buffer[k].el)) n++;
    }
    n += in_r->length;

    kv_resize(uint64_t, *ou, n); ou->n = 0;

    for (k = 0; k < in_e->length; k++) {
        if(!(in_e->buffer[k].el)) {
            ou->a[ou->n] = in_e->buffer[k].tn;
            ou->a[ou->n] <<= 1; ou->a[ou->n] |= in_e->buffer[k].rev;
            ou->n++;
        }
    }

    for (k = 0; k < in_r->length; k++) {
        ou->a[ou->n] = in_r->buffer[k].tn;
        ou->a[ou->n] <<= 1; ou->a[ou->n] |= in_r->buffer[k].rev;
        ou->n++;
    }

    radix_sort_ec64(ou->a, ou->a+ou->n);
}

void dbg_overlap_region_cigar(overlap_region *a, uint64_t a_n, char *qstr, All_reads *rref, UC_Read *tu)
{
    bit_extz_t ez; uint64_t i, k;
    for (i = 0; i < a_n; i++) {
        if(a[i].y_pos_strand) {
            recover_UC_Read_RC(tu, rref, a[i].y_id);
        } else {
            recover_UC_Read(tu, rref, a[i].y_id);
        }
        for (k = 0; k < a[i].w_list.n; k++) {
            if(is_ualn_win((a[i].w_list.a[k]))) continue;
            set_bit_extz_t(ez, a[i], k);
            if(!cigar_check(tu->seq, qstr, &ez)) {
                fprintf(stderr, "\n-0-[M::%s] x_id::%u, y_id::%u, x::[%u, %u), y::[%u, %u)\n", __func__, a[i].x_id, a[i].y_id, a[i].x_pos_s, a[i].x_pos_e + 1, a[i].y_pos_s, a[i].y_pos_e + 1);
                exit(1);
            } else {
                // fprintf(stderr, "\n-1-[M::%s] x_id::%u, y_id::%u, x::[%u, %u), y::[%u, %u)\n", __func__, a[i].x_id, a[i].y_id, a[i].x_pos_s, a[i].x_pos_e + 1, a[i].y_pos_s, a[i].y_pos_e + 1);
            }
        }
    }
}

overlap_region* h_ec_lchain_re(ha_abuf_t *ab, uint32_t rid, char* rs, uint64_t rl, UC_Read *tu, uint64_t mz_w, uint64_t mz_k, All_reads *rref, overlap_region_alloc *ol, Candidates_list *cl, bit_extz_t *exz, asg16_v* buf, double bw_thres, 
								 int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t is_accurate, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, ma_hit_t_alloc *in0, ma_hit_t_alloc *in1)
{
    // fprintf(stderr, "-mm-[M::%s]\tchain_cutoff::%u\n", __func__, chain_cutoff);
    uint64_t on = 0, k, ol0, wl = (asm_opt.is_ont)?(WINDOW_OHC):(WINDOW_HC), m; ma_hit_t *oa = NULL; overlap_region *z = NULL, *aux_o = NULL, t; Window_Pool w; double err = asm_opt.max_ov_diff_ec;

    int64_t max_skip, max_iter, max_dis, quick_check; double chn_pen_gap, chn_pen_skip; 
	set_lchain_dp_op(is_accurate, mz_k, &max_skip, &max_iter, &max_dis, &chn_pen_gap, &chn_pen_skip, &quick_check);

    init_Window_Pool(&w, rl, wl, (int)(1.0/err));

    on = in0->length + in1->length + 1;
    clear_overlap_region_alloc(ol); 
    clear_Candidates_list(cl); 
    if(on > ol->size) {
        REALLOC(ol->list, on);
        memset(ol->list+ol->size, 0, sizeof(overlap_region)*(on-ol->size));
        ol->size = on;   
    }
    on = in0->length + in1->length; aux_o = &(ol->list[on]);
    ol->length = 0; ol->mapped_overlaps_length = 0; m = 0;
    

    // get the list of anchors
    get_mz1(rs, rl, mz_w, mz_k, 0, !(asm_opt.flag & HA_F_NO_HPC), ab, ha_flt_tab, ha_idx, asm_opt.mz_sample_dist, k_flag, dbg_ct, NULL, -1, asm_opt.dp_min_len, -1, sp, asm_opt.mz_rewin, 0, NULL, 0);

    oa = in0->buffer; on = in0->length;
    for (k = 0; k < on; k++) {
        ol0 = ol->length;
        if(oa[k].el) {
            z = &(ol->list[ol->length++]);
            z->x_id = rid; z->y_id = oa[k].tn;
            z->x_pos_strand = 0; z->y_pos_strand = oa[k].rev;
            z->x_pos_s = (uint32_t)oa[k].qns;
            z->x_pos_e = oa[k].qe - 1;
            z->y_pos_s = oa[k].ts;
            z->y_pos_e = oa[k].te - 1;

            z->is_match = 1;
            z->align_length = z->overlapLen = z->shared_seed = z->x_pos_e + 1 - z->x_pos_s;
            z->non_homopolymer_errors = z->strong = 0;

            set_exact_exz(exz, z->x_pos_s, z->x_pos_e + 1, z->y_pos_s, z->y_pos_e + 1); push_alnw(z, exz);

            // if(oa[k].tn == 15382) fprintf(stderr, "-em-[M::%s]\tqn::%u\ttn::%u\terr::%u\n", __func__, rid, oa[k].tn, ol->list[ol->length-1].non_homopolymer_errors);
        } else {
            if(oa[k].rev) recover_UC_Read_RC(tu, rref, oa[k].tn);
            else recover_UC_Read(tu, rref, oa[k].tn);

            get_pi_ec_chain(ab, rid, rl, oa[k].tn, tu->seq, tu->length, mz_w, mz_k, ol, cl, bw_thres, apend_be, k_flag, dbg_ct, sp, high_occ, low_occ, gen_off, enable_mcopy, mcopy_rate, mcopy_khit_cut, max_skip, max_iter, max_dis, quick_check, chn_pen_gap, chn_pen_skip);
            assert(ol->length - ol0 <= 1);
            if((ol->length > ol0) && (gen_hc_r_alin_re(&(ol->list[ol0]), cl, rs, rl, tu->seq, tu->length, exz, aux_o, asm_opt.max_ov_diff_ec, w.window_length, rid, E_KHIT, 1, buf))) {
                ol->list[ol0].y_pos_strand = oa[k].rev;
                // if(oa[k].tn == 15382) fprintf(stderr, "-mm-[M::%s]\tqn::%u\ttn::%u\terr::%u\n", __func__, rid, oa[k].tn, ol->list[ol->length-1].non_homopolymer_errors);
            } 
        }

        if((ol->length <= ol0) || (ol->list[ol0].is_match != 1)) continue;

        if(m != ol0) {
            t = ol->list[m];
            ol->list[m] = ol->list[ol0];
            ol->list[ol0] = t;
        }

        m++;
    }


    oa = in1->buffer; on = in1->length;
    for (k = 0; k < on; k++) {
        ol0 = ol->length;

        if(oa[k].rev) recover_UC_Read_RC(tu, rref, oa[k].tn);
        else recover_UC_Read(tu, rref, oa[k].tn);

        get_pi_ec_chain(ab, rid, rl, oa[k].tn, tu->seq, tu->length, mz_w, mz_k, ol, cl, bw_thres, apend_be, k_flag, dbg_ct, sp, high_occ, low_occ, gen_off, enable_mcopy, mcopy_rate, mcopy_khit_cut, max_skip, max_iter, max_dis, quick_check, chn_pen_gap, chn_pen_skip);
        assert(ol->length - ol0 <= 1);
        if((ol->length > ol0) && (gen_hc_r_alin_re(&(ol->list[ol0]), cl, rs, rl, tu->seq, tu->length, exz, aux_o, asm_opt.max_ov_diff_ec, w.window_length, rid, E_KHIT, 1, buf))) {
            ol->list[ol0].y_pos_strand = oa[k].rev;
            // if(oa[k].tn == 15382) fprintf(stderr, "-mm-[M::%s]\tqn::%u\ttn::%u\terr::%u\n", __func__, rid, oa[k].tn, ol->list[ol->length-1].non_homopolymer_errors);
        } 

        if((ol->length <= ol0) || (ol->list[ol0].is_match != 1)) continue;

        if(m != ol0) {
            t = ol->list[m];
            ol->list[m] = ol->list[ol0];
            ol->list[ol0] = t;
        }
        
        m++;
    }

    ol->length = m;
    
    // fprintf(stderr, "[M::%s]\tnew_n::%lu\told_n::%u\n", __func__, ol->length, in0->length + in1->length);
    assert(ol->length <= (in0->length + in1->length));


    // dbg_overlap_region_cigar(ol->list, ol->length, rs, rref, tu);


    return aux_o;
}


overlap_region* h_ec_lchain_re1(ha_abuf_t *ab, uint32_t rid, UC_Read *qu, UC_Read *tu, uint64_t mz_w, uint64_t mz_k, All_reads *rref, overlap_region_alloc *ol, Candidates_list *cl, bit_extz_t *exz, asg16_v *buf, asg64_v *srt_i, double bw_thres, 
								 int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t is_accurate, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, ma_hit_t_alloc *in0, ma_hit_t_alloc *in1)
{
    // fprintf(stderr, "-mm-[M::%s]\tchain_cutoff::%u\n", __func__, chain_cutoff);
    uint64_t on = 0, k, one = 0, ol0, wl = (asm_opt.is_ont)?(WINDOW_OHC):(WINDOW_HC), m, m0, tid, trev; ma_hit_t *oa = NULL, *p = NULL; overlap_region *aux_o = NULL, *z = NULL, t; Window_Pool w; double err = asm_opt.max_ov_diff_ec;
    char* rs = qu->seq; uint64_t rl = qu->length; 
    int64_t max_skip, max_iter, max_dis, quick_check; double chn_pen_gap, chn_pen_skip; 
	set_lchain_dp_op(is_accurate, mz_k, &max_skip, &max_iter, &max_dis, &chn_pen_gap, &chn_pen_skip, &quick_check);
    init_Window_Pool(&w, rl, wl, (int)(1.0/err));

    srt_i->n = 0;
    oa = in0->buffer; on = in0->length; m0 = 0;
    for (k = 0; k < on; k++) {
        if(oa[k].el) {
            one++; continue;
        }
        m = oa[k].tn; m <<= 1; m |= ((uint64_t)oa[k].rev); m <<= 32; m |= (k<<1); m |= m0;
        kv_push(uint64_t, *srt_i, m);
    }
    oa = in1->buffer; on = in1->length; m0 = 1;
    for (k = 0; k < on; k++) {
        // if(oa[k].el) continue;
        m = oa[k].tn; m <<= 1; m |= ((uint64_t)oa[k].rev); m <<= 32; m |= (k<<1); m |= m0;
        kv_push(uint64_t, *srt_i, m);
    }
    radix_sort_ec64(srt_i->a, srt_i->a + srt_i->n);

    // get the list of anchors
    get_mz1(rs, rl, mz_w, mz_k, 0, !(asm_opt.flag & HA_F_NO_HPC), ab, ha_flt_tab, NULL/**ha_idx**/, asm_opt.mz_sample_dist, k_flag, dbg_ct, NULL, -1, asm_opt.dp_min_len, -1, sp, asm_opt.mz_rewin, 0, NULL, 0);

    h_ec_lchain_re_gen(ab, rid, rs, rl, mz_w, mz_k, ha_idx, rref, ol, cl, bw_thres, apend_be, k_flag, dbg_ct, sp, high_occ, low_occ, gen_off, enable_mcopy, mcopy_rate, mcopy_khit_cut, 
            max_skip, max_iter, max_dis, quick_check, chn_pen_gap, chn_pen_skip, tu, srt_i, scb.a);

    
    ///max size
    on = ol->length + one + srt_i->n + 1; m0 = on - 1;
    if(on > ol->size) {
        REALLOC(ol->list, on);
        memset(ol->list+ol->size, 0, sizeof(overlap_region)*(on-ol->size));
        ol->size = on;   
    }
    aux_o = &(ol->list[on-1]); ol->mapped_overlaps_length = 0; on = ol->length;

    // fprintf(stderr, "-0-[M::%s]\n", __func__);

    gen_hc_r_alin(ol, cl, rref, qu, tu, exz, aux_o, asm_opt.max_ov_diff_ec, w.window_length, rid, E_KHIT, 1, buf, 0); rs = qu->seq; 

    // fprintf(stderr, "-1-[M::%s]\n", __func__);

    ///handle unmatched chain
    for (k = m = ol->length; k < on; k++) {
        clear_fake_cigar(&(ol->list[k].f_cigar));
        clear_window_list_alloc(&(ol->list[k].w_list));
        clear_window_list_alloc(&(ol->list[k].boundary_cigars));

        ol0 = ol->length;

        tid = ol->list[k].y_id; trev = ol->list[k].y_pos_strand;
        if(trev) recover_UC_Read_RC(tu, rref, tid);
        else recover_UC_Read(tu, rref, tid);

        get_pi_ec_chain(ab, rid, rl, tid, tu->seq, tu->length, mz_w, mz_k, ol, cl, bw_thres, apend_be, k_flag, dbg_ct, sp, high_occ, low_occ, gen_off, enable_mcopy, mcopy_rate, mcopy_khit_cut, max_skip, max_iter, max_dis, quick_check, chn_pen_gap, chn_pen_skip);
        assert(ol->length - ol0 <= 1);
        if((ol->length > ol0) && (gen_hc_r_alin_re(&(ol->list[ol0]), cl, rs, rl, tu->seq, tu->length, exz, aux_o, asm_opt.max_ov_diff_ec, w.window_length, rid, E_KHIT, 1, buf))) {
            ol->list[ol0].y_pos_strand = trev;
            // if(oa[k].tn == 15382) fprintf(stderr, "-mm-[M::%s]\tqn::%u\ttn::%u\terr::%u\n", __func__, rid, oa[k].tn, ol->list[ol->length-1].non_homopolymer_errors);
        } 

        if((ol->length <= ol0) || (ol->list[ol0].is_match != 1)) continue;

        if(m != ol0) {
            t = ol->list[m];
            ol->list[m] = ol->list[ol0];
            ol->list[ol0] = t;
        }

        m++;
    }
    ol->length = m;

    for (k = ol->length; k < on; k++) {
        clear_fake_cigar(&(ol->list[k].f_cigar));
        clear_window_list_alloc(&(ol->list[k].w_list));
        clear_window_list_alloc(&(ol->list[k].boundary_cigars));
    }

    // fprintf(stderr, "[M::%s]\tnew::%lu\told::%lu\n", __func__, ol->length, ol0);
    
    oa = in0->buffer; on = in0->length;
    for (k = 0; k < on; k++) {
        if(oa[k].el) {
            z = &(ol->list[ol->length++]);
            z->x_id = rid; z->y_id = oa[k].tn;
            z->x_pos_strand = 0; z->y_pos_strand = oa[k].rev;
            z->x_pos_s = (uint32_t)oa[k].qns;
            z->x_pos_e = oa[k].qe - 1;
            z->y_pos_s = oa[k].ts;
            z->y_pos_e = oa[k].te - 1;

            z->is_match = 1;
            z->align_length = z->overlapLen = z->shared_seed = z->x_pos_e + 1 - z->x_pos_s;
            z->non_homopolymer_errors = z->strong = 0;

            set_exact_exz(exz, z->x_pos_s, z->x_pos_e + 1, z->y_pos_s, z->y_pos_e + 1); push_alnw(z, exz);

            // if(oa[k].tn == 1945) fprintf(stderr, "-em-[M::%s]\tqn::%u\ttn::%u\terr::%u\n", __func__, rid, oa[k].tn, ol->list[ol->length-1].non_homopolymer_errors);
        }
    }
    m = ol->length;

    for (k = 0; k < srt_i->n; k++) {
        ol0 = ol->length;
        if(srt_i->a[k]&1) {
            p = &(in1->buffer[((uint32_t)srt_i->a[k])>>1]);
        } else {
            p = &(in0->buffer[((uint32_t)srt_i->a[k])>>1]);
        }
        
        if(p->rev) recover_UC_Read_RC(tu, rref, p->tn);
        else recover_UC_Read(tu, rref, p->tn);

        get_pi_ec_chain(ab, rid, rl, p->tn, tu->seq, tu->length, mz_w, mz_k, ol, cl, bw_thres, apend_be, k_flag, dbg_ct, sp, high_occ, low_occ, gen_off, enable_mcopy, mcopy_rate, mcopy_khit_cut, max_skip, max_iter, max_dis, quick_check, chn_pen_gap, chn_pen_skip);
        assert(ol->length - ol0 <= 1);
        if((ol->length > ol0) && (gen_hc_r_alin_re(&(ol->list[ol0]), cl, rs, rl, tu->seq, tu->length, exz, aux_o, asm_opt.max_ov_diff_ec, w.window_length, rid, E_KHIT, 1, buf))) {
            ol->list[ol0].y_pos_strand = p->rev;
            // if(oa[k].tn == 15382) fprintf(stderr, "-mm-[M::%s]\tqn::%u\ttn::%u\terr::%u\n", __func__, rid, oa[k].tn, ol->list[ol->length-1].non_homopolymer_errors);
        } 

        if((ol->length <= ol0) || (ol->list[ol0].is_match != 1)) continue;

        if(m != ol0) {
            t = ol->list[m];
            ol->list[m] = ol->list[ol0];
            ol->list[ol0] = t;
        }

        m++;
    }
    
    ol->length = m;

    
    // fprintf(stderr, "[M::%s]\tnew_n::%lu\told_n::%lu\n", __func__, ol->length, m0);
    // fprintf(stderr, "[M::%s]\tnew_n::%lu\told_n::%u\n", __func__, ol->length, in0->length + in1->length);
    // assert(ol->length <= (in0->length + in1->length));


    // dbg_overlap_region_cigar(ol->list, ol->length, rs, rref, tu);

    return aux_o;
}

uint64_t direct_chain_cal(ha_abuf_t *ab, uint64_t qid, char *qs, uint64_t ql, uint64_t tid, char *ts, uint64_t tl, uint64_t trev, uint64_t mz_w, uint64_t mz_k, overlap_region_alloc *olst, Candidates_list *cl, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, double bw_thres, 
								 int apend_be, uint64_t max_cnt, uint64_t min_cnt, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, int64_t max_skip, int64_t max_iter, int64_t max_dis, int64_t quick_check, double chn_pen_gap, double chn_pen_skip,
                                 bit_extz_t *exz, overlap_region *aux_o, double e_rate, int64_t wl, int64_t khit, int64_t move_gap, asg16_v* buf)
{
    uint64_t ol0 = olst->length;
    get_pi_ec_chain(ab, qid, ql, tid, ts, tl, mz_w, mz_k, olst, cl, bw_thres, apend_be, k_flag, dbg_ct, sp, high_occ, low_occ, gen_off, enable_mcopy, mcopy_rate, mcopy_khit_cut, max_skip, max_iter, max_dis, quick_check, chn_pen_gap, chn_pen_skip);
    assert(olst->length - ol0 <= 1);
    if(olst->length > ol0) {
        if(gen_hc_r_alin_re(&(olst->list[ol0]), cl, qs, ql, ts, tl, exz, aux_o, e_rate, wl, qid, E_KHIT, 1, buf)) {
            olst->list[ol0].y_pos_strand = trev;
            return 1;
        } else {
            clear_fake_cigar(&(olst->list[ol0].f_cigar));
            clear_window_list_alloc(&(olst->list[ol0].w_list));
            clear_window_list_alloc(&(olst->list[ol0].boundary_cigars));
            olst->length--;
        }
    }
    return 0;
}


overlap_region* h_ec_lchain_re2(ha_abuf_t *ab, uint32_t rid, UC_Read *qu, UC_Read *tu, uint64_t mz_w, uint64_t mz_k, All_reads *rref, overlap_region_alloc *ol, Candidates_list *cl, bit_extz_t *exz, asg16_v *buf, asg64_v *srt_i, double bw_thres, 
								 int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t is_accurate, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, ma_hit_t_alloc *in0, ma_hit_t_alloc *in1)
{
    // fprintf(stderr, "-0-[M::%s]\tnew_n::%lu\told_n::%u\n", __func__, ol->length, in0->length + in1->length);
    // fprintf(stderr, "-mm-[M::%s]\tchain_cutoff::%u\n", __func__, chain_cutoff);
    uint64_t on = 0, k, l, i, one = 0, ol0, wl = (asm_opt.is_ont)?(WINDOW_OHC):(WINDOW_HC), m, m0, tid, trev, max_cnt = UINT32_MAX, min_cnt = 0; ma_hit_t *oa = NULL, *p = NULL; overlap_region *aux_o = NULL, *z = NULL, t; Window_Pool w; double err = asm_opt.max_ov_diff_ec; tiny_queue_t tq; memset(&tq, 0, sizeof(tiny_queue_t));
    char* rs = qu->seq; uint64_t rl = qu->length; int64_t n, zn, om;
    int64_t max_skip, max_iter, max_dis, quick_check; double chn_pen_gap, chn_pen_skip; 
	set_lchain_dp_op(is_accurate, mz_k, &max_skip, &max_iter, &max_dis, &chn_pen_gap, &chn_pen_skip, &quick_check);
    init_Window_Pool(&w, rl, wl, (int)(1.0/err));
    
    ///cutoff
    if(high_occ) {
        max_cnt = (*high_occ);
        if(max_cnt < 2) max_cnt = 2;
    }
    if(low_occ) {
        min_cnt = (*low_occ);
        if(min_cnt < 2) min_cnt = 2;
    }
    
    ///memory
    ol->length = 0; on = in0->length + in1->length + 1;
    if(on > ol->size) {
        REALLOC(ol->list, on);
        memset(ol->list+ol->size, 0, sizeof(overlap_region)*(on-ol->size));
        ol->size = on;   
    }
    aux_o = &(ol->list[on-1]); ol->mapped_overlaps_length = 0; on = in0->length + in1->length;
    
    ///overlap idx
    srt_i->n = 0;
    oa = in0->buffer; on = in0->length; m0 = 0;
    for (k = 0; k < on; k++) {
        if(oa[k].el) {
            one++; continue;
        }
        m = oa[k].tn; m <<= 1; m |= ((uint64_t)oa[k].rev); m <<= 32; m |= (k<<1); m |= m0;
        kv_push(uint64_t, *srt_i, m);
    }
    oa = in1->buffer; on = in1->length; m0 = 1;
    for (k = 0; k < on; k++) {
        // if(oa[k].el) continue;
        m = oa[k].tn; m <<= 1; m |= ((uint64_t)oa[k].rev); m <<= 32; m |= (k<<1); m |= m0;
        kv_push(uint64_t, *srt_i, m);
    }
    radix_sort_ec64(srt_i->a, srt_i->a + srt_i->n);

    // get the list of anchors
    get_mz1(rs, rl, mz_w, mz_k, 0, !(asm_opt.flag & HA_F_NO_HPC), ab, ha_flt_tab, NULL/**ha_idx**/, asm_opt.mz_sample_dist, k_flag, dbg_ct, NULL, -1, asm_opt.dp_min_len, -1, sp, asm_opt.mz_rewin, 0, NULL, 0);

    h_ec_lchain_re_gen_srt(ab, ha_idx, ol, cl);

    k = 1; l = 0; i = 0; n = zn = 0;
    while(h_ec_lchain_re_gen_qry(ab, &k, &l, &i, srt_i->a, srt_i->n, &tid, &trev)) {
        
        if(trev) recover_UC_Read_RC(tu, rref, tid);
        else recover_UC_Read(tu, rref, tid);

        ol0 = ol->length; om = 0;
        if(h_ec_lchain_re_chn(ab, l, k, rid, rs, rl, tid, tu->seq, tu->length, trev, mz_w, mz_k, ol, cl, bw_thres, apend_be, max_cnt, min_cnt, gen_off, enable_mcopy, mcopy_rate, mcopy_khit_cut, max_skip, max_iter, max_dis, quick_check, chn_pen_gap, chn_pen_skip, &tq, scb.a, &n, &zn)) {
            assert(ol->length - ol0 == 1);
            ol->list[ol0].y_pos_strand = 0;
            if(gen_hc_r_alin_re(&(ol->list[ol0]), cl, rs, rl, tu->seq, tu->length, exz, aux_o, asm_opt.max_ov_diff_ec, w.window_length, rid, E_KHIT, 1, buf)) {
                ol->list[ol0].y_pos_strand = trev; om = 1;
            } else {///unmatch
                clear_fake_cigar(&(ol->list[ol0].f_cigar));
                clear_window_list_alloc(&(ol->list[ol0].w_list));
                clear_window_list_alloc(&(ol->list[ol0].boundary_cigars));
                ol->length--;
            }
        }

        if(om) {
            srt_i->a[i] >>= 32; srt_i->a[i] <<= 32; srt_i->a[i] |= ((uint64_t)((uint32_t)-1));
        }

        l = k; k++; 
    }

    for (k = m = 0; k < srt_i->n; k++) {
        if(((uint32_t)srt_i->a[k]) == ((uint32_t)-1)) continue;
        srt_i->a[m++] = srt_i->a[k];
    }
    // fprintf(stderr, "[M::%s]\ttot::%lu\tremain::%lu\n", __func__, (uint64_t)srt_i->n, m);
    srt_i->n = m;
    


    oa = in0->buffer; on = in0->length;
    for (k = 0; k < on; k++) {
        if(oa[k].el) {
            z = &(ol->list[ol->length++]);
            z->x_id = rid; z->y_id = oa[k].tn;
            z->x_pos_strand = 0; z->y_pos_strand = oa[k].rev;
            z->x_pos_s = (uint32_t)oa[k].qns;
            z->x_pos_e = oa[k].qe - 1;
            z->y_pos_s = oa[k].ts;
            z->y_pos_e = oa[k].te - 1;

            z->is_match = 1;
            z->align_length = z->overlapLen = z->shared_seed = z->x_pos_e + 1 - z->x_pos_s;
            z->non_homopolymer_errors = z->strong = 0;

            set_exact_exz(exz, z->x_pos_s, z->x_pos_e + 1, z->y_pos_s, z->y_pos_e + 1); push_alnw(z, exz);

            // if(oa[k].tn == 1945) fprintf(stderr, "-em-[M::%s]\tqn::%u\ttn::%u\terr::%u\n", __func__, rid, oa[k].tn, ol->list[ol->length-1].non_homopolymer_errors);
        }
    }
    m = ol->length;


    clear_Candidates_list(cl);
    for (k = 0; k < srt_i->n; k++) {
        ol0 = ol->length;
        if(srt_i->a[k]&1) {
            p = &(in1->buffer[((uint32_t)srt_i->a[k])>>1]);
        } else {
            p = &(in0->buffer[((uint32_t)srt_i->a[k])>>1]);
        }
        
        if(p->rev) recover_UC_Read_RC(tu, rref, p->tn);
        else recover_UC_Read(tu, rref, p->tn);

        get_pi_ec_chain(ab, rid, rl, p->tn, tu->seq, tu->length, mz_w, mz_k, ol, cl, bw_thres, apend_be, k_flag, dbg_ct, sp, high_occ, low_occ, gen_off, enable_mcopy, mcopy_rate, mcopy_khit_cut, max_skip, max_iter, max_dis, quick_check, chn_pen_gap, chn_pen_skip);
        assert(ol->length - ol0 <= 1);
        if((ol->length > ol0) && (gen_hc_r_alin_re(&(ol->list[ol0]), cl, rs, rl, tu->seq, tu->length, exz, aux_o, asm_opt.max_ov_diff_ec, w.window_length, rid, E_KHIT, 1, buf))) {
            ol->list[ol0].y_pos_strand = p->rev;
            // if(oa[k].tn == 15382) fprintf(stderr, "-mm-[M::%s]\tqn::%u\ttn::%u\terr::%u\n", __func__, rid, oa[k].tn, ol->list[ol->length-1].non_homopolymer_errors);
        } 

        if((ol->length <= ol0) || (ol->list[ol0].is_match != 1)) continue;

        if(m != ol0) {
            t = ol->list[m];
            ol->list[m] = ol->list[ol0];
            ol->list[ol0] = t;
        }

        m++;
    }
    
    ol->length = m;
    
    // fprintf(stderr, "[M::%s]\tnew_n::%lu\told_n::%lu\n", __func__, ol->length, m0);
    // fprintf(stderr, "-1-[M::%s]\tnew_n::%lu\told_n::%u\n", __func__, ol->length, in0->length + in1->length);
    assert(ol->length <= (in0->length + in1->length));


    // dbg_overlap_region_cigar(ol->list, ol->length, rs, rref, tu);

    return aux_o;
}


overlap_region* h_ec_lchain_fast(ha_abuf_t *ab, uint32_t rid, UC_Read *qu, UC_Read *tu, uint64_t mz_w, uint64_t mz_k, All_reads *rref, overlap_region_alloc *ol, Candidates_list *cl, bit_extz_t *exz, asg16_v *buf, asg64_v *srt_i, double bw_thres, 
								 int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t is_accurate, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, ma_hit_t_alloc *in0, ma_hit_t_alloc *in1, double sh)
{
    // fprintf(stderr, "-0-[M::%s]\tnew_n::%lu\told_n::%u\n", __func__, ol->length, in0->length + in1->length);
    // fprintf(stderr, "-mm-[M::%s]\tchain_cutoff::%u\n", __func__, chain_cutoff);
    uint64_t on = 0, k, l, i, one = 0, ol0, wl = (asm_opt.is_ont)?(WINDOW_OHC):(WINDOW_HC), m, m0, tid, trev, max_cnt = UINT32_MAX, min_cnt = 0, is_match; ma_hit_t *oa = NULL, *p = NULL; overlap_region *aux_o = NULL, *z = NULL; Window_Pool w; double err = asm_opt.max_ov_diff_ec; tiny_queue_t tq; memset(&tq, 0, sizeof(tiny_queue_t));
    char* rs = qu->seq; uint64_t rl = qu->length; int64_t n, zn, om; uint64_t aq[2], at[2], bq[2], bt[2], ovlp, os, oe;
    int64_t max_skip, max_iter, max_dis, quick_check; double chn_pen_gap, chn_pen_skip; 
	set_lchain_dp_op(is_accurate, mz_k, &max_skip, &max_iter, &max_dis, &chn_pen_gap, &chn_pen_skip, &quick_check);
    init_Window_Pool(&w, rl, wl, (int)(1.0/err));
    
    ///cutoff
    if(high_occ) {
        max_cnt = (*high_occ);
        if(max_cnt < 2) max_cnt = 2;
    }
    if(low_occ) {
        min_cnt = (*low_occ);
        if(min_cnt < 2) min_cnt = 2;
    }
    
    ///memory
    ol->length = 0; on = in0->length + in1->length + 1;
    if(on > ol->size) {
        REALLOC(ol->list, on);
        memset(ol->list+ol->size, 0, sizeof(overlap_region)*(on-ol->size));
        ol->size = on;   
    }
    aux_o = &(ol->list[on-1]); ol->mapped_overlaps_length = 0; on = in0->length + in1->length;
    
    ///overlap idx
    srt_i->n = 0;
    oa = in0->buffer; on = in0->length; m0 = 0;
    for (k = 0; k < on; k++) {
        if(oa[k].el) {
            one++; continue;
        }
        m = oa[k].tn; m <<= 1; m |= ((uint64_t)oa[k].rev); m <<= 32; m |= (k<<1); m |= m0;
        kv_push(uint64_t, *srt_i, m);
    }
    oa = in1->buffer; on = in1->length; m0 = 1;
    for (k = 0; k < on; k++) {
        // if(oa[k].el) continue;
        m = oa[k].tn; m <<= 1; m |= ((uint64_t)oa[k].rev); m <<= 32; m |= (k<<1); m |= m0;
        kv_push(uint64_t, *srt_i, m);
    }
    radix_sort_ec64(srt_i->a, srt_i->a + srt_i->n);

    // get the list of anchors
    get_mz1(rs, rl, mz_w, mz_k, 0, !(asm_opt.flag & HA_F_NO_HPC), ab, ha_flt_tab, NULL/**ha_idx**/, asm_opt.mz_sample_dist, k_flag, dbg_ct, NULL, -1, asm_opt.dp_min_len, -1, sp, asm_opt.mz_rewin, 0, NULL, 0);

    h_ec_lchain_re_gen_srt(ab, ha_idx, ol, cl);

    k = 1; l = 0; i = 0; n = zn = 0;
    while(h_ec_lchain_re_gen_qry(ab, &k, &l, &i, srt_i->a, srt_i->n, &tid, &trev)) {
        
        if(trev) recover_UC_Read_RC(tu, rref, tid);
        else recover_UC_Read(tu, rref, tid);

        ol0 = ol->length; om = 0;
        if(h_ec_lchain_re_chn(ab, l, k, rid, rs, rl, tid, tu->seq, tu->length, trev, mz_w, mz_k, ol, cl, bw_thres, apend_be, max_cnt, min_cnt, gen_off, enable_mcopy, mcopy_rate, mcopy_khit_cut, max_skip, max_iter, max_dis, quick_check, chn_pen_gap, chn_pen_skip, &tq, scb.a, &n, &zn)) {
            assert(ol->length - ol0 == 1);
            ol->list[ol0].y_pos_strand = trev; om = 1; ol->list[ol0].shared_seed = 0; ol->list[ol0].is_match = 0;
            
            assert(((srt_i->a[i]>>33) == ol->list[ol0].y_id) && (((srt_i->a[i]>>32)&1) == ol->list[ol0].y_pos_strand));
            if(srt_i->a[i]&1) {
                p = &(in1->buffer[((uint32_t)srt_i->a[i])>>1]); is_match = 2;
            } else {
                p = &(in0->buffer[((uint32_t)srt_i->a[i])>>1]); is_match = 1;
            }
            ol->list[ol0].strong = p->ml; ol->list[ol0].without_large_indel = p->no_l_indel;

            aq[0] = (uint32_t)p->qns; aq[1] = p->qe; 
            at[0] = p->ts; at[1] = p->te;

            bq[0] = ol->list[ol0].x_pos_s; bq[1] = ol->list[ol0].x_pos_e + 1;
            bt[0] = ol->list[ol0].y_pos_s; bt[1] = ol->list[ol0].y_pos_e + 1;
            
            os = MAX(aq[0], bq[0]); oe = MIN(aq[1], bq[1]);
            ovlp = ((oe>os)? (oe-os):0);
            if(!((ovlp) && (ovlp >= ((aq[1] - aq[0])*sh)) && ((ovlp >= ((bq[1] - bq[0])*sh))))) om = 0;

            os = MAX(at[0], bt[0]); oe = MIN(at[1], bt[1]);
            ovlp = ((oe>os)? (oe-os):0);
            if(!((ovlp) && (ovlp >= ((at[1] - at[0])*sh)) && ((ovlp >= ((bt[1] - bt[0])*sh))))) om = 0;

            if(om) {
                if(exact_ec_check(rs, rl, tu->seq, tu->length, bq[0], bq[1], bt[0], bt[1])) {
                    if(is_match == 2) {
                        ol->list[ol0].strong = 0; ol->list[ol0].without_large_indel = 1;
                    }
                    is_match = 1; ol->list[ol0].shared_seed = 1;
                }
                ol->list[ol0].is_match = is_match;
            } else {///unmatch
                clear_fake_cigar(&(ol->list[ol0].f_cigar));
                clear_window_list_alloc(&(ol->list[ol0].w_list));
                clear_window_list_alloc(&(ol->list[ol0].boundary_cigars));
                ol->length--;
            }
        }

        if(om) {
            srt_i->a[i] >>= 32; srt_i->a[i] <<= 32; srt_i->a[i] |= ((uint64_t)((uint32_t)-1));
        }

        l = k; k++; 
    }

    for (k = m = 0; k < srt_i->n; k++) {
        if(((uint32_t)srt_i->a[k]) == ((uint32_t)-1)) continue;
        srt_i->a[m++] = srt_i->a[k];
    }
    // fprintf(stderr, "[M::%s]\ttot::%lu\tremain::%lu\n", __func__, (uint64_t)srt_i->n, m);
    srt_i->n = m;
    


    oa = in0->buffer; on = in0->length;
    for (k = 0; k < on; k++) {
        if(oa[k].el) {
            z = &(ol->list[ol->length++]);
            z->x_id = rid; z->y_id = oa[k].tn;
            z->x_pos_strand = 0; z->y_pos_strand = oa[k].rev;
            z->x_pos_s = (uint32_t)oa[k].qns;
            z->x_pos_e = oa[k].qe - 1;
            z->y_pos_s = oa[k].ts;
            z->y_pos_e = oa[k].te - 1;

            z->align_length = z->overlapLen = z->x_pos_e + 1 - z->x_pos_s;
            z->non_homopolymer_errors = 0;

            z->is_match = 1; z->shared_seed = 1;
            z->strong = oa[k].ml; z->without_large_indel = oa[k].no_l_indel;

            set_exact_exz(exz, z->x_pos_s, z->x_pos_e + 1, z->y_pos_s, z->y_pos_e + 1); push_alnw(z, exz);

            // if(oa[k].tn == 1945) fprintf(stderr, "-em-[M::%s]\tqn::%u\ttn::%u\terr::%u\n", __func__, rid, oa[k].tn, ol->list[ol->length-1].non_homopolymer_errors);
        }
    }
    /**
    m = ol->length;


    clear_Candidates_list(cl);
    for (k = 0; k < srt_i->n; k++) {
        ol0 = ol->length;
        if(srt_i->a[k]&1) {
            p = &(in1->buffer[((uint32_t)srt_i->a[k])>>1]);
        } else {
            p = &(in0->buffer[((uint32_t)srt_i->a[k])>>1]);
        }
        
        if(p->rev) recover_UC_Read_RC(tu, rref, p->tn);
        else recover_UC_Read(tu, rref, p->tn);

        get_pi_ec_chain(ab, rid, rl, p->tn, tu->seq, tu->length, mz_w, mz_k, ol, cl, bw_thres, apend_be, k_flag, dbg_ct, sp, high_occ, low_occ, gen_off, enable_mcopy, mcopy_rate, mcopy_khit_cut, max_skip, max_iter, max_dis, quick_check, chn_pen_gap, chn_pen_skip);
        assert(ol->length - ol0 <= 1);
        if((ol->length > ol0) && (gen_hc_r_alin_re(&(ol->list[ol0]), cl, rs, rl, tu->seq, tu->length, exz, aux_o, asm_opt.max_ov_diff_ec, w.window_length, rid, E_KHIT, 1, buf))) {
            ol->list[ol0].y_pos_strand = p->rev;
            // if(oa[k].tn == 15382) fprintf(stderr, "-mm-[M::%s]\tqn::%u\ttn::%u\terr::%u\n", __func__, rid, oa[k].tn, ol->list[ol->length-1].non_homopolymer_errors);
        } 

        if((ol->length <= ol0) || (ol->list[ol0].is_match != 1)) continue;

        if(m != ol0) {
            t = ol->list[m];
            ol->list[m] = ol->list[ol0];
            ol->list[ol0] = t;
        }

        m++;
    }
    
    ol->length = m;
    **/
    
    // fprintf(stderr, "[M::%s]\tnew_n::%lu\told_n::%lu\n", __func__, ol->length, m0);
    // fprintf(stderr, "-1-[M::%s]\tnew_n::%lu\told_n::%u\n", __func__, ol->length, in0->length + in1->length);
    assert(ol->length <= (in0->length + in1->length));


    // dbg_overlap_region_cigar(ol->list, ol->length, rs, rref, tu);

    return aux_o;
}

void h_ec_lchain_fast_new(ha_abuf_t *ab, uint32_t rid, UC_Read *qu, UC_Read *tu, All_reads *rref, overlap_region_alloc *ol, Candidates_list *cl, bit_extz_t *exz, asg16_v *buf, asg64_v *srt_i, ma_hit_t_alloc *in0, ma_hit_t_alloc *in1, double sh)
{
    // fprintf(stderr, "-0-[M::%s]\tnew_n::%lu\told_n::%u\n", __func__, ol->length, in0->length + in1->length);
    // fprintf(stderr, "-mm-[M::%s]\tchain_cutoff::%u\n", __func__, chain_cutoff);
    uint64_t on = 0, k, i, l, m, m0, tid, trev, is_match, is_usrt = 0; ma_hit_t *oa = NULL, *p = NULL; overlap_region *z = NULL, t; 
    char* rs = qu->seq; uint64_t rl = qu->length; int64_t om; uint64_t aq[2], at[2], bq[2], bt[2], ovlp, os, oe;
    

    
    ///overlap idx
    srt_i->n = 0;
    oa = in0->buffer; on = in0->length; m0 = 0;
    for (k = 0; k < on; k++) {
        // if(oa[k].el) continue;
        m = oa[k].tn; m <<= 1; m |= ((uint64_t)oa[k].rev); m <<= 32; m |= (k<<1); m |= m0;
        kv_push(uint64_t, *srt_i, m);
    }
    oa = in1->buffer; on = in1->length; m0 = 1;
    for (k = 0; k < on; k++) {
        // if(oa[k].el) continue;
        m = oa[k].tn; m <<= 1; m |= ((uint64_t)oa[k].rev); m <<= 32; m |= (k<<1); m |= m0;
        kv_push(uint64_t, *srt_i, m);
    }
    radix_sort_ec64(srt_i->a, srt_i->a + srt_i->n);

    k = 0; i = 0; 
    for (k = m = 0; k < ol->length; k++) {
        z = &(ol->list[k]); tid = z->y_id; trev = z->y_pos_strand;
        z->non_homopolymer_errors = 0;
        for (; (i < srt_i->n) && ((srt_i->a[i]>>32) < ((tid<<1)|trev)); i++);
        if((i < srt_i->n) && ((srt_i->a[i]>>32) == ((tid<<1)|trev))) {
            om = 1; z->shared_seed = 0; z->is_match = 0;
            if(srt_i->a[i]&1) {
                p = &(in1->buffer[((uint32_t)srt_i->a[i])>>1]); is_match = 2;
            } else {
                p = &(in0->buffer[((uint32_t)srt_i->a[i])>>1]); is_match = 1;
            }
            z->strong = p->ml; z->without_large_indel = p->no_l_indel;
            
            aq[0] = (uint32_t)p->qns; aq[1] = p->qe; 
            at[0] = p->ts; at[1] = p->te;

            bq[0] = z->x_pos_s; bq[1] = z->x_pos_e + 1;
            bt[0] = z->y_pos_s; bt[1] = z->y_pos_e + 1;
            
            os = MAX(aq[0], bq[0]); oe = MIN(aq[1], bq[1]);
            ovlp = ((oe>os)? (oe-os):0);
            if(!((ovlp) && (ovlp >= ((aq[1] - aq[0])*sh)) && ((ovlp >= ((bq[1] - bq[0])*sh))))) om = 0;
            z->non_homopolymer_errors += ((aq[1] - aq[0]) - ovlp);

            os = MAX(at[0], bt[0]); oe = MIN(at[1], bt[1]);
            ovlp = ((oe>os)? (oe-os):0);
            if(!((ovlp) && (ovlp >= ((at[1] - at[0])*sh)) && ((ovlp >= ((bt[1] - bt[0])*sh))))) om = 0;
            z->non_homopolymer_errors += ((at[1] - at[0]) - ovlp);

            if(om) {
                if(is_match == 1 && p->el == 1) p->el = 0;
                resize_UC_Read(tu, bt[1] - bt[0]);
                recover_UC_Read_sub_region(tu->seq, bt[0], bt[1] - bt[0], trev, rref, tid);
                if(exact_ec_check(rs, rl, tu->seq, bt[1] - bt[0], bq[0], bq[1], 0, bt[1] - bt[0])) {
                    if(is_match == 2) {
                        z->strong = 0; z->without_large_indel = 1;
                    }
                    is_match = 1; z->shared_seed = 1;
                }
                z->is_match = is_match;
            }
        } else {
            om = 0;
            bq[0] = z->x_pos_s; bq[1] = z->x_pos_e + 1;
            bt[0] = z->y_pos_s; bt[1] = z->y_pos_e + 1;
            resize_UC_Read(tu, bt[1] - bt[0]);
            recover_UC_Read_sub_region(tu->seq, bt[0], bt[1] - bt[0], trev, rref, tid);
            if(exact_ec_check(rs, rl, tu->seq, bt[1] - bt[0], bq[0], bq[1], 0, bt[1] - bt[0])) {
                z->strong = 0; z->without_large_indel = 1;
                z->shared_seed = 1; z->is_match = 1; om = 1;
            }
        }

        if(om) {
            if(m != k) {
                t = ol->list[m];
                ol->list[m] = ol->list[k];
                ol->list[k] = t;
            }
            m++;
        }
    }
    ol->length = m;


    oa = in0->buffer; on = in0->length;
    for (k = 0; k < on; k++) {
        if(oa[k].el) {
            // z = &(ol->list[ol->length++]);
            kv_pushp_ol(overlap_region, (*ol), &z);
            clear_fake_cigar(&(z->f_cigar));
            clear_window_list_alloc(&(z->w_list));
            clear_window_list_alloc(&(z->boundary_cigars));
            
            z->x_id = rid; z->y_id = oa[k].tn;
            z->x_pos_strand = 0; z->y_pos_strand = oa[k].rev;
            z->x_pos_s = (uint32_t)oa[k].qns;
            z->x_pos_e = oa[k].qe - 1;
            z->y_pos_s = oa[k].ts;
            z->y_pos_e = oa[k].te - 1;

            z->align_length = z->overlapLen = z->x_pos_e + 1 - z->x_pos_s;
            z->non_homopolymer_errors = 0;

            z->is_match = 1; z->shared_seed = 1;
            z->strong = oa[k].ml; z->without_large_indel = oa[k].no_l_indel;

            set_exact_exz(exz, z->x_pos_s, z->x_pos_e + 1, z->y_pos_s, z->y_pos_e + 1); push_alnw(z, exz);

            is_usrt = 1;
            // if(oa[k].tn == 1945) fprintf(stderr, "-em-[M::%s]\tqn::%u\ttn::%u\terr::%u\n", __func__, rid, oa[k].tn, ol->list[ol->length-1].non_homopolymer_errors);
        }
    }

    if(is_usrt) overlap_region_sort_y_id(ol->list, ol->length);

    if(ol->length > 1) {///for duplicated chains
        uint64_t mm_k, s; int64_t mm_sc, sc;
        for (k = 1, l = m = 0; k <= ol->length; k++) {
            if(k == ol->length || ol->list[k].y_id != ol->list[l].y_id) {
                mm_k = l;
                if(k - l > 1) {
                    for (s = l, mm_sc = INT32_MIN, mm_k = ((uint64_t)-1); s < k; s++) {
                        z = &(ol->list[s]);
                        sc = z->non_homopolymer_errors; sc = - sc;
                        if((sc > mm_sc) || ((sc == mm_sc) && ((ol->list[mm_k].x_pos_e+1-ol->list[mm_k].x_pos_s) < (z->x_pos_e+1-z->x_pos_s)))) {
                            mm_sc = sc; mm_k = s;
                        }
                    }
                }
                if(mm_k != ((uint64_t)-1)) {
                    if(mm_k != m) {
                        t = ol->list[mm_k];
                        ol->list[mm_k] = ol->list[m];
                        ol->list[m] = t;
                    }
                    m++;
                }
                l = k;
            }
        }
        ol->length = m;
    }
}

overlap_region* h_ec_lchain_re3(ha_abuf_t *ab, uint32_t rid, UC_Read *qu, UC_Read *tu, uint64_t mz_w, uint64_t mz_k, All_reads *rref, overlap_region_alloc *ol, Candidates_list *cl, bit_extz_t *exz, asg16_v *buf, asg64_v *srt_i, double bw_thres, 
                                 int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t is_accurate, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, ma_hit_t_alloc *in0, ma_hit_t_alloc *in1)
{
    // fprintf(stderr, "-mm-[M::%s]\tchain_cutoff::%u\n", __func__, chain_cutoff);
    uint64_t on = 0, k, one = 0, ol0, wl = (asm_opt.is_ont)?(WINDOW_OHC):(WINDOW_HC), m, m0, tid, trev; ma_hit_t *oa = NULL, *p = NULL; overlap_region *aux_o = NULL, *z = NULL, t; Window_Pool w; double err = asm_opt.max_ov_diff_ec;
    char* rs = qu->seq; uint64_t rl = qu->length; 
    int64_t max_skip, max_iter, max_dis, quick_check; double chn_pen_gap, chn_pen_skip; 
    set_lchain_dp_op(is_accurate, mz_k, &max_skip, &max_iter, &max_dis, &chn_pen_gap, &chn_pen_skip, &quick_check);
    init_Window_Pool(&w, rl, wl, (int)(1.0/err));

    srt_i->n = 0;
    oa = in0->buffer; on = in0->length; m0 = 0;
    for (k = 0; k < on; k++) {
        m = oa[k].tn; m <<= 1; m |= ((uint64_t)oa[k].rev); m <<= 32; 
        if(oa[k].el) {
            one++; m |= ((uint32_t)-1);
        } else {
           m |= (k<<1); m |= m0;
        }
        kv_push(uint64_t, *srt_i, m);
    }
    oa = in1->buffer; on = in1->length; m0 = 1;
    for (k = 0; k < on; k++) {
        // if(oa[k].el) continue;
        m = oa[k].tn; m <<= 1; m |= ((uint64_t)oa[k].rev); m <<= 32; m |= (k<<1); m |= m0;
        kv_push(uint64_t, *srt_i, m);
    }
    radix_sort_ec64(srt_i->a, srt_i->a + srt_i->n);

    // get the list of anchors
    get_mz1(rs, rl, mz_w, mz_k, 0, !(asm_opt.flag & HA_F_NO_HPC), ab, ha_flt_tab, NULL/**ha_idx**/, asm_opt.mz_sample_dist, k_flag, dbg_ct, NULL, -1, asm_opt.dp_min_len, -1, sp, asm_opt.mz_rewin, 0, NULL, 0);

    h_ec_lchain_re_gen3(ab, rid, rs, rl, mz_w, mz_k, ha_idx, rref, ol, cl, bw_thres, apend_be, k_flag, dbg_ct, sp, high_occ, low_occ, gen_off, enable_mcopy, mcopy_rate, mcopy_khit_cut, 
            max_skip, max_iter, max_dis, quick_check, chn_pen_gap, chn_pen_skip, tu, srt_i, scb.a);

    
    ///max size
    on = ol->length + one + srt_i->n + 1; m0 = on - 1;
    if(on > ol->size) {
        REALLOC(ol->list, on);
        memset(ol->list+ol->size, 0, sizeof(overlap_region)*(on-ol->size));
        ol->size = on;   
    }
    aux_o = &(ol->list[on-1]); ol->mapped_overlaps_length = 0; on = ol->length;

    // fprintf(stderr, "-0-[M::%s]\n", __func__);

    gen_hc_r_alin(ol, cl, rref, qu, tu, exz, aux_o, asm_opt.max_ov_diff_ec, w.window_length, rid, E_KHIT, 1, buf, 0); rs = qu->seq; 

    // fprintf(stderr, "-1-[M::%s]\n", __func__);

    ///handle unmatched chain
    for (k = m = ol->length; k < on; k++) {
        clear_fake_cigar(&(ol->list[k].f_cigar));
        clear_window_list_alloc(&(ol->list[k].w_list));
        clear_window_list_alloc(&(ol->list[k].boundary_cigars));

        ol0 = ol->length;

        tid = ol->list[k].y_id; trev = ol->list[k].y_pos_strand;
        if(trev) recover_UC_Read_RC(tu, rref, tid);
        else recover_UC_Read(tu, rref, tid);

        get_pi_ec_chain(ab, rid, rl, tid, tu->seq, tu->length, mz_w, mz_k, ol, cl, bw_thres, apend_be, k_flag, dbg_ct, sp, high_occ, low_occ, gen_off, enable_mcopy, mcopy_rate, mcopy_khit_cut, max_skip, max_iter, max_dis, quick_check, chn_pen_gap, chn_pen_skip);
        assert(ol->length - ol0 <= 1);
        if((ol->length > ol0) && (gen_hc_r_alin_re(&(ol->list[ol0]), cl, rs, rl, tu->seq, tu->length, exz, aux_o, asm_opt.max_ov_diff_ec, w.window_length, rid, E_KHIT, 1, buf))) {
            ol->list[ol0].y_pos_strand = trev;
            // if(oa[k].tn == 15382) fprintf(stderr, "-mm-[M::%s]\tqn::%u\ttn::%u\terr::%u\n", __func__, rid, oa[k].tn, ol->list[ol->length-1].non_homopolymer_errors);
        } 

        if((ol->length <= ol0) || (ol->list[ol0].is_match != 1)) continue;

        if(m != ol0) {
            t = ol->list[m];
            ol->list[m] = ol->list[ol0];
            ol->list[ol0] = t;
        }

        m++;
    }
    ol->length = m;

    for (k = ol->length; k < on; k++) {
        clear_fake_cigar(&(ol->list[k].f_cigar));
        clear_window_list_alloc(&(ol->list[k].w_list));
        clear_window_list_alloc(&(ol->list[k].boundary_cigars));
    }

    // fprintf(stderr, "[M::%s]\tnew::%lu\told::%lu\n", __func__, ol->length, ol0);
    
    oa = in0->buffer; on = in0->length;
    for (k = 0; k < on; k++) {
        if(oa[k].el) {
            z = &(ol->list[ol->length++]);
            z->x_id = rid; z->y_id = oa[k].tn;
            z->x_pos_strand = 0; z->y_pos_strand = oa[k].rev;
            z->x_pos_s = (uint32_t)oa[k].qns;
            z->x_pos_e = oa[k].qe - 1;
            z->y_pos_s = oa[k].ts;
            z->y_pos_e = oa[k].te - 1;

            z->is_match = 1;
            z->align_length = z->overlapLen = z->shared_seed = z->x_pos_e + 1 - z->x_pos_s;
            z->non_homopolymer_errors = z->strong = 0;

            set_exact_exz(exz, z->x_pos_s, z->x_pos_e + 1, z->y_pos_s, z->y_pos_e + 1); push_alnw(z, exz);

            // if(oa[k].tn == 1945) fprintf(stderr, "-em-[M::%s]\tqn::%u\ttn::%u\terr::%u\n", __func__, rid, oa[k].tn, ol->list[ol->length-1].non_homopolymer_errors);
        }
    }
    m = ol->length;

    for (k = 0; k < srt_i->n; k++) {
        ol0 = ol->length;
        if(srt_i->a[k]&1) {
            p = &(in1->buffer[((uint32_t)srt_i->a[k])>>1]);
        } else {
            p = &(in0->buffer[((uint32_t)srt_i->a[k])>>1]);
        }
        
        if(p->rev) recover_UC_Read_RC(tu, rref, p->tn);
        else recover_UC_Read(tu, rref, p->tn);

        get_pi_ec_chain(ab, rid, rl, p->tn, tu->seq, tu->length, mz_w, mz_k, ol, cl, bw_thres, apend_be, k_flag, dbg_ct, sp, high_occ, low_occ, gen_off, enable_mcopy, mcopy_rate, mcopy_khit_cut, max_skip, max_iter, max_dis, quick_check, chn_pen_gap, chn_pen_skip);
        assert(ol->length - ol0 <= 1);
        if((ol->length > ol0) && (gen_hc_r_alin_re(&(ol->list[ol0]), cl, rs, rl, tu->seq, tu->length, exz, aux_o, asm_opt.max_ov_diff_ec, w.window_length, rid, E_KHIT, 1, buf))) {
            ol->list[ol0].y_pos_strand = p->rev;
            // if(oa[k].tn == 15382) fprintf(stderr, "-mm-[M::%s]\tqn::%u\ttn::%u\terr::%u\n", __func__, rid, oa[k].tn, ol->list[ol->length-1].non_homopolymer_errors);
        } 

        if((ol->length <= ol0) || (ol->list[ol0].is_match != 1)) continue;

        if(m != ol0) {
            t = ol->list[m];
            ol->list[m] = ol->list[ol0];
            ol->list[ol0] = t;
        }

        m++;
    }
    
    ol->length = m;

    
    // fprintf(stderr, "[M::%s]\tnew_n::%lu\told_n::%lu\n", __func__, ol->length, m0);
    // fprintf(stderr, "[M::%s]\tnew_n::%lu\told_n::%u\n", __func__, ol->length, in0->length + in1->length);
    // assert(ol->length <= (in0->length + in1->length));


    // dbg_overlap_region_cigar(ol->list, ol->length, rs, rref, tu);

    return aux_o;
}


void gen_ori_seq0(char *tstr, uint64_t tl, UC_Read *qu, asg16_v *sc, uint64_t rid)
{
    uint64_t ck, qk, tk, k, wq[2], wt[2]; uint32_t len; uint16_t c, bq, bt; char *qstr = NULL;
    
    ck = qk = tk = 0;
    while (ck < sc->n) {
        wq[0] = qk; wt[0] = tk;
        ck = pop_trace_bp_f(sc, ck, &c, &bq, &bt, &len);
        if(c != 2) qk += len;
        if(c != 3) tk += len;
        wq[1] = qk; wt[1] = tk;
    }
    if(!(tk == tl)) {
        fprintf(stderr, "[M::%s] rid::%lu, tk::%lu, tl::%lu\n", __func__, rid, tk, tl);
    }
    assert(tk == tl);

    resize_UC_Read(qu, qk); qstr = qu->seq; qu->length = qk;
    ck = qk = tk = 0;
    while (ck < sc->n) {
        wq[0] = qk; wt[0] = tk;
        ck = pop_trace_bp_f(sc, ck, &c, &bq, &bt, &len);
        if(c != 2) qk += len;
        if(c != 3) tk += len;
        wq[1] = qk; wt[1] = tk;

        if(c == 0) {
            memcpy(qstr + wq[0], tstr + wt[0], (wq[1]-wq[0])*sizeof((*qstr)));
        } else if(c == 1 || c == 3) {
            for (k = wq[0]; k < wq[1]; k++) qstr[k] = s_H[bq]; 
        }
        // fprintf(stderr, "%u%c(%c)(x::[%lu,%ld))(y::[%lu,%ld))\n", len, cm[c], ((c==1)||(c==2))?(cc[bt]):('*'), wx[0], wx[1], wy[0], wy[1]); // s_H
    }
}

void gen_cc_fly(asg16_v *sc, char *qstr, uint64_t ql, char *tstr, uint64_t tl, bit_extz_t *exz, double e_rate, uint64_t maxn, uint64_t maxe)
{
    // fprintf(stderr, "[M::%s] ql::%lu, tl::%lu\n", __func__, ql, tl);
    if(ql == 0 && tl == 0) return;
    uint64_t k, ck, qk, tk, wq[2], wt[2], maxl, minl, diff, f, diff0 = 0;
    if(ql > 0 && tl == 0) {
        for (k = 0; k < ql; k++) {
            push_trace_bp_f(sc, 3, seq_nt6_table[(uint32_t)(qstr[k])], (uint16_t)-1, 1, 1); 
        }
        return;
    }
    if(ql == 0 && tl > 0) {
        for (k = 0; k < tl; k++) {
            push_trace_bp_f(sc, 2, (uint16_t)-1, seq_nt6_table[(uint32_t)(tstr[k])], 1, 1); 
        }
        return;
    }
    if(ql == tl && ql == 1) {
        if(qstr[0] == tstr[0]) push_trace_bp_f(sc, 0, (uint16_t)-1, (uint16_t)-1, 1, 1); 
        else push_trace_bp_f(sc, 1, seq_nt6_table[(uint32_t)(qstr[0])], seq_nt6_table[(uint32_t)(tstr[0])], 1, 1); 
        return;
    }



    if(ql >= tl) {
        maxl = ql; minl = tl;
    } else {
        maxl = tl; minl = ql;
    }
    f = 0;
            
    diff = 31;
    if(diff > (maxl - minl)) {
        if(diff > maxl) diff = maxl; 
        diff0 = diff; clear_align(*exz);
        cal_exz_global(tstr, tl, qstr, ql, diff, exz);
        if(is_align(*exz)) f = 1; 
    }

    if(!f) {
        diff = 63;
        if(diff > (maxl - minl)) {
            if(diff > maxl) diff = maxl; 
            if(diff > diff0) {
                diff0 = diff; clear_align(*exz);
                cal_exz_global(tstr, tl, qstr, ql, diff, exz);
                if(is_align(*exz)) f = 1; 
            }
        }
    }

    if(!f) {
        if((maxn > maxl) && (maxe > (maxl - minl))) {
            diff = maxl * e_rate;
            if(diff < 1) diff = 1;
            if(diff > maxe) diff = maxe;
            if(diff > diff0) {
                diff0 = diff; clear_align(*exz);
                cal_exz_global(tstr, tl, qstr, ql, diff, exz);
                if(is_align(*exz)) f = 1; 
            }
        }
    }

    // fprintf(stderr, "[M::%s] f::%lu, err::%d\n", __func__, f, exz->err);

    uint32_t on; uint16_t op;
    if(f) {
        
        for (ck = qk = tk = 0; ck < exz->cigar.n;) {
            wq[0] = qk; wt[0] = tk;
            ck = pop_trace(&(exz->cigar), ck, &op, &on);
            if(op!=2) qk += on;
            if(op!=3) tk += on; 
            wq[1] = qk; wt[1] = tk;

            if(op == 0) {
                push_trace_bp_f(sc, op, (uint16_t)-1, (uint16_t)-1, on, 1); 
            } else if(op == 1) {
                for (k = 0; k < on; k++) {
                    push_trace_bp_f(sc, op, seq_nt6_table[(uint32_t)(qstr[wq[0]+k])], seq_nt6_table[(uint32_t)(tstr[wt[0]+k])], 1, 1); 
                }
            } else if(op == 2) {
                for (k = 0; k < on; k++) {
                    push_trace_bp_f(sc, op, (uint16_t)-1, seq_nt6_table[(uint32_t)(tstr[wt[0]+k])], 1, 1); 
                }
            } else if(op == 3) {
                for (k = 0; k < on; k++) {
                    push_trace_bp_f(sc, op, seq_nt6_table[(uint32_t)(qstr[wq[0]+k])], (uint16_t)-1, 1, 1); 
                }
            }
        }
    } else {
        if(ql > 0) {
            op = 3; on = ql;
            for (k = 0; k < on; k++) {
                push_trace_bp_f(sc, op, seq_nt6_table[(uint32_t)(qstr[k])], (uint16_t)-1, 1, 1); 
            }
        }

        if(tl > 0) {
            op = 2; on = tl;
            for (k = 0; k < on; k++) {
                push_trace_bp_f(sc, op, (uint16_t)-1, seq_nt6_table[(uint32_t)(tstr[k])], 1, 1); 
            }
        }
        
    }
}

void cal_updated_trace_len(asg16_v *sc, uint64_t *ql, uint64_t *tl)
{
    uint64_t ck = 0, qk = 0, tk = 0; uint32_t len; uint16_t c, bq, bt;
    while (ck < sc->n) {
        ck = pop_trace_bp_f(sc, ck, &c, &bq, &bt, &len);
        if(c != 2) qk += len;
        if(c != 3) tk += len;
    }
    *ql = qk; *tl = tk;
}

void gen_updated_trace(asg16_v *qcc, asg16_v *tcc, asg16_v *tcc_res, char *qstr, uint64_t ql, char *tstr, uint64_t tl, asg64_v *srt, bit_extz_t *exz, uint64_t rid)
{
    uint64_t k, ck, qk, tk, wq[2], wt[2], old_dp, dp, s, e, srt_n, si, ei, so, os, oe, *qd, *td, qs, qe, ts, te, q0, t0; 
    asg16_v *cc; uint32_t len; uint16_t c, bq, bt;

    srt->n = 0;

    cc = qcc;
    ck = qk = tk = wq[0] = wq[1] = wt[0] = wt[1] = 0;
    while (ck < cc->n) {
        wq[0] = qk; wt[0] = tk;
        ck = pop_trace_bp_f(cc, ck, &c, &bq, &bt, &len);
        if(c != 2) qk += len;
        if(c != 3) tk += len;
        wq[1] = qk; wt[1] = tk;
        if(c != 0) continue;

        kv_push(uint64_t, (*srt), (wq[0]<<1));
        kv_push(uint64_t, (*srt), ((wq[1]<<1)|1));
    }

    cc = tcc;
    ck = qk = tk = wq[0] = wq[1] = wt[0] = wt[1] = 0;
    while (ck < cc->n) {
        wq[0] = qk; wt[0] = tk;
        ck = pop_trace_bp_f(cc, ck, &c, &bq, &bt, &len);
        if(c != 2) qk += len;
        if(c != 3) tk += len;
        wq[1] = qk; wt[1] = tk;
        if(c != 0) continue;

        kv_push(uint64_t, (*srt), (wt[0]<<1));
        kv_push(uint64_t, (*srt), ((wt[1]<<1)|1));
    }

    // fprintf(stderr, "[M::%s] rid::%lu, ql::%lu, tl::%lu\n", __func__, rid, ql, tl);

    radix_sort_ec64(srt->a, srt->a + srt->n);
    for (k = 0, dp = e = srt_n = 0, s = (uint64_t)-1; k < srt->n; k++) {
        old_dp = dp;
        //if a[k] is qe
        if (srt->a[k]&1) --dp;
        else ++dp;

        if(old_dp >= 2 && s != (uint64_t)-1) {
            e = srt->a[k]>>1;
            if(e > s) srt->a[srt_n++] = ((s<<32)|(e));
        }

        s = (uint64_t)-1;
        if(dp >= 2) s = srt->a[k]>>1;
        // if (old_dp < 2 && dp >= 2) {///old_dp < dp, a[k] is qs
        //     s = srt->a[k]>>1;
        // } else if (old_dp >= 2 && dp < 2) {///old_dp > dp, a[k] is qe
        //     e = srt->a[k]>>1;
        //     if(e > s) {
        //         srt->a[srt_n++] = ((s<<32)|(e));
        //         fprintf(stderr, "[M::%s] x::[%lu,\t%lu)\n", __func__, s, e);
        //     }
        // }
    }

    
    // for (k = 0; k < srt_n; k++) {
    //     fprintf(stderr, "[M::%s] x[%lu]::[%lu,\t%u)\n", __func__, k, srt->a[k]>>32, (uint32_t)srt->a[k]);
    // }
    

    if(srt_n > 0) {
        srt->n = srt_n;

        cc = qcc;
        k = ck = qk = tk = wq[0] = wq[1] = wt[0] = wt[1] = 0;
        while (ck < cc->n) {
            wq[0] = qk; wt[0] = tk;
            ck = pop_trace_bp_f(cc, ck, &c, &bq, &bt, &len);
            if(c != 2) qk += len;
            if(c != 3) tk += len;
            wq[1] = qk; wt[1] = tk;
            if(c != 0) continue;
            si = wq[0]; ei = wq[1];
            so = wt[0]; ///eo = wt[1];

            for (; (k > 0) && ((k >= srt_n) || (((uint32_t)(srt->a[k])) > si)); k--);
            for (; k < srt_n; k++) {
                s = srt->a[k]>>32; e = (uint32_t)(srt->a[k]);
                if(s >= ei) break;
                if(s >= si && e <= ei) {
                    os = so + s - si;
                    oe = so + e - si; 
                    kv_push(uint64_t, (*srt), ((os<<32)|(oe)));
                }
            }
        }
        assert(srt->n == (srt_n<<1));

        cc = tcc;
        k = ck = qk = tk = wq[0] = wq[1] = wt[0] = wt[1] = 0;
        while (ck < cc->n) {
            wq[0] = qk; wt[0] = tk;
            ck = pop_trace_bp_f(cc, ck, &c, &bq, &bt, &len);
            if(c != 2) qk += len;
            if(c != 3) tk += len;
            wq[1] = qk; wt[1] = tk;
            if(c != 0) continue;
            si = wt[0]; ei = wt[1];
            so = wq[0]; ///eo = wq[1];

            for (; (k > 0) && ((k >= srt_n) || (((uint32_t)(srt->a[k])) > si)); k--);
            for (; k < srt_n; k++) {
                s = srt->a[k]>>32; e = (uint32_t)(srt->a[k]);
                // fprintf(stderr, "######[M::%s] t[%lu]::[%lu,\t%lu) i::[%lu,\t%lu)\n", __func__, k, s, e, si, ei);
                if(s >= ei) break;
                if(s >= si && e <= ei) {
                    os = so + s - si;
                    oe = so + e - si; 
                    kv_push(uint64_t, (*srt), ((os<<32)|(oe)));
                    // fprintf(stderr, "[M::%s] ******\n", __func__);
                }
            }
        }
        // if(!(srt->n == (srt_n*3))) {
        //     fprintf(stderr, "[M::%s] rid::%lu, srt_n::%lu, srt->n::%lu\n", __func__, rid, srt_n, (uint64_t)srt->n);
        // }
        assert(srt->n == (srt_n*3));

        
        tcc_res->n = 0; ///reset tcc
        qd = srt->a + srt_n; td = srt->a + srt_n + srt_n; uint64_t nl = 0;///, dbg_ql, dbg_tl;
        for (k = q0 = t0 = 0; k < srt_n; k++) {
            qs = qd[k]>>32; qe = (uint32_t)qd[k];
            ts = td[k]>>32; te = (uint32_t)td[k];
            nl += qs - qe;
            // assert((qe - qs) == (te - ts));
            // assert(!memcmp(qstr + qs, tstr + ts, sizeof((*qstr))*(qe - qs)));

            // if(t0 > ts || q0 > qs) {
            //     fprintf(stderr, "[M::%s] rid::%lu, ql::%lu, tl::%lu\n", __func__, rid, ql, tl);
            // }
            // fprintf(stderr, "[M::%s] qseq::[%lu,%lu), ql::%lu, tseq::[%lu,%lu), tl::%lu\n", __func__, t0, ts, tl, q0, qs, ql);

            gen_cc_fly(tcc_res, tstr + t0, ts - t0, qstr + q0, qs - q0, exz, 0.25, MAX_SIN_L, MAX_SIN_E);

            // cal_updated_trace_len(tcc_res, &dbg_ql, &dbg_tl);
            // assert(dbg_ql == ts && dbg_tl == qs);

            // fprintf(stderr, "******\n");
            
            push_trace_bp_f(tcc_res, 0, (uint16_t)-1, (uint16_t)-1,  qe - qs, 1);

            // cal_updated_trace_len(tcc_res, &dbg_ql, &dbg_tl);
            // assert(dbg_ql == te && dbg_tl == qe);

            q0 = qe; t0 = te;
        }

        qs = ql; ts = tl;
        // fprintf(stderr, "[M::%s] qseq::[%lu,%lu), ql::%lu, tseq::[%lu,%lu), tl::%lu\n", __func__, t0, ts, tl, q0, qs, ql);
        gen_cc_fly(tcc_res, tstr + t0, ts - t0, qstr + q0, qs - q0, exz, 0.25, MAX_SIN_L, MAX_SIN_E);
        // if(!(dbg_ql == ts && dbg_tl == qs)) {
        //     fprintf(stderr, "[M::%s] rid::%lu, qseq::[%lu,%lu), ql::%lu, tseq::[%lu,%lu), tl::%lu\n", __func__, rid, t0, ts, tl, q0, qs, ql);
        // }
        // cal_updated_trace_len(tcc_res, &dbg_ql, &dbg_tl);
        // assert(dbg_ql == ts && dbg_tl == qs);

        // fprintf(stderr, "[M::%s] srt_n::%lu, nl::%lu, ql::%lu, tl::%lu\n", __func__, srt_n, nl, ql, tl);
    } else {
        gen_cc_fly(tcc_res, tstr, tl, qstr, ql, exz, 0.25, MAX_SIN_L, MAX_SIN_E);
    }
}

void update_scb(All_reads *R_INF, asg16_v *scc, asg16_v *scb, asg16_v *scb_res, UC_Read *qu, UC_Read *tu, asg64_v *srt, bit_extz_t *exz, uint64_t rid)
{
    char *qstr = NULL, *tstr = NULL; uint64_t ql = 0, tl = 0;
    uint64_t ck, qk, tk, k, wq[2], wt[2]; uint32_t len; uint16_t c, bq, bt;
    gen_ori_seq0(qu->seq, qu->length, tu, scb, rid); ///tstr = tu->seq; tl = tu->length;

    ck = qk = tk = 0; ql = qu->length;
    while (ck < scc->n) {
        wq[0] = qk; wt[0] = tk;
        ck = pop_trace_bp_f(scc, ck, &c, &bq, &bt, &len);
        if(c != 2) qk += len;
        if(c != 3) tk += len;
        wq[1] = qk; wt[1] = tk;
    }
    assert(qk == ql);
    tl = tk; resize_UC_Read(qu, ql + tl);
    qstr = qu->seq; tstr = qu->seq + ql;

    ck = 0; qk = tk = 0; 
    while (ck < scc->n) {
        wq[0] = qk; wt[0] = tk;
        ck = pop_trace_bp_f(scc, ck, &c, &bq, &bt, &len);
        if(c != 2) qk += len;
        if(c != 3) tk += len;
        wq[1] = qk; wt[1] = tk;
        // if(xk > (uint32_t)p->z.length) fprintf(stderr, "[M::%s] xk::%u, len::%u, c::%u, rid::%ld\n", __func__, xk, (uint32_t)p->z.length, c, i);
        if(c == 0) {
            memcpy(tstr + wt[0], qstr + wq[0], (wq[1]-wq[0])*sizeof((*qstr)));
        } else if(c == 1 || c == 2) {
            for (k = wt[0]; k < wt[1]; k++) tstr[k] = s_H[bt]; 
        }
        // if(i == 700) fprintf(stderr, "|%u%c(%c)(x::%u)(y::%u)", len, cm[c], ((c==1)||(c==2))?(cc[b]):('*'), wx[1], wy[1]); // s_H
    }

    qstr = tstr; ql = tl;
    tstr = tu->seq; tl = tu->length;

    // fprintf(stderr, "\n[M::%s] ql::%lu, tl::%lu, rid::%lu\n", __func__, ql, tl, rid);

    gen_updated_trace(scc, scb, scb_res, qstr, ql, tstr, tl, srt, exz, rid);

    

    ///debug
    // resize_UC_Read(tu, ql + tl);
    // memcpy(tu->seq + tl, qstr, ql); tstr = tu->seq; qstr = tu->seq + tl;

    // resize_UC_Read(qu, ql + tl);
    // memcpy(qu->seq, tu->seq, ql + tl); tstr = qu->seq; qstr = qu->seq + tl;

    // gen_ori_seq0(qstr, ql, tu, scb_res, rid);
    // assert(memcmp(tstr, tu->seq, tl) == 0);
    
}

uint32_t is_well_cal(asg64_v *idx, ma_hit_t_alloc *in0, ma_hit_t_alloc *in1, int64_t ql, int64_t occ_exact)
{
    ma_hit_t_alloc *paf = NULL; uint64_t k, s, e, vn; ma_hit_t *z; idx->n = 0;
    int64_t dp, old_dp, st = 0, ed;

    paf = in0;
    for (k = 0; k < paf->length; k++) {
        z = &(paf->buffer[k]);
        s = ((uint32_t)(z->qns)); e = z->qe;
        kv_push(uint64_t, (*idx), (s<<1));
        kv_push(uint64_t, (*idx), (e<<1)|1);
    }

    paf = in1;
    for (k = 0; k < paf->length; k++) {
        z = &(paf->buffer[k]);
        s = ((uint32_t)(z->qns)); e = z->qe;
        kv_push(uint64_t, (*idx), (s<<1));
        kv_push(uint64_t, (*idx), (e<<1)|1);
    }

    radix_sort_ec64(idx->a, idx->a + idx->n); vn = idx->n;
    for (k = 0, dp = 0, st = ed = 0; k < vn; ++k) {
        old_dp = dp;
        ///if a[j] is qe
        if (idx->a[k]&1) --dp;
        else ++dp;

        ed = idx->a[k]>>1;
        if((ed > st) && ((old_dp + 1) < occ_exact)) return 0;///+1 for self

        st = ed;
    }


    ed = ql; old_dp = dp;
    if((ed > st) && ((old_dp + 1) < occ_exact)) return 0;///+1 for self
    
    return 1;
}

static void worker_hap_dc_ec0(void *data, long i, int tid)
{
    // if(i == 6) fprintf(stderr, "-mm-[M::%s]\tqn::%u::%.*s\tf[i]::%u\n", __func__, (uint32_t)(i), (int)Get_NAME_LENGTH(R_INF, i), Get_NAME((R_INF), i), scc.f[i]);
    if(scc.f[i]) {
        scc.a[i].n = 0; sca.a[i].n = 0;
        // push_trace_bp(&(scc.a[i]), 0, (uint16_t)-1, Get_READ_LENGTH(R_INF, i), 0);
        push_trace_bp_f(&(scc.a[i]), 0, (uint16_t)-1, (uint16_t)-1, Get_READ_LENGTH(R_INF, i), 0);
        return;
    }
    ec_ovec_buf_t0 *b = &(((ec_ovec_buf_t*)data)->a[tid]);
    uint32_t high_occ = asm_opt.hom_cov * (2.0 - HA_KMER_GOOD_RATIO);
    uint32_t low_occ = asm_opt.hom_cov * HA_KMER_GOOD_RATIO;
    asg64_v buf0; overlap_region *aux_o = NULL; uint32_t qlen = 0;
    // overlap_region *aux_o = NULL; asg64_v buf0;

    // gen_ovlst_paf(&(R_INF.paf[i]), &(R_INF.reverse_paf[i]), &(b->v64));
    // if(i != 181) return;

    // fprintf(stderr, "-mm-[M::%s]\tqn::%u::%.*s\n", __func__, (uint32_t)(i), (int)Get_NAME_LENGTH(R_INF, i), Get_NAME((R_INF), i));

    recover_UC_Read(&b->self_read, &R_INF, i); qlen = b->self_read.length;


    /**
    if(is_well_cal(&b->v64, &(R_INF.paf[i]), &(R_INF.reverse_paf[i]), b->self_read.length, 4)) {
        // aux_o = h_ec_lchain_re1(b->ab, i, &b->self_read, &b->ovlp_read, asm_opt.mz_win, asm_opt.k_mer_length, &R_INF, &b->olist, &b->clist, &b->exz, &b->v16, &b->v64, 0.02, 1, NULL, NULL, &(b->sp), &high_occ, &low_occ, 1, 1, 0, 2, UINT32_MAX, &(R_INF.paf[i]), &(R_INF.reverse_paf[i]));
        aux_o = h_ec_lchain_re2(b->ab, i, &b->self_read, &b->ovlp_read, asm_opt.mz_win, asm_opt.k_mer_length, &R_INF, &b->olist, &b->clist, &b->exz, &b->v16, &b->v64, 0.02, 1, NULL, NULL, &(b->sp), &high_occ, &low_occ, 1, 1, 0, 2, UINT32_MAX, &(R_INF.paf[i]), &(R_INF.reverse_paf[i]));
    } else {
        aux_o = h_ec_lchain_re3(b->ab, i, &b->self_read, &b->ovlp_read, asm_opt.mz_win, asm_opt.k_mer_length, &R_INF, &b->olist, &b->clist, &b->exz, &b->v16, &b->v64, 0.02, 1, NULL, NULL, &(b->sp), &high_occ, &low_occ, 1, 1, 0, 2, UINT32_MAX, &(R_INF.paf[i]), &(R_INF.reverse_paf[i]));
        // fprintf(stderr, "[M::%s]\tqn::%u::%.*s\n", __func__, (uint32_t)(i), (int)Get_NAME_LENGTH(R_INF, i), Get_NAME((R_INF), i));
    }
    **/
   aux_o = h_ec_lchain_re2(b->ab, i, &b->self_read, &b->ovlp_read, asm_opt.mz_win, asm_opt.k_mer_length, &R_INF, &b->olist, &b->clist, &b->exz, &b->v16, &b->v64, 0.02, 1, NULL, NULL, &(b->sp), &high_occ, &low_occ, 1, 1, 0, 2, UINT32_MAX, &(R_INF.paf[i]), &(R_INF.reverse_paf[i]));
    
    ////for debug
    // scc.a[i].n = 0;
    // push_trace_bp(&(scc.a[i]), 0, (uint16_t)-1, Get_READ_LENGTH(R_INF, i), 0);

    // return;

    // aux_o = h_ec_lchain_re(b->ab, i, b->self_read.seq, b->self_read.length, &b->ovlp_read, asm_opt.mz_win, asm_opt.k_mer_length, &R_INF, &b->olist, &b->clist, &b->exz, &b->v16, 0.02, 1, NULL, NULL, &(b->sp), &high_occ, &low_occ, 1, 1, 0, 2, UINT32_MAX, &(R_INF.paf[i]), &(R_INF.reverse_paf[i]));
    
    b->cnt[0] += b->self_read.length;

    copy_asg_arr(buf0, b->sp); 
    rphase_hc(&b->olist, &R_INF, &b->hap, &b->self_read, &b->ovlp_read, &b->pidx, &b->v64, &buf0, 0, WINDOW_MAX_SIZE, b->self_read.length, 1/**, 1**/, i, (asm_opt.is_ont)?HPC_PL:0, asm_opt.is_ont, ((asm_opt.is_ont)?&(b->clist.chainDP):NULL), ((asm_opt.is_sc)?&(b->v8q):NULL), ((asm_opt.is_sc)?&(b->v8t):NULL));
    copy_asg_arr(b->sp, buf0); 

    copy_asg_arr(buf0, b->sp);
    b->cnt[1] += wcns_gen(&b->olist, &R_INF, &b->self_read, &b->ovlp_read, &b->exz, &b->pidx, &b->v64, &buf0, 0, 512, b->self_read.length, 3, 0.500001, aux_o, &b->v32, &b->cns, 256, i);
    copy_asg_arr(b->sp, buf0);

    push_nec_re(aux_o, &(scc.a[i]));
    update_scb(&R_INF, &(scc.a[i]), &(scb.a[i]), &(sca.a[i]), &b->self_read, &b->ovlp_read, &b->v64, &b->exz, i);
    
    push_ne_ovlp(&(R_INF.paf[i]), &b->olist, 1, &R_INF, &(scc.a[i])/**, i, &b->self_read, &b->ovlp_read**/);
    push_ne_ovlp(&(R_INF.reverse_paf[i]), &b->olist, 2, &R_INF, NULL/**, i, NULL, NULL**/);

    check_well_cal(&(scc.a[i]), &b->v64, &(R_INF.paf[i].is_fully_corrected), &(R_INF.paf[i].is_abnormal), qlen, (MIN_COVERAGE_THRESHOLD*2), &(R_INF.paf[i]));
    R_INF.trio_flag[i] = AMBIGU;

    refresh_ec_ovec_buf_t0(b, REFRESH_N);
}

void get_origin_ec_coor(asg16_v *ec, uint64_t *ts, uint64_t *te)
{
    uint64_t ts0 = *ts, te0 = *te, qk = 0, tk = 0, ck = 0, wq[2], wt[2]; uint16_t op, bq, bt, f = 0; uint32_t cl;
    while (ck < ec->n) {
        wq[0] = qk; wt[0] = tk; 
        ck = pop_trace_bp_f(ec, ck, &op, &bq, &bt, &cl);
        if(op != 2) qk += cl;
        if(op != 3) tk += cl;
        wq[1] = qk; wt[1] = tk;  
        if((op == 0) && (wt[0] <= ts0) && (wt[1] >= te0)) {
            (*ts) = wq[0] + ts0 - wt[0]; 
            (*te) = wq[0] + te0 - wt[0]; 
            f = 1;
            break;
        }
    }
    assert(f);
}

static void update_scb0(void *data, long i, int tid)
{
    if(sca.a[i].n) {
        kv_resize(uint16_t, scb.a[i], sca.a[i].n); scb.a[i].n = sca.a[i].n;
        memcpy(scb.a[i].a, sca.a[i].a, scb.a[i].n*sizeof((*(sca.a[i].a))));
    }

    // return;

    if(!scc.f[i]) return;

    ma_hit_t_alloc *ov, *os; uint64_t k, kr, qn, tn, ql, tl, qs, qe, ts, te; ma_hit_t *z, *r;
    uint64_t ck; uint16_t op, bq, bt; uint32_t cl;
    
    ov = &(R_INF.paf[i]);
    for (k = 0; k < ov->length; k++) {
        z = &(ov->buffer[k]);
        qn = z->qns>>32; tn = z->tn;
        if(scc.f[tn]) continue;

        os = &(R_INF.paf[tn]);
        for (kr = 0; kr < os->length; kr++) {
            if(os->buffer[kr].tn == qn) {
                r = &(os->buffer[kr]);
                break;
            }
        }
        if(kr >= os->length) continue;
        // if(!(r->el)) continue;

        qs = r->ts; qe = r->te;
        ts = (uint32_t)r->qns; te = r->qe;
        z->el = r->el; z->rev = r->rev; z->ml = r->ml; z->no_l_indel = r->no_l_indel;
        if(z->el) {
            get_origin_ec_coor(&(scc.a[tn]), &ts, &te);
        } 

        z->qns = qn; 
        z->qns = z->qns << 32;
        if(z->rev) {
            ql = Get_READ_LENGTH(R_INF, qn);
            tl = ck = 0;
            while (ck < scc.a[tn].n) {
                ck = pop_trace_bp_f(&(scc.a[tn]), ck, &op, &bq, &bt, &cl);
                if(op != 2) tl += cl;
            }
            z->qns = z->qns | (ql - qe);
            z->qe = ql - qs;
            z->ts = tl - te;
            z->te = tl - ts;
        } else {
            z->qns = z->qns | qs;
            z->qe = qe;
            z->ts = ts;
            z->te = te;
        }
    }
}

void dbg_rsc(char *str0, uint64_t l0, char *str1, uint64_t l1, asg16_v *sc, char *real, uint32_t id)
{
    uint64_t ck, xk, yk, k, wx[2], wy[2]; uint32_t len; uint16_t c, bq, bt; 
    
    ck = xk = yk = 0;
    while (ck < sc->n) {
        wx[0] = xk; wy[0] = yk;
        ck = pop_trace_bp_f(sc, ck, &c, &bq, &bt, &len);
        if(c != 2) xk += len;
        if(c != 3) yk += len;
        wx[1] = xk; wy[1] = yk;
    }
    assert(xk == l0);

    // char cm[4], cc[4]; 
    // cm[0] = 'M'; cm[1] = 'S'; cm[2] = 'I'; cm[3] = 'D';
    // cc[0] = 'A'; cc[1] = 'C'; cc[2] = 'G'; cc[3] = 'T';

    ck = xk = yk = 0;
    while (ck < sc->n) {
        wx[0] = xk; wy[0] = yk;
        ck = pop_trace_bp_f(sc, ck, &c, &bq, &bt, &len);
        if(c != 2) xk += len;
        if(c != 3) yk += len;
        wx[1] = xk; wy[1] = yk;

        if(c == 0) {
            memcpy(str0 + wx[0], str1 + wy[0], (wx[1]-wx[0])*sizeof((*str1)));
        } else if(c == 1 || c == 3) {
            for (k = wx[0]; k < wx[1]; k++) str0[k] = s_H[bq]; 
        }
        // fprintf(stderr, "%u%c(%c)(x::[%lu,%ld))(y::[%lu,%ld))\n", len, cm[c], ((c==1)||(c==2))?(cc[bt]):('*'), wx[0], wx[1], wy[0], wy[1]); // s_H
    }

    if(memcmp(str0, real, l0*sizeof((*str0)))) {
        fprintf(stderr, "-0-[M::%s]\tid::%u\n", __func__, id);
        // for (k = 0; k < l0 && str0[k] == real[k]; k++);
        
        // fprintf(stderr, "-0-[M::%s]\tid::%u\tk::%lu\tl0::%lu\tnc::%c\toc::%c\n", __func__, id, k, l0, str0[k], real[k]);
        // exit(1);
    } else {
        // fprintf(stderr, "-1-[M::%s]\tid::%u\n", __func__, id);
    }

}

static void worker_sl_ec(void *data, long i, int tid)
{
    // if(i != 0) return;

	sl_v *p = &(((sl_v*)data)[tid]); uint8_t *oa = NULL; char *na = NULL; uint64_t tqual, wqual;
    uint32_t ci = 0, len, xk, yk, wx[2], wy[2], k, Nn, yn = 0, tot_e; uint16_t c, bq, bt; 


    ci = 0; xk = yk = 0; tot_e = 0;
    while (ci < scc.a[i].n) {
        // ci = pop_trace_bp(&scc.a[i], ci, &c, &b, &len);
        ci = pop_trace_bp_f(&scc.a[i], ci, &c, &bq, &bt, &len);
        if(c != 3) yk += len;
        if(c != 0) tot_e += len;
        // fprintf(stderr, "|%u%c(%c)", len, cm[c], ((c==1)||(c==2))?(cc[b]):('*')); // s_H
    }
    if(tot_e == 0) return;///no change

    yn = yk; yk++; kv_resize(char, (*p), yk); p->a[yn] = '\0';
    recover_UC_Read(&p->z, &R_INF, i); ///b->z.length


    // char cm[4], cc[4]; 
    // cm[0] = 'M'; cm[1] = 'S'; cm[2] = 'I'; cm[3] = 'D';
    // cc[0] = 'A'; cc[1] = 'C'; cc[2] = 'G'; cc[3] = 'T';

    ci = 0; xk = yk = 0; Nn = 0;
    while (ci < scc.a[i].n) {
        wx[0] = xk; wy[0] = yk;
        // ci = pop_trace_bp(&scc.a[i], ci, &c, &b, &len);
        ci = pop_trace_bp_f(&scc.a[i], ci, &c, &bq, &bt, &len);
        if(c != 2) xk += len;
        if(c != 3) yk += len;
        wx[1] = xk; wy[1] = yk;
        // if(xk > (uint32_t)p->z.length) fprintf(stderr, "[M::%s] xk::%u, len::%u, c::%u, rid::%ld\n", __func__, xk, (uint32_t)p->z.length, c, i);
        if(c == 0) {
            // memcpy(p->a + wy[0], p->z.seq + wx[0], (wx[1]-wx[0])*sizeof((*(p->a))));
            for (; wx[0] < wx[1]; wx[0]++, wy[0]++) {
                p->a[wy[0]] = p->z.seq[wx[0]];
                if(p->a[wy[0]] == 'N') Nn++;
            }
        } else if(c == 1 || c == 2) {
            for (k = wy[0]; k < wy[1]; k++) {
                p->a[k] = s_H[bt]; 
                if(p->a[k] == 'N') Nn++;
            }
        }

        // if(i == 700) fprintf(stderr, "|%u%c(%c)(x::%u)(y::%u)", len, cm[c], ((c==1)||(c==2))?(cc[b]):('*'), wx[1], wy[1]); // s_H
    }

    // if(i == 700) fprintf(stderr, "|\n");
    if(asm_opt.is_sc) retrive_bqual(&(p->q), NULL, i, -1, -1, 0, sc_bn);


    if (R_INF.read_size[i] < yn) {
		R_INF.read_size[i] = yn;
		REALLOC(R_INF.read_sperate[i], R_INF.read_size[i]/4+1);
        if(asm_opt.is_sc) REALLOC(R_INF.rsc[i], ((R_INF.read_size[i]/sc_bn) + ((R_INF.read_size[i]%sc_bn)?1:0)));
	}
	R_INF.read_length[i] = yn;
    // if(Nn > 0) fprintf(stderr, "[M::%s] Nn->%u\n", __func__, Nn);

    // for (k = 0; k < yn; k++) {
    //     c = seq_nt6_table[(uint8_t)p->a[k]]; 
    //     if (c >= 4) { 
    //         fprintf(stderr, "[M::%s] Nn->%u, yn::%u, xn::%lld, k::%u, str::%c, c::%u, rid::%ld\n", __func__, Nn, yn, p->z.length, k, p->a[k], c, i);
    //     }
    // }


    ///debug
    // resize_UC_Read(&p->z, p->z.length * 2); 
    // dbg_rsc(p->z.seq + p->z.length, p->z.length, p->a, yn, &(scc.a[i]), p->z.seq, i);


	ha_compress_base(Get_READ(R_INF, i), p->a, yn, &R_INF.N_site[i], Nn);
    if(asm_opt.is_sc) {
        oa = p->q.a; na = p->a; 
        ci = 0; xk = yk = 0; Nn = 0;
        while (ci < scc.a[i].n) {
            wx[0] = xk; wy[0] = yk;
            ci = pop_trace_bp_f(&scc.a[i], ci, &c, &bq, &bt, &len);
            if(c != 2) xk += len;
            if(c != 3) yk += len;
            wx[1] = xk; wy[1] = yk;
            if(c == 0 || c == 1) {
                memcpy(na + wy[0], oa + wx[0], (wx[1]-wx[0])*sizeof((*oa)));
            } else if(c == 2) {
                get_wqual(i, wx[0], 0, NULL, oa, sc_wn, &tqual, &wqual);
                for (k = wy[0]; k < wy[1]; k++) na[k] = wqual;
            }
        }
        assert(yk == yn);
        ha_compress_qual_bit(Get_QUAL(R_INF, i), na, yn, sc_bn);
    }
}

uint64_t cal_ec_multiple(ec_ovec_buf_t *b, uint64_t n_thre, uint64_t n_a, uint64_t *r_base)
{
    double tt0 = yak_realtime_0();
    uint64_t k, num_base = 0, num_correct = 0; (*r_base) = 0;

    if(!(scc.a)) {
        scc.n = scc.m = n_a; CALLOC(scc.a, n_a); CALLOC(scc.f, n_a);
    }

    if(!(scb.a)) {
        scb.n = scb.m = n_a; CALLOC(scb.a, n_a);
    }

    for (k = 0; k < n_thre; ++k) b->a[k].cnt[0] = b->a[k].cnt[1] = 0;

    kt_for(n_thre, worker_hap_ec, b, n_a);///debug_for_fix

    for (k = 0; k < n_thre; ++k) {
        num_base += b->a[k].cnt[0];
        num_correct += b->a[k].cnt[1];
    }

    // fprintf(stderr, "\n[M::%s] # reads->%lu\n", __func__, n_a);
    // fprintf(stderr, "[M::%s] # input bases->%lu\n", __func__, num_base);
    // fprintf(stderr, "[M::%s] # corrected bases->%lu\n", __func__, num_correct);
    // fprintf(stderr, "[M::%s::%.3f] running time\n", __func__, yak_realtime_0()-tt0);
    fprintf(stderr, "[M::pec::%.3f] # bases: %lu; # corrected bases: %lu\n", yak_realtime_0()-tt0, num_base, num_correct);

    (*r_base) = num_base;
    return num_correct;
}


void ha_print_ovlp_stat_1(ec_ovec_buf_t *b, uint64_t n_thre, uint64_t n_a)
{
    double tt0 = yak_realtime_0();
    uint64_t k, forward, reverse, strong, weak, exact, no_l_indel;

    forward = reverse = strong = weak = exact = no_l_indel = 0;

    ///calculate overlaps
    for (k = 0; k < n_thre; ++k) {
        b->a[k].cnt[0] = b->a[k].cnt[1] = b->a[k].cnt[2] = b->a[k].cnt[3] = b->a[k].cnt[4] = b->a[k].cnt[5] = 0;
    }

    kt_for(n_thre, worker_hap_dc_ec_gen, b, n_a);

    for (k = 0; k < n_thre; ++k) {
        forward += b->a[k].cnt[0];
        reverse += b->a[k].cnt[1];
        strong += b->a[k].cnt[2];
        weak += b->a[k].cnt[3];
        exact += b->a[k].cnt[4];
        no_l_indel += b->a[k].cnt[5];
    }

    fprintf(stderr, "[M::%s] # overlaps: %lu\n", __func__, forward);
	fprintf(stderr, "[M::%s] # strong overlaps: %lu\n", __func__, strong);
	fprintf(stderr, "[M::%s] # weak overlaps: %lu\n", __func__, weak);
	fprintf(stderr, "[M::%s] # exact overlaps: %lu\n", __func__, exact); // this seems not right
	fprintf(stderr, "[M::%s] # inexact overlaps: %lu\n", __func__, forward - exact);
	fprintf(stderr, "[M::%s] # overlaps without large indels: %lu\n", __func__, no_l_indel);
	fprintf(stderr, "[M::%s] # reverse overlaps: %lu\n", __func__, reverse);
    fprintf(stderr, "[M::%s] # running time: %.3f\n", __func__, yak_realtime_0()-tt0);

    // fprintf(stderr, "\n[M::%s] # reads->%lu\n", __func__, n_a);
    // fprintf(stderr, "[M::%s] # corrected reads->%lu\n", __func__, rb);
    // fprintf(stderr, "[M::%s] # uncorrected reads->%lu\n", __func__, urb);
    // fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime_0()-tt0);
}

void ha_print_ovlp_stat_0(ec_ovec_buf_t *b, uint64_t n_thre, uint64_t n_a)
{
    double tt0 = yak_realtime_0();
    uint64_t k, forward, reverse, strong, weak, exact, no_l_indel;

    forward = reverse = strong = weak = exact = no_l_indel = 0;

    ///calculate overlaps
    for (k = 0; k < n_thre; ++k) {
        b->a[k].cnt[0] = b->a[k].cnt[1] = b->a[k].cnt[2] = b->a[k].cnt[3] = b->a[k].cnt[4] = b->a[k].cnt[5] = 0;
    }

    kt_for(n_thre, worker_hap_dc_ec_gen_new_idx, b, n_a);

    for (k = 0; k < n_thre; ++k) {
        forward += b->a[k].cnt[0];
        reverse += b->a[k].cnt[1];
        strong += b->a[k].cnt[2];
        weak += b->a[k].cnt[3];
        exact += b->a[k].cnt[4];
        no_l_indel += b->a[k].cnt[5];
    }

    fprintf(stderr, "[M::%s] # overlaps: %lu\n", __func__, forward);
	fprintf(stderr, "[M::%s] # strong overlaps: %lu\n", __func__, strong);
	fprintf(stderr, "[M::%s] # weak overlaps: %lu\n", __func__, weak);
	fprintf(stderr, "[M::%s] # exact overlaps: %lu\n", __func__, exact); // this seems not right
	fprintf(stderr, "[M::%s] # inexact overlaps: %lu\n", __func__, forward - exact);
	fprintf(stderr, "[M::%s] # overlaps without large indels: %lu\n", __func__, no_l_indel);
	fprintf(stderr, "[M::%s] # reverse overlaps: %lu\n", __func__, reverse);
    fprintf(stderr, "[M::%s] # running time: %.3f\n", __func__, yak_realtime_0()-tt0);

    // fprintf(stderr, "\n[M::%s] # reads->%lu\n", __func__, n_a);
    // fprintf(stderr, "[M::%s] # corrected reads->%lu\n", __func__, rb);
    // fprintf(stderr, "[M::%s] # uncorrected reads->%lu\n", __func__, urb);
    // fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime_0()-tt0);
}

uint64_t cal_sec_ec_multiple(ec_ovec_buf_t *b, uint64_t n_thre, uint64_t n_a, int64_t round)
{
    double tt0 = yak_realtime_0();
    uint64_t k, num_base, num_correct, rb, urb;
    num_base = num_correct = 0;
    
    ////counting
    rb = urb = 0;
    for (k = 0; k < n_thre; ++k) b->a[k].cnt[0] = b->a[k].cnt[1] = 0;

    kt_for(n_thre, worker_hap_dc_ec, b, n_a);///debug_for_fix
    
    for (k = 0; k < n_thre; ++k) {
        rb += b->a[k].cnt[0]; urb += b->a[k].cnt[1];
    }

    if(round >= 0) {
        if(!(sca.a)) {
            sca.n = sca.m = n_a; CALLOC(sca.a, n_a);
        }
        ////correct
        
        for (k = 0; k < n_thre; ++k) b->a[k].cnt[0] = b->a[k].cnt[1] = 0;
        
        kt_for(n_thre, worker_hap_dc_ec0, b, n_a);///debug_for_fix

        for (k = 0; k < n_thre; ++k) {
            num_base += b->a[k].cnt[0];
            num_correct += b->a[k].cnt[1];
        }

        kt_for(n_thre, update_scb0, b, n_a);
    }

    if(round >= 0) {
        fprintf(stderr, "[M::sec::%.3f] # bases: %lu; # corrected bases: %lu; # reads: %lu; # corrected reads: %lu\n", yak_realtime_0()-tt0, num_base, num_correct, rb, urb);
    } else {
        fprintf(stderr, "[M::sec::%.3f] # reads: %lu; # corrected reads: %lu\n", yak_realtime_0()-tt0, rb, urb);
    }

    // fprintf(stderr, "\n[M::%s] # reads->%lu\n", __func__, n_a);
    // fprintf(stderr, "[M::%s] # corrected reads->%lu\n", __func__, rb);
    // fprintf(stderr, "[M::%s] # uncorrected reads->%lu\n", __func__, urb);
    // if(round >= 0) {
    //     fprintf(stderr, "[M::%s] # input bases->%lu\n", __func__, num_base);
    //     fprintf(stderr, "[M::%s] # corrected bases->%lu\n", __func__, num_correct);
    //     fprintf(stderr, "[M::%s::%.3f] ==> round %ld\n", __func__, yak_realtime_0()-tt0, round);
    // }
    return num_correct;
}


void write_ec_reads(const char *suffix_ou)
{
    uint64_t k, strl; UC_Read qstr, tstr; char *nn = NULL, *str = NULL;
    init_UC_Read(&qstr); init_UC_Read(&tstr);
    MALLOC(nn, strlen(suffix_ou) + strlen(asm_opt.output_file_name) + 36);
    sprintf(nn, "%s.%s", asm_opt.output_file_name, suffix_ou);
    FILE *ou = fopen(nn, "w");
    free(nn);

    for (k = 0; k < R_INF.total_reads; k++) {
        recover_UC_Read(&qstr, &R_INF, k);
        if(scb.a) {
            gen_ori_seq0(qstr.seq, qstr.length, &tstr, &(scb.a[k]), k); str = tstr.seq; strl = tstr.length;
        } else {
            str = qstr.seq; strl = qstr.length;
        }

        fwrite(">", 1, 1, ou);
        fwrite(Get_NAME(R_INF, k), 1, Get_NAME_LENGTH(R_INF, k), ou);
        fwrite("\n", 1, 1, ou);
        fwrite(str, 1, strl, ou);
        fwrite("\n", 1, 1, ou);
    }

    fclose(ou); destory_UC_Read(&qstr); destory_UC_Read(&tstr);
}


void cal_ec_r(uint64_t n_thre, uint64_t round, uint64_t n_round, uint64_t n_a, uint64_t is_sv, uint64_t *tot_b, uint64_t *tot_e)
{
    // write_ec_reads("ec0.fa");

    ec_ovec_buf_t *b = NULL; uint64_t k, is_cr = (round&1);
    (*tot_b) = (*tot_e) = 0;


    b = gen_ec_ovec_buf_t(n_thre);
    (*tot_e) += cal_ec_multiple(b, n_thre, n_a, tot_b); ///exit(1);
    sl_ec_r(n_thre, n_a);

    for (k = 0; k < n_round; k++) {
        (*tot_e) += cal_sec_ec_multiple(b, n_thre, n_a, k);
        sl_ec_r(n_thre, n_a);
    }

    if(is_sv) kt_for(n_thre, worker_hap_dc_ec, b, n_a);///update overlaps

    if((!is_sv) || (is_sv && is_cr)) {
        kt_for(n_thre, worker_hap_post_rev, b, n_a);
    }

    // cal_sec_ec_multiple(b, n_thre, n_a, -1);

    // gen_sec_ec_multiple(b, n_thre, n_a);

    destroy_ec_ovec_buf_t(b);

    // write_ec_reads("ec16.fa");

    // uint64_t z;
    // for (z = 0; z < scc.n; z++) {
    //     if(scc.f[z]) continue;
    //     fprintf(stderr, "[M::%s]\tid::%lu::%.*s\n", __func__, z, (int)Get_NAME_LENGTH(R_INF, z), Get_NAME((R_INF), z));
    // }
}

void print_ov_dbg_paf(FILE *fp, char *ref_str, char *ref_id, int32_t ref_id_n, char *qry_str, char *qry_id, int32_t qry_id_n, uint64_t rs, uint64_t re, uint64_t rl, uint64_t qs, uint64_t qe, uint64_t ql, uint64_t rev, bit_extz_t *ez, char *ezh)
{
    uint64_t ci = 0; uint16_t c; uint32_t cl;
    fprintf(fp, "%.*s\t%lu\t%lu\t%lu\t", qry_id_n, qry_id, ql, qs, qe);
    fprintf(fp, "%c\t", "+-"[rev]);
    fprintf(fp, "%.*s\t%lu\t%lu\t%lu\t", ref_id_n, ref_id, rl, rs, re);
    fprintf(fp, "255\tcg:Z:");
    while (ci < ez->cigar.n) {
        ci = pop_trace(&(ez->cigar), ci, &c, &cl);
        fprintf(fp, "%u%c", cl, ezh[c]);
    }
    fprintf(fp, "\n");
}

static void *worker_ov_dbg_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
    cal_ec_r_dbg_t *p = (cal_ec_r_dbg_t*)data; char cm[4]; cm[0] = 'M'; cm[1] = 'M'; cm[2] = 'I'; cm[3] = 'D'; 
    // cal_ec_r_dbg_step_t
    if (step == 0) { // step 1: read a block of sequences
        cal_ec_r_dbg_step_t *s; CALLOC(s, 1);
        s->si = p->cn; p->cn += p->chunk_size; 
        if(p->cn > p->n_a) p->cn = p->n_a; s->ei = p->cn;
        if(s->si >= s->ei) free(s);
        else return s;
    } else if (step == 1) { // step 2: alignment
        cal_ec_r_dbg_step_t *s = (cal_ec_r_dbg_step_t*)in;
        s->buf = gen_ec_ovec_buf_t(p->n_thread); CALLOC(s->res, p->n_thread);  
        kt_for(p->n_thread, worker_hap_ec_dbg_paf, s, (s->ei - s->si));
        destroy_ec_ovec_buf_t(s->buf);
        return s;
    } else if (step == 2) { // step 3: dump
        cal_ec_r_dbg_step_t *s = (cal_ec_r_dbg_step_t*)in; uint64_t k, z; bit_extz_t ez; memset(&ez, 0, sizeof(ez));
        // UC_Read qu; UC_Read tu; init_UC_Read(&qu); init_UC_Read(&tu); 
        for (k = 0; k < p->n_thread; k++) {
            for (z = 0; z < s->res[k].n; z++) {
                ez.cigar.a = s->res[k].ec.a + s->res[k].a[z].bl; ez.cigar.n = ez.cigar.m = s->res[k].a[z].cc;

                // UC_Read_resize(qu, (s->res[k].a[z].qe - ((uint32_t)s->res[k].a[z].qns))); 
                // recover_UC_Read_sub_region(qu.seq, ((uint32_t)s->res[k].a[z].qns), (s->res[k].a[z].qe - ((uint32_t)s->res[k].a[z].qns)), 0, &R_INF, (s->res[k].a[z].qns>>32));

                // UC_Read_resize(tu, (s->res[k].a[z].te - s->res[k].a[z].ts));
                // recover_UC_Read_sub_region(tu.seq, s->res[k].a[z].ts, s->res[k].a[z].te - s->res[k].a[z].ts, 0, &R_INF, s->res[k].a[z].tn); 

                print_ov_dbg_paf(p->fp, NULL/**qu.seq**/, Get_NAME(R_INF, (s->res[k].a[z].qns>>32)), Get_NAME_LENGTH(R_INF, (s->res[k].a[z].qns>>32)), 
                            NULL/**tu.seq**/, Get_NAME(R_INF, (s->res[k].a[z].tn)), Get_NAME_LENGTH(R_INF, (s->res[k].a[z].tn)), (uint32_t)s->res[k].a[z].qns, s->res[k].a[z].qe, Get_READ_LENGTH(R_INF, (s->res[k].a[z].qns>>32)), s->res[k].a[z].ts, s->res[k].a[z].te, Get_READ_LENGTH(R_INF, (s->res[k].a[z].tn)), s->res[k].a[z].rev, &ez, cm);
            }
            free(s->res[k].a); free(s->res[k].ec.a);
        }
        free(s->res); free(s); ///destory_UC_Read(&qu); destory_UC_Read(&tu);
    }
    return 0;
}

void cal_ec_r_dbg(uint64_t n_thre, uint64_t n_a)
{
    char *paf = NULL; MALLOC(paf, (strlen(asm_opt.output_file_name)+64));
    sprintf(paf, "%s.ovlp.paf", asm_opt.output_file_name); 

    cal_ec_r_dbg_t sl; memset(&sl, 0, sizeof(sl));
    sl.n_thread = n_thre; sl.n_a = n_a; sl.chunk_size = 2000; sl.cn = 0; sl.fp = fopen(paf, "w");

    kt_pipeline(3, worker_ov_dbg_pipeline, &sl, 3);

    fclose(sl.fp); free(paf);
}

void destroy_cc_v(cc_v *z)
{
    uint64_t k;
    for (k = 0; k < z->m; k++) free(z->a[k].a);
    free(z->f); free(z->a);
    z->n = z->m = 0; z->f = NULL; z->a = NULL;
}

void cal_ov_r(uint64_t n_thre, uint64_t n_a, uint64_t new_idx)
{
    ec_ovec_buf_t *b = NULL;
    b = gen_ec_ovec_buf_t(n_thre);
    if(new_idx) {
        // kt_for(n_thre, worker_hap_dc_ec, b, n_a);///update overlaps
        destroy_cc_v(&scc); destroy_cc_v(&scb); destroy_cc_v(&sca);
        
        ha_print_ovlp_stat_0(b, n_thre, n_a);
    } else {
        ha_print_ovlp_stat_1(b, n_thre, n_a);
        destroy_cc_v(&scc); destroy_cc_v(&scb); destroy_cc_v(&sca);
    }

    destroy_ec_ovec_buf_t(b);
}

void sl_ec_r(uint64_t n_thre, uint64_t n_a)
{
    sl_v *b = NULL; uint64_t k; MALLOC(b, n_thre);
    for (k = 0; k < n_thre; k++) {
        b[k].a = NULL; b[k].n = b[k].m = 0;
        init_UC_Read(&b[k].z); kv_init(b[k].q);
    }

    kt_for(n_thre, worker_sl_ec, b, n_a);///debug_for_fix

    for (k = 0; k < n_thre; k++) {
        free(b[k].a); destory_UC_Read(&b[k].z); kv_destroy(b[k].q);
    }
    free(b);
}

void handle_chemical_r(uint64_t n_thre, uint64_t n_a)
{
    ec_ovec_buf_t *b = NULL; uint64_t k, chem_n = 0, dedup = 0;
    b = gen_ec_ovec_buf_t(n_thre);
    for (k = 0; k < n_thre; ++k) {
        b->a[k].cnt[0] = 0; b->a[k].cnt[1] = 0;
    }

    kt_for(n_thre, worker_hap_dc_ec_chemical_r, b, n_a);

    for (k = 0; k < n_thre; ++k) {
        chem_n += b->a[k].cnt[0]; 
        b->a[k].cnt[0] = 0; b->a[k].cnt[1] = 1;
    }

    kt_for(n_thre, worker_hap_dc_ec_chemical_r, b, n_a);
    
    for (k = 0; k < n_thre; ++k) {
        dedup += b->a[k].cnt[0]; 
        b->a[k].cnt[1] = 2;
    }

    kt_for(n_thre, worker_hap_dc_ec_chemical_r, b, n_a);

    fprintf(stderr, "[M::%s] # chimeric reads: %lu, # arcs:: %lu\n", __func__, chem_n, dedup);

    destroy_ec_ovec_buf_t(b);
}

void handle_chemical_arc(uint64_t n_thre, uint64_t n_a)
{
    ec_ovec_buf_t *b = NULL; uint64_t k, chem_n = 0;
    b = gen_ec_ovec_buf_t(n_thre);
    for (k = 0; k < n_thre; ++k) {
        b->a[k].cnt[0] = 0; b->a[k].cnt[1] = 0;
    }

    kt_for(n_thre, worker_hap_dc_ec_chemical_arc, b, n_a);

    for (k = 0; k < n_thre; ++k) {
        chem_n += b->a[k].cnt[0]; 
        b->a[k].cnt[0] = 0; b->a[k].cnt[1] = 1;
    }

    kt_for(n_thre, worker_hap_dc_ec_chemical_arc, b, n_a);

    fprintf(stderr, "[M::%s] # chimeric reads: %lu\n", __func__, chem_n);

    destroy_ec_ovec_buf_t(b);
}

uint8_t* gen_chemical_arc_rf(uint64_t n_thre, uint64_t n_a)
{
    ec_ovec_buf_t *b = NULL; uint64_t k, chem_n = 0; uint8_t *ra = NULL;
    b = gen_ec_ovec_buf_t(n_thre);
    for (k = 0; k < n_thre; ++k) {
        b->a[k].cnt[0] = 0; b->a[k].cnt[1] = 0;
    }
    MALLOC(ra, n_a); ///memset(ra, -1, sizeof((*ra))*n_a); 
    b->cr = ra;

    kt_for(n_thre, worker_hap_dc_ec_chemical_arc_mark, b, n_a);

    for (k = 0; k < n_thre; ++k) {
        chem_n += b->a[k].cnt[0]; 
        b->a[k].cnt[0] = 0; b->a[k].cnt[1] = 1;
    }

    kt_for(n_thre, worker_hap_dc_ec_chemical_arc_mark, b, n_a);

    fprintf(stderr, "[M::%s] # chimeric reads: %lu\n", __func__, chem_n);

    b->cr = NULL; destroy_ec_ovec_buf_t(b);
    return ra;
}