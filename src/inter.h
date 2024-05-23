#ifndef __INTER__
#define __INTER__
#include "Overlaps.h"
#include "Process_Read.h"
#include "hic.h"

#define G_CHAIN_BW 16//128
#define FLANK_M (0x7fffU)
#define P_CHAIN_COV 0.985
#define P_FRAGEMENT_CHAIN_COV 0.20
#define P_FRAGEMENT_PRIMARY_CHAIN_COV 0.70
#define P_FRAGEMENT_PRIMARY_SECOND_COV 0.25
#define P_CHAIN_SCORE 0.6
#define G_CHAIN_GAP 0.1
#define UG_SKIP 5
#define RG_SKIP 25
#define UG_SKIP_GRAPH_N 72
#define UG_SKIP_N 100
#define UG_ITER_N 5000
#define UG_DIS_N 50000
// #define UG_TRANS_W 2
#define UG_TRANS_W 2
// #define UG_TRANS_ERR_W 512
#define UG_TRANS_ERR_W 64
#define G_CHAIN_TRANS_RATE 0.25
#define G_CHAIN_TRANS_WEIGHT -1
#define G_CHAIN_INDEL 128
#define W_CHN_PEN_GAP 0.1
#define N_GCHAIN_RATE 0.04
#define PRIMARY_UL_CHAIN_MIN 75000

typedef struct {
	int w, k, bw, max_gap, is_HPC, hap_n, occ_weight, max_gap_pre, max_gc_seq_ext, seed;
	int max_lc_skip, max_lc_iter, min_lc_cnt, min_lc_score, max_gc_skip, ref_bonus;
	int min_gc_cnt, min_gc_score, sub_diff, best_n;
    float chn_pen_gap, mask_level, pri_ratio;
	///base-alignment
	double bw_thres, diff_ec_ul, diff_ec_ul_low, diff_ec_ul_hpc; int max_n_chain, ec_ul_round;
} mg_idxopt_t;

struct mg_tbuf_s {
	void *km;
	int frag_gap;
};
typedef struct mg_tbuf_s mg_tbuf_t;


mg_tbuf_t *mg_tbuf_init(void);

void mg_tbuf_destroy(mg_tbuf_t *b);

void *mg_tbuf_get_km(mg_tbuf_t *b);

typedef struct {
	FILE *fp;
	ul_vec_t u;
	uint64_t flag;
} ucr_file_t;

typedef struct {
	int32_t off, cnt;
	uint32_t v;
	int32_t score;
} mg_llchain_t;

typedef struct {
	int32_t id, parent;
	int32_t off, cnt;
	int32_t n_anchor, score;
	int32_t qs, qe;
	int32_t plen, ps, pe;
	int32_t blen, mlen;
	float div;
	uint32_t hash;
	int32_t subsc, n_sub;
	uint32_t mapq:8, flt:1, dummy:23;
} mg_gchain_t;

typedef struct {
	size_t n,m;
    uint64_t *a, tl;
	kvec_t(char) cc;
} mg_dbn_t;

typedef struct {
	int32_t cnt;
	uint32_t v;
	int32_t score;
	uint32_t qs, qe, ts, te;
} mg_lres_t;

typedef struct {
    int32_t n_gc, n_lc;
    mg_gchain_t *gc;///g_chain; idx in l_chains
    mg_lres_t *lc;///l_chain
    uint64_t qid, qlen;
} mg_gres_t;

typedef struct {
	size_t n,m;
    mg_gres_t *a;
	uint64_t total_base;
    uint64_t total_pair;
} mg_gres_a;

void push_uc_block_t(const ug_opt_t *uopt, kv_ul_ov_t *z, char **seq, uint64_t *len, uint64_t b_id);
void ul_resolve(ma_ug_t *ug, const asg_t *rg, const ug_opt_t *uopt, int hap_n);
void ul_load(const ug_opt_t *uopt);
uint64_t* get_hifi2ul_list(all_ul_t *x, uint64_t hid, uint64_t* a_n);
uint64_t ul_refine_alignment(const ug_opt_t *uopt, asg_t *sg);
ma_ug_t *ul_realignment(const ug_opt_t *uopt, asg_t *sg, uint32_t double_check_cache, const char *bin_file);
int32_t write_all_ul_t(all_ul_t *x, char* file_name, ma_ug_t *ug);
int32_t load_all_ul_t(all_ul_t *x, char* file_name, All_reads *hR, ma_ug_t *ug);
uint32_t ugl_cover_check(uint64_t is, uint64_t ie, ma_utg_t *u);
void filter_ul_ug(ma_ug_t *ug);
void gen_ul_vec_rid_t(all_ul_t *x, All_reads *rdb, ma_ug_t *ug);
void update_ug_arch_ul_mul(ma_ug_t *ug);
void print_ul_alignment(ma_ug_t *ug, all_ul_t *aln, uint32_t id, const char* cmd);
void clear_all_ul_t(all_ul_t *x);
void trans_base_infer(ma_ug_t *ug, asg_t *sg, ug_opt_t *uopt, kv_u_trans_t *res, bubble_type *bub);
hpc_re_t *gen_hpc_re_t(ma_ug_t *ug);
idx_emask_t* graph_ovlp_binning(ma_ug_t *ug, asg_t *sg, const ug_opt_t *uopt);
uint32_t gen_src_shared_interval_simple(uint32_t src, ma_ug_t *ug, uint64_t *flt, uint64_t flt_n, kv_ul_ov_t *res);
uint64_t check_ul_ov_t_consist(ul_ov_t *x, ul_ov_t *y, int64_t ql, int64_t tl, double diff);
uint32_t infer_se(uint32_t qs, uint32_t qe, uint32_t ts, uint32_t te, uint32_t rev, 
uint32_t rqs, uint32_t rqe, uint32_t *rts, uint32_t *rte);
uint32_t clean_contain_g(const ug_opt_t *uopt, asg_t *sg, uint32_t push_trans);
void dedup_contain_g(const ug_opt_t *uopt, asg_t *sg);
void trans_base_mmhap_infer(ma_ug_t *ug, asg_t *sg, ug_opt_t *uopt, kv_u_trans_t *res);
scaf_res_t *gen_contig_path(const ug_opt_t *uopt, asg_t *sg, ma_ug_t *ctg, ma_ug_t *ref);
void gen_contig_trans(const ug_opt_t *uopt, asg_t *sg, ma_ug_t *qry, scaf_res_t *qry_sc, ma_ug_t *ref, scaf_res_t *ref_sc, ma_ug_t *gfa, kv_u_trans_t *ta, uint32_t qoff, uint32_t toff, bubble_type *bu, kv_u_trans_t *res);
void gen_contig_self(const ug_opt_t *uopt, asg_t *sg, ma_ug_t *db, scaf_res_t *db_sc, ma_ug_t *gfa, kv_u_trans_t *ta, uint64_t soff, bubble_type *bu, kv_u_trans_t *res, uint32_t is_exact);
void order_contig_trans(kv_u_trans_t *in);
void sort_uc_block_qe(uc_block_t* a, uint64_t a_n);

#endif
