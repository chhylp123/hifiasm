#ifndef __COMMAND_LINE_PARSER__
#define __COMMAND_LINE_PARSER__

#define __STDC_LIMIT_MACROS
#include <pthread.h>
#include <stdint.h>

#define HA_VERSION "0.19.9-r616"

#define VERBOSE 0

#define HA_F_NO_HPC          0x1
#define HA_F_NO_KMER_FLT     0x2
#define HA_F_VERBOSE_GFA     0x4
#define HA_F_WRITE_EC        0x8
#define HA_F_WRITE_PAF       0x10
#define HA_F_SKIP_TRIOBIN    0x20
#define HA_F_PURGE_CONTAIN   0x40
#define HA_F_PURGE_JOIN      0x80
#define HA_F_BAN_POST_JOIN   0x100
#define HA_F_BAN_ASSEMBLY    0x200
#define HA_F_HIGH_HET        0x400
#define HA_F_PARTITION       0x800
#define HA_F_FAST            0x1000
#define HA_F_USKEW           0x2000

#define HA_MIN_OV_DIFF       0.02 // min sequence divergence in an overlap
#define MIN_N_CHAIN          100

typedef struct{
    int *l, n; 
    char **a;
}enzyme;

typedef struct {
	int flag;
    int num_reads;
    char** read_file_names;
    char* output_file_name;
    char* required_read_name;
	char *fn_bin_yak[2];
	char *fn_bin_list[2];
    char *fn_bin_poy;
	char *extract_list;
    enzyme *hic_reads[2];
    enzyme *hic_enzymes;
    enzyme *ar;
    enzyme *sec_in;
	int extract_iter;
    int thread_num;
    int k_mer_length;
    int hic_mer_length;
    int ul_mer_length;
    int trans_mer_length;
    int bub_mer_length;
	int mz_win;
    int ul_mz_win;
    int trans_win;
    int mz_rewin;
    int ul_mz_rewin;
	int mz_sample_dist;
	int bf_shift;
	int max_kmer_cnt;
	double high_factor; // coverage cutoff set to high_factor*hom_cov
	double max_ov_diff_ec;
	double max_ov_diff_final;
	int hom_cov;
    int het_cov;
    int b_low_cov;
    int b_high_cov;
    double m_rate;
	int max_n_chain; // fall-back max number of chains to consider
	int min_hist_kmer_cnt;
    int load_index_from_disk;
    int write_index_to_disk;
    int number_of_round;
    int adapterLen;
    int clean_round;
    int roundID;
    int max_hang_Len;
    int gap_fuzz;
    int min_overlap_Len;
    int min_overlap_coverage;
    int max_short_tip;
    int max_short_ul_tip;
    int max_contig_tip;
    int min_cnt;
    int mid_cnt;
    int purge_level_primary;
    int purge_level_trio;
    int purge_overlap_len;
    ///int purge_overlap_len_hic;
    int recover_atg_cov_min;
    int recover_atg_cov_max;
    int hom_global_coverage;
    int hom_global_coverage_set;
    int pur_global_coverage;
    int bed_inconsist_rate;
    int hic_inconsist_rate;

    float max_hang_rate;
    float min_drop_rate;
    float max_drop_rate;
    float purge_simi_rate_l2;
    float purge_simi_rate_l3;
    float purge_simi_thres;
    float trans_base_rate;
    float trans_base_rate_sec;
    float min_path_drop_rate;
    float max_path_drop_rate;
    // uint64_t path_clean_round;

    ///float purge_simi_rate_hic;

    long long small_pop_bubble_size;
    long long large_pop_bubble_size;
    long long num_bases;
    long long num_corrected_bases;
    long long num_recorrected_bases;
	long long mem_buf;
    long long coverage;
    int hap_occ;
    int polyploidy;
    int trio_flag_occ_thres;
    uint64_t seed;
    int32_t n_perturb;
    double f_perturb;
    int32_t n_weight;
    uint32_t is_alt;
    uint64_t misjoin_len;
    uint64_t scffold;
    int32_t dp_min_len;
    float dp_e;
    int64_t hg_size;
    float kpt_rate;
    int64_t infor_cov, s_hap_cov, trio_cov_het_ovlp;
    double ul_error_rate, ul_error_rate_low, ul_error_rate_hpc;
    int32_t ul_ec_round;
    uint8_t is_dbg_het_cnt;
    uint8_t is_low_het_ul;
    uint8_t is_base_trans;
    uint8_t is_read_trans;
    uint8_t is_topo_trans;
    uint8_t is_bub_trans;
    uint8_t bin_only;
    int32_t ul_clean_round;
    int32_t prt_dbg_gfa;
    int32_t integer_correct_round;
    uint8_t dbg_ovec_cal;
    uint8_t hifi_pst_join, ul_pst_join;
    uint32_t ul_min_base;
    uint8_t self_scaf;
    uint64_t self_scaf_min;
    uint64_t self_scaf_reliable_min;
    int64_t self_scaf_gap_max;
    int64_t somatic_cov;

    char *telo_motif;
    int64_t telo_pen;
    int64_t telo_drop;
    int64_t telo_mic_sc;
} hifiasm_opt_t;

extern hifiasm_opt_t asm_opt;

void init_opt(hifiasm_opt_t* asm_opt);
void destory_opt(hifiasm_opt_t* asm_opt);
void ha_opt_reset_to_round(hifiasm_opt_t* asm_opt, int round);
void ha_opt_update_cov(hifiasm_opt_t *opt, int hom_cov);
void ha_opt_update_cov_min(hifiasm_opt_t *opt, int hom_cov, int min_chain);
int CommandLine_process(int argc, char *argv[], hifiasm_opt_t* asm_opt);
double Get_T(void);

static inline int ha_opt_triobin(const hifiasm_opt_t *opt)
{
	return ((opt->fn_bin_yak[0] && opt->fn_bin_yak[1]) || (opt->fn_bin_list[0] && opt->fn_bin_list[1]));
}

static inline int ha_opt_hic(const hifiasm_opt_t *opt)
{
    return ((opt->hic_reads[0] && opt->hic_reads[1]));
}

#endif
