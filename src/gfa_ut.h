#ifndef __GFA_UT__
#define __GFA_UT__
#include "Overlaps.h"
#include "hic.h"

#define is_contain_r(ri, z) (((z)<(ri).len)&&((ri).index[(z)]!=(uint32_t)(-1))&&(!((ri).index[(z)]>>31)))

typedef struct {
	asg_t *g; 
	ma_hit_t_alloc *src;
    ma_sub_t *cov;
    R_to_U* ruIndex;
    int64_t max_hang; 
    int64_t min_ovlp; 
    int64_t ul_occ;
} sset_aux;

void ul_clean_gfa(ug_opt_t *uopt, asg_t *sg, ma_hit_t_alloc *src, ma_hit_t_alloc *rev, R_to_U* rI, int64_t clean_round, double min_ovlp_drop_ratio, double max_ovlp_drop_ratio, 
double ou_drop_rate, int64_t max_tip, int64_t gap_fuzz, bub_label_t *b_mask_t, int32_t is_ou, int32_t is_trio, uint32_t ou_thres, char *o_file);
uint32_t asg_arc_cut_tips(asg_t *g, uint32_t max_ext, asg64_v *in, uint32_t is_ou, R_to_U *ru, telo_end_t *te);
void asg_iterative_semi_circ(asg_t *g, ma_hit_t_alloc* src, asg64_v *in, uint32_t normal_len, uint32_t pop_chimer, asg64_v *dbg, telo_end_t *te);
void asg_arc_cut_chimeric(asg_t *g, ma_hit_t_alloc* src, asg64_v *in, uint32_t ou_thres, telo_end_t *te);
void asg_arc_cut_inexact(asg_t *g, ma_hit_t_alloc* src, asg64_v *in, int32_t max_ext, uint32_t is_ou, uint32_t is_trio, uint32_t min_diff, float ou_rat/**, asg64_v *dbg**/);
void asg_arc_cut_length(asg_t *g, asg64_v *in, int32_t max_ext, float len_rat, float ou_rat, uint32_t is_ou, uint32_t is_trio, 
uint32_t is_topo, uint32_t min_diff, uint32_t min_ou, ma_hit_t_alloc *rev, R_to_U* rI, uint32_t *max_drop_len);
void asg_arc_cut_bub_links(asg_t *g, asg64_v *in, float len_rat, float sec_len_rat, float ou_rat, uint32_t is_ou, uint64_t check_dist, ma_hit_t_alloc *rev, R_to_U* rI, int32_t max_ext);
void asg_arc_cut_complex_bub_links(asg_t *g, asg64_v *in, float len_rat, float ou_rat, uint32_t is_ou, bub_label_t *b_mask_t);
uint32_t asg_cut_large_indel(asg_t *g, asg64_v *in, int32_t max_ext, float ou_rat, uint32_t is_ou, uint32_t min_diff);
uint32_t asg_cut_semi_circ(asg_t *g, uint32_t lim_len, uint32_t is_clean);
void ul_realignment_gfa(ug_opt_t *uopt, asg_t *sg, int64_t clean_round, double min_ovlp_drop_ratio, 
double max_ovlp_drop_ratio, int64_t max_tip, int64_t max_ul_tip, bub_label_t *b_mask_t, uint32_t is_trio, char *o_file, ul_renew_t *ropt, const char *bin_file, uint64_t free_uld, uint64_t is_bridg, uint64_t deep_clean);
void recover_contain_g(asg_t *g, ma_hit_t_alloc *src, R_to_U* ruIndex, int64_t max_hang, int64_t min_ovlp, int64_t ul_occ);
void normalize_gou(asg_t *g);
void prt_specfic_sge(asg_t *g, uint32_t src, uint32_t dst, const char* cmd);
asg_t *gen_ng(ma_ug_t *ug, asg_t *sg, ug_opt_t *uopt, ma_sub_t **cov, R_to_U *ruI, uint64_t scaffold_len);
void post_rescue(ug_opt_t *uopt, asg_t *sg, ma_hit_t_alloc *src, ma_hit_t_alloc *rev, R_to_U* rI, bub_label_t *b_mask_t, long long no_trio_recover);
// void print_raw_u2rgfa_seq(all_ul_t *aln, R_to_U* rI, uint32_t is_detail);
bubble_type *gen_bubble_chain(asg_t *sg, ma_ug_t *ug, ug_opt_t *uopt, uint8_t **ir_het, uint8_t avoid_het);
void filter_sg_by_ug(asg_t *rg, ma_ug_t *ug, ug_opt_t *uopt);
void ug_ext_gfa(ug_opt_t *uopt, asg_t *sg, uint32_t max_len);
void update_sg_uo(asg_t *g, ma_hit_t_alloc *src);
uint32_t get_arcs(asg_t *g, uint32_t v, uint32_t* idx, uint32_t idx_n);
uint64_t ug_occ_w(uint64_t is, uint64_t ie, ma_utg_t *u);

#endif
