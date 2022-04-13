#ifndef __GFA_UT__
#define __GFA_UT__
#include "Overlaps.h"

void ul_clean_gfa(ug_opt_t *uopt, asg_t *sg, ma_hit_t_alloc *src, ma_hit_t_alloc *rev, R_to_U* rI, int64_t clean_round, double min_ovlp_drop_ratio, double max_ovlp_drop_ratio, 
double ou_drop_rate, int64_t max_tip, bub_label_t *b_mask_t, int32_t is_ou, int32_t is_trio, uint32_t ou_thres);
uint32_t asg_arc_cut_tips(asg_t *g, uint32_t max_ext, asg64_v *in, uint32_t is_ou);
void asg_iterative_semi_circ(asg_t *g, ma_hit_t_alloc* src, asg64_v *in, uint32_t normal_len, uint32_t pop_chimer, asg64_v *dbg);
void asg_arc_cut_chimeric(asg_t *g, ma_hit_t_alloc* src, asg64_v *in, uint32_t ou_thres);
void asg_arc_cut_inexact(asg_t *g, ma_hit_t_alloc* src, asg64_v *in, int32_t max_ext, uint32_t is_ou, uint32_t is_trio);
void asg_arc_cut_length(asg_t *g, asg64_v *in, int32_t max_ext, float len_rat, float ou_rat, uint32_t is_ou, uint32_t is_trio, 
uint32_t is_topo, ma_hit_t_alloc *rev, R_to_U* rI, uint32_t *max_drop_len);
void asg_arc_cut_bub_links(asg_t *g, asg64_v *in, float len_rat, float sec_len_rat, float ou_rat, uint32_t is_ou, uint64_t check_dist, ma_hit_t_alloc *rev, R_to_U* rI, int32_t max_ext);
void asg_arc_cut_complex_bub_links(asg_t *g, asg64_v *in, float len_rat, float ou_rat, uint32_t is_ou, bub_label_t *b_mask_t);
uint32_t asg_cut_large_indel(asg_t *g, asg64_v *in, int32_t max_ext, float ou_rat, uint32_t is_ou);
uint32_t asg_cut_semi_circ(asg_t *g, uint32_t lim_len, uint32_t is_clean);

#endif
