#ifndef __TOVLP__
#define __TOVLP__
#include <stdint.h>
#include "Overlaps.h"

utg_trans_t *init_utg_trans_t(ma_ug_t *ug, ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, R_to_U* ruIndex, asg_t *read_g, int max_hang, int min_ovlp);
void destroy_utg_trans_t(utg_trans_t **o);
void asg_bub_collect_ovlp(ma_ug_t *ug, uint32_t v0, buf_t *b, utg_trans_t *o);
void collect_trans_ovlp(const char* cmd, buf_t* pri, uint64_t pri_offset, buf_t* aux, uint64_t aux_offset,
ma_ug_t *ug, utg_trans_t *o);
int asg_arc_decompress(asg_t *g, ma_ug_t *ug, asg_t *read_sg, ma_hit_t_alloc* reverse_sources,
R_to_U* ruIndex, utg_trans_t *o);
int asg_arc_decompress_mul(asg_t *g, ma_ug_t *ug, asg_t *read_sg, uint32_t positive_flag, uint32_t negative_flag,
ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex, utg_trans_t *o);
#endif
