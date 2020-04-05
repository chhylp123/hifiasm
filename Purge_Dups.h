#ifndef __PURGEDUPS__
#define __PURGEDUPS__
#include <stdio.h>
#include <stdint.h>
#include "kvec.h"
#include "kdq.h"
#include "Overlaps.h"
#include "Hash_Table.h"

void purge_dups(ma_ug_t *ug, asg_t *read_g, ma_sub_t* coverage_cut, ma_hit_t_alloc* reverse_sources, 
R_to_U* ruIndex, kvec_asg_arc_t_warp* edge, float density, uint32_t bi_graph_Len, uint32_t long_hap_overlap, 
float lable_match_rate, int max_hang, int min_ovlp, long long bubble_dist, float drop_ratio,
uint32_t just_contain);
void fill_unitig(uint64_t* buffer, uint32_t bufferLen, asg_t* read_g, kvec_asg_arc_t_warp* edge,
uint32_t is_circle, uint64_t* rLen);

void enable_debug_mode();

#endif