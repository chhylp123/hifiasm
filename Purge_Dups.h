#ifndef __PURGEDUPS__
#define __PURGEDUPS__
#include <stdio.h>
#include <stdint.h>
#include "kvec.h"
#include "kdq.h"
#include "Overlaps.h"
#include "Hash_Table.h"
#define COV_COUNT 1024
#define HOM_PEAK_RATE 1.25
#define HET_PEAK_RATE (HOM_PEAK_RATE*2)
#define ALTER_COV_THRES 0.9
#define REAL_ALTER_THRES 0.1
#define CHAIN_FILTER_RATE 0.7

void purge_dups(ma_ug_t *ug, asg_t *read_g, ma_sub_t* coverage_cut, ma_hit_t_alloc* sources, 
ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex, kvec_asg_arc_t_warp* edge, float density, 
uint32_t purege_minLen, int max_hang, int min_ovlp, long long bubble_dist, float drop_ratio, 
uint32_t just_contain, uint32_t just_coverage, hc_links* link, hap_cov_t *cov);
void fill_unitig(uint64_t* buffer, uint32_t bufferLen, asg_t* read_g, kvec_asg_arc_t_warp* edge,
uint32_t is_circle, uint64_t* rLen);
void get_contig_length(ma_ug_t *ug, asg_t *g, uint64_t* primaryLen, uint64_t* alterLen);
void enable_debug_mode(uint32_t mode);
hap_cov_t* init_hap_cov_t(ma_ug_t *ug, asg_t* read_g, ma_hit_t_alloc* sources, R_to_U* ruIndex, 
ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, int max_hang, int min_ovlp);
void destory_hap_cov_t(hap_cov_t **x);
void chain_trans_ovlp(hap_cov_t *cov, ma_ug_t *ug, asg_t *read_sg, buf_t* xReads, uint32_t targetBaseLen, uint32_t* xEnd);


#endif