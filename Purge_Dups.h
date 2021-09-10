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
#define REAL_ALTER_THRES 0.25
#define CHAIN_FILTER_RATE 0.7
#define REV_W 8

#define SELF_EXIST 0
#define REVE_EXIST 1
#define DELETE 2
#define MIXED 3
#define FLIP 4

#define X2Y 0
#define Y2X 1
#define XCY 2
#define YCX 3

#define Cal_Off(OFF) ((long long)((uint32_t)((OFF)>>32)) - (long long)((uint32_t)((OFF))))
#define Get_xOff(OFF) ((long long)((uint32_t)((OFF)>>32)))
#define Get_yOff(OFF) ((long long)((uint32_t)((OFF))))
#define Get_match(x) ((x).weight)
#define Get_total(x) ((x).index_beg)
#define Get_type(x) ((x).index_end)
#define Get_x_beg(x) ((x).x_beg_pos)
#define Get_x_end(x) ((x).x_end_pos)
#define Get_y_beg(x) ((x).y_beg_pos)
#define Get_y_end(x) ((x).y_end_pos)
#define Get_rev(x) ((x).rev)

typedef struct {
    uint8_t rev;
    uint8_t type;
    uint8_t status;
    uint32_t x_beg_pos;
    uint32_t x_end_pos;
    uint32_t y_beg_pos;
    uint32_t y_end_pos;
    uint32_t x_beg_id;
    uint32_t x_end_id;
    uint32_t y_beg_id;
    uint32_t y_end_id;
    uint32_t xUid;
    uint32_t yUid;
    uint32_t weight;
    long long score;
    float s;
}hap_overlaps;

typedef struct {
    kvec_t(hap_overlaps) a;
}kvec_hap_overlaps;

typedef struct {
    kvec_hap_overlaps* x;
    uint32_t num;
}hap_overlaps_list;

void purge_dups(ma_ug_t *ug, asg_t *read_g, ma_sub_t* coverage_cut, ma_hit_t_alloc* sources, 
ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex, kvec_asg_arc_t_warp* edge, float density, 
uint32_t purege_minLen, int max_hang, int min_ovlp, float drop_ratio, uint32_t just_contain, 
uint32_t just_coverage, hap_cov_t *cov, uint32_t collect_p_trans, uint32_t collect_p_trans_f);
void fill_unitig(uint64_t* buffer, uint32_t bufferLen, asg_t* read_g, kvec_asg_arc_t_warp* edge,
uint32_t is_circle, uint64_t* rLen);
void get_contig_length(ma_ug_t *ug, asg_t *g, uint64_t* primaryLen, uint64_t* alterLen);
void enable_debug_mode(uint32_t mode);
hap_cov_t* init_hap_cov_t(ma_ug_t *ug, asg_t* read_g, ma_hit_t_alloc* sources, R_to_U* ruIndex, ma_hit_t_alloc* reverse_sources, 
ma_sub_t *coverage_cut, int max_hang, int min_ovlp, uint32_t is_collect_trans);
void destory_hap_cov_t(hap_cov_t **x);
void chain_trans_ovlp(hap_cov_t *cov, utg_trans_t *o, ma_ug_t *ug, asg_t *read_sg, buf_t* xReads, uint32_t targetBaseLen, uint32_t* xEnd);
int get_specific_hap_overlap(kvec_hap_overlaps* x, uint32_t qn, uint32_t tn);
void set_reverse_hap_overlap(hap_overlaps* dest, hap_overlaps* source, uint32_t* types);
void print_hap_paf(ma_ug_t *ug, hap_overlaps* ovlp);
uint64_t get_xy_pos_by_pos(asg_t *read_g, asg_arc_t* t, uint32_t v_in_unitig, uint32_t w_in_unitig, 
uint32_t v_in_pos, uint32_t w_in_pos, uint32_t xUnitigLen, uint32_t yUnitigLen, uint8_t* rev);
void quick_LIS(asg_arc_t_offset* x, uint32_t n, kvec_t_i32_warp* tailIndex, kvec_t_i32_warp* prevIndex);
uint32_t classify_hap_overlap(long long xBeg, long long xEnd, long long xLen,
long long yBeg, long long yEnd, long long yLen, long long* r_xBeg, long long* r_xEnd, 
long long* r_yBeg, long long* r_yEnd);
int cmp_hap_alignment_chaining(const void * a, const void * b);
uint32_t classify_hap_overlap(long long xBeg, long long xEnd, long long xLen,
long long yBeg, long long yEnd, long long yLen, long long* r_xBeg, long long* r_xEnd, 
long long* r_yBeg, long long* r_yEnd);

#endif