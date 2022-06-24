#ifndef __INTER__
#define __INTER__
#include "Overlaps.h"
#include "Process_Read.h"

void ul_resolve(ma_ug_t *ug, const asg_t *rg, const ug_opt_t *uopt, int hap_n);
void ul_load(const ug_opt_t *uopt);
uint64_t* get_hifi2ul_list(all_ul_t *x, uint64_t hid, uint64_t* a_n);
uint64_t ul_refine_alignment(const ug_opt_t *uopt, asg_t *sg);
ma_ug_t *ul_realignment(const ug_opt_t *uopt, asg_t *sg);
int32_t write_all_ul_t(all_ul_t *x, char* file_name, ma_ug_t *ug);
int32_t load_all_ul_t(all_ul_t *x, char* file_name, All_reads *hR, ma_ug_t *ug);
uint32_t ugl_cover_check(uint64_t is, uint64_t ie, ma_utg_t *u);
void filter_ul_ug(ma_ug_t *ug);
void gen_ul_vec_rid_t(all_ul_t *x, All_reads *rdb, ma_ug_t *ug);
void update_ug_arch_ul_mul(ma_ug_t *ug);
void print_ul_alignment(ma_ug_t *ug, all_ul_t *aln, uint32_t id, const char* cmd);

#endif
