#ifndef __ECOVLP_PARSER__
#define __ECOVLP_PARSER__

#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include "Hash_Table.h"
#include "Process_Read.h"
#include "kdq.h"


void prt_chain(overlap_region_alloc *o);
void cal_ec_r(uint64_t n_thre, uint64_t round, uint64_t n_round, uint64_t n_a, uint64_t is_sv, uint64_t *tot_b, uint64_t *tot_e);
void sl_ec_r(uint64_t n_thre, uint64_t n_a);
void cal_ov_r(uint64_t n_thre, uint64_t n_a, uint64_t new_idx);
void handle_chemical_r(uint64_t n_thre, uint64_t n_a);
void handle_chemical_arc(uint64_t n_thre, uint64_t n_a);
uint8_t* gen_chemical_arc_rf(uint64_t n_thre, uint64_t n_a);
void cal_ec_r_dbg(uint64_t n_thre, uint64_t n_a);

#endif