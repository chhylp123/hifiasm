#ifndef __PURTIG__
#define __PURTIG__
#include <stdio.h>
#include <stdint.h>
#include "kvec.h"
#include "Overlaps.h"
#include "Purge_Dups.h"

typedef struct {
    uint32_t bS, bE;
    uint32_t nS, nE;
    uint32_t uID;
    uint8_t hs;
}pt_node_t;

typedef struct {
	///sid[0]: query id
	///sid[1]: target id
	uint32_t sid[2];
	uint32_t w;
} pt_match1_t;

typedef struct {
    kvec_t(uint64_t) idx;
    kvec_t(pt_match1_t) ma;
    uint64_t* cc;
    uint32_t n_seq;
} pt_match_t;

typedef struct {
	///cnt1: how many unique minimizers
	///cnt2: how many non-unique minimizers
	///uint32_t cnt2, cnt1;
	uint64_t m[2];
	int8_t s;
} pt_uinfo_t;

typedef struct {
    kvec_t(pt_node_t) p;
    ma_ug_t *ug;
    asg_t *rg;
    trans_chain* t_ch;
    ///kvec_t(int8_t) s; ///status
    kvec_t(pt_uinfo_t) info; ///status
    kvec_t(uint32_t) p_idx;
    pt_match_t* e;
}pt_g_t;

void pt_solve(hap_overlaps_list* ovlp, trans_chain* t_ch, ma_ug_t *ug, asg_t *read_g, double f_rate, uint8_t* trio_flag);

#endif