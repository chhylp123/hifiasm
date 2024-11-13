#ifndef __ECOVLP_PARSER__
#define __ECOVLP_PARSER__

#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include "Hash_Table.h"
#include "Process_Read.h"
#include "kdq.h"

KDQ_INIT(uint32_t)

typedef struct {
	uint32_t v:31, f:1;
	uint32_t sc;
} cns_arc;
typedef struct {size_t n, m, nou; cns_arc *a; } cns_arc_v;

typedef struct {
	// uint16_t c:2, t:2, f:1, sc:3;
	uint32_t c:2, f:1, sc:29;
	cns_arc_v arc;
}cns_t;

typedef struct {
	size_t n, m; 
	cns_t *a;
	uint32_t si, ei, off, bn, bb0, bb1, cns_g_wl;
	kdq_t(uint32_t) *q;
}cns_gfa;

typedef struct {
	// chaining and overlapping related buffers
	UC_Read self_read, ovlp_read;
	Candidates_list clist;
	overlap_region_alloc olist;
	ha_abuf_t *ab;
	// int64_t num_read_base, num_correct_base, num_recorrect_base;
	uint64_t cnt[6], rr;
	haplotype_evdience_alloc hap;
	bit_extz_t exz;
	kv_ul_ov_t pidx;
	asg64_v v64;
	asg32_v v32;
	asg16_v v16;
	asg8_v v8q, v8t;

    kvec_t_u8_warp k_flag;
	st_mt_t sp;
	cns_gfa cns;	
} ec_ovec_buf_t0;

typedef struct {
	ec_ovec_buf_t0 *a;
	uint32_t n, rev;
} ec_ovec_buf_t;

ec_ovec_buf_t* gen_ec_ovec_buf_t(uint32_t n);
void destroy_ec_ovec_buf_t(ec_ovec_buf_t *p);
void prt_chain(overlap_region_alloc *o);
void cal_ec_r(uint64_t n_thre, uint64_t round, uint64_t n_round, uint64_t n_a, uint64_t is_sv, uint64_t *tot_b, uint64_t *tot_e);
void sl_ec_r(uint64_t n_thre, uint64_t n_a);
void cal_ov_r(uint64_t n_thre, uint64_t n_a, uint64_t new_idx);

#endif