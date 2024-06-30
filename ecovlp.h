#ifndef __ECOVLP_PARSER__
#define __ECOVLP_PARSER__

#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include "Hash_Table.h"
#include "Process_Read.h"

typedef struct {
	uint32_t v:31, f:1;
	uint32_t sc;
} cns_arc;
typedef struct {size_t n, m; cns_arc *a;} cns_arc_v;

typedef struct {
	// uint16_t c:2, t:2, f:1, sc:3;
	uint32_t c:2, f:1, sc:29;
	cns_arc_v in, ou;
}cns_t;

typedef struct {
	size_t n, m; 
	cns_t *a;
	uint32_t si, ei;
}cns_gfa;

typedef struct {
	int is_final, save_ov;
	// chaining and overlapping related buffers
	UC_Read self_read, ovlp_read;
	Candidates_list clist;
	overlap_region_alloc olist;
    overlap_region tmp;
	ha_abuf_t *ab;
	// error correction related buffers
	int64_t num_read_base, num_correct_base, num_recorrect_base;
	Cigar_record cigar;
	Correct_dumy correct;
	haplotype_evdience_alloc hap;
	bit_extz_t exz;
	// asg32_v v32;
	kv_ul_ov_t pidx;
	asg64_v v64;
	asg32_v v32;
	asg16_v v16;

    kvec_t_u64_warp r_buf;
    kvec_t_u8_warp k_flag;
	st_mt_t sp;
	cns_gfa cns;
	
} ec_ovec_buf_t0;

typedef struct {
	ec_ovec_buf_t0 *a;
	uint32_t n;
} ec_ovec_buf_t;

ec_ovec_buf_t* gen_ec_ovec_buf_t(uint32_t n, uint32_t is_final, uint32_t save_ov);
void destroy_ec_ovec_buf_t(ec_ovec_buf_t *p);
void cal_ec_multiple(ec_ovec_buf_t *b, uint64_t n_thre, uint64_t n_a);
void prt_chain(overlap_region_alloc *o);

#endif