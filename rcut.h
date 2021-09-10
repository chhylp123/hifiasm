#ifndef __RCUT__
#define __RCUT__
#include <stdio.h>
#include <stdint.h>
#include "kvec.h"
#include "Overlaps.h"
#include "Purge_Dups.h"
#include "hic.h"

typedef struct {
    uint32_t bS, bE;
    uint32_t nS, nE;
    uint32_t uID;
    uint8_t hs;
}mc_interval_t;

#define mc_node_t int8_t
#define mcg_node_t uint32_t
// #define w_t int64_t
// #define t_w_t int64_t
// #define w_cast(x) ((t_w_t)((x) < 0 ? (x) - 0.5 : (x) + 0.5))

#define w_t double
#define t_w_t double
#define w_cast(x) ((t_w_t)((x)))
#define MC_NAME "debug_mc.bin"

typedef struct {
	uint64_t x; ///(uint64_t)nid1 << 32 | nid2;
	w_t w;  ///might be negative or positive
} mc_edge_t;

typedef struct {
    kvec_t(uint64_t) idx;
    kvec_t(mc_edge_t) ma;
    uint64_t* cc;
    uint32_t n_seq;
} mc_match_t;

typedef struct {
    kvec_t(mc_node_t) s;
    ma_ug_t *ug;
    asg_t *rg;
    mc_match_t* e;
}mc_g_t;

typedef struct {
    uint32_t a[2], occ[2]; 
    mc_node_t s[2];
    t_w_t z[4];
}mb_node_t;

typedef struct {
    kvec_t(uint32_t) bid;
    kvec_t(uint32_t) idx;
    kvec_t(mb_node_t) u;
}mb_nodes_t;

typedef struct {
	uint64_t x; ///(uint64_t)nid1 << 32 | nid2;
	t_w_t w[4];  ///might be negative or positive
} mb_edge_t;

typedef struct {
    kvec_t(uint64_t) idx;
    kvec_t(mb_edge_t) ma;
    uint64_t* cc;
    uint32_t n_seq;
} mb_match_t;

typedef struct {
    mb_nodes_t* u; 
    mb_match_t* e;
}mb_g_t;

typedef struct {
	mcg_node_t s;
	uint16_t h[2], hc;
	double hw[2];
}mc_gg_status;

typedef struct {
	mc_gg_status *a;
	size_t n, m;
}kv_gg_status;

typedef struct {
	mcg_node_t *a;
	size_t n, m;
}mcb_t;


typedef struct {
	kv_gg_status *s;
    // ma_ug_t *ug;
    // asg_t *rg;
    uint32_t un;
    mc_match_t* e;
	kvec_t(mcb_t) m;
	mcg_node_t mask;
	uint16_t hN;
}mc_gg_t;


static inline uint64_t kr_splitmix64(uint64_t x)
{
	uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
	return z ^ (z >> 31);
}

static inline double kr_drand_r(uint64_t *x)
{
    union { uint64_t i; double d; } u;
	*x = kr_splitmix64(*x);
    u.i = 0x3FFULL << 52 | (*x) >> 12;
    return u.d - 1.0;
}

void mc_solve(hap_overlaps_list* ovlp, trans_chain* t_ch, kv_u_trans_t *ta, ma_ug_t *ug, asg_t *read_g, double f_rate, uint8_t* trio_flag, uint32_t renew_s, int8_t *s, uint32_t is_sys, bubble_type* bub, kv_u_trans_t *ref, int clean_ov);
void debug_mc_g_t(const char* name);
void mc_solve_general(kv_u_trans_t *ta, uint32_t un, kv_gg_status *s, uint16_t hapN, uint16_t update_ta, uint16_t write_dump);
kv_gg_status *init_mc_gg_status(ma_ug_t *ug, asg_t *read_g, ma_sub_t* coverage_cut, 
ma_hit_t_alloc* sources, R_to_U* ruIndex, uint64_t t_cov, uint16_t hapN);
void destory_mc_gg_t(mc_gg_t **p);
void debug_mc_gg_t(const char* fn, uint32_t update_ta, uint32_t convert_mc_g_t);
#endif