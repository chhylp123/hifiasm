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
// #define w_t int32_t
// #define t_w_t int64_t
// #define w_cast(x) ((t_w_t)((x) < 0 ? (x) - 0.5 : (x) + 0.5))

#define w_t double
#define t_w_t double
#define w_cast(x) ((t_w_t)((x)))


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
void mc_solve(hap_overlaps_list* ovlp, trans_chain* t_ch, kv_u_trans_t *ta, ma_ug_t *ug, asg_t *read_g, double f_rate, uint8_t* trio_flag, uint32_t renew_s, int8_t *s, uint32_t is_sys, bubble_type* bub);
#endif