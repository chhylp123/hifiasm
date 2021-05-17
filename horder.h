#ifndef __HORDER__
#define __HORDER__
#include <stdint.h>
#include "hic.h"

typedef struct {
	uint32_t v;
    uint32_t u;
	uint32_t occ:31, del:1;
	double w, nw;
} osg_arc_t;

typedef struct {
	double mw[2], ez[2];
    uint8_t del;
} osg_seq_t;

typedef struct {
	uint32_t m_arc, n_arc:31, is_srt:1;
	osg_arc_t *arc;

	uint32_t m_seq, n_seq:31, is_symm:1;
	osg_seq_t *seq;

	uint64_t *idx;
} osg_t;

typedef struct {
    osg_t *g;
}scg_t;
typedef struct {
	kvec_t(uint64_t) avoid;
    kvec_pe_hit r_hits, u_hits;
    ma_ug_t *ug;
    asg_t *r_g;
    scg_t sg;
}horder_t;

horder_t *init_horder_t(kvec_pe_hit *i_hits, uint64_t i_hits_uid_bits, uint64_t i_hits_pos_mode, 
asg_t *i_rg, ma_ug_t* i_ug, bubble_type* bub, ug_opt_t *opt, uint32_t round);
void destory_horder_t(horder_t **h);
#endif
