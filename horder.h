#ifndef __HORDER__
#define __HORDER__
#include <stdint.h>
#include "hic.h"

#define get_hit_srev(x, k) ((x).a.a[(k)].s>>63)
#define get_hit_slen(x, k) ((x).a.a[(k)].len>>32)
#define get_hit_suid(x, k) (((x).a.a[(k)].s<<1)>>(64 - (x).uID_bits))
#define get_hit_spos(x, k) ((x).a.a[(k)].s & (x).pos_mode)
#define get_hit_spos_e(x, k) (get_hit_srev((x),(k))?\
        ((get_hit_spos((x),(k))+1>=get_hit_slen((x),(k)))?\
            (get_hit_spos((x),(k))+1-get_hit_slen((x),(k))):0)\
                    :(get_hit_spos((x),(k))+get_hit_slen((x),(k))-1))

#define get_hit_erev(x, k) ((x).a.a[(k)].e>>63)
#define get_hit_elen(x, k) ((uint32_t)((x).a.a[(k)].len))
#define get_hit_euid(x, k) (((x).a.a[(k)].e<<1)>>(64 - (x).uID_bits))
#define get_hit_epos(x, k) ((x).a.a[(k)].e & (x).pos_mode)
#define get_hit_epos_e(x, k) (get_hit_erev((x),(k))?\
        ((get_hit_epos((x),(k))+1>=get_hit_elen((x),(k)))?\
            (get_hit_epos((x),(k))+1-get_hit_elen((x),(k))):0)\
                    :(get_hit_epos((x),(k))+get_hit_elen((x),(k))-1))

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
    // kvec_t(uint64_t) occ;
    // kvec_t(uint8_t) hf;
    kvec_pe_hit r_hits, u_hits;
    ma_ug_t *ug;
    asg_t *r_g;
    scg_t sg;
}horder_t;

horder_t *init_horder_t(kvec_pe_hit *i_hits, uint64_t i_hits_uid_bits, uint64_t i_hits_pos_mode, 
asg_t *i_rg, ma_ug_t* i_ug, bubble_type* bub, kv_u_trans_t *ref, ug_opt_t *opt, uint32_t round);
void destory_horder_t(horder_t **h);
void horder_clean_sg_by_utg(asg_t *sg, ma_ug_t *ug);
kvec_pe_hit *get_r_hits_for_trio(kvec_pe_hit *u_hits, asg_t* r_g, ma_ug_t* ug, bubble_type* bub, uint64_t uID_bits, uint64_t pos_mode);
void update_switch_unitig(ma_ug_t *ug, asg_t *rg, kvec_pe_hit *hits, kv_u_trans_t *k_trans, uint64_t cutoff_s, uint64_t cutoff_e,
uint64_t min_ulen, double boundaryRate);
kvec_pe_hit *get_r_hits_order(kvec_pe_hit *uhits, uint64_t hits_uid_bits, uint64_t hits_pos_mode, 
asg_t *rg, ma_ug_t* ug, bubble_type* bub);
void ha_aware_order(kvec_pe_hit *r_hits, asg_t *rg, ma_ug_t *ug_fa, ma_ug_t *ug_mo, kv_u_trans_t *ref, 
ug_opt_t *opt, uint32_t round);
spg_t *horder_utg(kvec_pe_hit *i_hits, uint64_t i_hits_uid_bits, uint64_t i_hits_pos_mode, 
asg_t *i_rg, ma_ug_t* i_ug, bubble_type* bub, ug_opt_t *opt);
#endif
