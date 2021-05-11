#ifndef __HORDER__
#define __HORDER__
#include <stdint.h>
#include "hic.h"
typedef struct {
    kvec_pe_hit r_hits, u_hits;
    ma_ug_t* ug;
    asg_t* r_g;
    kvec_t(uint8_t) hp;
}horder_t;

horder_t *init_horder_t(kvec_pe_hit *i_hits, uint64_t i_hits_uid_bits, uint64_t i_hits_pos_mode, 
asg_t *i_rg, ma_ug_t* i_ug, bubble_type* bub, ug_opt_t *opt);
void destory_horder_t(horder_t **h);
#endif
