#ifndef __HIC__
#define __HIC__
#include <stdint.h>
#include "Overlaps.h"

#define kdq_clear(q) ((q)->count = (q)->front = 0)
#define kv_malloc(v, s) ((v).n = 0, (v).m = (s), MALLOC((v).a, (s)))
#define RC_0 0
#define RC_1 1
#define RC_2 2

hc_edge* get_hc_edge(hc_links* link, uint64_t src, uint64_t dest, uint64_t dir);
void push_hc_edge(hc_linkeage* x, uint64_t uID, int weight, int dir, uint64_t* d);
void hic_analysis(ma_ug_t *ug, asg_t* read_g, hc_links* link);
void hic_benchmark(ma_ug_t *ug, asg_t* read_g);

typedef struct {
    uint32_t* index;
    ma_ug_t* ug;
    kvec_t(uint32_t) list;
    kvec_t(uint32_t) num;
    kvec_t(uint64_t) pathLen;
    kvec_t(uint64_t) b_s_idx;
    uint64_t f_bub, b_bub;
    asg_t *b_g;
} bubble_type;
#define P_het(B) ((B).num.n)
#define M_het(B) ((B).num.n + 1)
// #define IF_BUB(ID, B) ((B).index[(ID)] < (B).num.n)
// #define IF_HET(ID, B) ((B).index[(ID)] == (B).num.n)
// #define IF_HOM(ID, B) ((B).index[(ID)] > (B).num.n)
#define IF_BUB(ID, B) ((B).index[(ID)] < (B).f_bub+1)
#define IF_HET(ID, B) ((B).index[(ID)] == (B).f_bub+1)
#define IF_HOM(ID, B) ((B).index[(ID)] > (B).f_bub+1)
#define Get_bub_num(RECORD) ((RECORD).num.n-1)
void get_bubbles(bubble_type* bub, uint64_t id, uint32_t* beg, uint32_t* sink, uint32_t** a, uint32_t* n, uint64_t* pathBase);

#endif
