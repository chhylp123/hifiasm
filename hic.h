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
void push_hc_edge(hc_linkeage* x, uint64_t uID, double weight, int dir, uint64_t* d);
void hic_analysis(ma_ug_t *ug, asg_t* read_g, hc_links* link);
void hic_benchmark(ma_ug_t *ug, asg_t* read_g);

typedef struct {
    double w;
    uint32_t id, occ;
    ///uint32_t *bid, bid_n;
    ma_utg_t *u;
    uint64_t l_d, r_d;
}chain_hic_w_type;

typedef struct {
    size_t n, m;
    chain_hic_w_type* a;
    uint32_t max_bub_id;
    uint32_t *chain_idx, u_n;
}chain_hic_warp;

typedef struct {
    long long g_occ, b_occ;
    uint64_t id;
    uint8_t del;
}chain_w_type;

typedef struct {
    uint32_t* index, round_id, n_round;
    ma_ug_t* ug;
    kvec_t(uint32_t) list;
    kvec_t(uint32_t) num;
    kvec_t(uint64_t) pathLen;
    kvec_t(uint64_t) b_s_idx;
    uint64_t s_bub, f_bub, b_bub, b_end_bub, tangle_bub, cross_bub, mess_bub;
    uint32_t check_het;
    asg_t *b_g;
    ma_ug_t* b_ug;
    kvec_t(chain_w_type) chain_weight;
    chain_hic_warp c_w;
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
int load_hc_links(hc_links* link, const char *fn);
void write_hc_links(hc_links* link, const char *fn);
void destory_bubbles(bubble_type* bub);
void identify_bubbles(ma_ug_t* ug, bubble_type* bub, hc_links* link);
void resolve_bubble_chain_tangle(ma_ug_t* ug, bubble_type* bub);
uint32_t connect_bub_occ(bubble_type* bub, uint32_t root_id, uint32_t check_het);
void get_bub_id(bubble_type* bub, uint32_t root, uint64_t* id0, uint64_t* id1, uint32_t check_het);
void update_bubble_chain(ma_ug_t* ug, bubble_type* bub, uint32_t is_middle, uint32_t is_end);
void set_b_utg_weight_flag(bubble_type* bub, buf_t* b, uint32_t v, uint8_t* vis_flag, uint32_t flag, uint32_t* occ);

#endif
