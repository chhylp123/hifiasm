#ifndef __HIC__
#define __HIC__
#include <stdint.h>
#include "Overlaps.h"

#define kdq_clear(q) ((q)->count = (q)->front = 0)
#define kv_malloc(v, s) ((v).n = 0, (v).m = (s), MALLOC((v).a, (s)))

hc_edge* get_hc_edge(hc_links* link, uint64_t src, uint64_t dest, uint64_t dir);
void push_hc_edge(hc_linkeage* x, uint64_t uID, int weight, int dir, uint64_t* d);
void hic_analysis(ma_ug_t *ug, asg_t* read_g, hc_links* link);
void hic_benchmark(ma_ug_t *ug, asg_t* read_g);

#endif
