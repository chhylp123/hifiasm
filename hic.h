#ifndef __HIC__
#define __HIC__
#include <stdint.h>
#include "Overlaps.h"

void push_hc_edge(hc_linkeage* x, uint64_t uID, int weight, int dir, uint64_t* d);
void hic_analysis(ma_ug_t *ug);

#endif
