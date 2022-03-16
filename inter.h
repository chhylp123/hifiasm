#ifndef __INTER__
#define __INTER__
#include "Overlaps.h"
#include "Process_Read.h"

void ul_resolve(ma_ug_t *ug, const asg_t *rg, const ug_opt_t *uopt, int hap_n);
void ul_load(const ug_opt_t *uopt);
uint64_t* get_hifi2ul_list(all_ul_t *x, uint64_t hid, uint64_t* a_n);

#endif
