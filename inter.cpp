#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <zlib.h>
#include "inter.h"
#include "Overlaps.h"
#include "CommandLines.h"
#include "htab.h"

void uidx_build(ma_ug_t *ug, int is_HPC, int k, int w, int hap_n)
{
    int flag = asm_opt.flag;
    asm_opt.flag |= HA_F_NO_HPC;
    ug->u.h = hap_n;
    ha_flt_tab = ha_ft_ug_gen(&asm_opt, &(ug->u), is_HPC, k, w, hap_n, hap_n*10);
    ha_idx = ha_pt_ug_gen(&asm_opt, ha_flt_tab, &(ug->u), is_HPC, k, w, hap_n);
    asm_opt.flag = flag;
}

void uidx_destory()
{
    ha_ft_destroy(ha_flt_tab); 
    ha_pt_destroy(ha_idx);
}


void ul_resolve(ma_ug_t *ug, int hap_n)
{
    uidx_build(ug, 1, 63, 63, hap_n);
    
    uidx_destory();
}