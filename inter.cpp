#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <zlib.h>
#include "inter.h"
#include "Overlaps.h"
#include "CommandLines.h"
#include "htab.h"

void uidx_build(ma_ug_t *ug, int hap_n)
{
    int flag = asm_opt.flag;
    asm_opt.flag |= HA_F_NO_HPC;
    ha_flt_tab = ha_ft_ug_gen(&asm_opt, &(ug->u), hap_n);
    ha_idx = ha_pt_ug_gen(&asm_opt, ha_flt_tab, &(ug->u), hap_n);






    asm_opt.flag = flag;
    ha_ft_destroy(ha_flt_tab); 
    ha_pt_destroy(ha_idx);
}
