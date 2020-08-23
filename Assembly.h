#ifndef __ASSEMBLY__
#define __ASSEMBLY__
#include "CommandLines.h"

#define FORWARD 0
#define REVERSE_COMPLEMENT (0x8000000000000000)

#define Get_Cigar_Type(RECORD) (RECORD&3)
#define Get_Cigar_Length(RECORD) (RECORD>>2)

#define RESEED_DP 4
#define RESEED_PEAK_RATE 0.15
#define RESEED_LEN 2000
#define RESEED_HP_RATE 0.9

int ha_assemble(void);

#endif
