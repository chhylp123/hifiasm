#ifndef __ASSEMBLY__
#define __ASSEMBLY__
#include "CommandLines.h"

#define FORWARD 0
#define REVERSE_COMPLEMENT (0x8000000000000000)

#define Get_Cigar_Type(RECORD) (RECORD&3)
#define Get_Cigar_Length(RECORD) (RECORD>>2)

int ha_assemble(void);

#endif
