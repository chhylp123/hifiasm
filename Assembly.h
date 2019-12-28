#ifndef __ASSEMBLY__
#define __ASSEMBLY__
#include "CommandLines.h"

#define FORWARD 0
#define REVERSE_COMPLEMENT (0x8000000000000000)

#define Get_Cigar_Type(RECORD) (RECORD&3)
#define Get_Cigar_Length(RECORD) (RECORD>>2)

void Counting_multiple_thr();
void Build_hash_table_multiple_thr();
///int load_pre_cauculated_index();
void Overlap_calculate_multipe_thr();
void Correct_Reads(int last_round);
#endif
