#ifndef __ASSEMBLY__
#define __ASSEMBLY__

#define FORWARD 0
#define REVERSE_COMPLEMENT (0x8000000000000000)



void Counting_multiple_thr();
void Build_hash_table_multiple_thr();
int load_pre_cauculated_index();
void Overlap_calculate_multipe_thr();
/********************************for debug***************************************/
void Verify_Counting();

/********************************for debug*****************************************/
void verify_Position_hash_table();
#endif
