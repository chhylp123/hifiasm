#ifndef __CORRECT__
#define __CORRECT__
#include <stdint.h>
#include "Hash_Table.h"
#include "Levenshtein_distance.h"




#define MAX(x, y) ((x >= y)?x:y)
#define MIN(x, y) ((x <= y)?x:y)
#define OVERLAP(x_start, x_end, y_start, y_end) (MIN(x_end, y_end) - MAX(x_start, y_start) + 1) 
///#define OVERLAP(x_start, x_end, y_start, y_end) MIN(x_end, y_end) - MAX(x_start, y_start) + 1 

typedef struct
{
    uint64_t* overlapID;
    uint64_t length;
    uint64_t lengthNT;
    uint64_t size;
    uint64_t start_i;
    char overlap_region[WINDOW + THRESHOLD*2 + 10];
    char overlap_region_group[GROUP_SIZE][WINDOW + THRESHOLD*2 + 10];
    char path[WINDOW + THRESHOLD*2 + 10];
    int path_length;
    Word matrix_bit[((WINDOW + 10)<<3)];
    __m128i Peq_SSE[256];
} Correct_dumy;


void correct_overlap(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read, 
                        long long* matched_overlap_0, long long* matched_overlap_1, 
                        long long* potiental_matched_overlap_0, long long* potiental_matched_overlap_1);
void init_Correct_dumy(Correct_dumy* list);
void destory_Correct_dumy(Correct_dumy* list);
void clear_Correct_dumy(Correct_dumy* list, overlap_region_alloc* overlap_list);
void pre_filter_by_nearby(k_mer_pos* new_n_list, k_mer_pos* old_n_list, uint64_t n_length, uint64_t n_end_pos, UC_Read* g_read, 
All_reads* R_INF, Correct_dumy* dumy, uint64_t* new_n_length);
void pre_filter_by_nearby_single(k_mer_pos* new_n_list, k_mer_pos* old_n_list, uint64_t n_length, uint64_t n_end_pos, UC_Read* g_read, 
All_reads* R_INF, Correct_dumy* dumy, uint64_t* new_n_length);

/**********************for prefilter************************ */
void destory_k_mer_pos_list_alloc_prefilter(k_mer_pos_list_alloc* list);
void append_k_mer_pos_list_alloc_prefilter(k_mer_pos_list_alloc* list, k_mer_pos* n_list, uint64_t n_length, 
uint64_t n_end_pos, uint8_t n_direction, UC_Read* g_read, All_reads* R_INF, Correct_dumy* dumy);
/**********************for prefilter************************ */

#endif