#ifndef __CORRECT__
#define __CORRECT__
#include <stdint.h>
#include "Hash_Table.h"
#include "Levenshtein_distance.h"

#define WINDOW 350
#define THRESHOLD  15
#define GROUP_SIZE 4

typedef struct
{
    uint64_t* overlapID;
    uint64_t length;
    uint64_t lengthNT;
    uint64_t size;
    uint64_t start_i;
    char overlap_region[WINDOW + THRESHOLD*2 + 10];
    char overlap_region_group[GROUP_SIZE][WINDOW + THRESHOLD*2 + 10];
    __m128i Peq_SSE[256];
} Correct_dumy;


void correct_overlap(overlap_region_alloc* overlap_list, All_reads* R_INF, UC_Read* g_read, Correct_dumy* dumy);

void init_Correct_dumy(Correct_dumy* list);
void destory_Correct_dumy(Correct_dumy* list);
void clear_Correct_dumy(Correct_dumy* list, overlap_region_alloc* overlap_list);


#endif