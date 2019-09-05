#ifndef __CORRECT__
#define __CORRECT__
#include <stdint.h>
#include "Hash_Table.h"
#include "Levenshtein_distance.h"
#include "POA.h"
#include "Process_Read.h"

#define CORRECT_THRESHOLD 0.70
#define MIN_COVERAGE_THRESHOLD 4
#define CORRECT_INDEL_LENGTH 2
#define MISMATCH 1
#define INSERTION 2
#define DELETION 3

#define FLAG_THRE 1

#define MAX(x, y) ((x >= y)?x:y)
#define MIN(x, y) ((x <= y)?x:y)
#define OVERLAP(x_start, x_end, y_start, y_end) (MIN(x_end, y_end) - MAX(x_start, y_start) + 1) 
///#define OVERLAP(x_start, x_end, y_start, y_end) MIN(x_end, y_end) - MAX(x_start, y_start) + 1 

#define Get_MisMatch_Base(RECORD) (s_H[(RECORD>>3)])
#define Get_Match_Base(RECORD) (s_H[(RECORD&7)])



typedef struct
{
    /**[0-1] bits are type:**/
    /**[2-31] bits are length**/
    char current_operation;
    int current_operation_length;
    uint32_t* record;
    uint64_t size;
    uint64_t length;
    uint32_t new_read_length;


    char* lost_base;
    uint64_t lost_base_size;
    uint64_t lost_base_length;
    

}Cigar_record;


typedef struct
{
    ////the position of snp in read itself
    uint32_t site;
    ////the overlapID
    uint32_t overlapID;
    ////the position of snp in that overlap
    uint32_t overlapSite;
    ///there are several types: 0: equal to read 1: not equal to read, but it is a mismatch 2: is a gap
    uint8_t type;
    ///misbase
    char misBase;
}haplotype_evdience;



typedef struct
{
    ///the id of this snp
    uint32_t id;
    uint32_t overlap_num;
    uint32_t occ_0;
    uint32_t occ_1;
    uint32_t occ_2;
    int score;
    ////the position of snp in read itself
    uint32_t site;
}
SnpStats;


#define Get_SNP_Martix_Size(matrix) (matrix.snp * matrix.overlap)
#define Get_SNP_Vector(matrix, i) (matrix.snp_matrix + matrix.overlap * i)
#define Get_SNP_Vector_Length(matrix) (matrix.overlap)
#define Get_Result_SNP_Vector(matrix) (matrix.snp_matrix + matrix.overlap*matrix.snp)

typedef struct
{
    haplotype_evdience* list;
    uint32_t sub_list_start;
    uint32_t sub_list_length;
    uint32_t length;
    uint32_t size;
    
    uint32_t flag[WINDOW];

    uint32_t available_snp;
    uint32_t core_snp;
    uint32_t snp;
    uint32_t overlap;
    int8_t* snp_matrix;
    uint32_t snp_matrix_size;
    SnpStats* snp_stat; 
    SnpStats result_stat;
    uint32_t snp_stat_size;
}
haplotype_evdience_alloc;

inline int filter_snp(int x, int y, int total)
{
    double available;

    if(x <= y)
    {
        available = x;
    }
    else
    {
        available = y;
    }
    double threshold = 0.30;
    available = available/((double)(total));
    if(available <= threshold && available < 6)
    {
        return 0;
    }
    return 1;
}


inline int filter_one_snp(int occ_0, int occ_1, int total)
{
    

    double available;

    if(occ_0 <= occ_1)
    {
        available = occ_0;
    }
    else
    {
        available = occ_1;
    }
    double threshold = 0.35;
    available = available/((double)(total));

    

    if(available < threshold || occ_0 < MIN_COVERAGE_THRESHOLD + 1 || total < 10)
    {
        return 0;
    }
    return 1;
}

inline void InsertSNPVector(haplotype_evdience_alloc* h, haplotype_evdience* sub_list, long long sub_length, char misBase)
{
    if(sub_length <= 0)
        return;
    long long i = 0;
    h->snp_stat[h->available_snp].id = h->available_snp;
    h->snp_stat[h->available_snp].occ_0 = 0;
    h->snp_stat[h->available_snp].occ_1 = 0;
    h->snp_stat[h->available_snp].occ_2 = 0;
    h->snp_stat[h->available_snp].overlap_num = 0;
    
    h->snp_stat[h->available_snp].site = sub_list[0].site;

    int8_t* vector = Get_SNP_Vector((*h), h->available_snp);
    for (i = 0; i < sub_length; i++)
    {
        if(sub_list[i].type == 0)
        {
            vector[sub_list[i].overlapID] = 0;
            h->snp_stat[h->available_snp].occ_0++;
        }
        else if(sub_list[i].type == 1 && sub_list[i].misBase == misBase)
        {
            vector[sub_list[i].overlapID] = 1;
            h->snp_stat[h->available_snp].occ_1++;
        }
        else
        {
            vector[sub_list[i].overlapID] = 2;
            h->snp_stat[h->available_snp].occ_2++;
        }
        h->snp_stat[h->available_snp].overlap_num++;
    }


    int new_occ_0 = h->snp_stat[h->available_snp].occ_0 + 1;
    int new_occ_1 = h->snp_stat[h->available_snp].occ_1;

    if(filter_snp(new_occ_0, new_occ_1, new_occ_0 + new_occ_1) == 0)
    {
        h->snp_stat[h->available_snp].score = -1;
    }
    else
    {
        h->core_snp++;


        double consensus = new_occ_0 + new_occ_1 - abs(new_occ_0 - new_occ_1);

        consensus = consensus /((double)(new_occ_0 + new_occ_1));

        ///50% vs 50%
        if(new_occ_0 == new_occ_1)
        {
            consensus = consensus + 0.25;
        }
        else if(consensus >= 0.8)
        {
            consensus = consensus + 0.2;
        }
        else if(consensus >= 0.6)
        {
            consensus = consensus + 0.15;
        }
        else if(consensus >= 0.4)
        {
            consensus = consensus + 0.1;
        }
        else if(consensus >= 0.2)
        {
            consensus = consensus + 0.05;
        }
        

        consensus= consensus*((double)(new_occ_0 + new_occ_1));

        h->snp_stat[h->available_snp].score = consensus;
    }


    
    h->available_snp++;
}

inline void SetSnpMatrix(haplotype_evdience_alloc* h, long long snp_num, long long overlap_num)
{
    long long new_size = (snp_num + 1)* overlap_num;

    if(h->snp_matrix_size < new_size)
    {
        h->snp_matrix_size = new_size;
        h->snp_matrix = (int8_t*)realloc(h->snp_matrix, h->snp_matrix_size);
    }

    if(h->snp_stat_size < snp_num)
    {
        h->snp_stat_size = snp_num;
        h->snp_stat = (SnpStats*)realloc(h->snp_stat, h->snp_stat_size * sizeof(SnpStats));
    }

    ///h->snp may be different with the number of snp vector
    ///since some snps have been filtered
    h->snp = snp_num;
    h->overlap = overlap_num;
    h->available_snp = 0;
    h->core_snp = 0;

    memset(h->snp_matrix, -1, h->snp * h->overlap);

    
}

inline void InitHaplotypeEvdience(haplotype_evdience_alloc* h)
{
    h->snp = 0;
    h->available_snp = 0;
    h->overlap = 0;
    h->snp_matrix_size = 0;
    h->snp_stat_size = 0;
    h->snp_matrix = NULL;
    h->snp_stat = NULL;


    h->sub_list_start = 0;
    h->sub_list_length = 0;
    h->length = 0;
    h->size = 100;
    h->list = (haplotype_evdience*)calloc(h->size, sizeof(haplotype_evdience));
    memset(h->flag, 0, WINDOW * sizeof(uint32_t));
}

inline void StarSubListHaplotypeEvdience(haplotype_evdience_alloc* h)
{
    h->sub_list_start = h->length;
}

inline void EndSubListHaplotypeEvdience(haplotype_evdience_alloc* h)
{
    h->sub_list_length = h->length - h->sub_list_start;
}


inline void destoryHaplotypeEvdience(haplotype_evdience_alloc* h)
{
    free(h->list);
    free(h->snp_stat);
    free(h->snp_matrix);
}

inline void ResizeInitHaplotypeEvdience(haplotype_evdience_alloc* h)
{
    h->snp = 0;
    h->length = 0;
    h->sub_list_start = 0;
    h->sub_list_length = 0;
    memset(h->flag, 0, WINDOW * sizeof(uint32_t));
}

inline void RsetInitHaplotypeEvdienceFlag(haplotype_evdience_alloc* h)
{
    memset(h->flag, 0, WINDOW * sizeof(uint32_t));
}

inline void addHaplotypeEvdience(haplotype_evdience_alloc* h, haplotype_evdience* ev)
{
    uint32_t new_length = h->length + 1;
    if(new_length > h->size)
    {
        h->size = h->size * 2;
        if(h->size < new_length)
        {
            h->size = new_length;
        }
        h->list = (haplotype_evdience*)realloc(h->list, sizeof(haplotype_evdience)*h->size);
    }

    h->list[h->length] = (*ev);
    h->length++;
}

typedef struct
{
    char* corrected_read;
    long long corrected_read_length;
    long long last_boundary_length;
    long long corrected_read_size;
    long long corrected_base;

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

typedef struct
{
    Correct_dumy dumy;
    Cigar_record cigar;
    Cigar_record tmp_cigar;
    long long obtained_cigar_length;
}
Round2_alignment;

void init_Round2_alignment(Round2_alignment* h);
void destory_Round2_alignment(Round2_alignment* h);
void clear_Round2_alignment(Round2_alignment* h);



void correct_overlap(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read, Graph* g,
                        long long* matched_overlap_0, long long* matched_overlap_1, 
                        long long* potiental_matched_overlap_0, long long* potiental_matched_overlap_1,
                        Cigar_record* current_cigar, haplotype_evdience_alloc* hap,
                        Round2_alignment* second_round);
void init_Correct_dumy(Correct_dumy* list);
void destory_Correct_dumy(Correct_dumy* list);
void clear_Correct_dumy(Correct_dumy* list, overlap_region_alloc* overlap_list);
void clear_Correct_dumy_pure(Correct_dumy* list);
void pre_filter_by_nearby(k_mer_pos* new_n_list, k_mer_pos* old_n_list, uint64_t n_length, uint64_t n_end_pos, UC_Read* g_read, 
All_reads* R_INF, Correct_dumy* dumy, uint64_t* new_n_length);
void pre_filter_by_nearby_single(k_mer_pos* new_n_list, k_mer_pos* old_n_list, uint64_t n_length, uint64_t n_end_pos, UC_Read* g_read, 
All_reads* R_INF, Correct_dumy* dumy, uint64_t* new_n_length);
void get_seq_from_Graph(Graph* backbone, Correct_dumy* dumy);

void init_Cigar_record(Cigar_record* dummy);
void destory_Cigar_record(Cigar_record* dummy);
void clear_Cigar_record(Cigar_record* dummy);




inline void add_new_cell_to_cigar_record(Cigar_record* dummy, uint32_t len, uint32_t type)
{
    uint32_t tmp;
    tmp = len;
    tmp = tmp << 2;
    tmp = tmp | type;

    dummy->length++;


    if(dummy->length > dummy->size)
    {
        dummy->size = dummy->size * 2;
        dummy->record = (uint32_t*)realloc(dummy->record, dummy->size*sizeof(uint32_t));
    }

    dummy->record[dummy->length - 1] = tmp;
}

inline void add_existing_cell_to_cigar_record(Cigar_record* dummy, uint32_t len, uint32_t type)
{
    uint32_t tmp;

    tmp = dummy->record[dummy->length - 1] >> 2;
    tmp = tmp + len;
    tmp = tmp << 2;
    tmp = tmp | type;
    dummy->record[dummy->length - 1] = tmp;
}


inline void add_new_cell_to_cigar_record_with_different_base(Cigar_record* dummy, uint32_t len, uint32_t type, char* seq)
{
    uint32_t tmp;
    tmp = len;
    tmp = tmp << 2;
    tmp = tmp | type;
   

    dummy->length++;

    if(dummy->length > dummy->size)
    {
        dummy->size = dummy->size * 2;
        dummy->record = (uint32_t*)realloc(dummy->record, dummy->size*sizeof(uint32_t));
    }

    dummy->record[dummy->length - 1] = tmp;

    

    if (dummy->lost_base_length + len> dummy->lost_base_size)
    {
        dummy->lost_base_size = (dummy->lost_base_length + len) * 2;
        dummy->lost_base = (char*)realloc(dummy->lost_base, dummy->lost_base_size*sizeof(char));
    }

    int i = 0;
    for (i = 0; i < len; i++, dummy->lost_base_length++)
    {
        dummy->lost_base[dummy->lost_base_length] = seq[i];
    }
}

inline void add_existing_cell_to_cigar_record_with_different_base(Cigar_record* dummy, uint32_t len, uint32_t type, char* seq)
{
    uint32_t tmp;

    tmp = dummy->record[dummy->length - 1] >> 2;
    tmp = tmp + len;
    tmp = tmp << 2;
    tmp = tmp | type;
    dummy->record[dummy->length - 1] = tmp;

    if (dummy->lost_base_length + len> dummy->lost_base_size)
    {
        dummy->lost_base_size = (dummy->lost_base_length + len) * 2;
        dummy->lost_base = (char*)realloc(dummy->lost_base, dummy->lost_base_size*sizeof(char));
    }

    int i = 0;
    for (i = 0; i < len; i++, dummy->lost_base_length++)
    {
        dummy->lost_base[dummy->lost_base_length] = seq[i];
    }
}


/***
 type:
 0. match
 1. mismatch
 2. insertion
 3. deletion
 ***/
inline void add_cigar_record(char* seq, uint32_t len, Cigar_record* dummy, uint32_t type)
{

    
    uint32_t tmp;
    
    if(type == 0)///match
    {
        ///add to existing cell, just increase length
        if(dummy->current_operation == type)
        {
            add_existing_cell_to_cigar_record(dummy, len, type);
        }
        else  ///add to new cell
        {
            add_new_cell_to_cigar_record(dummy, len, type);
        }

        dummy->new_read_length += len;
    }
    else if(type == 1)///mismatch
    {
        ///add to existing cell, just increase length
        ///and add different bases
        if(dummy->current_operation == type)
        {
            add_existing_cell_to_cigar_record_with_different_base(dummy, len, type, seq);
        }
        else
        {
            add_new_cell_to_cigar_record_with_different_base(dummy, len, type, seq);
        }

        dummy->new_read_length += len;
    }
    else if(type == 3)///deletion, the bases in previous read will be removed
    {
        ///add to existing cell, just increase length
        ///and add different bases
        if(dummy->current_operation == type)
        {
            add_existing_cell_to_cigar_record_with_different_base(dummy, len, type, seq);
        }
        else
        {
            add_new_cell_to_cigar_record_with_different_base(dummy, len, type, seq);
        }
    }
    else if(type == 2)///insertion
    {
        /**
        ///add to existing cell, just increase length
        if(dummy->current_operation == type)
        {
            add_existing_cell_to_cigar_record(dummy, len, type);
        }
        else  ///add to new cell
        {
            add_new_cell_to_cigar_record(dummy, len, type);
        }
        **/
       ///add to existing cell, just increase length
        ///and add different bases
        if(dummy->current_operation == type)
        {
            add_existing_cell_to_cigar_record_with_different_base(dummy, len, type, seq);
        }
        else
        {
            add_new_cell_to_cigar_record_with_different_base(dummy, len, type, seq);
        }

        dummy->new_read_length += len;
    }
    
    dummy->current_operation = type;





    
}


/**********************for prefilter************************ */
void destory_k_mer_pos_list_alloc_prefilter(k_mer_pos_list_alloc* list);
void append_k_mer_pos_list_alloc_prefilter(k_mer_pos_list_alloc* list, k_mer_pos* n_list, uint64_t n_length, 
uint64_t n_end_pos, uint8_t n_direction, UC_Read* g_read, All_reads* R_INF, Correct_dumy* dumy);
/**********************for prefilter************************ */

#endif