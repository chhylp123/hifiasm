#ifndef __CORRECT__
#define __CORRECT__
#include <stdint.h>
#include "Hash_Table.h"
#include "Levenshtein_distance.h"
#include "POA.h"
#include "Process_Read.h"
#include "Correct.h"

//#define CORRECT_THRESHOLD 0.70
#define CORRECT_THRESHOLD 0.60
///#define CORRECT_THRESHOLD_SECOND 0.55
#define CORRECT_THRESHOLD_HOMOPOLYMER 0.515
#define MIN_COVERAGE_THRESHOLD 3
#define CORRECT_INDEL_LENGTH 2
#define MISMATCH 1
#define INSERTION 2
#define DELETION 3

#define WINDOW_MAX_SIZE (WINDOW + (int)(1.0 / HA_MIN_OV_DIFF) + 3) // TODO: why 1/max_ov_diff?

///#define FLAG_THRE 0

#define MAX(x, y) (((x) >= (y))?(x):(y))
#define MIN(x, y) (((x) <= (y))?(x):(y))
#define DIFF(x, y) ((MAX((x), (y))) - (MIN((x), (y))))
#define OVERLAP(x_start, x_end, y_start, y_end) (MIN(x_end, y_end) - MAX(x_start, y_start) + 1) 
///#define OVERLAP(x_start, x_end, y_start, y_end) MIN(x_end, y_end) - MAX(x_start, y_start) + 1 

#define Get_MisMatch_Base(RECORD) (s_H[(RECORD>>3)])
#define Get_Match_Base(RECORD) (s_H[(RECORD&7)])
#define Coverage_Threshold(coverage, r_len) (coverage*r_len*1.1)

 

#define Get_Max_DP_Value(RECORD) (RECORD>>32)
#define Get_Max_DP_ID(RECORD) (RECORD&(uint64_t)0xffffffff)

#define Adjust_Threshold(threshold, x_len) ((threshold == 0 && x_len >= 4)? 1: threshold)

typedef struct
{
    long long read_length;
    long long window_length;
    long long window_num;
    long long window_start;
    long long window_end;
    long long tail_length;
    int terminal;
}Window_Pool;

inline void init_Window_Pool(Window_Pool* dumy, long long read_length, long long window_length, long long tail_length)
{
    dumy->terminal = 0;
    dumy->read_length = read_length;
    dumy->window_length = window_length;
    dumy->tail_length = tail_length;

    dumy->window_start = 0;
    dumy->window_end = dumy->window_length - 1;
    if (dumy->window_end >= dumy->read_length)
    {
        dumy->window_end = dumy->read_length - 1;
    }

    dumy->window_num = (dumy->read_length + dumy->window_length - 1) / dumy->window_length;
}

inline int get_Window(Window_Pool* dumy, long long* w_beg, long long* w_end)
{   
    (*w_beg) = dumy->window_start;
    (*w_end) = dumy->window_end;

    if(dumy->window_end == dumy->read_length - 1)
    {
        if(dumy->terminal == 1 || dumy->read_length == 0)
        {
            return 0;
        }
        else if(dumy->terminal == 0)
        {
            dumy->terminal = 1;
        }
    }


    dumy->window_start = dumy->window_start + dumy->window_length;
    dumy->window_end = dumy->window_end + dumy->window_length;
    if (dumy->window_end >= dumy->read_length)
    {
        dumy->window_end = dumy->read_length - 1;
    }
    return 1;
}


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
  long long length;
  long long size;
  Cigar_record* buffer;
}Cigar_record_alloc;


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
    uint32_t homopolymer_num;
    uint32_t non_homopolymer_num;
    int score;
    ////the position of snp in read itself
    uint32_t site;
    uint8_t is_homopolymer;
}
SnpStats;



typedef struct
{
    uint32_t beg;
    uint32_t end;
    uint32_t occ_0;
    uint32_t occ_1;
    uint32_t homopolymer_num;
    uint32_t non_homopolymer_num;
    uint32_t is_remove;
}
Snp_ID_Vector;


typedef struct
{
    long long IDs_size;
    long long IDs_length;
    long long max_snp_id;
    Snp_ID_Vector* IDs;

    long long buffer_size;
    long long buffer_length;
    uint32_t* buffer;
}
Snp_ID_Vector_Alloc;

#define Get_DP_Backtrack_Column(matrix, i) (matrix.backtrack + matrix.snp_num * i)
#define Get_DP_Backtrack_Column_Length(matrix, i) (matrix.snp_num)


typedef struct
{
    // uint32_t snp_size;
    // uint32_t snp_num;
    // uint32_t* max;
    // uint32_t* colum_len;
    // uint32_t* colum;
    // uint32_t matrix_size;

    uint32_t snp_num;
    
    uint8_t* visit;
    uint32_t* max;
    uint64_t* max_for_sort;



    uint32_t snp_size;
    
    uint32_t* backtrack_length;
    uint32_t* backtrack;
    uint32_t backtrack_size;


    uint32_t* buffer;
    uint32_t* max_buffer;


    int max_snp_num;
    ///int max_snp_ID;
    int max_score;
    int current_snp_num;


    Snp_ID_Vector_Alloc SNP_IDs;
}
DP_matrix;


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
    
    /****************************may have bugs********************************/
    uint8_t flag[WINDOW_MAX_SIZE];
    /****************************may have bugs********************************/

    uint32_t available_snp;
    uint32_t core_snp;
    uint32_t snp;
    uint32_t overlap;
    int8_t* snp_matrix;
    uint32_t snp_matrix_size;
    SnpStats* snp_stat; 
    SnpStats result_stat;
    uint32_t snp_stat_size;

    DP_matrix dp;
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
    ///if(available < threshold || total < 10)
    {
        return 0;
    }
    return 1;
}

inline void count_nearby_snps(haplotype_evdience_alloc* hap, uint32_t* SNPs, int SNPsLen, int* nearsnp, int* non_nearsnps)
{
    long long i, current_id, large_id, small_id;
    long long distance = 10;

    

    if(SNPsLen == 1)
    {
        (*non_nearsnps) = 1;
        (*nearsnp) = 0;
        return;
    }

    if(SNPsLen == 0)
    {
        (*nearsnp) = 0;
        (*non_nearsnps) = 0;
        return;
    }

    (*nearsnp) = 0;
    (*non_nearsnps) = 0;

    for (i = 0; i < SNPsLen; i++)
    {
        if(i > 0 && i < SNPsLen - 1)
        {
            current_id= SNPs[i];

            ///since SNPs[i - 1].site is larger than SNPs[i]
            large_id = SNPs[i - 1];
            small_id = SNPs[i + 1];
            

            if(hap->snp_stat[large_id].site - hap->snp_stat[current_id].site < distance
              ||
              hap->snp_stat[current_id].site - hap->snp_stat[small_id].site < distance)
            {
                (*nearsnp)++;
            }
            else
            {
                (*non_nearsnps)++;
            }

            if(hap->snp_stat[current_id].site > hap->snp_stat[large_id].site || 
               hap->snp_stat[current_id].site < hap->snp_stat[small_id].site)
            {
                fprintf(stderr, "error\n");
            }
        }
        else if(i == 0)
        {
            current_id= SNPs[i];
            small_id = SNPs[i + 1];
            if(hap->snp_stat[current_id].site - hap->snp_stat[small_id].site < distance)
            {
                (*nearsnp)++;
            }
            else
            {
                (*non_nearsnps)++;
            }

            if(hap->snp_stat[current_id].site < hap->snp_stat[small_id].site)
            {
                fprintf(stderr, "error\n");
            }
        }
        else
        {
            large_id = SNPs[i - 1];
            current_id= SNPs[i];
            if(hap->snp_stat[large_id].site - hap->snp_stat[current_id].site < distance)
            {
                (*nearsnp)++;
            }
            else
            {
                (*non_nearsnps)++;
            }

            if(hap->snp_stat[current_id].site > hap->snp_stat[large_id].site)
            {
                fprintf(stderr, "error\n");
            }
            
        }
        
    }

    if((*nearsnp) + (*non_nearsnps) != SNPsLen)
    {
        fprintf(stderr, "(*nearsnp): %d, (*non_nearsnps): %d, SNPsLen: %d\n", (*nearsnp), (*non_nearsnps), SNPsLen);
    }
}

inline int filter_one_snp_advance_nearby(haplotype_evdience_alloc* hap, int occ_0, int occ_1, int total,
long long homopolymer_num, long long non_homopolymer_num,
uint32_t* SNPs, int SNPsLen)
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
    int min = available;

    double threshold1 = 0.35;
    double threshold2 = 0.24;
    available = available/((double)(total));
    


    ///if((non_homopolymer_num > 0 && min >= 5) || (min >= 6))
    if(min >= 5)
    {
        if(available < threshold2 || total < 10)
        {
            return 0;
        }
    }
    else if(available < threshold1 || occ_0 < MIN_COVERAGE_THRESHOLD + 1 || total < 10)
    {
        return 0;
    }
    return 1;
}


inline int if_is_homopolymer_strict(long long site, char* read, long long read_length)
{
    long long beg, end, i;
    long long threshold = 3;

    beg = site - threshold;
    if(beg < 0)
    {
        beg = 0;
    }

    end = site + threshold;

    if(end >= read_length)
    {
        end = read_length - 1;
    }

    char f_homopolymer_ch = 0;
    long long f_homopolymer_len = 0;

    for (i = site + 1; i <= end; i++)
    {
        if(f_homopolymer_ch == 0)
        {
            f_homopolymer_ch = read[i];
            f_homopolymer_len = 1;
        }
        else 
        {
            if(read[i] != f_homopolymer_ch)
            {
                break;
            }
            else
            {
                f_homopolymer_len++;
            }
        }
    }

    char b_homopolymer_ch = 0;
    long long b_homopolymer_len = 0;

    for (i = site - 1; i >= beg; i--)
    {
        if(b_homopolymer_ch == 0)
        {
            b_homopolymer_ch = read[i];
            b_homopolymer_len = 1;
        }
        else
        {
            if(read[i] != b_homopolymer_ch)
            {
                break;
            }
            else
            {
                b_homopolymer_len++;
            }
        }
    }

    if(f_homopolymer_ch == read[site])
    {
        f_homopolymer_len++;
    }
    else if(b_homopolymer_ch == read[site])
    {
        b_homopolymer_len++;
    }

    if(f_homopolymer_len >= threshold || b_homopolymer_len >= threshold)
    {
        return 1;
    }

    if (read[site] == f_homopolymer_ch 
        && 
        b_homopolymer_ch == f_homopolymer_ch 
        && 
        (f_homopolymer_len + b_homopolymer_len >= threshold))
    {
        return 1;
    }
    
    
    return 0;
}


inline int if_is_homopolymer_repeat(long long site, char* read, long long read_length)
{
    long long beg, end, i;
    long long threshold = 3;

    beg = site - threshold;
    if(beg < 0)
    {
        beg = 0;
    }

    end = site + threshold;

    if(end >= read_length)
    {
        end = read_length - 1;
    }

    char f_homopolymer_ch = 0;
    long long f_homopolymer_len = 0;

    for (i = site + 1; i <= end; i++)
    {
        if(f_homopolymer_ch == 0)
        {
            f_homopolymer_ch = read[i];
            f_homopolymer_len = 1;
        }
        else 
        {
            if(read[i] != f_homopolymer_ch)
            {
                break;
            }
            else
            {
                f_homopolymer_len++;
            }
        }
    }

    char b_homopolymer_ch = 0;
    long long b_homopolymer_len = 0;

    for (i = site - 1; i >= beg; i--)
    {
        if(b_homopolymer_ch == 0)
        {
            b_homopolymer_ch = read[i];
            b_homopolymer_len = 1;
        }
        else
        {
            if(read[i] != b_homopolymer_ch)
            {
                break;
            }
            else
            {
                b_homopolymer_len++;
            }
        }
    }

    if(f_homopolymer_ch == read[site])
    {
        f_homopolymer_len++;
    }
    else if(b_homopolymer_ch == read[site])
    {
        b_homopolymer_len++;
    }

    if(f_homopolymer_len >= threshold || b_homopolymer_len >= threshold)
    {
        return 1;
    }

    if (read[site] == f_homopolymer_ch 
        && 
        b_homopolymer_ch == f_homopolymer_ch 
        && 
        (f_homopolymer_len + b_homopolymer_len >= threshold))
    {
        return 1;
    }
    
    
    return 0;
}

inline void InsertSNPVector(haplotype_evdience_alloc* h, haplotype_evdience* sub_list, long long sub_length, char misBase, 
UC_Read* g_read)
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

    h->snp_stat[h->available_snp].is_homopolymer = 
    if_is_homopolymer_strict(h->snp_stat[h->available_snp].site, g_read->seq, g_read->length);

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


inline int calculate_score(int new_occ_0, int new_occ_1)
{
    if(new_occ_0 + new_occ_1 == 0)
    {
        return -1;
    }

    if(filter_snp(new_occ_0, new_occ_1, new_occ_0 + new_occ_1) == 0)
    {
        return -1;
    }

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

    return consensus;
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


inline void init_SNP_IDs(Snp_ID_Vector_Alloc* SNP_IDs)
{
    SNP_IDs->max_snp_id = -1;
    SNP_IDs->buffer_length = 0;
    SNP_IDs->buffer_size = 1000;
    SNP_IDs->buffer = (uint32_t*)malloc(sizeof(uint32_t) * SNP_IDs->buffer_size);

    SNP_IDs->IDs_length = 0;
    SNP_IDs->IDs_size = 10;
    SNP_IDs->IDs = (Snp_ID_Vector*)calloc(SNP_IDs->IDs_size, sizeof(Snp_ID_Vector));

}

inline void clear_SNP_IDs(Snp_ID_Vector_Alloc* SNP_IDs)
{
    SNP_IDs->IDs_length = 0;
    SNP_IDs->buffer_length = 0;
    SNP_IDs->max_snp_id = -1;
}

inline void destory_SNP_IDs(Snp_ID_Vector_Alloc* SNP_IDs)
{
    free(SNP_IDs->buffer);
    free(SNP_IDs->IDs);
}

inline void insert_SNP_IDs_addition(Snp_ID_Vector_Alloc* SNP_IDs, uint32_t* IDs_vec, int IDs_vec_length,
 uint32_t occ_0, uint32_t occ_1, uint32_t homopolymer_num, uint32_t non_homopolymer_num)
{
    if(SNP_IDs->IDs_length + 1 > SNP_IDs->IDs_size)
    {
        SNP_IDs->IDs_size = SNP_IDs->IDs_size * 2;
        SNP_IDs->IDs = (Snp_ID_Vector*)realloc(SNP_IDs->IDs, SNP_IDs->IDs_size * sizeof(Snp_ID_Vector));
    }
    SNP_IDs->IDs[SNP_IDs->IDs_length].beg = SNP_IDs->buffer_length;
    SNP_IDs->IDs[SNP_IDs->IDs_length].end = SNP_IDs->buffer_length + IDs_vec_length - 1;
    
    SNP_IDs->IDs[SNP_IDs->IDs_length].occ_0 = occ_0;
    SNP_IDs->IDs[SNP_IDs->IDs_length].occ_1 = occ_1;
    SNP_IDs->IDs[SNP_IDs->IDs_length].homopolymer_num = homopolymer_num;
    SNP_IDs->IDs[SNP_IDs->IDs_length].non_homopolymer_num = non_homopolymer_num;
    


    if(SNP_IDs->buffer_length + IDs_vec_length > SNP_IDs->buffer_size)
    {
        SNP_IDs->buffer_size = (SNP_IDs->buffer_length + IDs_vec_length) * 2;
        SNP_IDs->buffer = (uint32_t*)realloc(SNP_IDs->buffer, SNP_IDs->buffer_size * sizeof(uint32_t));
    }
    memcpy(SNP_IDs->buffer + SNP_IDs->buffer_length, IDs_vec, IDs_vec_length*sizeof(uint32_t));


    SNP_IDs->IDs_length++;
    SNP_IDs->buffer_length += IDs_vec_length;
}


inline void insert_SNP_IDs_addition(Snp_ID_Vector_Alloc* SNP_IDs, uint32_t* IDs_vec, int IDs_vec_length)
{
    if(SNP_IDs->IDs_length + 1 > SNP_IDs->IDs_size)
    {
        SNP_IDs->IDs_size = SNP_IDs->IDs_size * 2;
        SNP_IDs->IDs = (Snp_ID_Vector*)realloc(SNP_IDs->IDs, SNP_IDs->IDs_size * sizeof(Snp_ID_Vector));
    }
    SNP_IDs->IDs[SNP_IDs->IDs_length].beg = SNP_IDs->buffer_length;
    SNP_IDs->IDs[SNP_IDs->IDs_length].end = SNP_IDs->buffer_length + IDs_vec_length - 1;


    if(SNP_IDs->buffer_length + IDs_vec_length > SNP_IDs->buffer_size)
    {
        SNP_IDs->buffer_size = (SNP_IDs->buffer_length + IDs_vec_length) * 2;
        SNP_IDs->buffer = (uint32_t*)realloc(SNP_IDs->buffer, SNP_IDs->buffer_size * sizeof(uint32_t));
    }
    memcpy(SNP_IDs->buffer + SNP_IDs->buffer_length, IDs_vec, IDs_vec_length*sizeof(uint32_t));


    SNP_IDs->IDs_length++;
    SNP_IDs->buffer_length += IDs_vec_length;
}


inline void init_DP_matrix(DP_matrix* dp, uint32_t snp_num)
{

    if(snp_num > dp->snp_size)
    {
        dp->snp_size = snp_num;
        dp->max = (uint32_t*)realloc(dp->max, dp->snp_size * sizeof(uint32_t));
        dp->max_for_sort = (uint64_t*)realloc(dp->max_for_sort, dp->snp_size * sizeof(uint64_t));
        dp->visit = (uint8_t*)realloc(dp->visit, dp->snp_size * sizeof(uint8_t));
        dp->buffer = (uint32_t*)realloc(dp->buffer, dp->snp_size * sizeof(uint32_t));
        dp->max_buffer = (uint32_t*)realloc(dp->max_buffer, dp->snp_size * sizeof(uint32_t));


        
        dp->backtrack_length = (uint32_t*)realloc(dp->backtrack_length, dp->snp_size * sizeof(uint32_t));

        dp->backtrack_size = snp_num*snp_num;
        dp->backtrack = (uint32_t*)realloc(dp->backtrack, dp->backtrack_size * sizeof(uint32_t));
    }

    dp->snp_num = snp_num;

    clear_SNP_IDs(&(dp->SNP_IDs));


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

    /****************************may have bugs********************************/
    memset(h->flag, 0, WINDOW_MAX_SIZE * sizeof(uint8_t));
    /****************************may have bugs********************************/


    // h->dp.max = NULL;
    // h->dp.colum = NULL;
    // h->dp.colum_len = NULL;
    // h->dp.matrix_size = 0;
    // h->dp.snp_num = 0;
    // h->dp.snp_size = 0;
    h->dp.snp_num = 0;
    h->dp.max = NULL;
    h->dp.max_for_sort = NULL;
    h->dp.visit = NULL;
    h->dp.snp_size = 0;
    h->dp.backtrack = NULL;
    h->dp.backtrack_size = 0;
    h->dp.backtrack_length = NULL;
    h->dp.buffer = NULL;
    h->dp.max_buffer = NULL;

    init_SNP_IDs(&(h->dp.SNP_IDs));

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
    free(h->dp.backtrack);
    free(h->dp.max);
    free(h->dp.max_for_sort);
    free(h->dp.visit);
    free(h->dp.backtrack_length);
    free(h->dp.buffer);
    free(h->dp.max_buffer);
    destory_SNP_IDs(&(h->dp.SNP_IDs));
}

inline void ResizeInitHaplotypeEvdience(haplotype_evdience_alloc* h)
{
    h->snp = 0;
    h->length = 0;
    h->sub_list_start = 0;
    h->sub_list_length = 0;
    /****************************may have bugs********************************/
    memset(h->flag, 0, WINDOW_MAX_SIZE * sizeof(uint8_t));
    /****************************may have bugs********************************/
}

inline void RsetInitHaplotypeEvdienceFlag(haplotype_evdience_alloc* h, long long beg, long long useful_length)
{
    /****************************may have bugs********************************/
    memset(h->flag + beg, 0, useful_length * sizeof(uint8_t));
    /****************************may have bugs********************************/
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

    /****************************may have bugs********************************/
    // char overlap_region[WINDOW + THRESHOLD*2 + 10];
    // char overlap_region_group[GROUP_SIZE][WINDOW + THRESHOLD*2 + 10];
    // char path[WINDOW + THRESHOLD*2 + 10];
    // Word matrix_bit[((WINDOW + 10)<<3)];

    char overlap_region[WINDOW_MAX_SIZE + THRESHOLD_MAX_SIZE*2 + 10];
    char overlap_region_group[GROUP_SIZE][WINDOW_MAX_SIZE + THRESHOLD_MAX_SIZE*2 + 10];
    char path[WINDOW_MAX_SIZE + THRESHOLD_MAX_SIZE*2 + 10];

    char path_fix[WINDOW_MAX_SIZE + THRESHOLD_MAX_SIZE*2 + 10];
    char overlap_region_fix[WINDOW_MAX_SIZE + THRESHOLD_MAX_SIZE*2 + 10];
    Word matrix_bit[((WINDOW_MAX_SIZE + 10)<<3)];
    /****************************may have bugs********************************/

    int path_length;
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
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read, Graph* g, Graph* DAGCon, 
                        Cigar_record* current_cigar, haplotype_evdience_alloc* hap,
                        Round2_alignment* second_round, int force_repeat, int is_consensus,
                        int* fully_cov, int* abnormal);
void init_Correct_dumy(Correct_dumy* list);
void destory_Correct_dumy(Correct_dumy* list);
void clear_Correct_dumy(Correct_dumy* list, overlap_region_alloc* overlap_list);
void clear_Correct_dumy_pure(Correct_dumy* list);
void get_seq_from_Graph(Graph* backbone, Graph* DAGCon, Correct_dumy* dumy, Cigar_record* current_cigar, char* self_string,
char* r_string, long long r_string_length, long long r_string_site);
void init_Cigar_record(Cigar_record* dummy);
void destory_Cigar_record(Cigar_record* dummy);
void clear_Cigar_record(Cigar_record* dummy);
void add_new_cell_to_cigar_record(Cigar_record* dummy, uint32_t len, uint32_t type);
void add_existing_cell_to_cigar_record(Cigar_record* dummy, uint32_t len, uint32_t type);
void add_new_cell_to_cigar_record_with_different_base(Cigar_record* dummy, uint32_t len, uint32_t type, char* seq);
void add_existing_cell_to_cigar_record_with_different_base(Cigar_record* dummy, uint32_t len, uint32_t type, char* seq);


/***
 type:
 0. match
 1. mismatch
 2. insertion
 3. deletion
 ***/
inline void add_cigar_record(char* seq, uint32_t len, Cigar_record* dummy, uint32_t type)
{    
    if(type == 0)///match
    {
        ///add to existing cell, just increase length
        if((uint32_t)dummy->current_operation == type)
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
        if((uint32_t)dummy->current_operation == type)
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
        if((uint32_t)dummy->current_operation == type)
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
        if((uint32_t)dummy->current_operation == type)
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

int verify_cigar_2(char* x, int x_len, char* y, int y_len, Cigar_record* cigar, int error);


int verify_single_window(long long x_start, long long x_end, 
long long overlap_x_s, long long overlap_y_s, int x_id,
int y_id, int y_strand, char* x_buffer, char* y_buffer,
All_reads* R_INF);

void init_Cigar_record_alloc(Cigar_record_alloc* x);
void resize_Cigar_record_alloc(Cigar_record_alloc* x, long long new_size);
void destory_Cigar_record_alloc(Cigar_record_alloc* x);

void afine_gap_alignment(const char *qseq, uint8_t* qnum, const int ql, 
const char *tseq, uint8_t* tnum, const int tl, const uint8_t *c2n, const int strand,
int sc_mch, int sc_mis, int gapo, int gape, int bandLen, int zdrop, int end_bonus,
long long* max_q_pos, long long* max_t_pos, long long* global_score, 
long long* extention_score, long long* q_boundary_score, long long* q_boundary_t_coordinate,
long long* t_boundary_score, long long* t_boundary_q_coordinate,
long long* droped, int mode);
void correct_overlap_high_het(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read);
long long get_affine_gap_score(overlap_region* ovc, UC_Read* g_read, UC_Read* overlap_read, uint8_t* x_num,
uint8_t* y_num, uint64_t EstimateXOlen, uint64_t EstimateYOlen);
int collect_hp_regions(overlap_region_alloc* olist, All_reads* R_INF, kvec_t_u8_warp* k_flag, float hp_rate, int rlen, FILE* fp);

inline int if_exact_match(char* x, long long xLen, char* y, long long yLen, long long xBeg, long long xEnd, long long yBeg, long long yEnd)
{
    long long overlapLen = xEnd - xBeg + 1;

    if(yEnd - yBeg + 1 == overlapLen)
    {
        if(memcmp(x+xBeg, y+yBeg, overlapLen)==0) return 1;
    }

    return 0;
}



#define FORWARD_KSW 0
#define BACKWARD_KSW 1
#define MATCH_SCORE_KSW 2
#define MISMATCH_SCORE_KSW 4
#define GAP_OPEN_KSW 4
#define GAP_EXT_KSW 2
#define Z_DROP_KSW 400
#define BAND_KSW 500
#endif
