#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "Hash_Table.h"
#include "Process_Read.h"
#include "Correct.h"
#include "CommandLines.h"
#include "kmer.h"
#include <pthread.h>
#include "ksort.h"
pthread_mutex_t output_mutex;

#define overlap_region_key(a) ((a).y_id)
KRADIX_SORT_INIT(overlap_region_sort, overlap_region, overlap_region_key, member_size(overlap_region, y_id))

void overlap_region_sort_y_id(overlap_region *a, long long n)
{
	radix_sort_overlap_region_sort(a, a + n);
}

void Init_Heap(HeapSq* HBT)
{
    HBT->MaxSize = 1000;
    HBT->heap = (ElemType*)malloc(HBT->MaxSize*sizeof(ElemType));
    HBT->index_i = (uint64_t*)malloc(HBT->MaxSize*sizeof(uint64_t));
    HBT->len = 0;
}

void destory_Heap(HeapSq* HBT)
{
    free(HBT->heap);
    free(HBT->index_i);
}

void clear_Heap(HeapSq* HBT)
{
    HBT->len = 0;
}



inline int cmp_ElemType(ElemType* x, ElemType* y)
{

    if (x->node.readID < y->node.readID)
    {
        return 1;
    }
    else if (x->node.readID > y->node.readID)
    {
        return 2;
    }
    else
    {   if (x->node.strand < y->node.strand)
        {
            return 1;
        }
        else if (x->node.strand > y->node.strand)
        {
            return 2;
        }
        else
        {

            if(x->node.offset < y->node.offset)
            {
                return 1;
            }
            else if(x->node.offset > y->node.offset)
            {
                return 2;
            }
            else
            {
                if (x->node.self_offset < y->node.self_offset)
                {
                    return 1;
                }
                else  if (x->node.self_offset > y->node.self_offset)
                {
                    return 2;
                }
                else  ///if both r_pos_x and self_offset are equal, offset must be equal
                {
                    return 0;
                }
                
            }   
        }    
    }
}


inline void Insert_Heap(HeapSq* HBT, ElemType* x)
{
    long long i, j;
    if (HBT->len == HBT->MaxSize) 
    {
        HBT->MaxSize = 2*HBT->MaxSize;
        HBT->heap = (ElemType*)realloc(HBT->heap, HBT->MaxSize*sizeof(ElemType));
        HBT->index_i = (uint64_t*)realloc(HBT->index_i,HBT->MaxSize*sizeof(uint64_t));
    }
    HBT->heap[HBT->len] = *x; //add element to tail
    HBT->len++; 
    i = HBT->len - 1; 
    while (i != 0)
    {
        j = (i - 1) / 2; 
        ///if (x >= HBT->heap[j]) 
        ///1: x<y
        if (cmp_ElemType(x, &HBT->heap[j])!=1)
            break;
        HBT->heap[i] = HBT->heap[j]; 
        i = j; 
    }
    HBT->heap[i] = *x;
}

inline int DeleteHeap(HeapSq* HBT, ElemType* get)
{
    ElemType temp, x;
    int i, j;
    if (HBT->len == 0)
    {
        return 0;
    }
    temp = HBT->heap[0]; 
    HBT->len--;
    if (HBT->len == 0) 
    {
        *get = temp;
        return 2;
    }

    x = HBT->heap[HBT->len]; 
    i = 0; 
    j = 2 * i + 1;
    while (j <= HBT->len - 1)
    {
        ///if (j < HBT->len - 1 && HBT->heap[j] > HBT->heap[j+1])
        if (j < HBT->len - 1 && cmp_ElemType(&HBT->heap[j], &HBT->heap[j + 1]) == 2)
            j++;
        ///if (x <= HBT->heap[j]) 
        if (cmp_ElemType(&x, &HBT->heap[j])!=2)
            break;
        HBT->heap[i] = HBT->heap[j];
        i = j; 
        j = 2 * i + 1;
    }
    HBT->heap[i] = x; 

    *get = temp;
    return 1;
}

void init_overlap_region_alloc(overlap_region_alloc* list)
{
    list->size = 1000;
    list->length = 0;
    ///list->list = (overlap_region*)malloc(sizeof(overlap_region)*list->size);
    list->list = (overlap_region*)calloc(list->size, sizeof(overlap_region));
    uint64_t i;
    for (i = 0; i < list->size; i++)
    { 
        init_fake_cigar(&(list->list[i].f_cigar));
        init_window_list_alloc(&(list->list[i].boundary_cigars));
    }
}
void clear_overlap_region_alloc(overlap_region_alloc* list)
{
    list->length = 0;
    list->mapped_overlaps_length = 0;
    uint64_t i = 0;
    for (i = 0; i < list->size; i++)
    {   
        list->list[i].w_list_length = 0;
        clear_fake_cigar(&(list->list[i].f_cigar));
        clear_window_list_alloc(&(list->list[i].boundary_cigars));
    }
}

void destory_overlap_region_alloc(overlap_region_alloc* list)
{
    uint64_t i = 0;
    for (i = 0; i < list->size; i++)
    {
        if (list->list[i].w_list_size != 0)
        {
            free(list->list[i].w_list);
        }
        destory_fake_cigar(&(list->list[i].f_cigar));
        destory_window_list_alloc(&(list->list[i].boundary_cigars));
    }
    free(list->list);
}


int get_fake_gap_pos(Fake_Cigar* x, int index)
{
    return (x->buffer[index]>>32);
}

int get_fake_gap_shift(Fake_Cigar* x, int index)
{
    uint32_t tmp = ((uint32_t)(x->buffer[index]));
    int result;
    if(tmp & ((uint32_t)1))
    {
        tmp = tmp >> 1;
        result = tmp;
        result = result * -1;
    }
    else
    {
        tmp = tmp >> 1;
        result = tmp;
    }
    
    return result;
}



int append_inexact_overlap_region_alloc(overlap_region_alloc* list, overlap_region* tmp, 
All_reads* R_INF, int add_beg_end)
{
   
    if (list->length + 1 > list->size)
    {
        list->size = list->size * 2;
        list->list = (overlap_region*)realloc(list->list, sizeof(overlap_region)*list->size);
        /// need to set new space to be 0
        memset(list->list + (list->size/2), 0, sizeof(overlap_region)*(list->size/2));
    }

    if (list->length!=0 && list->list[list->length - 1].y_id==tmp->y_id)
    {    
        ///if(list->list[list->length - 1].shared_seed >= tmp->shared_seed)
        if((list->list[list->length - 1].shared_seed > tmp->shared_seed)
           ||
           ((list->list[list->length - 1].shared_seed == tmp->shared_seed) && 
           (list->list[list->length - 1].overlapLen <= tmp->overlapLen)))
        {
            return 0;
        }
        else
        {
            list->length--;
        }
    }

    if(tmp->x_pos_s <= tmp->y_pos_s)
    {
        tmp->y_pos_s = tmp->y_pos_s - tmp->x_pos_s;
        tmp->x_pos_s = 0;
    }
    else
    {
        tmp->x_pos_s = tmp->x_pos_s - tmp->y_pos_s;
        tmp->y_pos_s = 0;
    }


    long long x_right_length = Get_READ_LENGTH((*R_INF), tmp->x_id) - tmp->x_pos_e - 1;
    long long y_right_length = Get_READ_LENGTH((*R_INF), tmp->y_id) - tmp->y_pos_e - 1;



    if(x_right_length <= y_right_length)
    {
        tmp->x_pos_e = Get_READ_LENGTH((*R_INF), tmp->x_id) - 1;
        tmp->y_pos_e = tmp->y_pos_e + x_right_length;        
    }
    else
    {
        tmp->x_pos_e = tmp->x_pos_e + y_right_length;
        tmp->y_pos_e = Get_READ_LENGTH((*R_INF), tmp->y_id) - 1;
    }

    
    if (tmp->x_pos_strand == 1)
    {
        list->list[list->length].x_id = tmp->x_id;
        list->list[list->length].x_pos_e = Get_READ_LENGTH((*R_INF), tmp->x_id) - tmp->x_pos_s - 1;
        list->list[list->length].x_pos_s = Get_READ_LENGTH((*R_INF), tmp->x_id) - tmp->x_pos_e - 1;
        list->list[list->length].x_pos_strand = 0;

        list->list[list->length].y_id = tmp->y_id;
        list->list[list->length].y_pos_e = Get_READ_LENGTH((*R_INF), tmp->y_id) - tmp->y_pos_s - 1;
        list->list[list->length].y_pos_s = Get_READ_LENGTH((*R_INF), tmp->y_id) - tmp->y_pos_e - 1;
        list->list[list->length].y_pos_strand = 1;




        resize_fake_cigar(&(list->list[list->length].f_cigar), (tmp->f_cigar.length + 2));
        if(add_beg_end == 1)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), list->list[list->length].x_pos_s, 0);
        }
        
        long long distance_gap;
        /****************************may have bugs********************************/
        ///long long pre_distance_gap = 0;
        long long pre_distance_gap = 0xfffffffffffffff;
        /****************************may have bugs********************************/
        long long i = 0;
        for (i = 0; i < (long long)tmp->f_cigar.length; i++)
        {
            distance_gap = get_fake_gap_shift(&(tmp->f_cigar), i);
            if(distance_gap != pre_distance_gap)
            {
                pre_distance_gap = distance_gap;
                add_fake_cigar(&(list->list[list->length].f_cigar), 
                Get_READ_LENGTH((*R_INF), tmp->x_id) - get_fake_gap_pos(&(tmp->f_cigar), i) - 1, 
                pre_distance_gap);
            }
        }

        if(add_beg_end == 1 && get_fake_gap_pos(&(list->list[list->length].f_cigar), 
            list->list[list->length].f_cigar.length - 1) != (long long)list->list[list->length].x_pos_e)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), 
            list->list[list->length].x_pos_e, 
            get_fake_gap_shift(&(list->list[list->length].f_cigar), 
            list->list[list->length].f_cigar.length - 1));
        }
    }
    else
    {
        list->list[list->length].x_id = tmp->x_id;
        list->list[list->length].x_pos_e = tmp->x_pos_e;
        list->list[list->length].x_pos_s = tmp->x_pos_s;
        list->list[list->length].x_pos_strand = tmp->x_pos_strand;

        list->list[list->length].y_id = tmp->y_id;
        list->list[list->length].y_pos_e = tmp->y_pos_e;
        list->list[list->length].y_pos_s = tmp->y_pos_s;
        list->list[list->length].y_pos_strand = tmp->y_pos_strand;



        resize_fake_cigar(&(list->list[list->length].f_cigar), (tmp->f_cigar.length + 2));
        if(add_beg_end == 1)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), list->list[list->length].x_pos_s, 0);
        }
        
        long long distance_self_pos = tmp->x_pos_e - tmp->x_pos_s;
        long long distance_pos = tmp->y_pos_e - tmp->y_pos_s;
        long long init_distance_gap = distance_pos - distance_self_pos;
        /****************************may have bugs********************************/
        ///long long pre_distance_gap = init_distance_gap;
        long long pre_distance_gap = 0xfffffffffffffff;
        /****************************may have bugs********************************/
        long long distance_gap;
        long long i = 0;
        for (i = tmp->f_cigar.length - 1; i >= 0; i--)
        {
            distance_gap = get_fake_gap_shift(&(tmp->f_cigar), i);
            if(distance_gap != pre_distance_gap)
            {
                pre_distance_gap = distance_gap;

                add_fake_cigar(&(list->list[list->length].f_cigar), 
                get_fake_gap_pos(&(tmp->f_cigar), i), init_distance_gap - pre_distance_gap);
            }
        }

        if(add_beg_end == 1 && get_fake_gap_pos(&(list->list[list->length].f_cigar), 
        list->list[list->length].f_cigar.length - 1) != (long long)list->list[list->length].x_pos_e)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), 
            list->list[list->length].x_pos_e, 
            get_fake_gap_shift(&(list->list[list->length].f_cigar), 
            list->list[list->length].f_cigar.length - 1));
        } 
    }

    list->list[list->length].shared_seed = tmp->shared_seed;
    list->list[list->length].align_length = 0;
    list->list[list->length].is_match = 0;
    list->list[list->length].non_homopolymer_errors = 0;
    list->list[list->length].strong = 0;

    list->length++;

    return 1;
}


void append_overlap_region_alloc_debug(overlap_region_alloc* list, overlap_region* tmp)
{
   
    if (list->length + 1 > list->size)
    {
        list->size = list->size * 2;
        list->list = (overlap_region*)realloc(list->list, sizeof(overlap_region)*list->size);
    }

    
    list->list[list->length].x_id = tmp->x_id;
    list->list[list->length].x_pos_e = tmp->x_pos_e;
    list->list[list->length].x_pos_s = tmp->x_pos_s;
    list->list[list->length].x_pos_strand = tmp->x_pos_strand;

    list->list[list->length].y_id = tmp->y_id;
    list->list[list->length].y_pos_e = tmp->y_pos_e;
    list->list[list->length].y_pos_s = tmp->y_pos_s;
    list->list[list->length].y_pos_strand = tmp->y_pos_strand;

    list->list[list->length].shared_seed = tmp->shared_seed;

    list->length++;
}


int cmp_by_x_pos_s(const void * a, const void * b)
{
    if ((*(overlap_region*)a).x_pos_s > (*(overlap_region*)b).x_pos_s)
    {
        return 1;
    }
    else if ((*(overlap_region*)a).x_pos_s < (*(overlap_region*)b).x_pos_s)
    {
        return -1;
    }
    else
    {

        if ((*(overlap_region*)a).x_pos_e > (*(overlap_region*)b).x_pos_e)
        {
            return 1;
        }
        else if ((*(overlap_region*)a).x_pos_e < (*(overlap_region*)b).x_pos_e)
        {
            return -1;
        }
        else
        {
            return 0;
        }
        
    }
}


int cmp_by_x_pos_e(const void * a, const void * b)
{
    if ((*(overlap_region*)a).x_pos_e > (*(overlap_region*)b).x_pos_e)
    {
        return 1;
    }
    else if ((*(overlap_region*)a).x_pos_e < (*(overlap_region*)b).x_pos_e)
    {
        return -1;
    }
    else
    {

        if ((*(overlap_region*)a).x_pos_s > (*(overlap_region*)b).x_pos_s)
        {
            return 1;
        }
        else if ((*(overlap_region*)a).x_pos_s < (*(overlap_region*)b).x_pos_s)
        {
            return -1;
        }
        else
        {
            return 0;
        }
        
    }
}


void debug_chain(k_mer_hit* a, long long a_n, Chain_Data* dp)
{
    long long i, j, current_j;
    long long selfLen, indels;
    long long distance_self_pos, distance_pos, distance_gap;
    for (i = 0; i < a_n; ++i) 
    {
        selfLen = indels = 0;
        j = i;
        while (j >= 0)
        {
            current_j = j;

            j = dp->pre[j];

            if(j != -1)
            {
                distance_self_pos = a[current_j].self_offset - a[j].self_offset;
                distance_pos = a[current_j].offset - a[j].offset;
                distance_gap = distance_pos > distance_self_pos? distance_pos - distance_self_pos : distance_self_pos - distance_pos;

                indels += distance_gap; 
                selfLen += distance_self_pos;
            }
        }

        if(indels != dp->indels[i])
        {
            fprintf(stderr, "indels: %lld, dp->indels[i]: %lld\n", 
            indels, dp->indels[i]);
        }

        if(selfLen != dp->self_length[i])
        {
            fprintf(stderr, "selfLen: %lld, dp->self_length[i]: %lld\n", 
            selfLen, dp->self_length[i]);
        }

    }
}

long long get_chainLen(long long x_beg, long long x_end, long long xLen, 
long long y_beg, long long y_end, long long yLen)
{
    if(x_beg <= y_beg)
    {
        y_beg = y_beg - x_beg;
        x_beg = 0;
    }
    else
    {
        x_beg = x_beg - y_beg;
        y_beg = 0;
    }

    long long x_right_length = xLen - x_end - 1;
    long long y_right_length = yLen - y_end - 1;


    if(x_right_length <= y_right_length)
    {
        x_end = xLen - 1;
        y_end = y_end + x_right_length;        
    }
    else
    {
        x_end = x_end + y_right_length;
        y_end = yLen - 1;
    }

    return x_end - x_beg + 1;
}


///double band_width_threshold = 0.05;
void chain_DP(k_mer_hit* a, long long a_n, Chain_Data* dp, overlap_region* result, 
double band_width_threshold, int max_skip, int x_readLen, int y_readLen)
{
    long long i, j;
    long long self_pos, pos, max_j, max_i, max_score, score, n_skip;
    long long distance_pos, distance_self_pos, distance_gap,  distance_min;
    ///double band_width_threshold = 0.05;
    double band_width_penalty = 1 / band_width_threshold;
    long long min_score = asm_opt.k_mer_length;
    long long max_indels, max_self_length;
    double gap_rate;
    long long total_indels, total_self_length;
    
    resize_Chain_Data(dp, a_n);
    // fill the score and backtrack arrays
	for (i = 0; i < a_n; ++i) 
    {
        pos = a[i].offset;
        self_pos = a[i].self_offset;
        max_j = -1;
        max_score = min_score;
        n_skip = 0;
        max_indels = 0;
        max_self_length = 0;
        

        ///may have a pre-cut condition for j
        for (j = i - 1; j >= 0; --j) 
        {
            distance_pos = pos - a[j].offset;
            distance_self_pos = self_pos - a[j].self_offset;
            ///a has been sorted by a[].offset
            ///note for a, we do not have any two elements that have both equal offsets and self_offsets
            ///but there maybe two elements that have equal offsets or equal self_offsets
            if(distance_pos == 0 || distance_self_pos <= 0)
            {
                continue;
            }

            distance_gap = distance_pos > distance_self_pos? distance_pos - distance_self_pos : distance_self_pos - distance_pos;
            
            total_indels = dp->indels[j] + distance_gap;
            total_self_length = dp->self_length[j] + distance_self_pos;
            if(total_indels > band_width_threshold * total_self_length)
            {
                continue;
            }

            ///min distance
            distance_min = distance_pos < distance_self_pos? distance_pos:distance_self_pos;
            score = distance_min < min_score? distance_min : min_score;

            gap_rate = (double)((double)(total_indels)/(double)(total_self_length));
            ///if the gap rate > 0.06, score will be negative
            score -= (long long)(gap_rate * score * band_width_penalty);

            score += dp->score[j];

            ///find a new max score
            if(score > max_score)
            {
                max_score = score;
                max_j = j;
                max_indels = total_indels;
                max_self_length = total_self_length;
                /****************************may have bugs********************************/
                n_skip = 0;
                /****************************may have bugs********************************/
            }/****************************may have bugs********************************/
            else
            {
                n_skip++;
                if(n_skip > max_skip)
                {
                    break;
                }
            }
            /****************************may have bugs********************************/
        }

        dp->score[i] = max_score;
        dp->pre[i] = max_j;
        dp->indels[i] = max_indels;
        dp->self_length[i] = max_self_length;
    }


    ///debug_chain(a, a_n, dp);


    

    max_score = -1;
    max_i = -1;
    long long mini_xLen = x_readLen * 2 + 2, tmp_xLen;
    for (i = 0; i < a_n; ++i) 
    {
        if(dp->score[i] > max_score)
        {
            max_score = dp->score[i];
            max_i = i;
            mini_xLen = get_chainLen(a[i].self_offset, a[i].self_offset, x_readLen, 
            a[i].offset, a[i].offset, y_readLen);
        }
        else if(dp->score[i] == max_score)
        {
            tmp_xLen = get_chainLen(a[i].self_offset, a[i].self_offset, x_readLen, 
            a[i].offset, a[i].offset, y_readLen);

            if(tmp_xLen < mini_xLen)
            {
                max_score = dp->score[i];
                max_i = i;

                mini_xLen = tmp_xLen;
            }
        }
        
    }


    clear_fake_cigar(&(result->f_cigar));
    ///note a has been sorted by offset, that means has been sorted by query offset
    i = max_i;
    result->x_pos_e = a[i].self_offset;
    result->y_pos_e = a[i].offset;
    result->shared_seed = max_score;
    result->overlapLen = mini_xLen;

    distance_self_pos = result->x_pos_e - a[i].self_offset;
    distance_pos = result->y_pos_e - a[i].offset;
    long long pre_distance_gap = distance_pos - distance_self_pos;
    ///record first site
    ///the length of f_cigar should be at least 1
    ///record the offset of reference
    add_fake_cigar(&(result->f_cigar), a[i].self_offset, pre_distance_gap);
    long long chainLen = 0;
    if(result->x_pos_strand == 1)
    {
        while (i >= 0)
        {
            distance_self_pos = result->x_pos_e - a[i].self_offset;
            distance_pos = result->y_pos_e - a[i].offset;
            distance_gap = distance_pos - distance_self_pos;
            if(distance_gap != pre_distance_gap)
            {
                pre_distance_gap = distance_gap;
                ///record this site
                add_fake_cigar(&(result->f_cigar), a[i].self_offset, pre_distance_gap);
            }

            chainLen++;
            result->x_pos_s = a[i].self_offset;
            result->y_pos_s = a[i].offset;
            i = dp->pre[i];
        }
    }
    else
    {
    
        while (i >= 0)
        {
            distance_self_pos = result->x_pos_e - a[i].self_offset;
            distance_pos = result->y_pos_e - a[i].offset;
            distance_gap = distance_pos - distance_self_pos;
            if(distance_gap == pre_distance_gap)
            {
                result->f_cigar.length--;
                add_fake_cigar(&(result->f_cigar), a[i].self_offset, pre_distance_gap);
            }
            else
            {
                pre_distance_gap = distance_gap;
                add_fake_cigar(&(result->f_cigar), a[i].self_offset, pre_distance_gap);
            }

            chainLen++;
            result->x_pos_s = a[i].self_offset;
            result->y_pos_s = a[i].offset;
            i = dp->pre[i];
        }
    }
}



void calculate_overlap_region_by_chaining(Candidates_list* candidates, overlap_region_alloc* overlap_list, 
uint64_t readID, uint64_t readLength, All_reads* R_INF, double band_width_threshold, int add_beg_end)
{
    overlap_region tmp_region;
    long long i = 0;
    uint64_t current_ID;
    uint64_t current_stand;

    if (candidates->length == 0)
    {
        return;
    }

    long long sub_region_beg;
    long long sub_region_end;

    init_fake_cigar(&(tmp_region.f_cigar));

    i = 0;
    while (i < candidates->length)
    {
        current_ID = candidates->list[i].readID;
        current_stand = candidates->list[i].strand;

        ///reference read
        tmp_region.x_id = readID;
        tmp_region.x_pos_strand = current_stand;
        ///query read
        tmp_region.y_id = current_ID;
        ///here the strand of query is always 0
        tmp_region.y_pos_strand = 0;  



        sub_region_beg = i;
        sub_region_end = i;
        i++;

        while (i < candidates->length 
        && 
        current_ID == candidates->list[i].readID
        &&
        current_stand == candidates->list[i].strand)
        {
            sub_region_end = i;
            i++;
        }

        if (tmp_region.x_id == tmp_region.y_id)
        {
            continue;
        }

    
        chain_DP(candidates->list + sub_region_beg, 
        sub_region_end - sub_region_beg + 1, &(candidates->chainDP), &tmp_region, band_width_threshold,
        50, Get_READ_LENGTH((*R_INF), tmp_region.x_id), Get_READ_LENGTH((*R_INF), tmp_region.y_id));

        // chain_DP_back(candidates->list + sub_region_beg, 
        // sub_region_end - sub_region_beg + 1, &(candidates->chainDP), &tmp_region, band_width_threshold);
        
        ///自己和自己重叠的要排除
        ///if (tmp_region.x_id != tmp_region.y_id && tmp_region.shared_seed > 1)
        if (tmp_region.x_id != tmp_region.y_id)
        {
            append_inexact_overlap_region_alloc(overlap_list, &tmp_region, R_INF, add_beg_end);
            ///append_inexact_overlap_region_alloc_back(overlap_list, &tmp_region, R_INF);
        }
    }

    destory_fake_cigar(&(tmp_region.f_cigar));

    qsort(overlap_list->list, overlap_list->length, sizeof(overlap_region), cmp_by_x_pos_s);
}



void append_window_list(overlap_region* region, uint64_t x_start, uint64_t x_end, int y_start, int y_end, int error,
int extra_begin, int extra_end, int error_threshold)
{
    
    long long length = region->x_pos_e - region->x_pos_s + 1;
    ///the length of window may large or small than WINDOW
    /****************************may have bugs********************************/
    uint64_t num_windows = length / WINDOW + 4;
    /****************************may have bugs********************************/

    ///w_list_length has alredy set to be 0 at clear_overlap_region_alloc
    if (num_windows > region->w_list_size)
    {
        region->w_list_size = num_windows;
        region->w_list = (window_list*)realloc(region->w_list, region->w_list_size*sizeof(window_list));
    }
    

    region->w_list[region->w_list_length].x_start = x_start;
    region->w_list[region->w_list_length].x_end = x_end;
    region->w_list[region->w_list_length].y_start = y_start;
    region->w_list[region->w_list_length].y_end = y_end;
    region->w_list[region->w_list_length].error = error;
    region->w_list[region->w_list_length].cigar.length = -1;
    region->w_list[region->w_list_length].extra_begin = extra_begin;
    region->w_list[region->w_list_length].extra_end = extra_end;
    region->w_list[region->w_list_length].error_threshold = error_threshold;
    region->w_list_length++;
}



void init_k_mer_pos_list_alloc(k_mer_pos_list_alloc* list)
{
    list->size = 1000;
    list->length = 0;
    //list->list = (k_mer_pos_list*)malloc(sizeof(k_mer_pos_list)*list->size);
    list->list = (k_mer_pos_list*)calloc(list->size, sizeof(k_mer_pos_list)); 
}

void clear_k_mer_pos_list_alloc(k_mer_pos_list_alloc* list)
{
    list->length = 0;
}

void destory_k_mer_pos_list_alloc(k_mer_pos_list_alloc* list)
{
    free(list->list);
}

void append_k_mer_pos_list_alloc(k_mer_pos_list_alloc* list, k_mer_pos* n_list, uint64_t n_length, 
uint64_t n_end_pos, uint8_t n_direction)
{
   
    if (list->length + 1 > list->size)
    {
        list->size = list->size * 2;
        list->list = (k_mer_pos_list*)realloc(list->list, sizeof(k_mer_pos_list)*list->size);
    }
    
    list->list[list->length].list = n_list;
    list->list[list->length].length = n_length;
    list->list[list->length].direction = n_direction;
    list->list[list->length].end_pos = n_end_pos;

    list->length++;
}




int cmp_k_mer_pos_list(const void * a, const void * b)
{
    if ((*(k_mer_pos_list*)a).length > (*(k_mer_pos_list*)b).length)
    {
        return 1;
    }
    else if ((*(k_mer_pos_list*)a).length < (*(k_mer_pos_list*)b).length)
    {
        return -1;
    }
    else
    {
        return 0;
    }
    
}

inline void append_pos_to_Candidates_list(Candidates_list* candidates, ElemType* x)
{
    candidates->list[candidates->length] = x->node;
    candidates->length++;
}




void test_single_list(Candidates_list* candidates, k_mer_pos* n_list, uint64_t n_lengh, uint64_t end_pos, uint64_t strand)
{
    uint64_t i;
    long long j = 0;
    for (i = 0; i < n_lengh; i++)
    {
        
        for (; j < candidates->length; j++)
        {
            if (
                n_list[i].offset == (uint64_t)candidates->list[j].offset
                &&
                n_list[i].readID == candidates->list[j].readID
                &&
                end_pos == (uint64_t)candidates->list[j].self_offset
                &&
                strand == candidates->list[j].strand
            )
            {
                break;
            }   
        }

        if (j == candidates->length)
        {
            fprintf(stderr, "ERROR 4\n");
        }
    }

}


void merge_k_mer_pos_list_alloc_heap_sort(k_mer_pos_list_alloc* list, Candidates_list* candidates, HeapSq* HBT)
{
    clear_Heap(HBT);

    uint64_t total_length = 0;
    uint64_t i;
    ElemType x, y;
    //add the first element of each list to stack
    for (i = 0; i < list->length; i++)
    {
        x.ID = i;
        x.node.offset = list->list[i].list[0].offset;
        x.node.readID = list->list[i].list[0].readID;
        x.node.self_offset = list->list[i].end_pos;
        x.node.strand = list->list[i].direction;

        Insert_Heap(HBT, &x);

        HBT->index_i[i] = 1;
        
        total_length = total_length + list->list[i].length;
    }


    candidates->length = 0;
    if(total_length > (uint64_t)candidates->size)
    {
        candidates->size = total_length;
        candidates->list = (k_mer_hit*)realloc(candidates->list, sizeof(k_mer_hit)*candidates->size);
        candidates->tmp = (k_mer_hit*)realloc(candidates->tmp, sizeof(k_mer_hit)*candidates->size);
    }

    uint64_t ID;
    int flag;
    ///while (flag = DeleteHeap(HBT, &x))
    while ((flag = DeleteHeap(HBT, &x)))
    {
        append_pos_to_Candidates_list(candidates, &x);
        
        i = HBT->index_i[x.ID];
        ID = x.ID;

        if (flag == 2)
        {
            for (; i < list->list[x.ID].length; i++)
            {
                y.ID = ID;
                y.node.offset = list->list[x.ID].list[i].offset;
                y.node.readID = list->list[x.ID].list[i].readID;
                y.node.self_offset = list->list[x.ID].end_pos;
                y.node.strand = list->list[x.ID].direction;
                append_pos_to_Candidates_list(candidates, &y);
            }
            
            break;
        }

        if (i < list->list[x.ID].length)
        {
            y.ID = ID;
            y.node.offset = list->list[x.ID].list[i].offset;
            y.node.readID = list->list[x.ID].list[i].readID;
            y.node.self_offset = list->list[x.ID].end_pos;
            y.node.strand = list->list[x.ID].direction;
            Insert_Heap(HBT, &y);
            HBT->index_i[ID]++;
        }
    }     
}



void init_Count_Table(Count_Table** table)
{
    *table = kh_init(COUNT64);
}

void init_Pos_Table(Pos_Table** table)
{
    *table = kh_init(POS64);
}

void init_Total_Count_Table(int k, Total_Count_Table* TCB)
{
    if(k>64)
    {
        fprintf(stderr, "k-mer is too long. The length of k-mer must <= 64.");
        fflush(stderr);
        exit(0);
    }

    int total_bits = k * 2; 
    TCB->prefix_bits = PREFIX_BITS;
    TCB->suffix_bits = total_bits - TCB->prefix_bits;
    if (TCB->suffix_bits > MAX_SUFFIX_BITS)
    {
        TCB->suffix_bits = MAX_SUFFIX_BITS;
        TCB->prefix_bits = total_bits - TCB->suffix_bits;
    }
    ///TCB->suffix_mode = (1ULL<<TCB->suffix_bits) - 1;
    ///right shift is safe, since TCB->suffix_bits cannot be 0
    TCB->suffix_mode = ALL >> (64 - TCB->suffix_bits);
    
    ///number of small hash table
    TCB->size = (1ULL<<TCB->prefix_bits);
    TCB->sub_h = (Count_Table**)malloc(sizeof(Count_Table*)*TCB->size);
    TCB->sub_h_lock = (Hash_table_spin_lock*)malloc(sizeof(Hash_table_spin_lock)*TCB->size);
    memset(TCB->sub_h_lock, 0, sizeof(Hash_table_spin_lock)*TCB->size);

    int i = 0;
    for (i = 0; i < TCB->size; i++)
    {
        init_Count_Table(&(TCB->sub_h[i]));
        TCB->sub_h_lock[i].lock = 0;
    }
    TCB->non_unique_k_mer = 0;
}



void init_Total_Pos_Table(Total_Pos_Table* TCB, Total_Count_Table* pre_TCB)
{

    TCB->prefix_bits = pre_TCB->prefix_bits;
    TCB->suffix_bits = pre_TCB->suffix_bits;
    TCB->suffix_mode = pre_TCB->suffix_mode;
    TCB->size = pre_TCB->size;
    TCB->useful_k_mer = 0;
    TCB->total_occ = 0;
    TCB->k_mer_index = NULL;
    TCB->sub_h_lock = (Hash_table_spin_lock*)malloc(sizeof(Hash_table_spin_lock)*TCB->size);
    memset(TCB->sub_h_lock, 0, sizeof(Hash_table_spin_lock)*TCB->size);
    TCB->sub_h = (Pos_Table**)malloc(sizeof(Pos_Table*)*TCB->size);
    TCB->pos = NULL;

    int i = 0;
    for (i = 0; i < TCB->size; i++)
    {
        init_Pos_Table(&(TCB->sub_h[i]));
        TCB->sub_h_lock[i].lock = 0;
    }
}


void destory_Total_Count_Table(Total_Count_Table* TCB)
{
    int i;
    for (i = 0; i < TCB->size; i++)
    {
        kh_destroy(COUNT64, TCB->sub_h[i]);
    }
    free(TCB->sub_h);
    free(TCB->sub_h_lock);
}


void destory_Total_Pos_Table(Total_Pos_Table* TCB)
{
    free(TCB->k_mer_index);
    free(TCB->sub_h_lock);
    free(TCB->pos);

    int i;
    for (i = 0; i < TCB->size; i++)
    {
        kh_destroy(POS64, TCB->sub_h[i]);
    }
    free(TCB->sub_h);
}


void write_Total_Pos_Table(Total_Pos_Table* TCB, char* read_file_name)
{
    fprintf(stderr, "Writing index to disk... \n");
    char* index_name = (char*)malloc(strlen(read_file_name)+5);
    sprintf(index_name, "%s.idx", read_file_name);
    FILE* fp = fopen(index_name, "w");
    fwrite(&asm_opt.adapterLen, sizeof(asm_opt.adapterLen), 1, fp);
    fwrite(&asm_opt.k_mer_min_freq, sizeof(asm_opt.k_mer_min_freq), 1, fp);
    fwrite(&asm_opt.k_mer_max_freq, sizeof(asm_opt.k_mer_max_freq), 1, fp);
    fwrite(&TCB->prefix_bits, sizeof(TCB->prefix_bits), 1, fp);
    fwrite(&TCB->suffix_bits, sizeof(TCB->suffix_bits), 1, fp);
    fwrite(&TCB->suffix_mode, sizeof(TCB->suffix_mode), 1, fp);
    fwrite(&TCB->size, sizeof(TCB->size), 1, fp);
    fwrite(&TCB->useful_k_mer, sizeof(TCB->useful_k_mer), 1, fp);
    fwrite(&TCB->total_occ, sizeof(TCB->total_occ), 1, fp);
    fwrite(TCB->k_mer_index, sizeof(uint64_t), TCB->useful_k_mer+1, fp);
    fwrite(TCB->pos, sizeof(k_mer_pos), TCB->total_occ, fp);


    int i;
    for (i = 0; i < TCB->size; i++)
    {
        kh_write(POS64, TCB->sub_h[i], fp);
    }

    free(index_name);    
    fclose(fp);
    fprintf(stderr, "Index has been written.\n");
}


int load_Total_Pos_Table(Total_Pos_Table* TCB, char* read_file_name)
{
    fprintf(stderr, "Loading index from disk... \n");
    char* index_name = (char*)malloc(strlen(read_file_name)+5);
    sprintf(index_name, "%s.idx", read_file_name);
    FILE* fp = fopen(index_name, "r");
    if (!fp)
    {
        return 0;
    }
    int f_flag;
    int local_adapterLen;
    f_flag = fread(&local_adapterLen, sizeof(local_adapterLen), 1, fp);
    if(local_adapterLen != asm_opt.adapterLen)
    {
        fprintf(stderr, "the adapterLen of index is: %d, but the adapterLen set by user is: %d\n", 
        local_adapterLen, asm_opt.adapterLen);
        exit(1);
    }
    f_flag += fread(&asm_opt.k_mer_min_freq, sizeof(asm_opt.k_mer_min_freq), 1, fp);
    f_flag += fread(&asm_opt.k_mer_max_freq, sizeof(asm_opt.k_mer_max_freq), 1, fp);
    f_flag += fread(&TCB->prefix_bits, sizeof(TCB->prefix_bits), 1, fp);
    f_flag += fread(&TCB->suffix_bits, sizeof(TCB->suffix_bits), 1, fp);
    f_flag += fread(&TCB->suffix_mode, sizeof(TCB->suffix_mode), 1, fp);
    f_flag += fread(&TCB->size, sizeof(TCB->size), 1, fp);
    f_flag += fread(&TCB->useful_k_mer, sizeof(TCB->useful_k_mer), 1, fp);
    f_flag += fread(&TCB->total_occ, sizeof(TCB->total_occ), 1, fp);

    if (TCB->useful_k_mer+1)
    {
        TCB->k_mer_index = (uint64_t*)malloc(sizeof(uint64_t)*(TCB->useful_k_mer+1));
        f_flag += fread(TCB->k_mer_index, sizeof(uint64_t), TCB->useful_k_mer+1, fp);
    }
    else
    {
        TCB->k_mer_index = NULL;
    }
    
    
    if (TCB->total_occ)
    {
        TCB->pos = (k_mer_pos*)malloc(sizeof(k_mer_pos)*TCB->total_occ);
        f_flag += fread(TCB->pos, sizeof(k_mer_pos), TCB->total_occ, fp);
    }
    else
    {
        TCB->pos = NULL;
    }
    
    TCB->sub_h_lock = (Hash_table_spin_lock*)malloc(sizeof(Hash_table_spin_lock)*TCB->size);
    memset(TCB->sub_h_lock, 0, sizeof(Hash_table_spin_lock)*TCB->size);

    TCB->sub_h = (Pos_Table**)malloc(sizeof(Pos_Table*)*TCB->size);

    int i;
    for (i = 0; i < TCB->size; i++)
    {
        init_Pos_Table(&(TCB->sub_h[i]));
        TCB->sub_h_lock[i].lock = 0;
        kh_load(POS64, TCB->sub_h[i], fp);
    }

    free(index_name);    
    fclose(fp);
    fprintf(stderr, "Index has been loaded.\n");

    return 1;
}

typedef struct
{
    long long* list;
    uint64_t length;
} H_peaks;

void insert_H_peaks(H_peaks* h, long long index, long long value)
{
    if(h->length <= (uint64_t)index)
    {
        long long newLen = index + 1;
        h->list = (long long*)realloc(h->list, newLen*sizeof(long long));
        memset(h->list + h->length, 0, sizeof(long long) * (newLen - h->length));
        h->length = newLen;
    }

    h->list[index] += value;
}

inline void RC_Hash_code(Hash_code* code, Hash_code* rc_code, int k)
{
    rc_code->x[0] = 0;
    rc_code->x[1] = 0;
    int i;
    for (i = 0; i < k; i++)
    {
        rc_code->x[0] = rc_code->x[0] << 1;
        rc_code->x[1] = rc_code->x[1] << 1;
        rc_code->x[0] |= (((uint64_t)((code->x[0] >> i) & 1))^((uint64_t)1));
        rc_code->x[1] |= (((uint64_t)((code->x[1] >> i) & 1))^((uint64_t)1));
    }
}



void get_peak_debug(Total_Count_Table* TCB, long long* min, long long* max)
{
    int i;
    Count_Table* h;
    khint_t k;
    long long c_count;
    H_peaks LH;
    LH.list = NULL;
    LH.length = 0;
    uint64_t sub_ID;
    uint64_t sub_key;
    Hash_code code, rc_code, debug_code;
    char str[100];
    char rc_str[100];


    
    for (i = 0; i < TCB->size; i++)
    {
        h = TCB->sub_h[i];
        for (k = kh_begin(h); k != kh_end(h); ++k)
        {
            if (kh_exist(h, k))            // test if a bucket contains data
            {
                sub_ID = i;
                sub_key = kh_key(h, k);

                recover_hash_code(sub_ID, sub_key, &code, TCB->suffix_mode, 
                TCB->suffix_bits, asm_opt.k_mer_length);
                RC_Hash_code(&code, &rc_code, asm_opt.k_mer_length);
                RC_Hash_code(&rc_code, &debug_code, asm_opt.k_mer_length);
                if(code.x[0] != debug_code.x[0] || code.x[1] != debug_code.x[1])
                {
                    fprintf(stderr, "error\n");
                }

                Hashcode_to_string(&code, str, asm_opt.k_mer_length);
                Hashcode_to_string(&rc_code, rc_str, asm_opt.k_mer_length);
                reverse_complement(str, asm_opt.k_mer_length);
                if(memcmp(str, rc_str, asm_opt.k_mer_length) != 0)
                {
                    fprintf(stderr, "error\n");
                    int j;
                    for (j = 0; j < asm_opt.k_mer_length; j++)
                    {
                        fprintf(stderr, "%c",str[j]);
                    }
                    fprintf(stderr, "\n");

                    for (j = 0; j < asm_opt.k_mer_length; j++)
                    {
                        fprintf(stderr, "%c",rc_str[j]);
                    }
                    fprintf(stderr, "\n");
                    
                }


                ///get_Total_Count_Table(&TCB, &k_code, k_mer_length);

                c_count = kh_value(h, k);

                if(get_Total_Count_Table(TCB, &code, asm_opt.k_mer_length) != c_count)
                {
                    fprintf(stderr, "error\n");
                }

                insert_H_peaks(&LH, c_count, c_count);
            }
        }
    }

    (*max) = -1;
    (*min) = -1;
    long long max_value = -1;
    for (i = 0; i < (long long)LH.length; i++)
    {
        if(LH.list[i] >= max_value)
        {
            max_value = LH.list[i];
            (*max) = i;
        }
    }

    long long min_value = max_value;
    for (i = 0; i < (long long)LH.length; i++)
    {
        if(LH.list[i] < min_value && LH.list[i] != 0)
        {
            min_value = LH.list[i];
            (*min) = i;
        }
    }

    for (i = 0; i < (long long)LH.length; i++)
    {
        ///fprintf(stderr, "%d, %d\n", i, LH.list[i]);
        fprintf(stderr, "%lld\n", LH.list[i]);
    }
    
    

    free(LH.list);
}

///1: a > b; -1: a < b; 0: a=b
int cmp_Hash_code(Hash_code* a, Hash_code* b)
{
    if(a->x[1] > b->x[1])
    {
        return 1;
    }
    if(a->x[1] < b->x[1])
    {
        return -1;
    }
    ///a->x[1] == b->x[1]
    if(a->x[0] > b->x[0])
    {
        return 1;
    }
    if(a->x[0] < b->x[0])
    {
        return -1;
    }

    return 0;
}

int get_total_freq(Total_Count_Table* TCB, uint64_t sub_ID, uint64_t sub_key, long long* T_count)
{
    Hash_code code, rc_code;
    long long count, rc_count;

    recover_hash_code(sub_ID, sub_key, &code, TCB->suffix_mode, 
                TCB->suffix_bits, asm_opt.k_mer_length);
    RC_Hash_code(&code, &rc_code, asm_opt.k_mer_length);

    count = get_Total_Count_Table(TCB, &code, asm_opt.k_mer_length);
    rc_count = get_Total_Count_Table(TCB, &rc_code, asm_opt.k_mer_length);
    (*T_count) = count + rc_count;

    if(count == 0)
    {
        return 0;
    }///count > 0 && rc_count == 0
    else if(rc_count == 0)
    {
        return 1;
    }///count > 0 && rc_count > 0
    else
    {
        int flag = cmp_Hash_code(&code, &rc_code);

        ///code > rc_code
        if(flag > 0)
        {
            return 1;
        }///code < rc_code
        else if(flag < 0)
        {
            return 0;
        }
        else
        {
            (*T_count) = (*T_count)/2;
            return 1;
        }
    }
}

void get_peak(Total_Count_Table* TCB, long long* min, long long* max, long long* up_boundary)
{
    int i;
    Count_Table* h;
    khint_t k;
    long long count;
    H_peaks LH;
    LH.list = NULL;
    LH.length = 0;
    uint64_t sub_ID;
    uint64_t sub_key;
    
    
    
    for (i = 0; i < TCB->size; i++)
    {
        h = TCB->sub_h[i];
        for (k = kh_begin(h); k != kh_end(h); ++k)
        {
            if (kh_exist(h, k))            // test if a bucket contains data
            {
                sub_ID = i;
                sub_key = kh_key(h, k);

                if(get_total_freq(TCB, sub_ID, sub_key, &count)==1)
                {
                    insert_H_peaks(&LH, count, count);
                }
            }
        }
    }

    (*max) = -1;
    (*min) = -1;
    long long max_value = -1;
    //// seed with freq 1 is useless
    for (i = 2; i < (long long)LH.length; i++)
    {
        if(LH.list[i] >= max_value)
        {
            max_value = LH.list[i];
            (*max) = i;
        }
    }

    long long opt = 4;
    (*up_boundary) = -1;
    for (i = (*max) + opt; i < (long long)LH.length; i++)
    {
        if(LH.list[i] > LH.list[i-opt])
        {
            long long j = i-opt;
            for (; j < i; j++)
            {
                if(LH.list[j] < LH.list[j+1])
                {
                    (*up_boundary) = j;
                    goto end_opt;
                }
            }
            
            (*up_boundary) = i;
            goto end_opt;
        }
    }

    end_opt:
    if((*up_boundary) == -1 || (*up_boundary) > (*max) * 10)
    {
        (*up_boundary) = (*max) * 10;
    }




    long long min_value = max_value;
    //// seed with freq 1 is useless
    for (i = 2; i < (long long)LH.length && i < (*max); i++)
    {
        if(LH.list[i] < min_value && LH.list[i] != 0)
        {
            min_value = LH.list[i];
            (*min) = i;
        }
    }

    free(LH.list);
}



void Traverse_Counting_Table(Total_Count_Table* TCB, Total_Pos_Table* PCB, int k_mer_min_freq, int k_mer_max_freq)
{
    int i;
    Count_Table* h;
    khint_t k;
    uint64_t sub_key;
    uint64_t sub_ID;
    PCB->useful_k_mer = 0;
    PCB->total_occ = 0;

    long long freq_min, max, freq_up;
    ///get_peak_debug(TCB, &freq_min, &freq_max);
    get_peak(TCB, &freq_min, &max, &freq_up);
    // fprintf(stdout, "freq_min: %d, freq_max: %d, freq_up:%d\n", 
    // freq_min, max, freq_up);
    if(freq_min < k_mer_min_freq)
    {
        k_mer_min_freq = freq_min;
    }
    if(freq_up > k_mer_max_freq)
    {
        k_mer_max_freq = freq_up;
    }

    // fprintf(stdout, "k_mer_min_freq: %d, k_mer_max_freq: %d\n", 
    // k_mer_min_freq, k_mer_max_freq);


    khint_t t;  
    int absent;
    long long count;

    /********************************************
     hash_table(key) ----> PCB->k_mer_index ------> PCB->pos
     ********************************************/
    for (i = 0; i < TCB->size; i++)
    {
        h = TCB->sub_h[i];
        for (k = kh_begin(h); k != kh_end(h); ++k)
        {
            if (kh_exist(h, k))            // test if a bucket contains data
            {
                sub_ID = i;
                sub_key = kh_key(h, k);
                get_total_freq(TCB, sub_ID, sub_key, &count);

                if (count>=k_mer_min_freq && count<=k_mer_max_freq)
                {
                    t = kh_put(POS64, PCB->sub_h[sub_ID], sub_key, &absent);
    
                    if (absent)
                    {
                        ///kh_value(PCB->sub_h[sub_ID], t) = useful_k_mer + total_occ;
                        kh_value(PCB->sub_h[sub_ID], t) = PCB->useful_k_mer;
                    }
                    else   
                    {
                        ///kh_value(PCB->sub_h[sub_ID], t)++;
                        fprintf(stderr, "ERROR\n");
                    }

                   
                    PCB->useful_k_mer++;
                    PCB->total_occ = PCB->total_occ + kh_value(h, k);
                }
            }
        }
    }


    // fprintf(stdout, "useful_k_mer: %lld\n",PCB->useful_k_mer);
    // fprintf(stdout, "total_occ: %lld\n",PCB->total_occ);

    PCB->k_mer_index = (uint64_t*)malloc(sizeof(uint64_t)*(PCB->useful_k_mer+1));

    PCB->k_mer_index[0] = 0;

    PCB->total_occ = 0;
    PCB->useful_k_mer = 0;

    for (i = 0; i < TCB->size; i++)
    {
        h = TCB->sub_h[i];
        for (k = kh_begin(h); k != kh_end(h); ++k)
        {
            if (kh_exist(h, k))            // test if a bucket contains data
            {
                sub_ID = i;
                sub_key = kh_key(h, k);
                get_total_freq(TCB, sub_ID, sub_key, &count);

                ///if (kh_value(h, k)>=k_mer_min_freq && kh_value(h, k)<=k_mer_max_freq)
                if (count>=k_mer_min_freq && count<=k_mer_max_freq)
                {
                    PCB->useful_k_mer++;
                    PCB->total_occ = PCB->total_occ + kh_value(h, k);
                    PCB->k_mer_index[PCB->useful_k_mer] = PCB->total_occ;
                }
            }
        }
    }

    PCB->pos = (k_mer_pos*)malloc(sizeof(k_mer_pos)*PCB->total_occ);
    memset(PCB->pos, 0, sizeof(k_mer_pos)*PCB->total_occ);
    
    
}




int cmp_k_mer_pos(const void * a, const void * b)
{
    if ((*(k_mer_pos*)a).readID != (*(k_mer_pos*)b).readID)
    {
        return (*(k_mer_pos*)a).readID > (*(k_mer_pos*)b).readID ? 1 : -1; 
    }
    else
    {
        if ((*(k_mer_pos*)a).offset != (*(k_mer_pos*)b).offset)
        {
            return (*(k_mer_pos*)a).offset > (*(k_mer_pos*)b).offset ? 1 : -1; 
        }
        else
        {
            return 0;
        }
    }
}


void init_Chain_Data(Chain_Data* x)
{
    x->length = 0;
    x->size = 0;
    x->score = NULL;
    x->pre = NULL;
    x->indels = NULL;
    x->self_length = NULL;
}

void clear_Chain_Data(Chain_Data* x)
{
    x->length = 0;
}


void destory_Chain_Data(Chain_Data* x)
{
    free(x->score);
    free(x->pre);
    free(x->indels);
    free(x->self_length);
}


void resize_Chain_Data(Chain_Data* x, long long size)
{
    if(size > x->size)
    {
        x->size = size;
        x->score = (long long*)realloc(x->score, x->size*sizeof(long long));   
        x->pre = (long long*)realloc(x->pre, x->size*sizeof(long long));
        x->indels = (long long*)realloc(x->indels, x->size*sizeof(long long));
        x->self_length = (long long*)realloc(x->self_length, x->size*sizeof(long long));
    }
}


void init_Candidates_list(Candidates_list* l)
{
    l->length = 0;
    l->size = 0;
    l->list = NULL;
    l->tmp = NULL;
    l->foward_pos = 0;
    l->rc_pos = 0;
    init_Chain_Data(&(l->chainDP));
}

void clear_Candidates_list(Candidates_list* l)
{
    l->length = 0;
    l->foward_pos = 0;
    l->rc_pos = 0;
    clear_Chain_Data(&(l->chainDP));
}

void destory_Candidates_list(Candidates_list* l)
{
    free(l->list);
    free(l->tmp);
    destory_Chain_Data(&(l->chainDP));
}



int cmp_candidates_list(const void * a, const void * b)
{
    long long r_pos_a, r_pos_b;

    if((*((k_mer_hit*)a)).strand != (*((k_mer_hit*)b)).strand)
    {
        return (*((k_mer_hit*)a)).strand > (*((k_mer_hit*)b)).strand ? 1: -1;
    }
    else
    {
        if((*((k_mer_hit*)a)).readID != (*((k_mer_hit*)b)).readID)
        {
            return (*((k_mer_hit*)a)).readID > (*((k_mer_hit*)b)).readID ? 1: -1;
        }
        else
        {
            r_pos_a = (*((k_mer_hit*)a)).offset - (*((k_mer_hit*)a)).self_offset;
            r_pos_b = (*((k_mer_hit*)b)).offset - (*((k_mer_hit*)b)).self_offset;

            if(r_pos_a != r_pos_b)
            {
                return r_pos_a > r_pos_b ? 1: -1;
            }
            else
            {
                if((*((k_mer_hit*)a)).self_offset != (*((k_mer_hit*)b)).self_offset)
                {
                    return (*((k_mer_hit*)a)).self_offset > (*((k_mer_hit*)b)).self_offset ? 1: -1;
                }
                else
                {
                    return 0;
                }
                
            }
            
        }
        
    }
}




void init_fake_cigar(Fake_Cigar* x)
{
    x->buffer = NULL;
    x->length = 0;
    x->size = 0;
}

void destory_fake_cigar(Fake_Cigar* x)
{
    if(x->size > 0)
    {
        free(x->buffer);
    }
}

void clear_fake_cigar(Fake_Cigar* x)
{
    x->length = 0;
}

void add_fake_cigar(Fake_Cigar* x, uint32_t gap_site, int32_t gap_shift)
{
    if(x->length + 1 > x->size)
    {
        x->size = (x->length + 1) * 2;
        x->buffer = (uint64_t*)realloc(x->buffer, sizeof(uint64_t) * x->size);
    }

    x->buffer[x->length] = gap_site;
    x->buffer[x->length] = x->buffer[x->length] << 32;

    if(gap_shift < 0)
    {
        gap_shift = gap_shift * -1;
        gap_site = gap_shift;
        gap_site = gap_site << 1;
        gap_site = gap_site | ((uint32_t)1);
    }
    else
    {
        gap_site = gap_shift;
        gap_site = gap_site << 1;
    }

    x->buffer[x->length] = x->buffer[x->length] | ((uint32_t)gap_site);

    x->length++;
}


void resize_fake_cigar(Fake_Cigar* x, uint64_t size)
{
    if(size > x->size)
    {
        x->size = size;
        x->buffer = (uint64_t*)realloc(x->buffer, sizeof(uint64_t) * x->size);
    }

    x->length = 0;
}


void init_window_list_alloc(window_list_alloc* x)
{
    x->buffer = NULL;
    x->length = 0;
    x->size = 0;
}

void clear_window_list_alloc(window_list_alloc* x)
{
    x->length = 0;
}

void destory_window_list_alloc(window_list_alloc* x)
{
    if(x->size != 0)
    {
        free((x->buffer));
    }
}

void resize_window_list_alloc(window_list_alloc* x, long long size)
{
    if(size > x->size)
    {
        x->size = size;
        x->buffer = (window_list*)realloc(x->buffer, sizeof(window_list) * x->size);
    }

    long long i;
    for (i = 0; i < x->size; i++)
    {
        x->buffer[i].error = -1; 
    }
    

    x->length = 0;
}