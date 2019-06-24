#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "Hash_Table.h"
#include "Process_Read.h"
#include <pthread.h>
pthread_mutex_t output_mutex;


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

inline int cmp_ElemType_back(ElemType* x, ElemType* y)
{
    if (x->node.strand < y->node.strand)
    {
        return 1;
    }
    else if (x->node.strand > y->node.strand)
    {
        return 2;
    }
    else
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
        {
            if (x->node.offset < y->node.offset)
            {
                return 1;
            }
            else if (x->node.offset > y->node.offset)
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
                else
                {
                    return 0;
                }
                
            }   
        }    
    }
}


inline int cmp_ElemType(ElemType* x, ElemType* y)
{
    long long r_pos_x, r_pos_y;
    if (x->node.strand < y->node.strand)
    {
        return 1;
    }
    else if (x->node.strand > y->node.strand)
    {
        return 2;
    }
    else
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
        {

            r_pos_x = x->node.offset - x->node.self_offset;
            r_pos_y = y->node.offset - y->node.self_offset;

            ///if (x->node.offset < y->node.offset)
            if (r_pos_x < r_pos_y)
            {
                return 1;
            }
            ///else if (x->node.offset > y->node.offset)
            else if (r_pos_x > r_pos_y)
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
                else  ///如果r_pos_x和self_offset都相等，那么offset肯定也相等，也没必要再比了
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
    if (HBT->len == HBT->MaxSize) //若堆满，将数组空间扩展为原来的2倍
    {
        HBT->MaxSize = 2*HBT->MaxSize;
        HBT->heap = (ElemType*)realloc(HBT->heap, HBT->MaxSize*sizeof(ElemType));
        HBT->index_i = (uint64_t*)realloc(HBT->index_i,HBT->MaxSize*sizeof(uint64_t));
    }
    HBT->heap[HBT->len] = *x; //向堆尾添加新元素
    HBT->len++; //堆长度加1
    i = HBT->len - 1; //i指向待调整元素的位置，即其数组下标，初始指向新元素所在的堆尾位置
    while (i != 0)
    {
        j = (i - 1) / 2; //j指向下标为i的元素的双亲
        ///if (x >= HBT->heap[j]) //若新元素大于待调整元素的双亲，则比较调整结束，退出循环
        ///1: x<y
        if (cmp_ElemType(x, &HBT->heap[j])!=1)
            break;
        HBT->heap[i] = HBT->heap[j]; //将双亲元素下移到待调整元素的位置
        i = j; //使待调整位置变为其双亲位置，进行下一次循环
    }
    HBT->heap[i] = *x;//把新元素调整到最终位置
}


inline int DeleteHeap(HeapSq* HBT, ElemType* get)
{
    ElemType temp, x;
    int i, j;
    if (HBT->len == 0)
    {
        return 0;
    }
    temp = HBT->heap[0]; //暂存堆顶元素
    HBT->len--;
    if (HBT->len == 0) //若删除操作后堆为空则返回
    {
        *get = temp;
        return 2;
    }

    x = HBT->heap[HBT->len]; //将待调整的原堆尾元素暂存x中，以便放入最终位置
    i = 0; //用i指向待调整元素的位置，初始指向堆顶位置
    j = 2 * i + 1;//用j指向i的左孩子位置，初始指向下标为1的位置
    while (j <= HBT->len - 1)//寻找待调整元素的最终位置，每次使孩子元素上移一层，调整到孩子为空时止
    {
        ///if (j < HBT->len - 1 && HBT->heap[j] > HBT->heap[j+1])//若存在右孩子且较小，使j指向右孩子
        if (j < HBT->len - 1 && cmp_ElemType(&HBT->heap[j], &HBT->heap[j + 1]) == 2)
            j++;
        ///if (x <= HBT->heap[j]) //若x比其较小的孩子还小，则调整结束，退出循环
        if (cmp_ElemType(&x, &HBT->heap[j])!=2)
            break;
        HBT->heap[i] = HBT->heap[j];//否则，将孩子元素移到双亲位置
        i = j; //将待调整位置变为其较小的孩子位置
        j = 2 * i + 1;//将j变为新的待调整位置的左孩子位置，继续下一次循环
    }
    HBT->heap[i] = x; //把x放到最终位置

    //返回原堆顶元素
    *get = temp;
    return 1;
}





void init_overlap_region_alloc(overlap_region_alloc* list)
{
    list->size = 1000;
    list->length = 0;
    list->list = (overlap_region*)malloc(sizeof(overlap_region)*list->size);
}
void clear_overlap_region_alloc(overlap_region_alloc* list)
{
    list->length = 0;
}

void destory_overlap_region_alloc(overlap_region_alloc* list)
{
    free(list->list);
}




///r->length = Get_READ_LENGTH((*R_INF), ID);
void append_overlap_region_alloc(overlap_region_alloc* list, overlap_region* tmp, All_reads* R_INF)
{
   
    if (list->length + 1 > list->size)
    {
        list->size = list->size * 2;
        list->list = (overlap_region*)realloc(list->list, sizeof(overlap_region)*list->size);
    }

    if (list->length!=0 && 
    list->list[list->length - 1].y_id==tmp->y_id
    )
    {    
        if(list->list[list->length - 1].shared_seed >= tmp->shared_seed)
        {
            return;
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
    
    tmp->x_pos_e = Get_READ_LENGTH((*R_INF), tmp->x_id) - tmp->x_pos_s - 1;
    tmp->y_pos_e = Get_READ_LENGTH((*R_INF), tmp->y_id) - tmp->y_pos_s - 1;

    if(tmp->x_pos_e <= tmp->y_pos_e)
    {
        tmp->y_pos_e = tmp->y_pos_s + tmp->x_pos_e;
        tmp->x_pos_e = tmp->x_pos_s + tmp->x_pos_e;
    }
    else
    {
        tmp->x_pos_e = tmp->x_pos_s + tmp->y_pos_e;
        tmp->y_pos_e = tmp->y_pos_s + tmp->y_pos_e;
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


void debug_overlap_region(Candidates_list* candidates, overlap_region_alloc* overlap_list, uint64_t readID)
{
    uint64_t i = 0;
    uint64_t total_length = 0;
    uint64_t pre_total_length = 0;
    uint64_t current_ID;
    uint64_t current_stand;
    long long current_pos_diff;
    long long current_self_pos;
    long long tmp_pos_distance;
    long long tmp_self_pos_distance;
    long long constant_distance = 5;
    double error_rate = 0.05;



    for ( i = 0; i < overlap_list->length; i++)
    {
        pre_total_length = total_length;
        total_length = total_length + overlap_list->list[i].shared_seed;

        uint64_t j = 0;
        current_ID = candidates->list[pre_total_length].readID;
        current_stand = candidates->list[pre_total_length].strand;
        current_pos_diff = candidates->list[pre_total_length].offset - candidates->list[pre_total_length].self_offset;
        current_self_pos = candidates->list[pre_total_length].self_offset;
        for (j = pre_total_length; j < total_length; j++)
        {
            if(current_ID != candidates->list[j].readID || 
            current_stand != candidates->list[j].strand)
            {
                fprintf(stderr, "ERROR Overlap 1!\n");
            }


            ///这个一定是正值
            tmp_pos_distance = candidates->list[j].offset - candidates->list[j].self_offset - current_pos_diff;
            ///这个不一定是正值
            tmp_self_pos_distance = candidates->list[j].self_offset - current_self_pos;
            if (tmp_self_pos_distance < 0)
            {
                fprintf(stderr, "ERROR Overlap 3!\n");
            }

            if (tmp_self_pos_distance >= 0)
            {
                tmp_self_pos_distance = tmp_self_pos_distance * error_rate + constant_distance;
                if (tmp_pos_distance >= tmp_self_pos_distance)
                {
                    fprintf(stderr, "ERROR Overlap 4!\n");
                }
            }
            
            
        }

        if (total_length < candidates->length)
        {
            j = total_length;

            if (current_ID == candidates->list[j].readID && 
            current_stand == candidates->list[j].strand)
            {
                ///这个一定是正值
                tmp_pos_distance = candidates->list[j].offset - candidates->list[j].self_offset - current_pos_diff;
                ///这个不一定是正值
                tmp_self_pos_distance = candidates->list[j].self_offset - current_self_pos;
                if (tmp_self_pos_distance >= 0)
                {
                    tmp_self_pos_distance = tmp_self_pos_distance * error_rate + constant_distance;
                    if (tmp_pos_distance < tmp_self_pos_distance)
                    {
                       fprintf(stderr, "ERROR Overlap 5!\n");
                    }
                }
            }
        }
        
        
    }

    if (total_length != candidates->length)
    {
        fprintf(stderr, "ERROR Overlap 2!\n");
        fprintf(stderr, "total_length: %llu\n", total_length);
        fprintf(stderr, "candidates->length: %llu\n", candidates->length);
    }
}



void print_overlap_region(Candidates_list* candidates, overlap_region_alloc* overlap_list, All_reads* R_INF)
{

    pthread_mutex_lock(&output_mutex);


    uint64_t i = 0;
    for ( i = 0; i < overlap_list->length; i++)
    {
        fprintf(stderr, "************i: %llu***********\n", i);
        fprintf(stderr, "x: %llu, %llu, %llu, %llu, %llu\n",
        overlap_list->list[i].x_id,overlap_list->list[i].x_pos_s, 
        overlap_list->list[i].x_pos_e, overlap_list->list[i].x_pos_strand,
        Get_READ_LENGTH((*R_INF), overlap_list->list[i].x_id)
        );

        fprintf(stderr, "y: %llu, %llu, %llu, %llu, %llu\n",
        overlap_list->list[i].y_id,overlap_list->list[i].y_pos_s, 
        overlap_list->list[i].y_pos_e, overlap_list->list[i].y_pos_strand,
        Get_READ_LENGTH((*R_INF), overlap_list->list[i].y_id)
        );
        
    }

		
		

    fprintf(stderr, "#######################\nLength: %llu\n", candidates->length);

    
    for (i = 0; i < candidates->length; i++)
    {
        fprintf(stderr, "\ni: %llu\n", i);
        fprintf(stderr, "strand: %llu\n", candidates->list[i].strand);
        fprintf(stderr, "readID: %llu\n", candidates->list[i].readID);
        fprintf(stderr, "diff: %lld\n", candidates->list[i].offset - candidates->list[i].self_offset);
        fprintf(stderr, "offset: %lld\n", candidates->list[i].offset);
        fprintf(stderr, "self_offset: %lld\n", candidates->list[i].self_offset);

    }

    pthread_mutex_unlock(&output_mutex);
}



	///r->length = Get_READ_LENGTH((*R_INF), ID);

void calculate_overlap_region(Candidates_list* candidates, overlap_region_alloc* overlap_list, 
uint64_t readID, uint64_t readLength, All_reads* R_INF)
{
    overlap_region tmp_region;
    uint64_t i = 0;
    long long current_pos_diff;
    long long current_self_pos;
    uint64_t current_ID;
    uint64_t current_stand;
    long long tmp_pos_distance;
    long long tmp_self_pos_distance;
    long long constant_distance = 5;
    double error_rate = 0.05;

    if (candidates->length == 0)
    {
        return;
    }
    

    i = 0;
    while (i < candidates->length)
    {
        current_pos_diff = candidates->list[i].offset - candidates->list[i].self_offset;
        current_self_pos = candidates->list[i].self_offset;
        current_ID = candidates->list[i].readID;
        current_stand = candidates->list[i].strand;

        tmp_region.shared_seed = 1;
        ///这个是查询read的信息
        tmp_region.x_id = readID;
        tmp_region.x_pos_s = candidates->list[i].self_offset;
        tmp_region.x_pos_e = candidates->list[i].self_offset;
        tmp_region.x_pos_strand = current_stand;

        ///这个是被查询的read的信息
        tmp_region.y_id = current_ID;
        tmp_region.y_pos_s = candidates->list[i].offset;
        tmp_region.y_pos_e = candidates->list[i].offset;
        tmp_region.y_pos_strand = 0;  ///永远是0

        i++;

        while (i < candidates->length)
        {
            if (current_ID == candidates->list[i].readID && 
            current_stand == candidates->list[i].strand)
            {
                ///这个一定是正值
                tmp_pos_distance = candidates->list[i].offset - candidates->list[i].self_offset - current_pos_diff;
                ///这个不一定是正值
                tmp_self_pos_distance = candidates->list[i].self_offset - current_self_pos;
                if (tmp_self_pos_distance >= 0)
                {
                    tmp_self_pos_distance = tmp_self_pos_distance * error_rate + constant_distance;
                    if (tmp_pos_distance < tmp_self_pos_distance)
                    {
                        i++;
                        tmp_region.shared_seed++;
                        tmp_region.x_pos_e = candidates->list[i].self_offset;
                        tmp_region.y_pos_e = candidates->list[i].offset;
                        continue;
                    }
                }
            }
            ///i++;
            break;
        }

        append_overlap_region_alloc(overlap_list, &tmp_region, R_INF);
    }


    ///debug_overlap_region(candidates, overlap_list, readID);
    ///print_overlap_region(candidates, overlap_list, R_INF);
}
























void init_k_mer_pos_list_alloc(k_mer_pos_list_alloc* list)
{
    list->size = 1000;
    list->length = 0;
    list->list = (k_mer_pos_list*)malloc(sizeof(k_mer_pos_list)*list->size);
}

void destory_k_mer_pos_list_alloc(k_mer_pos_list_alloc* list)
{
    free(list->list);
}

void clear_k_mer_pos_list_alloc(k_mer_pos_list_alloc* list)
{
    list->length = 0;
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

    int j = 0;

    for (i = 0; i < n_lengh; i++)
    {
        
        for (; j < candidates->length; j++)
        {
            if (
                n_list[i].offset == candidates->list[j].offset
                &&
                n_list[i].readID == candidates->list[j].readID
                &&
                end_pos == candidates->list[j].self_offset
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

void verify_merge_result(k_mer_pos_list_alloc* list, Candidates_list* candidates)
{
    uint64_t total_length = 0;
    uint64_t i = 0;
    long long r_dis_x, r_dis_y;
    for (i = 0; i < list->length; i++)
    {
        total_length = total_length + list->list[i].length;
    }

    if (total_length!=candidates->length)
    {
        fprintf(stderr, "ERROR length & size.\n");
    }

    for (i = 1; i < candidates->length; i++)
    {
        if (candidates->list[i].strand < candidates->list[i-1].strand)
        {
            fprintf(stderr, "ERROR -1\n");
        }
        else if (candidates->list[i].strand == candidates->list[i-1].strand)
        {
            if (candidates->list[i].readID < candidates->list[i-1].readID)
            {
                fprintf(stderr, "ERROR 0\n");
            }
            else if (candidates->list[i].readID == candidates->list[i-1].readID)
            {
                r_dis_x = candidates->list[i].offset - candidates->list[i].self_offset;
                r_dis_y = candidates->list[i-1].offset - candidates->list[i-1].self_offset;

                ///if (candidates->list[i].offset < candidates->list[i-1].offset)
                if (r_dis_x < r_dis_y)
                {
                    fprintf(stderr, "ERROR 1\n");
                }
                ///else if (candidates->list[i].offset == candidates->list[i-1].offset)
                else if (r_dis_x == r_dis_y)
                {
                    if (candidates->list[i].self_offset < candidates->list[i-1].self_offset)
                    {
                        fprintf(stderr, "ERROR 2\n");
                    }
                }
            } 
        }
    }

    
    uint64_t j = 0;
    for (i = 0; i < list->length; i++)
    {
        test_single_list(candidates, list->list[i].list, list->list[i].length, list->list[i].end_pos, list->list[i].direction);
    }
    
    
    
}


void output_to_stderr(Candidates_list* candidates)
{
    uint64_t i = 0;

    pthread_mutex_lock(&output_mutex);
		
		

    fprintf(stderr, "************************\nLength: %llu\n", candidates->length);

    
    for (i = 0; i < candidates->length; i++)
    {
        fprintf(stderr, "\ni: %llu\n", i);
        fprintf(stderr, "strand: %llu\n", candidates->list[i].strand);
        fprintf(stderr, "readID: %llu\n", candidates->list[i].readID);
        fprintf(stderr, "diff: %lld\n", candidates->list[i].offset - candidates->list[i].self_offset);
        fprintf(stderr, "offset: %lld\n", candidates->list[i].offset);
        fprintf(stderr, "self_offset: %lld\n", candidates->list[i].self_offset);

    }

    pthread_mutex_unlock(&output_mutex);

    
}


void merge_k_mer_pos_list_alloc_heap_sort_back(k_mer_pos_list_alloc* list, Candidates_list* candidates, HeapSq* HBT)
{
    clear_Heap(HBT);

    uint64_t total_length = 0;
    uint64_t i;
    ElemType x, y;
    ///所有list的长度都不是0
    ///把各个表第一个元素加入到堆中
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
    if(total_length > candidates->size)
    {
        candidates->size = total_length;
        candidates->list = (k_mer_hit*)realloc(candidates->list, sizeof(k_mer_hit)*candidates->size);
        candidates->tmp = (k_mer_hit*)realloc(candidates->tmp, sizeof(k_mer_hit)*candidates->size);
    }

    uint64_t ID;
    while (DeleteHeap(HBT, &x))
    {
        append_pos_to_Candidates_list(candidates, &x);
        ///x.ID说明是从第x.ID个列表中的这个节点已经从堆里出来了
        ///HBT->index_i[x.ID]是第x.ID个列表的当前元素的下标
        i = HBT->index_i[x.ID];
        ID = x.ID;
        ///fprintf(stderr, "x.ID: %llu, i: %llu, length: %llu\n", x.ID, i, list->list[x.ID].length);


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

    ///verify_merge_result(list, candidates);
    
    
     
}


void merge_k_mer_pos_list_alloc_heap_sort(k_mer_pos_list_alloc* list, Candidates_list* candidates, HeapSq* HBT)
{
    clear_Heap(HBT);

    uint64_t total_length = 0;
    uint64_t i;
    ElemType x, y;
    ///所有list的长度都不是0
    ///把各个表第一个元素加入到堆中
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
    if(total_length > candidates->size)
    {
        candidates->size = total_length;
        candidates->list = (k_mer_hit*)realloc(candidates->list, sizeof(k_mer_hit)*candidates->size);
        candidates->tmp = (k_mer_hit*)realloc(candidates->tmp, sizeof(k_mer_hit)*candidates->size);
    }

    uint64_t ID;
    int flag;
    while (flag = DeleteHeap(HBT, &x))
    {
        append_pos_to_Candidates_list(candidates, &x);
        ///x.ID说明是从第x.ID个列表中的这个节点已经从堆里出来了
        ///HBT->index_i[x.ID]是第x.ID个列表的当前元素的下标
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


    ///verify_merge_result(list, candidates);
    ///output_to_stderr(candidates);
     
}


///有bug，找时间排一下
void merge_k_mer_pos_list_alloc_heap_sort_advance(k_mer_pos_list_alloc* list, Candidates_list* candidates, HeapSq* HBT)
{

    uint64_t total_length = 0;
    uint64_t i;
    ElemType x, y;

    /************************************forward*************************************************/
    clear_Heap(HBT);

    ///所有list的长度都不是0
    ///把各个表第一个元素加入到堆中
    for (i = 0; i < list->length; i++)
    {
        if (list->list[i].direction == 0)
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
    }


    candidates->length = 0;
    if(total_length > candidates->size)
    {
        candidates->size = total_length;
        candidates->list = (k_mer_hit*)realloc(candidates->list, sizeof(k_mer_hit)*candidates->size);
        candidates->tmp = (k_mer_hit*)realloc(candidates->tmp, sizeof(k_mer_hit)*candidates->size);
    }


    uint64_t ID;
    int flag;
    while (flag = DeleteHeap(HBT, &x))
    {
        append_pos_to_Candidates_list(candidates, &x);
        ///x.ID说明是从第x.ID个列表中的这个节点已经从堆里出来了
        ///HBT->index_i[x.ID]是第x.ID个列表的当前元素的下标
        i = HBT->index_i[x.ID];
        ID = x.ID;
        ///fprintf(stderr, "x.ID: %llu, i: %llu, length: %llu\n", x.ID, i, list->list[x.ID].length);
        /**
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
        **/
        
        
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



    /************************************reverse complement*************************************************/

    clear_Heap(HBT);
    ///所有list的长度都不是0
    ///把各个表第一个元素加入到堆中
    for (i = 0; i < list->length; i++)
    {
        if (list->list[i].direction == 1)
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
    }

    if(total_length > candidates->size)
    {
        candidates->size = total_length;
        candidates->list = (k_mer_hit*)realloc(candidates->list, sizeof(k_mer_hit)*candidates->size);
        candidates->tmp = (k_mer_hit*)realloc(candidates->tmp, sizeof(k_mer_hit)*candidates->size);
    }

    while (flag = DeleteHeap(HBT, &x))
    {
        append_pos_to_Candidates_list(candidates, &x);
        ///x.ID说明是从第x.ID个列表中的这个节点已经从堆里出来了
        ///HBT->index_i[x.ID]是第x.ID个列表的当前元素的下标
        i = HBT->index_i[x.ID];
        ID = x.ID;
        ///fprintf(stderr, "x.ID: %llu, i: %llu, length: %llu\n", x.ID, i, list->list[x.ID].length);
        /**
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
        **/
        

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


void merge_k_mer_pos_list_alloc(k_mer_pos_list_alloc* list, Candidates_list* candidates)
{
    qsort(list->list, list->length, sizeof(k_mer_pos_list), cmp_k_mer_pos_list);

    uint64_t i;
    for (i = 0; i < list->length; i++)
    {
        /**
        if (i > 0 && (list->list[i].length < list->list[i-1].length || list->list[i].length == 0 || list->list[i-1].length == 0))
        {
            fprintf(stderr, "ERROR\n");
        }
        **/
        
        merge_Candidates_list(candidates, list->list[i].list, list->list[i].length, 
        list->list[i].end_pos, list->list[i].direction);
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
        fprintf(stdout, "k-mer is too long. The length of k-mer must <= 64.");
        fflush(stdout);
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
    ///这个右移是安全的，因为TCB->suffix_bits不可能是0
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
    fprintf(stdout, "Writing index to disk ...... \n");
    char* index_name = (char*)malloc(strlen(read_file_name)+5);
    sprintf(index_name, "%s.idx", read_file_name);
    FILE* fp = fopen(index_name, "w");
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
    fprintf(stdout, "Index has been written.\n");
}


int load_Total_Pos_Table(Total_Pos_Table* TCB, char* read_file_name)
{
    fprintf(stdout, "Loading index to disk ...... \n");
    char* index_name = (char*)malloc(strlen(read_file_name)+5);
    sprintf(index_name, "%s.idx", read_file_name);
    FILE* fp = fopen(index_name, "r");
    if (!fp)
    {
        return 0;
    }
    
    
    fread(&TCB->prefix_bits, sizeof(TCB->prefix_bits), 1, fp);
    fread(&TCB->suffix_bits, sizeof(TCB->suffix_bits), 1, fp);
    fread(&TCB->suffix_mode, sizeof(TCB->suffix_mode), 1, fp);
    fread(&TCB->size, sizeof(TCB->size), 1, fp);
    fread(&TCB->useful_k_mer, sizeof(TCB->useful_k_mer), 1, fp);
    fread(&TCB->total_occ, sizeof(TCB->total_occ), 1, fp);


    if (TCB->useful_k_mer+1)
    {
        TCB->k_mer_index = (uint64_t*)malloc(sizeof(uint64_t)*(TCB->useful_k_mer+1));
        fread(TCB->k_mer_index, sizeof(uint64_t), TCB->useful_k_mer+1, fp);
    }
    else
    {
        TCB->k_mer_index = NULL;
    }
    
    
    if (TCB->total_occ)
    {
        TCB->pos = (k_mer_pos*)malloc(sizeof(k_mer_pos)*TCB->total_occ);
        fread(TCB->pos, sizeof(k_mer_pos), TCB->total_occ, fp);
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
    fprintf(stdout, "Index has been loaded.\n");

    return 1;
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

    ///init_Total_Pos_Table(PCB, TCB);

    khint_t t;  ///这就是个迭代器
    int absent;

    for (i = 0; i < TCB->size; i++)
    {
        h = TCB->sub_h[i];
        for (k = kh_begin(h); k != kh_end(h); ++k)
        {
            if (kh_exist(h, k))            // test if a bucket contains data
            {
                ///只有符合频率范围要求的k-mer，才会被加入到pos table中
                if (kh_value(h, k)>=k_mer_min_freq && kh_value(h, k)<=k_mer_max_freq)
                {

                    sub_ID = i;
                    sub_key = kh_key(h, k);

                    t = kh_put(POS64, PCB->sub_h[sub_ID], sub_key, &absent);
    
                    if (absent)
                    {
                        ///kh_value(PCB->sub_h[sub_ID], t) = useful_k_mer + total_occ;
                        kh_value(PCB->sub_h[sub_ID], t) = PCB->useful_k_mer;
                    }
                    else   ///哈希表中已有的元素
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


    fprintf(stdout, "useful_k_mer: %lld\n",PCB->useful_k_mer);
    fprintf(stdout, "total_occ: %lld\n",PCB->total_occ);

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
                if (kh_value(h, k)>=k_mer_min_freq && kh_value(h, k)<=k_mer_max_freq)
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


void init_Candidates_list(Candidates_list* l)
{
    l->length = 0;
    l->size = 0;
    l->list = NULL;
    l->tmp = NULL;
    l->foward_pos = 0;
    l->rc_pos = 0;
}

void clear_Candidates_list(Candidates_list* l)
{
    l->length = 0;
    l->foward_pos = 0;
    l->rc_pos = 0;
}

void destory_Candidates_list(Candidates_list* l)
{
    free(l->list);
    free(l->tmp);
}

///1是x小，2是y小，0是相等
inline int cmp_pos(k_mer_hit* x, k_mer_pos* y, uint64_t y_pos, uint8_t strand)
{
    if (x->strand < strand)
    {
        return 1;
    }
    else if (x->strand > strand)
    {
        return 2;
    }
    else
    {
        if (x->readID < y->readID)
        {
            return 1;
        }
        else if (x->readID > y->readID)
        {
            return 2;
        }
        else
        {
            if (x->offset < y->offset)
            {
                return 1;
            }
            else if (x->offset > y->offset)
            {
                return 2;
            }
            else
            {
                if (x->self_offset < y_pos)
                {
                    return 1;
                }
                else  if (x->self_offset > y_pos)
                {
                    return 2;
                }
                else
                {
                    return 0;
                }
                
            }   
        }    
    }
    
}

void debug_merge_Candidates_list(Candidates_list* l, k_mer_pos* n_list, uint64_t n_lengh, uint64_t end_pos, uint64_t i,
uint64_t i_1, uint64_t i_2, uint64_t len1, uint64_t len2, uint64_t strand)
{
    if (i != l->length)
    {
        fprintf(stderr, "ERROR\n");
        fprintf(stderr, "i: %llu, l->length: %llu\n", i , l->length);
        fprintf(stderr, "i_1: %llu, len1: %llu\n", i_1 , len1);
        fprintf(stderr, "i_2: %llu, len2: %llu\n", i_2 , len2);
        fprintf(stderr, "strand: %llu\n", strand);
        
    }

    for (i = 1; i < l->length; i++)
    {
        if (l->list[i].strand < l->list[i-1].strand)
        {
            fprintf(stderr, "ERROR -1\n");
        }
        else if (l->list[i].strand == l->list[i-1].strand)
        {
            if (l->list[i].readID < l->list[i-1].readID)
            {
                fprintf(stderr, "ERROR 0\n");
            }
            else if (l->list[i].readID == l->list[i-1].readID)
            {
                if (l->list[i].offset < l->list[i-1].offset)
                {
                    fprintf(stderr, "ERROR 1\n");
                }
                else if (l->list[i].offset == l->list[i-1].offset)
                {
                    if (l->list[i].self_offset < l->list[i-1].self_offset)
                    {
                        fprintf(stderr, "ERROR 2\n");
                    }
                }
            } 
        }
    }

    for (i = 0; i < n_lengh; i++)
    {
        int j = 0;
        for (j = 0; j < l->length; j++)
        {
            if (
                n_list[i].offset == l->list[j].offset
                &&
                n_list[i].readID == l->list[j].readID
                &&
                end_pos == l->list[j].self_offset
                &&
                strand == l->list[j].strand
            )
            {
                break;
            }   
        }

        if (j == l->length)
        {
            fprintf(stderr, "ERROR 4\n");
        }
    }


    for (i = 0; i < l->length - n_lengh; i++)
    {
        int j = 0;
        for (j = 0; j < l->length; j++)
        {
            if (
                l->tmp[i].offset == l->list[j].offset
                &&
                l->tmp[i].readID == l->list[j].readID
                &&
                l->tmp[i].self_offset == l->list[j].self_offset
                &&
                l->tmp[i].strand == l->list[j].strand
            )
            {
                break;
            }   
        }

        if (j == l->length)
        {
            fprintf(stderr, "ERROR 5\n");
        }
    }
}

void merge_Candidates_list(Candidates_list* l, k_mer_pos* n_list, uint64_t n_lengh, uint64_t end_pos, int strand)
{
    if (n_lengh == 0) ///不加这个就会巨慢
    {
        return;
    }
    

    if(l->length + n_lengh > l->size)
    {
        l->size = l->length + n_lengh;
        l->list = (k_mer_hit*)realloc(l->list, sizeof(k_mer_hit)*l->size);
        l->tmp = (k_mer_hit*)realloc(l->tmp, sizeof(k_mer_hit)*l->size);
    }

    int flag;
    uint64_t len1, len2, i_1, i_2, i;
    k_mer_hit *l_1;
    k_mer_pos *l_2;

    l_1 = l->list;
    if (strand == 0)
    {
        ///由于标号全是0，则标号为1的位置不受影响(数组后半部分)，可以直接复制过去
        ///需要从l_1 + l->rc_pos开始，复制l->length-l->rc_pos个元素，到l->tmp + l->rc_pos + n_lengh处
        memcpy(l->tmp + l->rc_pos + n_lengh, l_1 + l->rc_pos, sizeof(k_mer_hit)*(l->length-l->rc_pos));
        ///l_1只需要从0扫描到l->rc_pos就好了
        len1 = l->rc_pos; 
        i_1 = 0;
        ///此时目标数组没数据，
        i = 0;
        ///更新l->rc_pos
        l->rc_pos = l->rc_pos + n_lengh;
    }
    else
    {
        ///如果strand==1，则标号为0的位置不受影响(数组前半部分)，可以直接复制过去
        ///从l_1开始，复制l->rc_pos个元素，到l->tmp
        memcpy(l->tmp, l_1, sizeof(k_mer_hit)*l->rc_pos);
        ///l_1只需要从l->rc_pos扫描到末尾就好了
        i_1 = l->rc_pos;
        len1 = l->length; 
        ///此时目标数组已经有了l->rc_pos
        i = i_1;
        ///l->rc_pos不用更新
    }
    
    
    

    l_2 = n_list;
    len2 = n_lengh;
    i_2 = 0;


    while (i_1<len1 && i_2<len2)
    {
        flag = cmp_pos(&l_1[i_1], &l_2[i_2], end_pos, strand); 
        ///l_1小
        if (flag == 1)
        {
            l->tmp[i].readID = l_1[i_1].readID;
            l->tmp[i].offset = l_1[i_1].offset;
            l->tmp[i].self_offset = l_1[i_1].self_offset;
            l->tmp[i].strand = l_1[i_1].strand;
            i_1++;
            i++;
        }
        else  ///l_2小或者l_1 == l_2
        {
            l->tmp[i].readID = l_2[i_2].readID;
            l->tmp[i].offset = l_2[i_2].offset;
            l->tmp[i].self_offset = end_pos;
            l->tmp[i].strand = strand;
            i_2++;
            i++;
        }        
    }


    while (i_1<len1)
    {
        memcpy(l->tmp + i, l_1 + i_1, sizeof(k_mer_hit)*(len1-i_1));
        i = i + len1-i_1;
        i_1 = len1;
    }

    while (i_2<len2)
    {
        l->tmp[i].readID = l_2[i_2].readID;
        l->tmp[i].offset = l_2[i_2].offset;
        l->tmp[i].self_offset = end_pos;
        l->tmp[i].strand = strand;
        i_2++;
        i++;
    }
    
    

    k_mer_hit* k;
    k = l->list;
    l->list = l->tmp;
    l->tmp = k;
    l->length = l->length + n_lengh;

    
    /**
    if (strand == 0)
    {
        i = l->length;
        len1 = l->length - n_lengh;
    }
    debug_merge_Candidates_list(l, n_list, n_lengh, end_pos, i, i_1, i_2, len1, len2, strand);
    **/
    
}







/********************************for debug***************************************/
void merge_Candidates_list_version(Candidates_list* l, k_mer_pos* n_list, uint64_t n_lengh, uint64_t end_pos, int strand)
{
    if (n_lengh == 0) ///不加这个就会巨慢
    {
        return;
    }
    

    if(l->length + n_lengh > l->size)
    {
        l->size = l->length + n_lengh;
        l->list = (k_mer_hit*)realloc(l->list, sizeof(k_mer_hit)*l->size);
        l->tmp = (k_mer_hit*)realloc(l->tmp, sizeof(k_mer_hit)*l->size);
    }

    int flag;
    uint64_t len1, len2, i_1, i_2, i;
    k_mer_hit *l_1;
    k_mer_pos *l_2;

    l_1 = l->list;
    len1 = l->length; 
    i_1 = 0;

    l_2 = n_list;
    len2 = n_lengh;
    i_2 = 0;

    i = 0;

    while (i_1<len1 && i_2<len2)
    {
        flag = cmp_pos(&l_1[i_1], &l_2[i_2], end_pos, strand); 
        ///l_1小
        if (flag == 1)
        {
            l->tmp[i].readID = l_1[i_1].readID;
            l->tmp[i].offset = l_1[i_1].offset;
            l->tmp[i].self_offset = l_1[i_1].self_offset;
            l->tmp[i].strand = l_1[i_1].strand;
            i_1++;
            i++;
        }
        else  ///l_2小或者l_1 == l_2
        {
            l->tmp[i].readID = l_2[i_2].readID;
            l->tmp[i].offset = l_2[i_2].offset;
            l->tmp[i].self_offset = end_pos;
            l->tmp[i].strand = strand;
            i_2++;
            i++;
        }        
    }


    while (i_1<len1)
    {
        l->tmp[i].readID = l_1[i_1].readID;
        l->tmp[i].offset = l_1[i_1].offset;
        l->tmp[i].self_offset = l_1[i_1].self_offset;
        l->tmp[i].strand = l_1[i_1].strand;
        i_1++;
        i++;
    }

    while (i_2<len2)
    {
        l->tmp[i].readID = l_2[i_2].readID;
        l->tmp[i].offset = l_2[i_2].offset;
        l->tmp[i].self_offset = end_pos;
        l->tmp[i].strand = strand;
        i_2++;
        i++;
    }
    
    k_mer_hit* k;
    k = l->list;
    l->list = l->tmp;
    l->tmp = k;
    l->length = l->length + n_lengh;

    ///debug_merge_Candidates_list(l, n_list, n_lengh, end_pos, i, i_1, i_2, len1, len2, strand);
}







/********************************for debug***************************************/
void debug_mode(uint64_t d, uint64_t thread_ID, uint64_t thread_num)
{
    fprintf(stderr, "#%llu: Start debug_mode...\n", thread_ID);

    uint64_t i;

    for (i = thread_ID; i < 4,294,967,296; i = i + thread_num)
    {
        if(i%d != mod_d(0, i, d))
        {
            fprintf(stderr, "ERROR mode, i: %llu\n", i);
            fprintf(stderr, "i mod d: %llu\n", i%d);
            fprintf(stderr, "mod_d(0, i, d): %llu\n", mod_d(0, i, d));
        }

    }

    fprintf(stderr, "#%llu: Finish debug_mode\n", thread_ID);
    
}


/********************************for debug***************************************/
void test_COUNT64()
{

    uint64_t sb[7] = {5000000000, 6000000000, 7000000000, 8000000000, 9000000000, 10000000000 ,5000000000};
    Count_Table* h;
    ///构建哈希表
    init_Count_Table(&h);

    int i;
    khint_t k;  ///这就是个迭代器
    int absent;

    for (i = 0; i < 7; i++)
    {
        ///将5作为key插入到哈希表COUNT64中
        ///absent为0代表哈希表里面已经有这个key了
        k = kh_put(COUNT64, h, sb[i], &absent);
        ///代表插入了一个哈希表中没有的元素
        ///则直接插入
        if (absent)
        {
            kh_value(h, k) = i;
        }
        else   ///哈希表中已有的元素
        {
            kh_value(h, k) = i;
        }
    }

    for (i = 0; i < 7; i++)
    {
        ///查询哈希表，key为k
        k = kh_get(COUNT64, h, sb[i]);

        if (k != kh_end(h))
        {
            fprintf(stderr, "value: %d\n", kh_value(h, k));
        }
        else
        {
            fprintf(stderr, "not found!\n");
        }
    }
    

    kh_destroy(COUNT64, h);
    

}