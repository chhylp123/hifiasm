#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "Hash_Table.h"




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