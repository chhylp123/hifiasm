#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "Hash_Table.h"


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