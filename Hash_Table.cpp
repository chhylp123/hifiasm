#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "Hash_Table.h"


void init_Count_Table(Count_Table** table)
{
    *table = kh_init(COUNT64);
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

    int i = 0;
    for (i = 0; i < TCB->size; i++)
    {
        init_Count_Table(&(TCB->sub_h[i]));
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
