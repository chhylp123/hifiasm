#ifndef __HASHTABLE__
#define __HASHTABLE__
#include "khash.h"
#include "kmer.h"

KHASH_MAP_INIT_INT64(COUNT64, int) 

typedef khash_t(COUNT64) Count_Table;

#define PREFIX_BITS 16
#define MAX_SUFFIX_BITS 64
#define MODE_VALUE 101

typedef struct
{
    volatile int lock;

}Hash_table_spin_lock;

typedef struct
{
	Count_Table** sub_h;
    Hash_table_spin_lock* sub_h_lock;
    int prefix_bits;
    int suffix_bits;
	int size;
    uint64_t suffix_mode;
} Total_Count_Table;

/********************************for debug***************************************/
inline void print_64bit(uint64_t x)
{
    int i;
    for(i = 63; i >= 0; i--)
	{
		if(x & ((1ULL<<i)))
			fprintf(stderr, "1");
		else
			fprintf(stderr, "0");
	}

    fprintf(stderr, "\n");
}

inline uint64_t mod_d(uint64_t h_key, uint64_t low_key, uint64_t d)
{
    uint64_t result = (h_key >> 32) % d;
    result = ((result << 32) + (h_key & (uint64_t)0xffffffff)) % d;
    result = ((result << 32) + (low_key >> 32)) % d;
    result = ((result << 32) + (low_key & (uint64_t)0xffffffff)) % d;

    return result;
}



inline int get_sub_table(uint64_t* get_sub_ID, uint64_t* get_sub_key, Total_Count_Table* TCB, Hash_code* code, int k)
{
    uint64_t h_key, low_key;
    ///k有可能是64，所以可能会有问题
    ///low_key = code->x[0] | (code->x[1] << k);
    low_key = code->x[0] | (code->x[1] << SAFE_SHIFT(k));
    //k不可能为0, 所以这个右移不会有问题
    h_key = code->x[1] >> (64 - k);

    if(mod_d(h_key, low_key, MODE_VALUE) > 3)
    {
        return 0;
    }


    ///注意suffix_bits最大就是64
    ///前一个右移不安全，因为TCB->suffix_bits有可能为64
    ///后一个左移安全，因为TCB->suffix_bits不可能为0
    //uint64_t sub_ID = (low_key >> TCB->suffix_bits) | (h_key << (64 - TCB->suffix_bits));
    uint64_t sub_ID = (low_key >> SAFE_SHIFT(TCB->suffix_bits)) | (h_key << (64 - TCB->suffix_bits));
    uint64_t sub_key = (low_key & TCB->suffix_mode);

    *get_sub_ID = sub_ID;
    *get_sub_key = sub_key;

     return 1;

}


inline int insert_Total_Count_Table(Total_Count_Table* TCB, Hash_code* code, int k)
{
    uint64_t sub_ID, sub_key;
    if(!get_sub_table(&sub_ID, &sub_key, TCB, code, k))
    {
        return 0;
    }

    khint_t t;  ///这就是个迭代器
    int absent;


    while (__sync_lock_test_and_set(&TCB->sub_h_lock[sub_ID].lock, 1))
    {
        while (TCB->sub_h_lock[sub_ID].lock);
    }

    t = kh_put(COUNT64, TCB->sub_h[sub_ID], sub_key, &absent);
    if (absent)
    {
        kh_value(TCB->sub_h[sub_ID], t) = 1;
    }
    else   ///哈希表中已有的元素
    {
        //kh_value(TCB->sub_h[sub_ID], t) = kh_value(TCB->sub_h[sub_ID], t) + 1;
        kh_value(TCB->sub_h[sub_ID], t)++;
    }

    __sync_lock_release(&TCB->sub_h_lock[sub_ID].lock);

    return 1;
}

inline int get_Total_Count_Table(Total_Count_Table* TCB, Hash_code* code, int k)
{

    uint64_t sub_ID, sub_key;
    if(!get_sub_table(&sub_ID, &sub_key, TCB, code, k))
    {
        return 0;
    }

    khint_t t;  ///这就是个迭代器
    int absent;

    ///查询哈希表，key为k
    t = kh_get(COUNT64, TCB->sub_h[sub_ID], sub_key);

    if (t != kh_end(TCB->sub_h[sub_ID]))
    {
        return kh_value(TCB->sub_h[sub_ID], t);
    }
    else
    {
        return -1;
    }

    return 1;

}

/********************************for debug***************************************/
inline int verify_Total_Count_Table(Total_Count_Table* TCB, Hash_code* code, int k)
{
    uint64_t sub_ID, sub_key;
    if(!get_sub_table(&sub_ID, &sub_key, TCB, code, k))
    {
        return 0;
    }

    khint_t t;  ///这就是个迭代器
    int absent;

    ///查询哈希表，key为k
    t = kh_get(COUNT64, TCB->sub_h[sub_ID], sub_key);

    if (t != kh_end(TCB->sub_h[sub_ID]))
    {
        kh_value(TCB->sub_h[sub_ID], t)--;
        if (kh_value(TCB->sub_h[sub_ID], t)<0)
        {
           return -1;
        }
        else
        {
            return 1;
        }
    }
    else
    {
        return -1;
    }
}

/********************************for debug***************************************/
inline int Traverse_Total_Count_Table(Total_Count_Table* TCB)
{
    int i;
    Count_Table* h;
    khint_t k;

    long long non_empty_k_mer = 0;

    for (i = 0; i < TCB->size; i++)
    {
        h = TCB->sub_h[i];
        for (k = kh_begin(h); k != kh_end(h); ++k)
        {
            if (kh_exist(h, k))            // test if a bucket contains data
            {
                non_empty_k_mer++;

                if (kh_value(h, k)!= 0)
                {
                    fprintf(stderr, "ERROR when Traversing!\n");
                }
            }
        }
    }
    
    fprintf(stdout, "non_empty_k_mer: %lld\n", non_empty_k_mer);
}


void init_Total_Count_Table(int k, Total_Count_Table* TCB);
void destory_Total_Count_Table(Total_Count_Table* TCB);

void init_Count_Table(Count_Table** table);

/********************************for debug***************************************/
void test_COUNT64();

/********************************for debug***************************************/
void debug_mode(uint64_t d, uint64_t thread_ID, uint64_t thread_num);


#endif