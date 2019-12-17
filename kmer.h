#ifndef __KMER__
#define __KMER__
#include "Process_Read.h"

///#define ALL (0xffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff)
#define ALL (0xffffffffffffffff)
/****************************may have bugs********************************/
#define SAFE_SHIFT(k) k & ((k < 64)?ALL:0)
/****************************may have bugs********************************/






typedef struct
{
    ///最大64-mer
    ///x[0]低位
    ///x[1]高位
    uint64_t x[2];

} Hash_code;

typedef struct {
    Hash_code key; ///k-mer itself
    uint64_t value;  ///offset
} k_v;

typedef struct {
	k_v* buffer;
    uint32_t size;
    uint32_t length;
} small_hash_table;

void init_small_hash_table(small_hash_table* x);
void clear_small_hash_table(small_hash_table* x);
void resize_small_hash_table(small_hash_table* x, uint64_t size);
void destory_small_hash_table(small_hash_table* x);
void add_small_hash_table(small_hash_table* x, k_v* element);
void sort_small_hash_table(small_hash_table* x);
int compare_k_mer(k_v* x, k_v* y);
int query_small_hash_table(small_hash_table* target, k_v* query, long long* l_end, long long* r_end);


typedef struct
{
	char* str;
	long long l;
    long long i;
    long long N_occ;

} HPC_seq;


inline uint64_t get_HPC_code(HPC_seq* seq, uint64_t* end_pos)
{

    if(seq->i < seq ->l)
    {
        uint8_t code = seq_nt6_table[(uint8_t)seq->str[seq->i]];

        (*end_pos) = seq->i;

        for (; seq->i < seq->l; seq->i++)
        {
            ///统计N的个数
            if (seq_nt6_table[(uint8_t)seq->str[seq->i]] >= 4)
            {
                seq->N_occ++;
            }
            
            if (seq_nt6_table[(uint8_t)seq->str[seq->i]] != code)
            {
                break;
            }
        }

        return (uint64_t)code;
    }
    else
    {
        ///end
        return 6;
    }
    
}

inline void k_mer_append(Hash_code* code, uint64_t c, int k)
{

    uint64_t mask = ALL >> (64 -k);

    code->x[0] = ((code->x[0]<<1) | (c&1))  & mask;
    code->x[1] = ((code->x[1]<<1) | (c>>1)) & mask;
}

inline void Hashcode_to_string(Hash_code* code, char* str, int k)
{
    uint8_t c;
    int i;
    for (i = 0; i < k; i++)
    {
        c = (code->x[1] >> (k - i - 1)) & ((uint64_t)1);
        c = c << 1;
        c = c | ((code->x[0] >> (k - i - 1)) & ((uint64_t)1));

        str[i] = s_H[c];
    }
    
}

void init_HPC_seq(HPC_seq* seq, char* str, long long l);
void init_Hash_code(Hash_code* code);


#endif