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
    //can represent at most 64-mer
    uint64_t x[4];
} Hash_code;

typedef struct
{
	char* str;
	long long l;
    long long i;
    long long N_occ;

} HPC_seq;

inline void init_HPC_seq(HPC_seq* seq, char* str, long long l)
{
    seq->i = 0;
    seq->l = l;
    seq->N_occ = 0;
    seq->str = str;
}

inline uint64_t get_HPC_code(HPC_seq* seq, uint64_t* end_pos)
{

    if(seq->i < seq ->l)
    {
        uint8_t code = seq_nt6_table[(uint8_t)seq->str[seq->i]];

        (*end_pos) = seq->i;

        for (; seq->i < seq->l; seq->i++)
        {
            ///number of Ns
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

inline void init_Hash_code(Hash_code* code)
{
    code->x[0] = code->x[1] = code->x[2] = code->x[3] = 0;
}

inline void k_mer_append(Hash_code* code, uint64_t c, int k)
{
    uint64_t mask = ALL >> (64 - k), shift = k - 1;
    code->x[0] = ((code->x[0]<<1) | (c&1))  & mask;
    code->x[1] = ((code->x[1]<<1) | (c>>1)) & mask;
	code->x[2] = code->x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
	code->x[3] = code->x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
}

#endif
