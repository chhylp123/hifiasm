#include <stdio.h>
#include <stdlib.h>
#include "kmer.h"

void init_HPC_seq(HPC_seq* seq, char* str, long long l)
{
    seq->i = 0;
    seq->l = l;
    seq->N_occ = 0;
    seq->str = str;
}

void init_Hash_code(Hash_code* code)
{
    code->x[0] = 0;
    code->x[1] = 0;
}