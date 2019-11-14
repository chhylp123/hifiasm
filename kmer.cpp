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



void init_small_hash_table(small_hash_table* x)
{
    x->size = 0;
    x->buffer = NULL;
    x->length = 0;
}

void clear_small_hash_table(small_hash_table* x)
{
    x->length = 0;
}

void resize_small_hash_table(small_hash_table* x, uint64_t size)
{
    if(size > x->size)
    {
        x->size = size;
        x->buffer = (k_v*)realloc(x->buffer, x->size*sizeof(k_v));
    }
}

void destory_small_hash_table(small_hash_table* x)
{
    free(x->buffer);
}


void add_small_hash_table(small_hash_table* x, k_v* element)
{
    if(x->length + 1 > x->size)
    {
        x->size = (x->length + 1) * 2;
        x->buffer = (k_v*)realloc(x->buffer, x->size*sizeof(k_v));
    }

    x->buffer[x->length] = (*element);
    x->length++;
}

//x > y, return 1; x < y, return -1, x == y, return 0
int compare_k_mer(k_v* x, k_v* y)
{
    if(x->key.x[1] != y->key.x[1])
    {
        return x->key.x[1] > y->key.x[1] ? 1: -1;
    }
    else
    {
        if(x->key.x[0] != y->key.x[0])
        {
            return x->key.x[0] > y->key.x[0] ? 1: -1;
        }
        else
        {
            return 0;
        }
    }
    
}

int cmp_k_mer_kv(const void * a, const void * b)
{
    int flag = compare_k_mer((k_v*)a, (k_v*)b);

    if(flag == 0)
    {
        if ((*(k_v*)a).value != (*(k_v*)b).value)
        {
            return (*(k_v*)a).value > (*(k_v*)b).value ? 1: -1;
        }
        else
        {
            return 0;
        }
        
    }
    else
    {
        return flag;
    }
}

void sort_small_hash_table(small_hash_table* x)
{
    qsort(x->buffer, x->length, sizeof(k_v), cmp_k_mer_kv);
}


inline long long firstEqual(k_v* arr, long long arrLen, k_v* key)
{
    long long L = 0, R = arrLen - 1; //[L, R]
    long long mid;
    int flag;
    while( L <= R)
    {
        mid = L + (R - L)/2;

        flag = compare_k_mer(&(arr[mid]), key);

        ///arr[mid] >= key
        if(flag >= 0)
        {
            R = mid - 1;
        }
        else
        {
            L = mid + 1;
        }
    }

    
    if(L < arrLen && (flag = compare_k_mer(&(arr[L]), key) == 0))
    {
        return L;
    }
        
    return -1;
}

inline long long lastEqual(k_v* arr, long long arrLen, k_v* key)
{
    long long L = 0, R = arrLen - 1; //[L, R]
    long long mid;
    int flag;
    while( L <= R)
    {
        mid = L + (R - L)/2;
        flag = compare_k_mer(&(arr[mid]), key);
        ///arr[mid] <= key
        if(flag <= 0)
        {
            L = mid + 1;
        }
        else
        {
            R = mid - 1;
        }
    }
    
    if(R >= 0 && ((flag = compare_k_mer(&(arr[R]), key)) == 0))
    {
        return R;
    }
        
    return -1;
}

int query_small_hash_table(small_hash_table* target, k_v* query, long long* l_end, long long* r_end)
{
    (*l_end) = -1;
    (*r_end) = -1;
    long long left_end;
    long long right_end;

    left_end = firstEqual(target->buffer, target->length, query);

    if(left_end != -1)
    {
        right_end = lastEqual(target->buffer + left_end, target->length - left_end, query) + left_end;

        (*l_end) = left_end;
        (*r_end) = right_end;
        
        if(right_end == -1)
        {
            fprintf(stderr, "error\n");
        }
        

        return right_end - left_end + 1;
    }

    return 0;
}

