#ifndef __READ__
#define __READ__

#include<stdint.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
#include "kseq.h"

#define READ_BLOCK_SIZE 64
#define READ_BLOCK_NUM_PRE_THR 32

#define IS_FULL(buffer) ((buffer.num >= buffer.size)?1:0)
#define IS_EMPTY(buffer) ((buffer.num == 0)?1:0)



KSEQ_INIT(gzFile, gzread)

void init_kseq(char* file);
void destory_kseq();
int get_read(kseq_t *s);



typedef struct
{
	kseq_t* read;
	long long num;

} R_buffer_block;


typedef struct
{
	R_buffer_block* sub_block;
	long long block_inner_size;
	long long size;
	long long num;
	int all_read_end;
} R_buffer;

void init_R_buffer(int thread_num);
void* input_reads_muti_threads(void*);
void init_R_buffer_block(R_buffer_block* curr_sub_block);
int get_reads_mul_thread(R_buffer_block* curr_sub_block);

void Counting_block();

#endif
