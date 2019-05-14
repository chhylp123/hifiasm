#ifndef __READ__
#define __READ__

#include<stdint.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
#include "kseq.h"

#define READ_INIT_NUMBER 1000

#define READ_BLOCK_SIZE 64
#define READ_BLOCK_NUM_PRE_THR 100

#define IS_FULL(buffer) ((buffer.num >= buffer.size)?1:0)
#define IS_EMPTY(buffer) ((buffer.num == 0)?1:0)
#define Get_READ_LENGTH(R_INF, ID) (R_INF.index[ID+1] - R_INF.index[ID])
#define Get_NAME_LENGTH(R_INF, ID) (R_INF.name_index[ID+1] - R_INF.name_index[ID])
#define Get_READ(R_INF, ID) (R_INF.read + R_INF.index[ID])
#define Get_NAME(R_INF, ID) (R_INF.name + R_INF.name_index[ID])



KSEQ_INIT(gzFile, gzread)

void init_kseq(char* file);
void destory_kseq();
int get_read(kseq_t *s);



typedef struct
{
	char* read;
	char* name;
	uint64_t* index;
	uint64_t index_size;
	uint64_t* name_index;
	uint64_t name_index_size;
	uint64_t total_reads;
	uint64_t total_reads_bases;
	uint64_t total_name_length;

} All_reads;

extern All_reads R_INF;

void malloc_All_reads(All_reads* r);

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
void init_All_reads(All_reads* r);
void* input_reads_muti_threads(void*);
void init_R_buffer_block(R_buffer_block* curr_sub_block);
int get_reads_mul_thread(R_buffer_block* curr_sub_block);

void Counting_block();
void destory_R_buffer_block(R_buffer_block* curr_sub_block);
void destory_R_buffer();
void clear_R_buffer();


#endif
