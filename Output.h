#ifndef __OUTPUT__
#define __OUTPUT__

#include<stdint.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>

typedef struct
{
	char* buffer;
	long long size;
	long long length;
} Output_buffer_sub_block;

typedef struct Output_buffer
{
	Output_buffer_sub_block* sub_buffer;
	long long sub_block_size;
	long long sub_block_number;
	int all_buffer_end;
} Output_buffer;

#define OUTPUT_BUFFER_SIZE 100
#define SUB_BLOCK_INIT_SIZE 10000

void init_buffer_sub_block(Output_buffer_sub_block* sub_block);
void* pop_buffer(void*);
void add_segment_to_sub_buffer(Output_buffer_sub_block* current_sub_buffer, char* seg, long long segLen);
void add_base_to_sub_buffer(Output_buffer_sub_block* current_sub_buffer, char base);
void push_results_to_buffer(Output_buffer_sub_block* sub_block);
void finish_output_buffer();
void destory_buffer_sub_block(Output_buffer_sub_block* sub_block);
void destory_output_buffer();
void init_output_buffer(int thread_number);

#endif
