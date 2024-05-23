#include "Output.h"
#include "CommandLines.h"
#include <stdio.h>

pthread_mutex_t o_queueMutex;
pthread_cond_t o_flushCond;
pthread_cond_t o_stallCond;
pthread_mutex_t o_doneMutex;


Output_buffer buffer_out;
Output_buffer_sub_block tmp_buffer_sub_block;

void init_buffer_sub_block(Output_buffer_sub_block* sub_block)
{
	sub_block->length = 0;
	sub_block->size = SUB_BLOCK_INIT_SIZE;
	sub_block->buffer = (char*)malloc(sub_block->size);
}

void destory_buffer_sub_block(Output_buffer_sub_block* sub_block)
{
    free(sub_block->buffer);
}

void destory_output_buffer()
{
	for (int i = 0; i < buffer_out.sub_block_size; i++)
	{
		destory_buffer_sub_block(&(buffer_out.sub_buffer[i]));
	}

    free(buffer_out.sub_buffer);
}

void init_output_buffer(int thread_number)
{
	buffer_out.sub_block_size = OUTPUT_BUFFER_SIZE * thread_number;
	buffer_out.sub_block_number = 0;

	buffer_out.sub_buffer = (Output_buffer_sub_block*)malloc(sizeof(Output_buffer_sub_block)*buffer_out.sub_block_size);



	for (int i = 0; i < buffer_out.sub_block_size; i++)
	{
		init_buffer_sub_block(&(buffer_out.sub_buffer[i]));
	}

	buffer_out.all_buffer_end = 0;
}

inline int if_empty_buffer()
{

	if (buffer_out.sub_block_number == 0)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}


inline int if_full_buffer()
{

	if (buffer_out.sub_block_number >= buffer_out.sub_block_size)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

inline void pop_single_buffer(Output_buffer_sub_block* curr_sub_block)
{
	buffer_out.sub_block_number--;
	char *k;
	k = buffer_out.sub_buffer[buffer_out.sub_block_number].buffer;
	buffer_out.sub_buffer[buffer_out.sub_block_number].buffer = curr_sub_block->buffer;
	curr_sub_block->buffer = k;

	
	long long tmp_size;
	tmp_size = curr_sub_block->size;
	curr_sub_block->size = buffer_out.sub_buffer[buffer_out.sub_block_number].size;
	buffer_out.sub_buffer[buffer_out.sub_block_number].size = tmp_size;


	curr_sub_block->length = buffer_out.sub_buffer[buffer_out.sub_block_number].length;
	buffer_out.sub_buffer[buffer_out.sub_block_number].length = 0;
}



void add_segment_to_sub_buffer(Output_buffer_sub_block* current_sub_buffer, char* seg, long long segLen)
{
    if(current_sub_buffer->length + segLen + 2 > current_sub_buffer->size)
    {
        current_sub_buffer->size = current_sub_buffer->length + segLen + 2;
        current_sub_buffer->buffer = (char*)realloc(current_sub_buffer->buffer, current_sub_buffer->size);
    }

    memcpy(current_sub_buffer->buffer + current_sub_buffer->length, seg, segLen);
    current_sub_buffer->length += segLen;
    current_sub_buffer->buffer[current_sub_buffer->length] = '\0';
}


void add_base_to_sub_buffer(Output_buffer_sub_block* current_sub_buffer, char base)
{
    if(current_sub_buffer->length + 2 > current_sub_buffer->size)
    {
        current_sub_buffer->size = current_sub_buffer->length + 2;
        current_sub_buffer->buffer = (char*)realloc(current_sub_buffer->buffer, current_sub_buffer->size);
    }

    current_sub_buffer->buffer[current_sub_buffer->length] = base;
    current_sub_buffer->length++;
    current_sub_buffer->buffer[current_sub_buffer->length] = '\0';

}


inline void push_single_buffer(Output_buffer_sub_block* curr_sub_block)
{

	char *k;
	k = buffer_out.sub_buffer[buffer_out.sub_block_number].buffer;
	buffer_out.sub_buffer[buffer_out.sub_block_number].buffer = curr_sub_block->buffer;
	curr_sub_block->buffer = k;

	long long tmp_size;
	tmp_size = curr_sub_block->size;
	curr_sub_block->size = buffer_out.sub_buffer[buffer_out.sub_block_number].size;
	buffer_out.sub_buffer[buffer_out.sub_block_number].size = tmp_size;

	buffer_out.sub_buffer[buffer_out.sub_block_number].length = curr_sub_block->length;
	curr_sub_block->length = 0;

	buffer_out.sub_block_number++;
}


void* pop_buffer(void*)
{

	FILE* output_file = fopen(asm_opt.output_file_name, "w");

	init_buffer_sub_block(&tmp_buffer_sub_block);

	

	while (buffer_out.all_buffer_end < asm_opt.thread_num)
	{

		pthread_mutex_lock(&o_queueMutex);

		while (if_empty_buffer() && (buffer_out.all_buffer_end < asm_opt.thread_num))
		{
			pthread_cond_signal(&o_stallCond);
			pthread_cond_wait(&o_flushCond, &o_queueMutex);
		}

		if (!if_empty_buffer())
		{
			pop_single_buffer(&tmp_buffer_sub_block);
		}

		pthread_cond_signal(&o_stallCond);
		pthread_mutex_unlock(&o_queueMutex);

        

		if (tmp_buffer_sub_block.length != 0)
		{
			fprintf(output_file, "%s", tmp_buffer_sub_block.buffer);
		}

	}


	while (buffer_out.sub_block_number>0)
	{
		buffer_out.sub_block_number--;
		fprintf(output_file, "%s", buffer_out.sub_buffer[buffer_out.sub_block_number].buffer);
	}


    destory_buffer_sub_block(&tmp_buffer_sub_block);

    fclose(output_file);

	return NULL;
}


void push_results_to_buffer(Output_buffer_sub_block* sub_block)
{
	
	pthread_mutex_lock(&o_queueMutex);

	while (if_full_buffer())
	{
		pthread_cond_signal(&o_flushCond);
		pthread_cond_wait(&o_stallCond, &o_queueMutex);
	}

	push_single_buffer(sub_block);
	pthread_cond_signal(&o_flushCond);
	pthread_mutex_unlock(&o_queueMutex);
	
}



void finish_output_buffer()
{

    pthread_mutex_lock(&o_doneMutex);
    buffer_out.all_buffer_end++;


    if (buffer_out.all_buffer_end == asm_opt.thread_num)
    {
        pthread_cond_signal(&o_flushCond);
    }
        
    pthread_mutex_unlock(&o_doneMutex);
}