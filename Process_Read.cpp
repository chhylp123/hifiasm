#include "Process_Read.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <pthread.h>


gzFile fp;
kseq_t *seq;
R_buffer RDB;

pthread_mutex_t i_readinputMutex;
pthread_mutex_t i_queueMutex;
pthread_mutex_t i_terminateMutex;
pthread_cond_t i_flushCond;
pthread_cond_t i_readinputflushCond;
pthread_cond_t i_stallCond;
pthread_cond_t i_readinputstallCond;
pthread_mutex_t i_doneMutex;


void init_kseq(char* file)
{
	fp = gzopen(file, "r");
  	seq = kseq_init(fp);
}

void destory_kseq()
{
	kseq_destroy(seq);
  	gzclose(fp);
}


inline void exchage_kstring_t(kstring_t* a, kstring_t* b)
{
	kstring_t tmp;
	tmp = *a;
	*a = *b;
	*b = tmp;
}

int get_read(kseq_t *s)
{
	int l;

	if ((l = kseq_read(seq)) >= 0)
	{
		
		exchage_kstring_t(&seq->comment, &s->comment);
		exchage_kstring_t(&seq->name, &s->name);
		exchage_kstring_t(&seq->qual, &s->qual);
		exchage_kstring_t(&seq->seq, &s->seq);
		
		return 1;
	}
	else
	{
		return 0;
	}
	
	
}

void init_R_buffer_block(R_buffer_block* curr_sub_block)
{
	curr_sub_block->read = (kseq_t*)calloc(RDB.block_inner_size, sizeof(kseq_t));	
	curr_sub_block->num = 0;
}


void init_R_buffer(int thread_num)
{
	RDB.all_read_end = 0;
	RDB.num = 0;
	RDB.block_inner_size = READ_BLOCK_SIZE;
	RDB.size = thread_num*READ_BLOCK_NUM_PRE_THR;

	RDB.sub_block = (R_buffer_block*)malloc(sizeof(R_buffer_block)*RDB.size);

	int i = 0;

	for (i = 0; i < RDB.size; i++)
	{
		init_R_buffer_block(&RDB.sub_block[i]);
	}
	
}


void destory_R_buffer_block(R_buffer_block* curr_sub_block)
{

	free(curr_sub_block->read);
}


void destory_R_buffer()
{
	int i = 0;

	for (i = 0; i < RDB.size; i++)
	{
		destory_R_buffer_block(&RDB.sub_block[i]);
	}

	free(RDB.sub_block);
	
}


inline void load_read_block(R_buffer_block* read_batch, int batch_read_size,
	int* return_file_flag)
{
	int inner_i = 0;
	int file_flag = 1;




	while (inner_i<batch_read_size)
	{

		file_flag = get_read(&read_batch->read[inner_i]);

		if (file_flag == 1)
		{
			inner_i++;
		}
		else if (file_flag == 0)
		{
			break;
		}
	}

	if (inner_i || file_flag)
	{
		file_flag = 1;
	}

	*return_file_flag = file_flag;
	read_batch->num = inner_i;

}


inline void push_R_block(R_buffer_block* tmp_sub_block)
{
	

	///only exchange pointers
	kseq_t *k1;
	k1 = RDB.sub_block[RDB.num].read;

	RDB.sub_block[RDB.num].read = tmp_sub_block->read;

	tmp_sub_block->read = k1;

	RDB.sub_block[RDB.num].num = tmp_sub_block->num;
	tmp_sub_block->num = 0;

	RDB.num++;
}


inline void pop_R_block(R_buffer_block* curr_sub_block)
{
	RDB.num--;
	
	///only exchange pointers
	kseq_t *k1;
	k1 = RDB.sub_block[RDB.num].read;

	RDB.sub_block[RDB.num].read = curr_sub_block->read;

	curr_sub_block->read = k1;

	curr_sub_block->num = RDB.sub_block[RDB.num].num;
	RDB.sub_block[RDB.num].num = 0;



}


void* input_reads_muti_threads(void*)
{


	int i = 0;
	int file_flag = 1;

	R_buffer_block tmp_buf;

	init_R_buffer_block(&tmp_buf);



	while (1)
	{



		load_read_block(&tmp_buf, RDB.block_inner_size, &file_flag);

		if (file_flag == 0)
		{
			break;
		}


		pthread_mutex_lock(&i_readinputMutex);
		while (IS_FULL(RDB))
		{

			pthread_cond_signal(&i_readinputstallCond);
			pthread_cond_wait(&i_readinputflushCond, &i_readinputMutex);
		}


		push_R_block(&tmp_buf);

		pthread_cond_signal(&i_readinputstallCond);
		pthread_mutex_unlock(&i_readinputMutex);
	}


	pthread_mutex_lock(&i_readinputMutex);
	RDB.all_read_end = 1;
	pthread_cond_signal(&i_readinputstallCond);  //important
	pthread_mutex_unlock(&i_readinputMutex);

	destory_R_buffer_block(&tmp_buf);

}



int get_reads_mul_thread(R_buffer_block* curr_sub_block)
{


	pthread_mutex_lock(&i_readinputMutex);


	while (IS_EMPTY(RDB) && RDB.all_read_end == 0)
	{

		pthread_cond_signal(&i_readinputflushCond);
		pthread_cond_wait(&i_readinputstallCond, &i_readinputMutex);
	}


	if (!IS_EMPTY(RDB))
	{
		pop_R_block(curr_sub_block);
		pthread_cond_signal(&i_readinputflushCond);
		pthread_mutex_unlock(&i_readinputMutex);


		return 1;
	}
	else
	{
		curr_sub_block->num = 0;

		pthread_cond_signal(&i_readinputstallCond);   //important

		pthread_mutex_unlock(&i_readinputMutex);

		return 0;
	}
	

}



void Counting_block()
{

	long long read_number = 0;
    int i = 0;
	int file_flag = 1;

	R_buffer_block tmp_buf;

	init_R_buffer_block(&tmp_buf);

	while (1)
	{


		load_read_block(&tmp_buf, RDB.block_inner_size,
			&file_flag);


		if (file_flag == 0)
		{
			break;
		}

		for (i = 0; i < tmp_buf.num; i++)
		{
			fprintf(stderr,"@%s\n", tmp_buf.read[i].name.s);
			fprintf(stderr,"%s\n",tmp_buf.read[i].seq.s);
			fprintf(stderr,"+\n");
			fprintf(stderr,"%s\n",tmp_buf.read[i].qual.s);

			read_number++;
		}
		
	}

    fprintf(stdout, "read_number: %lld\n",read_number);



}