#include "Assembly.h"
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "Process_Read.h"
#include "CommandLines.h"




void Counting()
{


    kseq_t *seq = (kseq_t*)calloc(1, sizeof(kseq_t));	


    long long read_number = 0;

    while (get_read(seq))
    {

        fprintf(stderr,"@%s\n",seq->name.s);
        fprintf(stderr,"%s\n",seq->seq.s);
        fprintf(stderr,"+\n");
        fprintf(stderr,"%s\n",seq->qual.s);

        read_number++;
    }
    
    fprintf(stdout, "read_number: %lld\n",read_number);



}















void* Perform_Counting(void* arg)
{
    R_buffer_block curr_sub_block;

    init_R_buffer_block(&curr_sub_block);

    long long read_number = 0;
    
	int file_flag = 1;
	while (file_flag != 0)
	{

        file_flag = get_reads_mul_thread(&curr_sub_block);

        read_number = read_number + curr_sub_block.num; 

    }

    fprintf(stdout, "#########read_number: %lld\n",read_number);
    fflush(stdout);

}

















void Counting_multiple_thr()
{
    pthread_t inputReadsHandle;
    
    init_R_buffer(thread_num);

    pthread_create(&inputReadsHandle, NULL, input_reads_muti_threads, NULL);


    pthread_t *_r_threads;

	_r_threads = (pthread_t *)malloc(sizeof(pthread_t)*thread_num);

    int i = 0;

	for (i = 0; i < thread_num; i++)
	{
		int *arg = (int*)malloc(sizeof(*arg));
		*arg = i;

		pthread_create(_r_threads + i, NULL, Perform_Counting, (void*)arg);

	}



    pthread_join(inputReadsHandle, NULL);

    for (i = 0; i<thread_num; i++)
		pthread_join(_r_threads[i], NULL);

    free(_r_threads);


}















