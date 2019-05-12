#include "Assembly.h"
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "Process_Read.h"
#include "CommandLines.h"
#include "kmer.h"
#include "Hash_Table.h"

Total_Count_Table TCB;



/********************************for debug***************************************/
void Verify_Counting()
{


    kseq_t *seq = (kseq_t*)calloc(1, sizeof(kseq_t));	


    long long read_number = 0;

    HPC_seq HPC_read;
    uint64_t code;
    Hash_code k_code;
    int avalible_k = 0;

    fprintf(stdout, "Start Verifying ...\n");

    while (get_read(seq))
    {

        init_HPC_seq(&HPC_read, seq->seq.s, seq->seq.l);
        init_Hash_code(&k_code);

        avalible_k = 0;

        while ((code = get_HPC_code(&HPC_read)) != 6)
        {
            if(code < 4)
            {
                k_mer_append(&k_code,code,k_mer_length);
                avalible_k++;
                if (avalible_k>=k_mer_length)
                {
                    if(verify_Total_Count_Table(&TCB, &k_code, k_mer_length) == -1)
                    {
                        fprintf(stderr, "ERROR when subtracting!\n");
                    }
                }
                
            }
            else
            {
                avalible_k = 0;
                init_Hash_code(&k_code);
            }
            
        }




        read_number++;
    }
    
    fprintf(stdout, "read_number: %lld\n",read_number);

    fprintf(stdout, "Start Traversing ...\n");

    Traverse_Total_Count_Table(&TCB);

    fprintf(stdout, "Finish Traversing!\n");


}




void* Perform_Counting(void* arg)
{
    int thr_ID = *((int*)arg);

    ///debug_mode(101, thr_ID, thread_num);


    int i = 0;
    HPC_seq HPC_read;

    R_buffer_block curr_sub_block;

    init_R_buffer_block(&curr_sub_block);

    long long read_number = 0;
    long long select_k_mer_number = 0 ;
    long long k_mer_number = 0 ;

	int file_flag = 1;

    uint64_t code;

    Hash_code k_code;

    int avalible_k = 0;


	while (file_flag != 0)
	{

        file_flag = get_reads_mul_thread(&curr_sub_block);

        read_number = read_number + curr_sub_block.num; 

        for (i = 0; i < curr_sub_block.num; i++)
        {

            init_HPC_seq(&HPC_read, curr_sub_block.read[i].seq.s, curr_sub_block.read[i].seq.l);
            init_Hash_code(&k_code);

            avalible_k = 0;

            while ((code = get_HPC_code(&HPC_read)) != 6)
            {
                if(code < 4)
                {
                    k_mer_append(&k_code,code,k_mer_length);
                    avalible_k++;
                    if (avalible_k>=k_mer_length)
                    {
                        ///插入
                        if(insert_Total_Count_Table(&TCB, &k_code, k_mer_length))
                        {
                            select_k_mer_number++;
                        }

                        k_mer_number++;
                        
                    }
                    
                }
                else
                {
                    avalible_k = 0;
                    init_Hash_code(&k_code);
                }
                
            }
            

        }
        

    }

    destory_R_buffer_block(&curr_sub_block);
    free(arg);

    fprintf(stdout, "thr_ID: %d, read_number: %lld, select_k_mer_number: %lld, k_mer_number: %lld\n",
    thr_ID, read_number, select_k_mer_number, k_mer_number);

}
















void Counting_multiple_thr()
{

    init_Total_Count_Table(k_mer_length, &TCB);

    fprintf(stdout, "TCB.prefix_bits: %d\n", TCB.prefix_bits);
    fprintf(stdout, "TCB.suffix_bits: %d\n", TCB.suffix_bits);
    fprintf(stdout, "TCB.size: %d\n", TCB.size);



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

    destory_R_buffer();


    destory_Total_Count_Table(&TCB);
    

}























