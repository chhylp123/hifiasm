#include "Assembly.h"
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "Process_Read.h"
#include "CommandLines.h"
#include "kmer.h"
#include "Hash_Table.h"

Total_Count_Table TCB;
Total_Pos_Table PCB;
All_reads R_INF;




void* Perform_Counting(void* arg)
{
    int thr_ID = *((int*)arg);

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

            ///forward strand
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


            /**
            ///reverse complement strand
            reverse_complement(curr_sub_block.read[i].seq.s, curr_sub_block.read[i].seq.l);
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
            **/


        }
        

    }

    destory_R_buffer_block(&curr_sub_block);
    free(arg);

}






void* Build_hash_table(void* arg)
{
    int thr_ID = *((int*)arg);

    int i = 0;
    HPC_seq HPC_read;

    R_buffer_block curr_sub_block;

    init_R_buffer_block(&curr_sub_block);

	int file_flag = 1;

    uint64_t code;

    long long HPC_base;

    Hash_code k_code;

    int avalible_k = 0;


	while (file_flag != 0)
	{

        file_flag = get_reads_mul_thread(&curr_sub_block);


        for (i = 0; i < curr_sub_block.num; i++)
        {
 
            ///forward strand
            init_HPC_seq(&HPC_read, curr_sub_block.read[i].seq.s, curr_sub_block.read[i].seq.l);
            init_Hash_code(&k_code);

            avalible_k = 0;

            HPC_base = 0;

            while ((code = get_HPC_code(&HPC_read)) != 6)
            {
                if(code < 4)
                {
                    k_mer_append(&k_code,code,k_mer_length);
                    avalible_k++;
                    if (avalible_k>=k_mer_length)
                    {

                       ///选取的k-mer满足两个要求
                       ///1. hash(k-mer) % 101 <= 3
                       ///2. occ(k-mer)要满足范围
                       ///TCB表中的元素仅满足第一个要求，而PCB表中的元素满足两个要求
                       ///所以如果当前k-mer在PCB表中存在，则他的位置一定要加入到候选位置中去
                        ///insert_Total_Pos_Table(&PCB, &k_code, k_mer_length, curr_sub_block.read[i].ID, HPC_base - k_mer_length + 1, FORWARD);
                        insert_Total_Pos_Table(&PCB, &k_code, k_mer_length, 
                            curr_sub_block.read[i].ID, HPC_base - k_mer_length + 1);
                    }
                    
                }
                else
                {
                    avalible_k = 0;
                    init_Hash_code(&k_code);
                }

                HPC_base++;
            }

            
            ///load read
            compress_base(Get_READ(R_INF, curr_sub_block.read[i].ID),
            curr_sub_block.read[i].seq.s, curr_sub_block.read[i].seq.l, 
            &R_INF.N_site[curr_sub_block.read[i].ID], HPC_read.N_occ);

            memcpy(R_INF.name+R_INF.name_index[curr_sub_block.read[i].ID], 
            curr_sub_block.read[i].name.s, curr_sub_block.read[i].name.l);

            /**
            ///reverse complement strand
            reverse_complement(curr_sub_block.read[i].seq.s, curr_sub_block.read[i].seq.l);

            init_HPC_seq(&HPC_read, curr_sub_block.read[i].seq.s, curr_sub_block.read[i].seq.l);
            init_Hash_code(&k_code);

            avalible_k = 0;

            HPC_base = 0;

            while ((code = get_HPC_code(&HPC_read)) != 6)
            {
                if(code < 4)
                {
                    k_mer_append(&k_code,code,k_mer_length);
                    avalible_k++;
                    if (avalible_k>=k_mer_length)
                    {

                       ///选取的k-mer满足两个要求
                       ///1. hash(k-mer) % 101 <= 3
                       ///2. occ(k-mer)要满足范围
                       ///TCB表中的元素仅满足第一个要求，而PCB表中的元素满足两个要求
                       ///所以如果当前k-mer在PCB表中存在，则他的位置一定要加入到候选位置中去
                        insert_Total_Pos_Table(&PCB, &k_code, k_mer_length, 
                            curr_sub_block.read[i].ID, HPC_base - k_mer_length + 1, REVERSE_COMPLEMENT);
                    }
                    
                }
                else
                {
                    avalible_k = 0;
                    init_Hash_code(&k_code);
                }

                HPC_base++;
            }
            **/

        }
        

    }

    destory_R_buffer_block(&curr_sub_block);
    free(arg);



}









void Counting_multiple_thr()
{

    double start_time = Get_T();

    fprintf(stdout, "Begin Counting ...... \n");

    init_kseq(read_file_name);

    init_All_reads(&R_INF);

    init_Total_Count_Table(k_mer_length, &TCB);

    pthread_t inputReadsHandle;
    
    init_R_buffer(thread_num);

    int *is_insert = (int*)malloc(sizeof(*is_insert));
	*is_insert = 1;

    pthread_create(&inputReadsHandle, NULL, input_reads_muti_threads, (void*)is_insert);

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

    ///destory_R_buffer();
    
    ///destory_Total_Count_Table(&TCB);

    destory_kseq();

    fprintf(stdout, "Finish Counting ...... \n");

    fprintf(stdout, "%-30s%18.2f\n\n", "Counting time:", Get_T() - start_time);
   

}




void Build_hash_table_multiple_thr()
{
    double start_time = Get_T();


    fprintf(stdout, "Begin Building hash table ...... \n");

    init_Total_Pos_Table(&PCB, &TCB);

    double T_start_time = Get_T();

    Traverse_Counting_Table(&TCB, &PCB, k_mer_min_freq, k_mer_max_freq);

    fprintf(stdout, "%-30s%18.2f\n\n", "Traverse time:", Get_T() - T_start_time);
    fflush(stdout);
   


    init_kseq(read_file_name);

    clear_R_buffer();

    malloc_All_reads(&R_INF);

    pthread_t inputReadsHandle;

    int *is_insert = (int*)malloc(sizeof(*is_insert));
	*is_insert = 0;

    pthread_create(&inputReadsHandle, NULL, input_reads_muti_threads, (void*)is_insert);




    pthread_t *_r_threads;

	_r_threads = (pthread_t *)malloc(sizeof(pthread_t)*thread_num);

    int i = 0;

	for (i = 0; i < thread_num; i++)
	{
		int *arg = (int*)malloc(sizeof(*arg));
		*arg = i;

		pthread_create(_r_threads + i, NULL, Build_hash_table, (void*)arg);

	}

    pthread_join(inputReadsHandle, NULL);

    for (i = 0; i<thread_num; i++)
		pthread_join(_r_threads[i], NULL);

    free(_r_threads);

    destory_kseq();

    fprintf(stdout, "Finish Building hash table ...... \n");

    fprintf(stdout, "%-30s%18.2f\n\n", "Build hash table time:", Get_T() - start_time);

    ///destory_Total_Count_Table(&TCB);
    if (write_index_to_disk)
    {
        write_Total_Pos_Table(&PCB, read_file_name);
        destory_Total_Pos_Table(&PCB);
        load_Total_Pos_Table(&PCB, read_file_name);

        write_All_reads(&R_INF, read_file_name);
        destory_All_reads(&R_INF);
        load_All_reads(&R_INF, read_file_name);
    }
    

    
}




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



/********************************for debug*****************************************/
void verify_Position_hash_table()
{

    init_kseq(read_file_name);

    kseq_t *seq = (kseq_t*)calloc(1, sizeof(kseq_t));	


    long long read_number = 0;
    long long HPC_base;

    HPC_seq HPC_read;
    uint64_t code;
    Hash_code k_code;
    int avalible_k = 0;
    k_mer_pos* list;
    uint64_t sub_ID;

    UC_Read g_read;
    init_UC_Read(&g_read);

    fprintf(stdout, "Start Verifying Position Table...\n");

    while (get_read(seq))
    {

        if (seq->seq.l
            != Get_READ_LENGTH(R_INF, read_number))
        {
            fprintf(stderr, "seq error\n");
        }

        recover_UC_Read(&g_read, &R_INF, read_number);


        if(memcmp(seq->seq.s, g_read.seq, seq->seq.l))
        {
            fprintf(stderr, "\nseq error ID: %llu, length: %llu\n",read_number, seq->seq.l);
            
        }

        if (seq->name.l
            != Get_NAME_LENGTH(R_INF, read_number))
        {
            fprintf(stderr, "name error\n");
        }


         if(memcmp(seq->name.s, Get_NAME(R_INF, read_number), seq->name.l))
        {
            fprintf(stderr, "name error\n");
        }



        init_HPC_seq(&HPC_read, seq->seq.s, seq->seq.l);
        init_Hash_code(&k_code);

        avalible_k = 0;
        HPC_base = 0;

        while ((code = get_HPC_code(&HPC_read)) != 6)
        {
            if(code < 4)
            {
                k_mer_append(&k_code,code,k_mer_length);
                avalible_k++;
                if (avalible_k>=k_mer_length)
                {

                    uint64_t count1 = get_Total_Count_Table(&TCB, &k_code, k_mer_length);
                    uint64_t count2 = count_Total_Pos_Table(&PCB, &k_code, k_mer_length);
                    if(count1>=k_mer_min_freq && count1<= k_mer_max_freq && count1!=count2)
                    {
                        fprintf(stderr, "count1: %lld\n",count1);
                        fprintf(stderr, "count2: %lld\n",count2);
                    }
                    

                    if(locate_Total_Pos_Table(&PCB, &k_code, &list, k_mer_length, &sub_ID) == count2)
                    {
                        int j = 0;
                        for(j=1; j<count2; j++)
                        {
                            if(list[j].readID < list[j - 1].readID)
                            {
                                fprintf(stderr, "locate error\n");
                            }
                            else if(list[j].readID == list[j - 1].readID)
                            {
                                if(list[j].offset < list[j - 1].offset)
                                {
                                    fprintf(stderr, "locate error\n");
                                }
                            }
                        }

                        if(count2 != 0)
                        {
                            for(j=0; j<count2; j++)
                            {   
                                if(list[j].readID == read_number && list[j].offset == HPC_base - k_mer_length + 1)
                                {
                                    break;
                                }
                            }

                            if(j == count2)
                            {
                                fprintf(stderr, "locate error, read_number: %llu, pos: %llu\n", 
                                read_number, HPC_base - k_mer_length + 1);

                                for(j=0; j<count2; j++)
                                {   

                                    fprintf(stderr, "j: %llu, readID: %llu, offset: %llu\n", j, list[j].readID, list[j].offset);
                            
                                }
                            }
                        }
                        
                        
                    }
                    else
                    {
                        fprintf(stderr, "locate error\n");
                    }

                }
                
            }
            else
            {
                avalible_k = 0;
                init_Hash_code(&k_code);
            }

            HPC_base++;
            
        }

        read_number++;
    }

    destory_kseq();

    fprintf(stdout, "Finish Verifying Position Table!\n");

}












