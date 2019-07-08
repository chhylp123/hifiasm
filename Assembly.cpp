#include "Assembly.h"
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "Process_Read.h"
#include "CommandLines.h"
#include "kmer.h"
#include "Hash_Table.h"
#include "POA.h"
#include "Correct.h"

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

    uint64_t end_pos;

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

            while ((code = get_HPC_code(&HPC_read, &end_pos)) != 6)
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

            while ((code = get_HPC_code(&HPC_read, &end_pos)) != 6)
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
    uint64_t end_pos;

    ///long long HPC_base;

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

            ///HPC_base = 0;

            while ((code = get_HPC_code(&HPC_read, &end_pos)) != 6)
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
                            curr_sub_block.read[i].ID, end_pos);
                    }
                    
                }
                else
                {
                    avalible_k = 0;
                    init_Hash_code(&k_code);
                }

                ///HPC_base++;
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

            while ((code = get_HPC_code(&HPC_read, &end_pos)) != 6)
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
        ///destory_Total_Pos_Table(&PCB);
        ///load_Total_Pos_Table(&PCB, read_file_name);

        write_All_reads(&R_INF, read_file_name);
        ///destory_All_reads(&R_INF);
        ///load_All_reads(&R_INF, read_file_name);
    }
    

    
}


void* Overlap_calculate_version1(void* arg)
{

    int thr_ID = *((int*)arg);

    long long i = 0;
    int avalible_k = 0;

    UC_Read g_read;
    init_UC_Read(&g_read);
    HPC_seq HPC_read;
    Hash_code k_code;
    uint64_t code;
    uint64_t end_pos;
    k_mer_pos* list;
    uint64_t list_length;
    uint64_t sub_ID;
    Candidates_list l;

    k_mer_pos_list* merge_list;

    init_Candidates_list(&l);

    for (i = thr_ID; i < R_INF.total_reads; i = i + thread_num)
    ///for (i = thr_ID; i < R_INF.total_reads/50; i = i + thread_num)
    {
        /**
        if (i % 1000 == 0)
        {
            fprintf(stderr, "i: %llu\n", i);
        }
        **/
        
        

        clear_Candidates_list(&l);
        recover_UC_Read(&g_read, &R_INF, i);


        ///forward strand
        init_HPC_seq(&HPC_read, g_read.seq, g_read.length);
        init_Hash_code(&k_code);
        avalible_k = 0;

        while ((code = get_HPC_code(&HPC_read, &end_pos)) != 6)
        {
            if(code < 4)
            {
                k_mer_append(&k_code,code,k_mer_length);
                avalible_k++;
                if (avalible_k>=k_mer_length)
                {
                    list_length = locate_Total_Pos_Table(&PCB, &k_code, &list, k_mer_length, &sub_ID);
                    merge_Candidates_list(&l, list, list_length, end_pos, 0);
                    //merge_Candidates_list_version(&l, list, list_length, end_pos, 0);
                }
            }
            else
            {
                avalible_k = 0;
                init_Hash_code(&k_code);
            }

            ///HPC_base++;
        }



        ///reverse complement strand
        reverse_complement(g_read.seq, g_read.length);
        init_HPC_seq(&HPC_read, g_read.seq, g_read.length);
        init_Hash_code(&k_code);
        avalible_k = 0;

        while ((code = get_HPC_code(&HPC_read, &end_pos)) != 6)
        {
            if(code < 4)
            {
                k_mer_append(&k_code,code,k_mer_length);
                avalible_k++;
                if (avalible_k>=k_mer_length)
                {
                    list_length = locate_Total_Pos_Table(&PCB, &k_code, &list, k_mer_length, &sub_ID);
                    merge_Candidates_list(&l, list, list_length, end_pos, 1);
                    //merge_Candidates_list_version(&l, list, list_length, end_pos, 1);
                }
            }
            else
            {
                avalible_k = 0;
                init_Hash_code(&k_code);
            }

            ///HPC_base++;
        }

    }
    
    destory_Candidates_list(&l);
}




void debug_merge_result(Candidates_list* x, Candidates_list* y)
{
    uint64_t i;
    if (x->length != y->length)
    {
        fprintf(stderr, "ERROR: different list\n");
        fprintf(stderr, "x->length: %llu, y->length: %llu\n", x->length, y->length);
    }
    
    for (i = 0; i < x->length; i++)
    {
        if (x->list[i].offset != y->list[i].offset 
            ||
            x->list[i].readID != y->list[i].readID
            ||
            x->list[i].self_offset != y->list[i].self_offset
            ||
            x->list[i].strand != y->list[i].strand
        )
        {
            fprintf(stderr, "ERROR: different pos\n");
        }
        
    }
    
}



void* Overlap_calculate(void* arg)
{

    int thr_ID = *((int*)arg);

    long long i = 0;
    int avalible_k = 0;

    UC_Read g_read;
    init_UC_Read(&g_read);
    HPC_seq HPC_read;
    Hash_code k_code;
    uint64_t code;
    uint64_t end_pos;
    k_mer_pos* list;
    uint64_t list_length;
    uint64_t sub_ID;

    Candidates_list l;
    //Candidates_list debug_l;

    init_Candidates_list(&l);
    //init_Candidates_list(&debug_l);

    k_mer_pos_list_alloc array_list;
    init_k_mer_pos_list_alloc(&array_list);


    for (i = thr_ID; i < R_INF.total_reads; i = i + thread_num)
    ///for (i = thr_ID; i < R_INF.total_reads/50; i = i + thread_num)
    {
        /**
        if (i % 1000 == 0)
        {
            fprintf(stderr, "i: %llu\n", i);
        }
        **/
        
        

        clear_Candidates_list(&l);
        ///clear_Candidates_list(&debug_l);

        clear_k_mer_pos_list_alloc(&array_list);
        recover_UC_Read(&g_read, &R_INF, i);


        ///forward strand
        init_HPC_seq(&HPC_read, g_read.seq, g_read.length);
        init_Hash_code(&k_code);
        avalible_k = 0;

        while ((code = get_HPC_code(&HPC_read, &end_pos)) != 6)
        {
            if(code < 4)
            {
                k_mer_append(&k_code,code,k_mer_length);
                avalible_k++;
                if (avalible_k>=k_mer_length)
                {
                    list_length = locate_Total_Pos_Table(&PCB, &k_code, &list, k_mer_length, &sub_ID);

                    if (list_length != 0)
                    {
                        append_k_mer_pos_list_alloc(&array_list, list, list_length, end_pos, 0);
                    }
                    ///merge_Candidates_list(&l, list, list_length, end_pos, 0);
                    ///merge_Candidates_list_version(&debug_l, list, list_length, end_pos, 0);
                }
            }
            else
            {
                avalible_k = 0;
                init_Hash_code(&k_code);
            }

            ///HPC_base++;
        }



        ///reverse complement strand
        reverse_complement(g_read.seq, g_read.length);
        init_HPC_seq(&HPC_read, g_read.seq, g_read.length);
        init_Hash_code(&k_code);
        avalible_k = 0;

        while ((code = get_HPC_code(&HPC_read, &end_pos)) != 6)
        {
            if(code < 4)
            {
                k_mer_append(&k_code,code,k_mer_length);
                avalible_k++;
                if (avalible_k>=k_mer_length)
                {
                    list_length = locate_Total_Pos_Table(&PCB, &k_code, &list, k_mer_length, &sub_ID);
                    if (list_length != 0)
                    {
                        append_k_mer_pos_list_alloc(&array_list, list, list_length, end_pos, 1);
                    }
                    ///merge_Candidates_list(&l, list, list_length, end_pos, 1);
                    ///merge_Candidates_list_version(&debug_l, list, list_length, end_pos, 1);
                }
            }
            else
            {
                avalible_k = 0;
                init_Hash_code(&k_code);
            }

            ///HPC_base++;
        }

        merge_k_mer_pos_list_alloc(&array_list, &l);

        ///debug_merge_result(&l, &debug_l);

    }
    
    destory_Candidates_list(&l);
    ///destory_Candidates_list(&debug_l);


    destory_k_mer_pos_list_alloc(&array_list);
}


void* Overlap_calculate_heap_merge(void* arg)
{

    int thr_ID = *((int*)arg);
    uint64_t POA_i;
    long long i = 0;
    int avalible_k = 0;

    UC_Read g_read;
    init_UC_Read(&g_read);
    HPC_seq HPC_read;
    Hash_code k_code;
    uint64_t code;
    uint64_t end_pos;
    k_mer_pos* list;
    uint64_t list_length;
    uint64_t sub_ID;
    long long total_shared_seed = 0;
    long long candidate_overlap_reads = 0;
    

    Candidates_list l;
    //Candidates_list debug_l;
    Graph POA_Graph;

    init_Graph(&POA_Graph);
    init_Candidates_list(&l);
    //init_Candidates_list(&debug_l);

    k_mer_pos_list_alloc array_list;
    init_k_mer_pos_list_alloc(&array_list);

    overlap_region_alloc overlap_list;
    init_overlap_region_alloc(&overlap_list);

    HeapSq heap;

    Init_Heap(&heap);

    Correct_dumy correct;
    init_Correct_dumy(&correct);

    for (i = thr_ID; i < R_INF.total_reads; i = i + thread_num)
    {

        clear_Heap(&heap);
        clear_Candidates_list(&l);
        ///clear_Candidates_list(&debug_l);

        clear_k_mer_pos_list_alloc(&array_list);
        clear_overlap_region_alloc(&overlap_list);

        recover_UC_Read(&g_read, &R_INF, i);


        ///forward strand
        init_HPC_seq(&HPC_read, g_read.seq, g_read.length);
        init_Hash_code(&k_code);
        avalible_k = 0;

        while ((code = get_HPC_code(&HPC_read, &end_pos)) != 6)
        {
            if(code < 4)
            {
                k_mer_append(&k_code,code,k_mer_length);
                avalible_k++;
                if (avalible_k>=k_mer_length)
                {
                    list_length = locate_Total_Pos_Table(&PCB, &k_code, &list, k_mer_length, &sub_ID);

                    if (list_length != 0)
                    {
                        append_k_mer_pos_list_alloc(&array_list, list, list_length, end_pos, 0);
                    }
                    ///merge_Candidates_list(&l, list, list_length, end_pos, 0);
                    //merge_Candidates_list_version(&debug_l, list, list_length, end_pos, 0);
                }
            }
            else
            {
                avalible_k = 0;
                init_Hash_code(&k_code);
            }

            ///HPC_base++;
        }




        



        ///reverse complement strand
        reverse_complement(g_read.seq, g_read.length);
        /**
        UC_Read rc_read;
        init_UC_Read(&rc_read);
        recover_UC_Read_RC(&rc_read, &R_INF, i);
        uint64_t j = 0;
        for (j = 0; j < g_read.length; j++)
        {
            if (g_read.seq[j] != rc_read.seq[j])
            {
                fprintf(stderr, "j error: %llu, i: %llu\n", j, i);
                fprintf(stderr, "g_read.seq[j]: %c\n", g_read.seq[j]);
                fprintf(stderr, "rc_read.seq[j]: %c\n", rc_read.seq[j]);
            }
            
        }
        destory_UC_Read(&rc_read);
        **/
        init_HPC_seq(&HPC_read, g_read.seq, g_read.length);
        init_Hash_code(&k_code);
        avalible_k = 0;

        while ((code = get_HPC_code(&HPC_read, &end_pos)) != 6)
        {
            if(code < 4)
            {
                k_mer_append(&k_code,code,k_mer_length);
                avalible_k++;
                if (avalible_k>=k_mer_length)
                {
                    list_length = locate_Total_Pos_Table(&PCB, &k_code, &list, k_mer_length, &sub_ID);
                    if (list_length != 0)
                    {
                        append_k_mer_pos_list_alloc(&array_list, list, list_length, end_pos, 1);
                    }
                    ///merge_Candidates_list(&l, list, list_length, end_pos, 1);
                    //merge_Candidates_list_version(&debug_l, list, list_length, end_pos, 1);
                }
            }
            else
            {
                avalible_k = 0;
                init_Hash_code(&k_code);
            }

            ///HPC_base++;
        }

        ///merge_k_mer_pos_list_alloc(&array_list, &l);
        merge_k_mer_pos_list_alloc_heap_sort(&array_list, &l, &heap);
        /**
        if (array_list.length < 3)
        {
            merge_k_mer_pos_list_alloc(&array_list, &l);
        }
        else
        {
            merge_k_mer_pos_list_alloc_heap_sort_advance(&array_list, &l, &heap);
        }
        **/

       ///以x_pos_e，即结束位置为主元排序
       calculate_overlap_region(&l, &overlap_list, i, g_read.length, &R_INF);


       correct_overlap(&overlap_list, &R_INF, &g_read, &correct);

        
        
        /**
        POA_i = 0;
        fprintf(stderr, "\n\n**************\ni: %u\n", i);

        for (POA_i = 0; POA_i < overlap_list.length; POA_i++)
        {
            fprintf(stderr, "x_id: %u, y_id: %u\n", overlap_list.list[POA_i].x_id, overlap_list.list[POA_i].y_id);
            fprintf(stderr, "x_strand: %u, y_strand: %u\n", overlap_list.list[POA_i].x_pos_strand, overlap_list.list[POA_i].y_pos_strand);
            fprintf(stderr, "x_pos_s: %u\n", overlap_list.list[POA_i].x_pos_s);

        }
        **/
        
       /**
       POA_i = 0;

        for (POA_i = 1; POA_i < overlap_list.length; POA_i++)
        {
            if(overlap_list.list[POA_i].x_pos_s < overlap_list.list[POA_i - 1].x_pos_s)
            {
                fprintf(stderr, "1 sbsbsbsbs\n");
            }
            else if(overlap_list.list[POA_i].x_pos_s == overlap_list.list[POA_i - 1].x_pos_s &&
            overlap_list.list[POA_i].x_pos_e > overlap_list.list[POA_i - 1].x_pos_e)
            {
                 fprintf(stderr, "2 sbsbsbsbs\n");
            }
            

        }
        **/



        

        ///merge_k_mer_pos_list_alloc_heap_sort_advance(&array_list, &l, &heap);
        
        

        ///debug_merge_result(&l, &debug_l);

        /**
        clear_Graph(&POA_Graph);


        Perform_POA(&POA_Graph, &overlap_list, &R_INF, &g_read);
        **/

        ///fprintf(stderr, "i: %u\n", i);


        /**
        POA_i = 0;

        if (overlap_list.length > 0)
        {
            
        }
        

        for (POA_i = 1; POA_i < overlap_list.length; POA_i++)
        {
            if(overlap_list.list[POA_i].x_pos_strand == 1)
            {
                fprintf(stderr, "sbsbsbsbs\n");
            }

        }
        **/
        




    }

    /**
    fprintf(stderr, "candidate_overlap_reads: %llu\n", candidate_overlap_reads);
    fprintf(stderr, "total_shared_seed: %llu\n", total_shared_seed);
    **/
    destory_Candidates_list(&l);
    destory_overlap_region_alloc(&overlap_list);
    //destory_Candidates_list(&debug_l);

    destory_Heap(&heap);
    destory_k_mer_pos_list_alloc(&array_list);

    destory_Graph(&POA_Graph);
    destory_UC_Read(&g_read);

    destory_Correct_dumy(&correct);
}








void Overlap_calculate_multipe_thr()
{
    double start_time = Get_T();

    fprintf(stdout, "R_INF.total_reads: %llu\n", R_INF.total_reads);


    fprintf(stdout, "Begin Overlap Calculate ...... \n");

    pthread_t *_r_threads;

	_r_threads = (pthread_t *)malloc(sizeof(pthread_t)*thread_num);

    int i = 0;

	for (i = 0; i < thread_num; i++)
	{
		int *arg = (int*)malloc(sizeof(*arg));
		*arg = i;

		pthread_create(_r_threads + i, NULL, Overlap_calculate_heap_merge, (void*)arg);
        //pthread_create(_r_threads + i, NULL, Overlap_calculate, (void*)arg);

	}
    

    for (i = 0; i<thread_num; i++)
		pthread_join(_r_threads[i], NULL);

    free(_r_threads);

    

    fprintf(stdout, "Finish Overlap Calculate.\n");

    fprintf(stdout, "%-30s%18.2f\n\n", "Calculate Overlap time:", Get_T() - start_time);    
}


int load_pre_cauculated_index()
{
    if(load_Total_Pos_Table(&PCB, read_file_name) && load_All_reads(&R_INF, read_file_name))
    {
        return 1;
    }
    else
    {
        return 0;
    }
    
}


/********************************for debug***************************************/
void Verify_Counting()
{


    kseq_t *seq = (kseq_t*)calloc(1, sizeof(kseq_t));	


    long long read_number = 0;

    HPC_seq HPC_read;
    uint64_t code;
    uint64_t end_pos;
    Hash_code k_code;
    int avalible_k = 0;

    fprintf(stdout, "Start Verifying ...\n");

    while (get_read(seq))
    {

        init_HPC_seq(&HPC_read, seq->seq.s, seq->seq.l);
        init_Hash_code(&k_code);

        avalible_k = 0;

        while ((code = get_HPC_code(&HPC_read, &end_pos)) != 6)
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
    ///long long HPC_base;

    HPC_seq HPC_read;
    uint64_t code;
    uint64_t end_pos;
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
        ///HPC_base = 0;

        while ((code = get_HPC_code(&HPC_read, &end_pos)) != 6)
        {


            if(code < 4)
            {
                k_mer_append(&k_code,code,k_mer_length);
                avalible_k++;
                if (avalible_k>=k_mer_length)
                {

                    ///uint64_t count1 = get_Total_Count_Table(&TCB, &k_code, k_mer_length);
                    uint64_t count2 = count_Total_Pos_Table(&PCB, &k_code, k_mer_length);
                    /**
                    if(count1>=k_mer_min_freq && count1<= k_mer_max_freq && count1!=count2)
                    {
                        fprintf(stderr, "count1: %lld\n",count1);
                        fprintf(stderr, "count2: %lld\n",count2);
                    }
                    **/



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
                                if(list[j].readID == read_number && list[j].offset == end_pos)
                                {
                                    break;
                                }
                            }

                            if(j == count2)
                            {
                                fprintf(stderr, "locate error, read_number: %llu, pos: %llu\n", 
                                read_number, end_pos);

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

            ///HPC_base++;
            
        }

        read_number++;
    }

    destory_kseq();

    fprintf(stdout, "Finish Verifying Position Table!\n");

}












