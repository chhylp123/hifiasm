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
#include "Output.h"

Total_Count_Table TCB;
Total_Pos_Table PCB;
All_reads R_INF;
Assembly_Graph assembly;

pthread_mutex_t statistics;
long long total_matched_overlap_0 = 0;
long long total_matched_overlap_1 = 0;
long long total_potiental_matched_overlap_0 = 0;
long long total_potiental_matched_overlap_1 = 0;
long long total_num_read_base = 0;
long long total_num_correct_base = 0;
long long total_second_num_correct_base = 0;
int roundID = 0;


long long complete_threads = 0;


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
    ///free(arg);

}



void* Perform_Counting_non_first(void* arg)
{
    int thr_ID = *((int*)arg);

    int i = 0;
    HPC_seq HPC_read;


    long long read_number = 0;
    long long select_k_mer_number = 0 ;
    long long k_mer_number = 0 ;

	int file_flag = 1;

    uint64_t code;

    uint64_t end_pos;

    Hash_code k_code;

    int avalible_k = 0;

    UC_Read g_read;
    init_UC_Read(&g_read);

    for (i = thr_ID; i < R_INF.total_reads; i = i + thread_num)
    {
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

    destory_UC_Read(&g_read);
}



void* Build_hash_table_non_first(void* arg)
{
    int thr_ID = *((int*)arg);

    int i = 0;
    HPC_seq HPC_read;

	int file_flag = 1;

    uint64_t code;
    uint64_t end_pos;

    ///long long HPC_base;

    Hash_code k_code;

    int avalible_k = 0;

    UC_Read g_read;
    init_UC_Read(&g_read);

    for (i = thr_ID; i < R_INF.total_reads; i = i + thread_num)
    {
        recover_UC_Read(&g_read, &R_INF, i);
        ///forward strand
        init_HPC_seq(&HPC_read, g_read.seq, g_read.length);
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
                        i, end_pos);
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

    destory_UC_Read(&g_read);
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
    ///free(arg);



}









void Counting_multiple_thr()
{

    double start_time = Get_T();

    fprintf(stdout, "Begin Counting ...... \n");

    init_Total_Count_Table(k_mer_length, &TCB);

    pthread_t inputReadsHandle;

    int *is_insert = (int*)malloc(sizeof(*is_insert));
	*is_insert = 1;

    if (roundID == 0)
    {
        init_kseq(read_file_name);
        init_All_reads(&R_INF);
        init_R_buffer(thread_num);
        pthread_create(&inputReadsHandle, NULL, input_reads_muti_threads, (void*)is_insert);
    }
    

    pthread_t *_r_threads;

	_r_threads = (pthread_t *)malloc(sizeof(pthread_t)*thread_num);

    int i = 0;

	for (i = 0; i < thread_num; i++)
	{
		int *arg = (int*)malloc(sizeof(*arg));
		*arg = i;
        if (roundID == 0)
        {
            pthread_create(_r_threads + i, NULL, Perform_Counting, (void*)arg);
        }
        else
        {
            pthread_create(_r_threads + i, NULL, Perform_Counting_non_first, (void*)arg);
        }
	}



    

    for (i = 0; i<thread_num; i++)
		pthread_join(_r_threads[i], NULL);

    free(_r_threads);

    ///destory_R_buffer();
    
    ///destory_Total_Count_Table(&TCB);
    if (roundID == 0)
    {
        pthread_join(inputReadsHandle, NULL);
        destory_kseq();
    }

    fprintf(stdout, "Finish Counting ...... \n");

    fprintf(stdout, "%-30s%18.2f\n\n", "Counting time:", Get_T() - start_time);
   
    free(is_insert);

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
   
    ///at this moment, TCB can be free
    destory_Total_Count_Table(&TCB);

    pthread_t inputReadsHandle;

    int *is_insert = (int*)malloc(sizeof(*is_insert));
	*is_insert = 0;

    if (roundID == 0)
    {
        init_kseq(read_file_name);

        clear_R_buffer();

        malloc_All_reads(&R_INF);

        pthread_create(&inputReadsHandle, NULL, input_reads_muti_threads, (void*)is_insert);
    }



    pthread_t *_r_threads;

	_r_threads = (pthread_t *)malloc(sizeof(pthread_t)*thread_num);

    int i = 0;

	for (i = 0; i < thread_num; i++)
	{
		int *arg = (int*)malloc(sizeof(*arg));
		*arg = i;

        if (roundID == 0)
        {
            pthread_create(_r_threads + i, NULL, Build_hash_table, (void*)arg);
        }
        else
        {
            pthread_create(_r_threads + i, NULL, Build_hash_table_non_first, (void*)arg);
        }
	}

    if (roundID == 0)
    {
        pthread_join(inputReadsHandle, NULL);
    }

    for (i = 0; i<thread_num; i++)
		pthread_join(_r_threads[i], NULL);

    free(_r_threads);

    

    fprintf(stdout, "Finish Building hash table ...... \n");

    fprintf(stdout, "%-30s%18.2f\n\n", "Build hash table time:", Get_T() - start_time);

    if (roundID == 0)
    {
        destory_kseq();
        destory_R_buffer();
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

    ///destory_Total_Count_Table(&TCB);
    
    free(is_insert);
    
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

inline void output_read_to_buffer(long long readID, All_reads* R_INF, char* corrected_read, long long correct_read_length,
Output_buffer_sub_block* current_sub_buffer)
{
    /**
    if (seq->name.l
            != Get_NAME_LENGTH(R_INF, read_number))
        {
            fprintf(stderr, "name error\n");
        }


         if(memcmp(seq->name.s, Get_NAME(R_INF, read_number), seq->name.l))
        {
            fprintf(stderr, "name error\n");
        }
        **/


    ///先清空
    current_sub_buffer->length = 0;
    
    ///一个是>一个是\n
    add_base_to_sub_buffer(current_sub_buffer, '>');
    add_segment_to_sub_buffer(current_sub_buffer, Get_NAME((*R_INF), readID), Get_NAME_LENGTH((*R_INF), readID));
    add_base_to_sub_buffer(current_sub_buffer, '\n');
    add_segment_to_sub_buffer(current_sub_buffer, corrected_read, correct_read_length);
    add_base_to_sub_buffer(current_sub_buffer, '\n');
    push_results_to_buffer(current_sub_buffer);
    
}


void get_corrected_read_from_cigar(Cigar_record* cigar, char* pre_read, int pre_length,
char* new_read, int* new_length)
{
    int i, j;
    int pre_i, new_i;
    int operation, operation_length;
    pre_i = new_i = 0;
    int diff_char_i = 0;


    for (i = 0; i < cigar->length; i++)
    {
        operation = Get_Cigar_Type(cigar->record[i]);
        operation_length = Get_Cigar_Length(cigar->record[i]);

        if (operation == 0)
        {
            ///fprintf(stderr, "0 new_i: %d\n", new_i);
            memcpy(new_read + new_i, pre_read + pre_i, operation_length);
            pre_i = pre_i + operation_length;
            new_i = new_i + operation_length;
            ///fprintf(stderr, "0 new_i: %d\n", new_i);
        }
        else if (operation == 1)
        {
            
            for (j = 0; j < operation_length; j++)
            {
                // fprintf(stderr, "1 new_i: %d, diff_char_i: %d, lost_base_length: %d, lost_base: %d\n", 
                // new_i, diff_char_i, cigar->lost_base_length, cigar->lost_base[diff_char_i]);

                new_read[new_i] = Get_MisMatch_Base(cigar->lost_base[diff_char_i]);
                new_i++;
                diff_char_i++;
                ///fprintf(stderr, "1 new_i: %d\n", new_i);
            }
            pre_i = pre_i + operation_length;
        }
        else if (operation == 3)
        {
            pre_i = pre_i + operation_length;
            diff_char_i = diff_char_i + operation_length;
        }
        else if (operation == 2)
        {
            ///fprintf(stderr, "2 new_i: %d\n", new_i);
            memcpy(new_read + new_i, cigar->lost_base + diff_char_i, operation_length);
            new_i = new_i + operation_length;
            diff_char_i = diff_char_i + operation_length;
            ///fprintf(stderr, "2 new_i: %d\n", new_i);
        }
    }
    *new_length = new_i;

    // if(pre_i != pre_length)
    // {
    //     fprintf(stderr, "error\n");
    // }
    ///0xffffffff;
}


void get_uncorrected_read_from_cigar(Cigar_record* cigar, char* new_read, int new_length, char* pre_read, int* pre_length)
{
    int i, j;
    int pre_i, new_i;
    int operation, operation_length;
    pre_i = new_i = 0;
    int diff_char_i = 0;


    for (i = 0; i < cigar->length; i++)
    {
        operation = Get_Cigar_Type(cigar->record[i]);
        operation_length = Get_Cigar_Length(cigar->record[i]);


        if (operation == 0)
        {
            memcpy(pre_read + pre_i, new_read + new_i, operation_length);
            pre_i = pre_i + operation_length;
            new_i = new_i + operation_length;
        }
        else if (operation == 1)
        {
            
            for (j = 0; j < operation_length; j++)
            {
                pre_read[pre_i] = Get_Match_Base(cigar->lost_base[diff_char_i]);
                pre_i++;
                diff_char_i++;
            }
            new_i = new_i + operation_length;
        }
        else if (operation == 3)
        {
            memcpy(pre_read + pre_i, cigar->lost_base + diff_char_i, operation_length);
            pre_i = pre_i + operation_length;
            diff_char_i = diff_char_i + operation_length;
        }
        else if (operation == 2)
        {
            new_i = new_i + operation_length;
            diff_char_i = diff_char_i + operation_length;
        }
    }

    *pre_length = pre_i;
}


inline int get_cigar_errors(Cigar_record* cigar)
{
    int i;
    int total_errors = 0;
    for (i = 0; i < cigar->length; i++)
    {
        if (Get_Cigar_Type(cigar->record[i]) > 0)
        {
            total_errors = total_errors + Get_Cigar_Length(cigar->record[i]);
        }
    }

    return total_errors;
}

int debug_cigar(Cigar_record* cigar, char* pre_read, int pre_length,
char* new_read, int new_length, int correct_base)
{
    int i;
    int total_errors = 0;
    for (i = 0; i < cigar->length; i++)
    {
        if (Get_Cigar_Type(cigar->record[i]) > 0)
        {
            total_errors = total_errors + Get_Cigar_Length(cigar->record[i]);
        }
    }

    if(total_errors!=correct_base)
    {
        fprintf(stderr, "total_errors: %d, correct_base: %d\n", total_errors, correct_base);
    }

    

    int pre_i, new_i;
    int operation, operation_length;
    pre_i = new_i = 0;


    for (i = 0; i < cigar->length; i++)
    {
        operation = Get_Cigar_Type(cigar->record[i]);
        operation_length = Get_Cigar_Length(cigar->record[i]);

        if (operation == 0)
        {
            pre_i = pre_i + operation_length;
            new_i = new_i + operation_length;
        }

        if (operation == 1)
        {
            pre_i = pre_i + operation_length;
            new_i = new_i + operation_length;
        }

        if (operation == 3)
        {
            pre_i = pre_i + operation_length;
        }

        if (operation == 2)
        {
            new_i = new_i + operation_length;
        }
        
    }



    if (pre_i != pre_length)
    {
        fprintf(stderr, "pre_i: %d, pre_length: %d\n", pre_i, pre_length);
    }
    
    
    if(new_i != new_length)
    {
        fprintf(stderr, "new_i: %d, new_length: %d\n", new_i, new_length);
    }

    return 1;
    

    char* tmp_seq = (char*)malloc(new_length + pre_length);
    int tmp_length;

    get_corrected_read_from_cigar(cigar, pre_read, pre_length, tmp_seq, &tmp_length);

    if(tmp_length != new_length)
    {
        fprintf(stderr, "tmp_length: %d, new_length: %d\n", tmp_length, new_length);
    }

    if(memcmp(new_read, tmp_seq, new_length)!=0)
    {
        fprintf(stderr, "error new string\n");
    }




    
    get_uncorrected_read_from_cigar(cigar, new_read, new_length, tmp_seq, &tmp_length);

    
    if(tmp_length != pre_length)
    {
        fprintf(stderr, "tmp_length: %d, pre_length: %d\n", tmp_length, pre_length);
    }

    if(memcmp(pre_read, tmp_seq, pre_length)!=0)
    {
        fprintf(stderr, "error pre string\n");
    }
    

    free(tmp_seq);

    if(cigar->new_read_length != new_length)
    {
        fprintf(stderr, "cigar->new_read_length: %d, new_length: %d\n", cigar->new_read_length, new_length);
    }
    
    
    
}



inline void push_cigar(Compressed_Cigar_record* records, long long ID, Cigar_record* input)
{

    if (input->length > records[ID].size)
    {
        records[ID].size = input->length;
        records[ID].record = (uint32_t*)realloc(records[ID].record, records[ID].size*sizeof(uint32_t));     
    }
    records[ID].length = input->length;
    memcpy(records[ID].record, input->record, input->length*sizeof(uint32_t));

    if (input->lost_base_length > records[ID].lost_base_size)
    {
        records[ID].lost_base_size = input->lost_base_length;
        records[ID].lost_base = (char*)realloc(records[ID].lost_base, records[ID].lost_base_size);
    }
    records[ID].lost_base_length = input->lost_base_length;
    memcpy(records[ID].lost_base, input->lost_base, input->lost_base_length);

    records[ID].new_length = input->new_read_length;
    
    
}


void push_overlaps(ma_hit_t_alloc* paf, overlap_region_alloc* overlap_list, int flag)
{
    long long i = 0;
    ma_hit_t tmp;
    clear_ma_hit_t_alloc(paf);
    for (i = 0; i < overlap_list->length; i++)
    {
        if (overlap_list->list[i].is_match == flag)
        {
            tmp.qns = overlap_list->list[i].x_id;
            tmp.qns = tmp.qns << 32;
            tmp.qns = tmp.qns | (uint64_t)(overlap_list->list[i].x_pos_s);

            tmp.qe = overlap_list->list[i].x_pos_e;

            tmp.tn = overlap_list->list[i].y_id;
            tmp.ts = overlap_list->list[i].y_pos_s;
            tmp.te = overlap_list->list[i].y_pos_e;

            ///for overlap_list, the x_strand of all overlaps are 0, so the tmp.rev is the same as the y_strand
            tmp.rev = overlap_list->list[i].y_pos_strand;

            tmp.bl = R_INF.read_length[overlap_list->list[i].y_id];
            tmp.ml = overlap_list->list[i].strong;
            tmp.no_l_indel = overlap_list->list[i].without_large_indel;

            add_ma_hit_t_alloc(paf, &tmp);
        }
    }
    
}

int check_weak_overlap(ma_hit_t_alloc* reverse_paf_list, 
overlap_region_alloc* overlap_list, long long weakID)
{
    long long i = 0;
    long long strongID, index;
    for (i = 0; i < overlap_list->length; i++)
    {
        ///if this is a matched strong overlap
        if (overlap_list->list[i].is_match == 1 && overlap_list->list[i].strong == 1)
        {
            strongID = overlap_list->list[i].y_id;
            index = get_specific_overlap(&(reverse_paf_list[strongID]), strongID, weakID);
            if(index != -1)
            {
                return 0;
            }
        }
    }

    return 1;
}

int if_exact_match(char* x, long long xLen, char* y, long long yLen, 
long long xBeg, long long xEnd, long long yBeg, long long yEnd)
{
    long long overlapLen = xEnd - xBeg + 1;

    if(yEnd - yBeg + 1 == overlapLen)
    {
        long long i;

        for (i = 0; i < overlapLen; i++)
        {
            if(x[xBeg + i] != y[yBeg + i])
            {
                break;
            }
        }

        if(i == overlapLen)
        {
            return 1;
        }
    }
    
    return 0;
}

long long push_final_overlaps(ma_hit_t_alloc* paf, ma_hit_t_alloc* reverse_paf_list, 
overlap_region_alloc* overlap_list, UC_Read* x_read, UC_Read* y_read)
{
    long long i = 0;
    long long available_overlaps = 0;
    ma_hit_t tmp;
    clear_ma_hit_t_alloc(paf);
    for (i = 0; i < overlap_list->length; i++)
    {
        if (overlap_list->list[i].is_match == 1)
        {
            available_overlaps++;
            /**********************query***************************/
            //the interval of overlap is half-open [start, end) 
            tmp.qns = overlap_list->list[i].x_id;
            tmp.qns = tmp.qns << 32;
            tmp.qns = tmp.qns | (uint64_t)(overlap_list->list[i].x_pos_s);
            ///the end pos is open
            tmp.qe = overlap_list->list[i].x_pos_e + 1;
            /**********************query***************************/



            ///for overlap_list, the x_strand of all overlaps are 0, so the tmp.rev is the same as the y_strand
            tmp.rev = overlap_list->list[i].y_pos_strand;


            /**********************target***************************/
            tmp.tn = overlap_list->list[i].y_id;
            if(tmp.rev == 1)
            {
                long long y_readLen = R_INF.read_length[overlap_list->list[i].y_id];
                tmp.ts = y_readLen - overlap_list->list[i].y_pos_e - 1;
                tmp.te = y_readLen - overlap_list->list[i].y_pos_s - 1;
            }
            else
            {
                tmp.ts = overlap_list->list[i].y_pos_s;
                tmp.te = overlap_list->list[i].y_pos_e;
            }
            ///the end pos is open
            tmp.te++;
            /**********************target***************************/
            
            tmp.bl = R_INF.read_length[overlap_list->list[i].y_id];
            tmp.ml = overlap_list->list[i].strong;
            tmp.no_l_indel = overlap_list->list[i].without_large_indel;





            if(overlap_list->list[i].y_pos_strand == 0)
            {
                recover_UC_Read(y_read, &R_INF, overlap_list->list[i].y_id);
            }
            else
            {
                recover_UC_Read_RC(y_read, &R_INF, overlap_list->list[i].y_id);
            }
            
            tmp.el = if_exact_match(x_read->seq, x_read->length, y_read->seq, y_read->length, 
            overlap_list->list[i].x_pos_s, overlap_list->list[i].x_pos_e, 
            overlap_list->list[i].y_pos_s, overlap_list->list[i].y_pos_e);




            
            add_ma_hit_t_alloc(paf, &tmp);
        }
    }


    return available_overlaps;
    
}

int fix_overlap_region_by_cigar(long long* r_beg, long long* r_end, Cigar_record* cigar, long long new_read_length,
long long pre_read_length, long long y_ID)
{
    long long pre_r_beg = (*r_beg);
    long long pre_r_end = (*r_end);
    (*r_beg) = (*r_end) = -1;
    long long pre_i, new_i, cigar_i;
    pre_i = new_i = cigar_i = 0;
    int operation, operationLen;

    for (cigar_i = 0; cigar_i < cigar->length; cigar_i++)
    {
        operation = Get_Cigar_Type(cigar->record[cigar_i]);
        operationLen = Get_Cigar_Length(cigar->record[cigar_i]);
        if (operation == 0 || operation == 1)
        {
            if(pre_r_beg >= pre_i && pre_r_beg < pre_i + operationLen)
            {
                (*r_beg) = new_i + (pre_r_beg - pre_i);
            }


            if(pre_r_end >= pre_i && pre_r_end < pre_i + operationLen)
            {
                (*r_end) = new_i + (pre_r_end - pre_i);
            }

            if((*r_beg) != -1 && (*r_end) != -1)
            {
                return 1;
            }


            new_i += operationLen;
            pre_i += operationLen;
        }
        else if (operation == 2) ///2是x缺字符（y多字符）
        { 
            new_i += operationLen;
        }///3是y缺字符（x多字符）
        else if (operation == 3)
        {

            if(pre_r_beg >= pre_i && pre_r_beg < pre_i + operationLen)
            {
                (*r_beg) = new_i - 1;
                if((*r_beg) < 0)
                {
                    (*r_beg) = 0;
                }
            }


            if(pre_r_end >= pre_i && pre_r_end < pre_i + operationLen)
            {
                (*r_end) = new_i;
                if((*r_end) >= new_read_length)
                {
                    (*r_end) = new_read_length - 1;
                }
            }

            if((*r_beg) != -1 && (*r_end) != -1)
            {
                return 1;
            }

            pre_i += operationLen;
        }
    }

    if(pre_i != pre_read_length)
    {
                if(memcmp("m64013_190412_043951/172426111/ccs", Get_NAME(R_INF, y_ID), 
        Get_NAME_LENGTH(R_INF, y_ID)) == 0)
        {
            fprintf(stderr, "error, pre_i: %d, new_i:%d, pre_read_length: %d, new_read_length:%d\n", 
            pre_i, new_i, pre_read_length, new_read_length);
            for (cigar_i = 0; cigar_i < cigar->length; cigar_i++)
            {
                operation = Get_Cigar_Type(cigar->record[cigar_i]);
                operationLen = Get_Cigar_Length(cigar->record[cigar_i]);
                fprintf(stderr, "operation: %d, operationLen: %d\n", operation, operationLen);
            }
            
        }
        ///fprintf(stderr, "error, pre_i: %d, pre_read_length: %d\n", pre_i, pre_read_length);
        return 0;
    }

    ///return;

    if((*r_end) == -1 && pre_r_end >= new_read_length)
    {
        (*r_end) = new_read_length - 1;
    }

    if((*r_beg) == -1 && pre_r_beg >= new_read_length)
    {
        (*r_beg) = new_read_length - 1;
    }

    return 1;

    // if((*r_beg) == -1 || (*r_end) == -1)
    // {
    //     fprintf(stderr, "pre_r_beg: %d, pre_r_end: %d, pre_read_length: %d, new_read_length: %d\n",
    //     pre_r_beg, pre_r_end, pre_read_length, new_read_length);
    //     for (cigar_i = 0; cigar_i < cigar->length; cigar_i++)
    //     {
    //         operation = Get_Cigar_Type(cigar->record[cigar_i]);
    //         operationLen = Get_Cigar_Length(cigar->record[cigar_i]);
    //         fprintf(stderr, "operation: %d, operationLen: %d\n", operation, operationLen);
    //     }
    // }
}

void convert_kmer(k_v* kv, Hash_code* k_code, uint64_t end_pos)
{
    kv->key.x[0] = k_code->x[0];
    kv->key.x[1] = k_code->x[1];
    kv->value = end_pos;
}


void debug_sort_small_hash(small_hash_table* Table)
{
    long long i;
    for (i = 1; i < Table->length; i++)
    {
        if(compare_k_mer(&Table->buffer[i], &Table->buffer[i-1]) < 0)
        {
            fprintf(stderr, "error\n");
        }

        if(compare_k_mer(&Table->buffer[i], &Table->buffer[i-1]) == 0)
        {
            if(Table->buffer[i].value < Table->buffer[i-1].value)
            {
                fprintf(stderr, "error\n");
            }
        }
    }

    for (i = 0; i < Table->length; i++)
    {
        long long left, right;

        if(query_small_hash_table(Table, &Table->buffer[i], &left, &right) == 0)
        {
            fprintf(stderr, "Table->length: %d, i: %d, x[1]: %llu, x[0]: %llu, value: %llu\n", 
            Table->length, i,
            Table->buffer[i].key.x[1], Table->buffer[i].key.x[0],
            Table->buffer[i].value);
        }

        long long single_count = 0;
        long long first_i = -1;

        for (long long j = 0; j < Table->length; j++)
        {
            if(compare_k_mer(&Table->buffer[i], &Table->buffer[j]) == 0)
            {
                single_count++;
                if(first_i == -1)
                {
                    first_i = j;
                }
            }
            else if(compare_k_mer(&Table->buffer[i], &Table->buffer[j]) < 0)
            {
                break;
            }
            
        }

        if(left != first_i)
        {
            fprintf(stderr, "error left\n");
        }

        if(right - left + 1 != single_count)
        {
            fprintf(stderr, "error right\n");
        }
        
    }
}

void get_candidates_from_existing_overlaps(long long readID, UC_Read* g_read, UC_Read* overlap_read, 
overlap_region_alloc* overlap_list, k_mer_pos_list_alloc* array_list,
HeapSq* heap, Candidates_list* l, small_hash_table* forward, small_hash_table* reverse)
{
    HPC_seq HPC_read;
    Hash_code k_code;
    long long avalible_k;
    uint64_t code;
    uint64_t end_pos;
    k_mer_pos* list;
    uint64_t list_length;
    uint64_t sub_ID;
    k_v k_mer_kv;

    clear_Heap(heap);
    clear_Candidates_list(l);

    clear_k_mer_pos_list_alloc(array_list);
    clear_overlap_region_alloc(overlap_list);

    clear_small_hash_table(forward);
    clear_small_hash_table(reverse);

    recover_UC_Read(g_read, &R_INF, readID);

    ///forward strand
    init_HPC_seq(&HPC_read, g_read->seq, g_read->length);
    init_Hash_code(&k_code);
    avalible_k = 0;

    while ((code = get_HPC_code(&HPC_read, &end_pos)) != 6)
    {
        if(code < 4)
        {
            k_mer_append(&k_code,code, k_mer_length);
            avalible_k++;
            if (avalible_k>= k_mer_length)
            {
                if(if_k_mer_available(&k_code, k_mer_length))
                {
                    convert_kmer(&k_mer_kv, &k_code, end_pos);
                    add_small_hash_table(forward, &k_mer_kv);
                }
            }
        }
        else
        {
            avalible_k = 0;
            init_Hash_code(&k_code);
        }
    }

    sort_small_hash_table(forward);

    ///debug_sort_small_hash(forward);



    /**
    ///reverse complement strand
    reverse_complement(g_read->seq, g_read->length);
    init_HPC_seq(&HPC_read, g_read->seq, g_read->length);
    init_Hash_code(&k_code);
    avalible_k = 0;

    while ((code = get_HPC_code(&HPC_read, &end_pos)) != 6)
    {
        if(code < 4)
        {
            k_mer_append(&k_code,code, k_mer_length);
            avalible_k++;
            if (avalible_k>= k_mer_length)
            {
                if(if_k_mer_available(&k_code, k_mer_length))
                {
                    convert_kmer(&k_mer_kv, &k_code, end_pos);
                    add_small_hash_table(reverse, &k_mer_kv);
                }
            }
        }
        else
        {
            avalible_k = 0;
            init_Hash_code(&k_code);
        }
    }

    sort_small_hash_table(reverse);

    ///debug_sort_small_hash(reverse);
    **/

    long long i;
    long long y_start, y_end, y_strand, y_ID, y_overlapLen, extraLen, y_readLen, tmp;
    char* new_y_string;
    long long new_y_length;
    long long result_left, result_right, result_occ;
    long long total_result_occ = 0;
    Cigar_record cigar;
    for (i = 0; i < R_INF.paf[readID].length; i++)
    {
        y_ID = R_INF.paf[readID].buffer[i].tn;
        y_start = R_INF.paf[readID].buffer[i].ts;
        y_end = R_INF.paf[readID].buffer[i].te;
        y_strand = R_INF.paf[readID].buffer[i].rev;
        y_readLen = R_INF.read_length[y_ID];



        y_overlapLen = y_end - y_start + 1;
        extraLen = y_overlapLen * 0.1;
        y_start = y_start - extraLen;
        if(y_start < 0)
        {
            y_start = 0;
        }
        y_end = y_end + extraLen;
        if(extraLen >= y_readLen)
        {
            y_end = y_readLen - 1;
        }


        y_start = y_readLen - y_start - 1;
        y_end = y_readLen - y_end - 1;
        tmp = y_start;
        y_start = y_end;
        y_end = tmp;


        

        if(y_strand == 0)
        {
            recover_UC_Read(overlap_read, &R_INF, y_ID);
        }
        else
        {
            recover_UC_Read_RC(overlap_read, &R_INF, y_ID);
        }
        
        new_y_string = overlap_read->seq + y_start;
        new_y_length = y_end - y_start + 1;
        
        /**
        new_y_string = overlap_read->seq;
        new_y_length = overlap_read->length;
        **/
        total_result_occ = 0;
        clear_Candidates_list(l);


        init_HPC_seq(&HPC_read, new_y_string, new_y_length);
        init_Hash_code(&k_code);
        avalible_k = 0;
        while ((code = get_HPC_code(&HPC_read, &end_pos)) != 6)
        {
            if(code < 4)
            {
                k_mer_append(&k_code,code, k_mer_length);
                avalible_k++;
                if (avalible_k>= k_mer_length)
                {
                    if(if_k_mer_available(&k_code, k_mer_length))
                    {
                        convert_kmer(&k_mer_kv, &k_code, end_pos);
                        result_occ = query_small_hash_table(forward, &k_mer_kv, &result_left, &result_right);
                        if(result_occ > 0)
                        {
                            total_result_occ += result_occ;
                            
                            insert_kv_list_to_candidates(forward->buffer + result_left, result_occ, y_ID, 
                            end_pos + y_start, y_strand, l);
                            
                           /**
                           insert_kv_list_to_candidates(forward->buffer + result_left, result_occ, y_ID, 
                           end_pos, y_strand, l);
                           **/                            
                        }

                    }
                }
            }
            else
            {
                avalible_k = 0;
                init_Hash_code(&k_code);
            }
        }

        ///要在这里整理出一个candidate位置，然后插入进overlap_list
        sort_candidates(l, readID, overlap_list, &R_INF);

        // fprintf(stderr, "x_pos_s: (pre: %d), x_pos_e: (pre: %d), y_pos_s: (pre: %d), y_pos_e: (pre: %d), y_strand: (pre: %d)\n\n", 
        // (uint32_t)(R_INF.paf[readID].buffer[i].qns), R_INF.paf[readID].buffer[i].qe, R_INF.paf[readID].buffer[i].ts, R_INF.paf[readID].buffer[i].te,
        // R_INF.paf[readID].buffer[i].rev);

        // fprintf(stderr, "x_id: %d, y_ID:%d, x_len: %d, y_len: %d, overlap_list->length: %d\n", readID, y_ID, 
        // R_INF.read_length[readID], R_INF.read_length[y_ID], overlap_list->length);
        // fprintf(stderr, "x_pos_s: %d (pre: %d), x_pos_e: %d (pre: %d), y_pos_s: %d (pre: %d), y_pos_e: %d (pre: %d), x_strand: %d, y_strand: %d (pre: %d)\n\n\n", 
        // overlap_list->list[overlap_list->length - 1].x_pos_s, (uint32_t)(R_INF.paf[readID].buffer[i].qns),
        // overlap_list->list[overlap_list->length - 1].x_pos_e, R_INF.paf[readID].buffer[i].qe,
        // overlap_list->list[overlap_list->length - 1].y_pos_s, R_INF.paf[readID].buffer[i].ts,
        // overlap_list->list[overlap_list->length - 1].y_pos_e, R_INF.paf[readID].buffer[i].te,
        // overlap_list->list[overlap_list->length - 1].x_pos_strand, overlap_list->list[overlap_list->length - 1].y_pos_strand,
        // R_INF.paf[readID].buffer[i].rev);

        
        
        // if(l->length != total_result_occ)
        // {
        //     fprintf(stderr, "error\n");
        // }
        ///fprintf(stderr, "i:%d, total_result_occ: %d, y_strand: %d, new_y_length: %d\n", i, total_result_occ, y_strand, new_y_length);
    }


    qsort(overlap_list->list, overlap_list->length, sizeof(overlap_region), cmp_by_x_pos_s);

    reverse_complement(g_read->seq, g_read->length);


    // for (long long i = 0; i < overlap_list->length; i++)
    // {
    //     fprintf(stderr, "****x_pos_s: %d, x_pos_e: %d, y_pos_s: %d, y_pos_e: %d, x_strand: %d, y_strand: %d, x_id: %d, y_id: %d\n", 
    //         overlap_list->list[i].x_pos_s, overlap_list->list[i].x_pos_e, 
    //         overlap_list->list[i].y_pos_s, overlap_list->list[i].y_pos_e,
    //         overlap_list->list[i].x_pos_strand, overlap_list->list[i].y_pos_strand,
    //         overlap_list->list[i].x_id, overlap_list->list[i].y_id);
    // }
    
}

void get_new_candidates(long long readID, UC_Read* g_read, overlap_region_alloc* overlap_list, k_mer_pos_list_alloc* array_list,
HeapSq* heap, Candidates_list* l, double band_width_threshold)
{
    HPC_seq HPC_read;
    Hash_code k_code;
    long long avalible_k;
    uint64_t code;
    uint64_t end_pos;
    k_mer_pos* list;
    uint64_t list_length;
    uint64_t sub_ID;

    clear_Heap(heap);
    clear_Candidates_list(l);

    clear_k_mer_pos_list_alloc(array_list);
    clear_overlap_region_alloc(overlap_list);

    recover_UC_Read(g_read, &R_INF, readID);

    ///forward strand
    init_HPC_seq(&HPC_read, g_read->seq, g_read->length);
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
                    append_k_mer_pos_list_alloc(array_list, list, list_length, end_pos, 0);
                }
            }
        }
        else
        {
            avalible_k = 0;
            init_Hash_code(&k_code);
        }
    }



    ///reverse complement strand
    reverse_complement(g_read->seq, g_read->length);
    init_HPC_seq(&HPC_read, g_read->seq, g_read->length);
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
                    append_k_mer_pos_list_alloc(array_list, list, list_length, end_pos, 1);
                }
            }
        }
        else
        {
            avalible_k = 0;
            init_Hash_code(&k_code);
        }
    }

    

    merge_k_mer_pos_list_alloc_heap_sort(array_list, l, heap);
    
    ///以x_pos_e，即结束位置为主元排序
    ///calculate_overlap_region(l, overlap_list, readID, g_read->length, &R_INF);
    calculate_overlap_region_by_chaining(l, overlap_list, readID, g_read->length, &R_INF, band_width_threshold);
}




void* Overlap_calculate_heap_merge(void* arg)
{
    /************需要注释掉**********/
    // long long debug_overlap = 0;
    // long long filtered_debug_overlap = 0;
    /************需要注释掉**********/

    long long matched_overlap_0 = 0;
    long long matched_overlap_1 = 0;
    long long potiental_matched_overlap_0 = 0;
    long long potiental_matched_overlap_1 = 0;
    long long num_read_base = 0;
    long long num_correct_base = 0;
    long long num_second_correct_base = 0;
    long long j;
    int fully_cov, abnormal;

    int thr_ID = *((int*)arg);
    uint64_t POA_i;
    long long i = 0;
    int avalible_k = 0;

    UC_Read g_read;
    init_UC_Read(&g_read);

    UC_Read overlap_read;
    init_UC_Read(&overlap_read);

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
    Graph DAGCon;
    init_Graph(&DAGCon);
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


    Output_buffer_sub_block current_sub_buffer;

    init_buffer_sub_block(&current_sub_buffer);

    Cigar_record current_cigar;
    init_Cigar_record(&current_cigar);

    haplotype_evdience_alloc hap;
    InitHaplotypeEvdience(&hap);


    Round2_alignment second_round;
    init_Round2_alignment(&second_round);

    small_hash_table forward, reverse;
    init_small_hash_table(&forward);
    init_small_hash_table(&reverse);


    for (i = thr_ID; i < R_INF.total_reads; i = i + thread_num)
    {
        get_new_candidates(i, &g_read, &overlap_list, &array_list, &heap, &l, 0.02);
        
        clear_Cigar_record(&current_cigar);
        clear_Round2_alignment(&second_round);

        correct_overlap(&overlap_list, &R_INF, &g_read, &correct, &overlap_read, &POA_Graph, &DAGCon,
        &matched_overlap_0, &matched_overlap_1, &potiental_matched_overlap_0, &potiental_matched_overlap_1,
        &current_cigar, &hap, &second_round, 0, 1, &fully_cov, &abnormal);

        num_read_base = num_read_base + g_read.length;
        num_correct_base = num_correct_base + correct.corrected_base;
        num_second_correct_base = num_second_correct_base + second_round.dumy.corrected_base;

        push_cigar(R_INF.cigars, i, &current_cigar);
        push_cigar(R_INF.second_round_cigar, i, &(second_round.cigar));



        R_INF.paf[i].is_fully_corrected = 0;
        if(fully_cov)
        {
           if(
            get_cigar_errors(&current_cigar) == 0 
            && 
            get_cigar_errors(&second_round.cigar) == 0)
            {
                R_INF.paf[i].is_fully_corrected = 1;
            }     
        }
        R_INF.paf[i].is_abnormal = abnormal;





        push_overlaps(&(R_INF.paf[i]), &overlap_list, 1);
        push_overlaps(&(R_INF.reverse_paf[i]), &overlap_list, 2);
    }


    /************需要注释掉**********/
    // fprintf(stderr, "debug_overlap: %llu\n", debug_overlap);
    // fprintf(stderr, "filtered_debug_overlap: %llu\n", filtered_debug_overlap);
    /************需要注释掉**********/

    finish_output_buffer();

    destory_buffer_sub_block(&current_sub_buffer);
    
        

    /**
    fprintf(stderr, "candidate_overlap_reads: %llu\n", candidate_overlap_reads);
    fprintf(stderr, "total_shared_seed: %llu\n", total_shared_seed);
    **/
    destory_Candidates_list(&l);
    destory_overlap_region_alloc(&overlap_list);
    //destory_Candidates_list(&debug_l);

    destory_Heap(&heap);
    destory_k_mer_pos_list_alloc(&array_list);
    ///destory_k_mer_pos_list_alloc_prefilter(&array_list);

    destory_Graph(&POA_Graph);
    destory_Graph(&DAGCon);

    
    destory_UC_Read(&g_read);
    destory_UC_Read(&overlap_read);
    destory_Cigar_record(&current_cigar);

    destory_Correct_dumy(&correct);

    destoryHaplotypeEvdience(&hap);

    destory_Round2_alignment(&second_round);

    destory_small_hash_table(&forward);
    destory_small_hash_table(&reverse);



    pthread_mutex_lock(&statistics);
    total_matched_overlap_0 += matched_overlap_0;
    total_matched_overlap_1 += matched_overlap_1;
    total_potiental_matched_overlap_0 += potiental_matched_overlap_0;
    total_potiental_matched_overlap_1 += potiental_matched_overlap_1;
    total_num_read_base += num_read_base;
    total_num_correct_base += num_correct_base;
    total_second_num_correct_base +=num_second_correct_base;

    complete_threads++;
    if(complete_threads == thread_num)
    {
        fprintf(stderr, "total_matched_overlap_0: %llu\n", total_matched_overlap_0);
        fprintf(stderr, "total_matched_overlap_1: %llu\n", total_matched_overlap_1);
        fprintf(stderr, "total_potiental_matched_overlap_0: %llu\n", total_potiental_matched_overlap_0);
        fprintf(stderr, "total_potiental_matched_overlap_1: %llu\n", total_potiental_matched_overlap_1);
        fprintf(stderr, "total_num_read_base: %llu\n", total_num_read_base);
        fprintf(stderr, "total_num_correct_base: %llu\n", total_num_correct_base);
        fprintf(stderr, "total_second_num_correct_base: %llu\n", total_second_num_correct_base);
        
    }
	pthread_mutex_unlock(&statistics);

    free(arg);
}



void* Output_related_reads(void* arg)
{
    /************需要注释掉**********/
    // long long debug_overlap = 0;
    // long long filtered_debug_overlap = 0;
    /************需要注释掉**********/

    long long matched_overlap_0 = 0;
    long long matched_overlap_1 = 0;
    long long potiental_matched_overlap_0 = 0;
    long long potiental_matched_overlap_1 = 0;
    long long num_read_base = 0;
    long long num_correct_base = 0;
    long long num_second_correct_base = 0;
    long long j;
    int fully_cov;

    int thr_ID = *((int*)arg);
    uint64_t POA_i;
    long long i = 0;
    int avalible_k = 0;

    UC_Read g_read;
    init_UC_Read(&g_read);

    UC_Read overlap_read;
    init_UC_Read(&overlap_read);

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
    Graph DAGCon;
    init_Graph(&DAGCon);
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


    Output_buffer_sub_block current_sub_buffer;

    init_buffer_sub_block(&current_sub_buffer);

    Cigar_record current_cigar;
    init_Cigar_record(&current_cigar);

    haplotype_evdience_alloc hap;
    InitHaplotypeEvdience(&hap);


    Round2_alignment second_round;
    init_Round2_alignment(&second_round);

    small_hash_table forward, reverse;
    init_small_hash_table(&forward);
    init_small_hash_table(&reverse);

    long long required_read_name_length = strlen(required_read_name);
    for (i = thr_ID; i < R_INF.total_reads; i = i + thread_num)
    {
        
        if(required_read_name_length == Get_NAME_LENGTH((R_INF),i)
            &&
           memcmp(required_read_name, Get_NAME((R_INF), i), Get_NAME_LENGTH((R_INF),i)) == 0)
        {
            get_new_candidates(i, &g_read, &overlap_list, &array_list, &heap, &l, 0.02);


            fprintf(stderr, ">%.*s\n", Get_NAME_LENGTH((R_INF), i),
            Get_NAME((R_INF), i));
            recover_UC_Read(&g_read, &R_INF, i);
            fprintf(stderr, "%.*s\n", g_read.length, g_read.seq);


            long long k;
            for (k = 0; k < overlap_list.length; k++)
            {
                fprintf(stderr, ">%.*s\n", Get_NAME_LENGTH((R_INF),overlap_list.list[k].y_id),
                Get_NAME((R_INF),overlap_list.list[k].y_id));
                recover_UC_Read(&g_read, &R_INF, overlap_list.list[k].y_id);
                fprintf(stderr, "%.*s\n", g_read.length, g_read.seq);
            }
            
        }
    }

    finish_output_buffer();

    destory_buffer_sub_block(&current_sub_buffer);
    
        


    destory_Candidates_list(&l);
    destory_overlap_region_alloc(&overlap_list);
    //destory_Candidates_list(&debug_l);

    destory_Heap(&heap);
    destory_k_mer_pos_list_alloc(&array_list);
    ///destory_k_mer_pos_list_alloc_prefilter(&array_list);

    destory_Graph(&POA_Graph);
    destory_Graph(&DAGCon);

    
    destory_UC_Read(&g_read);
    destory_UC_Read(&overlap_read);
    destory_Cigar_record(&current_cigar);

    destory_Correct_dumy(&correct);

    destoryHaplotypeEvdience(&hap);

    destory_Round2_alignment(&second_round);

    destory_small_hash_table(&forward);
    destory_small_hash_table(&reverse);

    free(arg);
}



inline long long get_N_occ(char* seq, long long length)
{
    long long N_occ = 0;
    long long j;
    for (j = 0; j < length; j++)
    {
        if(seq_nt6_table[seq[j]] >= 4)
        {
            ///fprintf(stderr, "seq[%d]: %c\n", j, seq[j]);
            N_occ++;
        }
    }

    return N_occ;
}


void* Save_corrected_reads(void* arg)
{
    int thr_ID = *((int*)arg);
    long long i, j;
    UC_Read g_read;
    init_UC_Read(&g_read);

    int first_round_read_size = 10000;
    char* first_round_read = (char*)malloc(first_round_read_size);

    int second_round_read_size = 10000;
    char* second_round_read = (char*)malloc(second_round_read_size);

    Cigar_record cigar;
    int first_round_read_length;
    int second_round_read_length;
    uint64_t N_occ;

    char* new_read;
    int new_read_length;
    
    for (i = thr_ID; i < R_INF.total_reads; i = i + thread_num)
    {
        recover_UC_Read(&g_read, &R_INF, i);



        /********************************1 round******************************/
        if(R_INF.cigars[i].new_length > first_round_read_size)
        {
            first_round_read_size = R_INF.cigars[i].new_length;
            first_round_read = (char*)realloc(first_round_read, first_round_read_size);
        }

        

        


        cigar.length = R_INF.cigars[i].length;
        cigar.lost_base_length = R_INF.cigars[i].lost_base_length;
        cigar.record = R_INF.cigars[i].record;
        cigar.lost_base = R_INF.cigars[i].lost_base;

        get_corrected_read_from_cigar(&cigar, g_read.seq, g_read.length, first_round_read, &first_round_read_length);
        
        



        /********************debug********************/
        // N_occ = get_N_occ(first_round_read, first_round_read_length);
        // fprintf(stderr, "i: %d, first N_occ: %d\n", i, N_occ);
        // fflush(stderr);
        /********************debug********************/
        /********************************1 round******************************/

        /********************************2 round******************************/
        if(R_INF.second_round_cigar[i].new_length > second_round_read_size)
        {
            second_round_read_size = R_INF.second_round_cigar[i].new_length;
            second_round_read = (char*)realloc(second_round_read, second_round_read_size);
        }
        cigar.length = R_INF.second_round_cigar[i].length;
        cigar.lost_base_length = R_INF.second_round_cigar[i].lost_base_length;
        cigar.record = R_INF.second_round_cigar[i].record;
        cigar.lost_base = R_INF.second_round_cigar[i].lost_base;
        get_corrected_read_from_cigar(&cigar, first_round_read, first_round_read_length, 
        second_round_read, &second_round_read_length);
        

        /********************debug********************/
        // N_occ = get_N_occ(second_round_read, second_round_read_length);
        // fprintf(stderr, "i: %d, second N_occ: %d\n\n\n", i, N_occ);
        // fflush(stderr);
        /********************debug********************/
        /********************************2 round******************************/


        

        new_read = second_round_read;
        new_read_length = second_round_read_length;

        
        if (roundID != number_of_round - 1)
        {
            ///need modification
            reverse_complement(new_read, new_read_length);
        }
        else if(number_of_round % 2 == 0)
        {
            ///need modification
            reverse_complement(new_read, new_read_length);
        }

        
        N_occ = get_N_occ(new_read, new_read_length);

        
        if(R_INF.read_size[i] < new_read_length)
        {
            R_INF.read_size[i] = new_read_length;
            R_INF.read_sperate[i] = (uint8_t*)realloc(R_INF.read_sperate[i], R_INF.read_size[i]/4+1);
        }

        R_INF.read_length[i] = new_read_length;
        

        compress_base(Get_READ(R_INF, i),
            new_read, new_read_length, 
            &R_INF.N_site[i], N_occ);
    }

    destory_UC_Read(&g_read);
    free(first_round_read);
    free(second_round_read);
    free(arg);
}


void Output_corrected_reads()
{
    long long i, j;
    UC_Read g_read;
    init_UC_Read(&g_read);
    FILE* output_file = fopen(output_file_name, "w");


    for (i = 0; i < R_INF.total_reads; i++)
    {
        recover_UC_Read(&g_read, &R_INF, i);
        fwrite(">", 1, 1, output_file);
        fwrite(Get_NAME(R_INF, i), 1, Get_NAME_LENGTH(R_INF, i), output_file);
        fwrite("\n", 1, 1, output_file);
        fwrite(g_read.seq, 1, g_read.length, output_file);
        fwrite("\n", 1, 1, output_file);
    }

    /**
    if(number_of_round % 2 == 0)
    {
        for (i = 0; i < R_INF.total_reads; i++)
        {
            recover_UC_Read(&g_read, &R_INF, i);
            fwrite(">", 1, 1, output_file);
            fwrite(Get_NAME(R_INF, i), 1, Get_NAME_LENGTH(R_INF, i), output_file);
            fwrite("\n", 1, 1, output_file);
            fwrite(g_read.seq, 1, g_read.length, output_file);
            fwrite("\n", 1, 1, output_file);
        }
    }
    else
    {
        for (i = 0; i < R_INF.total_reads; i++)
        {
            recover_UC_Read_RC(&g_read, &R_INF, i);
            fwrite(">", 1, 1, output_file);
            fwrite(Get_NAME(R_INF, i), 1, Get_NAME_LENGTH(R_INF, i), output_file);
            fwrite("\n", 1, 1, output_file);
            fwrite(g_read.seq, 1, g_read.length, output_file);
            fwrite("\n", 1, 1, output_file);
        }
    }
    **/
    
    


    destory_UC_Read(&g_read);
    fclose(output_file);
}




void Overlap_calculate_multipe_thr()
{


    
    double start_time = Get_T();
    pthread_t outputResultSinkHandle;
    
    fprintf(stdout, "Begin Overlap Calculate ...... \n");

    pthread_t *_r_threads;

	_r_threads = (pthread_t *)malloc(sizeof(pthread_t)*thread_num);

    int i = 0;

	for (i = 0; i < thread_num; i++)
	{
		int *arg = (int*)malloc(sizeof(*arg));
		*arg = i;
        if(!required_read_name)
        {
            pthread_create(_r_threads + i, NULL, Overlap_calculate_heap_merge, (void*)arg);
        }
        else
        {
            pthread_create(_r_threads + i, NULL, Output_related_reads, (void*)arg);
        }
        
		
	}
    

    for (i = 0; i<thread_num; i++)
		pthread_join(_r_threads[i], NULL);


    free(_r_threads);

    if(required_read_name)
    {
        exit(0);
    }


    destory_Total_Pos_Table(&PCB);

    fprintf(stdout, "Finish Overlap Calculate.\n");

    fprintf(stdout, "%-30s%18.2f\n\n", "Calculate Overlap time:", Get_T() - start_time);    

    start_time = Get_T();

    _r_threads = (pthread_t *)malloc(sizeof(pthread_t)*thread_num);
    for (i = 0; i < thread_num; i++)
	{
		int *arg = (int*)malloc(sizeof(*arg));
		*arg = i;
		pthread_create(_r_threads + i, NULL, Save_corrected_reads, (void*)arg);
	}


    for (i = 0; i<thread_num; i++)
		pthread_join(_r_threads[i], NULL);
    free(_r_threads);
    fprintf(stdout, "%-30s%18.2f\n\n", "Save corrected read time:", Get_T() - start_time);   
    ///fflush(stdout);


    /*********************debug************************/
    
    ///only the last round can output read to disk
    if (roundID == number_of_round - 1)
    {
        start_time = Get_T();
        
        Output_corrected_reads();
        
        fprintf(stdout, "%-30s%18.2f\n\n", "Output time:", Get_T() - start_time);    
    }
    
   /*********************debug************************/
    

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

    while (get_read(seq, adapterLen))
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

    while (get_read(seq, adapterLen))
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



void get_inexact_overlaps(long long readID, UC_Read* g_read, overlap_region_alloc* overlap_list, k_mer_pos_list_alloc* array_list,
HeapSq* heap, Candidates_list* l)
{
    HPC_seq HPC_read;
    Hash_code k_code;
    long long avalible_k;
    uint64_t code;
    uint64_t end_pos;
    k_mer_pos* list;
    uint64_t list_length;
    uint64_t sub_ID;

    clear_Heap(heap);
    clear_Candidates_list(l);

    clear_k_mer_pos_list_alloc(array_list);
    clear_overlap_region_alloc(overlap_list);

    recover_UC_Read(g_read, &R_INF, readID);

    ///forward strand
    init_HPC_seq(&HPC_read, g_read->seq, g_read->length);
    init_Hash_code(&k_code);
    avalible_k = 0;

    // if(g_read->length != R_INF.read_length[readID])
    // {
    //     fprintf(stderr, "error\n");
    // }

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
                    append_k_mer_pos_list_alloc(array_list, list, list_length, end_pos, 0);
                }
            }
        }
        else
        {
            avalible_k = 0;
            init_Hash_code(&k_code);
        }
    }



    ///reverse complement strand
    reverse_complement(g_read->seq, g_read->length);
    init_HPC_seq(&HPC_read, g_read->seq, g_read->length);
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
                    append_k_mer_pos_list_alloc(array_list, list, list_length, end_pos, 1);
                }
            }
        }
        else
        {
            avalible_k = 0;
            init_Hash_code(&k_code);
        }
    }

    

    merge_k_mer_pos_list_alloc_heap_sort(array_list, l, heap);
    calculate_inexact_overlap_region(l, overlap_list, readID, g_read->length, &R_INF);
    ///calculate_overlap_region(l, overlap_list, readID, g_read->length, &R_INF);




    overlap_region_sort_y_id(overlap_list->list, overlap_list->length);
    ma_hit_sort_tn(R_INF.paf[readID].buffer, R_INF.paf[readID].length);
    overlap_list->mapped_overlaps_length = 0;
    long long inner_j = 0;
    long long j = 0;
    while (j < overlap_list->length && inner_j < R_INF.paf[readID].length)
    {
        if(overlap_list->list[j].y_id < R_INF.paf[readID].buffer[inner_j].tn)
        {
            j++;
        }
        else if(overlap_list->list[j].y_id > R_INF.paf[readID].buffer[inner_j].tn)
        {
            inner_j++;
        }
        else
        {
            if(overlap_list->list[j].y_pos_strand == R_INF.paf[readID].buffer[inner_j].rev)
            {
                overlap_list->list[j].is_match = 1;
                overlap_list->mapped_overlaps_length++;
            }
            j++;
            inner_j++;
        }
    }
}


void* Final_overlap_calculate_heap_merge(void* arg)
{
    long long matched_overlap_0 = 0;
    long long matched_overlap_1 = 0;
    long long potiental_matched_overlap_0 = 0;
    long long potiental_matched_overlap_1 = 0;
    long long num_correct_base = 0;
    long long num_read_base = 0;
    long long num_second_correct_base = 0;
    long long j, inner_j;

    int thr_ID = *((int*)arg);
    uint64_t POA_i;
    long long i = 0;
    int avalible_k = 0;

    UC_Read g_read;
    init_UC_Read(&g_read);

    UC_Read overlap_read;
    init_UC_Read(&overlap_read);

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
    Graph DAGCon;
    init_Graph(&DAGCon);
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


    Output_buffer_sub_block current_sub_buffer;

    init_buffer_sub_block(&current_sub_buffer);

    Cigar_record current_cigar;
    init_Cigar_record(&current_cigar);

    haplotype_evdience_alloc hap;
    InitHaplotypeEvdience(&hap);


    Round2_alignment second_round;
    init_Round2_alignment(&second_round);

    small_hash_table forward, reverse;
    init_small_hash_table(&forward);
    init_small_hash_table(&reverse);


    for (i = thr_ID; i < R_INF.total_reads; i = i + thread_num)
    {
        ///0.1%
        get_new_candidates(i, &g_read, &overlap_list, &array_list, &heap, &l, 0.001);
        /**
        correct_overlap(&overlap_list, &R_INF, &g_read, &correct, &overlap_read, &POA_Graph, &DAGCon,
        &matched_overlap_0, &matched_overlap_1, &potiental_matched_overlap_0, &potiental_matched_overlap_1,
        &current_cigar, &hap, &second_round, 0, 0);

        push_final_overlaps(&(R_INF.paf[i]), &overlap_list);
        **/

        overlap_region_sort_y_id(overlap_list.list, overlap_list.length);
        ma_hit_sort_tn(R_INF.paf[i].buffer, R_INF.paf[i].length);




        overlap_list.mapped_overlaps_length = 0;
        inner_j = 0;
        j = 0;
        while (j < overlap_list.length && inner_j < R_INF.paf[i].length)
        {
            if(overlap_list.list[j].y_id < R_INF.paf[i].buffer[inner_j].tn)
            {
                j++;
            }
            else if(overlap_list.list[j].y_id > R_INF.paf[i].buffer[inner_j].tn)
            {
                inner_j++;
            }
            else
            {
                if(overlap_list.list[j].y_pos_strand == R_INF.paf[i].buffer[inner_j].rev)
                {

                    overlap_list.list[j].is_match = 1;
                    overlap_list.list[j].strong = R_INF.paf[i].buffer[inner_j].ml;
                    overlap_list.list[j].without_large_indel = R_INF.paf[i].buffer[inner_j].no_l_indel;
                    overlap_list.mapped_overlaps_length++;

                    if(overlap_list.list[j].strong == 1)
                    {
                        matched_overlap_1++;
                    }
                    else if(overlap_list.list[j].strong == 0)
                    {
                        matched_overlap_0++;
                    }
                    else
                    {
                        fprintf(stderr, "error\n");
                    }


                    if(overlap_list.list[j].without_large_indel == 0)
                    {
                        num_second_correct_base++;
                    }



                }
                j++;
                inner_j++;
            }
        }


        ///recover missing exact overlaps 
        reverse_complement(g_read.seq, g_read.length);
        for (j = 0; j < overlap_list.length; j++)
        {
            if (overlap_list.list[j].is_match != 1)
            {
                if(overlap_list.list[j].y_pos_strand == 0)
                {
                    recover_UC_Read(&overlap_read, &R_INF, overlap_list.list[j].y_id);
                }
                else
                {
                    recover_UC_Read_RC(&overlap_read, &R_INF, overlap_list.list[j].y_id);
                }

                if(if_exact_match(g_read.seq, g_read.length, overlap_read.seq, overlap_read.length, 
                overlap_list.list[j].x_pos_s, overlap_list.list[j].x_pos_e, 
                overlap_list.list[j].y_pos_s, overlap_list.list[j].y_pos_e))
                {
                    overlap_list.list[j].is_match = 1;
                    overlap_list.list[j].strong = 0;
                    overlap_list.list[j].without_large_indel = 1;
                    overlap_list.mapped_overlaps_length++;
                    potiental_matched_overlap_0++;
                }
                
            }
        }


        if(R_INF.paf[i].is_fully_corrected)
        {
            potiental_matched_overlap_1++;
        }
        



        num_correct_base += 
        push_final_overlaps(&(R_INF.paf[i]), R_INF.reverse_paf, &overlap_list, &g_read, &overlap_read);
    }

    finish_output_buffer();

    destory_buffer_sub_block(&current_sub_buffer);
    destory_Candidates_list(&l);
    destory_overlap_region_alloc(&overlap_list);
    destory_Heap(&heap);
    destory_k_mer_pos_list_alloc(&array_list);
    destory_Graph(&POA_Graph);
    destory_Graph(&DAGCon);    
    destory_UC_Read(&g_read);
    destory_UC_Read(&overlap_read);
    destory_Cigar_record(&current_cigar);
    destory_Correct_dumy(&correct);
    destoryHaplotypeEvdience(&hap);
    destory_Round2_alignment(&second_round);
    destory_small_hash_table(&forward);
    destory_small_hash_table(&reverse);



    pthread_mutex_lock(&statistics);
    total_matched_overlap_0 += matched_overlap_0;
    total_matched_overlap_1 += matched_overlap_1;
    total_potiental_matched_overlap_0 += potiental_matched_overlap_0;
    total_potiental_matched_overlap_1 += potiental_matched_overlap_1;
    total_num_correct_base += num_correct_base;
    total_second_num_correct_base += num_second_correct_base;


    complete_threads++;
    if(complete_threads == thread_num)
    {


        fprintf(stderr, "overlaps with large indels: %llu\n", total_second_num_correct_base);
        fprintf(stderr, "weak overlaps: %llu\n", total_matched_overlap_0);
        fprintf(stderr, "strong overlaps: %llu\n", total_matched_overlap_1); 
        fprintf(stderr, "recover weak overlaps: %llu\n", total_potiental_matched_overlap_0);
        fprintf(stderr, "final available overlaps: %llu\n", total_num_correct_base); 

        fprintf(stderr, "fully corrected reads: %llu\n", total_potiental_matched_overlap_1);

        
              
    }
    pthread_mutex_unlock(&statistics);

    free(arg);
}


void Output_PAF()
{

    fprintf(stdout, "Writing PAF to disk ...... \n");
    char* paf_name = (char*)malloc(strlen(output_file_name)+5);
    sprintf(paf_name, "%s.paf", output_file_name);
    FILE* output_file = fopen(paf_name, "w");
    long long i, j;
    ma_hit_t_alloc* sources = R_INF.paf;



    for (i = 0; i < R_INF.total_reads; i++)
    {
        for (j = 0; j < sources[i].length; j++)
        {
            fwrite(Get_NAME(R_INF, Get_qn(sources[i].buffer[j])), 1, 
            Get_NAME_LENGTH(R_INF, Get_qn(sources[i].buffer[j])), output_file);
            fwrite("\t", 1, 1, output_file);
            fprintf(output_file, "%d\t", Get_READ_LENGTH(R_INF, Get_qn(sources[i].buffer[j])));
            fprintf(output_file, "%d\t", Get_qs(sources[i].buffer[j]));
            fprintf(output_file, "%d\t", Get_qe(sources[i].buffer[j]));
            if(sources[i].buffer[j].rev)
            {
                fprintf(output_file, "-\t");
            }
            else
            {
                fprintf(output_file, "+\t");
            }
            fwrite(Get_NAME(R_INF, Get_tn(sources[i].buffer[j])), 1, 
            Get_NAME_LENGTH(R_INF, Get_tn(sources[i].buffer[j])), output_file);
            fwrite("\t", 1, 1, output_file);
            fprintf(output_file, "%d\t", Get_READ_LENGTH(R_INF, Get_tn(sources[i].buffer[j])));
            fprintf(output_file, "%d\t", Get_ts(sources[i].buffer[j]));
            fprintf(output_file, "%d\t", Get_te(sources[i].buffer[j]));
            fprintf(output_file, "%d\t", sources[i].buffer[j].ml);
            fprintf(output_file, "%d\t", sources[i].buffer[j].bl);
            fprintf(output_file, "255\n");
            
        }
    }

    free(paf_name);
    fclose(output_file);
}



void generate_overlaps(int last_round)
{
    double start_time = Get_T();
    roundID = number_of_round - last_round;
    fprintf(stdout, "Calculting final overlaps ...\n");

    Counting_multiple_thr();
    Build_hash_table_multiple_thr();

    ///thread_num = 1;
    pthread_t *_r_threads;

    _r_threads = (pthread_t *)malloc(sizeof(pthread_t)*thread_num);

    int i = 0;

    for (i = 0; i < thread_num; i++)
    {
        int *arg = (int*)malloc(sizeof(*arg));
        *arg = i;
        pthread_create(_r_threads + i, NULL, Final_overlap_calculate_heap_merge, (void*)arg);
    }
    

    for (i = 0; i<thread_num; i++)
        pthread_join(_r_threads[i], NULL);
    free(_r_threads);

    fprintf(stdout, "Final overlaps have been calculated.\n");
    fprintf(stdout, "%-30s%18.2f\n\n", "Final overlaps calculation time:", Get_T() - start_time); 

    ///Output_PAF();
    /**
    build_string_graph(MIN_OVERLAP_COVERAGE, R_INF.paf, R_INF.reverse_paf, R_INF.total_reads, R_INF.read_length, 
    MIN_OVERLAP_LEN, MAX_HANG_LEN, 3, 0.5, 0.7, 0.8, output_file_name, MAX_BUBBLE_DIST);
    **/

    
    build_string_graph_without_clean(MIN_OVERLAP_COVERAGE, R_INF.paf, R_INF.reverse_paf, R_INF.total_reads, R_INF.read_length, 
    MIN_OVERLAP_LEN, MAX_HANG_LEN, c_round, 0.2, 0.8, 0.8, output_file_name, MAX_BUBBLE_DIST, 0, 1);

}





void Correct_Reads(int last_round)
{
    
    if(load_index_from_disk && load_all_data_from_disk(&R_INF.paf, &R_INF.reverse_paf, 
    output_file_name))
    {
        build_string_graph_without_clean(MIN_OVERLAP_COVERAGE, R_INF.paf, R_INF.reverse_paf, 
        R_INF.total_reads, R_INF.read_length, MIN_OVERLAP_LEN, MAX_HANG_LEN, c_round, 0.2, 
        0.8, 0.8, output_file_name, MAX_BUBBLE_DIST, 0, 0);
        exit(1);
    }
    else
    {
        fprintf(stdout, "Cannot find overlap file.\n");
    }
    

    complete_threads = 0;
    total_matched_overlap_0 = 0;
    total_matched_overlap_1 = 0;
    total_potiental_matched_overlap_0 = 0;
    total_potiental_matched_overlap_1 = 0;
    total_num_read_base = 0;
    total_num_correct_base = 0;
    total_second_num_correct_base = 0;


    if(last_round == 0)
    {
        generate_overlaps(last_round);
        return;
    }
        

    roundID = number_of_round - last_round;
    fprintf(stdout, "Error correction: start the %d-th round ...\n", roundID);
    ///only the first round correction can load index from disk
    if (roundID == 0 && load_index_from_disk && load_pre_cauculated_index())
    {
        ;
    }
    else
    {
        Counting_multiple_thr();

        Build_hash_table_multiple_thr();
    }

    
    fprintf(stdout, "Total pos in hash tabe: %d\n", PCB.total_occ);
    fprintf(stdout, "k_mer_min_freq in hashtable: %d\n", k_mer_min_freq);
    fprintf(stdout, "k_mer_max_freq in hashtable: %d\n", k_mer_max_freq);
    
    
    
    ///verify_Position_hash_table();
    Overlap_calculate_multipe_thr();

    fprintf(stdout, "Error correction: the %d-th round has been completed.\n", roundID);

    Correct_Reads(last_round - 1);
}












