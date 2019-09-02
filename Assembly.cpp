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

pthread_mutex_t statistics;
long long total_matched_overlap_0 = 0;
long long total_matched_overlap_1 = 0;
long long total_potiental_matched_overlap_0 = 0;
long long total_potiental_matched_overlap_1 = 0;
long long total_num_read_base = 0;
long long total_num_correct_base = 0;
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

int get_required_read(const char *required_name, long long RID, All_reads* R_INF)
{
    


    int required_name_length = strlen(required_name);
    char* debug_name = Get_NAME((*R_INF), RID);
    int debug_name_length = Get_NAME_LENGTH((*R_INF), RID);
    int i;


    if (required_name_length == debug_name_length)
    {
        for (i = 0; i < debug_name_length; i++)
        {
            if (required_name[i] != debug_name[i])
            {
                break;
            }
        }

        if (i == debug_name_length)
        {
            fprintf(stderr, "required_name: %s\n", required_name);

            return 1;
            ///fprintf(stderr, "read_length: %d\n", R_INF->g_read->length);
            
            /**
            int aviable_overlap_name = 0;
            for (i = 0; i < overlap_list->length; i++)
            {
                long long Len_x = overlap_list->list[i].x_pos_e -  overlap_list->list[i].x_pos_s + 1;

                if (Len_x * OVERLAP_THRESHOLD <=  overlap_list->list[i].align_length)
                {
                    fprintf(stderr, "a_i: %d\n", aviable_overlap_name);
                    fprintf(stderr, "x_pos_s: %d, x_pos_e: %d\n", 
                    overlap_list->list[i].x_pos_s, overlap_list->list[i].x_pos_e);
                    aviable_overlap_name++;

                    debug_name = Get_NAME((*R_INF), overlap_list->list[i].y_id);
                    debug_name_length = Get_NAME_LENGTH((*R_INF), overlap_list->list[i].y_id);

                    int j = 0;
                    for (j = 0; j < debug_name_length; j++)
                    {
                        fprintf(stderr, "%c", debug_name[j]);
                    }
                    fprintf(stderr, "\n");
                    
                }
            }
            **/
        }
        
    }


    return 0;

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
            memcpy(new_read + new_i, pre_read + pre_i, operation_length);
            pre_i = pre_i + operation_length;
            new_i = new_i + operation_length;
        }
        else if (operation == 1)
        {
            
            for (j = 0; j < operation_length; j++)
            {
                new_read[new_i] = Get_MisMatch_Base(cigar->lost_base[diff_char_i]);
                new_i++;
                diff_char_i++;
            }
            pre_i = pre_i + operation_length;
            
            /**
            pre_i = pre_i + operation_length;
            memcpy(new_read + new_i, cigar->lost_base + diff_char_i, operation_length);
            new_i = new_i + operation_length;
            diff_char_i = diff_char_i + operation_length;
            **/
            
        }
        else if (operation == 3)
        {
            pre_i = pre_i + operation_length;
            diff_char_i = diff_char_i + operation_length;
        }
        else if (operation == 2)
        {
            memcpy(new_read + new_i, cigar->lost_base + diff_char_i, operation_length);
            new_i = new_i + operation_length;
            diff_char_i = diff_char_i + operation_length;
        }
    }
    *new_length = new_i;
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

    /**
    fprintf(stderr, "total_errors: %d, correct_base: %d, length: %d, lost_base_length: %d\n", 
    total_errors, correct_base, cigar->length, cigar->lost_base_length);
    **/

    if (pre_i != pre_length || new_i != new_length)
    {
        fprintf(stderr, "pre_i: %d, pre_length: %d\n", pre_i, pre_length);
        fprintf(stderr, "new_i: %d, new_length: %d\n", new_i, new_length);
    }

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

void just_debug(overlap_region_alloc* overlap_list, All_reads* R_INF)
{
    int j;
    int i;
    int y_id;
    int y_strand;
    int y_readLen;
    int overlap_length;
    int high_quality_overlaps = 0;
    UC_Read g_read;
    init_UC_Read(&g_read);

    for (j = 0; j < overlap_list->length; j++)
    {
        y_id = overlap_list->list[j].y_id;
        y_strand = overlap_list->list[j].y_pos_strand;
        y_readLen = Get_READ_LENGTH((*R_INF), y_id);
        overlap_length = overlap_list->list[j].x_pos_e - overlap_list->list[j].x_pos_s + 1;
        if (overlap_length * OVERLAP_THRESHOLD <=  overlap_list->list[j].align_length)
        {
            high_quality_overlaps++;

            fprintf(stderr, "i: %d, overlap_length: %d, align_length: %d, y_strand: %d\n", 
                    high_quality_overlaps, overlap_length, overlap_list->list[j].align_length, y_strand);
            fprintf(stderr, "x_pos_s: %d, x_pos_e: %d\n", 
                    overlap_list->list[j].x_pos_s, overlap_list->list[j].x_pos_e);

            fprintf(stderr, "%.*s\n\n", Get_NAME_LENGTH((*R_INF), y_id), Get_NAME((*R_INF), y_id));

        }
    }


    recover_UC_Read_RC(&g_read, R_INF, overlap_list->list[0].x_id);

    for (i = 0; i < g_read.length && i < 81; i++)
    {
        fprintf(stderr, "%c", g_read.seq[i]);
    }

    fprintf(stderr, "\n");

    fprintf(stderr, "x_length: %d\n", g_read.length);


    

    fprintf(stderr, "high_quality_overlaps: %d, overlap_list->length: %d\n", high_quality_overlaps,overlap_list->length);

    destory_UC_Read(&g_read);

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
    long long j;

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

    for (i = thr_ID; i < R_INF.total_reads; i = i + thread_num)
    {
        
        
        

        clear_Cigar_record(&current_cigar);
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
                        /************需要注释掉**********/
                        // debug_overlap = debug_overlap + list_length;
                        /************需要注释掉**********/
                        append_k_mer_pos_list_alloc(&array_list, list, list_length, end_pos, 0);
                        ///append_k_mer_pos_list_alloc_prefilter(&array_list, list, list_length, end_pos, 0, &g_read, &R_INF, &correct);
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
                        /************需要注释掉**********/
                        // debug_overlap = debug_overlap + list_length;
                        /************需要注释掉**********/
                        append_k_mer_pos_list_alloc(&array_list, list, list_length, end_pos, 1);
                        ///append_k_mer_pos_list_alloc_prefilter(&array_list, list, list_length, end_pos, 1, &g_read, &R_INF, &correct);
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

        


        /************需要注释掉**********/
        // for (int ijk = 0; ijk < array_list.length; ijk++)
        // {
        //     filtered_debug_overlap = filtered_debug_overlap + array_list.list[ijk].length;
        // }
        /************需要注释掉**********/
        
        ///merge_k_mer_pos_list_alloc(&array_list, &l);
        merge_k_mer_pos_list_alloc_heap_sort(&array_list, &l, &heap);
        
        // if (array_list.length < 3)
        // {
        //     merge_k_mer_pos_list_alloc(&array_list, &l);
        // }
        // else
        // {
        //     merge_k_mer_pos_list_alloc_heap_sort_advance(&array_list, &l, &heap);
        // }
        

       ///以x_pos_e，即结束位置为主元排序
       calculate_overlap_region(&l, &overlap_list, i, g_read.length, &R_INF);

       ///clear_Graph(&POA_Graph);

       correct_overlap(&overlap_list, &R_INF, &g_read, &correct, &overlap_read, &POA_Graph,
       &matched_overlap_0, &matched_overlap_1, &potiental_matched_overlap_0, &potiental_matched_overlap_1,
       &current_cigar, &hap);
      
        num_read_base = num_read_base + g_read.length;
        num_correct_base = num_correct_base + correct.corrected_base;

        push_cigar(R_INF.cigars, i, &current_cigar);

        /**
        if(memcmp("m54334_180926_225337/39780640/ccs", Get_NAME(R_INF, i), Get_NAME_LENGTH(R_INF, i)) == 0)
        {
            fprintf(stderr, "i: %d\n", i);
            ///fprintf(stderr, "overlap_list.length: %d\n", overlap_list.length);
            just_debug(&overlap_list, &R_INF);
        }
        **/
        
        



        ///output_read_to_buffer(i, &R_INF, correct.corrected_read, correct.corrected_read_length, &current_sub_buffer);
        /**
        debug_cigar(&current_cigar, g_read.seq, g_read.length, correct.corrected_read, correct.corrected_read_length, 
        correct.corrected_base);
        **/
        
        
        


        
        
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
    destory_UC_Read(&g_read);
    destory_UC_Read(&overlap_read);
    destory_Cigar_record(&current_cigar);

    destory_Correct_dumy(&correct);

    destoryHaplotypeEvdience(&hap);


    pthread_mutex_lock(&statistics);
    total_matched_overlap_0 += matched_overlap_0;
    total_matched_overlap_1 += matched_overlap_1;
    total_potiental_matched_overlap_0 += potiental_matched_overlap_0;
    total_potiental_matched_overlap_1 += potiental_matched_overlap_1;
    total_num_read_base += num_read_base;
    total_num_correct_base += num_correct_base;

    complete_threads++;
    if(complete_threads == thread_num)
    {
        fprintf(stderr, "total_matched_overlap_0: %llu\n", total_matched_overlap_0);
        fprintf(stderr, "total_matched_overlap_1: %llu\n", total_matched_overlap_1);
        fprintf(stderr, "total_potiental_matched_overlap_0: %llu\n", total_potiental_matched_overlap_0);
        fprintf(stderr, "total_potiental_matched_overlap_1: %llu\n", total_potiental_matched_overlap_1);
        fprintf(stderr, "total_num_read_base: %llu\n", total_num_read_base);
        fprintf(stderr, "total_num_correct_base: %llu\n", total_num_correct_base);
    }
	pthread_mutex_unlock(&statistics);

    free(arg);
}


void* Save_corrected_reads(void* arg)
{
    int thr_ID = *((int*)arg);
    long long i, j;
    UC_Read g_read;
    init_UC_Read(&g_read);

    int new_read_size = 10000;
    char* new_read = (char*)malloc(new_read_size);
    Cigar_record cigar;
    int new_read_length;
    uint64_t N_occ;
    
    for (i = thr_ID; i < R_INF.total_reads; i = i + thread_num)
    {
        recover_UC_Read(&g_read, &R_INF, i);

        if(R_INF.cigars[i].new_length>new_read_size)
        {
            new_read_size = R_INF.cigars[i].new_length;
            new_read = (char*)realloc(new_read, new_read_size);
        }

        cigar.length = R_INF.cigars[i].length;
        cigar.lost_base_length = R_INF.cigars[i].lost_base_length;
        cigar.record = R_INF.cigars[i].record;
        cigar.lost_base = R_INF.cigars[i].lost_base;

        get_corrected_read_from_cigar(&cigar, g_read.seq, g_read.length, new_read, &new_read_length);

        ///need modification
        reverse_complement(new_read, new_read_length);


        N_occ = 0;
        for (j = 0; j < new_read_length; j++)
        {
            if(new_read[j] == 'N')
            {
                N_occ++;
            }
        }
        
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
    free(new_read);
    free(arg);
}


void Output_corrected_reads()
{
    long long i, j;
    UC_Read g_read;
    init_UC_Read(&g_read);
    FILE* output_file = fopen(output_file_name, "w");


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
    
    


    destory_UC_Read(&g_read);
    fclose(output_file);
}






void Overlap_calculate_multipe_thr()
{


    
    double start_time = Get_T();
    pthread_t outputResultSinkHandle;
    /**
    if (roundID == number_of_round - 1)
    {
        init_output_buffer(thread_num);
        pthread_create(&outputResultSinkHandle, NULL, pop_buffer, NULL);
    }
    **/

    

    fprintf(stdout, "Begin Overlap Calculate ...... \n");

    pthread_t *_r_threads;

	_r_threads = (pthread_t *)malloc(sizeof(pthread_t)*thread_num);

    int i = 0;

	for (i = 0; i < thread_num; i++)
	{
		int *arg = (int*)malloc(sizeof(*arg));
		*arg = i;

		pthread_create(_r_threads + i, NULL, Overlap_calculate_heap_merge, (void*)arg);

	}
    

    for (i = 0; i<thread_num; i++)
		pthread_join(_r_threads[i], NULL);

    /**
    if (roundID == number_of_round - 1)
    {
        pthread_join(outputResultSinkHandle, NULL);
        destory_output_buffer();
        ///destory_All_reads(&R_INF);
    }
    **/
    

    free(_r_threads);


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
        //pthread_create(_r_threads + i, NULL, Overlap_calculate, (void*)arg);

	}
    for (i = 0; i<thread_num; i++)
		pthread_join(_r_threads[i], NULL);
    free(_r_threads);
    fprintf(stdout, "%-30s%18.2f\n\n", "Save corrected read time:", Get_T() - start_time);   
    ///fflush(stdout);

    
    ///only the last round can output read to disk
    if (roundID == number_of_round - 1)
    {
        start_time = Get_T();
        
        Output_corrected_reads();
        
        fprintf(stdout, "%-30s%18.2f\n\n", "Output time:", Get_T() - start_time);    
    }
    

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


void Correct_Reads(int last_round)
{
    complete_threads = 0;
    total_matched_overlap_0 = 0;
    total_matched_overlap_1 = 0;
    total_potiental_matched_overlap_0 = 0;
    total_potiental_matched_overlap_1 = 0;
    total_num_read_base = 0;
    total_num_correct_base = 0;


    if(last_round == 0)
        return;

    roundID = number_of_round - last_round;
    fprintf(stdout, "Error correction: start the %d-th round ...\n", roundID);
    ///only the first round correction can load index from disk
    if (roundID == 0 & load_index_from_disk && load_pre_cauculated_index())
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












