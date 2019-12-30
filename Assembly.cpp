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

void* Perform_Counting(void* arg)
{
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
                    k_mer_append(&k_code,code, asm_opt.k_mer_length);
                    avalible_k++;
                    if (avalible_k >= asm_opt.k_mer_length)
                    {
                        if(insert_Total_Count_Table(&TCB, &k_code, asm_opt.k_mer_length))
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

    return NULL;
}

void* Perform_Counting_non_first(void* arg)
{
    int thr_ID = *((int*)arg);

    uint64_t i = 0;
    HPC_seq HPC_read;

    long long select_k_mer_number = 0 ;
    long long k_mer_number = 0 ;

    uint64_t code;

    uint64_t end_pos;

    Hash_code k_code;

    int avalible_k = 0;

    UC_Read g_read;
    init_UC_Read(&g_read);

    for (i = thr_ID; i < R_INF.total_reads; i = i + asm_opt.thread_num)
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
                k_mer_append(&k_code,code, asm_opt.k_mer_length);
                avalible_k++;
                if (avalible_k >= asm_opt.k_mer_length)
                {
                    if(insert_Total_Count_Table(&TCB, &k_code, asm_opt.k_mer_length))
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
    free(arg);

    return NULL;
}

void* Build_hash_table_non_first(void* arg)
{
    int thr_ID = *((int*)arg);

    uint64_t i = 0;
    HPC_seq HPC_read;
    uint64_t code;
    uint64_t end_pos;

    ///long long HPC_base;

    Hash_code k_code;

    int avalible_k = 0;

    UC_Read g_read;
    init_UC_Read(&g_read);

    for (i = thr_ID; i < R_INF.total_reads; i = i + asm_opt.thread_num)
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
                k_mer_append(&k_code,code, asm_opt.k_mer_length);
                avalible_k++;
                if (avalible_k >= asm_opt.k_mer_length)
                {
                    ///there are two requirements
                    ///1. hash(k-mer) 
                    ///2. occ(k-mer) 
                    ///TCB just meet the first requirement，while PCB needs to meet both of them
                    insert_Total_Pos_Table(&PCB, &k_code, asm_opt.k_mer_length, i, end_pos);
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
    free(arg);

    return NULL;
}


void* Build_hash_table(void* arg)
{

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
                    k_mer_append(&k_code,code, asm_opt.k_mer_length);
                    avalible_k++;
                    if (avalible_k >= asm_opt.k_mer_length)
                    {
                       ///there are two requirements
                       ///1. hash(k-mer) 
                       ///2. occ(k-mer) 
                       ///TCB just meet the first requirement，while PCB needs to meet both of them
                        insert_Total_Pos_Table(&PCB, &k_code, asm_opt.k_mer_length, curr_sub_block.read[i].ID, end_pos);
                    }
                }
                else
                {
                    avalible_k = 0;
                    init_Hash_code(&k_code);
                }
            }

            ///load read
            compress_base(Get_READ(R_INF, curr_sub_block.read[i].ID),
            curr_sub_block.read[i].seq.s, curr_sub_block.read[i].seq.l, 
            &R_INF.N_site[curr_sub_block.read[i].ID], HPC_read.N_occ);
            
            memcpy(R_INF.name+R_INF.name_index[curr_sub_block.read[i].ID], 
            curr_sub_block.read[i].name.s, curr_sub_block.read[i].name.l);
        }
    }

    destory_R_buffer_block(&curr_sub_block);
    free(arg);

    return NULL;
}


void Counting_multiple_thr()
{

    double start_time = Get_T();

    fprintf(stderr, "Begin Counting... \n");

    init_Total_Count_Table(asm_opt.k_mer_length, &TCB);

    pthread_t inputReadsHandle;

    int *is_insert = (int*)malloc(sizeof(*is_insert));
	*is_insert = 1;

    if (asm_opt.roundID == 0)
    {
        init_gz_files(&asm_opt);
        init_All_reads(&R_INF);
        init_R_buffer(asm_opt.thread_num);
        pthread_create(&inputReadsHandle, NULL, input_reads_muti_threads, (void*)is_insert);
    }
    

    pthread_t *_r_threads;

	_r_threads = (pthread_t *)malloc(sizeof(pthread_t)*asm_opt.thread_num);

    int i = 0;

	for (i = 0; i < asm_opt.thread_num; i++)
	{
		int *arg = (int*)malloc(sizeof(*arg));
		*arg = i;
        if (asm_opt.roundID == 0)
        {
            pthread_create(_r_threads + i, NULL, Perform_Counting, (void*)arg);
        }
        else
        {
            pthread_create(_r_threads + i, NULL, Perform_Counting_non_first, (void*)arg);
        }
	}


    for (i = 0; i < asm_opt.thread_num; i++)
		pthread_join(_r_threads[i], NULL);

    free(_r_threads);

    ///destory_R_buffer();
    
    ///destory_Total_Count_Table(&TCB);
    if (asm_opt.roundID == 0)
    {
        pthread_join(inputReadsHandle, NULL);
        ///destory_kseq();
        destory_gz_files();
    }

    fprintf(stderr, "Counting has been completed.\n");

    fprintf(stderr, "%-30s%18.2f\n\n", "Counting time:", Get_T() - start_time);
   
    free(is_insert);

}


void Build_hash_table_multiple_thr()
{
    double start_time = Get_T();

    fprintf(stderr, "Begin building hash table... \n");

    init_Total_Pos_Table(&PCB, &TCB);

    double T_start_time = Get_T();

    Traverse_Counting_Table(&TCB, &PCB, asm_opt.k_mer_min_freq, asm_opt.k_mer_max_freq);

    fprintf(stderr, "%-30s%18.2f\n\n", "Traverse time:", Get_T() - T_start_time);
   
    ///at this moment, TCB can be free
    destory_Total_Count_Table(&TCB);

    pthread_t inputReadsHandle;

    int *is_insert = (int*)malloc(sizeof(*is_insert));
	*is_insert = 0;

    if (asm_opt.roundID == 0)
    {
        init_gz_files(&asm_opt);
        clear_R_buffer();
        malloc_All_reads(&R_INF);
        pthread_create(&inputReadsHandle, NULL, input_reads_muti_threads, (void*)is_insert);
    }



    pthread_t *_r_threads;

	_r_threads = (pthread_t *)malloc(sizeof(pthread_t) * asm_opt.thread_num);

    int i = 0;

	for (i = 0; i < asm_opt.thread_num; i++)
	{
		int *arg = (int*)malloc(sizeof(*arg));
		*arg = i;

        if (asm_opt.roundID == 0)
        {
            pthread_create(_r_threads + i, NULL, Build_hash_table, (void*)arg);
        }
        else
        {
            pthread_create(_r_threads + i, NULL, Build_hash_table_non_first, (void*)arg);
        }
	}

    if (asm_opt.roundID == 0)
    {
        pthread_join(inputReadsHandle, NULL);
    }

    for (i = 0; i < asm_opt.thread_num; i++)
		pthread_join(_r_threads[i], NULL);

    free(_r_threads);

    

    fprintf(stderr, "Hash table has been built.\n");

    fprintf(stderr, "%-30s%18.2f\n\n", "Build hash table time:", Get_T() - start_time);

    if (asm_opt.roundID == 0)
    {
        ///destory_kseq();
        destory_gz_files();
        destory_R_buffer();
    }
    
    free(is_insert);
}


void get_corrected_read_from_cigar(Cigar_record* cigar, char* pre_read, int pre_length,
char* new_read, int* new_length)
{
    int i, j;
    int pre_i, new_i;
    int operation, operation_length;
    pre_i = new_i = 0;
    int diff_char_i = 0;


    for (i = 0; i < (long long)cigar->length; i++)
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


    for (i = 0; i < (long long)cigar->length; i++)
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
    for (i = 0; i < (long long)cigar->length; i++)
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
    for (i = 0; i < (long long)cigar->length; i++)
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


    for (i = 0; i < (long long)cigar->length; i++)
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

    if((int)cigar->new_read_length != new_length)
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


void push_overlaps(ma_hit_t_alloc* paf, overlap_region_alloc* overlap_list, int flag, 
All_reads* R_INF, int if_reverse)
{
    long long i = 0, xLen, yLen;
    ma_hit_t tmp;
    clear_ma_hit_t_alloc(paf);
    for (i = 0; i < (long long)overlap_list->length; i++)
    {
        if (overlap_list->list[i].is_match == flag)
        {
            xLen = Get_READ_LENGTH((*R_INF), overlap_list->list[i].x_id);
            yLen = Get_READ_LENGTH((*R_INF), overlap_list->list[i].y_id);

            tmp.qns = overlap_list->list[i].x_id;
            tmp.qns = tmp.qns << 32;
            tmp.tn = overlap_list->list[i].y_id;
        
            if(if_reverse != 0)
            {
                tmp.qns = tmp.qns | (uint64_t)(xLen - overlap_list->list[i].x_pos_s - 1);
                tmp.qe = xLen - overlap_list->list[i].x_pos_e - 1;
                tmp.ts = yLen - overlap_list->list[i].y_pos_s - 1;
                tmp.te = yLen - overlap_list->list[i].y_pos_e - 1;
            }
            else
            {
                tmp.qns = tmp.qns | (uint64_t)(overlap_list->list[i].x_pos_s);
                tmp.qe = overlap_list->list[i].x_pos_e;
                tmp.ts = overlap_list->list[i].y_pos_s;
                tmp.te = overlap_list->list[i].y_pos_e;
            }
            
    
            ///for overlap_list, the x_strand of all overlaps are 0, so the tmp.rev is the same as the y_strand
            tmp.rev = overlap_list->list[i].y_pos_strand;

            ///tmp.bl = R_INF.read_length[overlap_list->list[i].y_id];
            tmp.bl = Get_READ_LENGTH((*R_INF), overlap_list->list[i].y_id);
            tmp.ml = overlap_list->list[i].strong;
            tmp.no_l_indel = overlap_list->list[i].without_large_indel;

            add_ma_hit_t_alloc(paf, &tmp);
        }
    }
    
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
overlap_region_alloc* overlap_list, int flag)
{
    long long i = 0;
    long long available_overlaps = 0;
    ma_hit_t tmp;
    clear_ma_hit_t_alloc(paf);
    for (i = 0; i < (long long)overlap_list->length; i++)
    {
        if (overlap_list->list[i].is_match == flag)
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

            tmp.el = overlap_list->list[i].shared_seed;

            add_ma_hit_t_alloc(paf, &tmp);
        }
    }


    return available_overlaps;
    
}

void get_new_candidates(long long readID, UC_Read* g_read, overlap_region_alloc* overlap_list, k_mer_pos_list_alloc* array_list,
HeapSq* heap, Candidates_list* l, double band_width_threshold, int keep_whole_chain)
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
            k_mer_append(&k_code,code, asm_opt.k_mer_length);
            avalible_k++;
            if (avalible_k >= asm_opt.k_mer_length)
            {
                list_length = locate_Total_Pos_Table(&PCB, &k_code, &list, asm_opt.k_mer_length, &sub_ID);

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
            k_mer_append(&k_code,code, asm_opt.k_mer_length);
            avalible_k++;
            if (avalible_k >= asm_opt.k_mer_length)
            {
                list_length = locate_Total_Pos_Table(&PCB, &k_code, &list, asm_opt.k_mer_length, &sub_ID);
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
    
    calculate_overlap_region_by_chaining(l, overlap_list, readID, g_read->length, &R_INF, 
    band_width_threshold, keep_whole_chain);
}


void* Overlap_calculate_heap_merge(void* arg)
{
    long long num_read_base = 0;
    long long num_correct_base = 0;
    long long num_recorrect_base = 0;
    int fully_cov, abnormal;

    int thr_ID = *((int*)arg);
    long long i = 0;

    UC_Read g_read;
    init_UC_Read(&g_read);

    UC_Read overlap_read;
    init_UC_Read(&overlap_read);

    Candidates_list l;
    Graph POA_Graph;
    Graph DAGCon;
    init_Graph(&DAGCon);
    init_Graph(&POA_Graph);
    init_Candidates_list(&l);
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

    for (i = thr_ID; i < (long long)R_INF.total_reads; i = i + asm_opt.thread_num)
    {
        ///get_new_candidates(i, &g_read, &overlap_list, &array_list, &heap, &l, THRESHOLD_RATE*1.5);
        get_new_candidates(i, &g_read, &overlap_list, &array_list, &heap, &l, 0.02, 1);

        clear_Cigar_record(&current_cigar);
        clear_Round2_alignment(&second_round);

        correct_overlap(&overlap_list, &R_INF, &g_read, &correct, &overlap_read, &POA_Graph, &DAGCon,
        &current_cigar, &hap, &second_round, 0, 1, &fully_cov, &abnormal);

        num_read_base += g_read.length;
        num_correct_base += correct.corrected_base;
        num_recorrect_base += second_round.dumy.corrected_base;

        push_cigar(R_INF.cigars, i, &current_cigar);
        push_cigar(R_INF.second_round_cigar, i, &(second_round.cigar));


        R_INF.paf[i].is_fully_corrected = 0;
        if(fully_cov)
        {
           if(get_cigar_errors(&current_cigar) == 0 && 
            get_cigar_errors(&second_round.cigar) == 0)
            {
                R_INF.paf[i].is_fully_corrected = 1;
            }     
        }
        R_INF.paf[i].is_abnormal = abnormal;

        push_overlaps(&(R_INF.paf[i]), &overlap_list, 1, &R_INF, asm_opt.roundID%2);
        push_overlaps(&(R_INF.reverse_paf[i]), &overlap_list, 2, &R_INF, asm_opt.roundID%2);
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
    


    pthread_mutex_lock(&statistics);
    asm_opt.num_bases += num_read_base;
    asm_opt.num_corrected_bases += num_correct_base;
    asm_opt.num_recorrected_bases += num_recorrect_base;

    asm_opt.complete_threads++;
    if(asm_opt.complete_threads == asm_opt.thread_num)
    {
        fprintf(stderr, "total bases #: %lld\n", asm_opt.num_bases);
        fprintf(stderr, "total corrected bases: %lld\n", asm_opt.num_corrected_bases);
        fprintf(stderr, "total recorrected bases: %lld\n", asm_opt.num_recorrected_bases);
    }
	pthread_mutex_unlock(&statistics);
    free(arg);

    return NULL;
}

void* Output_related_reads(void* arg)
{
    int thr_ID = *((int*)arg);
    long long i = 0;

    UC_Read g_read;
    init_UC_Read(&g_read);

    UC_Read overlap_read;
    init_UC_Read(&overlap_read);

    Candidates_list l;
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


    long long required_read_name_length = strlen(asm_opt.required_read_name);
    for (i = thr_ID; i < (long long)R_INF.total_reads; i = i + asm_opt.thread_num)
    {
        
        if(required_read_name_length == (long long)Get_NAME_LENGTH((R_INF),i)
            &&
           memcmp(asm_opt.required_read_name, Get_NAME((R_INF), i), Get_NAME_LENGTH((R_INF),i)) == 0)
        {
            ////get_new_candidates(i, &g_read, &overlap_list, &array_list, &heap, &l, THRESHOLD_RATE*1.5);
            get_new_candidates(i, &g_read, &overlap_list, &array_list, &heap, &l, 0.02, 1);

            fprintf(stderr, ">%.*s\n", (int)Get_NAME_LENGTH((R_INF), i),
            Get_NAME((R_INF), i));
            recover_UC_Read(&g_read, &R_INF, i);
            fprintf(stderr, "%.*s\n", (int)g_read.length, g_read.seq);


            uint64_t k;
            for (k = 0; k < overlap_list.length; k++)
            {
                fprintf(stderr, ">%.*s\n", (int)Get_NAME_LENGTH((R_INF),overlap_list.list[k].y_id),
                Get_NAME((R_INF),overlap_list.list[k].y_id));
                recover_UC_Read(&g_read, &R_INF, overlap_list.list[k].y_id);
                fprintf(stderr, "%.*s\n", (int)g_read.length, g_read.seq);
            }
            
        }
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
    

    free(arg);

    return NULL;
}

inline long long get_N_occ(char* seq, long long length)
{
    long long N_occ = 0;
    long long j;
    for (j = 0; j < length; j++)
    {
        if(seq_nt6_table[(uint8_t)seq[j]] >= 4)
        {
            N_occ++;
        }
    }
    return N_occ;
}


void* Save_corrected_reads(void* arg)
{
    int thr_ID = *((int*)arg);
    long long i;
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
    
    for (i = thr_ID; i < (long long)R_INF.total_reads; i = i + asm_opt.thread_num)
    {
        recover_UC_Read(&g_read, &R_INF, i);

        /********************************1 round******************************/
        if((long long)R_INF.cigars[i].new_length > first_round_read_size)
        {
            first_round_read_size = R_INF.cigars[i].new_length;
            first_round_read = (char*)realloc(first_round_read, first_round_read_size);
        }

        cigar.length = R_INF.cigars[i].length;
        cigar.lost_base_length = R_INF.cigars[i].lost_base_length;
        cigar.record = R_INF.cigars[i].record;
        cigar.lost_base = R_INF.cigars[i].lost_base;

        get_corrected_read_from_cigar(&cigar, g_read.seq, g_read.length, first_round_read, &first_round_read_length);
        
        /********************************1 round******************************/

        /********************************2 round******************************/
        if((long long)R_INF.second_round_cigar[i].new_length > second_round_read_size)
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
        
        /********************************2 round******************************/


        

        new_read = second_round_read;
        new_read_length = second_round_read_length;

        
        if (asm_opt.roundID != asm_opt.number_of_round - 1)
        {
            ///need modification
            reverse_complement(new_read, new_read_length);
        }
        else if(asm_opt.number_of_round % 2 == 0)
        {
            ///need modification
            reverse_complement(new_read, new_read_length);
        }

        
        N_occ = get_N_occ(new_read, new_read_length);

        
        if((long long)R_INF.read_size[i] < new_read_length)
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

    return NULL;
}

void Output_corrected_reads()
{
    long long i;
    UC_Read g_read;
    init_UC_Read(&g_read);
    char* gfa_name = (char*)malloc(strlen(asm_opt.output_file_name)+35);
    sprintf(gfa_name, "%s.ec.fa", asm_opt.output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    free(gfa_name);

    for (i = 0; i < (long long)R_INF.total_reads; i++)
    {
        recover_UC_Read(&g_read, &R_INF, i);
        fwrite(">", 1, 1, output_file);
        fwrite(Get_NAME(R_INF, i), 1, Get_NAME_LENGTH(R_INF, i), output_file);
        fwrite("\n", 1, 1, output_file);
        fwrite(g_read.seq, 1, g_read.length, output_file);
        fwrite("\n", 1, 1, output_file);
    }
    destory_UC_Read(&g_read);
    fclose(output_file);
}


void Overlap_calculate_multipe_thr()
{
    double start_time = Get_T();
    
    fprintf(stderr, "Begin calculating overlaps... \n");

    pthread_t *_r_threads;

	_r_threads = (pthread_t *)malloc(sizeof(pthread_t)*asm_opt.thread_num);

    int i = 0;

	for (i = 0; i < asm_opt.thread_num; i++)
	{
		int *arg = (int*)malloc(sizeof(*arg));
		*arg = i;
        if(!asm_opt.required_read_name)
        {
            pthread_create(_r_threads + i, NULL, Overlap_calculate_heap_merge, (void*)arg);
        }
        else
        {
            pthread_create(_r_threads + i, NULL, Output_related_reads, (void*)arg);
        }
	}
    

    for (i = 0; i < asm_opt.thread_num; i++)
		pthread_join(_r_threads[i], NULL);


    free(_r_threads);

    if(asm_opt.required_read_name)
    {
        exit(1);
    }


    destory_Total_Pos_Table(&PCB);

    fprintf(stderr, "All overlaps have been calculated.\n");

    fprintf(stderr, "%-30s%18.2f\n\n", "Overlap calculation time:", Get_T() - start_time);    

    start_time = Get_T();

    _r_threads = (pthread_t *)malloc(sizeof(pthread_t)*asm_opt.thread_num);
    for (i = 0; i < asm_opt.thread_num; i++)
	{
		int *arg = (int*)malloc(sizeof(*arg));
		*arg = i;
		pthread_create(_r_threads + i, NULL, Save_corrected_reads, (void*)arg);
	}


    for (i = 0; i < asm_opt.thread_num; i++)
		pthread_join(_r_threads[i], NULL);
    free(_r_threads);

    fprintf(stderr, "%-30s%18.2f\n\n", "Corrected read saving time:", Get_T() - start_time);   

    
    ///only the last round can output read to disk
    if (asm_opt.roundID == asm_opt.number_of_round - 1)
    {
        start_time = Get_T();
        
        Output_corrected_reads();
        
        fprintf(stderr, "%-30s%18.2f\n\n", "Output time:", Get_T() - start_time);    
    }
}


// int load_pre_cauculated_index()
// {
//     if(load_Total_Pos_Table(&PCB, asm_opt.read_file_name) && load_All_reads(&R_INF, asm_opt.read_file_name))
//     {
//         return 1;
//     }
//     else
//     {
//         return 0;
//     }
    
// }


void update_overlaps(overlap_region_alloc* overlap_list, ma_hit_t_alloc* paf, 
UC_Read* g_read, UC_Read* overlap_read, int is_match, int is_exact)
{

    uint64_t inner_j = 0;
    uint64_t j = 0;
    long long x_overlapLen, y_overlapLen;
    while (j < overlap_list->length && inner_j < paf->length)
    {
        if(overlap_list->list[j].y_id < paf->buffer[inner_j].tn)
        {
            j++;
        }
        else if(overlap_list->list[j].y_id > paf->buffer[inner_j].tn)
        {
            inner_j++;
        }
        else
        {
            if(overlap_list->list[j].y_pos_strand == paf->buffer[inner_j].rev)
            {
                x_overlapLen = Get_qe(paf->buffer[inner_j]) - Get_qs(paf->buffer[inner_j]) + 1;
                y_overlapLen = Get_te(paf->buffer[inner_j]) - Get_ts(paf->buffer[inner_j]) + 1;
                if(x_overlapLen < y_overlapLen) x_overlapLen = y_overlapLen;
                x_overlapLen = x_overlapLen * 0.1;

                // if(
                // ((DIFF(overlap_list->list[j].x_pos_s, Get_qs(paf->buffer[inner_j])) < x_overlapLen)
                // && (DIFF(overlap_list->list[j].x_pos_e, Get_qe(paf->buffer[inner_j])) < x_overlapLen))
                // ||
                // ((DIFF(overlap_list->list[j].y_pos_s, Get_ts(paf->buffer[inner_j])) < x_overlapLen)
                // && (DIFF(overlap_list->list[j].y_pos_e, Get_te(paf->buffer[inner_j])) < x_overlapLen)))
                if(
                ((DIFF(overlap_list->list[j].x_pos_s, Get_qs(paf->buffer[inner_j])) < (uint64_t)x_overlapLen)
                && (DIFF(overlap_list->list[j].x_pos_e, Get_qe(paf->buffer[inner_j])) < (uint64_t)x_overlapLen))
                || 
                ((DIFF(overlap_list->list[j].y_pos_s, Get_ts(paf->buffer[inner_j])) < (uint64_t)x_overlapLen)
                && (DIFF(overlap_list->list[j].y_pos_e, Get_te(paf->buffer[inner_j])) < (uint64_t)x_overlapLen))
                )
                {
                    overlap_list->list[j].is_match = is_match;
                    overlap_list->list[j].strong = paf->buffer[inner_j].ml;
                    overlap_list->list[j].without_large_indel = paf->buffer[inner_j].no_l_indel;
                    if(is_exact == 1)
                    {
                        if(overlap_list->list[j].y_pos_strand == 0)
                        {
                            recover_UC_Read(overlap_read, &R_INF, overlap_list->list[j].y_id);
                        }
                        else
                        {
                            recover_UC_Read_RC(overlap_read, &R_INF, overlap_list->list[j].y_id);
                        }
                        if(if_exact_match(g_read->seq, g_read->length, overlap_read->seq, overlap_read->length, 
                        overlap_list->list[j].x_pos_s, overlap_list->list[j].x_pos_e, 
                        overlap_list->list[j].y_pos_s, overlap_list->list[j].y_pos_e))
                        {
                            overlap_list->list[j].shared_seed = 1;
                        }
                        else
                        {
                            overlap_list->list[j].shared_seed = 0;
                        }
                    }
                }
                else
                {
                    overlap_list->list[j].is_match = 3;
                }
            }
            else
            {
                overlap_list->list[j].is_match = 3;
            }
            
            j++;
            inner_j++;
        }
    }
}

void update_exact_overlaps(overlap_region_alloc* overlap_list, UC_Read* g_read, UC_Read* overlap_read)
{
    uint64_t j;
    for (j = 0; j < overlap_list->length; j++)
    {
        if (overlap_list->list[j].is_match != 1)
        {
            if(overlap_list->list[j].y_pos_strand == 0)
            {
                recover_UC_Read(overlap_read, &R_INF, overlap_list->list[j].y_id);
            }
            else
            {
                recover_UC_Read_RC(overlap_read, &R_INF, overlap_list->list[j].y_id);
            }

            if(if_exact_match(g_read->seq, g_read->length, overlap_read->seq, overlap_read->length, 
            overlap_list->list[j].x_pos_s, overlap_list->list[j].x_pos_e, 
            overlap_list->list[j].y_pos_s, overlap_list->list[j].y_pos_e))
            {
                overlap_list->list[j].is_match = 1;
                overlap_list->list[j].strong = 0;
                overlap_list->list[j].without_large_indel = 1;
                overlap_list->list[j].shared_seed = 1;
            }
        }
    }
}


void statistic(ma_hit_t_alloc* paf, ma_hit_t_alloc* rev_paf, long long readNum)
{

    long long forward, reverse, strong, weak, exact, no_l_indel;
    no_l_indel = forward = reverse = exact = strong = weak = 0;
    long long i, j;
    
    for (i = 0; i < readNum; i++)
    {
        forward += paf[i].length;
        reverse += rev_paf[i].length;
        for (j = 0; j < paf[i].length; j++)
        {
            if(paf[i].buffer[j].el == 1) exact++;
            if(paf[i].buffer[j].ml == 1) strong++;
            if(paf[i].buffer[j].ml == 0) weak++;
            if(paf[i].buffer[j].no_l_indel == 1) no_l_indel++;
        }
    }


    fprintf(stderr, "****************statistic for overlaps****************\n");
    fprintf(stderr, "overlaps #: %lld\n", forward);
    fprintf(stderr, "strong overlaps #: %lld\n", strong);
    fprintf(stderr, "weak overlaps #: %lld\n", weak);
    fprintf(stderr, "exact overlaps #: %lld\n", exact);
    fprintf(stderr, "inexact overlaps #: %lld\n", forward - exact);
    fprintf(stderr, "overlaps without large indels#: %lld\n", no_l_indel);
    fprintf(stderr, "reverse overlaps #: %lld\n", reverse);
    fprintf(stderr, "****************statistic for overlaps****************\n");
}


void fill_chain(Fake_Cigar* chain, char* x_string, char* y_string, long long xBeg, long long yBeg,
long long x_readLen, long long y_readLen, Cigar_record* cigar, uint8_t* c2n)
{
    long long i, xOffset, yOffset, xRegionLen, yRegionLen, /**bandLen,**/ maxXpos, maxYpos, mapScore, zdroped;
    ///float band_rate = 0.08;
    int endbouns;
    if(chain->length <= 0) return;

    kvec_t(uint8_t) x_num;
    kvec_t(uint8_t) y_num;
    kv_init(x_num);
    kv_init(y_num);

    ///deal with region 0 backward
    i = 0;
    endbouns = 0;

    xOffset = get_fake_gap_pos(chain, 0);
    xOffset = xOffset - 1;
    yOffset = (xOffset - xBeg) + yBeg + get_fake_gap_shift(chain, 0);
    if(xOffset >= 0 && yOffset >= 0)
    {
        xRegionLen = xOffset + 1;
        yRegionLen = yOffset + 1;
        //note here cannot use DIFF(xRegionLen, yRegionLen)
        // bandLen = (MIN(xRegionLen, yRegionLen))*band_rate;
        // if(bandLen == 0) bandLen = MIN(xRegionLen, yRegionLen);

        ///do alignment backward
        kv_resize(uint8_t, x_num, (uint64_t)xRegionLen);
        kv_resize(uint8_t, y_num, (uint64_t)yRegionLen);
        ///text is x, query is y
        afine_gap_alignment(x_string, x_num.a, xRegionLen, y_string, y_num.a, yRegionLen, 
        c2n, BACKWARD_KSW, MATCH_SCORE_KSW, MISMATCH_SCORE_KSW, GAP_OPEN_KSW, GAP_EXT_KSW,
        /**bandLen,**/BAND_KSW, Z_DROP_KSW, endbouns, &maxXpos, &maxYpos, &mapScore, &zdroped);
        // fprintf(stderr, "* xOffset: %lld, yOffset: %lld, xRegionLen: %lld, yRegionLen: %lld, bandLen: %lld, maxXpos: %lld, maxYpos: %lld, zdroped: %lld\n",
        // xOffset, yOffset, xRegionLen, yRegionLen, BAND_KSW, maxXpos, maxYpos, zdroped);
    }

    ///align forward
    for (i = 0; i < (long long)chain->length; i++)
    {
        // xOffset = get_fake_gap_pos(chain, i);
        // yOffset = xOffset + get_fake_gap_shift(chain, i);
        xOffset = get_fake_gap_pos(chain, i);
        yOffset = (xOffset - xBeg) + yBeg + get_fake_gap_shift(chain, i);
        ///last region
        if(i == (long long)(chain->length - 1))
        {
            endbouns = 0;
            xRegionLen = x_readLen - xOffset;
            yRegionLen = y_readLen - yOffset;
            //note here cannot use DIFF(xRegionLen, yRegionLen)
            // bandLen = (MIN(xRegionLen, yRegionLen))*band_rate;
            // if(bandLen == 0) bandLen = MIN(xRegionLen, yRegionLen);
        }
        else
        {
            ///higher endbouns for middle regions
            endbouns = MATCH_SCORE_KSW;
            xRegionLen = get_fake_gap_pos(chain, i+1) - xOffset;
            yRegionLen = (get_fake_gap_pos(chain, i+1) + get_fake_gap_shift(chain, i+1)) - 
                         (get_fake_gap_pos(chain, i) + get_fake_gap_shift(chain, i));

            // bandLen = MAX((MIN(xRegionLen, yRegionLen))*band_rate, DIFF(xRegionLen, yRegionLen));
            // if(bandLen == 0) bandLen = MIN(xRegionLen, yRegionLen);
        }

        
        ///do alignment forward
        kv_resize(uint8_t, x_num, (uint64_t)xRegionLen);
        kv_resize(uint8_t, y_num, (uint64_t)yRegionLen);
        ///text is x, query is y
        afine_gap_alignment(x_string+xOffset, x_num.a, xRegionLen, y_string+yOffset, y_num.a, yRegionLen, 
        c2n, FORWARD_KSW, MATCH_SCORE_KSW, MISMATCH_SCORE_KSW, GAP_OPEN_KSW, GAP_EXT_KSW,
        /**bandLen,**/BAND_KSW, Z_DROP_KSW, endbouns, &maxXpos, &maxYpos, &mapScore, &zdroped);
        // fprintf(stderr, "# xOffset: %lld, yOffset: %lld, xRegionLen: %lld, yRegionLen: %lld, bandLen: %lld, maxXpos: %lld, maxYpos: %lld, zdroped: %lld\n",
        // xOffset, yOffset, xRegionLen, yRegionLen, BAND_KSW, maxXpos, maxYpos, zdroped);
    }


    kv_destroy(x_num);
    kv_destroy(y_num);
}
void Final_phasing(overlap_region_alloc* overlap_list, Cigar_record_alloc* cigarline,
UC_Read* g_read, UC_Read* overlap_read, uint8_t* c2n)
{
    uint64_t i, xLen, yStrand;
    char* x_string;
    char* y_string;
    Cigar_record* cigar;
    resize_Cigar_record_alloc(cigarline, overlap_list->length);
    

    
    for (i = 0; i < overlap_list->length; i++)
    {
        if(overlap_list->list[i].is_match == 1 ||
           overlap_list->list[i].is_match == 2 ||
           overlap_list->list[i].is_match == 3)
        {
            xLen = overlap_list->list[i].x_pos_e - overlap_list->list[i].x_pos_s + 1;
            yStrand = overlap_list->list[i].y_pos_strand;
            cigar = &(cigarline->buffer[i]);
            ///has already been matched exactly
            if(overlap_list->list[i].is_match == 1 && overlap_list->list[i].shared_seed == 1)
            {
                add_cigar_record(g_read->seq + overlap_list->list[i].x_pos_s, xLen, cigar, 0);
            }
            else 
            {
                if(yStrand == 0)
                {
                    recover_UC_Read(overlap_read, &R_INF, overlap_list->list[i].y_id);
                }
                else
                {
                    recover_UC_Read_RC(overlap_read, &R_INF, overlap_list->list[i].y_id);
                }
                x_string = g_read->seq;
                y_string = overlap_read->seq;

                fill_chain(&(overlap_list->list[i].f_cigar), x_string, y_string, 
                overlap_list->list[i].x_pos_s, overlap_list->list[i].y_pos_s,
                Get_READ_LENGTH(R_INF, overlap_list->list[i].x_id), 
                Get_READ_LENGTH(R_INF, overlap_list->list[i].y_id), cigar, c2n);
            }
            
        }
    }
    
}

void* Final_overlap_calculate_heap_merge(void* arg)
{
    int thr_ID = *((int*)arg);
    uint64_t i = 0;

    UC_Read g_read;
    init_UC_Read(&g_read);

    UC_Read overlap_read;
    init_UC_Read(&overlap_read);

    Candidates_list l;    
    init_Candidates_list(&l);

    k_mer_pos_list_alloc array_list;
    init_k_mer_pos_list_alloc(&array_list);

    overlap_region_alloc overlap_list;
    init_overlap_region_alloc(&overlap_list);

    HeapSq heap;
    Init_Heap(&heap);

    Cigar_record_alloc cigarline;
    init_Cigar_record_alloc(&cigarline);

    uint8_t c2n[256];
    memset(c2n, 4, 256);
    c2n[(uint8_t)'A'] = c2n[(uint8_t)'a'] = 0; c2n[(uint8_t)'C'] = c2n[(uint8_t)'c'] = 1;
	c2n[(uint8_t)'G'] = c2n[(uint8_t)'g'] = 2; c2n[(uint8_t)'T'] = c2n[(uint8_t)'t'] = 3; // build the encoding table
    




    for (i = thr_ID; i < R_INF.total_reads; i = i + asm_opt.thread_num)
    {

        get_new_candidates(i, &g_read, &overlap_list, &array_list, &heap, &l, 0.001, 0);
        /**
        correct_overlap(&overlap_list, &R_INF, &g_read, &correct, &overlap_read, &POA_Graph, &DAGCon,
        &matched_overlap_0, &matched_overlap_1, &potiental_matched_overlap_0, &potiental_matched_overlap_1,
        &current_cigar, &hap, &second_round, 0, 0);
        push_final_overlaps(&(R_INF.paf[i]), &overlap_list);
        **/

        overlap_region_sort_y_id(overlap_list.list, overlap_list.length);
        ma_hit_sort_tn(R_INF.paf[i].buffer, R_INF.paf[i].length);
        ma_hit_sort_tn(R_INF.reverse_paf[i].buffer, R_INF.reverse_paf[i].length);
        reverse_complement(g_read.seq, g_read.length);


        update_overlaps(&overlap_list, &(R_INF.paf[i]), &g_read, &overlap_read, 1, 1);
        update_overlaps(&overlap_list, &(R_INF.reverse_paf[i]), &g_read, &overlap_read, 2, 0);
        ///recover missing exact overlaps 
        update_exact_overlaps(&overlap_list, &g_read, &overlap_read);

        ///Final_phasing(&overlap_list, &cigarline, &g_read, &overlap_read, c2n);

        push_final_overlaps(&(R_INF.paf[i]), R_INF.reverse_paf, 
        &overlap_list, 1);
        push_final_overlaps(&(R_INF.reverse_paf[i]), R_INF.reverse_paf, 
        &overlap_list, 2);

        
    }

    finish_output_buffer();

    destory_Candidates_list(&l);
    destory_overlap_region_alloc(&overlap_list);
    destory_Heap(&heap);
    destory_k_mer_pos_list_alloc(&array_list);  
    destory_UC_Read(&g_read);
    destory_UC_Read(&overlap_read);
    destory_Cigar_record_alloc(&cigarline);
    

    pthread_mutex_lock(&statistics);
    asm_opt.complete_threads++;
    if(asm_opt.complete_threads == asm_opt.thread_num)
    {
        if(VERBOSE >= 1)
        {
            statistic(R_INF.paf, R_INF.reverse_paf, R_INF.total_reads);
        }
    }
    pthread_mutex_unlock(&statistics);
    free(arg);

    return NULL;
}


void Output_PAF()
{

    fprintf(stderr, "Writing PAF to disk ...... \n");
    char* paf_name = (char*)malloc(strlen(asm_opt.output_file_name)+5);
    sprintf(paf_name, "%s.ovlp.paf", asm_opt.output_file_name);
    FILE* output_file = fopen(paf_name, "w");
    uint64_t i, j;
    ma_hit_t_alloc* sources = R_INF.paf;



    for (i = 0; i < R_INF.total_reads; i++)
    {
        for (j = 0; j < sources[i].length; j++)
        {
            fwrite(Get_NAME(R_INF, Get_qn(sources[i].buffer[j])), 1, 
            Get_NAME_LENGTH(R_INF, Get_qn(sources[i].buffer[j])), output_file);
            fwrite("\t", 1, 1, output_file);
            fprintf(output_file, "%lu\t", (unsigned long)Get_READ_LENGTH(R_INF, Get_qn(sources[i].buffer[j])));
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
            fprintf(output_file, "%lu\t", (unsigned long)Get_READ_LENGTH(R_INF, Get_tn(sources[i].buffer[j])));
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


int check_cluster(uint64_t* list, long long listLen, ma_hit_t_alloc* paf, float threshold)
{
    long long i, k;
    uint32_t qn, tn;
    long long T_edges, A_edges;
    T_edges = A_edges = 0;
    for (i = 0; i < listLen; i++)
    {
        qn = (uint32_t)list[i];
        for (k = i + 1; k < listLen; k++)
        {
            tn = (uint32_t)list[k];
            if(get_specific_overlap(&(paf[qn]), qn, tn) != -1)
            {
                A_edges++;
            }
            
            if(get_specific_overlap(&(paf[tn]), tn, qn) != -1)
            {
                A_edges++;
            }

            T_edges = T_edges + 2;
        }
        
    }

    if(A_edges >= (T_edges*threshold))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}


void rescue_edges(ma_hit_t_alloc* paf, ma_hit_t_alloc* rev_paf, 
long long readNum, long long rescue_threshold, float cluster_threshold)
{
    double startTime = Get_T();
    long long i, j, revises = 0;
    uint32_t qn, tn;

    kvec_t(uint64_t) edge_vector;
    kv_init(edge_vector);
    kvec_t(uint64_t) edge_vector_index;
    kv_init(edge_vector_index);
    uint64_t flag;
    int index;
    
    for (i = 0; i < readNum; i++)
    {
        edge_vector.n = 0;
        edge_vector_index.n = 0;
        for (j = 0; j < paf[i].length; j++)
        {
            qn = Get_qn(paf[i].buffer[j]);
            tn = Get_tn(paf[i].buffer[j]);
            index = get_specific_overlap(&(rev_paf[tn]), tn, qn);
            if(index != -1)
            {
                flag = tn;
                flag = flag << 32;
                flag = flag | (uint64_t)(index);
                kv_push(uint64_t, edge_vector, flag);
                kv_push(uint64_t, edge_vector_index, j);
            }
        }

        ///the read itself has these overlaps, but all related reads do not have 
        ///we need to remove all overlaps from paf[i], and then add all overlaps to rev_paf[i]
        if((long long)edge_vector.n >= rescue_threshold && 
        check_cluster(edge_vector.a, edge_vector.n, paf, cluster_threshold) == 1)
        {
            // fprintf(stderr,"\nremove following %u edges...\n", edge_vector.n);
            // print_revise_edges(&(paf[i]), edge_vector_index.a, edge_vector_index.n);

            add_overlaps(&(paf[i]), &(rev_paf[i]), edge_vector_index.a, edge_vector_index.n);
            remove_overlaps(&(paf[i]), edge_vector_index.a, edge_vector_index.n);
            revises = revises + edge_vector.n;
        }

        edge_vector.n = 0;
        edge_vector_index.n = 0;
        for (j = 0; j < rev_paf[i].length; j++)
        {
            qn = Get_qn(rev_paf[i].buffer[j]);
            tn = Get_tn(rev_paf[i].buffer[j]);
            index = get_specific_overlap(&(paf[tn]), tn, qn);
            if(index != -1)
            {
                flag = tn;
                flag = flag << 32;
                flag = flag | (uint64_t)(index);
                kv_push(uint64_t, edge_vector, flag);
                kv_push(uint64_t, edge_vector_index, j);
            }
        }

        ///the read itself do not have these overlaps, but all related reads have
        ///we need to remove all overlaps from rev_paf[i], and then add all overlaps to paf[i]
        if((long long)edge_vector.n >= rescue_threshold && 
        check_cluster(edge_vector.a, edge_vector.n, paf, cluster_threshold) == 1)
        {
            // fprintf(stderr,"\nadd following %u edges...\n", edge_vector.n);
            // print_revise_edges(&(rev_paf[i]), edge_vector_index.a, edge_vector_index.n);


            remove_overlaps(&(rev_paf[i]), edge_vector_index.a, edge_vector_index.n);
            add_overlaps_from_different_sources(paf, &(paf[i]), edge_vector.a, edge_vector.n);
            revises = revises + edge_vector.n;
        }
    }


    kv_destroy(edge_vector);
    kv_destroy(edge_vector_index);

    fprintf(stderr, "[M::%s] took %0.2fs, revise edges #: %lld\n\n", __func__, Get_T()-startTime, revises);
}


void generate_overlaps(int last_round)
{
    double start_time = Get_T();
    asm_opt.roundID = asm_opt.number_of_round - last_round;
    fprintf(stderr, "Begin calculting final overlaps ...\n");

    Counting_multiple_thr();
    Build_hash_table_multiple_thr();

    pthread_t *_r_threads;

    _r_threads = (pthread_t *)malloc(sizeof(pthread_t) * asm_opt.thread_num);

    int i = 0;

    for (i = 0; i < asm_opt.thread_num; i++)
    {
        int *arg = (int*)malloc(sizeof(*arg));
        *arg = i;
        pthread_create(_r_threads + i, NULL, Final_overlap_calculate_heap_merge, (void*)arg);
    }
    

    for (i = 0; i < asm_opt.thread_num; i++)
        pthread_join(_r_threads[i], NULL);
    free(_r_threads);

    ///rescue_edges(R_INF.paf, R_INF.reverse_paf, R_INF.total_reads, 4, 0.985);

    fprintf(stderr, "Final overlaps have been calculated.\n");
    fprintf(stderr, "%-30s%18.2f\n\n", "Final overlaps calculation time:", Get_T() - start_time); 

    Output_PAF();
    
    build_string_graph_without_clean(MIN_OVERLAP_COVERAGE, R_INF.paf, R_INF.reverse_paf, 
    R_INF.total_reads, R_INF.read_length, MIN_OVERLAP_LEN, MAX_HANG_LEN, asm_opt.clean_round, 
    asm_opt.pop_bubble_size, asm_opt.min_drop_rate, asm_opt.max_drop_rate, asm_opt.output_file_name, 
    MAX_BUBBLE_DIST, 0, 1);
}


void Correct_Reads(int last_round)
{
    
    if(asm_opt.load_index_from_disk && load_all_data_from_disk(&R_INF.paf, &R_INF.reverse_paf, 
    asm_opt.output_file_name))
    {
        build_string_graph_without_clean(MIN_OVERLAP_COVERAGE, R_INF.paf, R_INF.reverse_paf, 
        R_INF.total_reads, R_INF.read_length, MIN_OVERLAP_LEN, MAX_HANG_LEN, asm_opt.clean_round, 
        asm_opt.pop_bubble_size, asm_opt.min_drop_rate, asm_opt.max_drop_rate, 
        asm_opt.output_file_name, MAX_BUBBLE_DIST, 0, 0);
        exit(1);
    }
    else
    {
        ///fprintf(stderr, "Cannot find overlap file. Please run the whole hifiasm.\n");
    }
    
    clear_opt(&asm_opt, last_round);

    if(last_round == 0)
    {
        generate_overlaps(last_round);
        return;
    }
        
    fprintf(stderr, "Error correction: Start the %d-th round ...\n", asm_opt.roundID);
    
    Counting_multiple_thr();
    Build_hash_table_multiple_thr();
    Overlap_calculate_multipe_thr();

    fprintf(stderr, "Error correction: The %d-th round has been completed.\n", asm_opt.roundID);

    Correct_Reads(last_round - 1);
}












