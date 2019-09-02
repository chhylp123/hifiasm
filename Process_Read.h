#ifndef __READ__
#define __READ__

#include<stdint.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
#include "kseq.h"

#define READ_INIT_NUMBER 1000

#define READ_BLOCK_SIZE 64
#define READ_BLOCK_NUM_PRE_THR 100

#define IS_FULL(buffer) ((buffer.num >= buffer.size)?1:0)
#define IS_EMPTY(buffer) ((buffer.num == 0)?1:0)
///#define Get_READ_LENGTH(R_INF, ID) (R_INF.index[ID+1] - R_INF.index[ID])
#define Get_READ_LENGTH(R_INF, ID) R_INF.read_length[ID]
#define Get_NAME_LENGTH(R_INF, ID) (R_INF.name_index[ID+1] - R_INF.name_index[ID])
///#define Get_READ(R_INF, ID) R_INF.read + (R_INF.index[ID]>>2) + ID
#define Get_READ(R_INF, ID) R_INF.read_sperate[ID]
#define Get_NAME(R_INF, ID) R_INF.name + R_INF.name_index[ID]



KSEQ_INIT(gzFile, gzread)



static uint8_t seq_nt6_table[256] = {
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 0, 5, 1,  5, 5, 5, 2,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  3, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 0, 5, 1,  5, 5, 5, 2,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  3, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

static char bit_t_seq_table[256][4] = {0};
static char bit_t_seq_table_rc[256][4] = {0};
static char s_H[5] = {'A', 'C', 'G', 'T', 'N'};
static char rc_Table[5] = {'T', 'G', 'C', 'A', 'N'};

#define RC_CHAR(x) rc_Table[seq_nt6_table[(uint8_t)x]]


void init_kseq(char* file);
void destory_kseq();
int get_read(kseq_t *s);

typedef struct
{
    /**[0-1] bits are type:**/
    /**[2-31] bits are length**/
    uint32_t* record;
    uint32_t length;
	uint32_t size;

    char* lost_base;
    uint32_t lost_base_length;
	uint32_t lost_base_size;
	uint32_t new_length;
}Compressed_Cigar_record;


typedef struct
{
	uint64_t** N_site;
	///uint8_t* read;
	char* name;


	uint8_t** read_sperate;
	uint64_t* read_length;
	uint64_t* read_size;

	///seq start pos in uint8_t* read
	///do not need it
	///uint64_t* index;
	uint64_t index_size;


    ///name start pos in char* name
	uint64_t* name_index;
	uint64_t name_index_size;
	uint64_t total_reads;
	uint64_t total_reads_bases;
	uint64_t total_name_length;

	Compressed_Cigar_record* cigars; 

} All_reads;

extern All_reads R_INF;

void malloc_All_reads(All_reads* r);

typedef struct
{
	kseq_t* read;
	long long num;

} R_buffer_block;


typedef struct
{
	R_buffer_block* sub_block;
	long long block_inner_size;
	long long size;
	long long num;
	int all_read_end;
} R_buffer;


typedef struct
{
	char* seq;
	long long length;
	long long size;
	long long RID;
} UC_Read;

void init_R_buffer(int thread_num);
void init_All_reads(All_reads* r);
void* input_reads_muti_threads(void*);
void init_R_buffer_block(R_buffer_block* curr_sub_block);
int get_reads_mul_thread(R_buffer_block* curr_sub_block);
void compress_base(uint8_t* dest, char* src, uint64_t src_l, uint64_t** N_site_lis, uint64_t N_site_occ);
void init_UC_Read(UC_Read* r);
void recover_UC_Read(UC_Read* r, All_reads* R_INF, uint64_t ID);
void recover_UC_Read_RC(UC_Read* r, All_reads* R_INF, uint64_t ID);
void recover_UC_Read_sub_region(char* r, long long start_pos, long long length, uint8_t strand, All_reads* R_INF, long long ID);
void destory_UC_Read(UC_Read* r);
void reverse_complement(char* pattern, uint64_t length);
void write_All_reads(All_reads* r, char* read_file_name);
int load_All_reads(All_reads* r, char* read_file_name);
void destory_All_reads(All_reads* r);

void Counting_block();
void destory_R_buffer_block(R_buffer_block* curr_sub_block);
void destory_R_buffer();
void clear_R_buffer();


#endif
