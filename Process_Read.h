#ifndef __READ__
#define __READ__

#include<stdint.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
#include "Overlaps.h"
#include "CommandLines.h"
///#include "Hash_Table.h"

#define READ_INIT_NUMBER 1000

#define READ_BLOCK_SIZE 64
#define READ_BLOCK_NUM_PRE_THR 100

#define IS_FULL(buffer) ((buffer.num >= buffer.size)?1:0)
#define IS_EMPTY(buffer) ((buffer.num == 0)?1:0)
///#define Get_READ_LENGTH(R_INF, ID) (R_INF.index[ID+1] - R_INF.index[ID])
#define Get_READ_LENGTH(R_INF, ID) (R_INF).read_length[(ID)]
#define Get_NAME_LENGTH(R_INF, ID) ((R_INF).name_index[(ID)+1] - (R_INF).name_index[(ID)])
///#define Get_READ(R_INF, ID) R_INF.read + (R_INF.index[ID]>>2) + ID
#define Get_READ(R_INF, ID) (R_INF).read_sperate[(ID)]
#define Get_NAME(R_INF, ID) ((R_INF).name + (R_INF).name_index[(ID)])
#define CHECK_BY_NAME(R_INF, NAME, ID) (Get_NAME_LENGTH((R_INF),(ID))==strlen((NAME)) && \
                                        memcmp((NAME), Get_NAME((R_INF), (ID)), Get_NAME_LENGTH((R_INF),(ID))) == 0)

extern uint8_t seq_nt6_table[256];
extern char bit_t_seq_table[256][4];
extern char bit_t_seq_table_rc[256][4];
extern char s_H[5];
extern char rc_Table[5];


#define RC_CHAR(x) rc_Table[seq_nt6_table[(uint8_t)x]]

void init_aux_table();

typedef struct
{
    uint64_t x_id;
    uint64_t x_pos_s;
    uint64_t x_pos_e;
    uint8_t x_pos_strand;

    uint64_t y_id;
    uint64_t y_pos_s;
    uint64_t y_pos_e;
    uint8_t y_pos_strand;

	uint64_t matchLen;
	uint64_t totalLen;

} PAF;

typedef struct
{
    PAF* list;
    uint64_t size;
    uint64_t length;
} PAF_alloc;


inline void init_PAF_alloc(PAF_alloc* list)
{
	list->size = 15;
    list->length = 0;
    list->list = (PAF*)malloc(sizeof(PAF)*list->size);
}

inline void append_PAF_alloc(PAF_alloc* list, PAF* e)
{
    if(list->length+1 > list->size)
    {
        list->size = list->size * 2;
        list->list = (PAF*)realloc(list->list, sizeof(PAF)*list->size);
    }

    list->list[list->length] = (*e);
    list->length++;
}

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
} Compressed_Cigar_record;


#define AMBIGU 0
#define FATHER 1
#define MOTHER 2
#define MIX_TRIO 3
#define NON_TRIO 4
#define DROP 5
#define SET_TRIO 8
#define CHAIN_MATCH 1
#define CHAIN_UNMATCH 0.334

typedef struct
{
	uint64_t** N_site;
	///uint8_t* read;
	char* name;

	uint8_t** read_sperate;
	uint64_t* read_length;
	uint64_t* read_size;
	uint8_t* trio_flag;

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
	Compressed_Cigar_record* second_round_cigar;

    ma_hit_t_alloc* paf;
    ma_hit_t_alloc* reverse_paf;

    ///kvec_t_u64_warp* pb_regions;
} All_reads;

extern All_reads R_INF;

typedef struct
{
	char* seq;
	long long length;
	long long size;
	long long RID;
} UC_Read;

typedef struct
{
    char** read_name;
    uint64_t query_num;
    kvec_t_u64_warp* candidate_count;
    FILE* fp;
    pthread_mutex_t OutputMutex;
} Debug_reads;

void init_All_reads(All_reads* r);
void malloc_All_reads(All_reads* r);
void ha_insert_read_len(All_reads *r, int read_len, int name_len);
void ha_compress_base(uint8_t* dest, char* src, uint64_t src_l, uint64_t** N_site_lis, uint64_t N_site_occ);
void init_UC_Read(UC_Read* r);
void recover_UC_Read(UC_Read* r, const All_reads *R_INF, uint64_t ID);
void recover_UC_Read_RC(UC_Read* r, All_reads* R_INF, uint64_t ID);
void recover_UC_Read_sub_region(char* r, long long start_pos, long long length, uint8_t strand, All_reads* R_INF, long long ID);
void destory_UC_Read(UC_Read* r);
void reverse_complement(char* pattern, uint64_t length);
void write_All_reads(All_reads* r, char* read_file_name);
int load_All_reads(All_reads* r, char* read_file_name);
void destory_All_reads(All_reads* r);
int destory_read_bin(All_reads* r);
void init_Debug_reads(Debug_reads* x, const char* file);
void destory_Debug_reads(Debug_reads* x);
void recover_UC_sub_Read(UC_Read* i_r, long long start_pos, long long length, uint8_t strand, All_reads* R_INF, long long ID);
#endif
