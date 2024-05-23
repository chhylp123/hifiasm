#ifndef __READ__
#define __READ__

#define __STDC_LIMIT_MACROS
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
#define IS_SCAF_READ(R_INF, ID) ((R_INF).read_sperate[(ID)] == NULL)

extern uint8_t seq_nt6_table[256];
extern char bit_t_seq_table[256][4];
extern char bit_t_seq_table_rc[256][4];
extern char s_H[5];
extern char rc_Table[6];


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
    uint64_t *read_id;
    uint64_t query_num;
    kvec_t_u64_warp* candidate_count;
    FILE *fp, *fp_r0, *fp_r1;
    pthread_mutex_t OutputMutex;
} Debug_reads;


typedef struct
{
    uint32_t hid;
    uint32_t qs, qe, ts, te; uint32_t pidx, pdis, aidx;///TODO: enable pdis
    uint8_t pchain:5, rev:1, base:1, el:1;
} uc_block_t;

typedef struct
{
    uint32_t *a;
    size_t n, m;
} N_t;

typedef struct
{
    char *a; uint32_t n;
} nid_t;

typedef struct
{
    kvec_t(uint8_t) r_base;
    uint32_t rlen; 

    kvec_t(uc_block_t) bb;
    N_t N_site;

    uint8_t dd;
} ul_vec_t;

typedef struct{
    kvec_t(uint32_t) idx;
    kvec_t(uint64_t) occ;
} ul_vec_rid_t;

typedef struct
{
    kvec_t(nid_t) nid;
    ul_vec_rid_t ridx;
    ul_vec_t *a;
    size_t n, m;
    All_reads *hR;
    // idx_emask_t *mm;
    // uint32_t mm;
} all_ul_t;


typedef struct {
    ul_vec_t *a;
    size_t n, m;
    uint8_t dd;
} scaf_res_t;

extern all_ul_t UL_INF;
extern all_ul_t ULG_INF;
// extern uint32_t *het_cnt;
// extern uint32_t debug_out;

void init_All_reads(All_reads* r);
void malloc_All_reads(All_reads* r);
void ha_insert_read_len(All_reads *r, int read_len, int name_len);
void ha_compress_base(uint8_t* dest, char* src, uint64_t src_l, uint64_t** N_site_lis, uint64_t N_site_occ);
void init_UC_Read(UC_Read* r);
void recover_UC_Read(UC_Read* r, const All_reads *R_INF, uint64_t ID);
void recover_UC_Read_RC(UC_Read* r, All_reads* R_INF, uint64_t ID);
void recover_UC_Read_sub_region(char* r, int64_t start_pos, int64_t length, uint8_t strand, All_reads* R_INF, int64_t ID);
void destory_UC_Read(UC_Read* r);
void reverse_complement(char* pattern, uint64_t length);
void write_All_reads(All_reads* r, char* read_file_name);
int load_All_reads(All_reads* r, char* read_file_name);
int append_All_reads(All_reads* r, char *idx, uint32_t id);
void destory_All_reads(All_reads* r);
int destory_read_bin(All_reads* r);
void init_Debug_reads(Debug_reads* x, const char* file);
void destory_Debug_reads(Debug_reads* x);
void recover_UC_sub_Read(UC_Read* i_r, long long start_pos, long long length, uint8_t strand, All_reads* R_INF, long long ID);

void init_all_ul_t(all_ul_t *x, All_reads *hR);
void destory_all_ul_t(all_ul_t *x);
void append_ul_t(all_ul_t *x, uint64_t *rid, char* id, int64_t id_l, char* str, int64_t str_l, ul_ov_t *o, int64_t on, float p_chain_rate, const ug_opt_t *uopt, uint32_t save_bases);
void retrieve_ul_t(UC_Read* i_r, char *i_s, all_ul_t *ref, uint64_t ID, uint8_t strand, int64_t s, int64_t l);
void retrieve_u_seq(UC_Read* i_r, char* i_s, ma_utg_t *u, uint8_t strand, int64_t s, int64_t l, void *km);
void debug_retrieve_rc_sub(const ug_opt_t *uopt, all_ul_t *ref, const All_reads *R_INF, ul_idx_t *ul, uint32_t n_step);
uint32_t retrieve_u_cov(const ul_idx_t *ul, uint64_t id, uint8_t strand, uint64_t pos, uint8_t dir, int64_t *pi);
uint64_t retrieve_u_cov_region(const ul_idx_t *ul, uint64_t id, uint8_t strand, uint64_t s, uint64_t e, int64_t *pi);
uint64_t retrieve_r_cov_region(const ul_idx_t *ul, uint64_t id, uint8_t strand, uint64_t s, uint64_t e, int64_t *pi);
void append_ul_t_back(all_ul_t *x, uint64_t *rid, char* id, int64_t id_l, char* str, int64_t str_l, ul_ov_t *o, int64_t on, float p_chain_rate);
void write_compress_base_disk(FILE *fp, uint64_t ul_rid, char *str, uint32_t len, ul_vec_t *buf);
int64_t load_compress_base_disk(FILE *fp, uint64_t *ul_rid, char *dest, uint32_t *len, ul_vec_t *buf);
scaf_res_t *init_scaf_res_t(uint32_t n);
void destroy_scaf_res_t(scaf_res_t *p);
void read_ma(ma_hit_t* x, FILE* fp);

#endif
