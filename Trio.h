#ifndef __TRIO__
#define __TRIO__
#include <stdint.h>

#define YAK_MAX_KMER     31
#define YAK_COUNTER_BITS 10
#define YAK_N_COUNTS     (1<<YAK_COUNTER_BITS)
#define YAK_MAX_COUNT    ((1<<YAK_COUNTER_BITS)-1)

#define YAK_BLK_SHIFT  9 // 64 bytes, the size of a cache line
#define YAK_BLK_MASK   ((1<<(YAK_BLK_SHIFT)) - 1)

#define YAK_LOAD_ALL       1
#define YAK_LOAD_TRIOBIN1  2
#define YAK_LOAD_TRIOBIN2  3

#define YAK_MAGIC "YAK\2"

typedef struct {
	int n_shift, n_hashes;
	uint8_t *b;
} yak_bf_t;

struct yak_ht_t;

typedef struct {
	struct yak_ht_t *h;
	yak_bf_t *b;
} yak_ch1_t;

typedef struct {
	int k, pre, n_hash, n_shift;
	uint64_t tot;
	yak_ch1_t *h;
} yak_ch_t;


typedef struct {
	int max;
	uint32_t *s;
} tb_buf_t;


typedef struct {
	int k, n_threads, print_diff;
	double ratio_thres;
	///bseq_file_t *fp;
	const yak_ch_t *ch;
	tb_buf_t *buf;
	UC_Read *bseq;
	All_reads* seq; 
} tb_shared_t;

typedef struct {
	int c[16];
	int sc[2];
	int nk;
} tb_cnt_t;


typedef struct {
	int n_seq;
	tb_shared_t *aux;
	///bseq1_t *seq;
	///All_reads* seq;
	///tb_cnt_t *cnt;
} tb_step_t;


void trio_partition();
#endif