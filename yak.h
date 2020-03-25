#ifndef __YAK_H__
#define __YAK_H__
#define __STDC_LIMIT_MACROS
#include <stdint.h>

#define YAK_BLK_SHIFT  9 // 64 bytes, the size of a cache line
#define YAK_BLK_MASK   ((1<<(YAK_BLK_SHIFT)) - 1)

typedef struct {
	uint64_t x;
	uint64_t rid:28, pos:27, rev:1, span:8;
} ha_mz1_t;

typedef struct { uint32_t n, m; ha_mz1_t *a; } ha_mz1_v;

typedef struct {
	int n_shift, n_hashes;
	uint8_t *b;
} yak_bf_t;

extern const unsigned char seq_nt4_table[256];

int ha_hf_isflt(const void *hh, uint64_t y);
void ha_hf_destroy(void *h);

void trio_partition(void);

double yak_cputime(void);
void yak_reset_realtime(void);
double yak_realtime(void);
long yak_peakrss(void);

void ha_sketch(const char *str, int len, int w, int k, uint32_t rid, int is_hpc, ha_mz1_v *p, const void *hf);
int yak_analyze_count(int n_cnt, const int64_t *cnt, int *peak_het);

yak_bf_t *yak_bf_init(int n_shift, int n_hashes);
void yak_bf_destroy(yak_bf_t *b);
int yak_bf_insert(yak_bf_t *b, uint64_t hash);

static inline uint64_t yak_hash64(uint64_t key, uint64_t mask) // invertible integer hash function
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

static inline uint64_t yak_hash64_64(uint64_t key)
{
	key = ~key + (key << 21);
	key = key ^ key >> 24;
	key = (key + (key << 3)) + (key << 8);
	key = key ^ key >> 14;
	key = (key + (key << 2)) + (key << 4);
	key = key ^ key >> 28;
	key = key + (key << 31);
	return key;
}

static inline uint64_t yak_hash_long(uint64_t x[4])
{
	int j = x[1] < x[3]? 0 : 1;
	return yak_hash64_64(x[j<<1|0]) + yak_hash64_64(x[j<<1|1]);
}

#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

#ifndef klib_unused
#if (defined __clang__ && __clang_major__ >= 3) || (defined __GNUC__ && __GNUC__ >= 3)
#define klib_unused __attribute__ ((__unused__))
#else
#define klib_unused
#endif
#endif /* klib_unused */

#endif // __YAK_H__
