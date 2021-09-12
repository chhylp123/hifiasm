#ifndef __HA_HTAB_H__
#define __HA_HTAB_H__
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include "Process_Read.h"
#include "CommandLines.h"

typedef struct {
	int n, m;
	uint64_t *a;
} st_mt_t;

typedef struct {
	uint64_t x; ///x is the hash key
	///rid is the read id, pos is the end pos of this minimizer, rev is the direction
	///span is the length of this k-mer. For non-HPC k-mer, span may not be equal to k
	uint64_t rid:28, pos:27, rev:1, span:8;
} ha_mz1_t;

typedef struct {
	uint64_t rid:28, pos:27, rev:1, span:8; // actually it is not necessary to keep span in the index
} ha_idxpos_t;

typedef struct { uint32_t n, m; ha_mz1_t *a; } ha_mz1_v;

typedef struct {
	uint64_t x; ///x is the hash key
	uint64_t rid:31, rev:1, pos:32;
	uint8_t span;
} ha_mzl_t;

typedef struct {
	uint64_t rid:31, rev:1, pos:32;
	uint8_t span;
} ha_idxposl_t;

typedef struct { uint32_t n, m; ha_mzl_t *a; } ha_mzl_v;

typedef struct { // a simplified version of kdq
	int front, count;
	int a[64];
} tiny_queue_t;

static inline void tq_push(tiny_queue_t *q, int x)
{
	q->a[((q->count++) + q->front) & 0x3f] = x;
}

static inline int tq_shift(tiny_queue_t *q)
{
	int x;
	if (q->count == 0) return -1;
	x = q->a[q->front++];
	q->front &= 0x3f;
	--q->count;
	return x;
}

struct ha_pt_s;
typedef struct ha_pt_s ha_pt_t;

struct ha_abuf_s;
typedef struct ha_abuf_s ha_abuf_t;

struct ha_abufl_s;
typedef struct ha_abufl_s ha_abufl_t;

extern const unsigned char seq_nt4_table[256];
extern void *ha_flt_tab;
extern ha_pt_t *ha_idx;
extern void *ha_flt_tab_hp;
extern ha_pt_t *ha_idx_hp;
extern void *ha_ct_table;

void *ha_ft_ug_gen(const hifiasm_opt_t *asm_opt, ma_utg_v *us, int is_HPC, int k, int w, int min_freq, int max_freq);
void *ha_ft_gen(const hifiasm_opt_t *asm_opt, All_reads *rs, int *hom_cov, int is_hp_mode);
int32_t ha_ft_cnt(const void *hh, uint64_t y);
void ha_ft_destroy(void *h);

ha_pt_t *ha_pt_ug_gen(const hifiasm_opt_t *asm_opt, const void *flt_tab, ma_utg_v *us, int is_HPC, int k, int w, int min_freq);
ha_pt_t *ha_pt_gen(const hifiasm_opt_t *asm_opt, const void *flt_tab, int read_from_store, int is_hp_mode, All_reads *rs, int *hom_cov, int *het_cov);
void ha_pt_destroy(ha_pt_t *h);
const ha_idxpos_t *ha_pt_get(const ha_pt_t *h, uint64_t hash, int *n);
const ha_idxposl_t *ha_ptl_get(const ha_pt_t *h, uint64_t hash, int *n);
const int ha_pt_cnt(const ha_pt_t *h, uint64_t hash);

int write_pt_index(void *flt_tab, ha_pt_t *ha_idx, All_reads* r, hifiasm_opt_t* opt, char* file_name);
int load_pt_index(void **r_flt_tab, ha_pt_t **r_ha_idx, All_reads* r, hifiasm_opt_t* opt, char* file_name);
int write_ct_index(void *ct_idx, char* file_name);
int load_ct_index(void **ct_idx, char* file_name);
int query_ct_index(void* ct_idx, uint64_t hash);

ha_abuf_t *ha_abuf_init(void);
void ha_abuf_destroy(ha_abuf_t *ab);
uint64_t ha_abuf_mem(const ha_abuf_t *ab);
ha_abufl_t *ha_abufl_init(void);
void ha_abufl_destroy(ha_abufl_t *ab);
uint64_t ha_abufl_mem(const ha_abufl_t *ab);

double yak_cputime(void);
void yak_reset_realtime(void);
double yak_realtime(void);
long yak_peakrss(void);
double yak_peakrss_in_gb(void);
double yak_cpu_usage(void);

void ha_triobin(const hifiasm_opt_t *opt);

void mz1_ha_sketch(const char *str, int len, int w, int k, uint32_t rid, int is_hpc, ha_mz1_v *p, const void *hf, int sample_dist, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, ha_pt_t *pt, int min_freq, int32_t dp_min_len, float dp_e, st_mt_t *mt, int32_t ws, int32_t is_unique);
void mz2_ha_sketch(const char *str, int len, int w, int k, uint32_t rid, int is_hpc, ha_mzl_v *p, const void *hf, int sample_dist, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, ha_pt_t *pt, int min_freq, int32_t dp_min_len, float dp_e, st_mt_t *mt, int32_t ws, int32_t is_unique);
int ha_analyze_count(int n_cnt, int start_cnt, int m_peak_hom, const int64_t *cnt, int *peak_het);
int adj_m_peak_hom(int m_peak_hom, int max_i, int max2_i, int max3_i, int *peak_het);
void print_hist_lines(int n_cnt, int start_cnt, const int64_t *cnt);
void debug_adapter(const hifiasm_opt_t *asm_opt, All_reads *rs);

inline int mz_low_b(int peak_hom, int peak_het)
{
	int low_freq = 2;
	if(peak_het > 0) low_freq = peak_het/2;
	else if(peak_hom > 0) low_freq = peak_hom/4;
	if(low_freq < 2) low_freq = 2;
	return low_freq;
}

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
	///compare forward k-mer and reverse complementary strand
	int j = x[1] < x[3]? 0 : 1;
	return yak_hash64_64(x[j<<1|0]) + yak_hash64_64(x[j<<1|1]);
}

#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#ifndef kroundup64
#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, x|=(x)>>32, ++(x))
#endif

#ifndef klib_unused
#if (defined __clang__ && __clang_major__ >= 3) || (defined __GNUC__ && __GNUC__ >= 3)
#define klib_unused __attribute__ ((__unused__))
#else
#define klib_unused
#endif
#endif /* klib_unused */

#endif // __YAK_H__
