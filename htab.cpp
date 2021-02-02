#include <stdint.h>
#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "kthread.h"
#include "khashl.h"
#include "kseq.h"
#include "ksort.h"
#include "htab.h"

#define YAK_COUNTER_BITS 12
#define YAK_N_COUNTS     (1<<YAK_COUNTER_BITS)
#define YAK_MAX_COUNT    ((1<<YAK_COUNTER_BITS)-1)

const unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

void *ha_flt_tab;
ha_pt_t *ha_idx;
void *ha_flt_tab_hp;
ha_pt_t *ha_idx_hp;
void *ha_ct_table;

/***************************
 * Yak specific parameters *
 ***************************/

typedef struct {
	int32_t bf_shift, bf_n_hash;
	int32_t k, w, is_HPC;
	int32_t pre;
	int32_t n_thread;
	int64_t chunk_size;
	int adaLen;
} yak_copt_t;

void yak_copt_init(yak_copt_t *o)
{
	memset(o, 0, sizeof(yak_copt_t));
	o->bf_shift = 0;
	o->bf_n_hash = 4;
	o->k = 31;
	o->w = 1;
	o->pre = YAK_COUNTER_BITS;
	o->n_thread = 4;
	o->chunk_size = 20000000;
}

/************************
 * Blocked bloom filter *
 ************************/

#define YAK_BLK_SHIFT  9 // 64 bytes, the size of a cache line
#define YAK_BLK_MASK   ((1<<(YAK_BLK_SHIFT)) - 1)

typedef struct {
	int n_shift, n_hashes;
	uint8_t *b;
} yak_bf_t;
///in most cases, n_shift = 25, n_hashes = 4
yak_bf_t *yak_bf_init(int n_shift, int n_hashes)
{
	yak_bf_t *b;
	void *ptr = 0;
	if (n_shift + YAK_BLK_SHIFT > 64 || n_shift < YAK_BLK_SHIFT) return 0;
	CALLOC(b, 1);
	b->n_shift = n_shift;
	b->n_hashes = n_hashes;
	posix_memalign(&ptr, 1<<(YAK_BLK_SHIFT-3), 1ULL<<(n_shift-3));
	b->b = (uint8_t*)ptr;
	bzero(b->b, 1ULL<<(n_shift-3));
	return b;
}

void yak_bf_destroy(yak_bf_t *b)
{
	if (b == 0) return;
	free(b->b); free(b);
}

int yak_bf_insert(yak_bf_t *b, uint64_t hash)
{
	int x = b->n_shift - YAK_BLK_SHIFT;
	uint64_t y = hash & ((1ULL<<x) - 1);
	int h1 = hash >> x & YAK_BLK_MASK;
	int h2 = hash >> b->n_shift & YAK_BLK_MASK;
	uint8_t *p = &b->b[y<<(YAK_BLK_SHIFT-3)];
	int i, z = h1, cnt = 0;
	if ((h2&31) == 0) h2 = (h2 + 1) & YAK_BLK_MASK; // otherwise we may repeatedly use a few bits
	for (i = 0; i < b->n_hashes; z = (z + h2) & YAK_BLK_MASK) {
		uint8_t *q = &p[z>>3], u;
		u = 1<<(z&7);
		cnt += !!(*q & u);
		*q |= u;
		++i;
	}
	return cnt;
}

/********************
 * Count hash table *
 ********************/

#define yak_ct_eq(a, b) ((a)>>YAK_COUNTER_BITS == (b)>>YAK_COUNTER_BITS) // lower 8 bits for counts; higher bits for k-mer
#define yak_ct_hash(a) ((a)>>YAK_COUNTER_BITS)
KHASHL_SET_INIT(static klib_unused, yak_ct_t, yak_ct, uint64_t, yak_ct_hash, yak_ct_eq)

typedef struct {
	yak_ct_t *h;
	yak_bf_t *b;
} ha_ct1_t;

typedef struct {
	int k, pre, n_hash, n_shift;
	uint64_t tot;  ///number of distinct k-mers
	ha_ct1_t *h;
} ha_ct_t;

///for 0-th counting, k = 51, pre = 12, n_hash = 4, n_shift = 37
///for 1-th counting, opt.k = 51, opt->pre = 12, opt->bf_n_hash = 4, opt.bf_shift = 0
static ha_ct_t *ha_ct_init(int k, int pre, int n_hash, int n_shift)
{
	ha_ct_t *h;
	int i;
	if (pre < YAK_COUNTER_BITS) return 0;
	CALLOC(h, 1);
	h->k = k, h->pre = pre;
	CALLOC(h->h, 1<<h->pre);
	///i<h->pre = 4096
	///it seems there is a large hash table h, consisting 4096 small hash tables
	for (i = 0; i < 1<<h->pre; ++i)
		h->h[i].h = yak_ct_init();
	///for 0-th counting, enter here; used for bloom filter
	if (n_hash > 0 && n_shift > h->pre) {
		h->n_hash = n_hash, h->n_shift = n_shift;
		for (i = 0; i < 1<<h->pre; ++i)
			h->h[i].b = yak_bf_init(h->n_shift - h->pre, h->n_hash); ///h->n_shift = 37, h->pre = 12, h->n_hash = 4 
	}
	return h;
}

static void ha_ct_destroy_bf(ha_ct_t *h)
{
	int i;
	for (i = 0; i < 1<<h->pre; ++i) {
		if (h->h[i].b)
			yak_bf_destroy(h->h[i].b);
		h->h[i].b = 0;
	}
}

static void ha_ct_destroy(ha_ct_t *h)
{
	int i;
	if (h == 0) return;
	ha_ct_destroy_bf(h);
	for (i = 0; i < 1<<h->pre; ++i)
		yak_ct_destroy(h->h[i].h);
	free(h->h); free(h);
}

static int ha_ct_insert_list(ha_ct_t *h, int create_new, int n, const uint64_t *a)
{
	int j, mask = (1<<h->pre) - 1, n_ins = 0;
	ha_ct1_t *g;
	if (n == 0) return 0;
	///corresponding small hash index
	g = &h->h[a[0]&mask];
	for (j = 0; j < n; ++j) {
		int ins = 1, absent;
		///x is a 64-bit word, h->pre=12
		///all elements at a have the same low 12 bits
		///so low 12 bits are not useful
		uint64_t x = a[j] >> h->pre;
		khint_t k;
		if ((a[j]&mask) != (a[0]&mask)) continue;
		if (create_new) {
			///for 0-th counting, g->b = NULL
			if (g->b)
				ins = (yak_bf_insert(g->b, x) == h->n_hash);
			///for 0-th counting, g->b = NULL
			///x = the high 52 bits of a[j] + low 12 bits 0
			///the low 12 bits are used for counting
			if (ins) {
				k = yak_ct_put(g->h, x << YAK_COUNTER_BITS | (g->b? 1 : 0), &absent);
				if (absent) ++n_ins;
				if ((kh_key(g->h, k)&YAK_MAX_COUNT) < YAK_MAX_COUNT)
					++kh_key(g->h, k);
			}
		} else {
			k = yak_ct_get(g->h, x<<YAK_COUNTER_BITS);
			if (k != kh_end(g->h) && (kh_key(g->h, k)&YAK_MAX_COUNT) < YAK_MAX_COUNT)
				++kh_key(g->h, k);
		}
	}
	return n_ins;
}

/*** generate histogram ***/

typedef struct {
	uint64_t c[YAK_N_COUNTS];
} buf_cnt_t;

typedef struct {
	const ha_ct_t *h;
	buf_cnt_t *cnt;
} hist_aux_t;

static void worker_ct_hist(void *data, long i, int tid) // callback for kt_for()
{
	hist_aux_t *a = (hist_aux_t*)data;
	uint64_t *cnt = a->cnt[tid].c;
	yak_ct_t *g = a->h->h[i].h;
	khint_t k;
	for (k = 0; k < kh_end(g); ++k)
		if (kh_exist(g, k))
			++cnt[kh_key(g, k)&YAK_MAX_COUNT];
}

///YAK_N_COUNTS is also 4096
///used for calculating k-mer histogram
static void ha_ct_hist(const ha_ct_t *h, int64_t cnt[YAK_N_COUNTS], int n_thread)
{
	hist_aux_t a;
	int i, j;
	a.h = h;
	memset(cnt, 0, YAK_N_COUNTS * sizeof(uint64_t));
	CALLOC(a.cnt, n_thread);
	///start 4096 threads
	kt_for(n_thread, worker_ct_hist, &a, 1<<h->pre);
	for (i = 0; i < YAK_N_COUNTS; ++i) cnt[i] = 0;
	for (j = 0; j < n_thread; ++j)
		for (i = 0; i < YAK_N_COUNTS; ++i)
			cnt[i] += a.cnt[j].c[i];
	free(a.cnt);
}

/*** shrink a hash table ***/

typedef struct {
	int min, max;
	ha_ct_t *h;
} shrink_aux_t;

static void worker_ct_shrink(void *data, long i, int tid) // callback for kt_for()
{
	shrink_aux_t *a = (shrink_aux_t*)data;
	ha_ct_t *h = a->h;
	yak_ct_t *g = h->h[i].h, *f;
	khint_t k;
	f = yak_ct_init();
	yak_ct_resize(f, kh_size(g));
	for (k = 0; k < kh_end(g); ++k) {
		if (kh_exist(g, k)) {
			int absent, c = kh_key(g, k) & YAK_MAX_COUNT;
			if (c >= a->min && c <= a->max)
				yak_ct_put(f, kh_key(g, k), &absent);
		}
	}
	yak_ct_destroy(g);
	h->h[i].h = f;
}

static void ha_ct_shrink(ha_ct_t *h, int min, int max, int n_thread)
{
	int i;
	shrink_aux_t a;
	a.h = h, a.min = min, a.max = max;
	///still start 4096 threads
	kt_for(n_thread, worker_ct_shrink, &a, 1<<h->pre);
	for (i = 0, h->tot = 0; i < 1<<h->pre; ++i)
		h->tot += kh_size(h->h[i].h);
}

/***********************
 * Position hash table *
 ***********************/

KHASHL_MAP_INIT(static klib_unused, yak_pt_t, yak_pt, uint64_t, uint64_t, yak_ct_hash, yak_ct_eq)
#define generic_key(x) (x)
KRADIX_SORT_INIT(ha64, uint64_t, generic_key, 8)

typedef struct {
	yak_pt_t *h;
	uint64_t n;
	ha_idxpos_t *a;
} ha_pt1_t;

struct ha_pt_s {
	int k, pre;
	uint64_t tot, tot_pos;
	ha_pt1_t *h;
};

typedef struct {
	const ha_ct_t *ct;
	ha_pt_t *pt;
} pt_gen_aux_t;

static void worker_pt_gen(void *data, long i, int tid) // callback for kt_for()
{
	pt_gen_aux_t *a = (pt_gen_aux_t*)data;
	ha_pt1_t *b = &a->pt->h[i];
	yak_ct_t *g = a->ct->h[i].h;
	khint_t k;
	for (k = 0, b->n = 0; k != kh_end(g); ++k) {
		if (kh_exist(g, k)) {
			int absent;
			khint_t l;
			l = yak_pt_put(b->h, kh_key(g, k) >> a->ct->pre << YAK_COUNTER_BITS, &absent);
			///this should be the start index of kh_key's corresponding pos at ha_idxpos_t* a
			kh_val(b->h, l) = b->n;
			b->n += kh_key(g, k) & YAK_MAX_COUNT;
		}
	}
	yak_ct_destroy(g);
	a->ct->h[i].h = 0;
	CALLOC(b->a, b->n);
}

ha_pt_t *ha_pt_gen(ha_ct_t *ct, int n_thread)
{
	pt_gen_aux_t a;
	int i;
	ha_pt_t *pt;
	ha_ct_destroy_bf(ct);
	CALLOC(pt, 1);
	pt->k = ct->k, pt->pre = ct->pre, pt->tot = ct->tot;
	CALLOC(pt->h, 1<<pt->pre);
	for (i = 0; i < 1<<pt->pre; ++i) {
		pt->h[i].h = yak_pt_init();
		yak_pt_resize(pt->h[i].h, kh_size(ct->h[i].h));
	}
	a.ct = ct, a.pt = pt;
	kt_for(n_thread, worker_pt_gen, &a, 1<<pt->pre);
	free(ct->h); free(ct);
	return pt;
}

int ha_pt_insert_list(ha_pt_t *h, int n, const ha_mz1_t *a)
{
	int j, mask = (1<<h->pre) - 1, n_ins = 0;
	ha_pt1_t *g;
	if (n == 0) return 0;
	g = &h->h[a[0].x&mask];
	for (j = 0; j < n; ++j) {
		uint64_t x = a[j].x >> h->pre;
		khint_t k;
		int n;
		ha_idxpos_t *p;
		if ((a[j].x&mask) != (a[0].x&mask)) continue;
		k = yak_pt_get(g->h, x<<YAK_COUNTER_BITS);
		if (k == kh_end(g->h)) continue;
		n = kh_key(g->h, k) & YAK_MAX_COUNT;
		assert(n < YAK_MAX_COUNT);
		p = &g->a[kh_val(g->h, k) + n];
		p->rid = a[j].rid, p->rev = a[j].rev, p->pos = a[j].pos, p->span = a[j].span;
		//(uint64_t)a[j].rid<<36 | (uint64_t)a[j].rev<<35 | (uint64_t)a[j].pos<<8 | (uint64_t)a[j].span;
		++kh_key(g->h, k);
		++n_ins;
	}
	return n_ins;
}
/*
static void worker_pt_sort(void *data, long i, int tid)
{
	ha_pt_t *h = (ha_pt_t*)data;
	ha_pt1_t *g = &h->h[i];
	khint_t k;
	for (k = 0; k < kh_end(g->h); ++k) {
		int n;
		uint64_t *p;
		if (!kh_exist(g->h, k)) continue;
		n = kh_key(g->h, k) & YAK_MAX_COUNT;
		p = &g->a[kh_val(g->h, k)];
		radix_sort_ha64(p, p + n);
	}
}

void ha_pt_sort(ha_pt_t *h, int n_thread)
{
	kt_for(n_thread, worker_pt_sort, h, 1<<h->pre);
}
*/
void ha_pt_destroy(ha_pt_t *h)
{
	int i;
	if (h == 0) return;
	for (i = 0; i < 1<<h->pre; ++i) {
		yak_pt_destroy(h->h[i].h);
		free(h->h[i].a);
	}
	free(h->h); free(h);
}

const ha_idxpos_t *ha_pt_get(const ha_pt_t *h, uint64_t hash, int *n)
{
	khint_t k;
	const ha_pt1_t *g = &h->h[hash & ((1ULL<<h->pre) - 1)];
	*n = 0;
	k = yak_pt_get(g->h, hash >> h->pre << YAK_COUNTER_BITS);
	if (k == kh_end(g->h)) return 0;
	*n = kh_key(g->h, k) & YAK_MAX_COUNT;
	return &g->a[kh_val(g->h, k)];
}

/**********************************
 * Buffer for counting all k-mers *
 **********************************/

typedef struct {
	int n, m;
	uint64_t n_ins;
	uint64_t *a;
	ha_mz1_t *b;
} ch_buf_t;

///p = 12
static inline void ct_insert_buf(ch_buf_t *buf, int p, uint64_t y) // insert a k-mer $y to a linear buffer
{
	///assign k-mer to one of the 4096 bins
	///using low 12 bits for assigning
	///so all elements at b have the same low 12 bits
	int pre = y & ((1<<p) - 1);
	ch_buf_t *b = &buf[pre];
	if (b->n == b->m) {
		b->m = b->m < 8? 8 : b->m + (b->m>>1);
		REALLOC(b->a, b->m);
	}
	b->a[b->n++] = y;
}

static inline void pt_insert_buf(ch_buf_t *buf, int p, const ha_mz1_t *y)
{
	///assign minimizer to one of 4096 bins by low 12 bits
	int pre = y->x & ((1<<p) - 1);
	ch_buf_t *b = &buf[pre];
	if (b->n == b->m) {
		b->m = b->m < 8? 8 : b->m + (b->m>>1);
		REALLOC(b->b, b->m);
	}
	b->b[b->n++] = *y;
}

///buf is the read block, k is the k-mer length, p = 12, len is the read length, seq is the read
static void count_seq_buf(ch_buf_t *buf, int k, int p, int len, const char *seq) // insert k-mers in $seq to linear buffer $buf
{
	int i, l;
	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		///c = 00, 01, 10, 11
		if (c < 4) { // not an "N" base
			///x[0] & x[1] are the forward k-mer
			///x[2] & x[3] are the reverse complementary k-mer
			x[0] = (x[0] << 1 | (c&1))  & mask;
			x[1] = (x[1] << 1 | (c>>1)) & mask;
			x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
			x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
			if (++l >= k)
				ct_insert_buf(buf, p, yak_hash_long(x));
		} else l = 0, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
	}
}

static void count_seq_buf_HPC(ch_buf_t *buf, int k, int p, int len, const char *seq) // insert k-mers in $seq to linear buffer $buf
{
	int i, l, last = -1;
	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			if (c != last) {
				x[0] = (x[0] << 1 | (c&1))  & mask;
				x[1] = (x[1] << 1 | (c>>1)) & mask;
				x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
				x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
				if (++l >= k)
					ct_insert_buf(buf, p, yak_hash_long(x));
				last = c;
			}
		} else l = 0, last = -1, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
	}
}

/******************
 * K-mer counting *
 ******************/

KSEQ_INIT(gzFile, gzread)

#define HAF_COUNT_EXACT  0x1
#define HAF_COUNT_ALL    0x2
#define HAF_RS_WRITE_LEN 0x4
#define HAF_RS_WRITE_SEQ 0x8
#define HAF_RS_READ      0x10
#define HAF_CREATE_NEW   0x20
#define HAF_SKIP_READ    0x40

typedef struct { // global data structure for kt_pipeline()
	const yak_copt_t *opt;
	const void *flt_tab;
	int flag, create_new, is_store;
	uint64_t n_seq; ///number of total reads
	kseq_t *ks;
	UC_Read ucr;
	ha_ct_t *ct;
	ha_pt_t *pt;
	const All_reads *rs_in;
	All_reads *rs_out;
} pl_data_t;

typedef struct { // data structure for each step in kt_pipeline()
	pl_data_t *p;
	uint64_t n_seq0;  ///the start index of current buffer block at R_INF
	///sum_len = total bases, nk = number of k-mers
	int n_seq, m_seq, sum_len, nk;
	int *len;
	char **seq;
	ha_mz1_v *mz_buf;
	ha_mz1_v *mz;
	ch_buf_t *buf;
} st_data_t;

static void worker_for_insert(void *data, long i, int tid) // callback for kt_for()
{
	st_data_t *s = (st_data_t*)data;
	ch_buf_t *b = &s->buf[i];
	if (s->p->pt)
		b->n_ins += ha_pt_insert_list(s->p->pt, b->n, b->b);
	else///for 0-th count, go into here
		b->n_ins += ha_ct_insert_list(s->p->ct, s->p->create_new, b->n, b->a);
}

static void worker_for_mz(void *data, long i, int tid)
{
	st_data_t *s = (st_data_t*)data;
	///get the corresponding minimzer vector of this read
	ha_mz1_v *b = &s->mz_buf[tid];
	s->mz_buf[tid].n = 0;
	///s->p->opt->w = 51, s->p->opt->k
	ha_sketch(s->seq[i], s->len[i], s->p->opt->w, s->p->opt->k, s->n_seq0 + i, s->p->opt->is_HPC, b, s->p->flt_tab);
	s->mz[i].n = s->mz[i].m = b->n;
	MALLOC(s->mz[i].a, b->n);
	memcpy(s->mz[i].a, b->a, b->n * sizeof(ha_mz1_t));
}

static void *worker_count(void *data, int step, void *in) // callback for kt_pipeline()
{
	pl_data_t *p = (pl_data_t*)data;
	if (step == 0) { // step 1: read a block of sequences
		int ret;
		st_data_t *s;
		CALLOC(s, 1);
		s->p = p;
		s->n_seq0 = p->n_seq;
		if (p->rs_in && (p->flag & HAF_RS_READ)) {
			while (p->n_seq < p->rs_in->total_reads) {
				if((p->flag & HAF_SKIP_READ) && p->rs_in->trio_flag[p->n_seq] != AMBIGU)
				{
					++p->n_seq;
					continue;
				}
				int l;
				recover_UC_Read(&p->ucr, p->rs_in, p->n_seq);
				l = p->ucr.length;
				if (s->n_seq == s->m_seq) {
					s->m_seq = s->m_seq < 16? 16 : s->m_seq + (s->m_seq>>1);
					REALLOC(s->len, s->m_seq);
					REALLOC(s->seq, s->m_seq);
				}
				MALLOC(s->seq[s->n_seq], l);
				memcpy(s->seq[s->n_seq], p->ucr.seq, l);
				s->len[s->n_seq++] = l;
				++p->n_seq;
				s->sum_len += l;
				s->nk += l >= p->opt->k? l - p->opt->k + 1 : 0;
				if (s->sum_len >= p->opt->chunk_size)
					break;
			}
		} else {
			while ((ret = kseq_read(p->ks)) >= 0) {
				int l = (int)(p->ks->seq.l) - (int)(p->opt->adaLen) - (int)(p->opt->adaLen);
				if(l <= 0) continue;

				if (p->n_seq >= 1<<28) {
					fprintf(stderr, "ERROR: this implementation supports no more than %d reads\n", 1<<28);
					exit(1);
				}
				if (p->rs_out) {
					///for 0-th count, just insert read length to R_INF, instead of read
					if (p->flag & HAF_RS_WRITE_LEN) {
						assert(p->n_seq == p->rs_out->total_reads);
						ha_insert_read_len(p->rs_out, l, p->ks->name.l);
					} else if (p->flag & HAF_RS_WRITE_SEQ) {
						int i, n_N;
						assert(l == (int)p->rs_out->read_length[p->n_seq]);
						for (i = n_N = 0; i < l; ++i) // count number of ambiguous bases
							if (seq_nt4_table[(uint8_t)p->ks->seq.s[i+p->opt->adaLen]] >= 4)
								++n_N;
						ha_compress_base(Get_READ(*p->rs_out, p->n_seq), p->ks->seq.s+p->opt->adaLen, l, &p->rs_out->N_site[p->n_seq], n_N);
						memcpy(&p->rs_out->name[p->rs_out->name_index[p->n_seq]], p->ks->name.s, p->ks->name.l);
					}
				}
				///for 0-th count, insert both seq and length to local block
				if (s->n_seq == s->m_seq) {
					s->m_seq = s->m_seq < 16? 16 : s->m_seq + (s->m_seq>>1);
					REALLOC(s->len, s->m_seq);
					REALLOC(s->seq, s->m_seq);
				}
				MALLOC(s->seq[s->n_seq], l);
				memcpy(s->seq[s->n_seq], p->ks->seq.s+p->opt->adaLen, l);
				s->len[s->n_seq++] = l;
				++p->n_seq;
				s->sum_len += l;
				s->nk += l >= p->opt->k? l - p->opt->k + 1 : 0;
				///p->opt->chunk_size is the block max size
				if (s->sum_len >= p->opt->chunk_size)
					break;
			}
		}
		if (s->sum_len == 0) free(s);
		else return s;
	} else if (step == 1) { // step 2: extract k-mers
		///s is the block of reads
		st_data_t *s = (st_data_t*)in;
		///for 0-th counting, n_pre = 4096
		int i, n_pre = 1<<p->opt->pre, m;
		// allocate the k-mer buffer
		CALLOC(s->buf, n_pre);
		m = (int)(s->nk * 1.2 / n_pre) + 1;
		//pre-allocate memory for each of 4096 buffer
		for (i = 0; i < n_pre; ++i) {
			s->buf[i].m = m;
			///for 0-th counting, p->pt = NULL
			if (p->pt) MALLOC(s->buf[i].b, m);
			else MALLOC(s->buf[i].a, m);
		}
		// fill the buffer
		///for 0-th counting, p->opt->w == 1
		if (p->opt->w == 1) { // enumerate all k-mers
			///scan all reads
			for (i = 0; i < s->n_seq; ++i) {
				if (p->opt->is_HPC)
					count_seq_buf_HPC(s->buf, p->opt->k, p->opt->pre, s->len[i], s->seq[i]);
				else
					count_seq_buf(s->buf, p->opt->k, p->opt->pre, s->len[i], s->seq[i]);
				if (!p->is_store) free(s->seq[i]);
			}
		} else { // minimizers only
			uint32_t j;
			// compute minimizers
			// s->n_seq is how many reads at this buffer
			// s->mz && s->mz_buf are lists of minimzer vectors
			CALLOC(s->mz, s->n_seq);
			CALLOC(s->mz_buf, p->opt->n_thread);
			///calculate minimzers for each read, each read corresponds to one thread
			kt_for(p->opt->n_thread, worker_for_mz, s, s->n_seq);
			for (i = 0; i < p->opt->n_thread; ++i)
				free(s->mz_buf[i].a);
			free(s->mz_buf);
			// insert minimizers
			if (p->pt) {///insert whole minimizer
				for (i = 0; i < s->n_seq; ++i)
					for (j = 0; j < s->mz[i].n; ++j)
						pt_insert_buf(s->buf, p->opt->pre, &s->mz[i].a[j]);
			} else {///just insert the hash key of minimizer
				for (i = 0; i < s->n_seq; ++i)
					for (j = 0; j < s->mz[i].n; ++j)
						ct_insert_buf(s->buf, p->opt->pre, s->mz[i].a[j].x);
			}
			for (i = 0; i < s->n_seq; ++i) {
				free(s->mz[i].a);
				if (!p->is_store) free(s->seq[i]);
			}
			free(s->mz);
		}
		///just clean seq
		free(s->seq); free(s->len);
		s->seq = 0, s->len = 0;
		return s;
	} else if (step == 2) { // step 3: insert k-mers to hash table
		st_data_t *s = (st_data_t*)in;
		int i, n = 1<<p->opt->pre;
		uint64_t n_ins = 0;
		///for 0-th counting, p->pt = NULL
		kt_for(p->opt->n_thread, worker_for_insert, s, n);
		///n_ins is number of distinct k-mers
		for (i = 0; i < n; ++i) {
			n_ins += s->buf[i].n_ins;
			if (p->pt) free(s->buf[i].b);
			else free(s->buf[i].a);
		}
		if (p->ct) p->ct->tot += n_ins;
		if (p->pt) p->pt->tot_pos += n_ins;
		free(s->buf);
		#if 0
		fprintf(stderr, "[M::%s::%.3f*%.2f] processed %ld sequences; %ld %s in the hash table\n", __func__,
				yak_realtime(), yak_cpu_usage(), (long)s->n_seq0 + s->n_seq,
				(long)(p->pt? p->pt->tot_pos : p->ct->tot), p->pt? "positions" : "distinct k-mers");
		#endif
		free(s);
	}
	return 0;
}

void debug_adapter(const hifiasm_opt_t *asm_opt, All_reads *rs)
{
	int ret;
    uint32_t i, m, pass, unpass;
    gzFile fp = 0;
    kseq_t *ks = NULL;
	UC_Read ucr;
	init_UC_Read(&ucr);

    for (i = m = pass = unpass = 0; i < (uint32_t)asm_opt->num_reads; ++i)
    {
		if ((fp = gzopen(asm_opt->read_file_names[i], "r")) == 0) continue;
		ks = kseq_init(fp);
		while ((ret = kseq_read(ks)) >= 0) 
		{
			int l = ks->seq.l;
			if((l - asm_opt->adapterLen*2) <= 0) continue;
			recover_UC_Read(&ucr, rs, m);
			fprintf(stderr, "l: %d, ucr.length: %lld, asm_opt->adapterLen: %d\n", 
			l, ucr.length, asm_opt->adapterLen);
			if(memcmp(ucr.seq, ks->seq.s+asm_opt->adapterLen, ucr.length) == 0)
			{
				pass++;
			}
			else
			{
				unpass++;
			}
			m++;
		}
        kseq_destroy(ks);
        gzclose(fp);
		ks = NULL;
		fp = 0;
    }

	destory_UC_Read(&ucr);

	fprintf(stderr, "[M::%s::# reads: %u, # pass: %u, # unpass: %u\n]", __func__, m, pass, unpass);
	exit(1);
}

static ha_ct_t *yak_count(const yak_copt_t *opt, const char *fn, int flag, ha_pt_t *p0, ha_ct_t *c0, const void *flt_tab, All_reads *rs, int64_t *n_seq)
{
	///for 0-th counting, flag = HAF_COUNT_ALL|HAF_RS_WRITE_LEN|HAF_CREATE_NEW
	int read_rs = (rs && (flag & HAF_RS_READ));
	pl_data_t pl;
	gzFile fp = 0;
	memset(&pl, 0, sizeof(pl_data_t));
	pl.n_seq = *n_seq;
	if (read_rs) {
		pl.rs_in = rs;
		init_UC_Read(&pl.ucr);
	} else {///for 0-th counting, go into here
		if ((fp = gzopen(fn, "r")) == 0) return 0;
		pl.ks = kseq_init(fp);
	}
	///for 0-th counting, read all reads into pl.rs_out
	if (rs && (flag & (HAF_RS_WRITE_LEN|HAF_RS_WRITE_SEQ)))
		pl.rs_out = rs;
	///for 0-th counting, flt_tab = NULL
	///for 1-th counting, flt_tab = NULL
	pl.flt_tab = flt_tab;
	pl.opt = opt;
	pl.flag = flag;
	if (p0) {///for 1-th counting, p0 = NULL
		pl.pt = p0, pl.create_new = 0; // never create new elements in a position table
		assert(p0->k == opt->k && p0->pre == opt->pre);
	} else if (c0) {
		pl.ct = c0, pl.create_new = !!(flag&HAF_CREATE_NEW);
		assert(c0->k == opt->k && c0->pre == opt->pre);
	} else {///for ft-th counting and 1-th counting, go into here
		pl.create_new = 1; // alware create new elements if the count table is empty
		///for 0-th counting, opt.k = 51, opt->pre = 12, opt->bf_n_hash = 4, opt.bf_shift = 37
		///for 1-th counting, opt.k = 51, opt->pre = 12, opt->bf_n_hash = 4, opt.bf_shift = 0
		///building a large hash table consisting of 4096 small hash tables
		pl.ct = ha_ct_init(opt->k, opt->pre, opt->bf_n_hash, opt->bf_shift);
	}
	kt_pipeline(3, worker_count, &pl, 3);
	if (read_rs) {
		destory_UC_Read(&pl.ucr);
	} else {
		kseq_destroy(pl.ks);
		gzclose(fp);
	}
	*n_seq = pl.n_seq;
	return pl.ct;
}

ha_ct_t *ha_count(const hifiasm_opt_t *asm_opt, int flag, ha_pt_t *p0, const void *flt_tab, All_reads *rs)
{
	int i;
	int64_t n_seq = 0;
	yak_copt_t opt;
	ha_ct_t *h = 0;
	assert(!(flag & HAF_RS_WRITE_LEN) || !(flag & HAF_RS_WRITE_SEQ)); // not both
	///for 0-th counting, flag = HAF_COUNT_ALL|HAF_RS_WRITE_LEN
	if (rs) {
		if (flag & HAF_RS_WRITE_LEN)
			init_All_reads(rs);
		else if (flag & HAF_RS_WRITE_SEQ)
			malloc_All_reads(rs);
	}
	yak_copt_init(&opt);
	opt.k = asm_opt->k_mer_length;
	///always 0
	opt.is_HPC = !(asm_opt->flag&HA_F_NO_HPC);
	///for ft-counting, shoud be 1
	opt.w = flag & HAF_COUNT_ALL? 1 : asm_opt->mz_win;
	///for ft-counting, shoud be 37
	///for ha_pt_gen, shoud be 0
	opt.bf_shift = flag & HAF_COUNT_EXACT? 0 : asm_opt->bf_shift;
	opt.n_thread = asm_opt->thread_num;
	opt.adaLen = asm_opt->adapterLen;
	///asm_opt->num_reads is the number of fastq files
	for (i = 0; i < asm_opt->num_reads; ++i)
		h = yak_count(&opt, asm_opt->read_file_names[i], flag|HAF_CREATE_NEW, p0, h, flt_tab, rs, &n_seq);
	if (h && opt.bf_shift > 0)
		ha_ct_destroy_bf(h);
	return h;
}

/***************************
 * High count filter table *
 ***************************/

KHASHL_SET_INIT(static klib_unused, yak_ft_t, yak_ft, uint64_t, kh_hash_dummy, kh_eq_generic)

static yak_ft_t *gen_hh(const ha_ct_t *h)
{
	int i;
	yak_ft_t *hh;
	hh = yak_ft_init();
	yak_ft_resize(hh, h->tot * 2);
	for (i = 0; i < 1<<h->pre; ++i) {
		yak_ct_t *ht = h->h[i].h;
		khint_t k;
		for (k = 0; k < kh_end(ht); ++k) {
			if (kh_exist(ht, k)) {
				uint64_t y = kh_key(ht, k) >> h->pre << YAK_COUNTER_BITS | i;
				int absent;
				yak_ft_put(hh, y, &absent);
			}
		}
	}
	return hh;
}

int ha_ft_isflt(const void *hh, uint64_t y)
{
	yak_ft_t *h = (yak_ft_t*)hh;
	khint_t k;
	k = yak_ft_get(h, y);
	return k == kh_end(h)? 0 : 1;
}

void ha_ft_destroy(void *h)
{
	if (h) yak_ft_destroy((yak_ft_t*)h);
}


void debug_ct_index(void* q_ct_idx, void* r_ct_idx)
{
	ha_ct_t* ct_idx = (ha_ct_t*)q_ct_idx;
	yak_ct_t *g = NULL;
	uint64_t i;
	khint_t k;
	for (i = 0; (int)i < 1<<ct_idx->pre; i++)
	{
		g = ct_idx->h[i].h;
		for (k = 0; k < kh_end(g); ++k) 
		{
			if (kh_exist(g, k)) 
			{
				int c = kh_key(g, k) & YAK_MAX_COUNT;
				uint64_t hash = ((kh_key(g, k) >> ct_idx->pre)<<ct_idx->pre) | i;
				int q = query_ct_index(r_ct_idx, hash);
				if(q!=c)
				{
					fprintf(stderr, "ERROR:c: %d, q: %d\n", c, q);
				}
			}
		}
	}
}

/*************************
 * High-level interfaces *
 *************************/

void *ha_ft_gen(const hifiasm_opt_t *asm_opt, All_reads *rs, int *hom_cov, int is_hp_mode)
{
	yak_ft_t *flt_tab;
	int64_t cnt[YAK_N_COUNTS];
	int peak_hom, peak_het, cutoff = YAK_MAX_COUNT - 1, ex_flag = 0;
	if(is_hp_mode) ex_flag = HAF_RS_READ|HAF_SKIP_READ;
	ha_ct_t *h;
	h = ha_count(asm_opt, HAF_COUNT_ALL|HAF_RS_WRITE_LEN|ex_flag, NULL, NULL, rs);
	if((asm_opt->flag & HA_F_VERBOSE_GFA))
	{
		write_ct_index((void*)h, asm_opt->output_file_name);
		// load_ct_index(&ha_ct_table, asm_opt->output_file_name);
		// debug_ct_index((void*)h, ha_ct_table);
		// debug_ct_index(ha_ct_table, (void*)h);
		// ha_ct_destroy((ha_ct_t *)ha_ct_table);
	}
	
	if(!(ex_flag & HAF_SKIP_READ))
	{
		ha_ct_hist(h, cnt, asm_opt->thread_num);
		peak_hom = ha_analyze_count(YAK_N_COUNTS, asm_opt->min_hist_kmer_cnt, cnt, &peak_het);
		if (hom_cov) *hom_cov = peak_hom;
		if (peak_hom > 0) fprintf(stderr, "[M::%s] peak_hom: %d; peak_het: %d\n", __func__, peak_hom, peak_het);
		///in default, asm_opt->high_factor = 5.0
		cutoff = (int)(peak_hom * asm_opt->high_factor);
		if (cutoff > YAK_MAX_COUNT - 1) cutoff = YAK_MAX_COUNT - 1;
	}
	ha_ct_shrink(h, cutoff, YAK_MAX_COUNT, asm_opt->thread_num);
	flt_tab = gen_hh(h);
	ha_ct_destroy(h);
	fprintf(stderr, "[M::%s::%.3f*%.2f@%.3fGB] ==> filtered out %ld k-mers occurring %d or more times\n", __func__,
			yak_realtime(), yak_cpu_usage(), yak_peakrss_in_gb(), (long)kh_size(flt_tab), cutoff);
	return (void*)flt_tab;
}

ha_pt_t *ha_pt_gen(const hifiasm_opt_t *asm_opt, const void *flt_tab, int read_from_store, int is_hp_mode, All_reads *rs, int *hom_cov, int *het_cov)
{
	int64_t cnt[YAK_N_COUNTS], tot_cnt;
	int peak_hom, peak_het, i, extra_flag1, extra_flag2;
	ha_ct_t *ct;
	ha_pt_t *pt;
	if (read_from_store) {///if reads have already been read
		extra_flag1 = extra_flag2 = HAF_RS_READ;
	} else if (rs->total_reads == 0) {///if reads & length have not been scanned
		extra_flag1 = HAF_RS_WRITE_LEN;
		extra_flag2 = HAF_RS_WRITE_SEQ;
	} else {///if length has been loaded but reads have not
		extra_flag1 = HAF_RS_WRITE_SEQ;
		extra_flag2 = HAF_RS_READ;
	}
	if(is_hp_mode) extra_flag1 |= HAF_SKIP_READ, extra_flag2 |= HAF_SKIP_READ;

	ct = ha_count(asm_opt, HAF_COUNT_EXACT|extra_flag1, NULL, flt_tab, rs);
	fprintf(stderr, "[M::%s::%.3f*%.2f] ==> counted %ld distinct minimizer k-mers\n", __func__,
			yak_realtime(), yak_cpu_usage(), (long)ct->tot);
	ha_ct_hist(ct, cnt, asm_opt->thread_num);
	fprintf(stderr, "[M::%s] count[%d] = %ld (for sanity check)\n", __func__, YAK_MAX_COUNT, (long)cnt[YAK_MAX_COUNT]);
	peak_hom = ha_analyze_count(YAK_N_COUNTS, asm_opt->min_hist_kmer_cnt, cnt, &peak_het);
	if (hom_cov) *hom_cov = peak_hom;
	if (het_cov) *het_cov = peak_het;
	if (peak_hom > 0) fprintf(stderr, "[M::%s] peak_hom: %d; peak_het: %d\n", __func__, peak_hom, peak_het);
	///here ha_ct_shrink is mostly used to remove k-mer appearing only 1 time
	if (flt_tab == 0) {
		int cutoff = (int)(peak_hom * asm_opt->high_factor);
		if (cutoff > YAK_MAX_COUNT - 1) cutoff = YAK_MAX_COUNT - 1;
		if((extra_flag1 & HAF_SKIP_READ) && (extra_flag2 & HAF_SKIP_READ)) cutoff = YAK_MAX_COUNT - 1;
		ha_ct_shrink(ct, 2, cutoff, asm_opt->thread_num);
		for (i = 2, tot_cnt = 0; i <= cutoff; ++i) tot_cnt += cnt[i] * i;
	} else {
		///Note: here is just to remove minimizer appearing YAK_MAX_COUNT times
		///minimizer with YAK_MAX_COUNT occ may apper > YAK_MAX_COUNT times, so it may lead to overflow at ha_pt_gen
		ha_ct_shrink(ct, 2, YAK_MAX_COUNT - 1, asm_opt->thread_num);
		for (i = 2, tot_cnt = 0; i <= YAK_MAX_COUNT - 1; ++i) tot_cnt += cnt[i] * i;
	}
	pt = ha_pt_gen(ct, asm_opt->thread_num);
	ha_count(asm_opt, HAF_COUNT_EXACT|extra_flag2, pt, flt_tab, rs);
	assert((uint64_t)tot_cnt == pt->tot_pos);
	//ha_pt_sort(pt, asm_opt->thread_num);
	fprintf(stderr, "[M::%s::%.3f*%.2f] ==> indexed %ld positions\n", __func__,
			yak_realtime(), yak_cpu_usage(), (long)pt->tot_pos);
	return pt;
}

int query_ct_index(void* ct_idx, uint64_t hash)
{
	ha_ct1_t *g = &(((ha_ct_t*)ct_idx)->h[hash & ((1ULL<<((ha_ct_t*)ct_idx)->pre) - 1)]);
	khint_t k;
	k = yak_ct_get(g->h, hash);
	if (k == kh_end(g->h)) return 0;
	return kh_key(g->h, k)&YAK_MAX_COUNT;
}


int write_ct_index(void *i_ct_idx, char* file_name)
{
	char* gfa_name = (char*)malloc(strlen(file_name)+25);
    sprintf(gfa_name, "%s.ct_flt", file_name);
    FILE* fp = fopen(gfa_name, "w");
	if (!fp) {
		free(gfa_name);
        return 0;
    }
	ha_ct_t* ct_idx = (ha_ct_t*)i_ct_idx;
	int i;
	ha_ct1_t *g;
	fwrite(&ct_idx->k, sizeof(ct_idx->k), 1, fp);
	fwrite(&ct_idx->pre, sizeof(ct_idx->pre), 1, fp);
	fwrite(&ct_idx->n_hash, sizeof(ct_idx->n_hash), 1, fp);
	fwrite(&ct_idx->n_shift, sizeof(ct_idx->n_shift), 1, fp);
	fwrite(&ct_idx->tot, sizeof(ct_idx->tot), 1, fp);
	for (i = 0; i < 1<<ct_idx->pre; i++)
	{
		g = &(ct_idx->h[i]);
		yak_ct_save(g->h, fp);
	}


	fprintf(stderr, "[M::%s] Index has been written.\n", __func__);
    free(gfa_name);
	fclose(fp);
	return 1;
}

int load_ct_index(void **i_ct_idx, char* file_name)
{
	char* gfa_name = (char*)malloc(strlen(file_name)+25);
    sprintf(gfa_name, "%s.ct_flt", file_name);
    FILE* fp = fopen(gfa_name, "r");
	if (!fp) {
		free(gfa_name);
        return 0;
    }
	ha_ct_t** ct_idx = (ha_ct_t**)i_ct_idx;
	double index_time = 0;
	uint64_t flag = 0;
	int i;
	ha_ct_t *h = 0;
	ha_ct1_t *g;
	CALLOC(h, 1);

	flag += fread(&h->k, sizeof(h->k), 1, fp);
	flag += fread(&h->pre, sizeof(h->pre), 1, fp);
	flag += fread(&h->n_hash, sizeof(h->n_hash), 1, fp);
	flag += fread(&h->n_shift, sizeof(h->n_shift), 1, fp);
	flag += fread(&h->tot, sizeof(h->tot), 1, fp);
	CALLOC(h->h, 1<<h->pre);


	index_time = yak_realtime();
	for (i = 0; i < 1<<h->pre; ++i) 
	{
		g = &(h->h[i]);
		yak_ct_load(&(g->h), fp);
	}

	(*ct_idx) = h;
	fprintf(stderr, "[M::%s::%.3f] ==> Loaded count table\n", __func__, yak_realtime() - index_time);
	fprintf(stderr, "[M::%s] Index has been loaded.\n", __func__);
	free(gfa_name);
	return 1;
}

int write_pt_index(void *flt_tab, ha_pt_t *ha_idx, All_reads* r, hifiasm_opt_t* opt, char* file_name)
{
    char* gfa_name = (char*)malloc(strlen(file_name)+25);
    sprintf(gfa_name, "%s.pt_flt", file_name);
    FILE* fp = fopen(gfa_name, "w");
	if (!fp) {
		free(gfa_name);
        return 0;
    }
	yak_ft_t *ha_flt_tab = (yak_ft_t*)flt_tab;

	if(ha_flt_tab)
	{
		fwrite("f", 1, 1, fp);
		yak_ft_save(ha_flt_tab, fp);
	}


	if(ha_idx)
	{
		int i;
		ha_pt1_t *g;
		fwrite("h", 1, 1, fp);
		fwrite(&ha_idx->k, sizeof(ha_idx->k), 1, fp);
		fwrite(&ha_idx->pre, sizeof(ha_idx->pre), 1, fp);
		fwrite(&ha_idx->tot, sizeof(ha_idx->tot), 1, fp);
		fwrite(&ha_idx->tot_pos, sizeof(ha_idx->tot_pos), 1, fp);

		for (i = 0; i < 1<<ha_idx->pre; ++i) 
		{
			g = &(ha_idx->h[i]);
			yak_pt_save(g->h, fp);
			fwrite(&g->n, sizeof(g->n), 1, fp);
			fwrite(g->a, sizeof(ha_idxpos_t), g->n, fp);
		}
	}

	fwrite(&opt->number_of_round, sizeof(opt->number_of_round), 1, fp);
	fwrite(&opt->hom_cov, sizeof(opt->hom_cov), 1, fp);
	fwrite(&opt->het_cov, sizeof(opt->het_cov), 1, fp);
	fwrite(&opt->max_n_chain, sizeof(opt->max_n_chain), 1, fp);


	write_All_reads(r, gfa_name);

	fprintf(stderr, "[M::%s] Index has been written.\n", __func__);
    free(gfa_name);
	fclose(fp);
	return 1;
}

int load_pt_index(void **r_flt_tab, ha_pt_t **r_ha_idx, All_reads* r, hifiasm_opt_t* opt, char* file_name)
{
	char* gfa_name = (char*)malloc(strlen(file_name)+25);
    sprintf(gfa_name, "%s.pt_flt", file_name);
    FILE* fp = fopen(gfa_name, "r");
	if (!fp) {
		free(gfa_name);
        return 0;
    }

	ha_pt_t *ha_idx = NULL;
	char mode = 0;
	int f_flag, absent, i;
	double index_time, index_s_time, pos_time, pos_s_time;

	
	
	f_flag += fread(&mode, 1, 1, fp);
	if(mode == 'f')
	{
		index_time = yak_realtime();

		yak_ft_load((yak_ft_t **)r_flt_tab, fp);
		
	    f_flag += fread(&mode, 1, 1, fp);

		fprintf(stderr, "[M::%s::%.3f] ==> Loaded flt table\n", __func__, yak_realtime()-index_time);
	}
	///insert using multiple threads???
	if(mode == 'h')
	{
		pos_time = index_time = 0;

		CALLOC(ha_idx, 1);
		ha_pt1_t *g;
		f_flag += fread(&ha_idx->k, sizeof(ha_idx->k), 1, fp);
		f_flag += fread(&ha_idx->pre, sizeof(ha_idx->pre), 1, fp);
		f_flag += fread(&ha_idx->tot, sizeof(ha_idx->tot), 1, fp);
		f_flag += fread(&ha_idx->tot_pos, sizeof(ha_idx->tot_pos), 1, fp);
		CALLOC(ha_idx->h, 1<<ha_idx->pre);
		for (i = 0; i < 1<<ha_idx->pre; ++i) 
		{
			index_s_time = yak_realtime();

			g = &(ha_idx->h[i]);
			yak_pt_load(&(g->h), fp);

			index_time += yak_realtime() - index_s_time;
			
			pos_s_time = yak_realtime();

			f_flag += fread(&g->n, sizeof(g->n), 1, fp);
			MALLOC(g->a, g->n);
			f_flag += fread(g->a, sizeof(ha_idxpos_t), g->n, fp);

			pos_time += yak_realtime() - pos_s_time;
		}
		(*r_ha_idx) = ha_idx;

		fprintf(stderr, "[M::%s::%.3f(index)/%.3f(pos)] ==> Loaded pos table\n", __func__, index_time, pos_time);
	}

	if(mode != 'h' && mode != 'f')
	{
		free(gfa_name);
		fclose(fp);
		return 0;
	}


	f_flag += fread(&absent, sizeof(absent), 1, fp);
	if(absent != opt->number_of_round)
	{
		fprintf(stderr, "ERROR: different number of rounds!\n");
		exit(1);
	}

	f_flag += fread(&opt->hom_cov, sizeof(opt->hom_cov), 1, fp);
	f_flag += fread(&opt->het_cov, sizeof(opt->het_cov), 1, fp);
	f_flag += fread(&opt->max_n_chain, sizeof(opt->max_n_chain), 1, fp);


	fclose(fp);
	
	if(!load_All_reads(r, gfa_name))
	{
		free(gfa_name);
		return 0;
	}


	memset(r->trio_flag, AMBIGU, r->total_reads*sizeof(uint8_t));
	r->paf = (ma_hit_t_alloc*)malloc(sizeof(ma_hit_t_alloc)*r->total_reads);
    r->reverse_paf = (ma_hit_t_alloc*)malloc(sizeof(ma_hit_t_alloc)*r->total_reads);
	for (i = 0; i < (long long)r->total_reads; i++)
    {
        init_ma_hit_t_alloc(&(r->paf[i]));
        init_ma_hit_t_alloc(&(r->reverse_paf[i]));
    }

	fprintf(stderr, "[M::%s] Index has been loaded.\n", __func__);

	free(gfa_name);
	return 1;
}
