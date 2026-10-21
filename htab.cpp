#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <pthread.h>
#include <sys/mman.h>
#include "kthread.h"
#include "khashl.h"
#include "kseq.h"
#include "ksort.h"
#include "htab.h"
#include "Process_Read.h"

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
	int32_t k, w, is_HPC, uq;
	int32_t pre;
	int32_t n_thread;
	int64_t chunk_size;
	int adaLen, min_rcnt;
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
	///the large filters (16 GB of random inserts in the filter pass) aligned to 2 MB: with THP enabled
	///they get 2 MB pages when the kernel has free ones. No MADV_HUGEPAGE: with defrag=madvise it makes
	///the kernel compact memory for them, which costs more than it gives on a fragmented host.
	if (n_shift - 3 >= 21) posix_memalign(&ptr, 1ULL<<21, 1ULL<<(n_shift-3));
	else posix_memalign(&ptr, 1<<(YAK_BLK_SHIFT-3), 1ULL<<(n_shift-3));
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

///prefetch the block of yak_bf_insert(b, hash): no change of the filter
static inline void yak_bf_prefetch(const yak_bf_t *b, uint64_t hash)
{
	int x = b->n_shift - YAK_BLK_SHIFT;
	uint64_t y = hash & ((1ULL<<x) - 1);
	__builtin_prefetch(&b->b[y<<(YAK_BLK_SHIFT-3)], 1);
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
	uint64_t bs;
	ha_ct1_t *h;
	void *bf_map; uint64_t bf_map_n; ///the one mapping of all bloom filters (ha_ct_init()), or NULL
} ha_ct_t;

///for 0-th counting, k = 51, pre = 12, n_hash = 4, n_shift = 37
///for 1-th counting, opt.k = 51, opt->pre = 12, opt->bf_n_hash = 4, opt.bf_shift = 0
///the 4096 bloom filters of a count table take 2^(n_shift-3) bytes together (16 GB for n_shift = 37).
///They are slices of one anonymous mapping, 2 MB-aligned: a mmap() or munmap() per filter takes the
///address-space lock for writing and stalls the page faults of all other threads (T1: 0.8 s at about
///4 busy cores). The mapping is zero; the workers fault the pages in parallel with one write per 4 KB,
///and at the end free them in parallel (MADV_DONTNEED) before the one munmap(). The contents are the
///same as with yak_bf_init().
static void ha_ct_bf_init_worker(void *data, long i, int tid)///callback for kt_for()
{
	ha_ct_t *h = (ha_ct_t*)data;
	uint64_t j, n = 1ULL<<(h->n_shift - h->pre - 3);
	yak_bf_t *b;
	CALLOC(b, 1);
	b->n_shift = h->n_shift - h->pre;
	b->n_hashes = h->n_hash;
	b->b = (uint8_t*)(((uintptr_t)h->bf_map + (1ULL<<21) - 1) & ~(uintptr_t)((1ULL<<21) - 1)) + (uint64_t)i * n;
	for (j = 0; j < n; j += 4096) b->b[j] = 0;
	h->h[i].b = b;
}

static void ha_ct_bf_destroy_worker(void *data, long i, int tid)///callback for kt_for()
{
	ha_ct_t *h = (ha_ct_t*)data;
	if (h->h[i].b) {
		if (h->bf_map) {///a slice of the mapping: free its pages now, the mapping after the kt_for()
			madvise(h->h[i].b->b, 1ULL<<(h->n_shift - h->pre - 3), MADV_DONTNEED);
			free(h->h[i].b);
		} else yak_bf_destroy(h->h[i].b);
	}
	h->h[i].b = 0;
}

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
		///h->n_shift = 37, h->pre = 12, h->n_hash = 4 
		///the bounds of yak_bf_init(): else there are no filters, as there
		if (n_shift - h->pre + YAK_BLK_SHIFT <= 64 && n_shift - h->pre >= YAK_BLK_SHIFT) {
			h->bf_map_n = (1ULL<<(n_shift - 3)) + (1ULL<<21);
			h->bf_map = mmap(0, h->bf_map_n, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
			if (h->bf_map == MAP_FAILED) {
				fprintf(stderr, "[E::%s] cannot map %llu bytes for the bloom filters\n", __func__, (unsigned long long)h->bf_map_n);
				exit(1);
			}
			kt_for(asm_opt.thread_num > 0? asm_opt.thread_num : 1, ha_ct_bf_init_worker, h, 1<<h->pre);
		}
	}
	return h;
}

static void ha_ct_destroy_bf(ha_ct_t *h)
{
	kt_for(asm_opt.thread_num > 0? asm_opt.thread_num : 1, ha_ct_bf_destroy_worker, h, 1<<h->pre);
	if (h->bf_map) munmap(h->bf_map, h->bf_map_n);
	h->bf_map = 0; h->bf_map_n = 0;
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

#define HA_BF_PREFETCH 16 ///the bloom block of the k-mer HA_BF_PREFETCH ahead
static int ha_ct_insert_list(ha_ct_t *h, int create_new, int n, const uint64_t *a)
{
	int j, mask = (1<<h->pre) - 1, n_ins = 0;
	ha_ct1_t *g;
	if (n == 0) return 0;
	///corresponding small hash index
	g = &h->h[a[0]&mask];
	for (j = 0; j < n; ++j) {
		int ins = 1, absent;
		if (create_new && g->b && j + HA_BF_PREFETCH < n) yak_bf_prefetch(g->b, a[j + HA_BF_PREFETCH] >> h->pre);///a random block for each k-mer
		///x is a 64-bit word, h->pre=12
		///all elements at a have the same low 12 bits
		///so low 12 bits are not useful
		uint64_t x = a[j] >> h->pre;
		khint_t k;
		assert((a[j]&mask) == (a[0]&mask));
		if (create_new) {
			if (g->b)
				ins = (yak_bf_insert(g->b, x) == h->n_hash);
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
	fprintf(stderr, "[M::%s::%.3f*%.2f] ==> counted %ld distinct minimizer k-mers\n", __func__,
			yak_realtime(), yak_cpu_usage(), (long)h->tot);
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
	ha_idxposl_t *al;
} ha_pt1_t;

struct ha_pt_s {
	int k, pre;
	uint64_t tot, tot_pos;
	ha_pt1_t *h;
};

typedef struct {
	const ha_ct_t *ct;
	ha_pt_t *pt;
	int is_l;
} pt_gen_aux_t;


static void worker_pt_shrink(void *data, long i, int tid) // callback for kt_for()
{
    ha_pt_t *h = (ha_pt_t*)data;
	ha_pt1_t *b = &h->h[i];
    yak_pt_t *f = NULL;
    khint_t k;
    f = yak_pt_init();
    for (k = 0, b->n = 0; k < kh_end(b->h); ++k) {
        if (kh_exist(b->h, k)) {
			if(kh_val(b->h, k) <= 0) continue;
			int absent; khint_t l;
			l = yak_pt_put(f, (kh_key(b->h, k) >> h->pre) << YAK_COUNTER_BITS, &absent);
			kh_val(f, l) = b->n; 
			b->n += kh_key(b->h, k) & YAK_MAX_COUNT;
        }
    }
    yak_pt_destroy(b->h);
    h->h[i].h = f;
	CALLOC(b->a, b->n);///need fix
}

static uint64_t ha_pt_shrink(ha_pt_t *h, int n_thread)
{
    int i;
	uint64_t occ;
    ///still start 4096 threads
    kt_for(n_thread, worker_pt_shrink, h, 1<<h->pre);
    for (i = 0, occ = 0, h->tot = 0; i < 1<<h->pre; ++i)
	{
		h->tot += kh_size(h->h[i].h);
		occ += h->h[i].n;
	}
	return occ;
}

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
	if(a->is_l) CALLOC(b->al, b->n);
	else CALLOC(b->a, b->n);
}

ha_pt_t *ha_pt_gen(ha_ct_t *ct, int n_thread, int is_l)
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
	a.ct = ct, a.pt = pt, a.is_l = is_l;
	kt_for(n_thread, worker_pt_gen, &a, 1<<pt->pre);
	free(ct->h); free(ct);
	return pt;
}

static void worker_pt_gen_count(void *data, long i, int tid) // callback for kt_for()
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
			kh_val(b->h, l) = 0; kh_key(b->h, l) |= kh_key(g, k) & YAK_MAX_COUNT;
		}
	}
	yak_ct_destroy(g);
	a->ct->h[i].h = 0;
}

ha_pt_t *ha_pt_gen_count(ha_ct_t *ct, int n_thread)
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
	kt_for(n_thread, worker_pt_gen_count, &a, 1<<pt->pre);
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
		assert((a[j].x&mask) == (a[0].x&mask));
		k = yak_pt_get(g->h, x<<YAK_COUNTER_BITS);
		if (k == kh_end(g->h)) continue; // TODO: understand why we sometimes come here
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


int ha_pt_cnt_insert_list(ha_pt_t *h, int n, const uint64_t *a)
{
	int j, mask = (1<<h->pre) - 1, n_ins = 0;
	ha_pt1_t *g;
	if (n == 0) return 0;
	g = &h->h[a[0]&mask];
	for (j = 0; j < n; ++j) {
		uint64_t x = a[j] >> h->pre;
		khint_t k;
		assert((a[j]&mask) == (a[0]&mask));
		k = yak_pt_get(g->h, x<<YAK_COUNTER_BITS);
		if (k == kh_end(g->h)) continue; // TODO: understand why we sometimes come here
		++kh_val(g->h, k);
		++n_ins;
	}
	// fprintf(stderr, "n: %d, n_ins: %d\n", n, n_ins);
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
		if(h->h[i].a){
			free(h->h[i].a); h->h[i].a = NULL;
		}
		if(h->h[i].al){
			free(h->h[i].al); h->h[i].al = NULL;
		}	
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

const ha_idxposl_t *ha_ptl_get(const ha_pt_t *h, uint64_t hash, int *n)
{
	khint_t k;
	const ha_pt1_t *g = &h->h[hash & ((1ULL<<h->pre) - 1)];
	*n = 0;
	k = yak_pt_get(g->h, hash >> h->pre << YAK_COUNTER_BITS);
	if (k == kh_end(g->h)) return 0;
	*n = kh_key(g->h, k) & YAK_MAX_COUNT;
	return &g->al[kh_val(g->h, k)];
}

const int ha_pt_cnt(const ha_pt_t *h, uint64_t hash)
{
	khint_t k;
	const ha_pt1_t *g = &h->h[hash & ((1ULL<<h->pre) - 1)];
	k = yak_pt_get(g->h, hash >> h->pre << YAK_COUNTER_BITS);
	if (k == kh_end(g->h)) return 0;
	return kh_key(g->h, k) & YAK_MAX_COUNT;
}

inline uint64_t flt_quals(char *sc_a, uint64_t sc_l, uint64_t sc_off, int64_t sc_cut)
{
	int64_t sc_min = sc_l * sc_cut, sc_tot; uint64_t k;
	for (k = sc_tot = 0; (k < sc_l) && (sc_tot < sc_min); k++) {
		sc_tot += (((uint8_t)sc_a[k]) - sc_off);
	}

	// if(sc_tot < sc_min) {
	// 	fprintf(stderr, "[M::%s] sc_tot::%ld, sc_min::%ld, sc_l::%lu\n", __func__, sc_tot, sc_min, sc_l);
	// }

	if(sc_tot < sc_min) return 0;
	return 1;
}

/**********************************
 * Buffer for counting all k-mers *
 **********************************/
KSEQ_INIT(gzFile, gzread)
#define HAF_COUNT_EXACT  0x1
#define HAF_COUNT_ALL    0x2
#define HAF_RS_WRITE_LEN 0x4
#define HAF_RS_WRITE_SEQ 0x8
#define HAF_RS_READ      0x10
#define HAF_CREATE_NEW   0x20
#define HAF_SKIP_READ    0x40
#define HAF_UG_READ      0x80
#define HAF_COUNT_REFINE 0x100
#define HAF_RS_STASH     0x200 ///with HAF_RS_WRITE_LEN: keep the compressed reads and the names (ha_rss)

typedef struct { // global data structure for kt_pipeline()
	const yak_copt_t *opt;
	const void *flt_tab;
	int flag, create_new, is_store, uq;
	uint64_t n_mz, n_seq; ///number of total reads
	kseq_t *ks;
	UC_Read ucr;
	ha_ct_t *ct;
	ha_pt_t *pt;
	const All_reads *rs_in;
	All_reads *rs_out;
	const ma_utg_v *us_in;
	pthread_mutex_t bin_mutex;
	int n_bin_spare;
	void *bin_spare[3];
} pl_data_t;

/**parallel k-mer extraction for w == 1 (all k-mers): the reads of a chunk are split into nb
 * contiguous blocks; each block hashes its k-mers into its own list and counts them per bin; then
 * the k-mers are scattered to the bins. Within each bin, the k-mers keep the order of the serial
 * loop (block 0 first, then block 1, ...), so the counts with the bloom filter are the same.**/
typedef struct {
	int nb, n_pre, k, pre, is_HPC;
	const int *len; char **seq;
	int *bs;///block b holds the reads [bs[b], bs[b+1])
	uint64_t **a; int64_t *n, *m;///the k-mers of each block
	int64_t *cnt;///nb*n_pre: the count of block b in bin i; then the offset of block b in bin i
	uint64_t **dst;///the bin buffers
} ha_kblk_t;

///the k-mer list of one block in local variables: the n[] and m[] of adjacent blocks share cache lines
///(a store for each k-mer to t->n[b] was false sharing between the threads)
typedef struct {
	uint64_t *a; int64_t n, m, *cnt; uint64_t pmask;
} ha_kblk_buf_t;

static inline void ha_kblk_push(ha_kblk_buf_t *u, uint64_t y)
{
	if (u->n == u->m) {
		u->m = u->m < 1024? 1024 : u->m + (u->m>>1);
		REALLOC(u->a, u->m);
	}
	u->a[u->n++] = y;
	++u->cnt[y & u->pmask];
}

static void ha_kblk_worker_extract(void *data, long b, int tid)///callback for kt_for()
{
	ha_kblk_t *t = (ha_kblk_t*)data;
	int r, i, l, last, k = t->k;
	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;
	ha_kblk_buf_t u;
	u.a = t->a[b]; u.n = t->n[b]; u.m = t->m[b]; u.cnt = t->cnt + (int64_t)b*t->n_pre; u.pmask = (1<<t->pre) - 1;
	for (r = t->bs[b]; r < t->bs[b+1]; ++r) {
		const char *seq = t->seq[r]; int len = t->len[r];
		for (i = l = 0, last = -1, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i) {
			int c = seq_nt4_table[(uint8_t)seq[i]];
			if (c < 4) {
				if (!t->is_HPC || c != last) {
					x[0] = (x[0] << 1 | (c&1))  & mask;
					x[1] = (x[1] << 1 | (c>>1)) & mask;
					x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
					x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
					if (++l >= k) ha_kblk_push(&u, yak_hash_long(x));
					last = c;
				}
			} else l = 0, last = -1, x[0] = x[1] = x[2] = x[3] = 0;
		}
	}
	t->a[b] = u.a; t->n[b] = u.n; t->m[b] = u.m;
}

static void ha_kblk_worker_scatter(void *data, long b, int tid)///callback for kt_for()
{
	ha_kblk_t *t = (ha_kblk_t*)data;
	int64_t j, *off = t->cnt + (int64_t)b*t->n_pre;
	uint64_t pmask = (1<<t->pre) - 1;
	for (j = 0; j < t->n[b]; ++j) {
		uint64_t y = t->a[b][j];
		t->dst[y&pmask][off[y&pmask]++] = y;
	}
}

///extract the k-mers of the reads [0, n_seq) in parallel; ha_kblk_scatter() writes them to the bins
static void ha_kblk_extract(ha_kblk_t *t, int n_thread, int n_seq, const int *len, char **seq, int k, int pre, int is_HPC)
{
	int b, r; int64_t tot = 0, acc = 0;
	memset(t, 0, sizeof(*t));
	///2 blocks per thread: the inserts of step 2 run on the same cores, and a thread that waits for a core
	///must not hold back the whole step
	t->nb = 2*n_thread < n_seq? 2*n_thread : n_seq; if (t->nb < 1) t->nb = 1;
	t->n_pre = 1<<pre; t->k = k; t->pre = pre; t->is_HPC = is_HPC; t->len = len; t->seq = seq;
	CALLOC(t->bs, t->nb + 1); CALLOC(t->a, t->nb); CALLOC(t->n, t->nb); CALLOC(t->m, t->nb);
	CALLOC(t->cnt, (int64_t)t->nb * t->n_pre);
	for (r = 0; r < n_seq; ++r) tot += len[r];
	///blocks of about the same count of bases
	for (r = 0, b = 1; r < n_seq && b < t->nb; ++r) {
		acc += len[r];
		if (acc * t->nb >= tot * b) t->bs[b++] = r + 1;
	}
	for (; b <= t->nb; ++b) t->bs[b] = n_seq;
	kt_for(n_thread, ha_kblk_worker_extract, t, t->nb);
}

///bin_n[i]: the count of k-mers in bin i; then the bins dst[i] must hold bin_n[i] k-mers
static void ha_kblk_bin_count(ha_kblk_t *t, int64_t *bin_n)
{
	int b, i;
	///block by block (the rows of cnt in memory order); bin_n[i] keeps the sum of the blocks before b
	for (i = 0; i < t->n_pre; ++i) bin_n[i] = 0;
	for (b = 0; b < t->nb; ++b) {
		int64_t *c = t->cnt + (int64_t)b*t->n_pre;
		for (i = 0; i < t->n_pre; ++i) {
			int64_t x = c[i];
			c[i] = bin_n[i]; bin_n[i] += x;
		}
	}
}

static void ha_kblk_scatter(ha_kblk_t *t, int n_thread, uint64_t **dst)
{
	int b;
	t->dst = dst;
	kt_for(n_thread, ha_kblk_worker_scatter, t, t->nb);
	for (b = 0; b < t->nb; ++b) free(t->a[b]);
	free(t->a); free(t->n); free(t->m); free(t->cnt); free(t->bs);
}

/**the minimizers of one read in the cache: n minimizers at b, 11 bytes each (the hash x; pos minus the
 * pos before, 15 bits, and rev, 1 bit; span) with the rid of the read, or 16 bytes each (raw: ha_mz1_t)
 * if a pos difference does not fit or the rids differ.**/
typedef struct {
	uint32_t n, rid:31, raw:1;
	uint8_t *b;
} ha_mzk_t;

///a mapping that holds the packed minimizers of reads [.., r_end): its own mmap(), so a free gives the
///memory back to the system at once (the freed small malloc blocks of 5.5 M reads stay in the arenas)
typedef struct {
	uint8_t *p;
	uint64_t sz, used, r_end;
} ha_mzb_t;

#define HA_MZK_PACKED 11
#define HA_MZB_SIZE (64ULL<<20)

/**the minimizers of the first ha_count() pass of ha_pt_gen() (the count), kept per read in input
 * order and reused by the second pass (the index): both passes sketch the same reads with the same
 * options, so the second pass does not decode and sketch them again. mode: 0 off, 1 record, 2 replay.
 * a holds the packed minimizers of each read (ha_mzk_t), in the mappings blk[0..n_blk).**/
typedef struct {
	int mode, replayed;
	int query;///1: the cache holds the minimizers of each stored read (ids 0..n-1), for ha_mzc_query()
	uint64_t n, m, pos;
	ha_mzk_t *a;
	uint64_t n_blk, m_blk, blk_free;///blk[0..blk_free) are unmapped
	ha_mzb_t *blk;
	uint64_t f_m; ha_mz1_t *f;///the decoded minimizers of the chunk of a replay (one array for all chunks)
} ha_mzc_t;
static ha_mzc_t ha_mzc;

///the bytes of the packed minimizers of v, and the rid and form of k
static uint64_t ha_mzk_size(const ha_mz1_v *v, ha_mzk_t *k)
{
	uint64_t j, pos = 0;
	k->raw = 0; k->rid = v->n? v->a[0].rid : 0;
	for (j = 0; j < v->n; j++) {
		if (v->a[j].pos < pos || v->a[j].pos - pos >= (1U<<15) || v->a[j].rid != k->rid) {
			k->raw = 1;
			return (uint64_t)v->n * sizeof(ha_mz1_t);
		}
		pos = v->a[j].pos;
	}
	return (uint64_t)v->n * HA_MZK_PACKED;
}

static void ha_mzk_put(uint8_t *p, const ha_mz1_v *v, uint32_t raw)
{
	uint64_t j, pos = 0;
	if (raw) {
		memcpy(p, v->a, (uint64_t)v->n * sizeof(ha_mz1_t));
		return;
	}
	for (j = 0; j < v->n; j++, p += HA_MZK_PACKED) {
		const ha_mz1_t *z = &v->a[j];
		uint16_t d = (uint16_t)((z->pos - pos) << 1 | z->rev);
		memcpy(p, &z->x, 8); memcpy(p + 8, &d, 2); p[10] = z->span;
		pos = z->pos;
	}
}

///the n minimizers of k into o; a query (rid 0) if is_q
static void ha_mzk_get(const ha_mzk_t *k, int is_q, ha_mz1_t *o)
{
	const uint8_t *p = k->b;
	uint64_t j, pos = 0, rid = is_q? 0 : k->rid;
	if (k->raw) {
		memcpy(o, p, (uint64_t)k->n * sizeof(ha_mz1_t));
		if (is_q) for (j = 0; j < k->n; j++) o[j].rid = 0;
		return;
	}
	for (j = 0; j < k->n; j++, p += HA_MZK_PACKED) {
		uint64_t x; uint16_t d;
		memcpy(&x, p, 8); memcpy(&d, p + 8, 2);
		pos += d >> 1;
		o[j].x = x; o[j].rid = rid; o[j].pos = pos; o[j].rev = d & 1; o[j].span = p[10];
	}
}

/**the free mappings of HA_MZB_SIZE bytes: the next cache takes them (its pages are there: no munmap() of
 * 27 GB and no new page faults at each index rebuild). ha_mzc_release() unmaps them in a thread.**/
static uint8_t **ha_mzb_pool; static uint64_t ha_mzb_pool_n, ha_mzb_pool_m;

/**the record of a chunk (ha_count() step 1, mode 1) has three parts: ha_mzc_grow() gives the entries of its
 * n reads; the first kt_for() of sf##_mzbin_dist() sets each entry and its size (ha_mzk_size()); then
 * ha_mzc_reserve() gives the bytes and the second kt_for() packs each read (ha_mzk_put()); ha_mzc_commit().
 * Thus the cache adds no kt_for() to the pipeline (a kt_for() for each chunk of 20 Mb wakes all threads).**/
static ha_mzk_t *ha_mzc_grow(ha_mzc_t *z, int64_t n)
{
	if (z->n + n > z->m) {
		z->m = z->n + n; z->m += z->m>>1;
		REALLOC(z->a, z->m);
	}
	return z->a + z->n;
}

///tot bytes in the last mapping of z (a new mapping if they do not fit)
static uint8_t *ha_mzc_reserve(ha_mzc_t *z, uint64_t tot)
{
	ha_mzb_t *b = z->n_blk? &z->blk[z->n_blk - 1] : NULL;
	uint8_t *p;
	if (tot == 0) return NULL;
	if (b == NULL || b->p == NULL || b->used + tot > b->sz) {
		if (z->n_blk == z->m_blk) {
			z->m_blk = z->m_blk? z->m_blk<<1 : 64;
			REALLOC(z->blk, z->m_blk);
		}
		b = &z->blk[z->n_blk++];
		b->sz = tot > HA_MZB_SIZE? tot : HA_MZB_SIZE; b->used = 0; b->r_end = z->n;
		if (b->sz == HA_MZB_SIZE && ha_mzb_pool_n) b->p = ha_mzb_pool[--ha_mzb_pool_n];
		else b->p = (uint8_t*)mmap(NULL, b->sz, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
		if (b->p == MAP_FAILED) {
			fprintf(stderr, "[E::%s] mmap of %lu bytes failed\n", __func__, (unsigned long)b->sz);
			exit(1);
		}
	}
	p = b->p + b->used; b->used += tot;
	return p;
}

static void ha_mzc_commit(ha_mzc_t *z, int64_t n)
{
	z->n += n;
	if (z->n_blk) z->blk[z->n_blk - 1].r_end = z->n;
}

///the arrays of reads [r0, r0+n) of z: mz[i].a points into z->f, one array for the chunks of the replay (a
///new array for each chunk of 20 Mb was a new mmap() with new pages); the first kt_for() of
///sf##_mzbin_dist() decodes each read (ha_mzk_get())
static void ha_mzc_get(ha_mzc_t *z, uint64_t r0, int64_t n, ha_mz1_v *mz)
{
	int64_t i; uint64_t tot;
	for (i = 0, tot = 0; i < n; i++) tot += z->a[r0 + i].n;
	if (tot > z->f_m) {
		z->f_m = tot + (tot>>2);
		free(z->f); MALLOC(z->f, z->f_m);
	}
	for (i = 0, tot = 0; i < n; i++) {
		mz[i].a = z->f + tot; mz[i].n = mz[i].m = z->a[r0 + i].n;
		tot += mz[i].n;
	}
}

static pthread_t ha_mzb_unmap_t; static int ha_mzb_unmap_on;
typedef struct { uint8_t **p; uint64_t n; } ha_mzb_list_t;

static void *ha_mzb_unmap_run(void *data)
{
	ha_mzb_list_t *l = (ha_mzb_list_t*)data; uint64_t i;
	for (i = 0; i < l->n; i++) munmap(l->p[i], HA_MZB_SIZE);
	free(l->p); free(l);
	return NULL;
}

static void ha_mzb_unmap_join(void)
{
	if (ha_mzb_unmap_on) pthread_join(ha_mzb_unmap_t, NULL);
	ha_mzb_unmap_on = 0;
}

///unmap the pool in a thread
static void ha_mzb_pool_unmap(void)
{
	static int exit_set = 0;
	ha_mzb_list_t *l;
	ha_mzb_unmap_join();
	if (ha_mzb_pool_n == 0) return;
	if (!exit_set) {
		atexit(ha_mzb_unmap_join); exit_set = 1;
	}
	MALLOC(l, 1); l->p = ha_mzb_pool; l->n = ha_mzb_pool_n;
	ha_mzb_pool = NULL; ha_mzb_pool_n = ha_mzb_pool_m = 0;
	if (pthread_create(&ha_mzb_unmap_t, NULL, ha_mzb_unmap_run, l) == 0) ha_mzb_unmap_on = 1;
	else ha_mzb_unmap_run(l);
}

///the mappings of z that hold only reads before r go to the pool (or are unmapped if larger)
static void ha_mzc_free_before(ha_mzc_t *z, uint64_t r)
{
	for (; z->blk_free < z->n_blk && z->blk[z->blk_free].r_end <= r; z->blk_free++) {
		ha_mzb_t *b = &z->blk[z->blk_free];
		if (b->p && b->sz == HA_MZB_SIZE) {
			if (ha_mzb_pool_n == ha_mzb_pool_m) {
				ha_mzb_pool_m = ha_mzb_pool_m? ha_mzb_pool_m<<1 : 64;
				REALLOC(ha_mzb_pool, ha_mzb_pool_m);
			}
			ha_mzb_pool[ha_mzb_pool_n++] = b->p;
		} else if (b->p) munmap(b->p, b->sz);
		b->p = NULL;
	}
}

/**the reads of the filter count pass (ha_ft_gen()), compressed as R_INF keeps them (ha_compress_base()),
 * and their names: the first ha_pt_gen() adopts them (ha_rss_adopt()) and does not parse the input
 * files a second time. seq[i] and N[i] are the read_sperate[i] and N_site[i] of read i.**/
typedef struct {
	uint64_t n, m, name_l, name_m;
	uint8_t **seq; uint64_t **N; char *name;
} ha_rss_t;
static ha_rss_t ha_rss;

static void ha_rss_free(void)
{
	uint64_t i;
	for (i = 0; i < ha_rss.n; i++) {
		free(ha_rss.seq[i]); free(ha_rss.N[i]);
	}
	free(ha_rss.seq); free(ha_rss.N); free(ha_rss.name);
	memset(&ha_rss, 0, sizeof(ha_rss));
}

static void ha_rss_add_name(const char *name, uint64_t l)
{
	if (ha_rss.name_l + l > ha_rss.name_m) {
		ha_rss.name_m = ha_rss.name_l + l; ha_rss.name_m += ha_rss.name_m>>1;
		REALLOC(ha_rss.name, ha_rss.name_m);
	}
	memcpy(ha_rss.name + ha_rss.name_l, name, l); ha_rss.name_l += l;
}

///the reads [n0, n0+n) of a chunk, in the order of the chunks
static void ha_rss_add_reads(uint64_t n0, int n, uint8_t **seq, uint64_t **N)
{
	assert(ha_rss.n == n0);
	if (ha_rss.n + n > ha_rss.m) {
		ha_rss.m = ha_rss.n + n; ha_rss.m += ha_rss.m>>1;
		REALLOC(ha_rss.seq, ha_rss.m); REALLOC(ha_rss.N, ha_rss.m);
	}
	memcpy(ha_rss.seq + ha_rss.n, seq, n * sizeof(*seq));
	memcpy(ha_rss.N + ha_rss.n, N, n * sizeof(*N));
	ha_rss.n += n;
}

///give the kept reads to rs (its lengths are loaded, its reads are not): the same bytes as a second parse
///with HAF_RS_WRITE_SEQ. Returns 0 (and frees them) if they do not match rs.
static int ha_rss_adopt(All_reads *rs)
{
	uint64_t i;
	if (ha_rss.seq == NULL || ha_rss.n != rs->total_reads || ha_rss.name_l != rs->total_name_length || rs->rsc) {
		ha_rss_free();
		return 0;
	}
	malloc_All_reads2(rs, 0);///no buffer for each read: the kept reads take their place
	for (i = 0; i < rs->total_reads; i++) {
		rs->read_sperate[i] = ha_rss.seq[i];
		rs->N_site[i] = ha_rss.N[i];
	}
	memcpy(rs->name, ha_rss.name, ha_rss.name_l);
	free(ha_rss.seq); free(ha_rss.N); free(ha_rss.name);
	memset(&ha_rss, 0, sizeof(ha_rss));
	return 1;
}

static void ha_mzv_free_all(ha_mzc_t *z)
{
	ha_mzc_free_before(z, UINT64_MAX);
	free(z->a); free(z->blk); free(z->f); memset(z, 0, sizeof(*z));
}

static void ha_mzc_free(void)
{
	ha_mzv_free_all(&ha_mzc);
}

/**the reads that changed since the cache of the queries (ha_mzc.query) was made: each writer of a
 * stored read marks it (ha_rd_mark(): the correction, the strand flip). The next ha_pt_gen() moves the
 * cache to ha_mzp, and its count pass takes the cached minimizers of each read that did not change:
 * the sketch of the count pass depends only on the read, its id and the options (no index, no
 * k_flag), and the cache holds the minimizers of that pass.**/
static uint8_t *ha_rd_dirty; static uint64_t ha_rd_n;
static ha_mzc_t ha_mzp; static uint8_t *ha_mzp_dirty;

void ha_rd_mark(uint64_t rid)
{
	if (ha_rd_dirty && rid < ha_rd_n) ha_rd_dirty[rid] = 1;
}

static void ha_mzp_free(void)
{
	ha_mzv_free_all(&ha_mzp);
	free(ha_mzp_dirty); ha_mzp_dirty = NULL;
}

///1 if the count pass can take the cached minimizers of read rid
static inline int ha_mzp_hit(uint64_t rid)
{
	return ha_mzp.a && rid < ha_mzp.n && !ha_mzp_dirty[rid];
}

///if all reads changed, free the cache before the count pass, which sketches them all again: else it stays
///next to the new cache (after rounds 1 and 2 all reads change). Else the count pass unmaps each mapping
///of ha_mzp when it is past its reads.
static void ha_mzp_drop_dirty(void)
{
	uint64_t i;
	for (i = 0; i < ha_mzp.n && ha_mzp_dirty[i]; i++);
	if (i == ha_mzp.n) ha_mzp_free();
}

void ha_mzc_release(void)
{
	ha_mzc_free();
	free(ha_rd_dirty); ha_rd_dirty = NULL; ha_rd_n = 0;
	ha_mzb_pool_unmap();
}

#ifdef HA_MZR_CHECK
static uint64_t ha_mzr_n_check = 0, ha_mzr_n_bad = 0;
static void ha_mzr_report(void) { fprintf(stderr, "[M::ha_mzr_check] cached minimizers: %lu reads; different: %lu\n", (unsigned long)ha_mzr_n_check, (unsigned long)ha_mzr_n_bad); }
static int ha_mzr_report_set = atexit(ha_mzr_report);
#endif

///the minimizers of stored read rid, as mz1_ha_sketch() gives them for a query (rid 0, no k_flag, no
///dbg_ct, no refine): from the index of this round, if it holds the read with length rl
int ha_mzc_query(uint64_t rid, uint64_t rl, const All_reads *rs, ha_mz1_v *out)
{
	const ha_mzk_t *c;
	if (!ha_mzc.query || rid >= ha_mzc.n || rl != Get_READ_LENGTH((*rs), rid)) return 0;
	c = &ha_mzc.a[rid];
	if (out->m < c->n) {
		out->m = c->n; REALLOC(out->a, out->m);
	}
	out->n = c->n;
	ha_mzk_get(c, 1, out->a);
	return 1;
}

#ifdef HA_MZR_CHECK ///decode and sketch the cached reads too, and compare
#define HA_MZR_DECODE 1
#define HA_MZR_CHECK_READ do {\
		mz1_ha_sketch(s->seq[i], s->len[i], s->p->opt->w, s->p->opt->k, s->n_seq0 + i, s->p->opt->is_HPC, (ha_mz1_v*)b, s->p->flt_tab, asm_opt.mz_sample_dist, 0, 0, \
		(s->p->pt&&(s->p->flag&HAF_COUNT_REFINE))?s->p->pt:NULL, s->p->opt->min_rcnt, asm_opt.dp_min_len, asm_opt.dp_e, &mtl, asm_opt.mz_rewin, s->uq, NULL);\
		__sync_fetch_and_add(&ha_mzr_n_check, 1);\
		if (b->n != s->mz[i].n || (b->n && memcmp(b->a, s->mz[i].a, b->n * sizeof(*b->a)) != 0)) __sync_fetch_and_add(&ha_mzr_n_bad, 1);\
	} while (0)
#else
#define HA_MZR_DECODE 0
#define HA_MZR_CHECK_READ do {} while (0)
#endif

#ifdef HA_BINREUSE_CHECK
static uint64_t ha_binreuse_chunks, ha_binreuse_elements, ha_binreuse_bad, ha_binreuse_reused;
static void ha_binreuse_report(void)
{
	fprintf(stderr, "[M::HA_BINREUSE_CHECK] chunks=%lu reused=%lu elements=%lu different_bins=%lu\n", (unsigned long)ha_binreuse_chunks, (unsigned long)ha_binreuse_reused, (unsigned long)ha_binreuse_elements, (unsigned long)ha_binreuse_bad);
}
#define HA_BINREUSE_USED() __sync_fetch_and_add(&ha_binreuse_reused, 1)
#define HA_BINREUSE_VERIFY(sf, s, is_pt) do {\
	sf##_ch_buf_t *old; int r, bi, np = 1<<s->p->opt->pre; uint32_t mj; uint64_t ne = 0, bad = 0;\
	CALLOC(old, np);\
	for (r = 0; r < s->n_seq; ++r) {\
		if (s->p->opt->w == 1) {\
			if (s->p->opt->is_HPC) sf##_count_seq_buf_HPC(old, s->p->opt->k, s->p->opt->pre, s->len[r], s->seq[r]);\
			else sf##_count_seq_buf(old, s->p->opt->k, s->p->opt->pre, s->len[r], s->seq[r]);\
		} else for (mj = 0; mj < s->mz[r].n; ++mj) {\
			if (is_pt) sf##_pt_insert_buf(old, s->p->opt->pre, &s->mz[r].a[mj]);\
			else sf##_ct_insert_buf(old, s->p->opt->pre, s->mz[r].a[mj].x);\
		}\
	}\
	for (bi = 0; bi < np; ++bi) {\
		sf##_ch_buf_t *c = &old[bi], *b = &s->buf[bi];\
		ne += c->n;\
		if (c->n != b->n || (c->n && (is_pt? memcmp(c->b, b->b, c->n * sizeof(*c->b)) : memcmp(c->a, b->a, c->n * sizeof(*c->a))))) ++bad;\
		free(c->a); free(c->b);\
	}\
	free(old);\
	__sync_fetch_and_add(&ha_binreuse_chunks, 1); __sync_fetch_and_add(&ha_binreuse_elements, ne); __sync_fetch_and_add(&ha_binreuse_bad, bad);\
} while (0)
#else
#define HA_BINREUSE_VERIFY(sf, s, is_pt) do {} while (0)
#define HA_BINREUSE_USED() do {} while (0)
#endif

#define MZ_TEST_INIT(sf, HType, VType, IType, Ia) \
typedef struct {int n, m; uint64_t n_ins; uint64_t *a; HType *b;} sf##_ch_buf_t;\
/**a buffer belongs to one chunk until its insert finishes; the pipeline has at most three chunks**/\
static sf##_ch_buf_t *sf##_bin_take(pl_data_t *p, int n_pre)\
{\
	sf##_ch_buf_t *b = NULL; int i;\
	pthread_mutex_lock(&p->bin_mutex);\
	if (p->n_bin_spare) b = (sf##_ch_buf_t*)p->bin_spare[--p->n_bin_spare];\
	pthread_mutex_unlock(&p->bin_mutex);\
	if (b == NULL) CALLOC(b, n_pre);\
	else {\
		HA_BINREUSE_USED();\
		for (i = 0; i < n_pre; ++i) b[i].n = 0, b[i].n_ins = 0;\
	}\
	return b;\
}\
static void sf##_bin_free(pl_data_t *p)\
{\
	int i, j, n_pre = 1<<p->opt->pre;\
	for (j = 0; j < p->n_bin_spare; ++j) {\
		sf##_ch_buf_t *b = (sf##_ch_buf_t*)p->bin_spare[j];\
		for (i = 0; i < n_pre; ++i) free(b[i].a), free(b[i].b);\
		free(b);\
	}\
}\
static inline void sf##_ct_insert_buf(sf##_ch_buf_t *buf, int p, uint64_t y) /** insert a k-mer $y to a linear buffer**/\
{\
	/**assign k-mer to one of the 4096 bins**/\
	/**using low 12 bits for assigning**/\
	/**so all elements at b have the same low 12 bits**/\
	int pre = y & ((1<<p) - 1);\
	sf##_ch_buf_t *b = &buf[pre];\
	if (b->n == b->m) {\
		b->m = b->m < 8? 8 : b->m + (b->m>>1);\
		REALLOC(b->a, b->m);\
	}\
	b->a[b->n++] = y;\
}\
/**buf is the read block, k is the k-mer length, p = 12, len is the read length, seq is the read**/\
static void sf##_count_seq_buf(sf##_ch_buf_t *buf, int k, int p, int len, const char *seq) /**insert k-mers in $seq to linear buffer $buf**/\
{\
	int i, l;\
	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;\
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i) {\
		int c = seq_nt4_table[(uint8_t)seq[i]];\
		/**c = 00, 01, 10, 11**/\
		if (c < 4) { /** not an "N" base**/\
			/**x[0] & x[1] are the forward k-mer**/\
			/**x[2] & x[3] are the reverse complementary k-mer**/\
			x[0] = (x[0] << 1 | (c&1))  & mask;\
			x[1] = (x[1] << 1 | (c>>1)) & mask;\
			x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;\
			x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;\
			if (++l >= k)\
				sf##_ct_insert_buf(buf, p, yak_hash_long(x));\
		} else l = 0, x[0] = x[1] = x[2] = x[3] = 0; /** if there is an "N", restart**/\
	}\
}\
static void sf##_count_seq_buf_HPC(sf##_ch_buf_t *buf, int k, int p, int len, const char *seq) /**insert k-mers in $seq to linear buffer $buf**/\
{\
	int i, l, last = -1;\
	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;\
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i) {\
		int c = seq_nt4_table[(uint8_t)seq[i]];\
		if (c < 4) { /** not an "N" base**/\
			if (c != last) {\
				x[0] = (x[0] << 1 | (c&1))  & mask;\
				x[1] = (x[1] << 1 | (c>>1)) & mask;\
				x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;\
				x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;\
				if (++l >= k)\
					sf##_ct_insert_buf(buf, p, yak_hash_long(x));\
				last = c;\
			}\
		} else l = 0, last = -1, x[0] = x[1] = x[2] = x[3] = 0; /**if there is an "N", restart**/\
	}\
}\
int sf##_ha_pt_insert_list(ha_pt_t *h, int n, const HType *a)\
{\
	int j, mask = (1<<h->pre) - 1, n_ins = 0;\
	ha_pt1_t *g;\
	if (n == 0) return 0;\
	g = &h->h[a[0].x&mask];\
	/**fprintf(stderr, "a[0].x&mask: %lu, a[0].x&mask: %lu, n: %d\n", a[0].x&mask, a[0].x, n);**/\
	for (j = 0; j < n; ++j) {\
		uint64_t x = a[j].x >> h->pre;\
		khint_t k;\
		int n;\
		IType *p;\
		assert((a[j].x&mask) == (a[0].x&mask));\
		k = yak_pt_get(g->h, x<<YAK_COUNTER_BITS);\
		if (k == kh_end(g->h)) continue; \
		n = kh_key(g->h, k) & YAK_MAX_COUNT;\
		/**fprintf(stderr, "j: %d, n: %d\n", j, n);**/\
		assert(n < YAK_MAX_COUNT);\
		p = &g->Ia[kh_val(g->h, k) + n];\
		p->rid = a[j].rid, p->rev = a[j].rev, p->pos = a[j].pos, p->span = a[j].span;\
		/**(uint64_t)a[j].rid<<36 | (uint64_t)a[j].rev<<35 | (uint64_t)a[j].pos<<8 | (uint64_t)a[j].span;**/\
		++kh_key(g->h, k);\
		++n_ins;\
	}\
	return n_ins;\
}\
/** data structure for each step in kt_pipeline()**/\
typedef struct {pl_data_t *p;uint64_t n_seq0; int n_seq, m_seq, sum_len, nk, uq, *len; char **seq; VType *mz_buf; VType *mz;sf##_ch_buf_t *buf;st_mt_t *mt; uint8_t **rss_seq; uint64_t **rss_N; uint8_t *hit; void *mzf; const ha_mzk_t *kdec; ha_mzk_t *kenc; uint64_t *koff; uint8_t *kb;} sf##_st_data_t;\
static void sf##_worker_for_insert(void *data, long i, int tid) /** callback for kt_for()**/\
{\
	sf##_st_data_t *s = (sf##_st_data_t*)data;\
	sf##_ch_buf_t *b = &s->buf[i];\
	if (s->p->pt){\
		if(s->p->flag&HAF_COUNT_REFINE) b->n_ins += ha_pt_cnt_insert_list(s->p->pt, b->n, b->a);\
		else b->n_ins += sf##_ha_pt_insert_list(s->p->pt, b->n, b->b);\
	}else{\
		b->n_ins += ha_ct_insert_list(s->p->ct, s->p->create_new, b->n, b->a);\
	}\
}\
static void sf##_worker_for_mz(void *data, long i, int tid)\
{\
	sf##_st_data_t *s = (sf##_st_data_t*)data;\
	/**get the corresponding minimzer vector of this read; the buffers of the thread are local here: the\
	 * mz_buf[] and mt[] of adjacent threads share cache lines (a push to them was false sharing)**/\
	VType bl = s->mz_buf[tid], *b = &bl; st_mt_t mtl = s->mt[tid];\
	bl.n = 0;\
	if (s->hit && s->hit[i]) { /**the cached minimizers of an unchanged read (ha_mzp)**/\
		const ha_mzk_t *c = &ha_mzp.a[s->n_seq0 + i];\
		s->mz[i].n = s->mz[i].m = c->n;\
		MALLOC(s->mz[i].a, c->n);\
		ha_mzk_get(c, 0, (ha_mz1_t*)s->mz[i].a);\
		HA_MZR_CHECK_READ;\
		s->mz_buf[tid] = bl; s->mt[tid] = mtl;\
		return;\
	}\
	sf##_ha_sketch(s->seq[i], s->len[i], s->p->opt->w, s->p->opt->k, s->n_seq0 + i, s->p->opt->is_HPC, b, s->p->flt_tab, asm_opt.mz_sample_dist, 0, 0, \
	(s->p->pt&&(s->p->flag&HAF_COUNT_REFINE))?s->p->pt:NULL, s->p->opt->min_rcnt, asm_opt.dp_min_len, asm_opt.dp_e, &mtl, asm_opt.mz_rewin, s->uq, NULL);\
	s->mz[i].n = s->mz[i].m = b->n;\
	MALLOC(s->mz[i].a, b->n);\
	MEMCPY(s->mz[i].a, b->a, b->n);\
	s->mz_buf[tid] = bl; s->mt[tid] = mtl;\
}\
static inline void sf##_pt_insert_buf(sf##_ch_buf_t *buf, int p, const HType *y){\
	/**assign minimizer to one of 4096 bins by low 12 bits**/\
	int pre = y->x & ((1<<p) - 1);\
	sf##_ch_buf_t *b = &buf[pre];\
	if (b->n == b->m) {\
		b->m = b->m < 8? 8 : b->m + (b->m>>1);\
		REALLOC(b->b, b->m);\
	}\
	b->b[b->n++] = *y;\
}\
/**parallel distribution of the minimizers of a chunk to the 4096 bins: the reads are split into nb\
 * contiguous blocks; each block counts its minimizers per bin, then writes them at its offset. In each\
 * bin, the minimizers keep the order of the serial loop (read order), so the tables are the same.**/\
typedef struct {sf##_st_data_t *s; int nb, n_pre, is_pt; int *bs; int64_t *cnt;} sf##_mzbin_t;\
static void sf##_mzbin_worker_count(void *data, long b, int tid)\
{\
	sf##_mzbin_t *t = (sf##_mzbin_t*)data;\
	int64_t *cnt = t->cnt + (int64_t)b*t->n_pre; uint64_t pmask = t->n_pre - 1;\
	int i; uint32_t j;\
	for (i = t->bs[b]; i < t->bs[b+1]; ++i) {\
		if (t->s->kdec) ha_mzk_get(&t->s->kdec[i], 0, (ha_mz1_t*)t->s->mz[i].a);/**a replay of the cache (ha_mzc)**/\
		for (j = 0; j < t->s->mz[i].n; ++j) ++cnt[t->s->mz[i].a[j].x & pmask];\
		if (t->s->kenc) t->s->koff[i] = ha_mzk_size((const ha_mz1_v*)&t->s->mz[i], &t->s->kenc[i]);/**a record**/\
	}\
}\
static void sf##_mzbin_worker_scatter(void *data, long b, int tid)\
{\
	sf##_mzbin_t *t = (sf##_mzbin_t*)data;\
	int64_t *off = t->cnt + (int64_t)b*t->n_pre; uint64_t pmask = t->n_pre - 1;\
	sf##_ch_buf_t *buf = t->s->buf;\
	int i; uint32_t j;\
	for (i = t->bs[b]; i < t->bs[b+1]; ++i) {\
		const HType *a = t->s->mz[i].a;\
		if (t->s->kenc) {/**pack read i into the cache (ha_mzc)**/\
			ha_mzk_t *k = &t->s->kenc[i];\
			k->n = t->s->mz[i].n; k->b = k->n? t->s->kb + t->s->koff[i] : NULL;\
			if (k->n) ha_mzk_put(k->b, (const ha_mz1_v*)&t->s->mz[i], k->raw);\
		}\
		if (t->is_pt) {\
			for (j = 0; j < t->s->mz[i].n; ++j) buf[a[j].x & pmask].b[off[a[j].x & pmask]++] = a[j];\
		} else {\
			for (j = 0; j < t->s->mz[i].n; ++j) buf[a[j].x & pmask].a[off[a[j].x & pmask]++] = a[j].x;\
		}\
	}\
}\
static void sf##_mzbin_dist(sf##_st_data_t *s, int n_thread, int pre, int is_pt)\
{\
	sf##_mzbin_t t; int b, i, trim; int64_t tot = 0, acc = 0, cap = 0;\
	t.s = s; t.n_pre = 1<<pre; t.is_pt = is_pt;\
	t.nb = n_thread < s->n_seq? n_thread : s->n_seq; if (t.nb < 1) t.nb = 1;\
	CALLOC(t.bs, t.nb + 1); CALLOC(t.cnt, (int64_t)t.nb * t.n_pre);\
	for (i = 0; i < s->n_seq; ++i) tot += s->mz[i].n;\
	/**blocks of about the same count of minimizers**/\
	for (i = 0, b = 1; i < s->n_seq && b < t.nb; ++i) {\
		acc += s->mz[i].n;\
		if (acc * t.nb >= tot * b) t.bs[b++] = i + 1;\
	}\
	for (; b <= t.nb; ++b) t.bs[b] = s->n_seq;\
	kt_for(n_thread, sf##_mzbin_worker_count, &t, t.nb);\
	if (s->kenc) {/**the offsets of the packed reads, and their bytes in the cache**/\
		uint64_t ktot = 0, ksz;\
		for (i = 0; i < s->n_seq; ++i) ksz = s->koff[i], s->koff[i] = ktot, ktot += ksz;\
		s->kb = ha_mzc_reserve(&ha_mzc, ktot);\
	}\
	{/**counts to offsets, after the elements already in the bin; block by block (the rows of cnt in\
	  * memory order): run[i] keeps the end of bin i so far**/\
		int64_t *run; MALLOC(run, t.n_pre);\
		for (i = 0; i < t.n_pre; ++i) run[i] = s->buf[i].n;\
		for (b = 0; b < t.nb; ++b) {\
			int64_t *cb = t.cnt + (int64_t)b*t.n_pre;\
			for (i = 0; i < t.n_pre; ++i) {\
				int64_t x = cb[i];\
				cb[i] = run[i]; run[i] += x;\
			}\
		}\
		for (i = 0; i < t.n_pre; ++i) s->buf[i].n = run[i], cap += s->buf[i].m;/**the new size, read below**/\
		free(run);\
	}\
	/**a sequence of skewed chunks must not retain a large buffer in every bin**/\
	trim = cap > 2*(tot + t.n_pre);\
	for (i = 0; i < t.n_pre; ++i) {\
		sf##_ch_buf_t *c = &s->buf[i];\
		int64_t n = c->n;\
		if (n > c->m || c->m == 0 || (trim && n < c->m)) {\
			int64_t m = trim? n : (int64_t)c->m + (c->m>>1);\
			if (m > INT_MAX) m = n;\
			c->m = n > m? n : m; if (c->m == 0) c->m = 1;\
			if (is_pt) REALLOC(c->b, c->m);\
			else REALLOC(c->a, c->m);\
		}\
		c->n = n;\
	}\
	kt_for(n_thread, sf##_mzbin_worker_scatter, &t, t.nb);\
	free(t.bs); free(t.cnt);\
}\
static void sf##_worker_for_rs(void *data, long i, int tid) /** compress read i to R_INF**/\
{\
	sf##_st_data_t *s = (sf##_st_data_t*)data;\
	All_reads *rs = s->p->rs_out;\
	uint64_t rid = s->n_seq0 + i;\
	int j, l = s->len[i], n_N;\
	for (j = n_N = 0; j < l; ++j) /** count number of ambiguous bases**/\
		if (seq_nt4_table[(uint8_t)s->seq[i][j]] >= 4)\
			++n_N;\
	ha_compress_base(Get_READ(*rs, rid), s->seq[i], l, &rs->N_site[rid], n_N);\
}\
static void sf##_worker_for_rss(void *data, long i, int tid) /** compress read i for ha_rss (HAF_RS_STASH)**/\
{\
	sf##_st_data_t *s = (sf##_st_data_t*)data;\
	int j, l = s->len[i], n_N;\
	for (j = n_N = 0; j < l; ++j)\
		if (seq_nt4_table[(uint8_t)s->seq[i][j]] >= 4)\
			++n_N;\
	MALLOC(s->rss_seq[i], l/4+1);\
	ha_compress_base(s->rss_seq[i], s->seq[i], l, &s->rss_N[i], n_N);\
}\
static void *sf##_worker_count(void *data, int step, void *in) /** callback for kt_pipeline()**/\
{\
	pl_data_t *p = (pl_data_t*)data;\
	/**uint8_t src_a[1000000], des_a[1000000];**/\
	if (step == 0) { /** step 1: read a block of sequences**/\
		int ret;\
		sf##_st_data_t *s;\
		CALLOC(s, 1);\
		s->p = p;\
		s->n_seq0 = p->n_seq;\
		s->uq = p->opt->uq;\
		if (p->rs_in && (p->flag & HAF_RS_READ)) {\
			while (p->n_seq < p->rs_in->total_reads) {\
				if ((p->flag & HAF_SKIP_READ) && p->rs_in->trio_flag[p->n_seq] != AMBIGU) {\
					++p->n_seq;\
					continue;\
				}\
				int l, hit = (ha_mzc.mode == 1 && ha_mzp_hit(p->n_seq));\
				if (ha_mzc.mode == 2 || (hit && !HA_MZR_DECODE)) {/**replay or a cached read (ha_mzp): the minimizers are kept, the read is not needed**/\
					l = Get_READ_LENGTH((*p->rs_in), p->n_seq);\
				} else {\
					recover_UC_Read(&p->ucr, p->rs_in, p->n_seq);\
					l = p->ucr.length;\
				}\
				if (s->n_seq == s->m_seq) {\
					s->m_seq = s->m_seq < 16? 16 : s->m_seq + (s->m_seq>>1);\
					REALLOC(s->len, s->m_seq);\
					REALLOC(s->seq, s->m_seq);\
					if (ha_mzp.a) REALLOC(s->hit, s->m_seq);\
				}\
				if (ha_mzp.a) s->hit[s->n_seq] = hit;\
				if (ha_mzc.mode == 2 || (hit && !HA_MZR_DECODE)) {\
					s->seq[s->n_seq] = NULL;\
				} else {\
					MALLOC(s->seq[s->n_seq], l);\
					memcpy(s->seq[s->n_seq], p->ucr.seq, l);\
				}\
				s->len[s->n_seq++] = l;\
				++p->n_seq;\
				s->sum_len += l;\
				s->nk += l >= p->opt->k? l - p->opt->k + 1 : 0;\
				if (s->sum_len >= p->opt->chunk_size)\
					break;\
			}\
		} else if(p->us_in) {\
			ma_utg_t *u;\
			while (p->n_seq < p->us_in->n) {\
				u = &(p->us_in->a[p->n_seq]);\
				if (s->n_seq == s->m_seq) {\
					s->m_seq = s->m_seq < 16? 16 : s->m_seq + (s->m_seq>>1);\
					REALLOC(s->len, s->m_seq);\
					REALLOC(s->seq, s->m_seq);\
				}\
				MALLOC(s->seq[s->n_seq], u->len);\
				if(u->s) memcpy(s->seq[s->n_seq], u->s, u->len);\
				else retrieve_u_seq(NULL, s->seq[s->n_seq], u, 0, 0, -1, NULL);\
				s->len[s->n_seq++] = u->len;\
				++p->n_seq;\
				s->sum_len += u->len;\
				s->nk += u->len >= p->opt->k? u->len - p->opt->k + 1 : 0;\
				if (s->sum_len >= p->opt->chunk_size)\
					break;\
			}\
		} else {\
			while ((ret = kseq_read(p->ks)) >= 0) {\
				int l = (int)(p->ks->seq.l) - (int)(p->opt->adaLen) - (int)(p->opt->adaLen);\
				if((l <= 0) || (l < asm_opt.rl_cut)) continue;\
				if((asm_opt.is_sc) && (asm_opt.sc_cut > 0) && (!flt_quals(p->ks->qual.s+p->opt->adaLen, l, 33, asm_opt.sc_cut))) continue;\
				if (p->n_seq >= 1<<28) {\
					fprintf(stderr, "ERROR: this implementation supports no more than %d reads\n", 1<<28);\
					exit(1);\
				}\
				if (p->rs_out) {\
					/**for 0-th count, just insert read length to R_INF, instead of read**/\
					if (p->flag & HAF_RS_WRITE_LEN) {\
						assert(p->n_seq == p->rs_out->total_reads);\
						ha_insert_read_len(p->rs_out, l, p->ks->name.l);\
						if (p->flag & HAF_RS_STASH) ha_rss_add_name(p->ks->name.s, p->ks->name.l);\
					} else if (p->flag & HAF_RS_WRITE_SEQ) {\
						/**the bases are compressed in parallel in step 1 (sf##_worker_for_rs)**/\
						assert(l == (int)p->rs_out->read_length[p->n_seq]);\
						memcpy(&p->rs_out->name[p->rs_out->name_index[p->n_seq]], p->ks->name.s, p->ks->name.l);\
						if(p->rs_out->rsc) {\
							ha_compress_qual(Get_QUAL(*p->rs_out, p->n_seq), p->ks->qual.s+p->opt->adaLen, l, sc_bn, 33);\
							/**print_fastq(NULL, p->ks->name.s, p->ks->seq.s, p->ks->qual.s, (1<<sc_bn), 33);**/\
							/**if(l <= 1000000) {\
								convert_qual(src_a, p->ks->qual.s+p->opt->adaLen, l, (1<<sc_bn), 0, 33);\
								retrive_bqual(NULL, des_a, p->n_seq, -1, -1, 0, sc_bn);\
								if(memcmp(src_a, des_a, l)!=0) fprintf(stderr, "ERROR: incorrect qual values\n");\
								else fprintf(stderr, "Correct: correct qual values\n");\
							}**/\
						}\
					}\
				}\
				if (s->n_seq == s->m_seq) {\
					s->m_seq = s->m_seq < 16? 16 : s->m_seq + (s->m_seq>>1);\
					REALLOC(s->len, s->m_seq);\
					REALLOC(s->seq, s->m_seq);\
				}\
				MALLOC(s->seq[s->n_seq], l);\
				memcpy(s->seq[s->n_seq], p->ks->seq.s+p->opt->adaLen, l);\
				s->len[s->n_seq++] = l;\
				++p->n_seq;\
				s->sum_len += l;\
				s->nk += l >= p->opt->k? l - p->opt->k + 1 : 0;\
				/**p->opt->chunk_size is the block max size**/\
				if (s->sum_len >= p->opt->chunk_size)\
					break;\
			}\
		}\
		if (s->sum_len == 0) free(s);\
		else return s;\
	} else if (step == 1) { /** step 2: extract k-mers**/\
		/**s is the block of reads**/\
		sf##_st_data_t *s = (sf##_st_data_t*)in;\
		int i, n_pre = 1<<p->opt->pre, m;\
		if (p->ks && p->rs_out && (p->flag & HAF_RS_WRITE_LEN) && (p->flag & HAF_RS_STASH)) {\
			CALLOC(s->rss_seq, s->n_seq); CALLOC(s->rss_N, s->n_seq);\
			kt_for(p->opt->n_thread, sf##_worker_for_rss, s, s->n_seq);\
		}\
		if (p->ks && p->rs_out && (p->flag & HAF_RS_WRITE_SEQ) && !(p->flag & HAF_RS_WRITE_LEN))\
			kt_for(p->opt->n_thread, sf##_worker_for_rs, s, s->n_seq);\
		/**allocate the k-mer buffer**/\
		s->buf = sf##_bin_take(p, n_pre);\
		m = (int)(s->nk * 1.2 / n_pre) + 1;\
		/**pre-allocate memory for each of 4096 buffer; for minimizers, sf##_mzbin_dist() allocates\
		 * the exact size (m is for all k-mers: large mallocs and frees for each chunk)**/\
		for (i = 0; i < n_pre && p->opt->w == 1; ++i) {\
			if (s->buf[i].m > (int64_t)m*2) {\
				free(s->buf[i].a); free(s->buf[i].b);\
				s->buf[i].a = NULL; s->buf[i].b = NULL; s->buf[i].m = 0;\
			}\
			if (s->buf[i].m) continue;\
			s->buf[i].m = m;\
			/**for 0-th counting, p->pt = NULL**/\
			if (p->pt && !(p->flag&HAF_COUNT_REFINE)) MALLOC(s->buf[i].b, m);\
			else MALLOC(s->buf[i].a, m);\
		}\
		if (p->opt->w == 1) { /** enumerate all k-mers, in parallel (ha_kblk_*)**/\
			int i; ha_kblk_t kb; int64_t *bin_n; uint64_t **dst;\
			ha_kblk_extract(&kb, p->opt->n_thread, s->n_seq, s->len, s->seq, p->opt->k, p->opt->pre, p->opt->is_HPC);\
			MALLOC(bin_n, n_pre); MALLOC(dst, n_pre);\
			ha_kblk_bin_count(&kb, bin_n);\
			for (i = 0; i < n_pre; ++i) {\
				if (bin_n[i] > s->buf[i].m) {\
					s->buf[i].m = bin_n[i];\
					REALLOC(s->buf[i].a, s->buf[i].m);\
				}\
				s->buf[i].n = bin_n[i]; dst[i] = s->buf[i].a;\
			}\
			ha_kblk_scatter(&kb, p->opt->n_thread, dst);\
			HA_BINREUSE_VERIFY(sf, s, 0);\
			free(bin_n); free(dst);\
			for (i = 0; i < s->n_seq; ++i)\
				if (!p->is_store) free(s->seq[i]);\
		} else { /** minimizers only**/\
			/**s->n_seq is how many reads at this buffer**/\
			/**s->mz && s->mz_buf are lists of minimzer vectors**/\
			CALLOC(s->mz, s->n_seq);\
			if (ha_mzc.mode == 2) {/**replay the minimizers of the first pass (ha_mzc), into one array s->mzf**/\
				assert(ha_mzc.pos + s->n_seq <= ha_mzc.n);\
				ha_mzc_get(&ha_mzc, ha_mzc.pos, s->n_seq, (ha_mz1_v*)s->mz);\
				s->mzf = ha_mzc.f;\
				s->kdec = ha_mzc.a + ha_mzc.pos;\
				ha_mzc.pos += s->n_seq;\
				ha_mzc.replayed = 1;\
			} else {\
				CALLOC(s->mz_buf, p->opt->n_thread), CALLOC(s->mt, p->opt->n_thread);\
				/**calculate minimzers for each read, each read corresponds to one thread**/\
				kt_for(p->opt->n_thread, sf##_worker_for_mz, s, s->n_seq);\
				for (i = 0; i < p->opt->n_thread; ++i) free(s->mt[i].a), free(s->mz_buf[i].a);\
				free(s->mt), free(s->mz_buf);\
				if (ha_mzp.a) ha_mzc_free_before(&ha_mzp, s->n_seq0 + s->n_seq);/**the mappings of the reads before**/\
			}\
			if (ha_mzc.mode == 1) {/**record the minimizers for the second pass (ha_mzc), in sf##_mzbin_dist()**/\
				s->kenc = ha_mzc_grow(&ha_mzc, s->n_seq);\
				MALLOC(s->koff, s->n_seq);\
			}\
			/**insert minimizers: whole minimizers to the position table, else the hash keys only**/\
			sf##_mzbin_dist(s, p->opt->n_thread, p->opt->pre, p->pt && !(p->flag&HAF_COUNT_REFINE));\
			HA_BINREUSE_VERIFY(sf, s, p->pt && !(p->flag&HAF_COUNT_REFINE));\
			if (s->kenc) {\
				ha_mzc_commit(&ha_mzc, s->n_seq);\
				free(s->koff); s->koff = NULL; s->kenc = NULL; s->kb = NULL;\
			}\
			s->kdec = NULL;\
			for (i = 0; i < s->n_seq; ++i) {\
				p->n_mz += s->mz[i].n;\
				if (s->mzf == NULL) free(s->mz[i].a);\
				if (!p->is_store) free(s->seq[i]);\
			}\
			free(s->mz); s->mzf = NULL;\
		}\
		/**just clean seq**/\
		free(s->seq); free(s->len); free(s->hit);\
		s->seq = 0, s->len = 0, s->hit = 0;\
		return s;\
	} else if (step == 2) { /** step 3: insert k-mers to hash table**/\
		sf##_st_data_t *s = (sf##_st_data_t*)in;\
		int i, n = 1<<p->opt->pre;uint64_t n_ins = 0;\
		/**for 0-th counting, p->pt = NULL**/\
		kt_for(p->opt->n_thread, sf##_worker_for_insert, s, n);\
		if (s->rss_seq) {\
			ha_rss_add_reads(s->n_seq0, s->n_seq, s->rss_seq, s->rss_N);\
			free(s->rss_seq); free(s->rss_N);\
		}\
		/**n_ins is number of distinct k-mers**/\
		for (i = 0; i < n; ++i) {\
			n_ins += s->buf[i].n_ins;\
		}\
		if (p->ct) p->ct->tot += n_ins, p->ct->bs += s->sum_len;\
		if (p->pt) p->pt->tot_pos += n_ins;\
		pthread_mutex_lock(&p->bin_mutex);\
		assert(p->n_bin_spare < 3);\
		p->bin_spare[p->n_bin_spare++] = s->buf;\
		pthread_mutex_unlock(&p->bin_mutex);\
		free(s);\
	}\
	return 0;\
}

MZ_TEST_INIT(mz1, ha_mz1_t, ha_mz1_v, ha_idxpos_t, a)
MZ_TEST_INIT(mz2, ha_mzl_t, ha_mzl_v, ha_idxposl_t, al)


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

#ifdef HA_PL_TIME
///the time of each step of the ha_count() pipeline, to stderr
static void *ha_plt_worker(void *data, int step, void *in)
{
	double t0 = yak_realtime(); void *r;
	r = mz1_worker_count(data, step, in);
	fprintf(stderr, "[PLT] %d %.4f %.4f\n", step, t0, yak_realtime());
	return r;
}
#endif

static ha_ct_t *yak_count(const yak_copt_t *opt, const char *fn, int flag, ha_pt_t *p0, ha_ct_t *c0, const void *flt_tab, All_reads *rs, ma_utg_v *us, int64_t *n_seq)
{
	///for 0-th counting, flag = HAF_COUNT_ALL|HAF_RS_WRITE_LEN|HAF_CREATE_NEW
	int read_rs = (rs && (flag & HAF_RS_READ));
	int ug_rs = (us && (flag & HAF_UG_READ));
	pl_data_t pl;
	gzFile fp = 0;
	memset(&pl, 0, sizeof(pl_data_t));
	pl.n_seq = *n_seq; 
	if(ug_rs) {
		pl.us_in = us;
	} else if (read_rs) {
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
	if(pl.ct) pl.ct->bs = 0;
	pthread_mutex_init(&pl.bin_mutex, NULL);
	if(ug_rs) kt_pipeline(3, mz2_worker_count, &pl, 3);
#ifdef HA_PL_TIME
	else kt_pipeline(3, ha_plt_worker, &pl, 3);
#else
	else kt_pipeline(3, mz1_worker_count, &pl, 3);
#endif
	if (ug_rs) mz2_bin_free(&pl);
	else mz1_bin_free(&pl);
	pthread_mutex_destroy(&pl.bin_mutex);
	if (read_rs) {
		destory_UC_Read(&pl.ucr);
	} else if(!read_rs && !ug_rs) {
		kseq_destroy(pl.ks);
		gzclose(fp);
	}
	*n_seq = pl.n_seq;
	if (pl.opt->w > 1) fprintf(stderr, "[M::%s] collected %ld minimizers\n", __func__, (long)pl.n_mz);
	return pl.ct;
}

ha_ct_t *ha_count(const hifiasm_opt_t *asm_o, int flag, int HPC, int k, int w, ha_pt_t *p0, const void *flt_tab, All_reads *rs, ma_utg_v *us, int keep_adapter, int *low_freq, int unique_only)
{
#ifdef HA_BINREUSE_CHECK
	static int reuse_registered = 0;
	if (!reuse_registered) atexit(ha_binreuse_report), reuse_registered = 1;
#endif
	int i;
	int64_t n_seq = 0;
	uint64_t n_bs = 0;
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
	opt.k = k;
	opt.is_HPC = HPC;
	///for ft-counting, shoud be 1
	opt.w = flag & HAF_COUNT_ALL? 1 : w;
	///for ft-counting, shoud be 37
	///for ha_pt_gen, shoud be 0
	opt.bf_shift = flag & HAF_COUNT_EXACT? 0 : asm_o->bf_shift;
	opt.n_thread = asm_o->thread_num;
	opt.adaLen = (keep_adapter? asm_o->adapterLen : 0);
	opt.min_rcnt = (low_freq?*low_freq:-1);
	opt.uq = (unique_only?1:0);
	
	n_bs = 0;
	///pending for integration
	/**
	if(rs && asm_o->ar && asm_o->ul_mod) {
		for (i = 0; i < asm_o->ar->n; ++i){
			h = yak_count(&opt, asm_o->ar->a[i], flag|HAF_CREATE_NEW, p0, h, flt_tab, rs, us, &n_seq);
			if(h) n_bs += h->bs;
		}
	}**/
	///asm_opt->num_reads is the number of fastq files
	for (i = 0; i < (us?1:asm_o->num_reads); ++i){
		h = yak_count(&opt, asm_o->read_file_names[i], flag|HAF_CREATE_NEW, p0, h, flt_tab, rs, us, &n_seq);
		if(h) n_bs += h->bs;
	}
	if(h) h->bs = n_bs;	
	if (h && opt.bf_shift > 0)
		ha_ct_destroy_bf(h);
	return h;
}

/***************************
 * High count filter table *
 ***************************/

// Warning: the max count is 32767
KHASHL_MAP_INIT(static klib_unused, yak_ft_t, yak_ft, uint64_t, int16_t, kh_hash_dummy, kh_eq_generic)

static yak_ft_t *gen_hh(const ha_ct_t *h, int max_cnt)
{
	int i;
	yak_ft_t *hh;
	if (max_cnt > YAK_MAX_COUNT - 1) max_cnt = YAK_MAX_COUNT - 1;
	if (max_cnt > INT16_MAX - 1) max_cnt = INT16_MAX - 1;
	hh = yak_ft_init();
	yak_ft_resize(hh, h->tot * 2);
	for (i = 0; i < 1<<h->pre; ++i) {
		yak_ct_t *ht = h->h[i].h;
		khint_t k, l;
		for (k = 0; k < kh_end(ht); ++k) {
			if (kh_exist(ht, k)) {
				uint64_t y = kh_key(ht, k) >> h->pre << YAK_COUNTER_BITS | i;
				int absent;
				l = yak_ft_put(hh, y, &absent);
				if (absent) {
					int cnt = kh_key(ht, k) & YAK_MAX_COUNT;
					kh_val(hh, l) = cnt > max_cnt? INT16_MAX : cnt;
				}
			}
		}
	}
	return hh;
}

int32_t ha_ft_cnt(const void *hh, uint64_t y)
{
	yak_ft_t *h = (yak_ft_t*)hh;
	khint_t k;
	k = yak_ft_get(h, y);
	return k == kh_end(h)? 0 : kh_val(h, k) == INT16_MAX? INT32_MAX : kh_val(h, k);
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
void *ha_ft_ul_gen(const hifiasm_opt_t *asm_opt, ma_utg_v *us, int k, int w, int cutoff)
{
	yak_ft_t *flt_tab;
    ha_ct_t *h;
	///HAF_COUNT_EXACT ---> no bf; HAF_COUNT_ALL ---> no minimizer
	h = ha_count(asm_opt, HAF_COUNT_ALL|HAF_UG_READ, !(asm_opt->flag&HA_F_NO_HPC), k, w, NULL, NULL, NULL, us, 0, NULL, 0);
	
	// cutoff = (int)(asm_opt->hom_cov * asm_opt->high_factor);
    if (cutoff > YAK_MAX_COUNT - 1) cutoff = YAK_MAX_COUNT - 1;
	// fprintf(stderr, "[M::%s::] cutoff->%d\n\n", __func__, cutoff);
	ha_ct_shrink(h, cutoff, YAK_MAX_COUNT, asm_opt->thread_num);
	flt_tab = gen_hh(h, asm_opt->max_kmer_cnt);
	ha_ct_destroy(h);
	fprintf(stderr, "[M::%s::%.3f*%.2f@%.3fGB] ==> filtered out %ld k-mers occurring %d or more times\n", __func__,
            yak_realtime(), yak_cpu_usage(), yak_peakrss_in_gb(), (long)kh_size(flt_tab), cutoff);
	HA_MEMSTAT_AT("filter count end");
	return (void*)flt_tab;
}

void *ha_ft_ug_gen(const hifiasm_opt_t *asm_opt, ma_utg_v *us, int is_HPC, int k, int w, int min_freq, int max_freq)
{
	yak_ft_t *flt_tab;
	ha_ct_t *h;
	///HAF_COUNT_EXACT ---> no bf; HAF_COUNT_ALL ---> no minimizer
	h = ha_count(asm_opt, HAF_COUNT_ALL|HAF_UG_READ|HAF_COUNT_EXACT, is_HPC, k, w, NULL, NULL, NULL, us, 0, NULL, 0);
	ha_ct_shrink(h, min_freq, max_freq>YAK_MAX_COUNT-1?YAK_MAX_COUNT-1:max_freq, asm_opt->thread_num);
	flt_tab = gen_hh(h, YAK_MAX_COUNT);
	ha_ct_destroy(h);
	return (void*)flt_tab;
}

void *ha_ft_gen(const hifiasm_opt_t *asm_opt, All_reads *rs, int *hom_cov, int is_hp_mode, int read_from_store)
{
	yak_ft_t *flt_tab;
	int64_t cnt[YAK_N_COUNTS];
	int peak_hom, peak_het, cutoff = YAK_MAX_COUNT - 1, ex_flag = 0;
	if(is_hp_mode) ex_flag = HAF_RS_READ|HAF_SKIP_READ;
	ha_rss_free();
	if(!is_hp_mode && !read_from_store && !asm_opt->is_sc) ex_flag |= HAF_RS_STASH;///keep the reads for the first ha_pt_gen()
	ha_ct_t *h; 
	h = ha_count(asm_opt, HAF_COUNT_ALL|ex_flag|((read_from_store)?(HAF_RS_READ):(HAF_RS_WRITE_LEN)), !(asm_opt->flag&HA_F_NO_HPC), asm_opt->k_mer_length, asm_opt->mz_win, NULL, NULL, rs, NULL, 1, NULL, 0);
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
		peak_hom = ha_analyze_count(YAK_N_COUNTS, asm_opt->min_hist_kmer_cnt, asm_opt->hg_size>0?(h->bs/asm_opt->hg_size):(-1), cnt, &peak_het);
		if (hom_cov) *hom_cov = peak_hom;
		if (peak_hom > 0) fprintf(stderr, "[M::%s] peak_hom: %d; peak_het: %d\n", __func__, peak_hom, peak_het);
		///in default, asm_opt->high_factor = 5.0
		cutoff = (int)(peak_hom * asm_opt->high_factor);
		if (cutoff > YAK_MAX_COUNT - 1) cutoff = YAK_MAX_COUNT - 1;
	}
	ha_ct_shrink(h, cutoff, YAK_MAX_COUNT, asm_opt->thread_num);
	flt_tab = gen_hh(h, asm_opt->max_kmer_cnt);
	ha_ct_destroy(h);
	fprintf(stderr, "[M::%s::%.3f*%.2f@%.3fGB] ==> filtered out %ld k-mers occurring %d or more times\n", __func__,
			yak_realtime(), yak_cpu_usage(), yak_peakrss_in_gb(), (long)kh_size(flt_tab), cutoff);
	HA_MEMSTAT_AT("filter count end");
	return (void*)flt_tab;
}

ha_pt_t *ha_pt_ul_gen(const hifiasm_opt_t *asm_opt, const void *flt_tab, ma_utg_v *us, int k, int w, int cutoff)
{
	ha_ct_t *ct;
	ha_pt_t *pt;
	///HAF_COUNT_EXACT: no bf
	ct = ha_count(asm_opt, HAF_COUNT_EXACT|HAF_UG_READ, !(asm_opt->flag&HA_F_NO_HPC), k, w, NULL, flt_tab, NULL, us, 0, NULL, 0);
	fprintf(stderr, "[M::%s::%.3f*%.2f] ==> counted %ld distinct minimizer k-mers\n", __func__,
			yak_realtime(), yak_cpu_usage(), (long)ct->tot);
	///minimizer with YAK_MAX_COUNT occ may apper > YAK_MAX_COUNT times, so it may lead to overflow at ha_pt_gen
	if (flt_tab == 0) {
        if (cutoff > YAK_MAX_COUNT - 1) cutoff = YAK_MAX_COUNT - 1;
        ha_ct_shrink(ct, /**2**/1, cutoff, asm_opt->thread_num);
    } else {
        ///Note: here is just to remove minimizer appearing YAK_MAX_COUNT times
        ///minimizer with YAK_MAX_COUNT occ may apper > YAK_MAX_COUNT times, so it may lead to overflow at ha_pt_gen
        ha_ct_shrink(ct, /**2**/1, YAK_MAX_COUNT - 1, asm_opt->thread_num);
    }
	
	pt = ha_pt_gen(ct, asm_opt->thread_num, 1);
	ha_count(asm_opt, HAF_COUNT_EXACT|HAF_UG_READ, !(asm_opt->flag&HA_F_NO_HPC), k, w, pt, flt_tab, NULL, us, 0, NULL, 0);
	//ha_pt_sort(pt, asm_opt->thread_num);
	fprintf(stderr, "[M::%s::%.3f*%.2f] ==> indexed %ld positions\n", __func__,
			yak_realtime(), yak_cpu_usage(), (long)pt->tot_pos);
	return pt;
}

ha_pt_t *ha_pt_ug_gen(const hifiasm_opt_t *asm_opt, const void *flt_tab, ma_utg_v *us, int is_HPC, int k, int w, int min_freq)
{
	ha_ct_t *ct;
	ha_pt_t *pt;
	///HAF_COUNT_EXACT: no bf
	ct = ha_count(asm_opt, HAF_COUNT_EXACT|HAF_UG_READ, is_HPC, k, w, NULL, flt_tab, NULL, us, 0, NULL, 1);
	fprintf(stderr, "[M::%s::%.3f*%.2f] ==> counted %ld distinct minimizer k-mers\n", __func__,
			yak_realtime(), yak_cpu_usage(), (long)ct->tot);
	///minimizer with YAK_MAX_COUNT occ may apper > YAK_MAX_COUNT times, so it may lead to overflow at ha_pt_gen
	ha_ct_shrink(ct, min_freq, YAK_MAX_COUNT - 1, asm_opt->thread_num);

	pt = ha_pt_gen(ct, asm_opt->thread_num, 1);
	ha_count(asm_opt, HAF_COUNT_EXACT|HAF_UG_READ, is_HPC, k, w, pt, flt_tab, NULL, us, 0, NULL, 1);
	//ha_pt_sort(pt, asm_opt->thread_num);
	fprintf(stderr, "[M::%s::%.3f*%.2f] ==> indexed %ld positions\n", __func__,
			yak_realtime(), yak_cpu_usage(), (long)pt->tot_pos);
	return pt;
}

ha_pt_t *ha_pt_gen_dp(const hifiasm_opt_t *asm_opt, ha_ct_t *ct, int flag, int n_thread, const void *flt_tab, All_reads *rs, int peak_hom, int peak_het)
{
	int low_freq = mz_low_b(peak_hom, peak_het);
	ha_pt_t *pt = ha_pt_gen_count(ct, n_thread); ///key = cnt, val = 0
	ha_count(asm_opt, HAF_COUNT_EXACT|HAF_COUNT_REFINE|flag, !(asm_opt->flag&HA_F_NO_HPC), asm_opt->k_mer_length, asm_opt->mz_win, pt, flt_tab, rs, NULL, 1, &low_freq, 0);
	uint64_t occ = ha_pt_shrink(pt, n_thread);
	if(flag&HAF_RS_WRITE_LEN) flag -= HAF_RS_WRITE_LEN;
	if(flag&HAF_RS_WRITE_SEQ) flag -= HAF_RS_WRITE_SEQ;
	flag |= HAF_RS_READ; pt->tot_pos = 0;
	ha_count(asm_opt, HAF_COUNT_EXACT|flag, !(asm_opt->flag&HA_F_NO_HPC), asm_opt->k_mer_length, asm_opt->mz_win, pt, flt_tab, rs, NULL, 1, NULL, 0);
	// fprintf(stderr, "[M::%s::] counted %lu distinct minimizer k-mers\n", __func__, pt->tot);
	// fprintf(stderr, "[M::%s::] collected %lu minimizers\n\n\n", __func__, pt->tot_pos);
	assert(occ == pt->tot_pos);
	return pt;
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
		if (!is_hp_mode && ha_rss_adopt(rs)) extra_flag1 = HAF_RS_READ;///the reads of the filter count pass
	}
	ha_rss_free();///reads that are not adopted
	HA_MEMSTAT_AT("index start");
	if(is_hp_mode) extra_flag1 |= HAF_SKIP_READ, extra_flag2 |= HAF_SKIP_READ;

	///keep the minimizers of the count for the index (ha_mzc), in the normal mode; the cache of the
	///queries of the last round gives the minimizers of the reads that did not change (ha_mzp)
	ha_mzp_free();
	if (ha_mzc.query && ha_rd_dirty && ha_rd_n == rs->total_reads && ha_mzc.n == rs->total_reads && read_from_store && !is_hp_mode && !(asm_opt->flag & HA_F_FAST)) {
		ha_mzp = ha_mzc; memset(&ha_mzc, 0, sizeof(ha_mzc));
		ha_mzp_dirty = ha_rd_dirty; ha_rd_dirty = NULL;
		ha_mzp_drop_dirty();
	}
	free(ha_rd_dirty); ha_rd_dirty = NULL; ha_rd_n = 0;
	ha_mzc_free();
	ha_mem_trim();///the free blocks of the filter pass or of the last round
	HA_MEMSTAT_AT("index trim");
	if(!(asm_opt->flag & HA_F_FAST)) ha_mzc.mode = 1;
	ct = ha_count(asm_opt, HAF_COUNT_EXACT|extra_flag1,  !(asm_opt->flag&HA_F_NO_HPC), asm_opt->k_mer_length, asm_opt->mz_win, NULL, flt_tab, rs, NULL, 1, NULL, 0);
	ha_mzc.mode = 0;
	fprintf(stderr, "[M::%s::%.3f*%.2f] ==> counted %ld distinct minimizer k-mers\n", __func__,
			yak_realtime(), yak_cpu_usage(), (long)ct->tot);
	HA_MEMSTAT_AT("index count end");
	ha_ct_hist(ct, cnt, asm_opt->thread_num);
	fprintf(stderr, "[M::%s] count[%d] = %ld (for sanity check)\n", __func__, YAK_MAX_COUNT, (long)cnt[YAK_MAX_COUNT]);
	peak_hom = ha_analyze_count(YAK_N_COUNTS, asm_opt->min_hist_kmer_cnt, asm_opt->hg_size>0?(ct->bs/asm_opt->hg_size):(-1), cnt, &peak_het);
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
	if(!(asm_opt->flag & HA_F_FAST))
	{
		fprintf(stderr, "[M::%s::] counting in normal mode\n", __func__);
		pt = ha_pt_gen(ct, asm_opt->thread_num, 0);
		///replay only if the second pass reads the stored reads (HAF_RS_READ), as the first pass did
		if((extra_flag2 & HAF_RS_READ) && ha_mzc.n > 0) ha_mzc.mode = 2;
		ha_count(asm_opt, HAF_COUNT_EXACT|extra_flag2,  !(asm_opt->flag&HA_F_NO_HPC), asm_opt->k_mer_length, asm_opt->mz_win, pt, flt_tab, rs, NULL, 1, NULL, 0);
		if(ha_mzc.mode == 2) assert(ha_mzc.pos == ha_mzc.n);
		assert((uint64_t)tot_cnt == pt->tot_pos);
	}
	else
	{
		fprintf(stderr, "[M::%s::] counting in fast mode\n", __func__);
		pt = ha_pt_gen_dp(asm_opt, ct, HAF_COUNT_EXACT|extra_flag2, asm_opt->thread_num, flt_tab, rs, peak_hom, peak_het);
	}
	///after a full replay of all stored reads (no skipped read), keep the minimizers for the queries of
	///the next correction round (ha_mzc_query); else free them
	ha_mzp_free();
	if(ha_mzc.replayed && ha_mzc.n == rs->total_reads && !(extra_flag2 & HAF_SKIP_READ)) {
		ha_mzc.mode = 0; ha_mzc.query = 1;
		ha_rd_n = rs->total_reads; CALLOC(ha_rd_dirty, ha_rd_n);///from now on, the writers mark the changed reads
	} else {
		ha_mzc_free();
	}
	//ha_pt_sort(pt, asm_opt->thread_num);
	fprintf(stderr, "[M::%s::%.3f*%.2f] ==> indexed %ld positions, counted %ld distinct minimizer k-mers\n", __func__,
			yak_realtime(), yak_cpu_usage(), (long)pt->tot_pos, (long)pt->tot);
	HA_MEMSTAT_AT("index end");
	return pt;
}

int query_ct_index(void* ct_idx, uint64_t hash)
{
	ha_ct1_t *g = &(((ha_ct_t*)ct_idx)->h[hash & ((1ULL<<((ha_ct_t*)ct_idx)->pre) - 1)]);
	khint_t k;
	k = yak_ct_get(g->h, hash);
	if (k == kh_end(g->h)) return 0;
	return ((kh_key(g->h, k)&YAK_MAX_COUNT)==YAK_MAX_COUNT)?-1:(kh_key(g->h, k)&YAK_MAX_COUNT);
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
    char* gfa_name = (char*)malloc(strlen(file_name)+64);
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

	sprintf(gfa_name, "%s.pt_flt.paf.bin", file_name); 
	fclose(fp); fp = fopen(gfa_name, "w"); uint64_t k;
	if (!fp) {
		free(gfa_name);
        return 0;
    }
	fwrite(&(r->total_reads), sizeof(r->total_reads), 1, fp);
	for (k = 0; k < r->total_reads; k++) {
		fwrite(&(r->paf[k].is_fully_corrected), sizeof(r->paf[k].is_fully_corrected), 1, fp);
        fwrite(&(r->paf[k].is_abnormal), sizeof(r->paf[k].is_abnormal), 1, fp);
        fwrite(&(r->paf[k].length), sizeof(r->paf[k].length), 1, fp);
        fwrite(r->paf[k].buffer, sizeof((*(r->paf[k].buffer))), r->paf[k].length, fp);
	}

	fprintf(stderr, "[M::%s] Index has been written.\n", __func__);
    free(gfa_name);
	fclose(fp);
	return 1;
}

int load_pt_index(void **r_flt_tab, ha_pt_t **r_ha_idx, All_reads *r, hifiasm_opt_t* opt, char* file_name)
{
	char* gfa_name = (char*)malloc(strlen(file_name)+64);
    sprintf(gfa_name, "%s.pt_flt", file_name);
    FILE* fp = fopen(gfa_name, "r");
	if (!fp) {
		free(gfa_name);
        return 0;
    }

	ha_pt_t *ha_idx = NULL;
	char mode = 0;
	int f_flag = 0, absent, i;
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


	// fclose(fp);
	
	if(!load_All_reads(r, gfa_name))
	{
		free(gfa_name);
		return 0;
	}

	memset(r->trio_flag, AMBIGU, r->total_reads*sizeof(uint8_t));

	sprintf(gfa_name, "%s.pt_flt.paf.bin", file_name);
	fclose(fp); fp = fopen(gfa_name, "r"); uint64_t k;
    if (!fp) {
        free(gfa_name);
        return 0;
    }
	f_flag += fread(&(r->total_reads), sizeof(r->total_reads), 1, fp);
	r->paf = (ma_hit_t_alloc*)malloc(sizeof(ma_hit_t_alloc)*r->total_reads);
    r->reverse_paf = (ma_hit_t_alloc*)malloc(sizeof(ma_hit_t_alloc)*r->total_reads);
	for (k = 0; k < r->total_reads; k++) {
        // init_ma_hit_t_alloc(&(r->paf[k]));
        init_ma_hit_t_alloc(&(r->reverse_paf[k]));

		f_flag += fread(&(r->paf[k].is_fully_corrected), sizeof(r->paf[k].is_fully_corrected), 1, fp);
        f_flag += fread(&(r->paf[k].is_abnormal), sizeof(r->paf[k].is_abnormal), 1, fp);
        f_flag += fread(&(r->paf[k].length), sizeof(r->paf[k].length), 1, fp);
        r->paf[k].size = r->paf[k].length;

        r->paf[k].buffer = NULL;
        if(r->paf[k].length == 0) continue;
        
        r->paf[k].buffer = (ma_hit_t*)malloc(sizeof(ma_hit_t)*r->paf[k].length);
        fread(r->paf[k].buffer, sizeof((*(r->paf[k].buffer))), r->paf[k].length, fp);
    }
	fclose(fp);

	fprintf(stderr, "[M::%s] Index has been loaded.\n", __func__);

	free(gfa_name);
	return 1;
}



int uidx_write(void *flt_tab, ha_pt_t *ha_idx, char* file_name, ma_ug_t *ug)
{
	char* gfa_name = (char*)malloc(strlen(file_name)+25);
    sprintf(gfa_name, "%s.uidx.bin", file_name);
    FILE* fp = fopen(gfa_name, "w");
    if (!fp) {
        free(gfa_name);
        return 0;
    }

	if(ug) write_dbug(ug, fp);

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
            fwrite(g->al, sizeof(ha_idxposl_t), g->n, fp);
        }
    }
	fprintf(stderr, "[M::%s] Index has been written.\n", __func__);
    free(gfa_name);
    fclose(fp);
    return 1;
}


int uidx_load(void **r_flt_tab, ha_pt_t **r_ha_idx, char* file_name, ma_ug_t *ug)
{
    char* gfa_name = (char*)malloc(strlen(file_name)+25);
    sprintf(gfa_name, "%s.uidx.bin", file_name);
    FILE* fp = fopen(gfa_name, "r");
    if (!fp) {
        free(gfa_name);
        return 0;
    }

	if(ug && (!test_dbug(ug, fp))) {
        fprintf(stderr, "[M::%s] Renew UL Index\n", __func__);
        free(gfa_name); fclose(fp);
        return 0;
    }
	// fprintf(stderr, "[M::%s]\t%s\tftell::%ld\n", __func__, file_name, ftell(fp));

    ha_pt_t *ha_idx = NULL;
    char mode = 0;
    int f_flag = 0, i;
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
            MALLOC(g->al, g->n);
            f_flag += fread(g->al, sizeof(ha_idxposl_t), g->n, fp);

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

    fprintf(stderr, "[M::%s] Index has been loaded.\n", __func__);

	fclose(fp);
    free(gfa_name);
    return 1;
}
