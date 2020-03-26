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

/***************************
 * Yak specific parameters *
 ***************************/

typedef struct {
	int32_t bf_shift, bf_n_hash;
	int32_t k, w, is_HPC;
	int32_t pre;
	int32_t n_thread;
	int64_t chunk_size;
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
	o->chunk_size = 10000000;
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
	uint64_t tot;
	ha_ct1_t *h;
} ha_ct_t;

static ha_ct_t *ha_ct_init(int k, int pre, int n_hash, int n_shift)
{
	ha_ct_t *h;
	int i;
	if (pre < YAK_COUNTER_BITS) return 0;
	CALLOC(h, 1);
	h->k = k, h->pre = pre;
	CALLOC(h->h, 1<<h->pre);
	for (i = 0; i < 1<<h->pre; ++i)
		h->h[i].h = yak_ct_init();
	if (n_hash > 0 && n_shift > h->pre) {
		h->n_hash = n_hash, h->n_shift = n_shift;
		for (i = 0; i < 1<<h->pre; ++i)
			h->h[i].b = yak_bf_init(h->n_shift - h->pre, h->n_hash);
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
	g = &h->h[a[0]&mask];
	for (j = 0; j < n; ++j) {
		int ins = 1, absent;
		uint64_t x = a[j] >> h->pre;
		khint_t k;
		if ((a[j]&mask) != (a[0]&mask)) continue;
		if (create_new) {
			if (g->b)
				ins = (yak_bf_insert(g->b, x) == h->n_hash);
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

static void ha_ct_hist(const ha_ct_t *h, int64_t cnt[YAK_N_COUNTS], int n_thread)
{
	hist_aux_t a;
	int i, j;
	a.h = h;
	memset(cnt, 0, YAK_N_COUNTS * sizeof(uint64_t));
	CALLOC(a.cnt, n_thread);
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
	uint64_t *a;
} ha_pt1_t;

typedef struct {
	int k, pre;
	uint64_t tot, tot_pos;
	ha_pt1_t *h;
} ha_pt_t;

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
			l = yak_pt_put(b->h, kh_key(g, k) >> YAK_COUNTER_BITS << YAK_COUNTER_BITS, &absent);
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

static int ha_pt_insert_list(ha_pt_t *h, int n, const ha_mz1_t *a)
{
	int j, mask = (1<<h->pre) - 1, n_ins = 0;
	ha_pt1_t *g;
	if (n == 0) return 0;
	g = &h->h[a[0].x&mask];
	for (j = 0; j < n; ++j) {
		uint64_t x = a[j].x >> h->pre;
		khint_t k;
		if ((a[j].x&mask) != (a[0].x&mask)) continue;
		k = yak_pt_get(g->h, x<<YAK_COUNTER_BITS);
		if (k == kh_end(g->h)) continue;
		assert((kh_key(g->h, k)&YAK_MAX_COUNT) < YAK_MAX_COUNT);
		g->a[kh_val(g->h, k) + (kh_key(g->h, k)&YAK_MAX_COUNT)]
			= (uint64_t)a[j].rid<<36 | (uint64_t)a[j].rev<<35 | (uint64_t)a[j].pos<<8 | (uint64_t)a[j].span;
		++kh_key(g->h, k);
		++n_ins;
	}
	return n_ins;
}

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

static void ha_pt_sort(ha_pt_t *h, int n_thread)
{
	kt_for(n_thread, worker_pt_sort, h, 1<<h->pre);
}

static void ha_pt_destroy(ha_pt_t *h)
{
	int i;
	if (h == 0) return;
	for (i = 0; i < 1<<h->pre; ++i) {
		yak_pt_destroy(h->h[i].h);
		free(h->h[i].a);
	}
	free(h->h); free(h);
}

void ha_idx_destroy(void *h)
{
	ha_pt_destroy((ha_pt_t*)h);
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

static inline void ct_insert_buf(ch_buf_t *buf, int p, uint64_t y) // insert a k-mer $y to a linear buffer
{
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
	int pre = y->x & ((1<<p) - 1);
	ch_buf_t *b = &buf[pre];
	if (b->n == b->m) {
		b->m = b->m < 8? 8 : b->m + (b->m>>1);
		REALLOC(b->b, b->m);
	}
	b->b[b->n++] = *y;
}

static void count_seq_buf(ch_buf_t *buf, int k, int p, int len, const char *seq) // insert k-mers in $seq to linear buffer $buf
{
	int i, l;
	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
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

typedef struct { // global data structure for kt_pipeline()
	const yak_copt_t *opt;
	const void *flt_tab;
	int flag, create_new, is_store;
	uint64_t n_seq;
	kseq_t *ks;
	UC_Read ucr;
	ha_ct_t *ct;
	ha_pt_t *pt;
	const All_reads *rs_in;
	All_reads *rs_out;
} pl_data_t;

typedef struct { // data structure for each step in kt_pipeline()
	pl_data_t *p;
	uint64_t n_seq0;
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
	else
		b->n_ins += ha_ct_insert_list(s->p->ct, s->p->create_new, b->n, b->a);
}

static void worker_for_mz(void *data, long i, int tid)
{
	st_data_t *s = (st_data_t*)data;
	ha_mz1_v *b = &s->mz_buf[tid];
	s->mz_buf[tid].n = 0;
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
				int l = p->ks->seq.l;
				if (p->n_seq >= 1<<28) {
					fprintf(stderr, "ERROR: this implementation supports no more than %d reads\n", 1<<28);
					exit(1);
				}
				if (p->rs_out) {
					if (p->flag & HAF_RS_WRITE_LEN) {
						assert(p->n_seq == p->rs_out->total_reads);
						ha_insert_read_len(p->rs_out, l, p->ks->name.l);
					} else if (p->flag & HAF_RS_WRITE_SEQ) {
						int i, n_N;
						assert(l == (int)p->rs_out->read_length[p->n_seq]);
						for (i = n_N = 0; i < l; ++i) // count number of ambiguous bases
							if (seq_nt4_table[(uint8_t)p->ks->seq.s[i]] >= 4)
								++n_N;
						ha_compress_base(Get_READ(*p->rs_out, p->n_seq), p->ks->seq.s, l, &p->rs_out->N_site[p->n_seq], n_N);
					}
				}
				if (s->n_seq == s->m_seq) {
					s->m_seq = s->m_seq < 16? 16 : s->m_seq + (s->m_seq>>1);
					REALLOC(s->len, s->m_seq);
					REALLOC(s->seq, s->m_seq);
				}
				MALLOC(s->seq[s->n_seq], l);
				memcpy(s->seq[s->n_seq], p->ks->seq.s, l);
				s->len[s->n_seq++] = l;
				++p->n_seq;
				s->sum_len += l;
				s->nk += l >= p->opt->k? l - p->opt->k + 1 : 0;
				if (s->sum_len >= p->opt->chunk_size)
					break;
			}
		}
		if (s->sum_len == 0) free(s);
		else return s;
	} else if (step == 1) { // step 2: extract k-mers
		st_data_t *s = (st_data_t*)in;
		int i, n_pre = 1<<p->opt->pre, m;
		// allocate the k-mer buffer
		CALLOC(s->buf, n_pre);
		m = (int)(s->nk * 1.2 / n_pre) + 1;
		for (i = 0; i < n_pre; ++i) {
			s->buf[i].m = m;
			if (p->pt) MALLOC(s->buf[i].b, m);
			else MALLOC(s->buf[i].a, m);
		}
		// fill the buffer
		if (p->opt->w == 1) { // enumerate all k-mers
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
			CALLOC(s->mz, s->n_seq);
			CALLOC(s->mz_buf, p->opt->n_thread);
			kt_for(p->opt->n_thread, worker_for_mz, s, s->n_seq);
			for (i = 0; i < p->opt->n_thread; ++i)
				free(s->mz_buf[i].a);
			free(s->mz_buf);
			// insert minimizers
			if (p->pt) {
				for (i = 0; i < s->n_seq; ++i)
					for (j = 0; j < s->mz[i].n; ++j)
						pt_insert_buf(s->buf, p->opt->pre, &s->mz[i].a[j]);
			} else {
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
		free(s->seq); free(s->len);
		s->seq = 0, s->len = 0;
		return s;
	} else if (step == 2) { // step 3: insert k-mers to hash table
		st_data_t *s = (st_data_t*)in;
		int i, n = 1<<p->opt->pre;
		uint64_t n_ins = 0;
		kt_for(p->opt->n_thread, worker_for_insert, s, n);
		for (i = 0; i < n; ++i) {
			n_ins += s->buf[i].n_ins;
			if (p->pt) free(s->buf[i].b);
			else free(s->buf[i].a);
		}
		if (p->ct) p->ct->tot += n_ins;
		if (p->pt) p->pt->tot_pos += n_ins;
		free(s->buf);
		fprintf(stderr, "[M::%s::%.3f*%.2f] processed %ld sequences; %ld %s in the hash table\n", __func__,
				yak_realtime(), yak_cputime() / yak_realtime(), (long)s->n_seq0 + s->n_seq,
				(long)(p->pt? p->pt->tot_pos : p->ct->tot), p->pt? "positions" : "distinct k-mers");
		free(s);
	}
	return 0;
}

static ha_ct_t *yak_count(const yak_copt_t *opt, const char *fn, int flag, ha_pt_t *p0, ha_ct_t *c0, const void *flt_tab, All_reads *rs)
{
	int read_rs = (rs && (flag & HAF_RS_READ));
	pl_data_t pl;
	gzFile fp = 0;
	memset(&pl, 0, sizeof(pl_data_t));
	if (read_rs) {
		pl.rs_in = rs;
		init_UC_Read(&pl.ucr);
	} else {
		if ((fp = gzopen(fn, "r")) == 0) return 0;
		pl.ks = kseq_init(fp);
	}
	if (rs && (flag & (HAF_RS_WRITE_LEN|HAF_RS_WRITE_SEQ)))
		pl.rs_out = rs;
	pl.flt_tab = flt_tab;
	pl.opt = opt;
	pl.flag = flag;
	if (p0) {
		pl.pt = p0, pl.create_new = 0;
		assert(p0->k == opt->k && p0->pre == opt->pre);
	} else if (c0) {
		pl.ct = c0, pl.create_new = 0;
		assert(c0->k == opt->k && c0->pre == opt->pre);
	} else {
		pl.create_new = 1;
		pl.ct = ha_ct_init(opt->k, opt->pre, opt->bf_n_hash, opt->bf_shift);
	}
	kt_pipeline(3, worker_count, &pl, 3);
	if (read_rs) {
		destory_UC_Read(&pl.ucr);
	} else {
		kseq_destroy(pl.ks);
		gzclose(fp);
	}
	return pl.ct;
}

ha_ct_t *ha_count(const hifiasm_opt_t *asm_opt, int flag, ha_pt_t *p0, const void *flt_tab, All_reads *rs)
{
	int i;
	yak_copt_t opt;
	ha_ct_t *h = 0;
	assert(!(flag & HAF_RS_WRITE_LEN) || !(flag & HAF_RS_WRITE_SEQ)); // not both
	if (rs) {
		if (flag & HAF_RS_WRITE_LEN)
			init_All_reads(rs);
		else if (flag & HAF_RS_WRITE_SEQ)
			malloc_All_reads(rs);
	}
	yak_copt_init(&opt);
	opt.k = asm_opt->k_mer_length;
	opt.is_HPC = !asm_opt->no_HPC;
	opt.w = flag & HAF_COUNT_ALL? 1 : asm_opt->mz_win;
	opt.bf_shift = flag & HAF_COUNT_EXACT? 0 : asm_opt->bf_shift;
	opt.n_thread = asm_opt->thread_num;
	for (i = 0; i < asm_opt->num_reads; ++i)
		h = yak_count(&opt, asm_opt->read_file_names[i], flag, p0, h, flt_tab, rs);
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
				uint64_t y = kh_key(ht, k) >> YAK_COUNTER_BITS << h->pre | i;
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
	yak_ft_destroy((yak_ft_t*)h);
}

/*************************
 * High-level interfaces *
 *************************/

void *ha_gen_flt_tab(const hifiasm_opt_t *asm_opt, All_reads *rs)
{
	yak_ft_t *flt_tab;
	int64_t cnt[YAK_N_COUNTS];
	int peak_hom, peak_het, cutoff;
	ha_ct_t *h;
	h = ha_count(asm_opt, HAF_COUNT_ALL|HAF_RS_WRITE_LEN, NULL, NULL, rs);
	ha_ct_hist(h, cnt, asm_opt->thread_num);
	peak_hom = yak_analyze_count(YAK_N_COUNTS, cnt, &peak_het);
	if (peak_hom > 0) fprintf(stderr, "[M::%s] peak_hom: %d; peak_het: %d\n", __func__, peak_hom, peak_het);
	cutoff = (int)(peak_hom * asm_opt->high_factor);
	if (cutoff > YAK_MAX_COUNT - 1) cutoff = YAK_MAX_COUNT - 1;
	ha_ct_shrink(h, cutoff, YAK_MAX_COUNT, asm_opt->thread_num);
	flt_tab = gen_hh(h);
	ha_ct_destroy(h);
	fprintf(stderr, "[M::%s::%.3f*%.2f] ==> filtered out %ld k-mers occurring %d or more times\n", __func__,
			yak_realtime(), yak_cputime() / yak_realtime(), (long)kh_size(flt_tab), cutoff);
	return (void*)flt_tab;
}

void *ha_gen_mzidx(const hifiasm_opt_t *asm_opt, const void *flt_tab, int read_from_store, All_reads *rs)
{
	int64_t cnt[YAK_N_COUNTS], tot_cnt;
	int peak_hom, peak_het, i, extra_flag = read_from_store? HAF_RS_READ : HAF_RS_WRITE_SEQ;
	ha_ct_t *ct;
	ha_pt_t *pt;
	ct = ha_count(asm_opt, HAF_COUNT_EXACT|extra_flag, NULL, flt_tab, rs);
	fprintf(stderr, "[M::%s::%.3f*%.2f] ==> counted %ld distinct minimizer k-mers\n", __func__,
			yak_realtime(), yak_cputime() / yak_realtime(), (long)ct->tot);
	ha_ct_hist(ct, cnt, asm_opt->thread_num);
	fprintf(stderr, "[M::%s] count[%d] = %ld (for sanity check)\n", __func__, YAK_MAX_COUNT, (long)cnt[YAK_MAX_COUNT]);
	peak_hom = yak_analyze_count(YAK_N_COUNTS, cnt, &peak_het);
	if (peak_hom > 0) fprintf(stderr, "[M::%s] peak_hom: %d; peak_het: %d\n", __func__, peak_hom, peak_het);
	ha_ct_shrink(ct, 2, YAK_MAX_COUNT - 1, asm_opt->thread_num);
	for (i = 2, tot_cnt = 0; i <= YAK_MAX_COUNT - 1; ++i) tot_cnt += cnt[i] * i;
	pt = ha_pt_gen(ct, asm_opt->thread_num);
	ha_count(asm_opt, HAF_COUNT_EXACT|HAF_RS_READ, pt, flt_tab, rs);
	assert((uint64_t)tot_cnt == pt->tot_pos);
	ha_pt_sort(pt, asm_opt->thread_num);
	fprintf(stderr, "[M::%s::%.3f*%.2f] ==> indexed %ld positions\n", __func__,
			yak_realtime(), yak_cputime() / yak_realtime(), (long)pt->tot_pos);
	return pt;
}
