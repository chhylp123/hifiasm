#include <stdint.h>
#include "CommandLines.h"
#include "yak.h"
#include "khashl.h"

#define YAK_COUNTER_BITS 12
#define YAK_N_COUNTS     (1<<YAK_COUNTER_BITS)
#define YAK_MAX_COUNT    ((1<<YAK_COUNTER_BITS)-1)

#define yak_ch_eq(a, b) ((a)>>YAK_COUNTER_BITS == (b)>>YAK_COUNTER_BITS) // lower 8 bits for counts; higher bits for k-mer
#define yak_ch_hash(a) ((a)>>YAK_COUNTER_BITS)
KHASHL_SET_INIT(static klib_unused, yak_ht_t, yak_ht, uint64_t, yak_ch_hash, yak_ch_eq)

KHASHL_SET_INIT(static klib_unused, yak_hh_t, yak_hh, uint64_t, kh_hash_dummy, kh_eq_generic)

typedef struct {
	int32_t bf_shift, bf_n_hash;
	int32_t k, w, is_HPC;
	int32_t pre;
	int32_t n_thread;
	int64_t chunk_size;
} yak_copt_t;

typedef struct {
	yak_ht_t *h;
	yak_bf_t *b;
} yak_ch1_t;

typedef struct {
	int k, pre, n_hash, n_shift;
	uint64_t tot;
	yak_ch1_t *h;
} yak_ch_t;

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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "kthread.h"

/*** hash table ***/

static yak_ch_t *yak_ch_init(int k, int pre, int n_hash, int n_shift)
{
	yak_ch_t *h;
	int i;
	if (pre < YAK_COUNTER_BITS) return 0;
	CALLOC(h, 1);
	h->k = k, h->pre = pre;
	CALLOC(h->h, 1<<h->pre);
	for (i = 0; i < 1<<h->pre; ++i)
		h->h[i].h = yak_ht_init();
	if (n_hash > 0 && n_shift > h->pre) {
		h->n_hash = n_hash, h->n_shift = n_shift;
		for (i = 0; i < 1<<h->pre; ++i)
			h->h[i].b = yak_bf_init(h->n_shift - h->pre, h->n_hash);
	}
	return h;
}

static void yak_ch_destroy_bf(yak_ch_t *h)
{
	int i;
	for (i = 0; i < 1<<h->pre; ++i) {
		if (h->h[i].b)
			yak_bf_destroy(h->h[i].b);
		h->h[i].b = 0;
	}
}

static void yak_ch_destroy(yak_ch_t *h)
{
	int i;
	if (h == 0) return;
	yak_ch_destroy_bf(h);
	for (i = 0; i < 1<<h->pre; ++i)
		yak_ht_destroy(h->h[i].h);
	free(h->h); free(h);
}

static int yak_ch_insert_list(yak_ch_t *h, int create_new, int n, const uint64_t *a)
{
	int j, mask = (1<<h->pre) - 1, n_ins = 0;
	yak_ch1_t *g;
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
				k = yak_ht_put(g->h, x<<YAK_COUNTER_BITS, &absent);
				if (absent) ++n_ins;
				if ((kh_key(g->h, k)&YAK_MAX_COUNT) < YAK_MAX_COUNT)
					++kh_key(g->h, k);
			}
		} else {
			k = yak_ht_get(g->h, x<<YAK_COUNTER_BITS);
			if (k != kh_end(g->h) && (kh_key(g->h, k)&YAK_MAX_COUNT) < YAK_MAX_COUNT)
				++kh_key(g->h, k);
		}
	}
	return n_ins;
}

static int yak_ch_get(const yak_ch_t *h, uint64_t x)
{
	int mask = (1<<h->pre) - 1;
	yak_ht_t *g = h->h[x&mask].h;
	khint_t k;
	k = yak_ht_get(g, x >> h->pre << YAK_COUNTER_BITS);
	return k == kh_end(g)? -1 : kh_key(g, k)&YAK_MAX_COUNT;
}

/*** generate histogram ***/

typedef struct {
	uint64_t c[YAK_N_COUNTS];
} buf_cnt_t;

typedef struct {
	const yak_ch_t *h;
	buf_cnt_t *cnt;
} hist_aux_t;

static void worker_hist(void *data, long i, int tid) // callback for kt_for()
{
	hist_aux_t *a = (hist_aux_t*)data;
	uint64_t *cnt = a->cnt[tid].c;
	yak_ht_t *g = a->h->h[i].h;
	khint_t k;
	for (k = 0; k < kh_end(g); ++k)
		if (kh_exist(g, k))
			++cnt[kh_key(g, k)&YAK_MAX_COUNT];
}

static void yak_ch_hist(const yak_ch_t *h, int64_t cnt[YAK_N_COUNTS], int n_thread)
{
	hist_aux_t a;
	int i, j;
	a.h = h;
	memset(cnt, 0, YAK_N_COUNTS * sizeof(uint64_t));
	CALLOC(a.cnt, n_thread);
	kt_for(n_thread, worker_hist, &a, 1<<h->pre);
	for (i = 0; i < YAK_N_COUNTS; ++i) cnt[i] = 0;
	for (j = 0; j < n_thread; ++j)
		for (i = 0; i < YAK_N_COUNTS; ++i)
			cnt[i] += a.cnt[j].c[i];
	free(a.cnt);
}

/*** shrink a hash table ***/

typedef struct {
	int min, max;
	yak_ch_t *h;
} shrink_aux_t;

static void worker_shrink(void *data, long i, int tid) // callback for kt_for()
{
	shrink_aux_t *a = (shrink_aux_t*)data;
	yak_ch_t *h = a->h;
	yak_ht_t *g = h->h[i].h, *f;
	khint_t k;
	f = yak_ht_init();
	yak_ht_resize(f, kh_size(g));
	for (k = 0; k < kh_end(g); ++k) {
		if (kh_exist(g, k)) {
			int absent, c = kh_key(g, k) & YAK_MAX_COUNT;
			if (c >= a->min && c <= a->max)
				yak_ht_put(f, kh_key(g, k), &absent);
		}
	}
	yak_ht_destroy(g);
	h->h[i].h = f;
}

static void yak_ch_shrink(yak_ch_t *h, int min, int max, int n_thread)
{
	int i;
	shrink_aux_t a;
	a.h = h, a.min = min, a.max = max;
	kt_for(n_thread, worker_shrink, &a, 1<<h->pre);
	for (i = 0, h->tot = 0; i < 1<<h->pre; ++i)
		h->tot += kh_size(h->h[i].h);
}

#include <zlib.h>
#include <string.h>
#include "kseq.h" // FASTA/Q parser
KSEQ_INIT(gzFile, gzread)

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

typedef struct {
	int n, m;
	uint64_t n_ins;
	uint64_t *a;
} ch_buf_t;

static inline void ch_insert_buf(ch_buf_t *buf, int p, uint64_t y) // insert a k-mer $y to a linear buffer
{
	int pre = y & ((1<<p) - 1);
	ch_buf_t *b = &buf[pre];
	if (b->n == b->m) {
		b->m = b->m < 8? 8 : b->m + (b->m>>1);
		REALLOC(b->a, b->m);
	}
	b->a[b->n++] = y;
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
				ch_insert_buf(buf, p, yak_hash_long(x));
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
					ch_insert_buf(buf, p, yak_hash_long(x));
				last = c;
			}
		} else l = 0, last = -1, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
	}
}

/******************
 * K-mer counting *
 ******************/

typedef struct { // global data structure for kt_pipeline()
	const yak_copt_t *opt;
	const void *flt_tab;
	int create_new, is_store;
	kseq_t *ks;
	yak_ch_t *h;
} pl_data_t;

typedef struct { // data structure for each step in kt_pipeline()
	pl_data_t *p;
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
	yak_ch_t *h = s->p->h;
	b->n_ins += yak_ch_insert_list(h, s->p->create_new, b->n, b->a);
}

static void worker_for_mz(void *data, long i, int tid)
{
	st_data_t *s = (st_data_t*)data;
	ha_mz1_v *b = &s->mz_buf[tid];
	s->mz_buf[tid].n = 0;
	ha_sketch(s->seq[i], s->len[i], s->p->opt->w, s->p->opt->k, 0, s->p->opt->is_HPC, b, s->p->flt_tab);
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
		while ((ret = kseq_read(p->ks)) >= 0) {
			int l = p->ks->seq.l;
			if (l < p->opt->k) continue;
			if (s->n_seq == s->m_seq) {
				s->m_seq = s->m_seq < 16? 16 : s->m_seq + (s->m_seq>>1);
				REALLOC(s->len, s->m_seq);
				REALLOC(s->seq, s->m_seq);
			}
			MALLOC(s->seq[s->n_seq], l);
			memcpy(s->seq[s->n_seq], p->ks->seq.s, l);
			s->len[s->n_seq++] = l;
			s->sum_len += l;
			s->nk += l - p->opt->k + 1;
			if (s->sum_len >= p->opt->chunk_size)
				break;
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
			MALLOC(s->buf[i].a, m);
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
			// compute minimizers
			CALLOC(s->mz, s->n_seq);
			CALLOC(s->mz_buf, p->opt->n_thread);
			kt_for(p->opt->n_thread, worker_for_mz, s, s->n_seq);
			for (i = 0; i < p->opt->n_thread; ++i)
				free(s->mz_buf[i].a);
			free(s->mz_buf);
			// insert minimizers
			for (i = 0; i < s->n_seq; ++i) {
				uint32_t j;
				for (j = 0; j < s->mz[i].n; ++j)
					ch_insert_buf(s->buf, p->opt->pre, s->mz[i].a[j].x);
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
			free(s->buf[i].a);
		}
		p->h->tot += n_ins;
		free(s->buf);
		fprintf(stderr, "[M::%s::%.3f*%.2f] processed %d sequences; %ld distinct k-mers in the hash table\n", __func__,
				yak_realtime(), yak_cputime() / yak_realtime(), s->n_seq, (long)p->h->tot);
		free(s);
	}
	return 0;
}

static yak_ch_t *yak_count(const char *fn, const yak_copt_t *opt, yak_ch_t *h0, const void *flt_tab)
{
	pl_data_t pl;
	gzFile fp;
	if ((fp = gzopen(fn, "r")) == 0) return 0;
	memset(&pl, 0, sizeof(pl_data_t));
	pl.ks = kseq_init(fp);
	pl.flt_tab = flt_tab;
	pl.opt = opt;
	if (h0) {
		pl.h = h0, pl.create_new = 0;
		assert(h0->k == opt->k && h0->pre == opt->pre);
	} else {
		pl.create_new = 1;
		pl.h = yak_ch_init(opt->k, opt->pre, opt->bf_n_hash, opt->bf_shift);
	}
	kt_pipeline(3, worker_count, &pl, 3);
	kseq_destroy(pl.ks);
	gzclose(fp);
	return pl.h;
}

static yak_ch_t *yak_count_file(const yak_copt_t *opt, int n_fn, char **fn, const void *flt_tab)
{
	int i;
	yak_ch_t *h = 0;
	for (i = 0; i < n_fn; ++i)
		h = yak_count(fn[i], opt, h, flt_tab);
	if (opt->bf_shift > 0)
		yak_ch_destroy_bf(h);
	return h;
}

static yak_hh_t *gen_hh(const yak_ch_t *h)
{
	int i;
	yak_hh_t *hh;
	hh = yak_hh_init();
	yak_hh_resize(hh, h->tot * 2);
	for (i = 0; i < 1<<h->pre; ++i) {
		yak_ht_t *ht = h->h[i].h;
		khint_t k;
		for (k = 0; k < kh_end(ht); ++k) {
			if (kh_exist(ht, k)) {
				uint64_t y = kh_key(ht, k) >> YAK_COUNTER_BITS << h->pre | i;
				int absent;
				yak_hh_put(hh, y, &absent);
			}
		}
	}
	return hh;
}

void *ha_gen_flt_tab(const hifiasm_opt_t *asm_opt)
{
	yak_hh_t *flt_tab;
	int64_t cnt[YAK_N_COUNTS];
	int peak_hom, peak_het, cutoff;
	yak_copt_t opt;
	yak_ch_t *h;
	yak_copt_init(&opt);
	opt.is_HPC = !asm_opt->no_HPC;
	opt.k = asm_opt->k_mer_length;
	opt.n_thread = asm_opt->thread_num;
	opt.bf_shift = asm_opt->bf_shift;
	h = yak_count_file(&opt, asm_opt->num_reads, asm_opt->read_file_names, 0);
	yak_ch_hist(h, cnt, opt.n_thread);
	peak_hom = yak_analyze_count(YAK_N_COUNTS, cnt, &peak_het);
	if (peak_hom > 0) fprintf(stderr, "[M::%s] peak_hom: %d; peak_het: %d\n", __func__, peak_hom, peak_het);
	cutoff = (int)(peak_hom * asm_opt->high_factor);
	if (cutoff > YAK_MAX_COUNT - 1) cutoff = YAK_MAX_COUNT - 1;
	yak_ch_shrink(h, cutoff, YAK_MAX_COUNT, opt.n_thread);
	flt_tab = gen_hh(h);
	yak_ch_destroy(h);
	fprintf(stderr, "[M::%s] filtered out %ld k-mers occurring %d or more times\n",
			__func__, (long)kh_size(flt_tab), cutoff);
	return (void*)flt_tab;
}

int ha_hf_isflt(const void *hh, uint64_t y)
{
	yak_hh_t *h = (yak_hh_t*)hh;
	khint_t k;
	k = yak_hh_get(h, y);
	return k == kh_end(h)? 0 : 1;
}

void ha_hf_destroy(void *h)
{
	yak_hh_destroy((yak_hh_t*)h);
}
