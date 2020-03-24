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

typedef struct {
	int32_t bf_shift, bf_n_hash;
	int32_t k;
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

typedef struct { // global data structure for kt_pipeline()
	const yak_copt_t *opt;
	int create_new;
	kseq_t *ks;
	yak_ch_t *h;
} pldat_t;

typedef struct { // data structure for each step in kt_pipeline()
	pldat_t *p;
	int n, m, sum_len, nk;
	int *len;
	char **seq;
	ch_buf_t *buf;
} stepdat_t;

static void worker_for(void *data, long i, int tid) // callback for kt_for()
{
	stepdat_t *s = (stepdat_t*)data;
	ch_buf_t *b = &s->buf[i];
	yak_ch_t *h = s->p->h;
	b->n_ins += yak_ch_insert_list(h, s->p->create_new, b->n, b->a);
}

static void *worker_count_all(void *data, int step, void *in) // callback for kt_pipeline()
{
	pldat_t *p = (pldat_t*)data;
	if (step == 0) { // step 1: read a block of sequences
		int ret;
		stepdat_t *s;
		CALLOC(s, 1);
		s->p = p;
		while ((ret = kseq_read(p->ks)) >= 0) {
			int l = p->ks->seq.l;
			if (l < p->opt->k) continue;
			if (s->n == s->m) {
				s->m = s->m < 16? 16 : s->m + (s->n>>1);
				REALLOC(s->len, s->m);
				REALLOC(s->seq, s->m);
			}
			MALLOC(s->seq[s->n], l);
			memcpy(s->seq[s->n], p->ks->seq.s, l);
			s->len[s->n++] = l;
			s->sum_len += l;
			s->nk += l - p->opt->k + 1;
			if (s->sum_len >= p->opt->chunk_size)
				break;
		}
		if (s->sum_len == 0) free(s);
		else return s;
	} else if (step == 1) { // step 2: extract k-mers
		stepdat_t *s = (stepdat_t*)in;
		int i, n = 1<<p->opt->pre, m;
		CALLOC(s->buf, n);
		m = (int)(s->nk * 1.2 / n) + 1;
		for (i = 0; i < n; ++i) {
			s->buf[i].m = m;
			MALLOC(s->buf[i].a, m);
		}
		for (i = 0; i < s->n; ++i) {
			count_seq_buf(s->buf, p->opt->k, p->opt->pre, s->len[i], s->seq[i]);
			free(s->seq[i]);
		}
		free(s->seq); free(s->len);
		return s;
	} else if (step == 2) { // step 3: insert k-mers to hash table
		stepdat_t *s = (stepdat_t*)in;
		int i, n = 1<<p->opt->pre;
		uint64_t n_ins = 0;
		kt_for(p->opt->n_thread, worker_for, s, n);
		for (i = 0; i < n; ++i) {
			n_ins += s->buf[i].n_ins;
			free(s->buf[i].a);
		}
		p->h->tot += n_ins;
		free(s->buf);
		fprintf(stderr, "[M::%s::%.3f*%.2f] processed %d sequences; %ld distinct k-mers in the hash table\n", __func__,
				yak_realtime(), yak_cputime() / yak_realtime(), s->n, (long)p->h->tot);
		free(s);
	}
	return 0;
}

static yak_ch_t *yak_count(const char *fn, const yak_copt_t *opt, yak_ch_t *h0)
{
	pldat_t pl;
	gzFile fp;
	if ((fp = gzopen(fn, "r")) == 0) return 0;
	pl.ks = kseq_init(fp);
	pl.opt = opt;
	if (h0) {
		pl.h = h0, pl.create_new = 0;
		assert(h0->k == opt->k && h0->pre == opt->pre);
	} else {
		pl.create_new = 1;
		pl.h = yak_ch_init(opt->k, opt->pre, opt->bf_n_hash, opt->bf_shift);
	}
	kt_pipeline(3, worker_count_all, &pl, 3);
	kseq_destroy(pl.ks);
	gzclose(fp);
	return pl.h;
}

static yak_ch_t *yak_count_file(const yak_copt_t *opt, int n_fn, char **fn)
{
	int i;
	yak_ch_t *h = 0;
	for (i = 0; i < n_fn; ++i)
		h = yak_count(fn[i], opt, h);
	if (opt->bf_shift > 0)
		yak_ch_destroy_bf(h);
	return h;
}

static void yak_hist_line(int c, int x, int exceed, int64_t cnt)
{
	int j;
	if (c >= 0) fprintf(stderr, "[M::%s] %5d: ", __func__, c);
	else fprintf(stderr, "[M::%s] %5s: ", __func__, "rest");
	for (j = 0; j < x; ++j) fputc('*', stderr);
	if (exceed) fputc('>', stderr);
	fprintf(stderr, " %lld\n", (long long)cnt);
}

int yak_analyze_count(int n_cnt, const int64_t *cnt, int *peak_het)
{
	const int hist_max = 100;
	int i, low_i, max_i, max2_i, max3_i;
	int64_t max, max2, max3, min;

	// find the low point from the left
	*peak_het = -1;
	low_i = 2;
	for (i = 3; i < n_cnt; ++i)
		if (cnt[i] > cnt[i-1]) break;
	low_i = i - 1;
	fprintf(stderr, "[M::%s] lowest: count[%d] = %ld\n", __func__, low_i, (long)cnt[low_i]);
	if (low_i == n_cnt - 1) return -1; // low coverage

	// find the highest peak
	max_i = low_i + 1, max = cnt[max_i];
	for (i = low_i + 1; i < n_cnt; ++i)
		if (cnt[i] > max)
			max = cnt[i], max_i = i;
	fprintf(stderr, "[M::%s] highest: count[%d] = %ld\n", __func__, max_i, (long)cnt[max_i]);

	// print histogram
	for (i = 2; i < n_cnt; ++i) {
		int x, exceed = 0;
		x = (int)((double)hist_max * cnt[i] / cnt[max_i] + .499);
		if (x > hist_max) exceed = 1, x = hist_max; // may happen if cnt[2] is higher
		if (i > max_i && x == 0) break;
		yak_hist_line(i, x, exceed, cnt[i]);
	}
	{
		int x, exceed = 0;
		int64_t rest = 0;
		for (; i < n_cnt; ++i) rest += cnt[i];
		x = (int)((double)hist_max * rest / cnt[max_i] + .499);
		if (x > hist_max) exceed = 1, x = hist_max;
		yak_hist_line(-1, x, exceed, rest);
	}

	// look for smaller peak on the low end
	max2 = -1; max2_i = -1;
	for (i = max_i - 1; i > low_i; --i) {
		if (cnt[i] >= cnt[i-1] && cnt[i] >= cnt[i+1]) {
			if (cnt[i] > max2) max2 = cnt[i], max2_i = i;
		}
	}
	if (max2_i > low_i && max2_i < max_i) {
		for (i = max2_i + 1, min = max; i < max_i; ++i)
			if (cnt[i] < min) min = cnt[i];
		if (max2 < max * 0.05 || min > max2 * 0.95)
			max2 = -1, max2_i = -1;
	}
	if (max2 > 0) fprintf(stderr, "[M::%s] left: count[%d] = %ld\n", __func__, max2_i, (long)cnt[max2_i]);
	else fprintf(stderr, "[M::%s] left: none\n", __func__);

	// look for smaller peak on the high end
	max3 = -1; max3_i = -1;
	for (i = max_i + 1; i < n_cnt - 1; ++i) {
		if (cnt[i] >= cnt[i-1] && cnt[i] >= cnt[i+1]) {
			if (cnt[i] > max3) max3 = cnt[i], max3_i = i;
		}
	}
	if (max3_i > max_i) {
		for (i = max_i + 1, min = max; i < max3_i; ++i)
			if (cnt[i] < min) min = cnt[i];
		if (max3 < max * 0.05 || min > max3 * 0.95 || max3_i > max_i * 2.5)
			max3 = -1, max3_i = -1;
	}
	if (max3 > 0) fprintf(stderr, "[M::%s] right: count[%d] = %ld\n", __func__, max3_i, (long)cnt[max3_i]);
	else fprintf(stderr, "[M::%s] right: none\n", __func__);
	if (max3_i > 0) {
		*peak_het = max_i;
		return max3_i;
	} else {
		if (max2_i > 0) *peak_het = max2_i;
		return max_i;
	}
}

void ha_count_high(const hifiasm_opt_t *asm_opt)
{
	int64_t cnt[YAK_N_COUNTS];
	int peak_hom, peak_het;
	yak_copt_t opt;
	yak_ch_t *h;
	yak_copt_init(&opt);
	opt.k = asm_opt->k_mer_length;
	opt.n_thread = asm_opt->thread_num;
	opt.bf_shift = asm_opt->bf_shift;
	h = yak_count_file(&opt, asm_opt->num_reads, asm_opt->read_file_names);
	yak_ch_hist(h, cnt, opt.n_thread);
	peak_hom = yak_analyze_count(YAK_N_COUNTS, cnt, &peak_het);
	if (peak_hom > 0) fprintf(stderr, "[M::%s] peak_hom: %d; peak_het: %d\n", __func__, peak_hom, peak_het);
	yak_ch_destroy(h);
}
