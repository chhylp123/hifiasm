#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <zlib.h>
#include "khashl.h" // hash table
#include "kthread.h"
#include "kseq.h"
#include "Process_Read.h"
#include "htab.h"
#include "CommandLines.h"

#define YAK_MAX_KMER     31
#define YAK_COUNTER_BITS 10 // yak uses 10, but hifiasm uses 12; we have to copy over some yak code here due to this
#define YAK_N_COUNTS     (1<<YAK_COUNTER_BITS)
#define YAK_MAX_COUNT    ((1<<YAK_COUNTER_BITS)-1)

#define YAK_LOAD_ALL       1
#define YAK_LOAD_TRIOBIN1  2
#define YAK_LOAD_TRIOBIN2  3

#define YAK_MAGIC "YAK\2"

#define yak_ch_eq(a, b) ((a)>>YAK_COUNTER_BITS == (b)>>YAK_COUNTER_BITS) // lower 8 bits for counts; higher bits for k-mer
#define yak_ch_hash(a) ((a)>>YAK_COUNTER_BITS)
KHASHL_SET_INIT(static klib_unused, yak_ht_t, yak_ht, uint64_t, yak_ch_hash, yak_ch_eq)

typedef const char *ha_cstr_t;
KHASHL_MAP_INIT(static klib_unused, cstr_ht_t, cstr_ht, ha_cstr_t, int64_t, kh_hash_str, kh_eq_str)

KSTREAM_INIT(gzFile, gzread, 65536)

typedef struct {
	struct yak_ht_t *h;
} yak_ch1_t;

typedef struct {
	int k, pre, n_hash, n_shift;
	uint64_t tot;
	yak_ch1_t *h;
} yak_ch_t;

static int yak_ch_get(const yak_ch_t *h, uint64_t x)
{
	int mask = (1<<h->pre) - 1;
	yak_ht_t *g = h->h[x&mask].h;
	khint_t k;
	k = yak_ht_get(g, x >> h->pre << YAK_COUNTER_BITS);
	return k == kh_end(g)? -1 : kh_key(g, k)&YAK_MAX_COUNT;
}

static yak_ch_t *yak_ch_init(int k, int pre)
{
	yak_ch_t *h;
	int i;
	if (pre < YAK_COUNTER_BITS) return 0;
	CALLOC(h, 1);
	h->k = k, h->pre = pre;
	CALLOC(h->h, 1<<h->pre);
	for (i = 0; i < 1<<h->pre; ++i)
		h->h[i].h = yak_ht_init();
	return h;
}

static yak_ch_t *yak_ch_restore_core(yak_ch_t *ch0, const char *fn, int mode, ...)
{
	va_list ap;
	FILE *fp;
	uint32_t t[3], f_tmp = 0;
	char magic[4];
	int i, j, absent, min_cnt = 0, mid_cnt = 0, mode_err = 0;
	uint64_t mask = (1ULL<<YAK_COUNTER_BITS) - 1, n_ins = 0, n_new = 0;
	yak_ch_t *ch;

	va_start(ap, mode);
	if (mode == YAK_LOAD_ALL) { // do nothing
	} else if (mode == YAK_LOAD_TRIOBIN1 || mode == YAK_LOAD_TRIOBIN2) {
		assert(YAK_COUNTER_BITS >= 4);
		min_cnt = va_arg(ap, int);
		mid_cnt = va_arg(ap, int);
		if (ch0 == 0 && mode == YAK_LOAD_TRIOBIN2)
			mode_err = 1;
	} else mode_err = 1;
	va_end(ap);
	if (mode_err) return 0;

	if ((fp = fopen(fn, "rb")) == 0) return 0;
	if (fread(magic, 1, 4, fp) != 4) return 0;
	if (strncmp(magic, YAK_MAGIC, 4) != 0) {
		fprintf(stderr, "ERROR: wrong file magic.\n");
		fclose(fp);
		return 0;
	}
	f_tmp += fread(t, 4, 3, fp);
	if (t[2] != YAK_COUNTER_BITS) {
		fprintf(stderr, "ERROR: saved counter bits: %d; compile-time counter bits: %d\n", t[2], YAK_COUNTER_BITS);
		fclose(fp);
		return 0;
	}
	///t[0] = k; t[1] = pre, t[2] = YAK_COUNTER_BITS;
	ch = ch0 == 0? yak_ch_init(t[0], t[1]) : ch0;
	assert((int)t[0] == ch->k && (int)t[1] == ch->pre);
	for (i = 0; i < 1<<ch->pre; ++i) {
		yak_ht_t *h = ch->h[i].h;
		f_tmp += fread(t, 4, 2, fp);
		///t[0] = kh_capacity(h), t[1] = kh_size(h);
		if (ch0 == 0) yak_ht_resize(h, t[0]);
		for (j = 0; j < (int)t[1]; ++j) {
			uint64_t key;
			f_tmp += fread(&key, 8, 1, fp);
			if (mode == YAK_LOAD_ALL) {
				++n_ins;
				yak_ht_put(h, key, &absent);
				if (absent) ++n_new;
			} else if (mode == YAK_LOAD_TRIOBIN1 || mode == YAK_LOAD_TRIOBIN2) {
				int cnt = key & mask, x, shift = mode == YAK_LOAD_TRIOBIN1? 0 : 2;
				//1. filter singleton k-mer; 2. label non-repeat and repeat
				if (cnt >= mid_cnt) x = 2<<shift;
				else if (cnt >= min_cnt) x = 1<<shift;
				else x = -1;
				if (x >= 0) {
					khint_t k;
					///no need cnt at all
					key = (key & ~mask) | x;
					++n_ins;
					k = yak_ht_put(h, key, &absent);
					if (absent) ++n_new;
					else kh_key(h, k) = kh_key(h, k) | x;
				}
			}
		}
	}
	fclose(fp);
	///fprintf(stderr, "[M::%s] inserted %ld k-mers, of which %ld are new\n", __func__, (long)n_ins, (long)n_new);
	return ch;
}

static void yak_ch_destroy(yak_ch_t *h)
{
	int i;
	if (h == 0) return;
	for (i = 0; i < 1<<h->pre; ++i)
		yak_ht_destroy(h->h[i].h);
	free(h->h); free(h);
}

typedef struct {
	int max;
	uint32_t *s;
} tb_buf_t;

typedef struct {
	int k, n_threads, print_diff;
	double ratio_thres;
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
} tb_step_t;

static char tb_classify(const int sc[2], const int *c, int k, double ratio_thres)
{
	char type;
	if (sc[0] == 0 && sc[1] == 0) {
		if (c[0<<2|2] == c[2<<2|0]) type = '0';
		else if (c[0<<2|2] >= k - 4 + c[2<<2|0] && (c[2<<2|0] <= 1 || c[0<<2|2] * 0.05 > c[2<<2|0])) type = 'p';
		else if (c[2<<2|0] >= k - 4 + c[0<<2|2] && (c[0<<2|2] <= 1 || c[2<<2|0] * 0.05 > c[0<<2|2])) type = 'm';
		else type = '0';
	} else if (sc[0] > k && sc[1] > k) {
		type = 'a';
	} else if (sc[0] >= k - 4 + sc[1] && sc[0] * 0.05 >= sc[1] && c[0<<2|2] * ratio_thres > c[2<<2|0]) {
		type = 'p';
	} else if (sc[1] >= k - 4 + sc[0] && sc[1] * 0.05 >= sc[0] && c[2<<2|0] * ratio_thres > c[0<<2|2]) {
		type = 'm';
	} else {
		type = 'a';
	}
	return type;
}

static void tb_worker(void *_data, long k, int tid)
{
	tb_shared_t *aux = (tb_shared_t*)_data;
	UC_Read *s = &aux->bseq[tid]; 
	recover_UC_Read(s, aux->seq, k);
	tb_buf_t *b = &aux->buf[tid];
	tb_cnt_t cnt; memset(&cnt, 0, sizeof(tb_cnt_t));
	uint64_t x[4], mask;
	int i, l, shift;
	if (aux->ch->k < 32) {
		mask = (1ULL<<2*aux->ch->k) - 1;
		shift = 2 * (aux->ch->k - 1);
	} else {
		mask = (1ULL<<aux->ch->k) - 1;
		shift = aux->ch->k - 1;
	}
	if (s->length > b->max) {
		b->max = s->length;
		kroundup32(b->max);
		b->s = (uint32_t*)realloc(b->s, b->max * sizeof(uint32_t));
	}
	memset(b->s, 0, s->length * sizeof(uint32_t));
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < s->length; ++i) {
		int flag, c = seq_nt4_table[(uint8_t)s->seq[i]];
		if (c < 4) {
			if (aux->ch->k < 32) {
				x[0] = (x[0] << 2 | c) & mask;
				x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
			} else {
				x[0] = (x[0] << 1 | (c&1))  & mask;
				x[1] = (x[1] << 1 | (c>>1)) & mask;
				x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
				x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
			}
			if (++l >= aux->k) {
				int type = 0, c1, c2;
				uint64_t y;
				++cnt.nk;
				if (aux->ch->k < 32)
					y = yak_hash64(x[0] < x[1]? x[0] : x[1], mask);
				else
					y = yak_hash_long(x);
				flag = yak_ch_get(aux->ch, y);
				if (flag < 0) flag = 0;
				c1 = flag&3, c2 = flag>>2&3;
				if (c1 == 2 && c2 == 0) type = 1;
				else if (c2 == 2 && c1 == 0) type = 2;
				b->s[i] = type;
				++cnt.c[flag];
			}
		} else l = 0, x[0] = x[1] = x[2] = x[3] = 0;
	}
	for (l = 0, i = 1; i <= s->length; ++i) {
		if (i == s->length || b->s[i] != b->s[l]) {
			if (b->s[l] > 0 && i - l >= aux->k - 4)
				cnt.sc[b->s[l] - 1] += i - l;
			l = i;
		}
	}

	int *c = cnt.c;
	char type;
	type = tb_classify(cnt.sc, c, aux->k, aux->ratio_thres);
	aux->seq->trio_flag[k] = AMBIGU;
	if(type == 'p') aux->seq->trio_flag[k] = FATHER;
	if(type == 'm') aux->seq->trio_flag[k] = MOTHER;
}

static void ha_triobin_yak(const hifiasm_opt_t *opt)
{
    yak_ch_t *ch;
    int i /**, min_cnt = 2, mid_cnt = 5**/;
    tb_shared_t aux;
    memset(&aux, 0, sizeof(tb_shared_t));
	aux.n_threads = opt->thread_num, aux.print_diff = 0;
	aux.ratio_thres = 0.33;
	aux.seq = &R_INF;

    ch = yak_ch_restore_core(0,  opt->fn_bin_yak[0], YAK_LOAD_TRIOBIN1, opt->min_cnt, opt->mid_cnt);
	ch = yak_ch_restore_core(ch, opt->fn_bin_yak[1], YAK_LOAD_TRIOBIN2, opt->min_cnt, opt->mid_cnt);

    aux.k = ch->k;
    aux.ch = ch;
	aux.buf = (tb_buf_t*)calloc(aux.n_threads, sizeof(tb_buf_t));
	aux.bseq = (UC_Read*)calloc(aux.n_threads, sizeof(UC_Read));
	for (i = 0; i < aux.n_threads; ++i)
		init_UC_Read(&aux.bseq[i]);
	
	kt_for(aux.n_threads, tb_worker, &aux, aux.seq->total_reads);
	
    for (i = 0; i < aux.n_threads; ++i) {
		free(aux.buf[i].s);
		destory_UC_Read(&aux.bseq[i]);
	} 
	free(aux.buf);
	free(aux.bseq);
    yak_ch_destroy(ch);

	fprintf(stderr, "[M::%s::%.3f*%.2f] ==> partitioned reads using yak dumps\n", __func__, yak_realtime(), yak_cpu_usage());
}

static int ha_triobin_set_list(const cstr_ht_t *h, const char *fn, int flag)
{
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	int dret;
	int64_t n_tot = 0, n_bin = 0;
	fp = gzopen(fn, "r");
	if (fp == 0) {
		fprintf(stderr, "ERROR: failed to open file '%s'\n", fn);
		return -1;
	}
	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		char *p;
		khint_t k;
		++n_tot;
		for (p = str.s; *p; ++p)
			if (*p == '\t' || *p == ' ')
				*p = 0;
		k = cstr_ht_get(h, str.s);
		if (k != kh_end(h)) {
			R_INF.trio_flag[kh_val(h, k)] = flag;
			++n_bin;
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	fprintf(stderr, "[M::%s::%.3f*%.2f] flagged %ld reads, out of %ld lines in file '%s'\n",
			__func__, yak_realtime(), yak_cpu_usage(), (long)n_bin, (long)n_tot, fn);
	return 0;
}

static void ha_triobin_list(const hifiasm_opt_t *opt)
{
	int64_t i;
	khint_t k;
	cstr_ht_t *h;
	assert(R_INF.total_reads < (uint32_t)-1);
	h = cstr_ht_init();
	for (i = 0; i < (int64_t)R_INF.total_reads; ++i) {
		int absent;
		char *str = (char*)calloc(Get_NAME_LENGTH(R_INF, i) + 1, 1);
		strncpy(str, Get_NAME(R_INF, i), Get_NAME_LENGTH(R_INF, i));
		k = cstr_ht_put(h, str, &absent);
		if (absent) kh_val(h, k) = i;
	}
	fprintf(stderr, "[M::%s::%.3f*%.2f] created the hash table for read names\n", __func__, yak_realtime(), yak_cpu_usage());
	ha_triobin_set_list(h, opt->fn_bin_list[0], FATHER);
	ha_triobin_set_list(h, opt->fn_bin_list[1], MOTHER);
	for (k = 0; k < kh_end(h); ++k)
		if (kh_exist(h, k))
			free((char*)kh_key(h, k));
	cstr_ht_destroy(h);
	fprintf(stderr, "[M::%s::%.3f*%.2f] ==> partitioned reads with external lists\n", __func__, yak_realtime(), yak_cpu_usage());
}

void ha_triobin(const hifiasm_opt_t *opt)
{
	memset(R_INF.trio_flag, AMBIGU, R_INF.total_reads * sizeof(uint8_t));
	if (opt->fn_bin_list[0] && opt->fn_bin_list[1])
		ha_triobin_list(opt);
	if (opt->fn_bin_yak[0] && opt->fn_bin_yak[1])
		ha_triobin_yak(opt);
}
