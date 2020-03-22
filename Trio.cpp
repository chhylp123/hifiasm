#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include "khashl.h" // hash table
#include "kthread.h"
#include "Process_Read.h"
#include "Trio.h"
#include "CommandLines.h"

#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))
#define yak_ch_eq(a, b) ((a)>>YAK_COUNTER_BITS == (b)>>YAK_COUNTER_BITS) // lower 8 bits for counts; higher bits for k-mer
#define yak_ch_hash(a) ((a)>>YAK_COUNTER_BITS)
KHASHL_SET_INIT(, yak_ht_t, yak_ht, uint64_t, yak_ch_hash, yak_ch_eq)

///#define CHUNK_SIZE 200000

unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
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

int yak_ch_get(const yak_ch_t *h, uint64_t x)
{
	int mask = (1<<h->pre) - 1;
	yak_ht_t *g = h->h[x&mask].h;
	khint_t k;
	k = yak_ht_get(g, x >> h->pre << YAK_COUNTER_BITS);
	return k == kh_end(g)? -1 : kh_key(g, k)&YAK_MAX_COUNT;
}

yak_bf_t *yak_bf_init(int n_shift, int n_hashes)
{
	yak_bf_t *b;
	void *ptr = 0;
	if (n_shift + YAK_BLK_SHIFT > 64 || n_shift < YAK_BLK_SHIFT) return 0;
	b = (yak_bf_t*)calloc(1, sizeof(yak_bf_t));
	b->n_shift = n_shift;
	b->n_hashes = n_hashes;
	posix_memalign(&ptr, 1<<(YAK_BLK_SHIFT-3), 1ULL<<(n_shift-3));
	b->b = (uint8_t*)ptr;
	bzero(b->b, 1ULL<<(n_shift-3));
	return b;
}

yak_ch_t *yak_ch_init(int k, int pre, int n_hash, int n_shift)
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

void yak_bf_destroy(yak_bf_t *b)
{
	if (b == 0) return;
	free(b->b); free(b);
}

void yak_ch_destroy_bf(yak_ch_t *h)
{
	int i;
	for (i = 0; i < 1<<h->pre; ++i) {
		if (h->h[i].b)
			yak_bf_destroy(h->h[i].b);
		h->h[i].b = 0;
	}
}


yak_ch_t *yak_ch_restore_core(yak_ch_t *ch0, const char *fn, int mode, ...)
{
	va_list ap;
	FILE *fp;
	uint32_t t[3];
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
	fread(t, 4, 3, fp);
	if (t[2] != YAK_COUNTER_BITS) {
		fprintf(stderr, "ERROR: saved counter bits: %d; compile-time counter bits: %d\n", t[2], YAK_COUNTER_BITS);
		fclose(fp);
		return 0;
	}

	ch = ch0 == 0? yak_ch_init(t[0], t[1], 0, 0) : ch0;
	assert((int)t[0] == ch->k && (int)t[1] == ch->pre);
	for (i = 0; i < 1<<ch->pre; ++i) {
		yak_ht_t *h = ch->h[i].h;
		fread(t, 4, 2, fp);
		if (ch0 == 0) yak_ht_resize(h, t[0]);
		for (j = 0; j < t[1]; ++j) {
			uint64_t key;
			fread(&key, 8, 1, fp);
			if (mode == YAK_LOAD_ALL) {
				++n_ins;
				yak_ht_put(h, key, &absent);
				if (absent) ++n_new;
			} else if (mode == YAK_LOAD_TRIOBIN1 || mode == YAK_LOAD_TRIOBIN2) {
				int cnt = key & mask, x, shift = mode == YAK_LOAD_TRIOBIN1? 0 : 2;
				if (cnt >= mid_cnt) x = 2<<shift;
				else if (cnt >= min_cnt) x = 1<<shift;
				else x = -1;
				if (x >= 0) {
					khint_t k;
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


void yak_ch_destroy(yak_ch_t *h)
{
	int i;
	if (h == 0) return;
	yak_ch_destroy_bf(h);
	for (i = 0; i < 1<<h->pre; ++i)
		yak_ht_destroy(h->h[i].h);
	free(h->h); free(h);
}


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
	///tb_step_t *t = (tb_step_t*)_data;
	///tb_shared_t *aux = t->aux;
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


				//++t->cnt[k].nk;
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


				///++t->cnt[k].c[flag];
				++cnt.c[flag];

				// if (aux->print_diff && (flag>>2&3) != (flag&3))
				// 	printf("D\t%s\t%d\t%d\t%d\n", s->name, i, flag&3, flag>>2&3);
			}
		} else l = 0, x[0] = x[1] = x[2] = x[3] = 0;
	}
	for (l = 0, i = 1; i <= s->length; ++i) {
		if (i == s->length || b->s[i] != b->s[l]) {
			if (b->s[l] > 0 && i - l >= aux->k - 4)
			{
				///t->cnt[k].sc[b->s[l] - 1] += i - l;
				cnt.sc[b->s[l] - 1] += i - l;
			}
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



void trio_partition()
{
    if(asm_opt.pat_index == NULL || asm_opt.mat_index == NULL)
	{
		memset(R_INF.trio_flag, AMBIGU, R_INF.total_reads*sizeof(uint8_t));
		return;
	}

	double start_time = Get_T();
	fprintf(stderr, "Start trio binning ...... \n");

    yak_ch_t *ch;
    int i/**, min_cnt = 2, mid_cnt = 5**/;
    tb_shared_t aux;
    memset(&aux, 0, sizeof(tb_shared_t));
	aux.n_threads = asm_opt.thread_num, aux.print_diff = 0;
	aux.ratio_thres = 0.33;
	aux.seq = &R_INF;


    ch = yak_ch_restore_core(0,  asm_opt.pat_index, YAK_LOAD_TRIOBIN1, asm_opt.min_cnt, asm_opt.mid_cnt);
	ch = yak_ch_restore_core(ch, asm_opt.mat_index, YAK_LOAD_TRIOBIN2, asm_opt.min_cnt, asm_opt.mid_cnt);



    aux.k = ch->k;
    aux.ch = ch;
	aux.buf = (tb_buf_t*)calloc(aux.n_threads, sizeof(tb_buf_t));
	aux.bseq = (UC_Read*)calloc(aux.n_threads, sizeof(UC_Read));
	for (i = 0; i < aux.n_threads; ++i)
	{
		init_UC_Read(&aux.bseq[i]);
	} 
	
	kt_for(aux.n_threads, tb_worker, &aux, aux.seq->total_reads);
	
    for (i = 0; i < aux.n_threads; ++i)
	{
		free(aux.buf[i].s);
		destory_UC_Read(&aux.bseq[i]);
	} 
	free(aux.buf);
	free(aux.bseq);
    yak_ch_destroy(ch);

	fprintf(stderr, "Trio binning has been done.\n");
	fprintf(stderr, "%-30s%18.2f\n\n", "Trio binning time:", Get_T() - start_time); 
}