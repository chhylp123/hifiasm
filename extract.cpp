#include <zlib.h>
#include <string.h>
#include "Process_Read.h"
#include "khashl.h"
#include "kseq.h"

typedef const char *cstr_t;
KHASHL_CSET_INIT(KH_LOCAL, strset_t, ss, cstr_t, kh_hash_str, kh_eq_str)
KHASHL_MAP_INIT(KH_LOCAL, hm64_t, h64, uint64_t, int, kh_hash_uint64, kh_eq_generic)
KSTREAM_INIT(gzFile, gzread, 65536)

#define GFA_MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define GFA_REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

char *gfa_strdup(const char *src)
{
	int32_t len;
	char *dst;
	len = strlen(src);
	GFA_MALLOC(dst, len + 1);
	memcpy(dst, src, len + 1);
	return dst;
}

char *gfa_strndup(const char *src, size_t n)
{
	char *dst;
	GFA_MALLOC(dst, n + 1);
	strncpy(dst, src, n);
	dst[n] = 0;
	return dst;
}

char **gv_read_list(const char *o, int *n_)
{
	int n = 0, m = 0;
	char **s = 0;
	*n_ = 0;
	if (*o != '@') {
		const char *q = o, *p;
		for (p = q;; ++p) {
			if (*p == ',' || *p == 0) {
				if (n == m) {
					m = m? m<<1 : 16;
					GFA_REALLOC(s, m);
				}
				s[n++] = gfa_strndup(q, p - q);
				if (*p == 0) break;
				q = p + 1;
			}
		}
	} else {
		gzFile fp;
		kstream_t *ks;
		kstring_t str = {0,0,0};
		int dret;

		fp = gzopen(o + 1, "r");
		if (fp == 0) return 0;
		ks = ks_init(fp);
		while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
			char *p;
			for (p = str.s; *p && !isspace(*p); ++p);
			if (n == m) {
				m = m? m<<1 : 16;
				GFA_REALLOC(s, m);
			}
			s[n++] = gfa_strndup(str.s, p - str.s);
		}
		ks_destroy(ks);
		gzclose(fp);
	}
	if (s) s = (char**)realloc(s, n * sizeof(char*));
	*n_ = n;
	return s;
}

void ha_extract_print(const All_reads *rs, int n_rounds, int n, char **list)
{
	hm64_t *h;
	khint_t k;
	int i, absent, m, l;
	uint64_t j;
	const ma_hit_t_alloc *ov[2] = { rs->paf, rs->reverse_paf };
	FILE *fp = stdout;

	if (n > 0) {
		int max_len = 0;
		char *s = 0;
		strset_t *ss;
		ss = ss_init();
		for (i = 0; i < n; ++i)
			ss_put(ss, list[i], &absent);
		for (j = 0; j < rs->total_reads; ++j)
			if (max_len < (int)Get_NAME_LENGTH(*rs, j))
				max_len = Get_NAME_LENGTH(*rs, j);
		GFA_MALLOC(s, max_len + 1);
		h = h64_init();
		for (j = 0; j < rs->total_reads; ++j) {
			strncpy(s, Get_NAME(*rs, j), Get_NAME_LENGTH(*rs, j));
			s[Get_NAME_LENGTH(*rs, j)] = 0;
			if (ss_get(ss, s) != kh_end(ss)) {
				k = h64_put(h, j, &absent);
				kh_val(h, k) = -1;
			}
		}
		free(s);
		ss_destroy(ss);
	} else return;

	for (m = 0; m < n_rounds; ++m) {
		for (j = 0; j < rs->total_reads; ++j) {
			for (l = 0; l < 2; ++l) {
				const ma_hit_t_alloc *o = &ov[l][j];
				for (i = 0; i < (int)o->length; ++i) {
					uint64_t q = Get_qn(o->buffer[i]);
					uint64_t t = Get_tn(o->buffer[i]);
					int q_hit = 0, t_hit = 0;
					k = h64_get(h, q);
					q_hit = (k < kh_end(h) && kh_val(h, k) < m);
					k = h64_get(h, t);
					t_hit = (k < kh_end(h) && kh_val(h, k) < m);
					if ((!q_hit && !t_hit) || (q_hit && t_hit)) continue;
					if (!q_hit) {
						k = h64_put(h, q, &absent);
						if (absent) kh_val(h, k) = m;
					}
					if (!t_hit) {
						k = h64_put(h, t, &absent);
						if (absent) kh_val(h, k) = m;
					}
				}
			}
		}
	}

	for (j = 0; j < rs->total_reads; ++j) {
		for (l = 0; l < 2; ++l) {
			const ma_hit_t_alloc *o = &ov[l][j];
			for (i = 0; i < (int)o->length; ++i) {
				uint64_t q = Get_qn(o->buffer[i]);
				uint64_t t = Get_tn(o->buffer[i]);
				int q_hit = 0, t_hit = 0;
				q_hit = (h64_get(h, q) < kh_end(h));
				t_hit = (h64_get(h, t) < kh_end(h));
				if (!q_hit && !t_hit) continue;
				fwrite(Get_NAME(*rs, q), 1, Get_NAME_LENGTH(*rs, q), fp);
				fwrite("\t", 1, 1, fp);
				fprintf(fp, "%lu\t", (unsigned long)Get_READ_LENGTH(*rs, q));
				fprintf(fp, "%d\t", Get_qs(o->buffer[i]));
				fprintf(fp, "%d\t", Get_qe(o->buffer[i]));
				fputs(o->buffer[i].rev? "-\t" : "+\t", fp);
				fwrite(Get_NAME(*rs, t), 1, Get_NAME_LENGTH(*rs, t), fp);
				fwrite("\t", 1, 1, fp);
				fprintf(fp, "%lu\t", (unsigned long)Get_READ_LENGTH(*rs, t));
				fprintf(fp, "%d\t", Get_ts(o->buffer[i]));
				fprintf(fp, "%d\t%d\t%d\t%d\n", Get_te(o->buffer[i]), o->buffer[i].ml, o->buffer[i].bl, !l);
			}
		}
	}

	h64_destroy(h);
}

void ha_extract_print_list(const All_reads *rs, int n_rounds, const char *o)
{
	int i, n;
	char **list;
	list = gv_read_list(o, &n);
	ha_extract_print(rs, n_rounds, n, list);
	for (i = 0; i < n; ++i) free(list[i]);
	free(list);
}
