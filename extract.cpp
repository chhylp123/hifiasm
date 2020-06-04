#include <zlib.h>
#include <string.h>
#include "Process_Read.h"
#include "khashl.h"
#include "kseq.h"

typedef const char *cstr_t;
KHASHL_CMAP_INIT(KH_LOCAL, strmap_t, ss, cstr_t, int, kh_hash_str, kh_eq_str)
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
	strmap_t *h;
	khint_t k;
	int i, absent, m, l, max_len = 0;
	uint64_t j;
	char *s = 0;
	const ma_hit_t_alloc *ov[2] = { rs->paf, rs->reverse_paf };
	FILE *fp = stdout;

	for (j = 0; j < rs->total_reads; ++j) {
		if (max_len < (int)Get_NAME_LENGTH(*rs, j))
			max_len = Get_NAME_LENGTH(*rs, j);
	}
	GFA_MALLOC(s, max_len + 1);

	h = ss_init();
	for (i = 0; i < n; ++i) {
		k = ss_put(h, gfa_strdup(list[i]), &absent);
		kh_val(h, k) = -1;
	}

	for (m = 0; m < n_rounds; ++m) {
		for (j = 0; j < rs->total_reads; ++j) {
			for (l = 0; l < 2; ++l) {
				const ma_hit_t_alloc *o = &ov[l][j];
				for (i = 0; i < (int)o->length; ++i) {
					uint64_t q = Get_qn(o->buffer[i]);
					uint64_t t = Get_tn(o->buffer[i]);
					int q_hit = 0, t_hit = 0;
					strncpy(s, Get_NAME(*rs, q), Get_NAME_LENGTH(*rs, q)); s[Get_NAME_LENGTH(*rs, q)] = 0;
					k = ss_get(h, s);
					q_hit = (k < kh_end(h) && kh_val(h, k) < m);
					strncpy(s, Get_NAME(*rs, t), Get_NAME_LENGTH(*rs, t)); s[Get_NAME_LENGTH(*rs, t)] = 0;
					k = ss_get(h, s);
					t_hit = (k < kh_end(h) && kh_val(h, k) < m);
					if ((!q_hit && !t_hit) || (q_hit && t_hit)) continue;
					if (!q_hit) {
						char *tmp = gfa_strndup(Get_NAME(*rs, q), Get_NAME_LENGTH(*rs, q));
						k = ss_put(h, tmp, &absent);
						if (absent) kh_val(h, k) = m;
						else free(tmp);
					}
					if (!t_hit) {
						char *tmp = gfa_strndup(Get_NAME(*rs, t), Get_NAME_LENGTH(*rs, t));
						k = ss_put(h, tmp, &absent);
						if (absent) kh_val(h, k) = m;
						else free(tmp);
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
				strncpy(s, Get_NAME(*rs, q), Get_NAME_LENGTH(*rs, q)); s[Get_NAME_LENGTH(*rs, q)] = 0;
				k = ss_get(h, s);
				q_hit = (k < kh_end(h));
				strncpy(s, Get_NAME(*rs, t), Get_NAME_LENGTH(*rs, t)); s[Get_NAME_LENGTH(*rs, t)] = 0;
				k = ss_get(h, s);
				t_hit = (k < kh_end(h));
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

	for (k = 0; k != kh_end(h); ++k)
		if (kh_exist(h, k))
			free((char*)kh_key(h, k));
	ss_destroy(h);
	free(s);
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
