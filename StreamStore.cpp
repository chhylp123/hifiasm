#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#define _FILE_OFFSET_BITS 64

#include "StreamStore.h"
#include "kthread.h"

#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#define HA_STB_MAGIC       0x53544231U
#define HA_QUAL_SLICE_SIZE (1U << 20)
#define HA_STB_F_ALL_QUAL  0x1U
#define HA_STB_F_ANY_QUAL  0x2U

typedef struct {
    uint32_t read_idx;
    uint32_t pos;
    uint8_t ch;
} ha_seq_exc_t;

typedef struct {
    uint8_t *ram_blob;
    uint64_t blob_size;
    uint64_t disk_offset;
    uint8_t on_disk;
} ha_stream_batch_meta_t;

typedef struct {
    ha_stream_batch_meta_t *batches;
    size_t n_batches;
    size_t m_batches;
    size_t read_cursor;
    int captured;
} ha_stream_file_slot_t;

struct ha_stream_store_s {
    int num_files;
    ha_stream_file_slot_t *files;
    uint64_t max_ram_bytes;
    uint64_t ram_used_bytes;
    uint64_t peak_ram_bytes;
    uint64_t disk_used_bytes;
    uint64_t next_disk_offset;
    uint64_t total_records;
    uint64_t total_bases;
    uint64_t total_comp_bytes;
    int spill_fd;
    char *spill_dir;
    pthread_mutex_t mu;
};

static uint8_t g_dna_decode_init = 0;
static char g_dna_decode_tab[256][4];

static void ensure_dna_decode_tab(void)
{
    if (__builtin_expect(g_dna_decode_init, 1)) return;
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    for (int b = 0; b < 256; ++b) {
        g_dna_decode_tab[b][0] = bases[(b >> 6) & 3];
        g_dna_decode_tab[b][1] = bases[(b >> 4) & 3];
        g_dna_decode_tab[b][2] = bases[(b >> 2) & 3];
        g_dna_decode_tab[b][3] = bases[b & 3];
    }
    g_dna_decode_init = 1;
}

static inline void put_u32(uint8_t *p, uint32_t v) { memcpy(p, &v, 4); }
static inline void put_u64(uint8_t *p, uint64_t v) { memcpy(p, &v, 8); }
static inline uint32_t get_u32(const uint8_t *p) { uint32_t v; memcpy(&v, p, 4); return v; }
static inline uint64_t get_u64(const uint8_t *p) { uint64_t v; memcpy(&v, p, 8); return v; }

static int parse_dev_fd(const char *fn)
{
    const char *num_str = NULL;
    if (strncmp(fn, "/dev/fd/", 8) == 0) num_str = fn + 8;
    else if (strncmp(fn, "/proc/self/fd/", 14) == 0) num_str = fn + 14;
    if (!num_str || !num_str[0]) return -1;
    char *endp = NULL;
    long v = strtol(num_str, &endp, 10);
    if (!endp || *endp != '\0' || v < 0 || v > INT_MAX) return -1;
    return (int)v;
}

int ha_is_stream_path(const char *fn)
{
    if (!fn || !fn[0]) return 0;
    if (strcmp(fn, "-") == 0 || strcmp(fn, "/dev/stdin") == 0) return 1;
    if (strncmp(fn, "/dev/fd/", 8) == 0 || strncmp(fn, "/proc/", 6) == 0) return 1;
    struct stat st;
    if (stat(fn, &st) == 0) {
        if (!S_ISREG(st.st_mode)) return 1;
    }
    return 0;
}

int ha_check_input_readable(const char *fn, const char *opt_label)
{
    if (!fn || !fn[0]) {
        fprintf(stderr, "[ERROR] file does not exist (-%s)\n", opt_label ? opt_label : "input");
        return 0;
    }
    if (strcmp(fn, "-") == 0 || strcmp(fn, "/dev/stdin") == 0) {
        if (fcntl(STDIN_FILENO, F_GETFD) == -1) {
            fprintf(stderr, "[ERROR] stdin is not open (-%s)\n", opt_label ? opt_label : "stdin");
            return 0;
        }
        return 1;
    }
    int fd_num = parse_dev_fd(fn);
    if (fd_num >= 0) {
        if (fcntl(fd_num, F_GETFD) == -1) {
            fprintf(stderr, "[ERROR] %s is not an open descriptor (-%s)\n", fn, opt_label ? opt_label : "fd");
            return 0;
        }
        return 1;
    }
    struct stat st;
    if (stat(fn, &st) != 0) {
        fprintf(stderr, "[ERROR] %s does not exist (-%s)\n", fn, opt_label ? opt_label : "input");
        return 0;
    }
    if (!S_ISREG(st.st_mode)) {
        if (access(fn, R_OK) != 0) {
            fprintf(stderr, "[ERROR] %s is not readable (-%s)\n", fn, opt_label ? opt_label : "input");
            return 0;
        }
        return 1;
    }
    FILE *fp = fopen(fn, "r");
    if (!fp) {
        fprintf(stderr, "[ERROR] %s does not exist (-%s)\n", fn, opt_label ? opt_label : "input");
        return 0;
    }
    fclose(fp);
    return 1;
}

gzFile ha_gzopen_read(const char *fn)
{
    if (!fn || !fn[0]) return 0;
    gzFile fp = 0;
    if (strcmp(fn, "-") == 0 || strcmp(fn, "/dev/stdin") == 0) {
        int dup_fd = dup(STDIN_FILENO);
        if (dup_fd < 0) return 0;
        fp = gzdopen(dup_fd, "r");
        if (!fp) {
            close(dup_fd);
            return 0;
        }
    } else {
        int fd_num = parse_dev_fd(fn);
        if (fd_num >= 0) {
            int dup_fd = dup(fd_num);
            if (dup_fd >= 0) {
                fp = gzdopen(dup_fd, "r");
                if (!fp) {
                    close(dup_fd);
                    return 0;
                }
            }
        } else {
            fp = gzopen(fn, "r");
        }
    }
#if ZLIB_VERNUM >= 0x1240
    if (fp) gzbuffer(fp, 131072);
#endif
    return fp;
}

static char *resolve_spill_dir(const char *fallback_path)
{
    const char *env_tmp = getenv("TMPDIR");
    if (env_tmp && env_tmp[0]) {
        return strdup(env_tmp);
    }
    if (fallback_path && fallback_path[0]) {
        const char *slash = strrchr(fallback_path, '/');
        if (slash && slash > fallback_path) {
            size_t dlen = (size_t)(slash - fallback_path);
            char *d = (char *)malloc(dlen + 1);
            memcpy(d, fallback_path, dlen);
            d[dlen] = '\0';
            return d;
        }
    }
    return strdup(".");
}

static int open_unlinked_spill_fd(const char *dir)
{
    int fd = -1;
#ifdef O_TMPFILE
    fd = open(dir, O_TMPFILE | O_RDWR | O_CLOEXEC, 0600);
    if (fd >= 0) return fd;
#endif
    size_t dlen = strlen(dir);
    const char *suffix = "/.hifiasm_stream_XXXXXX";
    char *tmpl = (char *)malloc(dlen + strlen(suffix) + 1);
    if (!tmpl) return -1;
    memcpy(tmpl, dir, dlen);
    strcpy(tmpl + dlen, suffix);
    fd = mkostemp(tmpl, O_CLOEXEC);
    if (fd >= 0) {
        unlink(tmpl);
    }
    free(tmpl);
    return fd;
}

static void pwrite_all(int fd, const uint8_t *buf, uint64_t len, uint64_t offset)
{
    uint64_t done = 0;
    while (done < len) {
        size_t chunk = (len - done > (1ULL << 30)) ? (size_t)(1ULL << 30) : (size_t)(len - done);
        ssize_t nw = pwrite(fd, buf + done, chunk, (off_t)(offset + done));
        if (nw < 0) {
            if (errno == EINTR) continue;
            fprintf(stderr, "[ERROR] StreamStore pwrite failed: %s\n", strerror(errno));
            exit(1);
        }
        if (nw == 0) {
            fprintf(stderr, "[ERROR] StreamStore short pwrite (disk full in TMPDIR)\n");
            exit(1);
        }
        done += (uint64_t)nw;
    }
}

static void pread_all(int fd, uint8_t *buf, uint64_t len, uint64_t offset)
{
    uint64_t done = 0;
    while (done < len) {
        size_t chunk = (len - done > (1ULL << 30)) ? (size_t)(1ULL << 30) : (size_t)(len - done);
        ssize_t nr = pread(fd, buf + done, chunk, (off_t)(offset + done));
        if (nr < 0) {
            if (errno == EINTR) continue;
            fprintf(stderr, "[ERROR] StreamStore pread failed: %s\n", strerror(errno));
            exit(1);
        }
        if (nr == 0) {
            fprintf(stderr, "[ERROR] StreamStore unexpected EOF on pread\n");
            exit(1);
        }
        done += (uint64_t)nr;
    }
}

static uint8_t *zlib_compress_buf(const uint8_t *src, uint64_t src_len, uint64_t *dst_len_out, int level)
{
    if (src_len == 0) {
        *dst_len_out = 0;
        return NULL;
    }
    uLongf bound = compressBound((uLong)src_len);
    uint8_t *dst = (uint8_t *)malloc((size_t)bound);
    if (!dst) {
        fprintf(stderr, "[ERROR] StreamStore malloc failed during compression\n");
        exit(1);
    }
    uLongf dest_len = bound;
    int rc = compress2(dst, &dest_len, src, (uLong)src_len, level);
    if (rc != Z_OK) {
        fprintf(stderr, "[ERROR] StreamStore zlib compress2 failed (%d)\n", rc);
        exit(1);
    }
    *dst_len_out = (uint64_t)dest_len;
    return dst;
}

static void zlib_decompress_buf(const uint8_t *src, uint64_t src_len, uint8_t *dst, uint64_t dst_len)
{
    if (dst_len == 0) return;
    uLongf out_len = (uLongf)dst_len;
    int rc = uncompress(dst, &out_len, src, (uLong)src_len);
    if (rc != Z_OK || (uint64_t)out_len != dst_len) {
        fprintf(stderr, "[ERROR] StreamStore zlib uncompress failed (%d)\n", rc);
        exit(1);
    }
}

typedef struct {
    const uint8_t *raw_qual;
    uint64_t total_qual_len;
    uint8_t **slice_bufs;
    uint32_t *slice_comp_lens;
    uint32_t *slice_raw_lens;
} qual_comp_ctx_t;

static void worker_qual_compress(void *data, long idx, int tid)
{
    (void)tid;
    qual_comp_ctx_t *ctx = (qual_comp_ctx_t *)data;
    uint64_t start = (uint64_t)idx * HA_QUAL_SLICE_SIZE;
    uint64_t rem = ctx->total_qual_len - start;
    uint32_t slen = (rem > HA_QUAL_SLICE_SIZE) ? HA_QUAL_SLICE_SIZE : (uint32_t)rem;
    uint64_t clen = 0;
    ctx->slice_raw_lens[idx] = slen;
    ctx->slice_bufs[idx] = zlib_compress_buf(ctx->raw_qual + start, slen, &clen, 1);
    ctx->slice_comp_lens[idx] = (uint32_t)clen;
}

typedef struct {
    const uint8_t *comp_base;
    const uint64_t *slice_comp_offs;
    const uint32_t *slice_comp_lens;
    const uint32_t *slice_raw_lens;
    uint8_t *dst_qual;
} qual_decomp_ctx_t;

static void worker_qual_decompress(void *data, long idx, int tid)
{
    (void)tid;
    qual_decomp_ctx_t *ctx = (qual_decomp_ctx_t *)data;
    uint64_t dst_off = (uint64_t)idx * HA_QUAL_SLICE_SIZE;
    zlib_decompress_buf(
        ctx->comp_base + ctx->slice_comp_offs[idx],
        ctx->slice_comp_lens[idx],
        ctx->dst_qual + dst_off,
        ctx->slice_raw_lens[idx]
    );
}

ha_stream_store_t *ha_stream_store_create(int num_files, double max_ram_gb, const char *fallback_path)
{
    ensure_dna_decode_tab();
    ha_stream_store_t *st = (ha_stream_store_t *)calloc(1, sizeof(ha_stream_store_t));
    st->num_files = (num_files > 0) ? num_files : 1;
    st->files = (ha_stream_file_slot_t *)calloc((size_t)st->num_files, sizeof(ha_stream_file_slot_t));
    if (max_ram_gb <= 0.0) {
        st->max_ram_bytes = 0;
    } else {
        double bytes = max_ram_gb * 1073741824.0;
        st->max_ram_bytes = (bytes >= (double)UINT64_MAX) ? UINT64_MAX : (uint64_t)bytes;
    }
    st->spill_fd = -1;
    st->spill_dir = resolve_spill_dir(fallback_path);
    pthread_mutex_init(&st->mu, NULL);
    return st;
}

void ha_stream_store_destroy(ha_stream_store_t **store_ptr)
{
    if (!store_ptr || !(*store_ptr)) return;
    ha_stream_store_t *st = *store_ptr;
    for (int f = 0; f < st->num_files; ++f) {
        ha_stream_file_slot_t *slot = &st->files[f];
        for (size_t b = 0; b < slot->n_batches; ++b) {
            free(slot->batches[b].ram_blob);
        }
        free(slot->batches);
    }
    free(st->files);
    if (st->spill_fd >= 0) close(st->spill_fd);
    free(st->spill_dir);
    pthread_mutex_destroy(&st->mu);
    free(st);
    *store_ptr = NULL;
}

int ha_stream_store_is_captured(const ha_stream_store_t *store, int file_idx)
{
    if (!store || file_idx < 0 || file_idx >= store->num_files) return 0;
    return store->files[file_idx].captured;
}

void ha_stream_store_mark_captured(ha_stream_store_t *store, int file_idx)
{
    if (!store || file_idx < 0 || file_idx >= store->num_files) return;
    pthread_mutex_lock(&store->mu);
    store->files[file_idx].captured = 1;
    store->files[file_idx].read_cursor = 0;
    pthread_mutex_unlock(&store->mu);
}

void ha_stream_store_rewind_file(ha_stream_store_t *store, int file_idx)
{
    if (!store || file_idx < 0 || file_idx >= store->num_files) return;
    pthread_mutex_lock(&store->mu);
    store->files[file_idx].read_cursor = 0;
    pthread_mutex_unlock(&store->mu);
}

static inline int encode_dna_char(uint8_t ch, uint8_t *two_bit, int *is_n, int *is_lower, int *is_exc)
{
    *is_n = 0;
    *is_lower = 0;
    *is_exc = 0;
    switch (ch) {
        case 'A': *two_bit = 0; return 0;
        case 'C': *two_bit = 1; return 0;
        case 'G': *two_bit = 2; return 0;
        case 'T': *two_bit = 3; return 0;
        case 'N': *two_bit = 0; *is_n = 1; return 0;
        case 'a': *two_bit = 0; *is_lower = 1; return 0;
        case 'c': *two_bit = 1; *is_lower = 1; return 0;
        case 'g': *two_bit = 2; *is_lower = 1; return 0;
        case 't': *two_bit = 3; *is_lower = 1; return 0;
        case 'n': *two_bit = 0; *is_n = 1; *is_lower = 1; return 0;
        default:  *two_bit = 0; *is_exc = 1; return 0;
    }
}

/*
 * Serialized StreamBatch Header (96 bytes):
 *  [0..4)   magic (uint32)
 *  [4..8)   n_seq (uint32)
 *  [8..12)  flags (uint32)
 *  [12..16) n_qual_slices (uint32)
 *  [16..24) sum_seq_len (uint64)
 *  [24..32) sum_dna2_bytes (uint64)
 *  [32..40) sum_name_len (uint64)
 *  [40..48) comp_name_bytes (uint64)
 *  [48..56) sum_qual_len (uint64)
 *  [56..64) comp_qual_bytes (uint64)
 *  [64..72) total_n_pos (uint64)
 *  [72..80) comp_n_bytes (uint64)
 *  [80..88) comp_case_bytes (uint64)
 *  [88..96) comp_exc_bytes (uint64)
 *  [96..104) total_exc (uint64)
 */
#define HA_STB_HDR_BYTES 104U

void ha_stream_store_append_batch(
    ha_stream_store_t *store,
    int file_idx,
    int n_seq_in,
    char **seqs,
    const int *seq_lens,
    const ha_stream_raw_rec_t *raw_recs,
    int n_threads
)
{
    if (!store || n_seq_in <= 0) return;
    ensure_dna_decode_tab();

    uint32_t n_seq = (uint32_t)n_seq_in;
    uint64_t sum_seq_len = 0;
    uint64_t sum_dna2_bytes = 0;
    uint64_t sum_name_len = 0;
    uint64_t sum_qual_len = 0;
    uint32_t n_with_qual = 0;

    for (uint32_t i = 0; i < n_seq; ++i) {
        uint32_t sl = (uint32_t)seq_lens[i];
        sum_seq_len += sl;
        sum_dna2_bytes += (sl + 3U) >> 2;
        sum_name_len += raw_recs[i].name_len;
        if (raw_recs[i].has_qual) {
            ++n_with_qual;
            sum_qual_len += sl;
        }
    }

    uint32_t flags = 0;
    if (n_with_qual == n_seq) flags |= HA_STB_F_ALL_QUAL | HA_STB_F_ANY_QUAL;
    else if (n_with_qual > 0) flags |= HA_STB_F_ANY_QUAL;

    uint8_t *dna2 = (sum_dna2_bytes > 0) ? (uint8_t *)malloc((size_t)sum_dna2_bytes) : NULL;
    uint32_t *n_counts = (uint32_t *)calloc(n_seq, sizeof(uint32_t));
    size_t n_pos_cap = 256;
    uint64_t total_n_pos = 0;
    uint32_t *n_positions = (uint32_t *)malloc(n_pos_cap * sizeof(uint32_t));

    uint64_t case_bm_bytes = (sum_seq_len + 7ULL) >> 3;
    uint8_t *case_bm = NULL;
    int any_lower = 0;

    size_t exc_cap = 0;
    uint64_t total_exc = 0;
    ha_seq_exc_t *excs = NULL;

    uint64_t dna2_off = 0;
    uint64_t global_base_idx = 0;

    for (uint32_t r = 0; r < n_seq; ++r) {
        const uint8_t *s = (const uint8_t *)seqs[r];
        uint32_t sl = (uint32_t)seq_lens[r];
        uint32_t p = 0;
        while (p + 4U <= sl) {
            uint8_t b0, b1, b2, b3;
            int n0, n1, n2, n3, lw0, lw1, lw2, lw3, ex0, ex1, ex2, ex3;
            encode_dna_char(s[p],     &b0, &n0, &lw0, &ex0);
            encode_dna_char(s[p + 1], &b1, &n1, &lw1, &ex1);
            encode_dna_char(s[p + 2], &b2, &n2, &lw2, &ex2);
            encode_dna_char(s[p + 3], &b3, &n3, &lw3, &ex3);
            dna2[dna2_off++] = (uint8_t)((b0 << 6) | (b1 << 4) | (b2 << 2) | b3);

            if (__builtin_expect(n0 | n1 | n2 | n3 | lw0 | lw1 | lw2 | lw3 | ex0 | ex1 | ex2 | ex3, 0)) {
                const int ns[4] = {n0, n1, n2, n3};
                const int lws[4] = {lw0, lw1, lw2, lw3};
                const int exs[4] = {ex0, ex1, ex2, ex3};
                for (uint32_t k = 0; k < 4; ++k) {
                    if (ns[k]) {
                        if (total_n_pos == n_pos_cap) {
                            n_pos_cap <<= 1;
                            n_positions = (uint32_t *)realloc(n_positions, n_pos_cap * sizeof(uint32_t));
                        }
                        n_positions[total_n_pos++] = p + k;
                        n_counts[r]++;
                    }
                    if (lws[k]) {
                        if (!any_lower) {
                            case_bm = (uint8_t *)calloc((size_t)case_bm_bytes, 1);
                            any_lower = 1;
                        }
                        uint64_t gbi = global_base_idx + p + k;
                        case_bm[gbi >> 3] |= (uint8_t)(1U << (gbi & 7U));
                    }
                    if (exs[k]) {
                        if (total_exc == exc_cap) {
                            exc_cap = exc_cap ? (exc_cap << 1) : 64;
                            excs = (ha_seq_exc_t *)realloc(excs, exc_cap * sizeof(ha_seq_exc_t));
                        }
                        excs[total_exc].read_idx = r;
                        excs[total_exc].pos = p + k;
                        excs[total_exc].ch = s[p + k];
                        total_exc++;
                    }
                }
            }
            p += 4U;
        }
        if (p < sl) {
            uint8_t packed = 0;
            uint32_t shift = 6;
            for (; p < sl; ++p, shift -= 2) {
                uint8_t b0;
                int n0, lw0, ex0;
                encode_dna_char(s[p], &b0, &n0, &lw0, &ex0);
                packed |= (uint8_t)(b0 << shift);
                if (n0) {
                    if (total_n_pos == n_pos_cap) {
                        n_pos_cap <<= 1;
                        n_positions = (uint32_t *)realloc(n_positions, n_pos_cap * sizeof(uint32_t));
                    }
                    n_positions[total_n_pos++] = p;
                    n_counts[r]++;
                }
                if (lw0) {
                    if (!any_lower) {
                        case_bm = (uint8_t *)calloc((size_t)case_bm_bytes, 1);
                        any_lower = 1;
                    }
                    uint64_t gbi = global_base_idx + p;
                    case_bm[gbi >> 3] |= (uint8_t)(1U << (gbi & 7U));
                }
                if (ex0) {
                    if (total_exc == exc_cap) {
                        exc_cap = exc_cap ? (exc_cap << 1) : 64;
                        excs = (ha_seq_exc_t *)realloc(excs, exc_cap * sizeof(ha_seq_exc_t));
                    }
                    excs[total_exc].read_idx = r;
                    excs[total_exc].pos = p;
                    excs[total_exc].ch = s[p];
                    total_exc++;
                }
            }
            dna2[dna2_off++] = packed;
        }
        global_base_idx += sl;
    }
    assert(dna2_off == sum_dna2_bytes);

    /* Compress concatenated read names */
    uint8_t *raw_names = (sum_name_len > 0) ? (uint8_t *)malloc((size_t)sum_name_len) : NULL;
    uint64_t noff = 0;
    for (uint32_t r = 0; r < n_seq; ++r) {
        if (raw_recs[r].name_len > 0) {
            memcpy(raw_names + noff, raw_recs[r].name, raw_recs[r].name_len);
            noff += raw_recs[r].name_len;
        }
    }
    uint64_t comp_name_bytes = 0;
    uint8_t *comp_names = zlib_compress_buf(raw_names, sum_name_len, &comp_name_bytes, 1);
    free(raw_names);

    /* Compress exact N positions if present */
    uint64_t comp_n_bytes = 0;
    uint8_t *comp_n = NULL;
    if (total_n_pos > 0) {
        uint64_t raw_n_bytes = ((uint64_t)n_seq + total_n_pos) * sizeof(uint32_t);
        uint8_t *raw_n_buf = (uint8_t *)malloc((size_t)raw_n_bytes);
        memcpy(raw_n_buf, n_counts, (size_t)n_seq * sizeof(uint32_t));
        memcpy(raw_n_buf + (size_t)n_seq * sizeof(uint32_t), n_positions, (size_t)total_n_pos * sizeof(uint32_t));
        comp_n = zlib_compress_buf(raw_n_buf, raw_n_bytes, &comp_n_bytes, 1);
        free(raw_n_buf);
    }
    free(n_counts);
    free(n_positions);

    /* Compress lowercase bitmap if present */
    uint64_t comp_case_bytes = 0;
    uint8_t *comp_case = NULL;
    if (any_lower && case_bm) {
        comp_case = zlib_compress_buf(case_bm, case_bm_bytes, &comp_case_bytes, 1);
    }
    free(case_bm);

    /* Compress IUPAC exceptions if present */
    uint64_t comp_exc_bytes = 0;
    uint8_t *comp_exc = NULL;
    if (total_exc > 0) {
        uint64_t raw_exc_bytes = total_exc * 9ULL;
        uint8_t *raw_exc_buf = (uint8_t *)malloc((size_t)raw_exc_bytes);
        for (uint64_t e = 0; e < total_exc; ++e) {
            put_u32(raw_exc_buf + e * 9ULL, excs[e].read_idx);
            put_u32(raw_exc_buf + e * 9ULL + 4ULL, excs[e].pos);
            raw_exc_buf[e * 9ULL + 8ULL] = excs[e].ch;
        }
        comp_exc = zlib_compress_buf(raw_exc_buf, raw_exc_bytes, &comp_exc_bytes, 1);
        free(raw_exc_buf);
    }
    free(excs);

    /* Multi-threaded block compression of concatenated base qualities */
    uint32_t n_qual_slices = 0;
    uint64_t comp_qual_bytes = 0;
    uint8_t **qual_slice_bufs = NULL;
    uint32_t *qual_slice_clens = NULL;
    uint32_t *qual_slice_rlens = NULL;

    if (sum_qual_len > 0) {
        uint8_t *raw_qual = (uint8_t *)malloc((size_t)sum_qual_len);
        uint64_t qoff = 0;
        for (uint32_t r = 0; r < n_seq; ++r) {
            if (raw_recs[r].has_qual && seq_lens[r] > 0) {
                memcpy(raw_qual + qoff, raw_recs[r].qual, (size_t)seq_lens[r]);
                qoff += (uint64_t)seq_lens[r];
            }
        }
        assert(qoff == sum_qual_len);

        n_qual_slices = (uint32_t)((sum_qual_len + HA_QUAL_SLICE_SIZE - 1ULL) / HA_QUAL_SLICE_SIZE);
        qual_slice_bufs = (uint8_t **)calloc(n_qual_slices, sizeof(uint8_t *));
        qual_slice_clens = (uint32_t *)calloc(n_qual_slices, sizeof(uint32_t));
        qual_slice_rlens = (uint32_t *)calloc(n_qual_slices, sizeof(uint32_t));

        qual_comp_ctx_t qctx;
        qctx.raw_qual = raw_qual;
        qctx.total_qual_len = sum_qual_len;
        qctx.slice_bufs = qual_slice_bufs;
        qctx.slice_comp_lens = qual_slice_clens;
        qctx.slice_raw_lens = qual_slice_rlens;

        int nth = (n_threads > 1) ? n_threads : 1;
        if ((uint32_t)nth > n_qual_slices) nth = (int)n_qual_slices;
        kt_for(nth, worker_qual_compress, &qctx, (long)n_qual_slices);
        free(raw_qual);

        for (uint32_t s = 0; s < n_qual_slices; ++s) {
            comp_qual_bytes += qual_slice_clens[s];
        }
    }

    /* Assemble single contiguous serialized batch blob */
    uint64_t qual_Bm_bytes = ((flags & HA_STB_F_ANY_QUAL) && !(flags & HA_STB_F_ALL_QUAL))
                                 ? ((n_seq + 7ULL) >> 3)
                                 : 0ULL;
    uint64_t slice_tbl_bytes = (uint64_t)n_qual_slices * 8ULL;
    uint64_t blob_size = HA_STB_HDR_BYTES
                       + (uint64_t)n_seq * 8ULL
                       + qual_Bm_bytes
                       + sum_dna2_bytes
                       + comp_name_bytes
                       + comp_n_bytes
                       + comp_case_bytes
                       + comp_exc_bytes
                       + slice_tbl_bytes
                       + comp_qual_bytes;

    uint8_t *blob = (uint8_t *)malloc((size_t)blob_size);
    if (!blob) {
        fprintf(stderr, "[ERROR] StreamStore failed to allocate batch blob (%lu bytes)\n", (unsigned long)blob_size);
        exit(1);
    }

    put_u32(blob + 0, HA_STB_MAGIC);
    put_u32(blob + 4, n_seq);
    put_u32(blob + 8, flags);
    put_u32(blob + 12, n_qual_slices);
    put_u64(blob + 16, sum_seq_len);
    put_u64(blob + 24, sum_dna2_bytes);
    put_u64(blob + 32, sum_name_len);
    put_u64(blob + 40, comp_name_bytes);
    put_u64(blob + 48, sum_qual_len);
    put_u64(blob + 56, comp_qual_bytes);
    put_u64(blob + 64, total_n_pos);
    put_u64(blob + 72, comp_n_bytes);
    put_u64(blob + 80, comp_case_bytes);
    put_u64(blob + 88, comp_exc_bytes);
    put_u64(blob + 96, total_exc);

    uint64_t cur = HA_STB_HDR_BYTES;
    for (uint32_t r = 0; r < n_seq; ++r) {
        put_u32(blob + cur + (uint64_t)r * 4ULL, (uint32_t)seq_lens[r]);
    }
    cur += (uint64_t)n_seq * 4ULL;

    for (uint32_t r = 0; r < n_seq; ++r) {
        put_u32(blob + cur + (uint64_t)r * 4ULL, raw_recs[r].name_len);
    }
    cur += (uint64_t)n_seq * 4ULL;

    if (qual_Bm_bytes > 0) {
        memset(blob + cur, 0, (size_t)qual_Bm_bytes);
        for (uint32_t r = 0; r < n_seq; ++r) {
            if (raw_recs[r].has_qual) {
                blob[cur + (r >> 3)] |= (uint8_t)(1U << (r & 7U));
            }
        }
        cur += qual_Bm_bytes;
    }

    if (sum_dna2_bytes > 0) {
        memcpy(blob + cur, dna2, (size_t)sum_dna2_bytes);
        cur += sum_dna2_bytes;
    }
    free(dna2);

    if (comp_name_bytes > 0) {
        memcpy(blob + cur, comp_names, (size_t)comp_name_bytes);
        cur += comp_name_bytes;
    }
    free(comp_names);

    if (comp_n_bytes > 0) {
        memcpy(blob + cur, comp_n, (size_t)comp_n_bytes);
        cur += comp_n_bytes;
    }
    free(comp_n);

    if (comp_case_bytes > 0) {
        memcpy(blob + cur, comp_case, (size_t)comp_case_bytes);
        cur += comp_case_bytes;
    }
    free(comp_case);

    if (comp_exc_bytes > 0) {
        memcpy(blob + cur, comp_exc, (size_t)comp_exc_bytes);
        cur += comp_exc_bytes;
    }
    free(comp_exc);

    if (n_qual_slices > 0) {
        for (uint32_t s = 0; s < n_qual_slices; ++s) {
            put_u32(blob + cur + (uint64_t)s * 8ULL, qual_slice_clens[s]);
            put_u32(blob + cur + (uint64_t)s * 8ULL + 4ULL, qual_slice_rlens[s]);
        }
        cur += slice_tbl_bytes;
        for (uint32_t s = 0; s < n_qual_slices; ++s) {
            if (qual_slice_clens[s] > 0) {
                memcpy(blob + cur, qual_slice_bufs[s], qual_slice_clens[s]);
                cur += qual_slice_clens[s];
            }
            free(qual_slice_bufs[s]);
        }
        free(qual_slice_bufs);
        free(qual_slice_clens);
        free(qual_slice_rlens);
    }
    assert(cur == blob_size);

    /* Commit batch to RAM or spill to unlinked O_TMPFILE in TMPDIR */
    pthread_mutex_lock(&store->mu);
    ha_stream_file_slot_t *slot = &store->files[file_idx];
    if (slot->n_batches == slot->m_batches) {
        slot->m_batches = slot->m_batches ? (slot->m_batches << 1) : 16;
        slot->batches = (ha_stream_batch_meta_t *)realloc(
            slot->batches, slot->m_batches * sizeof(ha_stream_batch_meta_t)
        );
    }
    ha_stream_batch_meta_t *meta = &slot->batches[slot->n_batches++];
    meta->blob_size = blob_size;

    if (store->ram_used_bytes + blob_size <= store->max_ram_bytes) {
        meta->ram_blob = blob;
        meta->disk_offset = 0;
        meta->on_disk = 0;
        store->ram_used_bytes += blob_size;
        if (store->ram_used_bytes > store->peak_ram_bytes) {
            store->peak_ram_bytes = store->ram_used_bytes;
        }
    } else {
        if (store->spill_fd < 0) {
            store->spill_fd = open_unlinked_spill_fd(store->spill_dir);
            if (store->spill_fd < 0) {
                fprintf(stderr, "[ERROR] StreamStore failed to open unlinked temp file in %s: %s\n",
                        store->spill_dir, strerror(errno));
                exit(1);
            }
        }
        uint64_t off = store->next_disk_offset;
        store->next_disk_offset += blob_size;
        store->disk_used_bytes += blob_size;
        meta->ram_blob = NULL;
        meta->disk_offset = off;
        meta->on_disk = 1;
        pthread_mutex_unlock(&store->mu);
        pwrite_all(store->spill_fd, blob, blob_size, off);
        free(blob);
        pthread_mutex_lock(&store->mu);
    }

    store->total_records += n_seq;
    store->total_bases += sum_seq_len;
    store->total_comp_bytes += blob_size;
    pthread_mutex_unlock(&store->mu);
}

int ha_stream_store_replay_batch(
    ha_stream_store_t *store,
    int file_idx,
    int consume_ram,
    int n_threads,
    ha_stream_decoded_batch_t *out
)
{
    memset(out, 0, sizeof(*out));
    if (!store || file_idx < 0 || file_idx >= store->num_files) return 0;
    ensure_dna_decode_tab();

    uint8_t *blob = NULL;
    uint64_t blob_size = 0;
    int free_blob = 0;

    pthread_mutex_lock(&store->mu);
    ha_stream_file_slot_t *slot = &store->files[file_idx];
    if (slot->read_cursor >= slot->n_batches) {
        pthread_mutex_unlock(&store->mu);
        return 0;
    }
    ha_stream_batch_meta_t *meta = &slot->batches[slot->read_cursor++];
    blob_size = meta->blob_size;
    if (!meta->on_disk) {
        blob = meta->ram_blob;
        if (consume_ram) {
            meta->ram_blob = NULL;
            if (store->ram_used_bytes >= blob_size) store->ram_used_bytes -= blob_size;
            else store->ram_used_bytes = 0;
            free_blob = 1;
        }
        pthread_mutex_unlock(&store->mu);
    } else {
        uint64_t off = meta->disk_offset;
        int fd = store->spill_fd;
        pthread_mutex_unlock(&store->mu);
        blob = (uint8_t *)malloc((size_t)blob_size);
        if (!blob) {
            fprintf(stderr, "[ERROR] StreamStore failed to allocate replay buffer\n");
            exit(1);
        }
        pread_all(fd, blob, blob_size, off);
        free_blob = 1;
    }

    if (!blob || get_u32(blob) != HA_STB_MAGIC) {
        fprintf(stderr, "[ERROR] StreamStore corrupt batch header\n");
        exit(1);
    }

    uint32_t n_seq = get_u32(blob + 4);
    uint32_t flags = get_u32(blob + 8);
    uint32_t n_qual_slices = get_u32(blob + 12);
    uint64_t sum_seq_len = get_u64(blob + 16);
    uint64_t sum_dna2_bytes = get_u64(blob + 24);
    uint64_t sum_name_len = get_u64(blob + 32);
    uint64_t comp_name_bytes = get_u64(blob + 40);
    uint64_t sum_qual_len = get_u64(blob + 48);
    uint64_t comp_qual_bytes = get_u64(blob + 56);
    uint64_t total_n_pos = get_u64(blob + 64);
    uint64_t comp_n_bytes = get_u64(blob + 72);
    uint64_t comp_case_bytes = get_u64(blob + 80);
    uint64_t comp_exc_bytes = get_u64(blob + 88);
    uint64_t total_exc = get_u64(blob + 96);

    out->n_seq = n_seq;
    out->seq_lens = (int *)malloc((size_t)n_seq * sizeof(int));
    out->seqs = (char **)malloc((size_t)n_seq * sizeof(char *));
    out->name_lens = (uint32_t *)malloc((size_t)n_seq * sizeof(uint32_t));
    out->name_offs = (uint32_t *)malloc((size_t)n_seq * sizeof(uint32_t));
    out->has_qual = (uint8_t *)calloc((size_t)n_seq, 1);
    out->qual_offs = (uint32_t *)calloc((size_t)n_seq, sizeof(uint32_t));

    uint64_t cur = HA_STB_HDR_BYTES;
    for (uint32_t r = 0; r < n_seq; ++r) {
        out->seq_lens[r] = (int)get_u32(blob + cur + (uint64_t)r * 4ULL);
    }
    cur += (uint64_t)n_seq * 4ULL;

    uint32_t running_noff = 0;
    for (uint32_t r = 0; r < n_seq; ++r) {
        uint32_t nl = get_u32(blob + cur + (uint64_t)r * 4ULL);
        out->name_lens[r] = nl;
        out->name_offs[r] = running_noff;
        running_noff += nl + 1U;
    }
    cur += (uint64_t)n_seq * 4ULL;

    if (flags & HA_STB_F_ALL_QUAL) {
        memset(out->has_qual, 1, (size_t)n_seq);
    } else if (flags & HA_STB_F_ANY_QUAL) {
        uint64_t qual_bm_bytes = (n_seq + 7ULL) >> 3;
        for (uint32_t r = 0; r < n_seq; ++r) {
            out->has_qual[r] = (uint8_t)((blob[cur + (r >> 3)] >> (r & 7U)) & 1U);
        }
        cur += qual_bm_bytes;
    }

    const uint8_t *dna2 = blob + cur;
    cur += sum_dna2_bytes;

    const uint8_t *comp_names = blob + cur;
    cur += comp_name_bytes;

    const uint8_t *comp_n = blob + cur;
    cur += comp_n_bytes;

    const uint8_t *comp_case = blob + cur;
    cur += comp_case_bytes;

    const uint8_t *comp_exc = blob + cur;
    cur += comp_exc_bytes;

    /* Unpack 2-bit DNA into per-read buffers */
    uint64_t dna2_off = 0;
    for (uint32_t r = 0; r < n_seq; ++r) {
        uint32_t sl = (uint32_t)out->seq_lens[r];
        char *dst = (char *)malloc((size_t)sl + 4U);
        uint32_t nbytes = (sl + 3U) >> 2;
        for (uint32_t b = 0; b < nbytes; ++b) {
            memcpy(dst + ((size_t)b << 2), g_dna_decode_tab[dna2[dna2_off + b]], 4);
        }
        dst[sl] = '\0';
        out->seqs[r] = dst;
        dna2_off += nbytes;
    }
    assert(dna2_off == sum_dna2_bytes);

    /* Restore exact N positions */
    if (total_n_pos > 0) {
        uint64_t raw_n_bytes = ((uint64_t)n_seq + total_n_pos) * sizeof(uint32_t);
        uint8_t *raw_n_buf = (uint8_t *)malloc((size_t)raw_n_bytes);
        zlib_decompress_buf(comp_n, comp_n_bytes, raw_n_buf, raw_n_bytes);
        const uint8_t *ncnt_ptr = raw_n_buf;
        const uint8_t *npos_ptr = raw_n_buf + (size_t)n_seq * sizeof(uint32_t);
        uint64_t pidx = 0;
        for (uint32_t r = 0; r < n_seq; ++r) {
            uint32_t nc = get_u32(ncnt_ptr + (size_t)r * 4U);
            char *dst = out->seqs[r];
            for (uint32_t k = 0; k < nc; ++k) {
                uint32_t pos = get_u32(npos_ptr + (size_t)(pidx++) * 4U);
                dst[pos] = 'N';
            }
        }
        assert(pidx == total_n_pos);
        free(raw_n_buf);
    }

    /* Restore lowercase bases if present */
    if (comp_case_bytes > 0) {
        uint64_t case_bm_bytes = (sum_seq_len + 7ULL) >> 3;
        uint8_t *case_bm = (uint8_t *)malloc((size_t)case_bm_bytes);
        zlib_decompress_buf(comp_case, comp_case_bytes, case_bm, case_bm_bytes);
        uint64_t gbi = 0;
        for (uint32_t r = 0; r < n_seq; ++r) {
            uint32_t sl = (uint32_t)out->seq_lens[r];
            char *dst = out->seqs[r];
            for (uint32_t p = 0; p < sl; ++p, ++gbi) {
                if ((case_bm[gbi >> 3] >> (gbi & 7U)) & 1U) {
                    dst[p] = (char)((uint8_t)dst[p] | 0x20U);
                }
            }
        }
        free(case_bm);
    }

    /* Restore rare IUPAC/non-ACGTN bytes if present */
    if (total_exc > 0) {
        uint64_t raw_exc_bytes = total_exc * 9ULL;
        uint8_t *raw_exc_buf = (uint8_t *)malloc((size_t)raw_exc_bytes);
        zlib_decompress_buf(comp_exc, comp_exc_bytes, raw_exc_buf, raw_exc_bytes);
        for (uint64_t e = 0; e < total_exc; ++e) {
            uint32_t r = get_u32(raw_exc_buf + e * 9ULL);
            uint32_t pos = get_u32(raw_exc_buf + e * 9ULL + 4ULL);
            uint8_t ch = raw_exc_buf[e * 9ULL + 8ULL];
            out->seqs[r][pos] = (char)ch;
        }
        free(raw_exc_buf);
    }

    /* Decompress read names into null-terminated arena */
    out->names_buf = (char *)malloc((size_t)running_noff + 1U);
    if (sum_name_len > 0) {
        uint8_t *raw_names = (uint8_t *)malloc((size_t)sum_name_len);
        zlib_decompress_buf(comp_names, comp_name_bytes, raw_names, sum_name_len);
        uint64_t src_off = 0;
        for (uint32_t r = 0; r < n_seq; ++r) {
            uint32_t nl = out->name_lens[r];
            if (nl > 0) {
                memcpy(out->names_buf + out->name_offs[r], raw_names + src_off, nl);
                src_off += nl;
            }
            out->names_buf[out->name_offs[r] + nl] = '\0';
        }
        free(raw_names);
    } else {
        for (uint32_t r = 0; r < n_seq; ++r) {
            out->names_buf[out->name_offs[r]] = '\0';
        }
    }

    /* Multi-threaded decompression of base qualities if present */
    if (n_qual_slices > 0 && sum_qual_len > 0) {
        uint32_t *slice_clens = (uint32_t *)malloc((size_t)n_qual_slices * sizeof(uint32_t));
        uint32_t *slice_rlens = (uint32_t *)malloc((size_t)n_qual_slices * sizeof(uint32_t));
        uint64_t *slice_coffs = (uint64_t *)malloc((size_t)n_qual_slices * sizeof(uint64_t));
        uint64_t coff = 0;
        for (uint32_t s = 0; s < n_qual_slices; ++s) {
            slice_clens[s] = get_u32(blob + cur + (uint64_t)s * 8ULL);
            slice_rlens[s] = get_u32(blob + cur + (uint64_t)s * 8ULL + 4ULL);
            slice_coffs[s] = coff;
            coff += slice_clens[s];
        }
        assert(coff == comp_qual_bytes);
        cur += (uint64_t)n_qual_slices * 8ULL;

        uint8_t *raw_qual = (uint8_t *)malloc((size_t)sum_qual_len + 1U);
        qual_decomp_ctx_t dctx;
        dctx.comp_base = blob + cur;
        dctx.slice_comp_offs = slice_coffs;
        dctx.slice_comp_lens = slice_clens;
        dctx.slice_raw_lens = slice_rlens;
        dctx.dst_qual = raw_qual;

        int nth = (n_threads > 1) ? n_threads : 1;
        if ((uint32_t)nth > n_qual_slices) nth = (int)n_qual_slices;
        kt_for(nth, worker_qual_decompress, &dctx, (long)n_qual_slices);

        free(slice_clens);
        free(slice_rlens);
        free(slice_coffs);

        uint32_t total_qarena = (uint32_t)(sum_qual_len + n_seq + 1ULL);
        out->quals_buf = (char *)malloc((size_t)total_qarena);
        uint64_t qsrc = 0;
        uint32_t qdst = 0;
        for (uint32_t r = 0; r < n_seq; ++r) {
            out->qual_offs[r] = qdst;
            if (out->has_qual[r]) {
                uint32_t ql = (uint32_t)out->seq_lens[r];
                if (ql > 0) {
                    memcpy(out->quals_buf + qdst, raw_qual + qsrc, ql);
                    qsrc += ql;
                    qdst += ql;
                }
            }
            out->quals_buf[qdst++] = '\0';
        }
        assert(qsrc == sum_qual_len);
        free(raw_qual);
    }

    if (free_blob) free(blob);
    return 1;
}

void ha_stream_decoded_batch_free_aux(ha_stream_decoded_batch_t *b)
{
    if (!b) return;
    free(b->name_lens);
    free(b->name_offs);
    free(b->names_buf);
    free(b->has_qual);
    free(b->qual_offs);
    free(b->quals_buf);
    b->name_lens = NULL;
    b->name_offs = NULL;
    b->names_buf = NULL;
    b->has_qual = NULL;
    b->qual_offs = NULL;
    b->quals_buf = NULL;
}

uint64_t ha_stream_store_ram_bytes(const ha_stream_store_t *store)
{
    return store ? store->ram_used_bytes : 0;
}

uint64_t ha_stream_store_peak_ram_bytes(const ha_stream_store_t *store)
{
    return store ? store->peak_ram_bytes : 0;
}

uint64_t ha_stream_store_disk_bytes(const ha_stream_store_t *store)
{
    return store ? store->disk_used_bytes : 0;
}

uint64_t ha_stream_store_total_records(const ha_stream_store_t *store)
{
    return store ? store->total_records : 0;
}

uint64_t ha_stream_store_total_bases(const ha_stream_store_t *store)
{
    return store ? store->total_bases : 0;
}
