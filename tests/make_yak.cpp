#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "../khashl.h"
#include "../kseq.h"

KSEQ_INIT(gzFile, gzread)

#define YAK_K            31
#define YAK_PRE          10
#define YAK_COUNTER_BITS 10
#define YAK_MAX_COUNT    ((1 << YAK_COUNTER_BITS) - 1)
#define YAK_MAGIC        "YAK\2"

#define yak_ch_eq(a, b) ((a) >> YAK_COUNTER_BITS == (b) >> YAK_COUNTER_BITS)
#define yak_ch_hash(a)  ((a) >> YAK_COUNTER_BITS)
KHASHL_SET_INIT(static klib_unused, yak_ht_t, yak_ht, uint64_t, yak_ch_hash, yak_ch_eq)

static const unsigned char nt4_tab[256] = {
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 3,3,4,4, 4,4,4,4, 4,4,4,4,
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 3,3,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4
};

static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
    key = (~key + (key << 21)) & mask;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "Usage: make_yak <out.yak> <in1.fa/fq[.gz]> [in2 ...]\n");
        return 1;
    }
    const int n_buckets = 1 << YAK_PRE;
    const uint64_t pre_mask = (uint64_t)n_buckets - 1ULL;
    const uint64_t kmask = (1ULL << (2 * YAK_K)) - 1ULL;
    const int shift = 2 * (YAK_K - 1);

    yak_ht_t **buckets = (yak_ht_t **)calloc((size_t)n_buckets, sizeof(yak_ht_t *));
    for (int i = 0; i < n_buckets; ++i) buckets[i] = yak_ht_init();

    for (int f = 2; f < argc; ++f) {
        gzFile fp = gzopen(argv[f], "r");
        if (!fp) {
            fprintf(stderr, "Failed to open %s\n", argv[f]);
            return 1;
        }
        kseq_t *ks = kseq_init(fp);
        while (kseq_read(ks) >= 0) {
            uint64_t x0 = 0, x1 = 0;
            int l = 0;
            for (size_t i = 0; i < ks->seq.l; ++i) {
                unsigned char c = nt4_tab[(uint8_t)ks->seq.s[i]];
                if (c < 4) {
                    x0 = ((x0 << 2) | c) & kmask;
                    x1 = (x1 >> 2) | ((uint64_t)(3U - c) << shift);
                    if (++l >= YAK_K) {
                        uint64_t canon = (x0 < x1) ? x0 : x1;
                        uint64_t h = hash64(canon, kmask);
                        yak_ht_t *g = buckets[h & pre_mask];
                        uint64_t base_key = (h >> YAK_PRE) << YAK_COUNTER_BITS;
                        int absent = 0;
                        khint_t k = yak_ht_put(g, base_key, &absent);
                        if ((kh_key(g, k) & YAK_MAX_COUNT) < YAK_MAX_COUNT) {
                            ++kh_key(g, k);
                        }
                    }
                } else {
                    l = 0;
                    x0 = x1 = 0;
                }
            }
        }
        kseq_destroy(ks);
        gzclose(fp);
    }

    FILE *out = fopen(argv[1], "wb");
    if (!out) {
        fprintf(stderr, "Failed to open output %s\n", argv[1]);
        return 1;
    }
    fwrite(YAK_MAGIC, 1, 4, out);
    uint32_t hdr[3] = { YAK_K, YAK_PRE, YAK_COUNTER_BITS };
    fwrite(hdr, 4, 3, out);
    for (int i = 0; i < n_buckets; ++i) {
        yak_ht_t *g = buckets[i];
        uint32_t bt[2] = { kh_capacity(g), kh_size(g) };
        fwrite(bt, 4, 2, out);
        for (khint_t k = 0; k < kh_end(g); ++k) {
            if (kh_exist(g, k)) {
                uint64_t key = kh_key(g, k);
                fwrite(&key, 8, 1, out);
            }
        }
        yak_ht_destroy(g);
    }
    free(buckets);
    fclose(out);
    return 0;
}
