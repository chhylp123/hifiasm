#ifndef __HA_STREAM_STORE_H__
#define __HA_STREAM_STORE_H__

#include <stdint.h>
#include <stddef.h>
#include <zlib.h>

typedef struct {
    char *name;
    uint32_t name_len;
    char *qual;
    uint32_t qual_len;
    uint8_t has_qual;
} ha_stream_raw_rec_t;

typedef struct {
    uint32_t n_seq;
    int *seq_lens;
    char **seqs;
    uint32_t *name_lens;
    uint32_t *name_offs;
    char *names_buf;
    uint8_t *has_qual;
    uint32_t *qual_offs;
    char *quals_buf;
} ha_stream_decoded_batch_t;

struct ha_stream_store_s;
typedef struct ha_stream_store_s ha_stream_store_t;

int ha_is_stream_path(const char *fn);
int ha_check_input_readable(const char *fn, const char *opt_label);
gzFile ha_gzopen_read(const char *fn);

ha_stream_store_t *ha_stream_store_create(int num_files, double max_ram_gb, const char *fallback_path);
void ha_stream_store_destroy(ha_stream_store_t **store_ptr);

int ha_stream_store_is_captured(const ha_stream_store_t *store, int file_idx);
void ha_stream_store_mark_captured(ha_stream_store_t *store, int file_idx);
void ha_stream_store_rewind_file(ha_stream_store_t *store, int file_idx);

void ha_stream_store_append_batch(
    ha_stream_store_t *store,
    int file_idx,
    int n_seq,
    char **seqs,
    const int *seq_lens,
    const ha_stream_raw_rec_t *raw_recs,
    int n_threads
);

int ha_stream_store_replay_batch(
    ha_stream_store_t *store,
    int file_idx,
    int consume_ram,
    int n_threads,
    ha_stream_decoded_batch_t *out
);

void ha_stream_decoded_batch_free_aux(ha_stream_decoded_batch_t *b);

uint64_t ha_stream_store_ram_bytes(const ha_stream_store_t *store);
uint64_t ha_stream_store_peak_ram_bytes(const ha_stream_store_t *store);
uint64_t ha_stream_store_disk_bytes(const ha_stream_store_t *store);
uint64_t ha_stream_store_total_records(const ha_stream_store_t *store);
uint64_t ha_stream_store_total_bases(const ha_stream_store_t *store);

void ha_stream_global_destroy(void);

#endif
