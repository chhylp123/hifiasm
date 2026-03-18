#ifndef HAF_SEQ_READER_H
#define HAF_SEQ_READER_H

#include <stddef.h>
#include <stdint.h>
#include <zlib.h>

typedef struct {
	char *s;
	size_t l, m;
} ha_sr_buf_t;

typedef struct {
	ha_sr_buf_t name, seq, qual;
} ha_sr_rec_t;

typedef enum {
	HA_SR_FASTX = 0,
	HA_SR_BAM = 1
} ha_sr_mode_t;

typedef struct {
	gzFile fp;
	void *ks;
	ha_sr_rec_t rec;
	ha_sr_mode_t mode;
	char *bam_buf;
	size_t bam_buf_cap;
	char err[256];
	int bam_header_ready;
} ha_seq_reader_t;

int ha_seq_open(ha_seq_reader_t *r, const char *fn);
int ha_seq_read(ha_seq_reader_t *r);
void ha_seq_close(ha_seq_reader_t *r);
const char *ha_seq_error(const ha_seq_reader_t *r);

#endif
