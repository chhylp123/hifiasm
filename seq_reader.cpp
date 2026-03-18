#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include "kseq.h"
#include "seq_reader.h"

KSEQ_INIT(gzFile, gzread)

static void sr_set_error(ha_seq_reader_t *r, const char *msg)
{
	snprintf(r->err, sizeof(r->err), "%s", msg);
}

static void sr_set_error2(ha_seq_reader_t *r, const char *prefix, const char *fn)
{
	snprintf(r->err, sizeof(r->err), "%s%s", prefix, fn);
}

static int sr_gzread_exact(ha_seq_reader_t *r, void *buf, size_t len)
{
	size_t off = 0;
	while (off < len) {
		int ret = gzread(r->fp, (char *)buf + off, (unsigned int)(len - off));
		if (ret <= 0) return 0;
		off += ret;
	}
	return 1;
}

static uint32_t sr_le_u32(const uint8_t *p)
{
	return ((uint32_t)p[0]) | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}

static int32_t sr_le_i32(const uint8_t *p)
{
	return (int32_t)sr_le_u32(p);
}

static int sr_reserve(ha_sr_buf_t *b, size_t need)
{
	char *tmp;
	if (b->m >= need) return 1;
	b->m = need;
	tmp = (char *)realloc(b->s, b->m);
	if (!tmp) return 0;
	b->s = tmp;
	return 1;
}

static int sr_prepare_bam(ha_seq_reader_t *r)
{
	uint8_t hdr[8];
	int32_t l_text, n_ref, i;
	if (r->bam_header_ready) return 1;
	if (!sr_gzread_exact(r, hdr, 4)) {
		sr_set_error(r, "failed to read BAM magic");
		return 0;
	}
	if (memcmp(hdr, "BAM\1", 4) != 0) {
		sr_set_error(r, "input is not a BAM file");
		return 0;
	}
	if (!sr_gzread_exact(r, hdr, 4)) {
		sr_set_error(r, "failed to read BAM header length");
		return 0;
	}
	l_text = sr_le_i32(hdr);
	if (l_text < 0) {
		sr_set_error(r, "corrupted BAM header");
		return 0;
	}
	if (l_text > 0) {
		if (r->bam_buf_cap < (size_t)l_text) {
			char *tmp = (char *)realloc(r->bam_buf, (size_t)l_text);
			if (!tmp) {
				sr_set_error(r, "out of memory while reading BAM header");
				return 0;
			}
			r->bam_buf = tmp;
			r->bam_buf_cap = (size_t)l_text;
		}
		if (!sr_gzread_exact(r, r->bam_buf, (size_t)l_text)) {
			sr_set_error(r, "failed to read BAM text header");
			return 0;
		}
	}
	if (!sr_gzread_exact(r, hdr, 4)) {
		sr_set_error(r, "failed to read BAM reference count");
		return 0;
	}
	n_ref = sr_le_i32(hdr);
	if (n_ref < 0) {
		sr_set_error(r, "corrupted BAM reference count");
		return 0;
	}
	for (i = 0; i < n_ref; ++i) {
		int32_t l_name;
		if (!sr_gzread_exact(r, hdr, 4)) {
			sr_set_error(r, "failed to read BAM reference header");
			return 0;
		}
		l_name = sr_le_i32(hdr);
		if (l_name <= 0) {
			sr_set_error(r, "corrupted BAM reference name");
			return 0;
		}
		if (r->bam_buf_cap < (size_t)l_name) {
			char *tmp = (char *)realloc(r->bam_buf, (size_t)l_name);
			if (!tmp) {
				sr_set_error(r, "out of memory while reading BAM references");
				return 0;
			}
			r->bam_buf = tmp;
			r->bam_buf_cap = (size_t)l_name;
		}
		if (!sr_gzread_exact(r, r->bam_buf, (size_t)l_name) || !sr_gzread_exact(r, hdr, 4)) {
			sr_set_error(r, "failed to read BAM reference entry");
			return 0;
		}
	}
	r->bam_header_ready = 1;
	return 1;
}

static int sr_read_bam(ha_seq_reader_t *r)
{
	static const char bam_nt16_rev_table[16] = {
		'N', 'A', 'C', 'N', 'G', 'N', 'N', 'N', 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
	};
	uint8_t core[36], *seq, *qual;
	int32_t block_size, ref_id, pos, l_seq;
	uint32_t bin_mq_nl, flag_nc, l_read_name, n_cigar, flag;
	size_t seq_bytes, data_len, qual_off, i;
	int ret = gzread(r->fp, core, 4);
	if (ret == 0) return -1;
	if (ret != 4 || !sr_gzread_exact(r, core + 4, 32)) {
		sr_set_error(r, "truncated BAM record");
		return -2;
	}
	block_size = sr_le_i32(core);
	ref_id = sr_le_i32(core + 4);
	pos = sr_le_i32(core + 8);
	bin_mq_nl = sr_le_u32(core + 12);
	flag_nc = sr_le_u32(core + 16);
	l_seq = sr_le_i32(core + 20);
	if (block_size < 32 || l_seq < 0) {
		sr_set_error(r, "corrupted BAM record");
		return -2;
	}
	l_read_name = bin_mq_nl & 0xffU;
	n_cigar = flag_nc & 0xffffU;
	flag = flag_nc >> 16;
	data_len = (size_t)block_size - 32;
	if (r->bam_buf_cap < data_len) {
		char *tmp = (char *)realloc(r->bam_buf, data_len);
		if (!tmp) {
			sr_set_error(r, "out of memory while reading BAM record");
			return -2;
		}
		r->bam_buf = tmp;
		r->bam_buf_cap = data_len;
	}
	if (!sr_gzread_exact(r, r->bam_buf, data_len)) {
		sr_set_error(r, "truncated BAM payload");
		return -2;
	}
	if (!(flag & 0x4U) || ref_id != -1 || pos != -1 || n_cigar != 0) {
		sr_set_error(r, "only unaligned BAM input is supported");
		return -2;
	}
	seq_bytes = ((size_t)l_seq + 1) >> 1;
	qual_off = (size_t)l_read_name + ((size_t)n_cigar << 2) + seq_bytes;
	if (data_len < qual_off + (size_t)l_seq) {
		sr_set_error(r, "corrupted BAM read layout");
		return -2;
	}
	if (!sr_reserve(&r->rec.name, l_read_name ? (size_t)l_read_name : 1) ||
		!sr_reserve(&r->rec.seq, (size_t)l_seq + 1) ||
		!sr_reserve(&r->rec.qual, (size_t)l_seq + 1)) {
		sr_set_error(r, "out of memory while decoding BAM read");
		return -2;
	}
	r->rec.name.l = l_read_name ? (size_t)l_read_name - 1 : 0;
	memcpy(r->rec.name.s, r->bam_buf, r->rec.name.l);
	r->rec.name.s[r->rec.name.l] = 0;
	seq = (uint8_t *)r->bam_buf + l_read_name + ((size_t)n_cigar << 2);
	for (i = 0; i < (size_t)l_seq; ++i) {
		uint8_t c = seq[i >> 1];
		r->rec.seq.s[i] = bam_nt16_rev_table[(i & 1)? (c & 0xf) : (c >> 4)];
	}
	r->rec.seq.l = (size_t)l_seq;
	r->rec.seq.s[r->rec.seq.l] = 0;
	qual = (uint8_t *)r->bam_buf + qual_off;
	r->rec.qual.l = (size_t)l_seq;
	for (i = 0; i < (size_t)l_seq; ++i) {
		if (qual[i] == 0xff) {
			r->rec.qual.l = 0;
			break;
		}
		r->rec.qual.s[i] = (char)(qual[i] + 33);
	}
	if (r->rec.qual.l) r->rec.qual.s[r->rec.qual.l] = 0;
	return (int)r->rec.seq.l;
}

int ha_seq_open(ha_seq_reader_t *r, const char *fn)
{
	uint8_t magic[4];
	memset(r, 0, sizeof(*r));
	r->fp = gzopen(fn, "r");
	if (!r->fp) {
		sr_set_error2(r, "failed to open input file: ", fn);
		return 0;
	}
	if (gzread(r->fp, magic, 4) != 4) {
		sr_set_error2(r, "failed to read input file: ", fn);
		ha_seq_close(r);
		return 0;
	}
	if (gzrewind(r->fp) != 0) {
		sr_set_error2(r, "failed to rewind input file: ", fn);
		ha_seq_close(r);
		return 0;
	}
	if (memcmp(magic, "BAM\1", 4) == 0) {
		r->mode = HA_SR_BAM;
		if (!sr_prepare_bam(r)) {
			if (r->fp) gzclose(r->fp);
			r->fp = 0;
			return 0;
		}
	} else {
		r->mode = HA_SR_FASTX;
		r->ks = kseq_init(r->fp);
		if (!r->ks) {
			sr_set_error(r, "failed to initialize FASTA/FASTQ reader");
			ha_seq_close(r);
			return 0;
		}
	}
	return 1;
}

int ha_seq_read(ha_seq_reader_t *r)
{
	if (r->mode == HA_SR_BAM) return sr_read_bam(r);
	if (r->ks) {
		kseq_t *ks = (kseq_t *)r->ks;
		int ret = kseq_read(ks);
		if (ret < 0) return ret;
		r->rec.name.s = ks->name.s;
		r->rec.name.l = ks->name.l;
		r->rec.name.m = ks->name.m;
		r->rec.seq.s = ks->seq.s;
		r->rec.seq.l = ks->seq.l;
		r->rec.seq.m = ks->seq.m;
		r->rec.qual.s = ks->qual.s;
		r->rec.qual.l = ks->qual.l;
		r->rec.qual.m = ks->qual.m;
		return ret;
	}
	sr_set_error(r, "reader is not initialized");
	return -2;
}

void ha_seq_close(ha_seq_reader_t *r)
{
	if (r->ks) kseq_destroy((kseq_t *)r->ks);
	if (r->fp) gzclose(r->fp);
	if (r->mode == HA_SR_BAM) {
		free(r->rec.name.s);
		free(r->rec.seq.s);
		free(r->rec.qual.s);
	}
	free(r->bam_buf);
	r->fp = 0;
	r->ks = 0;
	r->bam_buf = 0;
	r->bam_buf_cap = 0;
	r->bam_header_ready = 0;
	r->rec.name.s = r->rec.seq.s = r->rec.qual.s = 0;
	r->rec.name.l = r->rec.name.m = 0;
	r->rec.seq.l = r->rec.seq.m = 0;
	r->rec.qual.l = r->rec.qual.m = 0;
}

const char *ha_seq_error(const ha_seq_reader_t *r)
{
	return r->err;
}
