#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "kvec.h"
#include "yak.h"

typedef struct { // a simplified version of kdq
	int front, count;
	int a[64];
} tiny_queue_t;

static inline void tq_push(tiny_queue_t *q, int x)
{
	q->a[((q->count++) + q->front) & 0x3f] = x;
}

static inline int tq_shift(tiny_queue_t *q)
{
	int x;
	if (q->count == 0) return -1;
	x = q->a[q->front++];
	q->front &= 0x3f;
	--q->count;
	return x;
}

/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      minimizers
 */
void ha_sketch(const char *str, int len, int w, int k, uint32_t rid, int is_hpc, ha_mz1_v *p, void *hf)
{
	static const ha_mz1_t dummy = { UINT64_MAX, 0, 0, 0 };
	uint64_t shift1 = k - 1, mask = (1ULL<<k) - 1, kmer[4] = {0,0,0,0};
	int i, j, l, buf_pos, min_pos, kmer_span = 0;
	ha_mz1_t buf[256], min = dummy;
	tiny_queue_t tq;

	assert(len > 0 && len < 1<<27 && rid < 1<<28 && (w > 0 && w < 256) && (k > 0 && k <= 63));
	memset(buf, 0xff, w * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(ha_mz1_t, *p, p->n + len/w);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		ha_mz1_t info = dummy;
		if (c < 4) { // not an ambiguous base
			int z;
			if (is_hpc) {
				int skip_len = 1;
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
					for (skip_len = 2; i + skip_len < len; ++skip_len)
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
							break;
					i += skip_len - 1; // put $i at the end of the current homopolymer run
				}
				tq_push(&tq, skip_len);
				kmer_span += skip_len;
				if (tq.count > k) kmer_span -= tq_shift(&tq);
			} else kmer_span = l + 1 < k? l + 1 : k;
			kmer[0] = (kmer[0] << 1 | (c&1))  & mask;                  // forward k-mer
			kmer[1] = (kmer[1] << 1 | (c>>1)) & mask;
			kmer[2] = kmer[2] >> 1 | (uint64_t)(1 - (c&1))  << shift1; // reverse k-mer
			kmer[3] = kmer[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift1;
			if (kmer[1] == kmer[3]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[1] < kmer[3]? 0 : 1; // strand
			++l;
			if (l >= k && kmer_span < 256) {
				uint64_t y;
				y = yak_hash64_64(kmer[z<<1|0]) + yak_hash64_64(kmer[z<<1|1]);
				if (hf == 0 || ha_hf_isflt(hf, y) == 0)
					info.x = y, info.rid = rid, info.pos = i, info.rev = z, info.span = kmer_span;
			}
		} else l = 0, tq.count = tq.front = 0, kmer_span = 0;
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].pos != min.pos) kv_push(ha_mz1_t, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].pos != min.pos) kv_push(ha_mz1_t, *p, buf[j]);
		}
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k && min.x != UINT64_MAX) kv_push(ha_mz1_t, *p, min);
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1 && min.x != UINT64_MAX) kv_push(ha_mz1_t, *p, min);
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.pos != buf[j].pos) kv_push(ha_mz1_t, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.pos != buf[j].pos) kv_push(ha_mz1_t, *p, buf[j]);
			}
		}
		if (++buf_pos == w) buf_pos = 0;
	}
	if (min.x != UINT64_MAX)
		kv_push(ha_mz1_t, *p, min);
}
