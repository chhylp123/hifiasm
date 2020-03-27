#include <stdio.h>
#include "htab.h"
#include "ksort.h"
#include "Hash_Table.h"

typedef struct { // this struct is not strictly necessary; we can use k_mer_pos instead, with modifications
	uint64_t srt;
	uint32_t self_off;
	uint32_t other_off;
} anchor1_t;

#define an_key1(a) ((a).srt)
#define an_key2(a) ((a).self_off)
KRADIX_SORT_INIT(ha_an1, anchor1_t, an_key1, 8)
KRADIX_SORT_INIT(ha_an2, anchor1_t, an_key2, 4)

typedef struct {
	int n;
	const ha_idxpos_t *a;
} seed1_t;

struct ha_abuf_s {
	uint64_t n_a, m_a;
	uint32_t old_mz_m;
	ha_mz1_v mz;
	seed1_t *seed;
	anchor1_t *a;
};

ha_abuf_t *ha_abuf_init(void)
{
	return (ha_abuf_t*)calloc(1, sizeof(ha_abuf_t));
}

void ha_abuf_destroy(ha_abuf_t *ab)
{
	free(ab->seed); free(ab->a); free(ab->mz.a); free(ab);
}

void ha_get_new_candidates(ha_abuf_t *ab, int64_t rid, UC_Read *ucr, overlap_region_alloc *overlap_list, Candidates_list *cl, double band_width_threshold, int keep_whole_chain)
{
	extern void *ha_flt_tab;
	extern ha_pt_t *ha_idx;
	uint32_t i;
	uint64_t k, l;

	// prepare
    clear_Candidates_list(cl);
    clear_overlap_region_alloc(overlap_list);
	recover_UC_Read(ucr, &R_INF, rid);
	ab->mz.n = 0, ab->n_a = 0;

	// get the list of anchors
	ha_sketch(ucr->seq, ucr->length, asm_opt.mz_win, asm_opt.k_mer_length, 0, !asm_opt.no_HPC, &ab->mz, ha_flt_tab);
	if (ab->mz.m > ab->old_mz_m) {
		ab->old_mz_m = ab->mz.m;
		REALLOC(ab->seed, ab->old_mz_m);
	}
	for (i = 0, ab->n_a = 0; i < ab->mz.n; ++i) {
		ab->seed[i].a = ha_pt_get(ha_idx, ab->mz.a[i].x, &ab->seed[i].n);
		ab->n_a += ab->seed[i].n;
	}
	if (ab->n_a > ab->m_a) {
		ab->m_a = ab->n_a;
		ab->m_a = ab->m_a > 16? ab->m_a + (ab->m_a>>1) : 16;
		REALLOC(ab->a, ab->m_a);
	}
	for (i = 0, k = 0; i < ab->mz.n; ++i) {
		int j;
		ha_mz1_t *z = &ab->mz.a[i];
		seed1_t *s = &ab->seed[i];
		for (j = 0; j < s->n; ++j) {
			const ha_idxpos_t *y = &s->a[j];
			anchor1_t *an = &ab->a[k++];
			uint8_t rev = z->rev == y->rev? 0 : 1;
			an->other_off = y->pos;
			an->self_off = rev? ucr->length - 1 - (z->pos + 1 - z->span) : z->pos;
			an->srt = (uint64_t)y->rid<<33 | (uint64_t)rev<<32 | an->other_off;
		}
	}

	// sort anchors
	radix_sort_ha_an1(ab->a, ab->a + ab->n_a);
	for (k = 1, l = 0; k <= ab->n_a; ++k) {
		if (k == ab->n_a || ab->a[k].srt != ab->a[l].srt) {
			if (k - l > 1)
				radix_sort_ha_an2(ab->a + l, ab->a + k);
			l = k;
		}
	}

	// copy over to _cl_
	if (ab->m_a >= (uint64_t)cl->size) {
		cl->size = ab->m_a;
		REALLOC(cl->list, cl->size);
	}
	for (k = 0; k < ab->n_a; ++k) {
		k_mer_hit *p = &cl->list[k];
		p->readID = ab->a[k].srt >> 33;
		p->strand = ab->a[k].srt >> 32 & 1;
		p->offset = ab->a[k].other_off;
		p->self_offset = ab->a[k].self_off;
	}
	cl->length = ab->n_a;

	calculate_overlap_region_by_chaining(cl, overlap_list, rid, ucr->length, &R_INF, band_width_threshold, keep_whole_chain);
	#if 0
	fprintf(stderr, "B\t%ld\t%ld\n", (long)overlap_list->length, (long)overlap_list->mapped_overlaps_length);
	for (int i = 0; i < (int)overlap_list->length; ++i) {
		overlap_region *r = &overlap_list->list[i];
		fprintf(stderr, "C\t%d\t%d\t%d\t%c\t%d\t%d\t%d\t%c\t%d\t%d\n", (int)r->x_id, (int)r->x_pos_s, (int)r->x_pos_e, "+-"[r->x_pos_strand],
				(int)r->y_id, (int)r->y_pos_s, (int)r->y_pos_e, "+-"[r->y_pos_strand], (int)r->shared_seed, r->is_match);
	}
	#endif
}
