#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "htab.h"
#include "ksort.h"
#include "Hash_Table.h"
#include "kalloc.h"
#include "Overlaps.h"
#include "Levenshtein_distance.h"

#define HA_KMER_GOOD_RATIO 0.333
#define OFL 0.95
#define CH_OCC 4
#define CH_SC 16

typedef struct { // this struct is not strictly necessary; we can use k_mer_pos instead, with modifications
	uint64_t srt;
	uint32_t self_off;
	uint32_t other_off;
	uint32_t cnt;
} anchor1_t;

#define an_key1(a) ((a).srt)
#define an_key2(a) ((a).self_off)
#define an_key3(a) ((a).other_off)
KRADIX_SORT_INIT(ha_an1, anchor1_t, an_key1, 8)
KRADIX_SORT_INIT(ha_an2, anchor1_t, an_key2, 4)
KRADIX_SORT_INIT(ha_an3, anchor1_t, an_key3, 4)
#define generic_key(x) (x)
KRADIX_SORT_INIT(anc64, uint64_t, generic_key, 8)

#define oreg_xs_lt(a, b) (((uint64_t)(a).x_pos_s<<32|(a).x_pos_e) < ((uint64_t)(b).x_pos_s<<32|(b).x_pos_e))
KSORT_INIT(or_xs, overlap_region, oreg_xs_lt)

#define oreg_ss_lt(a, b) ((a).shared_seed > (b).shared_seed) // in the decending order
KSORT_INIT(or_ss, overlap_region, oreg_ss_lt)

#define oreg_occ_lt(a, b) ((a).align_length > (b).align_length) // in the decending order
KSORT_INIT(or_occ, overlap_region, oreg_occ_lt)

#define oreg_id_lt(a, b) ((a).y_id < (b).y_id)
KSORT_INIT(or_id, overlap_region, oreg_id_lt)

#define ha_mz1_t_key(p) ((p).x)
KRADIX_SORT_INIT(ha_mz1_v_srt, ha_mz1_t, ha_mz1_t_key, member_size(ha_mz1_t, x))

typedef struct {
	int n;
	const ha_idxpos_t *a;
} seed1_t;

typedef struct {
	int n, cnt;
	const ha_idxposl_t *a;
} seedl_t;

struct ha_abuf_s {
	uint64_t n_a, m_a;///number of anchors (seed positions)
	uint32_t old_mz_m;///number of seeds
	ha_mz1_v mz;
	seed1_t *seed;
	anchor1_t *a;
};

struct ha_abufl_s {
	uint64_t n_a, m_a;///number of anchors (seed positions)
	uint32_t old_mz_m;///number of seeds
	ha_mzl_v mz;
	seedl_t *seed;
	anchor1_t *a;
};

#define HA_ABUF_INIT(HType, MZType, SDType, sf) \
HType *sf##_init_buf(void *km){HType *b = NULL; KCALLOC((km), b, 1); return b;}\
HType *sf##_init(void){return (HType*)calloc(1, sizeof(HType));}\
void sf##_free_buf(void *km, HType *ab, int is_z){if(ab){kfree(km, ab->seed); kfree(km, ab->a); kfree(km, ab->mz.a); if((is_z)){memset(ab, 0, sizeof(*ab));}}}\
void sf##_destroy_buf(void *km, HType *ab){if(ab){kfree(km, ab->seed); kfree(km, ab->a); kfree(km, ab->mz.a); kfree(km, ab);}}\
void sf##_destroy(HType *ab){if(ab){free(ab->seed); free(ab->a); free(ab->mz.a); free(ab);}}\
uint64_t sf##_mem(const HType *ab){\
	return ab->m_a * sizeof(anchor1_t) + ab->mz.m * (sizeof(MZType) + sizeof(SDType)) + sizeof(HType);\
}

HA_ABUF_INIT(ha_abuf_s, ha_mz1_t, seed1_t, ha_abuf)
HA_ABUF_INIT(ha_abufl_s, ha_mzl_t, seedl_t, ha_abufl)

int ha_ov_type(const overlap_region *r, uint32_t len)
{
	if (r->x_pos_s == 0 && r->x_pos_e == len - 1) return 2; // contained in a longer read
	else if (r->x_pos_s > 0 && r->x_pos_e < len - 1) return 3; // containing a shorter read
	else return r->x_pos_s == 0? 0 : 1;
}

void ha_get_new_candidates(ha_abuf_t *ab, int64_t rid, UC_Read *ucr, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, int max_n_chain, int keep_whole_chain,
						   kvec_t_u8_warp* k_flag, kvec_t_u64_warp* chain_idx, void *ha_flt_tab, ha_pt_t *ha_idx, overlap_region* f_cigar, kvec_t_u64_warp* dbg_ct, st_mt_t *sp)
{
	uint32_t i, rlen;
	uint64_t k, l;
	uint32_t low_occ = asm_opt.hom_cov * HA_KMER_GOOD_RATIO;
	uint32_t high_occ = asm_opt.hom_cov * (2.0 - HA_KMER_GOOD_RATIO);
	if(low_occ < 2) low_occ = 2;

	// prepare
    clear_Candidates_list(cl);
    clear_overlap_region_alloc(overlap_list);
	recover_UC_Read(ucr, &R_INF, rid);
	ab->mz.n = 0, ab->n_a = 0;
	rlen = Get_READ_LENGTH(R_INF, rid); // read length

	// get the list of anchors
	mz1_ha_sketch(ucr->seq, ucr->length, asm_opt.mz_win, asm_opt.k_mer_length, 0, !(asm_opt.flag & HA_F_NO_HPC), &ab->mz, ha_flt_tab, asm_opt.mz_sample_dist, k_flag, dbg_ct, NULL, -1, asm_opt.dp_min_len, -1, sp, asm_opt.mz_rewin, 0, NULL);
	// minimizer of queried read
	if (ab->mz.m > ab->old_mz_m) {
		ab->old_mz_m = ab->mz.m;
		REALLOC(ab->seed, ab->old_mz_m);
	}
	for (i = 0, ab->n_a = 0; i < ab->mz.n; ++i) {
		int n;
		ab->seed[i].a = ha_pt_get(ha_idx, ab->mz.a[i].x, &n);
		ab->seed[i].n = n;
		ab->n_a += n;
	}
	if (ab->n_a > ab->m_a) {
		ab->m_a = ab->n_a;
		kroundup64(ab->m_a);
		REALLOC(ab->a, ab->m_a);
	}
	for (i = 0, k = 0; i < ab->mz.n; ++i) {
		int j;
		///z is one of the minimizer
		ha_mz1_t *z = &ab->mz.a[i];
		seed1_t *s = &ab->seed[i];
		for (j = 0; j < s->n; ++j) {
			const ha_idxpos_t *y = &s->a[j];
			anchor1_t *an = &ab->a[k++];
			uint8_t rev = z->rev == y->rev? 0 : 1;
			an->other_off = y->pos;
			an->self_off = rev? ucr->length - 1 - (z->pos + 1 - z->span) : z->pos;
			an->cnt = s->n;
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
		if(ab->a[k].cnt > low_occ && ab->a[k].cnt < high_occ){
			p->cnt = 1;
		}
		else if(ab->a[k].cnt <= low_occ){
			p->cnt = 2;
		}
		else{
			p->cnt = 1 + ((ab->a[k].cnt + (high_occ<<1) - 1)/(high_occ<<1));
			p->cnt = pow(p->cnt, 1.1);
		}
	}
	cl->length = ab->n_a;

	calculate_overlap_region_by_chaining(cl, overlap_list, chain_idx, rid, ucr->length, &R_INF, NULL, bw_thres, keep_whole_chain, f_cigar, NULL);

	#if 0
	if (overlap_list->length > 0) {
		fprintf(stderr, "B\t%ld\t%ld\t%d\n", (long)rid, (long)overlap_list->length, rlen);
		for (int i = 0; i < (int)overlap_list->length; ++i) {
			overlap_region *r = &overlap_list->list[i];
			fprintf(stderr, "C\t%d\t%d\t%d\t%c\t%d\t%ld\t%d\t%d\t%c\t%d\t%d\n", (int)r->x_id, (int)r->x_pos_s, (int)r->x_pos_e, "+-"[r->x_pos_strand],
					(int)r->y_id, (long)Get_READ_LENGTH(R_INF, r->y_id), (int)r->y_pos_s, (int)r->y_pos_e, "+-"[r->y_pos_strand], (int)r->shared_seed, ha_ov_type(r, rlen));
		}
	}
	#endif
	
	if ((int)overlap_list->length > max_n_chain) {
		int32_t w, n[4], s[4];
		n[0] = n[1] = n[2] = n[3] = 0, s[0] = s[1] = s[2] = s[3] = 0;
		ks_introsort_or_ss(overlap_list->length, overlap_list->list);
		for (i = 0; i < (uint32_t)overlap_list->length; ++i) {
			const overlap_region *r = &overlap_list->list[i];
			w = ha_ov_type(r, rlen);
			++n[w];
			if ((int)n[w] == max_n_chain) s[w] = r->shared_seed;
		}
		if (s[0] > 0 || s[1] > 0 || s[2] > 0 || s[3] > 0) {
			// n[0] = n[1] = n[2] = n[3] = 0;
			for (i = 0, k = 0; i < (uint32_t)overlap_list->length; ++i) {
				overlap_region *r = &overlap_list->list[i];
				w = ha_ov_type(r, rlen);
				// ++n[w];
				// if (((int)n[w] <= max_n_chain) || (r->shared_seed >= s[w] && s[w] >= (asm_opt.k_mer_length<<1))) {
				if (r->shared_seed >= s[w]) {
					if ((uint32_t)k != i) {
						overlap_region t;
						t = overlap_list->list[k];
						overlap_list->list[k] = overlap_list->list[i];
						overlap_list->list[i] = t;
					}
					++k;
				}
			}
			overlap_list->length = k;
		}
	}

	///ks_introsort_or_xs(overlap_list->length, overlap_list->list);
}

void ha_get_new_ul_candidates(ha_abufl_t *ab, int64_t rid, char* rs, int64_t rl, uint64_t mz_w, uint64_t mz_k, const ul_idx_t *uref, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, int max_n_chain, int keep_whole_chain,
						   kvec_t_u8_warp* k_flag, kvec_t_u64_warp* chain_idx, void *ha_flt_tab, ha_pt_t *ha_idx, overlap_region* f_cigar, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t high_occ, void *km)
{
	uint32_t i;
	uint64_t k, l;
	if(high_occ < 1) high_occ = 1;
	// uint32_t high_occ = asm_opt.hom_cov >= 1?asm_opt.hom_cov:1;

	// prepare
    clear_Candidates_list(cl);
    clear_overlap_region_alloc(overlap_list);
	ab->mz.n = 0, ab->n_a = 0;

	// get the list of anchors
	mz2_ha_sketch(rs, rl, mz_w, mz_k, 0, !(asm_opt.flag & HA_F_NO_HPC), &ab->mz, ha_flt_tab, asm_opt.mz_sample_dist, k_flag, dbg_ct, NULL, -1, asm_opt.dp_min_len, -1, sp, asm_opt.mz_rewin, 0, km);
	
	// minimizer of queried read
	if (ab->mz.m > ab->old_mz_m) {
		ab->old_mz_m = ab->mz.m;
		KREALLOC(km, ab->seed, ab->old_mz_m);
	}

	for (i = 0, ab->n_a = 0; i < ab->mz.n; ++i) {
		int n;
		ab->seed[i].a = ha_ptl_get(ha_idx, ab->mz.a[i].x, &n);
		ab->seed[i].n = n;
		ab->n_a += n;
	}
	if (ab->n_a > ab->m_a) {
		ab->m_a = ab->n_a;
		KREALLOC(km, ab->a, ab->m_a);
	}

	for (i = 0, k = 0; i < ab->mz.n; ++i) {
		int j;
		///z is one of the minimizer
		ha_mzl_t *z = &ab->mz.a[i];
		seedl_t *s = &ab->seed[i];
		for (j = 0; j < s->n; ++j) {
			const ha_idxposl_t *y = &s->a[j];
			anchor1_t *an = &ab->a[k++];
			uint8_t rev = z->rev == y->rev? 0 : 1;
			an->other_off = y->pos;
			an->self_off = rev? rl - 1 - (z->pos + 1 - z->span) : z->pos;
			an->cnt = s->n;
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
		KREALLOC(km, cl->list, cl->size);
	}

	for (k = 0; k < ab->n_a; ++k) {
		k_mer_hit *p = &cl->list[k];
		p->readID = ab->a[k].srt >> 33;
		p->strand = ab->a[k].srt >> 32 & 1;
		p->offset = ab->a[k].other_off;
		p->self_offset = ab->a[k].self_off;
		if(ab->a[k].cnt <= high_occ){
			p->cnt = 1;
		}
		else{
			p->cnt = 1 + ((ab->a[k].cnt + (high_occ<<1) - 1)/(high_occ<<1));
			p->cnt = pow(p->cnt, 1.1);
		}
	}
	cl->length = ab->n_a;

	calculate_overlap_region_by_chaining(cl, overlap_list, chain_idx, rid, rl, NULL, uref, bw_thres, keep_whole_chain, f_cigar, km);
	
	#if 0
	if (overlap_list->length > 0) {
		fprintf(stderr, "B\t%ld\t%ld\t%d\n", (long)rid, (long)overlap_list->length, rlen);
		for (int i = 0; i < (int)overlap_list->length; ++i) {
			overlap_region *r = &overlap_list->list[i];
			fprintf(stderr, "C\t%d\t%d\t%d\t%c\t%d\t%ld\t%d\t%d\t%c\t%d\t%d\n", (int)r->x_id, (int)r->x_pos_s, (int)r->x_pos_e, "+-"[r->x_pos_strand],
					(int)r->y_id, (long)Get_READ_LENGTH(R_INF, r->y_id), (int)r->y_pos_s, (int)r->y_pos_e, "+-"[r->y_pos_strand], (int)r->shared_seed, ha_ov_type(r, rlen));
		}
	}
	#endif
	
	if ((int)overlap_list->length > max_n_chain) {
		int32_t w, n[4], s[4];
		n[0] = n[1] = n[2] = n[3] = 0, s[0] = s[1] = s[2] = s[3] = 0;
		ks_introsort_or_ss(overlap_list->length, overlap_list->list);
		for (i = 0; i < (uint32_t)overlap_list->length; ++i) {
			const overlap_region *r = &overlap_list->list[i];
			w = ha_ov_type(r, rl);
			++n[w];
			if ((int)n[w] == max_n_chain) s[w] = r->shared_seed;
		}
		if (s[0] > 0 || s[1] > 0 || s[2] > 0 || s[3] > 0) {
			// n[0] = n[1] = n[2] = n[3] = 0;
			for (i = 0, k = 0; i < (uint32_t)overlap_list->length; ++i) {
				overlap_region *r = &overlap_list->list[i];
				w = ha_ov_type(r, rl);
				// ++n[w];
				// if (((int)n[w] <= max_n_chain) || (r->shared_seed >= s[w] && s[w] >= (asm_opt.k_mer_length<<1))) {
				if (r->shared_seed >= s[w]) {
					if ((uint32_t)k != i) {
						overlap_region t;
						t = overlap_list->list[k];
						overlap_list->list[k] = overlap_list->list[i];
						overlap_list->list[i] = t;
					}
					++k;
				}
			}
			overlap_list->length = k;
		}
	}

	///ks_introsort_or_xs(overlap_list->length, overlap_list->list);
}


void calculate_ug_chaining(Candidates_list* candidates, overlap_region_alloc* overlap_list, kvec_t_u64_warp* chain_idx,
                                          uint64_t readID, ma_utg_v *ua, double band_width_threshold, int add_beg_end, overlap_region* f_cigar, long long mz_occ, double mz_rate)
{
    long long i = 0;
    uint64_t current_ID;
    uint64_t current_stand;

    if (candidates->length == 0)
    {
        return;
    }

    long long sub_region_beg;
    long long sub_region_end;
	long long chain_len;

    clear_fake_cigar(&((*f_cigar).f_cigar));

    i = 0;
    while (i < candidates->length)
    {
        chain_idx->a.n = 0;
        current_ID = candidates->list[i].readID;
        current_stand = candidates->list[i].strand;

        ///reference read
        (*f_cigar).x_id = readID;
        (*f_cigar).x_pos_strand = current_stand;
        ///query read
        (*f_cigar).y_id = current_ID;
        ///here the strand of query is always 0
        (*f_cigar).y_pos_strand = 0;  

        sub_region_beg = i;
        sub_region_end = i;
        i++;

        while (i < candidates->length 
        && 
        current_ID == candidates->list[i].readID
        &&
        current_stand == candidates->list[i].strand)
        {
            sub_region_end = i;
            i++;
        }

        if ((*f_cigar).x_id == (*f_cigar).y_id)
        {
            continue;
        }

		chain_len = chain_DP(candidates->list + sub_region_beg,
				sub_region_end - sub_region_beg + 1, &(candidates->chainDP), f_cigar, band_width_threshold,
				50, ua->a[(*f_cigar).x_id].len, ua->a[(*f_cigar).y_id].len, NULL);
		
		
		// if ((*f_cigar).x_id != (*f_cigar).y_id)
		if ((*f_cigar).x_id != (*f_cigar).y_id && chain_len > mz_occ*mz_rate)
        {
			append_utg_inexact_overlap_region_alloc(overlap_list, f_cigar, ua, add_beg_end, NULL);
        }
    }
}


void ha_get_inter_candidates(ha_abufl_t *ab, uint64_t id, char* r, uint64_t rlen, uint64_t rw, uint64_t rk, uint64_t is_hpc,
        overlap_region_alloc *ol, Candidates_list *cl, double bw_thres, int max_n_chain, int keep_whole_chain,
						    kvec_t_u8_warp* k_flag, kvec_t_u64_warp* chain_idx, void *ha_flt_tab, ha_pt_t *ha_idx, 
                            overlap_region* f_cigar, kvec_t_u64_warp* dbg_ct, st_mt_t *sp)
{
	uint64_t i, k, l;
	// prepare
    clear_Candidates_list(cl);
    clear_overlap_region_alloc(ol);
	ab->mz.n = 0, ab->n_a = 0;

	// get the list of anchors
    mz2_ha_sketch(r, rlen, rw, rk, 0, is_hpc, &ab->mz, ha_flt_tab, asm_opt.mz_sample_dist, k_flag, dbg_ct,   
    NULL, -1, asm_opt.dp_min_len, -1, sp, asm_opt.mz_rewin, 1, NULL);
	
	// minimizer of queried read
	if (ab->mz.m > ab->old_mz_m) {
		ab->old_mz_m = ab->mz.m;
		REALLOC(ab->seed, ab->old_mz_m);
	}
	for (i = 0, ab->n_a = 0; i < ab->mz.n; ++i) {
		int n;
		ab->seed[i].a = ha_ptl_get(ha_idx, ab->mz.a[i].x, &n);
		ab->seed[i].n = n;
		ab->seed[i].cnt = ha_ft_cnt(ha_flt_tab, ab->mz.a[i].x);
		ab->n_a += n;
	}
	if (ab->n_a > ab->m_a) {
		ab->m_a = ab->n_a;
		kroundup64(ab->m_a);
		REALLOC(ab->a, ab->m_a);
	}
	for (i = 0, k = 0; i < ab->mz.n; ++i) {
		int j;
		///z is one of the minimizer
		ha_mzl_t *z = &ab->mz.a[i];
		seedl_t *s = &ab->seed[i];
		for (j = 0; j < s->n; ++j) {
			const ha_idxposl_t *y = &s->a[j];
			anchor1_t *an = &ab->a[k++];
			uint8_t rev = z->rev == y->rev? 0 : 1;
			an->other_off = y->pos;
			an->self_off = rev? rlen - 1 - (z->pos + 1 - z->span) : z->pos;
			an->cnt = s->n;
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
		p->cnt = (ab->a[k].cnt == 1? 1 : 16);
	}
	cl->length = ab->n_a;

	calculate_overlap_region_by_chaining(cl, ol, chain_idx, id, rlen, /**&R_INF**/NULL, NULL, bw_thres, keep_whole_chain, f_cigar, NULL);

	#if 0
	if (ol->length > 0) {
		fprintf(stderr, "B\t%ld\t%ld\t%d\n", (long)rid, (long)ol->length, rlen);
		for (int i = 0; i < (int)ol->length; ++i) {
			overlap_region *r = &ol->list[i];
			fprintf(stderr, "C\t%d\t%d\t%d\t%c\t%d\t%ld\t%d\t%d\t%c\t%d\t%d\n", (int)r->x_id, (int)r->x_pos_s, (int)r->x_pos_e, "+-"[r->x_pos_strand],
					(int)r->y_id, (long)Get_READ_LENGTH(R_INF, r->y_id), (int)r->y_pos_s, (int)r->y_pos_e, "+-"[r->y_pos_strand], (int)r->shared_seed, ha_ov_type(r, rlen));
		}
	}
	#endif
	
	if ((int)ol->length > max_n_chain) {
		int32_t w, n[4], s[4];
		n[0] = n[1] = n[2] = n[3] = 0, s[0] = s[1] = s[2] = s[3] = 0;
		ks_introsort_or_ss(ol->length, ol->list);
		for (i = 0; i < (uint32_t)ol->length; ++i) {
			const overlap_region *r = &ol->list[i];
			w = ha_ov_type(r, rlen);
			++n[w];
			if ((int)n[w] == max_n_chain) s[w] = r->shared_seed;
		}
		if (s[0] > 0 || s[1] > 0 || s[2] > 0 || s[3] > 0) {
			// n[0] = n[1] = n[2] = n[3] = 0;
			for (i = 0, k = 0; i < (uint32_t)ol->length; ++i) {
				overlap_region *r = &ol->list[i];
				w = ha_ov_type(r, rlen);
				// ++n[w];
				// if (((int)n[w] <= max_n_chain) || (r->shared_seed >= s[w] && s[w] >= (asm_opt.k_mer_length<<1))) {
				if (r->shared_seed >= s[w]) {
					if ((uint32_t)k != i) {
						overlap_region t;
						t = ol->list[k];
						ol->list[k] = ol->list[i];
						ol->list[i] = t;
					}
					++k;
				}
			}
			ol->length = k;
		}
	}

	///ks_introsort_or_xs(overlap_list->length, overlap_list->list);
}

void ha_get_ug_candidates(ha_abuf_t *ab, int64_t rid, ma_utg_t *u, ma_utg_v *ua, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, int max_n_chain, int keep_whole_chain, kvec_t_u8_warp* k_flag,
kvec_t_u64_warp* chain_idx, void *ha_flt_tab, ha_pt_t *ha_idx, overlap_region* f_cigar, kvec_t_u64_warp* dbg_ct, double chain_match_rate)
{
    uint32_t i;
    uint64_t k, l;

    // prepare
    clear_Candidates_list(cl);
    clear_overlap_region_alloc(overlap_list);
    ab->mz.n = 0, ab->n_a = 0;

    // get the list of anchors
    //should use the new version...
    ///ha_sketch_query(u->s, u->len, asm_opt.mz_win, asm_opt.k_mer_length, 0, !(asm_opt.flag & HA_F_NO_HPC), &ab->mz, ha_flt_tab, k_flag, dbg_ct);
    // minimizer of queried read
    if (ab->mz.m > ab->old_mz_m) {
        ab->old_mz_m = ab->mz.m;
        REALLOC(ab->seed, ab->old_mz_m);
    }
    for (i = 0, ab->n_a = 0; i < ab->mz.n; ++i) {
        int n;
        ab->seed[i].a = ha_pt_get(ha_idx, ab->mz.a[i].x, &n);
        ab->seed[i].n = n;
        ab->n_a += n;
    }
    if (ab->n_a > ab->m_a) {
        ab->m_a = ab->n_a;
        kroundup64(ab->m_a);
        REALLOC(ab->a, ab->m_a);
    }
    for (i = 0, k = 0; i < ab->mz.n; ++i) {
        int j;
        ///z is one of the minimizer
        ha_mz1_t *z = &ab->mz.a[i];
        seed1_t *s = &ab->seed[i];
        for (j = 0; j < s->n; ++j) {
            const ha_idxpos_t *y = &s->a[j];
            anchor1_t *an = &ab->a[k++];
            uint8_t rev = z->rev == y->rev? 0 : 1;
            an->other_off = y->pos;
            an->self_off = rev? u->len - 1 - (z->pos + 1 - z->span) : z->pos;
            an->cnt = 1;
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
        p->cnt = 1;
    }
    cl->length = ab->n_a;

    calculate_ug_chaining(cl, overlap_list, chain_idx, rid, ua, bw_thres, keep_whole_chain, f_cigar, ab->mz.n, chain_match_rate);

    #if 0
    if (overlap_list->length > 0) {
        fprintf(stderr, "B\t%ld\t%ld\t%d\n", (long)rid, (long)overlap_list->length, rlen);
        for (int i = 0; i < (int)overlap_list->length; ++i) {
            overlap_region *r = &overlap_list->list[i];
            fprintf(stderr, "C\t%d\t%d\t%d\t%c\t%d\t%ld\t%d\t%d\t%c\t%d\t%d\n", (int)r->x_id, (int)r->x_pos_s, (int)r->x_pos_e, "+-"[r->x_pos_strand],
                    (int)r->y_id, (long)Get_READ_LENGTH(R_INF, r->y_id), (int)r->y_pos_s, (int)r->y_pos_e, "+-"[r->y_pos_strand], (int)r->shared_seed, ha_ov_type(r, rlen));
        }
    }
    #endif

    if ((int)overlap_list->length > max_n_chain) {
        int32_t w, n[4], s[4];
        n[0] = n[1] = n[2] = n[3] = 0, s[0] = s[1] = s[2] = s[3] = 0;
        ks_introsort_or_ss(overlap_list->length, overlap_list->list);
        for (i = 0; i < (uint32_t)overlap_list->length; ++i) {
            const overlap_region *r = &overlap_list->list[i];
            w = ha_ov_type(r, u->len);
            ++n[w];
            if ((int)n[w] == max_n_chain) s[w] = r->shared_seed;
        }
        if (s[0] > 0 || s[1] > 0 || s[2] > 0 || s[3] > 0) {
            for (i = 0, k = 0; i < (uint32_t)overlap_list->length; ++i) {
                overlap_region *r = &overlap_list->list[i];
                w = ha_ov_type(r, u->len);
                if (r->shared_seed >= s[w]) {
                    if ((uint32_t)k != i) {
                        overlap_region t;
                        t = overlap_list->list[k];
                        overlap_list->list[k] = overlap_list->list[i];
                        overlap_list->list[i] = t;
                    }
                    ++k;
                }
            }
            overlap_list->length = k;
        }
    }

    ///ks_introsort_or_xs(overlap_list->length, overlap_list->list);
}


void lable_matched_ovlp(overlap_region_alloc* overlap_list, ma_hit_t_alloc* paf)
{
	uint64_t j = 0, inner_j = 0;
	while (j < overlap_list->length && inner_j < paf->length)
    {
        if(overlap_list->list[j].y_id < paf->buffer[inner_j].tn)
        {
            j++;
        }
        else if(overlap_list->list[j].y_id > paf->buffer[inner_j].tn)
        {
            inner_j++;
        }
        else
        {
            if(overlap_list->list[j].y_pos_strand == paf->buffer[inner_j].rev)
            {
				overlap_list->list[j].is_match = 1;
            }
            j++;
            inner_j++;
        }
    }
}


void ha_get_candidates_interface(ha_abuf_t *ab, int64_t rid, UC_Read *ucr, overlap_region_alloc *overlap_list, overlap_region_alloc *overlap_list_hp, Candidates_list *cl, double bw_thres, 
								 int max_n_chain, int keep_whole_chain, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* chain_idx, ma_hit_t_alloc* paf, ma_hit_t_alloc* rev_paf, overlap_region* f_cigar,
								 kvec_t_u64_warp* dbg_ct, st_mt_t *sp)
{
	extern void *ha_flt_tab;
	extern ha_pt_t *ha_idx;
	extern void *ha_flt_tab_hp;
	extern ha_pt_t *ha_idx_hp;

	ha_get_new_candidates(ab, rid, ucr, overlap_list, cl, bw_thres, max_n_chain, keep_whole_chain, k_flag, chain_idx, ha_flt_tab, ha_idx, f_cigar, dbg_ct, sp);

	if(ha_idx_hp)
	{
		uint32_t i, k, y_id, overlapLen, max_i;
		int shared_seed;
		overlap_region t;
		overlap_region_sort_y_id(overlap_list->list, overlap_list->length);
		ma_hit_sort_tn(paf->buffer, paf->length);
		ma_hit_sort_tn(rev_paf->buffer, rev_paf->length);
		lable_matched_ovlp(overlap_list, paf);
		lable_matched_ovlp(overlap_list, rev_paf);

		for (i = 0, k = 0; i < overlap_list->length; ++i) 
		{
			if(overlap_list->list[i].is_match == 1)
			{
				if(k != i)
				{
					t = overlap_list->list[k];
					overlap_list->list[k] = overlap_list->list[i];
					overlap_list->list[i] = t;
					overlap_list->list[k].is_match = 0;
				}
				k++;
			}
		}
		overlap_list->length = k;


		ha_get_new_candidates(ab, rid, ucr, overlap_list_hp, cl, bw_thres, max_n_chain, keep_whole_chain, k_flag, chain_idx, ha_flt_tab_hp, ha_idx_hp, f_cigar, dbg_ct, sp);	
		
		if(overlap_list->length + overlap_list_hp->length > overlap_list->size)
		{
			overlap_list->list = (overlap_region*)realloc(overlap_list->list, 
				sizeof(overlap_region)*(overlap_list->length + overlap_list_hp->length));
			memset(overlap_list->list + overlap_list->size, 0, sizeof(overlap_region)*
				(overlap_list->length + overlap_list_hp->length - overlap_list->size));
			overlap_list->size = overlap_list->length + overlap_list_hp->length;
		}
		
		for (i = 0, k = overlap_list->length; i < overlap_list_hp->length; i++, k++)
		{
			t = overlap_list->list[k];
			overlap_list->list[k] = overlap_list_hp->list[i];
			overlap_list_hp->list[i] = t;
		}
		overlap_list->length = k;
		
		overlap_region_sort_y_id(overlap_list->list, overlap_list->length);

		i = k = 0;
		while (i < overlap_list->length)
		{
			y_id = overlap_list->list[i].y_id;
			shared_seed = overlap_list->list[i].shared_seed;
			overlapLen = overlap_list->list[i].overlapLen;
			max_i = i;
			i++;
			while (i < overlap_list->length && overlap_list->list[i].y_id == y_id)
			{
				if((overlap_list->list[i].shared_seed > shared_seed) || 
				  ((overlap_list->list[i].shared_seed == shared_seed) && (overlap_list->list[i].overlapLen <= overlapLen)))
				{
					y_id = overlap_list->list[i].y_id;
					shared_seed = overlap_list->list[i].shared_seed;
					overlapLen = overlap_list->list[i].overlapLen;
					max_i = i;
				}
				i++;
			}

			if(k != max_i)
			{
				t = overlap_list->list[k];
				overlap_list->list[k] = overlap_list->list[max_i];
				overlap_list->list[max_i] = t;
			}
			k++;
		}

		overlap_list->length = k;
	}
	
	ks_introsort_or_xs(overlap_list->length, overlap_list->list);
}


void ha_get_ul_candidates_interface(ha_abufl_t *ab, int64_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k, const ul_idx_t *uref, overlap_region_alloc *overlap_list, overlap_region_alloc *overlap_list_hp, Candidates_list *cl, double bw_thres, 
								 int max_n_chain, int keep_whole_chain, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* chain_idx, overlap_region* f_cigar, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t high_occ, void *km)
{
	extern void *ha_flt_tab;
	extern ha_pt_t *ha_idx;

	ha_get_new_ul_candidates(ab, rid, rs, rl, mz_w, mz_k, uref, overlap_list, cl, bw_thres, max_n_chain, keep_whole_chain, k_flag, chain_idx, ha_flt_tab, ha_idx, f_cigar, dbg_ct, sp, high_occ, km);
	if(km) {
		ha_abufl_free_buf(km, ab, 1); 
		destory_Candidates_list_buf(km, cl, 1); 
	}
	ks_introsort_or_xs(overlap_list->length, overlap_list->list);
}



void ha_sort_list_by_anchor(overlap_region_alloc *overlap_list)
{
	ks_introsort_or_xs(overlap_list->length, overlap_list->list);
}


void minimizers_gen(ha_abufl_t *ab, char* rs, int64_t rl, uint64_t mz_w, uint64_t mz_k, Candidates_list *cl, kvec_t_u8_warp* k_flag, 
void *ha_flt_tab, ha_pt_t *ha_idx, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ)
{
	// fprintf(stderr, "+[M::%s]\n", __func__);
	uint64_t i, k, l, max_cnt = UINT32_MAX, min_cnt = 0; int n, j; ha_mzl_t *z; seedl_t *s; 
	if(high_occ) {
		max_cnt = (*high_occ);
		if(max_cnt < 2) max_cnt = 2;
	}
	if(low_occ) {
		min_cnt = (*low_occ);
		if(min_cnt < 2) min_cnt = 2;
	}
	clear_Candidates_list(cl); ab->mz.n = 0, ab->n_a = 0;
	
	// get the list of anchors
	mz2_ha_sketch(rs, rl, mz_w, mz_k, 0, !(asm_opt.flag & HA_F_NO_HPC), &ab->mz, ha_flt_tab, asm_opt.mz_sample_dist, k_flag, dbg_ct, NULL, -1, asm_opt.dp_min_len, -1, sp, asm_opt.mz_rewin, 0, NULL);

	// minimizer of queried read
	if (ab->mz.m > ab->old_mz_m) {
		ab->old_mz_m = ab->mz.m;
		REALLOC(ab->seed, ab->old_mz_m);
	}

	for (i = 0, ab->n_a = 0; i < ab->mz.n; ++i) {
		ab->seed[i].a = ha_ptl_get(ha_idx, ab->mz.a[i].x, &n);
		ab->seed[i].n = n;
		ab->n_a += n;
	}

	if (ab->n_a > ab->m_a) {
		ab->m_a = ab->n_a;
		REALLOC(ab->a, ab->m_a);
	}

	for (i = 0, k = 0; i < ab->mz.n; ++i) {
		///z is one of the minimizer
		z = &ab->mz.a[i]; s = &ab->seed[i];
		for (j = 0; j < s->n; ++j) {
			const ha_idxposl_t *y = &s->a[j];
			anchor1_t *an = &ab->a[k++];
			uint8_t rev = z->rev == y->rev? 0 : 1;
			an->other_off = y->pos;
			an->self_off = rev? rl - 1 - (z->pos + 1 - z->span) : z->pos;
			///an->cnt: cnt<<8|span
			an->cnt = s->n; if(an->cnt > ((uint32_t)(0xffffffu))) an->cnt = 0xffffffu;
			an->cnt <<= 8; an->cnt |= ((z->span <= ((uint32_t)(0xffu)))?z->span:((uint32_t)(0xffu)));
			an->srt = (uint64_t)y->rid<<33 | (uint64_t)rev<<32 | an->other_off;
		}
	}

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
		if(((ab->a[k].cnt>>8) < max_cnt) && ((ab->a[k].cnt>>8) > min_cnt)){
			p->cnt = 1;
		} else if((ab->a[k].cnt>>8) <= min_cnt) {
			p->cnt = 2;
		} else{
			p->cnt = 1 + (((ab->a[k].cnt>>8) + (max_cnt<<1) - 1)/(max_cnt<<1));
			p->cnt = pow(p->cnt, 1.1);
		}
		if(p->cnt > ((uint32_t)(0xffffffu))) p->cnt = 0xffffffu;
		p->cnt <<= 8; p->cnt |= (((uint32_t)(0xffu))&(ab->a[k].cnt));
	}
	cl->length = ab->n_a;
}


void minimizers_qgen(ha_abufl_t *ab, char* rs, int64_t rl, uint64_t mz_w, uint64_t mz_k, Candidates_list *cl, kvec_t_u8_warp* k_flag, 
void *ha_flt_tab, ha_pt_t *ha_idx, All_reads* rdb, const ul_idx_t *udb, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, 
uint32_t *low_occ)
{
	// fprintf(stderr, "+[M::%s]\n", __func__);
	uint64_t i, k, l, max_cnt = UINT32_MAX, min_cnt = 0; int n, j; ha_mzl_t *z; seedl_t *s; 
	if(high_occ) {
		max_cnt = (*high_occ);
		if(max_cnt < 2) max_cnt = 2;
	}
	if(low_occ) {
		min_cnt = (*low_occ);
		if(min_cnt < 2) min_cnt = 2;
	}
	clear_Candidates_list(cl); ab->mz.n = 0, ab->n_a = 0;
	
	// get the list of anchors
	mz2_ha_sketch(rs, rl, mz_w, mz_k, 0, !(asm_opt.flag & HA_F_NO_HPC), &ab->mz, ha_flt_tab, asm_opt.mz_sample_dist, k_flag, dbg_ct, NULL, -1, asm_opt.dp_min_len, -1, sp, asm_opt.mz_rewin, 0, NULL);

	// minimizer of queried read
	if (ab->mz.m > ab->old_mz_m) {
		ab->old_mz_m = ab->mz.m;
		REALLOC(ab->seed, ab->old_mz_m);
	}

	for (i = 0, ab->n_a = 0; i < ab->mz.n; ++i) {
		ab->seed[i].a = ha_ptl_get(ha_idx, ab->mz.a[i].x, &n);
		ab->seed[i].n = n;
		ab->n_a += n;
	}

	if (ab->n_a > ab->m_a) {
		ab->m_a = ab->n_a;
		REALLOC(ab->a, ab->m_a);
	}

	for (i = 0, k = 0; i < ab->mz.n; ++i) {
		///z is one of the minimizer
		z = &ab->mz.a[i]; s = &ab->seed[i];
		for (j = 0; j < s->n; ++j) {
			const ha_idxposl_t *y = &s->a[j];
			anchor1_t *an = &ab->a[k++];
			uint8_t rev = z->rev == y->rev? 0 : 1;
			an->other_off = rev?((uint32_t)-1)-1-(y->pos+1-y->span):y->pos;
			an->self_off = z->pos;
			///an->cnt: cnt<<8|span
			an->cnt = s->n; if(an->cnt > ((uint32_t)(0xffffffu))) an->cnt = 0xffffffu;
			an->cnt <<= 8; an->cnt |= ((z->span <= ((uint32_t)(0xffu)))?z->span:((uint32_t)(0xffu)));
			an->srt = (uint64_t)y->rid<<33 | (uint64_t)rev<<32 | an->self_off;
		}
	}

	// copy over to _cl_
	if (ab->m_a >= (uint64_t)cl->size) {
		cl->size = ab->m_a;
		REALLOC(cl->list, cl->size);
	}

	k_mer_hit *p; uint64_t tid = (uint64_t)-1, tl = (uint64_t)-1;
	radix_sort_ha_an1(ab->a, ab->a + ab->n_a);
	for (k = 1, l = 0; k <= ab->n_a; ++k) {
		if (k == ab->n_a || ab->a[k].srt != ab->a[l].srt) {
			if (k-l>1) radix_sort_ha_an3(ab->a+l, ab->a+k);
			if((ab->a[l].srt>>33)!=tid) {
				tid = ab->a[l].srt>>33;
				tl = rdb?Get_READ_LENGTH((*rdb), tid):udb->ug->u.a[tid].len;
			} 
			for (i = l; i < k; i++) {
				p = &cl->list[i];
				p->readID = ab->a[i].srt>>33;
				p->strand = (ab->a[i].srt>>32)&1;
				if(!(p->strand)) {
					p->offset = ab->a[i].other_off;
				} else {
					p->offset = ((uint32_t)-1)-ab->a[i].other_off;
					p->offset = tl-p->offset;
				}
				p->self_offset = ab->a[i].self_off;
				if(((ab->a[i].cnt>>8) < max_cnt) && ((ab->a[i].cnt>>8) > min_cnt)){
					p->cnt = 1;
				} else if((ab->a[i].cnt>>8) <= min_cnt) {
					p->cnt = 2;
				} else{
					p->cnt = 1 + (((ab->a[i].cnt>>8) + (max_cnt<<1) - 1)/(max_cnt<<1));
					p->cnt = pow(p->cnt, 1.1);
				}
				if(p->cnt > ((uint32_t)(0xffffffu))) p->cnt = 0xffffffu;
				p->cnt <<= 8; p->cnt |= (((uint32_t)(0xffu))&(ab->a[i].cnt));
			}
			l = k;
		}
	}
	cl->length = ab->n_a;
}


void minimizers_qgen0(ha_abuf_t *ab, char* rs, int64_t rl, uint64_t mz_w, uint64_t mz_k, Candidates_list *cl, kvec_t_u8_warp* k_flag, 
void *ha_flt_tab, ha_pt_t *ha_idx, All_reads* rdb, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ)
{
	// fprintf(stderr, "+[M::%s]\n", __func__);
	uint64_t i, k, l, max_cnt = UINT32_MAX, min_cnt = 0; int n, j; ha_mz1_t *z; seed1_t *s; 
	if(high_occ) {
		max_cnt = (*high_occ);
		if(max_cnt < 2) max_cnt = 2;
	}
	if(low_occ) {
		min_cnt = (*low_occ);
		if(min_cnt < 2) min_cnt = 2;
	}
	clear_Candidates_list(cl); ab->mz.n = 0, ab->n_a = 0;
	
	// get the list of anchors
	mz1_ha_sketch(rs, rl, mz_w, mz_k, 0, !(asm_opt.flag & HA_F_NO_HPC), &ab->mz, ha_flt_tab, asm_opt.mz_sample_dist, k_flag, dbg_ct, NULL, -1, asm_opt.dp_min_len, -1, sp, asm_opt.mz_rewin, 0, NULL);

	// minimizer of queried read
	if (ab->mz.m > ab->old_mz_m) {
		ab->old_mz_m = ab->mz.m;
		REALLOC(ab->seed, ab->old_mz_m);
	}

	for (i = 0, ab->n_a = 0; i < ab->mz.n; ++i) {

		ab->seed[i].a = ha_pt_get(ha_idx, ab->mz.a[i].x, &n);
		ab->seed[i].n = n;
		ab->n_a += n;
	}

	if (ab->n_a > ab->m_a) {
		ab->m_a = ab->n_a;
		REALLOC(ab->a, ab->m_a);
	}

	for (i = 0, k = 0; i < ab->mz.n; ++i) {
		///z is one of the minimizer
		z = &ab->mz.a[i]; s = &ab->seed[i];
		for (j = 0; j < s->n; ++j) {
			const ha_idxpos_t *y = &s->a[j];
			anchor1_t *an = &ab->a[k++];
			uint8_t rev = z->rev == y->rev? 0 : 1;
			an->other_off = rev?((uint32_t)-1)-1-(y->pos+1-y->span):y->pos;
			an->self_off = z->pos;
			///an->cnt: cnt<<8|span
			an->cnt = s->n; if(an->cnt > ((uint32_t)(0xffffffu))) an->cnt = 0xffffffu;
			an->cnt <<= 8; an->cnt |= ((z->span <= ((uint32_t)(0xffu)))?z->span:((uint32_t)(0xffu)));
			an->srt = (uint64_t)y->rid<<33 | (uint64_t)rev<<32 | an->self_off;
		}
	}

	// copy over to _cl_
	if (ab->m_a >= (uint64_t)cl->size) {
		cl->size = ab->m_a;
		REALLOC(cl->list, cl->size);
	}

	k_mer_hit *p; uint64_t tid = (uint64_t)-1, tl = (uint64_t)-1;
	radix_sort_ha_an1(ab->a, ab->a + ab->n_a);
	for (k = 1, l = 0; k <= ab->n_a; ++k) {
		if (k == ab->n_a || ab->a[k].srt != ab->a[l].srt) {
			if (k-l>1) radix_sort_ha_an3(ab->a+l, ab->a+k);
			if((ab->a[l].srt>>33)!=tid) {
				tid = ab->a[l].srt>>33; 
				tl = Get_READ_LENGTH((*rdb), tid);
				// tl = rdb?Get_READ_LENGTH((*rdb), tid):udb->ug->u.a[tid].len;
			} 
			for (i = l; i < k; i++) {
				p = &cl->list[i];
				p->readID = ab->a[i].srt>>33;
				p->strand = (ab->a[i].srt>>32)&1;
				if(!(p->strand)) {
					p->offset = ab->a[i].other_off;
				} else {
					p->offset = ((uint32_t)-1)-ab->a[i].other_off;
					p->offset = tl-p->offset;
				}
				p->self_offset = ab->a[i].self_off;
				if(((ab->a[i].cnt>>8) < max_cnt) && ((ab->a[i].cnt>>8) > min_cnt)){
					p->cnt = 1;
				} else if((ab->a[i].cnt>>8) <= min_cnt) {
					p->cnt = 2;
				} else{
					p->cnt = 1 + (((ab->a[i].cnt>>8) + (max_cnt<<1) - 1)/(max_cnt<<1));
					p->cnt = pow(p->cnt, 1.1);
				}
				if(p->cnt > ((uint32_t)(0xffffffu))) p->cnt = 0xffffffu;
				p->cnt <<= 8; p->cnt |= (((uint32_t)(0xffu))&(ab->a[i].cnt));
			}
			l = k;
		}
	}
	cl->length = ab->n_a;
}

void minimizers_qgen0_amz(ha_abuf_t *ab, char* rs, int64_t rl, uint64_t mz_w, uint64_t mz_k, Candidates_list *cl, kvec_t_u8_warp* k_flag, 
void *ha_flt_tab, ha_pt_t *ha_idx, All_reads* rdb, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ)
{
	// fprintf(stderr, "+[M::%s]\n", __func__);
	uint64_t i, k, l, max_cnt = UINT32_MAX, min_cnt = 0; int n, j; ha_mz1_t *z; seed1_t *s; 
	if(high_occ) {
		max_cnt = (*high_occ);
		if(max_cnt < 2) max_cnt = 2;
	}
	if(low_occ) {
		min_cnt = (*low_occ);
		if(min_cnt < 2) min_cnt = 2;
	}
	clear_Candidates_list(cl); ab->mz.n = 0, ab->n_a = 0;
	
	// get the list of anchors
	mz1_ha_sketch(rs, rl, mz_w, mz_k, 0, !(asm_opt.flag & HA_F_NO_HPC), &ab->mz, ha_flt_tab, asm_opt.mz_sample_dist, k_flag, dbg_ct, NULL, -1, asm_opt.dp_min_len, -1, sp, asm_opt.mz_rewin, 0, NULL);

	// minimizer of queried read
	if (ab->mz.m > ab->old_mz_m) {
		ab->old_mz_m = ab->mz.m;
		REALLOC(ab->seed, ab->old_mz_m);
	}

	for (i = 0, ab->n_a = 0; i < ab->mz.n; ++i) {

		ab->seed[i].a = ha_pt_get(ha_idx, ab->mz.a[i].x, &n);
		ab->seed[i].n = n;
		ab->n_a += n;
	}

	if (ab->n_a > ab->m_a) {
		ab->m_a = ab->n_a;
		REALLOC(ab->a, ab->m_a);
	}

	for (i = 0, k = 0; i < ab->mz.n; ++i) {
		///z is one of the minimizer
		z = &ab->mz.a[i]; s = &ab->seed[i];
		for (j = 0; j < s->n; ++j) {
			const ha_idxpos_t *y = &s->a[j];
			anchor1_t *an = &ab->a[k++];
			uint8_t rev = z->rev == y->rev? 0 : 1;
			an->other_off = rev?((uint32_t)-1)-1-(y->pos+1-y->span):y->pos;
			an->self_off = z->pos;
			///an->cnt: cnt<<8|span
			an->cnt = s->n; if(an->cnt > ((uint32_t)(0xffffffu))) an->cnt = 0xffffffu;
			an->cnt <<= 8; an->cnt |= ((z->span <= ((uint32_t)(0xffu)))?z->span:((uint32_t)(0xffu)));
			an->srt = (uint64_t)y->rid<<33 | (uint64_t)rev<<32 | an->self_off;
		}
	}

	// copy over to _cl_
	if (ab->m_a >= (uint64_t)cl->size) {
		cl->size = ab->m_a;
		REALLOC(cl->list, cl->size);
	}

	k_mer_hit *p; uint64_t tid = (uint64_t)-1, tl = (uint64_t)-1;
	radix_sort_ha_an1(ab->a, ab->a + ab->n_a);
	for (k = 1, l = 0; k <= ab->n_a; ++k) {
		if (k == ab->n_a || ab->a[k].srt != ab->a[l].srt) {
			if (k-l>1) radix_sort_ha_an3(ab->a+l, ab->a+k);
			if((ab->a[l].srt>>33)!=tid) {
				tid = ab->a[l].srt>>33; 
				tl = Get_READ_LENGTH((*rdb), tid);
				// tl = rdb?Get_READ_LENGTH((*rdb), tid):udb->ug->u.a[tid].len;
			} 
			for (i = l; i < k; i++) {
				p = &cl->list[i];
				p->readID = ab->a[i].srt>>33;
				p->strand = (ab->a[i].srt>>32)&1;
				if(!(p->strand)) {
					p->offset = ab->a[i].other_off;
				} else {
					p->offset = ((uint32_t)-1)-ab->a[i].other_off;
					p->offset = tl-p->offset;
				}
				p->self_offset = ab->a[i].self_off;
				if(((ab->a[i].cnt>>8) < max_cnt) && ((ab->a[i].cnt>>8) > min_cnt)){
					p->cnt = 1;
				} else if((ab->a[i].cnt>>8) <= min_cnt) {
					p->cnt = 2;
				} else{
					p->cnt = 1 + (((ab->a[i].cnt>>8) + (max_cnt<<1) - 1)/(max_cnt<<1));
					p->cnt = pow(p->cnt, 1.1);
				}
				if(p->cnt > ((uint32_t)(0xffffffu))) p->cnt = 0xffffffu;
				p->cnt <<= 8; p->cnt |= (((uint32_t)(0xffu))&(ab->a[i].cnt));
			}
			l = k;
		}
	}
	cl->length = ab->n_a;
}


uint64_t lchain_qgen_mcopy_fast_re0(ha_abuf_t *ab, ha_pt_t *ha_idx, ha_mz1_t *ra, uint64_t rn, ha_mz1_t *qa, uint64_t qn, uint64_t qid, Candidates_list *cl, uint32_t *high_occ, uint32_t *low_occ)
{
	// fprintf(stderr, "+[M::%s]\n", __func__);
	uint64_t i, k, l, max_cnt = UINT32_MAX, min_cnt = 0, ri, qi, sn; anchor1_t *an; k_mer_hit *p;
	if(high_occ) {
		max_cnt = (*high_occ);
		if(max_cnt < 2) max_cnt = 2;
	}
	if(low_occ) {
		min_cnt = (*low_occ);
		if(min_cnt < 2) min_cnt = 2;
	}
	// clear_Candidates_list(cl); 

	///first try
	for (k = 1, l = i = ab->n_a = 0; k <= rn; ++k) {
		if (k == rn || ra[k].x != ra[l].x) {
			for (; i < qn && qa[i].x < ra[l].x; i++);
			if(i < qn && qa[i].x == ra[l].x) {
				sn = 0;
				if(ab->n_a < ab->m_a) sn = ab->seed[l].n;///ha_pt_cnt(ha_idx, ra[l].x);
							
				for (qi = i; qi < qn && qa[qi].x == ra[l].x; qi++) {
					for (ri = l; ri < k; ri++) {
						if(qa[qi].rev != ra[ri].rev) continue;
						if(ab->n_a < ab->m_a) {
							an = &(ab->a[ab->n_a++]); 

							an->other_off = qa[qi].pos;
							an->self_off = ra[ri].pos;
							///an->cnt: cnt<<8|span
							an->cnt = sn; if(an->cnt > ((uint32_t)(0xffffffu))) an->cnt = 0xffffffu;
							an->cnt <<= 8; an->cnt |= ((ra[ri].span <= ((uint32_t)(0xffu)))?ra[ri].span:((uint32_t)(0xffu)));
							an->srt = (((uint64_t)(an->self_off))<<32)|((uint64_t)(an->other_off));
						} else {
							ab->n_a++;
						}
					}
				}
			}
			l = k;
		}
	}

	
	if (ab->n_a > ab->m_a) {
		ab->m_a = ab->n_a; REALLOC(ab->a, ab->m_a);

		for (k = 1, l = i = ab->n_a = 0; k <= rn; ++k) {
			if (k == rn || ra[k].x != ra[l].x) {
				for (; i < qn && qa[i].x < ra[l].x; i++);
				if(i < qn && qa[i].x == ra[l].x) {
					sn = ab->seed[l].n;///ha_pt_cnt(ha_idx, ra[l].x);
					
					for (qi = i; qi < qn && qa[qi].x == ra[l].x; qi++) {
						for (ri = l; ri < k; ri++) {
							if(qa[qi].rev != ra[ri].rev) continue;

							an = &(ab->a[ab->n_a++]); 

							an->other_off = qa[qi].pos;
							an->self_off = ra[ri].pos;
							///an->cnt: cnt<<8|span
							an->cnt = sn; if(an->cnt > ((uint32_t)(0xffffffu))) an->cnt = 0xffffffu;
							an->cnt <<= 8; an->cnt |= ((ra[ri].span <= ((uint32_t)(0xffu)))?ra[ri].span:((uint32_t)(0xffu)));
							an->srt = (((uint64_t)(an->self_off))<<32)|((uint64_t)(an->other_off));
						}
					}
				}
				l = k;
			}
		}
	}

	// copy over to _cl_
	sn = ab->n_a + cl->length;
	if (sn > (uint64_t)cl->size) {
		cl->size = sn;
		REALLOC(cl->list, cl->size);
	}

	radix_sort_ha_an1(ab->a, ab->a + ab->n_a);
	for (i = 0; i < ab->n_a; i++) {
		p = &cl->list[cl->length++];
		p->readID = qid;
		p->strand = 0;
		p->offset = ab->a[i].other_off;
		p->self_offset = ab->a[i].self_off;

		if(((ab->a[i].cnt>>8) < max_cnt) && ((ab->a[i].cnt>>8) > min_cnt)){
			p->cnt = 1;
		} else if((ab->a[i].cnt>>8) <= min_cnt) {
			p->cnt = 2;
		} else{
			p->cnt = 1 + (((ab->a[i].cnt>>8) + (max_cnt<<1) - 1)/(max_cnt<<1));
			p->cnt = pow(p->cnt, 1.1);
		}
		if(p->cnt > ((uint32_t)(0xffffffu))) p->cnt = 0xffffffu;
		p->cnt <<= 8; p->cnt |= (((uint32_t)(0xffu))&(ab->a[i].cnt));
	}
	// cl->length = ab->n_a;
	return ab->n_a;
}

void gen_pair_chain(ha_abufl_t *ab, uint64_t rid, st_mt_t *tid, uint64_t tid_n, ha_mzl_t *in, uint64_t in_n, ha_mzl_t *idx, int64_t idx_n, uint64_t mzl_cutoff)
{
	if(!tid_n) return;
	uint64_t i, l, k, n, m, rev, x, tn, kn; ha_mzl_t *z, *p, *y; int64_t zi, ns, ne; anchor1_t *an;
	for (i = tn = 0; i < in_n; ++i) {
		///z is one of the minimizer
		z = &in[i]; p = &(idx[z->x]); n = 0;
		assert(z->pos == p->pos && z->rid == p->rid && z->span == p->span && z->rev == p->rev);
		for (zi = z->x+1; zi < idx_n && idx[zi].x == p->x && n < mzl_cutoff; zi++) {
			if(idx[zi].rid == rid) continue; 
			x = idx[zi].rid; x <<= 1; x |= ((uint64_t)((z->rev == idx[zi].rev)?0:1));
			for (m = 0; m < tid_n; m++) {
				if(tid->a[m] == x) {
					n++; break;
				}
			}
		}
		for (zi = ((int64_t)z->x)-1; zi >= 0 && idx[zi].x == p->x && n < mzl_cutoff; zi--) {
			if(idx[zi].rid == rid) continue; 
			x = idx[zi].rid; x <<= 1; x |= ((uint64_t)((z->rev == idx[zi].rev)?0:1));
			for (m = 0; m < tid_n; m++) {
				if(tid->a[m] == x) {
					n++; break;
				}
			}
		}
		if((!n) || (n >= mzl_cutoff)) continue;
		tn += n;
	}
	if(!tn) return;
	ab->n_a += tn;
	if (ab->n_a > ab->m_a) {
		ab->m_a = ab->n_a;
		REALLOC(ab->a, ab->m_a);
	}


	for (i = 0, k = ab->n_a - tn; i < in_n; ++i) {
		///z is one of the minimizer
		z = &in[i]; p = &(idx[z->x]); n = 0; ns = 0; ne = -1;
		for (zi = z->x+1; zi < idx_n && idx[zi].x == p->x && n < mzl_cutoff; zi++) {
			if(idx[zi].rid == rid) continue; 
			x = idx[zi].rid; x <<= 1; x |= ((uint64_t)((z->rev == idx[zi].rev)?0:1));
			for (m = 0; m < tid_n; m++) {
				if(tid->a[m] == x) {
					n++; break;
				}
			}
		}
		ne = zi;
		for (zi = ((int64_t)z->x)-1; zi >= 0 && idx[zi].x == p->x && n < mzl_cutoff; zi--) {
			if(idx[zi].rid == rid) continue; 
			x = idx[zi].rid; x <<= 1; x |= ((uint64_t)((z->rev == idx[zi].rev)?0:1));
			for (m = 0; m < tid_n; m++) {
				if(tid->a[m] == x) {
					n++; break;
				}
			}
		}
		ns = zi + 1; 
		if((!n) || (n >= mzl_cutoff)) continue;
		n = ne - ns; if(n > ((uint32_t)(0xffffffu))) n = 0xffffffu; n<<=8;

		for (zi = z->x+1; zi < idx_n && idx[zi].x == p->x; zi++) {
			if(idx[zi].rid == rid) continue; 
			x = idx[zi].rid; x <<= 1; x |= ((uint64_t)((z->rev == idx[zi].rev)?0:1));
			for (m = 0; m < tid_n; m++) {
				if(tid->a[m] == x) break;
			}
			if(m >= tid_n) continue; 


			y = &(idx[zi]); an = &ab->a[k++];
			rev = z->rev == y->rev? 0 : 1;
			an->other_off = rev?((uint32_t)-1)-1-(y->pos+1-y->span):y->pos;
			an->self_off = z->pos;
			///an->cnt: cnt<<8|span
			an->cnt = ((z->span <= ((uint32_t)(0xffu)))?z->span:((uint32_t)(0xffu))); an->cnt |= n;
			an->srt = (uint64_t)y->rid<<33 | (uint64_t)rev<<32 | an->self_off;
		}
		for (zi = ((int64_t)z->x)-1; zi >= 0 && idx[zi].x == p->x; zi--) {
			if(idx[zi].rid == rid) continue; 
			x = idx[zi].rid; x <<= 1; x |= ((uint64_t)((z->rev == idx[zi].rev)?0:1));
			for (m = 0; m < tid_n; m++) {
				if(tid->a[m] == x) break;
			}
			if(m >= tid_n) continue; 



			y = &(idx[zi]); an = &ab->a[k++];
			rev = z->rev == y->rev? 0 : 1;
			an->other_off = rev?((uint32_t)-1)-1-(y->pos+1-y->span):y->pos;
			an->self_off = z->pos;
			///an->cnt: cnt<<8|span
			an->cnt = ((z->span <= ((uint32_t)(0xffu)))?z->span:((uint32_t)(0xffu))); an->cnt |= n;
			an->srt = (uint64_t)y->rid<<33 | (uint64_t)rev<<32 | an->self_off;
		}
	}
	assert(k == ab->n_a);
	radix_sort_ha_an1(ab->a + ab->n_a - tn, ab->a + ab->n_a); kn = 0;
	for (k = ab->n_a - tn + 1, l = ab->n_a - tn; k <= ab->n_a; ++k) {
		if (k == ab->n_a || (ab->a[k].srt>>32) != (ab->a[l].srt>>32)) {
			for (; kn < tid_n && tid->a[kn] < (ab->a[l].srt>>32); kn++);
			if(kn < tid_n && tid->a[kn] == (ab->a[l].srt>>32)) kv_push(uint64_t, *tid, l);
			for (; kn < tid_n && tid->a[kn] == (ab->a[l].srt>>32); kn++);
			l = k;
		}
	}
}

void minimizers_qgen_input(ha_abufl_t *ab, uint64_t rid, char* rs, int64_t rl, uint64_t mz_w, uint64_t mz_k, Candidates_list *cl, kvec_t_u8_warp* k_flag, 
void *ha_flt_tab, ha_pt_t *ha_idx, All_reads* rdb, const ul_idx_t *udb, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, 
uint32_t *low_occ, ha_mzl_t *in, uint64_t in_n, ha_mzl_t *idx, int64_t idx_n, uint64_t mzl_cutoff, uint64_t chain_cutoff, kv_u_trans_t *kov)
{
	// fprintf(stderr, "+[M::%s]\n", __func__);
	uint64_t i, k, l, max_cnt = UINT32_MAX, min_cnt = 0, n; ha_mzl_t *z, *p, *y; anchor1_t *an;
	int64_t zi, ns, ne; uint8_t rev; k_mer_hit *s;
	if(high_occ) {
		max_cnt = (*high_occ);
		if(max_cnt < 2) max_cnt = 2;
	}
	if(low_occ) {
		min_cnt = (*low_occ);
		if(min_cnt < 2) min_cnt = 2;
	}
	clear_Candidates_list(cl); 

	for (i = 0, ab->n_a = 0; i < in_n; ++i) {
		///z is one of the minimizer
		z = &in[i]; p = &(idx[z->x]); n = 0;
		assert(z->pos == p->pos && z->rid == p->rid && z->span == p->span && z->rev == p->rev);
		for (zi = z->x+1; zi < idx_n && idx[zi].x == p->x && n < mzl_cutoff; zi++) {
			if(idx[zi].rid == rid) continue; n++;
		}
		for (zi = ((int64_t)z->x)-1; zi >= 0 && idx[zi].x == p->x && n < mzl_cutoff; zi--) {
			if(idx[zi].rid == rid) continue; n++;
		}
		if((!n) || (n >= mzl_cutoff)) continue;
		ab->n_a += n;
	}
	// if(rid == 0) {
	// 	fprintf(stderr, "-0-[M::%s] ab->n_a::%lu, in_n::%lu\n", __func__, ab->n_a, in_n);
	// }

	if (ab->n_a > ab->m_a) {
		ab->m_a = ab->n_a;
		REALLOC(ab->a, ab->m_a);
	}

	for (i = 0, k = 0; i < in_n; ++i) {
		///z is one of the minimizer
		z = &in[i]; p = &(idx[z->x]); n = 0; ns = 0; ne = -1;
		for (zi = z->x+1; zi < idx_n && idx[zi].x == p->x && n < mzl_cutoff; zi++) {
			if(idx[zi].rid == rid) continue; n++;
		}
		ne = zi;
		for (zi = ((int64_t)z->x)-1; zi >= 0 && idx[zi].x == p->x && n < mzl_cutoff; zi--) {
			if(idx[zi].rid == rid) continue; n++;
		}
		ns = zi + 1; 
		if((!n) || (n >= mzl_cutoff)) continue;
		n = ne - ns; if(n > ((uint32_t)(0xffffffu))) n = 0xffffffu; n<<=8;

		for (zi = z->x+1; zi < idx_n && idx[zi].x == p->x; zi++) {
			if(idx[zi].rid == rid) continue; 
			y = &(idx[zi]); an = &ab->a[k++];
			rev = z->rev == y->rev? 0 : 1;
			an->other_off = rev?((uint32_t)-1)-1-(y->pos+1-y->span):y->pos;
			an->self_off = z->pos;
			///an->cnt: cnt<<8|span
			an->cnt = ((z->span <= ((uint32_t)(0xffu)))?z->span:((uint32_t)(0xffu))); an->cnt |= n;
			an->srt = (uint64_t)y->rid<<33 | (uint64_t)rev<<32 | an->self_off;
		}
		for (zi = ((int64_t)z->x)-1; zi >= 0 && idx[zi].x == p->x; zi--) {
			if(idx[zi].rid == rid) continue; 
			y = &(idx[zi]); an = &ab->a[k++];
			rev = z->rev == y->rev? 0 : 1;
			an->other_off = rev?((uint32_t)-1)-1-(y->pos+1-y->span):y->pos;
			an->self_off = z->pos;
			///an->cnt: cnt<<8|span
			an->cnt = ((z->span <= ((uint32_t)(0xffu)))?z->span:((uint32_t)(0xffu))); an->cnt |= n;
			an->srt = (uint64_t)y->rid<<33 | (uint64_t)rev<<32 | an->self_off;
		}
	}
	assert(k == ab->n_a);

	uint64_t kn, spn; u_trans_t *ka; 
	kn = spn = 0; ka = NULL; sp->n = 0;
	if(kov) {
		kn = u_trans_n(*kov, rid); ka = u_trans_a(*kov, rid); 
	}
	if(kn > 0) {
		kv_resize(uint64_t, *sp, kn);
		for (k = 0; k < kn; ++k) {
			if(ka[k].del) continue;
			l = ka[k].tn; l <<= 1; l |= ka[k].rev; 
			kv_push(uint64_t, *sp, l);
		}
		if(sp->n > 1) radix_sort_anc64(sp->a, sp->a + sp->n);
	}

	radix_sort_ha_an1(ab->a, ab->a + ab->n_a); kn = 0;
	for (k = 1, l = n = 0, spn = sp->n; k <= ab->n_a; ++k) {
		if (k == ab->n_a || (ab->a[k].srt>>32) != (ab->a[l].srt>>32)) {
			if(k - l >= chain_cutoff) {
				for (; kn < spn && sp->a[kn] < (ab->a[l].srt>>32); kn++);
				if(kn < spn && sp->a[kn] == (ab->a[l].srt>>32)) kv_push(uint64_t, *sp, n);
				for (; kn < spn && sp->a[kn] == (ab->a[l].srt>>32); kn++) sp->a[kn] = (uint64_t)-1;
				for (i = l; i < k; i++) ab->a[n++] = ab->a[i];
			}
			l = k;
		}
	}
	ab->n_a = n;
	for (k = kn = 0; k < spn; k++) {
		if(sp->a[k] == (uint64_t)-1) continue;
		sp->a[kn++] = sp->a[k];
	}
	// sp->n = kn;
	if(kn > 0) gen_pair_chain(ab, rid, sp, kn, in, in_n, idx, idx_n, mzl_cutoff);
	if(spn > 0) {
		for (k = spn, kn = 0; k < sp->n; k++) sp->a[kn++] = sp->a[k];
		sp->n = kn;
	}
	

	// copy over to _cl_
	if (ab->n_a >= (uint64_t)cl->size) {
		cl->size = ab->n_a;
		REALLOC(cl->list, cl->size);
	}

    uint64_t tid = (uint64_t)-1, tl = (uint64_t)-1;
    for (k = 1, l = 0; k <= ab->n_a; ++k) {
        if (k == ab->n_a || ab->a[k].srt != ab->a[l].srt) {
            if (k-l>1) radix_sort_ha_an3(ab->a+l, ab->a+k);
            if((ab->a[l].srt>>33)!=tid) {
                tid = ab->a[l].srt>>33;
                tl = rdb?Get_READ_LENGTH((*rdb), tid):udb->ug->u.a[tid].len;
            } 
            for (i = l; i < k; i++) {
                s = &cl->list[i];
                s->readID = ab->a[i].srt>>33;
                s->strand = (ab->a[i].srt>>32)&1;
                if(!(s->strand)) {
                    s->offset = ab->a[i].other_off;
                } else {
                    s->offset = ((uint32_t)-1)-ab->a[i].other_off;
                    s->offset = tl-s->offset;
                }
                s->self_offset = ab->a[i].self_off;
                if(((ab->a[i].cnt>>8) < max_cnt) && ((ab->a[i].cnt>>8) > min_cnt)){
                    s->cnt = 1;
                } else if((ab->a[i].cnt>>8) <= min_cnt) {
                    s->cnt = 2;
                } else{
                    s->cnt = 1 + (((ab->a[i].cnt>>8) + (max_cnt<<1) - 1)/(max_cnt<<1));
                    s->cnt = pow(s->cnt, 1.1);
                }
                if(s->cnt > ((uint32_t)(0xffffffu))) s->cnt = 0xffffffu;
                s->cnt <<= 8; s->cnt |= (((uint32_t)(0xffu))&(ab->a[i].cnt));
            }
            l = k;
        }
    }
	cl->length = ab->n_a;
}

void inline reverse_k_mer_hit(k_mer_hit *a, uint64_t a_n, uint64_t xl, uint64_t yl)
{
	uint64_t z, han = a_n>>1; k_mer_hit *ai, *aj, ka;
	for (z = 0; z < han; z++) {
		ai = &(a[z]); aj = &(a[a_n-z-1]);
		ka = (*ai); (*ai) = (*aj); (*aj) = ka;

		ai->self_offset = xl-ai->self_offset-1;
		ai->offset = yl-ai->offset-1;

		aj->self_offset = xl-aj->self_offset-1;
		aj->offset = yl-aj->offset-1;
	}
	if(a_n&1) {
		a[z].self_offset = xl-a[z].self_offset-1;
		a[z].offset = yl-a[z].offset-1;
	}
}


void inline reset_k_mer_hit(k_mer_hit *a, uint64_t a_n, uint64_t xl, uint64_t yl, uint64_t rev, uint64_t *nid)
{
	uint64_t z, han = a_n>>1; k_mer_hit *ai, *aj, ka;
	if(rev) {
		for (z = 0; z < han; z++) {
			ai = &(a[z]); aj = &(a[a_n-z-1]); 
			ka = (*ai); (*ai) = (*aj); (*aj) = ka;

			ai->self_offset = xl-ai->self_offset-1;
			ai->offset = yl-ai->offset-1;

			aj->self_offset = xl-aj->self_offset-1;
			aj->offset = yl-aj->offset-1;
			if(nid) ai->readID = aj->readID = (*nid);
		}
		if(a_n&1) {
			a[z].self_offset = xl-a[z].self_offset-1;
			a[z].offset = yl-a[z].offset-1;
			if(nid) a[z].readID = (*nid);
		}
	} else if(nid) {
		for (z = 0; z < a_n; z++) a[z].readID = (*nid);
	}
}

void lchain_gen(Candidates_list* cl, overlap_region_alloc* ol, uint32_t rid, uint64_t rl, All_reads* rdb, 
				const ul_idx_t *udb, uint32_t apend_be, overlap_region* tf, uint64_t max_n_chain,
				int64_t max_skip, int64_t max_iter, int64_t max_dis, double chn_pen_gap, double chn_pen_skip, double bw_rate, int64_t quick_check, uint32_t gen_off)
{
	// fprintf(stderr, "+[M::%s]\n", __func__);
	uint64_t i, k, l, m, sm, cn = cl->length; overlap_region *r; ///srt = 0
	clear_overlap_region_alloc(ol);
	clear_fake_cigar(&(tf->f_cigar));

	for (l = 0, k = 1, m = 0; k <= cn; k++) {
		if((k == cn) || (cl->list[k].readID != cl->list[l].readID) 
													|| (cl->list[k].strand != cl->list[l].strand)) {
			if(cl->list[l].readID != rid) {
				tf->x_id = rid; 
				tf->x_pos_strand = cl->list[l].strand;
				tf->y_id = cl->list[l].readID; 
				tf->y_pos_strand = 0;///always 0
				// fprintf(stderr, "+[M::%s] l::%lu, k::%lu\n", __func__, l, k);
				sm = lchain_dp(cl->list+l, k-l, cl->list+m, &(cl->chainDP), tf, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_rate, 
              	rl, rdb?Get_READ_LENGTH((*rdb), (*tf).y_id):udb->ug->u.a[(*tf).y_id].len, quick_check);
				// assert(sm > 0);
				if(ovlp_chain_gen(ol, tf, rl, rdb?Get_READ_LENGTH((*rdb), (*tf).y_id):udb->ug->u.a[(*tf).y_id].len, apend_be, cl->list+m, sm)) {
					r = &(ol->list[ol->length-1]); r->non_homopolymer_errors = m; 
					// if(r->y_pos_strand) {
					// 	reverse_k_mer_hit(cl->list+m, sm, rl, rdb?Get_READ_LENGTH((*rdb), r->y_id):udb->ug->u.a[r->y_id].len);
					// }
					reset_k_mer_hit(cl->list+m, sm, rl, rdb?Get_READ_LENGTH((*rdb), r->y_id):udb->ug->u.a[r->y_id].len, r->y_pos_strand, &(ol->length));
					if(gen_off) gen_fake_cigar(&(r->f_cigar), r, apend_be, cl->list+m, sm);
					m += sm;
				}
			}
			l = k;
		}
	}
	cl->length = m;


	k = ol->length;
	if (ol->length > max_n_chain) {
		int32_t w, n[4], s[4]; overlap_region t;
		n[0] = n[1] = n[2] = n[3] = 0, s[0] = s[1] = s[2] = s[3] = 0;
		ks_introsort_or_ss(ol->length, ol->list); ///srt = 1;
		for (i = 0; i < ol->length; ++i) {
			r = &(ol->list[i]);
			w = ha_ov_type(r, rl);
			++n[w];
			if (((uint64_t)n[w]) == max_n_chain) s[w] = r->shared_seed;
		}
		if (s[0] > 0 || s[1] > 0 || s[2] > 0 || s[3] > 0) {
			// n[0] = n[1] = n[2] = n[3] = 0;
			for (i = 0, k = 0; i < ol->length; ++i) {
				r = &(ol->list[i]);
				w = ha_ov_type(r, rl);
				// ++n[w];
				// if (((int)n[w] <= max_n_chain) || (r->shared_seed >= s[w] && s[w] >= (asm_opt.k_mer_length<<1))) {
				if (r->shared_seed >= s[w]) {
					if (k != i) {
						t = ol->list[k];
						ol->list[k] = ol->list[i];
						ol->list[i] = t;
					}
					++k;
				}
			}
			ol->length = k;
		}
	}
	/**
	if(!gen_off) {
		if(srt) ks_introsort_or_id(ol->length, ol->list);
		uint64_t cln = cl->length;
		for (i = k = m = 0; i < ol->length; ++i) {
			r = &(ol->list[i]);
			for (;(k<cln)&&((cl->list[k].readID!=r->y_id)||(cl->list[k].strand!=r->y_pos_strand)); k++);
			// assert(k<cln);
			for (cn = m; k < cln; k++, m++) {
				if((cl->list[k].readID!=r->y_id)||(cl->list[k].strand!=r->y_pos_strand)) break;
				if(m != k) cl->list[m] = cl->list[k];
			}
			// assert(m - cn > 0);
			if(r->y_pos_strand) {
				reverse_k_mer_hit(cl->list+cn, m-cn, rl, rdb?Get_READ_LENGTH((*rdb), r->y_id):udb->ug->u.a[r->y_id].len);
			}
			// gen_fake_cigar(&(r->f_cigar), r, apend_be, cl->list+cn, m-cn);
		}
		cl->length = m;
	}
	**/

	ks_introsort_or_xs(ol->length, ol->list);
}

void lchain_qgen(Candidates_list* cl, overlap_region_alloc* ol, uint32_t rid, uint64_t rl, All_reads* rdb, 
				const ul_idx_t *udb, uint32_t apend_be, overlap_region* tf, uint64_t max_n_chain,
				int64_t max_skip, int64_t max_iter, int64_t max_dis, double chn_pen_gap, double chn_pen_skip, double bw_rate, int64_t quick_check, uint32_t gen_off)
{
	uint64_t i, k, l, m, sm, cn = cl->length; overlap_region *r; ///srt = 0
	clear_overlap_region_alloc(ol);
	clear_fake_cigar(&(tf->f_cigar));

	for (l = 0, k = 1, m = 0; k <= cn; k++) {
		if((k == cn) || (cl->list[k].readID != cl->list[l].readID) 
													|| (cl->list[k].strand != cl->list[l].strand)) {
			if(cl->list[l].readID != rid) {
				tf->x_id = rid; 
				tf->x_pos_strand = cl->list[l].strand;
				tf->y_id = cl->list[l].readID; 
				tf->y_pos_strand = 0;///always 0
				// fprintf(stderr, "+[M::%s] l::%lu, k::%lu\n", __func__, l, k);
				sm = lchain_qdp(cl->list+l, k-l, cl->list+m, &(cl->chainDP), tf, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_rate, 
              	rl, rdb?Get_READ_LENGTH((*rdb), (*tf).y_id):udb->ug->u.a[(*tf).y_id].len, quick_check);
				// assert(sm > 0);
				if(ovlp_chain_qgen(ol, tf, rl, rdb?Get_READ_LENGTH((*rdb), (*tf).y_id):udb->ug->u.a[(*tf).y_id].len, apend_be, cl->list+m, sm)) {
					r = &(ol->list[ol->length-1]); r->non_homopolymer_errors = m; 
					// if(tf->y_id == 66 || tf->y_id == 66) {
					// 	fprintf(stderr, "\n[M::%s::] utg%.6dl(%c), i::%lu\n", 
        			// 			__func__, (int32_t)tf->y_id+1, "+-"[tf->x_pos_strand], m);
					// }
					// reset_k_mer_hit(cl->list+m, sm, rl, rdb?Get_READ_LENGTH((*rdb), r->y_id):udb->ug->u.a[r->y_id].len, r->y_pos_strand, &(ol->length));
					for (i = 0; i < sm; i++) {
						cl->list[m+i].readID = ol->length;
						// if(tf->y_id == 126) fprintf(stderr, "[M::%s::qoff->%u::toff->%u]\n", __func__, cl->list[m+i].self_offset, cl->list[m+i].offset);
					}
					if(gen_off) gen_fake_cigar(&(r->f_cigar), r, apend_be, cl->list+m, sm);
					m += sm;
				}
			}
			l = k;
		}
	}
	cl->length = m;
	// fprintf(stderr, "+[M::%s] cn::%lu, m::%lu, ol->length::%lu\n", __func__, cn, m, ol->length);
	// for (k = 0; k < ol->length; k++) {
	// 	fprintf(stderr, "---[M::%s::utg%.6dl] q[%d, %d), t[%d, %d), khit_off::%u\n", __func__, 
    //     (int32_t)ol->list[k].y_id+1, ol->list[k].x_pos_s, ol->list[k].x_pos_e+1,
	// 	ol->list[k].y_pos_s, ol->list[k].y_pos_e+1, ol->list[k].non_homopolymer_errors);
	// }


	k = ol->length;
	if (ol->length > max_n_chain) {
		int32_t w, n[4], s[4]; overlap_region t;
		n[0] = n[1] = n[2] = n[3] = 0, s[0] = s[1] = s[2] = s[3] = 0;
		ks_introsort_or_ss(ol->length, ol->list); ///srt = 1;
		for (i = 0; i < ol->length; ++i) {
			r = &(ol->list[i]);
			w = ha_ov_type(r, rl);
			++n[w];
			if (((uint64_t)n[w]) == max_n_chain) s[w] = r->shared_seed;
		}
		if (s[0] > 0 || s[1] > 0 || s[2] > 0 || s[3] > 0) {
			// n[0] = n[1] = n[2] = n[3] = 0;
			for (i = 0, k = 0; i < ol->length; ++i) {
				r = &(ol->list[i]);
				w = ha_ov_type(r, rl);
				// ++n[w];
				// if (((int)n[w] <= max_n_chain) || (r->shared_seed >= s[w] && s[w] >= (asm_opt.k_mer_length<<1))) {
				if (r->shared_seed >= s[w]) {
					if (k != i) {
						t = ol->list[k];
						ol->list[k] = ol->list[i];
						ol->list[i] = t;
					}
					++k;
				}
			}
			ol->length = k;
		}
	}
	ks_introsort_or_xs(ol->length, ol->list);
}


void lchain_qgen_mcopy(Candidates_list* cl, overlap_region_alloc* ol, uint32_t rid, uint64_t rl, All_reads* rdb, 
				const ul_idx_t *udb, uint32_t apend_be, uint64_t max_n_chain, int64_t max_skip, int64_t max_iter, 
				int64_t max_dis, double chn_pen_gap, double chn_pen_skip, double bw_rate, int64_t quick_check, 
				uint32_t gen_off, double mcopy_rate, uint32_t chain_cutoff, uint32_t mcopy_khit_cut, st_mt_t *sp)
{
	// fprintf(stderr, "+[M::%s]\n", __func__);
	uint64_t i, k, l, m, cn = cl->length, yid, ol0, lch; overlap_region *r, t; ///srt = 0
	clear_overlap_region_alloc(ol);

	for (l = 0, k = 1, m = 0, lch = 0; k <= cn; k++) {
		if((k == cn) || (cl->list[k].readID != cl->list[l].readID)) {
			if(cl->list[l].readID != rid) {
				yid = cl->list[l].readID; ol0 = ol->length;
				m += lchain_qdp_mcopy(cl, l, k-l, m, &(cl->chainDP), ol, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_rate, 
				rid, rl, rdb?Get_READ_LENGTH((*rdb), yid):udb->ug->u.a[yid].len, quick_check, apend_be, gen_off, 1, mcopy_rate, mcopy_khit_cut, 1);
				if((chain_cutoff >= 2) && (!lch)) {
					for (i = ol0; (i<ol->length) && (!lch); i++) {
						if(ol->list[i].align_length < chain_cutoff) lch = 1;
					}
				}
			}
			l = k;
		}
	}
	cl->length = m;

	// for (k = 0; k < ol->length; k++) {
	// 	fprintf(stderr, "---[M::%s::utg%.6dl] q[%d, %d), t[%d, %d), khit_off::%u\n", __func__, 
    //     (int32_t)ol->list[k].y_id+1, ol->list[k].x_pos_s, ol->list[k].x_pos_e+1,
	// 	ol->list[k].y_pos_s, ol->list[k].y_pos_e+1, ol->list[k].non_homopolymer_errors);
	// }

	k = ol->length;
	if (ol->length > max_n_chain) {
		int32_t w, n[4], s[4]; 
		n[0] = n[1] = n[2] = n[3] = 0, s[0] = s[1] = s[2] = s[3] = 0; 
		ks_introsort_or_ss(ol->length, ol->list); 
		for (i = 0; i < ol->length; ++i) {
			r = &(ol->list[i]);
			w = ha_ov_type(r, rl);
			++n[w];
			if (((uint64_t)n[w]) == max_n_chain) s[w] = r->shared_seed;
		}
		if (s[0] > 0 || s[1] > 0 || s[2] > 0 || s[3] > 0) {
			// n[0] = n[1] = n[2] = n[3] = 0;
			for (i = 0, k = 0, lch = 0; i < ol->length; ++i) {
				r = &(ol->list[i]);
				w = ha_ov_type(r, rl);
				// ++n[w];
				// if (((int)n[w] <= max_n_chain) || (r->shared_seed >= s[w] && s[w] >= (asm_opt.k_mer_length<<1))) {
				if (r->shared_seed >= s[w]) {
					if (k != i) {
						t = ol->list[k];
						ol->list[k] = ol->list[i];
						ol->list[i] = t;
					}
					if(ol->list[k].align_length < chain_cutoff) lch = 1;
					++k;
				}
			}
			ol->length = k;
		}
	} 
	
	ks_introsort_or_xs(ol->length, ol->list);
	if(lch) {
		//@brief r485
		uint64_t zs, ze, rs, re, ob, os, oe, ocn, pp, kn, ms, me; int64_t osc;
		for (i = l = 0, cn = cl->length; i < ol->length; ++i) {
			if(ol->list[i].align_length < chain_cutoff) {
				zs = ol->list[i].x_pos_s; ze = ol->list[i].x_pos_e + 1;
				ob = (ze - zs)*OFL; if(ob < 16) ob = 16;
				osc = ol->list[i].shared_seed*CH_SC;
				ocn = ol->list[i].align_length<<CH_OCC;
				for (k = 0; (k < ol->length) && (ze > ol->list[k].x_pos_s); k++) {
					if(ol->list[k].align_length < chain_cutoff) continue;
					if(ol->list[k].align_length < ocn) continue;
					if(ol->list[k].shared_seed < osc) continue;
					rs = ol->list[k].x_pos_s; re = ol->list[k].x_pos_e + 1;
					os = ((rs>=zs)?rs:zs); oe = ((re<=ze)?re:ze);
					if((oe > os) && (oe - os) >= ob) {
						m = ol->list[k].non_homopolymer_errors;
						pp = cl->list[m].readID; kn = 0;
						for (; (m < cn) && (cl->list[m].readID == pp) && (kn < ocn); m++) {
							me = cl->list[m].self_offset; ms = me - (cl->list[m].cnt&(0xffu)); 
							if((ms >= os) && (me <= oe)) kn++;
						}
						if(kn >= ocn) break;
					}
				}
				if((k < ol->length) && (ze > ol->list[k].x_pos_s)) continue;
			}
			if (l != i) {
				t = ol->list[l];
				ol->list[l] = ol->list[i];
				ol->list[i] = t;
			}
			l++;
		}
		// fprintf(stderr, "+[M::%s] rid::%u, ol->length0::%lu, ol->length1::%lu\n", __func__, rid, ol->length, l);
		ol->length = l;


		/**
		//@brief r484
		for (i = sp->n = 0; i < ol->length; ++i) {
			if(ol->list[i].align_length < chain_cutoff) continue;
			os = ol->list[i].x_pos_s; oe = ol->list[i].x_pos_e + 1;
			if((sp->n) && (((uint32_t)sp->a[sp->n-1]) >= os)) {
				if(oe > ((uint32_t)sp->a[sp->n-1])) {
					oe = oe - ((uint32_t)sp->a[sp->n-1]);
					sp->a[sp->n-1] += oe; 
				}
			} else {
				os = (os<<32)|oe; kv_push(uint64_t, *sp, os);
			}
		}

		for (i = k = 0; i < ol->length; ++i) {
			if(ol->list[i].align_length < chain_cutoff) {///ol has been sorted by x_pos_s
				r = &(ol->list[i]); rs = r->x_pos_s; re = r->x_pos_e + 1;
				rl = re - rs; ovl = 0;
				for (m = 0; (m < sp->n) && (re > (sp->a[m]>>32)); m++) {
					os = ((rs>=(sp->a[m]>>32))? rs:(sp->a[m]>>32));
            		oe = ((re<=((uint32_t)sp->a[m]))? re:((uint32_t)sp->a[m]));
					if(oe > os) {
						ovl += (oe - os); if(ovl >= (rl*0.95)) break;
					}
				}
				if(ovl >= (rl*0.95)) continue;
			}
			if (k != i) {
				t = ol->list[k];
				ol->list[k] = ol->list[i];
				ol->list[i] = t;
			}
			ol->list[k++].align_length = 0;
		}
		// fprintf(stderr, "+[M::%s] ol->length0::%lu, ol->length1::%lu\n", __func__, ol->length, k);
		ol->length = k;
		**/
	} 
	/**else {
		for (i = 0; i < ol->length; ++i) ol->list[i].align_length = 0;
	}
	**/
	for (i = 0; i < ol->length; ++i) ol->list[i].align_length = 0;
}

void lchain_qgen_mcopy_fast(Candidates_list* cl, overlap_region_alloc* ol, uint32_t rid, uint64_t rl, All_reads* rdb, 
				uint32_t apend_be, uint64_t max_n_chain, int64_t max_skip, int64_t max_iter, 
				int64_t max_dis, double chn_pen_gap, double chn_pen_skip, double bw_rate, int64_t quick_check, 
				uint32_t gen_off, int64_t mcopy_num, double mcopy_rate, uint32_t chain_cutoff, uint32_t mcopy_khit_cut, st_mt_t *sp, uint64_t ocv_w)
{
	// fprintf(stderr, "+[M::%s] chain_cutoff::%u\n", __func__, chain_cutoff);
	uint64_t i, k, l, m, cn = cl->length, yid, ol0, lch, *cc = NULL, cwn = 0, cws, cwe, os, oe, rs, re, cw0, cw1/**, dbgn = 0**/; overlap_region *r, t; ///srt = 0
	clear_overlap_region_alloc(ol);

	for (l = 0, k = 1, m = 0, lch = 0; k <= cn; k++) {
		if((k == cn) || (cl->list[k].readID != cl->list[l].readID)) {
			if(cl->list[l].readID != rid) {
				yid = cl->list[l].readID; ol0 = ol->length;
				m += lchain_qdp_mcopy_fast(cl, l, k-l, m, &(cl->chainDP), ol, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_rate, 
				rid, rl, Get_READ_LENGTH((*rdb), yid), quick_check, apend_be, gen_off, mcopy_num, mcopy_rate, mcopy_khit_cut, 1);
				if((chain_cutoff >= 2) && (!lch)) {
					for (i = ol0; (i<ol->length) && (!lch); i++) {
						if(ol->list[i].align_length < chain_cutoff) lch = 1;
					}
				}
			}
			l = k;
		}
	}
	cl->length = m;

	// fprintf(stderr, "[M::%s::] rn::%lu\tmax_n_chain::%lu\n", __func__, ol->length, max_n_chain);
	// for (k = 0; k < ol->length; k++) {
	// 	fprintf(stderr, "---[M::%s::%.*s(qid::%u)] q[%d, %d), t[%d, %d), sc::%d, type::%d\n", __func__, 
	// 	(int32_t)Get_NAME_LENGTH(R_INF, ol->list[k].y_id), Get_NAME(R_INF, ol->list[k].y_id), ol->list[k].y_id,
	// 	ol->list[k].x_pos_s, ol->list[k].x_pos_e+1, ol->list[k].y_pos_s, ol->list[k].y_pos_e+1, ol->list[k].shared_seed, ha_ov_type(&(ol->list[k]), rl));
	// }

	k = ol->length;
	if (ol->length > max_n_chain) {
		int32_t w, n[4], s[4]; 
		n[0] = n[1] = n[2] = n[3] = 0, s[0] = s[1] = s[2] = s[3] = 0; 
		ks_introsort_or_ss(ol->length, ol->list); 
		for (i = 0; i < ol->length; ++i) {
			r = &(ol->list[i]);
			w = ha_ov_type(r, rl);
			++n[w];
			if (((uint64_t)n[w]) == max_n_chain) s[w] = r->shared_seed;
		}
		// fprintf(stderr, "[M::%s::] s[0]::%d, s[1]::%d, s[2]::%d, s[3]::%d\n", __func__, s[0], s[1], s[2], s[3]);
		if (s[0] > 0 || s[1] > 0 || s[2] > 0 || s[3] > 0) {
			if((((uint64_t)n[3]) >= max_n_chain) && (rl >= ocv_w)) {
				cwn = (rl/ocv_w) + ((rl % ocv_w)?(1):(0));
				kv_resize(uint64_t, (*sp), (cwn));
				cc = sp->a;
				for (i = cws = cwe = 0; i < cwn; i++) {
					cwe = cws + ocv_w; if(cwe > rl) cwe = rl;
					assert(cwe > cws);
					cc[i] = (cwe - cws)*(max_n_chain>>1); 
					// fprintf(stderr, "[M::%s::] cw::[%lu, %lu), cc::%lu\n", __func__, cws, cwe, cc[i]);
					if(cc[i] > UINT32_MAX) cc[i] = UINT32_MAX; cc[i] <<= 32;
					cws += ocv_w;
				}
			}
			for (i = 0, k = 0, lch = 0; i < ol->length; ++i) {
				r = &(ol->list[i]);
				w = ha_ov_type(r, rl);
				// ++n[w];
				// if (((int)n[w] <= max_n_chain) || (r->shared_seed >= s[w] && s[w] >= (asm_opt.k_mer_length<<1))) {
				if (r->shared_seed >= s[w]) {
					if(cwn) {
						m = (ol->list[i].x_pos_s/ocv_w); 
						rs = ol->list[i].x_pos_s; re = ol->list[i].x_pos_e + 1;
						for (cws = m*ocv_w; m < cwn; m++) {
							cwe = cws + ocv_w; if(cwe > rl) cwe = rl;
							os = ((rs >= cws)? rs : cws);
        					oe = ((re <= cwe)? re : cwe);
							if(oe <= os) break;
							if(((uint32_t)cc[m]) + (oe - os) < UINT32_MAX) {
								cc[m] += (oe - os);
							} else {
								cc[m] >>= 32; cc[m] <<= 32; cc[m] |= UINT32_MAX;
							}
							cws += ocv_w;
						}
					}
					if (k != i) {
						t = ol->list[k];
						ol->list[k] = ol->list[i];
						ol->list[i] = t;
					}
					if(ol->list[k].align_length < chain_cutoff) lch = 1;
					++k;
				} else if(w == 3 && cwn > 0) {
					m = (ol->list[i].x_pos_s/ocv_w); cw0 = cw1 = 0;
					rs = ol->list[i].x_pos_s; re = ol->list[i].x_pos_e + 1;
					for (cws = m*ocv_w; m < cwn; m++) {
						cwe = cws + ocv_w; if(cwe > rl) cwe = rl;
						os = ((rs >= cws)? rs : cws);
						oe = ((re <= cwe)? re : cwe);
						if(oe <= os) break;
						// fprintf(stderr, "+++[M::%s::%.*s(qid::%u)] o[%lu, %lu), cur::%lu, max::%lu\n", __func__, 
						// 	(int32_t)Get_NAME_LENGTH(R_INF, ol->list[i].y_id), Get_NAME(R_INF, ol->list[i].y_id), ol->list[i].y_id,
						// 	os, oe, ((uint64_t)((uint32_t)cc[m])), (cc[m]>>32));
						if((oe - os) + ((uint64_t)((uint32_t)cc[m])) >= (cc[m]>>32)) {
							cw1 += (oe - os);
						} else {
							cw0 += (oe - os);
						}
						cws += ocv_w;
					}

					if(cw0 >= ((cw0 + cw1)*0.7)) {
						m = (ol->list[i].x_pos_s/ocv_w); 
						rs = ol->list[i].x_pos_s; re = ol->list[i].x_pos_e + 1;
						for (cws = m*ocv_w; m < cwn; m++) {
							cwe = cws + ocv_w; if(cwe > rl) cwe = rl;
							os = ((rs >= cws)? rs : cws);
        					oe = ((re <= cwe)? re : cwe);
							if(oe <= os) break;
							if(((uint32_t)cc[m]) + (oe - os) < UINT32_MAX) {
								cc[m] += (oe - os);
							} else {
								cc[m] >>= 32; cc[m] <<= 32; cc[m] |= UINT32_MAX;
							}
							cws += ocv_w;
						}
						if (k != i) {
							t = ol->list[k];
							ol->list[k] = ol->list[i];
							ol->list[i] = t;
						}
						if(ol->list[k].align_length < chain_cutoff) lch = 1;
						++k;
						// dbgn++;
						// fprintf(stderr, "+++[M::%s::%.*s(qid::%u)] q[%d, %d), t[%d, %d), sc::%d, type::%d, cw0::%lu, cw1::%lu\n", __func__, 
						// 	(int32_t)Get_NAME_LENGTH(R_INF, ol->list[k].y_id), Get_NAME(R_INF, ol->list[k].y_id), ol->list[k].y_id,
						// 	ol->list[k].x_pos_s, ol->list[k].x_pos_e+1, ol->list[k].y_pos_s, ol->list[k].y_pos_e+1, ol->list[k].shared_seed, ha_ov_type(&(ol->list[k]), rl), cw0, cw1);
					}
				}
			}
			ol->length = k;
		}
	} 
	
	ks_introsort_or_xs(ol->length, ol->list);
	if(lch) {
		//@brief r485
		uint64_t zs, ze, rs, re, ob, os, oe, ocn, pp, kn, ms, me; int64_t osc;
		for (i = l = 0, cn = cl->length; i < ol->length; ++i) {
			if(ol->list[i].align_length < chain_cutoff) {
				zs = ol->list[i].x_pos_s; ze = ol->list[i].x_pos_e + 1;
				ob = (ze - zs)*OFL; if(ob < 16) ob = 16;
				osc = ol->list[i].shared_seed*CH_SC;
				ocn = ol->list[i].align_length<<CH_OCC;
				for (k = 0; (k < ol->length) && (ze > ol->list[k].x_pos_s); k++) {
					if(ol->list[k].align_length < chain_cutoff) continue;
					if(ol->list[k].align_length < ocn) continue;
					if(ol->list[k].shared_seed < osc) continue;
					rs = ol->list[k].x_pos_s; re = ol->list[k].x_pos_e + 1;
					os = ((rs>=zs)?rs:zs); oe = ((re<=ze)?re:ze);
					if((oe > os) && (oe - os) >= ob) {
						m = ol->list[k].non_homopolymer_errors;
						pp = cl->list[m].readID; kn = 0;
						for (; (m < cn) && (cl->list[m].readID == pp) && (kn < ocn); m++) {
							me = cl->list[m].self_offset; ms = me - (cl->list[m].cnt&(0xffu)); 
							if((ms >= os) && (me <= oe)) kn++;
						}
						if(kn >= ocn) break;
					}
				}
				if((k < ol->length) && (ze > ol->list[k].x_pos_s)) continue;
			}
			if (l != i) {
				t = ol->list[l];
				ol->list[l] = ol->list[i];
				ol->list[i] = t;
			}
			l++;
		}
		ol->length = l;
	} 

	for (i = 0; i < ol->length; ++i) ol->list[i].align_length = 0;
	// fprintf(stderr, "+[M::%s] rid::%u, ol->length0::%lu, dbgn::%lu\n", __func__, rid, ol->length, dbgn);
}

void lchain_qgen_mcopy_fast_re1(Candidates_list* cl, uint32_t cl_beg, overlap_region_alloc* ol, uint32_t rid, uint64_t rl, uint64_t tl, 
				uint32_t apend_be, int64_t max_skip, int64_t max_iter, 
				int64_t max_dis, double chn_pen_gap, double chn_pen_skip, double bw_rate, int64_t quick_check, 
				uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut)
{
	uint64_t cn = cl->length, m = cl_beg;
	if(cl_beg >= cn) return;

	m += lchain_qdp_mcopy_fast(cl, cl_beg, cn - cl_beg, m, &(cl->chainDP), ol, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_rate, 
                rid, rl, tl, quick_check, apend_be, gen_off, enable_mcopy, mcopy_rate, mcopy_khit_cut, 1);
	cl->length = m;
}

inline uint64_t special_lchain(Candidates_list* cl, overlap_region_alloc* ol, uint32_t rid, uint64_t rl, All_reads* rdb, 
				const ul_idx_t *udb, uint32_t apend_be, int64_t max_skip, int64_t max_iter, int64_t max_dis, double chn_pen_gap, 
				double chn_pen_skip, double bw_rate, int64_t quick_check, uint32_t gen_off, double mcopy_rate, uint32_t mcopy_khit_cut, 
				st_mt_t *sp, uint64_t *si, uint64_t m, uint64_t l, uint64_t k)
{
	uint64_t z, s, e, yid, ol0 = ol->length, ol1, m0 = m; 
	if(l >= k) return m;
	yid = cl->list[l].readID;
	for (z = l; z < k && cl->list[z].strand == cl->list[l].strand; z++);
	if(z > l) {
		s = l; e = z; ol1 = ol->length; // m0 = m;
		m += lchain_qdp_mcopy(cl, s, e-s, m, &(cl->chainDP), ol, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_rate, 
						rid, rl, rdb?Get_READ_LENGTH((*rdb), yid):udb->ug->u.a[yid].len, quick_check, apend_be, gen_off, 0, mcopy_rate, mcopy_khit_cut, 0);
		// if(rid == 7) {
		// 	fprintf(stderr, "+[M::%s::utg%.6ul] inn::[%lu, %lu), (*si)::%lu, ol1::%lu, ol->length::%lu\n", 
		// 	__func__, rid+1, s, e, (*si), ol1, ol->length);
		// }
		if(((*si) < sp->n) && (sp->a[(*si)] == s) && (ol->length > ol1)) {
			ol->list[ol->length-1].x_pos_strand = 1;
			// fprintf(stderr, "+[M::%s] utg%.6u%c -> utg%.6u%c, n_khits0::%lu, n_khits::%lu\n", __func__, rid+1, "lc"[udb->ug->u.a[rid].circ], 
			// ol->list[ol->length-1].y_id+1, "lc"[udb->ug->u.a[ol->list[ol->length-1].y_id].circ], e-s, m-m0);
		}
		for (; (*si) < sp->n && sp->a[(*si)] == s; (*si)++);
	}

	if(z < k) {
		s = z; e = k; ol1 = ol->length;	// m0 = m;
		m += lchain_qdp_mcopy(cl, s, e-s, m, &(cl->chainDP), ol, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_rate, 
						rid, rl, rdb?Get_READ_LENGTH((*rdb), yid):udb->ug->u.a[yid].len, quick_check, apend_be, gen_off, 0, mcopy_rate, mcopy_khit_cut, 0);
		// if(rid == 7) {
		// 	fprintf(stderr, "-[M::%s::utg%.6ul] inn::[%lu, %lu), (*si)::%lu, ol1::%lu, ol->length::%lu\n", 
		// 	__func__, rid+1, s, e, (*si), ol1, ol->length);
		// }
		if(((*si) < sp->n) && (sp->a[(*si)] == s) && (ol->length > ol1)) {
			ol->list[ol->length-1].x_pos_strand = 1;
			// fprintf(stderr, "-[M::%s] utg%.6u%c -> utg%.6u%c, n_khits0::%lu, n_khits::%lu\n", __func__, rid+1, "lc"[udb->ug->u.a[rid].circ], 
			// ol->list[ol->length-1].y_id+1, "lc"[udb->ug->u.a[ol->list[ol->length-1].y_id].circ], e-s, m-m0);
		}
		for (; (*si) < sp->n && sp->a[(*si)] == s; (*si)++);
	}
	if(ol->length > ol0) {
		overlap_region *p = &(ol->list[ol0]); uint64_t pi = ol0; overlap_region t;
		for (z = ol0; z < ol->length; z++) {
			if((p->shared_seed < ol->list[z].shared_seed) || 
				((p->shared_seed == ol->list[z].shared_seed) && (ol->list[z].x_pos_strand == 1))) {
				p = &(ol->list[z]); pi = z;
			}
		}

		// if(rid == 7) {
		// 	fprintf(stderr, ">[M::%s::] utg%.6ul\t->\tutg%.6ul\tflag::%u\n", 
		// 	__func__, p->x_id+1, p->y_id+1, p->x_pos_strand);
		// }

		for (z = s = ol0; z < ol->length; z++) {
			if((ol->list[z].x_pos_strand == 0) && (z != pi)) continue;
			if(s != z) {
				t = ol->list[s];
    			ol->list[s] = ol->list[z];
    			ol->list[z] = t;
			}
			s++;
		}

		if(s != ol->length) {
			ol->length = s;
			for (z = ol0; z < ol->length; z++) {
				s = ol->list[z].non_homopolymer_errors;
				pi = cl->list[s].readID;
				ol->list[z].non_homopolymer_errors = m0;
				for (; s < m && cl->list[s].readID == pi; s++) {
					cl->list[m0] = cl->list[s];
					cl->list[m0].readID = z; m0++;
				}
			}
			m = m0;
		}
		
	}
	return m;
}

void lchain_qgen_mcopy_input(Candidates_list* cl, overlap_region_alloc* ol, uint32_t rid, uint64_t rl, All_reads* rdb, 
				const ul_idx_t *udb, uint32_t apend_be, uint64_t max_n_chain, int64_t max_skip, int64_t max_iter, 
				int64_t max_dis, double chn_pen_gap, double chn_pen_skip, double bw_rate, double bw_thres_sec, 
				int64_t quick_check, uint32_t gen_off, double mcopy_rate, uint32_t mcopy_khit_cut, st_mt_t *sp)
{
	// fprintf(stderr, "+[M::%s]\n", __func__);
	uint64_t i, k, l, m, cn = cl->length, yid, si; overlap_region *r; ///srt = 0
	clear_overlap_region_alloc(ol);

	// if(rid == 160) {
	// 	fprintf(stderr, "[M::%s::utg%.6ul] sp->n::%u\n", __func__, rid+1, (uint32_t)sp->n);
	// }
	for (l = 0, k = 1, m = 0, si = 0; k <= cn; k++) {
		if((k == cn) || (cl->list[k].readID != cl->list[l].readID)) {
			if(cl->list[l].readID != rid) {
				yid = cl->list[l].readID;
				// if(rid == 160) {
				// 	fprintf(stderr, "+[M::%s::utg%.6lul] [%lu, %lu), m::%lu, ol->length::%lu\n", __func__, yid+1, l, k, m, ol->length);
				// }
				for (; si < sp->n && sp->a[si] < l; si++);
				if(si >= sp->n || sp->a[si] >= k) {///might be unmatched
					m += lchain_qdp_mcopy(cl, l, k-l, m, &(cl->chainDP), ol, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_rate, 
						rid, rl, rdb?Get_READ_LENGTH((*rdb), yid):udb->ug->u.a[yid].len, quick_check, apend_be, gen_off, 0, mcopy_rate, mcopy_khit_cut, 0);
				} else {
					m = special_lchain(cl, ol, rid, rl, rdb, udb, apend_be, max_skip, max_iter, max_dis, chn_pen_gap, 
					chn_pen_skip, bw_thres_sec, quick_check, gen_off, mcopy_rate, mcopy_khit_cut, sp, &si, m, l, k);
				}
				// if(rid == 160) {
				// 	fprintf(stderr, "-[M::%s::utg%.6lul] [%lu, %lu), m::%lu, ol->length::%lu\n", __func__, yid+1, l, k, m, ol->length);
				// }
			}
			l = k;
		}
	}
	cl->length = m;

	// for (k = 0; k < ol->length; k++) {
	// 	fprintf(stderr, "---[M::%s::utg%.6dl] q[%d, %d), t[%d, %d), khit_off::%u\n", __func__, 
    //     (int32_t)ol->list[k].y_id+1, ol->list[k].x_pos_s, ol->list[k].x_pos_e+1,
	// 	ol->list[k].y_pos_s, ol->list[k].y_pos_e+1, ol->list[k].non_homopolymer_errors);
	// }

	k = ol->length;
	if (ol->length > max_n_chain) {
		int32_t w, n[4], s[4]; overlap_region t;
		n[0] = n[1] = n[2] = n[3] = 0, s[0] = s[1] = s[2] = s[3] = 0;
		ks_introsort_or_ss(ol->length, ol->list); ///srt = 1;
		for (i = 0; i < ol->length; ++i) {
			r = &(ol->list[i]);
			w = ha_ov_type(r, rl);
			++n[w];
			if (((uint64_t)n[w]) == max_n_chain) s[w] = r->shared_seed;
		}
		if (s[0] > 0 || s[1] > 0 || s[2] > 0 || s[3] > 0) {
			// n[0] = n[1] = n[2] = n[3] = 0;
			for (i = 0, k = 0; i < ol->length; ++i) {
				r = &(ol->list[i]);
				w = ha_ov_type(r, rl);
				// ++n[w];
				// if (((int)n[w] <= max_n_chain) || (r->shared_seed >= s[w] && s[w] >= (asm_opt.k_mer_length<<1))) {
				if ((r->shared_seed >= s[w]) || (r->x_pos_strand == 1)) {
					if (k != i) {
						t = ol->list[k];
						ol->list[k] = ol->list[i];
						ol->list[i] = t;
					}
					++k;
				}
			}
			ol->length = k;
		}
	}
	ks_introsort_or_xs(ol->length, ol->list);
}

void set_lchain_dp_op(uint32_t is_accurate, uint32_t mz_k, int64_t *max_skip, int64_t *max_iter, int64_t *max_dis, double *chn_pen_gap, double *chn_pen_skip, int64_t *quick_check)
{
	double div, pen_gap, pen_skip, tmp;
	if(is_accurate) {
		(*quick_check) = 1; (*max_skip) = 25; (*max_iter) = 5000; (*max_dis) = 5000; div = 0.01; pen_gap = 0.5f; pen_skip = 0.0005f; 
	} else {
		(*quick_check) = 0; (*max_skip) = 25; (*max_iter) = 5000; (*max_dis) = 5000; div = 0.1; pen_gap = 0.5f; pen_skip = 0.0005f;
	}
	tmp = expf(-div * (double)mz_k);///0.60049557881 -> HiFi; 0.18268352405 -> ont
	*chn_pen_gap = pen_gap * tmp;///0.300247789405 -> HiFi; 0.091341762025 -> ont
	///0.000300247789405 -> HiFi (>3330 will be negative);
	//0.000091341762025 -> ont (>10947 will be negative);
	*chn_pen_skip = pen_skip * tmp; 
}

void ul_map_lchain(ha_abufl_t *ab, uint32_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k, const ul_idx_t *uref, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, 
								 int max_n_chain, int apend_be, kvec_t_u8_warp* k_flag, overlap_region* f_cigar, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t is_accurate, uint32_t gen_off, double mcopy_rate, uint32_t chain_cutoff, uint32_t mcopy_khit_cut)
{
	extern void *ha_flt_tab;
	extern ha_pt_t *ha_idx;
	int64_t max_skip, max_iter, max_dis, quick_check; double chn_pen_gap, chn_pen_skip;
	set_lchain_dp_op(is_accurate, mz_k, &max_skip, &max_iter, &max_dis, &chn_pen_gap, &chn_pen_skip, &quick_check);
	// minimizers_gen(ab, rs, rl, mz_w, mz_k, cl, k_flag, ha_flt_tab, ha_idx, dbg_ct, sp, high_occ, low_occ);
	minimizers_qgen(ab, rs, rl, mz_w, mz_k, cl, k_flag, ha_flt_tab, ha_idx, NULL, uref, dbg_ct, sp, high_occ, low_occ);
	// lchain_gen(cl, overlap_list, rid, rl, NULL, uref, apend_be, f_cigar, max_n_chain, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_thres, quick_check, gen_off);
	// lchain_qgen(cl, overlap_list, rid, rl, NULL, uref, apend_be, f_cigar, max_n_chain, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_thres, quick_check, gen_off);
	///no need to sort here, overlap_list has been sorted at lchain_gen
	lchain_qgen_mcopy(cl, overlap_list, rid, rl, NULL, uref, apend_be, max_n_chain, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_thres, quick_check, gen_off, mcopy_rate, chain_cutoff, mcopy_khit_cut, sp);
}

void h_ec_lchain(ha_abuf_t *ab, uint32_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k, All_reads *rref, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, 
								 int max_n_chain, int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t is_accurate, uint32_t gen_off, int64_t mcopy_num, double mcopy_rate, uint32_t chain_cutoff, uint32_t mcopy_khit_cut, uint64_t ocv_w)
{
	extern void *ha_flt_tab;
	extern ha_pt_t *ha_idx;
	int64_t max_skip, max_iter, max_dis, quick_check; double chn_pen_gap, chn_pen_skip;
	set_lchain_dp_op(is_accurate, mz_k, &max_skip, &max_iter, &max_dis, &chn_pen_gap, &chn_pen_skip, &quick_check);
	// minimizers_gen(ab, rs, rl, mz_w, mz_k, cl, k_flag, ha_flt_tab, ha_idx, dbg_ct, sp, high_occ, low_occ);
	minimizers_qgen0(ab, rs, rl, mz_w, mz_k, cl, k_flag, ha_flt_tab, ha_idx, rref, dbg_ct, sp, high_occ, low_occ);
	// lchain_gen(cl, overlap_list, rid, rl, NULL, uref, apend_be, f_cigar, max_n_chain, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_thres, quick_check, gen_off);
	// lchain_qgen(cl, overlap_list, rid, rl, NULL, uref, apend_be, f_cigar, max_n_chain, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_thres, quick_check, gen_off);
	///no need to sort here, overlap_list has been sorted at lchain_gen
	lchain_qgen_mcopy_fast(cl, overlap_list, rid, rl, rref, apend_be, max_n_chain, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_thres, quick_check, gen_off, mcopy_num, mcopy_rate, chain_cutoff, mcopy_khit_cut, sp, ocv_w);
}

void h_ec_lchain_amz(ha_abuf_t *ab, uint32_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k, All_reads *rref, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, 
								 int max_n_chain, int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t is_accurate, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t chain_cutoff, uint32_t mcopy_khit_cut, uint64_t ocv_w)
{
	extern void *ha_flt_tab;
	extern ha_pt_t *ha_idx;
	int64_t max_skip, max_iter, max_dis, quick_check; double chn_pen_gap, chn_pen_skip;
	set_lchain_dp_op(is_accurate, mz_k, &max_skip, &max_iter, &max_dis, &chn_pen_gap, &chn_pen_skip, &quick_check);
	// minimizers_gen(ab, rs, rl, mz_w, mz_k, cl, k_flag, ha_flt_tab, ha_idx, dbg_ct, sp, high_occ, low_occ);
	minimizers_qgen0_amz(ab, rs, rl, mz_w, mz_k, cl, k_flag, ha_flt_tab, ha_idx, rref, dbg_ct, sp, high_occ, low_occ);
	// lchain_gen(cl, overlap_list, rid, rl, NULL, uref, apend_be, f_cigar, max_n_chain, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_thres, quick_check, gen_off);
	// lchain_qgen(cl, overlap_list, rid, rl, NULL, uref, apend_be, f_cigar, max_n_chain, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_thres, quick_check, gen_off);
	///no need to sort here, overlap_list has been sorted at lchain_gen
	lchain_qgen_mcopy_fast(cl, overlap_list, rid, rl, rref, apend_be, max_n_chain, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_thres, quick_check, gen_off, enable_mcopy, mcopy_rate, chain_cutoff, mcopy_khit_cut, sp, ocv_w);
}

uint64_t recalu_minimizer0(char *s, uint64_t len, uint64_t is_hpc, int64_t mz_k, uint64_t mz_h, tiny_queue_t *tq, uint64_t *rpos, uint64_t *rspan)
{
	uint64_t k, l, shift1 = mz_k - 1, mask = (1ULL<<mz_k) - 1, kmer[4] = {0,0,0,0}, c, z, hs; int64_t mz_l, mz_span;
	(*rpos) = (*rspan) = (uint64_t)-1;
	if(is_hpc) {
		tq->front = tq->count = 0;
		for (k = 1, l = 0, mz_l = mz_span = 0; k <= len; ++k) {
        	if (k == len || seq_nt4_table[(uint8_t)s[k]] != seq_nt4_table[(uint8_t)s[l]]) {
				c = seq_nt4_table[(uint8_t)s[l]];
				if(c < 4) {
					kmer[0] = (kmer[0] << 1 | (c&1))  & mask;/**forward k-mer**/
					kmer[1] = (kmer[1] << 1 | (c>>1)) & mask;
					kmer[2] = kmer[2] >> 1 | (uint64_t)(1 - (c&1))  << shift1; /**reverse k-mer**/
					kmer[3] = kmer[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift1;
					if (kmer[1] == kmer[3]) continue; /** skip "symmetric k-mers" as we don't know it strand**/
					z = kmer[1] < kmer[3]? 0 : 1; /** strand**/
					mz_l++; mz_span += (k - l);
					tq_push(tq, k - l); if (tq->count > mz_k) mz_span -= tq_shift(tq);
					if(mz_l >= mz_k && mz_span < 256) {
						hs = yak_hash64_64(kmer[z<<1|0]) + yak_hash64_64(kmer[z<<1|1]);
						if(mz_h == hs) {
							(*rpos) = k - 1; (*rspan) = mz_span;
							return 1;
						}
					}
				} else {
					mz_l = mz_span = 0;
				}

				l = k;
			}
		}
	} else {
		for (k = 0, mz_l = mz_span = 0; k < len; ++k) {
			c = seq_nt4_table[(uint8_t)s[k]];
			if(c < 4) {
				kmer[0] = (kmer[0] << 1 | (c&1))  & mask;/**forward k-mer**/
				kmer[1] = (kmer[1] << 1 | (c>>1)) & mask;
				kmer[2] = kmer[2] >> 1 | (uint64_t)(1 - (c&1))  << shift1; /**reverse k-mer**/
				kmer[3] = kmer[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift1;
				if (kmer[1] == kmer[3]) continue; /** skip "symmetric k-mers" as we don't know it strand**/
				z = kmer[1] < kmer[3]? 0 : 1; /** strand**/
				mz_l++; mz_span++; if(mz_span > mz_k) mz_span = mz_k;
				if(mz_l >= mz_k) {
					hs = yak_hash64_64(kmer[z<<1|0]) + yak_hash64_64(kmer[z<<1|1]);
					if(mz_h == hs) {
						(*rpos) = k; (*rspan) = mz_span;
						return 1;
					}
				}
			} else {
				mz_l = mz_span = 0;
			}
		}
	}
	
	return 0;
}

uint64_t recalu_minimizer0_adv(char *s, uint64_t len, uint64_t is_hpc, int64_t mz_k, uint64_t mz_h, uint64_t rev, tiny_queue_t *tq, uint64_t *rpos, uint64_t *rspan)
{
	uint64_t k, l, shift1 = mz_k - 1, mask = (1ULL<<mz_k) - 1, kmer[4] = {0,0,0,0}, c, z, hs; int64_t mz_l, mz_span;
	(*rpos) = (*rspan) = (uint64_t)-1;
	if(!rev) {
		if(is_hpc) {
			tq->front = tq->count = 0;
			for (k = 1, l = 0, mz_l = mz_span = 0; k <= len; ++k) {
				if (k == len || seq_nt4_table[(uint8_t)s[k]] != seq_nt4_table[(uint8_t)s[l]]) {
					c = seq_nt4_table[(uint8_t)s[l]];
					if(c < 4) {
						kmer[0] = (kmer[0] << 1 | (c&1))  & mask;/**forward k-mer**/
						kmer[1] = (kmer[1] << 1 | (c>>1)) & mask;
						kmer[2] = kmer[2] >> 1 | (uint64_t)(1 - (c&1))  << shift1; /**reverse k-mer**/
						kmer[3] = kmer[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift1;
						if (kmer[1] == kmer[3]) continue; /** skip "symmetric k-mers" as we don't know it strand**/
						z = kmer[1] < kmer[3]? 0 : 1; /** strand**/
						mz_l++; mz_span += (k - l);
						tq_push(tq, k - l); if (tq->count > mz_k) mz_span -= tq_shift(tq);
						if(mz_l >= mz_k && mz_span < 256) {
							hs = yak_hash64_64(kmer[z<<1|0]) + yak_hash64_64(kmer[z<<1|1]);
							if(mz_h == hs) {
								(*rpos) = k - 1; (*rspan) = mz_span;
								return 1;
							}
						}
					} else {
						mz_l = mz_span = 0;
					}

					l = k;
				}
			}
		} else {
			for (k = 0, mz_l = mz_span = 0; k < len; ++k) {
				c = seq_nt4_table[(uint8_t)s[k]];
				if(c < 4) {
					kmer[0] = (kmer[0] << 1 | (c&1))  & mask;/**forward k-mer**/
					kmer[1] = (kmer[1] << 1 | (c>>1)) & mask;
					kmer[2] = kmer[2] >> 1 | (uint64_t)(1 - (c&1))  << shift1; /**reverse k-mer**/
					kmer[3] = kmer[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift1;
					if (kmer[1] == kmer[3]) continue; /** skip "symmetric k-mers" as we don't know it strand**/
					z = kmer[1] < kmer[3]? 0 : 1; /** strand**/
					mz_l++; mz_span++; if(mz_span > mz_k) mz_span = mz_k;
					if(mz_l >= mz_k) {
						hs = yak_hash64_64(kmer[z<<1|0]) + yak_hash64_64(kmer[z<<1|1]);
						if(mz_h == hs) {
							(*rpos) = k; (*rspan) = mz_span;
							return 1;
						}
					}
				} else {
					mz_l = mz_span = 0;
				}
			}
		}
	} else {
		uint8_t ch[5] = {3, 2, 1, 0, 5};
		if(is_hpc) {
			tq->front = tq->count = 0;
			for (k = 1, l = 0, mz_l = mz_span = 0; k <= len; ++k) {
				if (k == len || seq_nt4_table[(uint8_t)s[len-k-1]] != seq_nt4_table[(uint8_t)s[len-l-1]]) {
					c = ch[seq_nt4_table[(uint8_t)s[len-l-1]]];
					if(c < 4) {
						kmer[0] = (kmer[0] << 1 | (c&1))  & mask;/**forward k-mer**/
						kmer[1] = (kmer[1] << 1 | (c>>1)) & mask;
						kmer[2] = kmer[2] >> 1 | (uint64_t)(1 - (c&1))  << shift1; /**reverse k-mer**/
						kmer[3] = kmer[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift1;
						if (kmer[1] == kmer[3]) continue; /** skip "symmetric k-mers" as we don't know it strand**/
						z = kmer[1] < kmer[3]? 0 : 1; /** strand**/
						mz_l++; mz_span += (k - l);
						tq_push(tq, k - l); if (tq->count > mz_k) mz_span -= tq_shift(tq);
						if(mz_l >= mz_k && mz_span < 256) {
							hs = yak_hash64_64(kmer[z<<1|0]) + yak_hash64_64(kmer[z<<1|1]);
							if(mz_h == hs) {
								(*rpos) = k - 1; (*rspan) = mz_span;
								return 1;
							}
						}
					} else {
						mz_l = mz_span = 0;
					}

					l = k;
				}
			}
		} else {
			for (k = 0, mz_l = mz_span = 0; k < len; ++k) {
				c = ch[seq_nt4_table[(uint8_t)s[len-k-1]]];
				if(c < 4) {
					kmer[0] = (kmer[0] << 1 | (c&1))  & mask;/**forward k-mer**/
					kmer[1] = (kmer[1] << 1 | (c>>1)) & mask;
					kmer[2] = kmer[2] >> 1 | (uint64_t)(1 - (c&1))  << shift1; /**reverse k-mer**/
					kmer[3] = kmer[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift1;
					if (kmer[1] == kmer[3]) continue; /** skip "symmetric k-mers" as we don't know it strand**/
					z = kmer[1] < kmer[3]? 0 : 1; /** strand**/
					mz_l++; mz_span++; if(mz_span > mz_k) mz_span = mz_k;
					if(mz_l >= mz_k) {
						hs = yak_hash64_64(kmer[z<<1|0]) + yak_hash64_64(kmer[z<<1|1]);
						if(mz_h == hs) {
							(*rpos) = k; (*rspan) = mz_span;
							return 1;
						}
					}
				} else {
					mz_l = mz_span = 0;
				}
			}
		}
	}
	
	return 0;
}


uint64_t recalu_minimizer(uint64_t rid, anchor1_t *z, asg16_v *sc, int64_t *iok, int64_t *ink, int64_t *ick, int64_t *str_s, int64_t *str_e, uint64_t mz_k, uint64_t mz_h, tiny_queue_t *tq, All_reads *rref, UC_Read *tu, char *qstr, uint64_t qs, uint64_t qe)
{
	int64_t id = z->srt>>33; char *str = NULL; uint64_t rpos, rspan;
	int64_t ok = *iok, nk = *ink, ck = *ick, cn = sc->n, s0, e0, s1, e1, os, oe, ots, ote, ol, wo[2], wn[2], weo[2], wen[2], ovlp, len = Get_READ_LENGTH((*rref), id); uint16_t op, bq, bt; uint32_t cl;
	e0 = ((uint32_t)z->srt) + 1; s0 = e0 - z->other_off; s1 = e1 = -1; weo[0] = weo[1] = wen[0] = wen[1] = -1;
	if(e0 <= s0) return 0;

	///debug
	// ok = nk = ck = 0;

	if((ck < 0) || (ck > cn)) {//(*ck) == cn is allowed
        ck = ok = nk = 0;
    }

	// if(id == 1/**s0 == 10260 && e0 == 10334**/) {
	// 	fprintf(stderr, "\n\n-0-qk::%ld,\ttk::%ld,\tck::%ld,\ts0::%ld,\te0::%ld\n", ok, nk, ck, s0, e0);
	// }

	while (ck > 0 && ok >= s0) {///x -> t; y -> p; first insertion and then match/mismatch
        --ck;
        op = sc->a[ck]>>14;
        // ol = (((op == 1) || (op == 2))?(sc->a[ck]&(0xfff)):(sc->a[ck]&(0x3fff)));
		if((op == 2) || (op == 3)) {
            ol = sc->a[ck]&(0xfff);
        } else if(op == 1) {
            ol = sc->a[ck]&(0x3ff);
        } else {
            ol = sc->a[ck]&(0x3fff);
        }
        if(op != 2) ok -= ol;
        if(op != 3) nk -= ol;
    }

	// char cm[4]; cm[0] = 'M'; cm[1] = 'S'; cm[2] = 'I'; cm[3] = 'D';  
	// if(id == 1/**s0 == 10260 && e0 == 10334**/) {
	// 	fprintf(stderr, "-1-qk::%ld,\ttk::%ld,\tck::%ld\n", ok, nk, ck);
	// }

	while (ck < cn && ok < e0) { ///[s0, e0)
        wo[0] = ok; wn[0] = nk; 
        // ck = pop_trace_bp(sc, ck, &op, &b, &cl);
		ck = pop_trace_bp_f(sc, ck, &op, &bq, &bt, &cl);
        if(op != 2) ok += cl;
        if(op != 3) nk += cl;
        wo[1] = ok; wn[1] = nk;  

		os = ((s0 >= wo[0])? s0 : wo[0]);
		oe = ((e0 <= wo[1])? e0 : wo[1]);
		// os = MAX(s0, wo[0]); oe = MIN(e0, wo[1]);
        ovlp = ((oe>os)? (oe-os):0);
		// if(id == 1/**s0 == 10260 && e0 == 10334**/) {
		// 	fprintf(stderr, "%u%c(q::[%ld,%ld))(t::[%ld,%ld))(ck::%ld)\n", cl, cm[op], wo[0], wo[1], wn[0], wn[1], ck);
		// }

		if(op != 2) {
			if(!ovlp) continue;
		} else {///wo[0] == wo[1]
			if(wo[0] < s0 || wo[0] >= e0) continue;
		}

		if(op < 2) {
            ots = os - wo[0] + wn[0]; ote = oe - wo[0] + wn[0];
        } else {///op == 2: more y; p == 3: more x
            ots = wn[0]; ote = wn[1];
        } 

		if(s1 == -1) s1 = ots;
		e1 = ote;

        if((op == 0) && (ovlp > 0) && (ovlp > weo[1] - weo[0])) {
			weo[0] = os; weo[1] = oe;
			wen[0] = ots; wen[1] = ote;
        }
    }

	while (ck < cn && ok <= e0) { ///[s0, e0)
        wo[0] = ok; wn[0] = nk; 
        // pop_trace_bp(sc, ck, &op, &b, &cl);
		pop_trace_bp_f(sc, ck, &op, &bq, &bt, &cl);
		if(op != 2) break;
		// ck = pop_trace_bp(sc, ck, &op, &b, &cl);
		ck = pop_trace_bp_f(sc, ck, &op, &bq, &bt, &cl);
        // if(op != 2) ok += cl;
        // if(op != 3) nk += cl;
		nk += cl;
        wo[1] = ok; wn[1] = nk;  

		// if(id == 1/**s0 == 10260 && e0 == 10334**/) {
		// 	fprintf(stderr, "%u%c(q::[%ld,%ld))(t::[%ld,%ld))(ck::%ld)\n", cl, cm[op], wo[0], wo[1], wn[0], wn[1], ck);
		// }

		if(wo[0] >= s0 && wo[0] <= e0) {
			ots = wn[0]; ote = wn[1];
			if(s1 == -1) s1 = ots;
			e1 = ote;
		}
    }

	assert(wen[1] <= len);
	assert(e1 <= len);

	*iok = ok; *ink = nk; *ick = ck; 

	if(weo[0] == s0 && weo[1] == e0) {
		z->srt >>= 32; z->srt <<= 32; z->srt |= ((uint64_t)(wen[1]-1));
		///debug
		// char sstr[256]; recover_UC_Read_sub_region(sstr, wen[0], wen[1] - wen[0], 0, rref, id);
		// if(recalu_minimizer0(sstr, wen[1] - wen[0], !(asm_opt.flag & HA_F_NO_HPC), mz_k, mz_h, tq, &rpos, &rspan) && (rpos + 1 == ((uint64_t)(wen[1] - wen[0]))) && (((uint64_t)(wen[1] - wen[0])) == rspan)) {
		// 	// fprintf(stderr, "-0-[M::%s]\n", __func__);
		// } else {
		// 	// if(((z->srt>>32)&1) == 0) {
		// 		fprintf(stderr, "-1-[M::%s]\trid::%lu\ttid::%ld\t%c\trg1::[%ld,%ld)\trg0::[%ld,%ld)\tspan::%ld\n", __func__, rid, id, "+-"[(z->srt>>32)&1], wen[0], wen[1], s0, e0, e0 - s0);
		// 		fprintf(stderr, "tstr::%.*s\n", ((uint32_t)(wen[1] - wen[0])), sstr);
		// 		fprintf(stderr, "qstr::%.*s\n", ((uint32_t)(qe - qs)), qstr + qs);
		// 		exit(1);
		// 	// }
		// }
		// if(rid == 3196 && id == 3199) fprintf(stderr, "-full-[M::%s]\n", __func__);
		return 1;
	}

	if(e1 <= s1) return 0;

	

	if(s1 >= (*str_s) && e1 <= (*str_e)) {
		str = tu->seq + s1 - (*str_s);
	} else {
		if(s1 >= (*str_s) && (s1 < (*str_e)) && (e1 > (*str_e))) {
			UC_Read_resize((*tu), (e1 - (*str_s))); 
			recover_UC_Read_sub_region(tu->seq + (*str_e) - (*str_s), (*str_e), e1 - (*str_e), 0, rref, id);
			str = tu->seq + s1 - (*str_s); (*str_e) = e1;
		} else {
			UC_Read_resize((*tu), (e1 - s1)); 
			recover_UC_Read_sub_region(tu->seq, s1, e1 - s1, 0, rref, id);
			str = tu->seq; (*str_s) = s1; (*str_e) = e1;
		}
	}

	if(recalu_minimizer0(str, e1 - s1, !(asm_opt.flag & HA_F_NO_HPC), mz_k, mz_h, tq, &rpos, &rspan)) {
		rpos += s1;
		z->srt >>= 32; z->srt <<= 32; z->srt |= rpos; z->other_off = rspan;
		// fprintf(stderr, "-1-[M::%s]\trid::%lu\ttid::%ld\t%c\trg1::[%ld,%ld)\trg0::[%ld,%ld)\tspan::%ld\n", __func__, rid, id, "+-"[(z->srt>>32)&1], s1, e1, s0, e0, e0 - s0);
		// if(rid == 3196 && id == 3199) fprintf(stderr, "-part-[M::%s]\n", __func__);
		return 1;
	}
	// fprintf(stderr, "-0-[M::%s]\trid::%lu\ttid::%ld\t%c\trg1::[%ld,%ld)\trg0::[%ld,%ld)\tspan::%ld\n", __func__, rid, id, "+-"[(z->srt>>32)&1], s1, e1, s0, e0, e0 - s0);

	return 0;
}

int64_t hpc_minimizer_test(int64_t id, All_reads *rref, UC_Read *bu, int64_t str_s, int64_t str_e, int64_t z, int64_t len, int64_t rev, int64_t step)
{
	if(z < 0 || z >= len) return 0;
	if(z == 0 && rev == 1) return 0;
	if(z == len - 1 && rev == 0) return 0;

	int64_t an, sc, ec, k, tn = 0; char c = 0, *a = NULL;
	if(z >= str_s && z < str_e) c = bu->seq[z-str_s];

	if(rev) {
		if(c != 0) {
			if(z > str_s) {
				a = bu->seq; an = z - str_s;
				for (k = an - 1; k >= 0 && a[k] == c; k--);
				tn += an - k - 1;
				if((k >= 0) || (z - tn == 0)) return tn;
				ec = str_s; 
			} else {
				ec = z; 
			}
		} else {
			ec = z + 1; 
		}

		sc = ec - step; if(sc < 0) sc = 0;
		if(ec <= 0 || ec > len || ec <= sc) return tn;

		UC_Read_resize((*bu), (str_e - str_s) + (step)); 
		a = bu->seq + str_e - str_s;

		while (1) {
			recover_UC_Read_sub_region(a, sc, ec - sc, 0, rref, id); an = ec - sc; 
			if(c == 0) c = a[--an];
			for (k = an - 1; k >= 0 && a[k] == c; k--);
			tn += an - k - 1;
			if((k >= 0) || (z - tn == 0)) return tn;
			ec = sc;
			sc = ec - step; if(sc < 0) sc = 0;
			if(ec <= 0 || ec > len || ec <= sc) return tn;
		}
	} else {
		if(c != 0) {
			if(z + 1 < str_e) {
				a = bu->seq + z + 1 - str_s; an = str_e - z - 1;
				for (k = 0; k < an && a[k] == c; k++);
				tn += k;
				if((k < an) || (z + tn + 1 == len)) return tn;
				sc = str_e; 
			} else {
				sc = z + 1; 
			}
		} else {
			sc = z; 
		}

		ec = sc + step; if(ec > len) ec = len;
		if(sc < 0 || sc >= len || sc >= ec) return tn;

		UC_Read_resize((*bu), (str_e - str_s) + (step)); 
		a = bu->seq + str_e - str_s;

		while (1) {
			// if(bu->size < (a - bu->seq) + (ec - sc)) {
			// 	fprintf(stderr, "-1-[M::%s]\tid::%ld\n", __func__, id);
			// }
			recover_UC_Read_sub_region(a, sc, ec - sc, 0, rref, id); an = ec - sc; 
			a = bu->seq + str_e - str_s;
			if(c == 0) {
				c = a[0]; a = a + 1; an--; 
			} 
			for (k = 0; k < an && a[k] == c; k++);
			tn += k;
			if((k < an) || (z + tn + 1 == len)) return tn;
			sc = ec;
			ec = sc + step; if(ec > len) ec = len;
			if(sc < 0 || sc >= len || sc >= ec) return tn;
		}
	}
}

uint64_t recalu_minimizer_bd(uint64_t rid, anchor1_t *z, asg16_v *sc, int64_t *iok, int64_t *ink, int64_t *ick, int64_t *str_s, int64_t *str_e, uint64_t mz_k, uint64_t mz_h, tiny_queue_t *tq, All_reads *rref, UC_Read *tu, char *qstr, uint64_t qs, uint64_t qe)
{
	int64_t id = z->srt>>33; char *str = NULL; uint64_t rpos, rspan;
	int64_t ok = *iok, nk = *ink, ck = *ick, cn = sc->n, s0, e0, s1, e1, si, ei, os, oe, ots, ote, ol, wo[2], wn[2], weo[2], wen[2], ovlp, len = Get_READ_LENGTH((*rref), id); uint16_t op, bq, bt; uint32_t cl;
	e0 = ((uint32_t)z->srt) + 1; s0 = e0 - z->other_off; s1 = e1 = si = ei = -1; weo[0] = weo[1] = wen[0] = wen[1] = -1; 
	if(e0 <= s0) return 0;

	///debug
	// ok = nk = ck = 0;

	if((ck < 0) || (ck > cn)) {//(*ck) == cn is allowed
        ck = ok = nk = 0;
    }

	// if(id == 1/**s0 == 10260 && e0 == 10334**/) {
	// 	fprintf(stderr, "\n\n-0-qk::%ld,\ttk::%ld,\tck::%ld,\ts0::%ld,\te0::%ld\n", ok, nk, ck, s0, e0);
	// }

	while (ck > 0 && ok >= s0) {///x -> t; y -> p; first insertion and then match/mismatch
        --ck;
        op = sc->a[ck]>>14;
        // ol = (((op == 1) || (op == 2))?(sc->a[ck]&(0xfff)):(sc->a[ck]&(0x3fff)));
		if((op == 2) || (op == 3)) {
            ol = sc->a[ck]&(0xfff);
        } else if(op == 1) {
            ol = sc->a[ck]&(0x3ff);
        } else {
            ol = sc->a[ck]&(0x3fff);
        }
        if(op != 2) ok -= ol;
        if(op != 3) nk -= ol;
    }

	// char cm[4]; cm[0] = 'M'; cm[1] = 'S'; cm[2] = 'I'; cm[3] = 'D';  
	// if(id == 1/**s0 == 10260 && e0 == 10334**/) {
	// 	fprintf(stderr, "-1-qk::%ld,\ttk::%ld,\tck::%ld\n", ok, nk, ck);
	// }
	while (ck < cn && ok < e0) { ///[s0, e0)
        wo[0] = ok; wn[0] = nk; 
        // ck = pop_trace_bp(sc, ck, &op, &b, &cl);
		ck = pop_trace_bp_f(sc, ck, &op, &bq, &bt, &cl);
        if(op != 2) ok += cl;
        if(op != 3) nk += cl;
        wo[1] = ok; wn[1] = nk;  

		os = ((s0 >= wo[0])? s0 : wo[0]);
		oe = ((e0 <= wo[1])? e0 : wo[1]);
		// os = MAX(s0, wo[0]); oe = MIN(e0, wo[1]);
        ovlp = ((oe>os)? (oe-os):0);
		// if(id == 1/**s0 == 10260 && e0 == 10334**/) {
		// 	fprintf(stderr, "%u%c(q::[%ld,%ld))(t::[%ld,%ld))(ck::%ld)\n", cl, cm[op], wo[0], wo[1], wn[0], wn[1], ck);
		// }

		if(op != 2) {
			if(!ovlp) continue;
		} else {///wo[0] == wo[1]
			if(wo[0] < s0 || wo[0] >= e0) continue;
		}

		if(op < 2) {
            ots = os - wo[0] + wn[0]; ote = oe - wo[0] + wn[0];
        } else {///op == 2: more y; p == 3: more x
            ots = wn[0]; ote = wn[1];
        } 

		if(s1 == -1) {
			s1 = ots;
			// if((op == 0) && (wn[1] > ots) && ((wn[0] < ots) || (ots == 0))) si = 1;
			if((op == 0) && (wo[1] > s0) && ((wo[0] < s0) || (s1 == 0))) si = 1;
		}
		
		e1 = ote;
		if((op == 0) && (wo[0] < e0) && ((wo[1] > e0) || (e1 == len))) ei = 1;

        if((op == 0) && (ovlp > 0) && (ovlp > weo[1] - weo[0])) {
			weo[0] = os; weo[1] = oe;
			wen[0] = ots; wen[1] = ote;
        }
    }

	while (ck < cn && ok <= e0) { ///[s0, e0)
        wo[0] = ok; wn[0] = nk; 
        // pop_trace_bp(sc, ck, &op, &b, &cl);
		pop_trace_bp_f(sc, ck, &op, &bq, &bt, &cl);
		if(op != 2) break;
		// ck = pop_trace_bp(sc, ck, &op, &b, &cl);
		ck = pop_trace_bp_f(sc, ck, &op, &bq, &bt, &cl);
        // if(op != 2) ok += cl;
        // if(op != 3) nk += cl;
		nk += cl;
        wo[1] = ok; wn[1] = nk;  

		// if(id == 1/**s0 == 10260 && e0 == 10334**/) {
		// 	fprintf(stderr, "%u%c(q::[%ld,%ld))(t::[%ld,%ld))(ck::%ld)\n", cl, cm[op], wo[0], wo[1], wn[0], wn[1], ck);
		// }

		if(wo[0] >= s0 && wo[0] <= e0) {
			ots = wn[0]; ote = wn[1];
			if(s1 == -1) s1 = ots;
			e1 = ote;
		}
    }

	assert(wen[1] <= len);
	assert(e1 <= len);

	*iok = ok; *ink = nk; *ick = ck; 
	if(e1 <= s1) return 0;
	rpos = rspan = ((uint64_t)-1);

	if(weo[0] == s0 && weo[1] == e0) {
		rpos = wen[1]-1; rspan = e0 - s0;
		/** 
		if(si == 1 && ei == 1) {
			z->srt >>= 32; z->srt <<= 32; z->srt |= ((uint64_t)(wen[1]-1));
			///debug
			// char sstr[256]; recover_UC_Read_sub_region(sstr, wen[0], wen[1] - wen[0], 0, rref, id);
			// if(recalu_minimizer0(sstr, wen[1] - wen[0], !(asm_opt.flag & HA_F_NO_HPC), mz_k, mz_h, tq, &rpos, &rspan) && (rpos + 1 == ((uint64_t)(wen[1] - wen[0]))) && (((uint64_t)(wen[1] - wen[0])) == rspan)) {
			// 	// fprintf(stderr, "-0-[M::%s]\n", __func__);
			// } else {
			// 	// if(((z->srt>>32)&1) == 0) {
			// 		fprintf(stderr, "-1-[M::%s]\trid::%lu\ttid::%ld\t%c\trg1::[%ld,%ld)\trg0::[%ld,%ld)\tspan::%ld\n", __func__, rid, id, "+-"[(z->srt>>32)&1], wen[0], wen[1], s0, e0, e0 - s0);
			// 		fprintf(stderr, "tstr::%.*s\n", ((uint32_t)(wen[1] - wen[0])), sstr);
			// 		fprintf(stderr, "qstr::%.*s\n", ((uint32_t)(qe - qs)), qstr + qs);
			// 		exit(1);
			// 	// }
			// }
			// if(rid == 3196 && id == 3199) fprintf(stderr, "-full-[M::%s]\n", __func__);
			return 1;
		} 
		**/
	} else {
		if(s1 >= (*str_s) && e1 <= (*str_e)) {
			str = tu->seq + s1 - (*str_s);
		} else {
			if(s1 >= (*str_s) && (s1 < (*str_e)) && (e1 > (*str_e))) {
				UC_Read_resize((*tu), (e1 - (*str_s))); 
				recover_UC_Read_sub_region(tu->seq + (*str_e) - (*str_s), (*str_e), e1 - (*str_e), 0, rref, id);
				str = tu->seq + s1 - (*str_s); (*str_e) = e1;
			} else {
				UC_Read_resize((*tu), (e1 - s1)); 
				recover_UC_Read_sub_region(tu->seq, s1, e1 - s1, 0, rref, id);
				str = tu->seq; (*str_s) = s1; (*str_e) = e1;
			}
		}

		if(recalu_minimizer0(str, e1 - s1, !(asm_opt.flag & HA_F_NO_HPC), mz_k, mz_h, tq, &rpos, &rspan)) {
			rpos += s1;
		}
	}

	if(rpos != ((uint64_t)-1) && rspan != ((uint64_t)-1)) {
        e0 = rpos + 1; s0 = e0 - rspan;
        if(!(asm_opt.flag & HA_F_NO_HPC)) {
			if(si == -1) {
				s0 -= hpc_minimizer_test(id, rref, tu, *str_s, *str_e, s0, len, 1, 8);
				assert(s0 >= 0 && s0 < len && s0 < e0);
			}

			if(ei == -1) {
				e0 += hpc_minimizer_test(id, rref, tu, *str_s, *str_e, e0-1, len, 0, 8);
				assert(e0 >= 0 && e0 < len && s0 < e0);
			}
        }
        
        rpos = e0 - 1; rspan = e0 - s0;
        if(rspan < 256) {
            z->srt >>= 32; z->srt <<= 32; z->srt |= rpos; z->other_off = rspan;
            return 1;
        }
    }
	/**
	if(recalu_minimizer0(str, e1 - s1, !(asm_opt.flag & HA_F_NO_HPC), mz_k, mz_h, tq, &rpos, &rspan)) {
		rpos += s1;
		z->srt >>= 32; z->srt <<= 32; z->srt |= rpos; z->other_off = rspan;
		// fprintf(stderr, "-1-[M::%s]\trid::%lu\ttid::%ld\t%c\trg1::[%ld,%ld)\trg0::[%ld,%ld)\tspan::%ld\n", __func__, rid, id, "+-"[(z->srt>>32)&1], s1, e1, s0, e0, e0 - s0);
		// if(rid == 3196 && id == 3199) fprintf(stderr, "-part-[M::%s]\n", __func__);
		return 1;
	}
	// fprintf(stderr, "-0-[M::%s]\trid::%lu\ttid::%ld\t%c\trg1::[%ld,%ld)\trg0::[%ld,%ld)\tspan::%ld\n", __func__, rid, id, "+-"[(z->srt>>32)&1], s1, e1, s0, e0, e0 - s0);
	**/

	return 0;
}


uint64_t recalu_minimizer_non_retrieve(uint64_t rid, anchor1_t *z, asg16_v *sc, int64_t *iok, int64_t *ink, int64_t *ick, uint64_t mz_k, uint64_t mz_h, tiny_queue_t *tq, char *tstr, int64_t tl, uint64_t trev/** , char *qstr, uint64_t qs, uint64_t qe**/)
{
	char *str = NULL; uint64_t rpos, rspan; char c;
	int64_t ok = *iok, nk = *ink, ck = *ick, cn = sc->n, s0, e0, s1, e1, os, oe, ots, ote, ol, wo[2], wn[2], weo[2], wen[2], ovlp, k; uint16_t op, bq, bt; uint32_t cl;
	e0 = ((uint32_t)z->srt) + 1; s0 = e0 - z->other_off; s1 = e1 = -1; weo[0] = weo[1] = wen[0] = wen[1] = -1;
	if(e0 <= s0) return 0;


	if((ck < 0) || (ck > cn)) {//(*ck) == cn is allowed
        ck = ok = nk = 0;
    }

	while (ck > 0 && ok >= s0) {///x -> t; y -> p; first insertion and then match/mismatch
        --ck;
        op = sc->a[ck]>>14;
        // ol = (((op == 1) || (op == 2))?(sc->a[ck]&(0xfff)):(sc->a[ck]&(0x3fff)));
		if((op == 2) || (op == 3)) {
            ol = sc->a[ck]&(0xfff);
        } else if(op == 1) {
            ol = sc->a[ck]&(0x3ff);
        } else {
            ol = sc->a[ck]&(0x3fff);
        }
        if(op != 2) ok -= ol;
        if(op != 3) nk -= ol;
    }

	while (ck < cn && ok < e0) { ///[s0, e0)
        wo[0] = ok; wn[0] = nk; 
        // ck = pop_trace_bp(sc, ck, &op, &b, &cl);
		ck = pop_trace_bp_f(sc, ck, &op, &bq, &bt, &cl);
        if(op != 2) ok += cl;
        if(op != 3) nk += cl;
        wo[1] = ok; wn[1] = nk;  

		os = ((s0 >= wo[0])? s0 : wo[0]);
		oe = ((e0 <= wo[1])? e0 : wo[1]);
		// os = MAX(s0, wo[0]); oe = MIN(e0, wo[1]);
        ovlp = ((oe>os)? (oe-os):0);

		if(op != 2) {
			if(!ovlp) continue;
		} else {///wo[0] == wo[1]
			if(wo[0] < s0 || wo[0] >= e0) continue;
		}

		if(op < 2) {
            ots = os - wo[0] + wn[0]; ote = oe - wo[0] + wn[0];
        } else {///op == 2: more y; p == 3: more x
            ots = wn[0]; ote = wn[1];
        } 

		if(s1 == -1) s1 = ots;
		e1 = ote;

        if((op == 0) && (ovlp > 0) && (ovlp > weo[1] - weo[0])) {
			weo[0] = os; weo[1] = oe;
			wen[0] = ots; wen[1] = ote;
        }
    }

	while (ck < cn && ok <= e0) { ///[s0, e0)
        wo[0] = ok; wn[0] = nk; 
        // pop_trace_bp(sc, ck, &op, &b, &cl);
		pop_trace_bp_f(sc, ck, &op, &bq, &bt, &cl);
		if(op != 2) break;
		// ck = pop_trace_bp(sc, ck, &op, &b, &cl);
		ck = pop_trace_bp_f(sc, ck, &op, &bq, &bt, &cl);
        // if(op != 2) ok += cl;
        // if(op != 3) nk += cl;
		nk += cl;
        wo[1] = ok; wn[1] = nk;  


		if(wo[0] >= s0 && wo[0] <= e0) {
			ots = wn[0]; ote = wn[1];
			if(s1 == -1) s1 = ots;
			e1 = ote;
		}
    }

	assert(wen[1] <= tl);
	assert(e1 <= tl);

	*iok = ok; *ink = nk; *ick = ck; 

	rpos = rspan = ((uint64_t)-1);
	if(weo[0] == s0 && weo[1] == e0) {
		rpos = wen[1]-1; rspan = e0 - s0;
		// z->srt >>= 32; z->srt <<= 32; z->srt |= ((uint64_t)(wen[1]-1));
		///debug
		// char sstr[256]; recover_UC_Read_sub_region(sstr, wen[0], wen[1] - wen[0], 0, rref, id);
		// if(recalu_minimizer0(sstr, wen[1] - wen[0], !(asm_opt.flag & HA_F_NO_HPC), mz_k, mz_h, tq, &rpos, &rspan) && (rpos + 1 == ((uint64_t)(wen[1] - wen[0]))) && (((uint64_t)(wen[1] - wen[0])) == rspan)) {
		// 	// fprintf(stderr, "-0-[M::%s]\n", __func__);
		// } else {
		// 	// if(((z->srt>>32)&1) == 0) {
		// 		fprintf(stderr, "-1-[M::%s]\trid::%lu\ttid::%ld\t%c\trg1::[%ld,%ld)\trg0::[%ld,%ld)\tspan::%ld\n", __func__, rid, id, "+-"[(z->srt>>32)&1], wen[0], wen[1], s0, e0, e0 - s0);
		// 		fprintf(stderr, "tstr::%.*s\n", ((uint32_t)(wen[1] - wen[0])), sstr);
		// 		fprintf(stderr, "qstr::%.*s\n", ((uint32_t)(qe - qs)), qstr + qs);
		// 		exit(1);
		// 	// }
		// }
		// if(rid == 3196 && id == 3199) fprintf(stderr, "-full-[M::%s]\n", __func__);
		// return 1;
	} else {
		if(e1 <= s1) return 0;
		// if(s1 >= tl || e1 >= tl) return 0;
		
		if(trev) str = tstr + tl - e1;
		else str = tstr + s1; 

		if(recalu_minimizer0_adv(str, e1 - s1, !(asm_opt.flag & HA_F_NO_HPC), mz_k, mz_h, trev, tq, &rpos, &rspan)) {
			rpos += s1;
			// z->srt >>= 32; z->srt <<= 32; z->srt |= rpos; z->other_off = rspan;
			// fprintf(stderr, "-1-[M::%s]\trid::%lu\ttid::%ld\t%c\trg1::[%ld,%ld)\trg0::[%ld,%ld)\tspan::%ld\n", __func__, rid, id, "+-"[(z->srt>>32)&1], s1, e1, s0, e0, e0 - s0);
			// if(rid == 3196 && id == 3199) fprintf(stderr, "-part-[M::%s]\n", __func__);
			// return 1;
		}
	}

	if(rpos != ((uint64_t)-1) && rspan != ((uint64_t)-1)) {
		e0 = rpos + 1; s0 = e0 - rspan;
		if(!(asm_opt.flag & HA_F_NO_HPC)) {
			if(!trev) {
				s1 = s0; e1 = e0; 
			} else {
				s1 = tl - e0; e1 = tl - s0;
			}

			c = tstr[s1];
			for (k = s1 - 1; k >= 0 && tstr[k] == c; k--); s1 = k + 1;

			c = tstr[e1-1];
			for (k = e1; k < tl && tstr[k] == c; k++); e1 = k;

			if(!trev) {
				s0 = s1; e0 = e1; 
			} else {
				s0 = tl - e1; e0 = tl - s1;
			}
		}
		
		rpos = e0 - 1; rspan = e0 - s0;
		if(rspan < 256) {
			z->srt >>= 32; z->srt <<= 32; z->srt |= rpos; z->other_off = rspan;
			return 1;
		}
	}
	// fprintf(stderr, "-0-[M::%s]\trid::%lu\ttid::%ld\t%c\trg1::[%ld,%ld)\trg0::[%ld,%ld)\tspan::%ld\n", __func__, rid, id, "+-"[(z->srt>>32)&1], s1, e1, s0, e0, e0 - s0);

	return 0;
}

void hpc_ext_check(All_reads *rref, uint32_t id, int64_t s0, int64_t e0, int64_t l, uint64_t rev, char *buf)
{
	if(s0 < 0) s0 = 0; if(e0 > l) e0 = l;
	int64_t s = s0 - 256, e = e0 + 256, n, k, os0, oe0, os1, oe1; char c; if(s < 0) s = 0; if(e > l) e = l;
	recover_UC_Read_sub_region(buf, s, e - s, rev, rref, id); 
	os0 = s0 - s; oe0 = e0 - s; n = e - s; os1 = os0; oe1 = oe0;
	c = buf[os0];
	for (k = os0 - 1; k >= 0 && buf[k] == c; k--); os1 = k + 1;

	c = buf[oe0-1];
	for (k = oe0; k < n && buf[k] == c; k++); oe1 = k;

	if(os1 != os0 || oe1 != oe0) {
		fprintf(stderr, "[M::%s]\to0::[%ld,%ld)\to1::[%ld,%ld)\n", __func__, os0 + s0, oe0 + s0, os1 + s0, oe1 + s0);
	}

}

void h_ec_lchain_re_gen(ha_abuf_t *ab, uint32_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k, ha_pt_t *ha_idx, All_reads *rref, overlap_region_alloc *olst, Candidates_list *cl, double bw_thres, 
								 int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, 
								 int64_t max_skip, int64_t max_iter, int64_t max_dis, int64_t quick_check, double chn_pen_gap, double chn_pen_skip, UC_Read *tu, asg64_v *oidx, asg16_v *scc)
{
	uint64_t i, k, l, m, max_cnt = UINT32_MAX, min_cnt = 0; int n, n0, j; ha_mz1_t *z; seed1_t *s; tiny_queue_t tq; memset(&tq, 0, sizeof(tiny_queue_t));
    if(high_occ) {
        max_cnt = (*high_occ);
        if(max_cnt < 2) max_cnt = 2;
    }
    if(low_occ) {
        min_cnt = (*low_occ);
        if(min_cnt < 2) min_cnt = 2;
    }
    clear_Candidates_list(cl); ab->n_a = 0;

	// minimizer of queried read
    if (ab->mz.m > ab->old_mz_m) {
        ab->old_mz_m = ab->mz.m;
        REALLOC(ab->seed, ab->old_mz_m);
    }

    for (i = 0, ab->n_a = 0; i < ab->mz.n; ++i) {
        ab->seed[i].a = ha_pt_get(ha_idx, ab->mz.a[i].x, &n);
        ab->seed[i].n = n;
        ab->n_a += n;
    }

	if (ab->n_a > ab->m_a) {
        ab->m_a = ab->n_a;
        REALLOC(ab->a, ab->m_a);
    }

	for (i = 0, k = 0; i < ab->mz.n; ++i) {
        ///z is one of the minimizer
        z = &ab->mz.a[i]; s = &ab->seed[i];

		// uint64_t rpos, rspan;
		// if(!recalu_minimizer0(rs + (z->pos+1-z->span), z->span, !(asm_opt.flag & HA_F_NO_HPC), mz_k, z->x, &tq, &rpos, &rspan)) {
		// 	fprintf(stderr, "-1-[M::%s]\t%c\trg0::[%u,%u)\trid::%u\n", __func__, "+-"[z->rev], z->pos+1-z->span, z->pos+1, rid);
		// } 
		// else {
		// 	fprintf(stderr, "-0-[M::%s]\t%c\trg0::[%u,%u)\trid::%u\n", __func__, "+-"[z->rev], z->pos+1-z->span, z->pos+1, rid);
		// }

        for (j = 0; j < s->n; ++j) {
            const ha_idxpos_t *y = &s->a[j];
            anchor1_t *an = &ab->a[k++];
            uint8_t rev = z->rev == y->rev? 0 : 1;

            // an->other_off = rev?((uint32_t)-1)-1-(y->pos+1-y->span):y->pos;
			an->other_off = y->span;

            // an->self_off = z->pos;
			an->self_off = i;
			
            ///an->cnt: cnt<<8|span
            an->cnt = s->n; if(an->cnt > ((uint32_t)(0xffffffu))) an->cnt = 0xffffffu;
            an->cnt <<= 8; an->cnt |= ((z->span <= ((uint32_t)(0xffu)))?z->span:((uint32_t)(0xffu)));

            // an->srt = (uint64_t)y->rid<<33 | (uint64_t)rev<<32 | an->self_off;
			an->srt = (uint64_t)y->rid<<33 | (uint64_t)rev<<32 | y->pos;
        }
    }

	// copy over to _cl_
    if (ab->m_a >= (uint64_t)cl->size) {
        cl->size = ab->m_a;
        REALLOC(cl->list, cl->size);
    }
	clear_overlap_region_alloc(olst);

	// char dbg[256]; uint64_t rpos, rspan, thash;
	// char dbg[768];
	k_mer_hit *p; uint64_t tid = (uint64_t)-1, tl = (uint64_t)-1, trev, tspan, ol, zn, olst_n; int64_t ok, nk, ck, str_s, str_e;
    radix_sort_ha_an1(ab->a, ab->a + ab->n_a);
    for (k = 1, l = ol = n = zn = 0; k <= ab->n_a; ++k) {
        if (k == ab->n_a || (ab->a[k].srt>>32) != (ab->a[l].srt>>32)) {
			for (; (ol < oidx->n) && ((oidx->a[ol]>>32) < (ab->a[l].srt>>32)); ol++);
			if((ol < oidx->n) && ((ab->a[l].srt>>32) == (oidx->a[ol]>>32))) {
				tl = Get_READ_LENGTH((*rref), ab->a[l].srt>>33); tid = ab->a[l].srt>>33; trev = (ab->a[l].srt>>32) & 1; olst_n = olst->length;
				for (i = l, m = 0, ok = nk = ck = str_s = str_e = 0; i < k; i++) {
					if(!recalu_minimizer(rid, &(ab->a[i]), &(scc[tid]), &ok, &nk, &ck, &str_s, &str_e, mz_k, ab->mz.a[ab->a[i].self_off].x, &tq, rref, tu, rs, ab->mz.a[ab->a[i].self_off].pos+1-ab->mz.a[ab->a[i].self_off].span, ab->mz.a[ab->a[i].self_off].pos+1)) continue;
					///debug
					// recover_UC_Read_sub_region(dbg, ((uint32_t)ab->a[i].srt) + 1 - ab->a[i].other_off, ab->a[i].other_off, 0, rref, tid);
					// if(recalu_minimizer0(dbg, ab->a[i].other_off, !(asm_opt.flag & HA_F_NO_HPC), mz_k, ab->mz.a[ab->a[i].self_off].x, &tq, &rpos, &rspan) && (rpos + 1 == ab->a[i].other_off) && (ab->a[i].other_off == rspan)) {
					// 	// fprintf(stderr, "-0-[M::%s]\n", __func__);
					// } else {
					// 	fprintf(stderr, "-1-[M::%s]\n", __func__);
					// }
					// thash = ab->mz.a[ab->a[i].self_off].x;
					ab->a[m] = ab->a[i]; tspan = ab->a[m].other_off;
					ab->a[m].other_off = (uint32_t)ab->a[m].srt;
					if(trev) ab->a[m].other_off = tl - (ab->a[m].other_off+1-tspan) - 1;///looks like a bug
					ab->a[m].self_off = ab->mz.a[ab->a[m].self_off].pos;
					ab->a[m].srt = ab->a[m].self_off; ab->a[m].srt <<= 32; ab->a[m].srt |= ab->a[m].other_off;
					///debug
					// recover_UC_Read_sub_region(dbg, ab->a[m].other_off + 1 - tspan, tspan, trev, rref, tid);
					// if(recalu_minimizer0(dbg, tspan, !(asm_opt.flag & HA_F_NO_HPC), mz_k, thash, &tq, &rpos, &rspan) && (rpos + 1 == tspan) && (tspan == rspan)) {
					// 	fprintf(stderr, "-0-[M::%s]\n", __func__);
					// } else {
					// 	fprintf(stderr, "-1-[M::%s]\ttrev::%lu\n", __func__, trev);
					// }
					///debug
					// if(rid == 3196 && tid == 3199) hpc_ext_check(rref, tid, ab->a[m].other_off + 1 - tspan, ab->a[m].other_off + 1, tl, trev, dbg);


					m++;
				}
				if(m > 1) radix_sort_ha_an1(ab->a, ab->a + m);
				for (i = 0, n0 = n; i < m; i++) {
					p = &cl->list[n++];
					p->readID = tid;
                	p->strand = trev;
                    p->offset = ab->a[i].other_off;
                	p->self_offset = ab->a[i].self_off;
					if(((ab->a[i].cnt>>8) < max_cnt) && ((ab->a[i].cnt>>8) > min_cnt)){
						p->cnt = 1;
					} else if((ab->a[i].cnt>>8) <= min_cnt) {
						p->cnt = 2;
					} else{
						p->cnt = 1 + (((ab->a[i].cnt>>8) + (max_cnt<<1) - 1)/(max_cnt<<1));
						p->cnt = pow(p->cnt, 1.1);
					}
					if(p->cnt > ((uint32_t)(0xffffffu))) p->cnt = 0xffffffu;
					p->cnt <<= 8; p->cnt |= (((uint32_t)(0xffu))&(ab->a[i].cnt));
				}

				// if(rid == 3196 && tid == 3199) fprintf(stderr, "[M::%s]\ttid::%lu\ttrev::%lu\told::%lu\tnew::%lu\n", __func__, tid, trev, k - l, m);

				if(m > 0 && tid != rid) {
					zn += lchain_qdp_mcopy_fast(cl, n0, n-n0, zn, &(cl->chainDP), olst, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_thres, 
                													rid, rl, tl, quick_check, apend_be, gen_off, enable_mcopy, mcopy_rate, mcopy_khit_cut, 1);
				}

				if(olst->length > olst_n) {
					oidx->a[ol] >>= 32; oidx->a[ol] <<= 32; oidx->a[ol] |= ((uint64_t)((uint32_t)-1));
				}
			}
			
            l = k;
        }
    }
    cl->length = zn;

	for (k = m = 0; k < oidx->n; k++) {
		if(((uint32_t)oidx->a[k]) == ((uint32_t)-1)) continue;
		oidx->a[m++] = oidx->a[k];
	}
	// fprintf(stderr, "[M::%s]\ttot::%lu\tremain::%lu\n", __func__, (uint64_t)oidx->n, m);
	oidx->n = m;

	for (i = 0; i < olst->length; ++i) olst->list[i].align_length = 0;

	// minimizers_qgen0(ab, rs, rl, mz_w, mz_k, cl, k_flag, ha_flt_tab, ha_idx, rref, dbg_ct, sp, high_occ, low_occ);
	
	// lchain_qgen_mcopy_fast(cl, overlap_list, rid, rl, rref, apend_be, max_n_chain, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_thres, quick_check, gen_off, enable_mcopy, mcopy_rate, chain_cutoff, mcopy_khit_cut, sp);
}


void h_ec_lchain_re_gen3(ha_abuf_t *ab, uint32_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k, ha_pt_t *ha_idx, All_reads *rref, overlap_region_alloc *olst, Candidates_list *cl, double bw_thres, 
								 int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, 
								 int64_t max_skip, int64_t max_iter, int64_t max_dis, int64_t quick_check, double chn_pen_gap, double chn_pen_skip, UC_Read *tu, asg64_v *oidx, asg16_v *scc)
{
	uint64_t i, k, l, m, max_cnt = UINT32_MAX, min_cnt = 0; int n, n0, j; ha_mz1_t *z; seed1_t *s; tiny_queue_t tq; memset(&tq, 0, sizeof(tiny_queue_t));
    if(high_occ) {
        max_cnt = (*high_occ);
        if(max_cnt < 2) max_cnt = 2;
    }
    if(low_occ) {
        min_cnt = (*low_occ);
        if(min_cnt < 2) min_cnt = 2;
    }
    clear_Candidates_list(cl); ab->n_a = 0;

	// minimizer of queried read
    if (ab->mz.m > ab->old_mz_m) {
        ab->old_mz_m = ab->mz.m;
        REALLOC(ab->seed, ab->old_mz_m);
    }

    for (i = 0, ab->n_a = 0; i < ab->mz.n; ++i) {
        ab->seed[i].a = ha_pt_get(ha_idx, ab->mz.a[i].x, &n);
        ab->seed[i].n = n;
        ab->n_a += n;
    }

	if (ab->n_a > ab->m_a) {
        ab->m_a = ab->n_a;
        REALLOC(ab->a, ab->m_a);
    }

	for (i = 0, k = 0; i < ab->mz.n; ++i) {
        ///z is one of the minimizer
        z = &ab->mz.a[i]; s = &ab->seed[i];

		// uint64_t rpos, rspan;
		// if(!recalu_minimizer0(rs + (z->pos+1-z->span), z->span, !(asm_opt.flag & HA_F_NO_HPC), mz_k, z->x, &tq, &rpos, &rspan)) {
		// 	fprintf(stderr, "-1-[M::%s]\t%c\trg0::[%u,%u)\trid::%u\n", __func__, "+-"[z->rev], z->pos+1-z->span, z->pos+1, rid);
		// } 
		// else {
		// 	fprintf(stderr, "-0-[M::%s]\t%c\trg0::[%u,%u)\trid::%u\n", __func__, "+-"[z->rev], z->pos+1-z->span, z->pos+1, rid);
		// }

        for (j = 0; j < s->n; ++j) {
            const ha_idxpos_t *y = &s->a[j];
            anchor1_t *an = &ab->a[k++];
            uint8_t rev = z->rev == y->rev? 0 : 1;

            // an->other_off = rev?((uint32_t)-1)-1-(y->pos+1-y->span):y->pos;
			an->other_off = y->span;

            // an->self_off = z->pos;
			an->self_off = i;
			
            ///an->cnt: cnt<<8|span
            an->cnt = s->n; if(an->cnt > ((uint32_t)(0xffffffu))) an->cnt = 0xffffffu;
            an->cnt <<= 8; an->cnt |= ((z->span <= ((uint32_t)(0xffu)))?z->span:((uint32_t)(0xffu)));

            // an->srt = (uint64_t)y->rid<<33 | (uint64_t)rev<<32 | an->self_off;
			an->srt = (uint64_t)y->rid<<33 | (uint64_t)rev<<32 | y->pos;
        }
    }

	// copy over to _cl_
    if (ab->m_a >= (uint64_t)cl->size) {
        cl->size = ab->m_a;
        REALLOC(cl->list, cl->size);
    }
	clear_overlap_region_alloc(olst);

	// char dbg[256]; uint64_t rpos, rspan, thash;
	// char dbg[768]; uint64_t rpos, rspan, thash;
	k_mer_hit *p; uint64_t tid = (uint64_t)-1, tl = (uint64_t)-1, trev, tspan, ol, zn, olst_n, sk, sv; int64_t ok, nk, ck, str_s, str_e;
    radix_sort_ha_an1(ab->a, ab->a + ab->n_a);
    for (k = 1, l = ol = n = zn = 0; k <= ab->n_a; ++k) {
        if (k == ab->n_a || (ab->a[k].srt>>32) != (ab->a[l].srt>>32)) {
			
			for (sk = sv = 0; (ol < oidx->n) && ((oidx->a[ol]>>33) < (ab->a[l].srt>>33)); ol++);
			if((ol < oidx->n) && ((ab->a[l].srt>>33) == (oidx->a[ol]>>33))) {
				for (; (ol < oidx->n) && ((oidx->a[ol]>>32) < (ab->a[l].srt>>32)); ol++);
				if((ol < oidx->n) && ((ab->a[l].srt>>32) == (oidx->a[ol]>>32))) {
					if(((uint32_t)oidx->a[ol]) == ((uint32_t)-1)) sk = 1;///exact match
					else sv = 1;
				} else {///match in rev
					sk = 1;
				}
			}

			if(!sk) {
				tl = Get_READ_LENGTH((*rref), ab->a[l].srt>>33); tid = ab->a[l].srt>>33; trev = (ab->a[l].srt>>32) & 1; olst_n = olst->length;
				if(tid != rid) {
					for (i = l, m = 0, ok = nk = ck = str_s = str_e = 0; i < k; i++) {
						if(!recalu_minimizer_bd(rid, &(ab->a[i]), &(scc[tid]), &ok, &nk, &ck, &str_s, &str_e, mz_k, ab->mz.a[ab->a[i].self_off].x, &tq, rref, tu, rs, ab->mz.a[ab->a[i].self_off].pos+1-ab->mz.a[ab->a[i].self_off].span, ab->mz.a[ab->a[i].self_off].pos+1)) continue;
						///debug
						// recover_UC_Read_sub_region(dbg, ((uint32_t)ab->a[i].srt) + 1 - ab->a[i].other_off, ab->a[i].other_off, 0, rref, tid);
						// if(recalu_minimizer0(dbg, ab->a[i].other_off, !(asm_opt.flag & HA_F_NO_HPC), mz_k, ab->mz.a[ab->a[i].self_off].x, &tq, &rpos, &rspan) && (rpos + 1 == ab->a[i].other_off) && (ab->a[i].other_off == rspan)) {
						// 	// fprintf(stderr, "-0-[M::%s]\n", __func__);
						// } else {
						// 	fprintf(stderr, "-1-[M::%s]\n", __func__);
						// }
						// thash = ab->mz.a[ab->a[i].self_off].x;
						ab->a[m] = ab->a[i]; tspan = ab->a[m].other_off;
						ab->a[m].other_off = (uint32_t)ab->a[m].srt;
						if(trev) ab->a[m].other_off = tl - (ab->a[m].other_off+1-tspan) - 1;///looks like a bug
						ab->a[m].self_off = ab->mz.a[ab->a[m].self_off].pos;
						ab->a[m].srt = ab->a[m].self_off; ab->a[m].srt <<= 32; ab->a[m].srt |= ab->a[m].other_off;
						///debug
						/**
						recover_UC_Read_sub_region(dbg, ab->a[m].other_off + 1 - tspan, tspan, trev, rref, tid);
						if(recalu_minimizer0(dbg, tspan, !(asm_opt.flag & HA_F_NO_HPC), mz_k, thash, &tq, &rpos, &rspan) && (rpos + 1 == tspan) && (tspan == rspan)) {
							// fprintf(stderr, "-0-[M::%s]\n", __func__);
						} else {
							fprintf(stderr, "-1-[M::%s]\ttrev::%lu\n", __func__, trev);
						}
						///debug
						hpc_ext_check(rref, tid, ab->a[m].other_off + 1 - tspan, ab->a[m].other_off + 1, tl, trev, dbg);
						**/

						m++;
					}
					if(m > 1) radix_sort_ha_an1(ab->a, ab->a + m);
					for (i = 0, n0 = n; i < m; i++) {
						p = &cl->list[n++];
						p->readID = tid;
						p->strand = trev;
						p->offset = ab->a[i].other_off;
						p->self_offset = ab->a[i].self_off;
						if(((ab->a[i].cnt>>8) < max_cnt) && ((ab->a[i].cnt>>8) > min_cnt)){
							p->cnt = 1;
						} else if((ab->a[i].cnt>>8) <= min_cnt) {
							p->cnt = 2;
						} else{
							p->cnt = 1 + (((ab->a[i].cnt>>8) + (max_cnt<<1) - 1)/(max_cnt<<1));
							p->cnt = pow(p->cnt, 1.1);
						}
						if(p->cnt > ((uint32_t)(0xffffffu))) p->cnt = 0xffffffu;
						p->cnt <<= 8; p->cnt |= (((uint32_t)(0xffu))&(ab->a[i].cnt));
					}

					if(m > 0) {
						zn += lchain_qdp_mcopy_fast(cl, n0, n-n0, zn, &(cl->chainDP), olst, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_thres, 
																		rid, rl, tl, quick_check, apend_be, gen_off, enable_mcopy, mcopy_rate, mcopy_khit_cut, 1);
					}

					if(sv && olst->length > olst_n) {
						oidx->a[ol] >>= 32; oidx->a[ol] <<= 32; oidx->a[ol] |= ((uint64_t)((uint32_t)-1));
					}
				}
			}

            l = k;
        }
    }
    cl->length = zn;

	for (k = m = 0; k < oidx->n; k++) {
		if(((uint32_t)oidx->a[k]) == ((uint32_t)-1)) continue;
		oidx->a[m++] = oidx->a[k];
	}
	// fprintf(stderr, "[M::%s]\ttot::%lu\tremain::%lu\n", __func__, (uint64_t)oidx->n, m);
	oidx->n = m;

	for (i = 0; i < olst->length; ++i) olst->list[i].align_length = 0;

	// minimizers_qgen0(ab, rs, rl, mz_w, mz_k, cl, k_flag, ha_flt_tab, ha_idx, rref, dbg_ct, sp, high_occ, low_occ);
	
	// lchain_qgen_mcopy_fast(cl, overlap_list, rid, rl, rref, apend_be, max_n_chain, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_thres, quick_check, gen_off, enable_mcopy, mcopy_rate, chain_cutoff, mcopy_khit_cut, sp);
}


void h_ec_lchain_re_gen_srt(ha_abuf_t *ab, ha_pt_t *ha_idx, overlap_region_alloc *olst, Candidates_list *cl)
{
	uint64_t i, k; int j, n; ha_mz1_t *z; seed1_t *s; 
    
    clear_Candidates_list(cl); ab->n_a = 0;

	// minimizer of queried read
    if (ab->mz.m > ab->old_mz_m) {
        ab->old_mz_m = ab->mz.m;
        REALLOC(ab->seed, ab->old_mz_m);
    }

    for (i = 0, ab->n_a = 0; i < ab->mz.n; ++i) {
        ab->seed[i].a = ha_pt_get(ha_idx, ab->mz.a[i].x, &n);
        ab->seed[i].n = n;
        ab->n_a += n;
    }

	if (ab->n_a > ab->m_a) {
        ab->m_a = ab->n_a;
        REALLOC(ab->a, ab->m_a);
    }

	for (i = 0, k = 0; i < ab->mz.n; ++i) {
        ///z is one of the minimizer
        z = &ab->mz.a[i]; s = &ab->seed[i];

        for (j = 0; j < s->n; ++j) {
            const ha_idxpos_t *y = &s->a[j];
            anchor1_t *an = &ab->a[k++];
            uint8_t rev = z->rev == y->rev? 0 : 1;

			an->other_off = y->span;
			an->self_off = i;
            ///an->cnt: cnt<<8|span
            an->cnt = s->n; if(an->cnt > ((uint32_t)(0xffffffu))) an->cnt = 0xffffffu;
            an->cnt <<= 8; an->cnt |= ((z->span <= ((uint32_t)(0xffu)))?z->span:((uint32_t)(0xffu)));
			an->srt = (uint64_t)y->rid<<33 | (uint64_t)rev<<32 | y->pos;
        }
    }

	// copy over to _cl_
    if (ab->m_a >= (uint64_t)cl->size) {
        cl->size = ab->m_a;
        REALLOC(cl->list, cl->size);
    }
	clear_overlap_region_alloc(olst);
	radix_sort_ha_an1(ab->a, ab->a + ab->n_a);
}

uint64_t h_ec_lchain_re_gen_qry(ha_abuf_t *ab, uint64_t *k, uint64_t *l, uint64_t *i, uint64_t *idx_a, uint64_t idx_n, uint64_t *tid, uint64_t *trev)
{
	while ((*k) <= ab->n_a) {
        if ((*k) == ab->n_a || (ab->a[*k].srt>>32) != (ab->a[*l].srt>>32)) {
			for (; ((*i) < idx_n) && ((idx_a[*i]>>32) < (ab->a[*l].srt>>32)); (*i)++);
			if(((*i) < idx_n) && ((ab->a[*l].srt>>32) == (idx_a[*i]>>32))) {
				(*tid) = ab->a[*l].srt>>33; (*trev) = (ab->a[*l].srt>>32) & 1; 
				
				return 1;
			}
            (*l) = (*k);
        }
		++(*k);
    }
	return 0;
}

uint64_t h_ec_lchain_re_chn(ha_abuf_t *ab, uint64_t si, uint64_t ei, uint32_t rid, char* rs, uint64_t rl, uint64_t tid, char* ts, uint64_t tl, uint64_t trev, uint64_t mz_w, uint64_t mz_k, overlap_region_alloc *olst, Candidates_list *cl, double bw_thres, 
								 int apend_be, uint64_t max_cnt, uint64_t min_cnt, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, int64_t max_skip, int64_t max_iter, int64_t max_dis, int64_t quick_check, double chn_pen_gap, double chn_pen_skip, tiny_queue_t *tq, asg16_v *scc, int64_t *n, int64_t *zn)
{
	// char dbg[768]; uint64_t rpos, rspan, thash;

	int64_t ok, nk, ck, n0; uint64_t i, m, tspan, on0; k_mer_hit *p;
	for (i = si, m = ok = nk = ck = 0; i < ei; i++) {
		if(!recalu_minimizer_non_retrieve(rid, &(ab->a[i]), &(scc[tid]), &ok, &nk, &ck, mz_k, ab->mz.a[ab->a[i].self_off].x, tq, ts, tl, trev)) continue;
		///debug
		// thash = ab->mz.a[ab->a[i].self_off].x;
		
		ab->a[m] = ab->a[i]; tspan = ab->a[m].other_off;
		ab->a[m].other_off = (uint32_t)ab->a[m].srt;
		if(trev) ab->a[m].other_off = tl - (ab->a[m].other_off+1-tspan) - 1;///looks like a bug
		ab->a[m].self_off = ab->mz.a[ab->a[m].self_off].pos;
		ab->a[m].srt = ab->a[m].self_off; ab->a[m].srt <<= 32; ab->a[m].srt |= ab->a[m].other_off;


		///debug
		// recover_UC_Read_sub_region(dbg, ab->a[m].other_off + 1 - tspan, tspan, trev, &R_INF, tid);
		// if(recalu_minimizer0(dbg, tspan, !(asm_opt.flag & HA_F_NO_HPC), mz_k, thash, tq, &rpos, &rspan) && (rpos + 1 == tspan) && (tspan == rspan)) {
		//  	// fprintf(stderr, "-0-[M::%s]\n", __func__);
		// } else {
		//  	fprintf(stderr, "-1-[M::%s]\ttrev::%lu\n", __func__, trev);
		// }
		///debug
		// hpc_ext_check(&R_INF, tid, ab->a[m].other_off + 1 - tspan, ab->a[m].other_off + 1, tl, trev, dbg);

		m++;
	}

	if(m > 1) radix_sort_ha_an1(ab->a, ab->a + m);

	for (i = 0, n0 = (*n); i < m; i++) {
		p = &cl->list[(*n)++];
		p->readID = tid;
		p->strand = trev;
		p->offset = ab->a[i].other_off;
		p->self_offset = ab->a[i].self_off;
		if(((ab->a[i].cnt>>8) < max_cnt) && ((ab->a[i].cnt>>8) > min_cnt)){
			p->cnt = 1;
		} else if((ab->a[i].cnt>>8) <= min_cnt) {
			p->cnt = 2;
		} else{
			p->cnt = 1 + (((ab->a[i].cnt>>8) + (max_cnt<<1) - 1)/(max_cnt<<1));
			p->cnt = pow(p->cnt, 1.1);
		}
		if(p->cnt > ((uint32_t)(0xffffffu))) p->cnt = 0xffffffu;
		p->cnt <<= 8; p->cnt |= (((uint32_t)(0xffu))&(ab->a[i].cnt));
	}

	// if(rid == 3196 && tid == 3199) fprintf(stderr, "[M::%s]\ttid::%lu\ttrev::%lu\told::%lu\tnew::%lu\n", __func__, tid, trev, k - l, m);
	// fprintf(stderr, "[M::%s]\ttid::%lu\ttrev::%lu\told::%lu\tnew::%lu\n", __func__, tid, trev, ei - si, m);

	if(m > 0 && tid != rid) {
		on0 = olst->length;
		*zn += lchain_qdp_mcopy_fast(cl, n0, (*n)-n0, *zn, &(cl->chainDP), olst, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_thres, 
																	rid, rl, tl, quick_check, apend_be, gen_off, enable_mcopy, mcopy_rate, mcopy_khit_cut, 1);
		cl->length = *zn;
		if(olst->length > on0) return 1;
	}

	cl->length = *zn;
	return 0;
}

uint64_t get_mz1(const char *str, int len, int w, int k, uint32_t rid, int is_hpc, ha_abuf_t *ab, const void *hf, ha_pt_t *ha_idx, int sample_dist, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, ha_pt_t *pt, int min_freq, int32_t dp_min_len, float dp_e, st_mt_t *mt, int32_t ws, int32_t is_unique, void *km, uint64_t beg_i)
{
	ab->mz.n = beg_i;
	// get the list of anchors
	mz1_ha_sketch(str, len, w, k, 0, !(asm_opt.flag & HA_F_NO_HPC), &ab->mz, hf, asm_opt.mz_sample_dist, k_flag, dbg_ct, NULL, -1, asm_opt.dp_min_len, -1, mt, asm_opt.mz_rewin, 0, NULL);

	radix_sort_ha_mz1_v_srt(ab->mz.a + beg_i, ab->mz.a + ab->mz.n);

	if(ha_idx) {
		uint64_t i;

		if (ab->mz.m > ab->old_mz_m) {
			ab->old_mz_m = ab->mz.m;
			REALLOC(ab->seed, ab->old_mz_m);
		}

		for (i = 0; i < ab->mz.n; ++i) {
			ab->seed[i].a = NULL;
			ab->seed[i].n = ha_pt_cnt(ha_idx, ab->mz.a[i].x);
		}
	}

	return ab->mz.n;
}

void get_pi_ec_chain(ha_abuf_t *ab, uint64_t rid, uint64_t rl, uint32_t tid, char* ts, uint64_t tl, uint64_t mz_w, uint64_t mz_k, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, 
								int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, /**uint32_t is_accurate,**/ uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, 
								int64_t max_skip, int64_t max_iter, int64_t max_dis, int64_t quick_check, double chn_pen_gap, double chn_pen_skip)
{
	extern void *ha_flt_tab;
	extern ha_pt_t *ha_idx;
	uint64_t rn = ab->mz.n, tn, cn = 0;

	tn = get_mz1(ts, tl, mz_w, mz_k, 0, !(asm_opt.flag & HA_F_NO_HPC), ab, ha_flt_tab, NULL, asm_opt.mz_sample_dist, k_flag, dbg_ct, NULL, -1, asm_opt.dp_min_len, -1, sp, asm_opt.mz_rewin, 0, NULL, rn);

	cn = lchain_qgen_mcopy_fast_re0(ab, ha_idx, ab->mz.a, rn, ab->mz.a + rn, tn - rn, tid, cl, high_occ, low_occ);

	if(cn) {
		lchain_qgen_mcopy_fast_re1(cl, cl->length - cn, overlap_list, rid, rl, tl, apend_be, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_thres, quick_check, gen_off, enable_mcopy, mcopy_rate, mcopy_khit_cut);
	}

	ab->mz.n = rn;
}


int64_t ug_map_lchain(ha_abufl_t *ab, uint32_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k, const ul_idx_t *uref, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, double bw_thres_sec,
								 int max_n_chain, int apend_be, kvec_t_u8_warp* k_flag, overlap_region* f_cigar, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t is_accurate, 
								 uint32_t gen_off, double mcopy_rate, uint32_t mcopy_khit_cut, uint32_t is_hpc, ha_mzl_t *res, uint64_t res_n, ha_mzl_t *idx, uint64_t idx_n, uint64_t mzl_cutoff, uint64_t chain_cutoff, kv_u_trans_t *kov)
{
	int64_t max_skip, max_iter, max_dis, quick_check; double chn_pen_gap, chn_pen_skip;
	set_lchain_dp_op(is_accurate, mz_k, &max_skip, &max_iter, &max_dis, &chn_pen_gap, &chn_pen_skip, &quick_check);
	if((!overlap_list) || (!cl)) {
		ab->mz.n = 0;
		mz2_ha_sketch(rs, rl, mz_w, mz_k, rid, is_hpc, &ab->mz, NULL, asm_opt.mz_sample_dist, k_flag, dbg_ct, NULL, -1, asm_opt.dp_min_len, -1, sp, asm_opt.mz_rewin, 0, NULL);
		if(res) memcpy(res, ab->mz.a, ab->mz.n * (sizeof((*(ab->mz.a)))));
		return ab->mz.n;
	} else {
		sp->n = 0;
		minimizers_qgen_input(ab, rid, rs, rl, mz_w, mz_k, cl, k_flag, ha_flt_tab, ha_idx, NULL, uref, dbg_ct, sp, high_occ, low_occ, res, res_n, idx, idx_n, mzl_cutoff, chain_cutoff, kov);
		lchain_qgen_mcopy_input(cl, overlap_list, rid, rl, NULL, uref, apend_be, max_n_chain, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_thres, bw_thres_sec, quick_check, gen_off, mcopy_rate, mcopy_khit_cut, sp);
		return 0;
	}
}

int64_t ug_map_lchain_simple(ha_abufl_t *ab, uint32_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k, const ul_idx_t *uref, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres,
								 int max_n_chain, int apend_be, kvec_t_u8_warp* k_flag, overlap_region* f_cigar, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t is_accurate, 
								 uint32_t gen_off, double mcopy_rate, uint32_t mcopy_khit_cut, uint32_t is_hpc, ha_mzl_t *res, uint64_t res_n, ha_mzl_t *idx, uint64_t idx_n, uint64_t mzl_cutoff, uint64_t chain_cutoff)
{
	int64_t max_skip, max_iter, max_dis, quick_check; double chn_pen_gap, chn_pen_skip;
	set_lchain_dp_op(is_accurate, mz_k, &max_skip, &max_iter, &max_dis, &chn_pen_gap, &chn_pen_skip, &quick_check);
	if((!overlap_list) || (!cl)) {
		ab->mz.n = 0;
		mz2_ha_sketch(rs, rl, mz_w, mz_k, rid, is_hpc, &ab->mz, NULL, asm_opt.mz_sample_dist, k_flag, dbg_ct, NULL, -1, asm_opt.dp_min_len, -1, sp, asm_opt.mz_rewin, 0, NULL);
		if(res) memcpy(res, ab->mz.a, ab->mz.n * (sizeof((*(ab->mz.a)))));
		return ab->mz.n;
	} else {
		sp->n = 0;
		minimizers_qgen_input(ab, rid, rs, rl, mz_w, mz_k, cl, k_flag, ha_flt_tab, ha_idx, NULL, uref, dbg_ct, sp, high_occ, low_occ, res, res_n, idx, idx_n, mzl_cutoff, chain_cutoff, NULL);
		// lchain_qgen_mcopy_input(cl, overlap_list, rid, rl, NULL, uref, apend_be, max_n_chain, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_thres, bw_thres_sec, quick_check, gen_off, mcopy_rate, mcopy_khit_cut, sp);
		lchain_qgen_mcopy(cl, overlap_list, rid, rl, NULL, uref, apend_be, max_n_chain, max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, bw_thres, quick_check, gen_off, mcopy_rate, chain_cutoff, mcopy_khit_cut, sp);
		return 0;
	}
}