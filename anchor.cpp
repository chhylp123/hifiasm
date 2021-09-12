#include <stdio.h>
#include <math.h>
#include "htab.h"
#include "ksort.h"
#include "Hash_Table.h"

#define HA_KMER_GOOD_RATIO 0.333

typedef struct { // this struct is not strictly necessary; we can use k_mer_pos instead, with modifications
	uint64_t srt;
	uint32_t self_off;
	uint32_t other_off;
	uint32_t cnt;
} anchor1_t;

#define an_key1(a) ((a).srt)
#define an_key2(a) ((a).self_off)
KRADIX_SORT_INIT(ha_an1, anchor1_t, an_key1, 8)
KRADIX_SORT_INIT(ha_an2, anchor1_t, an_key2, 4)

#define oreg_xs_lt(a, b) (((uint64_t)(a).x_pos_s<<32|(a).x_pos_e) < ((uint64_t)(b).x_pos_s<<32|(b).x_pos_e))
KSORT_INIT(or_xs, overlap_region, oreg_xs_lt)

#define oreg_ss_lt(a, b) ((a).shared_seed > (b).shared_seed) // in the decending order
KSORT_INIT(or_ss, overlap_region, oreg_ss_lt)

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
HType *sf##_init(void){return (HType*)calloc(1, sizeof(HType));}\
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
	mz1_ha_sketch(ucr->seq, ucr->length, asm_opt.mz_win, asm_opt.k_mer_length, 0, !(asm_opt.flag & HA_F_NO_HPC), &ab->mz, ha_flt_tab, asm_opt.mz_sample_dist, k_flag, dbg_ct, NULL, -1, asm_opt.dp_min_len, -1, sp, asm_opt.mz_rewin, 0);
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

	calculate_overlap_region_by_chaining(cl, overlap_list, chain_idx, rid, ucr->length, &R_INF, bw_thres, keep_whole_chain, f_cigar);

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
				50, ua->a[(*f_cigar).x_id].len, ua->a[(*f_cigar).y_id].len);
		
		
		// if ((*f_cigar).x_id != (*f_cigar).y_id)
		if ((*f_cigar).x_id != (*f_cigar).y_id && chain_len > mz_occ*mz_rate)
        {
			append_utg_inexact_overlap_region_alloc(overlap_list, f_cigar, ua, add_beg_end);
        }
    }
}

/**
void ha_get_inter_candidates(ha_abufl_t *ab, int64_t uid, ma_utg_v *ua, overlap_region_alloc *ovlp, 
Candidates_list *cl, double bw_thres, int max_n_chain, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* chain_idx, 
void *ha_flt_tab, ha_pt_t *ha_idx, overlap_region* f_cigar, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, 
double chain_match_rate, uint32_t uk, uint32_t uw, uint32_t is_hpc, uint32_t h_occ)
{

	uint32_t i, rlen;
    uint64_t k, l;
	ma_utg_t *u = &(ua->a[uid]);
	// prepare
    clear_Candidates_list(cl);
    clear_overlap_region_alloc(ovlp);
	ab->mz.n = 0, ab->n_a = 0;

	mz2_ha_sketch(u->s, u->len, uw, uk, uid, is_hpc, &ab->mz, ha_flt_tab, asm_opt.mz_sample_dist, k_flag, dbg_ct,	
	NULL, -1, asm_opt.dp_min_len, -1, sp, asm_opt.mz_rewin, 1);
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
            an->self_off = rev? u->len - 1 - (z->pos + 1 - z->span) : z->pos;
            an->cnt = s->cnt;
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
    // calculate_overlap_region_by_chaining(cl, ovlp, chain_idx, rid, ucr->length, &R_INF, bw_thres, keep_whole_chain, f_cigar);
	// calculate_ug_chaining(cl, overlap_list, chain_idx, rid, ua, bw_thres, keep_whole_chain, f_cigar, ab->mz.n, chain_match_rate);


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
    
    if ((int)ovlp->length > max_n_chain) {
        int32_t w, n[4], s[4];
        n[0] = n[1] = n[2] = n[3] = 0, s[0] = s[1] = s[2] = s[3] = 0;
        ks_introsort_or_ss(ovlp->length, ovlp->list);
        for (i = 0; i < (uint32_t)ovlp->length; ++i) {
            const overlap_region *r = &ovlp->list[i];
            w = ha_ov_type(r, rlen);
            ++n[w];
            if ((int)n[w] == max_n_chain) s[w] = r->shared_seed;
        }
        if (s[0] > 0 || s[1] > 0 || s[2] > 0 || s[3] > 0) {
            // n[0] = n[1] = n[2] = n[3] = 0;
            for (i = 0, k = 0; i < (uint32_t)ovlp->length; ++i) {
                overlap_region *r = &ovlp->list[i];
                w = ha_ov_type(r, rlen);
                // ++n[w];
                // if (((int)n[w] <= max_n_chain) || (r->shared_seed >= s[w] && s[w] >= (asm_opt.k_mer_length<<1))) {
                if (r->shared_seed >= s[w]) {
                    if ((uint32_t)k != i) {
                        overlap_region t;
                        t = ovlp->list[k];
                        ovlp->list[k] = ovlp->list[i];
                        ovlp->list[i] = t;
                    }
                    ++k;
                }
            }
            ovlp->length = k;
        }
    }
    ///ks_introsort_or_xs(overlap_list->length, overlap_list->list);
}
**/


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

void ha_sort_list_by_anchor(overlap_region_alloc *overlap_list)
{
	ks_introsort_or_xs(overlap_list->length, overlap_list->list);
}
