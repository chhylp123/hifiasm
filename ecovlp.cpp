#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "Correct.h"
#include "Process_Read.h"
#include "ecovlp.h"
#include "kthread.h"
#include "htab.h"
#define HA_KMER_GOOD_RATIO 0.333
#define E_KHIT 31

#define generic_key(x) (x)
KRADIX_SORT_INIT(ec16, uint16_t, generic_key, 2)
KRADIX_SORT_INIT(ec32, uint32_t, generic_key, 4)
KRADIX_SORT_INIT(ec64, uint64_t, generic_key, 8)

void h_ec_lchain(ha_abuf_t *ab, uint32_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k, All_reads *rref, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, 
								 int max_n_chain, int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t is_accurate, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t chain_cutoff, uint32_t mcopy_khit_cut);



ec_ovec_buf_t* gen_ec_ovec_buf_t(uint32_t n, uint32_t is_final, uint32_t save_ov)
{
    uint32_t k; ec_ovec_buf_t0 *z = NULL;
    ec_ovec_buf_t *p = NULL; CALLOC(p, 1);
    p->n = n; CALLOC(p->a, p->n);
    for (k = 0; k < p->n; k++) {
        z = &(p->a[k]);
        z->is_final = !!is_final; z->save_ov = !!save_ov;
        init_UC_Read(&z->self_read);
	    init_UC_Read(&z->ovlp_read);
	    init_Candidates_list(&z->clist);
	    init_overlap_region_alloc(&z->olist);

        init_fake_cigar(&(z->tmp.f_cigar));
        memset(&(z->tmp.w_list), 0, sizeof(z->tmp.w_list));
        CALLOC(z->tmp.w_list.a, 1); z->tmp.w_list.n = z->tmp.w_list.m = 1;

        // kv_init(z->b_buf.a);
        kv_init(z->r_buf.a);
        kv_init(z->k_flag.a);
        kv_init(z->sp);
        kv_init(z->pidx);
	    kv_init(z->v64);
        kv_init(z->v32);
        kv_init(z->v16);
        init_bit_extz_t(&(z->exz), 31);

        z->ab = ha_abuf_init();
        if (!z->is_final) {
            init_Cigar_record(&z->cigar);
            // init_Graph(&b->POA_Graph);
            // init_Graph(&b->DAGCon);
            init_Correct_dumy(&z->correct);
            InitHaplotypeEvdience(&z->hap);
        }
    }
    
    return p;
}

void destroy_cns_gfa(cns_gfa *p)
{
    size_t k;
    for (k = 0; k < p->m; k++) {
        kv_destroy(p->a[k].in);
        kv_destroy(p->a[k].ou);
    }
    free(p->a);
}

void destroy_ec_ovec_buf_t(ec_ovec_buf_t *p)
{
    uint32_t k; ec_ovec_buf_t0 *z = NULL;
    for (k = 0; k < p->n; k++) {
        z = &(p->a[k]);
        destory_UC_Read(&z->self_read);
        destory_UC_Read(&z->ovlp_read);
        destory_Candidates_list(&z->clist);
	    destory_overlap_region_alloc(&z->olist);

        destory_fake_cigar(&(z->tmp.f_cigar));
        free(z->tmp.w_list.a); free(z->tmp.w_list.c.a);

        kv_destroy(z->r_buf.a);
        kv_destroy(z->k_flag.a);
        kv_destroy(z->sp);
        kv_destroy(z->pidx);
	    kv_destroy(z->v64);
        kv_destroy(z->v32);
        kv_destroy(z->v16);
        destroy_bit_extz_t(&(z->exz));

        ha_abuf_destroy(z->ab);
        if (!z->is_final) {
            destory_Cigar_record(&z->cigar);
            // destory_Graph(&b->POA_Graph);
            // destory_Graph(&b->DAGCon);
            destory_Correct_dumy(&z->correct);
            destoryHaplotypeEvdience(&z->hap);
        }
        destroy_cns_gfa(&(z->cns));

        asm_opt.num_bases += z->num_read_base;
		asm_opt.num_corrected_bases += z->num_correct_base;
		asm_opt.num_recorrected_bases += z->num_recorrect_base;
		// asm_opt.mem_buf += ha_ovec_mem(b[i], NULL);
    }
    free(p->a); free(p);

    fprintf(stderr, "[M::%s-chains] #->%lld\n", __func__, asm_opt.num_bases);
    fprintf(stderr, "[M::%s-passed-chains-0] #->%lld\n", __func__, asm_opt.num_corrected_bases);
    fprintf(stderr, "[M::%s-cis-chains-1] #->%lld\n", __func__, asm_opt.num_recorrected_bases);
}

void prt_chain(overlap_region_alloc *o)
{
    uint64_t k;
    for (k = 0; k < o->length; k++) {
        fprintf(stderr, "[M::%s]\t#%u\tlen::%lu\t%u\t%u\t%c\t#%u\tlen::%lu\t%u\t%u\tsc::%d\taln::%u\terr::%u\n", __func__, o->list[k].x_id, Get_READ_LENGTH(R_INF, o->list[k].x_id), o->list[k].x_pos_s, o->list[k].x_pos_e+1, "+-"[o->list[k].y_pos_strand], 
                                                                                        o->list[k].y_id, Get_READ_LENGTH(R_INF, o->list[k].y_id), o->list[k].y_pos_s, o->list[k].y_pos_e+1, o->list[k].shared_seed, o->list[k].align_length, o->list[k].non_homopolymer_errors);
    }    
}

overlap_region *fetch_aux_ovlp(overlap_region_alloc* ol) /// exactly same to gen_aux_ovlp
{
	if (ol->length + 1 > ol->size) {
		uint64_t sl = ol->size;
        ol->size = ol->length + 1;
        kroundup64(ol->size);
        REALLOC(ol->list, ol->size);
        /// need to set new space to be 0
        memset(ol->list + sl, 0, sizeof(overlap_region)*(ol->size - sl));
	}
	return &(ol->list[ol->length+1]);
}

typedef struct {
	ul_ov_t *c_idx;
    asg64_v *idx;
    int64_t i, i0, srt_n, rr;
    uint64_t mms, mme;
} cc_idx_t;


///[s, e)
int64_t extract_sub_cigar_mm(overlap_region *z, int64_t s, int64_t e, ul_ov_t *p, uint64_t *ct)
{
    int64_t wk = ovlp_cur_wid(*p), xk = ovlp_cur_xoff(*p), yk = ovlp_cur_yoff(*p), ck = ovlp_cur_coff(*p), os, oe, t;
    bit_extz_t ez; int64_t bd = ovlp_bd(*p), s0, e0; 
    s0 = ((int64_t)(z->w_list.a[wk].x_start)) + bd; 
    e0 = ((int64_t)(z->w_list.a[wk].x_end)) + 1 - bd;
    if(s < s0) s = s0; if(e > e0) e = e0;///exclude boundary
    if(s >= e) return -1;
    os = MAX(s, s0); oe = MIN(e, e0);
    if(oe <= os) return -1;

    set_bit_extz_t(ez, (*z), wk);
    if(!ez.cigar.n) return -1;
    int64_t cn = ez.cigar.n, op; int64_t ws, we, ovlp; 
    if((ck < 0) || (ck > cn)) {//(*ck) == cn is allowed
        ck = 0; xk = ez.ts; yk = ez.ps;
    }
    
    while (ck > 0 && xk >= s) {///x -> t; y -> p; first insertion and then match/mismatch
        --ck;
        op = ez.cigar.a[ck]>>14;
        if(op!=2) xk -= (ez.cigar.a[ck]&(0x3fff));
        if(op!=3) yk -= (ez.cigar.a[ck]&(0x3fff)); 
    }

     //some cigar will span s or e
    while (ck < cn && xk < e) {//[s, e)
        ws = xk; 
        op = ez.cigar.a[ck]>>14;
        ///op == 3: -> x; op == 2: -> y;
        if(op!=2) xk += (ez.cigar.a[ck]&(0x3fff));
        if(op!=3) yk += (ez.cigar.a[ck]&(0x3fff)); 
        ck++; we = xk;

        os = MAX(s, ws); oe = MIN(e, we);
        ovlp = ((oe>os)? (oe-os):0);
        if(op != 2) {
            if(!ovlp) continue;
        } else {///ws == we
            if(ws < s || ws >= e) continue;
        }
        

        if(op == 0) {
            for (t = os + 1; t < oe; t++) {
                ct[(t-s)<<1]++; ct[(t-s)<<1] += ((uint64_t)(0x100000000));
                ct[((t-s)<<1)+1]++; ct[((t-s)<<1)+1] += ((uint64_t)(0x100000000));
            }

            t = os;
            if(t < oe) {
                ct[(t-s)<<1]++; ct[(t-s)<<1] += ((uint64_t)(0x100000000));
                if(os > ws) {
                    ct[((t-s)<<1)+1]++; ct[((t-s)<<1)+1] += ((uint64_t)(0x100000000));
                }
            }
        } else if(op!=2) {
            for (t = os + 1; t < oe; t++) {
                ct[(t-s)<<1]++; 
                ct[((t-s)<<1)+1]++;
            }

            t = os;
            if(t < oe) {
                ct[(t-s)<<1]++; 
                if(os > ws) {
                    ct[((t-s)<<1)+1]++;
                }
            }
        } else {
            ct[((ws-s)<<1)+1]++; ///ct[((ws-s)<<1)+1] += ((uint64_t)(0x100000000));
        }
    }

    ovlp_cur_xoff(*p) = xk; ovlp_cur_yoff(*p) = yk; ovlp_cur_coff(*p) = ck; ovlp_cur_ylen(*p) = 0;
    
    return 1;
}

#define simp_vote_len 6

///[s, e)
uint32_t extract_sub_cigar_ii(overlap_region *z, All_reads *rref, int64_t s, int64_t e, UC_Read* tu, ul_ov_t *p)
{
    int64_t wk = ovlp_cur_wid(*p), xk = ovlp_cur_xoff(*p), yk = ovlp_cur_yoff(*p), ck = ovlp_cur_coff(*p), os, oe, ol;
    bit_extz_t ez; int64_t bd = ovlp_bd(*p), s0, e0, ii[2], it[2]; uint32_t res = (uint32_t)-1;
    s0 = ((int64_t)(z->w_list.a[wk].x_start)) + bd; 
    e0 = ((int64_t)(z->w_list.a[wk].x_end)) + 1 - bd;
    if(s < s0) s = s0; if(e > e0) e = e0;///exclude boundary
    if(s > e) return -1;///it is possible s == e
    os = MAX(s, s0); oe = MIN(e, e0);
    if(oe < os) return -1;///it is possible os == oe

    set_bit_extz_t(ez, (*z), wk);
    if(!ez.cigar.n) return -1;
    int64_t cn = ez.cigar.n; uint16_t op; int64_t ws, we, wts, wte, ovlp, cc = 0, cci; 
    if((ck < 0) || (ck > cn)) {//(*ck) == cn is allowed
        ck = 0; xk = ez.ts; yk = ez.ps;
    }
    
    while (ck > 0 && xk >= s) {///x -> t; y -> p; first insertion and then match/mismatch
        --ck;
        op = ez.cigar.a[ck]>>14;
        if(op!=2) xk -= (ez.cigar.a[ck]&(0x3fff));
        if(op!=3) yk -= (ez.cigar.a[ck]&(0x3fff)); 
    }

    char cm[4]; cm[0] = 'M'; cm[1] = 'S'; cm[2] = 'I'; cm[3] = 'D';  
    //some cigar will span s or e
    ii[0] = ii[1] = it[0] = it[1] = -1; res = cc = 0;
    while (ck < cn && xk < e) {//[s, e)
        ws = xk; wts = yk;
        op = ez.cigar.a[ck]>>14; ol = (ez.cigar.a[ck]&(0x3fff));
        ///op == 3: -> x; op == 2: -> y;
        if(op!=2) xk += ol;
        if(op!=3) yk += ol; 
        ck++; we = xk; wte = yk;

        // if(s == 10480) {
        //     fprintf(stderr, "[%ld, %ld)\t%c\n", ws, we, cm[op]);
        // }

        os = MAX(s, ws); oe = MIN(e, we);
        ovlp = ((oe>os)? (oe-os):0);
        if(op != 2) {
            if(!ovlp) continue;
        } else {///ws == we
            if(ws < s || ws >= e) continue;
        }
        



        if(ii[0] == -1) {
            ii[0] = os;
            if(op < 2) {
                it[0] = os - ws + wts;
            } else {///op == 2: more y; p == 3: more x
                it[0] = wts;
            } 
        }

        ii[1] = oe;
        if(op < 2) {
            it[1] = oe - ws + wts;
        } else {///op == 2: more y; p == 3: more x
            it[1] = wte;
        } 




        if(op != 2) ol = oe-os;
        cc += ol;
        fprintf(stderr, "%ld%c", ol, cm[op]);
        if(cc <= simp_vote_len) {
            for (cci = 0; cci < ol; cci++) {
                res <<= 2; res |= op;
            }
        }
    }

    while (ck < cn && xk <= e) {//[s, e)
        ws = xk; wts = yk;
        op = ez.cigar.a[ck]>>14; ol = (ez.cigar.a[ck]&(0x3fff));
        if(op != 2) break;
        yk += (ez.cigar.a[ck]&(0x3fff)); 
        ck++; we = xk; wte = yk;
        if(ws >= s && ws <= e) {
            
            if(ii[0] == -1) {
                ii[0] = os; it[0] = wts;
            }
            ii[1] = oe; it[1] = wte;
            
            cc += ol;
            fprintf(stderr, "%ld%c", ol, cm[op]);
            if(cc <= simp_vote_len) {
                for (cci = 0; cci < ol; cci++) {
                    res <<= 2; res |= op;
                }
            }
        }
    }

    fprintf(stderr, "\tx::[%ld, %ld)\ty::[%ld, %ld)\tcc::%ld\n", ii[0], ii[1], it[0], it[1], cc);
    if((cc <= simp_vote_len) 
            && (ii[1] >= ii[0]) && (ii[1] - ii[0] <= simp_vote_len) 
            && (it[1] >= it[0]) && (it[1] - it[0] <= simp_vote_len)) {
        // ii[0] = ii[0] - s; ii[1] = e - ii[1];
        if((ii[0] == s) && (ii[1] == e)) {
            op = cc; op <<= 12; res |= op;

            char *ystr = NULL; res <<= 16; cc = it[1] - it[0]; op = 0;
            if(cc > 0) {
                UC_Read_resize(*tu, (it[1] - it[0])); ystr = tu->seq;
                recover_UC_Read_sub_region(ystr, it[0], (it[1] - it[0]), z->y_pos_strand, rref, z->y_id);
                
                for (cci = 0; cci < cc; cci++) {
                    op <<= 2; op |= seq_nt6_table[(uint32_t)(ystr[cci])];
                }
            } 
            res |= op;

            op = it[1] - it[0]; op <<= 12; res |= op;
        } else {
            res = (uint32_t)-1;
        }
    } else {
        res = (uint32_t)-1;
    }
    

    ovlp_cur_xoff(*p) = xk; ovlp_cur_yoff(*p) = yk; ovlp_cur_coff(*p) = ck; ovlp_cur_ylen(*p) = 0;
    
    return res;
}


uint64_t iter_cc_idx_t(overlap_region* ol, cc_idx_t *z, int64_t s, int64_t e, uint64_t is_reduce, uint64_t is_insert, uint64_t **ra)
{
    int64_t rm_n, q[2], os, oe; ul_ov_t *cp; uint64_t m; *ra = NULL;

    // if(s == 15816 && e == 15819) {
    //     fprintf(stderr, "[M::%s] is_reduce::%lu\n", __func__, is_reduce);
    // }

    if(is_reduce) {
        for (m = rm_n = z->srt_n; m < z->idx->n; m++) {
            cp = &(z->c_idx[z->idx->a[m]]);
            // if(s == 15816 && e == 15819) {
            //     fprintf(stderr, "-0-[M::%s] ii::%lu, ii0::%ld\n", __func__, z->idx->a[m], z->i0);
            // }
            
            q[0] = ol[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_start+ovlp_bd(*cp);
            q[1] = ol[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_end+1-ovlp_bd(*cp);
            os = MAX(q[0], s); oe = MIN(q[1], e);
            if((oe > os) || ((is_insert) && (s == e) && (s >= q[0]) && (s <= q[1]))) {
                z->idx->a[rm_n++] = z->idx->a[m];
            }
        }
        z->idx->n = rm_n;
    }

    for (; z->i < z->srt_n; ++z->i) {
        cp = &(z->c_idx[(uint32_t)z->idx->a[z->i]]);
        // if(s == 15816 && e == 15819) {
        //     fprintf(stderr, "-1-[M::%s] ii::%u, ii0::%ld\n", __func__, (uint32_t)z->idx->a[z->i], z->i0);
        // }
        q[0] = ol[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_start+ovlp_bd(*cp);
        q[1] = ol[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_end+1-ovlp_bd(*cp);
        if(q[0] > e) break;
        if((!is_insert) && (q[0] >= e)) break;
        os = MAX(q[0], s); oe = MIN(q[1], e);
        if((oe > os) || ((is_insert) && (s == e) && (s >= q[0]) && (s <= q[1]))) {
            kv_push(uint64_t, *(z->idx), ((uint32_t)z->idx->a[z->i]));
        }
    }

    (*ra) = z->idx->a + z->srt_n; 
    return z->idx->n - z->srt_n;
}

void debug_inter0(overlap_region* ol, ul_ov_t *c_idx, uint64_t *idx, int64_t idx_n, uint64_t *res, int64_t res_n, int64_t s, int64_t e, uint64_t is_insert, uint64_t is_hard_check, const char *cmd)
{
    ul_ov_t *cp; int64_t q[2], a_n = 0, i, k = 0, os, oe; 
    for (i = 0; i < idx_n; i++) {
        cp = &(c_idx[(uint32_t)idx[i]]);
        q[0] = ol[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_start+ovlp_bd(*cp);
        q[1] = ol[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_end+1-ovlp_bd(*cp);

        // fprintf(stderr, "%s[M::%s] tid::%u\t%.*s\twid::%u\tq::[%u, %u)\terr::%d\toerr::%u\n", cmd, __func__, ol[ovlp_id(*cp)].y_id, (int)Get_NAME_LENGTH(R_INF, ol[ovlp_id(*cp)].y_id), Get_NAME(R_INF, ol[ovlp_id(*cp)].y_id),
        // ovlp_cur_wid(*cp), ol[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_start, ol[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_end+1, ol[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].error, ol[ovlp_id(*cp)].non_homopolymer_errors);

        os = MAX(q[0], s); oe = MIN(q[1], e);
        if((oe > os) || ((is_insert) && (s == e) && (s >= q[0]) && (s <= q[1]))) {
            a_n++; 
            // if(!(((uint32_t)idx[i]) == res[k])) {
            //     fprintf(stderr, "[M::%s] a_n::%ld\tres_n::%ld\ts::%ld\te::%ld\ti::%ld\tk::%ld\n", __func__, a_n, res_n, s, e, i, k);
            // }
            if(is_hard_check) {
                assert(((uint32_t)idx[i]) == res[k++]);
            } else {
                for (; (k < res_n) && (((uint32_t)idx[i]) != res[k]); k++);
                assert(k < res_n);
            }
        }
    }
    // if(a_n != res_n) {
    //     fprintf(stderr, "[M::%s] a_n::%ld\tres_n::%ld\ts::%ld\te::%ld\tidx_n::%ld\n", __func__, a_n, res_n, s, e, idx_n);
    // }
    if(is_hard_check) {
        assert(a_n == res_n);
    } else {
        assert(a_n <= res_n);
    }
}

void prt_cigar0(uint64_t in, int64_t len)
{
    int64_t k; uint64_t mp;
    char cm[4]; cm[0] = 'M'; cm[1] = 'S'; cm[2] = 'I'; cm[3] = 'D';  
    for (k = 0; k < len; k++) {
        mp = len - 1 - k; mp <<= 1;
        fprintf(stderr, "%c", cm[(in >> mp)&3]);
    }
    fprintf(stderr, "\n");
}

void prt_bp0(uint64_t in, int64_t len)
{
    int64_t k; uint64_t mp;
    char cm[4]; cm[0] = 'A'; cm[1] = 'C'; cm[2] = 'G'; cm[3] = 'T';  
    for (k = 0; k < len; k++) {
        mp = len - 1 - k; mp <<= 1;
        fprintf(stderr, "%c", cm[(in >> mp)&3]);
    }
    fprintf(stderr, "\n");
}

uint64_t cns_gen0(overlap_region* ol, All_reads *rref, uint64_t s, uint64_t e, UC_Read* tu, cc_idx_t *idx, uint64_t occ_tot, double occ_max, asg32_v* b32, uint32_t *rc)
{
    if(e > s + simp_vote_len) return 0;///too long

    uint64_t *id_a = NULL, id_n, an = 0, oc[2]; b32->n = 0; uint32_t m, *a = NULL;
    id_n = iter_cc_idx_t(ol, idx, s, e, idx->rr, ((s==e)?1:0), &id_a);
    // debug_inter0(ol, idx->c_idx, idx->idx->a + idx->i0, idx->srt_n - idx->i0, id_a, id_n, s, e, ((s==e)?1:0), 0, "-1-");
    uint64_t k, l, q[2], os, oe; ul_ov_t *p; overlap_region *z; idx->rr = 0;
    fprintf(stderr, "[M::%s] [%lu, %lu) id_n::%lu\n", __func__, s, e, id_n);
    for (k = 0; k < id_n; k++) {
        p = &(idx->c_idx[id_a[k]]); z = &(ol[ovlp_id(*p)]); 
        q[0] = z->w_list.a[ovlp_cur_wid(*p)].x_start+ovlp_bd(*p);
        q[1] = z->w_list.a[ovlp_cur_wid(*p)].x_end+1-ovlp_bd(*p);

        // fprintf(stderr, "[M::%s] tid::%u\t%.*s\twid::%u\tq::[%u, %u)\terr::%d\toerr::%u\n", __func__, ol[ovlp_id(*p)].y_id, (int)Get_NAME_LENGTH(R_INF, ol[ovlp_id(*p)].y_id), Get_NAME(R_INF, ol[ovlp_id(*p)].y_id),
        // ovlp_cur_wid(*p), ol[ovlp_id(*p)].w_list.a[ovlp_cur_wid(*p)].x_start, ol[ovlp_id(*p)].w_list.a[ovlp_cur_wid(*p)].x_end+1, ol[ovlp_id(*p)].w_list.a[ovlp_cur_wid(*p)].error, ol[ovlp_id(*p)].non_homopolymer_errors);

        if(q[1] <= e) idx->rr = 1;
        os = MAX(q[0], s); oe = MIN(q[1], e);
        if((oe > os) || ((s == e) && (s >= q[0]) && (s <= q[1]))) {
        // if(oe >= os) {
            ///[-4-][-12-][-4-][-12-]
            ///[cigar_len][cigar][base_len][base]
            m = extract_sub_cigar_ii(z, rref, os, oe, tu, p); an++;
            if(m != ((uint32_t)-1)) {///no gap in both sides
                kv_push(uint32_t, *b32, m);
            }
        }
    }

    oc[0] = b32->n; oc[1] = an + 1; //+1 for the reference read
    fprintf(stderr, "-0-[M::%s] oc[0]::%lu, oc[1]::%lu\n", __func__, oc[0], oc[1]);
    if(((oc[0] > (oc[1]*occ_max)) && (oc[0] > (oc[1]-oc[0])) && (oc[1] >= occ_tot) && (oc[0] > 1))) {
        radix_sort_ec32(b32->a, b32->a+b32->n); an = 0;
        for (k = 1, l = 0; k <= b32->n; ++k) {
            if (k == b32->n || b32->a[k] != b32->a[l]) {
                if(k - l > an) {
                    an = k - l; a = b32->a + l;
                }
                l = k;
            }
        }      
        oc[0] = an;
        fprintf(stderr, "-1-[M::%s] oc[0]::%lu, oc[1]::%lu\n", __func__, oc[0], oc[1]);
        if(((oc[0] > (oc[1]*occ_max)) && (oc[0] > (oc[1]-oc[0])) && (oc[1] >= occ_tot) && (oc[0] > 1))) {
            (*rc) = a[0];
            // prt_cigar0((a[0]<<4)>>20, a[0]>>28);
            // prt_bp0((a[0]<<20)>>20, (a[0]<<16)>>28);
            return 1;
        }
    }

    return 0;
}

void push_correct0(window_list *idx, window_list_alloc *res, uint32_t len0, uint32_t rc)
{
    if(len0 != ((uint32_t)-1)) {
        ;
    } else if(rc != ((uint32_t)-1)) {
        uint32_t cc = (rc<<4)>>20, cn = rc>>28, ck = 0, cs, cp;
        uint32_t bc = (rc<<20)>>20, bn = (rc<<16)>>28, bk = 0, bs, bp;
        for (ck = 0; ck < cn; ck++) {
            cs = (cn-1-ck)<<1; cp = (cc>>cs)&3;

            bp = (uint32_t)-1;
            if(cp != 3) {
                bs = (bn-1-bk)<<1; bp = (bc>>bs)&3; 
                bk++;
            }
            push_trace_bp(((asg16_v *)(&(res->c))), cp, bp, 1, ((idx->clen>0)?1:0));

            fprintf(stderr, "%c", cm[(in >> mp)&3]);
        }
    }
}

void push_cns_anchor(overlap_region* ol, All_reads *rref, uint64_t s, uint64_t e, UC_Read* tu, cc_idx_t *idx, overlap_region *aux_o, uint64_t is_tail, uint64_t occ_tot, double occ_max, asg32_v* b32)
{
    if((!is_tail) && (s >= e)) return;//if s >= e && is_tail = 1, gen the cns of the last a few bases -> s = e = ql
    fprintf(stderr, "\n[M::%s] [%lu, %lu)\n", __func__, s, e);
    window_list *p = NULL; uint64_t e0 = 0; uint32_t rc;
    if(aux_o->w_list.n > 0) {
        p = &(aux_o->w_list.a[aux_o->w_list.n-1]);
        e0 = p->x_end+1;
        ///make sure e > s
    }
    assert(s >= e0);
    if((((!is_tail) && (s > 0)) || ((is_tail) && (s > e0)))
        && (cns_gen0(ol, rref, e0, s, tu, idx, occ_tot, occ_max, b32, &rc))) {///CNS in between
        push_correct0(p, &(aux_o->w_list), (uint32_t)-1, rc);
    } else {

    }

    kv_pushp(window_list, aux_o->w_list, &p);
    p->x_start = s; p->x_end = e-1;
    // p->y_start = exz->ps; p->y_end = exz->pe;
}

uint64_t wcns_vote(overlap_region* ol, All_reads *rref, char* qstr, UC_Read* tu, uint64_t *id_a, uint64_t id_n, uint64_t s, uint64_t e, ul_ov_t *c_idx, cc_idx_t *occ, uint64_t occ_tot, double occ_exact, overlap_region *aux_o, asg32_v* b32)
{
    uint64_t k, q[2], rr = 0, os, oe, wl, oc[2], fI; ul_ov_t *p, *gp; overlap_region *z;
    // uint64_t *ct = occ->idx->a;///occ->idx->a[0, wl<<1)
    for (k = 0; k < id_n; k++) {
        p = &(c_idx[id_a[k]]); z = &(ol[ovlp_id(*p)]); 
        q[0] = z->w_list.a[ovlp_cur_wid(*p)].x_start+ovlp_bd(*p);
        q[1] = z->w_list.a[ovlp_cur_wid(*p)].x_end+1-ovlp_bd(*p);
        if(q[1] <= e) rr = 1;
        os = MAX(q[0], s); oe = MIN(q[1], e);
        if(oe > os) {
            ///prepare for CNS
            gp = &(occ->c_idx[id_a[k]]);
            ovlp_cur_xoff(*gp) = ovlp_cur_xoff(*p); ovlp_cur_yoff(*gp) = ovlp_cur_yoff(*p); ovlp_cur_coff(*gp) = ovlp_cur_coff(*p); ovlp_cur_ylen(*gp) = ovlp_cur_ylen(*p);
            // assert(ovlp_cur_wid(*p) == ovlp_cur_wid(*gp));
            // assert(ovlp_id(*p) == ovlp_id(*gp));
            // fprintf(stderr, "[M::%s] tid::%u\t%.*s\twid::%u\tq::[%u, %u)\tos::%lu\toe::%lu\n", __func__, ol[ovlp_id(*p)].y_id, (int)Get_NAME_LENGTH(R_INF, ol[ovlp_id(*p)].y_id), Get_NAME(R_INF, ol[ovlp_id(*p)].y_id),
            // ovlp_cur_wid(*p), ol[ovlp_id(*p)].w_list.a[ovlp_cur_wid(*p)].x_start, ol[ovlp_id(*p)].w_list.a[ovlp_cur_wid(*p)].x_end+1, os, oe);
            extract_sub_cigar_mm(z, os, oe, p, occ->idx->a + os - s);
        }
    }

    wl = e - s;
    os = occ->mms; oe = occ->mme;

    // fprintf(stderr, "[M::%s] s::%lu\te::%lu\n", __func__, s, e);
    
    for (k = 0; k < wl; k++) {
        //+1 for the reference read
        oc[0] = (occ->idx->a[(k<<1)]>>32) + 1;
        oc[1] = ((uint32_t)occ->idx->a[(k<<1)]) + 1;
        // fprintf(stderr, "-0-p::%lu\toc[0]::%lu\toc[1]::%lu\tgoc[0]::%lu\tgoc[1]::%u\n", s + k, oc[0], oc[1], (ct[(k<<1)+1]>>32) + 1, ((uint32_t)ct[(k<<1)+1]) + 1);
        // if(oc[1] < occ_tot || oc[0] <= 1) {
        //     ct[(k<<1)] = ct[(k<<1)+1] = 0;
        //     continue;
        // }
        // fprintf(stderr, "-1-p::%lu\toc[0]::%lu\toc[1]::%lu\n", s + k, oc[0], oc[1]);

        if((oc[0] > (oc[1]*occ_exact)) && (oc[0] > (oc[1]-oc[0])) && (oc[1] >= occ_tot) && (oc[0] > 1)) {
            ///note: there might be insertions at q[k-1, k], insead if q[k, k+1]
            fI = 1;
            ///make sure there is no insertion
            //+1 for the reference read
            oc[0] = (occ->idx->a[(k<<1)+1]>>32) + 1;
            oc[1] = ((uint32_t)occ->idx->a[(k<<1)+1]) + 1;
            if(((oc[0] > (oc[1]*occ_exact)) && (oc[0] > (oc[1]-oc[0])) && (oc[1] >= occ_tot) && (oc[0] > 1))) fI = 0;
            if(fI) {
                // fprintf(stderr, "-1-p::%lu\toc[0]::%lu\toc[1]::%u\tgoc[0]::%lu\tgoc[1]::%u\n", s + k, (occ->idx->a[(k<<1)]>>32) + 1, ((uint32_t)occ->idx->a[(k<<1)]) + 1, (occ->idx->a[(k<<1)+1]>>32) + 1, ((uint32_t)occ->idx->a[(k<<1)+1]) + 1);
                if(oe > os && os != ((uint64_t)-1)) {///push previous intervals
                    push_cns_anchor(ol, rref, os, oe, tu, occ, aux_o, 0, occ_tot, occ_exact, b32);
                }
                os = oe = (uint64_t)-1;
            }

            //+1 for the reference read
            oc[0] = (occ->idx->a[(k<<1)]>>32) + 1;
            oc[1] = ((uint32_t)occ->idx->a[(k<<1)]) + 1;
            if((s+k) == oe) {
                oe++;
            } else {
                if(oe > os && os != ((uint64_t)-1)) {///push previous intervals
                    push_cns_anchor(ol, rref, os, oe, tu, occ, aux_o, 0, occ_tot, occ_exact, b32);
                }
                os = s+k; oe = s+k+1;
            }
        } else {
            // fprintf(stderr, "-2-p::%lu\toc[0]::%lu\toc[1]::%u\tgoc[0]::%lu\tgoc[1]::%u\n", s + k, (occ->idx->a[(k<<1)]>>32) + 1, ((uint32_t)occ->idx->a[(k<<1)]) + 1, (occ->idx->a[(k<<1)+1]>>32) + 1, ((uint32_t)occ->idx->a[(k<<1)+1]) + 1);
            if(oe > os && os != ((uint64_t)-1)) {///push previous intervals
                push_cns_anchor(ol, rref, os, oe, tu, occ, aux_o, 0, occ_tot, occ_exact, b32);
            }
            os = oe = (uint64_t)-1;
        }
        occ->idx->a[(k<<1)] = occ->idx->a[(k<<1)+1] = 0;
    }

    occ->mms = occ->mme = (uint64_t)-1;
    if(oe > os && os != ((uint64_t)-1)) {
        occ->mms = os; occ->mme = oe;
    }
    return rr;
}



void print_debug_ovlp_cigar(overlap_region_alloc* ol, asg64_v* idx, kv_ul_ov_t *c_idx)
{
    uint64_t k, ci; uint32_t cl; ul_ov_t *cp; bit_extz_t ez; uint16_t c; char cm[4]; 
    cm[0] = 'M'; cm[1] = 'S'; cm[2] = 'I'; cm[3] = 'D';  
    for (k = 0; k < idx->n; k++) {
        cp = &(c_idx->a[(uint32_t)idx->a[k]]);
        fprintf(stderr, "**********[M::%s] tid::%u\t%.*s\twid::%u\tq::[%u, %u)\terr::%d\toerr::%u**********\n", __func__, ol->list[ovlp_id(*cp)].y_id, (int)Get_NAME_LENGTH(R_INF, ol->list[ovlp_id(*cp)].y_id), Get_NAME(R_INF, ol->list[ovlp_id(*cp)].y_id),
        ovlp_cur_wid(*cp), ol->list[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_start, ol->list[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_end+1, ol->list[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].error, ol->list[ovlp_id(*cp)].non_homopolymer_errors);
        set_bit_extz_t(ez, ol->list[ovlp_id(*cp)], ovlp_cur_wid(*cp)); ci = 0;
        while (ci < ez.cigar.n) {
            ci = pop_trace(&(ez.cigar), ci, &c, &cl);
            fprintf(stderr, "%u%c", cl, cm[c]);
        }
        fprintf(stderr, "\n");
    }    
}

void wcns_gen(overlap_region_alloc* ol, All_reads *rref, UC_Read* qu, UC_Read* tu, kv_ul_ov_t *c_idx, asg64_v* idx, asg64_v* buf, int64_t bd, uint64_t wl, int64_t ql, uint64_t occ_tot, double occ_exact, overlap_region *aux_o, asg32_v* b32)
{
    int64_t on = ol->length, k, i, zwn, q[2]; 
    uint64_t m, *ra, rn; overlap_region *z; ul_ov_t *cp;

    for (k = idx->n = c_idx->n = 0; k < on; k++) {
        z = &(ol->list[k]); zwn = z->w_list.n; 
        if((!zwn) || (z->is_match != 1)) continue;
        for (i = 0; i < zwn; i++) {
            if(is_ualn_win(z->w_list.a[i])) continue;
            q[0] = z->w_list.a[i].x_start; q[1] = z->w_list.a[i].x_end;
            q[0] += bd; q[1] -= bd;
            if(q[1] >= q[0]) {
                m = ((uint64_t)q[0]); m <<= 32; 
                m += c_idx->n; kv_push(uint64_t, *idx, m);

                kv_pushp(ul_ov_t, *c_idx, &cp);
                ovlp_id(*cp) = k; ///ovlp id
                // ovlp_min_wid(*cp) = i; ///beg id of windows
                // ovlp_max_wid(*cp) = i; ///end id of windows
                ovlp_cur_wid(*cp) = i; ///cur id of windows
                ovlp_cur_xoff(*cp) = z->w_list.a[i].x_start; ///cur xpos
                ovlp_cur_yoff(*cp) = z->w_list.a[i].y_start; ///cur xpos
                ovlp_cur_ylen(*cp) = 0;
                ovlp_cur_coff(*cp) = 0; ///cur cigar off in cur window
                ovlp_bd(*cp) = bd;
            }
        }
    }

    int64_t srt_n = idx->n, s, e, t, rr; i = 0;
    radix_sort_ec64(idx->a, idx->a+idx->n);
    for (k = 1, i = 0; k < srt_n; k++) {
        if (k == srt_n || (idx->a[k]>>32) != (idx->a[i]>>32)) {
            if(k - i > 1) {
                for (t = i; t < k; t++) {
                    cp = &(c_idx->a[(uint32_t)idx->a[t]]);
                    // s = ol->list[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_start+ovlp_bd(*cp);
                    // assert(s == (int64_t)(idx->a[i]>>32));
                    m = ol->list[ovlp_id(*cp)].w_list.a[ovlp_cur_wid(*cp)].x_end+1-ovlp_bd(*cp);
                    m <<= 32; m += ((uint32_t)idx->a[t]); idx->a[t] = m;
                    // fprintf(stderr, "[M::%s] s::%ld\tsi::%lu\n", __func__, s, (idx->a[i]>>32));
                }
                radix_sort_ec64(idx->a + i, idx->a + k);
            }
            i = k;
        }
    }

    print_debug_ovlp_cigar(ol, idx, c_idx);

    ///second index
    kv_resize(ul_ov_t, *c_idx, (c_idx->n<<1));
    ul_ov_t *idx_a = NULL, *idx_b = NULL;
    idx_a = c_idx->a; idx_b = c_idx->a;
    memcpy(idx_b, idx_a, c_idx->n * (sizeof((*(idx_a)))));

    kv_resize(uint64_t, *buf, ((wl<<1) + idx->n)); buf->n = ((wl<<1) + idx->n);
    memcpy(buf->a + (wl<<1), idx->a, idx->n * (sizeof((*(idx->a)))));
    memset(buf->a, 0, (wl<<1)*(sizeof((*(idx->a)))));

    cc_idx_t ii_a, ii_b; memset(&ii_a, 0, sizeof(ii_a)); memset(&ii_b, 0, sizeof(ii_b));

    ii_a.c_idx = idx_a; ii_a.idx = idx; ii_a.i = ii_a.i0 = 0; ii_a.srt_n = ii_a.idx->n; ii_a.mms = ii_a.mme = (uint64_t)-1;
    ii_b.c_idx = idx_b; ii_b.idx = buf; ii_b.i = ii_b.i0 = (wl<<1); ii_b.srt_n = ii_b.idx->n; ii_b.mms = ii_b.mme = (uint64_t)-1;

    s = 0; e = wl; e = ((e<=ql)?e:ql); rr = 0; 
    aux_o->w_list.n = aux_o->w_list.c.n = 0; ///for cigar
    for (; s < ql; ) {
        rn = iter_cc_idx_t(ol->list, &ii_a, s, e, rr, 0, &ra);
        // debug_inter0(ol->list, ii_a.c_idx, ii_a.idx->a + ii_a.i0, ii_a.srt_n - ii_a.i0, ra, rn, s, e, 0, 1, "-0-");
        rr = wcns_vote(ol->list, rref, qu->seq, tu, ra, rn, s, e, ii_a.c_idx, &ii_b, occ_tot, occ_exact, aux_o, b32);
        s += wl; e += wl; e = ((e<=ql)?e:ql);
    }

    if(ii_b.mme > ii_b.mms && ii_b.mms != (uint64_t)-1) {
        push_cns_anchor(ol->list, rref, ii_b.mms, ii_b.mme, tu, &ii_b, aux_o, 0, occ_tot, occ_exact, b32);
    }

    push_cns_anchor(ol->list, rref, ql, ql, tu, &ii_b, aux_o, 1, occ_tot, occ_exact, b32);
}

static void worker_hap_ec(void *data, long i, int tid)
{
	ec_ovec_buf_t0 *b = &(((ec_ovec_buf_t*)data)->a[tid]);
    uint32_t high_occ = asm_opt.hom_cov * (2.0 - HA_KMER_GOOD_RATIO);
    uint32_t low_occ = asm_opt.hom_cov * HA_KMER_GOOD_RATIO;
    overlap_region *aux_o = NULL; asg64_v buf0;

    // if (memcmp("m64012_190920_173625/7210046/ccs", Get_NAME((R_INF), i), Get_NAME_LENGTH((R_INF),i)) == 0) {
    //     fprintf(stderr, "[M::%s-beg] rid->%ld\n", __func__, i);
    // } else {
    //     return;
    // }

    if(i != 596/**1024**/) return;

    recover_UC_Read(&b->self_read, &R_INF, i);

    h_ec_lchain(b->ab, i, b->self_read.seq, b->self_read.length, asm_opt.mz_win, asm_opt.k_mer_length, &R_INF, &b->olist, &b->clist, 0.02, asm_opt.max_n_chain, 1, NULL, NULL, &(b->sp), &high_occ, &low_occ, 1, 1, 0, 2, 2, UINT32_MAX);

    b->num_read_base += b->olist.length;
    aux_o = fetch_aux_ovlp(&b->olist);///must be here

    gen_hc_r_alin(&b->olist, &b->clist, &R_INF, &b->self_read, &b->ovlp_read, &b->exz, aux_o, asm_opt.max_ov_diff_ec, WINDOW_HC, i, E_KHIT/**asm_opt.k_mer_length**/, 1, &b->v16);

    // fprintf(stderr, "\n[M::%s] rid::%ld\t%.*s\tlen::%lld\tocc::%lu\n", __func__, i, (int)Get_NAME_LENGTH(R_INF, i), 
    //             Get_NAME(R_INF, i), b->self_read.length, b->olist.length);

    b->num_correct_base += b->olist.length;

    copy_asg_arr(buf0, b->sp);
    rphase_hc(&b->olist, &R_INF, &b->hap, &b->self_read, &b->ovlp_read, &b->pidx, &b->v64, &buf0, 0, WINDOW_MAX_SIZE, b->self_read.length, 1);
    copy_asg_arr(b->sp, buf0);

    copy_asg_arr(buf0, b->sp);
    wcns_gen(&b->olist, &R_INF, &b->self_read, &b->ovlp_read, &b->pidx, &b->v64, &buf0, 0, 512, b->self_read.length, 3, 0.500001, aux_o, &b->v32);
    copy_asg_arr(b->sp, buf0);


    uint32_t k;
    for (k = 0; k < b->olist.length; k++) {
        if(b->olist.list[k].is_match == 1) b->num_recorrect_base++;
    }

    // exit(1);
    
    

    // prt_chain(&b->olist);

    // ul_map_lchain(b->abl, (uint32_t)-1, s->seq[i], s->len[i], s->opt->w, s->opt->k, s->uu, &b->olist, &b->clist, s->opt->bw_thres, 
    //         s->opt->max_n_chain, 1, NULL, &(b->tmp_region), NULL, &(b->sp), &high_occ, NULL, 0, 1, 0.2/**0.75**/, 2, 3);

    /**
	int fully_cov, abnormal;
    // if(i != 12578) return;
    // fprintf(stderr, "[M::%s-beg] rid->%ld\n", __func__, i);
    // if (memcmp("7897e875-76e5-42c8-bc37-94b370c4cc8d", Get_NAME((R_INF), i), Get_NAME_LENGTH((R_INF),i)) == 0) {
    //     fprintf(stderr, "[M::%s-beg] rid->%ld\n", __func__, i);
    // } else {
    //     return;
    // }

    ha_get_candidates_interface(b->ab, i, &b->self_read, &b->olist, &b->olist_hp, &b->clist, 
    0.02, asm_opt.max_n_chain, 1, NULL, &b->r_buf, &(R_INF.paf[i]), &(R_INF.reverse_paf[i]), &(b->tmp_region), NULL, &(b->sp));

	clear_Cigar_record(&b->cigar1);
	clear_Round2_alignment(&b->round2);

	correct_overlap(&b->olist, &R_INF, &b->self_read, &b->correct, &b->ovlp_read, &b->POA_Graph, &b->DAGCon,
			&b->cigar1, &b->hap, &b->round2, &b->r_buf, &(b->tmp_region.w_list), 0, 1, &fully_cov, &abnormal);

	b->num_read_base += b->self_read.length;
	b->num_correct_base += b->correct.corrected_base;
	b->num_recorrect_base += b->round2.dumy.corrected_base;

	push_cigar(R_INF.cigars, i, &b->cigar1);
	push_cigar(R_INF.second_round_cigar, i, &b->round2.cigar);

	R_INF.paf[i].is_fully_corrected = 0;
	if (fully_cov) {
		if (get_cigar_errors(&b->cigar1) == 0 && get_cigar_errors(&b->round2.cigar) == 0)
			R_INF.paf[i].is_fully_corrected = 1;
	}
	R_INF.paf[i].is_abnormal = abnormal;

    R_INF.trio_flag[i] = AMBIGU;
    
    ///need to be fixed in r305
    // if(ha_idx_hp == NULL)
    // {
    //     R_INF.trio_flag[i] += collect_hp_regions(&b->olist, &R_INF, &(b->k_flag), RESEED_HP_RATE, Get_READ_LENGTH(R_INF, i), NULL);
    // }

    if (R_INF.trio_flag[i] != AMBIGU || b->save_ov) {
		int is_rev = (asm_opt.number_of_round % 2 == 0);
		push_overlaps(&(R_INF.paf[i]), &b->olist, 1, &R_INF, is_rev);
		push_overlaps(&(R_INF.reverse_paf[i]), &b->olist, 2, &R_INF, is_rev);
	}

    if(het_cnt) het_cnt[i] = get_het_cnt(&b->hap);
    // fprintf(stderr, "[M::%s-end] rid->%ld\n", __func__, i);
    **/
}

void cal_ec_multiple(ec_ovec_buf_t *b, uint64_t n_thre, uint64_t n_a)
{
    double tt0 = yak_realtime_0();
    kt_for(n_thre, worker_hap_ec, b, n_a);///debug_for_fix
    fprintf(stderr, "[M::%s-reads] #->%lu\n", __func__, n_a);
    fprintf(stderr, "[M::%s::%.3f] ==> chaining\n", __func__, yak_realtime_0()-tt0);
}