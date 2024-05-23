#ifndef __HASHTABLE__
#define __HASHTABLE__
#include "htab.h"

#define PREFIX_BITS 16
#define MAX_SUFFIX_BITS 64
#define MODE_VALUE 101

#define WINDOW 375
#define WINDOW_BOUNDARY 375
///for one side, the first or last WINDOW_UNCORRECT_SINGLE_SIDE_BOUNDARY bases should not be corrected
#define WINDOW_UNCORRECT_SINGLE_SIDE_BOUNDARY 25
#define THRESHOLD  15
#define OVERLAP_THRESHOLD_HIFI_FILTER 0.9
#define OVERLAP_THRESHOLD_NOSI_FILTER 0.7
#define OVERLAP_THRESHOLD_FILTER_HPC 0.75
#define HIGH_HET_OVERLAP_THRESHOLD_FILTER 0.3
#define HIGH_HET_ERROR_RATE 0.08
#define THRESHOLD_MAX_SIZE  31
#define THRESHOLD_UL_MAX 0.2
#define WINDOW_UL 75
#define WINDOW_UL_H 200
// #define WINDOW_UL_H 150
#define MIN_UL_ALIN_RATE 0.5
#define MIN_UL_ALIN_LEN (WINDOW_UL*6)
#define WINDOW_UL_BOUND 48
#define WINDOW_UL_BOUND_RATE 0.55

#define GROUP_SIZE 4
///the max cigar likes 10M10D10M10D10M
///#define CIGAR_MAX_LENGTH THRESHOLD*2+2
#define CIGAR_MAX_LENGTH 31*2+4

typedef struct
{
    uint32_t offset;
    uint32_t readID:31, rev:1;
} k_mer_pos;

typedef struct
{
    k_mer_pos* list;
    uint64_t length;
    uint64_t size;
    uint8_t direction;
    uint64_t end_pos;
} k_mer_pos_list;

typedef struct
{
  ///the begining and end of a window, instead of the whole overlap
  int32_t x_start, x_end;
  int32_t y_start, y_end;
  int16_t extra_begin, extra_end;
  int16_t error, error_threshold;
  uint32_t cidx, clen;
} window_list;

typedef struct
{
    size_t n, m; 
    window_list *a;
    kvec_t(uint16_t) c;
} window_list_alloc;

typedef struct
{
    uint64_t* buffer;
	uint32_t length;
	uint32_t size;
} Fake_Cigar;

typedef struct
{
    uint32_t x_id;
    ///the begining and end of the whole overlap
    uint32_t x_pos_s;
    uint32_t x_pos_e;
    uint32_t x_pos_strand;

    uint32_t y_id;
    uint32_t y_pos_s;
    uint32_t y_pos_e;
    uint32_t y_pos_strand;

    uint32_t overlapLen;
    int32_t shared_seed;
    uint32_t align_length;
    uint8_t is_match;
    uint8_t without_large_indel;
    int8_t strong;
    uint32_t non_homopolymer_errors;

    // window_list* w_list;
    // uint32_t w_list_size;
    // uint32_t w_list_length;
    Fake_Cigar f_cigar;

    window_list_alloc w_list;
    window_list_alloc boundary_cigars;
} overlap_region;

typedef struct
{
    overlap_region* list;
    uint64_t size;
    uint64_t length;
    int64_t mapped_overlaps_length;
} overlap_region_alloc;

typedef struct
{
	uint32_t readID:31, strand:1;
	uint32_t offset, self_offset, cnt;
} k_mer_hit;

typedef struct {
	int32_t *score;
	int64_t *pre;
	int32_t *indels;
	int32_t *self_length;
    int32_t *occ;
	int64_t *tmp; // MUST BE 64-bit integer
	int64_t length;
	int64_t size;
} Chain_Data;

typedef struct
{
    k_mer_hit* list;
    long long length;
    long long size;
    Chain_Data chainDP;
} Candidates_list;

void init_Candidates_list(Candidates_list* l);
void clear_Candidates_list(Candidates_list* l);
void destory_Candidates_list(Candidates_list* l);
void destory_Candidates_list_buf(void *km, Candidates_list* l, int is_z);

void init_overlap_region_alloc(overlap_region_alloc* list);
void clear_overlap_region_alloc(overlap_region_alloc* list);
void destory_overlap_region_alloc(overlap_region_alloc* list);
void append_window_list(overlap_region* region, uint64_t x_start, uint64_t x_end, int y_start, int y_end, int error,
int extra_begin, int extra_end, int error_threshold, int blockLen, void *km);

void overlap_region_sort_y_id(overlap_region *a, long long n);

void calculate_overlap_region_by_chaining(Candidates_list* candidates, overlap_region_alloc* overlap_list, kvec_t_u64_warp* chain_idx,
uint64_t readID, uint64_t readLength, All_reads* R_INF, const ul_idx_t *uref, double band_width_threshold, int add_beg_end, overlap_region* f_cigar, void *km);

void init_fake_cigar(Fake_Cigar* x);
void destory_fake_cigar(Fake_Cigar* x);
void clear_fake_cigar(Fake_Cigar* x);
void add_fake_cigar(Fake_Cigar* x, uint32_t gap_site, int32_t gap_shift, void *km);
void resize_fake_cigar(Fake_Cigar* x, uint64_t size, void *km);
int get_fake_gap_pos(Fake_Cigar* x, int index);
int get_fake_gap_shift(Fake_Cigar* x, int index);

static inline long long y_start_offset(long long x_start, Fake_Cigar* o)
{
    if(x_start == get_fake_gap_pos(o, o->length - 1))
    {
        return get_fake_gap_shift(o, o->length - 1);
    }
    
    long long i;
    for (i = 0; i < (long long)o->length; i++)
    {
        if(x_start < get_fake_gap_pos(o, i))
        {
            break;
        }
    }

    if(i == 0 || i == (long long)o->length)
    {
        fprintf(stderr, "ERROR at %s:%d\n", __FILE__, __LINE__);
        exit(0);
    }

    ///note here return i - 1
    return get_fake_gap_shift(o, i - 1);
}

void resize_Chain_Data(Chain_Data* x, long long size, void *km);
void init_window_list_alloc(window_list_alloc* x);
void clear_window_list_alloc(window_list_alloc* x);
void destory_window_list_alloc(window_list_alloc* x);
void resize_window_list_alloc(window_list_alloc* x, uint64_t size);
long long chain_DP(k_mer_hit* a, long long a_n, Chain_Data* dp, overlap_region* result, double band_width_threshold, int max_skip, int x_readLen, int y_readLen, void *km);
uint64_t lchain_dp(k_mer_hit* a, int64_t a_n, k_mer_hit* des, Chain_Data* dp, overlap_region* res, 
              int64_t max_skip, int64_t max_iter, int64_t max_dis, double chn_pen_gap, double chn_pen_skip, double bw_rate, 
              int64_t xl, int64_t yl, int64_t quick_check);
int ovlp_chain_gen(overlap_region_alloc* ol, overlap_region* t, int64_t xl, int64_t yl, int64_t apend_be, k_mer_hit* hit, int64_t n_hit);
void gen_fake_cigar(Fake_Cigar* z, overlap_region *o, int64_t apend_be, k_mer_hit* hit, int64_t n_hit);
int append_utg_inexact_overlap_region_alloc(overlap_region_alloc* list, overlap_region* tmp, 
                                        ma_utg_v *ua, int add_beg_end, void *km);

#define kv_pushp_cl(type, v, p) do {									\
		if ((v).length == (v).size) {										\
			(v).size = (v).size? (v).size<<1 : 2;							\
			(v).list = (type*)realloc((v).list, sizeof(type) * (v).size);	\
		}															\
		*(p) = &((v).list[(v).length++]); \
	} while (0)

#define kv_resize_cl(type, v, s) do { \
		if ((v).size < (s)) { \
			(v).size = (s); \
			kv_roundup32((v).size); \
			(v).list = (type*)realloc((v).list, sizeof(type) * (v).size); \
		} \
	} while (0)

#define is_alnw(a) (((a).readID) == ((uint32_t)(0x7fffffff)))
#define is_pri_aln(a) ((((a).readID) == ((uint32_t)(0x7fffffff)))||((a).cnt >= (a).readID))

uint64_t lchain_dp_trace(k_mer_hit* a, int64_t a_n, int64_t max_lgap, double sgap_rate, int64_t sgap);
uint64_t lchain_qdp(k_mer_hit* a, int64_t a_n, k_mer_hit* des, Chain_Data* dp, overlap_region* res, 
              int64_t max_skip, int64_t max_iter, int64_t max_dis, double chn_pen_gap, double chn_pen_skip, double bw_rate, 
              int64_t xl, int64_t yl, int64_t quick_check);
int ovlp_chain_qgen(overlap_region_alloc* ol, overlap_region* t, int64_t xl, int64_t yl, int64_t apend_be, k_mer_hit* hit, int64_t n_hit);
uint64_t lchain_refine(k_mer_hit* a, int64_t a_n, k_mer_hit* des, Chain_Data* dp, 
                                    int64_t max_skip, int64_t max_iter, int64_t max_dis, int64_t long_gap);
uint64_t lchain_qdp_fix(k_mer_hit* a, int64_t a_n, Chain_Data* dp, int64_t max_skip, 
                int64_t max_iter, int64_t max_dis, double chn_pen_gap, double chn_pen_skip, 
                double bw_rate, int64_t xl, int64_t yl, int64_t quick_check, 
                int64_t left_fix, int64_t right_fix);
uint64_t lchain_simple(k_mer_hit* a, int64_t a_n, k_mer_hit* des, Chain_Data* dp, 
                                                        int64_t max_skip, int64_t max_iter);
uint64_t lchain_qdp_mcopy(Candidates_list *cl, int64_t a_idx, int64_t a_n, int64_t des_idx, 
              Chain_Data* dp, overlap_region_alloc* res, int64_t max_skip, int64_t max_iter, 
              int64_t max_dis, double chn_pen_gap, double chn_pen_skip, double bw_rate, 
              uint32_t xid, int64_t xl, int64_t yl, int64_t quick_check, uint32_t apend_be, 
              int64_t gen_cigar, int64_t enable_mcopy, double mcopy_rate, int64_t mcopy_khit_cutoff, 
              int64_t khit_n);

#endif
