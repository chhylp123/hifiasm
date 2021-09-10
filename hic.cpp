#define __STDC_LIMIT_MACROS
#include "float.h"
#include <math.h>
#include "hic.h"
#include "htab.h"
#include "assert.h"
#include "Overlaps.h"
#include "Hash_Table.h"
#include "Correct.h"
#include "Purge_Dups.h"
#include "rcut.h"
#include "khashl.h"
#include "kthread.h"
#include "ksort.h"
#include "kseq.h" // FASTA/Q parser
#include "kdq.h"
#include "horder.h"
KSEQ_INIT(gzFile, gzread)
KDQ_INIT(uint64_t)


#define OFFSET_RATE 0.000000001
#define OFFSET_SECOND_RATE 0.0000000001
#define SCALL 10000
#define OFFSET_RATE_MAX_W 6.90675477865*SCALL
#define OFFSET_RATE_MIN_W 4.0000003e-10*SCALL

#define HIC_COUNTER_BITS 12
#define HIC_MAX_COUNT    ((1<<HIC_COUNTER_BITS)-1)
#define HIC_KEY_MODE ((uint64_t)(((uint64_t)-1)-HIC_MAX_COUNT))
#define HIC_R_E_RATE 0.01

const unsigned char b2rc[5] = {'T', 'G', 'C', 'A', 'N'};

#define hic_ct_eq(a, b) ((a)>>HIC_COUNTER_BITS == (b)>>HIC_COUNTER_BITS) 
#define hic_ct_hash(a) ((a)>>HIC_COUNTER_BITS)
KHASHL_MAP_INIT(static klib_unused, hc_pt_t, hc_pt, uint64_t, uint64_t, hic_ct_hash, hic_ct_eq)

#define u_trans_m_key(a) (((uint64_t)((a).qn)<<32) | ((uint64_t)((a).tn)))
KRADIX_SORT_INIT(u_trans_m, u_trans_t, u_trans_m_key, 8)

#define u_trans_occ_key(a) ((a).occ)
KRADIX_SORT_INIT(u_trans_occ, u_trans_t, u_trans_occ_key, member_size(u_trans_t, occ))

#define is_hom_hit(a) ((a).id == (uint64_t)-1)
#define HC_PT_MA 65

typedef struct {
    kv_gg_status sg;
    uint64_t xs;
} psg_t;

typedef struct{
    kvec_t(char) name;
    kvec_t(uint64_t) name_Len;
    kvec_t(char) r;
    kvec_t(uint64_t) r_Len;
    uint64_t idx;
} reads_t;


typedef struct{
    kvec_t(hc_edge_warp) rGraph;
    kvec_t(uint64_t) order;
    pdq pq;
    kvec_t(uint8_t) rGraphSet;
    kvec_t(uint8_t) rGraphVis;
    kvec_t(uint8_t) utgVis;
    kvec_t(uint8_t) bmerVis;
    kdq_t(uint64_t) *q;
    kvec_t(uint32_t) parent;
    kvec_t(double) p_weight;
    const uint64_t* enzymes;
    uint64_t uID_mode, uID_shift, n, src, dest, n_e, c_e;
    int p_mer, a_mer, b_mer;
} min_cut_t;


typedef struct{
    kvec_t(uint32_t) a;
    uint32_t h[2];
    uint8_t full_bub;
    int status[2];
    double weight[2], weight_convex;
}partition_warp;

typedef struct{
    size_t n, m; 
    partition_warp* a;
    uint32_t* index;
}G_partition;

typedef struct{
    kvec_t(uint8_t) vis;
    double weight;
    long long bid, uid, chainID;
}block_phase_type;

typedef struct{
    uint64_t n; 
    uint8_t* lock;
    uint32_t* hap;
    uint32_t m[3];
    uint32_t label, label_add, label_shift;
    hc_links* link;
    G_partition g_p;
    G_partition group_g_p;
    kvec_t(double) label_buffer;
    block_phase_type b;
}H_partition;

typedef struct {
	uint32_t p; // the optimal parent vertex
	uint32_t d; // the shortest distance from the initial vertex
	uint32_t nc; // max count of reads, no matter positive or negative
    double nh, w[2];
    uint32_t uc, ac; // used vertex/allowed vertex
	uint32_t r:31, s:1; // r: the number of remaining incoming arc; s: state
	//s: state, s=0, this edge has not been visited, otherwise, s=1
} bub_p_t;

typedef struct {
	///all information for each node
	bub_p_t *a;
	kvec_t(uint32_t) S; // set of vertices without parents, nodes with all incoming edges visited
	kvec_t(uint32_t) T; // set of tips
	kvec_t(uint32_t) b; // visited vertices
	kvec_t(uint32_t) e; // visited edges/arcs
    uint32_t exist_hap_label;
} bub_p_t_warp;


typedef struct {
	hc_pt_t *h;
	uint64_t n;
	uint64_t *a;
    khint_t end;///end of total idx
} hc_pt1_t;

typedef struct {
	ma_ug_t* ug;
    asg_t* read_g;
    ///hc_links* link;
    trans_chain* t_ch;
    uint64_t uID_bits;
    uint64_t uID_mode;
    uint64_t pos_bits;
    uint64_t pos_mode;
    uint64_t rev_mode;
    uint64_t k;
    uint64_t hap_cnt;


    uint64_t pre;
    uint64_t tot;
    uint64_t tot_pos;
    uint64_t up_bound, low_bound;
    hc_pt1_t* idx_buf;
    long double a, b, frac, max_d;
} ha_ug_index;

typedef struct { // data structure for each step in kt_pipeline()
    uint64_t key, pos;
} ch_buf_t;

typedef struct {
	kvec_t(uint64_t) a;
} kvec_cnt;

typedef struct {
	kvec_t(ch_buf_t) a;
} kvec_pos;

typedef struct { // global data structure for kt_pipeline()
	int is_cnt;
    uint64_t buf_bytes;
	ha_ug_index *h;
    kvec_cnt* cnt;
	kvec_pos* buf;
    uint64_t n_thread;
} pldat_t;

typedef struct {
	uint64_t *a, id;
    uint16_t occ1, occ2;
} pe_hit_hap;

typedef struct {
    pe_hit_hap* a;
    size_t n, m;
    uint64_t n_u;
} kvec_pe_hit_hap;

typedef struct {
    kvec_t(hc_edge) a;
}kvec_hc_edge;

#define pe_hit_an1_key(x) ((x).s)
KRADIX_SORT_INIT(pe_hit_an1, pe_hit, pe_hit_an1_key, member_size(pe_hit, s))
#define pe_hit_an2_key(x) ((x).e)
KRADIX_SORT_INIT(pe_hit_an2, pe_hit, pe_hit_an2_key, member_size(pe_hit, e))

#define pe_hit_an1_idx_key(x) ((x).s<<1)
KRADIX_SORT_INIT(pe_hit_idx_an1, pe_hit, pe_hit_an1_idx_key, member_size(pe_hit, s))
#define pe_hit_an2_idx_key(x) ((x).e<<1)
KRADIX_SORT_INIT(pe_hit_idx_an2, pe_hit, pe_hit_an2_idx_key, member_size(pe_hit, e))

#define generic_key(x) (x)
KRADIX_SORT_INIT(hc64, uint64_t, generic_key, 8)
KRADIX_SORT_INIT(u32, uint32_t, generic_key, 4)
#define g_partition_key(x) (((x)>>1)+((x)<<63))
KRADIX_SORT_INIT(g_partition, uint64_t, g_partition_key, 8)

#define get_pe_s(x) ((x).a[0])
#define get_pe_e(x) ((x).a[(x).occ1])
KRADIX_SORT_INIT(pe_an1, pe_hit_hap, get_pe_s, 8)
KRADIX_SORT_INIT(pe_an2, pe_hit_hap, get_pe_e, 8)

#define pe_occ_key_1(x) ((x).occ1)
KRADIX_SORT_INIT(pe_occ1, pe_hit_hap, pe_occ_key_1, member_size(pe_hit_hap, occ1))
#define pe_occ_key_2(x) ((x).occ2)
KRADIX_SORT_INIT(pe_occ2, pe_hit_hap, pe_occ_key_2, member_size(pe_hit_hap, occ2))
#define pe_occ_key_t(x) (((uint64_t)((x).occ1))+((uint64_t)((x).occ2)))
KRADIX_SORT_INIT(pe_occ_t, pe_hit_hap, pe_occ_key_t, 8)

#define asg_arc_key(a) ((a).ul)
KRADIX_SORT_INIT(asg_e, asg_arc_t, asg_arc_key, 8)

typedef struct { // global data structure for kt_pipeline()
	const ha_ug_index* idx;
	kseq_t *ks1, *ks2;
    int64_t chunk_size;
    uint64_t n_thread;
    uint64_t total_base;
    uint64_t total_pair;
    kvec_pe_hit hits;
    ///kvec_pe_hit_hap hits;
    trans_chain* t_ch;
} sldat_t;

typedef struct {
	uint64_t ref;
    uint64_t off_cnt; 
} s_hit;

typedef struct {
    kvec_t(s_hit) a;
} kvec_vote;

typedef struct { // data structure for each step in kt_pipeline()
    const ha_ug_index* idx;
	int n, m, sum_len;
	uint64_t *len, id;
	char **seq;
	ch_buf_t *buf;
    kvec_vote* pos_buf;
    pe_hit* pos;
    ///pe_hit_hap* pos;
    trans_chain* t_ch;
} stepdat_t;

#define generic_key(x) (x)
KRADIX_SORT_INIT(b64, uint64_t, generic_key, 8)
#define ch_buf_t_key(a) ((a).key)
KRADIX_SORT_INIT(ch_buf, ch_buf_t, ch_buf_t_key, member_size(ch_buf_t, key))
#define hc_pos_key(x) ((x)<<1)
KRADIX_SORT_INIT(hc_pos, uint64_t, hc_pos_key, 8)
#define hc_s_hit_an1_key(a) ((a).ref)
KRADIX_SORT_INIT(hc_s_hit_an1, s_hit, hc_s_hit_an1_key, 8)
#define hc_s_hit_an2_key(a) ((uint32_t)(a).off_cnt)
KRADIX_SORT_INIT(hc_s_hit_an2, s_hit, hc_s_hit_an2_key, 8)
#define hc_s_hit_off_cnt_key(a) ((a).off_cnt)
KRADIX_SORT_INIT(hc_s_hit_off_cnt, s_hit, hc_s_hit_off_cnt_key, 8)
#define hc_edge_key_u(a) ((a).uID)
KRADIX_SORT_INIT(hc_edge_u, hc_edge, hc_edge_key_u, 4)
#define hc_edge_key_d(a) ((a).dis)
KRADIX_SORT_INIT(hc_edge_d, hc_edge, hc_edge_key_d, member_size(hc_edge, dis))

#define k_trans_qs_key(a) ((a).qs)
KRADIX_SORT_INIT(k_trans_qs, u_trans_t, k_trans_qs_key, member_size(u_trans_t, qs))

#define get_hit_suid(x, k) (((x).a.a[(k)].s<<1)>>(64 - (x).uID_bits))
#define get_hit_spos(x, k) ((x).a.a[(k)].s & (x).pos_mode)
#define get_hit_euid(x, k) (((x).a.a[(k)].e<<1)>>(64 - (x).uID_bits))
#define get_hit_epos(x, k) ((x).a.a[(k)].e & (x).pos_mode)


typedef struct {
    kvec_t(kvec_t_u64_warp) matrix;
    uint64_t uID_shift, dis_mode;
} MT;

typedef struct{
    uint64_t beg, end, dis, cnt_0, cnt_1;
} trans_p_t;

typedef struct{
    trans_p_t* a;
    size_t n, m;
    uint64_t max, med;
} trans_idx;

reads_t R1, R2;
ha_ug_index* ug_index;

void print_debug_bubble_graph(bubble_type* bub, ma_ug_t* ug, const char *fn);

void build_bub_graph(ma_ug_t* ug, bubble_type* bub);

void init_ha_ug_index_opt(ha_ug_index* idx, ma_ug_t *ug, int k, pldat_t* p, uint64_t up_occ, 
uint64_t low_occ, uint64_t thread_num)
{
    uint64_t i, n;
    for (idx->uID_bits=1; (uint64_t)(1<<idx->uID_bits)<(uint64_t)ug->u.n; idx->uID_bits++);
    idx->pos_bits = 64 - idx->uID_bits - 1;
    idx->uID_mode = (((uint64_t)-1) << (64-idx->uID_bits))>>1;
    idx->pos_mode = ((uint64_t)-1) >> (64-idx->pos_bits);
    idx->rev_mode = ((uint64_t)1) << 63;
    idx->ug = ug;
    idx->k = k;
    idx->pre = HIC_COUNTER_BITS;
    idx->tot = 1 << idx->pre;
    idx->tot_pos = 0;
    ///idx->up_bound = 1;
    idx->up_bound = up_occ;
    idx->low_bound = low_occ;
    CALLOC(idx->idx_buf, idx->tot);
    for (i = 0; i < idx->tot; i++)
    {
        idx->idx_buf[i].h = hc_pt_init();
    }
    for (i = n = 0; i < ug->u.n; i++)
    {
        n += ug->u.a[i].len;
    }
    n = n << 3;
    p->h = idx;
    p->buf_bytes = n>>7;
    CALLOC(p->cnt, idx->tot);
    CALLOC(p->buf, idx->tot);
    for (i = 0; i < idx->tot; i++)
    {
        kv_init(p->cnt[i].a);
        kv_init(p->buf[i].a);
    }
    p->n_thread = thread_num;
}

inline uint64_t get_k_direction(uint64_t x[4])
{
    if(x[1] != x[3]) 
    {
        return x[1] < x[3]? 0 : 1;
    }
    else if(x[0] != x[2])
    {
        return x[0] < x[2]? 0 : 1;
    }
    else
    {
        return (uint64_t)-1;
    }
}

inline uint64_t hc_hash_long(uint64_t x[4], uint64_t* skip, uint64_t k)
{
	///compare forward k-mer and reverse complementary strand
	(*skip) = get_k_direction(x);
    if((*skip) == (uint64_t)-1) return (*skip);
    if (k <= 32) return ((x[(*skip)<<1|0]<<32)|(x[(*skip)<<1|1]));
	return yak_hash64_64(x[(*skip)<<1|0]) + yak_hash64_64(x[(*skip)<<1|1]);
}

inline uint64_t get_hc_pt1_count(ha_ug_index* index, uint64_t key, uint64_t** pos_list)
{
    uint64_t bucket_mask = (1ULL<<index->pre) - 1;
    hc_pt1_t* h = &(index->idx_buf[key & bucket_mask]);
    uint64_t beg;
    khint_t k;
    k = hc_pt_get(h->h, key);
    if (k == kh_end(h->h))
    {
        return 0;
    }
    beg = kh_val(h->h, k);
    if(pos_list) *pos_list = h->a + beg;
    if((kh_key(h->h, k)&HIC_MAX_COUNT)<HIC_MAX_COUNT) return kh_key(h->h, k)&HIC_MAX_COUNT;
    if(k == h->end) return h->n - beg;
    for (k++; k != kh_end(h->h); ++k)
    {
        if (kh_exist(h->h, k)) 
        {
            return kh_val(h->h, k) - beg;
        }
    }
    return h->n - beg;
}

void test_hc_pt1(char* seq, uint64_t len, uint64_t uID, ha_ug_index* idx)
{
    uint64_t i, l, k, pos, *pos_list = NULL, cnt;
    uint64_t x[4], mask = (1ULL<<idx->k) - 1, shift = idx->k - 1, hash, skip;
    for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i) {
        int c = seq_nt4_table[(uint8_t)seq[i]];
        ///c = 00, 01, 10, 11
        if (c < 4) { // not an "N" base
            ///x[0] & x[1] are the forward k-mer
            ///x[2] & x[3] are the reverse complementary k-mer
            x[0] = (x[0] << 1 | (c&1))  & mask;
            x[1] = (x[1] << 1 | (c>>1)) & mask;
            x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
            x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
            if (++l >= idx->k)
            {
                hash = hc_hash_long(x, &skip, idx->k);
                if(skip == (uint64_t)-1) continue;
                pos = (skip << 63) | ((uID << (64-idx->uID_bits))>>1) | (i & idx->pos_mode);
                cnt = get_hc_pt1_count(idx, hash, &pos_list);
                if(cnt == 0) fprintf(stderr, "ERROR cnt, uID: %lu\n", uID);
                for (k = 0; k < cnt; k++)
                {
                    if(pos_list[k]==pos)
                    {
                        pos_list[k] = (uint64_t)-1;
                        break;
                    }
                }
                if(k == cnt) fprintf(stderr, "ERROR k\n");
                
            }
        } else l = 0, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
    }
}

void test_unitig_index(ha_ug_index* idx, ma_ug_t *ug)
{
    double index_time = yak_realtime();
    uint32_t i, j;
    ma_utg_t *u = NULL;
    hc_pt1_t *h = NULL;
    idx->ug = ug;
    for (i = 0; i < idx->ug->u.n; i++)
    {
        u = &(idx->ug->u.a[i]);
        if(u->m == 0) continue;
        test_hc_pt1(u->s, u->len, i, idx);
    }

    for (i = 0; i < idx->tot; i++)
    {
        h = &(idx->idx_buf[i]);
        for (j = 0; j < h->n; j++)
        {
            if(h->a[j] != (uint64_t)-1)
            {
                fprintf(stderr, "ERROR j\n");
            }
        }
        
    }

    fprintf(stderr, "[M::%s::%.3f] ==> Test has been passed\n", __func__, yak_realtime()-index_time);
}

void hc_pt_t_gen_single(hc_pt1_t* pt, uint64_t* up_bound, uint64_t* low_bound)
{
    khint_t k;
    uint64_t c;

    if(up_bound || low_bound)
    {
        for (k = 0; k != kh_end(pt->h); ++k) {
            if (kh_exist(pt->h, k)) {
                if((up_bound && kh_val(pt->h, k) > (*up_bound)) || 
                            (low_bound && kh_val(pt->h, k) < (*low_bound)))
                {
                    kh_val(pt->h, k) = 0;
                    kh_key(pt->h, k) = (kh_key(pt->h, k)&HIC_KEY_MODE)|
                                (kh_val(pt->h, k)<HIC_MAX_COUNT?kh_val(pt->h, k):HIC_MAX_COUNT);
                }
            }
        }
    }
    

    for (k = 0, pt->n = 0; k != kh_end(pt->h); ++k) {
        if (kh_exist(pt->h, k)) {
            c = kh_val(pt->h, k);
            kh_val(pt->h, k) = pt->n;
            pt->n += c;
            pt->end = k;
        }
    }
    CALLOC(pt->a, pt->n);
}

void write_dbug(ma_ug_t* ug, FILE* fp)
{
    ma_utg_t *u = NULL;
    uint32_t t, i;
    fwrite(&(ug->u.n), sizeof(ug->u.n), 1, fp);
    for (i = 0; i < ug->u.n; i++)
    {
        u = &(ug->u.a[i]);
        t = u->len;
        fwrite(&t, sizeof(t), 1, fp);
        t = u->circ;
        fwrite(&t, sizeof(t), 1, fp);
        fwrite(&(u->start), sizeof(u->start), 1, fp);
        fwrite(&(u->end), sizeof(u->end), 1, fp);
        fwrite(&(u->n), sizeof(u->n), 1, fp);
        fwrite(u->a, sizeof(uint64_t), u->n, fp);
    }
}

uint32_t test_dbug(ma_ug_t* ug, FILE* fp)
{
    uint32_t f_flag = 0, t, i, r_flag = 0;
    size_t tt;
    ma_utg_t ua, *ub = NULL; memset(&ua, 0, sizeof(ua));
    f_flag = fread(&tt, sizeof(tt), 1, fp);
    if(f_flag == 0 || tt != ug->u.n) goto DES;
    for (i = 0; i < tt; i++)
    {
        ub = &(ug->u.a[i]);
        f_flag = fread(&t, sizeof(t), 1, fp);
        if(f_flag == 0 || t != ub->len) goto DES;
        f_flag = fread(&t, sizeof(t), 1, fp);
        if(f_flag == 0 || t != ub->circ) goto DES;
        f_flag = fread(&(ua.start), sizeof(ua.start), 1, fp);
        if(f_flag == 0 || ua.start != ub->start) goto DES;
        f_flag = fread(&(ua.end), sizeof(ua.end), 1, fp);
        if(f_flag == 0 || ua.end != ub->end) goto DES;
        f_flag = fread(&(ua.n), sizeof(ua.n), 1, fp);
        if(f_flag == 0 || ua.n != ub->n) goto DES;
        t = ua.n;
        ua.n = 0;
        kv_resize(uint64_t, ua, t);
        ua.n = t;
        f_flag = fread(ua.a, sizeof(uint64_t), ua.n, fp);
        if(f_flag == 0 || memcmp(ua.a, ub->a, ua.n)) goto DES;
    }
    r_flag = 1;

    DES:
    free(ua.a);
    return r_flag;
}

int write_hc_pt_index(ha_ug_index* idx, char* file_name)
{
    char* gfa_name = (char*)malloc(strlen(file_name)+25);
    sprintf(gfa_name, "%s.hic.tlb.bin", file_name);
    FILE* fp = fopen(gfa_name, "w");
    if (!fp) {
        free(gfa_name);
        return 0;
    }
    uint64_t i = HC_PT_MA;

    fwrite(&i, sizeof(i), 1, fp);
    fwrite(&idx->uID_bits, sizeof(idx->uID_bits), 1, fp);
    fwrite(&idx->uID_mode, sizeof(idx->uID_mode), 1, fp);
    fwrite(&idx->pos_bits, sizeof(idx->pos_bits), 1, fp);
    fwrite(&idx->pos_mode, sizeof(idx->pos_mode), 1, fp);
    fwrite(&idx->rev_mode, sizeof(idx->rev_mode), 1, fp);
    fwrite(&idx->k, sizeof(idx->k), 1, fp);
    fwrite(&idx->pre, sizeof(idx->pre), 1, fp);
    fwrite(&idx->tot, sizeof(idx->tot), 1, fp);
    fwrite(&idx->tot_pos, sizeof(idx->tot_pos), 1, fp);

    for (i = 0; i < idx->tot; i++)
    {
        fwrite(&idx->idx_buf[i].n, sizeof(idx->idx_buf[i].n), 1, fp);
        fwrite(&idx->idx_buf[i].end, sizeof(idx->idx_buf[i].end), 1, fp);
        fwrite(idx->idx_buf[i].a, sizeof(uint64_t), idx->idx_buf[i].n, fp);
        hc_pt_save(idx->idx_buf[i].h, fp);
    }

    write_dbug(idx->ug, fp);

    fprintf(stderr, "[M::%s] Index has been written.\n", __func__);
    free(gfa_name);
    fclose(fp);
    return 1;
}

void destory_hc_pt_index(ha_ug_index* idx);
int load_hc_pt_index(ha_ug_index** r_idx, ma_ug_t *ug, char* file_name)
{
    uint64_t flag = 0;
    // double index_time = yak_realtime();
    char* gfa_name = (char*)malloc(strlen(file_name)+25);
    sprintf(gfa_name, "%s.hic.tlb.bin", file_name);
    FILE* fp = fopen(gfa_name, "r");
    if (!fp) {
        free(gfa_name);
        return 0;
    }
    ha_ug_index* idx = NULL; CALLOC(idx, 1);
    uint64_t i;

    flag += fread(&i, sizeof(i), 1, fp);
    if(i != HC_PT_MA)
    {
        free(gfa_name);
        destory_hc_pt_index(idx);
        free(idx);
        (*r_idx) = NULL;
        fclose(fp);
        fprintf(stderr, "[M::%s::] ==> Renew Hi-C index\n", __func__);
        return 0;
    }
    flag += fread(&idx->uID_bits, sizeof(idx->uID_bits), 1, fp);
    flag += fread(&idx->uID_mode, sizeof(idx->uID_mode), 1, fp);
    flag += fread(&idx->pos_bits, sizeof(idx->pos_bits), 1, fp);
    flag += fread(&idx->pos_mode, sizeof(idx->pos_mode), 1, fp);
    flag += fread(&idx->rev_mode, sizeof(idx->rev_mode), 1, fp);
    flag += fread(&idx->k, sizeof(idx->k), 1, fp);
    flag += fread(&idx->pre, sizeof(idx->pre), 1, fp);
    flag += fread(&idx->tot, sizeof(idx->tot), 1, fp);
    flag += fread(&idx->tot_pos, sizeof(idx->tot_pos), 1, fp);
    MALLOC(idx->idx_buf, idx->tot);

    for (i = 0; i < idx->tot; i++)
    {
        flag += fread(&idx->idx_buf[i].n, sizeof(idx->idx_buf[i].n), 1, fp);
        flag += fread(&idx->idx_buf[i].end, sizeof(idx->idx_buf[i].end), 1, fp);
        MALLOC(idx->idx_buf[i].a, idx->idx_buf[i].n);
        flag += fread(idx->idx_buf[i].a, sizeof(uint64_t), idx->idx_buf[i].n, fp);
        hc_pt_load(&(idx->idx_buf[i].h), fp);
    }


    (*r_idx) = idx;

    free(gfa_name);
    if(!test_dbug(ug, fp))
    {
        destory_hc_pt_index(idx);
        free(idx);
        (*r_idx) = NULL;
        fclose(fp);
        fprintf(stderr, "[M::%s::] ==> Renew Hi-C index\n", __func__);
        return 0;
    }

    fclose(fp);
    // fprintf(stderr, "[M::%s::%.3f] ==> HiC index has been loaded\n", __func__, yak_realtime()-index_time);
    return 1;
}

static void worker_for_sort(void *data, long i, int tid) // callback for kt_for()
{
    pldat_t *pl = (pldat_t*)data;
    hc_pt1_t *h = &(pl->h->idx_buf[i]);
    khint_t k;
    uint64_t beg, cnt = 0;
    uint64_t* pos_list;
    for (k = 0; k != kh_end(h->h); ++k) {
        if (kh_exist(h->h, k)) {
            beg = kh_val(h->h, k);
            pos_list = h->a + beg;
            if((kh_key(h->h, k)&HIC_MAX_COUNT)<HIC_MAX_COUNT)
            {
                cnt = kh_key(h->h, k)&HIC_MAX_COUNT;
            }
            else if(k == h->end)
            {
                cnt = h->n - beg;
            }
            else
            {
                for (k++; k != kh_end(h->h); ++k)
                {
                    if (kh_exist(h->h, k)) 
                    {
                        cnt = kh_val(h->h, k) - beg;
                        break;
                    }
                }
            }
            if(cnt > 0) radix_sort_hc_pos(pos_list, pos_list+cnt);
        }
    }

}

void hc_pt_t_gen(ha_ug_index* idx, pldat_t* pl)
{
    if(pl == NULL)
    {
        uint64_t i;
        for (i = 0; i < idx->tot; i++)
        {
            hc_pt_t_gen_single(&(idx->idx_buf[i]), &(idx->up_bound), &(idx->low_bound));
        }
    }
    else
    {
        kt_for(pl->n_thread, worker_for_sort, pl, pl->h->tot);
    }
}

static void worker_for(void *data, long i, int tid) // callback for kt_for()
{
	pldat_t *pl = (pldat_t*)data;
    hc_pt1_t *h = &(pl->h->idx_buf[i]);
    uint64_t m = 0, beg, end, occ;
    khint_t key;
    int absent;
    
    
    if(pl->is_cnt)
    {
        uint64_t* cnt = NULL;
        if(pl->cnt[i].a.n > 2) radix_sort_b64(pl->cnt[i].a.a, pl->cnt[i].a.a + pl->cnt[i].a.n);
        cnt = pl->cnt[i].a.a;
        occ = pl->cnt[i].a.n;
        for (m = beg = end = 0; m < occ; m++)
        {
            if(cnt[beg] == cnt[m])
            {
                end = m;
            }
            else
            {
                key = hc_pt_put(h->h, cnt[beg], &absent);
                if(absent) kh_val(h->h, key) = 0;
                kh_val(h->h, key) += (end - beg + 1);
                kh_key(h->h, key) = (kh_key(h->h, key)&HIC_KEY_MODE)|
                        (kh_val(h->h, key)<HIC_MAX_COUNT?kh_val(h->h, key):HIC_MAX_COUNT);
                beg = end = m;
            }
        }
        if(occ > 0)
        {
            key = hc_pt_put(h->h, cnt[beg], &absent);
            if(absent) kh_val(h->h, key) = 0;
            kh_val(h->h, key) += (end - beg + 1);
            kh_key(h->h, key) = (kh_key(h->h, key)&HIC_KEY_MODE)|
                        (kh_val(h->h, key)<HIC_MAX_COUNT?kh_val(h->h, key):HIC_MAX_COUNT);
        }
        pl->cnt[i].a.n = 0;
    } 
    
    if(!pl->is_cnt) 
    {
        ch_buf_t* pos = NULL;
        uint64_t num, *pos_list = NULL, k, k_n, pos_k;
        if(pl->buf[i].a.n > 2) radix_sort_ch_buf(pl->buf[i].a.a, pl->buf[i].a.a + pl->buf[i].a.n);
        pos = pl->buf[i].a.a;
        occ = pl->buf[i].a.n;
        for (m = beg = end = 0; m < occ; m++)
        {
            if(pos[beg].key == pos[m].key)
            {
                end = m;
            }
            else
            {
                num = get_hc_pt1_count(pl->h, pos[beg].key, &pos_list);
                if(num > 0)
                {
                    k_n=(end-beg+1);pos_k=pos_list[num-1];pos_list[num-1]+=k_n;
                    for (k = 0; k < k_n; k++)
                    {
                        pos_list[pos_k+k] = pos[beg+k].pos;
                    }
                }
                
                beg = end = m;
            }
            
        }
        if(occ > 0)
        {
            num = get_hc_pt1_count(pl->h, pos[beg].key, &pos_list);
            if(num > 0)
            {
                k_n=(end-beg+1);pos_k=pos_list[num-1];pos_list[num-1]+=k_n;
                for (k = 0; k < k_n; k++)
                {
                    pos_list[pos_k+k] = pos[beg+k].pos;
                }
            }
        }
        pl->buf[i].a.n = 0;
    }
}

void parallel_count_hc_pt1(pldat_t* pl)
{
    uint64_t i, l = 0, uID, num_pos = 0, pos_thre;
    uint64_t x[4], mask = (1ULL<<pl->h->k) - 1, shift = pl->h->k - 1, hash, pos, skip, bucket_mask = (1ULL<<pl->h->pre) - 1;
    ma_utg_t *u = NULL;
    ch_buf_t k_pos;
    if(pl->is_cnt) l = ((pl->buf_bytes>>3)/pl->h->tot) + 1, pos_thre = pl->buf_bytes>>3;
    if(!pl->is_cnt) l = ((pl->buf_bytes>>4)/pl->h->tot) + 1, pos_thre = pl->buf_bytes>>4;
    for (i = 0; i < pl->h->tot; i++)
    {
        if(pl->is_cnt)
        {
            kv_resize(uint64_t, pl->cnt[i].a, l);
            pl->cnt[i].a.n = 0;
        } 
        
        if(!pl->is_cnt)
        {
            kv_resize(ch_buf_t, pl->buf[i].a, l);
            pl->buf[i].a.n = 0;
        } 
    }
        
     
    for (uID = 0; uID < pl->h->ug->u.n; uID++)
    {
        u = &(pl->h->ug->u.a[uID]);
        if(u->m == 0) continue;

        for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < u->len; ++i) {
        int c = seq_nt4_table[(uint8_t)u->s[i]];
        ///c = 00, 01, 10, 11
        if (c < 4) { // not an "N" base
            ///x[0] & x[1] are the forward k-mer
            ///x[2] & x[3] are the reverse complementary k-mer
            x[0] = (x[0] << 1 | (c&1))  & mask;
            x[1] = (x[1] << 1 | (c>>1)) & mask;
            x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
            x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
            if (++l >= pl->h->k)
            {
                hash = hc_hash_long(x, &skip, pl->h->k);
                if(skip == (uint64_t)-1) continue;
                if(pl->is_cnt)
                {
                    kv_push(uint64_t, pl->cnt[hash & bucket_mask].a, hash);
                }
                else
                {
                    pos = (skip << 63) | ((uID << (64-pl->h->uID_bits))>>1) | (i & pl->h->pos_mode);
                    k_pos.key = hash; k_pos.pos = pos;
                    kv_push(ch_buf_t, pl->buf[hash & bucket_mask].a, k_pos);
                }
                num_pos++;

                if(num_pos >= pos_thre)
                {
                    num_pos = 0;
                    kt_for(pl->n_thread, worker_for, pl, pl->h->tot);
                }
            }
                
        } else l = 0, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
    }
    }

    if(num_pos > 0) kt_for(pl->n_thread, worker_for, pl, pl->h->tot);

    for (i = 0; i < pl->h->tot; i++)
    {
        if(pl->cnt[i].a.m > 0) kv_destroy(pl->cnt[i].a), kv_init(pl->cnt[i].a);
        if(pl->buf[i].a.m > 0) kv_destroy(pl->buf[i].a), kv_init(pl->buf[i].a);
    }
}

ha_ug_index* build_unitig_index(ma_ug_t *ug, int k, uint64_t up_occ, uint64_t low_occ, uint64_t thread_num)
{
    ha_ug_index* idx = NULL; CALLOC(idx, 1);
    pldat_t pl; pl.h = idx; pl.is_cnt = 1;
    double index_time = yak_realtime(), beg_time;
    init_ha_ug_index_opt(idx, ug, k, &pl, up_occ, low_occ, thread_num);

    beg_time = yak_realtime();
    pl.is_cnt = 1;
    parallel_count_hc_pt1(&pl);
    fprintf(stderr, "[M::%s::%.3f] ==> Counting\n", __func__, yak_realtime()-beg_time);

    beg_time = yak_realtime();
    hc_pt_t_gen(pl.h, NULL);
    fprintf(stderr, "[M::%s::%.3f] ==> Memory allocating\n", __func__, yak_realtime()-beg_time);

    beg_time = yak_realtime();
    pl.is_cnt = 0;
    parallel_count_hc_pt1(&pl);
    fprintf(stderr, "[M::%s::%.3f] ==> Filling pos\n", __func__, yak_realtime()-beg_time);

    beg_time = yak_realtime();
    hc_pt_t_gen(pl.h, &pl);
    fprintf(stderr, "[M::%s::%.3f] ==> Sorting pos\n", __func__, yak_realtime()-beg_time);

    fprintf(stderr, "[M::%s::%.3f] ==> HiC index has been built\n", __func__, yak_realtime()-index_time);
    
    uint64_t i;
    for (i = 0; i < idx->tot; i++)
    {
        kv_destroy(pl.cnt[i].a);
        kv_destroy(pl.buf[i].a);
    }
    free(pl.cnt); free(pl.buf);

    return idx;
}

void destory_hc_pt_index(ha_ug_index* idx)
{
    if(idx->idx_buf)
    {
        uint64_t i = 0;
        for (i = 0; i < idx->tot; i++)
        {
            if(idx->idx_buf[i].a) free(idx->idx_buf[i].a);
            if(idx->idx_buf[i].h) hc_pt_destroy(idx->idx_buf[i].h);
        }
        free(idx->idx_buf);
    }
}

inline void interpret_pos(const ha_ug_index* idx, s_hit *p, uint64_t* rev, uint64_t* uID, 
uint64_t* ref_p, uint64_t* self_p, uint64_t* exact_len, uint64_t* total_len)
{
    (*rev) = p->ref>>63;
    (*uID) = (p->ref << 1) >> (64 - idx->uID_bits);
    (*self_p) = (uint32_t)p->off_cnt;
    ///(*exact_len) = p->off_cnt >> 32;
    (*exact_len) = (p->off_cnt>>32) & ((uint64_t)65535);
    if(total_len != NULL)
    {
        ///(*exact_len) = (p->off_cnt>>32) & ((uint64_t)65535);
        (*total_len) = (p->off_cnt>>48) + (*exact_len);
    }
    if((p->ref & idx->pos_mode)>>(idx->pos_bits - 1))
    {
        (*ref_p) = (*self_p) - (p->ref&(idx->pos_mode>>1));
    }
    else
    {
        (*ref_p) = (*self_p) + (p->ref&(idx->pos_mode));
    }
}


inline uint64_t check_exact_match(char* a, long long a_beg, long long a_total, char* b, long long b_beg, 
long long b_total, long long Len, uint64_t rev, uint64_t dir)
{
    long long i = 0;
    if(rev == 0)
    {
        if(dir == 0)
        {
            for (i = 0; i < Len && a_beg < a_total && b_beg < b_total; i++)
            {
                if(a[a_beg++] != b[b_beg++]) return i;
            }
        }
        else
        {
            for (i = 0; i < Len && a_beg >= 0 && b_beg >= 0; i++)
            {
                if(a[a_beg--] != b[b_beg--]) return i;
            }
        }
    }
    else
    {

        if(dir == 0)
        {
            for (i = 0; i < Len && a_beg < a_total && b_beg < b_total; i++)
            {
                if(a[a_beg] != b2rc[seq_nt4_table[(uint8_t)b[b_total - b_beg - 1]]]) return i;
                a_beg++; b_beg++;
            }
        }
        else
        {
            for (i = 0; i < Len && a_beg >= 0 && b_beg >= 0; i++)
            {
                if(a[a_beg] != b2rc[seq_nt4_table[(uint8_t)b[b_total - b_beg - 1]]]) return i;
                a_beg--; b_beg--;
            }
        }
    }

    return i;
}

uint64_t debug_hash_value(char *r, uint64_t end, uint64_t k_mer)
{
    uint64_t i;
    uint64_t x[4], mask = (1ULL<<k_mer) - 1, shift = k_mer - 1, skip;
    for (i = end + 1 - k_mer, x[0] = x[1] = x[2] = x[3] = 0; i <= end; i++)
    {
        int c = seq_nt4_table[(uint8_t)r[i]];
        ///c = 00, 01, 10, 11
        if (c < 4) { // not an "N" base
            ///x[0] & x[1] are the forward k-mer
            ///x[2] & x[3] are the reverse complementary k-mer
            x[0] = (x[0] << 1 | (c&1))  & mask;
            x[1] = (x[1] << 1 | (c>>1)) & mask;
            x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
            x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
        }
    }

    return hc_hash_long(x, &skip, k_mer);
    
}

inline uint64_t collect_votes(s_hit* a, uint64_t n)
{
    if(n == 0) return 0;
    if(n == 1) return (a[0].off_cnt>>32); //seed length, is right
    long long i = 0;
    uint64_t cur_beg, cur_end, beg, end, ovlp = 0, tLen = 0; 
    cur_end = (uint32_t)a[n-1].off_cnt;
    cur_beg = cur_end + 1 - (a[n-1].off_cnt>>32);
    
    if(n >= 2)
    {
        for (i = n - 2; i >= 0; i--)
        {
            end = (uint32_t)a[i].off_cnt;
            beg = end + 1 - (a[i].off_cnt>>32);
            if(MAX(cur_beg, beg) <= MIN(cur_end, end))
            {
                cur_beg = MIN(cur_beg, beg);
                ///cur_end = MAX(cur_end, end);
            }
            else
            {
                ovlp += (cur_end + 1 - cur_beg);
                cur_beg = beg;
                cur_end = end;
            }
        }
    }
    
    ovlp += (cur_end + 1 - cur_beg);
    tLen = (uint32_t)a[n-1].off_cnt + 1 - cur_beg;
    tLen = tLen - ovlp;
    tLen = tLen << 16;
    return ovlp | tLen;
}

inline void compress_mapped_pos(const ha_ug_index* idx, kvec_vote* buf, uint64_t buf_iter, uint64_t max_i, uint64_t thres)
{
    if(buf_iter >= buf->a.n)
    {
        buf->a.n = buf_iter;
        return;
    }
    uint64_t rev, uID, ref_p, self_p, eLen, tLen, i, max_beg, max_end, cur_beg, cur_end, ovlp, max_eLen;
    uint64_t secondLen = 0, second_i = (uint64_t)-1;
    interpret_pos((ha_ug_index*)idx, &buf->a.a[max_i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
    max_end = self_p;
    max_beg = self_p + 1 - tLen;
    max_eLen = eLen;
    for (i = buf_iter; i < buf->a.n; i++)
    {
        if(i == max_i) continue;
        interpret_pos(idx, &buf->a.a[i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
        cur_end = self_p;
        cur_beg = self_p + 1 - tLen;
        if(MAX(cur_beg, max_beg) <= MIN(cur_end, max_end))
        {
            ovlp = MIN(cur_end, max_end) - MAX(cur_beg, max_beg) + 1;
            if(ovlp > thres)
            {
                if(eLen >= max_eLen * 0.8)
                {
                    buf->a.n = buf_iter;
                    return;
                }
                continue;
            } 
        }
        if(secondLen < eLen) secondLen = eLen, second_i = i;
    }

    if(second_i == (uint64_t)-1)
    {
        buf->a.a[buf_iter] = buf->a.a[max_i];
        buf->a.n = buf_iter + 1;
    }
    else
    {
        buf->a.a[buf_iter] = buf->a.a[MIN(max_i, second_i)];
        buf->a.a[buf_iter+1] = buf->a.a[MAX(max_i, second_i)];
        buf->a.n = buf_iter + 2;
    }
}

inline void print_pos_list(const ha_ug_index* idx, s_hit *l, uint64_t occ, uint64_t rid, uint64_t r1)
{
    if(rid == 33045391 || rid == 4239289 || rid == 5267597 || rid == 34474764 || rid == 35016489
            || rid == 36002255 || rid == 37811694 || rid == 46805824)
    {
        uint64_t i, rev, uID, ref_p, self_p, cnt;
        for (i = 0; i < occ; i++)
        {
            interpret_pos(idx, &l[i], &rev, &uID, &ref_p, &self_p, &cnt, NULL);
            fprintf(stderr, "(r%lu) rid: %lu, i: %lu, rev: %lu, uID: %lu, ref_p: %lu, self_p: %lu\n", 
                                        r1, rid, i, rev, uID, ref_p, self_p);
        }
        
    }
}

uint64_t get_longest_hit(char *r, uint64_t len, uint64_t k_mer, uint64_t self_p, uint64_t self_rev, kvec_vote* buf, const ha_ug_index* idx, 
uint64_t *pos_list, uint64_t cnt, uint64_t* c_sfx)
{
    uint64_t max_p, map_p_occ, i, j, m, rev, ref_p, u_len, uID, k_len;
    s_hit *p = NULL;
    ///each k-mer at different unitigs
    ///rev:uID:pos
    if(c_sfx) (*c_sfx) = (uint64_t)-1;
    for (j = 0; j < cnt; j++)
    {
        ///get
        kv_pushp(s_hit, buf->a, &p);
        rev = (pos_list[j]>>63) != self_rev;
        ref_p = pos_list[j] & idx->pos_mode;
        uID = (pos_list[j] << 1) >> (64 - idx->uID_bits);
        u_len = idx->ug->u.a[uID].len;
        if(rev) ref_p = u_len - 1 - (ref_p + 1 - k_mer);
        p->off_cnt = self_p | ((uint64_t)k_mer << 32); ///high bits should be the legnth

        p->ref = ref_p >= self_p? (ref_p-self_p) 
                        : (self_p-ref_p) + ((uint64_t)1 << (idx->pos_bits - 1));
        p->ref = (rev << 63)|(pos_list[j] & idx->uID_mode)|(p->ref&idx->pos_mode);


        ///extend
        k_len = check_exact_match(r, self_p + 1, len, idx->ug->u.a[uID].s, ref_p + 1, u_len, len, rev, 0);
        if(c_sfx && cnt == idx->hap_cnt && k_len < (*c_sfx)) (*c_sfx) = k_len;
        
        p->off_cnt += ((uint64_t)k_len << 32) + k_len;

        if(self_p >= k_mer && ref_p >= k_mer)
        {
            k_len = check_exact_match(r, self_p - k_mer, len, idx->ug->u.a[uID].s, 
                                                    ref_p - k_mer, u_len, len, rev, 1);
            p->off_cnt += ((uint64_t)k_len << 32);
        }
        // if(cnt > 0) fprintf(stderr, "inner j: %lu, rev: %lu, uID: %lu, ref_p: %lu, self_p: %u, len: %lu\n", j, rev, uID, ref_p, (uint32_t)p->off_cnt, p->off_cnt>>32);
    }

    p = buf->a.a + buf->a.n - cnt;
    if(cnt > 1) radix_sort_hc_s_hit_off_cnt(p, p + cnt);
    max_p = map_p_occ = 0; 
    for (j = 1, i = 0; j <= cnt; ++j)
    {
        if(j == cnt || p[j].off_cnt != p[i].off_cnt)
        {
            if((max_p>>32) < (p[i].off_cnt>>32))
            {
                max_p = p[i].off_cnt;
                map_p_occ = j - i;
            }
            else if(((max_p>>32) == (p[i].off_cnt>>32)) && ((j - i) > map_p_occ))
            {
                max_p = p[i].off_cnt;
                map_p_occ = j - i;
            }
            i = j;///must
        }
    }

    buf->a.n -= cnt;
    for (j = m = 0; j < cnt; j++)
    {
        if(p[j].off_cnt == max_p)
        {
            p[m] = p[j];
            m++;
        }
    }
    cnt = m;
    buf->a.n += cnt;
    
    // if(cnt > 0) fprintf(stderr, "max_p_offset: %u, max_p_len: %lu, map_p_occ: %lu\n", (uint32_t)max_p, max_p>>32, map_p_occ);
    return max_p;
}

#define is_update_hit(mL, mR, cL, cR) (((mL)<(cL))||((mL)==(cL)&&(mR)<(cR))) 
inline void compress_mapped_pos_advance(const ha_ug_index* idx, kvec_vote* buf, uint64_t buf_iter, uint64_t ovlp_thre)
{
    if(buf_iter >= buf->a.n)
    {
        buf->a.n = buf_iter;
        return;
    }
    s_hit *p = NULL;
    uint64_t rev, uID, ref_p, self_p, eLen, tLen, i, j, cnt;
    uint64_t max_beg = 0, max_end = 0, max_i, max_occ, cur_beg, cur_end, ovlp;
    uint64_t second_i = (uint64_t)-1, second_occ;
    uint64_t max_eLen, sec_eLen;
    double max_eRate, sec_eRate, eRate;
    p = buf->a.a + buf_iter;
    cnt = buf->a.n - buf_iter;
    if(cnt > 1) radix_sort_hc_s_hit_off_cnt(p, p + cnt); ///buf save all hits, here sort by offset in reads 
    max_eLen = 0; max_i = (uint64_t)-1; max_occ = 0; max_eRate = -1;
    for (j = 1, i = 0; j <= cnt; ++j)
    {
        if(j == cnt || p[j].off_cnt != p[i].off_cnt)
        {
            ///occ = j - i;
            interpret_pos((ha_ug_index*)idx, &p[i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
            eRate = (double)(eLen)/(double)(tLen);
            if(is_update_hit(max_eLen, max_eRate, eLen, eRate))
            {
                max_eLen = eLen; max_eRate = eRate;
                max_end = self_p; max_beg = self_p + 1 - tLen;
                max_i = i; max_occ = j - i;
            }
            // fprintf(stderr, "\n++++++[%lu, %lu] uID: %lu, ref_p: %lu, self_p: %lu, eLen: %lu, tLen: %lu, max_i: %lu\n", 
            // i, j, uID, ref_p, self_p, eLen, tLen, max_i);
            i = j;///must
        }
    }


    sec_eLen = 0; second_i = (uint64_t)-1; second_occ = 0; sec_eRate = -1;
    for (j = 1, i = 0; j <= cnt; ++j)
    {
        if(j == cnt || p[j].off_cnt != p[i].off_cnt)
        {
            if(i != max_i)
            {
                interpret_pos((ha_ug_index*)idx, &p[i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
                eRate = (double)(eLen)/(double)(tLen);
                cur_end = self_p;
                cur_beg = self_p + 1 - tLen;
                // fprintf(stderr, "\n----[%lu, %lu] uID: %lu, ref_p: %lu, self_p: %lu, eLen: %lu, tLen: %lu, max_i: %lu\n", 
                // i, j, uID, ref_p, self_p, eLen, tLen, max_i);
                // fprintf(stderr, "max_beg: %lu, max_end: %lu, cur_beg: %lu, cur_end: %lu\n", 
                // max_beg, max_end, cur_beg, cur_end);
                ///overlap with max interval
                if(MAX(cur_beg, max_beg) <= MIN(cur_end, max_end))
                {
                    ovlp = MIN(cur_end, max_end) - MAX(cur_beg, max_beg) + 1;
                    /*******************************for debug************************************/
                    if(ovlp == MIN(max_end+1-max_end, tLen))///for non-unique k-mer
                    {
                        i = j;///must
                        continue;///fully contain
                    }
                    if(ovlp > ((max_end+1-max_end)*0.8) && eLen > (max_eLen*0.8))///best is not unique
                    {
                        buf->a.n = buf_iter;
                        return;
                    }
                    if(ovlp > ((max_end+1-max_end)*0.15) + 1)
                    {
                        i = j;///must
                        continue;///fully contain
                    }
                    
                    // if(ovlp > (MIN((max_end+1-max_end), (cur_end+1-cur_end))*0.15) + 1)
                    // {
                    //     if(eLen > (max_eLen*0.8))///best is not unique
                    //     {
                    //         buf->a.n = buf_iter;
                    //         return;
                    //     }
                    //     i = j;///must
                    //     continue;
                    // }
                    /*******************************for debug************************************/
                }

                if(is_update_hit(sec_eLen, sec_eRate, eLen, eRate))
                {
                    sec_eLen = eLen; sec_eRate = eRate;
                    second_i = i; second_occ = j - i;
                }
            }
            i = j;///must
        }
    }

    // fprintf(stderr, "max_i: %lu, max_occ: %lu, second_i: %lu, second_occ: %lu\n", 
    //             max_i, max_occ, second_i, second_occ);

    if(second_i == (uint64_t)-1)
    {
        i = 0;
        for (j = max_i; j < max_i + max_occ; j++, i++) p[i] = p[j];
    }
    else ///be carful about overwritten
    {
        i = 0;
        if(max_i <= second_i)
        {
            for (j = max_i; j < max_i + max_occ; j++, i++) p[i] = p[j];
            for (j = second_i; j < second_i + second_occ; j++, i++) p[i] = p[j];
        }
        else
        {
            for (j = second_i; j < second_i + second_occ; j++, i++) p[i] = p[j];
            for (j = max_i; j < max_i + max_occ; j++, i++) p[i] = p[j];
        } 
    }

    buf->a.n = buf_iter + max_occ + second_occ;
}

void get_alignment(char *r, uint64_t len, uint64_t k_mer, kvec_vote* buf, const ha_ug_index* idx, uint64_t buf_iter, uint64_t rid)
{
    uint64_t i, j, k, l = 0, k_len, c_sfx, m, skip, *pos_list = NULL, cnt, rev, self_p, ref_p, uID;
    uint64_t x[4], mask = (1ULL<<k_mer) - 1, shift = k_mer - 1, hash;
    ///buf->a.n = 0;
    for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i) {
        int c = seq_nt4_table[(uint8_t)r[i]];
        ///c = 00, 01, 10, 11
        if (c < 4) { // not an "N" base
            ///x[0] & x[1] are the forward k-mer
            ///x[2] & x[3] are the reverse complementary k-mer
            x[0] = (x[0] << 1 | (c&1))  & mask;
            x[1] = (x[1] << 1 | (c>>1)) & mask;
            x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
            x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
            if (++l >= k_mer)
            {
                hash = hc_hash_long(x, &skip, k_mer);
                if(skip == (uint64_t)-1) continue;
                cnt = get_hc_pt1_count((ha_ug_index*)idx, hash, &pos_list);
                if(cnt > idx->hap_cnt || cnt <= 0) continue;

                if(cnt > 1)
                {
                    for (j = 0; j < cnt; j++)
                    {
                        uID = (pos_list[j] << 1) >> (64 - idx->uID_bits);
                        for (k = j + 1; k < cnt; k++)
                        {
                            if(uID == ((pos_list[k] << 1) >> (64 - idx->uID_bits))) break;
                        }
                        if(k < cnt) break;
                    }

                    if(j < cnt) continue;
                }
                // if(cnt > 0) fprintf(stderr, "+i: %lu, l: %lu, cnt: %lu\n", i, l, cnt);
                get_longest_hit(r, len, k_mer, i, skip, buf, idx, pos_list, cnt, &c_sfx);
                // if(cnt > 0) fprintf(stderr, "c_sfx: %lu\n", c_sfx);
                if(c_sfx != (uint64_t)-1)
                {
                    k_len = c_sfx;
                    if((k_len + 1) >= k_mer)
                    {
                        l = 0, x[0] = x[1] = x[2] = x[3] = 0;
                        i = i + k_len - (k_mer - 1);
                    }
                    else
                    {
                        ///l = i - (i + k_len - (k_mer - 1));
                        l = k_mer - k_len - 1;
                    }  
                }
                // if(cnt > 0) fprintf(stderr, "-i: %lu, l: %lu\n", i, l);
            }
            
        } else l = 0, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
    }
    
    if(buf->a.n - buf_iter == 0) return;
    if(buf->a.n - buf_iter > 1) radix_sort_hc_s_hit_an1(buf->a.a + buf_iter, buf->a.a + buf->a.n);


    uint64_t cur_ref_p, thres = (len * HIC_R_E_RATE) + 1, index_beg, ovlp;
    i = m = buf_iter;
    while (i < buf->a.n)
    {
        interpret_pos(idx, &buf->a.a[i], &rev, &uID, &ref_p, &self_p, &cnt, NULL);
        ///fprintf(stderr, "after-i: %lu, uID: %lu, ref_p: %lu, self_p: %lu\n", i, uID, ref_p, self_p);
        cur_ref_p = buf->a.a[i].ref;
        index_beg = i;
        ///ref>>(idx->pos_bits-1) = (rev:1):(uID:uID-bits):(ref_pos>=self_pos:1)
        while ((i < buf->a.n) && 
                ((buf->a.a[i].ref>>(idx->pos_bits-1)) == (cur_ref_p>>(idx->pos_bits-1))) && 
               (buf->a.a[i].ref - cur_ref_p <= thres))
        {
            i++;
        }
        if(i - index_beg > 1)
        {
            radix_sort_hc_s_hit_an2(buf->a.a + index_beg, buf->a.a + i);//sort by self_p
        }
        ovlp = collect_votes(buf->a.a + index_beg, i - index_beg);
        ///fprintf(stderr, "i-1: %lu, self_p: %u\n", i-1, (uint32_t)buf->a.a[i - 1].off_cnt);
        buf->a.a[m] = buf->a.a[i - 1];
        buf->a.a[m].off_cnt = (buf->a.a[m].off_cnt << 32)>>32;
        buf->a.a[m].off_cnt += ((uint64_t)ovlp<<32);
        ///fprintf(stderr, "m: %lu, self_p: %u\n", m, (uint32_t)buf->a.a[m].off_cnt);
        m++;
    }
    buf->a.n = m; 

    /*******************************for debug************************************/
    // for (i = buf_iter; i < buf->a.n; i++)
    // {
    //     uint64_t eLen, tLen;
    //     interpret_pos(idx, &buf->a.a[i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
    //     if(maxLen < eLen) fprintf(stderr, "ERROR1\n");
    //     if(i == max_i && maxLen != eLen) fprintf(stderr, "ERROR2\n");
    // }
    /*******************************for debug************************************/
    ///select the best alignment at [buf_iter, m)
    
    /*******************************for debug************************************/
    // fprintf(stderr, "len1:%lu, max_i: %lu\n", buf->a.n - buf_iter, max_i);
    // for (i = buf_iter; i < buf->a.n; i++)
    // {
    //     uint64_t eLen, tLen;
    //     interpret_pos(idx, &buf->a.a[i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
    //     fprintf(stderr, "(%lu) rev: %lu, uID: %lu, ref_p: %lu, self_p: %lu, eLen: %lu, tLen: %lu\n", 
    //                                                         i, rev, uID, ref_p, self_p, eLen, tLen);
    // }
    /*******************************for debug************************************/
    compress_mapped_pos_advance(idx, buf, buf_iter, (k_mer * 0.1) > 0? (k_mer * 0.1) : 1);
    /*******************************for debug************************************/
    // fprintf(stderr, "len2:%lu, max_i: %lu\n", buf->a.n - buf_iter, max_i);
    // for (i = buf_iter; i < buf->a.n; i++)
    // {
    //     uint64_t eLen, tLen;
    //     interpret_pos(idx, &buf->a.a[i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
    //     fprintf(stderr, "(%lu) rev: %lu, uID: %lu, ref_p: %lu, self_p: %lu, eLen: %lu, tLen: %lu\n", 
    //                                                         i, rev, uID, ref_p, self_p, eLen, tLen);
    // }
    // if(buf->a.n != m) fprintf(stderr, "Changed\n");
    // fprintf(stderr, "\n");
    /*******************************for debug************************************/
}

inline void compress_mapped_pos_debug(const ha_ug_index* idx, kvec_vote* buf, uint64_t buf_iter, uint64_t max_i, uint64_t thres)
{
    if(buf_iter >= buf->a.n)
    {
        buf->a.n = buf_iter;
        return;
    }
    uint64_t rev, uID, ref_p, self_p, eLen, tLen, i, max_beg, max_end, cur_beg, cur_end, ovlp, max_eLen;
    uint64_t secondLen = 0, second_i = (uint64_t)-1;
    interpret_pos((ha_ug_index*)idx, &buf->a.a[max_i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
    max_end = self_p;
    max_beg = self_p + 1 - tLen;
    max_eLen = eLen;
    for (i = buf_iter; i < buf->a.n; i++)
    {
        if(i == max_i) continue;
        interpret_pos(idx, &buf->a.a[i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
        cur_end = self_p;
        cur_beg = self_p + 1 - tLen;
        if(MAX(cur_beg, max_beg) <= MIN(cur_end, max_end))
        {
            ovlp = MIN(cur_end, max_end) - MAX(cur_beg, max_beg) + 1;
            if(ovlp > thres)
            {
                if(eLen >= max_eLen * 0.8)
                {
                    buf->a.n = buf_iter;
                    return;
                }
                continue;
            } 
        }
        if(secondLen < eLen) secondLen = eLen, second_i = i;
    }

    if(second_i == (uint64_t)-1)
    {
        buf->a.a[buf_iter] = buf->a.a[max_i];
        buf->a.n = buf_iter + 1;
    }
    else
    {
        buf->a.a[buf_iter] = buf->a.a[MIN(max_i, second_i)];
        buf->a.a[buf_iter+1] = buf->a.a[MAX(max_i, second_i)];
        buf->a.n = buf_iter + 2;
    }
}

void get_alignment_debug(char *r, uint64_t len, uint64_t k_mer, kvec_vote* buf, const ha_ug_index* idx, uint64_t buf_iter, uint64_t rid)
{
    uint64_t i, j, l = 0,  m, skip, *pos_list = NULL, cnt, rev, self_p, ref_p, uID;
    uint64_t x[4], mask = (1ULL<<k_mer) - 1, shift = k_mer - 1, hash;
    /****************************may have bugs********************************/
    /**
    ///buf->a.n = 0;
    uint64_t k_len, c_sfx, k;
    for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i) {
        int c = seq_nt4_table[(uint8_t)r[i]];
        ///c = 00, 01, 10, 11
        if (c < 4) { // not an "N" base
            ///x[0] & x[1] are the forward k-mer
            ///x[2] & x[3] are the reverse complementary k-mer
            x[0] = (x[0] << 1 | (c&1))  & mask;
            x[1] = (x[1] << 1 | (c>>1)) & mask;
            x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
            x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
            if (++l >= k_mer)
            {
                hash = hc_hash_long(x, &skip, k_mer);
                if(skip == (uint64_t)-1) continue;
                cnt = get_hc_pt1_count((ha_ug_index*)idx, hash, &pos_list);
                if(cnt > idx->hap_cnt || cnt <= 0) continue;

                if(cnt > 1)
                {
                    for (j = 0; j < cnt; j++)
                    {
                        uID = (pos_list[j] << 1) >> (64 - idx->uID_bits);
                        for (k = j + 1; k < cnt; k++)
                        {
                            if(uID == ((pos_list[k] << 1) >> (64 - idx->uID_bits))) break;
                        }
                        if(k < cnt) break;
                    }

                    if(j < cnt) continue;
                }
                // if(cnt > 0) fprintf(stderr, "+i: %lu, l: %lu, cnt: %lu\n", i, l, cnt);
                get_longest_hit(r, len, k_mer, i, skip, buf, idx, pos_list, cnt, &c_sfx);
                // if(cnt > 0) fprintf(stderr, "c_sfx: %lu\n", c_sfx);
                if(c_sfx != (uint64_t)-1)
                {
                    k_len = c_sfx;
                    if((k_len + 1) >= k_mer)
                    {
                        l = 0, x[0] = x[1] = x[2] = x[3] = 0;
                        i = i + k_len - (k_mer - 1);
                    }
                    else
                    {
                        ///l = i - (i + k_len - (k_mer - 1));
                        l = k_mer - k_len - 1;
                    }  
                }
                // if(cnt > 0) fprintf(stderr, "-i: %lu, l: %lu\n", i, l);
            }
            
        } else l = 0, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
    }
    **/

    s_hit *p = NULL; uint64_t u_len;
    ///buf->a.n = 0;
    for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i) {
        int c = seq_nt4_table[(uint8_t)r[i]];
        ///c = 00, 01, 10, 11
        if (c < 4) { // not an "N" base
            ///x[0] & x[1] are the forward k-mer
            ///x[2] & x[3] are the reverse complementary k-mer
            x[0] = (x[0] << 1 | (c&1))  & mask;
            x[1] = (x[1] << 1 | (c>>1)) & mask;
            x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
            x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
            if (++l >= k_mer)
            {
                hash = hc_hash_long(x, &skip, k_mer);
                if(skip == (uint64_t)-1) continue;
                /*******************************for debug************************************/
                // if(debug_hash_value(r, i, k_mer) != hash)
                // {
                //     fprintf(stderr, "ERROR\n");
                // }
                /*******************************for debug************************************/
                cnt = get_hc_pt1_count((ha_ug_index*)idx, hash, &pos_list);
                if(cnt > idx->hap_cnt) continue;
                if(cnt != 1) continue; ///might be able to be disabled in future

                
                for (j = 0; j < cnt; j++)
                {
                    kv_pushp(s_hit, buf->a, &p);
                    rev = (pos_list[j]>>63) != skip;
                    self_p = i;
                    ref_p = pos_list[j] & idx->pos_mode;
                    uID = (pos_list[j] << 1) >> (64 - idx->uID_bits);
                    u_len = idx->ug->u.a[uID].len;
                    if(rev) ref_p = u_len - 1 - (ref_p + 1 - k_mer);
                    p->off_cnt = self_p | ((uint64_t)k_mer << 32); ///high bits should be the legnth

                    p->ref = ref_p >= self_p? (ref_p-self_p) 
                                    : (self_p-ref_p) + ((uint64_t)1 << (idx->pos_bits - 1));
                    p->ref = (rev << 63)|(pos_list[j] & idx->uID_mode)|(p->ref&idx->pos_mode);


                    /*******************************for debug************************************/
                    // if(check_exact_match(r, i + 1 - k_mer, len,
                    //                     idx->ug->u.a[uID].s, ref_p  + 1 - k_mer, u_len, k_mer, rev, 0) != k_mer
                    //     ||
                    //    check_exact_match(r, i, len,
                    //                     idx->ug->u.a[uID].s, ref_p, u_len, k_mer, rev, 1) != k_mer)
                    // {
                    //     fprintf(stderr, "ERROR\n");
                    // }
                    /*******************************for debug************************************/
                }
                
                if(cnt == 1)
                {
                    ///uint64_t debug_right = 0, debug_left = 0, debug_len;

                    j = check_exact_match(r, self_p + 1, len, idx->ug->u.a[uID].s, ref_p + 1, u_len, len, rev, 0);
                    
                    ///debug_right = j;
                    ///if(j == 0) continue;
                    if((j + 1) >= k_mer)
                    {
                        l = 0, x[0] = x[1] = x[2] = x[3] = 0;
                        i = i + j - (k_mer - 1);
                    }
                    else
                    {
                        ///l = i - (i + j - (k_mer - 1));
                        l = k_mer - j -1;
                    }
                    buf->a.a[buf->a.n-1].off_cnt += ((uint64_t)j << 32) + j;

                    if(self_p >= k_mer && ref_p >= k_mer)
                    {
                        j = check_exact_match(r, self_p - k_mer, len, idx->ug->u.a[uID].s, 
                                                                ref_p - k_mer, u_len, len, rev, 1);
                        buf->a.a[buf->a.n-1].off_cnt += ((uint64_t)j << 32);
                        ///debug_left = j;
                    }


                    // debug_len = check_exact_match(r, self_p + debug_right, len, idx->ug->u.a[uID].s, 
                    // ref_p + debug_right, u_len, len, rev, 1);
                    // if(debug_len!= (debug_left + debug_right + k_mer))
                    // {
                    //     fprintf(stderr, "debug_len: %lu, debug_left: %lu, debug_right: %lu\n",
                    //     debug_len, debug_left, debug_right);
                    // }
                }
                
            }
            
        } else l = 0, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
    }
    /****************************may have bugs********************************/
    
    if(buf->a.n - buf_iter == 0) return;
    if(buf->a.n - buf_iter > 1) radix_sort_hc_s_hit_an1(buf->a.a + buf_iter, buf->a.a + buf->a.n);


    uint64_t cur_ref_p, thres = (len * HIC_R_E_RATE) + 1, index_beg, ovlp;
    /****************************may have bugs********************************/
    uint64_t maxLen = 0, max_i = (uint64_t)-1;
    /****************************may have bugs********************************/
    i = m = buf_iter;
    while (i < buf->a.n)
    {
        interpret_pos(idx, &buf->a.a[i], &rev, &uID, &ref_p, &self_p, &cnt, NULL);
        cur_ref_p = buf->a.a[i].ref;
        index_beg = i;
        while ((i < buf->a.n) && 
                ((buf->a.a[i].ref>>(idx->pos_bits-1)) == (cur_ref_p>>(idx->pos_bits-1))) && 
               (buf->a.a[i].ref - cur_ref_p <= thres))
        {
            i++;
        }
        if(i - index_beg > 1)
        {
            radix_sort_hc_s_hit_an2(buf->a.a + index_beg, buf->a.a + i);//sort by self_p
        }
        ovlp = collect_votes(buf->a.a + index_beg, i - index_beg);
        buf->a.a[m] = buf->a.a[i - 1];
        buf->a.a[m].off_cnt = (buf->a.a[m].off_cnt << 32)>>32;
        buf->a.a[m].off_cnt += ((uint64_t)ovlp<<32);
        m++;
        /****************************may have bugs********************************/
        if(maxLen < (ovlp&((uint64_t)65535))) maxLen = (ovlp&((uint64_t)65535)), max_i = m;
        /****************************may have bugs********************************/
    }
    buf->a.n = m; 
    /****************************may have bugs********************************/
    ///compress_mapped_pos_advance(idx, buf, buf_iter, (k_mer * 0.1) > 0? (k_mer * 0.1) : 1);
    compress_mapped_pos_debug(idx, buf, buf_iter, max_i, thres);
    /****************************may have bugs********************************/
}

inline int is_unreliable_hits(long long rev, long long ref_p, long long tLen, uint64_t uID, trans_chain* t_ch)
{
    uint64_t i;
    long long p_beg, p_end;
    bed_in* p = NULL;
    if(rev)
    {
        p_end = ref_p;
        p_beg = p_end + 1 - tLen;
    }
    else
    {
        p_beg = ref_p; 
        p_end = p_beg + tLen - 1;
    }
    if(p_beg < 0) p_beg = 0;
    if(p_end < 0) p_end = 0;

    p = &(t_ch->bed.a[uID]);
    for (i = 0; i < p->n; i++)
    {
        if(inter_interval(p_beg, p_end, p->a[i].beg, p->a[i].end, NULL, NULL)) break;
    }
    if(p->n > 0 && i < p->n) return 1;

    return 0;
}


void get_5_3_list(ha_ug_index* idx, s_hit* p, uint64_t cnt, s_hit** l5, uint64_t* l5_occ, 
s_hit** l3, uint64_t* l3_occ)
{
    (*l5) = (*l3) = NULL;
    (*l5_occ) = (*l3_occ) = 0;
    uint64_t i, j, rev, uID, ref_p, self_p, eLen, tLen, cur_beg, num;
    uint64_t beg_5 = (uint64_t)-1;
    for (j = 1, i = 0, num = 0; j <= cnt; ++j)
    {
        if(j == cnt || p[j].off_cnt != p[i].off_cnt)
        {
            interpret_pos((ha_ug_index*)idx, &p[i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
            ///cur_end = self_p;
            cur_beg = self_p + 1 - tLen;
            num++;
            if(cur_beg <= beg_5)
            {
                (*l3_occ) = (*l5_occ); (*l3) = (*l5);
                beg_5 = cur_beg; (*l5_occ) = j - i; (*l5) = p + i;
            }
            else
            {
                (*l3_occ) = j - i; (*l3) = p + i;
            }
            i = j;///must
        }
    }

    ///if(num > 2) fprintf(stderr, "ERROR: get_5_3_list\n");
}

inline void set_pe_pos_hap(ha_ug_index* idx, s_hit *l1, uint64_t occ1, s_hit *l2, uint64_t occ2, 
pe_hit_hap* x, uint64_t rid, trans_chain *t_ch)
{
    if(occ1 == 0 || occ2 == 0) return;
    uint64_t rev, uID, ref_p, self_p, eLen, tLen, i, is_unreliable = 0;
    s_hit *l1_5 = NULL, *l1_3 = NULL, *l2_5 = NULL, *l2_3 = NULL; 
    uint64_t l1_5_occ = 0, l1_3_occ = 0, l2_5_occ = 0, l2_3_occ = 0;

    /***************************for debug******************************/
    // fprintf(stderr, "\nrid: %lu, occ1: %lu, occ2: %lu\n", rid, occ1, occ2);
    // for (i = 0; i < occ1; i++)
    // {
    //     interpret_pos(idx, &l1[i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
    //     fprintf(stderr, "***-1-rev: %lu, uID: %lu, ref_p: %lu, self_p: %lu\n", 
    //     rev, uID, ref_p, self_p);
    // }

    // for (i = 0; i < occ2; i++)
    // {
    //     interpret_pos(idx, &l2[i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
    //     fprintf(stderr, "***-2-rev: %lu, uID: %lu, ref_p: %lu, self_p: %lu\n", 
    //     rev, uID, ref_p, self_p);
    // }
    /***************************for debug******************************/


    get_5_3_list(idx, l1, occ1, &l1_5, &l1_5_occ, &l1_3, &l1_3_occ);
    get_5_3_list(idx, l2, occ2, &l2_5, &l2_5_occ, &l2_3, &l2_3_occ);
    if(l1_5_occ == 0 || l2_5_occ == 0) return;
    x->id = rid;
    MALLOC(x->a, l1_5_occ + l2_5_occ);

    x->occ1 = 0; 
    for (i = 0; i < l1_5_occ; i++)
    {
        interpret_pos(idx, &l1_5[i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
        if(ref_p < self_p) continue;
        ref_p -= self_p; 
        if(rev) ref_p = idx->ug->u.a[uID].len - 1 - ref_p;
        if(t_ch && (is_unreliable_hits(rev, ref_p, tLen, uID, t_ch)))
        {
            is_unreliable = 1;
            continue;
        }
        x->a[x->occ1++] = (rev<<63) | ((uID << (64-idx->uID_bits))>>1) | (ref_p & idx->pos_mode);
    }
    
    x->occ2 = x->occ1; 
    for (i = 0; i < l2_5_occ; i++)
    {
        interpret_pos(idx, &l2_5[i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
        if(ref_p < self_p) continue;
        ref_p -= self_p; 
        if(rev) ref_p = idx->ug->u.a[uID].len - 1 - ref_p;
        if(t_ch && (is_unreliable_hits(rev, ref_p, tLen, uID, t_ch)))
        {
            is_unreliable = 1;
            continue;
        }
        x->a[x->occ2++] = (rev<<63) | ((uID << (64-idx->uID_bits))>>1) | (ref_p & idx->pos_mode);
    }
    x->occ2 -= x->occ1; 

    if(x->occ1 == 0 || x->occ2 == 0 || is_unreliable)
    {
        free(x->a); x->occ1 = x->occ2 = 0; x->a = NULL; x->id = (uint64_t)-1;
        return;
    }


    if(x->occ1 > 1) radix_sort_hc64(x->a, x->a + x->occ1);
    if(x->occ2 > 1) radix_sort_hc64(x->a + x->occ1, x->a + x->occ1 + x->occ2);


    /***************************for debug******************************/
    // fprintf(stderr, "-------------saved: x->occ1: %u, x->occ2: %u-------------\n", x->occ1, x->occ2);
    // for (i = 0; i < x->occ1; i++)
    // {
    //     fprintf(stderr, "###-1-rev: %lu, uID: %lu, ref_p: %lu\n", 
    //     x->a[i]>>63, (x->a[i]<<1)>>(64-idx->uID_bits), x->a[i] & idx->pos_mode);
    // }

    // for (i = 0; i < x->occ2; i++)
    // {
    //     fprintf(stderr, "###-2-rev: %lu, uID: %lu, ref_p: %lu\n", 
    //     x->a[i+x->occ1]>>63, (x->a[i+x->occ1]<<1)>>(64-idx->uID_bits), x->a[i+x->occ1] & idx->pos_mode);
    // }
    // fprintf(stderr, "-------------get_pe_s-rev: %lu, uID: %lu, ref_p: %lu-------------\n", 
    //     get_pe_s(*x)>>63, (get_pe_s(*x)<<1)>>(64-idx->uID_bits), get_pe_s(*x) & idx->pos_mode);
    // fprintf(stderr, "-------------get_pe_e-rev: %lu, uID: %lu, ref_p: %lu-------------\n", 
    //     get_pe_e(*x)>>63, (get_pe_e(*x)<<1)>>(64-idx->uID_bits), get_pe_e(*x) & idx->pos_mode);
    /***************************for debug******************************/
}

inline void set_pe_pos(ha_ug_index* idx, s_hit *l1, uint64_t occ1, s_hit *l2, uint64_t occ2, 
pe_hit* x, uint64_t rid, trans_chain* t_ch)
{
    if(occ1 == 0 || occ2 == 0) return;
    uint64_t rev, uID, ref_p, self_p, eLen, tLen, i, is_unreliable = 0;
    s_hit *l1_5 = NULL, *l1_3 = NULL, *l2_5 = NULL, *l2_3 = NULL; 
    uint64_t l1_5_occ = 0, l1_3_occ = 0, l2_5_occ = 0, l2_3_occ = 0;

    /***************************for debug******************************/
    // fprintf(stderr, "\nrid: %lu, occ1: %lu, occ2: %lu\n", rid, occ1, occ2);
    // for (i = 0; i < occ1; i++)
    // {
    //     interpret_pos(idx, &l1[i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
    //     fprintf(stderr, "***-1-rev: %lu, uID: %lu, ref_p: %lu, self_p: %lu\n", 
    //     rev, uID, ref_p, self_p);
    // }

    // for (i = 0; i < occ2; i++)
    // {
    //     interpret_pos(idx, &l2[i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
    //     fprintf(stderr, "***-2-rev: %lu, uID: %lu, ref_p: %lu, self_p: %lu\n", 
    //     rev, uID, ref_p, self_p);
    // }
    /***************************for debug******************************/


    get_5_3_list(idx, l1, occ1, &l1_5, &l1_5_occ, &l1_3, &l1_3_occ);
    get_5_3_list(idx, l2, occ2, &l2_5, &l2_5_occ, &l2_3, &l2_3_occ);
    if(l1_5_occ == 0 || l2_5_occ == 0) return;
    x->id = rid; x->len = 0;

    ///if(l1_5_occ != 1 || l2_5_occ != 1) fprintf(stderr, "ERROR\n");

    for (i = 0; i < l1_5_occ; i++)
    {
        interpret_pos(idx, &l1_5[i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
        if((ref_p + 1) < tLen) continue;
        ref_p = ref_p + 1 - tLen; 
        if(rev) ref_p = idx->ug->u.a[uID].len - 1 - ref_p;

        if(t_ch && (is_unreliable_hits(rev, ref_p, tLen, uID, t_ch)))
        {
            is_unreliable = 1;
            continue;
        }
        x->s = (rev<<63) | ((uID << (64-idx->uID_bits))>>1) | (ref_p & idx->pos_mode);
        x->len = tLen; x->len <<= 32;
    }
    
    for (i = 0; i < l2_5_occ; i++)
    {
        interpret_pos(idx, &l2_5[i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
        if((ref_p + 1) < tLen) continue;
        ref_p = ref_p + 1 - tLen;  
        if(rev) ref_p = idx->ug->u.a[uID].len - 1 - ref_p;

        if(t_ch && (is_unreliable_hits(rev, ref_p, tLen, uID, t_ch)))
        {
            is_unreliable = 1;
            continue;
        }
        x->e = (rev<<63) | ((uID << (64-idx->uID_bits))>>1) | (ref_p & idx->pos_mode);
        x->len |= tLen;
    }

    if(is_unreliable || x->s == (uint64_t)-1 || x->e == (uint64_t)-1)
    {
        x->id = x->s = x->e = x->len = (uint64_t)-1;
        return;
    }


    /****************************may have bugs********************************/
    // for (i = 0; i < l1_3_occ; i++)
    // {
    //     interpret_pos(idx, &l1_3[i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
    //     if(uID != ((x->s << 1) >> (64 - idx->uID_bits)) && 
    //                         uID != ((x->e << 1) >> (64 - idx->uID_bits)))
    //     {
    //         x->id = x->s = x->e = x->len = (uint64_t)-1;
    //         return;
    //     }
    // }

    // for (i = 0; i < l2_3_occ; i++)
    // {
    //     interpret_pos(idx, &l2_3[i], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
    //     if(uID != ((x->s << 1) >> (64 - idx->uID_bits)) && 
    //                         uID != ((x->e << 1) >> (64 - idx->uID_bits)))
    //     {
    //         x->id = x->s = x->e = x->len = (uint64_t)-1;
    //         return;
    //     }
    // }
    /****************************may have bugs********************************/



     
    /***************************for debug******************************/
    // fprintf(stderr, "-------------saved: x->occ1: %u, x->occ2: %u-------------\n", x->occ1, x->occ2);
    // for (i = 0; i < x->occ1; i++)
    // {
    //     fprintf(stderr, "###-1-rev: %lu, uID: %lu, ref_p: %lu\n", 
    //     x->a[i]>>63, (x->a[i]<<1)>>(64-idx->uID_bits), x->a[i] & idx->pos_mode);
    // }

    // for (i = 0; i < x->occ2; i++)
    // {
    //     fprintf(stderr, "###-2-rev: %lu, uID: %lu, ref_p: %lu\n", 
    //     x->a[i+x->occ1]>>63, (x->a[i+x->occ1]<<1)>>(64-idx->uID_bits), x->a[i+x->occ1] & idx->pos_mode);
    // }
    // fprintf(stderr, "-------------get_pe_s-rev: %lu, uID: %lu, ref_p: %lu-------------\n", 
    //     get_pe_s(*x)>>63, (get_pe_s(*x)<<1)>>(64-idx->uID_bits), get_pe_s(*x) & idx->pos_mode);
    // fprintf(stderr, "-------------get_pe_e-rev: %lu, uID: %lu, ref_p: %lu-------------\n", 
    //     get_pe_e(*x)>>63, (get_pe_e(*x)<<1)>>(64-idx->uID_bits), get_pe_e(*x) & idx->pos_mode);
    /***************************for debug******************************/
}


uint64_t if_debug_read(uint64_t rid)
{
    if(rid == 1169718 || rid == 2665829 || rid == 4239289)
    {
        return 1;
    }
    return 0;
}

static void worker_for_alignment(void *data, long i, int tid) // callback for kt_for()
{
    stepdat_t *s = (stepdat_t*)data;
    ///s->pos[i].id = (uint64_t)-1; s->pos[i].occ1 = s->pos[i].occ2 = 0; s->pos[i].a = NULL;
    s->pos[i].id = s->pos[i].s = s->pos[i].e = s->pos[i].len = (uint64_t)-1;

    /*******************************for debug************************************/
    // if(!if_debug_read(s->id+i)) return;
    // fprintf(stderr, "work-rid: %lu\n", (uint64_t)(s->id+i));
    /*******************************for debug************************************/

    uint64_t len1 = s->len[i]>>32, len2 = (uint32_t)s->len[i], occ1, occ2;
    char *r1 = s->seq[i], *r2 = s->seq[i] + len1;


    // fprintf(stderr, "**********R1**********\n");
    s->pos_buf[tid].a.n = 0;
    get_alignment(r1, len1, s->idx->k, &s->pos_buf[tid], s->idx, 0, s->id+i);
    occ1 = s->pos_buf[tid].a.n;
    if(occ1 == 0) return;


    // fprintf(stderr, "**********R2**********\n");
    get_alignment(r2, len2, s->idx->k, &s->pos_buf[tid], s->idx, occ1, s->id+i);
    occ2 = s->pos_buf[tid].a.n - occ1;
    if(occ2 == 0) return;   
    
    set_pe_pos((ha_ug_index*)s->idx, s->pos_buf[tid].a.a, occ1, s->pos_buf[tid].a.a + occ1, occ2, &(s->pos[i]), s->id+i, s->t_ch);

    /*******************************for debug************************************/
    // if(memcmp(r1, R1.r.a + R1.r_Len.a[s->id+i], len1) != 0)
    // {
    //     fprintf(stderr, "haha1\n");
    // }
    // if(memcmp(r2, R2.r.a + R2.r_Len.a[s->id+i], len2) != 0)
    // {
    //     fprintf(stderr, "haha2\n");
    // }
    // uint64_t j, rev, uID, ref_p, self_p, eLen, tLen;
    // char dir[2] = {'+', '-'};
    // fprintf(stderr, "(R1) %.*s\n", (int)(R1.name_Len.a[s->id + i + 1] - R1.name_Len.a[s->id+i]), 
    // R1.name.a + R1.name_Len.a[s->id+i]);
    // for (j = 0; j < occ1; j++)
    // {
    //     interpret_pos(s->idx, &s->pos_buf[tid].a.a[j], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
    //     fprintf(stderr, "utg%.6lu\t%c\t%lu\t%lu-%lu\n", uID+1, dir[rev], ref_p, self_p + 1 - tLen, self_p);
    // }
    // fprintf(stderr, "(R2) %.*s\n", (int)(R2.name_Len.a[s->id + i + 1] - R2.name_Len.a[s->id+i]), 
    // R2.name.a + R2.name_Len.a[s->id+i]);
    // for (j = 0; j < occ2; j++)
    // {
    //     interpret_pos(s->idx, &s->pos_buf[tid].a.a[j+occ1], &rev, &uID, &ref_p, &self_p, &eLen, &tLen);
    //     fprintf(stderr, "utg%.6lu\t%c\t%lu\t%lu-%lu\n", uID+1, dir[rev], ref_p, self_p + 1 - tLen, self_p);
    // }
    // fprintf(stderr, "\n");
    /*******************************for debug************************************/
}

static void *worker_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
    sldat_t *p = (sldat_t*)data;
    ///uint64_t total_base = 0, total_pair = 0;
    if (step == 0) { // step 1: read a block of sequences
        int ret1, ret2;
        uint64_t l1, l2;
		stepdat_t *s;
		CALLOC(s, 1);
        s->idx = p->idx; s->id = p->total_pair; s->t_ch = p->t_ch;
        while (((ret1 = kseq_read(p->ks1)) >= 0)&&((ret2 = kseq_read(p->ks2)) >= 0)) 
        {
            if (p->ks1->seq.l < p->idx->k || p->ks2->seq.l < p->idx->k) continue;
            if (s->n == s->m) {
                s->m = s->m < 16? 16 : s->m + (s->n>>1);
                REALLOC(s->len, s->m);
                REALLOC(s->seq, s->m);
            }

            l1 = p->ks1->seq.l; l2 = p->ks2->seq.l;
            MALLOC(s->seq[s->n], l1+l2);
            s->sum_len += l1+l2;
            memcpy(s->seq[s->n], p->ks1->seq.s, l1);
            memcpy(s->seq[s->n]+l1, p->ks2->seq.s, l2);
            s->len[s->n++] = (uint64_t)(l1<<32)|(uint64_t)l2;

            if (s->sum_len >= p->chunk_size) break;            
        }
        p->total_pair += s->n;
        if (s->sum_len == 0) free(s);
		else return s;
    }
    else if (step == 1) { // step 2: alignment
        stepdat_t *s = (stepdat_t*)in;
        CALLOC(s->pos_buf, p->n_thread);
        CALLOC(s->pos, s->n);
        int i;
        kt_for(p->n_thread, worker_for_alignment, s, s->n);
        for (i = 0; i < s->n; ++i) {
            free(s->seq[i]);
            p->total_base += (s->len[i]>>32) + (uint32_t)s->len[i];
        }
        
        free(s->seq); free(s->len);
        for (i = 0; i < (int)p->n_thread; ++i) {
            free(s->pos_buf[i].a.a);
        }
        free(s->pos_buf);
		return s;
    }
    else if (step == 2) { // step 3: dump
        stepdat_t *s = (stepdat_t*)in;
        int i;
        for (i = 0; i < s->n; ++i) {
            // if(s->pos[i].a == NULL) continue;
            // kv_push(pe_hit_hap, p->hits, s->pos[i]);
            if(s->pos[i].s == (uint64_t)-1) continue;
            kv_push(pe_hit, p->hits.a, s->pos[i]);
        }
        free(s->pos);
        free(s);
    }
    return 0;
}


int load_reads(reads_t* x, const enzyme *fn1, const enzyme *fn2)
{
    kv_init(x->name);
    kv_init(x->name_Len);
    kv_init(x->r);
    kv_init(x->r_Len);
    
    int ret;
    uint64_t name_tot, base_total;
    int i;
    name_tot = base_total = 0;

    for (i = 0; i < fn1->n && i < fn2->n; i++)
    {
        gzFile fp;
        if ((fp = gzopen(fn1->a[i], "r")) == 0)
        {
            kv_destroy(x->name);
            kv_destroy(x->name_Len);
            kv_destroy(x->r);
            kv_destroy(x->r_Len);
            return 0;
        } 
        
        kseq_t *ks;
        ks = kseq_init(fp);

        while (((ret = kseq_read(ks)) >= 0))
        {
            kv_push(uint64_t, x->name_Len, name_tot);
            kv_resize(char, x->name, name_tot + ks->name.l);
            memcpy(x->name.a + name_tot, ks->name.s, ks->name.l);
            name_tot += ks->name.l;

            kv_push(uint64_t, x->r_Len, base_total);
            kv_resize(char, x->r, base_total + ks->seq.l);
            memcpy(x->r.a + base_total, ks->seq.s, ks->seq.l);
            base_total += ks->seq.l;
        }
        
        kseq_destroy(ks);
        gzclose(fp);
    }
    kv_push(uint64_t, x->name_Len, name_tot);
    kv_push(uint64_t, x->r_Len, base_total);
    x->idx = 0;
    return 1;
}


void test_reads(reads_t* x, const char *fn)
{
    gzFile fp;
    kseq_t *ks;
    int ret, i = 0;

    if ((fp = gzopen(fn, "r")) == 0) return;
    ks = kseq_init(fp);
    while (((ret = kseq_read(ks)) >= 0))
    {
        if(memcmp(ks->name.s, x->name.a + x->name_Len.a[i], ks->name.l) != 0)
        {
            fprintf(stderr, "ERROR222: i: %d, len: %lu\n", i, x->name_Len.a[i]);
        }
        i++;
    }

    kseq_destroy(ks);
    gzclose(fp);
}

void destory_reads(reads_t* x)
{
    kv_destroy(x->name);
    kv_destroy(x->name_Len);
    kv_destroy(x->r);
    kv_destroy(x->r_Len);
}

void print_hits(ha_ug_index* idx, kvec_pe_hit* hits, const enzyme *fn1, const enzyme *fn2)
{
    uint64_t k, shif = 64 - idx->uID_bits;
    reads_t r1;
    load_reads(&r1, fn1, fn2);
    char dir[2] = {'+', '-'};
    for (k = 0; k < hits->a.n; ++k) 
    { 
        fprintf(stderr, "%.*s\t%c\ts-utg%.6dl\t%lu\t%c\te-utg%.6dl\t%lu\ti:%lu\n", 
        (int)(r1.name_Len.a[hits->a.a[k].id + 1] - r1.name_Len.a[hits->a.a[k].id]), 
        r1.name.a + r1.name_Len.a[hits->a.a[k].id],
        dir[hits->a.a[k].s>>63], (int)((hits->a.a[k].s<<1)>>shif)+1, hits->a.a[k].s&idx->pos_mode,
        dir[hits->a.a[k].e>>63], (int)((hits->a.a[k].e<<1)>>shif)+1, hits->a.a[k].e&idx->pos_mode,
        hits->a.a[k].id);        
    }
    destory_reads(&r1);
}

void print_hits_simp(ha_ug_index* idx, kvec_pe_hit* hits)
{
    uint64_t k, shif = 64 - idx->uID_bits;
    char dir[2] = {'+', '-'};
    for (k = 0; k < hits->a.n; ++k) 
    { 
        fprintf(stderr, "r-%lu-th\t%c\trs-utg%.6dl\t%lu\t%c\tre-utg%.6dl\t%lu\n", 
        hits->a.a[k].id,
        dir[hits->a.a[k].s>>63], (int)((hits->a.a[k].s<<1)>>shif)+1, hits->a.a[k].s&idx->pos_mode,
        dir[hits->a.a[k].e>>63], (int)((hits->a.a[k].e<<1)>>shif)+1, hits->a.a[k].e&idx->pos_mode);        
    }
}

inline void swap_pe_hit_hap(pe_hit_hap* x, pe_hit_hap* y)
{
    pe_hit_hap tmp;
    tmp = (*x); (*x) = (*y); (*y) = tmp;
}

void dedup_hits(kvec_pe_hit* hits, uint64_t is_dup)
{
    double index_time = yak_realtime();
    uint64_t k, l, m = 0, cur;
    radix_sort_pe_hit_an1(hits->a.a, hits->a.a + hits->a.n);
    for (k = 1, l = 0; k <= hits->a.n; ++k) 
    {   
        if (k == hits->a.n || hits->a.a[k].s != hits->a.a[l].s) 
        {
            if (k - l > 1) radix_sort_pe_hit_an2(hits->a.a + l, hits->a.a + k);
            if(is_dup)
            {
                cur = (uint64_t)-1;
                while (l < k)
                {
                    if(hits->a.a[l].e != cur)
                    {
                        cur = hits->a.a[l].e;
                        hits->a.a[m++] = hits->a.a[l];
                    }
                    l++;
                }
            }
            l = k;
        }
    }
    if(is_dup) hits->a.n = m;
    fprintf(stderr, "[M::%s::%.3f] ==> Dedup\n", __func__, yak_realtime()-index_time);
}

void sort_hits(kvec_pe_hit* hits)
{
    double index_time = yak_realtime();
    uint64_t k, l;
    radix_sort_pe_hit_an1(hits->a.a, hits->a.a + hits->a.n);
    for (k = 1, l = 0; k <= hits->a.n; ++k) 
    {   
        if (k == hits->a.n || (hits->a.a[k].s<<1) != (hits->a.a[l].s<<1)) 
        {
            if (k - l > 1) radix_sort_pe_hit_an2(hits->a.a + l, hits->a.a + k);
            l = k;
        }
    }
    fprintf(stderr, "[M::%s::%.3f] ==> Sort\n", __func__, yak_realtime()-index_time);
}

void destory_bubbles(bubble_type* bub)
{
    if(bub->index) free(bub->index);
    kv_destroy(bub->list);
    kv_destroy(bub->num);
    kv_destroy(bub->pathLen);
    kv_destroy(bub->b_s_idx);
    kv_destroy(bub->chain_weight);
    asg_destroy(bub->b_g);
    ma_ug_destroy(bub->b_ug);
}

void get_bubbles(bubble_type* bub, uint64_t id, uint32_t* beg, uint32_t* sink, uint32_t** a, uint32_t* n, uint64_t* pathBase)
{
    if(a) (*a) = bub->list.a + bub->num.a[id] + 2;
    if(n) (*n) = bub->num.a[id+1] - bub->num.a[id] - 2;
    if(beg) (*beg) = bub->list.a[bub->num.a[id]];
    if(sink) (*sink) = bub->list.a[bub->num.a[id] + 1];
    if(pathBase) (*pathBase) = bub->pathLen.a[id];
}

void dfs_bubble_broken(asg_t *g, kvec_t_u32_warp* stack, kvec_t_u32_warp* result, uint8_t* vis_flag, 
uint32_t vis_flag_n, uint32_t v_d, uint32_t beg_d, uint32_t sink_d)
{
    memset(vis_flag, 0, vis_flag_n);
    asg_arc_t *acur = NULL;
    uint32_t cur, ncur, i, p_beg = (uint32_t)-1, p_sink = (uint32_t)-1, v;
    stack->a.n = result->a.n = 0;
    v = v_d;
    if(v != (beg_d^1) && v != (sink_d^1)) kv_push(uint32_t, stack->a, v);
    while (stack->a.n > 0)
    {
        stack->a.n--;
        cur = stack->a.a[stack->a.n];
        vis_flag[cur] = 1;
        if((v>>1) != (cur>>1)) kv_push(uint32_t, result->a, cur>>1);
        ncur = asg_arc_n(g, cur);
        acur = asg_arc_a(g, cur);
        for (i = 0; i < ncur; i++)
        {
            if(acur[i].del) continue;
            if(vis_flag[acur[i].v]) continue;
            if((acur[i].v>>1) == (beg_d>>1) || (acur[i].v>>1) == (sink_d>>1))
            {
                if((acur[i].v>>1) == (beg_d>>1)) p_beg = acur[i].v;
                if((acur[i].v>>1) == (sink_d>>1)) p_sink = acur[i].v;
                continue;
            } 
            kv_push(uint32_t, stack->a, acur[i].v);
        }
    }
    
    memset(vis_flag, 0, vis_flag_n);
    v ^= 1;
    if(v != (beg_d^1) && v != (sink_d^1)) kv_push(uint32_t, stack->a, v);
    while (stack->a.n > 0)
    {
        stack->a.n--;
        cur = stack->a.a[stack->a.n];
        vis_flag[cur] = 1;
        if((v>>1) != (cur>>1)) kv_push(uint32_t, result->a, cur>>1);
        ncur = asg_arc_n(g, cur);
        acur = asg_arc_a(g, cur);
        for (i = 0; i < ncur; i++)
        {
            if(acur[i].del) continue;
            if(vis_flag[acur[i].v]) continue;
            if((acur[i].v>>1) == (beg_d>>1) || (acur[i].v>>1) == (sink_d>>1))
            {
                if((acur[i].v>>1) == (beg_d>>1)) p_beg = acur[i].v;
                if((acur[i].v>>1) == (sink_d>>1)) p_sink = acur[i].v;
                continue;
            } 
            kv_push(uint32_t, stack->a, acur[i].v);
        }
    }

    if(p_beg != (uint32_t)-1) kv_push(uint32_t, result->a, beg_d>>1);
    if(p_sink != (uint32_t)-1) kv_push(uint32_t, result->a, sink_d>>1);
}

void dfs_bubble(asg_t *g, kvec_t_u32_warp* stack, kvec_t_u32_warp* result, uint32_t v, uint32_t beg, uint32_t sink)
{
    asg_arc_t *acur = NULL;
    uint32_t cur, ncur, i;
    stack->a.n = result->a.n = 0;
    v = v << 1;
    kv_push(uint32_t, stack->a, v);
    while (stack->a.n > 0)
    {
        stack->a.n--;
        cur = stack->a.a[stack->a.n];
        if((v>>1) != (cur>>1)) kv_push(uint32_t, result->a, cur>>1);
        ncur = asg_arc_n(g, cur);
        acur = asg_arc_a(g, cur);
        for (i = 0; i < ncur; i++)
        {
            if(acur[i].del) continue;
            if((acur[i].v>>1) == beg || (acur[i].v>>1) == sink) continue;
            kv_push(uint32_t, stack->a, acur[i].v);
        }
    }
    

    v = v + 1;
    kv_push(uint32_t, stack->a, v);
    while (stack->a.n > 0)
    {
        stack->a.n--;
        cur = stack->a.a[stack->a.n];
        if((v>>1) != (cur>>1)) kv_push(uint32_t, result->a, cur>>1);
        ncur = asg_arc_n(g, cur);
        acur = asg_arc_a(g, cur);
        for (i = 0; i < ncur; i++)
        {
            if(acur[i].del) continue;
            if((acur[i].v>>1) == beg || (acur[i].v>>1) == sink) continue;
            kv_push(uint32_t, stack->a, acur[i].v);
        }
    }
}

uint32_t get_unitig_het_arb(ma_ug_t* ug, uint32_t uid, uint8_t *r_het_flag, kv_u_trans_t *ref, uint32_t m_het_occ,
uint32_t m_het_label, uint32_t p_het_label, uint32_t n_het_label)
{
    if(ref && ref->idx.n > 0 && u_trans_n(*ref, uid) > 0) return m_het_label;
    ma_utg_t *u = &(ug->u.a[uid]);
    uint32_t k, rId;
    uint32_t het_occ, hom_occ;
    for (k = 0, het_occ = hom_occ = 0; k < u->n; k++)
    {
        rId = u->a[k]>>33;
        if((r_het_flag[rId] & P_HET)) return m_het_label;
        if((r_het_flag[rId] & C_HET) || (r_het_flag[rId] & P_HET))
        {
            het_occ++;
        }
        else
        {
            hom_occ++;
        }
    }

    if((het_occ+hom_occ) == 0) return n_het_label; ///hom
    if((het_occ > ((het_occ+hom_occ)*0.8)) && ((het_occ+hom_occ) > m_het_occ)) return m_het_label; ///must het
    if(het_occ >= hom_occ) return p_het_label; ///potential het
    return n_het_label; ///hom
}
void update_bub_b_s_idx(bubble_type* bub);
void identify_bubbles(ma_ug_t* ug, bubble_type* bub, uint8_t *r_het_flag, kv_u_trans_t *ref)
{
    asg_cleanup(ug->g);
    if (!ug->g->is_symm) asg_symm(ug->g);
    uint32_t v, n_vtx = ug->g->n_seq * 2, i, k, mode = (((uint32_t)-1)<<2);
    uint32_t beg, sink, n, *a, n_occ;
    uint64_t pathLen;
    bub->ug = ug; 
    bub->b_bub = bub->b_end_bub = bub->tangle_bub = bub->cross_bub = bub->mess_bub = 0;
    if(bub->round_id == 0)
    {
        buf_t b; memset(&b, 0, sizeof(buf_t)); b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
        uint64_t tLen = get_bub_pop_max_dist_advance(ug->g, &b);
        kv_init(bub->list); kv_init(bub->num); kv_init(bub->pathLen);
        kv_init(bub->b_s_idx); kv_malloc(bub->b_s_idx, ug->g->n_seq); 
        bub->b_ug = NULL; kv_init(bub->chain_weight);
        bub->b_s_idx.n = ug->g->n_seq;
        memset(bub->b_s_idx.a, -1, bub->b_s_idx.n * sizeof(uint64_t));
        CALLOC(bub->index, n_vtx);
        for (i = 0; i < ug->g->n_seq; i++) ug->g->seq[i].c = 0;
        for (v = 0; v < n_vtx; ++v) 
        {
            if(ug->g->seq[v>>1].del) continue;
            if(asg_arc_n(ug->g, v) < 2) continue;
            if((bub->index[v]&(uint32_t)3) != 0) continue;
            if(asg_bub_pop1_primary_trio(ug->g, NULL, v, tLen, &b, (uint32_t)-1, (uint32_t)-1, 0, NULL, NULL, NULL, 0, 0, NULL))
            {
                //beg is v, end is b.S.a[0]
                //note b.b include end, does not include beg
                for (i = 0; i < b.b.n; i++)
                {
                    if(b.b.a[i]==v || b.b.a[i]==b.S.a[0]) continue;
                    bub->index[b.b.a[i]] &= mode; bub->index[b.b.a[i]] += 1;
                    bub->index[b.b.a[i]^1] &= mode; bub->index[b.b.a[i]^1] += 1;
                }
                bub->index[v] &= mode; bub->index[v] += 2;
                bub->index[b.S.a[0]^1] &= mode; bub->index[b.S.a[0]^1] += 3;
            }
        }

        kvec_t_u32_warp stack, result;
        kv_init(stack.a); kv_init(result.a);
        for (v = 0; v < n_vtx; ++v) 
        {
            if((bub->index[v]&(uint32_t)3) !=2) continue;
            if(asg_bub_pop1_primary_trio(ug->g, ug, v, tLen, &b, (uint32_t)-1, (uint32_t)-1, 0, &pathLen, NULL, NULL, 0, 0, NULL))
            {   
                //note b.b include end, does not include beg
                i = b.b.n + 1;
                if(b.b.n == 2 || b.b.n == 3 || b.b.n == 5)
                {
                    for (i = 0; i < b.b.n; i++)
                    {
                        if(b.b.a[i]==v || b.b.a[i]==b.S.a[0]) continue;
                        dfs_bubble(ug->g, &stack, &result, b.b.a[i]>>1, v>>1, b.S.a[0]>>1);
                        if((result.a.n + 3) != b.b.n && (result.a.n + 2) != b.b.n) break;
                    }
                }

                if(i == b.b.n)
                {
                    kv_push(uint32_t, bub->num, v);
                }
                else
                {
                    kv_push(uint32_t, bub->num, v + (1<<31));
                }
            }
        }
        kv_destroy(stack.a); kv_destroy(result.a);
        radix_sort_u32(bub->num.a, bub->num.a + bub->num.n);
        bub->s_bub = 0;
        for (k = 0; k < bub->num.n; k++)
        {
            if((bub->num.a[k]>>31) == 0) bub->s_bub++;
            v = (bub->num.a[k]<<1)>>1;
            bub->num.a[k] = bub->list.n;
            if(asg_bub_pop1_primary_trio(ug->g, ug, v, tLen, &b, (uint32_t)-1, (uint32_t)-1, 0, &pathLen, NULL, NULL, 0, 0, NULL))
            {
                kv_push(uint64_t, bub->pathLen, pathLen);
                //beg is v, end is b.S.a[0]
                kv_push(uint32_t, bub->list, v);
                kv_push(uint32_t, bub->list, b.S.a[0]^1);
                
                //note b.b include end, does not include beg
                for (i = 0; i < b.b.n; i++)
                {
                    if(b.b.a[i]==v || b.b.a[i]==b.S.a[0]) continue;
                    kv_push(uint32_t, bub->list, b.b.a[i]);
                }
            }
        }
        kv_push(uint32_t, bub->num, bub->list.n);
        free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
        bub->f_bub = bub->num.n - 1; ///bub->s_bub = bub->num.n - 1;
        
        for (i = 0; i < ug->g->n_seq; i++)
        {
            bub->index[i] = get_unitig_het_arb(ug, i, r_het_flag, ref, 20, M_het(*bub), P_het(*bub), (uint32_t)-1);
        }
        for (i = 0; i < bub->f_bub; i++)
        {
            get_bubbles(bub, i, &beg, &sink, &a, &n, &pathLen);
            for (v = n_occ = 0; v < n; v++)
            {
                bub->index[(a[v]>>1)] = i;
                n_occ += ug->u.a[a[v]>>1].n;
            }

            // if((pathLen*2) >= ug->g->seq[beg>>1].len && (pathLen*2) >= ug->g->seq[sink>>1].len)
            // {
            //     bub->index[(beg>>1)] = (uint32_t)-1;
            //     bub->index[(sink>>1)] = (uint32_t)-1;
            // }

            if(n_occ > 3)
            {
                if(bub->index[(beg>>1)] != M_het(*bub)) bub->index[(beg>>1)] = (uint32_t)-1;
                if(bub->index[(sink>>1)] != M_het(*bub)) bub->index[(sink>>1)] = (uint32_t)-1;
            }
            
            
            v = beg>>1;
            if(bub->b_s_idx.a[v] == (uint64_t)-1)
            {
                bub->b_s_idx.a[v] <<= 32;
                bub->b_s_idx.a[v] |= i;
            }
            else if((bub->b_s_idx.a[v] & 0xffffffff00000000) == 0xffffffff00000000)
            {
                bub->b_s_idx.a[v] <<= 32;
                bub->b_s_idx.a[v] |= i;
            }


            v = sink>>1;
            if(bub->b_s_idx.a[v] == (uint64_t)-1)
            {
                bub->b_s_idx.a[v] <<= 32;
                bub->b_s_idx.a[v] |= i;
            }
            else if((bub->b_s_idx.a[v] & 0xffffffff00000000) == 0xffffffff00000000)
            {
                bub->b_s_idx.a[v] <<= 32;
                bub->b_s_idx.a[v] |= i;
            }        
        }
        for (i = 0; i < ug->g->n_seq; i++)
        {
            if(bub->index[i] == M_het(*bub)) bub->index[i] = P_het(*bub);
        }
    }
    else
    {
        bub->num.n = bub->f_bub + 1; 
        bub->pathLen.n = bub->f_bub;
        bub->list.n = bub->num.a[bub->num.n-1];
        update_bub_b_s_idx(bub);
        bub->check_het = 0;
        asg_destroy(bub->b_g); bub->b_g = NULL;
        ma_ug_destroy(bub->b_ug); bub->b_ug = NULL;
        kv_destroy(bub->chain_weight); kv_init(bub->chain_weight);
    }
    bub->b_g = NULL;
    bub->b_ug = NULL;
    build_bub_graph(ug, bub);
    // fprintf(stderr, "-bub->index[18759]: %u, bub->num.n: %u\n",  (uint32_t)bub->index[18759], bub->num.n);
}

void print_bubbles(ma_ug_t* ug, bubble_type* bub, kvec_pe_hit* hits, hc_links* link, ha_ug_index* idx)
{
    uint64_t tLen, t_utg, i, k;
    uint32_t beg, sink, n, *a;
    for (i = 0, tLen = 0; i < bub->ug->u.n; i++) tLen += bub->ug->u.a[i].len;
    fprintf(stderr, "[M::%s] # unitigs: %lu, # bases: %lu\n",  __func__, bub->ug->u.n, tLen);
    for (i = 0, tLen = 0, t_utg = 0; i < bub->f_bub; i++)
    {
        get_bubbles(bub, i, &beg, &sink, &a, &n, NULL);
        t_utg += n;
        for (k = 0; k < n; k++)
        {
            tLen +=bub->ug->u.a[(a[k]>>1)].len;
        }
    }
    fprintf(stderr, "[M::%s] # bubbles: %lu, # unitigs: %lu, # bases: %lu\n",  __func__, 
    Get_bub_num(*bub), t_utg, tLen);


    for (i = 0, tLen = 0, t_utg = 0; i < ug->g->n_seq; i++)
    {
        if(IF_BUB(i, *bub))
        {
            t_utg++;
            tLen +=bub->ug->u.a[i].len;
        }
    }
    fprintf(stderr, "[M::%s] # bubbles: %lu, # unitigs: %lu, # bases: %lu\n",  __func__, 
    Get_bub_num(*bub), t_utg, tLen);

    for (i = 0, tLen = 0, t_utg = 0; i < ug->g->n_seq; i++)
    {
        if(IF_HET(i, *bub))
        {
            t_utg++;
            tLen +=bub->ug->u.a[i].len;
        }
    }
    fprintf(stderr, "[M::%s] # het unitigs: %lu, # het bases: %lu\n",  __func__, t_utg, tLen);

    uint8_t* flag; CALLOC(flag, ug->g->n_seq);
    uint64_t s_uid, e_uid, shif = 64 - idx->uID_bits;
    if(hits)
    {
        for (k = 0; k < hits->a.n; ++k) 
        {
            s_uid = ((hits->a.a[k].s<<1)>>shif);
            e_uid = ((hits->a.a[k].e<<1)>>shif);
            if(bub->index[s_uid] == (uint32_t)-1 || bub->index[e_uid] == (uint32_t)-1) continue;
            if(IF_BUB(s_uid, *bub) && IF_BUB(e_uid, *bub))
            {
                flag[s_uid] |= 1;
                flag[e_uid] |= 1;
                continue;
            }
            if(IF_HET(s_uid, *bub) && IF_HET(e_uid, *bub))
            {
                flag[s_uid] |= 4;
                flag[e_uid] |= 4;
                continue;
            }
            if(IF_BUB(s_uid, *bub)) flag[s_uid] |= 2, flag[e_uid] |= 2;
            if(IF_BUB(e_uid, *bub)) flag[e_uid] |= 2, flag[s_uid] |= 2;
        }
    }
    else if(link)
    {
        for (i = 0; i < link->a.n; i++)
        {
            for (k = 0; k < link->a.a[i].e.n; ++k) 
            {
                if(link->a.a[i].e.a[k].del) continue;
                s_uid = i;
                e_uid = link->a.a[i].e.a[k].uID;
                if(bub->index[s_uid] == (uint32_t)-1 || bub->index[e_uid] == (uint32_t)-1) continue;
                if(IF_BUB(s_uid, *bub) && IF_BUB(e_uid, *bub))
                {
                    flag[s_uid] |= 1;
                    flag[e_uid] |= 1;
                    continue;
                }
                if(IF_HET(s_uid, *bub) && IF_HET(e_uid, *bub))
                {
                    flag[s_uid] |= 4;
                    flag[e_uid] |= 4;
                    continue;
                }
                if(IF_BUB(s_uid, *bub)) flag[s_uid] |= 2, flag[e_uid] |= 2;
                if(IF_BUB(e_uid, *bub)) flag[e_uid] |= 2, flag[s_uid] |= 2;
            }
        }
    }
    



    for (i = 0, tLen = 0, t_utg = 0; i < ug->g->n_seq; i++)
    {
        if(flag[i] & (uint32_t)1)
        {
            t_utg++;
            tLen +=bub->ug->u.a[i].len;
        } 
    }
    fprintf(stderr, "[M::%s] # bubble-chained unitigs: %lu, # bubble-chained bases: %lu\n",  
                    __func__, t_utg, tLen);

    for (i = 0, tLen = 0, t_utg = 0; i < ug->g->n_seq; i++)
    {
        if((flag[i] & (uint32_t)1) || (flag[i] & (uint32_t)2))
        {
            t_utg++;
            tLen +=bub->ug->u.a[i].len;
        } 
    }
    fprintf(stderr, "[M::%s] # (bubble && het)-chained unitigs: %lu, # (bubble && het)-chained bases: %lu\n",  
                    __func__, t_utg, tLen);


    for (i = 0, tLen = 0, t_utg = 0; i < ug->g->n_seq; i++)
    {
        if((flag[i] & (uint32_t)1) || (flag[i] & (uint32_t)2) || (flag[i] & (uint32_t)4))
        {
            t_utg++;
            tLen +=bub->ug->u.a[i].len;
        } 
    }
    fprintf(stderr, "[M::%s] # (bubble || het)-chained unitigs: %lu, # (bubble || het)-chained bases: %lu\n",  
                    __func__, t_utg, tLen);
    free(flag);


    fprintf(stderr, "************bubble utgs************\n");
    uint64_t pathLen;
    for (i = 0, tLen = 0, t_utg = 0; i < bub->f_bub; i++)
    {
        get_bubbles(bub, i, &beg, &sink, &a, &n, &pathLen);
        t_utg += n;
        fprintf(stderr, "(full-%lu)\tbeg:utg%.6u\tsink:utg%.6u\tpathLen:%lu\t%s\n", 
        i, (beg>>1)+1, (sink>>1)+1, pathLen, i < bub->s_bub? "s-bub":(i<bub->f_bub?"f-bub":"b-bub"));
        for (k = 0; k < n; k++)
        {
            tLen +=bub->ug->u.a[(a[k]>>1)].len;
            fprintf(stderr, "utg%.6u,", (a[k]>>1)+1);
        }
        fprintf(stderr, "\n");
        ///if(i < bub->s_bub && (n != 4 && n != 2 && n != 1)) fprintf(stderr, "weird\n");
    }

    for (i = bub->f_bub, tLen = 0, t_utg = 0; i < bub->f_bub + bub->b_bub; i++)
    {
        get_bubbles(bub, i, &beg, &sink, &a, &n, &pathLen);
        t_utg += n;
        fprintf(stderr, "(broken-%lu)\tbeg:utg%.6u\tsink:utg%.6u\tpathLen:%lu\t%s\n", 
        i, (beg>>1)+1, (sink>>1)+1, pathLen, i < bub->s_bub? "s-bub":(i<bub->f_bub?"f-bub":"b-bub"));
        for (k = 0; k < n; k++)
        {
            tLen +=bub->ug->u.a[(a[k]>>1)].len;
            fprintf(stderr, "utg%.6u,", (a[k]>>1)+1);
        }
        fprintf(stderr, "\n");
    }

    // fprintf(stderr, "************het utgs************\n");
    // for (i = 0; i < ug->g->n_seq; i++)
    // {
    //     if(IF_HET(i, *bub)) fprintf(stderr, "utg%.6lu\n", i+1);
    // }
    // fprintf(stderr, "************het utgs************\n");
}

hc_edge* push_hc_edge(hc_linkeage* x, uint64_t uID, double weight, int dir, uint64_t* d)
{
    uint64_t k, n;
    hc_edge* a = NULL;
    hc_edge* p = NULL;
    if(dir == 0)
    {
        a = x->e.a;
        n = x->e.n;
    }
    else
    {
        a = x->f.a;
        n = x->f.n;
    }
    
    for (k = 0; k < n; k++)
    {
        if(a[k].del) continue;
        if(a[k].uID == uID)
        {
            a[k].weight += weight;
            if(d) a[k].dis = (*d);
            return &(a[k]);
        }
    }

    if(dir == 0)
    {
        kv_pushp(hc_edge, x->e, &p);
    }
    else
    {
        kv_pushp(hc_edge, x->f, &p);
    }
    
    ///p->del = p->enzyme = 0;
    p->del = 0;
    p->uID = uID;
    p->weight = weight;
    if(d) p->dis = (*d);
    return p;
}

long long get_enzyme_occ_debug(char* t, long long tlen, char* p, long long plen)
{
    long long s = 0, j, occ = 0;
    while(s <= (tlen - plen)) 
    { 
        j = plen-1; 

        while(j >= 0)
        {
            if(seq_nt4_table[(uint8_t)t[s+j]] >= 4) break;
            if((p[j] != t[s+j]) && seq_nt4_table[(uint8_t)p[j]] < 4) break;
            j--; 
        } 
        
        if (j < 0) occ++;
        s++;
    } 

    return occ;
}

int check_exact_match(char* x, long long xlen, char* y, long long ylen)
{
    long long i;
    if(xlen != ylen) return 0;
    for (i = 0; i < xlen; i++)
    {
        if(seq_nt4_table[(uint8_t)x[i]] >= 4) return 0;
        if((x[i] != y[i]) && seq_nt4_table[(uint8_t)y[i]] < 4) return 0;
    } 

    return 1;
}

long long get_enzyme_occ(char* t, long long tlen, char* p, long long plen)
{
    long long i, c, s = 0, j, occ = 0;
    int badchar[5]; badchar[0] = badchar[1] = badchar[2] = badchar[3] = badchar[4] = -1;
    for (i = 0; i < plen; i++)
    {
        c = seq_nt4_table[(uint8_t)p[i]];
        badchar[c] = i; 
        if(c == 4) badchar[0] = badchar[1] = badchar[2] = badchar[3] = i; 
    } 
    badchar[4] = -1;

    while(s <= (tlen - plen)) 
    { 
        j = plen-1; 

        while(j >= 0)
        {
            if(seq_nt4_table[(uint8_t)t[s+j]] >= 4) break;
            if((p[j] != t[s+j]) && seq_nt4_table[(uint8_t)p[j]] < 4) break;
            j--; 
        } 
        
            
        if (j < 0) 
        { 
            occ++;
            ///s += (s+m < n)? m-badchar[txt[s+m]] : 1; 
            s++;
        } 
        else
        {
            /*******************************for debug************************************/
            // long long f, end = s + MAX(1, j - badchar[seq_nt4_table[(uint8_t)t[s+j]]]);
            // for (f = s+1; f < end; f++)
            // {
            //     if(check_exact_match(t+f, plen, p, plen))
            //     {
            //         fprintf(stderr, "s: %lld, end: %lld, s+j: %lld, t[s+j]: %c, badchar: %d, j: %lld\n", 
            //                     s, end, s+j, t[s+j], badchar[seq_nt4_table[(uint8_t)t[s+j]]], j);
            //     }
            // }            
            /*******************************for debug************************************/
            
            s += MAX(1, j - badchar[seq_nt4_table[(uint8_t)t[s+j]]]); 
        }
    } 

    return occ;
}

#define pdq_cnt(q) ((q).x.a[0])

void init_pdq(pdq* q, uint64_t utg_num)
{
    kv_init(q->x); kv_push(uint64_t, q->x, 0);
    kv_malloc(q->dis, utg_num); q->dis.n = utg_num;
    kv_malloc(q->vis, utg_num); q->vis.n = utg_num;

    uint64_t i;
    for (i = 1; (uint64_t)(1<<i) < utg_num; i++);
    q->uID_mode = ((uint64_t)-1) >> (64-i); 
    q->uID_shift = i;
}

void destory_pdq(pdq* q)
{
    kv_destroy(q->x);
    kv_destroy(q->dis);
    kv_destroy(q->vis);
}

void reset_pdq(pdq* q)
{
    q->x.n = 1; pdq_cnt(*q) = 0;
    memset(q->dis.a, -1, sizeof(uint64_t)*q->dis.n);
    memset(q->vis.a, 0, sizeof(uint8_t)*q->vis.n);
}

void swap_pdq(uint64_t* i, uint64_t* j)
{
    uint64_t k;
    k = (*i);
    (*i) = (*j);
    (*j) = k;
}

#define weight(q, i) (get_dv_adv((q).x.a[i], (q).uID_mode, (q).uID_shift, &(q).tmp_v, &(q).tmp_d))

uint64_t inline set_dv_adv(uint64_t v, uint64_t dis, uint64_t v_mode, uint64_t v_shift)
{
    dis <<= v_shift; dis |= (v&v_mode);
    return dis; 
}

uint64_t inline get_dv_adv(uint64_t x, uint64_t v_mode, uint64_t v_shift, uint64_t* v, uint64_t* dis)
{
    (*v) = x & v_mode;
    (*dis) = x >> v_shift;
    return (*dis);
}

void push_pdq(pdq* q, uint64_t v, uint64_t dis)
{
    kv_push(uint64_t, q->x, set_dv_adv(v, dis, q->uID_mode, q->uID_shift));
    pdq_cnt(*q)++;
    int c_i = pdq_cnt(*q), p_i = c_i>>1;

    while ((p_i > 0) && (weight(*q, c_i) < weight(*q, p_i)))
    {
        swap_pdq(&(q->x.a[c_i]), &(q->x.a[p_i]));
        c_i  = p_i;
        p_i = c_i >> 1;
    }
}

void pop_pdq(pdq* q, uint64_t* min_v, uint64_t* min_dis)
{
    (*min_v) = (*min_dis) = (uint64_t)-1;
    if(pdq_cnt(*q) == 0) return;
    get_dv_adv((*q).x.a[1], (*q).uID_mode, (*q).uID_shift, min_v, min_dis);
    /*******************************for debug************************************/
    // uint64_t i;
    // for (i = 1; i < q->x.n; i++)
    // {
    //     if(weight(*q, i) < (*min_dis)) fprintf(stderr, "ERROR\n");
    // }
    /*******************************for debug************************************/
    ///min = q->x.a[1];
    swap_pdq(&(q->x.a[1]), &(q->x.a[pdq_cnt(*q)]));
    pdq_cnt(*q)--;
    q->x.n--;

    int c_i = 1, left_i, right_i, min_i, flag = 1;
    while(flag == 1) 
    {
        flag = 0;
        left_i = c_i << 1;
        right_i = left_i + 1;
        if(left_i > (int)(pdq_cnt(*q)))
        {
            break; // both children are null
        } 
        else if(right_i > (int)(pdq_cnt(*q)))
        {
            min_i = left_i; // right children is null
        }
        else
        {
            min_i = (weight(*q, left_i) < weight(*q, right_i))? left_i : right_i;
        }

        if(weight(*q, c_i) > weight(*q, min_i))
        {
            swap_pdq(&(q->x.a[c_i]), &(q->x.a[min_i]));
            c_i = min_i;
            flag = 1;
        }
    }
}

void get_shortest_path(uint32_t src, pdq* pq, asg_t *sg, uint32_t* pre)
{
    uint64_t v, u, i, nv, w;
    asg_arc_t *av = NULL;
    reset_pdq(pq);
    pq->dis.a[src] = 0;
    if(pre) pre[src] = (uint32_t)-1;
    push_pdq(pq, src, 0);
    while (pdq_cnt(*pq) > 0)
    {
        pop_pdq(pq, &v, &w);
        pq->vis.a[v] = 1;

        av = asg_arc_a(sg, v);
        nv = asg_arc_n(sg, v);

        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            u = av[i].v;
            w = (uint32_t)av[i].ul;

            if(pq->vis.a[u] == 0 && pq->dis.a[u] > pq->dis.a[v] + w)
            {
                pq->dis.a[u] = pq->dis.a[v] + w;
                push_pdq(pq, u, pq->dis.a[u]);
                if(pre) pre[u] = v;
            }
        }
    }
}

void all_pair_shortest_path(asg_t *sg, hc_links* link, MT* M)
{
    hc_linkeage* t = NULL;
    pdq pq;
    init_pdq(&pq, sg->n_seq<<1);
    uint32_t n_vtx = sg->n_seq<<1, v;
    uint64_t k, *p = NULL;

    for (v = 0; v < n_vtx; ++v)
    {
        if (sg->seq[v>>1].del) continue;
        t = &(link->a.a[v>>1]);
        if (t->e.n == 0) continue;
        get_shortest_path(v, &pq, sg, NULL);
        for (k = 0; k < pq.dis.n; k++)
        {
            if(pq.dis.a[k] == (uint64_t)-1) continue;
            kv_pushp(uint64_t, M->matrix.a[v].a, &p);
            (*p) = k << M->uID_shift;
            (*p) = (*p) | pq.dis.a[k];
        }
    }
    
    destory_pdq(&pq);
}

typedef struct{
    pdq* pq;
    uint32_t src;
    asg_t *sg;
    uint8_t *dest, flag;
    uint64_t occ, df_occ, v;
    uint32_t* pre;
}pdq_spec;

uint32_t get_specific_shortest_path(pdq_spec *p)
{
    uint64_t u, i, nv, w;
    asg_arc_t *av = NULL;
    if(p->v == (uint64_t)-1)
    {
        // reset_pdq(p->pq);
        p->pq->dis.a[p->src] = 0;
        if(p->pre) p->pre[p->src] = (uint32_t)-1;
        push_pdq(p->pq, p->src, 0);
        p->occ = 0;
    }
    if(p->occ >= p->df_occ) return 0;
    
    while (pdq_cnt(*(p->pq)) > 0)
    {
        pop_pdq(p->pq, &(p->v), &w);
        p->pq->vis.a[p->v] = 1;
        if(p->dest && p->dest[p->v] == p->flag) p->occ++;
        if(p->occ > p->df_occ) return 0;

        av = asg_arc_a(p->sg, p->v);
        nv = asg_arc_n(p->sg, p->v);

        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            u = av[i].v;
            w = (uint32_t)av[i].ul;

            if(p->pq->vis.a[u] == 0 && p->pq->dis.a[u] > p->pq->dis.a[p->v] + w)
            {
                p->pq->dis.a[u] = p->pq->dis.a[p->v] + w;
                push_pdq(p->pq, u, p->pq->dis.a[u]);
                if(p->pre) p->pre[u] = p->v;
            }
        }
        if(p->dest && p->dest[p->v] == p->flag) return 1;
        if(!p->dest) return 1;
    }

    return 0;
}
void get_utg_path(uint64_t s, uint64_t e, uint32_t *path, buf_t *res)
{
    uint64_t p = path[e], i;
    res->b.n = 0;
    kv_push(uint32_t, res->b, e);
    while (p != s)
    {
        kv_push(uint32_t, res->b, p);
        p = path[p];
    }
    kv_push(uint32_t, res->b, s);
    for (i = 0; i < (res->b.n>>1); i++)
    {
        p = res->b.a[i]; 
        res->b.a[i] = res->b.a[res->b.n-i-1];
        res->b.a[res->b.n-i-1] = p;
    }
    res->b.n--;
}
uint32_t check_trans_relation_by_path(uint32_t v, uint32_t w, pdq* pqv, uint32_t* path_v, buf_t *resv,
pdq* pqw, uint32_t* path_w, buf_t *resw, asg_t *sg, uint8_t *dest, uint8_t df, uint32_t df_occ, 
double rate, long long *dis)
{
    uint64_t r1, r2, d1, d2;
    pdq_spec p1, p2, *p = NULL;
    p1.v = (uint64_t)-1; p1.src = v; p1.pq = pqv; p1.sg = sg; p1.df_occ = df_occ;
    p1.occ = df_occ; p1.flag = df; p1.dest = dest; p1.pre = path_v; if(resv) resv->b.n = 0;
    reset_pdq(p1.pq);

    p2.v = (uint64_t)-1; p2.src = w; p2.pq = pqw; p2.sg = sg; p2.df_occ = df_occ;
    p2.occ = df_occ; p2.flag = df; p2.dest = dest; p2.pre = path_w; if(resw) resw->b.n = 0;
    reset_pdq(p2.pq);

    if(dis) (*dis) = 0;

    while (1)
    {
        p = &p1;
        r1 = get_specific_shortest_path(p);
        if(r1 && p1.pq->vis.a[p->v] && p2.pq->vis.a[p->v])
        {
            d1 = p1.pq->dis.a[p->v]; d2 = p2.pq->dis.a[p->v];
            if(d1 != (uint64_t)-1 && d2 != (uint64_t)-1)
            {
                if(d2 <= d1*(1+rate) && d2 >= d1*(1-rate))
                {
                    if(resv) get_utg_path(p1.src, p->v, p1.pre, resv);
                    if(resw) get_utg_path(p2.src, p->v, p2.pre, resw);
                    if(dis)
                    {
                        (*dis) = d1 + d2;
                        (*dis) -= (sg->seq[p1.src>>1].len + sg->seq[p2.src>>1].len);
                    }
                    
                    // fprintf(stderr, "p-utg%.6ul(d1-%lu), a-utg%.6ul(d2-%lu), conver-utg%.6lul\n", 
                    // (v>>1) + 1, d1, (w>>1) + 1, d2, (p->v>>1) + 1);


                    return 1;
                }
            }
        }
        p = &p2;
        r2 = get_specific_shortest_path(p);
        if(r2 && p1.pq->vis.a[p->v] && p2.pq->vis.a[p->v])
        {
            d1 = p1.pq->dis.a[p->v]; d2 = p2.pq->dis.a[p->v];
            if(d1 != (uint64_t)-1 && d2 != (uint64_t)-1)
            {
                if(d2 <= d1*(1+rate) && d2 >= d1*(1-rate))
                {
                    if(resv) get_utg_path(p1.src, p->v, p1.pre, resv);
                    if(resw) get_utg_path(p2.src, p->v, p2.pre, resw);
                    if(dis)
                    {
                        (*dis) = d1 + d2;
                        (*dis) -= (sg->seq[p1.src>>1].len + sg->seq[p2.src>>1].len);
                    }
                    // fprintf(stderr, "p-utg%.6ul(d1-%lu), a-utg%.6ul(d2-%lu), conver-utg%.6lul\n", 
                    // (v>>1) + 1, d1, (w>>1) + 1, d2, (p->v>>1) + 1);
                    return 1;
                }
            }
        }
        if(!r1 && !r2) return 0;
    }
    
}

void set_utg_by_dis(uint32_t v, pdq* pq, asg_t *g, kvec_t_u32_warp *res, uint32_t dis)
{
    uint64_t r;
    pdq_spec a;
    a.v = (uint64_t)-1; a.src = v; a.pq = pq; a.sg = g; a.df_occ = (uint64_t)-1;
    a.occ = 0; a.flag = (uint8_t)-1; a.dest = NULL; a.pre = NULL;
    reset_pdq(a.pq);
    while (1)
    {
        r = get_specific_shortest_path(&a);
        if(!r) return;
        if(a.pq->dis.a[a.v] <= dis) kv_push(uint32_t, res->a, a.v);
        else return;
    }
}

uint64_t LCA_distance(long long d_x, long long d_y, long long xLen, long long yLen, uint8_t* rev)
{
    (*rev) = 0;
    long long x_beg, x_end, y_beg, y_end, t_beg, t_end;
    x_end = d_x; x_beg = x_end - xLen + 1;
    y_end = d_y; y_beg = y_end - yLen + 1;
    if(x_end >= y_end)
    {
        t_end = x_end; (*rev) = 0;
        t_beg = y_beg;
    }
    else
    {
        t_end = y_end; (*rev) = 1;
        t_beg = x_beg;
    }

    return t_end + 1 - t_beg;
}

uint64_t get_LCA_bubble(uint32_t x, uint64_t xLen, uint32_t y, uint64_t yLen, uint8_t* dis, uint64_t n, MT* M, bubble_type* bub, uint64_t* min_rev)
{
    uint32_t j, v, k;
    uint64_t u, d = (uint64_t)-1, tmp;
    uint8_t rev;
    uint32_t root[2], a_n, *a;
    get_bubbles(bub, bub->index[x>>1], &root[0], &root[1], &a, &a_n, NULL);
    root[0] ^= 1; root[1] ^= 1;
    if(root[0] > root[1])
    {
        k = root[0];
        root[0] = root[1];
        root[1] = k;
    }

    dis[root[0]] = (uint8_t)-1;
    dis[root[1]] = (uint8_t)-1;

    v = x;
    for (j = 0; j < M->matrix.a[v].a.n; j++)
    {
        u = M->matrix.a[v].a.a[j] >> M->uID_shift;
        d = M->matrix.a[v].a.a[j] & M->dis_mode;
        dis[u] = dis[u] >> 4;
    }

    v = y;
    for (j = 0; j < M->matrix.a[v].a.n; j++)
    {
        u = M->matrix.a[v].a.a[j] >> M->uID_shift;
        d = M->matrix.a[v].a.a[j] & M->dis_mode;
        dis[u] = dis[u] >> 4;
    }
    uint64_t x_i = 0, y_i = 0, d_x, d_y, min_d = (uint64_t)-1;
    uint32_t min_j = (uint32_t)-1;
    (*min_rev) = (uint64_t)-1;


    for (k = 0; k < 2; k++)
    {
        j = root[k];
        if(dis[j] != 0)
        {
            dis[j] = (uint8_t)-1;
            continue;
        } 
        
        for (; x_i < M->matrix.a[x].a.n; x_i++)
        {
            u = M->matrix.a[x].a.a[x_i] >> M->uID_shift;
            d = M->matrix.a[x].a.a[x_i] & M->dis_mode;
            if(u == j) break;
        }
        if(x_i == M->matrix.a[x].a.n && M->matrix.a[x].a.n != 0) fprintf(stderr, "ERROR X\n");
        d_x = d;

        for (; y_i < M->matrix.a[y].a.n; y_i++)
        {
            u = M->matrix.a[y].a.a[y_i] >> M->uID_shift;
            d = M->matrix.a[y].a.a[y_i] & M->dis_mode;
            if(u == j) break;
        }
        if(y_i == M->matrix.a[y].a.n && M->matrix.a[y].a.n != 0) fprintf(stderr, "ERROR Y\n");
        d_y = d;

        tmp = LCA_distance(d_x, d_y, xLen, yLen, &rev);
        if(tmp < min_d) min_d = tmp, (*min_rev) = rev, min_j = j;
    }

    if(min_j == x || min_j == y) return (uint64_t)-1;

    return min_d;


}

uint64_t get_LCA(uint32_t x, uint64_t xLen, uint32_t y, uint64_t yLen, uint8_t* dis, uint64_t n, MT* M, bubble_type* bub, uint64_t* min_rev)
{
    if(IF_BUB(x>>1, *bub) && IF_BUB(y>>1, *bub) && bub->index[x>>1] == bub->index[y>>1])
    {
        return get_LCA_bubble(x, xLen, y, yLen, dis, n, M, bub, min_rev);
    }
    else
    {
        memset(dis, -1, sizeof(uint8_t)*n);
    }

    uint32_t j, v;
    uint64_t u, d = (uint64_t)-1, tmp;
    uint8_t rev;
    
    v = x;
    for (j = 0; j < M->matrix.a[v].a.n; j++)
    {
        u = M->matrix.a[v].a.a[j] >> M->uID_shift;
        d = M->matrix.a[v].a.a[j] & M->dis_mode;
        dis[u] = dis[u] >> 4;
    }

    v = y;
    for (j = 0; j < M->matrix.a[v].a.n; j++)
    {
        u = M->matrix.a[v].a.a[j] >> M->uID_shift;
        d = M->matrix.a[v].a.a[j] & M->dis_mode;
        dis[u] = dis[u] >> 4;
    }

    uint64_t x_i = 0, y_i = 0, d_x, d_y, min_d = (uint64_t)-1;
    uint32_t min_j = (uint32_t)-1;
    (*min_rev) = (uint64_t)-1;
    for (j = 0; j < n; j++)
    {
        if(dis[j] != 0)
        {
            dis[j] = (uint8_t)-1;
            continue;
        }

        for (; x_i < M->matrix.a[x].a.n; x_i++)
        {
            u = M->matrix.a[x].a.a[x_i] >> M->uID_shift;
            d = M->matrix.a[x].a.a[x_i] & M->dis_mode;
            if(u == j) break;
        }
        if(x_i == M->matrix.a[x].a.n && M->matrix.a[x].a.n != 0) fprintf(stderr, "ERROR X\n");
        d_x = d;

        for (; y_i < M->matrix.a[y].a.n; y_i++)
        {
            u = M->matrix.a[y].a.a[y_i] >> M->uID_shift;
            d = M->matrix.a[y].a.a[y_i] & M->dis_mode;
            if(u == j) break;
        }
        if(y_i == M->matrix.a[y].a.n && M->matrix.a[y].a.n != 0) fprintf(stderr, "ERROR Y\n");
        d_y = d;


        tmp = LCA_distance(d_x, d_y, xLen, yLen, &rev);
        if(tmp < min_d) min_d = tmp, (*min_rev) = rev, min_j = j;
    }

    if(min_j == x || min_j == y) return (uint64_t)-1;

    return min_d;
}

typedef struct { // data structure for each step in kt_pipeline()
    asg_t *sg;
    hc_links* link; 
    MT* M; 
    bubble_type* bub;
    uint8_t** dis_buf;
} utg_d_t;

static void worker_for_dis(void *data, long i, int tid) 
{
    utg_d_t* s = (utg_d_t*)data;
    hc_links* link = s->link;
    MT* M = s->M; 
    bubble_type* bub = s->bub;
    uint8_t* dis_buf = s->dis_buf[tid];
    asg_t *sg = s->sg;
    hc_linkeage* t = NULL;
    uint32_t n_vtx = sg->n_seq<<1, v, u, k, j;
    uint64_t d[2], db[2], q_u, min, min_i, min_b, rev[2], min_rev;

    if (sg->seq[i].del) return;
    t = &(link->a.a[i]);
    if (t->e.n == 0) return;

    for (k = 0; k < t->e.n; k++)
    {
        if(t->e.a[k].del) continue;
        if(i == t->e.a[k].uID) continue;
        u = t->e.a[k].uID;

        for (v = ((uint64_t)(i)<<1); v < ((uint64_t)(i+1)<<1); v++)///two directions
        {
            d[0] = d[1] = db[0] = db[1] = (uint64_t)-1;            
            for (j = 0; j < M->matrix.a[v].a.n; j++)
            {
                q_u = M->matrix.a[v].a.a[j] >> M->uID_shift;
                if((q_u>>1) == u) d[q_u&1] = (M->matrix.a[v].a.a[j] & M->dis_mode) + sg->seq[q_u>>1].len;
                if((q_u>>1) > u) break;///just for speeding up, doesn't affect results
            }

            min = min_i = min_b = (uint64_t)-1;
            if(t->e.a[k].dis != (uint64_t)-1) min = t->e.a[k].dis >> 3;
            
            if(d[0] < min) min = d[0], min_i = 0, min_b = 0;
            if(d[1] < min) min = d[1], min_i = 1, min_b = 0;

            if(min_i != (uint64_t)-1 && min != (uint64_t)-1)
            {
                t->e.a[k].dis = min<<1;
                t->e.a[k].dis += min_b;
                t->e.a[k].dis <<=1;
                t->e.a[k].dis += v&1;//s-direction
                t->e.a[k].dis <<=1;
                t->e.a[k].dis += min_i;//e-direction
            }
        }


        ///might be wrong
        if(IF_BUB(i, *bub) && IF_BUB(u, *bub)
            && bub->index[i] != bub->index[u] && t->e.a[k].dis != (uint64_t)-1)
        {
            continue;
        }

        for (v = ((uint64_t)(i)<<1); v < ((uint64_t)(i+1)<<1); v++)
        {
            d[0] = d[1] = db[0] = db[1] = (uint64_t)-1;  
            db[0] = get_LCA(v, sg->seq[v>>1].len, u<<1, sg->seq[u].len, 
                                                        dis_buf, n_vtx, M, bub, &rev[0]);
            db[1] = get_LCA(v, sg->seq[v>>1].len, (u<<1) + 1, sg->seq[u].len, 
                                                        dis_buf, n_vtx, M, bub, &rev[1]);

            min = min_i = min_b = min_rev = (uint64_t)-1;
            if(t->e.a[k].dis != (uint64_t)-1) min = t->e.a[k].dis >> 3;
            
            if(db[0] < min) min = db[0], min_i = 0, min_b = 1, min_rev = rev[0];
            if(db[1] < min) min = db[1], min_i = 1, min_b = 1, min_rev = rev[1];

            if(min_i != (uint64_t)-1 && min != (uint64_t)-1)
            {
                t->e.a[k].dis = min<<1;
                t->e.a[k].dis += min_b;
                t->e.a[k].dis <<=1;
                t->e.a[k].dis += ((v&1)^min_rev);//s-direction
                t->e.a[k].dis <<=1;
                t->e.a[k].dis += (min_i^min_rev);//e-direction
            }
        }
    }
}

void fill_utg_distance_multi(asg_t *sg, hc_links* link, MT* M, bubble_type* bub)
{
    // double index_time = yak_realtime();
    uint32_t i;
    utg_d_t s;
    s.sg = sg; s.link = link; s.M = M; s.bub = bub;
    s.dis_buf = (uint8_t**)malloc(sizeof(uint8_t*)*asm_opt.thread_num);
    for (i = 0; i < (uint32_t)asm_opt.thread_num; i++)
    {
        s.dis_buf[i] = (uint8_t*)malloc(sizeof(uint8_t)*(s.sg->n_seq<<1));
    }
    
    kt_for(asm_opt.thread_num, worker_for_dis, &s, s.sg->n_seq);


    for (i = 0; i < (uint32_t)asm_opt.thread_num; i++)
    {
        free(s.dis_buf[i]);
    }
    free(s.dis_buf);
    // fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
}

void init_MT(MT* M, uint32_t n_vtx)
{
    uint32_t v;
    kv_init(M->matrix); kv_malloc(M->matrix, n_vtx); M->matrix.n = n_vtx;
    for (v = 0; v < n_vtx; ++v) kv_init(M->matrix.a[v].a);
    for (v = 1; (uint64_t)(1<<v) < n_vtx; v++);
    M->uID_shift = 64 - v; M->dis_mode = ((uint64_t)-1) >> v; 
}

void destory_MT(MT* M)
{
    uint32_t v;
    for (v = 0; v < M->matrix.n; ++v) kv_destroy(M->matrix.a[v].a);
    kv_destroy(M->matrix);
}

int get_trans_ug_arch(uint32_t qn, uint32_t qs, uint32_t qe, uint32_t qLen, 
uint32_t tn, uint32_t ts, uint32_t te, uint32_t tLen, uint32_t rev, asg_arc_t* t)
{
    ma_hit_t h;
    h.qns = qn;
    h.qns = h.qns << 32;
    h.qns = h.qns | qs;
    h.qe = qe;
    h.tn = tn;
    h.ts = ts;
    h.te = te;
    h.rev = rev;
    h.del = 0;
    h.bl = h.el = h.ml = h.no_l_indel = 0;
    return ma_hit2arc(&h, qLen, tLen, MAX(qLen, tLen) + 1, 0, 0, t);
}

void update_ug_by_trans(asg_t *sg, kv_u_trans_t *ta)
{
    u_trans_t *a = NULL, *p = NULL;
    uint32_t i, k, st, occ, n, m;
    uint32_t qn, tn, qs, qe, ts, te, rev;
    asg_arc_t t, *e = NULL;
    long long r_qs, r_qe, r_ts, r_te;
    int r;

    for (k = occ = 0; k < ta->idx.n; k++)
    {
        a = u_trans_a(*ta, k);
        n = u_trans_n(*ta, k);
        for (st = 0, i = 1; i <= n; ++i)
        {
            if (i == n || a[i].tn != a[st].tn)
            {
                for (m = st, p = &(a[st]); m < i; m++)
                {
                    if(a[m].nw > p->nw) p = &(a[m]);
                }

                if(p->f == RC_2)///dis-connected
                {
                    rev = p->rev;
                    qn = p->qn; 
                    qs = p->qs; 
                    qe = p->qe - 1;
                    if(rev)
                    {
                        tn = p->tn;
                        ts = sg->seq[tn].len - (p->te - 1) - 1;
                        te = sg->seq[tn].len - p->ts - 1;
                    }
                    else
                    {
                        tn = p->tn; 
                        ts = p->ts; 
                        te = p->te - 1;
                    }

                    classify_hap_overlap(qs, qe, sg->seq[qn].len, ts, te, sg->seq[tn].len, 
                                                                            &r_qs, &r_qe, &r_ts, &r_te);
                    
                    qs = r_qs; qe = r_qe + 1;
                    if(rev)
                    {
                        ts = sg->seq[tn].len - r_te - 1;
                        te = sg->seq[tn].len - r_ts - 1 + 1;
                    }
                    else
                    {
                        ts = r_ts; te = r_te + 1;
                    }

                    r = get_trans_ug_arch(qn, qs, qe, sg->seq[qn].len, 
                                        tn, ts, te, sg->seq[tn].len, rev, &t);
                    if(r >= 0)
                    {
                        e = asg_arc_pushp(sg);
                        *e = t;
                        occ++;
                    }
                }
                st = i;
            }
        }
    }

    if(occ > 0)
    {
        free(sg->idx);
        sg->idx = 0;
        sg->is_srt = 0;
        asg_cleanup(sg);
    }
    asg_arc_del_trans(sg, asm_opt.gap_fuzz);

}

void push_LCA_edges(long long d_x, long long d_y, long long xLen, long long yLen, 
uint32_t v, uint32_t w, uint64_t *e0, uint64_t *e1)
{
    long long x_beg, x_end, y_beg, y_end;
    uint64_t d, rev;
    x_end = d_x; x_beg = x_end - xLen + 1;
    y_end = d_y; y_beg = y_end - yLen + 1;

    long long ovlp = ((MIN(x_end, y_end) >= MAX(x_beg, y_beg))? 
                                    MIN(x_end, y_end) - MAX(x_beg, y_beg) + 1 : 0);
    if(ovlp != xLen && ovlp != yLen)
    {
        d = MAX(x_end, y_end) - MIN(x_beg, y_beg) + 1;
        if(x_end >= y_end) rev = 1;
        else rev = 0;

        (*e0) = d;
        (*e0) <<= 1;
        (*e0) += 1;
        (*e0) <<= 1;
        (*e0) += ((v&1)^rev);//s-direction
        (*e0) <<= 1;
        (*e0) += ((w&1)^rev);//e-direction

        rev ^= 1;
        (*e1) = d;
        (*e1) <<= 1;
        (*e1) += 1;
        (*e1) <<= 1;
        (*e1) += ((w&1)^rev);//s-direction
        (*e1) <<= 1;
        (*e1) += ((v&1)^rev);//e-direction
    }
    else
    {
        if(xLen >= yLen) d = (x_beg + 1) + (y_end - y_beg + 1);
        else d = (x_end - x_beg + 1) + (yLen - y_end - 1);

        (*e0) = d;
        (*e0) <<= 1;
        (*e0) += 1;
        (*e0) <<= 1;
        (*e0) += (v&1);//s-direction
        (*e0) <<= 1;
        (*e0) += (w&1);//e-direction

        if(xLen >= yLen) d = (y_end - y_beg + 1) + (xLen - x_end - 1);
        else d = (y_beg + 1) + (x_end - x_beg + 1);
        (*e1) = d;
        (*e1) <<= 1;
        (*e1) += 1;
        (*e1) <<= 1;
        (*e1) += (w&1);//s-direction
        (*e1) <<= 1;
        (*e1) += (v&1);//e-direction
    }
}

void push_LCA_edges_rev(hc_edge *hx, hc_edge *hy, 
long long xLen, long long yLen, long long rLen, 
uint32_t v, uint32_t w, uint64_t *e0, uint64_t *e1)
{
    long long x_beg, x_end, y_beg, y_end, dx, dy;
    uint64_t d, rev;
    dx = (hx->dis>>3); 
    if((hx->dis&(uint64_t)2))
    {
        x_beg = rLen - dx - 1; x_end = x_beg + xLen - 1;
    }
    else
    {
        x_end = dx; x_beg = x_end - xLen + 1;
    } 
    
    dy = (hy->dis>>3);
    if((hy->dis&(uint64_t)2))
    {
        y_beg = rLen - dy - 1; y_end = y_beg + yLen - 1;
    } 
    else
    {
        y_end = dy; y_beg = y_end - yLen + 1;
    }

    long long ovlp = ((MIN(x_end, y_end) >= MAX(x_beg, y_beg))? 
                                    MIN(x_end, y_end) - MAX(x_beg, y_beg) + 1 : 0);
    if(ovlp != xLen && ovlp != yLen)
    {
        d = MAX(x_end, y_end) - MIN(x_beg, y_beg) + 1;
        if(x_end >= y_end) rev = 1;
        else rev = 0;

        (*e0) = d;
        (*e0) <<= 1;
        (*e0) += 1;
        (*e0) <<= 1;
        (*e0) += ((v&1)^rev);//s-direction
        (*e0) <<= 1;
        (*e0) += ((w&1)^rev);//e-direction

        rev ^= 1;
        (*e1) = d;
        (*e1) <<= 1;
        (*e1) += 1;
        (*e1) <<= 1;
        (*e1) += ((w&1)^rev);//s-direction
        (*e1) <<= 1;
        (*e1) += ((v&1)^rev);//e-direction
    }
    else
    {
        if(xLen >= yLen) d = (x_beg + 1) + (y_end - y_beg + 1);
        else d = (x_end - x_beg + 1) + (yLen - y_end - 1);
        (*e0) = d;
        (*e0) <<= 1;
        (*e0) += 1;
        (*e0) <<= 1;
        (*e0) += (v&1);//s-direction
        (*e0) <<= 1;
        (*e0) += (w&1);//e-direction

        if(xLen >= yLen) d = (y_end - y_beg + 1) + (xLen - x_end - 1);
        else d = (y_beg + 1) + (x_end - x_beg + 1);
        (*e1) = d;
        (*e1) <<= 1;
        (*e1) += 1;
        (*e1) <<= 1;
        (*e1) += (w&1);//s-direction
        (*e1) <<= 1;
        (*e1) += (v&1);//e-direction
    }
}

uint32_t up_contain(kv_u_trans_t *ta, hc_links* link, uint8_t *uc_idx, asg_t *sg, kvec_t_u64_warp *buf)
{
    hc_edge *he = NULL, *ht = NULL, *hx = NULL;
    uint64_t t_d = (uint64_t)-1, hd;
    uint32_t qn, tn, qs, qe, ts, te, rev, is_c;
    u_trans_t *a = NULL, *p = NULL;
    asg_arc_t t;
    uint32_t i, k, st, occ, n, m;
    long long r_qs, r_qe, r_ts, r_te;
    int r;

    for (k = 0, buf->a.n = 0; k < ta->idx.n; k++)
    {
        if((uc_idx[k]&2)&&(!(uc_idx[k]&1)))///contain others
        {
            if(link->a.a[k].e.n == 0)
            {
                uc_idx[k] -= 2;
                continue;
            }
            a = u_trans_a(*ta, k);
            n = u_trans_n(*ta, k);
            for (st = 0, i = 1; i <= n; ++i)
            {
                if (i == n || a[i].tn != a[st].tn)
                {
                    for (m = st, p = &(a[st]); m < i; m++)
                    {
                        if(a[m].nw > p->nw) p = &(a[m]);
                    }

                    if(p->f == RC_2)///dis-connected
                    {
                        rev = p->rev;
                        qn = p->qn; 
                        qs = p->qs; 
                        qe = p->qe - 1;
                        if(rev)
                        {
                            tn = p->tn;
                            ts = sg->seq[tn].len - (p->te - 1) - 1;
                            te = sg->seq[tn].len - p->ts - 1;
                        }
                        else
                        {
                            tn = p->tn; 
                            ts = p->ts; 
                            te = p->te - 1;
                        }

                        classify_hap_overlap(qs, qe, sg->seq[qn].len, ts, te, sg->seq[tn].len, 
                                                                                &r_qs, &r_qe, &r_ts, &r_te);
                        
                        qs = r_qs; qe = r_qe + 1;
                        if(rev)
                        {
                            ts = sg->seq[tn].len - r_te - 1;
                            te = sg->seq[tn].len - r_ts - 1 + 1;
                        }
                        else
                        {
                            ts = r_ts; te = r_te + 1;
                        }

                        r = get_trans_ug_arch(qn, qs, qe, sg->seq[qn].len, 
                                            tn, ts, te, sg->seq[tn].len, rev, &t);
                        if(r == MA_HT_TCONT)//q contains t
                        {
                            hd = tn; hd <<= 32; hd |= qn;
                            kv_push(uint64_t, buf->a, hd);
                            ///qn->tn
                            he = get_hc_edge(link, qn, tn, 0);
                            if(!he)
                            {
                                push_hc_edge(&(link->a.a[qn]), tn, 0, 0, &t_d);
                                he = get_hc_edge(link, qn, tn, 0);
                            } 
                            hd = (te - ts) + qs;
                            if(hd < he->dis)
                            {
                                he->dis = hd << 1;
                                he->dis += 1;
                                he->dis <<= 1;
                                he->dis += 0;//s-direction
                                he->dis <<= 1;
                                he->dis += rev;//e-direction
                            }

                            ///tn->qn
                            he = get_hc_edge(link, tn, qn, 0);
                            if(!he)
                            {
                                push_hc_edge(&(link->a.a[tn]), qn, 0, 0, &t_d);
                                he = get_hc_edge(link, tn, qn, 0);
                            } 
                            hd = (te - ts) + sg->seq[qn].len - qe;
                            if(hd < he->dis)
                            {
                                he->dis = hd << 1;
                                he->dis += 1;
                                he->dis <<= 1;
                                he->dis += 0;//s-direction
                                he->dis <<= 1;
                                he->dis += rev;//e-direction
                            }
                        }
                    }
                    st = i;
                }
            }
            uc_idx[k] -= 2;
        }
    }
    
    uint64_t e0, e1;
    for (k = 0; k < buf->a.n; k++)
    {
        qn = buf->a.a[k] >> 32;
        tn = (uint32_t)buf->a.a[k];

        he = get_hc_edge(link, tn, qn, 0); //tn contains qn
        for (i = 0; i < link->a.a[tn].e.n; i++)
        {
            if(link->a.a[tn].e.a[i].del) continue;
            if(link->a.a[tn].e.a[i].uID == qn) continue;
            ht = &(link->a.a[tn].e.a[i]);
            if(ht->dis == (uint64_t)-1) continue;

            if((he->dis&(uint64_t)2) == (ht->dis&(uint64_t)2))///s in same direction
            {
                push_LCA_edges(he->dis>>3, ht->dis>>3, sg->seq[he->uID].len, sg->seq[ht->uID].len, 
                he->dis&1, ht->dis&1, &e0, &e1);

                ///forward
                hx = get_hc_edge(link, he->uID, ht->uID, 0);
                if(!hx)
                {
                    push_hc_edge(&(link->a.a[he->uID]), ht->uID, 0, 0, &t_d);
                    hx = get_hc_edge(link, he->uID, ht->uID, 0);
                }
                if ((e0>>3) < (hx->dis>>3)) hx->dis = e0;

                ///backward
                hx = get_hc_edge(link, ht->uID, he->uID, 0);
                if(!hx)
                {
                    push_hc_edge(&(link->a.a[ht->uID]), he->uID, 0, 0, &t_d);
                    hx = get_hc_edge(link, ht->uID, he->uID, 0);
                }
                if ((e1>>3) < (hx->dis>>3)) hx->dis = e1;
            }
            else
            {
                push_LCA_edges_rev(he, ht, sg->seq[he->uID].len, sg->seq[ht->uID].len, 
                sg->seq[tn].len, he->dis&1, ht->dis&1, &e0, &e1);

                ///forward
                hx = get_hc_edge(link, he->uID, ht->uID, 0);
                if(!hx)
                {
                    push_hc_edge(&(link->a.a[he->uID]), ht->uID, 0, 0, &t_d);
                    hx = get_hc_edge(link, he->uID, ht->uID, 0);
                }
                if ((e0>>3) < (hx->dis>>3)) hx->dis = e0;

                ///backward
                hx = get_hc_edge(link, ht->uID, he->uID, 0);
                if(!hx)
                {
                    push_hc_edge(&(link->a.a[ht->uID]), he->uID, 0, 0, &t_d);
                    hx = get_hc_edge(link, ht->uID, he->uID, 0);
                }
                if ((e1>>3) < (hx->dis>>3)) hx->dis = e1;
            }
        }
    }


    for (k = 0, occ = 0; k < ta->idx.n; k++)
    {
        if((uc_idx[k]&1) && (uc_idx[k]&2))
        {
            is_c = 0;
            a = u_trans_a(*ta, k);
            n = u_trans_n(*ta, k);
            for (st = 0, i = 1; i <= n; ++i)
            {
                if (i == n || a[i].tn != a[st].tn)
                {
                    for (m = st, p = &(a[st]); m < i; m++)
                    {
                        if(a[m].nw > p->nw) p = &(a[m]);
                    }

                    if(p->f == RC_2)///dis-connected
                    {
                        rev = p->rev;
                        qn = p->qn; 
                        qs = p->qs; 
                        qe = p->qe - 1;
                        if(rev)
                        {
                            tn = p->tn;
                            ts = sg->seq[tn].len - (p->te - 1) - 1;
                            te = sg->seq[tn].len - p->ts - 1;
                        }
                        else
                        {
                            tn = p->tn; 
                            ts = p->ts; 
                            te = p->te - 1;
                        }

                        classify_hap_overlap(qs, qe, sg->seq[qn].len, ts, te, sg->seq[tn].len, 
                                                                                &r_qs, &r_qe, &r_ts, &r_te);
                        
                        qs = r_qs; qe = r_qe + 1;
                        if(rev)
                        {
                            ts = sg->seq[tn].len - r_te - 1;
                            te = sg->seq[tn].len - r_ts - 1 + 1;
                        }
                        else
                        {
                            ts = r_ts; te = r_te + 1;
                        }

                        r = get_trans_ug_arch(qn, qs, qe, sg->seq[qn].len, 
                                            tn, ts, te, sg->seq[tn].len, rev, &t);
                        if(r == MA_HT_QCONT && (uc_idx[tn]&2))//t contains q
                        {
                            is_c = 1;
                        }
                    }
                    st = i;
                }
            }
            if(is_c == 0) uc_idx[k] -= 1;
        }
        if(uc_idx[k]&2) occ++;
    }

    return occ;
}

void print_u_trans_t(u_trans_t *p)
{
    fprintf(stderr, "q-utg%.6ul\tqs(%u)\tqe(%u)\tt-utg%.6ul\tts(%u)\tte(%u)\trev(%u)\tw(%f)\n", p->qn+1, p->qs, p->qe, p->tn+1, p->ts, p->te, p->rev, p->nw);
}

void update_containment_distance(asg_t *sg, kv_u_trans_t *ta, hc_links* link)
{
    uint8_t *uc_idx = NULL;
    CALLOC(uc_idx, sg->n_seq);

    u_trans_t *a = NULL, *p = NULL;
    uint32_t i, k, st, n, m;
    uint32_t qn, tn, qs, qe, ts, te, rev;
    asg_arc_t t;
    long long r_qs, r_qe, r_ts, r_te;
    int r;

    for (k = 0; k < ta->idx.n; k++)
    {
        a = u_trans_a(*ta, k);
        n = u_trans_n(*ta, k);
        for (st = 0, i = 1; i <= n; ++i)
        {
            if (i == n || a[i].tn != a[st].tn)
            {
                for (m = st, p = &(a[st]); m < i; m++)
                {
                    if(a[m].nw > p->nw) p = &(a[m]);
                }

                if(p->f == RC_2)///dis-connected
                {
                    rev = p->rev;
                    qn = p->qn; 
                    qs = p->qs; 
                    qe = p->qe - 1;
                    if(rev)
                    {
                        tn = p->tn;
                        ts = sg->seq[tn].len - (p->te - 1) - 1;
                        te = sg->seq[tn].len - p->ts - 1;
                    }
                    else
                    {
                        tn = p->tn; 
                        ts = p->ts; 
                        te = p->te - 1;
                    }

                    classify_hap_overlap(qs, qe, sg->seq[qn].len, ts, te, sg->seq[tn].len, 
                                                                            &r_qs, &r_qe, &r_ts, &r_te);
                    
                    qs = r_qs; qe = r_qe + 1;
                    if(rev)
                    {
                        ts = sg->seq[tn].len - r_te - 1;
                        te = sg->seq[tn].len - r_ts - 1 + 1;
                    }
                    else
                    {
                        ts = r_ts; te = r_te + 1;
                    }

                    r = get_trans_ug_arch(qn, qs, qe, sg->seq[qn].len, 
                                        tn, ts, te, sg->seq[tn].len, rev, &t);
                    if(r == MA_HT_QCONT) uc_idx[qn] |= 1;
                    else if(r == MA_HT_TCONT) uc_idx[qn] |= 2;

                    // if(r < 0) print_u_trans_t(p);
                }
                st = i;
            }
        }
    }

    kvec_t_u64_warp buf; kv_init(buf.a);
    while(up_contain(ta, link, uc_idx, sg, &buf))
    {
        if(buf.a.n != 0) continue;
        for (k = 0; k < ta->idx.n; k++)
        {
            if((uc_idx[k]&1) && (uc_idx[k]&2))
            {
                uc_idx[k] = 2;
                break;
            }
        }
    }

    kv_destroy(buf.a);
    free(uc_idx);
}


void collect_hc_links(const ha_ug_index* idx, kvec_pe_hit* hits, hc_links* link, bubble_type* bub, MT* M)
{
    double index_time = yak_realtime();
    uint64_t k, i, shif = 64 - idx->uID_bits, beg, end, t_d;
    for (k = 0; k < hits->a.n; ++k) 
    {
        beg = ((hits->a.a[k].s<<1)>>shif);
        end = ((hits->a.a[k].e<<1)>>shif);

        if(beg == end) continue;
        if(IF_HOM(beg, *bub)) continue;
        if(IF_HOM(end, *bub)) continue;

        t_d = (uint64_t)-1;
        push_hc_edge(&(link->a.a[beg]), end, 0, 0, &t_d);
        push_hc_edge(&(link->a.a[end]), beg, 0, 0, &t_d);
    }
    asg_t *copy_sg = copy_read_graph(idx->ug->g);


    update_ug_by_trans(copy_sg, &(idx->t_ch->k_trans));
    all_pair_shortest_path(copy_sg, link, M);
    fill_utg_distance_multi(copy_sg, link, M, bub);
    update_containment_distance(copy_sg, &(idx->t_ch->k_trans), link);
    asg_destroy(copy_sg);

    fprintf(stderr, "[M::%s::%.3f] ==> Hi-C linkages have been counted\n", __func__, yak_realtime()-index_time);
    return;





    index_time = yak_realtime();
    for (k = 0; k < link->enzymes.n; k++)
    {
        link->enzymes.a[k] = 0;
        for (i = 0; i < (uint64_t)asm_opt.hic_enzymes->n; i++)
        {
            link->enzymes.a[k] += get_enzyme_occ(idx->ug->u.a[k].s, idx->ug->u.a[k].len, 
                                                asm_opt.hic_enzymes->a[i], asm_opt.hic_enzymes->l[i]);
        }
    }
    fprintf(stderr, "[M::%s::%.3f] ==> Enzymes have been counted\n", __func__, yak_realtime()-index_time);
}


void update_dis_connected_gfa(asg_t *sg, hc_links *link, MT *M)
{
    uint32_t i, m, v, nv, x, y;
    asg_arc_t *av = NULL;
	kvec_t(uint32_t) stack; kv_init(stack);
    uint32_t *flag = NULL; MALLOC(flag, sg->n_seq);
    uint64_t *group = NULL; MALLOC(group, sg->n_seq);

    for (i = 0; i < link->a.n; i++)
    {
        for (m = 0; m < link->a.a[i].e.n; m++)
        {
            if(link->a.a[i].e.a[m].dis == (uint64_t)-1)
            {
                if(link->a.a[i].e.a[m].is_cc != 0)
                {
                    fprintf(stderr, "ERROR 1\n");
                }
            }

            if(link->a.a[i].e.a[m].occ == (uint64_t)-1)
            {
                if(link->a.a[i].e.a[m].is_cc != 0 || link->a.a[i].e.a[m].dis == (uint64_t)-1)
                {
                    fprintf(stderr, "ERROR 2\n");
                }
            }
        }
    }


    memset(flag, -1, sg->n_seq*sizeof(uint32_t));
	// connected componets
	for (i = 0; i < sg->n_seq; ++i) {
		if (flag[i] != (uint32_t)-1) continue;
		stack.n = 0;
		kv_push(uint32_t, stack, i);
		while (stack.n > 0) {
			stack.n--;
			flag[stack.a[stack.n]] = i;///group id

            v = (stack.a[stack.n])<<1;
            av = asg_arc_a(sg, v);
            nv = asg_arc_n(sg, v);
			for (m = 0; m < nv; ++m) {
				if (flag[av[m].v>>1] != (uint32_t)-1) continue;
				kv_push(uint32_t, stack, av[m].v>>1);
			}

            v++;
            av = asg_arc_a(sg, v);
            nv = asg_arc_n(sg, v);
			for (m = 0; m < nv; ++m) {
				if (flag[av[m].v>>1] != (uint32_t)-1) continue;
				kv_push(uint32_t, stack, av[m].v>>1);
			}
		}
	}
	kv_destroy(stack);



	// precalculate the size of each group
	for (i = 0; i < sg->n_seq; ++i)
		group[i] = (uint64_t)flag[i] << 32 | i;
    radix_sort_hc64(group, group + sg->n_seq);
	for (i = 1, x = y = 0; i <= sg->n_seq; ++i) {
		if (i == sg->n_seq || (group[i]>>32) != (group[x]>>32)) {
			uint32_t j;
			for (j = x; j < i; ++j)
				group[j] = (uint64_t)y << 32 | (uint32_t)group[j];///(group id)|first element in this group
			++y, x = i;
		}
	}

    memset(flag, 0, sg->n_seq*sizeof(uint32_t));
    for (i = 1, x = y = 0; i <= sg->n_seq; ++i) 
    {
		if (i == sg->n_seq || (group[i]>>32) != (group[x]>>32)) 
        {
			x = i;
		}
	}


	free(flag);
}

void append_tig_link(uint64_t *cc, uint32_t cc_off, uint32_t cc_size, hc_links *link, asg_t *sg, kvec_t_u32_warp *buf)
{
    uint32_t i, k, id, tig_occ, b_cc_occ, *tig, *b_cc, qn, tn;
    uint64_t t_d = (uint64_t)-1;
    hc_edge *p = NULL;
    for (i = 0; i < cc_size; ++i) ///how many nodes
    {
        id = (uint32_t)cc[cc_off + i];///node id
        for (k = 0; k < link->a.a[id].e.n; k++)
        {
            if(link->a.a[id].e.a[k].is_cc == 0) break;
        }
        if(k < link->a.a[id].e.n) break;
    }
    if(i >= cc_size) return;

    buf->a.n = 0; tig_occ = 0;
    for (i = 0; i < cc_size; ++i) ///how many nodes
    {
        id = (uint32_t)cc[cc_off + i];///node id
        
        if(asg_arc_n(sg, (id<<1)) == 0 || asg_arc_n(sg, ((id<<1)+1)) == 0)///tig
        {
            kv_push(uint32_t, buf->a, id);
            tig_occ++;
        }
    }

    if(tig_occ == 0) return;

    b_cc_occ = 0;
    for (i = 0; i < cc_size; ++i) ///how many nodes
    {
        id = (uint32_t)cc[cc_off + i];///node id
        if(link->a.a[id].e.n == 0) continue;

        for (k = 0; k < link->a.a[id].e.n; k++)
        {
            if(link->a.a[id].e.a[k].is_cc == 0) break;
        }
        if(k < link->a.a[id].e.n)///disconnected
        {
            kv_push(uint32_t, buf->a, id);
            b_cc_occ++;
        }
    }

    if(b_cc_occ == 0) return;

    tig = buf->a.a;
    b_cc = buf->a.a + tig_occ;
    for (i = 0; i < tig_occ; i++)
    {
        qn = tig[i];
        for (k = 0; k < b_cc_occ; k++)
        {
            tn = b_cc[k];
            if(qn == tn) continue;
            p = get_hc_edge(link, qn, tn, 0);
            if(p) continue;

            t_d = (uint64_t)-1;
            p = push_hc_edge(&(link->a.a[qn]), tn, 0, 0, &t_d);
            p->is_cc = 0; p->occ = (uint64_t)-1;
            p = push_hc_edge(&(link->a.a[tn]), qn, 0, 0, &t_d);
            p->is_cc = 0; p->occ = (uint64_t)-1;
        }
    }
}

void update_ug_by_tigs(asg_t *sg, hc_links *link)
{
    uint32_t i, m, v, nv, x, y, qn, tn;
    asg_arc_t *av = NULL;
	kvec_t(uint32_t) stack; kv_init(stack);
    uint32_t *flag = NULL; MALLOC(flag, sg->n_seq);
    uint64_t *group = NULL; MALLOC(group, sg->n_seq);
    kvec_t_u32_warp buf; kv_init(buf.a);


    memset(flag, -1, sg->n_seq*sizeof(uint32_t));
	// connected componets
	for (i = 0; i < sg->n_seq; ++i) {
		if (flag[i] != (uint32_t)-1) continue;
		stack.n = 0;
		kv_push(uint32_t, stack, i);
		while (stack.n > 0) {
			stack.n--;
			flag[stack.a[stack.n]] = i;///group id

            v = (stack.a[stack.n])<<1;
            av = asg_arc_a(sg, v);
            nv = asg_arc_n(sg, v);
			for (m = 0; m < nv; ++m) {
				if (flag[av[m].v>>1] != (uint32_t)-1) continue;
				kv_push(uint32_t, stack, av[m].v>>1);
			}

            v++;
            av = asg_arc_a(sg, v);
            nv = asg_arc_n(sg, v);
			for (m = 0; m < nv; ++m) {
				if (flag[av[m].v>>1] != (uint32_t)-1) continue;
				kv_push(uint32_t, stack, av[m].v>>1);
			}
		}
	}
	kv_destroy(stack);



	// precalculate the size of each group
	for (i = 0; i < sg->n_seq; ++i)
		group[i] = (uint64_t)flag[i] << 32 | i;
    radix_sort_hc64(group, group + sg->n_seq);
	for (i = 1, x = y = 0; i <= sg->n_seq; ++i) {
		if (i == sg->n_seq || (group[i]>>32) != (group[x]>>32)) {
			uint32_t j;
			for (j = x; j < i; ++j)
				group[j] = (uint64_t)y << 32 | (uint32_t)group[j];///(group id)|first element in this group
			++y, x = i;
		}
	}


    for (i = 0; i < link->a.n; i++)
    {
        qn = i; 
        for (m = 0; m < link->a.a[i].e.n; m++)
        {
            link->a.a[i].e.a[m].occ = 0;
            link->a.a[i].e.a[m].is_cc = 0;
            tn = link->a.a[i].e.a[m].uID;
            if((group[qn]>>32) == (group[tn]>>32))
            {
                link->a.a[i].e.a[m].is_cc = 1;
            }
        }
    }


    for (i = 1, x = 0; i <= sg->n_seq; ++i) 
    {
		if (i == sg->n_seq || (group[i]>>32) != (group[x]>>32)) 
        {
            append_tig_link(group, x, i - x, link, sg, &buf);
			x = i;
		}
	}

    free(flag); kv_destroy(buf.a);
}

void idx_hc_links(kvec_pe_hit* hits, ha_ug_index* idx, bubble_type* bub);
void filter_disconnect_edges(ha_ug_index* idx, kvec_pe_hit *hits, hc_links *link, bubble_type *bub, uint32_t thres, double rate)
{
    uint32_t k, l, i, m, h_occ, *occ = NULL;
    uint64_t shif = 64 - idx->uID_bits, qn, tn, u_dis;
    pe_hit *h_a = NULL;
    hc_linkeage *t = NULL;
    u_trans_t *p = NULL;
    hc_edge *e = NULL;
    kv_u_trans_t k_trans; 
    kv_init(k_trans);
    
    if(hits->idx.n == 0) idx_hc_links(hits, idx, bub);
    CALLOC(occ, hits->idx.n);
    for (qn = 0; qn < hits->idx.n; qn++)
    {
        if(IF_HOM(qn, *bub)) continue;
        h_a = hits->a.a + (hits->idx.a[qn]>>32);
        h_occ = (uint32_t)(hits->idx.a[qn]);

        for (k = 1, l = 0; k <= h_occ; ++k) ///same qn
        {
            if (k == h_occ || ((h_a[k].e<<1)>>shif) != ((h_a[l].e<<1)>>shif)) //same qn and tn
            {
                tn = ((h_a[l].e<<1)>>shif);
                if(!IF_HOM(tn, *bub) && tn != qn)
                {
                    t = &(link->a.a[qn]);
                    for (i = 0, u_dis = (uint64_t)-1; i < t->e.n; i++)
                    {
                        if(t->e.a[i].del || t->e.a[i].uID != tn) continue;
                        u_dis = (t->e.a[i].dis ==(uint64_t)-1? (uint64_t)-1 : t->e.a[i].dis>>3);
                        break;
                    }
                    
                    if(u_dis == (uint64_t)-1)
                    {
                        kv_pushp(u_trans_t, k_trans, &p);
                        p->qn = qn; p->tn = tn; p->occ = (k-l);
                        kv_pushp(u_trans_t, k_trans, &p);
                        p->qn = tn; p->tn = qn; p->occ = (k-l);
                    }
                    else
                    {
                        occ[qn] += (k-l); occ[tn] += (k-l);
                    }
                }
                l = k;
            }
        }
    }


    radix_sort_u_trans_m(k_trans.a, k_trans.a + k_trans.n);

    for (k = 1, l = 0, m = 0; k <= k_trans.n; ++k)
    {
        if (k == k_trans.n || k_trans.a[l].qn != k_trans.a[k].qn || k_trans.a[l].tn != k_trans.a[k].tn) //same qn and tn
        {
            if(k - l > 2) fprintf(stderr, "ERROR-3\n");
            for (i = l, h_occ = 0; i < k; i++)
            {
                h_occ += k_trans.a[i].occ;
            }

            k_trans.a[m] = k_trans.a[l];
            // k_trans.a[m].occ = ((uint32_t)-1) - h_occ;
            k_trans.a[m].occ = h_occ;

            if(h_occ > thres || h_occ >= (occ[k_trans.a[m].qn]*rate) || h_occ >= (occ[k_trans.a[m].tn]*rate))
            {
                e = get_hc_edge(link, k_trans.a[m].qn, k_trans.a[m].tn, 0);
                if(e->dis != (uint64_t)-1) fprintf(stderr, "ERROR-3-0\n");
                e->is_cc = 1;

                e = get_hc_edge(link, k_trans.a[m].tn, k_trans.a[m].qn, 0);
                if(e->dis != (uint64_t)-1) fprintf(stderr, "ERROR-3-0\n");
                e->is_cc = 1;
            }
            m++;
            l = k;
        }
    }
    k_trans.n = m;

    kv_destroy(k_trans); free(occ);
}

void measure_distance(ha_ug_index* idx, const ma_ug_t* ug, kvec_pe_hit* hits, hc_links* link, bubble_type* bub, kv_u_trans_t *ta)
{
    // double index_time = yak_realtime();
    MT M;
    init_MT(&M, ug->g->n_seq<<1);
    uint64_t uID_bits;
    for (uID_bits=1; (uint64_t)(1<<uID_bits)<(uint64_t)ug->u.n; uID_bits++);
    uint64_t k, i, shif = 64 - uID_bits, beg, end, t_d;
    if(hits)
    {
        for (k = 0; k < hits->a.n; ++k) 
        {
            beg = ((hits->a.a[k].s<<1)>>shif);
            end = ((hits->a.a[k].e<<1)>>shif);

            if(beg == end) continue;
            if(IF_HOM(beg, *bub)) continue;
            if(IF_HOM(end, *bub)) continue;

            t_d = (uint64_t)-1;
            push_hc_edge(&(link->a.a[beg]), end, 0, 0, &t_d);
            push_hc_edge(&(link->a.a[end]), beg, 0, 0, &t_d);
        }
    }
    else
    {
        for (k = 0; k < ug->g->n_seq; ++k)
        {
            if(IF_HOM(k, *bub)) continue;
            for (i = 0; i < ug->g->n_seq; ++i)
            {
                if(i == k || IF_HOM(i, *bub)) continue;
                t_d = (uint64_t)-1;
                push_hc_edge(&(link->a.a[i]), k, 0, 0, &t_d);
                push_hc_edge(&(link->a.a[k]), i, 0, 0, &t_d);
            }
        }
    }

    asg_t *copy_sg = copy_read_graph(ug->g);

    update_ug_by_trans(copy_sg, ta);
    // update_ug_by_tigs(copy_sg, link);
    all_pair_shortest_path(copy_sg, link, &M);
    fill_utg_distance_multi(copy_sg, link, &M, bub);
    update_containment_distance(copy_sg, ta, link);
    // update_dis_connected_gfa(copy_sg, link, &M);
    asg_destroy(copy_sg);
    destory_MT(&M);

    for (i = 0; i < link->a.n; i++)
    {
        for (k = 0; k < link->a.a[i].e.n; k++)
        {
            link->a.a[i].e.a[k].is_cc = 0;
            if(link->a.a[i].e.a[k].del) continue;
            if(link->a.a[i].e.a[k].dis == (uint64_t)-1) continue;
            link->a.a[i].e.a[k].is_cc = 1;
        }
    }

    // filter_disconnect_edges(idx, hits, link, bub, 1, 0.01);
    // fprintf(stderr, "[M::%s::%.3f] ==> Hi-C linkages have been counted\n", __func__, yak_realtime()-index_time);
    return;
}

void set_reverse_links(uint32_t* bub, uint32_t n, kvec_t_u32_warp* reach, uint32_t root, hc_links* link)
{
    uint64_t i, k, d = RC_0;
    uint32_t v;
    for (i = 0; i < n; i++)
    {
        v = bub[i]>>1;
        if(v == root) continue;
        for (k = 0; k < reach->a.n; k++)
        {
            if(v == reach->a.a[k]) break;
        }

        ///if(k == reach->a.n && reach->a.n > 0)
        if(k == reach->a.n)
        {
            push_hc_edge(&(link->a.a[root]), v, 1, 1, &d);
            push_hc_edge(&(link->a.a[v]), root, 1, 1, &d);
        }
    }
    
}

void collect_hc_reverse_links(hc_links* link, ma_ug_t* ug, bubble_type* bub)
{
    uint64_t i, j, k, d = RC_0, m, pre;
    uint32_t beg, sink, n, v, *a = NULL;
    kvec_t_u32_warp stack, result;
    hc_edge *e = NULL;
    kv_init(stack.a); kv_init(result.a);
    ///clean all reverse overlaps within bubbles
    ///might be wrong
    for (i = 0; i < bub->f_bub; i++)
    {
        get_bubbles(bub, i, &beg, &sink, &a, &n, NULL);
        for (k = 0; k < n; k++)
        {
            v = a[k]>>1;
            for (j = 0; j < link->a.a[v].f.n; j++)
            {
                if(link->a.a[v].f.a[j].del) continue;
                e = get_hc_edge(link, link->a.a[v].f.a[j].uID, v, 1);
                e->del = 1;
            }
            link->a.a[v].f.n = 0;
        }

        v = beg>>1; 
        if(IF_HOM(v, *bub))
        {
            for (j = 0; j < link->a.a[v].f.n; j++)
            {
                if(link->a.a[v].f.a[j].del) continue;
                e = get_hc_edge(link, link->a.a[v].f.a[j].uID, v, 1);
                e->del = 1;
            }
            link->a.a[v].f.n = 0;
        }

        v = sink>>1;
        if(IF_HOM(v, *bub))
        {
            for (j = 0; j < link->a.a[v].f.n; j++)
            {
                if(link->a.a[v].f.a[j].del) continue;
                e = get_hc_edge(link, link->a.a[v].f.a[j].uID, v, 1);
                e->del = 1;
            }
            link->a.a[v].f.n = 0;
        }
    }

    for (i = 0; i < bub->f_bub; i++)
    {
        get_bubbles(bub, i, &beg, &sink, &a, &n, NULL);
        if(n == 2)
        {
            push_hc_edge(&(link->a.a[a[0]>>1]), a[1]>>1, 1, 1, &d);
            push_hc_edge(&(link->a.a[a[1]>>1]), a[0]>>1, 1, 1, &d);
            continue;
        }
        ///for complex bubbles, shouldn't have any assumption
        ///if(i >= bub->s_bub) continue;

        beg = beg>>1; sink = sink>>1;
        for (k = 0; k < n; k++)
        {
            v = a[k]>>1;
            dfs_bubble(ug->g, &stack, &result, v, beg, sink);
            set_reverse_links(a, n, &result, v, link);
        }
    }
    
    uint8_t* vis_flag = NULL;
    MALLOC(vis_flag, ug->g->n_seq*2);
    ///for broken bubbles
    for (i = bub->f_bub; i < bub->f_bub + bub->b_bub; i++)
    {
        get_bubbles(bub, i, &beg, &sink, &a, &n, NULL);

        for (k = 0; k < n; k++)
        {
            v = a[k];
            dfs_bubble_broken(ug->g, &stack, &result, vis_flag, ug->g->n_seq*2, v, beg, sink);
            set_reverse_links(a, n, &result, v>>1, link);
        }
    }
    
    kv_destroy(stack.a); kv_destroy(result.a); free(vis_flag);


    for (i = 0; i < link->a.n; i++)
    {
        for (k = m = 0; k < link->a.a[i].f.n; k++)
        {
            if(link->a.a[i].f.a[k].del) continue;
            link->a.a[i].f.a[m] = link->a.a[i].f.a[k];
            m++;
        }
        link->a.a[i].f.n = m;
        radix_sort_hc_edge_u(link->a.a[i].f.a, link->a.a[i].f.a + link->a.a[i].f.n);

        for (k = m = 0, pre = (uint64_t)-1; k < link->a.a[i].f.n; k++)
        {
            if(link->a.a[i].f.a[k].del) continue;
            if(link->a.a[i].f.a[k].uID == pre)
            {
                if(link->a.a[i].f.a[k].dis == RC_0) link->a.a[i].f.a[m-1].dis = RC_0;
                continue;
            }
             
            pre = link->a.a[i].f.a[k].uID;
            link->a.a[i].f.a[m] = link->a.a[i].f.a[k];
            m++;
        }
        link->a.a[i].f.n = m;
        radix_sort_hc_edge_d(link->a.a[i].f.a, link->a.a[i].f.a + link->a.a[i].f.n);
    }


    // hc_edge *e = NULL;
    // for (i = 0; i < link->a.n; i++)
    // {
    //     for (k = 0; k < link->a.a[i].f.n; k++)
    //     {
    //         if(link->a.a[i].f.a[k].del) continue;
    //         e = get_hc_edge(link, link->a.a[i].f.a[k].uID, i, 1);
    //         if(e == NULL) fprintf(stderr, "ERROR\n");
    //     }
    // }
    
}

void write_hc_links(hc_links* link, const char *fn)
{
    uint64_t k;
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
	sprintf(buf, "%s.hic.link.bin", fn);
    FILE* fp = fopen(buf, "w");

    fwrite(&link->a.n, sizeof(link->a.n), 1, fp);
    for (k = 0; k < link->a.n; k++)
    {
        fwrite(&link->a.a[k].e.n, sizeof(link->a.a[k].e.n), 1, fp);
        fwrite(link->a.a[k].e.a, sizeof(hc_edge), link->a.a[k].e.n, fp);

        fwrite(&link->a.a[k].f.n, sizeof(link->a.a[k].f.n), 1, fp);
        fwrite(link->a.a[k].f.a, sizeof(hc_edge), link->a.a[k].f.n, fp);
    }
    
    fwrite(&link->enzymes.n, sizeof(link->enzymes.n), 1, fp);
    fwrite(link->enzymes.a, sizeof(uint64_t), link->enzymes.n, fp);
    
    fclose(fp);
    free(buf);
    fprintf(stderr, "[M::%s::] ==> Hi-C linkages have been written\n", __func__);
}

int load_hc_links(hc_links* link, const char *fn)
{
    uint64_t k, flag = 0;
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
	sprintf(buf, "%s.hic.link.bin", fn);

    FILE* fp = NULL; 
    fp = fopen(buf, "r"); 
    if(!fp)
    {
        free(buf);
        return 0;
    } 
    

    kv_init(link->a);
    flag += fread(&link->a.n, sizeof(link->a.n), 1, fp);
    link->a.m = link->a.n; CALLOC(link->a.a, link->a.n);
    for (k = 0; k < link->a.n; k++)
    {
        flag += fread(&link->a.a[k].e.n, sizeof(link->a.a[k].e.n), 1, fp);
        link->a.a[k].e.m = link->a.a[k].e.n; MALLOC(link->a.a[k].e.a, link->a.a[k].e.n);
        flag += fread(link->a.a[k].e.a, sizeof(hc_edge), link->a.a[k].e.n, fp);

        flag += fread(&link->a.a[k].f.n, sizeof(link->a.a[k].f.n), 1, fp);
        link->a.a[k].f.m = link->a.a[k].f.n; MALLOC(link->a.a[k].f.a, link->a.a[k].f.n);
        flag += fread(link->a.a[k].f.a, sizeof(hc_edge), link->a.a[k].f.n, fp);
    }

    kv_init(link->enzymes);
    flag += fread(&link->enzymes.n, sizeof(link->enzymes.n), 1, fp);
    link->enzymes.m = link->enzymes.n; MALLOC(link->enzymes.a, link->enzymes.n);
    flag += fread(link->enzymes.a, sizeof(uint64_t), link->enzymes.n, fp);

    fclose(fp);
    free(buf);
    fprintf(stderr, "[M::%s::] ==> Hi-C linkages have been loaded\n", __func__);
    return 1;
}

void write_hc_hits(kvec_pe_hit* hits, ma_ug_t* ug, const char *fn)
{
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
    sprintf(buf, "%s.hic.lk.bin", fn);
    FILE* fp = fopen(buf, "w");

    fwrite(&hits->a.n, sizeof(hits->a.n), 1, fp);
    fwrite(hits->a.a, sizeof(pe_hit), hits->a.n, fp);
    write_dbug(ug, fp);

    fclose(fp);
    free(buf);
}

void write_hc_hits_v14(kvec_pe_hit_hap* i_hits, const char *fn)
{
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
	sprintf(buf, "%s.v14.hic.lk.bin", fn);
    FILE* fp = fopen(buf, "w");
    kvec_pe_hit hits; 
    kv_init(hits.a);
    uint64_t i, m_u = (uint64_t)-1, m_m = (uint64_t)-1;
    pe_hit* p = NULL;
    for (i = 0; i < i_hits->n; i++)
    {
        if(i_hits->a[i].occ1 == 1 && i_hits->a[i].occ2 == 1) 
        {
            kv_pushp(pe_hit, hits.a, &p);
            p->id = i_hits->a[i].id;
            p->s = i_hits->a[i].a[0];
            p->e = i_hits->a[i].a[1];
            m_u = i;
        }
        else
        {
            if(m_m == (uint64_t)-1) m_m = i;
        }
    }
    fprintf(stderr, "m_u: %lu, m_m: %lu, n_u: %lu\n", m_u, m_m, i_hits->n_u);

    fwrite(&hits.a.n, sizeof(hits.a.n), 1, fp);
    fwrite(hits.a.a, sizeof(pe_hit), hits.a.n, fp);

    kv_destroy(hits.a);
    fclose(fp);
    free(buf);
    exit(1);
}

#define pe_hit_hap_id_key(x) ((x).id)
KRADIX_SORT_INIT(pe_hit_hap_id, pe_hit_hap, pe_hit_hap_id_key, member_size(pe_hit_hap, id))

#define pe_hit_id_key(x) ((x).id)
KRADIX_SORT_INIT(pe_hit_id, pe_hit, pe_hit_id_key, member_size(pe_hit, id))

void debug_hc_hits_v14(kvec_pe_hit_hap* i_hits, const char *fn, const ha_ug_index* idx)
{
    uint64_t flag = 0;
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
	sprintf(buf, "%s.v14.hic.lk.bin", fn);
    kvec_pe_hit hits; 
    kv_init(hits.a);
    FILE* fp = NULL; 
    fp = fopen(buf, "r"); 

    kv_init(hits.a);
    flag += fread(&hits.a.n, sizeof(hits.a.n), 1, fp);
    hits.a.m = hits.a.n; MALLOC(hits.a.a, hits.a.n);
    flag += fread(hits.a.a, sizeof(pe_hit), hits.a.n, fp);

    radix_sort_pe_hit_id(hits.a.a, hits.a.a + hits.a.n);
    radix_sort_pe_hit_hap_id(i_hits->a, i_hits->a + i_hits->n_u);

    fprintf(stderr, "i_hits->n_u: %lu, hits.a.n: %lu\n", (uint64_t)i_hits->n_u, (uint64_t)hits.a.n);

    uint64_t i, k;
    uint64_t i_beg_utg, i_beg_pos, i_beg_rev;
    uint64_t i_end_utg, i_end_pos, i_end_rev;
    uint64_t k_beg_utg, k_beg_pos, k_beg_rev;
    uint64_t k_end_utg, k_end_pos, k_end_rev;
    uint64_t i_id, k_id;
    uint64_t same_occ = 0, diff_occ = 0, miss_occ = 0;
    for (i = 0, k = 0; i < i_hits->n_u; i++)
    {
        i_beg_rev = get_pe_s(i_hits->a[i])>>63;
        i_beg_utg = ((get_pe_s(i_hits->a[i])<<1)>>(64 - idx->uID_bits));
        i_beg_pos = get_pe_s(i_hits->a[i]) & idx->pos_mode;

        i_end_rev = get_pe_e(i_hits->a[i])>>63;
        i_end_utg = ((get_pe_e(i_hits->a[i])<<1)>>(64 - idx->uID_bits));
        i_end_pos = get_pe_e(i_hits->a[i]) & idx->pos_mode;

        i_id = i_hits->a[i].id;
        for (; k < hits.a.n; k++)
        {
            k_beg_rev = hits.a.a[k].s>>63;
            k_beg_utg = ((hits.a.a[k].s<<1)>>(64 - idx->uID_bits));
            k_beg_pos = hits.a.a[k].s & idx->pos_mode;

            k_end_rev = hits.a.a[k].e>>63;
            k_end_utg = ((hits.a.a[k].e<<1)>>(64 - idx->uID_bits));
            k_end_pos = hits.a.a[k].e & idx->pos_mode;

            k_id = hits.a.a[k].id;

            if(k_id > i_id)
            {
                miss_occ++;
                fprintf(stderr, "\n[MISS]rid=%lu\n", i_id);
                fprintf(stderr, "********v0.15********\n");
                fprintf(stderr, "beg_rev: %lu, beg_utg: %lu, beg_pos: %lu\n", 
                i_beg_rev, i_beg_utg, i_beg_pos);
                fprintf(stderr, "end_rev: %lu, end_utg: %lu, end_pos: %lu\n", 
                i_end_rev, i_end_utg, i_end_pos);
                break;
            } 
            
            if(k_id == i_id)
            {
                if(get_pe_s(i_hits->a[i]) == hits.a.a[k].s && get_pe_e(i_hits->a[i]) == hits.a.a[k].e)
                {
                    same_occ++;
                    // fprintf(stderr, "\n[SAME]rid=%lu\n", i_id);
                    // fprintf(stderr, "********v0.15********\n");
                    // fprintf(stderr, "beg_rev: %lu, beg_utg: %lu, beg_pos: %lu\n", 
                    // i_beg_rev, i_beg_utg, i_beg_pos);
                    // fprintf(stderr, "end_rev: %lu, end_utg: %lu, end_pos: %lu\n", 
                    // i_end_rev, i_end_utg, i_end_pos);
                    // fprintf(stderr, "********v0.14********\n");
                    // fprintf(stderr, "beg_rev: %lu, beg_utg: %lu, beg_pos: %lu\n", 
                    // k_beg_rev, k_beg_utg, k_beg_pos);
                    // fprintf(stderr, "end_rev: %lu, end_utg: %lu, end_pos: %lu\n", 
                    // k_end_rev, k_end_utg, k_end_pos);
                }
                else
                {
                    diff_occ++;
                    fprintf(stderr, "\n[DIFF]rid=%lu\n", i_id);
                    fprintf(stderr, "********v0.15********\n");
                    fprintf(stderr, "beg_rev: %lu, beg_utg: %lu, beg_pos: %lu\n", 
                    i_beg_rev, i_beg_utg, i_beg_pos);
                    fprintf(stderr, "end_rev: %lu, end_utg: %lu, end_pos: %lu\n", 
                    i_end_rev, i_end_utg, i_end_pos);
                    fprintf(stderr, "********v0.14********\n");
                    fprintf(stderr, "beg_rev: %lu, beg_utg: %lu, beg_pos: %lu\n", 
                    k_beg_rev, k_beg_utg, k_beg_pos);
                    fprintf(stderr, "end_rev: %lu, end_utg: %lu, end_pos: %lu\n", 
                    k_end_rev, k_end_utg, k_end_pos);
                }
                break;
            }
        }
    }
    
    fprintf(stderr, "same_occ: %lu, diff_occ: %lu, miss_occ: %lu", same_occ, diff_occ, miss_occ);


    kv_destroy(hits.a);
    fclose(fp);
    free(buf);
    exit(1);
}

int load_hc_hits(kvec_pe_hit* hits, ma_ug_t* ug, const char *fn)
{
    uint64_t flag = 0;
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
    sprintf(buf, "%s.hic.lk.bin", fn);

    FILE* fp = NULL; 
    fp = fopen(buf, "r"); 
    if(!fp) 
    {
        free(buf);
        return 0;
    }
    

    kv_init(hits->a);
    flag += fread(&hits->a.n, sizeof(hits->a.n), 1, fp);
    hits->a.m = hits->a.n; MALLOC(hits->a.a, hits->a.n);
    flag += fread(hits->a.a, sizeof(pe_hit), hits->a.n, fp);
    free(buf);

    if(!test_dbug(ug, fp))
    {
        free(hits->a.a);
        kv_init(hits->a);
        fclose(fp);
        fprintf(stderr, "[M::%s::] ==> Renew Hi-C linkages\n", __func__);
        return 0;
    }

    fclose(fp);
    fprintf(stderr, "[M::%s::] ==> Hi-C linkages have been loaded\n", __func__);
    return 1;
}

inline int get_phase_status(H_partition* hap, uint32_t uID)
{
    int d = -2;
    if(hap->hap[uID] & hap->m[0]) d = 1;
    if(hap->hap[uID] & hap->m[1]) d = -1;
    if(hap->hap[uID] & hap->m[2]) d = 0;
    return d;
}

inline uint32_t get_phase_group(H_partition* hap, uint32_t uID)
{
    return hap->hap[uID]>>hap->label_shift;
}

void print_hc_links(hc_links* link, int dir, H_partition* hap)
{
    uint64_t i, k;
    if(dir == 0)
    {
        double f_w, r_w;
        for (i = 0; i < link->a.n; ++i) 
        { 
            f_w = r_w = 0;
            for (k = 0; k < link->a.a[i].e.n; k++)
            {
                if(link->a.a[i].e.a[k].del) continue;
                fprintf(stderr, "s-utg%.6dl(%c)\tCLU:%d:%u\td-utg%.6dl(%c)\tCLU:%d:%u\t%lu\t%c\t%f\te\n", 
                (int)(i+1), "01"[!!(link->a.a[i].e.a[k].dis&(uint64_t)2)],
                get_phase_status(hap, i), hap->hap[i]>>3, 
                (int)(link->a.a[i].e.a[k].uID+1), "01"[!!(link->a.a[i].e.a[k].dis&(uint64_t)1)],
                get_phase_status(hap, link->a.a[i].e.a[k].uID), hap->hap[link->a.a[i].e.a[k].uID]>>3, 
                link->a.a[i].e.a[k].dis == (uint64_t)-1? (uint64_t)-1 : link->a.a[i].e.a[k].dis>>3,
                "fb"[!!(link->a.a[i].e.a[k].dis&(uint64_t)4)], link->a.a[i].e.a[k].weight);  
                if(get_phase_status(hap, i) == get_phase_status(hap, link->a.a[i].e.a[k].uID))
                {
                    f_w += link->a.a[i].e.a[k].weight;
                }  
                else
                {
                    r_w += link->a.a[i].e.a[k].weight;
                }
            }

            // fprintf(stderr, "self-utg%.6dl\tFW:%f\tRW:%f\tRT:%f\n**************************************************\n", 
            // (int)(i+1), f_w, r_w, (f_w+r_w) != 0? r_w/(f_w+r_w):0);
        }
    }
    

    if(dir == 1)
    {
        for (i = 0; i < link->a.n; ++i) 
        { 
            for (k = 0; k < link->a.a[i].f.n; k++)
            {
                if(link->a.a[i].f.a[k].del) continue;
                fprintf(stderr, "s-utg%.6d\td-utg%.6d\t%lu\te\n", 
                (int)(i+1), (int)(link->a.a[i].f.a[k].uID+1), link->a.a[i].f.a[k].dis);    
            }
        }
    }
}

void normalize_hc_links(hc_links* link)
{
    uint64_t i, k;
    for (i = 0; i < link->a.n; ++i) 
    { 
        for (k = 0; k < link->a.a[i].e.n; k++)
        {
            if(link->a.a[i].e.a[k].del) continue;
            link->a.a[i].e.a[k].weight *= 100;
            link->a.a[i].e.a[k].weight /= (double)(MIN(link->enzymes.a[i], link->enzymes.a[link->a.a[i].e.a[k].uID])); 
            ///link->a.a[i].e.a[k].weight /= (double)(link->enzymes.a[i] + link->enzymes.a[link->a.a[i].e.a[k].uID]); 
        }
    }
}

hc_edge* get_rGraph_edge(min_cut_t* x, uint64_t src, uint64_t dest)
{
    if(src >= x->rGraph.n) return NULL;
    uint64_t i;
    for (i = 0; i < x->rGraph.a[src].n; i++)
    {
        if(x->rGraph.a[src].a[i].del) continue;
        if(x->rGraph.a[src].a[i].uID == dest) return &(x->rGraph.a[src].a[i]);
    }
    return NULL;
}

void init_min_cut_t(min_cut_t* x, hc_links* link, const bubble_type* bub, const ma_ug_t *ug)
{
    uint64_t utg_num = link->a.n, i, k, u, v;
    x->n = utg_num;
    x->n_e = x->c_e = 0;

    kv_malloc(x->rGraphSet, utg_num); x->rGraphSet.n = utg_num;
    ///must utg_num<<1)
    kv_malloc(x->rGraphVis, utg_num); x->rGraphVis.n = utg_num;
    kv_malloc(x->utgVis, utg_num); x->utgVis.n = utg_num;
    kv_malloc(x->bmerVis, utg_num); x->bmerVis.n = utg_num;
    kv_malloc(x->order, utg_num); x->order.n = utg_num;
    kv_malloc(x->parent, utg_num); x->parent.n = utg_num;
    kv_malloc(x->p_weight, utg_num); x->p_weight.n = utg_num;
    ///uresolved BUGs, if use kv_resize segfault; if use kv_malloc, work?????
    kv_malloc(x->rGraph, utg_num); x->rGraph.n = utg_num;
    // kv_init(x->rGraph); kv_resize(hc_edge_warp, x->rGraph, utg_num); x->rGraph.n = utg_num;
    x->enzymes = link->enzymes.a;
    init_pdq(&(x->pq), utg_num<<1);

    //must be utg_num + 2 since we may need to add fake nodes
    for (i = 1; (uint64_t)(1<<i) < (utg_num + 2); i++);
    x->uID_mode = ((uint64_t)-1) >> (64-i); 
    x->uID_shift = i;
    for (i = 0; i < utg_num; i++)
    {
        ///x->order.a[i] = link->a.a[i].f.n;
        ///x->order.a[i] = ug->u.a[i].len;
        x->order.a[i] = x->enzymes[i];
        x->order.a[i] <<= x->uID_shift; 
        x->order.a[i] |= (uint64_t)(i & x->uID_mode);

        x->rGraphSet.a[i] = 0;
        x->rGraphVis.a[i] = 0;
        x->utgVis.a[i] = 0;
        x->bmerVis.a[i] = 0;
        x->parent.a[i] = (uint32_t)-1;
        
        ///uresolved BUGs, if use kv_resize segfault; if use kv_malloc, work?????
        // kv_init(x->rGraph.a[i]); kv_resize(hc_edge, x->rGraph.a[i], link->a.a[i].e.n);
        kv_malloc(x->rGraph.a[i], link->a.a[i].e.n); 
        x->rGraph.a[i].n = link->a.a[i].e.n;

        if(x->rGraph.a[i].n)
        {
            for (k = 0; k < x->rGraph.a[i].n; k++)
            {
                ///kv_push(hc_edge, x->rGraph.a[i], link->a.a[i].e.a[k]);
                x->rGraph.a[i].a[k] = link->a.a[i].e.a[k];
                x->n_e++;
                
                if((x->rGraph.a[i].a[k].weight == 0) || IF_HOM(x->rGraph.a[i].a[k].uID, *bub) 
                   || IF_HOM(i, *bub) || (x->rGraph.a[i].a[k].del))
                {
                    x->rGraph.a[i].a[k].del = 1;
                    x->n_e--;
                }
            }
        }
    }

    hc_edge *p = NULL;
    for (i = 0; i < utg_num; i++)
    {
        v = i;
        for (k = 0; k < link->a.a[v].f.n; k++)
        {
            if(link->a.a[v].f.a[k].del) continue;
            u = link->a.a[v].f.a[k].uID;


            p = get_rGraph_edge(x, v, u); 
            if(p)
            {
                p->del = 1;
                x->n_e--;
            }


            p = get_rGraph_edge(x, u, v);
            if(p)
            {
                p->del = 1;
                x->n_e--;
            } 
        }
    }

    x->q = kdq_init(uint64_t);
    radix_sort_hc64(x->order.a, x->order.a + x->order.n);

    x->b_mer = asm_opt.bub_mer_length;
    ///fprintf(stderr, "[M::%s]\n",  __func__);
    ///exit(0);
}

void destory_min_cut_t(min_cut_t* x)
{
    kv_destroy(x->order);
    kv_destroy(x->parent);
    kv_destroy(x->p_weight);
    kv_destroy(x->rGraphSet);
    kv_destroy(x->rGraphVis);
    kv_destroy(x->utgVis);
    kv_destroy(x->bmerVis);
    destory_pdq(&(x->pq));
    uint64_t i;
    for (i = 0; i < x->rGraph.m; i++)
    {
        kv_destroy(x->rGraph.a[i]);
    }
    kv_destroy(x->rGraph);
    kdq_destroy(uint64_t, x->q);
}

void reset_min_cut_t(min_cut_t* x, hc_links* link)
{
    ///no need to reset parent[] and q
    uint64_t i, j;
    ///important to have this line
    x->bmerVis.n = x->parent.n = x->p_weight.n = x->order.n = x->rGraph.n = x->rGraphVis.n = x->rGraphSet.n = link->a.n;
    kdq_clear(x->q);

    for (i = 0; i < x->rGraphSet.n; i++)
    {
        x->rGraphVis.a[i] = 0;
        ///x->bmerVis.a[i] = 0;
        ///important to have this line
        x->rGraph.a[i].n = link->a.a[i].e.n;

        if(x->rGraphSet.a[i] == 0) continue;
        for (j = 0; j < x->rGraph.a[i].n; j++)
        {
            x->rGraph.a[i].a[j].weight = link->a.a[i].e.a[j].weight;
        }
        x->rGraphSet.a[i] = 0;
    }
}

void update_link_by_min_cut_t(min_cut_t* x, hc_links* link)
{
    uint64_t i, j;

    for (i = 0; i < link->a.n; i++)
    {
        for (j = 0; j < link->a.a[i].e.n; j++)
        {
            link->a.a[i].e.a[j].del = x->rGraph.a[i].a[j].del;
        }
    }
}

uint64_t add_mul_convex(min_cut_t* x, uint64_t* a, uint64_t n)
{
    if(n == 0) return (uint64_t)-1;
    if(n == 1) return a[0];
    kv_push(uint8_t, x->rGraphSet, 0);
    kv_push(uint8_t, x->rGraphVis, 0);
    kv_push(uint8_t, x->bmerVis, 0);
    kv_push(uint32_t, x->parent, 0);
    kv_push(double, x->p_weight, 0);
    kv_resize(hc_edge_warp, x->rGraph, x->rGraph.n+1); 
    kv_init(x->rGraph.a[x->rGraph.n]);
    uint64_t i, k;
    hc_edge t;
    for (i = 0; i < n; i++)
    {
        ///t.uID = a[i]; t.del = t.enzyme = t.weight = 0;
        t.uID = a[i]; t.del = t.weight = 0;
        for (k = 0; k < x->rGraph.a[a[i]].n; k++)
        {
            if(x->rGraph.a[a[i]].a[k].del) continue;
            t.weight += x->rGraph.a[a[i]].a[k].weight;
        }
        kv_push(hc_edge, x->rGraph.a[x->rGraph.n], t);
        t.uID = x->rGraph.n;
        kv_push(hc_edge, x->rGraph.a[a[i]], t);
    }
    
    x->rGraph.n++;
    return x->rGraph.n - 1;
}

void get_s_t(min_cut_t* x, hc_links* link, uint64_t uID, uint64_t* src, uint64_t* dest, kvec_t_u64_warp* buff)
{
    buff->a.n = 0; (*src) = (*dest) = (uint64_t)-1;
    if(link->a.a[uID].f.n == 0) return;
    (*src) = uID; 

    uint64_t i, n;
    for (i = 0, n = 0; i < link->a.a[uID].f.n; i++)
    {
        if(link->a.a[uID].f.a[i].del) continue;
        kv_push(uint64_t, buff->a, link->a.a[uID].f.a[i].uID);
        (*dest) = link->a.a[uID].f.a[i].uID;
        n++;
    }
    if(n == 1 || n == 0) return;
    (*dest) = add_mul_convex(x, buff->a.a, buff->a.n);
}

uint64_t bfs_flow(uint64_t src, uint64_t dest, min_cut_t* x, kvec_t_u64_warp* buff)
{
    uint64_t *p = NULL, v, u, i;
    if(dest != (uint64_t)-1) memset(x->rGraphVis.a, 0, x->rGraphVis.n);

    kdq_push(uint64_t, x->q, src); 
    if(buff) kv_push(uint64_t, buff->a, src);

    x->rGraphVis.a[src] = 1;
    x->parent.a[src] = (uint32_t)-1;

    while (1)
    {
        p = kdq_shift(uint64_t, x->q);
        if(!p) break;
        v = *p;
        if(v == dest) return 1;
        for (i = 0; i < x->rGraph.a[v].n; i++)
        {
            if(x->rGraph.a[v].a[i].del) continue;
            if(x->rGraph.a[v].a[i].weight == 0) continue;
            u = x->rGraph.a[v].a[i].uID;
            if(x->rGraphVis.a[u]) continue;
            if(!x->bmerVis.a[u]) continue;

            x->parent.a[u] = v;
            x->p_weight.a[u] = x->rGraph.a[v].a[i].weight;

            kdq_push(uint64_t, x->q, u); 
            if(buff) kv_push(uint64_t, buff->a, u);
            ///set u or v to be 1? doesn't matter
            x->rGraphVis.a[u] = 1;
        }
    }

    return 0;
}

uint64_t maxFlow(uint64_t src, uint64_t dest, min_cut_t* x)
{
    double flow = 0, max_flow = 0;
    uint64_t v, u;
    hc_edge *p;

    while (bfs_flow(src, dest, x, NULL))
    {
        kdq_clear(x->q);
        flow = DBL_MAX;
        for (v = dest; v != src; v = x->parent.a[v])
        {
            flow = MIN(flow, x->p_weight.a[v]);
        }

        /*******************************for debug************************************/
        // if(src == 26818) fprintf(stderr, "***********flow: %f*********\n", flow);
        /*******************************for debug************************************/

        for (v = dest; v != src; v = x->parent.a[v])
        {
            u = x->parent.a[v];
            p = get_rGraph_edge(x, u, v);

            /*******************************for debug************************************/
            // if(src == 26818) fprintf(stderr, "utg%.6lul (%f)\n", u+1, p->weight);
            /*******************************for debug************************************/

            p->weight -= flow;
            p = get_rGraph_edge(x, v, u);
            p->weight += flow;
            x->rGraphSet.a[u] = x->rGraphSet.a[v] = 1;
        }

        max_flow += flow;
    }
    
    return (max_flow != 0);
}

uint64_t print_path(uint64_t src, uint64_t dest, min_cut_t* x)
{
    double flow = 0, max_flow = 0;
    uint64_t v, u;
    hc_edge *p;

    if(bfs_flow(src, dest, x, NULL))
    {
        kdq_clear(x->q);
        flow = DBL_MAX;
        for (v = dest; v != src; v = x->parent.a[v])
        {
            flow = MIN(flow, x->p_weight.a[v]);
        }

        /*******************************for debug************************************/
        fprintf(stderr, "***********flow: %f*********\n", flow);
        /*******************************for debug************************************/

        for (v = dest; v != src; v = x->parent.a[v])
        {
            u = x->parent.a[v];
            p = get_rGraph_edge(x, u, v);

            /*******************************for debug************************************/
            fprintf(stderr, "utg%.6lul (%f)\n", u+1, p->weight);
            /*******************************for debug************************************/
        }

        max_flow += flow;
    }
    
    return (max_flow != 0);
}

void print_src_dest(uint64_t src, min_cut_t* x, const char* command)
{
    uint64_t i;
    fprintf(stderr, "********************\n%s\n", command);
    if(src >= x->n)
    {
        for (i = 0; i < x->rGraph.a[src].n; i++)
        {
            if(x->rGraph.a[src].a[i].del) continue;
            fprintf(stderr, "utg%.6ul\n", x->rGraph.a[src].a[i].uID + 1);
        }
        
    }
    else
    {
        fprintf(stderr, "utg%.6lul\n", src+1);
    }
    fprintf(stderr, "!!!!!!!!!!!!!!!!!!!!\n");
    
}

void print_debug_rGraph(min_cut_t* x)
{
    fprintf(stderr, "******rGraph******\n");
    uint64_t i, j, u;
    for (i = 0; i < x->rGraphVis.n; i++)
    {
        if(!x->bmerVis.a[i]) continue;
        for (j = 0; j < x->rGraph.a[i].n; j++)
        {
            if(x->rGraph.a[i].a[j].del) continue;
            u = x->rGraph.a[i].a[j].uID;
            if(!x->bmerVis.a[u]) continue;
            fprintf(stderr, "***utg%.6lul\tutg%.6lul\t%f\n", i+1, u+1, x->rGraph.a[i].a[j].weight);    
        }
    }
    fprintf(stderr, "******rGraph******\n");
}

void graph_cut(uint64_t src, uint64_t dest, min_cut_t* x)
{
    /*******************************for debug************************************/
    ///if(src == 45179) print_debug_rGraph(x);
    /*******************************for debug************************************/
    if(maxFlow(src, dest, x))
    {
        ///in the last time bfs of maxFlow, rGraphVis has already been set
        uint64_t i, j, v, u;
        hc_edge *p;
        /*******************************for debug************************************/
        if(src == 45179)
        ///if(src == 26818)
        {
            ///print_debug_rGraph(x);
            print_src_dest(src, x, "src utg:");
            print_src_dest(dest, x, "dest utg:");
        }
        /*******************************for debug************************************/
        for (i = 0; i < x->rGraphVis.n; i++)
        {
            if(x->rGraphVis.a[i] == 0) continue;
            if(!x->bmerVis.a[i]) continue;
            v = i;
            for (j = 0; j < x->rGraph.a[i].n; j++)
            {
                if(x->rGraph.a[i].a[j].del) continue;
                u = x->rGraph.a[i].a[j].uID;
                if(x->rGraphVis.a[u]) continue;
                if(!x->bmerVis.a[u]) continue;
                /*******************************for debug************************************/
                if(src == 45179) fprintf(stderr, "utg%.6lul\tutg%.6lul\t%f\n", v+1, u+1, x->rGraph.a[i].a[j].weight);    
                /*******************************for debug************************************/
                ///delete <v, u>
                x->rGraph.a[i].a[j].del = 1;
                ///delete <u, v>
                p = get_rGraph_edge(x, u, v);
                p->del = 1;
                x->c_e += 2;
            }
            
        }

        /*******************************for debug************************************/
        ///if(src == 45179 || src == 31635)
        // if(src == 26818)
        // {
        //     fprintf(stderr, "hahahaha\n");
        //     print_path(26818, 1143, x);
        // }
        /*******************************for debug************************************/
    }
    /*******************************for debug************************************/
    ///if(src == 45179 || src == 31635)
    // {
    //     print_src_dest(src, x, "++++++src utg:");
    //     uint64_t m;
    //     for (m = 0; m < x->rGraph.a[src].n; m++)
    //     {
    //         if(x->rGraph.a[src].a[m].del) continue;
    //         fprintf(stderr, "src(utg%.6dl, enz:%lu)\tdes(utg%.6dl, enz:%lu)\t%f\n", 
    //         (int)(src+1), x->enzymes[src], 
    //         (int)(x->rGraph.a[src].a[m].uID+1), x->enzymes[x->rGraph.a[src].a[m].uID], 
    //         x->rGraph.a[src].a[m].weight);   
    //     }
    // }
    /*******************************for debug************************************/
}

void check_connective(min_cut_t* x, hc_links* link)
{
    double index_time = yak_realtime();
    kvec_t_u64_warp buff;
    kv_init(buff.a);
    uint64_t i, k, uID;
    for (i = 0; i < x->n; i++)
    {
        uID = x->order.a[i] & x->uID_mode;
        if(link->a.a[uID].f.n == 0) continue;
        for (k = 0; k < link->a.a[uID].f.n; k++)
        {
            if(link->a.a[uID].f.a[k].del) continue;
            if(x->utgVis.a[link->a.a[uID].f.a[k].uID] == 0) break;
        }
        if(k == link->a.a[uID].f.n) continue;
        reset_min_cut_t(x, link);
        get_s_t(x, link, uID, &(x->src), &(x->dest), &buff);

        bfs_flow(x->src, x->dest, x, NULL);

        x->utgVis.a[uID] = 1;
    }

    //reset x.utgVis
    memset(x->utgVis.a, 0, x->utgVis.n);
    kv_destroy(buff.a);
    fprintf(stderr, "[M::%s::%.3f] \n", __func__, yak_realtime()-index_time);
}

void get_Connected_Components(min_cut_t* x)
{
    double index_time = yak_realtime();
    uint64_t i, j, k = 0, uID, e;
    kvec_t_u64_warp buff;
    kv_init(buff.a);
    while (1)
    {
        for (i = 0; i < x->n; i++)
        {
            uID = x->order.a[i] & x->uID_mode;
            if(x->rGraphVis.a[uID] == 0) break;
        }
        if(i < x->n)
        {
            e = buff.a.n = 0;
            bfs_flow(uID, (uint64_t)-1, x, &buff);
            for (i = 0; i < buff.a.n; i++)
            {
                for (j = 0; j < x->rGraph.a[buff.a.a[i]].n; j++)
                {
                    if(x->rGraph.a[buff.a.a[i]].a[j].del == 0) e++;
                }
            }
            e >>= 1;
            if(buff.a.n > 1)
            {
                fprintf(stderr, "(%lu) Component: # nodes: %lu, # edges: %lu\n", 
                                                                k, (uint64_t)buff.a.n, e);
            }
            k++;
        }
        else
        {
            break;
        }
    }

    kv_destroy(buff.a);
    fprintf(stderr, "[M::%s::%.3f] # Connected Components: %lu\n", 
    __func__, yak_realtime()-index_time, k);
}

void print_rGraph(min_cut_t* x)
{
    uint64_t i, k;
    for (i = 0; i < x->rGraph.n; ++i) 
    { 
        for (k = 0; k < x->rGraph.a[i].n; k++)
        {
            if(x->rGraph.a[i].a[k].del) continue;
            fprintf(stderr, "src(utg%.6dl, enz:%lu)\tdes(utg%.6dl, enz:%lu)\t%f\n", 
            (int)(i+1), x->enzymes[i], 
            (int)(x->rGraph.a[i].a[k].uID+1), x->enzymes[x->rGraph.a[i].a[k].uID], 
            x->rGraph.a[i].a[k].weight);    
        }
    }
}

int select_large_node(const ma_ug_t *ug, min_cut_t* x, 
uint64_t src, uint64_t dest, uint64_t utg_thres, int weight_thres)
{
    if(src >= ug->u.n || dest >= ug->u.n) return 0;
    if(ug->u.a[src].n < utg_thres || ug->u.a[dest].n < utg_thres) return 0;

    uint64_t k;
    for (k = 0; k < x->rGraph.a[src].n; k++)
    {
        if(x->rGraph.a[src].a[k].del) continue;
        if(x->rGraph.a[src].a[k].weight >= weight_thres) break;
    }
    if(k == x->rGraph.a[src].n) return 0;

    src = dest;
    for (k = 0; k < x->rGraph.a[src].n; k++)
    {
        if(x->rGraph.a[src].a[k].del) continue;
        if(x->rGraph.a[src].a[k].weight >= weight_thres) break;
    }
    if(k == x->rGraph.a[src].n) return 0;

    return 1;
}

uint64_t inline set_dv(uint64_t v, uint64_t dis)
{
    dis <<= 32; dis |= v;
    return dis; 
}

uint64_t select_bmer(uint32_t src, uint64_t k, const bubble_type* bub, min_cut_t* x, uint32_t bub_only)
{
    uint32_t beg, sink, n, *a;
    uint32_t v, d, u, i, nv, b_mer_d, j;
    asg_t *sg = bub->ug->g;
    uint64_t *p = NULL;
    asg_arc_t *av = NULL;

    memset(x->rGraphVis.a, 0, x->rGraphVis.n);
    kdq_push(uint64_t, x->q, set_dv(src , 0));
    b_mer_d = 0; 

    x->rGraphVis.a[src] = 1;
    x->bmerVis.a[src] = 1;

    while (1)
    {
        p = kdq_shift(uint64_t, x->q);
        if(!p) break;
        v = (uint32_t)(*p); d = ((uint64_t)(*p))>>32;

        v = v<<1;
        av = asg_arc_a(sg, v);
        nv = asg_arc_n(sg, v);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            u = av[i].v>>1;
            if(x->rGraphVis.a[u]) continue;
            x->rGraphVis.a[u] = 1;
            if(IF_HOM(u, *bub))
            {
                if(d < k) kdq_push(uint64_t, x->q, set_dv(u, d+1));
            }
            else
            {
                kdq_push(uint64_t, x->q, set_dv(u , d));
                b_mer_d = d;
                if(IF_BUB(u, *bub) && x->bmerVis.a[u] == 0)
                {
                    get_bubbles((bubble_type*)bub, bub->index[u], &beg, &sink, &a, &n, NULL);
                    for (j = 0; j < n; j++) x->bmerVis.a[(a[j]>>1)] = 1;
                }
                //must be here
                if(bub_only == 0) x->bmerVis.a[u] = 1;
            }            
        }


        v = v + 1;
        av = asg_arc_a(sg, v);
        nv = asg_arc_n(sg, v);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            u = av[i].v>>1;
            if(x->rGraphVis.a[u]) continue;
            x->rGraphVis.a[u] = 1;
            if(IF_HOM(u, *bub))
            {
                if(d < k) kdq_push(uint64_t, x->q, set_dv(u, d+1));
            }
            else
            {
                kdq_push(uint64_t, x->q, set_dv(u , d));
                b_mer_d = d;
                if(IF_BUB(u, *bub) && x->bmerVis.a[u] == 0)
                {
                    get_bubbles((bubble_type*)bub, bub->index[u], &beg, &sink, &a, &n, NULL);
                    for (j = 0; j < n; j++) x->bmerVis.a[(a[j]>>1)] = 1;
                }
                //must be here
                if(bub_only == 0) x->bmerVis.a[u] = 1;
            }            
        }
    }

    return b_mer_d;
}

void select_bmer_distance(uint32_t src, uint64_t k, const bubble_type* bub, min_cut_t* x, 
uint32_t bub_only, uint32_t bub_extend)
{
    uint32_t beg, sink, n, *a;
    asg_t *sg = bub->ug->g;
    uint64_t v, u, i, j, nv, w, first = 1;
    asg_arc_t *av = NULL;
    reset_pdq(&(x->pq));

    x->bmerVis.a[src>>1] = 1;
    x->pq.dis.a[src] = 0;
    push_pdq(&(x->pq), src, 0);

    while (pdq_cnt(x->pq) > 0)
    {
        pop_pdq(&(x->pq), &v, &w);
        x->pq.vis.a[v] = 1;
        if(x->pq.dis.a[v] > k) break;

        ///fprintf(stderr, "******utg%.6dl, dis: %lu\n", (int)((v>>1)+1), x->pq.dis.a[v]); 

        if(IF_BUB(v>>1, *bub))
        {
            if(bub_extend && x->bmerVis.a[v>>1] == 0)
            {
                get_bubbles((bubble_type*)bub, bub->index[v>>1], &beg, &sink, &a, &n, NULL);
                for (j = 0; j < n; j++) x->bmerVis.a[(a[j]>>1)] = 1;
            }
            x->bmerVis.a[v>>1] = 1;
        }
         
        if(IF_HET(v>>1, *bub) && bub_only == 0) x->bmerVis.a[v>>1] = 1;

        av = asg_arc_a(sg, v);
        nv = asg_arc_n(sg, v);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            u = av[i].v;
            w = (uint32_t)av[i].ul;
            if(first) w = 0;

            if(x->pq.vis.a[u] == 0 && x->pq.dis.a[u] > x->pq.dis.a[v] + w)
            {
                x->pq.dis.a[u] = x->pq.dis.a[v] + w;
                push_pdq(&(x->pq), u, x->pq.dis.a[u]);
            }
        }

        first = 0;
    }
}

void get_bmer_unitgs(min_cut_t* x, const bubble_type* bub, uint64_t k, uint64_t src)
{
    uint32_t beg, sink, n, *a;
    if(!IF_BUB(src, *bub)) return;
    get_bubbles((bubble_type*)bub, bub->index[src], &beg, &sink, &a, &n, NULL);
    memset(x->bmerVis.a, 0, x->bmerVis.n);
    ///select_bmer(src, k, bub, x, 1);
    select_bmer_distance(beg^1, k, bub, x, 1, 1);
    select_bmer_distance(sink^1, k, bub, x, 1, 1);
}

min_cut_t* clean_hap(hc_links* link, bubble_type* bub, const ma_ug_t *ug)
{
    double index_time = yak_realtime();
    min_cut_t* x; CALLOC(x, 1);
    kvec_t_u64_warp buff;
    kv_init(buff.a);
    init_min_cut_t(x, link, (const bubble_type*)bub, ug);

    // get_Connected_Components(&x);
    // check_connective(&x, link);
    // print_rGraph(&x);

    long long i;
    uint64_t k, uID;
    
    ///for (i = 0; (uint64_t)i < x.n; i++)
    for (i = x->n - 1; i >= 0; i--)
    {
        uID = x->order.a[i] & x->uID_mode;
        ///fprintf(stderr, "uID: %lu, f.n: %lu\n", uID, (uint64_t)link->a.a[uID].f.n);
        if(link->a.a[uID].f.n == 0) continue;
        for (k = 0; k < link->a.a[uID].f.n; k++)
        {
            if(link->a.a[uID].f.a[k].del) continue;
            if(x->utgVis.a[link->a.a[uID].f.a[k].uID] == 0) break;
        }
        ///fprintf(stderr, "k: %lu\n", k);
        if(k == link->a.a[uID].f.n) continue;
        reset_min_cut_t(x, link);
        ///fprintf(stderr, "reset\n");
        get_s_t(x, link, uID, &(x->src), &(x->dest), &buff);
        ///fprintf(stderr, "x.src: %lu, x.dest: %lu\n", x.src, x.dest);
        ///Note: should only consider edges betweem bubbles, ignore edges to homo untigs

        /*******************************for debug************************************/
        ///if(!select_large_node(ug, &x, x.src, x.dest, 10, 0)) continue; 
        ///if(uID != 26818) continue;
        //if(uID != 45179) continue;
        ///memset(x.bmerVis.a, 1, x.bmerVis.n);
        get_bmer_unitgs(x, bub, x->b_mer, x->src);
        x->bmerVis.a[x->src] = x->bmerVis.a[x->dest] = 1;
        /*******************************for debug************************************/

        graph_cut(x->src, x->dest, x);
        ///fprintf(stderr, "graph_cut\n");
        x->utgVis.a[uID] = 1;
        ///exit(0);
    }

    reset_min_cut_t(x, link);

    fprintf(stderr, "[M::%s::%.3f] # edges: %lu, # cutted edges: %lu\n", 
                                    __func__, yak_realtime()-index_time, x->n_e, x->c_e);
    update_link_by_min_cut_t(x, link);
    ///destory_min_cut_t(x);
    kv_destroy(buff.a);
    return x;
}

void init_G_partition(G_partition* x, uint64_t n_utg)
{
    uint64_t i;
    kv_init(*x);
    MALLOC(x->index, n_utg);
    for (i = 0; i < n_utg; i++)
    {
        x->index[i] = (uint32_t)-1;
    }
}

void reset_G_partition(G_partition* x, uint64_t n_utg)
{
    uint64_t i;
    x->n = 0;
    for (i = 0; i < n_utg; i++)
    {
        x->index[i] = (uint32_t)-1;
    }
}

void destory_G_partition(G_partition* x)
{
    uint64_t i;
    for (i = 0; i < x->n; i++)
    {
        kv_destroy(x->a[i].a);
    }
    kv_destroy(*x);
    free(x->index);
}

double get_hc_weight(uint32_t query, uint32_t v0, uint32_t root, bub_p_t_warp *b, min_cut_t* x)
{
    if(v0 == root) return 0;
    uint32_t v, u;
    hc_edge *p = NULL;
    double weight = 0;
    v = v0;
    do {
		u = b->a[v].p; // u->v
        p = get_rGraph_edge(x, query>>1, v>>1);
        if(p) weight += p->weight;
		v = u;
	} while (v != root);

    return weight;
}

void set_path(bub_p_t_warp *b, uint32_t root, uint8_t* flag, uint8_t label)
{
    uint32_t v, u;
    ///v is the sink of this bubble
	v = b->S.a[0];
    do {
		u = b->a[v].p; // u->v
        flag[v>>1] |= label;
		v = u;
	} while (v != root);
    flag[b->S.a[0]>>1] = 0;
}

uint64_t trace_phase_path(ma_ug_t *ug, uint32_t s, uint32_t d, bub_p_t_warp *b, min_cut_t* x, uint8_t* flag, uint8_t label)
{
    asg_t *g = ug->g;
    if(g->seq[s>>1].del) return 0; // already deleted
    if(get_real_length(g, s, NULL)<2) return 0;
    uint32_t i, n_pending, is_first, to_replace, cur_nc, cur_uc, cur_ac, n_tips, tip_end, n_pop;
    double cur_nh, cur_rate, max_rate;
    ///S saves nodes with all incoming edges visited
    b->S.n = b->T.n = b->b.n = b->e.n = 0;
    ///for each node, b->a saves all related information
    b->a[s].d = b->a[s].nc = b->a[s].ac = b->a[s].uc = 0; b->a[s].nh = 0;
    ///b->S is the nodes with all incoming edges visited
    kv_push(uint32_t, b->S, s);
    n_pop = n_tips = n_pending = 0;
    tip_end = (uint32_t)-1;
    is_first = 1;

    do {
        ///v is a node that all incoming edges have been visited
        ///d is the distance from v0 to v
        uint32_t v = kv_pop(b->S);
        uint32_t d = b->a[v].d, nc = b->a[v].nc, uc = b->a[v].uc, ac = b->a[v].ac; 
        double nh = b->a[v].nh;

        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        for (i = 0; i < nv; ++i) {
            uint32_t w = av[i].v, l = (uint32_t)av[i].ul; // v->w with length l, not overlap length
            bub_p_t *t = &b->a[w];
            //got a circle
            if ((w>>1) == (s>>1)) goto pop_reset;
            //important when poping at long untig graph
            if(is_first) l = 0;
            if (av[i].del) continue;
            ///push the edge
            kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);

            if (t->s == 0) 
            { // this vertex has never been visited
                kv_push(uint32_t, b->b, w); // save it for revert
                ///t->p is the parent node of 
                ///t->s = 1 means w has been visited
                ///d is len(v0->v), l is len(v->w), so t->d is len(v0->w)
                t->p = v, t->s = 1, t->d = d + l, t->nc = nc + ug->u.a[(w>>1)].n;
                t->r = get_real_length(g, w^1, NULL);
                /**need fix**/
                t->nh = nh + get_hc_weight(w, v, s, b, x);
                t->ac = ac + (flag[(w>>1)] == 0? ug->u.a[(w>>1)].n : 0);
                t->uc = uc + (flag[(w>>1)] != 0? ug->u.a[(w>>1)].n : 0);
                ++n_pending;
            }
            else {
                to_replace = 0;

                cur_nc = nc + ug->u.a[(w>>1)].n;
                /**need fix**/
                cur_nh = nh + get_hc_weight(w, v, s, b, x);
                cur_ac = ac + (flag[(w>>1)] == 0? ug->u.a[(w>>1)].n : 0);
                cur_uc = uc + (flag[(w>>1)] != 0? ug->u.a[(w>>1)].n : 0);
                cur_rate = ((double)(cur_ac)/(double)(cur_ac+cur_uc));
                max_rate = ((double)(t->ac)/(double)(t->ac+t->uc));

                if(cur_rate > max_rate)
                {
                    to_replace = 1;
                }
                else if(cur_rate == max_rate)
                {
                    if(cur_nh > t->nh)
                    {
                        to_replace = 1;
                    }
                    else if(cur_nh == t->nh)
                    {
                        if(cur_nc > t->nc)
                        {
                            to_replace = 1;
                        }
                        else if(cur_nc == t->nc)
                        {
                            if(d + l > t->d)
                            {
                                to_replace = 1;
                            }
                        }
                    }
                }


                if(to_replace)
                {
                    t->p = v;
                    t->nc = cur_nc;
                    t->nh = cur_nh;
                    t->ac = cur_ac;
                    t->uc = cur_uc;
                }


                if (d + l < t->d) t->d = d + l; // update dist
            }

            if (--(t->r) == 0) {
                uint32_t x = get_real_length(g, w, NULL);
                if(x > 0)
                {
                    kv_push(uint32_t, b->S, w);
                }
                else
                {
                    ///at most one tip
                    if(n_tips != 0) goto pop_reset;
                    n_tips++;
                    tip_end = w;
                }
                --n_pending;
            }
        }
        is_first = 0;


        if(n_tips == 1)
        {
            if(tip_end != (uint32_t)-1 && n_pending == 0 && b->S.n == 0)
            {
                ///sink is b.S.a[0]
                kv_push(uint32_t, b->S, tip_end);
                break;
            }
            else
            {
                goto pop_reset;
            }
        }

        if (i < nv || b->S.n == 0) goto pop_reset;
    }while (b->S.n > 1 || n_pending);


    n_pop = 1;
    /**need fix**/
    set_path(b, s, flag, label);
    pop_reset:

    for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
        bub_p_t *t = &b->a[b->b.a[i]];
        t->p = t->d = t->nc = t->ac = t->uc = t->r = t->s = 0;
        t->nh = 0;
    }

    return n_pop;
}

inline void get_phased_block(G_partition* x, bubble_type* bub, uint64_t id, 
uint32_t* beg, uint32_t* sink, uint32_t** h0, uint32_t* h0_n, uint32_t** h1, uint32_t* h1_n,
uint32_t* phased, uint32_t* bub_id)
{
    if(bub && beg && sink && bub_id)
    {
        (*bub_id) = (*beg) = (*sink) = (uint32_t)-1;
        if(x->a[id].a.n > 0)
        {
            (*bub_id) = bub->index[x->a[id].a.a[0]];
            if(IF_BUB(x->a[id].a.a[0], *bub))
            {
                (*beg) = bub->list.a[bub->num.a[(*bub_id)]];
                (*sink) = bub->list.a[bub->num.a[(*bub_id)] + 1];
            }
        }
    }
    
    (*h0) = x->a[id].a.a; 
    (*h0_n) = x->a[id].h[0];

    (*h1) = x->a[id].a.a + x->a[id].h[0];
    (*h1_n) = x->a[id].h[1];
    if(phased) (*phased) = x->a[id].full_bub;
    if((*h0_n) == 0) (*h0) = NULL;
    if((*h1_n) == 0) (*h1) = NULL;
}

double get_co_weight(uint32_t *query, uint32_t query_n, uint32_t *target, uint32_t target_n, min_cut_t* m)
{
    double weight = 0;
    hc_edge *p = NULL;
    uint32_t i, k;
    for (i = 0; i < query_n; i++)
    {
        for (k = 0; k < target_n; k++)
        {
            p = get_rGraph_edge(m, query[i], target[k]);
            if(p) weight += p->weight;
        }
    }

    return weight;
}

void phase_bubble(uint64_t bid, bub_p_t_warp *b, bubble_type* bub, uint8_t* flag, const ma_ug_t *ug, 
min_cut_t* m, hc_links* link, G_partition* x)
{
    #define HAP1_LAB 1
    #define HAP2_LAB 2

    partition_warp* res = NULL;
    kv_pushp(partition_warp, *x, &res);
    memset(flag, 0, ug->g->n_seq);
    uint32_t beg, sink, n, *a, i, k;
    get_bubbles(bub, bid, &beg, &sink, &a, &n, NULL);
    res->full_bub = 0;
    trace_phase_path((ma_ug_t *)ug, beg, sink, b, m, flag, HAP1_LAB);
    trace_phase_path((ma_ug_t *)ug, beg, sink, b, m, flag, HAP2_LAB);
    kv_init(res->a);
    for (i = 0; i < ug->g->n_seq; i++)
    {
        if(flag[i] & (uint8_t)HAP1_LAB) kv_push(uint32_t, res->a, i);
    }
    res->h[0] = res->a.n;
    for (i = 0; i < ug->g->n_seq; i++)
    {
        if(flag[i] & (uint8_t)HAP2_LAB) kv_push(uint32_t, res->a, i);
    }
    res->h[1] = res->a.n - res->h[0];
    if(n == 2) res->full_bub = 1;
    if(res->full_bub == 0)
    {
        double self = 0, intersec = 0;
        uint32_t *h0 = NULL, *h1 = NULL;
        h0 = res->a.a; h1 = res->a.a + res->h[0];
        self += get_co_weight(h0, res->h[0], h0, res->h[0], m);
        self += get_co_weight(h1, res->h[1], h1, res->h[1], m);

        intersec += get_co_weight(h0, res->h[0], h1, res->h[1], m);
        intersec = intersec * 2;

        if(self > intersec) res->full_bub = 1;
    }

    if(res->full_bub == 0)
    {
        res->a.n = 0;
        uint32_t v, u = 0, uv, k_n, pre_n = x->n;
        hc_linkeage* t = NULL;
        x->n--;
        for (i = 0; i < n; i++)
        {
            v = a[i]>>1;
            t = &(link->a.a[v]);
            for (k = k_n = 0; k < t->f.n; k++)
            {
                if(t->f.a[k].del) continue;
                k_n++;
                u = t->f.a[k].uID;
            }
            if(k_n != 1) continue;

            t = &(link->a.a[u]);
            for (k = k_n = 0; k < t->f.n; k++)
            {
                if(t->f.a[k].del) continue;
                k_n++;
                uv = t->f.a[k].uID;
            }
            if(k_n != 1) continue;
            if(uv != v) continue;
            
            ///avoid dups
            for (k = 0; k < i; k++)
            {
                if((a[k]>>1) == u) break;
            }
            if(k < i) continue;
            

            kv_pushp(partition_warp, *x, &res);
            if(x->n > pre_n) kv_init(res->a);
            res->full_bub = 0;
            res->h[0] = res->h[1] = 1;
            kv_push(uint32_t, res->a, v);
            kv_push(uint32_t, res->a, u);

            for (k = 0; k < res->a.n; k++)
            {
                x->index[res->a.a[k]] = x->n-1;
            }
        }
    }
    else
    {
        for (k = 0; k < res->a.n; k++)
        {
            x->index[res->a.a[k]] = x->n-1;
        }
    }
    
    /*******************************for debug************************************/
    // for (i = 0; i < ug->g->n_seq; i++)
    // {
    //     if(flag[i] & (uint8_t)3)
    //     {
    //         uint32_t k;
    //         for (k = 0; k < n; k++)
    //         {
    //             if((a[k]>>1) == i)
    //             {
    //                 break;
    //             }
    //         }

    //         if(k == n) fprintf(stderr, "ERROR5\n");
    //     }
    // }
    /*******************************for debug************************************/
}

void print_phased_bubble(G_partition* x, bubble_type* bub, uint32_t utg_n)
{
    uint64_t i, k;
    uint32_t beg = 0, sink = 0, h0_n, h1_n, *h0, *h1, full_bub = 0, bubID = 0;

    for (i = 0; i < x->n; i++)
    {
        get_phased_block(x, bub, i, &beg, &sink, &h0, &h0_n, &h1, &h1_n, &full_bub, &bubID);

        fprintf(stderr, "\n[%lu]\tbeg:utg%.6ul\tsink:utg%.6ul\tphased=%u\n", i, (beg>>1)+1, (sink>>1)+1, full_bub);
        for (k = 0; k < h0_n; k++)
        {
            fprintf(stderr, "(0) utg%.6ul\n", h0[k] + 1);
        }

        for (k = 0; k < h1_n; k++)
        {
            fprintf(stderr, "(1) utg%.6ul\n", h1[k] + 1);
        }

        uint32_t n, *a;
        get_bubbles(bub, bubID, &beg, &sink, &a, &n, NULL);
        if(n > 2) fprintf(stderr, "complex\n");
    }

    /*******************************for debug************************************/
    for (i = 0; i < utg_n; i++)
    {
        if(x->index[i] == (uint32_t)-1) continue;
        partition_warp* p = &(x->a[x->index[i]]);
        for (k = 0; k < p->a.n; k++)
        {
            if(p->a.a[k] != i) break;
        }

        if(k == p->a.n) fprintf(stderr, "ERROR\n");
    }
    /*******************************for debug************************************/
}

G_partition* clean_bubbles(hc_links* link, bubble_type* bub, min_cut_t* m, const ma_ug_t *ug)
{
    double index_time = yak_realtime();
    uint64_t i;
    bub_p_t_warp b; 
    memset(&b, 0, sizeof(bub_p_t_warp));
    CALLOC(b.a, ug->g->n_seq*2); 
    uint8_t* flag = NULL;
    CALLOC(flag, ug->g->n_seq);
    G_partition* x; CALLOC(x, 1);
    init_G_partition(x, ug->g->n_seq);

    for (i = 0; i < bub->f_bub; i++)
    {
        phase_bubble(i, &b, bub, flag, ug, m, link, x);
    }

    free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
    free(flag);
    fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
    ///print_phased_bubble(x, bub, ug->g->n_seq);

    return x;
}

uint64_t get_hic_distance(pe_hit* hit, hc_links* link, const ha_ug_index* idx, uint32_t *is_cc)
{
    uint64_t s_uid, s_dir, e_uid, e_dir, u_dis, k;
    long long s_pos, e_pos;
    s_uid = ((hit->s<<1)>>(64 - idx->uID_bits)); s_pos = hit->s & idx->pos_mode;
    e_uid = ((hit->e<<1)>>(64 - idx->uID_bits)); e_pos = hit->e & idx->pos_mode;
    if(s_uid == e_uid)
    {
        if(is_cc) (*is_cc) = 1;
        return MAX(s_pos, e_pos) - MIN(s_pos, e_pos);
    } 
    hc_linkeage* t = &(link->a.a[s_uid]);
    for (k = 0; k < t->e.n; k++)
    {
        if(t->e.a[k].del || t->e.a[k].uID != e_uid) continue;
        if(is_cc) (*is_cc) = t->e.a[k].is_cc;
        s_dir = (!!(t->e.a[k].dis&(uint64_t)2));
        e_dir = (!!(t->e.a[k].dis&(uint64_t)1));
        u_dis = (t->e.a[k].dis ==(uint64_t)-1? (uint64_t)-1 : t->e.a[k].dis>>3);
        // if(s_uid == 24684 && s_pos == 124953 && e_uid == 16950 && e_pos == 93039)
        // {
        //     fprintf(stderr, "*****************s_dir: %lu, e_dir: %lu, u_dis: %lu\n", s_dir, e_dir, u_dis);
        // }
        if(u_dis == (uint64_t)-1) return (uint64_t)-1;
        if(s_dir == 1) s_pos = (long long)idx->ug->g->seq[s_uid].len - s_pos - 1;
        if(e_dir == 1) e_pos = (long long)idx->ug->g->seq[e_uid].len - e_pos - 1;
        e_pos = e_pos + u_dis - (long long)idx->ug->g->seq[e_uid].len;
        return MAX(s_pos, e_pos) - MIN(s_pos, e_pos);
    }

    return (uint64_t)-1;
}

hc_edge* get_hc_edge(hc_links* link, uint64_t src, uint64_t dest, uint64_t dir)
{
    if(src >= link->a.n) return NULL;
    uint64_t i, n;
    hc_edge* a = NULL;
    if(dir == 0)
    {
        n = link->a.a[src].e.n;
        a = link->a.a[src].e.a;
    }
    else
    {
        n = link->a.a[src].f.n;
        a = link->a.a[src].f.a;
    }
    
    for (i = 0; i < n; i++)
    {
        if(a[i].del) continue;
        if(a[i].uID == dest) return &(a[i]);
    }

    return NULL;
}

inline double get_trans(const ha_ug_index* idx, uint64_t x)
{
    return idx->a*(x/idx->frac) + idx->b;
}

inline double get_trans_weight_advance(const ha_ug_index* idx, uint64_t x, trans_idx* dis)
{
    long double rate = 0;

    ///if(x == (uint64_t)-1) x = dis->med;
    if(x != (uint64_t)-1)
    {
        if(x < dis->max)
        {
            uint64_t i;
            for (i = 0; i < dis->n; i++)
            {
                if(x < dis->a[i].end && x >= dis->a[i].beg) break;
            }
            if(i < dis->n)
            {
                rate = ((double)(dis->a[i].cnt_1))/((double)(dis->a[i].cnt_0 + dis->a[i].cnt_1));
            }
            else
            {
                rate = get_trans(idx, x);
            } 
        }
        else
        {
            rate = get_trans(idx, x);
        }
    }
    else
    {
        // rate = 0.2;
        rate = 0.4;
    }

    
    if(rate < 0) rate = 0;
    rate += OFFSET_RATE;
    if(rate > 0.5) rate = 0.5;
    rate -= OFFSET_SECOND_RATE; //[OFFSET_RATE - OFFSET_SECOND_RATE, 0.5 - OFFSET_SECOND_RATE]

    long double w = logl((1/rate)-1)*SCALL;
    if(w < OFFSET_RATE_MIN_W) w = OFFSET_RATE_MIN_W;
    if(w > OFFSET_RATE_MAX_W) w = OFFSET_RATE_MAX_W;
    return w;
}

void LeastSquare_advance(trans_idx* dis, ha_ug_index* idx, uint64_t med)
{
    #define SCAL_RATE 1000
    long double t1=0, t2=0, t3=0, t4=0, x, y;
    uint64_t i, m, ava_size;

    for (i = m = 0; i < dis->n; i++)
    {
        x = ((double)(dis->a[i].beg + dis->a[i].end))/2;
        y = ((double)(dis->a[i].cnt_1))/((double)(dis->a[i].cnt_0 + dis->a[i].cnt_1));
        if(dis->a[i].beg >= med) break;

        t1 += x*x;
        t2 += x;
        t3 += x*y;
        t4 += y;
        m++;
    }

    if(i < dis->n)
    {
        uint64_t beg, end, cnt_0, cnt_1;
        for (beg = dis->a[i].beg, end = dis->a[i].end, cnt_0 = cnt_1 = 0; i < dis->n; i++)
        {
            cnt_0 += dis->a[i].cnt_0;
            cnt_1 += dis->a[i].cnt_1;
            beg = MIN(beg, dis->a[i].beg);
            end = MAX(end, dis->a[i].end);
        }

        x = ((double)(beg + end))/2;
        y = ((double)(cnt_1))/((double)(cnt_0 + cnt_1));

        t1 += x*x;
        t2 += x;
        t3 += x*y;
        t4 += y;
        m++;
    }



    if(t2 > t4)
    {
        idx->frac = t2/t4;
        if(idx->frac > SCAL_RATE) idx->frac = idx->frac / SCAL_RATE;
        if(idx->frac < 1) idx->frac = 1;
    }

    ava_size = m;
    t1 /= (idx->frac*idx->frac);
    t2 /= idx->frac;
    t3 /= idx->frac;
    if((t1*ava_size - t2*t2) != 0)
    {
        idx->a = (t3*ava_size - t2*t4) / (t1*ava_size - t2*t2);
    }
    if((t1*ava_size - t2*t2) != 0)
    {
        idx->b = (t1*t4 - t2*t3) / (t1*ava_size - t2*t2);
    }
}

void weight_edges_advance(ha_ug_index* idx, kvec_pe_hit* hits, hc_links* link, bubble_type* bub, trans_idx* dis)
{
    uint64_t k, i, shif = 64 - idx->uID_bits, beg, end, t_d;
    hc_edge *e1 = NULL, *e2 = NULL;
    long double weight;

    for (i = 0; i < link->a.n; i++)
    {
        for (k = 0; k < link->a.a[i].e.n; k++)
        {
            if(link->a.a[i].e.a[k].del) continue;
            link->a.a[i].e.a[k].weight = 0;
        }
    }

    for (k = 0; k < hits->a.n; ++k) 
    {
        beg = ((hits->a.a[k].s<<1)>>shif);
        end = ((hits->a.a[k].e<<1)>>shif);

        if(beg == end) continue;
        if(IF_HOM(beg, *bub)) continue;
        if(IF_HOM(end, *bub)) continue;
        
        t_d = get_hic_distance(&(hits->a.a[k]), link, idx, NULL);
        if(t_d == (uint64_t)-1) continue;

        e1 = get_hc_edge(link, beg, end, 0);
        e2 = get_hc_edge(link, end, beg, 0);
        if(e1 == NULL || e2 == NULL) continue;
        weight = 1;
        if(dis)
        {
            weight = get_trans_weight_advance(idx, t_d, dis);
        }
        
        e1->weight += weight; e1->occ++;
        e2->weight += weight; e2->occ++;
    }
}

void get_bub_id(bubble_type* bub, uint32_t root, uint64_t* id0, uint64_t* id1, uint32_t check_het)
{
    if(id0) (*id0) = (uint64_t)-1;
    if(id1) (*id1) = (uint64_t)-1;
    uint64_t b_id0 = (uint64_t)-1, b_id1 = (uint64_t)-1;
    uint32_t beg, sink;

    if((bub->b_s_idx.a[root]&0xffffffff) != 0xffffffff)
    {
        b_id0 = bub->b_s_idx.a[root]&0xffffffff;
        if(check_het)
        {
            get_bubbles(bub, b_id0, &beg, &sink, NULL, NULL, NULL);
            if(IF_HET(beg>>1, *bub) && IF_HET(sink>>1, *bub)) b_id0 = (uint64_t)-1;
        } 
    }

    if((bub->b_s_idx.a[root]&0xffffffff00000000) != 0xffffffff00000000)
    {
        b_id1 = bub->b_s_idx.a[root]&0xffffffff00000000; b_id1 >>= 32;
        if(check_het)
        {
            get_bubbles(bub, b_id1, &beg, &sink, NULL, NULL, NULL);
            if(IF_HET(beg>>1, *bub) && IF_HET(sink>>1, *bub)) b_id1 = (uint64_t)-1;
        }        
    }

    if(b_id0 == (uint64_t)-1 && b_id1 != (uint64_t)-1)
    {
        b_id0 = b_id1;
        b_id1 = (uint64_t)-1;
    }

    if(id0) (*id0) = b_id0;
    if(id1) (*id1) = b_id1;
}

///return how many bubbles linked by this node
uint32_t connect_bub_occ(bubble_type* bub, uint32_t root_id, uint32_t check_het)
{
    uint64_t id0, id1, occ = 2;  
    get_bub_id(bub, root_id, &id0, &id1, check_het);
    if(id0 == (uint64_t)-1) occ--;
    if(id1 == (uint64_t)-1) occ--;
    return occ;
}
///x_0 and x_1 are the ids of unitigs; 
///x_0_b_id and x_1_b_id are the ids of bubble graph; 
int ma_2_bub_arc(bubble_type* bub, uint32_t x_0, uint32_t* x_0_b_id, uint32_t x_1, uint32_t* x_1_b_id,
asg_arc_t *p, uint32_t check_het)
{
    uint64_t id0, ori_0, id1, ori_1, tmp_id;
    uint32_t beg, sink, n, *a, x;
    uint32_t beg_0, sink_0, beg_1, sink_1;
    if(x_0_b_id) id0 = (*x_0_b_id);
    if(x_1_b_id) id1 = (*x_1_b_id);
    if((x_0 != (uint32_t)-1) && (x_1 != (uint32_t)-1))
    {
        if(((x_0>>1) == (x_1>>1)))
        {
            if(x_0_b_id == NULL && x_1_b_id == NULL)
            {
                get_bub_id(bub, x_0>>1, &id0, &id1, check_het);
            }
            
            get_bubbles(bub, id0, &beg_0, &sink_0, &a, &n, NULL);
            get_bubbles(bub, id1, &beg_1, &sink_1, &a, &n, NULL);


            ori_0 = (uint64_t)-1;
            if(x_0 == (beg_0^1))
            {
                ori_0 = 1;
            }
            else if(x_0 == (sink_0^1))
            {
                ori_0 = 0;
            }
            else if(x_0 == (beg_1^1))
            {
                ori_0 = 1+2;
            }
            else if(x_0 == (sink_1^1))
            {
                ori_0 = 0+2;
            }
            else
            {
                fprintf(stderr, "error 0\n");
                return 0;
            }


            ori_1 = (uint64_t)-1;
            if(x_1 == (beg_0^1))
            {
                ori_1 = 1;
            }
            else if(x_1 == (sink_0^1))
            {
                ori_1 = 0;
            }
            else if(x_1 == (beg_1^1))
            {
                ori_1 = 1 + 2;
            }
            else if(x_1 == (sink_1^1))
            {
                ori_1 = 0 + 2;
            }
            else
            {
                fprintf(stderr, "error 1\n");
                return 0;
            }

            

            if((((ori_0>>1)^(ori_1>>1))&1) != 1)
            {
                fprintf(stderr, "error 10\n");
                fprintf(stderr, "x_0: %u, id0: %lu, beg_0: %u, sink_0: %u, ori_0: %lu\n", 
                                                            x_0, id0, beg_0, sink_0, ori_0);
                fprintf(stderr, "x_1: %u, id1: %lu, beg_1: %u, sink_1: %u, ori_1: %lu\n", 
                                                            x_1, id1, beg_1, sink_1, ori_1);
                return 0;
            }
             
            if(ori_0 & 2)
            {
                tmp_id = id0; id0 = id1; id1 = tmp_id;
            }
            
            ori_0 &= 1; ori_1 &= 1; ori_1 ^= 1;


            
            p->ul = (id0<<1) | ori_0; p->ul <<= 32; p->ul += 0;
            p->v = (id1<<1) | ori_1; 
            p->ol = 0; p->del = 0; p->el = p->no_l_indel = p->strong = 1;
        }
        else
        {
            if(x_0_b_id == NULL) get_bub_id(bub, x_0>>1, &id0, NULL, check_het);
            if(x_1_b_id == NULL) get_bub_id(bub, x_1>>1, &id1, NULL, check_het);


            get_bubbles(bub, id0, &beg, &sink, &a, &n, NULL);
            ori_0 = (uint64_t)-1;
            if(x_0 == (beg^1))
            {
                ori_0 = 1;
            }
            else if(x_0 == (sink^1))
            {
                ori_0 = 0;
            }
            else
            {
                fprintf(stderr, "error 0\n");
                return 0;
            }
            
            

            get_bubbles(bub, id1, &beg, &sink, &a, &n, NULL);
            ori_1 = (uint64_t)-1;
            if(x_1 == (beg^1))
            {
                ori_1 = 1;
            }
            else if(x_1 == (sink^1))
            {
                ori_1 = 0;
            }
            else
            {
                fprintf(stderr, "error 1\n");
                return 0;
            }
            ori_0 &= 1; ori_1 &= 1; ori_1 ^= 1;

            p->ul = (id0<<1) | ori_0; p->ul <<= 32; p->ul += 0;
            p->v = (id1<<1) | ori_1; 
            p->ol = 0; p->del = 0; p->el = p->no_l_indel = p->strong = 1;
        }
    }
    else
    {
        x = (uint32_t)-1;
        if(x_0 != (uint32_t)-1) x = x_0;
        if(x_1 != (uint32_t)-1) x = x_1;
        if(x == (uint32_t)-1) return 0;
        if(x_0_b_id == NULL && x_1_b_id == NULL)
        {
            get_bub_id(bub, x>>1, &id0, &id1, check_het);
        }
        
        if(id0 != (uint64_t)-1)
        {
            get_bubbles(bub, id0, &beg, &sink, &a, &n, NULL);
            if(x == (beg^1))
            {
                return 1;
            }
            else if(x == (sink^1))
            {
                return 1;
            }

            return 0;
        }
        
        if(id1 != (uint64_t)-1)
        {
            get_bubbles(bub, id1, &beg, &sink, &a, &n, NULL);
            if(x == (beg^1))
            {
                return 1;
            }
            else if(x == (sink^1))
            {
                return 1;
            }

            return 0;
        }
    }
    
    return 1;
}

#define arc_first(g, v) ((g)->arc[(g)->idx[(v)]>>32])
#define arc_cnt(g, v) ((uint32_t)(g)->idx[(v)])
void debug_bub_utg(bubble_type* bub, ma_ug_t *bug, asg_t *bsg, uint32_t check_het)
{
    uint32_t i, k, rId, rId_next, ori, ori_next, root, beg, end;
    uint64_t id0, id1;
    ma_utg_t *u = NULL;
    asg_arc_t *t = NULL;
    for (i = 0; i < bug->u.n; i++)
    {
        u = &(bug->u.a[i]);
        if(u->n == 0) continue;
        for (k = 0; k < u->n; k++)
        {
            if(k+1 >= u->n) continue;
            rId = u->a[k]>>33;
            ori = u->a[k]>>32&1;
            get_bubbles(bub, rId, ori == 1?&root:NULL, ori == 0?&root:NULL, NULL, NULL, NULL);

            t = &(arc_first(bsg, u->a[k]>>32));
            get_bub_id(bub, root>>1, &id0, &id1, check_het);
            if(id0 == (uint64_t)-1 || (t->el == 1 && id1 == (uint64_t)-1) || (t->el == 0 && id1 != (uint64_t)-1))
            {
                fprintf(stderr, "sbsbsb0sbsbsb-utg%.6d, check_het: %u\n", (int)((root>>1)+1), check_het);
                fprintf(stderr, "id0: %lu, id1: %lu, t->el: %u\n", id0, id1, t->el);
                continue;
            }

            ///fprintf(stderr, "aaaaaaaa10aaaaaaaa-utg%.6d\n", (int)((root>>1)+1));

            rId_next = u->a[k+1]>>33;
            ori_next = u->a[k+1]>>32&1;


            get_bubbles(bub, rId, &beg, &end, NULL, NULL, NULL);
            if(ori == 1)
            {
                if(root != beg) fprintf(stderr, "sbsbsb1sbsbsb, root: %u, beg: %u, end: %u\n", root, beg, end);
            }
            else
            {
                if(root != end) fprintf(stderr, "sbsbsb2sbsbsb, root: %u, beg: %u, end: %u\n", root, beg, end);
            }

            if(t->el == 1)
            {
                get_bubbles(bub, rId_next, &beg, &end, NULL, NULL, NULL);
                if(ori_next == 0)
                {
                    if(root != (beg^1)) fprintf(stderr, "sbsbsb3sbsbsb, root: %u, beg: %u, end: %u\n", root, beg, end);
                }
                else
                {
                    if(root != (end^1)) fprintf(stderr, "sbsbsb4sbsbsb, root: %u, beg: %u, end: %u\n", root, beg, end);
                }
            }
        }
    }

    fprintf(stderr, "[M::%s]\n", __func__);
}

///just change the hap status of beg/sink, but they are are still at a chain of bubble
///might be ok 
inline void set_bub_idx(bubble_type* bub, ma_utg_t *bu, asg_t *untig_sg, int beg_idx, int end_idx, 
uint32_t is_to_hom, uint32_t check_het)
{
    int k;
    uint32_t rId, ori, root;
    uint64_t id0, id1, len0, len1;
    for (k = beg_idx; k <= end_idx; k++)
    {
        rId = bu->a[k]>>33;
        ori = bu->a[k]>>32&1;
        get_bubbles(bub, rId, ori == 1?&root:NULL, ori == 0?&root:NULL, NULL, NULL, NULL);

        if(is_to_hom && IF_HOM(root>>1, *bub)) continue;
        if(!is_to_hom && IF_HET(root>>1, *bub)) continue;

        
        get_bub_id(bub, root>>1, &id0, &id1, check_het);
        if(id0 == (uint64_t)-1 || id1 == (uint64_t)-1) continue;
        get_bubbles(bub, id0, NULL, NULL, NULL, NULL, &len0);
        get_bubbles(bub, id1, NULL, NULL, NULL, NULL, &len1);

        if(is_to_hom)
        {
            if(untig_sg->seq[root>>1].len > (MIN(len0, len1)*3)) continue;
            bub->index[root>>1] = (uint32_t)-1;
        }
        else
        {
            bub->index[root>>1] = bub->f_bub+1;
        }
    
    }
}

void determine_bub_idx(bubble_type* bub, ma_utg_t *bu, asg_t *untig_sg, uint64_t pLen, 
uint64_t rLEN, uint64_t r_hetLen, int beg_idx, int end_idx, uint32_t check_het)
{
    if(beg_idx > end_idx) return;

    uint64_t r_homLen = rLEN - r_hetLen;
    ///pLen: total length, rLEN: beg/sink length
    if(pLen > 0 && rLEN > 0 && r_hetLen > 0 && rLEN < pLen*0.5 && r_hetLen < rLEN * 0.2) ///set het to hom
    {
        set_bub_idx(bub, bu, untig_sg, beg_idx, end_idx, 1, bub->check_het);
    }
    else if(pLen > 0 && rLEN > 0 && r_homLen > 0 && rLEN > pLen*0.9 && r_homLen < rLEN * 0.1) ///set hom to het
    {
        set_bub_idx(bub, bu, untig_sg, beg_idx, end_idx, 0, bub->check_het);
    }
}

void detect_bub_graph(bubble_type* bub, asg_t *untig_sg)
{
    asg_t *bg = bub->b_g;
    ma_ug_t *ug = NULL;
    ug = ma_ug_gen(bub->b_g);
    ///debug_bub_utg(bub, ug, bg, bub->check_het);
    uint32_t i, k, rId, ori, root, r_root;
    int beg_idx, end_idx;
    uint64_t pLen, rLEN, r_hetLen;
    ma_utg_t *u = NULL;
    asg_arc_t *t = NULL;
    for (i = 0; i < ug->u.n; i++)
    {
        u = &(ug->u.a[i]);
        if(u->n == 0) continue;
        for (k = pLen = rLEN = r_hetLen = beg_idx = 0, end_idx = -1; k < u->n; k++)
        {
            rId = u->a[k]>>33;
            ori = u->a[k]>>32&1;
            get_bubbles(bub, rId, ori == 1?&root:&r_root, ori == 0?&root:&r_root, NULL, NULL, NULL);
            
            t = NULL;
            if(k+1 < u->n) t = &(arc_first(bg, u->a[k]>>32));

            pLen += bg->seq[rId].len;///path length in bubble
            if(end_idx < beg_idx) ///first bubble
            {
                pLen += untig_sg->seq[r_root>>1].len;
                rLEN += untig_sg->seq[r_root>>1].len;
                if(IF_HET(r_root>>1, *bub)) r_hetLen += untig_sg->seq[r_root>>1].len;
            }

            if(t)
            {
                if(t->el == 0)
                {
                    if(end_idx >= beg_idx)
                    {
                        pLen += untig_sg->seq[root>>1].len;
                        rLEN += untig_sg->seq[root>>1].len;
                        if(IF_HET(root>>1, *bub)) r_hetLen += untig_sg->seq[root>>1].len;
                        determine_bub_idx(bub, u, untig_sg, pLen, rLEN, r_hetLen, beg_idx, end_idx, bub->check_het);
                    } 
                    
                    pLen = rLEN = r_hetLen = 0;
                    beg_idx = k + 1; end_idx = k;
                }
                else
                {
                    pLen += t->ol;
                    rLEN += t->ol;
                    if(IF_HET(root>>1, *bub)) r_hetLen += t->ol;
                    end_idx = k;
                }
            }
        }
        
        if(end_idx >= beg_idx)
        {
            pLen += untig_sg->seq[root>>1].len;
            rLEN += untig_sg->seq[root>>1].len;
            if(IF_HET(root>>1, *bub)) r_hetLen += untig_sg->seq[root>>1].len;
            determine_bub_idx(bub, u, untig_sg, pLen, rLEN, r_hetLen, beg_idx, end_idx, bub->check_het);
        }
    }
    
    ma_ug_destroy(ug);
}

void get_bub_graph(ma_ug_t* ug, bubble_type* bub)
{
    asg_t *sg = ug->g;
    asg_arc_t t, *p = NULL;
    pdq pq;
    init_pdq(&pq, sg->n_seq<<1);
    uint32_t n_vtx = sg->n_seq<<1, v, k;
    uint32_t *pre = NULL; MALLOC(pre, n_vtx);
    uint32_t pre_id, adjecent, bub_occ;
    asg_t *bub_g = asg_init();
    for (v = 0; v < bub->f_bub; v++)
    {
        uint64_t pathbase;
        uint32_t beg, sink;
        get_bubbles(bub, v, &beg, &sink, NULL, NULL, &pathbase);
        asg_seq_set(bub_g, v, pathbase, (bub->check_het && IF_HET(beg>>1, *bub) && IF_HET(sink>>1, *bub))?1:0);
        bub_g->seq[v].c = PRIMARY_LABLE;
    }

    //check all unitigs
    for (v = 0; v < n_vtx; ++v)
    {
        if(sg->seq[v>>1].del) continue;
        if(bub->b_s_idx.a[v>>1] == (uint64_t)-1) continue; ///if (v>>1) is not a beg or sink of bubbles
        bub_occ = connect_bub_occ(bub, v>>1, bub->check_het);
        if(bub_occ == 0) continue;
        if(bub_occ == 2)
        {
            if(ma_2_bub_arc(bub, v, NULL, v^1, NULL, &t, bub->check_het))
            {
                t.ol = sg->seq[v>>1].len;
                p = asg_arc_pushp(bub_g);
                *p = t;
            }
            continue;
        }
        if(ma_2_bub_arc(bub, v, NULL, (uint32_t)-1, NULL, &t, bub->check_het) == 0) continue;

        get_shortest_path(v, &pq, sg, pre);
        for (k = 0; k < pq.dis.n; k++)
        {
            if(pq.dis.a[k] == (uint64_t)-1) continue;
            if(bub->b_s_idx.a[k>>1] == (uint64_t)-1) continue;
            if(connect_bub_occ(bub, k>>1, bub->check_het) == 0) continue;
            if((k>>1) == (v>>1)) continue;
            pre_id = pre[k];
            adjecent = 0;
            while (pre_id != v)
            {
                if(connect_bub_occ(bub, pre_id>>1, bub->check_het) > 0)
                {
                    adjecent = 1;
                    break;
                } 
                pre_id = pre[pre_id];
            }

            if(adjecent == 0)
            {
                if(ma_2_bub_arc(bub, v, NULL, k^1, NULL, &t, bub->check_het))
                {
                    t.el = 0; t.ol = pq.dis.a[k] + sg->seq[k>>1].len;
                    p = asg_arc_pushp(bub_g);
                    *p = t;
                }
            }
        }
    }

    free(pre);
    destory_pdq(&pq);

    asg_cleanup(bub_g);
    bub_g->r_seq = bub_g->n_seq;
    bub->b_g = bub_g;
}

void print_bubble_chain(bubble_type* bub, const char* command)
{
    ma_ug_t *ug = NULL;
    ug = ma_ug_gen(bub->b_g);
    uint32_t i, k, j, rId, beg, sink, *a, n;
    ma_utg_t *u = NULL;
    asg_arc_t *t = NULL;
    for (i = 0; i < ug->u.n; i++)
    {
        u = &(ug->u.a[i]);
        if(u->n == 0) continue;
        fprintf(stderr,"\n%s: chain-%u\n", command, i);
        for (k = 0; k < u->n; k++)
        {
            rId = u->a[k]>>33;
            get_bubbles(bub, rId, &beg, &sink, &a, &n, NULL);
            t = NULL;
            if(k+1 < u->n) t = &(arc_first(bub->b_g, u->a[k]>>32));
            fprintf(stderr, "[utg%.6dl, utg%.6dl] el=%u no_long_indel=%u, rId=%u, nv: %u, nv^: %u\n", 
                (int)((beg>>1)+1), (int)((sink>>1)+1), t?t->el:16, t?t->no_l_indel:16, rId, arc_cnt(bub->b_g, u->a[k]>>32), arc_cnt(bub->b_g, (u->a[k]>>32)^1));
            // if((u->a[k]>>33) == 9658)
            // {
            //     asg_arc_t *av;
            //     uint32_t nv, nv_i;
            //     av = asg_arc_a(bub->b_g, u->a[k]>>32);
            //     nv = asg_arc_n(bub->b_g, u->a[k]>>32);
            //     for (nv_i = 0; nv_i < nv; nv_i++)
            //     {
            //         if(av[nv_i].del) continue;
            //         fprintf(stderr, "v--->%u\n", av[nv_i].v>>1);
            //     }
                


            //     av = asg_arc_a(bub->b_g, (u->a[k]>>32)^1);
            //     nv = asg_arc_n(bub->b_g, (u->a[k]>>32)^1);
            //     for (nv_i = 0; nv_i < nv; nv_i++)
            //     {
            //         if(av[nv_i].del) continue;
            //         fprintf(stderr, "v^1--->%u\n", av[nv_i].v>>1);
            //     }

            // }
            if(bub->b_g->seq[rId].c == HAP_LABLE)
            {
                for (j = 0; j < n; j++)
                {
                    fprintf(stderr, ">>>utg%.6dl\n", (int)((a[j]>>1)+1));
                }
            }
        }
    }

    ma_ug_destroy(ug);
}

int is_simple_broken_bubble(ma_ug_t *unitig_ug, uint32_t x, uint32_t beg, uint32_t sink, uint32_t* new_het)
{
    uint32_t nv, v = (uint32_t)-1, u_s = (uint32_t)-1, u_e = (uint32_t)-1, i;
    asg_arc_t *av = NULL;
    (*new_het) = (uint32_t)-1;

    if((asg_arc_n(unitig_ug->g, x) == 1) 
                        && (asg_arc_n(unitig_ug->g, x^1) == 0))
    {
        v = x;
    }

    if((asg_arc_n(unitig_ug->g, x^1) == 1) 
                        && (asg_arc_n(unitig_ug->g, x) == 0))
    {
        v = x^1;
    }

    if(v == (uint32_t)-1) return 0;

    
    av = asg_arc_a(unitig_ug->g, v);
    nv = asg_arc_n(unitig_ug->g, v);
    for (i = 0; i < nv; i++)
    {
        if(av[i].del) continue;
        if((av[i].v>>1) == (beg>>1)) u_s = beg, u_e = sink;
        if((av[i].v>>1) == (sink>>1)) u_s = sink, u_e = beg;
    }

    if(u_s == (uint32_t)-1 || u_e == (uint32_t)-1) return 0;

    av = asg_arc_a(unitig_ug->g, u_s);
    nv = asg_arc_n(unitig_ug->g, u_s);
    if(nv != 2) return 0;
    for (i = 0; i < nv; i++)
    {
        if(av[i].del) continue;
        if(av[i].v == (v^1)) continue;
        if(av[i].v == (u_e^1))
        {
            (*new_het) = u_e;
            return 1;
        }
    }

    return 0;
}

int double_check_broken_bubble(asg_t *g, kvec_t_u32_warp* broken, uint32_t beg, uint32_t sink, 
uint8_t* vis_flag, uint32_t vis_flag_n, kvec_t_u32_warp* stack, asg_t *bsg, asg_arc_t *p_t)
{
    uint32_t cur, ncur, i, n, pre, occ;
    radix_sort_u32(broken->a.a, broken->a.a + broken->a.n);
    for (i = n = 0, pre = (uint32_t)-1; i < broken->a.n; i++)
    {
        if((broken->a.a[i]>>1) == (pre>>1)) continue;
        pre = broken->a.a[i];
        broken->a.a[n] = pre;
        n++;
    }
    broken->a.n = n;

    asg_arc_t *acur = NULL;
    memset(vis_flag, 0, vis_flag_n);
    stack->a.n = 0;
    kv_push(uint32_t, stack->a, beg);
    occ = 0;
    while (stack->a.n > 0)
    {
        stack->a.n--;
        cur = stack->a.a[stack->a.n];
        if(vis_flag[cur] == 0 && vis_flag[cur^1] == 0) occ++;
        if(vis_flag[cur] == 0 && cur != (beg^1) && (sink == (uint32_t)-1 || cur != (sink^1)))
        {
            vis_flag[cur] = 1;
            ncur = asg_arc_n(g, cur);
            acur = asg_arc_a(g, cur);
            for (i = 0; i < ncur; i++)
            {
                if(acur[i].del) continue;
                if(vis_flag[acur[i].v]) continue;
                kv_push(uint32_t, stack->a, acur[i].v);
            }
        }
        vis_flag[cur] = 1;
        

        cur^=1;
        if(vis_flag[cur] == 0 && cur != (beg^1) && (sink == (uint32_t)-1 || cur != (sink^1)))
        {
            vis_flag[cur] = 1;
            ncur = asg_arc_n(g, cur);
            acur = asg_arc_a(g, cur);
            for (i = 0; i < ncur; i++)
            {
                if(acur[i].del) continue;
                if(vis_flag[acur[i].v]) continue;
                kv_push(uint32_t, stack->a, acur[i].v);
            }
        }
        vis_flag[cur] = 1;
    }

    n = broken->a.n;
    if(beg != (uint32_t)-1) n++;
    if(sink != (uint32_t)-1) n++;

    if(occ > n)
    {
        ///fprintf(stderr, "\n+++++sb+++++beg-utg%.6ul, sink-utg%.6ul, occ: %u, n: %u\n", (beg>>1)+1, (sink>>1)+1, occ, n);
        if(bsg && p_t && beg != (uint32_t)-1 && sink != (uint32_t)-1)
        {
            p_t->del = 1;
            asg_arc_del(bsg, (p_t->v)^1, (p_t->ul>>32)^1, 1);
        }
        /*******************************for debug************************************/
        // memset(vis_flag, 0, vis_flag_n);
        // stack->a.n = 0;
        // kv_push(uint32_t, stack->a, beg);
        // occ = 0;
        // while (stack->a.n > 0)
        // {
        //     occ++;
        //     stack->a.n--;
        //     cur = stack->a.a[stack->a.n];

        //     fprintf(stderr, "cur-utg%.6ul\n", (cur>>1)+1);

        //     vis_flag[cur] = 1;
        //     if(cur == (beg^1) || cur == (sink^1)) continue;
        //     ncur = asg_arc_n(g, cur);
        //     acur = asg_arc_a(g, cur);
        //     for (i = 0; i < ncur; i++)
        //     {
        //         if(acur[i].del) continue;
        //         if(vis_flag[acur[i].v]) continue;
        //         kv_push(uint32_t, stack->a, acur[i].v);
        //     }

        //     cur^=1;
        //     if(vis_flag[cur]) continue;
        //     vis_flag[cur] = 1;
        //     if(cur == (beg^1) || cur == (sink^1)) continue;
        //     ncur = asg_arc_n(g, cur);
        //     acur = asg_arc_a(g, cur);
        //     for (i = 0; i < ncur; i++)
        //     {
        //         if(acur[i].del) continue;
        //         if(vis_flag[acur[i].v]) continue;
        //         kv_push(uint32_t, stack->a, acur[i].v);
        //     }
        // }

        // for (i = 0; i < broken->a.n; i++)
        // {
        //     fprintf(stderr, "*****cur-utg%.6ul\n", (broken->a.a[i]>>1)+1);
        // }
        /*******************************for debug************************************/
        return 0;
    }
     
    return 1;
}

int is_local_simple_circle(asg_t *g, uint32_t v)
{
    if(asg_arc_n(g, v) != asg_arc_n(g, v^1)) return 0;
    if(asg_arc_n(g, v) == 1) v = arc_first(g, v).v;
    if(asg_arc_n(g, v) != asg_arc_n(g, v^1)) return 0;
    if(asg_arc_n(g, v) != 2) return 0;
    uint32_t ncur, i, u;
    asg_arc_t *acur = NULL;
    ncur = asg_arc_n(g, v);
    acur = asg_arc_a(g, v);
    for (i = 0; i < ncur; i++)
    {
        if(acur[i].del) continue;
        u = acur[i].v;
        if(asg_arc_n(g, u) != 1 || asg_arc_n(g, u^1) != 1) continue;
        if(arc_first(g, u).v != v) continue;
        return 1;
    }
    return 0;
}

///actually not useful, and may have bug when one bubble at multipe chains
void update_bub_b_s_idx(bubble_type* bub) 
{
    memset(bub->b_s_idx.a, -1, bub->b_s_idx.n * sizeof(uint64_t));
    uint32_t i, v, beg, sink, n_bub = bub->num.n - 1;
    for (i = 0; i < n_bub; i++)
    {
        get_bubbles(bub, i, &beg, &sink, NULL, NULL, NULL);

        if(beg != (uint32_t)-1)
        {
            v = beg>>1;
            if(bub->b_s_idx.a[v] == (uint64_t)-1)
            {
                bub->b_s_idx.a[v] <<= 32;
                bub->b_s_idx.a[v] |= i;
            }
            else if((bub->b_s_idx.a[v] & 0xffffffff00000000) == 0xffffffff00000000)
            {
                bub->b_s_idx.a[v] <<= 32;
                bub->b_s_idx.a[v] |= i;
            }
        }

        

        if(sink != (uint32_t)-1)
        {
            v = sink>>1;
            if(bub->b_s_idx.a[v] == (uint64_t)-1)
            {
                bub->b_s_idx.a[v] <<= 32;
                bub->b_s_idx.a[v] |= i;
            }
            else if((bub->b_s_idx.a[v] & 0xffffffff00000000) == 0xffffffff00000000)
            {
                bub->b_s_idx.a[v] <<= 32;
                bub->b_s_idx.a[v] |= i;
            }
        }
    }
}

void update_bubble_graph(kvec_t_u32_warp* broken, uint32_t beg, uint32_t beg_bub_id, 
uint32_t sink, uint32_t sink_bub_id, bubble_type* bub, kvec_asg_arc_t_warp* edges, asg_t *bsg, 
asg_arc_t *p_t, uint8_t *bsg_idx, ma_ug_t *unitig_ug, uint64_t* occ_thres, uint64_t is_b_bub)
{
    uint32_t i, pre, n, bub_id, v;
    uint64_t occ;
    asg_arc_t t_f, t_r;
    
    radix_sort_u32(broken->a.a, broken->a.a + broken->a.n);
    for (i = n = occ = 0, pre = (uint32_t)-1; i < broken->a.n; i++)
    {
        if((broken->a.a[i]>>1) == (pre>>1)) continue;
        if(IF_HOM((broken->a.a[i]>>1), *bub))
        {
            if(is_local_simple_circle(unitig_ug->g, broken->a.a[i]))
            {
                bub->index[broken->a.a[i]>>1] = bub->f_bub+1;
            }
            else
            {
                continue;
            }
        }
         
        pre = broken->a.a[i];
        broken->a.a[n] = pre;
        occ += unitig_ug->u.a[broken->a.a[n]>>1].n;
        n++;
    }
    broken->a.n = n;
    ///if(broken->a.n == 0) return;
    if(broken->a.n == 1)
    {
        ///fprintf(stderr, "+++++sb+++++utg%.6ul\n", (broken->a.a[0]>>1)+1);
        if(beg != (uint32_t)-1 && sink != (uint32_t)-1 && 
                            is_simple_broken_bubble(unitig_ug, broken->a.a[0], beg, sink, &v))
        {
            if((v>>1) != (broken->a.a[0]>>1))
            {
                ///fprintf(stderr, "-----sb-----utg%.6ul\n", (v>>1)+1);
                kv_push(uint32_t, broken->a, v);
                occ += unitig_ug->u.a[v>>1].n;
                if(occ_thres && occ > (*occ_thres)) return;
                bub->index[v>>1] = bub->f_bub+1; ///set to het
            }
        }
    } 
    
    
    if(occ_thres && occ > (*occ_thres)) return;
    /********************push graph node********************/
    bub_id = bub->b_g->n_seq;
    asg_seq_set(bub->b_g, bub_id, 0, 0);
    bub->b_g->seq[bub_id].c = HAP_LABLE;
    if(is_b_bub) bub->b_bub++;
    /********************push graph node********************/

    /********************push bubble********************/
    kv_push(uint32_t, bub->num, bub->list.n);
    kv_push(uint64_t, bub->pathLen, 0);
    kv_push(uint32_t, bub->list, beg);
    kv_push(uint32_t, bub->list, sink);

    for (i = 0; i < broken->a.n; i++)
    {
        kv_push(uint32_t, bub->list, broken->a.a[i]);
        if(bsg_idx) bsg_idx[broken->a.a[i]>>1] = 1;
    }
    /********************push bubble********************/
    if(beg != (uint32_t)-1) ///beg_bub_id ----> bub_id
    {
        if(ma_2_bub_arc(bub, beg, &beg_bub_id, beg^1, &bub_id, &t_f, bub->check_het) && 
         ma_2_bub_arc(bub, beg^1, &bub_id, beg, &beg_bub_id, &t_r, bub->check_het))
         {
            t_f.el = 0; t_f.no_l_indel = 0; t_f.del = 0;
            kv_push(asg_arc_t, edges->a, t_f);

            t_r.el = 0; t_r.no_l_indel = 0; t_r.del = 0;
            kv_push(asg_arc_t, edges->a, t_r);
        }
    }

    if(sink != (uint32_t)-1) ///bub_id ----> sink_bub_id
    {
        if(ma_2_bub_arc(bub, sink^1, &bub_id, sink, &sink_bub_id, &t_f, bub->check_het) && 
         ma_2_bub_arc(bub, sink, &sink_bub_id, sink^1, &bub_id, &t_r, bub->check_het))
        {
            t_f.el = 0; t_f.no_l_indel = 0; t_f.del = 0;
            kv_push(asg_arc_t, edges->a, t_f);

            t_r.el = 0; t_r.no_l_indel = 0; t_r.del = 0;
            kv_push(asg_arc_t, edges->a, t_r);
        }
    }

    if(beg != (uint32_t)-1 && sink != (uint32_t)-1 && p_t)
    {
        p_t->del = 1;
        asg_arc_del(bsg, (p_t->v)^1, (p_t->ul>>32)^1, 1);
    }
}

void get_related_bub_nodes(kvec_t_u32_warp* broken, bubble_type* bub, pdq* pq, asg_t *unitig_g, 
                                        uint32_t *pre, uint32_t src, uint32_t dest, uint8_t *bsg_idx)
{
    uint32_t j_i, pre_id, adjecent;
    src ^= 1;
    get_shortest_path(src, pq, unitig_g, pre);
            
    for (j_i = 0; j_i < pq->dis.n; j_i++)
    {
        if(pq->dis.a[j_i] == (uint64_t)-1) continue;
        ///if(IF_HOM(j_i>>1, *bub)) continue;
        if((j_i>>1) == (src>>1)) continue;
        if((dest != (uint32_t)-1) && ((j_i>>1) == (dest>>1))) continue;

        pre_id = pre[j_i];
        adjecent = 0;
        while (pre_id != src)
        { 
            if(((dest != (uint32_t)-1) && ((pre_id>>1) == (dest>>1))) 
                                                || ((pre_id>>1) == (src>>1)))
            {
                adjecent = 1;
                break;
            }
            pre_id = pre[pre_id];
        }
        
        if(adjecent == 0)
        {
            if(broken->a.n == 0 || (broken->a.n > 0 && (j_i>>1) != (broken->a.a[broken->a.n - 1]>>1)))
            {
                if(bsg_idx && bsg_idx[(j_i>>1)])
                {
                    broken->a.n = 0;
                    return;
                }
                kv_push(uint32_t, broken->a, j_i);
            }
        }
    }
}

uint64_t calculate_chain_weight(ma_utg_t *u, bubble_type* bub, ma_ug_t *unitig_ug, chain_w_type* x)
{
    x->b_occ = x->g_occ = 0;
    uint32_t i, j, *a, n;
    uint64_t occ, occ_n, thres;
    for (i = occ = occ_n = 0; i < u->n; i++)
    {
        if(bub->b_g->seq[u->a[i]>>33].c != HAP_LABLE)
        {
            get_bubbles(bub, u->a[i]>>33, NULL, NULL, &a, &n, NULL);
            for (j = 0; j < n; j++)
            {
                occ += unitig_ug->u.a[a[j]>>1].n;
            }
            occ_n++;
        }
    }

    thres = (uint64_t)-1;
    if(occ_n > 0) thres = (occ*6)/occ_n;

    for (i = 0; i < u->n; i++)
    {
        occ = 0;
        get_bubbles(bub, u->a[i]>>33, NULL, NULL, &a, &n, NULL);
        for (j = 0; j < n; j++)
        {
            occ += unitig_ug->u.a[a[j]>>1].n;
        }

        if(bub->b_g->seq[u->a[i]>>33].c != HAP_LABLE || occ < thres)
        {
            x->g_occ += occ;
        }
        else
        {
            x->b_occ += occ;
        }
    }

    return thres;
}

int cmp_chain_weight(const void * a, const void * b)
{
    if((*(chain_w_type*)a).del != (*(chain_w_type*)b).del)
    {
        return (*(chain_w_type*)a).del > (*(chain_w_type*)b).del? 1 : -1;
    }
    else
    {
        long long a_occ = (*(chain_w_type*)a).g_occ - (*(chain_w_type*)a).b_occ;
        long long b_occ = (*(chain_w_type*)b).g_occ - (*(chain_w_type*)b).b_occ;
        if(a_occ != b_occ)
        {
            return a_occ > b_occ? -1 : 1;
        }
        else
        {
            return 0;
        }
    }
}

void resolve_bubble_chain_tangle_back(ma_ug_t* ug, bubble_type* bub, hc_links* link)
{
    ma_ug_t *copy_ug = copy_untig_graph(bub->b_ug);
    asg_arc_t *av = NULL;
    uint32_t i, j, k, v, w, w1, w2, nw1, nw2, nv, occ_e_1, occ_e_2, occ_c;
    ma_ug_t *bub_ug = copy_ug;
    ///ma_utg_t *u = NULL;
    buf_t b; memset(&b, 0, sizeof(buf_t));
    kvec_t_u32_warp stack, result;
    kv_init(stack.a); kv_init(result.a);
    uint8_t *vis = NULL; CALLOC(vis, ug->g->n_seq<<1);
    uint8_t *is_vis = NULL; CALLOC(is_vis, ug->g->n_seq<<1);
    kvec_t(uint64_t) occ_idx; kv_init(occ_idx); uint64_t tmp, *p = NULL;

    for (k = occ_idx.n = 0; k < bub_ug->g->n_seq; k++)
    {
        v = (k<<1);
        av = asg_arc_a(bub_ug->g, v);
        nv = asg_arc_n(bub_ug->g, v);
        for (i = 0, w = (uint32_t)-1, nw1 = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            nw1++;
            if((av[i].v>>1) == (v>>1)) continue;
            if(w != (uint32_t)-1) break;
            w = av[i].v;
        }
        if(i < nv) continue;
        w1 = w;


        v = (k<<1)+1;
        av = asg_arc_a(bub_ug->g, v);
        nv = asg_arc_n(bub_ug->g, v);
        for (i = 0, w = (uint32_t)-1, nw2 = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            nw2++;
            if((av[i].v>>1) == (v>>1)) continue;
            if(w != (uint32_t)-1) break;
            w = av[i].v;
        }
        if(i < nv) continue;
        w2 = w;

        if(nw1 <= 1 && nw2 <= 1) continue;
        if(w1 == (uint32_t)-1 && w2 == (uint32_t)-1) continue;
        if(w1 != (uint32_t)-1) w1 ^=1;
        if(w2 != (uint32_t)-1) w2 ^=1;

        if(w1 != (uint32_t)-1)
        {   
            w = (uint32_t)-1;
            if(w2 != (uint32_t)-1) w = w2^1;

            av = asg_arc_a(bub_ug->g, w1);
            nv = asg_arc_n(bub_ug->g, w1);
            for (i = 0; i < nv; i++)
            {
                if(av[i].del) continue;
                if((av[i].v>>1) == k) continue;
                if(av[i].v == w) continue;
                break;
            }
            if(i < nv) continue;
        }


        if(w2 != (uint32_t)-1)
        {   
            w = (uint32_t)-1;
            if(w1 != (uint32_t)-1) w = w1^1;

            av = asg_arc_a(bub_ug->g, w2);
            nv = asg_arc_n(bub_ug->g, w2);
            for (i = 0; i < nv; i++)
            {
                if(av[i].del) continue;
                if((av[i].v>>1) == k) continue;
                if(av[i].v == w) continue;
                break;
            }
            if(i < nv) continue;
        }

        
        if(w1 == (uint32_t)-1 && w2 != (uint32_t)-1) w1 = w2;
        if(w1 == w2) w2 = (uint32_t)-1;

        occ_c = occ_e_1 = occ_e_2 = (uint32_t)-1;
        set_b_utg_weight_flag(bub, &b, k<<1, NULL, 0, &occ_c);
        if(w1 != (uint32_t)-1) set_b_utg_weight_flag(bub, &b, w1^1, NULL, 0, &occ_e_1);
        if(w2 != (uint32_t)-1) set_b_utg_weight_flag(bub, &b, w2^1, NULL, 0, &occ_e_2);

        fprintf(stderr, "\n>>>>>>k=btg%.6ul (n=%u), w1=btg%.6ul (n=%u), w2=utg%.6ul (n=%u)\n", k+1, occ_c, 
                (w1>>1)+1, occ_e_1, (w2>>1)+1, occ_e_2);

        if(occ_c*5 >= occ_e_1) continue;
        if(occ_c*5 >= occ_e_2) continue;
        if(occ_c*10 >= (occ_e_1 + occ_e_2)) continue;
        kv_pushp(uint64_t, occ_idx, &p);
        (*p) = occ_e_1 + occ_e_2 - occ_c;
        (*p) <<= 32; (*p) += k;
        fprintf(stderr, "passed\n");
    }

    radix_sort_hc64(occ_idx.a, occ_idx.a + occ_idx.n);
    
    for (k = 0; k < occ_idx.n; ++k) 
    {
        tmp = occ_idx.a[k];
        occ_idx.a[k] = occ_idx.a[occ_idx.n - k - 1];
        occ_idx.a[occ_idx.n - k - 1] = tmp;
    }

    for (j = 0; j < occ_idx.n; j++)
    {
        k = (uint32_t)occ_idx.a[j];

        v = (k<<1);
        av = asg_arc_a(bub_ug->g, v);
        nv = asg_arc_n(bub_ug->g, v);
        for (i = 0, w = (uint32_t)-1, nw1 = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            nw1++;
            if((av[i].v>>1) == (v>>1)) continue;
            if(w != (uint32_t)-1) break;
            w = av[i].v;
        }
        if(i < nv) continue;
        w1 = w;


        v = (k<<1)+1;
        av = asg_arc_a(bub_ug->g, v);
        nv = asg_arc_n(bub_ug->g, v);
        for (i = 0, w = (uint32_t)-1, nw2 = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            nw2++;
            if((av[i].v>>1) == (v>>1)) continue;
            if(w != (uint32_t)-1) break;
            w = av[i].v;
        }
        if(i < nv) continue;
        w2 = w;

        if(nw1 <= 1 && nw2 <= 1) continue;
        if(w1 == (uint32_t)-1 && w2 == (uint32_t)-1) continue;
        if(w1 != (uint32_t)-1) w1 ^=1;
        if(w2 != (uint32_t)-1) w2 ^=1;



    }

    free(vis); free(is_vis); free(b.b.a); kv_destroy(occ_idx); kv_destroy(stack.a); kv_destroy(result.a);
    ma_ug_destroy(copy_ug);
}

uint32_t get_btg_occ(bubble_type* bub, uint32_t v)
{
    ma_ug_t *bub_ug = bub->b_ug;
    ma_utg_t *u = NULL;
    uint32_t k_i, k_j, *a = NULL, n, tan_occ = 0;

    u = &(bub_ug->u.a[v]);
    for (k_i = 0; k_i < u->n; k_i++)
    {
        get_bubbles(bub, u->a[k_i]>>33, NULL, NULL, &a, &n, NULL);
        for (k_j = 0; k_j < n; k_j++)
        {
            tan_occ += bub->ug->u.a[a[k_j]>>1].n;
        }
    }
    return tan_occ;
}

int check_bubble_tangle(bubble_type* bub, ma_ug_t* ug, uint32_t beg, uint32_t sink, 
double side_rate, double total_rate, uint32_t beg_occ, uint32_t sink_occ, 
uint8_t* is_vis, kvec_t_u32_warp* stack, kvec_t_u32_warp* res, uint8_t* chain_flag,
uint32_t* extra_check)
{
    if(extra_check) (*extra_check) = 1;
    uint32_t cur, tan_occ = 0, ncur, i, no_first = 0;
    asg_arc_t *acur = NULL;
    memset(is_vis, 0, ug->g->n_seq<<1);
    stack->a.n = 0;
    kv_push(uint32_t, stack->a, beg);
    if(res) res->a.n = 0;
    while (stack->a.n > 0)
    {
        stack->a.n--;
        cur = stack->a.a[stack->a.n];
        
        if(no_first && cur == beg) return 0;
        if(sink != (uint32_t)-1 && cur == sink) return 0;


        if(is_vis[cur] == 0 && is_vis[cur^1] == 0)
        {
            if((cur>>1) != (beg>>1) && (sink == (uint32_t)-1 || (cur>>1) != (sink>>1)))
            {
                if(res) kv_push(uint32_t, res->a, cur);
                if(chain_flag && chain_flag[cur>>1] != 0 && extra_check)
                {
                    (*extra_check) = 0;
                } 
                
                if(bub)
                {
                    tan_occ += get_btg_occ(bub, cur>>1);
                    if(tan_occ*side_rate >= beg_occ) return 0;
                    if(sink != (uint32_t)-1 && (tan_occ*side_rate >= sink_occ)) return 0;
                    if(tan_occ*total_rate >= (beg_occ + ((sink != (uint32_t)-1)?sink_occ : 0))) return 0;
                }
            }
        }
        

        if(is_vis[cur] == 0 && cur != (beg^1) && (sink == (uint32_t)-1 || cur != (sink^1)))
        {
            is_vis[cur] = 1;
            ncur = asg_arc_n(ug->g, cur);
            acur = asg_arc_a(ug->g, cur);
            for (i = 0; i < ncur; i++)
            {
                if(acur[i].del) continue;

                if(acur[i].v == beg) return 0;
                if(sink != (uint32_t)-1 && acur[i].v == sink) return 0;

                if(is_vis[acur[i].v]) continue;
                kv_push(uint32_t, stack->a, acur[i].v);
            }
        }
        is_vis[cur] = 1;
        

        cur^=1;
        if(is_vis[cur] == 0 && cur != (beg^1) && (sink == (uint32_t)-1 || cur != (sink^1)))
        {
            is_vis[cur] = 1;
            ncur = asg_arc_n(ug->g, cur);
            acur = asg_arc_a(ug->g, cur);
            for (i = 0; i < ncur; i++)
            {
                if(acur[i].del) continue;

                if(acur[i].v == beg) return 0;
                if(sink != (uint32_t)-1 && acur[i].v == sink) return 0;

                if(is_vis[acur[i].v]) continue;
                kv_push(uint32_t, stack->a, acur[i].v);
            }
        }
        is_vis[cur] = 1;
        no_first = 1;
    }

    if(bub)
    {
        if(tan_occ*side_rate >= beg_occ) return 0;
        if(sink != (uint32_t)-1 && (tan_occ*side_rate >= sink_occ)) return 0;
        if(tan_occ*total_rate >= (beg_occ + ((sink != (uint32_t)-1)?sink_occ : 0))) return 0;
    }
    
    return 1;
}


int find_bubble_tangle(bubble_type* bub, ma_ug_t* ug, uint8_t* is_vis, uint8_t* is_vis2, 
uint32_t v, double side_rate, double total_rate, kvec_t_u32_warp* stack, 
kvec_t_u32_warp* stack2, kvec_t_u32_warp* res_btg, kvec_t_u32_warp* res_utg, uint8_t* chain_flag, 
uint32_t* r_b_utg_beg, uint32_t* r_b_utg_sink, uint32_t* r_b_tg_beg, uint32_t* r_b_tg_sink,
uint32_t* r_utg_beg, uint32_t* r_utg_sink)
{
    (*r_b_utg_beg) = (*r_b_utg_sink) = (*r_b_tg_beg) = (*r_b_tg_sink) = (*r_utg_beg) = (*r_utg_sink) = (uint32_t)-1;
    ma_ug_t *bub_ug = bub->b_ug;
    ma_utg_t *u = NULL;
    uint32_t tan_occ = 0, cur, ncur, i, k, no_root = 0, v_occ, c_occ, utg_occ, w, btg_beg, btg_sink, utg_beg, utg_sink, is_t, extra_check;
    stack->a.n = 0;
    asg_arc_t *acur = NULL;

    memset(is_vis, 0, bub_ug->g->n_seq<<1);
    stack->a.n = 0;
    kv_push(uint32_t, stack->a, v);
    v_occ = get_btg_occ(bub, v>>1);
    while (stack->a.n > 0)
    {
        stack->a.n--;
        cur = stack->a.a[stack->a.n];
        if(is_vis[cur]) continue;
    
        c_occ = 0;
        if(no_root && cur == v) return 0;
        if(no_root && is_vis[cur] == 0 && is_vis[cur^1] == 0)
        {
            c_occ = get_btg_occ(bub, cur>>1);
            ///assume v_occ is the beg node, c_occ is the end node, which means tan_occ cannot be too large
            if((tan_occ*side_rate) < c_occ && (tan_occ*total_rate) < (c_occ + v_occ))
            {
                if(check_bubble_tangle(bub, bub->b_ug, v, cur^1, side_rate, total_rate, v_occ, c_occ, 
                is_vis2, stack2, NULL, NULL, NULL))
                {
                    check_bubble_tangle(bub, bub->b_ug, v, cur^1, side_rate, total_rate, v_occ, c_occ, 
                    is_vis2, stack2, res_btg, NULL, NULL);

                    for (k = 0; k < res_btg->a.n; k++) 
                    {
                        set_b_utg_weight_flag(bub, NULL, res_btg->a.a[k], chain_flag, 0, NULL);
                    }
                    btg_beg = v; btg_sink = cur^1;///b_utg id

                    u = &(bub_ug->u.a[btg_beg>>1]);
                    if((btg_beg&1)==1)
                    {
                        get_bubbles(bub, (u->a[0]>>32)>>1, (((u->a[0]>>32)&1)^1)==1?&w:NULL, 
                                                    (((u->a[0]>>32)&1)^1) == 0?&w:NULL, NULL, NULL, NULL);
                        (*r_b_tg_beg) = u->a[0]>>32;
                    } 
                    else
                    {
                        get_bubbles(bub, (u->a[u->n-1]>>32)>>1, ((u->a[u->n-1]>>32)&1)==1?&w:NULL, 
                                                    ((u->a[u->n-1]>>32)&1) == 0?&w:NULL, NULL, NULL, NULL);
                        (*r_b_tg_beg) = u->a[u->n-1]>>32;
                    }
                    utg_beg = w^1; ///ug id


                    u = &(bub_ug->u.a[btg_sink>>1]);
                    if((btg_sink&1)==1)
                    {
                        get_bubbles(bub, (u->a[0]>>32)>>1, (((u->a[0]>>32)&1)^1)==1?&w:NULL, 
                                                    (((u->a[0]>>32)&1)^1) == 0?&w:NULL, NULL, NULL, NULL);
                        (*r_b_tg_sink) = u->a[0]>>32;
                    } 
                    else
                    {
                        get_bubbles(bub, (u->a[u->n-1]>>32)>>1, ((u->a[u->n-1]>>32)&1)==1?&w:NULL, 
                                                    ((u->a[u->n-1]>>32)&1) == 0?&w:NULL, NULL, NULL, NULL);
                        (*r_b_tg_sink) = u->a[u->n-1]>>32;
                    }
                    utg_sink = w^1; ///ug id

                    is_t = check_bubble_tangle(NULL, ug, utg_beg, utg_sink, side_rate, total_rate, 
                    (uint32_t)-1, (uint32_t)-1, is_vis2, stack2, res_utg, chain_flag, &extra_check);

                    if(is_t == 1 && extra_check == 0)
                    {
                        for (k = utg_occ = 0; k < res_utg->a.n; k++) 
                        {
                            if(IF_HOM((res_utg->a.a[k]>>1), *bub)) continue;
                            utg_occ += ug->u.a[res_utg->a.a[k]>>1].n;
                        }

                        if(utg_occ*total_rate >= (v_occ+c_occ)) is_t = 0;
                    }

                    for (k = 0; k < res_btg->a.n; k++) 
                    {
                        set_b_utg_weight_flag(bub, NULL, res_btg->a.a[k], chain_flag, 1, NULL);
                    }
                    
                    if(is_t)
                    {
                        (*r_b_utg_beg) = btg_beg;
                        (*r_b_utg_sink) = btg_sink;
                        (*r_utg_beg) = utg_beg;
                        (*r_utg_sink) = utg_sink;
                        return is_t;
                    } 
                }
            }
        }
        
        is_vis[cur] = 1; 
        if(cur != (v^1))
        {
            ncur = asg_arc_n(bub_ug->g, cur);
            acur = asg_arc_a(bub_ug->g, cur);
            for (i = 0; i < ncur; i++)
            {
                if(acur[i].del) continue;
                if(acur[i].v == v) return 0;
                if(is_vis[acur[i].v]) continue;
                kv_push(uint32_t, stack->a, acur[i].v);
            }
        }
        

        if(no_root) tan_occ += c_occ;
        if((tan_occ*side_rate) >= v_occ) return 0;
        no_root = 1;
    }

    if(tan_occ*side_rate >= v_occ) return 0;
    if(tan_occ*total_rate >= v_occ) return 0;
    //let one end as a tangle
    if(check_bubble_tangle(bub, bub->b_ug, v, (uint32_t)-1, side_rate, total_rate, v_occ, (uint32_t)-1, 
                    is_vis2, stack2, res_btg, NULL, NULL) == 0)
    {
        return 0;
    }

    for (k = 0; k < res_btg->a.n; k++) 
    {
        set_b_utg_weight_flag(bub, NULL, res_btg->a.a[k], chain_flag, 0, NULL);
    }

    btg_beg = v;
    u = &(bub_ug->u.a[btg_beg>>1]);
    if((btg_beg&1)==1)
    {
        get_bubbles(bub, (u->a[0]>>32)>>1, (((u->a[0]>>32)&1)^1)==1?&w:NULL, 
                                    (((u->a[0]>>32)&1)^1) == 0?&w:NULL, NULL, NULL, NULL);
        (*r_b_tg_beg) = u->a[0]>>32;
    } 
    else
    {
        get_bubbles(bub, (u->a[u->n-1]>>32)>>1, ((u->a[u->n-1]>>32)&1)==1?&w:NULL, 
                                    ((u->a[u->n-1]>>32)&1) == 0?&w:NULL, NULL, NULL, NULL);
        (*r_b_tg_beg) = u->a[u->n-1]>>32;
    }
    utg_beg = w^1;

    is_t = check_bubble_tangle(NULL, ug, utg_beg, (uint32_t)-1, side_rate, total_rate, 
                    (uint32_t)-1, (uint32_t)-1, is_vis2, stack2, res_utg, chain_flag, &extra_check);
    
    if(is_t == 1 && extra_check == 0)
    {
        for (k = utg_occ = 0; k < res_utg->a.n; k++) 
        {
            if(IF_HOM((res_utg->a.a[k]>>1), *bub)) continue;
            utg_occ += ug->u.a[res_utg->a.a[k]>>1].n;
        }

        if(utg_occ*total_rate >= v_occ) is_t = 0;
    }
    
    for (k = 0; k < res_btg->a.n; k++) 
    {
        set_b_utg_weight_flag(bub, NULL, res_btg->a.a[k], chain_flag, 1, NULL);
    }

    if(is_t)
    {
        (*r_b_utg_beg) = btg_beg;
        (*r_utg_beg) = utg_beg;
    }
    return is_t;
}

uint32_t get_utg_end_from_btg(bubble_type* bub, ma_ug_t *bub_ug, uint32_t v)
{
    ma_utg_t *u = &(bub_ug->u.a[v>>1]);
    
    if((v&1)==1)
    {
        return (u->a[0]>>32)^1;
    }
    else
    {
        return u->a[u->n-1]>>32;
    } 
}
 
void drop_g_edges_by_utg(bubble_type* bub, asg_t *bsg, ma_ug_t *bub_ug, kvec_t_u32_warp* res_btg, 
uint32_t b_utg_beg, uint32_t b_utg_sink)
{
    uint32_t i, k, v, root, nv;
    asg_arc_t *av = NULL;
    if(b_utg_beg != (uint32_t)-1)
    {
        root = b_utg_beg;
        v = get_utg_end_from_btg(bub, bub_ug, root);
        nv = asg_arc_n(bsg, v);
        av = asg_arc_a(bsg, v);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            av[i].del = 1;
            asg_arc_del(bsg, (av[i].v)^1, (av[i].ul>>32)^1, 1);
        }
    }
    
    if(b_utg_sink != (uint32_t)-1)
    {
        root = b_utg_sink;
        v = get_utg_end_from_btg(bub, bub_ug, root);
        nv = asg_arc_n(bsg, v);
        av = asg_arc_a(bsg, v);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            av[i].del = 1;
            asg_arc_del(bsg, (av[i].v)^1, (av[i].ul>>32)^1, 1);
        }
    }
    

    if(res_btg == NULL) return;

    for (k = 0; k < res_btg->a.n; k++)
    {
        root = res_btg->a.a[k];
        v = get_utg_end_from_btg(bub, bub_ug, root);
        nv = asg_arc_n(bsg, v);
        av = asg_arc_a(bsg, v);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            av[i].del = 1;
            asg_arc_del(bsg, (av[i].v)^1, (av[i].ul>>32)^1, 1);
        }


        root = res_btg->a.a[k]^1;
        v = get_utg_end_from_btg(bub, bub_ug, root);
        nv = asg_arc_n(bsg, v);
        av = asg_arc_a(bsg, v);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            av[i].del = 1;
            asg_arc_del(bsg, (av[i].v)^1, (av[i].ul>>32)^1, 1);
        }
    }
}

void debug_tangle_bubble(bubble_type* bub, long long beg_idx, long long end_idx, const char* command)
{
    // long long beg_idx = (long long)bub->b_g->n_seq - bub->tangle_bub;
    // long long end_idx = (long long)bub->b_g->n_seq - 1;
    long long i, j, k;
    ma_utg_t *u = NULL;
    uint32_t beg_utg, sink_utg, *a = NULL, n, btg_left, ori_left, btg_right, ori_right, root_0, root_1;
    for (i = beg_idx; i <= end_idx; i++)
    {
        get_bubbles(bub, i, &beg_utg, &sink_utg, &a, &n, NULL);
        fprintf(stderr, "\n(%lld) %s: beg=utg%.6ul, sink=utg%.6ul, n: %u\n", i, command, (beg_utg>>1)+1, (sink_utg>>1)+1, n);
        for (k = 0; k < n; k++)
        {
            fprintf(stderr, "mid=utg%.6ul\n", (a[k]>>1)+1);
        }
        
        for (j = 0; j < bub->b_ug->g->n_seq; j++)
        {
            u = &(bub->b_ug->u.a[j]);
            if(u->n) continue;
            for (k = 0; k < u->n; k++)
            {
                if((long long)(u->a[k]>>33) != i) continue;
                fprintf(stderr, "is the %lld-th bubble at btg%.6lldl\n", k, j+1);
                if(k > 0)
                {
                    btg_left = u->a[k-1]>>33;
                    ori_left = u->a[k-1]>>32&1;
                    get_bubbles(bub, btg_left, ori_left == 1?&root_0:NULL, ori_left == 0?&root_0:NULL, NULL, NULL, NULL);
                    fprintf(stderr, "left-utg%.6ul\n", (ori_left>>1)+1);
                }

                if(k + 1 < u->n)
                {
                    btg_right = u->a[k+1]>>33;
                    ori_right = (u->a[k+1]>>32&1)^1;
                    get_bubbles(bub, btg_right, ori_right == 1?&root_1:NULL, ori_right == 0?&root_1:NULL, NULL, NULL, NULL);
                    fprintf(stderr, "right-utg%.6ul\n", (ori_right>>1)+1);
                }
            }
        }
        
    }
}

uint32_t print_b_utg_occ(bubble_type* bub, uint32_t v)
{
    ma_ug_t *bub_ug = bub->b_ug;
    ma_utg_t *u = NULL;
    uint32_t k_i, k_j, *a = NULL, n, tan_occ = 0, beg, sink;

    
    u = &(bub_ug->u.a[v]);
    fprintf(stderr, "\nstart: %u-th bubble-utg-start (# bubbles: %u)\n", v, (uint32_t)u->n);
    for (k_i = 0; k_i < u->n; k_i++)
    {
        get_bubbles(bub, u->a[k_i]>>33, &beg, &sink, &a, &n, NULL);
        for (k_j = 0; k_j < n; k_j++)
        {
            tan_occ += bub->ug->u.a[a[k_j]>>1].n;
        }

        fprintf(stderr, "bid: %lu, n: %u, beg-utg%.6dl(%u), sink-utg%.6dl(%u)\n", u->a[k_i]>>33, n, (beg>>1)+1, beg&1, (sink>>1)+1, sink&1);
    }

    fprintf(stderr, "end: %u-th bubble-utg-end\n\n", v);
    return tan_occ;
}

void update_bsg(asg_t *bsg, kvec_asg_arc_t_warp* edges)
{
    asg_arc_t *t = NULL;
    uint32_t k, l, i, convex, max_i;
    long long max, nodeLen, baseLen, max_stop_nodeLen, max_stop_baseLen;

    for (k = 0; k < edges->a.n; k++)
    {   
        t = asg_arc_pushp(bsg);
        *t = edges->a.a[k];
    }
    bsg->is_srt = 0; free(bsg->idx); bsg->idx = 0;
    asg_cleanup(bsg);



    radix_sort_asg_e(edges->a.a, edges->a.a + edges->a.n);
    for (k = 1, l = 0; k <= edges->a.n; ++k)
    {
        if (k == edges->a.n || (edges->a.a[k].ul>>32) != (edges->a.a[l].ul>>32))
        {
            if(k - l > 1)
            {
                for (i = l, max = -1, max_i = (uint32_t)-1; i < k; i++)
                {
                    get_unitig(bsg, NULL, edges->a.a[i].v, &convex, &nodeLen, &baseLen, &max_stop_nodeLen, 
                                                                    &max_stop_baseLen, 1, NULL);
                    if(max < nodeLen) max = nodeLen, max_i = i;
                }

                ///fprintf(stderr, "k - l: %u, max_i: %u\n", k - l, max_i);
                for (i = l; i < k; i++)
                {
                    // fprintf(stderr, "i: %u, +t->ul>>32: %lu, t->v: %u\n", 
                    //                     i, edges->a.a[i].ul>>32, edges->a.a[i].v);
                    if(i == max_i) continue;
                    asg_arc_del(bsg, (edges->a.a[i].ul>>32), (edges->a.a[i].v), 1);
                    asg_arc_del(bsg, (edges->a.a[i].v)^1, (edges->a.a[i].ul>>32)^1, 1);
                    ///edges->a.a[i].del = 1;
                }
            }
            l = k;
        }
    }
    
    asg_cleanup(bsg);
    
}

void resolve_bubble_chain_tangle(ma_ug_t* ug, bubble_type* bub)
{
    // double index_time = yak_realtime();
    ma_ug_t *bub_ug = bub->b_ug;
    asg_t *bsg = bub->b_g;
    uint32_t k, i, v, n_vx, new_bub;
    n_vx = MAX((MAX(ug->g->n_seq<<1, bub->b_ug->g->n_seq<<1)), bub->b_g->n_seq<<1);
    buf_t b; memset(&b, 0, sizeof(buf_t));
    kvec_t_u32_warp stack, stack2, res_btg, res_utg;
    kv_init(stack.a); kv_init(stack2.a); kv_init(res_btg.a); kv_init(res_utg.a);
    kvec_asg_arc_t_warp edges; kv_init(edges.a);
    uint8_t *is_vis = NULL; CALLOC(is_vis, n_vx);
    uint8_t *is_vis2 = NULL; CALLOC(is_vis2, n_vx);
    uint8_t *is_used = NULL; CALLOC(is_used, n_vx);
    uint8_t *chain_flag = NULL; CALLOC(chain_flag, n_vx);
    kvec_t(uint64_t) occ_idx; kv_init(occ_idx); uint64_t tmp, *p = NULL;
    double side_rate = 2.5, total_rate = 8;
    uint32_t b_utg_beg, b_utg_sink, b_tg_beg, b_tg_sink, utg_beg, utg_sink;
    

    while(1)
    {
        occ_idx.n = 0; edges.a.n = 0;
        if(n_vx < (uint32_t)(MAX((MAX(ug->g->n_seq<<1, bub->b_ug->g->n_seq<<1)), bub->b_g->n_seq<<1)))
        {
            n_vx = MAX((MAX(ug->g->n_seq<<1, bub->b_ug->g->n_seq<<1)), bub->b_g->n_seq<<1);
            is_vis = (uint8_t*)realloc(is_vis, n_vx);
            is_vis2 = (uint8_t*)realloc(is_vis2, n_vx);
            is_used = (uint8_t*)realloc(is_used, n_vx);
            chain_flag = (uint8_t*)realloc(chain_flag, n_vx);
        }
        memset(is_vis, 0, n_vx);
        memset(is_vis2, 0, n_vx);
        memset(is_used, 0, n_vx);
        memset(chain_flag, 0, n_vx);
        if(bub->num.n > 0) bub->num.n--;
        new_bub = bub->b_g->n_seq;
        //label all unitigs in bubble chain
        for (k = 0; k < bub_ug->g->n_seq; k++)
        {
            kv_pushp(uint64_t, occ_idx, &p);
            (*p) = get_btg_occ(bub, k);
            (*p) <<= 32; (*p) += k;
            set_b_utg_weight_flag(bub, NULL, k<<1, chain_flag, 1, NULL);
        }
        radix_sort_hc64(occ_idx.a, occ_idx.a + occ_idx.n);
        for (k = 0; k < occ_idx.n>>1; ++k) 
        {
            tmp = occ_idx.a[k];
            occ_idx.a[k] = occ_idx.a[occ_idx.n - k - 1];
            occ_idx.a[occ_idx.n - k - 1] = tmp;
        }

        for (k = 0; k < bub_ug->g->n_seq; k++)
        {
            v = ((uint32_t)(occ_idx.a[k]))<<1;
            if(is_used[v] == 0 && asg_arc_n(bub_ug->g, v) > 0)
            {
                if(find_bubble_tangle(bub, ug, is_vis, is_vis2, v, side_rate, total_rate, &stack, &stack2, 
                        &res_btg, &res_utg, chain_flag, &b_utg_beg, &b_utg_sink, &b_tg_beg, &b_tg_sink, 
                        &utg_beg, &utg_sink))
                {
                    if(utg_beg != (uint32_t)-1 && (!IF_HOM(utg_beg>>1, *bub)))
                    {
                        kv_push(uint32_t, res_utg.a, utg_beg);
                    } 
                    if(utg_sink != (uint32_t)-1 && (!IF_HOM(utg_sink>>1, *bub)))
                    {
                        kv_push(uint32_t, res_utg.a, utg_sink);
                    }

                    for (i = 0; i < res_btg.a.n; i++)
                    {
                        is_used[res_btg.a.a[i]] = 1; 
                        is_used[res_btg.a.a[i]^1] = 1; 
                    }
                    if(b_utg_beg != (uint32_t)-1) is_used[b_utg_beg] = 1;
                    if(b_utg_sink != (uint32_t)-1) is_used[b_utg_sink] = 1;
                    if(b_tg_beg != (uint32_t)-1) b_tg_beg>>=1;
                    if(b_tg_sink != (uint32_t)-1) b_tg_sink>>=1;    
                    update_bubble_graph(&res_utg, utg_beg, b_tg_beg, utg_sink, b_tg_sink, 
                    bub, &edges, bsg, NULL, NULL, ug, NULL, 0);
                    drop_g_edges_by_utg(bub, bsg, bub_ug, &res_btg, b_utg_beg, b_utg_sink);
                }
            }
            

            v ^= 1;
            if(is_used[v] == 0 && asg_arc_n(bub_ug->g, v) > 0)
            {

                if(find_bubble_tangle(bub, ug, is_vis, is_vis2, v, side_rate, total_rate, &stack, &stack2, 
                        &res_btg, &res_utg, chain_flag, &b_utg_beg, &b_utg_sink, &b_tg_beg, &b_tg_sink, 
                        &utg_beg, &utg_sink))
                {
                    if(utg_beg != (uint32_t)-1 && (!IF_HOM(utg_beg>>1, *bub)))
                    {
                        kv_push(uint32_t, res_utg.a, utg_beg);
                    } 
                    if(utg_sink != (uint32_t)-1 && (!IF_HOM(utg_sink>>1, *bub)))
                    {
                        kv_push(uint32_t, res_utg.a, utg_sink);
                    }

                    for (i = 0; i < res_btg.a.n; i++)
                    {
                        is_used[res_btg.a.a[i]] = 1; 
                        is_used[res_btg.a.a[i]^1] = 1; 
                    }

                    if(b_utg_beg != (uint32_t)-1) is_used[b_utg_beg] = 1;
                    if(b_utg_sink != (uint32_t)-1) is_used[b_utg_sink] = 1;
                    if(b_tg_beg != (uint32_t)-1) b_tg_beg>>=1;
                    if(b_tg_sink != (uint32_t)-1) b_tg_sink>>=1;
                    /*******************************for debug************************************/
                    // if(utg_beg == (utg_sink^1))
                    // {
                    //     print_b_utg_occ(bub, b_utg_beg>>1);
                    //     print_b_utg_occ(bub, b_utg_sink>>1);
                    //     print_b_utg_occ(bub, 42);
                    //     ///print_debug_bubble_graph(bub, ug, asm_opt.output_file_name);
                    // }
                    /*******************************for debug************************************/
                    update_bubble_graph(&res_utg, utg_beg, b_tg_beg, utg_sink, b_tg_sink, bub, &edges, bsg, NULL, NULL, ug, NULL, 0);
                    drop_g_edges_by_utg(bub, bsg, bub_ug, &res_btg, b_utg_beg, b_utg_sink);
                    ///fprintf(stderr, "->>>>>>beg=btg%.6ul, sink=btg%.6ul\n", (b_utg_beg>>1)+1, (b_utg_sink>>1)+1);
                }
            }
        }
        kv_push(uint32_t, bub->num, bub->list.n);
        new_bub = bub->b_g->n_seq - new_bub;
        bub->tangle_bub += new_bub;
        ///actually not useful, and may have bug when one bubble at multipe chains
        if(new_bub) update_bub_b_s_idx(bub);
        update_bsg(bsg, &edges);
        
        ma_ug_destroy(bub_ug);
        bub_ug = ma_ug_gen(bub->b_g);
        bub->b_ug = bub_ug;
        ///fprintf(stderr, "new_bub: %u, bub->tangle_bub: %lu\n", new_bub, bub->tangle_bub);
        if(new_bub == 0) break;
    }

    kv_destroy(bub->chain_weight);
    ma_utg_t *u = NULL;
    bub_ug = bub->b_ug;
    kv_malloc(bub->chain_weight, bub_ug->u.n); bub->chain_weight.n = bub_ug->u.n;
    for (i = 0; i < bub_ug->u.n; i++)
    {
        u = &(bub_ug->u.a[i]);
        bub->chain_weight.a[i].id = i;
        // if(u->n <= 1) ///not a chain
        // {
        //     bub->chain_weight.a[i].b_occ = bub->chain_weight.a[i].g_occ = 0;
        //     bub->chain_weight.a[i].del = 1;
        // } 
        // else
        {
            bub->chain_weight.a[i].del = 0;
            calculate_chain_weight(u, bub, ug, &(bub->chain_weight.a[i]));
        }
    }
    qsort(bub->chain_weight.a, bub->chain_weight.n, sizeof(chain_w_type), cmp_chain_weight);

    ///debug_tangle_bubble(bub);

    free(is_vis); free(is_vis2); free(is_used); free(chain_flag); free(b.b.a); 
    kv_destroy(occ_idx); kv_destroy(stack.a); kv_destroy(stack2.a); 
    kv_destroy(res_btg.a); kv_destroy(res_utg.a); kv_destroy(edges.a);

    ///print_debug_bubble_graph(bub, ug, asm_opt.output_file_name);
    // fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
}

void update_bubble_chain(ma_ug_t* ug, bubble_type* bub, uint32_t is_middle, uint32_t is_end)
{
    // double index_time = yak_realtime();
    if(bub->b_ug) ma_ug_destroy(bub->b_ug);
    if(bub->chain_weight.a) kv_destroy(bub->chain_weight);
    kvec_t_u32_warp broken;
    kv_init(broken.a);
    kvec_asg_arc_t_warp edges;
    kv_init(edges.a);
    ma_utg_t *u = NULL;
    asg_arc_t *t = NULL;
    asg_t *sg = ug->g;
    pdq pq;
    init_pdq(&pq, sg->n_seq<<1);
    asg_t *bsg = bub->b_g;
    ma_ug_t *bub_ug = NULL;
    bub_ug = ma_ug_gen(bub->b_g);
    uint32_t i, j, k_i, rId_0, ori_0, root_0, rId_1, ori_1, root_1, n_vtx = sg->n_seq<<1, new_bub;
    uint32_t *pre = NULL; MALLOC(pre, n_vtx); 
    uint8_t* vis_flag = NULL; MALLOC(vis_flag, ug->g->n_seq*2);
    kvec_t_u32_warp stack; kv_init(stack.a);
    ///chain_w_type x;
    ///uint64_t end_thres;

    uint8_t *bsg_idx = NULL; CALLOC(bsg_idx, n_vtx>>1); 
    for (i = 0; i < bub_ug->u.n; i++)
    {
        u = &(bub_ug->u.a[i]);
        if(u->n == 0) continue;
        for (k_i = 0; k_i < u->n; k_i++)
        {
            uint32_t *a, n;
            get_bubbles(bub, u->a[k_i]>>33, &root_0, &root_1, &a, &n, NULL);
            for (j = 0; j < n; j++)
            {
                bsg_idx[a[j]>>1] = 1;
            }
            bsg_idx[root_0>>1] = 1;
            bsg_idx[root_1>>1] = 1;
        }
    }

    if(bub->num.n > 0) bub->num.n--;
    new_bub = bub->b_g->n_seq;
    for (i = 0; i < bub_ug->u.n; i++)
    {
        u = &(bub_ug->u.a[i]);
        if(u->n == 0) continue;
        ///end_thres = calculate_chain_weight(u, bub, ug, &x);
        if(is_middle)
        {
            for (k_i = 0; k_i < u->n; k_i++)
            {
                if(k_i+1 >= u->n) continue;
                ///note: must igore .del here, since bsg might be changed
                t = &(arc_first(bsg, u->a[k_i]>>32));
                if(t->el == 1) continue;

                rId_0 = u->a[k_i]>>33;
                ori_0 = u->a[k_i]>>32&1;
                get_bubbles(bub, rId_0, ori_0 == 1?&root_0:NULL, ori_0 == 0?&root_0:NULL, NULL, NULL, NULL);
                
                rId_1 = u->a[k_i+1]>>33;
                ori_1 = (u->a[k_i+1]>>32&1)^1;
                get_bubbles(bub, rId_1, ori_1 == 1?&root_1:NULL, ori_1 == 0?&root_1:NULL, NULL, NULL, NULL);

                broken.a.n = 0;
                get_related_bub_nodes(&broken, bub, &pq, sg, pre, root_0, root_1, NULL);
                get_related_bub_nodes(&broken, bub, &pq, sg, pre, root_1, root_0, NULL);
                ///no need to cut the edge, we still have chance to flip by chain
                if(double_check_broken_bubble(ug->g, &broken, root_0^1, root_1^1, vis_flag, 
                                                                    ug->g->n_seq*2, &stack, NULL, NULL/**bsg, t**/) == 0)
                {
                    continue;
                }
                if(!IF_HOM(root_0>>1, *bub)) kv_push(uint32_t, broken.a, root_0);
                if(!IF_HOM(root_1>>1, *bub)) kv_push(uint32_t, broken.a, root_1);
                if(broken.a.n > 0)
                {
                    update_bubble_graph(&broken, root_0^1, rId_0, root_1^1, rId_1, bub, &edges, bsg, t, bsg_idx, ug, NULL, 1);
                }
            }
        }


        if(is_end)
        {
            if(u->n >0 && arc_cnt(bub_ug->g, (i<<1)+1) == 0)
            {
                rId_0 = u->a[0]>>33;
                ori_0 = (u->a[0]>>32&1)^1;
                get_bubbles(bub, rId_0, ori_0 == 1?&root_0:NULL, ori_0 == 0?&root_0:NULL, NULL, NULL, NULL);
                broken.a.n = 0;
                get_related_bub_nodes(&broken, bub, &pq, sg, pre, root_0, (uint32_t)-1, NULL);

                if(double_check_broken_bubble(ug->g, &broken, root_0^1, (uint32_t)-1, vis_flag, 
                                                                    ug->g->n_seq*2, &stack, NULL, NULL))
                {
                    if(!IF_HOM(root_0>>1, *bub)) kv_push(uint32_t, broken.a, root_0);
                    if(broken.a.n > 0)
                    {
                        ///fprintf(stderr, "root_0: utg%.6ul, broken.a.n: %u\n", (root_0>>1)+1, (uint32_t)broken.a.n);
                        update_bubble_graph(&broken, root_0^1, rId_0, (uint32_t)-1, (uint32_t)-1, bub, &edges, bsg, NULL, bsg_idx, ug, NULL, 0);
                    }
                }   
            }

        
            if(u->n >0 && arc_cnt(bub_ug->g, i<<1) == 0)
            {
                rId_1 = u->a[u->n-1]>>33;
                ori_1 = u->a[u->n-1]>>32&1;
                get_bubbles(bub, rId_1, ori_1 == 1?&root_1:NULL, ori_1 == 0?&root_1:NULL, NULL, NULL, NULL);
                broken.a.n = 0;
                get_related_bub_nodes(&broken, bub, &pq, sg, pre, root_1, (uint32_t)-1, bsg_idx);

                if(double_check_broken_bubble(ug->g, &broken, root_1^1, (uint32_t)-1, vis_flag, 
                                                                    ug->g->n_seq*2, &stack, NULL, NULL))
                {
                    if(!IF_HOM(root_1>>1, *bub)) kv_push(uint32_t, broken.a, root_1);
                    if(broken.a.n > 0)
                    {
                        ///fprintf(stderr, "root_1: utg%.6ul, broken.a.n: %u\n", (root_1>>1)+1, (uint32_t)broken.a.n);
                        update_bubble_graph(&broken, (uint32_t)-1, (uint32_t)-1, root_1^1, rId_1, bub, &edges, bsg, NULL, bsg_idx, ug, NULL, 0);
                    }

                }
            }
        }
    }
    kv_push(uint32_t, bub->num, bub->list.n);
    new_bub = bub->b_g->n_seq - new_bub;
    if(is_end) bub->b_end_bub += new_bub;
    ///actually not useful, and may have bug when one bubble at multipe chains
    if(new_bub) update_bub_b_s_idx(bub);


    for (i = 0; i < edges.a.n; i++)
    {        
        t = asg_arc_pushp(bsg);
        *t = edges.a.a[i];
    }

    bsg->is_srt = 0; free(bsg->idx); bsg->idx = 0;
    asg_cleanup(bsg);
    ma_ug_destroy(bub_ug);
    destory_pdq(&pq);
    free(pre);
    kv_destroy(broken.a);
    kv_destroy(edges.a);
    free(bsg_idx);
    kv_destroy(stack.a); 
    free(vis_flag);

    bub->b_ug = ma_ug_gen(bub->b_g);
    bub_ug = bub->b_ug;
    kv_malloc(bub->chain_weight, bub_ug->u.n); bub->chain_weight.n = bub_ug->u.n;
    for (i = 0; i < bub_ug->u.n; i++)
    {
        u = &(bub_ug->u.a[i]);
        bub->chain_weight.a[i].id = i;
        // if(u->n <= 1) ///not a chain
        // {
        //     bub->chain_weight.a[i].b_occ = bub->chain_weight.a[i].g_occ = 0;
        //     bub->chain_weight.a[i].del = 1;
        // } 
        // else
        {
            bub->chain_weight.a[i].del = 0;
            calculate_chain_weight(u, bub, ug, &(bub->chain_weight.a[i]));
        }
    }

    qsort(bub->chain_weight.a, bub->chain_weight.n, sizeof(chain_w_type), cmp_chain_weight);



    /**
    uint32_t d_v, d_u, v;
    for (i = 0; i < bsg->n_arc; i++)
    {
        d_v = (uint32_t)(bsg->arc[i].ul>>32);
        d_u = bsg->arc[i].v;
        for (v = 0; v < bsg->n_arc; v++)
        {
            if(((bsg->arc[v].ul>>32) == (d_u^1)) && (bsg->arc[v].v == (d_v^1))) break;      
        }

        if(v == bsg->n_arc)
        {
            fprintf(stderr, "hahaha, el: %u, ul>>33: %lu, ul&1: %lu, v>>1: %u, v&1: %u\n", 
                    bsg->arc[i].el, bsg->arc[i].ul>>33, (bsg->arc[i].ul>>32)&1, bsg->arc[i].v>>1, bsg->arc[i].v&1);
        }
        // else
        // {
        //     fprintf(stderr, "hehehe, el: %u, ul>>33: %lu, ul&1: %lu, v>>1: %u, v&1: %u\n", 
        //             bsg->arc[i].el, bsg->arc[i].ul>>33, (bsg->arc[i].ul>>32)&1, bsg->arc[i].v>>1, bsg->arc[i].v&1);
        // }
    }

    asg_arc_t *av = NULL, *au = NULL;
    uint32_t nv, nu;
    for (v = 0; v < (uint32_t)(bsg->n_seq<<1); v++)
    {
        av = asg_arc_a(bsg, v);
        nv = asg_arc_n(bsg, v);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;

            au = asg_arc_a(bsg, av[i].v^1);
            nu = asg_arc_n(bsg, av[i].v^1);
            for (k_i = 0; k_i < nu; k_i++)
            {
                if(au[k_i].del) continue;
                if(au[k_i].v == (v^1)) break;
            }
            if(k_i == nu) fprintf(stderr, "hahaha: v: %u, u: %u\n", v, av[i].v);
        }        
    }
    **/
    // fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
}

void set_b_utg_weight_flag(bubble_type* bub, buf_t* b, uint32_t v, uint8_t* vis_flag, uint32_t flag, uint32_t* occ)
{
    ma_ug_t *bub_ug = bub->b_ug;
    long long nodeLen, baseLen, max_stop_nodeLen, max_stop_baseLen;
    ma_utg_t *u = NULL;
    uint32_t convex, k, k_i, k_j, *a, n, beg, sink;
    if(b)
    {
        b->b.n = 0;
        get_unitig(bub_ug->g, NULL, v, &convex, &nodeLen, &baseLen, &max_stop_nodeLen, 
                                                                    &max_stop_baseLen, 1, b);
    }
    
    if(occ) (*occ) = 0;
    for (k = 0; k < (b?b->b.n:1); k++)
    {
        u = &(bub_ug->u.a[b?(b->b.a[k]>>1):(v>>1)]);
        if(u->n == 0) continue;
        for (k_i = 0; k_i < u->n; k_i++)
        {
            get_bubbles(bub, u->a[k_i]>>33, &beg, &sink, &a, &n, NULL);

            for (k_j = 0; k_j < n; k_j++)
            {
                if(vis_flag) vis_flag[a[k_j]>>1] = flag;
                if(occ) (*occ) += bub->ug->u.a[a[k_j]>>1].n;
            }
            if(beg != (uint32_t)-1 && vis_flag) vis_flag[beg>>1] = flag;
            if(sink != (uint32_t)-1 && vis_flag) vis_flag[sink>>1] = flag;
        }
    }
}

void set_b_utg_weight_flag_xor(bubble_type* bub, ma_ug_t *bub_ug, buf_t* b, uint32_t v, uint8_t* vis_flag, uint32_t flag, uint32_t* occ)
{
    long long nodeLen, baseLen, max_stop_nodeLen, max_stop_baseLen;
    ma_utg_t *u = NULL;
    uint32_t convex, k, k_i, k_j, *a, n, beg, sink;
    if(b)
    {
        b->b.n = 0;
        get_unitig(bub_ug->g, NULL, v, &convex, &nodeLen, &baseLen, &max_stop_nodeLen, 
                                                                    &max_stop_baseLen, 1, b);
    }
    
    if(occ) (*occ) = 0;
    for (k = 0; k < (b?b->b.n:1); k++)
    {
        u = &(bub_ug->u.a[b?(b->b.a[k]>>1):(v>>1)]);
        if(u->n == 0) continue;
        for (k_i = 0; k_i < u->n; k_i++)
        {
            get_bubbles(bub, u->a[k_i]>>33, &beg, &sink, &a, &n, NULL);

            for (k_j = 0; k_j < n; k_j++)
            {
                if(vis_flag) vis_flag[a[k_j]>>1] ^= flag;
                if(occ) (*occ) += bub->ug->u.a[a[k_j]>>1].n;
            }
            if(beg != (uint32_t)-1 && vis_flag) vis_flag[beg>>1] ^= flag;
            if(sink != (uint32_t)-1 && vis_flag) vis_flag[sink>>1] ^= flag;
        }
    }
}

double dfs_weight(uint32_t v, uint8_t* vis_flag, uint8_t* is_vis, hc_links* link, 
kvec_t_u32_warp* stack, kvec_t_u32_warp* result, uint32_t e_flag, uint32_t ava_flag,
uint32_t* link_occ)
{
    uint32_t cur, i, next = (uint32_t)-1;
    stack->a.n = 0;
    kv_push(uint32_t, stack->a, v);
    double w = 0;
    if(link_occ) (*link_occ) = 0;
    while (stack->a.n > 0)
    {
        stack->a.n--;
        cur = stack->a.a[stack->a.n];
        if(is_vis[cur]) continue;
        is_vis[cur] = 1;
        if(cur!=v && vis_flag[cur] != ava_flag) continue;
        for (i = 0; i < link->a.a[cur].e.n; i++)
        {
            if(link->a.a[cur].e.a[i].del) continue;
            next = link->a.a[cur].e.a[i].uID;
            ///if(vis_flag[next] == e_flag)
            if(vis_flag[next]&e_flag)
            {
                w += link->a.a[cur].e.a[i].weight;
                if(link_occ) (*link_occ) += link->a.a[cur].e.a[i].occ;
                continue;
            }
            if(is_vis[next]) continue;
            if(vis_flag[next] != ava_flag) continue;
            kv_push(uint32_t, stack->a, next);
        }
    }
    return w;
}

void if_conflict_utg(uint32_t root, uint32_t* aim_0, uint32_t* aim_1, ma_ug_t* ug, uint8_t* vis_flag,
uint8_t* is_vis_2, uint32_t ava_flag, kvec_t_u32_warp* stack)
{
    uint32_t n_vx = ug->g->n_seq<<1, k, cur, ncur;
    asg_arc_t *acur = NULL;
    memset(is_vis_2, 0, n_vx);

    stack->a.n = 0;
    kv_push(uint32_t, stack->a, root);
    while (stack->a.n > 0)
    {
        stack->a.n--;
        cur = stack->a.a[stack->a.n];
        if(is_vis_2[cur]) continue;
        is_vis_2[cur] = 1;
       
        ncur = asg_arc_n(ug->g, cur);
        acur = asg_arc_a(ug->g, cur);
        for (k = 0; k < ncur; k++)
        {
            if(acur[k].del) continue;
            if(is_vis_2[acur[k].v]) continue;
            if(vis_flag[acur[k].v>>1] != 0 && vis_flag[acur[k].v>>1] != ava_flag)
            {
                if(aim_0 && (acur[k].v>>1) == (*aim_0)) continue;
                if(aim_1 && (acur[k].v>>1) == (*aim_1)) continue;
                break;
            }
            kv_push(uint32_t, stack->a, acur[k].v);
        }

        if(k < ncur) return;
    }

    for (k = 0; k < n_vx; k++)
    {
        if(is_vis_2[k] && vis_flag[k>>1] == 0)
        {
            ///fprintf(stderr, "******************k=utg%.6ul, vis_flag: %u\n", (k>>1)+1, vis_flag[k>>1]);
            vis_flag[k>>1] = ava_flag;
        } 
        
    }
}

double get_chain_weight(bubble_type* bub, ma_ug_t *bub_ug, buf_t* b, uint32_t v, uint32_t convex_source, hc_links* link, 
uint8_t* vis_flag, uint8_t* is_vis, ma_ug_t* ug, kvec_t_u32_warp* stack, kvec_t_u32_warp* result, 
uint32_t e_flag, uint32_t ava_flag, kvec_t_u32_warp* res_utg, uint32_t* link_occ)
{
    long long nodeLen, baseLen, max_stop_nodeLen, max_stop_baseLen;
    ma_utg_t *u = NULL;
    uint32_t convex, k, k_i, k_j, *a, n, beg, sink, uID, root, cur, ncur, n_vx = ug->g->n_seq<<1, occ;
    asg_arc_t *acur = NULL;
    double w = 0;
    b->b.n = 0;
    get_unitig(bub_ug->g, NULL, v, &convex, &nodeLen, &baseLen, &max_stop_nodeLen, 
                                                                    &max_stop_baseLen, 1, b);
    memset(is_vis, 0, n_vx);

    for (k = 0; k < b->b.n; k++)
    {
        u = &(bub_ug->u.a[b->b.a[k]>>1]);
        if(u->n == 0) continue;

        for (k_i = 0; k_i < u->n; k_i++)
        {
            get_bubbles(bub, u->a[k_i]>>33, &beg, &sink, &a, &n, NULL);
            for(k_j = 0; k_j < n; k_j++) is_vis[a[k_j]] = is_vis[a[k_j]^1] = 1;
            if(beg != (uint32_t)-1) is_vis[beg] = is_vis[beg^1] = 1; 
            if(sink != (uint32_t)-1) is_vis[sink] = is_vis[sink^1] = 1; 
        }
    }

    u = &(bub_ug->u.a[v>>1]);
    if((v&1)==0)
    {
        get_bubbles(bub, (u->a[0]>>32)>>1, (((u->a[0]>>32)&1)^1)==1?&root:NULL, 
                                    (((u->a[0]>>32)&1)^1) == 0?&root:NULL, NULL, NULL, NULL);
    } 
    else
    {
        get_bubbles(bub, (u->a[u->n-1]>>32)>>1, ((u->a[u->n-1]>>32)&1)==1?&root:NULL, 
                                    ((u->a[u->n-1]>>32)&1) == 0?&root:NULL, NULL, NULL, NULL);
    }
    
    root ^= 1;
    ///fprintf(stderr, "root=utg%.6dl\n", (root>>1)+1);
    is_vis[root] = 0; 
    stack->a.n = 0;
    kv_push(uint32_t, stack->a, root);
    while (stack->a.n > 0)
    {
        stack->a.n--;
        cur = stack->a.a[stack->a.n];
        if(is_vis[cur]) continue;
        is_vis[cur] = 1;
        if(vis_flag[cur>>1] == 0) vis_flag[cur>>1] = ava_flag;
       
        ncur = asg_arc_n(ug->g, cur);
        acur = asg_arc_a(ug->g, cur);
        for (k = 0; k < ncur; k++)
        {
            if(acur[k].del) continue;
            if(is_vis[acur[k].v]) continue;
            if(vis_flag[acur[k].v>>1] != 0 && vis_flag[acur[k].v>>1] != ava_flag) continue;
            kv_push(uint32_t, stack->a, acur[k].v);
        }
    }



    uint32_t aim_0, aim_1, root_source;
    aim_0 = root>>1;
    u = &(bub_ug->u.a[convex_source>>1]);
    if((convex_source&1)==1)
    {
        get_bubbles(bub, (u->a[0]>>32)>>1, (((u->a[0]>>32)&1)^1)==1?&root_source:NULL, 
                                    (((u->a[0]>>32)&1)^1) == 0?&root_source:NULL, NULL, NULL, NULL);
    } 
    else
    {
        get_bubbles(bub, (u->a[u->n-1]>>32)>>1, ((u->a[u->n-1]>>32)&1)==1?&root_source:NULL, 
                                    ((u->a[u->n-1]>>32)&1) == 0?&root_source:NULL, NULL, NULL, NULL);
    }
    root_source ^= 1;
    aim_1 = root_source>>1;

    ///fprintf(stderr, "aim_0=utg%.6ul, aim_1=utg%.6ul\n", aim_0+1, aim_1+1);

    cur = root_source;
    ncur = asg_arc_n(ug->g, cur);
    acur = asg_arc_a(ug->g, cur);
    for (k_i = 0; k_i < ncur; k_i++)
    {
        if(acur[k_i].del) continue;
        if(vis_flag[acur[k_i].v>>1] != 0) continue;
        if_conflict_utg(acur[k_i].v, &aim_0, &aim_1, ug, vis_flag, is_vis, ava_flag, stack);
    }


    for (k = 0; k < ug->g->n_seq; k++)
    {
        if(vis_flag[k] == ava_flag)
        {
            cur = k<<1;
            ncur = asg_arc_n(ug->g, cur);
            acur = asg_arc_a(ug->g, cur);
            for (k_i = 0; k_i < ncur; k_i++)
            {
                if(acur[k_i].del) continue;
                if(vis_flag[acur[k_i].v>>1] != 0) continue;
                if_conflict_utg(acur[k_i].v, &aim_0, &aim_1, ug, vis_flag, is_vis, ava_flag, stack);
            }



            cur = (k<<1)+1;
            ncur = asg_arc_n(ug->g, cur);
            acur = asg_arc_a(ug->g, cur);
            for (k_i = 0; k_i < ncur; k_i++)
            {
                if(acur[k_i].del) continue;
                if(vis_flag[acur[k_i].v>>1] != 0) continue;
                if_conflict_utg(acur[k_i].v, &aim_0, &aim_1, ug, vis_flag, is_vis, ava_flag, stack);
            }
        }
    }
    






    memset(is_vis, 0, n_vx);
    if(link_occ) (*link_occ) = 0;
    for (k = result->a.n = 0, w = 0; k < b->b.n; k++)
    {
        u = &(bub_ug->u.a[b->b.a[k]>>1]);
        if(u->n == 0) continue;
        for (k_i = 0; k_i < u->n; k_i++)
        {
            get_bubbles(bub, u->a[k_i]>>33, NULL, NULL, &a, &n, NULL);

            for (k_j = 0; k_j < n; k_j++)
            {
                uID = a[k_j]>>1;
                w += dfs_weight(uID, vis_flag, is_vis, link, stack, result, e_flag, ava_flag, &occ);
                if(link_occ) (*link_occ) += occ;
            }
        }
    }


    for (k = 0; k < ug->g->n_seq; k++)
    {
        if(vis_flag[k] == ava_flag)
        {
            vis_flag[k] = 0;
            if(res_utg && (!IF_HOM(k, *bub)))
            {
                kv_push(uint32_t, res_utg->a, k<<1);
            }
        }
    }
    

    return w;
}

int double_check_bub_branch(asg_arc_t *t, ma_ug_t *bs_ug, double *e_w, uint32_t *e_occ, double cutoff, uint32_t max_w_occ)
{
    uint32_t v = t->v^1, w = (t->ul>>32)^1, i, nv, rv, max_i, w_i, *a_occ = NULL;
    asg_arc_t *av = NULL;
    double *aw = NULL, max_w = cutoff - 1, w_w = 1;
    av = asg_arc_a(bs_ug->g, v);
    nv = asg_arc_n(bs_ug->g, v);
    aw = (&e_w[bs_ug->g->idx[v]>>32]);
    a_occ = (&e_occ[bs_ug->g->idx[v]>>32]);

    if(nv <= 1) return 1;

    for (i = rv = 0, max_i = w_i = (uint32_t)-1; i < nv; i++)
    {
        if(av[i].del) continue;
        rv++;
        if(av[i].v == w)
        {
            w_i = i;
            w_w = aw[i];
            continue;
        }
        if(max_i == (uint32_t)-1)
        {
            max_i = i;
            max_w = aw[i];
        }
        else if(max_w < aw[i])
        {
            max_i = i;
            max_w = aw[i];
        }
    }

    if(rv <= 1) return 1; ///must be here

    if(max_i == (uint32_t)-1 || w_i == (uint32_t)-1) return 0;
    ///if(max_w <= max_w_cutoff) return 0; //must be <=
    if(a_occ[max_i] <= max_w_occ) return 0; //must be <=
    
    if(w_w*cutoff < max_w) return 1;
    return 0;
}

void clean_bubble_chain_by_HiC(ma_ug_t* ug, hc_links* link, bubble_type* bub)
{   
    double index_time = yak_realtime();
    ma_ug_t *bs_ug = bub->b_ug;
    uint32_t v, u, i, m, max_i, nv, rv, n_vx, root, flag_pri = 1, flag_aux = 2, flag_ava = 4, occ;
    double w, cutoff = 2/**, max_w_cutoff = MAX(MIN(100*OFFSET_RATE_MIN_W, OFFSET_RATE_MAX_W/100), OFFSET_RATE_MIN_W)**/;
    uint32_t max_w_occ = 4;
    asg_arc_t *av = NULL;
    n_vx = bs_ug->g->n_seq << 1;
    uint8_t *vis = NULL; CALLOC(vis, ug->g->n_seq<<1);
    uint8_t *is_vis = NULL; CALLOC(is_vis, ug->g->n_seq<<1);
    uint8_t *is_used = NULL; CALLOC(is_used, n_vx);
    uint8_t *dedup = NULL; CALLOC(dedup, ug->g->n_seq<<1);
    buf_t b; memset(&b, 0, sizeof(buf_t));
    kvec_t_u32_warp stack, result, res_utg;
    kv_init(stack.a); kv_init(result.a); kv_init(res_utg.a);
    double *e_w = NULL; MALLOC(e_w, bs_ug->g->n_arc);
    uint32_t *e_occ = NULL, *a_occ = NULL; CALLOC(e_occ, bs_ug->g->n_arc);
    double *aw = NULL, max_w = 0;
    kvec_asg_arc_t_warp edges; kv_init(edges.a);
    ma_ug_t *back_bs_ug = copy_untig_graph(bs_ug);

    for (i = 0; i < bs_ug->g->n_arc; i++)///weight of bs_ug's edges
    {
        e_w[i] = -1;
    }

    for (i = 0; i < bs_ug->g->n_seq; i++)///init all chain with flag_aux 
    {
        set_b_utg_weight_flag(bub, &b, i<<1, vis, flag_aux, NULL);
    }


    for (v = 0; v < n_vx; v++)
    {
        av = asg_arc_a(bs_ug->g, v);
        nv = asg_arc_n(bs_ug->g, v);
        aw = (&e_w[bs_ug->g->idx[v]>>32]);
        a_occ = (&e_occ[bs_ug->g->idx[v]>>32]);
        if(nv <= 1 || get_real_length(bs_ug->g, v, NULL) <= 1) continue;
        set_b_utg_weight_flag_xor(bub, bs_ug, &b, v^1, vis, flag_pri, NULL);

        ///fprintf(stderr, "\n******pri>btg%.6dl\n", (v>>1)+1);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            //fprintf(stderr, "aux>btg%.6dl\n", (av[i].v>>1)+1);
            w = get_chain_weight(bub, bs_ug, &b, av[i].v, v, link, vis, is_vis, ug, &stack, &result, flag_pri, flag_ava, NULL, &occ);
            ///fprintf(stderr, "aux>btg%.6dl, w: %f\n", (av[i].v>>1)+1, w);
            aw[i] = w;
            a_occ[i] = occ;
        }

        set_b_utg_weight_flag_xor(bub, bs_ug, &b, v^1, vis, flag_pri, NULL);
    }


    for (v = 0; v < n_vx; v++)
    {
        av = asg_arc_a(bs_ug->g, v);
        nv = asg_arc_n(bs_ug->g, v);
        aw = (&e_w[bs_ug->g->idx[v]>>32]);
        a_occ = (&e_occ[bs_ug->g->idx[v]>>32]);
        if(nv <= 1 || get_real_length(bs_ug->g, v, NULL) <= 1) continue;

        for (i = rv = 0, max_i = (uint32_t)-1; i < nv; i++)
        {
            if(av[i].del) continue;
            if(max_i == (uint32_t)-1)
            {
                max_i = i;
                max_w = aw[i];
            }
            else if(max_w < aw[i])
            {
                max_i = i;
                max_w = aw[i];
            }
            rv++;
        }

        if(max_i == (uint32_t)-1) continue;
        ///if(max_w <= max_w_cutoff) continue; //must be <=
        if(a_occ[max_i] <= max_w_occ) continue; //must be <=
        if(rv < 2) continue;

        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            if(i == max_i) continue;
            ///if((av[i].v>>1) == (v>>1) && aw[i] <= max_w_cutoff) continue; ///might be not reasonable
            if((av[i].v>>1) == (v>>1) && a_occ[i] <= max_w_occ) continue; ///might be not reasonable
            if(aw[i]*cutoff < max_w && double_check_bub_branch(&av[i], bs_ug, e_w, e_occ, cutoff, max_w_occ))
            {
                av[i].del = 1; asg_arc_del(bs_ug->g, (av[i].v)^1, (av[i].ul>>32)^1, 1);
            }
        }
    }

    uint32_t rId_0, ori_0, rId_1, ori_1, root_0, root_1, new_bub;
    if(bub->num.n > 0) bub->num.n--;
    new_bub = bub->b_g->n_seq;

    for (v = 0; v < n_vx; v++)
    {
        av = asg_arc_a(bs_ug->g, v);
        nv = asg_arc_n(bs_ug->g, v);
        rv = get_real_length(bs_ug->g, v, NULL);
        if(nv == rv) continue;
        if(rv != 1 || nv <= 1) continue;
        get_real_length(bs_ug->g, v, &u);
        u ^= 1;
        if(get_real_length(bs_ug->g, u, NULL) != 1) continue;
        drop_g_edges_by_utg(bub, bub->b_g, bs_ug, NULL, v, u);
        if(is_used[v] || is_used[u]) continue;

        is_used[v] = is_used[u] = 1;
        root = get_utg_end_from_btg(bub, bs_ug, v);
        rId_0 = root>>1;
        ori_0 = root&1;
        get_bubbles(bub, rId_0, ori_0 == 1?&root_0:NULL, ori_0 == 0?&root_0:NULL, NULL, NULL, NULL);

        root = get_utg_end_from_btg(bub, bs_ug, u);
        rId_1 = root>>1;
        ori_1 = root&1;
        get_bubbles(bub, rId_1, ori_1 == 1?&root_1:NULL, ori_1 == 0?&root_1:NULL, NULL, NULL, NULL);

        res_utg.a.n = 0;
        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, v^1, vis, flag_pri, NULL);
        get_chain_weight(bub, back_bs_ug, &b, u^1, v, link, vis, is_vis, ug, &stack, &result, flag_pri, flag_ava, &res_utg, NULL);
        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, v^1, vis, flag_pri, NULL);
        for (i = 0; i < res_utg.a.n; i++) dedup[res_utg.a.a[i]>>1] |= 1;

        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, u^1, vis, flag_pri, NULL);
        get_chain_weight(bub, back_bs_ug, &b, v^1, u, link, vis, is_vis, ug, &stack, &result, flag_pri, flag_ava, &res_utg, NULL);
        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, u^1, vis, flag_pri, NULL);
        for (; i < res_utg.a.n; i++) dedup[res_utg.a.a[i]>>1] |= 2;

        
        for (i = m = 0; i < res_utg.a.n; i++)
        {
            if(dedup[res_utg.a.a[i]>>1] == 3)
            {
                res_utg.a.a[m] = res_utg.a.a[i];
                m++;
            }
            dedup[res_utg.a.a[i]>>1] = 0;
        }
        res_utg.a.n = m;

        // fprintf(stderr, "res_utg.a.n: %u, m: %u, beg-utg%.6ul, sink-utg%.6ul\n", 
        //             res_utg.a.n, m, (root_0>>1)+1, (root_1>>1)+1);
        
        if(!IF_HOM(root_0>>1, *bub)) kv_push(uint32_t, res_utg.a, root_0);
        if(!IF_HOM(root_1>>1, *bub)) kv_push(uint32_t, res_utg.a, root_1);

        update_bubble_graph(&res_utg, root_0^1, rId_0, root_1^1, rId_1, bub, &edges, bub->b_g, NULL, NULL, ug, NULL, 0);

        ///fprintf(stderr, "\n******src-btg%.6ul------>dest-btg%.6ul\n", (v>>1)+1, (u>>1)+1);
    }

    kv_push(uint32_t, bub->num, bub->list.n);
    new_bub = bub->b_g->n_seq - new_bub;
    bub->cross_bub += new_bub;
    ///actually not useful, and may have bug when one bubble at multipe chains
    if(new_bub) update_bub_b_s_idx(bub);
    
    ///debug_tangle_bubble(bub, bub->b_g->n_seq - bub->cross_bub, bub->b_g->n_seq - 1, "Cross-tangle");
    update_bsg(bub->b_g, &edges);

    ma_ug_destroy(bs_ug);
    bs_ug = ma_ug_gen(bub->b_g);
    bub->b_ug = bs_ug;
    kv_destroy(bub->chain_weight);
    ma_utg_t *u_x = NULL;
    bs_ug = bub->b_ug;
    kv_malloc(bub->chain_weight, bs_ug->u.n); bub->chain_weight.n = bs_ug->u.n;
    for (i = 0; i < bs_ug->u.n; i++)
    {
        u_x = &(bs_ug->u.a[i]);
        bub->chain_weight.a[i].id = i;
        // if(u->n <= 1) ///not a chain
        // {
        //     bub->chain_weight.a[i].b_occ = bub->chain_weight.a[i].g_occ = 0;
        //     bub->chain_weight.a[i].del = 1;
        // } 
        // else
        {
            bub->chain_weight.a[i].del = 0;
            calculate_chain_weight(u_x, bub, ug, &(bub->chain_weight.a[i]));
        }
    }
    qsort(bub->chain_weight.a, bub->chain_weight.n, sizeof(chain_w_type), cmp_chain_weight);

    
    
    free(vis); free(is_vis); free(is_used); free(dedup); free(b.b.a); free(e_w); free(e_occ);
    kv_destroy(stack.a); kv_destroy(result.a); kv_destroy(res_utg.a); kv_destroy(edges.a);
    ma_ug_destroy(back_bs_ug);
    fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
}

void append_boundary_chain(ma_ug_t* ug, hc_links* link, bubble_type* bub)
{
    double index_time = yak_realtime();
    ma_ug_t *bs_ug = bub->b_ug;
    uint32_t v, u, i, k, beg_idx, m, nv, n_vx, flag_pri = 1, flag_aux = 2, flag_ava = 4;
    uint32_t root, rId_0, ori_0, root_0, new_bub;
    asg_arc_t *av = NULL;
    n_vx = bs_ug->g->n_seq << 1;
    uint8_t *vis = NULL; CALLOC(vis, ug->g->n_seq<<1);
    uint8_t *is_vis = NULL; CALLOC(is_vis, ug->g->n_seq<<1);
    uint8_t *is_used = NULL; CALLOC(is_used, n_vx);
    uint8_t *dedup = NULL; CALLOC(dedup, ug->g->n_seq<<1);
    buf_t b; memset(&b, 0, sizeof(buf_t));
    kvec_t_u32_warp stack, result, res_utg;
    kv_init(stack.a); kv_init(result.a); kv_init(res_utg.a);
    kvec_asg_arc_t_warp edges; kv_init(edges.a);

    for (i = 0; i < bs_ug->g->n_seq; i++)
    {
        set_b_utg_weight_flag(bub, &b, i<<1, vis, flag_aux, NULL);
    }

    if(bub->num.n > 0) bub->num.n--;
    new_bub = bub->b_g->n_seq;
    for (v = 0; v < n_vx; v++)
    {
        av = asg_arc_a(bs_ug->g, v);
        nv = asg_arc_n(bs_ug->g, v);
        if(nv == 0 || get_real_length(bs_ug->g, v, NULL) == 0) continue;
        res_utg.a.n = 0;
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            u = av[i].v^1;
            /**
            beg_idx = res_utg.a.n;
            set_b_utg_weight_flag_xor(bub, bs_ug, &b, v^1, vis, flag_pri, NULL);
            get_chain_weight(bub, bs_ug, &b, u^1, v, link, vis, is_vis, ug, &stack, &result, flag_pri, flag_ava, &res_utg, NULL);
            set_b_utg_weight_flag_xor(bub, bs_ug, &b, v^1, vis, flag_pri, NULL);
            for (k = m = beg_idx; k < res_utg.a.n; k++)
            {
                if(dedup[res_utg.a.a[k]>>1] != 0) continue;
                dedup[res_utg.a.a[k]>>1] = 1;
                res_utg.a.a[m] = res_utg.a.a[k];
                m++;
            }
            res_utg.a.n = m;
            **/
            beg_idx = res_utg.a.n;
            set_b_utg_weight_flag_xor(bub, bs_ug, &b, u^1, vis, flag_pri, NULL);
            get_chain_weight(bub, bs_ug, &b, v^1, u, link, vis, is_vis, ug, &stack, &result, flag_pri, flag_ava, &res_utg, NULL);
            set_b_utg_weight_flag_xor(bub, bs_ug, &b, u^1, vis, flag_pri, NULL);
            for (k = m = beg_idx; k < res_utg.a.n; k++)
            {
                if(dedup[res_utg.a.a[k]>>1] != 0) continue;
                dedup[res_utg.a.a[k]>>1] = 1;
                res_utg.a.a[m] = res_utg.a.a[k];
                m++;
            }
            res_utg.a.n = m;
        }

        for (k = 0; k < res_utg.a.n; k++) dedup[res_utg.a.a[k]>>1] = 0;
        /*******************************for debug************************************/
        // for (i = 0; i < res_utg.a.n; i++)
        // {
        //     for (k = 0; k < res_utg.a.n; k++)
        //     {
        //         if(k == i) continue;
        //         if((res_utg.a.a[i]>>1) == (res_utg.a.a[k]>>1)) fprintf(stderr, "ERROR\n");
        //     }
        // }
        /*******************************for debug************************************/

        root = get_utg_end_from_btg(bub, bs_ug, v);
        rId_0 = root>>1;
        ori_0 = root&1;
        get_bubbles(bub, rId_0, ori_0 == 1?&root_0:NULL, ori_0 == 0?&root_0:NULL, NULL, NULL, NULL);

        if(root_0 != (uint32_t)-1 && (!IF_HOM(root_0>>1, *bub))) kv_push(uint32_t, res_utg.a, root_0);

        if(v&1)
        {
            update_bubble_graph(&res_utg, root_0^1, rId_0, (uint32_t)-1, (uint32_t)-1, bub, &edges, bub->b_g, NULL, NULL, ug, NULL, 0);
        }
        else
        {
            update_bubble_graph(&res_utg, (uint32_t)-1, (uint32_t)-1, root_0^1, rId_0, bub, &edges, bub->b_g, NULL, NULL, ug, NULL, 0);
        }
    }
    kv_push(uint32_t, bub->num, bub->list.n);
    new_bub = bub->b_g->n_seq - new_bub;
    bub->mess_bub += new_bub;
    ///actually not useful, and may have bug when one bubble at multipe chains
    if(new_bub) update_bub_b_s_idx(bub);


    for (v = 0; v < n_vx; v++)
    {
        av = asg_arc_a(bs_ug->g, v);
        nv = asg_arc_n(bs_ug->g, v);
        if(nv == 0 || get_real_length(bs_ug->g, v, NULL) == 0) continue;
        drop_g_edges_by_utg(bub, bub->b_g, bs_ug, NULL, v, (uint32_t)-1);
    }

    update_bsg(bub->b_g, &edges);

    ma_ug_destroy(bs_ug);
    bs_ug = ma_ug_gen(bub->b_g);
    bub->b_ug = bs_ug;
    kv_destroy(bub->chain_weight);
    ma_utg_t *u_x = NULL;
    bs_ug = bub->b_ug;
    kv_malloc(bub->chain_weight, bs_ug->u.n); bub->chain_weight.n = bs_ug->u.n;
    for (i = 0; i < bs_ug->u.n; i++)
    {
        u_x = &(bs_ug->u.a[i]);
        bub->chain_weight.a[i].id = i;
        // if(u->n <= 1) ///not a chain
        // {
        //     bub->chain_weight.a[i].b_occ = bub->chain_weight.a[i].g_occ = 0;
        //     bub->chain_weight.a[i].del = 1;
        // } 
        // else
        {
            bub->chain_weight.a[i].del = 0;
            calculate_chain_weight(u_x, bub, ug, &(bub->chain_weight.a[i]));
        }
    }
    qsort(bub->chain_weight.a, bub->chain_weight.n, sizeof(chain_w_type), cmp_chain_weight);



    free(vis); free(is_vis); free(is_used); free(dedup); free(b.b.a);
    kv_destroy(stack.a); kv_destroy(result.a); kv_destroy(res_utg.a);
    kv_destroy(edges.a);


    /*******************************for debug************************************/
    // for (v = 0; v < (uint32_t)(bs_ug->g->n_seq<<1); v++)
    // {
    //     if(asg_arc_n(bs_ug->g, v) > 0) fprintf(stderr, "ERROR, btg%.6ul\n", (v>>1)+1);
    //     ma_utg_t *utg = &(bs_ug->u.a[v>>1]);
    //     for (i = 0; i < utg->n; i++)
    //     {
    //         if((utg->a[i]>>33) >= 
    //                 (bub->f_bub + bub->b_bub + bub->b_end_bub + bub->tangle_bub + bub->cross_bub))
    //         {
    //             if(i != 0 && i != utg->n - 1) fprintf(stderr, "ERROR, btg%.6ul, i: %u\n", (v>>1)+1, i);
    //         }
    //     }
    // }
    /*******************************for debug************************************/
    fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
}

int cmp_chain_hic_w_weight(const void * a, const void * b)
{
    if((*(chain_hic_w_type*)a).w != (*(chain_hic_w_type*)b).w)
    {
        return (*(chain_hic_w_type*)a).w > (*(chain_hic_w_type*)b).w? -1 : 1;
    }
    else
    {
        return 0;
    }
}


#define is_useful_bub(ID, B) (((ID)>=((B).f_bub + (B).b_bub + (B).b_end_bub + (B).tangle_bub + (B).cross_bub))\
                                   && ((ID)<((B).f_bub + (B).b_bub + (B).b_end_bub + (B).tangle_bub + (B).cross_bub + (B).mess_bub)))
void init_chain_hic_warp(ma_ug_t* ug, hc_links* link, bubble_type* bub, chain_hic_warp* c_w)
{
    ma_ug_t *bs_ug = bub->b_ug;
    uint32_t *a = NULL, n, occ, i, k_i, k_j, k_k, uID, is_del, m, bub_mess;
    double w;
    ma_utg_t *u_x = NULL;

    kv_init((*c_w)); 
    kv_malloc((*c_w), bs_ug->u.n); 
    (*c_w).n = bs_ug->u.n;
    (*c_w).max_bub_id = 0;
    (*c_w).u_n = ug->u.n;
    (*c_w).chain_idx = NULL;
    MALLOC((*c_w).chain_idx, ug->u.n); 
    memset((*c_w).chain_idx, -1, sizeof(uint32_t)*ug->u.n);

    for (i = bub_mess = 0; i < bs_ug->u.n; i++)
    {
        u_x = &(bs_ug->u.a[i]);
        (*c_w).a[i].id = i;
        (*c_w).a[i].w = 0;
        (*c_w).a[i].occ = 0;
        (*c_w).a[i].u = NULL;
        for (k_i = 0, w = 0, occ = 0; k_i < u_x->n; k_i++)
        {
            if(is_useful_bub(u_x->a[k_i]>>33, *bub))
            {
                bub_mess++;
                continue;
            }
            get_bubbles(bub, u_x->a[k_i]>>33, NULL, NULL, &a, &n, NULL);

            for (k_j = 0; k_j < n; k_j++)
            {
                uID = a[k_j]>>1;
                occ += ug->u.a[uID].n;
                for (k_k = 0; k_k < link->a.a[uID].e.n; k_k++)
                {
                    if(link->a.a[uID].e.a[k_k].del) continue;
                    w += link->a.a[uID].e.a[k_k].weight;
                }
            }
        }
        (*c_w).a[i].w = w;
        (*c_w).a[i].occ = occ;
    }
    if(bub_mess != bub->mess_bub) fprintf(stderr, "ERROR\n");
    ///fprintf(stderr, "bub_mess: %u, bub->mess_bub: %lu\n", bub_mess, bub->mess_bub);

    for (i = 0; i < bs_ug->u.n; i++)
    {
        u_x = &(bs_ug->u.a[i]);
        for (k_i = 0; k_i < u_x->n; k_i++)
        {
            if(is_useful_bub(u_x->a[k_i]>>33, *bub))
            {
                continue;
            }

            get_bubbles(bub, u_x->a[k_i]>>33, NULL, NULL, &a, &n, NULL);
            for (k_j = 0; k_j < n; k_j++)
            {
                uID = a[k_j]>>1;
                if((*c_w).chain_idx[uID] == (uint32_t)-1)
                {
                    (*c_w).chain_idx[uID] = i;
                }
                else
                {
                    if((*c_w).a[i].occ > (*c_w).a[(*c_w).chain_idx[uID]].occ)
                    {
                        (*c_w).chain_idx[uID] = i;
                    }
                }
            }
        }
    }

    for (i = m = 0; i < (*c_w).n; i++)
    {
        u_x = &(bs_ug->u.a[(*c_w).a[i].id]);
        is_del = 1;
        for (k_i = 0; k_i < u_x->n; k_i++)
        {
            if(is_useful_bub(u_x->a[k_i]>>33, *bub))
            {
                continue;
            }

            get_bubbles(bub, u_x->a[k_i]>>33, NULL, NULL, &a, &n, NULL);
            for (k_j = 0; k_j < n; k_j++)
            {
                uID = a[k_j]>>1;
                if((*c_w).chain_idx[uID] == (*c_w).a[i].id)
                {
                    is_del = 0;
                    break;
                } 
            }
            if(is_del == 0) break;
        }
        if(is_del == 0)
        {
            (*c_w).a[m] = (*c_w).a[i];
            m++;
        }
    }

    ///fprintf(stderr, "# chain: %u, # pre chain: %u\n", m, (uint32_t)(*c_w).n);
    (*c_w).n = m;
    for (i = 0; i < (*c_w).n; i++)
    {
        u_x = &(bs_ug->u.a[(*c_w).a[i].id]);
        CALLOC((*c_w).a[i].u, 1);
        for (k_i = (*c_w).a[i].u->n = 0; k_i < u_x->n; k_i++)
        {
            if(is_useful_bub(u_x->a[k_i]>>33, *bub)) continue;
            (*c_w).a[i].u->n++;
        }
        (*c_w).a[i].u->m = (*c_w).a[i].u->n;
        MALLOC((*c_w).a[i].u->a, (*c_w).a[i].u->m);
        for (k_i = (*c_w).a[i].u->n = 0; k_i < u_x->n; k_i++)
        {
            if(is_useful_bub(u_x->a[k_i]>>33, *bub)) continue;
            (*c_w).a[i].u->a[(*c_w).a[i].u->n] = u_x->a[k_i];
            (*c_w).a[i].u->n++;
        }
    }
    (*c_w).max_bub_id = (*c_w).n;
    
    if(bub->num.n > 0) bub->num.n--;
    chain_hic_w_type* p = NULL;
    for (i = 0; i < ug->u.n; i++)
    {
        uID = i;
        if(IF_HOM(uID, *bub)) continue;
        if((*c_w).chain_idx[uID] == (uint32_t)-1)
        {
            kv_pushp(chain_hic_w_type, (*c_w), &p);
            CALLOC(p->u, 1);
            p->u->n = p->u->m = 1;
            MALLOC(p->u->a, p->u->m);
            p->u->a[0] = (bub->pathLen.n)<<33;
            /********************push bubble********************/
            kv_push(uint32_t, bub->num, bub->list.n);
            kv_push(uint64_t, bub->pathLen, 0);
            kv_push(uint32_t, bub->list, uID<<1);
            kv_push(uint32_t, bub->list, uID<<1);
            kv_push(uint32_t, bub->list, uID<<1);
            /********************push bubble********************/
            p->occ = ug->u.a[uID].n;
            p->w = 0;
            for (k_k = 0; k_k < link->a.a[uID].e.n; k_k++)
            {
                if(link->a.a[uID].e.a[k_k].del) continue;
                p->w += link->a.a[uID].e.a[k_k].weight;
            }
            p->id = (*c_w).n - 1;
            (*c_w).chain_idx[uID] = p->id;
        }
    }
    kv_push(uint32_t, bub->num, bub->list.n);




    memset((*c_w).chain_idx, -1, sizeof(uint32_t)*ug->u.n);
    for (i = 0; i < (*c_w).n; i++)
    {
        (*c_w).a[i].id = i;
        u_x = (*c_w).a[i].u;
        for (k_i = 0; k_i < u_x->n; k_i++)
        {
            get_bubbles(bub, u_x->a[k_i]>>33, NULL, NULL, &a, &n, NULL);
            for (k_j = 0; k_j < n; k_j++)
            {
                uID = a[k_j]>>1;
                if((*c_w).chain_idx[uID] == (uint32_t)-1)
                {
                    (*c_w).chain_idx[uID] = (*c_w).a[i].id;
                }
                else
                {
                    if((*c_w).a[i].occ > (*c_w).a[(*c_w).chain_idx[uID]].occ)
                    {
                        (*c_w).chain_idx[uID] = (*c_w).a[i].id;
                    }
                }
            }
        }
    }

    /**
    uint32_t rId_0, ori_0, root_0, rId_1, ori_1, root_1;
    for (i = 0; i < (*c_w).n; i++)
    {
        (*c_w).a[i].l_d = (*c_w).a[i].r_d = (uint64_t)-1;
        u_x = (*c_w).a[i].u;
        if(u_x->n == 0) continue;

        rId_0 = u_x->a[0]>>33;
        ori_0 = (u_x->a[0]>>32&1)^1;
        get_bubbles(bub, rId_0, ori_0 == 1?&root_0:NULL, ori_0 == 0?&root_0:NULL, NULL, NULL, NULL);
        root_0 ^= 1;
    }
    **/

    ///fprintf(stderr, "# chain: %u, # c_w.max_bub_id: %u\n", (uint32_t)(*c_w).n, (*c_w).max_bub_id);
    ///qsort((*c_w).a, (*c_w).n, sizeof(chain_hic_w_type), cmp_chain_hic_w_weight);
}

void destory_chain_hic_warp(chain_hic_warp* c_w)
{
    uint32_t i;
    for (i = 0; i < c_w->n; i++)
    {
        free(c_w->a[i].u->a);
        free(c_w->a[i].u);
    }
    kv_destroy((*c_w)); 
    free((*c_w).chain_idx);
}

void build_bub_graph(ma_ug_t* ug, bubble_type* bub)
{
    bub->check_het = 0;
    get_bub_graph(ug, bub); ///just create nodes/edges from f_bub
    detect_bub_graph(bub, ug->g);
    asg_destroy(bub->b_g);
    bub->check_het = 1;
    get_bub_graph(ug, bub);
    ///print_bubble_chain(bub, "first round");
    // detect_bub_graph(bub, ug->g, 1);
    update_bubble_chain(ug, bub, 1, 0);
    ///print_bubble_chain(bub, "second round");
    update_bubble_chain(ug, bub, 0, 1);
    resolve_bubble_chain_tangle(ug, bub);
}

void get_forward_distance(uint32_t src, uint32_t dest, asg_t *sg, hc_links* link, MT* M)
{
    hc_edge *e = NULL;
    e = get_hc_edge(link, src, dest, 0);
    if(e == NULL) return;
    uint32_t v, j;
    uint64_t d[2], db[2], q_u, min, min_i, min_b;
    e->dis = (uint64_t)-1;

    for (v = ((uint64_t)(src)<<1); v < ((uint64_t)(src+1)<<1); v++)
    {
        d[0] = d[1] = db[0] = db[1] = (uint64_t)-1;
        for (j = 0; j < M->matrix.a[v].a.n; j++)
        {
            q_u = M->matrix.a[v].a.a[j] >> M->uID_shift;
            if((q_u>>1) == dest) d[q_u&1] = (M->matrix.a[v].a.a[j] & M->dis_mode) + sg->seq[q_u>>1].len;
            if((q_u>>1) > dest) break;///just for speeding up, doesn't affect results
        }

        min = min_i = min_b = (uint64_t)-1;
        if(e->dis != (uint64_t)-1) min = e->dis >> 3;

        if(d[0] < min) min = d[0], min_i = 0, min_b = 0;
        if(d[1] < min) min = d[1], min_i = 1, min_b = 0;
        if(min_i != (uint64_t)-1 && min != (uint64_t)-1)
        {
            e->dis = min<<1;
            e->dis += min_b;
            e->dis <<=1;
            e->dis += v&1;
            e->dis <<=1;
            e->dis += min_i;
        }
    }
    // fprintf(stderr, "%s\t%s\tdis(%lu)\n", e->dis == (uint64_t)-1? "unreach cur": "**reach cur", 
    //                                             ((e->dis>>2)&1)?"back":"forw", e->dis>>3);
}

uint32_t is_same_phase(uint64_t bid, uint64_t eid, H_partition* hap, int8_t *s, mc_gg_status *sa)
{
    if(bid == eid) return 1;
    if(hap || s)
    {
        int beg_status, end_status;
        beg_status = (hap? get_phase_status(hap, bid):s[bid]);
        if(beg_status != 1 && beg_status != -1) return (uint32_t)-1;
        end_status = (hap? get_phase_status(hap, eid):s[eid]);
        if(end_status != 1 && end_status != -1) return (uint32_t)-1;
        if(beg_status == end_status) return 1;
        return 0;
    }

    if(sa)
    {
        if(sa[bid].s == 0 || sa[eid].s == 0) return (uint32_t)-1;
        return !!(sa[bid].s&sa[eid].s);
    }
    return (uint32_t)-1;
}

int get_trans_rate_function_advance(ha_ug_index* idx, kvec_pe_hit* hits, hc_links* link, bubble_type* bub, 
H_partition* hap, int8_t *s, mc_gg_status *sa, trans_idx* dis)
{
    kvec_t(uint64_t) buf;
    kv_init(buf);
    uint64_t beg, end, cnt[2];
    uint64_t k, i, t_d, r_idx, f_idx, med = (uint64_t)-1;
    uint32_t is_s;
    // int beg_status, end_status;

    buf.n = 0;
    for (k = 0; k < hits->a.n; ++k) 
    {
        beg = ((hits->a.a[k].s<<1)>>(64 - idx->uID_bits));
        end = ((hits->a.a[k].e<<1)>>(64 - idx->uID_bits));

        if(IF_HOM(beg, *bub)) continue;
        if(IF_HOM(end, *bub)) continue;
        if(is_hom_hit(hits->a.a[k])) continue;


        t_d = get_hic_distance(&(hits->a.a[k]), link, idx, NULL);
        if(t_d == (uint64_t)-1) continue;
        // if(beg == end)
        // {
        //     t_d = (t_d << 1);
        // }
        // else
        // {
        //     beg_status = (hap? get_phase_status(hap, beg):s[beg]);
        //     if(beg_status != 1 && beg_status != -1) continue;
        //     end_status = (hap? get_phase_status(hap, end):s[end]);
        //     if(end_status != 1 && end_status != -1) continue;
        //     if(beg_status != end_status)
        //     {
        //         t_d = (t_d << 1) + 1; 
        //     }
        //     else
        //     {
        //         t_d = (t_d << 1);
        //     }            
        // }
        is_s = is_same_phase(beg, end, hap, s, sa);
        if(is_s == (uint32_t)-1) continue;
        t_d = (t_d << 1) + 1 - is_s;

        kv_push(uint64_t, buf, t_d); 
    }

    ///might have bias, we may not use right linkage larger than trans rc linkage
    radix_sort_hc64(buf.a, buf.a+buf.n);

    for (k = 0, r_idx = f_idx = (uint64_t)-1; k < buf.n; k++)
    {
        if((buf.a[k]&1) == 0) r_idx = k;
        if((buf.a[k]&1) == 1) f_idx = k;
    }
    buf.n = MIN(r_idx, f_idx);

    
    trans_p_t* p = NULL;
    dis->n = 0;
    uint64_t bin_size = MIN(2250, buf.n>>8), m;
    if(bin_size == 0)
    {
        for (i = 8; i > 0; i--)
        {
            bin_size = buf.n>>i;
            if(bin_size > 0) break;
        }
        if(bin_size == 0) bin_size = buf.n;
    }
    i = 0;
    while (i < buf.n)
    {
        kv_pushp(trans_p_t, *dis, &p);
        p->beg = i;

        cnt[0] = cnt[1] = 0;
        k = MIN(i+bin_size, buf.n);
        for (; i < k; i++)
        {
            cnt[buf.a[i]&1]++;
        }
        p->end = i;
        p->cnt_0 = cnt[0];
        p->cnt_1 = cnt[1];
    }

    i = m = 0; med = (uint64_t)-1;
    while(i < dis->n)
    {
        if(dis->a[i].cnt_0 > 0 && dis->a[i].cnt_1 > 0)
        {
            dis->a[m] = dis->a[i];
            m++;
            i++;
            continue;
        }

        if(med == (uint64_t)-1) med = buf.a[dis->a[i].beg]>>1;
        k = i; cnt[0] = cnt[1] = 0;
        for (; i < dis->n; i++)
        {
            cnt[0] += dis->a[i].cnt_0;
            cnt[1] += dis->a[i].cnt_1;
            if(cnt[0] > 0 && cnt[1] > 0) break;
        }

        if(i < dis->n)
        {
            dis->a[m].cnt_0 = cnt[0];
            dis->a[m].cnt_1 = cnt[1];
            dis->a[m].beg = dis->a[k].beg;
            dis->a[m].end = dis->a[i].end;
            m++;
            i++;
            continue;
        }

        cnt[0] -= dis->a[k].cnt_0;
        cnt[1] -= dis->a[k].cnt_1;
        while (1)
        {
            cnt[0] += dis->a[k].cnt_0;
            cnt[1] += dis->a[k].cnt_1;
            if(cnt[0] > 0 && cnt[1] > 0) break;

            if(k == 0)
            {
                k = (uint64_t)-1;
                break;
            } 
            k--;
        }

        if(k != (uint64_t)-1)
        {
            dis->a[m].cnt_0 = cnt[0];
            dis->a[m].cnt_1 = cnt[1];
            dis->a[m].end = dis->a[i-1].end;
            continue;
        }
        m = 0;
        break;
    }
    dis->n = m;
    if(dis->n == 0)
    {
        kv_destroy(buf);
        return 0;
    } 

    for (i = 0; i < dis->n; i++)
    {
        dis->a[i].beg = buf.a[dis->a[i].beg]>>1;
        dis->a[i].end = (buf.a[dis->a[i].end-1]>>1);
    }

    for (i = 0; i < dis->n - 1; i++)
    {
        dis->a[i].end += ((dis->a[i+1].beg - dis->a[i].end)/2);
        dis->a[i+1].beg = dis->a[i].end;
    }

    // for (i = 0; i < dis->n; i++)
    // {
    //     if(i > 0 && dis->a[i].beg != dis->a[i-1].end) fprintf(stderr, "ERROR: dis->a[i].beg: %lu, dis->a[i-1].end: %lu\n", dis->a[i].beg, dis->a[i-1].end);
    //     fprintf(stderr, "beg: %lu, end: %lu, cnt_0: %lu, cnt_1: %lu, error_rate: %f\n", 
    //     dis->a[i].beg, dis->a[i].end, dis->a[i].cnt_0, dis->a[i].cnt_1, (double)(dis->a[i].cnt_1)/(double)(dis->a[i].cnt_1 + dis->a[i].cnt_0));
    // }
    
    LeastSquare_advance(dis, idx, med);
    // fprintf(stderr, "idx->a: %f, idx->b: %f, idx->frac: %f, med: %lu\n", 
    // (double)idx->a, (double)idx->b, (double)idx->frac, med);

    dis->max = dis->a[dis->n-1].end;


    
    if(idx->a < 0) idx->a = 0;
    if(idx->a == 0)
    {
        idx->b = MAX((((double)(dis->a[dis->n-1].cnt_1))/((double)(dis->a[dis->n-1].cnt_0 + dis->a[dis->n-1].cnt_1))), idx->b);
    }
    if(idx->b < 0 && get_trans(idx, dis->max) < 0)
    {
        idx->b = ((double)(dis->a[dis->n-1].cnt_1))/((double)(dis->a[dis->n-1].cnt_0 + dis->a[dis->n-1].cnt_1));
    } 

    // fprintf(stderr, "idx->a: %f, idx->b: %f, idx->frac: %f, med: %lu\n", 
    // (double)idx->a, (double)idx->b, (double)idx->frac, med);





    buf.n = 0;
    for (k = 0; k < hits->a.n; ++k) 
    {
        beg = ((hits->a.a[k].s<<1)>>(64 - idx->uID_bits));
        end = ((hits->a.a[k].e<<1)>>(64 - idx->uID_bits));

        if(IF_HOM(beg, *bub)) continue;
        if(IF_HOM(end, *bub)) continue;
        if(is_hom_hit(hits->a.a[k])) continue;
        if(beg == end) continue;

        t_d = get_hic_distance(&(hits->a.a[k]), link, idx, NULL);
        if(t_d == (uint64_t)-1) continue;
        kv_push(uint64_t, buf, t_d); 
    }

    ///might have bias, we may not use right linkage larger than trans rc linkage
    radix_sort_hc64(buf.a, buf.a+buf.n);
    if(buf.n == 0)
    {
        dis->med = 0;
    }
    else
    {
        dis->med = ((buf.n&1)?buf.a[buf.n>>1]:((buf.a[buf.n>>1]+buf.a[(buf.n>>1)-1])/2));
    }

    // fprintf(stderr, "dis->med: %lu\n", dis->med);
    
    kv_destroy(buf);
    return 1;
}

void init_hic_advance(ha_ug_index* idx, kvec_pe_hit* hits, hc_links* link, bubble_type* bub, H_partition* hap, uint32_t ignore_dis)
{
    uint64_t k, i, m, is_comples_weight = 0;
    trans_idx dis;
    kv_init(dis);
    
    if(bub->round_id > 0 && ignore_dis == 0)
    {
        is_comples_weight = get_trans_rate_function_advance(idx, hits, link, bub, hap, NULL, NULL, &dis);
    }
    

    hc_edge *e = NULL;
    for (i = 0; i < link->a.n; i++)
    {
        for (k = 0; k < link->a.a[i].e.n; k++)
        {
            if(link->a.a[i].e.a[k].del) continue;
            if(link->a.a[i].e.a[k].dis == (uint64_t)-1)
            {
                e = get_hc_edge(link, link->a.a[i].e.a[k].uID, i, 0);
                if(!e) fprintf(stderr, "ERROR\n");
                e->del = link->a.a[i].e.a[k].del = 1;
            }
        }
    }

    for (i = 0; i < link->a.n; i++)
    {
        for (k = m = 0; k < link->a.a[i].e.n; k++)
        {
            if(link->a.a[i].e.a[k].del) continue;
            link->a.a[i].e.a[m] = link->a.a[i].e.a[k];
            link->a.a[i].e.a[m].weight = 0;
            link->a.a[i].e.a[m].occ = 0;
            m++;
        }
        link->a.a[i].e.n = m;
    }

    weight_edges_advance(idx, hits, link, bub, is_comples_weight == 1? &dis : NULL);

    for (i = 0; i < link->a.n; i++)
    {
        for (k = 0; k < link->a.a[i].e.n; k++)
        {
            if(link->a.a[i].e.a[k].del) continue;
            if(link->a.a[i].e.a[k].weight <= 0)
            {
                e = get_hc_edge(link, link->a.a[i].e.a[k].uID, i, 0);
                e->del = link->a.a[i].e.a[k].del = 1;
            }
        }
    }

    for (i = 0; i < link->a.n; i++)
    {
        for (k = m = 0; k < link->a.a[i].e.n; k++)
        {
            if(link->a.a[i].e.a[k].del) continue;
            link->a.a[i].e.a[m] = link->a.a[i].e.a[k];
            m++;
        }
        link->a.a[i].e.n = m;
    }
    kv_destroy(dis);
}


#define is_hap_set(i, Hap) (!!((Hap).hap[(i)]&((Hap).m[0]|(Hap).m[1]|(Hap).m[2])))
#define is_hap_set_label(i, Hap, label) (is_hap_set((i), (Hap))&&((Hap).hap[(i)]>>(Hap).label_shift)==((label)>>(Hap).label_shift))

double get_path_weight(uint32_t query, uint32_t v0, uint32_t root, bub_p_t_warp *b, hc_links* x)
{
    if(v0 == root) return 0;
    uint32_t v, u;
    hc_edge *p = NULL;
    double weight = 0;
    v = v0;
    do {
		u = b->a[v].p; // u->v
        p = get_hc_edge(x, query>>1, v>>1, 0);
        if(p) weight += p->weight;
		v = u;
	} while (v != root);

    return weight;
}

uint32_t get_related_weight(uint32_t x, H_partition* hap, double* w0, double* w1, uint32_t* hap_label)
{
    (*w0) = (*w1) = 0;
    if(x >= hap->link->a.n) return 0;
    uint32_t i, a_n = hap->link->a.a[x].e.n, occ;
    hc_edge* a = hap->link->a.a[x].e.a;
    for (i = occ = 0; i < a_n; i++)
    {
        if(a[i].del) continue;
        if(hap_label && (!is_hap_set_label(a[i].uID, *hap, *hap_label))) continue;
        if(is_hap_set(a[i].uID, *hap)) occ++;
        if((hap->hap[a[i].uID] & hap->m[0])) (*w0)+= a[i].weight;
        if((hap->hap[a[i].uID] & hap->m[1])) (*w1)+= a[i].weight;
    }
    return occ;
}

void set_path_hap(bub_p_t_warp *b, uint32_t root, H_partition* hap, uint32_t max_hap_label)
{
    uint32_t v, u, label;
    double w0 = 0, w1 = 0, cur_w0, cur_w1;
    ///v is the sink of this bubble
    v = b->S.a[0];
    do {
		u = b->a[v].p; // u->v
        if(v != b->S.a[0]) 
        {
            get_related_weight(v>>1, hap, &cur_w0, &cur_w1, &max_hap_label);
            w0 += cur_w0; w1 += cur_w1;
        }
		v = u;
	} while (v != root);

    if(w0 > w1)
    {
        label = max_hap_label | hap->m[0];
        b->exist_hap_label = hap->m[0];
    }
    else if(w0 < w1)
    {
        label = max_hap_label | hap->m[1];
        b->exist_hap_label = hap->m[1];
    }
    else
    {
        if(b->exist_hap_label == (uint32_t)-1)
        {
            label = max_hap_label | hap->m[0];
            b->exist_hap_label = hap->m[0];
        }
        else
        {
            if(b->exist_hap_label == hap->m[0])
            {
                label = max_hap_label | hap->m[1];
                b->exist_hap_label = hap->m[1];
            }
            else
            {
                label = max_hap_label | hap->m[0];
                b->exist_hap_label = hap->m[0];
            }
        }   
    }
    

	v = b->S.a[0];
    do {
		u = b->a[v].p; // u->v
        if(v != b->S.a[0]) hap->hap[v>>1] |= label;
		v = u;
	} while (v != root);
}


uint64_t get_phase_path(ma_ug_t *ug, uint32_t s, uint32_t d, bub_p_t_warp *b, H_partition* hap, uint32_t max_hap_label)
{
    asg_t *g = ug->g;
    if(g->seq[s>>1].del) return 0; // already deleted
    if(get_real_length(g, s, NULL)<2) return 0;
    uint32_t i, n_pending, is_first, to_replace, cur_nc, cur_uc, cur_ac, n_tips, tip_end, n_pop;
    double cur_nh, cur_w0, cur_w1, cur_rate, max_rate, cur_weight, max_weight;
    ///S saves nodes with all incoming edges visited
    b->S.n = b->T.n = b->b.n = b->e.n = 0;
    ///for each node, b->a saves all related information
    b->a[s].d = b->a[s].nc = b->a[s].ac = b->a[s].uc = 0; b->a[s].nh = b->a[s].w[0] = b->a[s].w[1] = 0;
    ///b->S is the nodes with all incoming edges visited
    kv_push(uint32_t, b->S, s);
    n_pop = n_tips = n_pending = 0;
    tip_end = (uint32_t)-1;
    is_first = 1;

    do {
        ///v is a node that all incoming edges have been visited
        ///d is the distance from v0 to v
        uint32_t v = kv_pop(b->S);
        uint32_t d = b->a[v].d, nc = b->a[v].nc, uc = b->a[v].uc, ac = b->a[v].ac; 
        double nh = b->a[v].nh;
        double nw_0 = b->a[v].w[0], nw_1 = b->a[v].w[1];

        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        for (i = 0; i < nv; ++i) {
            uint32_t w = av[i].v, l = (uint32_t)av[i].ul; // v->w with length l, not overlap length
            bub_p_t *t = &b->a[w];
            //got a circle
            if ((w>>1) == (s>>1)) goto pop_reset;
            //important when poping at long untig graph
            if(is_first) l = 0;
            if (av[i].del) continue;
            ///push the edge
            kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);

            if (t->s == 0) 
            { // this vertex has never been visited
                kv_push(uint32_t, b->b, w); // save it for revert
                ///t->p is the parent node of 
                ///t->s = 1 means w has been visited
                ///d is len(v0->v), l is len(v->w), so t->d is len(v0->w)
                t->p = v, t->s = 1, t->d = d + l, t->nc = nc + ug->u.a[(w>>1)].n;
                t->r = get_real_length(g, w^1, NULL);
                /**need fix**/
                t->nh = nh + get_path_weight(w, v, s, b, hap->link);

                get_related_weight(w>>1, hap, &(t->w[0]), &(t->w[1]), &max_hap_label);
                t->w[0] += nw_0; t->w[1] += nw_1; 

                t->ac = ac + ((!is_hap_set_label(w>>1, *hap, max_hap_label))?ug->u.a[(w>>1)].n : 0);
                t->uc = uc + ((is_hap_set_label(w>>1, *hap, max_hap_label))?ug->u.a[(w>>1)].n : 0);
                
                ++n_pending;
            }
            else {
                to_replace = 0;

                cur_nc = nc + ug->u.a[(w>>1)].n;
                /**need fix**/
                cur_nh = nh + get_path_weight(w, v, s, b, hap->link);
                get_related_weight(w>>1, hap, &cur_w0, &cur_w1, &max_hap_label);
                cur_w0 += nw_0; cur_w1 += nw_1;
                cur_weight = cur_nh + MAX(cur_w0, cur_w1) - MIN(cur_w0, cur_w1);
                max_weight = t->nh + MAX(t->w[0], t->w[1]) - MIN(t->w[0], t->w[1]);

                cur_ac = ac + ((!is_hap_set_label(w>>1, *hap, max_hap_label))?ug->u.a[(w>>1)].n : 0);;
                cur_uc = uc + ((is_hap_set_label(w>>1, *hap, max_hap_label))?ug->u.a[(w>>1)].n : 0);
                cur_rate = ((double)(cur_ac)/(double)(cur_ac+cur_uc));
                max_rate = ((double)(t->ac)/(double)(t->ac+t->uc));

                if(cur_rate > max_rate)
                {
                    to_replace = 1;
                }
                else if(cur_rate == max_rate)
                {
                    ///if(cur_nh > t->nh)
                    if(cur_weight > max_weight)
                    {
                        to_replace = 1;
                    }
                    else if(cur_weight == max_weight)///(cur_nh == t->nh)
                    {
                        if(cur_nc > t->nc)
                        {
                            to_replace = 1;
                        }
                        else if(cur_nc == t->nc)
                        {
                            if(d + l > t->d)
                            {
                                to_replace = 1;
                            }
                        }
                    }
                }


                if(to_replace)
                {
                    t->p = v;
                    t->nc = cur_nc;
                    t->nh = cur_nh;
                    t->ac = cur_ac;
                    t->uc = cur_uc;
                    t->w[0] = cur_w0;
                    t->w[1] = cur_w1;
                }


                if (d + l < t->d) t->d = d + l; // update dist
            }

            if (--(t->r) == 0) {
                uint32_t x = get_real_length(g, w, NULL);
                if(x > 0)
                {
                    kv_push(uint32_t, b->S, w);
                }
                else
                {
                    ///at most one tip
                    if(n_tips != 0) goto pop_reset;
                    n_tips++;
                    tip_end = w;
                }
                --n_pending;
            }
        }
        is_first = 0;


        if(n_tips == 1)
        {
            if(tip_end != (uint32_t)-1 && n_pending == 0 && b->S.n == 0)
            {
                ///sink is b.S.a[0]
                kv_push(uint32_t, b->S, tip_end);
                break;
            }
            else
            {
                goto pop_reset;
            }
        }

        if (i < nv || b->S.n == 0) goto pop_reset;
    }while (b->S.n > 1 || n_pending);


    n_pop = 1;
    /**need fix**/
    set_path_hap(b, s, hap, max_hap_label);
    pop_reset:

    for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
        bub_p_t *t = &b->a[b->b.a[i]];
        t->p = t->d = t->nc = t->ac = t->uc = t->r = t->s = 0;
        t->nh = t->w[0] = t->w[1] = 0;
    }

    return n_pop;
}

uint32_t get_weightest_hap_label_from_uid(uint64_t x, H_partition* hap, uint8_t* hap_label_flag, 
uint32_t* max_hap_label, double* max_hap_weight)
{
    (*max_hap_label) = (uint32_t)-1;
    kv_resize(double, hap->label_buffer, (hap->label>>hap->label_shift)+1);
    hap->label_buffer.n = (hap->label>>hap->label_shift)+1;
    uint32_t i, is_set_ava, is_unset_ava, a_n;
    hc_edge* a = NULL;
    for (i = 0; i < hap->label_buffer.n; i++)
    {
        hap->label_buffer.a[i] = 0;
    }

    is_set_ava = is_unset_ava = 0;
    a_n = hap->link->a.a[x].e.n;
    a = hap->link->a.a[x].e.a;
    for (i = 0; i < a_n; i++)
    {
        if(a[i].del) continue;
        if(is_hap_set(a[i].uID, *hap))
        {
            if(hap_label_flag && hap_label_flag[hap->hap[a[i].uID]] == 0) continue;
            hap->label_buffer.a[hap->hap[a[i].uID]>>hap->label_shift] += a[i].weight;
            is_set_ava = 1;
        }
        is_unset_ava = 1;
    }
    
    if(hap_label_flag && is_set_ava == 0) return 0;
    if(is_unset_ava == 0) return 0;
    if(is_set_ava == 0 && is_unset_ava > 0)
    {
        hap->label += hap->label_add;
        (*max_hap_label) = hap->label;
        return 1;
    }
    
    double max_weight;
    uint32_t max_i;
    for (i = 0, max_weight = -1, max_i = (uint32_t)-1; i < hap->label_buffer.n; i++)
    {
        if(hap->label_buffer.a[i] > max_weight)
        {
            max_weight = hap->label_buffer.a[i];
            max_i = i;
        }
    }

    (*max_hap_label) = max_i<<hap->label_shift;
    if(max_hap_weight) (*max_hap_weight) = max_weight;
    return 1;

}

uint32_t get_weightest_hap_label_from_bubble(uint64_t bid, H_partition* hap, bubble_type* bub, 
uint8_t* hap_label_flag, uint32_t* max_hap_label, double* max_hap_weight)
{
    (*max_hap_label) = (uint32_t)-1;

    kv_resize(double, hap->label_buffer, (hap->label>>hap->label_shift)+1);
    hap->label_buffer.n = (hap->label>>hap->label_shift)+1;
    uint32_t i, m, is_set_ava, is_unset_ava, a_n, *x_a, x_n, x;
    hc_edge* a = NULL;
    for (i = 0; i < hap->label_buffer.n; i++)
    {
        hap->label_buffer.a[i] = 0;
    }

    ///bid might be bubble or non-bubble
    get_bubbles(bub, bid, NULL, NULL, &x_a, &x_n, NULL);
    for (m = is_set_ava = is_unset_ava = 0; m < x_n; m++)
    {
        ///x is uid
        x = x_a[m]>>1;
        a_n = hap->link->a.a[x].e.n;
        a = hap->link->a.a[x].e.a;
        for (i = 0; i < a_n; i++)
        {
            if(a[i].del) continue;
            if(is_hap_set(a[i].uID, *hap))
            {
                if(hap_label_flag && hap_label_flag[hap->hap[a[i].uID]] == 0) continue; ///not at current chain
                hap->label_buffer.a[hap->hap[a[i].uID]>>hap->label_shift] += a[i].weight;
                is_set_ava = 1;
            }
            is_unset_ava = 1;
        }
    }
    if(hap_label_flag && is_set_ava == 0) return 0; ///no connection in current chain
    
    if(is_unset_ava == 0) return 0; ///no any connection
    if(is_set_ava == 0 && is_unset_ava > 0) ///update hap->label
    {
        hap->label += hap->label_add;
        (*max_hap_label) = hap->label;
        return 1;
    }
    
    double max_weight;
    uint32_t max_i;
    for (i = 0, max_weight = -1, max_i = (uint32_t)-1; i < hap->label_buffer.n; i++)
    {
        if(hap->label_buffer.a[i] > max_weight)
        {
            max_weight = hap->label_buffer.a[i];
            max_i = i;
        }
    }

    (*max_hap_label) = max_i<<hap->label_shift;
    if(max_hap_weight) (*max_hap_weight) = max_weight;
    return 1;
}


uint32_t get_available_com(H_partition* hap, bubble_type* bub, ma_ug_t *ug, uint32_t check_self, uint32_t check_others, 
uint8_t* hap_label_flag, uint32_t* max_hap_label)
{
    hc_links* link = hap->link;
    uint32_t beg, sink, n, *a, i, j, k, uID, max_bub_i, max_non_bub_i, max_i, is_ava;
    uint32_t hap_label, max_bub_label = (uint32_t)-1, max_non_bub_label = (uint32_t)-1;
    double w, max_bub_w, max_non_bub_w;
    max_i = (uint32_t)-1;

    for (i = 0, max_bub_w = -1, max_bub_i = (uint32_t)-1; i < bub->f_bub/**bub->s_bub**/; i++)
    {
        get_bubbles(bub, i, &beg, &sink, &a, &n, NULL);
        for (j = 0, w = 0; j < n; j++)
        {
            uID = a[j]>>1;
            if(check_self && is_hap_set(uID, *hap)) break;
        }
        if(j != n) continue;

        is_ava = 0;
        hap_label = (uint32_t)-1;
        if(check_others)
        {
            if(get_weightest_hap_label_from_bubble(i, hap, bub, 
                                                    hap_label_flag, &hap_label, &w)>0)
            {
                is_ava = 1;
            }
        }
        else
        {
            for (j = 0, w = 0, is_ava = 0; j < n; j++)
            {
                uID = a[j]>>1;
                for (k = 0; k < link->a.a[uID].e.n; k++)
                {
                    if(link->a.a[uID].e.a[k].del) continue;
                    w += link->a.a[uID].e.a[k].weight;
                    is_ava = 1;
                }
            }
        }
        
        if(is_ava == 0) continue;

        if(w > max_bub_w)
        {
            max_bub_w = w;
            max_bub_i = i;
            max_bub_label = hap_label;
        }
    }

    for (i = 0, max_non_bub_w = -1, max_non_bub_i = (uint32_t)-1; i < ug->u.n; i++)
    {
        if(IF_HET(i, *bub)/** || (bub->index[i] >= bub->s_bub && bub->index[i] < bub->f_bub)**/)
        {
            uID = i;
            is_ava = 0;
            if(check_self && is_hap_set(uID, *hap)) continue;
            hap_label = (uint32_t)-1;
            if(check_others)
            {
                if(get_weightest_hap_label_from_uid(uID, hap, hap_label_flag, &hap_label, &w)>0)
                {
                    is_ava = 1;
                }
            }
            else
            {
                for (k = 0, w = 0, is_ava = 0; k < link->a.a[uID].e.n; k++)
                {
                    if(link->a.a[uID].e.a[k].del) continue;
                    w += link->a.a[uID].e.a[k].weight;
                    is_ava = 1;
                }
            }
            
            if(is_ava == 0) continue;

            if(w > max_non_bub_w)
            {
                max_non_bub_w = w;
                max_non_bub_i = i;
                max_non_bub_label = hap_label;
            }
        }
    }

    if(max_bub_i != (uint32_t)-1 && max_non_bub_i != (uint32_t)-1)
    {
        if(max_non_bub_w > max_bub_w)
        {
            max_i = (max_non_bub_i << 1) + 1;
            w = max_non_bub_w;
            (*max_hap_label) = max_non_bub_label;
        }
        else
        {
            max_i = (max_bub_i << 1);
            w = max_bub_w;
            (*max_hap_label) = max_bub_label;
        }
    }
    else if(max_bub_i != (uint32_t)-1)
    {
        max_i = (max_bub_i << 1);
        w = max_bub_w;
        (*max_hap_label) = max_bub_label;
    }
    else if(max_non_bub_i != (uint32_t)-1)
    {
        max_i = (max_non_bub_i << 1) + 1;
        w = max_non_bub_w;
        (*max_hap_label) = max_non_bub_label;
    }

    if(max_i == (uint32_t)-1)
    {
        for (i = 0, max_non_bub_w = -1, max_non_bub_i = (uint32_t)-1; i < ug->u.n; i++)
        {
            ///if(bub->index[i] < bub->s_bub)
            if(IF_BUB(i, *bub))
            {

                uID = i;
                is_ava = 0;
                if(check_self && is_hap_set(uID, *hap)) continue;
                hap_label = (uint32_t)-1;
                if(check_others)
                {
                    if(get_weightest_hap_label_from_uid(uID, hap, hap_label_flag, &hap_label, &w)>0)
                    {
                        is_ava = 1;
                    }
                }
                else
                {
                    for (k = 0, w = 0, is_ava = 0; k < link->a.a[uID].e.n; k++)
                    {
                        if(link->a.a[uID].e.a[k].del) continue;
                        w += link->a.a[uID].e.a[k].weight;
                        is_ava = 1;
                    }
                }
                
                if(is_ava == 0) continue;

                if(w > max_non_bub_w)
                {
                    max_non_bub_w = w;
                    max_non_bub_i = i;
                    max_non_bub_label = hap_label;
                }

            }
        }

        if(max_non_bub_i != (uint32_t)-1)
        {
            max_i = (max_non_bub_i << 1) + 1;
            w = max_non_bub_w;
            (*max_hap_label) = max_non_bub_label;
        }
    }

    // if(max_i == (uint32_t)-1)
    // {
    //     fprintf(stderr, "-Cannot find!\n");
    // }
    // else if(max_i & 1)
    // {
    //     fprintf(stderr, "-utg-%uth, phasing ID: %u, w: %f, max_bub_i: %u, max_bub_w: %f, max_non_bub_i: %u, max_non_bub_w: %f\n", 
    //     max_i>>1, hap->label>>3, w, max_bub_i, max_bub_w, max_non_bub_i, max_non_bub_w);
    // }
    // else
    // {
    //     fprintf(stderr, "-bubble-%uth, phasing ID: %u, w: %f, max_bub_i: %u, max_bub_w: %f, max_non_bub_i: %u, max_non_bub_w: %f\n", 
    //     max_i>>1, hap->label>>3, w, max_bub_i, max_bub_w, max_non_bub_i, max_non_bub_w);
    // }

    return max_i;
}

void reset_ambiguous_label(H_partition* hap, uint8_t* hap_label_flag, uint32_t uID)
{
    uint32_t hap_label = (uint32_t)-1;
    if(get_weightest_hap_label_from_uid(uID, hap, hap_label_flag, &hap_label, NULL)>0)
    {
        double cur_w0, cur_w1;
        get_related_weight(uID, hap, &cur_w0, &cur_w1, &hap_label);
        if(cur_w0 >= cur_w1)
        {
            hap->hap[uID] |= (hap_label | hap->m[0]);
        }
        else
        {
            hap->hap[uID] |= (hap_label | hap->m[1]);
        }
    }
}

uint32_t get_unset_com(H_partition* hap, bubble_type* bub, ma_ug_t *ug, uint8_t* hap_label_flag, uint32_t* max_hap_label)
{
    uint32_t max_i = get_available_com(hap, bub, ug, 1, 1, hap_label_flag, max_hap_label);

    if(max_i == (uint32_t)-1)
    {
        max_i = get_available_com(hap, bub, ug, 1, 0, hap_label_flag, max_hap_label);
        if(max_i != (uint32_t)-1)
        {
            hap->label += hap->label_add;
            (*max_hap_label) = hap->label;
        }
    }

    return max_i;
}

void phase_com(H_partition* hap, ma_ug_t *ug, bub_p_t_warp* b, bubble_type* bub, uint32_t bid, uint32_t max_hap_label)
{
    
    if((bid & 1) == 0) ///bubble
    {
        uint32_t beg = (uint32_t)-1, sink = (uint32_t)-1, n, *a;
        get_bubbles(bub, bid>>1, &beg, &sink, &a, &n, NULL);
        ///fprintf(stderr, "+bubble-%uth, beg: %u, sink: %u, phasing ID: %u\n", bid>>1, beg>>1, sink>>1, hap->label>>3);
        b->exist_hap_label = (uint32_t)-1;
        get_phase_path(ug, beg, sink, b, hap, max_hap_label);
        get_phase_path(ug, beg, sink, b, hap, max_hap_label);
        ///fprintf(stderr, "-bubble-%uth, beg: %u, sink: %u, phasing ID: %u\n", bid>>1, beg>>1, sink>>1, hap->label>>3);
    }
    else
    {
        double cur_w0, cur_w1;
        get_related_weight(bid>>1, hap, &cur_w0, &cur_w1, &max_hap_label);
        ///fprintf(stderr, "utg-%uth, phasing ID: %u\n", bid>>1, hap->label>>3);
        if(cur_w0 >= cur_w1)
        {
            hap->hap[bid>>1] |= (max_hap_label | hap->m[0]);
        }
        else
        {
            hap->hap[bid>>1] |= (max_hap_label | hap->m[1]);
        }
    }
}


double get_cluster_weight(H_partition* hap, hc_links* link, uint32_t *h, uint32_t h_n)
{
    int o_d = 0;
    double weight = 0;
    uint32_t j, k, m;
    for (j = 0, weight = 0; j < h_n; j++)
    {
        for (k = 0; k < link->a.a[h[j]].e.n; k++)
        {
            if(link->a.a[h[j]].e.a[k].del) continue;
            for (m = 0; m < h_n; m++)
            {
                if(h[m] == link->a.a[h[j]].e.a[k].uID) break;
            }
            if(m < h_n) continue;

            o_d = get_phase_status(hap, link->a.a[h[j]].e.a[k].uID);
            if(o_d < -1) continue;
            ///if(o_d < -1) fprintf(stderr, "ERROR\n");
            weight += (o_d*link->a.a[h[j]].e.a[k].weight);
        }
    }

    return weight;
}


double get_cluster_inner_weight(H_partition* hap, hc_links* link, uint32_t *h0, uint32_t h0_n, 
uint32_t *h1, uint32_t h1_n)
{
    double weight = 0;
    uint32_t j, k, m;
    for (j = 0, weight = 0; j < h0_n; j++)
    {
        for (k = 0; k < link->a.a[h0[j]].e.n; k++)
        {
            if(link->a.a[h0[j]].e.a[k].del) continue;
            for (m = 0; m < h1_n; m++)
            {
                if(h1[m] == link->a.a[h0[j]].e.a[k].uID) break;
            }
            if(m == h1_n) continue;

            weight += link->a.a[h0[j]].e.a[k].weight;
        }
    }

    return weight * 2;
}

void update_partition_flag(H_partition* h, G_partition* g_p, hc_links* link, uint32_t id)
{
    uint32_t k, *h0, h0_n, *h1, h1_n, uID, flag = 0;
    int status;
    get_phased_block(g_p, NULL, id, NULL, NULL, &h0, &h0_n, &h1, &h1_n, NULL, NULL);

    status = g_p->a[id].status[0];
    if(status == 1) flag = h->m[0];
    else if(status == -1) flag = h->m[1];
    else if(status == 0) flag = h->m[2];
    else if(status == -2) flag = 0;
    for (k = 0; k < h0_n; k++)
    {
        uID = h0[k];
        h->hap[uID] >>= 3;
        h->hap[uID] <<= 3;
        h->hap[uID] |= flag;
    }

    status = g_p->a[id].status[1];
    if(status == 1) flag = h->m[0];
    else if(status == -1) flag = h->m[1];
    else if(status == 0) flag = h->m[2];
    else if(status == -2) flag = 0;
    for (k = 0; k < h1_n; k++)
    {
        uID = h1[k];
        h->hap[uID] >>= 3;
        h->hap[uID] <<= 3;
        h->hap[uID] |= flag;
    }
}


void update_partition_flag_debug(H_partition* h, G_partition* g_p, hc_links* link, uint32_t id)
{
    uint32_t k, *h0, h0_n, *h1, h1_n, uID, flag = 0;
    int status;
    get_phased_block(g_p, NULL, id, NULL, NULL, &h0, &h0_n, &h1, &h1_n, NULL, NULL);

    status = g_p->a[id].status[0];
    if(status == 1) flag = h->m[0];
    else if(status == -1) flag = h->m[1];
    else if(status == 0) flag = h->m[2];
    else if(status == -2) flag = 0;
    for (k = 0; k < h0_n; k++)
    {
        uID = h0[k];
        if(flag != (h->hap[uID]&7)) fprintf(stderr, "h0, id: %u, uID: %u, pre_flag: %u, cur_flag: %u\n", id, uID, (h->hap[uID]&7), flag);
        h->hap[uID] >>= 3;
        h->hap[uID] <<= 3;
        h->hap[uID] |= flag;
    }

    status = g_p->a[id].status[1];
    if(status == 1) flag = h->m[0];
    else if(status == -1) flag = h->m[1];
    else if(status == 0) flag = h->m[2];
    else if(status == -2) flag = 0;
    for (k = 0; k < h1_n; k++)
    {
        uID = h1[k];
        if(flag != (h->hap[uID]&7)) fprintf(stderr, "h1, id: %u, uID: %u, pre_flag: %u, cur_flag: %u\n", id, uID, (h->hap[uID]&7), flag);
        h->hap[uID] >>= 3;
        h->hap[uID] <<= 3;
        h->hap[uID] |= flag;
    }
}

void print_contig_partition(H_partition* hap, const char* debug)
{
    uint32_t i;
    int status;
    for (i = 0; i < hap->n; i++)
    {
        status = get_phase_status(hap, i);
        fprintf(stderr, "%s\tutg%.6d\tP:%u\tHG:A:%d\n", debug, (int)(i+1), hap->hap[i]>>3, status);
    }
}

void adjust_contig_partition(H_partition* hap, hc_links* link)
{
    double index_time = yak_realtime();
    uint32_t i, k, *h0, h0_n, *h1, h1_n;
    uint32_t h0_status[4], h1_status[4], h0_status_max;
    int h0_h, h1_h;
    for (i = 0; i < hap->g_p.n; i++)
    {

        get_phased_block(&(hap->g_p), NULL, i, NULL, NULL, &h0, &h0_n, &h1, &h1_n, NULL, NULL);
        hap->g_p.a[i].status[0] = hap->g_p.a[i].status[1] = -2;
        hap->g_p.a[i].weight[0] = hap->g_p.a[i].weight[1] = hap->g_p.a[i].weight_convex = 0;


        h0_status[0] = h0_status[1] = h0_status[2] = h0_status[3] = 0;
        for (k = 0; k < h0_n; k++)
        {
            h0_status[get_phase_status(hap, h0[k])+2]++;
        }

        hap->g_p.a[i].weight[0] = get_cluster_weight(hap, link, h0, h0_n);
        if(h1_n == 0)
        {
            if(h0_status[0] == h0_n) ///unset, flag = -2
            {
                hap->g_p.a[i].status[0] = -2;
            }
            else
            {
                if(h0_status[1] > 0 || h0_status[3] > 0) ///phased flag = 1/-1
                {
                    h0_status[0] = h0_status[2] = 0;

                    h0_h = -2;
                    h0_status_max = 0;
                    for (k = 0; k < 4; k++)
                    {
                        if(h0_status[k] > h0_status_max)
                        {
                            h0_status_max = h0_status[k];
                            h0_h = (int)(k) - 2;
                        }
                    }
                    hap->g_p.a[i].status[0] = h0_h;
                }
                else if(h0_status[2] > 0) ///hom flag
                {
                    hap->g_p.a[i].status[0] = 0;
                }
                else  //unset flag
                {
                    hap->g_p.a[i].status[0] = -2;
                }
            }
        }
        else
        {
            hap->g_p.a[i].weight[1] = get_cluster_weight( hap, link, h1, h1_n);
            h1_status[0] = h1_status[1] = h1_status[2] = h1_status[3] = 0;
            for (k = 0; k < h1_n; k++)
            {
                h1_status[get_phase_status(hap, h1[k])+2]++;
            }


            h0_h = h1_h = 0;
            for (k = 0; k < 4; k++)
            {
                if(h0_status[k] == h0_n) h0_h = (int)(k) - 2;
                if(h1_status[k] == h1_n) h1_h = (int)(k) - 2;
            }

            if(h0_h * h1_h == -1)
            {
                hap->g_p.a[i].status[0] = h0_h;
                hap->g_p.a[i].status[1] = h1_h;
            }
            else
            {
                if(hap->g_p.a[i].weight[0] >= hap->g_p.a[i].weight[1])
                {
                    hap->g_p.a[i].status[0] = 1;
                    hap->g_p.a[i].status[1] = -1;
                }
                else
                {
                    hap->g_p.a[i].status[0] = -1;
                    hap->g_p.a[i].status[1] = 1;
                }
            }
        }
        hap->g_p.a[i].weight_convex = get_cluster_inner_weight(hap, link, h0, h0_n, h1, h1_n);
        update_partition_flag(hap, &(hap->g_p), link, i);
        ///update_partition_flag_debug(hap, &(hap->g_p), link, i);
    }

    for (i = 0; i < hap->g_p.n; i++)
    {
        get_phased_block(&(hap->g_p), NULL, i, NULL, NULL, &h0, &h0_n, &h1, &h1_n, NULL, NULL);
        hap->g_p.a[i].weight[0] = get_cluster_weight(hap, link, h0, h0_n);
        hap->g_p.a[i].weight[1] = get_cluster_weight(hap, link, h1, h1_n);
        hap->g_p.a[i].weight_convex = get_cluster_inner_weight(hap, link, h0, h0_n, h1, h1_n);
    }
    fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
}

uint32_t get_weightest_uid(uint32_t* a, uint32_t n, H_partition* hap, uint8_t* hap_label_flag, 
uint32_t* max_hap_label)
{
    double cur_w0, cur_w1, max_weight;
    uint32_t j, max_i, is_ava, hap_label;

    for (j = is_ava = 0, max_i = (uint32_t)-1; j < n; j++)
    {
        if(is_hap_set(a[j]>>1, *hap)) ///avoid repeat phasing
        {
            continue;
        }
        
        get_weightest_hap_label_from_uid(a[j]>>1, hap, hap_label_flag, &hap_label, NULL);

        if(hap_label == (uint32_t)-1)
        {
            continue;
        }

        if(get_related_weight(a[j]>>1, hap, &cur_w0, &cur_w1, &hap_label) > 0)
        {
            if(max_i == (uint32_t)-1)
            {
                max_i = j;
                max_weight = cur_w0 + cur_w1;
                (*max_hap_label) = hap_label;
            }
            else if((cur_w0 + cur_w1) > max_weight)
            {
                max_i = j;
                max_weight = cur_w0 + cur_w1;
                (*max_hap_label) = hap_label;
            }
            is_ava++;
        }
    }

    ///no useful unitig
    if(is_ava == 0)
    {
        for (j = is_ava = 0, max_i = (uint32_t)-1; j < n; j++)
        {
            if(is_hap_set(a[j]>>1, *hap)) 
            {
                continue;
            }
            
            if(get_weightest_hap_label_from_uid(a[j]>>1, hap, NULL, &hap_label, NULL)==0) continue;

            if(get_related_weight(a[j]>>1, hap, &cur_w0, &cur_w1, &hap_label) > 0)
            {
                if(max_i == (uint32_t)-1)
                {
                    max_i = j;
                    max_weight = cur_w0 + cur_w1;
                    (*max_hap_label) = hap_label;
                }
                else if((cur_w0 + cur_w1) > max_weight)
                {
                    max_i = j;
                    max_weight = cur_w0 + cur_w1;
                    (*max_hap_label) = hap_label;
                }
                is_ava++;
            }
        }
    }
    if(max_i == (uint32_t)-1) return max_i;

    return a[max_i]>>1;
}


void get_weightest_hap_label_from_chain(ma_utg_t *u, H_partition* hap, bubble_type* bub, 
uint32_t* max_hap_label, uint32_t* max_bid_idx, uint32_t* is_forward_first)
{
    (*max_hap_label) = (uint32_t)-1;
    if(max_bid_idx) (*max_bid_idx) = 0;
    if(is_forward_first) (*is_forward_first) = 1;

    if(u->n == 0) return;
    kv_resize(double, hap->label_buffer, (hap->label>>hap->label_shift)+1);///how many group
    hap->label_buffer.n = (hap->label>>hap->label_shift)+1;
    uint32_t i, k, m, is_ava, a_n, *x_a, x_n, x;
    uint64_t bid;
    hc_edge* a = NULL;
    for (i = 0; i < hap->label_buffer.n; i++)
    {
        hap->label_buffer.a[i] = 0;///count weight for each group
    }

    for (k = is_ava = 0; k < u->n; k++)
    {
        bid = u->a[k]>>33;
        ///bid might be bubble or non-bubble
        get_bubbles(bub, bid, NULL, NULL, &x_a, &x_n, NULL);
        for (m = 0; m < x_n; m++)
        {
            ///x is uid
            x = x_a[m]>>1;
            a_n = hap->link->a.a[x].e.n;
            a = hap->link->a.a[x].e.a;
            for (i = 0; i < a_n; i++)
            {
                if(a[i].del) continue;
                if(is_hap_set(a[i].uID, *hap))
                {
                    hap->label_buffer.a[hap->hap[a[i].uID]>>hap->label_shift] += a[i].weight;
                    is_ava = 1;
                }
            }
        }
    }
    
    if(is_ava == 0) return; ///totally new chain
    double max_weight;
    uint32_t max_i;
    ///select the best exsiting group to u
    for (i = 0, max_weight = -1, max_i = (uint32_t)-1; i < hap->label_buffer.n; i++)
    {
        if(hap->label_buffer.a[i] > max_weight)
        {
            max_weight = hap->label_buffer.a[i];
            max_i = i;
        }
    }

    (*max_hap_label) = max_i<<hap->label_shift;
    if(max_bid_idx)
    {
        double current_weight, tot_w = 0, half_w = 0;
        for (k = 0, max_weight = -1, max_i = (uint32_t)-1; k < u->n; k++)
        {
            bid = u->a[k]>>33;
            ///bid might be bubble or non-bubble
            get_bubbles(bub, bid, NULL, NULL, &x_a, &x_n, NULL);
            for (m = 0, current_weight = 0; m < x_n; m++)
            {
                x = x_a[m]>>1;
                a_n = hap->link->a.a[x].e.n;
                a = hap->link->a.a[x].e.a;
                for (i = 0; i < a_n; i++)
                {
                    if(a[i].del) continue;
                    if(is_hap_set_label(a[i].uID, *hap, *max_hap_label))
                    {
                        current_weight += a[i].weight;
                    }
                }
            }

            if(current_weight > max_weight) /// the weight of each bubble
            {
                max_weight = current_weight;
                max_i = k;
                half_w = 0;
            }
            tot_w += current_weight;
            half_w += current_weight;
        }

        (*max_bid_idx) = max_i;
        if(is_forward_first)
        {
            if(half_w >= (tot_w - half_w))
            {   
                (*is_forward_first) = 1;
            }
            else
            {
                (*is_forward_first) = 0;
            }
        }
    }

    return;
}

void phase_bubble_chain_dir(H_partition* hap, ma_ug_t *ug, bub_p_t_warp* b, bubble_type* bub, 
ma_utg_t *u, uint8_t* hap_label_flag, uint32_t beg_idx, uint32_t end_idx, uint32_t is_forward)
{
    uint32_t i, j, beg, sink, *a, n, max_hap_label, max_uid;
    uint64_t bid;
    double cur_w0, cur_w1;
    for (i = beg_idx; i <= end_idx; i++)
    {
        bid = is_forward? u->a[i]>>33:u->a[u->n-i-1]>>33;
        get_bubbles(bub, bid, &beg, &sink, &a, &n, NULL);
        
        if(bub->b_g->seq[bid].c != HAP_LABLE/** && bid < bub->s_bub**/) ///simple bubble
        {   
            for (j = 0; j < n; j++) ///avoiding repeat phasing
            {
                if(is_hap_set(a[j]>>1, *hap)) break;
            }

            if(j < n) goto complete;

            ///get current max_hap_label from current chain
            get_weightest_hap_label_from_bubble(bid, hap, bub, hap_label_flag, &max_hap_label, NULL);

            
            ///three levels: 
            ///1. has setted weight with same hap label (using current max_hap_label)
            ///2. has setted weight but with different hap labels (using max max_hap_label from current weight)
            ///3. has unsetted weight (add hap->hap_label)
            ///4. skip, do nothing
            if(max_hap_label == (uint32_t)-1 && 
                    get_weightest_hap_label_from_bubble(bid, hap, bub, NULL, &max_hap_label, NULL) == 0) 
            {
                goto complete;
            }

            ///phase bubble
            b->exist_hap_label = (uint32_t)-1;
            get_phase_path(ug, beg, sink, b, hap, max_hap_label);
            get_phase_path(ug, beg, sink, b, hap, max_hap_label);
            for (j = 0; j < n; j++) ///set bubble as visited
            {
                hap_label_flag[a[j]>>1] = 1;
            }
        }
        else
        {
            ///select unitig with highest related weight at one time
            while (1)
            {
                max_uid = get_weightest_uid(a, n, hap, hap_label_flag, &max_hap_label);
                if(max_uid == (uint32_t)-1) break;

                get_related_weight(max_uid, hap, &cur_w0, &cur_w1, &max_hap_label);
                if(cur_w0 >= cur_w1)
                {
                    hap->hap[max_uid] |= (max_hap_label | hap->m[0]);
                }
                else
                {
                    hap->hap[max_uid] |= (max_hap_label | hap->m[1]);
                }

                hap_label_flag[max_uid] = 1;
            }
        }

        complete:;
    }
}
///ignore unitigs wihci have already been labeled in current chain (might happen)
void phase_bubble_chain(H_partition* hap, ma_ug_t *ug, bub_p_t_warp* b, bubble_type* bub, 
                                                    uint8_t* hap_label_flag, uint32_t chain_id)
{
    ma_utg_t *u = &(bub->b_ug->u.a[chain_id]);
    if(u->n == 0) return;
    uint32_t is_forward = 1, i, max_hap_label, max_bid_idx;
    memset(hap_label_flag, 0, ug->g->n_seq);

    get_weightest_hap_label_from_chain(u, hap, bub, &max_hap_label, &max_bid_idx, &is_forward);

    if(max_hap_label != (uint32_t)-1) ///means this is not a new chain
    {
        for (i = 0; i < hap->n; i++)
        {
            if(is_hap_set_label(i, *hap, max_hap_label)) hap_label_flag[i] = 1;
        }
    }

    if(max_hap_label == (uint32_t)-1) ///a totally new chain
    {   
        phase_bubble_chain_dir(hap, ug, b, bub, u, hap_label_flag, 0, u->n - 1, 1);
    }
    else
    {
        if(max_bid_idx == 0)
        {
            phase_bubble_chain_dir(hap, ug, b, bub, u, hap_label_flag, 0, u->n - 1, 1);
        }
        else
        {
            if(is_forward)
            {
                phase_bubble_chain_dir(hap, ug, b, bub, u, hap_label_flag, max_bid_idx, u->n - 1, 1);
                phase_bubble_chain_dir(hap, ug, b, bub, u, hap_label_flag, 0, max_bid_idx-1, 0);
            }
            else
            {
                phase_bubble_chain_dir(hap, ug, b, bub, u, hap_label_flag, 0, max_bid_idx-1, 0);
                phase_bubble_chain_dir(hap, ug, b, bub, u, hap_label_flag, max_bid_idx, u->n - 1, 1);
            }
        }
    }
}


uint32_t if_flip(H_partition* h, G_partition* g_p, hc_links* link,  
bubble_type* bub, uint32_t gid)
{
    double weight = 0;
    if(h->lock[gid]) return 0;
    if(g_p->a[gid].h[0] > 0 && (g_p->a[gid].status[0] == 1 || g_p->a[gid].status[0] == -1))
    {
        weight += (g_p->a[gid].weight[0] * g_p->a[gid].status[0]);
    }

    if(g_p->a[gid].h[1] > 0 && (g_p->a[gid].status[1] == 1 || g_p->a[gid].status[1] == -1))
    {
        weight += (g_p->a[gid].weight[1] * g_p->a[gid].status[1]);
    }
    weight += g_p->a[gid].weight_convex*2;
    if(weight >= 0) return 0;
    return 1;
}
void flip_unitig(G_partition* g_p, hc_links* link, bubble_type* bub, uint32_t id);
uint32_t phasing_improvement(H_partition* h, G_partition* g_p, ha_ug_index* idx, bubble_type* bub, hc_links* link);
uint32_t get_max_unitig(H_partition* h, G_partition* g_p, hc_links* link, bubble_type* bub);
double get_cluster_weight_debug(G_partition* g_p, hc_links* link, uint32_t *h, uint32_t h_n);

void debug_flip(G_partition* g_p, hc_links* link, bubble_type* bub, uint32_t id)
{
    fprintf(stderr, "33333333333\n");
    uint32_t *h0, h0_n, *h1, h1_n, k, wrong;
    double hw0, hw1;
    for (k = wrong = 0; k < g_p->n; k++)
    {
        hw0 = hw1 = 0;
        get_phased_block(g_p, NULL, k, NULL, NULL, &h0, &h0_n, &h1, &h1_n, NULL, NULL);
        if(h0_n >0) hw0 = get_cluster_weight_debug(g_p, link, h0, h0_n);
        if(h1_n >0) hw1 = get_cluster_weight_debug(g_p, link, h1, h1_n);
        if(hw0 != g_p->a[k].weight[0])
        {
            if((uint32_t)hw0 != (uint32_t)g_p->a[k].weight[0]) wrong = 1;
            fprintf(stderr, "k: %u, ERROR(id: %u): hw0: %f, weight[0]: %f\n", k, id, hw0, g_p->a[k].weight[0]);
        } 
        
        if(hw1 != g_p->a[k].weight[1])
        {
            if((uint32_t)hw1 != (uint32_t)g_p->a[k].weight[1]) wrong = 1;
            fprintf(stderr, "k: %u, ERROR(id: %u): hw1: %f, weight[1]: %f\n", k, id, hw1, g_p->a[k].weight[1]);
        }

        if(wrong) break;
    }
}
void merge_phase_group_by_chain(H_partition* hap, G_partition* g_p, bubble_type* bub, uint32_t chain_id)
{
    uint32_t i, k;
    uint32_t beg, sink, *a, n, pre_id, hap_label_id;
    uint64_t bid, uid;
    ma_utg_t *u = &(bub->b_ug->u.a[chain_id]);
    
    for (i = 0, pre_id = (uint32_t)-1; i < u->n; i++)
    {
        // fprintf(stderr, "inner i: %u, u->n: %u\n", i, (uint32_t)u->n);
        bid = u->a[i]>>33; ///here is a bubble
        // fprintf(stderr, "bid: %lu\n", bid);
        get_bubbles(bub, bid, &beg, &sink, &a, &n, NULL);
        for (k = 0; k < n; k++)
        {
            // fprintf(stderr, "k: %u, n: %u\n", k, n);
            uid = a[k]>>1;
            // fprintf(stderr, "uid: %lu\n", uid);
            if(g_p->index[uid] == (uint32_t)-1) ///mean this unitig doesn't have hap label
            {
                pre_id = (uint32_t)-1;
                continue;
            }
            hap_label_id = g_p->index[uid]>>1;
            // fprintf(stderr, "hap_label_id: %u, hap->n: %lu\n", hap_label_id, hap->n);
            if(hap_label_id == pre_id) continue;
            pre_id = hap_label_id;
            if(hap->lock[hap_label_id] == 1) continue;
            if(if_flip(hap, g_p, hap->link, bub, hap_label_id))
            {
                // fprintf(stderr, "2222222222\n");
                flip_unitig(g_p, hap->link, bub, hap_label_id);
                ///debug_flip(g_p, hap->link, bub, hap_label_id);
            }
        }
    }
}

void print_phase_group(G_partition* g_p, bubble_type* bub, const char* command)
{
    uint32_t i, k;
    partition_warp *res = NULL;
    for (i = 0; i < g_p->n; i++)
    {
        res = &(g_p->a[i]);
        fprintf(stderr, "\n%s: %u-th group: # %d = %u (weight: %f), # %d = %u (weight: %f), inner_weight: %f\n", command, i,
        res->status[0], res->h[0], res->weight[0], res->status[1], res->h[1], res->weight[1], res->weight_convex);
        for (k = 0; k < res->h[0]; k++)
        {
            fprintf(stderr, "%d: utg%.6ul\n", res->status[0], int(res->a.a[k]+1));
        }

        for (; k < res->a.n; k++)
        {
            fprintf(stderr, "%d: utg%.6ul\n", res->status[1], int(res->a.a[k]+1));
        }
    }
}

void set_bubble(H_partition* hap, G_partition* g_p, bubble_type* bub, block_phase_type* block, 
uint64_t bid)
{
    uint32_t beg, sink, k, uid, *a, n, gid;
    block->weight = 0;
    get_bubbles(bub, bid, &beg, &sink, &a, &n, NULL);
    for (k = 0; k < n; k++)
    {
        uid = a[k]>>1;
        if(g_p->index[uid] == (uint32_t)-1) continue;
        gid = g_p->index[uid]>>1;
        block->vis.a[gid] = 1;
    }
}

uint32_t next_hap_label_id(block_phase_type* b, G_partition* g_p, bubble_type* bub, ma_utg_t *u, 
int is_forward, long long* c_bid, long long* c_uid)
{
    uint32_t beg, sink, uid, *a, n, gid, pre_gid;
    while (1) ///while(b->bid < (long long)u->n)
    {
        if(is_forward == 1 && b->bid >= (long long)u->n) break;
        if(is_forward == 0 && b->bid < 0) break;

        get_bubbles(bub, u->a[b->bid]>>33, &beg, &sink, &a, &n, NULL);
        while (1) ///while (b->uid < (long long)n)
        {
            if(is_forward == 1 && b->uid >= (long long)n) break;
            if(is_forward == 0 && b->uid < 0) break;

            uid = a[b->uid]>>1;
            gid = (uint32_t)-1;
            if(g_p->index[uid] != (uint32_t)-1)
            {
                gid = g_p->index[uid]>>1;
            }
            if(c_bid) (*c_bid) = b->bid;
            if(c_uid) (*c_uid) = b->uid;

            if(is_forward == 1) b->uid++;
            if(is_forward == 0) b->uid--;
            if(gid == (uint32_t)-1) continue;


            pre_gid = uid = (uint32_t)-1;
            if(is_forward == 1 && (b->uid >= 2)) uid = a[b->uid - 2]>>1;
            if(is_forward == 0 && (b->uid + 2 < n)) uid = a[b->uid + 2]>>1;
            if(uid != (uint32_t)-1 && g_p->index[uid] != (uint32_t)-1) pre_gid = g_p->index[uid]>>1;

            if(pre_gid == gid) continue;

            return gid;
        }
        
        if(is_forward == 1) b->bid++, b->uid = 0;
        if(is_forward == 0) b->bid--, b->uid = (long long)n - (long long)1;
    }

    return (uint32_t)-1;
}

double get_new_weight(G_partition* g_p, uint8_t* flag, hc_links* link, uint32_t gid)
{
    double total_weight = 0, weight;
    int status, o_d;
    uint32_t *h0, h0_n, *h1, h1_n, *h, h_n, j, k, uID;

    if(g_p->a[gid].h[0] > 0 && (g_p->a[gid].status[0] == 1 || g_p->a[gid].status[0] == -1))
    {
        total_weight += (g_p->a[gid].weight[0] * g_p->a[gid].status[0]);
    }

    if(g_p->a[gid].h[1] > 0 && (g_p->a[gid].status[1] == 1 || g_p->a[gid].status[1] == -1))
    {
        total_weight += (g_p->a[gid].weight[1] * g_p->a[gid].status[1]);
    }
    total_weight += g_p->a[gid].weight_convex*2;



    get_phased_block(g_p, NULL, gid, NULL, NULL, &h0, &h0_n, &h1, &h1_n, NULL, NULL);

    h = h0; h_n = h0_n; status = g_p->a[gid].status[0];
    for (j = 0, weight = 0; j < h_n; j++)
    {
        for (k = 0; k < link->a.a[h[j]].e.n; k++)
        {
            if(link->a.a[h[j]].e.a[k].del) continue;
            uID = link->a.a[h[j]].e.a[k].uID;
            if(g_p->index[uID] == (uint32_t)-1) continue;
            if(flag[g_p->index[uID]>>1] == 0) continue;
            if((g_p->index[uID]>>1) == gid) continue;
            o_d = g_p->a[g_p->index[uID]>>1].status[g_p->index[uID]&1];
            weight += (status*o_d*link->a.a[h[j]].e.a[k].weight);
        }
    }
    total_weight -= (2*weight);




    h = h1; h_n = h1_n; status = g_p->a[gid].status[1];
    for (j = 0, weight = 0; j < h_n; j++)
    {
        for (k = 0; k < link->a.a[h[j]].e.n; k++)
        {
            if(link->a.a[h[j]].e.a[k].del) continue;
            uID = link->a.a[h[j]].e.a[k].uID;
            if(g_p->index[uID] == (uint32_t)-1) continue;
            if(flag[g_p->index[uID]>>1] == 0) continue;
            if((g_p->index[uID]>>1) == gid) continue;
            o_d = g_p->a[g_p->index[uID]>>1].status[g_p->index[uID]&1];
            weight += (status*o_d*link->a.a[h[j]].e.a[k].weight);
        }
    }
    total_weight -= (2*weight);



    return total_weight;
}

int identify_best_interval(block_phase_type* i_buf, uint8_t* lock, G_partition* g_p, bubble_type* bub, 
ma_utg_t *u, hc_links* link, long long f_bid, long long f_uid, long long* l_bid, long long* l_uid)
{
    long long c_bid, c_uid, min_bid, min_uid;
    double w = 0, min_w = 1;
    uint32_t gid, val = 0;
    block_phase_type b;
    b.bid = f_bid; b.uid = f_uid;
    memset(i_buf->vis.a, 0, g_p->n);
    i_buf->weight = min_w = 1;  min_bid = min_uid = -1;
    (*l_bid) = (*l_uid) = -1;

    while (1)
    {
        gid = next_hap_label_id(&b, g_p, bub, u, 1, &c_bid, &c_uid);
        if(gid == (uint32_t)-1) break;
        if(i_buf->vis.a[gid] == 1) continue;
        w += get_new_weight(g_p, i_buf->vis.a, link, gid);
        i_buf->vis.a[gid] = 1;
        if(lock[gid] == 0) val = 1;
        if(val == 0) continue;

        if(min_w > w)
        {
            min_w = w;
            min_bid = c_bid;
            min_uid = c_uid;
        }
    }

    if(min_w < 0 && min_bid != -1 && min_uid != -1)
    {
        (*l_bid) = min_bid; (*l_uid) = min_uid;
        i_buf->weight = min_w;
        ///fprintf(stderr, "+min_w: %f, min_bid: %lld, min_uid: %lld\n", min_w, min_bid, min_uid);
    }

    if(val == 0) return 1;

    return 0;
}

void identify_best_interval_debug(block_phase_type* i_buf, G_partition* g_p, bubble_type* bub, 
ma_utg_t *u, hc_links* link, long long f_bid, long long f_uid, long long* l_bid, long long* l_uid)
{
    long long c_bid, c_uid, min_bid, min_uid;
    double w = 0, min_w = 1;
    uint32_t gid;
    block_phase_type b;
    b.bid = f_bid; b.uid = f_uid;
    memset(i_buf->vis.a, 0, g_p->n);
    min_w = 1;  min_bid = min_uid = -1;
    (*l_bid) = (*l_uid) = -1;

    while (1)
    {
        gid = next_hap_label_id(&b, g_p, bub, u, 1, &c_bid, &c_uid);
        if(gid == (uint32_t)-1) break;
        if(i_buf->vis.a[gid] == 1) continue;
        w += get_new_weight(g_p, i_buf->vis.a, link, gid);
        if(f_bid == 689 && f_uid == 1)
        {
            fprintf(stderr, "gid: %u, w: %f\n", gid, w);
        }
        
        if(min_w > w)
        {
            min_w = w;
            min_bid = c_bid;
            min_uid = c_uid;
        }
        i_buf->vis.a[gid] = 1;
    }

    if(min_w < 0 && min_bid != -1 && min_uid != -1)
    {
        (*l_bid) = min_bid; (*l_uid) = min_uid;
        i_buf->weight = min_w;
        b.bid = min_bid; b.uid = min_uid;


        gid = next_hap_label_id(&b, g_p, bub, u, 1, &c_bid, &c_uid);
        fprintf(stderr, "-min_w: %f, f_bid: %lld, f_uid: %lld, min_bid: %lld, min_uid: %lld, gid: %u\n", 
                                                        min_w, f_bid, f_uid, min_bid, min_uid, gid);
        
    }
}

int flip_block(block_phase_type* i_buf, G_partition* g_p, bubble_type* bub, 
ma_utg_t *u, hc_links* link, uint8_t* lock, long long f_bid, long long f_uid, long long l_bid, long long l_uid)
{
    long long c_bid = l_bid, c_uid = l_uid;
    uint32_t gid, val = 0;
    block_phase_type b;
    b.bid = f_bid; b.uid = f_uid;
    memset(i_buf->vis.a, 0, g_p->n);
    val = 0;
    while (1)
    {
        gid = next_hap_label_id(&b, g_p, bub, u, 1, &c_bid, &c_uid);
        ////fprintf(stderr, "+gid: %u\n", gid);
        if(gid == (uint32_t)-1) break;
        ///fprintf(stderr, "lock[gid]: %u\n", lock[gid]);
        if(lock[gid] == 0)
        {
            val = 1;
            break;
        } 
        if(c_bid == l_bid && c_uid == l_uid) break;
    }

    if(val == 0) return 0;

    b.bid = f_bid; b.uid = f_uid;
    while (1)
    {
        gid = next_hap_label_id(&b, g_p, bub, u, 1, &c_bid, &c_uid);
        ///fprintf(stderr, "-gid: %u\n", gid);
        if(gid == (uint32_t)-1) break;
        ///fprintf(stderr, "vis[gid]: %u\n", i_buf->vis.a[gid]);
        if(i_buf->vis.a[gid] == 1) continue;
        lock[gid] = 1;
        i_buf->vis.a[gid] = 1;
        // if(f_bid == 689 && f_uid == 1)
        // {
        //     fprintf(stderr, "sssssssssssssssssssss\n");
        //     print_phase_group(g_p, bub, "Small-1");
        //     fprintf(stderr, "sbsbsbsb-gid: %u\n", gid);
        // }
        flip_unitig(g_p, link, bub, gid);
        ///flip_unitig_debug(g_p, link, bub, gid);
        // if(f_bid == 689 && f_uid == 1)
        // {
        //     fprintf(stderr, "sasasasa-gid: %u\n", gid);
        //     print_phase_group(g_p, bub, "Small-2");
        //     fprintf(stderr, "eeeeeeeeeeeeeeeeeeeeee\n");
        // }
        if(c_bid == l_bid && c_uid == l_uid) break;
    }

    return 1;
}

double get_total_weight(H_partition* h, G_partition* g_p)
{
    uint32_t i, k, uID;
    hc_links* link = h->link;
    int o_d = 0, o_f = 0;
    double w, t_w;
    for (i = 0, t_w = 0; i < h->n; i++)
    {
        if(g_p->index[i] == (uint32_t)-1) continue;
        o_f = g_p->a[g_p->index[i]>>1].status[g_p->index[i]&1];
        for (k = 0; k < link->a.a[i].e.n; k++)
        {
            if(link->a.a[i].e.a[k].del) continue;
            uID = link->a.a[i].e.a[k].uID;
            w = link->a.a[i].e.a[k].weight;
            if(g_p->index[uID] == (uint32_t)-1) continue;
            o_d = g_p->a[g_p->index[uID]>>1].status[g_p->index[uID]&1];
            t_w += (o_f*o_d*w);
        }
    }

    return t_w;
}

void hap_label_fliping(H_partition* hap, G_partition* g_p, bubble_type* bub, hc_links* link, uint32_t chain_id)
{
    long long c_bid, c_uid, l_bid, l_uid;
    uint32_t gid;
    ma_utg_t *u = &(bub->b_ug->u.a[chain_id]);
    hap->b.bid = hap->b.uid = 0;
    while (1)
    {
        gid = next_hap_label_id(&(hap->b), g_p, bub, u, 1, &c_bid, &c_uid);
        if(gid == (uint32_t)-1) break;
        identify_best_interval(&(hap->b), hap->lock, g_p, bub, u, link, c_bid, c_uid, &l_bid, &l_uid);
        if(l_bid == -1 || l_uid == -1) continue;
        
        // fprintf(stderr, "\nbefore weight: %f\n", get_total_weight(hap, g_p));
        // identify_best_interval_debug(&(hap->b), g_p, bub, u, link, c_bid, c_uid, &l_bid, &l_uid);
        if(flip_block(&(hap->b), g_p, bub, u, link, hap->lock, c_bid, c_uid, l_bid, l_uid))
        {
            ///fprintf(stderr, "after weight +: %f\n", get_total_weight(hap, g_p));
            hap->b.bid = l_bid;
            hap->b.uid = l_uid;
            gid = next_hap_label_id(&(hap->b), g_p, bub, u, 1, &c_bid, &c_uid);
            if(gid == (uint32_t)-1) break;
        }
        ///fprintf(stderr, "after weight -: %f\n", get_total_weight(hap, g_p));
        // exit(0);
    }
}

typedef struct{
    long long min_chain_id;
    long long min_f_bid;
    long long min_f_uid;
    long long min_l_bid; 
    long long min_l_uid;
    long long min_idx;
    double min_w;
}block_res_type;

typedef struct{
    block_phase_type* x;
    uint32_t n_thread;
    bubble_type* bub;
    uint64_t* chain_idx;///index of each chain
    uint64_t chain_idx_n;
    uint64_t chain_ele_occ;
    block_res_type* res;
    H_partition* h;
    G_partition* g_p;
}mul_block_phase_type;

uint32_t shift_block_phase_type(ma_utg_t *u, G_partition* g_p, bubble_type* bub, 
block_phase_type* b, uint32_t offset)
{
    long long c_bid, c_uid;
    uint32_t gid, occ = 0;
    b->bid = b->uid = 0;
    while (1)
    {
        if(occ == offset) break;
        gid = next_hap_label_id(b, g_p, bub, u, 1, &c_bid, &c_uid);
        if(gid == (uint32_t)-1) break;
        occ++;
    }
    return occ;
}

void get_block_phase_type(uint64_t* chain_idx, G_partition* g_p, bubble_type* bub, uint32_t id, block_phase_type* i_b)
{
    uint64_t i;
    ma_utg_t *u = NULL;
    for (i = 0; i < bub->chain_weight.n; i++)
    {
        if(id >= chain_idx[i] && id < chain_idx[i+1]) break;
    }

    u = &(bub->b_ug->u.a[bub->chain_weight.a[i].id]);
    shift_block_phase_type(u, g_p, bub, i_b, id - chain_idx[i]);
    i_b->chainID = bub->chain_weight.a[i].id;
}

void init_mul_block_phase_type(mul_block_phase_type* x, G_partition* g_p, bubble_type* bub, uint32_t n_thread, H_partition* hap)
{
    ma_utg_t *u = NULL;
    uint32_t i, n;
    block_phase_type b;
    x->bub = bub;
    x->n_thread = n_thread;
    CALLOC(x->res, x->n_thread);
    CALLOC(x->x, x->n_thread);
    for (i = 0; i < x->n_thread; i++)
    {
        kv_init(x->x[i].vis);
        kv_malloc(x->x[i].vis, hap->n); 
        x->x[i].vis.n = hap->n;
    }

    x->chain_idx_n = 0;
    MALLOC(x->chain_idx, bub->chain_weight.n+1);
    for (i = n = 0; i < bub->chain_weight.n; i++)
    {
        x->chain_idx[i] = n;
        if(bub->chain_weight.a[i].del) continue;
        u = &(bub->b_ug->u.a[bub->chain_weight.a[i].id]);
        n += shift_block_phase_type(u, g_p, bub, &b, (uint32_t)-1);
        x->chain_idx_n++;
    }
    x->chain_idx[i] = n;///index of chain
    x->chain_ele_occ = n;
}

void destory_mul_block_phase_type(mul_block_phase_type* x)
{
    uint32_t i;
    free(x->res);
    free(x->chain_idx);
    for (i = 0; i < x->n_thread; i++)
    {
        kv_destroy(x->x[i].vis);
    }
}

void select_max_block_by_utg_multi_thread(H_partition* h, G_partition* g_p, bubble_type* bub, 
hc_links* link, block_phase_type* i_b, uint64_t* chain_idx, uint32_t id, block_res_type* res)
{
    long long c_bid, c_uid, l_bid, l_uid;
    uint32_t gid;
    get_block_phase_type(chain_idx, g_p, bub, id, i_b);

    ma_utg_t *u = &(bub->b_ug->u.a[i_b->chainID]);
    
    gid = next_hap_label_id(i_b, g_p, bub, u, 1, &c_bid, &c_uid);

    if(gid == (uint32_t)-1) return;
    if(identify_best_interval(i_b, h->lock, g_p, bub, u, link, c_bid, c_uid, &l_bid, &l_uid)) return;
    if(l_bid == -1 || l_uid == -1) return;

    if((res->min_w > i_b->weight) || (res->min_w == i_b->weight && id < res->min_idx))
    {
        res->min_w = i_b->weight;
        res->min_f_bid = c_bid;
        res->min_f_uid = c_uid;
        res->min_l_bid = l_bid;
        res->min_l_uid = l_uid;
        res->min_chain_id = i_b->chainID;
        res->min_idx = id;
    }
}

static void worker_for_max_block(void *data, long i, int tid) // callback for kt_for()
{
    mul_block_phase_type* x = (mul_block_phase_type*)data;
    select_max_block_by_utg_multi_thread(x->h, x->g_p, x->bub, x->h->link, 
    &(x->x[tid]), x->chain_idx, i, &(x->res[tid]));
}

void select_max_block_by_utg_multi_thread_by_chain(H_partition* h, G_partition* g_p, bubble_type* bub, 
hc_links* link, block_phase_type* i_b, uint32_t id, block_res_type* res)
{
    long long c_bid, c_uid, l_bid, l_uid;
    uint32_t gid;
    ma_utg_t *u = &(bub->b_ug->u.a[id]);
    i_b->bid = i_b->uid = 0; i_b->chainID = id;
    while (1)
    {
        gid = next_hap_label_id(i_b, g_p, bub, u, 1, &c_bid, &c_uid);
        if(gid == (uint32_t)-1) break;

        if(identify_best_interval(i_b, h->lock, g_p, bub, u, link, c_bid, c_uid, &l_bid, &l_uid))
        {
            break;
        }
        if(l_bid == -1 || l_uid == -1) continue;
        if(res->min_w > i_b->weight)
        {
            res->min_w = i_b->weight;
            res->min_f_bid = c_bid;
            res->min_f_uid = c_uid;
            res->min_l_bid = l_bid;
            res->min_l_uid = l_uid;
            res->min_chain_id = i_b->chainID;
            res->min_idx = id;
        }
    }
}


int get_max_block_multi_thread(H_partition* h, G_partition* g_p, bubble_type* bub, mul_block_phase_type* x,
long long* min_u, long long* min_f_bid, long long* min_f_uid, long long* min_l_bid, long long* min_l_uid, 
double* min_w)
{
    uint32_t i;
    (*min_w) = 1;
    (*min_u) = (*min_f_bid) = (*min_f_uid) = (*min_l_bid) = (*min_l_uid) = -1;
    for (i = 0; i < x->n_thread; i++)
    {
        x->res[i].min_chain_id = x->res[i].min_f_bid = x->res[i].min_f_uid = -1;
        x->res[i].min_l_bid = x->res[i].min_l_uid = x->res[i].min_idx = -1; 
        x->res[i].min_w = 1;
    }
    x->g_p = g_p;
    x->h = h;
    kt_for(x->n_thread, worker_for_max_block, x, x->chain_ele_occ);
    ///kt_for(x->n_thread, worker_for_max_block_by_chain, x, x->chain_idx_n);

    long long min_idx = -1;
    for (i = 0; i < x->n_thread; i++)
    {
        if(x->res[i].min_chain_id == -1) continue;
        if(x->res[i].min_f_bid == -1 || x->res[i].min_f_uid == -1) continue;
        if(x->res[i].min_l_bid == -1 || x->res[i].min_l_uid == -1) continue;
        if(((*min_w) > x->res[i].min_w) || ((*min_w) == x->res[i].min_w && x->res[i].min_idx < min_idx))
        {
            (*min_w) = x->res[i].min_w;
            (*min_u) = x->res[i].min_chain_id;
            (*min_f_bid) = x->res[i].min_f_bid;
            (*min_f_uid) = x->res[i].min_f_uid;
            (*min_l_bid) = x->res[i].min_l_bid;
            (*min_l_uid) = x->res[i].min_l_uid;
            min_idx = x->res[i].min_idx;
        }
    }

    if((*min_u) != -1 && (*min_f_bid) != -1 && (*min_f_uid) != -1 && (*min_l_bid) != -1 && (*min_l_uid) != -1)
    {
        return 1;
    }

    return 0;
}

void select_max_block_by_utg(H_partition* hap, G_partition* g_p, bubble_type* bub, hc_links* link, uint32_t chain_id,
long long* min_f_bid, long long* min_f_uid, long long* min_l_bid, long long* min_l_uid, double* min_w)
{
    long long c_bid, c_uid, l_bid, l_uid;
    uint32_t gid;
    ma_utg_t *u = &(bub->b_ug->u.a[chain_id]);
    (*min_w) = 1;
    hap->b.bid = hap->b.uid = 0;
    (*min_f_bid) = (*min_f_uid) = (*min_l_bid) = (*min_l_uid) = -1;
    while (1)
    {
        gid = next_hap_label_id(&(hap->b), g_p, bub, u, 1, &c_bid, &c_uid);
        if(gid == (uint32_t)-1) break;

        if(identify_best_interval(&(hap->b), hap->lock, g_p, bub, u, link, c_bid, c_uid, &l_bid, &l_uid))
        {
            break;
        }
        if(l_bid == -1 || l_uid == -1) continue;
        if((*min_w) > hap->b.weight)
        {
            (*min_w) = hap->b.weight;
            (*min_f_bid) = c_bid;
            (*min_f_uid) = c_uid;
            (*min_l_bid) = l_bid;
            (*min_l_uid) = l_uid;
        }
    }
}


int get_max_block(H_partition* h, G_partition* g_p, bubble_type* bub, long long* min_u,
long long* min_f_bid, long long* min_f_uid, long long* min_l_bid, long long* min_l_uid, 
double* min_w)
{
    uint32_t i;
    long long f_bid, f_uid, l_bid, l_uid;
    double w;
    (*min_w) = 1;
    (*min_u) = (*min_f_bid) = (*min_f_uid) = (*min_l_bid) = (*min_l_uid) = -1;
    for (i = 0; i < bub->chain_weight.n; i++)
    {
        if(bub->chain_weight.a[i].del) continue;
        select_max_block_by_utg(h, g_p, bub, h->link, bub->chain_weight.a[i].id,
        &f_bid, &f_uid, &l_bid, &l_uid, &w);
        if(f_bid == -1 || f_uid == -1 || l_bid == -1 || l_uid == -1) continue;
        if((*min_w) > w)
        {
            (*min_w) = w;
            (*min_u) = bub->chain_weight.a[i].id;
            (*min_f_bid) = f_bid;
            (*min_f_uid) = f_uid;
            (*min_l_bid) = l_bid;
            (*min_l_uid) = l_uid;
        }
    }

    if((*min_u) != -1 && (*min_f_bid) != -1 && (*min_f_uid) != -1 && (*min_l_bid) != -1 && (*min_l_uid) != -1)
    {
        return 1;
    }

    return 0;
}

void phasing_improvement_by_block(H_partition* h, G_partition* g_p, bubble_type* bub, mul_block_phase_type* x)
{
    long long min_u, min_f_bid, min_f_uid, min_l_bid, min_l_uid; 
    double min_w;

    memset(h->lock, 0, sizeof(uint8_t)*h->n);
    while(get_max_block_multi_thread(h, g_p, bub, x, &min_u, &min_f_bid, &min_f_uid, &min_l_bid, &min_l_uid, &min_w))
    ///while(get_max_block(h, g_p, bub, &min_u, &min_f_bid, &min_f_uid, &min_l_bid, &min_l_uid, &min_w))
    {
        ///fprintf(stderr, "\nmin_w: %f, min_u: %lld, min_f_bid: %lld, min_f_uid: %lld, min_l_bid: %lld, min_l_uid: %lld\n", min_w, min_u, min_f_bid, min_f_uid, min_l_bid, min_l_uid);
        ///fprintf(stderr, "before weight: %f\n", get_total_weight(h, g_p));
        flip_block(&(h->b), g_p, bub, &(bub->b_ug->u.a[min_u]), h->link, h->lock, min_f_bid, 
        min_f_uid, min_l_bid, min_l_uid);
        ///fprintf(stderr, "after weight: %f\n", get_total_weight(h, g_p));
    }
}

void flip_by_chain(H_partition* h, G_partition* g_p, bubble_type* bub)
{
    uint32_t i;
    memset(h->lock, 0, sizeof(uint8_t)*h->n);
    for (i = 0; i < bub->chain_weight.n; i++)
    {
        if(bub->chain_weight.a[i].del) continue;
        merge_phase_group_by_chain(h, g_p, bub, bub->chain_weight.a[i].id);
    }

    double pre_w = get_total_weight(h, g_p), current_w;
    uint32_t round = 0;
    while (1)
    {
        memset(h->lock, 0, sizeof(uint8_t)*h->n);
        while (1)
        {
            i = get_max_unitig(h, g_p, h->link, bub);
            if(i == (uint32_t)-1) break;
            h->lock[i] = 1;
            flip_unitig(g_p, h->link, bub, i);
        }
        current_w = get_total_weight(h, g_p);
        ///fprintf(stderr, "[M::%s::round %u, pre_w: %f, current_w: %f]\n", __func__, round, pre_w, current_w);
        if(ceil(current_w) <= ceil(pre_w)) break;
        round++;
        pre_w = current_w;
    }


    ///print_phase_group(g_p, bub, "Large-pre");
    
    ///fprintf(stderr, "[M::%s::round %u, before block flipping: %f]\n", __func__, round, get_total_weight(h, g_p));

    mul_block_phase_type b_x;
    init_mul_block_phase_type(&b_x, g_p, bub, asm_opt.thread_num, h);

    pre_w = get_total_weight(h, g_p);
    while (1)
    {
        phasing_improvement_by_block(h, g_p, bub, &b_x);
        current_w = get_total_weight(h, g_p);
        ///fprintf(stderr, "[M::%s::round %u, after block flipping: %f]\n", __func__, round, get_total_weight(h, g_p));
        ///debug_flip(g_p, h->link, bub, 0);
        if(ceil(current_w) <= ceil(pre_w)) break;
        round++;
        pre_w = current_w;
    }
    destory_mul_block_phase_type(&b_x);

    for (i = 0; i < g_p->n; i++)
    {
        update_partition_flag(h, g_p, h->link, i);
    }
}

void flip_by_node(H_partition* h, G_partition* g_p, bubble_type* bub)
{
    uint32_t i;
    memset(h->lock, 0, sizeof(uint8_t)*h->n);
    

    double pre_w = get_total_weight(h, g_p), current_w;
    uint32_t round = 0;
    while (1)
    {
        memset(h->lock, 0, sizeof(uint8_t)*h->n);
        while (1)
        {
            i = get_max_unitig(h, g_p, h->link, bub);
            if(i == (uint32_t)-1) break;
            h->lock[i] = 1;
            flip_unitig(g_p, h->link, bub, i);
        }
        current_w = get_total_weight(h, g_p);
        ///fprintf(stderr, "[M::%s::round %u, pre_w: %f, current_w: %f]\n", __func__, round, pre_w, current_w);
        if(ceil(current_w) <= ceil(pre_w)) break;
        round++;
        pre_w = current_w;
    }
    
    ///fprintf(stderr, "[M::%s::round %u, before block flipping: %f]\n", __func__, round, get_total_weight(h, g_p));
    
    for (i = 0; i < g_p->n; i++)
    {
        update_partition_flag(h, g_p, h->link, i);
    }
}


void link_phase_group(H_partition* hap, bubble_type* bub)
{
    double index_time = yak_realtime();
    uint32_t i, k, n = (hap->label>>hap->label_shift)+1, *h0, h0_n, *h1, h1_n;;
    init_G_partition(&(hap->group_g_p), hap->n);
    partition_warp *res = NULL;
    for (i = 0; i < n; i++)///how many haplotype group
    {
        kv_pushp(partition_warp, hap->group_g_p, &res);
        kv_init(res->a);
        res->full_bub = 0;
        res->h[0] = res->h[1] = 0;
        res->status[0] = 1; res->status[1] = -1;
        res->weight[0] = res->weight[1] = res->weight_convex = 0;
        ///all unitigs
        for (k = 0; k < hap->n; k++)
        {
            if(get_phase_group(hap, k) == i && get_phase_status(hap, k) == 1)
            {
                kv_push(uint32_t, res->a, k);
                res->h[0]++;
            }
        }

        for (k = 0; k < hap->n; k++)
        {
            if(get_phase_group(hap, k) == i && get_phase_status(hap, k) == -1)
            {
                kv_push(uint32_t, res->a, k);
                res->h[1]++;
            }
        }

        for (k = 0; k < res->h[0]; k++)
        {
            ///if(hap->group_g_p.index[res->a.a[k]] != (uint32_t)-1) fprintf(stderr, "ERROR---00\n");
            hap->group_g_p.index[res->a.a[k]] = hap->group_g_p.n-1;
            hap->group_g_p.index[res->a.a[k]] = hap->group_g_p.index[res->a.a[k]] << 1;
        } 

        for (; k < res->a.n; k++)
        {
            ///if(hap->group_g_p.index[res->a.a[k]] != (uint32_t)-1) fprintf(stderr, "ERROR---11\n");
            hap->group_g_p.index[res->a.a[k]] = hap->group_g_p.n-1;
            hap->group_g_p.index[res->a.a[k]] = (hap->group_g_p.index[res->a.a[k]] << 1) + 1;
        }

        get_phased_block(&(hap->group_g_p), NULL, i, NULL, NULL, &h0, &h0_n, &h1, &h1_n, NULL, NULL);
        if(h0_n >0) res->weight[0] = get_cluster_weight(hap, hap->link, h0, h0_n); 
        if(h1_n >0) res->weight[1] = get_cluster_weight(hap, hap->link, h1, h1_n); 
        res->weight_convex = get_cluster_inner_weight(hap, hap->link, h0, h0_n, h1, h1_n);
    }

    fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
    // for (i = 0; i < n; i++)
    // {
    //     double w0 = 0, w1 = 0;
    //     res = &(hap->group_g_p.a[i]);
    //     fprintf(stderr, "%u-th group: # %d = %u, # %d = %u\n", i,
    //     res->status[0], res->h[0], res->status[1], res->h[1]);
    //     get_phased_block(&(hap->group_g_p), NULL, i, NULL, NULL, &h0, &h0_n, &h1, &h1_n, NULL, NULL);
    //     if(h0_n >0) w0 = get_cluster_weight_debug(&(hap->group_g_p), hap->link, h0, h0_n); 
    //     if(h1_n >0) w1 = get_cluster_weight_debug(&(hap->group_g_p), hap->link, h1, h1_n); 
    //     if(w0 != res->weight[0]) fprintf(stderr, "i: %u, ERROR: w0: %f, weight[0]: %f\n", i, w0, res->weight[0]);
    //     if(w1 != res->weight[1]) fprintf(stderr, "i: %u, ERROR: w1: %f, weight[1]: %f\n", i, w1, res->weight[1]);


    //     for (k = 0; k < res->h[0]; k++)
    //     {
    //         fprintf(stderr, "%d: utg%.6ul\n", res->status[0], int(res->a.a[k]+1));
    //     }

    //     for (; k < res->a.n; k++)
    //     {
    //         fprintf(stderr, "%d: utg%.6ul\n", res->status[1], int(res->a.a[k]+1));
    //     }
    // }

    /*******************************for debug************************************/
    // for (i = 0; i < hap->n; i++)
    // {
    //     if(hap->link->a.a[i].e.n == 0) continue;
    //     if(get_phase_status(hap, i) == -2)
    //     {
    //         fprintf(stderr, "ERROR+++: i: %u, group: %u, bub->index: %u\n", i, get_phase_group(hap, i), bub->index[i]);
    //         for (k = 0; k < hap->link->a.a[i].e.n; k++)
    //         {
    //             fprintf(stderr, "k: %u, uID: %u, weight: %f, del: %u\n", k, hap->link->a.a[i].e.a[k].uID, 
    //                         hap->link->a.a[i].e.a[k].weight, hap->link->a.a[i].e.a[k].del);   
    //         }
    //     } 
    // }
    /*******************************for debug************************************/
    

    flip_by_chain(hap, &(hap->group_g_p), bub);
    ///print_phase_group(&(hap->group_g_p), bub, "Large");
}

void print_chain_phasing(H_partition* hap, ma_ug_t *ug, bubble_type* bub, uint32_t chain_id)
{
    uint32_t i, k;
    uint32_t beg, sink, *a, n;
    uint64_t bid, uid;
    ma_utg_t *u = &(bub->b_ug->u.a[chain_id]);
    fprintf(stderr, "\n**********chain_id: %u**********\n", chain_id);
    for (i = 0; i < u->n; i++)
    {
        bid = u->a[i]>>33;
        fprintf(stderr, "(%u) chain_id: %u, u->n: %u, bid: %u\n", i, chain_id, (uint32_t)u->n, i);
        get_bubbles(bub, bid, &beg, &sink, &a, &n, NULL);
        for (k = 0; k < n; k++)
        {
            uid = a[k]>>1;
            fprintf(stderr, "utg%.6ul, hap: %u, group: %u, stats: %d\n", (int)(uid+1), hap->hap[uid], 
                        get_phase_group(hap, uid), get_phase_status(hap, uid));
        }
    }
}

int graph_bipartiteness(uint32_t* b_a, uint32_t b_a_n, uint8_t *color, hc_links* link, kvec_t_u32_warp* stack)
{
    if(b_a_n == 0) return 0;
    uint32_t i, uID, cur, occ = 0, sucess = 0, c;
    for (i = 0; i < b_a_n; i++) color[b_a[i]>>1] = 8;
    stack->a.n = 0;
    kv_push(uint32_t, stack->a, b_a[0]>>1);
    while (stack->a.n > 0)
    {
        stack->a.n--;
        cur = stack->a.a[stack->a.n];
        if((color[cur] & 1) == 0) occ++;
        color[cur] |= 1;
        for (i = 0; i < link->a.a[cur].f.n; i++)
        {
            if(link->a.a[cur].f.a[i].del) continue;
            if(link->a.a[cur].f.a[i].dis != RC_0) continue;
            uID = link->a.a[cur].f.a[i].uID;
            if((color[uID] & 8) == 0) continue;
            if((color[uID] & 1) == 1) continue;
            kv_push(uint32_t, stack->a, uID);
        }
    }
    if(occ != b_a_n) goto Failed;

    sucess = 1;
    for (i = 0; i < b_a_n; i++) color[b_a[i]>>1] = 8;
    stack->a.n = 0;
    kv_push(uint32_t, stack->a, b_a[0]>>1);
    color[b_a[0]>>1] |= 2;///colored
    while (stack->a.n > 0)
    {
        stack->a.n--;
        cur = stack->a.a[stack->a.n];
        color[cur] |= 1;
        c = color[cur] & 4; ///get color
        for (i = 0; i < link->a.a[cur].f.n; i++)
        {
            if(link->a.a[cur].f.a[i].del) continue;
            if(link->a.a[cur].f.a[i].dis != RC_0) continue;
            uID = link->a.a[cur].f.a[i].uID;
            if((color[uID] & 8) == 0) continue;
            if((color[uID] & 2) && ((color[uID] & 4) == c)) break; ///conflict
            if((color[uID] & 1) == 1) continue;
            kv_push(uint32_t, stack->a, uID);
            color[uID] |= 2; color[uID] |= (c^4);
        }

        if(i != link->a.a[cur].f.n)
        {
            sucess = -1;
            break;
        }
    }

    Failed: 
    if(sucess != 1)
    {
        for (i = 0; i < b_a_n; i++) color[b_a[i]>>1] = 0;
    }
    
    return sucess;
}

void assign_per_unitig_G_partition(G_partition* g_p, uint64_t hap_n, hc_links* link, bubble_type* bub,
uint32_t bubble_first)
{
    double index_time = yak_realtime();
    reset_G_partition(g_p, hap_n);

    partition_warp* res = NULL;
    hc_edge *a = NULL;
    uint32_t i, a_n, v, u, uv = (uint32_t)-1, k, k_n, k_nv, k_nu, beg, sink, *b_a = NULL, b_a_n;
    
    if(bubble_first)
    {
        int c;
        kvec_t_u32_warp stack; kv_init(stack.a);
        uint8_t *color = NULL; CALLOC(color, hap_n);
        uint32_t n_bub = bub->f_bub + bub->b_bub;
        for (i = 0; i < n_bub; i++)
        {
            get_bubbles(bub, i, &beg, &sink, &b_a, &b_a_n, NULL);
            if(b_a_n == 2 && i < bub->f_bub)
            {
                continue;
            }
            ///full bubble do not overlap with any others
            ///broken bubbles might be, but should do nothing
            c = graph_bipartiteness(b_a, b_a_n, color, link, &stack);

            if(c == 0)
            {
                fprintf(stderr, "too good: s-utg%.6ul && e-utg%.6ul && %s\n",(beg>>1)+1, (sink>>1)+1, b_a_n != 4? "abnormal" : "normal");
            }
            if(c == -1)
            {
                fprintf(stderr, "too bad: s-utg%.6ul && e-utg%.6ul\n",(beg>>1)+1, (sink>>1)+1);
            }
            if(c == 1)
            {
                fprintf(stderr, "\nprefect=%u: s-utg%.6ul && e-utg%.6ul\n", b_a_n, (beg>>1)+1, (sink>>1)+1);
                
                
                for (k = 0; k < b_a_n; k++)
                {
                    if((color[b_a[k]>>1] & 2) == 0) fprintf(stderr, "ERROR\n");
                    if((color[b_a[k]>>1] & 4) == 0) fprintf(stderr, "0: utg%.6ul\n", (b_a[k]>>1)+1);
                }

                for (k = 0; k < b_a_n; k++)
                {
                    if((color[b_a[k]>>1] & 2) == 0) fprintf(stderr, "ERROR\n");
                    if((color[b_a[k]>>1] & 4) != 0) fprintf(stderr, "1: utg%.6ul\n", (b_a[k]>>1)+1);
                }

                for (k = 0; k < b_a_n; k++) color[b_a[k]>>1] = 0;
            }
        }
        free(color);
        kv_destroy(stack.a);
    }

    for (i = 0; i < hap_n; i++)
    {
        v = i;
        a = link->a.a[v].f.a;
        a_n = link->a.a[v].f.n;
        for (k = k_n = 0; k < a_n; k++)
        {
            if(a[k].del) continue;
            if(a[k].dis != RC_0) break;
            u = a[k].uID;
            k_n++;
        }
        if(k_n != 1)
        {   
            u = (uint32_t)-1;
            goto push_uv;
        } 
        
        a = link->a.a[u].f.a;
        a_n = link->a.a[u].f.n;
        for (k = k_n = 0; k < a_n; k++)
        {
            if(a[k].del) continue;
            if(a[k].dis != RC_0) break;
            uv = a[k].uID;
            k_n++;
        }
        if(k_n != 1 || uv != v)
        {
            u = (uint32_t)-1;
            goto push_uv;
        } 

        push_uv:
        k_nv = 0;k_nu = 0;

        // not such easy. need to deal with here very carefully
        // if(g_p->index[v] != (uint32_t)-1) continue; 
        // if(u != (uint32_t)-1 && g_p->index[u] != (uint32_t)-1) u = (uint32_t)-1;

        a = link->a.a[v].e.a;
        a_n = link->a.a[v].e.n;
        for (k = 0; k < a_n; k++)
        {
            if(a[k].del) continue;
            k_nv++;
        }
                
        if(u != (uint32_t)-1)
        {
            a = link->a.a[u].e.a;
            a_n = link->a.a[u].e.n;
            for (k = 0; k < a_n; k++)
            {
                if(a[k].del) continue;
                k_nu++;
            }
        }
        if(k_nv == 0) continue;
        if(k_nv > 0 && k_nu > 0 && v > u) continue;

        kv_pushp(partition_warp, *g_p, &res);
        kv_init(res->a);
        res->full_bub = 0;
        res->h[0] = 1; res->h[1] = 0;
        kv_push(uint32_t, res->a, v);
        if(u != (uint32_t)-1)
        {
            res->h[1] = 1;
            kv_push(uint32_t, res->a, u);
        }

        for (k = 0; k < res->h[0]; k++)
        {
            g_p->index[res->a.a[k]] = g_p->n-1;
            g_p->index[res->a.a[k]] = g_p->index[res->a.a[k]] << 1;
        } 

        for (; k < res->a.n; k++)
        {
            g_p->index[res->a.a[k]] = g_p->n-1;
            g_p->index[res->a.a[k]] = (g_p->index[res->a.a[k]] << 1) + 1;
        } 
    }
    fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
}

typedef struct {
    double weight;
    uint64_t p_id, beg_idx, end_idx;
    uint8_t used;
}bub_sort_type;

typedef struct {
    bub_sort_type* a;
    size_t n, m;
}bub_sort_vec;

double get_specific_weight_by_chain(uint64_t* ids, uint64_t beg_idx, uint64_t end_idx, uint64_t p_id,
hc_links* link, uint8_t* vis, uint8_t flag)
{
    uint64_t x, k;
    uint32_t uid;
    double w;
    for (x = beg_idx, w = 0; x <= end_idx; x++)
    {
        uid = (uint32_t)((uint32_t)ids[x])>>1;
        for (k = 0; k < link->a.a[uid].e.n; k++)
        {
            if(link->a.a[uid].e.a[k].del) continue;
            if(vis[link->a.a[uid].e.a[k].uID] != flag) continue;
            w += link->a.a[uid].e.a[k].weight;
        }
    }

    return w;
}

int cmp_bubble_ele_by_chain(const void * a, const void * b)
{
    if((*(bub_sort_type*)a).weight == (*(bub_sort_type*)b).weight)
    {
        return (*(bub_sort_type*)a).weight > (*(bub_sort_type*)b).weight? -1 : 1;
    }

    return 0;
}

uint32_t get_max_hap_g(bub_sort_vec* w_stack, uint32_t* require_iso)
{
    uint32_t k, max_idx = (uint32_t)-1;
    double max_w;
    for (k = require_iso? (*require_iso)+1 : 0, max_idx = (uint32_t)-1; k < w_stack->n; k++)
    {
        if(w_stack->a[k].used) continue;
        if(!require_iso)
        {
            if(w_stack->a[k].p_id == (uint32_t)-1) continue;
            if((max_idx == (uint32_t)-1) || (max_idx != (uint32_t)-1 && max_w < w_stack->a[k].weight))
            {
                max_w = w_stack->a[k].weight;
                max_idx = k;
            }
        }
        else
        {
            if(w_stack->a[k].p_id != (uint32_t)-1) continue;
            return k;
        }
    }

    if(require_iso && (*require_iso) != 0)
    {
        for (k = 0; k < w_stack->n; k++)
        {
            if(w_stack->a[k].used) continue;
            if(w_stack->a[k].p_id != (uint32_t)-1) continue;
            return k;
        }
    }

    return max_idx;
}

void update_bub_sort_vec(uint64_t* ids, bub_sort_vec* w_stack, uint32_t max_idx, hc_links* link,
uint32_t* set_hap)
{
    w_stack->a[max_idx].used = 1;
    uint64_t i, k;
    uint32_t uid, pid_idx;
    for (i = w_stack->a[max_idx].beg_idx; i <= w_stack->a[max_idx].end_idx; i++)
    {
        uid = (uint32_t)((uint32_t)ids[i])>>1;
        for (k = 0; k < link->a.a[uid].e.n; k++)
        {
            if(link->a.a[uid].e.a[k].del) continue;
            pid_idx = set_hap[link->a.a[uid].e.a[k].uID];
            if(pid_idx == (uint32_t)-1) continue;
            if(w_stack->a[pid_idx].used) continue;
            w_stack->a[pid_idx].weight += link->a.a[uid].e.a[k].weight;
        }
    }
}

void sort_bubble_ele_by_chain(G_partition* g_p, hc_links* link, bubble_type* bub, kvec_t_u64_warp* stack, 
bub_sort_vec* w_stack, uint8_t* vis, uint32_t* set_hap, uint32_t n_utg, uint32_t chain_id)
{
    uint32_t max_idx, i, k, j, m, beg, sink, *a, n, flag_cur = 3, flag_right = 2, flag_left = 1, flag_unset = 0;
    uint64_t bid, uid, pid, pre_pid;
    ma_utg_t *u = &(bub->b_ug->u.a[chain_id]);
    bub_sort_type *p = NULL;
    memset(vis, flag_unset, n_utg);
    for (i = 0; i < u->n; i++)
    {
        bid = u->a[i]>>33;
        get_bubbles(bub, bid, &beg, &sink, &a, &n, NULL);
        for (k = 0; k < n; k++)
        {
            uid = a[k]>>1;
            vis[uid] = flag_right;
        }
    }

    for (i = 0; i < u->n; i++)
    {
        stack->a.n = 0; w_stack->n = 0;
        bid = u->a[i]>>33;
        get_bubbles(bub, bid, &beg, &sink, &a, &n, NULL);
        for (k = 0; k < n; k++)
        {
            uid = a[k]>>1;
            pid = g_p->index[uid];
            if(pid != (uint32_t)-1) pid >>= 1;
            kv_push(uint64_t, stack->a, (pid<<32)|a[k]);
            vis[uid] = flag_cur;
        }
        radix_sort_hc64(stack->a.a, stack->a.a + stack->a.n);///sort is to dedup pid

        for (k = 0, pre_pid = (uint64_t)-1; k < stack->a.n; k++)
        {
            if((stack->a.a[k]>>32) == pre_pid) continue;
            if(w_stack->n > 0) w_stack->a[w_stack->n-1].end_idx = k - 1;

            pre_pid = stack->a.a[k]>>32;
            kv_pushp(bub_sort_type, *w_stack, &p);
            p->weight = 0;
            p->p_id = pre_pid;
            p->beg_idx = k;
            p->end_idx = (uint64_t)-1;
            p->used = 0;
        }
        if(w_stack->n > 0) w_stack->a[w_stack->n-1].end_idx = k - 1;

        ///get each hap id
        for (k = 0; k < w_stack->n; k++)
        {
            w_stack->a[k].weight += get_specific_weight_by_chain(stack->a.a, w_stack->a[k].beg_idx, w_stack->a[k].end_idx, w_stack->a[k].p_id, 
            link, vis, flag_left);
            w_stack->a[k].weight -= get_specific_weight_by_chain(stack->a.a, w_stack->a[k].beg_idx, w_stack->a[k].end_idx, w_stack->a[k].p_id, 
            link, vis, flag_right);
            for (j = w_stack->a[k].beg_idx; j <= w_stack->a[k].end_idx; j++)
            {
                set_hap[((uint32_t)stack->a.a[j])>>1] = k;
            }
        }

        m = 0;
        while ((max_idx = get_max_hap_g(w_stack, NULL)) != (uint32_t)-1)
        {
            for (j = w_stack->a[max_idx].beg_idx; j <= w_stack->a[max_idx].end_idx; j++)
            {
                a[m] = (uint32_t)stack->a.a[j];
                m++;
            }
            update_bub_sort_vec(stack->a.a, w_stack, max_idx, link, set_hap);
        }

        while ((max_idx = get_max_hap_g(w_stack, &max_idx)) != (uint32_t)-1)
        {
            for (j = w_stack->a[max_idx].beg_idx; j <= w_stack->a[max_idx].end_idx; j++)
            {
                a[m] = (uint32_t)stack->a.a[j];
                m++;
            }
            update_bub_sort_vec(stack->a.a, w_stack, max_idx, link, set_hap);
        }
        

        /**
        qsort(w_stack->a, w_stack->n, sizeof(bub_sort_type), cmp_bubble_ele_by_chain);

        m = 0;

        for (k = 0; k < w_stack->n; k++)
        {
            if(w_stack->a[k].p_id == (uint32_t)-1) continue;
            for (j = w_stack->a[k].beg_idx; j <= w_stack->a[k].end_idx; j++)
            {
                a[m] = (uint32_t)stack->a.a[j];
                m++;
            }
        }

        for (k = 0; k < w_stack->n; k++)
        {
            if(w_stack->a[k].p_id != (uint32_t)-1) continue;
            for (j = w_stack->a[k].beg_idx; j <= w_stack->a[k].end_idx; j++)
            {
                a[m] = (uint32_t)stack->a.a[j];
                m++;
            }
        }
        **/

        for (k = 0; k < n; k++)
        {
            uid = a[k]>>1;
            vis[uid] = flag_left;
            set_hap[uid] = (uint32_t)-1;
        }
    }
}

void sort_bubble_ele(G_partition* g_p, hc_links* link, bubble_type* bub, uint32_t n_utg)
{
    double index_time = yak_realtime();
    kvec_t_u64_warp stack; kv_init(stack.a);
    bub_sort_vec w_stack; kv_init(w_stack);
    uint8_t* vis = NULL; MALLOC(vis, n_utg);
    uint32_t* set_hap = NULL; MALLOC(set_hap, n_utg); memset(set_hap, -1, sizeof(uint32_t)*n_utg);
    uint32_t i;
    
   for (i = 0; i < bub->chain_weight.n; i++)
   {
       if(bub->chain_weight.a[i].del) continue;
       sort_bubble_ele_by_chain(g_p, link, bub, &stack, &w_stack, vis, set_hap, n_utg, bub->chain_weight.a[i].id);
   }

   kv_destroy(stack.a); kv_destroy(w_stack); free(vis); free(set_hap);
   fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
}

uint32_t init_contig_partition(H_partition* hap, ha_ug_index* idx, bubble_type* bub, hc_links* link)
{
    ///hc_links* link = idx->link;
    ma_ug_t *ug = idx->ug;
    bub_p_t_warp b;
    memset(&b, 0, sizeof(bub_p_t_warp));
    CALLOC(b.a, ug->g->n_seq*2); 
    uint32_t i, nv = ug->g->n_seq * 2, max_i, max_hap_label;
    for (i = 0; i < nv; i++)
    {
        b.a[i].w[0] = b.a[i].w[1] = b.a[i].nh = 0;
        b.a[i].p =b.a[i].d = b.a[i].nc = b.a[i].uc = b.a[i].ac = b.a[i].r = b.a[i].s = 0;
    }
    uint8_t* hap_label_flag = NULL;
    CALLOC(hap_label_flag, ug->g->n_seq);
    

    hap->n = ug->u.n;
    MALLOC(hap->hap, hap->n);
    memset(hap->hap, 0, hap->n*sizeof(uint32_t));
    MALLOC(hap->lock, hap->n);
    memset(hap->lock, 0, hap->n);
    hap->m[0] = 1; hap->m[1] = 2; hap->m[2] = 4;
    hap->link = link;
    hap->label = 0;
    hap->label_add = 8;//1000, for hap group
    for(hap->label_shift=1; (uint64_t)(1<<hap->label_shift)<(uint64_t)hap->label_add; hap->label_shift++);
    kv_init(hap->label_buffer);
    kv_init(hap->b.vis); kv_malloc(hap->b.vis, hap->n); hap->b.vis.n = hap->n;

    ///sorted by weight
    for (i = 0; i < bub->chain_weight.n; i++)
    {
        if(bub->chain_weight.a[i].del) continue;
        phase_bubble_chain(hap, ug, &b, bub, hap_label_flag, bub->chain_weight.a[i].id);
        ///print_chain_phasing(hap, ug, bub, bub->chain_weight.a[i].id);
    }
    
    memset(hap_label_flag, 1, ug->g->n_seq);
    while (1)
    {
        max_i = get_unset_com(hap, bub, ug, hap_label_flag, &max_hap_label);
        if(max_i == (uint32_t)-1) break;
        phase_com(hap, ug, &b, bub, max_i, max_hap_label);
    }

    for (i = 0; i < hap->n; i++)
    {
        if((hap->hap[i]&hap->m[0])&&(hap->hap[i]&hap->m[1]))
        {
            hap->hap[i] >>= hap->label_shift;
            hap->hap[i] <<= hap->label_shift;
            hap->hap[i] |= hap->m[2];
            reset_ambiguous_label(hap, hap_label_flag, i);
        }
    }

    free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);

    init_G_partition(&(hap->g_p), hap->n);

    link_phase_group(hap, bub);

    assign_per_unitig_G_partition(&(hap->g_p), hap->n, link, bub, 0);///warp unitigs

    adjust_contig_partition(hap, link);

    update_bubble_chain(ug, bub, 0, 1);

    resolve_bubble_chain_tangle(ug, bub);

    clean_bubble_chain_by_HiC(ug, link, bub);

    append_boundary_chain(ug, link, bub);

    sort_bubble_ele(&(hap->g_p), link, bub, hap->n);

    free(hap_label_flag);
    return 1;
}

double get_path_phasing_weight(uint32_t query, uint32_t v0, uint32_t root, bub_p_t_warp *b, kv_u_trans_t *ta)
{
    if(v0 == root) return 0;
    uint32_t v, u;
    u_trans_t *p = NULL;
    double nw = 0;
    v = v0;
    do {
        u = b->a[v].p; // u->v
        get_u_trans_spec(ta, query>>1, v>>1, &p, NULL);
        if(p && (!p->del)) nw += p->nw;
        v = u;
    } while (v != root);

    return nw;
}

void get_related_phasing_weight(uint32_t x, kv_u_trans_t *ta, double* w0, double* w1, int8_t *s)
{
    (*w0) = (*w1) = 0;
    if(x >= ta->idx.n) return;
    uint32_t e_n, k;
    u_trans_t* e = u_trans_a(*ta, x);
    e_n = u_trans_n(*ta, x);
    for (k = 0; k < e_n; k++)
    {
        if(e[k].del) continue;
        if(s[e[k].tn] > 0) (*w0) += e[k].nw;
        else if(s[e[k].tn] < 0) (*w1) += e[k].nw;
    }
}

void set_phase_path(bub_p_t_warp *b, uint32_t root, kv_u_trans_t *ta, ps_t *s)
{
    int8_t f;
    uint32_t v, u;
    double z[2], cur_w[2];
    z[0] = z[1] = 0;

    v = b->S.a[0];
    do {
        u = b->a[v].p; // u->v
        if(v != b->S.a[0]) 
        {
            get_related_phasing_weight(v>>1, ta, &cur_w[0], &cur_w[1], s->s);
            z[0] += cur_w[0]; z[1] += cur_w[1];
        }
        v = u;
    } while (v != root);

    if(z[0] - z[1] < 0)
    {
        f = 1;
    }
    else if(z[0] - z[1] > 0)
    {
        f = -1;
    }
    else
    {
        s->xs = kr_splitmix64(s->xs);
        f = s->xs&1? 1 : -1;
    }

    v = b->S.a[0];
    do {
        u = b->a[v].p; // u->v
        if(v != b->S.a[0]) s->s[v>>1] = f;
        v = u;
    } while (v != root);
}

uint32_t bub_phase(ma_ug_t *ug, uint32_t beg, uint32_t end, bub_p_t_warp *b, kv_u_trans_t *ta, ps_t *s)
{
    asg_t *g = ug->g;
    if(g->seq[beg>>1].del) return 0; // already deleted
    if(get_real_length(g, beg, NULL)<2) return 0;
    uint32_t i, is_end, n_pending, to_replace, cur_nc, cur_uc, cur_ac, n_tips, tip_end, n_pop;
    double cur_nh, cur_w0, cur_w1, cur_rate, max_rate, cur_weight, min_weight;
    ///S saves nodes with all incoming edges visited
    b->S.n = b->T.n = b->b.n = b->e.n = 0;
    ///for each node, b->a saves all related information
    b->a[beg].d = b->a[beg].nc = b->a[beg].ac = b->a[beg].uc = 0; 
    b->a[beg].nh = b->a[beg].w[0] = b->a[beg].w[1] = 0;
    b->a[beg].p = (uint32_t)-1;
    ///b->S is the nodes with all incoming edges visited
    kv_push(uint32_t, b->S, beg);
    n_pop = n_tips = n_pending = 0;
    tip_end = (uint32_t)-1;

    do {
        ///v is a node that all incoming edges have been visited
        ///d is the distance from v0 to v
        uint32_t v = kv_pop(b->S);
        uint32_t d = b->a[v].d, nc = b->a[v].nc, uc = b->a[v].uc, ac = b->a[v].ac; 
        double nh = b->a[v].nh;///path weight
        double nw_0 = b->a[v].w[0], nw_1 = b->a[v].w[1];///weight to haplotype 1/2

        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        for (i = 0; i < nv; ++i) {
            uint32_t w = av[i].v, l = (uint32_t)av[i].ul; // v->w with length l, not overlap length
            bub_p_t *t = &b->a[w];
            is_end = 0;
            if((w>>1) == (end>>1)) is_end = 1;
            //got a circle
            if ((w>>1) == (beg>>1)) goto pop_reset;
            if (av[i].del) continue;
            ///push the edge
            kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);

            if (t->s == 0) 
            { // this vertex has never been visited
                kv_push(uint32_t, b->b, w); // save it for revert
                ///t->p is the parent node of 
                ///t->s = 1 means w has been visited
                ///d is len(v0->v), l is len(v->w), so t->d is len(v0->w)
                t->p = v, t->s = 1, t->d = d + l; 
                t->r = get_real_length(g, w^1, NULL);

                if(is_end == 0)
                {
                    t->nc = nc + ug->u.a[(w>>1)].n;
                    t->nh = nh + get_path_phasing_weight(w, v, beg, b, ta);
                    get_related_phasing_weight(w>>1, ta, &(t->w[0]), &(t->w[1]), s->s);
                    t->w[0] += nw_0; t->w[1] += nw_1; 
                    t->ac = ac + ((s->s[w>>1] == 0)? ug->u.a[(w>>1)].n : 0);
                    t->uc = uc + ((s->s[w>>1] != 0)? ug->u.a[(w>>1)].n : 0);
                }

                ++n_pending;
            }
            else {
                to_replace = 0;

                if(is_end)
                {
                    cur_nc = nc; cur_nh = nh; 
                    cur_w0 = nw_0; cur_w1= nw_1;
                    cur_ac = ac; cur_uc = uc;
                }
                else
                {
                    cur_nc = nc + ug->u.a[(w>>1)].n;
                    cur_nh = nh + get_path_phasing_weight(w, v, beg, b, ta);
                    get_related_phasing_weight(w>>1, ta, &cur_w0, &cur_w1, s->s);
                    cur_w0 += nw_0; cur_w1 += nw_1;
                    cur_ac = ac + ((s->s[w>>1] == 0)? ug->u.a[(w>>1)].n : 0);
                    cur_uc = uc + ((s->s[w>>1] != 0)? ug->u.a[(w>>1)].n : 0);
                }
                
                cur_weight = cur_nh + MIN(cur_w0, cur_w1) - MAX(cur_w0, cur_w1);
                min_weight = t->nh + MIN(t->w[0], t->w[1]) - MAX(t->w[0], t->w[1]);
                cur_rate = ((cur_ac+cur_uc == 0)? -1 : ((double)(cur_ac)/(double)(cur_ac+cur_uc)));
                max_rate = ((t->ac+t->uc == 0)? -1 : ((double)(t->ac)/(double)(t->ac+t->uc)));

                if(cur_rate > max_rate)
                {
                    to_replace = 1;
                }
                else if(cur_rate == max_rate)
                {
                    if(cur_weight < min_weight)
                    {
                        to_replace = 1;
                    }
                    else if(cur_weight == min_weight)
                    {
                        if(cur_nc > t->nc)
                        {
                            to_replace = 1;
                        }
                        else if(cur_nc == t->nc)
                        {
                            if(d + l > t->d)
                            {
                                to_replace = 1;
                            }
                        }
                    }
                }


                if(to_replace)
                {
                    t->p = v;
                    t->nc = cur_nc;
                    t->nh = cur_nh;
                    t->ac = cur_ac;
                    t->uc = cur_uc;
                    t->w[0] = cur_w0;
                    t->w[1] = cur_w1;
                }


                if (d + l < t->d) t->d = d + l; // update dist
            }

            if (--(t->r) == 0) {
                uint32_t x = get_real_length(g, w, NULL);
                if(x > 0)
                {
                    kv_push(uint32_t, b->S, w);
                }
                else
                {
                    ///at most one tip
                    if(n_tips != 0) goto pop_reset;
                    n_tips++;
                    tip_end = w;
                }
                --n_pending;
            }
        }

        if(n_tips == 1)
        {
            if(tip_end != (uint32_t)-1 && n_pending == 0 && b->S.n == 0)
            {
                ///sink is b.S.a[0]
                kv_push(uint32_t, b->S, tip_end);
                break;
            }
            else
            {
                goto pop_reset;
            }
        }

        if (i < nv || b->S.n == 0) goto pop_reset;
    }while (b->S.n > 1 || n_pending);

    n_pop = 1;
    /**need fix**/
    set_phase_path(b, beg, ta, s);

    pop_reset:

    for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
        bub_p_t *t = &b->a[b->b.a[i]];
        t->p = t->d = t->nc = t->ac = t->uc = t->r = t->s = 0;
        t->nh = t->w[0] = t->w[1] = 0;
    }
    return n_pop;
}

uint32_t get_weightest_node(kv_u_trans_t *ta, bubble_type* bub, ma_ug_t* ug, int8_t *s, uint8_t *vis)
{
    double w_a, w_n, max_w_a, max_w_n;
    uint32_t i, occ, k, m, v, *a = NULL, a_n, e_n, max_w_a_i, max_w_n_i;
    u_trans_t *e = NULL;
    max_w_a = max_w_n = -1; max_w_a_i = max_w_n_i = (uint32_t)-1;
    for (i = 0; i < bub->f_bub; i++)
    {
        get_bubbles(bub, i, NULL, NULL, &a, &a_n, NULL);
        if(vis[i]) continue;
        for (k = 0, w_a = w_n = 0; k < a_n; k++)
        {
            v = a[k]>>1;
            if(s[v] != 0) break; 
            e = u_trans_a(*ta, v);
            e_n = u_trans_n(*ta, v);
            for (m = 0; m < e_n; m++)
            {
                if(s[e[m].tn] != 0)
                {
                    w_a += (e[m].nw>=0?e[m].nw:-e[m].nw);
                }
                else
                {
                    w_n += (e[m].nw>=0?e[m].nw:-e[m].nw);
                }
            }
        }

        if(k >= a_n && a_n > 0)//unset whole bubble
        {
            if(w_a > max_w_a)
            {
                max_w_a_i = i<<1;
                max_w_a = w_a;
            }

            if(w_n > max_w_n)
            {
                max_w_n_i = i<<1;
                max_w_n = w_n;
            }
        }
        else
        {
            for (k = occ = 0; k < a_n; k++)
            {
                v = a[k]>>1;
                if(s[v] != 0)
                {
                    occ++;
                    continue;
                } 
                
                e = u_trans_a(*ta, v);
                e_n = u_trans_n(*ta, v);
                for (m = 0, w_a = w_n = 0; m < e_n; m++)
                {
                    if(s[e[m].tn] != 0)
                    {
                        w_a += (e[m].nw>=0?e[m].nw:-e[m].nw);
                    }
                    else
                    {
                        w_n += (e[m].nw>=0?e[m].nw:-e[m].nw);
                    }
                }

                if(w_a > max_w_a)
                {
                    max_w_a_i = (v<<1)+1;
                    max_w_a = w_a;
                }

                if(w_n > max_w_n)
                {
                    max_w_n_i = (v<<1)+1;
                    max_w_n = w_n;
                }
            }

            if(occ == a_n) vis[i] = 1;
        }
    }

    for (i = 0; i < ug->g->n_seq; i++)
    {
        if(s[i] != 0 || (IF_BUB(i, *bub))) continue;
        v = i;
        e = u_trans_a(*ta, v);
        e_n = u_trans_n(*ta, v);
        for (m = 0, w_a = w_n = 0; m < e_n; m++)
        {
            if(s[e[m].tn] != 0)
            {
                w_a += (e[m].nw>=0?e[m].nw:-e[m].nw);
            }
            else
            {
                w_n += (e[m].nw>=0?e[m].nw:-e[m].nw);
            }
        }

        if(w_a > max_w_a)
        {
            max_w_a_i = (i<<1)+1;
            max_w_a = w_a;
        }

        if(w_n > max_w_n)
        {
            max_w_n_i = (i<<1)+1;
            max_w_n = w_n;
        }
    }
    
    if(max_w_a_i != (uint32_t)-1) return max_w_a_i;
    return max_w_n_i;
}

void init_phase(ha_ug_index* idx, kv_u_trans_t *ta, bubble_type* bub, ps_t *st)
{
    double index_time = yak_realtime();
    uint8_t *vis = NULL; CALLOC(vis, idx->ug->g->n_seq);
    bub_p_t_warp b; memset(&b, 0, sizeof(bub_p_t_warp));
    CALLOC(b.a, idx->ug->g->n_seq*2);
    uint32_t i, k, *a = NULL, n, beg, end;
    memset(st->s, 0, sizeof(int8_t)*idx->ug->g->n_seq);
    
    double z[2];
    u_trans_t *e = NULL;
    uint32_t e_n;
    while(1)
    {
        i = get_weightest_node(ta, bub, idx->ug, st->s, vis);
        if(i == (uint32_t)-1) break;
        if(i&1)
        {
            i>>=1;
            z[0] = z[1] = 0;
            e = u_trans_a(*ta, i);
            e_n = u_trans_n(*ta, i);
            for (k = 0; k < e_n; k++)
            {
                if(e[k].del) continue;
                if(st->s[e[k].tn] > 0) z[0] += e[k].nw;
                else if(st->s[e[k].tn] < 0) z[1] += e[k].nw;
            }
            if(z[0] - z[1] < 0)
            {
                st->s[i] = 1;
            } 
            else if(z[0] - z[1] > 0)
            {
                st->s[i] = -1;
            }
            else
            {
                st->xs = kr_splitmix64(st->xs);
                st->s[i] = st->xs&1? 1 : -1;
            }
        }
        else
        {
            i>>=1;
            get_bubbles(bub, i, &beg, &end, &a, &n, NULL);
            bub_phase(idx->ug, beg, end, &b, ta, st);
            bub_phase(idx->ug, beg, end, &b, ta, st);
        }
    }
    free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
    free(vis);

    for (i = 0; i < idx->ug->g->n_seq; i++)
    {
        if(st->s[i] == 0) fprintf(stderr, "ERROR\n");
        if(u_trans_n(*ta, i) == 0) st->s[i] = 0;
    }

    fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
}


double dfs_weight_hic(uint32_t v, uint8_t* vis_flag, uint8_t* is_vis, kv_u_trans_t *ta, 
kvec_t_u32_warp* stack, kvec_t_u32_warp* result, uint32_t e_flag, uint32_t ava_flag, uint32_t* link_occ)
{
    u_trans_t *e = NULL;
    uint32_t cur, i, next = (uint32_t)-1, e_n;
    stack->a.n = 0;
    kv_push(uint32_t, stack->a, v);
    double w = 0;
    if(link_occ) (*link_occ) = 0;
    while (stack->a.n > 0)
    {
        stack->a.n--;
        cur = stack->a.a[stack->a.n];
        if(is_vis[cur]) continue;
        is_vis[cur] = 1;
        if(cur!=v && vis_flag[cur] != ava_flag) continue;
        e = u_trans_a(*ta, cur); e_n = u_trans_n(*ta, cur);
        for (i = 0; i < e_n; i++)
        {
            if(e[i].del) continue;
            next = e[i].tn;
            if(vis_flag[next]&e_flag)
            {
                w += (e[i].nw >= 0? e[i].nw : -e[i].nw);
                if(link_occ) (*link_occ) += e[i].occ;
                continue;
            }
            if(is_vis[next]) continue;
            if(vis_flag[next] != ava_flag) continue;
            kv_push(uint32_t, stack->a, next);
        }
    }
    return w;
}

double get_chain_weight_hic(bubble_type* bub, ma_ug_t *bub_ug, buf_t* b, uint32_t v, uint32_t convex_source, 
kv_u_trans_t *ta, uint8_t* vis_flag, uint8_t* is_vis, ma_ug_t* ug, kvec_t_u32_warp* stack, 
kvec_t_u32_warp* result, uint32_t e_flag, uint32_t ava_flag, kvec_t_u32_warp* res_utg, uint32_t* link_occ)
{
    long long nodeLen, baseLen, max_stop_nodeLen, max_stop_baseLen;
    ma_utg_t *u = NULL;
    uint32_t convex, k, k_i, k_j, *a, n, beg, sink, uID, root, cur, ncur, n_vx = ug->g->n_seq<<1, occ;
    asg_arc_t *acur = NULL;
    double w = 0;
    b->b.n = 0;
    get_unitig(bub_ug->g, NULL, v, &convex, &nodeLen, &baseLen, &max_stop_nodeLen, 
                                                                    &max_stop_baseLen, 1, b);
    memset(is_vis, 0, n_vx);

    for (k = 0; k < b->b.n; k++)
    {
        u = &(bub_ug->u.a[b->b.a[k]>>1]);
        if(u->n == 0) continue;

        for (k_i = 0; k_i < u->n; k_i++)
        {
            get_bubbles(bub, u->a[k_i]>>33, &beg, &sink, &a, &n, NULL);
            for(k_j = 0; k_j < n; k_j++) is_vis[a[k_j]] = is_vis[a[k_j]^1] = 1;
            if(beg != (uint32_t)-1) is_vis[beg] = is_vis[beg^1] = 1; 
            if(sink != (uint32_t)-1) is_vis[sink] = is_vis[sink^1] = 1; 
        }
    }

    u = &(bub_ug->u.a[v>>1]);
    if((v&1)==0)
    {
        get_bubbles(bub, (u->a[0]>>32)>>1, (((u->a[0]>>32)&1)^1)==1?&root:NULL, 
                                    (((u->a[0]>>32)&1)^1) == 0?&root:NULL, NULL, NULL, NULL);
    } 
    else
    {
        get_bubbles(bub, (u->a[u->n-1]>>32)>>1, ((u->a[u->n-1]>>32)&1)==1?&root:NULL, 
                                    ((u->a[u->n-1]>>32)&1) == 0?&root:NULL, NULL, NULL, NULL);
    }
    
    root ^= 1;
    is_vis[root] = 0; 
    stack->a.n = 0;
    kv_push(uint32_t, stack->a, root);
    while (stack->a.n > 0)///label all untigs not in any chain
    {
        stack->a.n--;
        cur = stack->a.a[stack->a.n];
        if(is_vis[cur]) continue;
        is_vis[cur] = 1;
        if(vis_flag[cur>>1] == 0) vis_flag[cur>>1] = ava_flag;///unitig not in any chain
       
        ncur = asg_arc_n(ug->g, cur);
        acur = asg_arc_a(ug->g, cur);
        for (k = 0; k < ncur; k++)
        {
            if(acur[k].del) continue;
            if(is_vis[acur[k].v]) continue;
            if(vis_flag[acur[k].v>>1] != 0 && vis_flag[acur[k].v>>1] != ava_flag) continue;
            kv_push(uint32_t, stack->a, acur[k].v);
        }
    }


    ///vis_flag keeps isloated nodes
    uint32_t aim_0, aim_1, root_source;
    aim_0 = root>>1;
    u = &(bub_ug->u.a[convex_source>>1]);
    if((convex_source&1)==1)
    {
        get_bubbles(bub, (u->a[0]>>32)>>1, (((u->a[0]>>32)&1)^1)==1?&root_source:NULL, 
                                    (((u->a[0]>>32)&1)^1) == 0?&root_source:NULL, NULL, NULL, NULL);
    } 
    else
    {
        get_bubbles(bub, (u->a[u->n-1]>>32)>>1, ((u->a[u->n-1]>>32)&1)==1?&root_source:NULL, 
                                    ((u->a[u->n-1]>>32)&1) == 0?&root_source:NULL, NULL, NULL, NULL);
    }
    root_source ^= 1;
    aim_1 = root_source>>1;


    cur = root_source;///scan nodes that cannot be reached from root but can be reached from root_source
    ncur = asg_arc_n(ug->g, cur);
    acur = asg_arc_a(ug->g, cur);
    for (k_i = 0; k_i < ncur; k_i++)
    {
        if(acur[k_i].del) continue;
        if(vis_flag[acur[k_i].v>>1] != 0) continue;///skip nodes that are already reachable
        ///don't label any path that can reack other chains
        if_conflict_utg(acur[k_i].v, &aim_0, &aim_1, ug, vis_flag, is_vis, ava_flag, stack);
    }


    for (k = 0; k < ug->g->n_seq; k++)
    {
        if(vis_flag[k] == ava_flag)
        {
            cur = k<<1;
            ncur = asg_arc_n(ug->g, cur);
            acur = asg_arc_a(ug->g, cur);
            for (k_i = 0; k_i < ncur; k_i++)
            {
                if(acur[k_i].del) continue;
                if(vis_flag[acur[k_i].v>>1] != 0) continue;
                if_conflict_utg(acur[k_i].v, &aim_0, &aim_1, ug, vis_flag, is_vis, ava_flag, stack);
            }



            cur = (k<<1)+1;
            ncur = asg_arc_n(ug->g, cur);
            acur = asg_arc_a(ug->g, cur);
            for (k_i = 0; k_i < ncur; k_i++)
            {
                if(acur[k_i].del) continue;
                if(vis_flag[acur[k_i].v>>1] != 0) continue;
                if_conflict_utg(acur[k_i].v, &aim_0, &aim_1, ug, vis_flag, is_vis, ava_flag, stack);
            }
        }
    }
    


    memset(is_vis, 0, n_vx);
    if(link_occ) (*link_occ) = 0;
    for (k = result->a.n = 0, w = 0; k < b->b.n; k++)
    {
        u = &(bub_ug->u.a[b->b.a[k]>>1]);
        if(u->n == 0) continue;
        for (k_i = 0; k_i < u->n; k_i++)
        {
            get_bubbles(bub, u->a[k_i]>>33, NULL, NULL, &a, &n, NULL);

            for (k_j = 0; k_j < n; k_j++)
            {
                uID = a[k_j]>>1;
                w += dfs_weight_hic(uID, vis_flag, is_vis, ta, stack, result, e_flag, ava_flag, &occ);
                if(link_occ) (*link_occ) += occ;
            }
        }
    }


    for (k = 0; k < ug->g->n_seq; k++)
    {
        if(vis_flag[k] == ava_flag)
        {
            vis_flag[k] = 0;
            if(res_utg && (!IF_HOM(k, *bub)))
            {
                kv_push(uint32_t, res_utg->a, k<<1);
            }
        }
    }

    return w;
}


void clean_bubble_chain_by_hic(ma_ug_t* ug, kv_u_trans_t *ta, bubble_type* bub)
{   
    double index_time = yak_realtime();
    ma_ug_t *bs_ug = bub->b_ug;
    uint32_t v, u, i, m, max_i, nv, rv, n_vx, root, flag_pri = 1, flag_aux = 2, flag_ava = 4, occ;
    double w, cutoff = 2;
    uint32_t max_w_occ = 4;
    asg_arc_t *av = NULL;
    n_vx = bs_ug->g->n_seq << 1;
    uint8_t *vis = NULL; CALLOC(vis, ug->g->n_seq<<1);
    uint8_t *is_vis = NULL; CALLOC(is_vis, ug->g->n_seq<<1);
    uint8_t *is_used = NULL; CALLOC(is_used, n_vx);
    uint8_t *dedup = NULL; CALLOC(dedup, ug->g->n_seq<<1);
    buf_t b; memset(&b, 0, sizeof(buf_t));
    kvec_t_u32_warp stack, result, res_utg;
    kv_init(stack.a); kv_init(result.a); kv_init(res_utg.a);
    double *e_w = NULL; MALLOC(e_w, bs_ug->g->n_arc);
    uint32_t *e_occ = NULL, *a_occ = NULL; CALLOC(e_occ, bs_ug->g->n_arc);
    double *aw = NULL, max_w = 0;
    kvec_asg_arc_t_warp edges; kv_init(edges.a);
    ma_ug_t *back_bs_ug = copy_untig_graph(bs_ug);

    for (i = 0; i < bs_ug->g->n_arc; i++)///weight of bs_ug's edges
    {
        e_w[i] = -1;
    }

    for (i = 0; i < bs_ug->g->n_seq; i++)///init all chain with flag_aux 
    {
        set_b_utg_weight_flag(bub, &b, i<<1, vis, flag_aux, NULL);
    }


    for (v = 0; v < n_vx; v++)
    {
        av = asg_arc_a(bs_ug->g, v);
        nv = asg_arc_n(bs_ug->g, v);
        aw = (&e_w[bs_ug->g->idx[v]>>32]);
        a_occ = (&e_occ[bs_ug->g->idx[v]>>32]);
        if(nv <= 1 || get_real_length(bs_ug->g, v, NULL) <= 1) continue;
        set_b_utg_weight_flag_xor(bub, bs_ug, &b, v^1, vis, flag_pri, NULL);

        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            w = get_chain_weight_hic(bub, bs_ug, &b, av[i].v, v, ta, vis, is_vis, ug, &stack, &result, flag_pri, flag_ava, NULL, &occ);
            aw[i] = w;
            a_occ[i] = occ;
        }

        set_b_utg_weight_flag_xor(bub, bs_ug, &b, v^1, vis, flag_pri, NULL);
    }


    for (v = 0; v < n_vx; v++)
    {
        av = asg_arc_a(bs_ug->g, v);
        nv = asg_arc_n(bs_ug->g, v);
        aw = (&e_w[bs_ug->g->idx[v]>>32]);
        a_occ = (&e_occ[bs_ug->g->idx[v]>>32]);
        if(nv <= 1 || get_real_length(bs_ug->g, v, NULL) <= 1) continue;

        for (i = rv = 0, max_i = (uint32_t)-1; i < nv; i++)
        {
            if(av[i].del) continue;
            if(max_i == (uint32_t)-1)
            {
                max_i = i;
                max_w = aw[i];
            }
            else if(max_w < aw[i])
            {
                max_i = i;
                max_w = aw[i];
            }
            rv++;
        }

        if(max_i == (uint32_t)-1) continue;
        ///if(max_w <= max_w_cutoff) continue; //must be <=
        if(a_occ[max_i] <= max_w_occ) continue; //must be <=
        if(rv < 2) continue;

        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            if(i == max_i) continue;
            ///if((av[i].v>>1) == (v>>1) && aw[i] <= max_w_cutoff) continue; ///might be not reasonable
            if((av[i].v>>1) == (v>>1) && a_occ[i] <= max_w_occ) continue; ///might be not reasonable
            if(aw[i]*cutoff < max_w && double_check_bub_branch(&av[i], bs_ug, e_w, e_occ, cutoff, max_w_occ))
            {
                av[i].del = 1; asg_arc_del(bs_ug->g, (av[i].v)^1, (av[i].ul>>32)^1, 1);
            }
        }
    }

    uint32_t rId_0, ori_0, rId_1, ori_1, root_0, root_1, new_bub;
    if(bub->num.n > 0) bub->num.n--;
    new_bub = bub->b_g->n_seq;

    for (v = 0; v < n_vx; v++)
    {
        av = asg_arc_a(bs_ug->g, v);
        nv = asg_arc_n(bs_ug->g, v);
        rv = get_real_length(bs_ug->g, v, NULL);
        if(nv == rv) continue;///no edge drop
        if(rv != 1 || nv <= 1) continue;
        get_real_length(bs_ug->g, v, &u);
        u ^= 1;
        if(get_real_length(bs_ug->g, u, NULL) != 1) continue;
        drop_g_edges_by_utg(bub, bub->b_g, bs_ug, NULL, v, u);
        if(is_used[v] || is_used[u]) continue;

        is_used[v] = is_used[u] = 1;
        root = get_utg_end_from_btg(bub, bs_ug, v);
        rId_0 = root>>1;
        ori_0 = root&1;
        get_bubbles(bub, rId_0, ori_0 == 1?&root_0:NULL, ori_0 == 0?&root_0:NULL, NULL, NULL, NULL);

        root = get_utg_end_from_btg(bub, bs_ug, u);
        rId_1 = root>>1;
        ori_1 = root&1;
        get_bubbles(bub, rId_1, ori_1 == 1?&root_1:NULL, ori_1 == 0?&root_1:NULL, NULL, NULL, NULL);

        res_utg.a.n = 0;
        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, v^1, vis, flag_pri, NULL);
        get_chain_weight_hic(bub, back_bs_ug, &b, u^1, v, ta, vis, is_vis, ug, &stack, &result, flag_pri, flag_ava, &res_utg, NULL);
        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, v^1, vis, flag_pri, NULL);
        for (i = 0; i < res_utg.a.n; i++) dedup[res_utg.a.a[i]>>1] |= 1;

        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, u^1, vis, flag_pri, NULL);
        get_chain_weight_hic(bub, back_bs_ug, &b, v^1, u, ta, vis, is_vis, ug, &stack, &result, flag_pri, flag_ava, &res_utg, NULL);
        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, u^1, vis, flag_pri, NULL);
        for (; i < res_utg.a.n; i++) dedup[res_utg.a.a[i]>>1] |= 2;

        for (i = m = 0; i < res_utg.a.n; i++)
        {
            if(dedup[res_utg.a.a[i]>>1] == 3)
            {
                res_utg.a.a[m] = res_utg.a.a[i];
                m++;
            }
            dedup[res_utg.a.a[i]>>1] = 0;
        }
        res_utg.a.n = m;

        if(!IF_HOM(root_0>>1, *bub)) kv_push(uint32_t, res_utg.a, root_0);
        if(!IF_HOM(root_1>>1, *bub)) kv_push(uint32_t, res_utg.a, root_1);

        update_bubble_graph(&res_utg, root_0^1, rId_0, root_1^1, rId_1, bub, &edges, bub->b_g, NULL, NULL, ug, NULL, 0);

        ///fprintf(stderr, "\n******src-btg%.6ul------>dest-btg%.6ul\n", (v>>1)+1, (u>>1)+1);
    }

    kv_push(uint32_t, bub->num, bub->list.n);
    new_bub = bub->b_g->n_seq - new_bub;
    bub->cross_bub += new_bub;
    ///actually not useful, and may have bug when one bubble at multipe chains
    if(new_bub) update_bub_b_s_idx(bub);
    
    ///debug_tangle_bubble(bub, bub->b_g->n_seq - bub->cross_bub, bub->b_g->n_seq - 1, "Cross-tangle");
    update_bsg(bub->b_g, &edges);

    ma_ug_destroy(bs_ug);
    bs_ug = ma_ug_gen(bub->b_g);
    bub->b_ug = bs_ug;
    kv_destroy(bub->chain_weight);
    ma_utg_t *u_x = NULL;
    bs_ug = bub->b_ug;
    kv_malloc(bub->chain_weight, bs_ug->u.n); bub->chain_weight.n = bs_ug->u.n;
    for (i = 0; i < bs_ug->u.n; i++)
    {
        u_x = &(bs_ug->u.a[i]);
        bub->chain_weight.a[i].id = i;
        // if(u->n <= 1) ///not a chain
        // {
        //     bub->chain_weight.a[i].b_occ = bub->chain_weight.a[i].g_occ = 0;
        //     bub->chain_weight.a[i].del = 1;
        // } 
        // else
        {
            bub->chain_weight.a[i].del = 0;
            calculate_chain_weight(u_x, bub, ug, &(bub->chain_weight.a[i]));
        }
    }
    qsort(bub->chain_weight.a, bub->chain_weight.n, sizeof(chain_w_type), cmp_chain_weight);

    
    
    free(vis); free(is_vis); free(is_used); free(dedup); free(b.b.a); free(e_w); free(e_occ);
    kv_destroy(stack.a); kv_destroy(result.a); kv_destroy(res_utg.a); kv_destroy(edges.a);
    ma_ug_destroy(back_bs_ug);
    fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
}

void clean_sub_tangle(bubble_type* bub, uint8_t* vis_flag, kvec_t_u32_warp *res_utg, 
uint8_t *dedup, uint8_t *is_tangle, kvec_asg_arc_t_warp* edges, uint32_t flag, uint32_t s, 
uint32_t e)
{
    uint32_t i, k_i, pv, v, occ, beg, sink, p_beg, p_sink, n, *a = NULL;
    ma_utg_t *u = NULL;

    for (i = 0, pv = (uint32_t)-1; i < res_utg->a.n; i++)
    {
        if(is_tangle[res_utg->a.a[i]>>1] == 1)
        {
            if(pv == (uint32_t)-1)
            {
                pv = res_utg->a.a[i]>>1;
            }
            else if(pv != (res_utg->a.a[i]>>1))
            {
                pv = (uint32_t)-1;
                break;
            }
        }
        else if(is_tangle[res_utg->a.a[i]>>1] == (uint8_t)-1)
        {
            pv = (uint32_t)-1;
            break;
        }
    }

    if(pv == (uint32_t)-1)
    {
        for (i = 0; i < res_utg->a.n; i++)
        {
            is_tangle[res_utg->a.a[i]>>1] = (uint8_t)-1;
        }
        return;
    } 

    occ = 0;
    u = &(bub->b_ug->u.a[s]);
    for (i = 0, p_beg = p_sink = (uint32_t)-1; i < u->n; i++)
    {
        get_bubbles(bub, u->a[i]>>33, &beg, &sink, &a, &n, NULL);

        for (k_i = 0; k_i < n; k_i++)
        {
            occ += bub->ug->u.a[a[k_i]>>1].n;
        }
        if(beg != (uint32_t)-1 && (beg>>1) != (p_beg>>1) && (beg>>1) != (p_sink>>1))
        {
            occ += bub->ug->u.a[beg>>1].n; 
        } 
        if(sink != (uint32_t)-1 && (sink>>1) != (p_beg>>1) && (sink>>1) != (p_sink>>1))
        {
            occ += bub->ug->u.a[sink>>1].n; 
        } 
        p_beg = beg; p_sink = sink;
    }

    u = &(bub->b_ug->u.a[e]);
    for (i = 0, p_beg = p_sink = (uint32_t)-1; i < u->n; i++)
    {
        get_bubbles(bub, u->a[i]>>33, &beg, &sink, &a, &n, NULL);

        for (k_i = 0; k_i < n; k_i++)
        {
            occ += bub->ug->u.a[a[k_i]>>1].n;
        }
        if(beg != (uint32_t)-1 && (beg>>1) != (p_beg>>1) && (beg>>1) != (p_sink>>1))
        {
            occ += bub->ug->u.a[beg>>1].n; 
        } 
        if(sink != (uint32_t)-1 && (sink>>1) != (p_beg>>1) && (sink>>1) != (p_sink>>1))
        {
            occ += bub->ug->u.a[sink>>1].n; 
        } 
        p_beg = beg; p_sink = sink;
    }

    if(occ <= bub->ug->u.a[pv].n*200)
    {
        for (i = 0; i < res_utg->a.n; i++)
        {
            is_tangle[res_utg->a.a[i]>>1] = (uint8_t)-1;
        }
        return;
    }

    if(!edges) return;

    u = &(bub->b_ug->u.a[s]);
    for (i = 0; i < u->n; i++)
    {
        get_bubbles(bub, u->a[i]>>33, &beg, &sink, &a, &n, NULL);

        for (k_i = 0; k_i < n; k_i++)
        {
            vis_flag[a[k_i]>>1] ^= flag;
        }
        if(beg != (uint32_t)-1)
        {
            vis_flag[beg>>1] ^= flag;
        } 
        if(sink != (uint32_t)-1)
        {
            vis_flag[sink>>1] ^= flag;
        } 
    }
    u = &(bub->b_ug->u.a[e]);
    for (i = 0; i < u->n; i++)
    {
        get_bubbles(bub, u->a[i]>>33, &beg, &sink, &a, &n, NULL);

        for (k_i = 0; k_i < n; k_i++)
        {
            vis_flag[a[k_i]>>1] ^= flag;
        }
        if(beg != (uint32_t)-1)
        {
            vis_flag[beg>>1] ^= flag;
        } 
        if(sink != (uint32_t)-1)
        {
            vis_flag[sink>>1] ^= flag;
        } 
    }
    for (i = 0; i < res_utg->a.n; i++)
    {
        if(dedup[res_utg->a.a[i]>>1] == 3)
        {
            vis_flag[res_utg->a.a[i]>>1] ^= flag;
        }
    }



    asg_arc_t *av = NULL, *p = NULL;
    uint32_t nv, c[2];

    c[0] = c[1] = 0;
    while (1)
    {
        c[0] = c[1] = 0; p = NULL;
        v = pv<<1;
        av = asg_arc_a(bub->ug->g, v);
        nv = asg_arc_n(bub->ug->g, v);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            c[(!!(vis_flag[av[i].v>>1]&flag))]++;
            if(p == NULL || p->ol > av[i].ol)
            {
                p = &(av[i]);
            }
        }
        



        v = (pv<<1) + 1;
        av = asg_arc_a(bub->ug->g, v);
        nv = asg_arc_n(bub->ug->g, v);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            c[(!!(vis_flag[av[i].v>>1]&flag))]++;
            if(p == NULL || p->ol > av[i].ol)
            {
                p = &(av[i]);
            }
        }

        if(p)
        {
            p->del = 1;
            c[(!!(vis_flag[p->v>>1]&flag))]--;
            if((p->v>>1) == pv)
            {
                asg_arc_del(bub->ug->g, (p->v)^1, (p->ul>>32)^1, 1);
                c[0] = c[1] = 0;
                v = pv<<1;
                av = asg_arc_a(bub->ug->g, v);
                nv = asg_arc_n(bub->ug->g, v);
                for (i = 0; i < nv; i++)
                {
                    if(av[i].del) continue;
                    c[(!!(vis_flag[av[i].v>>1]&flag))]++;
                }

                v = (pv<<1) + 1;
                av = asg_arc_a(bub->ug->g, v);
                nv = asg_arc_n(bub->ug->g, v);
                for (i = 0; i < nv; i++)
                {
                    if(av[i].del) continue;
                    c[(!!(vis_flag[av[i].v>>1]&flag))]++;
                }
            }

            
        }
        if(c[0] == 0 || c[1] == 0) break;
    }
    


    
    v = pv<<1;
    av = asg_arc_a(bub->ug->g, v);
    nv = asg_arc_n(bub->ug->g, v);
    for (i = 0; i < nv; i++)
    {
        if(!av[i].del) continue;
        av[i].del = 0;
        kv_push(asg_arc_t, edges->a, av[i]);
    }
    

    v = (pv<<1) + 1;
    av = asg_arc_a(bub->ug->g, v);
    nv = asg_arc_n(bub->ug->g, v);
    for (i = 0; i < nv; i++)
    {
        if(!av[i].del) continue;
        av[i].del = 0;
        kv_push(asg_arc_t, edges->a, av[i]);
    }

    is_tangle[pv] = (uint8_t)-1;


    u = &(bub->b_ug->u.a[s]);
    for (i = 0; i < u->n; i++)
    {
        get_bubbles(bub, u->a[i]>>33, &beg, &sink, &a, &n, NULL);

        for (k_i = 0; k_i < n; k_i++)
        {
            vis_flag[a[k_i]>>1] ^= flag;
        }
        if(beg != (uint32_t)-1)
        {
            vis_flag[beg>>1] ^= flag;
        } 
        if(sink != (uint32_t)-1)
        {
            vis_flag[sink>>1] ^= flag;
        } 
    }
    u = &(bub->b_ug->u.a[e]);
    for (i = 0; i < u->n; i++)
    {
        get_bubbles(bub, u->a[i]>>33, &beg, &sink, &a, &n, NULL);

        for (k_i = 0; k_i < n; k_i++)
        {
            vis_flag[a[k_i]>>1] ^= flag;
        }
        if(beg != (uint32_t)-1)
        {
            vis_flag[beg>>1] ^= flag;
        } 
        if(sink != (uint32_t)-1)
        {
            vis_flag[sink>>1] ^= flag;
        } 
    }
    for (i = 0; i < res_utg->a.n; i++)
    {
        if(dedup[res_utg->a.a[i]>>1] == 3)
        {
            vis_flag[res_utg->a.a[i]>>1] ^= flag;
        }
    }


}

void delete_sg_e_by_ug(asg_t* rg, ma_ug_t* ug, uint32_t v, uint32_t w)
{
    uint32_t vx, wx;
    vx = (v&1?((ug->u.a[v>>1].a[0]>>32)^1):(ug->u.a[v>>1].a[ug->u.a[v>>1].n-1]>>32));
    wx = (w&1?((ug->u.a[w>>1].a[ug->u.a[w>>1].n-1]>>32)^1):(ug->u.a[w>>1].a[0]>>32));
    asg_arc_del(rg, vx, wx, 1); asg_arc_del(rg, wx^1, vx^1, 1);
}

void resolve_bubble_chain_by_hic(ha_ug_index *idx, kv_u_trans_t *ta, bubble_type* bub)
{   
    // double index_time = yak_realtime();
    ma_ug_t* ug = idx->ug;
    ma_ug_t *bs_ug = bub->b_ug;
    uint32_t v, u, i, max_i, nv, rv, n_vx, root, flag_pri = 1, flag_aux = 2, flag_ava = 4, occ;
    double w, cutoff = 2;
    uint32_t max_w_occ = 4;
    asg_arc_t *av = NULL;
    n_vx = bs_ug->g->n_seq << 1;
    uint8_t *vis = NULL; CALLOC(vis, ug->g->n_seq<<1);
    uint8_t *is_vis = NULL; CALLOC(is_vis, ug->g->n_seq<<1);
    uint8_t *is_used = NULL; CALLOC(is_used, n_vx);
    uint8_t *dedup = NULL; CALLOC(dedup, ug->g->n_seq<<1);
    uint8_t *is_tangle = NULL; CALLOC(is_tangle, ug->g->n_seq);
    buf_t b; memset(&b, 0, sizeof(buf_t));
    kvec_t_u32_warp stack, result, res_utg;
    kv_init(stack.a); kv_init(result.a); kv_init(res_utg.a);
    double *e_w = NULL; MALLOC(e_w, bs_ug->g->n_arc);
    uint32_t *e_occ = NULL, *a_occ = NULL; CALLOC(e_occ, bs_ug->g->n_arc);
    double *aw = NULL, max_w = 0;
    kvec_asg_arc_t_warp edges; kv_init(edges.a);
    ma_ug_t *back_bs_ug = copy_untig_graph(bs_ug);

    for (i = 0; i < bs_ug->g->n_arc; i++)///weight of bs_ug's edges
    {
        e_w[i] = -1;
    }

    for (i = 0; i < bs_ug->g->n_seq; i++)///init all chain with flag_aux 
    {
        set_b_utg_weight_flag(bub, &b, i<<1, vis, flag_aux, NULL);
    }


    for (v = 0; v < n_vx; v++)
    {
        av = asg_arc_a(bs_ug->g, v);
        nv = asg_arc_n(bs_ug->g, v);
        aw = (&e_w[bs_ug->g->idx[v]>>32]);
        a_occ = (&e_occ[bs_ug->g->idx[v]>>32]);
        if(nv <= 1 || get_real_length(bs_ug->g, v, NULL) <= 1) continue;
        set_b_utg_weight_flag_xor(bub, bs_ug, &b, v^1, vis, flag_pri, NULL);

        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            w = get_chain_weight_hic(bub, bs_ug, &b, av[i].v, v, ta, vis, is_vis, ug, &stack, &result, flag_pri, flag_ava, NULL, &occ);
            aw[i] = w;
            a_occ[i] = occ;
        }

        set_b_utg_weight_flag_xor(bub, bs_ug, &b, v^1, vis, flag_pri, NULL);
    }


    for (v = 0; v < n_vx; v++)
    {
        av = asg_arc_a(bs_ug->g, v);
        nv = asg_arc_n(bs_ug->g, v);
        aw = (&e_w[bs_ug->g->idx[v]>>32]);
        a_occ = (&e_occ[bs_ug->g->idx[v]>>32]);
        if(nv <= 1 || get_real_length(bs_ug->g, v, NULL) <= 1) continue;

        for (i = rv = 0, max_i = (uint32_t)-1; i < nv; i++)
        {
            if(av[i].del) continue;
            if(max_i == (uint32_t)-1)
            {
                max_i = i;
                max_w = aw[i];
            }
            else if(max_w < aw[i])
            {
                max_i = i;
                max_w = aw[i];
            }
            rv++;
        }

        if(max_i == (uint32_t)-1) continue;
        ///if(max_w <= max_w_cutoff) continue; //must be <=
        if(a_occ[max_i] <= max_w_occ) continue; //must be <=
        if(rv < 2) continue;

        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            if(i == max_i) continue;
            ///if((av[i].v>>1) == (v>>1) && aw[i] <= max_w_cutoff) continue; ///might be not reasonable
            if((av[i].v>>1) == (v>>1) && a_occ[i] <= max_w_occ) continue; ///might be not reasonable
            if(aw[i]*cutoff < max_w && double_check_bub_branch(&av[i], bs_ug, e_w, e_occ, cutoff, max_w_occ))
            {
                av[i].del = 1; asg_arc_del(bs_ug->g, (av[i].v)^1, (av[i].ul>>32)^1, 1);
            }
        }
    }

    uint32_t rId_0, ori_0, rId_1, ori_1, root_0, root_1;
    for (v = 0; v < n_vx; v++)
    {
        av = asg_arc_a(bs_ug->g, v);
        nv = asg_arc_n(bs_ug->g, v);
        rv = get_real_length(bs_ug->g, v, NULL);
        if(nv == rv) continue;///no edge drop
        if(rv != 1 || nv <= 1) continue;
        get_real_length(bs_ug->g, v, &u);
        u ^= 1;
        if(get_real_length(bs_ug->g, u, NULL) != 1) continue;
        drop_g_edges_by_utg(bub, bub->b_g, bs_ug, NULL, v, u);
        if(is_used[v] || is_used[u]) continue;

        is_used[v] = is_used[u] = 1;
        root = get_utg_end_from_btg(bub, bs_ug, v);
        rId_0 = root>>1;
        ori_0 = root&1;
        get_bubbles(bub, rId_0, ori_0 == 1?&root_0:NULL, ori_0 == 0?&root_0:NULL, NULL, NULL, NULL);

        root = get_utg_end_from_btg(bub, bs_ug, u);
        rId_1 = root>>1;
        ori_1 = root&1;
        get_bubbles(bub, rId_1, ori_1 == 1?&root_1:NULL, ori_1 == 0?&root_1:NULL, NULL, NULL, NULL);

        res_utg.a.n = 0;
        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, v^1, vis, flag_pri, NULL);
        get_chain_weight_hic(bub, back_bs_ug, &b, u^1, v, ta, vis, is_vis, ug, &stack, &result, flag_pri, flag_ava, &res_utg, NULL);
        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, v^1, vis, flag_pri, NULL);
        for (i = 0; i < res_utg.a.n; i++) dedup[res_utg.a.a[i]>>1] |= 1;

        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, u^1, vis, flag_pri, NULL);
        get_chain_weight_hic(bub, back_bs_ug, &b, v^1, u, ta, vis, is_vis, ug, &stack, &result, flag_pri, flag_ava, &res_utg, NULL);
        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, u^1, vis, flag_pri, NULL);
        for (; i < res_utg.a.n; i++) dedup[res_utg.a.a[i]>>1] |= 2;

        for (i = 0; i < res_utg.a.n; i++)
        {
            if(dedup[res_utg.a.a[i]>>1] == 3)
            {
                is_tangle[res_utg.a.a[i]>>1] = 3;
            }
            else
            {
                if(is_tangle[res_utg.a.a[i]>>1] != 3)
                {
                    is_tangle[res_utg.a.a[i]>>1] = 1;
                }
            }
            dedup[res_utg.a.a[i]>>1] = 0;
        }
    }

    /*******************************for debug************************************/
    // for (i = 0; i < ug->g->n_seq; i++)
    // {
    //     if(is_tangle[i] == 1)
    //     {
    //         fprintf(stderr, "*****tangle-utg%.6ul\n", i+1);
    //     }
    // }
    /*******************************for debug************************************/
    memset(is_used, 0, n_vx);
    for (v = 0; v < n_vx; v++)
    {
        av = asg_arc_a(bs_ug->g, v);
        nv = asg_arc_n(bs_ug->g, v);
        rv = get_real_length(bs_ug->g, v, NULL);
        if(nv == rv) continue;///no edge drop
        if(rv != 1 || nv <= 1) continue;
        get_real_length(bs_ug->g, v, &u);
        u ^= 1;
        if(get_real_length(bs_ug->g, u, NULL) != 1) continue;
        drop_g_edges_by_utg(bub, bub->b_g, bs_ug, NULL, v, u);
        if(is_used[v] || is_used[u]) continue;

        is_used[v] = is_used[u] = 1;
        root = get_utg_end_from_btg(bub, bs_ug, v);
        rId_0 = root>>1;
        ori_0 = root&1;
        get_bubbles(bub, rId_0, ori_0 == 1?&root_0:NULL, ori_0 == 0?&root_0:NULL, NULL, NULL, NULL);

        root = get_utg_end_from_btg(bub, bs_ug, u);
        rId_1 = root>>1;
        ori_1 = root&1;
        get_bubbles(bub, rId_1, ori_1 == 1?&root_1:NULL, ori_1 == 0?&root_1:NULL, NULL, NULL, NULL);

        res_utg.a.n = 0;
        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, v^1, vis, flag_pri, NULL);
        get_chain_weight_hic(bub, back_bs_ug, &b, u^1, v, ta, vis, is_vis, ug, &stack, &result, flag_pri, flag_ava, &res_utg, NULL);
        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, v^1, vis, flag_pri, NULL);
        for (i = 0; i < res_utg.a.n; i++) dedup[res_utg.a.a[i]>>1] |= 1;

        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, u^1, vis, flag_pri, NULL);
        get_chain_weight_hic(bub, back_bs_ug, &b, v^1, u, ta, vis, is_vis, ug, &stack, &result, flag_pri, flag_ava, &res_utg, NULL);
        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, u^1, vis, flag_pri, NULL);
        for (; i < res_utg.a.n; i++) dedup[res_utg.a.a[i]>>1] |= 2;


        /*******************************for debug************************************/
        // fprintf(stderr, "\nres_utg.a.n: %u, beg-utg%.6ul, sink-utg%.6ul\n", 
        //             (uint32_t)res_utg.a.n, (root_0>>1)+1, (root_1>>1)+1);
        // for (i = 0; i < res_utg.a.n; i++)
        // {
        //     if(dedup[res_utg.a.a[i]>>1] == 3)
        //     {
        //         fprintf(stderr, "share-utg%.6ul\n", (res_utg.a.a[i]>>1)+1);
        //     }
        // }
        // for (i = 0; i < res_utg.a.n; i++)
        // {
        //     if(dedup[res_utg.a.a[i]>>1] == 1)
        //     {
        //         fprintf(stderr, "1-utg%.6ul\n", (res_utg.a.a[i]>>1)+1);
        //     }
        // }
        // for (i = 0; i < res_utg.a.n; i++)
        // {
        //     if(dedup[res_utg.a.a[i]>>1] == 2)
        //     {
        //         fprintf(stderr, "2-utg%.6ul\n", (res_utg.a.a[i]>>1)+1);
        //     }
        // }
        /*******************************for debug************************************/

        clean_sub_tangle(bub, vis, &res_utg, dedup, is_tangle, NULL, flag_pri, v>>1, u>>1);

        for (i = 0; i < res_utg.a.n; i++)
        {
            dedup[res_utg.a.a[i]>>1] = 0;
        }
    }


    edges.a.n = 0;
    memset(is_used, 0, n_vx);
    for (v = 0; v < n_vx; v++)
    {
        av = asg_arc_a(bs_ug->g, v);
        nv = asg_arc_n(bs_ug->g, v);
        rv = get_real_length(bs_ug->g, v, NULL);
        if(nv == rv) continue;///no edge drop
        if(rv != 1 || nv <= 1) continue;
        get_real_length(bs_ug->g, v, &u);
        u ^= 1;
        if(get_real_length(bs_ug->g, u, NULL) != 1) continue;
        drop_g_edges_by_utg(bub, bub->b_g, bs_ug, NULL, v, u);
        if(is_used[v] || is_used[u]) continue;

        is_used[v] = is_used[u] = 1;
        root = get_utg_end_from_btg(bub, bs_ug, v);
        rId_0 = root>>1;
        ori_0 = root&1;
        get_bubbles(bub, rId_0, ori_0 == 1?&root_0:NULL, ori_0 == 0?&root_0:NULL, NULL, NULL, NULL);

        root = get_utg_end_from_btg(bub, bs_ug, u);
        rId_1 = root>>1;
        ori_1 = root&1;
        get_bubbles(bub, rId_1, ori_1 == 1?&root_1:NULL, ori_1 == 0?&root_1:NULL, NULL, NULL, NULL);

        res_utg.a.n = 0;
        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, v^1, vis, flag_pri, NULL);
        get_chain_weight_hic(bub, back_bs_ug, &b, u^1, v, ta, vis, is_vis, ug, &stack, &result, flag_pri, flag_ava, &res_utg, NULL);
        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, v^1, vis, flag_pri, NULL);
        for (i = 0; i < res_utg.a.n; i++) dedup[res_utg.a.a[i]>>1] |= 1;

        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, u^1, vis, flag_pri, NULL);
        get_chain_weight_hic(bub, back_bs_ug, &b, v^1, u, ta, vis, is_vis, ug, &stack, &result, flag_pri, flag_ava, &res_utg, NULL);
        set_b_utg_weight_flag_xor(bub, back_bs_ug, &b, u^1, vis, flag_pri, NULL);
        for (; i < res_utg.a.n; i++) dedup[res_utg.a.a[i]>>1] |= 2;

        clean_sub_tangle(bub, vis, &res_utg, dedup, is_tangle, &edges, flag_pri, v>>1, u>>1);

        for (i = 0; i < res_utg.a.n; i++)
        {
            dedup[res_utg.a.a[i]>>1] = 0;
        }
    }


    /*******************************for debug************************************/
    // for (i = 0; i < ug->g->n_seq; i++)
    // {
    //     if(is_tangle[i] == 1)
    //     {
    //         fprintf(stderr, "####tangle-utg%.6ul\n", i+1);
    //     }
    // }

    // for (i = 0; i < edges.a.n; i++)
    // {
    //     fprintf(stderr, "s-utg%.6lul<------>d-utg%.6ul\n", (edges.a.a[i].ul>>33) + 1, (edges.a.a[i].v>>1) + 1);
        
    // }
    /*******************************for debug************************************/
    if(edges.a.n > 0)
    {
        for (i = 0; i < edges.a.n; i++)
        {
            asg_arc_del(bub->ug->g, edges.a.a[i].ul>>32, edges.a.a[i].v, 1);
            asg_arc_del(bub->ug->g, (edges.a.a[i].v)^1, (edges.a.a[i].ul>>32)^1, 1);
            delete_sg_e_by_ug(idx->read_g, idx->ug, edges.a.a[i].ul>>32, edges.a.a[i].v);
            delete_sg_e_by_ug(idx->read_g, idx->ug, (edges.a.a[i].v)^1, (edges.a.a[i].ul>>32)^1);
        }
        asg_cleanup(bub->ug->g);
        asg_cleanup(idx->read_g);
    }
    

    free(vis); free(is_vis); free(is_used); free(dedup); free(b.b.a); free(e_w); free(e_occ); free(is_tangle);
    kv_destroy(stack.a); kv_destroy(result.a); kv_destroy(res_utg.a); kv_destroy(edges.a);
    ma_ug_destroy(back_bs_ug);
    // fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
}

void append_boundary_chain_hic(ma_ug_t* ug, kv_u_trans_t *ta, bubble_type* bub)
{
    // double index_time = yak_realtime();
    ma_ug_t *bs_ug = bub->b_ug;
    uint32_t v, u, i, k, beg_idx, m, nv, n_vx, flag_pri = 1, flag_aux = 2, flag_ava = 4;
    uint32_t root, rId_0, ori_0, root_0, new_bub;
    asg_arc_t *av = NULL;
    n_vx = bs_ug->g->n_seq << 1;
    uint8_t *vis = NULL; CALLOC(vis, ug->g->n_seq<<1);
    uint8_t *is_vis = NULL; CALLOC(is_vis, ug->g->n_seq<<1);
    uint8_t *is_used = NULL; CALLOC(is_used, n_vx);
    uint8_t *dedup = NULL; CALLOC(dedup, ug->g->n_seq<<1);
    buf_t b; memset(&b, 0, sizeof(buf_t));
    kvec_t_u32_warp stack, result, res_utg;
    kv_init(stack.a); kv_init(result.a); kv_init(res_utg.a);
    kvec_asg_arc_t_warp edges; kv_init(edges.a);

    for (i = 0; i < bs_ug->g->n_seq; i++)
    {
        set_b_utg_weight_flag(bub, &b, i<<1, vis, flag_aux, NULL);
    }

    if(bub->num.n > 0) bub->num.n--;
    new_bub = bub->b_g->n_seq;
    for (v = 0; v < n_vx; v++)
    {
        av = asg_arc_a(bs_ug->g, v);
        nv = asg_arc_n(bs_ug->g, v);
        if(nv == 0 || get_real_length(bs_ug->g, v, NULL) == 0) continue;
        res_utg.a.n = 0;
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            u = av[i].v^1;
            beg_idx = res_utg.a.n;
            set_b_utg_weight_flag_xor(bub, bs_ug, &b, u^1, vis, flag_pri, NULL);
            get_chain_weight_hic(bub, bs_ug, &b, v^1, u, ta, vis, is_vis, ug, &stack, &result, flag_pri, flag_ava, &res_utg, NULL);
            set_b_utg_weight_flag_xor(bub, bs_ug, &b, u^1, vis, flag_pri, NULL);
            for (k = m = beg_idx; k < res_utg.a.n; k++)
            {
                if(dedup[res_utg.a.a[k]>>1] != 0) continue;
                dedup[res_utg.a.a[k]>>1] = 1;
                res_utg.a.a[m] = res_utg.a.a[k];
                m++;
            }
            res_utg.a.n = m;
        }

        for (k = 0; k < res_utg.a.n; k++) dedup[res_utg.a.a[k]>>1] = 0;

        root = get_utg_end_from_btg(bub, bs_ug, v);
        rId_0 = root>>1;
        ori_0 = root&1;
        get_bubbles(bub, rId_0, ori_0 == 1?&root_0:NULL, ori_0 == 0?&root_0:NULL, NULL, NULL, NULL);

        if(root_0 != (uint32_t)-1 && (!IF_HOM(root_0>>1, *bub))) kv_push(uint32_t, res_utg.a, root_0);

        if(v&1)
        {
            update_bubble_graph(&res_utg, root_0^1, rId_0, (uint32_t)-1, (uint32_t)-1, bub, &edges, bub->b_g, NULL, NULL, ug, NULL, 0);
        }
        else
        {
            update_bubble_graph(&res_utg, (uint32_t)-1, (uint32_t)-1, root_0^1, rId_0, bub, &edges, bub->b_g, NULL, NULL, ug, NULL, 0);
        }
    }
    kv_push(uint32_t, bub->num, bub->list.n);
    new_bub = bub->b_g->n_seq - new_bub;
    bub->mess_bub += new_bub;
    ///actually not useful, and may have bug when one bubble at multipe chains
    if(new_bub) update_bub_b_s_idx(bub);


    for (v = 0; v < n_vx; v++)
    {
        av = asg_arc_a(bs_ug->g, v);
        nv = asg_arc_n(bs_ug->g, v);
        if(nv == 0 || get_real_length(bs_ug->g, v, NULL) == 0) continue;
        drop_g_edges_by_utg(bub, bub->b_g, bs_ug, NULL, v, (uint32_t)-1);
    }

    update_bsg(bub->b_g, &edges);

    ma_ug_destroy(bs_ug);
    bs_ug = ma_ug_gen(bub->b_g);
    bub->b_ug = bs_ug;
    kv_destroy(bub->chain_weight);
    ma_utg_t *u_x = NULL;
    bs_ug = bub->b_ug;
    kv_malloc(bub->chain_weight, bs_ug->u.n); bub->chain_weight.n = bs_ug->u.n;
    for (i = 0; i < bs_ug->u.n; i++)
    {
        u_x = &(bs_ug->u.a[i]);
        bub->chain_weight.a[i].id = i;
        // if(u->n <= 1) ///not a chain
        // {
        //     bub->chain_weight.a[i].b_occ = bub->chain_weight.a[i].g_occ = 0;
        //     bub->chain_weight.a[i].del = 1;
        // } 
        // else
        {
            bub->chain_weight.a[i].del = 0;
            calculate_chain_weight(u_x, bub, ug, &(bub->chain_weight.a[i]));
        }
    }
    qsort(bub->chain_weight.a, bub->chain_weight.n, sizeof(chain_w_type), cmp_chain_weight);



    free(vis); free(is_vis); free(is_used); free(dedup); free(b.b.a);
    kv_destroy(stack.a); kv_destroy(result.a); kv_destroy(res_utg.a);
    kv_destroy(edges.a);
    // fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
}

double get_specific_hic_weight_by_chain(uint32_t uid, kv_u_trans_t *ta, uint8_t* vis, uint8_t flag)
{
    u_trans_t *e = u_trans_a(*ta, uid);
    uint64_t k, e_n = u_trans_n(*ta, uid);
    double w = 0;
    
    for (k = 0; k < e_n; k++)
    {
        if(e[k].del) continue;
        if(vis[e[k].tn] != flag) continue;
        w += (e[k].nw>=0?e[k].nw:-e[k].nw);
    }
    return w;
}

void update_bubble_weight(bub_sort_vec* w_stack, uint32_t idx, kv_u_trans_t *ta, uint8_t* vis, uint32_t flag_cur)
{
    w_stack->a[idx].used = 1;
    uint32_t k, m, uid = w_stack->a[idx].p_id>>1;
    u_trans_t *e = u_trans_a(*ta, uid);
    uint32_t e_n = u_trans_n(*ta, uid);
    for (k = 0; k < e_n; k++)
    {
        if(e[k].del) continue;
        if(vis[e[k].tn] != flag_cur) continue;
        for (m = 0; m < w_stack->n; m++)
        {
            if(e[k].tn == (w_stack->a[m].p_id>>1)) break;
        }
        if(m >= w_stack->n) continue;
        if(w_stack->a[m].used) continue;
        w_stack->a[m].weight += (e[k].nw>=0?e[k].nw:-e[k].nw);
    }
}

void reorder_bubbble_chain(kv_u_trans_t *ta, bubble_type* bub, bub_sort_vec* w_stack, 
uint8_t* vis, uint32_t n_utg, uint32_t chain_id)
{
    uint32_t i, k, m, max_idx, flag_cur = 3, flag_right = 2, flag_left = 1, flag_unset = 0, *a, n;
    uint64_t bid, uid;
    ma_utg_t *u = &(bub->b_ug->u.a[chain_id]);
    memset(vis, flag_unset, n_utg);
    for (i = 0; i < u->n; i++)
    {
        bid = u->a[i]>>33;
        get_bubbles(bub, bid, NULL, NULL, &a, &n, NULL);
        for (k = 0; k < n; k++)
        {
            uid = a[k]>>1;
            vis[uid] = flag_right;
        }
    }

    for (i = 0; i < u->n; i++)
    {
        w_stack->n = 0;
        bid = u->a[i]>>33;
        get_bubbles(bub, bid, NULL, NULL, &a, &n, NULL);
        kv_resize(bub_sort_type, *w_stack, n); 
        w_stack->n = n;
        for (k = 0; k < n; k++)
        {
            uid = a[k]>>1;
            vis[uid] = flag_cur;
            w_stack->a[k].p_id = a[k];
            w_stack->a[k].weight = 0;
            w_stack->a[k].used = 0;
        }
        
        for (k = 0; k < w_stack->n; k++)
        {
            w_stack->a[k].weight += get_specific_hic_weight_by_chain(w_stack->a[k].p_id>>1, 
            ta, vis, flag_left);
            w_stack->a[k].weight -= get_specific_hic_weight_by_chain(w_stack->a[k].p_id>>1, 
            ta, vis, flag_right);
        }
        m = 0;
        while ((max_idx = get_max_hap_g(w_stack, NULL)) != (uint32_t)-1)
        {
            a[m] = w_stack->a[max_idx].p_id;
            m++;
            update_bubble_weight(w_stack, max_idx, ta, vis, flag_cur);
        }

        while ((max_idx = get_max_hap_g(w_stack, &max_idx)) != (uint32_t)-1)
        {
            a[m] = w_stack->a[max_idx].p_id;
            m++;
            update_bubble_weight(w_stack, max_idx, ta, vis, flag_cur);
        }
        if(m != n) fprintf(stderr, "ERROR\n");


        for (k = 0; k < n; k++)
        {
            uid = a[k]>>1;
            vis[uid] = flag_left;
        }
    }
}

void reorder_bubbles(bubble_type* bub, kv_u_trans_t *ta, uint32_t n_utg)
{
    // double index_time = yak_realtime();
    uint8_t* vis = NULL; MALLOC(vis, n_utg);
    bub_sort_vec w_stack; kv_init(w_stack);
    uint32_t i;

    for (i = 0; i < bub->chain_weight.n; i++)
    {
        if(bub->chain_weight.a[i].del) continue;
        reorder_bubbble_chain(ta, bub, &w_stack, vis, n_utg, bub->chain_weight.a[i].id);
    }

    free(vis); kv_destroy(w_stack);
    // fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
}

void update_trans_g(ha_ug_index* idx, kv_u_trans_t *ta, bubble_type* bub)
{
    // double index_time = yak_realtime();
    // update_bubble_chain(idx->ug, bub, 0, 1);

    // resolve_bubble_chain_tangle(idx->ug, bub);

    clean_bubble_chain_by_hic(idx->ug, ta, bub);

    // print_debug_bubble_graph(bub, idx->ug, "bub-2");

    // append_boundary_chain_hic(idx->ug, ta, bub);

    // fprintf(stderr, "s_bub: %lu, f_bub: %lu, b_bub: %lu, b_end_bub: %lu, tangle_bub: %lu, cross_bub: %lu, mess_bub: %lu\n", 
    // bub->s_bub, bub->f_bub, bub->b_bub, bub->b_end_bub, bub->tangle_bub, bub->cross_bub, bub->mess_bub);

    ///reorder_bubbles(bub, ta, idx->ug->g->n_seq);
    // fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
}

uint32_t get_max_unitig(H_partition* h, G_partition* g_p, hc_links* link, bubble_type* bub)
{
    double min, weight;
    uint32_t i, min_i;
     
    for (i = 0, min = 1, min_i = (uint32_t)-1; i < g_p->n; i++)
    {
        if(h->lock[i]) continue;
        weight = 0;
        if(g_p->a[i].h[0] > 0 && (g_p->a[i].status[0] == 1 || g_p->a[i].status[0] == -1))
        {
            weight += (g_p->a[i].weight[0] * g_p->a[i].status[0]);
        }

        if(g_p->a[i].h[1] > 0 && (g_p->a[i].status[1] == 1 || g_p->a[i].status[1] == -1))
        {
            weight += (g_p->a[i].weight[1] * g_p->a[i].status[1]);
        }
        weight += g_p->a[i].weight_convex*2;

        if(weight >= 0) continue;
        if(weight < min)
        {
            min = weight;
            min_i = i;
        }
    }
    ///fprintf(stderr, "*****************min: %f\n", min);
    return min_i;
}

double get_cluster_weight_debug(G_partition* g_p, hc_links* link, uint32_t *h, uint32_t h_n)
{
    int o_d = 0;
    double weight = 0;
    uint32_t j, k, m, uID;
    for (j = 0, weight = 0; j < h_n; j++)
    {
        for (k = 0; k < link->a.a[h[j]].e.n; k++)
        {
            if(link->a.a[h[j]].e.a[k].del) continue;
            for (m = 0; m < h_n; m++)
            {
                if(h[m] == link->a.a[h[j]].e.a[k].uID) break;
            }
            if(m < h_n) continue;

            uID = link->a.a[h[j]].e.a[k].uID;
            ///o_d = get_phase_status(hap, link->a.a[h[j]].e.a[k].uID);
            o_d = g_p->a[g_p->index[uID]>>1].status[g_p->index[uID]&1];
            ///if(o_d < -1) fprintf(stderr, "ERROR\n");
            weight += (o_d*link->a.a[h[j]].e.a[k].weight);
        }
    }

    return weight;
}

void flip_unitig(G_partition* g_p, hc_links* link, bubble_type* bub, uint32_t id)
{
    if(g_p->a[id].h[0] > 0 && g_p->a[id].status[0] != 1 && g_p->a[id].status[0] != -1) return;
    if(g_p->a[id].h[1] > 0 && g_p->a[id].status[1] != 1 && g_p->a[id].status[1] != -1) return;
    uint32_t k, j, m, *h0, h0_n, *h1, h1_n, uID, *h = NULL, h_n;
    int status;
    double weight;
    get_phased_block(g_p, NULL, id, NULL, NULL, &h0, &h0_n, &h1, &h1_n, NULL, NULL);
    // fprintf(stderr, "h0_n: %u, h1_n: %u\n", h0_n, h1_n);
    if(h0_n > 0)
    {

        status = g_p->a[id].status[0];
        // fprintf(stderr, "+status: %d\n", status);
        h = h0; h_n = h0_n;
        for (j = 0; j < h_n; j++)
        {
            // fprintf(stderr, "+j: %u, h_n: %u\n", j, h_n);
            for (k = 0; k < link->a.a[h[j]].e.n; k++)
            {
                // fprintf(stderr, "+k: %u, e_n: %u\n", k, (uint32_t)link->a.a[h[j]].e.n);
                if(link->a.a[h[j]].e.a[k].del) continue;
                for (m = 0; m < h_n; m++)
                {
                    if(h[m] == link->a.a[h[j]].e.a[k].uID) break;
                }
                // fprintf(stderr, "+m: %u, h_n: %u\n", m, h_n);
                if(m < h_n) continue;

                uID = link->a.a[h[j]].e.a[k].uID;

                weight = link->a.a[h[j]].e.a[k].weight;

                if(g_p->index[uID] == (uint32_t)-1) continue;
                

                g_p->a[g_p->index[uID]>>1].weight[g_p->index[uID]&1] -= (2*status*weight);

                
            }
        }
        g_p->a[id].status[0] *= -1;
    }

    // fprintf(stderr, "hehehe\n");
    if(h1_n > 0)
    {
        status = g_p->a[id].status[1];
        // fprintf(stderr, "-status: %d\n", status);
        h = h1; h_n = h1_n;
        for (j = 0; j < h_n; j++)
        {
            // fprintf(stderr, "-j: %u, h_n: %u\n", j, h_n);
            for (k = 0; k < link->a.a[h[j]].e.n; k++)
            {
                // fprintf(stderr, "-k: %u, e_n: %u\n", k, (uint32_t)link->a.a[h[j]].e.n);
                if(link->a.a[h[j]].e.a[k].del) continue;
                for (m = 0; m < h_n; m++)
                {
                    if(h[m] == link->a.a[h[j]].e.a[k].uID) break;
                }
                // fprintf(stderr, "-m: %u, h_n: %u\n", m, h_n);
                if(m < h_n) continue;

                uID = link->a.a[h[j]].e.a[k].uID;
                // fprintf(stderr, "-uID: %u\n", uID);

                weight = link->a.a[h[j]].e.a[k].weight;

                if(g_p->index[uID] == (uint32_t)-1) continue;


                g_p->a[g_p->index[uID]>>1].weight[g_p->index[uID]&1] -= (2*status*weight);

            }
        }
        g_p->a[id].status[1] *= -1;
    }
}

void flip_unitig_debug(G_partition* g_p, hc_links* link, bubble_type* bub, uint32_t id)
{
    if(g_p->a[id].h[0] > 0 && g_p->a[id].status[0] != 1 && g_p->a[id].status[0] != -1) return;
    if(g_p->a[id].h[1] > 0 && g_p->a[id].status[1] != 1 && g_p->a[id].status[1] != -1) return;
    uint32_t k, j, m, *h0, h0_n, *h1, h1_n, uID, *h = NULL, h_n;
    int status;
    double weight;
    get_phased_block(g_p, NULL, id, NULL, NULL, &h0, &h0_n, &h1, &h1_n, NULL, NULL);
    // fprintf(stderr, "h0_n: %u, h1_n: %u\n", h0_n, h1_n);
    if(h0_n > 0)
    {

        status = g_p->a[id].status[0];
        // fprintf(stderr, "+status: %d\n", status);
        h = h0; h_n = h0_n;
        for (j = 0; j < h_n; j++)
        {
            // fprintf(stderr, "+j: %u, h_n: %u\n", j, h_n);
            for (k = 0; k < link->a.a[h[j]].e.n; k++)
            {
                // fprintf(stderr, "+k: %u, e_n: %u\n", k, (uint32_t)link->a.a[h[j]].e.n);
                if(link->a.a[h[j]].e.a[k].del) continue;
                for (m = 0; m < h_n; m++)
                {
                    if(h[m] == link->a.a[h[j]].e.a[k].uID) break;
                }
                // fprintf(stderr, "+m: %u, h_n: %u\n", m, h_n);
                if(m < h_n) continue;

                uID = link->a.a[h[j]].e.a[k].uID;

                weight = link->a.a[h[j]].e.a[k].weight;

                if(g_p->index[uID] == (uint32_t)-1) continue;

                if(id == 19675)
                {
                    fprintf(stderr, "+uID+: %u, e-weight: %f, g_p->index[uID]>>1: %u, status[0]: %d, pre_uID_weight: %f\n", 
                                             uID, weight, g_p->index[uID]>>1, status, g_p->a[g_p->index[uID]>>1].weight[g_p->index[uID]&1]);
                }
                

                g_p->a[g_p->index[uID]>>1].weight[g_p->index[uID]&1] -= (2*status*weight);

                if(id == 19675)
                {
                    fprintf(stderr, "+uID+: %u, new_uID_weight: %f\n", uID, 
                                     g_p->a[g_p->index[uID]>>1].weight[g_p->index[uID]&1]);
                }   
                
            }
        }
        g_p->a[id].status[0] *= -1;
    }

    // fprintf(stderr, "hehehe\n");
    if(h1_n > 0)
    {
        status = g_p->a[id].status[1];
        // fprintf(stderr, "-status: %d\n", status);
        h = h1; h_n = h1_n;
        for (j = 0; j < h_n; j++)
        {
            // fprintf(stderr, "-j: %u, h_n: %u\n", j, h_n);
            for (k = 0; k < link->a.a[h[j]].e.n; k++)
            {
                // fprintf(stderr, "-k: %u, e_n: %u\n", k, (uint32_t)link->a.a[h[j]].e.n);
                if(link->a.a[h[j]].e.a[k].del) continue;
                for (m = 0; m < h_n; m++)
                {
                    if(h[m] == link->a.a[h[j]].e.a[k].uID) break;
                }
                // fprintf(stderr, "-m: %u, h_n: %u\n", m, h_n);
                if(m < h_n) continue;

                uID = link->a.a[h[j]].e.a[k].uID;
                // fprintf(stderr, "-uID: %u\n", uID);

                weight = link->a.a[h[j]].e.a[k].weight;

                if(g_p->index[uID] == (uint32_t)-1) continue;

                if(id == 19675)
                {
                    fprintf(stderr, "-uID-: %u, e-weight: %f, g_p->index[uID]>>1: %u, status[1]: %d, pre_uID_weight: %f\n", 
                                             uID, weight, g_p->index[uID]>>1, status, g_p->a[g_p->index[uID]>>1].weight[g_p->index[uID]&1]);
                }

                g_p->a[g_p->index[uID]>>1].weight[g_p->index[uID]&1] -= (2*status*weight);

                if(id == 19675)
                {
                    fprintf(stderr, "-uID-: %u, new_uID_weight: %f\n", uID, 
                                     g_p->a[g_p->index[uID]>>1].weight[g_p->index[uID]&1]);
                }   
            }
        }
        g_p->a[id].status[1] *= -1;
    }
}


uint32_t phasing_improvement(H_partition* h, G_partition* g_p, ha_ug_index* idx, bubble_type* bub, hc_links* link)
{   
    uint32_t i, occ = 0, round = 0;
    double pre_w, pre_total, current_w;
    mul_block_phase_type b_x;
    init_mul_block_phase_type(&b_x, g_p, bub, asm_opt.thread_num, h);
    ///double index_time = yak_realtime();

    while(1)
    {
        pre_w = get_total_weight(h, g_p);
        pre_total = pre_w;

        while (1)
        {
            memset(h->lock, 0, sizeof(uint8_t)*g_p->n);
            while (1)
            {
                i = get_max_unitig(h, g_p, link, bub);
                if(i == (uint32_t)-1) break;
                h->lock[i] = 1;
                flip_unitig(g_p, link, bub, i);
                occ++;
            }
            current_w = get_total_weight(h, g_p);
            ///fprintf(stderr, "[M::%s::round single %u, pre_w: %f, current_w: %f]\n", __func__, round, pre_w, current_w);
            if(ceil(current_w) <= ceil(pre_w)) break;
            round++;
            pre_w = current_w;
        }

        
        pre_w = get_total_weight(h, g_p);
        round = 0;
        while (1)
        {
            ///fprintf(stderr, "[M::%s::round block %u, h->n: %lu]\n", __func__, round, h->n);
            phasing_improvement_by_block(h, g_p, bub, &b_x);
            current_w = get_total_weight(h, g_p);
            ///fprintf(stderr, "[M::%s::round block %u, pre_w: %f, current_w: %f]\n", __func__, round, pre_w, current_w);
            if(ceil(current_w) <= ceil(pre_w)) break;
            round++;
            pre_w = current_w;
        }
        
        if(ceil(current_w) <= ceil(pre_total)) break;
    }

    destory_mul_block_phase_type(&b_x);
    ///fprintf(stderr, "[M::%s:Flipping time:%.3f]\n", __func__, yak_realtime()-index_time);


    for (i = 0; i < g_p->n; i++)
    {
        update_partition_flag(h, g_p, link, i);
    }

    ///print_phase_group(g_p, bub, "Small");
    // double w0 = 0, w1 = 0;
    // uint32_t *h0, h0_n, *h1, h1_n;
    // get_phased_block(g_p, NULL, 2973, NULL, NULL, &h0, &h0_n, &h1, &h1_n, NULL, NULL);
    // w0 = get_cluster_weight_debug(g_p, h->link, h0, h0_n); 
    // w1 = get_cluster_weight_debug(g_p, h->link, h1, h1_n); 
    // fprintf(stderr, "debug-w0: %f, g_p->a[2973].weight[0]: %f\n", w0, g_p->a[2973].weight[0]);
    // fprintf(stderr, "debug-w1: %f, g_p->a[2973].weight[1]: %f\n", w1, g_p->a[2973].weight[1]);

    return !!occ;
} 

void destory_contig_partition(H_partition* hap)
{
    free(hap->lock);
    free(hap->hap);
    destory_G_partition(&(hap->g_p));
    destory_G_partition(&(hap->group_g_p));
    kv_destroy(hap->label_buffer);
    kv_destroy(hap->b.vis);
}

void label_unitigs(G_partition* g_p, ma_ug_t* ug)
{
    memset(R_INF.trio_flag, AMBIGU, R_INF.total_reads * sizeof(uint8_t));
    uint32_t i, k, j, *h0, h0_n, *h1, h1_n, uID, *h = NULL, h_n, flag = AMBIGU;
    int status;
    ma_utg_t *u = NULL;

    for (i = 0; i < g_p->n; i++)
    {
        if(g_p->a[i].h[0] > 0 && g_p->a[i].status[0] != 1 && g_p->a[i].status[0] != -1) continue;
        if(g_p->a[i].h[1] > 0 && g_p->a[i].status[1] != 1 && g_p->a[i].status[1] != -1) continue;
        get_phased_block(g_p, NULL, i, NULL, NULL, &h0, &h0_n, &h1, &h1_n, NULL, NULL);

        status = g_p->a[i].status[0];
        h = h0; h_n = h0_n;
        if(status == 1)
        {
            flag = FATHER;
        }
        else if (status == -1)
        {
            flag = MOTHER;
        }
        for (j = 0; j < h_n; j++)
        {
            uID = h[j];
            u = &ug->u.a[uID];
            if(u->m == 0) continue;
            for (k = 0; k < u->n; k++)
            {
                R_INF.trio_flag[u->a[k]>>33] = flag;
            }
        }



        status = g_p->a[i].status[1];
        h = h1; h_n = h1_n;
        if(status == 1)
        {
            flag = FATHER;
        }
        else if (status == -1)
        {
            flag = MOTHER;
        }
        for (j = 0; j < h_n; j++)
        {
            uID = h[j];
            u = &ug->u.a[uID];
            if(u->m == 0) continue;
            for (k = 0; k < u->n; k++)
            {
                R_INF.trio_flag[u->a[k]>>33] = flag;
            }
        }
    }


    uint64_t occ = 0; 
    for (i = 0; i < ug->u.n; i++)
    {
        occ += ug->u.a[i].n;
    }

    ///fprintf(stderr, "# reads: %lu\n", occ);

    for (i = occ = 0; i < R_INF.total_reads; i++)
    {
        if(R_INF.trio_flag[i] == FATHER) occ++;
    }

    ///fprintf(stderr, "# Father reads: %lu\n", occ);

    for (i = occ = 0; i < R_INF.total_reads; i++)
    {
        if(R_INF.trio_flag[i] == MOTHER) occ++;
    }
    
    ///fprintf(stderr, "# Mother reads: %lu\n", occ);
}

void label_unitigs_sm(int8_t *s, mc_gg_status *sa, ma_ug_t* ug)
{
    memset(R_INF.trio_flag, AMBIGU, R_INF.total_reads * sizeof(uint8_t));
    uint32_t i, k, flag = AMBIGU;
    ma_utg_t *u = NULL;

    for (i = 0; i < ug->g->n_seq; i++)
    {
        if(ug->g->seq[i].del) continue;
        flag = 0;
        if(s)
        {
            if(s[i] == 0) continue;
            flag = (s[i] > 0? FATHER:MOTHER);
        }

        if(sa)
        {
            if(sa[i].s != 1 && sa[i].s != 2) continue;
            flag = sa[i].s;
        }
        
        u = &ug->u.a[i];
        if(u->m == 0) continue;
        for (k = 0; k < u->n; k++)
        {
            R_INF.trio_flag[u->a[k]>>33] = flag;
        }
    }
}


void print_bubble_graph(bubble_type* bub, ma_ug_t* ug, const char* prefix, FILE *fp)
{
    uint32_t i, k, *a, n, beg, sink, x;
    asg_t *b_g = bub->b_g;
    char name[32];
    for (i = 0; i < b_g->n_seq; i++)
    {
        get_bubbles(bub, i, &beg, &sink, &a, &n, NULL);
        sprintf(name, "%s%.6d%c", prefix, i, "fb"[i<bub->f_bub?0:1]);
        fprintf(stderr, "S\t%s\t*\tLN:i:%d\n", name, n);

        if(beg != (uint32_t)-1) fprintf(stderr, "A\tutg%.6d%c\t%s\n", (beg>>1)+1, "lc"[ug->u.a[(beg>>1)].circ], "beg");
        if(sink != (uint32_t)-1) fprintf(stderr, "A\tutg%.6d%c\t%s\n", (sink>>1)+1, "lc"[ug->u.a[(sink>>1)].circ], "sink");
        for (k = 0; k < n; k++)
        {
            x = a[k]>>1;
            fprintf(stderr, "A\tutg%.6d%c\t%s\n", x+1, "lc"[ug->u.a[x].circ], "mid");
        }
    }

    asg_arc_t* au = NULL;
    uint32_t nu, u, v;
    for (i = 0; i < b_g->n_seq; i++)
    {
        u = i<<1;
        au = asg_arc_a(b_g, u);
        nu = asg_arc_n(b_g, u);
        for (k = 0; k < nu; k++)
        {
            if(au[k].del) continue;
            v = au[k].v;
            fprintf(stderr, "L\t%s%.6d%c\t%c\t%s%.6d%c\t%c\t%dM\tL1:i:%d\n", 
            prefix, u>>1, "fb"[(u>>1)<bub->f_bub?0:1], "+-"[u&1],
            prefix, v>>1, "fb"[(v>>1)<bub->f_bub?0:1], "+-"[v&1], 0, 0);


            asg_arc_t* av = asg_arc_a(b_g, v^1);
            uint32_t nv = asg_arc_n(b_g, v^1), m;
            for (m = 0; m < nv; m++)
            {
                if(av[m].del) continue;
                if(av[m].v == (u^1)) break;
            }

            if(m == nv) fprintf(stderr, "sb1sb, nv: %u, nu: %u\n", nv, nu);
        }


        u = (i<<1) + 1;
        au = asg_arc_a(ug->g, u);
        nu = asg_arc_n(ug->g, u);
        for (k = 0; k < nu; k++)
        {
            if(au[k].del) continue;
            v = au[k].v;
            fprintf(stderr, "L\t%s%.6d%c\t%c\t%s%.6d%c\t%c\t%dM\tL1:i:%d\n", 
            prefix, u>>1, "fb"[(u>>1)<bub->f_bub?0:1], "+-"[u&1],
            prefix, v>>1, "fb"[(v>>1)<bub->f_bub?0:1], "+-"[v&1], 0, 0);

            asg_arc_t* av = asg_arc_a(b_g, v^1);
            uint32_t nv = asg_arc_n(b_g, v^1), m;
            for (m = 0; m < nv; m++)
            {
                if(av[m].del) continue;
                if(av[m].v == (u^1)) break;
            }

            if(m == nv) fprintf(stderr, "sb2sb, nv: %u, nu: %u\n", nv, nu);
        }
    }
}



void print_bubble_utg(bubble_type* bub, ma_ug_t* unitig_ug, const char* prefix, FILE *fp)
{
    uint32_t i, k, *a, n, beg, sink, x, occ;
    ma_ug_t *b_ug = bub->b_ug;
    char name[32];
    for (i = 0; i < b_ug->u.n; i++)
    {
        ma_utg_t *p = &b_ug->u.a[i];
        if(p->n == 0) continue;
        for (k = occ = 0; k < p->n; k++)
        {
            x = p->a[k]>>33;
            get_bubbles(bub, x, &beg, &sink, &a, &n, NULL);
            occ += n;
        }
        sprintf(name, "%s%.6d%c", prefix, i + 1, "lc"[p->circ]);
        fprintf(fp, "S\t%s\t*\tLN:i:%u\n", name, occ);
        for (k = 0; k < p->n; k++)
        {
            x = p->a[k]>>33;
            get_bubbles(bub, x, &beg, &sink, &a, &n, NULL);
            if(beg != (uint32_t)-1) fprintf(fp, "A\tutg%.6d%c\t%u\t%s\n", (beg>>1)+1, "lc"[unitig_ug->u.a[(beg>>1)].circ], n, "beg");
            if(sink != (uint32_t)-1) fprintf(fp, "A\tutg%.6d%c\t%u\t%s\n", (sink>>1)+1, "lc"[unitig_ug->u.a[(sink>>1)].circ], n, "sink");
        }
    }

    asg_arc_t* au = NULL;
    uint32_t nu, u, v, j;
    for (i = 0; i < b_ug->u.n; ++i) {
        if(b_ug->u.a[i].m == 0) continue;
        if(b_ug->u.a[i].circ)
        {
            fprintf(fp, "L\t%s%.6dc\t+\t%s%.6dc\t+\t%dM\tL1:i:%d\n", 
            prefix, i+1, prefix, i+1, 0, 0);
            fprintf(fp, "L\t%s%.6dc\t-\t%s%.6dc\t-\t%dM\tL1:i:%d\n", 
            prefix, i+1, prefix, i+1, 0, 0);
        } 
        u = i<<1;
        au = asg_arc_a(b_ug->g, u);
        nu = asg_arc_n(b_ug->g, u);
        for (j = 0; j < nu; j++)
        {
            if(au[j].del) continue;
            v = au[j].v;
            fprintf(fp, "L\t%s%.6d%c\t%c\t%s%.6d%c\t%c\t%dM\tL1:i:%d\n", 
            prefix, (u>>1)+1, "lc"[b_ug->u.a[u>>1].circ], "+-"[u&1],
            prefix, (v>>1)+1, "lc"[b_ug->u.a[v>>1].circ], "+-"[v&1], 0, 0);
        }


        u = (i<<1) + 1;
        au = asg_arc_a(b_ug->g, u);
        nu = asg_arc_n(b_ug->g, u);
        for (j = 0; j < nu; j++)
        {
            if(au[j].del) continue;
            v = au[j].v;
            fprintf(fp, "L\t%s%.6d%c\t%c\t%s%.6d%c\t%c\t%dM\tL1:i:%d\n", 
            prefix, (u>>1)+1, "lc"[b_ug->u.a[u>>1].circ], "+-"[u&1],
            prefix, (v>>1)+1, "lc"[b_ug->u.a[v>>1].circ], "+-"[v&1], 0, 0);
        }
    }

    
}


void print_debug_bubble_graph(bubble_type* bub, ma_ug_t* ug, const char *fn)
{
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
    sprintf(buf, "%s.bub.gfa", fn);
    FILE* fp = fopen(buf, "w");

    print_bubble_utg(bub, ug, "btg", fp);

    fclose(fp);
    free(buf);
}

void print_bubble_chain(bubble_type* bub)
{
    uint32_t m, i, k;
    uint32_t beg, sink, *a = NULL, n;
    uint64_t bid;
    ma_utg_t *u = NULL; 
    for (m = 0; m < bub->chain_weight.n; m++)
    {
        if(bub->chain_weight.a[m].del) continue;
        u = &(bub->b_ug->u.a[bub->chain_weight.a[m].id]);
        fprintf(stderr, "\nChain_id=%lu\n", bub->chain_weight.a[m].id);
        for (i = 0; i < u->n; i++)
        {
            bid = u->a[i]>>33;
            get_bubbles(bub, bid, &beg, &sink, &a, &n, NULL);
            fprintf(stderr, "btg%.6lu%c, beg-utg%.6ul, sink-utg%.6ul\n", 
                        bid, "fb"[bid<bub->f_bub?0:1], (beg>>1)+1, (sink>>1)+1);
            for (k = 0; k < n; k++)
            {
                if(k != 0 && (k%5)==0) fprintf(stderr, "\n");
                fprintf(stderr, "m-utg%.6ul\t", (a[k]>>1)+1);
            }
            fprintf(stderr, "\n");
            
        }
    }
}

void init_contig_H_partition(bubble_type* bub, ha_ug_index* idx, H_partition* hap)
{
    uint32_t i, k_i, k_j, uID, *a = NULL, n, *h0, h0_n, *h1, h1_n;
    destory_G_partition(&(hap->group_g_p)); memset(&(hap->group_g_p), 0, sizeof(G_partition));
    init_G_partition(&(hap->group_g_p), hap->n);
    partition_warp *res = NULL;
    ma_utg_t *u_x = NULL;
    chain_hic_warp *c_w = &(bub->c_w);

    for (i = 0; i < bub->c_w.n; i++)
    {
        kv_pushp(partition_warp, hap->group_g_p, &res);
        kv_init(res->a);
        res->full_bub = 0;
        res->h[0] = res->h[1] = 0;

        res->status[0] = 1; res->status[1] = -1;
        res->weight[0] = res->weight[1] = res->weight_convex = 0;


        u_x = (*c_w).a[(*c_w).a[i].id].u;
        for (k_i = 0; k_i < u_x->n; k_i++)
        {
            get_bubbles(bub, u_x->a[k_i]>>33, NULL, NULL, &a, &n, NULL);
            for (k_j = 0; k_j < n; k_j++)
            {
                uID = a[k_j]>>1;
                if((*c_w).chain_idx[uID] != (*c_w).a[i].id) continue;
                if(get_phase_status(hap, uID)==1)
                {
                    kv_push(uint32_t, res->a, uID);
                    res->h[0]++;
                }
            }
        }

        for (k_i = 0; k_i < u_x->n; k_i++)
        {
            get_bubbles(bub, u_x->a[k_i]>>33, NULL, NULL, &a, &n, NULL);
            for (k_j = 0; k_j < n; k_j++)
            {
                uID = a[k_j]>>1;
                if((*c_w).chain_idx[uID] != (*c_w).a[i].id) continue;
                if(get_phase_status(hap, uID)==-1)
                {
                    kv_push(uint32_t, res->a, uID);
                    res->h[1]++;
                }
            }
        }


        for (k_i = 0; k_i < res->h[0]; k_i++)
        {
            ///if(hap->group_g_p.index[res->a.a[k_i]] != (uint32_t)-1) fprintf(stderr, "ERROR---00\n");
            hap->group_g_p.index[res->a.a[k_i]] = hap->group_g_p.n-1;
            hap->group_g_p.index[res->a.a[k_i]] = hap->group_g_p.index[res->a.a[k_i]] << 1;
        } 

        for (; k_i < res->a.n; k_i++)
        {
            ///if(hap->group_g_p.index[res->a.a[k_i]] != (uint32_t)-1) fprintf(stderr, "ERROR---11\n");
            hap->group_g_p.index[res->a.a[k_i]] = hap->group_g_p.n-1;
            hap->group_g_p.index[res->a.a[k_i]] = (hap->group_g_p.index[res->a.a[k_i]] << 1) + 1;
        }

        get_phased_block(&(hap->group_g_p), NULL, i, NULL, NULL, &h0, &h0_n, &h1, &h1_n, NULL, NULL);
        if(h0_n >0) res->weight[0] = get_cluster_weight(hap, hap->link, h0, h0_n); 
        if(h1_n >0) res->weight[1] = get_cluster_weight(hap, hap->link, h1, h1_n); 
        res->weight_convex = get_cluster_inner_weight(hap, hap->link, h0, h0_n, h1, h1_n);
    }

    flip_by_node(hap, &(hap->group_g_p), bub);
    label_unitigs(&(hap->group_g_p), idx->ug);
}

void reset_H_partition(H_partition* hap, uint32_t is_init)
{
    if(!is_init)
    {
        hap->n = 0;
        free(hap->lock);
        free(hap->hap);
        hap->m[0] = hap->m[1] = hap->m[2] = (uint32_t)-1;
        hap->label = hap->label_add = hap->label_shift = (uint32_t)-1;
        destory_G_partition(&(hap->g_p)); memset(&(hap->g_p), 0, sizeof(G_partition));
        destory_G_partition(&(hap->group_g_p)); memset(&(hap->group_g_p), 0, sizeof(G_partition));
        kv_destroy(hap->label_buffer); kv_init(hap->label_buffer);
        kv_destroy(hap->b.vis); kv_init(hap->b.vis); memset(&(hap->b), 0, sizeof(block_phase_type));
    }

    memset(hap, 0, sizeof(H_partition));
}


int alignment_worker_pipeline(sldat_t* sl, const enzyme *fn1, const enzyme *fn2)
{
    double index_time = yak_realtime();
    int i;
    for (i = 0; i < fn1->n && i < fn2->n; i++)
    {
        gzFile fp1, fp2;
        if ((fp1 = gzopen(fn1->a[i], "r")) == 0) return 0;
        if ((fp2 = gzopen(fn2->a[i], "r")) == 0) return 0;
        sl->ks1 = kseq_init(fp1);
        sl->ks2 = kseq_init(fp2);

        kt_pipeline(3, worker_pipeline, sl, 3);
        
        kseq_destroy(sl->ks1);
        kseq_destroy(sl->ks2);
        gzclose(fp1);
        gzclose(fp2);
    }
    fprintf(stderr, "[M::%s::%.3f] ==> Qualification\n", __func__, yak_realtime()-index_time);

    dedup_hits(&(sl->hits), 1);    
    return 1;
}

void debug_gfa_space(ha_ug_index* idx, ma_ug_t* ug, trans_chain* t_ch, kv_u_trans_t *ref)
{
    bubble_type bub; 
    memset(&bub, 0, sizeof(bubble_type));
    bub.round_id = 0; bub.n_round = 2;

    identify_bubbles(ug, &bub, t_ch->ir_het, ref);

    hc_links link;
    init_hc_links(&link, ug->g->n_seq, t_ch);

    measure_distance(idx, ug, NULL, &link, &bub, &(t_ch->k_trans));

    // uint32_t i, k;
    // for (i = 0; i < link.a.n; ++i) 
    // {
    //     for (k = 0; k < link.a.a[i].e.n; k++)
    //     {
    //         if(link.a.a[i].e.a[k].del || link.a.a[i].e.a[k].dis == (uint64_t)-1) continue;
    //         fprintf(stderr, "s-utg%.6dl\td-utg%.6dl\t%lu\n", 
    //             (int)(i+1), (int)(link.a.a[i].e.a[k].uID+1), 
    //             link.a.a[i].e.a[k].dis == (uint64_t)-1? (uint64_t)-1 : link.a.a[i].e.a[k].dis>>3);
    //     }
    // }



    destory_bubbles(&bub);
    destory_hc_links(&link);
}

void idx_hc_links(kvec_pe_hit* hits, ha_ug_index* idx, bubble_type* bub)
{
    uint64_t k, l;
    uint32_t qn, tn;

    kv_resize(uint64_t, hits->idx, idx->ug->g->n_seq);
    hits->idx.n = idx->ug->g->n_seq;
    memset(hits->idx.a, 0, hits->idx.n*sizeof(uint64_t));

    kv_resize(uint64_t, hits->occ, idx->ug->g->n_seq);
    hits->occ.n = idx->ug->g->n_seq;
    memset(hits->occ.a, 0, hits->occ.n*sizeof(uint64_t));

    radix_sort_pe_hit_idx_an1(hits->a.a, hits->a.a + hits->a.n);
    for (k = 1, l = 0; k <= hits->a.n; ++k) 
    {   
        if (k == hits->a.n || 
            ((hits->a.a[k].s<<1)>>(64 - idx->uID_bits)) != ((hits->a.a[l].s<<1)>>(64 - idx->uID_bits))) 
        {
            if (k - l > 1) radix_sort_pe_hit_idx_an2(hits->a.a + l, hits->a.a + k);
            
            hits->idx.a[((hits->a.a[l].s<<1)>>(64 - idx->uID_bits))] 
                                                = (uint64_t)l << 32 | (k - l);

            for (; l < k; l++)
            {
                qn = ((hits->a.a[l].s<<1)>>(64 - idx->uID_bits));
                tn = ((hits->a.a[l].e<<1)>>(64 - idx->uID_bits));
                if(bub && IF_HOM(qn, *bub)) continue;
                if(bub && IF_HOM(tn, *bub)) continue;
                hits->occ.a[qn]++;
                hits->occ.a[tn]++;
            }
            
            l = k;
        }
    }
}

inline uint32_t trans_checking_pass(bubble_type* bub, kv_u_trans_t *ref, uint32_t x, uint32_t y)
{
    if(u_trans_n(*ref, x) == 0 || u_trans_n(*ref, y) == 0) return 0;
    u_trans_t *a = NULL;
    uint32_t n, i, f[2], qn, tn;

    qn = x; tn = y;
    a = u_trans_a(*ref, qn); n = u_trans_n(*ref, qn);
    for (i = 0, f[0] = f[1] = 0; i < n; i++)
    {
        if(a[i].del) continue;
        if(a[i].f == RC_2) continue;
        if(IF_HOM(a[i].tn, *bub)) continue;
        f[(a[i].tn == tn && a[i].f != RC_2)]++;
    }
    if(f[0] != 0) return 0;
    if(f[1] == 0) return 0;

    qn = y; tn = x;
    a = u_trans_a(*ref, qn); n = u_trans_n(*ref, qn);
    for (i = 0, f[0] = f[1] = 0; i < n; i++)
    {
        if(a[i].del) continue;
        if(IF_HOM(a[i].tn, *bub)) continue;
        f[(a[i].tn == tn && a[i].f != RC_2)]++;
    }
    if(f[0] != 0) return 0;
    if(f[1] == 0) return 0;

    fprintf(stderr, "M::%s::s-utg%.6ul<----->d-utg%.6ul\n", __func__, x+1, y+1);
    return 1;
}

uint32_t get_u_trans_spec_idx(kv_u_trans_t *ta, uint32_t qn, uint32_t tn, u_trans_t **r_a, uint32_t *occ, uint32_t *idx)
{
    if(r_a) (*r_a) = NULL; 
    if(occ) (*occ) = 0;
    if(idx) (*idx) = (uint32_t)-1;
    u_trans_t *a = NULL;
    uint32_t n, st, i;
    a = u_trans_a(*ta, qn);
    n = u_trans_n(*ta, qn);
    for (st = 0, i = 1; i <= n; ++i)
    {
        if (i == n || a[i].tn != a[st].tn)
        {
            if(a[st].tn == tn)
            {
                if(r_a) (*r_a) = a + st;
                if(occ) (*occ) = i - st;
                if(idx) (*idx) = st + ((*ta).idx.a[(qn)]>>32);
                return 1;
            }
            st = i;
        }
    }
    return 0;
}

void interpr_hit(ha_ug_index* idx, uint64_t x, uint32_t rLen, uint32_t *uid, uint32_t *beg, uint32_t *end);
int hic_sc_type(ha_ug_index* idx, kvec_pe_hit* hits, uint64_t k)
{
    uint32_t s_uid, s_beg, s_end, e_uid, e_beg, e_end, slen, elen, x = 0;

    interpr_hit(idx, hits->a.a[k].s, hits->a.a[k].len>>32, &s_uid, &s_beg, &s_end); 
    s_beg = (s_beg+s_end)>>1; slen = idx->ug->u.a[s_uid].len;
    if(s_beg >= (slen>>1)) x+=1;

    interpr_hit(idx, hits->a.a[k].e, (uint32_t)hits->a.a[k].len, &e_uid, &e_beg, &e_end); 
    e_beg = (e_beg+e_end)>>1; elen = idx->ug->u.a[e_uid].len;
    if(e_beg >= (elen>>1)) x+=2;
    return x;
}

void weight_kv_u_trans(ha_ug_index* idx, kvec_pe_hit* hits, hc_links* link, bubble_type* bub, 
kv_u_trans_t *ta, trans_idx* dis, int sc_weight)
{
    uint64_t k, i, shif = 64 - idx->uID_bits, beg, end, t_d;
    u_trans_t *e1 = NULL, *e2 = NULL;
    long double weight;
    double *sw = NULL;
    u_trans_t *p = NULL;
    uint32_t is_cc, ii1, ii2;

    for (i = 0, ta->idx.n = ta->n = 0; i < link->a.n; i++)
    {
        for (k = 0; k < link->a.a[i].e.n; k++)
        {
            if(link->a.a[i].e.a[k].del) continue;
            if(IF_HOM(i, *bub)) continue;
            if(IF_HOM(link->a.a[i].e.a[k].uID, *bub)) continue;
            if(i == link->a.a[i].e.a[k].uID) continue;

            kv_pushp(u_trans_t, *ta, &p);
            memset(p, 0, sizeof(u_trans_t));
            p->qn = i; p->tn = link->a.a[i].e.a[k].uID;
            p->nw = 0; p->occ = 0;
        }
    }
    kt_u_trans_t_idx(ta, idx->ug->g->n_seq);
    if(sc_weight) {
        k = ta->n*3;
        MALLOC(sw, k);
        for (i = 0; i < k; i++) sw[i] = 0;
    }


    for (k = 0; k < hits->a.n; ++k) 
    {
        beg = ((hits->a.a[k].s<<1)>>shif);
        end = ((hits->a.a[k].e<<1)>>shif);

        if(beg == end) continue;
        if(IF_HOM(beg, *bub)) continue;
        if(IF_HOM(end, *bub)) continue;
        if(is_hom_hit(hits->a.a[k])) continue;
        
        t_d = get_hic_distance(&(hits->a.a[k]), link, idx, &is_cc);
        // if(t_d == (uint64_t)-1) continue;
        // if(t_d == (uint64_t)-1 && is_cc == 0) continue;
        // if(t_d == (uint64_t)-1 && !dis) continue;

        get_u_trans_spec_idx(ta, beg, end, &e1, NULL, &ii1);
        get_u_trans_spec_idx(ta, end, beg, &e2, NULL, &ii2);

        if(e1 == NULL || e2 == NULL) continue;
        weight = 1;
        if(dis) weight = get_trans_weight_advance(idx, t_d, dis);
        
        if(sc_weight){
            i = hic_sc_type(idx, hits, k);
            if(i == 0){
                e1->nw -= weight; e2->nw -= weight;
            }
            else{
                i--;
                sw[(ii1*3)+i] -= weight;
                sw[(ii2*3)+i] -= weight;
            }
        } else{
            e1->nw -= weight; e2->nw -= weight; 
        }
        e1->occ++; e2->occ++;
    }    

    if(sc_weight) {
        for (i = 0; i < ta->n; ++i){
            if(ta->a[i].nw > sw[(i*3)]) ta->a[i].nw = sw[(i*3)];
            if(ta->a[i].nw > sw[(i*3)+1]) ta->a[i].nw = sw[(i*3)+1];
            if(ta->a[i].nw > sw[(i*3)+2]) ta->a[i].nw = sw[(i*3)+2];
            ta->a[i].nw *= 2/**4**/;
        }
        free(sw);
    }
}

void interpr_hit(ha_ug_index* idx, uint64_t x, uint32_t rLen, uint32_t *uid, uint32_t *beg, uint32_t *end)
{
    if(uid) (*uid) = ((x<<1)>>(64 - idx->uID_bits)); 
    uint32_t rev = (x>>63);
    long long ref_p = x & idx->pos_mode;
    long long p_beg, p_end;

    if(rev)
    {
        p_end = ref_p;
        p_beg = p_end + 1 - rLen;
    }
    else
    {
        p_beg = ref_p; 
        p_end = p_beg + rLen - 1;
    }
    if(p_beg < 0) p_beg = 0;
    if(p_end < 0) p_end = 0;
    if(beg) (*beg) = p_beg; 
    if(end) (*end) = p_end + 1;
}
double get_interval_weight(ha_ug_index* idx, hc_links* link, trans_idx* dis, 
pe_hit *hits, uint32_t occ, uint32_t qid, uint32_t qs, uint32_t qe, uint32_t tid, uint32_t ts, uint32_t te)
{
    int64_t s_idx = 0, e_idx = (int64_t)occ - 1, m_idx = 0;
    uint32_t m_uid = (uint32_t)-1;
    while (s_idx <= e_idx)
    {
        m_idx = s_idx + (e_idx - s_idx)/2;
        m_uid = ((hits[m_idx].e<<1)>>(64 - idx->uID_bits)); 
        if (m_uid == tid)
            break;
        if (m_uid < tid)
            s_idx = m_idx + 1;
        else
            e_idx = m_idx - 1;
    }
    if(m_uid != tid) return 0;
    

    uint32_t k, s_uid, s_beg, s_end, e_uid, e_beg, e_end;
    uint64_t t_d;
    double w, weight;
    w = 0;
    for (k = m_idx; k < occ; k++)///all hits of qid
    {   
        interpr_hit(idx, hits[k].s, hits[k].len>>32, &s_uid, &s_beg, &s_end);
        if(s_uid != qid) continue;
        if(!(qs <= s_beg && qe >= s_end)) continue;

        interpr_hit(idx, hits[k].e, (uint32_t)hits[k].len, &e_uid, &e_beg, &e_end);
        if(e_uid != tid) break;
        if(!(ts <= e_beg && te >= e_end)) continue;

        t_d = get_hic_distance(&hits[k], link, idx, NULL);
        if(t_d == (uint64_t)-1) continue;

        weight = 1;
        if(dis) weight = get_trans_weight_advance(idx, t_d, dis);

        w += weight;
    }

    for (m_idx -= 1; m_idx >= 0; m_idx--)
    {
        k = m_idx;
        interpr_hit(idx, hits[k].s, hits[k].len>>32, &s_uid, &s_beg, &s_end);
        if(s_uid != qid) continue;
        if(!(qs <= s_beg && qe >= s_end)) continue;

        interpr_hit(idx, hits[k].e, (uint32_t)hits[k].len, &e_uid, &e_beg, &e_end);
        if(e_uid != tid) break;
        if(!(ts <= e_beg && te >= e_end)) continue;

        t_d = get_hic_distance(&hits[k], link, idx, NULL);
        if(t_d == (uint64_t)-1) continue;

        weight = 1;
        if(dis) weight = get_trans_weight_advance(idx, t_d, dis);

        w += weight;
    }
    return w;
}

double get_hits_weight(ha_ug_index* idx, bubble_type* bub, kvec_pe_hit* hits, hc_links* link, trans_idx* dis, 
u_trans_t *t_a, uint32_t t_n, uint32_t qid, kv_u_trans_t *ta_idx)
{
    /****************************may have bugs********************************/
    ///need to record self hits
    ///if(u_trans_n(*ta_idx, qid) == 0) return 0;///no hit bridging qid
    /****************************may have bugs********************************/
    uint32_t k, i, x_n, y;
    u_trans_t *x_a = NULL;
    double occ_q, occ_t;
    double w, i_w;
    for (k = 0, w = 0; k < t_n; k++)
    {
        if(IF_HOM(qid, *bub)) continue;
        if(IF_HOM(t_a[k].tn, *bub)) continue;
        /****************************may have bugs********************************/
        if(qid != t_a[k].tn)
        {
            x_n = u_trans_n(*ta_idx, qid); x_a = u_trans_a(*ta_idx, qid); y = t_a[k].tn;
            if(u_trans_n(*ta_idx, t_a[k].tn) < x_n)
            {
                x_n = u_trans_n(*ta_idx, t_a[k].tn); x_a = u_trans_a(*ta_idx, t_a[k].tn); y = qid;
            }
            if(x_n == 0) continue; 

            for (i = 0; i < x_n; i++)
            {
                if(x_a[i].tn == y) break;
            }
            if(i >= x_n) continue; ///no hit bridging tn and qid
        }
        /****************************may have bugs********************************/
        occ_q = hits->occ.a[qid];
        occ_t = hits->occ.a[t_a[k].tn] * ((double)(t_a[k].te - t_a[k].ts) / (double)(idx->ug->g->seq[t_a[k].tn].len));
        if(occ_t < 1) occ_t = 1;

        ///q--->t hits
        i_w = get_interval_weight(idx, link, dis, hits->a.a + (hits->idx.a[qid]>>32),
        (uint32_t)(hits->idx.a[qid]), qid, 0, idx->ug->g->seq[qid].len, t_a[k].tn, 
        t_a[k].ts, t_a[k].te);
        if(i_w != 0) i_w /= (MIN(occ_q, occ_t));
        w += i_w;
        if(t_a[k].tn == qid) continue;

        ///t--->q hits
        i_w = get_interval_weight(idx, link, dis, hits->a.a + (hits->idx.a[t_a[k].tn]>>32),
        (uint32_t)(hits->idx.a[t_a[k].tn]), t_a[k].tn, t_a[k].ts, t_a[k].te, 
        qid, 0, idx->ug->g->seq[qid].len);    
        if(i_w != 0) i_w /= (MIN(occ_q, occ_t));
        w += i_w;
    }
    return w;
}


void adjust_weight_kv_u_trans(ha_ug_index* idx, kvec_pe_hit* hits, hc_links* link, bubble_type* bub, 
kv_u_trans_t *ta, kv_u_trans_t *ref, trans_idx* dis)
{
    u_trans_t *a = NULL, *p = NULL;
    uint32_t n, k, i, m, qn, tn;
    uint8_t *vis = NULL; CALLOC(vis, idx->ug->g->n_seq);
    double w;

    for (i = m = 0; i < ta->n; i++)
    {
        if(ta->a[i].del) continue;
        if(IF_HOM(ta->a[i].qn, *bub)) continue;
        if(IF_HOM(ta->a[i].tn, *bub)) continue;
        ta->a[m] = ta->a[i];
        if(ta->a[m].nw != 0)
        {
            ta->a[m].nw /= (double)(MIN(hits->occ.a[ta->a[m].qn], hits->occ.a[ta->a[m].tn]));
        }
        /*******************************for debug************************************/
        // trans_checking_pass(bub, ref, ta->a[i].qn, ta->a[i].tn);
        // if(trans_checking_pass(bub, ref, ta->a[i].qn, ta->a[i].tn)) ta->a[m].nw = 0;
        /*******************************for debug************************************/
        
        m++;
    }
    ta->n = m;
    
    for (k = 0; k < ta->idx.n; k++)///all nodes
    {
        if(IF_HOM(k, *bub)) continue;
        a = u_trans_a(*ta, k);
        n = u_trans_n(*ta, k);
        ///for each pair qn, tn
        ///count hic pairs between (qn, tn^) and (qn^1, tn)
        for (i = 0, vis[k] = 1; i < n; i++)
        {
            vis[a[i].tn] = 1;
            if(a[i].qn == a[i].tn) continue;
            if(IF_HOM(a[i].qn, *bub)) continue;
            if(IF_HOM(a[i].tn, *bub)) continue;
            ///(qn, tn^)
            a[i].nw += get_hits_weight(idx, bub, hits, link, dis, 
                                    u_trans_a(*ref, a[i].tn), u_trans_n(*ref, a[i].tn), a[i].qn, ta);
            ///(qn^1, tn)
            a[i].nw += get_hits_weight(idx, bub, hits, link, dis, 
                                    u_trans_a(*ref, a[i].qn), u_trans_n(*ref, a[i].qn), a[i].tn, ta);
        }
        
        
        for (i = 0; i < ta->idx.n; i++)///all edges
        {
            qn = k; tn = i; w = 0;
            if(IF_HOM(qn, *bub)) continue;
            if(IF_HOM(tn, *bub)) continue;
            if(vis[tn]) continue;
            if(qn == tn) continue;
            ///(qn, tn^)
            w += get_hits_weight(idx, bub, hits, link, dis, 
                                u_trans_a(*ref, tn), u_trans_n(*ref, tn), qn, ta);
            ///(qn^1, tn)
            w += get_hits_weight(idx, bub, hits, link, dis, 
                                u_trans_a(*ref, qn), u_trans_n(*ref, qn), tn, ta);
            if(w == 0) continue;
            kv_pushp(u_trans_t, *ta, &p);
            memset(p, 0, sizeof(u_trans_t));///extra edges
            p->nw = w; p->qn = qn; p->tn = tn;
        }

        // for (i = 0; i < u_trans_n(*ref, k); i++)///only trans edges
        // {
        //     qn = k; tn = (u_trans_a(*ref, k))[i].tn; w = 0;
        //     if(IF_HOM(qn, *bub)) continue;
        //     if(IF_HOM(tn, *bub)) continue;
        //     if(vis[tn]) continue;
        //     if(qn == tn) continue;
        //     ///(qn, tn^)
        //     w += get_hits_weight(idx, hits, link, dis, 
        //                         u_trans_a(*ref, tn), u_trans_n(*ref, tn), qn, ta);
        //     ///(qn^1, tn)
        //     w += get_hits_weight(idx, hits, link, dis, 
        //                         u_trans_a(*ref, qn), u_trans_n(*ref, qn), tn, ta);
        //     if(w == 0) continue;
        //     kv_pushp(u_trans_t, *ta, &p);
        //     memset(p, 0, sizeof(u_trans_t));///extra edges
        //     p->nw = w; p->qn = qn; p->tn = tn;
        // }

        for (i = 0, vis[k] = 0; i < n; i++) vis[a[i].tn] = 0;
    }

    // fprintf(stderr, "------ta->n=%u\n", (uint32_t)ta->n);

    for (i = m = 0; i < ta->n; i++)
    {
        if(ta->a[i].nw == 0 || ta->a[i].del) continue;
        ta->a[m] = ta->a[i];
        m++;
    }
    ta->n = m;

    free(vis);
    kt_u_trans_t_idx(ta, idx->ug->g->n_seq);
}


inline uint32_t get_trans_interval_weight(ha_ug_index* idx, hc_links* link, bubble_type* bub, 
trans_idx* dis, pe_hit *hit, uint32_t hit_n, uint32_t qn, uint32_t qs, uint32_t qe, uint32_t tn, 
uint32_t ts, uint32_t te, double *w_a) 
{
    uint32_t i, s_uid, s_beg, s_end, e_uid, e_beg, e_end, found, is_cc;
    uint64_t t_d;
    double weight;
    (*w_a) = 0; found = 0;
    for (i = 0; i < hit_n; i++)///all hits already have the same qn and tn
    {
        if(is_hom_hit(hit[i])) continue;
        interpr_hit(idx, hit[i].s, hit[i].len>>32, &s_uid, &s_beg, &s_end);
        if(s_uid != qn) continue;
        if(!(qs <= s_beg && qe >= s_end)) continue;

        interpr_hit(idx, hit[i].e, (uint32_t)hit[i].len, &e_uid, &e_beg, &e_end);
        if(e_uid != tn) continue;
        if(!(ts <= e_beg && te >= e_end)) continue;

        t_d = get_hic_distance(&hit[i], link, idx, &is_cc);
        // if(t_d == (uint64_t)-1) continue;
        // if(t_d == (uint64_t)-1 && is_cc == 0) continue;
        // if(t_d == (uint64_t)-1 && !dis) continue;

        weight = 1;
        if(dis) weight = get_trans_weight_advance(idx, t_d, dis);

        (*w_a) += weight;
        found = 1;
    }

    return found;
} 


void append_trans_hits(ha_ug_index* idx, hc_links* link, bubble_type* bub, 
trans_idx* dis, pe_hit *hit, uint32_t hit_n, uint64_t *hit_occ, kv_u_trans_t *ref,
kv_u_trans_t *res, uint32_t qn, uint32_t tn)
{
    u_trans_t *r_a = NULL, *p = NULL, *q = NULL;
    uint32_t r_n, i;
    uint64_t occ_q, occ_t;
    double w;

    ///qn -> tn^
    {
        r_a = u_trans_a(*ref, tn);
        r_n = u_trans_n(*ref, tn);
        for (i = 0; i < r_n; i++)
        {
            if(IF_HOM(r_a[i].tn, *bub)) continue;
            if(r_a[i].tn == qn) continue;

            if(!get_trans_interval_weight(idx, link, bub, dis, hit, hit_n, 
            qn, 0, idx->ug->g->seq[qn].len, tn, r_a[i].qs, r_a[i].qe, &w))
            {
                continue;
            } 

            occ_q = hit_occ[qn];
            occ_t = (hit_occ[tn] * ((double)(r_a[i].qe - r_a[i].qs)/(double)(idx->ug->g->seq[tn].len))) + 0.5;
            if(occ_t < 1) occ_t = 1;

            kv_pushp(u_trans_t, *res, &p);
            memset(p, 0, sizeof(u_trans_t));
            p->qn = qn; p->tn = r_a[i].tn;
            p->nw = w/**r_a[i].occ**/;
            p->occ = MIN(occ_q, occ_t);

            kv_pushp(u_trans_t, *res, &q);
            (*q) = (*p); q->qn = p->tn; q->tn = p->qn;
        }
    }

    if(qn == tn) return;

   ///tn -> qn^
    {
        r_a = u_trans_a(*ref, qn);
        r_n = u_trans_n(*ref, qn);
        for (i = 0; i < r_n; i++)
        {
            if(IF_HOM(r_a[i].tn, *bub)) continue;
            if(r_a[i].tn == tn) continue;

            if(!get_trans_interval_weight(idx, link, bub, dis, hit, hit_n, 
            qn, r_a[i].qs, r_a[i].qe, tn, 0, idx->ug->g->seq[tn].len, &w))
            {
                continue;
            } 

            occ_q = (hit_occ[qn] * ((double)(r_a[i].qe - r_a[i].qs)/(double)(idx->ug->g->seq[qn].len))) + 0.5;
            occ_t = hit_occ[tn];
            if(occ_q < 1) occ_q = 1;
            kv_pushp(u_trans_t, *res, &p);
            memset(p, 0, sizeof(u_trans_t));
            p->qn = tn; p->tn = r_a[i].tn;
            p->nw = w/**r_a[i].occ**/;
            p->occ = MIN(occ_q, occ_t);

            kv_pushp(u_trans_t, *res, &q);
            (*q) = (*p); q->qn = p->tn; q->tn = p->qn;
        }
    }
    
}

double merge_u_trans_list(u_trans_t* a, uint32_t a_n)
{
    radix_sort_u_trans_occ(a, a + a_n);
    uint32_t k, l, i;
    double weight, w = 0;
    for (k = 1, l = 0; k <= a_n; ++k) 
    {   
        if (k == a_n || a[k].occ != a[l].occ) 
        {
            for (i = l, weight = 0; i < k; i++)
            {
                weight += a[i].nw;
            }
            
            /*******************************for debug************************************/
            // w += (weight/(double)(a[l].occ));
            w += weight;
            /*******************************for debug************************************/
            l = k;
        }
    }
    return w;
}

void adjust_weight_kv_u_trans_advance(ha_ug_index* idx, kvec_pe_hit* hits, hc_links* link, bubble_type* bub, 
kv_u_trans_t *ta, kv_u_trans_t *ref, trans_idx* dis)
{
    double index_time = yak_realtime();
    uint32_t k, l, m, h_occ;
    uint64_t shif = 64 - idx->uID_bits, qn, tn;
    double w;
    pe_hit *h_a = NULL;

    for (k = m = 0; k < ta->n; k++)
    {
        if(ta->a[k].nw == 0) continue;
        ta->a[k].occ = MIN(hits->occ.a[ta->a[k].qn], hits->occ.a[ta->a[k].tn]);
        if(ta->a[k].occ == 0) continue;
        ta->a[m] = ta->a[k];
        m++;
    }
    ta->n = m;

    for (qn = 0; qn < hits->idx.n; qn++)
    {
        if(IF_HOM(qn, *bub)) continue;
        h_a = hits->a.a + (hits->idx.a[qn]>>32);
        h_occ = (uint32_t)(hits->idx.a[qn]);

        for (k = 1, l = 0; k <= h_occ; ++k) ///same qn
        {
            if (k == h_occ || ((h_a[k].e<<1)>>shif) != ((h_a[l].e<<1)>>shif)) //same tn
            {
                tn = ((h_a[l].e<<1)>>shif);
                if(!IF_HOM(tn, *bub))
                {
                    append_trans_hits(idx, link, bub, dis, h_a+l, k-l, hits->occ.a, ref, ta, qn, tn);
                }
                l = k;
            }
        }
    }
    
    radix_sort_u_trans_m(ta->a, ta->a + ta->n);

    for (k = 1, l = 0, m = 0; k <= ta->n; ++k) 
    {   
        if (k == ta->n || (ta->a[k].qn != ta->a[l].qn || ta->a[k].tn != ta->a[l].tn)) 
        {
            w = merge_u_trans_list(ta->a + l, k - l);
            if(w != 0)
            {
                ta->a[m] = ta->a[l];
                ta->a[m].nw = w;
                ta->a[m].occ = 0;
                m++;
            }
            l = k;
        }
    }
    ta->n = m;
    
    kt_u_trans_t_idx(ta, idx->ug->g->n_seq);
    fprintf(stderr, "[M::%s::%.3f] \n", __func__, yak_realtime()-index_time);
}

void print_kv_weight(kv_u_trans_t *ta, int8_t *s)
{
    uint32_t i;
    // u_trans_t *e = NULL;
    fprintf(stderr, "\n[M::%s]\n", __func__);
    fprintf(stderr, "*********ta->n: %u\n", (uint32_t)ta->n);
    for (i = 0; i < ta->n; i++)
    {
        fprintf(stderr, "+s-utg%.6ul(s:%d)\td-utg%.6ul(s:%d)\tw-%f\n", 
                        ta->a[i].qn+1, s[ta->a[i].qn], ta->a[i].tn+1, s[ta->a[i].tn], ta->a[i].nw);
        /**
        get_u_trans_spec(ta, ta->a[i].tn, ta->a[i].qn, &e, NULL);
        if(e)
        {
            fprintf(stderr, "-d-utg%.6ul->s-utg%.6ul: %f\n", e->qn+1, e->tn+1, e->nw);
        }
        else
        {
            fprintf(stderr, "ERROR");
        }
        **/
    }
}

void print_debug_hc_links(ha_ug_index* idx, bubble_type* bub, hc_links* lk, kv_u_trans_t *ta, kvec_pe_hit* hits)
{
    uint64_t k, len = 0, occ = 0, shif = 64 - idx->uID_bits, beg, end;
    hc_edge *e = NULL;
    for (k = 0; k < hits->a.n; ++k) 
    {
        beg = ((hits->a.a[k].s<<1)>>shif);
        end = ((hits->a.a[k].e<<1)>>shif);
        if(beg == end) continue;
        if(IF_HOM(beg, *bub)) continue;
        if(IF_HOM(end, *bub)) continue;
        len += (hits->a.a[k].len>>32) + ((uint32_t)hits->a.a[k].len);
        occ++;
    }

    fprintf(stderr, "# total Hi-C aligned bases: %lu\n", len);
    fprintf(stderr, "# total Hi-C aligned pairs: %lu\n", occ);
    for (k = 0; k < ta->n; k++)
    {
        e = get_hc_edge(lk, ta->a[k].qn, ta->a[k].tn, 0);
        fprintf(stderr, "s-utg%.6ul\td-utg%.6ul\tD:%lu\tW:%f\n", 
        ta->a[k].qn+1, ta->a[k].tn+1, e->dis == (uint64_t)-1? (uint64_t)-1 : e->dis>>3, ta->a[k].nw);
    }
}
void renew_kv_u_trans(kv_u_trans_t *ta, hc_links *lk, kvec_pe_hit* hits, kv_u_trans_t *ref,
ha_ug_index* idx, bubble_type* bub, int8_t *s, mc_gg_status *sa, uint32_t ignore_dis)
{   
    uint64_t k, i, m, is_comples_weight = 0;
    trans_idx dis;
    kv_init(dis);
    if(bub->round_id > 0 && ignore_dis == 0)
    {
        is_comples_weight = get_trans_rate_function_advance(idx, hits, lk, bub, NULL, s, sa, &dis);
    }

    for (i = 0; i < lk->a.n; i++)
    {
        for (k = m = 0; k < lk->a.a[i].e.n; k++)
        {
            if(lk->a.a[i].e.a[k].del) continue;
            lk->a.a[i].e.a[m] = lk->a.a[i].e.a[k];
            lk->a.a[i].e.a[m].weight = 0;
            lk->a.a[i].e.a[m].occ = 0;
            m++;
        }
        lk->a.a[i].e.n = m;
    }

    if(hits->idx.n == 0) idx_hc_links(hits, idx, bub);

    weight_kv_u_trans(idx, hits, lk, bub, ta, is_comples_weight == 1? &dis : NULL, asm_opt.flag&HA_F_USKEW?0:1);
    // if(bub->round_id == bub->n_round-1)
    // {
    //     print_debug_hc_links(idx, bub, lk, ta, hits);
    // }
    // adjust_weight_kv_u_trans(idx, hits, lk, bub, ta, ref, is_comples_weight == 1? &dis : NULL);
    adjust_weight_kv_u_trans_advance(idx, hits, lk, bub, ta, ref, is_comples_weight == 1? &dis : NULL);
    kv_destroy(dis);
    // print_kv_weight(ta);
}

void print_kv_u_trans(kv_u_trans_t *ta, hc_links* lk, int8_t *s)
{
    uint32_t i;
    u_trans_t *p = NULL;
    hc_edge *e = NULL;

    for (i = 0; i < ta->n; i++)
    {
        p = &(ta->a[i]);
        e = get_hc_edge(lk, p->qn, p->tn, 0);
        fprintf(stderr, "s-utg%.6ul\tS(%d)\td-utg%.6ul\tS(%d)\trev(%u)\td(%lld)\ttw(%f)\n", 
        p->qn+1, s[p->qn], p->tn+1, s[p->tn], p->rev, 
        (e == NULL || e->dis == (uint64_t)-1)? -1 : (long long)(e->dis>>3), p->nw);
    }
    
}

ps_t* init_ps_t(uint64_t seed, uint64_t n)
{
    ps_t *s = NULL; CALLOC(s, 1);
    s->xs = seed;
    CALLOC(s->s, n);
    s->n = n;
    return s;
}

void destory_ps_t(ps_t **s)
{
    free((*s)->s);
    free((*s));
}

uint32_t is_hom_map(uint64_t x, uint32_t len, mc_interval_t *p, uint32_t *p_idx, ha_ug_index* idx)
{
    mc_interval_t *a = NULL;
    uint32_t uid, qs, qe, as, ae, occ, k;
    uint64_t oLen;
    interpr_hit(idx, x, len, &uid, &qs, &qe);

    a = p + p_idx[uid]; 
    occ = p_idx[uid+1] - p_idx[uid];
    for (k = 0; k < occ; k++)
    {
        as = a[k].bS;
        ae = a[k].bE;
        oLen = ((MIN(qe, ae) >= MAX(qs, as))? MIN(qe, ae) - MAX(qs, as) + 1 : 0);
        if(oLen == 0) continue;
        if(oLen > len*0.2) return 1;
    }
    return 0;
}

void update_hits(ha_ug_index* idx, kvec_pe_hit* hits, uint8_t *r_het)
{
    ma_ug_t *ug = idx->ug;
    asg_t *rg = idx->read_g;
    uint32_t k, v, l, offset, l_pos;
    asg_t* nsg = ug->g;
    ma_utg_t *u = NULL;
    mc_interval_t *t = NULL;
    kvec_t(mc_interval_t) p; kv_init(p);
    kvec_t(uint32_t) p_idx; kv_init(p_idx);

    kv_push(uint32_t, p_idx, 0);
    for (v = 0; v < nsg->n_seq; v++)
    {
        u = &(ug->u.a[v]);
        for (k = 1, l = 0, offset = 0, l_pos = 0; k <= u->n; ++k) 
        {   
            if (k == u->n || r_het[u->a[k]>>33] != r_het[u->a[l]>>33])
            {
                if(r_het[u->a[l]>>33] == N_HET)///only keep hom suregions
                {
                    kv_pushp(mc_interval_t, p, &t);
                    t->uID = v;
                    t->hs = r_het[u->a[l]>>33];

                    t->bS = l_pos;
                    t->bE = offset + rg->seq[u->a[k-1]>>33].len - 1;

                    t->nS = l;
                    t->nE = k - 1;
                }
                l = k;
                l_pos = offset + (uint32_t)u->a[k-1];
            }
            offset += (uint32_t)u->a[k-1];
        }
        kv_push(uint32_t, p_idx, p.n);
    }
    

    for (k = 0; k < hits->a.n; ++k)
    {
        hits->a.a[k].id = (uint64_t)-1;
        if(is_hom_map(hits->a.a[k].s, hits->a.a[k].len>>32, p.a, p_idx.a, idx) || 
           is_hom_map(hits->a.a[k].e, (uint32_t)hits->a.a[k].len, p.a, p_idx.a, idx))
        {
            continue;
        }
        hits->a.a[k].id = 0;
    }

    kv_destroy(p); kv_destroy(p_idx);
}

void verbose_het_stat(bubble_type *bub)
{
    uint64_t i, hetBase = 0, homBase = 0;
    for (i = 0; i < bub->ug->g->n_seq; i++)
    {
        if(IF_HOM(i, *bub)) homBase += bub->ug->g->seq[i].len;
        else hetBase += bub->ug->g->seq[i].len;
    }

    fprintf(stderr, "[M::stat] # heterozygous bases: %lu; # homozygous bases: %lu\n", hetBase, homBase);
}

void debug_output_disconnected_hits(ha_ug_index* idx, kv_u_trans_t *ta, kvec_pe_hit *hits, hc_links *link, bubble_type *bub, int8_t *s)
{
    uint32_t k, l, i, m, h_occ;
    uint64_t shif = 64 - idx->uID_bits, qn, tn, u_dis;
    pe_hit *h_a = NULL;
    hc_linkeage *t = NULL;
    u_trans_t *p = NULL;
    kv_u_trans_t k_trans; 
    kv_init(k_trans);


    for (qn = 0; qn < hits->idx.n; qn++)
    {
        if(IF_HOM(qn, *bub)) continue;
        h_a = hits->a.a + (hits->idx.a[qn]>>32);
        h_occ = (uint32_t)(hits->idx.a[qn]);

        for (k = 1, l = 0; k <= h_occ; ++k) ///same qn
        {
            if (k == h_occ || ((h_a[k].e<<1)>>shif) != ((h_a[l].e<<1)>>shif)) //same qn and tn
            {
                tn = ((h_a[l].e<<1)>>shif);
                if(!IF_HOM(tn, *bub) && tn != qn)
                {
                    // t = &(link->a.a[qn]);
                    // for (i = 0, u_dis = (uint64_t)-1; i < t->e.n; i++)
                    // {
                    //     if(t->e.a[i].del || t->e.a[i].uID != tn) continue;
                    //     u_dis = (t->e.a[i].dis ==(uint64_t)-1? (uint64_t)-1 : t->e.a[i].dis>>3);
                    //     break;
                    // }
                    
                    // if(u_dis == (uint64_t)-1)
                    {
                        kv_pushp(u_trans_t, k_trans, &p);
                        p->qn = qn; p->tn = tn; p->occ = (k-l);
                        kv_pushp(u_trans_t, k_trans, &p);
                        p->qn = tn; p->tn = qn; p->occ = (k-l);
                    }
                }
                l = k;
            }
        }
    }


    radix_sort_u_trans_m(k_trans.a, k_trans.a + k_trans.n);

    for (k = 1, l = 0, m = 0; k <= k_trans.n; ++k)
    {
        if (k == k_trans.n || k_trans.a[l].qn != k_trans.a[k].qn || k_trans.a[l].tn != k_trans.a[k].tn) //same qn and tn
        {
            for (i = l, h_occ = 0; i < k; i++)
            {
                h_occ += k_trans.a[i].occ;
            }

            k_trans.a[m] = k_trans.a[l];
            k_trans.a[m].occ = ((uint32_t)-1) - h_occ;
            m++;
            l = k;
        }
    }
    k_trans.n = m;

    radix_sort_u_trans_occ(k_trans.a, k_trans.a + k_trans.n);
    for (i = 0; i < k_trans.n; i++)
    {
        qn = k_trans.a[i].qn;
        tn = k_trans.a[i].tn;

        t = &(link->a.a[qn]);
        for (k = 0, u_dis = (uint64_t)-1; k < t->e.n; k++)
        {
            if(t->e.a[k].del || t->e.a[k].uID != tn) continue;
            u_dis = (t->e.a[k].dis ==(uint64_t)-1? (uint64_t)-1 : t->e.a[k].dis>>3);
            break;
        }
        get_u_trans_spec(ta, qn, tn, &p, NULL);

        fprintf(stderr, "s-utg%.6lul[hap-%d]<--->d-utg%.6lul[hap-%d](occ: %u, ",
        qn + 1, s[qn], tn + 1, s[tn], ((uint32_t)-1) - k_trans.a[i].occ);
        if(p) fprintf(stderr, "weight: %f", p->nw);
        else fprintf(stderr, "weight: NA");
        fprintf(stderr, "):(dis-%lu)\n", u_dis);
    }

    // fprintf(stderr, "########hits########\n");
    // char dir[2] = {'+', '-'};
    // for (k = 0; k < hits->a.n; ++k) 
    // { 
    //     fprintf(stderr, "%c\tutg%.6dl(len-%u)\t%lu\t%c\tutg%.6dl(len-%u)\t%lu\ti:%lu\n", 
    //     dir[hits->a.a[k].s>>63], (int)((hits->a.a[k].s<<1)>>shif)+1, 
    //     idx->ug->g->seq[((hits->a.a[k].s<<1)>>shif)].len, hits->a.a[k].s&idx->pos_mode,
    //     dir[hits->a.a[k].e>>63], (int)((hits->a.a[k].e<<1)>>shif)+1, 
    //     idx->ug->g->seq[((hits->a.a[k].e<<1)>>shif)].len, hits->a.a[k].e&idx->pos_mode,
    //     hits->a.a[k].id);        
    // }
    kv_destroy(k_trans);
}

void resolve_tangles_hic(ha_ug_index *idx, bubble_type *bub, kvec_pe_hit *hits, kv_u_trans_t *ta)
{
    uint32_t i, k, l, m, h_occ;
    uint64_t shif = 64 - idx->uID_bits, qn, tn;
    pe_hit *h_a = NULL;
    u_trans_t *p = NULL;
    identify_bubbles(idx->ug, bub, idx->t_ch->ir_het, &(idx->t_ch->k_trans));
    ta->idx.n = ta->n = 0;
    if(hits->idx.n == 0) idx_hc_links(hits, idx, NULL);
    for (qn = 0; qn < hits->idx.n; qn++)
    {
        h_a = hits->a.a + (hits->idx.a[qn]>>32);
        h_occ = (uint32_t)(hits->idx.a[qn]);

        for (k = 1, l = 0; k <= h_occ; ++k) ///same qn
        {
            if (k == h_occ || ((h_a[k].e<<1)>>shif) != ((h_a[l].e<<1)>>shif)) //same qn and tn
            {
                tn = ((h_a[l].e<<1)>>shif);
                if(tn != qn)
                {
                    kv_pushp(u_trans_t, *ta, &p);
                    p->qn = qn; p->tn = tn; p->occ = (k-l);
                    kv_pushp(u_trans_t, *ta, &p);
                    p->qn = tn; p->tn = qn; p->occ = (k-l);
                }
                l = k;
            }
        }
    }
    radix_sort_u_trans_m(ta->a, ta->a + ta->n);
    for (k = 1, l = 0, m = 0; k <= ta->n; ++k)
    {
        if (k == ta->n || ta->a[l].qn != ta->a[k].qn || ta->a[l].tn != ta->a[k].tn) //same qn and tn
        {
            for (i = l, h_occ = 0; i < k; i++)
            {
                h_occ += ta->a[i].occ;
            }

            qn = ta->a[l].qn;
            tn = ta->a[l].tn;
            ta->a[m] = ta->a[l];
            ta->a[m].occ = h_occ;
            ta->a[m].nw = (h_occ*SCALL)/(MIN(hits->occ.a[qn], hits->occ.a[tn]));
            m++;
            l = k;
        }
    }
    ta->n = m;
    kt_u_trans_t_idx(ta, idx->ug->g->n_seq);
    
    resolve_bubble_chain_by_hic(idx, ta, bub);
    ta->idx.n = ta->n = 0;
}


void print_kv_u_trans_t(kv_u_trans_t *ta, ma_ug_t* ug)
{
    uint32_t i;
    u_trans_t *p = NULL;
    for (i = 0; i < ta->n; i++)
    {
        p = &(ta->a[i]);
        fprintf(stderr, "q-utg%.6ul\tql(%u)\tqs(%u)\tqe(%u)\tt-utg%.6ul\ttl(%u)\tts(%u)\tte(%u)\trev(%u)\tw(%f)\tf(%u)\n", 
        p->qn+1, ug->u.a[p->qn].len, p->qs, p->qe, p->tn+1, ug->u.a[p->tn].len, p->ts, p->te, p->rev, p->nw, p->f);
    }
    fprintf(stderr, "[M::%s::] \n", __func__);
}

void write_ps_t(ps_t *s, const char *fn)
{
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
    sprintf(buf, "%s.hic.pst.bin", fn);
    FILE* fp = fopen(buf, "w");

    fwrite(&(s->xs), sizeof(s->xs), 1, fp);
    fwrite(&(s->n), sizeof(s->n), 1, fp);
    fwrite(s->s, sizeof(int8_t), s->n, fp);

    fclose(fp);
    free(buf);
}

int load_ps_t(ps_t **s, const char *fn)
{
    uint64_t flag = 0;
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
    sprintf(buf, "%s.hic.pst.bin", fn);
    FILE* fp = NULL; 
    fp = fopen(buf, "r"); 
    if(!fp) return 0;
    CALLOC(*s, 1);
    flag += fread(&((*s)->xs), sizeof((*s)->xs), 1, fp);
    flag += fread(&((*s)->n), sizeof((*s)->n), 1, fp);
    MALLOC((*s)->s, (*s)->n);
    flag += fread((*s)->s, sizeof(int8_t), (*s)->n, fp);

    fclose(fp);
    free(buf);
    return 1;
}

int cmp_kv_u_weight(const void * a, const void * b)
{
    if((*(u_trans_t*)a).nw == (*(u_trans_t*)b).nw) return 0;
    return ((*(u_trans_t*)a).nw) > ((*(u_trans_t*)b).nw)?-1:1;
}

typedef struct{
    uint64_t s, e;
} f_chain_t;

typedef struct{
    f_chain_t* a;
    size_t n, m;
    kvec_t(uint64_t) b;
    double tw;
    uint64_t tov;
} kv_f_chain;

void insert_dip_chain(kv_f_chain *x, u_trans_t *p, double LenRate)
{
    uint32_t i;
    uint64_t ovlp = 0, tLen = MIN(x->tov, p->qe-p->qs);
    int64_t dp, old_dp, start = 0;
    f_chain_t *t = NULL;
    double r;
    for (i = 0; i < x->n; i++)
    {
        ovlp += ((MIN(x->a[i].e, p->qe) > MAX(x->a[i].s, p->qs))?
                                        (MIN(x->a[i].e, p->qe) - MAX(x->a[i].s, p->qs)):0);
    }
    if(tLen > 0 && ovlp > tLen*LenRate)
    {
        r = ((double)ovlp)/((double)tLen);
        if(r >= 0.75)
        {
            if(p->nw <= x->tw*0.75) p->del = 1;
        } 
        else
        {
            if(p->nw <= x->tw*0.5) p->del = 1;
        }
    }

    if(!p->del)
    {
        x->tw += p->nw - (((double)ovlp)/((double)(p->qe-p->qs)))*p->nw;   

        kv_push(uint64_t, x->b, p->qs<<1);
        kv_push(uint64_t, x->b, p->qe<<1|1);
        if(x->b.n > 2 && x->b.a[x->b.n-2] < x->b.a[x->b.n-3])
        {
            radix_sort_hc64(x->b.a, x->b.a + x->b.n);
        }
        x->n = 0; x->tov = 0;
        for (i = 0, dp = 0, start = 0; i < x->b.n; ++i) 
        {
            old_dp = dp;
            ///if a[j] is qe
            if (x->b.a[i]&1) 
            {
                --dp;
            }
            else
            {
                ++dp;
            } 


            if (old_dp < 1 && dp >= 1) ///old_dp < dp, b.a[j] is qs
            { 
                ///case 2, a[j] is qs
                start = x->b.a[i]>>1;
            } 
            else if (old_dp >= 1 && dp < 1) ///old_dp > min_dp, b.a[j] is qe
            {
                kv_pushp(f_chain_t, *x, &t);
                t->s = start; t->e = x->b.a[i]>>1;
                x->tov += t->e - t->s;
            }
        }
    }
}

uint32_t get_MAPQ(u_trans_t *a, uint32_t a_n, uint32_t idx, double secondRate)
{
    uint32_t i;
    u_trans_t *p = &(a[idx]);
    uint64_t ovlp = 0;
    double aw = 0;
    for (i = 0; i < a_n; i++)
    {
        if(i == idx) continue;
        ovlp = ((MIN(a[i].qe, p->qe) > MAX(a[i].qs, p->qs))?
                                        (MIN(a[i].qe, p->qe) - MAX(a[i].qs, p->qs)):0);
        // if(p->qn == 447 && p->tn == 5495)
        // {
        //     fprintf(stderr, "ovlp: %lu, secondRate: %f, len: %u\n", ovlp, secondRate, MIN(a[i].qe - a[i].qs, p->qe - p->qs));
        // }
        if(ovlp == 0) continue;
        if(ovlp <= (secondRate*MIN(a[i].qe - a[i].qs, p->qe - p->qs))) continue;
        
        aw += a[i].nw;
    }
    // if(p->qn == 447 && p->tn == 5495) fprintf(stderr, "aw: %f\n", aw);
    aw = aw/p->nw;
    // if(p->qn == 447 && p->tn == 5495) fprintf(stderr, "aw: %f\n", aw);
    if(aw >= 1) return 16;
    if(aw >= 0.75) return 8;
    if(aw >= 0.5) return 4;
    if(aw >= 0.25) return 2;
    return 1;
}

void filter_kv_u_trans_t(kv_u_trans_t *ta, ma_ug_t* ug, double secondRate)
{
    kv_f_chain x; 
    memset(&x, 0, sizeof(kv_f_chain)); x.tw = 0;
    uint32_t k, i, n;
    u_trans_t *a = NULL;
    for (k = 0; k < ta->idx.n; k++)
    {
        a = u_trans_a(*ta, k);
        n = u_trans_n(*ta, k);
        if(n == 0) continue;
        qsort(a, n, sizeof(u_trans_t), cmp_kv_u_weight);
        x.n = 0; x.tov = 0; x.tw = 0; x.b.n = 0;
        for (i = 0; i < n; i++)
        {
            if(a[i].del) continue;
            insert_dip_chain(&x, &(a[i]), secondRate);
        }
    }
    free(x.a); free(x.b.a);
    kt_u_trans_t_idx(ta, ug->g->n_seq);
    kt_u_trans_t_simple_symm(ta, ug->g->n_seq, 0);


    /**
    for (k = 0; k < ta->idx.n; k++)
    {
        a = u_trans_a(*ta, k);
        n = u_trans_n(*ta, k);
        if(n == 0) continue;
        // radix_sort_k_trans_qs(a, a+n);
        for (i = 0; i < n; i++)
        {
            a[i].occ = get_MAPQ(a, n, i, secondRate);
            // if(a[i].occ != 1)
            // {
            //     uint32_t debug_i;
            //     fprintf(stderr, "\nq-utg%.6ul\tt-utg%.6ul\tocc-%u\n", a[i].qn+1, a[i].tn+1, a[i].occ);


            //     for (debug_i = 0; debug_i < n; debug_i++)
            //     {
            //         fprintf(stderr, "q-utg%.6ul\tqs(%u)\tqe(%u)\tt-utg%.6ul\tts(%u)\tte(%u)\trev(%u)\tw(%f)\tf(%u)\n", 
            //         a[debug_i].qn+1, a[debug_i].qs, a[debug_i].qe, a[debug_i].tn+1, a[debug_i].ts, 
            //         a[debug_i].te, a[debug_i].rev, a[debug_i].nw, a[debug_i].f);
            //     }
            // }
        }
    }

    for (k = 0; k < ta->n; k++)
    {
        if(ta->a[k].del || ta->a[k].qn > ta->a[k].tn) continue;
        get_u_trans_spec(ta, ta->a[k].tn, ta->a[k].qn, &r_a, NULL);
        a = &(ta->a[k]); 
        a->occ = r_a->occ = MAX(a->occ, r_a->occ);
    }
    **/
}

uint64_t check_ovlp(f_chain_t *o, uint64_t on, uint64_t qs, uint64_t qe)
{
    uint64_t i, ovlp;
    for (i = 0; i < on; i++)
    {
        ovlp = ((MIN(o[i].e, qe) > MAX(o[i].s, qs))?
                                        (MIN(o[i].e, qe) - MAX(o[i].s, qs)):0);
        if(ovlp) return 1;
    }
    return 0;
}

void flter_by_cov(ha_ug_index* idx, kvec_pe_hit *hits, int64_t min_dp)
{
    uint64_t i, m, k, l, pos_bits = 64 - idx->uID_bits - 1, uid, pM = ((uint64_t)-1)>>idx->uID_bits, *a = NULL, a_n;
    uint32_t s_uid, s_beg, s_end, e_uid, e_beg, e_end;
    int64_t dp, old_dp, start = 0;
    f_chain_t *t = NULL;
    kvec_t(uint64_t) b; kv_init(b);
    kvec_t(f_chain_t) x; kv_init(x);
    uint64_t *index = NULL; CALLOC(index, idx->ug->u.n);
    for (i = 0; i < hits->a.n; i++)
    {
        interpr_hit(idx, hits->a.a[i].s, hits->a.a[i].len>>32, &s_uid, &s_beg, &s_end);
        m = s_uid; m <<= pos_bits; m += s_beg; m <<= 1; kv_push(uint64_t, b, m);
        m = s_uid; m <<= pos_bits; m += s_end; m <<= 1; m += 1; kv_push(uint64_t, b, m);


        interpr_hit(idx, hits->a.a[i].e, (uint32_t)hits->a.a[i].len, &e_uid, &e_beg, &e_end);
        m = e_uid; m <<= pos_bits; m += e_beg; m <<= 1; kv_push(uint64_t, b, m);
        m = e_uid; m <<= pos_bits; m += e_end; m <<= 1; m += 1; kv_push(uint64_t, b, m);
    }
    radix_sort_hc64(b.a, b.a + b.n);
    for (k = 1, l = 0; k <= b.n; ++k) 
    {
        if (k == b.n || (b.a[k]>>(64 - idx->uID_bits)) != (b.a[l]>>(64 - idx->uID_bits)))
        {
            uid = (b.a[l]>>(64 - idx->uID_bits));
            a = b.a + l; a_n = k - l; index[uid] = x.n;
            for (i = 0, dp = 0, start = 0; i < a_n; ++i) 
            {
                old_dp = dp;
                ///if a[j] is qe
                if (a[i]&1) --dp;
                else ++dp;
    
                if (old_dp < min_dp && dp >= min_dp) ///old_dp < dp, b.a[j] is qs
                { 
                    ///case 2, a[j] is qs
                    start = (a[i]&pM)>>1;
                } 
                else if (old_dp >= min_dp && dp < min_dp) ///old_dp > min_dp, b.a[j] is qe
                {
                    kv_pushp(f_chain_t, x, &t);
                    t->s = start; t->e = (a[i]&pM)>>1;
                }
            }

            index[uid] |= ((uint64_t)(x.n - index[uid]))<<32;
            l = k;
        }
    }
    



    for (i = m = 0; i < hits->a.n; i++)
    {
        interpr_hit(idx, hits->a.a[i].s, hits->a.a[i].len>>32, &s_uid, &s_beg, &s_end);
        if(!check_ovlp(x.a+((uint32_t)index[s_uid]), index[s_uid]>>32, s_beg, s_end)) continue;

        interpr_hit(idx, hits->a.a[i].e, (uint32_t)hits->a.a[i].len, &e_uid, &e_beg, &e_end);
        if(!check_ovlp(x.a+((uint32_t)index[e_uid]), index[e_uid]>>32, e_beg, e_end)) continue;
        hits->a.a[m] = hits->a.a[i];
        m++;
    }

    fprintf(stderr, "[M::%s::] # old Hi-C pairs: %lu, # new Hi-C pairs: %lu\n", 
                                                                            __func__, hits->a.n, m);
    hits->a.n = m;

    kv_destroy(b); kv_destroy(x); free(index);
    
}

void set_tag_pre_read(uint64_t *r_tag, uint64_t id, uint64_t off, uint64_t pos, ma_hit_t_alloc* sources)
{
    uint64_t i, tn;
    ma_hit_t *h = NULL;
    if(r_tag[id] == (uint64_t)-1) r_tag[id] = 0;
    r_tag[id] += off;
    // if(id == 2531247) fprintf(stderr, "***id: %lu, off: %lu, pos: %lu\n", id, off, pos);
    for (i = 0; i < (uint64_t)(sources[id].length); i++)
    {
        h = &(sources[id].buffer[i]);
        if(!h->el) continue;
        tn = Get_tn((*h));
        if(pos >= Get_qs((*h)) && pos < Get_qe((*h)))
        {
            // if(tn == 2531247) fprintf(stderr, "###id: %lu, off: %lu, tn: %lu, pos: %lu\n", id, off, tn, pos);
            if(r_tag[tn] == (uint64_t)-1) r_tag[tn] = 0;
            r_tag[tn] += off;
        }
        
    }
}

void tag_reads(ha_ug_index* idx, kvec_pe_hit *u_hits, bubble_type* bub, int8_t *s, ma_hit_t_alloc* sources)
{
    kvec_pe_hit *r_hits = get_r_hits_for_trio(u_hits, idx->read_g, idx->ug, bub, idx->uID_bits, idx->pos_mode);
    uint64_t *r_tag = NULL, i, k, srid, erid, /**rid, occ, a_occ, **/flag = AMBIGU, cons = 0, incons = 0/**, max, min**/; 
    ma_utg_t *u = NULL;
    ma_ug_t* ug = idx->ug;
    int8_t ss, es, sr;
    CALLOC(r_tag, idx->read_g->n_seq);

    cons = 0, incons = 0; 
    memset(r_tag, -1, sizeof(uint64_t)*idx->read_g->n_seq);
    for (i = 0; i < r_hits->a.n; i++)
    {
        srid = get_hit_suid(*r_hits, i); 
        erid = get_hit_euid(*r_hits, i); 
        ss = s[r_hits->a.a[i].id>>32];
        es = s[(uint32_t)r_hits->a.a[i].id];
        if(ss == 0 || es == 0) continue;
        if((r_hits->a.a[i].id>>32) == ((uint32_t)r_hits->a.a[i].id)) continue;

        set_tag_pre_read(r_tag, srid, es < 0? 1 : ((uint64_t)1<<32), get_hit_spos(*r_hits, i), sources);
        set_tag_pre_read(r_tag, erid, ss < 0? 1 : ((uint64_t)1<<32), get_hit_epos(*r_hits, i), sources);
    
        if(ss == es) cons++;
        else incons++;
    }       
    
    /**
    for (k = 0; k < ug->u.n; k++)
    {
        if(ug->g->seq[k].del) continue;
        u = &(ug->u.a[k]); occ = a_occ = 0;
        for (i = 0; i < u->n; i++)
        {
            rid = u->a[i]>>33; 
            if(r_tag[rid] == (uint64_t)-1) continue;
            max = MAX((r_tag[rid]>>32), ((uint32_t)r_tag[rid]));
            min = MIN((r_tag[rid]>>32), ((uint32_t)r_tag[rid]));
            ///if(max <= (total*0.7) && max != total)
            if(max - min <= max *0.1)
            {
                fprintf(stderr, "max-%lu, min-%lu\n", max, min);
                a_occ++;
            }
            occ++;
        }
        if(occ == 0) continue;
        if(a_occ > occ*0.85) s[k] = 0, filter++, fprintf(stderr, "k: %lu, a_occ: %lu, occ: %lu\n", k, a_occ, occ);
    }

    fprintf(stderr, "[M::%s::] # consistent Hi-C pairs: %lu, # inconsistent Hi-C pairs: %lu, # filter: %lu\n", 
                                                                            __func__, cons, incons, filter);
    
    memset(r_tag, -1, sizeof(uint64_t)*idx->read_g->n_seq);
    for (i = 0; i < r_hits->a.n; i++)
    {
        srid = get_hit_suid(*r_hits, i); 
        erid = get_hit_euid(*r_hits, i); 
        ss = s[r_hits->a.a[i].id>>32];
        es = s[(uint32_t)r_hits->a.a[i].id];
        if(ss == 0 || es == 0) continue;
        if((r_hits->a.a[i].id>>32) == ((uint32_t)r_hits->a.a[i].id)) continue;

        set_tag_pre_read(r_tag, srid, es > 0? 1 : ((uint64_t)1<<32), get_hit_spos(*r_hits, i), sources);
        set_tag_pre_read(r_tag, erid, ss > 0? 1 : ((uint64_t)1<<32), get_hit_epos(*r_hits, i), sources);
    }
    **/  
    fprintf(stderr, "[M::%s::] # consistent Hi-C pairs: %lu, # inconsistent Hi-C pairs: %lu\n", 
                                                                            __func__, cons, incons);

    cons = 0, incons = 0; 
    memset(R_INF.trio_flag, AMBIGU, R_INF.total_reads * sizeof(uint8_t));
    for (i = 0; i < ug->g->n_seq; i++)
    {
        if(ug->g->seq[i].del) continue;

        flag = 0;
        if(s[i] == 0) continue;
        flag = (s[i] > 0? FATHER:MOTHER);

        u = &ug->u.a[i];
        if(u->m == 0) continue;
        for (k = 0; k < u->n; k++)
        {
            if(r_tag[u->a[k]>>33] == (uint64_t)-1) continue;
            if((r_tag[u->a[k]>>33]>>32) == ((uint32_t)r_tag[u->a[k]>>33])) continue;
            sr = (r_tag[u->a[k]>>33]>>32) > ((uint32_t)r_tag[u->a[k]>>33])?1:-1;
            if(sr == s[i])
            {
                R_INF.trio_flag[u->a[k]>>33] = flag;
                fprintf(stderr, "*rid: %lu, s[i]: %d, 1-occ: %u, (-1)-occ: %u\n", u->a[k]>>33, s[i], (uint32_t)(r_tag[u->a[k]>>33]>>32), ((uint32_t)r_tag[u->a[k]>>33]));
                // cons_occ += (r_tag[u->a[k]>>33]>>32) + ((uint32_t)r_tag[u->a[k]>>33]);
                cons++;
            }
            else
            {
                // R_INF.trio_flag[u->a[k]>>33] = -flag;
                fprintf(stderr, "#rid: %lu, s[i]: %d, 1-occ: %u, (-1)-occ: %u\n", u->a[k]>>33, s[i], (uint32_t)(r_tag[u->a[k]>>33]>>32), ((uint32_t)r_tag[u->a[k]>>33]));
                // incons_occ += (r_tag[u->a[k]>>33]>>32) + ((uint32_t)r_tag[u->a[k]>>33]);
                incons++;
            }
        }
    }

    free(r_tag);
    kv_destroy(r_hits->a);
    kv_destroy(r_hits->idx);
    kv_destroy(r_hits->occ);
    free(r_hits);
    fprintf(stderr, "[M::%s::] # consistent reads: %lu, # inconsistent reads: %lu\n", __func__, cons, incons);
}

void renew_idx_para(ha_ug_index* idx, ma_ug_t* ug)
{
    for (idx->uID_bits=1; (uint64_t)(1<<idx->uID_bits)<(uint64_t)ug->u.n; idx->uID_bits++);
    idx->pos_bits = 64 - idx->uID_bits - 1;
    idx->uID_mode = (((uint64_t)-1) << (64-idx->uID_bits))>>1;
    idx->pos_mode = ((uint64_t)-1) >> (64-idx->pos_bits);
    idx->rev_mode = ((uint64_t)1) << 63;
}

uint32_t get_oe_occ(uint32_t qn, uint32_t tn, kvec_pe_hit* hits, ha_ug_index* idx)
{
    uint64_t shif = 64 - idx->uID_bits, occ = 0;
    pe_hit *h_a = hits->a.a + (hits->idx.a[qn]>>32);
    uint32_t h_occ = (uint32_t)(hits->idx.a[qn]), i;
    for (i = 0; i < h_occ; i++) {
        if(((h_a[i].s<<1)>>shif)==qn && ((h_a[i].e<<1)>>shif)==tn) occ++;
    }
    return occ;
}

void optimize_u_trans(kv_u_trans_t *ovlp, kvec_pe_hit* hits, ha_ug_index* idx)
{
    if(hits->idx.n == 0) idx_hc_links(hits, idx, NULL);
    uint64_t i, m, occ;
    u_trans_t *x = NULL, *p = NULL;
    kv_u_trans_t k_trans; 
    kv_init(k_trans); kv_init(k_trans.idx);
    for (i = 0; i < ovlp->n; i++){
        x = &(ovlp->a[i]);
        if(x->qn > x->tn) continue;
        if(x->f != RC_2 || x->del) continue;
        occ = get_oe_occ(x->qn, x->tn, hits, idx) + get_oe_occ(x->tn, x->qn, hits, idx);
        kv_pushp(u_trans_t, k_trans, &p);
        (*p) = (*x); p->nw = (x->nw*(1-(((double)(occ<<1))/((double)(hits->occ.a[x->qn]+hits->occ.a[x->tn])))));
        if(p->nw < 0) fprintf(stderr, "ERROR-nw\n");
        if(p->nw == 0) p->nw = x->nw*0.005;
        if(p->nw == 0) {
            k_trans.n--;
        }
        else {
            kv_pushp(u_trans_t, k_trans, &p);
            (*p) = k_trans.a[k_trans.n-2];
            p->qn = k_trans.a[k_trans.n-2].tn; p->qs = k_trans.a[k_trans.n-2].ts; p->qe = k_trans.a[k_trans.n-2].te;
            p->tn = k_trans.a[k_trans.n-2].qn; p->ts = k_trans.a[k_trans.n-2].qs; p->te = k_trans.a[k_trans.n-2].qe;
        }
    }
    kt_u_trans_t_idx(&k_trans, idx->ug->g->n_seq);
    mc_solve(NULL, NULL, &k_trans, idx->ug, idx->read_g, 0.8, R_INF.trio_flag, 1, NULL, 1, NULL, NULL, 1);
    for (i = m = 0; i < ovlp->n; i++){
        x = &(ovlp->a[i]);
        if(x->del) continue;
        if(x->f == RC_2){
            get_u_trans_spec(&k_trans, x->qn, x->tn, &p, NULL);
            if(!p) continue;
        }
        ovlp->a[m++] = ovlp->a[i];
    }
    ovlp->n = m;
    kt_u_trans_t_idx(ovlp, idx->ug->g->n_seq);
    free(hits->idx.a); hits->idx.a = NULL; hits->idx.n = hits->idx.m = 0;
    free(hits->occ.a); hits->occ.a = NULL; hits->occ.n = hits->occ.m = 0;
    kv_destroy(k_trans); kv_destroy(k_trans.idx);
}

int hic_short_align(const enzyme *fn1, const enzyme *fn2, ha_ug_index* idx, ug_opt_t *opt, kvec_pe_hit **rhits)
{
    double index_time = yak_realtime();
    sldat_t sl;
    sl.idx = idx;
    sl.t_ch = idx->t_ch;
    sl.chunk_size = 20000000;
    sl.n_thread = asm_opt.thread_num;
    sl.total_base = sl.total_pair = 0;
    idx->hap_cnt = asm_opt.hap_occ;
    kv_init(sl.hits.a); kv_init(sl.hits.idx); kv_init(sl.hits.occ);


    if(!load_hc_hits(&sl.hits, idx->ug, asm_opt.output_file_name))
    {
        alignment_worker_pipeline(&sl, fn1, fn2);
        write_hc_hits(&sl.hits, idx->ug, asm_opt.output_file_name);
    }
    sl.hits.uID_bits = idx->uID_bits; sl.hits.pos_mode = idx->pos_mode;
    optimize_u_trans(&(idx->t_ch->k_trans), &sl.hits, idx);
    filter_kv_u_trans_t(&(idx->t_ch->k_trans), idx->ug, 0.5);
    // print_kv_u_trans_t(&(idx->t_ch->k_trans), idx->ug);
    // flter_by_cov(idx, &sl.hits, 2);
    // update_hits(idx, &sl.hits, idx->t_ch->is_r_het);
    ///debug_hc_hits_v14(&sl.hits, asm_opt.output_file_name, sl.idx);
    ////dedup_hits(&(sl.hits), sl.idx);   
    ///write_hc_hits_v14(&sl.hits, asm_opt.output_file_name);
    
    if(asm_opt.misjoin_len > 0)
    {
        update_switch_unitig(idx->ug, idx->read_g, &(sl.hits), &(idx->t_ch->k_trans), 10, 20, asm_opt.misjoin_len, 0.15);
        renew_idx_para(idx, idx->ug);
    }
    // print_hits_simp(idx, &sl.hits);
    // print_kv_u_trans_t(&(idx->t_ch->k_trans));

    hc_links link;
    init_hc_links(&link, idx->ug->g->n_seq, idx->t_ch);
    ///H_partition hap;
    bubble_type bub; 
    kv_u_trans_t k_trans; 
    kv_init(k_trans); kv_init(k_trans.idx);
    ps_t *s = NULL;
    mb_nodes_t u; 
    kv_init(u.bid); kv_init(u.idx); kv_init(u.u);
    memset(&bub, 0, sizeof(bubble_type));
    bub.round_id = 0; bub.n_round = asm_opt.n_weight;

    resolve_tangles_hic(idx, &bub, &sl.hits, &k_trans);
    measure_distance(idx, idx->ug, &sl.hits, &link, &bub, &(idx->t_ch->k_trans));
    // if((asm_opt.flag & HA_F_VERBOSE_GFA) && load_ps_t(&s, asm_opt.output_file_name))
    // {
    //     bub.round_id = bub.n_round;
    //     label_unitigs_sm(s->s, NULL, idx->ug);
    //     goto skip_flipping;
    // }
    s = init_ps_t(11, idx->ug->g->n_seq);
    for (bub.round_id = 0; bub.round_id < bub.n_round; bub.round_id++)
    {
        // identify_bubbles(idx->ug, &bub, idx->t_ch->is_r_het, &(idx->t_ch->k_trans));
        renew_kv_u_trans(&k_trans, &link, &sl.hits, &(idx->t_ch->k_trans), idx, &bub, s->s, NULL, 0);
        // if(bub.round_id == 0) init_phase(idx, &k_trans, &bub, s); 
        // update_trans_g(idx, &k_trans, &bub);
        /*******************************for debug************************************/
        mc_solve(NULL, NULL, &k_trans, idx->ug, idx->read_g, 0.8, R_INF.trio_flag, 
        (bub.round_id == 0? 1 : 0), s->s, 1, /**&bub**/NULL, &(idx->t_ch->k_trans), 0);
        /*******************************for debug************************************/
        label_unitigs_sm(s->s, NULL, idx->ug);

        /*******************************for debug************************************/
        // if(bub.round_id == bub.n_round - 1)
        // {
        //     debug_output_disconnected_hits(idx, &k_trans, &sl.hits, &link, &bub, s->s);
        // }
        /*******************************for debug************************************/
        /**
        init_hic_advance((ha_ug_index*)sl.idx, &sl.hits, &link, &bub, &hap, 0);
        reset_H_partition(&hap, (bub.round_id == 0? 1 : 0));
        init_contig_partition(&hap, idx, &bub, &link);
        phasing_improvement(&hap, &(hap.g_p), idx, &bub, &link);
        label_unitigs(&(hap.g_p), idx->ug);
        **/
    }
    // write_ps_t(s, asm_opt.output_file_name);
    // skip_flipping:
    verbose_het_stat(&bub);

    if(rhits) (*rhits) = get_r_hits_order(&sl.hits, idx->uID_bits, idx->pos_mode, idx->read_g, idx->ug, &bub);

    // tag_reads(idx, &sl.hits, &bub, s->s, opt->sources);

    // print_kv_weight(&k_trans, s->s);

    // horder_t *ho = init_horder_t(&sl.hits, idx->uID_bits, idx->pos_mode, idx->read_g, idx->ug, &bub, &(idx->t_ch->k_trans), opt, 3);

    ///print_hc_links(&link, 0, &hap);
    // print_kv_u_trans(&k_trans, &link, s->s);


    ///print_bubbles(idx->ug, &bub, sl.hits.a.n?&sl.hits:NULL, idx->link, idx);
    ///print_hits(idx, &sl.hits, fn1);
    

    ///print_debug_bubble_graph(&bub, idx->ug, asm_opt.output_file_name);
    // print_bubble_chain(&bub);
    // destory_contig_partition(&hap);
    // destory_horder_t(&ho);
    kv_destroy(sl.hits.a); kv_destroy(sl.hits.idx); kv_destroy(sl.hits.occ);
    destory_hc_links(&link);
    kv_destroy(k_trans); kv_destroy(k_trans.idx);
    destory_ps_t(&s);
    kv_destroy(u.bid); kv_destroy(u.idx); kv_destroy(u.u);
    destory_bubbles(&bub);
    return 1;
    

    print_bubbles(idx->ug, &bub, sl.hits.a.n?&sl.hits:NULL, &link, idx);
    collect_hc_reverse_links(&link, idx->ug, &bub);
    normalize_hc_links(&link);
    min_cut_t* cut = clean_hap(&link, &bub, idx->ug);
    ///print_bubbles(idx->ug, &bub, NULL, &link, idx);
    G_partition* gp = clean_bubbles(&link, &bub, cut, idx->ug);
    ///print_hc_links(&link);
    destory_min_cut_t(cut); free(cut);
    destory_G_partition(gp); free(gp);
    destory_bubbles(&bub);
    
    fprintf(stderr, "[M::%s::%.3f] processed %lu pairs; %lu bases\n", __func__, yak_realtime()-index_time, sl.total_pair, sl.total_base);
    return 1;
}

spg_t *hic_short_pre_align(const enzyme *fn1, const enzyme *fn2, ha_ug_index* idx, ug_opt_t *opt, kvec_pe_hit **rhits)
{
    // double index_time = yak_realtime();
    sldat_t sl;
    sl.idx = idx;
    sl.t_ch = idx->t_ch;
    sl.chunk_size = 20000000;
    sl.n_thread = asm_opt.thread_num;
    sl.total_base = sl.total_pair = 0;
    idx->hap_cnt = asm_opt.hap_occ;
    kv_init(sl.hits.a); kv_init(sl.hits.idx); kv_init(sl.hits.occ);

    if(!load_hc_hits(&sl.hits, idx->ug, asm_opt.output_file_name))
    {
        alignment_worker_pipeline(&sl, fn1, fn2);
        fprintf(stderr, "sb0sb\n");
        write_hc_hits(&sl.hits, idx->ug, asm_opt.output_file_name);
        fprintf(stderr, "sb1sb\n");
    }
    sl.hits.uID_bits = idx->uID_bits; sl.hits.pos_mode = idx->pos_mode;
    kv_u_trans_t k_trans; 
    kv_init(k_trans); kv_init(k_trans.idx);
    bubble_type bub; 
    memset(&bub, 0, sizeof(bubble_type));
    bub.round_id = 0; bub.n_round = asm_opt.n_weight;
    fprintf(stderr, "sb2sb\n");
    resolve_tangles_hic(idx, &bub, &sl.hits, &k_trans);
    fprintf(stderr, "sb3sb\n");
    spg_t *scg = horder_utg(&sl.hits, idx->uID_bits, idx->pos_mode, idx->read_g, idx->ug, &bub, opt);
    fprintf(stderr, "sb4sb\n");
    if(rhits){
        CALLOC(*rhits, 1);
        (**rhits) = sl.hits;
        sl.hits.a.a = NULL;
        sl.hits.idx.a = NULL;
        sl.hits.occ.a = NULL;
    }
    kv_destroy(sl.hits.a); kv_destroy(sl.hits.idx); kv_destroy(sl.hits.occ);
    kv_destroy(k_trans); kv_destroy(k_trans.idx);
    destory_bubbles(&bub);
    return scg;
}

int load_psg_t(psg_t **sg, const char *fn)
{
    uint64_t flag = 0;
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
    sprintf(buf, "%s.hic.pst.bin", fn);
    FILE* fp = NULL; 
    fp = fopen(buf, "r"); 
    if(!fp) return 0;
    CALLOC(*sg, 1);
    flag += fread(&((*sg)->xs), sizeof((*sg)->xs), 1, fp);
    flag += fread(&((*sg)->sg.n), sizeof((*sg)->sg.n), 1, fp);
    (*sg)->sg.m = (*sg)->sg.n;
    MALLOC((*sg)->sg.a, (*sg)->sg.n);
    flag += fread((*sg)->sg.a, sizeof(mc_gg_status), (*sg)->sg.n, fp);

    fclose(fp);
    free(buf);
    return 1;
}

void write_psg_t(psg_t *sg, const char *fn)
{
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
    sprintf(buf, "%s.hic.pst.bin", fn);
    FILE* fp = fopen(buf, "w");

    fwrite(&(sg->xs), sizeof(sg->xs), 1, fp);
    fwrite(&(sg->sg.n), sizeof(sg->sg.n), 1, fp);
    fwrite(sg->sg.a, sizeof(mc_gg_status), sg->sg.n, fp);

    fclose(fp);
    free(buf);
}

psg_t* init_psg_t(uint64_t seed, ma_ug_t* ug, asg_t* rg, ug_opt_t *opt)
{
    psg_t *s = NULL; CALLOC(s, 1);
    s->xs = seed;
    
    kv_gg_status *sa = init_mc_gg_status(ug, rg, opt->coverage_cut, opt->sources, opt->ruIndex, 
    asm_opt.hom_global_coverage_set?asm_opt.hom_global_coverage:((double)asm_opt.hom_global_coverage)/((double)HOM_PEAK_RATE), 
    asm_opt.polyploidy);

    s->sg = *sa;
    free(sa);
    return s;
}

void destory_psg_t(psg_t **s)
{
    free((*s)->sg.a);
    free((*s));
}

int hic_short_align_poy(const enzyme *fn1, const enzyme *fn2, ha_ug_index* idx, ug_opt_t *opt)
{
    double index_time = yak_realtime();
    sldat_t sl;
    kvec_hc_edge back_hc_edge;
    kv_init(back_hc_edge.a);
    sl.idx = idx;
    sl.t_ch = idx->t_ch;
    sl.chunk_size = 20000000;
    sl.n_thread = asm_opt.thread_num;
    sl.total_base = sl.total_pair = 0;
    idx->hap_cnt = asm_opt.hap_occ;
    kv_init(sl.hits.a); kv_init(sl.hits.idx); kv_init(sl.hits.occ);

    
    if(!load_hc_hits(&sl.hits, idx->ug, asm_opt.output_file_name))
    {
        alignment_worker_pipeline(&sl, fn1, fn2);
        write_hc_hits(&sl.hits, idx->ug, asm_opt.output_file_name);
    }

    hc_links link;
    init_hc_links(&link, idx->ug->g->n_seq, idx->t_ch);
    ///H_partition hap;
    bubble_type bub; 
    kv_u_trans_t k_trans; 
    kv_init(k_trans); kv_init(k_trans.idx);
    psg_t *s = NULL;
    mb_nodes_t u; 
    kv_init(u.bid); kv_init(u.idx); kv_init(u.u);
    memset(&bub, 0, sizeof(bubble_type));
    bub.round_id = 0; bub.n_round = asm_opt.n_weight;

    resolve_tangles_hic(idx, &bub, &sl.hits, &k_trans);
    measure_distance(idx, idx->ug, &sl.hits, &link, &bub, &(idx->t_ch->k_trans));
    if((asm_opt.flag & HA_F_VERBOSE_GFA) && load_psg_t(&s, asm_opt.output_file_name))
    {
        bub.round_id = bub.n_round;
        goto skip_flipping;
    }
    s = init_psg_t(11, idx->ug, idx->read_g, opt);
    for (bub.round_id = 0; bub.round_id < bub.n_round; bub.round_id++)
    {
        // identify_bubbles(idx->ug, &bub, idx->t_ch->is_r_het, &(idx->t_ch->k_trans));
        renew_kv_u_trans(&k_trans, &link, &sl.hits, &(idx->t_ch->k_trans), idx, &bub, NULL, s->sg.a, 0);
        // if(bub.round_id == 0) init_phase(idx, &k_trans, &bub, s); 
        // update_trans_g(idx, &k_trans, &bub);
        /*******************************for debug************************************/
        mc_solve_general(&k_trans, idx->ug->u.n, &(s->sg), asm_opt.polyploidy, 0, 1);
        /*******************************for debug************************************/

        /*******************************for debug************************************/
        // if(bub.round_id == bub.n_round - 1)
        // {
        //     debug_output_disconnected_hits(idx, &k_trans, &sl.hits, &link, &bub, s->s);
        // }
        /*******************************for debug************************************/
    }
    if((asm_opt.flag & HA_F_VERBOSE_GFA)) write_psg_t(s, asm_opt.output_file_name);

    skip_flipping:
    verbose_het_stat(&bub);
    
    if(asm_opt.polyploidy == 2) label_unitigs_sm(NULL, s->sg.a, idx->ug);

    // print_kv_weight(&k_trans);

    // horder_t *ho = init_horder_t(&sl.hits, idx->uID_bits, idx->pos_mode, idx->read_g, idx->ug, &bub, &(idx->t_ch->k_trans), opt, 3);

    ///print_hc_links(&link, 0, &hap);
    // print_kv_u_trans(&k_trans, &link, s->s);


    ///print_bubbles(idx->ug, &bub, sl.hits.a.n?&sl.hits:NULL, idx->link, idx);
    ///print_hits(idx, &sl.hits, fn1);
    

    ///print_debug_bubble_graph(&bub, idx->ug, asm_opt.output_file_name);
    // print_bubble_chain(&bub);
    // destory_contig_partition(&hap);
    // destory_horder_t(&ho);
    kv_destroy(back_hc_edge.a);
    kv_destroy(sl.hits.a);
    kv_destroy(sl.hits.idx);
    kv_destroy(sl.hits.occ);
    destory_hc_links(&link);
    kv_destroy(k_trans); 
    kv_destroy(k_trans.idx);
    destory_psg_t(&s);
    kv_destroy(u.bid); kv_destroy(u.idx); kv_destroy(u.u);
    fprintf(stderr, "[M::%s::%.3f] processed %lu pairs; %lu bases\n", __func__, yak_realtime()-index_time, sl.total_pair, sl.total_base);
    return 1;
}


void hic_analysis(ma_ug_t *ug, asg_t* read_g, trans_chain* t_ch, ug_opt_t *opt, uint32_t is_poy, kvec_pe_hit **rhits)
{
    ug_index = NULL;
    int exist = (asm_opt.load_index_from_disk? 
                    load_hc_pt_index(&ug_index, ug, asm_opt.output_file_name) : 0);
    if(exist == 0) ug_index = build_unitig_index(ug, asm_opt.hic_mer_length, asm_opt.hap_occ, 0, asm_opt.thread_num);
    if(exist == 0) write_hc_pt_index(ug_index, asm_opt.output_file_name);
    ug_index->ug = ug;
    ug_index->read_g = read_g;
    ug_index->t_ch = t_ch;
    ///test_unitig_index(ug_index, ug);
    if(!is_poy) hic_short_align(asm_opt.hic_reads[0], asm_opt.hic_reads[1], ug_index, opt, rhits);
    else hic_short_align_poy(asm_opt.hic_reads[0], asm_opt.hic_reads[1], ug_index, opt);
    
    
    destory_hc_pt_index(ug_index);free(ug_index);
}

spg_t *hic_pre_analysis(ma_ug_t *ug, asg_t* read_g, trans_chain* t_ch, ug_opt_t *opt, kvec_pe_hit **rhits)
{
    ug_index = NULL;
    int exist = (asm_opt.load_index_from_disk? 
                    load_hc_pt_index(&ug_index, ug, asm_opt.output_file_name) : 0);
    if(exist == 0) ug_index = build_unitig_index(ug, asm_opt.hic_mer_length, asm_opt.hap_occ, 0, asm_opt.thread_num);
    if(exist == 0) write_hc_pt_index(ug_index, asm_opt.output_file_name);
    ug_index->ug = ug;
    ug_index->read_g = read_g;
    ug_index->t_ch = t_ch;
    return hic_short_pre_align(asm_opt.hic_reads[0], asm_opt.hic_reads[1], ug_index, opt, rhits);    
}


void init_ug_idx(ma_ug_t *ug, uint64_t k, uint64_t up_bound, uint64_t low_bound, uint64_t build_idx)
{
    ug_index = NULL;
    if(build_idx)
    {
        ug_index = build_unitig_index(ug, k, up_bound, low_bound, asm_opt.thread_num);
    }
}

void des_ug_idx()
{
    destory_hc_pt_index(ug_index);
}

uint64_t count_unique_k_mers(char *r, uint64_t len, uint64_t query, uint64_t target, uint64_t *all, uint64_t *found)
{
    if(!ug_index) return 0;
    uint64_t i, j, l = 0, skip, *pos_list = NULL, cnt, uID, is_q, k_mer = ug_index->k;
    uint64_t x[4], mask = (1ULL<<k_mer) - 1, shift = k_mer - 1, hash;
    (*all) = (*found) = 0;
    for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i) {
        int c = seq_nt4_table[(uint8_t)r[i]];
        ///c = 00, 01, 10, 11
        if (c < 4) { // not an "N" base
            ///x[0] & x[1] are the forward k-mer
            ///x[2] & x[3] are the reverse complementary k-mer
            x[0] = (x[0] << 1 | (c&1))  & mask;
            x[1] = (x[1] << 1 | (c>>1)) & mask;
            x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
            x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
            if (++l >= k_mer)
            {
                hash = hc_hash_long(x, &skip, k_mer);
                if(skip == (uint64_t)-1) continue;
                cnt = get_hc_pt1_count((ha_ug_index*)ug_index, hash, &pos_list);
                if(cnt <= 0) continue;
                
                for (j = 0, is_q = 0; j < cnt; j++)
                {
                    uID = (pos_list[j] << 1) >> (64 - ug_index->uID_bits);
                    if(query == uID) is_q = 1;
                    if(target == uID) break;
                }

                if(is_q)
                {
                    (*found)++;
                    if(j < cnt) (*all)++;
                }
            }
            
        } else l = 0, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
    }

    return 1;
}


typedef struct{
    //[uID_start, uID_end)
    uint64_t uID_start;
    uint64_t uID_end;
    uint64_t u_n;
    uint64_t r_n;
    uint64_t* r_idx;
} bench_utg;

typedef struct{
    uint64_t s, e;
}homo_interval;

typedef struct{
    kvec_t(bench_utg) ug_idx; 
    uint64_t uID_bits;
    uint64_t pos_mode;
    hc_links link;
    kvec_t(homo_interval) regions;
}bench_idx;

uint64_t* set_bench_idx(ma_ug_t *ug, asg_t* read_g, uint64_t uID_start, uint64_t uID_end, uint64_t uID_bits, uint64_t r_n)
{
    uint64_t *idx = (uint64_t*)malloc(sizeof(uint64_t)*r_n), i, k;
    memset(idx, -1, sizeof(uint64_t)*r_n);
    uint64_t rId, ori, start, l;
    ma_utg_t *u = NULL;
    for (i = uID_start; i < uID_end; i++)
    {
        u = &(ug->u.a[i]);
        if(u->n == 0) continue;
        for (k = l = 0; k < u->n; k++)
        {
            rId = u->a[k]>>33;
            ori = u->a[k]>>32&1;
            start = l;
			l += (uint32_t)u->a[k];
            if(idx[rId] != (uint64_t)-1)
            {
                idx[rId] = (uint64_t)-1;
            }
            else
            {
                idx[rId] = (ori<<63) + ((i<<(64-uID_bits))>>1) + start;
                if(ori) idx[rId] = idx[rId] + read_g->seq[rId].len - 1;
            }
        }
    }

    return idx;
}

void get_r_utg_bench(uint64_t index, bench_idx* idx, ma_ug_t *ug)
{

    bench_utg* a_list = idx->ug_idx.a;
    uint64_t a_n = idx->ug_idx.n;
    bench_utg *x = &(a_list[index]), *y = NULL;
    uint64_t i, k, t, rev, x_uid, y_uid, y_pos, x_pos, d;
    uint64_t rId, ori;
    ma_utg_t *u = NULL;
    for (i = x->uID_start; i < x->uID_end; i++)
    {
        u = &(ug->u.a[i]);
        x_uid = i;
        if(u->n == 0) continue;
        for (k = 0; k < u->n; k++)
        {
            rId = u->a[k]>>33;
            ori = u->a[k]>>32&1;
            if(x->r_idx[rId] == (uint64_t)-1) continue;
            x_pos = x->r_idx[rId] & idx->pos_mode;

            for (t = 0; t < a_n; t++)
            {
                if(t == index) continue;
                y = &(a_list[t]);
                if(y->r_idx[rId] == (uint64_t)-1) continue;
                rev = 0;
                if((y->r_idx[rId]>>63) != ori) rev = 1;
                y_uid = (y->r_idx[rId]<<1)>>(64 - idx->uID_bits);
                y_pos = y->r_idx[rId] & idx->pos_mode;
                if(rev) y_pos = ug->u.a[y_uid].len - y_pos - 1;
                ///if(ori) x_pos = ug->u.a[x_uid].len - x_pos - 1, y_pos = ug->u.a[y_uid].len - y_pos - 1;
                d = MAX(x_pos, y_pos) - MIN(x_pos, y_pos);
                d = (d<<2) + (rev<<1);
                if(y_pos > x_pos) d = d + 1;
                push_hc_edge(&(idx->link.a.a[x_uid]), y_uid, 1, 0, &d);
                if(x_pos != y_pos) d = d ^ 1;
                push_hc_edge(&(idx->link.a.a[y_uid]), x_uid, 1, 0, &d);
            }
        }
    }
}

void hap_ID(bench_idx* idx, uint64_t ID, uint64_t* hapID, uint64_t* uID)
{
    uint64_t i;
    (*hapID) = (*uID) = (uint64_t)-1;
    for (i = 0; i < idx->ug_idx.n; i++)
    {
        if(ID >= idx->ug_idx.a[i].uID_start && ID < idx->ug_idx.a[i].uID_end)
        {
            (*hapID) = i;
            (*uID) = ID - idx->ug_idx.a[i].uID_start;
            return;
        }
    }
    return;
}

void print_bench_idx(bench_idx* idx, ma_ug_t *ug)
{
    uint64_t i, k, s_uID, s_hapID, d_uID, d_hapID;
    long long x[2] = {1, -1};
    for (i = 0; i < idx->link.a.n; i++)
    {
        for (k = 0; k < idx->link.a.a[i].e.n; k++)
        {
            if(idx->link.a.a[i].e.a[k].del) continue;
            hap_ID(idx, i, &s_hapID, &s_uID);
            hap_ID(idx, idx->link.a.a[i].e.a[k].uID, &d_hapID, &d_uID);
            fprintf(stderr, "s-hap%lu-utg%.6d\td-hap%lu-utg%.6d\t%c\t%lld\n", 
            s_hapID, (int)(s_uID+1), d_hapID, (int)(d_uID+1), 
            "+-"[!!(idx->link.a.a[i].e.a[k].dis&(uint64_t)2)], 
            ((long long)(idx->link.a.a[i].e.a[k].dis>>2))*x[idx->link.a.a[i].e.a[k].dis&(uint64_t)1]);
        }
    }
    
    
}

uint64_t get_hic_distance_bench_hap(pe_hit_hap* hit, hc_links* link, bench_idx* idx, ma_ug_t *ug, uint64_t* is_trans)
{
    (*is_trans) = (uint64_t)-1;
    uint64_t s_uid, e_uid;
    long long s_pos, e_pos;
    s_uid = ((get_pe_s(*hit)<<1)>>(64 - idx->uID_bits)); s_pos = get_pe_s(*hit) & idx->pos_mode;
    e_uid = ((get_pe_e(*hit)<<1)>>(64 - idx->uID_bits)); e_pos = get_pe_e(*hit) & idx->pos_mode;
    if(s_uid == e_uid)
    {
        (*is_trans) = 0;
        return MAX(s_pos, e_pos) - MIN(s_pos, e_pos);
    }

    uint64_t s_i, e_i, k, ori;
    for (s_i = 0; s_i < idx->ug_idx.n; s_i++)
    {
        if(s_uid >= idx->ug_idx.a[s_i].uID_start && s_uid < idx->ug_idx.a[s_i].uID_end) break;
    }
    for (e_i = 0; e_i < idx->ug_idx.n; e_i++)
    {
        if(e_uid >= idx->ug_idx.a[e_i].uID_start && e_uid < idx->ug_idx.a[e_i].uID_end) break;
    }
    if(s_i == idx->ug_idx.n || e_i == idx->ug_idx.n) return (uint64_t)-1;
    if(s_i == e_i)
    {
        (*is_trans) = 0;
        return (uint64_t)-1;
    }

    (*is_trans) = 1;
    hc_linkeage* t = &(link->a.a[s_uid]);
    long long m_x[2] = {1, -1}, dis;
    for (k = 0; k < t->e.n; k++)
    {
        if(t->e.a[k].del || t->e.a[k].uID != e_uid) continue;
        ori = !!(t->e.a[k].dis & (uint64_t)2);
        dis = (long long)(t->e.a[k].dis>>2) * m_x[t->e.a[k].dis & (uint64_t)1];
        if(ori) e_pos = ug->u.a[e_uid].len - e_pos - 1;
        e_pos = e_pos + dis;
        return MAX(s_pos, e_pos) - MIN(s_pos, e_pos);
    }

    return (uint64_t)-1;
}

uint64_t get_hic_distance_bench(pe_hit* hit, hc_links* link, bench_idx* idx, ma_ug_t *ug, uint64_t* is_trans)
{
    (*is_trans) = (uint64_t)-1;
    uint64_t s_uid, e_uid;
    long long s_pos, e_pos;
    s_uid = ((hit->s<<1)>>(64 - idx->uID_bits)); s_pos = hit->s & idx->pos_mode;
    e_uid = ((hit->e<<1)>>(64 - idx->uID_bits)); e_pos = hit->e & idx->pos_mode;
    if(s_uid == e_uid)
    {
        (*is_trans) = 0;
        return MAX(s_pos, e_pos) - MIN(s_pos, e_pos);
    }

    uint64_t s_i, e_i, k, ori;
    for (s_i = 0; s_i < idx->ug_idx.n; s_i++)
    {
        if(s_uid >= idx->ug_idx.a[s_i].uID_start && s_uid < idx->ug_idx.a[s_i].uID_end) break;
    }
    for (e_i = 0; e_i < idx->ug_idx.n; e_i++)
    {
        if(e_uid >= idx->ug_idx.a[e_i].uID_start && e_uid < idx->ug_idx.a[e_i].uID_end) break;
    }
    if(s_i == idx->ug_idx.n || e_i == idx->ug_idx.n) return (uint64_t)-1;
    if(s_i == e_i)
    {
        (*is_trans) = 0;
        return (uint64_t)-1;
    }

    (*is_trans) = 1;
    hc_linkeage* t = &(link->a.a[s_uid]);
    long long m_x[2] = {1, -1}, dis;
    for (k = 0; k < t->e.n; k++)
    {
        if(t->e.a[k].del || t->e.a[k].uID != e_uid) continue;
        ori = !!(t->e.a[k].dis & (uint64_t)2);
        dis = (long long)(t->e.a[k].dis>>2) * m_x[t->e.a[k].dis & (uint64_t)1];
        if(ori) e_pos = ug->u.a[e_uid].len - e_pos - 1;
        e_pos = e_pos + dis;
        return MAX(s_pos, e_pos) - MIN(s_pos, e_pos);
    }

    return (uint64_t)-1;
}


void init_bench_idx(bench_idx* idx, asg_t* read_g, ma_ug_t *ug)
{
    uint64_t i, occ;
    kv_init(idx->ug_idx);
    kv_init(idx->regions);
    kv_malloc(idx->ug_idx, ug->occ.n); idx->ug_idx.n = ug->occ.n;
    for (idx->uID_bits = 1; (uint64_t)(1<<idx->uID_bits)<(uint64_t)ug->u.n; idx->uID_bits++);
    idx->pos_mode = ((uint64_t)-1)>>(idx->uID_bits+1);
    for (i = occ = 0; i < ug->occ.n; i++)
    {
        idx->ug_idx.a[i].uID_start = occ;
        occ += ug->occ.a[i];
        idx->ug_idx.a[i].uID_end = occ;
        idx->ug_idx.a[i].u_n = ug->occ.a[i];

        idx->ug_idx.a[i].r_n = read_g->n_seq;
        idx->ug_idx.a[i].r_idx 
            = set_bench_idx(ug, read_g, idx->ug_idx.a[i].uID_start, idx->ug_idx.a[i].uID_end, 
                                                idx->uID_bits, idx->ug_idx.a[i].r_n);
    }

    init_hc_links(&(idx->link), ug->u.n, NULL);

    for (i = 0; i < idx->ug_idx.n; i++)
    {
        get_r_utg_bench(i, idx, ug);
    }
}

void evaluate_bench_idx_hap(bench_idx* idx, kvec_pe_hit_hap* hits, ma_ug_t *ug)
{
    uint64_t k, distance, is_trans, trans[2];
    kvec_t(uint64_t) buf;
    kv_init(buf);
    for (k = trans[0] = trans[1] = 0; k < hits->n_u; ++k) 
    {
        distance = get_hic_distance_bench_hap(&(hits->a[k]), &(idx->link), idx, ug, &is_trans);
        if(is_trans != (uint64_t)-1) trans[is_trans]++;
        if(distance == (uint64_t)-1 || is_trans == (uint64_t)-1) continue;
        distance = (distance << 1) + is_trans;
        kv_push(uint64_t, buf, distance); 
    }

    radix_sort_hc64(buf.a, buf.a+buf.n);

    for (k = 0; k < buf.n; k++)
    {
        fprintf(stderr, "%lu\t%lu\n", buf.a[k]>>1, buf.a[k]&1);
    }
    /**
    uint64_t up_dis = buf.a[(uint64_t)(buf.n*0.99)]>>1, step = 7240;
    uint64_t step_s = 0, step_e = step, cnt[2];
    for (k = cnt[0] = cnt[1] = 0; k < buf.n; k++)
    {
        if(step_s > up_dis) step_e = (buf.a[buf.n-1]>>1) + 1;
        if((buf.a[k]>>1) < step_e && (buf.a[k]>>1) >= step_s)
        {
            cnt[buf.a[k]&1]++;
        }
        if((buf.a[k]>>1) >= step_e)
        {
            while (!((buf.a[k]>>1) < step_e && (buf.a[k]>>1) >= step_s))
            {
                fprintf(stderr, "i: %lu, step_s: %lu, step_e: %lu, cnt[0]: %lu, cnt[1]: %lu, rate: %f\n", 
                step_s/step, step_s, step_e, cnt[0], cnt[1], ((double)cnt[1])/(double)(cnt[0] + cnt[1]));
                step_s += step;
                step_e += step;
                cnt[0] = cnt[1] = 0;
            }
        }
    }

    if(cnt[0] > 0 || cnt[1] > 0)
    {
        fprintf(stderr, "i: %lu, step_s: %lu, step_e: %lu, cnt[0]: %lu, cnt[1]: %lu, rate: %f\n", 
                step_s/step, step_s, step_e, cnt[0], cnt[1], ((double)cnt[1])/(double)(cnt[0] + cnt[1]));
    }
    **/
    kv_destroy(buf);
}


void evaluate_bench_idx(bench_idx* idx, kvec_pe_hit* hits, ma_ug_t *ug)
{
    uint64_t k, distance, is_trans, trans[2];
    kvec_t(uint64_t) buf;
    kv_init(buf);
    for (k = trans[0] = trans[1] = 0; k < hits->a.n; ++k) 
    {
        distance = get_hic_distance_bench(&(hits->a.a[k]), &(idx->link), idx, ug, &is_trans);
        if(is_trans != (uint64_t)-1) trans[is_trans]++;
        if(distance == (uint64_t)-1 || is_trans == (uint64_t)-1) continue;
        distance = (distance << 1) + is_trans;
        kv_push(uint64_t, buf, distance); 
    }

    radix_sort_hc64(buf.a, buf.a+buf.n);

    for (k = 0; k < buf.n; k++)
    {
        fprintf(stderr, "%lu\t%lu\n", buf.a[k]>>1, buf.a[k]&1);
    }
    /**
    uint64_t up_dis = buf.a[(uint64_t)(buf.n*0.99)]>>1, step = 7240;
    uint64_t step_s = 0, step_e = step, cnt[2];
    for (k = cnt[0] = cnt[1] = 0; k < buf.n; k++)
    {
        if(step_s > up_dis) step_e = (buf.a[buf.n-1]>>1) + 1;
        if((buf.a[k]>>1) < step_e && (buf.a[k]>>1) >= step_s)
        {
            cnt[buf.a[k]&1]++;
        }
        if((buf.a[k]>>1) >= step_e)
        {
            while (!((buf.a[k]>>1) < step_e && (buf.a[k]>>1) >= step_s))
            {
                fprintf(stderr, "i: %lu, step_s: %lu, step_e: %lu, cnt[0]: %lu, cnt[1]: %lu, rate: %f\n", 
                step_s/step, step_s, step_e, cnt[0], cnt[1], ((double)cnt[1])/(double)(cnt[0] + cnt[1]));
                step_s += step;
                step_e += step;
                cnt[0] = cnt[1] = 0;
            }
        }
    }

    if(cnt[0] > 0 || cnt[1] > 0)
    {
        fprintf(stderr, "i: %lu, step_s: %lu, step_e: %lu, cnt[0]: %lu, cnt[1]: %lu, rate: %f\n", 
                step_s/step, step_s, step_e, cnt[0], cnt[1], ((double)cnt[1])/(double)(cnt[0] + cnt[1]));
    }
    **/
    kv_destroy(buf);
}

void destory_bench_idx(bench_idx* idx)
{
    uint64_t i;
    for (i = 0; i < idx->ug_idx.n; i++)
    {
        free(idx->ug_idx.a[i].r_idx);
    }
    kv_destroy(idx->ug_idx);
    kv_destroy(idx->regions);
    destory_hc_links(&(idx->link));
}


int hic_short_align_bench(const enzyme *fn1, const enzyme *fn2, const char *output_file_name, ha_ug_index* idx)
{
    double index_time = yak_realtime();
    sldat_t sl;
    sl.idx = idx;
    ///sl.link = NULL;
    sl.chunk_size = 20000000;
    sl.n_thread = asm_opt.thread_num;
    sl.total_base = sl.total_pair = 0;
    idx->hap_cnt = asm_opt.hap_occ;
    kv_init(sl.hits.a);
    fprintf(stderr, "u.n: %d, uID_bits: %lu, pos_bits: %lu\n", (uint32_t)idx->ug->u.n, idx->uID_bits, idx->pos_bits);
    
    if(!load_hc_hits(&sl.hits, idx->ug, output_file_name))
    {
        // kt_pipeline(3, worker_pipeline, &sl, 3);
        // dedup_hits(&sl.hits);
        alignment_worker_pipeline(&sl, fn1, fn2);
        write_hc_hits(&sl.hits, idx->ug, output_file_name);
    }
    bench_idx bench;
    init_bench_idx(&bench, idx->read_g, idx->ug);
    ///print_bench_idx(&bench, idx->ug);
    evaluate_bench_idx(&bench, &sl.hits, idx->ug);

    destory_bench_idx(&bench);
    kv_destroy(sl.hits.a);
    fprintf(stderr, "[M::%s::%.3f] processed %lu pairs; %lu bases\n", __func__, yak_realtime()-index_time, sl.total_pair, sl.total_base);
    return 1;
}

void hic_benchmark(ma_ug_t *ug, asg_t* read_g)
{
    char *output_file_name = (char*)calloc(strlen(asm_opt.output_file_name) + 25, 1);
	sprintf(output_file_name, "%s.bench", asm_opt.output_file_name);
    ug_index = NULL;
    int exist = load_hc_pt_index(&ug_index, ug, output_file_name);
    if(exist == 0) ug_index = build_unitig_index(ug, asm_opt.hic_mer_length, asm_opt.hap_occ, 0, asm_opt.thread_num);
    if(exist == 0) write_hc_pt_index(ug_index, output_file_name);
    ug_index->ug = ug;
    ug_index->read_g = read_g;

    hic_short_align_bench(asm_opt.hic_reads[0], asm_opt.hic_reads[1], output_file_name, ug_index);

    free(output_file_name);
}