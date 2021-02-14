#define __STDC_LIMIT_MACROS
#include "float.h"
#include <math.h>
#include "hic.h"
#include "htab.h"
#include "assert.h"
#include "Overlaps.h"
#include "Hash_Table.h"
#include "Correct.h"
#include "khashl.h"
#include "kthread.h"
#include "ksort.h"
#include "kseq.h" // FASTA/Q parser
#include "kdq.h"
KSEQ_INIT(gzFile, gzread)
KDQ_INIT(uint64_t)


#define OFFSET_RATE 0.000000001
#define OFFSET_SECOND_RATE 0.0000000001
#define SCALL 10000
#define OFFSET_RATE_MAX_W 20.8286263517*SCALL
#define OFFSET_RATE_MIN_W 4.0000003e-10*SCALL

#define HIC_COUNTER_BITS 12
#define HIC_MAX_COUNT    ((1<<HIC_COUNTER_BITS)-1)
#define HIC_KEY_MODE ((uint64_t)(((uint64_t)-1)-HIC_MAX_COUNT))
#define HIC_R_E_RATE 0.01

const unsigned char b2rc[5] = {'T', 'G', 'C', 'A', 'N'};

#define hic_ct_eq(a, b) ((a)>>HIC_COUNTER_BITS == (b)>>HIC_COUNTER_BITS) 
#define hic_ct_hash(a) ((a)>>HIC_COUNTER_BITS)
KHASHL_MAP_INIT(static klib_unused, hc_pt_t, hc_pt, uint64_t, uint64_t, hic_ct_hash, hic_ct_eq)

typedef struct{
    kvec_t(char) name;
    kvec_t(uint64_t) name_Len;
    kvec_t(char) r;
    kvec_t(uint64_t) r_Len;
} reads_t;

typedef struct{
    kvec_t(uint8_t) vis;
    kvec_t(uint64_t) x;
    kvec_t(uint64_t) dis;
    uint64_t uID_mode, uID_shift, tmp_v, tmp_d;
}pdq;

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
    hc_links* link;
    uint64_t uID_bits;
    uint64_t uID_mode;
    uint64_t pos_bits;
    uint64_t pos_mode;
    uint64_t rev_mode;
    uint64_t k;
    uint64_t max_cnt;


    uint64_t pre;
    uint64_t tot;
    uint64_t tot_pos;
    uint64_t up_bound;
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
	uint64_t s, e, id;
} pe_hit;

typedef struct {
    kvec_t(pe_hit) a;
} kvec_pe_hit;

typedef struct {
    kvec_t(hc_edge) a;
}kvec_hc_edge;


#define pe_hit_an1_key(x) ((x).s)
KRADIX_SORT_INIT(pe_hit_an1, pe_hit, pe_hit_an1_key, 8)
#define pe_hit_an2_key(x) ((x).e)
KRADIX_SORT_INIT(pe_hit_an2, pe_hit, pe_hit_an2_key, 8)
#define generic_key(x) (x)
KRADIX_SORT_INIT(hc64, uint64_t, generic_key, 8)
KRADIX_SORT_INIT(u32, uint32_t, generic_key, 4)
#define g_partition_key(x) (((x)>>1)+((x)<<63))
KRADIX_SORT_INIT(g_partition, uint64_t, g_partition_key, 8)


typedef struct { // global data structure for kt_pipeline()
	const ha_ug_index* idx;
	kseq_t *ks1, *ks2;
    int64_t chunk_size;
    uint64_t n_thread;
    uint64_t total_base;
    uint64_t total_pair;
    kvec_pe_hit hits;
    hc_links* link;
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
    hc_links* link;
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
#define hc_edge_key_u(a) ((a).uID)
KRADIX_SORT_INIT(hc_edge_u, hc_edge, hc_edge_key_u, 4)
#define hc_edge_key_d(a) ((a).dis)
KRADIX_SORT_INIT(hc_edge_d, hc_edge, hc_edge_key_d, member_size(hc_edge, dis))

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
    uint64_t max;
} trans_idx;


reads_t R1, R2;
ha_ug_index* ug_index;

void build_bub_graph(ma_ug_t* ug, bubble_type* bub);

void init_ha_ug_index_opt(ha_ug_index* idx, ma_ug_t *ug, int k, pldat_t* p)
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
    idx->up_bound = 1;
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
    p->n_thread = asm_opt.thread_num;
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

void hc_pt_t_gen_single(hc_pt1_t* pt, uint64_t* up_bound)
{
    khint_t k;
    uint64_t c;

    if(up_bound)
    {
        for (k = 0; k != kh_end(pt->h); ++k) {
            if (kh_exist(pt->h, k)) {
                if(kh_val(pt->h, k) > (*up_bound))
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

int write_hc_pt_index(ha_ug_index* idx, char* file_name)
{
    char* gfa_name = (char*)malloc(strlen(file_name)+25);
    sprintf(gfa_name, "%s.hic.tlb.bin", file_name);
    FILE* fp = fopen(gfa_name, "w");
    if (!fp) {
        free(gfa_name);
        return 0;
    }

    fwrite(&idx->uID_bits, sizeof(idx->uID_bits), 1, fp);
    fwrite(&idx->uID_mode, sizeof(idx->uID_mode), 1, fp);
    fwrite(&idx->pos_bits, sizeof(idx->pos_bits), 1, fp);
    fwrite(&idx->pos_mode, sizeof(idx->pos_mode), 1, fp);
    fwrite(&idx->rev_mode, sizeof(idx->rev_mode), 1, fp);
    fwrite(&idx->k, sizeof(idx->k), 1, fp);
    fwrite(&idx->pre, sizeof(idx->pre), 1, fp);
    fwrite(&idx->tot, sizeof(idx->tot), 1, fp);
    fwrite(&idx->tot_pos, sizeof(idx->tot_pos), 1, fp);

    uint64_t i = 0;
    for (i = 0; i < idx->tot; i++)
    {
        fwrite(&idx->idx_buf[i].n, sizeof(idx->idx_buf[i].n), 1, fp);
        fwrite(&idx->idx_buf[i].end, sizeof(idx->idx_buf[i].end), 1, fp);
        fwrite(idx->idx_buf[i].a, sizeof(uint64_t), idx->idx_buf[i].n, fp);
        hc_pt_save(idx->idx_buf[i].h, fp);
    }

    fprintf(stderr, "[M::%s] Index has been written.\n", __func__);
    free(gfa_name);
    fclose(fp);
    return 1;
}

int load_hc_pt_index(ha_ug_index** r_idx, char* file_name)
{
    uint64_t flag = 0;
    double index_time = yak_realtime();
    char* gfa_name = (char*)malloc(strlen(file_name)+25);
    sprintf(gfa_name, "%s.hic.tlb.bin", file_name);
    FILE* fp = fopen(gfa_name, "r");
    if (!fp) {
        free(gfa_name);
        return 0;
    }
    ha_ug_index* idx = NULL; CALLOC(idx, 1);

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
    uint64_t i = 0;
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
    fclose(fp);
    fprintf(stderr, "[M::%s::%.3f] ==> HiC index has been loaded\n", __func__, yak_realtime()-index_time);
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
            hc_pt_t_gen_single(&(idx->idx_buf[i]), &(idx->up_bound));
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

ha_ug_index* build_unitig_index(ma_ug_t *ug, int k)
{
    ha_ug_index* idx = NULL; CALLOC(idx, 1);
    pldat_t pl; pl.h = idx; pl.is_cnt = 1;
    double index_time = yak_realtime(), beg_time;
    init_ha_ug_index_opt(idx, ug, k, &pl);

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
    (*exact_len) = p->off_cnt >> 32;
    if(total_len != NULL)
    {
        (*exact_len) = (p->off_cnt>>32) & ((uint64_t)65535);
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
    if(n == 1) return (a[0].off_cnt>>32);
    long long i = 0;
    uint64_t cur_beg, cur_end, beg, end, ovlp = 0, tLen = 0;
    cur_end = (uint32_t)a[n-1].off_cnt;
    cur_beg = cur_end + 1 - (a[n-1].off_cnt>>32);

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

void get_alignment(char *r, uint64_t len, uint64_t k_mer, kvec_vote* buf, 
const ha_ug_index* idx, uint64_t buf_iter, uint64_t rid)
{
    uint64_t i, j, l = 0, skip, *pos_list = NULL, cnt, rev, self_p, ref_p, u_len, uID;
    uint64_t x[4], mask = (1ULL<<k_mer) - 1, shift = k_mer - 1, hash;
    s_hit *p = NULL;
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
                if(cnt > idx->max_cnt) continue;
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
    
    ///if(buf->a.n - buf_iter <= 1) return;
    if(buf->a.n - buf_iter == 0) return;
    if(buf->a.n - buf_iter > 1) radix_sort_hc_s_hit_an1(buf->a.a + buf_iter, buf->a.a + buf->a.n);



    /*******************************for debug************************************/
    // print_pos_list(idx, buf->a.a+buf_iter, buf->a.n - buf_iter, rid, (buf_iter != 0));
    // fprintf(stderr, "len0:%lu\n", buf->a.n - buf_iter);
    // for (i = buf_iter; i < buf->a.n; i++)
    // {
    //     interpret_pos(idx, &buf->a.a[i], &rev, &uID, &ref_p, &self_p, &cnt, NULL);
    //     fprintf(stderr, "(%lu) rev: %lu, uID: %lu, ref_p: %lu, self_p: %lu, len: %lu\n", 
    //                                                         i, rev, uID, ref_p, self_p, cnt);
    // }
    /*******************************for debug************************************/





    uint64_t cur_ref_p, thres = (len * HIC_R_E_RATE) + 1, m, index_beg, ovlp, maxLen = 0, max_i = (uint64_t)-1;
    i = m = buf_iter;
    while (i < buf->a.n)
    {
        interpret_pos(idx, &buf->a.a[i], &rev, &uID, &ref_p, &self_p, &cnt, NULL);
        /*******************************for debug************************************/
        // if(check_exact_match(r, self_p, len, idx->ug->u.a[uID].s, 
        //                         ref_p, idx->ug->u.a[uID].len, cnt, rev, 1) != cnt)
        // {
        //     fprintf(stderr, "ERROR\n");
        // }
        /*******************************for debug************************************/
        // if(self_p > ref_p)
        // {
        //     i++; 
        //     continue; ///fix this in future
        // } 
        cur_ref_p = buf->a.a[i].ref;
        index_beg = i;
        while ((i < buf->a.n) && ((buf->a.a[i].ref>>idx->pos_bits) == (cur_ref_p>>idx->pos_bits)) && 
               (buf->a.a[i].ref - cur_ref_p <= thres))
        {
            i++;
        }
        if(i - index_beg > 1)
        {
            radix_sort_hc_s_hit_an2(buf->a.a + index_beg, buf->a.a + i);
        }
        ovlp = collect_votes(buf->a.a + index_beg, i - index_beg);
        buf->a.a[m] = buf->a.a[i - 1];
        buf->a.a[m].off_cnt = (buf->a.a[m].off_cnt << 32)>>32;
        buf->a.a[m].off_cnt += ((uint64_t)ovlp<<32);

        if(maxLen < (ovlp&((uint64_t)65535))) maxLen = (ovlp&((uint64_t)65535)), max_i = m;

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

    compress_mapped_pos(idx, buf, buf_iter, max_i, thres);
    
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
inline int is_unreliable_hits(long long rev, long long ref_p, long long tLen, uint64_t uID, hc_links* link)
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

    p = &(link->bed.a[uID]);
    for (i = 0; i < p->n; i++)
    {
        if(inter_interval(p_beg, p_end, p->a[i].beg, p->a[i].end, NULL, NULL)) break;
    }
    if(p->n > 0 && i < p->n) return 1;

    return 0;
}
inline void set_pe_pos(ha_ug_index* idx, s_hit *l1, uint64_t occ1, s_hit *l2, uint64_t occ2, 
pe_hit* x, uint64_t rid, hc_links* link)
{
    if(occ1 == 0 || occ2 == 0) return;
    uint64_t rev1, rev2, uID1, uID2, ref_p1, ref_p2, self_p1, self_p2, eLen1, eLen2, tLen1, tLen2;
    uint64_t rev_t, uID_t, ref_p_t, self_p_t, eLen_t, tLen_t;
    ///5' end of r1 and r2
    interpret_pos(idx, &l1[0], &rev1, &uID1, &ref_p1, &self_p1, &eLen1, &tLen1);
    interpret_pos(idx, &l2[0], &rev2, &uID2, &ref_p2, &self_p2, &eLen2, &tLen2);
    /*******************************for debug************************************/
    // if(rid == 33045391 || rid == 4239289 || rid == 5267597 || rid == 34474764 || rid == 35016489
    //         || rid == 36002255 || rid == 37811694 || rid == 46805824)
    // {
    //     fprintf(stderr, "rid: %lu, rev1: %lu, uID1: %lu, ref_p1: %lu, self_p1: %lu, rev2: %lu, uID2: %lu, ref_p2: %lu, self_p2: %lu\n", 
    //     rid, rev1, uID1, self_p1, ref_p1, rev2, uID2, ref_p2, self_p2);
    // }
    /*******************************for debug************************************/
    ///if(uID1 == uID2) return;
    if(ref_p1 < self_p1 || ref_p2 < self_p2) return;
    if(occ1 > 1)
    {
        interpret_pos(idx, &l1[1], &rev_t, &uID_t, &ref_p_t, &self_p_t, &eLen_t, &tLen_t);
        if(uID_t != uID1 && uID_t != uID2) return;
    }

    if(occ2 > 1)
    {
        interpret_pos(idx, &l2[1], &rev_t, &uID_t, &ref_p_t, &self_p_t, &eLen_t, &tLen_t);
        if(uID_t != uID1 && uID_t != uID2) return;
    }
    
    x->id = rid; 
    ref_p1 -= self_p1; 
    if(rev1) ref_p1 = idx->ug->u.a[uID1].len - 1 - ref_p1;
    x->s = (rev1<<63) | ((uID1 << (64-idx->uID_bits))>>1) | (ref_p1 & idx->pos_mode);
    
    
    ref_p2 -= self_p2; 
    if(rev2) ref_p2 = idx->ug->u.a[uID2].len - 1 - ref_p2;
    x->e = (rev2<<63) | ((uID2 << (64-idx->uID_bits))>>1) | (ref_p2 & idx->pos_mode); 

    if(link && (is_unreliable_hits(rev1, ref_p1, tLen1, uID1, link) || 
                    is_unreliable_hits(rev2, ref_p2, tLen2, uID2, link)))
    {
        x->id = x->s = x->e = (uint64_t)-1;
    }
}

static void worker_for_alignment(void *data, long i, int tid) // callback for kt_for()
{
    stepdat_t *s = (stepdat_t*)data;
    s->pos[i].id = s->pos[i].s = s->pos[i].e = (uint64_t)-1;
    uint64_t len1 = s->len[i]>>32, len2 = (uint32_t)s->len[i], occ1, occ2;
    char *r1 = s->seq[i], *r2 = s->seq[i] + len1;
    /*******************************for debug************************************/
    // if(memcmp(r1, R1.r.a + R1.r_Len.a[s->id+i], len1) != 0)
    // {
    //     fprintf(stderr, "haha1\n");
    // }
    // if(memcmp(r2, R2.r.a + R2.r_Len.a[s->id+i], len2) != 0)
    // {
    //     fprintf(stderr, "haha2\n");
    // }
    /*******************************for debug************************************/
    s->pos_buf[tid].a.n = 0;
    get_alignment(r1, len1, s->idx->k, &s->pos_buf[tid], s->idx, 0, s->id+i);
    occ1 = s->pos_buf[tid].a.n;
    if(occ1 == 0) return;
    get_alignment(r2, len2, s->idx->k, &s->pos_buf[tid], s->idx, occ1, s->id+i);
    occ2 = s->pos_buf[tid].a.n - occ1;
    if(occ2 == 0) return;   
    
    set_pe_pos((ha_ug_index*)s->idx, s->pos_buf[tid].a.a, occ1, s->pos_buf[tid].a.a + occ1, occ2, &(s->pos[i]), s->id+i, s->link);

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
        s->idx = p->idx; s->id = p->total_pair; s->link = p->link;
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
        MALLOC(s->pos, s->n);
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
            if(s->pos[i].s == (uint64_t)-1) continue;
            kv_push(pe_hit, p->hits.a, s->pos[i]);
        }
        free(s->pos);
        free(s);
    }
    return 0;
}


void load_reads(reads_t* x, const char *fn)
{
    kv_init(x->name);
    kv_init(x->name_Len);
    kv_init(x->r);
    kv_init(x->r_Len);
    gzFile fp;
    kseq_t *ks;
    int ret;
    uint64_t name_tot, base_total;

    name_tot = base_total = 0;
    if ((fp = gzopen(fn, "r")) == 0) return;
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

    kv_push(uint64_t, x->name_Len, name_tot);
    kv_push(uint64_t, x->r_Len, base_total);

    kseq_destroy(ks);
    gzclose(fp);
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

void print_hits(ha_ug_index* idx, kvec_pe_hit* hits, const char *fn)
{
    uint64_t k, shif = 64 - idx->uID_bits;
    reads_t r1;
    load_reads(&r1, fn);
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

void dedup_hits(kvec_pe_hit* hits)
{
    double index_time = yak_realtime();
    uint64_t k, l, m = 0, cur;
    radix_sort_pe_hit_an1(hits->a.a, hits->a.a + hits->a.n);
    for (k = 1, l = 0; k <= hits->a.n; ++k) 
    {   
        if (k == hits->a.n || hits->a.a[k].s != hits->a.a[l].s) 
        {
            if (k - l > 1) radix_sort_pe_hit_an2(hits->a.a + l, hits->a.a + k);
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
            l = k;
        }
    }
    hits->a.n = m;
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

void update_bub_b_s_idx(bubble_type* bub);
void identify_bubbles(ma_ug_t* ug, bubble_type* bub, hc_links* link)
{
    asg_cleanup(ug->g);
    if (!ug->g->is_symm) asg_symm(ug->g);
    uint32_t v, n_vtx = ug->g->n_seq * 2, i, k, mode = (((uint32_t)-1)<<2);
    uint32_t beg, sink, n, *a, n_occ;
    uint64_t pathLen, tLen;
    bub->ug = ug; 
    for (i = 0, tLen = 1; i < ug->u.n; i++) tLen += ug->u.a[i].len;
    bub->b_bub = bub->b_end_bub = bub->tangle_bub = bub->cross_bub = bub->mess_bub = 0;

    if(bub->round_id == 0)
    {
        buf_t b; memset(&b, 0, sizeof(buf_t)); b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
        kv_init(bub->list); kv_init(bub->num); kv_init(bub->pathLen);
        kv_init(bub->b_s_idx); kv_malloc(bub->b_s_idx, ug->g->n_seq); 
        bub->b_ug = NULL; kv_init(bub->chain_weight);
        bub->b_s_idx.n = ug->g->n_seq;
        memset(bub->b_s_idx.a, -1, bub->b_s_idx.n * sizeof(uint64_t));

        CALLOC(bub->index, n_vtx);
        for (i = 0; i < ug->g->n_seq; i++)
        {
            if(ug->g->seq[i].c > 0)
            {
                bub->index[i] = (ug->g->seq[i].c << 2);
                ug->g->seq[i].c = 0;
            } 
        }

        

        for (v = 0; v < n_vtx; ++v) 
        {
            if(ug->g->seq[v>>1].del) continue;
            if(asg_arc_n(ug->g, v) < 2) continue;
            if((bub->index[v]&(uint32_t)3) != 0) continue;
            if(asg_bub_pop1_primary_trio(ug->g, NULL, v, tLen, &b, (uint32_t)-1, (uint32_t)-1, 0, NULL, NULL))
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
            if(asg_bub_pop1_primary_trio(ug->g, NULL, v, tLen, &b, (uint32_t)-1, (uint32_t)-1, 0, &pathLen, NULL))
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
            if(asg_bub_pop1_primary_trio(ug->g, NULL, v, tLen, &b, (uint32_t)-1, (uint32_t)-1, 0, &pathLen, NULL))
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
            if((bub->index[i]>>2) == 0)
            {
                bub->index[i] = (uint32_t)-1; 
            }
            else
            {
                if((bub->index[i]>>2) == 1)
                {
                    bub->index[i] = P_het(*bub); ///potential het
                }
                else
                {
                    bub->index[i] = M_het(*bub); ///must het
                }
            }         
        }


        for (i = 0; i < bub->f_bub; i++)
        {
            get_bubbles(bub, i, &beg, &sink, &a, &n, &pathLen);
            for (v = n_occ = 0; v < n; v++)
            {
                bub->index[(a[v]>>1)] = i;
                n_occ += ug->u.a[a[v]>>1].n;
            }

            if((pathLen*2) >= ug->g->seq[beg>>1].len && (pathLen*2) >= ug->g->seq[sink>>1].len)
            {
                bub->index[(beg>>1)] = (uint32_t)-1;
                bub->index[(sink>>1)] = (uint32_t)-1;
            }

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
            if(bub->index[i] > P_het(*bub))
            {
                if(link)
                {
                    for (k = 0; k < link->a.a[i].f.n; k++)
                    {
                        if(link->a.a[i].f.a[k].del || link->a.a[i].f.a[k].dis != RC_1) continue;
                        bub->index[i] = P_het(*bub);
                        break;
                    }
                }
            }
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





void push_hc_edge(hc_linkeage* x, uint64_t uID, double weight, int dir, uint64_t* d)
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
            return;
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



void all_pair_shortest_path(const ha_ug_index* idx, hc_links* link, MT* M)
{
    double index_time = yak_realtime();
    asg_t *sg = idx->ug->g;
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
    fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
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
    const ha_ug_index* idx; 
    hc_links* link; 
    MT* M; 
    bubble_type* bub;
    uint8_t** dis_buf;
} utg_d_t;


static void worker_for_dis(void *data, long i, int tid) 
{
    utg_d_t* s = (utg_d_t*)data;
    const ha_ug_index* idx = s->idx;
    hc_links* link = s->link;
    MT* M = s->M; 
    bubble_type* bub = s->bub;
    uint8_t* dis_buf = s->dis_buf[tid];
    asg_t *sg = idx->ug->g;
    hc_linkeage* t = NULL;
    uint32_t n_vtx = sg->n_seq<<1, v, u, k, j;
    uint64_t d[2], db[2], q_u, min, min_i, min_b, rev[2], min_rev;

    if (sg->seq[i].del) return;
    t = &(link->a.a[i]);
    if (t->e.n == 0) return;

    for (k = 0; k < t->e.n; k++)
    {
        if(t->e.a[k].del) continue;
        u = t->e.a[k].uID;

        for (v = ((uint64_t)(i)<<1); v < ((uint64_t)(i+1)<<1); v++)
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
                t->e.a[k].dis += v&1;
                t->e.a[k].dis <<=1;
                t->e.a[k].dis += min_i;
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
                t->e.a[k].dis += ((v&1)^min_rev);
                t->e.a[k].dis <<=1;
                t->e.a[k].dis += (min_i^min_rev);
            }
        }
    }

}

void fill_utg_distance_multi(const ha_ug_index* idx, hc_links* link, MT* M, bubble_type* bub)
{
    double index_time = yak_realtime();
    uint32_t i;
    utg_d_t s;
    s.idx = idx; s.link = link; s.M = M; s.bub = bub;
    s.dis_buf = (uint8_t**)malloc(sizeof(uint8_t*)*asm_opt.thread_num);
    for (i = 0; i < (uint32_t)asm_opt.thread_num; i++)
    {
        s.dis_buf[i] = (uint8_t*)malloc(sizeof(uint8_t)*(s.idx->ug->g->n_seq<<1));
    }
    
    kt_for(asm_opt.thread_num, worker_for_dis, &s, s.idx->ug->g->n_seq);


    for (i = 0; i < (uint32_t)asm_opt.thread_num; i++)
    {
        free(s.dis_buf[i]);
    }
    free(s.dis_buf);
    fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
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

    all_pair_shortest_path(idx, link, M);
    fill_utg_distance_multi(idx, link, M, bub);

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
	sprintf(buf, "%s.hic.link", fn);
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
    fwrite(&link->r_num, sizeof(link->r_num), 1, fp);
    fwrite(link->u_idx, sizeof(uint32_t), 1, fp);
    
    fwrite(&(link->bed.n), sizeof(link->bed.n), 1, fp);
    for (k = 0; k < link->bed.n; k++)
    {
        fwrite(&(link->bed.a[k].n), sizeof(link->bed.a[k].n), 1, fp);
        fwrite(link->bed.a[k].a, sizeof(uint64_t)*link->bed.a[k].n, 1, fp);
    }
    


    fclose(fp);
    free(buf);
    fprintf(stderr, "[M::%s::] ==> Hi-C linkages have been written\n", __func__);
}

int load_hc_links(hc_links* link, const char *fn)
{
    uint64_t k, flag = 0;
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
	sprintf(buf, "%s.hic.link", fn);

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
    fread(&link->r_num, sizeof(link->r_num), 1, fp);
    MALLOC(link->u_idx, link->r_num);
    fread(link->u_idx, sizeof(uint32_t), 1, fp);



    kv_init(link->bed);
    flag += fread(&(link->bed.n), sizeof(link->bed.n), 1, fp);
    link->bed.m = link->bed.n; CALLOC(link->bed.a, link->bed.n);
    for (k = 0; k < link->bed.n; k++)
    {
        flag += fread(&(link->bed.a[k].n), sizeof(link->bed.a[k].n), 1, fp);
        link->bed.a[k].m = link->bed.a[k].n; MALLOC(link->bed.a[k].a, link->bed.a[k].n);
        flag += fread(link->bed.a[k].a, sizeof(uint64_t)*link->bed.a[k].n, 1, fp);
    }



    fclose(fp);
    free(buf);
    fprintf(stderr, "[M::%s::] ==> Hi-C linkages have been loaded\n", __func__);
    return 1;
}


void write_hc_hits(kvec_pe_hit* hits, const char *fn)
{
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
	sprintf(buf, "%s.hic.lk.bin", fn);
    FILE* fp = fopen(buf, "w");

    fwrite(&hits->a.n, sizeof(hits->a.n), 1, fp);
    fwrite(hits->a.a, sizeof(pe_hit), hits->a.n, fp);

    fclose(fp);
    free(buf);
}

int load_hc_hits(kvec_pe_hit* hits, const char *fn)
{
    uint64_t flag = 0;
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
	sprintf(buf, "%s.hic.lk.bin", fn);

    FILE* fp = NULL; 
    fp = fopen(buf, "r"); 
    if(!fp) return 0;

    kv_init(hits->a);
    flag += fread(&hits->a.n, sizeof(hits->a.n), 1, fp);
    hits->a.m = hits->a.n; MALLOC(hits->a.a, hits->a.n);
    flag += fread(hits->a.a, sizeof(pe_hit), hits->a.n, fp);

    fclose(fp);
    free(buf);
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
        for (i = 0, f_w = r_w = 0; i < link->a.n; ++i) 
        { 
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

            fprintf(stderr, "self-utg%.6dl\tFW:%f\tRW:%f\tRT:%f\n**************************************************\n", 
            (int)(i+1), f_w, r_w, r_w/f_w);
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



uint64_t get_hic_distance(pe_hit* hit, hc_links* link, const ha_ug_index* idx)
{
    uint64_t s_uid, s_dir, e_uid, e_dir, u_dis, k;
    long long s_pos, e_pos;
    s_uid = ((hit->s<<1)>>(64 - idx->uID_bits)); s_pos = hit->s & idx->pos_mode;
    e_uid = ((hit->e<<1)>>(64 - idx->uID_bits)); e_pos = hit->e & idx->pos_mode;
    if(s_uid == e_uid) return MAX(s_pos, e_pos) - MIN(s_pos, e_pos);
    hc_linkeage* t = &(link->a.a[s_uid]);
    for (k = 0; k < t->e.n; k++)
    {
        if(t->e.a[k].del || t->e.a[k].uID != e_uid) continue;
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

inline double get_trans_weight(const ha_ug_index* idx, uint64_t x)
{
    ///return 1.0;
    long double rate = get_trans(idx, x);
    if(rate < 0) rate = 0;
    rate += OFFSET_RATE;
    if(rate > 0.5) rate = 0.5;
    rate -= OFFSET_SECOND_RATE; //[OFFSET_RATE - OFFSET_SECOND_RATE, 0.5 - OFFSET_SECOND_RATE]

    long double w = logl((1/rate)-1)*SCALL;
    if(w < OFFSET_RATE_MIN_W) w = OFFSET_RATE_MIN_W;
    if(w > OFFSET_RATE_MAX_W) w = OFFSET_RATE_MAX_W;
    return w;
}

inline double get_trans_weight_advance(const ha_ug_index* idx, uint64_t x, trans_idx* dis)
{
    long double rate = 0;
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
    
    if(rate < 0) rate = 0;
    rate += OFFSET_RATE;
    if(rate > 0.5) rate = 0.5;
    rate -= OFFSET_SECOND_RATE; //[OFFSET_RATE - OFFSET_SECOND_RATE, 0.5 - OFFSET_SECOND_RATE]

    long double w = logl((1/rate)-1)*SCALL;
    if(w < OFFSET_RATE_MIN_W) w = OFFSET_RATE_MIN_W;
    if(w > OFFSET_RATE_MAX_W) w = OFFSET_RATE_MAX_W;
    return w;
}

void LeastSquare(uint64_t* vec, uint64_t len, ha_ug_index* idx, uint64_t med)
{
    #define SCAL_RATE 1000
    long double t1=0, t2=0, t3=0, t4=0, x, y, thres;
    uint64_t i, len_convince, m;


    for (i = 0; i < len; i += 4)
    {
        if(vec[i+1] > med) break;
        x = ((double)(vec[i] + vec[i+1]))/2;
        y = ((double)(vec[i+3]))/((double)(vec[i+2] + vec[i+3]));

        t1 += x*x;
        t2 += x;
        t3 += x*y;
        t4 += y;
    }
    len_convince = i;

    if(t2 > t4)
    {
        idx->frac = t2/t4;
        if(idx->frac > SCAL_RATE) idx->frac = idx->frac / SCAL_RATE;
    }

    t1 /= (idx->frac*idx->frac);
    t2 /= idx->frac;
    t3 /= idx->frac;
    idx->a = idx->b = 0;
    if((t1*(len_convince>>2) - t2*t2) != 0)
    {
        idx->a = (t3*(len_convince>>2) - t2*t4) / (t1*(len_convince>>2) - t2*t2);
    }
    if((t1*(len_convince>>2) - t2*t2) != 0)
    {
        idx->b = (t1*t4 - t2*t3) / (t1*(len_convince>>2) - t2*t2);
    }
    
    

    if(len > 0)
    {
        vec[len - 3] = vec[len - 4] + (vec[1] - vec[0]);
    }
    if(len_convince >= len) return;

    thres = get_trans(idx, vec[len_convince] + vec[len_convince+1]);
    fprintf(stderr, "len_convince: %lu, len: %lu, t1: %f, t2: %f, t3: %f, t4: %f, idx->a: %f, idx->b: %f, thres: %f\n", 
                len_convince, len, (double)t1, (double)t2, (double)t3, (double)t4, (double)idx->a, (double)idx->b, (double)thres);


    
    for (i = m = 0; i < len; i += 4)
    {
        x = ((double)(vec[i] + vec[i+1]))/2;
        y = ((double)(vec[i+3]))/((double)(vec[i+2] + vec[i+3]));
        if(vec[i+1] > med && y < thres) continue;

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
    }

    len = m;
    t1 /= (idx->frac*idx->frac);
    t2 /= idx->frac;
    t3 /= idx->frac;
    ///fprintf(stderr, "len: %lu, t1: %f, t2: %f, t3: %f, t4: %f\n", len, (double)t1, (double)t2, (double)t3, (double)t4);
    if((t1*(len>>2) - t2*t2) != 0)
    {
        idx->a = (t3*(len>>2) - t2*t4) / (t1*(len>>2) - t2*t2);
    }
    if((t1*(len>>2) - t2*t2) != 0)
    {
        idx->b = (t1*t4 - t2*t3) / (t1*(len>>2) - t2*t2);
    }
    fprintf(stderr, "len: %lu, t1: %f, t2: %f, t3: %f, t4: %f, idx->a: %f, idx->b: %f\n", 
                len, (double)t1, (double)t2, (double)t3, (double)t4, (double)idx->a, (double)idx->b);

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


void weight_edges(ha_ug_index* idx, kvec_pe_hit* hits, hc_links* link, bubble_type* bub)
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
        
        t_d = get_hic_distance(&(hits->a.a[k]), link, idx);
        if(t_d == (uint64_t)-1) continue;

        e1 = get_hc_edge(link, beg, end, 0);
        e2 = get_hc_edge(link, end, beg, 0);
        if(e1 == NULL || e2 == NULL) continue;
        weight = get_trans_weight(idx, t_d);
        /*******************************for distance debug************************************/
        weight = 1;
        /*******************************for distance debug************************************/

        e1->weight += weight; e1->occ++;
        e2->weight += weight; e2->occ++;
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
        
        t_d = get_hic_distance(&(hits->a.a[k]), link, idx);
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

            pLen += bg->seq[rId].len;
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
        // if((beg>>1) == 162)
        // {
        //     fprintf(stderr, ">>>###0###>>>beg=btg%.6ul, beg&1: %u, cur=btg%.6ul, cur&1: %u, sink=btg%.6ul, sink&1: %u, tan_occ: %u\n", 
        //     (beg>>1)+1, beg&1, (cur>>1)+1, cur&1, (sink>>1)+1, sink&1, tan_occ);
        // }


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
                    btg_beg = v; btg_sink = cur^1;

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
                    utg_sink = w^1;

                    // fprintf(stderr, "\nbtg_beg=btg%.6ul, utg_beg=utg%.6ul\n", 
                    // (btg_beg>>1)+1, (utg_beg>>1)+1);
                    // fprintf(stderr, "btg_sink=btg%.6ul, utg_sink=utg%.6ul\n", 
                    // (btg_sink>>1)+1, (utg_sink>>1)+1);

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


void resolve_bubble_chain_tangle(ma_ug_t* ug, bubble_type* bub)
{
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
                    ///fprintf(stderr, "+>>>>>>beg=btg%.6ul, sink=btg%.6ul\n", (b_utg_beg>>1)+1, (b_utg_sink>>1)+1);
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
                    
                    update_bubble_graph(&res_utg, utg_beg, b_tg_beg, utg_sink, b_tg_sink, 
                    bub, &edges, bsg, NULL, NULL, ug, NULL, 0);
                    drop_g_edges_by_utg(bub, bsg, bub_ug, &res_btg, b_utg_beg, b_utg_sink);
                    ///fprintf(stderr, "->>>>>>beg=btg%.6ul, sink=btg%.6ul\n", (b_utg_beg>>1)+1, (b_utg_sink>>1)+1);
                }
            }
        }
        kv_push(uint32_t, bub->num, bub->list.n);
        new_bub = bub->b_g->n_seq - new_bub;
        bub->tangle_bub += new_bub;
        if(new_bub) update_bub_b_s_idx(bub);
        asg_arc_t *t = NULL;
        for (k = 0; k < edges.a.n; k++)
        {        
            t = asg_arc_pushp(bsg);
            *t = edges.a.a[k];
        }

        bsg->is_srt = 0; free(bsg->idx); bsg->idx = 0;
        asg_cleanup(bsg);
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
}


void update_bubble_chain(ma_ug_t* ug, bubble_type* bub, uint32_t is_middle, uint32_t is_end)
{
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

    for (i = 0; i < bs_ug->g->n_arc; i++)
    {
        e_w[i] = -1;
    }

    for (i = 0; i < bs_ug->g->n_seq; i++)
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
    if(new_bub) update_bub_b_s_idx(bub);
    
    ///fprintf(stderr, "bub->cross_bub: %u\n", (uint32_t)bub->cross_bub);
    ///debug_tangle_bubble(bub, bub->b_g->n_seq - bub->cross_bub, bub->b_g->n_seq - 1, "Cross-tangle");

    asg_arc_t *t = NULL;
    for (i = 0; i < edges.a.n; i++)
    {        
        t = asg_arc_pushp(bub->b_g);
        *t = edges.a.a[i];
    }
    bub->b_g->is_srt = 0; 
    free(bub->b_g->idx); 
    bub->b_g->idx = 0;
    asg_cleanup(bub->b_g);
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
}


void append_boundary_chain(ma_ug_t* ug, hc_links* link, bubble_type* bub)
{
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
    if(new_bub) update_bub_b_s_idx(bub);

    ///fprintf(stderr, "bub->mess_bub: %lu\n", bub->mess_bub);

    for (v = 0; v < n_vx; v++)
    {
        av = asg_arc_a(bs_ug->g, v);
        nv = asg_arc_n(bs_ug->g, v);
        if(nv == 0 || get_real_length(bs_ug->g, v, NULL) == 0) continue;
        drop_g_edges_by_utg(bub, bub->b_g, bs_ug, NULL, v, (uint32_t)-1);
    }


    asg_arc_t *t = NULL;
    for (i = 0; i < edges.a.n; i++)
    {        
        t = asg_arc_pushp(bub->b_g);
        *t = edges.a.a[i];
    }
    bub->b_g->is_srt = 0; 
    free(bub->b_g->idx); 
    bub->b_g->idx = 0;
    asg_cleanup(bub->b_g);
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




int get_trans_rate_function(ha_ug_index* idx, kvec_pe_hit* hits, hc_links* link, bubble_type* bub, MT* M, H_partition* hap, trans_idx* dis)
{
    kvec_t(uint64_t) buf, buf_idx;
    kv_init(buf);
    kv_init(buf_idx);
    uint64_t beg, end, cnt[2];
    uint64_t k, i, t_d, r_idx, f_idx, med = (uint64_t)-1;
    int beg_status, end_status;
    for (i = 0; i < link->a.n; i++)
    {
        for (k = 0; k < link->a.a[i].e.n; k++)
        {
            link->a.a[i].e.a[k].dis = (uint64_t)-1;
        }
    }
    fill_utg_distance_multi(idx, link, M, bub);


    buf.n = 0;
    for (k = 0; k < hits->a.n; ++k) 
    {
        beg = ((hits->a.a[k].s<<1)>>(64 - idx->uID_bits));
        end = ((hits->a.a[k].e<<1)>>(64 - idx->uID_bits));

        if(IF_HOM(beg, *bub)) continue;
        if(IF_HOM(end, *bub)) continue;


        t_d = get_hic_distance(&(hits->a.a[k]), link, idx);
        if(t_d == (uint64_t)-1) continue;
        if(beg == end)
        {
            t_d = (t_d << 1);
        }
        else
        {
            beg_status = get_phase_status(hap, beg);
            if(beg_status != 1 && beg_status != -1) continue;
            end_status = get_phase_status(hap, end);
            if(end_status != 1 && end_status != -1) continue;
            if(beg_status != end_status)
            {
                t_d = (t_d << 1) + 1; 
            }
            else
            {
                t_d = (t_d << 1);
            }            
        }

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

    for (k = 0; k < buf.n; k++)
    {
        if((buf.a[k]&1) == 1)
        {
            kv_push(uint64_t, buf_idx, buf.a[k]>>1); 
        } 
    }

    uint64_t cutoff = buf_idx.n * 0.9, t = buf_idx.n * 0.005, pre, step;
    k = 0;
    if(cutoff >= t) k = cutoff - t;
    pre = 0;
    if(cutoff >= t + 1) pre = buf_idx.a[cutoff - t - 1];
    for (t_d = i = 0; k < cutoff + t; k++)
    {
        t_d += (buf_idx.a[k] - pre);
        pre = buf_idx.a[k];
        i++;
    }
    
    if(t_d == 0 || i == 0 || t == 0)
    {
        kv_destroy(buf);
        kv_destroy(buf_idx);
        return 0;
    }

    step = (t_d/i)*20;

    if(step == 0)
    {
        kv_destroy(buf);
        kv_destroy(buf_idx);
        return 0;
    }

    trans_p_t* p = NULL;
    dis->n = 0;
    uint64_t step_s = 0, step_e = step;
    if(buf.n>0) step_s = buf.a[0]>>1, step_e = (buf.a[0]>>1) + step;
    for (k = cnt[0] = cnt[1] = 0; k < buf.n; k++)
    {
        if((buf.a[k]>>1) < step_e && (buf.a[k]>>1) >= step_s)
        {
            cnt[buf.a[k]&1]++;
        }

        if((buf.a[k]>>1) >= step_e)
        {
            while (!((buf.a[k]>>1) < step_e && (buf.a[k]>>1) >= step_s))
            {
                kv_pushp(trans_p_t, *dis, &p);
                p->beg = step_s;
                p->end = step_e;
                p->cnt_0 = cnt[0];
                p->cnt_1 = cnt[1];
                step_s += step;
                step_e += step;
                cnt[0] = cnt[1] = 0;
            }
        }
    }

    if(cnt[0] > 0 || cnt[1] > 0)
    {
        kv_pushp(trans_p_t, *dis, &p);
        p->beg = step_s;
        p->end = step_e;
        p->cnt_0 = cnt[0];
        p->cnt_1 = cnt[1];
    }


    uint64_t smooth_step = 20, k_i, cnt_0;
    if(dis->n > 0) med = dis->a[dis->n-1].end;
    for (k = 0; k+smooth_step < dis->n; k++)
    {
        for (k_i = cnt_0 = 0; k_i < smooth_step; k_i++)
        {
            if(dis->a[k+k_i].cnt_0 == 0 || dis->a[k+k_i].cnt_1 == 0) cnt_0++;
        }

        if(cnt_0 >= smooth_step * 0.2)
        {
            med = dis->a[k].beg;
            break;
        }
    }

   
    long long b_k = 0, b_i = 0, b_j, pass = 0;
    ///for (b_k = b_i = 0; b_k < (long long)dis->n; b_k++)
    while(b_k < (long long)dis->n)
    {
        pass = 1;
        beg = dis->a[b_k].beg;
        end = dis->a[b_k].end;
        cnt[0] = dis->a[b_k].cnt_0;
        cnt[1] = dis->a[b_k].cnt_1;
        if(cnt[0] > 0 && cnt[1] > 0)
        {
            dis->a[b_i].beg = beg;
            dis->a[b_i].end = end;
            dis->a[b_i].cnt_0 = cnt[0];
            dis->a[b_i].cnt_1 = cnt[1];
            b_i++;
            b_k++;
            continue;
        }

        b_k++;
        for (b_j = b_k; b_j < (long long)dis->n; b_j++, b_k++)
        {
            end = dis->a[b_j].end;
            cnt[0] += dis->a[b_j].cnt_0;
            cnt[1] += dis->a[b_j].cnt_1;
            if(cnt[0] > 0 && cnt[1] > 0) break;
        }


        if(b_j < (long long)dis->n)
        {
            dis->a[b_i].beg = beg;
            dis->a[b_i].end = end;
            dis->a[b_i].cnt_0 = cnt[0];
            dis->a[b_i].cnt_1 = cnt[1];
            b_i++;
            b_k++;
            continue;
        }

        for(b_j = b_i-1; b_j >= 0; b_j--)
        {
            beg = dis->a[b_j].beg;
            cnt[0] += dis->a[b_j].cnt_0;
            cnt[1] += dis->a[b_j].cnt_1;
            if(cnt[0] > 0 && cnt[1] > 0) break;
        }

        if(b_j >= 0)
        {
            b_i = b_j;
            dis->a[b_i].beg = beg;
            dis->a[b_i].end = end;
            dis->a[b_i].cnt_0 = cnt[0];
            dis->a[b_i].cnt_1 = cnt[1];
            b_i++;
            b_k++;
            continue;
        }

        pass = 0;
        break;
    }

    dis->n = b_i;
    if(dis->n == 0 || pass == 0)
    {
        kv_destroy(buf);
        kv_destroy(buf_idx);
        return 0;
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

    
    kv_destroy(buf);
    kv_destroy(buf_idx);

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

    return 1;
}

void init_hic_p(ha_ug_index* idx, kvec_pe_hit* hits, hc_links* link, bubble_type* bub, 
kvec_hc_edge* back_hc_edge, MT* M, H_partition* hap, uint32_t ignore_dis)
{
    uint64_t k, i, m, uID, is_comples_weight = 0;
    trans_idx dis;
    kv_init(dis);

    if(bub->round_id > 0 && ignore_dis == 0)
    {
        is_comples_weight = get_trans_rate_function(idx, hits, link, bub, M, hap, &dis);
    }


    hc_edge *e = NULL;
    for (i = 0; i < link->a.n; i++)
    {
        for (k = 0; k < link->a.a[i].f.n; k++)
        {
            if(link->a.a[i].f.a[k].del) continue;
            if(link->a.a[i].f.a[k].dis == RC_0)
            {
                uID = link->a.a[i].f.a[k].uID;
                e = get_hc_edge(link, i, uID, 0);
                if(e) 
                {
                    if(back_hc_edge) kv_push(hc_edge, back_hc_edge->a, *e);
                    e->del = 1;
                }

                e = get_hc_edge(link, uID, i, 0);
                if(e) 
                {
                    if(back_hc_edge) kv_push(hc_edge, back_hc_edge->a, *e);
                    e->del = 1;
                }
            }
            else if(link->a.a[i].f.a[k].dis == RC_1)
            {
                uID = link->a.a[i].f.a[k].uID;
                get_forward_distance(i, uID, idx->ug->g, link, M);
                get_forward_distance(uID, i, idx->ug->g, link, M);
            }
        }
    }
    

    for (i = 0; i < link->a.n; i++)
    {
        for (k = 0; k < link->a.a[i].e.n; k++)
        {
            if(link->a.a[i].e.a[k].del) continue;
            if(link->a.a[i].e.a[k].dis == (uint64_t)-1)
            {
                e = get_hc_edge(link, link->a.a[i].e.a[k].uID, i, 0);
                if(back_hc_edge) kv_push(hc_edge, back_hc_edge->a, link->a.a[i].e.a[k]);
                if(back_hc_edge) kv_push(hc_edge, back_hc_edge->a, *e);
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
                if(back_hc_edge) kv_push(hc_edge, back_hc_edge->a, link->a.a[i].e.a[k]);
                if(back_hc_edge) kv_push(hc_edge, back_hc_edge->a, *e);
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
    kv_resize(double, hap->label_buffer, (hap->label>>hap->label_shift)+1);
    hap->label_buffer.n = (hap->label>>hap->label_shift)+1;
    uint32_t i, k, m, is_ava, a_n, *x_a, x_n, x;
    uint64_t bid;
    hc_edge* a = NULL;
    for (i = 0; i < hap->label_buffer.n; i++)
    {
        hap->label_buffer.a[i] = 0;
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
    
    if(is_ava == 0) return; ///this is a totally new chain
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

            if(current_weight > max_weight)
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
    ///fprintf(stderr, "\n######max_bid_idx: %u, is_forward: %u, max_hap_label: %u\n", max_bid_idx, is_forward, max_hap_label);
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
uint32_t phasing_improvement(H_partition* h, G_partition* g_p, ha_ug_index* idx, bubble_type* bub);
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
/**
double get_add_weight(H_partition* h, G_partition* g_p, hc_links* link, block_phase_type* block,
bubble_type* bub, uint32_t gid)
{
    double weight = 0;
    if(g_p->a[gid].h[0] > 0 && (g_p->a[gid].status[0] == 1 || g_p->a[gid].status[0] == -1))
    {
        weight += (g_p->a[gid].weight[0] * g_p->a[gid].status[0]);
    }

    if(g_p->a[gid].h[1] > 0 && (g_p->a[gid].status[1] == 1 || g_p->a[gid].status[1] == -1))
    {
        weight += (g_p->a[gid].weight[1] * g_p->a[gid].status[1]);
    }
}

void update_block_weight(H_partition* hap, G_partition* g_p, bubble_type* bub, block_phase_type* block, 
uint64_t bid)
{
    uint32_t beg, sink, k, uid, *a, n, gid;
    get_bubbles(bub, bid, &beg, &sink, &a, &n, NULL);
    for (k = 0; k < n; k++)
    {
        uid = a[k]>>1;
        if(g_p->index[uid] == (uint32_t)-1) continue;
        gid = g_p->index[uid]>>1;
        if(block->vis.a[gid]) continue;

    }

}

**/

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
    uint64_t* chain_idx;
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
    x->chain_idx[i] = n;
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


// static void worker_for_max_block_by_chain(void *data, long i, int tid) // callback for kt_for()
// {
//     mul_block_phase_type* x = (mul_block_phase_type*)data;
//     select_max_block_by_utg_multi_thread_by_chain(x->h, x->g_p, x->bub, x->h->link, 
//     &(x->x[tid]), x->bub->chain_weight.a[i].id, &(x->res[tid]));
// }

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
    for (i = 0; i < n; i++)
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
}

uint32_t init_contig_partition(H_partition* hap, ha_ug_index* idx, bubble_type* bub)
{
    hc_links* link = idx->link;
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
    hap->label_add = 8;
    for(hap->label_shift=1; (uint64_t)(1<<hap->label_shift)<(uint64_t)hap->label_add; hap->label_shift++);
    kv_init(hap->label_buffer);
    kv_init(hap->b.vis); kv_malloc(hap->b.vis, hap->n); hap->b.vis.n = hap->n;

    ///fprintf(stderr, "hap->label: %u, hap->label_add: %u, hap->label_shift: %u\n", hap->label, hap->label_add, hap->label_shift);


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

    assign_per_unitig_G_partition(&(hap->g_p), hap->n, link, bub, 0);

    adjust_contig_partition(hap, link);

    update_bubble_chain(ug, bub, 0, 1);

    resolve_bubble_chain_tangle(ug, bub);

    clean_bubble_chain_by_HiC(ug, link, bub);

    append_boundary_chain(ug, link, bub);

    sort_bubble_ele(&(hap->g_p), link, bub, hap->n);

    free(hap_label_flag);
    return 1;
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


uint32_t phasing_improvement(H_partition* h, G_partition* g_p, ha_ug_index* idx, bubble_type* bub)
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
                i = get_max_unitig(h, g_p, idx->link, bub);
                if(i == (uint32_t)-1) break;
                h->lock[i] = 1;
                flip_unitig(g_p, idx->link, bub, i);
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
        update_partition_flag(h, g_p, idx->link, i);
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
    uint32_t m, i;
    uint32_t beg, sink;
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
            get_bubbles(bub, bid, &beg, &sink, NULL, NULL, NULL);
            fprintf(stderr, "btg%.6lu%c, beg-utg%.6ul, sink-utg%.6ul\n", 
                        bid, "fb"[bid<bub->f_bub?0:1], (beg>>1)+1, (sink>>1)+1);
        }
    }
}

void init_contig_H_partition(bubble_type* bub, ha_ug_index* idx, kvec_pe_hit* hits, H_partition* hap)
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

void cluster_contigs(bubble_type* bub, ha_ug_index* idx, kvec_pe_hit* hits, MT* M, H_partition* hap)
{   
    uint64_t k, i, shif = 64 - idx->uID_bits, beg, end, t_d;
    hc_links* link = idx->link;
    for (i = 0; i < link->a.n; i++) link->a.a[i].e.n = 0;
    for (k = 0; k < hits->a.n; ++k) 
    {
        beg = ((hits->a.a[k].s<<1)>>shif);
        end = ((hits->a.a[k].e<<1)>>shif);

        if(beg == end) continue;
        if(IF_HOM(beg, *bub)) continue;
        if(IF_HOM(end, *bub)) continue;

        t_d = 1;
        push_hc_edge(&(link->a.a[beg]), end, 0, 0, &t_d);
        push_hc_edge(&(link->a.a[end]), beg, 0, 0, &t_d);
    }

    init_hic_p((ha_ug_index*)idx, hits, link, bub, NULL, M, NULL, 1);
    
    init_chain_hic_warp(idx->ug, link, bub, &bub->c_w);

    hap->link = link;
    hap->n = idx->ug->u.n;

    init_contig_H_partition(bub, idx, hits, hap);

    destory_chain_hic_warp(&bub->c_w);
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
    int i;
    for (i = 0; i < fn1->n && i < fn2->n; i++)
    {
        gzFile fp1, fp2;
        if ((fp1 = gzopen(fn1->a[i], "r")) == 0) return 0;
        if ((fp2 = gzopen(fn2->a[i], "r")) == 0) return 0;
        sl->ks1 = kseq_init(fp1);
        sl->ks2 = kseq_init(fp2);

        kt_pipeline(3, worker_pipeline, sl, 3);

        ///fprintf(stderr, "fn1->a[i]: %s, fn2->a[i]: %s, sl->hits.a.n: %u\n", fn1->a[i], fn2->a[i], (uint32_t)sl->hits.a.n);
        
        
        kseq_destroy(sl->ks1);
        kseq_destroy(sl->ks2);
        gzclose(fp1);
        gzclose(fp2);
    }

    ///fprintf(stderr, "+sl->hits.a.n: %u\n", (uint32_t)sl->hits.a.n);

    dedup_hits(&(sl->hits));

    ///fprintf(stderr, "-sl->hits.a.n: %u\n", (uint32_t)sl->hits.a.n);
    
    return 1;
}

int hic_short_align(const enzyme *fn1, const enzyme *fn2, ha_ug_index* idx)
{
    double index_time = yak_realtime();
    sldat_t sl;
    kvec_hc_edge back_hc_edge;
    kv_init(back_hc_edge.a);
    sl.idx = idx;
    sl.link = idx->link;
    sl.chunk_size = 20000000;
    sl.n_thread = asm_opt.thread_num;
    sl.total_base = sl.total_pair = 0;
    idx->max_cnt = 5;
    kv_init(sl.hits.a);
    
    if(!load_hc_hits(&sl.hits, asm_opt.output_file_name))
    {
        /*******************************for debug************************************/
        // load_reads(&R1, fn1);
        // test_reads(&R1, fn1);
        // load_reads(&R2, fn2);
        // test_reads(&R1, fn1);
        /*******************************for debug************************************/

        // kt_pipeline(3, worker_pipeline, &sl, 3);
        // dedup_hits(&sl.hits);
        alignment_worker_pipeline(&sl, fn1, fn2);
        /*******************************for debug************************************/
        // sort_hits(&sl.hits);
        // print_hits(idx, &sl.hits, fn1);
        /*******************************for debug************************************/
        
        write_hc_hits(&sl.hits, asm_opt.output_file_name);
    }

    ///fprintf(stderr, "u.n: %d, uID_bits: %lu, pos_bits: %lu, sl.hits.a.n: %u\n", (uint32_t)idx->ug->u.n, idx->uID_bits, idx->pos_bits, (uint32_t)sl.hits.a.n);

    H_partition hap;
    MT M;
    init_MT(&M, idx->ug->g->n_seq<<1);
    bubble_type bub; 
    memset(&bub, 0, sizeof(bubble_type));
    bub.round_id = 0; bub.n_round = 2;
    for (bub.round_id = 0; bub.round_id < bub.n_round; bub.round_id++)
    {
        identify_bubbles(idx->ug, &bub, idx->link);
        if(bub.round_id == 0)
        {
            collect_hc_links(sl.idx, &sl.hits, idx->link, &bub, &M);
            collect_hc_reverse_links(idx->link, idx->ug, &bub);
        }
        init_hic_p((ha_ug_index*)sl.idx, &sl.hits, idx->link, &bub, &back_hc_edge, &M, &hap, 0);
        ///init_hic_p_new((ha_ug_index*)sl.idx, &sl.hits, idx->link, &bub, &back_hc_edge, &M);
        reset_H_partition(&hap, (bub.round_id == 0? 1 : 0));
        init_contig_partition(&hap, idx, &bub);
        phasing_improvement(&hap, &(hap.g_p), idx, &bub);
        label_unitigs(&(hap.g_p), idx->ug);

        ///print_hc_links(idx->link, 0, &hap);
    }

    cluster_contigs(&bub, idx, &sl.hits, &M, &hap);

    destory_MT(&M);

    ///print_bubbles(idx->ug, &bub, sl.hits.a.n?&sl.hits:NULL, idx->link, idx);
    ///print_hits(idx, &sl.hits, fn1);
    

    ///print_debug_bubble_graph(&bub, idx->ug, asm_opt.output_file_name);
    // print_bubble_chain(&bub);
    // print_hc_links(idx->link, 0, &hap);
    
    ///print_contig_partition(&hap, "final");

    // uint32_t i;
    // for (i = 0; i < idx->ug->g->n_seq; i++)
    // {
    //     fprintf(stderr, "utg%.6ul, index: %u\n", (int)(i+1), bub.index[i]);
    // }
    



    destory_contig_partition(&hap);
    kv_destroy(back_hc_edge.a);
    return 1;
    
    /*******************************for debug************************************/
    // destory_reads(&R1);
    // destory_reads(&R2);
    /*******************************for debug************************************/
    print_bubbles(idx->ug, &bub, sl.hits.a.n?&sl.hits:NULL, idx->link, idx);
    collect_hc_reverse_links(idx->link, idx->ug, &bub);
    normalize_hc_links(idx->link);
    /*******************************for debug************************************/
    ///print_hc_links(&link);
    /*******************************for debug************************************/
    min_cut_t* cut = clean_hap(idx->link, &bub, idx->ug);
    ///print_bubbles(idx->ug, &bub, NULL, &link, idx);
    G_partition* gp = clean_bubbles(idx->link, &bub, cut, idx->ug);
    ///print_hc_links(&link);

    destory_min_cut_t(cut); free(cut);
    destory_G_partition(gp); free(gp);
    kv_destroy(sl.hits.a);
    destory_bubbles(&bub);
    
    fprintf(stderr, "[M::%s::%.3f] processed %lu pairs; %lu bases\n", __func__, yak_realtime()-index_time, sl.total_pair, sl.total_base);
    return 1;
}


void hic_analysis(ma_ug_t *ug, asg_t* read_g, hc_links* link)
{
    ug_index = NULL;
    int exist = load_hc_pt_index(&ug_index, asm_opt.output_file_name);
    if(exist == 0) ug_index = build_unitig_index(ug, asm_opt.hic_mer_length);
    if(exist == 0) write_hc_pt_index(ug_index, asm_opt.output_file_name);
    ug_index->ug = ug;
    ug_index->read_g = read_g;
    ug_index->link = link;
    ///test_unitig_index(ug_index, ug);
    hic_short_align(asm_opt.hic_reads[0], asm_opt.hic_reads[1], ug_index);
    
    destory_hc_pt_index(ug_index);
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

    init_hc_links(&(idx->link), ug->u.n, ug->g->n_seq);

    for (i = 0; i < idx->ug_idx.n; i++)
    {
        get_r_utg_bench(i, idx, ug);
    }
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
    sl.link = NULL;
    sl.chunk_size = 20000000;
    sl.n_thread = asm_opt.thread_num;
    sl.total_base = sl.total_pair = 0;
    idx->max_cnt = 5;
    kv_init(sl.hits.a);
    fprintf(stderr, "u.n: %d, uID_bits: %lu, pos_bits: %lu\n", (uint32_t)idx->ug->u.n, idx->uID_bits, idx->pos_bits);
    
    if(!load_hc_hits(&sl.hits, output_file_name))
    {
        // kt_pipeline(3, worker_pipeline, &sl, 3);
        // dedup_hits(&sl.hits);
        alignment_worker_pipeline(&sl, fn1, fn2);
        write_hc_hits(&sl.hits, output_file_name);
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
    int exist = load_hc_pt_index(&ug_index, output_file_name);
    if(exist == 0) ug_index = build_unitig_index(ug, asm_opt.hic_mer_length);
    if(exist == 0) write_hc_pt_index(ug_index, output_file_name);
    ug_index->ug = ug;
    ug_index->read_g = read_g;

    hic_short_align_bench(asm_opt.hic_reads[0], asm_opt.hic_reads[1], output_file_name, ug_index);

    free(output_file_name);
}