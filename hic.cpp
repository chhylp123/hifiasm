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
    double weight[2];
}partition_warp;

typedef struct{
    size_t n, m; 
    partition_warp* a;
    uint32_t* index;
}G_partition;

typedef struct{
    uint64_t n; 
    uint8_t* lock;
    uint32_t* hap;
    uint32_t m[3];
    uint32_t label, label_add;
    hc_links* link;
    G_partition g_p;
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
} bub_p_t_warp;


typedef struct {
	hc_pt_t *h;
	uint64_t n;
	uint64_t *a;
    khint_t end;
} hc_pt1_t;

typedef struct {
    uint32_t* index;
    ma_ug_t* ug;
    kvec_t(uint32_t) list;
    kvec_t(uint32_t) num;
    kvec_t(uint64_t) pathLen;
} bubble_type;

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
    hc_pt1_t* idx_buf;
    long double a,b, frac;
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

#define pe_hit_an1_key(x) ((x).s)
KRADIX_SORT_INIT(pe_hit_an1, pe_hit, pe_hit_an1_key, 8)
#define pe_hit_an2_key(x) ((x).e)
KRADIX_SORT_INIT(pe_hit_an2, pe_hit, pe_hit_an2_key, 8)
#define generic_key(x) (x)
KRADIX_SORT_INIT(hc64, uint64_t, generic_key, 8)
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

#define Get_bub_num(RECORD) ((RECORD).num.n-1)

reads_t R1, R2;
ha_ug_index* ug_index;

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

void hc_pt_t_gen_single(hc_pt1_t* pt)
{
    khint_t k;
    uint64_t c;
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
    sprintf(gfa_name, "%s.hc_tlb", file_name);
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
    sprintf(gfa_name, "%s.hc_tlb", file_name);
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
            radix_sort_hc_pos(pos_list, pos_list+cnt);
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
            hc_pt_t_gen_single(&(idx->idx_buf[i]));
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
                
                k_n=(end-beg+1);pos_k=pos_list[num-1];pos_list[num-1]+=k_n;
                for (k = 0; k < k_n; k++)
                {
                    pos_list[pos_k+k] = pos[beg+k].pos;
                }
                beg = end = m;
            }
            
        }
        if(occ > 0)
        {
            num = get_hc_pt1_count(pl->h, pos[beg].key, &pos_list);
            k_n=(end-beg+1);pos_k=pos_list[num-1];pos_list[num-1]+=k_n;
            for (k = 0; k < k_n; k++)
            {
                pos_list[pos_k+k] = pos[beg+k].pos;
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
                    ///p->off_cnt = self_p | ((uint64_t)cnt << 32);
                    p->off_cnt = self_p | ((uint64_t)k_mer << 32);

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

inline void set_pe_pos(ha_ug_index* idx, s_hit *l1, uint64_t occ1, s_hit *l2, uint64_t occ2, 
pe_hit* x, uint64_t rid)
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
    
    set_pe_pos((ha_ug_index*)s->idx, s->pos_buf[tid].a.a, occ1, s->pos_buf[tid].a.a + occ1, occ2, &(s->pos[i]), s->id+i);

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
        s->idx = p->idx; s->id = p->total_pair;
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
        fprintf(stderr, "%.*s\t%c\tutg%.6d\t%lu\t%c\tutg%.6d\t%lu\ti:%lu\n", 
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
}

inline void get_bubbles(bubble_type* bub, uint64_t id, uint32_t* beg, uint32_t* sink, uint32_t** a, uint32_t* n, uint64_t* pathBase)
{
    (*a) = bub->list.a + bub->num.a[id] + 2;
    (*n) = bub->num.a[id+1] - bub->num.a[id] - 2;
    (*beg) = bub->list.a[bub->num.a[id]];
    (*sink) = bub->list.a[bub->num.a[id] + 1];
    if(pathBase) (*pathBase) = bub->pathLen.a[id];
}

void identify_bubbles(ma_ug_t* ug, bubble_type* bub)
{
    asg_cleanup(ug->g);
    if (!ug->g->is_symm) asg_symm(ug->g);
    memset(bub, 0, sizeof(bubble_type));
    uint32_t v, n_vtx = ug->g->n_seq * 2, tLen, i, mode = (((uint32_t)-1)<<2);
    uint64_t pathLen;
    bub->ug = ug; 
    CALLOC(bub->index, n_vtx);
    for (i = 0; i < ug->g->n_seq; i++)
    {
        if(ug->g->seq[i].c > 0)
        {
            bub->index[i] = (ug->g->seq[i].c << 2);
            ug->g->seq[i].c = 0;
        } 
    }
    kv_init(bub->list); kv_init(bub->num); kv_init(bub->pathLen);
    buf_t b; memset(&b, 0, sizeof(buf_t)); b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    for (i = 0, tLen = 1; i < ug->u.n; i++) tLen += ug->u.a[i].len;
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(ug->g->seq[v>>1].del) continue;
        if(asg_arc_n(ug->g, v) < 2) continue;
        if((bub->index[v]&(uint32_t)3) != 0) continue;
        if(asg_bub_pop1_primary_trio(ug->g, NULL, v, tLen, &b, (uint32_t)-1, (uint32_t)-1, 0, NULL))
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

    
    for (v = 0; v < n_vtx; ++v) 
    {
        if((bub->index[v]&(uint32_t)3) !=2) continue;
        kv_push(uint32_t, bub->num, bub->list.n);
        if(asg_bub_pop1_primary_trio(ug->g, NULL, v, tLen, &b, (uint32_t)-1, (uint32_t)-1, 0, &pathLen))
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
                bub->index[i] = bub->num.n;
            }
            else
            {
                bub->index[i] = bub->num.n + 1;
            }
        }         
    }

    uint32_t beg, sink, n, *a;
    for (i = 0; i < bub->num.n-1; i++)
    {
        get_bubbles(bub, i, &beg, &sink, &a, &n, NULL);
        for (v = 0; v < n; v++)
        {
            bub->index[(a[v]>>1)] = i;
        }

        if(bub->index[(beg>>1)] != (bub->num.n + 1)) bub->index[(beg>>1)] = (uint32_t)-1;
        if(bub->index[(sink>>1)] != (bub->num.n + 1)) bub->index[(sink>>1)] = (uint32_t)-1;
        ///bub->index[(beg>>1)] = bub->index[(sink>>1)] = (uint32_t)-1;
    }

    for (i = 0; i < ug->g->n_seq; i++)
    {
        if(bub->index[i] == bub->num.n + 1) bub->index[i] = bub->num.n;
    }

    ///free(bub->index); bub->index = NULL;
}



void print_bubbles(ma_ug_t* ug, bubble_type* bub, kvec_pe_hit* hits, hc_links* link, ha_ug_index* idx)
{
    uint64_t tLen, t_utg, i, k;
    uint32_t beg, sink, n, *a;
    for (i = 0, tLen = 0; i < bub->ug->u.n; i++) tLen += bub->ug->u.a[i].len;
    fprintf(stderr, "[M::%s] # unitigs: %lu, # bases: %lu\n",  __func__, bub->ug->u.n, tLen);
    for (i = 0, tLen = 0, t_utg = 0; i < bub->num.n-1; i++)
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
        if(bub->index[i] < bub->num.n)
        {
            t_utg++;
            tLen +=bub->ug->u.a[i].len;
        }
    }
    fprintf(stderr, "[M::%s] # bubbles: %lu, # unitigs: %lu, # bases: %lu\n",  __func__, 
    Get_bub_num(*bub), t_utg, tLen);

    for (i = 0, tLen = 0, t_utg = 0; i < ug->g->n_seq; i++)
    {
        if(bub->index[i] == bub->num.n)
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
            if(bub->index[s_uid] < bub->num.n && bub->index[e_uid] < bub->num.n)
            {
                flag[s_uid] |= 1;
                flag[e_uid] |= 1;
                continue;
            }
            if(bub->index[s_uid] == bub->num.n && bub->index[e_uid] == bub->num.n)
            {
                flag[s_uid] |= 4;
                flag[e_uid] |= 4;
                continue;
            }
            if(bub->index[s_uid] < bub->num.n) flag[s_uid] |= 2, flag[e_uid] |= 2;
            if(bub->index[e_uid] < bub->num.n) flag[e_uid] |= 2, flag[s_uid] |= 2;
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
                if(bub->index[s_uid] < bub->num.n && bub->index[e_uid] < bub->num.n)
                {
                    flag[s_uid] |= 1;
                    flag[e_uid] |= 1;
                    continue;
                }
                if(bub->index[s_uid] == bub->num.n && bub->index[e_uid] == bub->num.n)
                {
                    flag[s_uid] |= 4;
                    flag[e_uid] |= 4;
                    continue;
                }
                if(bub->index[s_uid] < bub->num.n) flag[s_uid] |= 2, flag[e_uid] |= 2;
                if(bub->index[e_uid] < bub->num.n) flag[e_uid] |= 2, flag[s_uid] |= 2;
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
    for (i = 0, tLen = 0, t_utg = 0; i < bub->num.n-1; i++)
    {
        get_bubbles(bub, i, &beg, &sink, &a, &n, &pathLen);
        t_utg += n;
        fprintf(stderr, "(%lu)\tbeg:utg%.6u\tsink:utg%.6u\tpathLen:%lu\n", i, (beg>>1)+1, (sink>>1)+1, pathLen);
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
    //     if(bub->index[i] == bub->num.n) fprintf(stderr, "utg%.6lu\n", i+1);
    // }
    // fprintf(stderr, "************het utgs************\n");
}





void push_hc_edge(hc_linkeage* x, uint64_t uID, int weight, int dir, uint64_t* d)
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


void get_shortest_path(uint32_t src, pdq* pq, asg_t *sg)
{
    uint64_t v, u, i, nv, w;
    asg_arc_t *av = NULL;
    reset_pdq(pq);
    pq->dis.a[src] = 0;
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
        get_shortest_path(v, &pq, sg);
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
        if(x_i == M->matrix.a[x].a.n) fprintf(stderr, "ERROR X\n");
        d_x = d;

        for (; y_i < M->matrix.a[y].a.n; y_i++)
        {
            u = M->matrix.a[y].a.a[y_i] >> M->uID_shift;
            d = M->matrix.a[y].a.a[y_i] & M->dis_mode;
            if(u == j) break;
        }
        if(y_i == M->matrix.a[y].a.n) fprintf(stderr, "ERROR Y\n");
        d_y = d;

        tmp = LCA_distance(d_x, d_y, xLen, yLen, &rev);
        if(tmp < min_d) min_d = tmp, (*min_rev) = rev, min_j = j;
    }

    if(min_j == x || min_j == y) return (uint64_t)-1;

    return min_d;


}

uint64_t get_LCA(uint32_t x, uint64_t xLen, uint32_t y, uint64_t yLen, uint8_t* dis, uint64_t n, MT* M, bubble_type* bub, uint64_t* min_rev)
{
    

    if(bub->index[x>>1] < bub->num.n && bub->index[y>>1] < bub->num.n &&
       bub->index[x>>1] == bub->index[y>>1])
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
        if(x_i == M->matrix.a[x].a.n) fprintf(stderr, "ERROR X\n");
        d_x = d;

        for (; y_i < M->matrix.a[y].a.n; y_i++)
        {
            u = M->matrix.a[y].a.a[y_i] >> M->uID_shift;
            d = M->matrix.a[y].a.a[y_i] & M->dis_mode;
            if(u == j) break;
        }
        if(y_i == M->matrix.a[y].a.n) fprintf(stderr, "ERROR Y\n");
        d_y = d;

        tmp = LCA_distance(d_x, d_y, xLen, yLen, &rev);
        if(tmp < min_d) min_d = tmp, (*min_rev) = rev, min_j = j;
    }

    if(min_j == x || min_j == y) return (uint64_t)-1;

    return min_d;
}

void fill_utg_distance(const ha_ug_index* idx, hc_links* link, MT* M, bubble_type* bub)
{
    double index_time = yak_realtime();
    asg_t *sg = idx->ug->g;
    hc_linkeage* t = NULL;
    uint32_t n_vtx = sg->n_seq<<1, v, u, k, j, i;
    uint64_t d[2], db[2], q_u, min, min_i, min_b, rev[2], min_rev;
    kvec_t(uint8_t) dis_buf;
    kv_malloc(dis_buf, n_vtx); dis_buf.n = n_vtx;
    ///for (v = 0; v < n_vtx; ++v)
    for (i = 0; i < sg->n_seq; i++)
    {
        if (sg->seq[i].del) continue;
        t = &(link->a.a[i]);
        if (t->e.n == 0) continue;

        for (k = 0; k < t->e.n; k++)
        {
            if(t->e.a[k].del) continue;
            u = t->e.a[k].uID;

            for (v = (i<<1); v < ((i+1)<<1); v++)
            {
                d[0] = d[1] = db[0] = db[1] = (uint64_t)-1;
                
                for (j = 0; j < M->matrix.a[v].a.n; j++)
                {
                    q_u = M->matrix.a[v].a.a[j] >> M->uID_shift;
                    if((q_u>>1) == u) d[q_u&1] = (M->matrix.a[v].a.a[j] & M->dis_mode) + sg->seq[q_u>>1].len;
                    if((q_u>>1) > u) break;
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


            if(bub->index[i] >= bub->num.n && bub->index[u] >= bub->num.n && 
                bub->index[i] == bub->index[u] && t->e.a[k].dis != (uint64_t)-1)
            {
                continue;
            }



            for (v = (i<<1); v < ((i+1)<<1); v++)
            {
                d[0] = d[1] = db[0] = db[1] = (uint64_t)-1;  
                db[0] = get_LCA(v, sg->seq[v>>1].len, u<<1, sg->seq[u].len, 
                                                            dis_buf.a, dis_buf.n, M, bub, &rev[0]);
                db[1] = get_LCA(v, sg->seq[v>>1].len, (u<<1) + 1, sg->seq[u].len, 
                                                            dis_buf.a, dis_buf.n, M, bub, &rev[1]);
        

                min = min_i = min_b = min_rev = (uint64_t)-1;
                if(t->e.a[k].dis != (uint64_t)-1) min = t->e.a[k].dis >> 3;
                
                if(db[0] < min) min = db[0], min_i = 0, min_b = 1, min_rev = rev[0];
                if(db[1] < min) min = db[1], min_i = 1, min_b = 1, min_rev = rev[1];

                if(min_i != (uint64_t)-1 && min != (uint64_t)-1)
                {
                    t->e.a[k].dis = min<<1;
                    t->e.a[k].dis += (min_b^min_rev);
                    t->e.a[k].dis <<=1;
                    t->e.a[k].dis += ((v&1)^min_rev);
                    t->e.a[k].dis <<=1;
                    t->e.a[k].dis += min_i;
                }
            }

        }
    }
    kv_destroy(dis_buf);
    fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
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
                if((q_u>>1) > u) break;
            }

            min = min_i = min_b = (uint64_t)-1;
            if(t->e.a[k].dis != (uint64_t)-1) min = t->e.a[k].dis >> 3;
            
            if(d[0] < min) min = d[0], min_i = 0, min_b = 0;
            if(d[1] < min) min = d[1], min_i = 1, min_b = 0;

            // if((v>>1) == 8185 && u == 3845)
            // {
            //     fprintf(stderr, "***v>>1: %u, v&1: %u, u: %u, d[0]: %lu, d[1]: %lu, db[0]: %lu, db[1]: %lu\n", 
            //     v>>1, v&1, u, d[0], d[1], db[0], db[1]);
            // }

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
        if(bub->index[i] < bub->num.n && bub->index[u] < bub->num.n 
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

            // if((v>>1) == 8185 && u == 3845)
            // {
            //     fprintf(stderr, "###v>>1: %u, v&1: %u, u: %u, d[0]: %lu, d[1]: %lu, db[0]: %lu, db[1]: %lu\n", 
            //     v>>1, v&1, u, d[0], d[1], db[0], db[1]);
            // }

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


        // if(bub->index[i] < bub->num.n && bub->index[u] < bub->num.n && 
        //    bub->index[i] == bub->index[u])
        // {
        //     if(t->e.a[k].dis == (uint64_t)-1) fprintf(stderr, "hahahahahahahaha\n");
        //     fprintf(stderr, "sb-utg%.6d\tdb-utg%.6d\n", (int)(i+1), (int)(u+1));
        // }

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

void collect_hc_links(const ha_ug_index* idx, kvec_pe_hit* hits, hc_links* link, bubble_type* bub)
{
    double index_time = yak_realtime();
    uint64_t k, i, shif = 64 - idx->uID_bits, beg, end, t_d;
    for (k = 0; k < hits->a.n; ++k) 
    {
        beg = ((hits->a.a[k].s<<1)>>shif);
        end = ((hits->a.a[k].e<<1)>>shif);

        if(beg == end) continue;
        if(bub->index[beg] > bub->num.n) continue;
        if(bub->index[end] > bub->num.n) continue;

        t_d = (uint64_t)-1;
        push_hc_edge(&(link->a.a[beg]), end, 0, 0, &t_d);
        push_hc_edge(&(link->a.a[end]), beg, 0, 0, &t_d);
    }

    uint32_t n_vtx = idx->ug->g->n_seq<<1, v;
    MT M;
    kv_init(M.matrix); kv_malloc(M.matrix, n_vtx); M.matrix.n = n_vtx;
    for (v = 0; v < n_vtx; ++v) kv_init(M.matrix.a[v].a);
    for (v = 1; (uint64_t)(1<<v) < n_vtx; v++);
    M.uID_shift = 64 - v; M.dis_mode = ((uint64_t)-1) >> v; 

    all_pair_shortest_path(idx, link, &M);
    ///fill_utg_distance(idx, link, &M, bub);
    fill_utg_distance_multi(idx, link, &M, bub);

    for (v = 0; v < n_vtx; ++v) kv_destroy(M.matrix.a[v].a);
    kv_destroy(M.matrix);
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

void dfs_bubble(asg_t *g, kvec_t_u32_warp* stack, kvec_t_u32_warp* result, 
uint32_t v, uint32_t beg, uint32_t sink)
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
        kv_push(uint32_t, result->a, cur>>1);
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
        kv_push(uint32_t, result->a, cur>>1);
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

void set_reverse_links(uint32_t* bub, uint32_t n, kvec_t_u32_warp* reach, uint32_t root, hc_links* link)
{
    uint64_t i, k, d = 0;
    uint32_t v;
    for (i = 0; i < n; i++)
    {
        v = bub[i]>>1;
        for (k = 0; k < reach->a.n; k++)
        {
            if(v == reach->a.a[k]) break;
        }

        if(k == reach->a.n)
        {
            push_hc_edge(&(link->a.a[root]), v, 1, 1, &d);
            push_hc_edge(&(link->a.a[v]), root, 1, 1, &d);
        }
    }
    
}

void collect_hc_reverse_links(hc_links* link, ma_ug_t* ug, bubble_type* bub)
{
    uint64_t i, j, k, d = 0, m, pre;
    uint32_t beg, sink, n, v, *a = NULL;
    kvec_t_u32_warp stack, result;
    hc_edge *e = NULL;
    kv_init(stack.a); kv_init(result.a);
    for (i = 0; i < bub->num.n-1; i++)
    {
        get_bubbles(bub, i, &beg, &sink, &a, &n, NULL);
        for (k = 0; k < n; k++)
        {
            v = a[k]>>1;
            for (j = 0; j < link->a.a[v].f.n; j++)
            {
                if(link->a.a[v].f.a[j].del) continue;
                e = get_hc_edge(link, link->a.a[v].f.a[j].uID, v, 1);
                if(e == NULL) fprintf(stderr, "ERROR\n");
                e->del = 1;
            }
            link->a.a[v].f.n = 0;
        }

        v = beg>>1; 
        if(bub->index[v] > bub->num.n)
        {
            for (j = 0; j < link->a.a[v].f.n; j++)
            {
                if(link->a.a[v].f.a[j].del) continue;
                e = get_hc_edge(link, link->a.a[v].f.a[j].uID, v, 1);
                if(e == NULL) fprintf(stderr, "ERROR\n");
                e->del = 1;
            }
            link->a.a[v].f.n = 0;
        }

        v = sink>>1;
        if(bub->index[v] > bub->num.n)
        {
            for (j = 0; j < link->a.a[v].f.n; j++)
            {
                if(link->a.a[v].f.a[j].del) continue;
                e = get_hc_edge(link, link->a.a[v].f.a[j].uID, v, 1);
                if(e == NULL) fprintf(stderr, "ERROR\n");
                e->del = 1;
            }
            link->a.a[v].f.n = 0;
        }
    }

    for (i = 0; i < bub->num.n-1; i++)
    {
        get_bubbles(bub, i, &beg, &sink, &a, &n, NULL);
        if(n == 2)
        {
            push_hc_edge(&(link->a.a[a[0]>>1]), a[1]>>1, 1, 1, &d);
            push_hc_edge(&(link->a.a[a[1]>>1]), a[0]>>1, 1, 1, &d);
            continue;
        }
        beg = beg>>1; sink = sink>>1;
        for (k = 0; k < n; k++)
        {
            v = a[k]>>1;
            dfs_bubble(ug->g, &stack, &result, v, beg, sink);
            set_reverse_links(a, n, &result, v, link);
        }
    }
    kv_destroy(stack.a); kv_destroy(result.a);

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
                if(link->a.a[i].f.a[k].dis == 0) link->a.a[i].f.a[m-1].dis = 0;
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

void write_hc_links(hc_links* link, kvec_pe_hit* hits, const char *fn)
{
    uint64_t k;
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
	sprintf(buf, "%s.hic.lk", fn);
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


    fwrite(&hits->a.n, sizeof(hits->a.n), 1, fp);
    fwrite(hits->a.a, sizeof(pe_hit), hits->a.n, fp);

    fclose(fp);
    free(buf);
}


void write_hc_hits(kvec_pe_hit* hits, const char *fn)
{
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
	sprintf(buf, "%s.hic.lk", fn);
    FILE* fp = fopen(buf, "w");

    fwrite(&hits->a.n, sizeof(hits->a.n), 1, fp);
    fwrite(hits->a.a, sizeof(pe_hit), hits->a.n, fp);

    fclose(fp);
    free(buf);
}

int load_hc_links(hc_links* link, kvec_pe_hit* hits, const char *fn)
{
    uint64_t k, flag = 0;
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
	sprintf(buf, "%s.hic.lk", fn);

    FILE* fp = NULL; 
    fp = fopen(buf, "r"); 
    if(!fp) return 0;

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


    kv_init(hits->a);
    flag += fread(&hits->a.n, sizeof(hits->a.n), 1, fp);
    hits->a.m = hits->a.n; MALLOC(hits->a.a, hits->a.n);
    flag += fread(hits->a.a, sizeof(pe_hit), hits->a.n, fp);


    fclose(fp);
    free(buf);
    fprintf(stderr, "[M::%s::] ==> Hi-C linkages have been loaded\n", __func__);
    return 1;
}

int load_hc_hits(kvec_pe_hit* hits, const char *fn)
{
    uint64_t flag = 0;
    char *buf = (char*)calloc(strlen(fn) + 25, 1);
	sprintf(buf, "%s.hic.lk", fn);

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

void print_hc_links(hc_links* link, int dir, H_partition* hap)
{
    uint64_t i, k;
    if(dir == 0)
    {
        for (i = 0; i < link->a.n; ++i) 
        { 
            for (k = 0; k < link->a.a[i].e.n; k++)
            {
                if(link->a.a[i].e.a[k].del) continue;
                fprintf(stderr, "s-utg%.6d(%c)\tCLU:%u:%d\td-utg%.6d(%c)\tCLU:%u:%d\t%lu\t%c\t%f\te\n", 
                (int)(i+1), "01"[!!(link->a.a[i].e.a[k].dis&(uint64_t)2)],
                hap->hap[i]>>3, get_phase_status(hap, i),
                (int)(link->a.a[i].e.a[k].uID+1), "01"[!!(link->a.a[i].e.a[k].dis&(uint64_t)1)],
                hap->hap[link->a.a[i].e.a[k].uID]>>3, get_phase_status(hap, link->a.a[i].e.a[k].uID),
                link->a.a[i].e.a[k].dis == (uint64_t)-1? (uint64_t)-1 : link->a.a[i].e.a[k].dis>>3,
                "fb"[!!(link->a.a[i].e.a[k].dis&(uint64_t)4)], link->a.a[i].e.a[k].weight);    
            }
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
                
                if((x->rGraph.a[i].a[k].weight == 0) || 
                   (bub->index[x->rGraph.a[i].a[k].uID] > bub->num.n) || 
                   (bub->index[i] > bub->num.n) || (x->rGraph.a[i].a[k].del))
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
            if(bub->index[u] > bub->num.n)
            {
                if(d < k) kdq_push(uint64_t, x->q, set_dv(u, d+1));
            }
            else
            {
                kdq_push(uint64_t, x->q, set_dv(u , d));
                b_mer_d = d;
                if(bub->index[u] < bub->num.n && x->bmerVis.a[u] == 0)
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
            if(bub->index[u] > bub->num.n)
            {
                if(d < k) kdq_push(uint64_t, x->q, set_dv(u, d+1));
            }
            else
            {
                kdq_push(uint64_t, x->q, set_dv(u , d));
                b_mer_d = d;
                if(bub->index[u] < bub->num.n && x->bmerVis.a[u] == 0)
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

        if(bub->index[v>>1] < bub->num.n)
        {
            if(bub_extend && x->bmerVis.a[v>>1] == 0)
            {
                get_bubbles((bubble_type*)bub, bub->index[v>>1], &beg, &sink, &a, &n, NULL);
                for (j = 0; j < n; j++) x->bmerVis.a[(a[j]>>1)] = 1;
            }
            x->bmerVis.a[v>>1] = 1;
        }
         
        if(bub->index[v>>1] == bub->num.n && bub_only == 0) x->bmerVis.a[v>>1] = 1;

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
    if(bub->index[src] >= bub->num.n) return;
    get_bubbles((bubble_type*)bub, bub->index[src], &beg, &sink, &a, &n, NULL);
    memset(x->bmerVis.a, 0, x->bmerVis.n);
    ///select_bmer(src, k, bub, x, 1);
    select_bmer_distance(beg^1, k, bub, x, 1, 1);
    select_bmer_distance(sink^1, k, bub, x, 1, 1);
    /*******************************for debug************************************/
    // uint64_t i;
    // for (i = 0; i < x->bmerVis.n; ++i) 
    // { 
    //     if(x->bmerVis.a[i] == 0) continue;
    //     fprintf(stderr, "(k)utg%.6dl\n", (int)(i+1));    
    // }
    // for (i = 0; i < x->bmerVis.n; ++i) 
    // { 
    //     fprintf(stderr, "(label)utg%.6dl: %u, (num)%u\n", (int)(i+1), bub->index[i], bub->num.n);    
    // }
    /*******************************for debug************************************/
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
            if((*bub_id) < bub->num.n)
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

    for (i = 0; i < bub->num.n-1; i++)
    {
        phase_bubble(i, &b, bub, flag, ug, m, link, x);
    }

    free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
    free(flag);
    fprintf(stderr, "[M::%s::%.3f]\n", __func__, yak_realtime()-index_time);
    ///print_phased_bubble(x, bub, ug->g->n_seq);

    return x;
}


void debug_hc_links(ha_ug_index* idx, hc_links* link, sldat_t* sl, bubble_type* bub, const char *fn1)
{
    ///print_hits(idx, &sl->hits, fn1);
    destory_hc_links(link);
    init_hc_links(link, sl->idx->ug->g->n_seq, R_INF.total_reads);
    collect_hc_links(sl->idx, &sl->hits, link, bub);
    ///print_hc_links(link, 0);
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
    #define OFFSET_RATE 0.000000001
    #define OFFSET_RATE_THRES 20.7232658359
    long double rate = get_trans(idx, x);
    if(rate < 0) rate = 0;
    rate += OFFSET_RATE;
    if(rate > 0.5) rate = 0.5;

    double w = log((1-rate)/rate);
    if(w < 0) w = 0;
    if(w > OFFSET_RATE_THRES) w = OFFSET_RATE_THRES;
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



void weight_edges(ha_ug_index* idx, kvec_pe_hit* hits, hc_links* link, bubble_type* bub)
{
    uint64_t k, i, shif = 64 - idx->uID_bits, beg, end, t_d;
    hc_edge *e1 = NULL, *e2 = NULL;
    double weight;

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
        if(bub->index[beg] > bub->num.n) continue;
        if(bub->index[end] > bub->num.n) continue;

        t_d = get_hic_distance(&(hits->a.a[k]), link, idx);
        if(t_d == (uint64_t)-1) continue;

        e1 = get_hc_edge(link, beg, end, 0);
        e2 = get_hc_edge(link, end, beg, 0);
        if(e1 == NULL || e2 == NULL) continue;
        weight = get_trans_weight(idx, t_d);

        e1->weight += weight;
        e2->weight += weight;
    }
}

void init_hic_p(ha_ug_index* idx, kvec_pe_hit* hits, hc_links* link, bubble_type* bub)
{
    uint64_t k, i, m, beg, end, t_d, r_idx, f_idx, b_size, med = (uint64_t)-1, uID;
    uint32_t b_beg, b_end, n, *a, b_cnt;
    kvec_t(uint64_t) buf, buf_idx;
    kv_init(buf);
    kv_init(buf_idx);

    buf.n = 0;
    for (i = 0; i < bub->num.n-1; i++)
    {
        get_bubbles(bub, i, &b_beg, &b_end, &a, &n, &b_size);
        for (k = b_cnt = 0; k < n; k++)
        {
            b_cnt +=bub->ug->u.a[(a[k]>>1)].n;
        }
        ///too small
        if(b_cnt <= 3) continue;
        kv_push(uint64_t, buf, b_size); 
    }
    radix_sort_hc64(buf.a, buf.a+buf.n);
    med = buf.a[(uint64_t)(buf.n*0.75)];
    // if((buf.n&1) == 1) med = buf.a[buf.n>>1];
    // if((buf.n&1) == 0) med = (buf.a[buf.n>>1] + buf.a[(buf.n>>1)-1])/2;
    
    
    

    buf.n = 0;
    for (k = 0; k < hits->a.n; ++k) 
    {
        beg = ((hits->a.a[k].s<<1)>>(64 - idx->uID_bits));
        end = ((hits->a.a[k].e<<1)>>(64 - idx->uID_bits));

        if(bub->index[beg] > bub->num.n) continue;
        if(bub->index[end] > bub->num.n) continue;

        

        if(beg == end)
        {
            t_d = get_hic_distance(&(hits->a.a[k]), link, idx);
            if(t_d == (uint64_t)-1) continue;

            t_d = (t_d << 1);
            kv_push(uint64_t, buf, t_d); 
        }
        else
        {
            if(get_hc_edge(link, beg, end, 1))
            {
                t_d = get_hic_distance(&(hits->a.a[k]), link, idx);
                if(t_d == (uint64_t)-1) continue;

                t_d = (t_d << 1) + 1; 
                kv_push(uint64_t, buf, t_d); 
            }
        }
    }

    ///might have bias, we may not use right linkage larger than trans rc linkage
    radix_sort_hc64(buf.a, buf.a+buf.n);

    for (k = 0, r_idx = f_idx = (uint64_t)-1; k < buf.n; k++)
    {
        if((buf.a[k]&1) == 0) r_idx = k;
        if((buf.a[k]&1) == 1) f_idx = k;

        ///fprintf(stderr, "%lu\t%lu\n", buf.a[k]>>1, buf.a[k]&1);
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
    for (k = cutoff - t, pre = buf_idx.a[cutoff - t - 1], t_d = 0; k < cutoff + t; k++)
    {
        t_d += (buf_idx.a[k] - pre);
        pre = buf_idx.a[k];
    }
    step = (t_d/(t * 2))*100;

    buf_idx.n = 0;
    uint64_t step_s = 0, step_e = step, cnt[2], k_end;
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
                // fprintf(stderr, "i: %u, step_s: %lu, step_e: %lu, cnt[0]: %lu, cnt[1]: %lu, rate: %f\n", 
                // (uint32_t)(buf_idx.n>>2), step_s, step_e, cnt[0], cnt[1], ((double)cnt[1])/(double)(cnt[0] + cnt[1]));
                kv_push(uint64_t, buf_idx, step_s);
                kv_push(uint64_t, buf_idx, step_e);
                kv_push(uint64_t, buf_idx, cnt[0]);
                kv_push(uint64_t, buf_idx, cnt[1]);
                step_s += step;
                step_e += step;
                cnt[0] = cnt[1] = 0;
            }
        }
    }

    if(cnt[0] > 0 || cnt[1] > 0)
    {
        kv_push(uint64_t, buf_idx, step_s);
        kv_push(uint64_t, buf_idx, step_e);
        kv_push(uint64_t, buf_idx, cnt[0]);
        kv_push(uint64_t, buf_idx, cnt[1]);
    }

    uint64_t smooth_step = 20, k_i, cnt_0;
    for (k = 0; k+smooth_step < (buf_idx.n>>2); k++)
    {
        for (k_i = cnt_0 = 0; k_i < smooth_step; k_i++)
        {
            if(buf_idx.a[((k+k_i)<<2)+2] == 0 || 
               buf_idx.a[((k+k_i)<<2)+3] == 0)
            {
                cnt_0++;
            }
        }

        if(cnt_0 >= smooth_step * 0.3)
        {
            break;
        }
    }

    for (k_end = k; k < (buf_idx.n>>2); k++)
    {
        buf_idx.a[(k_end<<2)+1] = buf_idx.a[(k<<2)+1];
        buf_idx.a[(k_end<<2)+2] += buf_idx.a[(k<<2)+2];
        buf_idx.a[(k_end<<2)+3] += buf_idx.a[(k<<2)+3];
    }
    
    buf_idx.n = (k_end+1)<<2;
    if(k_end == 0) buf_idx.n = 0;

    for (k = i = 0; k < buf_idx.n; k += 4)
    {
        if(buf_idx.a[k+2] == 0 && buf_idx.a[k+3] == 0) continue;
        buf_idx.a[i] = buf_idx.a[k];
        buf_idx.a[i+1] = buf_idx.a[k+1];
        buf_idx.a[i+2] = buf_idx.a[k+2];
        buf_idx.a[i+3] = buf_idx.a[k+3];
        i += 4;
    }

    buf_idx.n = i;

    // for (k = 0; k < buf_idx.n; k += 4)
    // {
    //     fprintf(stderr, "step_s: %lu, step_e: %lu, rate: %f\n", buf_idx.a[k], buf_idx.a[k+1], 
    //     (double)(buf_idx.a[k+3])/(double)(buf_idx.a[k+2] + buf_idx.a[k+3]));
    // }
    ///idx->step = step;
    LeastSquare(buf_idx.a, buf_idx.n, idx, med);

    fprintf(stderr, "idx->a: %f, idx->b: %f, idx->frac: %f, med: %lu\n", 
    (double)idx->a, (double)idx->b, (double)idx->frac, med);
    

    hc_edge *e = NULL;
    for (i = 0; i < link->a.n; i++)
    {
        for (k = 0; k < link->a.a[i].f.n; k++)
        {
            if(link->a.a[i].f.a[k].del) continue;
            if(link->a.a[i].f.a[k].dis != 0) continue;
            uID = link->a.a[i].f.a[k].uID;
            e = get_hc_edge(link, i, uID, 0);
            if(e) e->del = 1;
            e = get_hc_edge(link, uID, i, 0);
            if(e) e->del = 1;
        }
    }

    for (i = 0; i < link->a.n; i++)
    {
        for (k = m = 0; k < link->a.a[i].e.n; k++)
        {
            if(link->a.a[i].e.a[k].del) continue;
            if(link->a.a[i].e.a[k].dis == (uint64_t)-1) continue;

            link->a.a[i].e.a[m] = link->a.a[i].e.a[k];
            link->a.a[i].e.a[m].weight = 0;
            m++;
        }
        link->a.a[i].e.n = m;
    }

    weight_edges(idx, hits, link, bub);


    kv_destroy(buf);
    kv_destroy(buf_idx);
}

#define is_hap_set(i, Hap) (!!((Hap).hap[i]&((Hap).m[0]|(Hap).m[1])))


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

void get_related_weight(uint32_t x, H_partition* hap, double* w0, double* w1)
{
    (*w0) = (*w1) = 0;
    if(x >= hap->link->a.n) return;
    uint32_t i, a_n = hap->link->a.a[x].e.n;
    hc_edge* a = hap->link->a.a[x].e.a;
    for (i = 0; i < a_n; i++)
    {
        if(a[i].del) continue;
        if((hap->hap[a[i].uID] & hap->m[0])) (*w0)+= a[i].weight;
        if((hap->hap[a[i].uID] & hap->m[1])) (*w1)+= a[i].weight;
    }
    return;
}

void set_path_hap(bub_p_t_warp *b, uint32_t root, H_partition* hap)
{
    uint32_t v, u, label;
    double w0 = 0, w1 = 0, cur_w0, cur_w1;
    ///v is the sink of this bubble
    v = b->S.a[0];
    do {
		u = b->a[v].p; // u->v
        if(v != b->S.a[0]) 
        {
            get_related_weight(v>>1, hap, &cur_w0, &cur_w1);
            w0 += cur_w0; w1 += cur_w1;
        }
		v = u;
	} while (v != root);

    if(w0 >= w1)
    {
        label = hap->label | hap->m[0];
    }
    else
    {
        label = hap->label | hap->m[1];
    }

	v = b->S.a[0];
    do {
		u = b->a[v].p; // u->v
        if(v != b->S.a[0]) hap->hap[v>>1] |= label;
		v = u;
	} while (v != root);
}


uint64_t get_phase_path(ma_ug_t *ug, uint32_t s, uint32_t d, bub_p_t_warp *b, H_partition* hap)
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

                get_related_weight(w>>1, hap, &(t->w[0]), &(t->w[1]));
                t->w[0] += nw_0; t->w[1] += nw_1; 

                t->ac = ac + ((!is_hap_set(w>>1, *hap))?ug->u.a[(w>>1)].n : 0);
                t->uc = uc + ((is_hap_set(w>>1, *hap))?ug->u.a[(w>>1)].n : 0);
                
                ++n_pending;
            }
            else {
                to_replace = 0;

                cur_nc = nc + ug->u.a[(w>>1)].n;
                /**need fix**/
                cur_nh = nh + get_path_weight(w, v, s, b, hap->link);
                get_related_weight(w>>1, hap, &cur_w0, &cur_w1);
                cur_w0 += nw_0; cur_w1 += nw_1;
                cur_weight = cur_nh + MAX(cur_w0, cur_w1) - MIN(cur_w0, cur_w1);
                max_weight = t->nh + MAX(t->w[0], t->w[1]) - MIN(t->w[0], t->w[1]);

                cur_ac = ac + ((!is_hap_set(w>>1, *hap))?ug->u.a[(w>>1)].n : 0);;
                cur_uc = uc + ((is_hap_set(w>>1, *hap))?ug->u.a[(w>>1)].n : 0);
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
    set_path_hap(b, s, hap);
    pop_reset:

    for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
        bub_p_t *t = &b->a[b->b.a[i]];
        t->p = t->d = t->nc = t->ac = t->uc = t->r = t->s = 0;
        t->nh = t->w[0] = t->w[1] = 0;
    }

    return n_pop;
}

uint32_t get_available_com(H_partition* hap, bubble_type* bub, ma_ug_t *ug, uint32_t check_self, uint32_t check_others)
{
    hc_links* link = hap->link;
    uint32_t beg, sink, n, *a, i, j, k, uID, max_bub_i, max_non_bub_i, max_i, is_ava;
    double w, cur_w0 = 0, cur_w1 = 0, max_bub_w, max_non_bub_w;
    max_i = (uint32_t)-1;

    for (i = 0, max_bub_w = -1, max_bub_i = (uint32_t)-1; i < bub->num.n-1; i++)
    {
        get_bubbles(bub, i, &beg, &sink, &a, &n, NULL);

        for (j = 0, w = 0, is_ava = 0; j < n; j++)
        {
            uID = a[j]>>1;
            if(check_self && is_hap_set(uID, *hap)) break;
            if(check_others)
            {
                get_related_weight(uID, hap, &cur_w0, &cur_w1);
                w += (MAX(cur_w0, cur_w1) - MIN(cur_w0, cur_w1));
                for (k = 0; k < link->a.a[uID].e.n; k++)
                {
                    if(link->a.a[uID].e.a[k].del) continue;
                    if(!is_hap_set(link->a.a[uID].e.a[k].uID, *hap)) continue;
                    is_ava = 1;
                }
            }
            else
            {
                for (k = 0; k < link->a.a[uID].e.n; k++)
                {
                    if(link->a.a[uID].e.a[k].del) continue;
                    w += link->a.a[uID].e.a[k].weight;
                    is_ava = 1;
                }
            }
        }
        if(j != n) continue;
        if(is_ava == 0) continue;

        if(w > max_bub_w)
        {
            max_bub_w = w;
            max_bub_i = i;
        }
    }

    for (i = 0, max_non_bub_w = -1, max_non_bub_i = (uint32_t)-1; i < ug->u.n; i++)
    {
        if(bub->index[i] == bub->num.n)
        {
            uID = i;
            if(check_self && is_hap_set(uID, *hap)) continue;
            if(check_others)
            {
                get_related_weight(uID, hap, &cur_w0, &cur_w1);
                w += (MAX(cur_w0, cur_w1) - MIN(cur_w0, cur_w1));
                for (k = 0, w = 0, is_ava = 0; k < link->a.a[uID].e.n; k++)
                {
                    if(link->a.a[uID].e.a[k].del) continue;
                    if(!is_hap_set(link->a.a[uID].e.a[k].uID, *hap)) continue;
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
            }
        }
    }

    if(max_bub_i != (uint32_t)-1 && max_non_bub_i != (uint32_t)-1)
    {
        if(max_non_bub_w > max_bub_w)
        {
            max_i = (max_non_bub_i << 1) + 1;
            w = max_non_bub_w;
        }
        else
        {
            max_i = (max_bub_i << 1);
            w = max_bub_w;
        }
    }
    else if(max_bub_i != (uint32_t)-1)
    {
        max_i = (max_bub_i << 1);
        w = max_bub_w;
    }
    else if(max_non_bub_i != (uint32_t)-1)
    {
        max_i = (max_non_bub_i << 1) + 1;
        w = max_non_bub_w;
    }

    if(max_i == (uint32_t)-1)
    {
        for (i = 0, max_non_bub_w = -1, max_non_bub_i = (uint32_t)-1; i < ug->u.n; i++)
        {
            if(bub->index[i] < bub->num.n)
            {
                uID = i;
                if(check_self && is_hap_set(uID, *hap)) continue;
                for (k = 0, w = 0, is_ava = 0; k < link->a.a[uID].e.n; k++)
                {
                    if(link->a.a[uID].e.a[k].del) continue;
                    if(check_others && (!is_hap_set(link->a.a[uID].e.a[k].uID, *hap))) continue;
                    w += link->a.a[uID].e.a[k].weight;
                    is_ava = 1;
                }
                if(is_ava == 0) continue;

                if(w > max_non_bub_w)
                {
                    max_non_bub_w = w;
                    max_non_bub_i = i;
                }
            }
        }

        if(max_non_bub_i != (uint32_t)-1)
        {
            max_i = (max_non_bub_i << 1) + 1;
            w = max_non_bub_w;
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

uint32_t get_unset_com(H_partition* hap, bubble_type* bub, ma_ug_t *ug)
{
    uint32_t max_i = get_available_com(hap, bub, ug, 1, 1);

    if(max_i == (uint32_t)-1)
    {
        max_i = get_available_com(hap, bub, ug, 1, 0);
        if(max_i != (uint32_t)-1) hap->label += hap->label_add;
    }

    return max_i;
}

void phase_com(H_partition* hap, ma_ug_t *ug, bub_p_t_warp* b, bubble_type* bub, uint32_t bid)
{
    
    if((bid & 1) == 0) ///bubble
    {
        uint32_t beg = (uint32_t)-1, sink = (uint32_t)-1, n, *a;
        get_bubbles(bub, bid>>1, &beg, &sink, &a, &n, NULL);
        ///fprintf(stderr, "+bubble-%uth, beg: %u, sink: %u, phasing ID: %u\n", bid>>1, beg>>1, sink>>1, hap->label>>3);
        get_phase_path(ug, beg, sink, b, hap);
        get_phase_path(ug, beg, sink, b, hap);
        ///fprintf(stderr, "-bubble-%uth, beg: %u, sink: %u, phasing ID: %u\n", bid>>1, beg>>1, sink>>1, hap->label>>3);
    }
    else
    {
        double cur_w0, cur_w1;
        get_related_weight(bid>>1, hap, &cur_w0, &cur_w1);
        ///fprintf(stderr, "utg-%uth, phasing ID: %u\n", bid>>1, hap->label>>3);
        if(cur_w0 >= cur_w1)
        {
            hap->hap[bid>>1] |= (hap->label | hap->m[0]);
        }
        else
        {
            hap->hap[bid>>1] |= (hap->label | hap->m[1]);
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
            ///if(o_d < -1) fprintf(stderr, "ERROR\n");
            weight += (o_d*link->a.a[h[j]].e.a[k].weight);
        }
    }

    return weight;
}

void update_partition_flag(H_partition* hap, hc_links* link, uint32_t id)
{
    uint32_t k, *h0, h0_n, *h1, h1_n, uID, flag = 0;
    int status;
    get_phased_block(&(hap->g_p), NULL, id, NULL, NULL, &h0, &h0_n, &h1, &h1_n, NULL, NULL);

    status = hap->g_p.a[id].status[0];
    if(status == 1) flag = hap->m[0];
    else if(status == -1) flag = hap->m[1];
    else if(status == 0) flag = hap->m[2];
    else if(status == -2) flag = 0;
    for (k = 0; k < h0_n; k++)
    {
        uID = h0[k];
        hap->hap[uID] >>= 3;
        hap->hap[uID] <<= 3;
        hap->hap[uID] |= flag;
    }

    status = hap->g_p.a[id].status[1];
    if(status == 1) flag = hap->m[0];
    else if(status == -1) flag = hap->m[1];
    else if(status == 0) flag = hap->m[2];
    else if(status == -2) flag = 0;
    for (k = 0; k < h1_n; k++)
    {
        uID = h1[k];
        hap->hap[uID] >>= 3;
        hap->hap[uID] <<= 3;
        hap->hap[uID] |= flag;
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
        hap->g_p.a[i].weight[0] = hap->g_p.a[i].weight[1] = 0;


        h0_status[0] = h0_status[1] = h0_status[2] = h0_status[3] = 0;
        for (k = 0; k < h0_n; k++)
        {
            ///fprintf(stderr, "(0) utg%.6ul\n", h0[k] + 1);
            ///if((hap->g_p.index[h0[k]]>>1) != i) fprintf(stderr, "ERROR\n");
            h0_status[get_phase_status(hap, h0[k])+2]++;
        }

        hap->g_p.a[i].weight[0] = get_cluster_weight(hap, link, h0, h0_n);
        if(h1_n == 0)
        {
            if(h0_status[0] == h0_n)
            {
                hap->g_p.a[i].status[0] = -2;
            }
            else
            {
                if(h0_status[1] > 0 || h0_status[3] > 0)
                {
                    h0_status[0] = h0_status[2] = 0;

                    h0_h = -2;
                    h0_status_max = 0;
                    for (k = 0; k < 4; k++)
                    {
                        if(h0_status[k] > h0_status_max) h0_status_max = h0_status[k], h0_h = (int)(k) - 2;
                    }
                    hap->g_p.a[i].status[0] = h0_h;
                }
                else if(h0_status[2] > 0)
                {
                    hap->g_p.a[i].status[0] = 0;
                }
                else
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
                ///fprintf(stderr, "(1) utg%.6ul\n", h1[k] + 1);
                ///if((hap->g_p.index[h1[k]]>>1) != i) fprintf(stderr, "ERROR\n");
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

        update_partition_flag(hap, link, i);
    }
}

uint32_t init_contig_partition(H_partition* hap, ha_ug_index* idx, bubble_type* bub)
{
    hc_links* link = idx->link;
    ma_ug_t *ug = idx->ug;
    bub_p_t_warp b;
    memset(&b, 0, sizeof(bub_p_t_warp));
    CALLOC(b.a, ug->g->n_seq*2); 
    uint32_t i, k, k_n, nv = ug->g->n_seq * 2;
    for (i = 0; i < nv; i++)
    {
        b.a[i].w[0] = b.a[i].w[1] = b.a[i].nh = 0;
        b.a[i].p =b.a[i].d = b.a[i].nc = b.a[i].uc = b.a[i].ac = b.a[i].r = b.a[i].s = 0;
    }
    

    hap->n = ug->u.n;
    MALLOC(hap->hap, hap->n);
    memset(hap->hap, 0, hap->n*sizeof(uint32_t));
    MALLOC(hap->lock, hap->n);
    memset(hap->lock, 0, hap->n);
    hap->m[0] = 1; hap->m[1] = 2; hap->m[2] = 4;
    hap->link = link;
    hap->label = 0;
    hap->label_add = 8;
    init_G_partition(&(hap->g_p), hap->n);

    uint32_t max_i = get_available_com(hap, bub, ug, 0, 0);
    
    if(max_i == (uint32_t)-1) return 0;

    hap->label = 0;
    while (1)
    {
        max_i = get_unset_com(hap, bub, ug);
        if(max_i == (uint32_t)-1) break;
        phase_com(hap, ug, &b, bub, max_i);
    }

    for (i = 0; i < hap->n; i++)
    {
        if((hap->hap[i]&hap->m[0])&&(hap->hap[i]&hap->m[1]))
        {
            hap->hap[i] >>= 3;
            hap->hap[i] <<= 3;
            hap->hap[i] |= hap->m[2];
        }
    }

    free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);

    partition_warp* res = NULL;
    hc_edge *a = NULL;
    uint32_t a_n, v, u, uv = (uint32_t)-1, k_nv, k_nu;
    for (i = 0; i < hap->n; i++)
    {
        v = i;
        a = link->a.a[v].f.a;
        a_n = link->a.a[v].f.n;
        for (k = k_n = 0; k < a_n; k++)
        {
            if(a[k].del) continue;
            if(a[k].dis != 0) break;
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
            if(a[k].dis != 0) break;
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

        kv_pushp(partition_warp, hap->g_p, &res);
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
            if(hap->g_p.index[res->a.a[k]] != (uint32_t)-1) fprintf(stderr, "ERROR\n");
            hap->g_p.index[res->a.a[k]] = hap->g_p.n-1;
            hap->g_p.index[res->a.a[k]] = hap->g_p.index[res->a.a[k]] << 1;
        } 

        for (; k < res->a.n; k++)
        {
            if(hap->g_p.index[res->a.a[k]] != (uint32_t)-1) fprintf(stderr, "ERROR\n");
            hap->g_p.index[res->a.a[k]] = hap->g_p.n-1;
            hap->g_p.index[res->a.a[k]] = (hap->g_p.index[res->a.a[k]] << 1) + 1;
        } 
    }

    ///print_contig_partition(hap, "first");

    adjust_contig_partition(hap, link);

    ///print_contig_partition(hap, "second");

    return 1;
}


uint32_t get_max_unitig(H_partition* hap, hc_links* link, ma_ug_t *ug, bubble_type* bub)
{
    double min, weight;
    uint32_t i, min_i;
     
    for (i = 0, min = 1, min_i = (uint32_t)-1; i < hap->g_p.n; i++)
    {
        if(hap->lock[i]) continue;
        weight = 0;
        if(hap->g_p.a[i].h[0] > 0 && (hap->g_p.a[i].status[0] == 1 || hap->g_p.a[i].status[0] == -1))
        {
            weight += (hap->g_p.a[i].weight[0] * hap->g_p.a[i].status[0]);
        }

        if(hap->g_p.a[i].h[1] > 0 && (hap->g_p.a[i].status[1] == 1 || hap->g_p.a[i].status[1] == -1))
        {
            weight += (hap->g_p.a[i].weight[1] * hap->g_p.a[i].status[1]);
        }

        if(weight >= 0) continue;
        if(weight < min)
        {
            min = weight;
            min_i = i;
        }
    }

    return min_i;
}

void flip_unitig(H_partition* hap, hc_links* link, ma_ug_t *ug, bubble_type* bub, uint32_t id)
{
    if(hap->g_p.a[id].h[0] > 0 && hap->g_p.a[id].status[0] != 1 && hap->g_p.a[id].status[0] != -1) return;
    if(hap->g_p.a[id].h[1] > 0 && hap->g_p.a[id].status[1] != 1 && hap->g_p.a[id].status[1] != -1) return;
    uint32_t k, j, m, *h0, h0_n, *h1, h1_n, uID, *h = NULL, h_n;
    int status;
    double weight;
    get_phased_block(&(hap->g_p), NULL, id, NULL, NULL, &h0, &h0_n, &h1, &h1_n, NULL, NULL);
    if(h0_n > 0)
    {
        status = hap->g_p.a[id].status[0];
        h = h0; h_n = h0_n;
        for (j = 0; j < h_n; j++)
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
                weight = link->a.a[h[j]].e.a[k].weight;
                hap->g_p.a[hap->g_p.index[uID]>>1].weight[hap->g_p.index[uID]&1] -= (2*status*weight);
            }
        }
        hap->g_p.a[id].status[0] *= -1;
    }

    if(h1_n > 0)
    {
        status = hap->g_p.a[id].status[1];
        h = h1; h_n = h1_n;
        for (j = 0; j < h_n; j++)
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
                weight = link->a.a[h[j]].e.a[k].weight;
                hap->g_p.a[hap->g_p.index[uID]>>1].weight[hap->g_p.index[uID]&1] -= (2*status*weight);
            }
        }
        hap->g_p.a[id].status[1] *= -1;
    }
}


uint32_t phasing_improvement(H_partition* hap, ha_ug_index* idx, bubble_type* bub)
{   
    uint32_t i, occ = 0;
    memset(hap->lock, 0, sizeof(uint8_t)*hap->g_p.n);
    while (1)
    {
        i = get_max_unitig(hap, idx->link, idx->ug, bub);
        if(i == (uint32_t)-1) break;
        hap->lock[i] = 1;
        flip_unitig(hap, idx->link, idx->ug, bub, i);
        occ++;
    }

    for (i = 0; i < hap->g_p.n; i++)
    {
        update_partition_flag(hap, idx->link, i);
    }
    
    return !!occ;
} 

void destory_contig_partition(H_partition* hap)
{
    free(hap->lock);
    free(hap->hap);
    destory_G_partition(&(hap->g_p));
}

void label_unitigs(H_partition* hap, ma_ug_t* ug)
{
    memset(R_INF.trio_flag, AMBIGU, R_INF.total_reads * sizeof(uint8_t));
    uint32_t i, k, j, *h0, h0_n, *h1, h1_n, uID, *h = NULL, h_n, flag = AMBIGU;
    int status;
    ma_utg_t *u = NULL;

    for (i = 0; i < hap->g_p.n; i++)
    {
        if(hap->g_p.a[i].h[0] > 0 && hap->g_p.a[i].status[0] != 1 && hap->g_p.a[i].status[0] != -1) continue;
        if(hap->g_p.a[i].h[1] > 0 && hap->g_p.a[i].status[1] != 1 && hap->g_p.a[i].status[1] != -1) continue;
        get_phased_block(&(hap->g_p), NULL, i, NULL, NULL, &h0, &h0_n, &h1, &h1_n, NULL, NULL);

        status = hap->g_p.a[i].status[0];
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



        status = hap->g_p.a[i].status[1];
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

    fprintf(stderr, "# reads: %lu\n", occ);

    for (i = occ = 0; i < R_INF.total_reads; i++)
    {
        if(R_INF.trio_flag[i] == FATHER) occ++;
    }

    fprintf(stderr, "# Father reads: %lu\n", occ);

    for (i = occ = 0; i < R_INF.total_reads; i++)
    {
        if(R_INF.trio_flag[i] == MOTHER) occ++;
    }
    
    fprintf(stderr, "# Mother reads: %lu\n", occ);
}

int hic_short_align(const char *fn1, const char *fn2, ha_ug_index* idx)
{
    double index_time = yak_realtime();
    sldat_t sl;
    gzFile fp1, fp2;
    if ((fp1 = gzopen(fn1, "r")) == 0) return 0;
    if ((fp2 = gzopen(fn2, "r")) == 0) return 0;
    sl.ks1 = kseq_init(fp1);
    sl.ks2 = kseq_init(fp2);
    sl.idx = idx;
    sl.chunk_size = 20000000;
    sl.n_thread = asm_opt.thread_num;
    sl.total_base = sl.total_pair = 0;
    idx->max_cnt = 5;
    kv_init(sl.hits.a);
    fprintf(stderr, "u.n: %d, uID_bits: %lu, pos_bits: %lu\n", (uint32_t)idx->ug->u.n, idx->uID_bits, idx->pos_bits);
    
    bubble_type bub;
    identify_bubbles(idx->ug, &bub);
    ///print_bubbles(idx->ug, &bub, sl.hits.a.n?&sl.hits:NULL, idx->link, idx);

    if(!load_hc_hits(&sl.hits, asm_opt.output_file_name))
    {
        /*******************************for debug************************************/
        // load_reads(&R1, fn1);
        // test_reads(&R1, fn1);
        // load_reads(&R2, fn2);
        // test_reads(&R1, fn1);
        /*******************************for debug************************************/

        kt_pipeline(3, worker_pipeline, &sl, 3);
        /*******************************for debug************************************/
        // sort_hits(&sl.hits);
        // print_hits(idx, &sl.hits, fn1);
        /*******************************for debug************************************/
        dedup_hits(&sl.hits);
        write_hc_hits(&sl.hits, asm_opt.output_file_name);
    }

    collect_hc_links(sl.idx, &sl.hits, idx->link, &bub);
    collect_hc_reverse_links(idx->link, idx->ug, &bub);
    
    ///debug_hc_links(idx, idx->link, &sl, &bub, fn1);
    
    init_hic_p((ha_ug_index*)sl.idx, &sl.hits, idx->link, &bub);
    
    H_partition hap;
    init_contig_partition(&hap, idx, &bub);
    phasing_improvement(&hap, idx, &bub);
    label_unitigs(&hap, idx->ug);
    
    ///print_hc_links(idx->link, 0, &hap);
    ///print_contig_partition(&hap, "final");

    destory_contig_partition(&hap);

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
    kseq_destroy(sl.ks1);
    kseq_destroy(sl.ks2);
	gzclose(fp1);
    gzclose(fp2);
    
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
    hic_short_align(asm_opt.hic_files[0], asm_opt.hic_files[1], ug_index);
    
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

    // uint64_t up_dis = buf.a[(uint64_t)(buf.n*0.99)]>>1, step = 1000;
    // uint64_t step_s = 0, step_e = step, cnt[2];
    // for (k = cnt[0] = cnt[1] = 0; k < buf.n; k++)
    // {
    //     if(step_s > up_dis) step_e = (buf.a[buf.n-1]>>1) + 1;
    //     if((buf.a[k]>>1) < step_e && (buf.a[k]>>1) >= step_s)
    //     {
    //         cnt[buf.a[k]&1]++;
    //     }
    //     if((buf.a[k]>>1) >= step_e)
    //     {
    //         while (!((buf.a[k]>>1) < step_e && (buf.a[k]>>1) >= step_s))
    //         {
    //             fprintf(stderr, "i: %lu, step_s: %lu, step_e: %lu, cnt[0]: %lu, cnt[1]: %lu, rate: %f\n", 
    //             step_s/step, step_s, step_e, cnt[0], cnt[1], ((double)cnt[1])/(double)(cnt[0] + cnt[1]));
    //             step_s += step;
    //             step_e += step;
    //             cnt[0] = cnt[1] = 0;
    //         }
    //     }
    // }

    // if(cnt[0] > 0 || cnt[1] > 0)
    // {
    //     fprintf(stderr, "i: %lu, step_s: %lu, step_e: %lu, cnt[0]: %lu, cnt[1]: %lu, rate: %f\n", 
    //             step_s/step, step_s, step_e, cnt[0], cnt[1], ((double)cnt[1])/(double)(cnt[0] + cnt[1]));
    // }

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


int hic_short_align_bench(const char *fn1, const char *fn2, const char *output_file_name, ha_ug_index* idx)
{
    double index_time = yak_realtime();
    sldat_t sl;
    gzFile fp1, fp2;
    if ((fp1 = gzopen(fn1, "r")) == 0) return 0;
    if ((fp2 = gzopen(fn2, "r")) == 0) return 0;
    sl.ks1 = kseq_init(fp1);
    sl.ks2 = kseq_init(fp2);
    sl.idx = idx;
    sl.chunk_size = 20000000;
    sl.n_thread = asm_opt.thread_num;
    sl.total_base = sl.total_pair = 0;
    idx->max_cnt = 5;
    kv_init(sl.hits.a);
    fprintf(stderr, "u.n: %d, uID_bits: %lu, pos_bits: %lu\n", (uint32_t)idx->ug->u.n, idx->uID_bits, idx->pos_bits);
    
    if(!load_hc_hits(&sl.hits, output_file_name))
    {
        kt_pipeline(3, worker_pipeline, &sl, 3);
        dedup_hits(&sl.hits);
        write_hc_hits(&sl.hits, output_file_name);
    }
    bench_idx bench;
    init_bench_idx(&bench, idx->read_g, idx->ug);
    ///print_bench_idx(&bench, idx->ug);
    evaluate_bench_idx(&bench, &sl.hits, idx->ug);

    destory_bench_idx(&bench);
    kv_destroy(sl.hits.a);
    kseq_destroy(sl.ks1);
    kseq_destroy(sl.ks2);
	gzclose(fp1);
    gzclose(fp2);
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

    hic_short_align_bench(asm_opt.hic_files[0], asm_opt.hic_files[1], output_file_name, ug_index);

    free(output_file_name);
}