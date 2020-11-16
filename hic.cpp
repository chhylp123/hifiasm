#define __STDC_LIMIT_MACROS
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
#define kdq_clear(q) ((q)->count = (q)->front = 0)
#define kv_malloc(v, s) ((v).n = 0, (v).m = (s), MALLOC((v).a, (s)))

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
    int weight;
    uint32_t uID:31, del:1;
    uint32_t enzyme;
} hc_edge;

typedef struct{
    kvec_t(hc_edge) e;
    kvec_t(hc_edge) f;//forbiden
} hc_linkeage;

typedef struct{
    kvec_t(hc_linkeage) a;
} hc_links;

typedef struct{
    ///kvec_t(hc_edge) a;
    size_t n, m; 
    hc_edge *a;
}hc_edge_warp;

typedef struct{
    kvec_t(hc_edge_warp) rGraph;
    kvec_t(uint64_t) order;
    kvec_t(uint8_t) rGraphSet;
    kvec_t(uint8_t) rGraphVis;
    kvec_t(uint8_t) utgVis;
    kvec_t(uint8_t) bmerVis;
    kdq_t(uint64_t) *q;
    kvec_t(uint64_t) parent;
    uint64_t uID_mode, uID_shift, n, src, dest;
    int p_mer, a_mer, b_mer;
} min_cut_t;

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
} bubble_type;

typedef struct {
	ma_ug_t* ug;
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

#define pe_hit_an1_key(x) ((x).s<<1)
KRADIX_SORT_INIT(pe_hit_an1, pe_hit, pe_hit_an1_key, 8)
#define pe_hit_an2_key(x) ((x).e<<1)
KRADIX_SORT_INIT(pe_hit_an2, pe_hit, pe_hit_an2_key, 8)
#define generic_key(x) (x)
KRADIX_SORT_INIT(hc64, uint64_t, generic_key, 8)

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
    if(uID1 == uID2) return;
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
        if (k == hits->a.n || (hits->a.a[k].s<<1) != (hits->a.a[l].s<<1)) 
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
}

inline void get_bubbles(bubble_type* bub, uint64_t id, uint32_t* beg, uint32_t* sink, uint32_t** a, uint32_t* n)
{
    (*a) = bub->list.a + bub->num.a[id] + 2;
    (*n) = bub->num.a[id+1] - bub->num.a[id] - 2;
    (*beg) = bub->list.a[bub->num.a[id]];
    (*sink) = bub->list.a[bub->num.a[id] + 1];
}

void identify_bubbles(ma_ug_t* ug, bubble_type* bub)
{
    asg_cleanup(ug->g);
    if (!ug->g->is_symm) asg_symm(ug->g);
    memset(bub, 0, sizeof(bubble_type));
    uint32_t v, n_vtx = ug->g->n_seq * 2, tLen, i, mode = (((uint32_t)-1)<<2);
    bub->ug = ug; 
    CALLOC(bub->index, n_vtx);
    for (i = 0; i < ug->g->n_seq; i++)
    {
        if(ug->g->seq[i].c == 1)
        {
            ug->g->seq[i].c = 0;
            bub->index[i] = 4;
        } 
    }
    kv_init(bub->list); kv_init(bub->num);
    buf_t b; memset(&b, 0, sizeof(buf_t)); b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    for (i = 0, tLen = 1; i < ug->u.n; i++) tLen += ug->u.a[i].len;
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(ug->g->seq[v>>1].del) continue;
        if(asg_arc_n(ug->g, v) < 2) continue;
        if((bub->index[v]&(uint32_t)3) != 0) continue;
        if(asg_bub_pop1_primary_trio(ug->g, NULL, v, tLen, &b, (uint32_t)-1, (uint32_t)-1, 0))
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
            bub->index[b.S.a[0]^1] &= mode;; bub->index[b.S.a[0]^1] += 3;
        }
    }

    
    for (v = 0; v < n_vtx; ++v) 
    {
        if((bub->index[v]&(uint32_t)3) !=2) continue;
        kv_push(uint32_t, bub->num, bub->list.n);
        if(asg_bub_pop1_primary_trio(ug->g, NULL, v, tLen, &b, (uint32_t)-1, (uint32_t)-1, 0))
        {
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
            bub->index[i] = bub->num.n;
        }         
    }

    uint32_t beg, sink, n, *a;
    for (i = 0; i < bub->num.n-1; i++)
    {
        get_bubbles(bub, i, &beg, &sink, &a, &n);
        for (v = 0; v < n; v++)
        {
            bub->index[(a[v]>>1)] = i;
        }

        bub->index[(beg>>1)] = bub->index[(sink>>1)] = (uint32_t)-1;
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
        get_bubbles(bub, i, &beg, &sink, &a, &n);
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


    // fprintf(stderr, "************bubble utgs************\n");
    // for (i = 0, tLen = 0, t_utg = 0; i < bub->num.n-1; i++)
    // {
    //     get_bubbles(bub, i, &beg, &sink, &a, &n);
    //     t_utg += n;
    //     fprintf(stderr, "(%lu)\tbeg:utg%.6u\tsink:utg%.6u\n", i, (beg>>1)+1, (sink>>1)+1);
    //     for (k = 0; k < n; k++)
    //     {
    //         tLen +=bub->ug->u.a[(a[k]>>1)].len;
    //         fprintf(stderr, "utg%.6u,", (a[k]>>1)+1);
    //     }
    //     fprintf(stderr, "\n");
    // }

    // fprintf(stderr, "************het utgs************\n");
    // for (i = 0; i < ug->g->n_seq; i++)
    // {
    //     if(bub->index[i] == bub->num.n) fprintf(stderr, "utg%.6lu\n", i+1);
    // }
    // fprintf(stderr, "************het utgs************\n");
}




void init_hc_links(hc_links* link, uint64_t ug_num)
{
    kv_init(link->a);
    kv_resize(hc_linkeage, link->a, ug_num);
    link->a.n = ug_num;
    uint64_t i;
    for (i = 0; i < link->a.n; i++)
    {
        kv_init(link->a.a[i].e);
        kv_init(link->a.a[i].f);
    }
}

void destory_hc_links(hc_links* link)
{
    uint64_t i;
    for (i = 0; i < link->a.n; i++)
    {
        kv_destroy(link->a.a[i].e);
        kv_destroy(link->a.a[i].f);
    }
    kv_destroy(link->a);
}

void push_hc_edge(hc_linkeage* x, uint64_t uID, int weight, int dir)
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
    
    p->del = p->enzyme = 0;
    p->uID = uID;
    p->weight = weight;
}


void collect_hc_links(const ha_ug_index* idx, kvec_pe_hit* hits, hc_links* link)
{
    uint64_t k, shif = 64 - idx->uID_bits, beg, end;
    for (k = 0; k < hits->a.n; ++k) 
    {
        beg = ((hits->a.a[k].s<<1)>>shif);
        end = ((hits->a.a[k].e<<1)>>shif);
        push_hc_edge(&(link->a.a[beg]), end, 1, 0);
        push_hc_edge(&(link->a.a[end]), beg, 1, 0);
    }
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
    uint64_t i, k;
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
            push_hc_edge(&(link->a.a[root]), v, 1, 1);
            push_hc_edge(&(link->a.a[v]), root, 1, 1);
        }
    }
    
}

void collect_hc_reverse_links(hc_links* link, ma_ug_t* ug, bubble_type* bub)
{
    uint64_t i, k;
    uint32_t beg, sink, n, v, *a = NULL;
    kvec_t_u32_warp stack, result;
    kv_init(stack.a); kv_init(result.a);
    for (i = 0; i < bub->num.n-1; i++)
    {
        get_bubbles(bub, i, &beg, &sink, &a, &n);
        if(n == 2)
        {
            push_hc_edge(&(link->a.a[a[0]>>1]), a[1]>>1, 1, 1);
            push_hc_edge(&(link->a.a[a[1]>>1]), a[0]>>1, 1, 1);
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
}

void write_hc_links(hc_links* link, const char *fn)
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
    
    fclose(fp);
    free(buf);
}

int load_hc_links(hc_links* link, const char *fn)
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


    fclose(fp);
    free(buf);
    return 1;
}

void print_hc_links(hc_links* link)
{
    uint64_t i, k;
    for (i = 0; i < link->a.n; ++i) 
    { 
        for (k = 0; k < link->a.a[i].e.n; k++)
        {
            if(link->a.a[i].e.a[k].del) continue;
            fprintf(stderr, "utg%.6d\tutg%.6d\t%d\t+\n", 
            (int)(i+1), (int)(link->a.a[i].e.a[k].uID+1), link->a.a[i].e.a[k].weight);    
        }
    }


    for (i = 0; i < link->a.n; ++i) 
    { 
        for (k = 0; k < link->a.a[i].f.n; k++)
        {
            if(link->a.a[i].f.a[k].del) continue;
            fprintf(stderr, "utg%.6d\tutg%.6d\t%d\t-", 
            (int)(i+1), (int)(link->a.a[i].f.a[k].uID+1), link->a.a[i].f.a[k].weight);
            if(link->a.a[i].f.a[k].weight > 2) fprintf(stderr,"\tcomplex");  
            fprintf(stderr,"\n");  
        }
    }
}

void init_min_cut_t(min_cut_t* x, hc_links* link, const bubble_type* bub, const ma_ug_t *ug)
{
    uint64_t utg_num = link->a.n, i, k;
    x->n = utg_num;

    kv_malloc(x->rGraphSet, utg_num); x->rGraphSet.n = utg_num;
    ///must utg_num<<1)
    kv_malloc(x->rGraphVis, utg_num); x->rGraphVis.n = utg_num;
    kv_malloc(x->utgVis, utg_num); x->utgVis.n = utg_num;
    kv_malloc(x->bmerVis, utg_num); x->bmerVis.n = utg_num;
    kv_malloc(x->order, utg_num); x->order.n = utg_num;
    kv_malloc(x->parent, utg_num); x->parent.n = utg_num;
    ///uresolved BUGs, if use kv_resize segfault; if use kv_malloc, work?????
    kv_malloc(x->rGraph, utg_num); x->rGraph.n = utg_num;
    // kv_init(x->rGraph); kv_resize(hc_edge_warp, x->rGraph, utg_num); x->rGraph.n = utg_num;
    
    //must be utg_num + 2 since we may need to add fake nodes
    for (i = 1; (uint64_t)(1<<i) < (utg_num + 2); i++);
    x->uID_mode = ((uint64_t)-1) >> (64-i); 
    x->uID_shift = i;
    for (i = 0; i < utg_num; i++)
    {
        ///x->order.a[i] = link->a.a[i].f.n;
        x->order.a[i] = ug->u.a[i].len;
        x->order.a[i] <<= x->uID_shift; 
        x->order.a[i] |= (uint64_t)(i & x->uID_mode);

        x->rGraphSet.a[i] = 0;
        x->rGraphVis.a[i] = 0;
        x->utgVis.a[i] = 0;
        x->bmerVis.a[i] = 0;
        x->parent.a[i] = (uint64_t)-1;
        
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
                if((bub->index[x->rGraph.a[i].a[k].uID] > bub->num.n) || 
                   (bub->index[i] > bub->num.n))
                {
                    x->rGraph.a[i].a[k].del = 1;
                }
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
    kv_destroy(x->rGraphSet);
    kv_destroy(x->rGraphVis);
    kv_destroy(x->utgVis);
    kv_destroy(x->bmerVis);
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
    x->bmerVis.n = x->parent.n = x->order.n = x->rGraph.n = x->rGraphVis.n = x->rGraphSet.n = link->a.n;
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

uint64_t add_mul_convex(min_cut_t* x, uint64_t* a, uint64_t n)
{
    if(n == 0) return (uint64_t)-1;
    if(n == 1) return a[0];
    kv_push(uint8_t, x->rGraphSet, 0);
    kv_push(uint8_t, x->rGraphVis, 0);
    kv_push(uint8_t, x->bmerVis, 0);
    kv_push(uint64_t, x->parent, 0);
    kv_resize(hc_edge_warp, x->rGraph, x->rGraph.n+1); 
    kv_init(x->rGraph.a[x->rGraph.n]);
    uint64_t i, k;
    hc_edge t;
    for (i = 0; i < n; i++)
    {
        t.uID = a[i]; t.del = t.enzyme = t.weight = 0;
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
    if(link->a.a[uID].f.n == 1)
    {
        (*dest) = link->a.a[uID].f.a[0].uID;
        return;
    }

    uint64_t i;
    for (i = 0; i < link->a.a[uID].f.n; i++)
    {
        if(link->a.a[uID].f.a[i].del) continue;
        kv_push(uint64_t, buff->a, link->a.a[uID].f.a[i].uID);
    }
    (*dest) = add_mul_convex(x, buff->a.a, buff->a.n);
}


uint64_t bfs_flow(uint64_t src, uint64_t dest, min_cut_t* x, kvec_t_u64_warp* buff)
{
    uint64_t *p = NULL, v, u, i;
    if(dest != (uint64_t)-1) memset(x->rGraphVis.a, 0, x->rGraphVis.n);

    kdq_push(uint64_t, x->q, src); 
    if(buff) kv_push(uint64_t, buff->a, src);

    x->rGraphVis.a[src] = 1;
    x->parent.a[src] = (uint64_t)-1;

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

            ///x->parent.a[u] = v;
            x->parent.a[u] = x->rGraph.a[v].a[i].weight;
            x->parent.a[u] <<= x->uID_shift;
            x->parent.a[u] |= v;

            kdq_push(uint64_t, x->q, u); 
            if(buff) kv_push(uint64_t, buff->a, u);
            ///set u or v to be 1? doesn't matter
            x->rGraphVis.a[u] = 1;
        }
    }

    return 0;
}

hc_edge* get_rGraph_edge(min_cut_t* x, uint64_t src, uint64_t dest)
{
    if(src >= x->rGraph.n) return NULL;
    uint64_t i;
    for (i = 0; i < x->rGraph.a[src].n; i++)
    {
        if(x->rGraph.a[src].a[i].uID == dest) return &(x->rGraph.a[src].a[i]);
    }
    return NULL;
}

uint64_t maxFlow(uint64_t src, uint64_t dest, min_cut_t* x)
{
    uint64_t flow = 0, max_flow = 0, v, u;
    hc_edge *p;

    while (bfs_flow(src, dest, x, NULL))
    {
        kdq_clear(x->q);
        flow = (uint64_t)-1;
        for (v = dest; v != src; v = x->parent.a[v]&x->uID_mode)
        {
            flow = MIN(flow, (x->parent.a[v]>>x->uID_shift));
        }

        fprintf(stderr, "***********flow: %lu*********\n", flow);

        for (v = dest; v != src; v = x->parent.a[v]&x->uID_mode)
        {
            u = x->parent.a[v]&x->uID_mode;
            p = get_rGraph_edge(x, u, v);

            fprintf(stderr, "utg%.6lul (%d)\n", u+1, p->weight);


            p->weight -= flow;
            p = get_rGraph_edge(x, v, u);
            p->weight += flow;
            x->rGraphSet.a[u] = x->rGraphSet.a[v] = 1;
        }

        max_flow += flow;
    }
    
    return max_flow;
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

void graph_cut(uint64_t src, uint64_t dest, min_cut_t* x)
{
    if(maxFlow(src, dest, x))
    {
        ///in the last time bfs of maxFlow, rGraphVis has already been set
        uint64_t i, j, v, u;
        hc_edge *p;
        /*******************************for debug************************************/
        print_src_dest(src, x, "src utg:");
        print_src_dest(dest, x, "dest utg:");
        /*******************************for debug************************************/
        for (i = 0; i < x->rGraphVis.n; i++)
        {
            if(x->rGraphVis.a[i] == 0) continue;
            v = i;
            for (j = 0; j < x->rGraph.a[i].n; j++)
            {
                if(x->rGraph.a[i].a[j].del) continue;
                u = x->rGraph.a[i].a[j].uID;
                if(x->rGraphVis.a[u]) continue;
                if(!x->bmerVis.a[u]) continue;
                /*******************************for debug************************************/
                fprintf(stderr, "utg%.6lul\tutg%.6lul\t%d\n", v+1, u+1, x->rGraph.a[i].a[j].weight);    
                /*******************************for debug************************************/
                ///delete <v, u>
                x->rGraph.a[i].a[j].del = 1;
                ///delete <u, v>
                p = get_rGraph_edge(x, u, v);
                p->del = 1;
            }
            
        }
    }
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
            fprintf(stderr, "src(utg%.6dl)\tdes(utg%.6dl)\t%d\n", 
            (int)(i+1), (int)(x->rGraph.a[i].a[k].uID+1), x->rGraph.a[i].a[k].weight);    
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



uint64_t select_bmer(uint32_t src, uint64_t k, const bubble_type* bub, min_cut_t* x)
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
                x->bmerVis.a[u] = 1;
                kdq_push(uint64_t, x->q, set_dv(u , d));
                b_mer_d = d;
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
                    get_bubbles((bubble_type*)bub, bub->index[u], &beg, &sink, &a, &n);
                    for (j = 0; j < n; j++) x->bmerVis.a[(a[j]>>1)] = 1;
                }
                //must be here
                x->bmerVis.a[u] = 1;
            }            
        }
    }

    return b_mer_d;
}

void get_bmer_unitgs(min_cut_t* x, const bubble_type* bub, uint64_t k, uint64_t src)
{
    memset(x->bmerVis.a, 0, x->bmerVis.n);
    select_bmer(src, k, bub, x);
    /*******************************for debug************************************/
    uint64_t i;
    for (i = 0; i < x->bmerVis.n; ++i) 
    { 
        if(x->bmerVis.a[i] == 0) continue;
        fprintf(stderr, "(k)utg%.6dl\n", (int)(i+1));    
    }
    // for (i = 0; i < x->bmerVis.n; ++i) 
    // { 
    //     fprintf(stderr, "(label)utg%.6dl: %u, (num)%u\n", (int)(i+1), bub->index[i], bub->num.n);    
    // }
    /*******************************for debug************************************/
}

void clean_hap(hc_links* link, bubble_type* bub, const ma_ug_t *ug)
{
    min_cut_t x;
    kvec_t_u64_warp buff;
    kv_init(buff.a);
    init_min_cut_t(&x, link, (const bubble_type*)bub, ug);

    // get_Connected_Components(&x);
    // check_connective(&x, link);
    // print_rGraph(&x);
    long long i;
    uint64_t k, uID;
    ///for (i = x.n - 1; i >= 0; i--)
    for (i = 0; (uint64_t)i < x.n; i++)
    {
        ///fprintf(stderr, "i: %lu\n", i);
        uID = x.order.a[i] & x.uID_mode;
        ///fprintf(stderr, "uID: %lu, f.n: %lu\n", uID, (uint64_t)link->a.a[uID].f.n);
        if(link->a.a[uID].f.n == 0) continue;
        for (k = 0; k < link->a.a[uID].f.n; k++)
        {
            if(link->a.a[uID].f.a[k].del) continue;
            if(x.utgVis.a[link->a.a[uID].f.a[k].uID] == 0) break;
        }
        ///fprintf(stderr, "k: %lu\n", k);
        if(k == link->a.a[uID].f.n) continue;
        reset_min_cut_t(&x, link);
        ///fprintf(stderr, "reset\n");
        get_s_t(&x, link, uID, &(x.src), &(x.dest), &buff);
        ///fprintf(stderr, "x.src: %lu, x.dest: %lu\n", x.src, x.dest);
        ///Note: should only consider edges betweem bubbles, ignore edges to homo untigs

        /*******************************for debug************************************/
        if(!select_large_node(ug, &x, x.src, x.dest, 10, 10)) continue; 


        get_bmer_unitgs(&x, bub, x.b_mer, x.src);
        x.bmerVis.a[x.src] = x.bmerVis.a[x.dest] = 1;
        /*******************************for debug************************************/

        graph_cut(x.src, x.dest, &x);
        ///fprintf(stderr, "graph_cut\n");
        x.utgVis.a[uID] = 1;
        exit(0);
    }
    destory_min_cut_t(&x);
    kv_destroy(buff.a);
}

int hic_short_align(const char *fn1, const char *fn2, ha_ug_index* idx)
{
    double index_time = yak_realtime();
    hc_links link;
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
    
    if(!load_hc_links(&link, asm_opt.output_file_name))
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

        init_hc_links(&link, sl.idx->ug->g->n_seq);
        collect_hc_links(sl.idx, &sl.hits, &link);
        write_hc_links(&link, asm_opt.output_file_name);
    }
    
    /*******************************for debug************************************/
    // destory_reads(&R1);
    // destory_reads(&R2);
    /*******************************for debug************************************/
    bubble_type bub;
    identify_bubbles(idx->ug, &bub);
    print_bubbles(idx->ug, &bub, sl.hits.a.n?&sl.hits:NULL, &link, idx);
    collect_hc_reverse_links(&link, idx->ug, &bub);
    /*******************************for debug************************************/
    ///print_hc_links(&link);
    /*******************************for debug************************************/
    clean_hap(&link, &bub, idx->ug);

    kv_destroy(sl.hits.a);
    destory_hc_links(&link);
    destory_bubbles(&bub);
    kseq_destroy(sl.ks1);
    kseq_destroy(sl.ks2);
	gzclose(fp1);
    gzclose(fp2);
    
    fprintf(stderr, "[M::%s::%.3f] processed %lu pairs; %lu bases\n", __func__, yak_realtime()-index_time, sl.total_pair, sl.total_base);
    return 1;
}

void hic_analysis(ma_ug_t *ug)
{
    
    ug_index = NULL;
    int exist = load_hc_pt_index(&ug_index, asm_opt.output_file_name);
    if(exist == 0) ug_index = build_unitig_index(ug, asm_opt.hic_mer_length);
    if(exist == 0) write_hc_pt_index(ug_index, asm_opt.output_file_name);
    ug_index->ug = ug;
    ///test_unitig_index(ug_index, ug);
    hic_short_align(asm_opt.hic_files[0], asm_opt.hic_files[1], ug_index);
    
    destory_hc_pt_index(ug_index);
}