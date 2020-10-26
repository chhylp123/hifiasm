#define __STDC_LIMIT_MACROS
#include "hic.h"
#include "htab.h"
#include "assert.h"
#include "Overlaps.h"
#include "khashl.h"

#define generic_hc_key(x) (x)
KHASHL_MAP_INIT(static klib_unused, hc_pt_t, hc_pt, uint64_t, uint64_t, generic_hc_key, kh_eq_generic)

typedef struct {
	hc_pt_t *h;
	uint64_t n;
	uint64_t *a;
    khint_t end;
} hc_pt1_t;

typedef struct {
	ma_ug_t* ug;
    uint64_t uID_bits;
    uint64_t uID_mode;
    uint64_t pos_bits;
    uint64_t pos_mode;
    uint64_t rev_mode;
    hc_pt1_t idx;
    uint64_t k;
} ha_ug_index;

ha_ug_index* ug_index;

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

inline uint64_t hc_hash_long(uint64_t x[4], uint64_t* skip)
{
	///compare forward k-mer and reverse complementary strand
	(*skip) = get_k_direction(x);
    if((*skip) == (uint64_t)-1) return (*skip);
	return yak_hash64_64(x[(*skip)<<1|0]) + yak_hash64_64(x[(*skip)<<1|1]);
}

inline uint64_t get_hc_pt1_count(ha_ug_index* idx, uint64_t key, uint64_t** pos_list)
{
    uint64_t beg;
    khint_t k;
    k = hc_pt_get(idx->idx.h, key);
    if (k == kh_end(idx->idx.h))
    {
        return 0;
    }
    beg = kh_val(idx->idx.h, k);
    if(pos_list) *pos_list = idx->idx.a + beg;
    if(k == idx->idx.end) return idx->idx.n - beg;
    for (k++; k != kh_end(idx->idx.h); ++k)
    {
        if (kh_exist(idx->idx.h, k)) 
        {
            return kh_val(idx->idx.h, k) - beg;
        }
    }
    return idx->idx.n - beg;
}

void count_hc_pt1(char* seq, uint64_t len, ha_ug_index* idx)
{
    uint64_t i, l;
    khint_t key;
    int absent;
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
                hash = hc_hash_long(x, &skip);
                if(skip == (uint64_t)-1) continue;
                key = hc_pt_put(idx->idx.h, hash, &absent);
                if(absent) kh_val(idx->idx.h, key) = 0;
                kh_val(idx->idx.h, key)++;
            }
                
        } else l = 0, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
    }
}

void fill_hc_pt1(char* seq, uint64_t len, uint64_t uID, ha_ug_index* idx)
{
    uint64_t i, l, pos, *pos_list = NULL, cnt;
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
                hash = hc_hash_long(x, &skip);
                if(skip == (uint64_t)-1) continue;
                pos = (skip << 63) | ((uID << (64-idx->uID_bits))>>1) | (i & idx->pos_mode);
                cnt = get_hc_pt1_count(idx, hash, &pos_list);
                ///assert(cnt != 0);
                if(pos_list[cnt-1]!=cnt-1)
                {
                    pos_list[pos_list[cnt-1]]=pos;
                    pos_list[cnt-1]++;
                }
                else
                {
                    pos_list[pos_list[cnt-1]]=pos;
                }
            }
        } else l = 0, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
    }
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
                hash = hc_hash_long(x, &skip);
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

void test_unitig_index(ha_ug_index* idx)
{
    uint32_t i;
    ma_utg_t *u = NULL;
    for (i = 0; i < idx->ug->u.n; i++)
    {
        ///fprintf(stderr, "i: %u, n: %u\n", i, (uint32_t)idx->ug->u.n);
        u = &(idx->ug->u.a[i]);
        if(u->m == 0) continue;
        test_hc_pt1(u->s, u->len, i, idx);
    }
    for (i = 0; i < idx->idx.n; i++)
    {
        if(idx->idx.a[i]!=(uint64_t)-1)
        {
            fprintf(stderr, "ERROR i\n");
        }
    }
}

void test_count1(char* seq, uint64_t len, uint64_t uID, ha_ug_index* idx, uint64_t direct_cnt)
{
    uint64_t i, l, cnt;
    khint_t k;
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
                hash = hc_hash_long(x, &skip);
                if(skip == (uint64_t)-1) continue;
                if(!direct_cnt)
                {
                    cnt = get_hc_pt1_count(idx, hash, NULL);
                }
                else
                {
                    k = hc_pt_get(idx->idx.h, hash);
                    if (k == kh_end(idx->idx.h))
                    {
                        cnt = 0;
                    }
                    else
                    {
                        cnt = kh_val(idx->idx.h, k);
                    }
                    
                }
            }
        } else l = 0, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
    }
}

void hc_pt_t_gen(hc_pt1_t* pt)
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
    fwrite(&idx->idx.n, sizeof(idx->idx.n), 1, fp);
    fwrite(&idx->idx.end, sizeof(idx->idx.end), 1, fp);
    fwrite(idx->idx.a, sizeof(uint64_t), idx->idx.n, fp);
    hc_pt_save(idx->idx.h, fp);


    fprintf(stderr, "[M::%s] Index has been written.\n", __func__);
    free(gfa_name);
    fclose(fp);
    return 1;
}


int load_hc_pt_index(ha_ug_index** r_idx, char* file_name)
{
    char* gfa_name = (char*)malloc(strlen(file_name)+25);
    sprintf(gfa_name, "%s.hc_tlb", file_name);
    FILE* fp = fopen(gfa_name, "r");
    if (!fp) {
        free(gfa_name);
        return 0;
    }
    ha_ug_index* idx = NULL; CALLOC(idx, 1);

    fread(&idx->uID_bits, sizeof(idx->uID_bits), 1, fp);
    fread(&idx->uID_mode, sizeof(idx->uID_mode), 1, fp);
    fread(&idx->pos_bits, sizeof(idx->pos_bits), 1, fp);
    fread(&idx->pos_mode, sizeof(idx->pos_mode), 1, fp);
    fread(&idx->rev_mode, sizeof(idx->rev_mode), 1, fp);
    fread(&idx->k, sizeof(idx->k), 1, fp);
    fread(&idx->idx.n, sizeof(idx->idx.n), 1, fp);
    fread(&idx->idx.end, sizeof(idx->idx.end), 1, fp);
    MALLOC(idx->idx.a, idx->idx.n);
    fread(idx->idx.a, sizeof(uint64_t), idx->idx.n, fp);
    hc_pt_load(&(idx->idx.h), fp);
    (*r_idx) = idx;

    fprintf(stderr, "[M::%s] Index has been loaded.\n", __func__);
    free(gfa_name);
    fclose(fp);
    return 1;
}


ha_ug_index* build_unitig_index(ma_ug_t *ug, int k)
{
    uint32_t i;
    ma_utg_t *u = NULL;
    ha_ug_index* idx = NULL; CALLOC(idx, 1);
    double index_time = yak_realtime(), beg_time;
    for (idx->uID_bits=1; (uint64_t)(1<<idx->uID_bits)<(uint64_t)ug->u.n; idx->uID_bits++);
    idx->pos_bits = 64 - idx->uID_bits - 1;
    idx->uID_mode = (((uint64_t)-1) << (64-idx->uID_bits))>>1;
    idx->pos_mode = ((uint64_t)-1) >> (64-idx->pos_bits);
    idx->rev_mode = ((uint64_t)1) << 63;
    idx->ug = ug;
    idx->k = k;
    idx->idx.h = hc_pt_init();

    beg_time = yak_realtime();
    for (i = 0; i < idx->ug->u.n; i++)
    {
        u = &(idx->ug->u.a[i]);
        if(u->m == 0) continue;
        count_hc_pt1(u->s, u->len, idx);
    }
    fprintf(stderr, "[M::%s::%.3f] ==> Counting\n", __func__, yak_realtime()-beg_time);
    // beg_time = yak_realtime();
    // for (i = 0; i < idx->ug->u.n; i++)
    // {
    //     u = &(idx->ug->u.a[i]);
    //     if(u->m == 0) continue;
    //     test_count1(u->s, u->len, i, idx, 1);
    // }
    // fprintf(stderr, "[M::%s::%.3f] ==> Counting 1\n", __func__, yak_realtime()-beg_time);


    beg_time = yak_realtime();
    hc_pt_t_gen(&(idx->idx));
    fprintf(stderr, "[M::%s::%.3f] ==> Memory allocating\n", __func__, yak_realtime()-beg_time);


    // beg_time = yak_realtime();
    // for (i = 0; i < idx->ug->u.n; i++)
    // {
    //     u = &(idx->ug->u.a[i]);
    //     if(u->m == 0) continue;
    //     test_count1(u->s, u->len, i, idx, 0);
    // }
    // fprintf(stderr, "[M::%s::%.3f] ==> Counting 2\n", __func__, yak_realtime()-beg_time);


    beg_time = yak_realtime();
    for (i = 0; i < idx->ug->u.n; i++)
    {
        u = &(idx->ug->u.a[i]);
        if(u->m == 0) continue;
        fill_hc_pt1(u->s, u->len, i, idx);
    }
    fprintf(stderr, "[M::%s::%.3f] ==> Filling pos\n", __func__, yak_realtime()-beg_time);

    fprintf(stderr, "[M::%s::%.3f] ==> HiC index has been built\n", __func__, yak_realtime()-index_time);

    return idx;
}

void destory_hc_pt_index(ha_ug_index* r_idx)
{
    if(r_idx->idx.h) hc_pt_destroy(r_idx->idx.h);
    if(r_idx->idx.a) free(r_idx->idx.a);
}

void hic_analysis(ma_ug_t *ug)
{
    ug_index = NULL;
    ///int exist = 0;//load_hc_pt_index(&ug_index, asm_opt.output_file_name);
    int exist = load_hc_pt_index(&ug_index, asm_opt.output_file_name);
    if(exist == 0) ug_index = build_unitig_index(ug, asm_opt.hic_mer_length);
    if(exist == 0) write_hc_pt_index(ug_index, asm_opt.output_file_name);
    ug_index->ug = ug;
    
    test_unitig_index(ug_index);
    
    destory_hc_pt_index(ug_index);
}