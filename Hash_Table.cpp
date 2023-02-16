#define __STDC_LIMIT_MACROS
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include "Hash_Table.h"
#include "ksort.h"
#include "kalloc.h"
pthread_mutex_t output_mutex;

#define overlap_region_key(a) ((a).y_id)
KRADIX_SORT_INIT(overlap_region_sort, overlap_region, overlap_region_key, member_size(overlap_region, y_id))

#define generic_key(x) (x)
KRADIX_SORT_INIT(hc64i, int64_t, generic_key, 8)

#define oreg_sss_lt(a, b) ((a).shared_seed > (b).shared_seed) // in the decending order
KSORT_INIT(or_sss, overlap_region, oreg_sss_lt)

#define normal_w(x, y) ((x)>=(y)?(x)/(y):1)

void overlap_region_sort_y_id(overlap_region *a, long long n)
{
	radix_sort_overlap_region_sort(a, a + n);
}

void init_overlap_region_alloc(overlap_region_alloc* list)
{
    list->size = 1000;
    list->length = 0;
    ///list->list = (overlap_region*)malloc(sizeof(overlap_region)*list->size);
    list->list = (overlap_region*)calloc(list->size, sizeof(overlap_region));
    uint64_t i;
    for (i = 0; i < list->size; i++) { 
        init_fake_cigar(&(list->list[i].f_cigar));
        init_window_list_alloc(&(list->list[i].w_list));
        init_window_list_alloc(&(list->list[i].boundary_cigars));
    }
}

void clear_overlap_region_alloc(overlap_region_alloc* list)
{
    list->length = 0;
    list->mapped_overlaps_length = 0;
    uint64_t i = 0;
    for (i = 0; i < list->size; i++) {   
        clear_fake_cigar(&(list->list[i].f_cigar));
        clear_window_list_alloc(&(list->list[i].w_list));
        clear_window_list_alloc(&(list->list[i].boundary_cigars));
    }
}

void destory_overlap_region_alloc(overlap_region_alloc* list)
{
    uint64_t i = 0;
    for (i = 0; i < list->size; i++) {
        destory_fake_cigar(&(list->list[i].f_cigar));
        destory_window_list_alloc(&(list->list[i].w_list));
        destory_window_list_alloc(&(list->list[i].boundary_cigars));
    }
    free(list->list);
}

int get_fake_gap_pos(Fake_Cigar* x, int index)
{
    return (x->buffer[index]>>32);
}

int get_fake_gap_shift(Fake_Cigar* x, int index)
{
    uint32_t tmp = ((uint32_t)(x->buffer[index]));
    int result;
    if(tmp & ((uint32_t)1))
    {
        tmp = tmp >> 1;
        result = tmp;
        result = result * -1;
    }
    else
    {
        tmp = tmp >> 1;
        result = tmp;
    }
    
    return result;
}

void gen_fake_cigar(Fake_Cigar* z, overlap_region *o, int64_t apend_be, k_mer_hit* hit, int64_t n_hit)
{
    int64_t k, dq, dr, dd, pdd; z->length = 0;
    if(apend_be == 1) add_fake_cigar(z, o->x_pos_s, 0, NULL);
    for (k = 0, pdd = INT32_MAX; k < n_hit; k++) {
        dq = hit[k].self_offset - o->x_pos_s; 
        dr = hit[k].offset - o->y_pos_s; 
        dd = dr - dq;
        // if(print) {
        //     fprintf(stderr, "[M::k->%lu] x::%u, y::%u, cnt::%u, dd::%ld, pdd::%ld, z->n::%u\n", 
        //         k, hit[k].self_offset, hit[k].offset, hit[k].cnt&(0xffu), dd, pdd, z->length);
        // }
        if(dd != pdd) {
            pdd = dd;
            add_fake_cigar(z, hit[k].self_offset, pdd, NULL);
        }
    }
    
    if((apend_be == 1) && (get_fake_gap_pos(z, z->length-1)!=((int64_t)o->x_pos_e))) {
        add_fake_cigar(z, o->x_pos_e, get_fake_gap_shift(z, z->length-1), NULL);
    }
}

void debug_cigar(Fake_Cigar* z, overlap_region *o, int64_t apend_be, k_mer_hit* hit, uint64_t n_hit)
{
    gen_fake_cigar(z, o, apend_be, hit, n_hit);
    if(!((z->length==o->f_cigar.length) && 
                    (!memcmp(z->buffer, o->f_cigar.buffer, sizeof((*(o->f_cigar.buffer)))*o->f_cigar.length)))) {
        uint64_t k;
        fprintf(stderr, "\n[M::%s] z->n::%u, o->n::%u, rev::%u\n", __func__, z->length, o->f_cigar.length, o->y_pos_strand);
        for (k = 0; k < z->length; k++) {
            fprintf(stderr, "[z::k->%lu] pos::%d, off::%d\n", k, 
                                    get_fake_gap_pos(z, k), get_fake_gap_shift(z, k));
        }
        for (k = 0; k < o->f_cigar.length; k++) {
            fprintf(stderr, "[o::k->%lu] pos::%d, off::%d\n", k, 
                                    get_fake_gap_pos(&(o->f_cigar), k), get_fake_gap_shift(&(o->f_cigar), k));
        }
        // gen_fake_cigar(z, o, apend_be, hit, n_hit, 1);
        // for (k = 0; k < o->f_cigar.length; k++) {
        // 	fprintf(stderr, "[M::k->%lu] x::%u, y::%u, cnt::%u\n", 
		// 			k, hit[k].self_offset, hit[k].offset, hit[k].cnt&(0xffu));
        // }
        
    }
}

///for backup
int ovlp_chain_gen(overlap_region_alloc* ol, overlap_region* t, int64_t xl, int64_t yl, int64_t apend_be, k_mer_hit* hit, int64_t n_hit)
{
    if (ol->length + 1 > ol->size) {
        uint64_t sl = ol->size;
        ol->size = ol->length + 1;
        kroundup64(ol->size);
        REALLOC(ol->list, ol->size);
        /// need to set new space to be 0
        memset(ol->list + sl, 0, sizeof(overlap_region)*(ol->size - sl));
    }

    if ((ol->length!=0) && (ol->list[ol->length-1].y_id==t->y_id)) {    
        if((ol->list[ol->length-1].shared_seed > t->shared_seed) ||
           ((ol->list[ol->length-1].shared_seed == t->shared_seed) && 
           (ol->list[ol->length-1].overlapLen <= t->overlapLen))) {
            return 0;
        } else {
            ol->length--;
        }
    }

    int64_t xr, yr; 
    if(t->x_pos_s <= t->y_pos_s) {
        t->y_pos_s -= t->x_pos_s; t->x_pos_s = 0;
    } else {
        t->x_pos_s -= t->y_pos_s; t->y_pos_s = 0;
    }

    xr = xl-t->x_pos_e-1; yr = yl-t->y_pos_e-1;
    if(xr <= yr) {
        t->x_pos_e = xl-1; t->y_pos_e += xr;        
    } else {
        t->y_pos_e = yl-1; t->x_pos_e += yr; 
    }

    overlap_region *o = &(ol->list[ol->length++]);
    o->shared_seed = t->shared_seed;
    o->align_length = 0;
    o->is_match = 0;
    o->non_homopolymer_errors = 0;
    o->strong = 0;
    o->x_id = t->x_id;
    o->y_id = t->y_id;
    o->x_pos_strand = 0;///always 0
    o->y_pos_strand = t->x_pos_strand;

    if (t->x_pos_strand == 1) {
        o->x_pos_e = xl-t->x_pos_s-1; o->x_pos_s = xl-t->x_pos_e-1;
        o->y_pos_e = yl-t->y_pos_s-1; o->y_pos_s = yl-t->y_pos_e-1;
    } else {
        o->x_pos_e = t->x_pos_e; o->x_pos_s = t->x_pos_s;
        o->y_pos_e = t->y_pos_e; o->y_pos_s = t->y_pos_s;
    }
    ///debug
    // debug_cigar(&(t->f_cigar), o, apend_be, hit, n_hit);

    return 1;
}

int ovlp_chain_qgen(overlap_region_alloc* ol, overlap_region* t, int64_t xl, int64_t yl, int64_t apend_be, k_mer_hit* hit, int64_t n_hit)
{
    if (ol->length + 1 > ol->size) {
        uint64_t sl = ol->size;
        ol->size = ol->length + 1;
        kroundup64(ol->size);
        REALLOC(ol->list, ol->size);
        /// need to set new space to be 0
        memset(ol->list + sl, 0, sizeof(overlap_region)*(ol->size - sl));
    }

    if ((ol->length!=0) && (ol->list[ol->length-1].y_id==t->y_id)) {    
        if((ol->list[ol->length-1].shared_seed > t->shared_seed) ||
           ((ol->list[ol->length-1].shared_seed == t->shared_seed) && 
           (ol->list[ol->length-1].overlapLen <= t->overlapLen))) {
            return 0;
        } else {
            ol->length--;
        }
    }

    int64_t xr, yr; 
    if(t->x_pos_s <= t->y_pos_s) {
        t->y_pos_s -= t->x_pos_s; t->x_pos_s = 0;
    } else {
        t->x_pos_s -= t->y_pos_s; t->y_pos_s = 0;
    }

    xr = xl-t->x_pos_e-1; yr = yl-t->y_pos_e-1;
    if(xr <= yr) {
        t->x_pos_e = xl-1; t->y_pos_e += xr;        
    } else {
        t->y_pos_e = yl-1; t->x_pos_e += yr; 
    }

    overlap_region *o = &(ol->list[ol->length++]);
    o->shared_seed = t->shared_seed;
    o->align_length = 0;
    o->is_match = 0;
    o->non_homopolymer_errors = 0;
    o->strong = 0;
    o->x_id = t->x_id;
    o->y_id = t->y_id;
    o->x_pos_strand = 0;///always 0
    o->y_pos_strand = t->x_pos_strand;

    o->x_pos_e = t->x_pos_e; o->x_pos_s = t->x_pos_s;
    o->y_pos_e = t->y_pos_e; o->y_pos_s = t->y_pos_s;
    ///debug
    // debug_cigar(&(t->f_cigar), o, apend_be, hit, n_hit);

    return 1;
}

int ovlp_chain_gen_fcigar(overlap_region_alloc* ol, overlap_region* t, int64_t xl, int64_t yl, int64_t apend_be, k_mer_hit* hit, int64_t n_hit)
{
    if (ol->length + 1 > ol->size) {
        uint64_t sl = ol->size;
        ol->size = ol->length + 1;
        kroundup64(ol->size);
        REALLOC(ol->list, ol->size);
        /// need to set new space to be 0
        memset(ol->list + sl, 0, sizeof(overlap_region)*(ol->size - sl));
    }

    if ((ol->length!=0) && (ol->list[ol->length-1].y_id==t->y_id)) {    
        if((ol->list[ol->length-1].shared_seed > t->shared_seed) ||
           ((ol->list[ol->length-1].shared_seed == t->shared_seed) && 
           (ol->list[ol->length-1].overlapLen <= t->overlapLen))) {
            return 0;
        } else {
            ol->length--;
        }
    }

    int64_t xr, yr, dd, pdd, dq, dr, id, i, fn; 
    if(t->x_pos_s <= t->y_pos_s) {
        t->y_pos_s -= t->x_pos_s; t->x_pos_s = 0;
    } else {
        t->x_pos_s -= t->y_pos_s; t->y_pos_s = 0;
    }

    xr = xl-t->x_pos_e-1; yr = yl-t->y_pos_e-1;
    if(xr <= yr) {
        t->x_pos_e = xl-1; t->y_pos_e += xr;        
    } else {
        t->y_pos_e = yl-1; t->x_pos_e += yr; 
    }

    overlap_region *o = &(ol->list[ol->length++]);
    o->shared_seed = t->shared_seed;
    o->align_length = 0;
    o->is_match = 0;
    o->non_homopolymer_errors = 0;
    o->strong = 0;
    o->x_id = t->x_id;
    o->y_id = t->y_id;
    o->x_pos_strand = 0;///always 0
    o->y_pos_strand = t->x_pos_strand;

    resize_fake_cigar(&(o->f_cigar), (t->f_cigar.length + 2), NULL);
    if(apend_be == 1) {
        add_fake_cigar(&(o->f_cigar), ((t->x_pos_strand)?(xl-t->x_pos_e-1):(t->x_pos_s)), 0, NULL);
    }

    if (t->x_pos_strand == 1) {
        o->x_pos_e = xl-t->x_pos_s-1; o->x_pos_s = xl-t->x_pos_e-1;
        o->y_pos_e = yl-t->y_pos_s-1; o->y_pos_s = yl-t->y_pos_e-1;
        
        pdd = INT32_MAX; fn = t->f_cigar.length;
        for (i = 0; i < fn; i++) {
            dd = get_fake_gap_shift(&(t->f_cigar), i);
            if(dd != pdd) {
                pdd = dd;
                add_fake_cigar(&(o->f_cigar), xl-get_fake_gap_pos(&(t->f_cigar), i)-1, pdd, NULL);
            }
        }
    } else {
        o->x_pos_e = t->x_pos_e; o->x_pos_s = t->x_pos_s;
        o->y_pos_e = t->y_pos_e; o->y_pos_s = t->y_pos_s;
    
        dq = t->x_pos_e - t->x_pos_s;
        dr = t->y_pos_e - t->y_pos_s;
        id = dr - dq;///indel size from left

        pdd = INT32_MAX; fn = t->f_cigar.length;
        for (i = fn-1; i >= 0; i--) {
            dd = get_fake_gap_shift(&(t->f_cigar), i);
            if(dd != pdd) {
                pdd = dd;
                add_fake_cigar(&(o->f_cigar), get_fake_gap_pos(&(t->f_cigar), i), id - pdd, NULL);
            }
        }
    }

    if((apend_be == 1) && (get_fake_gap_pos(&(o->f_cigar), o->f_cigar.length-1) != ((int64_t)o->x_pos_e))) {
        add_fake_cigar(&(o->f_cigar), o->x_pos_e, get_fake_gap_shift(&(o->f_cigar), o->f_cigar.length-1), NULL);
    }

    ///debug
    debug_cigar(&(t->f_cigar), o, apend_be, hit, n_hit);

    return 1;
}

int append_inexact_overlap_region_alloc(overlap_region_alloc* list, overlap_region* tmp, 
                                        long long xLen, long long yLen, int add_beg_end, void *km)
{
   
    if (list->length + 1 > list->size)
    {
        uint64_t sl = list->size;
        list->size = list->length + 1;
        kroundup64(list->size);
        if(!km) {
            REALLOC(list->list, list->size);
        } else {
            KREALLOC(km, list->list, list->size);
        }
        /// need to set new space to be 0
        memset(list->list + sl, 0, sizeof(overlap_region)*(list->size - sl));
    }

    if (list->length!=0 && list->list[list->length - 1].y_id==tmp->y_id)
    {    
        ///if(list->list[list->length - 1].shared_seed >= tmp->shared_seed)
        if((list->list[list->length - 1].shared_seed > tmp->shared_seed)
           ||
           ((list->list[list->length - 1].shared_seed == tmp->shared_seed) && 
           (list->list[list->length - 1].overlapLen <= tmp->overlapLen)))
        {
            return 0;
        }
        else
        {
            list->length--;
        }
    }

    if(tmp->x_pos_s <= tmp->y_pos_s)
    {
        tmp->y_pos_s = tmp->y_pos_s - tmp->x_pos_s;
        tmp->x_pos_s = 0;
    }
    else
    {
        tmp->x_pos_s = tmp->x_pos_s - tmp->y_pos_s;
        tmp->y_pos_s = 0;
    }


    long long x_right_length = xLen - tmp->x_pos_e - 1;
    long long y_right_length = yLen - tmp->y_pos_e - 1;

    if(x_right_length <= y_right_length)
    {
        tmp->x_pos_e = xLen - 1;
        tmp->y_pos_e = tmp->y_pos_e + x_right_length;        
    }
    else
    {
        tmp->x_pos_e = tmp->x_pos_e + y_right_length;
        tmp->y_pos_e = yLen - 1;
    }
    
    if (tmp->x_pos_strand == 1)
    {
        list->list[list->length].x_id = tmp->x_id;
        list->list[list->length].x_pos_e = xLen - tmp->x_pos_s - 1;
        list->list[list->length].x_pos_s = xLen - tmp->x_pos_e - 1;
        list->list[list->length].x_pos_strand = 0;

        list->list[list->length].y_id = tmp->y_id;
        list->list[list->length].y_pos_e = yLen - tmp->y_pos_s - 1;
        list->list[list->length].y_pos_s = yLen - tmp->y_pos_e - 1;
        list->list[list->length].y_pos_strand = 1;

        resize_fake_cigar(&(list->list[list->length].f_cigar), (tmp->f_cigar.length + 2), km);
        if(add_beg_end == 1)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), list->list[list->length].x_pos_s, 0, km);
        }
        
        long long distance_gap;
        /****************************may have bugs********************************/
        ///long long pre_distance_gap = 0;
        long long pre_distance_gap = 0xfffffffffffffff;
        /****************************may have bugs********************************/
        long long i = 0;
        for (i = 0; i < (long long)tmp->f_cigar.length; i++)
        {
            distance_gap = get_fake_gap_shift(&(tmp->f_cigar), i);
            if(distance_gap != pre_distance_gap)
            {
                pre_distance_gap = distance_gap;
                add_fake_cigar(&(list->list[list->length].f_cigar), xLen - get_fake_gap_pos(&(tmp->f_cigar), i) - 1, 
                pre_distance_gap, km);
            }
        }

        if(add_beg_end == 1 && get_fake_gap_pos(&(list->list[list->length].f_cigar), 
            list->list[list->length].f_cigar.length - 1) != (long long)list->list[list->length].x_pos_e)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), 
            list->list[list->length].x_pos_e, 
            get_fake_gap_shift(&(list->list[list->length].f_cigar), 
            list->list[list->length].f_cigar.length - 1), km);
        }
    }
    else
    {
        list->list[list->length].x_id = tmp->x_id;
        list->list[list->length].x_pos_e = tmp->x_pos_e;
        list->list[list->length].x_pos_s = tmp->x_pos_s;
        list->list[list->length].x_pos_strand = tmp->x_pos_strand;

        list->list[list->length].y_id = tmp->y_id;
        list->list[list->length].y_pos_e = tmp->y_pos_e;
        list->list[list->length].y_pos_s = tmp->y_pos_s;
        list->list[list->length].y_pos_strand = tmp->y_pos_strand;



        resize_fake_cigar(&(list->list[list->length].f_cigar), (tmp->f_cigar.length + 2), km);
        if(add_beg_end == 1)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), list->list[list->length].x_pos_s, 0, km);
        }
        
        long long distance_self_pos = tmp->x_pos_e - tmp->x_pos_s;
        long long distance_pos = tmp->y_pos_e - tmp->y_pos_s;
        long long init_distance_gap = distance_pos - distance_self_pos;
        /****************************may have bugs********************************/
        ///long long pre_distance_gap = init_distance_gap;
        long long pre_distance_gap = 0xfffffffffffffff;
        /****************************may have bugs********************************/
        long long distance_gap;
        long long i = 0;
        for (i = tmp->f_cigar.length - 1; i >= 0; i--)
        {
            distance_gap = get_fake_gap_shift(&(tmp->f_cigar), i);
            if(distance_gap != pre_distance_gap)
            {
                pre_distance_gap = distance_gap;

                add_fake_cigar(&(list->list[list->length].f_cigar), 
                get_fake_gap_pos(&(tmp->f_cigar), i), init_distance_gap - pre_distance_gap, km);
            }
        }

        if(add_beg_end == 1 && get_fake_gap_pos(&(list->list[list->length].f_cigar), 
        list->list[list->length].f_cigar.length - 1) != (long long)list->list[list->length].x_pos_e)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), 
            list->list[list->length].x_pos_e, 
            get_fake_gap_shift(&(list->list[list->length].f_cigar), 
            list->list[list->length].f_cigar.length - 1), km);
        } 
    }

    list->list[list->length].shared_seed = tmp->shared_seed;
    list->list[list->length].align_length = 0;
    list->list[list->length].is_match = 0;
    list->list[list->length].non_homopolymer_errors = 0;
    list->list[list->length].strong = 0;

    list->length++;

    return 1;
}


int append_utg_inexact_overlap_region_alloc(overlap_region_alloc* list, overlap_region* tmp, 
                                        ma_utg_v *ua, int add_beg_end, void *km)
{
    if (list->length + 1 > list->size)
    {
        uint64_t sl = list->size;
        list->size = list->length + 1;
        kroundup64(list->size);
        if(!km) {
            REALLOC(list->list, list->size);
        } else {
            KREALLOC(km, list->list, list->size);
        }
        /// need to set new space to be 0
        memset(list->list + sl, 0, sizeof(overlap_region)*(list->size - sl));
    }

    if (list->length!=0 && list->list[list->length - 1].y_id==tmp->y_id)
    {    
        ///if(list->list[list->length - 1].shared_seed >= tmp->shared_seed)
        if((list->list[list->length - 1].shared_seed > tmp->shared_seed)
           ||
           ((list->list[list->length - 1].shared_seed == tmp->shared_seed) && 
           (list->list[list->length - 1].overlapLen <= tmp->overlapLen)))
        {
            return 0;
        }
        else
        {
            list->length--;
        }
    }

    if(tmp->x_pos_s <= tmp->y_pos_s)
    {
        tmp->y_pos_s = tmp->y_pos_s - tmp->x_pos_s;
        tmp->x_pos_s = 0;
    }
    else
    {
        tmp->x_pos_s = tmp->x_pos_s - tmp->y_pos_s;
        tmp->y_pos_s = 0;
    }


    long long x_right_length = ua->a[tmp->x_id].len - tmp->x_pos_e - 1;
    long long y_right_length = ua->a[tmp->y_id].len - tmp->y_pos_e - 1;

    if(x_right_length <= y_right_length)
    {
        tmp->x_pos_e = ua->a[tmp->x_id].len - 1;
        tmp->y_pos_e = tmp->y_pos_e + x_right_length;        
    }
    else
    {
        tmp->x_pos_e = tmp->x_pos_e + y_right_length;
        tmp->y_pos_e = ua->a[tmp->y_id].len - 1;
    }
    
    if (tmp->x_pos_strand == 1)
    {
        list->list[list->length].x_id = tmp->x_id;
        list->list[list->length].x_pos_e = ua->a[tmp->x_id].len - tmp->x_pos_s - 1;
        list->list[list->length].x_pos_s = ua->a[tmp->x_id].len - tmp->x_pos_e - 1;
        list->list[list->length].x_pos_strand = 0;

        list->list[list->length].y_id = tmp->y_id;
        list->list[list->length].y_pos_e = ua->a[tmp->y_id].len - tmp->y_pos_s - 1;
        list->list[list->length].y_pos_s = ua->a[tmp->y_id].len - tmp->y_pos_e - 1;
        list->list[list->length].y_pos_strand = 1;

        resize_fake_cigar(&(list->list[list->length].f_cigar), (tmp->f_cigar.length + 2), km);
        if(add_beg_end == 1)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), list->list[list->length].x_pos_s, 0, km);
        }
        
        long long distance_gap;
        /****************************may have bugs********************************/
        ///long long pre_distance_gap = 0;
        long long pre_distance_gap = 0xfffffffffffffff;
        /****************************may have bugs********************************/
        long long i = 0;
        for (i = 0; i < (long long)tmp->f_cigar.length; i++)
        {
            distance_gap = get_fake_gap_shift(&(tmp->f_cigar), i);
            if(distance_gap != pre_distance_gap)
            {
                pre_distance_gap = distance_gap;
                add_fake_cigar(&(list->list[list->length].f_cigar), 
                ua->a[tmp->x_id].len - get_fake_gap_pos(&(tmp->f_cigar), i) - 1, 
                pre_distance_gap, km);
            }
        }

        if(add_beg_end == 1 && get_fake_gap_pos(&(list->list[list->length].f_cigar), 
            list->list[list->length].f_cigar.length - 1) != (long long)list->list[list->length].x_pos_e)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), 
            list->list[list->length].x_pos_e, 
            get_fake_gap_shift(&(list->list[list->length].f_cigar), 
            list->list[list->length].f_cigar.length - 1), km);
        }
    }
    else
    {
        list->list[list->length].x_id = tmp->x_id;
        list->list[list->length].x_pos_e = tmp->x_pos_e;
        list->list[list->length].x_pos_s = tmp->x_pos_s;
        list->list[list->length].x_pos_strand = tmp->x_pos_strand;

        list->list[list->length].y_id = tmp->y_id;
        list->list[list->length].y_pos_e = tmp->y_pos_e;
        list->list[list->length].y_pos_s = tmp->y_pos_s;
        list->list[list->length].y_pos_strand = tmp->y_pos_strand;



        resize_fake_cigar(&(list->list[list->length].f_cigar), (tmp->f_cigar.length + 2), km);
        if(add_beg_end == 1)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), list->list[list->length].x_pos_s, 0, km);
        }
        
        long long distance_self_pos = tmp->x_pos_e - tmp->x_pos_s;
        long long distance_pos = tmp->y_pos_e - tmp->y_pos_s;
        long long init_distance_gap = distance_pos - distance_self_pos;
        /****************************may have bugs********************************/
        ///long long pre_distance_gap = init_distance_gap;
        long long pre_distance_gap = 0xfffffffffffffff;
        /****************************may have bugs********************************/
        long long distance_gap;
        long long i = 0;
        for (i = tmp->f_cigar.length - 1; i >= 0; i--)
        {
            distance_gap = get_fake_gap_shift(&(tmp->f_cigar), i);
            if(distance_gap != pre_distance_gap)
            {
                pre_distance_gap = distance_gap;

                add_fake_cigar(&(list->list[list->length].f_cigar), 
                get_fake_gap_pos(&(tmp->f_cigar), i), init_distance_gap - pre_distance_gap, km);
            }
        }

        if(add_beg_end == 1 && get_fake_gap_pos(&(list->list[list->length].f_cigar), 
        list->list[list->length].f_cigar.length - 1) != (long long)list->list[list->length].x_pos_e)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), 
            list->list[list->length].x_pos_e, 
            get_fake_gap_shift(&(list->list[list->length].f_cigar), 
            list->list[list->length].f_cigar.length - 1), km);
        } 
    }

    list->list[list->length].shared_seed = tmp->shared_seed;
    list->list[list->length].align_length = 0;
    list->list[list->length].is_match = 0;
    list->list[list->length].non_homopolymer_errors = 0;
    list->list[list->length].strong = 0;

    list->length++;

    return 1;
}

void append_overlap_region_alloc_debug(overlap_region_alloc* list, overlap_region* tmp)
{
   
    if (list->length + 1 > list->size)
    {
        list->size = list->size * 2;
        list->list = (overlap_region*)realloc(list->list, sizeof(overlap_region)*list->size);
    }

    
    list->list[list->length].x_id = tmp->x_id;
    list->list[list->length].x_pos_e = tmp->x_pos_e;
    list->list[list->length].x_pos_s = tmp->x_pos_s;
    list->list[list->length].x_pos_strand = tmp->x_pos_strand;

    list->list[list->length].y_id = tmp->y_id;
    list->list[list->length].y_pos_e = tmp->y_pos_e;
    list->list[list->length].y_pos_s = tmp->y_pos_s;
    list->list[list->length].y_pos_strand = tmp->y_pos_strand;

    list->list[list->length].shared_seed = tmp->shared_seed;

    list->length++;
}

void debug_chain(k_mer_hit* a, long long a_n, Chain_Data* dp)
{
    long long i, j, current_j;
    long long selfLen, indels;
    long long distance_self_pos, distance_pos, distance_gap;
    for (i = 0; i < a_n; ++i) 
    {
        selfLen = indels = 0;
        j = i;
        while (j >= 0)
        {
            current_j = j;

            j = dp->pre[j];

            if(j != -1)
            {
                distance_self_pos = a[current_j].self_offset - a[j].self_offset;
                distance_pos = a[current_j].offset - a[j].offset;
                distance_gap = distance_pos > distance_self_pos? distance_pos - distance_self_pos : distance_self_pos - distance_pos;

                indels += distance_gap; 
                selfLen += distance_self_pos;
            }
        }

        if(indels != dp->indels[i])
        {
            fprintf(stderr, "indels: %lld, dp->indels[i]: %ld\n", indels, (long)dp->indels[i]);
        }

        if(selfLen != dp->self_length[i])
        {
            fprintf(stderr, "selfLen: %lld, dp->self_length[i]: %ld\n", selfLen, (long)dp->self_length[i]);
        }

    }
}

void print_chain(k_mer_hit* a, long long a_n, Chain_Data* dp, long long topN)
{
    fprintf(stderr, "topN: %lld\n", topN);
    long long max_score = -1, max_i = -1, max_n = 0;;
    long long ss, i, j, current_j;
    kvec_t(long long) si; kv_init(si);
    for (ss = 0; ss < topN && ss < a_n; ss++){
        for (i = 0, max_i = -1, max_score = -1; i < a_n; ++i) {
            for (j = 0; j < (long long)si.n; j++){
                if(i == si.a[j]) break;
            }
            if(j < (long long)si.n) continue;
            if(dp->score[i] > max_score) max_score = dp->score[i], max_i = i;
        }
        if(max_i < 0) continue;
        j = max_i; max_n = 0;
        while (j >= 0)
        {
            current_j = j;
            if(current_j == -1) continue;
            j = dp->pre[j];
            max_n++;
        }

        fprintf(stderr, "\nmax_i: %lld, max_score: %lld, max_n: %lld\n", max_i, max_score, max_n);
        
        j = max_i; 
        while (j >= 0)
        {
            current_j = j;
            if(current_j == -1) continue;

            kv_push(long long, si, current_j);
            j = dp->pre[j];
            fprintf(stderr, "self_offset: %u, offset: %u, cnt: %u, score: %d\n", 
            a[current_j].self_offset, a[current_j].offset, a[current_j].cnt, dp->score[current_j]);
        }
    }
    kv_destroy(si);
}


long long get_chainLen(long long x_beg, long long x_end, long long xLen, 
                       long long y_beg, long long y_end, long long yLen)
{
    if(x_beg <= y_beg)
    {
        y_beg = y_beg - x_beg;
        x_beg = 0;
    }
    else
    {
        x_beg = x_beg - y_beg;
        y_beg = 0;
    }

    long long x_right_length = xLen - x_end - 1;
    long long y_right_length = yLen - y_end - 1;


    if(x_right_length <= y_right_length)
    {
        x_end = xLen - 1;
        y_end = y_end + x_right_length;        
    }
    else
    {
        x_end = x_end + y_right_length;
        y_end = yLen - 1;
    }

    return x_end - x_beg + 1;
}


void debug_chain_single_site(k_mer_hit* a, long long a_n, Chain_Data* dp, int x_readLen, int y_readLen, int s_index)
{
    long long j, current_j = s_index;
    long long selfLen = 0, indels = 0;
    long long distance_self_pos, distance_pos, distance_gap;

    j = s_index;
    while (j >= 0)
    {
        current_j = j;

        j = dp->pre[j];

        if(j != -1)
        {
            distance_self_pos = a[current_j].self_offset - a[j].self_offset;
            distance_pos = a[current_j].offset - a[j].offset;
            distance_gap = distance_pos > distance_self_pos? distance_pos - distance_self_pos : distance_self_pos - distance_pos;

            indels += distance_gap; 
            selfLen += distance_self_pos;
        }
        fprintf(stderr, "j: %lld, score: %lld, occ: %d, pre_j: %lld\n", 
                    current_j, (long long)dp->score[current_j], dp->occ[current_j], j);
    }

    fprintf(stderr, "s_self_offset: %u, s_offset: %u, e_self_offset: %u, e_offset: %u, ovlp length: %lld, x_readLen: %d, y_readLen: %d\n", 
    a[s_index].self_offset, a[s_index].offset, a[current_j].self_offset, a[current_j].offset,
    get_chainLen(a[s_index].self_offset, a[current_j].self_offset, x_readLen, 
                                     a[s_index].offset, a[current_j].offset, y_readLen), x_readLen, y_readLen);

    if(indels != dp->indels[s_index])
    {
        fprintf(stderr, "indels: %lld, dp->indels[i]: %ld\n", indels, (long)dp->indels[s_index]);
    }

    if(selfLen != dp->self_length[s_index])
    {
        fprintf(stderr, "selfLen: %lld, dp->self_length[i]: %ld\n", selfLen, (long)dp->self_length[s_index]);
    }
    fprintf(stderr,"\n");
}

int32_t ha_chain_check(k_mer_hit *a, int32_t n_a, Chain_Data *dp, int32_t min_sc, double bw_thres)
{
	int32_t i, tot_indel = 0, tot_len = 0;
	double bw_pen;
	if (n_a == 0) return -1;
	for (i = 1; i < n_a; ++i)
		if (a[i-1].self_offset >= a[i].self_offset)
			break;
	if (i < n_a) return -1;
	bw_pen = 1.0 / bw_thres;
	// dp->score[0] = a[0].good? min_sc : min_sc>>1;
    dp->score[0] = normal_w(min_sc, (int64_t)a[0].cnt);
	dp->pre[0] = -1, dp->indels[0] = 0, dp->self_length[0] = 0, dp->occ[0] = 1;
	for (i = 1; i < n_a; ++i) {
		int32_t score, dg;
		int32_t dx = (int32_t)a[i].offset - (int32_t)a[i-1].offset;
		int32_t dy = (int32_t)a[i].self_offset - (int32_t)a[i-1].self_offset;
		int32_t dd = dx > dy? dx - dy : dy - dx;
		double gap_rate;
		tot_indel += dd;
		tot_len += dy;
		if (tot_indel > tot_len * bw_thres) break;
		dg = dx < dy? dx : dy;
		if (dd > THRESHOLD_MAX_SIZE && dd > dg * bw_thres) break;
		score = dg < min_sc? dg : min_sc;
        score = normal_w(score, (int64_t)a[i].cnt);
		gap_rate = (double)tot_indel / tot_len;
		score -= (int)(gap_rate * score * bw_pen);
		dp->score[i] = dp->score[i-1] + score;
		dp->pre[i] = i - 1;
		dp->indels[i] = tot_indel;
		dp->self_length[i] = tot_len;
        dp->occ[i] = i + 1;
	}
	if (i < n_a) return -1;
	return n_a;
}

///double band_width_threshold = 0.05;
long long chain_DP(k_mer_hit* a, long long a_n, Chain_Data* dp, overlap_region* result, 
              double band_width_threshold, int max_skip, int x_readLen, int y_readLen, void *km)
{
    long long i, j;
    long long self_pos, pos, max_j, max_i, max_score, score;
    long long distance_pos, distance_self_pos, distance_gap,  distance_min;
    ///double band_width_threshold = 0.05;
    double band_width_penalty = 1 / band_width_threshold;
    long long min_score = asm_opt.k_mer_length;
    long long max_indels, max_self_length;
    double gap_rate;
    long long total_indels, total_self_length;
	int32_t ret;
    
    resize_Chain_Data(dp, a_n, km);

	ret = ha_chain_check(a, a_n, dp, min_score, band_width_threshold);
	if (ret > 0) {
		a_n = ret;
		goto skip_dp;
	}

    // fill the score and backtrack arrays
	for (i = 0; i < a_n; ++i) dp->tmp[i] = -1;
	for (i = 0; i < a_n; ++i) 
    {
		int n_chn_skip = 0;
		int n_max_skip = 0;

        pos = a[i].offset;
        self_pos = a[i].self_offset;
        max_j = -1;
        max_score = normal_w(min_score, (int64_t)a[i].cnt);
        max_indels = 0;
        max_self_length = 0;

        ///may have a pre-cut condition for j
        for (j = i - 1; j >= 0; --j) 
        {
            distance_pos = pos - a[j].offset;
            distance_self_pos = self_pos - a[j].self_offset;
            ///a has been sorted by a[].offset
            ///note for a, we do not have any two elements that have both equal offsets and self_offsets
            ///but there maybe two elements that have equal offsets or equal self_offsets
            if(distance_pos == 0 || distance_self_pos <= 0)
            {
                continue;
            }

            distance_gap = distance_pos > distance_self_pos? distance_pos - distance_self_pos : distance_self_pos - distance_pos;
            
            total_indels = dp->indels[j] + distance_gap;
            total_self_length = dp->self_length[j] + distance_self_pos;
            if(total_indels > band_width_threshold * total_self_length)
            {
                continue;
            }

            ///min distance
            distance_min = distance_pos < distance_self_pos? distance_pos:distance_self_pos;
            score = distance_min < min_score? distance_min : min_score;
            ///need to be fixed in r305
			///if (!a[j].good) score = (score >> 1) + (score & 1);
            score = normal_w(score, (int64_t)a[j].cnt);

            gap_rate = (double)((double)(total_indels)/(double)(total_self_length));
            ///if the gap rate > 0.06, score will be negative
            score -= (long long)(gap_rate * score * band_width_penalty);

            score += dp->score[j];

			///find a new max score
			if (score > max_score) {///must use > instead of >=
				max_score = score;
				max_j = j;
				max_indels = total_indels;
				max_self_length = total_self_length;
				n_max_skip = 0;
				if (n_chn_skip > 0) --n_chn_skip;
			} else {
				if (++n_max_skip > max_skip)
					break;
				if (dp->tmp[j] == i) {
					if (++n_chn_skip > max_skip)
						break;
				}
			}
			if (dp->pre[j] >= 0) dp->tmp[dp->pre[j]] = i;
        }

        dp->score[i] = max_score;
        dp->pre[i] = max_j;
        dp->indels[i] = max_indels;
        dp->self_length[i] = max_self_length;
        dp->occ[i] = 1;
        if(max_j != -1) dp->occ[i] = dp->occ[max_j] + 1;
    }

    ///debug_chain(a, a_n, dp);
    // if((*result).x_id == 2162668 && (*result).y_id == 182804) print_chain(a, a_n, dp, 10);
    
skip_dp:

    max_score = -1;
    max_i = -1;
    long long mini_xLen = x_readLen * 2 + 2, tmp_xLen;
    for (i = 0; i < a_n; ++i) 
    {
        if(dp->score[i] > max_score)
        {
            max_score = dp->score[i];
            max_i = i;
            mini_xLen = get_chainLen(a[i].self_offset, a[i].self_offset, x_readLen, 
                                     a[i].offset, a[i].offset, y_readLen);
        }
        else if(dp->score[i] == max_score)
        {
            tmp_xLen = get_chainLen(a[i].self_offset, a[i].self_offset, x_readLen, 
                                    a[i].offset, a[i].offset, y_readLen);

            if(tmp_xLen < mini_xLen)
            {
                max_score = dp->score[i];
                max_i = i;

                mini_xLen = tmp_xLen;
            }
        }
        
    }

    clear_fake_cigar(&(result->f_cigar));
    ///note a has been sorted by offset, that means has been sorted by query offset
    i = max_i;
    result->x_pos_e = a[i].self_offset;
    result->y_pos_e = a[i].offset;
    result->shared_seed = max_score;
    result->overlapLen = mini_xLen;

    distance_self_pos = result->x_pos_e - a[i].self_offset;
    distance_pos = result->y_pos_e - a[i].offset;
    long long pre_distance_gap = distance_pos - distance_self_pos;
    ///record first site
    ///the length of f_cigar should be at least 1
    ///record the offset of reference
    add_fake_cigar(&(result->f_cigar), a[i].self_offset, pre_distance_gap, km);
    long long chainLen = 0;
    if(result->x_pos_strand == 1)
    {
        while (i >= 0)
        {
            distance_self_pos = result->x_pos_e - a[i].self_offset;
            distance_pos = result->y_pos_e - a[i].offset;
            distance_gap = distance_pos - distance_self_pos;
            if(distance_gap != pre_distance_gap)
            {
                pre_distance_gap = distance_gap;
                ///record this site
                add_fake_cigar(&(result->f_cigar), a[i].self_offset, pre_distance_gap, km);
            }

            chainLen++;
            result->x_pos_s = a[i].self_offset;
            result->y_pos_s = a[i].offset;
            i = dp->pre[i];
        }
    }
    else
    {
        while (i >= 0)
        {
            distance_self_pos = result->x_pos_e - a[i].self_offset;
            distance_pos = result->y_pos_e - a[i].offset;
            distance_gap = distance_pos - distance_self_pos;
            if(distance_gap == pre_distance_gap)
            {
                result->f_cigar.length--;
                add_fake_cigar(&(result->f_cigar), a[i].self_offset, pre_distance_gap, km);
            }
            else
            {
                pre_distance_gap = distance_gap;
                add_fake_cigar(&(result->f_cigar), a[i].self_offset, pre_distance_gap, km);
            }

            chainLen++;
            result->x_pos_s = a[i].self_offset;
            result->y_pos_s = a[i].offset;
            i = dp->pre[i];
        }
    }
    return chainLen;
}



void calculate_overlap_region_by_chaining(Candidates_list* candidates, overlap_region_alloc* overlap_list, kvec_t_u64_warp* chain_idx,
                                          uint64_t readID, uint64_t readLength, All_reads* R_INF, const ul_idx_t *uref, double band_width_threshold, int add_beg_end, overlap_region* f_cigar, void *km)
{
    long long i = 0;
    uint64_t current_ID;
    uint64_t current_stand;

    if (candidates->length == 0) return;

    long long sub_region_beg;
    long long sub_region_end;

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

		chain_DP(candidates->list + sub_region_beg,
				sub_region_end - sub_region_beg + 1, &(candidates->chainDP), f_cigar, band_width_threshold,
				25, /**Get_READ_LENGTH((*R_INF), (*f_cigar).x_id)**/readLength, 
                R_INF?Get_READ_LENGTH((*R_INF), (*f_cigar).y_id):uref->ug->u.a[(*f_cigar).y_id].len, km);

        ///if (tmp_region.x_id != tmp_region.y_id && tmp_region.shared_seed > 1)
        if ((*f_cigar).x_id != (*f_cigar).y_id)
        {
            append_inexact_overlap_region_alloc(overlap_list, f_cigar, readLength, R_INF?Get_READ_LENGTH((*R_INF), (*f_cigar).y_id):uref->ug->u.a[(*f_cigar).y_id].len, add_beg_end, km);
        }
    }
}


void append_window_list(overlap_region* region, uint64_t x_start, uint64_t x_end, int y_start, int y_end, int error,
                        int extra_begin, int extra_end, int error_threshold, int blockLen, void *km)
{
    window_list *p = NULL;
    kv_pushp(window_list, region->w_list, &p);

    p->x_start = x_start;
    p->x_end = x_end;
    p->y_start = y_start;
    p->y_end = y_end;
    p->error = error;
    p->extra_begin = extra_begin;
    p->extra_end = extra_end;
    p->error_threshold = error_threshold;
    p->cidx = p->clen = 0;
}

void test_single_list(Candidates_list* candidates, k_mer_pos* n_list, uint64_t n_lengh, uint64_t end_pos, uint64_t strand)
{
    uint64_t i;
    long long j = 0;
    for (i = 0; i < n_lengh; i++)
    {
        
        for (; j < candidates->length; j++)
        {
            if (
                n_list[i].offset == (uint64_t)candidates->list[j].offset
                &&
                n_list[i].readID == candidates->list[j].readID
                &&
                end_pos == (uint64_t)candidates->list[j].self_offset
                &&
                strand == candidates->list[j].strand
            )
            {
                break;
            }   
        }

        if (j == candidates->length)
        {
            fprintf(stderr, "ERROR 4\n");
        }
    }
}


void init_Chain_Data(Chain_Data* x)
{
	memset(x, 0, sizeof(Chain_Data));
}

void clear_Chain_Data(Chain_Data* x)
{
    x->length = 0;
}

void destory_Chain_Data(Chain_Data* x)
{
    free(x->score);
    free(x->pre);
    free(x->indels);
    free(x->self_length);
    free(x->occ);
	free(x->tmp);
}

void resize_Chain_Data(Chain_Data* x, long long size, void *km)
{
	if (size + 1 > x->size) {
		x->size = size + 1;
		kroundup64(x->size);
        if(!km) {
            REALLOC(x->score, x->size);
            REALLOC(x->pre, x->size);
            REALLOC(x->indels, x->size);
            REALLOC(x->self_length, x->size);
            REALLOC(x->occ, x->size);
            REALLOC(x->tmp, x->size);
        } else {
            KREALLOC(km, x->score, x->size);
            KREALLOC(km, x->pre, x->size);
            KREALLOC(km, x->indels, x->size);
            KREALLOC(km, x->self_length, x->size);
            KREALLOC(km, x->occ, x->size);
            KREALLOC(km, x->tmp, x->size);
        }
	}
}

void init_Candidates_list(Candidates_list* l)
{
    l->length = 0;
    l->size = 0;
    l->list = NULL;
    init_Chain_Data(&(l->chainDP));
}

void clear_Candidates_list(Candidates_list* l)
{
    l->length = 0;
    clear_Chain_Data(&(l->chainDP));
}

void destory_Candidates_list(Candidates_list* l)
{
    free(l->list);
    destory_Chain_Data(&(l->chainDP));
}

void destory_Candidates_list_buf(void *km, Candidates_list* l, int is_z)
{
    kfree(km, l->list);
    kfree(km, l->chainDP.score);
    kfree(km, l->chainDP.pre);
    kfree(km, l->chainDP.indels);
    kfree(km, l->chainDP.self_length);
    kfree(km, l->chainDP.occ);
    kfree(km, l->chainDP.tmp);
    if(is_z) memset(l, 0, sizeof(*l));
}

void init_fake_cigar(Fake_Cigar* x)
{
    x->buffer = NULL;
    x->length = 0;
    x->size = 0;
}

void destory_fake_cigar(Fake_Cigar* x)
{
    if(x->size > 0)
    {
        free(x->buffer);
    }
}

void clear_fake_cigar(Fake_Cigar* x)
{
    x->length = 0;
}

void add_fake_cigar(Fake_Cigar* x, uint32_t gap_site, int32_t gap_shift, void *km)
{
    if(x->length + 1 > x->size)
    {
        x->size = x->length + 1;
        kroundup32(x->size);
        if(!km) {
            REALLOC(x->buffer, x->size);
        } else {
            KREALLOC(km, x->buffer, x->size);
        }
        // x->buffer = (uint64_t*)realloc(x->buffer, sizeof(uint64_t) * x->size);
    }

    x->buffer[x->length] = gap_site;
    x->buffer[x->length] = x->buffer[x->length] << 32;

    if(gap_shift < 0)
    {
        gap_shift = gap_shift * -1;
        gap_site = gap_shift;
        gap_site = gap_site << 1;
        gap_site = gap_site | ((uint32_t)1);
    }
    else
    {
        gap_site = gap_shift;
        gap_site = gap_site << 1;
    }

    x->buffer[x->length] = x->buffer[x->length] | ((uint32_t)gap_site);

    x->length++;
}


void resize_fake_cigar(Fake_Cigar* x, uint64_t size, void *km)
{
    if (size > x->size) {
        x->size = size;
        if(!km) {
            REALLOC(x->buffer, x->size);
        }
        else {
            KREALLOC(km, x->buffer, x->size);
        }
		
    }
    x->length = 0;
}


void init_window_list_alloc(window_list_alloc* x)
{
    memset(x, 0, sizeof((*x)));
}

void clear_window_list_alloc(window_list_alloc* x)
{
    x->n = x->c.n = 0;
}

void destory_window_list_alloc(window_list_alloc* x)
{
    free(x->a); free(x->c.a);
}

void resize_window_list_alloc(window_list_alloc* x, uint64_t size)
{
    kv_resize(window_list, *x, size); x->n = x->c.n = 0;
    uint64_t k;
    for (k = 0; k < x->m; k++) x->a[k].error = -1; 
    x->c.n = 0;
}

#define normal_cw(x) (((x)&(0xffu))>=((x)>>8)?(((x)&(0xffu))/((x)>>8)):1)

int32_t lchain_check(k_mer_hit *a, int32_t n_a, Chain_Data *dp, double bw_thres)
{
	int32_t i, tot_g = 0, sc, dg, dq, dr, dd, span;
	double bw_pen;
	if (n_a == 0) return -1;
    if (n_a > 1) {
        if ((a[0].self_offset >= a[n_a-1].self_offset)||(a[0].offset == a[n_a-1].offset)) return -1;
        dq = (int32_t)a[n_a-1].self_offset - (int32_t)a[0].self_offset;
        dr = (int32_t)a[n_a-1].offset - (int32_t)a[0].offset;
        dd = ((dq>=dr)? (dq-dr): (dr-dq));//gap
        dg = ((dq>=dr)? (dr): (dq));///len
        if (dg == 0 || dd > (dg*bw_thres)) return -1;
    }

	for (i = 1; i < n_a; ++i) {///a[] is sorted by offset, instead of self_offset; but offset might be equal
        if(a[i-1].self_offset >= a[i].self_offset) break;
        if(a[i-1].offset == a[i].offset) break;
    }
	if (i < n_a) return -1;

	bw_pen = 1.0 / bw_thres;
    dp->score[0] = normal_w((a[0].cnt&(0xffu)), (a[0].cnt>>8)); dp->pre[0] = -1;
	for (i = 1; i < n_a; ++i) {
        dq = (int32_t)a[i].self_offset - (int32_t)a[i-1].self_offset;
		dr = (int32_t)a[i].offset - (int32_t)a[i-1].offset;
		dd = ((dq>=dr)? (dq-dr): (dr-dq));//gap
        dg = ((dq>=dr)? (dr): (dq));///len
        if(dg == 0) break;
		
		tot_g += dd;
		if (dd > THRESHOLD_MAX_SIZE && dd > (dg*bw_thres)) break;
        span = a[i].cnt&(0xffu);
		sc = dg < span? dg : span;
        sc = normal_w(sc, ((int32_t)(a[i].cnt>>8)));
		sc -= (int32_t)((((double)dd)/((double)dg))*bw_pen*((double)sc));///bw_pen is 20 for HiFi

		dp->score[i] = dp->score[i-1] + sc;
		dp->pre[i] = i - 1;
	}
	if (i < n_a) return -1;

    if(n_a > 1) {
        dq = (int32_t)a[n_a-1].self_offset - (int32_t)a[0].self_offset;
        dr = (int32_t)a[n_a-1].offset - (int32_t)a[0].offset;
        dg = ((dq>=dr)? (dr): (dq));///len
        dd = tot_g;///gap
        if (dd > (dg*bw_thres)) return -1;
    }
	return n_a;
}


int32_t lchain_qcheck(k_mer_hit *a, int32_t n_a, Chain_Data *dp, double bw_thres)
{
	int32_t i, tot_g = 0, sc, dg, dq, dr, dd, span;
	double bw_pen;
	if (n_a == 0) return -1;
    if (n_a > 1) {
        if ((a[0].self_offset >= a[n_a-1].self_offset)||(a[0].offset >= a[n_a-1].offset)) return -1;
        dq = (int32_t)a[n_a-1].self_offset - (int32_t)a[0].self_offset;
        dr = (int32_t)a[n_a-1].offset - (int32_t)a[0].offset;
        dd = ((dq>=dr)? (dq-dr): (dr-dq));//gap
        dg = ((dq>=dr)? (dr): (dq));///len
        if (dg == 0 || dd > (dg*bw_thres)) return -1;
    }

	for (i = 1; i < n_a; ++i) {///a[] is sorted by self_offset
        if(a[i-1].self_offset >= a[i].self_offset) break;
        if(a[i-1].offset >= a[i].offset) break;
    }
	if (i < n_a) return -1;

	bw_pen = 1.0 / bw_thres;
    dp->score[0] = normal_w((a[0].cnt&(0xffu)), (a[0].cnt>>8)); dp->pre[0] = -1;
	for (i = 1; i < n_a; ++i) {
        dq = (int32_t)a[i].self_offset - (int32_t)a[i-1].self_offset;
		dr = (int32_t)a[i].offset - (int32_t)a[i-1].offset;
		dd = ((dq>=dr)? (dq-dr): (dr-dq));//gap
        dg = ((dq>=dr)? (dr): (dq));///len
        if(dg == 0) break;
		
		tot_g += dd;
		if (dd > THRESHOLD_MAX_SIZE && dd > (dg*bw_thres)) break;
        span = a[i].cnt&(0xffu);
		sc = dg < span? dg : span;
        sc = normal_w(sc, ((int32_t)(a[i].cnt>>8)));
		sc -= (int32_t)((((double)dd)/((double)dg))*bw_pen*((double)sc));///bw_pen is 20 for HiFi

		dp->score[i] = dp->score[i-1] + sc;
		dp->pre[i] = i - 1;
	}
	if (i < n_a) return -1;

    if(n_a > 1) {
        dq = (int32_t)a[n_a-1].self_offset - (int32_t)a[0].self_offset;
        dr = (int32_t)a[n_a-1].offset - (int32_t)a[0].offset;
        dg = ((dq>=dr)? (dr): (dq));///len
        dd = tot_g;///gap
        if (dd > (dg*bw_thres)) return -1;
    }
	return n_a;
}

inline int32_t cal_bw(const k_mer_hit *ai, const k_mer_hit *aj, double bw_rate, int64_t sf_l, int64_t ot_l)
{
    ///ai is the suffix of aj
    int64_t sf_s = aj->self_offset, sf_e = ai->self_offset + 1;
    int64_t ot_s = aj->offset, ot_e = ai->offset + 1;
    int64_t sf_r = sf_l - sf_e, ot_r = ot_l - ot_e;
    if(sf_s <= ot_s) sf_s = 0;
    else sf_s -= ot_s;

    if(sf_r <= ot_r) sf_e = sf_l;
    else sf_e += ot_r;

    return (sf_e - sf_s)*bw_rate;
}

inline int32_t comput_sc_ch(const k_mer_hit *ai, const k_mer_hit *aj, double bw_rate, double chn_pen_gap, double chn_pen_skip, int64_t sl, int64_t ol)
{
    ///ai is the suffix of aj
    int32_t dq, dr, dd, dg, q_span, sc; 
    dq = (int64_t)(ai->self_offset) - (int64_t)(aj->self_offset);
    if(dq <= 0) return INT32_MIN;
    dr = (int64_t)(ai->offset) - (int64_t)(aj->offset);
    if(dr <= 0) return INT32_MIN;
    dd = dr > dq? dr - dq : dq - dr;//gap
    if((dd > 16) && (dd > cal_bw(ai, aj, bw_rate, sl, ol))) return INT32_MIN;
    dg = dr < dq? dr : dq;//len
    q_span = ai->cnt&(0xffu); 
    sc = q_span < dg? q_span : dg;
    sc = normal_w(sc, ((int32_t)(ai->cnt>>8)));
    if (dd || (dg > q_span && dg > 0)) {
        double lin_pen, a_pen;
        lin_pen = (chn_pen_gap*(double)dd);
        a_pen = ((double)(sc))*((((double)dd)/((double)dg))/bw_rate);
        if(lin_pen > a_pen) lin_pen = a_pen;
        lin_pen += (chn_pen_skip*(double)dg);
        sc -= (int32_t)lin_pen;
    }
    return sc;
}

inline int32_t comput_sc_ff(const k_mer_hit *ai, const k_mer_hit *aj, double bw_rate, double chn_pen_gap, double chn_pen_skip, int64_t sl, int64_t ol)
{
    ///ai is the suffix of aj
    int32_t dq, dr, dd, dg, q_span, sc; 
    dq = (int64_t)(ai->self_offset) - (int64_t)(aj->self_offset);
    if(dq < 0) return INT32_MIN;
    dr = (int64_t)(ai->offset) - (int64_t)(aj->offset);
    if(dr < 0) return INT32_MIN;
    dd = dr > dq? dr - dq : dq - dr;//gap
    // if((dd > 16) && (dd > cal_bw(ai, aj, bw_rate, sl, ol))) return INT32_MIN;
    dg = dr < dq? dr : dq;//len
    q_span = ai->cnt&(0xffu); 
    sc = q_span < dg? q_span : dg;
    sc = normal_w(sc, ((int32_t)(ai->cnt>>8)));
    if (dd || (dg > q_span && dg > 0)) {
        double lin_pen, a_pen;
        lin_pen = (chn_pen_gap*(double)dd);
        a_pen = ((double)(sc))*((((double)dd)/((double)dg))/bw_rate);
        if(lin_pen > a_pen) lin_pen = a_pen;
        lin_pen += (chn_pen_skip*(double)dg);
        sc -= (int32_t)lin_pen;
    }
    return sc;
}

uint64_t lchain_dp(k_mer_hit* a, int64_t a_n, k_mer_hit* des, Chain_Data* dp, overlap_region* res, 
              int64_t max_skip, int64_t max_iter, int64_t max_dis, double chn_pen_gap, double chn_pen_skip, double bw_rate, 
              int64_t xl, int64_t yl, int64_t quick_check)
{
    int64_t *p, *t, max_f, n_skip, st, max_j, end_j, sc, msc, msc_i, bw, max_ii, ovl, movl; 
    int32_t *f, max, tmp; int64_t i, j, ret, cL = 0;
    resize_Chain_Data(dp, a_n, NULL);
    t = dp->tmp; f = dp->score; p = dp->pre;
    bw = ((xl < yl)?xl:yl); bw *= bw_rate;
    msc = msc_i = -1; movl = INT32_MAX;

    if(quick_check) {
        ret = lchain_check(a, a_n, dp, bw_rate);
        if (ret > 0) {
            a_n = ret; msc_i = a_n-1; msc = f[msc_i];
            goto skip_ldp;
        }
    }

    memset(t, 0, (a_n*sizeof((*t))));
	for (i = st = 0, max_ii = -1; i < a_n; ++i) {
        max_f = a[i].cnt&(0xffu); 
        n_skip = 0; max_j = end_j = -1;
        if ((i-st) > max_iter) st = i-max_iter;

        for (j = i - 1; j >= st; --j) {
			sc = comput_sc_ch(&a[i], &a[j], bw_rate, chn_pen_gap, chn_pen_skip, xl, yl);
			if (sc == INT32_MIN) continue;
			sc += f[j];
			if (sc > max_f) {
				max_f = sc, max_j = j;
				if (n_skip > 0) --n_skip;
			} else if (t[j] == (int32_t)i) {
				if (++n_skip > max_skip)
					break;
			}
			if (p[j] >= 0) t[p[j]] = i;
		}
        end_j = j;

        if (max_ii < 0 || ((int64_t)a[i].offset) - ((int64_t)a[max_ii].offset) > max_dis) {
			max = INT32_MIN; max_ii = -1;
            for (j = i - 1; (j >= st) && ((((int64_t)a[i].offset)-((int64_t)a[j].offset))<=max_dis); --j) {
                if (max < f[j]) {
                    max = f[j], max_ii = j;
                }
            }
		}

        if (max_ii >= 0 && max_ii < end_j) {///just have a try with a[i]<->a[max_ii]
			tmp = comput_sc_ch(&a[i], &a[max_ii], bw_rate, chn_pen_gap, chn_pen_skip, xl, yl);
			if (tmp != INT32_MIN && max_f < tmp + f[max_ii])
				max_f = tmp + f[max_ii], max_j = max_ii;
		}
        f[i] = max_f; p[i] = max_j;
        if ((max_ii < 0) || (((((int64_t)a[i].offset)-((int64_t)a[max_ii].offset))<=max_dis) && (f[max_ii]<f[i]))) {
            max_ii = i;
        }
        if(f[i] >= msc) {
            ovl = get_chainLen(a[i].self_offset, a[i].self_offset, xl, a[i].offset, a[i].offset, yl);
            if(f[i] > msc || ovl < movl) {
                msc = f[i]; msc_i = i; movl = ovl;
            }
        }
    }
    
    skip_ldp:
    ///a[] has been sorted by offset
    i = msc_i; 
    res->x_pos_s = res->x_pos_e = a[i].self_offset;
    res->y_pos_s = res->y_pos_e = a[i].offset;
    res->shared_seed = msc;

    cL = 0; 
    while (i >= 0) {
        t[cL++] = i; msc_i = i; i = p[i];
    }

    res->x_pos_s = a[t[cL-1]].self_offset;
    res->y_pos_s = a[t[cL-1]].offset;
    res->overlapLen = get_chainLen(res->x_pos_s, res->x_pos_e, xl, res->y_pos_s, res->y_pos_e, yl);
    for (i = 0; i < cL; i++) des[i] = a[t[cL-i-1]];
    return cL;
}

uint64_t lchain_qdp(k_mer_hit* a, int64_t a_n, k_mer_hit* des, Chain_Data* dp, overlap_region* res, 
              int64_t max_skip, int64_t max_iter, int64_t max_dis, double chn_pen_gap, double chn_pen_skip, double bw_rate, 
              int64_t xl, int64_t yl, int64_t quick_check)
{
    int64_t *p, *t, max_f, n_skip, st, max_j, end_j, sc, msc, msc_i, bw, max_ii, ovl, movl; 
    int32_t *f, max, tmp; int64_t i, j, ret, cL = 0;
    resize_Chain_Data(dp, a_n, NULL);
    t = dp->tmp; f = dp->score; p = dp->pre;
    bw = ((xl < yl)?xl:yl); bw *= bw_rate;
    msc = msc_i = -1; movl = INT32_MAX;
    // if(a_n && a[0].readID == 0) {
    //     fprintf(stderr, "---[M::%s::utg%.6dl::%c]\n", 
    //                     __func__, (int32_t)a[0].readID+1, "+-"[a[0].strand]);
    // }
    if(quick_check) {
        ret = lchain_qcheck(a, a_n, dp, bw_rate);
        if (ret > 0) {
            a_n = ret; msc_i = a_n-1; msc = f[msc_i];
            goto skip_ldp;
        }
    }

    memset(t, 0, (a_n*sizeof((*t))));
    for (i = st = 0, max_ii = -1; i < a_n; ++i) {
        max_f = a[i].cnt&(0xffu); 
        n_skip = 0; max_j = end_j = -1;
        if ((i-st) > max_iter) st = i-max_iter;

        for (j = i - 1; j >= st; --j) {
            sc = comput_sc_ch(&a[i], &a[j], bw_rate, chn_pen_gap, chn_pen_skip, xl, yl);
            if (sc == INT32_MIN) continue;
            sc += f[j];
            if (sc > max_f) {
                max_f = sc, max_j = j;
                if (n_skip > 0) --n_skip;
            } else if (t[j] == (int32_t)i) {
                if (++n_skip > max_skip)
                    break;
            }
            if (p[j] >= 0) t[p[j]] = i;
        }
        end_j = j;

        if (max_ii < 0 || ((int64_t)a[i].self_offset) - ((int64_t)a[max_ii].self_offset) > max_dis) {
            max = INT32_MIN; max_ii = -1;
            for (j = i - 1; (j >= st) && ((((int64_t)a[i].self_offset)-((int64_t)a[j].self_offset))<=max_dis); --j) {
                if (max < f[j]) {
                    max = f[j], max_ii = j;
                }
            }
        }

        if (max_ii >= 0 && max_ii < end_j) {///just have a try with a[i]<->a[max_ii]
            tmp = comput_sc_ch(&a[i], &a[max_ii], bw_rate, chn_pen_gap, chn_pen_skip, xl, yl);
            if (tmp != INT32_MIN && max_f < tmp + f[max_ii])
                max_f = tmp + f[max_ii], max_j = max_ii;
        }
        f[i] = max_f; p[i] = max_j;
        if ((max_ii < 0) || (((((int64_t)a[i].self_offset)-((int64_t)a[max_ii].self_offset))<=max_dis) && (f[max_ii]<f[i]))) {
            max_ii = i;
        }
        if(f[i] >= msc) {
            ovl = get_chainLen(a[i].self_offset, a[i].self_offset, xl, a[i].offset, a[i].offset, yl);
            if(f[i] > msc || ovl < movl) {
                msc = f[i]; msc_i = i; movl = ovl;
            }
        }
        // if(a_n && a[0].readID == 0) {
        //     fprintf(stderr, "i::%ld[M::%s::utg%.6dl::%c] x::%u, y::%u, st::%ld, max_ii::%ld, f[i]::%d, p[i]::%ld, msc_i::%ld, msc::%ld, movl::%ld\n", 
        //                     i, __func__, (int32_t)a[i].readID+1, "+-"[a[i].strand], 
        //                     a[i].self_offset, a[i].offset, st, max_ii, f[i], p[i], msc_i, msc, movl);
        // }
    }
    
    skip_ldp:
    ///a[] has been sorted by self_offset
    i = msc_i; 
    res->x_pos_s = res->x_pos_e = a[i].self_offset;
    res->y_pos_s = res->y_pos_e = a[i].offset;
    res->shared_seed = msc;

    cL = 0; 
    while (i >= 0) {
        t[cL++] = i; msc_i = i; i = p[i];
    }

    res->x_pos_s = a[t[cL-1]].self_offset;
    res->y_pos_s = a[t[cL-1]].offset;
    res->overlapLen = get_chainLen(res->x_pos_s, res->x_pos_e, xl, res->y_pos_s, res->y_pos_e, yl);
    for (i = 0; i < cL; i++) {
        des[i] = a[t[cL-i-1]];
        // if(a_n && a[0].readID == 0) {
        //     fprintf(stderr, "i::%ld[M::%s::utg%.6dl::%c] x::%u, y::%u, cL::%ld\n", 
        //                     i, __func__, (int32_t)des[i].readID+1, "+-"[des[i].strand], des[i].self_offset, des[i].offset, cL);
        // }
    }
    return cL;
}


#define kv_pushp_ol(type, v, p) do {									\
		if ((v).length == (v).size) {										\
			(v).list = (type*)realloc((v).list, sizeof(type)*((v).size?((v).size<<1):(2)));	\
            memset((v).list+(v).size, 0, sizeof(overlap_region)*(((v).size?((v).size<<1):2)-(v).size));\
            (v).size = (v).size?((v).size<<1):(2);							\
		}															\
		*(p) = &((v).list[(v).length++]); \
	} while (0)

void push_ovlp_chain_qgen(overlap_region* o, uint32_t xid, int64_t xl, int64_t yl, int64_t sc, 
k_mer_hit *beg, k_mer_hit *end)
{
    int64_t xr, yr; 
    o->x_id = xid; o->y_id = beg->readID;
    o->x_pos_strand = 0; o->y_pos_strand = beg->strand;
    o->x_pos_s = beg->self_offset; o->y_pos_s = beg->offset;
    o->x_pos_e = end->self_offset; o->y_pos_e = end->offset;

    if(o->x_pos_s <= o->y_pos_s) {
        o->y_pos_s -= o->x_pos_s; o->x_pos_s = 0;
    } else {
        o->x_pos_s -= o->y_pos_s; o->y_pos_s = 0;
    }

    xr = xl-o->x_pos_e-1; yr = yl-o->y_pos_e-1;
    if(xr <= yr) {
        o->x_pos_e = xl-1; o->y_pos_e += xr;        
    } else {
        o->y_pos_e = yl-1; o->x_pos_e += yr; 
    }

    o->shared_seed = sc;
    o->align_length = 0;
    o->is_match = 0;
    o->non_homopolymer_errors = 0;
    o->strong = 0;
    o->overlapLen = 0;
}

int64_t filter_non_ovlp_chains(overlap_region *a, int64_t a_n, int64_t *n_v)
{
    int64_t k, i, n_mchain, omx, omy, opx, opy, ovx, ovy, os, oe; overlap_region *m, *p, t;
    for (k = n_mchain = (*n_v) = 0; k < a_n; k++) {
		m = &(a[k]); 
        omx = m->x_pos_e + 1 - m->x_pos_s;
        omy = m->y_pos_e + 1 - m->y_pos_s;
		for (i = 0; i < n_mchain; i++) {
			p = &(a[i]);
            opx = p->x_pos_e + 1 - p->x_pos_s;
            opy = p->y_pos_e + 1 - p->y_pos_s;

            os = ((m->x_pos_s>=p->x_pos_s)? m->x_pos_s:p->x_pos_s);
            oe = ((m->x_pos_e<=p->x_pos_e)? m->x_pos_e:p->x_pos_e) + 1;
			ovx = oe>os?oe-os:0;
            if((ovx > omx*0.1) || (ovx > opx*0.1)) break;

            os = ((m->y_pos_s>=p->y_pos_s)? m->y_pos_s:p->y_pos_s);
            oe = ((m->y_pos_e<=p->y_pos_e)? m->y_pos_e:p->y_pos_e) + 1;
			ovy = oe>os?oe-os:0;
			if((ovy > omy*0.1) || (ovy > opy*0.1)) break;
		}
		if(i < n_mchain) continue;

        if (n_mchain != k) {
            t = a[k]; a[k] = a[n_mchain]; a[n_mchain] = t;
        }
        (*n_v) += a[n_mchain].align_length;
		n_mchain++;
	}
	return n_mchain;
}

int64_t filter_non_ovlp_xchains(overlap_region *a, int64_t a_n, int64_t *n_v)
{
    int64_t k, i, n_mchain, omx, opx, ovx, os, oe; overlap_region *m, *p, t;
    for (k = n_mchain = (*n_v) = 0; k < a_n; k++) {
		m = &(a[k]); 
        omx = m->x_pos_e + 1 - m->x_pos_s;
		for (i = 0; i < n_mchain; i++) {
			p = &(a[i]);
            opx = p->x_pos_e + 1 - p->x_pos_s;

            os = ((m->x_pos_s>=p->x_pos_s)? m->x_pos_s:p->x_pos_s);
            oe = ((m->x_pos_e<=p->x_pos_e)? m->x_pos_e:p->x_pos_e) + 1;
			ovx = oe>os?oe-os:0;
            if((ovx > omx*0.1) || (ovx > opx*0.1)) break;
		}
		if(i < n_mchain) continue;

        if (n_mchain != k) {
            t = a[k]; a[k] = a[n_mchain]; a[n_mchain] = t;
        }
        (*n_v) += a[n_mchain].align_length;
		n_mchain++;
	}
	return n_mchain;
}

uint64_t lchain_qdp_mcopy(Candidates_list *cl, int64_t a_idx, int64_t a_n, int64_t des_idx, 
              Chain_Data* dp, overlap_region_alloc* res, int64_t max_skip, int64_t max_iter, 
              int64_t max_dis, double chn_pen_gap, double chn_pen_skip, double bw_rate, 
              uint32_t xid, int64_t xl, int64_t yl, int64_t quick_check, uint32_t apend_be, 
              int64_t gen_cigar, int64_t enable_mcopy, double mcopy_rate, int64_t mcopy_khit_cutoff, 
              int64_t khit_n)
{
    if(a_n <= 0) return 0;
    int64_t *p, *t, max_f, n_skip, st, max_j, end_j, sc, msc, msc_i, bw, max_ii, ovl, movl, plus = 0, min_sc, ch_n; 
    int32_t *f, max, tmp, *ii; int64_t i, k, j, cL = 0; k_mer_hit* a; k_mer_hit* des; k_mer_hit *swap; overlap_region *z;
    resize_Chain_Data(dp, a_n, NULL);
    t = dp->tmp; f = dp->score; p = dp->pre; ii = dp->occ;
    bw = ((xl < yl)?xl:yl); bw *= bw_rate;
    msc = msc_i = INT32_MIN; movl = INT32_MAX; ch_n = 1;
    a = cl->list + a_idx; des = cl->list + des_idx;
    // if(a_n && a[0].readID == 0) {
    //     fprintf(stderr, "---[M::%s::utg%.6dl::%c]\n", 
    //                     __func__, (int32_t)a[0].readID+1, "+-"[a[0].strand]);
    // }

    memset(t, 0, (a_n*sizeof((*t))));
    for (i = st = plus = 0, max_ii = -1; i < a_n; ++i) {
        max_f = a[i].cnt&(0xffu); 
        n_skip = 0; max_j = end_j = -1;
        if ((i-st) > max_iter) st = i-max_iter;
        while (a[i].strand != a[st].strand) ++st;

        for (j = i - 1; j >= st; --j) {
            sc = comput_sc_ch(&a[i], &a[j], bw_rate, chn_pen_gap, chn_pen_skip, xl, yl);
            if (sc == INT32_MIN) continue;
            sc += f[j];
            if (sc > max_f) {
                max_f = sc, max_j = j;
                if (n_skip > 0) --n_skip;
            } else if (t[j] == (int32_t)i) {
                if (++n_skip > max_skip)
                    break;
            }
            if (p[j] >= 0) t[p[j]] = i;
        }
        end_j = j;

        if ((max_ii<0) || (a[i].self_offset>a[max_ii].self_offset+max_dis) || (a[i].strand!=a[max_ii].strand)) {
            max = INT32_MIN; max_ii = -1;
            for (j=i-1; (j>=st) && (a[i].self_offset<=max_dis+a[j].self_offset)&&(a[i].strand==a[j].strand); --j) {
                if (max < f[j]) {
                    max = f[j], max_ii = j;
                }
            }
        }

        if (max_ii >= 0 && max_ii < end_j) {///just have a try with a[i]<->a[max_ii]
            tmp = comput_sc_ch(&a[i], &a[max_ii], bw_rate, chn_pen_gap, chn_pen_skip, xl, yl);
            if (tmp != INT32_MIN && max_f < tmp + f[max_ii])
                max_f = tmp + f[max_ii], max_j = max_ii;
        }
        f[i] = max_f; p[i] = max_j;
        if ((max_ii < 0) || ((a[i].self_offset<=max_dis+a[max_ii].self_offset)&&(a[i].strand==a[max_ii].strand)&&(f[max_ii]<f[i]))) {
            max_ii = i;
        }
        if(f[i] >= msc) {
            ovl = get_chainLen(a[i].self_offset, a[i].self_offset, xl, a[i].offset, a[i].offset, yl);
            if(f[i] > msc || ovl < movl) {
                msc = f[i]; msc_i = i; movl = ovl;
            }
        }
        if(f[i] < plus) plus = f[i];
        ii[i] = 0;///for mcopy, not here
        // if(a_n && a[0].readID == 0) {
        //     fprintf(stderr, "i::%ld[M::%s::utg%.6dl::%c] x::%u, y::%u, st::%ld, max_ii::%ld, f[i]::%d, p[i]::%ld, msc_i::%ld, msc::%ld, movl::%ld\n", 
        //                     i, __func__, (int32_t)a[i].readID+1, "+-"[a[i].strand], 
        //                     a[i].self_offset, a[i].offset, st, max_ii, f[i], p[i], msc_i, msc, movl);
        // }
    }

    for (i = msc_i, cL = 0; i >= 0; i = p[i]) { ii[i] = 1; t[cL++] = i;}///label the best chain
    if((movl < xl) && enable_mcopy/**(movl < yl)**/) {
        if(cL >= mcopy_khit_cutoff) {///if there are too few k-mers, disable mcopy
            msc -= plus; min_sc = msc*mcopy_rate/**0.2**/; ii[msc_i] = 0;
            for (i = ch_n = 0; i < a_n; ++i) {///make all f[] positive
                f[i] -= plus; if(i >= ch_n) t[i] = 0;
                if((!(ii[i])) && (f[i] >= min_sc)) {
                    t[ch_n] = ((uint64_t)f[i])<<32; t[ch_n] += (i<<1); ch_n++;
                }
            }
            if(ch_n > 1) {
                int64_t n_v, n_v0, ni, n_u, n_u0 = res->length; 
                radix_sort_hc64i(t, t + ch_n);
                for (k = ch_n-1, n_v = n_u = 0; k >= 0; --k) {
                    n_v0 = n_v;
                    for (i = ((uint32_t)t[k])>>1; i >= 0 && (t[i]&1) == 0; ) {
                        ii[n_v++] = i; t[i] |= 1; i = p[i];
                    }
                    if(n_v0 == n_v) continue;
                    sc = (i<0?(t[k]>>32):((t[k]>>32)-f[i]));
                    if(sc >= min_sc) {
                        kv_pushp_ol(overlap_region, (*res), &z);
                        push_ovlp_chain_qgen(z, xid, xl, yl, sc+plus, &(a[ii[n_v-1]]), &(a[ii[n_v0]]));
                        ///mcopy_khit_cutoff <= 1: disable the mcopy_khit_cutoff filtering, for the realignment
                        if((mcopy_khit_cutoff <= 1) || ((z->x_pos_e+1-z->x_pos_s) <= (movl<<2))) {
                            z->align_length = n_v-n_v0; z->x_id = n_v0;
                            n_u++;
                        } else {///non-best is too long
                            res->length--; n_v = n_v0;
                        }
                    } else {
                        n_v = n_v0;
                    }
                }

                if(n_u > 1) ks_introsort_or_sss(n_u, res->list + n_u0); 
                res->length = n_u0 + filter_non_ovlp_xchains(res->list + n_u0, n_u, &n_v);
                n_u = res->length;
                if(n_u > n_u0 + 1) {
                    kv_resize_cl(k_mer_hit, (*cl), (n_v+cl->length));
                    a = cl->list + a_idx; des = cl->list + des_idx; swap = cl->list + cl->length; 
                    for (k = n_u0, i = n_v0 = n_v = 0; k < n_u; k++) {
                        z = &(res->list[k]);
                        z->non_homopolymer_errors = des_idx + i;
                        n_v0 = z->x_id; ni = z->align_length;
                        for (j = 0; j < ni; j++, i++) {
                            ///k0 + (ni - j - 1)
                            swap[i] = a[ii[n_v0 + (ni- j - 1)]]; 
                            swap[i].readID = k; 
                        }
                        z->x_id = xid; 
                        if(gen_cigar) gen_fake_cigar(&(z->f_cigar), z, apend_be, swap+i-ni, ni);
                        if(!khit_n) z->align_length = 0; 
                    }
                    memcpy(des, swap, i*sizeof((*swap))); //assert(i == ch_n);
                    
                    // fprintf(stderr, "[M::%s::msc->%ld] msc_k_hits::%u, cL::%ld, min_sc::%ld, best_sc::%ld, n_u0_sc::%d, mcopy_rate::%f, # chains::%ld\n", 
                    // __func__, msc, res->list[n_u0].align_length, cL, min_sc, msc+plus, res->list[n_u0].shared_seed,
                    // mcopy_rate, n_u-n_u0);
                } else if(n_u == n_u0 + 1) {
                    z = &(res->list[n_u0]); k = n_u0; i = 0;
                    z->non_homopolymer_errors = des_idx + i;
                    n_v0 = z->x_id; ni = z->align_length;
                    for (j = 0; j < ni; j++, i++) {
                        ///k0 + (ni - j - 1)
                        des[i] = a[ii[n_v0 + (ni- j - 1)]]; 
                        des[i].readID = k; 
                    }
                    z->x_id = xid; 
                    if(gen_cigar) gen_fake_cigar(&(z->f_cigar), z, apend_be, des+i-ni, ni);
                    if(!khit_n) z->align_length = 0; 
                }
                return i;
            } else {
                msc += plus; i = msc_i; cL = 0; 
                while (i >= 0) {t[cL++] = i; i = p[i];}
            }
        }
    }
    ///a[] has been sorted by self_offset
    // i = msc_i; cL = 0; 
    // while (i >= 0) {t[cL++] = i; i = p[i];}
    kv_pushp_ol(overlap_region, (*res), &z);
    push_ovlp_chain_qgen(z, xid, xl, yl, msc, &(a[t[cL-1]]), &(a[t[0]]));
    for (i = 0; i < cL; i++) {des[i] = a[t[cL-i-1]]; des[i].readID = res->length-1;}
    z->non_homopolymer_errors = des_idx;
    if(gen_cigar) gen_fake_cigar(&(z->f_cigar), z, apend_be, des, cL);
    if(khit_n) z->align_length = cL;
    return cL;
}



#define rev_khit(an, xl, yl) do {	\
		(an).self_offset = (xl)-1-((an).self_offset+1-((an).cnt&((uint32_t)(0xffu))));	\
		(an).offset = (yl)-1-((an).offset+1-((an).cnt&((uint32_t)(0xffu))));\
	} while (0)	

///it is unpossibale that (left_fix == 0 && right_fix == 0)
uint64_t lchain_qdp_fix(k_mer_hit* a, int64_t a_n, Chain_Data* dp, int64_t max_skip, 
                int64_t max_iter, int64_t max_dis, double chn_pen_gap, double chn_pen_skip, 
                double bw_rate, int64_t xl, int64_t yl, int64_t quick_check, 
                int64_t left_fix, int64_t right_fix)
{
    int64_t *p, *t, max_f, n_skip, st, max_j, end_j, sc, msc, msc_i, bw, max_ii, ovl, movl; 
    int32_t *f, max, tmp; int64_t i, j, ret, cL = 0, must_p = 1, is_reorder = 0; k_mer_hit z;
    resize_Chain_Data(dp, a_n, NULL);
    t = dp->tmp; f = dp->score; p = dp->pre;
    bw = ((xl < yl)?xl:yl); bw *= bw_rate;
    msc = msc_i = -1; movl = INT32_MAX;

    if(quick_check) {
        ret = lchain_qcheck(a, a_n, dp, bw_rate);
        if (ret > 0) {
            a_n = ret; msc_i = a_n-1; msc = f[msc_i];
            goto skip_ldp;
        }
    }

    memset(t, 0, (a_n*sizeof((*t))));
    if((right_fix) && (!left_fix)) {
        // fprintf(stderr, "\n[M::%s::] a_n::%ld\n", __func__, a_n);
        // for (i = 0; i < a_n; ++i) {
        //     fprintf(stderr, "+[M::%s::] x::[%u, %u), y::[%u, %u)\n", __func__, 
        //     a[i].self_offset+1-(a[i].cnt&((uint32_t)(0xffu))), a[i].self_offset+1,
        //     a[i].offset+1-(a[i].cnt&((uint32_t)(0xffu))), a[i].offset+1);
        // }
        n_skip = a_n>>1;
        for (i=0; i<n_skip; ++i) {
            z = a[i]; a[i] = a[a_n-i-1]; a[a_n-i-1] = z; 
            rev_khit(a[i], xl, yl); rev_khit(a[a_n-i-1], xl, yl);
        }
        if(a_n&1) rev_khit(a[i], xl, yl);
        // for (i = 0; i < a_n; ++i) {
        //     fprintf(stderr, "-[M::%s::] x::[%u, %u), y::[%u, %u)\n", __func__, 
        //     a[i].self_offset+1-(a[i].cnt&((uint32_t)(0xffu))), a[i].self_offset+1,
        //     a[i].offset+1-(a[i].cnt&((uint32_t)(0xffu))), a[i].offset+1);
        // }
        is_reorder = 1;
    } 

    if((!right_fix) && (!left_fix)) must_p = 0;

    for (i = st = 0, max_ii = -1; i < a_n; ++i) {
        max_f = a[i].cnt&(0xffu); if(must_p) max_f = INT32_MIN;
        n_skip = 0; max_j = end_j = -1;
        if ((i-st) > max_iter) st = i-max_iter;
        // if(a[0].self_offset == 171728) {
        //     fprintf(stderr, "i::%ld[M::%s::] x::[%u, %u), y::[%u, %u)\n", i, __func__, 
        //         a[i].self_offset+1-(a[i].cnt&((uint32_t)(0xffu))), a[i].self_offset+1,
        //         a[i].offset+1-(a[i].cnt&((uint32_t)(0xffu))), a[i].offset+1);
        // }
        for (j = i - 1; j >= 0; --j) {
            sc = comput_sc_ff(&a[i], &a[j], bw_rate, chn_pen_gap, chn_pen_skip, xl, yl);
            // if(a[0].self_offset == 171728) {
            //     fprintf(stderr, "j::%ld[M::%s::sc->%ld] x::[%u, %u), y::[%u, %u)\n", j, __func__, sc,
            //     a[j].self_offset+1-(a[j].cnt&((uint32_t)(0xffu))), a[j].self_offset+1,
            //     a[j].offset+1-(a[j].cnt&((uint32_t)(0xffu))), a[j].offset+1);
            // }
            if (sc == INT32_MIN) continue;
            sc += f[j];
            if (sc > max_f) {
                max_f = sc, max_j = j;
                if (n_skip > 0) --n_skip;
            } else if (t[j] == (int32_t)i) {
                if ((++n_skip) > max_skip) {
                    if((max_j != -1) || (must_p == 0)) break;
                } 
            }
            if (p[j] >= 0) t[p[j]] = i;
            ///put it here will allow at least one prefix no matter max_dis
            ///this is special for gap filling, not for chaining
            if (a[i].self_offset > (max_dis + a[j].self_offset)) {
                if((max_j != -1)) break;
            }
            if (j < st) {
                if((max_j != -1) || (must_p == 0)) break;
            }
        }
        end_j = j;

        if (max_ii < 0 || ((int64_t)a[i].self_offset) - ((int64_t)a[max_ii].self_offset) > max_dis) {
            max = INT32_MIN; max_ii = -1;
            for (j = i - 1; (j >= st) && ((((int64_t)a[i].self_offset)-((int64_t)a[j].self_offset))<=max_dis); --j) {
                if (max < f[j]) {
                    max = f[j], max_ii = j;
                }
            }
        }

        if (max_ii >= 0 && max_ii < end_j) {///just have a try with a[i]<->a[max_ii]
            tmp = comput_sc_ff(&a[i], &a[max_ii], bw_rate, chn_pen_gap, chn_pen_skip, xl, yl);
            if (tmp != INT32_MIN && max_f < tmp + f[max_ii])
                max_f = tmp + f[max_ii], max_j = max_ii;
        }
        if(max_j == -1) {
            f[i] = 0; p[i] = max_j;
        } else {
            f[i] = max_f; p[i] = max_j;
        }
        
        if ((max_ii < 0) || (((((int64_t)a[i].self_offset)-((int64_t)a[max_ii].self_offset))<=max_dis) && (f[max_ii]<f[i]))) {
            max_ii = i;
        }
        // if(a[0].self_offset == 171728) {
        //     fprintf(stderr, "i::%ld[M::%s::] f::%d, p::%ld\n", i, __func__, f[i], p[i]);
        // }
        if(f[i] >= msc) {
            ovl = get_chainLen(a[i].self_offset, a[i].self_offset, xl, a[i].offset, a[i].offset, yl);
            if(f[i] > msc || ovl < movl) {
                msc = f[i]; msc_i = i; movl = ovl;
            }
        }
    }
    
    skip_ldp:
    if(right_fix && left_fix) msc_i = a_n-1;
    ///a[] has been sorted by self_offset
    i = msc_i; cL = 0; 
    while (i >= 0) {
        t[cL++] = i; msc_i = i; i = p[i];
    }
    // if((right_fix) && (!left_fix)) {
    //     fprintf(stderr, "[M::%s::] cL::%ld\n", __func__, cL);
    // }

    if(is_reorder) {
        n_skip = a_n>>1;
        for (i=0; i<n_skip; ++i) {
            z = a[i]; a[i] = a[a_n-i-1]; a[a_n-i-1] = z; 
            rev_khit(a[i], xl, yl); rev_khit(a[a_n-i-1], xl, yl);
        }
        if(a_n&1) rev_khit(a[i], xl, yl);
        for (i = 0; i < cL; i++) {
            t[i] = a_n - t[i] - 1;
            // fprintf(stderr, ">[M::%s::] t[%ld]::%ld\n", __func__, i, t[i]);
        }
        // for (i = 0; i < a_n; ++i) {
        //     fprintf(stderr, ">[M::%s::] x::[%u, %u), y::[%u, %u)\n", __func__, 
        //     a[i].self_offset+1-(a[i].cnt&((uint32_t)(0xffu))), a[i].self_offset+1,
        //     a[i].offset+1-(a[i].cnt&((uint32_t)(0xffu))), a[i].offset+1);
        // }
    } else {
        n_skip = cL>>1;
        for (i = 0; i < n_skip; i++) {
            msc_i = t[i]; t[i] = t[cL-i-1]; t[cL-i-1] = msc_i;
        }
    } 
    // if(cL != a_n) {
    //     fprintf(stderr, "\n[M::%s::] a_n::%ld, left_fix::%ld, right_fix::%ld\n", __func__, a_n, left_fix, right_fix);
    //     for (i = 0; i < a_n; ++i) {
    //         fprintf(stderr, "+[M::%s::] x::[%u, %u), y::[%u, %u)\n", __func__, 
    //         a[i].self_offset+1-(a[i].cnt&((uint32_t)(0xffu))), a[i].self_offset+1,
    //         a[i].offset+1-(a[i].cnt&((uint32_t)(0xffu))), a[i].offset+1);
    //     }
    //     for (i = 0; i < cL; i++) {
    //         fprintf(stderr, ">[M::%s::] t[%ld]::%ld\n", __func__, i, t[i]);
    //     }
    // }
    return cL;
}


uint64_t lchain_refine(k_mer_hit* a, int64_t a_n, k_mer_hit* des, Chain_Data* dp, 
                                    int64_t max_skip, int64_t max_iter, int64_t max_dis, int64_t long_gap)
{
    if(a_n <= 0) return 0;
    int64_t *p, *t, max_f, n_skip, st, max_j, sc, msc, msc_i, dq, dr, dd; 
    int32_t *f; int64_t i, j, cL = 0; 
    resize_Chain_Data(dp, a_n, NULL);
    t = dp->tmp; f = dp->score; p = dp->pre; msc = msc_i = -1;

    for (i = 1, f[0] = 0, p[0] = -1, msc_i = a_n - 1; i < a_n; i++) {
        j = i-1;
        dq = (int64_t)(a[i].self_offset) - (int64_t)(a[j].self_offset);
        dr = (int64_t)(a[i].offset) - (int64_t)(a[j].offset);
        dd = dr > dq? dr - dq : dq - dr;//gap
        if(dd <= long_gap || dq > max_dis) {
            p[i] = i - 1; f[i] = i;
        } else {
            break;
        }
    }
    if(i >= a_n) goto ss_kip;
    

    memset(t, 0, (a_n*sizeof((*t))));
    f[0] = 0; p[0] = -1;

    for (i = 1, st = 0; i < a_n; ++i) {
        max_f = INT32_MIN; n_skip = 0; max_j = -1;
        if ((i-st) > max_iter) st = i-max_iter;
        ///i-1
        j = i - 1;
        dq = (int64_t)(a[i].self_offset) - (int64_t)(a[j].self_offset);
        dr = (int64_t)(a[i].offset) - (int64_t)(a[j].offset);
        dd = dr > dq? dr - dq : dq - dr;//gap
        if(dd <= long_gap) dd = 0;
        sc = f[j] - dd;
        if (sc > max_f) {
            max_f = sc, max_j = j;
        } 
        if (p[j] >= 0) t[p[j]] = i;
        
        ///[st, i-2]
        for (--j; (j >= st) && (a[i].self_offset <= (max_dis + a[j].self_offset)); --j) {
            dq = (int64_t)(a[i].self_offset) - (int64_t)(a[j].self_offset);
            dr = (int64_t)(a[i].offset) - (int64_t)(a[j].offset);
            dd = dr > dq? dr - dq : dq - dr;//gap
            if(dd <= long_gap) dd = 0;
            sc = f[j] - dd;
            if (sc > max_f) {
                max_f = sc, max_j = j;
                if (n_skip > 0) --n_skip;
            } else if (t[j] == (int32_t)i) {
                if (++n_skip > max_skip)
                    break;
            }
            if (p[j] >= 0) t[p[j]] = i;
        }

        f[i] = max_f; p[i] = max_j;
    }

    i = a_n-1; msc = f[i]; msc_i = i; 
    for (j = i-1; (j >= 0) && (a[i].self_offset <= (max_dis + a[j].self_offset)); --j) {
        if(msc < f[j] && p[j] >= 0) {///hold at least two hits in th final chain
            msc = f[j]; msc_i = j;
        }
    }
    
    ss_kip:
    ///a[] has been sorted by self_offset
    i = msc_i; 
    cL = 0; 
    while (i >= 0) {
        t[cL++] = i; i = p[i];
    }

    n_skip = cL>>1;
    for (i = 0; i < n_skip; i++) {
        msc_i = t[i]; t[i] = t[cL-i-1]; t[cL-i-1] = msc_i;
    }
    if(des) {
        for (i = 0; i < cL; i++) des[i] = a[t[i]];
    }
    return cL;
}


uint64_t lchain_simple(k_mer_hit* a, int64_t a_n, k_mer_hit* des, Chain_Data* dp, 
                                                        int64_t max_skip, int64_t max_iter)
{
    if(a_n <= 0) return 0;
    int64_t *p, *t, max_f, n_skip, st, max_j, sc, msc, msc_i; 
    int32_t *f; int64_t i, j, cL = 0; 
    resize_Chain_Data(dp, a_n, NULL);
    t = dp->tmp; f = dp->score; p = dp->pre; msc = msc_i = -1;

    for (i=1, f[0]=a[0].cnt, p[0]=-1, msc_i=a_n-1; i<a_n; i++) {
        j = i-1;
        if((a[i].self_offset > a[j].self_offset)&&(a[i].offset > a[j].offset)) {
            p[i] = j; f[i] = f[j]+a[i].cnt;
        } else {
            break;
        }
    }
    if(i >= a_n) goto ss_kip;
    
    memset(t, 0, (a_n*sizeof((*t))));
    f[0]=a[0].cnt; p[0]=-1; msc = f[0]; msc_i = 0;

    for (i = 1, st = 0; i < a_n; ++i) {
        max_f = INT32_MIN; n_skip = 0; max_j = -1;
        if ((i-st) > max_iter) st = i-max_iter;
        ///[st, i-2]
        for (j=i-1; j >= st; --j) {
            if((a[i].self_offset > a[j].self_offset)&&(a[i].offset > a[j].offset)) {
                sc = f[j]+a[i].cnt;
                if (sc > max_f) {
                    max_f = sc, max_j = j;
                    if (n_skip > 0) --n_skip;
                } else if (t[j] == (int32_t)i) {
                    if (++n_skip > max_skip)
                        break;
                }
                if (p[j] >= 0) t[p[j]] = i;
            }
        }
        f[i] = max_f; p[i] = max_j;
        if(f[i] > msc) {
            msc = f[i]; msc_i = i;
        }
    }

    ss_kip:
    ///a[] has been sorted by self_offset
    i = msc_i; 
    cL = 0; 
    while (i >= 0) {
        t[cL++] = i; i = p[i];
    }

    n_skip = cL>>1;
    for (i = 0; i < n_skip; i++) {
        msc_i = t[i]; t[i] = t[cL-i-1]; t[cL-i-1] = msc_i;
    }
    if(des) {
        for (i = 0; i < cL; i++) des[i] = a[t[i]];
    }
    return cL;
}

inline int64_t hit_long_gap(k_mer_hit *a, k_mer_hit *b, int64_t max_lgap, double small_bw_rate, int64_t min_small_bw)
{
    int64_t dq, dr, dd, dm;
    dq = b->self_offset-a->self_offset;
    dr = b->offset-a->offset;
	dd = dq>=dr? ((dq)-(dr)): ((dr)-(dq));
    if(max_lgap>=0) {
        if(dd <= max_lgap) return 1;
        return 0;
    } else {
        dm = dq>=dr?dr:dq;
        if((dd > (dm*small_bw_rate)) && (dd > min_small_bw)) return 0;
        return 1;
    }
}

int64_t filter_bad_seed_dp(k_mer_hit *sk, k_mer_hit *ek, k_mer_hit* a, int64_t a_n, int64_t max_lgap, double small_bw_rate, int64_t min_small_bw)
{
    int64_t k = 0; k_mer_hit *z; double bw_r; int64_t bw, mmgap, occ = 0;
    bw_r = small_bw_rate; bw = min_small_bw; mmgap = max_lgap;
    if(sk) {
        for (k = 0; k < a_n; k++) {
            z = &(a[k]); 
            if(z->cnt < z->readID) continue;
            mmgap = max_lgap;
            if(z->cnt <= 1) mmgap = -1;
            if((sk && (!hit_long_gap(sk, z, mmgap, bw_r, bw))) || 
                                (ek && (!hit_long_gap(z, ek, mmgap, bw_r, bw)))) {
                z->offset = z->self_offset = (uint32_t)-1; occ++;
            } else {
                break;
            }
        }
    }

    if(ek && k < a_n) {
        for (k = a_n-1; k >= 0; k--) {
            z = &(a[k]); 
            if(z->cnt < z->readID) continue;
            mmgap = max_lgap;
            if(z->cnt <= 1) mmgap = -1;
            if((sk && (!hit_long_gap(sk, z, mmgap, bw_r, bw))) || 
                                (ek && (!hit_long_gap(z, ek, mmgap, bw_r, bw)))) {
                z->offset = z->self_offset = (uint32_t)-1; occ++;
            } else {
                break;
            }
        }
    }

    return occ;    
}

uint64_t lchain_dp_trace(k_mer_hit* a, int64_t a_n, int64_t max_lgap, double sgap_rate, int64_t sgap)
{
    // fprintf(stderr, "[M::%s::] a_n::%ld\n", __func__, a_n);
    if(a_n <= 0) return 0;
    int64_t i, st, occ = 0; 

    for (i = 1, st = 0; i <= a_n; ++i) {
        // fprintf(stderr, "[M::%s::i->%ld] q::%u, t::%u, cnt::%u, readID::%u\n", __func__, 
        // i-1, a[i-1].self_offset, a[i-1].offset, a[i-1].cnt, a[i-1].readID);
        if((i == a_n) || (is_alnw(a[i]))) {///[st, i)
            if(i > st) {
                occ += filter_bad_seed_dp((st>0)?&(a[st-1]):NULL, (i<a_n)?&(a[i]):NULL, a, a_n, max_lgap, sgap_rate, sgap);
            }
            st = i+1;
        }
    }

    if(occ) {
        for (i = occ = 0; i < a_n; ++i) {
            if(a[i].self_offset == ((uint32_t)-1)) continue;
            a[occ++] = a[i];
        }
        a_n = occ;
    }
    return a_n;
}