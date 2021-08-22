#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include "Hash_Table.h"
#include "ksort.h"
pthread_mutex_t output_mutex;

#define overlap_region_key(a) ((a).y_id)
KRADIX_SORT_INIT(overlap_region_sort, overlap_region, overlap_region_key, member_size(overlap_region, y_id))

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
    for (i = 0; i < list->size; i++)
    { 
        init_fake_cigar(&(list->list[i].f_cigar));
        init_window_list_alloc(&(list->list[i].boundary_cigars));
    }
}

void clear_overlap_region_alloc(overlap_region_alloc* list)
{
    list->length = 0;
    list->mapped_overlaps_length = 0;
    uint64_t i = 0;
    for (i = 0; i < list->size; i++)
    {   
        list->list[i].w_list_length = 0;
        clear_fake_cigar(&(list->list[i].f_cigar));
        clear_window_list_alloc(&(list->list[i].boundary_cigars));
    }
}

void destory_overlap_region_alloc(overlap_region_alloc* list)
{
    uint64_t i = 0;
    for (i = 0; i < list->size; i++)
    {
        if (list->list[i].w_list_size != 0)
        {
            free(list->list[i].w_list);
        }
        destory_fake_cigar(&(list->list[i].f_cigar));
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

int append_inexact_overlap_region_alloc(overlap_region_alloc* list, overlap_region* tmp, 
                                        All_reads* R_INF, int add_beg_end)
{
   
    if (list->length + 1 > list->size)
    {
        list->size = list->size * 2;
        list->list = (overlap_region*)realloc(list->list, sizeof(overlap_region)*list->size);
        /// need to set new space to be 0
        memset(list->list + (list->size/2), 0, sizeof(overlap_region)*(list->size/2));
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


    long long x_right_length = Get_READ_LENGTH((*R_INF), tmp->x_id) - tmp->x_pos_e - 1;
    long long y_right_length = Get_READ_LENGTH((*R_INF), tmp->y_id) - tmp->y_pos_e - 1;

    if(x_right_length <= y_right_length)
    {
        tmp->x_pos_e = Get_READ_LENGTH((*R_INF), tmp->x_id) - 1;
        tmp->y_pos_e = tmp->y_pos_e + x_right_length;        
    }
    else
    {
        tmp->x_pos_e = tmp->x_pos_e + y_right_length;
        tmp->y_pos_e = Get_READ_LENGTH((*R_INF), tmp->y_id) - 1;
    }
    
    if (tmp->x_pos_strand == 1)
    {
        list->list[list->length].x_id = tmp->x_id;
        list->list[list->length].x_pos_e = Get_READ_LENGTH((*R_INF), tmp->x_id) - tmp->x_pos_s - 1;
        list->list[list->length].x_pos_s = Get_READ_LENGTH((*R_INF), tmp->x_id) - tmp->x_pos_e - 1;
        list->list[list->length].x_pos_strand = 0;

        list->list[list->length].y_id = tmp->y_id;
        list->list[list->length].y_pos_e = Get_READ_LENGTH((*R_INF), tmp->y_id) - tmp->y_pos_s - 1;
        list->list[list->length].y_pos_s = Get_READ_LENGTH((*R_INF), tmp->y_id) - tmp->y_pos_e - 1;
        list->list[list->length].y_pos_strand = 1;

        resize_fake_cigar(&(list->list[list->length].f_cigar), (tmp->f_cigar.length + 2));
        if(add_beg_end == 1)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), list->list[list->length].x_pos_s, 0);
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
                Get_READ_LENGTH((*R_INF), tmp->x_id) - get_fake_gap_pos(&(tmp->f_cigar), i) - 1, 
                pre_distance_gap);
            }
        }

        if(add_beg_end == 1 && get_fake_gap_pos(&(list->list[list->length].f_cigar), 
            list->list[list->length].f_cigar.length - 1) != (long long)list->list[list->length].x_pos_e)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), 
            list->list[list->length].x_pos_e, 
            get_fake_gap_shift(&(list->list[list->length].f_cigar), 
            list->list[list->length].f_cigar.length - 1));
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



        resize_fake_cigar(&(list->list[list->length].f_cigar), (tmp->f_cigar.length + 2));
        if(add_beg_end == 1)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), list->list[list->length].x_pos_s, 0);
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
                get_fake_gap_pos(&(tmp->f_cigar), i), init_distance_gap - pre_distance_gap);
            }
        }

        if(add_beg_end == 1 && get_fake_gap_pos(&(list->list[list->length].f_cigar), 
        list->list[list->length].f_cigar.length - 1) != (long long)list->list[list->length].x_pos_e)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), 
            list->list[list->length].x_pos_e, 
            get_fake_gap_shift(&(list->list[list->length].f_cigar), 
            list->list[list->length].f_cigar.length - 1));
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
                                        ma_utg_v *ua, int add_beg_end)
{
   
    if (list->length + 1 > list->size)
    {
        list->size = list->size * 2;
        list->list = (overlap_region*)realloc(list->list, sizeof(overlap_region)*list->size);
        /// need to set new space to be 0
        memset(list->list + (list->size/2), 0, sizeof(overlap_region)*(list->size/2));
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

        resize_fake_cigar(&(list->list[list->length].f_cigar), (tmp->f_cigar.length + 2));
        if(add_beg_end == 1)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), list->list[list->length].x_pos_s, 0);
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
                pre_distance_gap);
            }
        }

        if(add_beg_end == 1 && get_fake_gap_pos(&(list->list[list->length].f_cigar), 
            list->list[list->length].f_cigar.length - 1) != (long long)list->list[list->length].x_pos_e)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), 
            list->list[list->length].x_pos_e, 
            get_fake_gap_shift(&(list->list[list->length].f_cigar), 
            list->list[list->length].f_cigar.length - 1));
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



        resize_fake_cigar(&(list->list[list->length].f_cigar), (tmp->f_cigar.length + 2));
        if(add_beg_end == 1)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), list->list[list->length].x_pos_s, 0);
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
                get_fake_gap_pos(&(tmp->f_cigar), i), init_distance_gap - pre_distance_gap);
            }
        }

        if(add_beg_end == 1 && get_fake_gap_pos(&(list->list[list->length].f_cigar), 
        list->list[list->length].f_cigar.length - 1) != (long long)list->list[list->length].x_pos_e)
        {
            add_fake_cigar(&(list->list[list->length].f_cigar), 
            list->list[list->length].x_pos_e, 
            get_fake_gap_shift(&(list->list[list->length].f_cigar), 
            list->list[list->length].f_cigar.length - 1));
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
              double band_width_threshold, int max_skip, int x_readLen, int y_readLen)
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
    
    resize_Chain_Data(dp, a_n);

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
    add_fake_cigar(&(result->f_cigar), a[i].self_offset, pre_distance_gap);
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
                add_fake_cigar(&(result->f_cigar), a[i].self_offset, pre_distance_gap);
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
                add_fake_cigar(&(result->f_cigar), a[i].self_offset, pre_distance_gap);
            }
            else
            {
                pre_distance_gap = distance_gap;
                add_fake_cigar(&(result->f_cigar), a[i].self_offset, pre_distance_gap);
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
                                          uint64_t readID, uint64_t readLength, All_reads* R_INF, double band_width_threshold, int add_beg_end, overlap_region* f_cigar)
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
				25, Get_READ_LENGTH((*R_INF), (*f_cigar).x_id), Get_READ_LENGTH((*R_INF), (*f_cigar).y_id));

        ///if (tmp_region.x_id != tmp_region.y_id && tmp_region.shared_seed > 1)
        if ((*f_cigar).x_id != (*f_cigar).y_id)
        {
            append_inexact_overlap_region_alloc(overlap_list, f_cigar, R_INF, add_beg_end);
        }
    }
}


void append_window_list(overlap_region* region, uint64_t x_start, uint64_t x_end, int y_start, int y_end, int error,
                        int extra_begin, int extra_end, int error_threshold)
{
    
    long long length = region->x_pos_e - region->x_pos_s + 1;
    ///the length of window may large or small than WINDOW
    /****************************may have bugs********************************/
    uint64_t num_windows = length / WINDOW + 4;
    /****************************may have bugs********************************/

    ///w_list_length has alredy set to be 0 at clear_overlap_region_alloc
    if (num_windows > region->w_list_size)
    {
        region->w_list_size = num_windows;
        region->w_list = (window_list*)realloc(region->w_list, region->w_list_size*sizeof(window_list));
    }
    

    region->w_list[region->w_list_length].x_start = x_start;
    region->w_list[region->w_list_length].x_end = x_end;
    region->w_list[region->w_list_length].y_start = y_start;
    region->w_list[region->w_list_length].y_end = y_end;
    region->w_list[region->w_list_length].error = error;
    region->w_list[region->w_list_length].cigar.length = -1;
    region->w_list[region->w_list_length].extra_begin = extra_begin;
    region->w_list[region->w_list_length].extra_end = extra_end;
    region->w_list[region->w_list_length].error_threshold = error_threshold;
    region->w_list_length++;
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

void resize_Chain_Data(Chain_Data* x, long long size)
{
	if (size + 1 > x->size) {
		x->size = size + 1;
		kroundup64(x->size);
		REALLOC(x->score, x->size);
		REALLOC(x->pre, x->size);
		REALLOC(x->indels, x->size);
		REALLOC(x->self_length, x->size);
        REALLOC(x->occ, x->size);
		REALLOC(x->tmp, x->size);
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

void add_fake_cigar(Fake_Cigar* x, uint32_t gap_site, int32_t gap_shift)
{
    if(x->length + 1 > x->size)
    {
        x->size = (x->length + 1) * 2;
        x->buffer = (uint64_t*)realloc(x->buffer, sizeof(uint64_t) * x->size);
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


void resize_fake_cigar(Fake_Cigar* x, uint64_t size)
{
    if (size > x->size) {
        x->size = size;
		REALLOC(x->buffer, x->size);
    }
    x->length = 0;
}


void init_window_list_alloc(window_list_alloc* x)
{
    x->buffer = NULL;
    x->length = 0;
    x->size = 0;
}

void clear_window_list_alloc(window_list_alloc* x)
{
    x->length = 0;
}

void destory_window_list_alloc(window_list_alloc* x)
{
    if(x->size != 0)
    {
        free((x->buffer));
    }
}

void resize_window_list_alloc(window_list_alloc* x, long long size)
{
    if(size > x->size)
    {
        x->size = size;
        x->buffer = (window_list*)realloc(x->buffer, sizeof(window_list) * x->size);
    }

    long long i;
    for (i = 0; i < x->size; i++)
    {
        x->buffer[i].error = -1; 
    }
    x->length = 0;
}
