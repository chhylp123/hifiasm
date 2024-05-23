#define __STDC_LIMIT_MACROS
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include <assert.h>
#include <float.h>
#include "Correct.h"
#include "Levenshtein_distance.h"
#include "Assembly.h"
#include "CommandLines.h"
// #include "ksw2.h"
#include "ksort.h"
#include "kalloc.h"
#include "htab.h"
#include "Overlaps.h"
#include "inter.h"
#define A_L 16
#define ext_w 6
#define r_simi_w 0.05
#define rphase_thres 4

#define generic_key(x) (x)
KRADIX_SORT_INIT(b32, uint32_t, generic_key, 4)
KRADIX_SORT_INIT(bc64, uint64_t, generic_key, 8)

#define haplotype_evdience_key(x) ((x).site)
KRADIX_SORT_INIT(haplotype_evdience_srt, haplotype_evdience, haplotype_evdience_key, member_size(haplotype_evdience, site))

#define haplotype_evdience_id_key(x) ((x).overlapID)
KRADIX_SORT_INIT(haplotype_evdience_id_srt, haplotype_evdience, haplotype_evdience_id_key, member_size(haplotype_evdience, overlapID))

#define overlap_region_dp_key(x) ((x).x_pos_e)
KRADIX_SORT_INIT(overlap_region_dp_srt, overlap_region, overlap_region_dp_key, member_size(overlap_region, x_pos_e))

#define window_list_xs_key(x) ((x).x_start)
KRADIX_SORT_INIT(window_list_xs_srt, window_list, window_list_xs_key, member_size(window_list, x_start))

#define uov_qs_key(p) ((p).qs)
KRADIX_SORT_INIT(uov_srt_qs, ul_ov_t, uov_qs_key, member_size(ul_ov_t, qs))

#define k_mer_hit_self_key(p) ((p).self_offset)
KRADIX_SORT_INIT(k_mer_hit_self, k_mer_hit, k_mer_hit_self_key, member_size(k_mer_hit, self_offset))

#define k_mer_hit_off_key(p) ((p).offset)
KRADIX_SORT_INIT(k_mer_hit_off, k_mer_hit, k_mer_hit_off_key, member_size(k_mer_hit, offset))

#define ul_ov_srt_qs1_key(p) ((p).qs)
KRADIX_SORT_INIT(ul_ov_srt_qs1, ul_ov_t, ul_ov_srt_qs1_key, member_size(ul_ov_t, qs))

#define ul_ov_srt_tn1_key(p) ((p).tn)
KRADIX_SORT_INIT(ul_ov_srt_tn1, ul_ov_t, ul_ov_srt_tn1_key, member_size(ul_ov_t, tn))

#define ul_ov_srt_qn1_key(p) ((p).qn)
KRADIX_SORT_INIT(ul_ov_srt_qn1, ul_ov_t, ul_ov_srt_qn1_key, member_size(ul_ov_t, qn))

#define ul_ov_srt_qe1_key(p) ((p).qe)
KRADIX_SORT_INIT(ul_ov_srt_qe1, ul_ov_t, ul_ov_srt_qe1_key, member_size(ul_ov_t, qe))

#define MAX_SEC_ERR (0x3fffffffU)


int ha_ov_type(const overlap_region *r, uint32_t len);
void set_lchain_dp_op(uint32_t is_accurate, uint32_t mz_k, int64_t *max_skip, int64_t *max_iter, int64_t *max_dis, double *chn_pen_gap, double *chn_pen_skip, int64_t *quick_check);

void clear_Round2_alignment(Round2_alignment* h)
{
    clear_Correct_dumy_pure(&(h->dumy));
    clear_Cigar_record(&(h->cigar));
    clear_Cigar_record(&(h->tmp_cigar));
    h->obtained_cigar_length = 0;
}

void init_Round2_alignment(Round2_alignment* h)
{
    init_Correct_dumy(&(h->dumy));
    init_Cigar_record(&(h->cigar));
    init_Cigar_record(&(h->tmp_cigar));
    h->obtained_cigar_length = 0;
}

void init_Round2_alignment_buf(Round2_alignment* h, void *km)
{
    init_Correct_dumy_buf(&(h->dumy), km);
    init_Cigar_record_buf(&(h->cigar), km);
    init_Cigar_record_buf(&(h->tmp_cigar), km);
    h->obtained_cigar_length = 0;
}

void destory_Round2_alignment(Round2_alignment* h)
{
    destory_Correct_dumy(&(h->dumy));
    destory_Cigar_record(&(h->cigar));
    destory_Cigar_record(&(h->tmp_cigar));
}



inline int get_interval(long long window_start, long long window_end, overlap_region_alloc* overlap_list, Correct_dumy* dumy, long long blockLen)
{
    uint64_t i, fud = 0;
    long long Len;
    if(window_start == 0) dumy->start_i = 0;
    for (i = dumy->start_i; i < overlap_list->length; i++)
    {
        ///this interval is smaller than all overlaps
        ///in this case, the next interval should start from 0
        if (window_end < (long long)overlap_list->list[i].x_pos_s)
        {
            dumy->start_i = 0;
            dumy->length = 0;
            dumy->lengthNT = 0;
            return 0;
        }
        else ///if window_end >= overlap_list->list[i].x_pos_s，this overlap might be overlapped with current interval
        {
            dumy->start_i = i;
            break;
        }
    }


    ///this interval is larger than all overlaps, so we don't need to scan next overlap
    if (i >= overlap_list->length)
    {
        dumy->start_i = overlap_list->length;
        dumy->length = 0;
        dumy->lengthNT = 0;
        return -2;
    }
    
    dumy->length = 0;
    dumy->lengthNT = 0;
    fud = 0;

    for (; i < overlap_list->length; i++)
    {
        if((Len = OVERLAP(window_start, window_end, (long long)overlap_list->list[i].x_pos_s, (long long)overlap_list->list[i].x_pos_e)) > 0)
        {
            ///sometimes the length of window > WINDOW, but overlap length == WINDOW
            // if (Len == WINDOW && window_end - window_start + 1 == WINDOW)
            if (Len == blockLen && window_end - window_start + 1 == blockLen)
            {
                dumy->overlapID[dumy->length] = i;
                dumy->length++; 
            }
            else
            {
                dumy->lengthNT++;
                dumy->overlapID[dumy->size - dumy->lengthNT] = i;
            }
            if(fud == 0) fud = 1, dumy->start_i = i;
        }
        
        if((long long)overlap_list->list[i].x_pos_s > window_end)
        {
            break;
        }
    }

    if ( dumy->length + dumy->lengthNT == 0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}


inline int get_available_interval(long long window_start, long long window_end, overlap_region_alloc* overlap_list, Correct_dumy* dumy)
{
    uint64_t i, fud = 0;
    long long Len;
    if(window_start == 0) dumy->start_i = 0;
    for (i = dumy->start_i; i < overlap_list->length; i++)
    {
        ///this interval is smaller than all overlaps
        ///in this case, the next interval should start from 0
        if (window_end < (long long)overlap_list->list[i].x_pos_s)
        {
            dumy->start_i = 0;
            dumy->length = 0;
            dumy->lengthNT = 0;
            return 0;
        }
        else ///if window_end >= overlap_list->list[i].x_pos_s，this overlap might be overlapped with current interval
        {
            dumy->start_i = i;
            break;
        }
    }

    ///this interval is larger than all overlaps, so we don't need to scan next overlap
    if (i >= overlap_list->length)
    {
        dumy->start_i = overlap_list->length;
        dumy->length = 0;
        dumy->lengthNT = 0;
        return -2;
    }
    
    dumy->length = 0;
    dumy->lengthNT = 0;
    fud = 0;

    
    long long fake_length = 0;

    for (; i < overlap_list->length; i++)
    {
        ///check if the interval is overlapped with current overlap   
        if((Len = OVERLAP(window_start, window_end, (long long)overlap_list->list[i].x_pos_s, (long long)overlap_list->list[i].x_pos_e)) > 0)
        {
            ///number of overlaps
            fake_length++;

            ///check if this overlap is available
            if (overlap_list->list[i].is_match == 1)
            {
                dumy->overlapID[dumy->length] = i;
                dumy->length++; 
            }
            if(fud == 0) fud = 1, dumy->start_i = i;
        }
        
        if((long long)overlap_list->list[i].x_pos_s > window_end)
        {
            break;
        }
    }

    ///fake_length is the number of overlaps, instead of the number of available overlaps
    if (fake_length == 0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

///Len = OVERLAP(window_start, window_end, overlap_list->list[i].x_pos_s, overlap_list->list[i].x_pos_e))

void print_string(char* s, int l)
{
    int i;
    for (i = 0; i < l; i++)
    {
        fprintf(stderr, "%c", s[i]);
    }

    fprintf(stderr, "\n");
    
}




void fill_subregion(char* r, long long start_pos, long long length, uint8_t strand, All_reads* R_INF, long long ID, 
int extra_begin, int extra_end)
{
    
    recover_UC_Read_sub_region(r+extra_begin, start_pos, length, strand, R_INF, ID);
    memset(r, 'N', extra_begin);
    memset(r+extra_begin+length, 'N', extra_end);
}

void fill_subregion_ul(char* r, long long start_pos, long long length, uint8_t strand, const ul_idx_t *uref, long long ID, 
int extra_begin, int extra_end)
{
    retrieve_u_seq(NULL, r+extra_begin, &(uref->ug->u.a[ID]), strand, start_pos, length, NULL);
    memset(r, 'N', extra_begin);
    memset(r+extra_begin+length, 'N', extra_end);
}

int determine_overlap_region(int threshold, long long y_start, long long y_ID, long long Window_Len, /**All_reads* R_INF**/long long y_len,
int* r_extra_begin, int* r_extra_end, long long* r_y_start, long long* r_y_length)
{
    int extra_begin;
    int extra_end;
    long long currentIDLen;
    long long o_len;

    ///the length of y
    // currentIDLen = Get_READ_LENGTH((*R_INF), y_ID);
    currentIDLen = y_len;
    
    ///since Window_Len == x_len + (threshold << 1)
    if(y_start < 0 || currentIDLen <= y_start || 
    currentIDLen - y_start + 2 * threshold + THRESHOLD_MAX_SIZE < Window_Len)
    {
        return 0;
    }
    
    extra_begin = extra_end = 0;
    ///y maybe less than 0
    y_start = y_start - threshold;
    o_len = MIN(Window_Len, currentIDLen - y_start);
    extra_end = Window_Len - o_len;

    if (y_start < 0)
    {
        extra_begin = -y_start;
        y_start = 0;
        o_len = o_len - extra_begin;
    }

    (*r_extra_begin) = extra_begin;
    (*r_extra_end) = extra_end;
    (*r_y_start) = y_start;
    (*r_y_length) = o_len;

    return 1;
}


int verify_single_window(long long x_start, long long x_end, 
long long overlap_x_s, long long overlap_y_s, int x_id,
int y_id, int x_strand, char* x_buffer, char* y_buffer,
All_reads* R_INF)
{
    char* x_string = NULL;
    char* y_string = NULL;
    int extra_begin, extra_end, x_len, threshold;
    long long y_start;
    long long Window_Len, o_len;
    unsigned int error;

    

    x_len = x_end - x_start + 1;
    threshold = x_len * asm_opt.max_ov_diff_ec;
    /****************************may have bugs********************************/
    threshold = Adjust_Threshold(threshold, x_len);
    /****************************may have bugs********************************/
        
    y_start = (x_start - overlap_x_s) + overlap_y_s;
    
    Window_Len = x_len + (threshold << 1);

    if(!determine_overlap_region(threshold, y_start, y_id, Window_Len, Get_READ_LENGTH((*R_INF), y_id),
    &extra_begin, &extra_end, &y_start, &o_len))
    {
        return 0;
    }

    ///use unusual direction here
    /**
    fill_subregion(y_buffer, y_start, o_len, y_strand, R_INF, y_id, extra_begin, extra_end);
    ///x is always forward strand
    recover_UC_Read_sub_region(x_buffer, x_start, x_len, 0, R_INF, x_id);
    **/
   ///use unusual direction here, here y is always forward strand
    fill_subregion(y_buffer, y_start, o_len, 0, R_INF, y_id, extra_begin, extra_end);
    recover_UC_Read_sub_region(x_buffer, x_start, x_len, x_strand, R_INF, x_id);


    x_string = x_buffer;
    y_string = y_buffer;

    Reserve_Banded_BPM(y_string, Window_Len, x_string, x_len, threshold, &error);

    if (error!=(unsigned int)-1)
    {
        return 1;
    }

    return 0;
}

void verify_window(long long window_start, long long window_end, overlap_region_alloc* overlap_list,Correct_dumy* dumy, All_reads* R_INF,
char* r_string)
{
    long long i;
    long long currentID;
    long long x_start, y_start, o_len;
    long long Window_Len = WINDOW + (THRESHOLD << 1);
    char* x_string = NULL;
    char* y_string = NULL;
    long long x_end, x_len;
    int end_site;
    unsigned int error;
    int groupLen = 0;
    int return_sites[GROUP_SIZE];
	unsigned int return_sites_error[GROUP_SIZE];
    uint64_t overlapID[GROUP_SIZE];
    uint64_t y_startGroup[GROUP_SIZE];
    int y_extra_begin[GROUP_SIZE];
    int y_extra_end[GROUP_SIZE];
    int error_threshold[GROUP_SIZE];
    int extra_begin;
    int extra_end;

    ///here are overlaps fully covered by WINDOW
    for (i = 0; i < (long long)dumy->length; i++)
    {
        extra_begin = extra_end = 0;
        ///if the window has been fully covered, the interval at x is [window_start, window_end]
        x_len = WINDOW;
        currentID = dumy->overlapID[i];
        x_start = window_start;
        ///offset of y
        y_start = (x_start - overlap_list->list[currentID].x_pos_s) + overlap_list->list[currentID].y_pos_s;
        /****************************may have bugs********************************/
        y_start += y_start_offset(x_start, &(overlap_list->list[currentID].f_cigar));
        /****************************may have bugs********************************/
        
  
        if(!determine_overlap_region(THRESHOLD, y_start, overlap_list->list[currentID].y_id, Window_Len, Get_READ_LENGTH((*R_INF), overlap_list->list[currentID].y_id),
        &extra_begin, &extra_end, &y_start, &o_len)) {
            continue;
        }

        fill_subregion(dumy->overlap_region_group[groupLen], y_start, o_len, overlap_list->list[currentID].y_pos_strand, 
        R_INF, overlap_list->list[currentID].y_id, extra_begin, extra_end);

        y_extra_begin[groupLen] = extra_begin;
        y_extra_end[groupLen] = extra_end;
        overlapID[groupLen] = currentID;
        y_startGroup[groupLen] = y_start;
        error_threshold[groupLen] = THRESHOLD;
        x_string = r_string + x_start;
        groupLen++;
        

        if (groupLen == GROUP_SIZE)
        {
            Reserve_Banded_BPM_4_SSE_only(dumy->overlap_region_group[0], dumy->overlap_region_group[1], 
            dumy->overlap_region_group[2], dumy->overlap_region_group[3], Window_Len, x_string, WINDOW,
                return_sites, return_sites_error, THRESHOLD, dumy->Peq_SSE);
            groupLen = 0;

            
            if (return_sites_error[0]!=(unsigned int)-1) {
                overlap_list->list[overlapID[0]].align_length += x_len;

                append_window_list(&overlap_list->list[overlapID[0]], window_start, window_end, 
                            y_startGroup[0], y_startGroup[0] + return_sites[0], (int)return_sites_error[0],
                            y_extra_begin[0], y_extra_end[0], error_threshold[0], WINDOW, NULL);
            }
            

            if (return_sites_error[1]!=(unsigned int)-1) {
                overlap_list->list[overlapID[1]].align_length += x_len;

                append_window_list(&overlap_list->list[overlapID[1]], window_start, window_end, 
                            y_startGroup[1], y_startGroup[1] + return_sites[1], (int)return_sites_error[1],
                            y_extra_begin[1], y_extra_end[1], error_threshold[1], WINDOW, NULL);
            }
            

            if (return_sites_error[2]!=(unsigned int)-1) {
                overlap_list->list[overlapID[2]].align_length += x_len;

                append_window_list(&overlap_list->list[overlapID[2]], window_start, window_end, 
                            y_startGroup[2], y_startGroup[2] + return_sites[2], (int)return_sites_error[2],
                            y_extra_begin[2], y_extra_end[2], error_threshold[2], WINDOW, NULL);
            }
            

            if (return_sites_error[3]!=(unsigned int)-1) {
                overlap_list->list[overlapID[3]].align_length += x_len;

                append_window_list(&overlap_list->list[overlapID[3]], window_start, window_end, 
                            y_startGroup[3], y_startGroup[3] + return_sites[3], (int)return_sites_error[3],
                            y_extra_begin[3], y_extra_end[3], error_threshold[3], WINDOW, NULL);
            }  
        }
    }

    if (groupLen == 1)
    {
        end_site = Reserve_Banded_BPM(dumy->overlap_region_group[0], Window_Len, x_string, WINDOW, THRESHOLD, &error);

        if (error!=(unsigned int)-1) {
            overlap_list->list[overlapID[0]].align_length += x_len;

            append_window_list(&overlap_list->list[overlapID[0]], window_start, window_end, 
                                y_startGroup[0], y_startGroup[0] + end_site, (int)error,
                                y_extra_begin[0], y_extra_end[0], error_threshold[0], WINDOW, NULL);
        }
    }
    else if (groupLen > 1)
    {
        Reserve_Banded_BPM_4_SSE_only(dumy->overlap_region_group[0], dumy->overlap_region_group[1], 
                dumy->overlap_region_group[2], dumy->overlap_region_group[3], Window_Len, x_string, WINDOW,
			        return_sites, return_sites_error, THRESHOLD, dumy->Peq_SSE);

        for (i = 0; i < groupLen; i++)
        {
            if (return_sites_error[i]!=(unsigned int)-1) {
                overlap_list->list[overlapID[i]].align_length += x_len;
                append_window_list(&overlap_list->list[overlapID[i]], window_start, window_end, 
                                y_startGroup[i], y_startGroup[i] + return_sites[i], (int)return_sites_error[i],
                                y_extra_begin[i], y_extra_end[i], error_threshold[i], WINDOW, NULL);
            }
        }

        groupLen = 0;
    }

    long long reverse_i = dumy->size - 1;
    int threshold;
    
    ///here are overlaps partially covered by WINDOW
    for (i = 0; i < (long long)dumy->lengthNT; i++)
    {
        extra_begin = extra_end = 0;
        currentID = dumy->overlapID[reverse_i--];
        x_start = MAX(window_start, (long long)overlap_list->list[currentID].x_pos_s);
        x_end = MIN(window_end, (long long)overlap_list->list[currentID].x_pos_e);

        ///overlap length between [window_start, window_end]
        x_len = x_end - x_start + 1;
        threshold = x_len * asm_opt.max_ov_diff_ec;
        /****************************may have bugs********************************/
        threshold = Adjust_Threshold(threshold, x_len);
        if(threshold > THRESHOLD_MAX_SIZE) threshold = THRESHOLD_MAX_SIZE;
        /****************************may have bugs********************************/
        
        ///offset of y
        y_start = (x_start - overlap_list->list[currentID].x_pos_s) + overlap_list->list[currentID].y_pos_s;
        /****************************may have bugs********************************/
        y_start += y_start_offset(x_start, &(overlap_list->list[currentID].f_cigar));
        /****************************may have bugs********************************/

        Window_Len = x_len + (threshold << 1);

        if(!determine_overlap_region(threshold, y_start, overlap_list->list[currentID].y_id, Window_Len, Get_READ_LENGTH((*R_INF), overlap_list->list[currentID].y_id),
        &extra_begin, &extra_end, &y_start, &o_len)) {
            continue;
        }
        
        fill_subregion(dumy->overlap_region, y_start, o_len, overlap_list->list[currentID].y_pos_strand, 
        R_INF, overlap_list->list[currentID].y_id, extra_begin, extra_end);        

        x_string = r_string + x_start;
        y_string = dumy->overlap_region;

        end_site = Reserve_Banded_BPM(y_string, Window_Len, x_string, x_len, threshold, &error);

        if (error!=(unsigned int)-1) {
            overlap_list->list[currentID].align_length += x_len;
            append_window_list(&overlap_list->list[currentID], x_start, x_end, y_start, y_start + end_site, (int)error,
            extra_begin, extra_end, threshold, WINDOW, NULL);
        }
    }
}


void verify_ul_window(long long window_start, long long window_end, overlap_region_alloc* overlap_list,Correct_dumy* dumy, const ul_idx_t *uref,
char* r_string, double max_ov_diff_ec, long long blockLen, long long max_error, void *km)
{
    long long i;
    long long currentID;
    long long x_start, y_start, o_len;
    long long Window_Len = blockLen + (max_error << 1);
    char* x_string = NULL;
    char* y_string = NULL;
    long long x_end, x_len;
    int end_site;
    unsigned int error;
    int groupLen = 0;
    int return_sites[GROUP_SIZE];
    unsigned int return_sites_error[GROUP_SIZE];
    uint64_t overlapID[GROUP_SIZE];
    uint64_t y_startGroup[GROUP_SIZE];
    int y_extra_begin[GROUP_SIZE];
    int y_extra_end[GROUP_SIZE];
    int error_threshold[GROUP_SIZE];
    int extra_begin;
    int extra_end;

    ///here are overlaps fully covered by blockLen
    for (i = 0; i < (long long)dumy->length; i++)
    {
        extra_begin = extra_end = 0;
        ///if the window has been fully covered, the interval at x is [window_start, window_end]
        x_len = blockLen;
        currentID = dumy->overlapID[i];
        x_start = window_start;
        ///offset of y
        y_start = (x_start - overlap_list->list[currentID].x_pos_s) + overlap_list->list[currentID].y_pos_s;
        /****************************may have bugs********************************/
        y_start += y_start_offset(x_start, &(overlap_list->list[currentID].f_cigar));
        /****************************may have bugs********************************/
        
  
        if(!determine_overlap_region(max_error, y_start, overlap_list->list[currentID].y_id, Window_Len, uref->ug->u.a[overlap_list->list[currentID].y_id].len,
        &extra_begin, &extra_end, &y_start, &o_len)) {
            continue;
        }

        // if(overlap_list->list[currentID].y_id == 4) {
        //     fprintf(stderr, "[M::%s] q_s::%lld, t_s::%lld, t_pri_l::%lld, aux_beg::%d, aux_end::%d, aln_l::%lld\n", __func__, 
        //         x_start, y_start, o_len, extra_begin, extra_end, Window_Len);
        // }

        fill_subregion_ul(dumy->overlap_region_group[groupLen], y_start, o_len, overlap_list->list[currentID].y_pos_strand, 
        uref, overlap_list->list[currentID].y_id, extra_begin, extra_end);

        y_extra_begin[groupLen] = extra_begin;
        y_extra_end[groupLen] = extra_end;
        overlapID[groupLen] = currentID;
        y_startGroup[groupLen] = y_start;
        error_threshold[groupLen] = max_error;
        x_string = r_string + x_start;
        groupLen++;
        

        if (groupLen == GROUP_SIZE)
        {
            // Reserve_Banded_BPM_4_SSE_only(dumy->overlap_region_group[0], dumy->overlap_region_group[1], 
            // dumy->overlap_region_group[2], dumy->overlap_region_group[3], Window_Len, x_string, blockLen,
            //     return_sites, return_sites_error, max_error, dumy->Peq_SSE);
            return_sites[0] = 
                    Reserve_Banded_BPM(dumy->overlap_region_group[0], Window_Len, x_string, blockLen, max_error, &return_sites_error[0]);
            return_sites[1] = 
                    Reserve_Banded_BPM(dumy->overlap_region_group[1], Window_Len, x_string, blockLen, max_error, &return_sites_error[1]);
            return_sites[2] = 
                    Reserve_Banded_BPM(dumy->overlap_region_group[2], Window_Len, x_string, blockLen, max_error, &return_sites_error[2]);
            return_sites[3] = 
                    Reserve_Banded_BPM(dumy->overlap_region_group[3], Window_Len, x_string, blockLen, max_error, &return_sites_error[3]);
            




            groupLen = 0;
            
            if (return_sites_error[0]!=(unsigned int)-1) {
                overlap_list->list[overlapID[0]].align_length += x_len;

                append_window_list(&overlap_list->list[overlapID[0]], window_start, window_end, 
                            y_startGroup[0], y_startGroup[0] + return_sites[0], (int)return_sites_error[0],
                            y_extra_begin[0], y_extra_end[0], error_threshold[0], blockLen, km);
            }
            

            if (return_sites_error[1]!=(unsigned int)-1) {
                overlap_list->list[overlapID[1]].align_length += x_len;
                append_window_list(&overlap_list->list[overlapID[1]], window_start, window_end, 
                            y_startGroup[1], y_startGroup[1] + return_sites[1], (int)return_sites_error[1],
                            y_extra_begin[1], y_extra_end[1], error_threshold[1], blockLen, km);
            }
            

            if (return_sites_error[2]!=(unsigned int)-1) {
                overlap_list->list[overlapID[2]].align_length += x_len;

                append_window_list(&overlap_list->list[overlapID[2]], window_start, window_end, 
                            y_startGroup[2], y_startGroup[2] + return_sites[2], (int)return_sites_error[2],
                            y_extra_begin[2], y_extra_end[2], error_threshold[2], blockLen, km);
            }
            

            if (return_sites_error[3]!=(unsigned int)-1) {
                overlap_list->list[overlapID[3]].align_length += x_len;

                append_window_list(&overlap_list->list[overlapID[3]], window_start, window_end, 
                            y_startGroup[3], y_startGroup[3] + return_sites[3], (int)return_sites_error[3],
                            y_extra_begin[3], y_extra_end[3], error_threshold[3], blockLen, km);
            } 
        }
    }

    if (groupLen == 1)
    {
        end_site = Reserve_Banded_BPM(dumy->overlap_region_group[0], Window_Len, x_string, blockLen, max_error, &error);

        if (error!=(unsigned int)-1) {
            overlap_list->list[overlapID[0]].align_length += x_len;

            append_window_list(&overlap_list->list[overlapID[0]], window_start, window_end, 
                                y_startGroup[0], y_startGroup[0] + end_site, (int)error,
                                y_extra_begin[0], y_extra_end[0], error_threshold[0], blockLen, km);
        }
    }
    else if (groupLen > 1)
    {
        // Reserve_Banded_BPM_4_SSE_only(dumy->overlap_region_group[0], dumy->overlap_region_group[1], 
        //         dumy->overlap_region_group[2], dumy->overlap_region_group[3], Window_Len, x_string, blockLen,
        //             return_sites, return_sites_error, max_error, dumy->Peq_SSE);

        for (i = 0; i < groupLen; i++)
        {
            return_sites[i] = Reserve_Banded_BPM(dumy->overlap_region_group[i], Window_Len, x_string, blockLen, max_error, &return_sites_error[i]);

            if (return_sites_error[i]!=(unsigned int)-1) {
                overlap_list->list[overlapID[i]].align_length += x_len;
                append_window_list(&overlap_list->list[overlapID[i]], window_start, window_end, 
                                y_startGroup[i], y_startGroup[i] + return_sites[i], (int)return_sites_error[i],
                                y_extra_begin[i], y_extra_end[i], error_threshold[i], blockLen, km);
            }
        }

        groupLen = 0;
    }

    long long reverse_i = dumy->size - 1;
    int threshold;
    
    ///here are overlaps partially covered by blockLen
    for (i = 0; i < (long long)dumy->lengthNT; i++)
    {
        extra_begin = extra_end = 0;
        currentID = dumy->overlapID[reverse_i--];
        x_start = MAX(window_start, (long long)overlap_list->list[currentID].x_pos_s);
        x_end = MIN(window_end, (long long)overlap_list->list[currentID].x_pos_e);

        ///overlap length between [window_start, window_end]
        x_len = x_end - x_start + 1;
        threshold = x_len * max_ov_diff_ec;
        /****************************may have bugs********************************/
        threshold = Adjust_Threshold(threshold, x_len);
        if(threshold > THRESHOLD_MAX_SIZE) threshold = THRESHOLD_MAX_SIZE;
        /****************************may have bugs********************************/
        
        ///offset of y
        y_start = (x_start - overlap_list->list[currentID].x_pos_s) + overlap_list->list[currentID].y_pos_s;
        /****************************may have bugs********************************/
        y_start += y_start_offset(x_start, &(overlap_list->list[currentID].f_cigar));
        /****************************may have bugs********************************/

        Window_Len = x_len + (threshold << 1);

        if(!determine_overlap_region(threshold, y_start, overlap_list->list[currentID].y_id, Window_Len, uref->ug->u.a[overlap_list->list[currentID].y_id].len,
        &extra_begin, &extra_end, &y_start, &o_len)) {
            continue;
        }

        // if(overlap_list->list[currentID].y_id == 4) {
        //     fprintf(stderr, "[M::%s] q_s::%lld, t_s::%lld, t_pri_l::%lld, aux_beg::%d, aux_end::%d, aln_l::%lld\n", __func__, 
        //         x_start, y_start, o_len, extra_begin, extra_end, Window_Len);
        // }
        
        fill_subregion_ul(dumy->overlap_region, y_start, o_len, overlap_list->list[currentID].y_pos_strand, 
        uref, overlap_list->list[currentID].y_id, extra_begin, extra_end);        

        x_string = r_string + x_start;
        y_string = dumy->overlap_region;

        end_site = Reserve_Banded_BPM(y_string, Window_Len, x_string, x_len, threshold, &error);

        if (error!=(unsigned int)-1) {
            overlap_list->list[currentID].align_length += x_len;
            append_window_list(&overlap_list->list[currentID], x_start, x_end, y_start, y_start + end_site, (int)error,
            extra_begin, extra_end, threshold, blockLen, km);
        }
    }
}

int32_t init_waln(int64_t err, int64_t s, int64_t l, int64_t w_l, int64_t* aux_beg, int64_t* aux_end, int64_t* r_s, int64_t* r_l)
{
    (*aux_beg) = (*aux_end) = (*r_s) = (*r_l) = -1;
    ///since w_l == x_len + (err << 1)
    if((s < 0) || (s >= l) || ((l-s+(2*err)+THRESHOLD_MAX_SIZE) < w_l)) return 0;
    (*aux_beg) = (*aux_end) = 0;
    ///s might be less than 0
    (*r_s) = s - err; 
    (*r_l) = l-(*r_s); if((*r_l) > w_l) (*r_l) = w_l;
    (*aux_end) = w_l - (*r_l);

    if ((*r_s) < 0) {
        (*aux_beg) = -(*r_s); (*r_s) = 0; (*r_l) -= (*aux_beg);
    }
    return 1;
}

///[s, e)
int64_t get_num_wins(int64_t s, int64_t e, int64_t block_s)
{
    int64_t nl = e - ((s/block_s)*block_s), nw;
    nw = (nl/block_s); if((nl%block_s)>0) nw++;
    return nw;
}

void gen_str_seq(char *dst, int64_t s, int64_t pri_l, uint8_t rev, const ul_idx_t *uref, long long id, int64_t aux_beg, int64_t aux_end)
{
    // int64_t l = pri_l + aux_beg + aux_end;    
    memset(dst, 'N', aux_beg);
    retrieve_u_seq(NULL, dst+aux_beg, &(uref->ug->u.a[id]), rev, s, pri_l, NULL);
    memset(dst+aux_beg+pri_l, 'N', aux_end);
}

void verify_ul_window_s(overlap_region *z, const ul_idx_t *uref, char* qstr, char *tstr,
double e_rate, int64_t w_l, int64_t e_max, void *km)
{
    int64_t q_s, q_e, nw, k, q_l;
    int64_t aux_beg, aux_end, t_s, thre, aln_l, t_pri_l, t_end;
    char *q_string, *t_string; unsigned int error;
    z->w_list.n = 0; nw = get_num_wins(z->x_pos_s, z->x_pos_e+1, w_l);
    get_win_se_by_normalize_xs(z, (z->x_pos_s/w_l)*w_l, w_l, &q_s, &q_e);
    // q_s = z->x_pos_s; get_win_id_by_s(z, q_s, w_l, &q_e); 
    for (k = 0; k < nw; k++) {
        aux_beg = aux_end = 0; q_l = 1 + q_e - q_s;
        thre = q_l*e_rate; thre = Adjust_Threshold(thre, q_l);
        if(thre > THRESHOLD_MAX_SIZE) thre = THRESHOLD_MAX_SIZE;
        ///offset of y
        t_s = (q_s - z->x_pos_s) + z->y_pos_s;
        t_s += y_start_offset(q_s, &(z->f_cigar));

        aln_l = q_l + (thre<<1);
        // if(z->y_id == 115) {
        //     fprintf(stderr, "+[M::] q_s::%ld, t_s::%ld, t_pri_l::%ld, aux_beg::%ld, aux_end::%ld, aln_l::%ld, t_end::%ld, error::%u\n", 
        //     q_s, t_s, t_pri_l, aux_beg, aux_end, aln_l, t_end, error);
        // }
        if(init_waln(thre, t_s, uref->ug->u.a[z->y_id].len, aln_l, &aux_beg, &aux_end, &t_s, &t_pri_l)) {
            gen_str_seq(tstr, t_s, t_pri_l, z->y_pos_strand, uref, z->y_id, aux_beg, aux_end);
            q_string = qstr+q_s; t_string = tstr;
            t_end = Reserve_Banded_BPM(t_string, aln_l, q_string, q_l, thre, &error);
            // if(z->y_id == 115) {
            //     fprintf(stderr, "-[M::] q_s::%ld, t_s::%ld, t_pri_l::%ld, aux_beg::%ld, aux_end::%ld, aln_l::%ld, t_end::%ld, error::%u, thre::%ld\n", 
            //     q_s, t_s, t_pri_l, aux_beg, aux_end, aln_l, t_end, error, thre);
            // }
            if (error!=((unsigned int)-1)) {
                z->align_length += q_l;
                ///t_s do not have aux_beg, while t_s + t_end (aka, te) has 
                append_window_list(z, q_s, q_e, t_s, t_s + t_end, error, aux_beg, aux_end, thre, w_l, km);
            }
        }
        q_s = q_e + 1; q_e = q_s + w_l - 1; 
        if(q_e >= (int64_t)z->x_pos_e) q_e = z->x_pos_e;
    }

    // if(q_e != (int64_t)z->x_pos_e) {
    //     fprintf(stderr, "[M::%s] q_e::%ld, z->x_pos_s::%u, z->x_pos_e::%u, w_l::%ld, nw::%ld\n", __func__, 
    //     q_e, z->x_pos_s, z->x_pos_e, w_l, nw);
    // }
    assert(q_e == (int64_t)z->x_pos_e);
}

///error_rate should be 30%
long long get_high_error(long long x_start, long long x_end, 
long long y_start, long long y_end, long long y_id, long long y_strand, long long pre_threshold,
long long n_steps, float error_rate, All_reads* R_INF, Correct_dumy* dumy,
UC_Read* g_read)
{
    long long stepLen = (x_end - x_start + 1) / n_steps;
    if((x_end - x_start + 1) % n_steps != 0)
    {
        stepLen++;
    }
    long long SubLen, SubWindowLen;
    long long SubThreshold = THRESHOLD_MAX_SIZE;
    int extra_begin, extra_end;
    long long o_len;
    long long T_error = 0;
    y_start = y_start + pre_threshold;

    while (x_start <= x_end)
    {
        SubLen = x_end - x_start + 1;
        if(SubLen > stepLen)
        {
            SubLen = stepLen;
        }


        SubThreshold = SubLen * error_rate;
        if(SubThreshold > THRESHOLD_MAX_SIZE)
        {
            SubThreshold = THRESHOLD_MAX_SIZE;
        }

        SubThreshold = Adjust_Threshold(SubThreshold, SubLen);


        SubWindowLen = SubLen + (SubThreshold << 1);
        if(determine_overlap_region(SubThreshold, y_start, y_id, SubWindowLen, Get_READ_LENGTH((*R_INF), y_id),
        &extra_begin, &extra_end, &y_start, &o_len) == 0)
        {
            T_error = T_error + (x_end - x_start + 1) * error_rate * 1.5;
            break;
        }

        fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, R_INF, y_id, 
        extra_begin, extra_end);

        char* x_string = g_read->seq + x_start;
        char* y_string = dumy->overlap_region;
        int end_site;
        unsigned int error;

        end_site = Reserve_Banded_BPM(y_string, SubWindowLen, x_string, SubLen, SubThreshold, &error);
                    
        ///if error = -1, unmatched
        if (error!=(unsigned int)-1)
        {
            T_error = T_error + error;
            y_start = y_start + end_site - extra_begin + 1;
        }
        else
        {
            T_error = T_error + SubLen * error_rate * 1.5;
            y_start = y_start + SubThreshold - extra_begin + SubLen;
        }
        
        x_start = x_start + SubLen;
    }


    return T_error;
}

inline int double_error_threshold(int pre_threshold, int x_len)
{

    pre_threshold = Adjust_Threshold(pre_threshold, x_len);
    int threshold = pre_threshold * 2;
    ///may have some bugs
    if(x_len >= 300 && threshold < THRESHOLD_MAX_SIZE)
    {
        threshold = THRESHOLD_MAX_SIZE;
    }
    
    if(threshold > THRESHOLD_MAX_SIZE)
    {
        threshold = THRESHOLD_MAX_SIZE;
    }

    return threshold;
}

inline int double_ul_error_threshold(int pre_threshold, int x_len)
{

    pre_threshold = Adjust_Threshold(pre_threshold, x_len);
    int threshold = THRESHOLD_UL_MAX * x_len;
    if(threshold < pre_threshold) threshold = pre_threshold;
    if(threshold > THRESHOLD_MAX_SIZE) threshold = THRESHOLD_MAX_SIZE;
    return threshold;
}


inline int verify_sub_window(All_reads* R_INF, Correct_dumy* dumy, UC_Read* g_read, 
long long x_beg, long long xLen, long long y_beg, long long yLen, uint64_t y_id, 
uint64_t y_pos_strand, int threshold, int alignment_strand,
unsigned int* get_error, int* get_y_end, int* get_x_end, int* get_aligned_xLen)
{
    (*get_aligned_xLen) = 0;
    (*get_y_end) = -1;
    (*get_x_end) = -1;
    (*get_error) = (unsigned int)-1;

    int extra_begin, extra_end, r_x_end, r_y_end, aligned_xLen;
    long long o_len;
    unsigned int r_error;
    if(!determine_overlap_region(threshold, y_beg, y_id, yLen, 
    Get_READ_LENGTH((*R_INF), y_id), &extra_begin, &extra_end, &y_beg, &o_len))
    {
        return 0;
    }

    fill_subregion(dumy->overlap_region, y_beg, o_len, y_pos_strand, R_INF, 
    y_id, extra_begin, extra_end);

    char* x_string = g_read->seq + x_beg;
    char* y_string = dumy->overlap_region;

    aligned_xLen = 0;

    alignment_extension(y_string, yLen, x_string, xLen, threshold, 
                    alignment_strand, &r_error, &r_y_end, &r_x_end, &aligned_xLen);
    
    (*get_error) = r_error;
    (*get_y_end) = r_y_end;
    (*get_x_end) = r_x_end;
    (*get_aligned_xLen) = aligned_xLen;

    if(aligned_xLen == 0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

inline int verify_ul_sub_window(const ul_idx_t *uref, Correct_dumy* dumy, UC_Read* g_read, 
long long x_beg, long long xLen, long long y_beg, long long yLen, uint64_t y_id, 
uint64_t y_pos_strand, int threshold, int alignment_strand,
unsigned int* get_error, int* get_y_end, int* get_x_end, int* get_aligned_xLen)
{
    (*get_aligned_xLen) = 0;
    (*get_y_end) = -1;
    (*get_x_end) = -1;
    (*get_error) = (unsigned int)-1;

    int extra_begin, extra_end, r_x_end, r_y_end, aligned_xLen;
    long long o_len;
    unsigned int r_error;
    if(!determine_overlap_region(threshold, y_beg, y_id, yLen, 
    uref->ug->u.a[y_id].len, &extra_begin, &extra_end, &y_beg, &o_len))
    {
        return 0;
    }

    fill_subregion_ul(dumy->overlap_region, y_beg, o_len, y_pos_strand, uref, 
    y_id, extra_begin, extra_end);

    // if(y_id == 6) {
    //     fprintf(stderr, "-[M::%s::aln_dir->%d] qs->%lld, ts->%lld, thres->%d, aux_beg->%d, aux_end->%d, t_pri_l->%lld\n", 
    //     __func__, alignment_strand, x_beg, y_beg, threshold, extra_begin, extra_end, o_len);
    // }

    char* x_string = g_read->seq + x_beg;
    char* y_string = dumy->overlap_region;

    aligned_xLen = 0;

    alignment_extension(y_string, yLen, x_string, xLen, threshold, 
                    alignment_strand, &r_error, &r_y_end, &r_x_end, &aligned_xLen);
    
    (*get_error) = r_error;
    (*get_y_end) = r_y_end;
    (*get_x_end) = r_x_end;
    (*get_aligned_xLen) = aligned_xLen;

    if(aligned_xLen == 0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

inline int64_t get_init_err_thres(int64_t len, double e_rate, int64_t block_s, int64_t block_err)
{
    if(len >= block_s) return block_err;
    int64_t thres = len * e_rate;
    thres = Adjust_Threshold(thres, len);
    if(thres > THRESHOLD_MAX_SIZE) thres = THRESHOLD_MAX_SIZE;
    return thres;
}


uint32_t get_init_paras(All_reads* rref, const ul_idx_t *uref, overlap_region *z, int64_t x_s, int64_t x_e, double e_rate, int64_t block_s, 
                        int64_t *r_ys, int64_t *r_ex_beg, int64_t *r_ex_end, int64_t *r_err_thre)
{
    int e, ex_beg, ex_end; long long y_s, o_len, Window_Len;
    e = get_init_err_thres(x_e+1-x_s, e_rate, block_s, rref?THRESHOLD:THRESHOLD_MAX_SIZE);
    y_s = (x_s-z->x_pos_s) + z->y_pos_s; y_s += y_start_offset(x_s, &(z->f_cigar));
    Window_Len = (x_e+1-x_s) + (e<<1);

    if(!determine_overlap_region(e, y_s, z->y_id, Window_Len, (rref?(Get_READ_LENGTH((*rref), z->y_id)):(uref->ug->u.a[z->y_id].len)),
    &ex_beg, &ex_end, &y_s, &o_len)) {
        return 0;
    }
    (*r_ys) = y_s; (*r_ex_beg) = ex_beg; (*r_ex_end) = ex_end; (*r_err_thre) = e;
    return 1;
}

int64_t check_coverage_gap(uint64_t *v_idx, uint64_t w_s, uint64_t w_e, int64_t block_s)
{
    int64_t wid = w_s/block_s, a_n = (uint32_t)(v_idx[wid]), k;
    uint64_t *a = v_idx + (v_idx[wid]>>32);
    for (k = 0; k < a_n; k++) {
        if(((a[k]>>32) == w_s) && (((uint32_t)(a[k])) == w_e)) return 1;
    }
    return 0;
}

inline double non_trim_error_rate(overlap_region *z, All_reads* rref, const ul_idx_t *uref, const kvec_t_u64_warp* v_idx, Correct_dumy* dumy, UC_Read* g_read, double e_rate, int64_t block_s)
{
    int64_t nw, aw = z->w_list.n, k, m, w_id, wn_id, w_s, w_e, idx_e, tErr = 0, tLen = 0, y_s, ex_beg, ex_end, err_thre, p_err_thre;
    int64_t x_len, Window_Len, y_beg_left, y_beg_right;
    unsigned int r_error_left, r_error_right; int32_t r_x_end_left, r_y_end_left, aligned_xLen_left, r_x_end_right, r_y_end_right, aligned_xLen_right;
    nw = get_num_wins(z->x_pos_s, z->x_pos_e+1, block_s);
    assert(nw >= aw && aw > 0);

    for (k = aw-1, idx_e = nw; k >= 0; k--) {
        w_id = get_win_id_by_e(z, z->w_list.a[k].x_end, block_s, &w_s);
        assert(w_s == z->w_list.a[k].x_start && w_id < idx_e && k <= w_id);
        tLen += z->w_list.a[k].x_end + 1 - z->w_list.a[k].x_start;
        tErr += z->w_list.a[k].error;///matched window
        // if(z->y_id == 1) {
        //     fprintf(stderr, "+[M::%s] ws->%d, we->%d, tot_l->%ld, tot_e->%ld\n", 
        //                         __func__, z->w_list.a[k].x_start, z->w_list.a[k].x_end, tLen, tErr);
        // }
        // if(k != w_id) z->w_list.a[w_id] = z->w_list.a[k];
        ///from mapped window w_list.a[k] to the following unmapped windows
        for (m = w_id+1, w_e = z->w_list.a[k].x_end; m < idx_e; m++) {
            w_s = w_e + 1;
            wn_id = get_win_id_by_s(z, w_s, block_s, &w_e);
            assert(wn_id == m); x_len = w_e + 1 - w_s; tLen += x_len;
            ///check if there are some windows that cannot be algined by any overlaps/unitigs
            ///if no, it is likely that the UL read itself has issues
            if(uref && v_idx && z->is_match == 4) {
                if(check_coverage_gap(v_idx->a.a, w_s, w_e, block_s)) {
                    tErr += THRESHOLD_MAX_SIZE;
                    // if(z->y_id == 1) {
                    //     fprintf(stderr, "-[M::%s] ws->%ld, we->%ld, tot_l->%ld, tot_e->%ld\n", __func__, w_s, w_e, tLen, tErr);
                    // }
                    continue;
                } 
            }
            if(!get_init_paras(rref, uref, z, w_s, w_e, e_rate, block_s, &y_s, &ex_beg, &ex_end, &err_thre)) {
                tErr += x_len;
                // if(z->y_id == 1) {
                //     fprintf(stderr, "-[M::%s] ws->%ld, we->%ld, tot_l->%ld, tot_e->%ld\n", __func__, w_s, w_e, tLen, tErr);
                // }
                continue;
            }
            p_err_thre = err_thre;

            if(rref) {
                err_thre = double_error_threshold(err_thre, x_len);
            } else {
                err_thre = double_ul_error_threshold(err_thre, x_len);
            }
            Window_Len = x_len + (err_thre << 1); 
            r_error_left = r_error_right = 0;
            aligned_xLen_left = aligned_xLen_right = 0;
            y_beg_left = y_beg_right = -1;

            if(m == w_id+1) { ///if the previous window is mapped
                y_beg_left = z->w_list.a[k].y_end + 1;///incorrect
            }

            if(m+1 == idx_e && k+1 < aw) { ///if the next window is mapped
                y_beg_right = z->w_list.a[k+1].y_start-x_len;///incorrect
            }

            if(y_beg_left == -1 && y_beg_right == -1) {
                y_beg_left = y_s;
                if(ex_beg >= 0) y_beg_left = y_beg_left + p_err_thre - ex_beg;
                y_beg_right = y_beg_left;
            }

            if(y_beg_left == -1 && y_beg_right != -1) y_beg_left = y_beg_right;
            if(y_beg_right == -1 && y_beg_left != -1) y_beg_right = y_beg_left;

            if(y_beg_left != -1) {///note: this function will change tstr/qstr
                if(rref) {
                    verify_sub_window(rref, dumy, g_read, w_s, x_len, y_beg_left, Window_Len, 
                    z->y_id, z->y_pos_strand, err_thre, 0, &r_error_left, &r_y_end_left, &r_x_end_left, &aligned_xLen_left);
                } else {
                    verify_ul_sub_window(uref, dumy, g_read, w_s, x_len, y_beg_left, Window_Len, 
                    z->y_id, z->y_pos_strand, err_thre, 0, &r_error_left, &r_y_end_left, &r_x_end_left, &aligned_xLen_left);
                }
            }

            if(y_beg_right != -1) {
                if(rref) {
                    verify_sub_window(rref, dumy, g_read, w_s, x_len, y_beg_right, Window_Len, z->y_id, z->y_pos_strand,
                    err_thre, 1, &r_error_right, &r_y_end_right, &r_x_end_right, &aligned_xLen_right);
                } else {
                    verify_ul_sub_window(uref, dumy, g_read, w_s, x_len, y_beg_right, Window_Len, z->y_id, z->y_pos_strand,
                    err_thre, 1, &r_error_right, &r_y_end_right, &r_x_end_right, &aligned_xLen_right);
                }
            }

            ///aligned in both directions
            if(aligned_xLen_left != 0 && aligned_xLen_right != 0) {
                if(aligned_xLen_left + aligned_xLen_right <= x_len) {
                    tErr += r_error_left + r_error_right + (x_len - aligned_xLen_left - aligned_xLen_right);
                } else {
                    float E_rate = (float)(x_len)/(float)(aligned_xLen_left + aligned_xLen_right);
                    tErr += (r_error_left + r_error_right)*E_rate;
                }
            }///not aligned in both directions
            else if(aligned_xLen_left == 0 && aligned_xLen_right == 0) {
                tErr += x_len;
            }///only aligned in left
            else if(aligned_xLen_left != 0) {
                tErr += r_error_left + (x_len - aligned_xLen_left);
            }///only aligned in right
            else if(aligned_xLen_right != 0) {
                tErr += r_error_right + (x_len - aligned_xLen_right);
            }
            // if(z->y_id == 1) {
            //     fprintf(stderr, "*[M::%s] qs->%ld, ts->%ld, tb[0]->%ld, tb[1]->%ld, di[0]->%u, di[1]->%u, al[0]->%d, al[1]->%d, err_thre->%ld\n", __func__, 
            //     w_s, y_s, y_beg_left, y_beg_right, r_error_left, r_error_right, aligned_xLen_left, aligned_xLen_right, err_thre);
            // }
            // if(z->y_id == 1) {
            //     fprintf(stderr, "-[M::%s] ws->%ld, we->%ld, tot_l->%ld, tot_e->%ld\n", __func__, w_s, w_e, tLen, tErr);
            // }
        }
        idx_e = w_id;
    }

    if(idx_e > 0) {
        for (m = 0, w_e = (int64_t)z->x_pos_s-1; m < idx_e; m++) {
            w_s = w_e + 1;
            wn_id = get_win_id_by_s(z, w_s, block_s, &w_e);
            assert(wn_id == m); x_len = w_e + 1 - w_s; tLen += x_len;
            ///check if there are some windows that cannot be algined by any overlaps/unitigs
            ///if no, it is likely that the UL read itself has issues
            if(uref && v_idx && z->is_match == 4) {
                if(check_coverage_gap(v_idx->a.a, w_s, w_e, block_s)) {
                    tErr += THRESHOLD_MAX_SIZE;
                    // if(z->y_id == 1) {
                    //     fprintf(stderr, "-[M::%s] ws->%ld, we->%ld, tot_l->%ld, tot_e->%ld\n", __func__, w_s, w_e, tLen, tErr);
                    // }
                    continue;
                } 
                // else {
                //     if(z->y_id == 575) {
                //         fprintf(stderr, "---[M::%s::] z::y_id->%u, w_s->%ld, w_e->%ld\n", __func__, z->y_id, w_s, w_e);
                //     } 
                // }
            }
            if(!get_init_paras(rref, uref, z, w_s, w_e, e_rate, block_s, &y_s, &ex_beg, &ex_end, &err_thre)) {
                tErr += x_len;
                // if(z->y_id == 1) {
                //     fprintf(stderr, "-[M::%s] ws->%ld, we->%ld, tot_l->%ld, tot_e->%ld\n", __func__, w_s, w_e, tLen, tErr);
                // }
                continue;
            }
            p_err_thre = err_thre;
            if(rref) {
                err_thre = double_error_threshold(err_thre, x_len);
            } else {
                err_thre = double_ul_error_threshold(err_thre, x_len);
            }
            Window_Len = x_len + (err_thre << 1); 
            r_error_left = r_error_right = 0;
            aligned_xLen_left = aligned_xLen_right = 0;
            y_beg_left = y_beg_right = -1;

            ///impossible that the previous window is mapped
            // if(m == w_id+1) { ///if the previous window is mapped
            //     y_beg_left = z->w_list.a[k].y_end + 1;
            // }

            if(m+1 == idx_e && k+1 < aw) { ///if the next window is mapped
                y_beg_right = z->w_list.a[k+1].y_start-x_len;
            }

            if(y_beg_left == -1 && y_beg_right == -1) {
                y_beg_left = y_s;
                if(ex_beg >= 0) y_beg_left = y_beg_left + p_err_thre - ex_beg;
                y_beg_right = y_beg_left;
            }

            if(y_beg_left == -1 && y_beg_right != -1) y_beg_left = y_beg_right;
            if(y_beg_right == -1 && y_beg_left != -1) y_beg_right = y_beg_left;

            if(y_beg_left != -1) {
                if(rref) {
                    verify_sub_window(rref, dumy, g_read, w_s, x_len, y_beg_left, Window_Len, 
                    z->y_id, z->y_pos_strand, err_thre, 0, &r_error_left, &r_y_end_left, &r_x_end_left, &aligned_xLen_left);
                } else {
                    verify_ul_sub_window(uref, dumy, g_read, w_s, x_len, y_beg_left, Window_Len, 
                    z->y_id, z->y_pos_strand, err_thre, 0, &r_error_left, &r_y_end_left, &r_x_end_left, &aligned_xLen_left);
                }
            }

            if(y_beg_right != -1) {
                if(rref) {
                    verify_sub_window(rref, dumy, g_read, w_s, x_len, y_beg_right, Window_Len, z->y_id, z->y_pos_strand,
                    err_thre, 1, &r_error_right, &r_y_end_right, &r_x_end_right, &aligned_xLen_right);
                } else {
                    verify_ul_sub_window(uref, dumy, g_read, w_s, x_len, y_beg_right, Window_Len, z->y_id, z->y_pos_strand,
                    err_thre, 1, &r_error_right, &r_y_end_right, &r_x_end_right, &aligned_xLen_right);
                }
            }

            ///aligned in both directions
            if(aligned_xLen_left != 0 && aligned_xLen_right != 0) {
                if(aligned_xLen_left + aligned_xLen_right <= x_len) {
                    tErr += r_error_left + r_error_right + (x_len - aligned_xLen_left - aligned_xLen_right);
                } else {
                    float E_rate = (float)(x_len)/(float)(aligned_xLen_left + aligned_xLen_right);
                    tErr += (r_error_left + r_error_right)*E_rate;
                }
            }///not aligned in both directions
            else if(aligned_xLen_left == 0 && aligned_xLen_right == 0) {
                tErr += x_len;
            }///only aligned in left
            else if(aligned_xLen_left != 0) {
                tErr += r_error_left + (x_len - aligned_xLen_left);
            }///only aligned in right
            else if(aligned_xLen_right != 0) {
                tErr += r_error_right + (x_len - aligned_xLen_right);
            }
            // if(z->y_id == 1) {
            //     fprintf(stderr, "-[M::%s] ws->%ld, we->%ld, tot_l->%ld, tot_e->%ld\n", __func__, w_s, w_e, tLen, tErr);
            // }
        }
    }
    
    assert(tLen == z->x_pos_e + 1 - z->x_pos_s);
    return (double)(tErr)/(double)(tLen);
}

void append_unmatched_wins(overlap_region *z, int64_t block_s)
{
    int64_t nw, aw = z->w_list.n, k, m, w_id, wn_id, w_s, w_e, idx_e;
    nw = get_num_wins(z->x_pos_s, z->x_pos_e+1, block_s);
    assert(nw >= aw && aw > 0);
    if(nw == aw) return;///done
    kv_resize(window_list, z->w_list, (uint64_t)nw); z->w_list.n = nw;
    for (k = aw-1, idx_e = nw; k >= 0; k--) {
        w_id = get_win_id_by_e(z, z->w_list.a[k].x_end, block_s, &w_s);
        assert(w_s == z->w_list.a[k].x_start && w_id < idx_e && k <= w_id);
        if(k != w_id) z->w_list.a[w_id] = z->w_list.a[k];
        for (m = w_id+1, w_e = z->w_list.a[k].x_end; m < idx_e; m++) {
            z->w_list.a[m].cidx = z->w_list.a[m].clen = 0;
            z->w_list.a[m].y_start = z->w_list.a[m].y_end = -1;
            z->w_list.a[m].error = z->w_list.a[m].error_threshold = -1;
            z->w_list.a[m].extra_begin = z->w_list.a[m].extra_end = 0;
            z->w_list.a[m].x_start = w_e + 1;
            wn_id = get_win_id_by_s(z, z->w_list.a[m].x_start, block_s, &w_e);
            z->w_list.a[m].x_end = w_e;
            assert(wn_id == m);
        }
        idx_e = w_id;
    }

    if(idx_e > 0) {
        for (m = 0, w_e = (int64_t)z->x_pos_s-1; m < idx_e; m++) {
            z->w_list.a[m].cidx = z->w_list.a[m].clen = 0;
            z->w_list.a[m].y_start = z->w_list.a[m].y_end = -1;
            z->w_list.a[m].error = z->w_list.a[m].error_threshold = -1;
            z->w_list.a[m].extra_begin = z->w_list.a[m].extra_end = 0;
            z->w_list.a[m].x_start = w_e + 1;
            wn_id = get_win_id_by_s(z, z->w_list.a[m].x_start, block_s, &w_e);
            z->w_list.a[m].x_end = w_e;
            assert(wn_id == m);
        }
    }
}

/**
inline double non_trim_ul_error_rate(overlap_region_alloc* overlap_list, long long ID,
const ul_idx_t *uref, Correct_dumy* dumy, UC_Read* g_read)
{
    long long tLen, tError,i, subWinLen, subWinNum;
    
    tLen = 0;
    tError = 0;

    subWinNum = overlap_list->list[ID].w_list_length;

    
    for (i = 0; i < subWinNum; i++)
    {
        subWinLen = overlap_list->list[ID].w_list[i].x_end - overlap_list->list[ID].w_list[i].x_start + 1;
        tLen += subWinLen;

        if(overlap_list->list[ID].w_list[i].y_end != -1)
        {
            tError += overlap_list->list[ID].w_list[i].error;
        }
        else
        {
            int x_len = subWinLen;
            int threshold = double_ul_error_threshold(overlap_list->list[ID].w_list[i].error_threshold, x_len);
            int Window_Len = x_len + (threshold << 1);
            unsigned int r_error_left = 0;
            int r_x_end_left, r_y_end_left, aligned_xLen_left;
            unsigned int r_error_right = 0;
            int r_x_end_right, r_y_end_right, aligned_xLen_right;
            long long y_beg_left, y_beg_right;
            
            aligned_xLen_left = aligned_xLen_right = 0;
            y_beg_left = y_beg_right = -1;

            if(overlap_list->list[ID].w_list[i].y_start == -1)
            {
                tError += x_len;
                continue;
            }

            ///if the previous window is mapped
            if(i > 0 && overlap_list->list[ID].w_list[i - 1].y_end != -1)
            {
                y_beg_left = overlap_list->list[ID].w_list[i - 1].y_end + 1;
            }

            ///if the next window is mapped
            if(i < (long long)(overlap_list->list[ID].w_list_length - 1) && overlap_list->list[ID].w_list[i + 1].y_end != -1)
            {
                y_beg_right = 1 + overlap_list->list[ID].w_list[i + 1].y_start - 1 - x_len;
            }


            if(y_beg_left == -1 && y_beg_right == -1)
            {
                y_beg_left = overlap_list->list[ID].w_list[i].y_start;
                if(overlap_list->list[ID].w_list[i].extra_begin >= 0)
                {
                    y_beg_left = y_beg_left + overlap_list->list[ID].w_list[i].error_threshold - 
                    overlap_list->list[ID].w_list[i].extra_begin;
                }
                y_beg_right = y_beg_left;
            }

            if(y_beg_left == -1 && y_beg_right != -1)
            {
                y_beg_left = y_beg_right;
            }

            if(y_beg_right == -1 && y_beg_left != -1)
            {
                y_beg_right = y_beg_left;
            }

            
            if(y_beg_left != -1)
            {
                verify_ul_sub_window(uref, dumy, g_read, overlap_list->list[ID].w_list[i].x_start, 
                x_len, y_beg_left, Window_Len, overlap_list->list[ID].y_id, overlap_list->list[ID].y_pos_strand, 
                threshold, 0, &r_error_left, &r_y_end_left, &r_x_end_left, &aligned_xLen_left);
            }

            if(y_beg_right != -1)
            {
                verify_ul_sub_window(uref, dumy, g_read, overlap_list->list[ID].w_list[i].x_start, 
                x_len, y_beg_right, Window_Len, overlap_list->list[ID].y_id, overlap_list->list[ID].y_pos_strand,
                threshold, 1, &r_error_right, &r_y_end_right, &r_x_end_right, &aligned_xLen_right);
            }

            ///aligned in both direction
            if(aligned_xLen_left != 0 && aligned_xLen_right != 0)
            {
                if(aligned_xLen_left + aligned_xLen_right <= x_len)
                {
                    tError = tError + r_error_left + r_error_right + 
                    (x_len - aligned_xLen_left - aligned_xLen_right);
                }
                else
                {
                    float E_rate = (float)(x_len)/(float)(aligned_xLen_left + aligned_xLen_right);
                    tError = tError + (r_error_left + r_error_right)*E_rate;
                }
            }///not aligned in both direction
            else if(aligned_xLen_left == 0 && aligned_xLen_right == 0)
            {
                tError += x_len;
            }///only aligned in left
            else if(aligned_xLen_left != 0)
            {
                tError = tError + r_error_left + (x_len - aligned_xLen_left);
            }///only aligned in right
            else if(aligned_xLen_right != 0)
            {
                tError = tError + r_error_right + (x_len - aligned_xLen_right);
            }
        }
    }

    double error_rate = (double)(tError)/(double)(tLen);

    return error_rate;
}

int calculate_hpm_errors(char* x, int x_len, char* y, int y_len, CIGAR* cigar, int error)
{
    int x_i, y_i, cigar_i;
    x_i = 0;
    y_i = 0;
    cigar_i = 0;
    int operation;
    int operationLen;
    int i;
    int cigar_error = 0;
    int hpm_error = 0;


    while (cigar_i < cigar->length)
    {
        operation = cigar->C_C[cigar_i];
        operationLen = cigar->C_L[cigar_i];

        if (operation == 0)
        {
            x_i = x_i + operationLen;
            y_i = y_i + operationLen;
        }
        else if (operation == 1)
        {
            cigar_error += operationLen;
            for (i = 0; i < operationLen; i++)
            {
                if(if_is_homopolymer_repeat(x_i, x, x_len) || if_is_homopolymer_repeat(y_i, y, y_len))
                {
                    hpm_error++;
                }

                x_i++;
                y_i++;
            }
        }
        else if (operation == 2)
        {

            if(if_is_homopolymer_repeat(x_i, x, x_len) || if_is_homopolymer_repeat(y_i, y, y_len))
            {
                hpm_error++;
            }
            cigar_error += operationLen;
            y_i += operationLen;
        }
        else if (operation == 3)
        {

            if(if_is_homopolymer_repeat(x_i, x, x_len) || if_is_homopolymer_repeat(y_i, y, y_len))
            {
                hpm_error++;
            }

            cigar_error += operationLen;
            x_i += operationLen;
        }
        
        cigar_i++;
    }
    return hpm_error;
}

int verify_cigar(char* x, int x_len, char* y, int y_len, CIGAR* cigar, int error)
{
    int x_i, y_i, cigar_i;
    x_i = 0;
    y_i = 0;
    cigar_i = 0;
    int operation;
    int operationLen;
    int i;
    int cigar_error = 0;
    int flag_error = 0;

    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2 means there are more y, 3 means there are more x
    while (cigar_i < cigar->length)
    {
        operation = cigar->C_C[cigar_i];
        operationLen = cigar->C_L[cigar_i];

        if (operation == 0)
        {
            for (i = 0; i < operationLen; i++)
            {

                if (x[x_i]!=y[y_i])
                {
                    ///fprintf(stderr, "error match\n");
                    flag_error = 1;
                }
                x_i++;
                y_i++;
            }
        }
        else if (operation == 1)
        {
            cigar_error += operationLen;
            for (i = 0; i < operationLen; i++)
            {

                if (x[x_i]==y[y_i])
                {
                    ///fprintf(stderr, "error mismatch, cigar_i: %d, x_i: %d, y_i: %d\n",cigar_i, x_i, y_i);
                    flag_error = 1;
                }
                x_i++;
                y_i++;
            }
        }
        else if (operation == 2)
        {
            cigar_error += operationLen;
            y_i += operationLen;
        }
        else if (operation == 3)
        {
            cigar_error += operationLen;
            x_i += operationLen;
        }
        
        cigar_i++;
    }

    
    if (cigar_error != error)
    {
        
        // fprintf(stderr, "error cigar_error: cigar_error: %d, error: %d\n", cigar_error, error);
        // for (i = 0; i < cigar->length; i++)
        // {
        //     fprintf(stderr, "%u: %u\n", cigar->C_L[i], cigar->C_C[i]);
        // }
        
        
       flag_error = 1;
        
    }
    
   
    if (flag_error == 1)
    {
        
        // print_string(x, x_len);
        // print_string(y, y_len);
        // fprintf(stderr, "x_len: %d, y_len: %d, cigar_len: %d, error: %d\n", x_len, y_len, cigar->length, error);
        // for (i = 0; i < cigar->length; i++)
        // {
        //     fprintf(stderr, "%u: %u\n", cigar->C_L[i], cigar->C_C[i]);
        // }
        
        
    }
    
    
    return flag_error;
    
}
**/

int32_t scan_cigar(window_list *idx, window_list_alloc *cc, int64_t* get_error, int64_t scanXLen, int64_t direction)
{
    uint8_t c = (uint8_t)-1; uint32_t cl = (uint32_t)-1;
    (*get_error) = -1;
    if(idx->clen == 1) {
        get_cigar_cell(idx, cc, 0, &c, &cl);
        if(c == 0) {
            (*get_error) = 0;
            return 1;
        }
    }
    int32_t x_i = 0, y_i = 0, c_i, c_n = idx->clen, c_err = 0;
    uint32_t i;

    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2: there are more bases at y, 3: there are more bases at x
    if(direction == 0) {
        for (c_i = 0; c_i < c_n; c_i++) {
            get_cigar_cell(idx, cc, c_i, &c, &cl);
            if (c == 0) { //match
                x_i += cl; y_i += cl;
                if(x_i >= scanXLen) {
                    (*get_error) = c_err;
                    return 1;
                }
            }
            else if (c == 1) {
                for (i = 0; i < cl; i++) {
                    x_i++; y_i++; c_err++;
                    if(x_i >= scanXLen) {
                        (*get_error) = c_err;
                        return 1;
                    }
                }
            }
            else if (c == 2) {///y has more bases than x
                c_err += cl; y_i += cl;
            }
            else if (c == 3) {///x has more bases than y
                for (i = 0; i < cl; i++) {
                    x_i++; c_err++;
                    if(x_i >= scanXLen) {
                        (*get_error) = c_err;
                        return 1;
                    }
                }
            }
        }
    } else {
        for (c_i = c_n-1; c_i >= 0; c_i--) {
            get_cigar_cell(idx, cc, c_i, &c, &cl);
            if (c == 0) { //match
                x_i += cl; y_i += cl;
                if(x_i >= scanXLen) {
                    (*get_error) = c_err;
                    return 1;
                }
            } else if (c == 1) { //mismatch
                for (i = 0; i < cl; i++) {
                    x_i++; y_i++; c_err++;
                    if(x_i >= scanXLen) {
                        (*get_error) = c_err;
                        return 1;
                    }
                }
            } else if (c == 2) {///y has more bases than x
                c_err += cl; y_i += cl;
            } else if (c == 3) {///x has more bases than y
                for (i = 0; i < cl; i++) {
                    x_i++; c_err++;
                    if(x_i >= scanXLen) {
                        (*get_error) = c_err;
                        return 1;
                    }
                }
            }
        }
    }

    (*get_error) = c_err;
    return 0;
}

///[scanXbeg, scanXend]
int scan_cigar_interval(window_list *idx, window_list_alloc *cc, int64_t* get_error, int64_t scanXbeg, int64_t scanXend)
{
    uint8_t c; uint32_t cl;
    (*get_error) = -1;
     if(idx->clen == 1) {
        get_cigar_cell(idx, cc, 0, &c, &cl);
        if(c == 0) {
            (*get_error) = 0;
            return 1;
        }
    }



    int32_t x_i = 0, y_i = 0, c_i, c_n = idx->clen, c_err = 0;
    uint32_t i;
    

    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2: there are more bases at y, 3: there are more bases at x
    for (c_i = 0; c_i < c_n; c_i++) {
        get_cigar_cell(idx, cc, c_i, &c, &cl);
        if (c == 0) {//match
            for (i = 0; i < cl; i++) {
                if(x_i == scanXbeg) c_err = 0;
                x_i++; y_i++;

                if(x_i == scanXend + 1) {
                    (*get_error) = c_err;
                    return 1;
                }
            }
        } else if (c == 1) {//mismatch
            for (i = 0; i < cl; i++) {
                if(x_i == scanXbeg) c_err = 0;
                x_i++; y_i++; c_err++;

                if(x_i == scanXend + 1) {
                    (*get_error) = c_err;
                    return 1;
                }
            }
        } else if (c == 2) {///y has more bases than x
            c_err += cl; y_i += cl;
        }
        else if (c == 3) {
            for (i = 0; i < cl; i++) {
                if(x_i == scanXbeg) c_err = 0;
                x_i++; c_err++;

                if(x_i == scanXend + 1) {
                    (*get_error) = c_err;
                    return 1;
                }
            }
        }
    }

    (*get_error) = c_err;
    return 0;
}

inline int move_gap_greedy(char* path, int path_i, int path_length, char* x, int x_i, char* y, int y_i, unsigned int* new_error)
{
    if(path[path_i] < 2)
    {
        return 0;
    }

    /**
         * 
    GGCG-TGTGCCTGT
        *
    GGCAATGTGCCTGT
        *
    00013000000000
    **/

    int flag = 0;

    char oper = path[path_i];
    

    if(oper == 3)///there are more x
    {
        path_i++;
        y_i--;
        for (; path_i < path_length && x_i >= 0 && y_i >= 0; path_i++, x_i--, y_i--)
        {
            if(path[path_i] == 2 || path[path_i] == 3 || (path[path_i] == 0 && x[x_i] != y[y_i]))
            {
                break;
            }
            else ///path[path_i] = 1 || path[path_i] = 0, exchange path[path_i] with path[path_i-1]
            {
                if(path[path_i] == 1 && x[x_i] == y[y_i])
                {
                    path[path_i - 1] = 0;
                    (*new_error)--;
                }
                else
                {
                    path[path_i - 1] = path[path_i];
                }

                path[path_i] = oper;
                

                flag = 1;
            }
            
        }
    }
    else if(oper == 2)///there are more y
    {
        path_i++;
        x_i--;
        for (; path_i < path_length && x_i >= 0 && y_i >= 0; path_i++, x_i--, y_i--)
        {
            if(path[path_i] == 2 || path[path_i] == 3 || (path[path_i] == 0 && x[x_i] != y[y_i]))
            {
                break;
            }
            else
            {

                if(path[path_i] == 1 && x[x_i] == y[y_i])
                {
                    path[path_i - 1] = 0;
                    (*new_error)--;
                }
                else
                {
                    path[path_i - 1] = path[path_i];
                }


                path[path_i] = oper;
                flag = 1;
            }
            
        }
    }

    return flag;
}

inline void generate_cigar(char* path, int path_length, window_list *idx, window_list_alloc *res, int* start, int* end, unsigned int* old_error,
    char* x, int x_len, char* y)
{
    // uint8_t debug_c; uint32_t debug_c_len;
    idx->cidx = res->c.n;
    if ((*old_error) == 0) {
        push_cigar_cell(res, 0, idx->x_end + 1 - idx->x_start);
        idx->clen = res->c.n - idx->cidx;
        // get_cigar_cell(idx, res, idx->clen-1, &debug_c, &debug_c_len);
        // assert(debug_c==0 && debug_c_len==(idx->x_end + 1 - idx->x_start));
        return;
    }

    ///0 is match, 1 is mismatch, 2 is up, 3 is left
	int32_t i = 0, pre_cl = 0, trem_p = -1; char pre_c = 5;    
    for (i = 0; i < path_length; i++) {
        if(path[i] == 1) {
            path[i] = 3;(*end)--; trem_p = i;
        } else {
            break;
        }
    }

    for (i = path_length - 1; i >= 0; i--) {
        if(path[i] == 1) {
            path[i] = 3; (*start)++;
        }
        else {
            break;
        }
    }


    // for (i = path_length - 1; i >= 0; i--)
    // {

    //     if (pre_ciga != path[i])
    //     {
    //         if (pre_ciga_length != 0)
    //         {
    //             result->cigar.C_L[result->cigar.length] = pre_ciga_length;
    //             result->cigar.C_C[result->cigar.length] = pre_ciga;
    //             result->cigar.length++;
    //         }

    //         pre_ciga = path[i];
    //         pre_ciga_length = 1;
    //     }
    //     else
    //     {
    //         pre_ciga_length++;
    //     }
    // }

    // if (pre_ciga_length != 0)
    // {
    //     result->cigar.C_L[result->cigar.length] = pre_ciga_length;
    //     result->cigar.C_C[result->cigar.length] = pre_ciga;
    //     result->cigar.length++;
    // }

    ///verify_cigar(x, x_len, y + (*start), (*end) - (*start) + 1, &(result->cigar), error);
    
    
    y = y + (*start);
    int32_t x_i = 0, y_i = 0;
    ///terminate_site = -1 in default
    for (i = path_length - 1; i > trem_p; i--) {
        if(path[i] == 0) {
            x_i++; y_i++;
        }
        else if(path[i] == 1) {
            x_i++; y_i++;
        }
        else if(path[i] == 2) {///there are more y
            move_gap_greedy(path, i, path_length, x, x_i, y, y_i, old_error);
            y_i++;
        }
        else if(path[i] == 3) {///there are more x
            move_gap_greedy(path, i, path_length, x, x_i, y, y_i, old_error);
            x_i++;
        }
    }



    pre_c = 5; pre_cl = 0;
    for (i = path_length - 1; i >= 0; i--) {
        if (pre_c != path[i]) {
            if (pre_cl != 0) {
                push_cigar_cell(res, pre_c, pre_cl);
                // get_cigar_cell(idx, res, res->c.n - idx->cidx - 1, &debug_c, &debug_c_len);
                // assert(debug_c==pre_c && debug_c_len==pre_cl);
            }
            pre_c = path[i]; pre_cl = 1;
        }
        else {
            pre_cl++;
        }
    }

    if (pre_cl != 0) {
        push_cigar_cell(res, pre_c, pre_cl);   
        // get_cigar_cell(idx, res, res->c.n - idx->cidx -1, &debug_c, &debug_c_len);
        // assert(debug_c==pre_c && debug_c_len==pre_cl);
    }

    idx->clen = res->c.n - idx->cidx;
    // if(verify_cigar(x, x_len, y, (*end) - (*start) + 1, &(result->cigar), *old_error))
    // {
    //     fprintf(stderr, "error\n");
    // }
}

int verify_cigar_2(char* x, int x_len, char* y, int y_len, Cigar_record* cigar, int error)
{
    int x_i, y_i, cigar_i;
    x_i = 0;
    y_i = 0;
    cigar_i = 0;
    int operation;
    int operationLen;
    int i;
    int cigar_error = 0;
    int flag_error = 0;
    int diff_i = 0;

    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2 means there are more y, 3 means there are more x
    while (cigar_i < (long long)cigar->length)
    {
        operation = Get_Cigar_Type(cigar->record[cigar_i]);
        operationLen = Get_Cigar_Length(cigar->record[cigar_i]);

        if (operation == 0)
        {
            for (i = 0; i < operationLen; i++)
            {

                if (x[x_i]!=y[y_i])
                {
                    ///fprintf(stderr, "error match\n");
                    flag_error = 1;
                }
                x_i++;
                y_i++;
            }
        }
        else if (operation == 1)
        {
            cigar_error += operationLen;
            for (i = 0; i < operationLen; i++)
            {

                if (x[x_i]==y[y_i])
                {
                    ///fprintf(stderr, "error mismatch, cigar_i: %d, x_i: %d, y_i: %d\n",cigar_i, x_i, y_i);
                    flag_error = 1;
                }

                if(Get_MisMatch_Base(cigar->lost_base[diff_i]) != y[y_i])
                {
                    // fprintf(stderr, "mismatch x: %c, y: %c, mis[%d]: %c\n", x[x_i],y[y_i],diff_i, 
                    // Get_MisMatch_Base(cigar->lost_base[diff_i]));
                }

                if(Get_Match_Base(cigar->lost_base[diff_i]) != x[x_i])
                {
                    // fprintf(stderr, "match x: %c, y: %c, deletion[%d]: %c\n", 
                    // x[x_i],y[y_i],diff_i, 
                    // Get_Match_Base(cigar->lost_base[diff_i]));
                }
               
                


                x_i++;
                y_i++;
                diff_i++;
            }
        }
        else if (operation == 2)
        {
            cigar_error += operationLen;
            for (i = 0; i < operationLen; i++)
            {
                if(cigar->lost_base[diff_i] != y[y_i])
                {
                    ///fprintf(stderr, "insertion x: %c, y: %c, insertion[%d]: %c\n", x[x_i],y[y_i],diff_i, 
                    ///cigar->lost_base[diff_i]);
                }
                y_i++;
                diff_i++;
            }
        }
        else if (operation == 3)
        {
            cigar_error += operationLen;

            for (i = 0; i < operationLen; i++)
            {
                if(cigar->lost_base[diff_i] != x[x_i])
                {
                    ///fprintf(stderr, "deletion x: %c, y: %c, deletion[%d]: %c\n", x[x_i],y[y_i],diff_i, 
                    ///cigar->lost_base[diff_i]);
                }
                x_i++;
                diff_i++;
            }
        }
        
        cigar_i++;
    }


    ///return;
    /**
    if (cigar_error != error)
    {
        fprintf(stderr, "error cigar_error: cigar_error: %d, error: %d\n", cigar_error, error);
        for (i = 0; i < cigar->length; i++)
        {
            operation = Get_Cigar_Type(cigar->record[i]);
            operationLen = Get_Cigar_Length(cigar->record[i]);
            fprintf(stderr, "%u: %u\n", operationLen, operation);
        }
        
    }
    **/
    
   
    if (flag_error == 1)
    {
        print_string(x, x_len);
        print_string(y, y_len);
        ///fprintf(stderr, "x_len: %d, y_len: %d, cigar_len: %d, error: %d\n", x_len, y_len, cigar->length, error);
        for (i = 0; i < (long long)cigar->length; i++)
        {
            operation = Get_Cigar_Type(cigar->record[i]);
            operationLen = Get_Cigar_Length(cigar->record[i]);
            ///fprintf(stderr, "%u: %u\n", operationLen, operation);
        }
    }
    
    
    return flag_error;
    
}

inline int fix_ul_boundary(char* x_string, long long x_len, int threshold,
long long total_y_start, long long local_y_start, long long local_y_end,
long long old_extra_begin, long long old_extra_end,
long long y_ID, long long Window_Len, const ul_idx_t *uref,
Correct_dumy* dumy, int y_strand, unsigned int old_error,
long long* r_total_y_start, int* r_start_site, int* r_end_site,
int* r_extra_begin, int* r_extra_end, unsigned int* r_error)
{

    
    int new_extra_begin, new_extra_end;
    long long new_y_start, new_y_length;
    int new_end_site, new_start_site;
    unsigned int new_error;
    char* y_string;


    int path_length;

    ///if the start pos at the left boundary 
    if(local_y_start == 0)
    {
        total_y_start = total_y_start + local_y_start;
        ///if local_y_start == 0 and old_extra_begin != 0
        ///this means total_y_start == 0, so shift to the left cannot get a new start pos
        if(old_extra_begin != 0)
        {
            return 0;
        }

        ///if the begining of alignment is 0, we should try to shift the window to find a better result
        ///shift to the left by threshold-1 bases
        if(!determine_overlap_region(threshold, total_y_start, y_ID, Window_Len, uref->ug->u.a[y_ID].len,
        &new_extra_begin, &new_extra_end, &new_y_start, &new_y_length))
        {
            return 0;
        }

        ///if new_y_start is equal to total_y_start, recalculate makes no sense
        if(new_y_start == total_y_start)
        {
            return 0;
        }

        fill_subregion_ul(dumy->overlap_region_fix, new_y_start, new_y_length, y_strand, uref, y_ID, 
        new_extra_begin, new_extra_end);

        y_string = dumy->overlap_region_fix;

        new_end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &new_error, &new_start_site,
                    &path_length, dumy->matrix_bit, dumy->path_fix, -1, -1);
        
        if (new_error != (unsigned int)-1 && new_error < old_error)
        {
            (*r_total_y_start) = new_y_start;
            (*r_start_site) = new_start_site;
            (*r_end_site) = new_end_site;
            (*r_extra_begin) = new_extra_begin;
            (*r_extra_end) = new_extra_end;
            (*r_error) = new_error;

            dumy->path_length = path_length;
            memcpy(dumy->path, dumy->path_fix, path_length);
            memcpy(dumy->overlap_region, dumy->overlap_region_fix, Window_Len);
            return 1;
        }
    }
    else if(local_y_end == Window_Len - 1)
    {
        ///if local_y_end == Window_Len - 1 and old_extra_end > 0
        ///this means local_y_end is the end of the y
        ///so shit to the right makes no sense
        if(old_extra_end != 0)
        {
            return 0;
        }
        long long total_y_end = total_y_start + local_y_end;

        total_y_start = total_y_end - x_len + 1;

        if(!determine_overlap_region(threshold, total_y_start, y_ID, Window_Len, uref->ug->u.a[y_ID].len,
        &new_extra_begin, &new_extra_end, &new_y_start, &new_y_length))
        {
            return 0;
        }

        if(new_y_start == total_y_end - local_y_end)
        {
            return 0;
        }

        fill_subregion_ul(dumy->overlap_region_fix, new_y_start, new_y_length, y_strand, uref, y_ID, 
        new_extra_begin, new_extra_end);

        y_string = dumy->overlap_region_fix;

        new_end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &new_error, &new_start_site,
                    &path_length, dumy->matrix_bit, dumy->path_fix, -1, -1);

        if (new_error != (unsigned int)-1 && new_error < old_error)
        {
            (*r_total_y_start) = new_y_start;
            (*r_start_site) = new_start_site;
            (*r_end_site) = new_end_site;
            (*r_extra_begin) = new_extra_begin;
            (*r_extra_end) = new_extra_end;
            (*r_error) = new_error;

            dumy->path_length = path_length;
            memcpy(dumy->path, dumy->path_fix, path_length);
            memcpy(dumy->overlap_region, dumy->overlap_region_fix, Window_Len);
            return 1;
        }


    }
    return 0;
}


inline int fix_boundary(char* x_string, long long x_len, int threshold,
long long total_y_start, long long local_y_start, long long local_y_end,
long long old_extra_begin, long long old_extra_end,
long long y_ID, long long Window_Len, All_reads* R_INF,
Correct_dumy* dumy, int y_strand, unsigned int old_error,
long long* r_total_y_start, int* r_start_site, int* r_end_site,
int* r_extra_begin, int* r_extra_end, unsigned int* r_error)
{

    
    int new_extra_begin, new_extra_end;
    long long new_y_start, new_y_length;
    int new_end_site, new_start_site;
    unsigned int new_error;
    char* y_string;


    int path_length;

    ///if the start pos at the left boundary 
    if(local_y_start == 0)
    {
        total_y_start = total_y_start + local_y_start;
        ///if local_y_start == 0 and old_extra_begin != 0
        ///this means total_y_start == 0, so shift to the left cannot get a new start pos
        if(old_extra_begin != 0)
        {
            return 0;
        }

        ///if the begining of alignment is 0, we should try to shift the window to find a better result
        ///shift to the left by threshold-1 bases
        if(!determine_overlap_region(threshold, total_y_start, y_ID, Window_Len, Get_READ_LENGTH((*R_INF), y_ID),
        &new_extra_begin, &new_extra_end, &new_y_start, &new_y_length))
        {
            return 0;
        }

        ///if new_y_start is equal to total_y_start, recalculate makes no sense
        if(new_y_start == total_y_start)
        {
            return 0;
        }

        fill_subregion(dumy->overlap_region_fix, new_y_start, new_y_length, y_strand, R_INF, y_ID, 
        new_extra_begin, new_extra_end);

        y_string = dumy->overlap_region_fix;

        new_end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &new_error, &new_start_site,
                    &path_length, dumy->matrix_bit, dumy->path_fix, -1, -1);
        
        if (new_error != (unsigned int)-1 && new_error < old_error)
        {
            (*r_total_y_start) = new_y_start;
            (*r_start_site) = new_start_site;
            (*r_end_site) = new_end_site;
            (*r_extra_begin) = new_extra_begin;
            (*r_extra_end) = new_extra_end;
            (*r_error) = new_error;

            dumy->path_length = path_length;
            memcpy(dumy->path, dumy->path_fix, path_length);
            memcpy(dumy->overlap_region, dumy->overlap_region_fix, Window_Len);
            return 1;
        }
    }
    else if(local_y_end == Window_Len - 1)
    {
        ///if local_y_end == Window_Len - 1 and old_extra_end > 0
        ///this means local_y_end is the end of the y
        ///so shit to the right makes no sense
        if(old_extra_end != 0)
        {
            return 0;
        }
        long long total_y_end = total_y_start + local_y_end;

        total_y_start = total_y_end - x_len + 1;

        if(!determine_overlap_region(threshold, total_y_start, y_ID, Window_Len, Get_READ_LENGTH((*R_INF), y_ID),
        &new_extra_begin, &new_extra_end, &new_y_start, &new_y_length))
        {
            return 0;
        }

        if(new_y_start == total_y_end - local_y_end)
        {
            return 0;
        }

        fill_subregion(dumy->overlap_region_fix, new_y_start, new_y_length, y_strand, R_INF, y_ID, 
        new_extra_begin, new_extra_end);

        y_string = dumy->overlap_region_fix;

        new_end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &new_error, &new_start_site,
                    &path_length, dumy->matrix_bit, dumy->path_fix, -1, -1);

        if (new_error != (unsigned int)-1 && new_error < old_error)
        {
            (*r_total_y_start) = new_y_start;
            (*r_start_site) = new_start_site;
            (*r_end_site) = new_end_site;
            (*r_extra_begin) = new_extra_begin;
            (*r_extra_end) = new_extra_end;
            (*r_error) = new_error;

            dumy->path_length = path_length;
            memcpy(dumy->path, dumy->path_fix, path_length);
            memcpy(dumy->overlap_region, dumy->overlap_region_fix, Window_Len);
            return 1;
        }


    }
    return 0;
}

inline char *return_str_seq(char *buf, int64_t s, int64_t pri_l, uint8_t rev, hpc_t *hpc_g, const ul_idx_t *uref, int64_t id, int64_t aux_beg, int64_t aux_end)
{
    if(!hpc_g) {
        memset(buf, 'N', aux_beg);
        retrieve_u_seq(NULL, buf+aux_beg, &(uref->ug->u.a[id]), rev, s, pri_l, NULL);
        memset(buf+aux_beg+pri_l, 'N', aux_end);
        return buf;
    } else {
        char *z = hpc_str(*hpc_g, id, rev);
        if((aux_beg == 0) && (aux_end == 0)) {
            return z+s;
        } else {
            memset(buf, 'N', aux_beg);
            memcpy(buf+aux_beg, z+s, pri_l);
            memset(buf+aux_beg+pri_l, 'N', aux_end);
            return buf;
        }
    }
}

inline char *return_str_seq_exz(char *buf, int64_t s, int64_t pri_l, uint8_t rev, hpc_t *hpc_g, const ul_idx_t *uref, int64_t id)
{
    if(!hpc_g) {
        retrieve_u_seq(NULL, buf, &(uref->ug->u.a[id]), rev, s, pri_l, NULL);
        return buf;
    } else {
        return hpc_str(*hpc_g, id, rev) + s;
    }
}

///cannot use tstr in-place
inline int recal_boundary(char* qstr, char* tstr1, int64_t ql, int64_t thres,
int64_t global_ts0, int64_t local_ts0, int64_t local_te0,
int64_t aux_beg0, int64_t aux_end0, unsigned int err0, 
int64_t tid, int64_t aln_l, uint32_t rev, 
Correct_dumy* dumy, All_reads* rref, hpc_t *hpc_g, const ul_idx_t *uref,
int64_t* global_ts1, int* local_ts1, int* local_te1,
int64_t* aux_beg1, int64_t* aux_end1, unsigned int* err1)
{
    int64_t ts, t_tot_l, aux_beg, aux_end, t_pri_l, t_end; 
    char *q_string = qstr, *t_string; unsigned int error = (unsigned int)-1; 
    int r_ts = 0, path_length = 0;
    if(hpc_g) t_tot_l = hpc_len(*hpc_g, tid);
    else if(uref) t_tot_l = uref->ug->u.a[tid].len;
    else t_tot_l = Get_READ_LENGTH((*rref), tid);

    if(local_ts0 == 0) {//left boundary
        if(aux_beg0 > 0) return 0;///shift to the left cannot get a new start pos
        ts = global_ts0;
    } else if((local_te0 + 1) == aln_l) {//right boundary
        if(aux_end0 > 0) return 0;///shift to the right cannot get a new start pos
        ts = global_ts0 + local_te0 - ql + 1;
    } else {
        return 0;
    }
    if(!init_waln(thres, ts, t_tot_l, aln_l, &aux_beg, &aux_end, &ts, &t_pri_l)) return 0;
    if(ts == global_ts0) return 0;//unchanged, make no sense

    if(rref) {
        fill_subregion(tstr1, ts, t_pri_l, rev, rref, tid, aux_beg, aux_end); t_string = tstr1;
    } else {
        t_string = return_str_seq(tstr1, ts, t_pri_l, rev, hpc_g, uref, tid, aux_beg, aux_end);
    }

    t_end = Reserve_Banded_BPM_PATH(t_string, aln_l, q_string, ql, thres, &error, &r_ts,
                &path_length, dumy->matrix_bit, dumy->path_fix, -1, -1);

    if (error != (unsigned int)-1 && error < err0) {
        (*global_ts1) = ts;
        (*local_ts1) = r_ts;
        (*local_te1) = t_end;
        (*aux_beg1) = aux_beg;
        (*aux_end1) = aux_end;
        (*err1) = error;

        dumy->path_length = path_length;
        memcpy(dumy->path, dumy->path_fix, path_length);
        // memcpy(tstr0, t_string, aln_l);
        return 1;
    }
    return 0;
}

///cannot use tstr in-place
inline int recal_boundary_exz(char* qstr, char* tstr, int64_t ql0, int64_t tl0, int64_t thres,
int64_t toff, int64_t ts0, int64_t te0, int64_t err0, 
int64_t tid, uint32_t rev, bit_extz_t *exz,
All_reads* rref, hpc_t *hpc_g, const ul_idx_t *uref, 
int64_t *ts_r, int64_t *aux_beg_r, int64_t *aux_end_r)
{
    int64_t ts, tl, t_tot_l, aux_beg, aux_end, t_pri_l, aln_l = ql0 + (thres << 1);
    char *q_string = qstr, *t_string; 
    if(hpc_g) t_tot_l = hpc_len(*hpc_g, tid);
    else if(uref) t_tot_l = uref->ug->u.a[tid].len;
    else t_tot_l = Get_READ_LENGTH((*rref), tid);

    if(ts0 == 0) {//left boundary
        ts = toff;
    } else if((te0 + 1) == tl0) {//right boundary
        ts = toff + te0 - ql0 + 1;
    } else {
        return 0;
    }
    if(!init_waln(thres, ts, t_tot_l, aln_l, &aux_beg, &aux_end, &ts, &t_pri_l)) return 0;
    if(ts == toff && tl0 == t_pri_l) return 0;//unchanged, make no sense

    tl = t_pri_l;
    if(rref) {
        recover_UC_Read_sub_region(tstr, ts, tl, rev, rref, tid); t_string = tstr;
    } else {
        t_string = return_str_seq_exz(tstr, ts, tl, rev, hpc_g, uref, tid);
    }

    clear_align(*exz); 
    ed_band_cal_semi_64_w_absent_diag_trace(t_string, tl, q_string, ql0, thres, aux_beg, exz);

    if(is_align(*exz) && exz->err < err0) {
        (*aux_beg_r) = aux_beg;
        (*aux_end_r) = aux_end;
        (*ts_r) = ts;
        return 1;
    }
    return 0;
}

inline char *update_des_str(char *des, int64_t s, int64_t pri_l, uint8_t rev, All_reads *rref, hpc_t *hpc_g, 
    const ul_idx_t *uref, int64_t id, int64_t aux_beg, int64_t aux_end, char *src)
{
    if(src) {
        // memcpy(des, src, (pri_l+aux_beg+aux_end));
        // return des;
        return src;
    } else {
        if(rref) {
            fill_subregion(des, s, pri_l, rev, rref, id, aux_beg, aux_end); 
            return des;
        } else {
            return return_str_seq(des, s, pri_l, rev, hpc_g, uref, id, aux_beg, aux_end);
        }
    }
}

/**
void debug_scan_cigar(overlap_region* sub_list)
{
    long long i;
    int f_err, b_err, fLen, xLen;
    for (i = 0; i < (long long)sub_list->w_list_length; i++)
    {
        if(sub_list->w_list[i].y_end == -1 || sub_list->w_list[i].cigar.length == -1)
        {
            continue;
        }
        xLen = sub_list->w_list[i].x_end - sub_list->w_list[i].x_start + 1;

        scan_cigar(&(sub_list->w_list[i].cigar), &b_err, 
        xLen, 1);
        scan_cigar(&(sub_list->w_list[i].cigar), &f_err, 
        xLen, 0);

        if(b_err != sub_list->w_list[i].error || f_err != sub_list->w_list[i].error)
        {
            fprintf(stderr, "error\n");
        }

        scan_cigar(&(sub_list->w_list[i].cigar), &b_err, 
        WINDOW, 1);
        scan_cigar(&(sub_list->w_list[i].cigar), &f_err, 
        WINDOW, 0);

        if(b_err != sub_list->w_list[i].error || f_err != sub_list->w_list[i].error)
        {
            fprintf(stderr, "error\n");
        }

        scan_cigar_interval(&(sub_list->w_list[i].cigar), &b_err, 0, xLen-1);
        if(b_err != sub_list->w_list[i].error)
        {
            fprintf(stderr, "error\n");
        }

        fLen = xLen / 3;
        scan_cigar(&(sub_list->w_list[i].cigar), &f_err, fLen, 0);
        scan_cigar_interval(&(sub_list->w_list[i].cigar), &b_err, 0, fLen-1);
        if(f_err != b_err)
        {
            fprintf(stderr, "error\n");
        }

        fLen = xLen / 3;
        scan_cigar(&(sub_list->w_list[i].cigar), &f_err, fLen, 1);
        scan_cigar_interval(&(sub_list->w_list[i].cigar), &b_err, xLen-fLen, xLen-1);
        if(f_err != b_err)
        {
            fprintf(stderr, "\nerror\n");
            fprintf(stderr, "b_err: %d, f_err: %d\n",b_err, f_err);
            long long j;
            for (j = 0; j < sub_list->w_list[i].cigar.length; j++)
            {
                fprintf(stderr, "len: %d, opera: %d\n", 
                sub_list->w_list[i].cigar.C_L[j], sub_list->w_list[i].cigar.C_C[j]);
            }
        }

        // bLen = xLen / 3;
        // fLen = xLen - bLen;

        // scan_cigar(&(sub_list->w_list[i].cigar), &b_err, 
        // bLen, 1);
        // scan_cigar(&(sub_list->w_list[i].cigar), &f_err, 
        // fLen, 0);

        // if(b_err + f_err != sub_list->w_list[i].error)
        // {
        //     fprintf(stderr, "\nsub_list->w_list[i].error: %d, bLen: %d, b_err: %d, fLen: %d, f_err: %d\n",
        //     sub_list->w_list[i].error, bLen, b_err, fLen, f_err);
        //     long long j;
        //     for (j = 0; j < sub_list->w_list[i].cigar.length; j++)
        //     {
        //         fprintf(stderr, "len: %d, opera: %d\n", 
        //         sub_list->w_list[i].cigar.C_L[j], sub_list->w_list[i].cigar.C_C[j]);
        //     }
            
        // }
    }
}
**/

void calculate_boundary_cigars(overlap_region* z, All_reads* R_INF, Correct_dumy* dumy, UC_Read* g_read, double e_rate)
{
    assert(z->w_list.n > 0);
    int64_t nw = z->w_list.n;
    resize_window_list_alloc(&(z->boundary_cigars), nw - 1);
    int64_t y_id = z->y_id, y_strand = z->y_pos_strand;
    int64_t y_readLen = Get_READ_LENGTH((*R_INF), y_id);
    int64_t i, y_distance, f_err = -1, b_err = -1, m_error;
    int64_t scanLen = 10, boundaryLen = 200;
    int64_t single_sideLen = boundaryLen/2;
    int64_t force_useless_side = single_sideLen/2;
    int64_t L_useless_side, R_useless_side, alpha = 1;
    long long y_start, x_start, x_end, yLen, xLen, leftLen, rightLen, threshold, o_len;
    char *x_string = NULL, *y_string = NULL;
    int end_site, real_y_start, extra_begin, extra_end;
    unsigned int error;
    z->boundary_cigars.n = nw - 1;
    ///the (i)-th boundary between the (i)-th window and the (i+1)-th window
    ///that means it includes (the tail of (i)-th window) and (the header of (i+1)-th window)
    ///note the (i)-th boundary is calculated at the (i)-th window
    for (i = 0; i + 1 < nw; i++) {
        ///if both of the two windows are not aligned
        ///it is not necessary to calculate the boundary
        if(z->w_list.a[i].y_end == -1 || z->w_list.a[i+1].y_end == -1) {
            z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
            continue;
        }
        ///y_distance can be less than 0, or larger than 0
        y_distance = (int64_t)z->w_list.a[i+1].y_start - (int64_t)z->w_list.a[i].y_end - 1;

        ///if two windows are aligned
        if(z->w_list.a[i].y_end != -1 && z->w_list.a[i+1].y_end != -1 && y_distance == 0) {
            ///scan backward
            scan_cigar(&(z->w_list.a[i]), &(z->w_list), &b_err, scanLen, 1);
            ///scan forward
            scan_cigar(&(z->w_list.a[i+1]), &(z->w_list), &f_err, scanLen, 0);
            if(b_err == 0 && f_err == 0) {
                z->boundary_cigars.a[i].error = -2; z->boundary_cigars.a[i].y_end = -1;
                continue;
            }
        }

        
        if(z->w_list.a[i].y_end != -1) {
            y_start = z->w_list.a[i].y_end; x_start = z->w_list.a[i].x_end;
        }///if the (i)-th window is not matched, have a look at the (i+1)-th window
        else if(z->w_list.a[i+1].y_end != -1) {
            y_start = z->w_list.a[i+1].y_start; x_start = z->w_list.a[i+1].x_start;
        }///if both of these two windows are not matched, directly skip
        else {
            z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
            continue;
        }

        ///it seems we don't need to record x_start and y_start
        z->boundary_cigars.a[i].extra_begin = x_start;
        z->boundary_cigars.a[i].extra_end = y_start;

        ///leftLen and rightLen are used for x        
        ///x should be at [sub_list->w_list[i].x_start, sub_list->w_list[i+1].x_end]
        ///y shouldn't have limitation
        ///note that the x_start and x_end should not be -1 in any case
        ///up to now, x_start and y_start are not -1
        ///leftLen does not include x_start itself, rightLen does
        ///gnerally speaking, rightLen should be always larger than leftLen
        leftLen = MIN(MIN((x_start - (int64_t)z->w_list.a[i].x_start), y_start), single_sideLen);
        rightLen = MIN(MIN(((int64_t)z->w_list.a[i+1].x_end + 1 - x_start), y_readLen - y_start), single_sideLen);

        ///xLen should be the sum length of two windows
        xLen = leftLen + rightLen;
        x_start = x_start - leftLen;
        x_end = x_start + xLen - 1;
        y_start = y_start - leftLen;

        ///if we don't have enough leftLen and rightLen
        // if(leftLen <= useless_side || rightLen <= useless_side)
        // {
        //     sub_list->boundary_cigars.buffer[i].error = -1;
        //     sub_list->boundary_cigars.buffer[i].y_end = -1;
        //     continue;
        // }


        threshold = xLen * e_rate/**asm_opt.max_ov_diff_ec**/;
        threshold = Adjust_Threshold(threshold, xLen);
        threshold = double_error_threshold(threshold, xLen);

        yLen = xLen + (threshold << 1);
        if(!determine_overlap_region(threshold, y_start, y_id, yLen, Get_READ_LENGTH((*R_INF), y_id), 
                    &extra_begin, &extra_end, &y_start, &o_len)) {
            z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
            continue;
        }

        if(o_len < xLen) {
            z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
            continue;
        }

        fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, R_INF, y_id, extra_begin, extra_end);
        
        x_string = g_read->seq + x_start;
        y_string = dumy->overlap_region;

        end_site = Reserve_Banded_BPM_PATH(y_string, yLen, x_string, xLen, threshold, &error, 
        &real_y_start, &(dumy->path_length), dumy->matrix_bit, dumy->path, -1, -1);

        ///means this window is matched
        if (error!=(unsigned int)-1) {
            z->boundary_cigars.a[i].x_start = x_start;
            z->boundary_cigars.a[i].x_end = x_end;

            generate_cigar(dumy->path, dumy->path_length, &(z->boundary_cigars.a[i]), &(z->boundary_cigars),
                        &real_y_start, &end_site, &error, x_string, xLen, y_string);
            ///should not adjust cigar here, adjust cigar may cause problem
            ///that is not what we want

            ///y_distance can be less than 0, or larger than 0
            ///please if one of the two windows is not matched, 
            ///y_distance may have potential problems 
            if(y_distance < 0) y_distance = y_distance * (-1);
            ///leftLen, rightLen
            // if(leftLen <= useless_side || rightLen <= useless_side)
            // {
            //     sub_list->boundary_cigars.buffer[i].error = -1;
            //     sub_list->boundary_cigars.buffer[i].y_end = -1;
            //     continue;
            // }
            L_useless_side = R_useless_side = force_useless_side;

            ///first window
            if((i == 0) && (x_start == (int64_t)z->w_list.a[0].x_start)) {
                L_useless_side = 0;
            }
            ///last window
            if((i == (int64_t)(z->w_list.n) - 2) && 
                    (x_end == (long long)(z->w_list.a[(int64_t)(z->w_list.n)-1].x_end))) {
                R_useless_side = 0;
            }

            if(leftLen <= L_useless_side || rightLen <= R_useless_side) {
                z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
                z->boundary_cigars.c.n = z->boundary_cigars.a[i].cidx;
                continue;
            }



            ///up to now, if we require (i)-th window and (i+1)-th window are matched
            ///boundary_cigars.buffer[i].cigar, w_list[i].cigar and w_list[i+1].cigar are avaiable
            ///get the error excluding the first and the last useless_side bases
            scan_cigar_interval(&(z->boundary_cigars.a[i]), &(z->boundary_cigars), &m_error, L_useless_side, xLen-R_useless_side-1);
            scan_cigar(&(z->w_list.a[i]), &(z->w_list), &b_err, leftLen-L_useless_side, 1);
            scan_cigar(&(z->w_list.a[i+1]), &(z->w_list), &f_err, rightLen-R_useless_side, 0);

            if(f_err + b_err + y_distance + alpha < m_error) {
                z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
                z->boundary_cigars.c.n = z->boundary_cigars.a[i].cidx;
                continue;
            }

            z->boundary_cigars.a[i].error = error;
            z->boundary_cigars.a[i].y_start = y_start + real_y_start - extra_begin;
            z->boundary_cigars.a[i].y_end = y_start + end_site - extra_begin;

            z->boundary_cigars.a[i].x_start = x_start;
            z->boundary_cigars.a[i].x_end = x_end;
            ///sub_list->boundary_cigars.buffer[i].error_threshold = useless_side;
            z->boundary_cigars.a[i].extra_begin = L_useless_side;
            z->boundary_cigars.a[i].extra_end = R_useless_side;
        }
        else {
            z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
            continue;
        }     
    }

}


void calculate_ul_boundary_cigars(overlap_region* z, const ul_idx_t *uref, Correct_dumy* dumy, 
UC_Read* g_read, double max_ov_diff_ec, long long blockLen)
{
    assert(z->w_list.n > 0);
    int64_t nw = z->w_list.n;
    resize_window_list_alloc(&(z->boundary_cigars), nw - 1);
    int64_t y_id = z->y_id;
    int64_t y_strand = z->y_pos_strand;
    int64_t y_readLen = uref->ug->u.a[y_id].len;
    int64_t i, y_distance;
    int64_t f_err, b_err, m_error, scanLen = 10;
    int64_t boundaryLen = WINDOW_UL_BOUND_RATE*blockLen;
    boundaryLen >>= 2; boundaryLen <<= 2;
    if(boundaryLen < WINDOW_UL_BOUND) boundaryLen = WINDOW_UL_BOUND;
    int64_t single_sideLen = boundaryLen/2;
    int64_t force_useless_side = single_sideLen/2;
    int64_t L_useless_side, R_useless_side;
    int64_t alpha = 1;
    long long y_start, x_start, x_end, yLen, xLen, leftLen, rightLen, threshold, o_len;
    int extra_begin, extra_end, end_site, real_y_start;
    char* x_string;
    char* y_string;
    unsigned int error;
    z->boundary_cigars.n = nw - 1;
    ///the (i)-th boundary between the (i)-th window and the (i+1)-th window
    ///that means it includes (the tail of (i)-th window) and (the header of (i+1)-th window)
    ///note the (i)-th boundary is calculated at the (i)-th window
    for (i = 0; i + 1 < nw; i++) {
        ///if both of the two windows are not aligned
        ///it is not necessary to calculate the boundary
        ///if(sub_list->w_list[i].y_end == -1 && sub_list->w_list[i+1].y_end == -1)
        if(z->w_list.a[i].y_end == -1 || z->w_list.a[i+1].y_end == -1) {
            z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
            continue;
        }

        ///note if w_list[i+1].y_start or sub_list->w_list[i].y_end is -1
        ///y_distance might have some problems at the last of this function
        ///we need to deal with it carefully
        y_distance = (int64_t)z->w_list.a[i+1].y_start - (int64_t)z->w_list.a[i].y_end - 1;

        ///if two windows are aligned
        if(z->w_list.a[i].y_end != -1 && z->w_list.a[i+1].y_end != -1 && y_distance == 0) {
            ///scan backward
            scan_cigar(&(z->w_list.a[i]), &(z->w_list), &b_err, scanLen, 1);
            ///scan forward
            scan_cigar(&(z->w_list.a[i+1]), &(z->w_list), &f_err, scanLen, 0);
            if(b_err == 0 && f_err == 0) {
                z->boundary_cigars.a[i].error = -2; z->boundary_cigars.a[i].y_end = -1;
                continue;
            }
        }


        ///y_distance can be less than 0, or larger than 0
        if(z->w_list.a[i].y_end != -1) {
            y_start = z->w_list.a[i].y_end; x_start = z->w_list.a[i].x_end;
        }///if the (i)-th window is not matched, have a look at the (i+1)-th window
        else if(z->w_list.a[i+1].y_end != -1) {
            y_start = z->w_list.a[i+1].y_start; x_start = z->w_list.a[i+1].x_start;
        }///if both of these two windows are not matched, directly skip
        else {
            z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
            continue;
        }

        ///it seems we don't need to record x_start and y_start
        z->boundary_cigars.a[i].extra_begin = x_start;
        z->boundary_cigars.a[i].extra_end = y_start;

        ///leftLen and rightLen are used for x        
        ///x should be at [sub_list->w_list[i].x_start, sub_list->w_list[i+1].x_end]
        ///y shouldn't have limitation
        ///note that the x_start and x_end should not be -1 in any case
        ///up to now, x_start and y_start are not -1
        ///leftLen does not include x_start itself, rightLen does
        ///gnerally speaking, rightLen should be always larger than leftLen
        leftLen = MIN(MIN((x_start - (long long)z->w_list.a[i].x_start), y_start), single_sideLen);
        rightLen = MIN(MIN(((long long)z->w_list.a[i+1].x_end + 1 - x_start), y_readLen - y_start), single_sideLen);

        ///xLen should be the sum length of two windows
        xLen = leftLen + rightLen;
        x_start = x_start - leftLen;
        x_end = x_start + xLen - 1;
        y_start = y_start - leftLen;

        ///if we don't have enough leftLen and rightLen
        // if(leftLen <= useless_side || rightLen <= useless_side)
        // {
        //     sub_list->boundary_cigars.buffer[i].error = -1;
        //     sub_list->boundary_cigars.buffer[i].y_end = -1;
        //     continue;
        // }


        threshold = xLen * max_ov_diff_ec;
        threshold = Adjust_Threshold(threshold, xLen);
        threshold = double_ul_error_threshold(threshold, xLen);

        yLen = xLen + (threshold << 1);
        if(!determine_overlap_region(threshold, y_start, y_id, yLen, uref->ug->u.a[y_id].len, 
                    &extra_begin, &extra_end, &y_start, &o_len)) {
            z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
            continue;
        }

        if(o_len < xLen) {
            z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
            continue;
        }

        fill_subregion_ul(dumy->overlap_region, y_start, o_len, y_strand, 
                    uref, y_id, extra_begin, extra_end);
        
        x_string = g_read->seq + x_start;
        y_string = dumy->overlap_region;

        end_site = Reserve_Banded_BPM_PATH(y_string, yLen, x_string, xLen, threshold, &error, 
        &real_y_start, &(dumy->path_length), dumy->matrix_bit, dumy->path, -1, -1);

        ///means this window is matched
        if (error!=(unsigned int)-1) {
            z->boundary_cigars.a[i].x_start = x_start; z->boundary_cigars.a[i].x_end = x_end;

            generate_cigar(dumy->path, dumy->path_length, &(z->boundary_cigars.a[i]), &(z->boundary_cigars),
                        &real_y_start, &end_site, &error, x_string, xLen, y_string);
            ///should not adjust cigar here, adjust cigar may cause problem
            ///that is not what we want

            ///y_distance can be less than 0, or larger than 0
            ///please if one of the two windows is not matched, 
            ///y_distance may have potential problems 
            if(y_distance < 0) y_distance = y_distance * (-1);
            ///leftLen, rightLen
            // if(leftLen <= useless_side || rightLen <= useless_side)
            // {
            //     sub_list->boundary_cigars.buffer[i].error = -1;
            //     sub_list->boundary_cigars.buffer[i].y_end = -1;
            //     continue;
            // }
            L_useless_side = R_useless_side = force_useless_side;

            ///first window
            if((i == 0) && (x_start == (long long)z->w_list.a[0].x_start)) {
                L_useless_side = 0;
            }
            ///last window
            if((i == (int64_t)(z->w_list.n) - 2) && (x_end == (z->w_list.a[(int64_t)z->w_list.n - 1].x_end))) {
                R_useless_side = 0;
            }

            if(leftLen <= L_useless_side || rightLen <= R_useless_side) {
                z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
                z->boundary_cigars.c.n = z->boundary_cigars.a[i].cidx;
                continue;
            }



            ///up to now, if we require (i)-th window and (i+1)-th window are matched
            ///boundary_cigars.buffer[i].cigar, w_list[i].cigar and w_list[i+1].cigar are avaiable
            ///get the error excluding the first and the last useless_side bases
            scan_cigar_interval(&(z->boundary_cigars.a[i]), &(z->boundary_cigars), &m_error, 
            L_useless_side, xLen-R_useless_side-1);
            scan_cigar(&(z->w_list.a[i]), &(z->w_list), &b_err, leftLen-L_useless_side, 1);
            scan_cigar(&(z->w_list.a[i+1]), &(z->w_list), &f_err, rightLen-R_useless_side, 0);

            if(f_err + b_err + y_distance + alpha < m_error) {
                z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
                z->boundary_cigars.c.n = z->boundary_cigars.a[i].cidx;
                continue;
            }

            z->boundary_cigars.a[i].error = error;
            z->boundary_cigars.a[i].y_start = y_start + real_y_start - extra_begin;
            z->boundary_cigars.a[i].y_end = y_start + end_site - extra_begin;

            z->boundary_cigars.a[i].x_start = x_start;
            z->boundary_cigars.a[i].x_end = x_end;
            ///sub_list->boundary_cigars.buffer[i].error_threshold = useless_side;
            z->boundary_cigars.a[i].extra_begin = L_useless_side;
            z->boundary_cigars.a[i].extra_end = R_useless_side;
        }
        else
        {
            z->boundary_cigars.a[i].error = -1;
            z->boundary_cigars.a[i].y_end = -1;
            continue;
        }
         
    }

}


/**
void debug_window_cigar(overlap_region_alloc* overlap_list, UC_Read* g_read, Correct_dumy* dumy,
All_reads* R_INF, int test_window, int test_boundary)
{
    uint64_t i, j, y_id, y_strand;
    char* x_string;
    char* y_string;
    long long x_start;
    long long x_end;
    long long x_len;
    long long y_start;
    long long y_end;
    long long y_len;

    for (j = 0; j < overlap_list->length; j++)
    {
        y_id = overlap_list->list[j].y_id;
        y_strand = overlap_list->list[j].y_pos_strand;
        if(overlap_list->list[j].is_match == 1)
        {

            if(test_window == 1)
            {
                for (i = 0; i < overlap_list->list[j].w_list_length; i++)
                {
                    if(overlap_list->list[j].w_list[i].y_end != -1)
                    {
                        ///there is no problem for x
                        x_start = overlap_list->list[j].w_list[i].x_start;
                        x_end = overlap_list->list[j].w_list[i].x_end;
                        x_len = x_end - x_start + 1;

                        x_string = g_read->seq + x_start;

                        y_start = overlap_list->list[j].w_list[i].y_start;
                        y_end = overlap_list->list[j].w_list[i].y_end;
                        y_len = y_end - y_start + 1;

                        recover_UC_Read_sub_region(dumy->overlap_region, y_start, y_len, y_strand, R_INF, y_id);
                        y_string = dumy->overlap_region;


                        if(verify_cigar(x_string, x_len, y_string, y_len, &overlap_list->list[j].w_list[i].cigar, 
                        overlap_list->list[j].w_list[i].error))
                        {
                            fprintf(stderr, "error\n");
                        }
                    }
                }
            }
            

            
            if(test_boundary == 1)
            {
                for (i = 0; i < (uint64_t)overlap_list->list[j].boundary_cigars.length; i++)
                {
                    if(overlap_list->list[j].boundary_cigars.buffer[i].y_end != -1)
                    {
                        x_start = overlap_list->list[j].boundary_cigars.buffer[i].x_start;
                        x_end = overlap_list->list[j].boundary_cigars.buffer[i].x_end;
                        x_len = x_end - x_start + 1;
                        x_string = g_read->seq + x_start;

                        y_start = overlap_list->list[j].boundary_cigars.buffer[i].y_start;
                        y_end = overlap_list->list[j].boundary_cigars.buffer[i].y_end;
                        y_len = y_end - y_start + 1;

                        recover_UC_Read_sub_region(dumy->overlap_region, y_start, y_len, y_strand, 
                        R_INF, y_id);
                        y_string = dumy->overlap_region;


                        if(verify_cigar(x_string, x_len, y_string, y_len, 
                        &overlap_list->list[j].boundary_cigars.buffer[i].cigar, 
                        overlap_list->list[j].boundary_cigars.buffer[i].error))
                        {
                            fprintf(stderr, "error\n");
                        }
                    }
                }
            }


            if(test_window == 1 && test_boundary == 1)
            {
                if(overlap_list->list[j].w_list_length != 
                        (uint64_t)(overlap_list->list[j].boundary_cigars.length + 1))
                {
                    fprintf(stderr, "error\n");
                }
            }
            
        }
    }
}
**/
int64_t get_adjust_winid(overlap_region *z, int64_t win_beg, int64_t win_len)
{   
    int64_t win_id, k;
    win_id = (win_beg-((z->x_pos_s/win_len)*win_len))/win_len;
    if((uint64_t)win_id < z->w_list.n && z->w_list.a[win_id].x_start == win_beg) return win_id;
    if(z->w_list.n == 0) return -1;
    // if(z->w_list.a[win_id].x_start <= win_beg) {
    //     fprintf(stderr, "z->w_list.n::%u, z->w_list.a[%ld].x_start::%d, win_beg::%ld\n", 
    //     (uint32_t)z->w_list.n, win_id, z->w_list.a[win_id].x_start, win_beg);
    // }
    if((uint64_t)win_id > z->w_list.n) win_id = z->w_list.n;
    // assert((z->w_list.a[win_id].x_start > win_beg);
    for (k = win_id - 1; k >= 0; k--) {
        // if(k < 0 || k >= (int64_t)z->w_list.n) fprintf(stderr, "win_id::%ld, k::%ld, z->w_list.n::%ld\n", win_id, k, (int64_t)z->w_list.n);
        if(z->w_list.a[k].x_start == win_beg) return k;
        if(z->w_list.a[k].x_start < win_beg) return -1;
    }
    return -1;
}

void set_herror_win(overlap_region_alloc* ovlp, Correct_dumy* du, kvec_t_u64_warp* v_idx, double max_ov_diff_ec, int64_t rLen, int64_t blockLen)
{
    Window_Pool w_inf; int32_t flag = 0; uint64_t cID, mm, fc, fw, idx_n, idx_i;
    init_Window_Pool(&w_inf, rLen, blockLen, (int)(1.0/max_ov_diff_ec));
    long long window_start, window_end; int64_t i, k, mLen, w_list_id, ws, we; 

    idx_n = get_num_wins(0, rLen, blockLen); idx_i = 0;
    kv_resize(uint64_t, v_idx->a, idx_n); v_idx->a.n = idx_n; 

    while(get_Window(&w_inf, &window_start, &window_end) && flag != -2) {
        du->length = du->lengthNT = 0;
        flag = get_interval(window_start, window_end, ovlp, du, w_inf.window_length);
        switch (flag) {
            case 1:    ///no match here
                break;
            case 0:    ///no match here
                break;
            case -2: ///if flag == -2, loop would be terminated
                break;
        }

        v_idx->a.a[idx_i++] = ((uint64_t)(v_idx->a.n))<<32;
        for (i = 0; i < (int64_t)du->length; i++) {
            cID = (uint32_t)du->overlapID[i];
            if(ovlp->list[cID].is_match!=3 && ovlp->list[cID].is_match!=4) continue;
            w_list_id = get_adjust_winid(&(ovlp->list[cID]), window_start, w_inf.window_length);
            if(w_list_id >= 0) break;///a matched window
        }
        if(i < (int64_t)du->length) continue;///if there is a matched window

        for (i = 0, mm = 0; i < (int64_t)du->length; i++) {///all windows are unmatched
            cID = (uint32_t)du->overlapID[i];
            if(ovlp->list[cID].is_match!=3 && ovlp->list[cID].is_match!=4) continue;
            ovlp->list[cID].is_match = 4; 
            ovlp->list[cID].align_length += window_end + 1 - window_start;
            mm++;
        }
        if(mm > 0) {
            kv_push(uint64_t, v_idx->a, (((uint64_t)window_start)<<32)|((uint64_t)window_end));
            v_idx->a.a[idx_i-1]++;
        }


        ///shorter than blockLen
        for (i = du->size-du->lengthNT, mLen = du->size-du->lengthNT, fc = 0; i < (int64_t)du->size; i++) {
            cID = (uint32_t)du->overlapID[i];
            if(ovlp->list[cID].is_match!=3 && ovlp->list[cID].is_match!=4) continue;
            get_win_se_by_normalize_xs(&(ovlp->list[cID]), window_start, blockLen, &ws, &we);
            w_list_id = get_adjust_winid(&(ovlp->list[cID]), ws, blockLen);
            if (w_list_id >= 0) {///matched
                cID = w_list_id; cID <<= 32; cID += (uint32_t)du->overlapID[i]; du->overlapID[i] = cID;
                if(mLen != i) {
                    mm = du->overlapID[i]; du->overlapID[i] = du->overlapID[mLen]; du->overlapID[mLen] = mm;
                }
                mLen++;
            } else {///unmatched
                cID = (uint32_t)-1; cID <<= 32; cID += (uint32_t)du->overlapID[i]; du->overlapID[i] = cID;
                fc++;
            }
        }
        // if(mLen == (int64_t)du->size) continue;///if all windows shorter than blockLen are matched
        if(fc == 0) continue;///no unmatched windows that are shorter than blockLen
        for (i = mLen; i < (int64_t)du->size; i++){///check the remaining unmatched windows that are shorter than blockLen
            cID = (uint32_t)du->overlapID[i];
            if(ovlp->list[cID].is_match!=3 && ovlp->list[cID].is_match!=4) continue;
            assert((du->overlapID[i]>>32)==(uint32_t)-1);
            get_win_se_by_normalize_xs(&(ovlp->list[cID]), window_start, blockLen, &ws, &we);
            for (k = du->size-du->lengthNT; k < mLen; k++) {///all matched windows
                fc = (uint32_t)du->overlapID[k]; fw = du->overlapID[k]>>32;
                assert(fw!=(uint32_t)-1); assert(ovlp->list[fc].is_match == 3 || ovlp->list[fc].is_match == 4);
                // if (ovlp->list[fc].w_list[fw].y_end == -1 || (ovlp->list[fc].is_match!=3 && ovlp->list[fc].is_match!=4)) fprintf(stderr, "ERROR\n");        
                ///if there is one matched window can cover the unmatched window
                if(ovlp->list[fc].w_list.a[fw].x_start<=ws && ovlp->list[fc].w_list.a[fw].x_end>=we) { 
                        break;
                }
            }

            if(k >= mLen) {///no matched window can cover the unmatched window
                ovlp->list[cID].is_match = 4; 
                ovlp->list[cID].align_length += we + 1 - ws;
                kv_push(uint64_t, v_idx->a, (((uint64_t)ws)<<32)|((uint64_t)we));
                v_idx->a.a[idx_i-1]++;
            }
        }
    }
}


inline void recalcate_window_advance(overlap_region_alloc* overlap_list, All_reads *rref, const ul_idx_t *uref, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read, kvec_t_u64_warp* v_idx, int64_t block_s, double e_rate, double e_rate_final)
{
    long long j, k, i;
    int threshold;
    long long y_id;
    int y_strand;
    long long y_readLen;
    long long x_start;
    long long x_end;
    long long x_len;
    long long total_y_start;
    long long total_y_end;
    long long y_start;
    long long Window_Len;
    char* x_string;
    char* y_string;
    int end_site;
    unsigned int error;
    int real_y_start;
    long long overlap_length;
    int extra_begin, extra_end;
    long long o_len;
    int64_t nw, a_nw, w_id, w_s, w_e, is_srt;
    double error_rate;
    uint64_t *w_idx;
    overlap_region *z;
    window_list *p = NULL;


    overlap_list->mapped_overlaps_length = 0;
    for (j = 0; j < (long long)overlap_list->length; j++) {
        z = &(overlap_list->list[j]); z->is_match = 0; is_srt = 1;
        if(z->w_list.n == 0) continue;///no alignment
        nw = get_num_wins(z->x_pos_s, z->x_pos_e+1, block_s); a_nw = z->w_list.n;
        kv_resize(uint64_t, v_idx->a, (uint64_t)nw); memset(v_idx->a.a, -1, sizeof((*v_idx->a.a))*nw); w_idx = v_idx->a.a;
        for (i = 0; i < a_nw; i++) {
            assert(z->w_list.a[i].y_end != -1);
            w_id = get_win_id_by_s(z, z->w_list.a[i].x_start, block_s, NULL);
            w_idx[w_id] = i;
        }
        // if(j == 248) {
        //     fprintf(stderr, "0-[M::%s] j::%lld, nw::%ld, a_nw::%ld, z->x_pos_s::%u, z->x_pos_e::%u, z->y_pos_s::%u, z->y_pos_e::%u, w_idx[0]::%lu\n", __func__, 
        //     j, nw, a_nw, z->x_pos_s, z->x_pos_e, z->y_pos_s, z->y_pos_e, w_idx[0]);
        // }

        y_id = z->y_id; y_strand = z->y_pos_strand; 
        y_readLen = (rref?(Get_READ_LENGTH((*rref), y_id)):(uref->ug->u.a[y_id].len));
        for (i = a_nw-1; i >= 0; i--) { //utilize the the end pos of pre-window in forward
            w_id = get_win_id_by_s(z, z->w_list.a[i].x_start, block_s, &w_e);
            // if(z->w_list.a[i].x_end != w_e) {
            //     fprintf(stderr, "[M::%s] block_s->%ld, w_id->%ld, z::x_pos_s->%u, z::x_pos_e->%u, x_start->%d, x_end->%d, w_e->%ld\n", __func__, block_s, w_id, z->x_pos_s, z->x_pos_e,
            //     z->w_list.a[i].x_start, z->w_list.a[i].x_end, w_e);
            // }
            assert(z->w_list.a[i].x_end == w_e);
            total_y_start = z->w_list.a[i].y_end + 1 - z->w_list.a[i].extra_begin;
            for (k = w_id + 1; k < nw; k++) {
                if(w_idx[k] != (uint64_t)-1) break;
                w_s = w_e + 1;
                w_id = get_win_id_by_s(z, w_s, block_s, &w_e);
                assert(w_id == k);
                extra_begin = extra_end = 0;
                if (total_y_start >= y_readLen) break;
                ///there is no problem for x
                x_start = w_s; x_end = w_e; x_len = x_end + 1 - x_start; y_start = total_y_start; 
                ///there are two potiential reasons for unmatched window:
                ///1. this window has a large number of differences
                ///2. DP does not start from the right offset
                if(rref) {
                    threshold = double_error_threshold(get_init_err_thres(x_len, e_rate, block_s, THRESHOLD), x_len);
                } else {
                    threshold = double_ul_error_threshold(get_init_err_thres(x_len, e_rate, block_s, THRESHOLD_MAX_SIZE), x_len);
                }
                
                Window_Len = x_len + (threshold << 1);

                if(!determine_overlap_region(threshold, y_start, y_id, Window_Len, (rref?(Get_READ_LENGTH((*rref), y_id)):(uref->ug->u.a[y_id].len)), 
                &extra_begin, &extra_end, &y_start, &o_len)) {
                    break;
                }
                if(o_len + threshold < x_len) break;
                
                if(rref) {
                    fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, rref, y_id, extra_begin, extra_end);
                } else {
                    fill_subregion_ul(dumy->overlap_region, y_start, o_len, y_strand, uref, y_id, extra_begin, extra_end);
                }

                x_string = g_read->seq + x_start; y_string = dumy->overlap_region;
                ///note!!! need notification
                end_site = Reserve_Banded_BPM(y_string, Window_Len, x_string, x_len, threshold, &error);                    
                if (error!=(unsigned int)-1) {///unmatched
                    kv_pushp(window_list, z->w_list, &p);
                    p->x_start = x_start;
                    p->x_end = x_end;
                    p->y_start = y_start;
                    p->y_end = y_start + end_site;
                    p->error = error;
                    p->extra_begin = extra_begin;
                    p->extra_end = extra_end;
                    p->error_threshold = threshold;
                    p->cidx = p->clen = 0;

                    z->align_length += x_len; w_idx[k] = z->w_list.n - 1;

                    if(is_srt && z->w_list.n > 1 && p->x_start < z->w_list.a[z->w_list.n-2].x_start) is_srt = 0;
                }
                else {
                    break;
                }

                total_y_start = y_start + end_site + 1 - extra_begin;
            }
        }
        // if(j == 248) {
        //     fprintf(stderr, "1-[M::%s] j::%lld, nw::%ld, a_nw::%ld, z->x_pos_s::%u, z->x_pos_e::%u, z->y_pos_s::%u, z->y_pos_e::%u, w_idx[0]::%lu\n", __func__, 
        //     j, nw, a_nw, z->x_pos_s, z->x_pos_e, z->y_pos_s, z->y_pos_e, w_idx[0]);
        // }
        for (i = 0; i < nw; i++) { //utilize the the start pos of next window in backward
            ///find the first matched window, which should not be the first window
            ///the pre-window of this matched window must be unmatched
            if(i > 0 && w_idx[i] != (uint64_t)-1 && w_idx[i-1] == (uint64_t)-1) {
                w_s = z->w_list.a[w_idx[i]].x_start;
                ///check if the start pos of this matched window has been calculated
                if(z->w_list.a[w_idx[i]].clen == 0) {
                    p = &(z->w_list.a[w_idx[i]]);
                    ///there is no problem for x
                    x_start = p->x_start; x_end = p->x_end; x_len = x_end + 1 - x_start; threshold = p->error_threshold;
                    /****************************may have bugs********************************/
                    ///should not adjust threshold, since this window can be matched by the old threshold
                    ///threshold = Adjust_Threshold(threshold, x_len);
                    /****************************may have bugs********************************/
                    Window_Len = x_len + (threshold << 1);
                    ///y_start is the real y_start
                    y_start = p->y_start; extra_begin = p->extra_begin; extra_end = p->extra_end;
                    o_len = Window_Len - extra_end - extra_begin;
                    if(rref) {
                        fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, rref, y_id, extra_begin, extra_end);
                    } else {
                        fill_subregion_ul(dumy->overlap_region, y_start, o_len, y_strand, uref, y_id, extra_begin, extra_end);
                    }
                    
                    x_string = g_read->seq + x_start;
                    y_string = dumy->overlap_region;

                    end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
                    &(dumy->path_length), dumy->matrix_bit, dumy->path, p->error, p->y_end - y_start);
                    assert(error != (unsigned int)-1);

                    {
                        ///this condition is always wrong
                        ///in best case, real_y_start = threshold, end_site = Window_Len - threshold - 1
                        if (end_site == Window_Len - 1 || real_y_start == 0) {
                            if(rref) {
                                if(fix_boundary(x_string, x_len, threshold, y_start, real_y_start, 
                                end_site, extra_begin, extra_end, y_id, Window_Len, rref, dumy, 
                                y_strand, error, &y_start, &real_y_start, &end_site, &extra_begin, 
                                &extra_end, &error)) {
                                    p->error = error; p->extra_begin = extra_begin; p->extra_end = extra_end;
                                }
                            } else {
                                if(fix_ul_boundary(x_string, x_len, threshold, y_start, real_y_start, 
                                end_site, extra_begin, extra_end, y_id, Window_Len, uref, dumy, 
                                y_strand, error, &y_start, &real_y_start, &end_site, &extra_begin, 
                                &extra_end, &error)) {
                                    p->error = error; p->extra_begin = extra_begin; p->extra_end = extra_end;
                                }
                            }
                        }
                                                 
                        generate_cigar(dumy->path, dumy->path_length, p, &(z->w_list), &real_y_start, &end_site, &error, x_string, x_len, y_string);   

                        ///note!!! need notification
                        real_y_start = y_start + real_y_start - extra_begin;
                        p->y_start = real_y_start; 
                        ///I forget why don't reduce the extra_begin for y_end
                        ///it seems extra_begin will be reduced at the end of this function 
                        p->y_end = y_start + end_site;
                        p->error = error;                         
                    }
                } else {
                    real_y_start = p->y_start;
                }

                ///the end pos for pre window is real_y_start - 1
                total_y_end = real_y_start - 1;
                ///find the unmatched window on the left of current matched window
                ///k starts from i - 1
                for (k = i - 1; k >= 0 && w_idx[k] == (uint64_t)-1; k--) {  
                    w_e = w_s - 1;
                    w_id = get_win_id_by_e(z, w_e, block_s, &w_s);
                    assert(w_id == k);
                    ///there is no problem in x
                    x_start = w_s; x_end = w_e; x_len = x_end + 1 - x_start;
                    ///there are two potiential reasons for unmatched window:
                    ///1. this window has a large number of differences
                    ///2. DP does not start from the right offset
                    if(rref) {
                        threshold = double_error_threshold(get_init_err_thres(x_len, e_rate, block_s, THRESHOLD), x_len);
                    } else {
                        threshold = double_ul_error_threshold(get_init_err_thres(x_len, e_rate, block_s, THRESHOLD_MAX_SIZE), x_len);
                    }
                    Window_Len = x_len + (threshold << 1);
                    if(total_y_end <= 0) break;
                    
                    ///y_start might be less than 0
                    y_start = total_y_end - x_len + 1;
                    if(!determine_overlap_region(threshold, y_start, y_id, Window_Len, (rref?(Get_READ_LENGTH((*rref), y_id)):(uref->ug->u.a[y_id].len)), 
                    &extra_begin, &extra_end, &y_start, &o_len)) {
                        break;
                    }

                    if(o_len + threshold < x_len) break;
                    
                    if(rref) {
                        fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, rref, y_id, extra_begin, extra_end);
                    } else {
                        fill_subregion_ul(dumy->overlap_region, y_start, o_len, y_strand, uref, y_id, extra_begin, extra_end);
                    }
                    x_string = g_read->seq + x_start;
                    y_string = dumy->overlap_region;

                    ///note!!! need notification
                    end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
                    &(dumy->path_length), dumy->matrix_bit, dumy->path, -1, -1);

                    if (error!=(unsigned int)-1) { 
                        ///this condition is always wrong
                        ///in best case, real_y_start = threshold, end_site = Window_Len - threshold - 1
                        if (end_site == Window_Len - 1 || real_y_start == 0) {
                            if(rref) {
                                fix_boundary(x_string, x_len, threshold, y_start, real_y_start, 
                                end_site, extra_begin, extra_end, y_id, Window_Len, rref, dumy, 
                                y_strand, error, &y_start, &real_y_start, &end_site, &extra_begin, 
                                &extra_end, &error);
                            } else {
                                fix_ul_boundary(x_string, x_len, threshold, y_start, real_y_start, 
                                end_site, extra_begin, extra_end, y_id, Window_Len, uref, dumy, 
                                y_strand, error, &y_start, &real_y_start, &end_site, &extra_begin, 
                                &extra_end, &error);
                            }
                        }

                        kv_pushp(window_list, z->w_list, &p);
                        p->x_start = x_start; p->x_end = x_end;///must set x_start/x_end here
                        generate_cigar(dumy->path, dumy->path_length, p, &(z->w_list), &real_y_start, &end_site, &error, x_string, x_len, y_string);  
                        
                        ///y_start has no shift, but y_end has shift               
                        p->y_start = y_start + real_y_start - extra_begin;
                        p->y_end = y_start + end_site;
                        p->error = error;
                        p->extra_begin = extra_begin;
                        p->extra_end = extra_end;
                        p->error_threshold = threshold;
                        z->align_length += x_len; w_idx[k] = z->w_list.n - 1;

                        if(is_srt && z->w_list.n > 1 && p->x_start < z->w_list.a[z->w_list.n-2].x_start) is_srt = 0;
                    }
                    else {
                        break;
                    }

                    total_y_end = y_start + real_y_start - 1 - extra_begin;
                }
            }
        }

        // if(j == 248) {
        //     fprintf(stderr, "2-[M::%s] j::%lld, nw::%ld, a_nw::%ld, z->x_pos_s::%u, z->x_pos_e::%u, z->y_pos_s::%u, z->y_pos_e::%u, w_idx[0]::%lu, w_idx[0]->cidx::%u, w_idx[0]->clen::%u, w_idx[0]->cigar[0]:%u\n", __func__, 
        //     j, nw, a_nw, z->x_pos_s, z->x_pos_e, z->y_pos_s, z->y_pos_e, w_idx[0], z->w_list.a[w_idx[0]].cidx, z->w_list.a[w_idx[0]].clen, z->w_list.c.a[z->w_list.a[w_idx[0]].cidx]);
        // }

        if(uref) {
            z->is_match = 0;
            if((((z->x_pos_e + 1 - z->x_pos_s)*MIN_UL_ALIN_RATE) <= z->align_length) && (z->align_length >= MIN_UL_ALIN_LEN)){
                z->is_match = 3; overlap_list->mapped_overlaps_length += z->align_length;
                ///sort for set_herror_win
                if(!is_srt) radix_sort_window_list_xs_srt(z->w_list.a, z->w_list.a + z->w_list.n);
            }
        }
    }

    // fprintf(stderr, "+++[M::%s::idx->%d::y_id->%u] z::align_length->%u, e_threshold->%f\n", 
    //         __func__, 27, overlap_list->list[27].y_id, overlap_list->list[27].align_length, e_rate);
    
    // fprintf(stderr, "+++[M::%s::idx->%d::y_id->%u] z::align_length->%u, e_threshold->%f\n", 
    //         __func__, 45, overlap_list->list[45].y_id, overlap_list->list[45].align_length, e_rate);

    // fprintf(stderr, "+++[M::%s::idx->%d::y_id->%u] z::align_length->%u, e_threshold->%f\n", 
    //         __func__, 277, overlap_list->list[277].y_id, overlap_list->list[277].align_length, e_rate);

    if(uref && overlap_list->mapped_overlaps_length > 0) {
        set_herror_win(overlap_list, dumy, v_idx, e_rate, g_read->length, block_s);
    }

    overlap_list->mapped_overlaps_length = 0;
    for (j = 0; j < (long long)overlap_list->length; j++) {
        z = &(overlap_list->list[j]);
        y_id = z->y_id; y_strand = z->y_pos_strand; 
        y_readLen = (rref?(Get_READ_LENGTH((*rref), y_id)):(uref->ug->u.a[y_id].len));
        overlap_length = z->x_pos_e + 1 - z->x_pos_s; //z->is_match = 0;
        // if(y_id == 0 || y_id == 1) {
        //     fprintf(stderr, "[M::%s::j->%lld] utg%.6dl(%c), align_length::%u, overlap_length::%lld\n", __func__, 
        //         j, (int32_t)z->y_id + 1, "+-"[z->y_pos_strand], z->align_length, overlap_length);
        // }
        // if(y_id == 4) {
        //     fprintf(stderr, "[M::%s::idx->%lld::] z::x_pos_s->%u, z::x_pos_e->%u, ovl->%lld, aln->%u\n", 
        //     __func__, j, z->x_pos_s, z->x_pos_e, overlap_length, z->align_length);
        // }

        ///debug_scan_cigar(&(overlap_list->list[j]));
        ///only calculate cigar for high quality overlaps
        if ((rref && (overlap_length*OVERLAP_THRESHOLD_HIFI_FILTER <= z->align_length)) || 
                                            (uref && (overlap_length*(1-e_rate) <= z->align_length))) {
            a_nw = z->w_list.n; 
            // int64_t tt = 0;
            for (i = 0, is_srt = 1; i < a_nw; i++) {
                p = &(z->w_list.a[i]);
                ///check if the cigar of this window has been got 
                if(p->clen == 0) {
                    ///there is no problem for x
                    x_start = p->x_start; x_end = p->x_end; x_len = x_end - x_start + 1;
                    /****************************may have bugs********************************/
                    ///threshold = x_len * asm_opt.max_ov_diff_ec;
                    threshold = p->error_threshold;
                    /****************************may have bugs********************************/
                    /****************************may have bugs********************************/
                    ///should not adjust threshold, since this window can be matched by the old threshold
                    ///threshold = Adjust_Threshold(threshold, x_len);
                    /****************************may have bugs********************************/
                    Window_Len = x_len + (threshold << 1);


                    ///y_start is the real y_start
                    ///for the window with cigar, y_start has already reduced extra_begin
                    y_start = p->y_start; extra_begin = p->extra_begin; extra_end = p->extra_end;
                    o_len = Window_Len - extra_end - extra_begin;
                    if(rref) {
                        fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, rref, y_id, extra_begin, extra_end);
                    } else {
                        fill_subregion_ul(dumy->overlap_region, y_start, o_len, y_strand, uref, y_id, extra_begin, extra_end);
                    }
                    x_string = g_read->seq + x_start; y_string = dumy->overlap_region;


                    ///note!!! need notification
                    end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
                    &(dumy->path_length), dumy->matrix_bit, dumy->path, p->error, p->y_end - y_start);
                    // if(!(error != (unsigned int)-1)) {
                    //     fprintf(stderr, "[M::%s::]\tqid::%u\tqlen::%lu\tq::[%d,\t%d)\ttid::%u\ttlen::%lu\tt::[%d,\t%d)\te_beg::%d\te_end::%d\terr::%d\n", __func__, 
                    //     overlap_list->list[j].x_id, Get_READ_LENGTH((*rref), overlap_list->list[j].x_id), 
                    //     p->x_start, p->x_end+1, 
                    //     overlap_list->list[j].y_id, Get_READ_LENGTH((*rref), overlap_list->list[j].y_id), 
                    //     p->y_start, p->y_end+1, 
                    //     p->extra_begin, p->extra_end, p->error);
                    //     fprintf(stderr, "[M::%s::]\tqcal_len::%lld\ttcal_len::%lld\tthres::%d\n", __func__, 
                    //     x_len, Window_Len, threshold);
                    //     fprintf(stderr, "qid::%u\nqname::%.*s\n\t%.*s\n", overlap_list->list[j].x_id, 
                    //     (int32_t)Get_NAME_LENGTH((*rref), overlap_list->list[j].x_id), 
                    //     Get_NAME((*rref), overlap_list->list[j].x_id), (int32_t)x_len, x_string);

                    //     fprintf(stderr, "tid::%u\ntname::%.*s\n\t%.*s\n", overlap_list->list[j].y_id, 
                    //     (int32_t)Get_NAME_LENGTH((*rref), overlap_list->list[j].y_id), 
                    //     Get_NAME((*rref), overlap_list->list[j].y_id), (int32_t)Window_Len, y_string);

                    // }
                    assert(error != (unsigned int)-1);

                    {
                        if (end_site == Window_Len - 1 || real_y_start == 0) {
                            if(rref) {
                                if(fix_boundary(x_string, x_len, threshold, y_start, real_y_start, end_site,
                                extra_begin, extra_end, y_id, Window_Len, rref, dumy, y_strand, error,
                                &y_start, &real_y_start, &end_site, &extra_begin, &extra_end, &error)) {
                                    p->error = error; p->extra_begin = extra_begin; p->extra_end = extra_end;
                                }
                            } else {
                                if(fix_ul_boundary(x_string, x_len, threshold, y_start, real_y_start, 
                                end_site, extra_begin, extra_end, y_id, Window_Len, uref, dumy, 
                                y_strand, error, &y_start, &real_y_start, &end_site, &extra_begin, 
                                &extra_end, &error)) {
                                    p->error = error; p->extra_begin = extra_begin; p->extra_end = extra_end;
                                }
                            }
                        }

                        generate_cigar(dumy->path, dumy->path_length, p, &(z->w_list), &real_y_start, &end_site, &error, x_string, x_len, y_string);    

                        ///note!!! need notification
                        real_y_start = y_start + real_y_start - extra_begin;
                        p->y_start = real_y_start;  
                        p->y_end = y_start + end_site - extra_begin;
                        p->error = error;                              
                    }

                    // if(y_id == 4) {
                    //     fprintf(stderr, "+[M::idx->%lld::] y_start->%d, y_end->%d, error->%d\n", 
                    //     j, p->y_start, p->y_end, p->error);
                    // }
                }
                else {
                    p->y_end -= p->extra_begin;
                    // if(y_id == 4) {
                    //     fprintf(stderr, "-[M::idx->%lld::] y_start->%d, y_end->%d, error->%d\n", 
                    //     j, p->y_start, p->y_end, p->error);
                    // }
                }
                // tt += p->error;
                if(is_srt && i > 0 && p->x_start < z->w_list.a[i-1].x_start) is_srt = 0;
            }

            if(!is_srt) radix_sort_window_list_xs_srt(z->w_list.a, z->w_list.a + z->w_list.n);
            error_rate = non_trim_error_rate(z, rref, uref, v_idx, dumy, g_read, e_rate, block_s);
            z->is_match = 0;
            // if(y_id == 4) {
            //     fprintf(stderr, "[M::%s::idx->%lld::] z::x_pos_s->%u, z::x_pos_e->%u, ovl->%lld, aln->%u, error_rate->%f, e_rate_final->%f, tt->%ld\n", 
            //     __func__, j, z->x_pos_s, z->x_pos_e, overlap_length, z->align_length, error_rate, e_rate_final, tt);
            // }

            if (error_rate <= e_rate_final/**asm_opt.max_ov_diff_final**/) {
                overlap_list->mapped_overlaps_length += overlap_length;
                z->is_match = 1; append_unmatched_wins(z, block_s);
                // if(j == 248) {
                //     fprintf(stderr, "3-[M::%s] j::%lld, nw::%ld, a_nw::%ld, z->x_pos_s::%u, z->x_pos_e::%u, z->y_pos_s::%u, z->y_pos_e::%u, w_idx[0]::%lu, w_idx[0]->cidx::%u, w_idx[0]->clen::%u, w_idx[0]->cigar[0]:%u\n", __func__, 
                //     j, nw, a_nw, z->x_pos_s, z->x_pos_e, z->y_pos_s, z->y_pos_e, w_idx[0], z->w_list.a[w_idx[0]].cidx, z->w_list.a[w_idx[0]].clen, z->w_list.c.a[z->w_list.a[w_idx[0]].cidx]);
                // }
                if(rref) {
                    calculate_boundary_cigars(z, rref, dumy, g_read, e_rate);
                } else {
                    calculate_ul_boundary_cigars(z, uref, dumy, g_read, e_rate, block_s);
                }
                // if(j == 248) {
                //     fprintf(stderr, "4-[M::%s] j::%lld, nw::%ld, a_nw::%ld, z->x_pos_s::%u, z->x_pos_e::%u, z->y_pos_s::%u, z->y_pos_e::%u, w_idx[0]::%lu, w_idx[0]->cidx::%u, w_idx[0]->clen::%u, w_idx[0]->cigar[0]:%u\n", __func__, 
                //     j, nw, a_nw, z->x_pos_s, z->x_pos_e, z->y_pos_s, z->y_pos_e, w_idx[0], z->w_list.a[w_idx[0]].cidx, z->w_list.a[w_idx[0]].clen, z->w_list.c.a[z->w_list.a[w_idx[0]].cidx]);
                // }
                // if((int64_t)z->x_pos_s!=z->w_list.a[0].x_start || 
                //         (int64_t)z->x_pos_e!=z->w_list.a[z->w_list.n-1].x_end) {
                //     fprintf(stderr, "[M::%s] z::x_pos_s->%u, z::x_pos_e->%u, (0)::x_start->%d, (wn-1)x_end->%d, z->w_list.n->%ld\n", __func__,
                //      z->x_pos_s, z->x_pos_e, z->w_list.a[0].x_start, z->w_list.a[z->w_list.n-1].x_end, (int64_t)z->w_list.n);
                // }

                // assert(get_num_wins(z->x_pos_s, z->x_pos_e+1, block_s)==(int64_t)z->w_list.n);
                // assert((int64_t)z->x_pos_s==z->w_list.a[0].x_start && 
                //                             (int64_t)z->x_pos_e==z->w_list.a[z->w_list.n-1].x_end);
            } else if (error_rate <= /**asm_opt.max_ov_diff_final**/e_rate_final * 1.5) {
                z->is_match = 3;
            }
            // fprintf(stderr, "[M::%s::idx->%lld::is_match->%u] z::y_id->%u, z::x_pos_s->%u, z::x_pos_e->%u, error_rate->%f, e_threshold->%f\n", 
            // __func__, j, z->is_match, z->y_id,  z->x_pos_s, z->x_pos_e, error_rate, e_rate);
        } else {///it impossible to be matched
            z->is_match = 0;
            // fprintf(stderr, "[M::%s::idx->%ld::is_match->%u] z::x_pos_s->%u, z::x_pos_e->%u, error_rate->-1, e_threshold->%f\n", 
            // __func__, j, z->is_match, z->x_pos_s, z->x_pos_e, e_rate);
        }

        
    }
    ///debug_window_cigar(overlap_list, g_read, dumy, rref, 1, 1);
}

uint32_t inline simi_pass(int64_t ol, int64_t aln_ol, uint32_t second_ck, double o_rate, double *e_rate)
{
    if(aln_ol == 0 || ol == 0) return 0;
    if((!second_ck) && (!e_rate)) {
        // if((ol*OVERLAP_THRESHOLD_FILTER) <= aln_ol) return 1;
        if((ol*o_rate) <= aln_ol) return 1;
    } else if(e_rate) {
        if((ol*((double)(((double)1.0)-(*e_rate)))) <= aln_ol) return 1;
    } else if(second_ck) {
        if(((ol*MIN_UL_ALIN_RATE) <= aln_ol) && (aln_ol >= MIN_UL_ALIN_LEN)) return 1;
    }

    // if(rref) {
    //     if((ol*OVERLAP_THRESHOLD_FILTER) <= aln_ol) return 1;
    // } else if(uref) {
    //     if(e_rate) {
    //         if((ol*((double)(((double)1.0)-(*e_rate)))) <= aln_ol) return 1;
    //     } else {
    //         if(((ol*MIN_UL_ALIN_RATE) <= aln_ol) && (aln_ol >= MIN_UL_ALIN_LEN)) return 1;
    //     }
    // }

    return 0;
}

inline uint32_t gen_backtrace(window_list *p, overlap_region *z, All_reads *rref, const ul_idx_t *uref, UC_Read* g_read, Correct_dumy* dumy,
int32_t y_strand, int32_t y_id)
{
    int64_t x_start, x_end, x_len, Window_Len, o_len; 
    int32_t threshold; long long y_start; 
    int real_y_start = 0, end_site, extra_begin, extra_end;
    char *x_string, *y_string; unsigned int error;
    ///there is no problem for x
    x_start = p->x_start; x_end = p->x_end; x_len = x_end - x_start + 1;
    threshold = p->error_threshold; Window_Len = x_len + (threshold << 1);


    ///y_start is the real y_start
    ///for the window with cigar, y_start has already reduced extra_begin
    y_start = p->y_start; extra_begin = p->extra_begin; extra_end = p->extra_end;
    o_len = Window_Len - extra_end - extra_begin;
    if(rref) {
        fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, rref, y_id, extra_begin, extra_end);
    } else {
        fill_subregion_ul(dumy->overlap_region, y_start, o_len, y_strand, uref, y_id, extra_begin, extra_end);
    }
    x_string = g_read->seq + x_start; y_string = dumy->overlap_region;


    ///note!!! need notification
    end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
    &(dumy->path_length), dumy->matrix_bit, dumy->path, p->error, p->y_end - y_start);
    // assert(error != (unsigned int)-1);
    if(error != (unsigned int)-1) {
        ///this condition is always wrong
        ///in best case, real_y_start = threshold, end_site = Window_Len - threshold - 1
        if (end_site == Window_Len - 1 || real_y_start == 0) {
            if(rref) {
                if(fix_boundary(x_string, x_len, threshold, y_start, real_y_start, end_site,
                extra_begin, extra_end, y_id, Window_Len, rref, dumy, y_strand, error,
                &y_start, &real_y_start, &end_site, &extra_begin, &extra_end, &error)) {
                    p->error = error; p->extra_begin = extra_begin; p->extra_end = extra_end;
                }
            } else {
                if(fix_ul_boundary(x_string, x_len, threshold, y_start, real_y_start, 
                end_site, extra_begin, extra_end, y_id, Window_Len, uref, dumy, 
                y_strand, error, &y_start, &real_y_start, &end_site, &extra_begin, 
                &extra_end, &error)) {
                    p->error = error; p->extra_begin = extra_begin; p->extra_end = extra_end;
                }
            }
        }

        generate_cigar(dumy->path, dumy->path_length, p, &(z->w_list), &real_y_start, &end_site, &error, x_string, x_len, y_string);    

        ///note!!! need notification
        real_y_start = y_start + real_y_start - extra_begin;
        p->y_start = real_y_start;  
        p->y_end = y_start + end_site - extra_begin;
        p->error = error;         
        return 1;                     
    }
    p->error = -1;
    return 0;
}

void gen_rev_str(char *in, char **out, uint32_t len)
{
    char *r; uint32_t k; MALLOC(r, len); (*out) = r;
    for (k = 0; k < len; k++) r[k] = in[len-k-1];
}

inline uint32_t gen_backtrace_adv(window_list *p, overlap_region *z, All_reads *rref, hpc_t *hpc_g, const ul_idx_t *uref, 
char *qstr, char *tstr, char *tstr1, Correct_dumy* dumy, uint32_t rev, uint32_t id)
{
    int64_t qs, qe, ql, aln_l, t_pri_l, thres, ts; 
    int r_ts = 0, t_end; int64_t aux_beg, aux_end;
    char *q_string, *t_string; unsigned int error;
    ///there is no problem for x
    qs = p->x_start; qe = p->x_end; ql = qe + 1 - qs;
    thres = p->error_threshold; aln_l = ql + (thres<<1);

    ///y_start is the real y_start
    ///for the window with cigar, y_start has already reduced extra_begin
    ts = p->y_start; aux_beg = p->extra_begin; aux_end = p->extra_end;
    t_pri_l = aln_l - aux_beg - aux_end;

    q_string = qstr + qs; 
    if(rref) {
        fill_subregion(tstr, ts, t_pri_l, rev, rref, id, aux_beg, aux_end); t_string = tstr;
    } else {
        t_string = return_str_seq(tstr, ts, t_pri_l, rev, hpc_g, uref, id, aux_beg, aux_end);
    }
    
    t_end = Reserve_Banded_BPM_PATH(t_string, aln_l, q_string, ql, thres, &error, &r_ts,
        &(dumy->path_length), dumy->matrix_bit, dumy->path, p->error, p->y_end - ts);

    // assert(error != (unsigned int)-1);
    if(error != (unsigned int)-1) {
        /**
        bit_extz_t exz, exz64; init_bit_extz_t(&exz, thres); init_bit_extz_t(&exz64, thres);
        // ed_band_cal_extension_64_0_w(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz64);
        // ed_band_cal_extension_64_0_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz64);
        // cigar_check(t_string+r_ts, q_string, &(exz64));
        // exz64.err = INT32_MAX;
        // ed_band_cal_extension_64_0_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz64);
        // cigar_check(t_string+r_ts, q_string, &(exz64));

        // ed_band_cal_extension_infi_0_w(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, NULL, &exz);
        // ed_band_cal_extension_infi_0_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, NULL, &exz);
        // cigar_check(t_string+r_ts, q_string, &exz);
        // exz.err = INT32_MAX;
        // ed_band_cal_extension_infi_0_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, NULL, &exz);
        // cigar_check(t_string+r_ts, q_string, &exz);
        // assert(exz.err <= (int32_t)error && exz.err >= 0); 
        // assert(exz.err == exz64.err && exz.ps == exz64.ps && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);
        
        // ed_band_cal_extension_256_0_w(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz);
        // ed_band_cal_extension_256_0_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz);
        // cigar_check(t_string+r_ts, q_string, &exz);
        // exz.err = INT32_MAX;
        // ed_band_cal_extension_256_0_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz);
        // cigar_check(t_string+r_ts, q_string, &exz);
        // assert(exz.err <= (int32_t)error && exz.err >= 0); 
        // assert(exz.err == exz64.err && exz.ps == exz64.ps && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);

        // char *qr, *tr;
        // gen_rev_str(q_string, &qr, ql); gen_rev_str(t_string+r_ts, &tr, t_end+1-r_ts); 

        // ed_band_cal_extension_64_1_w(tr, t_end+1-r_ts, qr, ql, thres, &exz64);
        // ed_band_cal_extension_64_1_w_trace(tr, t_end+1-r_ts, qr, ql, thres, &exz64);
        // cigar_check(tr, qr, &exz64);
        // exz64.err = INT32_MAX;
        // ed_band_cal_extension_64_1_w_trace(tr, t_end+1-r_ts, qr, ql, thres, &exz64);
        // cigar_check(tr, qr, &exz64);
        // assert(exz.err == exz64.err && exz.ps == (exz64.pl-exz64.pe-1) && exz.pe == (exz64.pl-exz64.ps-1) && exz.ts == (exz64.tl-exz64.te-1) && exz.te == (exz64.tl-exz64.ts-1));

        // ed_band_cal_extension_infi_1_w(tr, t_end+1-r_ts, qr, ql, thres, NULL, &exz);
        // ed_band_cal_extension_infi_1_w_trace(tr, t_end+1-r_ts, qr, ql, thres, NULL, &exz);
        // cigar_check(tr, qr, &exz);
        // exz.err = INT32_MAX;
        // ed_band_cal_extension_infi_1_w_trace(tr, t_end+1-r_ts, qr, ql, thres, NULL, &exz);
        // cigar_check(tr, qr, &exz);
        // assert(exz.err == exz64.err && exz.ps == exz64.ps && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);

        // ed_band_cal_extension_256_1_w(tr, t_end+1-r_ts, qr, ql, thres, &exz);
        // ed_band_cal_extension_256_1_w_trace(tr, t_end+1-r_ts, qr, ql, thres, &exz);
        // cigar_check(tr, qr, &exz);
        // exz.err = INT32_MAX;
        // ed_band_cal_extension_256_1_w_trace(tr, t_end+1-r_ts, qr, ql, thres, &exz);
        // cigar_check(tr, qr, &exz);
        // assert(exz.err == exz64.err && exz.ps == exz64.ps && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);

        // free(qr); free(tr);


        // ed_band_cal_global_64_w(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz64);
        // ed_band_cal_global_64_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz64);
        // cigar_check(t_string+r_ts, q_string, &(exz64));
        // exz64.err = INT32_MAX;
        // ed_band_cal_global_64_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz64);
        // cigar_check(t_string+r_ts, q_string, &(exz64));

        // ed_band_cal_global_infi_w(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, NULL, &exz);
        // ed_band_cal_global_infi_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, NULL, &exz);
        // cigar_check(t_string+r_ts, q_string, &exz);
        // exz.err = INT32_MAX;
        // ed_band_cal_global_infi_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, NULL, &exz);
        // cigar_check(t_string+r_ts, q_string, &exz);
        // assert(exz.err <= (int32_t)error && exz.err >= 0); 
        // assert(exz.err == exz64.err && exz.ps == exz64.ps && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);

        // ed_band_cal_global_256_w(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz);
        // ed_band_cal_global_256_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz);
        // cigar_check(t_string+r_ts, q_string, &exz);
        // exz.err = INT32_MAX;
        // ed_band_cal_global_256_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz);
        // cigar_check(t_string+r_ts, q_string, &exz);
        // assert(exz.err <= (int32_t)error && exz.err >= 0); 
        // assert(exz.err == exz64.err && exz.ps == exz64.ps && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);

        
        // ed_band_cal_semi_64_w(t_string, aln_l, q_string, ql, thres, &exz64);
        // ed_band_cal_semi_64_w_trace(t_string, aln_l, q_string, ql, thres, &exz64);
        // cigar_check(t_string, q_string, &(exz64));
        // exz64.err = INT32_MAX;
        // ed_band_cal_semi_64_w_trace(t_string, aln_l, q_string, ql, thres, &exz64);
        // cigar_check(t_string, q_string, &(exz64));

        // ed_band_cal_semi_infi_w(t_string, aln_l, q_string, ql, thres, NULL, &exz);
        // ed_band_cal_semi_infi_w_trace(t_string, aln_l, q_string, ql, thres, NULL, &exz);
        // cigar_check(t_string, q_string, &exz);
        // exz.err = INT32_MAX;
        // ed_band_cal_semi_infi_w_trace(t_string, aln_l, q_string, ql, thres, NULL, &exz);
        // cigar_check(t_string, q_string, &exz);
        // assert(exz.err <= (int32_t)error && exz.err >= 0); 
        // assert(exz.err == exz64.err && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);

        // ed_band_cal_semi_256_w(t_string, aln_l, q_string, ql, thres, &exz);
        // ed_band_cal_semi_256_w_trace(t_string, aln_l, q_string, ql, thres, &exz);
        // cigar_check(t_string, q_string, &exz);
        // exz.err = INT32_MAX;
        // ed_band_cal_semi_256_w_trace(t_string, aln_l, q_string, ql, thres, &exz);
        // cigar_check(t_string, q_string, &exz);
        // assert(exz.err <= (int32_t)error && exz.err >= 0); 
        // assert(exz.err == exz64.err && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);



        ed_band_cal_semi_64_w_absent_diag(t_string+aux_beg, aln_l-aux_beg-aux_end, q_string, ql, 
                                                                                thres, aux_beg, &exz64);
        // ed_band_cal_semi_64_w_absent_diag_trace(t_string+aux_beg, aln_l-aux_beg-aux_end, q_string, ql, 
        //                                                                         thres, aux_beg, &exz64);
        // cigar_check(t_string+aux_beg, q_string, &(exz64));
        // exz64.err = INT32_MAX;
        // ed_band_cal_semi_64_w_absent_diag_trace(t_string+aux_beg, aln_l-aux_beg-aux_end, q_string, ql, 
        //                                                                         thres, aux_beg, &exz64);
        // cigar_check(t_string+aux_beg, q_string, &(exz64));



        ed_band_cal_semi_infi_w_absent_diag(t_string+aux_beg, aln_l-aux_beg-aux_end, q_string, ql, 
                                                                                thres, aux_beg, NULL, &exz);
        ed_band_cal_semi_infi_w_absent_diag_trace(t_string+aux_beg, aln_l-aux_beg-aux_end, q_string, ql, 
                                                                                thres, aux_beg, NULL, &exz);
        cigar_check(t_string+aux_beg, q_string, &(exz));
        exz.err = INT32_MAX;
        ed_band_cal_semi_infi_w_absent_diag_trace(t_string+aux_beg, aln_l-aux_beg-aux_end, q_string, ql, 
                                                                                thres, aux_beg, NULL, &exz);
        cigar_check(t_string+aux_beg, q_string, &(exz));
        assert(exz.err <= (int32_t)error && exz.err >= 0); 
        assert(exz.err == exz64.err && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);
        exz64.err = exz.err; exz64.ps = exz.ps; exz64.pe = exz.pe; exz64.ts = exz.ts; exz64.te = exz.te;


        ed_band_cal_semi_256_w_absent_diag(t_string+aux_beg, aln_l-aux_beg-aux_end, q_string, ql, 
                                                                                thres, aux_beg, &exz);
        ed_band_cal_semi_256_w_absent_diag_trace(t_string+aux_beg, aln_l-aux_beg-aux_end, q_string, ql, 
                                                                                thres, aux_beg, &exz);
        cigar_check(t_string+aux_beg, q_string, &(exz));
        exz.err = INT32_MAX;
        ed_band_cal_semi_256_w_absent_diag_trace(t_string+aux_beg, aln_l-aux_beg-aux_end, q_string, ql, 
                                                                                thres, aux_beg, &exz);
        cigar_check(t_string+aux_beg, q_string, &(exz));
        assert(exz.err <= (int32_t)error && exz.err >= 0); 
        assert(exz.err == exz64.err && exz.ps == exz64.ps && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);
        // if((!(exz.err == exz64.err && exz.ps == exz64.ps && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te)) || (exz.err != (int32_t)error)) {
        //     fprintf(stderr, "\n[M::%s::semi] error::%u, ql::%ld, thres::%ld, exz.err::%d, exz64.err::%d, exz.ps::%d, exz64.ps::%d, exz.pe::%d, exz64.pe::%d, exz.ts::%d, exz64.ts::%d, exz.te::%d, exz64.te::%d\n", __func__, 
        //     error, ql, thres, exz.err, exz64.err, exz.ps, exz64.ps, exz.pe, exz64.pe, exz.ts, exz64.ts, exz.te, exz64.te);
        //     fprintf(stderr, "[tstr] %.*s\n", (int32_t)aln_l, t_string);
        //     fprintf(stderr, "[qstr] %.*s\n", (int32_t)ql, q_string);
        // }
        destroy_bit_extz_t(&exz); destroy_bit_extz_t(&exz64);

        // if(exz.err > (int32_t)error && ql == 1) {
        //     fprintf(stderr, "[M::%s::] error::%u, ed_extension::%d, ql::%ld, thres::%ld\n", 
        //     __func__, error, exz.err, ql, thres);
        //     fprintf(stderr, "[tstr] %.*s\n", t_end+1-r_ts, t_string+r_ts);
        //     fprintf(stderr, "[qstr] %.*s\n", (int32_t)ql, q_string);
        // }
        
        // assert(ed_band_cal_global(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres) == 
        //                 ed_band_cal_global_128bit(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres));
        **/
        
        ///this condition is always wrong
        ///in best case, r_ts = threshold, t_end = aln_l - thres - 1
        if (((t_end+1) == aln_l) || (r_ts == 0)) {
            if(recal_boundary(q_string, tstr1, ql, thres, ts, r_ts, t_end,
            aux_beg, aux_end, error, id, aln_l, rev, dumy, rref, hpc_g, uref, 
            &ts, &r_ts, &t_end, &aux_beg, &aux_end, &error)) {
                p->error = error; p->extra_begin = aux_beg; p->extra_end = aux_end;
                t_string = update_des_str(tstr, ts, aln_l-aux_beg-aux_end, rev, rref, hpc_g, uref,
                                                    id, aux_beg, aux_end, hpc_g?NULL:tstr1);
            }
        }

        generate_cigar(dumy->path, dumy->path_length, p, &(z->w_list), &r_ts, &t_end, &error, q_string, ql, t_string);    

        p->y_start = ts + r_ts - aux_beg; 
        p->y_end = ts + t_end - aux_beg;
        p->error = error;         
        return 1;                     
    }
    p->error = -1;
    return 0;
}

inline uint32_t aln_wlst_adv(overlap_region *z, All_reads *rref, hpc_t *hpc_g, 
const ul_idx_t *uref, char *qstr, char *tstr, char *tstr1, Correct_dumy* dumy,
uint32_t rev, uint32_t id, int64_t qs, int64_t qe, int64_t t_s, int64_t block_s, 
double e_rate, uint32_t is_cigar)
{
    int64_t ql, aln_l, t_tot_l; window_list *p = NULL; int r_ts = 0, t_end;
    int64_t aux_beg, aux_end, t_pri_l;
    int64_t thres; char *q_string, *t_string; unsigned int error;
    ql = qe + 1 - qs;
    ///there are two potiential reasons for unmatched window:
    ///1. this window has a large number of differences
    ///2. DP does not start from the right offset
    if(rref) {
        thres = double_error_threshold(get_init_err_thres(ql, e_rate, block_s, THRESHOLD), ql);
    } else {
        thres = double_ul_error_threshold(get_init_err_thres(ql, e_rate, block_s, THRESHOLD_MAX_SIZE), ql);
    }
    aln_l = ql + (thres << 1);
    if(hpc_g) t_tot_l = hpc_len(*hpc_g, id);
    else if(uref) t_tot_l = uref->ug->u.a[id].len;
    else t_tot_l = Get_READ_LENGTH((*rref), id);

    if(!init_waln(thres, t_s, t_tot_l, aln_l, &aux_beg, &aux_end, &t_s, &t_pri_l)) return 0;
    if(t_pri_l + thres < ql) return 0;

    q_string = qstr + qs; 
    if(rref) {
        fill_subregion(tstr, t_s, t_pri_l, rev, rref, id, aux_beg, aux_end); t_string = tstr;
    } else {
        t_string = return_str_seq(tstr, t_s, t_pri_l, rev, hpc_g, uref, id, aux_beg, aux_end);
    }

    if(is_cigar) {
        ///note!!! need notification
        t_end = Reserve_Banded_BPM_PATH(t_string, aln_l, q_string, ql, thres, &error, &r_ts,
        &(dumy->path_length), dumy->matrix_bit, dumy->path, -1, -1);
    } else {
        ///note!!! need notification
        t_end = Reserve_Banded_BPM(t_string, aln_l, q_string, ql, thres, &error);    
    }
    if(error!=(unsigned int)-1) {
        if(is_cigar) {
            ///this condition is always wrong
            ///in best case, r_ts = threshold, t_end = aln_l - thres - 1
            if (((t_end+1) == aln_l) || (r_ts == 0)) {
                if(recal_boundary(q_string, tstr1, ql, thres, t_s, r_ts, t_end,
                aux_beg, aux_end, error, id, aln_l, rev, dumy, rref, hpc_g, uref, 
                &t_s, &r_ts, &t_end, &aux_beg, &aux_end, &error)) {
                    t_string = update_des_str(tstr, t_s, aln_l-aux_beg-aux_end, rev, rref, hpc_g, uref,
                                                    id, aux_beg, aux_end, hpc_g?NULL:tstr1);
                }
            }
        }

        kv_pushp(window_list, z->w_list, &p);
        p->x_start = qs; p->x_end = qe; ///must set x_start/x_end here
        if(is_cigar) {
            generate_cigar(dumy->path, dumy->path_length, p, &(z->w_list), &r_ts, &t_end, &error, q_string, ql, t_string);  
        } else {
            p->cidx = p->clen = 0;
        }
        p->y_start = t_s + r_ts;///difference
        p->y_end = t_s + t_end;
        p->error = error;
        p->extra_begin = aux_beg;
        p->extra_end = aux_end;
        p->error_threshold = thres;
        z->align_length += ql; 
        
        return 1;
    }
    return 0;
}

void push_wcigar(window_list *idx, window_list_alloc *res, bit_extz_t *exz)
{
    idx->cidx = res->c.n; idx->clen = exz->cigar.n; res->c.n += exz->cigar.n;
    kv_resize(uint16_t, res->c, res->c.n);
    memcpy(res->c.a+idx->cidx, exz->cigar.a, exz->cigar.n*sizeof(*(res->c.a)));
}

inline uint32_t aln_wlst_adv_exz(overlap_region *z, All_reads *rref, hpc_t *hpc_g, 
const ul_idx_t *uref, char *qstr, char *tstr, bit_extz_t *exz,
uint32_t rev, uint32_t id, int64_t qs, int64_t qe, int64_t t_s, int64_t block_s, 
double e_rate, uint32_t is_cigar)
{
    int64_t ql, tl, aln_l, t_tot_l; window_list *p = NULL; ///int r_ts = 0, t_end;
    int64_t aux_beg, aux_end, t_pri_l; int64_t thres; char *q_string, *t_string; 
    ql = qe + 1 - qs;
    ///there are two potiential reasons for unmatched window:
    ///1. this window has a large number of differences
    ///2. DP does not start from the right offset
    if(rref) {
        thres = double_error_threshold(get_init_err_thres(ql, e_rate, block_s, THRESHOLD), ql);
    } else {
        thres = double_ul_error_threshold(get_init_err_thres(ql, e_rate, block_s, THRESHOLD_MAX_SIZE), ql);
    }
    aln_l = ql + (thres << 1);
    if(hpc_g) t_tot_l = hpc_len(*hpc_g, id);
    else if(uref) t_tot_l = uref->ug->u.a[id].len;
    else t_tot_l = Get_READ_LENGTH((*rref), id);

    if(!init_waln(thres, t_s, t_tot_l, aln_l, &aux_beg, &aux_end, &t_s, &t_pri_l)) return 0;
    if(t_pri_l + thres < ql) return 0;

    q_string = qstr + qs; 
    if(rref) {
        recover_UC_Read_sub_region(tstr, t_s, t_pri_l, rev, rref, id); t_string = tstr;
    } else {
        t_string = return_str_seq_exz(tstr, t_s, t_pri_l, rev, hpc_g, uref, id);
    }
    tl = t_pri_l;
    if(is_cigar) {
        clear_align(*exz);
        ed_band_cal_semi_64_w_absent_diag_trace(t_string, tl, q_string, ql, thres, aux_beg, exz);
    } else {
        ed_band_cal_semi_64_w_absent_diag(t_string, tl, q_string, ql, thres, aux_beg, exz); exz->ps = 0;
    }

    // if(id == 40 && qs == 79670 && qe == 79824) {
    //     fprintf(stderr, "\n[M::%s::semi] exz->ps::%d, exz->pe::%d, exz->ts::%d, exz->te::%d, exz->err::%d, exz->cigar.n::%d\n", 
    //     __func__, exz->ps, exz->pe, exz->ts, exz->te, exz->err, (int32_t)exz->cigar.n);
    // }

    if(is_align(*exz)) {
        kv_pushp(window_list, z->w_list, &p);
        p->x_start = qs; p->x_end = qe; ///must set x_start/x_end here
        p->y_start = t_s + exz->ps;///difference
        p->y_end = t_s + exz->pe;
        p->error = exz->err;
        p->cidx = p->clen = 0;
        if(is_cigar) {
            push_wcigar(p, &(z->w_list), exz);
            ///this condition is always wrong
            ///in best case, r_ts = threshold, t_end = aln_l - thres - 1
            if ((((exz->pe+1) == tl) || (exz->ps == 0)) && (exz->err > 0)) {
                if(recal_boundary_exz(q_string, tstr, ql, tl, thres, t_s, exz->ps, exz->pe,
                exz->err, id, rev, exz, rref, hpc_g, uref, &t_s, &aux_beg, &aux_end)) {
                    //update cigar
                    z->w_list.c.n = p->cidx; push_wcigar(p, &(z->w_list), exz);

                    p->y_start = t_s + exz->ps;///difference
                    p->y_end = t_s + exz->pe;
                    p->error = exz->err;
                }
            }
        } 

        p->extra_begin = aux_beg;
        p->extra_end = aux_end;
        p->error_threshold = thres;
        z->align_length += ql; 
        return 1;
    }
    return 0;
}

inline uint32_t aln_wlst(overlap_region *z, All_reads *rref, const ul_idx_t *uref, UC_Read* g_read, Correct_dumy* dumy,
int32_t y_strand, int32_t y_id, int64_t x_start, int64_t x_end, long long y_start, int64_t block_s, double e_rate, int32_t is_cigar)
{
    int64_t x_len, Window_Len; window_list *p = NULL; long long o_len;
    int32_t threshold; int real_y_start = 0, end_site, extra_begin, extra_end;
    char *x_string, *y_string; unsigned int error;
    x_len = x_end + 1 - x_start;
    ///there are two potiential reasons for unmatched window:
    ///1. this window has a large number of differences
    ///2. DP does not start from the right offset
    if(rref) {
        threshold = double_error_threshold(get_init_err_thres(x_len, e_rate, block_s, THRESHOLD), x_len);
    } else {
        threshold = double_ul_error_threshold(get_init_err_thres(x_len, e_rate, block_s, THRESHOLD_MAX_SIZE), x_len);
    }
    Window_Len = x_len + (threshold << 1);
    ///y_start might be less than 0
    if(!determine_overlap_region(threshold, y_start, y_id, Window_Len, (rref?(Get_READ_LENGTH((*rref), y_id)):(uref->ug->u.a[y_id].len)), 
    &extra_begin, &extra_end, &y_start, &o_len)) {
        return 0;
    }

    if(o_len + threshold < x_len) return 0;
    if(rref) {
        fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, rref, y_id, extra_begin, extra_end);
    } else {
        fill_subregion_ul(dumy->overlap_region, y_start, o_len, y_strand, uref, y_id, extra_begin, extra_end);
    }
    x_string = g_read->seq + x_start;
    y_string = dumy->overlap_region;

    if(is_cigar) {
        ///note!!! need notification
        end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
        &(dumy->path_length), dumy->matrix_bit, dumy->path, -1, -1);
    } else {
        ///note!!! need notification
        end_site = Reserve_Banded_BPM(y_string, Window_Len, x_string, x_len, threshold, &error);    
    }
    if(error!=(unsigned int)-1) {
        if(is_cigar) {
            ///this condition is always wrong
            ///in best case, real_y_start = threshold, end_site = Window_Len - threshold - 1
            if (end_site == Window_Len - 1 || real_y_start == 0) {
                if(rref) {
                    fix_boundary(x_string, x_len, threshold, y_start, real_y_start, 
                    end_site, extra_begin, extra_end, y_id, Window_Len, rref, dumy, 
                    y_strand, error, &y_start, &real_y_start, &end_site, &extra_begin, 
                    &extra_end, &error);
                } else {
                    fix_ul_boundary(x_string, x_len, threshold, y_start, real_y_start, 
                    end_site, extra_begin, extra_end, y_id, Window_Len, uref, dumy, 
                    y_strand, error, &y_start, &real_y_start, &end_site, &extra_begin, 
                    &extra_end, &error);
                }
            }
        }

        kv_pushp(window_list, z->w_list, &p);
        p->x_start = x_start; p->x_end = x_end; ///must set x_start/x_end here
        if(is_cigar) {
            generate_cigar(dumy->path, dumy->path_length, p, &(z->w_list), &real_y_start, &end_site, &error, x_string, x_len, y_string);  
        } else {
            p->cidx = p->clen = 0;
        }
        p->y_start = y_start + real_y_start;///difference
        p->y_end = y_start + end_site;
        p->error = error;
        p->extra_begin = extra_begin;
        p->extra_end = extra_end;
        p->error_threshold = threshold;
        z->align_length += x_len; 
        
        return 1;
    }
    return 0;
}
uint64_t realign_ed(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, 
char *tstr, char *tstr_1, Correct_dumy* dumy, kvec_t_u64_warp* v_idx, int64_t block_s, double e_rate, 
double *e_rate_final, uint32_t sec_check, int64_t *is_sort);
inline void refine_ed_aln(overlap_region_alloc* overlap_list, All_reads *rref, const ul_idx_t *uref, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read, kvec_t_u64_warp* v_idx, int64_t block_s, double e_rate, double e_rate_final)
{
    int64_t j, k, i, on, y_id, y_readLen, x_start, x_end, x_len, total_y_start, total_y_end;
    int32_t y_strand, real_y_start;
    int64_t nw, a_nw, w_id, w_s, w_e, is_srt, mm_we, mm_ws, mm_aln, ovl;
    double error_rate; uint64_t *w_idx; overlap_region *z; window_list *p = NULL;

    overlap_list->mapped_overlaps_length = 0; on = overlap_list->length;
    for (j = 0; j < on; ++j) {
        // z = &(overlap_list->list[j]); ovl = z->x_pos_e+1-z->x_pos_s;
        // if(!realign_ed(z, uref, NULL, rref, g_read->seq, 
        //         dumy->overlap_region, dumy->overlap_region_fix, dumy, v_idx, block_s, e_rate, NULL, 1, &is_srt)) {
        //     continue;
        // }
        z = &(overlap_list->list[j]); z->is_match = 0; is_srt = 1;
        if(z->w_list.n == 0) continue;///no alignment
        nw = get_num_wins(z->x_pos_s, z->x_pos_e+1, block_s); a_nw = z->w_list.n;
        kv_resize(uint64_t, v_idx->a, (uint64_t)nw); 
        memset(v_idx->a.a, -1, sizeof((*v_idx->a.a))*nw); w_idx = v_idx->a.a;
        for (i = 0; i < a_nw; i++) { ///w_idx[] == (uint64_t) if unmatched
            assert(z->w_list.a[i].y_end != -1);
            w_id = get_win_id_by_s(z, z->w_list.a[i].x_start, block_s, NULL);
            w_idx[w_id] = i;
        }

        y_id = z->y_id; y_strand = z->y_pos_strand; 
        ovl = z->x_pos_e+1-z->x_pos_s; mm_we = z->x_pos_s; mm_aln = 0; 
        y_readLen = (rref?(Get_READ_LENGTH((*rref), y_id)):(uref->ug->u.a[y_id].len));
        for (i = a_nw-1; i >= 0; i--) { //utilize the the end pos of pre-window in forward
            w_id = get_win_id_by_s(z, z->w_list.a[i].x_start, block_s, &w_e);
            assert(z->w_list.a[i].x_end == w_e);
            if(w_e > mm_we) mm_we = w_e;
            ///in most cases, extra_begin = 0
            total_y_start = z->w_list.a[i].y_end + 1 - z->w_list.a[i].extra_begin;
            for (k = w_id + 1; k < nw && total_y_start < y_readLen; k++) {
                if(w_idx[k] != (uint64_t)-1) break;
                w_s = w_e + 1;
                w_id = get_win_id_by_s(z, w_s, block_s, &w_e);
                assert(w_id == k);
                x_start = w_s; x_end = w_e;
                if(aln_wlst(z, rref, uref, g_read, dumy, y_strand, y_id, x_start, x_end, 
                                                            total_y_start, block_s, e_rate, 0)) {
                    p = &(z->w_list.a[z->w_list.n-1]);
                    w_idx[k] = z->w_list.n - 1;
                    if(x_end > mm_we) mm_we = x_end;
                    if(is_srt && z->w_list.n > 1 && p->x_start < z->w_list.a[z->w_list.n-2].x_start) is_srt = 0;
                } else {
                    break;
                }
                total_y_start = p->y_end + 1 - p->extra_begin;
            }
        }
        mm_ws = z->x_pos_s; mm_aln = mm_we+1-mm_ws; 
        if(!simi_pass(ovl, mm_aln, uref?1:0, OVERLAP_THRESHOLD_NOSI_FILTER, NULL)) continue;
        if(nw > 0 && w_idx[0] != (uint64_t)-1) mm_ws = z->w_list.a[w_idx[0]].x_end+1;
        
        for (i = 1; i < nw; i++) { //utilize the the start pos of next window in backward
            ///find the first matched window, which should not be the first window
            ///the pre-window of this matched window must be unmatched
            if(w_idx[i] != (uint64_t)-1 && w_idx[i-1] == (uint64_t)-1) {
                w_s = z->w_list.a[w_idx[i]].x_start; mm_aln -= (w_s-mm_ws);
                ///check if the start pos of this matched window has been calculated
                if(z->w_list.a[w_idx[i]].clen == 0) {
                    p = &(z->w_list.a[w_idx[i]]);
                    gen_backtrace(p, z, rref, uref, g_read, dumy, y_strand, y_id);
                    assert(p->error != -1);
                    p->y_end += p->extra_begin;
                } 
                real_y_start = p->y_start;

                ///the end pos for pre window is real_y_start - 1
                total_y_end = real_y_start - 1;
                ///find the unmatched window on the left of current matched window
                ///k starts from i - 1
                for (k = i - 1; k >= 0 && w_idx[k] == (uint64_t)-1 && total_y_end > 0; k--) {  
                    w_e = w_s - 1;
                    w_id = get_win_id_by_e(z, w_e, block_s, &w_s);
                    assert(w_id == k);
                    x_start = w_s; x_end = w_e; x_len = x_end + 1 - x_start;
                    if(aln_wlst(z, rref, uref, g_read, dumy, y_strand, y_id, x_start, x_end, total_y_end+1-x_len, block_s, e_rate, 1)) {
                        p = &(z->w_list.a[z->w_list.n-1]);
                        p->y_start -= p->extra_begin; ///y_start has no shift, but y_end has shift  
                        w_idx[k] = z->w_list.n - 1;
                        mm_aln += x_len;
                        if(is_srt && z->w_list.n > 1 && p->x_start < z->w_list.a[z->w_list.n-2].x_start) is_srt = 0;
                    } else {
                        break;
                    }
                    total_y_end = p->y_start - 1;
                }
                if(!simi_pass(ovl, mm_aln, uref?1:0, OVERLAP_THRESHOLD_NOSI_FILTER, NULL)) break;
            }
            if(w_idx[i] != (uint64_t)-1) mm_ws = z->w_list.a[w_idx[i]].x_end+1;
        }

        if(i < nw) continue;
        if(uref && simi_pass(ovl, z->align_length, uref?1:0, OVERLAP_THRESHOLD_NOSI_FILTER, NULL)) {
            z->is_match = 3; overlap_list->mapped_overlaps_length += z->align_length;
            ///sort for set_herror_win
            if(!is_srt) radix_sort_window_list_xs_srt(z->w_list.a, z->w_list.a + z->w_list.n);
        }
    }

    if(uref && overlap_list->mapped_overlaps_length > 0) {
        set_herror_win(overlap_list, dumy, v_idx, e_rate, g_read->length, block_s);
    }

    overlap_list->mapped_overlaps_length = 0;
    for (j = 0; j < (long long)overlap_list->length; j++) {
        z = &(overlap_list->list[j]);
        y_id = z->y_id; y_strand = z->y_pos_strand; 
        y_readLen = (rref?(Get_READ_LENGTH((*rref), y_id)):(uref->ug->u.a[y_id].len));
        ovl = z->x_pos_e + 1 - z->x_pos_s; //z->is_match = 0;
        // if(y_id == 4) {
        //     fprintf(stderr, "[M::%s::idx->%ld::] z::x_pos_s->%u, z::x_pos_e->%u, ovl->%ld, aln->%u\n", 
        //     __func__, j, z->x_pos_s, z->x_pos_e, ovl, z->align_length);
        // }
        ///debug_scan_cigar(&(overlap_list->list[j]));
        ///only calculate cigar for high quality overlaps
        // int64_t tt = 0;
        if(simi_pass(ovl, z->align_length, 0, OVERLAP_THRESHOLD_NOSI_FILTER, &e_rate)) {
            a_nw = z->w_list.n;
            for (i = 0, is_srt = 1; i < a_nw; i++) {
                p = &(z->w_list.a[i]);
                ///check if the cigar of this window has been got 
                if(p->clen == 0) {
                    gen_backtrace(p, z, rref, uref, g_read, dumy, y_strand, y_id);
                    assert(p->error != -1);
                    // if(y_id == 4) {
                    //     fprintf(stderr, "+[M::idx->%ld::] y_start->%d, y_end->%d, error->%d\n", 
                    //     j, p->y_start, p->y_end, p->error);
                    // }
                }
                else {
                    p->y_end -= p->extra_begin;
                    // if(y_id == 4) {
                    //     fprintf(stderr, "-[M::idx->%ld::] y_start->%d, y_end->%d, error->%d\n", 
                    //     j, p->y_start, p->y_end, p->error);
                    // }
                }
                // tt += p->error;
                if(is_srt && i > 0 && p->x_start < z->w_list.a[i-1].x_start) is_srt = 0;
            }
            if(!is_srt) radix_sort_window_list_xs_srt(z->w_list.a, z->w_list.a + z->w_list.n);
            error_rate = non_trim_error_rate(z, rref, uref, v_idx, dumy, g_read, e_rate, block_s);
            z->is_match = 0;///must be here;
            // if(y_id == 4) {
            //     fprintf(stderr, "[M::%s::idx->%ld::] block_s->%ld, z::x_pos_s->%u, z::x_pos_e->%u, ovl->%ld, aln->%u, error_rate->%f, e_rate_final->%f\n", 
            //     __func__, j, block_s, z->x_pos_s, z->x_pos_e, ovl, z->align_length, error_rate, e_rate_final);
            //     exit(1);
            // }
            if (error_rate <= e_rate_final) {
                overlap_list->mapped_overlaps_length += ovl;
                z->is_match = 1; append_unmatched_wins(z, block_s);
                if(rref) {
                    calculate_boundary_cigars(z, rref, dumy, g_read, e_rate);
                } else {
                    calculate_ul_boundary_cigars(z, uref, dumy, g_read, e_rate, block_s);
                }
                // assert(get_num_wins(z->x_pos_s, z->x_pos_e+1, block_s)==(int64_t)z->w_list.n);
                // assert((int64_t)z->x_pos_s==z->w_list.a[0].x_start && 
                //                             (int64_t)z->x_pos_e==z->w_list.a[z->w_list.n-1].x_end);
            } else if (error_rate <= e_rate_final * 1.5) {
                z->is_match = 3;
            }
        } else {///it impossible to be matched
            z->is_match = 0;
            // fprintf(stderr, "[M::%s::idx->%ld::is_match->%u] z::x_pos_s->%u, z::x_pos_e->%u, error_rate->-1, e_threshold->%f\n", 
            // __func__, j, z->is_match, z->x_pos_s, z->x_pos_e, e_rate);
        }
    }
}


uint32_t align_ul_ed_post(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, char* qstr, char *tstr, char *tstr_1, 
Correct_dumy* dumy, double e_rate, int64_t w_l, double ovlp_cut, void *km);
double gen_extend_err(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, 
char *tstr, char *tstr_1, Correct_dumy* dumy, uint64_t *v_idx, int64_t block_s, double ovlp_cut, double e_rate, double e_max, int64_t *r_e);
inline void refine_ed_aln_test(overlap_region_alloc* overlap_list, All_reads *rref, const ul_idx_t *uref, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read, kvec_t_u64_warp* v_idx, int64_t block_s, double e_rate, double e_rate_final)
{
    int64_t j, on, ovl; uint64_t k; double rr; overlap_region *z; 

    overlap_list->mapped_overlaps_length = 0; on = overlap_list->length;
    for (j = 0; j < on; ++j) {
        z = &(overlap_list->list[j]); ovl = z->x_pos_e+1-z->x_pos_s;
        if(!align_ul_ed_post(z, uref, NULL, g_read->seq, dumy->overlap_region, dumy->overlap_region_fix, 
                        dumy, e_rate, block_s, OVERLAP_THRESHOLD_NOSI_FILTER, NULL)) {
            continue;
        }
        if(uref && simi_pass(ovl, z->align_length, uref?1:0, OVERLAP_THRESHOLD_NOSI_FILTER, NULL)) {
            z->is_match = 3; overlap_list->mapped_overlaps_length += z->align_length;
        }
    }

    if(uref && overlap_list->mapped_overlaps_length > 0) {
        set_herror_win(overlap_list, dumy, v_idx, e_rate, g_read->length, block_s);
    }

    double e_max = e_rate_final * 1.5;
    overlap_list->mapped_overlaps_length = 0; on = overlap_list->length;
    for (j = 0; j < on; j++) {
        z = &(overlap_list->list[j]); ovl = z->x_pos_e + 1 - z->x_pos_s;
        rr = gen_extend_err(z, uref, NULL, rref, g_read->seq, dumy->overlap_region, dumy->overlap_region_fix,
                    dumy, v_idx?v_idx->a.a:NULL, block_s, -1, e_rate, (e_max+0.000001), NULL);
        z->is_match = 0;///must be here;
        if (rr <= e_rate_final) {
            for (k = 0; k < z->w_list.n; k++) {
                if(z->w_list.a[k].clen) continue;
                gen_backtrace_adv(&(z->w_list.a[k]), z, rref, NULL, uref, g_read->seq, dumy->overlap_region, dumy->overlap_region_fix, 
                dumy, z->y_pos_strand, z->y_id);
            }
            
            overlap_list->mapped_overlaps_length += ovl;
            z->is_match = 1; append_unmatched_wins(z, block_s);
            if(rref) {
                calculate_boundary_cigars(z, rref, dumy, g_read, e_rate);
            } else {
                calculate_ul_boundary_cigars(z, uref, dumy, g_read, e_rate, block_s);
            }
        } else if (rr <= e_max) {
            z->is_match = 3;
        }
    }
}



inline void add_base_to_correct_read_directly(Correct_dumy* dumy, char base)
{
  
    if (dumy->corrected_read_length + 2 > dumy->corrected_read_size)
    {
        dumy->corrected_read_size = dumy->corrected_read_size * 2;
        dumy->corrected_read = (char*)realloc(dumy->corrected_read, dumy->corrected_read_size);
    }
    
    dumy->corrected_read[dumy->corrected_read_length] = base;
    dumy->corrected_read_length++;
    dumy->corrected_read[dumy->corrected_read_length] = '\0';

}

inline void add_base_to_correct_read(Correct_dumy* dumy, char base, int is_error)
{
    ///don't need to deal with deletion
    if (base != 'D')
    {
        if (dumy->corrected_read_length + 2 > dumy->corrected_read_size)
        {
            dumy->corrected_read_size = dumy->corrected_read_size * 2;
            dumy->corrected_read = (char*)realloc(dumy->corrected_read, dumy->corrected_read_size);
        }
        
        dumy->corrected_read[dumy->corrected_read_length] = base;
        dumy->corrected_read_length++;
        dumy->corrected_read[dumy->corrected_read_length] = '\0';
    }

    if (is_error)
    {
        dumy->corrected_base++;
    }
    
}


inline void add_segment_to_correct_read(Correct_dumy* dumy, char* segment, long long segment_length)
{

    if (dumy->corrected_read_length + segment_length + 2 > dumy->corrected_read_size)
    {
        dumy->corrected_read_size = dumy->corrected_read_length + segment_length + 2;
        dumy->corrected_read = (char*)realloc(dumy->corrected_read, dumy->corrected_read_size);
    }

    memcpy(dumy->corrected_read + dumy->corrected_read_length, segment, segment_length);
    dumy->corrected_read_length += segment_length;
    dumy->corrected_read[dumy->corrected_read_length] = '\0';
}




///return the ID of next node at backbone
long long inline add_path_to_correct_read(Graph* backbone, Correct_dumy* dumy, long long currentNodeID, 
long long type, long long edgeID, Cigar_record* current_cigar, char* self_string)
{
    //long long i;
    long long nodeID;

    ///Note: currentNodeID must be a backbone node
    ///currentNodeID = 0 means a fake node
    ///currentNodeID = i means self_string[i - 1]
    ///include match/mismatch
    if (type == MISMATCH)
    {
        ///match
        if(backbone->g_nodes.list[currentNodeID].mismatch_edges.list[edgeID].length == 0)
        {
            nodeID = backbone->g_nodes.list[currentNodeID].mismatch_edges.list[edgeID].out_node;
            add_base_to_correct_read_directly(dumy, backbone->g_nodes.list[nodeID].base);
            ///nodeID = i means self_string[i - 1]
            ///add_cigar_record(self_string+nodeID-1, 1, current_cigar, 0);
            add_cigar_record(&(backbone->g_nodes.list[nodeID].base), 1, current_cigar, 0);
            return nodeID;
        }
        else ///mismatch
        {
            nodeID = backbone->g_nodes.list[currentNodeID].mismatch_edges.list[edgeID].out_node;
            add_base_to_correct_read_directly(dumy, backbone->g_nodes.list[nodeID].base);
            dumy->corrected_base++;
            char merge_base = 0;
            merge_base = seq_nt6_table[(uint8_t)backbone->g_nodes.list[nodeID].base];
            merge_base = merge_base << 3;
            nodeID = backbone->g_nodes.list[nodeID].mismatch_edges.list[0].out_node;
            merge_base = merge_base | seq_nt6_table[(uint8_t)backbone->g_nodes.list[nodeID].base];
            add_cigar_record(&merge_base, 1, current_cigar, 1);
            return nodeID;
        }
    }
    else if (type == DELETION)
    {
        nodeID = backbone->g_nodes.list[currentNodeID].deletion_edges.list[edgeID].out_node;
        dumy->corrected_base += nodeID - currentNodeID;
        add_cigar_record(self_string + currentNodeID, nodeID - currentNodeID, current_cigar, DELETION);
        return nodeID;
    }
    else if (type == INSERTION)
    {
        ///pay attention to this line 
        backbone->g_nodes.list[currentNodeID].num_insertions = 0;

        nodeID = backbone->g_nodes.list[currentNodeID].insertion_edges.list[edgeID].out_node;
        long long step = backbone->g_nodes.list[currentNodeID].insertion_edges.list[edgeID].length;
        long long i;
        for (i = 0; i < step; i++)
        {
            add_base_to_correct_read_directly(dumy, backbone->g_nodes.list[nodeID].base);
            add_cigar_record(&backbone->g_nodes.list[nodeID].base, 1, current_cigar, INSERTION);
            nodeID = backbone->g_nodes.list[nodeID].insertion_edges.list[0].out_node;
        }
        dumy->corrected_base += step;
        return nodeID;
    }
    else
    {
        fprintf(stderr, "error type\n");
        return -1;
    }
}


///return the ID of next node at backbone
long long inline add_path_to_correct_read_new(Graph* backbone, Graph* DAGCon, Correct_dumy* dumy, long long currentNodeID, 
long long type, long long edgeID, Cigar_record* current_cigar, char* self_string)
{
    //long long i;
    long long nodeID;

    ///Note: currentNodeID must be a backbone node
    ///currentNodeID = 0 means a fake node
    ///currentNodeID = i means self_string[i - 1]
    ///include match/mismatch
    if (type == MISMATCH)
    {
        ///match
        if(backbone->g_nodes.list[currentNodeID].mismatch_edges.list[edgeID].length == 0)
        {
            nodeID = backbone->g_nodes.list[currentNodeID].mismatch_edges.list[edgeID].out_node;
            add_base_to_correct_read_directly(dumy, backbone->g_nodes.list[nodeID].base);
            ///nodeID = i means self_string[i - 1]
            ///add_cigar_record(self_string+nodeID-1, 1, current_cigar, 0);
            add_cigar_record(&(backbone->g_nodes.list[nodeID].base), 1, current_cigar, 0);
            return nodeID;
        }
        else ///mismatch
        {
            nodeID = backbone->g_nodes.list[currentNodeID].mismatch_edges.list[edgeID].out_node;
            add_base_to_correct_read_directly(dumy, backbone->g_nodes.list[nodeID].base);
            dumy->corrected_base++;
            char merge_base = 0;
            merge_base = seq_nt6_table[(uint8_t)backbone->g_nodes.list[nodeID].base];
            merge_base = merge_base << 3;
            nodeID = backbone->g_nodes.list[nodeID].mismatch_edges.list[0].out_node;
            merge_base = merge_base | seq_nt6_table[(uint8_t)backbone->g_nodes.list[nodeID].base];
            add_cigar_record(&merge_base, 1, current_cigar, 1);
            return nodeID;
        }
    }
    else if (type == DELETION)
    {
        nodeID = backbone->g_nodes.list[currentNodeID].deletion_edges.list[edgeID].out_node;
        dumy->corrected_base += nodeID - currentNodeID;
        add_cigar_record(self_string + currentNodeID, nodeID - currentNodeID, current_cigar, DELETION);
        return nodeID;
    }
    else if (type == INSERTION)
    {
        ///pay attention to this line
        backbone->g_nodes.list[currentNodeID].num_insertions = 0;
        long long str;
        char str_c;
        while (pop_from_Queue(&(DAGCon->node_q), &str))
        {
            str_c = (char)str;
            add_base_to_correct_read_directly(dumy, str_c);
            add_cigar_record(&str_c, 1, current_cigar, INSERTION);
            dumy->corrected_base++;
        }

        return currentNodeID;
    }
    else
    {
        fprintf(stderr, "error type\n");
    }
    
    return -1;
}


void Merge_Out_Nodes(Graph* DAGCon, Node* currentNode)
{
    ///if this node does not have any output, directly return
    if(Real_Length(Output_Edges((*currentNode))) == 0)
    {
        return;
    }

    RSet buf, out_buf;
    char Bases[4] = {'A', 'C', 'G', 'T'};
    char base;
    long long base_i, weight;
    int flag = 0;
    Node* get_node_1 = NULL;
    Node* out_node_of_get_node_1 = NULL;
    Node* consensus_node_1 = NULL;
    Edge* e_forward_1 = NULL;
    Edge* e_backward_1 = NULL;

    ///merge all base for each base
    for (base_i = 0; base_i < 4; base_i++)
    {
        base = Bases[base_i];

        clear_RSet(&buf);
        flag = 0;
        weight = 0;
        ///should use getOutputEdges, instead of getOutputNodes
        ///check all out-nodes of currentNode
        while(getOutputNodes(&buf, DAGCon, currentNode, &get_node_1))
        {
            ///check the corresponding node, this node must only have one in-node
            ///note this is the Real_Length, instead of the Input_Edges.length
            if((*get_node_1).base == base && Real_Length(Input_Edges(*get_node_1)) == 1)
            {
                if(flag == 0)
                {
                    flag = 1;
                    ///add a new node to merge all out-node
                    consensus_node_1 = get_node_1;
                    ///link consensus_node to currentNode
                    ///set the new edge to be visited
                    if(get_bi_Edge(DAGCon, currentNode, consensus_node_1, &e_forward_1, &e_backward_1))
                    {
                        Visit(*e_forward_1) = 1;
                        Visit(*e_backward_1) = 1;
                    }
                    else
                    {
                        fprintf(stderr, "error\n");
                    }

                    weight = (*e_forward_1).weight;
                }
                else
                {
                    flag++;
                    ///add the weight of get_node->currentNode
                    if(get_bi_Edge(DAGCon, currentNode, get_node_1, &e_forward_1, &e_backward_1))
                    {
                        weight = weight + (*e_forward_1).weight;
                    }
                    else
                    {
                        fprintf(stderr, "error\n");
                    }

                    ///process the out-nodes of get_node
                    clear_RSet(&out_buf);
                    while(getOutputNodes(&out_buf, DAGCon, get_node_1, &out_node_of_get_node_1))
                    {
                        ///link consensus_node to the out-nodes of get_node
                        if(get_bi_Edge(DAGCon, consensus_node_1, out_node_of_get_node_1, &e_forward_1, &e_backward_1))
                        {
                            Visit(*e_forward_1) = 1;
                            Visit(*e_backward_1) = 1;
                            (*e_forward_1).weight += get_Edge_Weight(DAGCon, get_node_1, out_node_of_get_node_1);
                            (*e_backward_1).weight = (*e_forward_1).weight;
                        }
                        else
                        {
                            add_bi_direction_edge(DAGCon, consensus_node_1, out_node_of_get_node_1, 
                            get_Edge_Weight(DAGCon, get_node_1, out_node_of_get_node_1), 1);
                        }
                    }

                    delete_Node_DAGCon(DAGCon, get_node_1);
                }
            }
        }

        if(flag > 1)
        {
            get_bi_Edge(DAGCon, currentNode, consensus_node_1, &e_forward_1, &e_backward_1);
            (*e_forward_1).weight = weight;
            (*e_backward_1).weight = (*e_forward_1).weight;
        }

        if(flag > 0)
        {
            Merge_Out_Nodes(DAGCon, consensus_node_1);   
        }
    }
}


void Merge_In_Nodes(Graph* DAGCon, Node* currentNode)
{
    ///if this node does not have any input, directly return
    if(Real_Length(Input_Edges((*currentNode))) == 0)
    {
        return;
    }

    RSet buf, in_buf;
    char Bases[4] = {'A', 'C', 'G', 'T'};
    char base;
    long long base_i, weight;
    int flag = 0;
    Node* get_node = NULL;
    Node* in_node_of_get_node = NULL;
    Node* consensus_node = NULL;
    Edge* e_forward = NULL;
    Edge* e_backward = NULL;

    ///merge all base for each base
    for (base_i = 0; base_i < 4; base_i++)
    {
        base = Bases[base_i];

        clear_RSet(&buf);
        flag = 0;
        weight = 0;
        ///should use getInputEdges, instead of getInputNodes
        ///check all in-nodes of currentNode
        while(getInputNodes(&buf, DAGCon, currentNode, &get_node))
        {
            ///check the corresponding node, this node must only have one out-node
            ///note this is the Real_Length, instead of the Output_Edges.length
            if((*get_node).base == base && Real_Length(Output_Edges(*get_node)) == 1)
            {
                if(flag == 0)
                {
                    flag = 1;
                    ///add a new node to merge all in-node
                    consensus_node = get_node;
                    ///link consensus_node to currentNode
                    ///set the new edge to be visited
                    if(get_bi_Edge(DAGCon, consensus_node, currentNode, &e_forward, &e_backward))
                    {
                        Visit(*e_forward) = 1;
                        Visit(*e_backward) = 1;
                    }
                    else
                    {
                        fprintf(stderr, "error\n");
                    }

                    weight = (*e_forward).weight;
                }
                else
                {
                    flag++;
                    ///add the weight of get_node->currentNode
                    if(get_bi_Edge(DAGCon, get_node, currentNode, &e_forward, &e_backward))
                    {
                        weight = weight + (*e_forward).weight;
                    }
                    else
                    {
                        fprintf(stderr, "error\n");
                    }

                    ///process the in-nodes of get_node
                    clear_RSet(&in_buf);
                    while(getInputNodes(&in_buf, DAGCon, get_node, &in_node_of_get_node))
                    {
                        ///link in-nodes of get_node to consensus_node
                        if(get_bi_Edge(DAGCon, in_node_of_get_node, consensus_node, &e_forward, &e_backward))
                        {
                            Visit(*e_forward) = 1;
                            Visit(*e_backward) = 1;
                            (*e_forward).weight += get_Edge_Weight(DAGCon, in_node_of_get_node, get_node);
                            (*e_backward).weight = (*e_forward).weight;
                        }
                        else
                        {
                            add_bi_direction_edge(DAGCon, in_node_of_get_node, consensus_node, 
                            get_Edge_Weight(DAGCon, in_node_of_get_node, get_node), 1);
                        }
                    }

                    delete_Node_DAGCon(DAGCon, get_node);
                }
            }
        }

        if(flag > 1)
        {
            get_bi_Edge(DAGCon, consensus_node, currentNode, &e_forward, &e_backward);
            (*e_forward).weight = weight;
            (*e_backward).weight = (*e_forward).weight;
        }

        if(flag > 0)
        {
            Merge_In_Nodes(DAGCon, consensus_node);   
        }
    }
}

void print_graph(Graph* DAGCon)
{
    uint64_t i;
    for (i = 0; i < DAGCon->g_nodes.length; i++)
    {
        Node* currentStartNode = &(G_Node(*DAGCon, i));
        RSet iter_out;

        if(If_Node_Exist(*currentStartNode))
        {
            fprintf(stderr, "ID: %lu (%c) (w: %lu)\n", (unsigned long)(*currentStartNode).ID, (*currentStartNode).base, (unsigned long)(*currentStartNode).weight);
            clear_RSet(&iter_out);
            Edge* e;
            fprintf(stderr, "****Out-node: ");
            while(getOutputEdges(&iter_out, DAGCon, currentStartNode, &e))
            {
                //fprintf(stderr, "%d[%c], ", G_Node(*DAGCon, e->out_node).ID, G_Node(*DAGCon, e->out_node).base);
                fprintf(stderr, "%lu(w: %lu), ", (unsigned long)G_Node(*DAGCon, e->out_node).ID, (unsigned long)e->weight);
            }
            fprintf(stderr, "\n");


            // clear_RSet(&iter_out);
            // fprintf(stderr, "In-node: ");
            // while(getInputEdges(&iter_out, DAGCon, currentStartNode, &e))
            // {
            //     fprintf(stderr, "%d[%c], ", G_Node(*DAGCon, e->in_node).ID, G_Node(*DAGCon, e->in_node).base);
            // }
        }
    }

    fprintf(stderr, "*******\n");

}

void debug_DAGCon(Graph* DAGCon)
{
    uint64_t i = 0;
    for (i = 0; i < DAGCon->g_nodes.length; i++)
    {
        Node* currentStartNode = &(G_Node(*DAGCon, i));
        RSet iter_out;

        if(If_Node_Exist(*currentStartNode))
        {
            clear_RSet(&iter_out);
            Edge* e_self;
            Edge* e_reverse;
            while(getOutputEdges(&iter_out, DAGCon, currentStartNode, &e_self))
            {

                get_bi_direction_edges(DAGCon, e_self, &e_self, &e_reverse);


                if(Visit(*e_self) == 0)
                {
                    fprintf(stderr, "Visit(*e_self): %lu, error visit flag: in_node: %lu, out_node: %lu\n", 
                    (unsigned long)Visit(*e_self), (unsigned long)(*e_self).in_node, (unsigned long)(*e_self).out_node);
                }


                if(Visit(*e_reverse) == 0)
                {
                    fprintf(stderr, "Visit(*e_reverse): %lu, error visit flag: in_node: %lu, out_node: %lu\n", 
                    (unsigned long)Visit(*e_reverse), (unsigned long)(*e_reverse).in_node, (unsigned long)(*e_reverse).out_node);
                }


                if(e_self->in_node != e_reverse->in_node)
                {
                    fprintf(stderr, "different in-node\n");
                }
                if(e_self->out_node != e_reverse->out_node)
                {
                    fprintf(stderr, "different out-node\n");
                }
                if(e_self->weight != e_reverse->weight)
                {
                    fprintf(stderr, "different weight\n");
                }
            }

            clear_RSet(&iter_out);
            while(getInputEdges(&iter_out, DAGCon, currentStartNode, &e_self))
            {

                get_bi_direction_edges(DAGCon, e_self, &e_self, &e_reverse);


                if(Visit(*e_self) == 0)
                {
                    fprintf(stderr, "Visit(*e_self): %lu, error visit flag: in_node: %lu, out_node: %lu\n", 
                    (unsigned long)Visit(*e_self), (unsigned long)(*e_self).in_node, (unsigned long)(*e_self).out_node);
                }


                if(Visit(*e_reverse) == 0)
                {
                    fprintf(stderr, "Visit(*e_reverse): %lu, error visit flag: in_node: %lu, out_node: %lu\n", 
                    (unsigned long)Visit(*e_reverse), (unsigned long)(*e_reverse).in_node, (unsigned long)(*e_reverse).out_node);
                }


                if(e_self->in_node != e_reverse->in_node)
                {
                    fprintf(stderr, "different in-node\n");
                }
                if(e_self->out_node != e_reverse->out_node)
                {
                    fprintf(stderr, "different out-node\n");
                }
                if(e_self->weight != e_reverse->weight)
                {
                    fprintf(stderr, "different weight\n");
                }
            }
        }
    }
}

void Merge_DAGCon(Graph* DAGCon)
{
    ///using the length of edge representing if it has been visited
    ///in default, the length of edge is 0
    RSet iter_node, iter_edge;
    long long flag;

    Node* currentNode;
    Node* outNode;
    Edge* edge;
    Edge* e_forward;
    Edge* e_backward;

    // int num_way = Real_Length(Output_Edges(G_Node(*DAGCon, DAGCon->s_start_nodeID)));
    // if(num_way > 2)
    // {
    //     print_graph(DAGCon);
    // }
    

    
    ///at begining, only the start node has no in-node
    currentNode = &(G_Node(*DAGCon, DAGCon->s_start_nodeID));
    Push_Node(DAGCon, &currentNode);
    

    
    while (Pop_Node(DAGCon, &currentNode))
    {
        ///merge in-node
        Merge_In_Nodes(DAGCon, currentNode);
        ///merge out-node
        Merge_Out_Nodes(DAGCon, currentNode);

        clear_RSet(&iter_edge);
        ///for all out-edges of currentNode, set as visited
        while (getOutputEdges(&iter_edge, DAGCon, currentNode, &edge))
        {
            get_bi_direction_edges(DAGCon, edge, &e_forward, &e_backward);
            Visit(*e_forward) = 1;
            Visit(*e_backward) = 1;
        }
        


        ///check all out-node of currentNode
        clear_RSet(&iter_node);
        while(getOutputNodes(&iter_node, DAGCon, currentNode, &outNode))
        {

            ///for each outNode, check if all in-edges have been visited
            flag = 0;
            clear_RSet(&iter_edge);
            while (getInputEdges(&iter_edge, DAGCon, outNode, &edge))
            {
                if(Visit(*edge) == 0)
                {
                    flag = 1;
                    break;
                }
            }
            //if all in-edges of Out_node have already been visited, push it to queue
            if(flag == 0)
            {
                Push_Node(DAGCon, &outNode);
            }
        }
    }



    // if(num_way > 2)
    // {
    //     print_graph(DAGCon);
    //     fprintf(stderr, "****************************note*****************\n\n");
    // }


    ///debug_DAGCon(DAGCon);
}


inline void generate_seq_from_path(Graph* DAGCon, Node* node, int direction)
{
    clear_Queue(&(DAGCon->node_q));
    RSet iter;
    Edge* e = NULL;
    uint64_t max;
    Node* max_node = NULL;
    

    if(direction == 0)
    {
        while (node->ID != DAGCon->s_end_nodeID)
        {
            push_to_Queue(&(DAGCon->node_q), node->base);
            clear_RSet(&iter);
            max = 0;
            while(getOutputEdges(&iter, DAGCon, node, &e))
            {
                if(e->weight > max)
                {
                    max = e->weight;
                    max_node = &(G_Node(*DAGCon, e->out_node));
                }
            }
            node = max_node;
        }
    }
    else
    {
        while (node->ID != DAGCon->s_start_nodeID)
        {
            push_to_Queue(&(DAGCon->node_q), node->base);
            clear_RSet(&iter);
            max = 0;
            while(getInputEdges(&iter, DAGCon, node, &e))
            {
                if(e->weight > max)
                {
                    max = e->weight;
                    max_node = &(G_Node(*DAGCon, e->in_node));
                }
            }
            node = max_node;
        }

        long long i, k;
        long long length = (DAGCon->node_q.end - DAGCon->node_q.beg);
        long long length_ex = length/2;
        long long* array = DAGCon->node_q.buffer + DAGCon->node_q.beg;
        for (i = 0; i < length_ex; i++)
        {
            k = array[i];
            array[i] = array[length - i - 1];
            array[length - i - 1] = k;
        }
    }
}


long long generate_best_seq_from_edges(Graph* DAGCon)
{
    long long max_start = 0, max_end = 0, max_start_edge = 0, max_end_edge = 0;
    RSet iter;
    Edge* e = NULL;
    Node* newNode = NULL;
    long long max_count = 0;
    

    ///check the out-edges of start node
    ///must to be 0
    max_start = 0;
    newNode = &(G_Node(*DAGCon, DAGCon->s_start_nodeID));
    clear_RSet(&iter);
    while(getOutputEdges(&iter, DAGCon, newNode, &e))
    {
        if(e->weight > (uint64_t)max_start)
        {
            max_start = e->weight;
            max_start_edge = iter.index - 1;
        }
    }

    ///check the in-edges of end node
    ///must to be 0
    max_end = 0;
    newNode = &(G_Node(*DAGCon, DAGCon->s_end_nodeID));
    clear_RSet(&iter);
    while(getInputEdges(&iter, DAGCon, newNode, &e))
    {
        if(e->weight > (uint64_t)max_end)
        {
            max_end = e->weight;
            max_end_edge = iter.index - 1;
        }
    }

    if(max_start >= max_end)
    {
        max_count = max_start;
        generate_seq_from_path(DAGCon, 
        &G_Node(*DAGCon, Output_Edges(G_Node(*DAGCon, DAGCon->s_start_nodeID)).list[max_start_edge].out_node), 0);
    }
    else
    {
        max_count = max_end;
        generate_seq_from_path(DAGCon, 
        &G_Node(*DAGCon, Input_Edges(G_Node(*DAGCon, DAGCon->s_end_nodeID)).list[max_end_edge].in_node), 1);
    }

    return max_count;
}


inline void generate_seq_from_node(Graph* DAGCon, Node* node, int direction)
{
    clear_Queue(&(DAGCon->node_q));
    RSet iter;
    uint64_t max;
    Node* max_node = NULL;
    Node* getNodes = NULL;
    

    if(direction == 0)
    {
        while (node->ID != DAGCon->s_end_nodeID)
        {
            push_to_Queue(&(DAGCon->node_q), node->base);
            clear_RSet(&iter);
            max = 0;
            while(getOutputNodes(&iter, DAGCon, node, &getNodes))
            {
                if(getNodes->weight > max)
                {
                    max = getNodes->weight;
                    max_node = getNodes;
                }
            }
            node = max_node;
        }
    }
    else
    {
        while (node->ID != DAGCon->s_start_nodeID)
        {
            push_to_Queue(&(DAGCon->node_q), node->base);
            clear_RSet(&iter);
            max = 0;
            while(getInputNodes(&iter, DAGCon, node, &getNodes))
            {
                if(getNodes->weight > max)
                {
                    max = getNodes->weight;
                    max_node = getNodes;
                }
            }
            node = max_node;
        }

        long long i, k;
        long long length = (DAGCon->node_q.end - DAGCon->node_q.beg);
        long long length_ex = length/2;
        long long* array = DAGCon->node_q.buffer + DAGCon->node_q.beg;
        for (i = 0; i < length_ex; i++)
        {
            k = array[i];
            array[i] = array[length - i - 1];
            array[length - i - 1] = k;
        }
    }
}


long long generate_best_seq_from_nodes(Graph* DAGCon)
{
    long long max_start = 0, max_end = 0;
    RSet iter;
    Edge* e = NULL;
    Node* newNode = NULL;
    Node* getNode = NULL;
    Node* max_start_node = NULL;
    Node* max_end_node = NULL;
    long long max_count = 0;
    uint64_t i;

   for (i = 0; i < DAGCon->g_nodes.length; i++)
   {
       newNode = &(G_Node(*DAGCon, i));
       if(If_Node_Exist(*newNode))
        {
            newNode->weight = 0;
            clear_RSet(&iter);
            while(getOutputEdges(&iter, DAGCon, newNode, &e))
            {
                newNode->weight += e->weight;
            }
        }
   }
    
    newNode = &(G_Node(*DAGCon, DAGCon->s_start_nodeID));
    newNode->weight = 0;
    clear_RSet(&iter);
    while(getOutputEdges(&iter, DAGCon, newNode, &e))
    {
        newNode->weight += e->weight;
    }


    newNode = &(G_Node(*DAGCon, DAGCon->s_end_nodeID));
    newNode->weight = 0;
    clear_RSet(&iter);
    while(getInputEdges(&iter, DAGCon, newNode, &e))
    {
        newNode->weight += e->weight;
    }


    

    ///check the out-edges of start node
    ///must to be 0
    max_start = 0;
    newNode = &(G_Node(*DAGCon, DAGCon->s_start_nodeID));
    clear_RSet(&iter);
    while(getOutputNodes(&iter, DAGCon, newNode, &getNode))
    {
        if(getNode->weight > (uint64_t)max_start)
        {
            max_start = getNode->weight;
            max_start_node = getNode;
        }
    }

    ///check the in-edges of end node
    ///must to be 0
    max_end = 0;
    newNode = &(G_Node(*DAGCon, DAGCon->s_end_nodeID));
    clear_RSet(&iter);
    while(getInputNodes(&iter, DAGCon, newNode, &getNode))
    {
        if(getNode->weight > (uint64_t)max_end)
        {
            max_end = getNode->weight;
            max_end_node = getNode;
        }
    }

    if(max_start >= max_end)
    {
        max_count = max_start;
        generate_seq_from_node(DAGCon, max_start_node, 0);
    }
    else
    {
        max_count = max_end;
        generate_seq_from_node(DAGCon, max_end_node, 1);
    }

    return max_count;
}


void build_DAGCon(Graph* DAGCon, Graph* backbone, long long currentNodeID, long long* max_count)
{
    long long i, j, path_weight, nodeID, step;
    char base;
    clear_Graph(DAGCon);
    Node* newNode;
    Node* lastNode;

    ///add the start node and the end node
    newNode = add_Node_DAGCon(DAGCon, 'S');
    DAGCon->s_start_nodeID = newNode->ID;

    newNode = add_Node_DAGCon(DAGCon, 'E');
    DAGCon->s_end_nodeID = newNode->ID;


    for (i = 0; i < (long long)G_Node(*backbone, currentNodeID).insertion_edges.length; i++)
    {
        path_weight = G_Node(*backbone, currentNodeID).insertion_edges.list[i].weight;
        
        lastNode = &(G_Node(*DAGCon, DAGCon->s_start_nodeID));

        step = G_Node(*backbone, currentNodeID).insertion_edges.list[i].length;
        if(step != 0)
        {
            nodeID = G_Node(*backbone, currentNodeID).insertion_edges.list[i].out_node;
        
            for (j = 0; j < step; j++)
            {
                base = G_Node(*backbone, nodeID).base;
                newNode = add_Node_DAGCon(DAGCon, base);
                add_bi_direction_edge(DAGCon, lastNode, newNode, path_weight, 0);

                nodeID = G_Node(*backbone, nodeID).insertion_edges.list[0].out_node;

                lastNode = newNode;
            }

            if(lastNode->ID != DAGCon->s_start_nodeID)
            {
                add_bi_direction_edge(DAGCon, lastNode, &(G_Node(*DAGCon, DAGCon->s_end_nodeID)), path_weight, 0);

            }
        }


        
    }


    Merge_DAGCon(DAGCon);

    ///(*max_count) = generate_best_seq_from_edges(DAGCon);
    (*max_count) = generate_best_seq_from_nodes(DAGCon);
    
    ///very important
    backbone->g_nodes.list[currentNodeID].num_insertions = 0;

}

void debug_whole_graph(Graph* g)
{
    long long i, j, k;
    for (i = 0; i < (long long)g->g_nodes.length; i++)
    {
        if(g->g_nodes.list[i].deletion_edges.length!= 0 &&
        g->g_nodes.list[i].deletion_edges.length!= 1)
        {
            fprintf(stderr, "g->g_nodes.list[i].deletion_edges.length: %lu\n",
            (unsigned long)g->g_nodes.list[i].deletion_edges.length);
        }
    }

    for (i = g->s_start_nodeID; i < (long long)g->s_end_nodeID; i++)
    {
        if(g->g_nodes.list[i].mismatch_edges.length > 4 
        || 
        g->g_nodes.list[i].mismatch_edges.length < 1)
        {
            fprintf(stderr, "g->s_end_nodeID: %lu, g->g_nodes.list[%lld].mismatch_edges.length: %lu\n",
            (unsigned long)g->s_end_nodeID, i, (unsigned long)g->g_nodes.list[i].mismatch_edges.length);
        }
    }

    char current[1000];
    char compare[1000];
    long long total_weight = 0;
    for (i = g->s_start_nodeID; i < (long long)g->s_end_nodeID; i++)
    {
        total_weight = 0;
        for (j = 0; j < (long long)G_Node(*g, i).insertion_edges.length; j++)
        {
            
            total_weight = total_weight + G_Node(*g, i).insertion_edges.list[j].weight;

            extract_path(g, i, j, current);

            for (k = j + 1; k < (long long)G_Node(*g, i).insertion_edges.length; k++)
            {
                extract_path(g, i, k, compare);
                if(strcmp(current, compare)==0)
                {
                    fprintf(stderr,"error\n");
                }
            }
        }

        if(total_weight != (long long)G_Node(*g, i).num_insertions)
        {
            fprintf(stderr,"error\n");
        }
    }

        

}

void get_seq_from_Graph(Graph* backbone, Graph* DAGCon, Correct_dumy* dumy, Cigar_record* current_cigar, char* self_string,
char* r_string, long long r_string_length, long long r_string_site)
{
    ///debug_whole_graph(backbone);
    long long currentNodeID;
    long long i;
    // There are several cases: 
    // 1. match 2. mismatch (A, C, G, T, N) 3. deletion 4. insertion (A, C, G, T)
    // in fact, 1. the weight of node itself 2. weight of alignToNode 3. weight of insertion node
    long long max_count;
    int max_type;
    long long  max_edge;
    long long total_count;
    long long current_weight;

    long long max_insertion_count;

    currentNodeID = backbone->s_start_nodeID;

    while (currentNodeID != (long long)backbone->s_end_nodeID)
    {
        total_count = 0;
        max_count = -1;
        max_type = -1;
        max_edge = -1;


        ///if it is a backbone node, there are three types od out-edges
        ///1. mismatch_edges 2. insertion_edges 3. deletion_edges
        if (currentNodeID >= (long long)backbone->s_start_nodeID && currentNodeID <= (long long)backbone->s_end_nodeID)
        {
            ///mismatch_edges
            for (i = 0; i < (long long)backbone->g_nodes.list[currentNodeID].mismatch_edges.length; i++)
            {
                if (backbone->g_nodes.list[currentNodeID].num_insertions != 0)
                {
                    current_weight = backbone->g_nodes.list[currentNodeID].mismatch_edges.list[i].weight -
                        backbone->g_nodes.list[currentNodeID].mismatch_edges.list[i].num_insertions;
                }
                else
                {
                    current_weight = backbone->g_nodes.list[currentNodeID].mismatch_edges.list[i].weight;
                }
                
                
                total_count = total_count + current_weight;

                ///for match, it needs to deal with both match and insertion
                ///if there is a insertion, we need to check this node two times
                ///1. num_insertions > 0, 2. num_insertions=0
                if (current_weight > max_count)
                {
                    max_count = current_weight;
                    max_edge = i;
                    max_type = MISMATCH;
                }
            }

            ///insertion_edges
            if (backbone->g_nodes.list[currentNodeID].num_insertions != 0)
            {
                ///this line must be prior than the next line
                ///since build_DAGCon will set backbone->g_nodes.list[currentNodeID].num_insertions to be 0
                total_count = total_count + backbone->g_nodes.list[currentNodeID].num_insertions;
                
                build_DAGCon(DAGCon, backbone, currentNodeID, &max_insertion_count);

                if(max_insertion_count > max_count)
                {
                    max_count = max_insertion_count;
                    max_type = INSERTION;
                }
            }

            ///deletion_edges
            for (i = 0; i < (long long)backbone->g_nodes.list[currentNodeID].deletion_edges.length; i++)
            {
                total_count = total_count + backbone->g_nodes.list[currentNodeID].deletion_edges.list[i].weight;

                if ((long long)backbone->g_nodes.list[currentNodeID].deletion_edges.list[i].weight > max_count)
                {
                    max_count = backbone->g_nodes.list[currentNodeID].deletion_edges.list[i].weight;
                    max_edge = i;
                    max_type = DELETION;
                }
            }

            ///do correction
            if(max_count >= total_count*(CORRECT_THRESHOLD))
            {
                currentNodeID = add_path_to_correct_read_new(backbone, DAGCon, dumy, currentNodeID, max_type, max_edge, current_cigar,
                self_string);
            }
            else  
            {
                 ///NOTE: currentNodeID = 0 is a tmp node without any sense
                 if(currentNodeID > 0 && if_is_homopolymer_strict(r_string_site + currentNodeID - 1, r_string, r_string_length)
                  && max_count >= total_count*CORRECT_THRESHOLD_HOMOPOLYMER)
                {
                    currentNodeID = add_path_to_correct_read_new(backbone, DAGCon, dumy, currentNodeID, max_type, max_edge, current_cigar,
                    self_string);
                }
                else///don't do correction, directly use the base of next backbone node
                {
                    currentNodeID++;
                    add_base_to_correct_read_directly(dumy, backbone->g_nodes.list[currentNodeID].base);

                    add_cigar_record(&(backbone->g_nodes.list[currentNodeID].base), 1, current_cigar, 0);
                }
            }
            
        }
        else  ///if there is a non-backbone node
        {
            fprintf(stderr, "error\n");
        }
        
    }

}


void window_consensus(char* r_string, long long r_total_length, long long window_start, long long window_end, 
overlap_region_alloc* overlap_list, Correct_dumy* dumy, All_reads* R_INF, Graph* g, Graph* DAGCon, Cigar_record* current_cigar)
{
    clear_Graph(g);
    clear_Graph(DAGCon);

    long long x_start;
    long long x_length; 
    char* x_string;
    char* y_string;
    char* backbone;
    long long backbone_length;
    uint64_t i;
    long long y_start, y_length;
    long long windowID;
    long long startNodeID, endNodeID, currentNodeID;
    overlap_region *z;


    backbone = r_string + window_start;
    backbone_length = window_end + 1 - window_start;

    addUnmatchedSeqToGraph(g, backbone, backbone_length, &startNodeID, &endNodeID);

    
    long long correct_x_pos_s;
    for (i = 0; i < dumy->length; i++)
    {
        // assert(dumy->overlapID[i]<overlap_list->length);
        ///this is the overlap ID
        z = &(overlap_list->list[dumy->overlapID[i]]);

        correct_x_pos_s = (z->x_pos_s / WINDOW) * WINDOW;
        windowID = (window_start - correct_x_pos_s) / WINDOW;
        // assert(windowID<(int64_t)z->w_list.n);

        ///if this window is not matched
        if (z->w_list.a[windowID].y_end == -1) continue;
        
        x_start = z->w_list.a[windowID].x_start;
        x_length = z->w_list.a[windowID].x_end + 1 - z->w_list.a[windowID].x_start;


        y_start = z->w_list.a[windowID].y_start;
        y_length = z->w_list.a[windowID].y_end + 1 - z->w_list.a[windowID].y_start;
        // assert(y_start>=0 && y_start<(int64_t)Get_READ_LENGTH((*R_INF), z->y_id));
        // assert((y_start+y_length)>=0 && (y_start+y_length)<=(int64_t)Get_READ_LENGTH((*R_INF), z->y_id));
        recover_UC_Read_sub_region(dumy->overlap_region, y_start, y_length, z->y_pos_strand, R_INF, z->y_id);

        x_string = r_string + x_start;
        y_string = dumy->overlap_region;
        ///here is the offset of the start base, also is the node ID
        currentNodeID = x_start - window_start;

        // if(window_start == 4500) {
        //     fprintf(stderr, "[M::%s] window_start::%lld, ovlp_id::%lu, windowID::%lld, x_start::%lld, x_length::%lld, y_start::%lld, y_length::%lld, currentNodeID::%lld\n", __func__, 
        //     window_start, dumy->overlapID[i], windowID, x_start, x_length, y_start, y_length, currentNodeID);
        // }
        
        ///cigar: overlap_list->list[overlapID].w_list[windowID].cigar;
        addmatchedSeqToGraph(g, currentNodeID, x_string, x_length, 
                    y_string, y_length, &(z->w_list.a[windowID]), &(z->w_list), startNodeID, endNodeID);
    }

    get_seq_from_Graph(g, DAGCon, dumy, current_cigar, backbone, r_string, r_total_length, window_start);
}


void add_cigar_to_cigar(Correct_dumy* backbone_dumy, Cigar_record* backbone_cigar,
Round2_alignment* second_round,
long long back_bone_start, long long back_bone_length,
long long new_start, long long new_length)
{
    Correct_dumy* new_dumy = &(second_round->dumy);
    Cigar_record* new_cigar = &(second_round->tmp_cigar);
    Cigar_record* result_cigar = &(second_round->cigar);

    char* x_string = backbone_dumy->corrected_read + back_bone_start;
    char* y_string = new_dumy->corrected_read + new_start;
    /**
    if(verify_cigar_2(x_string, back_bone_length, y_string, new_length, new_cigar, -1))
    {
        fprintf(stderr, "error\n");
    }
    **/


   ///if type == 0, x_string here is not useful
   ///output matches to cigar
   add_cigar_record(x_string, back_bone_start - second_round->obtained_cigar_length, result_cigar, 0);
   second_round->obtained_cigar_length = back_bone_start + back_bone_length;

    long long i, cigar_i, x_i, y_i;
    int operation;
    int operationLen;
    x_i = y_i = 0;
    char merge_base;

    for (i = 0; i < (long long)new_cigar->length; i++)
    {
        operation = Get_Cigar_Type(new_cigar->record[i]);
        operationLen = Get_Cigar_Length(new_cigar->record[i]);
        if (operation == 0)
        {
            ///if type == 0, x_string here is not useful
            add_cigar_record(x_string, operationLen, result_cigar, 0);
            x_i += operationLen;
            y_i += operationLen;
        }
        else if (operation == 1)
        {
            for (cigar_i = 0; cigar_i < operationLen; cigar_i++)
            {
                merge_base = 0;
                merge_base = seq_nt6_table[(uint8_t)y_string[y_i]];
                merge_base = merge_base << 3;
                merge_base = merge_base | seq_nt6_table[(uint8_t)x_string[x_i]];

                add_cigar_record(&merge_base, 1, result_cigar, 1);
                x_i++;
                y_i++;
            }
        }
        else if (operation == INSERTION)///2是x缺字符（y多字符）
        {
            add_cigar_record(y_string+y_i, operationLen, result_cigar, INSERTION);
            y_i += operationLen;
        }
        else if (operation == DELETION)
        {
            add_cigar_record(x_string+x_i, operationLen, result_cigar, DELETION);
            x_i += operationLen;
        }
    }
}

///correct bases of current_dumy->corrected_read in [start_base, end_base]
int merge_cigars(Correct_dumy* current_dumy, Cigar_record* current_cigar, 
Round2_alignment* second_round, long long total_start_base, long long total_end_base,
long long total_window_start, long long total_window_end)
{
    Cigar_record* new_cigar = &(second_round->tmp_cigar);

    if(new_cigar->length == 1 && Get_Cigar_Type(new_cigar->record[0]) == 0)
    {
        return 1;
    }

    
    long long start_base = total_start_base - total_window_start;
    long long end_base = total_end_base - total_window_start;
    long long x_i, y_i, cigar_i, i;
    x_i = 0;
    y_i = 0;
    int operation;
    int operationLen;

    long long get_x_start, get_x_end, get_y_start, get_y_end;
    get_x_start = get_x_end = get_y_start = get_y_end = -1;

    int start_cigar = -1;
    int end_cigar = -1;

    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///obtained x_i may larger than start_base/end_base
    ///when operation == 3
    ///so for operation == 3, we need deal with carefully
    for (i = 0; i < (long long)new_cigar->length; i++)
    {
        operation = Get_Cigar_Type(new_cigar->record[i]);
        operationLen = Get_Cigar_Length(new_cigar->record[i]);
        if (operation == 0)
        {
            for (cigar_i = 0; cigar_i < operationLen; cigar_i++)
            {
                if(x_i >= start_base && get_x_start == -1)
                {
                    get_x_start = x_i;
                    get_y_start = y_i;
                    start_cigar = i;
                }



                if(x_i >= end_base && get_x_end == -1)
                {
                    get_x_end = x_i;
                    get_y_end = y_i;
                    end_cigar = i;
                    break;
                }

                

                x_i++;
                y_i++;
            }
        }
        else if (operation == 1)
        {
            for (cigar_i = 0; cigar_i < operationLen; cigar_i++)
            {
                if(x_i >= start_base && get_x_start == -1)
                {
                    get_x_start = x_i;
                    get_y_start = y_i;
                    start_cigar = i;
                }

                if(x_i >= end_base && get_x_end == -1)
                {
                    get_x_end = x_i;
                    get_y_end = y_i;
                    end_cigar = i;
                    break;
                }

                x_i++;
                y_i++;
            }
        }
        else if (operation == 2)
        {
            y_i += operationLen;
        }
        else if (operation == 3)
        {
            ///obtained x_i may larger than start_base/end_base
            ///when operation == 3
            ///so for operation == 3, we need deal with carefully
            x_i += operationLen;
        }
    }

    ///if there are some gap at the end of x, it very likely miscorrection
    if(get_x_end == -1 || get_x_start == -1)
    {
        return 0;
    }
    
    x_i = 0;
    y_i = 0;

    uint32_t single_record = 0;

    for (i = 0; i < (long long)new_cigar->length; i++)
    {
        operation = Get_Cigar_Type(new_cigar->record[i]);
        operationLen = Get_Cigar_Length(new_cigar->record[i]);

        if (i == start_cigar)
        {
            single_record = 0;
            single_record = operationLen - (get_x_start - x_i);
            single_record = single_record << 2;
            single_record = single_record | operation;
            new_cigar->record[i] = single_record;

            if(operation > 1)
            {
                fprintf(stderr, "error\n");
            }

            if (i == end_cigar)
            {
                x_i = get_x_start;
                single_record = 0;
                single_record = get_x_end - x_i + 1;
                single_record = single_record << 2;
                single_record = single_record | operation;
                new_cigar->record[i] = single_record;
                if(operation > 1)
                {
                    fprintf(stderr, "error\n");
                }
                break;
            }
        }
        else if (i == end_cigar)
        {
            single_record = 0;
            single_record = get_x_end - x_i + 1;
            single_record = single_record << 2;
            single_record = single_record | operation;
            new_cigar->record[i] = single_record;
            if(operation > 1)
            {
                fprintf(stderr, "error\n");
            }
            break;
        }


        if (operation == 0 || operation == 1)
        {
            x_i += operationLen;
            y_i += operationLen;
        }
        else if (operation == 2)
        {
            y_i += operationLen;
        }
        else if (operation == 3)
        {
            x_i += operationLen;
        }
    }

    new_cigar->length = end_cigar - start_cigar + 1;
    ///should be improved
    memmove(new_cigar->record, new_cigar->record + start_cigar, new_cigar->length*sizeof(uint32_t));
    
    long long total_x_start = total_window_start + get_x_start;
    long long x_length = get_x_end -get_x_start + 1;
    long long total_y_start = get_y_start;
    long long y_length = get_y_end -get_y_start + 1;

    add_cigar_to_cigar(current_dumy, current_cigar, second_round,
    total_x_start, x_length, total_y_start, y_length);

    return 1;
}

int process_boundary(overlap_region_alloc* overlap_list, All_reads* R_INF, Correct_dumy* dumy, Graph* g, Graph* DAGCon,
Cigar_record* current_cigar, long long uncorrected_window_start, Round2_alignment* second_round, window_list_alloc* win_ciagr_buf)
{
    char* r_string = dumy->corrected_read;
    long long r_total_length = current_cigar->new_read_length;
    long long corrected_window_start, corrected_window_end;
    int extra_begin;
    int extra_end;

    if(dumy->last_boundary_length == 0)
    {
        return 0;
    }

    corrected_window_start = dumy->last_boundary_length - WINDOW_BOUNDARY/2;
    corrected_window_end = dumy->last_boundary_length + WINDOW_BOUNDARY/2 - 1;

    if(corrected_window_start < 0)
    {
        corrected_window_start = 0;
    }

    if (corrected_window_end >= current_cigar->new_read_length)
    {
        corrected_window_end = current_cigar->new_read_length - 1;
    }

    clear_Graph(g);
    clear_Graph(DAGCon);

    long long x_start, x_end;
    long long x_length, x_len, o_len;
    int threshold; 
    long long Window_Len;
    char* x_string = NULL;
    char* y_string = NULL;
    char* backbone = NULL;
    long long backbone_length;
    long long i;
    long long y_start, y_length;
    long long overlapID, windowID;
    long long startNodeID, endNodeID, currentNodeID;
    int end_site;
    unsigned int error;
    int real_y_start;
    long long total_error = 0;

    backbone = r_string + corrected_window_start;
    backbone_length = corrected_window_end - corrected_window_start + 1;
    addUnmatchedSeqToGraph(g, backbone, backbone_length, &startNodeID, &endNodeID);


    long long correct_x_pos_s;
    long long matched_coverage = 0;
    for (i = 0; i < (long long)dumy->length; i++)
    {
        overlapID = dumy->overlapID[i];
        correct_x_pos_s = (overlap_list->list[overlapID].x_pos_s / WINDOW) * WINDOW;
        windowID = (uncorrected_window_start - correct_x_pos_s) / WINDOW;

        ///skip if window is unmatched
        if (overlap_list->list[overlapID].w_list.a[windowID].y_end == -1)
        {
            continue;
        }

        x_start = overlap_list->list[overlapID].w_list.a[windowID].x_start;
        y_start = overlap_list->list[overlapID].w_list.a[windowID].y_start;

        


        /**
         * There are total 3 cases:
         * 1.  this window of x is overlapped fully by y
         *     x: ------|------|---------
         *     y: ------|------|---------
         *     in this case, x_start == uncorrected_window_start, x_length == WINDOW
         * 2.  the suiffx of x's window is overlapped by the prefix of y
         *     x: ------|------|---------
         *               y: |--|-----------  
         *      in this case, x_start > uncorrected_window_start, x_length < WINDOW
         *      this overlap is useless
         * 3.  the prefix of x's window is overlapped by y (see last window)
         *      x: |------|------|-----|---
         *                 y: |--|-----|------  
         *      or
         *          x: |------|------|-----|----
         *        y: --|------|------|-----|--
         *      
         *      in this case, x_start == uncorrected_window_start, x_length < WINDOW
         * 
         *      case 1 and case 3 are useful, while case 2 is useless
         * **/

        ///case 1 and case 3 are useful
        if(x_start == uncorrected_window_start)
        {
            extra_begin = extra_end = 0;
            x_start = corrected_window_start;
            x_end = corrected_window_end;
            x_len = x_end - x_start + 1;
            threshold = x_len * asm_opt.max_ov_diff_ec;
            /****************************may have bugs********************************/
            threshold = Adjust_Threshold(threshold, x_len);
            /****************************may have bugs********************************/
            ///y_start may less than 0
            y_start = y_start - WINDOW_BOUNDARY/2;

            ///in fact, we don't need this line, just worry for bug
            if(y_start < 0)
            {
                continue;
            }

            Window_Len = x_len + (threshold << 1);

            error =(unsigned int)-1;
            if(determine_overlap_region(threshold, y_start, overlap_list->list[overlapID].y_id, Window_Len, Get_READ_LENGTH((*R_INF), overlap_list->list[overlapID].y_id),
            &extra_begin, &extra_end, &y_start, &o_len))
            {
                fill_subregion(dumy->overlap_region, y_start, o_len, overlap_list->list[overlapID].y_pos_strand, 
                R_INF, overlap_list->list[overlapID].y_id, extra_begin, extra_end);

                x_string = r_string + x_start;
                y_string = dumy->overlap_region;

                ///both end site and real_y_start have extra_begin
                ///should be improved, since most of overlaps are exact overlaps
                ///we can do it quickly
                end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
                        &(dumy->path_length), dumy->matrix_bit, dumy->path, -1, -1);
            }

            

            ///try to calculate using higher threshold
            if(error==(unsigned int)-1)
            {
                extra_begin = extra_end = 0;
                x_start = corrected_window_start;
                x_end = corrected_window_end;
                x_len = x_end - x_start + 1;
                threshold = threshold * 2;
                /****************************may have bugs********************************/
                threshold = Adjust_Threshold(threshold, x_len);
                /****************************may have bugs********************************/
                if(x_len >= 300 && threshold < THRESHOLD_MAX_SIZE)
                {
                    threshold = THRESHOLD_MAX_SIZE;
                }
                if(threshold > THRESHOLD_MAX_SIZE)
                {
                    threshold = THRESHOLD_MAX_SIZE;
                }
                Window_Len = x_len + (threshold << 1);
                y_start = overlap_list->list[overlapID].w_list.a[windowID].y_start - WINDOW_BOUNDARY/2;

                ///in fact, we don't need this line, just worry for bug
                if(y_start < 0)
                {
                    continue;
                }

                error =(unsigned int)-1;
                if(determine_overlap_region(threshold, y_start, overlap_list->list[overlapID].y_id, Window_Len, Get_READ_LENGTH((*R_INF), overlap_list->list[overlapID].y_id),
                &extra_begin, &extra_end, &y_start, &o_len))
                {
                    fill_subregion(dumy->overlap_region, y_start, o_len, overlap_list->list[overlapID].y_pos_strand, 
                    R_INF, overlap_list->list[overlapID].y_id, extra_begin, extra_end);

                    x_string = r_string + x_start;
                    y_string = dumy->overlap_region;

                    ///both end site and real_y_start have extra_begin
                    ///should be improved, since most of overlaps are exact overlaps
                    ///we can do it quickly
                    end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
                            &(dumy->path_length), dumy->matrix_bit, dumy->path, -1, -1);
                }
            }

            if (error!=(unsigned int)-1)
            {

                

                total_error = total_error + error;
                matched_coverage++;
                win_ciagr_buf->a[0].x_start = x_start; win_ciagr_buf->a[0].x_end = x_end;
                win_ciagr_buf->c.n = 0;
                generate_cigar(dumy->path, dumy->path_length, &(win_ciagr_buf->a[0]), win_ciagr_buf, &real_y_start, &end_site, &error, x_string, x_len, y_string);  
                ///both end site and real_y_start have extra_begin
                real_y_start -= extra_begin;
                end_site -= extra_begin;

                y_length = end_site - real_y_start + 1;
                y_start = y_start + real_y_start;
                

                x_start = corrected_window_start;
                x_length = corrected_window_end - x_start + 1;

                ///here can be improved, make y_string = dumy->overlap_region + real_y_start + extra_begin
                recover_UC_Read_sub_region(dumy->overlap_region, y_start, y_length, overlap_list->list[overlapID].y_pos_strand, 
                R_INF, overlap_list->list[overlapID].y_id);

                x_string = r_string + x_start;
                y_string = dumy->overlap_region;

                currentNodeID = x_start - corrected_window_start;
                
                addmatchedSeqToGraph(g, currentNodeID, x_string, x_length, y_string, y_length, &(win_ciagr_buf->a[0]), win_ciagr_buf, startNodeID, endNodeID);
            }

        }///case 2 is useless
        else if(x_start != uncorrected_window_start)
        {
            continue;
        }        
    }


    if(matched_coverage >= MIN_COVERAGE_THRESHOLD)
    {
        ///if there are no error, we do not need correction
        if(total_error == 0)
        {
            return 0;
        }
        
        clear_Cigar_record(&(second_round->tmp_cigar));
        clear_Correct_dumy_pure(&(second_round->dumy));

        ///correct bases in [start_base, end_base]
        long long start_base = corrected_window_start + WINDOW_UNCORRECT_SINGLE_SIDE_BOUNDARY;
        long long end_base = corrected_window_end - WINDOW_UNCORRECT_SINGLE_SIDE_BOUNDARY;

        if(end_base > start_base)
        {
            ///note there is an additional "S" node
            ///and start from i-th node, we can correct (i+1)-th base
            /// so the condition when traversing graph is 
            ///(node >= start_base - corrected_window_start && node <= end_base - corrected_window_start)
            get_seq_from_Graph(g, DAGCon, &(second_round->dumy), &(second_round->tmp_cigar), backbone, 
            r_string, r_total_length, corrected_window_start);
            
            /**
            if(verify_cigar_2(backbone, backbone_length, second_round->dumy.corrected_read, 
                second_round->dumy.corrected_read_length, &(second_round->tmp_cigar), -1))
            {
                fprintf(stderr, "hahah\n");
            }
            **/
            
            merge_cigars(dumy, current_cigar, second_round, start_base, end_base,
            corrected_window_start, corrected_window_end);
            


        }
        
    }
    else
    {
        return 0;
    }

    return 1;
    
    
    
}


void generate_consensus(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, Graph* g, Graph* DAGCon, Cigar_record* current_cigar,
                        Round2_alignment* second_round, window_list_alloc* win_ciagr_buf)
{
    clear_Cigar_record(current_cigar);
    
    long long window_start, window_end;

    long long num_availiable_win = 0;
    
    


    Window_Pool w_inf;
    init_Window_Pool(&w_inf, g_read->length, WINDOW, (int)(1.0/asm_opt.max_ov_diff_ec));

    int flag = 0;
    ///for last window
    dumy->last_boundary_length = 0;
    while(get_Window(&w_inf, &window_start, &window_end) && flag != -2)
    {

        

        dumy->length = 0;
        dumy->lengthNT = 0;

        ///return overlaps that are overlaped with [window_start, window_end]
        flag = get_available_interval(window_start, window_end, overlap_list, dumy);
        switch (flag)
        {
            case 1:    ///match
                break;
            case 0:    ///unmatch
                break;
            case -2: ///unmatch, and cannot match for next window
                break;
        }

        num_availiable_win = num_availiable_win + dumy->length;

        // fprintf(stderr, "[M::%s] window_start::%lld, window_end::%lld, g_read->length::%lld, dumy->length::%lu\n", __func__, 
        // window_start, window_end, g_read->length, dumy->length);
        
        ///number of overlaps, also be the coverage
        if(dumy->length >= MIN_COVERAGE_THRESHOLD)
        {
            window_consensus(g_read->seq, g_read->length, window_start, window_end, overlap_list, 
            dumy, R_INF, g, DAGCon, current_cigar);

            if(dumy->last_boundary_length != 0)
            {
                process_boundary(overlap_list, R_INF, dumy, g, DAGCon, current_cigar, 
                window_start, second_round, win_ciagr_buf);
            }
            
        }
        else
        {
            add_segment_to_correct_read(dumy, g_read->seq + window_start, window_end - window_start + 1);
            add_cigar_record(g_read->seq + window_start, window_end - window_start + 1, current_cigar, 0);
        }


        dumy->last_boundary_length = current_cigar->new_read_length;
    }

    if (window_start < g_read->length)
    {
        add_segment_to_correct_read(dumy, g_read->seq + window_start, g_read->length - window_start);
        add_cigar_record(g_read->seq + window_start, g_read->length - window_start, current_cigar, 0);        
    }

    ///if type == 0, x_string here is not useful
    ///output matches to cigar
    if (current_cigar->new_read_length != second_round->obtained_cigar_length)
    {
        add_cigar_record(dumy->corrected_read, current_cigar->new_read_length - second_round->obtained_cigar_length, 
        &(second_round->cigar), 0);
    }
}


inline int get_available_fully_covered_interval(long long window_start, long long window_end, 
overlap_region_alloc* overlap_list, Correct_dumy* dumy, long long* real_length, long long* real_length_100)
{
    long long i, fud = 0;
    long long Len;
    long long overlap_length;

    if(window_start == 0) dumy->start_i = 0;
    for (i = dumy->start_i; i < (long long)overlap_list->length; i++)
    {
        if (window_end < (long long)overlap_list->list[i].x_pos_s)
        {
            dumy->start_i = 0;
            return 0;
        }
        else 
        {
            dumy->start_i = i;
            break;
        }
    }


    if (i >= (long long)overlap_list->length)
    {
        dumy->start_i = overlap_list->length;
        return -2;
    }
    



    
    long long fake_length = 0;
    overlap_length = window_end - window_start + 1;
    (*real_length) = 0; fud = 0;

    for (; i < (long long)overlap_list->length; i++)
    {
        if((Len = OVERLAP(window_start, window_end, (long long)overlap_list->list[i].x_pos_s, (long long)overlap_list->list[i].x_pos_e)) > 0)
        {
            fake_length++;

            if (overlap_length == Len && overlap_list->list[i].is_match == 1)
            {
                (*real_length)++;
            }

            if (overlap_length == Len && overlap_list->list[i].is_match == 100)
            {
                (*real_length_100)++;
            }
            if(fud == 0) fud = 1, dumy->start_i = i;
        }
        
        if((long long)overlap_list->list[i].x_pos_s > window_end)
        {
            break;
        }
    }

    if (fake_length == 0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

int check_if_fully_covered(overlap_region_alloc* overlap_list, 
All_reads* R_INF, UC_Read* g_read, Correct_dumy* dumy, Graph* g, int* abnormal)
{
    long long window_start, window_end;
    int return_flag = 1;
    (*abnormal) = 0;
    
    Window_Pool w_inf;
    init_Window_Pool(&w_inf, g_read->length, WINDOW, (int)(1.0/asm_opt.max_ov_diff_ec));

    int flag = 0;
    long long realLen = 0, tmpLen = 0;
    while(get_Window(&w_inf, &window_start, &window_end) && flag != -2)
    {
        dumy->length = 0;
        dumy->lengthNT = 0;
        
        ///return overlaps that are overlaped with [window_start, window_end]
        flag = get_available_fully_covered_interval(window_start, window_end, 
        overlap_list, dumy, &realLen, &tmpLen);


        switch (flag)
        {
            case 1:    ///match
                break;
            case 0:    ///unmatch
                break;
            case -2: ///unmatch, and cannot match for next window
                break;
        }

        if(realLen < MIN_COVERAGE_THRESHOLD * 2)
        {
            return_flag = 0;
           //return 0;
        }

        if(realLen == 0)
        {
            ///that means this window is a middle window
            if(window_start != 0 && window_end != g_read->length - 1)
            {
                (*abnormal) = 1;
            }
            else if((*abnormal)==0)
            {
                (*abnormal) = 2;
            }
        }
    }

    return return_flag;
}

///mark SNPs at [xBeg, xEnd], note we need to deal with flag_offset carefully
void markSNP_detail(window_list *cigar_idx, window_list_alloc *cigar_s, uint8_t* flag, 
long long xBeg, long long xEnd, long long flag_offset, const ul_idx_t *uref, long long y_total_start, int y_strand, int yid)
{
    if(xBeg > xEnd) return;

    int64_t x_i, y_i, c_i, c_n = cigar_idx->clen, pi = 0, cc = 0;
    uint32_t i, operLen = (uint32_t)-1; uint8_t oper = (uint8_t)-1;
    i = c_i = x_i = y_i = 0;
    for (c_i = 0; c_i < c_n; c_i++) {
        get_cigar_cell(cigar_idx, cigar_s, c_i, &oper, &operLen);
        if(x_i > xEnd) break;

        if (oper == 0) {///match
            x_i += operLen; y_i += operLen;
        } 
        else if(oper == 1) {///mismatch
            for (i = 0; i < operLen; i++) {
                /// note we need to deal with flag_offset carefully
                ///if(flag[x_i - flag_offset] < 127 && x_i >= xBeg && x_i <= xEnd)
                if(x_i >= xBeg && x_i <= xEnd && flag[x_i - flag_offset] < 127) {///Fix-attention
                    if(uref) {
                        cc = retrieve_u_cov(uref, yid, y_strand, y_total_start + y_i, y_strand, &pi);
                        cc += flag[x_i - flag_offset];
                        flag[x_i - flag_offset] = (cc <= 127?cc:127);
                        // if(cc <= 127) flag[x_i - flag_offset] = cc;
                    } else {
                        flag[x_i - flag_offset]++;
                    }
                }
                x_i++; y_i++;
            }
        } else if (oper == 2) {///insertion, that means y has more bases than x
            y_i += operLen;
        }
        else if (oper == 3) {
            x_i += operLen;
        }
    }
}


///window_offset is still the x-based offset
///x_total_start and y_total_start are global positions, instead of local positions
void markSNP_advance(
long long window_offset,
long long x_total_start, long long x_length, 
long long y_total_start, long long y_length, 
window_list *current_cigar, window_list_alloc *current_cigar_s, 
window_list *beg_cigar, window_list_alloc *beg_cigar_s,
window_list *end_cigar, window_list_alloc *end_cigar_s,
haplotype_evdience_alloc* hap, const ul_idx_t *uref, int strand, int yid)
{   
    long long x_total_end = x_total_start + x_length - 1;
    ///mismatches based on the offset of x
    long long inner_offset = x_total_start - window_offset;
    ///long long useless_side, xleftLen, xrightLen, x_interval_beg, x_interval_end;
    long long L_useless_side, R_useless_side, xleftLen, xrightLen, x_interval_beg, x_interval_end;
    long long current_cigar_beg, current_cigar_end;

    ///for current_cigar, [current_cigar_beg, current_cigar_end]
    current_cigar_beg = 0;
    current_cigar_end = x_length - 1;
    ///if the beg_cigar is available
    if(beg_cigar != NULL && beg_cigar->y_end!=-1)
    {
        ///useless_side = beg_cigar->error_threshold;
        L_useless_side = beg_cigar->extra_begin;
        R_useless_side = beg_cigar->extra_end;
        ///again, xleftLen does not include x_total_start itself, but includes beg_cigar->x_start
        xleftLen = x_total_start - beg_cigar->x_start;
        ///xrightLen includes both x_total_start and beg_cigar->x_end
        xrightLen = beg_cigar->x_end - x_total_start + 1;
        ///actually xleftLen could be no larger than useless_side
        ///but such window has already been filtered out at calculate_boundary_cigars
        ///if(xleftLen > useless_side && xrightLen > useless_side)
        if(xleftLen > L_useless_side && xrightLen > R_useless_side)
        {
            ///[x_interval_beg, x_interval_end] are the offsets to beg_cigar->x_start
            ///they are local postions, instead of global positions

            x_interval_beg = xleftLen;
            //x_interval_end = x_interval_beg + (xrightLen - useless_side) - 1;
            x_interval_end = x_interval_beg + (xrightLen - R_useless_side) - 1;
            ///current_cigar_beg is the offset of the current cigar
            ///that is the beg of current_cigar_beg
            ///current_cigar_beg = xrightLen - useless_side;
            current_cigar_beg = xrightLen - R_useless_side;

            markSNP_detail(beg_cigar, beg_cigar_s, hap->flag + inner_offset, x_interval_beg, x_interval_end, 
            x_interval_beg, uref, beg_cigar->y_start, strand, yid);
        }
    }

    if(end_cigar!=NULL && end_cigar->y_end!=-1)
    {
        ///useless_side = end_cigar->error_threshold;
        L_useless_side = end_cigar->extra_begin;
        R_useless_side = end_cigar->extra_end;
        ///again, xleftLen does not include x_total_end, but includes end_cigar->x_start
        ///it seems to be not what we want
        xleftLen = x_total_end - end_cigar->x_start;
        ///xrightLen includes both x_total_end and end_cigar->x_end
        ///it is also not what we want
        xrightLen = end_cigar->x_end - x_total_end + 1;
        ///we hope that x_total_end should be included in xleftLen, instead of xrightLen
        ///that means xleftLen should + 1, while xrightLen should -1
        ///but it is fine here 

        ///actually xrightLen could be no larger than useless_side
        ///but such window has already been filtered out in calculate_boundary_cigars
        ///if(xleftLen > useless_side && xrightLen > useless_side)
        if(xleftLen > L_useless_side && xrightLen > R_useless_side)
        {
            ///[x_interval_beg, x_interval_end] are the offsets to beg_cigar->x_start
            ///they are local postions, instead of global positions
            x_interval_end = xleftLen;
            ///the real left part length is (xleftLen + 1)
            ///so the useful left part length is ((xleftLen + 1) - useless_side)
            ///x_interval_beg = x_interval_end - ((xleftLen + 1) - useless_side) + 1;
            x_interval_beg = x_interval_end - ((xleftLen + 1) - L_useless_side) + 1;
            ///current_cigar_end = (x_length - 1) - ((xleftLen + 1) - useless_side);
            current_cigar_end = (x_length - 1) - ((xleftLen + 1) - L_useless_side);
            
            markSNP_detail(end_cigar, end_cigar_s, hap->flag + end_cigar->x_start - window_offset, x_interval_beg, 
            x_interval_end, 0, uref, end_cigar->y_start, strand, yid);
        }
    }

    markSNP_detail(current_cigar, current_cigar_s, hap->flag + inner_offset, current_cigar_beg, 
    current_cigar_end, 0, uref, current_cigar->y_start, strand, yid);
}



/**
void addSNPtohaplotype(
long long window_offset, int overlapID,
char* x_string, long long x_total_start, long long x_length, 
char* y_string, long long y_total_start, long long y_length, 
CIGAR* cigar, haplotype_evdience_alloc* hap, int snp_threshold)
{
    
    int x_i, y_i, cigar_i;
    x_i = 0;
    y_i = 0;
    cigar_i = 0;
    int operation;
    int operationLen;
    int i;
    long long inner_offset = x_total_start - window_offset;
    haplotype_evdience ev;
    
    ///note that node 0 is the start node
    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2 represents thre are more bases at y
    ///3 represents thre are more bases at x
    while (cigar_i < cigar->length)
    {
        operation = cigar->C_C[cigar_i];
        operationLen = cigar->C_L[cigar_i];

        ///matches
        if (operation == 0)
        {
            for (i = 0; i < operationLen; i++)
            {
                ///should be at least 2 mismatches
                if(hap->flag[inner_offset] > snp_threshold)
                {
                    ev.misBase =  y_string[y_i];
                    ev.overlapID = overlapID;
                    ev.site = x_total_start + x_i;
                    ev.overlapSite = y_total_start + y_i;
                    ev.type = 0;
                    addHaplotypeEvdience(hap, &ev, NULL);
                }


                inner_offset++;
                x_i++;
                y_i++;
            }

        }
        else if(operation == 1)
        {
            for (i = 0; i < operationLen; i++)
            {

                if(hap->flag[inner_offset] > snp_threshold)
                {
                    ev.misBase =  y_string[y_i];
                    ev.overlapID = overlapID;
                    ev.site = x_total_start + x_i;
                    ev.overlapSite = y_total_start + y_i;
                    ev.type = 1;
                    addHaplotypeEvdience(hap, &ev, NULL);
                }

                inner_offset++;
                x_i++;
                y_i++;
            }
        }///insertion
        else if (operation == 2)
        {
            y_i += operationLen;
        }
        else if (operation == 3)
        {
            //may have bugs
            for (i = 0; i < operationLen; i++)
            {
                if(hap->flag[inner_offset] > snp_threshold)
                {
                    ev.misBase =  'N';
                    ev.overlapID = overlapID;
                    ev.site = x_total_start + x_i;
                    ev.overlapSite = y_total_start + y_i;
                    ev.type = 2;
                    addHaplotypeEvdience(hap, &ev, NULL);
                }

                inner_offset++;
                x_i++;
            }
            //may have bugs
        }
        
        cigar_i++;
    }
}
**/


///mark SNPs at [xBeg, xEnd], note we need to deal with flag_offset carefully
void addSNPtohaplotype_details(window_list *cigar_idx, window_list_alloc *cigar_s, uint8_t* flag, 
char* x_string, char* y_string, long long x_total_start, long long y_total_start,
long long xBeg, long long xEnd, int overlapID, long long flag_offset, 
haplotype_evdience_alloc* hap, long long snp_threshold, const ul_idx_t *uref, int y_strand, int yid, void *km)
{
    if(xBeg > xEnd) return;
    int64_t x_i, y_i, c_i, pi = 0, c_n = cigar_idx->clen; 
    uint32_t i, operLen = (uint32_t)-1; uint8_t oper = (uint8_t)-1; i = c_i = x_i = y_i = 0;
    haplotype_evdience ev;

    ///note that node 0 is the start node
    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2 represents thre are more bases at y
    ///3 represents thre are more bases at x
    for (c_i = 0; c_i < c_n; c_i++) {
        get_cigar_cell(cigar_idx, cigar_s, c_i, &oper, &operLen);
        if(x_i > xEnd) break;

        if (oper == 0) { ///matches
            for (i = 0; i < operLen; i++) {
                ///should be at least 2 mismatches
                /// note we need to deal with flag_offset carefully
                ///if(flag[x_i - flag_offset] > snp_threshold && x_i >= xBeg && x_i <= xEnd)
                if(x_i >= xBeg && x_i <= xEnd && flag[x_i - flag_offset] > snp_threshold) {
                    ev.misBase =  y_string[y_i];
                    ev.overlapID = overlapID;
                    ev.site = x_total_start + x_i;
                    ev.overlapSite = y_total_start + y_i;
                    ev.type = 0;
                    ev.cov = uref?retrieve_u_cov(uref, yid, y_strand, y_total_start + y_i, y_strand, &pi):1;
                    addHaplotypeEvdience(hap, &ev, km);
                }
                ///inner_offset++;
                x_i++; y_i++;
            }
        }
        else if(oper == 1) {
            for (i = 0; i < operLen; i++) {
                /// should be at least 2 mismatches
                /// note we need to deal with flag_offset carefully
                ///if(flag[x_i - flag_offset] > snp_threshold && x_i >= xBeg && x_i <= xEnd)
                if(x_i >= xBeg && x_i <= xEnd && flag[x_i - flag_offset] > snp_threshold) {
                    ev.misBase =  y_string[y_i];
                    ev.overlapID = overlapID;
                    ev.site = x_total_start + x_i;
                    ev.overlapSite = y_total_start + y_i;
                    ev.type = 1;
                    ev.cov = uref?retrieve_u_cov(uref, yid, y_strand, y_total_start + y_i, y_strand, &pi):1;
                    addHaplotypeEvdience(hap, &ev, km);
                }
                ///inner_offset++;
                x_i++; y_i++;
            }
        }///insertion, 2 represents thre are more bases at y
        else if (oper == 2) {
            y_i += operLen;
        }///3 represents thre are more bases at x
        else if (oper == 3) {
            /****************************may have bugs********************************/
            for (i = 0; i < operLen; i++) {
                ///if(hap->flag[inner_offset] > snp_threshold)
                /// should be at least 2 mismatches
                /// note we need to deal with flag_offset carefully
                ///if(flag[x_i - flag_offset] > snp_threshold && x_i >= xBeg && x_i <= xEnd)
                if(x_i >= xBeg && x_i <= xEnd && flag[x_i - flag_offset] > snp_threshold) {
                    ev.misBase =  'N';
                    ev.overlapID = overlapID;
                    ev.site = x_total_start + x_i;
                    ev.overlapSite = y_total_start + y_i;
                    ev.type = 2;
                    ev.cov = uref?retrieve_u_cov(uref, yid, y_strand, y_total_start + y_i, y_strand, &pi):1;
                    addHaplotypeEvdience(hap, &ev, km);
                }

                ///inner_offset++;
                x_i++;
            }
            /****************************may have bugs********************************/
        }
    }
}


void addSNPtohaplotype_advance(
long long window_offset, int overlapID,
long long x_total_start, long long x_length, 
long long y_total_start, long long y_length, 
window_list* current_cigar, window_list_alloc *current_cigar_s, 
window_list* beg_cigar, window_list_alloc *beg_cigar_s,
window_list* end_cigar, window_list_alloc *end_cigar_s,
haplotype_evdience_alloc* hap, int snp_threshold, char* x_T_string, char* y_T_string, const ul_idx_t *uref, int strand, int yid, void *km)
{
    long long x_total_end = x_total_start + x_length - 1;
    long long inner_offset = x_total_start - window_offset;
    ///long long useless_side, xleftLen, xrightLen, x_interval_beg, x_interval_end;
    long long L_useless_side, R_useless_side, xleftLen, xrightLen, x_interval_beg, x_interval_end;
    long long current_cigar_beg, current_cigar_end;

    ///for current_cigar, [current_cigar_beg, current_cigar_end]
    current_cigar_beg = 0;
    current_cigar_end = x_length - 1;
    ///if the beg_cigar is available
    if(beg_cigar != NULL && beg_cigar->y_end!=-1)
    {
        ///useless_side = beg_cigar->error_threshold;
        L_useless_side = beg_cigar->extra_begin;
        R_useless_side = beg_cigar->extra_end;
        ///again, xleftLen does not include x_total_start itself, but includes beg_cigar->x_start
        xleftLen = x_total_start - beg_cigar->x_start;
        ///xrightLen includes both x_total_start and beg_cigar->x_end
        xrightLen = beg_cigar->x_end - x_total_start + 1;
        ///actually xleftLen could be no larger than useless_side
        ///but such window has already been filtered out at calculate_boundary_cigars
        ///if(xleftLen > useless_side && xrightLen > useless_side)
        if(xleftLen > L_useless_side && xrightLen > R_useless_side) {
            ///[x_interval_beg, x_interval_end] are the offsets to beg_cigar->x_start
            ///they are local postions, instead of global positions

            x_interval_beg = xleftLen;
            ///x_interval_end = x_interval_beg + (xrightLen - useless_side) - 1;
            x_interval_end = x_interval_beg + (xrightLen - R_useless_side) - 1;
            ///current_cigar_beg is the offset of the current cigar
            ///that is the beg of current_cigar_beg
            ///current_cigar_beg = xrightLen - useless_side;
            current_cigar_beg = xrightLen - R_useless_side;

            // markSNP_detail(cigar_record, hap->flag + inner_offset, x_interval_beg, 
            // x_interval_end, x_interval_beg);
            addSNPtohaplotype_details(beg_cigar, beg_cigar_s, hap->flag + inner_offset, 
            x_T_string + beg_cigar->x_start, y_T_string + beg_cigar->y_start, 
            beg_cigar->x_start, beg_cigar->y_start, x_interval_beg, x_interval_end, 
            overlapID, x_interval_beg, hap, snp_threshold, uref, strand, yid, km);
        }
    }


    if(end_cigar!=NULL && end_cigar->y_end!=-1)
    {
        ///useless_side = end_cigar->error_threshold;
        L_useless_side = end_cigar->extra_begin;
        R_useless_side = end_cigar->extra_end;
        ///again, xleftLen does not include x_total_end, but includes end_cigar->x_start
        ///it seems to be not what we want
        xleftLen = x_total_end - end_cigar->x_start;
        ///xrightLen includes both x_total_end and end_cigar->x_end
        ///it is also not what we want
        xrightLen = end_cigar->x_end - x_total_end + 1;
        ///we hope that x_total_end should be included in xleftLen, instead of xrightLen
        ///that means xleftLen should + 1, while xrightLen should -1
        ///but it is fine here 

        ///actually xrightLen could be no larger than useless_side
        ///but such window has already been filtered out in calculate_boundary_cigars
        ///if(xleftLen > useless_side && xrightLen > useless_side)
        if(xleftLen > L_useless_side && xrightLen > R_useless_side)
        {
            ///[x_interval_beg, x_interval_end] are the offsets to beg_cigar->x_start
            ///they are local postions, instead of global positions
            x_interval_end = xleftLen;
            ///the real left part length is (xleftLen + 1)
            ///so the useful left part length is ((xleftLen + 1) - useless_side)
            ///x_interval_beg = x_interval_end - ((xleftLen + 1) - useless_side) + 1;
            x_interval_beg = x_interval_end - ((xleftLen + 1) - L_useless_side) + 1;
            ///current_cigar_end = (x_length - 1) - ((xleftLen + 1) - useless_side);
            current_cigar_end = (x_length - 1) - ((xleftLen + 1) - L_useless_side);
            
            // markSNP_detail(cigar_record, hap->flag + end_cigar->x_start - window_offset,
            // x_interval_beg, x_interval_end, 0);
            addSNPtohaplotype_details(end_cigar, end_cigar_s, hap->flag + end_cigar->x_start - window_offset,
            x_T_string + end_cigar->x_start, y_T_string + end_cigar->y_start, 
            end_cigar->x_start, end_cigar->y_start, x_interval_beg, x_interval_end, 
            overlapID, 0, hap, snp_threshold, uref, strand, yid, km);
        }
    }

    // markSNP_detail(&(current_cigar->cigar), hap->flag + inner_offset, current_cigar_beg, 
    // current_cigar_end, 0);
    addSNPtohaplotype_details(current_cigar, current_cigar_s, hap->flag + inner_offset, 
    x_T_string + current_cigar->x_start, y_T_string + current_cigar->y_start, 
    current_cigar->x_start, current_cigar->y_start, current_cigar_beg, 
    current_cigar_end, overlapID, 0, hap, snp_threshold, uref, strand, yid, km);
}

/**
void cluster(char* r_string, long long window_start, long long window_end, 
overlap_region_alloc* overlap_list, Correct_dumy* dumy, All_reads* R_INF, haplotype_evdience_alloc* hap)
{
    ///window_start, window_end, and useful_length correspond to x, instead of y
    long long useful_length = window_end - window_start + 1;
    long long x_start;
    long long x_length; 
    char* x_string;
    char* y_string;
    long long i;
    long long y_start, y_length;
    long long overlapID, windowID;

    long long correct_x_pos_s;
    int snp_threshold;
    snp_threshold = 1;

    ///all overlaps related to the current window [window_start, window_end]
    ///first mark all snp pos
    for (i = 0; i < (long long)dumy->length; i++)
    {
        ///overlap id, instead of the window id or the y id
        overlapID = dumy->overlapID[i];

        ///overlap_list->list[overlapID].x_pos_s is the begining of the whole overlap
        correct_x_pos_s = (overlap_list->list[overlapID].x_pos_s / WINDOW) * WINDOW;
        ///window_start is the begining of this window in the whole x_read
        windowID = (window_start - correct_x_pos_s) / WINDOW;

        ///skip if this window is not matched
        if (overlap_list->list[overlapID].w_list[windowID].y_end == -1)
        {
            continue;
        }
        
        ///both x_start and y_start are the offsets of the whole x_read and y_read
        ///instead of the offsets of window
        x_start = overlap_list->list[overlapID].w_list[windowID].x_start;
        x_length = overlap_list->list[overlapID].w_list[windowID].x_end 
                - overlap_list->list[overlapID].w_list[windowID].x_start + 1;

        y_start = overlap_list->list[overlapID].w_list[windowID].y_start;
        y_length = overlap_list->list[overlapID].w_list[windowID].y_end
                - overlap_list->list[overlapID].w_list[windowID].y_start + 1;

            
        markSNP(window_start, x_start, x_length, y_start, y_length, 
        &(overlap_list->list[overlapID].w_list[windowID].cigar), hap);
    }


    //may have bugs
    long long last_snp = -1;
    long long first_snp = -1;
    for (i = 0; i < useful_length; i++)
    {
        if(hap->flag[i] != 0)
        {
            last_snp = i;
            if(first_snp == -1)
            {
                first_snp = i;
            }
        }
        ///for a real snp, the coverage should be at least 2
        if(hap->flag[i] > snp_threshold)
        {
            // hap->snp++;
            hap->nn_snp++;
        }
    }
    ///if there are any >0 elements, both first_snp and last_snp should be != -1
    if(first_snp == -1 || last_snp == -1)
    {
        first_snp = 0;
        last_snp = -1;
    }
    //may have bugs



    ///add the information related to snp to haplotype_evdience_alloc
    for (i = 0; i < (long long)dumy->length; i++)
    {
        ///overlap ID, instead of the window ID
        overlapID = dumy->overlapID[i];

        ///overlap_list->list[overlapID].x_pos_s is the begining of the whole overlap
        correct_x_pos_s = (overlap_list->list[overlapID].x_pos_s / WINDOW) * WINDOW;
        ///window_start is the begining of this window in the whole x_read
        windowID = (window_start - correct_x_pos_s) / WINDOW;

        ///skip if this window is not matched
        if (overlap_list->list[overlapID].w_list[windowID].y_end == -1)
        {
            continue;
        }
        
        ///both x_start and y_start are the offsets of the whole x_read and y_read
        ///instead of the offsets of window
        x_start = overlap_list->list[overlapID].w_list[windowID].x_start;
        x_length = overlap_list->list[overlapID].w_list[windowID].x_end 
                - overlap_list->list[overlapID].w_list[windowID].x_start + 1;

        y_start = overlap_list->list[overlapID].w_list[windowID].y_start;
        y_length = overlap_list->list[overlapID].w_list[windowID].y_end
                - overlap_list->list[overlapID].w_list[windowID].y_start + 1;


        recover_UC_Read_sub_region(dumy->overlap_region, y_start, y_length, overlap_list->list[overlapID].y_pos_strand, 
                R_INF, overlap_list->list[overlapID].y_id);

        x_string = r_string + x_start;
        y_string = dumy->overlap_region;


        addSNPtohaplotype(window_start, overlapID, x_string, x_start, x_length, 
        y_string, y_start, y_length, &(overlap_list->list[overlapID].w_list[windowID].cigar), 
        hap, snp_threshold);        
    }

    RsetInitHaplotypeEvdienceFlag(hap, first_snp, last_snp + 1 - first_snp);
}
**/


void get_related_cigars(window_list_alloc* boundary_cigars, long long id, window_list** beg_cigar, 
window_list** end_cigar)
{
    (*beg_cigar) = &(boundary_cigars->a[id*2]);
    (*end_cigar) = &(boundary_cigars->a[id*2+1]);
}

int cmp_haplotype_evdience(const void * a, const void * b)
{
    if ((*(haplotype_evdience*)a).site != (*(haplotype_evdience*)b).site)
    {
        return (*(haplotype_evdience*)a).site > (*(haplotype_evdience*)b).site ? 1 : -1; 
    }
    else
    {
        if ((*(haplotype_evdience*)a).type != (*(haplotype_evdience*)b).type)
        {
            return (*(haplotype_evdience*)a).type > (*(haplotype_evdience*)b).type ? 1 : -1; 
        }
        else
        {
            if ((*(haplotype_evdience*)a).misBase != (*(haplotype_evdience*)b).misBase)
            {
                return (*(haplotype_evdience*)a).misBase > (*(haplotype_evdience*)b).misBase ? 1 : -1; 
            }
            else
            {
                return 0;
            }
            
        }
        
    }
    
    
}

void cluster_advance(char* r_string, long long window_start, long long window_end, 
overlap_region_alloc* overlap_list, Correct_dumy* dumy, All_reads* R_INF, 
haplotype_evdience_alloc* hap, UC_Read* overlap_read, int snp_threshold)
{
    window_list* beg_cigar;
    window_list* end_cigar;
    ///window_start, window_end, and useful_length correspond to x, instead of y
    long long useful_length = window_end - window_start + 1;
    long long x_start, x_length, ll = hap->length, lr;
    char* x_string;
    char* y_string;
    long long i;
    long long y_start, y_length;
    long long overlapID, windowID;

    long long correct_x_pos_s;


    ///all overlaps related to the current window [window_start, window_end]
    ///first mark all snp pos
    for (i = 0; i < (long long)dumy->length; i++)
    {
        ///overlap id, instead of the window id or the y id
        overlapID = dumy->overlapID[i];

        ///overlap_list->list[overlapID].x_pos_s is the begining of the whole overlap
        correct_x_pos_s = (overlap_list->list[overlapID].x_pos_s / WINDOW) * WINDOW;
        ///window_start is the begining of this window in the whole x_read
        windowID = (window_start - correct_x_pos_s) / WINDOW;

        ///skip if this window is not matched
        if (overlap_list->list[overlapID].w_list.a[windowID].y_end == -1)
        {
            continue;
        }
        
        ///both x_start and y_start are the offsets of the whole x_read and y_read
        ///instead of the offsets of window
        x_start = overlap_list->list[overlapID].w_list.a[windowID].x_start;
        x_length = overlap_list->list[overlapID].w_list.a[windowID].x_end + 1
                - overlap_list->list[overlapID].w_list.a[windowID].x_start;

        y_start = overlap_list->list[overlapID].w_list.a[windowID].y_start;
        y_length = overlap_list->list[overlapID].w_list.a[windowID].y_end + 1
                - overlap_list->list[overlapID].w_list.a[windowID].y_start;

        beg_cigar = end_cigar = NULL;
        
        if(windowID >= 1)
        {
            beg_cigar = &(overlap_list->list[overlapID].boundary_cigars.a[windowID-1]);
        }

        if(windowID < (long long)(overlap_list->list[overlapID].w_list.n) - 1)
        {
            end_cigar = &(overlap_list->list[overlapID].boundary_cigars.a[windowID]);
        }
            

        markSNP_advance(window_start, x_start, x_length, y_start, y_length, 
        &(overlap_list->list[overlapID].w_list.a[windowID]), &(overlap_list->list[overlapID].w_list),
        beg_cigar, &(overlap_list->list[overlapID].boundary_cigars),
        end_cigar, &(overlap_list->list[overlapID].boundary_cigars), hap, NULL, 
        overlap_list->list[overlapID].y_pos_strand, overlap_list->list[overlapID].y_id);
    }


    /****************************may have bugs********************************/
    long long last_snp = -1;
    long long first_snp = -1;
    for (i = 0, lr = 0; i < useful_length; i++)
    {
        if(hap->flag[i] != 0)
        {
            last_snp = i;
            if(first_snp == -1)
            {
                first_snp = i;
            }
        }
        ///for a real snp, the coverage should be at least 2
        if(hap->flag[i] > snp_threshold)
        {
            // hap->snp++; 
            hap->nn_snp++;
            lr++;
        }
    }
    ///if there are any >0 elements, both first_snp and last_snp should be != -1
    if(first_snp == -1 || last_snp == -1)
    {
        first_snp = 0;
        last_snp = -1;
    }
    /****************************may have bugs********************************/



    ///add the information related to snp to haplotype_evdience_alloc
    for (i = 0; i < (long long)dumy->length; i++)
    {
        ///overlap ID, instead of the window ID
        overlapID = dumy->overlapID[i];

        ///overlap_list->list[overlapID].x_pos_s is the begining of the whole overlap
        correct_x_pos_s = (overlap_list->list[overlapID].x_pos_s / WINDOW) * WINDOW;
        ///window_start is the begining of this window in the whole x_read
        windowID = (window_start - correct_x_pos_s) / WINDOW;

        ///skip if this window is not matched
        if (overlap_list->list[overlapID].w_list.a[windowID].y_end == -1)
        {
            continue;
        }
        
        ///both x_start and y_start are the offsets of the whole x_read and y_read
        ///instead of the offsets of window
        x_start = overlap_list->list[overlapID].w_list.a[windowID].x_start;
        x_length = overlap_list->list[overlapID].w_list.a[windowID].x_end 
                - overlap_list->list[overlapID].w_list.a[windowID].x_start + 1;

        y_start = overlap_list->list[overlapID].w_list.a[windowID].y_start;
        y_length = overlap_list->list[overlapID].w_list.a[windowID].y_end
                - overlap_list->list[overlapID].w_list.a[windowID].y_start + 1;


        if(overlap_list->list[overlapID].y_pos_strand == 0)
        {
            recover_UC_Read(overlap_read, R_INF, overlap_list->list[overlapID].y_id);
        }
        else
        {
            recover_UC_Read_RC(overlap_read, R_INF, overlap_list->list[overlapID].y_id);
        }
        
        x_string = r_string;
        y_string = overlap_read->seq;


        beg_cigar = end_cigar = NULL;
        if(windowID >= 1)
        {
            beg_cigar = &(overlap_list->list[overlapID].boundary_cigars.a[windowID-1]);
        }
        if(windowID < (long long)(overlap_list->list[overlapID].w_list.n) - 1)
        {
            end_cigar = &(overlap_list->list[overlapID].boundary_cigars.a[windowID]);
        }


        addSNPtohaplotype_advance(window_start, overlapID, x_start, x_length, y_start, y_length, 
        &(overlap_list->list[overlapID].w_list.a[windowID]), &(overlap_list->list[overlapID].w_list),
        beg_cigar, &(overlap_list->list[overlapID].boundary_cigars), 
        end_cigar, &(overlap_list->list[overlapID].boundary_cigars), 
        hap, snp_threshold, x_string, y_string, NULL, overlap_list->list[overlapID].y_pos_strand, 
        overlap_list->list[overlapID].y_id, NULL);                
    }

    RsetInitHaplotypeEvdienceFlag(hap, first_snp, last_snp + 1 - first_snp);

    if(hap->length - ll > 1 && lr > 1) radix_sort_haplotype_evdience_srt(hap->list+ll, hap->list + hap->length);
}

void cluster_ul_advance(char* r_string, long long window_start, long long window_end, 
overlap_region_alloc* overlap_list, Correct_dumy* dumy, const ul_idx_t *uref, 
haplotype_evdience_alloc* hap, UC_Read* overlap_read, int snp_threshold, long long blockLen, void *km)
{
    window_list* beg_cigar;
    window_list* end_cigar;
    ///window_start, window_end, and useful_length correspond to x, instead of y
    long long useful_length = window_end - window_start + 1;
    long long x_start, x_length, ll = hap->length, lr;
    char* x_string;
    char* y_string;
    long long i;
    long long y_start, y_length;
    long long overlapID, windowID;

    long long correct_x_pos_s;


    ///all overlaps related to the current window [window_start, window_end]
    ///first mark all snp pos
    for (i = 0; i < (long long)dumy->length; i++)
    {
        ///overlap id, instead of the window id or the y id
        overlapID = dumy->overlapID[i];

        ///overlap_list->list[overlapID].x_pos_s is the begining of the whole overlap
        correct_x_pos_s = (overlap_list->list[overlapID].x_pos_s / blockLen) * blockLen;
        ///window_start is the begining of this window in the whole x_read
        windowID = (window_start - correct_x_pos_s) / blockLen;

        ///skip if this window is not matched
        if (overlap_list->list[overlapID].w_list.a[windowID].y_end == -1)
        {
            continue;
        }
        
        ///both x_start and y_start are the offsets of the whole x_read and y_read
        ///instead of the offsets of window
        x_start = overlap_list->list[overlapID].w_list.a[windowID].x_start;
        x_length = overlap_list->list[overlapID].w_list.a[windowID].x_end 
                - overlap_list->list[overlapID].w_list.a[windowID].x_start + 1;

        y_start = overlap_list->list[overlapID].w_list.a[windowID].y_start;
        y_length = overlap_list->list[overlapID].w_list.a[windowID].y_end
                - overlap_list->list[overlapID].w_list.a[windowID].y_start + 1;

        beg_cigar = end_cigar = NULL;
        
        if(windowID >= 1)
        {
            beg_cigar = &(overlap_list->list[overlapID].boundary_cigars.a[windowID-1]);
        }

        if(windowID < (long long)(overlap_list->list[overlapID].w_list.n) - 1)
        {
            end_cigar = &(overlap_list->list[overlapID].boundary_cigars.a[windowID]);
        }
            

        markSNP_advance(window_start, x_start, x_length, y_start, y_length, 
        &(overlap_list->list[overlapID].w_list.a[windowID]), &(overlap_list->list[overlapID].w_list),
        beg_cigar, &(overlap_list->list[overlapID].boundary_cigars),
        end_cigar, &(overlap_list->list[overlapID].boundary_cigars), 
        hap, uref, overlap_list->list[overlapID].y_pos_strand, overlap_list->list[overlapID].y_id);
    }


    /****************************may have bugs********************************/
    long long last_snp = -1;
    long long first_snp = -1;
    for (i = 0, lr = 0; i < useful_length; i++)
    {
        if(hap->flag[i] != 0)
        {
            last_snp = i;
            if(first_snp == -1)
            {
                first_snp = i;
            }
        }
        ///for a real snp, the coverage should be at least 2
        if(hap->flag[i] > snp_threshold)
        {
            // hap->snp++; 
            hap->nn_snp++;
            lr++;
        }
    }
    ///if there are any >0 elements, both first_snp and last_snp should be != -1
    if(first_snp == -1 || last_snp == -1)
    {
        first_snp = 0;
        last_snp = -1;
    }
    /****************************may have bugs********************************/



    ///add the information related to snp to haplotype_evdience_alloc
    for (i = 0; i < (long long)dumy->length; i++)
    {
        ///overlap ID, instead of the window ID
        overlapID = dumy->overlapID[i];

        ///overlap_list->list[overlapID].x_pos_s is the begining of the whole overlap
        correct_x_pos_s = (overlap_list->list[overlapID].x_pos_s / blockLen) * blockLen;
        ///window_start is the begining of this window in the whole x_read
        windowID = (window_start - correct_x_pos_s) / blockLen;

        ///skip if this window is not matched
        if (overlap_list->list[overlapID].w_list.a[windowID].y_end == -1)
        {
            continue;
        }
        
        ///both x_start and y_start are the offsets of the whole x_read and y_read
        ///instead of the offsets of window
        x_start = overlap_list->list[overlapID].w_list.a[windowID].x_start;
        x_length = overlap_list->list[overlapID].w_list.a[windowID].x_end 
                - overlap_list->list[overlapID].w_list.a[windowID].x_start + 1;

        y_start = overlap_list->list[overlapID].w_list.a[windowID].y_start;
        y_length = overlap_list->list[overlapID].w_list.a[windowID].y_end
                - overlap_list->list[overlapID].w_list.a[windowID].y_start + 1;

        retrieve_u_seq(overlap_read, NULL, &uref->ug->u.a[overlap_list->list[overlapID].y_id], 
        overlap_list->list[overlapID].y_pos_strand, 0, -1, km);
        
        
        x_string = r_string;
        y_string = overlap_read->seq;


        beg_cigar = end_cigar = NULL;
        if(windowID >= 1)
        {
            beg_cigar = &(overlap_list->list[overlapID].boundary_cigars.a[windowID-1]);
        }
        if(windowID < (long long)(overlap_list->list[overlapID].w_list.n) - 1)
        {
            end_cigar = &(overlap_list->list[overlapID].boundary_cigars.a[windowID]);
        }


        addSNPtohaplotype_advance(window_start, overlapID, x_start, x_length, y_start, y_length, 
        &(overlap_list->list[overlapID].w_list.a[windowID]), &(overlap_list->list[overlapID].w_list),
        beg_cigar, &(overlap_list->list[overlapID].boundary_cigars),
        end_cigar, &(overlap_list->list[overlapID].boundary_cigars),
        hap, snp_threshold, x_string, y_string, uref, overlap_list->list[overlapID].y_pos_strand, overlap_list->list[overlapID].y_id, km);                
    }

    RsetInitHaplotypeEvdienceFlag(hap, first_snp, last_snp + 1 - first_snp);

    if(hap->length - ll > 1 && lr > 1) radix_sort_haplotype_evdience_srt(hap->list+ll, hap->list + hap->length);
}


int cmp_snp_stats(const void * a, const void * b)
{
    if ((*(SnpStats*)a).score != (*(SnpStats*)b).score)
    {
        return (*(SnpStats*)a).score < (*(SnpStats*)b).score ? 1 : -1; 
    }
    else
    {
        if ((*(SnpStats*)a).occ_2 != (*(SnpStats*)b).occ_2)
        {
            return (*(SnpStats*)a).occ_2 > (*(SnpStats*)b).occ_2 ? 1 : -1; 
        }
        else
        {
            return 0;
        }   
    }
}


int cmp_max_DP(const void * a, const void * b)
{
    if(Get_Max_DP_Value((*(uint64_t*)a))!=Get_Max_DP_Value((*(uint64_t*)b)))
    {
        return Get_Max_DP_Value((*(uint64_t*)a)) < Get_Max_DP_Value((*(uint64_t*)b))? 1 : -1;
    }
    else
    {
        return 0;
    }
}


int split_sub_list(haplotype_evdience_alloc* hap, 
haplotype_evdience* sub_list, long long sub_length, 
overlap_region_alloc* overlap_list, All_reads* R_INF, UC_Read* g_read)
{
    long long i = 0;
    long long occ_0 = 0;
    long long occ_1 = 0;
    long long occ_1_array[5];
    memset(occ_1_array, 0, sizeof(long long) * 5);
    long long occ_2 = 0;


    for (i = 0; i < sub_length; i++)
    {
        if(sub_list[i].type == 0)
        {
            occ_0++;
        }
        else if(sub_list[i].type == 1)
        {
            occ_1_array[seq_nt6_table[(uint8_t)(sub_list[i].misBase)]]++;
            occ_1++;
        }
        else if(sub_list[i].type == 2)
        {
            occ_2++;
        }
    }

    

    /**
     1. if occ_0 = 0, that means all overlaps are different with this read at this site
     2. it is not possible that occ_1 = 0,
     3. if occ_1 = 1, there are only one difference. It must be a sequencing error.
        (for repeat, it maybe a snp at repeat. but ...)
    **/
    ///if(occ_0 == 0 || occ_1 <= 1)
    if(occ_0 == 0 || occ_1 == 0)
    {
        return 0;
    }

    

    ///note: if the max value except type0 is type2
    ///that means this is no snp hapolyte
    long long max = occ_2;
    long long max_i = -1;

    for (i = 0; i < 5; i++)
    {
        if(occ_1_array[i] > max)
        {
            max = occ_1_array[i];
            max_i = i;
        }
    }



    if(max_i == -1)
    {
        return 0;
    }
    
    if(max <= 1)
    {
        return 0;
    }


    ///if we have two max
    for (i = 0; i < 5; i++)
    {
        if(occ_1_array[i] == max && i != max_i)
        {
            return 0;
        }
    }

    long long new_0 = occ_0 + 1;
    long long new_total = sub_length + 1;
    ///note: here occ_0++ since the read itself has a type0
    double available = new_0 + max;
    double threshold = 0.95;
    available = available/((double)(new_total));
    if(available < threshold)///Fix-attention: looks definitely wrong 
    {
        return 0;
    }

    ///new_total is the number of errors here
    new_total = new_total - new_0;
    ///available is the number of selected errors here
    available = max;
    threshold = 0.70;
    available = available/((double)(new_total));
    if(available < threshold)///Fix-attention: looks definitely wrong 
    {
        return 0;
    }
    
    InsertSNPVector(hap, sub_list, sub_length, s_H[max_i], g_read);

    return 1;
}


int calculate_distance_snp_vector(int8_t *vector1, int8_t *vector2, int Len)
{
    int i;
    for (i = 0; i < Len; i++)
    {
        if(vector1[i] != vector2[i])
        {
            if ((vector1[i] == 0 || vector1[i] == 1) && (vector2[i] == 0 || vector2[i] == 1))
            {
                return 1;
            }
        }
    }
    
    return 0;
}

void print_core_snp(haplotype_evdience_alloc* hap)
{
    uint64_t i, j;
    for (i = 0; i < hap->core_snp; i++)
    {
        fprintf(stderr, "core(i): %lu, site: %u, occ_0: %u, occ_1: %u, occ_2: %u, score: %d\n", 
    (unsigned long)i, hap->snp_stat.a[i].site, hap->snp_stat.a[i].occ_0, hap->snp_stat.a[i].occ_1,
    hap->snp_stat.a[i].occ_2,
    hap->snp_stat.a[i].score); 

        int vectorID = hap->snp_stat.a[i].id;
        int8_t* vector = Get_SNP_Vector((*hap), vectorID);

        for (j = 0; j < hap->overlap; j++)
        {
            if(vector[j] == 0)
            {
                fprintf(stderr, "type: %d, ID: %lu\n", vector[j], (unsigned long)j);
            }
        }


        for (j = 0; j < hap->overlap; j++)
        {
            if(vector[j] == 1)
            {
                fprintf(stderr, "type: %d, ID: %lu\n", vector[j], (unsigned long)j);
            }
        }


        for (j = 0; j < hap->overlap; j++)
        {
            if(vector[j] == 2)
            {
                fprintf(stderr, "type: %d, ID: %lu\n", vector[j], (unsigned long)j);
            }
        }
        
    }
}



void add_to_result_snp_vector(haplotype_evdience_alloc* hap, int8_t *new_vector, int Len)
{
    int8_t *r_vector = Get_Result_SNP_Vector((*hap));
    int j;
    for (j = 0; j < Len; j++)
    {
        if(r_vector[j] == -1)
        {
            if(new_vector[j] == 0)
            {
                hap->result_stat.occ_0++;
                r_vector[j] = new_vector[j];
            }
            else if(new_vector[j] == 1)
            {
                hap->result_stat.occ_1++;
                r_vector[j] = new_vector[j];
            }
        }
        ///can debug here
    }

    hap->result_stat.overlap_num = hap->result_stat.occ_0 + hap->result_stat.occ_1;
}


int debug_add_to_result_snp_vector(haplotype_evdience_alloc* hap, int8_t *new_vector, int Len)
{
    int8_t *r_vector = Get_Result_SNP_Vector((*hap));
    int j;
    for (j = 0; j < Len; j++)
    {
        if(r_vector[j] == -1)
        {
            if(new_vector[j] == 0)
            {
                hap->result_stat.occ_0++;
                r_vector[j] = new_vector[j];
            }
            else if(new_vector[j] == 1)
            {
                hap->result_stat.occ_1++;
                r_vector[j] = new_vector[j];
            }
        }
        else ///can debug here
        {
            ///if((new_vector[j] != -1 && new_vector[j] != 2 && new_vector[j] != r_vector[j]))
            if((new_vector[j] == 0 || new_vector[j] == 1) && new_vector[j] != r_vector[j])
            {
                return j;
            }
        }
    }

    return -1;

}



int merge_snp_vectors_and_test(haplotype_evdience_alloc* hap, int diff_vector_ID)
{
    int8_t *r_vector = Get_Result_SNP_Vector((*hap));
    int vectorLen = Get_SNP_Vector_Length((*hap));
    memset(r_vector, -1, vectorLen);
    hap->result_stat.occ_0 = 0;
    hap->result_stat.occ_1 = 0;

    int8_t* vector;
    int vectorID;
    int i, j;

    for (i = 0; i < (int)hap->core_snp; i++)
    {
        if(i == diff_vector_ID)
        {
            continue;
        }

        vectorID = hap->snp_stat.a[i].id;
        vector = Get_SNP_Vector((*hap), vectorID);

        for (j = 0; j < vectorLen; j++)
        {
            if(r_vector[j] == -1)
            {
                if(vector[j] == 0)
                {
                    hap->result_stat.occ_0++;
                    r_vector[j] = vector[j];
                }
                else if(vector[j] == 1)
                {
                    hap->result_stat.occ_1++;
                    r_vector[j] = vector[j];
                }
            }
            else ///can debug here
            {
                ///has confilict
                if(vector[j] != -1 && vector[j] != 2 && vector[j] != r_vector[j])
                {
                    return 0;
                }
            }
        }
    }

    hap->result_stat.overlap_num = hap->result_stat.occ_0 + hap->result_stat.occ_1;

    return 1;
}

int generate_haplotypes(haplotype_evdience_alloc* hap)
{
    int j;

    int vectorID2;
    int8_t *vector, *vector2;
    
    if(hap->core_snp == 0)
    {
        return 0;
    }

    ///sort by weight
    // qsort(hap->snp_stat, hap->available_snp, sizeof(SnpStats), cmp_snp_stats);
    qsort(hap->snp_stat.a, hap->snp_stat.n, sizeof(SnpStats), cmp_snp_stats);


    ///the hap->core_snp is used to find centriod    
    ///if there are <5 vectors in core_snp, we didn't allow different vector
    if (hap->core_snp < 5)
    {
        if(merge_snp_vectors_and_test(hap, -1) == 0)
        {
            return 0;
        }
    }
    else ///for vectors in core_snp, we allow at most one different vector when there are >= 5 vectors in core_snp
    {
        ///there are two condition: 1. vector 0 is the different one. 2. vector 0 is not the different one        
        ///first try to merge all vector together
        if(merge_snp_vectors_and_test(hap, -1) == 0)
        {
            for (j = hap->core_snp - 1; j >= 0; j--)
            {
                if(merge_snp_vectors_and_test(hap, j) == 1)
                {
                    break;
                }
            }

            if(j == -1)
            {
                return 0;
            }
        }
    }


    

    ///after merge, we get result vector
    vector = Get_Result_SNP_Vector((*hap));
    ///and for each non-core snp vector, if it has no conflict with result vector
    /// add it to result vector
    for (j = hap->core_snp; j < (int)hap->snp_stat.n/**hap->available_snp**/; j++)
    {
        vectorID2 = hap->snp_stat.a[j].id;
        vector2 = Get_SNP_Vector((*hap), vectorID2);
        if(calculate_distance_snp_vector(vector, vector2, Get_SNP_Vector_Length((*hap))) == 0)
        {
            add_to_result_snp_vector(hap, vector2, Get_SNP_Vector_Length((*hap)));
        }
    }


    ///for read only have 1 snp, we need a more strict condition
    if (hap->core_snp == 1 && 
    filter_one_snp(hap->result_stat.occ_0 + 1, hap->result_stat.occ_1, 
    hap->result_stat.overlap_num + 1) == 0)
    {
        return 0;
    }


    return 1;
    
}

void Preorder_Merge(uint32_t snpID, haplotype_evdience_alloc* hap, int is_merge) 
{ 
    int vectorID = hap->snp_stat.a[snpID].id;
    int8_t* vector = Get_SNP_Vector((*hap), vectorID);
    hap->dp.visit[snpID] = 1;
    

    if(is_merge)
    {
        if(hap->snp_stat.a[snpID].is_homopolymer)
        {
            hap->result_stat.homopolymer_num++;
        }
        else
        {
            hap->result_stat.non_homopolymer_num++;
        }
        
        hap->result_stat.score++;
        int flag;
        if((flag = debug_add_to_result_snp_vector(hap, vector, Get_SNP_Vector_Length((*hap))))!= -1)
        {
            fprintf(stderr, "incompatible snp vector....\n");
            exit(0);
        }
    }

    uint32_t* column;
    int j;

    if(hap->dp.backtrack_length[snpID] != 0)
    {
        column = Get_DP_Backtrack_Column(hap->dp, snpID);
        
        if(is_merge)
        {
            int add_ID = 0;
            for (j = 0; j < (int)hap->dp.backtrack_length[snpID]; j++)
            {
                if(hap->snp_stat.a[column[j]].is_homopolymer == 0)
                {
                    add_ID = j;
                }
            }

            for (j = 0; j < (int)hap->dp.backtrack_length[snpID]; j++)
            {
                if(j == add_ID)
                {
                    Preorder_Merge(column[j], hap, 1); 
                }
                else
                {
                    Preorder_Merge(column[j], hap, 0); 
                }
            }
        }
        else
        {
            for (j = 0; j < (int)hap->dp.backtrack_length[snpID]; j++)
            {
                Preorder_Merge(column[j], hap, 0); 
            }
        }
    }
}


void Preorder_Merge_Advance_Repeat(uint32_t snpID, haplotype_evdience_alloc* hap, int pathLen) 
{ 
    hap->dp.visit[snpID] = 1;
    hap->dp.buffer[pathLen] = snpID;
    pathLen++;

    if(hap->dp.backtrack_length[snpID] == 0)
    {
        insert_SNP_IDs_addition(&(hap->dp.SNP_IDs), hap->dp.buffer, pathLen);
        return;
    }
    else
    {
        uint32_t* column;
        int j;

        column = Get_DP_Backtrack_Column(hap->dp, snpID);
        for (j = 0; j < (int)hap->dp.backtrack_length[snpID]; j++)
        {
            Preorder_Merge_Advance_Repeat(column[j], hap, pathLen); 
        }
    }
}

void generate_result_vector(haplotype_evdience_alloc* hap, int pathLen)
{
    if(pathLen != hap->dp.current_snp_num)
    {
        fprintf(stderr, "error\n");
    }


    int8_t* vector = Get_Result_SNP_Vector((*hap));
    memset(vector, -1, Get_SNP_Vector_Length((*hap)));
    hap->result_stat.occ_0 = 0;
    hap->result_stat.occ_1 = 0;
    hap->result_stat.occ_2 = 0;
    hap->result_stat.score = pathLen;
    hap->result_stat.homopolymer_num = 0;
    hap->result_stat.non_homopolymer_num = 0;

    long long snpID1;
    long long j = 0;
    int flag, vectorID;
    int current_score;
    for (j = 0; j < pathLen; j++)
    {
        snpID1 = hap->dp.buffer[j];
        vectorID = hap->snp_stat.a[snpID1].id;
        vector = Get_SNP_Vector((*hap), vectorID);
        if(hap->snp_stat.a[snpID1].is_homopolymer)
        {
            hap->result_stat.homopolymer_num++;
        }
        else
        {
            hap->result_stat.non_homopolymer_num++;
        }

        if((flag = debug_add_to_result_snp_vector(hap, vector, Get_SNP_Vector_Length((*hap))))!= -1)
        {
            fprintf(stderr, "incompatible snp vector....\n");
            exit(0);
        }
    }
    hap->result_stat.overlap_num = hap->result_stat.occ_0 + hap->result_stat.occ_1;
    

    

    ///check if this is a useful snp vector
    if(hap->result_stat.overlap_num !=0 && filter_one_snp_advance_nearby(hap, hap->result_stat.occ_0 + 1, 
    hap->result_stat.occ_1, hap->result_stat.overlap_num + 1,
    hap->result_stat.homopolymer_num, hap->result_stat.non_homopolymer_num,
    hap->dp.buffer, pathLen))
    {
        current_score = calculate_score(hap->result_stat.occ_0 + 1, hap->result_stat.occ_1);
        ///first useful snp vector
        if(hap->dp.max_snp_num < pathLen)
        {
            hap->dp.max_snp_num = pathLen;
            hap->dp.max_score = current_score;
            memcpy(hap->dp.max_buffer, hap->dp.buffer, sizeof(uint32_t) * pathLen);
        }///if we have multiple single best snp vector, select the vector with max score
        else if(hap->dp.max_snp_num == pathLen)
        {
            if(current_score > hap->dp.max_score)
            {
                hap->dp.max_score = current_score;
                memcpy(hap->dp.max_buffer, hap->dp.buffer, sizeof(uint32_t) * pathLen);
            }
        }
    }

}


void Preorder_Merge_Advance(uint32_t snpID, haplotype_evdience_alloc* hap, int pathLen) 
{ 
    hap->dp.visit[snpID] = 1;
    hap->dp.buffer[pathLen] = snpID;
    pathLen++;

    if(hap->dp.backtrack_length[snpID] == 0)
    {
        generate_result_vector(hap, pathLen);
        return;
    }
    else
    {
        uint32_t* column;
        int j;

        column = Get_DP_Backtrack_Column(hap->dp, snpID);
        for (j = 0; j < (int)hap->dp.backtrack_length[snpID]; j++)
        {
            Preorder_Merge_Advance(column[j], hap, pathLen); 
        }
    }
}



int if_snp_vector_useful(haplotype_evdience_alloc* hap,
long long occ_0, long long occ_1, uint32_t* SNPs, long long SNPsLen)
{

    double occ_1_coverage_low = (occ_0 + occ_1) * 0.3;
    
    if(occ_1 == 0 || occ_0 == 0)
    {
        return 0;
    }

    ///Fix-attention
    if(occ_1 >= occ_1_coverage_low && occ_0 >= occ_1_coverage_low)
    {
        return 1;
    }
    else if(occ_1 >= 5 && occ_0 >= 5)
    {
        return 1;
    }
    else if(occ_1 >= 3 && occ_0 >= 3 && SNPsLen >= 2)
    {
        int nearsnp;
        int non_nearsnps;
        count_nearby_snps(hap, SNPs, SNPsLen, &nearsnp, &non_nearsnps);
        if(non_nearsnps > 0)
        {
            return 1;
        }
    }

    return 0;
}


void merge_SNP_Vectors(haplotype_evdience_alloc* hap, uint32_t* SNPs, long long SNPLen)
{
    
    int8_t* vector = Get_Result_SNP_Vector((*hap));
    memset(vector, -1, Get_SNP_Vector_Length((*hap)));
    hap->result_stat.occ_0 = 0;
    hap->result_stat.occ_1 = 0;
    hap->result_stat.occ_2 = 0;
    hap->result_stat.score = SNPLen;
    hap->result_stat.homopolymer_num = 0;
    hap->result_stat.non_homopolymer_num = 0;

    long long snpID1;
    long long j = 0;
    int flag, vectorID;
    for (j = 0; j < SNPLen; j++)
    {
        snpID1 = SNPs[j];
        vectorID = hap->snp_stat.a[snpID1].id;
        vector = Get_SNP_Vector((*hap), vectorID);
        if(hap->snp_stat.a[snpID1].is_homopolymer)
        {
            hap->result_stat.homopolymer_num++;
        }
        else
        {
            hap->result_stat.non_homopolymer_num++;
        }

        if((flag = debug_add_to_result_snp_vector(hap, vector, Get_SNP_Vector_Length((*hap))))!= -1)
        {
            fprintf(stderr, "incompatible snp vector....\n");
            exit(0);
        }
    }
    hap->result_stat.overlap_num = hap->result_stat.occ_0 + hap->result_stat.occ_1;
}


void remove_reads(haplotype_evdience_alloc* hap, uint32_t* SNPs, long long SNPsLen, overlap_region_alloc* overlap_list)
{
    long long i, j, snpID, vectorID, overlapLen;
    int8_t *vector;

    for (i = 0; i < SNPsLen; i++)
    {
        snpID = SNPs[i];
        vectorID = hap->snp_stat.a[snpID].id;
        vector = Get_SNP_Vector((*hap), vectorID);
        ///hap->snp_stat[snpID].site;

        for (j = 0; j < Get_SNP_Vector_Length((*hap)); j++)
        {
            
            if(vector[j] == 1 && overlap_list->list[j].is_match == 1)
            {
                //overlap_list->list[j].is_match = 0;
                overlap_list->list[j].is_match = 2;
                overlapLen = overlap_list->list[j].x_pos_e - overlap_list->list[j].x_pos_s + 1;
                ///overlap_list->mapped_overlaps--;
                overlap_list->mapped_overlaps_length -= overlapLen;
            }

            /****************************may have bugs********************************/
            if( hap->snp_stat.a[snpID].site >= overlap_list->list[j].x_pos_s
                &&
                hap->snp_stat.a[snpID].site <= overlap_list->list[j].x_pos_e)
              {
                  overlap_list->list[j].strong = 1;
              }
            /****************************may have bugs********************************/

        }
    }
}

void try_to_remove_reads(int8_t* vector, long long vectorLen, overlap_region_alloc* overlap_list,
uint32_t* SNPs, long long SNPLen, haplotype_evdience_alloc* hap)
{
    long long i, overlapLen;
    long long removed_num = 0;

    for (i = 0; i < vectorLen; i++)
    {
        if(vector[i] == 1 && overlap_list->list[i].is_match == 1)
        {
            
            ///overlap_list->list[i].is_match = 0;
            overlap_list->list[i].is_match = 2;
            overlapLen = overlap_list->list[i].x_pos_e - overlap_list->list[i].x_pos_s + 1;
            ///overlap_list->mapped_overlaps--;
            overlap_list->mapped_overlaps_length -= overlapLen;
            removed_num++;
        }
    }


    long long snpID, j;
    for (i = 0; i < SNPLen; i++)
    {
        snpID = SNPs[i];
        ///check all overlaps
        for (j = 0; j < Get_SNP_Vector_Length((*hap)); j++)
        {
            /****************************may have bugs********************************/
            if( hap->snp_stat.a[snpID].site >= overlap_list->list[j].x_pos_s
                &&
                hap->snp_stat.a[snpID].site <= overlap_list->list[j].x_pos_e)
              {
                  overlap_list->list[j].strong = 1;
              }
            /****************************may have bugs********************************/
        }
    }
}


void process_repeat_snps(haplotype_evdience_alloc* hap, overlap_region_alloc* overlap_list)
{
    int i;
    uint32_t* snp_ids;
    long long length;

    

    for (i = 0; i < hap->dp.SNP_IDs.IDs_length; i++)
    {
        snp_ids = hap->dp.SNP_IDs.buffer + hap->dp.SNP_IDs.IDs[i].beg;
        length = hap->dp.SNP_IDs.IDs[i].end -hap->dp.SNP_IDs.IDs[i].beg + 1;

        merge_SNP_Vectors(hap, snp_ids, length);

        if(if_snp_vector_useful(hap, hap->result_stat.occ_0, hap->result_stat.occ_1, 
        snp_ids, length))
        {
            try_to_remove_reads(Get_Result_SNP_Vector((*hap)), Get_SNP_Vector_Length((*hap)), 
            overlap_list, snp_ids, length, hap);

            hap->dp.SNP_IDs.IDs[i].is_remove = 1;
        }
        else
        {
            hap->dp.SNP_IDs.IDs[i].is_remove = 0;
        }
    }
}


void lable_large_indels(overlap_region_alloc* overlap_list, long long read_length, Correct_dumy* dumy, double max_ov_diff_ec)
{
    long long i, j, c_i, c_n; uint32_t operLen; uint8_t oper;
    int is_delete = 0; window_list *c_idx;
    for (i = 0; i < (long long)overlap_list->length; i++)
    {
        ///should has at least 3 windows for this overlap
        if (overlap_list->list[i].is_match == 1 && overlap_list->list[i].w_list.n >= 3)
        {
            ///here w_list_length >= 3
            ///skip the first and last window
            for (j = 1; j + 1 < (long long)(overlap_list->list[i].w_list.n); j++)
            {
                ///this window is not matched, it seems to have large difference
                if(overlap_list->list[i].w_list.a[j].y_end == -1)
                {
                    overlap_list->list[i].is_match = 100;
                    is_delete = 1;
                    goto end_rem;
                }

                c_idx = &(overlap_list->list[i].w_list.a[j]);
                c_n = c_idx->clen;
                ///if there are <=2 cigar elements, skip it
                if(c_n < 3) continue;

                ///skip the first and last cigar elements
                for (c_i = 1; c_i + 1 < c_n; c_i++)
                {
                    get_cigar_cell(c_idx, &(overlap_list->list[i].w_list), c_i, &oper, &operLen);
                    
                    
                    if(operLen <= 5)
                    {
                        continue;
                    }
                    ///>=6 bp deletion or insertion 
                    if(oper == 2 || oper == 3)
                    {
                        overlap_list->list[i].is_match = 100;
                        is_delete = 1;
                        goto end_rem;
                    }
                }
            }
        }

        end_rem:;
    }


    if(is_delete == 1)
    {
        long long window_start, window_end;
        Window_Pool w_inf;
        init_Window_Pool(&w_inf, read_length, WINDOW, (int)(1.0/max_ov_diff_ec));
        int flag = 0;
        long long realLen = 0, realLen_100 = 0;
        int to_recover = 0;
        while(get_Window(&w_inf, &window_start, &window_end) && flag != -2)
        {
            dumy->length = 0;
            dumy->lengthNT = 0;
            ///return overlaps that is overlaped with [window_start, window_end]
            flag = get_available_fully_covered_interval(window_start, window_end, 
            overlap_list, dumy, &realLen, &realLen_100);


            switch (flag)
            {
                case 1:    ///match
                    break;
                case 0:    ///unmatch
                    break;
                case -2: ///unmatch, and the next window also cannot match
                    break;
            }

            ///it seems there is a long indel at the reference read itself
            if(realLen == 0 && realLen_100 > 0)
            {
                to_recover = 1;
                break;
            }
        }

        if(to_recover == 1)
        {
            for (i = 0; i < (long long)overlap_list->length; i++)
            {
                if (overlap_list->list[i].is_match == 100)
                {
                    overlap_list->list[i].is_match = 1;
                }
            }
        }
    }


    for (i = is_delete = 0; i < (long long)(overlap_list->length); i++)
    {
        if (overlap_list->list[i].is_match == 1)
        {
            overlap_list->list[i].without_large_indel = 1;
            is_delete++;
        }

        if (overlap_list->list[i].is_match == 100)
        {
            overlap_list->list[i].is_match = 1;
            overlap_list->list[i].without_large_indel = 0;
            is_delete++;
        }
    }

    // if(is_delete) radix_sort_overlap_region_dp_srt(overlap_list->list, overlap_list->list+overlap_list->length);
}

int debug_print_snp_stat(char* name, haplotype_evdience_alloc* hap, overlap_region_alloc* overlap_list, All_reads* R_INF)
{
    if(overlap_list->length > 0 &&
       memcmp(name, Get_NAME((*R_INF), overlap_list->list[0].x_id), 
       Get_NAME_LENGTH((*R_INF), overlap_list->list[0].x_id)) == 0)
    {
        fprintf(stderr, "\n%s, available_snp: %d\n", name, (int)hap->snp_stat.n);
        int i;
        for (i = 0; i < (int)hap->snp_stat.n; i++)
        {
            fprintf(stderr, "site: %d, occ_0: %d, occ_1: %d, occ_2: %d\n", 
            hap->snp_stat.a[i].site, hap->snp_stat.a[i].occ_0,
            hap->snp_stat.a[i].occ_1, hap->snp_stat.a[i].occ_2);
        }
    }

    return 1;  
}

int generate_haplotypes_DP(haplotype_evdience_alloc* hap, overlap_region_alloc* overlap_list, All_reads* R_INF, long long rLen, 
int force_repeat)
{
    int j, i;
    int vectorID, vectorID2;
    int8_t *vector, *vector2;

    
    // if(hap->available_snp == 0)
    if(hap->snp_stat.n == 0)
    {
        return 0;
    }

    ///debug_print_snp_stat("m64013_190324_024932/23660629/ccs", hap, overlap_list, R_INF);

    
    ///if hap->available_snp == 1, the following codes would have bugs
    ///filter snps that are highly likly false 
    // if(hap->available_snp > 1)
    if(hap->snp_stat.n > 1)
    {
        i = 0;
        ///if a snp is very close to others, it should not be a real snp
        for (j = 0; j < (int)hap->snp_stat.n/**hap->available_snp**/; j++)
        {
            if(j > 0 && j + 1 < (int)hap->snp_stat.n)
            {
                if(hap->snp_stat.a[j].site != hap->snp_stat.a[j - 1].site + 1
                    &&
                hap->snp_stat.a[j].site + 1 != hap->snp_stat.a[j + 1].site)
                {
                    hap->snp_stat.a[i] = hap->snp_stat.a[j];
                    i++;
                }

            }
            else if(j == 0)
            {
                if(hap->snp_stat.a[j].site + 1 != hap->snp_stat.a[j + 1].site)
                {
                    hap->snp_stat.a[i] = hap->snp_stat.a[j];
                    i++;
                }
            }
            else
            {
                if(hap->snp_stat.a[j].site != hap->snp_stat.a[j - 1].site + 1)
                {
                    hap->snp_stat.a[i] = hap->snp_stat.a[j];
                    i++;
                }
            }
        }
        // hap->available_snp = i;
        hap->snp_stat.n = i;
    }


    
    


    int flag;
    long long overlap_length, total_read, unuseful_read;
    total_read = unuseful_read = 0;
    ///check if any read may be conflict with others
    for (i = 0; i < (long long)overlap_list->length; i++)
    {
        overlap_length = overlap_list->list[i].x_pos_e - overlap_list->list[i].x_pos_s + 1;
        if (overlap_list->list[i].is_match == 1)
        {
            total_read++;
            flag = -1;
            for (j = 0; j < (int)hap->snp_stat.n; j++)
            {
                vectorID = hap->snp_stat.a[j].id;
                vector = Get_SNP_Vector((*hap), vectorID);

                ///flag == -1 means there are no useful signals yet
                if (flag == -1)
                {
                    if((vector[i] == 0 || vector[i] == 1 ))
                    {
                        flag = 0;
                    }
                }///flag == 0 means there is at least one useful signal yet
                else if (flag == 0) 
                {
                    if(vector[i] != 0 && vector[i] != 1)
                    {
                        flag = 2;
                    }
                }///flag == 0 means there is at least one useful signal first, and another unuseful signal after that
                else if(flag == 2) 
                {
                    if((vector[i] == 0 || vector[i] == 1 ))
                    {
                        flag = 3;
                        break;
                    }
                }
            }


            if(flag == 3) ///Fix-attention: definitely wrong
            {
                unuseful_read++;
                for (j = 0; j < (int)hap->snp_stat.n; j++)
                {
                    vectorID = hap->snp_stat.a[j].id;
                    vector = Get_SNP_Vector((*hap), vectorID);
                    

                    
                    if(vector[i] == 0)
                    {
                        hap->snp_stat.a[j].occ_0--;
                        hap->snp_stat.a[j].occ_2++;
                    }
                    else if(vector[i] == 1)
                    {
                        hap->snp_stat.a[j].occ_1--;
                        hap->snp_stat.a[j].occ_2++;
                    }
                    else if(vector[i] != 2)
                    {
                        hap->snp_stat.a[j].occ_2++;
                    }
                    

                    vector[i] = 2;
                }
                
                ///this read may be unuseful
                ///overlap_list->list[i].is_match = 0;
                ///overlap_list->list[i].is_match = 2;
                overlap_list->list[i].is_match = 4;
                ///overlap_list->mapped_overlaps--;
                overlap_list->mapped_overlaps_length -= overlap_length;
            }
        }
    }
    

    /*******************************DP********************************/
    init_DP_matrix(&(hap->dp), hap->snp_stat.n);

    long long equal_best = 0;
    uint32_t* column;

    

    for (i = 0; i < (int)hap->snp_stat.n; i++)
    {
        ///vector of snp i
        vectorID = hap->snp_stat.a[i].id;
        vector = Get_SNP_Vector((*hap), vectorID);
        hap->dp.visit[i] = 0;
        hap->dp.max[i] = 1;
        hap->dp.backtrack_length[i] = 0;
        equal_best = 0;
        column = Get_DP_Backtrack_Column(hap->dp, i);

        for (j = 0; j < i; j++)
        {
            ///vector of snp j
            vectorID2 = hap->snp_stat.a[j].id;
            vector2 = Get_SNP_Vector((*hap), vectorID2);

            ///vector is compatible with vector2
            if(calculate_distance_snp_vector(vector, vector2, Get_SNP_Vector_Length((*hap))) == 0)
            {
                
                if(hap->dp.max[i] < hap->dp.max[j] + 1)
                {
                    hap->dp.max[i] = hap->dp.max[j] + 1;

                    column[0] = j;
                    equal_best = 1;
                }
                else if(hap->dp.max[i] == hap->dp.max[j] + 1)
                {
                    column[equal_best] = j;
                    equal_best++;   
                }
                
                
            }
        }

        hap->dp.backtrack_length[i] = equal_best;
    }

    /*******************************DP********************************/

    uint64_t tmp_mode = 0;

    for (i = 0; i < (int)hap->snp_stat.n; i++)
    {
        tmp_mode = hap->dp.max[i];
        tmp_mode = tmp_mode << 32;
        tmp_mode = tmp_mode | (uint64_t)(i);
        hap->dp.max_for_sort[i] = tmp_mode;
    }
    
    qsort(hap->dp.max_for_sort, hap->snp_stat.n, sizeof(uint64_t), cmp_max_DP);


    int snpID;
    ///the minmum snp_num is 1
    hap->dp.max_snp_num = 0;
    hap->dp.max_score = -2;
    

    for (i = 0; i < (int)hap->snp_stat.n; i++)
    {
        snpID = Get_Max_DP_ID(hap->dp.max_for_sort[i]);
        if(hap->dp.visit[snpID] == 0)
        {
            hap->dp.current_snp_num = Get_Max_DP_Value(hap->dp.max_for_sort[i]);
            Preorder_Merge_Advance_Repeat(snpID, hap, 0);
        }
    }

    ///debug_print_snp_stat("m64013_190324_024932/23660629/ccs", hap, overlap_list, R_INF);


    //if(hap->dp.max_snp_num > 0)
    if(hap->snp_stat.n > 0)
    {
        process_repeat_snps(hap, overlap_list);
        return 1;
    }
    else
    {
        return 0;
    }
}


inline int check_informative_site(haplotype_evdience_alloc* hap, SnpStats* snp)
{
    long long vectorID = snp->id;
    int8_t *vector = Get_SNP_Vector((*hap), vectorID);
    snp->occ_0 = 0;
    snp->occ_1 = 0;
    snp->occ_2 = 0;
    long long i;
    for (i = 0; i < Get_SNP_Vector_Length((*hap)); i++)
    {
        if(vector[i] == 0)
        {
            snp->occ_0++;
        }
        else if(vector[i] == 1)
        {
            snp->occ_1++;
        }
        else if(vector[i] == 2)
        {
            snp->occ_2++;
        }
    }

    if(snp->occ_0 >= 2 || snp->occ_1 >= 2)
    {
        return 1;
    }

    return 0;
}


int generate_haplotypes_naive(haplotype_evdience_alloc* hap, overlap_region_alloc* overlap_list, All_reads* R_INF, long long rLen, 
int force_repeat)
{
    int j, i;

    
    if(hap->snp_stat.n == 0)
    {
        return 0;
    }

    ///if hap->available_snp == 1, the following codes would have bugs
    ///filter snps that are highly likly false 
    if(hap->snp_stat.n > 1)
    {
        i = 0;
        ///if a snp is very close to others, it should not be a real snp
        for (j = 0; j < (int)hap->snp_stat.n; j++)
        {
            
            if(j > 0 && j + 1 < (int)hap->snp_stat.n)
            {
                if(hap->snp_stat.a[j].site != hap->snp_stat.a[j - 1].site + 1
                    &&
                hap->snp_stat.a[j].site + 1 != hap->snp_stat.a[j + 1].site)
                {
                    hap->snp_stat.a[i] = hap->snp_stat.a[j];
                    i++;
                }

            }
            else if(j == 0)
            {
                if(hap->snp_stat.a[j].site + 1 != hap->snp_stat.a[j + 1].site)
                {
                    hap->snp_stat.a[i] = hap->snp_stat.a[j];
                    i++;
                }
            }
            else
            {
                if(hap->snp_stat.a[j].site != hap->snp_stat.a[j - 1].site + 1)
                {
                    hap->snp_stat.a[i] = hap->snp_stat.a[j];
                    i++;
                }
            }
        }
        hap->snp_stat.n = i;
    }


    long long m;
    if(hap->snp_stat.n > 0)
    {
        ///************************debug**************************///
        m = 0;
        for (i = 0; i < (int)hap->snp_stat.n; i++)
        {
            if(check_informative_site(hap, &(hap->snp_stat.a[i])))
            {
                hap->snp_stat.a[m] = hap->snp_stat.a[i];
                m++;
            }
        }
        hap->snp_stat.n = m;
        ///************************debug**************************///



        init_DP_matrix(&(hap->dp), hap->snp_stat.n);

        for (i = 0; i < (int)hap->snp_stat.n; i++)
        {
            hap->dp.max_buffer[i] = i;
        }
        hap->dp.max_snp_num = hap->snp_stat.n;
        remove_reads(hap, hap->dp.max_buffer, hap->dp.max_snp_num, overlap_list);

        return 1;

    }
    else
    {
        return 0;
    }

}


void generate_haplotypes_naive_advance(haplotype_evdience_alloc* hap, overlap_region_alloc* overlap_list, void *km)
{
    if(hap->length == 0) return;
    uint64_t k, l, i, o, *a, ii;
    int64_t z;
    SnpStats *s = NULL, *t = NULL;
    hap->snp_srt.n = 0;
    radix_sort_haplotype_evdience_id_srt(hap->list, hap->list + hap->length);
    for (k = 1, l = 0; k <= hap->length; ++k) {  
        if (k == hap->length || hap->list[k].overlapID != hap->list[l].overlapID) {
            for (i = l, o = 0; i < k; i++) {
                // if(hap->list[i].overlapSite == 55) {
                //     fprintf(stderr, "5555555555[M::%s::] utg%.6dl(%c), x_site::%u, x::[%u, %u), misBase::%c, occ_0::%u, occ_1::%u, occ_2::%u, snp_idx::%u\n", __func__, 
                //     ((int32_t)(overlap_list->list[hap->list[i].overlapID].y_id)) + 1,
                //     "+-"[overlap_list->list[hap->list[i].overlapID].y_pos_strand], 
                //     hap->list[i].site, overlap_list->list[hap->list[i].overlapID].x_pos_s,
                //     overlap_list->list[hap->list[i].overlapID].x_pos_e + 1, hap->list[i].misBase, 
                //     s->occ_0, s->occ_1, s->occ_2, hap->list[i].overlapSite);
                // }
                if(hap->list[i].type!=1) continue;///mismatch
                s = &(hap->snp_stat.a[hap->list[i].overlapSite]);
                // if(hap->list[i].overlapID == 125 || hap->list[i].overlapID == 127) {
                //     fprintf(stderr, "[M::%s::] utg%.6dl(%c), x_site::%u, x::[%u, %u), misBase::%c, occ_0::%u, occ_1::%u, occ_2::%u, snp_idx::%u\n", __func__, 
                //     ((int32_t)(overlap_list->list[hap->list[i].overlapID].y_id)) + 1,
                //     "+-"[overlap_list->list[hap->list[i].overlapID].y_pos_strand], 
                //     hap->list[i].site, overlap_list->list[hap->list[i].overlapID].x_pos_s,
                //     overlap_list->list[hap->list[i].overlapID].x_pos_e + 1, hap->list[i].misBase, 
                //     s->occ_0, s->occ_1, s->occ_2, hap->list[i].overlapSite);
                // }
                if(s->occ_0 < 2 || s->occ_1 < 2) continue;
                if(s->occ_0 >= asm_opt.s_hap_cov && s->occ_1 >= asm_opt.infor_cov) o++;///allels must be real
            }

            // if(hap->list[l].overlapID == 125 || hap->list[l].overlapID == 127) {
            //     fprintf(stderr, "+++[M::%s::] utg%.6dl(%c), o::%lu\n", __func__, 
            //     ((int32_t)(overlap_list->list[hap->list[l].overlapID].y_id)) + 1,
            //     "+-"[overlap_list->list[hap->list[l].overlapID].y_pos_strand], o);
            // }
            if(o > 0) {
                o = ((uint32_t)-1) - o;
                o <<= 32; o += l; 
                if(!km) kv_push(uint64_t, hap->snp_srt, o);
                else kv_push_km(km, uint64_t, hap->snp_srt, o);
            }
            l = k;
        }
    }
    // fprintf(stderr, "\nhap->snp_srt.n: %u, overlap_list->length: %lu, x_id: %u\n", 
    // (uint32_t)hap->snp_srt.n, overlap_list->length, overlap_list->list[0].x_id);
    if (hap->snp_srt.n > 0) {
        radix_sort_bc64(hap->snp_srt.a, hap->snp_srt.a + hap->snp_srt.n);///sort by how many snps in one overlap
        for (k = 0; k < hap->snp_srt.n; k++) {
            o = 0; l = (uint32_t)hap->snp_srt.a[k];
            for (i = l; i < hap->length && hap->list[i].overlapID == hap->list[l].overlapID; i++) {
                if(hap->list[i].type!=1) continue;
                s = &(hap->snp_stat.a[hap->list[i].overlapSite]);
                if(s->occ_0 < 2 || s->occ_1 < 2) continue;
                if(s->occ_0 >= asm_opt.s_hap_cov && s->occ_1 >= asm_opt.infor_cov) o++;
            }
            // if(hap->list[l].overlapID == 125 || hap->list[l].overlapID == 127) {
            //     fprintf(stderr, "sbsbsb[M::%s::] utg%.6dl(%c), o::%lu\n", __func__, 
            //     ((int32_t)(overlap_list->list[hap->list[l].overlapID].y_id)) + 1,
            //     "+-"[overlap_list->list[hap->list[l].overlapID].y_pos_strand], o);
            // }
            if(o == 0) continue;

            ii = hap->list[l].overlapID;
            if(overlap_list->list[ii].is_match == 1) overlap_list->list[ii].is_match = 2;
            for (i = l; i < hap->length && hap->list[i].overlapID == hap->list[l].overlapID; i++) {
                if(hap->list[i].type==1){
                    s = &(hap->snp_stat.a[hap->list[i].overlapSite]);
                    s->score = 1; 
                } else if(hap->list[i].type==0) {
                    z = hap->list[i].overlapSite; s = &(hap->snp_stat.a[z]);
                    for (z = hap->list[i].overlapSite; z >= 0; z--) {
                        t = &(hap->snp_stat.a[z]);
                        if(s->site!=t->site) break;
                        t->occ_0 -= hap->list[i].cov;
                        assert(t->occ_0 >= 1);// if(t->occ_0 < 1) fprintf(stderr, "WRONG-CORRECTION\n");
                    }
                }
            }
        }

        for (k = 0; k < hap->snp_srt.n; k++) {
            o = 0; l = (uint32_t)hap->snp_srt.a[k];
            for (i = l; i < hap->length && hap->list[i].overlapID == hap->list[l].overlapID; i++) {
                if(hap->list[i].type!=1) continue;
                s = &(hap->snp_stat.a[hap->list[i].overlapSite]);
                if(s->occ_0 < 2 || s->occ_1 < 2) continue;
                if(s->score == 1) {
                    o++;
                    // if(hap->list[i].overlapID == 125 || hap->list[i].overlapID == 127) {
                    //     fprintf(stderr, "[M::%s::] utg%.6dl(%c), x_site::%u, x::[%u, %u), misBase::%c, occ_0::%u, occ_1::%u, occ_2::%u, snp_idx::%u\n", __func__, 
                    //     ((int32_t)(overlap_list->list[hap->list[i].overlapID].y_id)) + 1,
                    //     "+-"[overlap_list->list[hap->list[i].overlapID].y_pos_strand], 
                    //     hap->list[i].site, overlap_list->list[hap->list[i].overlapID].x_pos_s,
                    //     overlap_list->list[hap->list[i].overlapID].x_pos_e + 1, hap->list[i].misBase, 
                    //     s->occ_0, s->occ_1, s->occ_2, hap->list[i].overlapSite);
                    // }
                }
            }
            ii = hap->list[l].overlapID;
            // if(hap->list[l].overlapID == 125 || hap->list[l].overlapID == 127) {
            //     fprintf(stderr, ">>>[M::%s::] utg%.6dl(%c), o::%lu\n", __func__, 
            //     ((int32_t)(overlap_list->list[hap->list[l].overlapID].y_id)) + 1,
            //     "+-"[overlap_list->list[hap->list[l].overlapID].y_pos_strand], o);
            // }
            if(overlap_list->list[ii].is_match == 2 && o == 0) {
                overlap_list->list[ii].is_match = 1;
            }
            if(overlap_list->list[ii].is_match == 1 && o > 0) {
                overlap_list->list[ii].is_match = 2;
            }
        }
    }

    hap->snp_srt.n = 0;
    for (k = 1, l = 0; k <= hap->length; ++k) {  
        if (k == hap->length || hap->list[k].overlapID != hap->list[l].overlapID) {
            if(overlap_list->list[hap->list[l].overlapID].is_match == 2) {
                l = k;
                continue;
            }
            for (i = l, o = 0; i < k; i++) {
                if(hap->list[i].type!=1) continue;
                s = &(hap->snp_stat.a[hap->list[i].overlapSite]);
                if(s->occ_0 < 2 || s->occ_1 < 2) continue;
                if(s->score == 1) continue;
                 o++; 
                 if(!km) kv_push(uint64_t, hap->snp_srt, hap->list[i].overlapSite);
                 else kv_push_km(km, uint64_t, hap->snp_srt, hap->list[i].overlapSite);
            }
            // if(hap->list[l].overlapID == 125 || hap->list[l].overlapID == 127) {
            //     fprintf(stderr, "---[M::%s::] utg%.6dl(%c), o::%lu\n", __func__, 
            //     ((int32_t)(overlap_list->list[hap->list[l].overlapID].y_id)) + 1,
            //     "+-"[overlap_list->list[hap->list[l].overlapID].y_pos_strand], o);
            // }
            hap->snp_srt.n -= o;
            if(o >= 2) {///there are at least two variants at one read
                radix_sort_bc64(hap->snp_srt.a + hap->snp_srt.n, hap->snp_srt.a + hap->snp_srt.n + o);
                a = hap->snp_srt.a + hap->snp_srt.n;
                for (i = z = 0; i < o; i++) {
                    if(i > 0) s = &(hap->snp_stat.a[a[i-1]]);
                    if(i + 1 < o) t = &(hap->snp_stat.a[a[i+1]]);
                    if(s && s->site + 32 > hap->snp_stat.a[a[i]].site) continue;
                    if(t && hap->snp_stat.a[a[i]].site + 32 > t->site) continue;
                    a[z] = a[i];
                    z++;
                }
                if(z >= 2) hap->snp_srt.n += z; 
            }
            l = k;
        }
    }
    if (hap->snp_srt.n > 0) {
        radix_sort_bc64(hap->snp_srt.a, hap->snp_srt.a + hap->snp_srt.n);
        for (k = 1, l = 0; k <= hap->snp_srt.n; ++k) {  
            if(k == hap->snp_srt.n || hap->snp_srt.a[k] != hap->snp_srt.a[l]) {
                if(k - l >= 2) hap->snp_stat.a[hap->snp_srt.a[l]].score = 1;
            }
            l = k;
        }
    }
    
    for (k = 1, l = 0; k <= hap->length; ++k) {  
        if (k == hap->length || hap->list[k].overlapID != hap->list[l].overlapID) {
            ii = hap->list[l].overlapID;
            if(overlap_list->list[ii].is_match==2) overlap_list->list[ii].is_match = 1;
            if(overlap_list->list[ii].is_match==1) {
                for (i = l; i < k; i++) {
                    if(hap->list[i].type==1 || hap->list[i].type==0) {
                        s = &(hap->snp_stat.a[hap->list[i].overlapSite]);
                        if(s->score == 1 && (!(s->occ_0 < 2 || s->occ_1 < 2))) {
                            overlap_list->list[ii].strong = 1;
                            if(hap->list[i].type==1) {
                                overlap_list->list[ii].is_match = 2;
                                overlap_list->mapped_overlaps_length -= 
                                            overlap_list->list[ii].x_pos_e + 1 - overlap_list->list[ii].x_pos_s;
                                break;
                            }
                        }
                    }

                }
            }
            l = k;
        }
    }

    // for (i = k = 0; i < overlap_list->length; i++) {
    //     if(overlap_list->list[i].is_match == 2) k++;
    // }
    // for (i = l = 0; i < hap->snp_stat.n; i++) {
    //     s = &(hap->snp_stat.a[i]);
    //     if(s->score == 1 && (!(s->occ_0 < 2 || s->occ_1 < 2))) l++;
    // }
    // fprintf(stderr, "#trans ovlp: %lu, # snp:: %lu\n", k, l);
    
}


void generate_haplotypes_naive_HiFi(haplotype_evdience_alloc* hap, overlap_region_alloc* overlap_list, double up, void *km)
{
    // fprintf(stderr, "[M::%s::] Done\n", __func__);
    if(hap->length == 0) return;
    uint64_t k, l, i, o, *a, ii, m_snp_stat, m_list, m_off;
    int64_t z; 
    SnpStats *s = NULL, *t = NULL;

    for (k = 1, l = 0, i = m_snp_stat = m_list = 0; k <= hap->snp_stat.n; ++k) {///filter snps
        if(k == hap->snp_stat.n || hap->snp_stat.a[k].site != hap->snp_stat.a[l].site) {
            if((l > 0) && (hap->snp_stat.a[l].site == (hap->snp_stat.a[l-1].site + 1))) {
                l = k; continue;
            }
            if((k < hap->snp_stat.n) && ((hap->snp_stat.a[l].site+1) == hap->snp_stat.a[k].site)) {
                l = k; continue;
            }

            for (; i < hap->length && hap->list[i].site != hap->snp_stat.a[l].site; i++);
            assert(i < hap->length && hap->list[i].site == hap->snp_stat.a[l].site);
            m_off = l - m_snp_stat; 
            for (; i < hap->length && hap->list[i].site == hap->snp_stat.a[l].site; i++) {
                assert(hap->list[i].overlapSite>=l && hap->list[i].overlapSite < k);
                // assert(hap->snp_stat.a[hap->list[i].overlapSite].site==hap->list[i].site);
                hap->list[m_list] = hap->list[i]; hap->list[m_list++].overlapSite -= m_off;
            }
        
            for (; l < k; l++) hap->snp_stat.a[m_snp_stat++] = hap->snp_stat.a[l];
        }
    }
    hap->snp_stat.n = m_snp_stat; hap->length = m_list;
    if(hap->snp_stat.n == 0 || hap->length == 0) return;

    hap->snp_srt.n = 0;
    radix_sort_haplotype_evdience_id_srt(hap->list, hap->list + hap->length);
    for (k = 1, l = 0; k <= hap->length; ++k) {  
        if (k == hap->length || hap->list[k].overlapID != hap->list[l].overlapID) {
            for (i = l, o = 0; i < k; i++) {
                if(hap->list[i].type!=1) continue;///mismatch
                s = &(hap->snp_stat.a[hap->list[i].overlapSite]); 
                assert(s->site == hap->list[i].site);
                if(s->occ_0 < 2 || s->occ_1 < 2) continue;
                if(s->occ_0 >= asm_opt.s_hap_cov && s->occ_1 >= asm_opt.infor_cov) o++;///allels must be real
            }
            if(o > 0) {
                o = ((uint32_t)-1) - o;
                o <<= 32; o += l; 
                if(!km) kv_push(uint64_t, hap->snp_srt, o);
                else kv_push_km(km, uint64_t, hap->snp_srt, o);
            }
            l = k;
        }
    }
    
    if (hap->snp_srt.n > 0) {
        radix_sort_bc64(hap->snp_srt.a, hap->snp_srt.a + hap->snp_srt.n);///sort by how many snps in one overlap
        for (k = 0; k < hap->snp_srt.n; k++) {
            o = 0; l = (uint32_t)hap->snp_srt.a[k];
            for (i = l; i < hap->length && hap->list[i].overlapID == hap->list[l].overlapID; i++) {
                if(hap->list[i].type!=1) continue;
                s = &(hap->snp_stat.a[hap->list[i].overlapSite]);
                if(s->occ_0 < 2 || s->occ_1 < 2) continue;
                if(s->occ_0 >= asm_opt.s_hap_cov && s->occ_1 >= asm_opt.infor_cov) o++;
            }
            if(o == 0) continue;

            ii = hap->list[l].overlapID;
            if(overlap_list->list[ii].is_match == 1) overlap_list->list[ii].is_match = 2;
            for (i = l; i < hap->length && hap->list[i].overlapID == hap->list[l].overlapID; i++) {
                if(hap->list[i].type==1){
                    s = &(hap->snp_stat.a[hap->list[i].overlapSite]);
                    s->score = 1; 
                } ///else if((hap->list[i].type==0) && (o>=(overlap_list->list[ii].align_length*up))) {
                else if(hap->list[i].type==0) {
                    ///not real allels
                    z = hap->list[i].overlapSite; s = &(hap->snp_stat.a[z]);
                    for (z = hap->list[i].overlapSite; z >= 0; z--) {
                        t = &(hap->snp_stat.a[z]);
                        if(s->site!=t->site) break;
                        t->occ_0 -= hap->list[i].cov;
                        assert(t->occ_0 >= 1);
                    }
                }
            }
        }

        for (k = 0; k < hap->snp_srt.n; k++) {///sorted by how many allels in each overlap; more -> less
            o = 0; l = (uint32_t)hap->snp_srt.a[k];
            for (i = l; i < hap->length && hap->list[i].overlapID == hap->list[l].overlapID; i++) {
                if(hap->list[i].type!=1) continue;
                s = &(hap->snp_stat.a[hap->list[i].overlapSite]);
                if(s->occ_0 < 2 || s->occ_1 < 2) continue;
                if(s->score == 1) o++;
            }
            ii = hap->list[l].overlapID;
            ///for HiFi, do not flip trans to cis
            // if(overlap_list->list[ii].is_match == 2 && o == 0) {
            //     overlap_list->list[ii].is_match = 1;
            // }
            if(overlap_list->list[ii].is_match == 1 && o > 0) {
                overlap_list->list[ii].is_match = 2;
            }
        }


        for (k = 1, l = 0; k <= hap->length; ++k) {  ///reset snp_stat
            if (k == hap->length || hap->list[k].overlapID != hap->list[l].overlapID) {
                ii = hap->list[l].overlapID;
                if(overlap_list->list[ii].is_match==1) {
                    for (i = l; i < k; i++) {
                        if(hap->list[i].type==1) {
                            hap->snp_stat.a[hap->list[i].overlapSite].score = -1;
                        }
                    }
                }
                l = k;
            }
        } 
    }

    

    hap->snp_srt.n = 0;
    for (k = 1, l = 0; k <= hap->length; ++k) {  
        if (k == hap->length || hap->list[k].overlapID != hap->list[l].overlapID) {
            if(overlap_list->list[hap->list[l].overlapID].is_match == 2) {
                l = k;
                continue;
            }
            for (i = l, o = 0; i < k; i++) {
                if(hap->list[i].type!=1) continue;
                s = &(hap->snp_stat.a[hap->list[i].overlapSite]);
                if(s->occ_0 < 2 || s->occ_1 < 2) continue;
                if(s->occ_0 >= asm_opt.s_hap_cov && s->occ_1 >= asm_opt.infor_cov) continue;
                if(s->score == 1) continue;
                o++; 
                if(!km) kv_push(uint64_t, hap->snp_srt, hap->list[i].overlapSite);
                else kv_push_km(km, uint64_t, hap->snp_srt, hap->list[i].overlapSite);
            }
            hap->snp_srt.n -= o;
            ///there are at least two variants at one read
            if(o>=(overlap_list->list[hap->list[l].overlapID].align_length*up)) {
                radix_sort_bc64(hap->snp_srt.a + hap->snp_srt.n, hap->snp_srt.a + hap->snp_srt.n + o);
                a = hap->snp_srt.a + hap->snp_srt.n;
                for (i = z = 0; i < o; i++) {
                    if(i > 0) s = &(hap->snp_stat.a[a[i-1]]);
                    if(i + 1 < o) t = &(hap->snp_stat.a[a[i+1]]);
                    if(s && s->site + 32 > hap->snp_stat.a[a[i]].site) continue;
                    if(t && hap->snp_stat.a[a[i]].site + 32 > t->site) continue;
                    a[z] = a[i];
                    z++;
                }
                if(z >= 2) hap->snp_srt.n += z; 
            }
            l = k;
        }
    }
    if (hap->snp_srt.n > 0) {
        radix_sort_bc64(hap->snp_srt.a, hap->snp_srt.a + hap->snp_srt.n);
        for (k = 1, l = 0; k <= hap->snp_srt.n; ++k) {  
            if(k == hap->snp_srt.n || hap->snp_srt.a[k] != hap->snp_srt.a[l]) {
                if(k - l >= 2) hap->snp_stat.a[hap->snp_srt.a[l]].score = 1;
            }
            l = k;
        }
    }

    for (k = 1, l = 0; k <= hap->length; ++k) {  
        if (k == hap->length || hap->list[k].overlapID != hap->list[l].overlapID) {
            ii = hap->list[l].overlapID;
            if(overlap_list->list[ii].is_match==2) {
                overlap_list->list[ii].strong = 1;
                overlap_list->mapped_overlaps_length -= 
                                            overlap_list->list[ii].x_pos_e + 1 - overlap_list->list[ii].x_pos_s;
            } else if(overlap_list->list[ii].is_match==1) {
                for (i = l; i < k; i++) {
                    if(hap->list[i].type==1 || hap->list[i].type==0) {
                        s = &(hap->snp_stat.a[hap->list[i].overlapSite]);
                        if(s->score == 1 && (!(s->occ_0 < 2 || s->occ_1 < 2))) {
                            overlap_list->list[ii].strong = 1;
                            if(hap->list[i].type==1) {
                                overlap_list->list[ii].is_match = 2;
                                overlap_list->mapped_overlaps_length -= 
                                            overlap_list->list[ii].x_pos_e + 1 - overlap_list->list[ii].x_pos_s;
                                break;
                            }
                        }
                    }

                }
            }
            l = k;
        }
    }    
}


void generate_haplotypes_naive_UL(haplotype_evdience_alloc* hap, overlap_region_alloc* overlap_list, double up, void *km)
{
    if(hap->length == 0) return;
    uint64_t k, l, i, o, *a, ii, m_snp_stat, m_list, m_off;
    int64_t z; 
    SnpStats *s = NULL, *t = NULL;

    for (k = 1, l = 0, i = m_snp_stat = m_list = 0; k <= hap->snp_stat.n; ++k) {///filter snps
        if(k == hap->snp_stat.n || hap->snp_stat.a[k].site != hap->snp_stat.a[l].site) {
            if((l > 0) && (hap->snp_stat.a[l].site == (hap->snp_stat.a[l-1].site + 1))) {
                l = k; continue;
            }
            if((k < hap->snp_stat.n) && ((hap->snp_stat.a[l].site+1) == hap->snp_stat.a[k].site)) {
                l = k; continue;
            }

            for (; i < hap->length && hap->list[i].site != hap->snp_stat.a[l].site; i++);
            assert(i < hap->length && hap->list[i].site == hap->snp_stat.a[l].site);
            m_off = l - m_snp_stat; 
            for (; i < hap->length && hap->list[i].site == hap->snp_stat.a[l].site; i++) {
                assert(hap->list[i].overlapSite>=l && hap->list[i].overlapSite < k);
                // assert(hap->snp_stat.a[hap->list[i].overlapSite].site==hap->list[i].site);
                hap->list[m_list] = hap->list[i]; hap->list[m_list++].overlapSite -= m_off;
            }
        
            for (; l < k; l++) hap->snp_stat.a[m_snp_stat++] = hap->snp_stat.a[l];
        }
    }
    hap->snp_stat.n = m_snp_stat; hap->length = m_list;
    if(hap->snp_stat.n == 0 || hap->length == 0) return;

    hap->snp_srt.n = 0;
    radix_sort_haplotype_evdience_id_srt(hap->list, hap->list + hap->length);
    for (k = 1, l = 0; k <= hap->length; ++k) {  
        if (k == hap->length || hap->list[k].overlapID != hap->list[l].overlapID) {
            for (i = l, o = 0; i < k; i++) {
                if(hap->list[i].type!=1) continue;///mismatch
                s = &(hap->snp_stat.a[hap->list[i].overlapSite]); 
                assert(s->site == hap->list[i].site);
                if(s->occ_0 < 2 || s->occ_1 < 2) continue;
                if(s->occ_0 >= asm_opt.s_hap_cov && s->occ_1 >= asm_opt.infor_cov) o++;///allels must be real
            }
            if(o > 0) {
                o = ((uint32_t)-1) - o;
                o <<= 32; o += l; 
                if(!km) kv_push(uint64_t, hap->snp_srt, o);
                else kv_push_km(km, uint64_t, hap->snp_srt, o);
            }
            l = k;
        }
    }
    
    if (hap->snp_srt.n > 0) {
        radix_sort_bc64(hap->snp_srt.a, hap->snp_srt.a + hap->snp_srt.n);///sort by how many snps in one overlap
        for (k = 0; k < hap->snp_srt.n; k++) {
            o = 0; l = (uint32_t)hap->snp_srt.a[k];
            for (i = l; i < hap->length && hap->list[i].overlapID == hap->list[l].overlapID; i++) {
                if(hap->list[i].type!=1) continue;
                s = &(hap->snp_stat.a[hap->list[i].overlapSite]);
                if(s->occ_0 < 2 || s->occ_1 < 2) continue;
                if(s->occ_0 >= asm_opt.s_hap_cov && s->occ_1 >= asm_opt.infor_cov) o++;
            }
            if(o == 0) continue;

            ii = hap->list[l].overlapID;
            if(overlap_list->list[ii].is_match == 1) overlap_list->list[ii].is_match = 2;
            for (i = l; i < hap->length && hap->list[i].overlapID == hap->list[l].overlapID; i++) {
                if(hap->list[i].type==1){
                    s = &(hap->snp_stat.a[hap->list[i].overlapSite]);
                    s->score = 1; 
                } ///else if((hap->list[i].type==0) && (o>=(overlap_list->list[ii].align_length*up))) {
                else if(hap->list[i].type==0) {
                    ///not real allels
                    z = hap->list[i].overlapSite; s = &(hap->snp_stat.a[z]);
                    for (z = hap->list[i].overlapSite; z >= 0; z--) {
                        t = &(hap->snp_stat.a[z]);
                        if(s->site!=t->site) break;
                        t->occ_0 -= hap->list[i].cov;
                        assert(t->occ_0 >= 1);
                    }
                }
            }
        }

        for (k = 0; k < hap->snp_srt.n; k++) {///sorted by how many allels in each overlap; more -> less
            o = 0; l = (uint32_t)hap->snp_srt.a[k];
            for (i = l; i < hap->length && hap->list[i].overlapID == hap->list[l].overlapID; i++) {
                if(hap->list[i].type!=1) continue;
                s = &(hap->snp_stat.a[hap->list[i].overlapSite]);
                if(s->occ_0 < 2 || s->occ_1 < 2) continue;
                if(s->score == 1) o++;
            }
            ii = hap->list[l].overlapID;
            ///for HiFi, do not flip trans to cis
            // if(overlap_list->list[ii].is_match == 2 && o == 0) {
            //     overlap_list->list[ii].is_match = 1;
            // }
            if(overlap_list->list[ii].is_match == 1 && o > 0) {
                overlap_list->list[ii].is_match = 2;
            }
        }


        for (k = 1, l = 0; k <= hap->length; ++k) {  ///reset snp_stat
            if (k == hap->length || hap->list[k].overlapID != hap->list[l].overlapID) {
                ii = hap->list[l].overlapID;
                if(overlap_list->list[ii].is_match==1) {
                    for (i = l; i < k; i++) {
                        if(hap->list[i].type==1) {
                            hap->snp_stat.a[hap->list[i].overlapSite].score = -1;
                        }
                    }
                }
                l = k;
            }
        } 
    }

    

    hap->snp_srt.n = 0;
    for (k = 1, l = 0; k <= hap->length; ++k) {  
        if (k == hap->length || hap->list[k].overlapID != hap->list[l].overlapID) {
            if(overlap_list->list[hap->list[l].overlapID].is_match == 2) {
                l = k;
                continue;
            }
            for (i = l, o = 0; i < k; i++) {
                if(hap->list[i].type!=1) continue;
                s = &(hap->snp_stat.a[hap->list[i].overlapSite]);
                if(s->occ_0 < 2 || s->occ_1 < 2) continue;
                if(s->occ_0 >= asm_opt.s_hap_cov && s->occ_1 >= asm_opt.infor_cov) continue;
                if(s->score == 1) continue;
                o++; 
                if(!km) kv_push(uint64_t, hap->snp_srt, hap->list[i].overlapSite);
                else kv_push_km(km, uint64_t, hap->snp_srt, hap->list[i].overlapSite);
            }
            hap->snp_srt.n -= o;
            ///there are at least two variants at one read
            if(o>=(overlap_list->list[hap->list[l].overlapID].align_length*up)) {
                radix_sort_bc64(hap->snp_srt.a + hap->snp_srt.n, hap->snp_srt.a + hap->snp_srt.n + o);
                a = hap->snp_srt.a + hap->snp_srt.n;
                for (i = z = 0; i < o; i++) {
                    if(i > 0) s = &(hap->snp_stat.a[a[i-1]]);
                    if(i + 1 < o) t = &(hap->snp_stat.a[a[i+1]]);
                    if(s && s->site + 32 > hap->snp_stat.a[a[i]].site) continue;
                    if(t && hap->snp_stat.a[a[i]].site + 32 > t->site) continue;
                    a[z] = a[i];
                    z++;
                }
                if(z >= 2) hap->snp_srt.n += z; 
            }
            l = k;
        }
    }
    if (hap->snp_srt.n > 0) {
        radix_sort_bc64(hap->snp_srt.a, hap->snp_srt.a + hap->snp_srt.n);
        for (k = 1, l = 0; k <= hap->snp_srt.n; ++k) {  
            if(k == hap->snp_srt.n || hap->snp_srt.a[k] != hap->snp_srt.a[l]) {
                if(k - l >= 2) hap->snp_stat.a[hap->snp_srt.a[l]].score = 1;
            }
            l = k;
        }
    }

    for (k = 1, l = 0; k <= hap->length; ++k) {  
        if (k == hap->length || hap->list[k].overlapID != hap->list[l].overlapID) {
            ii = hap->list[l].overlapID;
            if(overlap_list->list[ii].is_match==2) {
                overlap_list->list[ii].strong = 1;
                overlap_list->mapped_overlaps_length -= 
                                            overlap_list->list[ii].x_pos_e + 1 - overlap_list->list[ii].x_pos_s;
            } else if(overlap_list->list[ii].is_match==1) {
                for (i = l; i < k; i++) {
                    if(hap->list[i].type==1 || hap->list[i].type==0) {
                        s = &(hap->snp_stat.a[hap->list[i].overlapSite]);
                        if(s->score == 1 && (!(s->occ_0 < 2 || s->occ_1 < 2))) {
                            overlap_list->list[ii].strong = 1;
                            if(hap->list[i].type==1) {
                                overlap_list->list[ii].is_match = 2;
                                overlap_list->mapped_overlaps_length -= 
                                            overlap_list->list[ii].x_pos_e + 1 - overlap_list->list[ii].x_pos_s;
                                break;
                            }
                        }
                    }

                }
            }
            l = k;
        }
    }    
}

/**
void partition_overlaps(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, haplotype_evdience_alloc* hap,
                        int force_repeat)
{
    ResizeInitHaplotypeEvdience(hap);

    long long i;
    long long window_start, window_end;

    long long num_availiable_win = 0;
    
    Window_Pool w_inf;
    init_Window_Pool(&w_inf, g_read->length, WINDOW, (int)(1.0/asm_opt.max_ov_diff_ec));

    int flag = 0;
    while(get_Window(&w_inf, &window_start, &window_end) && flag != -2)
    {
        dumy->length = 0;
        dumy->lengthNT = 0;

        ///return overlaps that is overlaped with [window_start, window_end]
        flag = get_available_interval(window_start, window_end, overlap_list, dumy);
        switch (flag)
        {
            case 1:    ///found matched overlaps
                break;
            case 0:    ///do not find any matched overlaps
                break;
            case -2: ///do not find any matched overlaps, and the next window also cannot match
                break;
        }
        
        num_availiable_win = num_availiable_win + dumy->length;

        cluster(g_read->seq, window_start, window_end, overlap_list, dumy, R_INF, hap);
    }

    ///very time-consuming
    qsort(hap->list, hap->length, sizeof(haplotype_evdience), cmp_haplotype_evdience);

    
    

    ///debug_hap_information(overlap_list, R_INF, g_read, hap, dumy);

    SetSnpMatrix(hap, &(hap->nn_snp), &(overlap_list->length), 1, NULL);


    uint64_t pre_site = (uint64_t)-1;
    uint64_t num_of_snps = 0;
    long long pre_i = -1;
    long long sub_length;
    haplotype_evdience* sub_list;

    ////split reads
    for (i = 0; i < hap->length; i++)
    {
        if(pre_site != hap->list[i].site)
        {
            if(i != 0)
            {
                sub_list = hap->list + pre_i;
                sub_length = i - pre_i;
                split_sub_list(hap, sub_list, sub_length, overlap_list, R_INF, g_read);
            }
            num_of_snps++;
            pre_site = hap->list[i].site;
            pre_i = i;
        }
    }

    if(pre_i != -1)
    {
        sub_list = hap->list + pre_i;
        sub_length = i - pre_i;
        split_sub_list(hap, sub_list, sub_length, overlap_list, R_INF, g_read);
    }

    ///debug_snp_matrix(hap);
    generate_haplotypes_DP(hap, overlap_list, R_INF, g_read->length, force_repeat);
    ///generate_haplotypes_naive(hap, overlap_list, R_INF, g_read->length, force_repeat);

    lable_large_indels(overlap_list, g_read->length, dumy, asm_opt.max_ov_diff_ec);


    ///debug_snp_matrix(hap);
}
**/

void partition_overlaps_advance_back(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, UC_Read* overlap_read, Correct_dumy* dumy, 
                        haplotype_evdience_alloc* hap, int force_repeat)
{
    ResizeInitHaplotypeEvdience(hap);

    long long i;
    long long window_start, window_end;

    long long num_availiable_win = 0;
    
    Window_Pool w_inf;
    init_Window_Pool(&w_inf, g_read->length, WINDOW, (int)(1.0/asm_opt.max_ov_diff_ec));

    int flag = 0;
    while(get_Window(&w_inf, &window_start, &window_end) && flag != -2)
    {
        dumy->length = 0;
        dumy->lengthNT = 0;
        ///return overlaps that is overlaped with [window_start, window_end]
        flag = get_available_interval(window_start, window_end, overlap_list, dumy);
        switch (flag)
        {
            case 1:    ///found matched overlaps
                break;
            case 0:    ///do not find any matched overlaps
                break;
            case -2: ///do not find any matched overlaps, and the next window also cannot match
                break;
        }
        
        num_availiable_win = num_availiable_win + dumy->length;
        
        ///need to deal with 
        cluster_advance(g_read->seq, window_start, window_end, overlap_list, 
        dumy, R_INF, hap, overlap_read, 1);
    }


    
    ///very time-consuming
    ///Fix-attention ---> able to be sorted locally
    // qsort(hap->list, hap->length, sizeof(haplotype_evdience), cmp_haplotype_evdience);
    SetSnpMatrix(hap, &(hap->nn_snp), &(overlap_list->length), 1, NULL);

    uint64_t pre_site = (uint64_t)-1;
    uint64_t num_of_snps = 0;
    long long pre_i = -1;
    long long sub_length;
    haplotype_evdience* sub_list;

    ////split reads
    for (i = 0; i < hap->length; i++)
    {
        if(pre_site != hap->list[i].site)
        {
            if(i != 0)
            {
                sub_list = hap->list + pre_i;
                sub_length = i - pre_i;
                split_sub_list(hap, sub_list, sub_length, overlap_list, R_INF, g_read);
            }
            num_of_snps++;
            pre_site = hap->list[i].site;
            pre_i = i;
        }
    }

    if(pre_i != -1)
    {
        sub_list = hap->list + pre_i;
        sub_length = i - pre_i;
        split_sub_list(hap, sub_list, sub_length, overlap_list, R_INF, g_read);
    }

    generate_haplotypes_DP(hap, overlap_list, R_INF, g_read->length, force_repeat);
    ///generate_haplotypes_naive(hap, overlap_list, R_INF, g_read->length, force_repeat);

    lable_large_indels(overlap_list, g_read->length, dumy, asm_opt.max_ov_diff_ec);
}

inline void insert_snp_vv(haplotype_evdience_alloc* h, haplotype_evdience* a, uint64_t a_n, char misBase, UC_Read* g_read, void *km)
{
    if(a_n == 0) return;
    SnpStats *p = NULL; uint64_t /**nn = 0,**/ i;
    if(!km) kv_pushp(SnpStats, h->snp_stat, &p);
    else kv_pushp_km(km, SnpStats, h->snp_stat, &p);
    p->id = h->snp_stat.n-1;
    p->occ_0 = 1;
    p->occ_1 = 0;
    p->occ_2 = 0;
    p->overlap_num = 0;
    p->site = a[0].site;
    p->is_homopolymer = if_is_homopolymer_strict(p->site, g_read->seq, g_read->length);
    for (i = 0; i < a_n; i++) {
        if(a[i].type == 0) {
            a[i].overlapSite = p->id;
            h->snp_stat.a[p->id].occ_0 += a[i].cov;
        }
        else if(a[i].type == 1 && a[i].misBase == misBase) {
            a[i].overlapSite = p->id;
            h->snp_stat.a[p->id].occ_1 += a[i].cov;
        }
        else {
            h->snp_stat.a[p->id].occ_2 += a[i].cov;
        }
        h->snp_stat.a[p->id].overlap_num += a[i].cov;
    }
    h->snp_stat.a[p->id].score = -1;
}


int insert_snp_ee(haplotype_evdience_alloc* h, haplotype_evdience* a, uint64_t a_n, haplotype_evdience* u_a, UC_Read* g_read, void *km)
{
    uint64_t i, m, occ_0, occ_1[6], occ_2, diff;
    occ_0 = occ_2 = diff = 0; memset(occ_1, 0, sizeof(uint64_t)*6);

    for (i = 0; i < a_n; i++) {
        if(a[i].type == 0){
            occ_0 += a[i].cov;
        }else if(a[i].type == 1){
            occ_1[seq_nt6_table[(uint8_t)(a[i].misBase)]] += a[i].cov;
            diff += a[i].cov;
        }
        // else if(a[i].type == 2){
        //     occ_2++;
        //     diff++;
        // }
        occ_2 += a[i].cov;
    }

    /**
     1. if occ_0 = 0, that means all overlaps are different with this read at this site
     2. it is not possible that occ_1 = 0,
     3. if occ_1 = 1, there are only one difference. It must be a sequencing error.
        (for repeat, it maybe a snp at repeat. but ...)
    **/
    SnpStats *p = NULL; 
    uint32_t is_homopolymer = (uint32_t)-1; 
    if(occ_0 == 0 || diff <= 1) return 0;
    for (i = m = 0; i < 4; i++) {
        if(occ_1[i] >= 2){
            if(!km) kv_pushp(SnpStats, h->snp_stat, &p);
            else kv_pushp_km(km, SnpStats, h->snp_stat, &p);
            p->id = h->snp_stat.n-1;
            p->occ_0 = 1 + occ_0;
            p->occ_1 = occ_1[i];
            p->occ_2 = occ_2 - p->occ_0 - p->occ_1;
            p->overlap_num = 0;
            p->site = a[0].site;
            p->score = -1;
            p->overlap_num = occ_2;
            if(is_homopolymer == (uint32_t)-1) {
                is_homopolymer = if_is_homopolymer_strict(p->site, g_read->seq, g_read->length);
            }
            p->is_homopolymer = is_homopolymer;
            occ_1[i] = p->id;
            m++;
        } else {
            occ_1[i] = (uint64_t)-1;
        }
    }
    occ_1[4] = occ_1[5] = (uint64_t)-1;
    if(m == 0) return 0;
    
    for (i = m = 0; i < a_n; i++) {
        // fprintf(stderr, "[M::%s] a[%lu].misBase->%c\n", __func__, i, a[i].misBase);
        if(a[i].type == 0) {
            a[i].overlapSite = h->snp_stat.n-1;
        } else if(occ_1[seq_nt6_table[(uint8_t)(a[i].misBase)]]!=(uint64_t)-1){
            a[i].cov = a[i].overlapSite;///note: only renew cov here!!!
            a[i].overlapSite = occ_1[seq_nt6_table[(uint8_t)(a[i].misBase)]];
        } else {
            continue;
        }
        u_a[m++] = a[i];
    }

    /**
    // if(c_snp && ovlp) {
    //     ;
    // }
    for (i = m = 0; i < 4; i++) {
        if(occ_1[i] >= 2) {
            insert_snp_vv(h, a, a_n, s_H[i], g_read, km);
            m++;
        }
    }
    if(m == 0) return 0;

    for (i = m = 0; i < a_n; i++) {
        if(a[i].type == 0){
            u_a[m++] = a[i];
        }else if(a[i].type == 1){
            if(occ_1[seq_nt6_table[(uint8_t)(a[i].misBase)]] >= 2) u_a[m++] = a[i];
        }
    }
    **/

    return m;
}


void partition_overlaps_advance(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, UC_Read* overlap_read, Correct_dumy* dumy, 
                        haplotype_evdience_alloc* hap, int force_repeat)
{
    ResizeInitHaplotypeEvdience(hap);

    uint64_t k, l, m;
    long long window_start, window_end;
    long long num_availiable_win = 0;
    
    Window_Pool w_inf;
    init_Window_Pool(&w_inf, g_read->length, WINDOW, (int)(1.0/asm_opt.max_ov_diff_ec));

    int flag = 0;
    while(get_Window(&w_inf, &window_start, &window_end) && flag != -2)
    {
        dumy->length = 0;
        dumy->lengthNT = 0;
        ///return overlaps that is overlaped with [window_start, window_end]
        flag = get_available_interval(window_start, window_end, overlap_list, dumy);
        switch (flag)
        {
            case 1:    ///found matched overlaps
                break;
            case 0:    ///do not find any matched overlaps
                break;
            case -2: ///do not find any matched overlaps, and the next window also cannot match
                break;
        }
        
        num_availiable_win = num_availiable_win + dumy->length;
        
        ///need to deal with 
        cluster_advance(g_read->seq, window_start, window_end, overlap_list, dumy, R_INF, hap, overlap_read, 1);
    }

    ///very time-consuming
    ///Fix-attention ---> able to be sorted locally
    // qsort(hap->list, hap->length, sizeof(haplotype_evdience), cmp_haplotype_evdience);
    SetSnpMatrix(hap, &(hap->nn_snp), &(overlap_list->length), 0, NULL);
    for (k = 1, l = m = 0; k <= hap->length; ++k) {  
        if (k == hap->length || hap->list[k].site != hap->list[l].site) {
            m += insert_snp_ee(hap, hap->list+l, k-l, hap->list+m, g_read, NULL);
            l = k;
        }
    }
    hap->length = m;

    // generate_haplotypes_naive_advance(hap, overlap_list, NULL);
    generate_haplotypes_naive_HiFi(hap, overlap_list, 0.04, NULL);
    // generate_haplotypes_DP(hap, overlap_list, R_INF, g_read->length, force_repeat);
    // generate_haplotypes_naive(hap, overlap_list, R_INF, g_read->length, force_repeat);

    lable_large_indels(overlap_list, g_read->length, dumy, asm_opt.max_ov_diff_ec);
}


void debug_phasing_snp_site_status(haplotype_evdience_alloc* h, haplotype_evdience* a, uint64_t a_n, overlap_region_alloc* olist)
{
    uint64_t i;
    for (i = 0; i < a_n; i++) {
        if(a[i].site == 76046) {
            fprintf(stderr, "[M::utg%.6dl::]x_site->%u, y_site->%u, type->%u, cov->%u, misBase->%c\n", 
            (int)olist->list[a[i].overlapID].y_id+1, a[i].site, 
            a[i].overlapSite, a[i].type, a[i].cov, a[i].misBase);
        }
    }
}


void partition_ul_overlaps_advance(overlap_region_alloc* overlap_list, const ul_idx_t *uref, 
                        UC_Read* g_read, UC_Read* overlap_read, Correct_dumy* dumy, haplotype_evdience_alloc* hap, 
                        int force_repeat, double max_ov_diff_ec, long long blockLen, void *km)
{
    ResizeInitHaplotypeEvdience(hap);

    uint64_t k, l, m;
    long long window_start, window_end;
    long long num_availiable_win = 0;
    
    Window_Pool w_inf;
    init_Window_Pool(&w_inf, g_read->length, blockLen, (int)(1.0/max_ov_diff_ec));

    int flag = 0;
    while(get_Window(&w_inf, &window_start, &window_end) && flag != -2)
    {
        dumy->length = 0;
        dumy->lengthNT = 0;
        ///return overlaps that is overlaped with [window_start, window_end]
        flag = get_available_interval(window_start, window_end, overlap_list, dumy);
        switch (flag)
        {
            case 1:    ///found matched overlaps
                break;
            case 0:    ///do not find any matched overlaps
                break;
            case -2: ///do not find any matched overlaps, and the next window also cannot match
                break;
        }
        
        num_availiable_win = num_availiable_win + dumy->length;
        
        ///need to deal with 
        cluster_ul_advance(g_read->seq, window_start, window_end, overlap_list, dumy, uref, hap, overlap_read, 1, w_inf.window_length, km);
    }


    
    ///very time-consuming
    ///Fix-attention ---> able to be sorted locally
    // qsort(hap->list, hap->length, sizeof(haplotype_evdience), cmp_haplotype_evdience);
    SetSnpMatrix(hap, &(hap->nn_snp), &(overlap_list->length), 0, km);
    for (k = 1, l = m = 0; k <= hap->length; ++k) {  
        if (k == hap->length || hap->list[k].site != hap->list[l].site) {
            // debug_phasing_snp_site_status(hap, hap->list+l, k-l, overlap_list);
            m += insert_snp_ee(hap, hap->list+l, k-l, hap->list+m, g_read, km);
            l = k;
        }
    }
    hap->length = m;

    generate_haplotypes_naive_advance(hap, overlap_list, km);

    // generate_haplotypes_DP(hap, overlap_list, R_INF, g_read->length, force_repeat);
    // generate_haplotypes_naive(hap, overlap_list, R_INF, g_read->length, force_repeat);

    lable_large_indels(overlap_list, g_read->length, dumy, max_ov_diff_ec);
}


void collect_no_cov_regions(overlap_region_alloc* overlap_list, All_reads* R_INF, 
kvec_t_u32_warp* b, kvec_t_u64_warp* r, int min_dp, int min_len)
{
    b->a.n = r->a.n = 0;
    ///if(overlap_list->length == 0) return;
    long long i = 0, xLen = Get_READ_LENGTH((*R_INF), overlap_list->list[0].x_id);
    uint32_t qs, qe;
    uint64_t tmp;
    int dp, old_dp, s_start = 0, s_end = 0;
    ///at least 1
    if(min_len < 1) min_len = 1;

    
    for (i = 0; i < (long long)overlap_list->length; i++)
    {
        if (overlap_list->list[i].is_match != 1 && overlap_list->list[i].is_match != 2) continue;

        qs = overlap_list->list[i].x_pos_s;
        qe = overlap_list->list[i].x_pos_e + 1;
        kv_push(uint32_t, b->a, qs<<1);
        kv_push(uint32_t, b->a, qe<<1|1);
    }


    ///we can identify the qs and qe by the 0-th bit
    radix_sort_b32(b->a.a, b->a.a + b->a.n);

    for (i = 0, dp = 0; i < (long long)b->a.n; ++i) 
    {
        old_dp = dp;
        //if a[j] is qe
        if (b->a.a[i]&1) --dp;
        else ++dp;
        /**
        min_dp is the coverage drop threshold
        there are two cases: 
            1. old_dp = dp + 1 (b.a[j] is qe); 2. old_dp = dp - 1 (b.a[j] is qs);
        **/
        if (old_dp < min_dp && dp >= min_dp) ///old_dp < dp, b.a[j] is qs
        { 
            ///case 2, a[j] is qs
            s_end = b->a.a[i]>>1;
            ///at least 1
            if(s_end-s_start >= min_len)
            {
                tmp = s_start; tmp = tmp << 32; tmp = tmp | (uint64_t)(s_end-1);
                kv_push(uint64_t, r->a, tmp);
            }
        }
        else if (old_dp >= min_dp && dp < min_dp) ///old_dp > min_dp, b.a[j] is qe
        {
            s_start = b->a.a[i]>>1;
        }
    }


    if(s_start < xLen && xLen-s_start >= min_len)
    {
        s_end = xLen;
        tmp = s_start; tmp = tmp << 32; tmp = tmp | (uint64_t)(s_end-1);
        kv_push(uint64_t, r->a, tmp);
    }
}

int collect_hp_regions_back(overlap_region_alloc* olist, All_reads* R_INF, kvec_t_u32_warp* b, kvec_t_u64_warp* r, kvec_t_u8_warp* k_flag, float hp_rate, FILE* fp)
{
    int i, k, qs, qe, ava_k_mer = 0, hp_k_mer = 0, min_dp;
    // min_dp = RESEED_DP;
    // if(asm_opt.hom_cov > 0) min_dp = asm_opt.hom_cov * RESEED_PEAK_RATE;
    // if(min_dp > RESEED_DP) min_dp = RESEED_DP;
    min_dp = RESEED_DP;
    if(asm_opt.hom_cov > 0) min_dp = asm_opt.hom_cov * RESEED_PEAK_RATE;
    if(asm_opt.het_cov > 0) min_dp = asm_opt.het_cov * RESEED_PEAK_RATE;
    collect_no_cov_regions(olist, R_INF, b, r, min_dp, RESEED_LEN);

    for (i = 0; i < (int)r->a.n; i++)
    {
        ///[qs, qe]
        qs = r->a.a[i]>>32;
        qe = (r->a.a[i]<<32)>>32;
        
        for (k = qs; k <= qe; k++)
        {
            if(k_flag->a.a[k] > 1) ava_k_mer++;
            if(k_flag->a.a[k] > 2) hp_k_mer++;
        }

        if(fp) fprintf(fp, "qs: %d, qe: %d, ava_k_mer: %d, hp_k_mer: %d\n", qs, qe, ava_k_mer, hp_k_mer);        
    }

    if(fp) fprintf(fp, "ava_k_mer: %d, hp_k_mer: %d, hp_rate: %f, min_dp: %d, a.n: %d\n", ava_k_mer, hp_k_mer, hp_rate, min_dp, (int)r->a.n);

    // if(fp)
    // {
    //     for (k = 0; k < (int)k_flag->a.n; k++)
    //     {
    //         if(k_flag->a.a[k] > 0) fprintf(fp, "(%d) %u\n", k, k_flag->a.a[k]);
    //     }
    // }
    
    if(hp_k_mer > ava_k_mer*hp_rate) return 1; ///must use '>' instead of '>='
    r->a.n = 0;
    return 0;
}


int collect_hp_regions(overlap_region_alloc* olist, All_reads* R_INF, kvec_t_u8_warp* k_flag, 
float hp_rate, int rlen, FILE* fp)
{
    int i, ava_k_mer = 0, hp_k_mer = 0, vLen, min_dp;
    int32_t w, n[4];
    n[0] = n[1] = n[2] = n[3] = 0;
    min_dp = RESEED_DP;
    if(asm_opt.hom_cov > 0) min_dp = asm_opt.hom_cov * RESEED_PEAK_RATE;
    if(min_dp > RESEED_DP) min_dp = RESEED_DP;
    overlap_region* ov = NULL;

    for (i = 0; i < (long long)olist->length; i++)
    {
        ov = &(olist->list[i]);
        if (ov->is_match != 1 && ov->is_match != 2) continue;

        w = ha_ov_type(ov, rlen);
		++n[w];
    }

    if(fp) fprintf(fp, "n[0]: %d, n[1]: %d, n[2]: %d, n[3]: %d\n", n[0], n[1], n[2], n[3]);

    // n[0] += n[2];
    // n[1] += n[2];

    if(n[0] < min_dp)
    {
        ava_k_mer = hp_k_mer = 0;
        vLen = MIN(k_flag->a.n, RESEED_LEN);
        for (i = 0; i < vLen; i++)
        {
            if(k_flag->a.a[i] > 1) ava_k_mer++;
            if(k_flag->a.a[i] > 2) hp_k_mer++;
        }
        if(hp_k_mer > ava_k_mer*hp_rate) return 1;
    }


    if(n[1] < min_dp)
    {
        ava_k_mer = hp_k_mer = 0;
        vLen = MIN(k_flag->a.n, RESEED_LEN);
        for (i = k_flag->a.n - vLen; i < (int)k_flag->a.n; i++)
        {
            if(k_flag->a.a[i] > 1) ava_k_mer++;
            if(k_flag->a.a[i] > 2) hp_k_mer++;
        }
        if(hp_k_mer > ava_k_mer*hp_rate) return 1;
    }


    if(fp) fprintf(fp, "ava_k_mer: %d, hp_k_mer: %d, hp_rate: %f, min_dp: %d\n", ava_k_mer, hp_k_mer, hp_rate, min_dp);

    // if(fp)
    // {
    //     for (k = 0; k < (int)k_flag->a.n; k++)
    //     {
    //         if(k_flag->a.a[k] > 0) fprintf(fp, "(%d) %u\n", k, k_flag->a.a[k]);
    //     }
    // }
    
    return 0;
}

uint64_t ovlp_occ(overlap_region_alloc* overlap_list, uint8_t is_match)
{
    uint64_t occ = 0, k;
    for (k = 0; k < overlap_list->length; k++) {
        if(overlap_list->list[k].is_match == is_match) occ++;
    }
    return occ;
}

void correct_overlap(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read, 
                        Graph* g, Graph* DAGCon, Cigar_record* current_cigar, 
                        haplotype_evdience_alloc* hap, Round2_alignment* second_round, 
                        kvec_t_u64_warp* v_idx, window_list_alloc* win_ciagr_buf, int force_repeat, int is_consensus, int* fully_cov, int* abnormal)
{
    clear_Correct_dumy(dumy, overlap_list, NULL);

    long long window_start, window_end;

    Window_Pool w_inf;

    init_Window_Pool(&w_inf, g_read->length, WINDOW, (int)(1.0/asm_opt.max_ov_diff_ec));

    int flag = 0;

    while(get_Window(&w_inf, &window_start, &window_end) && flag != -2)
    {
        dumy->length = 0;
        dumy->lengthNT = 0;
        flag = get_interval(window_start, window_end, overlap_list, dumy, w_inf.window_length);

        switch (flag)
        {
            case 1:    ///no match here
                break;
            case 0:    ///no match here
                break;
            case -2: ///if flag == -2, loop would be terminated
                break;
        }

        ///dumy->lengthNT represent how many overlaps that the length of them is not equal to WINDOW; may larger or less than WINDOW
        ///dumy->length represent how many overlaps that the length of them is WINDOW
        /****************************may improve**************************/
        ///now the windows which are larger than WINDOW are verified one-by-one, to improve it, we can do it group-bygroup
        verify_window(window_start, window_end, overlap_list, dumy, R_INF, g_read->seq);
    }
    // fprintf(stderr, "###dumy->start_i:%lu, overlap_list->length:%lu\n\n", dumy->start_i, overlap_list->length);
    
    // recalcate_window(overlap_list, R_INF, g_read, dumy, overlap_read);
    // partition_overlaps(overlap_list, R_INF, g_read, dumy, hap, force_repeat);
    recalcate_window_advance(overlap_list, R_INF, NULL, g_read, dumy, overlap_read, v_idx, w_inf.window_length, asm_opt.max_ov_diff_ec, asm_opt.max_ov_diff_final);
    // fprintf(stderr, "[M::%s-beg] occ[0]->%lu, occ[1]->%lu, occ[2]->%lu, occ[3]->%lu\n", __func__, 
    // ovlp_occ(overlap_list, 0), ovlp_occ(overlap_list, 1), ovlp_occ(overlap_list, 2), ovlp_occ(overlap_list, 3));
    partition_overlaps_advance(overlap_list, R_INF, g_read, overlap_read, dumy, hap, force_repeat);
    // fprintf(stderr, "[M::%s-after] occ[0]->%lu, occ[1]->%lu, occ[2]->%lu, occ[3]->%lu\n", __func__, 
    // ovlp_occ(overlap_list, 0), ovlp_occ(overlap_list, 1), ovlp_occ(overlap_list, 2), ovlp_occ(overlap_list, 3));
    if(is_consensus)
    {
        generate_consensus(overlap_list, R_INF, g_read, dumy, g, DAGCon, current_cigar, second_round, win_ciagr_buf);
    }

    (*fully_cov) = check_if_fully_covered(overlap_list, R_INF, g_read, dumy, g, abnormal);
}

/**
void debug_phasing_status(overlap_region_alloc *olist, ma_ug_t *ug, uint64_t print_w_list,
haplotype_evdience_alloc* hap, UC_Read* g_read, int64_t flanking, uint64_t yid)
{
    uint64_t i, k, l, ii;
    int64_t t;
    SnpStats *s = NULL;
    
        for (i = 0; i < olist->length; i++) {
            if(olist->list[i].y_id != yid) continue;
            fprintf(stderr, "\n[M::utg%.6d%c::is_match->%u] rev->%u, x->[%u, %u), y->[%u, %u)\n", 
            (int)olist->list[i].y_id+1, "lc"[ug->u.a[olist->list[i].y_id].circ], olist->list[i].is_match, 
                olist->list[i].y_pos_strand, olist->list[i].x_pos_s, olist->list[i].x_pos_e+1, olist->list[i].y_pos_s, olist->list[i].y_pos_e+1);
            if(print_w_list) {
                for (k = 0; k < olist->list[i].w_list.n; k++) {
                    if(olist->list[i].w_list.a[k].y_end != -1) {
                        fprintf(stderr, "x->[%lu, %lu), y->[%d, %d), e->%d\n", 
                        olist->list[i].w_list.a[k].x_start, olist->list[i].w_list.a[k].x_end+1, 
                            olist->list[i].w_list.a[k].y_start, olist->list[i].w_list.a[k].y_end+1, 
                            olist->list[i].w_list.a[k].error);
                    } else {
                        fprintf(stderr, "x->[-1, -1), y->[-1, -1), e->-1\n");
                    }
                }
            }
        }
    

    for (k = 1, l = 0; k <= hap->length; ++k) {  
        if (k == hap->length || hap->list[k].overlapID != hap->list[l].overlapID) {
            ii = hap->list[l].overlapID;
            if(olist->list[ii].y_id != yid) {
                l = k;
                continue;
            }
            for (i = l; i < k; i++) {
                if(hap->list[i].type!=1) continue; 
                s = &(hap->snp_stat.a[hap->list[i].overlapSite]);
                if(s->score == 1 && (!(s->occ_0 < 2 || s->occ_1 < 2))) {
                    fprintf(stderr, "s->site:%u, s->occ_0:%u, s->occ_1:%u, s->occ_2:%u\n", s->site, s->occ_0, s->occ_1, s->occ_2);
                    for (t = s->site>=flanking?s->site-flanking:0; t<g_read->length && t<=s->site+flanking; t++){
                        if(t == s->site) fprintf(stderr,"[");
                        fprintf(stderr,"%c", g_read->seq[t]);
                        if(t == s->site) fprintf(stderr,"]");
                    }
                    fprintf(stderr,"\n");
                }
            }
            l = k;
        }
    }
}
**/

void print_ovlp_occ_stat(overlap_region_alloc* overlap_list, uint32_t xlen, uint8_t is_match)
{
    uint64_t k;
    for (k = 0; k < overlap_list->length; k++) {
        if(overlap_list->list[k].is_match != is_match) continue;
        fprintf(stderr, "[M::%s::xlen::%u] utg%.6dl(%c), is_match::%u, x::[%u, %u)\n", __func__, xlen,
            (int32_t)overlap_list->list[k].y_id + 1, 
            "+-"[overlap_list->list[k].y_pos_strand], overlap_list->list[k].is_match,
            overlap_list->list[k].x_pos_s, overlap_list->list[k].x_pos_e+1);
    }
}
void align_ul_ed(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, char* qstr, char *tstr, double e_rate, int64_t w_l, void *km);
void correct_ul_overlap(overlap_region_alloc* overlap_list, const ul_idx_t *uref, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read, 
                        Graph* g, Graph* DAGCon, Cigar_record* current_cigar, 
                        haplotype_evdience_alloc* hap, Round2_alignment* second_round, 
                        kvec_t_u64_warp* v_idx, window_list_alloc* win_ciagr_buf, 
                        int force_repeat, int is_consensus, int* fully_cov, int* abnormal, 
                        double max_ov_diff_ec, long long winLen, void *km)
{
    
    clear_Correct_dumy(dumy, overlap_list, km);

    Window_Pool w_inf;
    init_Window_Pool(&w_inf, g_read->length, winLen, (int)(1.0/max_ov_diff_ec));
    /**
    uint64_t i;
    for (i = 0; i < overlap_list->length; i++) {
        verify_ul_window_s(&(overlap_list->list[i]), uref, g_read->seq, dumy->overlap_region, max_ov_diff_ec, w_inf.window_length, THRESHOLD_MAX_SIZE, km);
        // align_ul_ed(&(overlap_list->list[i]), uref, NULL, g_read->seq, dumy->overlap_region, max_ov_diff_ec, w_inf.window_length, km);
    }
    **/
   
    
    // recalcate_window(overlap_list, R_INF, g_read, dumy, overlap_read);
    // partition_overlaps(overlap_list, R_INF, g_read, dumy, hap, force_repeat);
    // recalcate_window_ul_advance(overlap_list, uref, g_read, dumy, overlap_read, max_ov_diff_ec, w_inf.window_length, km);
    // recalcate_window_advance(overlap_list, NULL, uref, g_read, dumy, overlap_read, v_idx, w_inf.window_length, max_ov_diff_ec, max_ov_diff_ec);
    /**
    refine_ed_aln(overlap_list, NULL, uref, g_read, dumy, overlap_read, v_idx, w_inf.window_length, max_ov_diff_ec, max_ov_diff_ec);
    **/
    refine_ed_aln_test(overlap_list, NULL, uref, g_read, dumy, overlap_read, v_idx, w_inf.window_length, max_ov_diff_ec, max_ov_diff_ec);
    // fprintf(stderr, "[M::%s-beg] occ[0]->%lu, occ[1]->%lu, occ[2]->%lu, occ[3]->%lu\n", __func__, 
    // ovlp_occ(overlap_list, 0), ovlp_occ(overlap_list, 1), ovlp_occ(overlap_list, 2), ovlp_occ(overlap_list, 3));
    ///after this function, overlap_list is sorted by x_pos_e; used for g_chain
    partition_ul_overlaps_advance(overlap_list, uref, g_read, overlap_read, dumy, hap, force_repeat, max_ov_diff_ec, w_inf.window_length, km);
    // print_ovlp_occ_stat(overlap_list, g_read->length, 1);
    // print_ovlp_occ_stat(overlap_list, g_read->length, 2);
    // fprintf(stderr, "[M::%s-end] occ[0]->%lu, occ[1]->%lu, occ[2]->%lu, occ[3]->%lu\n", __func__, 
    // ovlp_occ(overlap_list, 0), ovlp_occ(overlap_list, 1), ovlp_occ(overlap_list, 2), ovlp_occ(overlap_list, 3));
    // debug_phasing_status(overlap_list, uref->ug, 0, hap, g_read, 20, 1176);
    // debug_phasing_status(overlap_list, uref->ug, 0, hap, g_read, 20, 1167);
    // debug_phasing_status(overlap_list, uref->ug, 0, hap, g_read, 20, 1170);
    /**
    
    
    if(is_consensus)
    {
        generate_consensus(overlap_list, R_INF, g_read, dumy, g, DAGCon, current_cigar, second_round);
    }


    (*fully_cov) = check_if_fully_covered(overlap_list, R_INF, g_read, dumy, g, abnormal);
    **/
}




void lchain_align(overlap_region_alloc* overlap_list, const ul_idx_t *uref, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read, 
                        Graph* g, Graph* DAGCon, Cigar_record* current_cigar, 
                        haplotype_evdience_alloc* hap, Round2_alignment* second_round, 
                        kvec_t_u64_warp* v_idx, window_list_alloc* win_ciagr_buf, 
                        int force_repeat, int is_consensus, int* fully_cov, int* abnormal, 
                        double max_ov_diff_ec, long long winLen, void *km)
{
    
    clear_Correct_dumy(dumy, overlap_list, km);

    long long window_start, window_end;

    Window_Pool w_inf;

    init_Window_Pool(&w_inf, g_read->length, winLen, (int)(1.0/max_ov_diff_ec));
    int flag = 0;
    
    while(get_Window(&w_inf, &window_start, &window_end) && flag != -2) {
        dumy->length = 0; dumy->lengthNT = 0;
        flag = get_interval(window_start, window_end, overlap_list, dumy, w_inf.window_length);

        switch (flag) {
            case 1:    ///no match here
                break;
            case 0:    ///no match here
                break;
            case -2: ///if flag == -2, loop would be terminated
                break;
        }

        ///dumy->lengthNT represent how many overlaps that the length of them is not equal to WINDOW; may larger or less than WINDOW
        ///dumy->length represent how many overlaps that the length of them is WINDOW
        ///now the windows which are larger than WINDOW are verified one-by-one, to improve it, we can do it group-bygroup
        verify_ul_window(window_start, window_end, overlap_list, dumy, uref, g_read->seq, max_ov_diff_ec, w_inf.window_length, /**THRESHOLD**/THRESHOLD_MAX_SIZE, km);
    }
    
    // recalcate_window(overlap_list, R_INF, g_read, dumy, overlap_read);
    // partition_overlaps(overlap_list, R_INF, g_read, dumy, hap, force_repeat);
    // recalcate_window_ul_advance(overlap_list, uref, g_read, dumy, overlap_read, max_ov_diff_ec, w_inf.window_length, km);
    refine_ed_aln(overlap_list, NULL, uref, g_read, dumy, overlap_read, v_idx, w_inf.window_length, max_ov_diff_ec, max_ov_diff_ec);
    // fprintf(stderr, "[M::%s-beg] occ[0]->%lu, occ[1]->%lu, occ[2]->%lu, occ[3]->%lu\n", __func__, 
    // ovlp_occ(overlap_list, 0), ovlp_occ(overlap_list, 1), ovlp_occ(overlap_list, 2), ovlp_occ(overlap_list, 3));
    ///after this function, overlap_list is sorted by x_pos_e; used for g_chain
    partition_ul_overlaps_advance(overlap_list, uref, g_read, overlap_read, dumy, hap, force_repeat, max_ov_diff_ec, w_inf.window_length, km);
    // print_ovlp_occ_stat(overlap_list, g_read->length, 1);
    // print_ovlp_occ_stat(overlap_list, g_read->length, 2);
    // fprintf(stderr, "[M::%s-end] occ[0]->%lu, occ[1]->%lu, occ[2]->%lu, occ[3]->%lu\n", __func__, 
    // ovlp_occ(overlap_list, 0), ovlp_occ(overlap_list, 1), ovlp_occ(overlap_list, 2), ovlp_occ(overlap_list, 3));
    // debug_phasing_status(overlap_list, uref->ug, 0, hap, g_read, 20, 1176);
    // debug_phasing_status(overlap_list, uref->ug, 0, hap, g_read, 20, 1167);
    // debug_phasing_status(overlap_list, uref->ug, 0, hap, g_read, 20, 1170);
    /**
    
    
    if(is_consensus)
    {
        generate_consensus(overlap_list, R_INF, g_read, dumy, g, DAGCon, current_cigar, second_round);
    }


    (*fully_cov) = check_if_fully_covered(overlap_list, R_INF, g_read, dumy, g, abnormal);
    **/
}

void init_Cigar_record(Cigar_record* dummy)
{
    dummy->length = 0;
    dummy->size = 100;
    dummy->record = (uint32_t*)malloc(sizeof(uint32_t)*dummy->size);


    dummy->lost_base_length = 0;
    dummy->lost_base_size = 100;
    dummy->lost_base = (char*)malloc(sizeof(char)*dummy->lost_base_size);



    dummy->current_operation_length = 0;
    dummy->current_operation = 127;
}

void init_Cigar_record_buf(Cigar_record* dummy, void *km)
{
    memset(dummy, 0, sizeof(*dummy));
    dummy->current_operation = 127;
}


void destory_Cigar_record(Cigar_record* dummy)
{
    free(dummy->record);
    free(dummy->lost_base);
}

void clear_Cigar_record(Cigar_record* dummy)
{
    dummy->new_read_length = 0;
    dummy->length = 0;
    dummy->lost_base_length = 0;

    dummy->current_operation_length = 0;
    dummy->current_operation = 127;
}


void init_Correct_dumy_buf(Correct_dumy* list, void *km)
{
    memset(list, 0, sizeof(Correct_dumy));
    int i;
    for (i = 0; i < 256; i++){
		list->Peq_SSE[i] = _mm_setzero_si128();
	}
}


void init_Correct_dumy(Correct_dumy* list)
{
    list->size = 0;
    list->length = 0;
    list->lengthNT = 0;
    list->start_i = 0;
    list->overlapID = NULL;
    int i;
    for (i = 0; i < 256; i++)
	{
		list->Peq_SSE[i] = _mm_setzero_si128();
	}


    list->corrected_read_size = 1000;
    list->corrected_read_length = 0;
    list->corrected_read = (char*)malloc(sizeof(char)*list->corrected_read_size);
    list->corrected_base = 0;

}

void destory_Correct_dumy(Correct_dumy* list)
{
    free(list->overlapID);
    free(list->corrected_read);
}


void clear_Correct_dumy(Correct_dumy* list, overlap_region_alloc* overlap_list, void *km)
{
    list->length = 0;
    list->lengthNT = 0;
    list->start_i = 0;

    if (list->size < overlap_list->length){
        list->size = overlap_list->length;
        if(!km) REALLOC(list->overlapID, list->size);
        else KREALLOC(km, list->overlapID, list->size);
    }

    list->last_boundary_length = 0;
    list->corrected_read_length = 0;
    list->corrected_base = 0;
    
}

void clear_Correct_dumy_pure(Correct_dumy* list)
{
    list->length = 0;
    list->lengthNT = 0;
    list->start_i = 0;
    list->last_boundary_length = 0;
    list->corrected_read_length = 0;
    list->corrected_base = 0;
}


void init_Cigar_record_alloc(Cigar_record_alloc* x)
{
    x->length = 0;
    x->size = 0;
    x->buffer = NULL;
}

void resize_Cigar_record_alloc(Cigar_record_alloc* x, long long new_size)
{
    long long i;
    if(new_size > x->size)
    {
        x->buffer = (Cigar_record*)realloc(x->buffer, new_size*sizeof(Cigar_record));
        for (i = 0; i < x->size; i++)
        {
            clear_Cigar_record(&(x->buffer[i]));
        }
        for (; i < new_size; i++)
        {
            init_Cigar_record(&(x->buffer[i]));
            clear_Cigar_record(&(x->buffer[i]));
        }

        x->size = new_size;
    }
    else
    {
        for (i = 0; i < new_size; i++)
        {
            clear_Cigar_record(&(x->buffer[i]));
        }
    }

    x->length = 0;
}
void destory_Cigar_record_alloc(Cigar_record_alloc* x)
{
    long long i;
    for (i = 0; i < x->size; i++)
    {
        destory_Cigar_record(&(x->buffer[i]));
    }
    free(x->buffer);
}


void add_new_cell_to_cigar_record(Cigar_record* dummy, uint32_t len, uint32_t type)
{
    uint32_t tmp;
    tmp = len;
    tmp = tmp << 2;
    tmp = tmp | type;

    dummy->length++;


    if(dummy->length > dummy->size)
    {
        dummy->size = dummy->size * 2;
        dummy->record = (uint32_t*)realloc(dummy->record, dummy->size*sizeof(uint32_t));
    }

    dummy->record[dummy->length - 1] = tmp;
}

void add_existing_cell_to_cigar_record(Cigar_record* dummy, uint32_t len, uint32_t type)
{
    uint32_t tmp;

    tmp = dummy->record[dummy->length - 1] >> 2;
    tmp = tmp + len;
    tmp = tmp << 2;
    tmp = tmp | type;
    dummy->record[dummy->length - 1] = tmp;
}


void add_new_cell_to_cigar_record_with_different_base(Cigar_record* dummy, uint32_t len, uint32_t type, char* seq)
{
    uint32_t tmp;
    tmp = len;
    tmp = tmp << 2;
    tmp = tmp | type;
   

    dummy->length++;

    if(dummy->length > dummy->size)
    {
        dummy->size = dummy->size * 2;
        dummy->record = (uint32_t*)realloc(dummy->record, dummy->size*sizeof(uint32_t));
    }

    dummy->record[dummy->length - 1] = tmp;

    

    if (dummy->lost_base_length + len> dummy->lost_base_size)
    {
        dummy->lost_base_size = (dummy->lost_base_length + len) * 2;
        dummy->lost_base = (char*)realloc(dummy->lost_base, dummy->lost_base_size*sizeof(char));
    }

    uint32_t i = 0;
    for (i = 0; i < len; i++, dummy->lost_base_length++)
    {
        dummy->lost_base[dummy->lost_base_length] = seq[i];
    }
}

void add_existing_cell_to_cigar_record_with_different_base(Cigar_record* dummy, uint32_t len, uint32_t type, char* seq)
{
    uint32_t tmp;

    tmp = dummy->record[dummy->length - 1] >> 2;
    tmp = tmp + len;
    tmp = tmp << 2;
    tmp = tmp | type;
    dummy->record[dummy->length - 1] = tmp;

    if (dummy->lost_base_length + len> dummy->lost_base_size)
    {
        dummy->lost_base_size = (dummy->lost_base_length + len) * 2;
        dummy->lost_base = (char*)realloc(dummy->lost_base, dummy->lost_base_size*sizeof(char));
    }

    uint32_t i = 0;
    for (i = 0; i < len; i++, dummy->lost_base_length++)
    {
        dummy->lost_base[dummy->lost_base_length] = seq[i];
    }
}

void afine_gap_alignment(const char *qseq, uint8_t* qnum, const int ql, 
const char *tseq, uint8_t* tnum, const int tl, const uint8_t *c2n, const int strand,
int sc_mch, int sc_mis, int gapo, int gape, int bandLen, int zdrop, int end_bonus,
long long* max_q_pos, long long* max_t_pos, long long* global_score, 
long long* extention_score, long long* q_boundary_score, long long* q_boundary_t_coordinate,
long long* t_boundary_score, long long* t_boundary_q_coordinate,
long long* droped, int mode)
{
    /**
    // for ksw2
    (*max_t_pos) = (*max_q_pos) = -1;
	int i, a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
    
	int8_t mat[25] = {(int8_t)a,(int8_t)b,(int8_t)b,(int8_t)b,0, 
    (int8_t)b,(int8_t)a,(int8_t)b,(int8_t)b,0, (int8_t)b,(int8_t)b,(int8_t)a,(int8_t)b,0, 
    (int8_t)b,(int8_t)b,(int8_t)b,(int8_t)a,0, 0,0,0,0,0};
	ksw_extz_t ez;
	memset(&ez, 0, sizeof(ksw_extz_t));

    if(strand == FORWARD_KSW)
    {
         for (i = 0; i < tl; ++i) tnum[i] = c2n[(uint8_t)tseq[i]]; // encode to 0/1/2/3
	     for (i = 0; i < ql; ++i) qnum[i] = c2n[(uint8_t)qseq[i]];
    }
    else if(strand == BACKWARD_KSW)
    {
         for (i = 0; i < tl; ++i) tnum[i] = c2n[(uint8_t)tseq[tl - i - 1]]; // encode to 0/1/2/3
	     for (i = 0; i < ql; ++i) qnum[i] = c2n[(uint8_t)qseq[ql - i - 1]];
    }
    ksw_extz2_sse(0, ql, qnum, tl, tnum, 5, mat, gapo, gape, bandLen, zdrop, end_bonus, 
    mode, &ez);
    (*global_score) = ez.score;
    (*extention_score) = ez.max;
    (*q_boundary_score) = ez.mqe;
    (*q_boundary_t_coordinate) = ez.mqe_t;
    (*t_boundary_score) = ez.mte;
    (*t_boundary_q_coordinate) = ez.mte_q;
    (*max_t_pos) = ez.max_t;
    (*max_q_pos) = ez.max_q;
    (*droped) = ez.zdropped;
    free(ez.cigar);
    
    
	// for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
	// 	printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
	// putchar('\n');
    
   // for ksw2
   **/
}

int fill_chain_by_affine_gap_debug(Fake_Cigar* chain, char* x_string, char* y_string, overlap_region* ovc,
long long x_readLen, long long y_readLen, Cigar_record* cigar, uint8_t* c2n, uint8_t* x_num, uint8_t* y_num,
long long* minus_score_thres, long long* final_scores)
{
    /**
    long long i, xOffset, yOffset, xRegionLen, yRegionLen, maxXpos, maxYpos, zdroped;
    long long mapGlobalScore, mapExtentScore;
    long long xBuoundaryScore, xBuoundaryYcoordinate, yBuoundaryScore, yBuoundaryXcoordinate;
    ///float band_rate = 0.08;
    int endbouns, mode;
    long long xBeg, yBeg;
    xBeg = ovc->x_pos_s;
    yBeg = ovc->y_pos_s;
    if(chain->length <= 0) return 0;
    // long long minus_score_thres = (EstimateOlen*HIGH_HET_ERROR_RATE*(MATCH_SCORE_KSW+(MAX(MISMATCH_SCORE_KSW,GAP_EXT_KSW))));
    // long long total_score_thres = EstimateOlen*MATCH_SCORE_KSW - minus_score_thres;
    long long sum_score = 0, current_ovlp = 0, zdrop_occ = 0;
    long long new_xBeg, new_yBeg, new_xEnd, new_yEnd;
    new_xBeg = ovc->x_pos_s;
    new_yBeg = ovc->y_pos_s;
    new_xEnd = ovc->x_pos_e;
    new_yEnd = ovc->y_pos_e;
    ///long long sub_score_sum;
    ///deal with region 0 backward
    i = 0;
    endbouns = 0;

    xOffset = get_fake_gap_pos(chain, 0);
    xOffset = xOffset - 1;
    yOffset = (xOffset - xBeg) + yBeg + get_fake_gap_shift(chain, 0);
    if(xOffset >= 0 && yOffset >= 0)
    {
        xRegionLen = xOffset + 1;
        yRegionLen = yOffset + 1;
        //note here cannot use DIFF(xRegionLen, yRegionLen)
        // bandLen = (MIN(xRegionLen, yRegionLen))*band_rate;
        // if(bandLen == 0) bandLen = MIN(xRegionLen, yRegionLen);


        ///do alignment backward
        ///for beginning part and end part, must use exact mode
        mode = KSW_EZ_SCORE_ONLY;
        afine_gap_alignment(x_string, x_num, xRegionLen, y_string, y_num, yRegionLen, 
        c2n, BACKWARD_KSW, MATCH_SCORE_KSW, MISMATCH_SCORE_KSW, GAP_OPEN_KSW, GAP_EXT_KSW,
        BAND_KSW, Z_DROP_KSW, endbouns, &maxXpos, &maxYpos, &mapGlobalScore, 
        &mapExtentScore, &xBuoundaryScore, &xBuoundaryYcoordinate, 
        &yBuoundaryScore, &yBuoundaryXcoordinate, &zdroped, mode);

        
        

        if(!zdroped)
        {
            if(xRegionLen <= yRegionLen)
            {
                sum_score += xBuoundaryScore;
                new_yBeg = yRegionLen - xBuoundaryYcoordinate - 1;
            }
            else
            {
                sum_score += yBuoundaryScore;
                new_xBeg = xRegionLen - yBuoundaryXcoordinate - 1;
            }
        }
        else
        {   ///return 0;
            sum_score += mapExtentScore;
            if(xRegionLen <= yRegionLen)
            {
                
                sum_score -= (GAP_OPEN_KSW + (xRegionLen - maxXpos)*GAP_EXT_KSW);
            }
            else
            {
                sum_score -= (GAP_OPEN_KSW + (yRegionLen - maxYpos)*GAP_EXT_KSW);
            }
        }   
    }

    ///align forward
    for (i = 0; i < (long long)chain->length; i++)
    {
        // xOffset = get_fake_gap_pos(chain, i);
        // yOffset = xOffset + get_fake_gap_shift(chain, i);
        xOffset = get_fake_gap_pos(chain, i);
        yOffset = (xOffset - xBeg) + yBeg + get_fake_gap_shift(chain, i);
        ///last region
        if(i == (long long)(chain->length - 1))
        {
            endbouns = 0;
            xRegionLen = x_readLen - xOffset;
            yRegionLen = y_readLen - yOffset;
            ///for beginning part and end part, must use exact mode
            mode = KSW_EZ_SCORE_ONLY;
            //note here cannot use DIFF(xRegionLen, yRegionLen)
            // bandLen = (MIN(xRegionLen, yRegionLen))*band_rate;
            // if(bandLen == 0) bandLen = MIN(xRegionLen, yRegionLen);
        }
        else
        {
            ///higher endbouns for middle regions
            endbouns = MATCH_SCORE_KSW;
            xRegionLen = get_fake_gap_pos(chain, i+1) - xOffset;
            yRegionLen = (get_fake_gap_pos(chain, i+1) + get_fake_gap_shift(chain, i+1)) - 
                         (get_fake_gap_pos(chain, i) + get_fake_gap_shift(chain, i));
            mode = KSW_EZ_SCORE_ONLY | KSW_EZ_APPROX_MAX | KSW_EZ_APPROX_DROP;
            // bandLen = MAX((MIN(xRegionLen, yRegionLen))*band_rate, DIFF(xRegionLen, yRegionLen));
            // if(bandLen == 0) bandLen = MIN(xRegionLen, yRegionLen);
        }


        if(minus_score_thres)
        {
            current_ovlp = MIN((xOffset + 1 - xBeg), (yOffset + 1 - yBeg));
            current_ovlp = current_ovlp*MATCH_SCORE_KSW;
            if(current_ovlp - sum_score > (*minus_score_thres))
            {
                return 0;
            }
        }

        if(xOffset < 0) xOffset = 0;
        if(yOffset < 0) yOffset = 0;
        if(xRegionLen < 0) xRegionLen = 0;
        if(yRegionLen < 0) yRegionLen = 0;

        ///do alignment forward
        ///text is x, query is y
        afine_gap_alignment(x_string+xOffset, x_num, xRegionLen, y_string+yOffset, y_num, yRegionLen, 
        c2n, FORWARD_KSW, MATCH_SCORE_KSW, MISMATCH_SCORE_KSW, GAP_OPEN_KSW, GAP_EXT_KSW,
        BAND_KSW, Z_DROP_KSW, endbouns, &maxXpos, &maxYpos, &mapGlobalScore, 
        &mapExtentScore, &xBuoundaryScore, &xBuoundaryYcoordinate, 
        &yBuoundaryScore, &yBuoundaryXcoordinate, &zdroped, mode);
        // fprintf(stderr, "# xOffset: %lld, yOffset: %lld, xRegionLen: %lld, yRegionLen: %lld, bandLen: %lld, maxXpos: %lld, maxYpos: %lld, zdroped: %lld\n",
        // xOffset, yOffset, xRegionLen, yRegionLen, BAND_KSW, maxXpos, maxYpos, zdroped);

        if(!zdroped)
        {
            if(i != (long long)(chain->length - 1))
            {
                sum_score += mapGlobalScore;
            }
            else
            {
                if(xRegionLen <= yRegionLen)
                {
                    sum_score += xBuoundaryScore;
                    new_yEnd = yOffset + xBuoundaryYcoordinate;
                }
                else
                {
                    sum_score += yBuoundaryScore;
                    new_xEnd = xOffset + yBuoundaryXcoordinate;
                    // if(new_xEnd != (long long)ovc->x_pos_e)
                    // {
                    //     fprintf(stderr, "\n******direction: %u, new_xBeg: %lld, new_xEnd: %lld, new_yBeg: %lld, new_yEnd: %lld, old_xBeg: %u, old_xEnd: %u, old_yBeg: %u, old_yEnd: %u\n", 
                    //     ovc->y_pos_strand, new_xBeg, new_xEnd, new_yBeg, new_yEnd, ovc->x_pos_s, ovc->x_pos_e, ovc->y_pos_s, ovc->y_pos_e);
                    //     fprintf(stderr, "x_readLen: %lld, y_readLen: %lld\n", x_readLen, y_readLen);
                    //     fprintf(stderr, "xID: %lld, yID: %lld\n", ovc->x_id, ovc->y_id);
                    //     fprintf(stderr, "xRegionLen: %lld, yRegionLen: %lld\n", xRegionLen, yRegionLen);
                    //     fprintf(stderr, "xOffset: %lld, yOffset: %lld\n", xOffset, yOffset);
                    //     fprintf(stderr, "yBuoundaryXcoordinate: %lld\n", yBuoundaryXcoordinate);
                    // }
                }
            }
        }
        else
        {
            ///return 0;
            if(i != (long long)(chain->length - 1)) zdrop_occ++;
            if(zdrop_occ > 1) return 0;
            sum_score += mapExtentScore;
            if(xRegionLen <= yRegionLen)
            {
                sum_score -= (GAP_OPEN_KSW + (xRegionLen - maxXpos)*GAP_EXT_KSW);
            }
            else
            {
                sum_score -= (GAP_OPEN_KSW + (yRegionLen - maxYpos)*GAP_EXT_KSW);
            }
        }
    }

    (*final_scores) = sum_score;
    if(new_xBeg != (long long)ovc->x_pos_s || new_xEnd != (long long)ovc->x_pos_e || 
       new_yBeg != (long long)ovc->y_pos_s || new_yEnd != (long long)ovc->y_pos_e)
    {   
        // fprintf(stderr, "\ntttdirection: %u, new_xBeg: %lld, new_xEnd: %lld, new_yBeg: %lld, new_yEnd: %lld, old_xBeg: %u, old_xEnd: %u, old_yBeg: %u, old_yEnd: %u\n", 
        // ovc->y_pos_strand, new_xBeg, new_xEnd, new_yBeg, new_yEnd, ovc->x_pos_s, ovc->x_pos_e, ovc->y_pos_s, ovc->y_pos_e);
        // fprintf(stderr, "x_readLen: %lld, y_readLen: %lld\n", x_readLen, y_readLen);
        // fprintf(stderr, "xID: %lld, yID: %lld\n", ovc->x_id, ovc->y_id);
        // for (i = 0; i < (long long)chain->length; i++)
        // {
        //     fprintf(stderr,"i: %lld, x_pos: %d, offset: %d\n", 
        //     i, get_fake_gap_pos(chain, i), get_fake_gap_shift(chain, i));
        // }
    }
    else
    {
        // fprintf(stderr, "\nkkkdirection: %u, new_xBeg: %lld, new_xEnd: %lld, new_yBeg: %lld, new_yEnd: %lld, old_xBeg: %u, old_xEnd: %u, old_yBeg: %u, old_yEnd: %u\n", 
        // ovc->y_pos_strand, new_xBeg, new_xEnd, new_yBeg, new_yEnd, ovc->x_pos_s, ovc->x_pos_e, ovc->y_pos_s, ovc->y_pos_e);
        // fprintf(stderr, "x_readLen: %lld, y_readLen: %lld\n", x_readLen, y_readLen);
        // fprintf(stderr, "xID: %lld, yID: %lld\n", ovc->x_id, ovc->y_id);
        // for (i = 0; i < (long long)chain->length; i++)
        // {
        //     fprintf(stderr,"i: %lld, x_pos: %d, offset: %d\n", 
        //     i, get_fake_gap_pos(chain, i), get_fake_gap_shift(chain, i));
        // }
    }
    **/
    
    return 1;
}



int fill_chain_by_affine_gap(Fake_Cigar* chain, char* x_string, char* y_string, overlap_region* ovc,
long long x_readLen, long long y_readLen, Cigar_record* cigar, uint8_t* c2n, uint8_t* x_num, uint8_t* y_num,
long long* minus_score_thres, long long* final_scores)
{
    /**
    long long i, xOffset, yOffset, xRegionLen, yRegionLen, maxXpos, maxYpos, zdroped;
    long long mapGlobalScore, mapExtentScore;
    long long xBuoundaryScore, xBuoundaryYcoordinate, yBuoundaryScore, yBuoundaryXcoordinate;
    ///float band_rate = 0.08;
    int endbouns, mode;
    long long xBeg, yBeg;
    xBeg = ovc->x_pos_s;
    yBeg = ovc->y_pos_s;
    if(chain->length <= 0) return 0;
    long long sum_score = 0, current_ovlp = 0, zdrop_occ = 0;
    long long chain_num = (long long)chain->length - 1;
    
    ///align forward
    for (i = 0; i < chain_num; i++)
    {
        xOffset = get_fake_gap_pos(chain, i);
        yOffset = (xOffset - xBeg) + yBeg + get_fake_gap_shift(chain, i);
        ///last region
        
        ///higher endbouns for middle regions
        endbouns = MATCH_SCORE_KSW;
        xRegionLen = get_fake_gap_pos(chain, i+1) - xOffset;
        yRegionLen = (get_fake_gap_pos(chain, i+1) + get_fake_gap_shift(chain, i+1)) - 
                        (get_fake_gap_pos(chain, i) + get_fake_gap_shift(chain, i));
        ///last region
        if(i == chain_num - 1)
        {
            xRegionLen++;
            yRegionLen++;
        }

        mode = KSW_EZ_SCORE_ONLY | KSW_EZ_APPROX_MAX | KSW_EZ_APPROX_DROP;
            

        if(minus_score_thres)
        {
            current_ovlp = MIN((xOffset - xBeg), (yOffset - yBeg));
            current_ovlp = current_ovlp*MATCH_SCORE_KSW;
            if(current_ovlp - sum_score > (*minus_score_thres))
            {
                return 0;
            }
        }

        if(xOffset < 0) xOffset = 0;
        if(yOffset < 0) yOffset = 0;
        if(xRegionLen < 0) xRegionLen = 0;
        if(yRegionLen < 0) yRegionLen = 0;

        ///do alignment forward
        ///text is x, query is y
        afine_gap_alignment(x_string+xOffset, x_num, xRegionLen, y_string+yOffset, y_num, yRegionLen, 
        c2n, FORWARD_KSW, MATCH_SCORE_KSW, MISMATCH_SCORE_KSW, GAP_OPEN_KSW, GAP_EXT_KSW,
        BAND_KSW, Z_DROP_KSW, endbouns, &maxXpos, &maxYpos, &mapGlobalScore, 
        &mapExtentScore, &xBuoundaryScore, &xBuoundaryYcoordinate, 
        &yBuoundaryScore, &yBuoundaryXcoordinate, &zdroped, mode);
        

        if(!zdroped)
        {
            sum_score += mapGlobalScore;
        }
        else
        {
            ///return 0;
            zdrop_occ++;
            ///if(zdrop_occ > 1) return 0;
            sum_score += mapExtentScore;
            if(xRegionLen <= yRegionLen)
            {
                sum_score -= (GAP_OPEN_KSW + (xRegionLen - maxXpos)*GAP_EXT_KSW);
            }
            else
            {
                sum_score -= (GAP_OPEN_KSW + (yRegionLen - maxYpos)*GAP_EXT_KSW);
            }
        }
    }

    (*final_scores) = sum_score;
    **/
    return 1;
}



long long get_affine_gap_score(overlap_region* ovc, UC_Read* g_read, UC_Read* overlap_read, uint8_t* x_num,
uint8_t* y_num, uint64_t EstimateXOlen, uint64_t EstimateYOlen)
{
    char* x_string;
    char* y_string;
    uint64_t yStrand;
    long long minus_score_thres = (MAX(EstimateXOlen, EstimateYOlen)*HIGH_HET_ERROR_RATE*(MATCH_SCORE_KSW+(MAX(MISMATCH_SCORE_KSW,GAP_EXT_KSW))));
    long long total_score_thres = MAX(EstimateXOlen, EstimateYOlen)*MATCH_SCORE_KSW - minus_score_thres;

    yStrand = ovc->y_pos_strand;

    if(yStrand == 0)
    {
        recover_UC_Read(overlap_read, &R_INF, ovc->y_id);
    }
    else
    {
        recover_UC_Read_RC(overlap_read, &R_INF, ovc->y_id);
    }
    x_string = g_read->seq;
    y_string = overlap_read->seq;

    long long sum = 0;
    
    if(fill_chain_by_affine_gap(&(ovc->f_cigar), x_string, y_string, ovc, Get_READ_LENGTH(R_INF, ovc->x_id), 
    Get_READ_LENGTH(R_INF, ovc->y_id), NULL, seq_nt6_table, x_num, y_num, &minus_score_thres, &sum) == 0)
    {
        return 0;
    }

    if(sum >= total_score_thres) return 1;
    return 0;
}

/**
void recalcate_high_het_overlap(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read)
{
    long long j, k, i;
    int threshold;
    long long y_id;
    int y_strand;
    long long y_readLen;
    long long x_start;
    long long x_end;
    long long x_len;
    long long total_y_start;
    long long total_y_end;
    long long y_start;
    long long Window_Len;
    char* x_string;
    char* y_string;
    int end_site;
    unsigned int error;
    int real_y_start;
    long long overlap_length;
    int extra_begin, extra_end;
    long long o_len;
    kvec_t(uint8_t) x_num;
    kvec_t(uint8_t) y_num;
    kv_init(x_num);
    kv_init(y_num);

    for (j = 0; j < (long long)overlap_list->length; j++)
    {
        
        if(overlap_list->list[j].w_list_length == 0) continue;
        y_id = overlap_list->list[j].y_id;
        y_strand = overlap_list->list[j].y_pos_strand;
        y_readLen = Get_READ_LENGTH((*R_INF), y_id);

        //i corresponding to each window of a overlap
        //utilize the the end pos of pre-window in backwards
        for (i = overlap_list->list[j].w_list_length - 1; i >= 0; i--)
        {
            ///the first matched window
            if(overlap_list->list[j].w_list[i].y_end != -1)
            {
                ///note!!! need notification
                ///this is the actual end postion in ystring
                total_y_start = overlap_list->list[j].w_list[i].y_end 
                            - overlap_list->list[j].w_list[i].extra_begin + 1;

                ///k corresponding to all unmatched windows at the right side of overlap_list->list[j].w_list[i]
                ///so k starts from i + 1, and end to the first matched window
                for (k = i + 1; k < (long long)overlap_list->list[j].w_list_length && overlap_list->list[j].w_list[k].y_end == -1; k++)
                {
                    extra_begin = extra_end = 0;

                    ///if y_start > y_readLen, direct terminate
                    if (total_y_start >= y_readLen)
                    {
                        break;
                    }
                    
                    ///there is no problem for x
                    x_start = overlap_list->list[j].w_list[k].x_start;
                    x_end = overlap_list->list[j].w_list[k].x_end;
                    x_len = x_end - x_start + 1;
                    ///there are two potiential reasons for unmatched window:
                    ///1. this window has a large number of differences
                    ///2. DP does not start from the right offset
                    threshold = double_error_threshold(overlap_list->list[j].w_list[k].error_threshold, x_len);
                    
                    y_start = total_y_start;
                    Window_Len = x_len + (threshold << 1);

                    if(!determine_overlap_region(threshold, y_start, y_id, Window_Len, Get_READ_LENGTH((*R_INF), y_id), 
                    &extra_begin, &extra_end, &y_start, &o_len))
                    {
                        break;
                    }

                    if(o_len + threshold < x_len)
                    {
                        break;
                    }
                    
                    fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, 
                    R_INF, y_id, extra_begin, extra_end);

                    x_string = g_read->seq + x_start;
                    y_string = dumy->overlap_region;

                    ///note!!! need notification
                    end_site = Reserve_Banded_BPM(y_string, Window_Len, x_string, x_len, threshold, &error);
                    
                    ///if error==-1, unmatched
                    if (error!=(unsigned int)-1)
                    {
                        overlap_list->list[j].w_list[k].cigar.length = -1;
                        overlap_list->list[j].w_list[k].y_start = y_start;
                        overlap_list->list[j].w_list[k].y_end = y_start + end_site;
                        overlap_list->list[j].w_list[k].error = (int)error;
                        ///note!!! need notification
                        overlap_list->list[j].w_list[k].extra_begin = extra_begin;
                        overlap_list->list[j].w_list[k].extra_end = extra_end;
                        overlap_list->list[j].w_list[k].error_threshold = threshold;
                        
                        overlap_list->list[j].align_length += x_len;
                    }
                    else
                    {
                        break;
                    }

                    ///note!!! need notification
                    total_y_start = y_start + end_site - extra_begin + 1;
                }
                
            }
            
        }
        


        //i corresponding to each window of a overlap
        //utilize the the start pos of next window in forward
        for (i = 0; i < (long long)overlap_list->list[j].w_list_length; i++)
        {
            ///find the first matched window, which should not be the first window
            ///the pre-window of this matched window must be unmatched
            if(overlap_list->list[j].w_list[i].y_end != -1 && i != 0 && overlap_list->list[j].w_list[i - 1].y_end == -1)
            {
                ///check if the start pos of this matched window has been calculated
                if(overlap_list->list[j].w_list[i].cigar.length == -1)
                {
                    ///there is no problem for x
                    x_start = overlap_list->list[j].w_list[i].x_start;
                    x_end = overlap_list->list[j].w_list[i].x_end;
                    x_len = x_end - x_start + 1;
                    //may have bugs
                    threshold = overlap_list->list[j].w_list[i].error_threshold;
                    //may have bugs
                    //may have bugs
                    ///should not adjust threshold, since this window can be matched by the old threshold
                    ///threshold = Adjust_Threshold(threshold, x_len);
                    //may have bugs
                    Window_Len = x_len + (threshold << 1);


                    ///y_start is the real y_start
                    y_start = overlap_list->list[j].w_list[i].y_start;
                    extra_begin = overlap_list->list[j].w_list[i].extra_begin;
                    extra_end = overlap_list->list[j].w_list[i].extra_end;
                    o_len = Window_Len - extra_end - extra_begin;
                    fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, 
                    R_INF, y_id, extra_begin, extra_end);
                    x_string = g_read->seq + x_start;
                    y_string = dumy->overlap_region;

                    ///note!!! need notification
                    end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
                    &(dumy->path_length), dumy->matrix_bit, dumy->path, 
                    overlap_list->list[j].w_list[i].error, overlap_list->list[j].w_list[i].y_end - y_start);


                    ///y_start has already been calculated
                    if (error != (unsigned int)-1)
                    {
                        ///this condition is always wrong
                        ///in best case, real_y_start = threshold, end_site = Window_Len - threshold - 1
                        if (end_site == Window_Len - 1 || real_y_start == 0)
                        {
                            if(fix_boundary(x_string, x_len, threshold, y_start, real_y_start, 
                            end_site, extra_begin, extra_end, y_id, Window_Len, R_INF, dumy, 
                            y_strand, error, &y_start, &real_y_start, &end_site, &extra_begin, 
                            &extra_end, &error))
                            {
                                overlap_list->list[j].w_list[i].error = error;
                                overlap_list->list[j].w_list[i].extra_begin = extra_begin;
                                overlap_list->list[j].w_list[i].extra_end = extra_end;
                            }
                        }
                                                 
                        generate_cigar(dumy->path, dumy->path_length, &(overlap_list->list[j].w_list[i]),
                        &real_y_start, &end_site, &error, x_string, x_len, y_string);   

                        ///note!!! need notification
                        real_y_start = y_start + real_y_start - extra_begin;
                        overlap_list->list[j].w_list[i].y_start = real_y_start; 
                        ///I forget why don't reduce the extra_begin for y_end
                        ///it seems extra_begin will be reduced at the end of this function 
                        overlap_list->list[j].w_list[i].y_end = y_start + end_site;
                        overlap_list->list[j].w_list[i].error = error;                         
                    }
                    else
                    {
                        fprintf(stderr, "error\n");
                    }
                }
                else
                {
                    real_y_start = overlap_list->list[j].w_list[i].y_start;
                }

                
                ///the end pos for pre window is real_y_start - 1
                total_y_end = real_y_start - 1;
                ///find the unmatched window on the left of current matched window
                ///k starts from i - 1
                for (k = i - 1; k >= 0 && overlap_list->list[j].w_list[k].y_end == -1; k--)
                {  
                    ///there is no problem in x
                    x_start = overlap_list->list[j].w_list[k].x_start;
                    x_end = overlap_list->list[j].w_list[k].x_end;
                    x_len = x_end - x_start + 1;
                    ///there are two potiential reasons for unmatched window:
                    ///1. this window has a large number of differences
                    ///2. DP does not start from the right offset
                    threshold = double_error_threshold(overlap_list->list[j].w_list[k].error_threshold, x_len);

                    Window_Len = x_len + (threshold << 1);

                    if(total_y_end <= 0)
                    {
                        break;
                    }

                    ///y_start might be less than 0
                    y_start = total_y_end - x_len + 1;
                    if(!determine_overlap_region(threshold, y_start, y_id, Window_Len, Get_READ_LENGTH((*R_INF), y_id),
                    &extra_begin, &extra_end, &y_start, &o_len))
                    {
                        break;
                    }

                    if(o_len + threshold < x_len)
                    {
                        break;
                    }

                    fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, 
                    R_INF, y_id, extra_begin, extra_end);
                    x_string = g_read->seq + x_start;
                    y_string = dumy->overlap_region;

                    ///note!!! need notification
                    end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
                    &(dumy->path_length), dumy->matrix_bit, dumy->path, -1, -1);

                    if (error!=(unsigned int)-1)
                    { 
                        ///this condition is always wrong
                        ///in best case, real_y_start = threshold, end_site = Window_Len - threshold - 1
                        if (end_site == Window_Len - 1 || real_y_start == 0)
                        {
                            fix_boundary(x_string, x_len, threshold, y_start, real_y_start, end_site,
                            extra_begin, extra_end, y_id, Window_Len, R_INF, dumy, y_strand, error,
                            &y_start, &real_y_start, &end_site,
                            &extra_begin, &extra_end, &error);
                        }

                        generate_cigar(dumy->path, dumy->path_length, &(overlap_list->list[j].w_list[k]),
                        &real_y_start, &end_site, &error, x_string, x_len, y_string);  

                        ///y_start has no shift, but y_end has shift           
                        overlap_list->list[j].w_list[k].y_start = y_start + real_y_start - extra_begin;
                        overlap_list->list[j].w_list[k].y_end = y_start + end_site;
                        overlap_list->list[j].w_list[k].error = error;
                        overlap_list->list[j].align_length += x_len;
                        overlap_list->list[j].w_list[k].extra_begin = extra_begin;
                        overlap_list->list[j].w_list[k].extra_end = extra_end;
                        overlap_list->list[j].w_list[k].error_threshold = threshold;
                    }
                    else
                    {
                        break;
                    }

                    total_y_end = y_start + real_y_start - 1 - extra_begin;
                }
            }
        }
    }



    overlap_list->mapped_overlaps_length = 0;
    
    double error_rate;
    int is_update = 0;
    for (j = 0; j < (long long)overlap_list->length; j++)
    {
        y_id = overlap_list->list[j].y_id;
        y_strand = overlap_list->list[j].y_pos_strand;
        y_readLen = Get_READ_LENGTH((*R_INF), y_id);
        overlap_length = overlap_list->list[j].x_pos_e - overlap_list->list[j].x_pos_s + 1;
        overlap_list->list[j].is_match = 0;
        is_update = 0;
        ///debug_scan_cigar(&(overlap_list->list[j]));
        if(overlap_list->list[j].w_list_length == 0 || overlap_length == 0 || overlap_list->list[j].align_length == 0) continue;
        ///only calculate cigar for high quality overlaps
        if (overlap_length * OVERLAP_THRESHOLD_FILTER <=  overlap_list->list[j].align_length)
        {
            
            for (i = 0; i < (long long)overlap_list->list[j].w_list_length; i++)
            {
                ///first we need to check if this window is matched
                if(overlap_list->list[j].w_list[i].y_end != -1)
                {
                    ///second check if the cigar of this window has been got 
                    if(overlap_list->list[j].w_list[i].cigar.length == -1)
                    {
                        ///there is no problem for x
                        x_start = overlap_list->list[j].w_list[i].x_start;
                        x_end = overlap_list->list[j].w_list[i].x_end;
                        x_len = x_end - x_start + 1;
                        //may have bugs
                        ///threshold = x_len * asm_opt.max_ov_diff_ec;
                        threshold = overlap_list->list[j].w_list[i].error_threshold;
                        //may have bugs
                        //may have bugs
                        ///should not adjust threshold, since this window can be matched by the old threshold
                        ///threshold = Adjust_Threshold(threshold, x_len);
                        //may have bugs
                        Window_Len = x_len + (threshold << 1);


                        ///y_start is the real y_start
                        ///for the window with cigar, y_start has already reduced extra_begin
                        y_start = overlap_list->list[j].w_list[i].y_start;
                        extra_begin = overlap_list->list[j].w_list[i].extra_begin;
                        extra_end = overlap_list->list[j].w_list[i].extra_end;
                        o_len = Window_Len - extra_end - extra_begin;
                        fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, 
                        R_INF, y_id, extra_begin, extra_end);
                        x_string = g_read->seq + x_start;
                        y_string = dumy->overlap_region;


                        ///note!!! need notification
                        end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
                        &(dumy->path_length), dumy->matrix_bit, dumy->path, 
                        overlap_list->list[j].w_list[i].error, overlap_list->list[j].w_list[i].y_end - y_start);

                        if (error != (unsigned int)-1)
                        {
                            if (end_site == Window_Len - 1 || real_y_start == 0)
                            {
                                
                                if(fix_boundary(x_string, x_len, threshold, y_start, real_y_start, end_site,
                                extra_begin, extra_end, y_id, Window_Len, R_INF, dumy, y_strand, error,
                                &y_start, &real_y_start, &end_site,
                                &extra_begin, &extra_end, &error))
                                {
                                    overlap_list->list[j].w_list[i].error = error;
                                    overlap_list->list[j].w_list[i].extra_begin = extra_begin;
                                    overlap_list->list[j].w_list[i].extra_end = extra_end;
                                }
                                
                            }

                            generate_cigar(dumy->path, dumy->path_length, &(overlap_list->list[j].w_list[i]),
                            &real_y_start, &end_site, &error, x_string, x_len, y_string);    

                            ///note!!! need notification
                            real_y_start = y_start + real_y_start - extra_begin;
                            overlap_list->list[j].w_list[i].y_start = real_y_start;  
                            overlap_list->list[j].w_list[i].y_end = y_start + end_site - extra_begin;
                            overlap_list->list[j].w_list[i].error = error;                              
                        }
                        else
                        {
                            fprintf(stderr, "error\n");
                        }
                    }
                    else
                    {
                        overlap_list->list[j].w_list[i].y_end -= overlap_list->list[j].w_list[i].extra_begin;
                    }


                }
            }

            error_rate = non_trim_error_rate(overlap_list, j, R_INF, dumy, g_read);
            
            
            if (error_rate <= HIGH_HET_ERROR_RATE)
            {
                is_update = 1;
            }
        }

        if((is_update == 0) && (overlap_list->list[j].align_length >= WINDOW) && 
           (overlap_length * HIGH_HET_OVERLAP_THRESHOLD_FILTER <=  overlap_list->list[j].align_length))
        {
            kv_resize(uint8_t, x_num, (uint64_t)(Get_READ_LENGTH((*R_INF), overlap_list->list[j].x_id)));
            kv_resize(uint8_t, y_num, (uint64_t)(Get_READ_LENGTH((*R_INF), overlap_list->list[j].y_id)));
            is_update = get_affine_gap_score(&(overlap_list->list[j]), g_read, overlap_read, x_num.a, y_num.a,
            overlap_list->list[j].x_pos_e + 1 - overlap_list->list[j].x_pos_s, 
            overlap_list->list[j].y_pos_e + 1 - overlap_list->list[j].y_pos_s);
        }

        if(is_update)
        {
            overlap_list->list[j].is_match = 2;
        }
    }

    kv_destroy(x_num);
    kv_destroy(y_num);
}



void correct_overlap_high_het(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read)
{
    clear_Correct_dumy(dumy, overlap_list, NULL);

    long long window_start, window_end;

    Window_Pool w_inf;

    init_Window_Pool(&w_inf, g_read->length, WINDOW, (int)(1.0/asm_opt.max_ov_diff_ec));

    int flag = 0;

    while(get_Window(&w_inf, &window_start, &window_end) && flag != -2)
    {
        dumy->length = 0;
        dumy->lengthNT = 0;
        flag = get_interval(window_start, window_end, overlap_list, dumy, w_inf.window_length);

        switch (flag)
        {
            case 1:    ///no match here
                break;
            case 0:    ///no match here
                break;
            case -2: ///if flag == -2, loop would be terminated
                break;
        }

        ///dumy->lengthNT represent how many overlaps that the length of them is not equal to WINDOW; may larger or less than WINDOW
        ///dumy->length represent how many overlaps that the length of them is WINDOW
        //may improve
        ///now the windows which are larger than WINDOW are verified one-by-one, to improve it, we can do it group-bygroup
        verify_window(window_start, window_end, overlap_list, dumy, R_INF, g_read->seq);
    }

    recalcate_high_het_overlap(overlap_list, R_INF, g_read, dumy, overlap_read);
}
**/

uint64_t update_ov_track_0(Fake_Cigar* z, overlap_region *o, int64_t apend_be, int64_t xl, int64_t yl,
k_mer_hit* hit, int64_t n_hit)
{
    int64_t k, dq, dr, dd, pdd = INT32_MAX, xr, yr; 
    z->length = 0;
    if(hit[0].readID != o->y_id || hit[0].strand != o->y_pos_strand) return 0;
    ///update o->s
    o->x_pos_s = hit[0].self_offset; o->y_pos_s = hit[0].offset;
    if(o->x_pos_s <= o->y_pos_s) {
        o->y_pos_s -= o->x_pos_s; o->x_pos_s = 0;
    } else {
        o->x_pos_s -= o->y_pos_s; o->y_pos_s = 0;
    }

    if(apend_be == 1) add_fake_cigar(z, o->x_pos_s, 0, NULL);
    for (k = 0; (k < n_hit) && (hit[k].readID == o->y_id) && (hit[k].strand == o->y_pos_strand); k++) {
        dq = hit[k].self_offset - o->x_pos_s; 
        dr = hit[k].offset - o->y_pos_s; 
        dd = dr - dq;
        if(dd != pdd) {
            pdd = dd;
            add_fake_cigar(z, hit[k].self_offset, pdd, NULL);
        }
    }
    ///update o->s
    o->x_pos_e = hit[k-1].self_offset; o->y_pos_e = hit[k-1].offset;
    xr = xl-o->x_pos_e-1; yr = yl-o->y_pos_e-1;
    if(xr <= yr) {
        o->x_pos_e = xl-1; o->y_pos_e += xr;        
    } else {
        o->y_pos_e = yl-1; o->x_pos_e += yr; 
    }

    if((apend_be == 1) && (get_fake_gap_pos(z, z->length-1)!=((int64_t)o->x_pos_e))) {
        add_fake_cigar(z, o->x_pos_e, get_fake_gap_shift(z, z->length-1), NULL);
    }

    return k;
}

int64_t iter_hpc(uint8_t *m, int64_t mn, int64_t *mo, int64_t rev, int64_t *so, int64_t *ho, int64_t sc)
{
    if(sc == 0) return 0;
    if(sc <= (*so)) return (*ho)-1;
    if(!rev) {
        while((*mo) < mn) {
            for (; (*mo) < mn && m[(*mo)] == 255; (*mo)++) {
                (*so) += m[(*mo)];
            }
            (*so) += m[(*mo)]; (*mo)++; (*ho)++;
            if(sc <= (*so)) return (*ho)-1;
        }        
    } else {
        while((*mo) < mn) {
            (*so) += m[mn-(*mo)-1]; (*mo)++; (*ho)++;
            for (; (*mo) < mn && m[mn-(*mo)-1] == 255; (*mo)++) {
                (*so) += m[mn-(*mo)-1];
            }
            if(sc <= (*so)) return (*ho)-1;
        }
    }

    return -1;
}

uint64_t update_ov_track_hpc_0(Fake_Cigar* z, overlap_region *o, int64_t apend_be, int64_t xhl,
uint32_t *x_idx, int64_t y_idx_map_l, uint8_t *y_idx_map, hpc_t *hpc_g, k_mer_hit* hit, int64_t n_hit)
{
    int64_t k, dq, dr, dd, pdd = INT32_MAX, xr, yr, yhl, mo = 0, so = 0, ho = 0, x1, y1; 
    z->length = 0;
    if(hit[0].readID != o->y_id || hit[0].strand != o->y_pos_strand) return 0;
    yhl = hpc_len(*hpc_g, o->y_id);

    x1 = x_idx[hit[0].self_offset];
    y1 = iter_hpc(y_idx_map, y_idx_map_l, &mo, o->y_pos_strand, &so, &ho, hit[0].offset); 
    assert(y1 >= 0);
    o->x_pos_s = x1; o->y_pos_s = y1;
    
    ///update o->s
    if(o->x_pos_s <= o->y_pos_s) {
        o->y_pos_s -= o->x_pos_s; o->x_pos_s = 0;
    } else {
        o->x_pos_s -= o->y_pos_s; o->y_pos_s = 0;
    }
    

    if(apend_be == 1) add_fake_cigar(z, o->x_pos_s, 0, NULL);
    for (k = 0; (k < n_hit) && (hit[k].readID == o->y_id) && (hit[k].strand == o->y_pos_strand); k++) {
        x1 = x_idx[hit[k].self_offset];
        y1 = iter_hpc(y_idx_map, y_idx_map_l, &mo, o->y_pos_strand, &so, &ho, hit[k].offset); 
        assert(y1 >= 0);
        dq = x1 - o->x_pos_s; 
        dr = y1 - o->y_pos_s; 
        dd = dr - dq;
        if(dd != pdd) {
            pdd = dd;
            add_fake_cigar(z, x1, pdd, NULL);
        }
    }
    ///update o->s
    x1 = x_idx[hit[k-1].self_offset];
    y1 = iter_hpc(y_idx_map, y_idx_map_l, &mo, o->y_pos_strand, &so, &ho, hit[k-1].offset); 
    assert(y1 >= 0);
    o->x_pos_e = x1; o->y_pos_e = y1;
    xr = xhl-o->x_pos_e-1; yr = yhl-o->y_pos_e-1;
    if(xr <= yr) {
        o->x_pos_e = xhl-1; o->y_pos_e += xr;        
    } else {
        o->y_pos_e = yhl-1; o->x_pos_e += yr; 
    }

    if((apend_be == 1) && (get_fake_gap_pos(z, z->length-1)!=((int64_t)o->x_pos_e))) {
        add_fake_cigar(z, o->x_pos_e, get_fake_gap_shift(z, z->length-1), NULL);
    }

    return k;
}

///(char *qstr, kvec_t_u64_warp* q_idx) -> only used for hpc
uint64_t update_ol_track(overlap_region_alloc* ol, Candidates_list *cl, hpc_t *hpc_g, const ul_idx_t *udb, uint32_t apend_be, uint64_t qlen, 
char *qstr, kvec_t_u32_warp* q_idx)
{
    uint64_t cln = cl->length, i, k, l, m = 0; overlap_region *r; 
    if(hpc_g) {
        if(ol->length) {
            q_idx->a.n = 0; kv_resize(uint32_t, q_idx->a, qlen); m = 0;
            for (l = 0, k = 1; k <= qlen; k++) {
                if((k == qlen) || (qstr[k] != qstr[l]) || (seq_nt4_table[(uint8_t)qstr[l]] >= 4)) {
                    for (i = l; i < k; i++) q_idx->a.a[i] = m;
                    l = k; m++;
                }
            }
    
            for (i = k = 0; i < ol->length; ++i) {
                r = &(ol->list[i]); k = r->non_homopolymer_errors;
                update_ov_track_hpc_0(&(r->f_cigar), r, apend_be, m, q_idx->a.a, 
                (uint32_t)hpc_g->mm->idx[r->y_id], hpc_g->mm->a + (hpc_g->mm->idx[r->y_id]>>32), 
                hpc_g, cl->list+k, cln-k);
            }
        }
    } else {
        if(ol->length) {
            for (i = k = 0; i < ol->length; ++i) {
                r = &(ol->list[i]); k = r->non_homopolymer_errors;
                update_ov_track_0(&(r->f_cigar), r, apend_be, qlen, udb->ug->u.a[r->y_id].len, cl->list+k, cln-k);
            }
        }
    }
    return m;///hpc length
}

uint64_t gen_hpc_str(const char *in, uint32_t in_l, UC_Read *z, uint64_t *in_hl)
{
    uint64_t hl, k, l;
    if(in_hl) {
        hl = (*in_hl);
    } else {
        for (l = hl = 0, k = 1; k <= in_l; k++) {
            if((k == in_l) || (in[k] != in[l]) || (seq_nt4_table[(uint8_t)in[l]] >= 4)) {
                hl++; l = k;
            }
        }
    }
    resize_UC_Read(z, hl); z->length = 0;
    for (l = 0, k = 1; k <= in_l; k++) {
        if((k == in_l) || (in[k] != in[l]) || (seq_nt4_table[(uint8_t)in[l]] >= 4)) {
            z->seq[z->length++] = in[l]; l = k;
        }
    }
    return hl;
}

///ts do not have aux_beg, while te has 
uint32_t push_wlst(const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, overlap_region* ol, 
                        char* qstr, char *tstr, char *tstr_1, Correct_dumy* dumy,
                        int64_t qs, int64_t qe, int64_t ts, int64_t te, int64_t tl, 
                        int64_t error, int64_t aux_beg, int64_t aux_end, int64_t thres, double e_rate, 
                        int64_t block_s, uint32_t sec_check, double ovlp_cut, void *km)
{

    window_list p, t, *a; int64_t w_e, w_s, ce = qs - 1, cs = ol->x_pos_s, toff, ovl, ualn, aln; 
    uint64_t a_n, k;
    
    p.x_start = qs; p.x_end = qe; p.y_start = ts; p.y_end = te; p.error = error;
    p.extra_begin = aux_beg; p.extra_end = aux_end;
    p.error_threshold = thres; p.cidx = p.clen = 0;
    if(ol->w_list.n > 0) { //utilize the the end pos of pre-window in forward
        w_e = ol->w_list.a[ol->w_list.n-1].x_end; 
        toff = ol->w_list.a[ol->w_list.n-1].y_end + 1 - ol->w_list.a[ol->w_list.n-1].extra_begin;
        while ((w_e < ce) && (toff < tl)) {
            w_s = w_e + 1;
            get_win_id_by_s(ol, w_s, block_s, &w_e);
            // x_start = w_s; x_end = w_e;
            if(aln_wlst_adv(ol, rref, hpc_g, uref, qstr, tstr, tstr_1, dumy,
            ol->y_pos_strand, ol->y_id, w_s, w_e, toff, block_s, e_rate, 0)) {
                toff = ol->w_list.a[ol->w_list.n-1].y_end + 1 - ol->w_list.a[ol->w_list.n-1].extra_begin;
            } else {
                break;
            }
        }
        cs = ol->w_list.a[ol->w_list.n-1].x_end + 1;
    }
    ///utilize the the start pos of next window in backward
    a_n = ol->w_list.n; w_s = qs;
    if(w_s > cs) {
        gen_backtrace_adv(&p, ol, rref, hpc_g, uref, qstr, tstr, tstr_1, dumy, ol->y_pos_strand, ol->y_id);
        p.y_end += p.extra_begin;
        toff = p.y_start - 1;
        while ((w_s > cs) && (toff > 0)) {
            w_e = w_s - 1;
            get_win_id_by_e(ol, w_e, block_s, &w_s);
            // x_start = w_s; x_end = w_e; x_len = x_end + 1 - x_start;
            if(aln_wlst_adv(ol, rref, hpc_g, uref, qstr, tstr, tstr_1, dumy,
                ol->y_pos_strand, ol->y_id, w_s, w_e, toff+1-(w_e+1-w_s), block_s, e_rate, 1)) {
                ///y_start has no shift, but y_end has shift  
                ol->w_list.a[ol->w_list.n-1].y_start -= ol->w_list.a[ol->w_list.n-1].extra_begin; 
                toff = ol->w_list.a[ol->w_list.n-1].y_start - 1;
            } else {
                break;
            }
        }
    }

    ol->align_length += qe + 1 - qs;
    ovl = ol->x_pos_e+1-ol->x_pos_s; ualn = (qe + 1 - ol->x_pos_s) - ol->align_length; aln = ovl-ualn;
    if((!simi_pass(ovl, aln, 0, ovlp_cut, &e_rate)) && (!simi_pass(ovl, aln, sec_check, ovlp_cut, NULL))) {
        kv_push(window_list, ol->w_list, p); 
        return 0;
    }

    if(ol->w_list.n > a_n) {
        a = ol->w_list.a + a_n; a_n = ol->w_list.n - a_n; toff = a_n; a_n >>=1;
        for (k = 0; k < a_n; k++) {
            t = a[k]; a[k] = a[toff-1-k]; a[toff-1-k] = t;
        }
    }
    kv_push(window_list, ol->w_list, p); 
    return 1;
}

uint32_t align_ul_ed_post(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, char* qstr, char *tstr, char *tstr_1, 
Correct_dumy* dumy, double e_rate, int64_t w_l, double ovlp_cut, void *km)
{
    int64_t q_s, q_e, nw, k, q_l, t_tot_l, sec_check = (uref&&(!hpc_g))?1:0;
    int64_t aux_beg, aux_end, t_s, thre, aln_l, t_pri_l, t_end;
    char *q_string, *t_string; unsigned int error;
    z->w_list.n = 0; z->is_match = 0; z->align_length = 0;
    nw = get_num_wins(z->x_pos_s, z->x_pos_e+1, w_l);
    get_win_se_by_normalize_xs(z, (z->x_pos_s/w_l)*w_l, w_l, &q_s, &q_e);
    for (k = 0; k < nw; k++) {
        aux_beg = aux_end = 0; q_l = 1 + q_e - q_s;
        thre = q_l*e_rate; thre = Adjust_Threshold(thre, q_l);
        if(thre > THRESHOLD_MAX_SIZE) thre = THRESHOLD_MAX_SIZE;
        ///offset of y
        t_s = (q_s - z->x_pos_s) + z->y_pos_s;
        t_s += y_start_offset(q_s, &(z->f_cigar));

        aln_l = q_l + (thre<<1); t_tot_l = hpc_g?hpc_len(*hpc_g, z->y_id):uref->ug->u.a[z->y_id].len;
        if(init_waln(thre, t_s, t_tot_l, aln_l, &aux_beg, &aux_end, &t_s, &t_pri_l)) {
            q_string = qstr+q_s; 
            t_string = return_str_seq(tstr, t_s, t_pri_l, z->y_pos_strand, hpc_g, uref, z->y_id, aux_beg, aux_end);

            t_end = Reserve_Banded_BPM(t_string, aln_l, q_string, q_l, thre, &error);
            // int32_t debug_t_end, debug_error;
            // debug_t_end = ed_band_cal_semi(t_string, aln_l, q_string, q_l, thre, &debug_error);
            // if((t_end != debug_t_end) || (t_end >= 0 && debug_t_end >= 0 && debug_error != (int32_t)error)) {
            //     fprintf(stderr, "[M::%s] debug_error->%d, error->%d\n", __func__, debug_error, error);
            // }
            

            if (error!=((unsigned int)-1)) {
                ///t_s do not have aux_beg, while t_s + t_end (aka, te) has 
                if(!push_wlst(uref, hpc_g, NULL, z, qstr, tstr, tstr_1, dumy, q_s, q_e, t_s, t_s + t_end, 
                                        t_tot_l, error, aux_beg, aux_end, thre, e_rate, w_l, sec_check, ovlp_cut, km)) {
                    return 0;
                }
                // append_window_list(z, q_s, q_e, t_s, t_s + t_end, error, aux_beg, aux_end, thre, w_l, km);
            }
        }
        q_s = q_e + 1; q_e = q_s + w_l - 1; 
        if(q_e >= (int64_t)z->x_pos_e) q_e = z->x_pos_e;
    }

    if((!simi_pass(z->x_pos_e+1-z->x_pos_s, z->align_length, 0, ovlp_cut, &e_rate)) &&
        (!simi_pass(z->x_pos_e+1-z->x_pos_s, z->align_length, sec_check, ovlp_cut, NULL))) return 0;
    return 1;
}

void prt_cigar(uint16_t *ca, uint32_t cn)
{
    uint32_t k;
    for (k = 0; k < cn; k++) {
        fprintf(stderr, "%u%c", ca[k]&0x3fff, "EMDI"[ca[k]>>14]);
    }
    fprintf(stderr, "\n");
}


inline uint32_t gen_backtrace_adv_exz(window_list *p, overlap_region *z, All_reads *rref, hpc_t *hpc_g, const ul_idx_t *uref, 
char *qstr, char *tstr, bit_extz_t *exz, uint32_t rev, uint32_t id)
{
    if(p->error < 0 || p->y_end < 0) return 0;
    int64_t qs, qe, ql, tl, aln_l, t_pri_l, thres, ts, t_tot_l; 
    int64_t aux_beg, aux_end;
    char *q_string, *t_string; 
    ///there is no problem for x
    qs = p->x_start; qe = p->x_end; ql = qe + 1 - qs;
    thres = p->error_threshold; aln_l = ql + (thres<<1);

    ///y_start is the real y_start
    ///for the window with cigar, y_start has already reduced extra_begin
    ts = p->y_start; aux_beg = p->extra_begin; aux_end = p->extra_end;
    if(aux_end >= 0) {
        t_pri_l = aln_l - aux_beg - aux_end;
    } else {
        if(hpc_g) t_tot_l = hpc_len(*hpc_g, id);
        else if(uref) t_tot_l = uref->ug->u.a[id].len;
        else t_tot_l = Get_READ_LENGTH((*rref), id);

        t_pri_l = ts + aln_l - aux_beg; if(t_pri_l > t_tot_l) t_pri_l = t_tot_l;
        t_pri_l = t_pri_l - ts;
    }
    
    q_string = qstr + qs; tl = t_pri_l;
    if(rref) {
        recover_UC_Read_sub_region(tstr, ts, t_pri_l, rev, rref, id); t_string = tstr;
    } else {
        t_string = return_str_seq_exz(tstr, ts, t_pri_l, rev, hpc_g, uref, id);
    }

    exz->ts = 0; exz->te = p->x_end-p->x_start; exz->tl = ql;
    exz->ps = -1; exz->pe = p->y_end-p->y_start; exz->pl = tl;  
    exz->err = p->error; exz->thre = p->error_threshold; 
    // clear_align(*exz);
    ed_band_cal_semi_64_w_absent_diag_trace(t_string, tl, q_string, ql, thres, aux_beg, exz);
    // if(id == 178 && p->x_start == 86800 && p->x_end == 86807) {
    //     fprintf(stderr, "\n[M::%s::semi] exz->ps::%d, exz->pe::%d, exz->ts::%d, exz->te::%d, exz->err::%d, exz->cigar.n::%d\n", 
    //     __func__, exz->ps, exz->pe, exz->ts, exz->te, exz->err, (int32_t)exz->cigar.n);
    //     fprintf(stderr, "[M::%s::semi] p->y_start::%d, p->y_end::%d, p->x_start::%d, p->x_end::%d, p->error::%d\n", 
    //     __func__, p->y_start, p->y_end, p->x_start, p->x_end, p->error);
    //     if(is_align(*exz)) {
    //         prt_cigar(exz->cigar.a, exz->cigar.n);
    //         fprintf(stderr, "[tstr] %.*s\n", exz->pe+1-exz->ps, t_string+exz->ps);
    //         fprintf(stderr, "[qstr] %.*s\n", exz->te+1-exz->ts, q_string+exz->ts);
    //     }
    // }
    // assert(is_align(*exz));
    // assert(cigar_check(t_string, q_string, exz));
    
            
    if(is_align(*exz)) {
        p->y_start = ts + exz->ps;///difference
        p->y_end = ts + exz->pe;
        p->error = exz->err;
        push_wcigar(p, &(z->w_list), exz);
        ///this condition is always wrong
        ///in best case, r_ts = threshold, t_end = aln_l - thres - 1
        if ((((exz->pe+1) == tl) || (exz->ps == 0)) && (exz->err > 0)) {
            if(recal_boundary_exz(q_string, tstr, ql, tl, thres, ts, exz->ps, exz->pe,
            exz->err, id, rev, exz, rref, hpc_g, uref, &ts, &aux_beg, &aux_end)) {
                //update cigar
                z->w_list.c.n = p->cidx; push_wcigar(p, &(z->w_list), exz);

                p->y_start = ts + exz->ps;///difference
                p->y_end = ts + exz->pe;
                p->error = exz->err;
            }
        }
        p->extra_begin = aux_beg; 
        p->extra_end = aux_end;
        return 1;
    }
    p->error = -1;
    return 0;
}


///ts do not have aux_beg, while te has 
uint32_t push_wlst_exz(const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, overlap_region* ol, 
                        char* qstr, char *tstr, bit_extz_t *exz,
                        int64_t qs, int64_t qe, int64_t ts, int64_t te, int64_t tl, 
                        int64_t aux_beg, int64_t aux_end, double e_rate, int64_t block_s, uint32_t sec_check, double ovlp_cut, int64_t force_aln, void *km)
{

    window_list p, t, *a; int64_t w_e, w_s, ce = qs - 1, cs = ol->x_pos_s, toff, ovl, ualn, aln, ys; 
    uint64_t a_n, k;
    
    p.x_start = qs; p.x_end = qe; p.y_start = ts; p.y_end = te; p.error = exz->err;
    p.extra_begin = aux_beg; p.extra_end = aux_end; p.error_threshold = exz->thre; p.cidx = p.clen = 0;
    if(ol->w_list.n > 0) { //utilize the the end pos of pre-window in forward
        w_e = ol->w_list.a[ol->w_list.n-1].x_end; 
        toff = ol->w_list.a[ol->w_list.n-1].y_end + 1;
        while ((w_e < ce) && (toff < tl)) {
            w_s = w_e + 1;
            get_win_id_by_s(ol, w_s, block_s, &w_e);
            // x_start = w_s; x_end = w_e;
            if(aln_wlst_adv_exz(ol, rref, hpc_g, uref, qstr, tstr, exz,
            ol->y_pos_strand, ol->y_id, w_s, w_e, toff, block_s, e_rate, 0)) {
                toff = ol->w_list.a[ol->w_list.n-1].y_end + 1;
            } else {
                break;
            }
        }
        cs = ol->w_list.a[ol->w_list.n-1].x_end + 1;
    }
    ///utilize the the start pos of next window in backward
    a_n = ol->w_list.n; w_s = qs;
    if(w_s > cs) {
        gen_backtrace_adv_exz(&p, ol, rref, hpc_g, uref, qstr, tstr, exz, ol->y_pos_strand, ol->y_id);
        toff = p.y_start - 1;
        while (w_s > cs) {
            w_e = w_s - 1;
            get_win_id_by_e(ol, w_e, block_s, &w_s); ys = toff+1-(w_e+1-w_s);
            // x_start = w_s; x_end = w_e; x_len = x_end + 1 - x_start;
            if((ys >= 0) && aln_wlst_adv_exz(ol, rref, hpc_g, uref, qstr, tstr, exz,
                ol->y_pos_strand, ol->y_id, w_s, w_e, ys, block_s, e_rate, 1)) {
                toff = ol->w_list.a[ol->w_list.n-1].y_start - 1;
            } else {
                break;
            }
        }
    }

    ol->align_length += qe + 1 - qs;
    ovl = ol->x_pos_e+1-ol->x_pos_s; ualn = (qe + 1 - ol->x_pos_s) - ol->align_length; aln = ovl-ualn;
    if((!force_aln) && (!simi_pass(ovl, aln, 0, ovlp_cut, &e_rate)) && (!simi_pass(ovl, aln, sec_check, ovlp_cut, NULL))) {
        kv_push(window_list, ol->w_list, p); 
        return 0;
    }

    if(ol->w_list.n > a_n) {
        a = ol->w_list.a + a_n; a_n = ol->w_list.n - a_n; toff = a_n; a_n >>=1;
        for (k = 0; k < a_n; k++) {
            t = a[k]; a[k] = a[toff-1-k]; a[toff-1-k] = t;
        }
    }
    kv_push(window_list, ol->w_list, p); 
    return 1;
}


uint32_t align_ul_ed_post_extz(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, char* qstr, char *tstr, bit_extz_t *exz, double e_rate, int64_t w_l, double ovlp_cut, int64_t force_aln, void *km)
{
    int64_t q_s, q_e, nw, k, q_l, t_l, t_tot_l, sec_check = (uref&&(!hpc_g))?1:0;
    int64_t aux_beg, aux_end, t_s, thre, aln_l, t_pri_l;
    char *q_string, *t_string; 
    z->w_list.n = 0; z->is_match = 0; z->align_length = 0;
    nw = get_num_wins(z->x_pos_s, z->x_pos_e+1, w_l);
    get_win_se_by_normalize_xs(z, (z->x_pos_s/w_l)*w_l, w_l, &q_s, &q_e);
    for (k = 0; k < nw; k++) {
        aux_beg = aux_end = 0; q_l = 1 + q_e - q_s;
        thre = q_l*e_rate; thre = Adjust_Threshold(thre, q_l);
        if(thre > THRESHOLD_MAX_SIZE) thre = THRESHOLD_MAX_SIZE;
        ///offset of y
        t_s = (q_s - z->x_pos_s) + z->y_pos_s;
        t_s += y_start_offset(q_s, &(z->f_cigar));

        aln_l = q_l + (thre<<1); t_tot_l = hpc_g?hpc_len(*hpc_g, z->y_id):uref->ug->u.a[z->y_id].len;
        if(init_waln(thre, t_s, t_tot_l, aln_l, &aux_beg, &aux_end, &t_s, &t_pri_l)) {
            q_string = qstr+q_s; 
            t_string = return_str_seq_exz(tstr, t_s, t_pri_l, z->y_pos_strand, hpc_g, uref, z->y_id);
            t_l = t_pri_l;
            // t_end = Reserve_Banded_BPM(t_string, aln_l, q_string, q_l, thre, &error);
            ed_band_cal_semi_64_w_absent_diag(t_string, t_l, q_string, q_l, thre, aux_beg, exz);
            // if(z->y_id == 178 && q_s == 86800 && q_e == 86807) {
            //     fprintf(stderr, "\n[M::%s::semi::t_s->%ld::t_pri_l->%ld::aux_beg->%ld::aux_end->%ld::thre->%ld] exz->ps::%d, exz->pe::%d, exz->ts::%d, exz->te::%d, exz->err::%d, exz->cigar.n::%d\n", 
            //     __func__, t_s, t_pri_l, aux_beg, aux_end, thre, exz->ps, exz->pe, exz->ts, exz->te, exz->err, (int32_t)exz->cigar.n);
            //     fprintf(stderr, "[tstr::len->%ld] %.*s\n", t_l, (int32_t)t_l, t_string);
            //     fprintf(stderr, "[qstr::len->%ld] %.*s\n", q_l, (int32_t)q_l, q_string);
            // } 
            if (is_align(*exz)) {
                // ed_band_cal_semi_64_w(t_string, aln_l, q_string, q_l, thre, exz);
                // assert(exz->err <= exz->thre);
                // fprintf(stderr, "[M::%s] exz->err::%d\n", __func__, exz->err);
                ///t_s do not have aux_beg, while t_s + t_end (aka, te) has 
                if(!push_wlst_exz(uref, hpc_g, NULL, z, qstr, tstr, exz, q_s, q_e, t_s, t_s + exz->pe, 
                                        t_tot_l, aux_beg, aux_end, e_rate, w_l, sec_check, ovlp_cut, force_aln, km)) {
                    return 0;
                }
                // append_window_list(z, q_s, q_e, t_s, t_s + t_end, error, aux_beg, aux_end, thre, w_l, km);
            }
        }
        q_s = q_e + 1; q_e = q_s + w_l - 1; 
        if(q_e >= (int64_t)z->x_pos_e) q_e = z->x_pos_e;
    }

    if((!force_aln) && (!simi_pass(z->x_pos_e+1-z->x_pos_s, z->align_length, 0, ovlp_cut, &e_rate)) &&
        (!simi_pass(z->x_pos_e+1-z->x_pos_s, z->align_length, sec_check, ovlp_cut, NULL))) return 0;
    return 1;
}

inline uint32_t ed_cut(const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, 
char *qstr, char *tstr, uint32_t rev, uint32_t id, 
int64_t qs, int64_t qe, int64_t t_s, int64_t block_s, double e_rate,
uint32_t aln_dir, int64_t* r_err, int64_t* qoff, int64_t* toff, int64_t* aln_qlen)
{
    (*aln_qlen) = 0; (*r_err) = INT32_MAX;
    if(qoff) (*qoff) = -1; if(toff) (*toff) = -1; 
    int64_t ql, aln_l, t_tot_l, aux_beg, aux_end, t_pri_l, thres;
    char *q_string, *t_string; unsigned int error; int t_end, q_end;

    ql = qe + 1 - qs;
    ///there are two potiential reasons for unmatched window:
    ///1. this window has a large number of differences
    ///2. DP does not start from the right offset
    if(rref) {
        thres = double_error_threshold(get_init_err_thres(ql, e_rate, block_s, THRESHOLD), ql);
    } else {
        thres = double_ul_error_threshold(get_init_err_thres(ql, e_rate, block_s, THRESHOLD_MAX_SIZE), ql);
    }
    
    aln_l = ql + (thres << 1);
    if(hpc_g) t_tot_l = hpc_len(*hpc_g, id);
    else if(uref) t_tot_l = uref->ug->u.a[id].len;
    else t_tot_l = Get_READ_LENGTH((*rref), id);

    if(!init_waln(thres, t_s, t_tot_l, aln_l, &aux_beg, &aux_end, &t_s, &t_pri_l)) return 0;
    // if(t_pri_l + thres < ql) return 0;

    q_string = qstr + qs; 
    if(rref) {
        fill_subregion(tstr, t_s, t_pri_l, rev, rref, id, aux_beg, aux_end); t_string = tstr;
    } else {
        t_string = return_str_seq(tstr, t_s, t_pri_l, rev, hpc_g, uref, id, aux_beg, aux_end);
    }

    // if(id == 6) {
    //     fprintf(stderr, "-[M::%s::aln_dir->%u] qs->%ld, ts->%ld, thres->%ld, aux_beg->%ld, aux_end->%ld, t_pri_l->%ld\n", 
    //     __func__, aln_dir, qs, t_s, thres, aux_beg, aux_end, t_pri_l);
    // }
    if(aln_dir == 0) {
        Reserve_Banded_BPM_Extension(t_string, aln_l, q_string, ql, thres, &error, &t_end, &q_end);
    } else {
        Reserve_Banded_BPM_Extension_REV(t_string, aln_l, q_string, ql, thres, &error, &t_end, &q_end);
    }

    if(t_end != -1 && q_end != -1) (*aln_qlen) = (aln_dir?(ql-q_end):(q_end+1));
    if(qoff) (*qoff) = q_end; if(toff) (*toff) = t_end; (*r_err) = error;

    if((*aln_qlen) == 0) return 0;
    return 1;
}

int64_t gen_extend_err_0(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, 
char *tstr, char *tstr_1, Correct_dumy* dumy, uint64_t *v_idx, int64_t block_s, double e_rate, 
int64_t qs, int64_t qe, int64_t pk)
{
    int64_t tot_e = 0, ts, di[2], al[2], tb[2], an = z->w_list.n; double rr;
    int64_t id = z->y_id, rev = z->y_pos_strand, ql = qe + 1 - qs;
    ///check if there are some windows that cannot be algined by any overlaps/unitigs
    ///if no, it is likely that the UL read itself has issues
    if(uref && v_idx && z->is_match == 4) {
        if(check_coverage_gap(v_idx, qs, qe, block_s)) {
            tot_e += THRESHOLD_MAX_SIZE; return tot_e;
        } 
    }
    ts = (qs - z->x_pos_s) + z->y_pos_s; ts += y_start_offset(qs, &(z->f_cigar));

    di[0] = di[1] = al[0] = al[1] = 0; tb[0] = tb[1] = -1;
    if((pk > 0) && (qs == (z->w_list.a[pk].x_end + 1))) {
        if(z->w_list.a[pk].clen == 0) {///do not have cigar
            gen_backtrace_adv(&(z->w_list.a[pk]), z, rref, hpc_g, uref, qstr, tstr, tstr_1, dumy, rev, id);
        }
        tb[0] = z->w_list.a[pk].y_end + 1;
    }

    if(((pk+1) < an) && ((qe+1) == (z->w_list.a[pk+1].x_start))) {
        if(z->w_list.a[pk+1].clen == 0) {///do not have cigar
            gen_backtrace_adv(&(z->w_list.a[pk+1]), z, rref, hpc_g, uref, qstr, tstr, tstr_1, dumy, rev, id);
        }
        tb[1] = z->w_list.a[pk+1].y_start-ql;
    }

    if(tb[0] == -1 && tb[1] == -1) tb[0] = tb[1] = ts;
    else if(tb[0] == -1 && tb[1] != -1) tb[0] = tb[1];
    else if(tb[1] == -1 && tb[0] != -1) tb[1] = tb[0];

    if(tb[0] != -1) {
        if(!ed_cut(uref, hpc_g, rref, qstr, tstr, rev, id, qs, qe, tb[0], block_s, e_rate,
                                                                    0, &(di[0]), NULL, NULL, &(al[0]))) {
            di[0] = ql; al[0] = 0;
        }
    }

    if(tb[1] != -1) {
        if(!ed_cut(uref, hpc_g, rref, qstr, tstr, rev, id, qs, qe, tb[1], block_s, e_rate,
                                                                    1, &(di[1]), NULL, NULL, &(al[1]))) {
            di[1] = ql; al[1] = 0;
        }
    }

    if(al[0] && al[1]) {///matched in both sides
        if((al[0] + al[1]) <= ql) {
            tot_e += di[0] + di[1] + ql - (al[0] + al[1]);
        } else {
            rr = ((double)ql)/((double)(al[0] + al[1]));
            tot_e += (di[0] + di[1])*rr;
        }
    } else if((!al[0]) && (!al[1])) {//failed
        tot_e += ql;
    } else if(al[0]) {
        tot_e += di[0] + (ql - al[0]);
    }else if(al[1]) {
        tot_e += di[1] + (ql - al[1]);
    }
    // if(z->y_id == 6) {
    //     fprintf(stderr, "-[M::%s] qs->%ld, ts->%ld, tb[0]->%ld, tb[1]->%ld, di[0]->%ld, di[1]->%ld, al[0]->%ld, al[1]->%ld, block_s->%ld, e_rate->%f\n", __func__, 
    //     qs, ts, tb[0], tb[1], di[0], di[1], al[0], al[1], block_s, e_rate);
    // }
    return tot_e;
}


int64_t gen_extend_err_0_exz(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, 
char *tstr, bit_extz_t *exz, uint64_t *v_idx, int64_t block_s, double e_rate, 
int64_t qs, int64_t qe, int64_t pk)
{
    int64_t tot_e = 0, ts, di[2], al[2], tb[2], an = z->w_list.n; double rr;
    int64_t id = z->y_id, rev = z->y_pos_strand, ql = qe + 1 - qs;
    ///check if there are some windows that cannot be algined by any overlaps/unitigs
    ///if no, it is likely that the UL read itself has issues
    if(uref && v_idx && z->is_match == 4) {
        if(check_coverage_gap(v_idx, qs, qe, block_s)) {
            tot_e += THRESHOLD_MAX_SIZE; return tot_e;
        } 
    }
    ts = (qs - z->x_pos_s) + z->y_pos_s; ts += y_start_offset(qs, &(z->f_cigar));

    di[0] = di[1] = al[0] = al[1] = 0; tb[0] = tb[1] = -1;
    if((pk > 0) && (qs == (z->w_list.a[pk].x_end + 1))) {
        if(z->w_list.a[pk].clen == 0) {///do not have cigar
            gen_backtrace_adv_exz(&(z->w_list.a[pk]), z, rref, hpc_g, uref, qstr, tstr, exz, rev, id);
        }
        tb[0] = z->w_list.a[pk].y_end + 1;
    }

    if(((pk+1) < an) && ((qe+1) == (z->w_list.a[pk+1].x_start))) {
        if(z->w_list.a[pk+1].clen == 0) {///do not have cigar
            gen_backtrace_adv_exz(&(z->w_list.a[pk+1]), z, rref, hpc_g, uref, qstr, tstr, exz, rev, id);
        }
        tb[1] = z->w_list.a[pk+1].y_start-ql;
    }

    if(tb[0] == -1 && tb[1] == -1) tb[0] = tb[1] = ts;
    else if(tb[0] == -1 && tb[1] != -1) tb[0] = tb[1];
    else if(tb[1] == -1 && tb[0] != -1) tb[1] = tb[0];

    if(tb[0] != -1) {
        if(!ed_cut(uref, hpc_g, rref, qstr, tstr, rev, id, qs, qe, tb[0], block_s, e_rate,
                                                                    0, &(di[0]), NULL, NULL, &(al[0]))) {
            di[0] = ql; al[0] = 0;
        }
    }

    if(tb[1] != -1) {
        if(!ed_cut(uref, hpc_g, rref, qstr, tstr, rev, id, qs, qe, tb[1], block_s, e_rate,
                                                                    1, &(di[1]), NULL, NULL, &(al[1]))) {
            di[1] = ql; al[1] = 0;
        }
    }

    if(al[0] && al[1]) {///matched in both sides
        if((al[0] + al[1]) <= ql) {
            tot_e += di[0] + di[1] + ql - (al[0] + al[1]);
        } else {
            rr = ((double)ql)/((double)(al[0] + al[1]));
            tot_e += (di[0] + di[1])*rr;
        }
    } else if((!al[0]) && (!al[1])) {//failed
        tot_e += ql;
    } else if(al[0]) {
        tot_e += di[0] + (ql - al[0]);
    }else if(al[1]) {
        tot_e += di[1] + (ql - al[1]);
    }
    // if(z->y_id == 6) {
    //     fprintf(stderr, "-[M::%s] qs->%ld, ts->%ld, tb[0]->%ld, tb[1]->%ld, di[0]->%ld, di[1]->%ld, al[0]->%ld, al[1]->%ld, block_s->%ld, e_rate->%f\n", __func__, 
    //     qs, ts, tb[0], tb[1], di[0], di[1], al[0], al[1], block_s, e_rate);
    // }
    return tot_e;
}

double gen_extend_err(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, 
char *tstr, char *tstr_1, Correct_dumy* dumy, uint64_t *v_idx, int64_t block_s, double ovlp_cut, double e_rate, double e_max, int64_t *r_e)
{
    int64_t ovl, k, ce, an = z->w_list.n, tot_l, tot_e, ws, we, ql; 
    ovl = z->x_pos_e+1-z->x_pos_s; if(r_e) (*r_e) = INT64_MAX;
    if(!simi_pass(ovl, z->align_length, 0, ovlp_cut, &e_rate)) return DBL_MAX;
    // nw = get_num_wins(z->x_pos_s, z->x_pos_e+1, block_s);
    for (k = 0; k < an; k++) {
        if(z->w_list.a[k].clen) z->w_list.a[k].y_end -= z->w_list.a[k].extra_begin;
    }

    tot_l = tot_e = 0;
    for (k = an-1, ce = z->x_pos_e; k >= 0; k--) {
        // assert(k == 0 || z->w_list.a[k].x_end > z->w_list.a[k-1].x_start);//sorted
        tot_l += z->w_list.a[k].x_end + 1 - z->w_list.a[k].x_start;
        tot_e += z->w_list.a[k].error;///matched window

        we = z->w_list.a[k].x_end;
        while (we < ce) {
            ws = we+1;
            get_win_id_by_s(z, ws, block_s, &we);
            ql = we+1-ws; tot_l += ql; 
            tot_e += gen_extend_err_0(z, uref, hpc_g, rref, qstr, tstr, tstr_1, dumy, v_idx, block_s, e_rate, ws, we, k);
            if((e_max > 0) && (tot_e > (ovl*e_max))) return DBL_MAX;
        }
        ce = z->w_list.a[k].x_start-1;
        if((e_max > 0) && (tot_e > (ovl*e_max))) return DBL_MAX;
    }

    if(ce >= ((int64_t)z->x_pos_s)) {
        we = ((int64_t)z->x_pos_s)-1;
        while (we < ce) {
            ws = we+1;
            get_win_id_by_s(z, ws, block_s, &we);
            ql = we+1-ws; tot_l += ql; 
            tot_e += gen_extend_err_0(z, uref, hpc_g, rref, qstr, tstr, tstr_1, dumy, v_idx, block_s, e_rate, ws, we, k);
            if((e_max > 0) && (tot_e > (ovl*e_max))) return DBL_MAX;
        }
    }

    assert(tot_l == ovl); if(r_e) (*r_e) = tot_e;
    return (double)(tot_e)/(double)(tot_l);
}

double gen_extend_err_exz(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, 
char *tstr, bit_extz_t *exz, uint64_t *v_idx, int64_t block_s, double ovlp_cut, double e_rate, double e_max, int64_t *r_e)
{
    int64_t ovl, k, ce, an = z->w_list.n, tot_l, tot_e, ws, we, ql; 
    ovl = z->x_pos_e+1-z->x_pos_s; if(r_e) (*r_e) = INT64_MAX;
    if(!simi_pass(ovl, z->align_length, 0, ovlp_cut, &e_rate)) return DBL_MAX;
    // nw = get_num_wins(z->x_pos_s, z->x_pos_e+1, block_s);
    // for (k = 0; k < an; k++) {
    //     if(z->w_list.a[k].clen) z->w_list.a[k].y_end -= z->w_list.a[k].extra_begin;
    // }

    tot_l = tot_e = 0;
    for (k = an-1, ce = z->x_pos_e; k >= 0; k--) {
        // assert(k == 0 || z->w_list.a[k].x_end > z->w_list.a[k-1].x_start);//sorted
        tot_l += z->w_list.a[k].x_end + 1 - z->w_list.a[k].x_start;
        tot_e += z->w_list.a[k].error;///matched window

        we = z->w_list.a[k].x_end;
        while (we < ce) {
            ws = we+1;
            get_win_id_by_s(z, ws, block_s, &we);
            ql = we+1-ws; tot_l += ql; 
            tot_e += gen_extend_err_0_exz(z, uref, hpc_g, rref, qstr, tstr, exz, v_idx, block_s, e_rate, ws, we, k);
            if((e_max > 0) && (tot_e > (ovl*e_max))) return DBL_MAX;
        }
        ce = z->w_list.a[k].x_start-1;
        if((e_max > 0) && (tot_e > (ovl*e_max))) return DBL_MAX;
    }

    if(ce >= ((int64_t)z->x_pos_s)) {
        we = ((int64_t)z->x_pos_s)-1;
        while (we < ce) {
            ws = we+1;
            get_win_id_by_s(z, ws, block_s, &we);
            ql = we+1-ws; tot_l += ql; 
            tot_e += gen_extend_err_0_exz(z, uref, hpc_g, rref, qstr, tstr, exz, v_idx, block_s, e_rate, ws, we, k);
            if((e_max > 0) && (tot_e > (ovl*e_max))) return DBL_MAX;
        }
    }

    assert(tot_l == ovl); if(r_e) (*r_e) = tot_e;
    return (double)(tot_e)/(double)(tot_l);
}


void push_anchors(window_list *z, window_list_alloc *zidx, asg64_v *anchor, uint64_t *qhp, int64_t qhp_l, int64_t *qhp_k, uint32_t mcl)
{
    int64_t xi = 0, yi = 0, ci, cn = z->clen; uint8_t c = (uint8_t)-1; uint32_t cl = (uint32_t)-1;
    for (ci = 0; ci < cn; ci++) {
        get_cigar_cell(z, zidx, ci, &c, &cl);
        if (c == 0) { //match   
            if(cl >= mcl) {
                ;
                ;
                ;
                ;
            }
            xi += cl; yi += cl;
        } else if (c == 1) {
            xi += cl; yi += cl;
        } else if (c == 2) {///y has more bases than x
            yi += cl;
        } else if (c == 3) {///x has more bases than y
            xi += cl;
        }
    }
}

#define gen_hpc_max_len(x) ((x)+((x)>>1)+1)
///[off_s, off_e)
uint64_t extract_mm_hpc(char *in, int64_t len, int64_t off_s, int64_t off_e, int64_t w, uint64_t rev)
{
    int64_t i, o, l, trim, k, tl; uint64_t m, sf; uint8_t c;
    ///forward
    for (k = 1, trim = 0; k <= w; k++) {
        m = 0; o = gen_hpc_max_len(k); sf = k<<1;
        if(!rev) {
            ///[off_s, off_e)
            for (i = ((off_s>=o)?(off_s-o):(0)), l = 0; i < off_e; i++) {
                c = seq_nt4_table[(uint8_t)in[i]]; 
                if((c < 4) && (((l >= k) && (((m>>sf)&3) == c)) || (l < k))) {
                    if(l < k) m = (m<<2) + c; 
                    else sf = (sf?(sf):(k<<1))-2;
                    l++;
                } else {
                    if(i > off_s) {
                        tl = i-off_s;
                        if((l >= o) && (trim < tl)) trim = tl;
                        l = -1;
                        break;
                    }
                    l = 0; sf = k<<1;
                }
            }
            tl = i-off_s;
            if((l!=-1) && (i > off_s) && (l >= o) && (trim < tl)) {
                trim = tl;
                if(trim >= (off_e-off_s)) break;
            }
        } else {
            ///[off_s, off_e)
            for (i = (((len-off_e)>=o)?(off_e+o):(len))-1, l = 0; i >= off_s; i--) {
                c = seq_nt4_table[(uint8_t)in[i]]; 
                if((c < 4) && (((l >= k) && (((m>>sf)&3) == c)) || (l < k))) {
                    if(l < k) m = (m<<2) + c; 
                    else sf = (sf?(sf):(k<<1))-2;
                    l++;
                } else {
                    if(i+1 < off_e) {
                        tl = off_e-i-1;
                        if((l >= o) && (trim < tl)) trim = tl;
                        l = -1;
                        break;
                    }
                    l = 0; sf = k<<1;
                }
            }
            tl = off_e-i-1;
            if((l!=-1) && (i+1 < off_e) && (l >= o) && (trim < tl)) {
                trim = tl;
                if(trim >= (off_e-off_s)) break;
            }
        }
    }
    return trim;
}

uint64_t trim_hpc(const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, char *tstr, int64_t hpc_max,
int64_t ql, int64_t tl, int64_t tid, int64_t trev, int64_t *rqs, int64_t *rqe, int64_t *rts, int64_t *rte)
{
    ///[qs, qe); [ts, te)
    int64_t qs = *rqs, qe = *rqe, ts = *rts, te = *rte, trim[2], hl, aux_l, subl = qe-qs; char *ss;
    if(hpc_max > 32) hpc_max = 32; trim[0] = trim[1] = 0; aux_l = gen_hpc_max_len(hpc_max);

    qs -= aux_l; if(qs < 0) qs = 0;
    qe += aux_l; if(qe > ql) qe = ql;
    ss = qstr + qs;
    hl = extract_mm_hpc(ss, qe - qs, (*rqs)-qs, (*rqe)-qs, hpc_max, 0);
    if(hl >= subl) return 0; trim[0] = hl;
    hl = extract_mm_hpc(ss, qe - qs, (*rqs)-qs, (*rqe)-qs, hpc_max, 1);
    if(hl >= subl) return 0; trim[1] = hl;
    if(trim[0] + trim[1] >= subl) return 0;

    ts -= aux_l; if(ts < 0) ts = 0;
    te += aux_l; if(te > tl) te = tl;
    if(rref) {
        fill_subregion(tstr, ts, te-ts, trev, rref, tid, 0, 0); ss = tstr;
    } else {
        ss = return_str_seq(tstr, ts, te-ts, trev, hpc_g, uref, tid, 0, 0);
    }
    hl = extract_mm_hpc(ss, te - ts, (*rts)-ts, (*rte)-ts, hpc_max, 0);
    if(hl >= subl) return 0; if(hl > trim[0]) trim[0] = hl;
    hl = extract_mm_hpc(ss, te - ts, (*rts)-ts, (*rte)-ts, hpc_max, 1);
    if(hl >= subl) return 0; if(hl > trim[1]) trim[1] = hl;
    if(trim[0] + trim[1] >= subl) return 0;

    (*rqs) += trim[0]; (*rts) += trim[0];
    (*rqe) -= trim[1]; (*rte) -= trim[1];
    return 1;
}

#define cl_pushp(type, v, p) do {									\
		if ((v).length == (v).size) {										\
			(v).size = (v).size? (v).size<<1 : 2;							\
			(v).list = (type*)realloc((v).list, sizeof(type) * (v).size);	\
		}															\
		*(p) = &(v).list[(v).length++]; \
	} while (0)

///ai is the suffix of aj
int64_t inline traceback_sc(const k_mer_hit *ai, const k_mer_hit *aj)
{
    int64_t qsi = ai->self_offset-ai->cnt, qej = aj->self_offset;
    int64_t tsi = ai->offset-ai->cnt, tej = aj->offset;
    if(qsi >= qej && tsi >= tej) return ai->cnt;
    return INT32_MIN;
}

void split_long_anchors(Candidates_list *ac, int64_t block, int64_t block_n)
{
    int64_t i, m, an = ac->length;
    if(block_n < 0) {
        for (i = block_n = 0; i < an; i++) {
            if(ac->list[i].cnt <= block) block_n++;
            else block_n += (ac->list[i].cnt/block) + (((ac->list[i].cnt%block) > 0)?1:0);
        }
    }
    if(block_n <= an) return;
    if(block_n > ac->size) {
        ac->size = block_n; REALLOC(ac->list, ac->size);
    }
    for (i = an-1, m = block_n-1; i >= 0; i--) {
        if(ac->list[i].cnt <= block) {
            ac->list[m--] = ac->list[i];
        } else {
            while (ac->list[i].cnt > 0) {
                ac->list[m] = ac->list[i];
                if(ac->list[i].cnt >= block) {
                    ac->list[m].cnt = block;
                    ac->list[i].cnt -= block;
                    ac->list[i].self_offset -= block;
                    ac->list[i].offset -= block;
                } else {
                    ac->list[i].cnt = 0; 
                }
                m--;
            }
        }
    }
    ac->length = block_n; assert(m == -1);
}

int64_t gen_affine_traceback_dp(Candidates_list *ac, int64_t max_skip, int64_t max_iter, int64_t max_dis, int64_t block, int64_t block_n)
{
    if(ac->length < 1) return 0;
    int64_t i, j, *p, *t, max_f, n_skip, max_j, end_j, st, max_ii, sc, max, tmp, msc_i, msc; 
    int32_t *f, cL; k_mer_hit* a = ac->list; int64_t a_n = ac->length; Chain_Data* dp;

    for (i = 1; i < a_n; ++i) {
        sc = traceback_sc(&a[i], &a[i-1]);
        if(sc == INT32_MIN) break;
    }
    if(i >= a_n) return a_n;

    split_long_anchors(ac, block, block_n);
    a = ac->list; a_n = ac->length; dp = &(ac->chainDP);
    
    resize_Chain_Data(dp, a_n, NULL); t = dp->tmp; f = dp->score; p = dp->pre;
    t[0] = 0; p[0] = -1; f[0] = a[0].cnt;
    msc_i = msc = -1; i = 0;

    memset(t, 0, (a_n*sizeof((*t))));
    for (i = st = 0, max_ii = -1; i < a_n; ++i) {
        max_f = a[i].cnt; n_skip = 0; max_j = end_j = -1;
        if ((i-st) > max_iter) st = i-max_iter;

        for (j = i - 1; j >= st; --j) {
            sc = traceback_sc(&a[i], &a[j]);
            if(sc == INT32_MIN) break;
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
			tmp = traceback_sc(&a[i], &a[max_ii]);
			if (tmp != INT32_MIN && max_f < tmp + f[max_ii])
				max_f = tmp + f[max_ii], max_j = max_ii;
		}
        f[i] = max_f; p[i] = max_j;
        if ((max_ii < 0) || (((((int64_t)a[i].offset)-((int64_t)a[max_ii].offset))<=max_dis) && (f[max_ii]<f[i]))) {
            max_ii = i;
        }

        if(f[i] > msc) {
            msc = f[i]; msc_i = i;
        }
    }

    cL = 0; i = msc_i; 
    while (i >= 0) {
        t[cL++] = i; i = p[i];
    }
    for (i = 0; i < cL; i++) a[i] = a[t[cL-i-1]];
    return cL;
}

uint64_t gen_affine_traceback(overlap_region *o, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, char *tstr, uint64_t ql, 
uint64_t *qhp, uint64_t qhp_l, Candidates_list *ac, uint32_t hpc_max, uint32_t min_ach, uint32_t block)
{
    if(o->w_list.n <= 0) return 0;
    int64_t nw = o->w_list.n, snw, k, t, qi, ti, ci, p_qi, p_ti, cn, qs, qe, ts, te, tl; 
    window_list *z; uint8_t c; uint32_t cl, pcl, id = o->y_id, rev = o->y_pos_strand; k_mer_hit *p;
    uint64_t hm, aocc = 0, bocc = 0; qi = ti = 0; p_qi = p_ti = INT32_MIN; pcl = 0;
    clear_Candidates_list(ac); 
    if(o->w_list.n > (uint64_t)ac->size) {
        ac->size = o->w_list.n; REALLOC(ac->list, ac->size);
    }

    if(hpc_g) tl = hpc_len(*hpc_g, id);
    else if(uref) tl = uref->ug->u.a[id].len;
    else tl = Get_READ_LENGTH((*rref), id);
    for (k = 0; k < nw; k++) {
        z = &(o->w_list.a[k]); ci = 0; cn = z->clen;
        qi = z->x_start; ti = z->y_start;
        for (ci = 0; ci < cn; ci++) {
            get_cigar_cell(z, &(o->w_list), ci, &c, &cl);
            if (c == 0) { //match   
                if((p_qi == qi) && (p_ti == ti)) {
                    pcl += cl;
                } else {
                    ///push
                    if(pcl > 0) {
                        hm = 0; qs = qi - pcl; qe = qi; ts = ti - pcl; te = ti;
                        if(pcl > min_ach) hm = trim_hpc(uref, hpc_g, rref, qstr, tstr, hpc_max, ql, tl, id, rev, &qs, &qe, &ts, &te);
                        if(hm) aocc++;
                        if(hm || aocc == 0) {
                            cl_pushp(k_mer_hit, *ac, &p);
                            p->readID = p->strand = !!hm;
                            p->cnt = qe - qs; p->self_offset = qe; p->offset = te;
                            if(p->cnt > min_ach) bocc++;
                        }
                    }
                    ///push
                    pcl = 0;
                }
                qi += cl; ti += cl;
                p_qi = qi; p_ti = ti;
            } else {
                ///push
                if(pcl > 0) {
                    hm = 0; qs = qi - pcl; qe = qi; ts = ti - pcl; te = ti;
                    if(pcl > min_ach) hm = trim_hpc(uref, hpc_g, rref, qstr, tstr, hpc_max, ql, tl, id, rev, &qs, &qe, &ts, &te);
                    if(hm) aocc++;
                    if(hm || aocc == 0) {
                        cl_pushp(k_mer_hit, *ac, &p);
                        p->readID = p->strand = !!hm;
                        p->cnt = qe - qs; p->self_offset = qe; p->offset = te;
                        if(p->cnt > min_ach) bocc++;
                    }
                }
                pcl = 0;
                ///push
                if (c == 1) {
                    qi += cl; ti += cl; 
                } if (c == 2) {///t has more bases than p
                    ti += cl; 
                } if (c == 3) {///p has more bases than t
                    qi += cl; 
                }
            } 
        }
    }

    ///push
    if(pcl > 0) {
        hm = 0; qs = qi - pcl; qe = qi; ts = ti - pcl; te = ti;
        if(pcl > min_ach) hm = trim_hpc(uref, hpc_g, rref, qstr, tstr, hpc_max, ql, tl, id, rev, &qs, &qe, &ts, &te);
        if(hm) aocc++;
        if(hm || aocc == 0) {
            cl_pushp(k_mer_hit, *ac, &p);
            p->readID = p->strand = !!hm;
            p->cnt = qe - qs; p->self_offset = qe; p->offset = te;
            if(p->cnt > min_ach) bocc++;
        }
    }
    ///push

    nw = ac->length; snw = -1;
    if(aocc > 0) {
        if(aocc < (uint64_t)ac->length) {
            for (k = t = snw = 0; k < nw; k++) {
                if(ac->list[k].readID) {
                    ac->list[t] = ac->list[k]; 
                    if(ac->list[t].cnt <= block) snw++;
                    else snw += (ac->list[t].cnt/block) + (((ac->list[t].cnt%block) > 0)?1:0);
                    t++;
                }
            }
            ac->length = t;
        }
    } else if (bocc > 0) {
        if(bocc < (uint64_t)ac->length) {
            for (k = t = snw = 0; k < nw; k++) {
                if(ac->list[k].cnt > min_ach) {
                    ac->list[t] = ac->list[k]; 
                    if(ac->list[t].cnt <= block) snw++;
                    else snw += (ac->list[t].cnt/block) + (((ac->list[t].cnt%block) > 0)?1:0);
                    t++;
                }
            }
            ac->length = t;
        }
    }
    nw = ac->length;
    gen_affine_traceback_dp(ac, 25, 5000, 5000, block, snw);
    return 1;
}


void align_ul_ed(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, char* qstr, char *tstr, double e_rate, int64_t w_l, void *km)
{
    int64_t q_s, q_e, nw, k, q_l, t_tot_l;
    int64_t aux_beg, aux_end, t_s, thre, aln_l, t_pri_l, t_end;
    char *q_string, *t_string; unsigned int error;
    z->w_list.n = 0; z->is_match = 0;
    nw = get_num_wins(z->x_pos_s, z->x_pos_e+1, w_l);
    get_win_se_by_normalize_xs(z, (z->x_pos_s/w_l)*w_l, w_l, &q_s, &q_e);
    for (k = 0; k < nw; k++) {
        aux_beg = aux_end = 0; q_l = 1 + q_e - q_s;
        thre = q_l*e_rate; thre = Adjust_Threshold(thre, q_l);
        if(thre > THRESHOLD_MAX_SIZE) thre = THRESHOLD_MAX_SIZE;
        ///offset of y
        t_s = (q_s - z->x_pos_s) + z->y_pos_s;
        t_s += y_start_offset(q_s, &(z->f_cigar));

        aln_l = q_l + (thre<<1); t_tot_l = hpc_g?hpc_len(*hpc_g, z->y_id):uref->ug->u.a[z->y_id].len;
        if(init_waln(thre, t_s, t_tot_l, aln_l, &aux_beg, &aux_end, &t_s, &t_pri_l)) {
            q_string = qstr+q_s; 
            t_string = return_str_seq(tstr, t_s, t_pri_l, z->y_pos_strand, hpc_g, uref, z->y_id, aux_beg, aux_end);

            t_end = Reserve_Banded_BPM(t_string, aln_l, q_string, q_l, thre, &error);
            if (error!=((unsigned int)-1)) {
                z->align_length += q_l;
                ///t_s do not have aux_beg, while t_s + t_end (aka, te) has 
                append_window_list(z, q_s, q_e, t_s, t_s + t_end, error, aux_beg, aux_end, thre, w_l, km);
            }
        }
        q_s = q_e + 1; q_e = q_s + w_l - 1; 
        if(q_e >= (int64_t)z->x_pos_e) q_e = z->x_pos_e;
    }

    assert(q_e == (int64_t)z->x_pos_e);
}

uint64_t realign_ed(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, 
char *tstr, char *tstr_1, Correct_dumy* dumy, kvec_t_u64_warp* v_idx, int64_t block_s, double e_rate, 
double *e_rate_final, uint32_t sec_check, double ovlp_cut, int64_t *is_sort)
{
    int64_t i, k, nw, a_nw, w_id, y_id, y_strand, real_y_start, x_start, x_end, x_len, ce, cs;
    int64_t w_s, w_e, mm_we, mm_ws, mm_aln, ovl, y_readLen, total_y_start, total_y_end;
    uint64_t *w_idx = NULL, srt = 1; window_list *p = NULL; if(sec_check && (!uref)) sec_check = 0;
    z->is_match = 0; if(is_sort) (*is_sort) = 1;
    if(z->w_list.n == 0) return 0;
    nw = get_num_wins(z->x_pos_s, z->x_pos_e+1, block_s); a_nw = z->w_list.n;
    y_id = z->y_id; y_strand = z->y_pos_strand; 
    ovl = z->x_pos_e+1-z->x_pos_s; mm_ws = mm_we = z->x_pos_s; mm_aln = 0; 

    if(hpc_g) y_readLen = hpc_len(*hpc_g, y_id);
    else if(uref) y_readLen = uref->ug->u.a[y_id].len;
    else y_readLen = Get_READ_LENGTH((*rref), y_id);
    
    for (i = a_nw-1, ce = z->x_pos_e; i >= 0; i--) { //utilize the the end pos of pre-window in forward
        w_e = mm_we = z->w_list.a[i].x_end; 
        total_y_start = z->w_list.a[i].y_end + 1 - z->w_list.a[i].extra_begin;
        while ((w_e < ce) && (total_y_start < y_readLen)) {
            w_s = w_e + 1;
            w_id = get_win_id_by_s(z, w_s, block_s, &w_e);
            x_start = w_s; x_end = w_e;
            if(aln_wlst_adv(z, rref, hpc_g, uref, qstr, tstr, tstr_1, dumy,
            y_strand, y_id, x_start, x_end, total_y_start, block_s, e_rate, 0)) {
                p = &(z->w_list.a[z->w_list.n-1]);
                mm_we = x_end;
            } else {
                break;
            }
            total_y_start = p->y_end + 1 - p->extra_begin;
        }
        ce = z->w_list.a[i].x_start-1;
        if(i == a_nw-1) {///only possiblity with the largest end pos
            mm_aln = mm_we+1-mm_ws;
            if(!simi_pass(ovl, mm_aln, sec_check, ovlp_cut, NULL)) break;
        }
    }
     
    if(z->w_list.a[a_nw-1].x_end > z->w_list.a[z->w_list.n-1].x_end) {
        srt = 0; if(is_sort) (*is_sort) = 0;
    }
    if(i >= 0) return 0;
    if((!srt) && (z->w_list.n <= (nw*0.2))) {///if very few windows are mapped
        radix_sort_window_list_xs_srt(z->w_list.a, z->w_list.a + z->w_list.n);
        srt = 1; if(is_sort) (*is_sort) = 1;
    }
    if(!srt) {///need sort
        a_nw = z->w_list.n; 
        kv_resize(uint64_t, v_idx->a, (uint64_t)nw); w_idx = v_idx->a.a;
        memset(v_idx->a.a, -1, sizeof((*v_idx->a.a))*nw); 
        for (i = 0; i < a_nw; i++) { ///w_idx[] == (uint64_t) if unmatched
            assert(z->w_list.a[i].y_end != -1);
            w_id = get_win_id_by_s(z, z->w_list.a[i].x_start, block_s, NULL);
            w_idx[w_id] = i;
        }
        ///deal with first window
        mm_ws = z->x_pos_s; 
        if(w_idx[0] != (uint64_t)-1) {
            w_s = z->w_list.a[w_idx[0]].x_start; mm_aln -= (w_s-mm_ws);
            mm_ws = z->w_list.a[w_idx[0]].x_end+1;
        }        
        for (i = 1; i < nw; i++) { //utilize the the start pos of next window in backward
            ///find the first matched window, which should not be the first window
            ///the pre-window of this matched window must be unmatched
            if(w_idx[i] != (uint64_t)-1 && w_idx[i-1] == (uint64_t)-1) {
                w_s = z->w_list.a[w_idx[i]].x_start; mm_aln -= (w_s-mm_ws);
                ///check if the start pos of this matched window has been calculated
                if(z->w_list.a[w_idx[i]].clen == 0) {
                    p = &(z->w_list.a[w_idx[i]]);
                    gen_backtrace_adv(p, z, rref, hpc_g, uref, qstr, tstr, tstr_1, dumy, y_strand, y_id);
                    assert(p->error != -1);
                    p->y_end += p->extra_begin;
                } 
                real_y_start = p->y_start;

                ///the end pos for pre window is real_y_start - 1
                total_y_end = real_y_start - 1;
                ///find the unmatched window on the left of current matched window
                ///k starts from i - 1
                for (k = i - 1; k >= 0 && w_idx[k] == (uint64_t)-1 && total_y_end > 0; k--) {  
                    w_e = w_s - 1;
                    w_id = get_win_id_by_e(z, w_e, block_s, &w_s);
                    assert(w_id == k);
                    x_start = w_s; x_end = w_e; x_len = x_end + 1 - x_start;
                    if(aln_wlst_adv(z, rref, hpc_g, uref, qstr, tstr, tstr_1, dumy,
                        y_strand, y_id, x_start, x_end, total_y_end+1-x_len, block_s, e_rate, 1)) {
                        p = &(z->w_list.a[z->w_list.n-1]);
                        p->y_start -= p->extra_begin; ///y_start has no shift, but y_end has shift  
                        w_idx[k] = z->w_list.n - 1;
                        mm_aln += x_len;
                        // if(is_sort && (*is_sort) && z->w_list.n > 1 && p->x_start < z->w_list.a[z->w_list.n-2].x_start) (*is_sort) = 0;
                    } else {
                        break;
                    }
                    total_y_end = p->y_start - 1;
                }
                if(!simi_pass(ovl, mm_aln, sec_check, ovlp_cut, NULL)) break;
            }
            if(w_idx[i] != (uint64_t)-1) mm_ws = z->w_list.a[w_idx[i]].x_end+1;
        }
        if(i < nw) return 0;
    } else {//sorted
        a_nw = z->w_list.n; mm_ws = z->x_pos_s; 
        for (i = 0, cs = z->x_pos_s; i < a_nw; i++) {
            p = &(z->w_list.a[i]);
            w_s = p->x_start; mm_aln -= (w_s-mm_ws);
            ///check if the start pos of this matched window has been calculated
            if((w_s > cs) && (p->clen == 0)) {
                gen_backtrace_adv(p, z, rref, hpc_g, uref, qstr, tstr, tstr_1, dumy, y_strand, y_id);
                assert(p->error != -1);
                p->y_end += p->extra_begin;
            } 
            real_y_start = p->y_start;
            ///the end pos for pre window is real_y_start - 1
            total_y_end = real_y_start - 1;
            while ((w_s > cs) && (total_y_end > 0)) {
                w_e = w_s - 1;
                w_id = get_win_id_by_e(z, w_e, block_s, &w_s);
                x_start = w_s; x_end = w_e; x_len = x_end + 1 - x_start;
                if(aln_wlst_adv(z, rref, hpc_g, uref, qstr, tstr, tstr_1, dumy,
                    y_strand, y_id, x_start, x_end, total_y_end+1-x_len, block_s, e_rate, 1)) {
                    p = &(z->w_list.a[z->w_list.n-1]);
                    p->y_start -= p->extra_begin; ///y_start has no shift, but y_end has shift  
                    mm_aln += x_len;
                    // if(is_sort && (*is_sort) && z->w_list.n > 1 && p->x_start < z->w_list.a[z->w_list.n-2].x_start) (*is_sort) = 0;
                } else {
                    break;
                }
                total_y_end = p->y_start - 1;
            }
            if(!simi_pass(ovl, mm_aln, sec_check, ovlp_cut, NULL)) break;
            mm_ws = cs = z->w_list.a[i].x_end+1;
        }
        if(a_nw < (int64_t)z->w_list.n) {
            srt = 0; if(is_sort) (*is_sort) = 0;
        }
        if(i < a_nw) return 0;
    }

    if(e_rate_final) {
        /**
        if(simi_pass(ovl, z->align_length, 0, &e_rate)) {
            a_nw = z->w_list.n;
            for (i = 0; i < a_nw; i++) {
                p = &(z->w_list.a[i]);
                ///check if the cigar of this window has been got 
                if(p->clen == 0) {
                    gen_backtrace_adv(p, z, rref, hpc_g, uref, qstr, tstr, tstr_1, dumy, y_strand, y_id);
                    assert(p->error != -1);
                } else {
                    p->y_end -= p->extra_begin;
                }
            }
            if(!srt) radix_sort_window_list_xs_srt(z->w_list.a, z->w_list.a + z->w_list.n);
            ///note: this function will change tstr/qstr
            error_rate = non_trim_error_rate(z, rref, uref, v_idx, dumy, g_read, e_rate, block_s);
            z->is_match = 0;///must be here;

            if (error_rate <= e_rate_final) {
                overlap_list->mapped_overlaps_length += ovl;
                z->is_match = 1; append_unmatched_wins(z, block_s);
                if(rref) {
                    calculate_boundary_cigars(z, rref, dumy, g_read, e_rate);
                } else {
                    calculate_ul_boundary_cigars(z, uref, dumy, g_read, e_rate, block_s);
                }
                // assert(get_num_wins(z->x_pos_s, z->x_pos_e+1, block_s)==(int64_t)z->w_list.n);
                // assert((int64_t)z->x_pos_s==z->w_list.a[0].x_start && 
                //                             (int64_t)z->x_pos_e==z->w_list.a[z->w_list.n-1].x_end);
            } else if (error_rate <= e_rate_final * 1.5) {
                z->is_match = 3;
            }
        }
        **/
       return 1;
    } else {
        return 1;
    }
    return 0;
}

uint64_t col_errors(overlap_region *z)
{
    uint64_t i, e = 0;
    for (i = 0; i < z->w_list.n; i++) e += z->w_list.a[i].error;
    return e;
}


void ul_lalign_hpc(overlap_region_alloc* ol, Candidates_list *cl, const ul_idx_t *uref, char *qstr, 
                        uint64_t ql, UC_Read* qu, UC_Read* tu, Correct_dumy* dumy, 
                        haplotype_evdience_alloc* hap, kvec_t_u64_warp* v_idx, kvec_t_u32_warp* q_idx,  
                        double e_rate, double eh_rate, int64_t wl, void *km)
{
    uint64_t i, qhl, bs, k, ovl, whl; Window_Pool w; double err; overlap_region t; overlap_region *z;
    whl = MIN((((double)THRESHOLD_MAX_SIZE)/eh_rate), WINDOW);
    ol->mapped_overlaps_length = 0;
    if(ol->length <= 0) return;
    
    ///hpc alignment
    ///init hpc seq
    qhl = update_ol_track(ol, cl, uref->hpc_g, uref, 1, ql, qstr, q_idx);
    gen_hpc_str(qstr, ql, qu, &qhl);
    ///verify hpc seq
    clear_Correct_dumy(dumy, ol, km); err = eh_rate; 
    init_Window_Pool(&w, qhl, whl, (int)(1.0/err));
    bs = (w.window_length)+(THRESHOLD_MAX_SIZE<<1)+1;
    resize_UC_Read(tu, bs<<1);
    for (i = k = 0; i < ol->length; i++) {
        if(!align_ul_ed_post(&(ol->list[i]), uref, uref->hpc_g, qu->seq, tu->seq, tu->seq+bs, dumy, err, w.window_length, OVERLAP_THRESHOLD_FILTER_HPC, km)) {
            continue;
        }
        // fprintf(stderr, "+++[M::%s] yid::%u, x::[%u, %u), y::[%u, %u), aln::%u, err::%lu\n", __func__, ol->list[i].y_id, 
        // ol->list[i].x_pos_s, ol->list[i].x_pos_e+1, ol->list[i].y_pos_s, ol->list[i].y_pos_e+1,
        // ol->list[i].align_length, col_errors(&(ol->list[i])));
        if(k != i) {
            t = ol->list[k];
            ol->list[k] = ol->list[i];
            ol->list[i] = t;
        }
        k++;
    }
    ol->length = k;
    if(ol->length <= 0) return;


    ///base alignment
    update_ol_track(ol, cl, NULL, uref, 1, ql, NULL, NULL);
    resize_UC_Read(qu, ql); qu->length = ql; memcpy(qu->seq, qstr, ql);
    clear_Correct_dumy(dumy, ol, km); err = e_rate; 
    init_Window_Pool(&w, ql, wl, (int)(1.0/err));
    bs = (w.window_length)+(THRESHOLD_MAX_SIZE<<1)+1;
    resize_UC_Read(tu, bs<<1);
    for (i = 0; i < ol->length; i++) {
        z = &(ol->list[i]); ovl = z->x_pos_e+1-z->x_pos_s;
        if(!align_ul_ed_post(z, uref, NULL, qu->seq, tu->seq, tu->seq+bs, dumy, err, w.window_length, -1, km)) {
            continue;
        }
        // fprintf(stderr, "---[M::%s] yid::%u, x::[%u, %u), y::[%u, %u), aln::%u, err::%lu\n", __func__, ol->list[i].y_id, 
        // ol->list[i].x_pos_s, ol->list[i].x_pos_e+1, ol->list[i].y_pos_s, ol->list[i].y_pos_e+1,
        // ol->list[i].align_length, col_errors(&(ol->list[i])));
        if(uref && simi_pass(ovl, z->align_length, uref?1:0, -1, NULL)) {
            z->is_match = 3; ol->mapped_overlaps_length += z->align_length;
        }
    }

    if(uref && ol->mapped_overlaps_length > 0) {
        set_herror_win(ol, dumy, v_idx, err, ql, w.window_length);
    }

    double e_max = err*1.5, rr;
    for (i = k = 0; i < ol->length; i++) {
        z = &(ol->list[i]); ovl = z->x_pos_e + 1 - z->x_pos_s;
        rr = gen_extend_err(z, uref, NULL, NULL, qu->seq, tu->seq, tu->seq+bs,
                    dumy, v_idx?v_idx->a.a:NULL, w.window_length, -1, err, (e_max+0.000001), NULL);
        z->is_match = 0;///must be here;
        if (rr <= err) {
            for (k = 0; k < z->w_list.n; k++) {
                if(z->w_list.a[k].clen) continue;
                gen_backtrace_adv(&(z->w_list.a[k]), z, NULL, NULL, uref, qu->seq, tu->seq, tu->seq+bs,
                dumy, z->y_pos_strand, z->y_id);
            }
            
            ol->mapped_overlaps_length += ovl; k++;
            z->is_match = 1; append_unmatched_wins(z, w.window_length);
            calculate_ul_boundary_cigars(z, uref, dumy, qu, err, w.window_length);
        } else if (rr <= e_max) {
            z->is_match = 3;
        }
    }

    partition_ul_overlaps_advance(ol, uref, qu, tu, dumy, hap, 1, err, w.window_length, km);
    
    // recalcate_window(overlap_list, R_INF, g_read, dumy, overlap_read);
    // partition_overlaps(overlap_list, R_INF, g_read, dumy, hap, force_repeat);
    // recalcate_window_ul_advance(overlap_list, uref, g_read, dumy, overlap_read, max_ov_diff_ec, w_inf.window_length, km);
    // recalcate_window_advance(overlap_list, NULL, uref, g_read, dumy, overlap_read, v_idx, w_inf.window_length, max_ov_diff_ec, max_ov_diff_ec);
    /**
    refine_ed_aln(overlap_list, NULL, uref, g_read, dumy, overlap_read, v_idx, w_inf.window_length, max_ov_diff_ec, max_ov_diff_ec);
    **/
    // fprintf(stderr, "[M::%s-beg] occ[0]->%lu, occ[1]->%lu, occ[2]->%lu, occ[3]->%lu\n", __func__, 
    // ovlp_occ(overlap_list, 0), ovlp_occ(overlap_list, 1), ovlp_occ(overlap_list, 2), ovlp_occ(overlap_list, 3));
    ///after this function, overlap_list is sorted by x_pos_e; used for g_chain
    /**
    partition_ul_overlaps_advance(overlap_list, uref, g_read, overlap_read, dumy, hap, force_repeat, max_ov_diff_ec, w_inf.window_length, km);
    **/
    // print_ovlp_occ_stat(overlap_list, g_read->length, 1);
    // print_ovlp_occ_stat(overlap_list, g_read->length, 2);
    // fprintf(stderr, "[M::%s-end] occ[0]->%lu, occ[1]->%lu, occ[2]->%lu, occ[3]->%lu\n", __func__, 
    // ovlp_occ(overlap_list, 0), ovlp_occ(overlap_list, 1), ovlp_occ(overlap_list, 2), ovlp_occ(overlap_list, 3));
    // debug_phasing_status(overlap_list, uref->ug, 0, hap, g_read, 20, 1176);
    // debug_phasing_status(overlap_list, uref->ug, 0, hap, g_read, 20, 1167);
    // debug_phasing_status(overlap_list, uref->ug, 0, hap, g_read, 20, 1170);
    /**
    
    
    if(is_consensus)
    {
        generate_consensus(overlap_list, R_INF, g_read, dumy, g, DAGCon, current_cigar, second_round);
    }


    (*fully_cov) = check_if_fully_covered(overlap_list, R_INF, g_read, dumy, g, abnormal);
    **/
}

// void ul_phase(overlap_region *oa, int64_t on, uint64_t *idx, uint64_t *buf)
// {
//     int64_t i, oi; overlap_region *z;
//     for (i = 0; i < on; i++) {
//         oi = (uint32_t)idx[i]; z = &(oa[oi]);
//         if(z->is_match != 1) continue;
//     }
// }

inline uint32_t cigar_check_dbg(char *pstr, char *tstr, bit_extz_t *ez)
{
	int32_t pi = ez->ps, ti = ez->ts, err = 0; uint32_t ci = 0, cl, k; uint16_t c;
	while (ci < ez->cigar.n) {
		ci = pop_trace(&(ez->cigar), ci, &c, &cl);
		// fprintf(stderr, "# %u = %u, cigar_n::%u\n", c, cl, (uint32_t)ez->cigar.n); 
		if(c == 0) {
			for (k=0;(k<cl)&&(pstr[pi]==tstr[ti]);k++,pi++,ti++);
			if(k!=cl) {
				fprintf(stderr, "ERROR-d-0\n"); 
				return 0;
			}
		} else {
			err += cl;
			if(c == 1) {
				for (k=0;(k<cl)&&(pstr[pi]!=tstr[ti]);k++,pi++,ti++);
				if(k!=cl) {
					fprintf(stderr, "ERROR-d-1\n"); 
					return 0;
				}
			} else if(c == 2) {///more p
				pi+=cl;
			} else if(c == 3) {
				ti+=cl;
			}
		}
	}
	if(err != ez->err) {
		fprintf(stderr, "ERROR-err, err::%d, ez->err::%d\n", err, ez->err); 
		return 0;
	}
	return 1;
}

int64_t wcigar_check(window_list *p, window_list_alloc *z, char *qstr, char *tstr)
{
    bit_extz_t ez;
    ez.cigar.a = z->c.a+p->cidx; ez.cigar.n = ez.cigar.m = p->clen;
    ez.ps = p->y_start; ez.pe = p->y_end;
    ez.ts = p->x_start; ez.te = p->x_end;
    ez.err = p->error;
    // prt_cigar(ez.cigar.a, ez.cigar.n);
    return cigar_check_dbg(tstr, qstr, &ez);
}

void verify_aln(int32_t sid, overlap_region *z, UC_Read* qu, UC_Read* tu, All_reads *rref, hpc_t *hpc_g, const ul_idx_t *uref)
{
    uint64_t tl, tid = z->y_id, rev = z->y_pos_strand, i; char *tstr, *qstr; window_list *p;
    if(hpc_g) tl = hpc_len(*hpc_g, tid);
    else if(uref) tl = uref->ug->u.a[tid].len;
    else tl = Get_READ_LENGTH((*rref), tid);
    
    qstr = qu->seq; resize_UC_Read(tu, tl);
    if(rref) {
        recover_UC_Read_sub_region(tu->seq, 0, tl, rev, rref, tid); tstr = tu->seq;
    } else {
        tstr = return_str_seq_exz(tu->seq, 0, tl, rev, hpc_g, uref, tid);
    }
    
    for (i = 0; i < z->w_list.n; i++) {
        p = &(z->w_list.a[i]);
        if(p->y_end == -1) continue;
        if(!wcigar_check(p, &(z->w_list), qstr, tstr)) break;
    }
    if(i < z->w_list.n) {
        fprintf(stderr, "sid::%d, tid::%lu, i::%lu, qs::%d, qe::%d, ts::%d, te::%d, err::%d\n", sid, tid, i, 
        p->x_start, p->x_end, p->y_start, p->y_end, p->error); 
        prt_cigar(z->w_list.c.a+p->cidx, p->clen);
        fprintf(stderr, "[tstr::[%u, %u]] %.*s\n", z->y_pos_s, z->y_pos_e, p->y_end+1-p->y_start, tstr+p->y_start);
        fprintf(stderr, "[qstr::[%u, %u]] %.*s\n", z->x_pos_s, z->x_pos_e, p->x_end+1-p->x_start, qstr+p->x_start);
        exit(0);
    }
}


void ul_lalign_old_ed(overlap_region_alloc* ol, Candidates_list *cl, const ul_idx_t *uref, char *qstr, 
                        uint64_t ql, UC_Read* qu, UC_Read* tu, Correct_dumy* dumy, 
                        haplotype_evdience_alloc* hap, kvec_t_u64_warp* v_idx,   
                        double e_rate, int64_t wl, uint64_t is_base, void *km)
{
    uint64_t i, bs, k, ovl/**, on**/; Window_Pool w; double err; 
    /**int64_t sc;**/ overlap_region t; overlap_region *z; 
    ol->mapped_overlaps_length = 0;
    if(ol->length <= 0) return;

    ///base alignment
    clear_Correct_dumy(dumy, ol, km); err = e_rate; 
    init_Window_Pool(&w, ql, wl, (int)(1.0/err));
    bs = (w.window_length)+(THRESHOLD_MAX_SIZE<<1)+1;
    resize_UC_Read(tu, bs<<1);

    if(is_base) {    
        resize_UC_Read(qu, ql); qu->length = ql; memcpy(qu->seq, qstr, ql);
        for (i = 0; i < ol->length; i++) {
            z = &(ol->list[i]); ovl = z->x_pos_e+1-z->x_pos_s;
            if(!align_ul_ed_post(z, uref, NULL, qu->seq, tu->seq, tu->seq+bs, dumy, err, w.window_length, -1, km)) {
                continue;
            }
            if(uref && simi_pass(ovl, z->align_length, uref?1:0, -1, NULL)) {
                z->is_match = 3; ol->mapped_overlaps_length += z->align_length;
            }
        }

        if(uref && ol->mapped_overlaps_length > 0) {
            set_herror_win(ol, dumy, v_idx, err, ql, w.window_length);
        }

        double e_max = err*1.5, rr; int64_t re;
        for (i = k = 0; i < ol->length; i++) {
            z = &(ol->list[i]); ovl = z->x_pos_e + 1 - z->x_pos_s;
            rr = gen_extend_err(z, uref, NULL, NULL, qu->seq, tu->seq, tu->seq+bs,
                        dumy, v_idx?v_idx->a.a:NULL, w.window_length, -1, err, (e_max+0.000001), &re);
            z->is_match = 0;///must be here;
            if (rr <= err) {
                if(k != i) {
                    t = ol->list[k];
                    ol->list[k] = ol->list[i];
                    ol->list[i] = t;
                }
                ol->list[k].is_match = 1; ol->list[k].non_homopolymer_errors = re;
                k++;
            } 
        }

        ol->length = k;
        // fprintf(stderr, "+[M::%s] on::%lu\n", __func__, ol->length);
        if(ol->length <= 0) return;
    } else {
        // fprintf(stderr, "-[M::%s] on::%lu\n", __func__, ol->length);
        if(ol->length <= 1) return;
        // for (i = 0; (i < ol->length) && (ol->list[i].is_match == 1); i++); on = i;
        // if(on <= 1) return;
        // kv_resize(uint64_t, v_idx->a, (on<<1)); v_idx->a.n = on;
        // for (i = 0; i < on; i++) {
        //     sc = ol->list[i].x_pos_e+1-ol->list[i].x_pos_s;
        //     sc -= ((int64_t)(ol->list[i].non_homopolymer_errors*ERROR_RATE));
        //     if(sc < 0) sc = 0;
        //     v_idx->a.a[i] = sc; v_idx->a.a[i] <<= 32; v_idx->a.a[i] += i;
        // }
        // radix_sort_bc64(v_idx->a.a, v_idx->a.a+on);
        // for (i = 0; i < on; i++) {
        //     ol->list[i].non_homopolymer_errors = 0;
        //     k = ((uint32_t)v_idx->a.a[i]);
        //     v_idx->a.a[k]<<=32; v_idx->a.a[k]>>=32; v_idx->a.a[k]|=i;
        // }
        // ul_phase(ol->list, on, v_idx->a.a, v_idx->a.a+on);


        for (i = 0; i < ol->length; i++) {
            z = &(ol->list[i]); ovl = z->x_pos_e+1-z->x_pos_s; z->is_match = 1;
            for (k = 0; k < z->w_list.n; k++) {
                if(z->w_list.a[k].clen) continue;
                gen_backtrace_adv(&(z->w_list.a[k]), z, NULL, NULL, uref, qu->seq, tu->seq, tu->seq+bs,
                dumy, z->y_pos_strand, z->y_id);
            }
            ol->mapped_overlaps_length += ovl; 
            append_unmatched_wins(z, w.window_length);
            calculate_ul_boundary_cigars(z, uref, dumy, qu, err, w.window_length);
        }
        partition_ul_overlaps_advance(ol, uref, qu, tu, dumy, hap, 1, err, w.window_length, km);
    }
}

inline uint64_t scale_ed_thre(uint32_t err, uint32_t max_err)
{
    uint64_t bd = (err<<1)+1, w; 
    w = (bd>>bitw); w <<= bitw; if(w < bd) w += bitwbit;
    err = (w-1)>>1; if(err > max_err) err = max_err;
    return err;
}


///[qs, qe)
int64_t update_semi_coord(const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, overlap_region *z, 
int64_t qs, int64_t qe, int64_t thre, int64_t *ts, int64_t *te, int64_t *aux_beg)
{
    int64_t ql = qe - qs, aln_l, t_tot_l, id = z->y_id, aux_end, tl; 
    (*ts) = (qs - z->x_pos_s) + z->y_pos_s;
    (*ts) += y_start_offset(qs, &(z->f_cigar));
    aln_l = ql + (thre<<1);
    if(hpc_g) t_tot_l = hpc_len(*hpc_g, id);
    else if(uref) t_tot_l = uref->ug->u.a[id].len;
    else t_tot_l = Get_READ_LENGTH((*rref), id);
    if(!init_waln(thre, (*ts), t_tot_l, aln_l, aux_beg, &aux_end, ts, &tl)) {
        (*ts) = (*te) = (*aux_beg) = -1;
        return 0;
    }
    (*te) = (*ts) + tl; 
    return 1;
}


void adjust_ext_offset(int64_t *qs, int64_t *qe, int64_t *ts, int64_t *te, int64_t ql, int64_t tl, int64_t thre, int64_t mode)
{
    int64_t qoff, toff;
    if(mode == 1) {///forward extension
        qoff = ql - (*qs); toff = tl - (*ts);
        if(qoff <= toff) {
            (*qe) = ql; (*te) = (*ts) + qoff + thre;
        } else {
            (*te) = tl; (*qe) = (*qs) + toff + thre;
        }
    } else if(mode == 2) {///backward extension
        qoff = (*qe); toff = (*te);
        if(qoff <= toff) {
            (*qs) = 0; (*ts) = (*te) - qoff - thre;
        } else {
            (*ts) = 0; (*qs) = (*qe) - toff - thre;
        }
    }
    if((*qs) < 0) (*qs) = 0;
    if((*ts) < 0) (*ts) = 0;
    if((*qe) > ql) (*qe) = ql;
    if((*te) > tl) (*te) = tl;
}

void adjust_specific_ext_offset(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, bit_extz_t *exz, char* qstr, UC_Read *tu, 
int64_t *qs0, int64_t *qe0, int64_t *ts0, int64_t *te0, int64_t ql, int64_t tl, int64_t wl, int64_t *mode)
{   
    ///[qs, qe)
    if((*mode) != 1 && (*mode) != 2) return;
    int64_t qs = *qs0, qe = *qe0, ts = *ts0, te = *te0, k, update, wid, we, ws, sl, gq, gt, gg, min_g = INT32_MAX, min_id = -1, gl;
    adjust_ext_offset(&qs, &qe, &ts, &te, ql, tl, 0, *mode);
    if((*mode) == 1) {///forward extension
        //qs0 and ts0 are fixed; te0 = -1, qe0 is unreliable
        if(qe > (*qe0)) {
            we = (*qe0); we/=wl; we *= wl; we +=wl; we--; ///next window
            for (k = 0; we < qe && k < ext_w; we+=wl) {//[ws, we]; [qs, qe)
                wid = get_win_id_by_e(z, we, wl, NULL);
                if(z->w_list.a[wid].y_end == -1) continue;//unmapped
                if(z->w_list.a[wid].x_end < (*qs0)) continue;
                if(z->w_list.a[wid].y_end < (*ts0)) continue;
                sl = z->w_list.a[wid].x_end+1-z->w_list.a[wid].x_start;
                if(z->w_list.a[wid].error > (sl/A_L)) continue;
                gq = z->w_list.a[wid].x_end-(*qs0); 
                gt = z->w_list.a[wid].y_end-(*ts0); 
                gl = MIN(gq, gt); gl/=16; if(gl <= 0) gl = 1;
                gg = (gq>=gt)?(gq-gt):(gt-gq); gg /= gl;
                update = 0;
                if(min_g>gg) {
                    update = 1;
                } else if((min_g==gg)&&(z->w_list.a[min_id].error>z->w_list.a[wid].error)) {
                    update = 1;
                } 
                if(update) {
                    min_g = gg; min_id = wid;
                }
                k++;
            }
        }
        if(min_id >= 0) {
            (*qe0) = z->w_list.a[min_id].x_end+1; 
            if(((*qe0)+wl) >= qe) {
                (*qe0) = qe; (*te0) = te; ///still extension
                return;
            }
            (*te0) = z->w_list.a[min_id].y_end+1; (*mode) = 0;///global
        } else {
            (*qe0) = qe; (*te0) = te; ///still extension
        }
    } else if((*mode) == 2) {///backward extension
        //qe0 and te0 are fixed; ts0 = -1, qs0 is unreliable
        if(qs < (*qs0)) {
            ws = (*qs0)-1; ws/=wl; ws*=wl; 
            for (k = 0; ws >= qs && k < ext_w; ws-=wl) {//[ws, we]; [qs, qe)
                wid = get_win_id_by_s(z, ws, wl, NULL);
                if(z->w_list.a[wid].y_end == -1) continue;//unmapped
                if(z->w_list.a[wid].x_start >= (*qe0)) continue;
                if(z->w_list.a[wid].y_start >= (*te0)) continue;
                sl = z->w_list.a[wid].x_end+1-z->w_list.a[wid].x_start;
                if(z->w_list.a[wid].error > (sl/A_L)) continue;
                gq = (*qe0) - z->w_list.a[wid].x_start; 
                gt = (*te0) - z->w_list.a[wid].y_start; 
                gl = MIN(gq, gt); gl/=16; if(gl <= 0) gl = 1;
                gg = (gq>=gt)?(gq-gt):(gt-gq); gg /= gl;
                update = 0;
                if(min_g>gg) {
                    update = 1;
                } else if((min_g==gg)&&(z->w_list.a[min_id].error>z->w_list.a[wid].error)) {
                    update = 1;
                } 
                if(update) {
                    min_g = gg; min_id = wid;
                }
                k++;
            }
        }
        if(min_id >= 0) {
            (*qs0) = z->w_list.a[min_id].x_start; 
            if((qs+wl) >= (*qs0)) {
                (*qs0) = qs; (*ts0) = ts; ///still extension
                return;
            }
            (*ts0) = z->w_list.a[min_id].y_start; (*mode) = 0;///global
        } else {
            (*qs0) = qs; (*ts0) = ts; ///still extension
        }
    }
}

int64_t cal_exz_infi(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, bit_extz_t *exz, char* qstr, UC_Read *tu, int64_t qs, int64_t qe, int64_t ts, int64_t te, int64_t thre, int64_t q_tot_l, int64_t mode)
{
    clear_align(*exz);
    int64_t aux_beg = 0, bd = (((thre)<<1)+1), ql, tl, t_tot_l = -1; int32_t nword = ((bd>>bitw)+(!!(bd&bitz)));
    char *q_string, *t_string; int32_t rev = z->y_pos_strand, id = z->y_id; ql = qe - qs;
    if(hpc_g) t_tot_l = hpc_len(*hpc_g, id);
    else if(uref) t_tot_l = uref->ug->u.a[id].len;
    else t_tot_l = Get_READ_LENGTH((*rref), id);

    if(mode == 3) {
        update_semi_coord(uref, hpc_g, rref, z, qs, qe, thre, &ts, &te, &aux_beg);
    } else if(mode == 1 || mode == 2) {
        adjust_ext_offset(&qs, &qe, &ts, &te, q_tot_l, t_tot_l, thre, mode);
    }
    
    if((qe > qs) && (te > ts) && (ts != -1) && (te != -1)) {
        ql = qe - qs; q_string = qstr + qs;
        tl = te - ts; resize_UC_Read(tu, tl);
        // fprintf(stderr, "q::[%ld, %ld), t::[%ld, %ld), thre::%ld, t_tot_l::%ld\n", qs, qe, ts, te, thre, t_tot_l);
        if(rref) {
            recover_UC_Read_sub_region(tu->seq, ts, tl, rev, rref, id); t_string = tu->seq;
        } else {
            t_string = return_str_seq_exz(tu->seq, ts, tl, rev, hpc_g, uref, id);
        }
        // fprintf(stderr, ", nword::%d", nword);
        // if(ql < 0) fprintf(stderr, "qs::%ld, qe::%ld\n", qs, qe);
        // return 0;
        if(nword <= 1) {
            if(mode == 0) { //global
                ed_band_cal_global_64_w_trace(t_string, tl, q_string, ql, thre, exz);
            } else if(mode == 1) {///forward extension
                // fprintf(stderr, "q::[%ld, %ld), t::[%ld, %ld), thre::%ld\n", qs, qe, ts, te, thre);
                ed_band_cal_extension_64_0_w_trace(t_string, tl, q_string, ql, thre, exz);
            } else if(mode == 2) {///backward extension
                ed_band_cal_extension_64_1_w_trace(t_string, tl, q_string, ql, thre, exz);
            } else if(mode == 3) {//semi-global
                ed_band_cal_semi_64_w_absent_diag_trace(t_string, tl, q_string, ql, thre, aux_beg, exz);
            }
        } else if(nword == 2) {
            if(mode == 0) { //global
                ed_band_cal_global_128_w_trace(t_string, tl, q_string, ql, thre, exz);
            } else if(mode == 1) {///forward extension
                ed_band_cal_extension_128_0_w_trace(t_string, tl, q_string, ql, thre, exz);
            } else if(mode == 2) {///backward extension
                ed_band_cal_extension_128_1_w_trace(t_string, tl, q_string, ql, thre, exz);
            } else if(mode == 3) {//semi-global
                ed_band_cal_semi_128_w_absent_diag_trace(t_string, tl, q_string, ql, thre, aux_beg, exz);
            }
        } else {
            if(mode == 0) { //global
                ed_band_cal_global_infi_w_trace(t_string, tl, q_string, ql, thre, &nword, exz);
            } else if(mode == 1) {///forward extension
                ed_band_cal_extension_infi_0_w_trace(t_string, tl, q_string, ql, thre, &nword, exz);
            } else if(mode == 2) {///backward extension
                ed_band_cal_extension_infi_1_w_trace(t_string, tl, q_string, ql, thre, &nword, exz);
            } else if(mode == 3) {//semi-global
                ed_band_cal_semi_infi_w_absent_diag_trace(t_string, tl, q_string, ql, thre, aux_beg, &nword, exz);
            }
        }
        if(is_align(*exz)) {

            return 1;
        }
        return 0;
    }
    return 0;
}

void hc_aln_exz(overlap_region *z, Candidates_list *cl, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, 
char* qstr, UC_Read *tu, int64_t qs, int64_t qe, int64_t ts, int64_t te, int64_t estimate_err,
int64_t mode, int64_t wl, bit_extz_t *exz, int64_t q_tot, double e_rate, int64_t maxl, int64_t maxe)
{
    clear_align(*exz);
    int64_t thre, ql = qe - qs, thre0, t_tot, id = z->y_id; 
    if(((ts == -1) && (te == -1))) mode = 3;///set to semi-global
    
    if(mode == 1 || mode == 2) {
        if(hpc_g) t_tot = hpc_len(*hpc_g, id);
        else if(uref) t_tot = uref->ug->u.a[id].len;
        else t_tot = Get_READ_LENGTH((*rref), id);
        // fprintf(stderr, "\n+[M::%s::ql->%ld::tl->%ld::mode->%ld] q::[%ld, %ld), t::[%ld, %ld), zx[%d, %d), zy[%d, %d)\n", 
        // __func__, q_tot, t_tot, mode, qs, qe, ts, te, z->x_pos_s, z->x_pos_e+1, z->y_pos_s, z->y_pos_e+1);
        ///find a aligned window >= qe to convert extension to global; or do extension directly
        adjust_specific_ext_offset(z, uref, hpc_g, rref, exz, qstr, tu, &qs, &qe, &ts, &te, q_tot, t_tot, wl, &mode);
        // fprintf(stderr, "-[M::%s::ql->%ld::tl->%ld::mode->%ld] q::[%ld, %ld), t::[%ld, %ld), zx[%d, %d), zy[%d, %d)\n", 
        // __func__, q_tot, t_tot, mode, qs, qe, ts, te, z->x_pos_s, z->x_pos_e+1, z->y_pos_s, z->y_pos_e+1);
    }

    // fprintf(stderr, "[M::%s::ql::%ld] qs::%ld, qe::%ld, ts::%ld, te::%ld, mode::%ld, estimate_err::%ld, e_rate::%f", 
    //                                                     __func__, ql, qs, qe, ts, te, mode, estimate_err, e_rate);

    if(ql <= maxl && (estimate_err*1.2) <= maxe) {
        thre = scale_ed_thre(estimate_err, maxe); if(thre > ql) thre = ql;
        if(cal_exz_infi(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, thre, q_tot, mode)) {
            // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(+)\n", exz->err, exz->thre, thre);
            return;
        }

        thre0 = thre; thre = ql*e_rate;
        thre = scale_ed_thre(thre, maxe); if(thre > ql) thre = ql;
        if(thre > thre0) {
            if(cal_exz_infi(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, thre, q_tot, mode)) {
                // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(-)\n", exz->err, exz->thre, thre);
                return;
            }
        }

        thre0 = thre; thre <<= 1;
        thre = scale_ed_thre(thre, maxe); if(thre > ql) thre = ql;
        if(thre > thre0) {
            if(cal_exz_infi(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, thre, q_tot, mode)) {
                // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(-)\n", exz->err, exz->thre, thre);
                return;
            }
        }

        thre0 = thre; thre = ql*0.51;
        thre = scale_ed_thre(thre, maxe); if(thre > ql) thre = ql;
        if(thre > thre0) {
            if(cal_exz_infi(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, thre, q_tot, mode)) {
                // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(*)\n", exz->err, exz->thre, thre);
                return;
            }
        }
    }
    // fprintf(stderr, ", err::%d, thre::%d\n", INT32_MAX, exz->thre);
    // anchor_aln(z, cl, uref, hpc_g, rref, qstr, tu, qs, qe, ts, te, thre, mode, wl, exz, q_tot, A_L);
    
}

void sub_ciagar_gen(overlap_region *z, Candidates_list *cl, uint64_t s, uint64_t e, uint64_t wl, 
const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, UC_Read *tu, bit_extz_t *exz, double e_rate,
int64_t ql, uint64_t rid)
{
    uint64_t qs, qe, sid, eid, k, l, m, tot_e, c_e; int64_t q[2], t[2], o[2], mode, kocc;
    qs = (s/wl)*wl; if(qs < z->x_pos_s) qs = z->x_pos_s; if(qs > z->x_pos_e) return;
    qe = (e/wl)*wl; if(qe < e) qe += wl; if(qe > z->x_pos_e+1) qe = z->x_pos_e+1;
    if(qe <= 0) return;
    sid = get_win_id_by_s(z, qs, wl, NULL);
    eid = get_win_id_by_e(z, qe-1, wl, NULL) + 1;///must qe-1 instead of qe!!!!!!
    if(sid >= eid) return;
    // fprintf(stderr, "\n***[M::%s::rid->%lu] s::%lu, e::%lu, n_qs::%lu, n_qe::%lu, z::[%u, %u), sid::%lu, eid::%lu, w_list.n::%lu\n", 
    // __func__, rid, s, e, qs, qe, z->x_pos_s, z->x_pos_e+1, sid, eid, (uint64_t)z->w_list.n);
    for (k = sid+1, l = sid; k <= eid; k++) {//[sid, eid)
        if(k == eid || z->w_list.a[k].extra_end < 0) {
            if(k - l > 1 || z->w_list.a[l].extra_end >= 0) {
                q[0] = q[1] = t[0] = t[1] = -1; 
                o[0] = o[1] = -1; mode = -1; tot_e = 0;
                if(z->w_list.a[l].extra_end < 0) {
                    q[0] = z->w_list.a[l].x_end+1; 
                    if(z->w_list.a[l].y_end != -1) {
                        t[0] = z->w_list.a[l].y_end+1; 
                    }
                    if(z->w_list.a[l].extra_end != INT16_MIN) {
                        o[0] = -z->w_list.a[l].extra_end;
                    }
                } else {///first window
                    q[0] = qs; 
                    if(z->w_list.a[l].y_end != -1) {
                        c_e = z->w_list.a[l].error;
                    } else {
                        c_e = z->w_list.a[l].x_end + 1 - z->w_list.a[l].x_start;
                        if(c_e > THRESHOLD_MAX_SIZE) c_e = THRESHOLD_MAX_SIZE;
                    }
                    tot_e += c_e;
                }

                if(k > sid && k < eid && z->w_list.a[k].extra_end < 0) {
                    q[1] = z->w_list.a[k-1].x_end+1; 
                    if(z->w_list.a[k].y_end != -1) {
                        // if(z->w_list.a[k-1].y_end == -1) {
                        //     fprintf(stderr, "[M::%s::rid::%lu] k::%lu, sid::%lu, eid::%lu, k_y_end::%d, k-1_y_end::%d, xk[%d, %d)\n", 
                        //     __func__, rid, k, sid, eid, z->w_list.a[k].y_end, z->w_list.a[k-1].y_end,
                        //     z->w_list.a[k].x_start, z->w_list.a[k].x_end+1);
                        // }
                        assert(z->w_list.a[k-1].y_end != -1);
                        t[1] = z->w_list.a[k-1].y_end+1;
                    }
                    if(z->w_list.a[k].extra_end != INT16_MIN) {
                        o[1] = -z->w_list.a[k].extra_end;
                    }
                } else {///last window
                    q[1] = qe;
                }

                if((t[0] != -1) && (t[1] != -1)) {
                    mode = 0;//global
                } else if((t[0] != -1) && (t[1] == -1)) {
                    /**t[1] = z->y_pos_e+1;**/ mode = 1;///forward extension
                } else if((t[0] == -1) && (t[1] != -1)) {
                    /**t[0] = z->y_pos_s;**/ mode = 2;///backward extension
                } else {
                    mode = 3;//semi-global
                }

                for (m = l+1; m < k; m++) {
                    if(z->w_list.a[m].y_end != -1) {
                        c_e = z->w_list.a[m].error;
                    } else {
                        c_e = z->w_list.a[m].x_end + 1 - z->w_list.a[m].x_start;
                        if(c_e > THRESHOLD_MAX_SIZE) c_e = THRESHOLD_MAX_SIZE;
                    }
                    tot_e += c_e;
                }
                // if(q[1] < q[0]) {
                //     fprintf(stderr, "[M::%s::ql::%lu] qs::%lu, qe::%lu, ts::%lu, te::%lu, mode::%ld, tot_e::%lu\n", 
                //                                         __func__, q[1]-q[0], q[0], q[1], t[0], t[1], mode, tot_e);
                // }
                kocc = MAX(o[0], o[1]); if(kocc < 0) kocc = 1;
                hc_aln_exz(z, cl, uref, hpc_g, rref, qstr, tu, q[0], q[1], t[0], t[1], tot_e, mode, wl, exz, ql, e_rate, MAX_CNS_L, MAX_CNS_E); 
                
            }
            l = k;
        }
    }
}


uint64_t cigar_gen(overlap_region *z, Candidates_list *cl, ul_ov_t *ov, uint64_t on, uint64_t qn, uint64_t wl,  
const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, UC_Read *tu, bit_extz_t *exz, double e_rate, 
int64_t ql, uint64_t rid, ul_ov_t *des)
{
    if(on <= 0) return 0;
    uint64_t i; 
    // des[0] = ov[0]; 
    // for (i = m = 1; i < on; i++) {
    //     fusion_merge();
    // }
    

    for (i = 0; i < on; i++) {
        assert((i<=0)||(ov[i].qs > ov[i-1].qe));
        sub_ciagar_gen(z, cl, ov[i].qs, ov[i].qe, wl, uref, hpc_g, rref, qstr, tu, exz, e_rate, ql, rid);
    }
    return on;
    // uint64_t i, qs = (uint64_t)-1, qe = (uint64_t)-1;
    // for (i = 0; i < on && ov[i].qn == qn; i++) {
    //     assert((i<=0)||(ov[i].qs >= ov[i-1].qe));
    //     if(ov[i].qs <= qe && qe != (uint64_t)-1) {
    //         qe = ov[i].qe;
    //     } else {
    //         if(qs != (uint64_t)-1) sub_ciagar_gen(z, cl, qs, qe, wl, uref, hpc_g, rref, qstr, tu, exz, e_rate, ql, rid);
    //         qs = ov[i].qs; qe = ov[i].qe;
    //     }
    // }
    // if(qs != (uint64_t)-1) sub_ciagar_gen(z, cl, qs, qe, wl, uref, hpc_g, rref, qstr, tu, exz, e_rate, ql, rid);
}



int64_t push_adp_k_hits(Candidates_list *cl, int64_t cln, uint64_t qs, uint64_t qe, uint64_t ts, uint64_t te, uint64_t readID, int64_t ci, int64_t dbgid)
{
    int64_t k = ci, m = -1; k_mer_hit *p = NULL;
    if(cl->length > cln) {
        m = cl->length-1; assert(is_alnw(cl->list[m]));
    }
    for (; k < cln && cl->list[k].readID == readID && cl->list[k].self_offset < qs; k++) {
        p = &(cl->list[k]);
        // fprintf(stderr, "[M::%s::] p->q::%ld, p->t::%ld\n", __func__, p->self_offset, p->offset);
        if(p->offset >= ts) continue;
        if((m>=0) && ((cl->list[m].offset >= p->offset) || (cl->list[m].self_offset >= p->self_offset))) {
            continue;
        }
        kv_pushp_cl(k_mer_hit, (*cl), &p); *p = cl->list[k]; p->readID = p->cnt;
        // if(dbgid == 109111 || dbgid == 75436) {
        //     fprintf(stderr, "[M::%s::k->%ld] self_offset::%u, offset::%u, cl->length::%lld, cln::%ld, p->readID::%u\n", 
        //     __func__, k, cl->list[k].self_offset, cl->list[k].offset, cl->length, cln, p->readID);
        // }
    }
    if(qs == (uint64_t)-1 || qe == (uint64_t)-1) return k;
    // if(!(k >= cln || cl->list[k].self_offset >= qe)) {
    //      fprintf(stderr, "[M::%s::] k::%ld, self_offset::%u, offset::%u\n", 
    //                 __func__, k, cl->list[k].self_offset, cl->list[k].offset);
    // }
    assert(k >= cln || cl->list[k].readID != readID || cl->list[k].self_offset >= qe);
    ///push qs, ts
    kv_pushp_cl(k_mer_hit, (*cl), &p); 
    p->readID = ((uint32_t)(0x7fffffff));
    p->cnt = (uint32_t)-1; p->strand = 0; 
    p->self_offset = qs; p->offset = ts;
    
    for (; k < cln && cl->list[k].readID == readID && cl->list[k].self_offset < qe; k++);

    ///push qe, te
    kv_pushp_cl(k_mer_hit, (*cl), &p); 
    p->readID = ((uint32_t)(0x7fffffff));
    p->cnt = (uint32_t)-1; p->strand = 1; 
    p->self_offset = qe-1; p->offset = te-1;
    return k;
}

int64_t gen_weight_khits0(uint32_t qs, uint32_t qe, k_mer_hit *a, int64_t an, int64_t k, uint64_t dp)///[qs, qe)
{
    for (; k >= 0 && a[k].self_offset >= qs; k--); 
    for (k = ((k>=0)?k:0); k < an && a[k].self_offset < qe; k++) {
        if((a[k].self_offset >= qs) && (a[k].self_offset < qe) && (!is_alnw(a[k]))) {
            a[k].readID = dp;
        }
    }
    return k;
}

void gen_weight_khits(asg64_v* idx, k_mer_hit *a, int64_t an)
{
    int64_t i, idx_n = idx->n, dp, old_dp, beg, end, k;
    // fprintf(stderr, "[M::%s::] idx->n::%ld\n", __func__, (int64_t)idx->n);
    for (i = k = 0, dp = old_dp = 0, beg = 0, end = -1; i < idx_n; ++i) {///[beg, end) but coordinates in idx is [, ]
        ///if idx->a.a[] is qe
        old_dp = dp;
        if ((idx->a[i]>>32)&1) {
            --dp; end = (idx->a[i]>>33)+1;
        }else {
            //meet a new overlap; the overlaps are pushed by the x_pos_s
            ++dp; end = (idx->a[i]>>33);
        }
        // fprintf(stderr, "[M::%s::] beg::%ld, end::%ld, old_dp::%ld\n", __func__, beg, end, old_dp);
        if(end > beg) k = gen_weight_khits0(beg, end, a, an, k, old_dp);
        beg = end;
    }
}

void prt_khit(Candidates_list *cl, overlap_region_alloc* ol, overlap_region *z, uint64_t utg_id, const char *cmd)
{
    uint64_t k, cid; int64_t i; 
    if(!z) {
        for (k = 0; k < ol->length; k++) {
            z = &(ol->list[k]);
            if(z->y_id == utg_id) break;
        }
    }
    if((z) || (z->y_id == utg_id)) {
            fprintf(stderr, "[M::%s::cmd->%s] ******\n", __func__, (char *)cmd);
        for (i = z->shared_seed, cid = cl->list[i].readID; i < cl->length && cl->list[i].readID == cid; i++) {
            fprintf(stderr, "[M::%s::i->%lu] qoff::%u, toff::%u, cid::%lu, cl->length::%lld\n", __func__, i, 
                        cl->list[i].self_offset, cl->list[i].offset, cid, cl->length);
        }
    }
}

int64_t gen_cns_chain(overlap_region_alloc* ol, overlap_region *z, Candidates_list *cl, asg64_v* iidx, int64_t max_lgap, double sgap_rate, int64_t need_filter_khit)
{
    int64_t k, wn = z->w_list.n, aln_n, qs, qe, ts, te, ci, id, rcn = cl->length, kn; k_mer_hit *ka; 
    if(wn <= 0) return 0;
    qs = qe = ts = te = -1; ci = z->shared_seed; id = cl->list[ci].readID;
    for (k = aln_n = 0; k < wn; k++) {
        if(z->w_list.a[k].extra_end < 0 && z->w_list.a[k].y_end != -1) {///anchor
            if(qs == -1) {
                qs = z->w_list.a[k].x_start; ts = z->w_list.a[k-1].y_end+1;
            }
            qe = z->w_list.a[k].x_end+1; te = z->w_list.a[k].y_end+1;
        } else {
            if(qs != -1) {
                // if(z->y_id == 109111 || z->y_id == 75436) {
                //     fprintf(stderr, "\n[M::%s::] q::[%ld, %ld), t::[%ld, %ld), ci::%ld, cn::%lld\n", 
                //         __func__, qs, qe, ts, te, ci, cl->length);
                // }
                ci = push_adp_k_hits(cl, rcn, qs, qe, ts, te, id, ci, z->y_id);
                aln_n++;//[qs, qe); [ts, te)
            }
            qs = qe = ts = te = -1;
        }
    }
    if(qs != -1) {
        // if(z->y_id == 109111 || z->y_id == 75436) {
        //     fprintf(stderr, "\n[M::%s::] q::[%ld, %ld), t::[%ld, %ld), ci::%ld, cn::%lld\n", 
        //         __func__, qs, qe, ts, te, ci, cl->length);
        // }
        ci = push_adp_k_hits(cl, rcn, qs, qe, ts, te, id, ci, z->y_id);
        aln_n++;//[qs, qe); [ts, te)
    }
    // if(z->y_id == 109111 || z->y_id == 75436) {
    //     fprintf(stderr, "\n[M::%s::] q::[%ld, %ld), t::[%ld, %ld), ci::%ld, cn::%lld, rcn::%ld\n", 
    //                     __func__, qs, qe, ts, te, ci, cl->length, rcn);
    // }
    push_adp_k_hits(cl, rcn, (uint64_t)-1, (uint64_t)-1, (uint64_t)-1, (uint64_t)-1, id, ci, z->y_id);

    ka = cl->list+rcn; kn = cl->length-rcn;
    // if(z->y_id == 66) {
    //     fprintf(stderr, "[M::%s::] kn::%ld, q_pos::%u, t_pos::%u\n", __func__, kn, ka[kn-1].self_offset, ka[kn-1].offset);
    // }
    // if(z->y_id == 109111 || z->y_id == 75436) {
    //     prt_khit(cl, ol, NULL, 109111, "sa");
    // }
    if(iidx) gen_weight_khits(iidx, ka, kn);
    // if(z->y_id == 109111 || z->y_id == 75436) {
    //     prt_khit(cl, ol, NULL, 109111, "sb");
    // }
    if(need_filter_khit) {
        kn = lchain_dp_trace(ka, kn, max_lgap, sgap_rate, SGAP); cl->length = rcn + kn;
    }
    cl->length = rcn;
    return kn;
}


uint64_t cigar_gen_cns(overlap_region *z, Candidates_list *cl, ul_ov_t *ov, uint64_t on, uint64_t qn, uint64_t wl,  
const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, UC_Read *tu, bit_extz_t *exz, double e_rate, 
int64_t ql, uint64_t rid, asg64_v* iidx, ul_ov_t *des)
{
    if(on <= 0 || iidx->n <= 0) return 0;
    uint64_t i; 
    for (i = 0; i < on; i++) {
        assert((i<=0)||(ov[i].qs > ov[i-1].qe));
        sub_ciagar_gen(z, cl, ov[i].qs, ov[i].qe, wl, uref, hpc_g, rref, qstr, tu, exz, e_rate, ql, rid);
    }
    return on;
}

///[ys, ye)
uint64_t inline get_win_aln(overlap_region *z, uint64_t wid, int64_t *ys, int64_t *ye, int64_t *err)
{
	(*err) = -2;
	if((wid > 0) && (z->w_list.a[wid].y_end != -1) && (z->w_list.a[wid-1].y_end != -1) && (z->w_list.a[wid].y_end > z->w_list.a[wid-1].y_end)) {
		(*ys) = z->w_list.a[wid-1].y_end+1; 
		(*ye) = z->w_list.a[wid].y_end+1; 
		(*err) = z->w_list.a[wid].error;
		return 1;
	}
	return 0;
}

char* retrive_str_piece_exz(All_reads *rref, const ul_idx_t *uref, char *buf, int64_t s, int64_t l, int64_t rev, int64_t id)
{
	if(rref) recover_UC_Read_sub_region(buf, s, l, rev, rref, id); 
	else if(uref) retrieve_u_seq(NULL, buf, &(uref->ug->u.a[id]), rev, s, l, NULL); 
	else return NULL;
	return buf;
}

uint64_t k_hits_bcheck(All_reads *rref, const ul_idx_t *uref, overlap_region_alloc* ol, Candidates_list *cl, 
uint64_t khit, uint64_t *a, uint64_t a_n, char* qstr, char *str0, char *str1)
{
    if(a_n < 2) return 1;
    uint64_t k, e; char *ref0, *ref1; overlap_region *z;

    e = cl->list[(uint32_t)a[0]].offset; 
    if(e >= khit) e-=khit;
    else return 0;
    z = &(ol->list[cl->list[(uint32_t)a[0]].readID]);
    ref0 = retrive_str_piece_exz(rref, uref, str0, e, khit, z->y_pos_strand, z->y_id);
    // fprintf(stderr, "\n[M::%s::] qstr::%.*s\n", __func__, 
    // (int32_t)khit, qstr+cl->list[(uint32_t)a[0]].self_offset-khit);
    // fprintf(stderr, "[M::%s::qoff->%u::toff->%u::%c] tstr0::%.*s\n", __func__, 
    // cl->list[(uint32_t)a[0]].self_offset, cl->list[(uint32_t)a[0]].offset, "+-"[z->y_pos_strand], (int32_t)khit, ref0);

    for (k = 1; k < a_n; k++) {
        e = cl->list[(uint32_t)a[k]].offset; 
        if(e >= khit) e-=khit;
        else return 0;
        z = &(ol->list[cl->list[(uint32_t)a[k]].readID]);
        ref1 = retrive_str_piece_exz(rref, uref, str1, e, khit, z->y_pos_strand, z->y_id);
        // fprintf(stderr, "[M::%s::qoff->%u::toff->%u::%c] tstr1::%.*s\n", __func__, 
        // cl->list[(uint32_t)a[k]].self_offset, cl->list[(uint32_t)a[k]].offset, "+-"[z->y_pos_strand], (int32_t)khit, ref1);
        if(memcmp(ref0, ref1, khit)) return 0;
    }
    return 1;
}

inline int64_t khit_long_gap(k_mer_hit *a, k_mer_hit *b, double small_bw_rate, int64_t min_small_bw)
{
    int64_t dq, dr, dd, dm;
    dq = b->self_offset-a->self_offset;
    dr = b->offset-a->offset;
	dd = dq>=dr? ((dq)-(dr)): ((dr)-(dq));

    dm = dq>=dr?dr:dq;
    if((dd > (dm*small_bw_rate)) && (dd > min_small_bw)) return 0;
    return 1;
}

int64_t filter_bad_khits(k_mer_hit *sk, k_mer_hit *ek, k_mer_hit* a, int64_t a_n, double small_bw_rate, int64_t min_small_bw)
{
    int64_t k = 0; k_mer_hit *z; double bw_r; int64_t bw, occ = 0;
    bw_r = small_bw_rate; bw = min_small_bw; 
    if(sk) {
        for (k = 0; k < a_n; k++) {
            z = &(a[k]); 
            if((sk && (!khit_long_gap(sk, z, bw_r, bw))) || 
                                (ek && (!khit_long_gap(z, ek, bw_r, bw)))) {
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
            if((sk && (!khit_long_gap(sk, z, bw_r, bw))) || 
                                (ek && (!khit_long_gap(z, ek, bw_r, bw)))) {
                z->offset = z->self_offset = (uint32_t)-1; occ++;
            } else {
                break;
            }
        }
    }

    if(occ) {
        for (k = occ = 0; k < a_n; k++) {
            if(a[k].offset == (uint32_t)-1) continue;
            a[occ++] = a[k];
        }
        a_n = occ;
    } 
    return a_n;
}

int64_t filter_adp_k_hits(k_mer_hit *ka, int64_t kn, Chain_Data *dp, uint64_t qs, uint64_t qe, uint64_t ts, uint64_t te, uint64_t readID, int64_t ci, int64_t *cmi, k_mer_hit *m,
double small_bw_rate, int64_t min_small_bw)
{
    int64_t k = ci, cmi0 = (*cmi); k_mer_hit *p = NULL, n;
    for (; k < kn && ka[k].readID == readID && ka[k].self_offset < qs; k++) {
        p = &(ka[k]);
        // fprintf(stderr, "[M::%s::] p->q::%ld, p->t::%ld\n", __func__, p->self_offset, p->offset);
        if(p->offset >= ts) continue;
        if((m) && ((m->offset >= p->offset) || (m->self_offset >= p->self_offset))) continue;
        ka[(*cmi)++] = *p;
    }
    n.self_offset = qs; n.offset = ts; p = &n;
    if(qs == (uint64_t)-1 || qe == (uint64_t)-1) p = NULL;
    // fprintf(stderr, "[M::%s::] cmi0::%ld, cmi::%ld\n", __func__, cmi0, (*cmi)); back = (*cmi);
    if((*cmi) > cmi0) {
        (*cmi) = cmi0 + filter_bad_khits(m, p, ka+cmi0, (*cmi)-cmi0, small_bw_rate, min_small_bw);
        if((*cmi) > cmi0) {
            (*cmi) = cmi0 + lchain_refine(ka+cmi0, (*cmi)-cmi0, ka+cmi0, dp, 50, 5000, 512, 16); 
        }
    }
    // if(back != (*cmi)) {
    //     fprintf(stderr, "sbsbsbsb[M::%s::] cmi0::%ld, cmi::%ld\n", __func__, cmi0, (*cmi));
    // }
    
    for (; k < kn && ka[k].readID == readID && ka[k].self_offset < qe; k++);
    return k;
}


void refine_khits(overlap_region *z, Candidates_list *cl, Chain_Data *dp, double sgap_rate)
{
    int64_t k, wn = z->w_list.n, aln_n, qs, qe, ts, te, ci, cmi, id, mm_gap = 64; 
    k_mer_hit p; p.self_offset = (uint32_t)-1; p.offset = (uint32_t)-1;
    if(wn <= 0) return;
    qs = qe = ts = te = -1; ci = cmi = z->shared_seed; id = cl->list[ci].readID;
    for (k = aln_n = 0; k < wn; k++) {
        if(z->w_list.a[k].extra_end < 0 && z->w_list.a[k].y_end != -1) {///anchor
            if(qs == -1) {
                qs = z->w_list.a[k].x_start; ts = z->w_list.a[k-1].y_end+1;
            }
            qe = z->w_list.a[k].x_end+1; te = z->w_list.a[k].y_end+1;
        } else {
            if(qs != -1) {
                // fprintf(stderr, "\n[M::%s::] q::[%ld, %ld), t::[%ld, %ld), ci::%ld\n", 
                // __func__, qs, qe, ts, te, ci);
                ci = filter_adp_k_hits(cl->list, cl->length, dp, qs, qe, ts, te, id, ci, &cmi, 
                                                                aln_n?&p:NULL, sgap_rate, mm_gap);                
                p.self_offset = qe - 1; p.offset = te -1;
                aln_n++;//[qs, qe); [ts, te)
            }
            qs = qe = ts = te = -1;
        }
    }
    if(qs != -1) {
        // fprintf(stderr, "\n[M::%s::] q::[%ld, %ld), t::[%ld, %ld), ci::%ld\n", 
        //         __func__, qs, qe, ts, te, ci);
        ci = filter_adp_k_hits(cl->list, cl->length, dp, qs, qe, ts, te, id, ci, &cmi, 
                                                                aln_n?&p:NULL, sgap_rate, mm_gap); 
        p.self_offset = qe - 1; p.offset = te -1;
        aln_n++;//[qs, qe); [ts, te)
    }
    // fprintf(stderr, "\n[M::%s::] q::[%ld, %ld), t::[%ld, %ld), ci::%ld\n", 
    //             __func__, qs, qe, ts, te, ci);
    ci = filter_adp_k_hits(cl->list, cl->length, dp, (uint64_t)-1, (uint64_t)-1, (uint64_t)-1, (uint64_t)-1, id, ci, &cmi, 
                                                            aln_n?&p:NULL, sgap_rate, mm_gap);
    for (k = cmi; k < ci; k++) cl->list[k].readID = ((uint32_t)(0x7fffffff));
}


void count_k_hits_filter(overlap_region_alloc* ol, Candidates_list *cl, asg64_v* ii, Chain_Data *dp)
{
    int64_t k, on = ol->length, m = 0, i, cn = cl->length, scn; k_mer_hit *ca;
    overlap_region *z; uint64_t pid; 
    kv_resize(uint64_t, *ii, ol->length);
    for (k = 0, ii->n = 0; k < on; k++, ii->n++) {
        z = &(ol->list[k]);  
        ii->a[ii->n] = z->shared_seed;
        ii->a[ii->n] <<= 32; ii->a[ii->n] |= k;
        i = z->shared_seed; pid = cl->list[i].readID;
        for (; i < cn && cl->list[i].readID == pid && cl->list[i].readID != ((uint32_t)(0x7fffffff)); i++);
        scn = i - z->shared_seed; ca = cl->list+z->shared_seed;
        i = lchain_refine(ca, scn, ca, dp, 50, 5000, 512, 16); 
        for (; i < scn; i++) ca[i].readID = ((uint32_t)(0x7fffffff));
    }
    radix_sort_bc64(ii->a, ii->a+ii->n); 

    for (k = m = 0; k < on; k++) {
        z = &(ol->list[(uint32_t)ii->a[k]]); i = z->shared_seed; 
        pid = cl->list[i].readID; z->shared_seed = m;
        for (; i < cn && cl->list[i].readID == pid && cl->list[i].readID != ((uint32_t)(0x7fffffff)); i++) {
            cl->list[m] = cl->list[i]; 
            cl->list[m].readID = (uint32_t)ii->a[k]; 
            cl->list[m].cnt = 0;
            m++;
        }
    }
    cl->length = cn = m; 
}



void count_k_hits_adv(All_reads *rref, const ul_idx_t *uref, char* qstr, UC_Read *buf, 
overlap_region_alloc* ol, Candidates_list *cl, asg64_v* ii, Chain_Data *dp, double sgap_rate, 
uint64_t khit, uint64_t basec)
{
    int64_t k, l, on = ol->length, m = 0, i, cn = cl->length, srt_n; 
    overlap_region *z; uint64_t t, ff, pid, *srt; char *str0 = NULL, *str1 = NULL;
    kv_resize(uint64_t, *ii, ol->length);
    for (k = 0, ii->n = 0; k < on; k++, ii->n++) {
        ii->a[ii->n] = ol->list[k].shared_seed;
        ii->a[ii->n] <<= 32; ii->a[ii->n] |= k;
        refine_khits(&(ol->list[k]), cl, dp, sgap_rate);
    }
    radix_sort_bc64(ii->a, ii->a+ii->n); 
    // prt_khit(cl, ol, NULL, 109111, "a");
    for (k = m = 0; k < on; k++) {
        z = &(ol->list[(uint32_t)ii->a[k]]); i = z->shared_seed; 
        pid = cl->list[i].readID; z->shared_seed = m;
        for (; i < cn && cl->list[i].readID == pid && cl->list[i].readID != ((uint32_t)(0x7fffffff)); i++) {
            cl->list[m] = cl->list[i]; 
            cl->list[m].readID = (uint32_t)ii->a[k]; 
            cl->list[m].cnt = 0;
            t = cl->list[m].self_offset; t <<= 32; t |= m;
            kv_push(uint64_t, (*ii), t);
            // if(z->y_id == 109111) {
            //     fprintf(stderr, "[M::%s::m->%ld] qoff::%u, toff::%u\n", __func__, m, 
            //                         cl->list[m].self_offset, cl->list[m].offset);
            // }
            m++;
        }
    }
    cl->length = cn = m; srt = ii->a+on; srt_n = ii->n-on;
    // prt_khit(cl, ol, NULL, 109111, "b");
    radix_sort_bc64(srt, srt+srt_n);
    if(basec) {
        resize_UC_Read(buf, (khit<<1)); str0 = buf->seq; str1 = buf->seq + khit;
    }
    // prt_khit(cl, ol, NULL, 109111, "c");
    for (k = 1, l = 0; k <= srt_n; k++) {
        if(k == cn || (srt[l]>>32) != (srt[k]>>32)) {
            ff = k - l;
            if(basec && ff > 1) {
                if(!k_hits_bcheck(rref, uref, ol, cl, khit, srt+l, k-l, qstr, str0, str1)) ff = 1;
            }
            for (i = l; i < k; i++) {
                cl->list[(uint32_t)srt[i]].cnt = ff;
                // fprintf(stderr, "[M::%s::] pos::%u, cnt::%lu\n", __func__, 
                // cl->list[(uint32_t)srt[i]].self_offset, ff);
            }
            
            l = k;
        }
    }
    // prt_khit(cl, ol, NULL, 109111, "d");
}


void count_k_hits(All_reads *rref, const ul_idx_t *uref, char* qstr, UC_Read *buf, 
overlap_region_alloc* ol, Candidates_list *cl, asg64_v* ii, uint64_t khit, uint64_t basec)
{
    int64_t k, l, on = ol->length, m = 0, i, cn = cl->length, srt_n; 
    overlap_region *z; uint64_t t, ff, pid, *srt; char *str0 = NULL, *str1 = NULL;
    kv_resize(uint64_t, *ii, ol->length);
    for (k = 0, ii->n = 0; k < on; k++, ii->n++) {
        ii->a[ii->n] = ol->list[k].shared_seed;
        ii->a[ii->n] <<= 32; ii->a[ii->n] |= k;
    }
    radix_sort_bc64(ii->a, ii->a+ii->n); 

    for (k = m = 0; k < on; k++) {
        z = &(ol->list[(uint32_t)ii->a[k]]); i = z->shared_seed; 
        pid = cl->list[i].readID; z->shared_seed = m;
        for (; i < cn && cl->list[i].readID == pid; i++) {
            cl->list[m] = cl->list[i]; 
            cl->list[m].readID = (uint32_t)ii->a[k]; 
            cl->list[m].cnt = 0;
            t = cl->list[m].self_offset; t <<= 32; t |= m;
            kv_push(uint64_t, (*ii), t);
            m++;
        }
    }
    cl->length = cn = m; srt = ii->a+on; srt_n = ii->n-on;
    radix_sort_bc64(srt, srt+srt_n);
    if(basec) {
        resize_UC_Read(buf, (khit<<1)); str0 = buf->seq; str1 = buf->seq + khit;
    }

    for (k = 1, l = 0; k <= srt_n; k++) {
        if(k == cn || (srt[l]>>32) != (srt[k]>>32)) {
            ff = k - l;
            if(basec && ff > 1) {
                if(!k_hits_bcheck(rref, uref, ol, cl, khit, srt+l, k-l, qstr, str0, str1)) ff = 1;
            }
            for (i = l; i < k; i++) {
                cl->list[(uint32_t)srt[i]].cnt = ff;
                // fprintf(stderr, "[M::%s::] pos::%u, cnt::%lu\n", __func__, 
                // cl->list[(uint32_t)srt[i]].self_offset, ff);
            }
            
            l = k;
        }
    }

    
}


void tuning_ext_offset(overlap_region *z, k_mer_hit *ch_a, int64_t ch_n, int64_t ql, int64_t tl, int64_t wl, 
int64_t fusion_k_len, int64_t fusion_win_occ, int64_t *ch_s0, int64_t *ch_e0, int64_t *qs0, int64_t *qe0, int64_t *ts0, int64_t *te0, int64_t *mode)
{
    ///[qs, qe)
    if((*mode) != 1 && (*mode) != 2) return;
    int64_t qs = *qs0, qe = *qe0, ts = *ts0, te = *te0, ch_s = *ch_s0, ch_e = *ch_e0;
    int64_t ke, ks, p[3], we, ws, k, wid;
    int64_t sl, gq, gt, gl, gg, update, min_g = INT32_MAX, min_id = -1;
    adjust_ext_offset(&qs, &qe, &ts, &te, ql, tl, 0, *mode); p[0] = p[1] = p[2] = -1;
    if((*mode) == 1) {///forward extension
        //qs0 and ts0 are fixed; te0 = -1, qe0 is unreliable
        if(qe > (*qe0)) {
            ///find k-mer hit
            ke = ch_a[ch_s].self_offset;
            for (ch_e = ch_s; (ch_e < ch_n) && (ch_a[ch_e].self_offset < (*qe0)); ch_e++) {
                ke = ch_a[ch_e].self_offset;
            }
            for (ke += fusion_k_len; (ch_e < ch_n) && (ch_a[ch_e].self_offset < ke); ch_e++) {
                if(ch_a[ch_e].self_offset<(*qe0)) continue;
                if((ch_a[ch_e].offset<=(*ts0))||(ch_a[ch_e].self_offset<=(*qs0))) continue;//not co-linear
                if(is_pri_aln(ch_a[ch_e])) {
                    p[2] = ch_e; break;
                } else if(ch_a[ch_e].cnt > 1 && p[1] == -1) {
                    p[1] = ch_e;
                } else if(p[0] == -1) {
                    p[0] = ch_e;
                }
            }
            if(p[2] != -1) p[0] = p[2];
            else if(p[1] != -1) p[0] = p[1];
            if(p[0] != -1) {
                (*mode) = 0;///global
                (*qe0) = ch_a[p[0]].self_offset; (*te0) = ch_a[p[0]].offset; (*ch_e0) = p[0];
                return;
            }

            ///find aligned window
            we = (*qe0); we/=wl; we *= wl; we +=wl; we--; ///next window
            for (k = 0; we < qe && we <= z->x_pos_e && k < fusion_win_occ; we+=wl) {//[ws, we]; [qs, qe)
                wid = get_win_id_by_e(z, we, wl, NULL);
                if(z->w_list.a[wid].y_end == -1) continue;//unmapped
                if(z->w_list.a[wid].x_end < (*qs0)) continue;
                if(z->w_list.a[wid].y_end < (*ts0)) continue;
                sl = z->w_list.a[wid].x_end+1-z->w_list.a[wid].x_start;
                if(z->w_list.a[wid].error > (sl/A_L)) continue;
                gq = z->w_list.a[wid].x_end-(*qs0); 
                gt = z->w_list.a[wid].y_end-(*ts0); 
                gl = MIN(gq, gt); gl/=16; if(gl <= 0) gl = 1;
                gg = (gq>=gt)?(gq-gt):(gt-gq); gg /= gl;
                update = 0;
                if(min_g>gg) {
                    update = 1;
                } else if((min_g==gg)&&(z->w_list.a[min_id].error>z->w_list.a[wid].error)) {
                    update = 1;
                } 
                if(update) {
                    min_g = gg; min_id = wid;
                }
                k++;
            }
        }

        if(min_id != -1) {
            (*qe0) = z->w_list.a[min_id].x_end+1; 
            (*te0) = z->w_list.a[min_id].y_end+1; 
            (*mode) = 0;///global
        } else {
            (*qe0) = qe; (*te0) = te;//extension
        }
    } else if((*mode) == 2) {///backward extension
        //qe0 and te0 are fixed; ts0 = -1, qs0 is unreliable
        if(qs < (*qs0)) {
            ///find k-mer hit
            ks = ch_a[ch_e].self_offset;
            for (ch_s = ch_e; (ch_s >= 0) && (ch_a[ch_s].self_offset > (*qs0)); ch_s--) {
                ks = ch_a[ch_s].self_offset;
            }
            for (ks -= fusion_k_len; (ch_s >= 0) && (ch_a[ch_s].self_offset > ks); ch_s--) {
                if(ch_a[ch_s].self_offset>(*qs0)) continue;
                if((ch_a[ch_s].offset>=(*te0))||(ch_a[ch_s].self_offset>=(*qe0))) continue;//not co-linear
                if(is_pri_aln(ch_a[ch_s])) {
                    p[2] = ch_s; break;
                } else if(ch_a[ch_s].cnt > 1 && p[1] == -1) {
                    p[1] = ch_s;
                } else if(p[0] == -1) {
                    p[0] = ch_s;
                }
            }
            // fprintf(stderr, "-[M::%s::] utg%.6dl(%c), p[0]::%ld, p[1]::%ld, p[2]::%ld\n", 
            //             __func__, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], p[0], p[1], p[2]);
            if(p[2] != -1) p[0] = p[2];
            else if(p[1] != -1) p[0] = p[1];
            if(p[0] != -1) {
                (*mode) = 0;///global
                (*qs0) = ch_a[p[0]].self_offset; (*ts0) = ch_a[p[0]].offset; (*ch_s0) = p[0];
                return;
            }

            ///find aligned window
            ws = (*qs0)-1; ws/=wl; ws*=wl; 
            for (k = 0; ws >= qs && ws >= z->x_pos_s && k < fusion_win_occ; ws-=wl) {//[ws, we]; [qs, qe)
                wid = get_win_id_by_s(z, ws, wl, NULL);
                if(z->w_list.a[wid].y_end == -1) continue;//unmapped
                if(z->w_list.a[wid].x_start >= (*qe0)) continue;
                if(z->w_list.a[wid].y_start >= (*te0)) continue;
                sl = z->w_list.a[wid].x_end+1-z->w_list.a[wid].x_start;
                if(z->w_list.a[wid].error > (sl/A_L)) continue;
                gq = (*qe0) - z->w_list.a[wid].x_start; 
                gt = (*te0) - z->w_list.a[wid].y_start; 
                gl = MIN(gq, gt); gl/=16; if(gl <= 0) gl = 1;
                gg = (gq>=gt)?(gq-gt):(gt-gq); gg /= gl;
                update = 0;
                if(min_g>gg) {
                    update = 1;
                } else if((min_g==gg)&&(z->w_list.a[min_id].error>z->w_list.a[wid].error)) {
                    update = 1;
                } 
                if(update) {
                    min_g = gg; min_id = wid;
                }
                k++;
            }
        }

        if(min_id != -1) {
            (*qs0) = z->w_list.a[min_id].x_start; 
            (*ts0) = z->w_list.a[min_id].y_start; 
            (*mode) = 0;///global
        } else {
            (*qs0) = qs; (*ts0) = ts; ///still extension
        }
    }

    if((*ch_e0) == -1) {
        for ((*ch_e0)=(*ch_s0);((*ch_e0)<ch_n)&&(ch_a[(*ch_e0)].self_offset>=(*qs0))&&(ch_a[(*ch_e0)].self_offset<=(*qe0)); (*ch_e0)++);
    }

    if((*ch_s0) == -1) {
        for ((*ch_s0)=(*ch_e0);((*ch_s0)>=0)&&(ch_a[(*ch_s0)].self_offset>=(*qs0))&&(ch_a[(*ch_s0)].self_offset<=(*qe0)); (*ch_s0)--);
        (*ch_s0)++;
    }
    ///boundry
    for (;((*ch_e0)<ch_n)&&(ch_a[(*ch_e0)].self_offset>=(*qs0))&&(ch_a[(*ch_e0)].self_offset<=(*qe0)); (*ch_e0)++);
}

int64_t cal_estimate_err(overlap_region *z, int64_t wl, int64_t qs, int64_t qe, double e_rate)
{
    int64_t k, ws, we, wid, os, oe, ovlp, tot, cov_l, est = (qe-qs)*e_rate, wn = z->w_list.n; 
    if(!wn) return est;
    if(qs < z->x_pos_s) qs = z->x_pos_s;
    if(qe > z->x_pos_e+1) qe = z->x_pos_e+1;
    ws = qs/wl; ws *= wl; wid = get_win_id_by_s(z, ws, wl, NULL);
    if(wid>=wn) wid = wn-1;
    for (k = wid;k < wn && qs > z->w_list.a[k].x_end; k++); if(k == wn) return est;
    ///qs <= z->w_list.a[k].x_end
    for (;k>=0 && qs < z->w_list.a[k].x_start; k--); if(k < 0) k = 0;
    ///qs >= z->w_list.a[k].x_start

    for (tot = cov_l = 0; k < wn && z->w_list.a[k].x_start < qe; k++) {
        if(z->w_list.a[k].y_end == -1) continue;
        ws = z->w_list.a[k].x_start; we = z->w_list.a[k].x_end+1;
        os = MAX(qs, ws); oe = MIN(qe, we);
	    ovlp = ((oe>os)? (oe-os):0);
        if(!ovlp) continue;
        cov_l += ovlp;
        // if(!ovlp) {
        //     fprintf(stderr, "\n[M::%s::] utg%.6dl(%c), q::[%ld, %ld), w::[%ld, %ld), z::::[%d, %d)\n", 
        //         __func__, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], qs, qe, ws, we, z->x_pos_s, z->x_pos_e+1);
        // }
        if(ovlp == (we-ws)) {
            tot += z->w_list.a[k].error;
        } else {
            tot += ((double)z->w_list.a[k].error)*((double)ovlp)/((double)(we-ws));
        }
    }
    tot += ((qe-qs)-cov_l)*e_rate;
    return tot;
}

#define UC_Read_resize(v, s) do {\
		if ((v).size<(s)) {REALLOC((v).seq,(s));(v).size=(s);}\
	} while (0)

///[s, e); [ps, pe)
inline char *retrieve_str_seq_exz(UC_Read *tu, int64_t s, int64_t l, 
int64_t ps, int64_t pl, uint8_t rev, const ul_idx_t *uref, hpc_t *hpc_g, 
All_reads *rref, int64_t id)
{
    if(!hpc_g) {
        char *str; int64_t ss = s, sl = l; tu->length = l;
        UC_Read_resize(*tu, sl); str = tu->seq; 
        if(s == ps) {
            if(l <= pl) return tu->seq;
            str = tu->seq + pl; ss = ps + pl; sl = l - pl;
        }
        if(uref) {
            retrieve_u_seq(NULL, str, &(uref->ug->u.a[id]), rev, ss, sl, NULL);
        } else if(rref) {
            recover_UC_Read_sub_region(str, ss, sl, rev, rref, id);
        }
        return tu->seq;
    } else {
        return hpc_str(*hpc_g, id, rev) + s;
    }
}

void cal_exz_global(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{
    int32_t bd, nword;
    bd = (((thre)<<1)+1); nword = ((bd>>bitw)+(!!(bd&bitz)));

    if(nword <= 1) {
        ed_band_cal_global_64_w_trace(pstr, pn, tstr, tn, thre, ez);
    } else if(nword == 2) {
        ed_band_cal_global_128_w_trace(pstr, pn, tstr, tn, thre, ez);
    } else {
        ed_band_cal_global_infi_w_trace(pstr, pn, tstr, tn, thre, &nword, ez);
    }
}

void cal_exz_global_simi(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{
    int32_t bd, nword;
    bd = (((thre)<<1)+1); nword = ((bd>>bitw)+(!!(bd&bitz)));

    if(nword <= 1) {
        ed_band_cal_global_64_w(pstr, pn, tstr, tn, thre, ez);
    } else if(nword == 2) {
        ed_band_cal_global_128_w(pstr, pn, tstr, tn, thre, ez);
    } else {
        ed_band_cal_global_infi_w(pstr, pn, tstr, tn, thre, &nword, ez);
    }
}

void cal_exz_extension_0(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{
    int32_t bd, nword;
    bd = (((thre)<<1)+1); nword = ((bd>>bitw)+(!!(bd&bitz)));

    if(nword <= 1) {
        ed_band_cal_extension_64_0_w_trace(pstr, pn, tstr, tn, thre, ez);
    } else if(nword == 2) {
        ed_band_cal_extension_128_0_w_trace(pstr, pn, tstr, tn, thre, ez);
    } else {
        ed_band_cal_extension_infi_0_w_trace(pstr, pn, tstr, tn, thre, &nword, ez);
    }
}

void cal_exz_extension_0_simi(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{
    int32_t bd, nword;
    bd = (((thre)<<1)+1); nword = ((bd>>bitw)+(!!(bd&bitz)));

    if(nword <= 1) {
        ed_band_cal_extension_64_0_w(pstr, pn, tstr, tn, thre, ez);
    } else if(nword == 2) {
        ed_band_cal_extension_128_0_w(pstr, pn, tstr, tn, thre, ez);
    } else {
        ed_band_cal_extension_infi_0_w(pstr, pn, tstr, tn, thre, &nword, ez);
    }
}

void cal_exz_extension_1(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{
    int32_t bd, nword;
    bd = (((thre)<<1)+1); nword = ((bd>>bitw)+(!!(bd&bitz)));

    if(nword <= 1) {
        ed_band_cal_extension_64_1_w_trace(pstr, pn, tstr, tn, thre, ez);
    } else if(nword == 2) {
        ed_band_cal_extension_128_1_w_trace(pstr, pn, tstr, tn, thre, ez);
    } else {
        ed_band_cal_extension_infi_1_w_trace(pstr, pn, tstr, tn, thre, &nword, ez);
    }
}

void cal_exz_extension_1_simi(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{
    int32_t bd, nword;
    bd = (((thre)<<1)+1); nword = ((bd>>bitw)+(!!(bd&bitz)));

    if(nword <= 1) {
        ed_band_cal_extension_64_1_w(pstr, pn, tstr, tn, thre, ez);
    } else if(nword == 2) {
        ed_band_cal_extension_128_1_w(pstr, pn, tstr, tn, thre, ez);
    } else {
        ed_band_cal_extension_infi_1_w(pstr, pn, tstr, tn, thre, &nword, ez);
    }
}

void cal_exz_semi(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, int32_t aux_beg, bit_extz_t *ez)
{
    int32_t bd, nword;
    bd = (((thre)<<1)+1); nword = ((bd>>bitw)+(!!(bd&bitz)));

    if(nword <= 1) {
        ed_band_cal_semi_64_w_absent_diag_trace(pstr, pn, tstr, tn, thre, aux_beg, ez);
    } else if(nword == 2) {
        ed_band_cal_semi_128_w_absent_diag_trace(pstr, pn, tstr, tn, thre, aux_beg, ez);
    } else {
        ed_band_cal_semi_infi_w_absent_diag_trace(pstr, pn, tstr, tn, thre, aux_beg, &nword, ez);
    }
}

void cal_exz_semi_simi(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, int32_t aux_beg, bit_extz_t *ez)
{
    int32_t bd, nword;
    bd = (((thre)<<1)+1); nword = ((bd>>bitw)+(!!(bd&bitz)));

    if(nword <= 1) {
        ed_band_cal_semi_64_w_absent_diag(pstr, pn, tstr, tn, thre, aux_beg, ez);
    } else if(nword == 2) {
        ed_band_cal_semi_128_w_absent_diag(pstr, pn, tstr, tn, thre, aux_beg, ez);
    } else {
        ed_band_cal_semi_infi_w_absent_diag(pstr, pn, tstr, tn, thre, aux_beg, &nword, ez);
    }
}

void ref_cigar_check(char* qstr, UC_Read *tu, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, int64_t id, int64_t rev, bit_extz_t *ez)
{
    int64_t pts = -1, pte = -1, tl = ez->pe-ez->ps+1, ql = ez->te-ez->ts+1, bps, bts, k; 
    char *q, *t; bps = ez->ps; bts = ez->ts;

    q = qstr + ez->ts; 
    t = retrieve_str_seq_exz(tu, ez->ps, tl, pts, pte-pts, rev, uref, hpc_g, rref, id);
    ez->ps = ez->ts = 0; 
    if(!cigar_check(t, q, ez)){
        fprintf(stderr, "[M::%s::] cigar_n::%d\n", __func__, (int32_t)ez->cigar.n);
        for (k = 0; k < (int32_t)ez->cigar.n; k++) {
            fprintf(stderr, "[M::%s::%ld] cigar_len::%u, c::%u\n", __func__, k,
            ez->cigar.a[k]&(0x3fff), (ez->cigar.a[k]>>14));
        }
        
        fprintf(stderr, "[M::%s::l->%ld] s::%ld, pstr::%.*s\n", __func__, tl, bps, (int32_t)tl, t);
        fprintf(stderr, "[M::%s::l->%ld] s::%ld, tstr::%.*s\n", __func__, ql, bts, (int32_t)ql, q);
        exit(1);
    }
    ez->ps = bps; ez->ts = bts;
}

int64_t cal_exz_infi_adv(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, 
bit_extz_t *exz, char* qstr, UC_Read *tu, int64_t qs, int64_t qe, int64_t ts, int64_t te, 
int64_t *pts, int64_t *pte, int64_t thre, int64_t *pthre, int64_t q_tot_l, int64_t mode)
{
    clear_align(*exz);
    int64_t aux_beg = 0, ql, tl, t_tot_l = -1, dd; 
    char *q_string, *t_string; int32_t rev = z->y_pos_strand, id = z->y_id; 
    ql = qe - qs; tl = te - ts; dd = MAX(ql, tl); 
    if(hpc_g) t_tot_l = hpc_len(*hpc_g, id);
    else if(uref) t_tot_l = uref->ug->u.a[id].len;
    else t_tot_l = Get_READ_LENGTH((*rref), id);
   
    if(mode == 3) {
        update_semi_coord(uref, hpc_g, rref, z, qs, qe, ((thre>dd)?dd:thre), &ts, &te, &aux_beg);
    } else if(mode == 1 || mode == 2) {
        adjust_ext_offset(&qs, &qe, &ts, &te, q_tot_l, t_tot_l, ((thre>dd)?dd:thre), mode);
    }
    
    if((qe > qs) && (te > ts) && (ts != -1) && (te != -1)) {
        ql = qe - qs; tl = te - ts; 
        dd = MAX(ql, tl); 
        if(thre > dd) thre = dd; 
        if(thre <= (*pthre)) return 0;
        (*pthre) = thre;
        
        q_string = qstr + qs; 
        t_string = retrieve_str_seq_exz(tu, ts, tl, (*pts), (*pte)-(*pts), rev, uref, hpc_g, rref, id);
        (*pts) = ts; (*pte) = te;

        if(mode == 0) { //global
            cal_exz_global(t_string, tl, q_string, ql, thre, exz);
        } else if(mode == 1) {///forward extension
            cal_exz_extension_0(t_string, tl, q_string, ql, thre, exz);
        } else if(mode == 2) {///backward extension
            cal_exz_extension_1(t_string, tl, q_string, ql, thre, exz);
        } else if(mode == 3) {//semi-global
            cal_exz_semi(t_string, tl, q_string, ql, thre, aux_beg, exz);
        }

        if(is_align(*exz)) {
            // cigar_check(t_string, q_string, exz); 
            // if(mode == 1) {
            //     fprintf(stderr, "\n[M::%s::ql::%ld] qs::%ld, qe::%ld, ts::%ld, te::%ld, mode::%ld, err::%d, thre::%d, exz_q[%d, %d], exz_t[%d, %d]\n", 
            //                      __func__, ql, qs, qe, ts, te, mode, exz->err, exz->thre, exz->ts, exz->te, exz->ps, exz->pe);
            //     fprintf(stderr, "[M::%s::] pstr::%.*s\n", __func__, (int32_t)tl, t_string);
            //     fprintf(stderr, "[M::%s::] tstr::%.*s\n", __func__, (int32_t)ql, q_string);
            //     fprintf(stderr, "[M::%s::] exz->cigar.n::%d\n", __func__, (int32_t)exz->cigar.n);
            // }
            exz->ps += ts; exz->pe += ts;
            exz->ts += qs; exz->te += qs;
            return 1;
        }
        return 0;
    }
    return 0;
}


int64_t cal_exact_exz(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, 
bit_extz_t *exz, char* qstr, UC_Read *tu, int64_t qs, int64_t qe, int64_t ts, int64_t te, 
int64_t *pts, int64_t *pte, int64_t q_tot_l, int64_t mode)
{
    clear_align(*exz); exz->thre = 0; exz->cigar.n = 0;
    int64_t ql, tl, t_tot_l = -1; 
    char *q_string, *t_string; int32_t rev = z->y_pos_strand, id = z->y_id; ql = qe - qs;
    if(hpc_g) t_tot_l = hpc_len(*hpc_g, id);
    else if(uref) t_tot_l = uref->ug->u.a[id].len;
    else t_tot_l = Get_READ_LENGTH((*rref), id);

    if(mode == 3) {//semi
        ts = (qs - z->x_pos_s) + z->y_pos_s; ts += y_start_offset(qs, &(z->f_cigar));
        te = ts + ql; 
    } else if(mode == 1) {///forward extension
        te = ts + ql; 
    } else if(mode == 2) {///backward extension
        ts = te - ql; 
    }
    if(ts < 0) ts = 0;
    if(ts > t_tot_l) ts = t_tot_l;
    if(te > t_tot_l) te = t_tot_l;
    ql = qe - qs; tl = te - ts; 
    if(ql != tl) return 0;

    q_string = qstr + qs; 
    t_string = retrieve_str_seq_exz(tu, ts, tl, (*pts), (*pte)-(*pts), rev, uref, hpc_g, rref, id);
    (*pts) = ts; (*pte) = te;

    if(memcmp(q_string, t_string, ql)) return 0;
    exz->err = 0; push_trace(&(exz->cigar), 0, ql);
    exz->pl = tl; exz->ps = 0; exz->pe = tl-1;
    exz->tl = ql; exz->ts = 0; exz->te = ql-1;
    // cigar_check(t_string, q_string, exz);
    // if(!cigar_check(t_string, q_string, exz)) {
    //     fprintf(stderr, "[M::%s::] cigar_n::%d\n", __func__, (int32_t)exz->cigar.n);
    //     fprintf(stderr, "[M::%s::] pstr::%.*s\n", __func__, (int32_t)tl, t_string);
    //     fprintf(stderr, "[M::%s::] tstr::%.*s\n", __func__, (int32_t)ql, q_string);
    // }
    exz->ps += ts; exz->pe += ts;
    exz->ts += qs; exz->te += qs;
    return 1;
}

//[qmin, qmax) && [tmin, tmax)
void adjust_ext_offset_fixed_t(int64_t *qs, int64_t *qe, int64_t *ts, int64_t *te, 
int64_t qmin, int64_t qmax, int64_t tmin, int64_t tmax, int64_t thre, int64_t mode)
{
    int64_t qoff, toff;
    if(mode == 1) {///forward extension
        qoff = qmax - (*qs); toff = tmax - (*ts);
        if(qoff <= toff) {
            (*qe) = qmax; (*te) = (*ts) + qoff + thre;
        } else {
            (*te) = tmax; (*qe) = (*qs) + toff + thre;
        }
    } else if(mode == 2) {///backward extension
        qoff = (*qe) - qmin; toff = (*te) - tmin;
        if(qoff <= toff) {
            (*qs) = qmin; (*ts) = (*te) - qoff - thre;
        } else {
            (*ts) = tmin; (*qs) = (*qe) - toff - thre;
        }
    }
    if((*qs) < qmin) (*qs) = qmin;
    if((*ts) < tmin) (*ts) = tmin;
    if((*qe) > qmax) (*qe) = qmax;
    if((*te) > tmax) (*te) = tmax;
}



int64_t cal_exz_infi_simi_adv(const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, 
bit_extz_t *exz, char* qstr, UC_Read *tu, int64_t qs, int64_t qe, int64_t ts, int64_t te, 
int64_t *pts, int64_t *pte, int64_t thre, int64_t *pthre, int64_t qmin, int64_t qmax, 
int64_t tmin, int64_t tmax, int32_t rev, int32_t id, int64_t mode, overlap_region *z,
int64_t gen_trace)
{   ///mode cannot be 3
    clear_align(*exz);
    int64_t aux_beg = 0, ql, tl, dd; char *q_string, *t_string; 
    ql = qe - qs; tl = te - ts; dd = MAX(ql, tl); 
    if(mode == 3 && z) {
        update_semi_coord(uref, hpc_g, rref, z, qs, qe, ((thre>dd)?dd:thre), &ts, &te, &aux_beg);
    } else if(mode == 1 || mode == 2) {
        adjust_ext_offset_fixed_t(&qs, &qe, &ts, &te, qmin, qmax, tmin, tmax, ((thre>dd)?dd:thre), mode);
    }

    // if(ts == 18327 && te == 18601 && qs == 145990 && qe == 146191) {
    //         fprintf(stderr, "[M::%s::] mode::%ld, q::[%ld, %ld), t::[%ld, %ld), thre::%ld, pthre::%ld\n", 
    //         __func__, mode, qs, qe, ts, te, thre, *pthre);
    // }
    
    if((qe > qs) && (te > ts) && (ts != -1) && (te != -1)) {
        ql = qe - qs; tl = te - ts; 
        dd = MAX(ql, tl); 
        if(thre > dd) thre = dd; 
        // if(ts == 18327 && te == 18601 && qs == 145990 && qe == 146191) {
        //     fprintf(stderr, "[M::%s::] mode::%ld, q::[%ld, %ld), t::[%ld, %ld), thre::%ld, pthre::%ld\n", 
        //     __func__, mode, qs, qe, ts, te, thre, *pthre);
        // }
        if(thre <= (*pthre)) return 0;
        (*pthre) = thre;
        
        q_string = qstr + qs; 
        t_string = retrieve_str_seq_exz(tu, ts, tl, (*pts), (*pte)-(*pts), rev, uref, hpc_g, rref, id);
        (*pts) = ts; (*pte) = te;

        if(!gen_trace) {
            if(mode == 0) { //global
                cal_exz_global_simi(t_string, tl, q_string, ql, thre, exz);
            } else if(mode == 1) {///forward extension
                cal_exz_extension_0_simi(t_string, tl, q_string, ql, thre, exz);
            } else if(mode == 2) {///backward extension
                cal_exz_extension_1_simi(t_string, tl, q_string, ql, thre, exz);
            } else if(mode == 3) {//semi-global; mode cannot be 3
                cal_exz_semi_simi(t_string, tl, q_string, ql, thre, aux_beg, exz);
            }
        } else {
            if(mode == 0) { //global
                cal_exz_global(t_string, tl, q_string, ql, thre, exz);
            } else if(mode == 1) {///forward extension
                cal_exz_extension_0(t_string, tl, q_string, ql, thre, exz);
            } else if(mode == 2) {///backward extension
                cal_exz_extension_1(t_string, tl, q_string, ql, thre, exz);
            } else if(mode == 3) {//semi-global
                cal_exz_semi(t_string, tl, q_string, ql, thre, aux_beg, exz);
            }
        }
        

        if(is_align(*exz)) {
            exz->ps += ts; exz->pe += ts;
            exz->ts += qs; exz->te += qs;
            // if(exz->ps == 18327 && exz->pe + 1 == 18601 && exz->ts == 145990 && exz->te + 1 == 146191) {
            //     fprintf(stderr, "[M::%s::] mode::%ld, q::[%ld, %ld), t::[%ld, %ld)\n", 
            //     __func__, mode, qs, qe, ts, te);
            //     fprintf(stderr, "[M::%s::] q_string::%.*s\n", __func__, (int32_t)ql, q_string);
            //     fprintf(stderr, "[M::%s::] t_string::%.*s\n", __func__, (int32_t)tl, t_string);
            // }
            return 1;
        }
        // if(ts == 18327 && te == 18601 && qs == 145990 && qe == 146191) {
        //     fprintf(stderr, "[M::%s::] mode::%ld, q::[%ld, %ld), t::[%ld, %ld)\n", 
        //     __func__, mode, qs, qe, ts, te);
        //     fprintf(stderr, "[M::%s::] q_string::%.*s\n", __func__, (int32_t)ql, q_string);
        //     fprintf(stderr, "[M::%s::] t_string::%.*s\n", __func__, (int32_t)tl, t_string);
        // }
        return 0;
    }
    return 0;
}

int64_t cal_exact_simi_exz(const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, 
bit_extz_t *exz, char* qstr, UC_Read *tu, int64_t qs, int64_t qe, int64_t ts, int64_t te, 
int64_t *pts, int64_t *pte, int64_t q_tot_l, int64_t t_tot_l, int32_t rev, int32_t id, 
int64_t mode, overlap_region *z)
{
    clear_align(*exz); exz->thre = 0; exz->cigar.n = 0;
    int64_t ql, tl; char *q_string, *t_string; ql = qe - qs;

    if(mode == 3 && z) {//semi
        ts = (qs - z->x_pos_s) + z->y_pos_s; ts += y_start_offset(qs, &(z->f_cigar));
        te = ts + ql; 
    } else if(mode == 1) {///forward extension
        te = ts + ql; 
    } else if(mode == 2) {///backward extension
        ts = te - ql; 
    }
    if(ts < 0) ts = 0;
    if(ts > t_tot_l) ts = t_tot_l;
    if(te > t_tot_l) te = t_tot_l;
    ql = qe - qs; tl = te - ts; 
    if(ql != tl) return 0;

    q_string = qstr + qs; 
    t_string = retrieve_str_seq_exz(tu, ts, tl, (*pts), (*pte)-(*pts), rev, uref, hpc_g, rref, id);
    (*pts) = ts; (*pte) = te;

    if(memcmp(q_string, t_string, ql)) return 0;
    exz->err = 0; push_trace(&(exz->cigar), 0, ql);
    exz->pl = tl; exz->ps = 0; exz->pe = tl-1;
    exz->tl = ql; exz->ts = 0; exz->te = ql-1;
    
    exz->ps += ts; exz->pe += ts;
    exz->ts += qs; exz->te += qs;
    return 1;
}


void append_wcigar(window_list *idx, window_list_alloc *res, bit_extz_t *exz)
{
    // fprintf(stderr, "[M::%s::] idx->cidx::%u, idx->clen::%u, res->c.n_0::%u, exz->cigar.n::%u, ", 
    //         __func__, idx->cidx, idx->clen, (uint32_t)res->c.n, (uint32_t)exz->cigar.n);
    if(idx->clen > 0) {
        // if(!(idx->cidx+idx->clen == res->c.n)) {
        //     fprintf(stderr, "[M::%s::] idx->cidx::%u, idx->clen::%u, res->c.n::%u\n", 
        //     __func__, idx->cidx, idx->clen, (uint32_t)res->c.n);
        // }
        assert(idx->cidx+idx->clen == res->c.n);
        if(exz->cigar.n > 0) {
            uint16_t c0, c; uint32_t l0, l, ci = 0, cn; 
            ///last item of old cigar
            c0 = (res->c.a[res->c.n-1]>>14); l0 = (res->c.a[res->c.n-1]&(0x3fff));
            ///first item of new cigar
            ci = pop_trace(&(exz->cigar), ci, &c, &l);
            if(c0 == c) {l += l0; res->c.n--; idx->clen--;}
            push_trace(((asg16_v *)(&(res->c))), c, l); 
            idx->clen = res->c.n-idx->cidx;

            cn = exz->cigar.n-ci;
            if(cn > 0) {
                kv_resize(uint16_t, res->c, (res->c.n+cn));
                memcpy(res->c.a+res->c.n, exz->cigar.a+ci, cn*sizeof(*(res->c.a)));
                idx->clen += cn; res->c.n += cn;
            }
        }
    } else {
        push_wcigar(idx, res, exz);///if exz is the first item
    }
    // fprintf(stderr, "[M::%s::] res->c.n::%u\n", __func__, (uint32_t)res->c.n);
}

void push_alnw(overlap_region *aux_o, bit_extz_t *exz)
{
    window_list *p = NULL; int64_t t;
    if(aux_o->w_list.n > 0) {
        p = &(aux_o->w_list.a[aux_o->w_list.n-1]);
        // fprintf(stderr, "+[M::%s::wn->%d] px::[%d, %d], py::[%d, %d], pe::%d, exz->t::[%u, %u], exz->p::[%u, %u], exz->e::%d, clen::%u\n", 
        //     __func__, (int32_t)(aux_o->w_list.n), p->x_start, p->x_end, p->y_start, p->y_end, p->error, 
        //     exz->ts, exz->te, exz->ps, exz->pe, exz->err, p->clen);
        // assert((p->x_end<exz->ts)&&(p->y_end<exz->ps));
        
        if(p->clen > 0) {
            t = ((int64_t)p->error) + ((int64_t)exz->err);
            ///note: t cannot be equal to INT16_MAX; otherwise it is unable to distiguish unaligned regions
            if(((p->x_end+1) == exz->ts) && ((p->y_end+1) == exz->ps) && (t < INT16_MAX)) {
                p->x_end = exz->te; p->y_end = exz->pe; p->error += exz->err;
                append_wcigar(p, &(aux_o->w_list), exz);
            //     fprintf(stderr, "-[M::%s::wn->%d] px::[%d, %d], py::[%d, %d], pe::%d, exz->t::[%u, %u], exz->p::[%u, %u], exz->e::%d, clen::%u\n", 
            // __func__, (int32_t)(aux_o->w_list.n), p->x_start, p->x_end, p->y_start, p->y_end, p->error, 
            // exz->ts, exz->te, exz->ps, exz->pe, exz->err, p->clen);
                return;
            }
        }
    }
    // fprintf(stderr, "[M::%s::wn->%d] exz->t::[%u, %u], exz->p::[%u, %u], exz->e::%d\n", 
    //         __func__, (int32_t)(aux_o->w_list.n), exz->ts, exz->te, exz->ps, exz->pe, exz->err);
    kv_pushp(window_list, aux_o->w_list, &p);
    p->x_start = exz->ts; p->x_end = exz->te;
    p->y_start = exz->ps; p->y_end = exz->pe;
    p->error_threshold = 0; p->error = exz->err;///single round of alignment cannot have INT16_MAX errors
    push_wcigar(p, &(aux_o->w_list), exz);
}

///[qs, qe] && [ts, te]
void push_unmap_alnw(overlap_region *aux_o, int32_t qs, int64_t qe, int64_t ts, int64_t te, int64_t mode)
{
    window_list *p = NULL;
    kv_pushp(window_list, aux_o->w_list, &p);
    p->x_start = qs; p->x_end = qe;
    p->y_start = ts; p->y_end = te;
    p->error_threshold = mode; p->error = INT16_MAX;
    p->extra_begin = p->extra_end = -1;
    p->cidx = p->clen = 0;
}

///[qs, qe] && [ts, te]
void push_replace_alnw(overlap_region *aux_o, int32_t qs, int64_t qe, int64_t ts, int64_t te, int64_t mode)
{
    window_list *p = NULL;
    kv_pushp(window_list, aux_o->w_list, &p);
    p->x_start = qs; p->x_end = qe;
    p->y_start = ts; p->y_end = te;
    p->error_threshold = mode; p->error = 0;
    p->extra_begin = p->extra_end = 0;
    p->cidx = p->clen = 0;
}

///[qs, qe] && [ts, te]
void push_replace_alnw_adv(overlap_region *z, int64_t wl, overlap_region *aux_o, int32_t qs, int64_t qe, int64_t ts, int64_t te, int64_t mode)
{
    int64_t k, wsk, wek, err, t; window_list *p = NULL;
    wsk = get_win_id_by_s(z, qs, wl, NULL); assert(z->w_list.a[wsk].x_start == qs);
    wek = get_win_id_by_s(z, qe+1, wl, NULL); assert(z->w_list.a[wek].x_end == (qe+1));
    kv_pushp(window_list, aux_o->w_list, &p);
    p->x_start = p->x_end = qs;
    p->y_start = p->y_end = ts;
    p->error_threshold = mode; p->error = 0;
    p->extra_begin = p->extra_end = 0;
    p->cidx = p->clen = 0; err = 0;

    for (k = wsk; k <= wek; k++) {
        assert(z->w_list.a[k].y_end != -1);
        t = err + z->w_list.a[k].error;
        if(t < INT16_MAX) {///note: t cannot be equal to INT16_MAX; otherwise it is unable to distiguish unaligned regions
            err += z->w_list.a[k].error;
            p->x_end = z->w_list.a[k].x_end;
            p->y_end = z->w_list.a[k].y_end;
        } else {
            err = z->w_list.a[k].error;
            kv_pushp(window_list, aux_o->w_list, &p);
            p->x_start = z->w_list.a[k].x_start;
            p->x_end = z->w_list.a[k].x_end;
            p->y_start = z->w_list.a[k].y_start;
            p->y_end = z->w_list.a[k].y_end;
            p->error_threshold = mode;
            p->extra_begin = p->extra_end = 0;
            p->cidx = p->clen = 0;
        }
        p->error = err; 
    }
    p->x_end = qe; p->y_end = te;
}



int64_t hc_aln_exz_adv(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, 
char* qstr, UC_Read *tu, int64_t qs, int64_t qe, int64_t ts, int64_t te, int64_t mode, int64_t wl, 
bit_extz_t *exz, int64_t q_tot, double e_rate, int64_t maxl, int64_t maxe, int64_t force_l, 
int64_t estimate_err, overlap_region *aux_o)
{
    clear_align(*exz); exz->thre = 0; 
    if(((ts == -1) && (te == -1))) mode = 3;///set to semi-global
    int64_t thre, ql = qe - qs, thre0, pts = -1, pte = -1, pthre = -1;
    if(ql == 0 && (te-ts) == 0) return 1;
    if((ql <= 0) || (te-ts) <= 0) return 0;
    if(estimate_err < 0) {
        if(ql > wl) estimate_err = cal_estimate_err(z, wl, qs, qe, e_rate);
        else estimate_err = ql*e_rate;
    }
    

    // fprintf(stderr, "[M::%s::ql::%ld] qs::[%ld, %ld), ts::[%ld, %ld), mode::%ld, est_err::%ld, e_rate::%f, maxe::%ld", 
    //                                                     __func__, ql, qs, qe, ts, te, mode, estimate_err, e_rate, maxe);
    if(ql <= 16) {
        if(cal_exact_exz(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, &pts, &pte, q_tot, mode)) {
            // ref_cigar_check(qstr, tu, uref, hpc_g, rref, z->y_id, z->y_pos_strand, exz);
            // fprintf(stderr, ", err::%d, thre::%d, scale::0(+)\n", exz->err, exz->thre);
            push_alnw(aux_o, exz);
            return 1;
        }
    }

    if(ql <= maxl && (estimate_err>>1) <= maxe) {
        thre = scale_ed_thre(estimate_err, maxe); 
        if(cal_exz_infi_adv(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, &pts, &pte, thre, &pthre, q_tot, mode)) {
            // ref_cigar_check(qstr, tu, uref, hpc_g, rref, z->y_id, z->y_pos_strand, exz);
            // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(+)\n", exz->err, exz->thre, thre);
            push_alnw(aux_o, exz);
            return 1;
        }

        thre0 = thre; thre = ql*e_rate; thre = scale_ed_thre(thre, maxe); 
        if(thre > thre0) {
            if(cal_exz_infi_adv(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, &pts, &pte, thre, &pthre, q_tot, mode)) {
                // ref_cigar_check(qstr, tu, uref, hpc_g, rref, z->y_id, z->y_pos_strand, exz);
                // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(-)\n", exz->err, exz->thre, thre);
                push_alnw(aux_o, exz);
                return 1;
            }
        }

        thre0 = thre; thre <<= 1; thre = scale_ed_thre(thre, maxe); 
        if(thre > thre0) {
            if(cal_exz_infi_adv(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, &pts, &pte, thre, &pthre, q_tot, mode)) {
                // ref_cigar_check(qstr, tu, uref, hpc_g, rref, z->y_id, z->y_pos_strand, exz);
                // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(-)\n", exz->err, exz->thre, thre);
                push_alnw(aux_o, exz);
                return 1;
            }
        }

        thre0 = thre; thre = ql*0.51; thre = scale_ed_thre(thre, maxe); 
        if(thre > thre0) {
            if(cal_exz_infi_adv(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, &pts, &pte, thre, &pthre, q_tot, mode)) {
                // ref_cigar_check(qstr, tu, uref, hpc_g, rref, z->y_id, z->y_pos_strand, exz);
                // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(*)\n", exz->err, exz->thre, thre);
                push_alnw(aux_o, exz);
                return 1;
            }
        }

        if(ql <= force_l) {
            thre = maxe; 
            if(cal_exz_infi_adv(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, &pts, &pte, thre, &pthre, q_tot, mode)) {
                // ref_cigar_check(qstr, tu, uref, hpc_g, rref, z->y_id, z->y_pos_strand, exz);
                // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(*)\n", exz->err, exz->thre, thre);
                push_alnw(aux_o, exz);
                return 1;
            }
        }
    }
    // fprintf(stderr, ", err::%d, thre::%d\n", INT32_MAX, exz->thre);
    // if(mode == 0) {
    //     fprintf(stderr, "[M::%s::] pstr::%.*s\n", __func__, (int32_t)tu->length, tu->seq);
    //     fprintf(stderr, "[M::%s::] tstr::%.*s\n", __func__, (int32_t)(qe-qs), qstr+qs);
    // }
    return 0;
    
}

void prt_k_mer_hit(k_mer_hit *ch_a, int64_t ch_n)
{
    int64_t k;
    for (k = 0; k < ch_n; k++) {
        fprintf(stderr, "[M::%s::k->%ld] q_pos::%u, t_pos::%u, cnt::%u, cov::%u\n", 
        __func__, k, ch_a[k].self_offset, ch_a[k].offset, ch_a[k].cnt, ch_a[k].readID);
    }
    
}

void debug_iter_k_mer_hit(k_mer_hit *ch_a, int64_t ch_n, uint64_t s, uint64_t e, int64_t ibeg, int64_t iend)
{
    int64_t i, beg = -1, end = -1;
    for (i = 0; i < ch_n; i++) {
        if((ch_a[i].self_offset >= s) && (ch_a[i].self_offset < e)) {
            if(beg == -1) beg = i;
            end = i+1;
        }
    }
    assert(ibeg==beg && iend==end);
}

int64_t chain_aln(overlap_region *z, Chain_Data *dp, k_mer_hit *ch_a, int64_t ch_n, const ul_idx_t *uref, 
hpc_t *hpc_g, All_reads *rref, char* qstr, UC_Read *tu, int64_t qs, int64_t qe, int64_t ts, 
int64_t te, int64_t mode, int64_t wl, bit_extz_t *exz, int64_t q_tot, double e_rate, 
int64_t min_chain_aln, uint64_t rid)
{
    int64_t *m, mn = 0, k; int64_t q[2], t[2], is_chain_aln = 1; 
    if(mode == 0) {
        ///wrong
        // assert((ch_n >= 2) && (ch_a[0].self_offset == qs) && (ch_a[0].offset == ts) && 
        //                     (ch_a[ch_n-1].self_offset == qe) && (ch_a[ch_n-1].offset == te));
        //does not work with CNS alignment; it also could not work here
        if(ch_n == 2 && ((qe-qs)+128) < min_chain_aln) is_chain_aln = 0;
    } else if(mode == 1) {
        // assert((ch_n >= 1) && (ch_a[0].self_offset == qs) && (ch_a[0].offset == ts));
        //does not work with CNS alignment; it also could not work here
        if(ch_n == 1 && ((qe-qs)+128) < min_chain_aln) is_chain_aln = 0;
    } else if(mode == 2) {
        // assert((ch_n >= 1) && (ch_a[ch_n-1].self_offset == qe) && (ch_a[ch_n-1].offset == te));
        //does not work with CNS alignment; it also could not work here
        if(ch_n == 1 && ((qe-qs)+128) < min_chain_aln) is_chain_aln = 0;
    }

    if(is_chain_aln) {
        mn = lchain_refine(ch_a, ch_n, NULL, dp, 50, 5000, 512, 16); m = dp->tmp;
        q[0] = qs; t[0] = ts;
        for (k = 0; k < mn; k++) {
            q[1] = ch_a[m[k]].self_offset; t[1] = ch_a[m[k]].offset;
            if(q[1] > q[0]) {
                if((t[0] != -1) && (t[1] != -1)) {
                    mode = 0;//global
                } else if((t[0] != -1) && (t[1] == -1)) {
                    mode = 1;///forward extension
                } else if((t[0] == -1) && (t[1] != -1)) {
                    mode = 2;///backward extension
                } else {
                    mode = 3;//semi-global
                }
                if(!hc_aln_exz_adv(z, uref, hpc_g, rref, qstr, tu, q[0], q[1], t[0], t[1], mode, wl, exz, q_tot, e_rate, MAX_SIN_L, MAX_SIN_E, FORCE_SIN_L, -1, NULL)) {

                }
            }
            q[0] = q[1]; t[0] = t[1];
        }
        q[1] = qe; t[1] = te;
        if(q[1] > q[0]) {
            if((t[0] != -1) && (t[1] != -1)) {
                mode = 0;//global
            } else if((t[0] != -1) && (t[1] == -1)) {
                mode = 1;///forward extension
            } else if((t[0] == -1) && (t[1] != -1)) {
                mode = 2;///backward extension
            } else {
                mode = 3;//semi-global
            }
            if(!hc_aln_exz_adv(z, uref, hpc_g, rref, qstr, tu, q[0], q[1], t[0], t[1], mode, wl, exz, q_tot, e_rate, MAX_SIN_L, MAX_SIN_E, FORCE_SIN_L, -1, NULL)) {

            }
        }
        
    } else {

    } 

    // {
    //     q[0] = qs; t[0] = ts; q[1] = qe; t[1] = te;
    //     if(q[1] > q[0]) {
    //         hc_aln_exz_adv(z, uref, hpc_g, rref, qstr, tu, q[0], q[1], t[0], t[1], mode, wl, exz, q_tot, e_rate, MAX_SIN_L, MAX_SIN_E);
    //     }
    // }
    return 0;
}

int64_t sub_base_aln(overlap_region *z, Chain_Data *dp, k_mer_hit *ch_a, int64_t ch_n, uint64_t pre_e, 
uint64_t s, uint64_t e, int64_t wl, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, 
char* qstr, UC_Read *tu, bit_extz_t *exz, double e_rate, int64_t ql, int64_t tl, 
int64_t ch_i, uint64_t rid)
{
    int64_t i = ch_i, l, ibeg, iend, mode, q[2], t[2], ch_s, ch_e;
    for (; i >= 0 && ch_a[i].self_offset >= s; i--);
    if(i<0) i = 0; ibeg=iend=-1;
    for (; i < ch_n && ch_a[i].self_offset < e; i++) {
        if((ch_a[i].self_offset >= s) && (ch_a[i].self_offset < e)) {
            if(ibeg < 0) ibeg = i;
            iend = i+1;
        }
    }
    // debug_iter_k_mer_hit(ch_a, ch_n, s, e, ibeg, iend);
    // fprintf(stderr, "***[M::%s::rid->%lu] utg%.6dl(%c), s::%lu, e::%lu, z::[%u, %u), ibeg::%ld, iend::%ld\n", 
    // __func__, rid, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], s, e, z->x_pos_s, z->x_pos_e+1, ibeg, iend);
    ch_i = i;
    if(ibeg >= 0 && iend >= 0 && iend > ibeg) {///find some anchors[ibeg, iend)
        for (i = l = ibeg; i <= iend; i++) {
            // fprintf(stderr, "\n[M::%s::i->%ld] q::%u, t::%u, cnt::%u, readID::%u\n", __func__, 
            //                                 i, ch_a[i].self_offset, ch_a[i].offset, ch_a[i].cnt, ch_a[i].readID);
            if(i == iend || is_pri_aln(ch_a[i])) {
                q[0] = q[1] = t[0] = t[1] = -1; mode = ch_s = ch_e = -1; 
                if(l < i && l < iend && is_pri_aln(ch_a[l])) {
                    q[0] = ch_a[l].self_offset; t[0] = ch_a[l].offset; ch_s = l;
                } else {
                    q[0] = s;
                }

                if(i < iend && is_pri_aln(ch_a[i])) {
                    q[1] = ch_a[i].self_offset; t[1] = ch_a[i].offset; ch_e = i;
                } else {
                    q[1] = e;
                }

                if((t[0] != -1) && (t[1] != -1)) {
                    mode = 0;//global
                } else if((t[0] != -1) && (t[1] == -1)) {
                    mode = 1;///forward extension
                } else if((t[0] == -1) && (t[1] != -1)) {
                    mode = 2;///backward extension
                } else {
                    mode = 3;//semi-global
                }

                if((mode == 0) && is_alnw(ch_a[l]) && is_alnw(ch_a[i]) 
                                            && (ch_a[l].strand == 0) && (ch_a[i].strand == 1)) {
                        ;   
                } else {
                    // fprintf(stderr, "+[M::%s::] utg%.6dl(%c), q::[%ld, %ld), t::[%ld, %ld), mode::%ld\n", 
                    //     __func__, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], q[0], q[1], t[0], t[1], mode);

                    if(mode == 1 || mode == 2) {
                        tuning_ext_offset(z, ch_a, ch_n, ql, tl, wl, MAX_SIN_L, 4, &ch_s, &ch_e, &q[0], &q[1], &t[0], &t[1], &mode);
                    } else if(ch_e >= 0) {//global
                        ch_e++;
                    }
                    if(!hc_aln_exz_adv(z, uref, hpc_g, rref, qstr, tu, q[0], q[1], t[0], t[1], mode, wl, exz, ql, e_rate, MAX_CNS_L, MAX_CNS_E, FORCE_CNS_L, -1, NULL)) {
                        if(ch_s < 0) ch_s = ibeg; if(ch_e < 0) ch_e = iend;
                        assert(ch_e > ch_s);
                        // chain_aln(z, dp, ch_a+ch_s, ch_e-ch_s, uref, hpc_g, rref, qstr, tu, q[0], q[1], t[0], t[1], mode, wl, exz, ql, e_rate, MAX_CNS_L, rid);
                    }
                }
                l = i;
            }
        }
    } else {//totoally no anchor; probably semi-global
        // fprintf(stderr, "\n***[M::%s::rid->%lu] utg%.6dl(%c), s::%lu, e::%lu, z::[%u, %u), ibeg::%ld, iend::%ld\n", 
        // __func__, rid, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], s, e, z->x_pos_s, z->x_pos_e+1, ibeg, iend);
    }

    return ch_i;
}

void cigar_gen_by_chain(overlap_region *z, Chain_Data *dp, k_mer_hit *ch_a, int64_t ch_n, ul_ov_t *ov, int64_t on, uint64_t wl,  
const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, UC_Read *tu, bit_extz_t *exz, double e_rate, 
int64_t ql, uint64_t rid)
{
    if(on <= 0) return;
    int64_t i, ch_i, tl, id = z->y_id; uint64_t pe = (uint64_t)-1;
    if(hpc_g) tl = hpc_len(*hpc_g, id);
    else if(uref) tl = uref->ug->u.a[id].len;
    else tl = Get_READ_LENGTH((*rref), id);
    // for (i = 0; i < wn; i++) z->w_list.a[i].clen = 0;///clean cigar
    // if(on > 1) {
    //     fprintf(stderr, "[M::%s::] rid::%lu, on::%ld\n", __func__, rid, on);
    // }
    // if(z->y_id == 126) prt_k_mer_hit(ch_a, ch_n);
    for (i = ch_i = 0; i < on; i++) {
        assert((i<=0)||(ov[i].qs > ov[i-1].qe));
        ch_i = sub_base_aln(z, dp, ch_a, ch_n, pe, ov[i].qs, ov[i].qe, wl, uref, hpc_g, rref, qstr, tu, exz, e_rate, ql, tl, ch_i, rid);
        pe = ov[i].qe;
    }
}


int64_t adjust_base_coordinates(overlap_region *z, k_mer_hit *ch_a, int64_t ch_n, 
ul_ov_t *res, int64_t wl, int64_t ql, int64_t tl, int64_t ch_i)
{
    res->ts = res->te = (uint32_t)-1; 
    res->sec = 3; res->qn = res->tn = (uint32_t)-1; //semi-global
    if(ch_n == 0) return ch_i; 
    int64_t i = ch_i, ibeg, iend; uint64_t s = res->qs, e = res->qe; 
    for (; i >= 0 && ch_a[i].self_offset >= s; i--); if(i < 0) i = 0;
    for (; i < ch_n && ch_a[i].self_offset < s; i++); 
    if((i >= 0) && ((ch_a[i].self_offset > s) || (i >= ch_n))) i--; ibeg = i; ///if ch_a[i].self_offset == s, do nothing
    for (i = (i>=0?i:0); i < ch_n && ch_a[i].self_offset < e; i++); iend = i;
    ch_i = i;///ch_i must be here
    ///ibeg might be < 0, iend might be == ch_n
    ///1. [s, e) contain anchors
    ///2. anchors contain [s, e)
    // fprintf(stderr, "[M::%s::] ibeg::%ld, iend::%ld, ch_n::%ld\n", __func__, ibeg, iend, ch_n);
    // if(ibeg >= 0) {
    //     fprintf(stderr, "[M::%s::] ch_a[ibeg].self_offset::%u, ch_a[ibeg].offset::%u\n", 
    //                                     __func__, ch_a[ibeg].self_offset, ch_a[ibeg].offset);
    // }
    // if(ibeg+1 >= 0) {
    //     fprintf(stderr, "[M::%s::] ch_a[ibeg+1].self_offset::%u, ch_a[ibeg+1].offset::%u\n", 
    //                                     __func__, ch_a[ibeg+1].self_offset, ch_a[ibeg+1].offset);
    // }

    // if(iend < ch_n) {
    //     fprintf(stderr, "[M::%s::] ch_a[iend].self_offset::%u, ch_a[iend].offset::%u\n", 
    //                                     __func__, ch_a[iend].self_offset, ch_a[iend].offset);
    // } 
    // if(iend > 0) {
    //     fprintf(stderr, "[M::%s::] ch_a[iend-1].self_offset::%u, ch_a[iend-1].offset::%u\n", 
    //                                     __func__, ch_a[iend-1].self_offset, ch_a[iend-1].offset);
    // }

    // if(z->y_id == 109111) {
    //     int64_t di;
    //     fprintf(stderr, "[M::%s] ch_n::%ld, s::%lu, e::%lu, ibeg::%ld, iend::%ld\n", __func__, ch_n, s, e, ibeg, iend);
    //     for (di = 0; di < ch_n; di++) {
    //         fprintf(stderr, "[M::%s::di->%ld] qoff::%u, toff::%u\n", __func__, di, 
    //                                 ch_a[di].self_offset, ch_a[di].offset);
    //     }        
    // }

    if(ibeg >= 0) {
        res->qs = ch_a[ibeg].self_offset;
        res->ts = ch_a[ibeg].offset;
        res->qn = ibeg;
    } else {//extension to left
        res->qs = 0; res->qn = (uint32_t)-1;
    }

    if(iend < ch_n) {
        res->qe = ch_a[iend].self_offset;
        res->te = ch_a[iend].offset;
        res->tn = iend;
    } else {//extension to right
        res->qe = ql; res->tn = ch_n;
    }

    return ch_i; 
}

inline int64_t translate_double_mode(uint64_t double_mode, uint64_t is_backward)
{
    if(double_mode == 0) return 0;
    if(double_mode == 4) return 3;
    if(double_mode == 1 || double_mode == 2) return double_mode;
}

int64_t fusion_chain_ovlp(overlap_region *z, k_mer_hit *ch_a, int64_t ch_n, ul_ov_t *ov, int64_t on, uint64_t wl, int64_t ql, int64_t tl)
{
    int64_t i, srt, ch_i, m, os, oe, ovlp; ul_ov_t *p;
    for (i = ch_i = 0, srt = 1; i < on; i++) {
        // fprintf(stderr, "[M::%s::i->%ld] ovq::[%u, %u)\n", __func__, i, ov->qs, ov->qe);
        ov[i].sec = 6;///do not know the aln type
        ch_i = adjust_base_coordinates(z, ch_a, ch_n, &(ov[i]), wl, ql, tl, ch_i);
        if(i > 0 && ov[i].qs < ov[i-1].qe) srt = 0;
    }
    if(on <= 1) return on;

    if(!srt) radix_sort_uov_srt_qs(ov, ov+on);
    for (i = m = 1; i < on; i++) {
        p = &(ov[m-1]);
        os = MAX(p->qs, ov[i].qs); 
        oe = MIN(p->qe, ov[i].qe);
        ovlp = oe - os;
        if(ovlp >= 0) {//merge
            p->qe = MAX(p->qe, ov[i].qe);
            p->te = MAX(p->te, ov[i].te);
            p->tn = MAX(p->tn, ov[i].tn);
        } else {//new
            ov[m++] = ov[i]; 
        }
    }

    on = m;
    return on;
}

// void chain_win_aln(overlap_region *z, Chain_Data *dp, Candidates_list *cl, int64_t qs, int64_t qe, 
// int64_t ts, int64_t te, int64_t ql, int64_t tl, int64_t wl, int64_t mode, bit_extz_t *exz)
// {

// }

int64_t ovlp_base_aln_all(overlap_region *z, k_mer_hit *ch_a, int64_t ch_n, 
int64_t soff, int64_t eoff, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, 
char* qstr, UC_Read *tu, ul_ov_t *ov, int64_t ql, int64_t tl, int64_t wl, bit_extz_t *exz, 
overlap_region *aux_o, double e_rate)
{
    int64_t ibeg, iend, i, l, mode, q[2], t[2], is_done;
    ibeg = soff; iend = eoff;
    for (l = ibeg, i = ibeg + 1; i <= iend; i++) {
        l = i - 1;
        q[0] = q[1] = t[0] = t[1] = mode = -1; is_done = 0;
        if(l >= 0) {
            q[0] = ch_a[l].self_offset; t[0] = ch_a[l].offset;
        } else {
            q[0] = ov->qs;
        }

        if(i < ch_n) {
            q[1] = ch_a[i].self_offset; t[1] = ch_a[i].offset;
        } else {
            q[1] = ov->qe;
        }

        if((t[0] != -1) && (t[1] != -1)) {
            mode = 0;//global
        } else if((t[0] != -1) && (t[1] == -1)) {
            mode = 1;///forward extension
        } else if((t[0] == -1) && (t[1] != -1)) {
            mode = 2;///backward extension
        } else {
            mode = 3;///no primary hit within [ibeg, iend] 
        }
        assert(mode != 3);
        if(mode == 1 || mode == 2) adjust_ext_offset(&(q[0]), &(q[1]), &(t[0]), &(t[1]), ql, tl, 0, mode);
        ///at cns chain, the base alignment fails; there is no anchor between soff and eoff
        if((eoff-soff<=1) && (((q[1]-q[0])>>1) < MAX_CNS_E)) {
            is_done = 0;
        } else {
            is_done = hc_aln_exz_adv(z, uref, hpc_g, rref, qstr, tu, q[0], q[1], t[0], t[1], mode, wl, exz, ql, e_rate, MAX_SIN_L, MAX_SIN_E, FORCE_SIN_L, -1, aux_o);
        }

        if(!is_done) {
            push_unmap_alnw(aux_o, q[0], q[1]-1, t[0], t[1]-1, mode);
            // chain_win_aln(z, dp, cl, q[0], q[1], t[0], t[1], ql, tl, wl, exz);
        }
        // if(aux_o->y_id == 109111) {
        //     fprintf(stderr, "<[M::%s::] utg%.6dl(%c), q::[%ld, %ld), t::[%ld, %ld), mode::%ld, is_done::%ld\n", 
        //             __func__, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], q[0], q[1], t[0], t[1], mode, is_done);
        // }
    }
    return 0;
}

void ovlp_base_aln(overlap_region *z, k_mer_hit *ch_a, int64_t ch_n, 
ul_ov_t *ov, int64_t wl, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, UC_Read *tu, 
bit_extz_t *exz, overlap_region *aux_o, double e_rate, int64_t ql, int64_t tl, uint64_t rid)
{
    int64_t ibeg, iend, i, l, mode, q[2], t[2], is_done;
    if(ov->qn == ((uint32_t)-1)) ibeg = -1;
    else ibeg = ov->qn;
    iend = ov->tn;
    assert(iend>=ibeg+1);
    // fprintf(stderr, "\n***[M::%s::rid->%lu] utg%.6dl(%c), s::%u, e::%u, z::[%u, %u), ibeg::%ld, iend::%ld, ch_n::%ld\n", 
    // __func__, rid, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], ov->qs, ov->qe, z->x_pos_s, z->x_pos_e+1, ibeg, iend, ch_n);
    for (l = ibeg, i = ibeg + 1; i <= iend; i++) {
        if(i == iend || is_pri_aln(ch_a[i])) {
            q[0] = q[1] = t[0] = t[1] = mode = -1; is_done = 0;
            if(l >= 0) {
                q[0] = ch_a[l].self_offset; t[0] = ch_a[l].offset;
            } else {
                q[0] = ov->qs;
            }

            if(i < ch_n) {
                q[1] = ch_a[i].self_offset; t[1] = ch_a[i].offset;
            } else {
                q[1] = ov->qe;
            }

            if((t[0] != -1) && (t[1] != -1)) {
                mode = 0;//global
            } else if((t[0] != -1) && (t[1] == -1)) {
                mode = 1;///forward extension
            } else if((t[0] == -1) && (t[1] != -1)) {
                mode = 2;///backward extension
            } else {
                mode = 3;///no primary hit within [ibeg, iend] 
            }
            
            if((mode == 0) && is_alnw(ch_a[l]) && is_alnw(ch_a[i]) 
                                            && (ch_a[l].strand == 0) && (ch_a[i].strand == 1)) {
                is_done = 1;
                // if(aux_o->x_id == 29033 && aux_o->y_id == 21307){
                //     fprintf(stderr, "*[M::%s::] utg%.6dl(%c), q::[%ld, %ld), t::[%ld, %ld), mode::%ld\n", 
                //             __func__, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], q[0], q[1], t[0], t[1], mode);
                // }
                // push_replace_alnw(aux_o, q[0], q[1]-1, t[0], t[1]-1, mode);///no need this
                push_replace_alnw_adv(z, wl, aux_o, q[0], q[1]-1, t[0], t[1]-1, mode);
            } else if(mode != 3) {
                if(mode == 1 || mode == 2) adjust_ext_offset(&(q[0]), &(q[1]), &(t[0]), &(t[1]), ql, tl, 0, mode);
                // if(aux_o->x_id == 29033 && aux_o->y_id == 21307){
                //     fprintf(stderr, "#[M::%s::] utg%.6dl(%c), q::[%ld, %ld), t::[%ld, %ld), mode::%ld\n", 
                //             __func__, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], q[0], q[1], t[0], t[1], mode);
                // }
                is_done = hc_aln_exz_adv(z, uref, hpc_g, rref, qstr, tu, q[0], q[1], t[0], t[1], mode, wl, exz, ql, e_rate, MAX_CNS_L, MAX_CNS_E, FORCE_CNS_L, -1, aux_o);
            }
            
            if(!is_done) {///postprocess
                // if(aux_o->x_id == 29033 && aux_o->y_id == 21307){
                //     fprintf(stderr, ">[M::%s::] utg%.6dl(%c), q::[%ld, %ld), t::[%ld, %ld), mode::%ld\n", 
                //             __func__, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], q[0], q[1], t[0], t[1], mode);
                // }
                is_done = ovlp_base_aln_all(z, ch_a, ch_n, l, i, uref, hpc_g, rref, qstr, tu, ov, ql, tl, wl, exz, aux_o, e_rate);
            }
            // if(aux_o->x_id == 29033 && aux_o->y_id == 21307){
            //     fprintf(stderr, "-is_done::%ld[M::%s::] utg%.6dl(%c), q::[%ld, %ld), t::[%ld, %ld), mode::%ld, l::%ld, i::%ld, ch_n::%ld\n", 
            //             is_done, __func__, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], q[0], q[1], t[0], t[1], mode, l, i, ch_n);
            // }
            // if(rid == (uint64_t)-1) {
            //     fprintf(stderr, "+is_done::%ld[M::%s::] utg%.6dl(%c), q::[%ld, %ld), t::[%ld, %ld), mode::%ld, l::%ld, i::%ld, ch_n::%ld\n", 
            //             is_done, __func__, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], q[0], q[1], t[0], t[1], mode, l, i, ch_n);
            // }
            l = i;
        }
    }

}

inline void push_khit(Candidates_list *res, int32_t xs, int32_t ys, uint32_t len, uint32_t h_khit, uint32_t *ic)
{
    uint32_t p, c; k_mer_hit *z;
    c = ((len >= h_khit)?1:2); if(ic) c = *ic; c <<= 8;
    if(len > 0) {
        while (len >= (0xffu)) {
            p = (c + (0xffu)); 
            kv_pushp_cl(k_mer_hit, (*res), &z); 
            memset(z, 0, sizeof((*z)));
            z->self_offset = xs; z->offset = ys; z->cnt = p;
            len -= (0xffu);
        }
        if(len) {
            p = (c + len); 
            kv_pushp_cl(k_mer_hit, (*res), &z); 
            memset(z, 0, sizeof((*z)));
            z->self_offset = xs; z->offset = ys; z->cnt = p;
        }
    } else {
        p = (c + len); 
        kv_pushp_cl(k_mer_hit, (*res), &z); 
        memset(z, 0, sizeof((*z)));
        z->self_offset = xs; z->offset = ys; z->cnt = p;
    }
}

uint32_t extract_exact_cigar(asg16_v *ez, int32_t ps, int32_t ts, int32_t pmin, int32_t pmax, 
int32_t tmin, int32_t tmax, Candidates_list *res, int32_t minl, int64_t min_w_l, int64_t h_khit)
{
    uint32_t ci = 0, cl, occ = 0; uint16_t c; 
    int32_t pi = ps, ti = ts, p[2], t[2], poff, toff, maxl;
    int32_t pos, poe, tos, toe, l; poff = toff = maxl = -1;
    while (ci < ez->n && pi < pmax && ti < tmax) {
        ci = pop_trace(ez, ci, &c, &cl);
        if(c == 0) {
            p[0] = pi; p[1] = pi + cl;
            t[0] = ti; t[1] = ti + cl;
            pos = MAX(p[0], pmin); poe = MIN(p[1], pmax);
            tos = MAX(t[0], tmin); toe = MIN(t[1], tmax);
            if((poe > pos) && (toe > tos)) {
                l = poe - pos;
                if(l == (toe - tos)) {
                    poe--; toe--;
                    if(l > maxl) {
                        poff = poe; toff = toe; maxl = l;
                    }
                    if(l >= minl) {
                        push_khit(res, toe, poe, l, h_khit, NULL); occ++;
                    }
                }
            }
            pi+=cl; ti+=cl;
        } else if(c == 1) {
            pi+=cl; ti+=cl;
        } else if(c == 2) {///more p
            pi+=cl;
        } else if(c == 3) {
            ti+=cl;
        }
    }
    ///(ts >= tmin) && (ti >= (min_w_l + ts)): here is a whole window
    if(maxl > 0 && maxl < minl && (ts >= tmin) && (ti >= (min_w_l + ts))) {
        uint32_t w = 3;
        push_khit(res, toff, poff, maxl, h_khit, &w); occ++;
    }
    return occ;
}

int64_t debug_k_mer_hit_retrive(k_mer_hit *z, hpc_t *hpc_g, All_reads *rref, const ul_idx_t *uref, 
char* qstr, UC_Read *tu, int64_t id, int64_t rev)
{
    int64_t qs, qe, ts, te; char *q_string, *t_string;
    qe = z->self_offset+1; qs = qe - (z->cnt&(0xffu));
    te = z->offset+1; ts = te - (z->cnt&(0xffu));
    if(qe == qs) return 1;
    q_string = qstr + qs; 
    t_string = retrieve_str_seq_exz(tu, ts, te-ts, -1, -1, rev, uref, hpc_g, rref, id);
    fprintf(stderr, "[M::%s::] q_string::%.*s\n", __func__, (int32_t)(qe-qs), q_string);
    fprintf(stderr, "[M::%s::] t_string::%.*s\n", __func__, (int32_t)(te-ts), t_string);
    if(memcmp(q_string, t_string, qe-qs)) {
        fprintf(stderr, "[M::%s::] qsite::%u, tsite::%u\n", __func__, z->self_offset, z->offset);
        return 0;
    }
    return 1;
}

int64_t gen_single_khit(Candidates_list *cl, int64_t ch_n, int64_t h_khit, int64_t mode, int64_t qs, int64_t qe, int64_t ts, int64_t te, int64_t max_skip, int64_t max_iter, int64_t rid)
{
    // if(ch_n != 3 || mode != 0 || qs != 171728 || qe != 172258) return 0;
    k_mer_hit *ch_a = cl->list + cl->length; int64_t k, i, j, occ, m, ncn, prefix, suffix, srt = 1;
    prefix = suffix = 0;
    if(mode == 0 || mode == 2) suffix = 1;
    if(mode == 0 || mode == 1) prefix = 1;
    // if(ch_n == 2 && mode == 2 && qe - qs == 2419 && te - ts == 2419) {
        // fprintf(stderr, "[M::%s::mode->%ld] ch_n::%ld, q::[%ld, %ld), t::[%ld, %ld)\n", 
        //     __func__, mode, ch_n, qs, qe, ts, te);
    // }
    
    for (k = occ = m = 0; k < ch_n; k++) {
        // if(ch_n == 2 && mode == 2 && qe - qs == 2419 && te - ts == 2419) {
        // fprintf(stderr, "+i::%ld[M::%s::] x::[%u, %u), y::[%u, %u), cnt::%u\n", k, __func__, 
        //     ch_a[k].self_offset+1-(ch_a[k].cnt&((uint32_t)(0xffu))), ch_a[k].self_offset+1,
        //     ch_a[k].offset+1-(ch_a[k].cnt&((uint32_t)(0xffu))), ch_a[k].offset+1, (ch_a[k].cnt&(0xffu)));
        // }
        if(!(ch_a[k].cnt&(0xffu))) continue;
        occ++;
        if((ch_a[k].cnt&(0xffu)) > 1) occ++;
        ch_a[m++] = ch_a[k];
    }
    ch_n = m; if(!ch_n) return ch_n;
    occ += prefix + suffix;
    // if(ch_n == 2 && mode == 2 && qe - qs == 2419 && te - ts == 2419) {
        // fprintf(stderr, "+[M::%s::] occ::%ld\n", __func__, occ);
    // }

    ncn = occ + cl->length;
    if(cl->size < ncn) {
        cl->size = ncn;
        REALLOC(cl->list, cl->size);
        // cl->list = (k_mer_hit*)realloc(cl->list, (sizeof((*(cl->list)))*cl->length));
    }
    ch_a = cl->list + cl->length; assert((cl->length+occ)<= cl->size);

    k_mer_hit cht;
    ///global or backward
    if(suffix) {
        cht.self_offset = qe;
        cht.offset = te;
        cht.cnt = 1; cht.readID = 1;//make it as primary chain
        cht.strand = 0;
        ch_a[--occ] = cht;
        // fprintf(stderr, "occ::%ld[M::%s::] x::%u, y::%u, cnt::%u, cov::%u\n", occ, __func__, 
        //     cht.self_offset, cht.offset, cht.cnt, cht.readID);
    }
    for (k = ch_n-1; k >= 0; k--) {
        if(!(ch_a[k].cnt&(0xffu))) continue;
        ///end
        cht.self_offset = ch_a[k].self_offset+1;
        cht.offset = ch_a[k].offset+1;
        cht.strand = 0;
        cht.cnt = cht.readID = (ch_a[k].cnt&(0xffu));
        //make it as non-primary chain
        if((ch_a[k].cnt&(0xffu)) < h_khit) cht.readID = cht.cnt + 1;
        ch_a[--occ] = cht;
        // fprintf(stderr, "occ::%ld[M::%s::] x::%u, y::%u, cnt::%u, cov::%u\n", occ, __func__, 
        //     cht.self_offset, cht.offset, cht.cnt, cht.readID);

        if((ch_a[k].cnt&(0xffu)) > 1) {
            ///start
            cht.self_offset = ch_a[k].self_offset+1-(ch_a[k].cnt&(0xffu));
            cht.offset = ch_a[k].offset+1-(ch_a[k].cnt&(0xffu));
            cht.strand = 0;
            cht.cnt = cht.readID = (ch_a[k].cnt&(0xffu));
            //make it as non-primary chain
            if((ch_a[k].cnt&(0xffu)) < h_khit) cht.readID = cht.cnt + 1;
            ch_a[--occ] = cht;
            // fprintf(stderr, "occ::%ld[M::%s::] x::%u, y::%u, cnt::%u, cov::%u\n", occ, __func__, 
            // cht.self_offset, cht.offset, cht.cnt, cht.readID);
        }
    }
    
    if(prefix) { ///global or forward
        cht.self_offset = qs;
        cht.offset = ts;
        cht.cnt = 1; cht.readID = 1;//make it as primary chain
        cht.strand = 0;
        ch_a[--occ] = cht;
        // fprintf(stderr, "occ::%ld[M::%s::] x::%u, y::%u, cnt::%u, cov::%u\n", occ, __func__, 
        //     cht.self_offset, cht.offset, cht.cnt, cht.readID);
    }
    // if(ch_n == 2 && mode == 2 && qe - qs == 2419 && te - ts == 2419) {
        // fprintf(stderr, "-[M::%s::] occ::%ld\n", __func__, occ);
    // }
    // if(!(occ == 0)) {
    //     fprintf(stderr, "[M::%s] rid::%ld, name::%.*s\n", __func__, rid,
    //      (int32_t)UL_INF.nid.a[rid].n, UL_INF.nid.a[rid].a);
    // }
    assert(occ == 0);
    ch_n = occ = ncn - cl->length;
    uint64_t q[2], t[2]; 
    q[0] = q[1] = t[0] = t[1] = (uint64_t)-1;
    if(prefix) {
        q[0] = qs; t[0] = ts;
    }
    if(suffix) {
        q[1] = qe; t[1] = te;
    }

    // for (k = 0; k < ch_n; k++) {
    //     fprintf(stderr, "<i::%ld[M::%s::] x::%u, y::%u, cnt::%u, cov::%u\n", k, __func__, 
    //         ch_a[k].self_offset, ch_a[k].offset, ch_a[k].cnt, ch_a[k].readID);
    // }
    for (k = m = occ = 0; k < ch_n; k++) {
        if((k>0) && (ch_a[k].self_offset==q[0]) && (ch_a[k].offset=t[0])) continue;
        if(((k+1)<ch_n) && (ch_a[k].self_offset==q[1]) && (ch_a[k].offset=t[1])) continue;
        if(m>0) {
            if((ch_a[k].self_offset>ch_a[m-1].self_offset) && (ch_a[k].offset>ch_a[m-1].offset)) {
                occ++;
            }
            if(ch_a[k].self_offset<=ch_a[m-1].self_offset) srt = 0;
        } else {
            occ++;
        }
        ch_a[m++] = ch_a[k];
    }
    ch_n = m;
    if(occ == ch_n) return ch_n;///already colinear
    if(!srt) {
        radix_sort_k_mer_hit_self(ch_a, ch_a + ch_n);
        for (i = 1, j = 0; i <= ch_n; i++) {
            if (i == ch_n || ch_a[i].self_offset != ch_a[j].self_offset) {
                if(i - j > 1) radix_sort_k_mer_hit_off(ch_a+j, ch_a+i);
                j = i;
            }
        }
    }

    // for (k = 0; k < ch_n; k++) {
    //     fprintf(stderr, ">i::%ld[M::%s::] x::%u, y::%u, cnt::%u, cov::%u\n", k, __func__, 
    //         ch_a[k].self_offset, ch_a[k].offset, ch_a[k].cnt, ch_a[k].readID);
    // }
    occ = ch_n;
    ch_n = lchain_simple(ch_a+prefix, ch_n-prefix-suffix, ch_a+prefix, &(cl->chainDP), max_skip, max_iter);
    ch_n += prefix + suffix; if(suffix) ch_a[ch_n-1] = ch_a[occ-1];
    // for (k = 0; k < ch_n; k++) {
    //     fprintf(stderr, "-i::%ld[M::%s::] x::%u, y::%u, cnt::%u, cov::%u\n", k, __func__, 
    //         ch_a[k].self_offset, ch_a[k].offset, ch_a[k].cnt, ch_a[k].readID);
    // }
    
    return ch_n;
}

///[qs, qe) && [ts, te)
int64_t gen_win_chain(overlap_region *z, Candidates_list *cl, int64_t qs, int64_t qe, int64_t ts, int64_t te, 
int64_t wl, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, UC_Read *tu, bit_extz_t *exz, 
int64_t ql, int64_t tl, double e_rate, int64_t h_khit, int64_t mode, int64_t rid)
{
    assert(mode < 3);
    int64_t k, ws, we, os, oe, wsk, rcn = cl->length, ncn, occ = 0, ovlp, wn = z->w_list.n; asg16_v ez; uint32_t w = 1;
    ws = qs; if(ws < z->x_pos_s) ws = z->x_pos_s; 
    we = qe-1; if(we > z->x_pos_e) we = z->x_pos_e; 
    wsk = get_win_id_by_s(z, ((ws/wl)*wl), wl, NULL); 
    // if(rid == 7) {
    //     fprintf(stderr, "[M::%s::]\tutg%.6u%c\t%u\t%u\t%c\tutg%.6u%c\t%u\t%u\tq::[%ld,%ld)\tw::[%ld,%ld]\twn::%ld\twsk::%ld\n", __func__, 
	// 			z->x_id+1, "lc"[uref->ug->u.a[z->x_id].circ], 
	// 			z->x_pos_s, z->x_pos_e+1, "+-"[z->y_pos_strand],
	// 			z->y_id+1, "lc"[uref->ug->u.a[z->y_id].circ], 
	// 			z->y_pos_s, z->y_pos_e+1, qs, qe, ws, we, wn, wsk);
    // }
    // assert((ws>=z->w_list.a[wsk].x_start) && (ws<=z->w_list.a[wsk].x_end));
    // wek = get_win_id_by_e(z, ((we/wl)*wl), wl, NULL);
    // assert((we>=z->w_list.a[wek].x_start) && (we<=z->w_list.a[wek].x_end));
    for(wsk=((wsk<wn)?(wsk):(wn-1)); wsk<wn && qs>z->w_list.a[wsk].x_end; wsk++);
    for(wsk=((wsk<wn)?(wsk):(wn-1)); wsk>=0 && qs<z->w_list.a[wsk].x_start; wsk--);
    if(wsk < 0) wsk = 0; ///qs >= z->w_list.a[wsk].x_start && qs <= z->w_list.a[wsk].x_end  
    ///global or forward
    if(mode == 0 || mode == 1) push_khit(cl, qs, ts, 0, 0, &w);
    //[ws, we] && [wsk, wek]; [qs, qe) && [ts, te)
    for (k = wsk; k<wn && z->w_list.a[k].x_start<qe; k++) {
        if(z->w_list.a[k].y_end == -1) continue;
        ws = z->w_list.a[k].x_start; we = z->w_list.a[k].x_end+1;
        os = MAX(qs, ws); oe = MIN(qe, we); ovlp = ((oe>os)? (oe-os):0);
        if(!ovlp) continue;
        if(!(z->w_list.a[k].clen)) {
            gen_backtrace_adv_exz(&(z->w_list.a[k]), z, NULL, NULL, uref, qstr, tu->seq, exz, z->y_pos_strand, z->y_id);
        }
        ez.a = z->w_list.c.a + z->w_list.a[k].cidx; 
        ez.n = ez.m = z->w_list.a[k].clen; 
        occ += extract_exact_cigar(&ez, z->w_list.a[k].y_start, z->w_list.a[k].x_start, ts, te, qs, qe, cl, 10, wl, h_khit);
    }
    ///global or backward
    if(mode == 0 || mode == 2) push_khit(cl, qe-1, te-1, 0, 0, &w);
    // fprintf(stderr, "[M::%s::] rcn::%ld, cl->length::%lld\n", __func__, rcn, cl->length);
    if(!occ) {
        cl->length = rcn; return 0;
    }
    ncn = cl->length; cl->length = rcn;
    k_mer_hit *ch_a = cl->list + rcn; int64_t ch_n0 = ncn - rcn, ch_n;
    int64_t max_skip, max_iter, max_dis, quick_check; double chn_pen_gap, chn_pen_skip;
	set_lchain_dp_op(0, h_khit, &max_skip, &max_iter, &max_dis, &chn_pen_gap, &chn_pen_skip, &quick_check);
    max_dis = MAX_SIN_L>>1;
    // for (k = 0; k < ch_n0; k++) {
    //     assert(debug_k_mer_hit_retrive(&(ch_a[k]), hpc_g, rref, uref, qstr, tu, z->y_id, z->y_pos_strand));
    // }
    ch_n = lchain_qdp_fix(ch_a, ch_n0, &(cl->chainDP), max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, 
                e_rate, ql, tl, 1, ((mode==0)||(mode==1))?1:0,  ((mode==0)||(mode==2))?1:0);
    for (k = occ = 0; k < ch_n; k++) {
        ch_a[k] = ch_a[cl->chainDP.tmp[k]];
        if((ch_a[k].cnt&(0xffu))) occ++;
        // assert(debug_k_mer_hit_retrive(&(ch_a[k]), hpc_g, rref, uref, qstr, tu, z->y_id, z->y_pos_strand));
    }
    // fprintf(stderr, "[M::%s::] ch_n0::%ld, ch_n::%ld, mode::%ld, ql::%ld, tl::%ld, occ::%ld\n", 
    // __func__, ch_n0, ch_n, mode, qe-qs, te-ts, occ);
    if(occ <= 0) return 0;
    ch_n = gen_single_khit(cl, ch_n, h_khit, mode, qs, qe, ts, te, max_skip, max_iter, rid);
    return ch_n;
}


void rechain_aln(overlap_region *z, Candidates_list *cl, overlap_region *aux_o, int64_t aux_i, int64_t wl, 
const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, UC_Read *tu, bit_extz_t *exz, double e_rate, 
int64_t ql, int64_t tl, int64_t h_khit, int64_t rid)
{
    int64_t rcn = cl->length, ch_n, qs, qe, ts, te, mode, an0, an, todo; 
    k_mer_hit *ch_a; ul_ov_t idx; uint8_t q[2], t[2];
    ///[qs, qe) && [ts, te)
    qs = aux_o->w_list.a[aux_i].x_start; qe = aux_o->w_list.a[aux_i].x_end+1;
    ts = aux_o->w_list.a[aux_i].y_start; te = aux_o->w_list.a[aux_i].y_end+1;
    if(qe - qs < FORCE_SIN_L || te - ts < FORCE_SIN_L) return;
    mode = aux_o->w_list.a[aux_i].error_threshold;
    ch_n = gen_win_chain(z, cl, qs, qe, ts, te, wl, uref, hpc_g, rref, qstr, tu, exz, ql, tl, e_rate, h_khit, mode, rid);
    ch_a = cl->list + rcn;
    if(ch_n) {
        idx.ts = idx.te = (uint32_t)-1; idx.qs = 0; idx.qe = ql; todo = 1;
        if(mode == 0) {//global
            idx.qn = 0; idx.tn = ch_n - 1;
            idx.qs = ch_a[idx.qn].self_offset;
            idx.ts = ch_a[idx.qn].offset;
            idx.qe = ch_a[idx.tn].self_offset;
            idx.te = ch_a[idx.tn].offset;
            assert(ch_a[0].self_offset == qs && ch_a[0].offset == ts);
            assert(ch_a[ch_n-1].self_offset == qe && ch_a[ch_n-1].offset == te);
            if(ch_n <= 2) todo = 0;
        } else if(mode == 1) {//forward ext
            idx.qn = 0; idx.tn = ch_n;
            idx.qs = ch_a[idx.qn].self_offset;
            idx.ts = ch_a[idx.qn].offset;
            idx.qe = ql; 
            assert(ch_a[0].self_offset == qs && ch_a[0].offset == ts);
            if(ch_n <= 1) todo = 0;
        } else if(mode == 2) {///backward ext
            idx.qn = (uint32_t)-1; idx.tn = ch_n-1;
            idx.qs = 0; 
            idx.qe = ch_a[idx.tn].self_offset;
            idx.te = ch_a[idx.tn].offset;
            assert(ch_a[ch_n-1].self_offset == qe && ch_a[ch_n-1].offset == te);
            if(ch_n <= 1) todo = 0;
        }
        if(todo) {
            an0 = aux_o->w_list.n;
            // if(z->x_id == 29033 && z->y_id == 21307) {
            //     fprintf(stderr, "[M::%s]\tan0::%ld\tq::[%u,\t%u)\tt::[%u,\t%u)\tlw::%u\trw::%u\n", __func__, an0,
            //     idx.qs, idx.qe, idx.ts, idx.te, idx.qn, idx.tn);
            // }
            ovlp_base_aln(z, ch_a, ch_n, &idx, wl, uref, hpc_g, rref, qstr, tu, exz, aux_o, e_rate, ql, tl, (uint64_t)-1);
            an = aux_o->w_list.n; q[0] = q[1] = t[0] = t[1] = 0; todo = 0;
            // if(z->x_id == 29033 && z->y_id == 21307) {
            //     fprintf(stderr, "[M::%s]\tan::%ld\n", __func__, an);
            // }
            // fprintf(stderr, "[M::%s::] awn0::%ld, awn::%lu\n", __func__, an0, an);
            ///old unaligned window could be replaced by the new aligned window
            if((an == (an0 + 1)) && (!(is_ualn_win(aux_o->w_list.a[an-1])))) {
                if(aux_o->w_list.a[aux_i].x_start == aux_o->w_list.a[an-1].x_start) q[0] = 1;
                if(aux_o->w_list.a[aux_i].x_end == aux_o->w_list.a[an-1].x_end) q[1] = 1;
                if(aux_o->w_list.a[aux_i].y_start == aux_o->w_list.a[an-1].y_start) t[0] = 1;
                if(aux_o->w_list.a[aux_i].y_end == aux_o->w_list.a[an-1].y_end) t[1] = 1;
                if((mode == 0) && q[0] && q[1] && t[0] && t[1]) todo = 1;
                if((mode == 1) && q[0] && t[0]) todo = 1;
                if((mode == 2) && q[1] && t[1]) todo = 1;
                if(todo) {
                    aux_o->w_list.a[aux_i] = aux_o->w_list.a[an-1]; aux_o->w_list.n--;
                }
            }
            // if(an > an0) {///should always > 0 as there are unmapped windows
            // }
            // aux_o->w_list.n = an0;
        }
    }
    cl->length = rcn;///must reset!!!!
}

void debug_overlap_region(overlap_region *au, char* qstr, UC_Read *tu, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref)
{
    int64_t wn = au->w_list.n, k; bit_extz_t ez;
    for (k = 0; k < wn; k++) {
        assert((k<=0)||((au->w_list.a[k].x_start>au->w_list.a[k-1].x_end)
                                        &&(au->w_list.a[k].y_start>au->w_list.a[k-1].y_end)));
        assert(au->w_list.a[k].x_end>au->w_list.a[k].x_start);
        assert(au->w_list.a[k].y_end>au->w_list.a[k].y_start);
        if(is_ualn_win(au->w_list.a[k]) || is_est_aln(au->w_list.a[k])) continue;
        ez.cigar.a = au->w_list.c.a + au->w_list.a[k].cidx; 
        ez.cigar.n = ez.cigar.m = au->w_list.a[k].clen;
        ez.ts = au->w_list.a[k].x_start; ez.te = au->w_list.a[k].x_end;
        ez.ps = au->w_list.a[k].y_start; ez.pe = au->w_list.a[k].y_end;
        ez.err = au->w_list.a[k].error;
        ref_cigar_check(qstr, tu, uref, hpc_g, rref, au->y_id, au->y_pos_strand, &ez);
    }
}

void update_overlap_region(overlap_region *des, overlap_region *src, int64_t xl, int64_t yl)
{
    kv_resize(uint16_t, des->w_list.c, src->w_list.c.n);
    des->w_list.c.n = src->w_list.c.n;
    memcpy(des->w_list.c.a, src->w_list.c.a, src->w_list.c.n*(sizeof((*(src->w_list.c.a)))));

    kv_resize(window_list, des->w_list, src->w_list.n);
    des->w_list.n = src->w_list.n;
    memcpy(des->w_list.a, src->w_list.a, src->w_list.n*(sizeof((*(src->w_list.a)))));

    if(src->w_list.n) {
        des->x_pos_s = src->w_list.a[0].x_start; des->x_pos_e = src->w_list.a[src->w_list.n-1].x_end; 
        des->y_pos_s = src->w_list.a[0].y_start; des->y_pos_e = src->w_list.a[src->w_list.n-1].y_end;
    }

    int64_t xr, yr; 
    if(des->x_pos_s <= des->y_pos_s) {
        des->y_pos_s -= des->x_pos_s; des->x_pos_s = 0;
    } else {
        des->x_pos_s -= des->y_pos_s; des->y_pos_s = 0;
    }

    xr = xl-des->x_pos_e-1; yr = yl-des->y_pos_e-1;
    if(xr <= yr) {
        des->x_pos_e = xl-1; des->y_pos_e += xr;        
    } else {
        des->y_pos_e = yl-1; des->x_pos_e += yr; 
    }
}

void cigar_gen_by_chain_adv(overlap_region *z, Candidates_list *cl, int64_t ch_idx, int64_t ch_n,
ul_ov_t *ov, int64_t on, uint64_t wl, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, 
UC_Read *tu, bit_extz_t *exz, overlap_region *aux_o, double e_rate, int64_t ql, uint64_t rid, int64_t h_khit)
{
    if(on <= 0) return;
    int64_t i, tl, id = z->y_id, m; 
    k_mer_hit *ch_a = cl->list + ch_idx;
    if(hpc_g) tl = hpc_len(*hpc_g, id);
    else if(uref) tl = uref->ug->u.a[id].len;
    else tl = Get_READ_LENGTH((*rref), id);
    // fprintf(stderr, "\n[M::%s::rid->%ld] utg%.6dl(%c), z::[%u, %u)\n", 
    // __func__, rid, (int32_t)z->y_id+1, "+-"[z->y_pos_strand],  z->x_pos_s, z->x_pos_e+1);
    on = fusion_chain_ovlp(z, ch_a, ch_n, ov, on, wl, ql, tl);
    aux_o->w_list.n = aux_o->w_list.c.n = 0; 
    aux_o->y_id = z->y_id; aux_o->y_pos_strand = z->y_pos_strand;
    aux_o->x_pos_s = z->x_pos_s; aux_o->x_pos_e = z->x_pos_e;
    aux_o->y_pos_s = z->y_pos_s; aux_o->y_pos_e = z->y_pos_e;

    
    for (i = 0; i < on; i++) {
        // if(aux_o->y_id == 109111) {
        //     fprintf(stderr, "[M::%s::i->%ld] ovq::[%u, %u), ovt::[%u, %u), hits::[%d, %d)\n", __func__, i, 
        //                                 ov->qs, ov->qe, ov->ts, ov->te,
        //                                 (ov->qn!=((uint32_t)-1))?(int32_t)ov->qn:-1, (int32_t)ov->tn);
        // }
        assert((i<=0)||(ov[i].qs>ov[i-1].qe));
        ovlp_base_aln(z, ch_a, ch_n, &(ov[i]), wl, uref, hpc_g, rref, qstr, tu, exz, aux_o, e_rate, ql, tl, rid);
    }

    int64_t aux_n = aux_o->w_list.n;
    ///for debug
    // if(aux_o->y_id == 109111) {
    //     for (i = 0; i < ((int64_t)aux_o->w_list.n); i++) {
    //         fprintf(stderr, "-0-[aln::i->%ld::ql->%d] q::[%d, %d), t::[%d, %d), err::%d, clen::%u, extra_end::%d, mode::%d\n", i, 
    //                     aux_o->w_list.a[i].x_end+1-aux_o->w_list.a[i].x_start,
    //                     aux_o->w_list.a[i].x_start, aux_o->w_list.a[i].x_end+1, 
    //                     aux_o->w_list.a[i].y_start, aux_o->w_list.a[i].y_end+1, 
    //                     aux_o->w_list.a[i].error, aux_o->w_list.a[i].clen, 
    //                     aux_o->w_list.a[i].extra_end, aux_o->w_list.a[i].error_threshold);
    //     }
    // }

    for (i = 0; i < aux_n; i++) {
        if(!(is_ualn_win(aux_o->w_list.a[i]))) continue;
        //will overwrite ch_a; does not matter
        rechain_aln(z, cl, aux_o, i, wl, uref, hpc_g, rref, qstr, tu, exz, e_rate, ql, tl, h_khit, rid);
    }

    ///for debug
    // if(aux_o->y_id == 109111) {
    //     for (i = 0; i < ((int64_t)aux_o->w_list.n); i++) {
    //         fprintf(stderr, "-2-[aln::i->%ld::ql->%d] q::[%d, %d), t::[%d, %d), err::%d, clen::%u, mode::%d\n", i, 
    //                     aux_o->w_list.a[i].x_end+1-aux_o->w_list.a[i].x_start,
    //                     aux_o->w_list.a[i].x_start, aux_o->w_list.a[i].x_end+1, 
    //                     aux_o->w_list.a[i].y_start, aux_o->w_list.a[i].y_end+1, 
    //                     aux_o->w_list.a[i].error, aux_o->w_list.a[i].clen, aux_o->w_list.a[i].error_threshold);
    //     }
    // }

    if(((int64_t)aux_o->w_list.n) > aux_n) {
        for (i = m = 0; i < ((int64_t)aux_o->w_list.n); i++) {
            if((i < aux_n) && (is_ualn_win(aux_o->w_list.a[i]))) continue;
            aux_o->w_list.a[m++] = aux_o->w_list.a[i];
        }
        aux_o->w_list.n = m;
        radix_sort_window_list_xs_srt(aux_o->w_list.a, aux_o->w_list.a+aux_o->w_list.n);
    }

    ///update z by aux_o
    update_overlap_region(z, aux_o, ql, tl);

    // debug_overlap_region(aux_o, qstr, tu, uref, hpc_g, rref);


    // ch_a = cl->list + ch_idx; //update
    // for (i = 0; i < wn; i++) z->w_list.a[i].clen = 0;///clean cigar
    // if(on > 1) {
    //     fprintf(stderr, "[M::%s::] rid::%lu, on::%ld\n", __func__, rid, on);
    // }
    // if(z->y_id == 126) prt_k_mer_hit(ch_a, ch_n);
    // for (i = ch_i = 0; i < on; i++) {
    //     assert((i<=0)||(ov[i].qs > ov[i-1].qe));
    //     ov[i].sec = 16;///do not know the aln type
    //     ch_i = sub_base_aln(z, dp, ch_a, ch_n, pe, ov[i].qs, ov[i].qe, wl, uref, hpc_g, rref, qstr, tu, exz, e_rate, ql, tl, ch_i, rid);
    //     pe = ov[i].qe;
    // }
}

void ovlp_base_direct(overlap_region *z, k_mer_hit *ch_a, int64_t ch_n, 
ul_ov_t *ov, int64_t wl, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, UC_Read *tu, 
bit_extz_t *exz, overlap_region *aux_o, double e_rate, int64_t ql, int64_t tl, uint64_t rid)
{
    int64_t ibeg, iend, i, l, mode, q[2], t[2], is_done;
    if(ov->qn == ((uint32_t)-1)) ibeg = -1;
    else ibeg = ov->qn;
    iend = ov->tn;
    assert(iend>=ibeg+1);
    // if(z->x_id == 57 && z->y_id == 2175) {
    //     fprintf(stderr, "\n***[M::%s::rid->%lu] utg%.6dl(%c), s::%u, e::%u, z::[%u, %u), ibeg::%ld, iend::%ld, ch_n::%ld\n", 
    //     __func__, rid, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], ov->qs, ov->qe, z->x_pos_s, z->x_pos_e+1, ibeg, iend, ch_n);
    // }
    for (l = ibeg, i = ibeg + 1; i <= iend; i++) {
        q[0] = q[1] = t[0] = t[1] = mode = -1; is_done = 0;
        if(l >= 0) {
            q[0] = ch_a[l].self_offset; t[0] = ch_a[l].offset;
        } else {
            q[0] = ov->qs;
        }

        if(i < ch_n) {
            q[1] = ch_a[i].self_offset; t[1] = ch_a[i].offset;
        } else {
            q[1] = ov->qe;
        }

        if((t[0] != -1) && (t[1] != -1)) {
            mode = 0;//global
        } else if((t[0] != -1) && (t[1] == -1)) {
            mode = 1;///forward extension
        } else if((t[0] == -1) && (t[1] != -1)) {
            mode = 2;///backward extension
        } else {
            mode = 3;///no primary hit within [ibeg, iend] 
        }

        if(mode == 1 || mode == 2) adjust_ext_offset(&(q[0]), &(q[1]), &(t[0]), &(t[1]), ql, tl, 0, mode);
        // if(z->x_id == 57 && z->y_id == 2175) {
        //     fprintf(stderr, "#[M::%s::] utg%.6dl(%c), q::[%ld, %ld), t::[%ld, %ld), mode::%ld\n", 
        //             __func__, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], q[0], q[1], t[0], t[1], mode);
        // }
        is_done = hc_aln_exz_adv(z, uref, hpc_g, rref, qstr, tu, q[0], q[1], t[0], t[1], mode, wl, exz, ql, e_rate, 
        MAX_SIN_L, MAX_SIN_E, FORCE_SIN_L, -1, aux_o);

        // if(z->x_id == 57 && z->y_id == 2175) {
        //     fprintf(stderr, "-is_done::%ld[M::%s::] utg%.6dl(%c), q::[%ld, %ld), t::[%ld, %ld), mode::%ld\n", 
        //             is_done, __func__, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], q[0], q[1], t[0], t[1], mode);
        // }
                
        if(!is_done) {///postprocess
            push_unmap_alnw(aux_o, q[0], q[1]-1, t[0], t[1]-1, mode);
        }
        l = i;
    }
}

void cigar_gen_by_chain_adv_local(overlap_region *z, Candidates_list *cl, ul_ov_t *ov, int64_t on, uint64_t wl, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, 
UC_Read *tu, bit_extz_t *exz, overlap_region *aux_o, double e_rate, int64_t ql, uint64_t rid, int64_t h_khit)
{
    if(on <= 0) return;
    int64_t ch_idx = z->shared_seed, ch_n;
    int64_t i, tl, id = z->y_id, m; 
    k_mer_hit *ch_a = cl->list + ch_idx;
    if(hpc_g) tl = hpc_len(*hpc_g, id);
    else if(uref) tl = uref->ug->u.a[id].len;
    else tl = Get_READ_LENGTH((*rref), id);
    for (i = ch_idx; i < cl->length && cl->list[i].readID == cl->list[ch_idx].readID; i++); ch_n = i-ch_idx;
    
    // fprintf(stderr, "[M::%s::rid->%ld] utg%.6dl(%c), z::[%u, %u)\n", 
    // __func__, rid, (int32_t)z->y_id+1, "+-"[z->y_pos_strand],  z->x_pos_s, z->x_pos_e+1);
    on = fusion_chain_ovlp(z, ch_a, ch_n, ov, on, wl, ql, tl);
    aux_o->w_list.n = aux_o->w_list.c.n = 0; 
    aux_o->y_id = z->y_id; aux_o->y_pos_strand = z->y_pos_strand;
    aux_o->x_pos_s = z->x_pos_s; aux_o->x_pos_e = z->x_pos_e;
    aux_o->y_pos_s = z->y_pos_s; aux_o->y_pos_e = z->y_pos_e;

    
    for (i = 0; i < on; i++) {
        // fprintf(stderr, "[M::%s::+i->%ld] ovq::[%u, %u), ovt::[%u, %u), hits::[%d, %d)\n", __func__, i, 
        //                             ov->qs, ov->qe, ov->ts, ov->te,
        //                             (ov->qn!=((uint32_t)-1))?(int32_t)ov->qn:-1, (int32_t)ov->tn);
        assert((i<=0)||(ov[i].qs>ov[i-1].qe));
        ovlp_base_direct(z, ch_a, ch_n, &(ov[i]), wl, uref, hpc_g, rref, qstr, tu, exz, aux_o, e_rate, ql, tl, rid);
    }

    int64_t aux_n = aux_o->w_list.n;
    for (i = 0; i < aux_n; i++) {
        if(!(is_ualn_win(aux_o->w_list.a[i]))) continue;
        // if((aux_o->w_list.a[i].x_end+1-aux_o->w_list.a[i].x_start) <= FORCE_CNS_L) {
            // fprintf(stderr, "[aln::-i->%ld::ql->%d] q::[%d, %d), t::[%d, %d), err::%d, clen::%u, mode::%d\n", i, 
            //             aux_o->w_list.a[i].x_end+1-aux_o->w_list.a[i].x_start,
            //             aux_o->w_list.a[i].x_start, aux_o->w_list.a[i].x_end+1, 
            //             aux_o->w_list.a[i].y_start, aux_o->w_list.a[i].y_end+1, 
            //             aux_o->w_list.a[i].error, aux_o->w_list.a[i].clen, aux_o->w_list.a[i].error_threshold);
        // }
        //will overwrite ch_a; does not matter
        rechain_aln(z, cl, aux_o, i, wl, uref, hpc_g, rref, qstr, tu, exz, e_rate, ql, tl, h_khit, rid);
    }
    if(((int64_t)aux_o->w_list.n) > aux_n) {
        for (i = m = 0; i < ((int64_t)aux_o->w_list.n); i++) {
            if((i < aux_n) && (is_ualn_win(aux_o->w_list.a[i]))) continue;
            aux_o->w_list.a[m++] = aux_o->w_list.a[i];
        }
        aux_o->w_list.n = m;
        radix_sort_window_list_xs_srt(aux_o->w_list.a, aux_o->w_list.a+aux_o->w_list.n);
    }

    ///update z by aux_o
    update_overlap_region(z, aux_o, ql, tl);

    // debug_overlap_region(aux_o, qstr, tu, uref, hpc_g, rref);


    // ch_a = cl->list + ch_idx; //update
    // for (i = 0; i < wn; i++) z->w_list.a[i].clen = 0;///clean cigar
    // if(on > 1) {
    //     fprintf(stderr, "[M::%s::] rid::%lu, on::%ld\n", __func__, rid, on);
    // }
    // if(z->y_id == 126) prt_k_mer_hit(ch_a, ch_n);
    // for (i = ch_i = 0; i < on; i++) {
    //     assert((i<=0)||(ov[i].qs > ov[i-1].qe));
    //     ov[i].sec = 16;///do not know the aln type
    //     ch_i = sub_base_aln(z, dp, ch_a, ch_n, pe, ov[i].qs, ov[i].qe, wl, uref, hpc_g, rref, qstr, tu, exz, e_rate, ql, tl, ch_i, rid);
    //     pe = ov[i].qe;
    // }
}



#define gen_err_unaligned(xl, yl) (((xl)<=FORCE_SIN_L)?(MAX((xl), (yl))):MAX((MIN((xl), (yl))), ((xl*0.51)+1)))

#define ovlp_id(x) ((x).tn)
#define ovlp_min_wid(x) ((x).ts)
#define ovlp_max_wid(x) ((x).te)
#define ovlp_cur_wid(x) ((x).qn)
#define ovlp_cur_xoff(x) ((x).qs)
#define ovlp_cur_coff(x) ((x).qe)
#define ovlp_bd(x) ((x).sec)

int64_t retrieve_cigar_err_debug(bit_extz_t *ez, int64_t s, int64_t e)
{
    if(!ez->cigar.n) return 0;
    int64_t err = 0, xk = ez->ts; int64_t ws, we, os, oe, ovlp;
    uint32_t ck = 0, cl; uint16_t op;

    os = MAX(s, ez->ts); oe = MIN(e, ez->te+1);
    ovlp = ((oe>os)? (oe-os):0); 
    if(!ovlp) {
        fprintf(stderr, "[M::%s::] s::%ld, e::%ld, ez->ts::%u, ez->te::%u\n", 
        __func__, s, e, ez->ts, ez->te);
    }
    assert(ovlp);

    //some cigar will span s or e
    while (ck < ez->cigar.n && xk < e) {//[s, e)
        ws = xk; 
        ck = pop_trace(&(ez->cigar), ck, &op, &cl);
        if(op!=2) xk += cl;
        we = xk;
        os = MAX(s, ws); oe = MIN(e, we);
        ovlp = ((oe>os)? (oe-os):0);
        if((op==2) && (ws>=s) && (ws<e)) {
            ovlp = cl;
        }
        if((!ovlp) || (!op)) continue;
        err += ovlp;
    }
    return err;
} 

///[s, e)
int64_t extract_sub_cigar_err_debug(overlap_region *z, int64_t s, int64_t e)
{
    int64_t wk = 0, wn = z->w_list.n, ws, we, os, oe, ovlp; 
    window_list *m; int64_t xl, yl, werr, err; bit_extz_t ez;
    for (wk = err = 0; wk < wn; wk++) {
        m = &(z->w_list.a[wk]); 
        ws = m->x_start; we = m->x_end+1;
        if(ws >= e) break;
        os = MAX(s, ws); oe = MIN(e, we);
        ovlp = ((oe>os)? (oe-os):0);
        if(ovlp) {
            xl = m->x_end+1-m->x_start;
            yl = m->y_end+1-m->y_start;
            if((is_ualn_win((*m))) || (is_est_aln((*m)))) {
                if(is_ualn_win((*m))) { //unmapped
                    werr = gen_err_unaligned(xl, yl);
                } else {
                    werr = m->error;//shared window
                }
                if(ovlp < xl) {
                    werr = (((double)ovlp)/((double)xl))*((double)werr);
                }
                //skip the whole window
                err += werr; 
            } else {
                if(ovlp == xl) {
                    //skip the whole window
                    err += m->error; 
                } else {
                    set_bit_extz_t(ez, (*z), wk); 
                    err += retrieve_cigar_err_debug(&ez, os, oe);
                }                
            }
        }
    }
    return err;
}

int64_t retrieve_cigar_err(bit_extz_t *ez, int64_t s, int64_t e, int64_t *xk, int64_t *ck)
{
    if(!ez->cigar.n) return 0;
    int64_t cn = ez->cigar.n, op, err = 0; int64_t ws, we, os, oe, ovlp;
    if(((*ck) < 0) || ((*ck) > cn)) {//(*ck) == cn is allowed
        (*ck) = 0; (*xk) = ez->ts; 
    }

    while ((*ck) > 0 && (*xk) > s) {
        --(*ck);
        op = ez->cigar.a[(*ck)]>>14;
        if(op!=2) (*xk) -= (ez->cigar.a[(*ck)]&(0x3fff));
    }

    //some cigar will span s or e
    while ((*ck) < cn && (*xk) < e) {//[s, e)
        ws = (*xk); 
        op = ez->cigar.a[(*ck)]>>14;
        if(op!=2) (*xk) += (ez->cigar.a[(*ck)]&(0x3fff));
        we = (*xk);
        os = MAX(s, ws); oe = MIN(e, we);
        ovlp = ((oe>os)? (oe-os):0);
        if((op==2) && (ws>=s) && (ws<e)) {
            ovlp = (ez->cigar.a[(*ck)]&(0x3fff));
        }

        // if(s == 22694 && e == 38018) {
        //     fprintf(stderr, "[M::%s]\tw::[%ld,\t%ld)\tc::%ld\top::%ld\terr::%ld\n", __func__, ws, we, (*ck), op, err + (op?ovlp:0));
        // } 
        (*ck)++;
        if((!ovlp) || (!op)) continue;
        err += ovlp;
    }

    // int64_t debug_err = retrieve_cigar_err_debug(ez, s, e);
    // if(!(err == debug_err)) {
    //     fprintf(stderr, "[M::%s::] err::%ld, debug_err::%ld, s::%ld, e::%ld, ez->ts::%u, ez->te::%u\n", 
    //     __func__, err, debug_err, s, e, ez->ts, ez->te);
    // }
    // assert(err == debug_err);
    return err;
} 
///[s, e)
int64_t extract_sub_cigar_err(overlap_region *z, int64_t s, int64_t e, ul_ov_t *p)
{
    int64_t wk = ovlp_cur_wid(*p), xk = ovlp_cur_xoff(*p), ck = ovlp_cur_coff(*p);
    int64_t min_w = ovlp_min_wid(*p), max_w = ovlp_max_wid(*p);//[min_w, max_w]
    bit_extz_t ez; window_list *m; 
    int64_t ws, we, os, oe, ovlp, err = 0, xl, yl, werr, tot = e - s;
    if(wk < min_w || wk > max_w) wk = min_w;
    for (; wk >= min_w && z->w_list.a[wk].x_start > s; wk--); 
    if(wk < min_w || wk > max_w) return -1;
    for (; wk <= max_w && z->w_list.a[wk].x_end < s; wk++); 
    if(wk < min_w || wk > max_w) return -1;
    //s >= w_list.a[wk].x_start && s <= w_list.a[wk].x_end
    if(wk != ovlp_cur_wid(*p)) {//xk is global, while ck is local
        xk = z->w_list.a[wk].x_start; ck = 0;
    }
    // fprintf(stderr, "[M::%s] wk::%ld, ck::%ld, xk::%ld\n", __func__, wk, ck, xk);
    // fprintf(stderr, "+[M::%s] wk::%ld, ck::%ld, xk::%ld, w::[%d, %d), bound::[%ld, %ld)\n", 
    //     __func__, wk, ck, xk, z->w_list.a[wk].x_start, z->w_list.a[wk].x_end+1, s, e);

    while(wk <= max_w && z->w_list.a[wk].x_start < e) {///[s, e)
        m = &(z->w_list.a[wk]); 
        ws = m->x_start; we = m->x_end+1;
        os = MAX(s, ws); oe = MIN(e, we);
        ovlp = ((oe>os)? (oe-os):0);
        // fprintf(stderr, "-[M::%s] ovlp::%ld, wk::%ld, ck::%ld, xk::%ld, w::[%ld, %ld), bound::[%ld, %ld)\n", 
        // __func__, ovlp, wk, ck, xk, ws, we, s, e);
        if(ovlp) {
            xl = m->x_end+1-m->x_start;
            yl = m->y_end+1-m->y_start;
            if((is_ualn_win((*m))) || (is_est_aln((*m)))) {
                if(is_ualn_win((*m))) { //unmapped
                    werr = gen_err_unaligned(xl, yl);
                } else {
                    werr = m->error;//shared window
                }
                if(ovlp < xl) {
                    werr = (((double)ovlp)/((double)xl))*((double)werr);
                }
                //skip the whole window
                err += werr; xk = m->x_end+1; ck = m->clen;
            } else {
                if(ovlp == xl) {
                    //skip the whole window
                    err += m->error; xk = m->x_end+1; ck = m->clen;
                } else {
                    // if(os == 22694 && oe == 38018) {
                    //     fprintf(stderr, "\n-[M::%s]\tos::%ld\toe::%ld\n", __func__, os, oe);
                    // } 
                    set_bit_extz_t(ez, (*z), wk); 
                    err += retrieve_cigar_err(&ez, os, oe, &xk, &ck);
                }                
            }
        }
        tot -= ovlp; 
        if(xk >= e) break;//[min_w, max_w] && [s, e)
        wk++; if(wk > max_w) break;
        xk = z->w_list.a[wk].x_start; ck = 0;//reset
    }
    assert(!tot);
    ovlp_cur_wid(*p) = wk; ovlp_cur_xoff(*p) = xk; ovlp_cur_coff(*p) = ck;
    return err;
}

#define bst_ov(x) ((x).misBase)
#define ov_dif(x) ((x).cov)
#define ov_id(x) ((x).overlapID)
#define ov_xoff(x) ((x).site)
#define var_id(x) ((x).overlapSite)

#define var_s(x) ((x).site)
#define var_l(x) ((x).overlap_num)
#define var_occ(x) ((x).occ_0)
#define var_min_dif(x) ((x).score)
#define var_min_ovid(x) ((x).id)
#define var_h_idx(x) ((x).occ_1)


uint64_t query_gen_gov_idx(asg64_v *ovidx, uint64_t v, uint64_t w)
{
    uint64_t m, s, e;
    if(v > w) {
        m = v; v = w; w = m;
    }
    s = ovidx->a[v]>>32; e = s + (uint32_t)ovidx->a[v];
    for (m = s; m < e; m++) {
        if(ovidx->a[m] == w) return 1;
    }
    return 0;
}

void push_sec_aln(overlap_region *z, int64_t s, int64_t e, int64_t sec_err)
{
    window_list *p;
    if(z->align_length > 0) {
        p = z->w_list.a + z->w_list.n + z->align_length - 1;
        if((p->x_end == s) && ((!!(p->clen)) == (!!sec_err))) {
            p->x_end = e; p->clen += sec_err;
            return;
        }
    }
    if((z->w_list.n+z->align_length)==z->w_list.m) {
        z->w_list.m = z->w_list.m? z->w_list.m<<1 : 2; 
        z->w_list.a = (window_list*)realloc(z->w_list.a, sizeof(window_list)*z->w_list.m); 
    }
    p = &(z->w_list.a[z->w_list.n+z->align_length]); z->align_length++; memset(p, 0, sizeof((*p)));
    p->x_start = s; p->x_end = e; p->clen += sec_err;
}

void push_sec_aln_robust(overlap_region *z, int64_t s, int64_t e, int64_t sec_err)
{
    window_list *p;
    if(z->align_length > 0) {
        p = z->w_list.a + z->w_list.n + z->align_length - 1;
        if((p->x_end == s) && (p->clen == 0) && ((!!(p->clen)) == (!!sec_err))) {
            p->x_end = e; p->clen += sec_err;
            return;
        }
    }
    if((z->w_list.n+z->align_length)==z->w_list.m) {
        z->w_list.m = z->w_list.m? z->w_list.m<<1 : 2; 
        z->w_list.a = (window_list*)realloc(z->w_list.a, sizeof(window_list)*z->w_list.m); 
    }
    p = &(z->w_list.a[z->w_list.n+z->align_length]); z->align_length++; memset(p, 0, sizeof((*p)));
    p->x_start = s; p->x_end = e; p->clen += sec_err;
}

// #define id_mm ((uint64_t)0x7fffffffffffffff)
#define id_set ((uint64_t)0x8000000000000000)
#define id_get(a) ((uint32_t)(a))
#define err_get(a) (((a)&((uint64_t)0x7fffffffffffffff))>>32)

void reassign_sec_err(overlap_region* ol, asg64_v *ovidx, asg64_v *buf, uint64_t bid0)
{
    uint64_t bn = buf->n, bid = bid0, oid, mid, m, s, e;
    if(buf->a[bid]&id_set) return;
    kv_push(uint64_t, *buf, bid);
    while (buf->n > bn) {
        bid = buf->a[--buf->n];
        if(buf->a[bid]&id_set) continue;
        buf->a[bid]|=id_set;

        oid = id_get(buf->a[bid]);
        s = ovidx->a[oid]>>32; e = (uint32_t)ovidx->a[oid];
        for (m = s; m < e; m++) {
            mid = ol[(uint32_t)ovidx->a[m]].overlapLen;
            if(mid == (uint32_t)-1) continue;
            if(buf->a[mid]&id_set) continue;
            kv_push(uint64_t, *buf, mid);
        }
    }
}
///[s, e)
uint64_t gen_region_phase(overlap_region* ol, uint64_t *id_a, uint64_t id_n, uint64_t s, uint64_t e, uint64_t dp, ul_ov_t *c_idx, asg64_v *buf, asg64_v *ovidx)
{
    if(!id_n) return id_n;
    uint64_t k, m, mn, q[2], buf_n, rm_n, oid; int64_t err, msc, msc_k, msc_n;
    overlap_region *z; ul_ov_t *p; buf->n = 0; kv_resize(uint64_t, *buf, dp);
    for (k = buf_n = rm_n = 0; k < id_n; k++) {
        p = &(c_idx[id_a[k]]);
        q[0] = ol[ovlp_id(*p)].w_list.a[ovlp_min_wid(*p)].x_start;
        q[1] = ol[ovlp_id(*p)].w_list.a[ovlp_max_wid(*p)].x_end+1;
        if(q[0]<=s && q[1]>=e) {
            kv_push(uint64_t, *buf, id_a[k]);
            // buf[buf_n++] = id_a[k];
        }
        if(q[1] < e) rm_n++;
    }
    buf_n = buf->n;
    assert(buf_n == dp);//not right

    if(buf_n > 0) {
        for (k = 0, msc = INT32_MAX, msc_k = -1, msc_n = 0; k < buf_n; k++) {
            p = &(c_idx[(uint32_t)buf->a[k]]); z = &(ol[ovlp_id(*p)]); 
            // fprintf(stderr, "+++[M::%s::utg%.6dl] wid::%u, xoff::%u, coff::%u\n", __func__, 
            //     (int32_t)ol[ovlp_id(*p)].y_id+1, ovlp_cur_wid(*p), ovlp_cur_xoff(*p), ovlp_cur_coff(*p));
            err = extract_sub_cigar_err(z, s, e, p);
            // int64_t debug_err = extract_sub_cigar_err_debug(z, s, e);
            // assert(err == debug_err);
            // fprintf(stderr, "[M::%s::] err::%ld, debug_err::%ld\n", __func__, err, debug_err);

            // fprintf(stderr, "---[M::%s::utg%.6dl] wid::%u, xoff::%u, coff::%u, err::%ld\n", __func__, 
            //     (int32_t)ol[ovlp_id(*p)].y_id+1, ovlp_cur_wid(*p), ovlp_cur_xoff(*p), ovlp_cur_coff(*p), err);
            fprintf(stderr, "---[M::%s::utg%.6dl] xoff::[%lu, %lu), err::%ld\n", __func__, 
                (int32_t)ol[ovlp_id(*p)].y_id+1, s, e, err);
            assert(err >= 0);  
            if(err < msc) {
                msc = err; msc_k = k; msc_n = 1;
            } else if(err == msc) {
                msc_n++;
            }
            buf->a[k] |= (((uint64_t)err)<<32);       
        }

        if(msc_n == 1) {
            p = &(c_idx[(uint32_t)buf->a[msc_k]]);
            z = &(ol[ovlp_id(*p)]); mn = 1;
            if(msc_k != 0) {
                m = buf->a[msc_k]; 
                buf->a[msc_k] = buf->a[0]; 
                buf->a[0] = m;
            }
        } else {
            for (k = mn = 0; k < buf_n && (int64_t)mn < msc_n; k++) {
                p = &(c_idx[(uint32_t)buf->a[k]]); 
                z = &(ol[ovlp_id(*p)]); 
                if((buf->a[k]>>32) == (uint64_t)msc) {
                    if(mn != k) {
                        m = buf->a[k]; 
                        buf->a[k] = buf->a[mn]; 
                        buf->a[mn] = m;
                    }
                    mn++;
                }
            }
        }
        // fprintf(stderr, "[M::%s] buf_n::%ld, msc_n::%ld, mn::%lu\n", __func__, buf_n, msc_n, mn);
        for (k = 0; k < buf_n; k++) {
            oid = ovlp_id((c_idx[(uint32_t)buf->a[k]]));
            buf->a[k] >>= 32; buf->a[k] <<= 32; buf->a[k] |= oid;
            ol[oid].overlapLen = k;
            // if(s == 158482) fprintf(stderr, "k->%ld::oid->%ld[M::%s::utg%.6dl] pos::[%lu, %lu)\n", k, oid, __func__,  
            //                         (int32_t)ol[oid].y_id+1, s, e);
        }

        for (k = 0; k < mn; k++) {///best alignment
            z = &(ol[id_get(buf->a[k])]); 
            reassign_sec_err(ol, ovidx, buf, k);
            push_sec_aln(z, s, e, 0);
        }

        for (k = mn; k < buf_n; k++) {
            z = &(ol[id_get(buf->a[k])]); 
            push_sec_aln(z, s, e, ((buf->a[k]&id_set)?(0):(err_get(buf->a[k])-msc)));
        }

        for (k = 0; k < buf_n; k++) {
            ol[id_get(buf->a[k])].overlapLen = (uint32_t)-1;
        }
    }


    if(rm_n) {
        for (k = m = 0; k < id_n; k++) {
            p = &(c_idx[id_a[k]]);
            q[1] = ol[ovlp_id(*p)].w_list.a[ovlp_max_wid(*p)].x_end+1;
            if(q[1] < e) continue;
            id_a[m++] = id_a[k];
        }
        id_n = m;
    }
    return id_n;
}


///[s, e)
uint64_t gen_region_phase_robust(overlap_region* ol, uint64_t *id_a, uint64_t id_n, uint64_t s, uint64_t e, uint64_t dp, ul_ov_t *c_idx, asg64_v *buf)
{
    if(!id_n) return id_n;
    uint64_t k, m, mn, q[2], buf_n, rm_n, oid; int64_t err, msc, msc_k, msc_n;
    overlap_region *z; ul_ov_t *p; buf->n = 0; kv_resize(uint64_t, *buf, dp);
    for (k = buf_n = rm_n = 0; k < id_n; k++) {
        p = &(c_idx[id_a[k]]);
        q[0] = ol[ovlp_id(*p)].w_list.a[ovlp_min_wid(*p)].x_start;
        q[1] = ol[ovlp_id(*p)].w_list.a[ovlp_max_wid(*p)].x_end+1;
        if(q[0]<=s && q[1]>=e) {
            kv_push(uint64_t, *buf, id_a[k]);
            // buf[buf_n++] = id_a[k];
        }
        if(q[1] < e) rm_n++;
    }
    buf_n = buf->n;
    assert(buf_n == dp);//not right

    if(buf_n > 0) {
        for (k = 0, msc = INT32_MAX, msc_k = -1, msc_n = 0; k < buf_n; k++) {
            p = &(c_idx[(uint32_t)buf->a[k]]); z = &(ol[ovlp_id(*p)]); 
            // fprintf(stderr, "+++[M::%s::utg%.6dl] wid::%u, xoff::%u, coff::%u\n", __func__, 
            //     (int32_t)ol[ovlp_id(*p)].y_id+1, ovlp_cur_wid(*p), ovlp_cur_xoff(*p), ovlp_cur_coff(*p));
            err = extract_sub_cigar_err(z, s, e, p);
            // int64_t debug_err = extract_sub_cigar_err_debug(z, s, e);
            // assert(err == debug_err);
            // fprintf(stderr, "[M::%s::] err::%ld, debug_err::%ld\n", __func__, err, debug_err);

            // fprintf(stderr, "---[M::%s::utg%.6dl] wid::%u, xoff::%u, coff::%u, err::%ld\n", __func__, 
            //     (int32_t)ol[ovlp_id(*p)].y_id+1, ovlp_cur_wid(*p), ovlp_cur_xoff(*p), ovlp_cur_coff(*p), err);
            // fprintf(stderr, "---[M::%s::utg%.6dl] xoff::[%lu, %lu), err::%ld\n", __func__, 
            //     (int32_t)ol[ovlp_id(*p)].y_id+1, s, e, err);
            assert(err >= 0);  
            if(err < msc) {
                msc = err; msc_k = k; msc_n = 1;
            } else if(err == msc) {
                msc_n++;
            }
            buf->a[k] |= (((uint64_t)err)<<32);       
        }

        if(msc_n == 1) {
            p = &(c_idx[(uint32_t)buf->a[msc_k]]);
            z = &(ol[ovlp_id(*p)]); mn = 1;
            if(msc_k != 0) {
                m = buf->a[msc_k]; 
                buf->a[msc_k] = buf->a[0]; 
                buf->a[0] = m;
            }
        } else {
            for (k = mn = 0; k < buf_n && (int64_t)mn < msc_n; k++) {
                p = &(c_idx[(uint32_t)buf->a[k]]); 
                z = &(ol[ovlp_id(*p)]); 
                if((buf->a[k]>>32) == (uint64_t)msc) {
                    if(mn != k) {
                        m = buf->a[k]; 
                        buf->a[k] = buf->a[mn]; 
                        buf->a[mn] = m;
                    }
                    mn++;
                }
            }
        }
        // fprintf(stderr, "[M::%s] buf_n::%ld, msc_n::%ld, mn::%lu\n", __func__, buf_n, msc_n, mn);
        for (k = 0; k < buf_n; k++) {
            oid = ovlp_id((c_idx[(uint32_t)buf->a[k]]));
            buf->a[k] >>= 32; buf->a[k] <<= 32; buf->a[k] |= oid;
            // ol[oid].overlapLen = k;//no need to set 
            // if(s == 158482) fprintf(stderr, "k->%ld::oid->%ld[M::%s::utg%.6dl] pos::[%lu, %lu)\n", k, oid, __func__,  
            //                         (int32_t)ol[oid].y_id+1, s, e);
        }

        for (k = 0; k < mn; k++) {///best alignment
            z = &(ol[id_get(buf->a[k])]); 
            // reassign_sec_err(ol, ovidx, buf, k);
            // push_sec_aln(z, s, e, 0);
            push_sec_aln_robust(z, s, e, 0);
        }

        for (k = mn; k < buf_n; k++) {
            z = &(ol[id_get(buf->a[k])]); 
            // push_sec_aln(z, s, e, ((buf->a[k]&id_set)?(0):(err_get(buf->a[k])-msc)));
            push_sec_aln_robust(z, s, e, (err_get(buf->a[k])-msc));
        }

        // for (k = 0; k < buf_n; k++) {
        //     ol[id_get(buf->a[k])].overlapLen = (uint32_t)-1;
        // }
    }


    if(rm_n) {
        for (k = m = 0; k < id_n; k++) {
            p = &(c_idx[id_a[k]]);
            q[1] = ol[ovlp_id(*p)].w_list.a[ovlp_max_wid(*p)].x_end+1;
            if(q[1] < e) continue;
            id_a[m++] = id_a[k];
        }
        id_n = m;
    }
    return id_n;
}

///[s, e)
int64_t extract_sub_cigar_err_rr(overlap_region *z, int64_t s, int64_t e, ul_ov_t *p)
{
    int64_t wk = ovlp_cur_wid(*p), xk = ovlp_cur_xoff(*p), ck = ovlp_cur_coff(*p);
    int64_t min_w = ovlp_min_wid(*p), max_w = ovlp_max_wid(*p);//[min_w, max_w]
    bit_extz_t ez; window_list *m; int64_t bd = ovlp_bd(*p), s0, e0;
    s0 = ((int64_t)(z->w_list.a[min_w].x_start)) + bd; 
    e0 = ((int64_t)(z->w_list.a[max_w].x_end))+1-bd;
    if(s < s0) s = s0; if(e > e0) e = e0;///exclude boundary
    if(s >= e) return -1;

    int64_t ws, we, os, oe, ovlp, err = 0, xl, yl, werr, tot = e - s;
    if(wk < min_w || wk > max_w) wk = min_w;
    for (; wk >= min_w && z->w_list.a[wk].x_start > s; wk--); 
    if(wk < min_w || wk > max_w) return -1;
    for (; wk <= max_w && z->w_list.a[wk].x_end < s; wk++); 
    if(wk < min_w || wk > max_w) return -1;
    //s >= w_list.a[wk].x_start && s <= w_list.a[wk].x_end
    if(wk != ovlp_cur_wid(*p)) {//xk is global, while ck is local
        xk = z->w_list.a[wk].x_start; ck = 0;
    }

    while(wk <= max_w && z->w_list.a[wk].x_start < e) {///[s, e); [min_w, max_w]
        m = &(z->w_list.a[wk]); 
        ws = m->x_start; we = m->x_end+1;
        os = MAX(s, ws); oe = MIN(e, we);
        ovlp = ((oe>os)? (oe-os):0);
        
        if(ovlp) {
            xl = m->x_end+1-m->x_start;
            yl = m->y_end+1-m->y_start;
            if((is_ualn_win((*m))) || (is_est_aln((*m)))) {
                if(is_ualn_win((*m))) { //unmapped
                    werr = gen_err_unaligned(xl, yl);
                } else {
                    werr = m->error;//shared window
                }
                if(ovlp < xl) {
                    werr = (((double)ovlp)/((double)xl))*((double)werr);
                }
                //skip the whole window
                err += werr; xk = m->x_end+1; ck = m->clen;
            } else {
                if(ovlp == xl) {
                    //skip the whole window
                    err += m->error; xk = m->x_end+1; ck = m->clen;
                } else {
                    set_bit_extz_t(ez, (*z), wk); 
                    err += retrieve_cigar_err(&ez, os, oe, &xk, &ck);
                }                
            }
        }
        tot -= ovlp; 
        if(xk >= e) break;//[min_w, max_w] && [s, e)
        wk++; if(wk > max_w) break;
        xk = z->w_list.a[wk].x_start; ck = 0;//reset
    }
    assert(!tot);
    ovlp_cur_wid(*p) = wk; ovlp_cur_xoff(*p) = xk; ovlp_cur_coff(*p) = ck;
    return err;
}


uint64_t is_mask_ov(mask_ul_ov_t *mk, uint64_t *bes_id, uint64_t bes_n, uint64_t sec_id)
{
    uint64_t s = mk->idx.a[sec_id]>>32, e = (uint32_t)(mk->idx.a[sec_id]), bk, si;
    bk = 0; si = s;
    while (bk < bes_n && si < e) {
        if (bes_id[bk] < mk->srt.a[si].tn) {
            bk++;
        } else if(mk->srt.a[si].tn < bes_id[bk]) {
            si++;
        } else {///bes_id[bk] == mk->srt.a[si].tn
            if(!(mk->srt.a[si].qs)) return 1;
            return 0;
        }
    }
    return 0;
}

///[s, e)
uint64_t gen_region_phase_robust_rr(overlap_region* ol, uint64_t *id_a, uint64_t id_n, uint64_t s, uint64_t e, uint64_t dp, ul_ov_t *c_idx, asg64_v *buf, mask_ul_ov_t *mk)
{
    if(!id_n) return id_n;
    uint64_t k, m, mn, q[2], buf_n, rm_n, oid; int64_t err, msc, msc_k, msc_n;
    overlap_region *z; ul_ov_t *p; buf->n = 0; kv_resize(uint64_t, *buf, dp);
    for (k = buf_n = rm_n = 0; k < id_n; k++) {
        p = &(c_idx[id_a[k]]);
        q[0] = ol[ovlp_id(*p)].w_list.a[ovlp_min_wid(*p)].x_start+ovlp_bd(*p);
        q[1] = ol[ovlp_id(*p)].w_list.a[ovlp_max_wid(*p)].x_end+1-ovlp_bd(*p);
        if(q[0]<=s && q[1]>=e) {
            kv_push(uint64_t, *buf, id_a[k]);
        }
        if(q[1] < e) rm_n++;
    }
    buf_n = buf->n;
    // fprintf(stderr, "[M::%s] buf_n::%lu, dp::%lu\n", __func__, buf_n, dp);
    assert(buf_n == dp);//not right

    if(buf_n > 0) {
        for (k = 0, msc = INT32_MAX, msc_k = -1, msc_n = 0; k < buf_n; k++) {
            p = &(c_idx[(uint32_t)buf->a[k]]); z = &(ol[ovlp_id(*p)]); 
            err = extract_sub_cigar_err_rr(z, s, e, p);
            assert(err >= 0);  
            if(err < msc) {
                msc = err; msc_k = k; msc_n = 1;
            } else if(err == msc) {
                msc_n++;
            }
            buf->a[k] |= (((uint64_t)err)<<32);       
        }

        if(msc_n == 1) {
            p = &(c_idx[(uint32_t)buf->a[msc_k]]);
            z = &(ol[ovlp_id(*p)]); mn = 1;
            if(msc_k != 0) {
                m = buf->a[msc_k]; 
                buf->a[msc_k] = buf->a[0]; 
                buf->a[0] = m;
            }
        } else {
            for (k = mn = 0; k < buf_n && (int64_t)mn < msc_n; k++) {
                p = &(c_idx[(uint32_t)buf->a[k]]); 
                z = &(ol[ovlp_id(*p)]); 
                if((buf->a[k]>>32) == (uint64_t)msc) {
                    if(mn != k) {
                        m = buf->a[k]; 
                        buf->a[k] = buf->a[mn]; 
                        buf->a[mn] = m;
                    }
                    mn++;
                }
            }
        }
        // fprintf(stderr, "[M::%s] buf_n::%ld, msc_n::%ld, mn::%lu\n", __func__, buf_n, msc_n, mn);
        for (k = 0; k < buf_n; k++) {
            oid = ovlp_id((c_idx[(uint32_t)buf->a[k]]));
            buf->a[k] >>= 32; buf->a[k] <<= 32; buf->a[k] |= oid;
            // ol[oid].overlapLen = k;//no need to set 
            // if(s == 158482) fprintf(stderr, "k->%ld::oid->%ld[M::%s::utg%.6dl] pos::[%lu, %lu)\n", k, oid, __func__,  
            //                         (int32_t)ol[oid].y_id+1, s, e);
        }

        for (k = 0; k < mn; k++) {///best alignment
            z = &(ol[id_get(buf->a[k])]); 
            // reassign_sec_err(ol, ovidx, buf, k);
            // push_sec_aln(z, s, e, 0);
            push_sec_aln_robust(z, s, e, 0);
            buf->a[k] = id_get(buf->a[k]);///all equally best overlap pieces
        }
        if(!mk) {
            for (k = mn; k < buf_n; k++) {
                err = err_get(buf->a[k])-msc;
                z = &(ol[id_get(buf->a[k])]); 
                push_sec_aln_robust(z, s, e, err);
            }
        } else {
            if(mn > 1) radix_sort_bc64(buf->a, buf->a + mn);
            for (k = mn; k < buf_n; k++) {
                err = err_get(buf->a[k])-msc;
                if(is_mask_ov(mk, buf->a, mn, id_get(buf->a[k]))) err = 0;
                z = &(ol[id_get(buf->a[k])]); 
                push_sec_aln_robust(z, s, e, err);
            }
        }
    }


    if(rm_n) {
        for (k = m = 0; k < id_n; k++) {
            p = &(c_idx[id_a[k]]);
            q[1] = ol[ovlp_id(*p)].w_list.a[ovlp_max_wid(*p)].x_end+1-ovlp_bd(*p);
            if(q[1] < e) continue;
            id_a[m++] = id_a[k];
        }
        id_n = m;
    }
    return id_n;
}


int64_t infer_rovlp(ul_ov_t *li, ul_ov_t *lj, uc_block_t *bi, uc_block_t *bj, All_reads *ridx, ma_ug_t *ug)
{
	int64_t in, is, ie, irev, iqs, iqe, jn, js, je, jrev, jqs, jqe, ir, jr, ts, te, max_s, min_e, s_shift, e_shift;
	
	if(li) {
		in = ug?ug->u.a[li->tn].len:Get_READ_LENGTH(R_INF, li->tn); 
		is = li->ts; ie = li->te; irev = li->rev; iqs = li->qs; iqe = li->qe;
	} else if(bi) {
		in = ug?ug->u.a[bi->hid].len:Get_READ_LENGTH(R_INF, bi->hid); 
		is = bi->ts; ie = bi->te; irev = bi->rev; iqs = bi->qs; iqe = bi->qe;
	} else {
		return 0;
	}

	if(lj) {
		jn = ug?ug->u.a[lj->tn].len:Get_READ_LENGTH(R_INF, lj->tn); 
		js = lj->ts; je = lj->te; jrev = lj->rev; jqs = lj->qs; jqe = lj->qe;
	} else if(bj) {
		jn = ug?ug->u.a[bj->hid].len:Get_READ_LENGTH(R_INF, bj->hid); 
		js = bj->ts; je = bj->te; jrev = bj->rev; jqs = bj->qs; jqe = bj->qe;
	} else {
		return 0;
	}

	max_s = MAX(iqs, jqs); min_e = MIN(iqe, jqe);
	if(min_e <= max_s) return 0;
	s_shift = get_offset_adjust(max_s - iqs, iqe-iqs, ie-is);
	e_shift = get_offset_adjust(iqe - min_e, iqe-iqs, ie-is);
	if(irev) {
		ts = s_shift; s_shift = e_shift; e_shift = ts;
	}
	is += s_shift; ie-= e_shift;

	// if(li && lj && li->tn == 324 && lj->tn == 319 && li->qs == 63841) {
	// 	fprintf(stderr, "+++in:%ld, is:%ld, ie:%ld, irev:%ld, jn:%ld, js:%ld, je:%ld, jrev:%ld\n", in, is, ie, irev, jn, js, je, jrev);
	// }

	s_shift = get_offset_adjust(max_s - jqs, jqe-jqs, je-js);
    e_shift = get_offset_adjust(jqe - min_e, jqe-jqs, je-js);
    if(jrev) {
        ts = s_shift; s_shift = e_shift; e_shift = ts;
    }
    js += s_shift; je-= e_shift;

	if(irev) {
		ts = in - ie; te = in - is;
		is = ts; ie = te;
	}

	if(jrev) {
		ts = jn - je; te = jn - js;
		js = ts; je = te;
	}
	
	// if(li && lj && li->tn == 324 && lj->tn == 319 && li->qs == 63841) {
	// 	fprintf(stderr, "---in:%ld, is:%ld, ie:%ld, irev:%ld, jn:%ld, js:%ld, je:%ld, jrev:%ld\n", in, is, ie, irev, jn, js, je, jrev);
	// }

	if(is <= js) {
		js -= is; is = 0;
	} else {
		is -= js; js = 0;
	}

	ir = in - ie; jr = jn - je;

	if(ir <= jr){
        ie = in; je += ir;        
    }
    else {
		je = jn; ie += jr;
    }

	ir = ie - is; jr = je - js;
	return MAX(ir, jr);
}

void convert_ul_ov_t(ul_ov_t *des, overlap_region *src, ma_ug_t *ug)
{
	des->qn = (uint32_t)-1; des->qs = src->x_pos_s; des->qe = src->x_pos_e+1;
    des->tn = src->y_id; des->el = 1; des->rev = src->y_pos_strand;
    des->sec = src->non_homopolymer_errors;
	if(des->rev) {
		des->ts = ug->u.a[des->tn].len - (src->y_pos_e+1); 
		des->te = ug->u.a[des->tn].len - src->y_pos_s;
	} else {
		des->ts = src->y_pos_s; 
		des->te = src->y_pos_e+1;
	}
}

uint64_t check_connect_ug(const ul_idx_t *uref, uint32_t v, uint32_t w, int64_t bw, double diff_ec_ul, int64_t dq)
{
	const asg_t *g = uref?uref->ug->g:NULL; int64_t dt = -1;
	uint32_t nv = asg_arc_n(g, v), i; asg_arc_t *av = asg_arc_a(g, v);
	for (i = 0; i < nv; i++) {
		if(av[i].del || av[i].v != w) continue;
		dt = av[i].ol; 
		break;
	}
	if(dt < 0) return 0;
	int64_t diff = (dq>dt? dq-dt:dt-dq), mm = MAX(dq, dt);
	mm *= diff_ec_ul; if(mm < bw) mm = bw;
	if(diff <= mm) return 1;
	return 0;
}

uint64_t check_connect_rg(const ul_idx_t *uref, const ug_opt_t *uopt, uint32_t uv, uint32_t uw, int64_t bw, double diff_ec_ul, int64_t dq)
{
	int64_t dt = -1;
	if(uref->ug->u.a[uv>>1].circ || uref->ug->u.a[uw>>1].circ) return 0;
	uint32_t rv = (uref->ug->u.a[uv>>1].a[(uv&1)?(0):(uref->ug->u.a[uv>>1].n-1)]>>32)^(uv&1);
	uint32_t rw = (uref->ug->u.a[uw>>1].a[(uw&1)?(uref->ug->u.a[uw>>1].n-1):(0)]>>32)^(uw&1);
	ma_hit_t_alloc* src = uopt->sources;
	int64_t min_ovlp = uopt->min_ovlp;
	int64_t max_hang = uopt->max_hang;
	uint64_t z, qn, tn, x = rv>>1; int32_t r = 1; asg_arc_t e;
	for (z = 0; z < src[x].length; z++) {
		qn = Get_qn(src[x].buffer[z]); tn = Get_tn(src[x].buffer[z]);
		if(tn != (rw>>1)) continue;
		r = ma_hit2arc(&(src[x].buffer[z]), Get_READ_LENGTH(R_INF, qn), Get_READ_LENGTH(R_INF, tn), max_hang, asm_opt.max_hang_rate, min_ovlp, &e);
		if(r < 0) continue;
		if((e.ul>>32) != rv || e.v != rw) continue;
		dt = e.ol; 
		break;
	}
	if(dt < 0) return 0;
	int64_t diff = (dq>dt? dq-dt:dt-dq), mm = MAX(dq, dt);
	mm *= diff_ec_ul; if(mm < bw) mm = bw;
	if(diff <= mm) return 1;
	return 0;
}

uint32_t govlp_check(const ul_idx_t *uref, const ug_opt_t *uopt, int64_t bw, double diff_ec_ul, ul_ov_t *li, ul_ov_t *lj)
{
	int64_t qo = infer_rovlp(li, lj, NULL, NULL, /**ridx**/NULL, uref->ug); ///overlap length in query (UL read)
	
	// fprintf(stderr, "+++[M::%s::utg%.6dl->utg%.6dl] qo::%ld\n", __func__, (int32_t)li->tn+1, (int32_t)lj->tn+1, qo);
	if(check_connect_ug(uref, ((li->tn<<1)|li->rev)^1, ((lj->tn<<1)|lj->rev)^1, bw, diff_ec_ul, qo)) return 1;
	// fprintf(stderr, "[M::%s::] check_connect_ug fail\n", __func__);
	if(check_connect_rg(uref, uopt, ((li->tn<<1)|li->rev)^1, ((lj->tn<<1)|lj->rev)^1, bw, diff_ec_ul, qo)) return 1;
	// fprintf(stderr, "[M::%s::] check_connect_rg fail\n", __func__);
	return 0;
}

void gen_gov_idx(overlap_region_alloc* ol, const ul_idx_t *uref, const ug_opt_t *uopt, int64_t bw, double diff_ec_ul, asg64_v* idx)
{
    int64_t on = ol->length, k, i; uint64_t os, oe, ovlp; ul_ov_t p, q, *li, *lj;
    kv_resize(uint64_t, *idx, (uint64_t)on); memset(idx->a, 0, sizeof(*(idx->a))*on);
    for (k = 0, idx->n = on; k < on; k++) {
        convert_ul_ov_t(&p, &(ol->list[k]), uref->ug); p.qn = k; 
        // idx->a[k] = idx->n; idx->a[k] <<= 32;
        for (i = on - 1; i >= 0 && i > k && ol->list[i].x_pos_e >= ol->list[k].x_pos_s; i--) {
            // if(k >= i) continue;
            convert_ul_ov_t(&q, &(ol->list[i]), uref->ug); q.qn = i;

            if(p.qe > q.qe) li = &p, lj = &q;
			else if(p.qe == q.qe && p.qs >= q.qs) li = &p, lj = &q;
			else lj = &p, li = &q;
			os = MAX(li->qs, lj->qs), oe = MIN(li->qe, lj->qe);
			ovlp = ((oe > os)? (oe - os):0);
			if(!ovlp) continue;//no overlap

            if(lj->qs <= li->qs+G_CHAIN_INDEL) {
				if(govlp_check(uref, uopt, bw, diff_ec_ul, li, lj)) {
                    idx->a[k]++; idx->a[i]++; 
                    kv_push(uint64_t, *idx, (((uint64_t)k)<<32)|((uint64_t)i));
                    kv_push(uint64_t, *idx, (((uint64_t)i)<<32)|((uint64_t)k));
                }
			} else if((lj->qe+G_CHAIN_INDEL>=li->qe) && (lj->qs+G_CHAIN_INDEL>=li->qs)) {
				if(govlp_check(uref, uopt, bw, diff_ec_ul, lj, li)) {
                    idx->a[k]++; idx->a[i]++; 
                    kv_push(uint64_t, *idx, (((uint64_t)k)<<32)|((uint64_t)i));
                    kv_push(uint64_t, *idx, (((uint64_t)i)<<32)|((uint64_t)k));
                }
			}
        }
    }
    radix_sort_bc64(idx->a + on, idx->a + idx->n);
    for (k = 0, os = oe = on; k < on; k++) {
        oe = os + idx->a[k];
        idx->a[k] = (os<<32)|oe;
        os = oe;
    }

    // for (k = 0; k < on; k++) {
    //     int64_t s, e;
    //     s = idx->a[k]>>32; e = (uint32_t)idx->a[k];
    //     for (i = s; i < e; i++) {
    //         assert((idx->a[i]>>32) == (uint32_t)k);
    //         fprintf(stderr, "k::%ld[M::%s::utg%.6dl] utg%.6dl\n", k, __func__, 
    //             (int32_t)ol->list[k].y_id+1, (int32_t)ol->list[(uint32_t)idx->a[i]].y_id+1);
    //     }
    // }   
}

void prt_overlap_region_stat(overlap_region *z, int64_t sid)
{
    uint64_t k = 0, aln = 0, ualn = 0, err = 0;
    for (k = 0; k < z->w_list.n; k++) {
        if(is_ualn_win(z->w_list.a[k])) {
            ualn += z->w_list.a[k].x_end+1-z->w_list.a[k].x_start;
        } else {
            aln += z->w_list.a[k].x_end+1-z->w_list.a[k].x_start;
            err += z->w_list.a[k].error;
        }
    }
    fprintf(stderr, "[M::%s::utg%.6dl::%c] sid::%ld, q::[%d, %d), t::[%d, %d), aln::%lu, ualn::%lu, err::%lu\n", 
        __func__, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], sid, z->x_pos_s, z->x_pos_e+1, 
        z->y_pos_s, z->y_pos_e+1, aln, ualn, err);
    
}


void prt_overlap_region_phase_stat(overlap_region *z, int64_t sid)
{
    uint64_t k = 0;
    fprintf(stderr, "[M::%s::utg%.6dl::%c] sid::%ld, q::[%d, %d), t::[%d, %d)\n", __func__, 
        (int32_t)z->y_id+1, "+-"[z->y_pos_strand], sid, z->x_pos_s, z->x_pos_e+1, z->y_pos_s, z->y_pos_e+1);
    for (k = 0; k < z->w_list.n; k++) {
        fprintf(stderr, "[k::%lu] sid::%ld, q::[%d, %d), sec::%u\n", k, sid,
        z->w_list.a[k].x_start, z->w_list.a[k].x_end, z->w_list.a[k].clen);
    }
}

void region_phase(overlap_region_alloc* ol, const ul_idx_t *uref, const ug_opt_t *uopt, kv_ul_ov_t *c_idx, asg64_v* idx, asg64_v* buf, asg64_v* buf1, int64_t ulid)
{
    int64_t on = ol->length, k, i, zwn, q[2], t[2], w[2]; 
    uint64_t m; overlap_region *z; ul_ov_t *cp;
    kv_resize(uint64_t, *idx, (ol->length<<1)); 
    kv_resize(ul_ov_t, *c_idx, ol->length);
    for (k = idx->n = c_idx->n = 0; k < on; k++) {
        z = &(ol->list[k]); zwn = z->w_list.n; 
        // prt_overlap_region_stat(z, ulid);
        // if(/**ulid == 35437 &&**/ z->y_id == 15199 || z->y_id == 31315) {
        //     fprintf(stderr, "+0+[M::%s::utg%.6dl::%c] ulid::%ld, all::%u, non-best::%ld, best::%u\n", __func__, 
        //     (int32_t)z->y_id+1, "+-"[z->y_pos_strand], ulid, z->overlapLen, zwn, z->align_length);
        //     prt_overlap_region_stat(&(ol->list[k]), ulid);
        //     // prt_overlap_region_phase_stat(&(ol->list[k]), ulid);
        // }
        // z->align_length = z->overlapLen = z->x_pos_e+1-z->x_pos_s; 
        z->align_length = 0; z->overlapLen = (uint32_t)-1;
        z->non_homopolymer_errors = 0;
        if(!zwn) continue;
        q[0] = q[1] = t[0] = t[1] = w[0] = w[1] = INT32_MIN;
        for (i = 0; i < zwn; i++) {
            if((z->w_list.a[i].x_start==(q[1]+1)) && ((z->w_list.a[i].y_start==(t[1]+1)))) {
                q[1] = z->w_list.a[i].x_end;
                t[1] = z->w_list.a[i].y_end;
                w[1] = i;
            } else {
                if(q[0] != INT32_MIN) {
                    m = ((uint64_t)q[0])<<1; m <<= 32; 
                    m += c_idx->n; kv_push(uint64_t, *idx, m);
                    m = (((uint64_t)q[1])<<1)+1; m <<= 32; 
                    m += c_idx->n; kv_push(uint64_t, *idx, m);

                    kv_pushp(ul_ov_t, *c_idx, &cp);
                    ovlp_id(*cp) = k; ///ovlp id
                    ovlp_min_wid(*cp) = w[0]; ///beg id of windows
                    ovlp_max_wid(*cp) = w[1]; ///end id of windows
                    ovlp_cur_wid(*cp) = w[0]; ///cur id of windows
                    ovlp_cur_xoff(*cp) = z->w_list.a[w[0]].x_start; ///cur xpos
                    ovlp_cur_coff(*cp) = 0; ///cur cigar off in cur window
                    ovlp_bd(*cp) = 0;
                }

                q[0] = z->w_list.a[i].x_start; q[1] = z->w_list.a[i].x_end;
                t[0] = z->w_list.a[i].y_start; t[1] = z->w_list.a[i].y_end;
                w[0] = i; w[1] = i;
            }
        }
        if(q[0] != INT32_MIN) {
            m = ((uint64_t)q[0])<<1; m <<= 32; 
            m += c_idx->n; kv_push(uint64_t, *idx, m);
            m = (((uint64_t)q[1])<<1)+1; m <<= 32; 
            m += c_idx->n; kv_push(uint64_t, *idx, m);

            kv_pushp(ul_ov_t, *c_idx, &cp);
            ovlp_id(*cp) = k; ///ovlp id
            ovlp_min_wid(*cp) = w[0]; ///beg id of windows
            ovlp_max_wid(*cp) = w[1]; ///end id of windows
            ovlp_cur_wid(*cp) = w[0]; ///cur id of windows
            ovlp_cur_xoff(*cp) = z->w_list.a[w[0]].x_start; ///cur xpos
            ovlp_cur_coff(*cp) = 0; ///cur cigar off in cur window
            ovlp_bd(*cp) = 0;
        }
        // if(ulid == 35437 && z->y_id == 109111) {
        //     fprintf(stderr, "+1+[M::%s::utg%.6dl::%c] ulid::%ld, all::%u, non-best::%ld, best::%u\n", __func__, 
        //     (int32_t)z->y_id+1, "+-"[z->y_pos_strand], ulid, z->overlapLen, zwn, z->align_length);
        //     prt_overlap_region_stat(&(ol->list[k]));
        //     prt_overlap_region_phase_stat(&(ol->list[k]));
        // }
    }
    radix_sort_bc64(idx->a, idx->a+idx->n);
    //this is used with gen_region_phase, now give up
    //gen_gov_idx(ol, uref, uopt, G_CHAIN_BW, N_GCHAIN_RATE, buf1);
    // for (m = 0; m < c_idx->n; m++) {
    //     fprintf(stderr, "+++[M::%s::utg%.6dl] q[%d, %d), t[%d, %d), wn::%d\n", __func__, 
    //     (int32_t)ol->list[ovlp_id(c_idx->a[m])].y_id+1,
    //     ol->list[ovlp_id(c_idx->a[m])].w_list.a[ovlp_min_wid(c_idx->a[m])].x_start,
    //     ol->list[ovlp_id(c_idx->a[m])].w_list.a[ovlp_max_wid(c_idx->a[m])].x_end+1,
    //     ol->list[ovlp_id(c_idx->a[m])].w_list.a[ovlp_min_wid(c_idx->a[m])].y_start,
    //     ol->list[ovlp_id(c_idx->a[m])].w_list.a[ovlp_max_wid(c_idx->a[m])].y_end+1,
    //     ovlp_max_wid(c_idx->a[m])+1-ovlp_min_wid(c_idx->a[m]));
    // }
    

    int64_t srt_n = idx->n, dp, old_dp, beg, end;
    for (i = k = 0, dp = old_dp = 0, beg = 0, end = -1; i < srt_n; ++i) {///[beg, end) but coordinates in idx is [, ]
        ///if idx->a.a[] is qe
        old_dp = dp;
        if ((idx->a[i]>>32)&1) {
            --dp; end = (idx->a[i]>>33)+1;
        }else {
            //meet a new overlap; the overlaps are pushed by the x_pos_s
            ++dp; end = (idx->a[i]>>33);
            kv_push(uint64_t, *idx, ((uint32_t)idx->a[i]));
        }
        // fprintf(stderr, "\n[M::%s::] beg::%ld, end::%ld, old_dp::%ld\n", __func__, beg, end, old_dp);
        if((end > beg) && (old_dp >= 2)) {
            // fprintf(stderr, "\n[M::%s::] beg::%ld, end::%ld, old_dp::%ld\n", __func__, beg, end, old_dp);
            // kv_resize(uint64_t, *buf, ((uint32_t)old_dp)<<1);
            // idx->n = srt_n + gen_region_phase(ol->list, idx->a+srt_n, idx->n-srt_n, beg, end, old_dp, c_idx->a, buf, buf1);   
            idx->n = srt_n + gen_region_phase_robust(ol->list, idx->a+srt_n, idx->n-srt_n, beg, end, old_dp, c_idx->a, buf);
        }
        beg = end;
    }
    ///hap->length

    for (k = 0; k < on; k++) {
        z = &(ol->list[k]); 
        // if(ulid == 35437 && z->y_id == 109111) {
        //     fprintf(stderr, "+2+[M::%s::utg%.6dl::%c] ulid::%ld, all::%u, non-best::%ld, best::%u\n", __func__, 
        //     (int32_t)z->y_id+1, "+-"[z->y_pos_strand], ulid, z->overlapLen, zwn, z->align_length);
        //     prt_overlap_region_stat(&(ol->list[k]));
        //     prt_overlap_region_phase_stat(&(ol->list[k]));
        // }



        z->overlapLen = z->x_pos_e+1-z->x_pos_s; 
        z->non_homopolymer_errors = 0; zwn = 0;
        for (i = m = 0; i < z->align_length; i++) {
            z->w_list.a[m] = z->w_list.a[z->w_list.n+i];
            if(z->w_list.a[m].clen > 0) {
                z->non_homopolymer_errors += z->w_list.a[m].clen;
                zwn += z->w_list.a[m].x_end-z->w_list.a[m].x_start;
            } 
            m++;
        }
        z->w_list.n = m; 
        // if(ulid == 35437 && z->y_id == 109111/**zwn > z->overlapLen**/) {
        //     fprintf(stderr, "+3+[M::%s::utg%.6dl::%c] ulid::%ld, all::%u, non-best::%ld, best::%u\n", __func__, 
        //     (int32_t)z->y_id+1, "+-"[z->y_pos_strand], ulid, z->overlapLen, zwn, z->align_length);
        //     prt_overlap_region_stat(&(ol->list[k]));
        //     prt_overlap_region_phase_stat(&(ol->list[k]));
        // }
        assert(zwn <= z->overlapLen);
        z->align_length = z->overlapLen - zwn;
        // /**if(ulid == 35437 && z->y_id == 15199 || z->y_id == 31315)**/ {
        //     prt_overlap_region_phase_stat(z, ulid);
        //     fprintf(stderr, "[M::%s::utg%.6dl::%c] all::%u, non-best::%ld, best::%u\n\n", __func__, 
        //     (int32_t)z->y_id+1, "+-"[z->y_pos_strand], z->overlapLen, zwn, z->align_length);
        // }
        // prt_overlap_region_phase_stat(&(ol->list[k]));
        // z = &(ol->list[k]); 
        // if(z->align_length == z->overlapLen) {///prefer alignments without any trans hit
        //     z->align_length = z->overlapLen = z->x_pos_e+1-z->x_pos_s;
        // }
    }
}


void rphase_rr(overlap_region_alloc* ol, const ul_idx_t *uref, const ug_opt_t *uopt, kv_ul_ov_t *c_idx, asg64_v* idx, asg64_v* buf, int64_t ulid, int64_t bd, mask_ul_ov_t *mk)
{
    int64_t on = ol->length, k, i, zwn, q[2], t[2], w[2]; 
    uint64_t m; overlap_region *z; ul_ov_t *cp;
    kv_resize(uint64_t, *idx, (ol->length<<1)); 
    kv_resize(ul_ov_t, *c_idx, ol->length);
    for (k = idx->n = c_idx->n = 0; k < on; k++) {
        z = &(ol->list[k]); zwn = z->w_list.n; 
        z->align_length = 0; z->overlapLen = (uint32_t)-1;
        z->non_homopolymer_errors = 0;
        if(!zwn) continue;
        q[0] = q[1] = t[0] = t[1] = w[0] = w[1] = INT32_MIN;
        for (i = 0; i < zwn; i++) {
            if((z->w_list.a[i].x_start==(q[1]+1)) && ((z->w_list.a[i].y_start==(t[1]+1)))) {
                q[1] = z->w_list.a[i].x_end;
                t[1] = z->w_list.a[i].y_end;
                w[1] = i;
            } else {
                if(q[0] != INT32_MIN) {
                    q[0] += bd; q[1] -= bd;
                    if(q[1] >= q[0]) {
                        m = ((uint64_t)q[0])<<1; m <<= 32; 
                        m += c_idx->n; kv_push(uint64_t, *idx, m);
                        m = (((uint64_t)q[1])<<1)+1; m <<= 32; 
                        m += c_idx->n; kv_push(uint64_t, *idx, m);

                        kv_pushp(ul_ov_t, *c_idx, &cp);
                        ovlp_id(*cp) = k; ///ovlp id
                        ovlp_min_wid(*cp) = w[0]; ///beg id of windows
                        ovlp_max_wid(*cp) = w[1]; ///end id of windows
                        ovlp_cur_wid(*cp) = w[0]; ///cur id of windows
                        ovlp_cur_xoff(*cp) = z->w_list.a[w[0]].x_start; ///cur xpos
                        ovlp_cur_coff(*cp) = 0; ///cur cigar off in cur window
                        ovlp_bd(*cp) = bd;
                    }
                }

                q[0] = z->w_list.a[i].x_start; q[1] = z->w_list.a[i].x_end;
                t[0] = z->w_list.a[i].y_start; t[1] = z->w_list.a[i].y_end;
                w[0] = i; w[1] = i;
            }
        }
        if(q[0] != INT32_MIN) {
            q[0] += bd; q[1] -= bd;
            if(q[1] >= q[0]) {
                m = ((uint64_t)q[0])<<1; m <<= 32; 
                m += c_idx->n; kv_push(uint64_t, *idx, m);
                m = (((uint64_t)q[1])<<1)+1; m <<= 32; 
                m += c_idx->n; kv_push(uint64_t, *idx, m);

                kv_pushp(ul_ov_t, *c_idx, &cp);
                ovlp_id(*cp) = k; ///ovlp id
                ovlp_min_wid(*cp) = w[0]; ///beg id of windows
                ovlp_max_wid(*cp) = w[1]; ///end id of windows
                ovlp_cur_wid(*cp) = w[0]; ///cur id of windows
                ovlp_cur_xoff(*cp) = z->w_list.a[w[0]].x_start; ///cur xpos
                ovlp_cur_coff(*cp) = 0; ///cur cigar off in cur window
                ovlp_bd(*cp) = bd;
            }
        }
    }
    radix_sort_bc64(idx->a, idx->a+idx->n);

    int64_t srt_n = idx->n, dp, old_dp, beg, end;
    for (i = k = 0, dp = old_dp = 0, beg = 0, end = -1; i < srt_n; ++i) {///[beg, end) but coordinates in idx is [, ]
        ///if idx->a.a[] is qe
        old_dp = dp;
        if ((idx->a[i]>>32)&1) {
            --dp; end = (idx->a[i]>>33)+1;
        }else {
            //meet a new overlap; the overlaps are pushed by the x_pos_s
            ++dp; end = (idx->a[i]>>33);
            kv_push(uint64_t, *idx, ((uint32_t)idx->a[i]));
        }
        if((end > beg) && (old_dp >= 2)) {
            idx->n = srt_n + gen_region_phase_robust_rr(ol->list, idx->a+srt_n, idx->n-srt_n, beg, end, old_dp, c_idx->a, buf, mk);
        }
        beg = end;
    }
}

int64_t gen_aln_ul_ov_t(int64_t in_id, int64_t tl, overlap_region *in, ul_ov_t *ou)
{
    int64_t zwn; memset(ou, 0, sizeof((*ou)));
    zwn = in->w_list.n; if(!zwn) return 0;

    ///in_id -> ovlp id
    ou->qn = in_id; ou->tn = in->y_id; 
    ou->rev = in->y_pos_strand;
    ou->qs = in->w_list.a[0].x_start;
    ou->qe = in->w_list.a[zwn-1].x_end+1;
    
    if(!(ou->rev)) {
        ou->ts = in->w_list.a[0].y_start;
        ou->te = in->w_list.a[zwn-1].y_end+1;
    } else {
        ou->ts = tl-(in->w_list.a[zwn-1].y_end+1);
        ou->te = tl-in->w_list.a[0].y_start;
    }
    return 1;
}

uint32_t push_emask_flt(kv_emask_t *in, uint64_t *flt, uint64_t flt_n, uint64_t qn, kv_ul_ov_t *res)
{
    uint32_t k = 0, z = 0, tn, rn0 = res->n; ul_ov_t *p;
    while (k < in->n && z < flt_n) {
        if (in->a[k].tn < (flt[z]>>32)) {
            k++;
        } else if((flt[z]>>32) < in->a[k].tn) {
            z++;
        } else {///in->a[k].tn == (flt[z]>>32)
            for (tn = in->a[k].tn; k < in->n && in->a[k].tn == tn; k++) {
                kv_pushp(ul_ov_t, *res, &p);
                p->qn = qn; p->qs = in->a[k].qs; p->qe = in->a[k].qe; 
                p->tn = in->a[k].tn; p->ts = in->a[k].ts; p->te = in->a[k].te; 
                p->rev = in->a[k].rev; p->el = in->a[k].full; p->sec = k;
            }
            for (; z < flt_n && (flt[z]>>32) == tn; z++);
        }
    }
    return res->n-rn0;
}

#define is_ul_ov_pe(z, zn, i, mm) ((((i)+1)<(zn))&&((z)[(i)+1].qn==(z)[(i)].qn)&&((z)[(i)+1].tn==(z)[(i)].tn)&&((z)[(i)+1].rev==(z)[(i)].rev)\
            &&((z)[(i)+1].sec==(z)[(i)].sec+1)&&(((z)[(i)].sec)<(mm).a[(z)[(i)].qn].n)&&(((z)[(i)+1].sec)<(mm).a[(z)[(i)+1].qn].n)\
            &&((mm).a[(z)[(i)].qn].a[((z)[(i)].sec)].pe)&&((mm).a[(z)[(i)+1].qn].a[((z)[(i)+1].sec)].pe)\
            &&((mm).a[(z)[(i)].qn].a[((z)[(i)].sec)].dir==0)&&((mm).a[(z)[(i)+1].qn].a[((z)[(i)+1].sec)].dir==1))

#define is_exact_ov(z, mm)  (((((z).sec)<(mm).a[(z).qn].n))&&(!((mm).a[(z).qn].a[((z).sec)].pe))\
            &&((mm).a[(z).qn].a[((z).sec)].el==(uint32_t)-1))


///ref_n <= 2
uint64_t is_cover_ul_ov_t(ma_ug_t *ug, double diff, ul_ov_t *ref, uint64_t ref_n, ul_ov_t *in)
{
    if((!check_ul_ov_t_consist(in, &(ref[0]), ug->g->seq[in->qn].len, ug->g->seq[in->tn].len, diff))) return 0;
    if((ref_n > 1) && (!check_ul_ov_t_consist(in, &(ref[ref_n-1]), ug->g->seq[in->qn].len, ug->g->seq[in->tn].len, 0.06))) return 0;
    uint64_t i; int64_t sql, stl, os, oe, ovlp;
    sql = in->qe - in->qs; stl = in->te - in->ts;   
    for (i = 0; i < ref_n && sql > 0 && stl > 0; i++) {
        os = MAX(in->qs, ref[i].qs); oe = MIN(in->qe, ref[i].qe);
		ovlp = ((oe>os)?(oe-os):(0)); sql -= ovlp;

        os = MAX(in->ts, ref[i].ts); oe = MIN(in->te, ref[i].te);
        ovlp = ((oe>os)?(oe-os):(0)); stl -= ovlp;
    }
    if((sql > 256) && (sql > ((in->qe - in->qs)*0.06))) return 0;
    if((sql <= 256) && (sql > ((in->qe - in->qs)*0.6))) return 0;
    if((stl > 256) && (stl > ((in->te - in->ts)*0.06))) return 0;
    if((stl <= 256) && (stl > ((in->te - in->ts)*0.6))) return 0;
    return 1;
}

uint64_t dedup_src_shared(ma_ug_t *ug, ul_ov_t *ta, uint64_t tn, ul_ov_t *qa, uint64_t qn, idx_emask_t *mm)
{
    uint32_t ti, qi, ts, te, id, k, l;
    ti = qi = 0;
    while (ti < tn && qi < qn) {
        if (ta[ti].tn < qa[qi].tn) {
            ti++;
        } else if(qa[qi].tn < ta[ti].tn) {
            qi++;
        } else {///ta[ti].tn == qa[qi].tn
            id = ta[ti].tn;
            for (ts = ti; ti < tn && ta[ti].tn == id; ti++); te = ti;
            for (; qi < qn && qa[qi].tn == id; qi++) {
                for (k = ts; k < te; k++) {
                    l = 1; 
                    if(is_ul_ov_pe(ta, te, k, (*mm))) {
                        l++; k++;
                    }
                    if(is_cover_ul_ov_t(ug, 0.06, ta+k+1-l, l, &(qa[qi]))) break;
                }
                if(k < te) qa[qi].qn = qa[qi].tn = (uint32_t)-1;
            }
        }
    }
    for (qi = l = 0; qi < qn; qi++) {
        if(qa[qi].tn == (uint32_t)-1) continue;
        qa[l++] = qa[qi];
    }
    return l;
}

void dedup_src_shared1(ma_ug_t *ug, kv_ul_ov_t *res, uint64_t tocc, uint64_t qocc, idx_emask_t *mm)
{
    if((!tocc) || (!qocc)) return;
    uint32_t ti, qi, tn, qn, rn = res->n, id, ts, te, k, l, ff = 0;
    // qi = res->n - qocc; fi = 0; m = qi; qocc = 0;
    // while (qi < res->n && fi < flt_n) {
    //     if(res->a[qi].tn < (flt[fi]>>32)) {
    //         qi++;
    //     } else if((flt[fi]>>32) < res->a[qi].tn) {
    //         fi++;
    //     } else {///res->a[qi].tn == (flt[fi]>>32)
    //         res->a[m++] = res->a[qi]; qocc++; qi++; fi++;
    //     }
    // }
    // res->n = m;
    // if((!tocc) || (!qocc)) return;

    rn = res->n;
    ti = res->n - tocc - qocc; qi = res->n - qocc;//t first, and then q 
    tn = ti + tocc; qn = qi + qocc; kv_resize(ul_ov_t, *res, (rn+tn));///the buf size should be at least tn
    while (ti < tn && qi < qn) {
        if (res->a[ti].tn < res->a[qi].tn) {
            kv_push(ul_ov_t, *res, res->a[ti]); ti++; 
        } else if(res->a[qi].tn < res->a[ti].tn) {
            kv_push(ul_ov_t, *res, res->a[qi]); qi++; 
        } else {///ta[ti].tn == qa[qi].tn
            id = res->a[ti].tn; ts = ti;
            for (; ti < tn && res->a[ti].tn == id; ti++) {
                kv_push(ul_ov_t, *res, res->a[ti]);
            }
            te = ti;
            
            for (; qi < qn && res->a[qi].tn == id; qi++) {
                for (k = ts; k < te; k++) {
                    l = k; 
                    if(is_ul_ov_pe(res->a, te, k, (*mm))) k++;
                    if(is_cover_ul_ov_t(ug, 0.06, res->a+l, k+1-l, &(res->a[qi]))) {
                        ff = 1; break;
                    }
                }
                if(k >= te) kv_push(ul_ov_t, *res, res->a[qi]);
            }
        }
    }
    while (ti < tn) {
        kv_push(ul_ov_t, *res, res->a[ti]); ti++; 
    }
    while (qi < qn) {
        kv_push(ul_ov_t, *res, res->a[qi]); qi++; 
    }
    
    // fprintf(stderr, "[M::%s]\tff::%u\n", __func__, ff);
    if(!ff) {///nothing has been removed as duplication
        res->n = rn;
    } else if(res->n > rn) {
        l = res->n; res->n = rn - tocc - qocc;
        for (k = rn; k < l; k++) res->a[res->n++] = res->a[k];
    }
}

#define aln2ov(in, ou, tl) \
        {(ou).qn=(in).x_id;(ou).tn=(in).y_id;(ou).rev=(in).y_pos_strand;\
        (ou).qs=(in).x_pos_s;(ou).qe=(in).x_pos_e+1;\
        if((in).y_pos_strand){(ou).ts=(tl)-(in).y_pos_e-1;(ou).te=(tl)-(in).y_pos_s;}\
        else {(ou).ts=(in).y_pos_s;(ou).te=(in).y_pos_e+1;}}

uint64_t cal_x_ul_ovlp(ma_ug_t *ug, overlap_region *q, overlap_region *t, ul_ov_t *res)
{
    ul_ov_t qo, to; uint32_t ql, tl, os, oe, rqs, rqe, rts, rte;
    ql = ug->g->seq[q->y_id].len; tl = ug->g->seq[t->y_id].len;
    aln2ov((*q), qo, ql); aln2ov((*t), to, tl); memset(res, 0, sizeof(*res));
    os = MAX(qo.qs, to.qs); oe = MIN(qo.qe, to.qe);
    if(oe <= os) return 0;
    if(infer_se(qo.qs, qo.qe, qo.ts, qo.te, qo.rev, os, oe, &rqs, &rqe) &&
			infer_se(to.qs, to.qe, to.ts, to.te, to.rev, os, oe, &rts, &rte)) {
        if(rqe > rqs && rte > rts) {
            res->qn = qo.tn; res->tn = to.tn; res->rev = ((qo.rev==to.rev)?0:1);
            res->qs = rqs; res->qe = rqe; res->ts = rts; res->te = rte; 
            return 1;
        }
    }
    return 0;
}

uint64_t check_mask_exist(ma_ug_t *ug, overlap_region *q, overlap_region *t, kv_ul_ov_t *ov_db, idx_emask_t *mm, double len_diff, uint64_t *is_exact, uint64_t *ovdb_idx)
{
    ul_ov_t rr, *a; uint64_t k, eid, an, ss; (*is_exact) = 0; (*ovdb_idx) = (uint64_t)-1;
    if(!cal_x_ul_ovlp(ug, q, t, &rr)) return 0;
    ss = q->overlapLen;
    if(q->x_pos_strand) {///if q does not have mask
        eid = rr.tn; a = ov_db->a + ss; an = q->x_pos_strand;
        for (k = 0; k < an && a[k].tn < eid; k++);
        for (; k < an && a[k].tn == eid; k++) {
            if((!check_ul_ov_t_consist(&rr, &(a[k]), ug->g->seq[rr.qn].len, ug->g->seq[rr.tn].len, len_diff))) continue;
            if(is_ul_ov_pe(a, an, k, (*mm))) {
                k++; if((!check_ul_ov_t_consist(&rr, &(a[k]), ug->g->seq[rr.qn].len, ug->g->seq[rr.tn].len, len_diff))) continue;
            }
            if(is_exact_ov(a[k],(*mm))) (*is_exact) = 1;
            (*ovdb_idx) = ss + k;
            return 1;
        }
    }

    k = rr.qn; rr.qn = rr.tn; rr.tn = k;
    k = rr.qs; rr.qs = rr.ts; rr.ts = k;
    k = rr.qe; rr.qe = rr.te; rr.te = k;
    ss = t->overlapLen;
    if(t->x_pos_strand) {///if t does not have mask
        eid = rr.tn; a = ov_db->a + ss; an = t->x_pos_strand;
        for (k = 0; k < an && a[k].tn < eid; k++);
        for (; k < an && a[k].tn == eid; k++) {
            if((!check_ul_ov_t_consist(&rr, &(a[k]), ug->g->seq[rr.qn].len, ug->g->seq[rr.tn].len, len_diff))) continue;
            if(is_ul_ov_pe(a, an, k, (*mm))) {
                k++; if((!check_ul_ov_t_consist(&rr, &(a[k]), ug->g->seq[rr.qn].len, ug->g->seq[rr.tn].len, len_diff))) continue;
            }
            if(is_exact_ov(a[k],(*mm))) (*is_exact) = 1;
            (*ovdb_idx) = ss + k;
            return 1;
        }
    }
    return 0;
}

void push_consist_ul_ovlps(ma_ug_t *ug, kv_ul_ov_t *ov_db, uint64_t s, uint64_t e, double len_diff, ul_ov_t *ref, idx_emask_t *mm, kv_ul_ov_t *out, uint64_t flip_res, uint64_t *is_exact)
{
    (*is_exact) = 0;
    if(s >= e) return;
    uint64_t k, z, ql = ug->g->seq[ref->qn].len, tl = ug->g->seq[ref->tn].len; ul_ov_t *p;
    for (k = s; k < e && ov_db->a[k].tn < ref->tn; k++);
    for (; k < e && ov_db->a[k].tn == ref->tn; k++) {
        if((!check_ul_ov_t_consist(ref, &(ov_db->a[k]), ql, tl, len_diff))) continue;
        z = k;
        if(is_ul_ov_pe(ov_db->a, ov_db->n, k, (*mm))) {
            k++; if((!check_ul_ov_t_consist(ref, &(ov_db->a[k]), ql, tl, len_diff))) continue;
        }
        if(is_exact_ov(ov_db->a[k], (*mm))) {
            (*is_exact) = 1; return;
        }

        for (; z <= k; z++) {
            kv_pushp(ul_ov_t, *out, &p);
            *p = ov_db->a[k];
            if(flip_res) {
                p->qn = ov_db->a[k].tn; p->tn = ov_db->a[k].qn;
                p->qs = ov_db->a[k].ts; p->qe = ov_db->a[k].te;
                p->ts = ov_db->a[k].qs; p->te = ov_db->a[k].qe;
            }
        }
    }
}

uint64_t cal_no_bd_coor(ul_ov_t *in, ul_ov_t *ou, idx_emask_t *mm, int64_t bd)
{
    int64_t qs, qe, ts, te; emask_t *z;
    qs = in->qs; qe = in->qe; ts = in->ts; te = in->te;
    if(in->sec != ((uint32_t)(0x3fffffff))) {///might be an end mask; no need shrink
        if((in->sec < mm->a[in->qn].n) && (in->tn == mm->a[in->qn].a[in->sec].tn) && (in->rev == mm->a[in->qn].a[in->sec].rev) && 
            (in->qs == mm->a[in->qn].a[in->sec].qs) && (in->qe == mm->a[in->qn].a[in->sec].qe) && 
                (in->ts == mm->a[in->qn].a[in->sec].ts) && (in->te == mm->a[in->qn].a[in->sec].te)) {
            z = &(mm->a[in->qn].a[in->sec]); assert(z->el != ((uint32_t)-1));
            if(z->dir == 0) {
                qe -= (z->el>>1); 
                if(in->rev) ts += (z->el>>1);
                else te -= (z->el>>1);
            } else {
                qs += (z->el>>1); 
                if(in->rev) te -= (z->el>>1);
                else ts += (z->el>>1);
            }
            // fprintf(stderr, "+[M::%s]\tz->dir::%u\tz->el::%u\tq::[%ld,\t%ld)\tt::[%ld,\t%ld)\n", __func__, z->dir, z->el, qs, qe, ts, te); 
        } else {
            assert((in->sec < mm->a[in->tn].n) && (in->qn == mm->a[in->tn].a[in->sec].tn) && (in->rev == mm->a[in->tn].a[in->sec].rev) && 
                (in->qs == mm->a[in->tn].a[in->sec].ts) && (in->qe == mm->a[in->tn].a[in->sec].te) && 
                (in->ts == mm->a[in->tn].a[in->sec].qs) && (in->te == mm->a[in->tn].a[in->sec].qe));
            z = &(mm->a[in->tn].a[in->sec]); assert(z->el != ((uint32_t)-1));
            if(z->dir == 1) {
                qe -= (z->el>>1); 
                if(in->rev) ts += (z->el>>1);
                else te -= (z->el>>1);
            } else {
                qs += (z->el>>1); 
                if(in->rev) te -= (z->el>>1);
                else ts += (z->el>>1);
            }
            // fprintf(stderr, "-[M::%s]\tz->dir::%u\tz->el::%u\tq::[%ld,\t%ld)\tt::[%ld,\t%ld)\n", __func__, z->dir, z->el, qs, qe, ts, te); 
        }
    } else {///shrink anyway
        qs += bd; qe -= bd; ts += bd; te -= bd;
    }
    if(qe > qs && te > ts) {
        (*ou) = (*in);
        ou->qs = qs; ou->qe = qe; ou->ts = ts; ou->te = te;
        return 1;
    }
    return 0;
}

typedef struct {
	int64_t qi, qs, qe;
	int64_t ti;
    int64_t qis, qie, tis, tie;
    int64_t wi, wn;
	int64_t tot, lc, cql, ctl, ci, lerr;
    overlap_region *z;
} pe_cigar_iter_t;

#define citer_end(m)  ((m).qi>=(m).qe)

/**
inline void pop_citer(pe_cigar_iter_t *it)
{
    it->lc = it->lerr = -1; it->cql = it->ctl = 0; 
    it->qis = it->qie; it->tis = it->tie;
    if(it->wi >= it->wn) return;
    int64_t ws, we, os, oe, ovlp, ql, tl, werr; 
    window_list *m; bit_extz_t ez;

    while (it->wi < it->wn) {
        m = &(it->z->w_list.a[it->wi]); 
        ws = m->x_start; we = m->x_end+1;
        os = MAX(it->qs, ws); oe = MIN(it->qe, we);
        ovlp = ((oe>os)? (oe-os):0);

        if(ovlp) {
            ql = m->x_end+1-m->x_start;
            tl = m->y_end+1-m->y_start;
            if((is_ualn_win((*m))) || (is_est_aln((*m)))) {
                it->lc = 4;//not an ordinary cigar
                if(is_ualn_win((*m))) { //unmapped
                    werr = gen_err_unaligned(ql, tl);
                } else {
                    werr = m->error;//shared window
                    if(!werr) it->lc = 0;///treat it as match
                }
                it->lerr = werr;
                if(ovlp < ql) {
                    werr = (((double)ovlp)/((double)ql))*((double)werr);
                }
                //skip the whole window
                it->tot += werr; 
                it->qis = it->qi; it->qi = it->qie = m->x_end+1; it->cql = it->qie - it->qis;
                it->tis = it->ti; it->ti = it->tie = m->x_start; it->ctl = it->tie - it->tis;
                it->ci = 0; it->wi++;
                return;
            } else {
                if(it->qi < m->x_start || it->ti < m->y_start) {
                    it->qi = m->x_start; it->ti = m->y_start; it->ci = 0;
                }
                set_bit_extz_t(ez, (*(it->z)), it->wi);
                // if(ovlp == ql) {
                //     //skip the whole window
                //     err += m->error; xk = m->x_end+1; ck = m->clen;
                // } else {
                //     set_bit_extz_t(ez, (*z), wk); 
                //     err += retrieve_cigar_err(&ez, os, oe, &xk, &ck);
                // }                
            }
        }
    }
}
**/

int64_t get_pe_diff(overlap_region *q, uint64_t *qmask, uint64_t qmask_n, overlap_region *t, uint64_t *tmask, uint64_t tmask_n, int64_t bd)
{
    int64_t s, e; overlap_region *z;
    pe_cigar_iter_t qi, ti;
    z = q; if(!(z->w_list.n)) return 0;
    s = z->w_list.a[0].x_start; e = z->w_list.a[z->w_list.n-1].x_end+1;
    s += bd; e -= bd; if(s >= e) return 0;
    memset(&qi, 0, sizeof(qi));
    qi.z = q; qi.qs = s; qi.qe = e; 
    qi.qi = q->w_list.a[0].x_start;
    qi.ti = q->w_list.a[0].y_start;
    qi.wi = 0; qi.wn = q->w_list.n;
    qi.tot = 0; qi.ci = 0; qi.lc = -1; qi.cql = qi.ctl = 0;
    qi.qis = qi.qie = qi.qi; 
    qi.tis = qi.tie = qi.ti; 
    qi.lerr = -1;

    z = t; if(!(z->w_list.n)) return 0;
    s = z->w_list.a[0].x_start; e = z->w_list.a[z->w_list.n-1].x_end+1;
    s += bd; e -= bd; if(s >= e) return 0;
    memset(&ti, 0, sizeof(ti)); 
    ti.z = q; ti.qs = s; ti.qe = e; 
    ti.qi = t->w_list.a[0].x_start;
    ti.ti = t->w_list.a[0].y_start;
    ti.wi = 0; ti.wn = t->w_list.n;
    ti.tot = 0; ti.ci = 0; ti.lc = -1; ti.cql = ti.ctl = 0;
    ti.qis = ti.qie = ti.qi; 
    ti.tis = ti.tie = ti.ti; 
    ti.lerr = -1;

    while((!citer_end(qi)) && (!citer_end(ti))) {
        if(qi.qi == ti.qi) {///check if it is been masked

        } else if((qi.lc == 0) && (ti.lc == 0)) {///check if it is been masked

        }

        if(qi.qi < ti.qi) {

        } if(ti.qi < qi.qi) {

        } else {//qi.qi == ti.qi

        }
    }
    return 0;
}

void retrieve_cigar_xcoord(bit_extz_t *ez, int64_t is, int64_t ie, int64_t *yk, int64_t *xk, int64_t *ck, int64_t *rxs, int64_t *rxe)
{
    if(!ez->cigar.n) return;
    int64_t cn = ez->cigar.n, op; int64_t ws, we, ws0, we0, os, oe;
    if(((*ck) < 0) || ((*ck) > cn)) {//(*ck) == cn is allowed
        (*ck) = 0; (*yk) = ez->ps; (*xk) = ez->ts; 
    }

    while ((*ck) > 0 && (*yk) > is) {
        --(*ck);
        op = ez->cigar.a[(*ck)]>>14;
        if(op!=3) (*yk) -= (ez->cigar.a[(*ck)]&(0x3fff));
        if(op!=2) (*xk) -= (ez->cigar.a[(*ck)]&(0x3fff));
    }

    //some cigar will span s or e
    while ((*ck) < cn && (*yk) < ie) {//[s, e)
        ws = (*yk); ws0 = (*xk); 
        op = ez->cigar.a[(*ck)]>>14;
        if(op!=3) (*yk) += (ez->cigar.a[(*ck)]&(0x3fff));
        if(op!=2) (*xk) += (ez->cigar.a[(*ck)]&(0x3fff));
        we = (*yk); we0 = (*xk); 
        if(ws >= is && we <= ie) {///the cigar is fully contained
            if((*rxs) > ws0) (*rxs) = ws0;
            if((*rxe) < we0) (*rxe) = we0;
        } else if(op == 0 || op == 1) {///overlap, it is hard to handle indels
            os = MAX(is, ws); oe = MIN(ie, we);
            if(oe > os) {
                ws0 += (os - ws);
                we0 = ws0 + (oe - os);
                if((*rxs) > ws0) (*rxs) = ws0;
                if((*rxe) < we0) (*rxe) = we0;
            }
        }
        (*ck)++;
    }
} 

uint64_t extract_xcoordates0(overlap_region *z, int64_t ys, int64_t ye, int64_t *rxs, int64_t *rxe, ul_ov_t *p)
{
    int64_t wk = ovlp_cur_wid(*p), yk = ovlp_cur_xoff(*p), xk = ovlp_id(*p), ck = ovlp_cur_coff(*p);
    int64_t min_w = ovlp_min_wid(*p), max_w = ovlp_max_wid(*p);//[min_w, max_w]
    bit_extz_t ez; window_list *m; (*rxs) = INT32_MAX; (*rxe) = -1;
    if(ys >= ye) return 0;
    ///[ys, ye) but [z->w_list.a[wk].y_start, z->w_list.a[wk].y_end]
    int64_t ws, we, os, oe, ovlp, yl, tot = ye - ys;
    if(wk < min_w || wk > max_w) wk = min_w;
    for (; wk >= min_w && z->w_list.a[wk].y_start > ys; wk--); 
    if(wk < min_w || wk > max_w) return 0;
    for (; wk <= max_w && z->w_list.a[wk].y_end < ys; wk++); 
    if(wk < min_w || wk > max_w) return 0;
    //s >= w_list.a[wk].x_start && s <= w_list.a[wk].x_end
    if(wk != ovlp_cur_wid(*p)) {//xk is global, while ck is local
        yk = z->w_list.a[wk].y_start; 
        xk = z->w_list.a[wk].x_start; 
        ck = 0;
    }

    ///[ys, ye) but [z->w_list.a[wk].y_start, z->w_list.a[wk].y_end]
    while(wk <= max_w && z->w_list.a[wk].y_start < ye) {
        m = &(z->w_list.a[wk]); 
        ws = m->y_start; we = m->y_end+1;
        os = MAX(ys, ws); oe = MIN(ye, we);
        ovlp = ((oe>os)? (oe-os):0);
        
        if(ovlp) {
            yl = m->y_end+1-m->y_start;
            if((is_ualn_win((*m))) || (is_est_aln((*m)))) {
                if(ovlp >= yl) {
                    if((*rxs) > m->x_start) (*rxs) = m->x_start;
                    if((*rxe) < (m->x_end+1)) (*rxe) = m->x_end+1;
                }
                //skip the whole window
                yk = m->y_end+1; xk = m->x_end+1; ck = m->clen;
            } else {
                if(ovlp == yl) {
                    //skip the whole window
                    yk = m->y_end+1; xk = m->x_end+1; ck = m->clen;
                    if((*rxs) > m->x_start) (*rxs) = m->x_start;
                    if((*rxe) < (m->x_end+1)) (*rxe) = m->x_end+1;
                } else {
                    set_bit_extz_t(ez, (*z), wk); 
                    retrieve_cigar_xcoord(&ez, os, oe, &yk, &xk, &ck, rxs, rxe);
                }                
            }
        }
        tot -= ovlp; 
        if(yk >= ye) break;//[min_w, max_w] && [ys, ye)
        wk++; if(wk > max_w) break;
        yk = z->w_list.a[wk].y_start; xk = z->w_list.a[wk].x_start; ck = 0;//reset
    }
    assert(!tot);
    ovlp_cur_wid(*p) = wk; ovlp_cur_xoff(*p) = yk; ovlp_id(*p) = xk; ovlp_cur_coff(*p) = ck;
    if((*rxe) > (*rxs)) return 1;
    return 0;
}

uint64_t extract_xcoordates(overlap_region *z, uint64_t *a, int64_t a_n, asg64_v *b)
{
    int64_t i, zwn, q[2], t[2], w[2], bn = b->n, a_i, a_z, as, ae, is, ie, os, oe, ovlp; 
    ul_ov_t m; int64_t qs, qe;
    zwn = z->w_list.n; if(!zwn) return 0;
    q[0] = q[1] = t[0] = t[1] = w[0] = w[1] = INT32_MIN; memset(&m, 0, sizeof(m));
    for (i = a_i = 0; i < zwn; i++) {
        if((z->w_list.a[i].x_start==(q[1]+1)) && ((z->w_list.a[i].y_start==(t[1]+1)))) {
            q[1] = z->w_list.a[i].x_end;
            t[1] = z->w_list.a[i].y_end;
            w[1] = i;
        } else {
            if(q[0] != INT32_MIN) {
                ovlp_id(m) = z->w_list.a[w[0]].x_start; ///cur xpos
                ovlp_min_wid(m) = w[0]; ///beg id of windows
                ovlp_max_wid(m) = w[1]; ///end id of windows
                ovlp_cur_wid(m) = w[0]; ///cur id of windows
                ovlp_cur_xoff(m) = z->w_list.a[w[0]].y_start; ///cur ypos
                ovlp_cur_coff(m) = 0; ///cur cigar off in cur window

                is = t[0]; ie = t[1]+1; ///[is, ie)
                for (a_z = a_i; a_z >= 0; a_z--) {
                    as = a[a_z]>>32; ae = (uint32_t)a[a_z];
                    if(ae <= is) break;
                }
                if(a_z < 0) a_z = 0;
                for (; a_z < a_n; a_z++) {
                    as = a[a_z]>>32; ae = (uint32_t)a[a_z];
                    os = MAX(is, as); oe = MIN(ie, ae);
                    ovlp = ((oe>os)? (oe-os):0);
                    if(ovlp) {
                        if(extract_xcoordates0(z, os, oe, &qs, &qe, &m)) {
                            kv_push(uint64_t, *b, (((uint64_t)qs)<<32)|((uint64_t)qe));
                        }
                    }
                    if(as >= ie) break;
                }
                a_i = a_z;
            }
            q[0] = z->w_list.a[i].x_start; q[1] = z->w_list.a[i].x_end;
            t[0] = z->w_list.a[i].y_start; t[1] = z->w_list.a[i].y_end;
            w[0] = i; w[1] = i;
        }
    }

    if(q[0] != INT32_MIN) {
        ovlp_id(m) = z->w_list.a[w[0]].x_start; ///cur xpos
        ovlp_min_wid(m) = w[0]; ///beg id of windows
        ovlp_max_wid(m) = w[1]; ///end id of windows
        ovlp_cur_wid(m) = w[0]; ///cur id of windows
        ovlp_cur_xoff(m) = z->w_list.a[w[0]].y_start; ///cur ypos
        ovlp_cur_coff(m) = 0; ///cur cigar off in cur window

        is = t[0]; ie = t[1]+1; 
        for (a_z = a_i; a_z >= 0; a_z--) {
            as = a[a_z]>>32; ae = (uint32_t)a[a_z];
            if(ae <= is) break;
        }
        if(a_z < 0) a_z = 0;
        for (; a_z < a_n; a_z++) {
            as = a[a_z]>>32; ae = (uint32_t)a[a_z];
            os = MAX(is, as); oe = MIN(ie, ae);
            ovlp = ((oe>os)? (oe-os):0);
            if(ovlp) {
                if(extract_xcoordates0(z, os, oe, &qs, &qe, &m)) {
                    kv_push(uint64_t, *b, (((uint64_t)qs)<<32)|((uint64_t)qe));
                }
            }
            if(as >= ie) break;
        }
        a_i = a_z;
    }
    return b->n - bn;
}

int64_t cal_xerr(overlap_region *z, uint64_t *a, int64_t a_n)
{
    int64_t i, zwn, q[2], t[2], w[2], a_i, a_z, as, ae, is, ie, os, oe, ovlp, err = 0; ul_ov_t m; 
    zwn = z->w_list.n; if(!zwn) return 0;
    q[0] = q[1] = t[0] = t[1] = w[0] = w[1] = INT32_MIN; memset(&m, 0, sizeof(m));
    for (i = a_i = 0; i < zwn; i++) {
        // fprintf(stderr, "[M::%s]\tw::[%d,\t%d)\n", __func__, z->w_list.a[i].x_start, z->w_list.a[i].x_end+1);
        if((z->w_list.a[i].x_start==(q[1]+1)) && ((z->w_list.a[i].y_start==(t[1]+1)))) {
            q[1] = z->w_list.a[i].x_end;
            t[1] = z->w_list.a[i].y_end;
            w[1] = i;
        } else {
            if(q[0] != INT32_MIN) {
                ovlp_id(m) = 0; ///ovlp id
                ovlp_min_wid(m) = w[0]; ///beg id of windows
                ovlp_max_wid(m) = w[1]; ///end id of windows
                ovlp_cur_wid(m) = w[0]; ///cur id of windows
                ovlp_cur_xoff(m) = z->w_list.a[w[0]].x_start; ///cur xpos
                ovlp_cur_coff(m) = 0; ///cur cigar off in cur window
                ovlp_bd(m) = 0;

                is = q[0]; ie = q[1]+1; ///[is, ie)
                // fprintf(stderr, "+[M::%s]\tis::%ld\tie::%ld\n", __func__, is, ie);
                for (a_z = a_i; a_z >= 0; a_z--) {
                    as = a[a_z]>>32; ae = (uint32_t)a[a_z];
                    if(ae <= is) break;
                }
                if(a_z < 0) a_z = 0;
                for (; a_z < a_n; a_z++) {
                    as = a[a_z]>>32; ae = (uint32_t)a[a_z];
                    os = MAX(is, as); oe = MIN(ie, ae);
                    ovlp = ((oe>os)? (oe-os):0);
                    if(ovlp) {
                        // fprintf(stderr, "+[M::%s]\tos::%ld\toe::%ld\n", __func__, os, oe);
                        err = err + extract_sub_cigar_err(z, os, oe, &m);
                    }
                    if(as >= ie) break;
                }
                a_i = a_z;
            }
            q[0] = z->w_list.a[i].x_start; q[1] = z->w_list.a[i].x_end;
            t[0] = z->w_list.a[i].y_start; t[1] = z->w_list.a[i].y_end;
            w[0] = i; w[1] = i;
        }
    }

    if(q[0] != INT32_MIN) {
        ovlp_id(m) = 0; ///ovlp id
        ovlp_min_wid(m) = w[0]; ///beg id of windows
        ovlp_max_wid(m) = w[1]; ///end id of windows
        ovlp_cur_wid(m) = w[0]; ///cur id of windows
        ovlp_cur_xoff(m) = z->w_list.a[w[0]].x_start; ///cur xpos
        ovlp_cur_coff(m) = 0; ///cur cigar off in cur window
        ovlp_bd(m) = 0;

        is = q[0]; ie = q[1]+1; 
        // fprintf(stderr, "-[M::%s]\tis::%ld\tie::%ld\n", __func__, is, ie);
        for (a_z = a_i; a_z >= 0; a_z--) {
            as = a[a_z]>>32; ae = (uint32_t)a[a_z];
            if(ae <= is) break;
        }
        if(a_z < 0) a_z = 0;
        for (; a_z < a_n; a_z++) {
            as = a[a_z]>>32; ae = (uint32_t)a[a_z];
            os = MAX(is, as); oe = MIN(ie, ae);
            ovlp = ((oe>os)? (oe-os):0);
            if(ovlp) {
                // if(os == 22758 && oe == 22764) {
                //     fprintf(stderr, "-[M::%s]\tos::%ld\toe::%ld\n", __func__, os, oe);
                // } 
                err = err + extract_sub_cigar_err(z, os, oe, &m);
            }
            if(as >= ie) break;
        }
        a_i = a_z;
    }
    return err;
}

void update_masks(asg64_v *in, overlap_region *qi, overlap_region *ti, int64_t bd)
{
    uint64_t in_n0 = in->n, in_n1, k, s, e; 
    int64_t zwn, i, q[2], t[2]; overlap_region *z;
    for (k = 0; k < in_n0; k++) {
        s = in->a[k]>>32; e = (uint32_t)in->a[k];
        // fprintf(stderr, "+[M::%s]\tk::%lu\tstr::[%lu,\t%lu)\n", __func__, k, s, e); 
        kv_push(uint64_t, *in, s<<1);
        kv_push(uint64_t, *in, (((e-1)<<1)|1));
    }
    in_n1 = in->n;

    uint64_t m, start, end; int64_t dp, old_dp, dp_mask, old_dp_mask;
    radix_sort_bc64(in->a+in_n0, in->a+in_n1);
    for (k = in_n0, dp = start = 0; k < in_n1; k++) {
        old_dp = dp;
        if (in->a[k]&1) --dp;
        else ++dp;

        if (old_dp < 2 && dp >= 2) {///old_dp < dp, b.a[j] is qs
            start = in->a[k]>>1;
        } else if (old_dp >= 2 && dp < 2){
            end = in->a[k]>>1;
            if(end >= start) {
                kv_push(uint64_t, *in, ((start<<32)|end));///[start, end]
                // fprintf(stderr, "-[M::%s]\tmsk::[%lu,\t%lu)\n", __func__, start, end+1); 
            }
        }
    }

    for (m = 0, k = in_n1; k < in->n; k++) in->a[m++] = in->a[k];
    in->n = in_n0 = m; ///in->a[0, in_n0) includes the masked regions

    z = qi;
    zwn = z->w_list.n;
    if(zwn > 0) {
        q[0] = q[1] = t[0] = t[1] = INT32_MIN;
        for (i = 0; i < zwn; i++) {
            if((z->w_list.a[i].x_start==(q[1]+1)) && ((z->w_list.a[i].y_start==(t[1]+1)))) {
                q[1] = z->w_list.a[i].x_end;
                t[1] = z->w_list.a[i].y_end;
            } else {
                if(q[0] != INT32_MIN) {
                    q[0] += bd; q[1] -= bd;
                    if(q[1] >= q[0]) {///unmasked regions
                        m = ((uint64_t)q[0])<<1; m <<= 32; 
                        kv_push(uint64_t, *in, m);
                        m = (((uint64_t)(q[1]))<<1)+1; m <<= 32; 
                        kv_push(uint64_t, *in, m);
                        // fprintf(stderr, "[M::%s]\tqstr::[%ld,\t%ld)\n", __func__, q[0], q[1]+1); 
                    }
                }
                q[0] = z->w_list.a[i].x_start; q[1] = z->w_list.a[i].x_end;
                t[0] = z->w_list.a[i].y_start; t[1] = z->w_list.a[i].y_end;
            }
        }
        if(q[0] != INT32_MIN) {
            q[0] += bd; q[1] -= bd;
            if(q[1] >= q[0]) {
                m = ((uint64_t)q[0])<<1; m <<= 32; 
                kv_push(uint64_t, *in, m);
                m = (((uint64_t)(q[1]))<<1)+1; m <<= 32; 
                kv_push(uint64_t, *in, m);
                // fprintf(stderr, "[M::%s]\tqstr::[%ld,\t%ld)\n", __func__, q[0], q[1]+1); 
            }
        }
    }

    z = ti;
    zwn = z->w_list.n;
    if(zwn > 0) {
        q[0] = q[1] = t[0] = t[1] = INT32_MIN;
        for (i = 0; i < zwn; i++) {
            if((z->w_list.a[i].x_start==(q[1]+1)) && ((z->w_list.a[i].y_start==(t[1]+1)))) {
                q[1] = z->w_list.a[i].x_end;
                t[1] = z->w_list.a[i].y_end;
            } else {
                if(q[0] != INT32_MIN) {
                    q[0] += bd; q[1] -= bd;
                    if(q[1] >= q[0]) {///unmasked regions
                        m = ((uint64_t)q[0])<<1; m <<= 32; 
                        kv_push(uint64_t, *in, m);
                        m = (((uint64_t)(q[1]))<<1)+1; m <<= 32; 
                        kv_push(uint64_t, *in, m);
                        // fprintf(stderr, "[M::%s]\ttstr::[%ld,\t%ld)\n", __func__, q[0], q[1]+1); 
                    }
                }
                q[0] = z->w_list.a[i].x_start; q[1] = z->w_list.a[i].x_end;
                t[0] = z->w_list.a[i].y_start; t[1] = z->w_list.a[i].y_end;
            }
        }
        if(q[0] != INT32_MIN) {
            q[0] += bd; q[1] -= bd;
            if(q[1] >= q[0]) {
                m = ((uint64_t)q[0])<<1; m <<= 32; 
                kv_push(uint64_t, *in, m);
                m = (((uint64_t)(q[1]))<<1)+1; m <<= 32; 
                kv_push(uint64_t, *in, m);
                // fprintf(stderr, "[M::%s]\ttstr::[%ld,\t%ld)\n", __func__, q[0], q[1]+1); 
            }
        }
    }


    for (k = 0; k < in_n0; k++) {
        start = in->a[k]>>32; end = (uint32_t)in->a[k];
        if(end >= start) {
            q[0] = start; q[1] = end;
            m = ((uint64_t)q[0])<<1; m <<= 32; m |= 1;
            kv_push(uint64_t, *in, m);
            m = (((uint64_t)(q[1]))<<1)+1; m <<= 32; m |= 1;
            kv_push(uint64_t, *in, m);
        }
    }
    in_n1 = in->n;

    radix_sort_bc64(in->a+in_n0, in->a+in_n1);///(uint32_t)in->a[]: 0-> original; 1-> mask
    for (k = in_n0, dp = old_dp = dp_mask = old_dp_mask = 0, start = 0, end = -1; k < in_n1; ++k) {///[beg, end) but coordinates in idx is [, ]
        ///if idx->a.a[] is qe
        old_dp = dp; old_dp_mask = dp_mask;
        if ((in->a[k]>>32)&1) {
            --dp; end = (in->a[k]>>33)+1; if((uint32_t)in->a[k]) dp_mask--;
        }else {
            //meet a new overlap; the overlaps are pushed by the x_pos_s
            ++dp; end = (in->a[k]>>33); if((uint32_t)in->a[k]) dp_mask++;
        }
        if((end > start) && (old_dp >= 2) && old_dp_mask <= 0) {
            kv_push(uint64_t, *in, ((start<<32)|end));///[start, end)
            if(qi->y_id == 4 && ti->y_id == 6) {
                fprintf(stderr, "[M::%s]\tout::[%ld,\t%ld)\n", __func__, start, end); 
            }
        }
        start = end; 
    }

    for (m = 0, k = in_n1; k < in->n; k++) in->a[m++] = in->a[k];
    in->n = in_n0 = m; 
}

int64_t cal_paired_distance(ma_ug_t *ug, overlap_region *q, overlap_region *t, ul_ov_t *a, uint32_t a_n, 
asg64_v *srt, asg64_v* buf1, idx_emask_t *mm, int64_t bd)
{
    uint64_t k, cn, *qa, *ta, *ca, m, qn, tn, rev, rev_n, tl, mt, qocc, tocc, s, e; 
    ul_ov_t ou; memset((&ou), 0, sizeof(ou)); int64_t eq, et;
    srt->n = 0; kv_resize(uint64_t, *srt, (a_n<<1)); qa = srt->a; ta = qa + a_n;

    // if(q->y_id == 4 && t->y_id == 6) {
    //     fprintf(stderr, "\n\n\n[M::%s]\tutg%.6u%c\tutg%.6u%c\n", __func__, q->y_id+1, "lc"[ug->u.a[q->y_id].circ], 
    //     t->y_id+1, "lc"[ug->u.a[t->y_id].circ]);
    // }

    for (k = cn = 0; k < a_n; k++) {
        // if(q->y_id == 4 && t->y_id == 6) {
        //     fprintf(stderr, "[M::%s]\tk::%lu\n", __func__, k); 
        //     fprintf(stderr, "au\tutg%.6u%c\txl::%u\tx::[%u,\t%u)\t%c\tutg%.6u%c\tyl::%u\ty::[%u,\t%u)\tfull::%u\tis_cal::%u\n", 
        //     a[k].qn + 1, "lc"[ug->u.a[a[k].qn].circ], ug->u.a[a[k].qn].len, a[k].qs, a[k].qe, "+-"[a[k].rev],
        //     a[k].tn + 1, "lc"[ug->u.a[a[k].tn].circ], ug->u.a[a[k].tn].len, a[k].ts, a[k].te,
        //     a[k].el, (a[k].sec!=((uint32_t)(0x3fffffff)))?1:0);
        // }
        if(!cal_no_bd_coor(&(a[k]), &ou, mm, bd)) continue;
        qa[cn] = (((uint64_t)ou.qs)<<32)|((uint64_t)ou.qe);
        ta[cn] = (((uint64_t)ou.ts)<<32)|((uint64_t)ou.te);
        cn++;
        // if(q->y_id == 4 && t->y_id == 6) {
        //     fprintf(stderr, "ou\tutg%.6u%c\txl::%u\tx::[%u,\t%u)\t%c\tutg%.6u%c\tyl::%u\ty::[%u,\t%u)\n", 
        //     ou.qn + 1, "lc"[ug->u.a[ou.qn].circ], ug->u.a[ou.qn].len, ou.qs, ou.qe, "+-"[ou.rev],
        //     ou.tn + 1, "lc"[ug->u.a[ou.tn].circ], ug->u.a[ou.tn].len, ou.ts, ou.te);
        // }
    }
    radix_sort_bc64(qa, qa+cn); ca = qa; rev = q->y_pos_strand; tl = ug->g->seq[q->y_id].len;
    for (k = m = 0; k < cn; k++) {
        if(m <= 0 || (((uint32_t)ca[m-1])) < (ca[k]>>32)) {
            ca[m++] = ca[k];
        } else if(((uint32_t)ca[m-1]) < ((uint32_t)ca[k])) {
            ca[m-1] += ((uint32_t)ca[k])-((uint32_t)ca[m-1]);
        }
    }
    if(rev) {
        rev_n = m>>1; ///tl = (tl<<32) + tl;
        for (k = 0; k < rev_n; k++) {
            mt = ca[k]; ca[k] = ca[m-k-1]; ca[m-k-1] = mt;

            s = tl-((uint32_t)ca[k]); e = tl-(ca[k]>>32);
            ca[k] = (s<<32)|e;
            s = tl-((uint32_t)ca[m-k-1]); e = tl-(ca[m-k-1]>>32);
            ca[m-k-1] = (s<<32)|e;
        }
        if(m&1) {
            s = tl-((uint32_t)ca[k]); e = tl-(ca[k]>>32);
            ca[k] = (s<<32)|e;
        }
    }
    qn = m;
    // if(q->y_id == 4 && t->y_id == 6) {
    // fprintf(stderr, "[M::%s]\tqn::%lu\trev::%lu\n", __func__, qn, rev); 
    // for (k = 0; k < qn; k++) {
    //     s = qa[k]>>32; e = ((uint32_t)qa[k]);
    //     fprintf(stderr, "[M::%s]\tqstr::[%lu,\t%lu)\n", __func__, s, e); 
    // }
    // }

    radix_sort_bc64(ta, ta+cn); ca = ta; rev = t->y_pos_strand; tl = ug->g->seq[t->y_id].len;
    for (k = m = 0; k < cn; k++) {
        if(m <= 0 || (((uint32_t)ca[m-1])) < (ca[k]>>32)) {
            ca[m++] = ca[k];
        } else if(((uint32_t)ca[m-1]) < ((uint32_t)ca[k])) {
            ca[m-1] += ((uint32_t)ca[k])-((uint32_t)ca[m-1]);
        }
    }
    if(rev) {
        rev_n = m>>1; ///tl = (tl<<32) + tl;
        for (k = 0; k < rev_n; k++) {
            mt = ca[k]; ca[k] = ca[m-k-1]; ca[m-k-1] = mt;

            s = tl-((uint32_t)ca[k]); e = tl-(ca[k]>>32);
            ca[k] = (s<<32)|e;
            s = tl-((uint32_t)ca[m-k-1]); e = tl-(ca[m-k-1]>>32);
            ca[m-k-1] = (s<<32)|e;
        }
        if(m&1) {
            s = tl-((uint32_t)ca[k]); e = tl-(ca[k]>>32);
            ca[k] = (s<<32)|e;
        }
    }
    tn = m;
    // if(q->y_id == 4 && t->y_id == 6) {
    // fprintf(stderr, "[M::%s]\ttn::%lu\trev::%lu\n", __func__, tn, rev); 
    // for (k = 0; k < tn; k++) {
    //     s = ta[k]>>32; e = ((uint32_t)ta[k]);
    //     fprintf(stderr, "[M::%s]\ttstr::[%lu,\t%lu)\n", __func__, s, e); 
    // }
    // }

    buf1->n = 0; qocc = tocc = 0; srt->n = 0;
    if(qn > 0 && tn > 0) {
        qocc = extract_xcoordates(q, qa, qn, buf1);
        tocc = extract_xcoordates(t, ta, tn, buf1);
        // if(q->y_id == 4 && t->y_id == 6) {
        //     fprintf(stderr, "[M::%s]\tqocc::%lu\ttocc::%lu\n", __func__, qocc, tocc);
        // }
        if((!qocc) || (!tocc)) qocc = tocc = buf1->n = 0;
    }
    if(!(buf1->n)) return INT32_MIN;
    update_masks(buf1, q, t, bd);   

    eq = et = 0;
    if(buf1->n) {
        eq = cal_xerr(q, buf1->a, buf1->n);
        et = cal_xerr(t, buf1->a, buf1->n);
    }   
    // if(q->y_id == 4 && t->y_id == 6) {
    //     s = MAX(q->x_pos_s, t->x_pos_s);
    //     e = MIN(q->x_pos_e, t->x_pos_e) + 1;
    //     buf1->n = 0;
    //     kv_push(uint64_t, *buf1, ((s<<32)|e));
    //     fprintf(stderr, "[M::%s]\teq::%ld\tet::%ld\tbuf1->n::%u\n", __func__, eq, et, (uint32_t)buf1->n);
    //     eq = cal_xerr(q, buf1->a, buf1->n); et = cal_xerr(t, buf1->a, buf1->n);
    //     fprintf(stderr, "[M::%s]\teq1::%ld\tet1::%ld\tbuf1->n::%u\ts::%lu\te::%lu\n", 
    //     __func__, eq, et, (uint32_t)buf1->n, s, e);
    // }

    return eq - et;    
    // return get_pe_diff(q, qa, qn, t, ta, tn, bd);
}

void prt_masks(ma_ug_t *ug, ul_ov_t *a, uint64_t a_n)
{
    uint64_t k;
    for (k = 0; k < a_n; k++) {
        fprintf(stderr, "[M::%s]\tutg%.6u%c\txl::%u\tx::[%u,\t%u)\t%c\tutg%.6u%c\tyl::%u\ty::[%u,\t%u)\tfull::%u\tis_cal::%u\n", 
        __func__, 
        a[k].qn + 1, "lc"[ug->u.a[a[k].qn].circ], ug->u.a[a[k].qn].len, a[k].qs, a[k].qe, "+-"[a[k].rev],
        a[k].tn + 1, "lc"[ug->u.a[a[k].tn].circ], ug->u.a[a[k].tn].len, a[k].ts, a[k].te,
        a[k].el, (a[k].sec!=((uint32_t)(0x3fffffff)))?1:0);
    }
}

uint32_t gen_mask_ovlp0(ma_ug_t *ug, overlap_region *a, uint32_t qi, uint32_t ti, kv_ul_ov_t *ov_db, kv_ul_ov_t *buf, idx_emask_t *mm, double len_diff, int64_t rlen, uint64_t ovdb_idx, asg64_v *coor_srt, 
asg64_v* buf1, int64_t bd, ul_ov_t *res)
{
    ul_ov_t *ref, r0, r1; uint32_t bn = buf->n, s, e; uint64_t is_exact = 0; 
    overlap_region *q = &(a[qi]), *t = &(a[ti]); int64_t dd = 0;
    assert(q->y_id < t->y_id); 
    assert(ovdb_idx < ov_db->n);///quickly jump to the related qn-tn pair ov_db[]
    cal_x_ul_ovlp(ug, q, t, &r0); 
    r1 = r0; 
    r1.qn = r0.tn; r1.tn = r0.qn;
    r1.qs = r0.ts; r1.qe = r0.te;
    r1.ts = r0.qs; r1.te = r0.qe;

    if((ov_db->a[ovdb_idx].qn==q->y_id) && (ov_db->a[ovdb_idx].tn==t->y_id)) {
        ref = &(r0); s = ovdb_idx; e = q->overlapLen+q->x_pos_strand;
        if((s!=((uint32_t)-1))&&(e > s)) push_consist_ul_ovlps(ug, ov_db, s, e, len_diff, ref, mm, buf, 0, &is_exact);

        if(!is_exact) {
            ref = &(r1); s = t->overlapLen; e = t->overlapLen+t->x_pos_strand;
            if((s!=((uint32_t)-1))&&(e > s)) push_consist_ul_ovlps(ug, ov_db, s, e, len_diff, ref, mm, buf, 1, &is_exact);
        }
    } else {
        assert((ov_db->a[ovdb_idx].tn==q->y_id) && (ov_db->a[ovdb_idx].qn==t->y_id));
        ref = &(r1); s = ovdb_idx; e = t->overlapLen+t->x_pos_strand;
        if((s!=((uint32_t)-1))&&(e > s)) push_consist_ul_ovlps(ug, ov_db, s, e, len_diff, ref, mm, buf, 1, &is_exact);

        if(!is_exact) {
            ref = &(r0); s = q->overlapLen; e = q->overlapLen+q->x_pos_strand;
            if((s!=((uint32_t)-1))&&(e > s)) push_consist_ul_ovlps(ug, ov_db, s, e, len_diff, ref, mm, buf, 0, &is_exact);
        }
    }
    // prt_masks(ug, buf->a + bn, buf->n - bn);
    ///p->qn/p->tn:: id within the ol->list
    ///p->el:: if equally best
    ///p->qs:: err(query)-err(target)
    ///p->ts:: err(target)-err(query)
    ref = &(r0); 
    res->qn = qi; res->tn = ti;
    res->el = 1; res->qs = res->ts = 0;
    if(!is_exact) {
        assert(buf->n > bn);
        dd = cal_paired_distance(ug, q, t, buf->a + bn, buf->n - bn, coor_srt, buf1, mm, bd);
        if(dd != 0) {
            res->el = 0;
            if(dd > 0) {///q has more error than t
                res->qs = dd; res->ts = 0;
            } else {
                res->qs = 0; res->ts = -dd;
            }
        }
    }
    buf->n = bn;
    if(dd!=INT32_MIN) return 1;
    return 0;
}

uint64_t if_direct_mask(const ul_idx_t *uref, const ug_opt_t *uopt, ul_ov_t *p, ul_ov_t *q, uint64_t bw, double len_diff)
{
    ul_ov_t *li, *lj; uint64_t os, oe, ovlp;
    if(p->qe > q->qe) li = p, lj = q;
    else if(p->qe == q->qe && p->qs >= q->qs) li = p, lj = q;
    else lj = p, li = q;
    os = MAX(li->qs, lj->qs), oe = MIN(li->qe, lj->qe);
    ovlp = ((oe > os)? (oe - os):0);
    if(!ovlp) return 0;//no overlap

    if(lj->qs <= li->qs+G_CHAIN_INDEL) {
        if(govlp_check(uref, uopt, bw, len_diff, li, lj)) return 1;
    } else if((lj->qe+G_CHAIN_INDEL>=li->qe) && (lj->qs+G_CHAIN_INDEL>=li->qs)) {
        if(govlp_check(uref, uopt, bw, len_diff, lj, li)) return 1;
    }
    return 0;
}

uint64_t mask_ovlps(const ul_idx_t *uref, const ug_opt_t *uopt, ma_ug_t *ug, overlap_region_alloc* ol, mask_ul_ov_t *mk, kv_ul_ov_t *ov_db, asg64_v *coor_srt, asg64_v* buf1, idx_emask_t *mm, double len_diff, uint64_t rlen, int64_t bd)
{
    mk->srt.n = mk->idx.n = 0;
    int64_t k, m, i, osrt_n, osrt_n1, on = ol->length; uint64_t wt, w[2], is_exact, sid, tot, tol, qi, ti; 
    overlap_region *z; ul_ov_t ou; kv_ul_ov_t *osrt = &(mk->srt); ul_ov_t p, q;
    for (k = osrt->n = 0; k < on; k++) {
        z = &(ol->list[k]); if(!(z->x_pos_strand)) continue;///x_pos_strand: how many masks
        if(gen_aln_ul_ov_t(k, ug->g->seq[z->y_id].len, z, &ou)) {
            kv_push(ul_ov_t, *osrt, ou);
        }
    }
    if(!(osrt->n)) return 0;

    osrt_n = osrt->n; memset(&ou, 0, sizeof(ou));
    radix_sort_uov_srt_qs(osrt->a, osrt->a+osrt->n); ///sort by qs; for quick filtering between overlaps
    for (k = tot = tol = 0; k < osrt_n; k++) {
        tol += osrt->a[k].qe-osrt->a[k].qs; 
        convert_ul_ov_t(&p, &(ol->list[osrt->a[k].qn]), ug); p.qn = osrt->a[k].qn;
        for (i = k+1; i < osrt_n && osrt->a[i].qs < osrt->a[k].qe; i++) {
            ///for a pair of overlap, only need to calculate it in one side
            if(ol->list[osrt->a[i].qn].y_id >= ol->list[osrt->a[k].qn].y_id) continue;
            is_exact = 0; sid = (uint32_t)-1;
            convert_ul_ov_t(&q, &(ol->list[osrt->a[i].qn]), ug); q.qn = osrt->a[i].qn;
            if(if_direct_mask(uref, uopt, &p, &q, 16, len_diff)) is_exact = 1;
            if((!is_exact) && (!check_mask_exist(ug, &(ol->list[osrt->a[i].qn]), &(ol->list[osrt->a[k].qn]), ov_db, mm, len_diff, &is_exact, &sid))) {
                continue;
            }
            
            if(is_exact) {///if two overlaps are exactly the same; no need to do anything
                ou.qs = 0; ou.qe = sid; ou.ts = ou.te = 0; ou.el = 1;
            } else {
                w[0] = ol->list[osrt->a[k].qn].non_homopolymer_errors;
                w[1] = ol->list[osrt->a[i].qn].non_homopolymer_errors;
                wt = w[0] + w[1] + ((w[0]>=w[1])?(w[0]-w[1]):(w[1]-w[0]));
                if(wt > (uint32_t)-1) wt = (uint32_t)-1;
                //qs: weight; qe: idx with in ov_db; el: exact match?
                ou.qs = wt; ou.qe = sid; ou.el = 0;
                ou.ts = osrt->a[i].qe-osrt->a[i].qs;
                ou.te = osrt->a[k].qe-osrt->a[k].qs;
            }
            ou.qn = osrt->a[i].qn; ou.tn = osrt->a[k].qn;///id of the ol->list, instead of real yid
            tot += ou.ts + ou.te;
            kv_push(ul_ov_t, *osrt, ou);
        }
        for (i = k-1; i >= 0; i--) {
            ///for a pair of overlap, only need to calculate it in one side
            if(ol->list[osrt->a[i].qn].y_id >= ol->list[osrt->a[k].qn].y_id) continue;
            is_exact = 0; sid = (uint32_t)-1;
            convert_ul_ov_t(&q, &(ol->list[osrt->a[i].qn]), ug); q.qn = osrt->a[i].qn;
            if(if_direct_mask(uref, uopt, &p, &q, 16, len_diff)) is_exact = 1;
            if((!is_exact) && (!check_mask_exist(ug, &(ol->list[osrt->a[i].qn]), &(ol->list[osrt->a[k].qn]), ov_db, mm, len_diff, &is_exact, &sid))) {
                continue;
            }
            
            if(is_exact) {
                ou.qs = 0; ou.qe = sid; ou.ts = ou.te = 0; ou.el = 1;
            } else {
                w[0] = ol->list[osrt->a[k].qn].non_homopolymer_errors;
                w[1] = ol->list[osrt->a[i].qn].non_homopolymer_errors;
                wt = w[0] + w[1] + ((w[0]>=w[1])?(w[0]-w[1]):(w[1]-w[0]));
                if(wt > (uint32_t)-1) wt = (uint32_t)-1;
                ou.qs = wt; ou.qe = sid; ou.el = 0;
                ou.ts = osrt->a[i].qe-osrt->a[i].qs;
                ou.te = osrt->a[k].qe-osrt->a[k].qs;
            }
            ou.qn = osrt->a[i].qn; ou.tn = osrt->a[k].qn;///id of the ol->list, instead of real yid
            tot += ou.ts + ou.te;
            kv_push(ul_ov_t, *osrt, ou);
        }
    }

    // fprintf(stderr, "0[M::%s]\tosrt_n::%u\tosrt->n::%u\n", __func__, (uint32_t)osrt_n, (uint32_t)osrt->n);

    if((tot > (rlen*256)) && (tot > (tol*32))) {
        tot = MAX((rlen*256), (tol*32)); osrt_n1 = osrt->n;
        radix_sort_uov_srt_qs(osrt->a+osrt_n, osrt->a+osrt->n); ///sort by weight (qs)
        for (k = osrt_n, wt = 0; k < osrt_n1 && wt <= tot; k++) {
            if(osrt->a[k].el) continue;///exact match, no additional work
            wt += (uint64_t)osrt->a[k].ts + (uint64_t)osrt->a[k].te;
        }
        osrt->n = k;
    }

    memset(&ou, 0, sizeof(ou)); osrt_n1 = osrt->n;
    ///note: osrt will be updated here, so do not use the address within osrt
    for (k = osrt_n, m = 0; k < osrt_n1; k++) {
        qi = osrt->a[k].qn; ti = osrt->a[k].tn; ///represent a pair of read ol->list[qi] <-> ol->list[ti]
        sid = osrt->a[k].qe; wt = 1;
        ///p->qn/p->tn:: id within the ol->list
        ///p->el:: if equally best
        ///p->qs:: err(query)-err(target)
        ///p->ts:: err(target)-err(query)
        if(osrt->a[k].el) {///exact match; no need base-check
            ou.qn = qi; ou.tn = ti; ou.el = 1; ou.qs = ou.ts = 0;
        } else {///do base-check
            wt = gen_mask_ovlp0(ug, ol->list, qi, ti, ov_db, osrt, mm, len_diff, rlen, sid, coor_srt, buf1, bd, &(ou));
        }
        if(wt) {
            osrt->a[m++] = ou;
            // fprintf(stderr, "[M::%s]\tutg%.6u%c\tqerr::%u\t\tutg%.6u%c\tterr::%u\tel::%u\n", 
            // __func__, ol->list[ou.qn].y_id+1, "lc"[ug->u.a[ol->list[ou.qn].y_id].circ], ou.qs, 
            // ol->list[ou.tn].y_id+1, "lc"[ug->u.a[ol->list[ou.tn].y_id].circ], ou.ts, ou.el);
        }
    }
    osrt->n = m;
    if(!(osrt->n)) return 0;

    ///double
    ///p->qn/p->tn:: id within the ol->list
    ///p->el:: if equally best
    ///p->qs:: err(query)-err(target)
    ///p->ts:: err(target)-err(query)
    kv_resize(ul_ov_t, *osrt, (osrt->n<<1)); 
    memcpy(osrt->a+osrt->n, osrt->a, (sizeof((*(osrt->a)))*osrt->n));
    osrt_n = osrt->n; osrt->n <<= 1; osrt_n1 = osrt->n;
    for (k = osrt_n; k < osrt_n1; k++) {
        m = osrt->a[k].qn; osrt->a[k].qn = osrt->a[k].tn; osrt->a[k].tn = m;
        m = osrt->a[k].qs; osrt->a[k].qs = osrt->a[k].ts; osrt->a[k].ts = m;
    }

    ///osrt = &(mk->srt);
    radix_sort_ul_ov_srt_qn1(osrt->a, osrt->a+osrt->n);
    kv_resize(uint64_t, mk->idx, ol->length); mk->idx.n = ol->length;
    memset(mk->idx.a, 0, sizeof((*(mk->idx.a)))*ol->length);
    for (k = 1, i = 0; k <= osrt_n1; k++) {
        if(k == osrt_n1 || osrt->a[i].qn != osrt->a[k].qn) {
            if(k - i > 1) radix_sort_ul_ov_srt_tn1(osrt->a+i, osrt->a+k);
            mk->idx.a[osrt->a[i].qn] = (((uint64_t)i)<<32)|((uint64_t)k);
            i = k;
        }
    }

    return 1;
}

void refine_rphase_back(overlap_region *za, uint64_t zid, ul_ov_t *a, uint64_t a_n, asg64_v *buf)
{
    uint64_t k, bn, qs, qe, err, zs, ze, zerr, m, zwn; overlap_region *z; window_list *p; 
    kv_resize(uint64_t, *buf, a_n);
    for (k = buf->n = 0; k < a_n; k++) {
        if(a[k].qs > 0) {//za[zid] has higher error rate
            // buf->a[buf->n++] = (((uint64_t)za[a[k].tn].x_pos_s)<<32)|((uint64_t)a[k].tn);
            buf->a[buf->n++] = (((uint64_t)za[a[k].tn].x_pos_s)<<32)|(k);
        } 
    }
    radix_sort_bc64(buf->a, buf->a+buf->n);

    bn = buf->n;
    for (k = 0; k < bn; k++) {
        qs = qe = err = 0;
        if(buf->n > bn) {
            qs = buf->a[buf->n-2]>>32; 
            qe = ((uint32_t)buf->a[buf->n-2]);
            err = buf->a[buf->n-1];
        }
        zs = MIN(za[a[((uint64_t)buf->a[k])].tn].x_pos_s, za[zid].x_pos_s);
        ze = MAX(za[a[((uint64_t)buf->a[k])].tn].x_pos_e, za[zid].x_pos_e) + 1;
        if(ze <= zs) continue; 
        zerr = a[((uint64_t)buf->a[k])].qs;
        if(qe <= zs) {
            kv_push(uint64_t, *buf, ((zs<<32)|ze));
            kv_push(uint64_t, *buf, zerr);
        } else {
            if(ze > qe) qe = ze;
            if(zerr > err) err = zerr;
            buf->a[buf->n-2] = ((qs<<32)|(qe));
            buf->a[buf->n-1] = err;
        }
    }

    for (k = bn, m = 0; k < buf->n; k++) buf->a[m++] = buf->a[k];
    buf->n = m;
    if(!(buf->n)) return;
    
    z = &(za[zid]);
    for (k = 0; k < z->align_length; k++) {
        p = &(z->w_list.a[z->w_list.n+k]);
        if(p->clen <= 0) continue;
        zs = p->x_start; ze = p->x_end; ///note here is [p->x_start, p->x_end)
        zerr = p->clen;
        
        kv_push(uint64_t, *buf, ((zs<<32)|ze));
        kv_push(uint64_t, *buf, zerr);
    }

    bn = buf->n; kv_resize(uint64_t, *buf, (bn + (bn>>1)));
    for (k = 0; k < bn; k+=2) {
        zs = buf->a[k]>>32; zs <<= 32; zs |= k;
        kv_push(uint64_t, *buf, zs);
    }

    uint64_t *wa = buf->a + bn, wan = buf->n - bn;
    radix_sort_bc64(wa, wa + wan); 
    z->x_pos_strand = (uint32_t)-1;
    z->w_list.n = 0;
    
    for (k = 0; k < wan; k++) {
        qs = qe = z->x_pos_s; err = 0;
        if(z->w_list.n) {
            qs = z->w_list.a[z->w_list.n-1].x_start;
            qe = z->w_list.a[z->w_list.n-1].x_end;
            err = z->w_list.a[z->w_list.n-1].clen;
        }

        zs = buf->a[((uint64_t)wa[k])]>>32;
        ze = (uint32_t)(buf->a[((uint64_t)wa[k])]);        
        zerr = buf->a[((uint64_t)wa[k])+1];

         if(qe <= zs) {
            if(qe < zs) {
                kv_pushp(window_list, z->w_list, &p); memset(p, 0, sizeof((*p)));
                p->clen = 0; p->x_start = qe; p->x_end = zs;
            }
            kv_pushp(window_list, z->w_list, &p); memset(p, 0, sizeof((*p)));
            p->clen = zerr; p->x_start = zs; p->x_end = ze;
        } else {
            if(ze > qe) qe = ze;
            if(zerr > err) err = zerr;
            p = &(z->w_list.a[z->w_list.n-1]); 
            p->clen = err;
            p->x_start = qs; 
            p->x_end = qe;
        }
    }
    

    qs = qe = z->x_pos_s; zs = z->x_pos_e + 1;
    if(z->w_list.n) {
        qs = z->w_list.a[z->w_list.n-1].x_start;
        qe = z->w_list.a[z->w_list.n-1].x_end;
    }
    if(qe < zs) {
        kv_pushp(window_list, z->w_list, &p); memset(p, 0, sizeof((*p)));
        p->clen = 0; p->x_start = qe; p->x_end = zs;
    }

    z->x_pos_strand = (uint32_t)-1;
    z->overlapLen = z->x_pos_e+1-z->x_pos_s; 
    z->non_homopolymer_errors = 0; zwn = 0;
    for (k = 0; k < z->w_list.n; k++) {
        if(z->w_list.a[k].clen > 0) {
            z->non_homopolymer_errors += z->w_list.a[k].clen;
            zwn += z->w_list.a[k].x_end-z->w_list.a[k].x_start;
        } 
    }
    assert(zwn <= z->overlapLen);///zwn is the length with secondary-best alignment
    z->align_length = z->overlapLen - zwn;
}

inline void push_rphase(overlap_region *z, uint64_t zs, uint64_t ze, uint64_t zerr, uint64_t is_pri)
{
    uint64_t pe, p_pri; window_list *p; 
    pe = z->x_pos_s; p_pri = 0;
    if(z->w_list.n) {
        pe = z->w_list.a[z->w_list.n-1].x_end;
        p_pri = z->w_list.a[z->w_list.n-1].cidx;
    }

    // fprintf(stderr, "+[M::%s] zs::%lu, ze::%lu, zerr::%lu, pe::%lu, is_pri::%lu, p_pri::%lu, wn::%u\n", 
    // __func__, zs, ze, zerr, pe, is_pri, p_pri, (uint32_t)z->w_list.n);

    if(pe <= zs) {///not overlap with [zs, ze)
        if(pe < zs) {
            kv_pushp(window_list, z->w_list, &p); memset(p, 0, sizeof((*p)));
            p->clen = 0; p->x_start = pe; p->x_end = zs; p->cidx = 0;
        }
        kv_pushp(window_list, z->w_list, &p); memset(p, 0, sizeof((*p)));
        p->clen = zerr; p->x_start = zs; p->x_end = ze; p->cidx = is_pri;
    } else {
        if(pe <= ze) {///prefix-suffix overlap
            if(is_pri) {
                if(p_pri) {
                    kv_pushp(window_list, z->w_list, &p); memset(p, 0, sizeof((*p)));
                    p->clen = zerr; p->cidx = is_pri;
                    p->x_start = zs; p->x_end = ze;
                    z->w_list.a[z->w_list.n-2].x_end = zs;///trim the previous primary window
                } else {
                    p = &(z->w_list.a[z->w_list.n-1]);
                    p->clen = zerr; p->cidx = is_pri;
                    p->x_end = ze;
                }
            } else {///if current is not primary
                p = &(z->w_list.a[z->w_list.n-1]);
                p->x_end = ze;
            }
        } else {//pe > ze; [zs, ze) is contained within [ps, pe)
            if(is_pri) {///is_pri == 0, no need to do anything
                if(!p_pri) {
                    p = &(z->w_list.a[z->w_list.n-1]);
                    p->clen = zerr; p->cidx = is_pri;
                } else {
                    kv_pushp(window_list, z->w_list, &p); memset(p, 0, sizeof((*p)));
                    p->clen = zerr; p->cidx = is_pri; 
                    p->x_start = zs; p->x_end = z->w_list.a[z->w_list.n-2].x_end;
                    z->w_list.a[z->w_list.n-2].x_end = zs;///trim the previous primary window
                }
            }
        }
    }
    // fprintf(stderr, "-[M::%s] zs::%lu, ze::%lu, zerr::%lu, pe::%lu, is_pri::%lu, p_pri::%lu, wn::%u, x::[%d, %d), sec::%u\n", 
    // __func__, zs, ze, zerr, pe, is_pri, p_pri, (uint32_t)z->w_list.n, 
    // z->w_list.a[z->w_list.n-1].x_start, z->w_list.a[z->w_list.n-1].x_end, z->w_list.a[z->w_list.n-1].clen);
}

void refine_rphase(ma_ug_t *ug, int64_t rlen, overlap_region *za, uint64_t zan, uint64_t zid, ul_ov_t *a, uint64_t a_n, asg64_v *buf)
{
    // fprintf(stderr, "\n[M::%s] zid::%lu, x::[%u, %u), zan::%lu, a_n::%lu\n", 
    // __func__, zid, za[zid].x_pos_s, za[zid].x_pos_e+1, zan, a_n);
    uint64_t k, bn, qs, qe, err, zs, ze, zerr, m, zwn; overlap_region *z; window_list *p; 
    kv_resize(uint64_t, *buf, a_n);
    for (k = buf->n = 0; k < a_n; k++) {
        if(a[k].qs > 0) {//za[zid] has higher error rate
            // buf->a[buf->n++] = (((uint64_t)za[a[k].tn].x_pos_s)<<32)|((uint64_t)a[k].tn);
            // fprintf(stderr, "+[M::%s] k::%lu, err::%u, x::[%u, %u)\n", __func__, k, a[k].qs, za[a[k].tn].x_pos_s, za[a[k].tn].x_pos_e+1);
            buf->a[buf->n++] = (((uint64_t)za[a[k].tn].x_pos_s)<<32)|(k);
        } 
    }
    radix_sort_bc64(buf->a, buf->a+buf->n);

    bn = buf->n;
    for (k = 0; k < bn; k++) {
        qs = qe = err = 0;
        if(buf->n > bn) {
            qs = buf->a[buf->n-2]>>32; 
            qe = ((uint32_t)buf->a[buf->n-2]);
            err = buf->a[buf->n-1];
        }
        zs = MAX(za[a[((uint32_t)buf->a[k])].tn].x_pos_s, za[zid].x_pos_s);
        ze = MIN(za[a[((uint32_t)buf->a[k])].tn].x_pos_e, za[zid].x_pos_e) + 1;
        if(ze <= zs) continue; 
        zerr = a[((uint32_t)buf->a[k])].qs;
        // fprintf(stderr, "-[M::%s] k::%lu, z::[%u, %u), zerr::%lu\n", __func__, k, zs, ze, zerr);
        if(qe <= zs) {
            kv_push(uint64_t, *buf, ((zs<<32)|ze));
            kv_push(uint64_t, *buf, zerr);
        } else {
            if(ze > qe) qe = ze;
            if(zerr > err) err = zerr;
            buf->a[buf->n-2] = ((qs<<32)|(qe));
            buf->a[buf->n-1] = err;
        }
    }

    for (k = bn, m = 0; k < buf->n; k++) buf->a[m++] = buf->a[k];
    buf->n = m;
    if(!(buf->n)) return;
    
    z = &(za[zid]);
    for (k = 0; k < z->align_length; k++) {
        p = &(z->w_list.a[z->w_list.n+k]);
        if(p->clen <= 0) continue;
        zs = p->x_start; ze = p->x_end; ///note here is [p->x_start, p->x_end)
        zerr = p->clen;
        // fprintf(stderr, ">[M::%s] k::%lu, z::[%u, %u), zerr::%lu\n", __func__, k, zs, ze, zerr);
        kv_push(uint64_t, *buf, ((zs<<32)|ze));
        kv_push(uint64_t, *buf, zerr);
    }

    uint64_t *ref, ref_n, ref_i, *qry, qry_n, qry_i;
    ref = buf->a; ref_n = m; 
    qry = ref + ref_n; qry_n = buf->n - ref_n;
    ref_i = qry_i = z->w_list.n = 0;
    // fprintf(stderr, "[M::%s] ref_n::%lu, qry_n::%lu\n", __func__, ref_n, qry_n);
    ///all ref and qry are regions with at least one error
    while (ref_i < ref_n && qry_i < qry_n) {
        if((ref[ref_i]>>32) < (qry[qry_i]>>32)) {
            push_rphase(z, ref[ref_i]>>32, (uint32_t)ref[ref_i], ref[ref_i+1], 0);
            ref_i += 2;
        } else if((qry[qry_i]>>32) < (ref[ref_i]>>32)) {
            push_rphase(z, qry[qry_i]>>32, (uint32_t)qry[qry_i], qry[qry_i+1], 1);
            qry_i += 2;
        } else {
            push_rphase(z, qry[qry_i]>>32, (uint32_t)qry[qry_i], qry[qry_i+1], 1);
            qry_i += 2;
        }
    }

    while(qry_i < qry_n) {
        push_rphase(z, qry[qry_i]>>32, (uint32_t)qry[qry_i], qry[qry_i+1], 1);
        qry_i += 2;
    }

    while(ref_i < ref_n) {
        push_rphase(z, ref[ref_i]>>32, (uint32_t)ref[ref_i], ref[ref_i+1], 0);
        ref_i += 2;
    }

    ///push remaining bases as the last window    
    qs = qe = z->x_pos_s; zs = z->x_pos_e + 1;
    if(z->w_list.n) {
        qs = z->w_list.a[z->w_list.n-1].x_start;
        qe = z->w_list.a[z->w_list.n-1].x_end;
    }
    if(qe < zs) {
        kv_pushp(window_list, z->w_list, &p); memset(p, 0, sizeof((*p)));
        p->clen = 0; p->x_start = qe; p->x_end = zs;
    }


    z->x_pos_strand = (uint32_t)-1;
    z->overlapLen = z->x_pos_e+1-z->x_pos_s; 
    z->non_homopolymer_errors = 0; zwn = 0;
    for (k = m = 0; k < z->w_list.n; k++) {
        if(z->w_list.a[k].clen > 0) {
            z->non_homopolymer_errors += z->w_list.a[k].clen;
            zwn += z->w_list.a[k].x_end-z->w_list.a[k].x_start;
        } 
        z->w_list.a[k].cidx = 0; m += z->w_list.a[k].x_end-z->w_list.a[k].x_start;
    }
    // fprintf(stderr, "[M::%s] m:%lu, z->overlapLen:%lu\n", __func__, m, z->overlapLen);
    assert(m == z->overlapLen);
    // assert(zwn <= z->overlapLen);///zwn is the length with secondary-best alignment
    z->align_length = z->overlapLen - zwn;

    // fprintf(stderr, ">[M::%s]\tutg%.6u%c\txl::%ld\tx::[%u,\t%u)\t%c\tyl::%u\ty::[%u,\t%u)\tsec::%u\tsec_len::%lu\n", 
    //     __func__, z->y_id + 1, "lc"[ug->u.a[z->y_id].circ], 
    //     rlen, z->x_pos_s, z->x_pos_e+1, "+-"[z->y_pos_strand], 
    //     ug->u.a[z->y_id].len, z->y_pos_s, z->y_pos_e+1, z->non_homopolymer_errors, zwn);
}

void rphase_hl(overlap_region_alloc* ol, const ul_idx_t *uref, const ug_opt_t *uopt, kv_ul_ov_t *c_idx, asg64_v* idx, asg64_v* buf, asg64_v* buf1, int64_t ulid, int64_t bd, 
int64_t rlen, mask_ul_ov_t *mk, idx_emask_t *mm, double len_diff)
{
    uint64_t on0 = ol->length, k, i, l, m0, m1, sec, zwn, m, mmov; 
    overlap_region *z, h; ma_ug_t *ug = uref->ug;

    for (k = i = 0; k < ol->length; k++) {
        z = &(ol->list[k]); 
        z->align_length = 0; 
        z->non_homopolymer_errors = 0;
        z->overlapLen = (uint32_t)-1;
        z->x_pos_strand = 0;
        // fprintf(stderr, "i[M::%s]\tutg%.6u%c\txl::%ld\tx::[%u,\t%u)\t%c\tyl::%u\ty::[%u,\t%u)\tsec::%u\n", 
        // __func__, z->y_id + 1, "lc"[ug->u.a[z->y_id].circ], 
        // rlen, z->x_pos_s, z->x_pos_e+1, "+-"[z->y_pos_strand], 
        // ug->u.a[z->y_id].len, z->y_pos_s, z->y_pos_e+1, z->non_homopolymer_errors);
        if(!(z->w_list.n)) continue;///no base-level; cis
        if(k != i) {
            h = ol->list[k];
            ol->list[k] = ol->list[i];
            ol->list[i] = h;
        }
        i++;
    }
    ol->length = i; ///[ol->length, on0) keeps cis overlaps without the base-level alignment
    
    rphase_rr(ol, uref, uopt, c_idx, idx, buf, ulid, bd, NULL); 
    ///could use idx && buf here
    for (k = idx->n = buf->n = 0; k < ol->length; k++) {
        z = &(ol->list[k]); 
        ///sort overlaps by sec error
        for (i = sec = 0; i < z->align_length; i++) {
            sec += z->w_list.a[z->w_list.n+i].clen;
        }
        z->non_homopolymer_errors = sec; sec = sec - (sec&3);///normalize by 4
        kv_push(uint64_t, *buf, (((uint64_t)sec)<<32)|((uint64_t)k));

        ///sort overlaps by yid for filtering
        kv_push(uint64_t, *idx, (((uint64_t)z->y_id)<<32)|((uint64_t)k));

        // fprintf(stderr, "+[M::%s]\tutg%.6u%c\txl::%ld\tx::[%u,\t%u)\t%c\tyl::%u\ty::[%u,\t%u)\tsec::%u\tn_sec::%lu\n", 
        // __func__, z->y_id + 1, "lc"[ug->u.a[z->y_id].circ], 
        // rlen, z->x_pos_s, z->x_pos_e+1, "+-"[z->y_pos_strand], 
        // ug->u.a[z->y_id].len, z->y_pos_s, z->y_pos_e+1, z->non_homopolymer_errors, sec);
    }
    radix_sort_bc64(idx->a, idx->a+idx->n); ///for filtering

    radix_sort_bc64(buf->a, buf->a+buf->n);///sort by error; smaller error first
    for (l = 0, k = 1; k <= buf->n; k++) {
        if(k == buf->n || (buf->a[l]>>32) == (buf->a[k]>>32)) {///with equal number of normalized errors
            if(((k - l) > 1) && (k < buf->n)) {
                for (i = l; i < k; i++) {
                    // fprintf(stderr, "[M::%s]\tl::%lu\tk::%lu\tbuf->a[i]>>32::%lu\n", __func__, l, k, buf->a[i]>>32);
                    z = &(ol->list[(uint32_t)buf->a[i]]);
                    m0 = z->x_pos_e+1-z->x_pos_s;
                    m1 = z->y_pos_e+1-z->y_pos_s;
                    sec = MIN(m0, m1); sec = ((uint32_t)-1)-((uint32_t)sec);
                    buf->a[i] = (((uint64_t)sec)<<32)|((uint32_t)buf->a[i]);
                }
                radix_sort_bc64(buf->a+l, buf->a+k);///sort by lenght; longer first
            }
            l = k;
        }
    }
    
    mk->idx.n = mk->srt.n = 0; c_idx->n = 0;
    mmov = 1000000; if(mmov > (ol->length*32)) mmov = (ol->length*32);
    for (k = 0; k < buf->n && c_idx->n <= mmov; k++) {///need to add overlap that directly connected in the graph
        z = &(ol->list[(uint32_t)buf->a[k]]); z->overlapLen = c_idx->n; 

        // fprintf(stderr, "\n-[M::%s]\tutg%.6u%c\txl::%ld\tx::[%u,\t%u)\t%c\tyl::%u\ty::[%u,\t%u)\tsec::%u\n", 
        // __func__, z->y_id + 1, "lc"[ug->u.a[z->y_id].circ], 
        // rlen, z->x_pos_s, z->x_pos_e+1, "+-"[z->y_pos_strand], 
        // ug->u.a[z->y_id].len, z->y_pos_s, z->y_pos_e+1, z->non_homopolymer_errors);

        m0 = push_emask_flt(&(mm->a[z->y_id]), idx->a, idx->n, z->y_id, c_idx);


        // fprintf(stderr, "-[M::%s]\tm0::%lu\n", __func__, m0);
        // prt_masks(ug, c_idx->a+c_idx->n-m0, m0);


        m1 = gen_src_shared_interval_simple(z->y_id, ug, idx->a, idx->n, c_idx);


        // fprintf(stderr, "-[M::%s]\tm1::%lu\n", __func__, m1);
        // prt_masks(ug, c_idx->a+c_idx->n-m1, m1);


        dedup_src_shared1(ug, c_idx, m0, m1, mm);
        z->x_pos_strand = c_idx->n - z->overlapLen;
        if(!(z->x_pos_strand)) z->overlapLen = (uint32_t)-1;

        // fprintf(stderr, "-[M::%s]\tutg%.6u%c\txl::%ld\tx::[%u,\t%u)\t%c\tyl::%u\ty::[%u,\t%u)\tsec::%u\tm0::%lu\tm1::%lu\tcan_n::%u\n", 
        // __func__, z->y_id + 1, "lc"[ug->u.a[z->y_id].circ], 
        // rlen, z->x_pos_s, z->x_pos_e+1, "+-"[z->y_pos_strand], 
        // ug->u.a[z->y_id].len, z->y_pos_s, z->y_pos_e+1, z->non_homopolymer_errors, m0, m1, z->x_pos_strand);
        // fprintf(stderr, "-[M::%s]\tcan_n::%u\n", __func__, z->x_pos_strand);
        // if(z->x_pos_strand) prt_masks(ug, c_idx->a+z->overlapLen, z->x_pos_strand);
    }

    ///if no mask overlaps, some overlaps may be still masked using the overlap within the graph
    if(((c_idx->n) || (ol->length > 1)) && mask_ovlps(uref, uopt, ug, ol, mk, c_idx, idx, buf1, mm, len_diff, rlen, bd)) {
        // for (k = 0; k < mk->srt.n; k++) {
        //     fprintf(stderr, "*[M::%s]\tutg%.6u%c\txl::%u\tqerr::%u\t%c\tutg%.6u%c\tyl::%u\tterr::%u\tel::%u\n", 
        //     __func__, 
        //     ol->list[mk->srt.a[k].qn].y_id + 1, "lc"[ug->u.a[ol->list[mk->srt.a[k].qn].y_id].circ], 
        //     ug->u.a[ol->list[mk->srt.a[k].qn].y_id].len, mk->srt.a[k].qs, 
        //     "+-"[mk->srt.a[k].rev],
        //     ol->list[mk->srt.a[k].tn].y_id + 1, "lc"[ug->u.a[ol->list[mk->srt.a[k].tn].y_id].circ], 
        //     ug->u.a[ol->list[mk->srt.a[k].tn].y_id].len, mk->srt.a[k].ts, 
        //     mk->srt.a[k].el);
        // }
        


        rphase_rr(ol, uref, uopt, c_idx, idx, buf, ulid, bd, mk); 
        for (k = 0; k < ol->length; k++) {
            if(((uint32_t)mk->idx.a[k]) > (mk->idx.a[k]>>32)) {///it has masks
                refine_rphase(ug, rlen, ol->list, ol->length, k, mk->srt.a+(mk->idx.a[k]>>32), ((uint32_t)mk->idx.a[k])-(mk->idx.a[k]>>32), idx);
            }
        }
    }
    // assert(on0 <= ol->length);
    ol->length = on0;

    for (k = 0; k < ol->length; k++) {
        z = &(ol->list[k]); 
        if(z->x_pos_strand == (uint32_t)-1) {
        // fprintf(stderr, ">[M::%s]\tutg%.6u%c\txl::%ld\tx::[%u,\t%u)\t%c\tyl::%u\ty::[%u,\t%u)\tsec::%u\n", 
        //     __func__, z->y_id + 1, "lc"[ug->u.a[z->y_id].circ], 
        //     rlen, z->x_pos_s, z->x_pos_e+1, "+-"[z->y_pos_strand], 
        //     ug->u.a[z->y_id].len, z->y_pos_s, z->y_pos_e+1, z->non_homopolymer_errors);

            z->x_pos_strand = 0; continue;
        }
        z->overlapLen = z->x_pos_e+1-z->x_pos_s; 
        z->non_homopolymer_errors = 0; zwn = 0;
        for (i = m = 0; i < z->align_length; i++) {
            z->w_list.a[m] = z->w_list.a[z->w_list.n+i];
            if(z->w_list.a[m].clen > 0) {
                z->non_homopolymer_errors += z->w_list.a[m].clen;
                zwn += z->w_list.a[m].x_end-z->w_list.a[m].x_start;
            } 
            m++;
        }
        z->w_list.n = m; z->x_pos_strand = 0;
        assert(zwn <= z->overlapLen);///zwn is the length with secondary-best alignment
        z->align_length = z->overlapLen - zwn;
        
        // fprintf(stderr, ">[M::%s]\tutg%.6u%c\txl::%ld\tx::[%u,\t%u)\t%c\tyl::%u\ty::[%u,\t%u)\tsec::%u\tsec_len::%lu\n", 
        // __func__, z->y_id + 1, "lc"[ug->u.a[z->y_id].circ], 
        // rlen, z->x_pos_s, z->x_pos_e+1, "+-"[z->y_pos_strand], 
        // ug->u.a[z->y_id].len, z->y_pos_s, z->y_pos_e+1, z->non_homopolymer_errors, zwn);
    }
}

uint64_t rphase_detect0(overlap_region_alloc* ol, kv_ul_ov_t *c_idx, asg64_v* idx, asg64_v* buf, asg64_v* suf)
{
    uint64_t srt_n = idx->n, i, k, h, qs, qe, ts, te, rr = 0, ovlp; overlap_region *z, *p;
    int64_t dp, old_dp, beg, end;
    for (i = k = 0, dp = old_dp = 0, beg = 0, end = -1; i < srt_n; ++i) {///[beg, end) but coordinates in idx is [, ]
        ///if idx->a.a[] is qe
        old_dp = dp;
        if ((idx->a[i]>>32)&1) {
            --dp; end = (idx->a[i]>>33)+1;
        }else {
            //meet a new overlap; the overlaps are pushed by the x_pos_s
            ++dp; end = (idx->a[i]>>33);
            kv_push(uint64_t, *idx, ((uint32_t)idx->a[i]));
        }
        // fprintf(stderr, "\n[M::%s::] beg::%ld, end::%ld, old_dp::%ld\n", __func__, beg, end, old_dp);
        if((end > beg) && (old_dp >= 2)) {
            // fprintf(stderr, "\n[M::%s::] beg::%ld, end::%ld, old_dp::%ld\n", __func__, beg, end, old_dp);
            // kv_resize(uint64_t, *buf, ((uint32_t)old_dp)<<1);
            // idx->n = srt_n + gen_region_phase(ol->list, idx->a+srt_n, idx->n-srt_n, beg, end, old_dp, c_idx->a, buf, buf1);   
            idx->n = srt_n + gen_region_phase_robust(ol->list, idx->a+srt_n, idx->n-srt_n, beg, end, old_dp, c_idx->a, buf);
        }
        beg = end;
    }
    idx->n = srt_n; 
    
    suf->n = 0; kv_resize(uint64_t, *suf, (ol->length<<1)); 
    uint64_t *sep = suf->a, *sid = suf->a + ol->length, sid_n; 
    for (k = h = sid_n = 0; k < ol->length; k++) {
        z = &(ol->list[k]); 
        z->overlapLen = (uint32_t)-1; z->non_homopolymer_errors = 0; 
        for (i = 0; i < z->align_length; i++) {
            z->non_homopolymer_errors += z->w_list.a[z->w_list.n+i].clen;
        }
        for (; h < c_idx->n && ovlp_id(c_idx->a[h]) < k; h++);
        if(h < c_idx->n && ovlp_id(c_idx->a[h]) == k) {
            qs = z->w_list.a[ovlp_min_wid(c_idx->a[h])].x_start;
            qe = z->w_list.a[ovlp_max_wid(c_idx->a[h])].x_end+1;
            for (; h < c_idx->n && ovlp_id(c_idx->a[h]) == k; h++) {
                qe = z->w_list.a[ovlp_max_wid(c_idx->a[h])].x_end+1;
            }
            sid[sid_n] = ((qs<<32)|sid_n);
            sep[sid_n] = ((k<<32)|qe);
            sid_n++;
            // kv_push(uint64_t, *suf, ((qs<<32)|h));
        }
    }
    radix_sort_bc64(sid, sid+sid_n);
    for (k = 0; k < sid_n; k++) {
        qs = sid[k]>>32; qe = ((uint32_t)sep[(uint32_t)sid[k]]); 
        z = &(ol->list[sep[(uint32_t)sid[k]]>>32]);
        if(z->non_homopolymer_errors == 0) continue;
        for (i = rr = 0; i < sid_n; i++) {
            ts = sid[i]>>32; te = ((uint32_t)sep[(uint32_t)sid[i]]); 
            p = &(ol->list[sep[(uint32_t)sid[i]]>>32]);
            if(ts >= qe) break;
            if(te <= qs) continue;
            ovlp = ((MIN(qe, te) > MAX(qs, ts))? (MIN(qe, te) - MAX(qs, ts)):0);
            if(ovlp < 512) continue;
            if(z->non_homopolymer_errors > p->non_homopolymer_errors + rphase_thres) continue;
            rr = 1; break;
        }
        if(rr) z->overlapLen = 1;
    }
    return rr;
}

void reacal_phase(overlap_region_alloc* ol, kv_ul_ov_t *c_idx, asg64_v* idx, asg64_v* buf, asg64_v* suf)
{
    uint64_t srt_n = idx->n, i, k, h, qs, qe, ts, te, rr = 0; overlap_region *z, *p;
    int64_t dp, old_dp, beg, end;
    for (i = k = 0, dp = old_dp = 0, beg = 0, end = -1; i < srt_n; ++i) {///[beg, end) but coordinates in idx is [, ]
        ///if idx->a.a[] is qe
        old_dp = dp;
        if ((idx->a[i]>>32)&1) {
            --dp; end = (idx->a[i]>>33)+1;
        }else {
            //meet a new overlap; the overlaps are pushed by the x_pos_s
            ++dp; end = (idx->a[i]>>33);
            kv_push(uint64_t, *idx, ((uint32_t)idx->a[i]));
        }
        // fprintf(stderr, "\n[M::%s::] beg::%ld, end::%ld, old_dp::%ld\n", __func__, beg, end, old_dp);
        if((end > beg) && (old_dp >= 2)) {
            // fprintf(stderr, "\n[M::%s::] beg::%ld, end::%ld, old_dp::%ld\n", __func__, beg, end, old_dp);
            // kv_resize(uint64_t, *buf, ((uint32_t)old_dp)<<1);
            // idx->n = srt_n + gen_region_phase(ol->list, idx->a+srt_n, idx->n-srt_n, beg, end, old_dp, c_idx->a, buf, buf1);   
            idx->n = srt_n + gen_region_phase_robust(ol->list, idx->a+srt_n, idx->n-srt_n, beg, end, old_dp, c_idx->a, buf);
        }
        beg = end;
    }
    idx->n = srt_n; 
    
    suf->n = 0; kv_resize(uint64_t, *suf, (ol->length<<1)); 
    uint64_t *sep = suf->a, *sid = suf->a + ol->length, sid_n; 
    for (k = h = sid_n = 0; k < ol->length; k++) {
        z = &(ol->list[k]); 
        z->overlapLen = (uint32_t)-1; z->non_homopolymer_errors = 0; 
        for (i = 0; i < z->align_length; i++) {
            z->non_homopolymer_errors += z->w_list.a[z->w_list.n+i].clen;
        }
        for (; h < c_idx->n && ovlp_id(c_idx->a[h]) < k; h++);
        if(h < c_idx->n && ovlp_id(c_idx->a[h]) == k) {
            qs = z->w_list.a[ovlp_min_wid(c_idx->a[h])].x_start;
            qe = z->w_list.a[ovlp_max_wid(c_idx->a[h])].x_end+1;
            for (; h < c_idx->n && ovlp_id(c_idx->a[h]) == k; h++) {
                qe = z->w_list.a[ovlp_max_wid(c_idx->a[h])].x_end+1;
            }
            sid[sid_n] = ((qs<<32)|sid_n);
            sep[sid_n] = ((k<<32)|qe);
            sid_n++;
            // kv_push(uint64_t, *suf, ((qs<<32)|h));
        }
    }
    radix_sort_bc64(sid, sid+sid_n);
    for (k = rr = 0; k < sid_n; k++) {
        qs = sid[k]>>32; qe = ((uint32_t)sep[(uint32_t)sid[k]]); 
        z = &(ol->list[sep[(uint32_t)sid[k]]>>32]);
        if(z->non_homopolymer_errors == 0) continue;
        for (i = 0; i < sid_n; i++) {
            ts = sid[i]>>32; te = ((uint32_t)sep[(uint32_t)sid[i]]); 
            p = &(ol->list[sep[(uint32_t)sid[i]]>>32]);
            if(ts >= qe) break;
            if(te <= qs) continue;
            if(p->non_homopolymer_errors > z->non_homopolymer_errors) continue;
            if(z->non_homopolymer_errors - p->non_homopolymer_errors <= rphase_thres) {
                z->overlapLen = p->overlapLen = 1; rr++;
            }
        }
    }
}

void region_phase_adv(overlap_region_alloc* ol, const ul_idx_t *uref, const ug_opt_t *uopt, kv_ul_ov_t *c_idx, asg64_v* idx, asg64_v* buf, asg64_v* buf1, int64_t ulid)
{
    int64_t on = ol->length, k, i, zwn, q[2], t[2], w[2]; 
    uint64_t m; overlap_region *z; ul_ov_t *cp;
    kv_resize(uint64_t, *idx, (ol->length<<1)); 
    kv_resize(ul_ov_t, *c_idx, ol->length);
    for (k = idx->n = c_idx->n = 0; k < on; k++) {
        z = &(ol->list[k]); zwn = z->w_list.n; 
        z->align_length = 0; z->overlapLen = (uint32_t)-1;
        z->non_homopolymer_errors = 0;
        if(!zwn) continue;
        q[0] = q[1] = t[0] = t[1] = w[0] = w[1] = INT32_MIN;
        for (i = 0; i < zwn; i++) {
            if((z->w_list.a[i].x_start==(q[1]+1)) && ((z->w_list.a[i].y_start==(t[1]+1)))) {
                q[1] = z->w_list.a[i].x_end;
                t[1] = z->w_list.a[i].y_end;
                w[1] = i;
            } else {
                if(q[0] != INT32_MIN) {
                    m = ((uint64_t)q[0])<<1; m <<= 32; 
                    m += c_idx->n; kv_push(uint64_t, *idx, m);
                    m = (((uint64_t)q[1])<<1)+1; m <<= 32; 
                    m += c_idx->n; kv_push(uint64_t, *idx, m);

                    kv_pushp(ul_ov_t, *c_idx, &cp);
                    ovlp_id(*cp) = k; ///ovlp id
                    ovlp_min_wid(*cp) = w[0]; ///beg id of windows
                    ovlp_max_wid(*cp) = w[1]; ///end id of windows
                    ovlp_cur_wid(*cp) = w[0]; ///cur id of windows
                    ovlp_cur_xoff(*cp) = z->w_list.a[w[0]].x_start; ///cur xpos
                    ovlp_cur_coff(*cp) = 0; ///cur cigar off in cur window
                }

                q[0] = z->w_list.a[i].x_start; q[1] = z->w_list.a[i].x_end;
                t[0] = z->w_list.a[i].y_start; t[1] = z->w_list.a[i].y_end;
                w[0] = i; w[1] = i;
            }
        }
        if(q[0] != INT32_MIN) {
            m = ((uint64_t)q[0])<<1; m <<= 32; 
            m += c_idx->n; kv_push(uint64_t, *idx, m);
            m = (((uint64_t)q[1])<<1)+1; m <<= 32; 
            m += c_idx->n; kv_push(uint64_t, *idx, m);

            kv_pushp(ul_ov_t, *c_idx, &cp);
            ovlp_id(*cp) = k; ///ovlp id
            ovlp_min_wid(*cp) = w[0]; ///beg id of windows
            ovlp_max_wid(*cp) = w[1]; ///end id of windows
            ovlp_cur_wid(*cp) = w[0]; ///cur id of windows
            ovlp_cur_xoff(*cp) = z->w_list.a[w[0]].x_start; ///cur xpos
            ovlp_cur_coff(*cp) = 0; ///cur cigar off in cur window
        }
    }
    radix_sort_bc64(idx->a, idx->a+idx->n);
    //this is used with gen_region_phase, now give up
    //gen_gov_idx(ol, uref, uopt, G_CHAIN_BW, N_GCHAIN_RATE, buf1);
    // for (m = 0; m < c_idx->n; m++) {
    //     fprintf(stderr, "+++[M::%s::utg%.6dl] q[%d, %d), t[%d, %d), wn::%d\n", __func__, 
    //     (int32_t)ol->list[ovlp_id(c_idx->a[m])].y_id+1,
    //     ol->list[ovlp_id(c_idx->a[m])].w_list.a[ovlp_min_wid(c_idx->a[m])].x_start,
    //     ol->list[ovlp_id(c_idx->a[m])].w_list.a[ovlp_max_wid(c_idx->a[m])].x_end+1,
    //     ol->list[ovlp_id(c_idx->a[m])].w_list.a[ovlp_min_wid(c_idx->a[m])].y_start,
    //     ol->list[ovlp_id(c_idx->a[m])].w_list.a[ovlp_max_wid(c_idx->a[m])].y_end+1,
    //     ovlp_max_wid(c_idx->a[m])+1-ovlp_min_wid(c_idx->a[m]));
    // }
    
    if(rphase_detect0(ol, c_idx, idx, buf, buf1)) {//refine phasing
        reacal_phase(ol, c_idx, idx, buf, buf1);
    }

    for (k = 0; k < on; k++) {
        z = &(ol->list[k]); 
        z->overlapLen = z->x_pos_e+1-z->x_pos_s; 
        z->non_homopolymer_errors = 0; zwn = 0;
        for (i = m = 0; i < z->align_length; i++) {
            z->w_list.a[m] = z->w_list.a[z->w_list.n+i];
            if(z->w_list.a[m].clen > 0) {
                z->non_homopolymer_errors += z->w_list.a[m].clen;
                zwn += z->w_list.a[m].x_end-z->w_list.a[m].x_start;
            } 
            m++;
        }
        z->w_list.n = m; 
        assert(zwn <= z->overlapLen);
        z->align_length = z->overlapLen - zwn;
    }
}

void ul_gap_filling_adv(overlap_region_alloc* ol, Candidates_list *cl, kv_ul_ov_t *aln, uint64_t wl, 
const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, UC_Read *tu, bit_extz_t *exz, overlap_region *aux_o,  
asg64_v* buf, asg64_v* iidx, double e_rate, int64_t ql, uint64_t rid, int64_t khit, int64_t base_chekc_k_hit,
int64_t max_lgap)
{
    int64_t k, l, ch_n, a_n = aln->n; uint64_t pqn, pk; overlap_region *z; //k_mer_hit *ch_a;
    // count_k_hits(rref, uref, qstr, tu, ol, cl, buf, khit, base_chekc_k_hit);
    count_k_hits_adv(rref, uref, qstr, tu, ol, cl, buf, &(cl->chainDP), e_rate, khit, base_chekc_k_hit);
    // prt_khit(cl, ol, NULL, 109111, "e");
    for (k = 1, l = 0, pqn = 0; k <= a_n; k++) {
        if(k == a_n || aln->a[l].qn != aln->a[k].qn) {
            z = &(ol->list[aln->a[l].qn]); assert(z->align_length == l);  
            // fprintf(stderr, "[M::%s::] oid::[%lu, %u)\n", __func__, pqn, aln->a[l].qn);
            for (pk = pqn; pk < aln->a[l].qn; pk++) ol->list[pk].w_list.n = 0;
            pqn = aln->a[l].qn+1;  
            // fprintf(stderr, "[M::%s::utg%.6dl::%c]\n", __func__, (int32_t)z->y_id+1, "+-"[z->y_pos_strand]);
            // prt_khit(cl, ol, NULL, 109111, "f");
            ch_n = gen_cns_chain(ol, z, cl, iidx, max_lgap, e_rate, 0);
            // prt_khit(cl, ol, NULL, 109111, "g");
            if(ch_n) {
                cigar_gen_by_chain_adv(z, cl, cl->length, ch_n, aln->a+l, k-l, wl, uref, hpc_g, rref, qstr, tu, exz, aux_o, e_rate, ql, rid, khit);
            }
            l = k;
        }
    }
    for (pk = pqn; pk < ol->length; pk++) ol->list[pk].w_list.n = 0;
}

void ul_gap_filling_local(overlap_region_alloc* ol, Candidates_list *cl, kv_ul_ov_t *aln, uint64_t wl, 
const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, UC_Read *tu, bit_extz_t *exz, overlap_region *aux_o,  
asg64_v* buf, asg64_v* iidx, double e_rate, int64_t ql, uint64_t rid, int64_t khit, int64_t base_chekc_k_hit,
int64_t max_lgap)
{
    int64_t k, l, a_n = aln->n; uint64_t pqn, pk; overlap_region *z; //k_mer_hit *ch_a;
    count_k_hits_filter(ol, cl, buf, &(cl->chainDP));
    for (k = 1, l = 0, pqn = 0; k <= a_n; k++) {
        if(k == a_n || aln->a[l].qn != aln->a[k].qn) {
            z = &(ol->list[aln->a[l].qn]); assert(z->align_length == l);  
            for (pk = pqn; pk < aln->a[l].qn; pk++) ol->list[pk].w_list.n = 0;
            pqn = aln->a[l].qn+1;  
            cigar_gen_by_chain_adv_local(z, cl, aln->a+l, k-l, wl, uref, hpc_g, rref, qstr, tu, exz, aux_o, e_rate, ql, rid, khit);
            // ch_n = gen_cns_chain(z, cl, iidx, max_lgap, e_rate, 0);
            // if(ch_n) {
            //     cigar_gen_by_chain_adv(z, cl, cl->length, ch_n, aln->a+l, k-l, wl, uref, hpc_g, rref, qstr, tu, exz, aux_o, e_rate, ql, rid, khit);
            // }
            l = k;
        }
    }
    for (pk = pqn; pk < ol->length; pk++) ol->list[pk].w_list.n = 0;
}

inline uint32_t ovlp_win_check(overlap_region *z, uint32_t id0, uint32_t id1, int64_t max_lgap, double small_bw_rate, int64_t min_small_bw)
{
    if(id0 == (uint32_t)-1 || id1 == (uint32_t)-1) return 1;
    int64_t qs0, qe0, ts0, te0, qs1, qe1, ts1, te1, err, dd, dm, dq, dr;
    if(!get_win_aln(z, id0, &ts0, &te0, &err)) return 1; 
    qs0 = z->w_list.a[id0].x_start; qe0 = z->w_list.a[id0].x_end+1;
    if(!get_win_aln(z, id1, &ts1, &te1, &err)) return 1; 
    qs1 = z->w_list.a[id1].x_start; qe1 = z->w_list.a[id1].x_end+1;

    if(qs1 < qs0 || qe1 < qe0) return 0; 
    if(ts1 < ts0 || te1 < te0) return 0; 
    dq = qe1 - qs0; dr = te1 - ts0; dd = dq>=dr? ((dq)-(dr)): ((dr)-(dq));
    if((ts1 < te0) && (qe0 == qs1)) {//has overlap in y
        dm = dq>=dr?dr:dq;
        if((dd > (dm*small_bw_rate)) && (dd > min_small_bw)) return 0;
    } else {
        if(dd > max_lgap) return 0;
    }
    return 1;
}

uint64_t dp_commen_sketch(kv_ul_ov_t *aln, overlap_region_alloc* ol, uint64_t *id_a, int64_t id_n, 
uint64_t *win_a, int64_t win_n, uint64_t *dp, int64_t n_skip, int64_t wl, int64_t cov, int64_t max_lgap, 
double sgap_rate, int64_t sgap)
{
    if(!win_n) return 0;
    int64_t k, ws, ws0, i, m, wid, wid0, sc, max_sc, p, long_sc, long_idx; overlap_region *z; 
    // fprintf(stderr, "\n[M::%s::] n_skip::%ld, win_n::%ld\n", __func__, n_skip, win_n);
    // for (k = 0; k < win_n; k++) {
    //     fprintf(stderr, "[M::%s::] win_a[%ld]::%lu\n", __func__, k, win_a[k]);
    // }
    if(n_skip < win_n) {
        long_sc = long_idx = -1;
        for (k = 0; k < n_skip; k++) {
            dp[k] = ((k>0)?(k-1):((uint32_t)-1)); dp[k] <<= 32; dp[k] |= k+1; max_sc = k+1;
            if(long_sc < max_sc) {
                long_sc = max_sc; long_idx = k;
            }
        }

        for (; k < win_n; k++) {
            ws = win_a[k]; sc = 1; p = -1; max_sc = sc;
            for (m = k-1; m >= 0; m--) {
                ws0 = win_a[m];
                for (i = 0; i < id_n; i++) {
                    z = &(ol->list[aln->a[(uint32_t)id_a[i]].qn]);
                    wid = get_win_id_by_s(z, ws, wl, NULL);
                    wid0 = get_win_id_by_s(z, ws0, wl, NULL);
                    if(!ovlp_win_check(z, wid0, wid, max_lgap, sgap_rate, sgap)) break;
                }
                if((i >= id_n) && ((sc + ((uint32_t)dp[m])) > max_sc)) {
                    max_sc = (sc + ((uint32_t)dp[m])); p = m;
                }
            }
            dp[k] = ((p>=0)?(p):((uint32_t)-1)); dp[k] <<= 32; dp[k] |= max_sc;
            if(long_sc < max_sc) {
                long_sc = max_sc; long_idx = k;
            }
        }
        // fprintf(stderr, "[M::%s::] long_idx::%lu, long_sc::%lu\n", __func__, long_idx, long_sc);
        k = long_idx; 
        while (k >= 0) {
            win_a[k] |= ((uint64_t)0x8000000000000000); 
            k = ((dp[k]>>32)!=((uint32_t)-1))?(dp[k]>>32):(-1);
        }

        for (k = 0, m = 0; k < win_n; k++) {
            if(!(win_a[k]&((uint64_t)0x8000000000000000))) continue;
            win_a[m++] = (win_a[k]<<1)>>1;
        }
        win_n = m;
        // fprintf(stderr, "[M::%s::] win_n::%lu\n", __func__, win_n);
    }

    for (k = 0, cov = -cov; k < win_n; k++) {
        ws = win_a[k];
        for (i = 0; i < id_n; i++) {
            z = &(ol->list[aln->a[(uint32_t)id_a[i]].qn]);
            wid = get_win_id_by_s(z, ws, wl, NULL);
            // if(wid >= z->w_list.n) {
            //     fprintf(stderr, "[M::%s::] ws::%ld, wl::%ld, wid::%ld, wn::%ld\n", 
            //         __func__, ws, wl, wid, ((int64_t)z->w_list.n));
            // }
            if(cov > INT16_MIN) z->w_list.a[wid].extra_end = cov; 
            else z->w_list.a[wid].extra_end = INT16_MIN;
            z->align_length = wid;///for conliner
        }
    }
    return win_n;
}


uint64_t gen_commen_sketch(All_reads *rref, const ul_idx_t *uref, overlap_region_alloc* ol, uint64_t *id_a, uint64_t id_n, uint64_t s, uint64_t e, uint64_t ql, uint64_t wl,
uint64_t *buf, uint64_t dp, char *str0, char *str1, kv_ul_ov_t *aln, asg64_v *trace,
int64_t max_lgap, double sgap_rate, int64_t sgap)///[s, e)
{
    if(!id_n) return id_n;
    uint64_t i, m, k, rm_n = 0, buf_n = 0, qs, qe, wid, co, occ; char *qstring, *tstring;
    overlap_region *z; uint64_t ws, we; int64_t r_y[2], r_err, p_y[2], p_err, rxl, ryl, pxl;
    ///shrink [qs, qe)
    qs = (s/wl)*wl; if(qs < s) qs += wl; if(qs >= ql) return id_n;
    qe = (e/wl)*wl; if(qe >= ql) qe = ql;
    if(qs >= qe) return id_n;
    //idx_a[] is sorted by aln[].qs
    for (k = 0; k < id_n; k++) {
        if(aln->a[id_a[k]].qs<=qs && aln->a[id_a[k]].qe>=qe) {
            buf[buf_n] = ol->list[aln->a[id_a[k]].qn].align_length; buf[buf_n] <<= 32; buf[buf_n] |= id_a[k];
            buf_n++;
        }
        if(aln->a[id_a[k]].qe < e) rm_n++;
    }
    assert(buf_n == dp && buf_n > 1);

    if(buf_n > 0) {
        ///fs = fe = (uint64_t)-1;
        co = 1; trace->n = occ = 0;
        for (k = qs; k < qe; k += wl) {
            ws = k; we = ws + wl; if(we > qe) we = qe;//[ws, we)
            ///first overlap
            z = &(ol->list[aln->a[(uint32_t)buf[0]].qn]);
            wid = get_win_id_by_s(z, ws, wl, NULL);
            if(!get_win_aln(z, wid, &(r_y[0]), &(r_y[1]), &r_err)) continue;
            if(!ovlp_win_check(z, z->align_length, wid, max_lgap, sgap_rate, sgap)) continue;///not co-linear
            rxl = z->w_list.a[wid].x_end+1-z->w_list.a[wid].x_start;
            ryl = r_y[1]-r_y[0];
            // fprintf(stderr, "+[M::%s::] buf_n::%lu, q::[%lu, %lu), w::[%lu, %lu), rxl::%ld, ryl::%ld, x::[%d, %d), y::[%ld, %ld)\n", 
            //         __func__, buf_n, qs, qe, ws, we, rxl, ryl, z->w_list.a[wid].x_start, z->w_list.a[wid].x_end+1, r_y[0], r_y[1]);
            if(rxl > 1 && ryl > 1) continue;///length of window should be longer than 1
            // if(rxl <= 1 || ryl <= 1) continue;///length of window should be longer than 1
            // fprintf(stderr, "-[M::%s::] buf_n::%lu, q::[%lu, %lu), w::[%lu, %lu), rxl::%ld, ryl::%ld, x::[%d, %d), y::[%ld, %ld)\n", 
            //         __func__, buf_n, qs, qe, ws, we, rxl, ryl, z->w_list.a[wid].x_start, z->w_list.a[wid].x_end+1, r_y[0], r_y[1]);
            qstring = tstring = NULL;
            for (i = 1; i < buf_n; i++) {
                z = &(ol->list[aln->a[(uint32_t)buf[i]].qn]);
                wid = get_win_id_by_s(z, ws, wl, NULL);
                if(!get_win_aln(z, wid, &(p_y[0]), &(p_y[1]), &p_err)) break;
                if(!ovlp_win_check(z, z->align_length, wid, max_lgap, sgap_rate, sgap)) break;///not co-linear
                pxl = z->w_list.a[wid].x_end+1-z->w_list.a[wid].x_start;
                if(pxl != rxl) break;
                if(((r_y[1]-r_y[0]) != (p_y[1]-p_y[0])) || (r_err != p_err)) break;
                if(r_err == 0) continue;
                if(!qstring) {
                    qstring = retrive_str_piece_exz(rref, uref, str0, r_y[0], r_y[1]-r_y[0], 
                                                ol->list[aln->a[(uint32_t)buf[0]].qn].y_pos_strand, ol->list[aln->a[(uint32_t)buf[0]].qn].y_id);
                }
                tstring = retrive_str_piece_exz(rref, uref, str1, p_y[0], p_y[1]-p_y[0], z->y_pos_strand, z->y_id);
                if(memcmp(qstring, tstring, (p_y[1]-p_y[0]))) {
                    // fprintf(stderr, "[M::%s::] qs::%ld, qe::%lu, ts::%ld, te::%lu, err::%ld, rts::%ld, rte::%lu, err::%ld\n", 
                    // __func__, ws, we, p_y[0], p_y[1], r_err, r_y[0], r_y[1], p_err);
                    // fprintf(stderr, "str0::%.*s\n", (int32_t)(r_y[1]-r_y[0]), str0);
                    // fprintf(stderr, "str1::%.*s\n", (int32_t)(p_y[1]-p_y[0]), str1);
                    break;
                }
            }
            if(i < buf_n) continue;
            if(co) {
                for (i = 0; i < buf_n; i++) {
                    z = &(ol->list[aln->a[(uint32_t)buf[i]].qn]);
                    wid = get_win_id_by_s(z, ws, wl, NULL);
                    co = ovlp_win_check(z, buf[i]>>32, wid, max_lgap, sgap_rate, sgap);
                    m = wid; m <<= 32; m |= (uint32_t)buf[i]; buf[i] = m;
                }
                if(co) occ++;
            }
            kv_push(uint64_t, *trace, ws);
        }

        if(trace->n) {
            if(!co) kv_resize(uint64_t, *trace, trace->n<<1);
            trace->n = dp_commen_sketch(aln, ol, buf, buf_n, trace->a, trace->n, trace->a + trace->n, occ, wl, dp, max_lgap, sgap_rate, sgap);
        }
    }
    if(rm_n) {
        for (i = m = 0; i < id_n; i++) {
            if(aln->a[id_a[i]].qe < e) continue;
            id_a[m++] = id_a[i];
        }
        id_n = m;
    }
    return id_n;
}

void update_sketch_trace(overlap_region_alloc* ol, const ul_idx_t *uref, const ug_opt_t *uopt, 
All_reads *rref, UC_Read* tu, asg64_v* idx, asg64_v *b0, asg64_v *b1, int64_t ql, int64_t wl, 
kv_ul_ov_t *aln, uint64_t rid, int64_t max_lgap, double sgap_rate)
{
    idx->n = 0;
	if(!aln->n) return;
	uint64_t i, k, own, srt_n; int64_t dp, old_dp, beg, end; overlap_region *z;
    for (i = 0; i < ol->length; i++) {
		z = &(ol->list[i]); append_unmatched_wins(z, wl); 
		own = z->w_list.n; z->align_length = (uint32_t)-1;
		for (k = 0; k < own; k++) {
			if(z->w_list.a[k].extra_end < 0) z->w_list.a[k].extra_end = 0;
		}		
    }

    kv_resize(uint64_t, *idx, (aln->n<<1)); kv_resize(uint64_t, *b0, aln->n);
    for (i = srt_n = 0; i < aln->n; i++) {
		// if(i == 0 || aln->a[i].qn != aln->a[i-1].qn) ol->list[aln->a[i].qn].align_length = i;
        ol->list[aln->a[i].qn].align_length = (uint32_t)-1;///for co-linear
        idx->a[srt_n] = aln->a[i].qs<<1; idx->a[srt_n] <<= 32; idx->a[srt_n] += i; srt_n++;
        idx->a[srt_n] = ((aln->a[i].qe-1)<<1)+1; idx->a[srt_n] <<= 32; idx->a[srt_n] += i; srt_n++;
		aln->a[i].el = 1;
    }

    radix_sort_bc64(idx->a, idx->a+srt_n); idx->n = srt_n; resize_UC_Read(tu, (wl<<1));
    for (i = 0, dp = 0, beg = 0, end = -1; i < srt_n; ++i) {///[beg, end), but the idx saves [qs, qe]
		old_dp = dp;
        ///if idx->a.a[] is qe
        if ((idx->a[i]>>32)&1) {
			--dp; end = (idx->a[i]>>33)+1;
		}else {
			//meet a new overlap; the overlaps are pushed by the x_pos_s
			++dp; end = (idx->a[i]>>33);
			kv_push(uint64_t, *idx, ((uint32_t)idx->a[i]));
		}
		
		if((end > beg) && (end - beg > wl) && (old_dp >= 2) ) {
            // fprintf(stderr, "+++[M::%s::] beg::%ld, end::%ld\n", __func__, beg, end);
			idx->n = srt_n + 
					gen_commen_sketch(rref, uref, ol, idx->a+srt_n, idx->n-srt_n, beg, end, ql, wl, b0->a, 
                    old_dp, tu->seq, tu->seq+wl, aln, b1, max_lgap, sgap_rate, SGAP);
		}
        beg = end;
	}

    for (i = 0; i < aln->n; i++) {
		if(i == 0 || aln->a[i].qn != aln->a[i-1].qn) ol->list[aln->a[i].qn].align_length = i;
    }
    idx->n = srt_n;
    // fprintf(stderr, "+++[M::%s::] idx->n::%ld\n", __func__, idx->n);
	return;
}



void ul_lalign(overlap_region_alloc* ol, Candidates_list *cl, const ul_idx_t *uref, const ug_opt_t *uopt, 
        char *qstr, uint64_t ql, UC_Read* qu, UC_Read* tu, Correct_dumy* dumy, bit_extz_t *exz, haplotype_evdience_alloc* hap, 
        kvec_t_u64_warp* v_idx, overlap_region *aux_o, double e_rate, int64_t wl, kv_ul_ov_t *aln, int64_t sid, uint64_t khit, 
        st_mt_t *stb, idx_emask_t *mm, mask_ul_ov_t *mk, void *km)
{
    uint64_t i, bs, k, ovl/**, on**/; Window_Pool w; double err; 
    /**int64_t sc;**/ overlap_region t; overlap_region *z; asg64_v iidx, buf, buf1;
    ol->mapped_overlaps_length = 0;
    if(ol->length <= 0) return;

    ///base alignment
    clear_Correct_dumy(dumy, ol, km); err = e_rate; 
    init_Window_Pool(&w, ql, wl, (int)(1.0/err));
    bs = (w.window_length)+(THRESHOLD_MAX_SIZE<<1)+1;
    resize_UC_Read(tu, bs<<1);

    if(!aln) {    
        resize_UC_Read(qu, ql); qu->length = ql; memcpy(qu->seq, qstr, ql);
        for (i = 0; i < ol->length; i++) {
            z = &(ol->list[i]); ovl = z->x_pos_e+1-z->x_pos_s; 
            z->shared_seed = z->non_homopolymer_errors;///for index
            if(!align_ul_ed_post_extz(z, uref, NULL, qu->seq, tu->seq, exz, err, w.window_length, -1, 0, km)) {
                continue;
            }
            if(uref && simi_pass(ovl, z->align_length, uref?1:0, -1, NULL)) {
                z->is_match = 3; ol->mapped_overlaps_length += z->align_length;
            }
        }

        if(uref && ol->mapped_overlaps_length > 0) {
            set_herror_win(ol, dumy, v_idx, err, ql, w.window_length);
        }

        double e_max = err*1.5, rr; int64_t re;
        for (i = k = 0; i < ol->length; i++) {
            z = &(ol->list[i]); ovl = z->x_pos_e + 1 - z->x_pos_s;
            rr = gen_extend_err_exz(z, uref, NULL, NULL, qu->seq, tu->seq, exz, v_idx?v_idx->a.a:NULL, w.window_length, -1, err, (e_max+0.000001), &re);
            z->is_match = 0;///must be here;


            // fprintf(stderr, "[M::%s::utg%.6dl::%c] sid::%ld, q::[%d, %d), t::[%d, %d), err::%ld\n", 
            //     __func__, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], sid, z->x_pos_s, z->x_pos_e+1, z->y_pos_s, z->y_pos_e+1,
            //     re);

            if (rr <= err) {
                if(k != i) {
                    t = ol->list[k];
                    ol->list[k] = ol->list[i];
                    ol->list[i] = t;
                }
                ol->list[k].is_match = 1; ol->list[k].non_homopolymer_errors = re;
                // fprintf(stderr, "+[M::%s] on::%lu\n", __func__, ol->length);
                /****for debug****/
                // z = &(ol->list[k]);
                // fprintf(stderr, "[M::%s::utg%.6dl::%c] sid::%ld, q::[%d, %d), t::[%d, %d), err::%u\n", 
                // __func__, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], sid, z->x_pos_s, z->x_pos_e+1, z->y_pos_s, z->y_pos_e+1,
                // z->non_homopolymer_errors);
                /****for debug****/
                k++;
            } 
        }

        ol->length = k;
        // fprintf(stderr, "+[M::%s] on::%lu\n", __func__, ol->length);
        if(ol->length <= 0) return;
    } else {
        // fprintf(stderr, "-[M::%s] on::%lu\n", __func__, ol->length);
        if(ol->length <= 1) {
            if(ol->length == 1) {
                window_list *p;
                z = &(ol->list[0]); z->w_list.n = 0;
                kv_pushp(window_list, z->w_list, &p);
                memset(p, 0, sizeof((*p)));
                p->x_start = z->x_pos_s; p->x_end = z->x_pos_e+1; p->clen = 0;
                z->align_length = z->overlapLen = z->x_pos_e+1-z->x_pos_s;
                z->non_homopolymer_errors = 0;
            }
            return;
        }
        ///coordinates for all intervals with cov > 1
        copy_asg_arr(iidx, hap->snp_srt); copy_asg_arr(buf, v_idx->a); copy_asg_arr(buf1, (*stb));
        // fprintf(stderr, "\n[M::%s] iidx_n::%ld\n", __func__, (int64_t)iidx.n);
        ul_gap_filling_adv(ol, cl, aln, wl, uref, NULL, NULL, qu->seq, tu, exz, aux_o, &buf, &iidx, err, ql, sid, khit, 1, MAX_LGAP(ql));
        copy_asg_arr(hap->snp_srt, iidx); copy_asg_arr(v_idx->a, buf); copy_asg_arr((*stb), buf1);

        copy_asg_arr(iidx, hap->snp_srt); copy_asg_arr(buf, v_idx->a); copy_asg_arr(buf1, (*stb));
        // region_phase(ol, uref, uopt, aln, &iidx, &buf, &buf1, sid);
        rphase_hl(ol, uref, uopt, aln, &iidx, &buf, &buf1, sid, 64, ql, mk, mm, err);
        copy_asg_arr(hap->snp_srt, iidx); copy_asg_arr(v_idx->a, buf); copy_asg_arr((*stb), buf1);
    }
}


int64_t get_chain_x_by_y(overlap_region* ot, int64_t q)
{
	int64_t x, y, off, i, lx = -1, ly = -1; Fake_Cigar* o = &(ot->f_cigar);
	x = get_fake_gap_pos(o, o->length - 1); 
	off = get_fake_gap_shift(o, o->length - 1);
	y = x - ot->x_pos_s + ot->y_pos_s + off;
	if(y == q) return x;

	for (i = 0; i < (int64_t)o->length; i++){
		x = get_fake_gap_pos(o, i); off = get_fake_gap_shift(o, i);
		y = x - ot->x_pos_s + ot->y_pos_s + off;
		if(q < y) {
			lx = x; ly = y;
			break;
		}
	}

    assert((i!=0)&&(i!=(int64_t)o->length));
	x = get_fake_gap_pos(o, i-1); off = get_fake_gap_shift(o, i-1);
	y = x - ot->x_pos_s + ot->y_pos_s + off;
	y = (((double)(q - y))/((double)(ly - y)))*((double)(lx -x)) + x;
	if(y < ot->x_pos_s) y = ot->x_pos_s;
	if(y > ot->x_pos_e) y = ot->x_pos_e;
	return y;
}

int64_t gen_contain_ov(const ul_idx_t *uref, utg_ct_t *p, overlap_region* o, kv_ul_ov_t *res)
{
	int64_t y_s, y_e, y_bs, y_be, x_s, x_e, q_s, q_e;
	if(o->y_pos_strand) {
		y_s = uref->ug->u.a[o->y_id].len - p->e;
		y_e = uref->ug->u.a[o->y_id].len - p->s - 1;
	} else {
		y_s = p->s; y_e = p->e - 1;
	}

	y_s = MAX(y_s, (int64_t)o->y_pos_s); y_e = MIN(y_e, (int64_t)o->y_pos_e);
	if(y_s > y_e) return 0;
	x_s = get_chain_x_by_y(o, y_s); x_e = get_chain_x_by_y(o, y_e) + 1;
    assert(x_s < x_e);
	// if(x_s >= x_e) fprintf(stderr, "+++y_s->%ld, y_e->%ld, x_s->%ld, x_e->%ld\n", y_s, y_e, x_s, x_e);
	if(o->y_pos_strand) {
		y_bs = uref->ug->u.a[o->y_id].len - (y_e+1);
		y_be = uref->ug->u.a[o->y_id].len - y_s;
	} else {
		y_bs = y_s; y_be = y_e + 1;
	}

	q_s = 0; q_e = p->e - p->s;
	if(p->x&1) {
		q_s += (p->e - y_be);
		q_e -= (y_bs - p->s);
	} else {
		q_s += (y_bs - p->s);
		q_e -= (p->e - y_be); 
	}

	ul_ov_t *x = NULL;
	kv_pushp(ul_ov_t, *res, &x);
	x->qn = o->x_id; x->qs = x_s; x->qe = x_e;
	x->tn = (uint32_t)(0x80000000); x->tn |= (p->x>>1);
	x->ts = q_s; x->te = q_e; x->el = 1;x->sec = 0; x->rev = ((o->y_pos_strand == (p->x&1))?0:1);
	return 1;
}


uint64_t gen_sub_ov(const ul_idx_t *udb, overlap_region* o, kv_ul_ov_t *res, uint64_t xs, uint64_t xe)
{
	uint64_t ts, te, i, l, rn = res->n, rev = o->y_pos_strand; ul_ov_t *z;
	utg_ct_t p; ma_utg_t *u = &(udb->ug->u.a[o->y_id]);
    if(!rev){
        ts = o->y_pos_s; te = o->y_pos_e + 1;
    } else {
        ts = u->len - (o->y_pos_e+1); te = u->len - o->y_pos_s;
    }

    for (i = l = 0; i < u->n; i++) {
        p.x = u->a[i]>>32; p.s = l; p.e = l + Get_READ_LENGTH(R_INF, (u->a[i]>>33));
        l += (uint32_t)u->a[i];
        if(p.e <= ts) continue; 
        if(p.s >= te) break;
        if(gen_contain_ov(udb, &p, o, res)) {
            z = &(res->a[res->n-1]);
            if(z->qs >= xe) break;
            if(z->qe <= xs) res->n--; 
        }
    }
    return res->n-rn;
}

uint64_t gen_conta_ov(const ul_idx_t *udb, overlap_region* o, utg_ct_t *ct_a, int64_t ct_n, kv_ul_ov_t *res, uint64_t xs, uint64_t xe)
{
	int64_t i; uint64_t ys, ye, rn = res->n; ma_utg_t *u = &(udb->ug->u.a[o->y_id]);
     ul_ov_t *z; utg_ct_t *p;

    if(o->y_pos_strand == 0){
        ys = o->y_pos_s; ye = o->y_pos_e + 1;
    } else {
        ys = u->len - (o->y_pos_e+1); ye = u->len - o->y_pos_s;
    }

    for (i = 0; i < ct_n; i++) {
        p = &(ct_a[i]);
        if(p->e <= ys) continue; 
        if(p->s >= ye) break;
        if(gen_contain_ov(udb, p, o, res)) {
            z = &(res->a[res->n-1]);
            if(z->qs >= xe) break;
            if(z->qe <= xs) res->n--; 
        }
    }
    return res->n-rn;
}


uint64_t gen_r_aln(const ul_idx_t *udb, overlap_region *z, uint64_t zid, asg64_v *idx, double o_rate, kv_ul_ov_t *aln, uint64_t min_ovlp)
{
    uint64_t ol = z->x_pos_e+1-z->x_pos_s, aln_ol = z->align_length, k, cn, is_srt; 
    ul_ov_t *p; ul_contain *ct = udb->ct; utg_ct_t *ca; uint64_t os, oe, ovlp, salnl, sol, raln = aln->n;
    int64_t sk, ek, wn, s, e, ws, we, minw, maxw;
    if((aln_ol < min_ovlp) || (!(z->w_list.n))) return 0;
    if((ol*o_rate) <= aln_ol) {
        kv_pushp(ul_ov_t, *aln, &p);
        p->qn = zid; p->tn = z->y_id; p->el = 1; p->rev = z->y_pos_strand; p->sec = 0;
        p->qs = z->x_pos_s; p->qe = z->x_pos_e+1; //[qs, qe)
        p->ts = 0; p->te = z->w_list.n; //[ts, te)
    } else {
        cn = ((uint32_t)(ct->idx.a[z->y_id])); 
        ca = ct->rids.a + ((ct->idx.a[z->y_id])>>32);
        wn = z->w_list.n; 
        ws = z->w_list.a[0].x_start; 
        we = z->w_list.a[z->w_list.n-1].x_end+1;

        gen_sub_ov(udb, z, aln, ws, we); gen_conta_ov(udb, z, ca, cn, aln, ws, we);

        idx->n = aln->n-raln; kv_resize(uint64_t, *idx, idx->n);
        for (k = raln, idx->n = 0, is_srt = 1; k < aln->n; k++) {
            if(k > raln && aln->a[k].qs < aln->a[k-1].qs) is_srt = 0;
            idx->a[idx->n++] = (((uint64_t)(aln->a[k].qs))<<32)|((uint64_t)(aln->a[k].qe));
        }
        if(!is_srt) radix_sort_bc64(idx->a, idx->a+idx->n);
        aln->n = raln;

        for (k = sk = ek = 0; k < idx->n; k++) {///idx->a is sorted by s, not by e
            s = idx->a[k]>>32; e = (uint32_t)idx->a[k]; 
            salnl = 0; sol = e - s; minw = INT32_MAX; maxw = 0;
            for (sk = ((sk<wn)?(sk):(wn-1)); sk < wn && s > z->w_list.a[sk].x_end; sk++); 
            for (sk = ((sk<wn)?(sk):(wn-1)); sk >= 0 && s < z->w_list.a[sk].x_start; sk--); 
            if(k < 0) k = 0; ///s >= z->w_list.a[sk].x_start && s <= z->w_list.a[sk].x_end            
            for (ek = sk; ek < wn; ek++) {
                ws = z->w_list.a[ek].x_start;
                we = z->w_list.a[ek].x_end + 1;
                if(ws >= e) break;
                if(z->w_list.a[ek].y_end == -1) continue;
                os = MAX(s, ws); oe = MIN(e, we);
                ovlp = ((oe>os)? (oe-os):0);
                if(!ovlp) continue;
                salnl += ovlp; 
                if(ek < minw) minw = ek;
                if(ek > maxw) maxw = ek;
            }
            p = NULL; 
            if(((sol*o_rate)<=salnl) && (maxw >= minw)) {
                ws = z->w_list.a[minw].x_start; we = z->w_list.a[maxw].x_end + 1; maxw++;
                if(aln->n > raln) {
                    p = &(aln->a[aln->n-1]);
                    os = MAX(s, p->qs); oe = MIN(e, p->qe);
                    if(oe>os) {
                        if(s < p->qs) p->qs = s;
                        if(e > p->qe) p->qe = e;
                        if(minw < p->ts) p->ts = minw;
                        if(maxw > p->te) p->te = maxw;
                    } else {
                        p = NULL;
                    }
                }
                if(!p) {
                    kv_pushp(ul_ov_t, *aln, &p);
                    p->qn = zid; p->tn = z->y_id; p->el = 0; 
                    p->rev = z->y_pos_strand; p->sec = 0;
                    p->qs = s; p->qe = e; ///[qs, qe)
                    p->ts = minw; p->te = maxw; ///[ts, te)
                }
            }
        }
    }

    return aln->n - raln;
}

///q[2], t[2]
int64_t get_win_yoff(overlap_region *o, int64_t toff, int64_t k, int64_t qs, int64_t qe, int64_t ts, int64_t te, int64_t *q, int64_t *t,
int64_t *rq, int64_t *rt)
{
    int64_t wn = o->w_list.n<<1, s, qdis, tdis, dis;
    if(k<0) k = 0; if(k>=wn) k = wn-1;
    for (; k < wn; k++) {
        s = ((k&1)?(o->w_list.a[k>>1].x_end):(o->w_list.a[k>>1].x_start));
        if(s >= toff) break;
    }
    for (k=(k>=wn?(wn-1):(k)); k >= 0; k--) {
        s = ((k&1)?(o->w_list.a[k>>1].x_end):(o->w_list.a[k>>1].x_start));
        if(s <= toff) break;
    }

    q[0] = qs; q[1] = qe; t[0] = ts; t[1] = te;
    if(k < 0) {//toff <= t[1]
        q[1] =  o->w_list.a[0].x_start;
        t[1] = o->w_list.a[0].y_start;
        tdis = t[1] - toff; qdis = q[1] - q[0]; dis = MIN(qdis, tdis);
        (*rq) = q[1] - dis; (*rt) = t[1] - dis;
    } else if(k == wn-1) {//toff >= t[0]
        q[0] = o->w_list.a[k>>1].x_end;
        t[0] = o->w_list.a[k>>1].y_end;
        tdis = toff - t[0]; qdis = q[1] - q[0]; dis = MIN(qdis, tdis);
        (*rq) = q[0] + dis; (*rt) = t[0] + dis;
    } else {//toff >= t[0] && toff <= t[1]
        q[0] = ((k&1)?(o->w_list.a[k>>1].x_end):(o->w_list.a[k>>1].x_start));
        t[0] = ((k&1)?(o->w_list.a[k>>1].y_end):(o->w_list.a[k>>1].y_start));
        k++;
        q[1] = ((k&1)?(o->w_list.a[k>>1].x_end):(o->w_list.a[k>>1].x_start));
        t[1] = ((k&1)?(o->w_list.a[k>>1].y_end):(o->w_list.a[k>>1].y_start));
        k--;
        (*rt) = toff;
        (*rq) = q[0] + get_offset_adjust(toff-t[0], t[1]-t[0], q[1]-q[0]);
    }
    return k;
}

int64_t get_win_off(int64_t ql, utg_ct_t *p, overlap_region* o, double o_rate, ul_ov_t *res)
{
    if(!o->w_list.n) return 0;
    int64_t q[2], t[2], wk[2], wq[2], wt[2], k, tot_l, qs, qe, ts, te, os, oe, ovlp, aln_l;
    qs = 0; qe = ql-1; ts = p->s; te = p->e-1; //[qs, qe] && [ps, pe]

    k = 0; memset(res, 0, sizeof((*res)));
    k = get_win_yoff(o, ts, k, qs, qe, ts, te, q, t, &(wq[0]), &(wt[0])); wk[0] = k; 
    k = get_win_yoff(o, te, k, qs, qe, ts, te, q, t, &(wq[1]), &(wt[1])); wk[1] = k;
    res->ts = ts; res->te = te + 1;///offset of t
    res->qs = wk[0]<0?(uint32_t)-1:wk[0]; res->qe = wk[1] + 1;///id of window

    wk[0] = ((wk[0]>=0)?(wk[0]>>1):0);
    wk[1] = ((wk[1]>=0)?(wk[1]>>1):0)+1;
    tot_l = wq[1]+1-wq[0]; aln_l = 0;
    for (k = wk[0]; k < wk[1]; k++) {
        if(is_ualn_win(o->w_list.a[k])) continue;
        os = MAX(wq[0], o->w_list.a[k].x_start); 
        oe = MIN(wq[1], o->w_list.a[k].x_end) + 1;
        ovlp = ((oe>os)? (oe-os):0); aln_l += ovlp;
    }
    if((aln_l < (tot_l*o_rate)) || (aln_l == 0)) return 0;
	return 1;
}

uint64_t gen_sub_ov_cigar(const ul_idx_t *udb, uint64_t ql, overlap_region* o, double o_rate, kv_ul_ov_t *res)
{
	uint64_t ts, te, i, l, tl, id = o->y_id, rn = res->n, s, e, rev = o->y_pos_strand; utg_ct_t p; 
    ma_utg_t *u = &(udb->ug->u.a[id]); ul_ov_t rr; int64_t wn, k;

    wn = o->w_list.n; tl = u->len;
    for (k = 0; (k < wn) && (is_ualn_win(o->w_list.a[k])); k++);
    if(k >= wn) return 0;
    if(!rev) {
        ts = o->w_list.a[k].x_start;
    } else {
        te = tl-o->w_list.a[k].x_start;
    }
    for (k = wn-1; (k >= 0) && (is_ualn_win(o->w_list.a[k])); k--);
    if(k < 0) return 0;
    if(!rev) {
        te = o->w_list.a[k].x_end+1;
    } else {
        ts = tl-o->w_list.a[k].x_end-1;
    }
    if(ts >= te) return 0;///[ts, te)
    

    for (i = l = 0; i < u->n; i++) {
        p.x = u->a[i]>>32; s = l; e = l + Get_READ_LENGTH(R_INF, (u->a[i]>>33));
        l += (uint32_t)u->a[i];
        if(e <= ts) continue; if(s >= te) break;
        if(!rev) {
            p.s = s; p.e = e;
        }
        else {
            p.s = tl - e; p.e = tl - s;
        }

        if(get_win_off(ql, &p, o, o_rate, &rr)) {
            rr.el = 0; kv_push(ul_ov_t, *res, rr);
        }
    }
    return res->n-rn;
}


uint64_t gen_conta_ov_cigar(const ul_idx_t *udb, uint64_t ql, overlap_region* o, utg_ct_t *ct_a, uint64_t ct_n, double o_rate, kv_ul_ov_t *res)
{
	uint64_t ts, te, i, tl, id = o->y_id, rn = res->n, rev = o->y_pos_strand; utg_ct_t p, *z; 
    ma_utg_t *u = &(udb->ug->u.a[id]); ul_ov_t rr; int64_t wn, k;

    wn = o->w_list.n; tl = u->len;
    for (k = 0; (k < wn) && (is_ualn_win(o->w_list.a[k])); k++);
    if(k >= wn) return 0;
    if(!rev) {
        ts = o->w_list.a[k].x_start;
    } else {
        te = tl-o->w_list.a[k].x_start;
    }
    for (k = wn-1; (k >= 0) && (is_ualn_win(o->w_list.a[k])); k--);
    if(k < 0) return 0;
    if(!rev) {
        te = o->w_list.a[k].x_end+1;
    } else {
        ts = tl-o->w_list.a[k].x_end-1;
    }
    if(ts >= te) return 0;///[ts, te)
    

    for (i = 0; i < ct_n; i++) {
        z = &(ct_a[i]);
        if(z->e <= ts) continue; 
        if(z->s >= te) break;
        p = *z;
        if(!rev) {
            p.s = z->s; p.e = z->e;
        } else {
            p.s = tl - z->e; p.e = tl - z->s;
        }
        
        if(get_win_off(ql, &p, o, o_rate, &rr)) {
            rr.el = 1; kv_push(ul_ov_t, *res, rr);
        }
    }
    return res->n-rn;
}


 /**
void cal_simi_ul_ov_t(overlap_region *z, ul_ov_t *o, kv_ul_ov_t *res, int64_t ql)
{
    int64_t beg_q[2], beg_t[2], end_q[2], end_t[2], wk[2], wq[2], wt[2], kbeg, kend; 
    int64_t k, tot_l, qs, qe, ts, te, os, oe, ovlp, aln_l, wts, wte, wtl, wn; ul_ov_t *p;
    qs = 0; qe = ql-1; ts = o->ts; te = o->te-1; //[qs, qe] && [ps, pe]
    k = (o->qs==(uint32_t)-1)?(-1):(o->qs);
    k = get_win_yoff(z, ts, k, qs, qe, ts, te, beg_q, beg_t, &(wq[0]), &(wt[0])); wk[0] = k; 
    k = o->qe;
    k = get_win_yoff(z, te, k, qs, qe, ts, te, end_q, end_t, &(wq[1]), &(wt[1])); wk[1] = k;
    o->sec = 0;
    kbeg = ((wk[0]>=0)?(wk[0]>>1):(0)); kend = ((wk[1]>=0)?(wk[1]>>1):(0)); aln_l = 0;
    for (k = kbeg; k <= kend; k++) {
        wts = z->w_list.a[k].y_start; wte = z->w_list.a[k].y_end; wtl = wte + 1 - wts;
        os = MAX(wts, ts); oe = MIN(wte, te) + 1;
        ovlp = ((oe>os)? (oe-os):0); aln_l += ovlp;
        assert(ovlp > 0);
        if(ovlp == wtl) {

        } else {
            
        }
    }
    if(aln_l < te+1-ts) {
        wn = z->w_list.n;
        assert(wk[0] == -1 || wk[1] == wn-1);
    }

    if(wk[0] == wk[1]) {///one window cover the whole [ts, te]
        kv_pushp(ul_ov_t, *res, &p); p->ts = ts; p->te = te;
    }
    if(wk[0] < 0) {
        
    }
}
**/

///[ts, te) ->  this is the reverse coordinates of t, not the original coordinates of t
int64_t extract_subov(int64_t ts, int64_t te, overlap_region *o, double o_rate, int64_t *in_k, ul_ov_t *res)
{
    int64_t rev = o->y_pos_strand, qs, qe, k = 0, wn = o->w_list.n, ws, we, os, oe, ovlp, salnl;
    if(ts < ((int64_t)o->y_pos_s)) ts = o->y_pos_s;
    if(te > ((int64_t)o->y_pos_e+1)) te = o->y_pos_e+1;
    if(ts >= te) return 0;
    qs = get_chain_x_by_y(o, ts); 
    qe = get_chain_x_by_y(o, te-1) + 1;
    assert(qs < qe);
    memset(res, 0, sizeof(*res)); res->rev = rev;
    res->qs = qs; res->qe = qe; res->ts = ts; res->te = te;

    if(o_rate >= 0) {
        if(wn <= 0) return 0;
        if(in_k) k = *in_k;

        if(k < 0) k = 0; if(k >= wn) k = wn-1;
        for(; k < wn && qs > o->w_list.a[k].x_end; k++);
        if(k < 0) k = 0; if(k >= wn) k = wn-1;
        for(; k >= 0 && qs < o->w_list.a[k].x_start; k--);
        ///qs <= o->w_list.a[k].x_end && qs >= o->w_list.a[k].x_start
        if(k < 0) k = 0;
        if(in_k) *in_k = k;

        for (salnl = 0; k < wn; k++) {
            ws = o->w_list.a[k].x_start;
            we = o->w_list.a[k].x_end + 1;
            if(ws >= qe) break;
            if((o->w_list.a[k].y_end == -1) || (is_ualn_win(o->w_list.a[k]))) continue;
            os = MAX(qs, ws); oe = MIN(qe, we);
            ovlp = ((oe>os)? (oe-os):0);
            if(!ovlp) continue;
            salnl += ovlp; 
        }
        // if(o->y_id == 700) {
        //     fprintf(stderr, "[M::%s] raw_t::[%ld, %ld), raw_q::[%ld, %ld), salnl::%ld\n", __func__, ts, te, qs, qe, salnl);
        // }
        if((((qe-qs)*o_rate)<=salnl) && (salnl > 0)) return 1;
        else return 0;
    }
	return 1;
}


uint64_t gen_sub_ov_adv(const ul_idx_t *udb, overlap_region* o, double o_rate, utg_ct_t *ct_a, uint64_t ct_n, kv_ul_ov_t *res)
{
	uint64_t ts, te, i, l, rn = res->n, rev = o->y_pos_strand, s, e, rid, t[2]; ul_ov_t z;
	ma_utg_t *u = &(udb->ug->u.a[o->y_id]); int64_t k;
    if(!rev){
        ts = o->y_pos_s; te = o->y_pos_e + 1; k = 0;
    } else {
        ts = u->len - (o->y_pos_e+1); te = u->len - o->y_pos_s; 
        k = ((int64_t)o->w_list.n)-1; if(k < 0) k = 0;
    }

    if(!ct_a) {    
        for (i = l = 0; i < u->n; i++) {
            rid = u->a[i]>>33; 
            s = l; e = l + Get_READ_LENGTH(R_INF, rid);
            l += (uint32_t)u->a[i];
            if(e <= ts) continue; 
            if(s >= te) break;
            t[0] = (rev?(u->len-e):(s)); t[1] = (rev?(u->len-s):(e));
            if(extract_subov(t[0], t[1], o, o_rate, &k, &z)) {
                ///[ts, te) -> whole interval rid at the unitig adjusted by the reverse 
                // z.ts = t[0]; z.te = t[1]; 
                ///the strand of rid at the unitig
                z.rev = rev; 
                ///non-contained read at the unitg
                z.el = 0; 
                ///rid
                z.tn = rid;
                ///i-th read at the unitig
                z.qn = i;
                kv_push(ul_ov_t, *res, z);
            }
        }
    } else {
        for (i = 0; i < ct_n; i++) {
            rid = ct_a[i].x>>1; s = ct_a[i].s; e = ct_a[i].e;
            if(e <= ts) continue; 
            if(s >= te) break;
            t[0] = (rev?(u->len-e):(s)); t[1] = (rev?(u->len-s):(e));
            if(extract_subov(t[0], t[1], o, o_rate, &k, &z)) {
                ///[ts, te) -> whole interval rid at the unitig adjusted by the reverse 
                // z.ts = t[0]; z.te = t[1]; 
                ///the strand of rid at the unitig
                z.rev = rev; 
                ///contained read at the unitg
                z.el = 1; 
                ///rid
                z.tn = rid;
                ///i-th read at the unitig
                z.qn = i;
                kv_push(ul_ov_t, *res, z);
            }
        }
    }
    return res->n-rn;
}

int64_t return_t_chain(overlap_region *z, Candidates_list *cl)
{
    int64_t i, cn = cl->length, scn; uint64_t pid; k_mer_hit *ca;

    // if(z->x_id == 57 && z->y_id == 2175) {
    //     fprintf(stderr, "\n-0-[M::%s]\tutg%.6ul\tx::[%u,\t%u)\t%c\tutg%.6ul\ty::[%u,\t%u)\n", 
    //         __func__, z->x_id+1, z->x_pos_s, z->x_pos_e+1,
	// 		"+-"[z->y_pos_strand], z->y_id+1, z->y_pos_s, z->y_pos_e+1);
    //     i = z->shared_seed; pid = cl->list[i].readID;
    //     for (; i < cn && cl->list[i].readID == pid && cl->list[i].readID != ((uint32_t)(0x7fffffff)); i++) {
    //         fprintf(stderr, "i::%ld[M::%s]\treadID::%u\tself_offset::%u\toffset::%u\t%c\n", 
    //             i, __func__, cl->list[i].readID, cl->list[i].self_offset, cl->list[i].offset, 
    //             "+-"[cl->list[i].strand]);
    //     }
    // }



    i = z->shared_seed; pid = cl->list[i].readID;
    for (; i < cn && cl->list[i].readID == pid && cl->list[i].readID != ((uint32_t)(0x7fffffff)); i++);
    scn = i - z->shared_seed; ca = cl->list+z->shared_seed;
    i = lchain_refine(ca, scn, ca, &(cl->chainDP), 50, 5000, 512, 16); cn = i;
    for (; i < scn; i++) ca[i].readID = ((uint32_t)(0x7fffffff));



    // if(z->x_id == 57 && z->y_id == 2175) {
    //     fprintf(stderr, "\n-1-[M::%s]\tutg%.6ul\tx::[%u,\t%u)\t%c\tutg%.6ul\ty::[%u,\t%u)\n", 
    //         __func__, z->x_id+1, z->x_pos_s, z->x_pos_e+1,
	// 		"+-"[z->y_pos_strand], z->y_id+1, z->y_pos_s, z->y_pos_e+1);
    //     i = z->shared_seed; pid = cl->list[i].readID;
    //     for (; i < cn && cl->list[i].readID == pid && cl->list[i].readID != ((uint32_t)(0x7fffffff)); i++) {
    //         fprintf(stderr, "i::%ld[M::%s]\treadID::%u\tself_offset::%u\toffset::%u\t%c\n", 
    //             i, __func__, cl->list[i].readID, cl->list[i].self_offset, cl->list[i].offset, 
    //             "+-"[cl->list[i].strand]);
    //     }
    // }
    return cn;
}

// int64_t gen_mix_tchain(Candidates_list *cl, int64_t kidx, int64_t kn, ul_ov_t *oa, int64_t on)
// {
//     int64_t rcn = cl->length, kk, ok;
//     kv_resize_cl(k_mer_hit, *cl, rcn+kn);
//     kk = ok = 0;
//     while(kk < kn && ok < on) {
//         if(cl->list[kidx+kk].offset)
//     }
// }


uint64_t gen_woff_idx(overlap_region *z, asg64_v *oidx)
{
    uint64_t mm, k, aln = 0;
    kv_resize(uint64_t, *oidx, (z->w_list.n<<1)+2); 
    mm = z->x_pos_s; mm <<= 32; mm += z->y_pos_s; kv_push(uint64_t, *oidx, mm);
    for (k = oidx->n = 0; k < z->w_list.n; k++) {
        mm = z->w_list.a[k].x_start; mm <<= 32; mm += z->w_list.a[k].y_start;
        if((oidx->n == 0) && (mm != oidx->a[oidx->n-1])) {
            kv_push(uint64_t, *oidx, mm);
        }

        mm = z->w_list.a[k].x_end; mm <<= 32; mm += z->w_list.a[k].y_end;
        if((oidx->n == 0) && (mm != oidx->a[oidx->n-1])) {
            kv_push(uint64_t, *oidx, mm);
        }
        if(!(is_ualn_win(z->w_list.a[k]))) aln += z->w_list.a[k].x_end - z->w_list.a[k].x_start;
    }
    mm = z->x_pos_e; mm <<= 32; mm += z->y_pos_e; 
    if((oidx->n == 0) && (mm != oidx->a[oidx->n-1])) {
        kv_push(uint64_t, *oidx, mm);
    }
    return aln;
}

//return [rq, rt]
void win_boundary_offset(window_list *a, int64_t w_n, int64_t wi, int64_t toff, int64_t ql, int64_t *rq, int64_t *rt)
{
    int64_t q[2], t[2], qdis, tdis, dis; q[0] = q[1] = t[0] = t[1] = -1;
    ///[q[0], q[1]] && [t[0], t[1]]
    if(toff >= a[wi].y_start && toff <= a[wi].y_end) {///within the window
        q[0] = a[wi].x_start; q[1] = a[wi].x_end;
        t[0] = a[wi].y_start; t[1] = a[wi].y_end;
    } else if(toff < a[wi].y_start) {///before the window
        if(wi > 0) {
            q[0] = a[wi-1].x_end+1; q[1] = a[wi].x_start-1;
            t[0] = a[wi-1].y_end+1; t[1] = a[wi].y_start-1;
        } else {
            qdis = a[wi].x_start; 
            tdis = a[wi].y_start - toff;  
            dis = MIN(qdis, tdis);
            (*rq) = a[wi].x_start - dis;
            (*rt) = a[wi].y_start - dis; 
        }
    } else if(toff > a[wi].y_end) {//after the window
        if(wi < w_n - 1) {
            q[0] = a[wi].x_end+1; q[1] = a[wi+1].x_start-1;
            t[0] = a[wi].y_end+1; t[1] = a[wi+1].y_start-1;
        } else {
            qdis = ql-1-a[wi].x_end; 
            tdis = toff-a[wi].y_end;  
            dis = MIN(qdis, tdis);
            (*rq) = a[wi].x_end + dis;
            (*rt) = a[wi].y_end + dis;
        }
    }

    if(q[0] >= 0 && q[1] >= 0 && t[0] >= 0 && t[1] >= 0) {
        (*rt) = toff;
        (*rq) = q[0] + get_offset_adjust(toff-t[0], t[1]-t[0], q[1]-q[0]);
    }
}

int64_t hc_aln_exz_simi_adv(int64_t id, int64_t rev, const ul_idx_t *uref, hpc_t *hpc_g, 
All_reads *rref, char* qstr, UC_Read *tu, int64_t qs, int64_t qe, int64_t ts, int64_t te, 
int64_t qmin, int64_t qmax, int64_t tmin, int64_t tmax, int64_t mode, bit_extz_t *exz, int64_t q_tot, 
double e_rate, int64_t maxl, int64_t maxe, int64_t force_l, int64_t estimate_err, overlap_region *z,
int64_t gen_trace)
{
    clear_align(*exz); exz->thre = 0; ///mode cannot be 3
    int64_t thre, ql = qe - qs, thre0, pts = -1, pte = -1, pthre = -1, t_tot;
    if(estimate_err < 0) estimate_err = ql*e_rate;
    if(ql <= 0) return 0;
    if(hpc_g) t_tot = hpc_len(*hpc_g, id);
    else if(uref) t_tot = uref->ug->u.a[id].len;
    else t_tot = Get_READ_LENGTH((*rref), id);

    // if(ts == 18327 && te == 18601 && qs == 145990 && qe == 146191) {
    //     fprintf(stderr, "-0-[M::%s::] mode::%ld, q::[%ld, %ld), t::[%ld, %ld)\n", 
    //     __func__, mode, qs, qe, ts, te);
    // }

    if(ql <= 16) {
        if(cal_exact_simi_exz(uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, &pts, &pte, q_tot,
        t_tot, rev, id, mode, z)) {
            // ref_cigar_check(qstr, tu, uref, hpc_g, rref, z->y_id, z->y_pos_strand, exz);
            // fprintf(stderr, ", err::%d, thre::%d, scale::0(+)\n", exz->err, exz->thre);
            return 1;
        }
    }

    // if(ts == 18327 && te == 18601 && qs == 145990 && qe == 146191) {
    //     fprintf(stderr, "-a-[M::%s::] mode::%ld, q::[%ld, %ld), t::[%ld, %ld), ql::%ld, maxl::%ld, estimate_err::%ld, maxe::%ld\n", 
    //     __func__, mode, qs, qe, ts, te, ql, maxl, estimate_err, maxe);
    // }

    if(ql <= maxl && (estimate_err>>1) <= maxe) {
        thre = scale_ed_thre(estimate_err, maxe); 
        // if(ts == 18327 && te == 18601 && qs == 145990 && qe == 146191) {
        //     fprintf(stderr, "-1-[M::%s::] mode::%ld, q::[%ld, %ld), t::[%ld, %ld), thre::%ld\n", 
        //     __func__, mode, qs, qe, ts, te, thre);
        // }
        if(cal_exz_infi_simi_adv(uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, &pts, &pte, thre, &pthre, 
                qmin, qmax, tmin, tmax, rev, id, mode, z, gen_trace)) {
            // ref_cigar_check(qstr, tu, uref, hpc_g, rref, z->y_id, z->y_pos_strand, exz);
            // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(+)\n", exz->err, exz->thre, thre);
            // push_alnw(aux_o, exz);
            return 1;
        }

        thre0 = thre; thre = ql*e_rate; thre = scale_ed_thre(thre, maxe); 
        // if(ts == 18327 && te == 18601 && qs == 145990 && qe == 146191) {
        //     fprintf(stderr, "-2-[M::%s::] mode::%ld, q::[%ld, %ld), t::[%ld, %ld), thre::%ld\n", 
        //     __func__, mode, qs, qe, ts, te, thre);
        // }
        if(thre > thre0) {
            if(cal_exz_infi_simi_adv(uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, &pts, &pte, thre, &pthre, 
                qmin, qmax, tmin, tmax, rev, id, mode, z, gen_trace)) {
                // ref_cigar_check(qstr, tu, uref, hpc_g, rref, z->y_id, z->y_pos_strand, exz);
                // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(-)\n", exz->err, exz->thre, thre);
                // push_alnw(aux_o, exz);
                return 1;
            }
        }

        thre0 = thre; thre <<= 1; thre = scale_ed_thre(thre, maxe); 
        // if(ts == 18327 && te == 18601 && qs == 145990 && qe == 146191) {
        //     fprintf(stderr, "-3-[M::%s::] mode::%ld, q::[%ld, %ld), t::[%ld, %ld), thre::%ld\n", 
        //     __func__, mode, qs, qe, ts, te, thre);
        // }
        if(thre > thre0) {
            if(cal_exz_infi_simi_adv(uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, &pts, &pte, thre, &pthre, 
                qmin, qmax, tmin, tmax, rev, id, mode, z, gen_trace)) {
                // ref_cigar_check(qstr, tu, uref, hpc_g, rref, z->y_id, z->y_pos_strand, exz);
                // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(-)\n", exz->err, exz->thre, thre);
                // push_alnw(aux_o, exz);
                return 1;
            }
        }

        thre0 = thre; thre = ql*0.51; thre = scale_ed_thre(thre, maxe); 
        // if(ts == 18327 && te == 18601 && qs == 145990 && qe == 146191) {
        //     fprintf(stderr, "-4-[M::%s::] mode::%ld, q::[%ld, %ld), t::[%ld, %ld), thre::%ld\n", 
        //     __func__, mode, qs, qe, ts, te, thre);
        // }
        if(thre > thre0) {
            if(cal_exz_infi_simi_adv(uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, &pts, &pte, thre, &pthre, 
                qmin, qmax, tmin, tmax, rev, id, mode, z, gen_trace)) {
                // ref_cigar_check(qstr, tu, uref, hpc_g, rref, z->y_id, z->y_pos_strand, exz);
                // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(*)\n", exz->err, exz->thre, thre);
                // push_alnw(aux_o, exz);
                return 1;
            }
        }

        if(ql <= force_l) {
            thre = maxe; 
            // if(ts == 18327 && te == 18601 && qs == 145990 && qe == 146191) {
            //     fprintf(stderr, "-5-[M::%s::] mode::%ld, q::[%ld, %ld), t::[%ld, %ld), thre::%ld\n", 
            //     __func__, mode, qs, qe, ts, te, thre);
            // }
            if(cal_exz_infi_simi_adv(uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, &pts, &pte, thre, &pthre, 
                qmin, qmax, tmin, tmax, rev, id, mode, z, gen_trace)) {
                // ref_cigar_check(qstr, tu, uref, hpc_g, rref, z->y_id, z->y_pos_strand, exz);
                // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(*)\n", exz->err, exz->thre, thre);
                // push_alnw(aux_o, exz);
                return 1;
            }
        }
    }

    // if(ts == 18327 && te == 18601 && qs == 145990 && qe == 146191) {
    //     fprintf(stderr, "-b-[M::%s::] mode::%ld, q::[%ld, %ld), t::[%ld, %ld)\n", 
    //     __func__, mode, qs, qe, ts, te);
    // }
    // fprintf(stderr, ", err::%d, thre::%d\n", INT32_MAX, exz->thre);
    // if(mode == 0) {
    //     fprintf(stderr, "[M::%s::] pstr::%.*s\n", __func__, (int32_t)tu->length, tu->seq);
    //     fprintf(stderr, "[M::%s::] tstr::%.*s\n", __func__, (int32_t)(qe-qs), qstr+qs);
    // }
    return 0;
    
}

#define update_ul_ov_t_coor(z, qbeg, qend, tbeg, tend) do {\
        if(((int64_t)(z).qs) > (qbeg)) (z).qs = (qbeg);\
        if(((int64_t)(z).qe) < (qend)) (z).qe = (qend);\
        if(((int64_t)(z).ts) > (tbeg)) (z).ts = (tbeg);\
        if(((int64_t)(z).te) < (tend)) (z).te = (tend);\
	} while (0)

///[ps, pe) && [ts, te])
int64_t hc_aln_exz_by_exist_cigar_with_p(bit_extz_t *ez, int64_t mode, int64_t ps, int64_t pe, 
int64_t *ts, int64_t *te, int64_t *cis, int64_t *cie, int64_t *cps, int64_t *cpe, int64_t *cts, int64_t *cte)
{
    int64_t ts0 = (*ts), te0 = (*te), ts1, te1;
    assert(ez->ps <= ps && ez->pe+1 >= pe && mode > 0); 
    if(mode == 1) (*te) = -1;
    else if(mode == 2) (*ts) = -1;
    else (*ts) = (*te) = -1;
    (*cis) = (*cie) = (*cps) = (*cpe) = (*cts) = (*cte) = -1;

    int32_t pi = ez->ps, ti = ez->ts, pi0, ti0, err[2], tot_err, ws, we; 
    uint32_t ci = 0, cl; uint16_t c, sset = 0, eset = 0;
    err[0] = err[1] = tot_err = 0;
    while (ci < ez->cigar.n) {
        ci = pop_trace(&(ez->cigar), ci, &c, &cl);
        ws = pi; pi0 = pi; ti0 = ti;
        if(c <= 1) {
            pi+=cl; ti+=cl;
        } else if(c == 2) {///more p
            pi+=cl;
        } else if(c == 3) {
            ti+=cl;
        }
        we = pi; 
        // if((ps == 5139) && (pe == 21733)) {
        //     fprintf(stderr, "[M::%s::ci->%u] wp::[%d, %d), pi::%d, ti::%d, c::%u, cl::%u, tot_err::%d, p::[%ld, %ld), t::[%ld, %ld)\n", 
        //     __func__, ci, ws, we, pi, ti, c, cl, tot_err, ps, pe, (*ts), (*te));
        // }
        if(we < ps) {///not we <= ps
            if(c != 0) tot_err += cl;
            continue;
        }
        if(ws > pe) {///not ws >= pe
            if(c != 0) tot_err += cl;
            break;
        }

        if(!sset) {
            if((ps>=ws) && (ps<we)) {
                // if((ps == 5139) && (pe == 21733)) {
                //     fprintf(stderr, "+[M::%s::ci->%u] c::%u, cl::%u\n", __func__, ci, c, cl);
                // }
                if(c <= 1) {
                    ts1 = ti - (we - ps); 
                    if((ts1 == ts0) || (mode != 1)) {
                        if(c == 1) err[0] = tot_err + (ps - ws); 
                        else err[0] = tot_err;
                        sset = 1; (*ts) = ts1; 
                        (*cis) = ci-1; (*cps) = pi0; (*cts) = ti0;
                    }
                } else if(c == 2) {///more p
                    // ts1 = ti - (we - ps); 
                    ts1 = ti;
                    if((ts1 == ts0) || (mode != 1)) {
                        err[0] = tot_err + (ps - ws); 
                        sset = 1; (*ts) = ts1; 
                        (*cis) = ci-1; (*cps) = pi0; (*cts) = ti0;
                    }
                } else if(c == 3) {
                    ts1 = ti; 
                    if((ts1 == ts0) || (mode != 1)) {
                        err[0] = tot_err + cl; 
                        sset = 1; (*ts) = ts1; 
                        (*cis) = ci-1; (*cps) = pi0; (*cts) = ti0;
                    }
                }
            } 

            if((c == 3) && (ps>=ws) && (ps<=we)) {
                // if((ps == 5139) && (pe == 21733)) {
                //     fprintf(stderr, "-[M::%s::ci->%u] c::%u, cl::%u\n", __func__, ci, c, cl);
                // }
                ts1 = ti - cl; 
                if(ts1 == ts0) {
                    err[0] = tot_err; sset = 1; (*ts) = ts1; 
                    (*cis) = ci-1; (*cps) = pi0; (*cts) = ti0;
                }
            }
        }

        if(!eset) {
            if((pe>ws) && (pe<=we)) {
                if(c <= 1) {
                    te1 = ti - (we - pe); 
                    if((te1 == te0) || (mode != 2)) {
                        if(c == 1) err[1] = tot_err + (pe - ws);
                        else err[1] = tot_err;
                        eset = 1; (*te) = te1; 
                        (*cie) = ci; (*cpe) = pi; (*cte) = ti;
                    }
                } else if(c == 2) {///more p
                    // te1 = ti - (we - pe); 
                    te1 = ti;
                    if((te1 == te0) || (mode != 2)) {
                        err[1] = tot_err + (pe - ws); 
                        eset = 1; (*te) = te1; 
                        (*cie) = ci; (*cpe) = pi; (*cte) = ti;
                    }
                } else if(c == 3) {
                    te1 = ti-cl; 
                    if((te1 == te0) || (mode != 2)) {
                        err[1] = tot_err; 
                        eset = 1; (*te) = te1; 
                        (*cie) = ci; (*cpe) = pi; (*cte) = ti;
                    }
                }
            }
            if((c == 3) && (pe>=ws) && (pe<=we)) {
                te1 = ti; 
                if(te1 == te0) {
                    err[1] = tot_err + cl; 
                    eset = 1; (*te) = te1; 
                    (*cie) = ci; (*cpe) = pi; (*cte) = ti;
                }
            }
            if(eset) break;
        }

        if(c != 0) tot_err += cl;
    }

    // if(!(((mode == 1)&&((*ts) == ts0))||((mode == 2)&&((*te) == te0))||(mode == 3))) {
    //     fprintf(stderr, "+[M::%s::mode->%ld] t0::[%ld, %ld), t::[%ld, %ld), err::%d\n", 
    //     __func__, mode, ts0, te0, (*ts), (*te), err[1] - err[0]);
    // }
    assert(((mode == 1)&&((*ts) == ts0))||((mode == 2)&&((*te) == te0))||(mode == 3));
    return err[1] - err[0];
}


inline void update_trace_idx(rtrace_t *tc, int64_t wid, int64_t wid_s, int64_t wid_e, 
int64_t qs, int64_t qe, int64_t ts, int64_t te)
{
    if(qs < tc->c_qs || ts < tc->c_ts) {
        tc->c_qs = qs; 
        tc->c_ts = ts;
        tc->c_wsid = wid;
        tc->c_wsii = wid_s;
    }

    if(qe > tc->c_qe || te > tc->c_te) {
        tc->c_qe = qe; 
        tc->c_te = te;
        tc->c_weid = wid;
        tc->c_weii = wid_e;
    }
}

int64_t scan_single_wcigar(bit_extz_t *ez, int64_t csi, int64_t cei, int64_t cps, int64_t cpe, int64_t cts, int64_t cte, 
int64_t tar_ps, int64_t tar_pe, int64_t tar_ts, int64_t tar_te)
{
    int64_t ci = csi, pi = cps, ti = cts, ws, we;
    int64_t tot_err[2] = {0}, err[2] = {0}, ts1, te1; uint32_t cl; 
    uint16_t c, sset = 0, eset = 0;
    while (ci < cei) {
        ci = pop_trace(&(ez->cigar), ci, &c, &cl);
        ws = pi; 
        if(c <= 1) {
            pi+=cl; ti+=cl;
        } else if(c == 2) {///more p
            pi+=cl;
        } else if(c == 3) {
            ti+=cl;
        }
        we = pi;
        // if(cps == 4489 && cpe == 16192) {
        //     fprintf(stderr, "[M::%s::ci->%ld] wp::[%ld, %ld), pi::%ld, ti::%ld, c::%u, cl::%u, tot_err::%ld\n", 
        //     __func__, ci, ws, we, pi, ti, c, cl, tot_err[0]); 
        // }
        if(we < tar_ps) {///not we <= ps
            if(c != 0) tot_err[0] += cl;
            continue;
        }
        if(ws > tar_pe) {///not ws >= pe
            if(c != 0) tot_err[0] += cl;
            break;
        }

        if(!sset) {
            if((tar_ps>=ws) && (tar_ps<we)) {
                if(c <= 1) {
                    ts1 = ti - (we - tar_ps); 
                    if(ts1 == tar_ts) {
                        if(c == 1) err[0] = tot_err[0] + (tar_ps - ws); 
                        else err[0] = tot_err[0];
                        sset = 1; 
                    }
                } else if(c == 2) {///more p
                    // ts1 = ti - (we - tar_ps); 
                    ts1 = ti;
                    if(ts1 == tar_ts) {
                        err[0] = tot_err[0] + (tar_ps - ws); 
                        sset = 1; 
                    }
                } else if(c == 3) {
                    ts1 = ti; 
                    if(ts1 == tar_ts) {
                        err[0] = tot_err[0] + cl; 
                        sset = 1; 
                    }
                }
            } 

            if((c == 3) && (tar_ps>=ws) && (tar_ps<=we)) {
                ts1 = ti - cl; 
                if(ts1 == tar_ts) {
                    err[0] = tot_err[0]; sset = 1; 
                }
            }
        }

        if(!eset) {
            if((tar_pe>ws) && (tar_pe<=we)) {
                if(c <= 1) {
                    te1 = ti - (we - tar_pe); 
                    if(te1 == tar_te) {
                        if(c == 1) err[1] = tot_err[0] + (tar_pe - ws);
                        else err[1] = tot_err[0];
                        eset = 1; 
                    }
                } else if(c == 2) {///more p
                    // te1 = ti - (we - tar_pe); 
                    te1 = ti;
                    if(te1 == tar_te) {
                        err[1] = tot_err[0] + (tar_pe - ws); 
                        eset = 1; 
                    }
                } else if(c == 3) {
                    te1 = ti-cl; 
                    if(te1 == tar_te) {
                        err[1] = tot_err[0]; 
                        eset = 1; 
                    }
                }
            }
            if((c == 3) && (tar_pe>=ws) && (tar_pe<=we)) {
                te1 = ti; 
                if(te1 == tar_te) {
                    err[1] = tot_err[0] + cl; 
                    eset = 1; 
                }
            }
            if(eset) break;
        }

        if(c != 0) tot_err[0] += cl;
        // fprintf(stderr, "[M::%s::ci->%ld] tot_err::%ld, c::%u\n", __func__, ci, tot_err[0], c); 
    }
    return err[1] - err[0];
}


int64_t scan_single_wcigar_toff_backward_backup(bit_extz_t *ez, int64_t csi, int64_t cei, int64_t cps, int64_t cpe, int64_t cts, int64_t cte, 
int64_t tar_ps, int64_t tar_pe, int64_t tar_ts, int64_t tar_te)
{
    int64_t ci = cei-1, pi = cpe, ti = cte, wts, wte;
    int64_t tot_err = 0, err[2] = {0}, ps1, pe1; 
    uint16_t c, sset = 0, eset = 0; uint32_t cl; 
    while (ci >= csi) {
        ci = pop_trace_back(&(ez->cigar), ci, &c, &cl);
        wte = ti; 
        if(c <= 1) {
            pi-=cl; ti-=cl;
        } else if(c == 2) {///more p
            pi-=cl;
        } else if(c == 3) {
            ti-=cl;
        }
        wts = ti;
        // if(tar_ts == 0 && tar_te == 14233 && cts == 0 && cte == 14233) {
        //     fprintf(stderr, "[M::%s::ci->%ld] wt::[%ld, %ld), pi::%ld, ti::%ld, c::%u, cl::%u, tot_err::%ld, err[0]::%ld, err[1]::%ld\n", 
        //     __func__, ci, wts, wte, pi, ti, c, cl, tot_err[0], err[0], err[1]); 
        // }
        if(wte < tar_ts) {///not we <= ps
            if(c != 0) tot_err += cl;
            break;
        }
        if(wts > tar_te) {///not ws >= pe
            if(c != 0) tot_err += cl;
            continue;
        }

        if(!sset) {
            // if(tar_ts == 19030 && tar_te == 31898 && cts == 19030 && cte == 31983) {
            //     fprintf(stderr, "+[M::%s::ci->%ld] wt::[%ld, %ld), tar_t::[%ld, %ld)\n", 
            //     __func__, ci, wts, wte, tar_ts, tar_te); 
            // }
            if((tar_ts>=wts) && (tar_ts<wte)) {
                assert(c != 2);
                if(c <= 1) {
                    ps1 = pi + (tar_ts - wts); 
                    if(ps1 == tar_ps) {
                        if(c == 1) err[0] = tot_err + (wte - tar_ts); 
                        else err[0] = tot_err;
                        sset = 1; 
                    }
                } else if(c == 3) {///more t
                    ps1 = pi; 
                    if(ps1 == tar_ps) {
                        err[0] = tot_err + (wte - tar_ts);
                        sset = 1; 
                    }
                }
            } 

            if((c == 2) && (tar_ts>=wts) && (tar_ts<=wte)) {///more p
                // ps1 = pi; 
                // if(ps1 == tar_ps) {
                //     err[0] = tot_err[0] + cl; sset = 1; 
                // }
                if(tar_ps >= pi && tar_ps < pi + cl) {
                    err[0] = tot_err + (pi + cl - tar_ps);
                    sset = 1; 
                }
            }
            if(sset) break;
        }

        if(!eset) {
            if((tar_te>wts) && (tar_te<=wte)) {
                assert(c != 2);
                if(c <= 1) {
                    pe1 = pi + (tar_te - wts); 
                    // if(tar_ts == 19030 && tar_te == 31898 && cts == 19030 && cte == 31983) {
                    //     fprintf(stderr, "-[M::%s::ci->%ld] wt::[%ld, %ld), tar_t::[%ld, %ld), pe1::%ld, tar_pe::%ld\n", 
                    //     __func__, ci, wts, wte, tar_ts, tar_te, pe1, tar_pe); 
                    // }
                    if(pe1 == tar_pe) {
                        if(c == 1) err[1] = tot_err + (wte - tar_te);
                        else err[1] = tot_err;
                        eset = 1; 
                    }
                } else if(c == 3) {///more t
                    pe1 = pi;
                    if(pe1 == tar_pe) {
                        err[1] = tot_err + (wte - tar_te);
                        eset = 1; 
                    }
                }
            }

            if((c == 2) && (tar_te>=wts) && (tar_te<=wte)) {///more p
                // pe1 = pi; 
                // if(pe1 == tar_pe) {
                //     err[1] = tot_err[0] + cl; eset = 1; 
                // }
                if(tar_pe >= pi && tar_pe < pi + cl) {
                    err[1] = tot_err + (pi + cl - tar_pe);
                    eset = 1; 
                }
            }
        }

        if(c != 0) tot_err += cl;
        // fprintf(stderr, "[M::%s::ci->%ld] tot_err::%ld, c::%u\n", __func__, ci, tot_err[0], c); 
    }
    return err[0] - err[1];
}


int64_t scan_single_wcigar_toff_backward(bit_extz_t *ez, int64_t csi, int64_t cei, int64_t cps, int64_t cpe, int64_t cts, int64_t cte, 
int64_t tar_ps, int64_t tar_pe, int64_t tar_ts, int64_t tar_te)
{
    int64_t ci = cei-1, pi = cpe, ti = cte, wts, wte;
    int64_t e = 0, ps1, pe1; 
    uint16_t c, sset = 0, eset = 0, ff[2]; uint32_t cl; 
    while (ci >= csi) {
        ci = pop_trace_back(&(ez->cigar), ci, &c, &cl);
        wte = ti; 
        if(c <= 1) {
            pi-=cl; ti-=cl;
        } else if(c == 2) {///more p
            pi-=cl;
        } else if(c == 3) {
            ti-=cl;
        }
        wts = ti; ff[0] = ff[1] = 0;
        // if(tar_ts == 0 && tar_te == 14233 && cts == 0 && cte == 14233) {
        //     fprintf(stderr, "[M::%s::ci->%ld] wt::[%ld, %ld), pi::%ld, ti::%ld, c::%u, cl::%u, tot_err::%ld, err[0]::%ld, err[1]::%ld\n", 
        //     __func__, ci, wts, wte, pi, ti, c, cl, tot_err[0], err[0], err[1]); 
        // }
        if(wte < tar_ts) {///not we <= ps
            if(c != 0) e += cl;
            break;
        }
        if(wts > tar_te) {///not ws >= pe
            if(c != 0) e += cl;
            continue;
        }

        if(!eset) {
            if((tar_te>wts) && (tar_te<=wte)) {
                assert(c != 2);
                if(c <= 1) {
                    pe1 = pi + (tar_te - wts); 
                    if(pe1 == tar_pe) {
                        if(c == 1) e = tar_te - wts;
                        else e = 0;
                        eset = 1; 
                    }
                } else if(c == 3) {///more t
                    pe1 = pi;
                    if(pe1 == tar_pe) {
                        e = tar_te - wts;
                        eset = 1; 
                    }
                }
            }

            if((c == 2) && (tar_te==wts) && (tar_te==wte)) {///more p
                if(tar_pe >= pi && tar_pe < pi + cl) {
                    e = tar_pe - pi;
                    eset = 1; 
                }
            }
            // if(eset) continue;
            if(eset) ff[0] = 1;
        }


        if(!sset) {
            // if(tar_ts == 19030 && tar_te == 31898 && cts == 19030 && cte == 31983) {
            //     fprintf(stderr, "+[M::%s::ci->%ld] wt::[%ld, %ld), tar_t::[%ld, %ld)\n", 
            //     __func__, ci, wts, wte, tar_ts, tar_te); 
            // }
            if((tar_ts>=wts) && (tar_ts<wte)) {
                assert(c != 2);
                if(c <= 1) {
                    ps1 = pi + (tar_ts - wts); 
                    if(ps1 == tar_ps) {
                        if(c == 1) e += (wte - tar_ts); 
                        sset = 1; 
                    }
                } else if(c == 3) {///more t
                    ps1 = pi; 
                    if(ps1 == tar_ps) {
                        e += (wte - tar_ts);
                        sset = 1; 
                    }
                }
            } 

            if((c == 2) && (tar_ts==wts) && (tar_ts==wte)) {///more p
                if(tar_ps >= pi && tar_ps < pi + cl) {
                    e += (pi + cl - tar_ps);
                    sset = 1; 
                }
            }
            // if(sset) break;
            if(sset) ff[1] = 1;
        }

        if(ff[1]) break;
        if(ff[0]) continue;
        if(c!=0) e += cl;
        // fprintf(stderr, "[M::%s::ci->%ld] tot_err::%ld, c::%u\n", __func__, ci, tot_err[0], c); 
    }
    return e;
}


int64_t debug_aln_err(overlap_region *o, int64_t qs, int64_t qe, int64_t ts, int64_t te, double e_rate, 
int64_t estz_err, const ul_idx_t *uref, char* qstr, UC_Read *tu, bit_extz_t *exz)
{
    int64_t ql = qe - qs, tl = te - ts, id = o->y_id, rev = o->y_pos_strand;
    if((!ql) && (!tl)) return 0;
    if((!ql) && (tl)) return tl;
    if((ql) && (!tl)) return ql;
    if(estz_err < 0) estz_err = -1;
    if(hc_aln_exz_simi_adv(id, rev, uref, NULL, NULL, qstr, tu, qs, qe, ts, te, 
                            qs, qe, ts, te, 0, exz, ql, e_rate, MAX_CNS_L, MAX_CNS_E, FORCE_CNS_L, estz_err, NULL, 0)) {
        return exz->err;
    } else {
        return gen_err_unaligned(ql, tl); 
    }
}

int64_t get_sub_cigar_err(ul_ov_t *aln, rtrace_t *tc, const ul_idx_t *uref, char* qstr, UC_Read *tu, overlap_region *o, bit_extz_t *exz, double e_rate, int64_t is_rev)
{
    int64_t k, q[2], t[2], cq[2], ct[2], ci[2], err, qwl, twl; 
    window_list *wa = o->w_list.a; bit_extz_t aux;
    if(!is_rev) {
        err = 0; q[0] = aln->qs; t[0] = aln->ts;
        for (k = tc->c_wsid; k <= tc->c_weid; k++) {
            q[1] = wa[k].x_start; t[1] = wa[k].y_start;
            if(q[0] <= q[1] && t[0] <= t[1]) {
                err += debug_aln_err(o, q[0], q[1], t[0], t[1], e_rate, k==tc->c_wsid?tc->pfx_e:-1, uref, qstr, tu, exz);
                // fprintf(stderr, "-0-[M::%s::k->%ld] q::[%ld, %ld), t::[%ld, %ld), err::%ld\n", 
                //         __func__, k, q[0], q[1], t[0], t[1], err);
            }

            cq[0] = q[0] = wa[k].x_start; cq[1] = q[1] = wa[k].x_end+1; 
            ct[0] = t[0] = wa[k].y_start; ct[1] = t[1] = wa[k].y_end+1; 
            ci[0] = 0; ci[1] = wa[k].clen; qwl = q[1] - q[0]; twl = t[1] - t[0];
            if(is_ualn_win(wa[k])) {
                err += gen_err_unaligned(qwl, twl); 
                // fprintf(stderr, "-1-[M::%s::k->%ld] q::[%ld, %ld), t::[%ld, %ld), err::%ld\n", 
                //         __func__, k, q[0], q[1], t[0], t[1], err);
            } else {
                set_bit_extz_t(aux, (*o), k); 
                if(k == tc->c_wsid) {
                    cq[0] = tc->c_qs; ct[0] = tc->c_ts; ci[0] = tc->c_wsii;
                    q[0] = MAX((int32_t)aln->qs, tc->c_qs);
                    t[0] = MAX((int32_t)aln->ts, tc->c_ts);
                }
                if(k == tc->c_weid) {
                    cq[1] = tc->c_qe; ct[1] = tc->c_te; ci[1] = tc->c_weii;
                    q[1] = MIN((int32_t)aln->qe, tc->c_qe);
                    t[1] = MIN((int32_t)aln->te, tc->c_te);
                }
                // err += scan_single_wcigar(&aux, ci[0], ci[1], ct[0], ct[1], cq[0], cq[1], 
                //     t[0], t[1], q[0], q[1]);
                err += scan_single_wcigar_toff_backward(&aux, ci[0], ci[1], ct[0], ct[1], cq[0], cq[1], 
                    t[0], t[1], q[0], q[1]);
                // fprintf(stderr, "-2-[M::%s::k->%ld] q::[%ld, %ld), t::[%ld, %ld), err::%ld, cq::[%ld, %ld), ct::[%ld, %ld)\n", 
                //         __func__, k, q[0], q[1], t[0], t[1], err, cq[0], cq[1], ct[0], ct[1]);
            }
            q[0] = q[1]; t[0] = t[1];
        }

        q[1] = aln->qe; t[1] = aln->te;
        if(q[0] <= q[1] && t[0] <= t[1]) {
            err += debug_aln_err(o, q[0], q[1], t[0], t[1], e_rate, tc->sfx_e, uref, qstr, tu, exz);
            // fprintf(stderr, "-3-[M::%s::k->%ld] q::[%ld, %ld), t::[%ld, %ld), err::%ld\n", 
            //             __func__, k, q[0], q[1], t[0], t[1], err);
        }
    } else {
        err = 0; q[1] = aln->qe; t[1] = aln->te;
        for (k = tc->c_weid; k >= tc->c_wsid; k--) {
            q[0] = wa[k].x_end+1; t[0] = wa[k].y_end+1;
            if(q[0] <= q[1] && t[0] <= t[1]) {
                err += debug_aln_err(o, q[0], q[1], t[0], t[1], e_rate, k==tc->c_weid?tc->sfx_e:-1, uref, qstr, tu, exz);
            }

            cq[0] = q[0] = wa[k].x_start; cq[1] = q[1] = wa[k].x_end+1; 
            ct[0] = t[0] = wa[k].y_start; ct[1] = t[1] = wa[k].y_end+1; 
            ci[0] = 0; ci[1] = wa[k].clen; qwl = q[1] - q[0]; twl = t[1] - t[0];
            if(is_ualn_win(wa[k])) {
                err += gen_err_unaligned(qwl, twl); 
            } else {
                set_bit_extz_t(aux, (*o), k); 
                if(k == tc->c_wsid) {
                    cq[0] = tc->c_qs; ct[0] = tc->c_ts; ci[0] = tc->c_wsii;
                    q[0] = MAX((int32_t)aln->qs, tc->c_qs);
                    t[0] = MAX((int32_t)aln->ts, tc->c_ts);
                }
                if(k == tc->c_weid) {
                    cq[1] = tc->c_qe; ct[1] = tc->c_te; ci[1] = tc->c_weii;
                    q[1] = MIN((int32_t)aln->qe, tc->c_qe);
                    t[1] = MIN((int32_t)aln->te, tc->c_te);
                }
                // err += scan_single_wcigar(&aux, ci[0], ci[1], ct[0], ct[1], cq[0], cq[1], 
                //     t[0], t[1], q[0], q[1]);
                err += scan_single_wcigar_toff_backward(&aux, ci[0], ci[1], ct[0], ct[1], cq[0], cq[1], 
                    t[0], t[1], q[0], q[1]);
            }
            q[1] = q[0]; t[1] = t[0];
        }

        q[0] = aln->qs; t[0] = aln->ts;
        if(q[0] <= q[1] && t[0] <= t[1]) {
            err += debug_aln_err(o, q[0], q[1], t[0], t[1], e_rate, tc->pfx_e, uref, qstr, tu, exz);
        }
    }
    // fprintf(stderr, "[M::%s::] err::%ld, aln_err::%u\n", __func__, err, aln->sec);
    return err;
}


void gen_clip_win_err(overlap_region *o, int64_t qs, int64_t qe, int64_t ts, int64_t te, double e_rate, 
int64_t estz_err, const ul_idx_t *uref, char* qstr, UC_Read *tu, bit_extz_t *exz)
{
    int64_t ql = qe - qs, tl = te - ts, id = o->y_id, rev = o->y_pos_strand;
    exz->ts = qs; exz->te = qe-1; exz->ps = ts; exz->pe = te-1; 
    exz->cigar.n = 0; exz->thre = exz->err = INT32_MAX; 
    if((!ql) && (!tl)) {
        exz->err = 0; return;
    }
    if((!ql) && (tl)) {
        exz->err = tl; return;
    }
    if((ql) && (!tl)) {
        exz->err = ql; return;
    }
    if(estz_err < 0) estz_err = -1;
    if(hc_aln_exz_simi_adv(id, rev, uref, NULL, NULL, qstr, tu, qs, qe, ts, te, 
                            qs, qe, ts, te, 0, exz, ql, e_rate, MAX_CNS_L, MAX_CNS_E, FORCE_CNS_L, estz_err, NULL, 1)) {
        return;
    } else {
        exz->ts = qs; exz->te = qe-1; exz->ps = ts; exz->pe = te-1; 
        exz->cigar.n = 0; exz->thre = INT32_MAX;
        exz->err = gen_err_unaligned(ql, tl); 
        return;
    }
}


int64_t extract_sub_werr(bit_extz_t *ez, int64_t csi, int64_t cei, int64_t cps, int64_t cpe, int64_t cts, int64_t cte, 
int64_t tar_ps, int64_t tar_pe, int64_t tar_ts, int64_t tar_te, int64_t spec_toff, rtrace_iter *idx)
{
    if(idx->toff >= tar_pe && idx->qoff >= tar_te) {
        idx->coff = cei-1;
        idx->toff = cpe;
        idx->qoff = cte;
        idx->cerr = 0;
    }
    if(!(is_align(*ez))) return INT32_MAX;
    if(!(ez->cigar.n)) return (((double)(tar_te-spec_toff))/((double)(tar_te-tar_ts)))*ez->err;

    int64_t ci = idx->coff, pi = idx->toff, ti = idx->qoff, err = idx->cerr;
    int64_t t[2], ci0, pi0, ti0, err0, re = INT32_MAX, wts, wte, ps1, pe1;
    uint16_t c, sset = 0, eset = 0, ff[2]; uint32_t cl;
    while (ci >= csi) {
        ci0 = ci; pi0 = pi; ti0 = ti; err0 = err;
        ci = pop_trace_back(&(ez->cigar), ci, &c, &cl);
        wte = ti; 
        if(c <= 1) {
            pi-=cl; ti-=cl;
        } else if(c == 2) {///more p
            pi-=cl;
        } else if(c == 3) {
            ti-=cl;
        }
        wts = ti; ff[0] = ff[1] = 0;
        t[0] = wts; t[1] = wte;
        // if(spec_toff == 118374) {
        //     fprintf(stderr, "+[M::%s::ci->%ld::cl->%u::c->%u] spec_toff::%ld, t::[%ld, %ld), p::[%ld, %ld), err0::%ld\n", 
        //     __func__, ci, cl, c, spec_toff, ti, ti0, pi, pi0, err0);
        // }

        if(wte < tar_ts) {///not we <= ps
            if(c != 0) err += cl;
            break;
        }
        if(wts > tar_te) {///not ws >= pe
            if(c != 0) err += cl;
            continue;
        }

        if(!eset) {
            if((tar_te>wts) && (tar_te<=wte)) {
                assert(c != 2);
                if(c <= 1) {
                    pe1 = pi + (tar_te - wts); 
                    if(pe1 == tar_pe) {
                        if(c == 1) err = tar_te - wts;
                        else err = 0;
                        eset = 1; 
                    }
                } else if(c == 3) {///more t
                    pe1 = pi;
                    if(pe1 == tar_pe) {
                        err = tar_te - wts;
                        eset = 1; 
                    }
                }
            }

            if((c == 2) && (tar_te>=wts) && (tar_te<=wte)) {///more p
                if(tar_pe >= pi && tar_pe < pi + cl) {
                    err = tar_pe - pi;
                    eset = 1; 
                }
            }
            if(eset) {
                ff[0] = 1; t[1] = tar_te;
            }
        }


        if(!sset) {
            if((tar_ts>=wts) && (tar_ts<wte)) {
                assert(c != 2);
                if(c <= 1) {
                    ps1 = pi + (tar_ts - wts); 
                    if(ps1 == tar_ps) {
                        if(c == 1) err += (wte - tar_ts); 
                        sset = 1; 
                    }
                } else if(c == 3) {///more t
                    ps1 = pi; 
                    if(ps1 == tar_ps) {
                        err += (wte - tar_ts);
                        sset = 1; 
                    }
                }
            } 

            if((c == 2) && (tar_ts>=wts) && (tar_ts<=wte)) {///more p
                if(tar_ps >= pi && tar_ps < pi + cl) {
                    err += (pi + cl - tar_ps);
                    sset = 1; 
                }
            }
            if(sset) {
                ff[1] = 1; t[0] = tar_ts;
            }
        }

        // if(spec_toff == 118374) {
        //     fprintf(stderr, "-[M::%s::ci->%ld::cl->%u::c->%u] spec_toff::%ld, t::[%ld, %ld), p::[%ld, %ld), err0::%ld\n", 
        //     __func__, ci, cl, c, spec_toff, ti, ti0, pi, pi0, err0);
        // }
        if(spec_toff > t[1]) break;

        if((spec_toff >= t[0] && spec_toff < t[1]) || (spec_toff == t[0] && spec_toff == t[1])) {
            re = err + (((!ff[1]) && (!ff[0]) && (c!=0))?cl:0);
            if(c == 1 || c == 3) re -= (spec_toff - t[0]);
            idx->coff = ci0;
            idx->toff = pi0;
            idx->qoff = ti0;
            idx->cerr = err0;
            // return re;
        }
        if(ff[1]) break;
        if(ff[0]) continue;
        if(c!=0) err += cl;
    }
    return re;
}

//get error within [qe, ql)
int64_t get_rid_backward_cigar_err_back(rtrace_iter *it, ul_ov_t *aln, kv_rtrace_t *trace, rtrace_t *tc,
const ul_idx_t *uref, char* qstr, UC_Read *tu, overlap_region_alloc *ol, overlap_region *o, 
bit_extz_t *exz, double e_rate, int64_t qs)
{
    if(qs == aln->qe) return 0;
    if(!tc) tc = &(trace->a[aln->qn]);
    if(!o) o = &(ol->list[tc->oid]);

	if(it->k == INT32_MAX) {
		it->k = tc->c_weid;
		it->q[1] = aln->qe; 
		it->t[1] = aln->te;
		it->werr = it->cerr = 0; 
		it->qoff = aln->qe;
		it->toff = aln->te;
		it->coff = INT32_MAX; 
		clear_align(*exz); 
		exz->ps = exz->pe = exz->ts = exz->te = INT32_MAX;
	}
	window_list *wa = o->w_list.a; int64_t qwl, twl, sub_err; bit_extz_t aux; 
    assert(qs <= it->qoff);

	for (; it->k >= tc->c_wsid; it->k--) {
		it->q[0] = wa[it->k].x_end+1; 
		it->t[0] = wa[it->k].y_end+1;
        // if(qs == 166327) {
        //     fprintf(stderr, "-0-[M::%s::k->%ld] qoff::%ld, toff::%ld, coff::%ld, inner_q::[%ld, %ld), inner_t::[%ld, %ld), werr::%ld, qs::%ld\n", 
        //         __func__, it->k, it->qoff, it->toff, it->coff, it->q[0], it->q[1], it->t[0], it->t[1], it->werr, qs);
        // }
        
		if(it->qoff > it->q[0]) {///[it->qoff, ql) has been calculated
			if(it->q[0] <= it->q[1] && it->t[0] <= it->t[1]) {
				if((exz->ps != it->t[0]) || (exz->pe != it->t[1]) || (exz->ts != it->q[0]) || (exz->te == it->q[1])) {
					///calculate on-the-fly
                    gen_clip_win_err(o, it->q[0], it->q[1], it->t[0], it->t[1], e_rate, it->k==tc->c_weid?tc->sfx_e:-1, uref, qstr, tu, exz);
				}
                // fprintf(stderr, "-1-[M::%s::k->%ld] qoff::%ld, toff::%ld, coff::%ld, inner_q::[%ld, %ld), inner_t::[%ld, %ld), werr::%ld, qs::%ld\n", 
                // __func__, it->k, it->qoff, it->toff, it->coff, it->q[0], it->q[1], it->t[0], it->t[1], it->werr, qs);
				if(qs >= it->q[0] && qs < it->q[1]) {
					///phrase cigar
                    sub_err = extract_sub_werr(exz, 0, exz->cigar.n, exz->ps, exz->pe+1, exz->ts, exz->te+1, 
                            exz->ps, exz->pe+1, exz->ts, exz->te+1, qs, it);
                    if(qs == it->q[0]) {
                        it->werr += sub_err; sub_err = 0;
                        it->qoff = it->q[0]; it->toff = it->t[0]; 
                    }
					return it->werr+sub_err;
				} 
                it->qoff = it->q[0]; 
                it->toff = it->t[0]; 
                it->werr += exz->err;
			}
		} 
        if(qs == it->qoff) return it->werr;

		it->cq[0] = it->q[0] = wa[it->k].x_start; it->cq[1] = it->q[1] = wa[it->k].x_end+1; 
		it->ct[0] = it->t[0] = wa[it->k].y_start; it->ct[1] = it->t[1] = wa[it->k].y_end+1; 
		it->ci[0] = 0; it->ci[1] = wa[it->k].clen; qwl = it->q[1] - it->q[0]; twl = it->t[1] - it->t[0];
		if(is_ualn_win(wa[it->k])) {
			if(it->qoff > it->q[0]) {
                // fprintf(stderr, "-2-[M::%s::k->%ld] qoff::%ld, toff::%ld, coff::%ld, inner_q::[%ld, %ld), inner_t::[%ld, %ld), werr::%ld, qs::%ld\n", 
                // __func__, it->k, it->qoff, it->toff, it->coff, it->q[0], it->q[1], it->t[0], it->t[1], it->werr, qs);
				if(qs /**>=**/> it->q[0] && qs < it->q[1]) {
					sub_err = ((((double)(it->q[1]-qs))/((double)(it->q[1]-it->q[0])))*
                                                                    (gen_err_unaligned(qwl, twl))); 
                    // if(qs == it->q[0]) {
                    //     it->werr += sub_err; sub_err = 0;
                    //     it->qoff = it->q[0]; it->toff = it->t[0]; 
                    // }
					return it->werr+sub_err;
				}
				it->qoff = it->q[0]; 
                it->toff = it->t[0]; 
                it->werr += gen_err_unaligned(qwl, twl); 
			}
		} else {
			set_bit_extz_t(aux, (*o), it->k); 
			if(it->k == tc->c_wsid) {
				it->cq[0] = tc->c_qs; it->ct[0] = tc->c_ts; it->ci[0] = tc->c_wsii;
				it->q[0] = MAX((int32_t)aln->qs, tc->c_qs);
				it->t[0] = MAX((int32_t)aln->ts, tc->c_ts);
			}
			if(it->k == tc->c_weid) {
				it->cq[1] = tc->c_qe; it->ct[1] = tc->c_te; it->ci[1] = tc->c_weii;
				it->q[1] = MIN((int32_t)aln->qe, tc->c_qe);
				it->t[1] = MIN((int32_t)aln->te, tc->c_te);
			}
            // fprintf(stderr, "-3-[M::%s::k->%ld] qoff::%ld, toff::%ld, coff::%ld, inner_q::[%ld, %ld), inner_t::[%ld, %ld), werr::%ld, qs::%ld\n", 
            //     __func__, it->k, it->qoff, it->toff, it->coff, it->q[0], it->q[1], it->t[0], it->t[1], it->werr, qs);
			if(it->qoff > it->q[0]) {
				if(qs /**>=**/> it->q[0] && qs < it->q[1]) {
					sub_err = extract_sub_werr(&aux, it->ci[0], it->ci[1], it->ct[0], it->ct[1], 
                        it->cq[0], it->cq[1], it->t[0], it->t[1], it->q[0], it->q[1], qs, it);
                    // if(qs == it->q[0]) {
                    //     it->werr += sub_err; sub_err = 0;
                    //     it->qoff = it->q[0]; it->toff = it->t[0]; 
                    // }
					return it->werr + sub_err; 
				}
                it->werr += extract_sub_werr(&aux, it->ci[0], it->ci[1], it->ct[0], it->ct[1], 
                        it->cq[0], it->cq[1], it->t[0], it->t[1], it->q[0], it->q[1], it->q[0], it);
				it->qoff = it->q[0]; 
                it->toff = it->t[0]; 
			}
		}

        // if(qs == 166327) {
        //     fprintf(stderr, "-1-[M::%s::k->%ld] qoff::%ld, toff::%ld, coff::%ld, inner_q::[%ld, %ld), inner_t::[%ld, %ld), werr::%ld, qs::%ld\n", 
        //         __func__, it->k, it->qoff, it->toff, it->coff, it->q[0], it->q[1], it->t[0], it->t[1], it->werr, qs);
        // }
		it->q[1] = it->q[0]; it->t[1] = it->t[0];
        if((qs == it->qoff))  {
            if(it->k > tc->c_wsid) {
                it->q[0] = wa[it->k-1].x_end+1; 
		        it->t[0] = wa[it->k-1].y_end+1;
            } else {
                it->q[0] = aln->qs; 
                it->t[0] = aln->ts;
            }

            // if(qs == 166327) {
            //     fprintf(stderr, "-2-[M::%s::k->%ld] qoff::%ld, toff::%ld, coff::%ld, inner_q::[%ld, %ld), inner_t::[%ld, %ld), werr::%ld, qs::%ld\n", 
            //         __func__, it->k, it->qoff, it->toff, it->coff, it->q[0], it->q[1], it->t[0], it->t[1], it->werr, qs);
            //     fprintf(stderr, "-2-[M::%s::k->%ld] pre_inner_q::[%d, %d), pre_inner_t::[%d, %d)\n", 
            //         __func__, it->k, wa[it->k-1].x_start, wa[it->k-1].x_end+1, wa[it->k-1].y_start, wa[it->k-1].y_end+1);
            // }
            ///otherwise there should be an indel at t
            if(it->q[0] == it->q[1] && it->t[0] < it->t[1]) {
                it->werr += it->t[1] - it->t[0]; continue;
            } else {
                return it->werr;
            }
        }
    }

	it->q[0] = aln->qs; it->t[0] = aln->ts;
	if(it->qoff > it->q[0]) {
		if(it->q[0] <= it->q[1] && it->t[0] <= it->t[1]) {
			if((exz->ps != it->t[0]) || (exz->pe != it->t[1]) || (exz->ts != it->q[0]) || (exz->te == it->q[1])) {
				///calculate on-the-fly
                gen_clip_win_err(o, it->q[0], it->q[1], it->t[0], it->t[1], e_rate, tc->pfx_e, uref, qstr, tu, exz);
            }
            // fprintf(stderr, "-4-[M::%s::k->%ld] qoff::%ld, toff::%ld, coff::%ld, inner_q::[%ld, %ld), inner_t::[%ld, %ld), werr::%ld, qs::%ld\n", 
            //     __func__, it->k, it->qoff, it->toff, it->coff, it->q[0], it->q[1], it->t[0], it->t[1], it->werr, qs);
			if(qs >= it->q[0] && qs < it->q[1]) {
				///phrase cigar
                sub_err = extract_sub_werr(exz, 0, exz->cigar.n, exz->ps, exz->pe+1, exz->ts, exz->te+1, 
                            exz->ps, exz->pe+1, exz->ts, exz->te+1, qs, it);
                if(qs == it->q[0]) {
                    it->werr += sub_err; sub_err = 0;
                    it->qoff = it->q[0]; it->toff = it->t[0]; 
                }
				return it->werr+sub_err; 
			} 

            it->qoff = it->q[0]; 
            it->toff = it->t[0]; 
            it->werr += exz->err;
		}
	} 
    if(it->k < tc->c_wsid) it->werr = aln->sec;
	return it->werr;
}


//get error within [qe, ql)
int64_t get_rid_backward_cigar_err(rtrace_iter *it, ul_ov_t *aln, kv_rtrace_t *trace, rtrace_t *tc,
const ul_idx_t *uref, char* qstr, UC_Read *tu, overlap_region_alloc *ol, overlap_region *o, 
bit_extz_t *exz, double e_rate, int64_t qs)
{
    //this is not right
    // if(qs == aln->qe) return 0;
    if(!tc) tc = &(trace->a[aln->qn]);
    if(!o) o = &(ol->list[tc->oid]);
    if(it->k < tc->c_wsid && qs > ((int64_t)aln->qs)) it->k = INT32_MAX;

    if(it->k == INT32_MAX) {
        it->k = tc->c_weid;
        it->q[1] = MAX(tc->c_qe, ((int64_t)aln->qe)); 
        it->t[1] = MAX(tc->c_te, ((int64_t)aln->te));
        it->qoff = MAX(tc->c_qe, ((int64_t)aln->qe));
        it->toff = MAX(tc->c_te, ((int64_t)aln->te));
        it->cur_qoff = MAX(tc->c_qe, ((int64_t)aln->qe));

        it->werr = it->werr0 = it->cerr = 0; 
        it->f = 0;
        it->coff = INT32_MAX; 
        clear_align(*exz); 
        exz->ps = exz->pe = exz->ts = exz->te = INT32_MAX;
    }
    // if(qs == 45766) {
    //     fprintf(stderr, "\n[M::%s::] q::[%u, %u), t::[%u, %u), tot_e::%u, cid::[%d, %d], cq::[%u, %u), ct::[%u, %u), pfx_e::%d, sfx_e::%d, mid_e::%d\n", 
    //                     __func__, aln->qs, aln->qe, aln->ts, aln->te, aln->sec, tc->c_wsid, tc->c_weid,
    //                     tc->c_qs, tc->c_qe, tc->c_ts, tc->c_te, tc->pfx_e, tc->sfx_e, tc->mid_e);
    // }
    
    window_list *wa = o->w_list.a; int64_t qwl, twl, sub_err; bit_extz_t aux; 
    assert(qs <= it->qoff);
    if(qs == it->cur_qoff) return it->werr0;

    // if(qs == 166327) {
        // fprintf(stderr, "\n-*-[M::%s::k->%ld] qoff::%ld, toff::%ld, coff::%ld, inner_q::[%ld, %ld), inner_t::[%ld, %ld), werr::%ld, qs::%ld, cur_qoff::%ld, werr0::%ld\n", 
        //     __func__, it->k, it->qoff, it->toff, it->coff, it->q[0], it->q[1], it->t[0], it->t[1], it->werr, qs, 
        //     it->cur_qoff, it->werr0);
    // }

    for (; it->k >= tc->c_wsid; it->k--) {
        it->q[0] = wa[it->k].x_end+1; 
        it->t[0] = wa[it->k].y_end+1;
        // if(qs == 166327) {
            // fprintf(stderr, "-0-[M::%s::k->%ld] qoff::%ld, toff::%ld, coff::%ld, inner_q::[%ld, %ld), inner_t::[%ld, %ld), werr::%ld, qs::%ld, cur_qoff::%ld, werr0::%ld, f::%ld\n", 
            //     __func__, it->k, it->qoff, it->toff, it->coff, it->q[0], it->q[1], it->t[0], it->t[1], it->werr, qs, 
            //     it->cur_qoff, it->werr0, it->f);
        // }
        ///now qs >= it->qoff -> qs >= q[0]
        if(it->f == 0) {
            if(it->qoff >= it->q[0]) {///[it->qoff, ql) has been calculated
                if(it->q[0] <= it->q[1] && it->t[0] <= it->t[1]) {
                    if((exz->ps != it->t[0]) || (exz->pe != it->t[1]) || (exz->ts != it->q[0]) || (exz->te == it->q[1])) {
                        ///calculate on-the-fly
                        gen_clip_win_err(o, it->q[0], it->q[1], it->t[0], it->t[1], e_rate, it->k==tc->c_weid?tc->sfx_e:-1, uref, qstr, tu, exz);
                    }
                    sub_err = INT32_MAX;
                    ///if there are indels at either ends of t, qs == it->q[0] || qs == it->q[1]
                    if(qs >= it->q[0] && qs <= it->q[1]) {
                        ///phrase cigar
                        sub_err = extract_sub_werr(exz, 0, exz->cigar.n, exz->ps, exz->pe+1, exz->ts, exz->te+1, 
                                exz->ps, exz->pe+1, exz->ts, exz->te+1, qs, it);
                        if(sub_err != INT32_MAX) {///find the coordinate for qs
                            if((it->cur_qoff != qs) || ((it->werr+sub_err) > it->werr0)) {
                                it->werr0 = it->werr+sub_err;
                            }
                            it->cur_qoff = qs;
                            if(qs > it->q[0] && qs <= it->q[1]) return it->werr0;
                        } else {///happen when qs == it->q[1] and no indels at the right end
                            it->qoff = it->q[1]; 
                            it->toff = it->t[1]; 
                            // assert(it->cur_qoff == qs);
                            it->cur_qoff = qs;
                            return it->werr0;
                        } 
                    } 

                    if(qs <= it->q[0]) {
                        it->qoff = it->q[0]; 
                        it->toff = it->t[0]; 
                        it->werr += exz->err;
                        it->f = 1;
                        if((it->cur_qoff != it->q[0]) || (it->werr > it->werr0)) {
                            it->werr0 = it->werr;
                        }
                        it->cur_qoff = it->q[0];
                    }
                } else {
                    it->f = 1;
                }
            } else {
                it->f = 1;
            }
        }
        
        it->cq[0] = it->q[0] = wa[it->k].x_start; it->cq[1] = it->q[1] = wa[it->k].x_end+1; 
        it->ct[0] = it->t[0] = wa[it->k].y_start; it->ct[1] = it->t[1] = wa[it->k].y_end+1; 
        it->ci[0] = 0; it->ci[1] = wa[it->k].clen; qwl = it->q[1] - it->q[0]; twl = it->t[1] - it->t[0];
        // if(qs == 166327) {
            // fprintf(stderr, "-1-[M::%s::k->%ld] qoff::%ld, toff::%ld, coff::%ld, inner_q::[%ld, %ld), inner_t::[%ld, %ld), werr::%ld, qs::%ld, cur_qoff::%ld, werr0::%ld, f::%ld\n", 
            //     __func__, it->k, it->qoff, it->toff, it->coff, it->q[0], it->q[1], it->t[0], it->t[1], 
            //     it->werr, qs, it->cur_qoff, it->werr0, it->f);
        // }
        if(is_ualn_win(wa[it->k])) {
            if(it->qoff >= it->q[0] && it->f == 1) {
                ///ignore indels of t at both ends; so ignore qs == it->q[0] and qs == it->q[1]
                if(qs == it->q[1]) {
                    // assert(it->cur_qoff == qs);
                    it->cur_qoff = qs;
                    return it->werr0;
                }
                if(qs > it->q[0] && qs < it->q[1]) {
                    sub_err = ((((double)(it->q[1]-qs))/((double)(it->q[1]-it->q[0])))*
                                                                    (gen_err_unaligned(qwl, twl))); 
                    return it->werr+sub_err;
                }
                ///qs <= it->q[0]
                it->qoff = it->q[0]; 
                it->toff = it->t[0]; 
                it->werr += gen_err_unaligned(qwl, twl); 
                it->f = 0;
                if((it->cur_qoff != it->q[0]) || (it->werr > it->werr0)) {
                    it->werr0 = it->werr;
                }
                it->cur_qoff = it->q[0];
            }
        } else {
            set_bit_extz_t(aux, (*o), it->k); 
            if(it->k == tc->c_wsid) {
                it->cq[0] = tc->c_qs; it->ct[0] = tc->c_ts; it->ci[0] = tc->c_wsii;
                it->q[0] = MAX((int32_t)aln->qs, tc->c_qs);
                it->t[0] = MAX((int32_t)aln->ts, tc->c_ts);
            }
            if(it->k == tc->c_weid) {
                it->cq[1] = tc->c_qe; it->ct[1] = tc->c_te; it->ci[1] = tc->c_weii;
                it->q[1] = MIN((int32_t)aln->qe, tc->c_qe);
                it->t[1] = MIN((int32_t)aln->te, tc->c_te);
            }

            if(it->qoff >= it->q[0] && it->f == 1) {
                sub_err = INT32_MAX;
                ///if there are indels at either ends of t, qs == it->q[0] || qs == it->q[1]
                if(qs >= it->q[0] && qs <= it->q[1]) {
                    sub_err = extract_sub_werr(&aux, it->ci[0], it->ci[1], it->ct[0], it->ct[1], 
                        it->cq[0], it->cq[1], it->t[0], it->t[1], it->q[0], it->q[1], qs, it);
                    // if(qs == 166327) {
                        // fprintf(stderr, "-3-[M::%s::k->%ld] qoff::%ld, toff::%ld, coff::%ld, inner_q::[%ld, %ld), inner_t::[%ld, %ld), werr::%ld, qs::%ld, cur_qoff::%ld, werr0::%ld, f::%ld, sub_err::%ld\n", 
                        //     __func__, it->k, it->qoff, it->toff, it->coff, it->q[0], it->q[1], it->t[0], it->t[1], it->werr, qs, 
                        //     it->cur_qoff, it->werr0, it->f, sub_err);
                    // }
                    if(sub_err != INT32_MAX) {///find the coordinate for qs
                        if((it->cur_qoff != qs) || ((it->werr+sub_err) > it->werr0)) {
                            it->werr0 = it->werr+sub_err;
                        }
                        it->cur_qoff = qs;
                        if(qs > it->q[0] && qs <= it->q[1]) return it->werr0;
                    } else {///happen when qs == it->q[1] and no indels at the right end
                        it->qoff = it->q[1]; 
                        it->toff = it->t[1]; 
                        // if(!(it->cur_qoff == qs)) {
                        //     fprintf(stderr, "\n[M::%s::] q::[%u, %u), t::[%u, %u), tot_e::%u, cid::[%d, %d], cq::[%u, %u), ct::[%u, %u), pfx_e::%d, sfx_e::%d, mid_e::%d\n", 
                        //     __func__, aln->qs, aln->qe, aln->ts, aln->te, aln->sec, tc->c_wsid, tc->c_weid,
                        //     tc->c_qs, tc->c_qe, tc->c_ts, tc->c_te, tc->pfx_e, tc->sfx_e, tc->mid_e);
                        // }
                        // assert(it->cur_qoff == qs);
                        it->cur_qoff = qs;
                        return it->werr0;
                    } 
                }
                // if(qs == 166327) {
                    // fprintf(stderr, "-4-[M::%s::k->%ld] qoff::%ld, toff::%ld, coff::%ld, inner_q::[%ld, %ld), inner_t::[%ld, %ld), werr::%ld, qs::%ld, cur_qoff::%ld, werr0::%ld, f::%ld, sub_err::%ld\n", 
                    //     __func__, it->k, it->qoff, it->toff, it->coff, it->q[0], it->q[1], it->t[0], it->t[1], it->werr, qs, 
                    //     it->cur_qoff, it->werr0, it->f, sub_err);
                // }
                ///qs <= it->q[0]
                if(qs <= it->q[0]) {
                    it->werr += ((qs==it->q[0])?(sub_err):(extract_sub_werr(&aux, it->ci[0], it->ci[1], it->ct[0], it->ct[1], 
                                    it->cq[0], it->cq[1], it->t[0], it->t[1], it->q[0], it->q[1], it->q[0], it)));
                    it->qoff = it->q[0]; 
                    it->toff = it->t[0]; 
                    it->f = 0;
                    if((it->cur_qoff != it->q[0]) || (it->werr > it->werr0)) {
                        it->werr0 = it->werr;
                    }
                    it->cur_qoff = it->q[0];
                }
            }
        }
        it->q[1] = it->q[0]; it->t[1] = it->t[0];
        // if(qs == 166327) {
            // fprintf(stderr, "-#-[M::%s::k->%ld] qoff::%ld, toff::%ld, coff::%ld, inner_q::[%ld, %ld), inner_t::[%ld, %ld), werr::%ld, qs::%ld, cur_qoff::%ld, werr0::%ld, f::%ld\n", 
            //     __func__, it->k, it->qoff, it->toff, it->coff, it->q[0], it->q[1], it->t[0], it->t[1], it->werr, qs, 
            //     it->cur_qoff, it->werr0, it->f);
        // }
    }

    it->q[0] = aln->qs; it->t[0] = aln->ts;
    if(it->qoff >= it->q[0] && it->f == 0) {
        if(it->q[0] <= it->q[1] && it->t[0] <= it->t[1]) {
            if((exz->ps != it->t[0]) || (exz->pe != it->t[1]) || (exz->ts != it->q[0]) || (exz->te == it->q[1])) {
                ///calculate on-the-fly
                gen_clip_win_err(o, it->q[0], it->q[1], it->t[0], it->t[1], e_rate, tc->pfx_e, uref, qstr, tu, exz);
            }
            sub_err = INT32_MAX;
            ///if there are indels at either ends of t, qs == it->q[0] || qs == it->q[1]
            if(qs >= it->q[0] && qs <= it->q[1]) {
                ///phrase cigar
                sub_err = extract_sub_werr(exz, 0, exz->cigar.n, exz->ps, exz->pe+1, exz->ts, exz->te+1, 
                            exz->ps, exz->pe+1, exz->ts, exz->te+1, qs, it);
                if(sub_err != INT32_MAX) {///find the coordinate for qs
                    if((it->cur_qoff != qs) || ((it->werr+sub_err) > it->werr0)) {
                        it->werr0 = it->werr+sub_err;
                    }
                    it->cur_qoff = qs;
                    if(qs > it->q[0] && qs <= it->q[1]) return it->werr0;
                } else {///happen when qs == it->q[1] and no indels at the right end
                    it->qoff = it->q[1]; 
                    it->toff = it->t[1]; 
                    // assert(it->cur_qoff == qs);
                    it->cur_qoff = qs;
                    return it->werr0;
                } 
            } 

            if(qs <= it->q[0]) {
                it->qoff = it->q[0]; 
                it->toff = it->t[0]; 
                it->werr += exz->err;
                it->f = 1;
                if((it->cur_qoff != it->q[0]) || (it->werr > it->werr0)) {
                    it->werr0 = it->werr;
                }
                it->cur_qoff = it->q[0];
            }
        }
    } 
    // if(it->k < tc->c_wsid) it->werr = aln->sec;
    return it->werr;
}

void debug_backtrace_step_err(uint64_t rid, ul_ov_t *aln, rtrace_t *tc, const ul_idx_t *uref, char* qstr, UC_Read *tu, overlap_region *o, bit_extz_t *exz, double e_rate)
{
    // if(aln->qs == 0 && aln->qe == 3013 && aln->ts == 28965 && aln->te == 31993) {
        // fprintf(stderr, "\n[M::%s::] q::[%u, %u), t::[%u, %u), tot_e::%u, cid::[%d, %d], cq::[%u, %u), ct::[%u, %u), pfx_e::%d, sfx_e::%d, mid_e::%d\n", 
        //                 __func__, aln->qs, aln->qe, aln->ts, aln->te, aln->sec, tc->c_wsid, tc->c_weid,
        //                 tc->c_qs, tc->c_qe, tc->c_ts, tc->c_te, tc->pfx_e, tc->sfx_e, tc->mid_e);
        int64_t k, err, qs = aln->qs, qe = aln->qe, err0; rtrace_iter it;
        k = qe; it.k = INT32_MAX;
        err0 = get_rid_backward_cigar_err(&it, aln, NULL, tc, uref, qstr, tu, NULL, o, exz, e_rate, k); 
        err = get_rid_backward_cigar_err(&it, aln, NULL, tc, uref, qstr, tu, NULL, o, exz, e_rate, k); 
        assert(err == err0); 
        for (k = qe, it.k = INT32_MAX, err0 = 0; k >= qs; k-=8) {
            err = get_rid_backward_cigar_err(&it, aln, NULL, tc, uref, qstr, tu, NULL, o, exz, e_rate, k);  
            // if(!(err >= 0 && err >= err0)) {
            //     fprintf(stderr, "[M::%s::rid->%lu] q::[%u, %u), t::[%u, %u), tot_e::%u, cid::[%d, %d], cq::[%u, %u), ct::[%u, %u)\n", 
            //             __func__, rid, aln->qs, aln->qe, aln->ts, aln->te, aln->sec, tc->c_wsid, tc->c_weid,
            //             tc->c_qs, tc->c_qe, tc->c_ts, tc->c_te);
            //     fprintf(stderr, "[M::%s::] q::[%ld, %ld), err::%ld, err0::%ld\n", __func__, k, qe, err, err0);   
            // }
            // fprintf(stderr, "[M::%s::] q::[%ld, %ld), err::%ld\n", __func__, k, qe, err); 
            assert(err >= 0 && err >= err0); 
            err0 = err;  
        }
        k = qs;
        err = get_rid_backward_cigar_err(&it, aln, NULL, tc, uref, qstr, tu, NULL, o, exz, e_rate, k); 
        // if(!(err >= 0 && err >= err0)) {
        //     fprintf(stderr, "[M::%s::rid->%lu] q::[%u, %u), t::[%u, %u), tot_e::%u, cid::[%d, %d], cq::[%u, %u), ct::[%u, %u)\n", 
        //             __func__, rid, aln->qs, aln->qe, aln->ts, aln->te, aln->sec, tc->c_wsid, tc->c_weid,
        //             tc->c_qs, tc->c_qe, tc->c_ts, tc->c_te);
        //     fprintf(stderr, "[M::%s::] q::[%ld, %ld), err::%ld, err0::%ld\n", __func__, k, qe, err, err0);   
        // }
        // fprintf(stderr, "[M::%s::] q::[%ld, %ld), err::%ld, exz->err::%d, exz->cigar.n::%d\n", 
        // __func__, qs, qe, err, exz->err, (int32_t)exz->cigar.n);    
        assert(err >= 0 && err >= err0); 
        if(!(err == (int64_t)aln->sec)) {
            fprintf(stderr, "[M::%s::rid->%lu] q::[%u, %u), t::[%u, %u), tot_e::%u, cid::[%d, %d], cq::[%u, %u), ct::[%u, %u)\n", 
                        __func__, rid, aln->qs, aln->qe, aln->ts, aln->te, aln->sec, tc->c_wsid, tc->c_weid,
                        tc->c_qs, tc->c_qe, tc->c_ts, tc->c_te);
            fprintf(stderr, "[M::%s::] aln->sec::%u, err::%ld\n", __func__, aln->sec, err);  
        }
        assert(err == (int64_t)aln->sec);  
    // }   
}

///[wsid, weid) && [ts, te)
void gen_raln(const ul_idx_t *uref, char* qstr, UC_Read *tu, overlap_region *o, bit_extz_t *exz, 
int64_t wsid, int64_t weid, int64_t ts, int64_t te, int64_t ql, int64_t id, int64_t rev, double e_rate, 
uint64_t rid, rtrace_t *tc, ul_ov_t *res)
{
    // fprintf(stderr, "\n[M::%s::] ii::[%ld, %ld), t::[%ld, %ld), ql::%ld, id::%ld\n", 
    //         __func__, wsid, weid, ts, te, ql, id);
    int64_t k, q[2], t[2], c[2], ct[2], cq[2], mode, qwl, twl, aln_e = 0, cur_e, is_aln; 
    bit_extz_t aux; window_list *wa = o->w_list.a; 
    memset(res, 0, sizeof((*res))); 
    res->qs = res->ts = UINT32_MAX; res->qe = res->te = 0;
    memset(tc, 0, sizeof((*tc)));
    tc->c_qs = tc->c_ts = INT32_MAX; tc->c_qe = tc->c_te = -1;
    tc->pfx_e = tc->mid_e = tc->sfx_e = 0;

    q[0] = q[1] = -1; t[0] = ts; t[1] = te;///[q[0], q[1]) && [t[0], t[1])
    for (k = wsid; k < weid; k++) {
        q[1] = wa[k].x_start; t[1] = wa[k].y_start; 
        mode = - 1; cur_e = 0; is_aln = 1;
        // if(rid == 53 && id == 6) {
        //     fprintf(stderr, "0-a[M::%s::k->%ld] q::[%ld, %ld), t::[%ld, %ld)\n", __func__, k, q[0], q[1], t[0], t[1]);
        // }
        //before window; there are gaps before the window that need to be filled
        if(t[0] < t[1] && q[0] < q[1]) {
            if(q[0] < 0) { ///backward extension
                mode = 2;
                adjust_ext_offset_fixed_t(&(q[0]), &(q[1]), &(t[0]), &(t[1]), 0, q[1], t[0], t[1], 0, mode);
            } else {
                mode = 0;
            }
            // if(rid == 53 && id == 6) {
            //     fprintf(stderr, "1-a[M::%s::k->%ld] q::[%ld, %ld), t::[%ld, %ld), mode::%ld\n", 
            //     __func__, k, q[0], q[1], t[0], t[1], mode);
            // }
            qwl = q[1] - q[0]; twl = t[1] - t[0];
            if(qwl || twl) {
                if(qwl == 0 && twl > 0) {
                    cur_e = twl;
                    // tot_e += twl;
                } else if(twl == 0 && qwl > 0) {
                    cur_e = qwl;
                    // tot_e += qwl;
                } else {
                    if(hc_aln_exz_simi_adv(id, rev, uref, NULL, NULL, qstr, tu, q[0], q[1], t[0], t[1], 
                            0, ql, ts, te, mode, exz, ql, e_rate, MAX_CNS_L, MAX_CNS_E, FORCE_CNS_L, -1, NULL, 0)) {
                        cur_e = exz->err; 
                        // tot_e += exz->err; 
                        q[0] = exz->ts; q[1] = exz->te + 1; t[0] = exz->ps; t[1] = exz->pe + 1;
                        // if(rid == 53 && id == 6) {
                        //     fprintf(stderr, "2-a[M::%s::k->%ld] q::[%ld, %ld), t::[%ld, %ld), mode::%ld, cur_e::%ld\n", 
                        //     __func__, k, q[0], q[1], t[0], t[1], mode, cur_e);
                        // }
                    } else {
                        cur_e = gen_err_unaligned(qwl, twl); is_aln = 0;
                        // tot_e += gen_err_unaligned(qwl, twl); 
                        // if(rid == 53 && id == 6) {
                        //     fprintf(stderr, "2-b[M::%s::k->%ld] q::[%ld, %ld), t::[%ld, %ld), mode::%ld, cur_e::%ld\n", 
                        //     __func__, k, q[0], q[1], t[0], t[1], mode, cur_e);
                        // }
                    }
                }
                aln_e += cur_e;
                update_ul_ov_t_coor((*res), q[0], q[1], t[0], t[1]);
                if(k == wsid) tc->pfx_e = cur_e*(is_aln?1:-1);
            }
        } 
        
        ///within window
        q[0] = wa[k].x_start; t[0] = wa[k].y_start;
        q[1] = wa[k].x_end+1; t[1] = wa[k].y_end+1;
        qwl = q[1] - q[0]; twl = t[1] - t[0]; cur_e = 0; is_aln = 1;
        // if(rid == 53 && id == 6) {
        //     fprintf(stderr, "3-a[M::%s::k->%ld] q::[%ld, %ld), t::[%ld, %ld)\n", 
        //     __func__, k, q[0], q[1], t[0], t[1]);
        // }
        if((ts <= t[0]) && (te >= t[1])) {///cover the whole window
            if(!is_ualn_win(wa[k])) cur_e = wa[k].error;
            else cur_e = gen_err_unaligned(qwl, twl); 
            update_trace_idx(tc, k, 0, wa[k].clen, q[0], q[1], t[0], t[1]);
            // if(rid == 53 && id == 6) {
            //     fprintf(stderr, "4-a[M::%s::k->%ld] q::[%ld, %ld), t::[%ld, %ld), cur_e::%ld\n", 
            //     __func__, k, q[0], q[1], t[0], t[1], cur_e);
            // }
        } else {///te < t[1]->cover a part of window
            if(ts > t[0]) {
                t[0] = ts; q[0] = -1;
            }
            if(te < t[1]) {
                t[1] = te; q[1] = -1;
            }

            if((q[0] != -1) && (q[1] != -1)) {
                mode = 0;//global
            } else if((q[0] != -1) && (q[1] == -1)) {
                mode = 1;///forward extension
            } else if((q[0] == -1) && (q[1] != -1)) {
                mode = 2;///backward extension
            } else {
                mode = 3;///no primary hit within [ibeg, iend] 
            }
            // if(rid == 53 && id == 6) {
            //     fprintf(stderr, "5-a[M::%s::k->%ld] q::[%ld, %ld), t::[%ld, %ld), mode::%ld\n", 
            //     __func__, k, q[0], q[1], t[0], t[1], mode);
            // }
            if(!is_ualn_win(wa[k])) {///scan cigar by the coordinates of y/t
                set_bit_extz_t(aux, (*o), k); 
                cur_e = hc_aln_exz_by_exist_cigar_with_p(&aux, mode, t[0], t[1], &(q[0]), &(q[1]), 
                &(c[0]), &(c[1]), &(ct[0]), &(ct[1]), &(cq[0]), &(cq[1]));
                update_trace_idx(tc, k, c[0], c[1], cq[0], cq[1], ct[0], ct[1]);
                // if(rid == 53 && id == 6) {
                //     fprintf(stderr, "6-a[M::%s::k->%ld] q::[%ld, %ld), t::[%ld, %ld), cur_e::%ld, c_q::[%ld, %ld), c_t::[%ld, %ld)\n", 
                //     __func__, k, q[0], q[1], t[0], t[1], cur_e, cq[0], cq[1], ct[0], ct[1]);
                // }
            } else {
                assert(mode == 1 || mode == 2);
                adjust_ext_offset_fixed_t(&(q[0]), &(q[1]), &(t[0]), &(t[1]), q[0], q[1], t[0], t[1], 0, mode);
                qwl = q[1] - q[0]; twl = t[1] - t[0];
                // if(rid == 53 && id == 6) {
                //     fprintf(stderr, "7-a[M::%s::k->%ld] q::[%ld, %ld), t::[%ld, %ld), cur_e::%ld\n", 
                //     __func__, k, q[0], q[1], t[0], t[1], cur_e);
                // }
                if(qwl || twl) {
                    if(qwl == 0 && twl > 0) {
                        cur_e = twl;
                    } else if(twl == 0 && qwl > 0) {
                        cur_e = qwl;
                    } else {
                        if(hc_aln_exz_simi_adv(id, rev, uref, NULL, NULL, qstr, tu, q[0], q[1], t[0], t[1], 
                                0, ql, ts, te, mode, exz, ql, e_rate, MAX_CNS_L, MAX_CNS_E, FORCE_CNS_L, -1, NULL, 0)) {
                            cur_e = exz->err; 
                            q[0] = exz->ts; q[1] = exz->te + 1; t[0] = exz->ps; t[1] = exz->pe + 1;
                            // if(rid == 53 && id == 6) {
                            //     fprintf(stderr, "7-b[M::%s::k->%ld] q::[%ld, %ld), t::[%ld, %ld), cur_e::%ld\n", 
                            //     __func__, k, q[0], q[1], t[0], t[1], cur_e);
                            // }
                        } else {
                            cur_e = gen_err_unaligned(qwl, twl); is_aln = 0;
                            // if(rid == 53 && id == 6) {
                            //     fprintf(stderr, "7-c[M::%s::k->%ld] q::[%ld, %ld), t::[%ld, %ld), cur_e::%ld\n", 
                            //     __func__, k, q[0], q[1], t[0], t[1], cur_e);
                            // }
                        }
                    }
                    tc->sfx_e = cur_e*(is_aln?1:-1);
                }
            }
        }
        aln_e += cur_e;
        update_ul_ov_t_coor((*res), q[0], q[1], t[0], t[1]);

        q[0] = q[1]; t[0] = t[1];
    }
    // if(ts == 8087 && te == 28530 && ql == 127662) {
    //     fprintf(stderr, "8-a[M::%s::k->%ld] q::[%ld, %ld), t::[%ld, %ld), cur_e::%ld, ql::%ld\n", 
    //     __func__, k, q[0], q[1], t[0], t[1], cur_e, ql);
    // }
    if(q[0] >= 0 && q[0] < ql) {///forward extension
        q[1] = ql; t[1] = te; cur_e = 0; mode = 1;///forward extension
        adjust_ext_offset_fixed_t(&(q[0]), &(q[1]), &(t[0]), &(t[1]), q[0], q[1], t[0], t[1], 0, mode);
        qwl = q[1] - q[0]; twl = t[1] - t[0]; is_aln = 1;
        // if(rid == 53 && id == 6) {
        //     fprintf(stderr, "9-a[M::%s::k->%ld] q::[%ld, %ld), t::[%ld, %ld), cur_e::%ld\n", 
        //     __func__, k, q[0], q[1], t[0], t[1], cur_e);
        // }
        if(qwl || twl) {
            if(qwl == 0 && twl > 0) {
                cur_e = twl;
            } else if(twl == 0 && qwl > 0) {
                cur_e = qwl;
            } else {
                if(hc_aln_exz_simi_adv(id, rev, uref, NULL, NULL, qstr, tu, q[0], q[1], t[0], t[1], 
                        0, ql, ts, te, mode, exz, ql, e_rate, MAX_CNS_L, MAX_CNS_E, FORCE_CNS_L, -1, NULL, 0)) {
                    cur_e = exz->err; 
                    q[0] = exz->ts; q[1] = exz->te + 1; t[0] = exz->ps; t[1] = exz->pe + 1;
                    // if(rid == 53 && id == 6) {
                    //     fprintf(stderr, "9-b[M::%s::k->%ld] q::[%ld, %ld), t::[%ld, %ld), cur_e::%ld, thre::%d\n", 
                    //     __func__, k, q[0], q[1], t[0], t[1], cur_e, exz->thre);
                    // }
                } else {
                    cur_e = gen_err_unaligned(qwl, twl); is_aln = 0;
                    // if(rid == 53 && id == 6) {
                    //     fprintf(stderr, "9-c[M::%s::k->%ld] q::[%ld, %ld), t::[%ld, %ld), cur_e::%ld\n", 
                    //     __func__, k, q[0], q[1], t[0], t[1], cur_e);
                    // }
                }
            }
            update_ul_ov_t_coor((*res), q[0], q[1], t[0], t[1]);
            tc->sfx_e = cur_e*(is_aln?1:-1);
        }
        aln_e += cur_e;
    }
    tc->mid_e = aln_e - abs(tc->pfx_e) - abs(tc->sfx_e); 
    if(aln_e <= MAX_SEC_ERR) res->sec = aln_e;
    else res->sec = MAX_SEC_ERR;
    // fprintf(stderr, "10-a[M::%s::k->%ld] q::[%u, %u), t::[%u, %u), tot_e::%u\n\n", 
    //                 __func__, k, res->qs, res->qe, res->ts, res->te, res->sec);


    // aln_e = get_sub_cigar_err(res, tc, uref, qstr, tu, o, exz, e_rate);
    // if(aln_e != (int64_t)res->sec) {
    //     fprintf(stderr, "[M::%s::rid->%lu] tot_e::%u, aln_e::%ld, id::%ld\n", 
    //     __func__, rid, res->sec, aln_e, id);
    //     exit(1);
    // }
    
    // assert(aln_e == get_sub_cigar_err(res, tc, uref, qstr, tu, o, exz, e_rate, 0));
    // assert(aln_e == get_sub_cigar_err(res, tc, uref, qstr, tu, o, exz, e_rate, 1));
    // debug_backtrace_step_err(rid, res, tc, uref, qstr, tu, o, exz, e_rate);
    // if(rid == 42 && res->sec == 706 && res->qs == 63686 && res->qe == 75773 
    //                     && res->ts == 9404 && res->te == 21663) {
    //     debug_backtrace_step_err(rid, res, tc, uref, qstr, tu, o, exz, e_rate);
    // }
}

void prt_aln_w(overlap_region *o)
{
    uint64_t k;
    for (k = 0; k < o->w_list.n; k++) {
        fprintf(stderr, "[M::%s::aln->%u] q::[%d, %d), t::[%d, %d), err::%d\n", 
            __func__, ((o->w_list.a[k].y_end != -1) && (!(is_ualn_win(o->w_list.a[k])))), 
            o->w_list.a[k].x_start, o->w_list.a[k].x_end+1, 
            o->w_list.a[k].y_start, o->w_list.a[k].y_end+1, o->w_list.a[k].error);
    }
}

///[ts, te) ->  this is the reverse coordinates of t, not the original coordinates of t
int64_t extract_subov_cigar(const ul_idx_t *uref, char* qstr, UC_Read *tu, bit_extz_t *exz,
int64_t ts0, int64_t te0, overlap_region *o, double o_rate, int64_t *in_k, int64_t ql, 
double e_rate, uint64_t rid, uint64_t dbg_id, rtrace_t *trace, ul_ov_t *res)
{
    int64_t rev = o->y_pos_strand, t[2], q[2], k = 0, wts, wte, ii[2];
    int64_t wn = o->w_list.n, os, oe, ovlp, salnl = 0;
    t[0] = ts0; t[1] = te0; if(wn <= 0) return 0;

    if(in_k) k = *in_k;
    if(k < 0) k = 0; if(k >= wn) k = wn-1;
    for(; k < wn && t[0] > o->w_list.a[k].y_end; k++);
    if(k < 0) k = 0; if(k >= wn) k = wn-1;
    for(; k >= 0 && t[0] < o->w_list.a[k].y_start; k--);
    ///qs <= o->w_list.a[k].x_end && qs >= o->w_list.a[k].x_start
    if(k < 0) k = 0;
    if(in_k) *in_k = k;

    for (ii[0] = INT32_MAX, ii[1] = -1; k < wn; k++) {
        wts = o->w_list.a[k].y_start;
        wte = o->w_list.a[k].y_end + 1;
        if(wts >= t[1]) break;
        if((o->w_list.a[k].y_end == -1) || (is_ualn_win(o->w_list.a[k]))) continue;
        os = MAX(t[0], wts); oe = MIN(t[1], wte);
        ovlp = ((oe>os)? (oe-os):0);
        if(!ovlp) continue;
        if(k < ii[0]) ii[0] = k;
        if(k > ii[1]) ii[1] = k;
        salnl += ovlp; 
    }
    // if(dbg_id == 5525) {
    //     prt_aln_w(o);
    //     fprintf(stderr, "[M::%s::aln->%ld] ii::[%ld, %ld), q_aln0::[%d, %d), t0::[%ld, %ld)\n", 
    //     __func__, salnl, ii[0], ii[1]+1, o->w_list.a[ii[0]].x_start, o->w_list.a[ii[0]].x_end+1, 
    //     ts0, te0);
    // }
    if((!salnl) || (ii[0] == INT32_MAX) || (ii[1] < 0)) return 0;
    if((ii[0] == ii[1]) && (is_ualn_win(o->w_list.a[ii[0]]))) return 0;
    //[ii[0], ii[1]]
    win_boundary_offset(o->w_list.a, o->w_list.n, ii[0], ts0, ql, &(q[0]), &(t[0]));
    win_boundary_offset(o->w_list.a, o->w_list.n, ii[1], te0-1, ql, &(q[1]), &(t[1]));
    q[1]++; t[1]++;
    // if(dbg_id == 5525) {
    //     fprintf(stderr, "[M::%s::aln->%ld] ii::[%ld, %ld), q::[%ld, %ld), t::[%ld, %ld)\n", 
    //             __func__, salnl, ii[0], ii[1]+1, q[0], q[1], t[0], t[1]);
    // }
    if(salnl < ((t[1]-t[0])*o_rate)) return 0;
    
    gen_raln(uref, qstr, tu, o, exz, ii[0], ii[1]+1, ts0, te0, ql, o->y_id, rev, e_rate, rid, trace, res);
    assert(res->ts >= ts0 && res->te <= te0 && res->qs >= 0 && res->qe <= ql);
    double simi_thre = e_rate + MIN(r_simi_w, (e_rate/2));
    if((res->sec <= ((res->qe-res->qs)*simi_thre)) && (res->sec <= ((res->te-res->ts)*simi_thre))) {
        return 1;
    } else {
        return 0;
    }
}


uint64_t gen_sub_ov_adv_cigar(const ul_idx_t *udb, overlap_region* o, char* qstr, UC_Read *tu, 
bit_extz_t *exz, int64_t ql, double o_rate, double e_rate, utg_ct_t *ct_a, uint64_t ct_n, 
uint64_t sid, uint64_t oid, kv_rtrace_t *trace, kv_ul_ov_t *res)
{
	uint64_t ts, te, i, l, rn = res->n, rev = o->y_pos_strand, s, e, rid, t[2]; 
	ma_utg_t *u = &(udb->ug->u.a[o->y_id]); ul_ov_t z; rtrace_t tz; int64_t k;
    if(!rev){
        ts = o->y_pos_s; te = o->y_pos_e + 1; k = 0;
    } else {
        ts = u->len - (o->y_pos_e+1); te = u->len - o->y_pos_s; 
        k = ((int64_t)o->w_list.n)-1; if(k < 0) k = 0;
    }

    if(!ct_a) {    
        for (i = l = 0; i < u->n; i++) {
            rid = u->a[i]>>33; 
            s = l; e = l + Get_READ_LENGTH(R_INF, rid);///note: [s, e) are pos of t
            l += (uint32_t)u->a[i];
            if(e <= ts) continue; 
            if(s >= te) break;
            t[0] = (rev?(u->len-e):(s)); t[1] = (rev?(u->len-s):(e));
            if(extract_subov_cigar(udb, qstr, tu, exz, t[0], t[1], o, o_rate, &k, ql, e_rate, sid, rid, &tz, &z)) {
                tz.oid = oid; 
                z.ts -= t[0]; z.te -= t[0];
                z.rev = ((o->y_pos_strand == ((u->a[i]>>32)&1))?0:1);
                if(z.rev) {
                    t[0] = z.ts; t[1] = z.te;
                    z.ts = Get_READ_LENGTH(R_INF, rid) - t[1];
                    z.te = Get_READ_LENGTH(R_INF, rid) - t[0];
                }
                ///non-contained read at the unitg
                z.el = 1; 
                ///rid
                z.tn = rid;
                ///i-th read at the unitig
                z.qn = trace->n;
                kv_push(ul_ov_t, *res, z);
                kv_push(rtrace_t, *trace, tz);
            //     fprintf(stderr, "+[M::%s::rid->%lu::rev->%lu] utg_t::[%ld, %ld), ql::%ld\n", 
            // __func__, rid, rev, t[0], t[1], ql);
                // fprintf(stderr, "+[M::%s::%.*s::%c] q::[%u, %u), ql::%ld, t::[%u, %u), tl::%lu, err::%u\n", 
                // __func__, (int)Get_NAME_LENGTH(R_INF, z.tn), Get_NAME(R_INF, z.tn), 
                // "+-"[z.rev], z.qs, z.qe, ql, z.ts, z.te, Get_READ_LENGTH(R_INF, z.tn), z.sec);
            }
        }
    } else {
        for (i = 0; i < ct_n; i++) {
            rid = ct_a[i].x>>1; s = ct_a[i].s; e = ct_a[i].e;
            if(e <= ts) continue; 
            if(s >= te) break;
            t[0] = (rev?(u->len-e):(s)); t[1] = (rev?(u->len-s):(e));
            if(extract_subov_cigar(udb, qstr, tu, exz, t[0], t[1], o, o_rate, &k, ql, e_rate, sid, rid, &tz, &z)) {
                tz.oid = oid;
                z.ts -= t[0]; z.te -= t[0];
                z.rev = ((o->y_pos_strand == (ct_a[i].x&1))?0:1);
                if(z.rev) {
                    t[0] = z.ts; t[1] = z.te;
                    z.ts = Get_READ_LENGTH(R_INF, rid) - t[1];
                    z.te = Get_READ_LENGTH(R_INF, rid) - t[0];
                }
                ///contained read at the unitg
                z.el = 0; 
                ///rid
                z.tn = rid;
                ///i-th read at the unitig
                z.qn = trace->n;
                kv_push(ul_ov_t, *res, z);
                kv_push(rtrace_t, *trace, tz);
            //     fprintf(stderr, "-[M::%s::rid->%lu::rev->%lu] utg_t::[%ld, %ld), ql::%ld\n", 
            // __func__, rid, rev, t[0], t[1], ql);
                // fprintf(stderr, "-[M::%s::%.*s::%c] q::[%u, %u), ql::%ld, t::[%u, %u), tl::%lu, err::%u\n", 
                // __func__, (int)Get_NAME_LENGTH(R_INF, z.tn), Get_NAME(R_INF, z.tn), 
                // "+-"[z.rev], z.qs, z.qe, ql, z.ts, z.te, Get_READ_LENGTH(R_INF, z.tn), z.sec);
            }
        }
    }
    return res->n-rn;
}

uint64_t win_cluster_fliter(window_list *wa, uint64_t wn, uint64_t min_ovlp, double o_rate)
{
    uint64_t k, ws = (uint64_t)-1, we = (uint64_t)-1, aln_ol, sk;
    for (k = 0; k < wn && wa[k].y_end == -1; k++);
    if(k >= wn) return 0; ws = wa[k].x_start;

    for (aln_ol = 0, sk = k; k < wn; k++) {
        we = wa[k].x_end+1; 
        if(wa[k].y_end != -1) {
            aln_ol += wa[k].x_end+1-wa[k].x_start;
        }
        if(we >= ws + min_ovlp) break;
    }
    // fprintf(stderr, "[M::%s::] aln_ol::%lu, w::[%lu, %lu), min_ovlp::%lu, o_rate::%f\n", 
    // __func__, aln_ol, ws, we, min_ovlp, o_rate);
    if((aln_ol >= (we-ws)*o_rate) && (we-ws >= min_ovlp)) return 1;

    for (k++; k < wn; k++) {
        we = wa[k].x_end+1; 
        if(wa[k].y_end != -1) {
            aln_ol += wa[k].x_end+1-wa[k].x_start;
        }
        for (; sk < wn && wa[sk].x_start + min_ovlp < we; sk++) {
            if(wa[sk].y_end != -1) aln_ol -= wa[sk].x_end+1-wa[sk].x_start;
        }
        ws = wa[sk].x_start; 
        // fprintf(stderr, "+[M::%s::] aln_ol::%lu, w::[%lu, %lu), min_ovlp::%lu, o_rate::%f, k::[%lu, %lu)\n", 
        //     __func__, aln_ol, ws, we, min_ovlp, o_rate, sk, k);
        if((aln_ol >= (we-ws)*o_rate) && (we-ws >= min_ovlp)) return 1;
        if((sk > 0) && (we-ws < min_ovlp)) {
            sk--; if(wa[sk].y_end != -1) aln_ol += wa[sk].x_end+1-wa[sk].x_start;
            ws = wa[sk].x_start;  
            // fprintf(stderr, "-[M::%s::] aln_ol::%lu, w::[%lu, %lu), min_ovlp::%lu, o_rate::%f, k::[%lu, %lu)\n", 
            // __func__, aln_ol, ws, we, min_ovlp, o_rate, sk, k);
            if((aln_ol >= (we-ws)*o_rate) && (we-ws >= min_ovlp)) return 1;
        }
    }
    return 0;
}


void print_aln_windows(overlap_region *z)
{
    uint64_t k;
    for (k = 0; k < z->w_list.n; k++) {
        fprintf(stderr, "[M::%s::] q::[%u, %u), t::[%u, %u), err::%d\n", __func__, 
                z->w_list.a[k].x_start, z->w_list.a[k].x_end+1, 
                z->w_list.a[k].y_start, z->w_list.a[k].y_end+1,
                z->w_list.a[k].error);
    }
}

uint64_t ul_local_aln(overlap_region *z, Candidates_list *cl, const ul_idx_t *udb, char* qstr, UC_Read *tu, 
bit_extz_t *exz, double e_rate, int64_t w_l, uint64_t min_ovlp, double o_rate, kv_ul_ov_t *rln, kv_ul_ov_t *cln,
kv_rtrace_t *trace, uint64_t ql, uint64_t rid, uint64_t oid, uint64_t khit, overlap_region *aux_o)
{   
    uint64_t ol, aln_ol, cn, k, mm, rln_0; ul_contain *ct = udb->ct; utg_ct_t *ca; ul_ov_t *p;
    if(!aux_o) {
        align_ul_ed_post_extz(z, udb, NULL, qstr, tu->seq, exz, e_rate, w_l, -1, 1, NULL); 
        ol = z->x_pos_e+1-z->x_pos_s, aln_ol = z->align_length;
        if(aln_ol >= (ol*o_rate)) return 1;
        if(aln_ol <= min_ovlp) return 0;
        // fprintf(stderr, "***[M::%s::aln_ol->%lu] utg%.6dl(%c), align::%u, q::[%u, %u), t::[%u, %u)\n", __func__, aln_ol,
        //         (int32_t)z->y_id + 1, "+-"[z->y_pos_strand], z->align_length,
        //         z->x_pos_s, z->x_pos_e+1, z->y_pos_s, z->y_pos_e+1);
        mm = win_cluster_fliter(z->w_list.a, z->w_list.n, min_ovlp, o_rate/2);
        // print_aln_windows(z);
        // fprintf(stderr, "***[M::%s::] mm::%lu\n", __func__, mm);
        return mm;
    } else {
        ol = z->x_pos_e+1-z->x_pos_s; aln_ol = z->align_length;
        if(aln_ol <= min_ovlp) return 0;
        // if(z->y_id == 700) {
        //     fprintf(stderr, "+[M::%s::aln_ol->%lu::o_rate->%f] utg%.6dl(%c), align::%u, q::[%u, %u), t::[%u, %u), ql::%lu\n", 
        //             __func__, aln_ol, o_rate, 
        //             (int32_t)z->y_id + 1, "+-"[z->y_pos_strand], z->align_length,
        //             z->x_pos_s, z->x_pos_e+1, z->y_pos_s, z->y_pos_e+1, ql);
        //     print_aln_windows(z);
        // }
        rln_0 = rln->n;
        if((ol*o_rate) <= aln_ol) {
            kv_pushp(ul_ov_t, *rln, &p); memset(p, 0, sizeof(*p));
            ///[ts, te) -> whole interval rid at the unitig adjusted by the reverse 
            p->ts = z->y_pos_s; p->te = z->y_pos_e+1; 
            p->qs = z->x_pos_s; p->qe = z->x_pos_e+1;
            p->rev = z->y_pos_strand; p->el = 0; p->tn = p->qn = (uint32_t)-1;
        } else {
            cn = ((uint32_t)(ct->idx.a[z->y_id])); 
            ca = ct->rids.a + ((ct->idx.a[z->y_id])>>32);
            gen_sub_ov_adv(udb, z, o_rate, NULL, 0, rln);
            gen_sub_ov_adv(udb, z, o_rate, ca, cn, rln);
        }

        if(rln->n <= rln_0) return 0;
        radix_sort_ul_ov_srt_qs1(rln->a+rln_0, rln->a+rln->n);
        // if(z->y_id == 700) {
        //     for (k = 0; k < rln->n; k++) {
        //         p = &(rln->a[k]);
        //         fprintf(stderr, "+[M::%s::k->%lu] candidate_q::[%u, %u)\n", __func__, k, p->qs, p->qe);
        //     }
        // }

        for (k = mm = rln_0, p = NULL; k < rln->n; k++) {
            if((!p) || (rln->a[k].qs >= p->qe)) p = NULL;
            if(p) {
                if((rln->a[k].qs<p->qe)&&(rln->a[k].qe>p->qe)) {
                    p->qe = rln->a[k].qe; 
                }
            } else {
                p = &(rln->a[mm++]); *p = rln->a[k];
            }
        }
        rln->n = mm;
        // if(z->y_id == 700) {
        //     for (k = 0; k < rln->n; k++) {
        //         p = &(rln->a[k]);
        //         fprintf(stderr, "-[M::%s::k->%lu] candidate_q::[%u, %u)\n", __func__, k, p->qs, p->qe);
        //     }
        // }

        return_t_chain(z, cl);
        cigar_gen_by_chain_adv_local(z, cl, rln->a+rln_0, rln->n-rln_0, w_l, udb, NULL, NULL, qstr, tu, exz, aux_o, e_rate, ql, rid, khit);
        
        // fprintf(stderr, "-[M::%s::aln_ol->%lu] utg%.6dl(%c), align::%u, q::[%u, %u), t::[%u, %u)\n", __func__, aln_ol,
        //         (int32_t)z->y_id + 1, "+-"[z->y_pos_strand], z->align_length,
        //         z->x_pos_s, z->x_pos_e+1, z->y_pos_s, z->y_pos_e+1);
        // print_aln_windows(z);


        rln->n = rln_0;
        cn = ((uint32_t)(ct->idx.a[z->y_id])); 
        ca = ct->rids.a + ((ct->idx.a[z->y_id])>>32);
        gen_sub_ov_adv_cigar(udb, z, qstr, tu, exz, ql, o_rate, e_rate, NULL, 0, rid, oid, trace, rln);
        gen_sub_ov_adv_cigar(udb, z, qstr, tu, exz, ql, o_rate, e_rate, ca, cn, rid, oid, trace, cln);
        return 1;
    }
}


int64_t ul_raw_aln(overlap_region *z, Candidates_list *cl, const ul_idx_t *udb, char* qstr, UC_Read *tu, 
bit_extz_t *exz, double e_rate, int64_t w_l, uint64_t ql, uint64_t rid, uint64_t khit, overlap_region *aux_o)
{   
    return_t_chain(z, cl);
    int64_t ch_idx = z->shared_seed, ch_n;
    int64_t i, tl, id = z->y_id, m; ul_ov_t ov;
    k_mer_hit *ch_a = cl->list + ch_idx;
    tl = udb->ug->u.a[id].len;
    for (i = ch_idx; i < cl->length && cl->list[i].readID == cl->list[ch_idx].readID; i++); ch_n = i-ch_idx;

    // if(z->x_id == 57 && z->y_id == 2175) {
    //     fprintf(stderr, "\n-*-[M::%s]\tutg%.6u%c\txl::%u\tx::[%u,\t%u)\t%c\tutg%.6u%c\tyl::%u\ty::[%u,\t%u)\tch_idx::%ld\tch_n::%ld\n", 
    //         __func__, 
	// 		z->x_id+1, "lc"[udb->ug->u.a[z->x_id].circ], udb->ug->u.a[z->x_id].len, z->x_pos_s, z->x_pos_e+1,
	// 		"+-"[z->y_pos_strand],
	// 		z->y_id+1, "lc"[udb->ug->u.a[z->y_id].circ], udb->ug->u.a[z->y_id].len, z->y_pos_s, z->y_pos_e+1,
    //         ch_idx, ch_n);
    //     for (i = ch_idx; i < cl->length && cl->list[i].readID == cl->list[ch_idx].readID; i++) {
    //         fprintf(stderr, "i::%ld[M::%s]\treadID::%u\tself_offset::%u\toffset::%u\t%c\n", 
    //             i, __func__, cl->list[i].readID, cl->list[i].self_offset, cl->list[i].offset, 
    //             "+-"[cl->list[i].strand]);
    //     }
    // }

    // on = fusion_chain_ovlp(z, ch_a, ch_n, ov, on, wl, ql, tl);
    aux_o->w_list.n = aux_o->w_list.c.n = 0; 
    aux_o->y_id = z->y_id; aux_o->y_pos_strand = z->y_pos_strand;
    aux_o->x_pos_s = z->x_pos_s; aux_o->x_pos_e = z->x_pos_e;
    aux_o->y_pos_s = z->y_pos_s; aux_o->y_pos_e = z->y_pos_e;

    memset(&ov, 0, sizeof(ov)); ov.sec = 3;
    ov.ts = ov.te = (uint32_t)-1;
    ov.qs = 0; ov.qn = (uint32_t)-1; //extension to left
    ov.qe = ql; ov.tn = ch_n; //extension to right

    ovlp_base_direct(z, ch_a, ch_n, &ov, w_l, udb, NULL, NULL, qstr, tu, exz, aux_o, e_rate, ql, tl, rid);

    int64_t aux_n = aux_o->w_list.n;


    // if(z->x_id == 57 && z->y_id == 2175) {
    //     fprintf(stderr, "\n-0-[M::%s]\tutg%.6u%c\txl::%u\tx::[%u,\t%u)\t%c\tutg%.6u%c\tyl::%u\ty::[%u,\t%u)\n", 
    //         __func__, 
	// 		z->x_id+1, "lc"[udb->ug->u.a[z->x_id].circ], udb->ug->u.a[z->x_id].len, z->x_pos_s, z->x_pos_e+1,
	// 		"+-"[z->y_pos_strand],
	// 		z->y_id+1, "lc"[udb->ug->u.a[z->y_id].circ], udb->ug->u.a[z->y_id].len, z->y_pos_s, z->y_pos_e+1);
    //     for (i = 0; i < aux_n; i++) {
    //         window_list *m = &(aux_o->w_list.a[i]);
    //         fprintf(stderr, "i::%ld[M::%s]\tutg%.6u%c\twx::[%u,\t%u)\t%c\tutg%.6u%c\twy::[%u,\t%u)\terr::%d\tualn::%u\test::%u\n", i, __func__, 
	// 			z->x_id+1, "lc"[udb->ug->u.a[z->x_id].circ], m->x_start, m->x_end+1,
	// 			"+-"[z->y_pos_strand], 
	// 			z->y_id+1, "lc"[udb->ug->u.a[z->y_id].circ], m->y_start, m->y_end+1, m->error, 
	// 			(is_ualn_win((*m))), (is_est_aln((*m))));
    //     }
    //     fprintf(stderr, "\n");
    // }
    

    for (i = 0; i < aux_n; i++) {
        if(!(is_ualn_win(aux_o->w_list.a[i]))) continue;
        //will overwrite ch_a; does not matter
        rechain_aln(z, cl, aux_o, i, w_l, udb, NULL, NULL, qstr, tu, exz, e_rate, ql, tl, khit, rid);
    }

    // if(z->x_id == 57 && z->y_id == 2175) {
    //     fprintf(stderr, "\n-1-[M::%s]\tutg%.6u%c\txl::%u\tx::[%u,\t%u)\t%c\tutg%.6u%c\tyl::%u\ty::[%u,\t%u)\n", 
    //         __func__, 
	// 		z->x_id+1, "lc"[udb->ug->u.a[z->x_id].circ], udb->ug->u.a[z->x_id].len, z->x_pos_s, z->x_pos_e+1,
	// 		"+-"[z->y_pos_strand],
	// 		z->y_id+1, "lc"[udb->ug->u.a[z->y_id].circ], udb->ug->u.a[z->y_id].len, z->y_pos_s, z->y_pos_e+1);
    //     for (i = 0; i < ((int64_t)aux_o->w_list.n); i++) {
    //         window_list *m = &(aux_o->w_list.a[i]);
    //         fprintf(stderr, "i::%ld[M::%s]\tutg%.6u%c\twx::[%u,\t%u)\t%c\tutg%.6u%c\twy::[%u,\t%u)\terr::%d\tualn::%u\test::%u\n", i, __func__, 
	// 			z->x_id+1, "lc"[udb->ug->u.a[z->x_id].circ], m->x_start, m->x_end+1,
	// 			"+-"[z->y_pos_strand], 
	// 			z->y_id+1, "lc"[udb->ug->u.a[z->y_id].circ], m->y_start, m->y_end+1, m->error, 
	// 			(is_ualn_win((*m))), (is_est_aln((*m))));
    //     }
    //     fprintf(stderr, "\n");
    // }
    
    if(((int64_t)aux_o->w_list.n) > aux_n) {
        for (i = m = 0; i < ((int64_t)aux_o->w_list.n); i++) {
            if((i < aux_n) && (is_ualn_win(aux_o->w_list.a[i]))) continue;
            aux_o->w_list.a[m++] = aux_o->w_list.a[i];
        }
        aux_o->w_list.n = m;
        radix_sort_window_list_xs_srt(aux_o->w_list.a, aux_o->w_list.a+aux_o->w_list.n);
    }

    ///update z by aux_o
    update_overlap_region(z, aux_o, ql, tl);
    assert((z->x_pos_e>=z->x_pos_s) && (z->y_pos_e>=z->y_pos_s));

    int64_t zwn = z->w_list.n, zerr = 0, zlen = z->x_pos_e+1-z->x_pos_s;
    for (i = 0; i < zwn; i++) {
        if(is_ualn_win(z->w_list.a[i])) {
            zerr += MAX((z->w_list.a[i].x_end+1-z->w_list.a[i].x_start), 
                                (z->w_list.a[i].y_end+1-z->w_list.a[i].y_start));
        } else {
            zerr += z->w_list.a[i].error;
        }
    }
    z->non_homopolymer_errors = zerr;

    // if(z->x_id == 57 && z->y_id == 2175) {
    //     fprintf(stderr, "\n-2-[M::%s]\tutg%.6u%c\txl::%u\tx::[%u,\t%u)\t%c\tutg%.6u%c\tyl::%u\ty::[%u,\t%u)\n", 
    //         __func__, 
	// 		z->x_id+1, "lc"[udb->ug->u.a[z->x_id].circ], udb->ug->u.a[z->x_id].len, z->x_pos_s, z->x_pos_e+1,
	// 		"+-"[z->y_pos_strand],
	// 		z->y_id+1, "lc"[udb->ug->u.a[z->y_id].circ], udb->ug->u.a[z->y_id].len, z->y_pos_s, z->y_pos_e+1);
    //     for (i = 0; i < zwn; i++) {
    //         window_list *m = &(z->w_list.a[i]);
    //         fprintf(stderr, "i::%ld[M::%s]\tutg%.6u%c\twx::[%u,\t%u)\t%c\tutg%.6u%c\twy::[%u,\t%u)\terr::%d\tualn::%u\test::%u\n", i, __func__, 
	// 			z->x_id+1, "lc"[udb->ug->u.a[z->x_id].circ], m->x_start, m->x_end+1,
	// 			"+-"[z->y_pos_strand], 
	// 			z->y_id+1, "lc"[udb->ug->u.a[z->y_id].circ], m->y_start, m->y_end+1, m->error, 
	// 			(is_ualn_win((*m))), (is_est_aln((*m))));
    //     }
    //     fprintf(stderr, "\n");
    // }

    // if(z->x_id == 57 && z->y_id == 2175) {
    //     fprintf(stderr, "\n-3-[M::%s]\tutg%.6u%c\txl::%u\tx::[%u,\t%u)\t%c\tutg%.6u%c\tyl::%u\ty::[%u,\t%u)\tzerr::%ld\tzlen::%ld\te_rate::%f\n", 
    //         __func__, 
	// 		z->x_id+1, "lc"[udb->ug->u.a[z->x_id].circ], udb->ug->u.a[z->x_id].len, z->x_pos_s, z->x_pos_e+1,
	// 		"+-"[z->y_pos_strand],
	// 		z->y_id+1, "lc"[udb->ug->u.a[z->y_id].circ], udb->ug->u.a[z->y_id].len, z->y_pos_s, z->y_pos_e+1,
    //         zerr, zlen, e_rate);
    // }

    if(zerr >= zlen) return 0;
    if(zerr <= 0 && zlen > 0) return 1;
    if(zerr > (zlen*e_rate)) return 0;
    return 1;
}

void dedup_ul_ov_t(kv_ul_ov_t *in)
{
    ul_ov_t *a = in->a; int64_t k, l, m, z, r, a_n = in->n; uint64_t qo, to; double rr = 0.95;
    for (k = 0; k < a_n; k++) a[k].tn = ((uint32_t)(a[k].tn<<1))|((uint32_t)(a[k].rev));
    radix_sort_ul_ov_srt_tn1(a, a + a_n);
    for (k = 1, l = m = 0; k <= a_n; k++) {
        if(k == a_n || a[k].tn != a[l].tn) {
            for (z = l; z < k; z++) {
                for (r=m-1; (r>=0) && (a[r].tn==(a[z].tn>>1)) && (a[r].rev==a[z].rev); r--) {
                    qo = ((MIN(a[z].qe, a[r].qe) > MAX(a[z].qs, a[r].qs))? 
															MIN(a[z].qe, a[r].qe) - MAX(a[z].qs, a[r].qs):0);
                    to = ((MIN(a[z].te, a[r].te) > MAX(a[z].ts, a[r].ts))?
                                                        MIN(a[z].te, a[r].te) - MAX(a[z].ts, a[r].ts):0);
                    if(qo >= ((a[r].qe - a[r].qs)*rr) && qo >= ((a[z].qe - a[z].qs)*rr) && 
						   to >= ((a[r].te - a[r].ts)*rr) && to >= ((a[z].te - a[z].ts)*rr)) {
                            break;
                    }
                }
                if(r >= 0 && (a[r].tn==(a[z].tn>>1)) && (a[r].rev==a[z].rev)) {
                    if(a[z].sec < a[r].sec) {
                        a[r] = a[z]; a[r].tn >>= 1; 
                    }
                    continue;
                }
                a[m] = a[z]; a[m].tn >>= 1; m++;
            }
            l = k;
        }
    }
    // fprintf(stderr, "[M::%s::] in->n0::%ld, in->n::%ld\n", __func__, (int64_t)in->n, m);
    in->n = m;
}

void ul_rid_lalign_adv(overlap_region_alloc* ol, Candidates_list *cl, const ul_idx_t *uref, const ug_opt_t *uopt, 
        char *qstr, uint64_t ql, UC_Read* qu, UC_Read* tu, bit_extz_t *exz, overlap_region *aux_o, double e_rate, 
        int64_t wl, kv_ul_ov_t *aln, kv_ul_ov_t *cln, kv_rtrace_t *trace, int64_t sid, uint64_t khit, void *km)
{
    uint64_t i, bs, k; Window_Pool w; double err; 
    overlap_region t; overlap_region *z; //asg64_v iidx, buf, buf1;
    ol->mapped_overlaps_length = 0;
    if(ol->length <= 0) return;

    ///base alignment
    err = e_rate; 
    init_Window_Pool(&w, ql, wl, (int)(1.0/err));
    bs = (w.window_length)+(THRESHOLD_MAX_SIZE<<1)+1;
    resize_UC_Read(tu, bs<<1); 

    if(!aux_o) {
        resize_UC_Read(qu, ql); qu->length = ql; memcpy(qu->seq, qstr, ql); 
        for (i = k = 0; i < ol->length; i++) {
            z = &(ol->list[i]); z->shared_seed = z->non_homopolymer_errors;///for index
            if(!ul_local_aln(z, cl, uref, qu->seq, tu, exz, err, w.window_length, 
            1000, OVERLAP_THRESHOLD_NOSI_FILTER, NULL, NULL, NULL, ql, sid, i, khit, NULL)) {
                continue;
            }
            if(k != i) {
                t = ol->list[k];
                ol->list[k] = ol->list[i];
                ol->list[i] = t;
            }
            z = &(ol->list[k++]); z->is_match = 1; 
        }
        ol->length = k;
        if(ol->length <= 0) return;
    } else {
        for (i = cln->n = trace->n = 0; i < ol->length; i++) {
            z = &(ol->list[i]); z->shared_seed = z->non_homopolymer_errors;///for index
            // fprintf(stderr, "\n[M::%s::utg%.6dl(%c)] i::%lu, aln_l::%u, q::[%u, %u), ql::%u\n", __func__, 
            // (int32_t)z->y_id + 1, "+-"[z->y_pos_strand], i, z->align_length, 
            // z->x_pos_s, z->x_pos_e+1, z->x_pos_e+1-z->x_pos_s);
            ul_local_aln(z, cl, uref, qu->seq, tu, exz, err, w.window_length, 1000, 
            OVERLAP_THRESHOLD_NOSI_FILTER, aln, cln, trace, ql, sid, i, khit, aux_o);
        }

        ///contained reads
        if(cln->n) dedup_ul_ov_t(cln);
        if(cln->n) {
            kv_resize(ul_ov_t, *aln, aln->n+cln->n);
            memcpy(aln->a+aln->n, cln->a, cln->n*sizeof(*(aln->a)));
            aln->n += cln->n;
        }

    }
}

/**
void ul_raw_lalign_adv(overlap_region_alloc* ol, Candidates_list *cl, const ul_idx_t *uref, All_reads *rdb, const ug_opt_t *uopt, 
        char *qstr, uint64_t ql, UC_Read* qu, UC_Read* tu, Correct_dumy* dumy, bit_extz_t *exz, haplotype_evdience_alloc* hap, 
        kvec_t_u64_warp* v_idx, overlap_region *aux_o, double e_rate, int64_t wl, kv_ul_ov_t *aln, kv_ul_ov_t *aln1,
        int64_t sid, uint64_t khit, st_mt_t *stb, void *km)
{
    uint64_t i, bs, k, aln_occ; Window_Pool w; double err; 
    overlap_region t; overlap_region *z; asg64_v iidx, buf, buf1;
    ol->mapped_overlaps_length = 0;
    if(ol->length <= 0) return;

    ///base alignment
    clear_Correct_dumy(dumy, ol, km); err = e_rate; 
    init_Window_Pool(&w, ql, wl, (int)(1.0/err));
    bs = (w.window_length)+(THRESHOLD_MAX_SIZE<<1)+1;
    resize_UC_Read(tu, bs<<1); 

    if(!aux_o) {
        resize_UC_Read(qu, ql); qu->length = ql; memcpy(qu->seq, qstr, ql); 
        copy_asg_arr(iidx, hap->snp_srt);
        for (i = k = 0, aln->n = 0; i < ol->length; i++) {
            z = &(ol->list[i]); z->shared_seed = z->non_homopolymer_errors;///for index
            align_ul_ed_post_extz(z, uref, NULL, qu->seq, tu->seq, exz, err, w.window_length, -1, 1, km);
            aln_occ = gen_r_aln(uref, z, k, &iidx, OVERLAP_THRESHOLD_FILTER, aln, 1000);
            if(!aln_occ) continue;
            if(k != i) {
                t = ol->list[k];
                ol->list[k] = ol->list[i];
                ol->list[i] = t;
            }
            z = &(ol->list[k++]); z->is_match = 1; 
        }
        copy_asg_arr(hap->snp_srt, iidx);
        ol->length = k;
        if(ol->length <= 0) return;
    } else {
        copy_asg_arr(iidx, hap->snp_srt); copy_asg_arr(buf, v_idx->a); copy_asg_arr(buf1, (*stb));
        ul_gap_filling_local(ol, cl, aln, wl, uref, NULL, NULL, qu->seq, tu, exz, aux_o, &buf, &iidx, err, ql, sid, khit, 1, MAX_LGAP(ql));
        copy_asg_arr(hap->snp_srt, iidx); copy_asg_arr(v_idx->a, buf); copy_asg_arr((*stb), buf1);

        copy_asg_arr(iidx, hap->snp_srt);
        gen_aln_local(ol, aln, aln1, &iidx, uref, qu->seq, tu, exz, aux_o, err, ql, OVERLAP_THRESHOLD_FILTER);
        copy_asg_arr(hap->snp_srt, iidx);
    }
    // } else {
    //     // fprintf(stderr, "-[M::%s] on::%lu\n", __func__, ol->length);
    //     if(ol->length <= 1) return;
    //     ///coordinates for all intervals with cov > 1
    //     copy_asg_arr(iidx, hap->snp_srt); copy_asg_arr(buf, v_idx->a); copy_asg_arr(buf1, (*stb));
    //     // fprintf(stderr, "\n[M::%s] iidx_n::%ld\n", __func__, (int64_t)iidx.n);
    //     ul_gap_filling_adv(ol, cl, aln, wl, uref, NULL, NULL, qu->seq, tu, exz, aux_o, &buf, &iidx, err, ql, sid, khit, 1, MAX_LGAP(ql));
    //     copy_asg_arr(hap->snp_srt, iidx); copy_asg_arr(v_idx->a, buf); copy_asg_arr((*stb), buf1);

    //     copy_asg_arr(iidx, hap->snp_srt); copy_asg_arr(buf, v_idx->a); copy_asg_arr(buf1, (*stb));
    //     region_phase(ol, uref, uopt, aln, &iidx, &buf, &buf1);
    //     copy_asg_arr(hap->snp_srt, iidx); copy_asg_arr(v_idx->a, buf); copy_asg_arr((*stb), buf1);
    // }
}
**/

uint64_t gen_nkhits(Candidates_list *cl, overlap_region *z)
{
    uint64_t pid, kn; int64_t k;
    pid = cl->list[z->shared_seed].readID;kn = 0;
    for (k = z->shared_seed; k < cl->length && cl->list[k].readID == pid; k++) kn++;
    return kn;
}

void ug_lalign(overlap_region_alloc* ol, Candidates_list *cl, const ul_idx_t *uref, const ug_opt_t *uopt, 
        char *qstr, uint64_t ql, UC_Read* qu, UC_Read* tu, Correct_dumy* dumy, bit_extz_t *exz,  
        overlap_region *aux_o, double e_rate, int64_t wl, int64_t sid, uint64_t khit, 
        uint64_t chain_cut, void *km)
{
    uint64_t i, bs, k, ovl, pid, kl, kn; Window_Pool w; double err; 
    overlap_region t; overlap_region *z; char *in = qstr;
    ol->mapped_overlaps_length = 0;
    if(ol->length <= 0) return;

    ///base alignment
    clear_Correct_dumy(dumy, ol, km); err = e_rate; 
    init_Window_Pool(&w, ql, wl, (int)(1.0/err));
    bs = (w.window_length)+(THRESHOLD_MAX_SIZE<<1)+1;
    resize_UC_Read(tu, bs<<1);

    if(!aux_o) {    
        if(qu) {
            resize_UC_Read(qu, ql); qu->length = ql; memcpy(qu->seq, qstr, ql); in = qu->seq;
        } 
        for (i = 0; i < ol->length; i++) {
            z = &(ol->list[i]); ovl = z->x_pos_e+1-z->x_pos_s; 
            z->shared_seed = z->non_homopolymer_errors;///for index
            pid = cl->list[z->shared_seed].readID; kl = cl->length; kn = 0;
            for (k = z->shared_seed; k < kl && cl->list[k].readID == pid && kn < chain_cut; k++) kn++;
            if(kn < chain_cut) continue;
            
            if(!align_ul_ed_post_extz(z, uref, NULL, in/**qu->seq**/, tu->seq, exz, err, w.window_length, -1, 0, km)) {
                continue;
            }
            if(uref && simi_pass(ovl, z->align_length, uref?1:0, -1, NULL)) {
                z->is_match = 3; ol->mapped_overlaps_length += z->align_length;
            }
        }

        // if(uref && ol->mapped_overlaps_length > 0) {
        //     set_herror_win(ol, dumy, v_idx, err, ql, w.window_length);
        // }

        double e_max = err*1.5, rr; int64_t re;
        for (i = k = 0; i < ol->length; i++) {
            z = &(ol->list[i]); ovl = z->x_pos_e + 1 - z->x_pos_s;
            rr = gen_extend_err_exz(z, uref, NULL, NULL, in/**qu->seq**/, tu->seq, exz, NULL/**v_idx?v_idx->a.a:NULL**/, w.window_length, -1, err, (e_max+0.000001), &re);
            z->is_match = 0;///must be here;
            
            // if(z->x_id == 57 && z->y_id == 2175) { 
            //     fprintf(stderr, "+utg%.6u%c\t%u\t%u\t%u\t%c\tutg%.6u%c\t%u\t%u\t%u\talign_length::%u\terr::%f\trr::%f\twn::%u\tk::%lu\ti::%lu\n", 
            //         z->x_id+1, "lc"[uref->ug->u.a[z->x_id].circ], uref->ug->u.a[z->x_id].len, 
            //         z->x_pos_s, z->x_pos_e+1, "+-"[z->y_pos_strand],
            //         z->y_id+1, "lc"[uref->ug->u.a[z->y_id].circ], uref->ug->u.a[z->y_id].len, 
            //         z->y_pos_s, z->y_pos_e+1, z->align_length, err, rr, (uint32_t)z->w_list.n/**gen_nkhits(cl, z)**/,
            //         k, i);
            // }
            if (rr <= err) {
                if(k != i) {
                    t = ol->list[k];
                    ol->list[k] = ol->list[i];
                    ol->list[i] = t;
                }
                ol->list[k].is_match = 1; ol->list[k].non_homopolymer_errors = re;
                k++;
            } 
        }

        ol->length = k;
        // if(sid == 57) { 
        //     fprintf(stderr, "+utg%.6ld%c\tol->length::%lu\n", sid+1, "lc"[uref->ug->u.a[sid].circ], ol->length);
        // }
        if(ol->length <= 0) return;
    } else {
        // if(ol->length <= 1) return;
        if(ol->length <= 0) return;
        if(qu) {
            resize_UC_Read(qu, ql); qu->length = ql; memcpy(qu->seq, qstr, ql); in = qu->seq;
        } 
        for (i = k = 0; i < ol->length; i++) {
            z = &(ol->list[i]); 
            // if(z->x_id == 57 && z->y_id == 2175) { 
            //     fprintf(stderr, "-1-utg%.6u%c\t%u\t%u\t%u\t%c\tutg%.6u%c\t%u\t%u\t%u\talign_length::%u\terr::%f\twn::%u\tk::%lu\ti::%lu\n", 
            //         z->x_id+1, "lc"[uref->ug->u.a[z->x_id].circ], uref->ug->u.a[z->x_id].len, 
            //         z->x_pos_s, z->x_pos_e+1, "+-"[z->y_pos_strand],
            //         z->y_id+1, "lc"[uref->ug->u.a[z->y_id].circ], uref->ug->u.a[z->y_id].len, 
            //         z->y_pos_s, z->y_pos_e+1, z->align_length, err, (uint32_t)z->w_list.n/**gen_nkhits(cl, z)**/,
            //         k, i);
            // }
            if(!ul_raw_aln(z, cl, uref, in/**qu->seq**/, tu, exz, err, wl, ql, sid, khit, aux_o)) {
                z->is_match = 0;
                continue;
            }

            // if(z->x_id == 57 && z->y_id == 2175) { 
            //     fprintf(stderr, "-2-utg%.6u%c\t%u\t%u\t%u\t%c\tutg%.6u%c\t%u\t%u\t%u\talign_length::%u\terr::%f\twn::%u\tk::%lu\ti::%lu\n", 
            //         z->x_id+1, "lc"[uref->ug->u.a[z->x_id].circ], uref->ug->u.a[z->x_id].len, 
            //         z->x_pos_s, z->x_pos_e+1, "+-"[z->y_pos_strand],
            //         z->y_id+1, "lc"[uref->ug->u.a[z->y_id].circ], uref->ug->u.a[z->y_id].len, 
            //         z->y_pos_s, z->y_pos_e+1, z->align_length, err, (uint32_t)z->w_list.n/**gen_nkhits(cl, z)**/,
            //         k, i);
            // }
            if(k != i) {
                t = ol->list[k];
                ol->list[k] = ol->list[i];
                ol->list[i] = t;
            }
            ol->list[k].is_match = 1; 
            k++;
        }
        ol->length = k;

        // if(sid == 57) { 
        //     fprintf(stderr, "-utg%.6ld%c\tol->length::%lu\n", sid+1, "lc"[uref->ug->u.a[sid].circ], ol->length);
        // }
    }
}