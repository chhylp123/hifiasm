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


int ha_ov_type(const overlap_region *r, uint32_t len);


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
            z->w_list.a[m].extra_begin = z->w_list.a[m].extra_end = -1;
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
        if ((rref && (overlap_length*OVERLAP_THRESHOLD_FILTER <= z->align_length)) || 
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
        if(!simi_pass(ovl, mm_aln, uref?1:0, OVERLAP_THRESHOLD_FILTER, NULL)) continue;
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
                if(!simi_pass(ovl, mm_aln, uref?1:0, OVERLAP_THRESHOLD_FILTER, NULL)) break;
            }
            if(w_idx[i] != (uint64_t)-1) mm_ws = z->w_list.a[w_idx[i]].x_end+1;
        }

        if(i < nw) continue;
        if(uref && simi_pass(ovl, z->align_length, uref?1:0, OVERLAP_THRESHOLD_FILTER, NULL)) {
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
        if(simi_pass(ovl, z->align_length, 0, OVERLAP_THRESHOLD_FILTER, &e_rate)) {
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

double test_err_rate(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, UC_Read* g_read, 
char* qstr, char *tstr, char *tstr_1, Correct_dumy* dumy, kvec_t_u64_warp* v_idx, int64_t block_s, double e_rate)
{
    int64_t ovl = z->x_pos_e+1-z->x_pos_s, a_nw = z->w_list.n, i; double error_rate;
    window_list *p; int64_t y_id = z->y_id, y_strand = z->y_pos_strand;
    if(!simi_pass(ovl, z->align_length, 0, OVERLAP_THRESHOLD_FILTER, &e_rate)) return DBL_MAX;

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
    error_rate = non_trim_error_rate(z, rref, uref, v_idx, dumy, g_read, e_rate, block_s);
    return error_rate;
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
                        dumy, e_rate, block_s, OVERLAP_THRESHOLD_FILTER, NULL)) {
            continue;
        }
        if(uref && simi_pass(ovl, z->align_length, uref?1:0, OVERLAP_THRESHOLD_FILTER, NULL)) {
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
                        int64_t aux_beg, int64_t aux_end, double e_rate, int64_t block_s, uint32_t sec_check, double ovlp_cut, void *km)
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


uint32_t align_ul_ed_post_extz(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, char* qstr, char *tstr, bit_extz_t *exz, double e_rate, int64_t w_l, double ovlp_cut, void *km)
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
                                        t_tot_l, aux_beg, aux_end, e_rate, w_l, sec_check, ovlp_cut, km)) {
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

inline uint64_t scale_ed_thre(uint32_t err)
{
    uint64_t bd = (err<<1)+1, w; 
    w = (bd>>bitw); w <<= bitw; if(w < bd) w += bitwbit;
    err = (w-1)>>1; if(err > MAX_E) err = MAX_E;
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

int64_t cal_exz_infi(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, bit_extz_t *exz, char* qstr, UC_Read *tu, int64_t qs, int64_t qe, int64_t ts, int64_t te, int64_t thre, int64_t q_tot_l, int64_t mode)
{
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
        clear_align(*exz);
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

void hc_aln_exz(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, 
char* qstr, UC_Read *tu, int64_t qs, int64_t qe, int64_t ts, int64_t te, int64_t estimate_err,
int64_t mode, int64_t wl, bit_extz_t *exz, int64_t q_tot, double e_rate)
{
    int64_t thre, ql = qe - qs, thre0; 
    if(((ts == -1) && (te == -1))) mode = 3;///set to semi-global
    // fprintf(stderr, "[M::%s::ql::%ld] qs::%ld, qe::%ld, ts::%ld, te::%ld, mode::%ld, estimate_err::%ld, e_rate::%f", 
    //                                                     __func__, ql, qs, qe, ts, te, mode, estimate_err, e_rate);
    
    if(ql <= MAX_L && (estimate_err*1.2) <= MAX_E) {
        thre = scale_ed_thre(estimate_err); if(thre > ql) thre = ql;
        if(cal_exz_infi(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, thre, q_tot, mode)) {
            // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(+)\n", exz->err, exz->thre, thre);
            return;
        }

        thre0 = thre; thre = ql*e_rate;
        thre = scale_ed_thre(thre); if(thre > ql) thre = ql;
        if(thre > thre0) {
            if(cal_exz_infi(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, thre, q_tot, mode)) {
                // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(-)\n", exz->err, exz->thre, thre);
                return;
            }
        }

        thre0 = thre; thre <<= 1;
        thre = scale_ed_thre(thre); if(thre > ql) thre = ql;
        if(thre > thre0) {
            if(cal_exz_infi(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, thre, q_tot, mode)) {
                // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(-)\n", exz->err, exz->thre, thre);
                return;
            }
        }

        thre0 = thre; thre = ql*0.51;
        thre = scale_ed_thre(thre); if(thre > ql) thre = ql;
        if(thre > thre0) {
            if(cal_exz_infi(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, thre, q_tot, mode)) {
                // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(*)\n", exz->err, exz->thre, thre);
                return;
            }
        }
    }
    // fprintf(stderr, ", err::%d, thre::%d\n", INT32_MAX, exz->thre);
    
}

void sub_ciagar_gen(overlap_region *z, uint64_t s, uint64_t e, uint64_t wl, 
const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, UC_Read *tu, bit_extz_t *exz, double e_rate,
int64_t ql, uint64_t rid)
{
    uint64_t qs, qe, sid, eid, k, l, m, tot_e, c_e; int64_t q[2], t[2], mode;
    qs = (s/wl)*wl; if(qs < z->x_pos_s) qs = z->x_pos_s; if(qs > z->x_pos_e) return;
    qe = (e/wl)*wl; if(qe < e) qe += wl; if(qe > z->x_pos_e+1) qe = z->x_pos_e+1;
    if(qe <= 0) return;
    sid = get_win_id_by_s(z, qs, wl, NULL);
    eid = get_win_id_by_e(z, qe-1, wl, NULL) + 1;///must qe-1 instead of qe!!!!!!
    if(sid >= eid) return;
    // fprintf(stderr, "***[M::%s] s::%lu, e::%lu, n_qs::%lu, n_qe::%lu, z::[%u, %u), sid::%lu, eid::%lu, w_list.n::%lu\n", 
    // __func__, s, e, qs, qe, z->x_pos_s, z->x_pos_e+1, sid, eid, (uint64_t)z->w_list.n);
    for (k = sid+1, l = sid; k <= eid; k++) {//[sid, eid)
        if(k == eid || z->w_list.a[k].extra_end < 0) {
            if(k - l > 1 || z->w_list.a[l].extra_end >= 0) {
                q[0] = q[1] = t[0] = t[1] = -1; mode = -1; tot_e = 0;
                if(z->w_list.a[l].extra_end < 0) {
                    q[0] = z->w_list.a[l].x_end+1; 
                    if(z->w_list.a[l].y_end != -1) {
                        t[0] = z->w_list.a[l].y_end+1; 
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
                hc_aln_exz(z, uref, hpc_g, rref, qstr, tu, q[0], q[1], t[0], t[1], tot_e, mode, wl, exz, ql, e_rate); 
                
            }
            l = k;
        }
    }
}

void cigar_gen(overlap_region *z, ul_ov_t *ov, uint64_t on, uint64_t qn, uint64_t wl,  
const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, UC_Read *tu, bit_extz_t *exz, double e_rate, int64_t ql, uint64_t rid)
{
    uint64_t i, qs = (uint64_t)-1, qe = (uint64_t)-1;
    for (i = 0; i < on && ov[i].qn == qn; i++) {
        assert((i<=0)||(ov[i].qs >= ov[i-1].qe));
        if(ov[i].qs <= qe && qe != (uint64_t)-1) {
            qe = ov[i].qe;
        } else {
            if(qs != (uint64_t)-1) sub_ciagar_gen(z, qs, qe, wl, uref, hpc_g, rref, qstr, tu, exz, e_rate, ql, rid);
            qs = ov[i].qs; qe = ov[i].qe;
        }
    }
    if(qs != (uint64_t)-1) sub_ciagar_gen(z, qs, qe, wl, uref, hpc_g, rref, qstr, tu, exz, e_rate, ql, rid);
}

void ul_gap_filling(overlap_region_alloc* ol, kv_ul_ov_t *aln, uint64_t wl, 
const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, UC_Read *tu, bit_extz_t *exz, 
double e_rate, int64_t ql, uint64_t rid)
{
    int64_t i, on = ol->length;
    for (i = 0; i < on; i++) {
        if(ol->list[i].align_length == (uint32_t)-1) continue;
        cigar_gen(&(ol->list[i]), aln->a+ol->list[i].align_length, aln->n-ol->list[i].align_length, i, wl, 
        uref, hpc_g, rref, qstr, tu, exz, e_rate, ql, rid);
    }
}


///[ys, ye)
uint64_t get_win_aln(overlap_region *z, uint64_t wid, int64_t *ys, int64_t *ye, int64_t *err)
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

uint64_t dp_commen_sketch(kv_ul_ov_t *aln, overlap_region_alloc* ol, uint64_t *id_a, int64_t id_n, 
uint64_t *win_a, int64_t win_n, uint64_t *dp, int64_t n_skip, int64_t wl)
{
    if(!win_n) return 0;
    int64_t k, ws, ws0, i, m, wid, wid0, s, e, err, s0, e0, err0, sc, max_sc, p, long_sc, long_idx; 
    overlap_region *z; s = e = err = s0 = e0 = err0 = -1;
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

                    get_win_aln(z, wid, &s, &e, &err);
                    get_win_aln(z, wid0, &s0, &e0, &err0);
                    if(e0 > s) break;
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

        k = long_idx; 
        while (k >= 0) {
            win_a[k] |= ((uint64_t)0x8000000000000000); 
            k = ((dp[k]>>32)!=((uint32_t)-1))?(dp[k]>>32):(-1);
        }

        for (k = 0, m = 0; k < win_n; k++) {
            if(!(win_a[k]&((uint64_t)0x8000000000000000))) continue;
            win_a[m++] = win_a[k];
        }
        win_n = m;
    }

    for (k = 0; k < win_n; k++) {
        ws = win_a[k];
        for (i = 0; i < id_n; i++) {
            z = &(ol->list[aln->a[(uint32_t)id_a[i]].qn]);
            wid = get_win_id_by_s(z, ws, wl, NULL);
            z->w_list.a[wid].extra_end = -1; 
            get_win_aln(z, wid, &s, &e, &err);
            if(z->align_length < (uint32_t)e) z->align_length = e;///for conliner
        }
    }
    return win_n;
}

uint64_t gen_commen_sketch(All_reads *rref, const ul_idx_t *uref, overlap_region_alloc* ol, uint64_t *id_a, uint64_t id_n, uint64_t s, uint64_t e, uint64_t ql, uint64_t wl,
uint64_t *buf, uint64_t dp, char *str0, char *str1, kv_ul_ov_t *aln, asg64_v *trace)///[s, e)
{
    if(!id_n) return id_n;
    uint64_t i, m, k, rm_n = 0, buf_n = 0, qs, qe, wid, co, occ; char *qstring, *tstring;
    overlap_region *z; uint64_t ws, we; int64_t r_y[2], r_err, p_y[2], p_err;
    ///shrink [qs, qe)
    qs = (s/wl)*wl; if(qs < s) qs += wl; if(qs >= ql) return id_n;
    qe = (e/wl)*wl; if(qe >= ql) qe = ql;
    if(qs >= qe) return id_n;
    //idx_a[] is sorted by aln[].qs
    for (k = 0; k < id_n; k++) {
        if(aln->a[id_a[k]].qs<=qs && aln->a[id_a[k]].qe>=qe) {
            buf[buf_n++] = id_a[k];
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
            if(r_y[0] < z->align_length) continue;///not co-linear
            qstring = tstring = NULL;
            for (i = 1; i < buf_n; i++) {
                z = &(ol->list[aln->a[(uint32_t)buf[i]].qn]);
                wid = get_win_id_by_s(z, ws, wl, NULL);
                if(!get_win_aln(z, wid, &(p_y[0]), &(p_y[1]), &p_err)) break;
                if(p_y[0] < z->align_length) continue;///not co-linear

                if(((r_y[1]-r_y[0]) != (p_y[1]-p_y[0])) || (r_err != p_err)) break;
                if(r_err == 0) continue;
                if(!qstring) {
                    qstring = retrive_str_piece_exz(rref, uref, str0, r_y[0], r_y[1]-r_y[0], 
                                                ol->list[aln->a[(uint32_t)buf[0]].qn].y_pos_strand, ol->list[aln->a[(uint32_t)buf[0]].qn].y_id);
                }
                tstring = retrive_str_piece_exz(rref, uref, str1, p_y[0], p_y[1]-p_y[0], z->y_pos_strand, z->y_id);
                if(memcmp(str0, str1, (p_y[1]-p_y[0]))) {
                    // fprintf(stderr, "[M::%s::] qs::%ld, qe::%lu, ts::%ld, te::%lu, err::%ld, rts::%ld, rte::%lu, err::%ld\n", 
                    // __func__, ws, we, p_y[0], p_y[1], r_err, r_y[0], r_y[1], p_err);
                    // fprintf(stderr, "str0::%.*s\n", (int32_t)(r_y[1]-r_y[0]), str0);
                    // fprintf(stderr, "str1::%.*s\n", (int32_t)(p_y[1]-p_y[0]), str1);
                    break;
                }
            }
            if(i < buf_n) continue;
            for (i = 0, p_y[0] = p_y[1] = -1; i < buf_n; i++) {
                z = &(ol->list[aln->a[(uint32_t)buf[i]].qn]);
                wid = get_win_id_by_s(z, ws, wl, NULL);

                if(co) {
                    get_win_aln(z, wid, &(p_y[0]), &(p_y[1]), &p_err);
                    if((buf[i]>>32) > (uint64_t)p_y[0]) co = 0;
                    m = p_y[1]; m <<= 32; m |= (uint32_t)buf[i]; buf[i] = m;
                }
            }
            if(co) occ++;
            kv_push(uint64_t, *trace, ws);
        }

        if(trace->n) {
            if(!co) kv_resize(uint64_t, *trace, trace->n<<1);
            trace->n = dp_commen_sketch(aln, ol, buf, buf_n, trace->a, trace->n, trace->a + trace->n, occ, wl);
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
All_reads *rref, UC_Read* tu, asg64_v* idx, asg64_v *b0, asg64_v *b1, int64_t ql, int64_t wl, kv_ul_ov_t *aln, uint64_t rid)
{
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
        ol->list[aln->a[i].qn].align_length = 0;///for co-linear
        idx->a[srt_n] = aln->a[i].qs<<1; idx->a[srt_n] <<= 32; idx->a[srt_n] += i; srt_n++;
        idx->a[srt_n] = ((aln->a[i].qe-1)<<1)+1; idx->a[srt_n] <<= 32; idx->a[srt_n] += i; srt_n++;
		aln->a[i].el = 1;
    }

    radix_sort_bc64(idx->a, idx->a+srt_n); idx->n = srt_n; resize_UC_Read(tu, (wl<<1));
    for (i = 0, dp = 0, beg = 0, end = -1; i < srt_n; ++i) {///[beg, end]
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
			idx->n = srt_n + 
					gen_commen_sketch(rref, uref, ol, idx->a+srt_n, idx->n-srt_n, beg, end, ql, wl, b0->a, old_dp, tu->seq, tu->seq+wl, aln, b1);
		}
        beg = end;
	}

    for (i = srt_n = 0; i < aln->n; i++) {
		if(i == 0 || aln->a[i].qn != aln->a[i-1].qn) ol->list[aln->a[i].qn].align_length = i;
    }

	return;
}


void ul_lalign(overlap_region_alloc* ol, Candidates_list *cl, const ul_idx_t *uref, char *qstr, 
                        uint64_t ql, UC_Read* qu, UC_Read* tu, Correct_dumy* dumy, bit_extz_t *exz,
                        haplotype_evdience_alloc* hap, kvec_t_u64_warp* v_idx,   
                        double e_rate, int64_t wl, kv_ul_ov_t *aln, int64_t sid, void *km)
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

    if(!aln) {    
        resize_UC_Read(qu, ql); qu->length = ql; memcpy(qu->seq, qstr, ql);
        for (i = 0; i < ol->length; i++) {
            z = &(ol->list[i]); ovl = z->x_pos_e+1-z->x_pos_s;
            if(!align_ul_ed_post_extz(z, uref, NULL, qu->seq, tu->seq, exz, err, w.window_length, -1, km)) {
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
        ul_gap_filling(ol, aln, wl, uref, NULL, NULL, qu->seq, tu, exz, err, ql, sid);
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
                gen_backtrace_adv_exz(&(z->w_list.a[k]), z, NULL, NULL, uref, qu->seq, tu->seq, exz, z->y_pos_strand, z->y_id);
            }
            // verify_aln(sid, z, qu, tu, NULL, NULL, uref);


            ol->mapped_overlaps_length += ovl; 
            append_unmatched_wins(z, w.window_length);
            calculate_ul_boundary_cigars(z, uref, dumy, qu, err, w.window_length);
        }
        partition_ul_overlaps_advance(ol, uref, qu, tu, dumy, hap, 1, err, w.window_length, km);
    }
}