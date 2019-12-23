#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include "Correct.h"
#include "Levenshtein_distance.h"
#include "edlib.h"
#include "Assembly.h"
#include "CommandLines.h"
#include "ksw2.h"

long long T_total_match=0;
long long T_total_unmatch=0;
long long T_total_mis=0;
pthread_mutex_t debug_statistics ;




void align(const char *tseq, const char *qseq, const int tl, const int ql, 
const uint8_t *c, int sc_mch, int sc_mis, int gapo, int gape, int bandLen, int zdrop,
int* max_q_pos, int* max_t, int* score)
{
	int i, a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
	int8_t mat[25] = {a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0};
	uint8_t *ts, *qs;
	ksw_extz_t ez;

	memset(&ez, 0, sizeof(ksw_extz_t));
	ts = (uint8_t*)malloc(tl);
	qs = (uint8_t*)malloc(ql);
	for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
	for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
	///ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);
    ksw_extz2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, bandLen, zdrop, sc_mch, 0, &ez);
    /**
	for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
		printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
	putchar('\n');
    **/
	free(ez.cigar); free(ts); free(qs);
    (*score) = ez.max;
}


void afine_gap_alignment(const char *tseq, uint8_t* tnum, const int tl, 
const char *qseq, uint8_t* qnum, const int ql, const uint8_t *c2n, const int strand,
int sc_mch, int sc_mis, int gapo, int gape, int bandLen, int zdrop, int end_bonus,
long long* max_t_pos, long long* max_q_pos, long long* score, long long* droped)
{
    (*max_t_pos) = (*max_q_pos) = -1;
	int i, a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
	int8_t mat[25] = {a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0};
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

    ksw_extz2_sse(0, ql, qnum, tl, tnum, 5, mat, gapo, gape, bandLen, zdrop, end_bonus, 0, &ez);

    (*score) = ez.max;
    (*max_t_pos) = ez.max_t;
    (*max_q_pos) = ez.max_q;
    (*droped) = ez.zdropped;
    
    /**
	for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
		printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
	putchar('\n');
    **/
	free(ez.cigar);
}


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

void destory_Round2_alignment(Round2_alignment* h)
{
    destory_Correct_dumy(&(h->dumy));
    destory_Cigar_record(&(h->cigar));
    destory_Cigar_record(&(h->tmp_cigar));
}



///y_length > x_length
unsigned int edit_distance_normal_test_banded(char* y, int y_length, char* x, int x_length, int error_cut, int matrix[1000][1000] )
{   memset(matrix, 0, sizeof(matrix));

    int i, j;
    for (i = 0; i <= x_length; i++)
    {
        matrix[i][0] = i;
    }

    int digonal, up, left;
    unsigned int min;

    ///一列列算的
    for (i = 0; i < x_length; i++)
    {
        for (j = 0; j < y_length; j++)
        {
            ///matrix[i + 1][j + 1]
            digonal = matrix[i][j] + (x[i] !=  y[j]);
            up = matrix[i + 1][j] + 1;
            left = matrix[i][j + 1] + 1;
            min = digonal;
            if (up < min)
            {
                min = up;
            }

            if (left< min)
            {
                min = left;
            }

            matrix[i + 1][j + 1] = min;
        }
    }

    min = (unsigned int)-1;
    for (j = x_length; j <= y_length; j++)
    {
        if (matrix[i][j] < min)
        {
            min = matrix[i][j];
        }
    }


    return min <= error_cut?min:(unsigned int)(-1);
}

void verify_get_interval(long long window_start, long long window_end, overlap_region_alloc* overlap_list,Correct_dumy* dumy)
{
    long long i;
    long long match_length = 0;
    long long match_lengthNT = 0;
    long long Len;

    for (i = 0; i < overlap_list->length; i++)
    {
        if((Len = OVERLAP(window_start, window_end, overlap_list->list[i].x_pos_s, overlap_list->list[i].x_pos_e)) > 0)
        {
            if (Len == WINDOW)
            {
                match_length++;
                long long j;
                for (j = 0; j < dumy->length; j++)
                {
                    if (i==dumy->overlapID[j])
                    {
                        break;
                    }
                }

                if (j >= dumy->length)
                {
                    fprintf(stderr, "+ERROR interval\n");
                
                    fprintf(stderr, "i: %u, window_start: %u, window_end: %u, x_pos_s: %u, x_pos_e: %u\n", 
                    i, window_start, window_end, overlap_list->list[i].x_pos_s, overlap_list->list[i].x_pos_e);
                }
            }
            else
            {
                match_lengthNT++;
                long long j;
                for (j = 0; j < dumy->lengthNT; j++)
                {
                    if (i==dumy->overlapID[dumy->size - j - 1])
                    {
                        break;
                    }
                }

                if (j >= dumy->lengthNT)
                {
                    fprintf(stderr, "-ERROR interval\n");
                
                    fprintf(stderr, "i: %u, window_start: %u, window_end: %u, x_pos_s: %u, x_pos_e: %u, dumy->lengthNT: %u\n", 
                    i, window_start, window_end, overlap_list->list[i].x_pos_s, overlap_list->list[i].x_pos_e, dumy->lengthNT);
                }
            }
            
            

           
        
        }
    }

    if (match_length != dumy->length || match_lengthNT != dumy->lengthNT)
    {
        fprintf(stderr, "****************ERROR interval length*******************\n");
        fprintf(stderr, "match_length: %u\n", match_length);
        fprintf(stderr, "dumy->length: %u\n", dumy->length);
        fprintf(stderr, "window_start: %u, window_end: %u\n", window_start, window_end);
    }
    
}


inline int get_interval_back(long long window_start, long long window_end, overlap_region_alloc* overlap_list, Correct_dumy* dumy)
{
    long long i;
    int flag = 0;

    for (i = dumy->start_i; i < overlap_list->length; i++)
    {
        ///只会发生在这个interval比list里所有元素都小的情况
        ///这种情况下一个interval需要从0开始
        if (window_start < overlap_list->list[i].x_pos_s)
        {
            dumy->start_i = 0;
            dumy->length = 0;
            return -1;
        }
        else if(window_start >= overlap_list->list[i].x_pos_s && window_start <= overlap_list->list[i].x_pos_e)
        {
            dumy->start_i = i;
            break;
        }
    }

    ///只会发生在这个window比list里所有元素都大的情况
    ///这种情况下一个window也无需遍历了
    if (i >= overlap_list->length)
    {
        dumy->start_i = overlap_list->length;
        dumy->length = 0;
        return -2;
    }
    
    ///走到这里的时候，至少window_start的要求是满足了
    dumy->length = 0;

    for (; i < overlap_list->length; i++)
    {
        if(overlap_list->list[i].x_pos_s <= window_start && overlap_list->list[i].x_pos_e >= window_end)
        {
            dumy->overlapID[dumy->length] = i;
            dumy->length++;
        }
        else if(overlap_list->list[i].x_pos_s > window_start)
        {
            break;
        }
    }

    if ( dumy->length == 0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}


inline int get_interval(long long window_start, long long window_end, overlap_region_alloc* overlap_list, Correct_dumy* dumy)
{
    long long i;
    int flag = 0;
    long long Begin, End, Len;

    // fprintf(stderr, "overlap_list->length: %d, dumy->size: %d, dumy->start_i: %d\n", 
    // overlap_list->length, dumy->size, dumy->start_i);
    // fflush(stderr);

    for (i = dumy->start_i; i < overlap_list->length; i++)
    {
        ///只会发生在这个interval比list里所有元素都小的情况
        ///这种情况下一个interval需要从0开始
        if (window_end < overlap_list->list[i].x_pos_s)
        {
            dumy->start_i = 0;
            dumy->length = 0;
            dumy->lengthNT = 0;
            return 0;
        }
        else ///只要window_end >= overlap_list->list[i].x_pos_s，就有可能重叠
        {
            dumy->start_i = i;
            break;
        }
    }

    

    ///只会发生在这个window比list里所有元素都大的情况
    ///这种情况下一个window也无需遍历了
    if (i >= overlap_list->length)
    {
        dumy->start_i = overlap_list->length;
        dumy->length = 0;
        dumy->lengthNT = 0;
        return -2;
    }
    
    dumy->length = 0;
    dumy->lengthNT = 0;

    for (; i < overlap_list->length; i++)
    {
        // fprintf(stderr, "inner i: %d, x_pos_s: %d, x_pos_e: %d, window_start: %d, window_end: %d\n",
        // i, overlap_list->list[i].x_pos_s, overlap_list->list[i].x_pos_e,
        // window_start, window_end);
        // fflush(stderr);

        // fprintf(stderr, "dumy->length: %d, dumy->lengthNT: %d, dumy->size: %d\n",
        // dumy->length, dumy->lengthNT, dumy->size);
        // fflush(stderr);
        
        if((Len = OVERLAP(window_start, window_end, overlap_list->list[i].x_pos_s, overlap_list->list[i].x_pos_e)) > 0)
        {
            ///sometimes the length of window > WINDOW, but overlap length == WINDOW
            if (Len == WINDOW && window_end - window_start + 1 == WINDOW)
            {
                dumy->overlapID[dumy->length] = i;
                dumy->length++; 
            }
            else
            {
                dumy->lengthNT++;
                dumy->overlapID[dumy->size - dumy->lengthNT] = i;
            }
        }
        
        if(overlap_list->list[i].x_pos_s > window_end)
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
    long long i;
    int flag = 0;
    long long Begin, End, Len;
    long long overlap_length;


    for (i = dumy->start_i; i < overlap_list->length; i++)
    {
        ///只会发生在这个interval比list里所有元素都小的情况
        ///这种情况下一个interval需要从0开始
        if (window_end < overlap_list->list[i].x_pos_s)
        {
            dumy->start_i = 0;
            dumy->length = 0;
            dumy->lengthNT = 0;
            return 0;
        }
        else ///只要window_end >= overlap_list->list[i].x_pos_s，就有可能重叠
        {
            dumy->start_i = i;
            break;
        }
    }

    ///只会发生在这个window比list里所有元素都大的情况
    ///这种情况下一个window也无需遍历了
    if (i >= overlap_list->length)
    {
        dumy->start_i = overlap_list->length;
        dumy->length = 0;
        dumy->lengthNT = 0;
        return -2;
    }
    
    dumy->length = 0;
    dumy->lengthNT = 0;


    
    long long fake_length = 0;

    for (; i < overlap_list->length; i++)
    {
        ///是否重叠   
        if((Len = OVERLAP(window_start, window_end, overlap_list->list[i].x_pos_s, overlap_list->list[i].x_pos_e)) > 0)
        {
            ///重叠数量
            fake_length++;

            ///重叠是否有效
            overlap_length = overlap_list->list[i].x_pos_e - overlap_list->list[i].x_pos_s + 1;
            if (overlap_list->list[i].is_match == 1)
            {
                dumy->overlapID[dumy->length] = i;
                dumy->length++; 
            }
        }
        
        if(overlap_list->list[i].x_pos_s > window_end)
        {
            break;
        }
    }

    ///fake_length是重叠的数量，而不是有效重叠的数量
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
    for (size_t i = 0; i < l; i++)
    {
        fprintf(stderr, "%c", s[i]);
    }

    fprintf(stderr, "\n");
    
}

void test_edit_distance_by_edlib(char* x_string, char* y_string, long long x_len, 
long long o_len, int threashold, int error, long long* total_mis)
{
    EdlibAlignResult result = edlibAlign(x_string, x_len, y_string, o_len, 
    edlibNewAlignConfig(threashold, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

    if (result.status == EDLIB_STATUS_OK) {
        if (result.editDistance != error)
        {   
            
            (*total_mis)++;


            if ((int)error != -1 && result.editDistance==-1)
            {
                fprintf(stderr, "ERROR1\n");
            }

            if ((int)error < result.editDistance==-1 && 
            (int)error != -1 && result.editDistance!=-1)
            {
                fprintf(stderr, "ERROR2\n");
            }

            int up_length = o_len - result.endLocations[0] - 1;
            int left_length = result.endLocations[0] - x_len;

            
            char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
            int cigar_length = strlen(cigar);

            int i = cigar_length - 1;
            int j = 0;
            char tmp;
           

            while (i >= 0)
            {
                if (up_length < 0 || left_length < 0)
                {
                    break;
                }

                if (cigar[i] == 'I' || cigar[i] == 'D')
                {
                    tmp = cigar[i];
                    cigar[i] = '\0';
                    j = i - 1;
                    while (cigar[j] <= '9' && cigar[j] >= '0' && j >= 0)
                    {
                        j--;
                    }
                    j++;
                    int Len = atoi(cigar + j);
                    cigar[i] = tmp;
                    i = j - 1;

                    if (tmp == 'I')
                    {
                        up_length = up_length - Len;
                        left_length = left_length + Len;
                    }
                    else
                    {
                        left_length = left_length - Len;
                        up_length = up_length + Len;
                    }
                }
                else
                {
                    i--;
                }
            }

            
            if (up_length >= 0 && left_length >= 0)
            {
                fprintf(stderr, "****\nedlib: %d, alignmentLength: %d, startLocations: %d, endLocations: %d\n", 
                result.editDistance, result.alignmentLength, result.startLocations[0], result.endLocations[0]);
                fprintf(stderr,"%s\n", cigar);
                print_string(x_string, x_len);
                print_string(y_string, o_len);
                fprintf(stderr, "BPM: %d\n", error);
            }

            free(cigar);
        }
    }
    edlibFreeAlignResult(result);
}


void fill_subregion(char* r, long long start_pos, long long length, uint8_t strand, All_reads* R_INF, long long ID, 
int extra_begin, int extra_end)
{
    
    recover_UC_Read_sub_region(r+extra_begin, start_pos, length, strand, R_INF, ID);
    memset(r, 'N', extra_begin);
    memset(r+extra_begin+length, 'N', extra_end);
}

int determine_overlap_region(int threshold, long long y_start, long long y_ID, long long Window_Len, All_reads* R_INF,
int* r_extra_begin, int* r_extra_end, long long* r_y_start, long long* r_y_length)
{
    int extra_begin;
    int extra_end;
    long long currentIDLen;
    long long o_len;

    ///the length of y
    currentIDLen = Get_READ_LENGTH((*R_INF), y_ID);
    
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
    int end_site;
    unsigned int error;

    

    x_len = x_end - x_start + 1;
    threshold = x_len * THRESHOLD_RATE;
    /****************************may have bugs********************************/
    threshold = Adjust_Threshold(threshold, x_len);
    /****************************may have bugs********************************/
        
    y_start = (x_start - overlap_x_s) + overlap_y_s;
    
    Window_Len = x_len + (threshold << 1);

    if(!determine_overlap_region(threshold, y_start, y_id, Window_Len, R_INF,
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

    end_site = Reserve_Banded_BPM(y_string, Window_Len, x_string, x_len, threshold, &error);

    // if(y_id == 439960 && (x_id == 439950 || x_id == 5598234))
    // {
    //     fprintf(stderr, "y_id: %d, x_id: %d, error: %llu", y_id, x_id, error);
    // }

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
    long long currentID, currentIDLen;
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

    ///这些是整个window被完全覆盖的
    for (i = 0; i < dumy->length; i++)
    {
        extra_begin = extra_end = 0;
        ///整个window被覆盖的话，read本身上的区间就是[window_start, window_end]
        x_len = WINDOW;
        currentID = dumy->overlapID[i];
        x_start = window_start;
        ///y上的相对位置
        y_start = (x_start - overlap_list->list[currentID].x_pos_s) + overlap_list->list[currentID].y_pos_s;
        /****************************may have bugs********************************/
        y_start += y_start_offset(x_start, &(overlap_list->list[currentID].f_cigar));
        /****************************may have bugs********************************/
        
  
        if(!determine_overlap_region(THRESHOLD, y_start, overlap_list->list[currentID].y_id, Window_Len, R_INF,
        &extra_begin, &extra_end, &y_start, &o_len))
        {
            append_window_list(&overlap_list->list[currentID], window_start, window_end, 
            -1, -1, -1, -1, -1, -1);
            continue;
        }

        fill_subregion(dumy->overlap_region_group[groupLen], y_start, o_len, overlap_list->list[currentID].y_pos_strand, 
        R_INF, overlap_list->list[currentID].y_id, extra_begin, extra_end);

        // fprintf(stderr, "!i: %d\n", i);
        // fflush(stderr);

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

            
            if (return_sites_error[0]!=(unsigned int)-1)
            {
                overlap_list->list[overlapID[0]].align_length += x_len;

                append_window_list(&overlap_list->list[overlapID[0]], window_start, window_end, 
                            y_startGroup[0], y_startGroup[0] + return_sites[0], (int)return_sites_error[0],
                            y_extra_begin[0], y_extra_end[0], error_threshold[0]);
            }
            else
            {
                append_window_list(&overlap_list->list[overlapID[0]], window_start, window_end, y_startGroup[0], -1, -1,
                y_extra_begin[0], y_extra_end[0], error_threshold[0]);
            }
            

            if (return_sites_error[1]!=(unsigned int)-1)
            {
                overlap_list->list[overlapID[1]].align_length += x_len;

                append_window_list(&overlap_list->list[overlapID[1]], window_start, window_end, 
                            y_startGroup[1], y_startGroup[1] + return_sites[1], (int)return_sites_error[1],
                            y_extra_begin[1], y_extra_end[1], error_threshold[1]);
            }
            else
            {
                append_window_list(&overlap_list->list[overlapID[1]], window_start, window_end, y_startGroup[1], -1, -1,
                y_extra_begin[1], y_extra_end[1], error_threshold[1]);
            }
            

            if (return_sites_error[2]!=(unsigned int)-1)
            {
                overlap_list->list[overlapID[2]].align_length += x_len;

                append_window_list(&overlap_list->list[overlapID[2]], window_start, window_end, 
                            y_startGroup[2], y_startGroup[2] + return_sites[2], (int)return_sites_error[2],
                            y_extra_begin[2], y_extra_end[2], error_threshold[2]);
            }
            else
            {
                append_window_list(&overlap_list->list[overlapID[2]], window_start, window_end, y_startGroup[2], -1, -1,
                y_extra_begin[2], y_extra_end[2], error_threshold[2]);
            }
            

            if (return_sites_error[3]!=(unsigned int)-1)
            {
                overlap_list->list[overlapID[3]].align_length += x_len;

                append_window_list(&overlap_list->list[overlapID[3]], window_start, window_end, 
                            y_startGroup[3], y_startGroup[3] + return_sites[3], (int)return_sites_error[3],
                            y_extra_begin[3], y_extra_end[3], error_threshold[3]);
            }
            else
            {
                append_window_list(&overlap_list->list[overlapID[3]], window_start, window_end, y_startGroup[3], -1, -1,
                y_extra_begin[3], y_extra_end[3], error_threshold[3]);
            }    
        }
    }

    // fprintf(stderr, "(1) dumy->size: %d\n", dumy->size);
    // fflush(stderr);
    

    if (groupLen == 1)
    {
        end_site = Reserve_Banded_BPM(dumy->overlap_region_group[0], Window_Len, x_string, WINDOW, THRESHOLD, &error);

        

        if (error!=(unsigned int)-1)
        {
            overlap_list->list[overlapID[0]].align_length += x_len;

            append_window_list(&overlap_list->list[overlapID[0]], window_start, window_end, 
                                y_startGroup[0], y_startGroup[0] + end_site, (int)error,
                                y_extra_begin[0], y_extra_end[0], error_threshold[0]);
        }
        else
        {
            append_window_list(&overlap_list->list[overlapID[0]], window_start, window_end, y_startGroup[0], -1, -1,
            y_extra_begin[0], y_extra_end[0], error_threshold[0]);
        }
    }
    else if (groupLen > 1)
    {
        Reserve_Banded_BPM_4_SSE_only(dumy->overlap_region_group[0], dumy->overlap_region_group[1], 
                dumy->overlap_region_group[2], dumy->overlap_region_group[3], Window_Len, x_string, WINDOW,
			        return_sites, return_sites_error, THRESHOLD, dumy->Peq_SSE);

        for (i = 0; i < groupLen; i++)
        {
            if (return_sites_error[i]!=(unsigned int)-1)
            {
                overlap_list->list[overlapID[i]].align_length += x_len;
                append_window_list(&overlap_list->list[overlapID[i]], window_start, window_end, 
                                y_startGroup[i], y_startGroup[i] + return_sites[i], (int)return_sites_error[i],
                                y_extra_begin[i], y_extra_end[i], error_threshold[i]);
            }
            else
            {
                append_window_list(&overlap_list->list[overlapID[i]], window_start, window_end, y_startGroup[i], -1, -1,
                y_extra_begin[i], y_extra_end[i], error_threshold[i]);
            }
            
        }

        groupLen = 0;
    }


    // fprintf(stderr, "(2) dumy->size: %d\n", dumy->size);
    // fflush(stderr);
    

    long long reverse_i = dumy->size - 1;
    int threshold;
    
    ///这些是整个window被部分覆盖的
    for (i = 0; i < dumy->lengthNT; i++)
    {
        
        extra_begin = extra_end = 0;
        currentID = dumy->overlapID[reverse_i--];
        x_start = MAX(window_start, overlap_list->list[currentID].x_pos_s);
        x_end = MIN(window_end, overlap_list->list[currentID].x_pos_e);


        ///这个是和当前窗口重叠的长度
        x_len = x_end - x_start + 1;
        threshold = x_len * THRESHOLD_RATE;
        /****************************may have bugs********************************/
        threshold = Adjust_Threshold(threshold, x_len);
        /****************************may have bugs********************************/
        

        ///y上的相对位置
        y_start = (x_start - overlap_list->list[currentID].x_pos_s) + overlap_list->list[currentID].y_pos_s;
        /****************************may have bugs********************************/
        y_start += y_start_offset(x_start, &(overlap_list->list[currentID].f_cigar));
        /****************************may have bugs********************************/


        



        // fprintf(stderr, "lengthNT: %d, i: %d, window_start: %d, window_end: %d, x_start: %d, x_end: %d\n", 
        // dumy->lengthNT, i, window_start, window_end, x_start, x_end);
        // fflush(stderr);

        // fprintf(stderr, "x_pos_s: %d, x_pos_e: %d, y_pos_s: %d, y_pos_e: %d\n", 
        // overlap_list->list[currentID].x_pos_s, 
        // overlap_list->list[currentID].x_pos_e,
        // overlap_list->list[currentID].y_pos_s,
        // overlap_list->list[currentID].y_pos_e);

        // fprintf(stderr, "x_id: %d, x_length: %d, y_id: %d, y_length: %d, y_pos_strand: %d\n",
        // overlap_list->list[currentID].x_id, 
        // Get_READ_LENGTH((*R_INF), overlap_list->list[currentID].x_id),
        // overlap_list->list[currentID].y_id,
        // Get_READ_LENGTH((*R_INF), overlap_list->list[currentID].y_id),
        // overlap_list->list[currentID].y_pos_strand);
        // fflush(stderr);

        // fprintf(stderr, "y_start: %d, y_start_offset: %d\n", y_start, 
        // y_start_offset(x_start, &(overlap_list->list[currentID].f_cigar)));
        // print_fake_gap(&(overlap_list->list[currentID].f_cigar));
        // fflush(stderr);
        
        // fprintf(stderr, "(31) i: %d, dumy->size: %d, y_start: %d, yLen: %d\n", i, dumy->size, y_start,
        // Get_READ_LENGTH((*R_INF), overlap_list->list[currentID].y_id));
        // fflush(stderr);
        
        Window_Len = x_len + (threshold << 1);

        if(!determine_overlap_region(threshold, y_start, overlap_list->list[currentID].y_id, Window_Len, R_INF,
        &extra_begin, &extra_end, &y_start, &o_len))
        {
            append_window_list(&overlap_list->list[currentID], x_start, x_end, 
            -1, -1, -1, -1, -1, -1);
            continue;
        }

        /**
        if(overlap_list->list[currentID].x_id == 18390 
        && overlap_list->list[currentID].y_id == 18419)
        {
            fprintf(stderr, "x_start: %d, x_len: %d, y_start: %d, y_offset: %d, extra_begin: %d, extra_end: %d, o_len: %d\n", 
                x_start, x_len, y_start, 
                y_start_offset(x_start, &(overlap_list->list[currentID].f_cigar)),
                extra_begin, extra_end, o_len);
            
        }
        **/


        // fprintf(stderr, "(32) i: %d, dumy->size: %d, y_start: %d, o_len: %d, extra_begin: %d, extra_end: %d\n", 
        // i, dumy->size, y_start, o_len, extra_begin, extra_end);
        // if(o_len == -13)
        // {
        //     print_fake_gap(&(overlap_list->list[currentID].f_cigar));
        // }
        // fflush(stderr);

        
        fill_subregion(dumy->overlap_region, y_start, o_len, overlap_list->list[currentID].y_pos_strand, 
        R_INF, overlap_list->list[currentID].y_id, extra_begin, extra_end);

        

        // fprintf(stderr, "(333332) i: %d, dumy->size: %d, y_start: %d, yLen: %d, threshold: %d, Window_Len: %d\n", i, dumy->size, y_start,
        // Get_READ_LENGTH((*R_INF), overlap_list->list[currentID].y_id), threshold, Window_Len);
        // fflush(stderr);
        

        x_string = r_string + x_start;
        y_string = dumy->overlap_region;


        end_site = Reserve_Banded_BPM(y_string, Window_Len, x_string, x_len, threshold, &error);


        // fprintf(stderr, "(33) dumy->size: %d\n", dumy->size);
        // fflush(stderr);

        if (error!=(unsigned int)-1)
        {
            overlap_list->list[currentID].align_length += x_len;
            append_window_list(&overlap_list->list[currentID], x_start, x_end, y_start, y_start + end_site, (int)error,
            extra_begin, extra_end, threshold);
        }
        else
        {
            append_window_list(&overlap_list->list[currentID], x_start, x_end, y_start, -1, -1,
            extra_begin, extra_end, threshold);
        }
    }

    // fprintf(stderr, "(3) dumy->size: %d\n", dumy->size);
    // fflush(stderr);


    // fprintf(stderr, "************groupLen: %d\n", groupLen);
    // fflush(stderr);
}

void debug_stats(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read, 
                        long long* matched_overlap_0, long long* matched_overlap_1)
{
    long long j;
    long long Len_x;
    int threshold;
    long long y_start;
    long long Len_y;
    long long currentIDLen;
    for (j = 0; j < overlap_list->length; j++)
    {
        Len_x = overlap_list->list[j].x_pos_e -  overlap_list->list[j].x_pos_s + 1;

        if (overlap_list->list[j].is_match == 1)
        {
            if (overlap_list->list[j].y_pos_strand == 0)
            {
                (*matched_overlap_0)++;
                ///(*matched_overlap_0) = (*matched_overlap_0) + overlap_list->list[j].align_length;
                ///(*matched_overlap_0) = (*matched_overlap_0) + Len_x;
            }
            else
            {
                (*matched_overlap_1)++;
                ///(*matched_overlap_1) = (*matched_overlap_1) + overlap_list->list[j].align_length;
                ///(*matched_overlap_1) = (*matched_overlap_1) + Len_x;
            }
        }
    }

    


}

inline double trim_error_rate(overlap_region_alloc* overlap_list, long long ID)
{
    long long tLen, tError,i, subWinLen, subWinNum;
    
    tLen = 0;
    tError = 0;

    subWinNum = overlap_list->list[ID].w_list_length;

    if(subWinNum < 5)
    {
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
                ///tError += (Adjust_Threshold(subWinLen*THRESHOLD_RATE, subWinLen) * 2);
                tError += (Adjust_Threshold(subWinLen*THRESHOLD_RATE, subWinLen) * 3);
            }
        }
    }
    else
    {
        for (i = 1; i < subWinNum - 1; i++)
        {
            subWinLen = overlap_list->list[ID].w_list[i].x_end - overlap_list->list[ID].w_list[i].x_start + 1;
            tLen += subWinLen;

            if(overlap_list->list[ID].w_list[i].y_end != -1)
            {
                tError += overlap_list->list[ID].w_list[i].error;
            }
            else
            {
                ///tError += (Adjust_Threshold(subWinLen*THRESHOLD_RATE, subWinLen) * 2);
                tError += (Adjust_Threshold(subWinLen*THRESHOLD_RATE, subWinLen) * 3);
            }
        }
    }
    

    

    double error_rate = (double)(tError)/(double)(tLen);

    return error_rate;
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
        if(determine_overlap_region(SubThreshold, y_start, y_id, SubWindowLen, R_INF,
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
                    
        ///error等于-1说明没匹配
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
    R_INF, &extra_begin, &extra_end, &y_beg, &o_len))
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

inline double non_trim_error_rate(overlap_region_alloc* overlap_list, long long ID,
All_reads* R_INF, Correct_dumy* dumy, UC_Read* g_read)
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
            int threshold = double_error_threshold(overlap_list->list[ID].w_list[i].error_threshold, x_len);
            int Window_Len = x_len + (threshold << 1);
            unsigned int r_error_left;
            int r_x_end_left, r_y_end_left, aligned_xLen_left;
            unsigned int r_error_right;
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
            if(i < overlap_list->list[ID].w_list_length - 1 && overlap_list->list[ID].w_list[i + 1].y_end != -1)
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
                verify_sub_window(R_INF, dumy, g_read, overlap_list->list[ID].w_list[i].x_start, 
                x_len, y_beg_left, Window_Len, overlap_list->list[ID].y_id, overlap_list->list[ID].y_pos_strand, 
                threshold, 0, &r_error_left, &r_y_end_left, &r_x_end_left, &aligned_xLen_left);
            }

            if(y_beg_right != -1)
            {
                verify_sub_window(R_INF, dumy, g_read, overlap_list->list[ID].w_list[i].x_start, 
                x_len, y_beg_right, Window_Len, overlap_list->list[ID].y_id, overlap_list->list[ID].y_pos_strand,
                threshold, 1, &r_error_right, &r_y_end_right, &r_x_end_right, &aligned_xLen_right);
            }

            /**
           if(i == 0 || i == subWinNum - 1)
           {
               if(((aligned_xLen_left + aligned_xLen_right) < x_len)
                 &&
                 (x_len - (aligned_xLen_left + aligned_xLen_right) > 20))
                {
                    return 1.0;
                }
           }
           **/


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



void mark_duplicate(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read)
{
    long long j, overlapLen;
    double rate;


    for (j = 0; j < overlap_list->length; j++)
    {
        if(overlap_list->list[j].is_match == 1)
        {
            rate = trim_error_rate(overlap_list, j);

            if(rate > 0.01)
            {
                overlapLen = overlap_list->list[j].x_pos_e - overlap_list->list[j].x_pos_s + 1;
                overlap_list->mapped_overlaps_length -= overlapLen;
                overlap_list->list[j].is_match = 0;
            }
        }
    }


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

    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2是x缺字符（y多字符），而3是y缺字符（x多字符）
    ///while (x_i < x_len && y_i < y_len && cigar_i < cigar->length)
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
        }///2是x缺字符（y多字符）
        else if (operation == 2)
        {

            if(if_is_homopolymer_repeat(x_i, x, x_len) || if_is_homopolymer_repeat(y_i, y, y_len))
            {
                hpm_error++;
            }/**
            else
            {
                if(x_i - 5 > 0 && x_i + 5 <= x_len
                &&
                y_i - 5 > 0 && y_i + 5 <= y_len)
                {
                    fprintf(stderr, "x: %.*s\ny: %.*s\n\n", 10, x + x_i - 5, 10, y + y_i - 5);
                }
                
            }
            **/
            

            cigar_error += operationLen;
            y_i += operationLen;
        }///3是y缺字符（x多字符）
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



void count_no_HPM_errors(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read, 
                        long long* total_errors, long long* total_hpm_errors)
{


    long long j, i;
    long long y_id, y_strand, y_readLen;
    long long x_start, x_end, x_len, y_start, y_end, y_len, error;
    char* x_string;
    char* y_string;
    CIGAR* cigar;
    int hpm_error;
    (*total_errors) = 0;
    (*total_hpm_errors) = 0;

    for (j = 0; j < overlap_list->length; j++)
    {
        y_id = overlap_list->list[j].y_id;
        y_strand = overlap_list->list[j].y_pos_strand;
        y_readLen = Get_READ_LENGTH((*R_INF), y_id);

        if (overlap_list->list[j].is_match == 1)
        {
            ///for (i = 0; i < overlap_list->list[j].w_list_length; i++)
            for (i = 1; i < overlap_list->list[j].w_list_length - 1; i++)
            {
                if(overlap_list->list[j].w_list[i].y_end != -1)
                {
                    x_start = overlap_list->list[j].w_list[i].x_start;
                    x_end = overlap_list->list[j].w_list[i].x_end;
                    x_len = x_end - x_start + 1;

                    x_string = g_read->seq + x_start;

                    y_start = overlap_list->list[j].w_list[i].y_start;
                    y_end = overlap_list->list[j].w_list[i].y_end;
                    y_len = y_end - y_start + 1;

                    recover_UC_Read_sub_region(dumy->overlap_region, y_start, y_len, y_strand, R_INF, y_id);
                    y_string = dumy->overlap_region;


                    cigar = &overlap_list->list[j].w_list[i].cigar;
                    error = overlap_list->list[j].w_list[i].error;

                    hpm_error = calculate_hpm_errors(x_string, x_len, y_string, y_len, cigar, error);

                    ///fprintf(stderr, "hpm_error: %d, error: %d\n", hpm_error, error);
                    (*total_errors) += error;
                    (*total_hpm_errors) += hpm_error;
                }
            }

        }
    }
    


}


void debug_output_overlaps(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read, 
                        long long* matched_overlap_0, long long* matched_overlap_1)
{
    long long j;
    long long Len_x;
    int threshold;
    long long y_start;
    long long Len_y;
    long long currentIDLen = 0;
    fprintf(stderr, "overlap_list->length: %d\n", overlap_list->length);
    for (j = 0; j < overlap_list->length; j++)
    {
        if(memcmp("m64013_190412_043951/108332093/ccs", Get_NAME((*R_INF),overlap_list->list[j].y_id), 
     Get_NAME_LENGTH((*R_INF), overlap_list->list[j].y_id)) == 0)
        {
            fprintf(stderr, "******************************x_id: %d, y_id: %d, y_name: %.*s, error_rate: %f, is_match: %d*******************************\n",
            overlap_list->list[j].x_id,
            overlap_list->list[j].y_id,
            Get_NAME_LENGTH((*R_INF), overlap_list->list[j].y_id), Get_NAME((*R_INF),overlap_list->list[j].y_id),
            trim_error_rate(overlap_list, j), overlap_list->list[j].is_match);

            print_fake_gap(&overlap_list->list[j].f_cigar);
        }

        Len_x = overlap_list->list[j].x_pos_e -  overlap_list->list[j].x_pos_s + 1;

        if (overlap_list->list[j].is_match == 1)
        {
            currentIDLen++;
            fprintf(stderr, "y_name: %.*s\n",
            Get_NAME_LENGTH((*R_INF), overlap_list->list[j].y_id), Get_NAME((*R_INF),overlap_list->list[j].y_id));
        }
        // else
        // {
        //     fprintf(stderr, "not match, y_name: %.*s\n",
        //     Get_NAME_LENGTH((*R_INF), overlap_list->list[j].y_id), Get_NAME((*R_INF),overlap_list->list[j].y_id));
        // }
        
    }

    fprintf(stderr, "currentIDLen: %d\n\n", currentIDLen);

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
    ///2是x缺字符（y多字符），而3是y缺字符（x多字符）
    ///while (x_i < x_len && y_i < y_len && cigar_i < cigar->length)
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
        }///2是x缺字符（y多字符）
        else if (operation == 2)
        {
            cigar_error += operationLen;
            y_i += operationLen;
        }///3是y缺字符（x多字符）
        else if (operation == 3)
        {
            cigar_error += operationLen;
            x_i += operationLen;
        }
        
        cigar_i++;
    }


    endloop: 
    ///return;
    
    if (cigar_error != error)
    {
        /**
        fprintf(stderr, "error cigar_error: cigar_error: %d, error: %d\n", cigar_error, error);
        for (i = 0; i < cigar->length; i++)
        {
            fprintf(stderr, "%u: %u\n", cigar->C_L[i], cigar->C_C[i]);
        }
        **/
        
       flag_error = 1;
        
    }
    
   
    if (flag_error == 1)
    {
        /**
        print_string(x, x_len);
        print_string(y, y_len);
        fprintf(stderr, "x_len: %d, y_len: %d, cigar_len: %d, error: %d\n", x_len, y_len, cigar->length, error);
        for (i = 0; i < cigar->length; i++)
        {
            fprintf(stderr, "%u: %u\n", cigar->C_L[i], cigar->C_C[i]);
        }
        **/
        
    }
    
    
    return flag_error;
    
}


int scan_cigar(CIGAR* cigar, int* get_error, int scanXLen, int direction)
{
    (*get_error) = -1;
    if(cigar->length == 1 && cigar->C_C[0] == 0)
    {
        (*get_error) = 0;
        return 1;
    }


    int x_i, y_i, cigar_i;
    x_i = 0;
    y_i = 0;
    
    int operation;
    int operationLen;
    int i;
    int cigar_error = 0;

    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2: there are more bases at y, 3: there are more bases at x
    if(direction == 0)
    {
        cigar_i = 0;
        while (cigar_i < cigar->length)
        {
            operation = cigar->C_C[cigar_i];
            operationLen = cigar->C_L[cigar_i];
            cigar_i++;

            if (operation == 0)
            {
                x_i = x_i + operationLen;
                y_i = y_i + operationLen;
                if(x_i >= scanXLen)
                {
                    (*get_error) = cigar_error;
                    return 1;
                }
            }
            else if (operation == 1)
            {
                for (i = 0; i < operationLen; i++)
                {
                    x_i++;
                    y_i++;
                    cigar_error++;
                    if(x_i >= scanXLen)
                    {
                        (*get_error) = cigar_error;
                        return 1;
                    }
                }
            }///2是x缺字符（y多字符）
            else if (operation == 2)
            {
                cigar_error += operationLen;
                y_i += operationLen;
            }///3是y缺字符（x多字符）
            else if (operation == 3)
            {
                for (i = 0; i < operationLen; i++)
                {
                    x_i++;
                    cigar_error++;
                    if(x_i >= scanXLen)
                    {
                        (*get_error) = cigar_error;
                        return 1;
                    }
                }
            }
        }
    }
    else
    {
        cigar_i = cigar->length - 1;
        while (cigar_i >= 0)
        {
            operation = cigar->C_C[cigar_i];
            operationLen = cigar->C_L[cigar_i];
            cigar_i--;
            if (operation == 0)
            {
                x_i = x_i + operationLen;
                y_i = y_i + operationLen;
                if(x_i >= scanXLen)
                {
                    (*get_error) = cigar_error;
                    return 1;
                }
            }
            else if (operation == 1)
            {
                for (i = 0; i < operationLen; i++)
                {
                    x_i++;
                    y_i++;
                    cigar_error++;
                    if(x_i >= scanXLen)
                    {
                        (*get_error) = cigar_error;
                        return 1;
                    }
                }
            }///2是x缺字符（y多字符）
            else if (operation == 2)
            {
                cigar_error += operationLen;
                y_i += operationLen;
            }///3是y缺字符（x多字符）
            else if (operation == 3)
            {
                for (i = 0; i < operationLen; i++)
                {
                    x_i++;
                    cigar_error++;
                    if(x_i >= scanXLen)
                    {
                        (*get_error) = cigar_error;
                        return 1;
                    }
                }
            }
        }
    }

    (*get_error) = cigar_error;
    return 0;
}

///[scanXbeg, scanXend]
int scan_cigar_interval(CIGAR* cigar, int* get_error, int scanXbeg, int scanXend)
{
    (*get_error) = -1;
    if(cigar->length == 1 && cigar->C_C[0] == 0)
    {
        (*get_error) = 0;
        return 1;
    }



    int x_i, y_i, cigar_i;
    x_i = 0;
    y_i = 0;
    
    int operation;
    int operationLen;
    int i;
    int cigar_error = 0;

    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2: there are more bases at y, 3: there are more bases at x
    cigar_i = 0;
    while (cigar_i < cigar->length)
    {
        operation = cigar->C_C[cigar_i];
        operationLen = cigar->C_L[cigar_i];
        cigar_i++;

        if (operation == 0)
        {
            for (i = 0; i < operationLen; i++)
            {
                if(x_i == scanXbeg)
                {
                    cigar_error = 0;
                }
                
                x_i++;
                y_i++;

                if(x_i == scanXend + 1)
                {
                    (*get_error) = cigar_error;
                    return 1;
                }
            }
        }
        else if (operation == 1)
        {
            for (i = 0; i < operationLen; i++)
            {
                if(x_i == scanXbeg)
                {
                    cigar_error = 0;
                }
                
                x_i++;
                y_i++;
                cigar_error++;

                if(x_i == scanXend + 1)
                {
                    (*get_error) = cigar_error;
                    return 1;
                }
            }
        }///2是x缺字符（y多字符）
        else if (operation == 2)
        {
            cigar_error += operationLen;
            y_i += operationLen;
        }///3是y缺字符（x多字符）
        else if (operation == 3)
        {
            for (i = 0; i < operationLen; i++)
            {
                if(x_i == scanXbeg)
                {
                    cigar_error = 0;
                }
                

                x_i++;
                cigar_error++;

                if(x_i == scanXend + 1)
                {
                    (*get_error) = cigar_error;
                    return 1;
                }
            }
        }
    }


    (*get_error) = cigar_error;
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
    

    if(oper == 3)
    {
        path_i++;
        y_i--;
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
    else if(oper == 2)
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

inline void generate_cigar(
	char* path, int path_length, window_list* result, int* start, int* end, unsigned int* old_error,
    char* x, int x_len, char* y)
{

    if ((*old_error) == 0)
    {
        result->cigar.C_L[0] = result->x_end - result->x_start + 1;
        result->cigar.C_C[0] = 0;
        result->cigar.length = 1;

        return;
    }
    

	int i = 0;
    result->cigar.length = 0;
    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    char pre_ciga = 5;
    int pre_ciga_length = 0;
    
    int terminate_site = -1;

    for (i = 0; i < path_length; i++)
    {
        if(path[i] == 1)
        {
            path[i] = 3;
            (*end)--;
            terminate_site = i;
        }
        else
        {
            break;
        }
    }

    for (i = path_length - 1; i >= 0; i--)
    {
        if(path[i] == 1)
        {
            path[i] = 3;
            (*start)++;
        }
        else
        {
            break;
        }
    }


    for (i = path_length - 1; i >= 0; i--)
    {

        if (pre_ciga != path[i])
        {
            if (pre_ciga_length != 0)
            {
                result->cigar.C_L[result->cigar.length] = pre_ciga_length;
                result->cigar.C_C[result->cigar.length] = pre_ciga;
                result->cigar.length++;
            }

            pre_ciga = path[i];
            pre_ciga_length = 1;
        }
        else
        {
            pre_ciga_length++;
        }
    }

    if (pre_ciga_length != 0)
    {
        result->cigar.C_L[result->cigar.length] = pre_ciga_length;
        result->cigar.C_C[result->cigar.length] = pre_ciga;
        result->cigar.length++;
    }

    ///verify_cigar(x, x_len, y + (*start), (*end) - (*start) + 1, &(result->cigar), error);
    
    
    y = y + (*start);

    int x_i, y_i;
    x_i = 0;
    y_i = 0;
    ///terminate_site = -1;
    for (i = path_length - 1; i > terminate_site; i--)
    {
        if(path[i] == 0)
        {
            x_i++;
            y_i++;
        }
        else if(path[i] == 1)
        {
            x_i++;
            y_i++;
        }
        else if(path[i] == 2)
        {
            move_gap_greedy(path, i, path_length, x, x_i, y, y_i, old_error);
            y_i++;
        }
        else if(path[i] == 3)
        {
            move_gap_greedy(path, i, path_length, x, x_i, y, y_i, old_error);
            x_i++;
        }
    }



    pre_ciga = 5;
    pre_ciga_length = 0;
    result->cigar.length = 0;
    for (i = path_length - 1; i >= 0; i--)
    {

        if (pre_ciga != path[i])
        {
            if (pre_ciga_length != 0)
            {
                result->cigar.C_L[result->cigar.length] = pre_ciga_length;
                result->cigar.C_C[result->cigar.length] = pre_ciga;
                result->cigar.length++;
            }

            pre_ciga = path[i];
            pre_ciga_length = 1;
        }
        else
        {
            pre_ciga_length++;
        }
    }

    if (pre_ciga_length != 0)
    {
        result->cigar.C_L[result->cigar.length] = pre_ciga_length;
        result->cigar.C_C[result->cigar.length] = pre_ciga;
        result->cigar.length++;
    }

    
    // if(verify_cigar(x, x_len, y, (*end) - (*start) + 1, &(result->cigar), *old_error))
    // {
    //     fprintf(stderr, "error\n");
    // }
    


    /**
    int x_i, y_i;
    x_i = 0;
    y_i = 0;
    int new_error = error;
    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2是x缺字符（y多字符），而3是y缺字符（x多字符）
    ///while (x_i < x_len && y_i < y_len && cigar_i < cigar->length)
    for (i = path_length - 1; i >= 0; i--)
    {
        if(path[i] == 0)
        {
            x_i++;
            y_i++;
        }
        else if(path[i] == 1)
        {
            x_i++;
            y_i++;
        }
        else if(path[i] == 2)
        {
            move_gap_greedy(path, i, path_length, x, x_i, y, y_i, &new_error);
            y_i++;
        }
        else if(path[i] == 3)
        {
            move_gap_greedy(path, i, path_length, x, x_i, y, y_i, &new_error);
            x_i++;
        }
    }






    
    CIGAR new_cigar;

    new_cigar.length = 0;
    pre_ciga = 5;
    pre_ciga_length = 0;
    

    for (i = path_length - 1; i >= 0; i--)
    {

        if (pre_ciga != path[i])
        {
            if (pre_ciga_length != 0)
            {
                new_cigar.C_L[new_cigar.length] = pre_ciga_length;
                new_cigar.C_C[new_cigar.length] = pre_ciga;
                new_cigar.length++;
            }

            pre_ciga = path[i];
            pre_ciga_length = 1;
        }
        else
        {
            pre_ciga_length++;
        }
    }

    if (pre_ciga_length != 0)
    {
        new_cigar.C_L[new_cigar.length] = pre_ciga_length;
        new_cigar.C_C[new_cigar.length] = pre_ciga;
        new_cigar.length++;
    }




    
    if(verify_cigar(x, x_len, y, (*end) - (*start) + 1, &new_cigar, new_error))
    {

        fprintf(stderr, "x_string: %.*s\n",  x_len, x);
                fprintf(stderr, "y_string: %.*s\n\n",  (*end) - (*start) + 1 , y);

        for (int j = 0; j < result->cigar.length; j++)
        {
            fprintf(stderr, "oper: %d, len: %d\n", result->cigar.C_C[j], result->cigar.C_L[j]);
        }

        for (int j = 0; j < new_cigar.length; j++)
        {
            fprintf(stderr, "new_cigar.oper: %d, new_cigar.len: %d\n", new_cigar.C_C[j], new_cigar.C_L[j]);
        }

    }
    **/
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
    ///2是x缺字符（y多字符），而3是y缺字符（x多字符）
    ///while (x_i < x_len && y_i < y_len && cigar_i < cigar->length)
    while (cigar_i < cigar->length)
    {
        operation = Get_Cigar_Type(cigar->record[cigar_i]);
        operationLen = Get_Cigar_Length(cigar->record[cigar_i]);

        if (operation == 0)
        {
            for (i = 0; i < operationLen; i++)
            {

                if (x[x_i]!=y[y_i])
                {
                    fprintf(stderr, "error match\n");
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
                    fprintf(stderr, "error mismatch, cigar_i: %d, x_i: %d, y_i: %d\n",cigar_i, x_i, y_i);
                    flag_error = 1;
                }

                if(Get_MisMatch_Base(cigar->lost_base[diff_i]) != y[y_i])
                {
                    fprintf(stderr, "mismatch x: %c, y: %c, mis[%d]: %c\n", x[x_i],y[y_i],diff_i, 
                    Get_MisMatch_Base(cigar->lost_base[diff_i]));
                }

                if(Get_Match_Base(cigar->lost_base[diff_i]) != x[x_i])
                {
                    fprintf(stderr, "match x: %c, y: %c, deletion[%d]: %c\n", 
                    x[x_i],y[y_i],diff_i, 
                    Get_Match_Base(cigar->lost_base[diff_i]));
                }
               
                


                x_i++;
                y_i++;
                diff_i++;
            }
        }///2是x缺字符（y多字符）
        else if (operation == 2)
        {
            cigar_error += operationLen;
            for (i = 0; i < operationLen; i++)
            {
                if(cigar->lost_base[diff_i] != y[y_i])
                {
                    fprintf(stderr, "insertion x: %c, y: %c, insertion[%d]: %c\n", x[x_i],y[y_i],diff_i, 
                    cigar->lost_base[diff_i]);
                }
                y_i++;
                diff_i++;
            }
        }///3是y缺字符（x多字符）
        else if (operation == 3)
        {
            cigar_error += operationLen;

            for (i = 0; i < operationLen; i++)
            {
                if(cigar->lost_base[diff_i] != x[x_i])
                {
                    fprintf(stderr, "deletion x: %c, y: %c, deletion[%d]: %c\n", x[x_i],y[y_i],diff_i, 
                    cigar->lost_base[diff_i]);
                }
                x_i++;
                diff_i++;
            }
            // x_i += operationLen;
            // diff_i += operationLen;
        }
        
        cigar_i++;
    }


    endloop: 
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
        fprintf(stderr, "x_len: %d, y_len: %d, cigar_len: %d, error: %d\n", x_len, y_len, cigar->length, error);
        for (i = 0; i < cigar->length; i++)
        {
            operation = Get_Cigar_Type(cigar->record[i]);
            operationLen = Get_Cigar_Length(cigar->record[i]);
            fprintf(stderr, "%u: %u\n", operationLen, operation);
        }
    }
    
    
    return flag_error;
    
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
        if(!determine_overlap_region(threshold, total_y_start, y_ID, Window_Len, R_INF,
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

        if(!determine_overlap_region(threshold, total_y_start, y_ID, Window_Len, R_INF,
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





inline void recalcate_window_back(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read)
{
    long long j, k, i;
    long long Len_x;
    int threshold;
    long long y_len;
    long long currentIDLen;
    long long matches;
    long long y_id;
    int y_strand;
    long long y_readLen;
    long long x_start;
    long long x_end;
    long long x_len;
    long long total_y_start;
    long long total_y_end;
    long long y_start;
    long long y_end;
    long long Window_Len;
    char* x_string;
    char* y_string;
    int end_site;
    unsigned int error;
    int real_y_start;
    long long overlap_length;
    int extra_begin, extra_end;
    long long o_len;


    ///j负责遍历整个overlap list
    for (j = 0; j < overlap_list->length; j++)
    {
        
        y_id = overlap_list->list[j].y_id;
        y_strand = overlap_list->list[j].y_pos_strand;
        y_readLen = Get_READ_LENGTH((*R_INF), y_id);


        
        
        
        //i负责每个overlap里面的window
        //倒着找
        //倒着用结束位置矫正
        for (i = overlap_list->list[j].w_list_length - 1; i >= 0; i--)
        {
            ///找到第一个匹配的window
            if(overlap_list->list[j].w_list[i].y_end != -1)
            {
                ///note!!! need notification
                ///total_y_start = overlap_list->list[j].w_list[i].y_end + 1;
                ///this is the actual end postion in ystring
                total_y_start = overlap_list->list[j].w_list[i].y_end - overlap_list->list[j].w_list[i].extra_begin + 1;

                

                ///k遍历匹配window右侧所有不匹配的window
                ///如果i匹配，则k从i+1开始
                ///知道第一个匹配的window结束
                for (k = i + 1; k < overlap_list->list[j].w_list_length && overlap_list->list[j].w_list[k].y_end == -1; k++)
                {
        /**
                    if(memcmp("m64013_190324_024932/92733922/ccs", Get_NAME((*R_INF),overlap_list->list[j].x_id), 
                Get_NAME_LENGTH((*R_INF), overlap_list->list[j].x_id)) == 0)
        {
                if(memcmp("m64013_190324_024932/123996697/ccs", Get_NAME((*R_INF),overlap_list->list[j].y_id), 
        Get_NAME_LENGTH((*R_INF), overlap_list->list[j].y_id)) == 0)
            {
                fprintf(stderr, "total_y_start: %d, x_start: %lld, x_end: %lld, y_start: %lld, y_end: %lld, match: %d, cigar.length: %d, error: %d\n", 
        total_y_start, overlap_list->list[j].w_list[i].x_start, overlap_list->list[j].w_list[i].x_end, 
        overlap_list->list[j].w_list[i].y_start, overlap_list->list[j].w_list[i].y_end,
        overlap_list->list[j].w_list[i].y_end - overlap_list->list[j].w_list[i].y_start + 1,
        overlap_list->list[j].w_list[i].cigar.length, overlap_list->list[j].w_list[i].error);
            }
        }
        **/


                    extra_begin = extra_end = 0;

                    ///y_start有可能大于y_readLen
                    ///这多发于最后一个window长度仅为几，而前面一个window的结束位置也超过了y_readLen-1
                    ///这个时候做动态规划会给超过的部分补N
                    if (total_y_start >= y_readLen)
                    {
                        break;
                    }
                    
                    ///there is no problem for x
                    x_start = overlap_list->list[j].w_list[k].x_start;
                    x_end = overlap_list->list[j].w_list[k].x_end;
                    x_len = x_end - x_start + 1;
                    // /****************************may have bugs********************************/
                    // ///threshold = x_len * THRESHOLD_RATE;
                    // threshold = overlap_list->list[j].w_list[k].error_threshold;
                    // /****************************may have bugs********************************/
                    // /****************************may have bugs********************************/
                    // threshold = Adjust_Threshold(threshold, x_len);
                    // /****************************may have bugs********************************/
                    threshold = double_error_threshold(overlap_list->list[j].w_list[k].error_threshold, x_len);
                    

                    y_start = total_y_start;
                    Window_Len = x_len + (threshold << 1);
                    determine_overlap_region(threshold, y_start, y_id, Window_Len, R_INF, 
                    &extra_begin, &extra_end, &y_start, &o_len);

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
                    
                    ///error等于-1说明没匹配
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
                    ///total_y_start = y_start + end_site + 1;
                    total_y_start = y_start + end_site - extra_begin + 1;

                }
                
            }
            
        }
        



        
        ///continue;
        ///i负责每个overlap里面的window
        ///正着找
        ///用起始位置矫正
        for (i = 0; i < overlap_list->list[j].w_list_length; i++)
        {
            ///找到第一个匹配的window
            ///首先这个window要匹配
            ///其次不要是第一个window，这没意义
            ///最后他之前的那个window必须是不匹配，如果之前那个window匹配，也没意义
            if(overlap_list->list[j].w_list[i].y_end != -1 && i != 0 && overlap_list->list[j].w_list[i - 1].y_end == -1)
            {
                
                ///判断这个匹配的window的起始位置有没有被计算出来
                ///如果没有，就需要重新计算
                if(overlap_list->list[j].w_list[i].cigar.length == -1)
                {
                    ///there is no problem for x
                    x_start = overlap_list->list[j].w_list[i].x_start;
                    x_end = overlap_list->list[j].w_list[i].x_end;
                    x_len = x_end - x_start + 1;
                    /****************************may have bugs********************************/
                    ///threshold = x_len * THRESHOLD_RATE;
                    threshold = overlap_list->list[j].w_list[i].error_threshold;
                    /****************************may have bugs********************************/
                    /****************************may have bugs********************************/
                    threshold = Adjust_Threshold(threshold, x_len);
                    /****************************may have bugs********************************/
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




                    
                    

                    ///到这里y_start已经被正确计算出来了
                    if (error != (unsigned int)-1)
                    {
                        
                        if (end_site == Window_Len - 1 || real_y_start == 0)
                        {
                            
                            if(fix_boundary(x_string, x_len, threshold, y_start, real_y_start, end_site,
                            extra_begin, extra_end, y_id, Window_Len, R_INF, dumy, y_strand, error,
                            &y_start, &real_y_start, &end_site,
                            &extra_begin, &extra_end, &error))
                            {
                                ///fprintf(stderr, "old_error: %d, new_error: %d\n", overlap_list->list[j].w_list[i].error, error);
                                overlap_list->list[j].w_list[i].error = error;
                                overlap_list->list[j].w_list[i].extra_begin = extra_begin;
                                overlap_list->list[j].w_list[i].extra_end = extra_end;
                            }
                            
                            
                        }
                        

                         
                        generate_cigar(dumy->path, dumy->path_length, &(overlap_list->list[j].w_list[i]),
                        &real_y_start, &end_site, &error, x_string, x_len, y_string);   

                        if(real_y_start < extra_begin || end_site >= Window_Len - extra_end)
                        {
                            fprintf(stderr, "\nreal_y_start: %d, extra_begin: %d\n", 
                            real_y_start, extra_begin);

                            fprintf(stderr, "end_site: %d, Window_Len: %d, extra_end: %d\n", 
                            end_site, Window_Len, extra_end);                            
                        }
                        
                        
                        ///note!!! need notification
                        ///real_y_start = y_start + real_y_start;
                        real_y_start = y_start + real_y_start - extra_begin;
                        overlap_list->list[j].w_list[i].y_start = real_y_start;  
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


                
                ///再次矫正的基础位置就是real_y_start
                total_y_end = real_y_start - 1;
                ///k遍历匹配window左侧所有不匹配的window
                ///如果i匹配，则k从i-1开始
                ///直到第一个匹配的window结束
                ///因为i!=0，所以k的大小不用担心
                for (k = i - 1; k >= 0 && overlap_list->list[j].w_list[k].y_end == -1; k--)
                {  
                    ///there is no problem in x
                    x_start = overlap_list->list[j].w_list[k].x_start;
                    x_end = overlap_list->list[j].w_list[k].x_end;
                    x_len = x_end - x_start + 1;
                    // /****************************may have bugs********************************/
                    // ///threshold = x_len * THRESHOLD_RATE;
                    // threshold = overlap_list->list[j].w_list[k].error_threshold;
                    // /****************************may have bugs********************************/
                    // /****************************may have bugs********************************/
                    // threshold = Adjust_Threshold(threshold, x_len);
                    // /****************************may have bugs********************************/
                    threshold = double_error_threshold(overlap_list->list[j].w_list[k].error_threshold, x_len);

                    Window_Len = x_len + (threshold << 1);

                    if(total_y_end <= 0)
                    {
                        break;
                    }

                    ///y_start may less than 0
                    y_start = total_y_end - x_len + 1;
                    determine_overlap_region(threshold, y_start, y_id, Window_Len, R_INF,
                    &extra_begin, &extra_end, &y_start, &o_len);

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


                    ///error等于-1说明没匹配
                    if (error!=(unsigned int)-1)
                    { 

                        if (end_site == Window_Len - 1 || real_y_start == 0)
                        {
                            fix_boundary(x_string, x_len, threshold, y_start, real_y_start, end_site,
                            extra_begin, extra_end, y_id, Window_Len, R_INF, dumy, y_strand, error,
                            &y_start, &real_y_start, &end_site,
                            &extra_begin, &extra_end, &error);
                        }

                        generate_cigar(dumy->path, dumy->path_length, &(overlap_list->list[j].w_list[k]),
                        &real_y_start, &end_site, &error, x_string, x_len, y_string);  

                        if(real_y_start < extra_begin || end_site >= Window_Len - extra_end)
                        {
                            fprintf(stderr, "\nreal_y_start: %d, extra_begin: %d\n", 
                            real_y_start, extra_begin);

                            fprintf(stderr, "end_site: %d, Window_Len: %d, extra_end: %d\n", 
                            end_site, Window_Len, extra_end);                            
                        }

                        ///y_start has no shift, but y_end has shift           
                        overlap_list->list[j].w_list[k].y_start = y_start + real_y_start - extra_begin;
                        overlap_list->list[j].w_list[k].y_end = y_start + end_site;
                        overlap_list->list[j].w_list[k].error = error;
                        overlap_list->list[j].align_length += x_len;
                        ///note!!! need notification
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
    
    int pre_threshold;
    long long tLen, tError;
    double error_rate;
    ///j负责遍历整个overlap list
    for (j = 0; j < overlap_list->length; j++)
    {
        y_id = overlap_list->list[j].y_id;
        y_strand = overlap_list->list[j].y_pos_strand;
        y_readLen = Get_READ_LENGTH((*R_INF), y_id);
        overlap_length = overlap_list->list[j].x_pos_e - overlap_list->list[j].x_pos_s + 1;
        overlap_list->list[j].is_match = 0;

    /**
        if(memcmp("m64013_190324_024932/92733922/ccs", Get_NAME((*R_INF),overlap_list->list[j].x_id), 
     Get_NAME_LENGTH((*R_INF), overlap_list->list[j].x_id)) == 0)
    {
            if(memcmp("m64013_190324_024932/123996697/ccs", Get_NAME((*R_INF),overlap_list->list[j].y_id), 
     Get_NAME_LENGTH((*R_INF), overlap_list->list[j].y_id)) == 0)
        {
            fprintf(stderr, "##############y_name: %.*s, error_rate: %f, is_match: %d, overlap_length: %d, align_length: %d##################\n",
            Get_NAME_LENGTH((*R_INF), overlap_list->list[j].y_id), Get_NAME((*R_INF),overlap_list->list[j].y_id),
            trim_error_rate(overlap_list, j), overlap_list->list[j].is_match, overlap_length, 
            overlap_list->list[j].align_length);
            for (i = 0; i < overlap_list->list[j].w_list_length; i++)
            {
                fprintf(stderr, "x_start: %lld, x_end: %lld, y_start: %lld, y_end: %lld, match: %d, cigar.length: %d, error: %d\n", 
                overlap_list->list[j].w_list[i].x_start, overlap_list->list[j].w_list[i].x_end, 
                overlap_list->list[j].w_list[i].y_start, overlap_list->list[j].w_list[i].y_end,
                overlap_list->list[j].w_list[i].y_end - overlap_list->list[j].w_list[i].y_start + 1,
                overlap_list->list[j].w_list[i].cigar.length, overlap_list->list[j].w_list[i].error);
                
            }
        }
    }
    **/



        ///only calculate cigar for high quality overlaps
        if (overlap_length * OVERLAP_THRESHOLD_FILTER <=  overlap_list->list[j].align_length)
        {
            
            for (i = 0; i < overlap_list->list[j].w_list_length; i++)
            {
                ///判断cigar是否被计算
                ///没被计算过就重算
                ///第一个条件是判断这个窗口是否匹配
                if(overlap_list->list[j].w_list[i].y_end != -1)
                {
                    if(overlap_list->list[j].w_list[i].cigar.length == -1)
                    {
                        ///there is no problem for x
                        x_start = overlap_list->list[j].w_list[i].x_start;
                        x_end = overlap_list->list[j].w_list[i].x_end;
                        x_len = x_end - x_start + 1;
                        /****************************may have bugs********************************/
                        ///threshold = x_len * THRESHOLD_RATE;
                        threshold = overlap_list->list[j].w_list[i].error_threshold;
                        /****************************may have bugs********************************/
                        /****************************may have bugs********************************/
                        threshold = Adjust_Threshold(threshold, x_len);
                        /****************************may have bugs********************************/
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


                        ///到这里y_start已经被正确计算出来了
                        if (error != (unsigned int)-1)
                        {

                            if (end_site == Window_Len - 1 || real_y_start == 0)
                            {
                                
                                if(fix_boundary(x_string, x_len, threshold, y_start, real_y_start, end_site,
                                extra_begin, extra_end, y_id, Window_Len, R_INF, dumy, y_strand, error,
                                &y_start, &real_y_start, &end_site,
                                &extra_begin, &extra_end, &error))
                                {
                                    ///fprintf(stderr, "old_error: %d, new_error: %d\n", overlap_list->list[j].w_list[i].error, error);
                                    overlap_list->list[j].w_list[i].error = error;
                                    overlap_list->list[j].w_list[i].extra_begin = extra_begin;
                                    overlap_list->list[j].w_list[i].extra_end = extra_end;
                                }
                                
                            }

                            generate_cigar(dumy->path, dumy->path_length, &(overlap_list->list[j].w_list[i]),
                            &real_y_start, &end_site, &error, x_string, x_len, y_string);    


                            if(real_y_start < extra_begin || end_site >= Window_Len - extra_end)
                            {
                                fprintf(stderr, "\nreal_y_start: %d, extra_begin: %d\n", 
                                real_y_start, extra_begin);

                                fprintf(stderr, "end_site: %d, Window_Len: %d, extra_end: %d\n", 
                                end_site, Window_Len, extra_end);                            
                            }

                            

                            ///note!!! need notification
                            ///real_y_start = y_start + real_y_start;
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
                else ///try to calculate using higher threshold
                {
                    ///there is no problem for x
                    ///there is no problem for x
                    x_start = overlap_list->list[j].w_list[i].x_start;
                    x_end = overlap_list->list[j].w_list[i].x_end;
                    x_len = x_end - x_start + 1;
                    ///double the threshold
                    pre_threshold = overlap_list->list[j].w_list[i].error_threshold;
                    threshold = double_error_threshold(pre_threshold, x_len);
                    // /****************************may have bugs********************************/
                    // ///pre_threshold = x_len * THRESHOLD_RATE;
                    // pre_threshold = overlap_list->list[j].w_list[i].error_threshold;
                    // /****************************may have bugs********************************/
                    // /****************************may have bugs********************************/
                    // pre_threshold = Adjust_Threshold(pre_threshold, x_len);
                    // /****************************may have bugs********************************/
                    // threshold = pre_threshold * 2;
                    // ///may have some bugs
                    // if(x_len >= 300 && threshold < THRESHOLD_MAX_SIZE)
                    // {
                    //     threshold = THRESHOLD_MAX_SIZE;
                    // }
                    // if(threshold > THRESHOLD_MAX_SIZE)
                    // {
                    //     threshold = THRESHOLD_MAX_SIZE;
                    // }
                    Window_Len = x_len + (threshold << 1);



                    ///if the previous window is mapped
                    if(i > 0 && overlap_list->list[j].w_list[i - 1].y_end != -1)
                    {
                        y_start = overlap_list->list[j].w_list[i - 1].y_end + 1;
                        determine_overlap_region(threshold, y_start, y_id, Window_Len, R_INF, 
                        &extra_begin, &extra_end, &y_start, &o_len);
                    }///if the next window is mapped
                    else if(i < overlap_list->list[j].w_list_length - 1 && overlap_list->list[j].w_list[i + 1].y_end != -1)
                    {
                        y_start = overlap_list->list[j].w_list[i + 1].y_start - 1 - x_len + 1;
                        determine_overlap_region(threshold, y_start, y_id, Window_Len, R_INF, 
                        &extra_begin, &extra_end, &y_start, &o_len);
                    }
                    else///if the previous window and next window are not mapped, using the y_start itself
                    {
                        ///y_start is the real y_start
                        y_start = overlap_list->list[j].w_list[i].y_start;
                        /// since y_start has already substacted pre_threshold
                        ///here we just need to substact threshold - pre_threshold
                        determine_overlap_region(threshold - pre_threshold, y_start, y_id, Window_Len, R_INF,
                        &extra_begin, &extra_end, &y_start, &o_len);
                    }
                    

                    fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, R_INF, y_id, extra_begin, extra_end);

                    x_string = g_read->seq + x_start;
                    y_string = dumy->overlap_region;

                    ///note!!! need notification
                    end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
                    &(dumy->path_length), dumy->matrix_bit, dumy->path, -1, -1);

                    if (error!=(unsigned int)-1)
                    { 
                        
                        if (end_site == Window_Len - 1 || real_y_start == 0)
                        {
                            fix_boundary(x_string, x_len, threshold, y_start, real_y_start, end_site,
                            extra_begin, extra_end, y_id, Window_Len, R_INF, dumy, y_strand, error,
                            &y_start, &real_y_start, &end_site,
                            &extra_begin, &extra_end, &error);
                        }

                        generate_cigar(dumy->path, dumy->path_length, &(overlap_list->list[j].w_list[i]),
                        &real_y_start, &end_site, &error, x_string, x_len, y_string);  


                        if(real_y_start < extra_begin || end_site >= Window_Len - extra_end)
                        {
                            fprintf(stderr, "\nreal_y_start: %d, extra_begin: %d, error: %d, Window_Len: %d, x_len: %d\n", 
                            real_y_start, extra_begin, error, Window_Len, x_len);

                            fprintf(stderr, "end_site: %d, Window_Len: %d, extra_end: %d\n", 
                            end_site, Window_Len, extra_end);                            
                        }

                        real_y_start = y_start + real_y_start - extra_begin;
                        overlap_list->list[j].w_list[i].y_start = real_y_start;  
                        overlap_list->list[j].w_list[i].y_end = y_start + end_site - extra_begin; 
                        overlap_list->list[j].w_list[i].error = error;
                        overlap_list->list[j].align_length += x_len;
                        overlap_list->list[j].w_list[i].extra_begin = extra_begin;
                        overlap_list->list[j].w_list[i].extra_end = extra_end;
                    }
                }
            }


            error_rate = trim_error_rate(overlap_list, j);

        

            ///if(error_rate <= 0.015)
            if(error_rate <= 0.025)
            {
                ///overlap_list->mapped_overlaps++;
                overlap_list->mapped_overlaps_length += overlap_length;
                overlap_list->list[j].is_match = 1;
            }
        }
    }
    
    



    /**
    ///j负责遍历整个overlap list
    for (j = 0; j < overlap_list->length; j++)
    {
        y_id = overlap_list->list[j].y_id;
        y_strand = overlap_list->list[j].y_pos_strand;
        y_readLen = Get_READ_LENGTH((*R_INF), y_id);
        overlap_length = overlap_list->list[j].x_pos_e - overlap_list->list[j].x_pos_s + 1;

        ///only calculate cigar for high quality overlaps
        if (overlap_length * OVERLAP_THRESHOLD_FILTER <=  overlap_list->list[j].align_length)
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
                        fprintf(stderr, "j: %d, i: %d, y_id: %d, y_start: %d, y_end: %d\n", j, i, y_id, y_start, y_end);
                    }
                }

            }
        }


        // for (i = 0; i < overlap_list->list[j].w_list_length; i++)
        // {

        //     if(overlap_list->list[j].w_list[i].x_end - overlap_list->list[j].w_list[i].x_start + 1 != WINDOW && 
        //     overlap_list->list[j].w_list[i].x_end != overlap_list->list[j].x_pos_e && 
        //     overlap_list->list[j].w_list[i].x_start != overlap_list->list[j].x_pos_s)
        //     {
        //         fprintf(stderr, "x_start:%d, x_end: %d, g_read->length: %d, x_pos_s: %d, x_pos_e: %d\n", 
        //         overlap_list->list[j].w_list[i].x_start, 
        //         overlap_list->list[j].w_list[i].x_end, 
        //         g_read->length,
        //         overlap_list->list[j].x_pos_s,
        //         overlap_list->list[j].x_pos_e);
        //     }
        // }
    }
    **/
    
    

    
    
}


inline void recalcate_window_simple(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read)
{
    long long j, k, i;
    long long Len_x;
    int threshold;
    long long y_len;
    long long currentIDLen;
    long long matches;
    long long y_id;
    int y_strand;
    long long y_readLen;
    long long x_start;
    long long x_end;
    long long x_len;
    long long total_y_start;
    long long total_y_end;
    long long y_start;
    long long y_end;
    long long Window_Len;
    char* x_string;
    char* y_string;
    int end_site;
    unsigned int error;
    int real_y_start;
    long long overlap_length;
    int extra_begin, extra_end;
    long long o_len;


    ///j负责遍历整个overlap list
    for (j = 0; j < overlap_list->length; j++)
    {
        
        y_id = overlap_list->list[j].y_id;
        y_strand = overlap_list->list[j].y_pos_strand;
        y_readLen = Get_READ_LENGTH((*R_INF), y_id);


        
        
        
        //i负责每个overlap里面的window
        //倒着找
        //倒着用结束位置矫正
        for (i = overlap_list->list[j].w_list_length - 1; i >= 0; i--)
        {
            ///找到第一个匹配的window
            if(overlap_list->list[j].w_list[i].y_end != -1)
            {
                ///note!!! need notification
                ///total_y_start = overlap_list->list[j].w_list[i].y_end + 1;
                ///this is the actual end postion in ystring
                total_y_start = overlap_list->list[j].w_list[i].y_end - overlap_list->list[j].w_list[i].extra_begin + 1;

                

                ///k遍历匹配window右侧所有不匹配的window
                ///如果i匹配，则k从i+1开始
                ///知道第一个匹配的window结束
                for (k = i + 1; k < overlap_list->list[j].w_list_length && overlap_list->list[j].w_list[k].y_end == -1; k++)
                {
                    extra_begin = extra_end = 0;

                    ///y_start有可能大于y_readLen
                    ///这多发于最后一个window长度仅为几，而前面一个window的结束位置也超过了y_readLen-1
                    ///这个时候做动态规划会给超过的部分补N
                    if (total_y_start >= y_readLen)
                    {
                        break;
                    }
                    
                    ///there is no problem for x
                    x_start = overlap_list->list[j].w_list[k].x_start;
                    x_end = overlap_list->list[j].w_list[k].x_end;
                    x_len = x_end - x_start + 1;
                    // /****************************may have bugs********************************/
                    // ///threshold = x_len * THRESHOLD_RATE;
                    // threshold = overlap_list->list[j].w_list[k].error_threshold;
                    // /****************************may have bugs********************************/
                    // /****************************may have bugs********************************/
                    // threshold = Adjust_Threshold(threshold, x_len);
                    // /****************************may have bugs********************************/
                    threshold = double_error_threshold(overlap_list->list[j].w_list[k].error_threshold, x_len);
                    

                    y_start = total_y_start;
                    Window_Len = x_len + (threshold << 1);
                    if(!determine_overlap_region(threshold, y_start, y_id, Window_Len, R_INF, 
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
                    
                    ///error等于-1说明没匹配
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
                    ///total_y_start = y_start + end_site + 1;
                    total_y_start = y_start + end_site - extra_begin + 1;

                }
                
            }
            
        }
        



        
        ///continue;
        ///i负责每个overlap里面的window
        ///正着找
        ///用起始位置矫正
        for (i = 0; i < overlap_list->list[j].w_list_length; i++)
        {
            ///找到第一个匹配的window
            ///首先这个window要匹配
            ///其次不要是第一个window，这没意义
            ///最后他之前的那个window必须是不匹配，如果之前那个window匹配，也没意义
            if(overlap_list->list[j].w_list[i].y_end != -1 && i != 0 && overlap_list->list[j].w_list[i - 1].y_end == -1)
            {
                
                ///判断这个匹配的window的起始位置有没有被计算出来
                ///如果没有，就需要重新计算
                if(overlap_list->list[j].w_list[i].cigar.length == -1)
                {
                    ///there is no problem for x
                    x_start = overlap_list->list[j].w_list[i].x_start;
                    x_end = overlap_list->list[j].w_list[i].x_end;
                    x_len = x_end - x_start + 1;
                    /****************************may have bugs********************************/
                    ///threshold = x_len * THRESHOLD_RATE;
                    threshold = overlap_list->list[j].w_list[i].error_threshold;
                    /****************************may have bugs********************************/
                    /****************************may have bugs********************************/
                    threshold = Adjust_Threshold(threshold, x_len);
                    /****************************may have bugs********************************/
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




                    
                    

                    ///到这里y_start已经被正确计算出来了
                    if (error != (unsigned int)-1)
                    {
                        
                        if (end_site == Window_Len - 1 || real_y_start == 0)
                        {
                            
                            if(fix_boundary(x_string, x_len, threshold, y_start, real_y_start, end_site,
                            extra_begin, extra_end, y_id, Window_Len, R_INF, dumy, y_strand, error,
                            &y_start, &real_y_start, &end_site,
                            &extra_begin, &extra_end, &error))
                            {
                                ///fprintf(stderr, "old_error: %d, new_error: %d\n", overlap_list->list[j].w_list[i].error, error);
                                overlap_list->list[j].w_list[i].error = error;
                                overlap_list->list[j].w_list[i].extra_begin = extra_begin;
                                overlap_list->list[j].w_list[i].extra_end = extra_end;
                            }
                            
                            
                        }
                        

                         
                        generate_cigar(dumy->path, dumy->path_length, &(overlap_list->list[j].w_list[i]),
                        &real_y_start, &end_site, &error, x_string, x_len, y_string);   

                        if(real_y_start < extra_begin || end_site >= Window_Len - extra_end)
                        {
                            fprintf(stderr, "\nreal_y_start: %d, extra_begin: %d\n", 
                            real_y_start, extra_begin);

                            fprintf(stderr, "end_site: %d, Window_Len: %d, extra_end: %d\n", 
                            end_site, Window_Len, extra_end);                            
                        }
                        
                        
                        ///note!!! need notification
                        ///real_y_start = y_start + real_y_start;
                        real_y_start = y_start + real_y_start - extra_begin;
                        overlap_list->list[j].w_list[i].y_start = real_y_start;  
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


                
                ///再次矫正的基础位置就是real_y_start
                total_y_end = real_y_start - 1;
                ///k遍历匹配window左侧所有不匹配的window
                ///如果i匹配，则k从i-1开始
                ///直到第一个匹配的window结束
                ///因为i!=0，所以k的大小不用担心
                for (k = i - 1; k >= 0 && overlap_list->list[j].w_list[k].y_end == -1; k--)
                {  
                    ///there is no problem in x
                    x_start = overlap_list->list[j].w_list[k].x_start;
                    x_end = overlap_list->list[j].w_list[k].x_end;
                    x_len = x_end - x_start + 1;
                    // /****************************may have bugs********************************/
                    // ///threshold = x_len * THRESHOLD_RATE;
                    // threshold = overlap_list->list[j].w_list[k].error_threshold;
                    // /****************************may have bugs********************************/
                    // /****************************may have bugs********************************/
                    // threshold = Adjust_Threshold(threshold, x_len);
                    // /****************************may have bugs********************************/
                    threshold = double_error_threshold(overlap_list->list[j].w_list[k].error_threshold, x_len);

                    Window_Len = x_len + (threshold << 1);

                    if(total_y_end <= 0)
                    {
                        break;
                    }

                    ///y_start may less than 0
                    y_start = total_y_end - x_len + 1;
                    if(!determine_overlap_region(threshold, y_start, y_id, Window_Len, R_INF,
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


                    ///error等于-1说明没匹配
                    if (error!=(unsigned int)-1)
                    { 

                        if (end_site == Window_Len - 1 || real_y_start == 0)
                        {
                            fix_boundary(x_string, x_len, threshold, y_start, real_y_start, end_site,
                            extra_begin, extra_end, y_id, Window_Len, R_INF, dumy, y_strand, error,
                            &y_start, &real_y_start, &end_site,
                            &extra_begin, &extra_end, &error);
                        }

                        generate_cigar(dumy->path, dumy->path_length, &(overlap_list->list[j].w_list[k]),
                        &real_y_start, &end_site, &error, x_string, x_len, y_string);  

                        if(real_y_start < extra_begin || end_site >= Window_Len - extra_end)
                        {
                            fprintf(stderr, "\nreal_y_start: %d, extra_begin: %d\n", 
                            real_y_start, extra_begin);

                            fprintf(stderr, "end_site: %d, Window_Len: %d, extra_end: %d\n", 
                            end_site, Window_Len, extra_end);                            
                        }

                        ///y_start has no shift, but y_end has shift           
                        overlap_list->list[j].w_list[k].y_start = y_start + real_y_start - extra_begin;
                        overlap_list->list[j].w_list[k].y_end = y_start + end_site;
                        overlap_list->list[j].w_list[k].error = error;
                        overlap_list->list[j].align_length += x_len;
                        ///note!!! need notification
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
    
    int pre_threshold;
    long long tLen, tError;
    double error_rate;
    ///j负责遍历整个overlap list
    for (j = 0; j < overlap_list->length; j++)
    {
        y_id = overlap_list->list[j].y_id;
        y_strand = overlap_list->list[j].y_pos_strand;
        y_readLen = Get_READ_LENGTH((*R_INF), y_id);
        overlap_length = overlap_list->list[j].x_pos_e - overlap_list->list[j].x_pos_s + 1;
        overlap_list->list[j].is_match = 0;


        ///only calculate cigar for high quality overlaps
        if (overlap_length * OVERLAP_THRESHOLD_FILTER <=  overlap_list->list[j].align_length)
        {
            
            for (i = 0; i < overlap_list->list[j].w_list_length; i++)
            {
                ///判断cigar是否被计算
                ///没被计算过就重算
                ///第一个条件是判断这个窗口是否匹配
                if(overlap_list->list[j].w_list[i].y_end != -1)
                {
                    if(overlap_list->list[j].w_list[i].cigar.length == -1)
                    {
                        ///there is no problem for x
                        x_start = overlap_list->list[j].w_list[i].x_start;
                        x_end = overlap_list->list[j].w_list[i].x_end;
                        x_len = x_end - x_start + 1;
                        /****************************may have bugs********************************/
                        ///threshold = x_len * THRESHOLD_RATE;
                        threshold = overlap_list->list[j].w_list[i].error_threshold;
                        /****************************may have bugs********************************/
                        /****************************may have bugs********************************/
                        threshold = Adjust_Threshold(threshold, x_len);
                        /****************************may have bugs********************************/
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

                        // if(error != overlap_list->list[j].w_list[i].error)
                        // {
                        //     fprintf(stderr, "error\n");
                        // }

                        ///到这里y_start已经被正确计算出来了
                        if (error != (unsigned int)-1)
                        {

                            if (end_site == Window_Len - 1 || real_y_start == 0)
                            {
                                
                                if(fix_boundary(x_string, x_len, threshold, y_start, real_y_start, end_site,
                                extra_begin, extra_end, y_id, Window_Len, R_INF, dumy, y_strand, error,
                                &y_start, &real_y_start, &end_site,
                                &extra_begin, &extra_end, &error))
                                {
                                    ///fprintf(stderr, "old_error: %d, new_error: %d\n", overlap_list->list[j].w_list[i].error, error);
                                    overlap_list->list[j].w_list[i].error = error;
                                    overlap_list->list[j].w_list[i].extra_begin = extra_begin;
                                    overlap_list->list[j].w_list[i].extra_end = extra_end;
                                }
                                
                            }

                            generate_cigar(dumy->path, dumy->path_length, &(overlap_list->list[j].w_list[i]),
                            &real_y_start, &end_site, &error, x_string, x_len, y_string);    


                            if(real_y_start < extra_begin || end_site >= Window_Len - extra_end)
                            {
                                fprintf(stderr, "\nreal_y_start: %d, extra_begin: %d\n", 
                                real_y_start, extra_begin);

                                fprintf(stderr, "end_site: %d, Window_Len: %d, extra_end: %d\n", 
                                end_site, Window_Len, extra_end);                            
                            }

                            

                            ///note!!! need notification
                            ///real_y_start = y_start + real_y_start;
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

            ///error_rate = trim_error_rate(overlap_list, j);
            error_rate = non_trim_error_rate(overlap_list, j, R_INF, dumy, g_read);
            

            ///if(error_rate <= 0.015)
            if(error_rate <= 0.03)
            {
                ///overlap_list->mapped_overlaps++;
                overlap_list->mapped_overlaps_length += overlap_length;
                overlap_list->list[j].is_match = 1;
            }
            else if(error_rate <= 0.045)
            {
                overlap_list->list[j].is_match = 3;
            }
        }
    }
    
    



    /**
    ///j负责遍历整个overlap list
    for (j = 0; j < overlap_list->length; j++)
    {
        y_id = overlap_list->list[j].y_id;
        y_strand = overlap_list->list[j].y_pos_strand;
        y_readLen = Get_READ_LENGTH((*R_INF), y_id);
        overlap_length = overlap_list->list[j].x_pos_e - overlap_list->list[j].x_pos_s + 1;

        ///only calculate cigar for high quality overlaps
        if (overlap_length * OVERLAP_THRESHOLD_FILTER <=  overlap_list->list[j].align_length)
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
                        fprintf(stderr, "j: %d, i: %d, y_id: %d, y_start: %d, y_end: %d\n", j, i, y_id, y_start, y_end);
                    }
                }

            }
        }


        // for (i = 0; i < overlap_list->list[j].w_list_length; i++)
        // {

        //     if(overlap_list->list[j].w_list[i].x_end - overlap_list->list[j].w_list[i].x_start + 1 != WINDOW && 
        //     overlap_list->list[j].w_list[i].x_end != overlap_list->list[j].x_pos_e && 
        //     overlap_list->list[j].w_list[i].x_start != overlap_list->list[j].x_pos_s)
        //     {
        //         fprintf(stderr, "x_start:%d, x_end: %d, g_read->length: %d, x_pos_s: %d, x_pos_e: %d\n", 
        //         overlap_list->list[j].w_list[i].x_start, 
        //         overlap_list->list[j].w_list[i].x_end, 
        //         g_read->length,
        //         overlap_list->list[j].x_pos_s,
        //         overlap_list->list[j].x_pos_e);
        //     }
        // }
    }
    **/
    
    

    
    
}




inline void recalcate_window(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read)
{
    long long j, k, i;
    long long Len_x;
    int threshold;
    long long y_len;
    long long currentIDLen;
    long long matches;
    long long y_id;
    int y_strand;
    long long y_readLen;
    long long x_start;
    long long x_end;
    long long x_len;
    long long total_y_start;
    long long total_y_end;
    long long y_start;
    long long y_end;
    long long Window_Len;
    char* x_string;
    char* y_string;
    int end_site;
    unsigned int error;
    int real_y_start;
    long long overlap_length;
    int extra_begin, extra_end;
    long long o_len;


    ///j负责遍历整个overlap list
    for (j = 0; j < overlap_list->length; j++)
    {
        
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
                for (k = i + 1; k < overlap_list->list[j].w_list_length && overlap_list->list[j].w_list[k].y_end == -1; k++)
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

                    if(!determine_overlap_region(threshold, y_start, y_id, Window_Len, R_INF, 
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
                    
                    ///error等于-1说明没匹配
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
        for (i = 0; i < overlap_list->list[j].w_list_length; i++)
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
                    /****************************may have bugs********************************/
                    threshold = overlap_list->list[j].w_list[i].error_threshold;
                    /****************************may have bugs********************************/
                    /****************************may have bugs********************************/
                    ///should not adjust threshold, since this window can be matched by the old threshold
                    ///threshold = Adjust_Threshold(threshold, x_len);
                    /****************************may have bugs********************************/
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

                        // if(real_y_start < extra_begin || end_site >= Window_Len - extra_end)
                        // {
                        //     fprintf(stderr, "\nreal_y_start: %d, extra_begin: %d\n", 
                        //     real_y_start, extra_begin);

                        //     fprintf(stderr, "end_site: %d, Window_Len: %d, extra_end: %d\n", 
                        //     end_site, Window_Len, extra_end);                            
                        // }
                        
                        
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
                ///k遍历匹配window左侧所有不匹配的window
                ///如果i匹配，则k从i-1开始
                ///直到第一个匹配的window结束
                ///因为i!=0，所以k的大小不用担心
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
                    if(!determine_overlap_region(threshold, y_start, y_id, Window_Len, R_INF,
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

                        // if(real_y_start < extra_begin || end_site >= Window_Len - extra_end)
                        // {
                        //     fprintf(stderr, "\nreal_y_start: %d, extra_begin: %d\n", 
                        //     real_y_start, extra_begin);

                        //     fprintf(stderr, "end_site: %d, Window_Len: %d, extra_end: %d\n", 
                        //     end_site, Window_Len, extra_end);                            
                        // }

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
    
    int pre_threshold;
    long long tLen, tError;
    double error_rate;
    for (j = 0; j < overlap_list->length; j++)
    {
        y_id = overlap_list->list[j].y_id;
        y_strand = overlap_list->list[j].y_pos_strand;
        y_readLen = Get_READ_LENGTH((*R_INF), y_id);
        overlap_length = overlap_list->list[j].x_pos_e - overlap_list->list[j].x_pos_s + 1;
        overlap_list->list[j].is_match = 0;


        ///only calculate cigar for high quality overlaps
        if (overlap_length * OVERLAP_THRESHOLD_FILTER <=  overlap_list->list[j].align_length)
        {
            
            for (i = 0; i < overlap_list->list[j].w_list_length; i++)
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
                        /****************************may have bugs********************************/
                        ///threshold = x_len * THRESHOLD_RATE;
                        threshold = overlap_list->list[j].w_list[i].error_threshold;
                        /****************************may have bugs********************************/
                        /****************************may have bugs********************************/
                        ///should not adjust threshold, since this window can be matched by the old threshold
                        ///threshold = Adjust_Threshold(threshold, x_len);
                        /****************************may have bugs********************************/
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

                        // if(error != overlap_list->list[j].w_list[i].error)
                        // {
                        //     fprintf(stderr, "error\n");
                        // }

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


                            // if(real_y_start < extra_begin || end_site >= Window_Len - extra_end)
                            // {
                            //     fprintf(stderr, "\nreal_y_start: %d, extra_begin: %d\n", 
                            //     real_y_start, extra_begin);

                            //     fprintf(stderr, "end_site: %d, Window_Len: %d, extra_end: %d\n", 
                            //     end_site, Window_Len, extra_end);                            
                            // }

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

            ///error_rate = trim_error_rate(overlap_list, j);
            error_rate = non_trim_error_rate(overlap_list, j, R_INF, dumy, g_read);
            

            ///if(error_rate <= 0.015)
            if(error_rate <= 0.03)
            {
                overlap_list->mapped_overlaps_length += overlap_length;
                overlap_list->list[j].is_match = 1;
            }
            else if(error_rate <= 0.045)
            {
                overlap_list->list[j].is_match = 3;
            }
        }
    }    
}


void debug_scan_cigar(overlap_region* sub_list)
{
    long long i;
    int f_err, b_err, fLen, bLen, xLen;
    for (i = 0; i < sub_list->w_list_length; i++)
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
            fprintf(stderr, "sbsbsbsbsb1\n");
        }

        scan_cigar(&(sub_list->w_list[i].cigar), &b_err, 
        WINDOW, 1);
        scan_cigar(&(sub_list->w_list[i].cigar), &f_err, 
        WINDOW, 0);

        if(b_err != sub_list->w_list[i].error || f_err != sub_list->w_list[i].error)
        {
            fprintf(stderr, "sbsbsbsbsb2\n");
        }

        scan_cigar_interval(&(sub_list->w_list[i].cigar), &b_err, 0, xLen-1);
        if(b_err != sub_list->w_list[i].error)
        {
            fprintf(stderr, "sbsbsbsbsb3\n");
        }

        fLen = xLen / 3;
        scan_cigar(&(sub_list->w_list[i].cigar), &f_err, fLen, 0);
        scan_cigar_interval(&(sub_list->w_list[i].cigar), &b_err, 0, fLen-1);
        if(f_err != b_err)
        {
            fprintf(stderr, "sbsbsbsbsb4\n");
        }

        fLen = xLen / 3;
        scan_cigar(&(sub_list->w_list[i].cigar), &f_err, fLen, 1);
        scan_cigar_interval(&(sub_list->w_list[i].cigar), &b_err, xLen-fLen, xLen-1);
        if(f_err != b_err)
        {
            fprintf(stderr, "\nsbsbsbsbsb5\n");
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

void calculate_boundary_cigars(overlap_region* sub_list, All_reads* R_INF, Correct_dumy* dumy, 
UC_Read* g_read)
{
    resize_window_list_alloc(&(sub_list->boundary_cigars), sub_list->w_list_length - 1);
    int x_id = sub_list->x_id;
    int y_id = sub_list->y_id;
    int y_strand = sub_list->y_pos_strand;
    long long x_readLen = Get_READ_LENGTH((*R_INF), x_id);
    long long y_readLen = Get_READ_LENGTH((*R_INF), y_id);
    long long i, y_distance;
    int f_err, b_err, m_error;
    int scanLen = 10;
    long long boundaryLen = 200;
    long long single_sideLen = boundaryLen/2;
    long long force_useless_side = single_sideLen/2;
    long long L_useless_side, R_useless_side;
    int alpha = 1;
    long long y_start, x_start;
    long long y_end, x_end;
    long long yLen, xLen;
    long long leftLen, rightLen;
    long long threshold;
    int extra_begin, extra_end;
    long long o_len;
    char* x_string;
    char* y_string;
    int end_site;
    unsigned int error;
    int real_y_start;
    sub_list->boundary_cigars.length = sub_list->w_list_length - 1;
    ///the (i)-th boundary between the (i)-th window and the (i+1)-th window
    ///that means it includes (the tail of (i)-th window) and (the header of (i+1)-th window)
    ///note the (i)-th boundary is calculated at the (i)-th window
    for (i = 0; i + 1 < sub_list->w_list_length; i++)
    {
        ///if both of the two windows are not aligned
        ///it is not necessary to calculate the boundary
        ///if(sub_list->w_list[i].y_end == -1 && sub_list->w_list[i+1].y_end == -1)
        if(sub_list->w_list[i].y_end == -1 || sub_list->w_list[i+1].y_end == -1)
        {
            sub_list->boundary_cigars.buffer[i].error = -1;
            sub_list->boundary_cigars.buffer[i].y_end = -1;
            continue;
        }

        ///note if w_list[i+1].y_start or sub_list->w_list[i].y_end is -1
        ///y_distance might have some problems at the last of this function
        ///we need to deal with it carefully
        y_distance = sub_list->w_list[i+1].y_start - sub_list->w_list[i].y_end - 1;

        ///if two windows are aligned
        if(sub_list->w_list[i].y_end != -1 && sub_list->w_list[i+1].y_end != -1 && y_distance == 0)
        {
            ///scan backward
            scan_cigar(&(sub_list->w_list[i].cigar), &b_err, scanLen, 1);
            ///scan forward
            scan_cigar(&(sub_list->w_list[i+1].cigar), &f_err, scanLen, 0);
            if(b_err == 0 && f_err == 0)
            {
                sub_list->boundary_cigars.buffer[i].error = -2;
                sub_list->boundary_cigars.buffer[i].y_end = -1;
                continue;
            }
        }


        ///y_distance can be less than 0, or larger than 0
        if(sub_list->w_list[i].y_end != -1)
        {
            y_start = sub_list->w_list[i].y_end;
            x_start = sub_list->w_list[i].x_end;
        }///if the (i)-th window is not matched, have a look at the (i+1)-th window
        else if(sub_list->w_list[i+1].y_end != -1)
        {
            y_start = sub_list->w_list[i+1].y_start;
            x_start = sub_list->w_list[i+1].x_start;
        }///if both of these two windows are not matched, directly skip
        else
        {
            sub_list->boundary_cigars.buffer[i].error = -1;
            sub_list->boundary_cigars.buffer[i].y_end = -1;
            continue;
        }

        ///it seems we don't need to record x_start and y_start
        sub_list->boundary_cigars.buffer[i].extra_begin = x_start;
        sub_list->boundary_cigars.buffer[i].extra_end = y_start;

        ///leftLen and rightLen are used for x        
        ///x should be at [sub_list->w_list[i].x_start, sub_list->w_list[i+1].x_end]
        ///y shouldn't have limitation
        ///note that the x_start and x_end should not be -1 in any case
        ///up to now, x_start and y_start are not -1
        ///leftLen does not include x_start itself, rightLen does
        ///gnerally speaking, rightLen should be always larger than leftLen
        leftLen = MIN(MIN((x_start - sub_list->w_list[i].x_start), y_start), single_sideLen);
        rightLen = MIN(MIN((sub_list->w_list[i+1].x_end + 1 - x_start), y_readLen - y_start), 
        single_sideLen);

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


        threshold = xLen * THRESHOLD_RATE;
        threshold = Adjust_Threshold(threshold, xLen);
        threshold = double_error_threshold(threshold, xLen);

        yLen = xLen + (threshold << 1);
        if(!determine_overlap_region(threshold, y_start, y_id, yLen, R_INF, 
                    &extra_begin, &extra_end, &y_start, &o_len))
        {
            sub_list->boundary_cigars.buffer[i].error = -1;
            sub_list->boundary_cigars.buffer[i].y_end = -1;
            continue;
        }

        if(o_len < xLen)
        {
            sub_list->boundary_cigars.buffer[i].error = -1;
            sub_list->boundary_cigars.buffer[i].y_end = -1;
            continue;
        }

        fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, 
                    R_INF, y_id, extra_begin, extra_end);
        
        x_string = g_read->seq + x_start;
        y_string = dumy->overlap_region;

        end_site = Reserve_Banded_BPM_PATH(y_string, yLen, x_string, xLen, threshold, &error, 
        &real_y_start, &(dumy->path_length), dumy->matrix_bit, dumy->path, -1, -1);

        ///means this window is matched
        if (error!=(unsigned int)-1)
        {
            sub_list->boundary_cigars.buffer[i].x_start = x_start;
            sub_list->boundary_cigars.buffer[i].x_end = x_end;

            generate_cigar(dumy->path, dumy->path_length, &(sub_list->boundary_cigars.buffer[i]),
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
            if((i == 0) && (x_start == sub_list->w_list[0].x_start))
            {
                L_useless_side = 0;
            }
            ///last window
            if((i == sub_list->w_list_length - 2) && 
                (x_end == sub_list->w_list[sub_list->w_list_length - 1].x_end))
            {
                R_useless_side = 0;
            }

            if(leftLen <= L_useless_side || rightLen <= R_useless_side)
            {
                sub_list->boundary_cigars.buffer[i].error = -1;
                sub_list->boundary_cigars.buffer[i].y_end = -1;
                continue;
            }



            ///up to now, if we require (i)-th window and (i+1)-th window are matched
            ///boundary_cigars.buffer[i].cigar, w_list[i].cigar and w_list[i+1].cigar are avaiable
            ///get the error excluding the first and the last useless_side bases
            scan_cigar_interval(&(sub_list->boundary_cigars.buffer[i].cigar), &m_error, 
            L_useless_side, xLen-R_useless_side-1);
            scan_cigar(&(sub_list->w_list[i].cigar), &b_err, leftLen-L_useless_side, 1);
            scan_cigar(&(sub_list->w_list[i+1].cigar), &f_err, rightLen-R_useless_side, 0);

            if(f_err + b_err + y_distance + alpha < m_error)
            {
                sub_list->boundary_cigars.buffer[i].error = -1;
                sub_list->boundary_cigars.buffer[i].y_end = -1;
                continue;
            }

            sub_list->boundary_cigars.buffer[i].error = error;
            sub_list->boundary_cigars.buffer[i].y_start = y_start + real_y_start - extra_begin;
            sub_list->boundary_cigars.buffer[i].y_end = y_start + end_site - extra_begin;

            sub_list->boundary_cigars.buffer[i].x_start = x_start;
            sub_list->boundary_cigars.buffer[i].x_end = x_end;
            ///sub_list->boundary_cigars.buffer[i].error_threshold = useless_side;
            sub_list->boundary_cigars.buffer[i].extra_begin = L_useless_side;
            sub_list->boundary_cigars.buffer[i].extra_end = R_useless_side;
        }
        else
        {
            sub_list->boundary_cigars.buffer[i].error = -1;
            sub_list->boundary_cigars.buffer[i].y_end = -1;
            continue;
        }
         
    }

}


void debug_window_cigar(overlap_region_alloc* overlap_list, UC_Read* g_read, Correct_dumy* dumy,
All_reads* R_INF, int test_window, int test_boundary)
{
    long long i, j, y_id, y_strand;
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
                            fprintf(stderr, "j: %d, i: %d, y_id: %d, y_start: %d, y_end: %d\n", j, i, y_id, y_start, y_end);
                        }
                    }
                }
            }
            

            
            if(test_boundary == 1)
            {
                for (i = 0; i < overlap_list->list[j].boundary_cigars.length; i++)
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
                            fprintf(stderr, "j: %d, i: %d, y_id: %d, y_start: %d, y_end: %d\n", j, i, y_id, y_start, y_end);
                        }
                    }
                }
            }


            if(test_window == 1 && test_boundary == 1)
            {
                if(overlap_list->list[j].w_list_length != 
                        overlap_list->list[j].boundary_cigars.length + 1)
                {
                    fprintf(stderr, "w_list_length: %d, boundary_cigars.length: %d\n", 
                    overlap_list->list[j].w_list_length, overlap_list->list[j].boundary_cigars.length);
                }
            }
            
        }
    }
}

inline void recalcate_window_advance(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read)
{
    long long j, k, i;
    long long Len_x;
    int threshold;
    long long y_len;
    long long currentIDLen;
    long long matches;
    long long y_id;
    int y_strand;
    long long y_readLen;
    long long x_start;
    long long x_end;
    long long x_len;
    long long total_y_start;
    long long total_y_end;
    long long y_start;
    long long y_end;
    long long Window_Len;
    char* x_string;
    char* y_string;
    int end_site;
    unsigned int error;
    int real_y_start;
    long long overlap_length;
    int extra_begin, extra_end;
    long long o_len;


    ///j负责遍历整个overlap list
    for (j = 0; j < overlap_list->length; j++)
    {
        
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
                for (k = i + 1; k < overlap_list->list[j].w_list_length && overlap_list->list[j].w_list[k].y_end == -1; k++)
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

                    if(!determine_overlap_region(threshold, y_start, y_id, Window_Len, R_INF, 
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
                    
                    ///error等于-1说明没匹配
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
        for (i = 0; i < overlap_list->list[j].w_list_length; i++)
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
                    /****************************may have bugs********************************/
                    threshold = overlap_list->list[j].w_list[i].error_threshold;
                    /****************************may have bugs********************************/
                    /****************************may have bugs********************************/
                    ///should not adjust threshold, since this window can be matched by the old threshold
                    ///threshold = Adjust_Threshold(threshold, x_len);
                    /****************************may have bugs********************************/
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

                        // if(real_y_start < extra_begin || end_site >= Window_Len - extra_end)
                        // {
                        //     fprintf(stderr, "\nreal_y_start: %d, extra_begin: %d\n", 
                        //     real_y_start, extra_begin);

                        //     fprintf(stderr, "end_site: %d, Window_Len: %d, extra_end: %d\n", 
                        //     end_site, Window_Len, extra_end);                            
                        // }
                        
                        
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
                ///k遍历匹配window左侧所有不匹配的window
                ///如果i匹配，则k从i-1开始
                ///直到第一个匹配的window结束
                ///因为i!=0，所以k的大小不用担心
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
                    if(!determine_overlap_region(threshold, y_start, y_id, Window_Len, R_INF,
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

                        // if(real_y_start < extra_begin || end_site >= Window_Len - extra_end)
                        // {
                        //     fprintf(stderr, "\nreal_y_start: %d, extra_begin: %d\n", 
                        //     real_y_start, extra_begin);

                        //     fprintf(stderr, "end_site: %d, Window_Len: %d, extra_end: %d\n", 
                        //     end_site, Window_Len, extra_end);                            
                        // }

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
    
    int pre_threshold;
    long long tLen, tError;
    double error_rate;
    for (j = 0; j < overlap_list->length; j++)
    {
        y_id = overlap_list->list[j].y_id;
        y_strand = overlap_list->list[j].y_pos_strand;
        y_readLen = Get_READ_LENGTH((*R_INF), y_id);
        overlap_length = overlap_list->list[j].x_pos_e - overlap_list->list[j].x_pos_s + 1;
        overlap_list->list[j].is_match = 0;

        ///debug_scan_cigar(&(overlap_list->list[j]));

        ///only calculate cigar for high quality overlaps
        if (overlap_length * OVERLAP_THRESHOLD_FILTER <=  overlap_list->list[j].align_length)
        {
            
            for (i = 0; i < overlap_list->list[j].w_list_length; i++)
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
                        /****************************may have bugs********************************/
                        ///threshold = x_len * THRESHOLD_RATE;
                        threshold = overlap_list->list[j].w_list[i].error_threshold;
                        /****************************may have bugs********************************/
                        /****************************may have bugs********************************/
                        ///should not adjust threshold, since this window can be matched by the old threshold
                        ///threshold = Adjust_Threshold(threshold, x_len);
                        /****************************may have bugs********************************/
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

                        // if(error != overlap_list->list[j].w_list[i].error)
                        // {
                        //     fprintf(stderr, "error\n");
                        // }

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


                            // if(real_y_start < extra_begin || end_site >= Window_Len - extra_end)
                            // {
                            //     fprintf(stderr, "\nreal_y_start: %d, extra_begin: %d\n", 
                            //     real_y_start, extra_begin);

                            //     fprintf(stderr, "end_site: %d, Window_Len: %d, extra_end: %d\n", 
                            //     end_site, Window_Len, extra_end);                            
                            // }

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

            ///error_rate = trim_error_rate(overlap_list, j);
            error_rate = non_trim_error_rate(overlap_list, j, R_INF, dumy, g_read);
            
            
            ///if(error_rate <= 0.015)
            ///if(error_rate <= 0.03)
            if(error_rate <= FINAL_OVERLAP_ERROR_RATE)
            {
                overlap_list->mapped_overlaps_length += overlap_length;
                overlap_list->list[j].is_match = 1;
                calculate_boundary_cigars(&(overlap_list->list[j]), R_INF, dumy, g_read);
            }
            else if(error_rate <= 0.045)
            {
                overlap_list->list[j].is_match = 3;
            }
        }
    }

    ///debug_window_cigar(overlap_list, g_read, dumy, R_INF, 1, 1);
}





/**
inline void adjust_alignment_windows(overlap_region* overlap, All_reads* R_INF)
{
    long long y_id = overlap->y_id;
    long long y_strand = overlap->y_pos_strand;
    long long y_readLen = Get_READ_LENGTH((*R_INF), y_id);
    long long overlap_length = overlap->x_pos_e - overlap->x_pos_s + 1;
    long long i;
    ///note here we start from i = 1, instead of i = 0
    for (i = 1; i < overlap->w_list_length; i++)
    {
        ///that means we have both start pos and end pos 
        if(overlap->w_list[i].y_end != -1 && overlap->w_list[i].cigar.length != -1)
        {
            ;
        }
    }
}
**/



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
    ///deletion就不要管
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




///返回下一个backbone节点上的ID
long long inline add_path_to_correct_read(Graph* backbone, Correct_dumy* dumy, long long currentNodeID, 
long long type, long long edgeID, Cigar_record* current_cigar, char* self_string)
{
    //long long i;
    long long nodeID;

    ///Note: currentNodeID must be a backbone node
    ///currentNodeID = 0 means a fake node
    ///currentNodeID = i means self_string[i - 1]
    ///包括匹配和误配两种情况
    if (type == MISMATCH)
    {
        ///这是match的情况
        if(backbone->g_nodes.list[currentNodeID].mismatch_edges.list[edgeID].length == 0)
        {
            nodeID = backbone->g_nodes.list[currentNodeID].mismatch_edges.list[edgeID].out_node;
            add_base_to_correct_read_directly(dumy, backbone->g_nodes.list[nodeID].base);
            ///match所以dumy->corrected_base不要+1

            ///nodeID = i means self_string[i - 1]
            ///add_cigar_record(self_string+nodeID-1, 1, current_cigar, 0);
            add_cigar_record(&(backbone->g_nodes.list[nodeID].base), 1, current_cigar, 0);

            /***********需要注释掉********* */
            if (nodeID != currentNodeID + 1)
            {
                fprintf(stderr, "error match\n");
            }
            /***********需要注释掉********* */

            return nodeID;
        }
        else ///这是mismatch的情况
        {
            

            nodeID = backbone->g_nodes.list[currentNodeID].mismatch_edges.list[edgeID].out_node;
            add_base_to_correct_read_directly(dumy, backbone->g_nodes.list[nodeID].base);
            dumy->corrected_base++;

            
            char merge_base = 0;
            merge_base = seq_nt6_table[(uint8_t)backbone->g_nodes.list[nodeID].base];
            merge_base = merge_base << 3;
            ///这种中间节点只有一个元素，所以直接list[0]
            nodeID = backbone->g_nodes.list[nodeID].mismatch_edges.list[0].out_node;
            merge_base = merge_base | seq_nt6_table[(uint8_t)backbone->g_nodes.list[nodeID].base];
            add_cigar_record(&merge_base, 1, current_cigar, 1);
            
            /**
            add_cigar_record(&(backbone->g_nodes.list[nodeID].base), 1, current_cigar, 1);
            nodeID = backbone->g_nodes.list[nodeID].mismatch_edges.list[0].out_node;
            **/

            /***********需要注释掉********* */
            if (nodeID != currentNodeID + 1)
            {
                fprintf(stderr, "error mismatch\n");
            }
            /***********需要注释掉********* */

            return nodeID;
        }
    }
    else if (type == DELETION)
    {
        nodeID = backbone->g_nodes.list[currentNodeID].deletion_edges.list[edgeID].out_node;
        dumy->corrected_base += nodeID - currentNodeID;

        // if(nodeID - currentNodeID != 1)
        // {
        //     fprintf(stderr, "error\n");
        // }

        ///currentNodeID = i means self_string[i - 1]
        add_cigar_record(self_string + currentNodeID, nodeID - currentNodeID, current_cigar, DELETION);



        /***********需要注释掉********* */
        if (!(nodeID >= backbone->s_start_nodeID && nodeID <= backbone->s_end_nodeID))
        {
            fprintf(stderr, "error deletion 1\n");
        }
        if (nodeID <= currentNodeID)
        {
            fprintf(stderr, "error deletion 2\n");
        }
        
        /***********需要注释掉********* */
        
        return nodeID;
    }
    else if (type == INSERTION)
    {
        ///这个一定要变成0
        backbone->g_nodes.list[currentNodeID].num_insertions = 0;

        nodeID = backbone->g_nodes.list[currentNodeID].insertion_edges.list[edgeID].out_node;
        long long step = backbone->g_nodes.list[currentNodeID].insertion_edges.list[edgeID].length;
        long long i;
        for (i = 0; i < step; i++)
        {
            add_base_to_correct_read_directly(dumy, backbone->g_nodes.list[nodeID].base);
            add_cigar_record(&backbone->g_nodes.list[nodeID].base, 1, current_cigar, INSERTION);
            ///只有一条边
            nodeID = backbone->g_nodes.list[nodeID].insertion_edges.list[0].out_node;
        }
        dumy->corrected_base += step;

        ///currentNodeID = i means self_string[i - 1]
        ///add_cigar_record(self_string + currentNodeID, step, current_cigar, INSERTION);

        /***********需要注释掉********* */
        if (nodeID != currentNodeID)
        {
            fprintf(stderr, "error insertion\n");
        }
        /***********需要注释掉********* */

        return nodeID;
    }
    else
    {
        fprintf(stderr, "error type\n");
    }
    
}


///返回下一个backbone节点上的ID
long long inline add_path_to_correct_read_new(Graph* backbone, Graph* DAGCon, Correct_dumy* dumy, long long currentNodeID, 
long long type, long long edgeID, Cigar_record* current_cigar, char* self_string)
{
    //long long i;
    long long nodeID;

    ///Note: currentNodeID must be a backbone node
    ///currentNodeID = 0 means a fake node
    ///currentNodeID = i means self_string[i - 1]
    ///包括匹配和误配两种情况
    if (type == MISMATCH)
    {
        ///这是match的情况
        if(backbone->g_nodes.list[currentNodeID].mismatch_edges.list[edgeID].length == 0)
        {
            nodeID = backbone->g_nodes.list[currentNodeID].mismatch_edges.list[edgeID].out_node;
            add_base_to_correct_read_directly(dumy, backbone->g_nodes.list[nodeID].base);
            ///match所以dumy->corrected_base不要+1

            ///nodeID = i means self_string[i - 1]
            ///add_cigar_record(self_string+nodeID-1, 1, current_cigar, 0);
            add_cigar_record(&(backbone->g_nodes.list[nodeID].base), 1, current_cigar, 0);

            /***********需要注释掉********* */
            if (nodeID != currentNodeID + 1)
            {
                fprintf(stderr, "error match\n");
            }
            /***********需要注释掉********* */

            return nodeID;
        }
        else ///这是mismatch的情况
        {
            

            nodeID = backbone->g_nodes.list[currentNodeID].mismatch_edges.list[edgeID].out_node;
            add_base_to_correct_read_directly(dumy, backbone->g_nodes.list[nodeID].base);
            dumy->corrected_base++;

            
            char merge_base = 0;
            merge_base = seq_nt6_table[(uint8_t)backbone->g_nodes.list[nodeID].base];
            merge_base = merge_base << 3;
            ///这种中间节点只有一个元素，所以直接list[0]
            nodeID = backbone->g_nodes.list[nodeID].mismatch_edges.list[0].out_node;
            merge_base = merge_base | seq_nt6_table[(uint8_t)backbone->g_nodes.list[nodeID].base];
            add_cigar_record(&merge_base, 1, current_cigar, 1);
            
            /**
            add_cigar_record(&(backbone->g_nodes.list[nodeID].base), 1, current_cigar, 1);
            nodeID = backbone->g_nodes.list[nodeID].mismatch_edges.list[0].out_node;
            **/

            /***********需要注释掉********* */
            if (nodeID != currentNodeID + 1)
            {
                fprintf(stderr, "error mismatch\n");
            }
            /***********需要注释掉********* */

            return nodeID;
        }
    }
    else if (type == DELETION)
    {
        nodeID = backbone->g_nodes.list[currentNodeID].deletion_edges.list[edgeID].out_node;
        dumy->corrected_base += nodeID - currentNodeID;

        // if(nodeID - currentNodeID != 1)
        // {
        //     fprintf(stderr, "error\n");
        // }

        ///currentNodeID = i means self_string[i - 1]
        add_cigar_record(self_string + currentNodeID, nodeID - currentNodeID, current_cigar, DELETION);



        /***********需要注释掉********* */
        if (!(nodeID >= backbone->s_start_nodeID && nodeID <= backbone->s_end_nodeID))
        {
            fprintf(stderr, "error deletion 1\n");
        }
        if (nodeID <= currentNodeID)
        {
            fprintf(stderr, "error deletion 2\n");
        }
        
        /***********需要注释掉********* */
        
        return nodeID;
    }
    else if (type == INSERTION)
    {
        ///这个一定要变成0
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

        // nodeID = backbone->g_nodes.list[currentNodeID].insertion_edges.list[edgeID].out_node;
        // long long step = backbone->g_nodes.list[currentNodeID].insertion_edges.list[edgeID].length;
        // long long i;
        // for (i = 0; i < step; i++)
        // {
        //     add_base_to_correct_read_directly(dumy, backbone->g_nodes.list[nodeID].base);
        //     add_cigar_record(&backbone->g_nodes.list[nodeID].base, 1, current_cigar, INSERTION);
        //     ///只有一条边
        //     nodeID = backbone->g_nodes.list[nodeID].insertion_edges.list[0].out_node;
        // }
        // dumy->corrected_base += step;
        // /***********需要注释掉********* */
        // if (nodeID != currentNodeID)
        // {
        //     fprintf(stderr, "error insertion\n");
        // }
        // /***********需要注释掉********* */
        // return nodeID;
    }
    else
    {
        fprintf(stderr, "error type\n");
    }
    
}



void test_single_path(Graph* DAGCon, Graph* backbone, int debug_node_in_backbone)
{
    char forward[1000];
    char reverse[1000];
    char pre[1000];
    long long i, j, outNode, inputNode, preNode, string_i, path_weight, step;

    if((Output_Edges(G_Node(*DAGCon, DAGCon->s_start_nodeID)).length != Input_Edges(G_Node(*DAGCon, DAGCon->s_end_nodeID)).length)
    || (G_Node(*backbone, debug_node_in_backbone).insertion_edges.length != Output_Edges(G_Node(*DAGCon, DAGCon->s_start_nodeID)).length))
    {
        fprintf(stderr, "s_out: %d, s_end: %d\n",
        Output_Edges(G_Node(*DAGCon, DAGCon->s_start_nodeID)).length,
        Input_Edges(G_Node(*DAGCon, DAGCon->s_end_nodeID)).length);
    }



    

    for (i = 0; i < Output_Edges(G_Node(*DAGCon, DAGCon->s_start_nodeID)).length; i++)
    {
        outNode = Output_Edges(G_Node(*DAGCon, DAGCon->s_start_nodeID)).list[i].out_node;
        string_i = 0;
        path_weight = Output_Edges(G_Node(*DAGCon, DAGCon->s_start_nodeID)).list[i].weight;

        while(outNode != DAGCon->s_end_nodeID)
        {
            forward[string_i++] = G_Node(*DAGCon, outNode).base;
            if(Output_Edges(G_Node(*DAGCon, outNode)).list[0].weight != path_weight)
            {
                fprintf(stderr, "error1\n");
            }
            outNode = Output_Edges(G_Node(*DAGCon, outNode)).list[0].out_node;
        }
        forward[string_i] = '\0';


        inputNode = Input_Edges(G_Node(*DAGCon, DAGCon->s_end_nodeID)).list[i].in_node;
        string_i = 0;
        while(inputNode != DAGCon->s_start_nodeID)
        {
            reverse[string_i++] = G_Node(*DAGCon, inputNode).base;



            if(Input_Edges(G_Node(*DAGCon, inputNode)).list[0].weight != path_weight)
            {
                fprintf(stderr, "error2\n");
            }
            inputNode = Input_Edges(G_Node(*DAGCon, inputNode)).list[0].in_node;
        }
        reverse[string_i] = '\0';

        for(j = 0; j < string_i/2; j ++)
        {
            char k = reverse[j];
            reverse[j] = reverse[string_i - j - 1];
            reverse[string_i - j - 1] = k;
        }

        if(strcmp(forward, reverse))
        {
            fprintf(stderr, "f: %s\n, r: %s\n\n", forward, reverse);
        }




        if(G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[i].weight != path_weight)
        {
            fprintf(stderr, "error3\n");
        }

        if(G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[i].length != string_i)
        {
            fprintf(stderr, "error4\n");
        }

        step = G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[i].length;
        if(step != 0)
        {
            string_i = 0;
            preNode = G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[i].out_node;
        
            for (j = 0; j < step; j++)
            {
                pre[string_i++] = G_Node(*backbone, preNode).base;
                preNode = G_Node(*backbone, preNode).insertion_edges.list[0].out_node;
            }
        }

        pre[string_i] = '\0';


        if(strcmp(forward, pre))
        {
            fprintf(stderr, "f: %s, r: %s, p: %s\n\n", forward, reverse, pre);
        }

    }
}




void test_single_path_new(Graph* DAGCon, Graph* backbone, int debug_node_in_backbone)
{
    char forward[1000];
    char reverse[1000];
    char pre[1000];
    long long i, j, preNode, string_i, path_weight, step;
    Node* outNode;
    Node* inputNode;
    Node* currentStartNode;
    Node* currentEndNode;

    if(
        (Real_Length(Output_Edges(G_Node(*DAGCon, DAGCon->s_start_nodeID))) 
        != 
        Real_Length(Input_Edges(G_Node(*DAGCon, DAGCon->s_end_nodeID))))
        || 
        (G_Node(*backbone, debug_node_in_backbone).insertion_edges.length
        != 
        Real_Length(Output_Edges(G_Node(*DAGCon, DAGCon->s_start_nodeID)))))
    {
        fprintf(stderr, "s_out: %d, s_end: %d, insertion_edges.length: %d\n",
        Real_Length(Output_Edges(G_Node(*DAGCon, DAGCon->s_start_nodeID))),
        Real_Length(Input_Edges(G_Node(*DAGCon, DAGCon->s_end_nodeID))),
        G_Node(*backbone, debug_node_in_backbone).insertion_edges.length);
    }

    RSet iter_out, iter_input;
    clear_RSet(&iter_out);
    clear_RSet(&iter_input);
    currentStartNode = &(G_Node(*DAGCon, DAGCon->s_start_nodeID));
    currentEndNode = &(G_Node(*DAGCon, DAGCon->s_end_nodeID));

    i = 0;



    while(getOutputNodes(&iter_out, DAGCon, currentStartNode, &outNode))
    {
        string_i = 0;
        path_weight = Input_Edges((*outNode)).list[0].weight;

        while(outNode != &(G_Node(*DAGCon, DAGCon->s_end_nodeID)))
        {
            RSet inner;
            Edge* e;
            forward[string_i++] = (*outNode).base;

            clear_RSet(&inner);
            while (getOutputEdges(&inner, DAGCon, outNode, &e))
            {
                if(e->weight != path_weight)
                {
                    fprintf(stderr, "error1\n");
                }
            }


            clear_RSet(&inner);
            while (getInputEdges(&inner, DAGCon, outNode, &e))
            {
                if(e->weight != path_weight)
                {
                    fprintf(stderr, "error1\n");
                }
            }
    
            clear_RSet(&inner);
            while(getOutputNodes(&inner, DAGCon, outNode, &outNode))
            {
                ;
            }
        }
        forward[string_i] = '\0';

        
        string_i = 0;
        if(!getInputNodes(&iter_input, DAGCon, currentEndNode, &inputNode))
        {
            fprintf(stderr, "sbsbsb\n");
        }

        while(inputNode != &(G_Node(*DAGCon,DAGCon->s_start_nodeID)))
        {
            reverse[string_i++] = (*inputNode).base;
            if(Input_Edges((*inputNode)).list[0].weight != path_weight)
            {
                fprintf(stderr, "error2\n");
            }
            RSet inner;
            Edge* e;
            clear_RSet(&inner);
            while (getOutputEdges(&inner, DAGCon, inputNode, &e))
            {
                if(e->weight != path_weight)
                {
                    fprintf(stderr, "error1\n");
                }
            }
            clear_RSet(&inner);
            while (getInputEdges(&inner, DAGCon, inputNode, &e))
            {
                if(e->weight != path_weight)
                {
                    fprintf(stderr, "error1\n");
                }
            }

            clear_RSet(&inner);
            while(getInputNodes(&inner, DAGCon, inputNode, &inputNode))
            {
                ;
            }

            ///inputNode = &(G_Node(*DAGCon, Input_Edges((*inputNode)).list[0].in_node));
        }
        reverse[string_i] = '\0';

        for(j = 0; j < string_i/2; j ++)
        {
            char k = reverse[j];
            reverse[j] = reverse[string_i - j - 1];
            reverse[string_i - j - 1] = k;
        }

        if(strcmp(forward, reverse)!=0)
        {
            fprintf(stderr, "f: %s, r: %s\n\n", forward, reverse);
        }





        if(G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[i].weight != path_weight)
        {
            fprintf(stderr, "error3\n");
        }

        if(G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[i].length != string_i)
        {
            fprintf(stderr, "error4\n");
        }

        /**
        step = G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[i].length;
        if(step != 0)
        {
            string_i = 0;
            preNode = G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[i].out_node;
        
            for (j = 0; j < step; j++)
            {
                pre[string_i++] = G_Node(*backbone, preNode).base;
                preNode = G_Node(*backbone, preNode).insertion_edges.list[0].out_node;
            }
        }

        pre[string_i] = '\0';
        **/

       extract_path(backbone, debug_node_in_backbone, i, pre);


        if(strcmp(forward, pre)!=0)
        {
            fprintf(stderr, "f: %s, r: %s, p: %s\n\n", forward, reverse, pre);
        }
        

        i++;
    }

    int pre_weight = 0;
    for (i = 0; i < G_Node(*backbone, debug_node_in_backbone).insertion_edges.length; i++)
    {
        if(i == 0)
        {
            extract_path(backbone, debug_node_in_backbone, i, pre);
            pre_weight = G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[i].weight;
        }
        else if(i > 0)
        {
            extract_path(backbone, debug_node_in_backbone, i, forward);

            if(strcmp(forward, pre)==0)
            {
                fprintf(stderr, "f: %s, f_weight: %d, p: %s, p_weight: %d\n\n", 
                forward, G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[i].weight,
                pre, pre_weight);
            }

            memcpy(pre, forward, strlen(forward) + 1);
            
        }
    }
    

    
    for (i = 0; i < DAGCon->g_nodes.length; i++)
    {
        currentStartNode = &(G_Node(*DAGCon, i));

        if(If_Node_Exist(*currentStartNode))
        {
            clear_RSet(&iter_out);
            Edge* e;
            Edge* e_self;
            Edge* e_reverse;
            while(getOutputEdges(&iter_out, DAGCon, currentStartNode, &e_self))
            {
                ///e_reverse = &(Input_Edges(G_Node(*DAGCon, e_self->out_node)).list[e_self->reverse_edge_ID]);
                e = e_self;

                get_bi_direction_edges(DAGCon, e_self, &e_self, &e_reverse);

                if(e != e_self)
                {
                    fprintf(stderr, "error0\n");
                }

                if(e_self->in_node != e_reverse->in_node || e_self->out_node != e_reverse->out_node || 
                e_self->weight != e_reverse->weight)
                {
                    fprintf(stderr, "error1\n");
                }


                if(e_self != &(Output_Edges(G_Node(*DAGCon, e_self->in_node)).list[e_self->self_edge_ID]))
                {
                    fprintf(stderr, "error2\n");
                }

                if(Visit(*e_self) == 1 || Visit(*e_reverse) == 1)
                {
                    fprintf(stderr, "error visit flag\n");
                }







                get_bi_direction_edges(DAGCon, e_reverse, &e_self, &e_reverse);

                if(e != e_self)
                {
                    fprintf(stderr, "error0\n");
                }

               if(e_self->in_node != e_reverse->in_node || e_self->out_node != e_reverse->out_node || 
                e_self->weight != e_reverse->weight)
                {
                    fprintf(stderr, "error1\n");
                }


                if(e_self != &(Output_Edges(G_Node(*DAGCon, e_self->in_node)).list[e_self->self_edge_ID]))
                {
                    fprintf(stderr, "error2\n");
                }
            }



            clear_RSet(&iter_input);
            while(getInputEdges(&iter_input, DAGCon, currentStartNode, &e_self))
            {

                e = e_self;

                get_bi_direction_edges(DAGCon, e_self, &e_self, &e_reverse);

                if(e != e_reverse)
                {
                    fprintf(stderr, "error0\n");
                }

                if(e_self->in_node != e_reverse->in_node || e_self->out_node != e_reverse->out_node || 
                e_self->weight != e_reverse->weight)
                {
                    fprintf(stderr, "error1\n");
                }


                if(e_self != &(Output_Edges(G_Node(*DAGCon, e_self->in_node)).list[e_self->self_edge_ID]))
                {
                    fprintf(stderr, "error2\n");
                }







                get_bi_direction_edges(DAGCon, e_reverse, &e_self, &e_reverse);

                if(e != e_reverse)
                {
                    fprintf(stderr, "error0\n");
                }

               if(e_self->in_node != e_reverse->in_node || e_self->out_node != e_reverse->out_node || 
                e_self->weight != e_reverse->weight)
                {
                    fprintf(stderr, "error1\n");
                }


                if(e_self != &(Output_Edges(G_Node(*DAGCon, e_self->in_node)).list[e_self->self_edge_ID]))
                {
                    fprintf(stderr, "error2\n");
                }
            }
        }
        // else
        // {
        //     fprintf(stderr, "node does not exist\n");
        // }
        
    }
}


void debug_Queue(Graph* DAGCon)
{
    long long* input;
    long long* output;

    srand((unsigned)time(0)); 
    long long array_length = rand() % 100;
    input = (long long*)malloc(sizeof(long long) * (array_length + 1));
    output = (long long*)malloc(sizeof(long long) * (array_length + 1));
    long long i;

    for (i = 0; i < array_length; i++)
    {
       input[i] = rand() % 1000000;
       push_to_Queue(&(DAGCon->node_q), input[i]);
    }

    i = 0;
    while (pop_from_Queue(&(DAGCon->node_q), &output[i]))
    {
        i++;
    }

    if(i != array_length)
    {
        fprintf(stderr, "array_length: %d\n", array_length);
    }
    else
    {
        for (i = 0; i < array_length; i++)
        {
            if(input[i] != output[i])
            {
                fprintf(stderr, "input[%d]:%d, output[%d]: %d\n", 
            i, input[i], i, output[i]);
            }
            
        }
    }
    
    

 
    
        
        
    

    srand((unsigned)time(0));
    if(array_length != 0)
    {
        array_length = rand() % array_length;
    }
    

    for (i = 0; i < array_length; i++)
    {
       input[i] = rand() % 1000000;
       push_to_Queue(&(DAGCon->node_q), input[i]);
    }

    i = 0;
    while (pop_from_Queue(&(DAGCon->node_q), &output[i]))
    {
        i++;
    }

    if(i != array_length)
    {
        fprintf(stderr, "array_length: %d\n", array_length);
    }
    else
    {
        for (i = 0; i < array_length; i++)
        {
            if(input[i] != output[i])
            {
                fprintf(stderr, "input[%d]:%d, output[%d]: %d\n", 
            i, input[i], i, output[i]);
            }
            
        }
    }


    free(input);
    free(output);
}


///return the in-edge ID of outNode
long long get_In_Edge_ID(Graph* DAGCon, long long inNode, long long outNode)
{
    long long i;

    for (i = 0; i < Input_Edges(G_Node(*DAGCon, outNode)).length; i++)
    {
        if (Input_Edges(G_Node(*DAGCon, outNode)).list[i].in_node == inNode)
        {
            return i;
        }
    }

    return -1;
}


///return the out-edge ID of inNode
long long get_Out_Edge_ID(Graph* DAGCon, long long inNode, long long outNode)
{
    long long i;
    for (i = 0; i < Output_Edges(G_Node(*DAGCon, inNode)).length; i++)
    {
        if(Output_Edges(G_Node(*DAGCon, inNode)).list[i].out_node == outNode)
        {
            return i;
        }
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
    long long base_i, i, weight;
    int flag = 0;
    Node* get_node_1;
    Node* out_node_of_get_node_1;
    Node* consensus_node_1;
    Edge* e_forward_1;
    Edge* e_backward_1;

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
    long long base_i, i, weight;
    int flag = 0;
    Node* get_node;
    Node* in_node_of_get_node;
    Node* consensus_node;
    Edge* e_forward;
    Edge* e_backward;

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
    long long i;
    for (long long i = 0; i < DAGCon->g_nodes.length; i++)
    {
        Node* currentStartNode = &(G_Node(*DAGCon, i));
        RSet iter_out;

        if(If_Node_Exist(*currentStartNode))
        {
            fprintf(stderr, "ID: %d (%c) (w: %d)\n", (*currentStartNode).ID, (*currentStartNode).base, (*currentStartNode).weight);
            clear_RSet(&iter_out);
            Edge* e;
            fprintf(stderr, "****Out-node: ");
            while(getOutputEdges(&iter_out, DAGCon, currentStartNode, &e))
            {
                //fprintf(stderr, "%d[%c], ", G_Node(*DAGCon, e->out_node).ID, G_Node(*DAGCon, e->out_node).base);
                fprintf(stderr, "%d(w: %d), ", G_Node(*DAGCon, e->out_node).ID, e->weight);
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
    for (long long i = 0; i < DAGCon->g_nodes.length; i++)
    {
        Node* currentStartNode = &(G_Node(*DAGCon, i));
        RSet iter_out;

        if(If_Node_Exist(*currentStartNode))
        {
            clear_RSet(&iter_out);
            Edge* e;
            Edge* e_self;
            Edge* e_reverse;
            while(getOutputEdges(&iter_out, DAGCon, currentStartNode, &e_self))
            {

                get_bi_direction_edges(DAGCon, e_self, &e_self, &e_reverse);


                if(Visit(*e_self) == 0)
                {
                    fprintf(stderr, "Visit(*e_self): %d, error visit flag: in_node: %d, out_node: %d\n", 
                    Visit(*e_self), (*e_self).in_node, (*e_self).out_node);
                }


                if(Visit(*e_reverse) == 0)
                {
                    fprintf(stderr, "Visit(*e_reverse): %d, error visit flag: in_node: %d, out_node: %d\n", 
                    Visit(*e_reverse), (*e_reverse).in_node, (*e_reverse).out_node);
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
                    fprintf(stderr, "Visit(*e_self): %d, error visit flag: in_node: %d, out_node: %d\n", 
                    Visit(*e_self), (*e_self).in_node, (*e_self).out_node);
                }


                if(Visit(*e_reverse) == 0)
                {
                    fprintf(stderr, "Visit(*e_reverse): %d, error visit flag: in_node: %d, out_node: %d\n", 
                    Visit(*e_reverse), (*e_reverse).in_node, (*e_reverse).out_node);
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
    Edge* e;
    long long max;
    Node* max_node;
    

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
    long long max_start, max_end, max_start_edge, max_end_edge;
    RSet iter;
    Edge* e;
    Node* newNode;
    long long max_count;
    

    ///check the out-edges of start node
    ///must to be 0
    max_start = 0;
    newNode = &(G_Node(*DAGCon, DAGCon->s_start_nodeID));
    clear_RSet(&iter);
    while(getOutputEdges(&iter, DAGCon, newNode, &e))
    {
        if(e->weight > max_start)
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
        if(e->weight > max_end)
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

    
    
    ///if((*direction) == 1)
    ///if(max_start < max_end && DAGCon->g_nodes.length > 5)
    ///if(DAGCon->g_nodes.length > 5)
    // if(max_start < max_end && DAGCon->node_q.end - DAGCon->node_q.beg > 1) 
    // {
    //     print_graph(DAGCon);
    //     long long str;
    //     while (pop_from_Queue(&(DAGCon->node_q), &str))
    //     {
    //         fprintf(stderr, "%c", (char)str);
    //     }
    //     fprintf(stderr, "\n");
    //     fprintf(stderr, "###################(*max_count): %d, (*max_edge): %d, (*direction): %d###################\n\n", 
    //     (*max_count), max_start >= max_end? max_start_edge:max_end_edge, max_start >= max_end? 0:1);
    // }
    // if(max_start >= max_end)
    // {
    //     newNode = &(G_Node(*DAGCon, DAGCon->s_start_nodeID));
    //     if(Output_Edges(*newNode).list[max_start_edge].weight != (*max_count))
    //     {
    //         fprintf(stderr, "ERROR\n");
    //     }
    // }
    // else
    // {
    //     newNode = &(G_Node(*DAGCon, DAGCon->s_end_nodeID));
    //     if(Input_Edges(*newNode).list[max_end_edge].weight != (*max_count))
    //     {
    //         fprintf(stderr, "ERROR\n");
    //     }
    // }

    return max_count;


}


inline void generate_seq_from_node(Graph* DAGCon, Node* node, int direction)
{
    clear_Queue(&(DAGCon->node_q));
    RSet iter;
    long long max;
    Node* max_node;
    Node* getNodes;
    

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
    long long max_start, max_end;
    RSet iter;
    Edge* e;
    Node* newNode;
    Node* getNode;
    Node* max_start_node;
    Node* max_end_node;
    long long max_count, i;

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
        if(getNode->weight > max_start)
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
        if(getNode->weight > max_end)
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

    
    /**
    if(DAGCon->g_nodes.length > 5)
    ///if(max_start < max_end && DAGCon->node_q.end - DAGCon->node_q.beg > 1) 
    ///if(max_start < max_end)
    {
        print_graph(DAGCon);
        long long str;
        while (pop_from_Queue(&(DAGCon->node_q), &str))
        {
            fprintf(stderr, "%c", (char)str);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "###################(*max_count): %d, (*max_node): %d, (*direction): %d###################\n\n", 
        max_count, max_start >= max_end? max_start_node->ID:max_end_node->ID, max_start >= max_end? 0:1);
    }
    if(max_start >= max_end)
    {
        newNode = &(G_Node(*DAGCon, DAGCon->s_start_nodeID));
        if(max_start_node->weight != max_count)
        {
            fprintf(stderr, "ERROR\n");
        }
    }
    else
    {
        newNode = &(G_Node(*DAGCon, DAGCon->s_end_nodeID));
        if(max_end_node->weight != max_count)
        {
            fprintf(stderr, "ERROR\n");
        }
    }
    **/

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


    for (i = 0; i < G_Node(*backbone, currentNodeID).insertion_edges.length; i++)
    {
        path_weight = G_Node(*backbone, currentNodeID).insertion_edges.list[i].weight;
        
        lastNode = &(G_Node(*DAGCon, DAGCon->s_start_nodeID));


        /*****************************debug***************************************/
        // newNode = add_Node_DAGCon(DAGCon, 'F');
        // add_bi_direction_edge(DAGCon, lastNode, newNode, path_weight, 0);
        // add_bi_direction_edge(DAGCon, newNode, &(G_Node(*DAGCon, DAGCon->s_end_nodeID)), path_weight, 0);
        // delete_Node_DAGCon(DAGCon, newNode);
        // if(remove_and_check_bi_direction_edge_from_nodes(DAGCon, &(G_Node(*DAGCon, DAGCon->s_start_nodeID)), newNode))
        // {
        //     fprintf(stderr, "step: %d, j: %d\n", step, j);
        // }
        // if(remove_and_check_bi_direction_edge_from_nodes(DAGCon, newNode, &(G_Node(*DAGCon, DAGCon->s_end_nodeID))))
        // {
        //     fprintf(stderr, "step: %d, j: %d\n", step, j);
        // }


        // Node* node0;
        // Node* node1;
        // Node* node2;
        // newNode = add_Node_DAGCon(DAGCon, 'T');
        // add_bi_direction_edge(DAGCon, lastNode, newNode, path_weight, 0);
        // node0 = newNode;
        // lastNode = newNode;
        // newNode = add_Node_DAGCon(DAGCon, 'T');
        // add_bi_direction_edge(DAGCon, lastNode, newNode, path_weight, 0);
        // node1 = newNode;
        // lastNode = newNode;
        // newNode = add_Node_DAGCon(DAGCon, 'T');
        // add_bi_direction_edge(DAGCon, lastNode, newNode, path_weight, 0);
        // node2 = newNode;


        // lastNode = &(G_Node(*DAGCon, DAGCon->s_start_nodeID));
        // delete_Node_DAGCon(DAGCon, node1);

        // Edge* e_forward;
        // Edge* e_backward;
        // if(get_bi_Edge(DAGCon, node2, &(G_Node(*DAGCon, DAGCon->s_end_nodeID)), &e_forward, &e_backward))
        // {
        //     remove_and_check_bi_direction_edge_from_edge(DAGCon, e_forward);
        // }

        // if(get_bi_Edge(DAGCon, node2, &(G_Node(*DAGCon, DAGCon->s_end_nodeID)), &e_forward, &e_backward))
        // {
        //     fprintf(stderr, "edge remove error\n");
        // }

        // remove_and_check_bi_direction_edge_from_nodes(DAGCon, lastNode, node0);
        // remove_and_check_bi_direction_edge_from_nodes(DAGCon, node0, node1);
        // remove_and_check_bi_direction_edge_from_nodes(DAGCon, node1, node2);
        // if(remove_and_check_bi_direction_edge_from_nodes(DAGCon, node2, &(G_Node(*DAGCon, DAGCon->s_end_nodeID))))
        // {
        //     fprintf(stderr, "edge remove error\n");
        // }
        /*****************************debug***************************************/


        step = G_Node(*backbone, currentNodeID).insertion_edges.list[i].length;
        if(step != 0)
        {
            nodeID = G_Node(*backbone, currentNodeID).insertion_edges.list[i].out_node;
        
            for (j = 0; j < step; j++)
            {
                base = G_Node(*backbone, nodeID).base;
                newNode = add_Node_DAGCon(DAGCon, base);
                add_bi_direction_edge(DAGCon, lastNode, newNode, path_weight, 0);

                /*****************************debug***************************************/
                // if(!add_and_check_bi_direction_edge(DAGCon, lastNode, newNode, path_weight, 0))
                // {
                //     fprintf(stderr, "haha\n");
                // }

                // if(add_and_check_bi_direction_edge(DAGCon, lastNode, newNode, path_weight, 0))
                // {
                //     fprintf(stderr, "haha\n");
                // }
                // if(j != 0)
                // {
                //     add_bi_direction_edge(DAGCon, &(G_Node(*DAGCon, DAGCon->s_start_nodeID)), newNode, path_weight, 0);
                //     if(!remove_and_check_bi_direction_edge_from_nodes(DAGCon, &(G_Node(*DAGCon, DAGCon->s_start_nodeID)), 
                //     newNode))
                //     {
                //         fprintf(stderr, "step: %d, j: %d\n", step, j);
                //     }
                // }
                /*****************************debug***************************************/

                nodeID = G_Node(*backbone, nodeID).insertion_edges.list[0].out_node;

                lastNode = newNode;
            }

            if(lastNode->ID != DAGCon->s_start_nodeID)
            {
                add_bi_direction_edge(DAGCon, lastNode, &(G_Node(*DAGCon, DAGCon->s_end_nodeID)), path_weight, 0);

                /*****************************debug***************************************/
                // if(!add_and_check_bi_direction_edge(DAGCon, lastNode, &(G_Node(*DAGCon, DAGCon->s_end_nodeID)), path_weight, 0))
                // {
                //     fprintf(stderr, "haha\n");
                // }

                // if(add_and_check_bi_direction_edge(DAGCon, lastNode, &(G_Node(*DAGCon, DAGCon->s_end_nodeID)), path_weight, 0))
                // {
                //     fprintf(stderr, "haha\n");
                // }
                // add_bi_direction_edge(DAGCon, &(G_Node(*DAGCon, DAGCon->s_start_nodeID)), 
                // &(G_Node(*DAGCon, DAGCon->s_end_nodeID)), path_weight, 0);
                // if(!remove_and_check_bi_direction_edge_from_nodes(DAGCon, &(G_Node(*DAGCon, DAGCon->s_start_nodeID)), 
                // &(G_Node(*DAGCon, DAGCon->s_end_nodeID))))
                // {
                //     fprintf(stderr, "step: %d, j: %d\n", step, j);
                // }
                /*****************************debug***************************************/
            }
        }


        
    }

    ///test_single_path_new(DAGCon, backbone, currentNodeID);


    Merge_DAGCon(DAGCon);

    ///(*max_count) = generate_best_seq_from_edges(DAGCon);
    (*max_count) = generate_best_seq_from_nodes(DAGCon);

    // long long k = 0;
    // for (i = 0; i < DAGCon->g_nodes.length; i++)
    // {
    //     if(i != DAGCon->s_start_nodeID && i != DAGCon->s_end_nodeID)
    //     {
    //         if(If_Node_Exist(G_Node(*DAGCon, i)) && G_Node(*DAGCon, i).weight > k)
    //         {
    //             k = G_Node(*DAGCon, i).weight;
    //         }
    //     }
    // }

    // if(k != (*max_count))
    // {
    //     fprintf(stderr, "k: %d, (*max_count): %d\n", k, (*max_count));
    // }
    
    
    ///very important
    backbone->g_nodes.list[currentNodeID].num_insertions = 0;

}

void debug_whole_graph(Graph* g)
{
    long long i, j, k;
    for (i = 0; i < g->g_nodes.length; i++)
    {
        if(g->g_nodes.list[i].deletion_edges.length!= 0 &&
        g->g_nodes.list[i].deletion_edges.length!= 1)
        {
            fprintf(stderr, "g->g_nodes.list[i].deletion_edges.length: %d\n",
            g->g_nodes.list[i].deletion_edges.length);
        }
    }

    for (i = g->s_start_nodeID; i < g->s_end_nodeID; i++)
    {
        if(g->g_nodes.list[i].mismatch_edges.length > 4 
        || 
        g->g_nodes.list[i].mismatch_edges.length < 1)
        {
            fprintf(stderr, "g->s_end_nodeID: %d, g->g_nodes.list[%d].mismatch_edges.length: %d\n",
            g->s_end_nodeID, i, g->g_nodes.list[i].mismatch_edges.length);
        }
    }

    char current[1000];
    char compare[1000];
    long long total_weight = 0;
    for (i = g->s_start_nodeID; i < g->s_end_nodeID; i++)
    {
        total_weight = 0;
        for (j = 0; j < G_Node(*g, i).insertion_edges.length; j++)
        {
            
            total_weight = total_weight + G_Node(*g, i).insertion_edges.list[j].weight;

            extract_path(g, i, j, current);

            for (k = j + 1; k < G_Node(*g, i).insertion_edges.length; k++)
            {
                extract_path(g, i, k, compare);
                if(strcmp(current, compare)==0)
                {
                    fprintf(stderr,"error\n");
                }
            }
        }

        if(total_weight != G_Node(*g, i).num_insertions)
        {
            fprintf(stderr,"error\n");
        }
    }

        

}

void get_seq_from_Graph(Graph* backbone, Graph* DAGCon, Correct_dumy* dumy, Cigar_record* current_cigar, char* self_string,
char* r_string, long long r_string_length, long long r_string_site)
{
    ///debug_whole_graph(backbone);

    // double threshold;
    // if(roundID > 0)
    // {
    //     threshold = CORRECT_THRESHOLD_SECOND;
    // }
    // else
    // {
    //     threshold = CORRECT_THRESHOLD;
    // }
    
    


    long long new_seq_length = 0;
    long long currentNodeID;
    long long i;
    // 总共有以下几种情况: 
    // 1. match 2. mismatch (A, C, G, T, N) 3. deletion 4. insertion (A, C, G, T)
    // 其实就是 1. 自己本身的weight 2. alignToNode的weight 3. insertion节点的weight
    long long max_count;
    int max_type;
    long long  max_edge;
    long long total_count;
    long long nodeID;
    char current_base;
    long long current_weight;

    long long max_insertion_count;

    currentNodeID = backbone->s_start_nodeID;

    ///fprintf(stderr, "currentNodeID: %d\n", currentNodeID);

    while (currentNodeID != backbone->s_end_nodeID)
    {
        total_count = 0;
        max_count = -1;
        max_type = -1;
        max_edge = -1;


        ///假如这是个backbone节点
        ///有三类出边
        ///1. mismatch_edges 2. insertion_edges 3. deletion_edges
        if (currentNodeID >= backbone->s_start_nodeID && currentNodeID <= backbone->s_end_nodeID)
        {
            
            
            
            

            ///mismatch_edges
            for (i = 0; i < backbone->g_nodes.list[currentNodeID].mismatch_edges.length; i++)
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

                ///match
                ///match要处理插入的情况
                ///如果这里有insertion, 这个节点会过两遍
                ///第一遍num_insertions > 0, 第二遍num_insertions=0
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
                // for (i = 0; i < backbone->g_nodes.list[currentNodeID].insertion_edges.length; i++)
                // {
                //     total_count = total_count + backbone->g_nodes.list[currentNodeID].insertion_edges.list[i].weight;

                //     if (backbone->g_nodes.list[currentNodeID].insertion_edges.list[i].weight > max_count)
                //     {
                //         max_count = backbone->g_nodes.list[currentNodeID].insertion_edges.list[i].weight;
                //         max_edge = i;
                //         max_type = INSERTION;
                //     }
                // }
            }

            
            
            

            ///deletion_edges
            for (i = 0; i < backbone->g_nodes.list[currentNodeID].deletion_edges.length; i++)
            {
                total_count = total_count + backbone->g_nodes.list[currentNodeID].deletion_edges.list[i].weight;

                if (backbone->g_nodes.list[currentNodeID].deletion_edges.list[i].weight > max_count)
                {
                    max_count = backbone->g_nodes.list[currentNodeID].deletion_edges.list[i].weight;
                    max_edge = i;
                    max_type = DELETION;
                }
            }



            
            
            ///这种情况下矫正
            if(max_count >= total_count*(CORRECT_THRESHOLD))
            ///if(max_count >= total_count * threshold)
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
                else///不矫正, 直接取下一个backbone节点
                {
                    currentNodeID++;
                    add_base_to_correct_read_directly(dumy, backbone->g_nodes.list[currentNodeID].base);

                    add_cigar_record(&(backbone->g_nodes.list[currentNodeID].base), 1, current_cigar, 0);
                }
            }
            
            ///fprintf(stderr, "currentNodeID: %d, max_type: %d\n", currentNodeID, max_type);
        }
        else  ///非backbone节点就会出错了
        {
            fprintf(stderr, "error\n");
        }
        
    }

}


void get_seq_from_Graph_print(Graph* backbone, Correct_dumy* dumy, Cigar_record* current_cigar, char* self_string,
char* r_string, long long r_string_length, long long r_string_site)
{
    long long new_seq_length = 0;
    long long currentNodeID;
    long long i;
    // 总共有以下几种情况: 
    // 1. match 2. mismatch (A, C, G, T, N) 3. deletion 4. insertion (A, C, G, T)
    // 其实就是 1. 自己本身的weight 2. alignToNode的weight 3. insertion节点的weight
    long long max_count;
    int max_type;
    long long  max_edge;
    long long total_count;
    long long nodeID;
    char current_base;
    long long current_weight;

    currentNodeID = backbone->s_start_nodeID;

    ///fprintf(stderr, "currentNodeID: %d\n", currentNodeID);

    while (currentNodeID != backbone->s_end_nodeID)
    {
        total_count = 0;
        max_count = -1;
        max_type = -1;
        max_edge = -1;


        ///假如这是个backbone节点
        ///有三类出边
        ///1. mismatch_edges 2. insertion_edges 3. deletion_edges
        if (currentNodeID >= backbone->s_start_nodeID && currentNodeID <= backbone->s_end_nodeID)
        {
            
            
            
            

            ///mismatch_edges
            for (i = 0; i < backbone->g_nodes.list[currentNodeID].mismatch_edges.length; i++)
            {
                if(currentNodeID == 187)
                {
                    fprintf(stderr, "backbone->g_nodes.list[currentNodeID].num_insertions: %d\n",
                    backbone->g_nodes.list[currentNodeID].num_insertions);
                }

                if (backbone->g_nodes.list[currentNodeID].num_insertions != 0)
                {
                    current_weight = backbone->g_nodes.list[currentNodeID].mismatch_edges.list[i].weight -
                        backbone->g_nodes.list[currentNodeID].mismatch_edges.list[i].num_insertions;
                }
                else
                {
                    current_weight = backbone->g_nodes.list[currentNodeID].mismatch_edges.list[i].weight;
                }

                if(currentNodeID == 187)
                {
                    fprintf(stderr, "current_weight: %d\n",
                    current_weight);
                }
                
                
                total_count = total_count + current_weight;

                ///match
                ///match要处理插入的情况
                ///如果这里有insertion, 这个节点会过两遍
                ///第一遍num_insertions > 0, 第二遍num_insertions=0
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
                for (i = 0; i < backbone->g_nodes.list[currentNodeID].insertion_edges.length; i++)
                {
                    total_count = total_count + backbone->g_nodes.list[currentNodeID].insertion_edges.list[i].weight;

                    if (backbone->g_nodes.list[currentNodeID].insertion_edges.list[i].weight > max_count)
                    {
                        max_count = backbone->g_nodes.list[currentNodeID].insertion_edges.list[i].weight;
                        max_edge = i;
                        max_type = INSERTION;
                    }

                    if(currentNodeID == 187)
                    {
                        fprintf(stderr, "backbone->g_nodes.list[currentNodeID].insertion_edges.list[i].weight: %d, length: %d\n",
                        backbone->g_nodes.list[currentNodeID].insertion_edges.list[i].weight,
                    backbone->g_nodes.list[currentNodeID].insertion_edges.list[i].length);

                    nodeID = backbone->g_nodes.list[currentNodeID].insertion_edges.list[i].out_node;
                    fprintf(stderr, "%c\n", backbone->g_nodes.list[nodeID].base);
                    }
                }
            }
            
            

            ///deletion_edges
            for (i = 0; i < backbone->g_nodes.list[currentNodeID].deletion_edges.length; i++)
            {
                total_count = total_count + backbone->g_nodes.list[currentNodeID].deletion_edges.list[i].weight;

                if (backbone->g_nodes.list[currentNodeID].deletion_edges.list[i].weight > max_count)
                {
                    max_count = backbone->g_nodes.list[currentNodeID].deletion_edges.list[i].weight;
                    max_edge = i;
                    max_type = DELETION;
                }

                if(currentNodeID == 187)
                {
                    fprintf(stderr, "backbone->g_nodes.list[currentNodeID].deletion_edges.list[i].weight: %d\n",
                    backbone->g_nodes.list[currentNodeID].deletion_edges.list[i].weight);
                }
            }

            /**
            if(currentNodeID > 0 && backbone->g_nodes.list[currentNodeID].base != r_string[r_string_site + currentNodeID - 1])
            {
                fprintf(stderr, "currentNodeID: %d\n", currentNodeID);
            }
            **/
            
            fprintf(stderr, "currentNodeID: %d, max_count: %d, max_type: %d, total_count: %d\n", 
            currentNodeID, max_count, max_type, total_count);



        if(currentNodeID == 187)
        {
            nodeID = backbone->g_nodes.list[currentNodeID].insertion_edges.list[max_edge].out_node;
            fprintf(stderr, "%c", backbone->g_nodes.list[nodeID].base);
            nodeID = backbone->g_nodes.list[nodeID].insertion_edges.list[0].out_node;
            fprintf(stderr, "%c\n", backbone->g_nodes.list[nodeID].base);
        }



            ///这种情况下矫正
            if(max_count >= total_count*CORRECT_THRESHOLD)
            {
                currentNodeID = add_path_to_correct_read(backbone, dumy, currentNodeID, max_type, max_edge, current_cigar,
                self_string);
            }
            else  
            {
                ///NOTE: currentNodeID = 0 is a tmp node without any sense
                if(currentNodeID > 0 && if_is_homopolymer_strict(r_string_site + currentNodeID - 1, r_string, r_string_length)
                  && max_count >= total_count*CORRECT_THRESHOLD_HOMOPOLYMER/** && max_type != MISMATCH**/)
                {
                    currentNodeID = add_path_to_correct_read(backbone, dumy, currentNodeID, max_type, max_edge, current_cigar,
                    self_string);
                }
                else///不矫正, 直接取下一个backbone节点
                {
                    currentNodeID++;
                    add_base_to_correct_read_directly(dumy, backbone->g_nodes.list[currentNodeID].base);

                    add_cigar_record(&(backbone->g_nodes.list[currentNodeID].base), 1, current_cigar, 0);
                }
            }
            
            ///fprintf(stderr, "currentNodeID: %d, max_type: %d\n", currentNodeID, max_type);
        }
        else  ///非backbone节点就会出错了
        {
            fprintf(stderr, "error\n");
        }
        
    }

}






/**
///从backbone_start遍历到backbone_end节点，生成出来的seq要接着放到dumy->corrected_read中
void get_seq_from_Graph(Graph* backbone, long long backbone_start, long long backbone_end, Correct_dumy* dumy)
{
    long long new_seq_length = 0;
    long long currentNodeID;
    long long i;
    // 总共有以下几种情况: 
    // 1. match 2. mismatch (A, C, G, T, N) 3. deletion 4. insertion (A, C, G, T)
    // 其实就是 1. 自己本身的weight 2. alignToNode的weight 3. insertion节点的weight
    long long max_count;
    int max_type;
    long long  max_node;
    long long total_count;
    long long nodeID;
    char current_base;

    currentNodeID = backbone_start;
    while (currentNodeID != backbone_end)
    {
        total_count = 0;
        max_count = -1;
        ///图上能够被遍历到的有两种节点
        ///1. backbone节点 2. insertion节点
        ///backbone节点才有match/mismatch/deletion
        ///insertion这些都没有，就是无脑看出边

        ///假如这是个backbone节点
        if (currentNodeID >= backbone_start && currentNodeID <= backbone_end)
        {
            ///match
            total_count += backbone->g_nodes.list[currentNodeID].weight;
            max_count = backbone->g_nodes.list[currentNodeID].weight;
            max_type = 0;
            max_node = currentNodeID;
            ///mismatch和deletion (A, C, G, T, N, D, 除了自己的那个字符)
            for (i = 0; i < backbone->g_nodes.list[currentNodeID].alignedTo_Nodes.length; i++)
            {
                nodeID = backbone->g_nodes.list[currentNodeID].alignedTo_Nodes.list[i].out_node;
                total_count += backbone->g_nodes.list[nodeID].weight;

                if (backbone->g_nodes.list[nodeID].weight > max_count)
                {
                    max_count = backbone->g_nodes.list[nodeID].weight;
                    max_type = 1;
                    max_node = nodeID;
                }
            }
            ///insertion (A, C, G, T, N)
            ///注意这个出边还得避开下一个backbone节点
            for (i = 0; i < backbone->g_nodes.list[currentNodeID].outcome_edges.length; i++)
            {
                if (backbone->g_nodes.list[currentNodeID].outcome_edges.list[i].weight == 2)
                {
                    
                    nodeID = backbone->g_nodes.list[currentNodeID].outcome_edges.list[i].out_node;
                    total_count += backbone->g_nodes.list[nodeID].weight;

                    if (backbone->g_nodes.list[nodeID].weight > max_count)
                    {
                        max_count = backbone->g_nodes.list[nodeID].weight;
                        max_type = 2;
                        max_node = nodeID;
                    }

                    ///拿到的nodeID应该一定不是backbone上的，如果是就错了
                    if (nodeID >= backbone_start && nodeID <= backbone_end)
                    {
                        fprintf(stderr, "error\n");
                    }
                    
                }
            }
        }
        else  ///如果是insertion节点，就无脑看出边
        {
            ///insertion (A, C, G, T, N)
            ///注意这个出边不用避开下一个backbone节点
            for (i = 0; i < backbone->g_nodes.list[currentNodeID].outcome_edges.length; i++)
            {
                if (backbone->g_nodes.list[currentNodeID].outcome_edges.list[i].weight == 2)
                {
                    
                    nodeID = backbone->g_nodes.list[currentNodeID].outcome_edges.list[i].out_node;
                    total_count += backbone->g_nodes.list[nodeID].weight;

                    if (backbone->g_nodes.list[nodeID].weight > max_count)
                    {
                        max_count = backbone->g_nodes.list[nodeID].weight;
                        max_type = 2;
                        max_node = nodeID;
                    }
                    
                }
                else  ///如果边不是2就不对了
                {
                    fprintf(stderr, "error\n");
                }
                
            }
        }
        
        





        if(max_count >= total_count*CORRECT_THRESHOLD)
        {
            current_base = backbone->g_nodes.list[max_node].base;
            ///说明是insertion
            if (max_type == 2)
            {
                currentNodeID = max_node;
            }
            else ///其他情况依然沿着backbone向前
            {
                currentNodeID++;
            }
        }
        else
        {
            ///假如这是个backbone节点, 不矫正
            if (currentNodeID >= backbone_start && currentNodeID <= backbone_end)
            {
                current_base = backbone->g_nodes.list[currentNodeID].base;
                currentNodeID++;
            }
            else///如果在insertion节点上不达标很麻烦...,只能选最大的了
            {
                current_base = backbone->g_nodes.list[max_node].base;
                ///说明是insertion
                if (max_type == 2)
                {
                    currentNodeID = max_node;
                }
                else ///insertion节点不可能出现这种情况
                {
                    fprintf(stderr, "error\n");
                }
            }
        }


        if (max_count <= 0)
        {
            fprintf(stderr, "error\n");
        }


        add_base_to_correct_read(dumy, current_base, max_type);
        
    }



    
    ///最后还要处理backbone_end这个节点
    total_count = 0;
    max_count = -1;
    ///这个节点肯定是backbone上的节点啊
    ///match
    total_count += backbone->g_nodes.list[currentNodeID].weight;
    max_count = backbone->g_nodes.list[currentNodeID].weight;
    max_type = 0;
    max_node = currentNodeID;
    ///mismatch和deletion (A, C, G, T, N, D, 除了自己的那个字符)
    for (i = 0; i < backbone->g_nodes.list[currentNodeID].alignedTo_Nodes.length; i++)
    {
        nodeID = backbone->g_nodes.list[currentNodeID].alignedTo_Nodes.list[i].out_node;
        total_count += backbone->g_nodes.list[nodeID].weight;

        if (backbone->g_nodes.list[nodeID].weight > max_count)
        {
            max_count = backbone->g_nodes.list[nodeID].weight;
            max_type = 1;
            max_node = nodeID;
        }
    }
    ///这个节点不应该有任何出边了
    if(backbone->g_nodes.list[currentNodeID].outcome_edges.length)
    {
        fprintf(stderr, "haha\n");
    }


    if(max_count >= total_count*CORRECT_THRESHOLD)
    {
        current_base = backbone->g_nodes.list[max_node].base;
    }
    else
    {
        current_base = backbone->g_nodes.list[currentNodeID].base;
    }


    if (max_count <= 0)
    {
        fprintf(stderr, "error\n");
    }
    
    add_base_to_correct_read(dumy, current_base, max_type);
    
}
**/










/**
///从backbone_start遍历到backbone_end节点，生成出来的seq要接着放到dumy->corrected_read中
void get_seq_from_Graph_Len2(Graph* backbone, long long backbone_start, long long backbone_end, Correct_dumy* dumy)
{
    long long new_seq_length = 0;
    long long currentNodeID;
    long long i;
    // 总共有以下几种情况: 
    // 1. match 2. mismatch (A, C, G, T, N) 3. deletion 4. insertion (A, C, G, T)
    // 其实就是 1. 自己本身的weight 2. alignToNode的weight 3. insertion节点的weight
    long long max_count;
    int max_type;
    long long  max_node;
    long long total_count;
    long long nodeID;
    char current_base;
    char buffer[2];

    currentNodeID = backbone_start;
    while (currentNodeID != backbone_end)
    {
        total_count = 0;
        max_count = -1;
        ///图上能够被遍历到的有两种节点
        ///1. backbone节点 2. insertion节点
        ///backbone节点才有match/mismatch/deletion
        ///insertion这些都没有，就是无脑看出边

        ///假如这是个backbone节点
        if (currentNodeID >= backbone_start && currentNodeID <= backbone_end)
        {
            ///match
            total_count += backbone->g_nodes.list[currentNodeID].weight;
            max_count = backbone->g_nodes.list[currentNodeID].weight;
            max_type = 0;
            max_node = currentNodeID;
            ///mismatch和deletion (A, C, G, T, N, D, 除了自己的那个字符)
            for (i = 0; i < backbone->g_nodes.list[currentNodeID].alignedTo_Nodes.length; i++)
            {
                nodeID = backbone->g_nodes.list[currentNodeID].alignedTo_Nodes.list[i].out_node;
                total_count += backbone->g_nodes.list[nodeID].weight;

                if (backbone->g_nodes.list[nodeID].weight > max_count)
                {
                    max_count = backbone->g_nodes.list[nodeID].weight;
                    max_type = 1;
                    max_node = nodeID;
                }
            }
            ///insertion (A, C, G, T, N)
            ///注意这个出边还得避开下一个backbone节点
            for (i = 0; i < backbone->g_nodes.list[currentNodeID].outcome_edges.length; i++)
            {
                if (backbone->g_nodes.list[currentNodeID].outcome_edges.list[i].weight == 2)
                {
                    ///fprintf(stderr, "error\n");
                    nodeID = backbone->g_nodes.list[currentNodeID].outcome_edges.list[i].out_node;
                    total_count += backbone->g_nodes.list[nodeID].weight;

                    if (backbone->g_nodes.list[nodeID].weight > max_count)
                    {
                        max_count = backbone->g_nodes.list[nodeID].weight;
                        max_type = 2;
                        max_node = nodeID;
                    }

                    ///拿到的nodeID应该一定不是backbone上的，如果是就错了
                    if (nodeID >= backbone_start && nodeID <= backbone_end)
                    {
                        fprintf(stderr, "error\n");
                    }
                    
                }
            }
        }
        else  ///如果是insertion节点，就无脑看出边
        {
            if (backbone->g_nodes.list[currentNodeID].outcome_edges.length!=1)
            {
                fprintf(stderr, "error000\n");
            }
            
            ///insertion (A, C, G, T, N)
            ///注意这个出边不用避开下一个backbone节点
            for (i = 0; i < backbone->g_nodes.list[currentNodeID].outcome_edges.length; i++)
            {
                if (backbone->g_nodes.list[currentNodeID].outcome_edges.list[i].weight == 2)
                {
                    
                    nodeID = backbone->g_nodes.list[currentNodeID].outcome_edges.list[i].out_node;
                    total_count += backbone->g_nodes.list[nodeID].weight;

                    if (backbone->g_nodes.list[nodeID].weight > max_count)
                    {
                        max_count = backbone->g_nodes.list[nodeID].weight;
                        max_type = 2;
                        max_node = nodeID;
                    }
                    
                }
                else  ///如果边不是2就不对了
                {
                    fprintf(stderr, "error\n");
                }
                
            }
        }
        
        





        if(max_count >= total_count*CORRECT_THRESHOLD)
        {
            current_base = backbone->g_nodes.list[max_node].base;
            ///说明是insertion
            if (max_type == 2)
            {
                currentNodeID = max_node;
            }
            else ///其他情况依然沿着backbone向前
            {
                currentNodeID++;
            }
        }
        else
        {
            ///假如这是个backbone节点, 不矫正
            if (currentNodeID >= backbone_start && currentNodeID <= backbone_end)
            {
                current_base = backbone->g_nodes.list[currentNodeID].base;
                currentNodeID++;
            }  ///应该不存在这个问题
            else///如果在insertion节点上不达标很麻烦...,只能选最大的了
            {
                ///因为现在每个insert节点只有一个出边，且这个出边到backbone
                fprintf(stderr, "error111\n");
            }
        }


        if (max_count <= 0)
        {
            fprintf(stderr, "error\n");
        }

        if (current_base < 'A')
        {
            buffer[0] = s_H[(current_base >> 2) & ((uint8_t)3)];
            buffer[1] = s_H[current_base & ((uint8_t)3)];
            add_base_to_correct_read(dumy, buffer[0], max_type);
            add_base_to_correct_read(dumy, buffer[1], max_type);
        }
        else
        {
            add_base_to_correct_read(dumy, current_base, max_type);
        }
        
    }



    
    ///最后还要处理backbone_end这个节点
    total_count = 0;
    max_count = -1;
    ///这个节点肯定是backbone上的节点啊
    ///match
    total_count += backbone->g_nodes.list[currentNodeID].weight;
    max_count = backbone->g_nodes.list[currentNodeID].weight;
    max_type = 0;
    max_node = currentNodeID;
    ///mismatch和deletion (A, C, G, T, N, D, 除了自己的那个字符)
    for (i = 0; i < backbone->g_nodes.list[currentNodeID].alignedTo_Nodes.length; i++)
    {
        nodeID = backbone->g_nodes.list[currentNodeID].alignedTo_Nodes.list[i].out_node;
        total_count += backbone->g_nodes.list[nodeID].weight;

        if (backbone->g_nodes.list[nodeID].weight > max_count)
        {
            max_count = backbone->g_nodes.list[nodeID].weight;
            max_type = 1;
            max_node = nodeID;
        }
    }
    ///这个节点不应该有任何出边了
    if(backbone->g_nodes.list[currentNodeID].outcome_edges.length)
    {
        fprintf(stderr, "haha\n");
    }


    if(max_count >= total_count*CORRECT_THRESHOLD)
    {
        current_base = backbone->g_nodes.list[max_node].base;
    }
    else
    {
        current_base = backbone->g_nodes.list[currentNodeID].base;
    }


    if (max_count <= 0)
    {
        fprintf(stderr, "error\n");
    }
    
    add_base_to_correct_read(dumy, current_base, max_type);
    
}
**/


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
    long long i;
    long long y_start, y_length;
    long long overlapID, windowID;
    long long startNodeID, endNodeID, currentNodeID;

    ///这个和前面算alignment还不一样
    ///那个时候x_start和x_end是当前窗口内的overlap的起始和结束位置
    ///这个window就是要做consensus啊，所以起始和结束就是window本身，固定的
    backbone = r_string + window_start;
    backbone_length = window_end - window_start + 1;

    addUnmatchedSeqToGraph(g, backbone, backbone_length, &startNodeID, &endNodeID);

    
    long long correct_x_pos_s;
    ///与当前window重叠的所有overlap
    for (i = 0; i < dumy->length; i++)
    {
        ///这个是那个overlap的ID，而不是overlap里对应窗口的ID
        overlapID = dumy->overlapID[i];

        correct_x_pos_s = (overlap_list->list[overlapID].x_pos_s / WINDOW) * WINDOW;
        windowID = (window_start - correct_x_pos_s) / WINDOW;

        ///如果这个window不匹配，跳过
        if (overlap_list->list[overlapID].w_list[windowID].y_end == -1)
        {
            continue;
        }
        

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
        ///这个是比对上的起始base在backbone上对应的位置，也就是节点ID
        currentNodeID = x_start - window_start;
        ///这个是要用的cigar: overlap_list->list[overlapID].w_list[windowID].cigar;


        /**
        if(memcmp("m54238_180909_174539/6947324/ccs", Get_NAME((*R_INF),overlap_list->list[overlapID].x_id), 
     Get_NAME_LENGTH((*R_INF), overlap_list->list[overlapID].x_id)) == 0 && window_start == 12750 && window_end == 13124)
        {
            fprintf(stderr, "********x_start: %d, window_start: %d, window_end: %d, dumy->length: %d, y_name: %.*s\n", 
            x_start, window_start, window_end, dumy->length, 
            Get_NAME_LENGTH((*R_INF), overlap_list->list[overlapID].y_id), Get_NAME((*R_INF),overlap_list->list[overlapID].y_id));


            for (int ijk = 0; ijk < overlap_list->list[overlapID].w_list[windowID].cigar.length; ijk++)
            {
                fprintf(stderr, "###### Oper: %d, Len: %d\n", 
                overlap_list->list[overlapID].w_list[windowID].cigar.C_C[ijk],
                overlap_list->list[overlapID].w_list[windowID].cigar.C_L[ijk]);
            }
        }
        

       if(memcmp("m64011_190326_191011/163906371/ccs", Get_NAME((*R_INF),overlap_list->list[0].x_id), 
     Get_NAME_LENGTH((*R_INF), overlap_list->list[0].x_id)) == 0 && window_start == 9375 && window_end == 9749)
        {
            fprintf(stderr, "y_name: %.*s\n", Get_NAME_LENGTH((*R_INF), overlap_list->list[overlapID].y_id), 
                Get_NAME((*R_INF),overlap_list->list[overlapID].y_id));

            reverse_complement(x_string, x_length);
            reverse_complement(y_string, y_length);

            fprintf(stderr, "x: %.*s\n", x_length, x_string);
            fprintf(stderr, "y: %.*s\n", y_length, y_string);
            for (int ijk = 0; ijk < overlap_list->list[overlapID].w_list[windowID].cigar.length; ijk++)
            {
                fprintf(stderr, "###### Oper: %d, Len: %d\n", 
                overlap_list->list[overlapID].w_list[windowID].cigar.C_C[ijk],
                overlap_list->list[overlapID].w_list[windowID].cigar.C_L[ijk]);
            }
            fprintf(stderr,"\n");

            reverse_complement(x_string, x_length);
            reverse_complement(y_string, y_length);
            
            
        }
        **/

        addmatchedSeqToGraph(g, currentNodeID, x_string, x_length, 
                    y_string, y_length, &(overlap_list->list[overlapID].w_list[windowID].cigar), startNodeID, endNodeID);

        
    }


    /**
    if(memcmp("m54238_180922_175520/52363405/ccs", Get_NAME((*R_INF),overlap_list->list[overlapID].x_id), 
     Get_NAME_LENGTH((*R_INF), overlap_list->list[overlapID].x_id)) == 0 && window_start == 6750 && window_end == 7124)
     {
         fprintf(stderr, "dumy->corrected_read_length: %d\n", dumy->corrected_read_length);
     }

    if(memcmp("m54238_180922_175520/52363405/ccs", Get_NAME((*R_INF),overlap_list->list[overlapID].x_id), 
     Get_NAME_LENGTH((*R_INF), overlap_list->list[overlapID].x_id)) == 0 && window_start == 6750 && window_end == 7124)
     {
         get_seq_from_Graph_print(g, dumy, current_cigar, backbone, r_string, r_total_length, window_start);
     }
     else
     {
         get_seq_from_Graph(g, dumy, current_cigar, backbone, r_string, r_total_length, window_start);
     }
     

    

    


    if(memcmp("m54238_180922_175520/52363405/ccs", Get_NAME((*R_INF),overlap_list->list[overlapID].x_id), 
     Get_NAME_LENGTH((*R_INF), overlap_list->list[overlapID].x_id)) == 0 && window_start == 6750 && window_end == 7124)
    {
        fprintf(stderr, "dumy->corrected_read_length: %d\n", dumy->corrected_read_length);
    }
    **/

    get_seq_from_Graph(g, DAGCon, dumy, current_cigar, backbone, r_string, r_total_length, window_start);

    ///get_seq_from_Graph(g, startNodeID, endNodeID, dumy);
    ///get_seq_from_Graph_Len2(g, startNodeID, endNodeID, dumy);


    
    ///debug_graph(g, backbone_length);
    
    

    /**
    for (i = 0; i < dumy->length; i++)
    {
        ///这个是那个overlap的ID，而不是overlap里对应窗口的ID
        overlapID = dumy->overlapID[i];

        correct_x_pos_s = (overlap_list->list[overlapID].x_pos_s / WINDOW) * WINDOW;
        windowID = (window_start - correct_x_pos_s) / WINDOW;

        ///如果这个window不匹配，跳过
        if (overlap_list->list[overlapID].w_list[windowID].y_end == -1)
        {
            continue;
        }
        

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
        ///这个是比对上的起始base在backbone上对应的位置，也就是节点ID
        currentNodeID = x_start - window_start;
        ///这个是要用的cigar: overlap_list->list[overlapID].w_list[windowID].cigar;

        Graph_debug(g, currentNodeID, x_string, x_length, 
                    y_string, y_length, &(overlap_list->list[overlapID].w_list[windowID].cigar), startNodeID, endNodeID);
    }

    for (i = 0; i < g->g_nodes.length; i++)
    {
        if (i >= startNodeID && i <= endNodeID)
        {
            if (g->g_nodes.list[i].weight != 1)
            {
                fprintf(stderr, "error 1\n");
            }
        }
        else
        {
            if (g->g_nodes.list[i].weight != 0)
            {
                fprintf(stderr, "error 2\n");
            }

            if (g->g_nodes.list[i].alignedTo_Nodes.length != 0)
            {
                fprintf(stderr, "error 3\n");
            }

            ///节点入边不为0，说明这个不是alignTO节点，而是insert节点
            if (g->g_nodes.list[i].income_edges.length != 0 && g->g_nodes.list[i].outcome_edges.length != 0)
            {
                fprintf(stderr, "error 4\n");
            }
        }
    }
    **/
    
    


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

    for (i = 0; i < new_cigar->length; i++)
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
    Correct_dumy* new_dumy = &(second_round->dumy);


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
    int cigar_error = 0;

    long long get_x_start, get_x_end, get_y_start, get_y_end;
    get_x_start = get_x_end = get_y_start = get_y_end = -1;

    int start_cigar = -1;
    int end_cigar = -1;
    char merge_base;

    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2是x缺字符（y多字符），而3是y缺字符（x多字符）
    ///obtained x_i may larger than start_base/end_base
    ///when operation == 3
    ///so for operation == 3, we need deal with carefully
    for (i = 0; i < new_cigar->length; i++)
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
        }///2是x缺字符（y多字符）
        else if (operation == 2)
        {
            y_i += operationLen;
        }///3是y缺字符（x多字符）
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
    
    /**
    fprintf(stderr, "get_x_start: %d, get_x_end: %d\n", get_x_start, get_x_end);
    fprintf(stderr, "get_y_start: %d, get_y_end: %d\n", get_y_start, get_y_end);
    for (int ijk = 0; ijk < new_cigar->length; ijk++)
    {
        fprintf(stderr, "Oper: %d, Len: %d\n", Get_Cigar_Type(new_cigar->record[ijk]), 
        Get_Cigar_Length(new_cigar->record[ijk]));
    }
    **/
    

    x_i = 0;
    y_i = 0;

    uint32_t single_record = 0;

    for (i = 0; i < new_cigar->length; i++)
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
        }///3是y缺字符（x多字符）
        else if (operation == 3)
        {
            x_i += operationLen;
        }
    }

    new_cigar->length = end_cigar - start_cigar + 1;
    ///可以优化
    memmove(new_cigar->record, new_cigar->record + start_cigar, new_cigar->length*sizeof(uint32_t));
    
    long long total_x_start = total_window_start + get_x_start;
    long long x_length = get_x_end -get_x_start + 1;
    long long total_y_start = get_y_start;
    long long y_length = get_y_end -get_y_start + 1;


    add_cigar_to_cigar(current_dumy, current_cigar, second_round,
    total_x_start, x_length, total_y_start, y_length);

}

int process_boundary(overlap_region_alloc* overlap_list, All_reads* R_INF, Correct_dumy* dumy, Graph* g, Graph* DAGCon,
Cigar_record* current_cigar, long long uncorrected_window_start, Round2_alignment* second_round)
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
    char* x_string;
    char* y_string;
    char* backbone;
    long long backbone_length;
    long long i;
    long long y_start, y_length;
    long long overlapID, windowID;
    long long startNodeID, endNodeID, currentNodeID;
    int end_site;
    unsigned int error;
    int real_y_start;
    window_list tmp_cigar;
    long long total_error = 0;

    backbone = r_string + corrected_window_start;
    backbone_length = corrected_window_end - corrected_window_start + 1;
    addUnmatchedSeqToGraph(g, backbone, backbone_length, &startNodeID, &endNodeID);


    long long correct_x_pos_s;
    long long matched_coverage = 0;
    for (i = 0; i < dumy->length; i++)
    {
        overlapID = dumy->overlapID[i];
        correct_x_pos_s = (overlap_list->list[overlapID].x_pos_s / WINDOW) * WINDOW;
        windowID = (uncorrected_window_start - correct_x_pos_s) / WINDOW;

        ///skip if window is unmatched
        if (overlap_list->list[overlapID].w_list[windowID].y_end == -1)
        {
            continue;
        }

        x_start = overlap_list->list[overlapID].w_list[windowID].x_start;
        y_start = overlap_list->list[overlapID].w_list[windowID].y_start;

        


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
            threshold = x_len * THRESHOLD_RATE;
            /****************************may have bugs********************************/
            threshold = Adjust_Threshold(threshold, x_len);
            /****************************may have bugs********************************/
            ///y_start may less than 0
            y_start = y_start - WINDOW_BOUNDARY/2;

            ///其实可以不加...怕出bug
            if(y_start < 0)
            {
                continue;
            }

            Window_Len = x_len + (threshold << 1);

            error =(unsigned int)-1;
            if(determine_overlap_region(threshold, y_start, overlap_list->list[overlapID].y_id, Window_Len, R_INF,
            &extra_begin, &extra_end, &y_start, &o_len))
            {
                fill_subregion(dumy->overlap_region, y_start, o_len, overlap_list->list[overlapID].y_pos_strand, 
                R_INF, overlap_list->list[overlapID].y_id, extra_begin, extra_end);

                x_string = r_string + x_start;
                y_string = dumy->overlap_region;

                ///both end site and real_y_start have extra_begin
                ///有很多是完全匹配，可以先快速判断是不是完全匹配
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
                y_start = overlap_list->list[overlapID].w_list[windowID].y_start - WINDOW_BOUNDARY/2;

                ///其实可以不加...怕出bug
                if(y_start < 0)
                {
                    continue;
                }

                error =(unsigned int)-1;
                if(determine_overlap_region(threshold, y_start, overlap_list->list[overlapID].y_id, Window_Len, R_INF,
                &extra_begin, &extra_end, &y_start, &o_len))
                {
                    fill_subregion(dumy->overlap_region, y_start, o_len, overlap_list->list[overlapID].y_pos_strand, 
                    R_INF, overlap_list->list[overlapID].y_id, extra_begin, extra_end);

                    x_string = r_string + x_start;
                    y_string = dumy->overlap_region;

                    ///both end site and real_y_start have extra_begin
                    ///有很多是完全匹配，可以先快速判断是不是完全匹配
                    end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
                            &(dumy->path_length), dumy->matrix_bit, dumy->path, -1, -1);
                }
            }

            if (error!=(unsigned int)-1)
            {

                

                total_error = total_error + error;
                matched_coverage++;
                tmp_cigar.x_start = x_start;
                tmp_cigar.x_end = x_end;
                generate_cigar(dumy->path, dumy->path_length, &tmp_cigar, &real_y_start, &end_site, &error, 
                x_string, x_len, y_string);  
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

    //             if(memcmp("m54238_180922_175520/52363405/ccs", Get_NAME((*R_INF),overlap_list->list[overlapID].x_id), 
    //  Get_NAME_LENGTH((*R_INF), overlap_list->list[overlapID].x_id)) == 0 && uncorrected_window_start == 6750)
    //     {
    //         fprintf(stderr, "********x_start: %d, uncorrected_window_start: %d, dumy->length: %d, y_name: %.*s\n", 
    //         x_start, uncorrected_window_start, dumy->length, 
    //         Get_NAME_LENGTH((*R_INF), overlap_list->list[overlapID].y_id), 
    //         Get_NAME((*R_INF),overlap_list->list[overlapID].y_id));


    //          fprintf(stderr, "*******error: %d****\n", error);
    //             for (int ijk = 0; ijk < tmp_cigar.cigar.length; ijk++)
    //             {
    //                 fprintf(stderr, "Oper: %d, Len: %d\n", tmp_cigar.cigar.C_C[ijk],
    //                 tmp_cigar.cigar.C_L[ijk]);
    //             }
                
    //     }

                
                // if(verify_cigar(x_string, x_length, y_string, y_length, &tmp_cigar.cigar, 
                //     error))
                // {
                //     fprintf(stderr, "*******error: %d****\n", error);
                //     for (int ijk = 0; ijk < tmp_cigar.cigar.length; ijk++)
                //     {
                //         fprintf(stderr, "Oper: %d, Len: %d\n", tmp_cigar.cigar.C_C[ijk],
                //         tmp_cigar.cigar.C_L[ijk]);
                //     }
                    
                //     fprintf(stderr, "*******dumy->path_length: %d\n****\n", dumy->path_length);
                // }
                

                currentNodeID = x_start - corrected_window_start;
                
                addmatchedSeqToGraph(g, currentNodeID, x_string, x_length, 
                    y_string, y_length, &(tmp_cigar.cigar), startNodeID, endNodeID);
            }

        }///case 2 is useless
        else if(x_start != uncorrected_window_start)
        {
            continue;
        }        
    }

    /**
    fprintf(stderr, "matched_coverage: %d, dumy->length: %d\n", 
    matched_coverage, dumy->length);
    **/

    if(matched_coverage >= MIN_COVERAGE_THRESHOLD)
    {
        ///if there are no error, we do not need correction
        if(total_error == 0)
        {
            return 0;
        }
        /**
        fprintf(stderr, "s_start_nodeID: %d, s_end_nodeID: %d, corrected_window_start: %d, corrected_window_end: %d\n", 
    g->s_start_nodeID, g->s_end_nodeID, corrected_window_start, corrected_window_end);
    **/
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
                        Round2_alignment* second_round)
{
    clear_Cigar_record(current_cigar);
    
    long long overlap_length;
    long long window_start, window_end;

    long long num_availiable_win = 0;
    
    


    Window_Pool w_inf;
    init_Window_Pool(&w_inf, g_read->length, WINDOW, TAIL_LENGTH);

    int flag = 0;
    ///for last window
    dumy->last_boundary_length = 0;
    while(get_Window(&w_inf, &window_start, &window_end) && flag != -2)
    {

        

        dumy->length = 0;
        dumy->lengthNT = 0;
        ///flag返回的是重叠数量
        ///dumy->length返回的是有效完全重叠的数量
        ///dumy->lengthNT返回的是有效不完全重叠的数量
        ///return overlaps that is overlaped with [window_start, window_end]
        flag = get_available_interval(window_start, window_end, overlap_list, dumy);
        switch (flag)
        {
            case 1:    ///找到匹配
                break;
            case 0:    ///没找到匹配
                break;
            case -2: ///下一个window也不会存在匹配, 直接跳出
                break;
        }


        ///这个是available overlap里所有window的数量...
        ///num_availiable_win = num_availiable_win + dumy->length + dumy->lengthNT;
        num_availiable_win = num_availiable_win + dumy->length;
        
        ///重叠窗口数，也就是coverage大小
        if(dumy->length >= MIN_COVERAGE_THRESHOLD)
        {

            window_consensus(g_read->seq, g_read->length, window_start, window_end, overlap_list, 
            dumy, R_INF, g, DAGCon, current_cigar);

            if(dumy->last_boundary_length != 0)
            {
                process_boundary(overlap_list, R_INF, dumy, g, DAGCon, current_cigar, 
                window_start, second_round);
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
    long long i;
    int flag = 0;
    long long Begin, End, Len;
    long long overlap_length;


    for (i = dumy->start_i; i < overlap_list->length; i++)
    {
        ///只会发生在这个interval比list里所有元素都小的情况
        ///这种情况下一个interval需要从0开始
        if (window_end < overlap_list->list[i].x_pos_s)
        {
            dumy->start_i = 0;
            return 0;
        }
        else ///只要window_end >= overlap_list->list[i].x_pos_s，就有可能重叠
        {
            dumy->start_i = i;
            break;
        }
    }

    ///只会发生在这个window比list里所有元素都大的情况
    ///这种情况下一个window也无需遍历了
    if (i >= overlap_list->length)
    {
        dumy->start_i = overlap_list->length;
        return -2;
    }
    



    
    long long fake_length = 0;
    overlap_length = window_end - window_start + 1;
    (*real_length) = 0;

    for (; i < overlap_list->length; i++)
    {
        ///是否重叠   
        if((Len = OVERLAP(window_start, window_end, overlap_list->list[i].x_pos_s, overlap_list->list[i].x_pos_e)) > 0)
        {
            ///重叠数量
            fake_length++;

            if (overlap_length == Len && overlap_list->list[i].is_match == 1)
            {
                (*real_length)++;
            }

            if (overlap_length == Len && overlap_list->list[i].is_match == 100)
            {
                (*real_length_100)++;
            }
        }
        
        if(overlap_list->list[i].x_pos_s > window_end)
        {
            break;
        }
    }

    ///fake_length是重叠的数量，而不是有效重叠的数量
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
    
    long long overlap_length;
    long long window_start, window_end;
    int return_flag = 1;
    (*abnormal) = 0;
    
    Window_Pool w_inf;
    init_Window_Pool(&w_inf, g_read->length, WINDOW, TAIL_LENGTH);

    int flag = 0;
    long long realLen, tmpLen;
    while(get_Window(&w_inf, &window_start, &window_end) && flag != -2)
    {
        dumy->length = 0;
        dumy->lengthNT = 0;
        ///flag返回的是重叠数量
        ///dumy->length返回的是有效完全重叠的数量
        ///dumy->lengthNT返回的是有效不完全重叠的数量
        ///return overlaps that is overlaped with [window_start, window_end]
        flag = get_available_fully_covered_interval(window_start, window_end, 
        overlap_list, dumy, &realLen, &tmpLen);


        switch (flag)
        {
            case 1:    ///找到匹配
                break;
            case 0:    ///没找到匹配
                break;
            case -2: ///下一个window也不会存在匹配, 直接跳出
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
    //return 1;
}







void markSNP(
long long window_offset,
long long x_total_start, long long x_length, 
long long y_total_start, long long y_length, 
CIGAR* cigar, haplotype_evdience_alloc* hap)
{
    
    int x_i, y_i, cigar_i;
    x_i = 0;
    y_i = 0;
    cigar_i = 0;
    int operation;
    int operationLen;
    int i;
    ///mismatches based on the offset of x
    long long inner_offset = x_total_start - window_offset;

    
    ///note that node 0 is the start node
    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2 represents thre are more bases at y
    ///3 represents thre are more bases at x
    while (cigar_i < cigar->length)
    {
        operation = cigar->C_C[cigar_i];
        operationLen = cigar->C_L[cigar_i];

        ///match
        if (operation == 0)
        {
            x_i += operationLen;
            y_i += operationLen;
        }
        else if(operation == 1) ///mismatch
        {
            for (i = 0; i < operationLen; i++)
            {
                if(hap->flag[inner_offset + x_i] < 127)
                {
                    hap->flag[inner_offset + x_i]++;
                }
                
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
            x_i += operationLen;
        }
        
        cigar_i++;
    }
}


///mark SNPs at [xBeg, xEnd], note we need to deal with flag_offset carefully
void markSNP_detail(CIGAR* cigar_record, uint8_t* flag, 
long long xBeg, long long xEnd, long long flag_offset)
{
    if(xBeg > xEnd) return;

    int operation;
    int operationLen;
    long long x_i, y_i, cigar_i, i;
    i = cigar_i = x_i = y_i = 0;

    while (cigar_i < cigar_record->length)
    {
        operation = cigar_record->C_C[cigar_i];
        operationLen = cigar_record->C_L[cigar_i];
        if(x_i > xEnd)
        {
            break;
        }
        
        ///match
        if (operation == 0)
        {
            x_i += operationLen;
            y_i += operationLen;
        }
        else if(operation == 1) ///mismatch
        {
            for (i = 0; i < operationLen; i++)
            {
                /// note we need to deal with flag_offset carefully
                if(flag[x_i - flag_offset] < 127 && x_i >= xBeg && x_i <= xEnd)
                {
                    flag[x_i - flag_offset]++;
                }
                
                x_i++;
                y_i++;
            }
        }///insertion, that means y has more bases than x
        else if (operation == 2)
        {
            y_i += operationLen;
        }
        else if (operation == 3)
        {
            x_i += operationLen;
        }
        
        cigar_i++;
    }
}


///window_offset is still the x-based offset
///x_total_start and y_total_start are global positions, instead of local positions
void markSNP_advance(
long long window_offset,
long long x_total_start, long long x_length, 
long long y_total_start, long long y_length, 
window_list* current_cigar, window_list* beg_cigar, window_list* end_cigar, 
haplotype_evdience_alloc* hap)
{   
    long long x_total_end = x_total_start + x_length - 1;
    ///mismatches based on the offset of x
    long long inner_offset = x_total_start - window_offset;
    CIGAR* cigar_record;
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
        cigar_record = &(beg_cigar->cigar);
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

            markSNP_detail(cigar_record, hap->flag + inner_offset, x_interval_beg, 
            x_interval_end, x_interval_beg);
        }
    }

    if(end_cigar!=NULL && end_cigar->y_end!=-1)
    {
        ///useless_side = end_cigar->error_threshold;
        L_useless_side = end_cigar->extra_begin;
        R_useless_side = end_cigar->extra_end;
        cigar_record = &(end_cigar->cigar);
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
            
            markSNP_detail(cigar_record, hap->flag + end_cigar->x_start - window_offset,
            x_interval_beg, x_interval_end, 0);
        }
    }

    markSNP_detail(&(current_cigar->cigar), hap->flag + inner_offset, current_cigar_beg, 
    current_cigar_end, 0);

    // end_mark:
    // if(x_length >= 4)
    // {
    //     markSNP_detail(&(current_cigar->cigar), hap->flag + inner_offset, 0, x_length/2, 0);
    //     markSNP_detail(&(current_cigar->cigar), hap->flag + inner_offset, x_length/2+1, x_length - 1, 0);
    // }
    // else
    // {
    //     markSNP_detail(&(current_cigar->cigar), hap->flag + inner_offset, 0, x_length - 1, 0);
    // }
}




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
                    addHaplotypeEvdience(hap, &ev);
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
                    addHaplotypeEvdience(hap, &ev);
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
            /****************************may have bugs********************************/
            for (i = 0; i < operationLen; i++)
            {
                if(hap->flag[inner_offset] > snp_threshold)
                {
                    ev.misBase =  'N';
                    ev.overlapID = overlapID;
                    ev.site = x_total_start + x_i;
                    ev.overlapSite = y_total_start + y_i;
                    ev.type = 2;
                    addHaplotypeEvdience(hap, &ev);
                }

                inner_offset++;
                x_i++;
            }
            /****************************may have bugs********************************/
        }
        
        cigar_i++;
    }
}







///mark SNPs at [xBeg, xEnd], note we need to deal with flag_offset carefully
void addSNPtohaplotype_details(CIGAR* cigar_record, uint8_t* flag, 
char* x_string, char* y_string, long long x_total_start, long long y_total_start,
long long xBeg, long long xEnd, int overlapID, long long flag_offset, 
haplotype_evdience_alloc* hap, long long snp_threshold)
{
    if(xBeg > xEnd) return;
    int operation;
    int operationLen;
    long long x_i, y_i, cigar_i, i;
    i = cigar_i = x_i = y_i = 0;

    haplotype_evdience ev;

    ///note that node 0 is the start node
    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2 represents thre are more bases at y
    ///3 represents thre are more bases at x
    while (cigar_i < cigar_record->length)
    {
        operation = cigar_record->C_C[cigar_i];
        operationLen = cigar_record->C_L[cigar_i];
        if(x_i > xEnd)
        {
            break;
        }

        ///matches
        if (operation == 0)
        {
            for (i = 0; i < operationLen; i++)
            {
                ///should be at least 2 mismatches
                /// note we need to deal with flag_offset carefully
                if(flag[x_i - flag_offset] > snp_threshold && x_i >= xBeg && x_i <= xEnd)
                {
                    ev.misBase =  y_string[y_i];
                    ev.overlapID = overlapID;
                    ev.site = x_total_start + x_i;
                    ev.overlapSite = y_total_start + y_i;
                    ev.type = 0;
                    addHaplotypeEvdience(hap, &ev);
                }
                ///inner_offset++;
                x_i++;
                y_i++;
            }
        }
        else if(operation == 1)
        {
            for (i = 0; i < operationLen; i++)
            {

                /// should be at least 2 mismatches
                /// note we need to deal with flag_offset carefully
                if(flag[x_i - flag_offset] > snp_threshold && x_i >= xBeg && x_i <= xEnd)
                {
                    ev.misBase =  y_string[y_i];
                    ev.overlapID = overlapID;
                    ev.site = x_total_start + x_i;
                    ev.overlapSite = y_total_start + y_i;
                    ev.type = 1;
                    addHaplotypeEvdience(hap, &ev);
                }

                ///inner_offset++;
                x_i++;
                y_i++;
            }
        }///insertion, 2 represents thre are more bases at y
        else if (operation == 2)
        {
            y_i += operationLen;
        }///3 represents thre are more bases at x
        else if (operation == 3)
        {
            /****************************may have bugs********************************/
            for (i = 0; i < operationLen; i++)
            {
                ///if(hap->flag[inner_offset] > snp_threshold)
                /// should be at least 2 mismatches
                /// note we need to deal with flag_offset carefully
                if(flag[x_i - flag_offset] > snp_threshold && x_i >= xBeg && x_i <= xEnd)
                {
                    ev.misBase =  'N';
                    ev.overlapID = overlapID;
                    ev.site = x_total_start + x_i;
                    ev.overlapSite = y_total_start + y_i;
                    ev.type = 2;
                    addHaplotypeEvdience(hap, &ev);
                }

                ///inner_offset++;
                x_i++;
            }
            /****************************may have bugs********************************/
        }
        cigar_i++;
    }
}


void addSNPtohaplotype_advance(
long long window_offset, int overlapID,
long long x_total_start, long long x_length, 
long long y_total_start, long long y_length, 
window_list* current_cigar, window_list* beg_cigar, window_list* end_cigar, 
haplotype_evdience_alloc* hap, int snp_threshold, char* x_T_string, char* y_T_string)
{
    long long x_total_end = x_total_start + x_length - 1;
    long long inner_offset = x_total_start - window_offset;
    CIGAR* cigar_record;
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
        cigar_record = &(beg_cigar->cigar);
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
            ///x_interval_end = x_interval_beg + (xrightLen - useless_side) - 1;
            x_interval_end = x_interval_beg + (xrightLen - R_useless_side) - 1;
            ///current_cigar_beg is the offset of the current cigar
            ///that is the beg of current_cigar_beg
            ///current_cigar_beg = xrightLen - useless_side;
            current_cigar_beg = xrightLen - R_useless_side;

            // markSNP_detail(cigar_record, hap->flag + inner_offset, x_interval_beg, 
            // x_interval_end, x_interval_beg);
            addSNPtohaplotype_details(cigar_record, hap->flag + inner_offset, 
            x_T_string + beg_cigar->x_start, y_T_string + beg_cigar->y_start, 
            beg_cigar->x_start, beg_cigar->y_start, x_interval_beg, x_interval_end, 
            overlapID, x_interval_beg, hap, snp_threshold);
        }
    }


    if(end_cigar!=NULL && end_cigar->y_end!=-1)
    {
        ///useless_side = end_cigar->error_threshold;
        L_useless_side = end_cigar->extra_begin;
        R_useless_side = end_cigar->extra_end;
        cigar_record = &(end_cigar->cigar);
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
            addSNPtohaplotype_details(cigar_record, hap->flag + end_cigar->x_start - window_offset,
            x_T_string + end_cigar->x_start, y_T_string + end_cigar->y_start, 
            end_cigar->x_start, end_cigar->y_start, x_interval_beg, x_interval_end, 
            overlapID, 0, hap, snp_threshold);
        }
    }

    // markSNP_detail(&(current_cigar->cigar), hap->flag + inner_offset, current_cigar_beg, 
    // current_cigar_end, 0);
    addSNPtohaplotype_details(&(current_cigar->cigar), hap->flag + inner_offset, 
    x_T_string + current_cigar->x_start, y_T_string + current_cigar->y_start, 
    current_cigar->x_start, current_cigar->y_start, current_cigar_beg, 
    current_cigar_end, overlapID, 0, hap, snp_threshold);

    /**
    end_add:
    if(x_length >= 4)
    {
        addSNPtohaplotype_details(&(current_cigar->cigar), hap->flag + inner_offset, 
        x_T_string + current_cigar->x_start, y_T_string + current_cigar->y_start, 
        current_cigar->x_start, current_cigar->y_start, 0, 
        x_length/2, overlapID, 0, hap, snp_threshold);

        addSNPtohaplotype_details(&(current_cigar->cigar), hap->flag + inner_offset, 
        x_T_string + current_cigar->x_start, y_T_string + current_cigar->y_start, 
        current_cigar->x_start, current_cigar->y_start, x_length/2 + 1, 
        x_length - 1, overlapID, 0, hap, snp_threshold);
    }
    else
    {
        addSNPtohaplotype_details(&(current_cigar->cigar), hap->flag + inner_offset, 
        x_T_string + current_cigar->x_start, y_T_string + current_cigar->y_start, 
        current_cigar->x_start, current_cigar->y_start, 0, 
        x_length - 1, overlapID, 0, hap, snp_threshold);
    }
    **/
}













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
    long long startNodeID, endNodeID, currentNodeID;

    long long correct_x_pos_s;
    long long inner_window_offset;
    int snp_threshold;
    snp_threshold = 1;

    ///all overlaps related to the current window [window_start, window_end]
    ///first mark all snp pos
    for (i = 0; i < dumy->length; i++)
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


    /****************************may have bugs********************************/
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
            hap->snp++;
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
    for (i = 0; i < dumy->length; i++)
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


void get_related_cigars(window_list_alloc* boundary_cigars, long long id, window_list** beg_cigar, 
window_list** end_cigar)
{
    (*beg_cigar) = &(boundary_cigars->buffer[id*2]);
    (*end_cigar) = &(boundary_cigars->buffer[id*2+1]);
}

void cluster_advance(char* r_string, long long window_start, long long window_end, 
overlap_region_alloc* overlap_list, Correct_dumy* dumy, All_reads* R_INF, 
haplotype_evdience_alloc* hap, UC_Read* overlap_read)
{
    window_list* beg_cigar;
    window_list* end_cigar;
    ///window_start, window_end, and useful_length correspond to x, instead of y
    long long useful_length = window_end - window_start + 1;
    long long x_start, x_length;
    char* x_string;
    char* y_string;
    long long i;
    long long y_start, y_length;
    long long overlapID, windowID;
    long long startNodeID, endNodeID, currentNodeID;

    long long correct_x_pos_s;
    long long inner_window_offset;
    int snp_threshold;
    snp_threshold = 1;


    ///all overlaps related to the current window [window_start, window_end]
    ///first mark all snp pos
    for (i = 0; i < dumy->length; i++)
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

        beg_cigar = end_cigar = NULL;
        
        if(windowID >= 1)
        {
            beg_cigar = &(overlap_list->list[overlapID].boundary_cigars.buffer[windowID-1]);
        }

        if(windowID < overlap_list->list[overlapID].w_list_length - 1)
        {
            end_cigar = &(overlap_list->list[overlapID].boundary_cigars.buffer[windowID]);
        }
            

        markSNP_advance(window_start, x_start, x_length, y_start, y_length, 
        &(overlap_list->list[overlapID].w_list[windowID]), beg_cigar, end_cigar, hap);
    }


    /****************************may have bugs********************************/
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
            hap->snp++;
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
    for (i = 0; i < dumy->length; i++)
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
            beg_cigar = &(overlap_list->list[overlapID].boundary_cigars.buffer[windowID-1]);
        }
        if(windowID < overlap_list->list[overlapID].w_list_length - 1)
        {
            end_cigar = &(overlap_list->list[overlapID].boundary_cigars.buffer[windowID]);
        }


        addSNPtohaplotype_advance(window_start, overlapID, x_start, x_length, 
        y_start, y_length, &(overlap_list->list[overlapID].w_list[windowID]),
        beg_cigar, end_cigar, hap, snp_threshold, x_string, y_string);                
    }

    RsetInitHaplotypeEvdienceFlag(hap, first_snp, last_snp + 1 - first_snp);
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

void debug_hap_information(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, haplotype_evdience_alloc* hap,
                        Correct_dumy* dumy)
{
    int i, overlapID, y_ID, y_Strand;
    long long x_start;
    long long x_length; 
    char* x_string;
    char* y_string;
    long long y_start;
    long long y_length;

    for (i = 0; i < hap->length; i++)
    {
        if(hap->list[i].type < 2)
        {
            overlapID = hap->list[i].overlapID;
            x_start = hap->list[i].site; 
            y_start = hap->list[i].overlapSite;


            y_ID = overlap_list->list[overlapID].y_id;
            y_Strand = overlap_list->list[overlapID].y_pos_strand;

            recover_UC_Read_sub_region(dumy->overlap_region, y_start, 1, y_Strand, R_INF, y_ID);

            x_string = g_read->seq + x_start;
            y_string = dumy->overlap_region;

            if(y_string[0] != hap->list[i].misBase)
            {
                fprintf(stderr, "y_string[0]: %c, hap->list[i].misBase: %c\n",
                    y_string[0], hap->list[i].misBase);
            }


            if(hap->list[i].type == 0)
            {
                if(x_string[0] != y_string[0])
                {
                    fprintf(stderr, "x_string[0]: %c, y_string[0]: %c\n",
                    x_string[0], y_string[0]);


                }
            }
            else if(hap->list[i].type == 1)
            {
                if(x_string[0] == y_string[0])
                {
                    fprintf(stderr, "x_string[0]: %c, y_string[0]: %c\n",
                    x_string[0], y_string[0]);
                }

            }
        }
    }


    for (i = 0; i < hap->length; i++)
    {
        if(i != 0 && hap->list[i].site  < hap->list[i-1].site)
        {
            fprintf(stderr, "wrong order\n");
        }
    }
}









int debug_split_sub_list(haplotype_evdience_alloc* hap, 
haplotype_evdience* sub_list, long long sub_length, long long num_haplotype)
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
    **/
    if(occ_0 == 0 || occ_1 <= 1)
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
    if(available < threshold)
    {
        return 0;
    }

    ///if we just have one snp, we need to phase it carefully
    if(num_haplotype == 1)
    {
        ///we must have just 1 match and 1 mismatch
        ///any other types are not good
        if(new_0 + max != new_total)
        {
            return 0;
        }
        
        if(filter_snp(new_0, max, new_total) == 0)
        {
            return 0;
        }
    }
    /**
    if(filter_snp(new_0, max, new_total) == 0)
    {
        return 0;
    }
    **/

    


    ///for each calculated snp, find if it is at snp matrix
    for (i = 0; i < hap->available_snp; i++)
    {
        if(hap->snp_stat[i].site == sub_list[0].site)
        {

            
            int j = 0;
            int vectorID = hap->snp_stat[i].id;
            int8_t* vector = Get_SNP_Vector((*hap), vectorID);

            if(hap->snp_stat[i].occ_0 != occ_0)
            {
                fprintf(stderr, "error occ0\n");
            }

            if(hap->snp_stat[i].occ_1 != occ_1_array[max_i])
            {
                fprintf(stderr, "error occ1\n");
            }

            if(hap->snp_stat[i].overlap_num != sub_length)
            {
                fprintf(stderr, "error overlap_num\n");
            }

            if(hap->snp_stat[i].overlap_num != hap->snp_stat[i].occ_0 + 
            hap->snp_stat[i].occ_1 + hap->snp_stat[i].occ_2)
            {
                fprintf(stderr, "error overlap_num\n");
            }

            ///for each element in snp vector, find if it is in calculated dataset
            for (j = 0; j < Get_SNP_Vector_Length((*hap)); j++)
            {
                if(vector[j] != -1)
                {
                    int x_i = 0;
                    for (x_i = 0; x_i < sub_length; x_i++)
                    {
                        if(j == sub_list[x_i].overlapID)
                        {
                            break;
                        }
                    }

                    if(x_i == sub_length)
                    {
                        fprintf(stderr, "error: j: %d\n",j);
                    }
                    else
                    {
                        if(vector[j] == 0 || sub_list[x_i].type == 0)
                        {
                            if(vector[j] != sub_list[x_i].type)
                            {
                                fprintf(stderr, "error: 0: %d\n",j);
                            }
                        }

                        if(vector[j] == 1)
                        {
                            if(sub_list[x_i].type != 1)
                            {
                                fprintf(stderr, "-error: 1: %d\n",j);
                            }


                            if(sub_list[x_i].type == 1 && sub_list[x_i].misBase != s_H[max_i])
                            {
                                fprintf(stderr, "+error: 1: %d\n",j);
                            }
                        }

                        if(vector[j] == 2)
                        {
                            if(sub_list[x_i].type != 2)
                            {
                                if(sub_list[x_i].type == 1 && sub_list[x_i].misBase != s_H[max_i])
                                {
                                    ;
                                }
                                else
                                {
                                    fprintf(stderr, "error: 2: %d\n",j);
                                }
                    
                            }

                        }
                        
         
                    }
                    
                }
            }


            ///for each calculated data, find if it is in snp vector
            for (j = 0; j < sub_length; j++)
            {
                if(vector[sub_list[j].overlapID] != sub_list[j].type)
                {

                    if(vector[sub_list[j].overlapID] == 2 && sub_list[j].type == 1 && sub_list[j].misBase != s_H[max_i])
                    {
                        ;
                    }
                    else
                    {
                            fprintf(stderr, "vector[sub_list[j].site]: %d, sub_list[j].type: %d\n",
                        vector[sub_list[j].overlapID], sub_list[j].type);
                    }
                }
            }
            
            
            break;
        }
    }

    
    if(i == hap->available_snp)
    {
        fprintf(stderr, "error\n");
    }    
    

    /**
    fprintf(stderr, "new_0: %d, occ_0: %d, max: %d, max_i: %d, sub_length: %d, new_total: %d, available: %lf\n", 
    new_0, occ_0, max, max_i, sub_length, new_total, available);
    for (i = 0; i < sub_length; i++)
    {
        
        fprintf(stderr, "i: %d, site: %d, type: %d, char: %c, ID: %d, name: %.*s\n", 
    i, sub_list[i].site, sub_list[i].type, sub_list[i].misBase, sub_list[i].overlapID, 
    Get_NAME_LENGTH((*R_INF), overlap_list->list[sub_list[i].overlapID].y_id), 
    Get_NAME((*R_INF),overlap_list->list[sub_list[i].overlapID].y_id)); 

    }
    fprintf(stderr, "\n");
    **/
    
    
    return 1;
    

    
}













int debug_snp_matrix(haplotype_evdience_alloc* hap)
{
    uint64_t pre_site = (uint64_t)-1;
    uint64_t num_of_snps = 0;
    long long pre_i = -1;
    long long sub_length;
    haplotype_evdience* sub_list;
    long long i;
    long long a_snp = 0;


    ////split reads
    for (i = 0; i < hap->length; i++)
    {
        if(pre_site != hap->list[i].site)
        {
            if(i != 0)
            {
                sub_list = hap->list + pre_i;
                sub_length = i - pre_i;
                ///debug_total_length = debug_total_length + sub_length;
                a_snp += debug_split_sub_list(hap, sub_list, sub_length, hap->snp);
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
        ///debug_total_length = debug_total_length + sub_length;
        a_snp += debug_split_sub_list(hap, sub_list, sub_length, hap->snp);
    }

    
    if(a_snp != hap->available_snp)
    {
        fprintf(stderr, "a_snp: %d, available_snp: %d\n", 
        a_snp, hap->available_snp);
    }
    

}

int split_sub_list(haplotype_evdience_alloc* hap, 
haplotype_evdience* sub_list, long long sub_length, long long num_haplotype, 
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
    if(available < threshold)
    {
        return 0;
    }

    ///new_total is the number of errors here
    new_total = new_total - new_0;
    ///available is the number of selected errors here
    available = max;
    threshold = 0.70;
    available = available/((double)(new_total));
    if(available < threshold)
    {
        return 0;
    }
    
   



    InsertSNPVector(hap, sub_list, sub_length, s_H[max_i], g_read);
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
    int i, j;
    for (i = 0; i < hap->core_snp; i++)
    {
        fprintf(stderr, "core(i): %d, site: %d, occ_0: %d, occ_1: %d, occ_2: %d, score: %d\n", 
    i, hap->snp_stat[i].site, hap->snp_stat[i].occ_0, hap->snp_stat[i].occ_1,
    hap->snp_stat[i].occ_2,
    hap->snp_stat[i].score); 

        int vectorID = hap->snp_stat[i].id;
        int8_t* vector = Get_SNP_Vector((*hap), vectorID);

        for (j = 0; j < hap->overlap; j++)
        {
            if(vector[j] == 0)
            {
                fprintf(stderr, "type: %d, ID: %d\n", vector[j], j);
            }
        }


        for (j = 0; j < hap->overlap; j++)
        {
            if(vector[j] == 1)
            {
                fprintf(stderr, "type: %d, ID: %d\n", vector[j], j);
            }
        }


        for (j = 0; j < hap->overlap; j++)
        {
            if(vector[j] == 2)
            {
                fprintf(stderr, "type: %d, ID: %d\n", vector[j], j);
            }
        }
        
    }
}

void merge_snp_vectors(haplotype_evdience_alloc* hap, int diff_vector_ID)
{
    int8_t *r_vector = Get_Result_SNP_Vector((*hap));
    int vectorLen = Get_SNP_Vector_Length((*hap));
    memset(r_vector, -1, vectorLen);
    hap->result_stat.occ_0 = 0;
    hap->result_stat.occ_1 = 0;

    int8_t* vector;
    int vectorID;
    int i, j;

    for (i = 0; i < hap->core_snp; i++)
    {
        if(i == diff_vector_ID)
        {
            continue;
        }

        vectorID = hap->snp_stat[i].id;
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
                if((vector[j] != -1 && vector[j] != 2 && vector[j] != r_vector[j]))
                {
                    fprintf(stderr, "j: %d, vector[j]: %d, r_vector[j]: %d, hap->core_snp: %d, diff_vector_ID: %d\n",
                    j, vector[j], r_vector[j], hap->core_snp, diff_vector_ID);

                    print_core_snp(hap);
                }
            }
            
            
            
        }
    }

    hap->result_stat.overlap_num = hap->result_stat.occ_0 + hap->result_stat.occ_1;
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
                ///fprintf(stderr, "j: %d\n", j);
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

    for (i = 0; i < hap->core_snp; i++)
    {
        if(i == diff_vector_ID)
        {
            continue;
        }

        vectorID = hap->snp_stat[i].id;
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

    int vectorID, vectorID2;
    int diff_core_vector = 0;
    int diff_vector_ID = -1;
    int8_t *vector, *vector2;
    
    if(hap->core_snp == 0)
    {
        return 0;
    }

    ///sort by weight
    qsort(hap->snp_stat, hap->available_snp, sizeof(SnpStats), cmp_snp_stats);


    // for (j = 0; j < hap->available_snp; j++)
    // {
    //     fprintf(stderr, "j: %d, score: %d\n",  j, hap->snp_stat[j].score);
    // }
    // fprintf(stderr, "\n\n");

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

        
        diff_vector_ID = -1;
        ///first try to merge all vector together
        if(merge_snp_vectors_and_test(hap, -1) == 0)
        {
            for (j = hap->core_snp - 1; j >= 0; j--)
            {
                if(merge_snp_vectors_and_test(hap, j) == 1)
                {
                    diff_vector_ID = j;
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
    for (j = hap->core_snp; j < hap->available_snp; j++)
    {
        vectorID2 = hap->snp_stat[j].id;
        vector2 = Get_SNP_Vector((*hap), vectorID2);
        if(calculate_distance_snp_vector(vector, vector2, Get_SNP_Vector_Length((*hap))) == 0)
        {
            add_to_result_snp_vector(hap, vector2, Get_SNP_Vector_Length((*hap)));
        }
    }

    ///merge_snp_vectors(hap, diff_vector_ID);

    ///for read only have 1 snp, we need a more strict condition
    if (hap->core_snp == 1 && 
    filter_one_snp(hap->result_stat.occ_0 + 1, hap->result_stat.occ_1, 
    hap->result_stat.overlap_num + 1) == 0)
    {
        return 0;
    }


    return 1;
    
}



void print_snp_in_line(haplotype_evdience_alloc* hap)
{
    int j, i;
    uint32_t* column;
        fprintf(stderr, "###########hap->available_snp: %d###########\n", hap->available_snp);
        for (j = 0; j < hap->available_snp; j++)
        {
            fprintf(stderr, "*********j: %d, site: %d, id: %d*********\n", j, hap->snp_stat[j].site, hap->snp_stat[j].id);

            fprintf(stderr, "type(0):\n");

            int vectorID = hap->snp_stat[j].id;
            int8_t* vector = Get_SNP_Vector((*hap), vectorID);
            for (i = 0; i < hap->overlap; i++)
            {
                if(vector[i] == 0)
                {
                    fprintf(stderr, "%3d, ", i);
                }
            }
            fprintf(stderr, "\n");

            fprintf(stderr, "type(1):\n");
            for (i = 0; i < hap->overlap; i++)
            {
                if(vector[i] == 1)
                {
                    fprintf(stderr, "%3d, ", i);
                }
            }
            fprintf(stderr, "\n");
        }

        fprintf(stderr, "***********************\n");
        for (i = 0; i < hap->dp.snp_num; i++)
        {
            fprintf(stderr, "hap->dp.max[%d]: %d, hap->dp.backtrack_length: %d\n", 
            i, hap->dp.max[i], hap->dp.backtrack_length[i]);

            if(hap->dp.backtrack_length[i] != 0)
            {
                column = Get_DP_Backtrack_Column(hap->dp, i);
                for (j = 0; j < hap->dp.backtrack_length[i]; j++)
                {
                    fprintf(stderr, "pre: %d,", column[j]);
                }
                fprintf(stderr, "\n");
            }
            
        }

        fprintf(stderr, "#########################\n\n\n");
}

///if j == -1, print result vector
void print_single_snp(haplotype_evdience_alloc* hap, int j)
{



    int i, vectorID;
    int8_t* vector;
    if(j != -1)
    {
        fprintf(stderr, "*********site: %d, j: %d, id: %d*********\n", 
        hap->snp_stat[j].site, j, hap->snp_stat[j].id);
        vectorID = hap->snp_stat[j].id;
        vector = Get_SNP_Vector((*hap), vectorID);
    }
    else
    {
        fprintf(stderr, "*********result snp*********\n");
        vector = Get_Result_SNP_Vector((*hap));
    }
    

    
    fprintf(stderr, "type(0):\n");
    for (i = 0; i < hap->overlap; i++)
    {
        if(vector[i] == 0)
        {
            fprintf(stderr, "%3d, ", i);
        }
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "type(1):\n");
    for (i = 0; i < hap->overlap; i++)
    {
        if(vector[i] == 1)
        {
            fprintf(stderr, "%3d, ", i);
        }
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "###############\n\n");
}



void Preorder_Merge(uint32_t snpID, haplotype_evdience_alloc* hap, int is_merge) 
{ 
    int vectorID = hap->snp_stat[snpID].id;
    int8_t* vector = Get_SNP_Vector((*hap), vectorID);
    hap->dp.visit[snpID] = 1;
    

    if(is_merge)
    {
        if(hap->snp_stat[snpID].is_homopolymer)
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
            for (j = 0; j < hap->dp.backtrack_length[snpID]; j++)
            {
                if(hap->snp_stat[column[j]].is_homopolymer == 0)
                {
                    add_ID = j;
                }
            }

            for (j = 0; j < hap->dp.backtrack_length[snpID]; j++)
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
            for (j = 0; j < hap->dp.backtrack_length[snpID]; j++)
            {
                Preorder_Merge(column[j], hap, 0); 
            }
        }
    }
}


void generate_result_vector_repeat(haplotype_evdience_alloc* hap, int pathLen)
{
    if(pathLen != hap->dp.current_snp_num)
    {
        fprintf(stderr, "hahah\n");
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
        vectorID = hap->snp_stat[snpID1].id;
        vector = Get_SNP_Vector((*hap), vectorID);
        if(hap->snp_stat[snpID1].is_homopolymer)
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
    

    insert_SNP_IDs_addition(&(hap->dp.SNP_IDs), hap->dp.buffer, pathLen, hap->result_stat.occ_0, hap->result_stat.occ_1,
    hap->result_stat.homopolymer_num, hap->result_stat.non_homopolymer_num);
}



void Preorder_Merge_Advance_Repeat(uint32_t snpID, haplotype_evdience_alloc* hap, int pathLen) 
{ 
    hap->dp.visit[snpID] = 1;
    hap->dp.buffer[pathLen] = snpID;
    pathLen++;

    if(hap->dp.backtrack_length[snpID] == 0)
    {
        ///generate_result_vector_repeat(hap, pathLen);
        insert_SNP_IDs_addition(&(hap->dp.SNP_IDs), hap->dp.buffer, pathLen);
        return;
    }
    else
    {
        uint32_t* column;
        int j;

        column = Get_DP_Backtrack_Column(hap->dp, snpID);
        for (j = 0; j < hap->dp.backtrack_length[snpID]; j++)
        {
            Preorder_Merge_Advance_Repeat(column[j], hap, pathLen); 
        }
    }
}



void generate_result_vector(haplotype_evdience_alloc* hap, int pathLen)
{
    if(pathLen != hap->dp.current_snp_num)
    {
        fprintf(stderr, "hahah\n");
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
        vectorID = hap->snp_stat[snpID1].id;
        vector = Get_SNP_Vector((*hap), vectorID);
        if(hap->snp_stat[snpID1].is_homopolymer)
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
        for (j = 0; j < hap->dp.backtrack_length[snpID]; j++)
        {
            Preorder_Merge_Advance(column[j], hap, pathLen); 
        }
    }
}



int if_snp_vector_useful_v2(haplotype_evdience_alloc* hap,
long long occ_0, long long occ_1, uint32_t* SNPs, long long SNPsLen)
{

    double occ_1_coverage_low = (occ_0 + occ_1) * 0.3;
    
    if(occ_1 == 0 || occ_0 == 0)
    {
        return 0;
    }


    if(occ_1 >= occ_1_coverage_low && occ_0 >= occ_1_coverage_low)
    {
        return 1;
    }
    else if(occ_1 >= 5 && occ_0 >= 5)
    {
        return 1;
    }
    else if(occ_1 >= 2 && occ_0 >= 2 && SNPsLen >= 2)
    {
        /**
        int nearsnp;
        int non_nearsnps;
        count_nearby_snps(hap, SNPs, SNPsLen, &nearsnp, &non_nearsnps);
        if(non_nearsnps > 0)
        {
            return 1;
        }
        **/
       return 1;
    }

    return 0;
}


int if_snp_vector_useful(haplotype_evdience_alloc* hap,
long long occ_0, long long occ_1, uint32_t* SNPs, long long SNPsLen)
{

    double occ_1_coverage_low = (occ_0 + occ_1) * 0.3;
    
    if(occ_1 == 0 || occ_0 == 0)
    {
        return 0;
    }


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
    int current_score;
    for (j = 0; j < SNPLen; j++)
    {
        snpID1 = SNPs[j];
        vectorID = hap->snp_stat[snpID1].id;
        vector = Get_SNP_Vector((*hap), vectorID);
        if(hap->snp_stat[snpID1].is_homopolymer)
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
        vectorID = hap->snp_stat[snpID].id;
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
            if( hap->snp_stat[snpID].site >= overlap_list->list[j].x_pos_s
                &&
                hap->snp_stat[snpID].site <= overlap_list->list[j].x_pos_e)
              {
                  overlap_list->list[j].strong = 1;
              }
            /****************************may have bugs********************************/

        }
    }
}


void output_reads_phase(haplotype_evdience_alloc* hap, uint32_t* SNPs, long long SNPsLen, 
overlap_region_alloc* overlap_list, All_reads* R_INF)
{
    long long i, j, snpID, vectorID, overlapLen;
    int8_t *vector;

    for (i = 0; i < SNPsLen; i++)
    {
        snpID = SNPs[i];
        vectorID = hap->snp_stat[snpID].id;
        vector = Get_SNP_Vector((*hap), vectorID);

        fprintf(stderr, "i: %d, site: %d, Get_SNP_Vector_Length((*hap)): %d\n", 
        i, hap->snp_stat[snpID].site, Get_SNP_Vector_Length((*hap)));
    
        fprintf(stderr, "flag 1\n");
        for (j = 0; j < Get_SNP_Vector_Length((*hap)); j++)
        {
            if(vector[j] == 1)
            {
            
                fprintf(stderr, "%.*s\n", Get_NAME_LENGTH((*R_INF), overlap_list->list[j].y_id),
                Get_NAME((*R_INF), overlap_list->list[j].y_id));
            
            }
        }


        fprintf(stderr, "flag 0\n");
        for (j = 0; j < Get_SNP_Vector_Length((*hap)); j++)
        {

            if(vector[j] == 0)
            {
            
                fprintf(stderr, "%.*s\n", Get_NAME_LENGTH((*R_INF), overlap_list->list[j].y_id),
                Get_NAME((*R_INF), overlap_list->list[j].y_id));
            
            }
        }

    }
}

void remove_reads_debug(haplotype_evdience_alloc* hap, uint32_t* SNPs, long long SNPsLen, overlap_region_alloc* overlap_list)
{
    fprintf(stderr, "SNPsLen: %d\n", SNPsLen);
    long long i, j, snpID, vectorID, overlapLen;
    int8_t *vector;

    for (i = 0; i < SNPsLen; i++)
    {
        snpID = SNPs[i];
        vectorID = hap->snp_stat[snpID].id;
        vector = Get_SNP_Vector((*hap), vectorID);
        fprintf(stderr, "i: %d, snpID:%d, SNPsLen: %d, available_snp: %d, snp_stat[snpID].site: %d\n", 
        i, snpID, SNPsLen, hap->available_snp, hap->snp_stat[snpID].site);

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
            if( hap->snp_stat[snpID].site >= overlap_list->list[j].x_pos_s
                &&
                hap->snp_stat[snpID].site <= overlap_list->list[j].x_pos_e)
              {
                  overlap_list->list[j].strong = 1;
              }

              fprintf(stderr, "j: %d, x_pos_s: %d, x_pos_e: %d, strong: %d, is_match: %d", j, overlap_list->list[j].x_pos_s, 
            overlap_list->list[j].x_pos_e, overlap_list->list[j].strong, 
            overlap_list->list[j].is_match);
            fprintf(stderr, "****************y: %.*s****************\n", 
            Get_NAME_LENGTH(R_INF, overlap_list->list[j].y_id), 
            Get_NAME(R_INF, overlap_list->list[j].y_id));
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
        
        
        ///if(overlap_list->list[0].x_id == 83735)
        ///if(overlap_list->list[0].x_id == 83739)
        // if(overlap_list->list[0].x_id == 1185538)
        // {
        //     fprintf(stderr, "SNPLen: %d, x_id: %d, hap->snp_stat[snpID].site: %d, occ_0: %d, occ_1: %d, occ_2: %d, overlap_num: %d\n", 
        //     SNPLen, overlap_list->list[0].x_id, hap->snp_stat[snpID].site,
        //     hap->snp_stat[snpID].occ_0, hap->snp_stat[snpID].occ_1, hap->snp_stat[snpID].occ_2,
        //     hap->snp_stat[snpID].overlap_num);
        // }
        
        

        ///check all overlaps
        for (j = 0; j < Get_SNP_Vector_Length((*hap)); j++)
        {
            /****************************may have bugs********************************/
            if( hap->snp_stat[snpID].site >= overlap_list->list[j].x_pos_s
                &&
                hap->snp_stat[snpID].site <= overlap_list->list[j].x_pos_e)
              {
                  overlap_list->list[j].strong = 1;
              }
            /****************************may have bugs********************************/
        }
    }
}


void process_repeat_snps(haplotype_evdience_alloc* hap, int coverage, overlap_region_alloc* overlap_list)
{
    int i, snpID, vectorID, flag;
    int8_t *vector;


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



    ///print_snp_in_line(hap);
    // fprintf(stderr, "-:overlap_list->mapped_overlaps: %d\n", overlap_list->mapped_overlaps);
    /**
    if(overlap_list->mapped_overlaps > coverage * 1.6)
    {
        for (i = 0; i < hap->dp.SNP_IDs.IDs_length; i++)
        {
            snp_ids = hap->dp.SNP_IDs.buffer + hap->dp.SNP_IDs.IDs[i].beg;
            length = hap->dp.SNP_IDs.IDs[i].end -hap->dp.SNP_IDs.IDs[i].beg + 1;
            if(hap->dp.SNP_IDs.IDs[i].is_remove == 0 && 
            if_snp_vector_useful(hap, hap->dp.SNP_IDs.IDs[i].occ_0, hap->dp.SNP_IDs.IDs[i].occ_1,
            occ_1_threshold_low, coverage, snp_ids, length, 1))
            {
                //fprintf(stderr, "i: %d \n", i);
                remove_reads(hap, snp_ids, length, overlap_list);
                
            }
        }
    }
    **/

    // fprintf(stderr, "-:overlap_list->mapped_overlaps: %d\n\n\n", overlap_list->mapped_overlaps);

    

    
}

void process_repeat_snps_debug(haplotype_evdience_alloc* hap, int coverage, 
overlap_region_alloc* overlap_list, All_reads* R_INF)
{
    int i, snpID, vectorID, flag;
    int8_t *vector;
    long long occ_1_threshold_low;
    long long occ_1_threshold_up;
    occ_1_threshold_low = 0;


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


            
            if(overlap_list->list[0].x_id == 5405)
            {
                fprintf(stderr, "snpid length: %d, occ_0: %d, occ_1: %d\n", 
                length,
                hap->result_stat.occ_0,
                hap->result_stat.occ_1);
                int k;
                for (k = 0; k < length; k++)
                {
                    fprintf(stderr, "i: %d, site: %d\n", 
                    i, hap->snp_stat[snp_ids[k]].site);
                }

                for (k = 0; k < Get_SNP_Vector_Length((*hap)); k++)
                {
                    fprintf(stderr, "flag 0\n");
                    if(Get_Result_SNP_Vector((*hap))[k] == 0)
                    {
                        fprintf(stderr, "%.*s\n", 
                        Get_NAME_LENGTH((*R_INF),overlap_list->list[k].y_id),
                        Get_NAME((*R_INF),overlap_list->list[k].y_id));
                    }

                    fprintf(stderr, "flag 1\n");
                    if(Get_Result_SNP_Vector((*hap))[k] == 1)
                    {
                        fprintf(stderr, "%.*s\n", 
                        Get_NAME_LENGTH((*R_INF),overlap_list->list[k].y_id),
                        Get_NAME((*R_INF),overlap_list->list[k].y_id));
                    }
                }

                
                
            }
            
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



void debug_repeat_vector(haplotype_evdience_alloc* hap)
{
    int j, i, snpID, vectorID, flag;
    int8_t *vector;
    

    // if(memcmp(hap->dp.SNP_IDs.buffer + hap->dp.SNP_IDs.IDs[hap->dp.SNP_IDs.max_snp_id].beg, 
    // hap->dp.max_buffer, 
    // sizeof(uint32_t) *(hap->dp.SNP_IDs.IDs[hap->dp.SNP_IDs.max_snp_id].end - 
    // hap->dp.SNP_IDs.IDs[hap->dp.SNP_IDs.max_snp_id].beg + 1)))
    // {
    //     fprintf(stderr, "error1\n");
    // }

    // if(hap->dp.SNP_IDs.IDs[hap->dp.SNP_IDs.max_snp_id].end - 
    // hap->dp.SNP_IDs.IDs[hap->dp.SNP_IDs.max_snp_id].beg + 1 != 
    // hap->dp.max_snp_num)
    // {
    //     fprintf(stderr, "error2\n");
    // }



    

    uint32_t* snp_ids;
    long long length;
    for (i = 0; i < hap->dp.SNP_IDs.IDs_length; i++)
    {
        snp_ids = hap->dp.SNP_IDs.buffer + hap->dp.SNP_IDs.IDs[i].beg;
        length = hap->dp.SNP_IDs.IDs[i].end -hap->dp.SNP_IDs.IDs[i].beg + 1;

        ////first clear result snp
        vector = Get_Result_SNP_Vector((*hap));
        memset(vector, -1, Get_SNP_Vector_Length((*hap)));

        long long non_hom = 0;
        long long hom = 0;
        for (j = 0; j < length; j++)
        {
            ///note here is hap->dp.max_buffer instead of hap->dp.buffer
            snpID = snp_ids[j];
            vectorID = hap->snp_stat[snpID].id;
            vector = Get_SNP_Vector((*hap), vectorID);


            if(hap->snp_stat[snpID].is_homopolymer)
            {
                hom++;
            }
            else
            {
                non_hom++;
            }
            

            if((flag = debug_add_to_result_snp_vector(hap, vector, Get_SNP_Vector_Length((*hap))))!= -1)
            {
                fprintf(stderr, "incompatible snp vector....\n");
                exit(0);
            }
        }

        vector = Get_Result_SNP_Vector((*hap));
        long long occ_0 = 0;
        long long occ_1 = 0;
        for (j = 0; j < Get_SNP_Vector_Length((*hap)); j++)
        {
            if(vector[j] == 0)
            {
                occ_0++;
            }

            if(vector[j] == 1)
            {
                occ_1++;
            }
        }

        if(hom != hap->dp.SNP_IDs.IDs[i].homopolymer_num)
        {
            fprintf(stderr, "error\n");
        }

        if(non_hom != hap->dp.SNP_IDs.IDs[i].non_homopolymer_num)
        {
            fprintf(stderr, "error\n");
        }

        if(occ_0 != hap->dp.SNP_IDs.IDs[i].occ_0)
        {
            fprintf(stderr, "error\n");
        }

        if(occ_1 != hap->dp.SNP_IDs.IDs[i].occ_1)
        {
            fprintf(stderr, "error\n");
        }
        


    }
}


int generate_haplotypes_DP_back(haplotype_evdience_alloc* hap, overlap_region_alloc* overlap_list, All_reads* R_INF, long long rLen, 
int force_repeat)
{
    int j, i;

    int vectorID, vectorID2;
    int diff_core_vector = 0;
    int diff_vector_ID = -1;
    int8_t *vector, *vector2;

    
    if(hap->available_snp == 0)
    {
        return 0;
    }

    

    ///if hap->available_snp == 1, the following codes would have bugs
    ///filter snps that are highly likly false 
    if(hap->available_snp > 1)
    {
        i = 0;
        ///if a snp is very near to others, it should not be a real snp
        for (j = 0; j < hap->available_snp; j++)
        {
            // if(hap->snp_stat[j].occ_1 == 1)
            // {
            //     fprintf(stderr, "***\n");
            // }
            if(j > 0 && j < hap->available_snp - 1)
            {
                if(hap->snp_stat[j].site != hap->snp_stat[j - 1].site + 1
                    &&
                hap->snp_stat[j].site + 1 != hap->snp_stat[j + 1].site)
                {
                    hap->snp_stat[i] = hap->snp_stat[j];
                    i++;
                }

            }
            else if(j == 0)
            {
                if(hap->snp_stat[j].site + 1 != hap->snp_stat[j + 1].site)
                {
                    hap->snp_stat[i] = hap->snp_stat[j];
                    i++;
                }
            }
            else
            {
                if(hap->snp_stat[j].site != hap->snp_stat[j - 1].site + 1)
                {
                    hap->snp_stat[i] = hap->snp_stat[j];
                    i++;
                }
            }
        }
        hap->available_snp = i;
    }


    
    


    int flag;
    long long overlap_length, total_read, unuseful_read, last_j, last_j_ID, last_j_flag;
    total_read = unuseful_read = 0;
    ///check if any read may be conflict with others
    for (i = 0; i < overlap_list->length; i++)
    {
        overlap_length = overlap_list->list[i].x_pos_e - overlap_list->list[i].x_pos_s + 1;
        if (overlap_list->list[i].is_match == 1)
        {
            total_read++;
            flag = -1;
            for (j = 0; j < hap->available_snp; j++)
            {
                vectorID = hap->snp_stat[j].id;
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
                        last_j = hap->snp_stat[j].site;
                        last_j_ID = j;
                        last_j_flag = vector[i];
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


            if(flag == 3)
            {
                unuseful_read++;
                for (j = 0; j < hap->available_snp; j++)
                {
                    vectorID = hap->snp_stat[j].id;
                    vector = Get_SNP_Vector((*hap), vectorID);
                    

                    
                    if(vector[i] == 0)
                    {
                        hap->snp_stat[j].occ_0--;
                        hap->snp_stat[j].occ_2++;
                    }
                    else if(vector[i] == 1)
                    {
                        hap->snp_stat[j].occ_1--;
                        hap->snp_stat[j].occ_2++;
                    }
                    else if(vector[i] != 2)
                    {
                        hap->snp_stat[j].occ_2++;
                    }
                    

                    vector[i] = 2;
                }
                // if(overlap_list->list[i].is_match == 0)
                // {
                //     fprintf(stderr, "error\n");
                // }
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
    init_DP_matrix(&(hap->dp), hap->available_snp);

    long long equal_best = 0;
    uint32_t* column;
    long long column_length;

    

    for (i = 0; i < hap->available_snp; i++)
    {
        ///vector of snp i
        vectorID = hap->snp_stat[i].id;
        vector = Get_SNP_Vector((*hap), vectorID);
        hap->dp.visit[i] = 0;
        hap->dp.max[i] = 1;
        hap->dp.backtrack_length[i] = 0;
        equal_best = 0;
        column = Get_DP_Backtrack_Column(hap->dp, i);
        column_length = Get_DP_Backtrack_Column_Length(hap->dp, i);

        for (j = 0; j < i; j++)
        {
            ///vector of snp j
            vectorID2 = hap->snp_stat[j].id;
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

    for (i = 0; i < hap->available_snp; i++)
    {
        tmp_mode = hap->dp.max[i];
        tmp_mode = tmp_mode << 32;
        tmp_mode = tmp_mode | (uint64_t)(i);
        hap->dp.max_for_sort[i] = tmp_mode;
    }
    
    qsort(hap->dp.max_for_sort, hap->available_snp, sizeof(uint64_t), cmp_max_DP);


    int snpID;
    int group_num = 0;
    ///the minmum snp_num is 1
    hap->dp.max_snp_num = 0;
    hap->dp.max_score = -2;
    
    
    
    

    //repeat
    if(overlap_list->mapped_overlaps_length > Coverage_Threshold(coverage, rLen) || force_repeat)
    {
        for (i = 0; i < hap->available_snp; i++)
        {
            snpID = Get_Max_DP_ID(hap->dp.max_for_sort[i]);
            if(hap->dp.visit[snpID] == 0)
            {
                hap->dp.current_snp_num = Get_Max_DP_Value(hap->dp.max_for_sort[i]);
                Preorder_Merge_Advance_Repeat(snpID, hap, 0);
            }
        }
    }
    else  //non-repeat
    {
        for (i = 0; i < hap->available_snp; i++)
        {
            snpID = Get_Max_DP_ID(hap->dp.max_for_sort[i]);
            if(hap->dp.visit[snpID] == 0)
            {
                hap->dp.current_snp_num = Get_Max_DP_Value(hap->dp.max_for_sort[i]);
                Preorder_Merge_Advance(snpID, hap, 0);
            }
        }
    }

    /**
    if(memcmp("m64016_190918_162737/49678749/ccs", 
    Get_NAME((*R_INF), overlap_list->list[0].x_id), 
    Get_NAME_LENGTH((*R_INF), overlap_list->list[0].x_id)) == 0)
    {
        fprintf(stderr, "x_id: %d, mapped_overlaps_length: %d\n",
        overlap_list->list[0].x_id,
        overlap_list->mapped_overlaps_length);
        fprintf(stderr, "coverage: %d\n",
        coverage);
        fprintf(stderr, "rLen: %d\n",
        rLen);
        fprintf(stderr, "Coverage_Threshold(coverage, rLen): %f\n",
        Coverage_Threshold(coverage, rLen));
        fprintf(stderr, "max_snp_num: %d\n",
        hap->dp.max_snp_num);
    }
    **/
    
    
    
    
    

    if(overlap_list->mapped_overlaps_length > Coverage_Threshold(coverage, rLen) || force_repeat)
    {
        ///debug_repeat_vector(hap);
        ///process_repeat_snps_debug(hap, coverage, overlap_list, R_INF);
        process_repeat_snps(hap, coverage, overlap_list);

        return 1;
        
    }
    else if(hap->dp.max_snp_num > 0)
    {
        /**
        if(memcmp("m64011_190329_072846/59507330/ccs", 
        Get_NAME((*R_INF), overlap_list->list[0].x_id), 
        Get_NAME_LENGTH((*R_INF), overlap_list->list[0].x_id)) == 0)
        {

            output_reads_phase(hap, hap->dp.max_buffer, hap->dp.max_snp_num, 
            overlap_list, R_INF);
        }
        **/
    


    
        remove_reads(hap, hap->dp.max_buffer, hap->dp.max_snp_num, overlap_list);
        return 1;
    }
    else
    {
        return 0;
    }
    

    

    /**
    vector = Get_Result_SNP_Vector((*hap));
    memset(vector, -1, Get_SNP_Vector_Length((*hap)));
    vector = Get_Result_SNP_Vector((*hap));
    for (i = 0; i < hap->available_snp; i++)
    {
        int debug_i = Get_Max_DP_ID(hap->dp.max_for_sort[i]);
        int round = hap->dp.max[debug_i];
        ///fprintf(stderr, "round: %d\n", round);
        if(round > 1)
        {
            vector = Get_Result_SNP_Vector((*hap));
            memset(vector, -1, Get_SNP_Vector_Length((*hap)));
            

            if(hap->dp.backtrack_length[debug_i] < 1)
            {
                fprintf(stderr, "error\n");
            }

        
            
            while (round > 0)
            {
                
                
                
                vectorID2 = hap->snp_stat[debug_i].id;
                vector2 = Get_SNP_Vector((*hap), vectorID2);
                int flag;
                if((flag = debug_add_to_result_snp_vector(hap, vector2, Get_SNP_Vector_Length((*hap))))!= -1)
                {
                    fprintf(stderr, "flag: %d, debug_i: %d, i: %d, x_name: %.*s\n", 
                    flag, debug_i, i, Get_NAME_LENGTH((*R_INF), overlap_list->list[0].x_id), Get_NAME((*R_INF),overlap_list->list[0].x_id));
                    print_snp_in_line(hap);
                }

                if(round != hap->dp.max[debug_i])
                {
                    fprintf(stderr, "error: %d\n", round);
                }


               if(hap->dp.backtrack_length[debug_i] != 0)
                {
                    column = Get_DP_Backtrack_Column(hap->dp, debug_i);
                    debug_i = column[0];
                }
                else if(round != 1)
                {
                    fprintf(stderr, "round: %d\n", round);
                }
                


                round--;
            }
            
        }
    }
    vector = Get_Result_SNP_Vector((*hap));
    memset(vector, -1, Get_SNP_Vector_Length((*hap)));
    

    vector = Get_Result_SNP_Vector((*hap));
    memset(vector, -1, Get_SNP_Vector_Length((*hap)));
    vector = Get_Result_SNP_Vector((*hap));
    for (i = 0; i < hap->available_snp; i++)
    {
        if(hap->dp.max[i] > 1)
        {
            vector = Get_Result_SNP_Vector((*hap));
            memset(vector, -1, Get_SNP_Vector_Length((*hap)));

            if(hap->dp.backtrack_length[i] < 1)
            {
                fprintf(stderr, "error\n");
            }

        
            int debug_i = i;
            int round = hap->dp.max[i];
            while (round > 0)
            {
                
                
                
                vectorID2 = hap->snp_stat[debug_i].id;
                vector2 = Get_SNP_Vector((*hap), vectorID2);
                int flag;
                if((flag = debug_add_to_result_snp_vector(hap, vector2, Get_SNP_Vector_Length((*hap))))!= -1)
                {
                    fprintf(stderr, "flag: %d, debug_i: %d, i: %d, x_name: %.*s\n", 
                    flag, debug_i, i, Get_NAME_LENGTH((*R_INF), overlap_list->list[0].x_id), Get_NAME((*R_INF),overlap_list->list[0].x_id));
                    print_snp_in_line(hap);
                }

                if(round != hap->dp.max[debug_i])
                {
                    fprintf(stderr, "error: %d\n", round);
                }


               if(hap->dp.backtrack_length[debug_i] != 0)
                {
                    column = Get_DP_Backtrack_Column(hap->dp, debug_i);
                    debug_i = column[0];
                }
                else if(round != 1)
                {
                    fprintf(stderr, "round: %d\n", round);
                }
                


                round--;
            }
            
        }
    }
    vector = Get_Result_SNP_Vector((*hap));
    memset(vector, -1, Get_SNP_Vector_Length((*hap)));
    **/
    /**
    if(hap->available_snp > 4)
    // if(memcmp("m54334_180924_221206/48759269/ccs", Get_NAME((*R_INF),overlap_list->list[0].x_id), 
    // Get_NAME_LENGTH((*R_INF), overlap_list->list[0].x_id)) == 0)
    // if(memcmp("m54328_180922_235017/65536381/ccs", Get_NAME((*R_INF),overlap_list->list[0].x_id), 
    // Get_NAME_LENGTH((*R_INF), overlap_list->list[0].x_id)) == 0)
    {
        fprintf(stderr, "\n\n\n###########hap->available_snp: %d###########\n", hap->available_snp);

        fprintf(stderr, "x_name: %.*s\n", 
             Get_NAME_LENGTH((*R_INF), overlap_list->list[0].x_id), Get_NAME((*R_INF),overlap_list->list[0].x_id));

        for (j = 0; j < hap->available_snp; j++)
        {
            fprintf(stderr, "*********site: %d, j: %d, id: %d*********\n", 
            hap->snp_stat[j].site, j, hap->snp_stat[j].id);

            fprintf(stderr, "type(0):\n");

            int vectorID = hap->snp_stat[j].id;
            int8_t* vector = Get_SNP_Vector((*hap), vectorID);
            for (i = 0; i < hap->overlap; i++)
            {
                if(vector[i] == 0)
                {
                    fprintf(stderr, "%3d, ", i);
                }
            }
            fprintf(stderr, "\n");

            fprintf(stderr, "type(1):\n");
            for (i = 0; i < hap->overlap; i++)
            {
                if(vector[i] == 1)
                {
                    fprintf(stderr, "%3d, ", i);
                }
            }
            fprintf(stderr, "\n");

            ///if(j == 21 || j == 23)
            // {
            //     for (i = 0; i < hap->overlap; i++)
            //     {
            //         if(vector[i] == 1)
            //         {
            //             fprintf(stderr, "1: i: %d, %.*s\n", i, Get_NAME_LENGTH((*R_INF), overlap_list->list[i].y_id), 
            //             Get_NAME((*R_INF),overlap_list->list[i].y_id));
            //         }
            //     }
            // }
        }

        fprintf(stderr, "***********************\n");
        for (i = 0; i < hap->dp.snp_num; i++)
        {
            fprintf(stderr, "hap->dp.max[%d]: %d, hap->dp.backtrack_length: %d\n", i, hap->dp.max[i], hap->dp.backtrack_length[i]);
            if(hap->dp.backtrack_length[i] != 0)
            {
                column = Get_DP_Backtrack_Column(hap->dp, i);
                for (j = 0; j < hap->dp.backtrack_length[i]; j++)
                {
                    fprintf(stderr, "pre: %d,", column[j]);
                }
                fprintf(stderr, "\n");
            }
            
        }
    }
    
    **/
    



    /**
    fprintf(stderr, "overlap_list->length: %u, total_read: %u, unuseful_read: %u\n", overlap_list->length, total_read, unuseful_read);
    for (i = 0; i < overlap_list->length; i++)
    {
        overlap_length = overlap_list->list[i].x_pos_e - overlap_list->list[i].x_pos_s + 1;
        if (overlap_length * OVERLAP_THRESHOLD <=  overlap_list->list[i].align_length)
        {
            total_read++;
            flag = -1;
            for (j = 0; j < hap->available_snp; j++)
            {
                vectorID = hap->snp_stat[j].id;
                vector = Get_SNP_Vector((*hap), vectorID);

                if(vector[i] != 2 && flag == 2)
                {
                    fprintf(stderr, "hahahah\n");
                    break;
                }

                if(vector[i] == 2)
                {
                    flag = 2;
                }
            }
        }
    }
    **/

    

}

void lable_large_indels(overlap_region_alloc* overlap_list, All_reads* R_INF, long long read_length,
Correct_dumy* dumy)
{
    long long i, j;
    long long cigar_i, operation, operationLen;
    int is_delete = 0;
    CIGAR* cigar;
    for (i = 0; i < overlap_list->length; i++)
    {
        ///should has at least 3 windows for this overlap
        if (overlap_list->list[i].is_match == 1 && overlap_list->list[i].w_list_length >= 3)
        {
            ///here w_list_length >= 3
            ///skip the first and last window
            for (j = 1; j < overlap_list->list[i].w_list_length - 1; j++)
            {
                ///this window is not matched, it seems to have large difference
                if(overlap_list->list[i].w_list[j].y_end == -1)
                {
                    overlap_list->list[i].is_match = 100;
                    is_delete = 1;
                    goto end_rem;
                }

                cigar = &(overlap_list->list[i].w_list[j].cigar);
                ///if there are <=2 cigar elements, skip it
                if(cigar->length < 3)
                {
                    continue;
                }
                ///skip the first and last cigar elements
                for (cigar_i = 1; cigar_i < cigar->length - 1; cigar_i++)
                {
                    operation = cigar->C_C[cigar_i];
                    operationLen = cigar->C_L[cigar_i];
                    
                    if(operationLen <= 5)
                    {
                        continue;
                    }
                    ///>=6 bp deletion or insertion 
                    if(operation == 2 || operation == 3)
                    {
                        overlap_list->list[i].is_match = 100;
                        is_delete = 1;
                        goto end_rem;
                    }
                }
            }
        }

        end_rem:
        overlap_list->list[i].w_list_length >= 3;

        
    }


    if(is_delete == 1)
    {
        long long overlap_length;
        long long window_start, window_end;
        Window_Pool w_inf;
        init_Window_Pool(&w_inf, read_length, WINDOW, TAIL_LENGTH);
        int flag = 0;
        long long realLen, realLen_100;
        int to_recover = 0;
        while(get_Window(&w_inf, &window_start, &window_end) && flag != -2)
        {
            dumy->length = 0;
            dumy->lengthNT = 0;
            ///flag返回的是重叠数量
            ///dumy->length返回的是有效完全重叠的数量
            ///dumy->lengthNT返回的是有效不完全重叠的数量
            ///return overlaps that is overlaped with [window_start, window_end]
            flag = get_available_fully_covered_interval(window_start, window_end, 
            overlap_list, dumy, &realLen, &realLen_100);


            switch (flag)
            {
                case 1:    ///找到匹配
                    break;
                case 0:    ///没找到匹配
                    break;
                case -2: ///下一个window也不会存在匹配, 直接跳出
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
            for (i = 0; i < overlap_list->length; i++)
            {
                if (overlap_list->list[i].is_match == 100)
                {
                    overlap_list->list[i].is_match = 1;
                }
            }
        }
    }


    for (i = 0; i < overlap_list->length; i++)
    {
        if (overlap_list->list[i].is_match == 1)
        {
            overlap_list->list[i].without_large_indel = 1;
        }

        if (overlap_list->list[i].is_match == 100)
        {
            overlap_list->list[i].is_match = 1;
            overlap_list->list[i].without_large_indel = 0;
        }
    }

}

int debug_print_snp_stat(char* name, haplotype_evdience_alloc* hap, overlap_region_alloc* overlap_list, All_reads* R_INF)
{
    if(overlap_list->length > 0 &&
       memcmp(name, Get_NAME((*R_INF), overlap_list->list[0].x_id), 
       Get_NAME_LENGTH((*R_INF), overlap_list->list[0].x_id)) == 0)
    {
        fprintf(stderr, "\n%s, available_snp: %d\n", name, hap->available_snp);
        int i;
        for (i = 0; i < hap->available_snp; i++)
        {
            fprintf(stderr, "site: %d, occ_0: %d, occ_1: %d, occ_2: %d\n", 
            hap->snp_stat[i].site, hap->snp_stat[i].occ_0,
            hap->snp_stat[i].occ_1, hap->snp_stat[i].occ_2);
        }
    }  
}

int generate_haplotypes_DP(haplotype_evdience_alloc* hap, overlap_region_alloc* overlap_list, All_reads* R_INF, long long rLen, 
int force_repeat)
{
    int j, i;

    int vectorID, vectorID2;
    int diff_core_vector = 0;
    int diff_vector_ID = -1;
    int8_t *vector, *vector2;

    
    if(hap->available_snp == 0)
    {
        return 0;
    }

    ///debug_print_snp_stat("m64013_190324_024932/23660629/ccs", hap, overlap_list, R_INF);

    
    ///if hap->available_snp == 1, the following codes would have bugs
    ///filter snps that are highly likly false 
    if(hap->available_snp > 1)
    {
        i = 0;
        ///if a snp is very near to others, it should not be a real snp
        for (j = 0; j < hap->available_snp; j++)
        {
            if(j > 0 && j < hap->available_snp - 1)
            {
                if(hap->snp_stat[j].site != hap->snp_stat[j - 1].site + 1
                    &&
                hap->snp_stat[j].site + 1 != hap->snp_stat[j + 1].site)
                {
                    hap->snp_stat[i] = hap->snp_stat[j];
                    i++;
                }

            }
            else if(j == 0)
            {
                if(hap->snp_stat[j].site + 1 != hap->snp_stat[j + 1].site)
                {
                    hap->snp_stat[i] = hap->snp_stat[j];
                    i++;
                }
            }
            else
            {
                if(hap->snp_stat[j].site != hap->snp_stat[j - 1].site + 1)
                {
                    hap->snp_stat[i] = hap->snp_stat[j];
                    i++;
                }
            }
        }
        hap->available_snp = i;
    }


    
    


    int flag;
    long long overlap_length, total_read, unuseful_read, last_j, last_j_ID, last_j_flag;
    total_read = unuseful_read = 0;
    ///check if any read may be conflict with others
    for (i = 0; i < overlap_list->length; i++)
    {
        overlap_length = overlap_list->list[i].x_pos_e - overlap_list->list[i].x_pos_s + 1;
        if (overlap_list->list[i].is_match == 1)
        {
            total_read++;
            flag = -1;
            for (j = 0; j < hap->available_snp; j++)
            {
                vectorID = hap->snp_stat[j].id;
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
                        last_j = hap->snp_stat[j].site;
                        last_j_ID = j;
                        last_j_flag = vector[i];
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


            if(flag == 3)
            {
                unuseful_read++;
                for (j = 0; j < hap->available_snp; j++)
                {
                    vectorID = hap->snp_stat[j].id;
                    vector = Get_SNP_Vector((*hap), vectorID);
                    

                    
                    if(vector[i] == 0)
                    {
                        hap->snp_stat[j].occ_0--;
                        hap->snp_stat[j].occ_2++;
                    }
                    else if(vector[i] == 1)
                    {
                        hap->snp_stat[j].occ_1--;
                        hap->snp_stat[j].occ_2++;
                    }
                    else if(vector[i] != 2)
                    {
                        hap->snp_stat[j].occ_2++;
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
    init_DP_matrix(&(hap->dp), hap->available_snp);

    long long equal_best = 0;
    uint32_t* column;
    long long column_length;

    

    for (i = 0; i < hap->available_snp; i++)
    {
        ///vector of snp i
        vectorID = hap->snp_stat[i].id;
        vector = Get_SNP_Vector((*hap), vectorID);
        hap->dp.visit[i] = 0;
        hap->dp.max[i] = 1;
        hap->dp.backtrack_length[i] = 0;
        equal_best = 0;
        column = Get_DP_Backtrack_Column(hap->dp, i);
        column_length = Get_DP_Backtrack_Column_Length(hap->dp, i);

        for (j = 0; j < i; j++)
        {
            ///vector of snp j
            vectorID2 = hap->snp_stat[j].id;
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

    for (i = 0; i < hap->available_snp; i++)
    {
        tmp_mode = hap->dp.max[i];
        tmp_mode = tmp_mode << 32;
        tmp_mode = tmp_mode | (uint64_t)(i);
        hap->dp.max_for_sort[i] = tmp_mode;
    }
    
    qsort(hap->dp.max_for_sort, hap->available_snp, sizeof(uint64_t), cmp_max_DP);


    int snpID;
    int group_num = 0;
    ///the minmum snp_num is 1
    hap->dp.max_snp_num = 0;
    hap->dp.max_score = -2;
    
    
    
    


    for (i = 0; i < hap->available_snp; i++)
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
    if(hap->available_snp > 0)
    {
        process_repeat_snps(hap, coverage, overlap_list);
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

inline long long snp_occ_in_one_read(haplotype_evdience_alloc* hap, overlap_region_alloc* overlap_list, 
long long readID)
{
    long long i;
    long long vectorID;
    int8_t *vector;

    long long snp_occ = 0;
    for (i = 0; i < hap->available_snp; i++)
    {
        vectorID = hap->snp_stat[i].id;
        vector = Get_SNP_Vector((*hap), vectorID);
        if(vector[readID] == 1)
        {
            snp_occ++;
        }
    }

    return snp_occ;
}

inline void remove_read_from_snps(haplotype_evdience_alloc* hap, overlap_region_alloc* overlap_list, 
long long readID)
{
    long long i;
    long long vectorID;
    int8_t *vector;

    for (i = 0; i < hap->available_snp; i++)
    {
        vectorID = hap->snp_stat[i].id;
        vector = Get_SNP_Vector((*hap), vectorID);
        vector[readID] = 2;
    }
    overlap_list->list[readID].is_match = 4;
}

int generate_haplotypes_naive(haplotype_evdience_alloc* hap, overlap_region_alloc* overlap_list, All_reads* R_INF, long long rLen, 
int force_repeat)
{
    int j, i;

    int vectorID, vectorID2;
    int diff_core_vector = 0;
    int diff_vector_ID = -1;
    int8_t *vector, *vector2;

    
    if(hap->available_snp == 0)
    {
        return 0;
    }

    

    ///if hap->available_snp == 1, the following codes would have bugs
    ///filter snps that are highly likly false 
    if(hap->available_snp > 1)
    {
        i = 0;
        ///if a snp is very close to others, it should not be a real snp
        for (j = 0; j < hap->available_snp; j++)
        {
            
            if(j > 0 && j < hap->available_snp - 1)
            {
                if(hap->snp_stat[j].site != hap->snp_stat[j - 1].site + 1
                    &&
                hap->snp_stat[j].site + 1 != hap->snp_stat[j + 1].site)
                {
                    hap->snp_stat[i] = hap->snp_stat[j];
                    i++;
                }

            }
            else if(j == 0)
            {
                if(hap->snp_stat[j].site + 1 != hap->snp_stat[j + 1].site)
                {
                    hap->snp_stat[i] = hap->snp_stat[j];
                    i++;
                }
            }
            else
            {
                if(hap->snp_stat[j].site != hap->snp_stat[j - 1].site + 1)
                {
                    hap->snp_stat[i] = hap->snp_stat[j];
                    i++;
                }
            }
        }
        hap->available_snp = i;
    }


    long long m, snp_occ;
    if(hap->available_snp > 0)
    {
        ///************************debug**************************///
        m = 0;
        for (i = 0; i < hap->available_snp; i++)
        {
            if(check_informative_site(hap, &(hap->snp_stat[i])))
            {
                hap->snp_stat[m] = hap->snp_stat[i];
                m++;
            }
        }
        hap->available_snp = m;
        ///************************debug**************************///



        init_DP_matrix(&(hap->dp), hap->available_snp);

        for (i = 0; i < hap->available_snp; i++)
        {
            hap->dp.max_buffer[i] = i;
        }
        hap->dp.max_snp_num = hap->available_snp;
        remove_reads(hap, hap->dp.max_buffer, hap->dp.max_snp_num, overlap_list);

        return 1;

    }
    else
    {
        return 0;
    }

}


void print_Haplotype(haplotype_evdience_alloc* hap, overlap_region_alloc* overlap_list, All_reads* R_INF)
{
        int j, i;

        
    fprintf(stderr, "\nhap->snp: %d, hap->length: %d, perc: %d, x_name: %.*s\n", 
    hap->snp, hap->length, (hap->snp == 0? 0: hap->length/hap->snp), 
    Get_NAME_LENGTH((*R_INF), overlap_list->list[0].x_id), Get_NAME((*R_INF),overlap_list->list[0].x_id)
    );
    fprintf(stderr, "hap->available_snp: %d, hap->core_snp:%d\n", 
    hap->available_snp, hap->core_snp);

    int Len_x, matched_overlap = 0;
        for (j = 0; j < overlap_list->length; j++)
        {
            Len_x = overlap_list->list[j].x_pos_e -  overlap_list->list[j].x_pos_s + 1;

            if (overlap_list->list[j].is_match == 1)
            {
                matched_overlap++;
            }
        }

        fprintf(stderr, "occ_0: %d, occ_1: %d, overlap_num: %d, matched_overlap: %d, overlap_list->length: %d\n", 
    hap->result_stat.occ_0, hap->result_stat.occ_1,
    hap->result_stat.overlap_num, matched_overlap, overlap_list->length); 
    

    return;

        fprintf(stderr, "Phaseing sucessfully!\n");
        fprintf(stderr, "occ_0: %d, occ_1: %d, overlap_num: %d\n", 
    hap->result_stat.occ_0, hap->result_stat.occ_1,
    hap->result_stat.overlap_num); 


        int8_t* vector = Get_Result_SNP_Vector((*hap));
        for (j = 0; j < hap->overlap; j++)
        {
            if(vector[j] == 0)
            {
                fprintf(stderr, "Ptype: %d, ID: %d\n", vector[j], j);
            }
        }


        for (j = 0; j < hap->overlap; j++)
        {
            if(vector[j] == 1)
            {
                fprintf(stderr, "Ptype: %d, ID: %d\n", vector[j], j);
            }
        }


        for (j = 0; j < hap->overlap; j++)
        {
            if(vector[j] == 2)
            {
                fprintf(stderr, "Ptype: %d, ID: %d\n", vector[j], j);
            }
        }


    for (i = 0; i < hap->core_snp; i++)
    {
        fprintf(stderr, "core(i): %d, site: %d, occ_0: %d, occ_1: %d, occ_2: %d, score: %d\n", 
    i, hap->snp_stat[i].site, hap->snp_stat[i].occ_0, hap->snp_stat[i].occ_1,
    hap->snp_stat[i].occ_2,
    hap->snp_stat[i].score); 

        int vectorID = hap->snp_stat[i].id;
        int8_t* vector = Get_SNP_Vector((*hap), vectorID);

        for (j = 0; j < hap->overlap; j++)
        {
            if(vector[j] == 0)
            {
                fprintf(stderr, "type: %d, ID: %d\n", vector[j], j);
            }
        }


        for (j = 0; j < hap->overlap; j++)
        {
            if(vector[j] == 1)
            {
                fprintf(stderr, "type: %d, ID: %d\n", vector[j], j);
            }
        }


        for (j = 0; j < hap->overlap; j++)
        {
            if(vector[j] == 2)
            {
                fprintf(stderr, "type: %d, ID: %d\n", vector[j], j);
            }
        }
        
    }



        for (; i < hap->available_snp; i++)
    {
        fprintf(stderr, "i: %d, site: %d, occ_0: %d, occ_1: %d, occ_2: %d, score: %d\n", 
    i, hap->snp_stat[i].site, hap->snp_stat[i].occ_0, hap->snp_stat[i].occ_1,
    hap->snp_stat[i].occ_2,
    hap->snp_stat[i].score); 

        int vectorID = hap->snp_stat[i].id;
        int8_t* vector = Get_SNP_Vector((*hap), vectorID);

        for (j = 0; j < hap->overlap; j++)
        {
            if(vector[j] == 0)
            {
                fprintf(stderr, "type: %d, ID: %d\n", vector[j], j);
            }
        }


        for (j = 0; j < hap->overlap; j++)
        {
            if(vector[j] == 1)
            {
                fprintf(stderr, "type: %d, ID: %d\n", vector[j], j);
            }
        }


        for (j = 0; j < hap->overlap; j++)
        {
            if(vector[j] == 2)
            {
                fprintf(stderr, "type: %d, ID: %d\n", vector[j], j);
            }
        }
        
    }
    
}



void debug_near_snp(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, haplotype_evdience_alloc* hap)
{
    int i, j, overlap_length, window_start, window_end, flag;
    long long window_num = (g_read->length + WINDOW - 1) / WINDOW;
    int total_read = 0;
    int unuseful_read = 0;
    for (i = 0; i < overlap_list->length; i++)
    {
        int flag = -1;
        overlap_length = overlap_list->list[i].x_pos_e - overlap_list->list[i].x_pos_s + 1;
        if (overlap_list->list[i].is_match == 1)
        {
            total_read++;
            for (j = 0; j < hap->available_snp; j++)
            {
                int vectorID = hap->snp_stat[j].id;
                int8_t* vector = Get_SNP_Vector((*hap), vectorID);

                if(vector[i] != 2 && flag == 2)
                {
                    unuseful_read++;
                    break;
                }

                if(vector[i] == 2)
                {
                    flag = 2;
                }
            }
        }
    }





    for (j = 1; j < hap->available_snp; j++)
    {
        if(hap->snp_stat[j].site <= hap->snp_stat[j - 1].site)
        {
            fprintf(stderr, "error\n");
        }
    }
    
    for (j = 0; j < hap->available_snp; j++)
    {
        fprintf(stderr, "\n\nsite: %d, score: %d\n", hap->snp_stat[j].site, hap->snp_stat[j].score );

        // if((j>0 && hap->snp_stat[j].site == hap->snp_stat[j - 1].site + 1)
        //     ||
        //    (j < hap->available_snp - 1 && hap->snp_stat[j].site + 1 == hap->snp_stat[j + 1].site))
        {

            fprintf(stderr, "x_name: %.*s\n", 
             Get_NAME_LENGTH((*R_INF), overlap_list->list[0].x_id), Get_NAME((*R_INF),overlap_list->list[0].x_id));

            int vectorID = hap->snp_stat[j].id;
            int8_t* vector = Get_SNP_Vector((*hap), vectorID);
            for (i = 0; i < hap->overlap; i++)
            {
                if(vector[i] == 0)
                {
                    fprintf(stderr, "type: %d, ID: %d\n", vector[i], i);
                }
            }

            for (i = 0; i < hap->overlap; i++)
            {
                if(vector[i] == 1)
                {
                    fprintf(stderr, "type: %d, ID: %d, %.*s\n", 
                    vector[i], i,
                    Get_NAME_LENGTH((*R_INF), overlap_list->list[i].y_id), Get_NAME((*R_INF),overlap_list->list[i].y_id));

                }
            }

            window_start = 0;
            window_end = WINDOW - 1;
            if (window_end >= g_read->length)
            {
                window_end = g_read->length - 1;
            }
            long long ijk;
            for (ijk = 0; ijk < window_num; ijk++)
            {


                if(hap->snp_stat[j].site <= window_end && hap->snp_stat[j].site >= window_start)
                {
                    fprintf(stderr, "winID: %d, winNum: %d, window_start: %d, window_end: %d\n",
                    ijk, window_num, window_start, window_end);
                    
                    dumy->length = 0;
                    dumy->lengthNT = 0;
                    ///flag返回的是重叠数量
                    ///dumy->length返回的是有效完全重叠的数量
                    ///dumy->lengthNT返回的是有效不完全重叠的数量
                    ///return overlaps that is overlaped with [window_start, window_end]
                    flag = get_available_interval(window_start, window_end, overlap_list, dumy);
                    switch (flag)
                    {
                        case 1:    ///找到匹配
                            break;
                        case 0:    ///没找到匹配
                            break;
                        case -2: ///下一个window也不会存在匹配, 直接跳出
                            i = window_num;
                            break;
                    }


                    for (i = 0; i < dumy->length; i++)
                    {
                        ///这个是那个overlap的ID，而不是overlap里对应窗口的ID
                        int overlapID = dumy->overlapID[i];

                        ///overlap_list->list[overlapID].x_pos_s is the begining of the whole overlap
                        int correct_x_pos_s = (overlap_list->list[overlapID].x_pos_s / WINDOW) * WINDOW;
                        ///window_start is the begining of this window in the whole x_read
                        int windowID = (window_start - correct_x_pos_s) / WINDOW;

                        ///如果这个window不匹配，跳过
                        if (overlap_list->list[overlapID].w_list[windowID].y_end == -1)
                        {
                            continue;
                        }

                        if(vector[overlapID] == 1)
                        {
                            fprintf(stderr, "overlapID: %d, x_strat: %d, x_end: %d, y_start: %d, y_end: %d\n",
                            overlapID,
                            overlap_list->list[overlapID].w_list[windowID].x_start,
                            overlap_list->list[overlapID].w_list[windowID].x_end,
                            overlap_list->list[overlapID].w_list[windowID].y_start,
                            overlap_list->list[overlapID].w_list[windowID].y_end);

                            recover_UC_Read_sub_region(dumy->overlap_region, overlap_list->list[overlapID].w_list[windowID].y_start, 
                            overlap_list->list[overlapID].w_list[windowID].y_end -overlap_list->list[overlapID].w_list[windowID].y_start + 1,
                            overlap_list->list[overlapID].y_pos_strand, 
                            R_INF, overlap_list->list[overlapID].y_id);

                            char* x_string = g_read->seq + overlap_list->list[overlapID].w_list[windowID].x_start;
                            char* y_string = dumy->overlap_region;

                            fprintf(stderr, "x_string: \n%.*s\n",  
                            overlap_list->list[overlapID].w_list[windowID].x_end - overlap_list->list[overlapID].w_list[windowID].x_start + 1,
                            x_string);

                            fprintf(stderr, "y_string: \n%.*s\n",  
                            overlap_list->list[overlapID].w_list[windowID].y_end - overlap_list->list[overlapID].w_list[windowID].y_start + 1,
                            y_string);


                            int haha_i = 0;

                            for (haha_i = 0; haha_i < overlap_list->list[overlapID].w_list[windowID].cigar.length; haha_i++)
                            {
                                fprintf(stderr, "oper: %d, len: %d\n", 
                                overlap_list->list[overlapID].w_list[windowID].cigar.C_C[haha_i],
                                overlap_list->list[overlapID].w_list[windowID].cigar.C_L[haha_i]
                                );
                            }

                        }

                        
                    }
                    
                }


                


                window_start = window_start + WINDOW;
                window_end = window_end + WINDOW;
                if (window_end >= g_read->length)
                {
                    window_end = g_read->length - 1;
                }
            }

            


        }

    }
    fprintf(stderr, "total_read: %d, unuseful_read: %d\n", total_read, unuseful_read);
    
}



void partition_overlaps(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, haplotype_evdience_alloc* hap,
                        int force_repeat)
{
    ResizeInitHaplotypeEvdience(hap);

    long long i, j, overlap_length;
    long long window_start, window_end;

    long long num_availiable_win = 0;
    
    Window_Pool w_inf;
    init_Window_Pool(&w_inf, g_read->length, WINDOW, TAIL_LENGTH);

    int flag = 0;
    while(get_Window(&w_inf, &window_start, &window_end) && flag != -2)
    {
        dumy->length = 0;
        dumy->lengthNT = 0;
        ///flag返回的是重叠数量
        ///dumy->length返回的是有效完全重叠的数量
        ///dumy->lengthNT返回的是有效不完全重叠的数量
        ///return overlaps that is overlaped with [window_start, window_end]
        flag = get_available_interval(window_start, window_end, overlap_list, dumy);
        switch (flag)
        {
            case 1:    ///找到匹配
                break;
            case 0:    ///没找到匹配
                break;
            case -2: ///下一个window也不会存在匹配, 直接跳出
                break;
        }
        
        ///这个是available overlap里所有window的数量...
        ///num_availiable_win = num_availiable_win + dumy->length + dumy->lengthNT;
        num_availiable_win = num_availiable_win + dumy->length;

        cluster(g_read->seq, window_start, window_end, overlap_list, dumy, R_INF, hap);
    }


    
    ///very time-consuming
    qsort(hap->list, hap->length, sizeof(haplotype_evdience), cmp_haplotype_evdience);

    
    

    ///debug_hap_information(overlap_list, R_INF, g_read, hap, dumy);

    SetSnpMatrix(hap, hap->snp, overlap_list->length);


    uint64_t pre_site = (uint64_t)-1;
    uint64_t num_of_snps = 0;
    long long pre_i = -1;
    long long sub_length;
    haplotype_evdience* sub_list;

    ///long long debug_total_length = 0;

    ////split reads
    for (i = 0; i < hap->length; i++)
    {
        if(pre_site != hap->list[i].site)
        {
            if(i != 0)
            {
                sub_list = hap->list + pre_i;
                sub_length = i - pre_i;
                ///debug_total_length = debug_total_length + sub_length;
                split_sub_list(hap, sub_list, sub_length, hap->snp, overlap_list, R_INF, g_read);
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
        ///debug_total_length = debug_total_length + sub_length;
        split_sub_list(hap, sub_list, sub_length, hap->snp, overlap_list, R_INF, g_read);
    }

    ///debug_snp_matrix(hap);
    generate_haplotypes_DP(hap, overlap_list, R_INF, g_read->length, force_repeat);
    ///generate_haplotypes_naive(hap, overlap_list, R_INF, g_read->length, force_repeat);

    lable_large_indels(overlap_list, R_INF, g_read->length, dumy);


    ///debug_snp_matrix(hap);
}


void partition_overlaps_advance(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, UC_Read* overlap_read, Correct_dumy* dumy, 
                        haplotype_evdience_alloc* hap, int force_repeat)
{
    ResizeInitHaplotypeEvdience(hap);

    long long i, j, overlap_length;
    long long window_start, window_end;

    long long num_availiable_win = 0;
    
    Window_Pool w_inf;
    init_Window_Pool(&w_inf, g_read->length, WINDOW, TAIL_LENGTH);

    int flag = 0;
    while(get_Window(&w_inf, &window_start, &window_end) && flag != -2)
    {
        dumy->length = 0;
        dumy->lengthNT = 0;
        ///flag返回的是重叠数量
        ///dumy->length返回的是有效完全重叠的数量
        ///dumy->lengthNT返回的是有效不完全重叠的数量
        ///return overlaps that is overlaped with [window_start, window_end]
        flag = get_available_interval(window_start, window_end, overlap_list, dumy);
        switch (flag)
        {
            case 1:    ///找到匹配
                break;
            case 0:    ///没找到匹配
                break;
            case -2: ///下一个window也不会存在匹配, 直接跳出
                break;
        }
        
        ///这个是available overlap里所有window的数量...
        num_availiable_win = num_availiable_win + dumy->length;
        
        ///need to deal with 
        cluster_advance(g_read->seq, window_start, window_end, overlap_list, 
        dumy, R_INF, hap, overlap_read);
    }


    
    ///very time-consuming
    qsort(hap->list, hap->length, sizeof(haplotype_evdience), cmp_haplotype_evdience);

    
    

    ///debug_hap_information(overlap_list, R_INF, g_read, hap, dumy);

    SetSnpMatrix(hap, hap->snp, overlap_list->length);


    uint64_t pre_site = (uint64_t)-1;
    uint64_t num_of_snps = 0;
    long long pre_i = -1;
    long long sub_length;
    haplotype_evdience* sub_list;

    ///long long debug_total_length = 0;

    ////split reads
    for (i = 0; i < hap->length; i++)
    {
        // if(overlap_list->list[0].x_id == 109837)
        // {
        //     fprintf(stderr, "hap->list[%d].site: %d, type: %d, misBase: %c, overlapID: %d, %.*s\n", 
        //     i, hap->list[i].site, hap->list[i].type, hap->list[i].misBase, hap->list[i].overlapID,
        //     Get_NAME_LENGTH((*R_INF), overlap_list->list[hap->list[i].overlapID].y_id),
        //     Get_NAME((*R_INF), overlap_list->list[hap->list[i].overlapID].y_id));
        // }
        if(pre_site != hap->list[i].site)
        {
            if(i != 0)
            {
                sub_list = hap->list + pre_i;
                sub_length = i - pre_i;
                ///debug_total_length = debug_total_length + sub_length;
                split_sub_list(hap, sub_list, sub_length, hap->snp, overlap_list, R_INF, g_read);
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
        ///debug_total_length = debug_total_length + sub_length;
        split_sub_list(hap, sub_list, sub_length, hap->snp, overlap_list, R_INF, g_read);
    }

    ///debug_snp_matrix(hap);
    generate_haplotypes_DP(hap, overlap_list, R_INF, g_read->length, force_repeat);
    ///generate_haplotypes_naive(hap, overlap_list, R_INF, g_read->length, force_repeat);

    lable_large_indels(overlap_list, R_INF, g_read->length, dumy);


    ///debug_snp_matrix(hap);
}




















void correct_overlap_back(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read, Graph* g, Graph* DAGCon,
                        long long* matched_overlap_0, long long* matched_overlap_1, 
                        long long* potiental_matched_overlap_0, long long* potiental_matched_overlap_1,
                        Cigar_record* current_cigar, haplotype_evdience_alloc* hap, 
                        Round2_alignment* second_round)
{

    
    reverse_complement(g_read->seq, g_read->length);

    clear_Correct_dumy(dumy, overlap_list);

    long long window_num = (g_read->length + WINDOW - 1) / WINDOW;

    long long i;

    long long window_start, window_end;

    window_start = 0;
    window_end = WINDOW - 1;
    if (window_end >= g_read->length)
    {
        window_end = g_read->length - 1;
    }



    int flag;

    for (i = 0; i < window_num; i++)
    {

        dumy->length = 0;
        dumy->lengthNT = 0;
        flag = get_interval(window_start, window_end, overlap_list, dumy);

        switch (flag)
        {
            case 1:    ///找到匹配
                break;
            case 0:    ///没找到匹配
                break;
            case -2: ///下一个window也不会存在匹配, 直接跳出
                i = window_num;
                break;
        }

        if(dumy->length + dumy->lengthNT>overlap_list->length)
        {
            fprintf(stderr, "error length\n");
        }
        
        
        ///verify_get_interval(window_start, window_end, overlap_list, dumy);

        verify_window(window_start, window_end, overlap_list, dumy, R_INF, g_read->seq);

        window_start = window_start + WINDOW;
        window_end = window_end + WINDOW;
        if (window_end >= g_read->length)
        {
            window_end = g_read->length - 1;
        }
    }

    debug_stats(overlap_list, R_INF, g_read, dumy, overlap_read, potiental_matched_overlap_0, potiental_matched_overlap_1);
    
    recalcate_window(overlap_list, R_INF, g_read, dumy, overlap_read);
    
    debug_stats(overlap_list, R_INF, g_read, dumy, overlap_read, matched_overlap_0, matched_overlap_1);
    
    partition_overlaps(overlap_list, R_INF, g_read, dumy, hap, 0);

    generate_consensus(overlap_list, R_INF, g_read, dumy, g, DAGCon, current_cigar, second_round);
}


void print_overlap(char* name, long long readID, 
overlap_region_alloc* overlap_list, All_reads* R_INF, int output_reads, int output_cigar)
{
    if(memcmp(name, Get_NAME((*R_INF), readID), 
    Get_NAME_LENGTH((*R_INF),readID)) == 0)
    {
        long long i, j;
        fprintf(stderr, "\n\n****************ref_read: %.*s, id: %d****************\n", 
            Get_NAME_LENGTH((*R_INF),readID), Get_NAME((*R_INF),readID), readID);

        fprintf(stderr, "\n###flag: 1\n");

        for (i = 0; i < overlap_list->length; i++)
        {
            if(overlap_list->list[i].is_match == 1)
            {
                fprintf(stderr, "%.*s\n", Get_NAME_LENGTH((*R_INF),overlap_list->list[i].y_id),
            Get_NAME((*R_INF),overlap_list->list[i].y_id));

            fprintf(stderr, "alignLen: %d, x_s: %d, x_e: %d, y_s: %d, y_e: %d, y_dir: %d, strong: %d\n", 
                overlap_list->list[i].align_length,
                overlap_list->list[i].x_pos_s,
                overlap_list->list[i].x_pos_e,
                overlap_list->list[i].y_pos_s,
                overlap_list->list[i].y_pos_e,
                overlap_list->list[i].y_pos_strand,
                overlap_list->list[i].strong);
            }
        }

        fprintf(stderr, "\n###flag: 2\n");

        for (i = 0; i < overlap_list->length; i++)
        {
            if(overlap_list->list[i].is_match == 2)
            {
                fprintf(stderr, "%.*s\n", Get_NAME_LENGTH((*R_INF),overlap_list->list[i].y_id),
            Get_NAME((*R_INF),overlap_list->list[i].y_id));

            fprintf(stderr, "alignLen: %d, x_s: %d, x_e: %d, y_s: %d, y_e: %d, y_dir: %d, strong: %d\n", 
                overlap_list->list[i].align_length,
                overlap_list->list[i].x_pos_s,
                overlap_list->list[i].x_pos_e,
                overlap_list->list[i].y_pos_s,
                overlap_list->list[i].y_pos_e,
                overlap_list->list[i].y_pos_strand,
                overlap_list->list[i].strong);
            }
        }

        fprintf(stderr, "\n###flag: 4\n");

        for (i = 0; i < overlap_list->length; i++)
        {
            if(overlap_list->list[i].is_match == 4)
            {
                fprintf(stderr, "%.*s\n", Get_NAME_LENGTH((*R_INF),overlap_list->list[i].y_id),
                Get_NAME((*R_INF),overlap_list->list[i].y_id));
                fprintf(stderr, "alignLen: %d, x_s: %d, x_e: %d, y_s: %d, y_e: %d, y_dir: %d, strong: %d\n", 
                overlap_list->list[i].align_length,
                overlap_list->list[i].x_pos_s,
                overlap_list->list[i].x_pos_e,
                overlap_list->list[i].y_pos_s,
                overlap_list->list[i].y_pos_e,
                overlap_list->list[i].y_pos_strand,
                overlap_list->list[i].strong);
            }
        }

        if(output_reads)
        {
            UC_Read g_read;
            init_UC_Read(&g_read);
            recover_UC_Read(&g_read, R_INF, readID);

            fprintf(stderr, "\n\nOutput all related reads\n");
            fprintf(stderr, "ref_read:\n");
            fprintf(stderr, ">%.*s\n", Get_NAME_LENGTH((*R_INF),readID), Get_NAME((*R_INF),readID));
            fprintf(stderr, "%.*s\n", g_read.length, g_read.seq);

            ///fprintf(stderr, "query_read:\n");
            for (i = 0; i < overlap_list->length; i++)
            {
                ///fprintf(stderr, "i: %d\n", i);
                recover_UC_Read(&g_read, R_INF, overlap_list->list[i].y_id);
                fprintf(stderr, ">%.*s\n", 
                Get_NAME_LENGTH((*R_INF),overlap_list->list[i].y_id),
                Get_NAME((*R_INF),overlap_list->list[i].y_id));
                fprintf(stderr, "%.*s\n", g_read.length, g_read.seq);
            }

            fprintf(stderr, "Has already output all related reads\n\n");

            destory_UC_Read(&g_read);
         }


        if(output_cigar)
        {
            for (i = 0; i < overlap_list->length; i++)
            {
                
                fprintf(stderr, "\ni: %d, %.*s, x_s: %d, x_e: %d, y_s: %d, y_end: %d, w_list_length: %d, dir: %d, strong: %d, is_match: %d\n", 
                i, Get_NAME_LENGTH((*R_INF),overlap_list->list[i].y_id),
                Get_NAME((*R_INF),overlap_list->list[i].y_id),
                overlap_list->list[i].x_pos_s,
                overlap_list->list[i].x_pos_e,
                overlap_list->list[i].y_pos_s,
                overlap_list->list[i].y_pos_e,
                overlap_list->list[i].w_list_length,
                overlap_list->list[i].y_pos_strand,
                overlap_list->list[i].strong,
                overlap_list->list[i].is_match);

                fprintf(stderr, "window cigar: \n");                
                for (j = 0; j < overlap_list->list[i].w_list_length; j++)
                {
                    fprintf(stderr, "************************\ncigar_j: %d, x_s: %d, x_e: %d, y_s: %d, y_end: %d\n", 
                    j, overlap_list->list[i].w_list[j].x_start,
                    overlap_list->list[i].w_list[j].x_end, 
                    overlap_list->list[i].w_list[j].y_start,
                    overlap_list->list[i].w_list[j].y_end);
                    if(overlap_list->list[i].w_list[j].y_end == -1)
                    {
                        fprintf(stderr, "not match\n");
                    }
                    else
                    {
                        int cigar_i, operation, operationLen;
                        CIGAR* cigar = &(overlap_list->list[i].w_list[j].cigar);
                        fprintf(stderr, "length: %d\n", cigar->length);
                        for (cigar_i = 0; cigar_i < cigar->length; cigar_i++)
                        {
                            operation = cigar->C_C[cigar_i];
                            operationLen = cigar->C_L[cigar_i];
                            fprintf(stderr, "oper: %d, Len: %d\n", operation, operationLen);
                        }
                    }
                }



                fprintf(stderr, "boundary cigar: \n");
                for (j = 0; j < overlap_list->list[i].boundary_cigars.length; j++)
                {
                    fprintf(stderr, "###################\ncigar_j: %d, x_s: %d, x_e: %d, y_s: %d, y_end: %d\n", 
                    j, overlap_list->list[i].boundary_cigars.buffer[j].x_start,
                    overlap_list->list[i].boundary_cigars.buffer[j].x_end,
                    overlap_list->list[i].boundary_cigars.buffer[j].y_start,
                    overlap_list->list[i].boundary_cigars.buffer[j].y_end);


                    if(overlap_list->list[i].boundary_cigars.buffer[j].y_end == -1)
                    {
                        fprintf(stderr, "not match\n");
                    }
                    else
                    {
                        int cigar_i, operation, operationLen;
                        CIGAR* cigar = &(overlap_list->list[i].boundary_cigars.buffer[j].cigar);
                        fprintf(stderr, "length: %d\n", cigar->length);
                        for (cigar_i = 0; cigar_i < cigar->length; cigar_i++)
                        {
                            operation = cigar->C_C[cigar_i];
                            operationLen = cigar->C_L[cigar_i];
                            fprintf(stderr, "oper: %d, Len: %d\n", operation, operationLen);
                        }
                    }
                }   
            }
        }


        fprintf(stderr, "End\n\n\n\n\n");
       
    }

}



void correct_overlap(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read, Graph* g, Graph* DAGCon,
                        long long* matched_overlap_0, long long* matched_overlap_1, 
                        long long* potiental_matched_overlap_0, long long* potiental_matched_overlap_1,
                        Cigar_record* current_cigar, haplotype_evdience_alloc* hap, 
                        Round2_alignment* second_round, int force_repeat, int is_consensus,
                        int* fully_cov, int* abnormal, uint8_t* c2n)
{
    reverse_complement(g_read->seq, g_read->length);

    clear_Correct_dumy(dumy, overlap_list);

    long long window_start, window_end;

    Window_Pool w_inf;

    init_Window_Pool(&w_inf, g_read->length, WINDOW, TAIL_LENGTH);

    int flag = 0;

    while(get_Window(&w_inf, &window_start, &window_end) && flag != -2)
    {
        dumy->length = 0;
        dumy->lengthNT = 0;
        flag = get_interval(window_start, window_end, overlap_list, dumy);

        switch (flag)
        {
            case 1:    ///no match here
                break;
            case 0:    ///no match here
                break;
            case -2: ///if flag == -2, loop would be terminated
                break;
        }

        if(dumy->length + dumy->lengthNT>overlap_list->length)
        {
            fprintf(stderr, "error length\n");
        }

        ///dumy->lengthNT represent how many overlaps that the length of them is not equal to WINDOW; may larger or less than WINDOW
        ///dumy->length represent how many overlaps that the length of them is WINDOW
        /****************************may improve**************************/
        ///now the windows which are larger than WINDOW are verified one-by-one, to improve it, we can do it group-bygroup
        verify_window(window_start, window_end, overlap_list, dumy, R_INF, g_read->seq);
    }

    /**
    recalcate_window(overlap_list, R_INF, g_read, dumy, overlap_read);
    partition_overlaps(overlap_list, R_INF, g_read, dumy, hap, force_repeat);
    **/

    recalcate_window_advance(overlap_list, R_INF, g_read, dumy, overlap_read);
    ///partition_overlaps(overlap_list, R_INF, g_read, dumy, hap, force_repeat);
    partition_overlaps_advance(overlap_list, R_INF, g_read, overlap_read, dumy, hap, force_repeat);






    // print_overlap("m64013_190410_223304/61081096/ccs", 
    // overlap_list->list[0].x_id, overlap_list, R_INF, 0, 1);

    // print_overlap("m64013_190412_043951/83756835/ccs", 
    // overlap_list->list[0].x_id, overlap_list, R_INF, 0, 1);

    // print_overlap("m64013_190322_203854/152832751/ccs", 
    // overlap_list->list[0].x_id, overlap_list, R_INF, 0, 1);

    // print_overlap("m64011_190329_072846/175047040/ccs", 
    // overlap_list->list[0].x_id, overlap_list, R_INF, 0, 1);

    // print_overlap("m64013_190322_203854/13174300/ccs", 
    // overlap_list->list[0].x_id, overlap_list, R_INF, 0, 1);

    // print_overlap("m64013_190322_203854/177145456/ccs", 
    // overlap_list->list[0].x_id, overlap_list, R_INF, 0, 1);

    // print_overlap("m64013_190322_203854/120717695/ccs", 
    // overlap_list->list[0].x_id, overlap_list, R_INF, 0, 1);

    // print_overlap("m64011_190329_072846/76548728/ccs", 
    // overlap_list->list[0].x_id, overlap_list, R_INF, 0, 1);

    // print_overlap("m64013_190324_024932/23660629/ccs", 
    // overlap_list->list[0].x_id, overlap_list, R_INF, 0, 1);

    


    if(is_consensus)
    {
        generate_consensus(overlap_list, R_INF, g_read, dumy, g, DAGCon, current_cigar, second_round);
    }


    (*fully_cov) = check_if_fully_covered(overlap_list, R_INF, g_read, dumy, g, abnormal);
    
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

void clear_Correct_dumy(Correct_dumy* list, overlap_region_alloc* overlap_list)
{
    list->length = 0;
    list->lengthNT = 0;
    list->start_i = 0;

    if (list->size < overlap_list->length)
    {
        list->size = overlap_list->length;
        list->overlapID = (uint64_t*)realloc(list->overlapID, list->size*sizeof(uint64_t));
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


void pre_filter_by_nearby_single(k_mer_pos* new_n_list, k_mer_pos* old_n_list, uint64_t n_length, uint64_t n_end_pos, UC_Read* g_read, 
All_reads* R_INF, Correct_dumy* dumy, uint64_t* new_n_length)
{
    (*new_n_length) = 0;
    ///这种就是0，new_n_list也不需要有数据
    if (n_length == 0)
    {
        return;
    }
    
    char* x_string = NULL;
    char* y_string = NULL;
    

    long long  x_offset = (long long)(n_end_pos) - WINDOW + 1;
    if (x_offset < 0)
    {
        x_offset = 0;
    }

    ///x_length是x上待验证区间的实际长度
    ///x_offset是x上待验证区间的起始位置
    ///如果是向前取待验证区间，那么x_length至少大于等于40 (compressed k-mer长度，不是实际长度)
    long long x_length = n_end_pos - x_offset + 1;
    ///这种也不要过滤了，直接把old_n_list全部赋值过来就好了
    if (x_length < (WINDOW/2))
    {
        (*new_n_length) = n_length;
        memcpy(new_n_list, old_n_list, sizeof(k_mer_pos)*n_length);
        return;
    }
    long long Window_Len = x_length + (THRESHOLD << 1);
    x_string = g_read->seq + x_offset;


    long long y_offset;
    long long y_length;

    long long i = 0;

    long long y_read_length;
    int end_site;
    unsigned int error;

    for (i = 0; i < n_length; i++)
    {
        ////old_n_list[i].offset是y上k-mer的结束位置
        ///n_end_pos是x上k-mer的结束位置
        ///x_length是x上区间长度
        y_offset = (long long)(old_n_list[i].offset) - x_length;
        ///这种情况下弃疗
        if (y_offset < 0)
        {
            new_n_list[(*new_n_length)].readID = old_n_list[i].readID;
            new_n_list[(*new_n_length)].offset = old_n_list[i].offset;
            (*new_n_length)++;
            continue;
        }

        y_offset = y_offset - THRESHOLD;
        ///还能抢救
        if (y_offset < 0)
        {
            y_offset = 0;
        }

        ///y的read的总长度
        y_read_length = Get_READ_LENGTH((*R_INF), old_n_list[i].readID);
        y_length = MIN(Window_Len, y_read_length - y_offset);

        ///如果y的长度比x长度还小，那就直接弃疗了
        if (y_length < x_length)
        {
            new_n_list[(*new_n_length)].readID = old_n_list[i].readID;
            new_n_list[(*new_n_length)].offset = old_n_list[i].offset;
            (*new_n_length)++;
            continue;
        }
        
        ///y的方向都是0，因为索引里都是0
        recover_UC_Read_sub_region(dumy->overlap_region, y_offset, y_length, 0, R_INF, old_n_list[i].readID);
        y_string = dumy->overlap_region;

        memset (y_string + y_length, 0, Window_Len - y_length);
        end_site = Reserve_Banded_BPM(y_string, y_length, x_string, x_length, THRESHOLD, &error);

        if (error!=(unsigned int)-1)
        {
            new_n_list[(*new_n_length)].readID = old_n_list[i].readID;
            new_n_list[(*new_n_length)].offset = old_n_list[i].offset;
            (*new_n_length)++;
        }  
    }
}



void pre_filter_by_nearby(k_mer_pos* new_n_list, k_mer_pos* old_n_list, uint64_t n_length, uint64_t n_end_pos, UC_Read* g_read, 
All_reads* R_INF, Correct_dumy* dumy, uint64_t* new_n_length)
{
    (*new_n_length) = 0;
    ///这种就是0，new_n_list也不需要有数据
    if (n_length == 0)
    {
        return;
    }
    
    char* x_string = NULL;
    char* y_string = NULL;
    

    long long  x_offset = (long long)(n_end_pos) - WINDOW + 1;
    if (x_offset < 0)
    {
        x_offset = 0;
    }

    ///x_length是x上待验证区间的实际长度
    ///x_offset是x上待验证区间的起始位置
    ///如果是向前取待验证区间，那么x_length至少大于等于40 (compressed k-mer长度，不是实际长度)
    long long x_length = n_end_pos - x_offset + 1;
    ///这种也不要过滤了，直接把old_n_list全部赋值过来就好了
    if (x_length < (WINDOW/2))
    {
        (*new_n_length) = n_length;
        memcpy(new_n_list, old_n_list, sizeof(k_mer_pos)*n_length);
        return;
    }
    long long Window_Len = x_length + (THRESHOLD << 1);
    x_string = g_read->seq + x_offset;


    long long y_offset;
    long long y_length;

    long long i = 0;

    long long y_read_length;
    int end_site;
    unsigned int error;

    int groupLen = 0;
    int return_sites[GROUP_SIZE];
    unsigned int return_sites_error[GROUP_SIZE];
    uint64_t readID[GROUP_SIZE];
    uint64_t offset[GROUP_SIZE];

    for (i = 0; i < n_length; i++)
    {
        ////old_n_list[i].offset是y上k-mer的结束位置
        ///n_end_pos是x上k-mer的结束位置
        ///x_length是x上区间长度
        y_offset = (long long)(old_n_list[i].offset) - x_length;
        ///这种情况下弃疗
        if (y_offset < 0)
        {
            new_n_list[(*new_n_length)].readID = old_n_list[i].readID;
            new_n_list[(*new_n_length)].offset = old_n_list[i].offset;
            (*new_n_length)++;
            continue;
        }

        y_offset = y_offset - THRESHOLD;
        ///还能抢救
        if (y_offset < 0)
        {
            y_offset = 0;
        }

        ///y的read的总长度
        y_read_length = Get_READ_LENGTH((*R_INF), old_n_list[i].readID);
        y_length = MIN(Window_Len, y_read_length - y_offset);

        ///如果y的长度比x长度还小，那就直接弃疗了
        if (y_length < x_length)
        {
            new_n_list[(*new_n_length)].readID = old_n_list[i].readID;
            new_n_list[(*new_n_length)].offset = old_n_list[i].offset;
            (*new_n_length)++;
            continue;
        }

        if(y_length == Window_Len)
        {
             ///y的方向都是0，因为索引里都是0
            recover_UC_Read_sub_region(dumy->overlap_region_group[groupLen], y_offset, y_length, 0, 
                            R_INF, old_n_list[i].readID);
            readID[groupLen] = old_n_list[i].readID;
            offset[groupLen] = old_n_list[i].offset;


            groupLen++;
            if (groupLen == GROUP_SIZE)
            {
                Reserve_Banded_BPM_4_SSE_only(dumy->overlap_region_group[0], dumy->overlap_region_group[1], 
                dumy->overlap_region_group[2], dumy->overlap_region_group[3], y_length, x_string, x_length,
                    return_sites, return_sites_error, THRESHOLD, dumy->Peq_SSE);
                groupLen = 0;

                if (return_sites_error[0]!=(unsigned int)-1)
                {
                    new_n_list[(*new_n_length)].readID = readID[0];
                    new_n_list[(*new_n_length)].offset = offset[0];
                    (*new_n_length)++;
                }

                if (return_sites_error[1]!=(unsigned int)-1)
                {
                    new_n_list[(*new_n_length)].readID = readID[1];
                    new_n_list[(*new_n_length)].offset = offset[1];
                    (*new_n_length)++;
                }

                if (return_sites_error[2]!=(unsigned int)-1)
                {
                    new_n_list[(*new_n_length)].readID = readID[2];
                    new_n_list[(*new_n_length)].offset = offset[2];
                    (*new_n_length)++;
                }

                if (return_sites_error[3]!=(unsigned int)-1)
                {
                    new_n_list[(*new_n_length)].readID = readID[3];
                    new_n_list[(*new_n_length)].offset = offset[3];
                    (*new_n_length)++;
                }
                
            }

        }
        else
        {
            ///y的方向都是0，因为索引里都是0
            recover_UC_Read_sub_region(dumy->overlap_region, y_offset, y_length, 0, R_INF, old_n_list[i].readID);
            y_string = dumy->overlap_region;

            memset (y_string + y_length, 0, Window_Len - y_length);
            end_site = Reserve_Banded_BPM(y_string, y_length, x_string, x_length, THRESHOLD, &error);

            if (error!=(unsigned int)-1)
            {
                new_n_list[(*new_n_length)].readID = old_n_list[i].readID;
                new_n_list[(*new_n_length)].offset = old_n_list[i].offset;
                (*new_n_length)++;
            }
        }
        
        
  
    }



    if (groupLen == 1)
    {
         end_site = Reserve_Banded_BPM(dumy->overlap_region_group[0], Window_Len, x_string, x_length, THRESHOLD, &error);

        if (error!=(unsigned int)-1)
        {
            new_n_list[(*new_n_length)].readID = readID[0];
            new_n_list[(*new_n_length)].offset = offset[0];
            (*new_n_length)++;
        }

    }
    else
    {
        Reserve_Banded_BPM_4_SSE_only(dumy->overlap_region_group[0], dumy->overlap_region_group[1], 
                dumy->overlap_region_group[2], dumy->overlap_region_group[3], Window_Len, x_string, x_length,
                    return_sites, return_sites_error, THRESHOLD, dumy->Peq_SSE);
        

        for (i = 0; i < groupLen; i++)
        {
            if (return_sites_error[i]!=(unsigned int)-1)
            {
                new_n_list[(*new_n_length)].readID = readID[i];
                new_n_list[(*new_n_length)].offset = offset[i];
                (*new_n_length)++;
            }
        }

        groupLen = 0;
        

    }

    ///要排序....
    if ((*new_n_length)>1)
    {
        qsort(new_n_list, (*new_n_length), sizeof(k_mer_pos), cmp_k_mer_pos);
    }
}



/**********************for prefilter************************ */

void destory_k_mer_pos_list_alloc_prefilter(k_mer_pos_list_alloc* list)
{
    long long i = 0;
    for (i = 0; i < list->size; i++)
    {
        if (list->list[i].size != 0)
        {
            free(list->list[i].list);
        }
    }

    free(list->list);
}


void append_k_mer_pos_list_alloc_prefilter(k_mer_pos_list_alloc* list, k_mer_pos* n_list, uint64_t n_length, 
uint64_t n_end_pos, uint8_t n_direction, UC_Read* g_read, All_reads* R_INF, Correct_dumy* dumy)
{
   
    if (list->length + 1 > list->size)
    {
        list->size = list->size * 2;
        list->list = (k_mer_pos_list*)realloc(list->list, sizeof(k_mer_pos_list)*list->size);
        ///新分配空间要初始化
        memset(list->list + (list->size/2), 0, sizeof(k_mer_pos_list)*(list->size/2));
    }

    if (list->list[list->length].size < n_length)
    {
        list->list[list->length].size = n_length;
        list->list[list->length].list = (k_mer_pos*)realloc(list->list[list->length].list, 
                    sizeof(k_mer_pos)*list->list[list->length].size);
    }
    
    
    ///list->list[list->length].list = n_list;
    ///memcpy(list->list[list->length].list, n_list, sizeof(k_mer_pos)*n_length);
    pre_filter_by_nearby(list->list[list->length].list, n_list, n_length, n_end_pos, 
            g_read, R_INF, dumy, &n_length);
    
            
    if(n_length > 0)
    {
        list->list[list->length].length = n_length;
        list->list[list->length].direction = n_direction;
        list->list[list->length].end_pos = n_end_pos;

        list->length++;
    }
}


/**********************for prefilter************************ */


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
