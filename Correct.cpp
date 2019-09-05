#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include "Correct.h"
#include "Levenshtein_distance.h"
#include "edlib.h"
#include "Assembly.h"

long long T_total_match=0;
long long T_total_unmatch=0;
long long T_total_mis=0;
pthread_mutex_t debug_statistics ;


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


    for (i = dumy->start_i; i < overlap_list->length; i++)
    {
        ///只会发生在这个interval比list里所有元素都小的情况
        ///这种情况下一个interval需要从0开始
        if (window_end < overlap_list->list[i].x_pos_s)
        {
            dumy->start_i = 0;
            dumy->length = 0;
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
        return -2;
    }
    
    dumy->length = 0;
    dumy->lengthNT = 0;

    for (; i < overlap_list->length; i++)
    {
        
        if((Len = OVERLAP(window_start, window_end, overlap_list->list[i].x_pos_s, overlap_list->list[i].x_pos_e)) > 0)
        {
            if (Len == WINDOW)
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
            if (overlap_length * OVERLAP_THRESHOLD <=  overlap_list->list[i].align_length)
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

void determine_overlap_region(int threshold, long long y_start, long long y_ID, long long Window_Len, All_reads* R_INF,
int* r_extra_begin, int* r_extra_end, long long* r_y_start, long long* r_y_length)
{
    int extra_begin;
    int extra_end;
    long long currentIDLen;
    long long o_len;

    extra_begin = extra_end = 0;
    ///y maybe less than 0
    y_start = y_start - threshold;
    ///the length of y
    currentIDLen = Get_READ_LENGTH((*R_INF), y_ID);
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


  
        determine_overlap_region(THRESHOLD, y_start, overlap_list->list[currentID].y_id, Window_Len, R_INF,
        &extra_begin, &extra_end, &y_start, &o_len);

        fill_subregion(dumy->overlap_region_group[groupLen], y_start, o_len, overlap_list->list[currentID].y_pos_strand, 
        R_INF, overlap_list->list[currentID].y_id, extra_begin, extra_end);

        y_extra_begin[groupLen] = extra_begin;
        y_extra_end[groupLen] = extra_end;
        overlapID[groupLen] = currentID;
        y_startGroup[groupLen] = y_start;
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
                            y_extra_begin[0], y_extra_end[0]);
            }
            else
            {
                append_window_list(&overlap_list->list[overlapID[0]], window_start, window_end, y_startGroup[0], -1, -1,
                y_extra_begin[0], y_extra_end[0]);
            }
            

            if (return_sites_error[1]!=(unsigned int)-1)
            {
                overlap_list->list[overlapID[1]].align_length += x_len;

                append_window_list(&overlap_list->list[overlapID[1]], window_start, window_end, 
                            y_startGroup[1], y_startGroup[1] + return_sites[1], (int)return_sites_error[1],
                            y_extra_begin[1], y_extra_end[1]);
            }
            else
            {
                append_window_list(&overlap_list->list[overlapID[1]], window_start, window_end, y_startGroup[1], -1, -1,
                y_extra_begin[1], y_extra_end[1]);
            }
            

            if (return_sites_error[2]!=(unsigned int)-1)
            {
                overlap_list->list[overlapID[2]].align_length += x_len;

                append_window_list(&overlap_list->list[overlapID[2]], window_start, window_end, 
                            y_startGroup[2], y_startGroup[2] + return_sites[2], (int)return_sites_error[2],
                            y_extra_begin[2], y_extra_end[2]);
            }
            else
            {
                append_window_list(&overlap_list->list[overlapID[2]], window_start, window_end, y_startGroup[2], -1, -1,
                y_extra_begin[2], y_extra_end[2]);
            }
            

            if (return_sites_error[3]!=(unsigned int)-1)
            {
                overlap_list->list[overlapID[3]].align_length += x_len;

                append_window_list(&overlap_list->list[overlapID[3]], window_start, window_end, 
                            y_startGroup[3], y_startGroup[3] + return_sites[3], (int)return_sites_error[3],
                            y_extra_begin[3], y_extra_end[3]);
            }
            else
            {
                append_window_list(&overlap_list->list[overlapID[3]], window_start, window_end, y_startGroup[3], -1, -1,
                y_extra_begin[3], y_extra_end[3]);
            }    
        }
    }


    

    if (groupLen == 1)
    {
        end_site = Reserve_Banded_BPM(dumy->overlap_region_group[0], Window_Len, x_string, WINDOW, THRESHOLD, &error);

        

        if (error!=(unsigned int)-1)
        {
            overlap_list->list[overlapID[0]].align_length += x_len;

            append_window_list(&overlap_list->list[overlapID[0]], window_start, window_end, 
                                y_startGroup[0], y_startGroup[0] + end_site, (int)error,
                                y_extra_begin[0], y_extra_end[0]);
        }
        else
        {
            append_window_list(&overlap_list->list[overlapID[0]], window_start, window_end, y_startGroup[0], -1, -1,
            y_extra_begin[0], y_extra_end[0]);
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
                                y_extra_begin[i], y_extra_end[i]);
            }
            else
            {
                append_window_list(&overlap_list->list[overlapID[i]], window_start, window_end, y_startGroup[i], -1, -1,
                y_extra_begin[i], y_extra_end[i]);
            }
            
        }

        groupLen = 0;
    }
    

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

        ///y上的相对位置
        y_start = (x_start - overlap_list->list[currentID].x_pos_s) + overlap_list->list[currentID].y_pos_s;

        
        Window_Len = x_len + (threshold << 1);
        determine_overlap_region(threshold, y_start, overlap_list->list[currentID].y_id, Window_Len, R_INF,
        &extra_begin, &extra_end, &y_start, &o_len);

        fill_subregion(dumy->overlap_region, y_start, o_len, overlap_list->list[currentID].y_pos_strand, 
        R_INF, overlap_list->list[currentID].y_id, extra_begin, extra_end);
        

        x_string = r_string + x_start;
        y_string = dumy->overlap_region;


        end_site = Reserve_Banded_BPM(y_string, Window_Len, x_string, x_len, threshold, &error);


        if (error!=(unsigned int)-1)
        {
            overlap_list->list[currentID].align_length += x_len;
            append_window_list(&overlap_list->list[currentID], x_start, x_end, y_start, y_start + end_site, (int)error,
            extra_begin, extra_end);
        }
        else
        {
            append_window_list(&overlap_list->list[currentID], x_start, x_end, y_start, -1, -1,
            extra_begin, extra_end);
        }
    }

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

        if (Len_x * OVERLAP_THRESHOLD <=  overlap_list->list[j].align_length)
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

            /**
            threshold = Len_x * 0.04;

            y_start = overlap_list->list[j].y_pos_s - threshold;
            if (y_start < 0)
            {
                y_start = 0;
            }

            
            ///当前y的长度
            currentIDLen = Get_READ_LENGTH((*R_INF), overlap_list->list[j].y_id);
            ///不能超过y的剩余长度
            Len_y = MIN(Len_x + 2 * threshold, currentIDLen - y_start);
            
            
            if (overlap_list->list[j].y_pos_strand == 0)
            {
                recover_UC_Read(overlap_read, R_INF, overlap_list->list[j].y_id);
                (*matched_overlap_0)++;
            }
            else
            {
                recover_UC_Read_RC(overlap_read, R_INF, overlap_list->list[j].y_id);
                (*matched_overlap_1)++;
            }
        

            
            
            EdlibAlignResult result = edlibAlign(g_read->seq + overlap_list->list[j].x_pos_s, Len_x, 
            overlap_read->seq + y_start, Len_y, 
            edlibNewAlignConfig(threshold, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
            

            if (result.status == EDLIB_STATUS_OK && result.editDistance != -1) {
                char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
                int cigar_length = strlen(cigar);
                free(cigar);

                
                if (overlap_list->list[j].y_pos_strand == 0)
                {
                    (*potiental_matched_overlap_0)++;
                }
                else
                {
                    (*potiental_matched_overlap_1)++;
                }
                
                // if (overlap_list->list[j].shared_seed == 1 && Len_x >= 1000)
                // {
                //     (*matched_overlap_0)++;
                // }

                // if (overlap_list->list[j].shared_seed == 2 && Len_x >= 1000)
                // {
                //     (*matched_overlap_1)++;
                // }
            
                
            }

            edlibFreeAlignResult(result);
            **/
        }
    }

}


void output_stats(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read)
{
    long long j;
    long long Len_x;
    int threshold;
    long long y_start;
    long long Len_y;
    long long currentIDLen;
    long long whole;
    for (j = 0; j < overlap_list->length; j++)
    {
        Len_x = overlap_list->list[j].x_pos_e -  overlap_list->list[j].x_pos_s + 1;

        if (Len_x * OVERLAP_THRESHOLD <=  overlap_list->list[j].align_length)
        {
            ;
        }
    }

}





inline void generate_cigar(
	char* path, int path_length, window_list* result, int* start, int* end, int error)
{

    if (error == 0)
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


    for (i = 0; i < path_length; i++)
    {
        if(path[i] == 1)
        {
            path[i] = 3;
            (*end)--;
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
        fprintf(stderr, "error cigar_error: cigar_error: %d, error: %d\n", cigar_error, error);
        for (i = 0; i < cigar->length; i++)
        {
            fprintf(stderr, "%u: %u\n", cigar->C_L[i], cigar->C_C[i]);
        }
        
    }
    
   
    if (flag_error == 1)
    {
        print_string(x, x_len);
        print_string(y, y_len);
        fprintf(stderr, "x_len: %d, y_len: %d, cigar_len: %d, error: %d\n", x_len, y_len, cigar->length, error);
        for (i = 0; i < cigar->length; i++)
        {
            fprintf(stderr, "%u: %u\n", cigar->C_L[i], cigar->C_C[i]);
        }
    }
    
    
    return flag_error;
    
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


/**
///记住等于0也要重算，虽然意义不那么大就是了
#define NEED_START_POS(x) (x.error<=0)
#define GET_ERROR(x) ((x.error>=0)? x.error : x.error*-1)
**/

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
                    threshold = x_len * THRESHOLD_RATE;


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
                    threshold = x_len * THRESHOLD_RATE;
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
                        


                         
                        generate_cigar(dumy->path, dumy->path_length, &(overlap_list->list[j].w_list[i]),
                        &real_y_start, &end_site, error);  

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
                    threshold = x_len * THRESHOLD_RATE;
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

                        generate_cigar(dumy->path, dumy->path_length, &(overlap_list->list[j].w_list[k]),
                        &real_y_start, &end_site, error);   

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




    
    ///j负责遍历整个overlap list
    for (j = 0; j < overlap_list->length; j++)
    {
        y_id = overlap_list->list[j].y_id;
        y_strand = overlap_list->list[j].y_pos_strand;
        y_readLen = Get_READ_LENGTH((*R_INF), y_id);
        overlap_length = overlap_list->list[j].x_pos_e - overlap_list->list[j].x_pos_s + 1;

        ///only calculate cigar for high quality overlaps
        if (overlap_length * OVERLAP_THRESHOLD <=  overlap_list->list[j].align_length)
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
                        threshold = x_len * THRESHOLD_RATE;
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

                            generate_cigar(dumy->path, dumy->path_length, &(overlap_list->list[j].w_list[i]),
                            &real_y_start, &end_site, error);    


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
        if (overlap_length * OVERLAP_THRESHOLD <=  overlap_list->list[j].align_length)
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
    }
    **/
    
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

void get_seq_from_Graph(Graph* backbone, Correct_dumy* dumy, Cigar_record* current_cigar, char* self_string)
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
                for (i = 0; i < backbone->g_nodes.list[currentNodeID].insertion_edges.length; i++)
                {
                    total_count = total_count + backbone->g_nodes.list[currentNodeID].insertion_edges.list[i].weight;

                    if (backbone->g_nodes.list[currentNodeID].insertion_edges.list[i].weight > max_count)
                    {
                        max_count = backbone->g_nodes.list[currentNodeID].insertion_edges.list[i].weight;
                        max_edge = i;
                        max_type = INSERTION;
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
            }
            
            ///这种情况下矫正
            if(max_count >= total_count*CORRECT_THRESHOLD)
            {
                currentNodeID = add_path_to_correct_read(backbone, dumy, currentNodeID, max_type, max_edge, current_cigar,
                self_string);
            }
            else  ///不矫正, 直接取下一个backbone节点
            {
                currentNodeID++;
                add_base_to_correct_read_directly(dumy, backbone->g_nodes.list[currentNodeID].base);

                add_cigar_record(&(backbone->g_nodes.list[currentNodeID].base), 1, current_cigar, 0);
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


void window_consensus(char* r_string, long long window_start, long long window_end, 
overlap_region_alloc* overlap_list, Correct_dumy* dumy, All_reads* R_INF, Graph* g, Cigar_record* current_cigar)
{
    clear_Graph(g);

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
        /**
        fprintf(stderr, "****window_start: %d, x_start: %d, x_end: %d, x_length: %d, y_start: %d, y_length: %d, dumy->last_boundary_length: %d\n", 
        window_start, x_start, overlap_list->list[overlapID].w_list[windowID].x_end , x_length, y_start, y_length,
        dumy->last_boundary_length);
        **/


        recover_UC_Read_sub_region(dumy->overlap_region, y_start, y_length, overlap_list->list[overlapID].y_pos_strand, 
                R_INF, overlap_list->list[overlapID].y_id);

        x_string = r_string + x_start;
        y_string = dumy->overlap_region;
        ///这个是比对上的起始base在backbone上对应的位置，也就是节点ID
        currentNodeID = x_start - window_start;
        ///这个是要用的cigar: overlap_list->list[overlapID].w_list[windowID].cigar;

        addmatchedSeqToGraph(g, currentNodeID, x_string, x_length, 
                    y_string, y_length, &(overlap_list->list[overlapID].w_list[windowID].cigar), startNodeID, endNodeID);
    }

    get_seq_from_Graph(g, dumy, current_cigar, backbone);


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
                merge_base = x_string[x_i];
                merge_base = merge_base << 3;
                merge_base = merge_base | y_string[y_i];
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
    /**
    fprintf(stderr, "start_cigar: %d, end_cigar: %d\n", start_cigar, end_cigar);
    for (int ijk = 0; ijk < new_cigar->length; ijk++)
    {
        fprintf(stderr, "new:Oper: %d, Len: %d\n", Get_Cigar_Type(new_cigar->record[ijk]), 
        Get_Cigar_Length(new_cigar->record[ijk]));
    }
    **/
    
    long long total_x_start = total_window_start + get_x_start;
    long long x_length = get_x_end -get_x_start + 1;
    long long total_y_start = get_y_start;
    long long y_length = get_y_end -get_y_start + 1;


    add_cigar_to_cigar(current_dumy, current_cigar, second_round,
    total_x_start, x_length, total_y_start, y_length);

}

int process_boundary(overlap_region_alloc* overlap_list, All_reads* R_INF, Correct_dumy* dumy, Graph* g, 
Cigar_record* current_cigar, long long uncorrected_window_start, Round2_alignment* second_round)
{
    char* r_string = dumy->corrected_read;
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

        ///如果这个window不匹配，跳过
        if (overlap_list->list[overlapID].w_list[windowID].y_end == -1)
        {
            continue;
        }

        x_start = overlap_list->list[overlapID].w_list[windowID].x_start;
        y_start = overlap_list->list[overlapID].w_list[windowID].y_start;


        /**
         * There are total 3 cases:
         * 1.  this window of x is overlapped totally by y
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
            ///y_start may less than 0
            y_start = y_start - WINDOW_BOUNDARY/2;

            ///其实可以不加...怕出bug
            if(y_start < 0)
            {
                continue;
            }

            Window_Len = x_len + (threshold << 1);
            determine_overlap_region(threshold, y_start, overlap_list->list[overlapID].y_id, Window_Len, R_INF,
            &extra_begin, &extra_end, &y_start, &o_len);
            fill_subregion(dumy->overlap_region, y_start, o_len, overlap_list->list[overlapID].y_pos_strand, 
            R_INF, overlap_list->list[overlapID].y_id, extra_begin, extra_end);

            x_string = r_string + x_start;
            y_string = dumy->overlap_region;

            ///both end site and real_y_start have extra_begin
            ///有很多是完全匹配，可以先快速判断是不是完全匹配
            end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
                    &(dumy->path_length), dumy->matrix_bit, dumy->path, -1, -1);

            if (error!=(unsigned int)-1)
            {
                total_error = total_error + error;
                matched_coverage++;
                tmp_cigar.x_start = x_start;
                tmp_cigar.x_end = x_end;
                generate_cigar(dumy->path, dumy->path_length, &tmp_cigar, &real_y_start, &end_site, error);
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
                /**
                if(memcmp("m54238_180914_183539/66650355/ccs", Get_NAME((*R_INF), overlap_list->list[overlapID].x_id), 
                Get_NAME_LENGTH((*R_INF), overlap_list->list[overlapID].x_id)) == 0 && uncorrected_window_start == 375)
                {
                    fprintf(stderr, "###### %.*s, x_string[124]: %c, y_string[124]: %c\n", Get_NAME_LENGTH((*R_INF), overlap_list->list[overlapID].y_id), 
                    Get_NAME((*R_INF), overlap_list->list[overlapID].y_id), x_string[124], y_string[124]);

                    fprintf(stderr, "######error: %d****\n", error);
                    for (int ijk = 0; ijk < tmp_cigar.cigar.length; ijk++)
                    {
                        fprintf(stderr, "###### Oper: %d, Len: %d\n", tmp_cigar.cigar.C_C[ijk],
                        tmp_cigar.cigar.C_L[ijk]);
                    }
                    
                    fprintf(stderr, "######dumy->path_length: %d\n****\n", dumy->path_length);
                }
                **/


                /**
                if(verify_cigar(x_string, x_length, y_string, y_length, &tmp_cigar.cigar, 
                    error))
                {
                    fprintf(stderr, "*******error: %d****\n", error);
                    for (int ijk = 0; ijk < tmp_cigar.cigar.length; ijk++)
                    {
                        fprintf(stderr, "Oper: %d, Len: %d\n", tmp_cigar.cigar.C_C[ijk],
                        tmp_cigar.cigar.C_L[ijk]);
                    }
                    
                    fprintf(stderr, "*******dumy->path_length: %d\n****\n", dumy->path_length);
                }
                **/

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
            get_seq_from_Graph(g, &(second_round->dumy), &(second_round->tmp_cigar), backbone);
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
                        UC_Read* g_read, Correct_dumy* dumy, Graph* g, Cigar_record* current_cigar,
                        Round2_alignment* second_round)
{
    clear_Cigar_record(current_cigar);
    
    long long window_num = (g_read->length + WINDOW - 1) / WINDOW;
    long long i, j, overlap_length;
    long long window_start, window_end;

    long long num_availiable_win = 0;
    
    
    window_start = 0;
    window_end = WINDOW - 1;
    if (window_end >= g_read->length)
    {
        window_end = g_read->length - 1;
    }


    int flag;
    ///for last window
    dumy->last_boundary_length = 0;
    for (i = 0; i < window_num; i++)
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
                i = window_num;
                break;
        }


        ///这个是available overlap里所有window的数量...
        ///num_availiable_win = num_availiable_win + dumy->length + dumy->lengthNT;
        num_availiable_win = num_availiable_win + dumy->length;

        ///重叠窗口数，也就是coverage大小
        if(dumy->length >= MIN_COVERAGE_THRESHOLD)
        {

            window_consensus(g_read->seq, window_start, window_end, overlap_list, dumy, R_INF, g, current_cigar);

            if(dumy->last_boundary_length != 0)
            {
                process_boundary(overlap_list, R_INF, dumy, g, current_cigar, window_start, second_round);
                /**
                if(memcmp("m54238_180914_183539/66650355/ccs", Get_NAME((*R_INF), overlap_list->list[0].x_id), 
                Get_NAME_LENGTH((*R_INF), overlap_list->list[0].x_id)) == 0)
                {
                    
                   ///fprintf(stderr, "%.*s\n\n", Get_NAME_LENGTH((*R_INF), overlap_list->list[0].x_id), Get_NAME((*R_INF), overlap_list->list[0].x_id));
                    
                   fprintf(stderr, "window_start: %d, window_end: %d, dumy->last_boundary_length: %d\n", 
                   window_start, window_end, dumy->last_boundary_length);
                   fprintf(stderr, "new_read_length: %d\n",second_round->cigar.new_read_length);
                   fprintf(stderr, "cigar.length: %d\n",second_round->cigar.length);
                   fprintf(stderr, "obtained_cigar_length: %d\n",second_round->obtained_cigar_length);

                    fprintf(stderr, "*******\n");
                   for (size_t ijk = 0; ijk < second_round->cigar.length; ijk++)
                   {
                       int operation = Get_Cigar_Type(second_round->cigar.record[ijk]);
                       int operation_length = Get_Cigar_Length(second_round->cigar.record[ijk]);
                       fprintf(stderr, "oper: %d, oper_len: %d\n", operation, operation_length);
                   }

                   fprintf(stderr, "*******\n");
                   
                }
                **/
            }
            
        }
        else
        {
            add_segment_to_correct_read(dumy, g_read->seq + window_start, window_end - window_start + 1);
            add_cigar_record(g_read->seq + window_start, window_end - window_start + 1, current_cigar, 0);
        }


        

        dumy->last_boundary_length = current_cigar->new_read_length;
        
        window_start = window_start + WINDOW;
        window_end = window_end + WINDOW;
        if (window_end >= g_read->length)
        {
            window_end = g_read->length - 1;
        }
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
    long long inner_offset = x_total_start - window_offset;

    
    ///note that node 0 is the start node
    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2是x缺字符（y多字符），而3是y缺字符（x多字符）
    while (cigar_i < cigar->length)
    {
        operation = cigar->C_C[cigar_i];
        operationLen = cigar->C_L[cigar_i];

        ///这种情况代表匹配和mismatch
        if (operation == 0)
        {
            x_i += operationLen;
            y_i += operationLen;
        }
        else if(operation == 1)
        {
            for (i = 0; i < operationLen; i++)
            {
                /**
                if(inner_offset + x_i >= WINDOW)
                {
                    fprintf(stderr, "error\n");
                }
                **/
                hap->flag[inner_offset + x_i]++;
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
    /**
    if(x_i != x_length || y_i != y_length)
    {
        fprintf(stderr, "x_i: %d, x_length: %d\n", x_i, x_length);
        fprintf(stderr, "y_i: %d, y_length: %d\n", y_i, y_length);
    }
    
   int no_zero = 0;
   for (i = 0; i < WINDOW; i++)
   {
       if(hap->flag[i] != 0)
       {
           no_zero++;
       }
   }

   fprintf(stderr, "no_zero: %d\n", no_zero);
   **/
   
}




void addSNPtohaplotype(
long long window_offset, int overlapID,
char* x_string, long long x_total_start, long long x_length, 
char* y_string, long long y_total_start, long long y_length, 
CIGAR* cigar, haplotype_evdience_alloc* hap)
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
    ///2是x缺字符（y多字符），而3是y缺字符（x多字符）
    while (cigar_i < cigar->length)
    {
        operation = cigar->C_C[cigar_i];
        operationLen = cigar->C_L[cigar_i];

        ///这种情况代表匹配和mismatch
        if (operation == 0)
        {
            for (i = 0; i < operationLen; i++)
            {
                if(hap->flag[inner_offset] > FLAG_THRE)
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

                if(hap->flag[inner_offset] > FLAG_THRE)
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

            for (i = 0; i < operationLen; i++)
            {
                if(hap->flag[inner_offset] > FLAG_THRE)
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
        }
        
        cigar_i++;
    }
}












void cluster(char* r_string, long long window_start, long long window_end, 
overlap_region_alloc* overlap_list, Correct_dumy* dumy, All_reads* R_INF, haplotype_evdience_alloc* hap)
{
    long long x_start;
    long long x_length; 
    char* x_string;
    char* y_string;
    long long i;
    long long y_start, y_length;
    long long overlapID, windowID;
    long long startNodeID, endNodeID, currentNodeID;

    RsetInitHaplotypeEvdienceFlag(hap);

    long long correct_x_pos_s;
    long long inner_window_offset;


    ///与当前window重叠的所有overlap
    ///first mark all snp pos
    for (i = 0; i < dumy->length; i++)
    {
        ///这个是那个overlap的ID，而不是overlap里对应窗口的ID
        overlapID = dumy->overlapID[i];

        ///overlap_list->list[overlapID].x_pos_s is the begining of the whole overlap
        correct_x_pos_s = (overlap_list->list[overlapID].x_pos_s / WINDOW) * WINDOW;
        ///window_start is the begining of this window in the whole x_read
        windowID = (window_start - correct_x_pos_s) / WINDOW;

        ///如果这个window不匹配，跳过
        if (overlap_list->list[overlapID].w_list[windowID].y_end == -1)
        {
            continue;
        }
        
        ///both x_start and y_start are the begining of the whole x_read and y_read
        x_start = overlap_list->list[overlapID].w_list[windowID].x_start;
        x_length = overlap_list->list[overlapID].w_list[windowID].x_end 
                - overlap_list->list[overlapID].w_list[windowID].x_start + 1;

        y_start = overlap_list->list[overlapID].w_list[windowID].y_start;
        y_length = overlap_list->list[overlapID].w_list[windowID].y_end
                - overlap_list->list[overlapID].w_list[windowID].y_start + 1;


        markSNP(window_start, x_start, x_length, y_start, y_length, &(overlap_list->list[overlapID].w_list[windowID].cigar), 
        hap);
    }


    for (i = 0; i < WINDOW; i++)
    {
        if(hap->flag[i] > FLAG_THRE)
        {
            hap->snp++;
        }
    }



    ///add the information related to snp to haplotype_evdience_alloc
    for (i = 0; i < dumy->length; i++)
    {
        ///这个是那个overlap的ID，而不是overlap里对应窗口的ID
        overlapID = dumy->overlapID[i];

        ///overlap_list->list[overlapID].x_pos_s is the begining of the whole overlap
        correct_x_pos_s = (overlap_list->list[overlapID].x_pos_s / WINDOW) * WINDOW;
        ///window_start is the begining of this window in the whole x_read
        windowID = (window_start - correct_x_pos_s) / WINDOW;

        ///如果这个window不匹配，跳过
        if (overlap_list->list[overlapID].w_list[windowID].y_end == -1)
        {
            continue;
        }
        
        ///both x_start and y_start are the begining of the whole x_read and y_read
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
        hap);        
    }



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
            else if(hap->list[i].type == 0)
            {
                if(x_string[0] == y_string[0])
                {
                    fprintf(stderr, "x_string[0]: %c, y_string[0]: %c\n",
                    x_string[0], y_string[0]);
                }

            }
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
overlap_region_alloc* overlap_list, All_reads* R_INF)
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
    
        


    

    
    InsertSNPVector(hap, sub_list, sub_length, s_H[max_i]);
    
    

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

            if (Len_x * OVERLAP_THRESHOLD <=  overlap_list->list[j].align_length)
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

void partition_overlaps(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, haplotype_evdience_alloc* hap)
{
    ResizeInitHaplotypeEvdience(hap);

    long long window_num = (g_read->length + WINDOW - 1) / WINDOW;
    long long i, j, overlap_length;
    long long window_start, window_end;

    long long num_availiable_win = 0;
    

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


        ///这个是available overlap里所有window的数量...
        ///num_availiable_win = num_availiable_win + dumy->length + dumy->lengthNT;
        num_availiable_win = num_availiable_win + dumy->length;


        cluster(g_read->seq, window_start, window_end, overlap_list, dumy, R_INF, hap);

        
        window_start = window_start + WINDOW;
        window_end = window_end + WINDOW;
        if (window_end >= g_read->length)
        {
            window_end = g_read->length - 1;
        }
        
    }


    
    ///very time-consuming
    ///it seems that hap->list has alread been sorted by site
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
                split_sub_list(hap, sub_list, sub_length, hap->snp, overlap_list, R_INF);
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
        split_sub_list(hap, sub_list, sub_length, hap->snp, overlap_list, R_INF);
    }

    ///debug_snp_matrix(hap);

    qsort(hap->snp_stat, hap->available_snp, sizeof(SnpStats), cmp_snp_stats);

    ///debug_snp_matrix(hap);

    if(generate_haplotypes(hap))
    {
        //print_Haplotype(hap, overlap_list, R_INF);
        ///if we can phase, we need to exclude the reads of different haplotype
        int8_t* vector = Get_Result_SNP_Vector((*hap));
        for (j = 0; j < Get_SNP_Vector_Length((*hap)); j++)
        {
            if(vector[j] == 1)
            {
                overlap_list->list[j].align_length = 0;
            }
        }
    }

}























void correct_overlap(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read, Graph* g,
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
    
    partition_overlaps(overlap_list, R_INF, g_read, dumy, hap);

    generate_consensus(overlap_list, R_INF, g_read, dumy, g, current_cigar, second_round);



    /**
    ///fprintf(stderr, "length: %lld, corrected_base: %lld\n", g_read->length, dumy->corrected_base);
    
    EdlibAlignResult result = edlibAlign(g_read->seq, g_read->length, dumy->corrected_read, dumy->corrected_read_length, 
            edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));

    

    if (result.status == EDLIB_STATUS_OK) {
        if (dumy->corrected_base == 0 && result.editDistance!= 0)
        {
            fprintf(stderr, "error 0\n");
        }


        if (dumy->corrected_base < result.editDistance)
        {
            fprintf(stderr, "error 1\n");
        }
        
    
        // fprintf(stderr, "****\n distance: %d, alignmentLength: %d, startLocations: %d, endLocations: %d, corrected_base: %d\n", 
        // result.editDistance, result.alignmentLength, result.startLocations[0], result.endLocations[0], dumy->corrected_base);
        
        
        // char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
        // fprintf(stderr,"%s\n", cigar);
        // free(cigar);
        
    }
    else
    {
        fprintf(stderr, "error\n");
    }
    
    edlibFreeAlignResult(result);

    for (i = 0; i < dumy->corrected_read_length; i++)
    {
        if (dumy->corrected_read[i] != 'A' && dumy->corrected_read[i] != 'C' && dumy->corrected_read[i] != 'G' && dumy->corrected_read[i] != 'T')
        {
            fprintf(stderr, "error\n");
        }
        
    }
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
