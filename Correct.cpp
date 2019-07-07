#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "Correct.h"
#include "Levenshtein_distance.h"
#include "edlib.h"

long long T_total_match=0;
long long T_total_unmatch=0;
long long T_total_mis=0;


#define MAX(x, y) ((x >= y)?x:y)
#define MIN(x, y) ((x <= y)?x:y)
#define OVERLAP(x_start, x_end, y_start, y_end) (MIN(x_end, y_end) - MAX(x_start, y_start) + 1) 
///#define OVERLAP(x_start, x_end, y_start, y_end) MIN(x_end, y_end) - MAX(x_start, y_start) + 1 



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
        ///if (result.editDistance != error && result.endLocations[0] > x_len)
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

            /**
            for (int i = cigar_length - 1; i >= 0; i--)
            {
                switch (cigar[i])
                {
                case 'I':
                    
                    break;
                
                default:
                    break;
                }
            }
            **/
            
            ///fprintf(stderr,"%s\n", cigar);
            free(cigar);
            
            /**
            fprintf(stderr, "****\ni: %u, edlib: %d, alignmentLength: %d, startLocations: %d, endLocations: %d\n", 
            i, result.editDistance, result.alignmentLength, result.startLocations[0], result.endLocations[0]);
            char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
            fprintf(stderr,"%s\n", cigar);
            free(cigar);
            print_string(x_string, x_len);
            print_string(y_string, o_len);

            fprintf(stderr, "BPM: %d\n", error);
            **/
            
        }
    }
    edlibFreeAlignResult(result);
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
    long long total_match=0;
    long long total_unmatch=0;
    long long total_mis=0;


    ///这些是整个window被完全覆盖的
    for (i = 0; i < dumy->length; i++)
    {
        ///整个window被覆盖的话，read本身上的区间就是[window_start, window_end]
        x_len = WINDOW;
        currentID = dumy->overlapID[i];
        x_start = window_start;
        ///y上的相对位置
        y_start = (x_start - overlap_list->list[currentID].x_pos_s) + overlap_list->list[currentID].y_pos_s;
        ///y上的起始
        y_start = y_start - THRESHOLD;
        if (y_start < 0)
        {
            y_start = 0;
        }
        ///当前y的长度
        currentIDLen = Get_READ_LENGTH((*R_INF), overlap_list->list[currentID].y_id);
        ///不能超过y的剩余长度
        o_len = MIN(Window_Len, currentIDLen - y_start);
        
        recover_UC_Read_sub_region(dumy->overlap_region, y_start, o_len, overlap_list->list[currentID].y_pos_strand, 
                R_INF, overlap_list->list[currentID].y_id);

        x_string = r_string + x_start;
        y_string = dumy->overlap_region;

        


        ///这个长度足够，可以用来一起比
        if (Window_Len == o_len)
        {
            
            end_site = Reserve_Banded_BPM(y_string, o_len, x_string, x_len, THRESHOLD, &error);


            if (error!=(unsigned int)-1)
            {
                total_match++;
            }
            else
            {
                total_unmatch++;
            }

            
            test_edit_distance_by_edlib(x_string, y_string, x_len, 
            o_len, THRESHOLD, error, &total_mis);

        }
        else  ///这个不够，只能单个比 ///不够的地方要置N
        {
            /**
            end_site = Reserve_Banded_BPM(y_string, o_len, x_string, x_len, THRESHOLD, &error);


            if (error!=-1)
            {
                total_match++;
            }
            else
            {
                total_unmatch++;
            }
            **/
        }
    }

    long long reverse_i = dumy->size - 1;
    
    ///这些是整个window被部分覆盖的
    for (i = 0; i < dumy->lengthNT; i++)
    {
        currentID = dumy->overlapID[reverse_i--];
        x_start = MAX(window_start, overlap_list->list[currentID].x_pos_s);
        x_end = MIN(window_end, overlap_list->list[currentID].x_pos_e);
        ///这个是和当前窗口重叠的长度
        x_len = x_end - x_start + 1;

        if (x_len <= 0)
        {
            fprintf(stderr, "ERROR\n");
        }

        ///y上的相对位置
        y_start = (x_start - overlap_list->list[currentID].x_pos_s) + overlap_list->list[currentID].y_pos_s;
        ///y上的起始
        y_start = y_start - THRESHOLD;
        if (y_start < 0)
        {
            y_start = 0;
        }
        
        ///当前y的长度
        currentIDLen = Get_READ_LENGTH((*R_INF), overlap_list->list[currentID].y_id);

        
        ///不能超过y的剩余长度
        Window_Len = x_len + (THRESHOLD << 1);
        o_len = MIN(Window_Len, currentIDLen - y_start);
        
        recover_UC_Read_sub_region(dumy->overlap_region, y_start, o_len, overlap_list->list[currentID].y_pos_strand, 
                R_INF, overlap_list->list[currentID].y_id);

        x_string = r_string + x_start;
        y_string = dumy->overlap_region;

        ///不够的地方要置N

        end_site = Reserve_Banded_BPM(y_string, o_len, x_string, x_len, THRESHOLD, &error);

        /**
        if (error!=-1)
        {
            total_match++;
        }
        else
        {
            total_unmatch++;
        }
        **/
    }
    
    /**
    if (total_match!=0)
    {
        
        fprintf(stderr, "total_match: %u\n", total_match);
        fprintf(stderr, "total_unmatch: %u\n", total_unmatch);
        fprintf(stderr, "total_mis: %u\n", total_mis);
        
    }
    **/

    T_total_match = T_total_match + total_match;
    T_total_unmatch = T_total_unmatch + total_unmatch;
    T_total_mis = T_total_mis + total_mis;
    
    


    
    
}

void correct_overlap(overlap_region_alloc* overlap_list, All_reads* R_INF, UC_Read* g_read, Correct_dumy* dumy)
{
    reverse_complement(g_read->seq, g_read->length);

    clear_Correct_dumy(dumy, overlap_list);

    long long window_num = (g_read->length + WINDOW - 1) / WINDOW;

    long long i;

    long long window_start, window_end;

    window_start = 0;
    window_end = WINDOW - 1;



    /**
    UC_Read debug_read;
    init_UC_Read(&debug_read);

    for (i = 0; i < overlap_list->length; i++)
    {
        
        if (overlap_list->list[i].y_pos_strand)
        {
            recover_UC_Read_RC(&debug_read, R_INF, overlap_list->list[i].y_id);
        }
        else
        {
            recover_UC_Read(&debug_read, R_INF, overlap_list->list[i].y_id);
        }
        
        
        long long x_length = overlap_list->list[i].x_pos_e - overlap_list->list[i].x_pos_s + 1;
        long long y_length = overlap_list->list[i].y_pos_e - overlap_list->list[i].y_pos_s + 1;

        EdlibAlignResult result = edlibAlign(g_read->seq+overlap_list->list[i].x_pos_s, 
        x_length, 
        debug_read.seq + overlap_list->list[i].y_pos_s, 
        y_length, 
            edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
            

        if (result.status == EDLIB_STATUS_OK) {
            if (result.editDistance>0 && result.editDistance < x_length*0.02)
            {
                fprintf(stderr, "inner_i: %u, x_length: %u, y_length: %u, editDistance: %u\n", 
                i, x_length, y_length, result.editDistance);
                fprintf(stderr, "x_pos_s: %u, x_pos_e: %u\n", 
                overlap_list->list[i].x_pos_s, overlap_list->list[i].x_pos_e);
                fprintf(stderr, "y_pos_s: %u, y_pos_e: %u, y_pos_strand: %u\n", 
                overlap_list->list[i].y_pos_s, overlap_list->list[i].y_pos_e, overlap_list->list[i].y_pos_strand);


                EdlibAlignResult n_result = edlibAlign(g_read->seq+overlap_list->list[i].x_pos_s, 
                350, 
                debug_read.seq + overlap_list->list[i].y_pos_s - 15, 
                380, 
                    edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
                
                fprintf(stderr, "n_editDistance: %u\n", 
                n_result.editDistance);


                

                edlibFreeAlignResult(n_result);
            }
            
            
        }
        edlibFreeAlignResult(result);
        

    }

    destory_UC_Read(&debug_read);
    **/
    ///return;


















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

        ///break;

    }

    
    fprintf(stderr, "total_match: %u, total_unmatch: %u, total_mis: %u\n", 
    T_total_match, T_total_unmatch, T_total_mis);
    



}


void init_Correct_dumy(Correct_dumy* list)
{
    list->size = 0;
    list->length = 0;
    list->lengthNT = 0;
    list->start_i = 0;
    list->overlapID = NULL;
}

void destory_Correct_dumy(Correct_dumy* list)
{
    free(list->overlapID);
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
    
}