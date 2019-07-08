#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include "Correct.h"
#include "Levenshtein_distance.h"
#include "edlib.h"

long long T_total_match=0;
long long T_total_unmatch=0;
long long T_total_mis=0;
pthread_mutex_t debug_statistics ;


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
    /************需要注释掉********* */ 
    // long long total_match=0;
    // long long total_unmatch=0;
    // long long total_mis=0;
    /************需要注释掉********* */
    int groupLen = 0;
    int return_sites[GROUP_SIZE];
	unsigned int return_sites_error[GROUP_SIZE];

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
        
        
        ///这个长度足够，可以用来一起比
        if (Window_Len == o_len)
        {

            recover_UC_Read_sub_region(dumy->overlap_region_group[groupLen], y_start, o_len, overlap_list->list[currentID].y_pos_strand, 
                R_INF, overlap_list->list[currentID].y_id);

            x_string = r_string + x_start;

            groupLen++;
            if (groupLen == GROUP_SIZE)
            {
            
                Reserve_Banded_BPM_4_SSE_only(dumy->overlap_region_group[0], dumy->overlap_region_group[1], 
                dumy->overlap_region_group[2], dumy->overlap_region_group[3], o_len, x_string, x_len,
			        return_sites, return_sites_error, THRESHOLD, dumy->Peq_SSE);
                groupLen = 0;


                // end_site = Reserve_Banded_BPM(dumy->overlap_region_group[0], o_len, x_string, x_len, THRESHOLD, &error);
                // if (error!=return_sites_error[0])
                // {
                //     fprintf(stderr, "error\n");
                // }
                // end_site = Reserve_Banded_BPM(dumy->overlap_region_group[1], o_len, x_string, x_len, THRESHOLD, &error); 
                // if (error!=return_sites_error[1])
                // {
                //     fprintf(stderr, "error\n");
                // }
                // end_site = Reserve_Banded_BPM(dumy->overlap_region_group[2], o_len, x_string, x_len, THRESHOLD, &error);
                // if (error!=return_sites_error[2])
                // {
                //     fprintf(stderr, "error\n");
                // } 
                // end_site = Reserve_Banded_BPM(dumy->overlap_region_group[3], o_len, x_string, x_len, THRESHOLD, &error);
                // if (error!=return_sites_error[3])
                // {
                //     fprintf(stderr, "error\n");
                // }  
                
               
            //    /************需要注释掉**********/
            //    for (size_t ijk = 0; ijk < GROUP_SIZE; ijk++)
            //    {
            //        if (return_sites_error[ijk]!=(unsigned int)-1)
            //         {
            //             total_match++;
            //         }
            //         else
            //         {
            //             total_unmatch++;
            //         }
            //         test_edit_distance_by_edlib(x_string, dumy->overlap_region_group[ijk], x_len, o_len, THRESHOLD, 
            //         return_sites_error[ijk], &total_mis);
            //    }
            //     /************需要注释掉********* */
                
            }

        }
        else  ///这个不够，只能单个比 ///不够的地方要置N
        {

            recover_UC_Read_sub_region(dumy->overlap_region, y_start, o_len, overlap_list->list[currentID].y_pos_strand, 
                R_INF, overlap_list->list[currentID].y_id);

            x_string = r_string + x_start;
            y_string = dumy->overlap_region;

            ///o_len < Window_Len, 说明y的长度不够，需要在y后面补N
            memset (y_string + o_len, 0, Window_Len - o_len);
            end_site = Reserve_Banded_BPM(y_string, o_len, x_string, x_len, THRESHOLD, &error);
            
            // /************需要注释掉**********/
            // if (error!=(unsigned int)-1)
            // {
            //     total_match++;
            // }
            // else
            // {
            //     total_unmatch++;
            // }
            // test_edit_distance_by_edlib(x_string, y_string, x_len, o_len, THRESHOLD, error, &total_mis);
            // /************需要注释掉********* */
        }
    }



    if (groupLen == 1)
    {
        end_site = Reserve_Banded_BPM(dumy->overlap_region_group[0], o_len, x_string, x_len, THRESHOLD, &error);
        // /************需要注释掉**********/
        // if (error!=(unsigned int)-1)
        // {
        //     total_match++;
        // }
        // else
        // {
        //     total_unmatch++;
        // }
        // test_edit_distance_by_edlib(x_string, dumy->overlap_region_group[0], x_len, o_len, THRESHOLD, error, &total_mis);
        // /************需要注释掉********* */
    }
    else if (groupLen > 1)
    {
        Reserve_Banded_BPM_4_SSE_only(dumy->overlap_region_group[0], dumy->overlap_region_group[1], 
                dumy->overlap_region_group[2], dumy->overlap_region_group[3], o_len, x_string, x_len,
			        return_sites, return_sites_error, THRESHOLD, dumy->Peq_SSE);
        // /************需要注释掉**********/
        // for (size_t ijk = 0; ijk < groupLen; ijk++)
        // {
        //     if (return_sites_error[ijk]!=(unsigned int)-1)
        //     {
        //         total_match++;
        //     }
        //     else
        //     {
        //         total_unmatch++;
        //     }
        //     test_edit_distance_by_edlib(x_string, dumy->overlap_region_group[ijk], x_len, o_len, THRESHOLD, 
        //     return_sites_error[ijk], &total_mis);
        // }
        // /************需要注释掉**********/
        groupLen = 0;
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
        if (Window_Len != o_len)
        {
            ///o_len < Window_Len, 说明y的长度不够，需要在y后面补N
            memset (y_string + o_len, 0, Window_Len - o_len);
        }
        end_site = Reserve_Banded_BPM(y_string, o_len, x_string, x_len, THRESHOLD, &error);


        // /************需要注释掉********* */
        // if (error!=(unsigned int)-1)
        // {
        //     total_match++;
        // }
        // else
        // {
        //     total_unmatch++;
        // }
        // test_edit_distance_by_edlib(x_string, y_string, x_len, o_len, THRESHOLD, error, &total_mis);
        // /************需要注释掉********* */
    }
    



    // /************需要注释掉********* */
    // pthread_mutex_lock(&debug_statistics);
    // T_total_match = T_total_match + total_match;
    // T_total_unmatch = T_total_unmatch + total_unmatch;
    // T_total_mis = T_total_mis + total_mis;
	// pthread_mutex_unlock(&debug_statistics);
    // /************需要注释掉********* */
        
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


    
    // /************需要注释掉********* */
    // pthread_mutex_lock(&debug_statistics);
    // fprintf(stderr, "total_match: %u, total_unmatch: %u, total_mis: %u\n", 
    // T_total_match, T_total_unmatch, T_total_mis);
	// pthread_mutex_unlock(&debug_statistics);
    // /************需要注释掉********* */
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