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
    uint64_t overlapID[GROUP_SIZE];
    uint64_t y_startGroup[GROUP_SIZE];

    ///这些是整个window被完全覆盖的
    for (i = 0; i < dumy->length; i++)
    {
        ///整个window被覆盖的话，read本身上的区间就是[window_start, window_end]
        x_len = WINDOW;
        currentID = dumy->overlapID[i];
        x_start = window_start;
        ///y上的相对位置
        y_start = (x_start - overlap_list->list[currentID].x_pos_s) + overlap_list->list[currentID].y_pos_s;

        // /************需要注释掉**********/
        // if(x_start < overlap_list->list[currentID].x_pos_s)
        // {
        //     fprintf(stderr, "ERROR\n");
        // }
        // /************需要注释掉**********/
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

        // /************需要注释掉**********/
        // if(o_len < x_len)
        // {
        //     fprintf(stderr, "ERROR\n");
        // }
        // /************需要注释掉**********/
        
        
        ///这个长度足够，可以用来一起比
        if (Window_Len == o_len)
        {

            recover_UC_Read_sub_region(dumy->overlap_region_group[groupLen], y_start, o_len, overlap_list->list[currentID].y_pos_strand, 
                R_INF, overlap_list->list[currentID].y_id);
            
            overlapID[groupLen] = currentID;
            y_startGroup[groupLen] = y_start;

            x_string = r_string + x_start;

            groupLen++;
            if (groupLen == GROUP_SIZE)
            {
            
                Reserve_Banded_BPM_4_SSE_only(dumy->overlap_region_group[0], dumy->overlap_region_group[1], 
                dumy->overlap_region_group[2], dumy->overlap_region_group[3], o_len, x_string, x_len,
			        return_sites, return_sites_error, THRESHOLD, dumy->Peq_SSE);
                groupLen = 0;

                
                if (return_sites_error[0]!=(unsigned int)-1)
                {
                    overlap_list->list[overlapID[0]].align_length += x_len;

                    append_window_list(&overlap_list->list[overlapID[0]], window_start, window_end, 
                                y_startGroup[0], y_startGroup[0] + return_sites[0], (int)return_sites_error[0]);
                }
                else
                {
                    append_window_list(&overlap_list->list[overlapID[0]], window_start, window_end, y_startGroup[0], -1, -1);
                }
                

                if (return_sites_error[1]!=(unsigned int)-1)
                {
                    overlap_list->list[overlapID[1]].align_length += x_len;

                    append_window_list(&overlap_list->list[overlapID[1]], window_start, window_end, 
                                y_startGroup[1], y_startGroup[1] + return_sites[1], (int)return_sites_error[1]);
                }
                else
                {
                    append_window_list(&overlap_list->list[overlapID[1]], window_start, window_end, y_startGroup[1], -1, -1);
                }
                

                if (return_sites_error[2]!=(unsigned int)-1)
                {
                    overlap_list->list[overlapID[2]].align_length += x_len;

                    append_window_list(&overlap_list->list[overlapID[2]], window_start, window_end, 
                                y_startGroup[2], y_startGroup[2] + return_sites[2], (int)return_sites_error[2]);
                }
                else
                {
                    append_window_list(&overlap_list->list[overlapID[2]], window_start, window_end, y_startGroup[2], -1, -1);
                }
                

                if (return_sites_error[3]!=(unsigned int)-1)
                {
                    overlap_list->list[overlapID[3]].align_length += x_len;

                    append_window_list(&overlap_list->list[overlapID[3]], window_start, window_end, 
                                y_startGroup[3], y_startGroup[3] + return_sites[3], (int)return_sites_error[3]);
                }
                else
                {
                    append_window_list(&overlap_list->list[overlapID[3]], window_start, window_end, y_startGroup[3], -1, -1);
                }
                

                
               
               /************需要注释掉**********/
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
                /************需要注释掉********* */
                
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

            


            if (error!=(unsigned int)-1)
            {
                overlap_list->list[currentID].align_length += x_len;
                append_window_list(&overlap_list->list[currentID], window_start, window_end, y_start, y_start + end_site, 
                (int)error);
            }
            else
            {
                append_window_list(&overlap_list->list[currentID], window_start, window_end, y_start, -1, -1);
            }
            
            
            /************需要注释掉**********/
            // if (error!=(unsigned int)-1)
            // {
            //     total_match++;
            // }
            // else
            // {
            //     total_unmatch++;
            // }
            // test_edit_distance_by_edlib(x_string, y_string, x_len, o_len, THRESHOLD, error, &total_mis);
            /************需要注释掉********* */
        }
    }



    if (groupLen == 1)
    {
        end_site = Reserve_Banded_BPM(dumy->overlap_region_group[0], o_len, x_string, x_len, THRESHOLD, &error);

        

        if (error!=(unsigned int)-1)
        {
            overlap_list->list[overlapID[0]].align_length += x_len;

            append_window_list(&overlap_list->list[overlapID[0]], window_start, window_end, 
                                y_startGroup[0], y_startGroup[0] + end_site, (int)error);
        }
        else
        {
            append_window_list(&overlap_list->list[overlapID[0]], window_start, window_end, y_startGroup[0], -1, -1);
        }
        

        /************需要注释掉**********/
        // if (error!=(unsigned int)-1)
        // {
        //     total_match++;
        // }
        // else
        // {
        //     total_unmatch++;
        // }
        // test_edit_distance_by_edlib(x_string, dumy->overlap_region_group[0], x_len, o_len, THRESHOLD, error, &total_mis);
        /************需要注释掉********* */
    }
    else if (groupLen > 1)
    {
        Reserve_Banded_BPM_4_SSE_only(dumy->overlap_region_group[0], dumy->overlap_region_group[1], 
                dumy->overlap_region_group[2], dumy->overlap_region_group[3], o_len, x_string, x_len,
			        return_sites, return_sites_error, THRESHOLD, dumy->Peq_SSE);

        for (i = 0; i < groupLen; i++)
        {
            if (return_sites_error[i]!=(unsigned int)-1)
            {
                overlap_list->list[overlapID[i]].align_length += x_len;
                append_window_list(&overlap_list->list[overlapID[i]], window_start, window_end, 
                                y_startGroup[i], y_startGroup[i] + return_sites[i], (int)return_sites_error[i]);
            }
            else
            {
                append_window_list(&overlap_list->list[overlapID[i]], window_start, window_end, y_startGroup[i], -1, -1);
            }
            
        }


        /************需要注释掉**********/
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
        /************需要注释掉**********/
        groupLen = 0;
    }
    
    

    long long reverse_i = dumy->size - 1;
    int threshold;
    
    ///这些是整个window被部分覆盖的
    for (i = 0; i < dumy->lengthNT; i++)
    {
        currentID = dumy->overlapID[reverse_i--];
        x_start = MAX(window_start, overlap_list->list[currentID].x_pos_s);
        x_end = MIN(window_end, overlap_list->list[currentID].x_pos_e);

        ///这个是和当前窗口重叠的长度
        x_len = x_end - x_start + 1;
        threshold = x_len * THRESHOLD_RATE;

        // /************需要注释掉**********/
        // if (x_len <= 0)
        // {
        //     fprintf(stderr, "ERROR\n");
        // }
        // /************需要注释掉**********/

        ///y上的相对位置
        y_start = (x_start - overlap_list->list[currentID].x_pos_s) + overlap_list->list[currentID].y_pos_s;

        // /************需要注释掉**********/
        // if(x_start < overlap_list->list[currentID].x_pos_s)
        // {
        //     fprintf(stderr, "ERROR\n");
        // }
        // /************需要注释掉**********/


        ///y上的起始
        ///y_start = y_start - THRESHOLD;
        y_start = y_start - threshold;
        if (y_start < 0)
        {
            y_start = 0;
        }
        
        ///当前y的长度
        currentIDLen = Get_READ_LENGTH((*R_INF), overlap_list->list[currentID].y_id);

        
        ///不能超过y的剩余长度
        ///Window_Len = x_len + (THRESHOLD << 1);
        Window_Len = x_len + (threshold << 1);
        o_len = MIN(Window_Len, currentIDLen - y_start);


        // /************需要注释掉**********/
        // if(o_len < x_len)
        // {
        //     fprintf(stderr, "ERROR\n");
        // }
        // /************需要注释掉**********/

        
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
        ///end_site = Reserve_Banded_BPM(y_string, o_len, x_string, x_len, THRESHOLD, &error);
        end_site = Reserve_Banded_BPM(y_string, o_len, x_string, x_len, threshold, &error);


        if (error!=(unsigned int)-1)
        {
            overlap_list->list[currentID].align_length += x_len;
            append_window_list(&overlap_list->list[currentID], x_start, x_end, y_start, y_start + end_site, (int)error);
        }
        else
        {
            append_window_list(&overlap_list->list[currentID], x_start, x_end, y_start, -1, -1);
        }




        /************需要注释掉********* */
        // if (error!=(unsigned int)-1)
        // {
        //     total_match++;
        // }
        // else
        // {
        //     total_unmatch++;
        // }
        // //test_edit_distance_by_edlib(x_string, y_string, x_len, o_len, THRESHOLD, error, &total_mis);
        // test_edit_distance_by_edlib(x_string, y_string, x_len, o_len, threshold, error, &total_mis);
        /************需要注释掉********* */
    }
    



    /************需要注释掉********* */
    // pthread_mutex_lock(&debug_statistics);
    // T_total_match = T_total_match + total_match;
    // T_total_unmatch = T_total_unmatch + total_unmatch;
    // T_total_mis = T_total_mis + total_mis;
	// pthread_mutex_unlock(&debug_statistics);
    /************需要注释掉********* */
        
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




inline void generate_cigar(
	char* path, int path_length, window_list* result)
{

    if (result->error == 0)
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
                    fprintf(stderr, "error mismatch, cigar_i: %d, x_i: %d, y_i: %d\n",
                    cigar_i, x_i, y_i);
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
                total_y_start = overlap_list->list[j].w_list[i].y_end + 1;
                ///k遍历匹配window右侧所有不匹配的window
                ///如果i匹配，则k从i+1开始
                ///知道第一个匹配的window结束
                for (k = i + 1; k < overlap_list->list[j].w_list_length && overlap_list->list[j].w_list[k].y_end == -1; k++)
                {   
                    ///y_start有可能大于y_readLen
                    ///这多发于最后一个window长度仅为几，而前面一个window的结束位置也超过了y_readLen-1
                    ///这个时候做动态规划会给超过的部分补N
                    if (total_y_start >= y_readLen)
                    {
                        break;
                    }
                    
                    x_start = overlap_list->list[j].w_list[k].x_start;
                    x_end = overlap_list->list[j].w_list[k].x_end;
                    x_len = x_end - x_start + 1;
                    ///if(x_len == )
                    threshold = x_len * THRESHOLD_RATE;
                    y_start = total_y_start - threshold;
                    if (y_start < 0)
                    {
                        y_start = 0;
                    }
                    Window_Len = x_len + (threshold << 1);
                    ///y_start有可能大于y_readLen
                    ///这多发于最后一个window长度仅为几，而前面一个window的结束位置也超过了y_readLen-1
                    ///这个时候做动态规划会给超过的部分补N
                    y_len = MIN(Window_Len, y_readLen - y_start);
                    ///这说明已经到y的结尾了
                    if (y_len < x_len)
                    {
                        break;
                    }
                    
                    

                    recover_UC_Read_sub_region(dumy->overlap_region, y_start, y_len, y_strand, R_INF, y_id);
                    
                    x_string = g_read->seq + x_start;
                    y_string = dumy->overlap_region;
                    if (Window_Len != y_len)
                    {
                        ///o_len < Window_Len, 说明y的长度不够，需要在y后面补N
                        memset (y_string + y_len, 0, Window_Len - y_len);
                    }
                    
                    end_site = Reserve_Banded_BPM(y_string, y_len, x_string, x_len, threshold, &error);
                    
                    ///error等于-1说明没匹配
                    if (error!=(unsigned int)-1)
                    {
                        overlap_list->list[j].w_list[k].cigar.length = -1;
                        overlap_list->list[j].w_list[k].y_start = y_start;
                        overlap_list->list[j].w_list[k].y_end = y_start + end_site;
                        overlap_list->list[j].w_list[k].error = (int)error;
                        overlap_list->list[j].align_length += x_len;
                    }
                    else
                    {
                        break;
                    }
                    
                    total_y_start = y_start + end_site + 1;

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
                    x_start = overlap_list->list[j].w_list[i].x_start;
                    x_end = overlap_list->list[j].w_list[i].x_end;
                    x_len = x_end - x_start + 1;
                    threshold = x_len * THRESHOLD_RATE;
                    y_start = overlap_list->list[j].w_list[i].y_start;
                    Window_Len = x_len + (threshold << 1);
                    y_len = MIN(Window_Len, y_readLen - y_start);

                    recover_UC_Read_sub_region(dumy->overlap_region, y_start, y_len, y_strand, R_INF, y_id);

                    x_string = g_read->seq + x_start;
                    y_string = dumy->overlap_region;
                    if (Window_Len != y_len)
                    {
                        ///o_len < Window_Len, 说明y的长度不够，需要在y后面补N
                        memset (y_string + y_len, 0, Window_Len - y_len);
                    }

                    end_site = Reserve_Banded_BPM_PATH(y_string, y_len, x_string, x_len, threshold, &error, &real_y_start,
                    &(dumy->path_length), dumy->matrix_bit, dumy->path, 
                    overlap_list->list[j].w_list[i].error, overlap_list->list[j].w_list[i].y_end - y_start);

                    
                    

                    ///到这里y_start已经被正确计算出来了
                    if (error != (unsigned int)-1)
                    {
                        real_y_start = y_start + real_y_start;
                        overlap_list->list[j].w_list[i].y_start = real_y_start;
                        generate_cigar(dumy->path, dumy->path_length, &(overlap_list->list[j].w_list[i]));                                
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
                    x_start = overlap_list->list[j].w_list[k].x_start;
                    x_end = overlap_list->list[j].w_list[k].x_end;
                    x_len = x_end - x_start + 1;
                    threshold = x_len * THRESHOLD_RATE;
                    ///这个和上面不同，先求y_end
                    y_end = total_y_end + threshold;
                    ///y_end不能大于y的总长度
                    if(y_end >= y_readLen)
                    {
                        y_end = y_readLen - 1;
                    }
                    Window_Len = x_len + (threshold << 1);
                    y_len = MIN(Window_Len, y_end + 1);
                    ///这说明已经到y的开始了，因为是倒着算的，没必要接着算了
                    if (y_len < x_len)
                    {
                        break;
                    }
                    y_start = y_end - y_len + 1;
                    recover_UC_Read_sub_region(dumy->overlap_region, y_start, y_len, y_strand, R_INF, y_id);
                    x_string = g_read->seq + x_start;
                    y_string = dumy->overlap_region;
                    
                    ///不要在前面补补，补了有可能触发剪纸
                    ///还是像上面一样在后面随便补点就好了
                    if (Window_Len != y_len)
                    {
                        ///o_len < Window_Len, 说明y的长度不够，需要在y后面补N
                        ///memset (y_string, 0, Window_Len - y_len);
                        memset (y_string + y_len, 0, Window_Len - y_len);
                    }
                    

                    end_site = Reserve_Banded_BPM_PATH(y_string, y_len, x_string, x_len, threshold, &error, &real_y_start,
                    &(dumy->path_length), dumy->matrix_bit, dumy->path, -1, -1);

                    ///error等于-1说明没匹配
                    if (error!=(unsigned int)-1)
                    {
                        
                        ///overlap_list->list[j].w_list[k].y_pre_start = overlap_list->list[j].w_list[k].y_start;
                        overlap_list->list[j].w_list[k].y_start = y_start + real_y_start;
                        overlap_list->list[j].w_list[k].y_end = y_start + end_site;
                        overlap_list->list[j].w_list[k].error = error;
                        overlap_list->list[j].align_length += x_len;
                        generate_cigar(dumy->path, dumy->path_length, &(overlap_list->list[j].w_list[k]));
                    }
                    else
                    {
                        break;
                    }

                    total_y_end = y_start + real_y_start - 1;
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

        if (overlap_length * OVERLAP_THRESHOLD <=  overlap_list->list[j].align_length)
        {
            for (i = 0; i < overlap_list->list[j].w_list_length; i++)
            {
                ///判断cigar是否被计算
                ///没被计算过就重算
                ///第一个条件是判断这个窗口是否匹配
                if(overlap_list->list[j].w_list[i].y_end != -1 && overlap_list->list[j].w_list[i].cigar.length == -1)
                {
                    x_start = overlap_list->list[j].w_list[i].x_start;
                    x_end = overlap_list->list[j].w_list[i].x_end;
                    x_len = x_end - x_start + 1;
                    threshold = x_len * THRESHOLD_RATE;
                    y_start = overlap_list->list[j].w_list[i].y_start;
                    Window_Len = x_len + (threshold << 1);
                    y_len = MIN(Window_Len, y_readLen - y_start);

                    recover_UC_Read_sub_region(dumy->overlap_region, y_start, y_len, y_strand, R_INF, y_id);

                    x_string = g_read->seq + x_start;
                    y_string = dumy->overlap_region;
                    if (Window_Len != y_len)
                    {
                        ///o_len < Window_Len, 说明y的长度不够，需要在y后面补N
                        memset (y_string + y_len, 0, Window_Len - y_len);
                    }

                    end_site = Reserve_Banded_BPM_PATH(y_string, y_len, x_string, x_len, threshold, &error, &real_y_start,
                    &(dumy->path_length), dumy->matrix_bit, dumy->path, 
                    overlap_list->list[j].w_list[i].error, overlap_list->list[j].w_list[i].y_end - y_start);

                    ///到这里y_start已经被正确计算出来了
                    if (error != (unsigned int)-1)
                    {
                        real_y_start = y_start + real_y_start;
                        overlap_list->list[j].w_list[i].y_start = real_y_start;
                        generate_cigar(dumy->path, dumy->path_length, &(overlap_list->list[j].w_list[i]));                                
                    }
                }
            }
        }
    }
    
    















    

    
    // for (j = 0; j < overlap_list->length; j++)
    // {
        
    //     y_id = overlap_list->list[j].y_id;
    //     y_strand = overlap_list->list[j].y_pos_strand;
    //     y_readLen = Get_READ_LENGTH((*R_INF), y_id);

    //     matches = 0;
    //     for (i = 0; i < overlap_list->list[j].w_list_length; i++)
    //     {
    //         x_start = overlap_list->list[j].w_list[i].x_start;
    //         x_end = overlap_list->list[j].w_list[i].x_end;
    //         x_len = x_end - x_start + 1;
    //         ///if(x_len == )
    //         threshold = x_len * THRESHOLD_RATE;   
    //         y_start = overlap_list->list[j].w_list[i].y_start;


            
    //         if (overlap_list->list[j].w_list[i].y_end != -1)
    //         {
    //             matches += x_len;
    //         }

    //         if (i + 1 < overlap_list->list[j].w_list_length)
    //         {
    //             if (
    //                 overlap_list->list[j].w_list[i + 1].x_end <= overlap_list->list[j].w_list[i].x_end
    //                 ||
    //                 overlap_list->list[j].w_list[i + 1].x_start <= overlap_list->list[j].w_list[i].x_start
    //                 )
    //             {
    //                 fprintf(stderr, "ERROR 1\n");
    //             }
                
    //         }


    //         ///当起始位置还没求出来的时候
    //         if(overlap_list->list[j].w_list[i].cigar.length == -1)
    //         {
                
    //             Window_Len = x_len + (threshold << 1);
    //             y_len = MIN(Window_Len, y_readLen - y_start);

    //             recover_UC_Read_sub_region(dumy->overlap_region, y_start, y_len, y_strand, R_INF, y_id);
                    
    //             x_string = g_read->seq + x_start;
    //             y_string = dumy->overlap_region;
    //             if (Window_Len != y_len)
    //             {
    //                 ///o_len < Window_Len, 说明y的长度不够，需要在y后面补N
    //                 memset (y_string + y_len, 0, Window_Len - y_len);
    //             }

    //             ///匹配
    //             if (overlap_list->list[j].w_list[i].y_end != -1)
    //             {
    //                 ///fprintf(stderr, "hahah\n");
    //                 end_site = Reserve_Banded_BPM(y_string, y_len, x_string, x_len, threshold, &error);

    //                 if (error == (unsigned int)-1)
    //                 {
    //                     fprintf(stderr, "1 ERROR, y_start: %u, y_end: %u\n", 
    //                     overlap_list->list[j].w_list[i].y_start, overlap_list->list[j].w_list[i].y_end);
    //                 }
                    
    //                 if (end_site + y_start != overlap_list->list[j].w_list[i].y_end)
    //                 {
    //                     fprintf(stderr, "ERROR, y_end: %u, end_site: %u, y_start: %u, +: %u\n", 
    //                     overlap_list->list[j].w_list[i].y_end, end_site, y_start, end_site + y_start);
    //                 }

    //                 int old_error = overlap_list->list[j].w_list[i].error;

    //                 if (error != old_error)
    //                 {
                        
    //                     fprintf(stderr, "ERROR new error: %d, old error: %d\n", error, overlap_list->list[j].w_list[i].error);
    //                 }
                    
                    

    //             }
    //             else
    //             {
    //                 end_site = Reserve_Banded_BPM(y_string, y_len, x_string, x_len, threshold, &error);

    //                 if (error != (unsigned int)-1)
    //                 {
    //                     fprintf(stderr, "2 ERROR\n");
    //                 }
    //             }
                

    //         }
    //         else
    //         {
                

    //             y_start = overlap_list->list[j].w_list[i].y_start;
    //             ///此时y_start是真正匹配的起始位置，所以要减去threshold
    //             y_start = y_start - threshold;
    //             if(y_start < 0)
    //             {
    //                 y_start = 0;
    //             }

    //             Window_Len = x_len + (threshold << 1);
    //             y_len = MIN(Window_Len, y_readLen - y_start);

    //             recover_UC_Read_sub_region(dumy->overlap_region, y_start, y_len, y_strand, R_INF, y_id);
                    
    //             x_string = g_read->seq + x_start;
    //             y_string = dumy->overlap_region;

    //             if (Window_Len != y_len)
    //             {
    //                 ///o_len < Window_Len, 说明y的长度不够，需要在y后面补N
    //                 memset (y_string + y_len, 0, Window_Len - y_len);
    //             }


                
    //             end_site = Reserve_Banded_BPM_PATH(y_string, y_len, x_string, x_len, threshold, &error, &real_y_start,
    //                     &(dumy->path_length), dumy->matrix_bit, dumy->path, 
    //                     overlap_list->list[j].w_list[i].error, overlap_list->list[j].w_list[i].y_end - y_start);
                        
    //             ///end_site = Reserve_Banded_BPM(y_string, y_len, x_string, x_len, threshold, &error);
                

    //             if (error == (unsigned int)-1)
    //             {
    //                 fprintf(stderr, "3 ERROR, y_start: %d, y_end: %d, x_len: %d, pre_error: %d, y_strand: %u\n", 
    //                 overlap_list->list[j].w_list[i].y_start, overlap_list->list[j].w_list[i].y_end, x_len, overlap_list->list[j].w_list[i].error, y_strand);

    //                 fprintf(stderr, "Window_Len: %d, y_len: %d, y_readLen: %d\n", 
    //                 Window_Len, y_len, y_readLen);

    //                 print_string(g_read->seq+overlap_list->list[j].w_list[i].x_start, 
    //                 overlap_list->list[j].w_list[i].x_end-overlap_list->list[j].w_list[i].x_start+1);

    //                 recover_UC_Read_sub_region(dumy->overlap_region, overlap_list->list[j].w_list[i].y_start, 
    //                 overlap_list->list[j].w_list[i].y_end-overlap_list->list[j].w_list[i].y_start + 1, y_strand, R_INF, y_id);

    //                 print_string(dumy->overlap_region, 
    //                 overlap_list->list[j].w_list[i].y_end-overlap_list->list[j].w_list[i].y_start + 1);
    //             }

    //             ////检验cigar是否正确
    //             x_start = overlap_list->list[j].w_list[i].x_start;
    //             x_end = overlap_list->list[j].w_list[i].x_end;
    //             x_len = x_end - x_start + 1;
    //             x_string = g_read->seq + x_start;

    //             y_start = overlap_list->list[j].w_list[i].y_start;
    //             y_end = overlap_list->list[j].w_list[i].y_end;
    //             y_len = y_end - y_start + 1;
    //             recover_UC_Read_sub_region(dumy->overlap_region, y_start, y_len, y_strand, R_INF, y_id);
    //             y_string = dumy->overlap_region;
                
                
    //             if(verify_cigar(x_string, x_len, y_string, y_len, &overlap_list->list[j].w_list[i].cigar, 
    //             overlap_list->list[j].w_list[i].error))
    //             {
    //                 fprintf(stderr, "j: %d, i: %d, y_id: %d, y_start: %d, y_end: %d\n", j, i, y_id, y_start, y_end);
    //             }
                
                
                
            
                

    //         }
            
    //     }

    //     if (matches != overlap_list->list[j].align_length)
    //     {
    //         fprintf(stderr, "ERROR 2: matches: %d, align_length: %d\n", 
    //                                 matches, overlap_list->list[j].align_length);
    //     }
    // }
    
    
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



void window_consensus(char* r_string, long long window_start, long long window_end, 
overlap_region_alloc* overlap_list, Correct_dumy* dumy, All_reads* R_INF, Graph* g)
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


    get_seq_from_Graph(g, startNodeID, endNodeID, dumy);


    /**
    if (startNodeID != 0 || endNodeID != backbone_length - 1)
    {
        fprintf(stderr, "error\n");
    }
    

    for (i = startNodeID; i <= backbone_length; i++)
    {
        if(g->g_nodes.list[i].alignedTo_Nodes.length > 5)
        {
            fprintf(stderr, "error 1\n");
        }

        if(g->g_nodes.list[i].outcome_edges.length > 4)
        {
            fprintf(stderr, "error 2\n");
        }
    }

    
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


void generate_consensus(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, Graph* g)
{

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
            window_consensus(g_read->seq, window_start, window_end, overlap_list, dumy, R_INF, g);
        }
        else
        {
            add_segment_to_correct_read(dumy, g_read->seq + window_start, window_end - window_start + 1);
        }
        
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
    }
    





    /***********************要注释掉*************************/
    // long long debug_num_availiable_win = 0;
    // for (j = 0; j < overlap_list->length; j++)
    // {
    //     overlap_length = overlap_list->list[j].x_pos_e - overlap_list->list[j].x_pos_s + 1;

    //     if (overlap_length * OVERLAP_THRESHOLD <=  overlap_list->list[j].align_length)
    //     {
    //         debug_num_availiable_win = debug_num_availiable_win + overlap_list->list[j].w_list_length;
    //     }
    // }

    // if (debug_num_availiable_win != num_availiable_win)
    // {
    //     fprintf(stderr, "error, debug_num_availiable_win: %d, num_availiable_win: %d\n",
    //     debug_num_availiable_win, num_availiable_win);
    // }
    /***********************要注释掉*************************/
}

void correct_overlap(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read, Graph* g,
                        long long* matched_overlap_0, long long* matched_overlap_1, 
                        long long* potiental_matched_overlap_0, long long* potiental_matched_overlap_1)
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


    generate_consensus(overlap_list, R_INF, g_read, dumy, g);


    ///fprintf(stderr, "length: %lld, corrected_base: %lld\n", g_read->length, dumy->corrected_base);
    /**
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
    **/
    
    /************需要注释掉********* */
    // pthread_mutex_lock(&debug_statistics);
    // fprintf(stderr, "total_match: %u, total_unmatch: %u, total_mis: %u\n", 
    // T_total_match, T_total_unmatch, T_total_mis);
	// pthread_mutex_unlock(&debug_statistics);
    /************需要注释掉********* */

    
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
