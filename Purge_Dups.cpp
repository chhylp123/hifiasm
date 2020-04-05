#include "Purge_Dups.h"
#include "Overlaps.h"
#include "Correct.h"
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <stdlib.h>
#include "ksort.h"

#define Cal_Off(OFF) ((long long)((uint32_t)((OFF)>>32)) - (long long)((uint32_t)((OFF))))
#define Get_match(x) ((x).weight)
#define Get_total(x) ((x).index_beg)
#define Get_type(x) ((x).index_end)
#define Get_x_beg(x) ((x).x_beg_pos)
#define Get_x_end(x) ((x).x_end_pos)
#define Get_y_beg(x) ((x).y_beg_pos)
#define Get_y_end(x) ((x).y_end_pos)
#define Get_rev(x) ((x).rev)

uint8_t debug_enable = 0;

typedef struct {
	asg_arc_t x;
    uint64_t Off;
    uint64_t weight;
}asg_arc_t_offset;

typedef struct {
	kvec_t(asg_arc_t_offset) a;
	uint64_t i;
}kvec_asg_arc_t_offset;


typedef struct {
    uint64_t weight;
    uint32_t x_beg_pos;
    uint32_t x_end_pos;
    uint32_t y_beg_pos;
    uint32_t y_end_pos;
    uint32_t index_beg;
    uint32_t index_end;
    uint8_t rev;
}hap_candidates;

typedef struct {
    kvec_t(hap_candidates) a;
    uint64_t i;
}kvec_hap_candidates;

#define SELF_EXIST 0
#define REVE_EXIST 1
#define DELETE 2

typedef struct {
    uint8_t rev;
    uint8_t type;
    uint8_t status;
    uint32_t x_beg_pos;
    uint32_t x_end_pos;
    uint32_t y_beg_pos;
    uint32_t y_end_pos;
    uint32_t x_beg_id;
    uint32_t x_end_id;
    uint32_t y_beg_id;
    uint32_t y_end_id;
    uint32_t xUid;
    uint32_t yUid;
    uint32_t weight;
}hap_overlaps;

typedef struct {
    kvec_t(hap_overlaps) a;
}kvec_hap_overlaps;

typedef struct {
    kvec_hap_overlaps* x;
    uint32_t num;
}hap_overlaps_list;

void init_hap_overlaps_list(hap_overlaps_list* x, uint32_t num)
{
    uint32_t i = 0;
    x->num = num;
    x->x = (kvec_hap_overlaps*)malloc(sizeof(kvec_hap_overlaps)*x->num);
    for (i = 0; i < x->num; i++)
    {
        kv_init(x->x[i].a);
    }
}

void enable_debug_mode()
{
    debug_enable = 1;
}

void destory_hap_overlaps_list(hap_overlaps_list* x)
{
    uint32_t i = 0;
    for (i = 0; i < x->num; i++)
    {
        kv_destroy(x->x[i].a);
    }
    free(x->x);
}



inline void clean_visit_flag(uint8_t* visit, asg_t *read_g, R_to_U* ruIndex, uint32_t contigNum, 
ma_hit_t_alloc* x)
{
    uint32_t k, rId, is_Unitig, Hap_cId;

    if(x->length*2 > contigNum)
    {
        memset(visit, 0, contigNum);
    }
    else
    {
        for (k = 0; k < x->length; k++)
        {
            rId = Get_tn(x->buffer[k]);
        
            if(read_g->seq[rId].del == 1)
            {
                ///get the id of read that contains it 
                get_R_to_U(ruIndex, rId, &rId, &is_Unitig);
                if(rId == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[rId].del == 1) continue;
            }

            ///there are two cases: 
            ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
            ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
            get_R_to_U(ruIndex, rId, &Hap_cId, &is_Unitig);
            if(is_Unitig == 0 || Hap_cId == (uint32_t)-1) continue;
            ///here rId is the id of the read coming from the different haplotype
            ///Hap_cId is the id of the corresponding contig (note here is the contig, instead of untig)
            visit[Hap_cId] = 0;
        }
    }
}

uint32_t prefilter(uint32_t x_pos, uint32_t y_pos, uint32_t xLen, uint32_t yLen, uint32_t dir, 
float Hap_rate, uint32_t seedOcc)
{
    uint32_t max_count = 0, min_count = 0;
    uint32_t /**xLeftBeg, **/xLeftLen, yLeftBeg, yLeftLen;
    uint32_t xRightBeg, xRightLen, yRightBeg, yRightLen;
    if(dir == 0)
    {
        /**xLeftBeg = 0;**/ xLeftLen = x_pos; xRightBeg = x_pos; xRightLen = xLen - xRightBeg;
        yLeftBeg = 0; yLeftLen = y_pos; yRightBeg = y_pos; yRightLen = yLen - yRightBeg;
    }
    else
    {
        /**xLeftBeg = 0;**/ xLeftLen = x_pos; xRightBeg = x_pos; xRightLen = xLen - xRightBeg;

        yLeftBeg = y_pos + 1; yLeftLen = yLen - yLeftBeg;
        yRightBeg = 0; yRightLen = y_pos + 1;
    }



    max_count = seedOcc;
    min_count = MIN(xLeftLen, yLeftLen) + MIN(xRightLen, yRightLen);
    if(min_count == 0) return NON_PLOID;
    if(max_count <= min_count*Hap_rate) return NON_PLOID;
    return PLOID;
}



inline uint64_t get_xy_pos(asg_t *read_g, asg_arc_t* t, uint32_t v_in_unitig, uint32_t w_in_unitig, 
uint32_t xUnitigLen, uint32_t yUnitigLen, uint64_t* position_index, uint8_t* rev)
{
    uint32_t x_pos, y_pos, x_dir = 0, y_dir = 0;
    uint64_t tmp;
    x_pos = y_pos = (uint32_t)-1;
    if((t->ul>>32)==v_in_unitig)///end pos
    {
        x_pos = (position_index[v_in_unitig>>1]>>32) + read_g->seq[v_in_unitig>>1].len - 1;
        x_dir = 0;
    }
    else if((t->ul>>32)==(v_in_unitig^1))///start pos
    {
        x_pos = (position_index[v_in_unitig>>1]>>32);
        x_dir = 1;
    }
    else
    {
        fprintf(stderr, "ERROR\n");
    }
    
    if(t->v == w_in_unitig)
    {
        y_pos = (position_index[w_in_unitig>>1]>>32) + t->ol - 1;
        y_dir = 0;
    }
    else if(t->v == (w_in_unitig^1))
    {
        y_pos = (position_index[w_in_unitig>>1]>>32) + read_g->seq[w_in_unitig>>1].len - t->ol;
        y_dir = 1;
    }
    else
    {
        fprintf(stderr, "ERROR\n");
    }

    (*rev) = x_dir^y_dir;
    if((*rev)) y_pos = yUnitigLen - y_pos - 1;
    
    tmp = x_pos; tmp = tmp << 32; tmp = tmp | y_pos;
    return tmp;
}


void print_debug_unitig(ma_utg_t *xReads, uint64_t* position_index, const char* infor)
{
    uint32_t k;
    fprintf(stderr, "\n%s: n = %u\n", infor, xReads->n);
    for (k = 0; k < xReads->n; k++)
    {
        fprintf(stderr, "(%u)v: %u, len: %u, index: %u, pos: %u\n", 
        k, (uint32_t)(xReads->a[k]>>32), (uint32_t)xReads->a[k], (uint32_t)(position_index[xReads->a[k]>>33]),
        (uint32_t)(position_index[xReads->a[k]>>33]>>32));
    }
}

void deduplicate_edge(kvec_asg_arc_t_offset* u_buffer)
{
    if(u_buffer->a.n == 0) return;
    long long i = u_buffer->a.n - 1, k, i_off;
    uint32_t v = u_buffer->a.a[i].x.ul>>33, m;

    for (; i >= 0; i--)
    {
        if((u_buffer->a.a[i].x.ul>>33) != v)
        {
            break;
        } 
    }

    ///fprintf(stderr, "u_buffer->a.n: %u, i: %lld\n", u_buffer->a.n, i);

    i = i + 1;
    for (m = i; i < (long long)u_buffer->a.n; i++)
    {
        if(u_buffer->a.a[i].x.del) continue;

        i_off = Cal_Off(u_buffer->a.a[i].Off);

        for (k = i + 1; k < (long long)u_buffer->a.n; k++)
        {
            if(u_buffer->a.a[k].x.del) continue;
            if(u_buffer->a.a[i].x.el != u_buffer->a.a[k].x.el) continue;
            if(i_off != Cal_Off(u_buffer->a.a[k].Off)) continue;
            u_buffer->a.a[k].x.del = 1;
            u_buffer->a.a[i].weight += u_buffer->a.a[k].weight;
        }
        u_buffer->a.a[m] = u_buffer->a.a[i];
        m++;
    }

    u_buffer->a.n = m;
    ///fprintf(stderr, "u_buffer->a.n: %u, i: %lld\n", u_buffer->a.n, i);
}

int cmp_hap_alignment(const void * a, const void * b)
{
    if((*(asg_arc_t_offset*)a).x.el > (*(asg_arc_t_offset*)b).x.el) return 1;
    if((*(asg_arc_t_offset*)a).x.el < (*(asg_arc_t_offset*)b).x.el) return -1;

    long long aOff = Cal_Off((*(asg_arc_t_offset*)a).Off);
    long long bOff = Cal_Off((*(asg_arc_t_offset*)b).Off);

    if(aOff > bOff) return 1;
    if(aOff < bOff) return -1;

    if(((*(asg_arc_t_offset*)a).Off>>32) > ((*(asg_arc_t_offset*)b).Off>>32)) return 1;
    if(((*(asg_arc_t_offset*)a).Off>>32) < ((*(asg_arc_t_offset*)b).Off>>32)) return -1;

    if((uint32_t)((*(asg_arc_t_offset*)a).Off) > (uint32_t)((*(asg_arc_t_offset*)b).Off)) return 1;
    if((uint32_t)((*(asg_arc_t_offset*)a).Off) < (uint32_t)((*(asg_arc_t_offset*)b).Off)) return -1;

    if((*(asg_arc_t_offset*)a).weight < (*(asg_arc_t_offset*)b).weight) return 1;
    if((*(asg_arc_t_offset*)a).weight > (*(asg_arc_t_offset*)b).weight) return -1;

    return 0;
}

int cmp_hap_candidates(const void * a, const void * b)
{
    if((*(hap_candidates*)a).weight < (*(hap_candidates*)b).weight) return 1;
    if((*(hap_candidates*)a).weight > (*(hap_candidates*)b).weight) return -1;

    if((*(hap_candidates*)a).index_beg > (*(hap_candidates*)b).index_beg) return 1;
    if((*(hap_candidates*)a).index_beg < (*(hap_candidates*)b).index_beg) return -1;
    
    return 0;
}
inline long long get_hap_overlapLen(long long x_beg, long long x_end, long long xLen, 
long long y_beg, long long y_end, long long yLen, long long* n_x_beg, long long* n_x_end,
long long* n_y_beg, long long* n_y_end)
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

    if(n_x_beg) (*n_x_beg) = x_beg;
    if(n_x_end) (*n_x_end) = x_end;
    if(n_y_beg) (*n_y_beg) = y_beg;
    if(n_y_end) (*n_y_end) = y_end;

    return x_end - x_beg + 1;
}


#define X2Y 0
#define Y2X 1
#define XCY 2
#define YCX 3
uint32_t classify_hap_overlap(long long xBeg, long long xEnd, long long xLen,
long long yBeg, long long yEnd, long long yLen, long long* r_xBeg, long long* r_xEnd, 
long long* r_yBeg, long long* r_yEnd)
{
    long long n_x_beg, n_x_end, n_y_beg, n_y_end;
    get_hap_overlapLen(xBeg, xEnd, xLen, yBeg, yEnd, yLen, &n_x_beg, &n_x_end, &n_y_beg, &n_y_end);
    if(r_xBeg) (*r_xBeg) = n_x_beg;
    if(r_xEnd) (*r_xEnd) = n_x_end;
    if(r_yBeg) (*r_yBeg) = n_y_beg;
    if(r_yEnd) (*r_yEnd) = n_y_end;
    if(n_x_beg == 0 && n_x_end == xLen - 1) return YCX;
    if(n_y_beg == 0 && n_y_end == yLen - 1) return XCY;
    if(n_y_beg == 0 && n_x_end == xLen - 1) return X2Y;
    if(n_x_beg == 0 && n_y_end == yLen - 1) return Y2X;
    return XCY;
}

void get_pair_hap_similarity(uint64_t* readIDs, uint32_t Len, uint32_t target_uId, 
ma_hit_t_alloc* reverse_sources, asg_t *read_g, R_to_U* ruIndex, double* Match, double* Total)
{
    #define CUTOFF_THRES 100
    uint32_t i, j, qn, tn, is_Unitig, uId, min_count = 0, max_count = 0, cutoff = 0;;
    for (i = 0; i < Len; i++)
    {
        if(cutoff > CUTOFF_THRES)
        {
            max_count = 0;
            min_count = Len;
        }
        qn = readIDs[i]>>33;
        if(reverse_sources[qn].length > 0) min_count++;
        for (j = 0; j < reverse_sources[qn].length; j++)
        {
            tn = Get_tn(reverse_sources[qn].buffer[j]);
            if(read_g->seq[tn].del == 1)
            {
                get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
                if(tn == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[tn].del == 1) continue;
            }


            get_R_to_U(ruIndex, tn, &uId, &is_Unitig);
            if(uId!=(uint32_t)-1 && is_Unitig == 1 && uId == target_uId)
            {
                max_count++;
                break;
            }
        }

        //means no match
        if(j == reverse_sources[qn].length)
        {
            cutoff++;
        }
        else
        {
            cutoff = 0;
        }
    }

    (*Match) = max_count;
    (*Total) = min_count;
}

inline void check_hap_match(uint32_t qn, uint32_t targetBeg, uint32_t targetEnd, uint32_t targetID, 
uint64_t* position_index, ma_hit_t_alloc* reverse_sources, asg_t *read_g, R_to_U* ruIndex, uint32_t* is_found, uint32_t* is_match)
{
    uint32_t j, tn, uId, is_Unitig, offset;
    (*is_found) = (*is_match) = 0;
    if(reverse_sources[qn].length > 0) (*is_found) = 1;
    for (j = 0; j < reverse_sources[qn].length; j++)
    {
        tn = Get_tn(reverse_sources[qn].buffer[j]);
        if(read_g->seq[tn].del == 1)
        {
            get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
            if(tn == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[tn].del == 1) continue;
        }


        get_R_to_U(ruIndex, tn, &uId, &is_Unitig);
        if(uId!=(uint32_t)-1 && is_Unitig == 1 && uId == targetID)
        {
            offset = (uint32_t)(position_index[tn]);
            if(offset >= targetBeg && offset <= targetEnd)
            {
                (*is_match) = 1;
                break;
            }                    
        }
    }

}


void determin_hap_alignment_boundary_single_side(uint64_t* readIDs, long long queryLen, long long targetBeg, 
long long targetEnd, long long targetID, long long eMatch, long long eTotal, long long dir, 
float Hap_rate, uint64_t* position_index, ma_hit_t_alloc* reverse_sources, asg_t *read_g, 
R_to_U* ruIndex, uint32_t* n_matchLen, uint32_t* n_max_count, uint32_t* n_min_count)
{
    long long i, maxId, min_count = eTotal, max_count = eMatch, matchLen = 0;
    uint32_t is_found, is_match;
    if(dir == 0)
    {
        for (i = 0, maxId = 0; i < queryLen; i++)
        {
            check_hap_match(readIDs[i]>>33, targetBeg, targetEnd, targetID, position_index, reverse_sources, 
            read_g, ruIndex, &is_found, &is_match);
            min_count += is_found;
            max_count += is_match;
            if(max_count > min_count*Hap_rate) maxId = i;
        }


        for (i = maxId; i >= 0; i--)
        {
            check_hap_match(readIDs[i]>>33, targetBeg, targetEnd, targetID, position_index, reverse_sources, 
            read_g, ruIndex, &is_found, &is_match);
            if(is_found == 1 && is_match == 1)
            {
                break;
            }
            min_count -= is_found;
            max_count -= is_match;
        }

        matchLen = i+1;
    }
    else
    {
        for (i = queryLen - 1, maxId = queryLen - 1; i >= 0; i--)
        {
            check_hap_match(readIDs[i]>>33, targetBeg, targetEnd, targetID, position_index, reverse_sources, 
            read_g, ruIndex, &is_found, &is_match);
            min_count += is_found;
            max_count += is_match;
            if(max_count > min_count*Hap_rate) maxId = i;
        }

        for (i = maxId; i < queryLen; i++)
        {
            check_hap_match(readIDs[i]>>33, targetBeg, targetEnd, targetID, position_index, reverse_sources, 
            read_g, ruIndex, &is_found, &is_match);
            if(is_found == 1 && is_match == 1)
            {
                break;
            }
            min_count -= is_found;
            max_count -= is_match;
        }

        matchLen = queryLen - i;
    }

    ///need to check if min_count == 0
    if(min_count == 0)
    {
        (*n_matchLen) = (*n_min_count) = (*n_max_count) = 0;
        return;
    }

    (*n_matchLen) = matchLen;
    (*n_min_count) = min_count;
    (*n_max_count) = max_count;
}

inline void modify_target_interval(long long beg, long long end, long long len,
long long* target_beg, long long* target_end)
{
    #define TARGET_SGIFT 3
    beg -= TARGET_SGIFT;
    end += TARGET_SGIFT;
    if(beg < 0) beg = 0;
    if(end >= len) end = len - 1;
    (*target_beg) = beg;
    (*target_end) = end;
}


void bi_direction_hap_alignment_extention(ma_utg_t* xReads, uint32_t xLeftBeg, uint32_t xLeftLen,
uint32_t xRightBeg, uint32_t xRightLen, uint32_t targetUid, uint32_t target_beg, uint32_t target_end,
float Hap_rate, uint64_t* position_index, ma_hit_t_alloc* reverse_sources, asg_t *read_g, R_to_U* ruIndex,
uint32_t rev, long long* x_interval_beg, long long* x_interval_end)
{
    if(rev)
    {
        uint32_t k;
        k = xLeftBeg; xLeftBeg = xRightBeg; xRightBeg = k;
        k = xLeftLen; xLeftLen = xRightLen; xRightLen = k;
    }
    uint32_t n_matchLenLeft, x_max_countLeft, x_min_countLeft;
    uint32_t n_matchLenRight, x_max_countRight, x_min_countRight;
    n_matchLenLeft = x_max_countLeft = x_min_countLeft = 0;
    determin_hap_alignment_boundary_single_side(xReads->a+xLeftBeg, xLeftLen,
    target_beg, target_end, targetUid, x_max_countLeft, x_min_countLeft, 1, Hap_rate, 
    position_index, reverse_sources, read_g, ruIndex, &n_matchLenLeft, &x_max_countLeft, 
    &x_min_countLeft);

    n_matchLenRight = x_max_countRight = x_min_countRight = 0;
    determin_hap_alignment_boundary_single_side(xReads->a+xRightBeg, xRightLen,
    target_beg, target_end, targetUid, x_max_countRight, x_min_countRight, 0, Hap_rate,
    position_index, reverse_sources, read_g, ruIndex, &n_matchLenRight, &x_max_countRight,
    &x_min_countRight);

    if(x_max_countLeft >= x_max_countRight)
    {
        determin_hap_alignment_boundary_single_side(xReads->a+xRightBeg, xRightLen,
        target_beg, target_end, targetUid, x_max_countLeft, x_min_countLeft, 0, Hap_rate,
        position_index, reverse_sources, read_g, ruIndex, &n_matchLenRight, &x_max_countRight,
        &x_min_countRight);
    }
    else
    {
        determin_hap_alignment_boundary_single_side(xReads->a+xLeftBeg, xLeftLen,
        target_beg, target_end, targetUid, x_max_countRight, x_min_countRight, 1, Hap_rate, 
        position_index, reverse_sources, read_g, ruIndex, &n_matchLenLeft, &x_max_countLeft, 
        &x_min_countLeft);
    }

    (*x_interval_beg) = xLeftBeg + xLeftLen; (*x_interval_beg) -= n_matchLenLeft;
    (*x_interval_end) = xRightBeg + n_matchLenRight; (*x_interval_end) -= 1;
}

void get_hap_alignment_boundary(ma_utg_t* xReads, ma_utg_t* yReads, uint32_t type, 
uint32_t xLeftMatch, uint32_t xLeftTotal, uint32_t yLeftMatch, uint32_t yLeftTotal, 
uint32_t xRightMatch, uint32_t xRightTotal, uint32_t yRightMatch, uint32_t yRightTotal,
uint32_t xLeftBeg, uint32_t xLeftLen, uint32_t yLeftBeg, uint32_t yLeftLen,
uint32_t xRightBeg, uint32_t xRightLen, uint32_t yRightBeg, uint32_t yRightLen,
uint32_t xUid, uint32_t yUid, float Hap_rate, uint64_t* position_index,
ma_hit_t_alloc* reverse_sources, asg_t *read_g, R_to_U* ruIndex, uint32_t rev,
long long* r_x_interval_beg, long long* r_x_interval_end, 
long long* r_y_interval_beg, long long* r_y_interval_end)
{
    
    uint32_t x_max_count, x_min_count, y_max_count, y_min_count, n_matchLen;
    long long x_interval_beg, x_interval_end, y_interval_beg, y_interval_end;
    long long target_beg, target_end;
    x_max_count = x_min_count = y_max_count = y_min_count = 0;

    if(type == X2Y)
    {
        /********************x*********************/
        x_max_count = xRightMatch;
        x_min_count = xRightTotal;

        modify_target_interval(yLeftBeg, yLeftBeg+yLeftLen-1, yReads->n, &target_beg, &target_end);
        determin_hap_alignment_boundary_single_side(xReads->a+xLeftBeg, xLeftLen, 
        /**yLeftBeg, yLeftBeg+yLeftLen-1,**/ target_beg, target_end, yUid, 
        x_max_count, x_min_count, 1, Hap_rate, position_index, reverse_sources, 
        read_g, ruIndex, &n_matchLen, &x_max_count, &x_min_count);

        x_interval_beg = xLeftBeg + xLeftLen; x_interval_beg -= n_matchLen;
        x_interval_end = xRightBeg + xRightLen; x_interval_end -= 1;
        /********************x*********************/

        /********************y*********************/
        y_max_count = yLeftMatch;
        y_min_count = yLeftTotal;
        modify_target_interval(xRightBeg, xRightBeg+xRightLen-1, xReads->n, &target_beg, &target_end);
        determin_hap_alignment_boundary_single_side(yReads->a+yRightBeg, yRightLen,
        /**xRightBeg, xRightBeg+xRightLen-1,**/ target_beg, target_end, xUid, 
        y_max_count, y_min_count, rev, Hap_rate, position_index, reverse_sources, 
        read_g, ruIndex, &n_matchLen, &y_max_count, &y_min_count);
        if(rev == 0)
        {
            y_interval_beg = yLeftBeg;
            y_interval_end = yRightBeg + n_matchLen; y_interval_end -= 1;
        }
        else
        {
            y_interval_beg = yRightBeg + yRightLen; y_interval_beg -= n_matchLen;
            y_interval_end = yLeftBeg + yLeftLen; y_interval_end -= 1;
        }
        
        /********************y*********************/
    }
    else if(type == Y2X)
    {
        /********************x*********************/
        x_max_count = xLeftMatch;
        x_min_count = xLeftTotal;
        modify_target_interval(yRightBeg, yRightBeg+yRightLen-1, yReads->n, &target_beg, &target_end);
        determin_hap_alignment_boundary_single_side(xReads->a+xRightBeg, xRightLen,
        /**yRightBeg, yRightBeg+yRightLen-1,**/ target_beg, target_end, yUid, 
        x_max_count, x_min_count, 0, Hap_rate, position_index, reverse_sources, 
        read_g, ruIndex, &n_matchLen, &x_max_count, &x_min_count);

        x_interval_beg = xLeftBeg;
        x_interval_end = xRightBeg + n_matchLen; x_interval_end -= 1;
        /********************x*********************/

        /********************y*********************/
        y_max_count = yRightMatch;
        y_min_count = yRightTotal;
        modify_target_interval(xLeftBeg, xLeftBeg+xLeftLen-1, xReads->n, &target_beg, &target_end);
        determin_hap_alignment_boundary_single_side(yReads->a+yLeftBeg, yLeftLen,
        /**xLeftBeg, xLeftBeg+xLeftLen-1,**/ target_beg, target_end, xUid, 
        y_max_count, y_min_count, 1-rev, Hap_rate, position_index, reverse_sources, 
        read_g, ruIndex, &n_matchLen, &y_max_count, &y_min_count);
        if(rev == 0)
        {
            y_interval_beg = yLeftBeg + yLeftLen; y_interval_beg -= n_matchLen;
            y_interval_end = yRightBeg + yRightLen; y_interval_end -= 1;
        }
        else
        {
            y_interval_beg = yRightBeg;
            y_interval_end = yLeftBeg + n_matchLen; y_interval_end -= 1;
        }
        /********************y*********************/
    }
    else if(type == XCY)
    {
        /********************x*********************/   
        bi_direction_hap_alignment_extention(xReads, xLeftBeg, xLeftLen, xRightBeg, xRightLen, 
        yUid, 0, yReads->n - 1, Hap_rate, position_index, reverse_sources, read_g, ruIndex, 0,
        &x_interval_beg, &x_interval_end);
        /********************x*********************/

        /********************y*********************/
        y_interval_beg = 0;
        y_interval_end = yReads->n; y_interval_end -= 1;
        /********************y*********************/
    }
    else if(type == YCX)
    {
        /********************x*********************/  
        x_interval_beg = 0;
        x_interval_end = xReads->n; x_interval_end -= 1;
        /********************x*********************/  

        /********************y*********************/
        bi_direction_hap_alignment_extention(yReads, yLeftBeg, yLeftLen, yRightBeg, yRightLen,
        xUid, 0, xReads->n - 1, Hap_rate, position_index, reverse_sources, read_g, ruIndex, rev,
        &y_interval_beg, &y_interval_end);
        /********************y*********************/
    }

    (*r_x_interval_beg) = x_interval_beg;
    (*r_x_interval_end) = x_interval_end;
    (*r_y_interval_beg) = y_interval_beg;
    (*r_y_interval_end) = y_interval_end;    
}

uint32_t vote_overlap_type(kvec_asg_arc_t_offset* u_buffer, hap_candidates* hap_can, 
uint64_t* position_index, ma_utg_t* xReads, ma_utg_t* yReads)
{
    uint32_t i, xBasePos, yBasePos;
    asg_arc_t_offset* arch = NULL;
    uint32_t flag[4];
    flag[X2Y] = flag[Y2X] = flag[XCY] = flag[YCX] = 0;


    for (i = hap_can->index_beg; i <= hap_can->index_end; i++)
    {
        arch = &(u_buffer->a.a[i]);
        xBasePos = (uint32_t)(arch->Off>>32);
        yBasePos = (uint32_t)(arch->Off);
        flag[classify_hap_overlap(xBasePos, xBasePos, xReads->len, yBasePos, yBasePos, yReads->len,
        NULL, NULL, NULL, NULL)]++;
    }
    
    uint32_t max_flag_i = 0;
    for (i = 0; i < 4; i++)
    {
        if(i == max_flag_i) continue;
        if(flag[i] > flag[max_flag_i])
        {
            max_flag_i = i;
        }
    }
    
    return max_flag_i;
}


void get_base_boundary(R_to_U* ruIndex, ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, 
asg_t *read_g, uint64_t* position_index, int max_hang, int min_ovlp, ma_utg_t *xReads, ma_utg_t *yReads,
uint32_t xUid, uint32_t yUid, long long begIndex, long long endIndex, uint32_t dir, uint32_t rev,
uint32_t* x_off, uint32_t* y_off)
{
    long long k, j;
    ma_hit_t_alloc *xR = NULL;
    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL, *st = NULL;
    int32_t r;
    asg_arc_t t;
    uint32_t rId, Hap_uId, is_Unitig, v, w, v_dir, w_dir, is_found = 0, oLen = 0;
    uint64_t tmp;
    (*x_off) = (*y_off) = (uint32_t)-1;
    if(dir == 1)
    {
        for (k = endIndex; k >= begIndex; k--)
        {
            xR = &(reverse_sources[xReads->a[k]>>33]);
            is_found = 0; oLen = 0;
            for (j = 0; j < xR->length; j++)
            {
                h = &(xR->buffer[j]);
                sq = &(coverage_cut[Get_qn(*h)]);
                st = &(coverage_cut[Get_tn(*h)]);
                if(st->del || read_g->seq[Get_tn(*h)].del) continue;

                r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, &t);
                ///if it is a contained overlap, skip
                if(r < 0) continue;

                rId = t.v>>1;
                if(read_g->seq[rId].del == 1) continue;
                ///there are two cases: 
                ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
                ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
                get_R_to_U(ruIndex, rId, &Hap_uId, &is_Unitig);
                if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
                if(Hap_uId != yUid) continue;

                v = xReads->a[k]>>32;
                get_R_to_U(ruIndex, v>>1, &Hap_uId, &is_Unitig);
                if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
                if(Hap_uId != xUid) continue;
                if((uint32_t)(position_index[v>>1]) != k) continue;

                w = (yReads->a[(uint32_t)(position_index[rId])])>>32;

                v_dir = ((t.ul>>32)==v)?1:0;
                w_dir = (t.v == w)?1:0;
                if(rev == 0 && v_dir != w_dir) continue;
                if(rev == 1 && v_dir == w_dir) continue;
                
                

                tmp = get_xy_pos(read_g, &t, v, w, xReads->len, yReads->len, position_index, &(t.el));
                ///if(is_found == 0 || ((uint32_t)(tmp>>32) > (*x_off) && ((uint32_t)tmp) > (*y_off)))
                if(is_found == 0 || t.ol > oLen)
                {
                    (*x_off) = tmp>>32;
                    (*y_off) = (uint32_t)tmp;
                    oLen = t.ol;
                } 

                is_found = 1; 
            }
            if(is_found) return;
        }   
    }
    else
    {
        for (k = begIndex; k <= endIndex; k++)
        {
            xR = &(reverse_sources[xReads->a[k]>>33]);
            is_found = 0; oLen = 0;
            for (j = 0; j < xR->length; j++)
            {
                h = &(xR->buffer[j]);
                sq = &(coverage_cut[Get_qn(*h)]);
                st = &(coverage_cut[Get_tn(*h)]);
                if(st->del || read_g->seq[Get_tn(*h)].del) continue;

                r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, &t);
                ///if it is a contained overlap, skip
                if(r < 0) continue;

                rId = t.v>>1;
                if(read_g->seq[rId].del == 1) continue;
                ///there are two cases: 
                ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
                ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
                get_R_to_U(ruIndex, rId, &Hap_uId, &is_Unitig);
                if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
                if(Hap_uId != yUid) continue;

                v = xReads->a[k]>>32;
                get_R_to_U(ruIndex, v>>1, &Hap_uId, &is_Unitig);
                if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
                if(Hap_uId != xUid) continue;
                if((uint32_t)(position_index[v>>1]) != k) continue;

                w = (yReads->a[(uint32_t)(position_index[rId])])>>32;

                v_dir = ((t.ul>>32)==v)?1:0;
                w_dir = (t.v == w)?1:0;
                if(rev == 0 && v_dir != w_dir) continue;
                if(rev == 1 && v_dir == w_dir) continue;
                
                

                tmp = get_xy_pos(read_g, &t, v, w, xReads->len, yReads->len, position_index, &(t.el));
                ///if(is_found == 0 || ((uint32_t)(tmp>>32) < (*x_off) && ((uint32_t)tmp) < (*y_off)))
                if(is_found == 0 || t.ol > oLen)
                {
                    (*x_off) = tmp>>32;
                    (*y_off) = (uint32_t)tmp;
                    oLen = t.ol;
                }     

                is_found = 1;         
            }
            if(is_found) return;
        }
    }


    (*x_off) = (*y_off) = (uint32_t)-1;
    
}

uint32_t determine_hap_overlap_type(hap_candidates* hap_can, ma_utg_t *xReads, ma_utg_t *yReads,
R_to_U* ruIndex, ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, asg_t *read_g, 
uint64_t* position_index, int max_hang, int min_ovlp, uint32_t xUid, uint32_t yUid,
long long* r_x_pos_beg, long long* r_x_pos_end, long long* r_y_pos_beg, long long* r_y_pos_end)
{
    uint32_t x_pos_beg, y_pos_beg, x_pos_end, y_pos_end;
    /*************************x***************************/
    get_base_boundary(ruIndex, reverse_sources, coverage_cut, read_g, position_index, max_hang, 
    min_ovlp, xReads, yReads, xUid, yUid, Get_x_beg(*hap_can), Get_x_end(*hap_can), 0, 
    Get_rev(*hap_can), &x_pos_beg, &y_pos_beg);
    get_base_boundary(ruIndex, reverse_sources, coverage_cut, read_g, position_index, max_hang, 
    min_ovlp, xReads, yReads, xUid, yUid, Get_x_beg(*hap_can), Get_x_end(*hap_can), 1, 
    Get_rev(*hap_can), &x_pos_end, &y_pos_end);
    /*************************x***************************/

    if(x_pos_beg == (uint32_t)-1 || y_pos_beg == (uint32_t)-1 
                    || x_pos_end == (uint32_t)-1 || y_pos_end == (uint32_t)-1)
    {
        return (uint32_t)-1;
    }
    if(x_pos_beg > x_pos_end || y_pos_beg > y_pos_end) return (uint32_t)-1;
    /**
    fprintf(stderr, "\nrev: %u, weight: %lu\n", Get_rev(*hap_can), Get_match(*hap_can));
        
    fprintf(stderr, "xUid: %u, xLen: %u, xBase: %u, x_interval_beg: %u, x_interval_end: %u, x_pos_beg: %u, x_pos_end: %u\n", 
    xUid, xReads->n, xReads->len, Get_x_beg(*hap_can), Get_x_end(*hap_can), x_pos_beg, x_pos_end);

    fprintf(stderr, "yUid: %u, yLen: %u, yBase: %u, y_interval_beg: %u, y_interval_end: %u, y_pos_beg: %u, y_pos_end: %u\n", 
    yUid, yReads->n, yReads->len, Get_y_beg(*hap_can), Get_y_end(*hap_can), y_pos_beg, y_pos_end);
    **/
    /**
    #define X2Y 0
    #define Y2X 1
    #define XCY 2
    #define YCX 3
    **/
    return classify_hap_overlap(x_pos_beg, x_pos_end, xReads->len, y_pos_beg, y_pos_end, yReads->len, 
    r_x_pos_beg, r_x_pos_end, r_y_pos_beg, r_y_pos_end);
}


uint32_t calculate_pair_hap_similarity(kvec_asg_arc_t_offset* u_buffer, hap_candidates* hap_can, 
uint64_t* position_index, uint32_t xUid, uint32_t yUid, ma_utg_t* xReads, ma_utg_t* yReads, 
ma_hit_t_alloc* reverse_sources, asg_t *read_g, R_to_U* ruIndex, ma_sub_t *coverage_cut, 
float Hap_rate, int max_hang, int min_ovlp, long long* r_x_pos_beg, long long* r_x_pos_end, 
long long* r_y_pos_beg, long long* r_y_pos_end)
{
    uint32_t max_count = 0, min_count = 0, i, flag;
    uint32_t xLen = xReads->n, xIndex/**, xBasePos**/;
    uint32_t yLen = yReads->n, yIndex/**, yBasePos**/;
    uint32_t xLeftBeg, xLeftLen, yLeftBeg, yLeftLen;
    uint32_t xRightBeg, xRightLen, yRightBeg, yRightLen;
    uint64_t totalWeigth;
    double xLeftMatch, xLeftTotal, yLeftMatch, yLeftTotal;
    double xRightMatch, xRightTotal, yRightMatch, yRightTotal;
    asg_arc_t_offset* arch = NULL;

    for (i = hap_can->index_beg, totalWeigth = 0; i <= hap_can->index_end; i++)
    {
        totalWeigth += u_buffer->a.a[i].weight;
        if(totalWeigth >= (hap_can->weight/2)) break;
    }
    if(i > hap_can->index_end) i = hap_can->index_end;

    arch = &(u_buffer->a.a[i]);
    xIndex = (uint32_t)(position_index[arch->x.ul>>33]);
    yIndex = (uint32_t)(position_index[arch->x.v>>1]);
    ///xBasePos = (uint32_t)(arch->Off>>32);
    ///yBasePos = (uint32_t)(arch->Off);
    
    if(hap_can->rev == 0)
    {
        xLeftBeg = 0; xLeftLen = xIndex; xRightBeg = xIndex; xRightLen = xLen - xRightBeg;
        yLeftBeg = 0; yLeftLen = yIndex; yRightBeg = yIndex; yRightLen = yLen - yRightBeg;
    }
    else
    {
        xLeftBeg = 0; xLeftLen = xIndex; xRightBeg = xIndex; xRightLen = xLen - xRightBeg;

        yLeftBeg = yIndex + 1; yLeftLen = yLen - yLeftBeg;
        yRightBeg = 0; yRightLen = yIndex + 1;
    }

    ///flag = classify_hap_overlap(xBasePos, xBasePos, xReads->len, yBasePos, yBasePos, yReads->len);
    flag = vote_overlap_type(u_buffer, hap_can, position_index, xReads, yReads);

    if(flag == XCY)
    {
        get_pair_hap_similarity(yReads->a, yLen, xUid, reverse_sources, read_g, ruIndex, 
        &yLeftMatch, &yLeftTotal);
        max_count = yLeftMatch;
        min_count = yLeftTotal;
    }
    else if(flag == YCX)
    {
        get_pair_hap_similarity(xReads->a, xLen, yUid, reverse_sources, read_g, ruIndex, 
        &xLeftMatch, &xLeftTotal);
        max_count = xLeftMatch;
        min_count = xLeftTotal;
    }
    else if(flag == X2Y)
    {
        get_pair_hap_similarity(yReads->a+yLeftBeg, yLeftLen, xUid, reverse_sources, read_g, ruIndex,
        &yLeftMatch, &yLeftTotal);
        get_pair_hap_similarity(xReads->a+xRightBeg, xRightLen, yUid, reverse_sources, read_g, ruIndex, 
        &xRightMatch, &xRightTotal);
        max_count = yLeftMatch + xRightMatch;
        min_count = yLeftTotal + xRightTotal;
    }
    else if(flag == Y2X)
    {
        get_pair_hap_similarity(xReads->a+xLeftBeg, xLeftLen, yUid, reverse_sources, read_g, ruIndex, 
        &xLeftMatch, &xLeftTotal);
        get_pair_hap_similarity(yReads->a+yRightBeg, yRightLen, xUid, reverse_sources, read_g, ruIndex,
        &yRightMatch, &yRightTotal);
        max_count = xLeftMatch + yRightMatch;
        min_count = xLeftTotal + yRightTotal;
    }
    

    hap_can->weight = hap_can->index_beg = 0;
    if(min_count == 0) return NON_PLOID;
    if(max_count > min_count*Hap_rate)
    {
        long long r_x_interval_beg, r_x_interval_end, r_y_interval_beg, r_y_interval_end;
        ///for containment, don't need to do anything
        get_hap_alignment_boundary(xReads, yReads, flag, xLeftMatch, xLeftTotal, 
        yLeftMatch, yLeftTotal, xRightMatch, xRightTotal, yRightMatch, yRightTotal,
        xLeftBeg, xLeftLen, yLeftBeg, yLeftLen, xRightBeg, xRightLen, yRightBeg, yRightLen,
        xUid, yUid, Hap_rate, position_index, reverse_sources, read_g, ruIndex, hap_can->rev,
        &r_x_interval_beg, &r_x_interval_end, &r_y_interval_beg, &r_y_interval_end);

        if(r_x_interval_beg < 0 || r_x_interval_end < 0 || r_y_interval_beg < 0 || r_y_interval_end < 0)
        {
            return NON_PLOID;
        }

        get_pair_hap_similarity(xReads->a + r_x_interval_beg, r_x_interval_end + 1 - r_x_interval_beg, 
        yUid, reverse_sources, read_g, ruIndex, &xLeftMatch, &xLeftTotal);
        if(xLeftMatch == 0 || xLeftTotal == 0) return NON_PLOID;
        
        hap_can->weight = xLeftMatch;
        hap_can->index_beg = xLeftTotal;
        hap_can->index_end = flag;
        hap_can->x_beg_pos = r_x_interval_beg;
        hap_can->x_end_pos = r_x_interval_end;
        hap_can->y_beg_pos = r_y_interval_beg;
        hap_can->y_end_pos = r_y_interval_end;

        hap_can->index_end = determine_hap_overlap_type(hap_can, xReads, yReads, ruIndex, 
        reverse_sources, coverage_cut, read_g, position_index, max_hang, min_ovlp, xUid, 
        yUid, r_x_pos_beg, r_x_pos_end, r_y_pos_beg, r_y_pos_end);
        if(hap_can->index_end == (uint32_t)-1) return NON_PLOID;

        return PLOID;
    } 
	return NON_PLOID;
}


void hap_alignment(ma_ug_t *ug, asg_t *read_g,  ma_hit_t_alloc* reverse_sources, 
R_to_U* ruIndex, ma_sub_t *coverage_cut, uint64_t* position_index, uint64_t* vote_counting, 
uint8_t* visit, kvec_t_u64_warp* u_vecs, kvec_asg_arc_t_offset* u_buffer, kvec_hap_candidates* u_can,
uint32_t Input_uId, float Hap_rate, int max_hang, int min_ovlp, float chain_rate, 
hap_overlaps_list* all_ovlp)
{
    ma_utg_t *xReads = NULL, *yReads = NULL;
    ma_hit_t_alloc *xR = NULL;
    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL, *st = NULL;
    asg_t* nsg = ug->g;
    uint32_t i, j, v, rId, k, is_Unitig, Hap_uId, xUid, yUid, seedOcc, xPos, yPos, is_update;
    uint64_t tmp;
    long long cur_offset, new_offset, interval_len;
    long long r_x_pos_beg, r_x_pos_end, r_y_pos_beg, r_y_pos_end;
    int32_t r;
    asg_arc_t t;
    asg_arc_t_offset t_offset;
    hap_candidates hap_can;
    hap_overlaps hap_align;
    xUid = Input_uId;
    if(nsg->seq[xUid].del || nsg->seq[xUid].c == ALTER_LABLE) return;
    memset(vote_counting, 0, sizeof(uint64_t)*nsg->n_seq);
    memset(visit, 0, nsg->n_seq);
    u_vecs->a.n = 0;
    u_can->a.n = 0;

    xReads = &(ug->u.a[xUid]);
    for (i = 0; i < xReads->n; i++)
    {
        xR = &(reverse_sources[xReads->a[i]>>33]);
        
        for (k = 0; k < xR->length; k++)
        {
            rId = Get_tn(xR->buffer[k]);
        
            if(read_g->seq[rId].del == 1)
            {
                ///get the id of read that contains it 
                get_R_to_U(ruIndex, rId, &rId, &is_Unitig);
                if(rId == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[rId].del == 1) continue;
            }

            ///there are two cases: 
            ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
            ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
            get_R_to_U(ruIndex, rId, &Hap_uId, &is_Unitig);
            if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
            ///here rId is the id of the read coming from the different haplotype
            ///Hap_cId is the id of the corresponding contig (note here is the contig, instead of untig)
            if(visit[Hap_uId]!=0) continue;
            visit[Hap_uId] = 1;
            if(vote_counting[Hap_uId] < UINT64_MAX) vote_counting[Hap_uId]++;
        }

        clean_visit_flag(visit, read_g, ruIndex, nsg->n_seq, xR);
    }



    u_vecs->a.n = 0;
    for (i = 0; i < nsg->n_seq; i++)
    {
        if(i == xUid) continue;
        if(vote_counting[i] == 0) continue;
        tmp = vote_counting[i]; tmp = tmp << 32; tmp = tmp | (uint64_t)i;
        kv_push(uint64_t, u_vecs->a, tmp);
    }

    if(u_vecs->a.n == 0) return;
    sort_kvec_t_u64_warp(u_vecs, 1);


    ///scan each candidate unitig
    for (i = 0; i < u_vecs->a.n; i++)
    {
        yUid = (uint32_t)u_vecs->a.a[i];
        seedOcc = u_vecs->a.a[i]>>32;
        xReads = &(ug->u.a[xUid]);
        yReads = &(ug->u.a[yUid]);
        u_buffer->a.n = 0;
        

        //if(xUid == 404 && yUid == 307)
        // if(xReads->n <= 100 && yReads->n <= 100)
        // {
        //     debug_enable = 1;
        // }
        // else
        // {
        //     debug_enable = 0;
        // }
        // debug_enable = 1;

        for (k = 0; k < xReads->n; k++)
        {
            xR = &(reverse_sources[xReads->a[k]>>33]);
            for (j = 0; j < xR->length; j++)
            {
                h = &(xR->buffer[j]);
                sq = &(coverage_cut[Get_qn(*h)]);
                st = &(coverage_cut[Get_tn(*h)]);
                if(st->del || read_g->seq[Get_tn(*h)].del) continue;

                r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, &t);
                ///if it is a contained overlap, skip
                if(r < 0) continue;

                rId = t.v>>1;
                if(read_g->seq[rId].del == 1) continue;
                ///there are two cases: 
                ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
                ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
                get_R_to_U(ruIndex, rId, &Hap_uId, &is_Unitig);
                if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
                if(Hap_uId != yUid) continue;

                v = xReads->a[k]>>32;
                get_R_to_U(ruIndex, v>>1, &Hap_uId, &is_Unitig);
                if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
                if(Hap_uId != xUid) continue;
                if((uint32_t)(position_index[v>>1]) != k) continue;


                if((prefilter((uint32_t)(position_index[v>>1]), (uint32_t)(position_index[rId]), 
                    xReads->n, yReads->n, 0, Hap_rate, seedOcc)==NON_PLOID) && 
                   (prefilter((uint32_t)(position_index[v>>1]), (uint32_t)(position_index[rId]), 
                    xReads->n, yReads->n, 1, Hap_rate, seedOcc)==NON_PLOID))
                {
                    continue;
                }
            
                t_offset.Off = get_xy_pos(read_g, &t, v, (yReads->a[(uint32_t)(position_index[rId])])>>32, 
                xReads->len, yReads->len, position_index, &(t.el));
                t_offset.x = t;
                t_offset.weight = 1;

                kv_push(asg_arc_t_offset, u_buffer->a, t_offset);
                // if(debug_enable)
                // {
                //     fprintf(stderr, "xUid: %u, yUid: %u, (%u), x_index: %u, y_index: %u, dis: %lld\n\n", 
                //     xUid, yUid, (uint32_t)(u_buffer->a.n-1), 
                //     (uint32_t)(position_index[u_buffer->a.a[u_buffer->a.n-1].x.ul>>33]),
                //     (uint32_t)(position_index[u_buffer->a.a[u_buffer->a.n-1].x.v>>1]),
                //     (long long)((uint32_t)(u_buffer->a.a[u_buffer->a.n-1].Off>>32)) - 
                //     (long long)((uint32_t)(u_buffer->a.a[u_buffer->a.n-1].Off)));
                // }    
            }

            deduplicate_edge(u_buffer);
        }

        if(u_buffer->a.n == 0) continue;

        // if(debug_enable)
        // {
        //     print_debug_unitig(xReads, position_index, "xReads");
        //     print_debug_unitig(yReads, position_index, "yReads");
        // }
        
        qsort(u_buffer->a.a, u_buffer->a.n, sizeof(asg_arc_t_offset), cmp_hap_alignment);
        k = 0;
        u_can->a.n = 0;
        while (k < u_buffer->a.n)
        {
            hap_can.rev = u_buffer->a.a[k].x.el;
            hap_can.index_beg = k;
            hap_can.index_end = k;
            hap_can.weight = u_buffer->a.a[k].weight;
            hap_can.x_beg_pos = hap_can.x_end_pos = (uint32_t)(u_buffer->a.a[k].Off>>32);
            hap_can.y_beg_pos = hap_can.y_end_pos = (uint32_t)(u_buffer->a.a[k].Off);
            cur_offset = Cal_Off(u_buffer->a.a[k].Off);
            interval_len = get_hap_overlapLen(hap_can.x_beg_pos, hap_can.x_end_pos, xReads->len, 
            hap_can.y_beg_pos, hap_can.y_end_pos, yReads->len, NULL, NULL, NULL, NULL);


            k++;
            while (k < u_buffer->a.n)
            {
                new_offset = Cal_Off(u_buffer->a.a[k].Off);
                if(u_buffer->a.a[k].x.el != hap_can.rev) break;
                if((new_offset - cur_offset)>(interval_len*chain_rate)) break;
                

                hap_can.index_end = k;
                hap_can.weight += u_buffer->a.a[k].weight;

                is_update = 0;
                xPos = (uint32_t)(u_buffer->a.a[k].Off>>32);
                yPos = (uint32_t)(u_buffer->a.a[k].Off);
                if(xPos < hap_can.x_beg_pos)
                {
                    hap_can.x_beg_pos = xPos;
                    is_update = 1;
                } 

                if(xPos > hap_can.x_end_pos)
                {
                    hap_can.x_end_pos = xPos;
                    is_update = 1;
                }

                if(yPos < hap_can.y_beg_pos)
                {
                    hap_can.y_beg_pos = yPos;
                    is_update = 1;
                }

                if(yPos > hap_can.y_end_pos)
                {
                    hap_can.y_end_pos = yPos;
                    is_update = 1;
                }

                if(new_offset == cur_offset) is_update = 0;

                if(is_update)
                {
                    interval_len = get_hap_overlapLen(hap_can.x_beg_pos, hap_can.x_end_pos, xReads->len, 
                    hap_can.y_beg_pos, hap_can.y_end_pos, yReads->len, NULL, NULL, NULL, NULL);
                }

                k++;
            }

            kv_push(hap_candidates, u_can->a, hap_can);
        }

        if(u_can->a.n == 0) continue;
        
        qsort(u_can->a.a, u_can->a.n, sizeof(hap_candidates), cmp_hap_candidates);

        Get_match(hap_can) = Get_total(hap_can) = 0;
        memset(&hap_align, 0, sizeof(hap_overlaps));
        
        for (k = 0; k < u_can->a.n; k++)
        {
            is_update = 0;
            if(u_can->a.a[k].weight < Get_match(hap_can)*Hap_rate) continue;

            if(calculate_pair_hap_similarity(u_buffer, &(u_can->a.a[k]), position_index, xUid, yUid, 
            xReads, yReads, reverse_sources, read_g, ruIndex, coverage_cut, Hap_rate, max_hang, 
            min_ovlp, &r_x_pos_beg, &r_x_pos_end, &r_y_pos_beg, &r_y_pos_end)!=PLOID)
            {
                continue;
            }
            
            if(Get_match(hap_can) < Get_match(u_can->a.a[k]))
            {
                is_update = 1;
            }
            else if(Get_match(hap_can) == Get_match(u_can->a.a[k]) && 
                                    Get_total(hap_can) > Get_total(u_can->a.a[k]))
            {
                is_update = 1;
            }

            if(is_update)
            {
                hap_can = u_can->a.a[k];
                hap_align.rev = Get_rev(hap_can);
                hap_align.type = Get_type(hap_can);
                hap_align.x_beg_id = Get_x_beg(hap_can);
                hap_align.x_end_id = Get_x_end(hap_can) + 1;
                hap_align.y_beg_id = Get_y_beg(hap_can);
                hap_align.y_end_id = Get_y_end(hap_can) + 1;
                hap_align.weight = Get_match(hap_can);
                hap_align.x_beg_pos = r_x_pos_beg;
                hap_align.x_end_pos = r_x_pos_end + 1;
                if(hap_align.rev == 0)
                {
                    hap_align.y_beg_pos = r_y_pos_beg;
                    hap_align.y_end_pos = r_y_pos_end + 1;
                }
                else
                {
                    hap_align.y_beg_pos = yReads->len - r_y_pos_end - 1;
                    hap_align.y_end_pos = yReads->len - r_y_pos_beg - 1 + 1;
                }
                hap_align.xUid = xUid;
                hap_align.yUid = yUid;
                hap_align.status = SELF_EXIST;
            }
        }

        if(Get_match(hap_can) == 0 || Get_total(hap_can) == 0) continue;
        
        kv_push(hap_overlaps, all_ovlp->x[hap_align.xUid].a, hap_align);

        /**
        fprintf(stderr, "\nsplit\nxUid: %u, yUid: %u, u_buffer->a.n: %u\n", xUid, yUid, (uint32_t)u_buffer->a.n);

        for (k = 0; k < u_buffer->a.n; k++)
        {
            v = u_buffer->a.a[k].x.ul>>33;
            w = u_buffer->a.a[k].x.v>>1;

            fprintf(stderr, "\n(%u), x_index: %u, x_pos: %u, x_real_pos: %u, dis: %lld, dir: %u, weight: %lu\n", 
            k, (uint32_t)(position_index[v]), (uint32_t)(position_index[v]>>32), 
            (uint32_t)(u_buffer->a.a[k].Off>>32), Cal_Off(u_buffer->a.a[k].Off),
            u_buffer->a.a[k].x.el, u_buffer->a.a[k].weight);

            fprintf(stderr, "(%u), y_index: %u, y_pos: %u, y_real_pos: %u\n", k,
            (uint32_t)(position_index[w]), (uint32_t)(position_index[w]>>32), 
            (uint32_t)(u_buffer->a.a[k].Off));
        }
        **/
    }

}

void print_hap_paf(ma_ug_t *ug, hap_overlaps* ovlp)
{
    fprintf(stderr, "utg%.6d%c\t%u\t%u\t%u\t%c\tutg%.6d%c\t%u\t%u\t%u\t%u\t%u\n", 
    ovlp->xUid+1, "lc"[ug->u.a[ovlp->xUid].circ], ug->u.a[ovlp->xUid].len, ovlp->x_beg_pos, ovlp->x_end_pos, "+-"[ovlp->rev], 
    ovlp->yUid+1, "lc"[ug->u.a[ovlp->yUid].circ], ug->u.a[ovlp->yUid].len, ovlp->y_beg_pos, ovlp->y_end_pos, ovlp->type, (uint32_t)ovlp->weight);
}


int inline get_specific_hap_overlap(kvec_hap_overlaps* x, uint32_t qn, uint32_t tn)
{
    uint32_t i;
    for (i = 0; i < x->a.n; i++)
    {
        if(x->a.a[i].xUid == qn && x->a.a[i].yUid == tn)
        {
            return i;
        }        
    }

    return -1;
}

void set_reverse_hap_overlap(hap_overlaps* dest, hap_overlaps* source, uint32_t* types)
{
    dest->status = REVE_EXIST;
    dest->rev = source->rev;
    dest->type = types[source->type];
    dest->weight = source->weight;
    dest->xUid = source->yUid;
    dest->yUid = source->xUid;
    dest->x_beg_pos = source->y_beg_pos;
    dest->x_end_pos = source->y_end_pos;
    dest->y_beg_pos = source->x_beg_pos;
    dest->y_end_pos = source->x_end_pos;
    dest->x_beg_id = source->y_beg_id;
    dest->x_end_id = source->y_end_id;
    dest->y_beg_id = source->x_beg_id;
    dest->y_end_id = source->x_end_id;
}

/**
#define X2Y 0
#define Y2X 1
#define XCY 2
#define YCX 3
**/
void normalize_hap_overlaps(hap_overlaps_list* all_ovlp, hap_overlaps_list* back_all_ovlp)
{
    hap_overlaps *x = NULL, *y = NULL;
    uint32_t v, i, uId, qn, tn;
    uint32_t types[4];
    types[X2Y] = Y2X; types[Y2X] = X2Y; types[XCY] = YCX; types[YCX] = XCY;
    int index;
    for (v = 0; v < all_ovlp->num; v++)
    {
        uId = v;
        for (i = 0; i < all_ovlp->x[uId].a.n; i++)
        {
            qn = all_ovlp->x[uId].a.a[i].xUid;
            tn = all_ovlp->x[uId].a.a[i].yUid;
            x = &(all_ovlp->x[uId].a.a[i]);
            index = get_specific_hap_overlap(&(all_ovlp->x[tn]), tn, qn);
            if(index != -1)
            {
                y = &(all_ovlp->x[tn].a.a[index]);
                if(x->rev == y->rev && types[x->type]==y->type) continue;
                if(x->weight >= y->weight)
                {
                    kv_push(hap_overlaps, back_all_ovlp->x[tn].a, (*y));
                    set_reverse_hap_overlap(y, x, types);
                }
                else
                {
                    kv_push(hap_overlaps, back_all_ovlp->x[qn].a, (*x));
                    set_reverse_hap_overlap(x, y, types);
                }
            }
            else
            {
                kv_pushp(hap_overlaps, all_ovlp->x[tn].a, &y);
                set_reverse_hap_overlap(y, x, types);
            }
        }
    }
}


void debug_hap_overlaps(hap_overlaps_list* all_ovlp, hap_overlaps_list* back_all_ovlp)
{
    hap_overlaps *x = NULL, *y = NULL;
    uint32_t v, i, uId, qn, tn;
    uint32_t types[4];
    types[X2Y] = Y2X; types[Y2X] = X2Y; types[XCY] = YCX; types[YCX] = XCY;
    int index;
    for (v = 0; v < all_ovlp->num; v++)
    {
        uId = v;
        for (i = 0; i < all_ovlp->x[uId].a.n; i++)
        {
            qn = all_ovlp->x[uId].a.a[i].xUid;
            tn = all_ovlp->x[uId].a.a[i].yUid;
            x = &(all_ovlp->x[uId].a.a[i]);
            index = get_specific_hap_overlap(&(all_ovlp->x[tn]), tn, qn);
            if(index == -1)
            {
                fprintf(stderr, "ERROR 0\n");
                continue;
            }

            y = &(all_ovlp->x[tn].a.a[index]);
            if(x->rev != y->rev || types[x->type] != y->type)
            {
                fprintf(stderr, "ERROR 1\n");
                continue;
            }

            if(x->status == REVE_EXIST && y->status != SELF_EXIST)
            {
                fprintf(stderr, "ERROR 2\n");
                continue;
            }

            if(x->status == REVE_EXIST)
            {
                if(x->weight != y->weight) fprintf(stderr, "ERROR 3\n"); 
                if(x->xUid != y->yUid) fprintf(stderr, "ERROR 4\n"); 
                if(x->yUid != y->xUid) fprintf(stderr, "ERROR 5\n"); 
                if(x->x_beg_pos != y->y_beg_pos) fprintf(stderr, "ERROR 6\n"); 
                if(x->x_end_pos != y->y_end_pos) fprintf(stderr, "ERROR 7\n"); 
                if(x->y_beg_pos != y->x_beg_pos) fprintf(stderr, "ERROR 8\n");
                if(x->y_end_pos != y->x_end_pos) fprintf(stderr, "ERROR 9\n");
                if(x->x_beg_id != y->y_beg_id) fprintf(stderr, "ERROR 10\n");
                if(x->x_end_id != y->y_end_id) fprintf(stderr, "ERROR 11\n");
                if(x->y_beg_id != y->x_beg_id) fprintf(stderr, "ERROR 12\n");
                if(x->y_beg_id != y->x_beg_id) fprintf(stderr, "ERROR 13\n");
                if(x->y_end_id != y->x_end_id) fprintf(stderr, "ERROR 14\n");

                index = get_specific_hap_overlap(&(back_all_ovlp->x[qn]), qn, tn);
                if(index != -1)
                {
                    if(back_all_ovlp->x[qn].a.a[index].weight > x->weight) fprintf(stderr, "ERROR 15\n");
                }
            }
           
        }
    }
}

void print_purge_gfa(ma_ug_t *ug, asg_t *purge_g)
{
    uint32_t v, i, n_vtx = purge_g->n_seq * 2;
    for (v = 0; v < n_vtx; v++)
    {
        if(v%2==0) fprintf(stderr, "\n");
        if(purge_g->seq[v>>1].del)
        {
            fprintf(stderr, "(D) v>>1: %u, v&1: %u, utg%.6d%c\n", v>>1, v&1, (v>>1)+1,
            "lc"[ug->u.a[v>>1].circ]);
            continue;
        } 

        fprintf(stderr, "(E) v>>1: %u, v&1: %u, utg%.6dl%c\n", v>>1, v&1, (v>>1)+1,
        "lc"[ug->u.a[v>>1].circ]);
        
        uint32_t nv = asg_arc_n(purge_g, v);
		asg_arc_t *av = asg_arc_a(purge_g, v);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            fprintf(stderr, "av[i].ul: %u (utg%.6d%c, dir: %u, len: %u), av[i].v: %u (utg%.6d%c, dir: %u, len: %u), ol: %u\n", 
                (uint32_t)(av[i].ul>>33), (uint32_t)(av[i].ul>>33)+1, "lc"[ug->u.a[av[i].ul>>33].circ], (uint32_t)(av[i].ul>>32)&1, ug->u.a[av[i].ul>>33].len, 
                av[i].v>>1, (av[i].v>>1)+1, "lc"[ug->u.a[av[i].v>>1].circ], av[i].v&1, ug->u.a[av[i].v>>1].len, av[i].ol);
        }

    }
    
}



// pop bubbles
int asg_pop_bubble_purge_graph(asg_t *purge_g, int max_dist)
{
	uint32_t v, n_vtx = purge_g->n_seq * 2;
	uint64_t n_pop = 0;
	buf_t b;
	if (!purge_g->is_symm) asg_symm(purge_g);
	memset(&b, 0, sizeof(buf_t));
	///set information for each node
	b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
	//traverse all node with two directions 
	for (v = 0; v < n_vtx; ++v) {
		uint32_t i, n_arc = 0, nv = asg_arc_n(purge_g, v);
		asg_arc_t *av = asg_arc_a(purge_g, v);
		///some node could be deleted
		if (nv < 2 || purge_g->seq[v>>1].del || purge_g->seq[v>>1].c == ALTER_LABLE) continue;
		///some edges could be deleted
		for (i = 0; i < nv; ++i) // asg_bub_pop1() may delete some edges/arcs
			if (!av[i].del) ++n_arc;
		if (n_arc > 1)
			n_pop += asg_bub_pop1_primary_trio(purge_g, NULL, v, max_dist, &b, (uint32_t)-1, DROP, 1);
	}
	free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
	if (n_pop) asg_cleanup(purge_g);

	if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] popped %lu bubbles\n", __func__, (unsigned long)n_pop);
    }
    return n_pop;
}

int get_hap_arch(hap_overlaps* hap, uint32_t qLen, uint32_t tLen, int max_hang, float max_hang_rate, 
int min_ovlp, asg_arc_t* t)
{
    int r;
    ma_hit_t h;
    h.qns = hap->xUid;
    h.qns = h.qns << 32;
    h.qns = h.qns | hap->x_beg_pos;
    h.qe = hap->x_end_pos;
    h.tn = hap->yUid;
    h.ts = hap->y_beg_pos;
    h.te = hap->y_end_pos;
    h.rev = hap->rev;
    h.del = 0;
    h.bl = h.el = h.ml = h.no_l_indel = 0;

    r = ma_hit2arc(&h, qLen, tLen, max_hang, max_hang_rate, min_ovlp, t);
    return r;
}



void clean_purge_graph(asg_t *purge_g, int max_dist, float drop_ratio)
{
    uint64_t operation = 1;
    while (operation > 0)
    {
        operation = 0;
        operation += asg_pop_bubble_purge_graph(purge_g, max_dist);
        operation += unitig_arc_del_short_diploid_by_length(purge_g, drop_ratio);
    }
    
    unitig_arc_del_short_diploid_by_length(purge_g, 1);
}


void get_node_boundary(R_to_U* ruIndex, ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, 
asg_t *read_g, uint64_t* position_index, int max_hang, int min_ovlp, ma_utg_t *xReads, ma_utg_t *yReads,
uint32_t xUid, uint32_t yUid, long long begIndex, long long endIndex, uint32_t dir, uint32_t rev,
asg_arc_t* reture_t_f, asg_arc_t* reture_t_r)
{
    long long k, j;
    ma_hit_t_alloc *xR = NULL;
    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL, *st = NULL;
    int r, index;
    asg_arc_t t_f, t_r;
    uint32_t rId, Hap_uId, is_Unitig, v, w, v_dir, w_dir, is_found = 0, oLen = 0;
    reture_t_f->del = reture_t_r->del = 1;
    if(dir == 1)
    {
        for (k = endIndex; k >= begIndex; k--)
        {
            xR = &(reverse_sources[xReads->a[k]>>33]);
            is_found = 0;
            for (j = 0; j < xR->length; j++)
            {
                h = &(xR->buffer[j]);
                sq = &(coverage_cut[Get_qn(*h)]);
                st = &(coverage_cut[Get_tn(*h)]);
                if(st->del || read_g->seq[Get_tn(*h)].del) continue;
                r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, &t_f);
                ///if it is a contained overlap, skip
                if(r < 0) continue;

                rId = t_f.v>>1;
                if(read_g->seq[rId].del == 1) continue;
                ///there are two cases: 
                ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
                ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
                get_R_to_U(ruIndex, rId, &Hap_uId, &is_Unitig);
                if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
                if(Hap_uId != yUid) continue;

                v = xReads->a[k]>>32;
                get_R_to_U(ruIndex, v>>1, &Hap_uId, &is_Unitig);
                if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
                if(Hap_uId != xUid) continue;
                if((uint32_t)(position_index[v>>1]) != k) continue;

                w = (yReads->a[(uint32_t)(position_index[rId])])>>32;
                v_dir = ((t_f.ul>>32)==v)?1:0;
                w_dir = (t_f.v == w)?1:0;

                if(rev == 0 && v_dir != w_dir) continue;
                if(rev == 1 && v_dir == w_dir) continue;
                if(v_dir == 1) continue;
                /************************get reverse edge*************************/
                index = get_specific_overlap(&(reverse_sources[Get_tn(*h)]), Get_tn(*h), Get_qn(*h));
                if(index == -1) continue;
                h = &(reverse_sources[Get_tn(*h)].buffer[index]);
                sq = &(coverage_cut[Get_qn(*h)]);
                st = &(coverage_cut[Get_tn(*h)]);
                if(st->del || read_g->seq[Get_tn(*h)].del) continue;
                r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, &t_r);
                if(r < 0) continue;
                /************************get reverse edge*************************/

                if(is_found == 0 || t_f.ol > oLen)
                {
                    (*reture_t_f) = t_f;
                    (*reture_t_r) = t_r;
                    oLen = t_f.ol;
                } 

                is_found = 1; 
            }
            if(is_found) return;
        }
    }
    else
    {
        for (k = begIndex; k <= endIndex; k++)
        {
            xR = &(reverse_sources[xReads->a[k]>>33]);
            is_found = 0; oLen = 0;
            for (j = 0; j < xR->length; j++)
            {
                h = &(xR->buffer[j]);
                sq = &(coverage_cut[Get_qn(*h)]);
                st = &(coverage_cut[Get_tn(*h)]);
                if(st->del || read_g->seq[Get_tn(*h)].del) continue;

                r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, &t_f);
                ///if it is a contained overlap, skip
                if(r < 0) continue;

                rId = t_f.v>>1;
                if(read_g->seq[rId].del == 1) continue;
                ///there are two cases: 
                ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
                ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
                get_R_to_U(ruIndex, rId, &Hap_uId, &is_Unitig);
                if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
                if(Hap_uId != yUid) continue;

                v = xReads->a[k]>>32;
                get_R_to_U(ruIndex, v>>1, &Hap_uId, &is_Unitig);
                if(is_Unitig == 0 || Hap_uId == (uint32_t)-1) continue;
                if(Hap_uId != xUid) continue;
                if((uint32_t)(position_index[v>>1]) != k) continue;

                w = (yReads->a[(uint32_t)(position_index[rId])])>>32;

                v_dir = ((t_f.ul>>32)==v)?1:0;
                w_dir = (t_f.v == w)?1:0;
                if(rev == 0 && v_dir != w_dir) continue;
                if(rev == 1 && v_dir == w_dir) continue;
                if(v_dir == 0) continue;
            
                /************************get reverse edge*************************/
                index = get_specific_overlap(&(reverse_sources[Get_tn(*h)]), Get_tn(*h), Get_qn(*h));
                if(index == -1) continue;
                h = &(reverse_sources[Get_tn(*h)].buffer[index]);
                sq = &(coverage_cut[Get_qn(*h)]);
                st = &(coverage_cut[Get_tn(*h)]);
                if(st->del || read_g->seq[Get_tn(*h)].del) continue;
                r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, &t_r);
                if(r < 0) continue;
                /************************get reverse edge*************************/

                if(is_found == 0 || t_f.ol > oLen)
                {
                    (*reture_t_f) = t_f;
                    (*reture_t_r) = t_r;
                    oLen = t_f.ol;
                } 

                is_found = 1;         
            }
            if(is_found) return;
        }
    }
    
}

void fill_unitig(uint64_t* buffer, uint32_t bufferLen, asg_t* read_g, kvec_asg_arc_t_warp* edge,
uint32_t is_circle, uint64_t* rLen)
{
    uint32_t i, k, totalLen, v, w, nv, l;
    asg_arc_t *av  = NULL;
    (*rLen) = totalLen = 0;
    for (i = 0; i < bufferLen - 1; i++)
    {
        v = (uint64_t)(buffer[i])>>32;
        w = (uint64_t)(buffer[i + 1])>>32; 
        av = asg_arc_a(read_g, v);
        nv = asg_arc_n(read_g, v);
        l = 0;
        
        for (k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            if(av[k].v == w) 
            {
                l = asg_arc_len(av[k]);
                break;
            }
        }

        if(k == nv)
        {
            for (k = 0; k < edge->a.n; k++)
            {
                if(edge->a.a[k].del) continue;
                if((edge->a.a[k].ul>>32) == v && edge->a.a[k].v == w)
                {
                    l = asg_arc_len(edge->a.a[k]);
                    break;
                }
            }

            if(k == edge->a.n)
            {
                fprintf(stderr, "####ERROR1: i: %u, v>>1: %u, v&1: %u, w>>1: %u, w&1: %u\n",
                i, v>>1, v&1, w>>1, w&1);
            }
        } 


        
        buffer[i] = v; buffer[i] = buffer[i]<<32; buffer[i] = buffer[i] | (uint64_t)(l);
        totalLen += l;
    }
    
    if(i < bufferLen)
    {
        if(is_circle)
        {
            v = (uint64_t)(buffer[i])>>32;
            w = (uint64_t)(buffer[0])>>32; 
            av = asg_arc_a(read_g, v);
            nv = asg_arc_n(read_g, v);
            l = 0;
            
            for (k = 0; k < nv; k++)
            {
                if(av[k].del) continue;
                if(av[k].v == w) 
                {
                    l = asg_arc_len(av[k]);
                    break;
                }
            }

            if(k == nv)
            {
                for (k = 0; k < edge->a.n; k++)
                {
                    if(edge->a.a[k].del) continue;
                    if((edge->a.a[k].ul>>32) == v && edge->a.a[k].v == w)
                    {
                        l = asg_arc_len(edge->a.a[k]);
                        break;
                    }
                }

                if(k == edge->a.n)
                {
                    fprintf(stderr, "####ERROR2: i: %u, v>>1: %u, v&1: %u, w>>1: %u, w&1: %u\n",
                    i, v>>1, v&1, w>>1, w&1);
                }
            } 

            buffer[i] = v; buffer[i] = buffer[i]<<32; buffer[i] = buffer[i] | (uint64_t)(l);
            totalLen += l;
        }
        else
        {
            v = (uint64_t)(buffer[i])>>32;
            l = read_g->seq[v>>1].len;
            buffer[i] = v;
            buffer[i] = buffer[i]<<32;
            buffer[i] = buffer[i] | (uint64_t)(l);
            totalLen += l;
        }
    }

    (*rLen) = totalLen;

}

void purge_merge(asg_t *purge_g, ma_ug_t *ug, hap_overlaps_list* all_ovlp, buf_t* b_0,
R_to_U* ruIndex, ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, asg_t *read_g, 
uint64_t* position_index, int max_hang, int min_ovlp, kvec_asg_arc_t_warp* edge, uint8_t* visit)
{
    uint32_t i, nv, k, v, w, x_beg_index, x_end_index, y_beg_index, y_end_index, cut_beg, cut_end, begIndex, endIndex, keepUid;
    hap_overlaps *x = NULL/**, *y = NULL**/;
    ma_utg_t *xReads = NULL, *yReads = NULL;
    asg_arc_t t_forward, t_backward;
    asg_arc_t *av = NULL;
    kvec_t(uint64_t) buffer;
    uint64_t totalLen;
    int index = 0;
    i = 0;
    while (i < b_0->b.n)
    {
        cut_beg = 0; cut_end = (uint32_t)-1;
        kv_init(buffer);
        /********************for the first node********************/
        v = b_0->b.a[i];
        keepUid = v>>1;
        xReads = &(ug->u.a[v>>1]);
        if(v&1)
        {
            for (k = 0; k < xReads->n; k++)
            {
                ///aim[query->n - j - 1] = (query->a[j])^(uint64_t)(0x100000000);
                kv_push(uint64_t, buffer, (xReads->a[xReads->n - k - 1])^(uint64_t)(0x100000000));
            }
        }
        else
        {
            for (k = 0; k < xReads->n; k++)
            {
                kv_push(uint64_t, buffer, xReads->a[k]);
            }
        }
        cut_beg = 0; cut_end = xReads->n - 1;
        i++;
        /********************for the first node********************/

        
        for (; i < b_0->b.n; i++)
        {
            ///x = y = NULL;
            x = NULL;

            v = b_0->b.a[i-1];
            w = b_0->b.a[i];
            

            index = get_specific_hap_overlap(&(all_ovlp->x[v>>1]), v>>1, w>>1);
            x = &(all_ovlp->x[v>>1].a.a[index]);


        
            xReads = &(ug->u.a[v>>1]);
            yReads = &(ug->u.a[w>>1]);

            begIndex = x->x_beg_id;
            if(cut_beg > begIndex) begIndex = cut_beg;

            endIndex = x->x_end_id-1;
            if(cut_end < endIndex) endIndex = cut_end;

            get_node_boundary(ruIndex, reverse_sources, coverage_cut, read_g, position_index, max_hang, 
            min_ovlp, xReads, yReads, v>>1, w>>1, begIndex, endIndex, v&1, x->rev, &t_forward, &t_backward);
            if(t_forward.del || t_backward.del) break;
            
            kv_push(asg_arc_t, edge->a, t_forward);
            kv_push(asg_arc_t, edge->a, t_backward);

            x_beg_index = 0; x_end_index = xReads->n - 1;
            y_beg_index = 0; y_end_index = yReads->n - 1;
            
            if((v&1) == 0)
            {
                x_end_index = (uint32_t)position_index[t_forward.ul>>33];
                buffer.n = buffer.n - (cut_end - x_end_index);
            }
            else
            {
                x_beg_index = (uint32_t)position_index[t_forward.ul>>33];
                buffer.n = buffer.n - (x_beg_index - cut_beg);
            }

            if((w&1) == 1)
            {
                y_end_index = (uint32_t)position_index[t_forward.v>>1];
            }
            else
            {
                y_beg_index = (uint32_t)position_index[t_forward.v>>1];
            }

            cut_beg = y_beg_index;
            cut_end = y_end_index;

            if((w&1) == 1)
            {
                for (k = y_end_index; k >= y_beg_index; k--)
                {
                    kv_push(uint64_t, buffer, (yReads->a[k])^(uint64_t)(0x100000000));
                    if(k==0) break;
                }
            }
            else
            {
                for (k = y_beg_index; k <= y_end_index; k++)
                {
                    kv_push(uint64_t, buffer, yReads->a[k]);
                }
            }

            purge_g->seq[w>>1].c = ALTER_LABLE;
        }


        fill_unitig(buffer.a, buffer.n, read_g, edge, 0, &totalLen);

        xReads = &(ug->u.a[keepUid]);
        free(xReads->a);
        xReads->a = buffer.a;
        xReads->n = buffer.n;
        xReads->m = buffer.m;
        xReads->len = totalLen;
        xReads->circ = 0;
        if(xReads->start != (xReads->a[0]>>32))
        {
            xReads->start = xReads->a[0]>>32;
            v = (keepUid<<1)+1;
            av = asg_arc_a(ug->g, v);
            nv = asg_arc_n(ug->g, v);
            for (k = 0; k < nv; k++)
            {
                if(av[k].del) continue;
                asg_arc_del(ug->g, av[k].ul>>32, av[k].v, 1);
                asg_arc_del(ug->g, av[k].v^1, av[k].ul>>32^1, 1);
            }
        }
        
        if(xReads->end != ((xReads->a[xReads->n-1]>>32)^1))
        {
            xReads->end = ((xReads->a[xReads->n-1]>>32)^1);
            v = (keepUid<<1);
            av = asg_arc_a(ug->g, v);
            nv = asg_arc_n(ug->g, v);
            for (k = 0; k < nv; k++)
            {
                if(av[k].del) continue;
                asg_arc_del(ug->g, av[k].ul>>32, av[k].v, 1);
                asg_arc_del(ug->g, av[k].v^1, av[k].ul>>32^1, 1);
            }
        }
    }

    for (i = 0; i < b_0->b.n; i++)
    {
        v = b_0->b.a[i];
        visit[v>>1] = 1;
        if(purge_g->seq[v>1].c != ALTER_LABLE) continue;
        asg_seq_drop(purge_g, v>1);
    }
    
}

void link_unitigs(asg_t *purge_g, ma_ug_t *ug, hap_overlaps_list* all_ovlp,
R_to_U* ruIndex, ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, asg_t *read_g, 
uint64_t* position_index, int max_hang, int min_ovlp, kvec_asg_arc_t_warp* edge, uint8_t* visit)
{
    uint32_t v, n_vtx = purge_g->n_seq * 2, beg, end;
    long long nodeLen, baseLen, max_stop_nodeLen, max_stop_baseLen;
    buf_t b_0;
    memset(&b_0, 0, sizeof(buf_t));
    memset(visit, 0, purge_g->n_seq);
    for (v = 0; v < n_vtx; ++v) 
    {
        if(purge_g->seq[v>>1].c == ALTER_LABLE || purge_g->seq[v>>1].del || visit[v>>1]) continue;
        if(get_real_length(purge_g, v, NULL) != 1) continue;
        if(get_real_length(purge_g, v^1, NULL) != 0) continue;

        beg = v;
        b_0.b.n = 0;
        if(get_unitig(purge_g, NULL, beg, &end, &nodeLen, &baseLen, &max_stop_nodeLen, 
                    &max_stop_baseLen, 1, &b_0) == LOOP)
        {
            continue;
        }
        purge_merge(purge_g, ug, all_ovlp, &b_0, ruIndex, reverse_sources, coverage_cut, 
        read_g, position_index, max_hang, min_ovlp, edge, visit);
    }
    free(b_0.b.a);
}
void purge_dups(ma_ug_t *ug, asg_t *read_g, ma_sub_t* coverage_cut, ma_hit_t_alloc* reverse_sources, 
R_to_U* ruIndex, kvec_asg_arc_t_warp* edge, float density, uint32_t bi_graph_Len, uint32_t long_hap_overlap, 
float lable_match_rate, int max_hang, int min_ovlp, long long bubble_dist, float drop_ratio,
uint32_t just_contain)
{
    asg_t *purge_g = NULL;
    purge_g = asg_init();
    kvec_t_u64_warp u_vecs;
    kv_init(u_vecs.a);
    asg_t* nsg = ug->g;
    uint32_t v, rId, uId, i, k, offset;
    ma_utg_t* reads = NULL;
    uint8_t* visit = NULL;
    visit = (uint8_t*)malloc(sizeof(uint8_t) * nsg->n_seq);
    ///int flag;
    memset(visit, 0, nsg->n_seq);
    uint64_t* position_index = (uint64_t*)malloc(sizeof(uint64_t)*read_g->n_seq);
    memset(position_index, -1, sizeof(uint64_t)*read_g->n_seq);
    uint64_t* vote_counting = (uint64_t*)malloc(sizeof(uint64_t)*nsg->n_seq);
    memset(vote_counting, 0, sizeof(uint64_t)*nsg->n_seq);
    kvec_asg_arc_t_offset u_buffer;
    kv_init(u_buffer.a);
    kvec_hap_candidates u_can;
    kv_init(u_can.a);
    hap_overlaps_list all_ovlp;
    init_hap_overlaps_list(&all_ovlp, nsg->n_seq);
    hap_overlaps_list back_all_ovlp;
    init_hap_overlaps_list(&back_all_ovlp, nsg->n_seq);
    ///uint32_t junk_cov, hap_cov, dip_cov, junk_occ, repeat_occ, single_cov;
    ma_sub_t* purge_cut = NULL;
    purge_cut = (ma_sub_t*)malloc(sizeof(ma_sub_t)*nsg->n_seq);
    asg_arc_t t;
    asg_arc_t* p = NULL;
    int r;

    /****************************may have bugs********************************/
    for (v = 0; v < nsg->n_seq; ++v) 
    {
        uId = v;
        if(nsg->seq[uId].c != HAP_LABLE) continue;
        reads = &(ug->u.a[uId]);
        for (k = 0; k < reads->n; k++)
        {
            rId = reads->a[k]>>33;
            read_g->seq[reads->a[k]>>33].c = HAP_LABLE;
        }
    }
    /****************************may have bugs********************************/

    for (v = 0; v < nsg->n_seq; v++)
    {
        uId = v;
        if(nsg->seq[uId].del || nsg->seq[uId].c == ALTER_LABLE)
        {
            asg_seq_set(purge_g, uId, 0, 1);
            purge_g->seq[uId].c = ALTER_LABLE;
            continue;
        } 
        reads = &(ug->u.a[uId]);
        for (i = 0, offset = 0; i < reads->n; i++)
        {
            rId = reads->a[i]>>33;
            set_R_to_U(ruIndex, rId, uId, 1);

            position_index[rId] = offset;
            position_index[rId] = position_index[rId] << 32;
            position_index[rId] = position_index[rId] | (uint64_t)i;

            offset += (uint32_t)reads->a[i];
        }

        purge_cut[uId].del = 0;
        purge_cut[uId].s = 0;
        purge_cut[uId].e = offset;
        purge_cut[uId].c = PRIMARY_LABLE;
        asg_seq_set(purge_g, uId, offset, 0);
        purge_g->seq[uId].c = PRIMARY_LABLE;
    }


    for (v = 0; v < nsg->n_seq; v++)
    {
        uId = v;
        if(nsg->seq[uId].del || nsg->seq[uId].c == ALTER_LABLE) continue;

        hap_alignment(ug, read_g,  reverse_sources, ruIndex, coverage_cut, position_index, 
        vote_counting, visit, &u_vecs, &u_buffer, &u_can, uId, density, max_hang, min_ovlp, 
        0.05, &all_ovlp);
    }

    normalize_hap_overlaps(&all_ovlp, &back_all_ovlp);
    ///debug_hap_overlaps(&all_ovlp, &back_all_ovlp);
    
    
    for (v = 0; v < all_ovlp.num; v++)
    {
        uId = v;
        for (i = 0; i < all_ovlp.x[uId].a.n; i++)
        {
            if(all_ovlp.x[uId].a.a[i].type == YCX)
            {
                nsg->seq[all_ovlp.x[uId].a.a[i].xUid].c = ALTER_LABLE;
                purge_g->seq[all_ovlp.x[uId].a.a[i].xUid].c = ALTER_LABLE;
                purge_g->seq[all_ovlp.x[uId].a.a[i].xUid].del = 1;
                all_ovlp.x[uId].a.a[i].status = DELETE;
            }

            if(all_ovlp.x[uId].a.a[i].type == XCY)
            {
                nsg->seq[all_ovlp.x[uId].a.a[i].yUid].c = ALTER_LABLE;
                purge_g->seq[all_ovlp.x[uId].a.a[i].yUid].c = ALTER_LABLE;
                purge_g->seq[all_ovlp.x[uId].a.a[i].yUid].del = 1;
                all_ovlp.x[uId].a.a[i].status = DELETE;
            }
            ///print_hap_paf(ug, &(all_ovlp.x[uId].a.a[i]));
        }
    }
    
    if(just_contain == 0)
    {
        for (v = 0; v < all_ovlp.num; v++)
        {
            uId = v;
            if(purge_g->seq[uId].del || purge_g->seq[uId].c == ALTER_LABLE) continue;
            for (i = 0; i < all_ovlp.x[uId].a.n; i++)
            {
                if(all_ovlp.x[uId].a.a[i].status == DELETE) continue;
                if(purge_g->seq[all_ovlp.x[uId].a.a[i].xUid].c == ALTER_LABLE||
                   purge_g->seq[all_ovlp.x[uId].a.a[i].xUid].del||
                   purge_g->seq[all_ovlp.x[uId].a.a[i].yUid].c == ALTER_LABLE||
                   purge_g->seq[all_ovlp.x[uId].a.a[i].yUid].del)
                   {
                       continue;
                   }
                ///print_hap_paf(ug, &(all_ovlp.x[uId].a.a[i]));
                r = get_hap_arch(&(all_ovlp.x[uId].a.a[i]), ug->u.a[all_ovlp.x[uId].a.a[i].xUid].len, 
                ug->u.a[all_ovlp.x[uId].a.a[i].yUid].len, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);

                if (r >= 0) 
                {
                    ///push node?
                    p = asg_arc_pushp(purge_g);
                    *p = t;
                }
                else
                {
                    print_hap_paf(ug, &(all_ovlp.x[uId].a.a[i]));
                    fprintf(stderr, "error\n");
                }
            }
        }

        asg_cleanup(purge_g);
		asg_symm(purge_g);

        clean_purge_graph(purge_g, bubble_dist, drop_ratio);
        ///if(debug_enable) print_purge_gfa(ug, purge_g);
        link_unitigs(purge_g, ug, &all_ovlp, ruIndex, reverse_sources, coverage_cut, read_g, position_index, 
        max_hang, min_ovlp, edge, visit);
    }

    for (v = 0; v < all_ovlp.num; v++)
    {
        uId = v;
        if(purge_g->seq[uId].c == ALTER_LABLE)
        {
            ug->g->seq[uId].c = ALTER_LABLE;
        }
    }

    uint32_t is_Unitig;
    for (v = 0; v < ruIndex->len; v++)
    {
        get_R_to_U(ruIndex, v, &uId, &is_Unitig);
        if(is_Unitig == 1) ruIndex->index[v] = (uint32_t)-1;
    }

    asg_cleanup(nsg);
    kv_destroy(u_vecs.a);
    kv_destroy(u_buffer.a);
    kv_destroy(u_can.a);
    destory_hap_overlaps_list(&all_ovlp);
    destory_hap_overlaps_list(&back_all_ovlp);
    asg_destroy(purge_g);
    free(purge_cut);
    free(position_index);
    free(vote_counting);
    free(visit);
}

