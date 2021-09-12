#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <zlib.h>
#include "Assembly.h"
#include "Process_Read.h"
#include "CommandLines.h"
#include "Hash_Table.h"
#include "POA.h"
#include "Correct.h"
#include "htab.h"
#include "kthread.h"
#include "rcut.h"

void ha_get_candidates_interface(ha_abuf_t *ab, int64_t rid, UC_Read *ucr, overlap_region_alloc *overlap_list, overlap_region_alloc *overlap_list_hp, Candidates_list *cl, double bw_thres, 
int max_n_chain, int keep_whole_chain, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* chain_idx, ma_hit_t_alloc* paf, ma_hit_t_alloc* rev_paf, overlap_region* f_cigar, kvec_t_u64_warp* dbg_ct, st_mt_t *sp);
void ha_get_ug_candidates(ha_abuf_t *ab, int64_t rid, ma_utg_t *u, ma_utg_v *ua, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, int max_n_chain, int keep_whole_chain, kvec_t_u8_warp* k_flag,
kvec_t_u64_warp* chain_idx, void *ha_flt_tab, ha_pt_t *ha_idx, overlap_region* f_cigar, kvec_t_u64_warp* dbg_ct, double chain_match_rate);
void ha_sort_list_by_anchor(overlap_region_alloc *overlap_list);

All_reads R_INF;
Debug_reads R_INF_FLAG;

void get_corrected_read_from_cigar(Cigar_record* cigar, char* pre_read, int pre_length, char* new_read, int* new_length)
{
    int i, j;
    int pre_i, new_i;
    int operation, operation_length;
    pre_i = new_i = 0;
    int diff_char_i = 0;


    for (i = 0; i < (long long)cigar->length; i++)
    {
        operation = Get_Cigar_Type(cigar->record[i]);
        operation_length = Get_Cigar_Length(cigar->record[i]);

        if (operation == 0)
        {
            memcpy(new_read + new_i, pre_read + pre_i, operation_length);
            pre_i = pre_i + operation_length;
            new_i = new_i + operation_length;
        }
        else if (operation == 1)
        {

            for (j = 0; j < operation_length; j++)
            {
                new_read[new_i] = Get_MisMatch_Base(cigar->lost_base[diff_char_i]);
                new_i++;
                diff_char_i++;
            }
            pre_i = pre_i + operation_length;
        }
        else if (operation == 3)
        {
            pre_i = pre_i + operation_length;
            diff_char_i = diff_char_i + operation_length;
        }
        else if (operation == 2)
        {
            memcpy(new_read + new_i, cigar->lost_base + diff_char_i, operation_length);
            new_i = new_i + operation_length;
            diff_char_i = diff_char_i + operation_length;
        }
    }
    *new_length = new_i;
}

void get_uncorrected_read_from_cigar(Cigar_record* cigar, char* new_read, int new_length, char* pre_read, int* pre_length)
{
    int i, j;
    int pre_i, new_i;
    int operation, operation_length;
    pre_i = new_i = 0;
    int diff_char_i = 0;


    for (i = 0; i < (long long)cigar->length; i++)
    {
        operation = Get_Cigar_Type(cigar->record[i]);
        operation_length = Get_Cigar_Length(cigar->record[i]);


        if (operation == 0)
        {
            memcpy(pre_read + pre_i, new_read + new_i, operation_length);
            pre_i = pre_i + operation_length;
            new_i = new_i + operation_length;
        }
        else if (operation == 1)
        {

            for (j = 0; j < operation_length; j++)
            {
                pre_read[pre_i] = Get_Match_Base(cigar->lost_base[diff_char_i]);
                pre_i++;
                diff_char_i++;
            }
            new_i = new_i + operation_length;
        }
        else if (operation == 3)
        {
            memcpy(pre_read + pre_i, cigar->lost_base + diff_char_i, operation_length);
            pre_i = pre_i + operation_length;
            diff_char_i = diff_char_i + operation_length;
        }
        else if (operation == 2)
        {
            new_i = new_i + operation_length;
            diff_char_i = diff_char_i + operation_length;
        }
    }

    *pre_length = pre_i;
}

inline int get_cigar_errors(Cigar_record* cigar)
{
    int i;
    int total_errors = 0;
    for (i = 0; i < (long long)cigar->length; i++)
    {
        if (Get_Cigar_Type(cigar->record[i]) > 0)
        {
            total_errors = total_errors + Get_Cigar_Length(cigar->record[i]);
        }
    }
    return total_errors;
}

int debug_cigar(Cigar_record* cigar, char* pre_read, int pre_length, char* new_read, int new_length, int correct_base)
{
    int i;
    int total_errors = 0;
    for (i = 0; i < (long long)cigar->length; i++)
    {
        if (Get_Cigar_Type(cigar->record[i]) > 0)
        {
            total_errors = total_errors + Get_Cigar_Length(cigar->record[i]);
        }
    }

    if(total_errors!=correct_base)
    {
        fprintf(stderr, "total_errors: %d, correct_base: %d\n", total_errors, correct_base);
    }

    int pre_i, new_i;
    int operation, operation_length;
    pre_i = new_i = 0;

    for (i = 0; i < (long long)cigar->length; i++)
    {
        operation = Get_Cigar_Type(cigar->record[i]);
        operation_length = Get_Cigar_Length(cigar->record[i]);

        if (operation == 0)
        {
            pre_i = pre_i + operation_length;
            new_i = new_i + operation_length;
        }

        if (operation == 1)
        {
            pre_i = pre_i + operation_length;
            new_i = new_i + operation_length;
        }

        if (operation == 3)
        {
            pre_i = pre_i + operation_length;
        }

        if (operation == 2)
        {
            new_i = new_i + operation_length;
        }

    }

    if (pre_i != pre_length)
    {
        fprintf(stderr, "pre_i: %d, pre_length: %d\n", pre_i, pre_length);
    }


    if(new_i != new_length)
    {
        fprintf(stderr, "new_i: %d, new_length: %d\n", new_i, new_length);
    }

    return 1;

    char* tmp_seq = (char*)malloc(new_length + pre_length);
    int tmp_length;

    get_corrected_read_from_cigar(cigar, pre_read, pre_length, tmp_seq, &tmp_length);

    if(tmp_length != new_length)
    {
        fprintf(stderr, "tmp_length: %d, new_length: %d\n", tmp_length, new_length);
    }

    if(memcmp(new_read, tmp_seq, new_length)!=0)
    {
        fprintf(stderr, "error new string\n");
    }


    get_uncorrected_read_from_cigar(cigar, new_read, new_length, tmp_seq, &tmp_length);


    if(tmp_length != pre_length)
    {
        fprintf(stderr, "tmp_length: %d, pre_length: %d\n", tmp_length, pre_length);
    }

    if(memcmp(pre_read, tmp_seq, pre_length)!=0)
    {
        fprintf(stderr, "error pre string\n");
    }


    free(tmp_seq);

    if((int)cigar->new_read_length != new_length)
    {
        fprintf(stderr, "cigar->new_read_length: %d, new_length: %d\n", cigar->new_read_length, new_length);
    }
}


inline void push_cigar(Compressed_Cigar_record* records, long long ID, Cigar_record* input)
{
    if (input->length > records[ID].size)
    {
        records[ID].size = input->length;
        records[ID].record = (uint32_t*)realloc(records[ID].record, records[ID].size*sizeof(uint32_t));     
    }
    records[ID].length = input->length;
    memcpy(records[ID].record, input->record, input->length*sizeof(uint32_t));

    if (input->lost_base_length > records[ID].lost_base_size)
    {
        records[ID].lost_base_size = input->lost_base_length;
        records[ID].lost_base = (char*)realloc(records[ID].lost_base, records[ID].lost_base_size);
    }
    records[ID].lost_base_length = input->lost_base_length;
    memcpy(records[ID].lost_base, input->lost_base, input->lost_base_length);

    records[ID].new_length = input->new_read_length;
}


void push_overlaps(ma_hit_t_alloc* paf, overlap_region_alloc* overlap_list, int flag, All_reads* R_INF, int if_reverse)
{
    long long i = 0, xLen, yLen;
	int32_t size = 0;
    ma_hit_t tmp;
	for (i = 0; i < (long long)overlap_list->length; ++i)
		if (overlap_list->list[i].is_match == flag)
			++size;
	resize_ma_hit_t_alloc(paf, size);
    clear_ma_hit_t_alloc(paf);
    for (i = 0; i < (long long)overlap_list->length; i++)
    {
        if (overlap_list->list[i].is_match == flag)
        {
            xLen = Get_READ_LENGTH((*R_INF), overlap_list->list[i].x_id);
            yLen = Get_READ_LENGTH((*R_INF), overlap_list->list[i].y_id);

            tmp.qns = overlap_list->list[i].x_id;
            tmp.qns = tmp.qns << 32;
            tmp.tn = overlap_list->list[i].y_id;

            if(if_reverse != 0)
            {
                tmp.qns = tmp.qns | (uint64_t)(xLen - overlap_list->list[i].x_pos_s - 1);
                tmp.qe = xLen - overlap_list->list[i].x_pos_e - 1;
                tmp.ts = yLen - overlap_list->list[i].y_pos_s - 1;
                tmp.te = yLen - overlap_list->list[i].y_pos_e - 1;
            }
            else
            {
                tmp.qns = tmp.qns | (uint64_t)(overlap_list->list[i].x_pos_s);
                tmp.qe = overlap_list->list[i].x_pos_e;
                tmp.ts = overlap_list->list[i].y_pos_s;
                tmp.te = overlap_list->list[i].y_pos_e;
            }

            ///for overlap_list, the x_strand of all overlaps are 0, so the tmp.rev is the same as the y_strand
            tmp.rev = overlap_list->list[i].y_pos_strand;

            ///tmp.bl = R_INF.read_length[overlap_list->list[i].y_id];
            tmp.bl = Get_READ_LENGTH((*R_INF), overlap_list->list[i].y_id);
            tmp.ml = overlap_list->list[i].strong;
            tmp.no_l_indel = overlap_list->list[i].without_large_indel;

            add_ma_hit_t_alloc(paf, &tmp);
        }
    }
}


long long push_final_overlaps(ma_hit_t_alloc* paf, ma_hit_t_alloc* reverse_paf_list, overlap_region_alloc* overlap_list, int flag)
{
    long long i = 0;
    long long available_overlaps = 0;
    ma_hit_t tmp;
    clear_ma_hit_t_alloc(paf); // paf has been preallocated, so we don't need preallocation
    for (i = 0; i < (long long)overlap_list->length; i++)
    {
        if (overlap_list->list[i].is_match == flag)
        {
            available_overlaps++;
            /**********************query***************************/
            //the interval of overlap is half-open [start, end) 
            tmp.qns = overlap_list->list[i].x_id;
            tmp.qns = tmp.qns << 32;
            tmp.qns = tmp.qns | (uint64_t)(overlap_list->list[i].x_pos_s);
            ///the end pos is open
            tmp.qe = overlap_list->list[i].x_pos_e + 1;
            /**********************query***************************/



            ///for overlap_list, the x_strand of all overlaps are 0, so the tmp.rev is the same as the y_strand
            tmp.rev = overlap_list->list[i].y_pos_strand;


            /**********************target***************************/
            tmp.tn = overlap_list->list[i].y_id;
            if(tmp.rev == 1)
            {
                long long y_readLen = R_INF.read_length[overlap_list->list[i].y_id];
                tmp.ts = y_readLen - overlap_list->list[i].y_pos_e - 1;
                tmp.te = y_readLen - overlap_list->list[i].y_pos_s - 1;
            }
            else
            {
                tmp.ts = overlap_list->list[i].y_pos_s;
                tmp.te = overlap_list->list[i].y_pos_e;
            }
            ///the end pos is open
            tmp.te++;
            /**********************target***************************/

            tmp.bl = R_INF.read_length[overlap_list->list[i].y_id];
            tmp.ml = overlap_list->list[i].strong;
            tmp.no_l_indel = overlap_list->list[i].without_large_indel;

            tmp.el = overlap_list->list[i].shared_seed;

            add_ma_hit_t_alloc(paf, &tmp);
        }
    }

    return available_overlaps;
}



long long push_final_overlaps_increment(ma_hit_t_alloc* paf, ma_hit_t_alloc* reverse_paf_list, overlap_region_alloc* overlap_list, int flag)
{
    long long i = 0;
    long long available_overlaps = paf->length;
    ma_hit_t tmp;
    ///clear_ma_hit_t_alloc(paf); // paf has been preallocated, so we don't need preallocation
    for (i = 0; i < (long long)overlap_list->length; i++)
    {
        if (overlap_list->list[i].is_match == flag)
        {
            available_overlaps++;
            /**********************query***************************/
            //the interval of overlap is half-open [start, end) 
            tmp.qns = overlap_list->list[i].x_id;
            tmp.qns = tmp.qns << 32;
            tmp.qns = tmp.qns | (uint64_t)(overlap_list->list[i].x_pos_s);
            ///the end pos is open
            tmp.qe = overlap_list->list[i].x_pos_e + 1;
            /**********************query***************************/



            ///for overlap_list, the x_strand of all overlaps are 0, so the tmp.rev is the same as the y_strand
            tmp.rev = overlap_list->list[i].y_pos_strand;


            /**********************target***************************/
            tmp.tn = overlap_list->list[i].y_id;
            if(tmp.rev == 1)
            {
                long long y_readLen = R_INF.read_length[overlap_list->list[i].y_id];
                tmp.ts = y_readLen - overlap_list->list[i].y_pos_e - 1;
                tmp.te = y_readLen - overlap_list->list[i].y_pos_s - 1;
            }
            else
            {
                tmp.ts = overlap_list->list[i].y_pos_s;
                tmp.te = overlap_list->list[i].y_pos_e;
            }
            ///the end pos is open
            tmp.te++;
            /**********************target***************************/

            tmp.bl = R_INF.read_length[overlap_list->list[i].y_id];
            tmp.ml = overlap_list->list[i].strong;
            tmp.no_l_indel = overlap_list->list[i].without_large_indel;

            tmp.el = overlap_list->list[i].shared_seed;

            add_ma_hit_t_alloc(paf, &tmp);
        }
    }

    return available_overlaps;
}

typedef struct {
	int is_final, save_ov;
	// chaining and overlapping related buffers
	UC_Read self_read, ovlp_read;
	Candidates_list clist;
	overlap_region_alloc olist;
    overlap_region_alloc olist_hp;
	ha_abuf_t *ab;
    ha_abufl_t *abl;
	// error correction related buffers
	int64_t num_read_base, num_correct_base, num_recorrect_base;
	Cigar_record cigar1;
	Graph POA_Graph;
	Graph DAGCon;
	Correct_dumy correct;
	haplotype_evdience_alloc hap;
	Round2_alignment round2;
    kvec_t_u32_warp b_buf;
    kvec_t_u64_warp r_buf;
    kvec_t_u8_warp k_flag;
    overlap_region tmp_region;
    ma_utg_v *ua;
    st_mt_t sp;
} ha_ovec_buf_t;

ha_ovec_buf_t *ha_ovec_init(int is_final, int save_ov, int is_ug)
{
	ha_ovec_buf_t *b;
	CALLOC(b, 1);
	b->is_final = !!is_final, b->save_ov = !!save_ov;
	init_UC_Read(&b->self_read);
	init_UC_Read(&b->ovlp_read);
	init_Candidates_list(&b->clist);
	init_overlap_region_alloc(&b->olist);
    init_overlap_region_alloc(&b->olist_hp);
    init_fake_cigar(&(b->tmp_region.f_cigar));
    kv_init(b->b_buf.a);
    kv_init(b->r_buf.a);
    kv_init(b->k_flag.a);
    kv_init(b->sp);
	if(!is_ug) b->ab = ha_abuf_init();
    else b->abl = ha_abufl_init();
	if (!b->is_final) {
		init_Cigar_record(&b->cigar1);
		init_Graph(&b->POA_Graph);
		init_Graph(&b->DAGCon);
		init_Correct_dumy(&b->correct);
		InitHaplotypeEvdience(&b->hap);
		init_Round2_alignment(&b->round2);
	}
	return b;
}

void ha_ovec_destroy(ha_ovec_buf_t *b)
{
	destory_UC_Read(&b->self_read);
	destory_UC_Read(&b->ovlp_read);
	destory_Candidates_list(&b->clist);
	destory_overlap_region_alloc(&b->olist);
    destory_overlap_region_alloc(&b->olist_hp);
	ha_abuf_destroy(b->ab);
    ha_abufl_destroy(b->abl);
    destory_fake_cigar(&(b->tmp_region.f_cigar));
    kv_destroy(b->b_buf.a);
    kv_destroy(b->r_buf.a);
    kv_destroy(b->k_flag.a);
    kv_destroy(b->sp);
	if (!b->is_final) {
		destory_Cigar_record(&b->cigar1);
		destory_Graph(&b->POA_Graph);
		destory_Graph(&b->DAGCon);
		destory_Correct_dumy(&b->correct);
		destoryHaplotypeEvdience(&b->hap);
		destory_Round2_alignment(&b->round2);
	}
	free(b);
}

static int64_t ha_Graph_mem(const Graph *g)
{
	int64_t i, mem = 0;
	mem = sizeof(Graph) + g->node_q.size * 8 + g->g_nodes.size * sizeof(Node);
	for (i = 0; i < (int64_t)g->g_nodes.size; ++i) {
		Node *n = &g->g_nodes.list[i];
		mem += n->mismatch_edges.size * sizeof(Edge);
		mem += n->deletion_edges.size * sizeof(Edge);
		mem += n->insertion_edges.size * sizeof(Edge);
	}
	mem += g->g_nodes.sort.size * 9;
	return mem;
}

int64_t ha_ovec_mem(const ha_ovec_buf_t *b)
{
	int64_t i, mem = 0, mem_clist, mem_olist;
	mem_clist = b->clist.size * sizeof(k_mer_hit) + b->clist.chainDP.size * 7 * 4;
    
	mem_olist = b->olist.size * sizeof(overlap_region);
	for (i = 0; i < (int64_t)b->olist.size; ++i) {
		const overlap_region *r = &b->olist.list[i];
		mem_olist += r->w_list_size * sizeof(window_list);
		mem_olist += r->f_cigar.size * 8;
		mem_olist += r->boundary_cigars.size * sizeof(window_list);
	}
    mem_olist += b->olist_hp.size * sizeof(overlap_region);
	for (i = 0; i < (int64_t)b->olist_hp.size; ++i) {
		const overlap_region *r = &b->olist_hp.list[i];
		mem_olist += r->w_list_size * sizeof(window_list);
		mem_olist += r->f_cigar.size * 8;
		mem_olist += r->boundary_cigars.size * sizeof(window_list);
	}

    if(b->ab) mem = ha_abuf_mem(b->ab) + mem_clist + mem_olist;
    if(b->abl) mem = ha_abufl_mem(b->abl) + mem_clist + mem_olist;
	if (!b->is_final) {
		mem += sizeof(Cigar_record) + b->cigar1.lost_base_size + b->cigar1.size * 4;
		mem += sizeof(Correct_dumy) + b->correct.size * 8;
		mem += sizeof(Round2_alignment) + b->round2.cigar.size * 4 + b->round2.tmp_cigar.size * 4;
		mem += sizeof(haplotype_evdience_alloc) + b->hap.size * sizeof(haplotype_evdience) + b->hap.snp_matrix_size + b->hap.snp_stat_size * sizeof(SnpStats);
		mem += ha_Graph_mem(&b->POA_Graph);
		mem += ha_Graph_mem(&b->DAGCon);
	}
	return mem;
}

static void worker_ovec(void *data, long i, int tid)
{
	ha_ovec_buf_t *b = ((ha_ovec_buf_t**)data)[tid];
	int fully_cov, abnormal;

    ha_get_candidates_interface(b->ab, i, &b->self_read, &b->olist, &b->olist_hp, &b->clist, 
    0.02, asm_opt.max_n_chain, 1, &(b->k_flag), &b->r_buf, &(R_INF.paf[i]), &(R_INF.reverse_paf[i]), &(b->tmp_region), NULL, &(b->sp));

	clear_Cigar_record(&b->cigar1);
	clear_Round2_alignment(&b->round2);

	correct_overlap(&b->olist, &R_INF, &b->self_read, &b->correct, &b->ovlp_read, &b->POA_Graph, &b->DAGCon,
			&b->cigar1, &b->hap, &b->round2, 0, 1, &fully_cov, &abnormal);

	b->num_read_base += b->self_read.length;
	b->num_correct_base += b->correct.corrected_base;
	b->num_recorrect_base += b->round2.dumy.corrected_base;

	push_cigar(R_INF.cigars, i, &b->cigar1);
	push_cigar(R_INF.second_round_cigar, i, &b->round2.cigar);

	R_INF.paf[i].is_fully_corrected = 0;
	if (fully_cov) {
		if (get_cigar_errors(&b->cigar1) == 0 && get_cigar_errors(&b->round2.cigar) == 0)
			R_INF.paf[i].is_fully_corrected = 1;
	}
	R_INF.paf[i].is_abnormal = abnormal;

    R_INF.trio_flag[i] = AMBIGU;
    
    ///need to be fixed in r305
    // if(ha_idx_hp == NULL)
    // {
    //     R_INF.trio_flag[i] += collect_hp_regions(&b->olist, &R_INF, &(b->k_flag), RESEED_HP_RATE, Get_READ_LENGTH(R_INF, i), NULL);
    // }

    if (R_INF.trio_flag[i] != AMBIGU || b->save_ov) {
		int is_rev = (asm_opt.number_of_round % 2 == 0);
		push_overlaps(&(R_INF.paf[i]), &b->olist, 1, &R_INF, is_rev);
		push_overlaps(&(R_INF.reverse_paf[i]), &b->olist, 2, &R_INF, is_rev);
	}
}


static void worker_ovec_related_reads(void *data, long i, int tid)
{
	ha_ovec_buf_t *b = ((ha_ovec_buf_t**)data)[tid];

    uint64_t k, queryNameLen;
    for (k = 0; k < R_INF_FLAG.query_num; k++)
    {
        queryNameLen = strlen(R_INF_FLAG.read_name[k]);
        if (queryNameLen != Get_NAME_LENGTH((R_INF),i)) continue;
        if (memcmp(R_INF_FLAG.read_name[k], Get_NAME((R_INF), i), Get_NAME_LENGTH((R_INF),i)) == 0)
        {
            break;
        }
    }

    if(k < R_INF_FLAG.query_num)
    {
        int fully_cov, abnormal, q_idx = k;

        ha_get_candidates_interface(b->ab, i, &b->self_read, &b->olist, &b->olist_hp, &b->clist, 
        0.02, asm_opt.max_n_chain, 1, &(b->k_flag), &b->r_buf, &(R_INF.paf[i]), &(R_INF.reverse_paf[i]), &(b->tmp_region), &(R_INF_FLAG.candidate_count[q_idx]), &(b->sp));

        clear_Cigar_record(&b->cigar1);
        clear_Round2_alignment(&b->round2);

        correct_overlap(&b->olist, &R_INF, &b->self_read, &b->correct, &b->ovlp_read, &b->POA_Graph, &b->DAGCon,
                &b->cigar1, &b->hap, &b->round2, 0, 1, &fully_cov, &abnormal);

        b->num_read_base += b->self_read.length;
        b->num_correct_base += b->correct.corrected_base;
        b->num_recorrect_base += b->round2.dumy.corrected_base;

        push_cigar(R_INF.cigars, i, &b->cigar1);
        push_cigar(R_INF.second_round_cigar, i, &b->round2.cigar);

        R_INF.paf[i].is_fully_corrected = 0;
        if (fully_cov) {
            if (get_cigar_errors(&b->cigar1) == 0 && get_cigar_errors(&b->round2.cigar) == 0)
                R_INF.paf[i].is_fully_corrected = 1;
        }
        R_INF.paf[i].is_abnormal = abnormal;

        

        pthread_mutex_lock(&R_INF_FLAG.OutputMutex);

        fprintf(R_INF_FLAG.fp, "\n>%.*s, rid: %ld\n", (int)Get_NAME_LENGTH((R_INF), i), Get_NAME((R_INF), i), i);
        fprintf(R_INF_FLAG.fp, "%d-th round, len: %lu, hom_cov: %d, max_n_chain: %d\n", 
            asm_opt.number_of_round, Get_READ_LENGTH(R_INF, i), asm_opt.hom_cov, asm_opt.max_n_chain);

        fprintf(R_INF_FLAG.fp, "***************************k-mer counts (%d)***************************\n", (int)(R_INF_FLAG.candidate_count[q_idx].a.n));
        
        sort_kvec_t_u64_warp(&(R_INF_FLAG.candidate_count[q_idx]), 0);
        for (k = 0; k < R_INF_FLAG.candidate_count[q_idx].a.n; k++) 
        {
            fprintf(R_INF_FLAG.fp, "[%lu] Count(%u): %lu, filtered: %lu\n", k, 
            (uint32_t)R_INF_FLAG.candidate_count[q_idx].a.a[k], R_INF_FLAG.candidate_count[q_idx].a.a[k]>>33, 
            (R_INF_FLAG.candidate_count[q_idx].a.a[k]>>32)&(uint64_t)1);
        }

        
        fprintf(R_INF_FLAG.fp, "***************************forward ovlp***************************\n");
        for (k = 0; k < b->olist.length; k++) 
        {
            if(b->olist.list[k].is_match != 1) continue;
			fprintf(R_INF_FLAG.fp, "%.*s\n", (int)Get_NAME_LENGTH((R_INF), b->olist.list[k].y_id), Get_NAME((R_INF), b->olist.list[k].y_id));
            fprintf(R_INF_FLAG.fp, "qs: %u, qe: %u, ts: %u, te: %u, rev: %u, strong: %u, no_l_indel: %u, len: %lu\n",
            b->olist.list[k].x_pos_s, b->olist.list[k].x_pos_e, b->olist.list[k].y_pos_s, b->olist.list[k].y_pos_e, 
            b->olist.list[k].y_pos_strand, b->olist.list[k].strong, b->olist.list[k].without_large_indel,
            Get_READ_LENGTH(R_INF, b->olist.list[k].y_id));
		}

        fprintf(R_INF_FLAG.fp, "***************************reverse ovlp***************************\n");
        for (k = 0; k < b->olist.length; k++) 
        {
            if(b->olist.list[k].is_match != 2) continue;
			fprintf(R_INF_FLAG.fp, "%.*s\n", (int)Get_NAME_LENGTH((R_INF), b->olist.list[k].y_id), Get_NAME((R_INF), b->olist.list[k].y_id));
            fprintf(R_INF_FLAG.fp, "qs: %u, qe: %u, ts: %u, te: %u, rev: %u, strong: %u, no_l_indel: %u, len: %lu\n",
            b->olist.list[k].x_pos_s, b->olist.list[k].x_pos_e, b->olist.list[k].y_pos_s, b->olist.list[k].y_pos_e, 
            b->olist.list[k].y_pos_strand, b->olist.list[k].strong, b->olist.list[k].without_large_indel,
            Get_READ_LENGTH(R_INF, b->olist.list[k].y_id));
		}

        fprintf(R_INF_FLAG.fp, "***************************unmatched ovlp***************************\n");
        for (k = 0; k < b->olist.length; k++) 
        {
            if(b->olist.list[k].is_match == 1) continue;
            if(b->olist.list[k].is_match == 2) continue;
			fprintf(R_INF_FLAG.fp, "%.*s\n", (int)Get_NAME_LENGTH((R_INF), b->olist.list[k].y_id), Get_NAME((R_INF), b->olist.list[k].y_id));
            fprintf(R_INF_FLAG.fp, "qs: %u, qe: %u, ts: %u, te: %u, rev: %u, strong: %u, no_l_indel: %u, len: %lu\n",
            b->olist.list[k].x_pos_s, b->olist.list[k].x_pos_e, b->olist.list[k].y_pos_s, b->olist.list[k].y_pos_e, 
            b->olist.list[k].y_pos_strand, b->olist.list[k].strong, b->olist.list[k].without_large_indel,
            Get_READ_LENGTH(R_INF, b->olist.list[k].y_id));
		}

        R_INF.trio_flag[i] = AMBIGU;

        ///need to be fixed in r305
        // if(ha_idx_hp == NULL)
        // {
        //     R_INF.trio_flag[i] += collect_hp_regions(&b->olist, &R_INF, &(b->k_flag), RESEED_HP_RATE, Get_READ_LENGTH(R_INF, i), R_INF_FLAG.fp);
        // }
        
        fprintf(R_INF_FLAG.fp, "R_INF.trio_flag[%ld]: %u\n", i, R_INF.trio_flag[i]);
    

	    pthread_mutex_unlock(&R_INF_FLAG.OutputMutex);
    }
}


static inline long long get_N_occ(char* seq, long long length)
{
	long long j, N_occ = 0;
	for (j = 0; j < length; j++)
		if(seq_nt6_table[(uint8_t)seq[j]] >= 4)
			N_occ++;
	return N_occ;
}

typedef struct {
	UC_Read g_read;
	int first_round_read_size;
	int second_round_read_size;
	char *first_round_read;
	char *second_round_read;
} ha_ecsave_buf_t;

static void worker_ec_save(void *data, long i, int tid)
{
	ha_ecsave_buf_t *e = (ha_ecsave_buf_t*)data + tid;

	Cigar_record cigar;
	int first_round_read_length;
	int second_round_read_length;
	uint64_t N_occ;

	char *new_read;
	int new_read_length;

	recover_UC_Read(&e->g_read, &R_INF, i);

	// round 1
	if ((long long)R_INF.cigars[i].new_length > e->first_round_read_size) {
		e->first_round_read_size = R_INF.cigars[i].new_length;
		REALLOC(e->first_round_read, e->first_round_read_size);
	}

	cigar.length = R_INF.cigars[i].length;
	cigar.lost_base_length = R_INF.cigars[i].lost_base_length;
	cigar.record = R_INF.cigars[i].record;
	cigar.lost_base = R_INF.cigars[i].lost_base;

	get_corrected_read_from_cigar(&cigar, e->g_read.seq, e->g_read.length, e->first_round_read, &first_round_read_length);

	// round 2
	if ((long long)R_INF.second_round_cigar[i].new_length > e->second_round_read_size) {
		e->second_round_read_size = R_INF.second_round_cigar[i].new_length;
		REALLOC(e->second_round_read, e->second_round_read_size);
	}
	cigar.length = R_INF.second_round_cigar[i].length;
	cigar.lost_base_length = R_INF.second_round_cigar[i].lost_base_length;
	cigar.record = R_INF.second_round_cigar[i].record;
	cigar.lost_base = R_INF.second_round_cigar[i].lost_base;

	get_corrected_read_from_cigar(&cigar, e->first_round_read, first_round_read_length, e->second_round_read, &second_round_read_length);

	new_read = e->second_round_read;
	new_read_length = second_round_read_length;

	if (asm_opt.roundID != asm_opt.number_of_round - 1)
	{
		///need modification
		reverse_complement(new_read, new_read_length);
	}
	else if(asm_opt.number_of_round % 2 == 0)
	{
		///need modification
		reverse_complement(new_read, new_read_length);
	}

	N_occ = get_N_occ(new_read, new_read_length);

	if ((long long)R_INF.read_size[i] < new_read_length) {
		R_INF.read_size[i] = new_read_length;
		REALLOC(R_INF.read_sperate[i], R_INF.read_size[i]/4+1);
	}
	R_INF.read_length[i] = new_read_length;
	ha_compress_base(Get_READ(R_INF, i), new_read, new_read_length, &R_INF.N_site[i], N_occ);
}

void Output_corrected_reads()
{
    long long i;
    UC_Read g_read;
    init_UC_Read(&g_read);
    char* gfa_name = (char*)malloc(strlen(asm_opt.output_file_name)+35);
    sprintf(gfa_name, "%s.ec.fa", asm_opt.output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    free(gfa_name);

    for (i = 0; i < (long long)R_INF.total_reads; i++)
    {
        recover_UC_Read(&g_read, &R_INF, i);
        fwrite(">", 1, 1, output_file);
        fwrite(Get_NAME(R_INF, i), 1, Get_NAME_LENGTH(R_INF, i), output_file);
        fwrite("\n", 1, 1, output_file);
        fwrite(g_read.seq, 1, g_read.length, output_file);
        fwrite("\n", 1, 1, output_file);
    }
    destory_UC_Read(&g_read);
    fclose(output_file);
}

void debug_print_pob_regions()
{
    uint64_t i, total = 0;
    for (i = 0; i < R_INF.total_reads; i++)
    {
        if(R_INF.trio_flag[i]!=AMBIGU)
        {
            total++;
            fprintf(stderr, "(%lu) %.*s\n", i, (int)Get_NAME_LENGTH(R_INF, i), Get_NAME(R_INF, i));
        }
    }
    fprintf(stderr, "total hp reads: %lu, R_INF.total_reads: %lu\n", total, R_INF.total_reads);
    exit(1);
}

void rescue_hp_reads(ha_ovec_buf_t **b)
{
    int hom_cov, het_cov;
    ha_flt_tab_hp = ha_idx_hp = NULL;
    if (!(asm_opt.flag & HA_F_NO_KMER_FLT)) {
       ha_flt_tab_hp = ha_ft_gen(&asm_opt, &R_INF, &hom_cov, 1);
    }
    ha_idx_hp = ha_pt_gen(&asm_opt, ha_flt_tab, 1, 1, &R_INF, &hom_cov, &het_cov);


    if (asm_opt.required_read_name)
		kt_for(asm_opt.thread_num, worker_ovec_related_reads, b, R_INF.total_reads);
	else
		kt_for(asm_opt.thread_num, worker_ovec, b, R_INF.total_reads);




    ha_ft_destroy(ha_flt_tab_hp); ha_flt_tab_hp = NULL;
    ha_pt_destroy(ha_idx_hp); ha_idx_hp = NULL;
}




void ha_overlap_and_correct(int round)
{
	int i, hom_cov, het_cov, r_out = 0;
	ha_ovec_buf_t **b;
	ha_ecsave_buf_t *e;
    ha_flt_tab_hp = ha_idx_hp = NULL;

    if((ha_idx == NULL)&&(asm_opt.flag & HA_F_VERBOSE_GFA)&&(round == asm_opt.number_of_round - 1))
    {
        r_out = 1;
    }

    if(asm_opt.required_read_name) init_Debug_reads(&R_INF_FLAG, asm_opt.required_read_name); // for debugging only
	// overlap and correct reads
	CALLOC(b, asm_opt.thread_num);
	for (i = 0; i < asm_opt.thread_num; ++i)
		b[i] = ha_ovec_init(0, (round == asm_opt.number_of_round - 1),0);
    if(ha_idx) hom_cov = asm_opt.hom_cov;
	if(ha_idx == NULL) ha_idx = ha_pt_gen(&asm_opt, ha_flt_tab, round == 0? 0 : 1, 0, &R_INF, &hom_cov, &het_cov); // build the index
	///debug_adapter(&asm_opt, &R_INF);
    if (round == 0 && ha_flt_tab == 0) // then asm_opt.hom_cov hasn't been updated
		ha_opt_update_cov(&asm_opt, hom_cov);
	if (asm_opt.required_read_name)
		kt_for(asm_opt.thread_num, worker_ovec_related_reads, b, R_INF.total_reads);
	else
		kt_for(asm_opt.thread_num, worker_ovec, b, R_INF.total_reads);

    if (r_out) write_pt_index(ha_flt_tab, ha_idx, &R_INF, &asm_opt, asm_opt.output_file_name);
	ha_pt_destroy(ha_idx);
	ha_idx = NULL;

	// collect statistics
	for (i = 0; i < asm_opt.thread_num; ++i) {
		asm_opt.num_bases += b[i]->num_read_base;
		asm_opt.num_corrected_bases += b[i]->num_correct_base;
		asm_opt.num_recorrected_bases += b[i]->num_recorrect_base;
		asm_opt.mem_buf += ha_ovec_mem(b[i]);
		ha_ovec_destroy(b[i]);
	}
	free(b);

	if (asm_opt.required_read_name) destory_Debug_reads(&R_INF_FLAG), exit(0); // for debugging only
    
	// save corrected reads to R_INF
	CALLOC(e, asm_opt.thread_num);
	for (i = 0; i < asm_opt.thread_num; ++i) {
		init_UC_Read(&e[i].g_read);
		e[i].first_round_read_size = e[i].second_round_read_size = 50000;
		CALLOC(e[i].first_round_read, e[i].first_round_read_size);
		CALLOC(e[i].second_round_read, e[i].second_round_read_size);
	}
	kt_for(asm_opt.thread_num, worker_ec_save, e, R_INF.total_reads);
	for (i = 0; i < asm_opt.thread_num; ++i) {
		destory_UC_Read(&e[i].g_read);
		free(e[i].first_round_read);
		free(e[i].second_round_read);
	}
	free(e);
    ///debug_print_pob_regions();
}


void update_overlaps(overlap_region_alloc* overlap_list, ma_hit_t_alloc* paf, 
UC_Read* g_read, UC_Read* overlap_read, int is_match, int is_exact)
{

    uint64_t inner_j = 0;
    uint64_t j = 0;
    long long x_overlapLen, y_overlapLen;
    while (j < overlap_list->length && inner_j < paf->length)
    {
        if(overlap_list->list[j].y_id < paf->buffer[inner_j].tn)
        {
            j++;
        }
        else if(overlap_list->list[j].y_id > paf->buffer[inner_j].tn)
        {
            inner_j++;
        }
        else
        {
            if(overlap_list->list[j].y_pos_strand == paf->buffer[inner_j].rev)
            {
                x_overlapLen = Get_qe(paf->buffer[inner_j]) - Get_qs(paf->buffer[inner_j]) + 1;
                y_overlapLen = Get_te(paf->buffer[inner_j]) - Get_ts(paf->buffer[inner_j]) + 1;
                if(x_overlapLen < y_overlapLen) x_overlapLen = y_overlapLen;
                x_overlapLen = x_overlapLen * 0.1;

                // if(
                // ((DIFF(overlap_list->list[j].x_pos_s, Get_qs(paf->buffer[inner_j])) < x_overlapLen)
                // && (DIFF(overlap_list->list[j].x_pos_e, Get_qe(paf->buffer[inner_j])) < x_overlapLen))
                // ||
                // ((DIFF(overlap_list->list[j].y_pos_s, Get_ts(paf->buffer[inner_j])) < x_overlapLen)
                // && (DIFF(overlap_list->list[j].y_pos_e, Get_te(paf->buffer[inner_j])) < x_overlapLen)))
                if(
                ((DIFF(overlap_list->list[j].x_pos_s, Get_qs(paf->buffer[inner_j])) < (uint64_t)x_overlapLen)
                && (DIFF(overlap_list->list[j].x_pos_e, Get_qe(paf->buffer[inner_j])) < (uint64_t)x_overlapLen))
                || 
                ((DIFF(overlap_list->list[j].y_pos_s, Get_ts(paf->buffer[inner_j])) < (uint64_t)x_overlapLen)
                && (DIFF(overlap_list->list[j].y_pos_e, Get_te(paf->buffer[inner_j])) < (uint64_t)x_overlapLen))
                )
                {
                    overlap_list->list[j].is_match = is_match;
                    overlap_list->list[j].strong = paf->buffer[inner_j].ml;
                    overlap_list->list[j].without_large_indel = paf->buffer[inner_j].no_l_indel;
                    if(is_exact == 1)
                    {
                        if(overlap_list->list[j].y_pos_strand == 0)
                        {
                            recover_UC_Read(overlap_read, &R_INF, overlap_list->list[j].y_id);
                        }
                        else
                        {
                            recover_UC_Read_RC(overlap_read, &R_INF, overlap_list->list[j].y_id);
                        }
                        if(if_exact_match(g_read->seq, g_read->length, overlap_read->seq, overlap_read->length, 
                        overlap_list->list[j].x_pos_s, overlap_list->list[j].x_pos_e, 
                        overlap_list->list[j].y_pos_s, overlap_list->list[j].y_pos_e))
                        {
                            overlap_list->list[j].shared_seed = 1;
                        }
                        else
                        {
                            overlap_list->list[j].shared_seed = 0;
                        }
                    }
                }
                else
                {
                    overlap_list->list[j].is_match = 3;
                }
            }
            else
            {
                overlap_list->list[j].is_match = 3;
            }

            j++;
            inner_j++;
        }
    }
}


int check_chain_indels(Fake_Cigar* chain, long long xBeg, long long xEnd, float indel_rate)
{
    uint64_t i = 0;
    long long indels = 0, xOffset;
    if(chain->length != 0)
    {
        indels += abs(get_fake_gap_shift(chain, 0));
        xOffset = get_fake_gap_pos(chain, 0);
        if(indels > (xOffset - xBeg + 1) * indel_rate) return 0;

        for (i = 1; i < chain->length; i++)
        {
            indels += abs((get_fake_gap_shift(chain, i) - get_fake_gap_shift(chain, i-1)));
            xOffset = get_fake_gap_pos(chain, i);
            if(indels > (xOffset - xBeg + 1) * indel_rate) return 0;
        }
    }

    if(indels > (xEnd - xBeg + 1) * indel_rate) return 0;
    return 1;
}

void update_overlaps_chain_width(overlap_region_alloc* overlap_list, ma_hit_t_alloc* paf, 
UC_Read* g_read, UC_Read* overlap_read, int is_match, int is_exact, float indel_rate)
{

    uint64_t inner_j = 0;
    uint64_t j = 0;
    long long x_overlapLen, y_overlapLen;
    while (j < overlap_list->length && inner_j < paf->length)
    {
        if(overlap_list->list[j].y_id < paf->buffer[inner_j].tn)
        {
            j++;
        }
        else if(overlap_list->list[j].y_id > paf->buffer[inner_j].tn)
        {
            inner_j++;
        }
        else
        {
            if(check_chain_indels(&(overlap_list->list[j].f_cigar), overlap_list->list[j].x_pos_s, 
            overlap_list->list[j].x_pos_e, indel_rate) == 1)
            {
                if(overlap_list->list[j].y_pos_strand == paf->buffer[inner_j].rev)
                {
                    x_overlapLen = Get_qe(paf->buffer[inner_j]) - Get_qs(paf->buffer[inner_j]) + 1;
                    y_overlapLen = Get_te(paf->buffer[inner_j]) - Get_ts(paf->buffer[inner_j]) + 1;
                    if(x_overlapLen < y_overlapLen) x_overlapLen = y_overlapLen;
                    x_overlapLen = x_overlapLen * 0.1;

                    // if(
                    // ((DIFF(overlap_list->list[j].x_pos_s, Get_qs(paf->buffer[inner_j])) < x_overlapLen)
                    // && (DIFF(overlap_list->list[j].x_pos_e, Get_qe(paf->buffer[inner_j])) < x_overlapLen))
                    // ||
                    // ((DIFF(overlap_list->list[j].y_pos_s, Get_ts(paf->buffer[inner_j])) < x_overlapLen)
                    // && (DIFF(overlap_list->list[j].y_pos_e, Get_te(paf->buffer[inner_j])) < x_overlapLen)))
                    if(
                    ((DIFF(overlap_list->list[j].x_pos_s, Get_qs(paf->buffer[inner_j])) < (uint64_t)x_overlapLen)
                    && (DIFF(overlap_list->list[j].x_pos_e, Get_qe(paf->buffer[inner_j])) < (uint64_t)x_overlapLen))
                    || 
                    ((DIFF(overlap_list->list[j].y_pos_s, Get_ts(paf->buffer[inner_j])) < (uint64_t)x_overlapLen)
                    && (DIFF(overlap_list->list[j].y_pos_e, Get_te(paf->buffer[inner_j])) < (uint64_t)x_overlapLen))
                    )
                    {
                        overlap_list->list[j].is_match = is_match;
                        overlap_list->list[j].strong = paf->buffer[inner_j].ml;
                        overlap_list->list[j].without_large_indel = paf->buffer[inner_j].no_l_indel;
                        if(is_exact == 1)
                        {
                            if(overlap_list->list[j].y_pos_strand == 0)
                            {
                                recover_UC_Read(overlap_read, &R_INF, overlap_list->list[j].y_id);
                            }
                            else
                            {
                                recover_UC_Read_RC(overlap_read, &R_INF, overlap_list->list[j].y_id);
                            }
                            if(if_exact_match(g_read->seq, g_read->length, overlap_read->seq, overlap_read->length, 
                            overlap_list->list[j].x_pos_s, overlap_list->list[j].x_pos_e, 
                            overlap_list->list[j].y_pos_s, overlap_list->list[j].y_pos_e))
                            {
                                overlap_list->list[j].shared_seed = 1;
                            }
                            else
                            {
                                overlap_list->list[j].shared_seed = 0;
                            }
                        }
                    }
                    else
                    {
                        overlap_list->list[j].is_match = 3;
                    }
                }
                else
                {
                    overlap_list->list[j].is_match = 3;
                }
            }

            j++;
            inner_j++;
        }
    }
}



void update_exact_overlaps(overlap_region_alloc* overlap_list, UC_Read* g_read, UC_Read* overlap_read)
{
    uint64_t j;
    for (j = 0; j < overlap_list->length; j++)
    {
        if (overlap_list->list[j].is_match != 1)
        {
            if((overlap_list->list[j].x_pos_e + 1 - overlap_list->list[j].x_pos_s) != 
               (overlap_list->list[j].y_pos_e + 1 - overlap_list->list[j].y_pos_s))
            {
                continue;
            }

            if(overlap_list->list[j].y_pos_strand == 0)
            {
                recover_UC_Read(overlap_read, &R_INF, overlap_list->list[j].y_id);
            }
            else
            {
                recover_UC_Read_RC(overlap_read, &R_INF, overlap_list->list[j].y_id);
            }

            if(if_exact_match(g_read->seq, g_read->length, overlap_read->seq, overlap_read->length, 
            overlap_list->list[j].x_pos_s, overlap_list->list[j].x_pos_e, 
            overlap_list->list[j].y_pos_s, overlap_list->list[j].y_pos_e))
            {
                overlap_list->list[j].is_match = 1;
                overlap_list->list[j].strong = 0;
                overlap_list->list[j].without_large_indel = 1;
                overlap_list->list[j].shared_seed = 1;
            }
        }
    }
}

void ha_print_ovlp_stat(ma_hit_t_alloc* paf, ma_hit_t_alloc* rev_paf, long long readNum)
{
	long long forward, reverse, strong, weak, exact, no_l_indel;
	long long i, j;

	no_l_indel = forward = reverse = exact = strong = weak = 0;
	for (i = 0; i < readNum; i++) {
		forward += paf[i].length;
		reverse += rev_paf[i].length;
		for (j = 0; j < paf[i].length; j++) {
			if (paf[i].buffer[j].el == 1) exact++;
			if (paf[i].buffer[j].ml == 1) strong++;
			if (paf[i].buffer[j].ml == 0) weak++;
			if (paf[i].buffer[j].no_l_indel == 1) no_l_indel++;
		}
	}
	fprintf(stderr, "[M::%s] # overlaps: %lld\n", __func__, forward);
	fprintf(stderr, "[M::%s] # strong overlaps: %lld\n", __func__, strong);
	fprintf(stderr, "[M::%s] # weak overlaps: %lld\n", __func__, weak);
	fprintf(stderr, "[M::%s] # exact overlaps: %lld\n", __func__, exact); // this seems not right
	fprintf(stderr, "[M::%s] # inexact overlaps: %lld\n", __func__, forward - exact);
	fprintf(stderr, "[M::%s] # overlaps without large indels: %lld\n", __func__, no_l_indel);
	fprintf(stderr, "[M::%s] # reverse overlaps: %lld\n", __func__, reverse);
}

void fill_chain(Fake_Cigar* chain, char* x_string, char* y_string, long long xBeg, long long yBeg,
long long x_readLen, long long y_readLen, Cigar_record* cigar, uint8_t* c2n)
{
    /**
    long long i, xOffset, yOffset, xRegionLen, yRegionLen, maxXpos, maxYpos, mapGlobalScore, mapExtentScore, zdroped;
    long long xBuoundaryScore, yBuoundaryScore;
    ///float band_rate = 0.08;
    int endbouns,mode;
    if(chain->length <= 0) return;

    kvec_t(uint8_t) x_num;
    kvec_t(uint8_t) y_num;
    kv_init(x_num);
    kv_init(y_num);

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
        kv_resize(uint8_t, x_num, (uint64_t)xRegionLen);
        kv_resize(uint8_t, y_num, (uint64_t)yRegionLen);
        ///text is x, query is y
        afine_gap_alignment(x_string, x_num.a, xRegionLen, y_string, y_num.a, yRegionLen, 
        c2n, BACKWARD_KSW, MATCH_SCORE_KSW, MISMATCH_SCORE_KSW, GAP_OPEN_KSW, GAP_EXT_KSW,
        BAND_KSW, Z_DROP_KSW, endbouns, &maxXpos, &maxYpos, &mapGlobalScore, 
        &mapExtentScore, &xBuoundaryScore, &yBuoundaryScore, &zdroped);

        // fprintf(stderr, "* xOffset: %lld, yOffset: %lld, xRegionLen: %lld, yRegionLen: %lld, bandLen: %lld, maxXpos: %lld, maxYpos: %lld, zdroped: %lld\n",
        // xOffset, yOffset, xRegionLen, yRegionLen, BAND_KSW, maxXpos, maxYpos, zdroped);
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

            // bandLen = MAX((MIN(xRegionLen, yRegionLen))*band_rate, DIFF(xRegionLen, yRegionLen));
            // if(bandLen == 0) bandLen = MIN(xRegionLen, yRegionLen);
        }


        ///do alignment forward
        kv_resize(uint8_t, x_num, (uint64_t)xRegionLen);
        kv_resize(uint8_t, y_num, (uint64_t)yRegionLen);
        ///text is x, query is y
        afine_gap_alignment(x_string+xOffset, x_num.a, xRegionLen, y_string+yOffset, y_num.a, yRegionLen, 
        c2n, FORWARD_KSW, MATCH_SCORE_KSW, MISMATCH_SCORE_KSW, GAP_OPEN_KSW, GAP_EXT_KSW,
        BAND_KSW, Z_DROP_KSW, endbouns, &maxXpos, &maxYpos, &mapGlobalScore, 
        &mapExtentScore, &xBuoundaryScore, &yBuoundaryScore, &zdroped);
        // fprintf(stderr, "# xOffset: %lld, yOffset: %lld, xRegionLen: %lld, yRegionLen: %lld, bandLen: %lld, maxXpos: %lld, maxYpos: %lld, zdroped: %lld\n",
        // xOffset, yOffset, xRegionLen, yRegionLen, BAND_KSW, maxXpos, maxYpos, zdroped);
    }


    kv_destroy(x_num);
    kv_destroy(y_num);
    **/
}
void Final_phasing(overlap_region_alloc* overlap_list, Cigar_record_alloc* cigarline,
UC_Read* g_read, UC_Read* overlap_read, uint8_t* c2n)
{
    uint64_t i, xLen, yStrand;
    char* x_string;
    char* y_string;
    Cigar_record* cigar;
    resize_Cigar_record_alloc(cigarline, overlap_list->length);



    for (i = 0; i < overlap_list->length; i++)
    {
        if(overlap_list->list[i].is_match == 1 ||
           overlap_list->list[i].is_match == 2 ||
           overlap_list->list[i].is_match == 3)
        {
            xLen = overlap_list->list[i].x_pos_e - overlap_list->list[i].x_pos_s + 1;
            yStrand = overlap_list->list[i].y_pos_strand;
            cigar = &(cigarline->buffer[i]);
            ///has already been matched exactly
            if(overlap_list->list[i].is_match == 1 && overlap_list->list[i].shared_seed == 1)
            {
                add_cigar_record(g_read->seq + overlap_list->list[i].x_pos_s, xLen, cigar, 0);
            }
            else 
            {
                if(yStrand == 0)
                {
                    recover_UC_Read(overlap_read, &R_INF, overlap_list->list[i].y_id);
                }
                else
                {
                    recover_UC_Read_RC(overlap_read, &R_INF, overlap_list->list[i].y_id);
                }
                x_string = g_read->seq;
                y_string = overlap_read->seq;

                fill_chain(&(overlap_list->list[i].f_cigar), x_string, y_string, 
                overlap_list->list[i].x_pos_s, overlap_list->list[i].y_pos_s,
                Get_READ_LENGTH(R_INF, overlap_list->list[i].x_id), 
                Get_READ_LENGTH(R_INF, overlap_list->list[i].y_id), cigar, c2n);
            }

        }
    }

}

static void worker_ov_final(void *data, long i, int tid)
{
	ha_ovec_buf_t *b = ((ha_ovec_buf_t**)data)[tid];

	//get_new_candidates(i, &g_read, &overlap_list, &array_list, &l, 0.001, 0);
    ha_get_candidates_interface(b->ab, i, &b->self_read, &b->olist, &b->olist_hp, &b->clist, 0.001, 
    asm_opt.max_n_chain, 0, &(b->k_flag), &b->r_buf, &(R_INF.paf[i]), &(R_INF.reverse_paf[i]), &(b->tmp_region), NULL, &(b->sp));

	overlap_region_sort_y_id(b->olist.list, b->olist.length);
	ma_hit_sort_tn(R_INF.paf[i].buffer, R_INF.paf[i].length);
	ma_hit_sort_tn(R_INF.reverse_paf[i].buffer, R_INF.reverse_paf[i].length);

	update_overlaps(&b->olist, &(R_INF.paf[i]), &b->self_read, &b->ovlp_read, 1, 1);
	update_overlaps(&b->olist, &(R_INF.reverse_paf[i]), &b->self_read, &b->ovlp_read, 2, 0);
	///recover missing exact overlaps
	update_exact_overlaps(&b->olist, &b->self_read, &b->ovlp_read);

	///Final_phasing(&overlap_list, &cigarline, &g_read, &overlap_read, c2n);
	push_final_overlaps(&(R_INF.paf[i]), R_INF.reverse_paf, &b->olist, 1);
	push_final_overlaps(&(R_INF.reverse_paf[i]), R_INF.reverse_paf, &b->olist, 2);
}




void reset_final_overlaps(overlap_region_alloc *overlap_list)
{
    
    uint64_t i;
    for (i = 0; i < overlap_list->length; i++)
    {
        if (overlap_list->list[i].is_match == 1 || overlap_list->list[i].is_match == 2)
        {
            overlap_list->list[i].x_pos_s = overlap_list->list[i].x_pos_e = (uint32_t)-1;
            overlap_list->list[i].y_pos_s = overlap_list->list[i].y_pos_e = (uint32_t)-1;
            overlap_list->list[i].is_match = 0;
        }
    }

    ha_sort_list_by_anchor(overlap_list);
}

void debug_affine_gap_alignment(overlap_region_alloc *overlap_list, UC_Read* g_read, UC_Read* overlap_read)
{
    uint64_t i;
    kvec_t(uint8_t) x_num;
    kvec_t(uint8_t) y_num;
    kv_init(x_num);
    kv_init(y_num);
    for (i = 0; i < overlap_list->length; i++)
    {
        if (overlap_list->list[i].is_match == 1 && overlap_list->list[i].shared_seed == 1)
        {

            kv_resize(uint8_t, x_num, (uint64_t)(Get_READ_LENGTH(R_INF, overlap_list->list[i].x_id)));
            kv_resize(uint8_t, y_num, (uint64_t)(Get_READ_LENGTH(R_INF, overlap_list->list[i].y_id)));

            get_affine_gap_score(&(overlap_list->list[i]), g_read, overlap_read, x_num.a, y_num.a,
            overlap_list->list[i].x_pos_e + 1 - overlap_list->list[i].x_pos_s, 
            overlap_list->list[i].y_pos_e + 1 - overlap_list->list[i].y_pos_s);
        }
    }

    kv_destroy(x_num);
    kv_destroy(y_num);
}

static void worker_ov_final_high_het(void *data, long i, int tid)
{
    ha_ovec_buf_t *b = ((ha_ovec_buf_t**)data)[tid];

    ha_get_candidates_interface(b->ab, i, &b->self_read, &b->olist, &b->olist_hp, &b->clist, HIGH_HET_ERROR_RATE, 
    asm_opt.max_n_chain, 1, &(b->k_flag), &b->r_buf, &(R_INF.paf[i]), &(R_INF.reverse_paf[i]), &(b->tmp_region), NULL, &(b->sp));

    overlap_region_sort_y_id(b->olist.list, b->olist.length);
    ma_hit_sort_tn(R_INF.paf[i].buffer, R_INF.paf[i].length);
    ma_hit_sort_tn(R_INF.reverse_paf[i].buffer, R_INF.reverse_paf[i].length);


    ///update_overlaps(&b->olist, &(R_INF.paf[i]), &b->self_read, &b->ovlp_read, 1, 1);
    update_overlaps_chain_width(&b->olist, &(R_INF.paf[i]), &b->self_read, &b->ovlp_read, 1, 1, 0.002);
    update_overlaps(&b->olist, &(R_INF.reverse_paf[i]), &b->self_read, &b->ovlp_read, 2, 0);
    ///recover missing exact overlaps
    update_exact_overlaps(&b->olist, &b->self_read, &b->ovlp_read);

    
    ///Final_phasing(&overlap_list, &cigarline, &g_read, &overlap_read, c2n);
    push_final_overlaps(&(R_INF.paf[i]), R_INF.reverse_paf, &b->olist, 1);
    push_final_overlaps(&(R_INF.reverse_paf[i]), R_INF.reverse_paf, &b->olist, 2);

    ///debug_affine_gap_alignment(&b->olist, &b->self_read, &b->ovlp_read);

    reset_final_overlaps(&b->olist);
    correct_overlap_high_het(&b->olist, &R_INF, &b->self_read, &b->correct, &b->ovlp_read);
    push_final_overlaps_increment(&(R_INF.reverse_paf[i]), R_INF.reverse_paf, &b->olist, 2);
}

void Output_PAF()
{
    fprintf(stderr, "Writing PAF to disk ...... \n");
    char* paf_name = (char*)malloc(strlen(asm_opt.output_file_name)+50);
    sprintf(paf_name, "%s.ovlp.paf", asm_opt.output_file_name);
    FILE* output_file = fopen(paf_name, "w");
    uint64_t i, j;
    ma_hit_t_alloc* sources = R_INF.paf;

    for (i = 0; i < R_INF.total_reads; i++)
    {
        for (j = 0; j < sources[i].length; j++)
        {
            fwrite(Get_NAME(R_INF, Get_qn(sources[i].buffer[j])), 1, 
            Get_NAME_LENGTH(R_INF, Get_qn(sources[i].buffer[j])), output_file);
            fwrite("\t", 1, 1, output_file);
            fprintf(output_file, "%lu\t", (unsigned long)Get_READ_LENGTH(R_INF, Get_qn(sources[i].buffer[j])));
            fprintf(output_file, "%d\t", Get_qs(sources[i].buffer[j]));
            fprintf(output_file, "%d\t", Get_qe(sources[i].buffer[j]));
            if(sources[i].buffer[j].rev)
            {
                fprintf(output_file, "-\t");
            }
            else
            {
                fprintf(output_file, "+\t");
            }
            fwrite(Get_NAME(R_INF, Get_tn(sources[i].buffer[j])), 1, 
            Get_NAME_LENGTH(R_INF, Get_tn(sources[i].buffer[j])), output_file);
            fwrite("\t", 1, 1, output_file);
            fprintf(output_file, "%lu\t", (unsigned long)Get_READ_LENGTH(R_INF, Get_tn(sources[i].buffer[j])));
            fprintf(output_file, "%d\t", Get_ts(sources[i].buffer[j]));
            fprintf(output_file, "%d\t", Get_te(sources[i].buffer[j]));
            fprintf(output_file, "%d\t", sources[i].buffer[j].ml);
            fprintf(output_file, "%d\t", sources[i].buffer[j].bl);
            fprintf(output_file, "255\n");

        }
    }

    free(paf_name);
    fclose(output_file);

    fprintf(stderr, "PAF has been written.\n");
}

int check_cluster(uint64_t* list, long long listLen, ma_hit_t_alloc* paf, float threshold)
{
    long long i, k;
    uint32_t qn, tn;
    long long T_edges, A_edges;
    T_edges = A_edges = 0;
    for (i = 0; i < listLen; i++)
    {
        qn = (uint32_t)list[i];
        for (k = i + 1; k < listLen; k++)
        {
            tn = (uint32_t)list[k];
            if(get_specific_overlap(&(paf[qn]), qn, tn) != -1)
            {
                A_edges++;
            }

            if(get_specific_overlap(&(paf[tn]), tn, qn) != -1)
            {
                A_edges++;
            }

            T_edges = T_edges + 2;
        }

    }

    if(A_edges >= (T_edges*threshold))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}


void rescue_edges(ma_hit_t_alloc* paf, ma_hit_t_alloc* rev_paf, 
long long readNum, long long rescue_threshold, float cluster_threshold)
{
    double startTime = Get_T();
    long long i, j, revises = 0;
    uint32_t qn, tn;

    kvec_t(uint64_t) edge_vector;
    kv_init(edge_vector);
    kvec_t(uint64_t) edge_vector_index;
    kv_init(edge_vector_index);
    uint64_t flag;
    int index;

    for (i = 0; i < readNum; i++)
    {
        edge_vector.n = 0;
        edge_vector_index.n = 0;
        for (j = 0; j < paf[i].length; j++)
        {
            qn = Get_qn(paf[i].buffer[j]);
            tn = Get_tn(paf[i].buffer[j]);
            index = get_specific_overlap(&(rev_paf[tn]), tn, qn);
            if(index != -1)
            {
                flag = tn;
                flag = flag << 32;
                flag = flag | (uint64_t)(index);
                kv_push(uint64_t, edge_vector, flag);
                kv_push(uint64_t, edge_vector_index, j);
            }
        }

        ///the read itself has these overlaps, but all related reads do not have 
        ///we need to remove all overlaps from paf[i], and then add all overlaps to rev_paf[i]
        if((long long)edge_vector.n >= rescue_threshold && 
        check_cluster(edge_vector.a, edge_vector.n, paf, cluster_threshold) == 1)
        {
            add_overlaps(&(paf[i]), &(rev_paf[i]), edge_vector_index.a, edge_vector_index.n);
            remove_overlaps(&(paf[i]), edge_vector_index.a, edge_vector_index.n);
            revises = revises + edge_vector.n;
        }

        edge_vector.n = 0;
        edge_vector_index.n = 0;
        for (j = 0; j < rev_paf[i].length; j++)
        {
            qn = Get_qn(rev_paf[i].buffer[j]);
            tn = Get_tn(rev_paf[i].buffer[j]);
            index = get_specific_overlap(&(paf[tn]), tn, qn);
            if(index != -1)
            {
                flag = tn;
                flag = flag << 32;
                flag = flag | (uint64_t)(index);
                kv_push(uint64_t, edge_vector, flag);
                kv_push(uint64_t, edge_vector_index, j);
            }
        }

        ///the read itself do not have these overlaps, but all related reads have
        ///we need to remove all overlaps from rev_paf[i], and then add all overlaps to paf[i]
        if((long long)edge_vector.n >= rescue_threshold && 
        check_cluster(edge_vector.a, edge_vector.n, paf, cluster_threshold) == 1)
        {
            remove_overlaps(&(rev_paf[i]), edge_vector_index.a, edge_vector_index.n);
            add_overlaps_from_different_sources(paf, &(paf[i]), edge_vector.a, edge_vector.n);
            revises = revises + edge_vector.n;
        }
    }


    kv_destroy(edge_vector);
    kv_destroy(edge_vector_index);

    fprintf(stderr, "[M::%s] took %0.2fs, revise edges #: %lld\n\n", __func__, Get_T()-startTime, revises);
}

void hap_recalculate_peaks(char* output_file_name)
{
    destory_read_bin(&R_INF);
    destory_ma_hit_t_alloc(R_INF.paf);
    destory_ma_hit_t_alloc(R_INF.reverse_paf);

    char* gfa_name = (char*)malloc(strlen(output_file_name)+25);
	sprintf(gfa_name, "%s.ec", output_file_name);

    int hom_cov, het_cov;
    // construct hash table for high occurrence k-mers
    if (!(asm_opt.flag & HA_F_NO_KMER_FLT)) {
        ha_flt_tab = ha_ft_gen(&asm_opt, &R_INF, &hom_cov, 0);
        ha_opt_update_cov(&asm_opt, hom_cov);
    }
    free(R_INF.read_length);
    free(R_INF.name_index);

    load_All_reads(&R_INF, gfa_name);

    ha_idx = ha_pt_gen(&asm_opt, ha_flt_tab, 1, 0, &R_INF, &hom_cov, &het_cov); // build the index
    asm_opt.hom_cov = hom_cov;
    asm_opt.het_cov = het_cov;
    ha_pt_destroy(ha_idx);
    ha_idx = 0;

    destory_read_bin(&R_INF);
    free(gfa_name);

    load_all_data_from_disk(&R_INF.paf, &R_INF.reverse_paf, asm_opt.output_file_name); 
    fprintf(stderr, "M::%s has done.\n", __func__);   
}

void ha_overlap_final(void)
{
	int i, hom_cov, het_cov;
	ha_ovec_buf_t **b;
    ha_flt_tab_hp = ha_idx_hp = NULL;

	CALLOC(b, asm_opt.thread_num);
	for (i = 0; i < asm_opt.thread_num; ++i)
		b[i] = ha_ovec_init(asm_opt.flag & HA_F_HIGH_HET, 1,0);///b[i] = ha_ovec_init(1, 1);
	ha_idx = ha_pt_gen(&asm_opt, ha_flt_tab, 1, 0, &R_INF, &hom_cov, &het_cov); // build the index
    if(asm_opt.flag & HA_F_HIGH_HET)
    {
        kt_for(asm_opt.thread_num, worker_ov_final_high_het, b, R_INF.total_reads);
    }
    else
    {
        kt_for(asm_opt.thread_num, worker_ov_final, b, R_INF.total_reads);
    }
    
	ha_pt_destroy(ha_idx);
	ha_idx = 0;
	for (i = 0; i < asm_opt.thread_num; ++i)
		ha_ovec_destroy(b[i]);
	free(b);
    asm_opt.hom_cov = hom_cov;
    asm_opt.het_cov = het_cov;
}

static void worker_ov_utg(void *data, long i, int tid)
{
	ha_ovec_buf_t *b = ((ha_ovec_buf_t**)data)[tid];
    if(b->ua->a[i].len == 0) return;

    ha_get_ug_candidates(b->ab, i, &(b->ua->a[i]), b->ua, &b->olist, &b->clist, 
    0.3, asm_opt.polyploidy*5, 0, &(b->k_flag), &b->r_buf, ha_flt_tab, ha_idx,
    &(b->tmp_region), NULL, /**0.3**/0);
    
	overlap_region_sort_y_id(b->olist.list, b->olist.length);
	ma_hit_sort_tn(R_INF.paf[i].buffer, R_INF.paf[i].length);
	ma_hit_sort_tn(R_INF.reverse_paf[i].buffer, R_INF.reverse_paf[i].length);

	update_overlaps(&b->olist, &(R_INF.paf[i]), &b->self_read, &b->ovlp_read, 1, 1);
	update_overlaps(&b->olist, &(R_INF.reverse_paf[i]), &b->self_read, &b->ovlp_read, 2, 0);
	///recover missing exact overlaps
	update_exact_overlaps(&b->olist, &b->self_read, &b->ovlp_read);

	///Final_phasing(&overlap_list, &cigarline, &g_read, &overlap_read, c2n);
	push_final_overlaps(&(R_INF.paf[i]), R_INF.reverse_paf, &b->olist, 1);
	push_final_overlaps(&(R_INF.reverse_paf[i]), R_INF.reverse_paf, &b->olist, 2);
}


void ug_idx_build(ma_ug_t *ug, int hap_n)
{
    int flag = asm_opt.flag&HA_F_NO_HPC, i;
    asm_opt.flag -= flag;
    // ha_flt_tab = ha_ft_ug_gen(&asm_opt, &(ug->u), hap_n, hap_n);
    // ha_idx = ha_pt_ug_gen(&asm_opt, ha_flt_tab, &(ug->u), hap_n);

    ha_ovec_buf_t **b = NULL;
    // overlap and correct reads
    CALLOC(b, asm_opt.thread_num);
    for (i = 0; i < asm_opt.thread_num; ++i)
    {
        b[i] = ha_ovec_init(1, 1, 0);
        b[i]->ua = &(ug->u);
    }
            
    kt_for(asm_opt.thread_num, worker_ov_utg, b, R_INF.total_reads);

    for (i = 0; i < asm_opt.thread_num; ++i)
        ha_ovec_destroy(b[i]);
    free(b);

    ha_ft_destroy(ha_flt_tab); 
    ha_pt_destroy(ha_idx);
    asm_opt.flag += flag;
    exit(1);
}

int ha_assemble(void)
{
    // debug_mc_g_t(MC_NAME);
    // debug_mc_gg_t(MC_NAME, 0, 0);
	extern void ha_extract_print_list(const All_reads *rs, int n_rounds, const char *o);
	int r, hom_cov = -1, ovlp_loaded = 0;
	if (asm_opt.load_index_from_disk && load_all_data_from_disk(&R_INF.paf, &R_INF.reverse_paf, asm_opt.output_file_name)) {
		ovlp_loaded = 1;
		fprintf(stderr, "[M::%s::%.3f*%.2f] ==> loaded corrected reads and overlaps from disk\n", __func__, yak_realtime(), yak_cpu_usage());
		if (asm_opt.extract_list) {
			ha_extract_print_list(&R_INF, asm_opt.extract_iter, asm_opt.extract_list);
			exit(0);
		}
		if (!(asm_opt.flag & HA_F_SKIP_TRIOBIN) && !(asm_opt.flag & HA_F_VERBOSE_GFA)) ha_triobin(&asm_opt);
        ///if (!(asm_opt.flag & HA_F_SKIP_TRIOBIN)) ha_triobin(&asm_opt), ovlp_loaded = 2;
		if (asm_opt.flag & HA_F_WRITE_EC) Output_corrected_reads();
		if (asm_opt.flag & HA_F_WRITE_PAF) Output_PAF();
        if (asm_opt.het_cov == -1024) hap_recalculate_peaks(asm_opt.output_file_name), ovlp_loaded = 2;
	}
	if (!ovlp_loaded) {
        ha_flt_tab = ha_idx = NULL;
        if((asm_opt.flag & HA_F_VERBOSE_GFA)) load_pt_index(&ha_flt_tab, &ha_idx, &R_INF, &asm_opt, asm_opt.output_file_name), load_ct_index(&ha_ct_table, asm_opt.output_file_name);

		// construct hash table for high occurrence k-mers
		if (!(asm_opt.flag & HA_F_NO_KMER_FLT) && ha_flt_tab == NULL) 
        {
			ha_flt_tab = ha_ft_gen(&asm_opt, &R_INF, &hom_cov, 0);
			ha_opt_update_cov(&asm_opt, hom_cov);
		}
		// error correction
		assert(asm_opt.number_of_round > 0);
		for (r = ha_idx?asm_opt.number_of_round-1:0; r < asm_opt.number_of_round; ++r) {
			ha_opt_reset_to_round(&asm_opt, r); // this update asm_opt.roundID and a few other fields
			ha_overlap_and_correct(r);
			fprintf(stderr, "[M::%s::%.3f*%.2f@%.3fGB] ==> corrected reads for round %d\n", __func__, yak_realtime(),
					yak_cpu_usage(), yak_peakrss_in_gb(), r + 1);
			fprintf(stderr, "[M::%s] # bases: %lld; # corrected bases: %lld; # recorrected bases: %lld\n", __func__,
					asm_opt.num_bases, asm_opt.num_corrected_bases, asm_opt.num_recorrected_bases);
			fprintf(stderr, "[M::%s] size of buffer: %.3fGB\n", __func__, asm_opt.mem_buf / 1073741824.0);
		}
		if (asm_opt.flag & HA_F_WRITE_EC) Output_corrected_reads();
		// overlap between corrected reads
		ha_opt_reset_to_round(&asm_opt, asm_opt.number_of_round);
		ha_overlap_final();
		fprintf(stderr, "[M::%s::%.3f*%.2f@%.3fGB] ==> found overlaps for the final round\n", __func__, yak_realtime(),
				yak_cpu_usage(), yak_peakrss_in_gb());
		ha_print_ovlp_stat(R_INF.paf, R_INF.reverse_paf, R_INF.total_reads);
		ha_ft_destroy(ha_flt_tab);
		if (asm_opt.flag & HA_F_WRITE_PAF) Output_PAF();
		ha_triobin(&asm_opt);
	}
    if(ovlp_loaded == 2) ovlp_loaded = 0;

    build_string_graph_without_clean(asm_opt.min_overlap_coverage, R_INF.paf, R_INF.reverse_paf, 
        R_INF.total_reads, R_INF.read_length, asm_opt.min_overlap_Len, asm_opt.max_hang_Len, asm_opt.clean_round, 
        asm_opt.gap_fuzz, asm_opt.min_drop_rate, asm_opt.max_drop_rate, asm_opt.output_file_name, asm_opt.large_pop_bubble_size, 0, !ovlp_loaded);
	destory_All_reads(&R_INF);
	return 0;
}
