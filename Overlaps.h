#ifndef __OVERLAPS__
#define __OVERLAPS__
#include <stdio.h>
#include <stdint.h>
#include "kvec.h"
#include "kdq.h"
#include "ksort.h"

///#define MIN_OVERLAP_LEN 2000
///#define MIN_OVERLAP_LEN 500
///#define MIN_OVERLAP_LEN 50
///#define MIN_OVERLAP_LEN 50
///#define MIN_OVERLAP_COVERAGE 1
///#define MIN_OVERLAP_COVERAGE 0
///#define MAX_HANG_LEN 1000
///#define MAX_HANG_PRE 0.8
///#define GAP_FUZZ 1000
///#define MAX_SHORT_TIPS 3
///#define MAX_BUBBLE_DIST 10000000
#define SMALL_BUBBLE_SIZE (uint32_t)-1
//#define SMALL_BUBBLE_SIZE 1000
#define PRIMARY_LABLE 0
#define ALTER_LABLE 1
#define HAP_LABLE 2
#define FAKE_LABLE 4
#define TRIO_THRES 0.9
#define DOUBLE_CHECK_THRES 0.1
#define FINAL_DOUBLE_CHECK_THRES 0.2
#define CHIMERIC_TRIM_THRES 4
#define GAP_LEN 100
// #define PRIMARY_LABLE 1
// #define ALTER_LABLE 2
// #define HAP_LABLE 4


#define Get_qn(RECORD) ((uint32_t)((RECORD).qns>>32))
#define Get_qs(RECORD) ((uint32_t)((RECORD).qns))
#define Get_qe(RECORD) ((RECORD).qe)
#define Get_tn(RECORD) ((RECORD).tn)
#define Get_ts(RECORD) ((RECORD).ts)
#define Get_te(RECORD) ((RECORD).te)

#define LONG_TIPS 0
#define TWO_INPUT 1
#define TWO_OUTPUT 2
#define MUL_INPUT 3
#define MUL_OUTPUT 4
#define END_TIPS 5
#define LONG_TIPS_UNDER_MAX_EXT 6
#define LOOP 7

#define TRIM 10
#define CUT 11
#define CUT_DIF_HAP 12




///query is the read itself
typedef struct {
	uint64_t qns;
	uint32_t qe, tn, ts, te;
	uint32_t ml:31, rev:1;
	uint32_t bl:31, del:1;
	uint8_t el;
	uint8_t no_l_indel;
} ma_hit_t;

typedef struct {
	ma_hit_t* buffer;
    uint32_t size;
    uint32_t length;
	uint8_t is_fully_corrected;
	uint8_t is_abnormal;
} ma_hit_t_alloc;


void init_ma_hit_t_alloc(ma_hit_t_alloc* x);
void clear_ma_hit_t_alloc(ma_hit_t_alloc* x);
void resize_ma_hit_t_alloc(ma_hit_t_alloc* x, uint32_t size);
void destory_ma_hit_t_alloc(ma_hit_t_alloc* x);
void add_ma_hit_t_alloc(ma_hit_t_alloc* x, ma_hit_t* element);
void ma_hit_sort_tn(ma_hit_t *a, long long n);
void ma_hit_sort_qns(ma_hit_t *a, long long n);

int load_all_data_from_disk(ma_hit_t_alloc **sources, ma_hit_t_alloc **reverse_sources, 
char* output_file_name);



typedef struct {
	uint32_t s:31, del:1, e;
	uint8_t c;
} ma_sub_t;

void ma_hit_sub(int min_dp, ma_hit_t_alloc* sources, long long n_read, uint64_t* readLen, 
long long mini_overlap_length, ma_sub_t** coverage_cut);
void ma_hit_cut(ma_hit_t_alloc* sources, long long n_read, uint64_t* readLen, 
long long mini_overlap_length, ma_sub_t** coverage_cut);
void ma_hit_flt(ma_hit_t_alloc* sources, long long n_read, const ma_sub_t *coverage_cut, 
int max_hang, int min_ovlp);
long long get_specific_overlap(ma_hit_t_alloc* x, uint32_t qn, uint32_t tn);


typedef struct {
	uint32_t qSpre, qEpre, qScur, qEcur, qn;///[qSp, qEp) && [qSn, qEn]
	uint32_t tSpre, tEpre, tScur, tEcur, tn;
} u_trans_hit_t;

typedef struct {
	size_t n, m;
	u_trans_hit_t* a;
} kv_u_trans_hit_t;



typedef struct {
	uint32_t qs, qe, qn;
	uint32_t ts, te, tn;
	uint32_t occ;
	double nw;
	uint8_t f:6, rev:1, del:1;
	///uint8_t qo:4, to:4;
} u_trans_t;

typedef struct {
	size_t n, m;
	u_trans_t* a;
	kvec_t(uint64_t) idx;
} kv_u_trans_t;

#define u_trans_a(x, id) ((x).a + ((x).idx.a[(id)]>>32))
#define u_trans_n(x, id) ((uint32_t)((x).idx.a[(id)]))

typedef struct {
	uint64_t ul;
	uint32_t v;
	uint32_t ol:31, del:1;
	uint8_t strong;
	uint8_t el;
	uint8_t no_l_indel;
} asg_arc_t;

typedef struct {
	size_t n, m;
	asg_arc_t* a;
} kv_asg_arc_t;


typedef struct {
	uint32_t len:31, circ:1; // len: length of the unitig; circ: circular if non-zero
	uint32_t start, end; // start: starting vertex in the string graph; end: ending vertex
	uint32_t m, n; // number of reads
	uint64_t *a; // list of reads
	char *s; // unitig sequence is not null
} ma_utg_t;



typedef struct {
	uint32_t len:31, del:1;
	uint8_t c;
} asg_seq_t;

typedef struct {
	uint32_t m_arc, n_arc:31, is_srt:1;
	asg_arc_t *arc;
	uint32_t m_seq, n_seq:31, is_symm:1;
	uint32_t r_seq;

	asg_seq_t *seq;
	uint64_t *idx;

	uint8_t* seq_vis;

	uint32_t n_F_seq;
	ma_utg_t* F_seq;
} asg_t;

asg_t *asg_init(void);
void asg_destroy(asg_t *g);
void asg_arc_sort(asg_t *g);
void asg_seq_set(asg_t *g, int sid, int len, int del);
void asg_arc_index(asg_t *g);
void asg_cleanup(asg_t *g);
void asg_symm(asg_t *g);
void print_gfa(asg_t *g);


typedef struct { size_t n, m; uint64_t *a; } asg64_v;


typedef struct { size_t n, m; ma_utg_t *a;} ma_utg_v;

typedef struct {
	ma_utg_v u;
	asg_t *g;
	kvec_t(uint64_t) occ;
} ma_ug_t;

typedef struct {
	uint32_t utg:31, ori:1, start, len;
} utg_intv_t;


#define MA_HT_INT        (-1)
#define MA_HT_QCONT      (-2)
#define MA_HT_TCONT      (-3)
#define MA_HT_SHORT_OVLP (-4)

///in default, max_hang = 1000, int_frac = 0.8, min_ovlp = 50
static inline int ma_hit2arc(const ma_hit_t *h, int ql, int tl, int max_hang, float int_frac, int min_ovlp, asg_arc_t *p)
{
	int32_t tl5, tl3, ext5, ext3, qs = (int32_t)h->qns;
	uint32_t u, v, l; // u: query end; v: target end; l: length from u to v

	///if query and target are in different strand
	if (h->rev) tl5 = tl - h->te, tl3 = h->ts; // tl5: 5'-end overhang (on the query strand); tl3: similar
	else tl5 = h->ts, tl3 = tl - h->te;

	///ext5 and ext3 is the hang on left side and right side, respectively
	ext5 = qs < tl5? qs : tl5;
	ext3 = ql - (int)h->qe < tl3? ql - (int)h->qe : tl3;


	/**
	if (ext5 > max_hang || ext3 > max_hang || h->qe - qs < (h->qe - qs + ext5 + ext3) * int_frac)
		return MA_HT_INT;
	**/
	///ext3 and ext5 should be always 0
	if (ext5 > max_hang || ext3 > max_hang 
	|| h->qe - qs < (h->qe - qs + ext5 + ext3) * int_frac
	|| h->te - h->ts < (h->te - h->ts + ext5 + ext3) * int_frac)
	{
		return MA_HT_INT;
	}
		
	/**
	********************************query-to-target overlap****************************
	case 1: u = 0, rev = 0                           in the view of target: direction is 1 
	query: CCCCCCCCTAATTAAAAT                        target: TAATTAAAATGGGGGG (use ex-target as query)
	               ||||||||||         <--->                  ||||||||||
	       target: TAATTAAAATGGGGGG           query: CCCCCCCCTAATTAAAAT (use ex-query as target)

	case 2: u = 0, rev = 1                           in the view of target: direction is 0
	query: CCCCCCCCTAATTAAAAT					     target: CCCCCCATTTTAATTA  (use ex-target as query)
                   ||||||||||        <--->                         ||||||||||
	       target: TAATTAAAATGGGGGG                         query: ATTTTAATTAGGGGGGGG  (use ex-query as target)
	********************************query-to-target overlap****************************

	********************************target-to-query overlap****************************
	case 3: u = 1, rev = 0                           in the view of target: direction is 0
			 query: AAATAATATCCCCCCGCG                target: GGGCCGGCAAATAATAT (use ex-target as query)
					|||||||||          <--->                          |||||||||
	target: GGGCCGGCAAATAATAT                                  query: AAATAATATCCCCCCGCG (use ex-query as target)

	case 4: u = 1, rev = 1                          in the view of target: direction is 1
	         query: AAATAATATCCCCCCGCG                        target: ATATTATTTGCCGGCCC (use ex-target as query)
                    |||||||||             <--->                       |||||||||
	target: GGGCCGGCAAATAATAT                          query: CGCGGGGGATATTATTT (use ex-query as target)
	********************************target-to-query overlap****************************
	**/

	if (qs <= tl5 && ql - (int)h->qe <= tl3) return MA_HT_QCONT; // query contained in target
	else if (qs >= tl5 && ql - (int)h->qe >= tl3) return MA_HT_TCONT; // target contained in query 
	else if (qs > tl5) u = 0, v = !!h->rev, l = qs - tl5; ///u = 0 means query-to-target overlap, l is the length of node in string graph (not the overlap length)
	else u = 1, v = !h->rev, l = (ql - h->qe) - tl3; ///u = 1 means target-to-query overlaps, l is the length of node in string graph (not the overlap length)
	if ((int)h->qe - qs + ext5 + ext3 < min_ovlp || (int)h->te - (int)h->ts + ext5 + ext3 < min_ovlp) return MA_HT_SHORT_OVLP; // short overlap
	///u = 0 / 1 means query-to-target / target-to-query overlaps, 
	///l is the length of node in string graph (not the overlap length between two reads)
	u |= h->qns>>32<<1, v |= h->tn<<1;
	/**
	p->ul: |____________31__________|__________1___________|______________32_____________|
	                    qn            direction of overlap       length of this node (not overlap length)
						                (in the view of query)
	p->v : |___________31___________|__________1___________|
				        tn             reverse direction of overlap
						              (in the view of target)
	p->ol: overlap length
	**/
	p->ul = (uint64_t)u<<32 | l, p->v = v, p->ol = ql - l, p->del = 0;
	///l is the length of node in string graph (not the overlap length)

	p->strong = h->ml;
	p->el = h->el;
	p->no_l_indel = h->no_l_indel;
	return l;
}



#define asg_arc_len(arc) ((uint32_t)(arc).ul)
#define asg_arc_n(g, v) ((uint32_t)(g)->idx[(v)])
#define asg_arc_a(g, v) (&(g)->arc[(g)->idx[(v)]>>32])

static inline uint32_t asg_get_arc(asg_t *g, uint32_t v, uint32_t w, asg_arc_t* t)
{
	uint32_t i, nv = asg_arc_n(g, v);
	asg_arc_t *av = asg_arc_a(g, v);
	for (i = 0; i < nv; ++i)
	{
		if(av[i].del) continue;
		if(av[i].v == w)
		{
			(*t) = av[i];
			return 1;
		}
	}

	return 0;
}

// append an arc
static inline asg_arc_t *asg_arc_pushp(asg_t *g)
{
	if (g->n_arc == g->m_arc) {
		g->m_arc = g->m_arc? g->m_arc<<1 : 16;
		g->arc = (asg_arc_t*)realloc(g->arc, g->m_arc * sizeof(asg_arc_t));
	}
	return &g->arc[g->n_arc++];
}

// set asg_arc_t::del for v->w
static inline void asg_arc_del(asg_t *g, uint32_t v, uint32_t w, int del)
{
	uint32_t i, nv = asg_arc_n(g, v);
	asg_arc_t *av = asg_arc_a(g, v);
	for (i = 0; i < nv; ++i)
		if (av[i].v == w) av[i].del = !!del;
}

// set asg_arc_t::del and asg_seq_t::del to 1 for sequence s and all its associated arcs
static inline void asg_seq_del(asg_t *g, uint32_t s)
{
	uint32_t k;
	g->seq[s].del = 1;
	for (k = 0; k < 2; ++k) {
		uint32_t i, v = s<<1 | k;
		uint32_t nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		for (i = 0; i < nv; ++i) {
			av[i].del = 1;
			asg_arc_del(g, av[i].v^1, v^1, 1);
		}
	}
}

static inline void asg_seq_drop(asg_t *g, uint32_t s)
{
	///s is not at primary
	if(g->seq[s].c == ALTER_LABLE)
	{
		uint32_t k;
		for (k = 0; k < 2; ++k) 
		{
			///two directions of this node
			uint32_t i, v = s<<1 | k;
			uint32_t nv = asg_arc_n(g, v);
			asg_arc_t *av = asg_arc_a(g, v);
			for (i = 0; i < nv; ++i) 
			{
				if(av[i].del) continue;
				///if output node is at primary
				/****************************may have hap bugs********************************/
				///if(g->seq[(av[i].v>>1)].c == PRIMARY_LABLE)
				///if(g->seq[(av[i].v>>1)].c == PRIMARY_LABLE || g->seq[(av[i].v>>1)].c == HAP_LABLE)
				if(g->seq[(av[i].v>>1)].c != ALTER_LABLE)
				{/****************************may have hap bugs********************************/
					av[i].del = 1;
					asg_arc_del(g, av[i].v^1, v^1, 1);
				}
			}
		}
	}
}





/******************
 * Bubble popping *
 ******************/

typedef struct {
	uint32_t p; // the optimal parent vertex
	uint32_t d; // the shortest distance from the initial vertex
	uint32_t c; // max count of positive reads
	uint32_t m; // max count of negative reads
	uint32_t np; // max count of non-positive reads
	uint32_t nc; // max count of reads, no matter positive or negative
	uint32_t r:31, s:1; // r: the number of remaining incoming arc; s: state
	//s: state, s=0, this edge has not been visited, otherwise, s=1
} binfo_t;

typedef struct {
	///all information for each node
	binfo_t *a;
	kvec_t(uint32_t) S; // set of vertices without parents, nodes with all incoming edges visited
	kvec_t(uint32_t) T; // set of tips
	kvec_t(uint32_t) b; // visited vertices
	kvec_t(uint32_t) e; // visited edges/arcs
} buf_t;

typedef struct {
	kvec_t(uint64_t) Nodes; 
	kvec_t(uint64_t) Edges; 
	uint32_t pre_n_seq, seqID; 
} C_graph;

typedef struct {
	kvec_t(uint8_t) a;
	uint32_t i;
} kvec_t_u8_warp;

typedef struct {
	kvec_t(uint32_t) a;
	uint32_t i;
} kvec_t_u32_warp;

typedef struct {
	kvec_t(int32_t) a;
	uint32_t i;
} kvec_t_i32_warp;

typedef struct {
	kvec_t(uint64_t) a;
	uint64_t i;
} kvec_t_u64_warp;

typedef struct {
	kvec_t(asg_arc_t) a;
	uint64_t i;
}kvec_asg_arc_t_warp;

void sort_kvec_t_u64_warp(kvec_t_u64_warp* u_vecs, uint32_t is_descend);
int asg_arc_del_multi(asg_t *g);
int asg_arc_del_asymm(asg_t *g);

typedef struct {
	uint32_t q_pos;
	uint32_t t_pos;
	uint32_t t_id;
	uint32_t is_color;
} Hap_Align;

typedef struct {
	kvec_t(Hap_Align) x;
	uint64_t i;
} Hap_Align_warp;

typedef struct {
	buf_t* b_0;
    uint32_t untigI;
    uint32_t readI;
    uint32_t offset;
} rIdContig;

// count the number of outgoing arcs, including reduced arcs
static inline int count_out_with_del(const asg_t *g, uint32_t v)
{
	uint32_t nv = asg_arc_n(g, v);
	return nv;
}


// count the number of outgoing arcs, including reduced arcs
static inline int count_out_without_del(const asg_t *g, uint32_t v)
{
	uint32_t i, n, nv = asg_arc_n(g, v);
	const asg_arc_t *av = asg_arc_a(g, v);

	for (i = n = 0; i < nv; ++i)
		if (!av[i].del) ++n;
	return n;
}


void build_string_graph_without_clean(
int min_dp, ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, 
long long n_read, uint64_t* readLen, long long mini_overlap_length, 
long long max_hang_length, long long clean_round, long long gap_fuzz,
float min_ovlp_drop_ratio, float max_ovlp_drop_ratio, char* output_file_name, 
long long bubble_dist, int read_graph, int write);

void debug_info_of_specfic_read(char* name, ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, int id, char* command);
void collect_abnormal_edges(ma_hit_t_alloc* paf, ma_hit_t_alloc* rev_paf, long long readNum);
void add_overlaps(ma_hit_t_alloc* source_paf, ma_hit_t_alloc* dest_paf, uint64_t* source_index, long long listLen);
void remove_overlaps(ma_hit_t_alloc* source_paf, uint64_t* source_index, long long listLen);
void add_overlaps_from_different_sources(ma_hit_t_alloc* source_paf_list, ma_hit_t_alloc* dest_paf, 
uint64_t* source_index, long long listLen);

#define EvaluateLen(U, id) ((U).a[(id)].start)
#define IsMerge(U, id) ((U).a[(id)].end)
#define kv_reuse(v, rn, rm, r) ((v).n = (rn), (v).m = (rm), (v).a = (r))
#define Get_vis(visit, v, d) (((visit)[(v)>>1])&(((((v)<<(d))&1)+1)))
#define Set_vis(visit, v, d) (((visit)[(v)>>1])|=(((((v)<<(d))&1)+1)))





typedef struct {
	uint64_t len;
	uint32_t* index;
	uint8_t* is_het;
} R_to_U;

void init_R_to_U(R_to_U* x, uint64_t len);
void destory_R_to_U(R_to_U* x);
void set_R_to_U(R_to_U* x, uint32_t rID, uint32_t uID, uint32_t is_Unitig, uint8_t* flag);
void get_R_to_U(R_to_U* x, uint32_t rID, uint32_t* uID, uint32_t* is_Unitig);
void transfor_R_to_U(R_to_U* x);
void debug_utg_graph(ma_ug_t *ug, asg_t* read_g, kvec_asg_arc_t_warp* edge, int require_equal_nv, int test_tangle);
long long asg_arc_del_simple_circle_untig(ma_hit_t_alloc* sources, ma_sub_t* coverage_cut, asg_t *g, long long circleLen, int is_drop);

typedef struct {
	asg_t* g;

	asg_arc_t *av;
	uint32_t nv;
	uint32_t av_i;

	asg_arc_t* new_edges;
	uint32_t new_edges_n;
	uint32_t new_edges_i;
} Edge_iter;

typedef struct {
	asg_arc_t x;
    uint64_t Off;
    uint64_t weight;
}asg_arc_t_offset;

typedef struct {
	kvec_t(asg_arc_t_offset) a;
	uint64_t i;
}kvec_asg_arc_t_offset;



void init_Edge_iter(asg_t* g, uint32_t v, asg_arc_t* new_edges, uint32_t new_edges_n, Edge_iter* x);
int get_arc_t(Edge_iter* x, asg_arc_t* get);


inline int get_real_length(asg_t *g, uint32_t v, uint32_t* v_s)
{
    uint32_t i, kv = 0;
    for (i = 0, kv = 0; i < asg_arc_n(g, v); i++)
    {
        if(!asg_arc_a(g, v)[i].del)
        {
            if(v_s) v_s[kv] = asg_arc_a(g, v)[i].v;
            kv++;
        }
    }

    return kv;
}

inline uint32_t check_tip(asg_t *sg, uint32_t begNode, uint32_t* endNode, buf_t* b, uint32_t max_ext)
{
    ///cut tip of length <= max_ext
    uint32_t v = begNode, w;
    uint32_t kv;
    uint32_t eLen = 0;
    (*endNode) = (uint32_t)-1;
    b->b.n = 0;
    while (1)
    {
        kv = get_real_length(sg, v, NULL);
        (*endNode) = v;
        eLen++;
        if(b) kv_push(uint32_t, b->b, v);
        if(kv == 0) return END_TIPS;
        if(kv > 1) return MUL_OUTPUT;
        ///if(eLen > max_ext) return LONG_TIPS;
        ///kv must be 1 here
        kv = get_real_length(sg, v, &w);
        ///here this value must be >= 1
        if(get_real_length(sg, w^1, NULL)!=1) return MUL_INPUT;
        v = w;
        if(v == begNode) return LOOP;
        if(eLen >= max_ext) return LONG_TIPS;
    }
}

inline uint32_t get_unitig(asg_t *sg, ma_ug_t *ug, uint32_t begNode, uint32_t* endNode, 
long long* nodeLen, long long* baseLen, long long* max_stop_nodeLen, long long* max_stop_baseLen, 
uint32_t stops_threshold, buf_t* b)
{
	ma_utg_v* u = NULL; 
    uint32_t v = begNode, w, k;
    uint32_t kv, return_flag, n_stops = 0;
	long long pre_baseLen = 0, pre_nodeLen = 0;
	long long cur_baseLen = 0, cur_nodeLen = 0;
	(*max_stop_nodeLen) = (*max_stop_baseLen) = (*nodeLen) = (*baseLen) = 0;
    (*endNode) = (uint32_t)-1;
	if(ug!=NULL) u = &(ug->u);

    while (1)
    {
        kv = get_real_length(sg, v, NULL);
        (*endNode) = v;
		if(u == NULL) 
		{
			(*nodeLen)++;
		}
		else
		{
			(*nodeLen) += EvaluateLen((*u), v>>1);
		}
        if(b) kv_push(uint32_t, b->b, v);
		///means reach the end of a unitig
		if(kv!=1) (*baseLen) += sg->seq[v>>1].len;
		if(kv==0)
		{
			return_flag = END_TIPS;
			break;
			///return END_TIPS;
		} 
		if(kv>1)
		{
			return_flag = MUL_OUTPUT;
			break;
			///return MUL_OUTPUT;
		}
        ///kv must be 1 here
        kv = get_real_length(sg, v, &w);
		///means reach the end of a unitig
        if(get_real_length(sg, w^1, NULL)!=1)
		{

			n_stops++;
			if(n_stops >= stops_threshold)
			{
				(*baseLen) += sg->seq[v>>1].len;
				return_flag = MUL_INPUT;
				break;
				///return MUL_INPUT;
			}
			else
			{
				for (k = 0; k < asg_arc_n(sg, v); k++)
				{
					if(asg_arc_a(sg, v)[k].del) continue;
					///here is just one undeleted edge
					(*baseLen) += asg_arc_len(asg_arc_a(sg, v)[k]);
					break;
				}
			}

			cur_baseLen = (*baseLen) - pre_baseLen;
			pre_baseLen = (*baseLen);
			if(cur_baseLen > (*max_stop_baseLen))
			{
				(*max_stop_baseLen) = cur_baseLen;
			}


			cur_nodeLen = (*nodeLen) - pre_nodeLen;
			pre_nodeLen = (*nodeLen);
			if(cur_nodeLen > (*max_stop_nodeLen))
			{
				(*max_stop_nodeLen) = cur_nodeLen;
			}
		}
		else
		{
			for (k = 0; k < asg_arc_n(sg, v); k++)
			{
				if(asg_arc_a(sg, v)[k].del) continue;
				///here is just one undeleted edge
				(*baseLen) += asg_arc_len(asg_arc_a(sg, v)[k]);
				break;
			}
		}
		

        v = w;
        if(v == begNode)
		{
			return_flag = LOOP;
			break;
			///return LOOP;
		} 
    }




	cur_baseLen = (*baseLen) - pre_baseLen;
	pre_baseLen = (*baseLen);
	if(cur_baseLen > (*max_stop_baseLen))
	{
		(*max_stop_baseLen) = cur_baseLen;
	}


	cur_nodeLen = (*nodeLen) - pre_nodeLen;
	pre_nodeLen = (*nodeLen);
	if(cur_nodeLen > (*max_stop_nodeLen))
	{
		(*max_stop_nodeLen) = cur_nodeLen;
	}

	return return_flag;
}

#define UNAVAILABLE (uint32_t)-1
#define PLOID 0
#define NON_PLOID 1
// #define DIFF_HAP_RATE 0.75
#define TRIO_DROP_THRES 0.9
#define TRIO_DROP_LENGTH_THRES 0.8
#define MAX_STOP_RATE 0.6
#define TANGLE_MISSED_THRES 0.6
#define HET_HOM_RATE 0.7

typedef struct {
	uint32_t father_occ;
	uint32_t mother_occ;
	uint32_t ambig_occ;
	uint32_t drop_occ;
	uint32_t total;
} Trio_counter;

typedef struct {
	uint32_t p; // the optimal parent vertex
	uint32_t d; // the shortest distance from the initial vertex
	uint32_t r:31, s:1; // r: the number of remaining incoming arc; s: state
} binfo_s_t;

typedef struct {
	///all information for each node
	binfo_s_t *a;
	kvec_t(uint32_t) S; // set of vertices without parents, nodes with all incoming edges visited
	kvec_t(uint32_t) b; // visited vertices
	kvec_t(uint32_t) e; // visited edges/arcs
} buf_s_t;

typedef struct{
    buf_s_t *b;
    uint32_t n_thres, n_reads;
    asg_t *g;
    uint32_t check_cross;
	uint64_t bub_dist;
} bub_label_t;

void resolve_tangles(ma_ug_t *src, asg_t *read_g, ma_hit_t_alloc* reverse_sources, long long minLongUntig, 
long long maxShortUntig, float l_untig_rate, float max_node_threshold, R_to_U* ruIndex, uint8_t* is_r_het, 
uint32_t trio_flag, float drop_ratio);
void adjust_utg_advance(asg_t *sg, ma_ug_t *ug, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex, bub_label_t* b_mask_t, uint8_t* is_r_het);
void rescue_contained_reads_aggressive(ma_ug_t *i_ug, asg_t *r_g, ma_hit_t_alloc* sources, ma_sub_t *coverage_cut,
R_to_U* ruIndex, int max_hang, int min_ovlp, uint32_t chainLenThres, uint32_t is_bubble_check, 
uint32_t is_primary_check, kvec_asg_arc_t_warp* new_rtg_edges, kvec_t_u32_warp* new_rtg_nodes, bub_label_t* b_mask_t);
void rescue_missing_overlaps_aggressive(ma_ug_t *i_ug, asg_t *r_g, ma_hit_t_alloc* sources, ma_sub_t *coverage_cut,
R_to_U* ruIndex, int max_hang, int min_ovlp, uint32_t is_bubble_check, uint32_t is_primary_check, kvec_asg_arc_t_warp* new_rtg_edges, bub_label_t* b_mask_t);
void all_to_all_deduplicate(ma_ug_t* ug, asg_t* read_g, ma_sub_t* coverage_cut, 
ma_hit_t_alloc* sources, uint8_t postive_flag, float drop_rate, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex, uint8_t* is_r_het, float double_check_rate, int non_tig_occ);
void drop_semi_circle(ma_ug_t *ug, asg_t* nsg, asg_t* read_g, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex, uint8_t* is_r_het);
void rescue_wrong_overlaps_to_unitigs(ma_ug_t *i_ug, asg_t *r_g,  ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, 
ma_sub_t *coverage_cut, R_to_U* ruIndex, int max_hang, int min_ovlp, long long bubble_dist, kvec_asg_arc_t_warp* keep_edges, bub_label_t* b_mask_t);
void get_unitig_trio_flag(ma_utg_t* nsu, uint32_t flag, uint32_t* require, uint32_t* non_require, uint32_t* ambigious);
void rescue_missing_overlaps_backward(ma_ug_t *i_ug, asg_t *r_g, ma_hit_t_alloc* sources, ma_sub_t *coverage_cut,
R_to_U* ruIndex, int max_hang, int min_ovlp, uint32_t backward_steps, uint32_t is_bubble_check, uint32_t is_primary_check, bub_label_t* b_mask_t);
uint32_t get_edge_from_source(ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, 
R_to_U* ruIndex, int max_hang, int min_ovlp, uint32_t query, uint32_t target, asg_arc_t* t);
int unitig_arc_del_short_diploid_by_length(asg_t *g, float drop_ratio);
void asg_bub_backtrack_primary(asg_t *g, uint32_t v0, buf_t *b);

typedef struct{
    double weight;
    uint32_t uID;
    uint64_t dis;
	uint8_t is_cc:7, del:1;
	uint64_t occ;
	///uint64_t occ:63, scaff:1;
    ///uint32_t enzyme;
} hc_edge;

typedef struct{
    kvec_t(hc_edge) e;
    kvec_t(hc_edge) f;//forbiden
} hc_linkeage;

typedef struct{
	uint64_t beg, end;
}bed_interval;

typedef struct{
	size_t n, m;
	bed_interval* a;
}bed_in;

typedef struct{
    kvec_t(hc_linkeage) a;
    kvec_t(uint64_t) enzymes;
} hc_links;

#define N_HET 0
#define C_HET 1
#define P_HET 2
#define S_HET 4

typedef struct {
	uint32_t p_x_p, p_y_p, p_x, p_y;
	uint32_t c_x_p, c_y_p;
	uint8_t c_rev;
} ca_buf_t;

typedef struct {
	size_t n, m;
	ca_buf_t* a;
} kv_ca_buf_t;

typedef struct {
	kvec_t(uint32_t) uIDs;
	kvec_t(uint32_t) iDXs;
	uint32_t chain_num;
} sub_tran_t;

typedef struct{
	uint32_t* rUidx;
	uint64_t* rUpos;
	uint8_t* ir_het;
	uint32_t r_num, u_num;
	kvec_t(bed_in) bed;
	kvec_t(uint32_t) topo_buf;
	kvec_t(uint32_t) topo_res;
	buf_t b_buf_0, b_buf_1;
	///uint32_t* uLen;
	kv_u_trans_t k_trans;
	kv_u_trans_hit_t k_t_b;
	kv_ca_buf_t c_buf;
	sub_tran_t st;
}trans_chain;

typedef struct {
    uint32_t n;
    uint32_t* cov;
    uint64_t* pos_idx;
    ma_hit_t_alloc* reverse_sources;
    ma_sub_t *coverage_cut;
    R_to_U* ruIndex;
	asg_t *read_g;
    int max_hang;
    int min_ovlp;
    kvec_asg_arc_t_offset u_buffer; 
    kvec_t_i32_warp tailIndex;
    kvec_t_i32_warp prevIndex;
	uint8_t* is_r_het;
	trans_chain* t_ch;
}hap_cov_t;

typedef struct{
    ///kvec_t(hc_edge) a;
    size_t n, m; 
    hc_edge *a;
}hc_edge_warp;

typedef struct {
	uint32_t qs, qe, qn, qus, que;
	uint32_t ts, te, tn, tus, tue;
} utg_thit_t;

typedef struct {
	size_t n, m;
	utg_thit_t* a;
} kv_utg_thit_t_t;

typedef struct {
	ma_hit_t_alloc* reverse_sources;
    ma_sub_t *coverage_cut;
    R_to_U* ruIndex;
    asg_t *read_g;	
	kvec_asg_arc_t_offset u_buffer; 
    kvec_t_i32_warp tailIndex;
    kvec_t_i32_warp prevIndex;
	kv_utg_thit_t_t k_t_b;
	kv_ca_buf_t c_buf;
	kv_u_trans_t k_trans;
	uint64_t *pos_idx, rn;
	kvec_t(uint32_t) topo_res;
	ma_ug_t *cug;
	int max_hang;
	int min_ovlp;

	ma_utg_v u;
	kv_u_trans_t t;
	buf_t b0, b1;
} utg_trans_t;

typedef struct {
	ma_ug_t *ug;
	kvec_t(uint64_t) idx;
	kvec_t(uint32_t) dst;
} spg_t;

void init_hc_links(hc_links* link, uint64_t ug_num, trans_chain* t_ch);
void destory_hc_links(hc_links* link);
uint64_t get_bub_pop_max_dist(asg_t *g, buf_t *b);
uint64_t get_bub_pop_max_dist_advance(asg_t *g, buf_t *b);
uint64_t asg_bub_pop1_primary_trio(asg_t *g, ma_ug_t *utg, uint32_t v0, uint64_t max_dist, buf_t *b, uint32_t positive_flag, 
uint32_t negative_flag, uint32_t is_pop, uint64_t* path_base_len, uint64_t* path_nodes, hap_cov_t *cov, uint32_t is_update_chain, uint32_t keep_d, utg_trans_t *o);

void adjust_utg_by_primary(ma_ug_t **ug, asg_t* read_g, float drop_rate,
ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, ma_sub_t* coverage_cut, 
long long tipsLen, float tip_drop_ratio, long long stops_threshold, 
R_to_U* ruIndex, float chimeric_rate, float drop_ratio, int max_hang, int min_ovlp,
kvec_asg_arc_t_warp* new_rtg_edges, hap_cov_t **i_cov, bub_label_t* b_mask_t, uint32_t collect_p_trans, uint32_t collect_p_trans_f);
ma_ug_t* copy_untig_graph(ma_ug_t *src);
ma_ug_t* output_trio_unitig_graph(asg_t *sg, ma_sub_t* coverage_cut, char* output_file_name, 
uint8_t flag, ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, 
long long tipsLen, float tip_drop_ratio, long long stops_threshold, R_to_U* ruIndex, 
float chimeric_rate, float drop_ratio, int max_hang, int min_ovlp, int is_bench, bub_label_t* b_mask_t);
asg_t* copy_read_graph(asg_t *src);
ma_ug_t *ma_ug_gen(asg_t *g);
void ma_ug_destroy(ma_ug_t *ug);

inline int inter_interval(int a_s, int a_e, int b_s, int b_e, int* i_s, int* i_e)
{
    if(a_s > b_e || b_s > a_e) return 0;
    if(i_s) (*i_s) = a_s >= b_s? a_s : b_s; ///MAX(a_s, b_s);
    if(i_e) (*i_e) = a_e <= b_e? a_e : b_e; ///MIN(a_e, b_e);
    return 1;
}

inline uint32_t get_origin_uid(uint32_t v, trans_chain* t_ch, uint32_t *off, uint32_t *idx)
{
	if(off) (*off) = (t_ch->rUpos[v>>1]>>32);
	if(idx) (*idx) = (uint32_t)(t_ch->rUpos[v>>1]);
    if(t_ch->rUpos[v>>1] == (uint64_t)-1) return (uint32_t)-1;
    return (uint32_t)(((t_ch->rUidx[v>>1]>>1)<<1) + ((t_ch->rUidx[v>>1]^v)&1));
}
void chain_origin_trans_uid_by_distance(hap_cov_t *cov, asg_t *read_sg, 
uint32_t *pri_a, uint32_t pri_n, uint32_t pri_beg, uint64_t *i_pri_len, 
uint32_t *aux_a, uint32_t aux_n, uint32_t aux_beg, uint64_t *i_aux_len,
ma_ug_t *ug, uint32_t flag, double overall_score, const char* cmd);
int asg_arc_del_trans(asg_t *g, int fuzz);
void kt_u_trans_t_idx(kv_u_trans_t *ta, uint32_t n);
void kt_u_trans_t_simple_symm(kv_u_trans_t *ta, uint32_t un, uint32_t symm_add);
uint32_t get_u_trans_spec(kv_u_trans_t *ta, uint32_t qn, uint32_t tn, u_trans_t **r_a, uint32_t *occ);
int ma_ug_seq(ma_ug_t *g, asg_t *read_g, ma_sub_t *coverage_cut, ma_hit_t_alloc* sources, 
kvec_asg_arc_t_warp* edge, int max_hang, int min_ovlp, kvec_asg_arc_t_warp *E, uint32_t is_polish);


typedef struct{
    ma_sub_t* coverage_cut;
    ma_hit_t_alloc* sources; 
    ma_hit_t_alloc* reverse_sources; 
    long long tipsLen; 
    float tip_drop_ratio; 
    long long stops_threshold; 
    R_to_U* ruIndex; 
    float chimeric_rate;
    float drop_ratio; 
    int max_hang;
    int min_ovlp;
    int is_bench;
	long long gap_fuzz;
    bub_label_t* b_mask_t;
}ug_opt_t;

void adjust_utg_by_trio(ma_ug_t **ug, asg_t* read_g, uint8_t flag, float drop_rate,
ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, ma_sub_t* coverage_cut, 
long long tipsLen, float tip_drop_ratio, long long stops_threshold, 
R_to_U* ruIndex, float chimeric_rate, float drop_ratio, int max_hang, int min_ovlp,
kvec_asg_arc_t_warp* new_rtg_edges, bub_label_t* b_mask_t);
uint32_t cmp_untig_graph(ma_ug_t *src, ma_ug_t *dest);
void reduce_hamming_error(asg_t *sg, ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, 
int max_hang, int min_ovlp, long long gap_fuzz);
int ma_ug_seq_scaffold(ma_ug_t *g, asg_t *read_g, ma_sub_t *coverage_cut, ma_hit_t_alloc* sources, 
kvec_asg_arc_t_warp* edge, int max_hang, int min_ovlp, kvec_asg_arc_t_warp *E, uint32_t is_polish);
void ma_ug_print(const ma_ug_t *ug, asg_t* read_g, const ma_sub_t *coverage_cut, 
ma_hit_t_alloc* sources, R_to_U* ruIndex, const char* prefix, FILE *fp);
void ma_ug_print_simple(const ma_ug_t *ug, asg_t* read_g, const ma_sub_t *coverage_cut, 
ma_hit_t_alloc* sources, R_to_U* ruIndex, const char* prefix, FILE *fp);
trans_chain* init_trans_chain(ma_ug_t *ug, uint64_t r_num);
void destory_trans_chain(trans_chain **x);

typedef struct {///[cBeg, cEnd)
	uint32_t u_i, r_i, len, s_pos_cur, s_pre_v, s_pre_w, p_v, p_idx, p_uId, cBeg, cEnd;
    ///buf_t* x;
    uint32_t *a, an;
    ma_ug_t *ug;
    asg_t *read_sg;
    trans_chain* t_ch;
} u_trans_hit_idx;
void reset_u_trans_hit_idx(u_trans_hit_idx *t, uint32_t* i_x_a, uint32_t i_x_n, ma_ug_t *i_ug, 
asg_t *i_read_sg, trans_chain* i_t_ch, uint32_t i_cBeg, uint32_t i_cEnd);
uint32_t get_u_trans_hit(u_trans_hit_idx *t, u_trans_hit_t *hit);
inline uint32_t get_offset_adjust(uint32_t offset, uint32_t offsetLen, uint32_t targetLen)
{
    return ((double)(offset)/(double)(offsetLen))*targetLen;
}

uint32_t set_utg_offset(uint32_t *a, uint32_t a_n, ma_ug_t *ug, asg_t *read_sg, uint64_t* pos_idx, uint32_t is_clear,
uint32_t only_len);
uint64_t get_utg_cov(ma_ug_t *ug, uint32_t uID, asg_t* read_g, 
const ma_sub_t* coverage_cut, ma_hit_t_alloc* sources, R_to_U* ruIndex, uint8_t* r_flag);
trans_chain* load_hc_trans(const char *fn);
char *get_outfile_name(char* output_file_name);
void reset_u_trans_hit_idx(u_trans_hit_idx *t, uint32_t* i_x_a, uint32_t i_x_n, ma_ug_t *i_ug, 
asg_t *i_read_sg, trans_chain* i_t_ch, uint32_t i_cBeg, uint32_t i_cEnd);
void extract_sub_overlaps(uint32_t i_tScur, uint32_t i_tEcur, uint32_t i_tSpre, uint32_t i_tEpre,
uint32_t tn, kv_u_trans_hit_t* ktb, uint32_t bn);
void clean_u_trans_t_idx(kv_u_trans_t *ta, ma_ug_t *ug, asg_t *read_g);


#define JUNK_COV 5
#define DISCARD_RATE 0.8

#endif
