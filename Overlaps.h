#ifndef __OVERLAPS__
#define __OVERLAPS__
#include <stdio.h>
#include <stdint.h>
#include "kvec.h"
#include "kdq.h"

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


typedef struct { size_t n, m; ma_utg_t *a; } ma_utg_v;

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
} R_to_U;

void init_R_to_U(R_to_U* x, uint64_t len);
void destory_R_to_U(R_to_U* x);
void set_R_to_U(R_to_U* x, uint32_t rID, uint32_t uID, uint32_t is_Unitig, uint8_t* flag);
void get_R_to_U(R_to_U* x, uint32_t rID, uint32_t* uID, uint32_t* is_Unitig);
void transfor_R_to_U(R_to_U* x);
void debug_utg_graph(ma_ug_t *ug, asg_t* read_g, kvec_asg_arc_t_warp* edge, int require_equal_nv, int test_tangle);
int asg_pop_bubble_primary(asg_t *g, int max_dist);
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

void init_Edge_iter(asg_t* g, uint32_t v, asg_arc_t* new_edges, uint32_t new_edges_n, Edge_iter* x);
int get_arc_t(Edge_iter* x, asg_arc_t* get);
int asg_pop_bubble_primary_trio(ma_ug_t *ug, int max_dist, uint32_t positive_flag, uint32_t negative_flag);


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

inline uint32_t get_unitig_back(asg_t *sg, ma_ug_t *ug, uint32_t begNode, uint32_t* endNode, 
long long* nodeLen, long long* baseLen, buf_t* b)
{
	ma_utg_v* u = NULL; 
    uint32_t v = begNode, w, k;
    uint32_t kv;
	(*nodeLen) = (*baseLen) = 0;
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
		if(kv==0) return END_TIPS;
        if(kv>1) return MUL_OUTPUT;
        ///kv must be 1 here
        kv = get_real_length(sg, v, &w);
		///means reach the end of a unitig
        if(get_real_length(sg, w^1, NULL)!=1)
		{
			(*baseLen) += sg->seq[v>>1].len;
			return MUL_INPUT;
		}

		for (k = 0; k < asg_arc_n(sg, v); k++)
        {
			if(asg_arc_a(sg, v)[k].del) continue;
			///here is just one undeleted edge
			(*baseLen) += asg_arc_len(asg_arc_a(sg, v)[k]);
			break;
		}

        v = w;
        if(v == begNode) return LOOP;
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
#define DIFF_HAP_RATE 0.75
#define TRIO_DROP_THRES 0.9
#define TRIO_DROP_LENGTH_THRES 0.8
#define MAX_STOP_RATE 0.6
#define TANGLE_MISSED_THRES 0.6
///if ug == NULL, nsg should be equal to read_sg
inline uint32_t check_different_haps(asg_t *nsg, ma_ug_t *ug, asg_t *read_sg, 
uint32_t v_0, uint32_t v_1, ma_hit_t_alloc* reverse_sources, buf_t* b_0, buf_t* b_1, 
R_to_U* ruIndex, uint32_t min_edge_length, uint32_t stops_threshold)
{
	uint32_t vEnd, qn, tn, j, is_Unitig, uId;
	long long ELen_0, ELen_1, tmp, max_stop_nodeLen, max_stop_baseLen;

	b_0->b.n = b_1->b.n = 0;
	if(get_unitig(nsg, ug, v_0, &vEnd, &ELen_0, &tmp, &max_stop_nodeLen, &max_stop_baseLen, 
	stops_threshold, b_0) == LOOP)
	{
		return UNAVAILABLE;
	}
	if(get_unitig(nsg, ug, v_1, &vEnd, &ELen_1, &tmp, &max_stop_nodeLen, &max_stop_baseLen, 
	stops_threshold, b_1) == LOOP)
	{
		return UNAVAILABLE;
	}
	if(ELen_0<=min_edge_length || ELen_1<=min_edge_length) return UNAVAILABLE;

	rIdContig b_max, b_min;
	b_max.b_0 = b_min.b_0 = NULL;
	b_max.offset = b_max.readI = b_max.untigI = 0;
	b_min.offset = b_min.readI = b_min.untigI = 0;

	if(ELen_0<=ELen_1)
	{
        b_min.b_0 = b_0;
        b_max.b_0 = b_1;
    }
    else
    {
        b_min.b_0 = b_1;
        b_max.b_0 = b_0;
    }

	uint32_t max_count = 0, min_count = 0;
	ma_utg_t *node_min = NULL, *node_max = NULL;
	if(ug != NULL)
	{
		/*****************************label all unitigs****************************************/
		for (b_max.untigI = 0; b_max.untigI < b_max.b_0->b.n; b_max.untigI++)
		{
			node_max = &(ug->u.a[b_max.b_0->b.a[b_max.untigI]>>1]);
			///each read
			for (b_max.readI = 0; b_max.readI < node_max->n; b_max.readI++)
			{
				qn = (node_max->a[b_max.readI]>>33);
				set_R_to_U(ruIndex, qn, (b_max.b_0->b.a[b_max.untigI]>>1), 1, &(read_sg->seq[qn].c));
			}
		}
		/*****************************label all unitigs****************************************/

		///each unitig
		for (b_min.untigI = 0; b_min.untigI < b_min.b_0->b.n; b_min.untigI++)
		{
			
			node_min = &(ug->u.a[(b_min.b_0->b.a[b_min.untigI]>>1)]);

			///each read
			for (b_min.readI = 0; b_min.readI < node_min->n; b_min.readI++)
        	{
				qn = node_min->a[b_min.readI]>>33;

				/************************BUG: don't forget****************************/
				if(reverse_sources[qn].length > 0) min_count++;
				///if(reverse_sources[qn].length >= 0) min_count++;
				/************************BUG: don't forget****************************/
				for (j = 0; j < (long long)reverse_sources[qn].length; j++)
				{
					tn = Get_tn(reverse_sources[qn].buffer[j]);
					if(read_sg->seq[tn].del == 1)
					{
						get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
						if(tn == (uint32_t)-1 || is_Unitig == 1 || read_sg->seq[tn].del == 1) continue;
					}

					get_R_to_U(ruIndex, tn, &uId, &is_Unitig);
					if(uId!=(uint32_t)-1 && is_Unitig == 1)
					{
						max_count++;
						break;
					}
				}
			}
		}
		/*****************************label all unitigs****************************************/
		for (b_max.untigI = 0; b_max.untigI < b_max.b_0->b.n; b_max.untigI++)
		{
			node_max = &(ug->u.a[b_max.b_0->b.a[b_max.untigI]>>1]);
			///each read
			for (b_max.readI = 0; b_max.readI < node_max->n; b_max.readI++)
			{
				qn = (node_max->a[b_max.readI]>>33);
				ruIndex->index[qn] = (uint32_t)-1;
			}
		}
		/*****************************label all unitigs****************************************/
	}
	else
	{
		/*****************************label all reads****************************************/
		for (b_max.untigI = 0; b_max.untigI < b_max.b_0->b.n; b_max.untigI++)
		{
			qn = (b_max.b_0->b.a[b_max.untigI]>>1);
			set_R_to_U(ruIndex, qn, 1, 1, &(read_sg->seq[qn].c));
		}
		/*****************************label all reads****************************************/

		///each read
		for (b_min.untigI = 0; b_min.untigI < b_min.b_0->b.n; b_min.untigI++)
		{
			qn = (b_min.b_0->b.a[b_min.untigI]>>1);

			/************************BUG: don't forget****************************/
			if(reverse_sources[qn].length > 0) min_count++;
			///if(reverse_sources[qn].length >= 0) min_count++;
			/************************BUG: don't forget****************************/

			for (j = 0; j < (long long)reverse_sources[qn].length; j++)
        	{
            	tn = Get_tn(reverse_sources[qn].buffer[j]);
				if(nsg->seq[tn].del == 1)
				{
					get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
					if(tn == (uint32_t)-1 || is_Unitig == 1 || nsg->seq[tn].del == 1) continue;
				}


				get_R_to_U(ruIndex, tn, &uId, &is_Unitig);
				if(uId!=(uint32_t)-1 && is_Unitig == 1)
				{
					max_count++;
					break;
				}
			}
		}

		/*****************************label all reads****************************************/
		for (b_max.untigI = 0; b_max.untigI < b_max.b_0->b.n; b_max.untigI++)
		{
			qn = (b_max.b_0->b.a[b_max.untigI]>>1);
			ruIndex->index[qn] = (uint32_t)-1;
		}
		/*****************************label all reads****************************************/
	}

	// if(((v_0==7707) && (v_1==26867))||((v_1==7707) && (v_0==26867)))
	// {
	// 	fprintf(stderr, "******\nv_0>>1: %u, v_0&1: %u, ELen_0: %u\n", v_0>>1, v_0&1, (uint32_t)ELen_0);
	// 	fprintf(stderr, "v_1>>1: %u, v_1&1: %u, ELen_1: %u\n", v_1>>1, v_1&1, (uint32_t)ELen_1);
	// 	fprintf(stderr, "min_count: %u, max_count: %u, DIFF_HAP_RATE: %f\n\n", 
	// 	min_count, max_count, DIFF_HAP_RATE);
	// }
	
	if(min_count == 0) return UNAVAILABLE;
	if(max_count > min_count*DIFF_HAP_RATE) return PLOID;
	return NON_PLOID;
}

inline uint32_t check_different_haps_naive(asg_t *nsg, ma_ug_t *ug, asg_t *read_sg, 
uint32_t v_0, uint32_t v_1, ma_hit_t_alloc* reverse_sources, buf_t* b_0, buf_t* b_1, 
R_to_U* ruIndex, uint32_t min_edge_length, uint32_t stops_threshold)
{
	uint32_t vEnd, qn, tn, j, is_Unitig;
	long long ELen_0, ELen_1, tmp, max_stop_nodeLen, max_stop_baseLen;

	b_0->b.n = b_1->b.n = 0;
	if(get_unitig(nsg, ug, v_0, &vEnd, &ELen_0, &tmp, &max_stop_nodeLen, &max_stop_baseLen, 
	stops_threshold, b_0) == LOOP)
	{
		return UNAVAILABLE;
	} 
	if(get_unitig(nsg, ug, v_1, &vEnd, &ELen_1, &tmp, &max_stop_nodeLen, &max_stop_baseLen, 
	stops_threshold, b_1) == LOOP)
	{
		return UNAVAILABLE;
	}

	if(ELen_0<=min_edge_length || ELen_1<=min_edge_length) return UNAVAILABLE;

	rIdContig b_max, b_min;
	b_max.b_0 = b_min.b_0 = NULL;
	b_max.offset = b_max.readI = b_max.untigI = 0;
	b_min.offset = b_min.readI = b_min.untigI = 0;

	if(ELen_0<=ELen_1)
	{
        b_min.b_0 = b_0;
        b_max.b_0 = b_1;
    }
    else
    {
        b_min.b_0 = b_1;
        b_max.b_0 = b_0;
    }

	uint32_t max_count = 0, min_count = 0;
	ma_utg_t *node_min = NULL, *node_max = NULL;

	if(ug != NULL)
	{
		///each unitig
		for (b_min.untigI = 0; b_min.untigI < b_min.b_0->b.n; b_min.untigI++)
		{
			
			node_min = &(ug->u.a[(b_min.b_0->b.a[b_min.untigI]>>1)]);

			///each read
			for (b_min.readI = 0; b_min.readI < node_min->n; b_min.readI++)
        	{
				qn = node_min->a[b_min.readI]>>33;

				if(reverse_sources[qn].length > 0) min_count++;
				for (j = 0; j < (long long)reverse_sources[qn].length; j++)
				{
					tn = Get_tn(reverse_sources[qn].buffer[j]);
					if(read_sg->seq[tn].del == 1)
					{
						get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
						if(tn == (uint32_t)-1 || is_Unitig == 1 || read_sg->seq[tn].del == 1) continue;
					}

					///each unitig
					for (b_max.untigI = 0; b_max.untigI < b_max.b_0->b.n; b_max.untigI++)
					{
						node_max = &(ug->u.a[b_max.b_0->b.a[b_max.untigI]>>1]);
						///each read
						for (b_max.readI = 0; b_max.readI < node_max->n; b_max.readI++)
						{
							if(tn == (node_max->a[b_max.readI]>>33))
							{
								max_count++;
								goto end_check_different_haps_ug;
							}
						}
					}
				}

				end_check_different_haps_ug:;
			}
		}
	}
	else
	{
		///each read
		for (b_min.untigI = 0; b_min.untigI < b_min.b_0->b.n; b_min.untigI++)
		{
			qn = (b_min.b_0->b.a[b_min.untigI]>>1);

			if(reverse_sources[qn].length > 0) min_count++;

			for (j = 0; j < (long long)reverse_sources[qn].length; j++)
        	{
            	tn = Get_tn(reverse_sources[qn].buffer[j]);
				if(nsg->seq[tn].del == 1)
				{
					get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
					if(tn == (uint32_t)-1 || is_Unitig == 1 || nsg->seq[tn].del == 1) continue;
				}

				///each read
				for (b_max.untigI = 0; b_max.untigI < b_max.b_0->b.n; b_max.untigI++)
                {
					if((b_max.b_0->b.a[b_max.untigI]>>1) == tn)
					{
						max_count++;
						goto end_check_different_haps_non_ug;
					}
				}
			}

			end_check_different_haps_non_ug:;
		}
	}

	// if(((v_0==7707) && (v_1==26867))||((v_1==7707) && (v_0==26867)))
	// {
	// 	fprintf(stderr, "******\nv_0>>1: %u, v_0&1: %u, ELen_0: %u\n", v_0>>1, v_0&1, (uint32_t)ELen_0);
	// 	fprintf(stderr, "v_1>>1: %u, v_1&1: %u, ELen_1: %u\n", v_1>>1, v_1&1, (uint32_t)ELen_1);
	// 	fprintf(stderr, "min_count: %u, max_count: %u, DIFF_HAP_RATE: %f\n\n", 
	// 	min_count, max_count, DIFF_HAP_RATE);
	// }
	
	if(min_count == 0) return UNAVAILABLE;
	if(max_count > min_count*DIFF_HAP_RATE) return PLOID;
	return NON_PLOID;
}



typedef struct {
	uint32_t father_occ;
	uint32_t mother_occ;
	uint32_t ambig_occ;
	uint32_t drop_occ;
	uint32_t total;
} Trio_counter;

void resolve_tangles(ma_ug_t *src, asg_t *read_g, ma_hit_t_alloc* reverse_sources, long long minLongUntig, 
long long maxShortUntig, float l_untig_rate, float max_node_threshold, R_to_U* ruIndex, uint32_t trio_flag, 
float drop_ratio);
void adjust_utg_advance(asg_t *sg, ma_ug_t *ug, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex);
void rescue_contained_reads_aggressive(ma_ug_t *i_ug, asg_t *r_g, ma_hit_t_alloc* sources, ma_sub_t *coverage_cut,
R_to_U* ruIndex, int max_hang, int min_ovlp, long long bubble_dist, uint32_t chainLenThres, uint32_t is_bubble_check, 
uint32_t is_primary_check, kvec_asg_arc_t_warp* new_rtg_edges, kvec_t_u32_warp* new_rtg_nodes);
void rescue_missing_overlaps_aggressive(ma_ug_t *i_ug, asg_t *r_g, ma_hit_t_alloc* sources, ma_sub_t *coverage_cut,
R_to_U* ruIndex, int max_hang, int min_ovlp, long long bubble_dist, uint32_t is_bubble_check, 
uint32_t is_primary_check, kvec_asg_arc_t_warp* new_rtg_edges);
void all_to_all_deduplicate(ma_ug_t* ug, asg_t* read_g, ma_sub_t* coverage_cut, 
ma_hit_t_alloc* sources, uint8_t postive_flag, float drop_rate, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex, float double_check_rate);
void drop_semi_circle(ma_ug_t *ug, asg_t* nsg, asg_t* read_g, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex);
void rescue_wrong_overlaps_to_unitigs(ma_ug_t *i_ug, asg_t *r_g,  ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, 
ma_sub_t *coverage_cut, R_to_U* ruIndex, int max_hang, int min_ovlp, long long bubble_dist, kvec_asg_arc_t_warp* keep_edges);
void get_unitig_trio_flag(ma_utg_t* nsu, uint32_t flag, uint32_t* require, uint32_t* non_require, uint32_t* ambigious);
void rescue_missing_overlaps_backward(ma_ug_t *i_ug, asg_t *r_g, ma_hit_t_alloc* sources, ma_sub_t *coverage_cut,
R_to_U* ruIndex, int max_hang, int min_ovlp, long long bubble_dist, uint32_t backward_steps, 
uint32_t is_bubble_check, uint32_t is_primary_check);
uint32_t get_edge_from_source(ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, 
R_to_U* ruIndex, int max_hang, int min_ovlp, uint32_t query, uint32_t target, asg_arc_t* t);
uint64_t asg_bub_pop1_primary_trio(asg_t *g, ma_ug_t *utg, uint32_t v0, int max_dist, buf_t *b, 
uint32_t positive_flag, uint32_t negative_flag, uint32_t is_pop, uint64_t* path_base_len, uint64_t* path_nodes);
int unitig_arc_del_short_diploid_by_length(asg_t *g, float drop_ratio);


typedef struct{
    double weight;
    uint32_t uID:31, del:1;
    uint64_t dis;
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
	kvec_t(bed_in) bed;
    uint32_t* u_idx;
	uint64_t r_num;
} hc_links;

typedef struct{
    ///kvec_t(hc_edge) a;
    size_t n, m; 
    hc_edge *a;
}hc_edge_warp;

void init_hc_links(hc_links* link, uint64_t ug_num, uint64_t r_num);
void destory_hc_links(hc_links* link);
void clean_primary_untig_graph(ma_ug_t *ug, asg_t *read_g, ma_hit_t_alloc* reverse_sources,
long long bubble_dist, long long tipsLen, float tip_drop_ratio, long long stops_threshold, 
R_to_U* ruIndex, buf_t* b_0, uint8_t* visit, float density, uint32_t miniHapLen, 
uint32_t miniBiGraph, float chimeric_rate, int is_final_clean, int just_bubble_pop, 
float drop_ratio, hc_links* link);
void adjust_utg_by_primary(ma_ug_t **ug, asg_t* read_g, float drop_rate,
ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, ma_sub_t* coverage_cut, 
long long bubble_dist, long long tipsLen, float tip_drop_ratio, long long stops_threshold, 
R_to_U* ruIndex, float chimeric_rate, float drop_ratio, int max_hang, int min_ovlp,
kvec_asg_arc_t_warp* new_rtg_edges, hc_links* link);
void collect_reverse_unitigs(buf_t* b_0, buf_t* b_1, hc_links* link, ma_ug_t *ug, asg_t *read_sg);
ma_ug_t* copy_untig_graph(ma_ug_t *src);
ma_ug_t* output_trio_unitig_graph(asg_t *sg, ma_sub_t* coverage_cut, char* output_file_name, 
uint8_t flag, ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, long long bubble_dist, 
long long tipsLen, float tip_drop_ratio, long long stops_threshold, R_to_U* ruIndex, 
float chimeric_rate, float drop_ratio, int max_hang, int min_ovlp, int is_bench);
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

#define JUNK_COV 5
#define DISCARD_RATE 0.8

#endif
