#ifndef __OVERLAPS__
#define __OVERLAPS__
#include <stdint.h>
#include "kvec.h"
#include "kdq.h"
///#include "Hash_Table.h"

///#define MIN_OVERLAP_LEN 2000
///#define MIN_OVERLAP_LEN 500
///#define MIN_OVERLAP_LEN 50
#define MIN_OVERLAP_LEN 50
///#define MIN_OVERLAP_COVERAGE 1
#define MIN_OVERLAP_COVERAGE 0
#define MAX_HANG_LEN 1000
#define MAX_HANG_PRE 0.8
#define GAP_FUZZ 1000
#define MAX_SHORT_TIPS 3
#define MAX_BUBBLE_DIST 10000000
#define SMALL_BUBBLE_SIZE (uint32_t)-1
//#define SMALL_BUBBLE_SIZE 1000


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
void resize_ma_hit_t_alloc(ma_hit_t_alloc* x, uint64_t size);
void destory_ma_hit_t_alloc(ma_hit_t_alloc* x);
void add_ma_hit_t_alloc(ma_hit_t_alloc* x, ma_hit_t* element);
void ma_hit_sort_tn(ma_hit_t *a, long long n);
void ma_hit_sort_qns(ma_hit_t *a, long long n);

int load_all_data_from_disk(ma_hit_t_alloc **sources, ma_hit_t_alloc **reverse_sources, 
char* output_file_name);


typedef struct {
	ma_hit_t_alloc overlaps;
} Assembly_Graph;

void init_Assembly_Graph(Assembly_Graph* x);
void destory_Assembly_Graph(Assembly_Graph* x);
void collect_ma_hit_t(ma_hit_t_alloc* dest, ma_hit_t_alloc* sources, long long num_sources);
void normalize_ma_hit_t(ma_hit_t_alloc* sources, long long num_sources);
void debug_normalize_ma_hit_t(ma_hit_t_alloc* sources, long long num_sources);


typedef struct {
	uint32_t s:31, del:1, e;
	uint8_t c;
} ma_sub_t;

void ma_hit_sub(int min_dp, ma_hit_t_alloc* sources, long long n_read, uint64_t* readLen, 
long long mini_overlap_length, ma_sub_t** coverage_cut);
void ma_hit_cut(int min_dp, ma_hit_t_alloc* sources, long long n_read, uint64_t* readLen, 
long long mini_overlap_length, ma_sub_t** coverage_cut);
void ma_hit_flt(ma_hit_t_alloc* sources, long long n_read, const ma_sub_t *coverage_cut, 
int max_hang, int min_ovlp);
long long get_specific_overlap(ma_hit_t_alloc* x, uint32_t qn, uint32_t tn);

void debug_cut_ma_hit_t(ma_hit_t_alloc* sources, long long num_sources, ma_sub_t *coverage_cut);

typedef struct {
	uint64_t ul;
	uint32_t v;
	uint32_t ol:31, del:1;
	uint8_t strong;
	uint8_t el;
	uint8_t no_l_indel;
} asg_arc_t;

typedef struct {
	uint32_t len:31, del:1;
	uint8_t c;
} asg_seq_t;

typedef struct {
	uint32_t m_arc, n_arc:31, is_srt:1;
	asg_arc_t *arc;
	uint32_t m_seq, n_seq:31, is_symm:1;
	asg_seq_t *seq;
	uint64_t *idx;

	uint8_t* seq_vis;
} asg_t;

typedef struct { size_t n, m; uint64_t *a; } asg64_v;


#define MA_HT_INT        (-1)
#define MA_HT_QCONT      (-2)
#define MA_HT_TCONT      (-3)
#define MA_HT_SHORT_OVLP (-4)

///in default, max_hang = 1000, int_frac = 0.05, min_ovlp = 2000
static inline int ma_hit2arc(const ma_hit_t *h, int ql, int tl, int max_hang, float int_frac, int min_ovlp, asg_arc_t *p)
{
	int32_t tl5, tl3, ext5, ext3, qs = (int32_t)h->qns;
	uint32_t u, v, l; // u: query end; v: target end; l: length from u to v

	///if query and target are in different strand
	if (h->rev) tl5 = tl - h->te, tl3 = h->ts; // tl5: 5'-end overhang (on the query strand); tl3: similar
	else tl5 = h->ts, tl3 = tl - h->te;

	///ext5 and ext3 is the hang on left side and right side, respectively
	ext5 = qs < tl5? qs : tl5;
	ext3 = ql - h->qe < tl3? ql - h->qe : tl3;


	/**
	if (ext5 > max_hang || ext3 > max_hang || h->qe - qs < (h->qe - qs + ext5 + ext3) * int_frac)
		return MA_HT_INT;
	**/
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

	if (qs <= tl5 && ql - h->qe <= tl3) return MA_HT_QCONT; // query contained in target
	else if (qs >= tl5 && ql - h->qe >= tl3) return MA_HT_TCONT; // target contained in query 
	else if (qs > tl5) u = 0, v = !!h->rev, l = qs - tl5; ///u = 0 means query-to-target overlap, l is the length of node in string graph (not the overlap length)
	else u = 1, v = !h->rev, l = (ql - h->qe) - tl3; ///u = 1 means target-to-query overlaps, l is the length of node in string graph (not the overlap length)
	if (h->qe - qs + ext5 + ext3 < min_ovlp || h->te - h->ts + ext5 + ext3 < min_ovlp) return MA_HT_SHORT_OVLP; // short overlap
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


void build_string_graph(int min_dp, ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, long long n_read, uint64_t* readLen, 
long long mini_overlap_length, long long max_hang_length,
long long clean_round, float min_ovlp_drop_ratio, float max_ovlp_drop_ratio,
float final_ovlp_drop_ratio, char* output_file_name, long long bubble_dist);


#define asg_arc_len(arc) ((uint32_t)(arc).ul)
#define asg_arc_n(g, v) ((uint32_t)(g)->idx[(v)])
#define asg_arc_a(g, v) (&(g)->arc[(g)->idx[(v)]>>32])

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


typedef struct {
	uint32_t len:31, circ:1; // len: length of the unitig; circ: circular if non-zero
	uint32_t start, end; // start: starting vertex in the string graph; end: ending vertex
	uint32_t m, n; // number of reads
	uint64_t *a; // list of reads
	char *s; // unitig sequence is not null
} ma_utg_t;

typedef struct { size_t n, m; ma_utg_t *a; } ma_utg_v;

typedef struct {
	ma_utg_v u;
	asg_t *g;
} ma_ug_t;

typedef struct {
	uint32_t utg:31, ori:1, start, len;
} utg_intv_t;


/******************
 * Bubble popping *
 ******************/

typedef struct {
	uint32_t p; // the optimal parent vertex
	uint32_t d; // the shortest distance from the initial vertex
	uint32_t c; // max count of reads
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

// count the number of outgoing arcs, including reduced arcs
static inline int count_out_with_del(const asg_t *g, uint32_t v)
{
	uint32_t i, n, nv = asg_arc_n(g, v);
	const asg_arc_t *av = asg_arc_a(g, v);
	/**
	for (i = n = 0; i < nv; ++i)
		if (!av[i].del) ++n;
	return n;
	**/
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

void debug_info_of_specfic_read(char* name, ma_hit_t_alloc* sources, 
ma_hit_t_alloc* reverse_sources, int id, char* fun);

void build_string_graph_without_clean(int min_dp, ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, long long n_read, uint64_t* readLen, 
long long mini_overlap_length, long long max_hang_length,
long long clean_round, float min_ovlp_drop_ratio, float max_ovlp_drop_ratio, 
float corase_ovlp_drop_ratio, char* output_file_name, long long bubble_dist, int read_graph,
int write);

void debug_info_of_specfic_read(char* name, ma_hit_t_alloc* sources, 
ma_hit_t_alloc* reverse_sources, int id, char* command);

void collect_abnormal_edges(ma_hit_t_alloc* paf, ma_hit_t_alloc* rev_paf, long long readNum);

void add_overlaps(ma_hit_t_alloc* source_paf, ma_hit_t_alloc* dest_paf, uint64_t* source_index, long long listLen);
void remove_overlaps(ma_hit_t_alloc* source_paf, uint64_t* source_index, long long listLen);
void add_overlaps_from_different_sources(ma_hit_t_alloc* source_paf_list, ma_hit_t_alloc* dest_paf, 
uint64_t* source_index, long long listLen);
void print_revise_edges(ma_hit_t_alloc* source_paf, uint64_t* source_index, long long listLen);
#endif