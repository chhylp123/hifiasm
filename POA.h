#ifndef __POA_PARSER__
#define __POA_PARSER__
#include <stdint.h>
#include "Hash_Table.h"
#include "Process_Read.h"

/**
 1. 单个节点信息
    (1) ID
    (2) base
    (3) 入边信息
    (4) 出边信息
    (5) 比对到什么节点
 2. 各个节点信息，用数组下标组织，数组下标就是节点ID; 还要存拓扑排序后的下标和节点ID的对应关系
 3. 边
    (1) 边的起始
    (2) 边的结束节点
    (3) 过这条边的序列的label，也就是名称
 4. 各个序列信息 
    (1) 这个序列本身 
    (2) 这个序列的name或者ID 
    (3) 这个序列的在图中对应的起始和结束节点ID
 5. 有两个回溯矩阵，一个是graph的，一个是seq的
 **/

typedef struct
{
    uint64_t in_node;
    uint64_t out_node;
    uint64_t weight;
} Edge;

typedef struct
{
    Edge* list;
    uint64_t size;
    uint64_t length;
} Edge_alloc;

typedef struct
{
    uint64_t ID;
    char base;
    Edge_alloc income_edges;
    Edge_alloc outcome_edges;
    Edge_alloc alignedTo_Nodes;
} Node;

typedef struct
{
    uint64_t* list;
    uint8_t* visit;
    uint64_t size;
    uint64_t length;

    uint64_t* iterative_buffer;
    uint8_t* iterative_buffer_visit;
    uint64_t iterative_i;
} topo_Sorting_buffer;

typedef struct
{
    Node* list;
    topo_Sorting_buffer sort;
    uint64_t size;
    uint64_t length;
} Node_alloc;

typedef struct
{
    uint64_t g_n_nodes;
    uint64_t g_n_edges;
    uint64_t g_next_nodeID;
    Node_alloc g_nodes;



    char* seq;
    uint64_t seqID;
    uint64_t s_start_nodeID;
    uint64_t s_end_nodeID;
} Graph;


void init_Edge_alloc(Edge_alloc* list);
void clear_Edge_alloc(Edge_alloc* list);
void destory_Edge_alloc(Edge_alloc* list);
void append_Edge_alloc(Edge_alloc* list,  uint64_t in_node, uint64_t out_node, uint64_t weight);

void init_Node_alloc(Node_alloc* list);
void destory_Node_alloc(Node_alloc* list);
void clear_Node_alloc(Node_alloc* list);
uint64_t append_Node_alloc(Node_alloc* list, char base);
uint64_t* get_Topo_Sort_Order(Node_alloc* list, int need_sort);


void init_Graph(Graph* g);
void addUnmatchedSeqToGraph(Graph* g, char* g_read_seq, long long g_read_length, long long* startID, long long* endID);
void destory_Graph(Graph* g);
void clear_Graph(Graph* g);
void Perform_POA(Graph* g, overlap_region_alloc* overlap_list, All_reads* R_INF, UC_Read* g_read);


#endif