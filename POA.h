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
    ///0是match，1是mismatch，2是x缺字符（y多字符），而3是y缺字符（x多字符）
    uint64_t weight;
    uint64_t num_insertions;
    ///这条路径上到backbone节点之前总共有多少节点
    uint64_t length;
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
    uint64_t weight;
    ///记录的是以当前节点为尾的deletion个数
    uint64_t num_insertions;
    char base;
    Edge_alloc mismatch_edges;
    Edge_alloc deletion_edges;
    Edge_alloc insertion_edges;

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
    ///has a indivial start node 0
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
void append_Edge_alloc(Edge_alloc* list,  uint64_t in_node, uint64_t out_node, uint64_t weight, uint64_t length);

void init_Node_alloc(Node_alloc* list);
void destory_Node_alloc(Node_alloc* list);
void clear_Node_alloc(Node_alloc* list);
uint64_t append_Node_alloc(Node_alloc* list, char base);
uint64_t* get_Topo_Sort_Order(Node_alloc* list, int need_sort);


void init_Graph(Graph* g);
void addUnmatchedSeqToGraph(Graph* g, char* g_read_seq, long long g_read_length, long long* startID, long long* endID);
void addmatchedSeqToGraph(Graph* backbone, long long currentNodeID, char* x_string, long long x_length, 
        char* y_string, long long y_length, CIGAR* cigar, long long backbone_start, long long backbone_end);
void destory_Graph(Graph* g);
void clear_Graph(Graph* g);
void Perform_POA(Graph* g, overlap_region_alloc* overlap_list, All_reads* R_INF, UC_Read* g_read);

void Graph_debug(Graph* backbone, long long currentNodeID, char* x_string, long long x_length, 
        char* y_string, long long y_length, CIGAR* cigar, long long backbone_start, long long backbone_end);

void debug_graph(Graph* g, long long backbone_length);

uint64_t inline add_Node_Graph(Graph* g, char base)
{
    return append_Node_alloc(&g->g_nodes, base);
}


///仅仅用于误配边
inline void add_mismatchEdge_weight(Graph* g, uint64_t in_node, char base, int last_operation)
{
    long long i = 0;
    long long nodeID;
    Edge_alloc* edge = &(g->g_nodes.list[in_node].mismatch_edges);

    for (i = 0; i < edge->length; i++)
    {
        nodeID = edge->list[i].out_node;
        if(g->g_nodes.list[nodeID].base == base)
        {
            edge->list[i].weight++;
            ///如果上一个操作是insertion
            if (last_operation == 2)
            {
                edge->list[i].num_insertions++;
            }
            
            break;
        }
    }

    ///说明不存在这么一条边
    if (i == edge->length)
    {
        nodeID  = add_Node_Graph(g, base);
        
        ///只有match边长度是0
        ///mismatch边长度都是1
        append_Edge_alloc(edge, in_node, nodeID, 1, 1);
        ///如果上一个操作是insertion
        if (last_operation == 2)
        {
            edge->list[edge->length - 1].num_insertions++;
        }

        ///将新节点的mismatch_edges连到backbone上
        append_Edge_alloc(&(g->g_nodes.list[nodeID].mismatch_edges), nodeID, in_node + 1, 1, 0);
    }
    ///获得节点的mismatch_edges长度为1，其他均为0

}



inline void add_single_deletionEdge_weight(Graph* g, long long alignNodeID, long long nextNodeID, uint64_t edge_length)
{
    long long i = 0;
    long long nodeID;
    Edge_alloc* edge = &(g->g_nodes.list[alignNodeID].deletion_edges);

    for (i = 0; i < edge->length; i++)
    {
        nodeID = edge->list[i].out_node;
        if(nodeID == nextNodeID)
        {
            edge->list[i].weight++;
            break;
        }
    }

    ///说明不存在这么一条边
    if (i == edge->length)
    {
        append_Edge_alloc(edge, alignNodeID, nextNodeID, 1, edge_length);
    }
}

inline void add_deletionEdge_weight(Graph* g, long long alignNodeID, long long deletion_length)
{
    if (deletion_length == 1)
    {
        add_single_deletionEdge_weight(g, alignNodeID, alignNodeID + 1, 0);
    }
    else if (deletion_length == 2)
    {
        add_single_deletionEdge_weight(g, alignNodeID, alignNodeID + 1, 0);
        add_single_deletionEdge_weight(g, alignNodeID, alignNodeID + 2, 0);
        add_single_deletionEdge_weight(g, alignNodeID + 1, alignNodeID + 2, 0);
    }
    else if (deletion_length > 2)
    {
        ///fprintf(stderr, "too long deletion!\n");
        add_single_deletionEdge_weight(g, alignNodeID, alignNodeID + deletion_length, 0);
    }
}


inline int getEdge(Graph* g, Edge_alloc* edge, uint64_t edge_length, char base)
{
    long long i = 0;
    long long nodeID;

    for (i = 0; i < edge->length; i++)
    {
        if (edge->list[i].length == edge_length)
        {
            nodeID = edge->list[i].out_node;
            if(g->g_nodes.list[nodeID].base == base)
            {
                return i;
            }
        }
    }
    
    return -1;
}


inline int get_insertion_Edges(Graph* g, Edge_alloc* edge, uint64_t edge_length, char* bases)
{
    long long i = 0;
    long long nodeID;
    long long edgeID;

    if (edge_length < 1)
    {
        return -1;
    }
    

    edgeID = getEdge(g, edge, edge_length, bases[0]);

    long long return_edgeID = edgeID;

    if(edgeID == -1)
    {
        return -1;
    }

    
    Edge_alloc* new_edge = edge;

    for (i = 1; i < edge_length; i++)
    {
        nodeID = new_edge->list[edgeID].out_node;
        new_edge = &(g->g_nodes.list[nodeID].insertion_edges);
        edgeID = getEdge(g, new_edge, edge_length - i, bases[i]);
        if(edgeID == -1)
        {
            return -1;
        }
    }

    return edgeID;
}



inline int create_insertion_Edges(Graph* g, long long alignNodeID, uint64_t edge_length, char* bases)
{
    long long i = 0;
    long long nodeID;
    ///最后应该连回原节点
    ///long long backboneID = alignNodeID + 1;
    long long backboneID = alignNodeID;


    if (edge_length < 1)
    {
        return -1;
    }


    nodeID  = add_Node_Graph(g, bases[0]);
    ///将新加入的节点通过insertion_edges接到alignNodeID上
    append_Edge_alloc(&(g->g_nodes.list[alignNodeID].insertion_edges), alignNodeID, nodeID, 1, edge_length);

    alignNodeID = nodeID;

    for (i = 1; i < edge_length; i++)
    {
        nodeID  = add_Node_Graph(g, bases[i]);
        ///将新加入的节点通过insertion_edges接到alignNodeID上
        append_Edge_alloc(&(g->g_nodes.list[alignNodeID].insertion_edges), alignNodeID, nodeID, 1, edge_length - i);
        alignNodeID = nodeID;        
    }

    append_Edge_alloc(&(g->g_nodes.list[alignNodeID].insertion_edges), alignNodeID, backboneID, 1, 0);

}

inline void add_insertionEdge_weight(Graph* g, long long alignNodeID, char* insert, long long insert_length)
{
    long long i = 0;
    long long nodeID;
    long long edgeID;
    Edge_alloc* edge = &(g->g_nodes.list[alignNodeID].insertion_edges);

    if (insert_length == 1)
    {
        edgeID = getEdge(g, edge, 1, insert[0]);
        if (edgeID != -1)
        {
            ///这条路均只有一个出度
            edge->list[edgeID].weight++;
        }
        else ///不存在这么一条边
        {
            nodeID  = add_Node_Graph(g, insert[0]);
            append_Edge_alloc(edge, alignNodeID, nodeID, 1, 1);
            ///将新加入的节点通过insertion_edges接回backbone上
            ///应该连回到原节点，而不是原节点的下一个节点
            ///append_Edge_alloc(&(g->g_nodes.list[nodeID].insertion_edges), nodeID, alignNodeID + 1, 1, 0);
            append_Edge_alloc(&(g->g_nodes.list[nodeID].insertion_edges), nodeID, alignNodeID, 1, 0);
        }
    }
    else if (insert_length == 2)
    {
        /*******************第0个字符********************* */
        edgeID = getEdge(g, edge, 1, insert[0]);
        if (edgeID != -1)
        {
            ///这条路均只有一个出度
            edge->list[edgeID].weight++;
        }
        else ///不存在这么一条边
        {
            nodeID  = add_Node_Graph(g, insert[0]);
            append_Edge_alloc(edge, alignNodeID, nodeID, 1, 1);
            ///将新加入的节点通过insertion_edges接回backbone上
            ///应该连回到原节点，而不是原节点的下一个节点
            ///append_Edge_alloc(&(g->g_nodes.list[nodeID].insertion_edges), nodeID, alignNodeID + 1, 1, 0);
            append_Edge_alloc(&(g->g_nodes.list[nodeID].insertion_edges), nodeID, alignNodeID, 1, 0);
        }
        /*******************第0个字符********************* */

        /*******************第1个字符********************* */
        if (insert[1] != insert[0])
        {
            edgeID = getEdge(g, edge, 1, insert[1]);
            if (edgeID != -1)
            {
                ///这条路均只有一个出度
                edge->list[edgeID].weight++;
            }
            else ///不存在这么一条边
            {
                nodeID  = add_Node_Graph(g, insert[1]);
                append_Edge_alloc(edge, alignNodeID, nodeID, 1, 1);
                ///将新加入的节点通过insertion_edges接回backbone上
                ///应该连回到原节点，而不是原节点的下一个节点
                ///append_Edge_alloc(&(g->g_nodes.list[nodeID].insertion_edges), nodeID, alignNodeID + 1, 1, 0);
                append_Edge_alloc(&(g->g_nodes.list[nodeID].insertion_edges), nodeID, alignNodeID, 1, 0);
            }
        }
        /*******************第1个字符********************* */

        /**********************两个字符******************* */
        
        edgeID = get_insertion_Edges(g, edge, 2, insert);
        if (edgeID != -1)
        {
            ///这条路均只有一个出度
            edge->list[edgeID].weight++;
        }
        else
        {
            create_insertion_Edges(g, alignNodeID, insert_length, insert);
        }
        
        /**********************两个字符******************* */
    }
    else if (insert_length > 2)
    {
        ////fprintf(stderr, "too long insertion\n");
        /*************************大于2个字符************************** */

        edgeID = get_insertion_Edges(g, edge, insert_length, insert);
        if (edgeID != -1)
        {
            ///这条路均只有一个出度
            edge->list[edgeID].weight++;
        }
        else
        {
            create_insertion_Edges(g, alignNodeID, insert_length, insert);
        }
    }
    
    
    
}


#endif