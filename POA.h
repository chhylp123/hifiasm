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
    long long beg;
    ///end is the index of next input data, instead of the index of last data
    long long end;
    long long length;
    long long size;
    long long* buffer;
} Queue;

inline void init_Queue(Queue* q)
{
    q->beg = 0;
    q->end = 0;
    q->length = 0;
    q->size = 20;
    q->buffer = (long long*)malloc(sizeof(long long)*q->size);
}

inline void clear_Queue(Queue* q)
{
    q->beg = 0;
    q->end = 0;
    q->length = 0;
}

inline void destory_Queue(Queue* q)
{
    free(q->buffer);
}

inline int is_empty_Queue(Queue* q)
{
    ///end is the index of next input data, instead of the index of last data
    if(q->beg == q->end)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

inline int is_full_Queue(Queue* q)
{
    ///end is the index of next input data, instead of the index of last data
    if(q->end < q->size)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

inline void push_to_Queue(Queue* q, long long nodeID)
{
    if(is_full_Queue(q))
    {
        long long move_length = q->beg;
        ///end is the index of next input data, instead of the index of last data
        long long current_length = q->end - q->beg;
        ///recalloc directly
        if(move_length == 0)
        {
            q->size = q->size * 2;
            q->buffer = (long long*)realloc(q->buffer, q->size*sizeof(long long));     
        }
        else
        {
            ///won't overlap
            if(current_length <= move_length)
            {
                memcpy(q->buffer, q->buffer+q->beg, sizeof(long long)*current_length);
            }
            else///may overlap
            {
                memmove(q->buffer, q->buffer+q->beg, sizeof(long long)*current_length);
            }

            q->beg = 0;
            q->end = current_length;
        }
    }
    

    q->buffer[q->end] = nodeID;
    q->end++;
}

inline int pop_from_Queue(Queue* q, long long* nodeID)
{
    if(is_empty_Queue(q))
    {
        (*nodeID) = -1;
        return 0;
    }
    else
    {
        (*nodeID) = q->buffer[q->beg];
        q->beg++;

        return 1;
    }
}





typedef struct
{
    uint64_t in_node;
    uint64_t out_node;
    ///0是match，1是mismatch，2是x缺字符（y多字符），而3是y缺字符（x多字符）
    uint64_t weight;
    uint64_t num_insertions;
    ///这条路径上到backbone节点之前总共有多少节点
    uint64_t length;
    uint64_t self_edge_ID;
    uint64_t reverse_edge_ID;
} Edge;

typedef struct
{
    Edge* list;
    uint64_t size;
    uint64_t length;
    uint64_t delete_length;
} Edge_alloc;

#define Real_Length(X) ((X).length - (X).delete_length)
#define Input_Edges(Node) ((Node).insertion_edges)
#define Output_Edges(Node) ((Node).deletion_edges)
#define G_Node(G, Node) ((G).g_nodes.list[(Node)])
#define If_Node_Exist(Node) ((Node).base != 'D')
#define If_Edge_Exist(E) ((E).out_node != (uint64_t)-1)
#define Visit(E) (E).length

typedef struct
{
    long long index;
} RSet;

inline void clear_RSet(RSet* set)
{
    set->index = 0;
}

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
    uint64_t delete_length;
} Node_alloc;

typedef struct
{
    uint64_t g_n_nodes;
    uint64_t g_n_edges;
    uint64_t g_next_nodeID;
    Node_alloc g_nodes;


    Queue node_q;
    char* seq;
    uint64_t seqID;
    uint64_t s_start_nodeID;
    uint64_t s_end_nodeID;
} Graph;

int add_and_check_bi_direction_edge(Graph* graph, Node* in_node, Node* out_node, uint64_t weight, uint64_t flag);
void add_bi_direction_edge(Graph* graph, Node* in_node, Node* out_node, uint64_t weight, uint64_t flag);
int remove_and_check_bi_direction_edge_from_nodes(Graph* graph, Node* in_node, Node* out_node);
int remove_and_check_bi_direction_edge_from_edge(Graph* graph, Edge* e);



inline int Pop_Node(Graph* DAGCon, Node** node)
{
    long long nodeID = 0;
    int return_flag = pop_from_Queue(&(DAGCon->node_q), &nodeID);

    (*node) = &(G_Node(*DAGCon, nodeID));
    return return_flag;
}

inline int Push_Node(Graph* DAGCon, Node** node)
{
    push_to_Queue(&(DAGCon->node_q), (**node).ID);
}

inline int getInputNodes(RSet* set, Graph* graph, Node* node, Node** get_Node)
{
    if(set->index >= Input_Edges(*node).length)
    {
        return 0;
    }

    ///skip all deleted edges
    while (
        set->index < Input_Edges(*node).length
        &&
        !(If_Edge_Exist(Input_Edges(*node).list[set->index]))
    )
    {
        set->index++;
    }


    if(
        set->index < Input_Edges(*node).length 
        && 
        If_Edge_Exist(Input_Edges(*node).list[set->index])
    )
    {
        (*get_Node) = &(G_Node((*graph), Input_Edges(*node).list[set->index].in_node));
        set->index++;
        return 1;
    }
    else
    {
        return 0;
    }    
}



inline int getInputEdges(RSet* set, Graph* graph, Node* node, Edge** get_Edge)
{
    if(set->index >= Input_Edges(*node).length)
    {
        return 0;
    }

    ///skip all deleted edges
    while (
        set->index < Input_Edges(*node).length
        &&
        !(If_Edge_Exist(Input_Edges(*node).list[set->index]))
    )
    {
        set->index++;
    }


    if(
        set->index < Input_Edges(*node).length 
        && 
        If_Edge_Exist(Input_Edges(*node).list[set->index])
    )
    {
        (*get_Edge) = &(Input_Edges(*node).list[set->index]);
        set->index++;
        return 1;
    }
    else
    {
        return 0;
    }    
}


inline int getOutputNodes(RSet* set, Graph* graph, Node* node, Node** get_Node)
{
    if(set->index >= Output_Edges(*node).length)
    {
        return 0;
    }

    ///skip all deleted edges
    while (
    set->index < Output_Edges(*node).length 
    && 
    !(If_Edge_Exist(Output_Edges(*node).list[set->index]))
    )
    {
        set->index++;
    }

    if(set->index < Output_Edges(*node).length && 
       If_Edge_Exist(Output_Edges(*node).list[set->index]))
    {
        (*get_Node) = &(G_Node((*graph), Output_Edges(*node).list[set->index].out_node));
        set->index++;
        return 1;
    }
    else
    {
        return 0;
    }
}



inline int getOutputEdges(RSet* set, Graph* graph, Node* node, Edge** get_Edge)
{
    if(set->index >= Output_Edges(*node).length)
    {
        return 0;
    }

    ///skip all deleted edges
    while (
    set->index < Output_Edges(*node).length 
    && 
    !(If_Edge_Exist(Output_Edges(*node).list[set->index]))
    )
    {
        set->index++;
    }

    if(set->index < Output_Edges(*node).length && 
       If_Edge_Exist(Output_Edges(*node).list[set->index]))
    {
        (*get_Edge) = &(Output_Edges(*node).list[set->index]);
        set->index++;
        return 1;
    }
    else
    {
        return 0;
    }
}


inline void get_bi_direction_edges(Graph* DAGCon, Edge* edge, Edge** e_forward, Edge** e_backward)
{
    long long in_node = edge->in_node;
    long long out_node = edge->out_node;

    if(
    edge->self_edge_ID < Output_Edges(G_Node(*DAGCon, in_node)).length
    &&
    Output_Edges(G_Node(*DAGCon, in_node)).list[edge->self_edge_ID].in_node == in_node
    &&
    Output_Edges(G_Node(*DAGCon, in_node)).list[edge->self_edge_ID].out_node == out_node
    )
    {
        (*e_forward) = &(Output_Edges(G_Node(*DAGCon, in_node)).list[edge->self_edge_ID]);
        (*e_backward) = &(Input_Edges(G_Node(*DAGCon, out_node)).list[edge->reverse_edge_ID]);
    }
    else
    {
        (*e_forward) = &(Output_Edges(G_Node(*DAGCon, in_node)).list[edge->reverse_edge_ID]);
        (*e_backward) = &(Input_Edges(G_Node(*DAGCon, out_node)).list[edge->self_edge_ID]);
    }
}


inline long long get_bi_Edge(Graph* DAGCon, Node* inNode, Node* outNode, Edge** e_forward, Edge** e_backward)
{
    Edge* e;
    RSet iter;
    clear_RSet(&iter);

    if(If_Node_Exist(*inNode) && If_Node_Exist(*outNode))
    {
        //find in-edge of outNode
        while(getInputEdges(&iter, DAGCon, outNode, &e))
        {
            if(e->in_node == inNode->ID)
            {
                get_bi_direction_edges(DAGCon, e, e_forward, e_backward);
                return 1;
            }
        }
    }

    return 0;
}


inline long long get_Edge_Weight(Graph* DAGCon, Node* inNode, Node* outNode)
{
    Edge* e_forward;
    Edge* e_backward;
    get_bi_Edge(DAGCon, inNode, outNode, &e_forward, &e_backward);
    return e_forward->weight;
}

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

inline Node* add_Node_DAGCon(Graph* g, char base)
{
    return &(G_Node(*g, append_Node_alloc(&g->g_nodes, base)));
}

///to delete a node
///1. set the corresponding base to be 'D'
///2. remove all related edges
///2. clear all related edges
///3. g_nodes.delete_length++, please do not substract g_nodes.length
uint64_t inline delete_Node_DAGCon(Graph* g, Node* node)
{
    g->g_nodes.delete_length++;
    g->g_nodes.list[(*node).ID].base = 'D';
    g->g_nodes.list[(*node).ID].num_insertions = (uint64_t)-1;
    g->g_nodes.list[(*node).ID].weight = (uint64_t)-1;


    RSet iter;
    Edge* e;
    clear_RSet(&iter);
    while (getOutputEdges(&iter, g, node, &e))
    {
        remove_and_check_bi_direction_edge_from_edge(g, e);
    }

    clear_RSet(&iter);
    while (getInputEdges(&iter, g, node, &e))
    {
        remove_and_check_bi_direction_edge_from_edge(g, e);
    }


    
    clear_Edge_alloc(&(g->g_nodes.list[(*node).ID].insertion_edges));
    clear_Edge_alloc(&(g->g_nodes.list[(*node).ID].mismatch_edges));
    clear_Edge_alloc(&(g->g_nodes.list[(*node).ID].deletion_edges));
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
    /**
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
    **/
   long long i;
   for (i = 0; i < deletion_length; i++)
   {
       add_single_deletionEdge_weight(g, alignNodeID + i, alignNodeID + i + 1, 0);
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
    /****************************may have bugs********************************/
    ///return edgeID;
    return return_edgeID;
    /****************************may have bugs********************************/
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


inline void extract_path(Graph* backbone, int debug_node_in_backbone, int path_i, char* pre)
{
    int step = G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[path_i].length;
    int string_i, preNode, j;
    if(step != 0)
    {
        string_i = 0;
        preNode = G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[path_i].out_node;
    
        for (j = 0; j < step; j++)
        {
            pre[string_i++] = G_Node(*backbone, preNode).base;
            preNode = G_Node(*backbone, preNode).insertion_edges.list[0].out_node;
        }
    }

    pre[string_i] = '\0';
}


inline void extract_path_debug(Graph* backbone, int debug_node_in_backbone, int path_i, char* pre)
{
    int step = G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[path_i].length;
    int string_i, preNode, preEdge, j;
    if(step != 0)
    {
        string_i = 0;
        preNode = G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[path_i].out_node;
        preEdge = G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[path_i].length;

        for (j = 0; j < step; j++)
        {
            ///pre[string_i++] = G_Node(*backbone, preNode).base;
            fprintf(stderr, "j: %d (%c%d), ", j, G_Node(*backbone, preNode).base, preEdge);
            preEdge = G_Node(*backbone, preNode).insertion_edges.list[0].length;
            preNode = G_Node(*backbone, preNode).insertion_edges.list[0].out_node;
        }
    }

    fprintf(stderr, "\n");

    ///pre[string_i] = '\0';
}


inline int getEdge_DEBUG(Graph* g, Edge_alloc* edge, uint64_t edge_length, char base)
{
    long long i = 0;
    long long nodeID;

    for (i = 0; i < edge->length; i++)
    {
        ///fprintf(stderr, "************i:%d, edge->list[i].length: %d, edge_length: %d\n",i, edge->list[i].length, edge_length);
        if (edge->list[i].length == edge_length)
        {
            nodeID = edge->list[i].out_node;
            fprintf(stderr, "########i:%d, edge->list[i].length: %d, edge_length: %d, nodeID: %d, list[nodeID].base: %c, base: %c\n",
        i, edge->list[i].length, edge_length, nodeID, g->g_nodes.list[nodeID].base, base);

            if(g->g_nodes.list[nodeID].base == base)
            {
                return i;
            }
        }
    }
    
    return -1;
}

inline int get_insertion_Edges_new(Graph* backbone, int debug_node_in_backbone, uint64_t edge_length, char* bases)
{
    int path_i, j, step, preNode;

    for (path_i = 0; path_i < G_Node(*backbone, debug_node_in_backbone).insertion_edges.length; path_i++)
    {
        step = G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[path_i].length;


        if(step != edge_length)
        {
            continue;
        }


        if(step != 0)
        {
            preNode = G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[path_i].out_node;

            for (j = 0; j < step; j++)
            {
                ///pre[string_i++] = G_Node(*backbone, preNode).base;
                ///fprintf(stderr, "path_i: %d, ID: %d\n", path_i, G_Node(*backbone, preNode).ID);
                if(G_Node(*backbone, preNode).base != bases[j])
                {
                    break;
                }

                preNode = G_Node(*backbone, preNode).insertion_edges.list[0].out_node;
            }

            if(j == step)
            {
                return path_i;
            }
        }
    }


    return -1;
}


inline int get_insertion_Edges_debug(Graph* g, Edge_alloc* edge, uint64_t edge_length, char* bases)
{
    long long i = 0;
    long long nodeID;
    long long edgeID;

    if (edge_length < 1)
    {
        return -1;
    }

    
    ///fprintf(stderr, "edge_length: %d, edge: %.*s\n", edge_length, edge_length, bases);
    

    edgeID = getEdge_DEBUG(g, edge, edge_length, bases[0]);
    fprintf(stderr, "i: %d, edgeID: %d, edge_length - i: %d\n", i, edgeID, edge_length);
    

    

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
        edgeID = getEdge_DEBUG(g, new_edge, edge_length - i, bases[i]);
        fprintf(stderr, "i: %d, edgeID: %d, edge_length - i: %d\n", i, edgeID, edge_length - i);
        if(edgeID == -1)
        {
            return -1;
        }
    }
    /****************************may have bugs********************************/
    ///return edgeID;
    return return_edgeID;
    /****************************may have bugs********************************/
}

inline void add_insertionEdge_weight(Graph* g, long long alignNodeID, char* insert, long long insert_length)
{
    
    long long nodeID;
    long long edgeID;
    Edge_alloc* edge = &(g->g_nodes.list[alignNodeID].insertion_edges);

    if (insert_length == 1)
    {
        edgeID = getEdge(g, edge, 1, insert[0]);
        // if(edgeID != get_insertion_Edges_new(g, alignNodeID, insert_length, insert))
        // {
        //     fprintf(stderr, "error\n");
        // }
        if (edgeID != -1)
        {
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
    else
    {
        ///edgeID = get_insertion_Edges(g, edge, insert_length, insert);
        edgeID = get_insertion_Edges_new(g, alignNodeID, insert_length, insert);
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

    // /******************************for homopolymer*************************/
    // long long i = 0;
    // char hom;
    // if (insert_length > 0)
    // {
    //     hom = insert[0];
    // }

    // for (i = 0; i < insert_length; i++)
    // {
    //     if(insert[i] != hom)
    //     {
    //         break;
    //     }
    // }

    // ///if it is a homopolymer
    // if(i == insert_length)
    // {
    //     ///single base
    //     edgeID = getEdge(g, edge, 1, insert[0]);
    //     if (edgeID != -1)
    //     {
    //         ///这条路均只有一个出度
    //         edge->list[edgeID].weight++;
    //     }
    //     else ///不存在这么一条边
    //     {
    //         nodeID  = add_Node_Graph(g, insert[0]);
    //         append_Edge_alloc(edge, alignNodeID, nodeID, 1, 1);
    //         ///将新加入的节点通过insertion_edges接回backbone上
    //         ///应该连回到原节点，而不是原节点的下一个节点
    //         ///append_Edge_alloc(&(g->g_nodes.list[nodeID].insertion_edges), nodeID, alignNodeID + 1, 1, 0);
    //         append_Edge_alloc(&(g->g_nodes.list[nodeID].insertion_edges), nodeID, alignNodeID, 1, 0);
    //     }

    //     ///multiple bases
    //     for (i = 1; i < insert_length; i++)
    //     {
    //         edgeID = get_insertion_Edges(g, edge, i + 1, insert);
    //         if (edgeID != -1)
    //         {
    //             ///这条路均只有一个出度
    //             edge->list[edgeID].weight++;
    //         }
    //         else
    //         {
    //             create_insertion_Edges(g, alignNodeID, i + 1, insert);
    //         }
    //     }

    //     return;
    // }
    // /******************************for homopolymer*************************/

    

    // if (insert_length == 1)
    // {
    //     edgeID = getEdge(g, edge, 1, insert[0]);
    //     if (edgeID != -1)
    //     {
    //         ///这条路均只有一个出度
    //         edge->list[edgeID].weight++;
    //     }
    //     else ///不存在这么一条边
    //     {
    //         nodeID  = add_Node_Graph(g, insert[0]);
    //         append_Edge_alloc(edge, alignNodeID, nodeID, 1, 1);
    //         ///将新加入的节点通过insertion_edges接回backbone上
    //         ///应该连回到原节点，而不是原节点的下一个节点
    //         ///append_Edge_alloc(&(g->g_nodes.list[nodeID].insertion_edges), nodeID, alignNodeID + 1, 1, 0);
    //         append_Edge_alloc(&(g->g_nodes.list[nodeID].insertion_edges), nodeID, alignNodeID, 1, 0);
    //     }
    // }
    // else if (insert_length == 2)
    // {
    //     /*******************第0个字符********************* */
    //     edgeID = getEdge(g, edge, 1, insert[0]);
    //     if (edgeID != -1)
    //     {
    //         ///这条路均只有一个出度
    //         edge->list[edgeID].weight++;
    //     }
    //     else ///不存在这么一条边
    //     {
    //         nodeID  = add_Node_Graph(g, insert[0]);
    //         append_Edge_alloc(edge, alignNodeID, nodeID, 1, 1);
    //         ///将新加入的节点通过insertion_edges接回backbone上
    //         ///应该连回到原节点，而不是原节点的下一个节点
    //         ///append_Edge_alloc(&(g->g_nodes.list[nodeID].insertion_edges), nodeID, alignNodeID + 1, 1, 0);
    //         append_Edge_alloc(&(g->g_nodes.list[nodeID].insertion_edges), nodeID, alignNodeID, 1, 0);
    //     }
    //     /*******************第0个字符********************* */

    //     /*******************第1个字符********************* */
    //     if (insert[1] != insert[0])
    //     {
    //         edgeID = getEdge(g, edge, 1, insert[1]);
    //         if (edgeID != -1)
    //         {
    //             ///这条路均只有一个出度
    //             edge->list[edgeID].weight++;
    //         }
    //         else ///不存在这么一条边
    //         {
    //             nodeID  = add_Node_Graph(g, insert[1]);
    //             append_Edge_alloc(edge, alignNodeID, nodeID, 1, 1);
    //             ///将新加入的节点通过insertion_edges接回backbone上
    //             ///应该连回到原节点，而不是原节点的下一个节点
    //             ///append_Edge_alloc(&(g->g_nodes.list[nodeID].insertion_edges), nodeID, alignNodeID + 1, 1, 0);
    //             append_Edge_alloc(&(g->g_nodes.list[nodeID].insertion_edges), nodeID, alignNodeID, 1, 0);
    //         }
    //     }
    //     /*******************第1个字符********************* */

    //     /**********************两个字符******************* */
        
    //     edgeID = get_insertion_Edges(g, edge, 2, insert);
    //     if (edgeID != -1)
    //     {
    //         ///这条路均只有一个出度
    //         edge->list[edgeID].weight++;
    //     }
    //     else
    //     {
    //         create_insertion_Edges(g, alignNodeID, insert_length, insert);
    //     }
        
    //     /**********************两个字符******************* */
    // }
    // else if (insert_length > 2)
    // {
    //     ////fprintf(stderr, "too long insertion\n");
    //     /*************************大于2个字符************************** */

    //     edgeID = get_insertion_Edges(g, edge, insert_length, insert);
    //     if (edgeID != -1)
    //     {
    //         ///这条路均只有一个出度
    //         edge->list[edgeID].weight++;
    //     }
    //     else
    //     {
    //         create_insertion_Edges(g, alignNodeID, insert_length, insert);
    //     }
    // }
    
    
    
}


void addmatchedSeqToGraph_print(Graph* backbone, long long currentNodeID, char* x_string, long long x_length, 
        char* y_string, long long y_length, CIGAR* cigar, long long backbone_start, long long backbone_end);



#endif