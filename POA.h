#ifndef __POA_PARSER__
#define __POA_PARSER__
#include <stdint.h>
#include "Hash_Table.h"
#include "Process_Read.h"

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
    ///0 is match，1 is mismatch，2 means y has more bases, 3 means x has more bases
    uint64_t weight;
    uint64_t num_insertions;
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
    ///number of deletion end with current node
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

    return 1;
}

inline int getInputNodes(RSet* set, Graph* graph, Node* node, Node** get_Node)
{
    if(set->index >= (long long)Input_Edges(*node).length)
    {
        return 0;
    }

    ///skip all deleted edges
    while (
        set->index < (long long)Input_Edges(*node).length
        &&
        !(If_Edge_Exist(Input_Edges(*node).list[set->index]))
    )
    {
        set->index++;
    }


    if(
        set->index < (long long)Input_Edges(*node).length 
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
    if(set->index >= (long long)Input_Edges(*node).length)
    {
        return 0;
    }

    ///skip all deleted edges
    while (
        set->index < (long long)Input_Edges(*node).length
        &&
        !(If_Edge_Exist(Input_Edges(*node).list[set->index]))
    )
    {
        set->index++;
    }


    if(
        set->index < (long long)Input_Edges(*node).length 
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
    if(set->index >= (long long)Output_Edges(*node).length)
    {
        return 0;
    }

    ///skip all deleted edges
    while (
    set->index < (long long)Output_Edges(*node).length 
    && 
    !(If_Edge_Exist(Output_Edges(*node).list[set->index]))
    )
    {
        set->index++;
    }

    if(set->index < (long long)Output_Edges(*node).length && 
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
    if(set->index >= (long long)Output_Edges(*node).length)
    {
        return 0;
    }

    ///skip all deleted edges
    while (
    set->index < (long long)Output_Edges(*node).length 
    && 
    !(If_Edge_Exist(Output_Edges(*node).list[set->index]))
    )
    {
        set->index++;
    }

    if(set->index < (long long)Output_Edges(*node).length && 
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
    (long long)Output_Edges(G_Node(*DAGCon, in_node)).list[edge->self_edge_ID].in_node == in_node
    &&
    (long long)Output_Edges(G_Node(*DAGCon, in_node)).list[edge->self_edge_ID].out_node == out_node
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
    Edge* e_forward = NULL;
    Edge* e_backward = NULL;
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

    return 1;
}





///just for mimatch edges
inline void add_mismatchEdge_weight(Graph* g, uint64_t in_node, char base, int last_operation)
{
    long long i = 0;
    long long nodeID;
    Edge_alloc* edge = &(g->g_nodes.list[in_node].mismatch_edges);

    for (i = 0; i < (long long)edge->length; i++)
    {
        nodeID = edge->list[i].out_node;
        if(g->g_nodes.list[nodeID].base == base)
        {
            edge->list[i].weight++;
            ///if last operation is insertion
            if (last_operation == 2)
            {
                edge->list[i].num_insertions++;
            }
            
            break;
        }
    }

    ///there are no such edge
    if (i == (long long)edge->length)
    {
        nodeID  = add_Node_Graph(g, base);
        
        ///the length of match edge is 0, while the length of mismatch edge is 1
        append_Edge_alloc(edge, in_node, nodeID, 1, 1);
        ///if last operation is insertion
        if (last_operation == 2)
        {
            edge->list[edge->length - 1].num_insertions++;
        }

        ///add the mismatch_edges of new node to the backbone
        append_Edge_alloc(&(g->g_nodes.list[nodeID].mismatch_edges), nodeID, in_node + 1, 1, 0);
    }
}



inline void add_single_deletionEdge_weight(Graph* g, long long alignNodeID, long long nextNodeID, uint64_t edge_length)
{
    long long i = 0;
    long long nodeID;
    Edge_alloc* edge = &(g->g_nodes.list[alignNodeID].deletion_edges);

    for (i = 0; i < (long long)edge->length; i++)
    {
        nodeID = edge->list[i].out_node;
        if(nodeID == nextNodeID)
        {
            edge->list[i].weight++;
            break;
        }
    }

    ///there are no such edge
    if (i == (long long)edge->length)
    {
        append_Edge_alloc(edge, alignNodeID, nextNodeID, 1, edge_length);
    }
}

inline void add_deletionEdge_weight(Graph* g, long long alignNodeID, long long deletion_length)
{
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

    for (i = 0; i < (long long)edge->length; i++)
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

    for (i = 1; i < (long long)edge_length; i++)
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
    return return_edgeID;
    /****************************may have bugs********************************/
}



inline int create_insertion_Edges(Graph* g, long long alignNodeID, uint64_t edge_length, char* bases)
{
    long long i = 0;
    long long nodeID;
    ///should link back to the intial node
    ///long long backboneID = alignNodeID + 1;
    long long backboneID = alignNodeID;


    if (edge_length < 1)
    {
        return -1;
    }


    nodeID  = add_Node_Graph(g, bases[0]);
    ///add the new node to alignNodeID by insertion_edges
    append_Edge_alloc(&(g->g_nodes.list[alignNodeID].insertion_edges), alignNodeID, nodeID, 1, edge_length);

    alignNodeID = nodeID;

    for (i = 1; i < (long long)edge_length; i++)
    {
        nodeID  = add_Node_Graph(g, bases[i]);
        ///add the new node to alignNodeID by insertion_edges
        append_Edge_alloc(&(g->g_nodes.list[alignNodeID].insertion_edges), alignNodeID, nodeID, 1, edge_length - i);
        alignNodeID = nodeID;        
    }

    append_Edge_alloc(&(g->g_nodes.list[alignNodeID].insertion_edges), alignNodeID, backboneID, 1, 0);

    return 1;
}


inline void extract_path(Graph* backbone, int debug_node_in_backbone, int path_i, char* pre)
{
    int step = G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[path_i].length;
    int string_i = 0, preNode = 0, j = 0;
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



inline int get_insertion_Edges_new(Graph* backbone, int debug_node_in_backbone, uint64_t edge_length, char* bases)
{
    int path_i, j, step, preNode;

    for (path_i = 0; path_i < (long long)G_Node(*backbone, debug_node_in_backbone).insertion_edges.length; path_i++)
    {
        step = G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[path_i].length;


        if(step != (long long)edge_length)
        {
            continue;
        }


        if(step != 0)
        {
            preNode = G_Node(*backbone, debug_node_in_backbone).insertion_edges.list[path_i].out_node;

            for (j = 0; j < step; j++)
            {
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




inline void add_insertionEdge_weight(Graph* g, long long alignNodeID, char* insert, long long insert_length)
{
    
    long long nodeID;
    long long edgeID;
    Edge_alloc* edge = &(g->g_nodes.list[alignNodeID].insertion_edges);

    if (insert_length == 1)
    {
        edgeID = getEdge(g, edge, 1, insert[0]);
        if (edgeID != -1)
        {
            edge->list[edgeID].weight++;
        }
        else ///there is no such edge
        {
            nodeID  = add_Node_Graph(g, insert[0]);
            append_Edge_alloc(edge, alignNodeID, nodeID, 1, 1);
            ///add the new node to alignNodeID by insertion_edges
            //should link to the initial node, instead of the next node of the initial node
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
            ///just one outdegree
            edge->list[edgeID].weight++;
        }
        else
        {
            create_insertion_Edges(g, alignNodeID, insert_length, insert);
        }
    }
}



#endif