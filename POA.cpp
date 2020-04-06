#include <stdlib.h>
#include <string.h>
#include "POA.h"
#include "Correct.h"
#include "Process_Read.h"
#define INIT_EDGE_SIZE 50
#define INCREASE_EDGE_SIZE 5
#define INIT_NODE_SIZE 16000

/********
 * Edge *
 ********/

void init_Edge_alloc(Edge_alloc* list)
{
	if (list->list == NULL) {
		list->size = INIT_EDGE_SIZE;
		list->length = 0;
		list->delete_length = 0;
		list->list = (Edge*)malloc(sizeof(Edge)*list->size);
	} else {
		list->length = 0;
		list->delete_length = 0;
	}
}

void clear_Edge_alloc(Edge_alloc* list)
{
    list->length = 0;
    list->delete_length = 0;
}

void destory_Edge_alloc(Edge_alloc* list)
{
	if (list && list->list)
	    free(list->list);
}

void append_Edge_alloc(Edge_alloc* list,  uint64_t in_node, uint64_t out_node, uint64_t weight, uint64_t length)
{
	if (list->length + 1 > list->size) {
		uint64_t old_size = list->size;
		list->size = list->size + INCREASE_EDGE_SIZE;
		list->list = (Edge*)realloc(list->list, sizeof(Edge)*list->size);
		memset(&list->list[old_size], 0, (list->size - old_size) * sizeof(Edge));
	}

	list->list[list->length].in_node = in_node;
	list->list[list->length].out_node = out_node;
	list->list[list->length].weight = weight;
	list->list[list->length].length = length;
	list->list[list->length].num_insertions = 0;
	list->list[list->length].self_edge_ID = list->length;

	list->length++;
}

int add_and_check_bi_direction_edge(Graph* graph, Node* in_node, Node* out_node, uint64_t weight, uint64_t flag)
{
    Edge* e_forward;
    Edge* e_backward;

    //if there are no edge from in_node to out_node
    if(!get_bi_Edge(graph, in_node, out_node, &e_forward, &e_backward))
    {
        append_Edge_alloc(&(Output_Edges((*in_node))), (*in_node).ID, (*out_node).ID, weight, flag);
        append_Edge_alloc(&(Input_Edges((*out_node))), (*in_node).ID, (*out_node).ID, weight, flag);

        Output_Edges((*in_node)).list[Output_Edges((*in_node)).length - 1].reverse_edge_ID 
        = Input_Edges((*out_node)).length - 1;

        Input_Edges((*out_node)).list[Input_Edges((*out_node)).length - 1].reverse_edge_ID
        = Output_Edges((*in_node)).length - 1;

        return 1;
    }
    else//if there is an edge from in_node to out_node, do nothing
    {
        return 0;
    }
}

void add_bi_direction_edge(Graph* graph, Node* in_node, Node* out_node, uint64_t weight, uint64_t flag)
{

    append_Edge_alloc(&(Output_Edges((*in_node))), (*in_node).ID, (*out_node).ID, weight, flag);
    append_Edge_alloc(&(Input_Edges((*out_node))), (*in_node).ID, (*out_node).ID, weight, flag);

    Output_Edges((*in_node)).list[Output_Edges((*in_node)).length - 1].reverse_edge_ID 
    = Input_Edges((*out_node)).length - 1;

    Input_Edges((*out_node)).list[Input_Edges((*out_node)).length - 1].reverse_edge_ID
    = Output_Edges((*in_node)).length - 1;
}

int remove_and_check_bi_direction_edge_from_nodes(Graph* graph, Node* in_node, Node* out_node)
{
    Edge* e_forward;
    Edge* e_backward;

    //if there are no edge from in_node to out_node
    //1. remove these two edges
    //2. increase the edge_list.delete_length in both in_node and out_node
    if(get_bi_Edge(graph, in_node, out_node, &e_forward, &e_backward))
    {
        e_forward->in_node = (uint64_t)-1;
        e_forward->out_node = (uint64_t)-1;
        e_forward->weight = (uint64_t)-1;
        e_forward->length = (uint64_t)-1;
        e_forward->num_insertions = (uint64_t)-1;
        e_forward->self_edge_ID = (uint64_t)-1;
        e_forward->reverse_edge_ID = (uint64_t)-1;


        e_backward->in_node = (uint64_t)-1;
        e_backward->out_node = (uint64_t)-1;
        e_backward->weight = (uint64_t)-1;
        e_backward->length = (uint64_t)-1;
        e_backward->num_insertions = (uint64_t)-1;
        e_backward->self_edge_ID = (uint64_t)-1;
        e_backward->reverse_edge_ID = (uint64_t)-1;

        Output_Edges(*in_node).delete_length++;
        Input_Edges((*out_node)).delete_length++;

        return 1;
    }
    else//if there is an edge from in_node to out_node, do nothing
    {
        return 0;
    }
}

int remove_and_check_bi_direction_edge_from_edge(Graph* graph, Edge* e)
{
    Edge* e_forward;
    Edge* e_backward;

    if(If_Edge_Exist(*e))
    {
        get_bi_direction_edges(graph, e, &e_forward, &e_backward);
        Output_Edges(G_Node(*graph, e_forward->in_node)).delete_length++;
        Input_Edges(G_Node(*graph, e_forward->out_node)).delete_length++;

        e_forward->in_node = (uint64_t)-1;
        e_forward->out_node = (uint64_t)-1;
        e_forward->weight = (uint64_t)-1;
        e_forward->length = (uint64_t)-1;
        e_forward->num_insertions = (uint64_t)-1;
        e_forward->self_edge_ID = (uint64_t)-1;
        e_forward->reverse_edge_ID = (uint64_t)-1;


        e_backward->in_node = (uint64_t)-1;
        e_backward->out_node = (uint64_t)-1;
        e_backward->weight = (uint64_t)-1;
        e_backward->length = (uint64_t)-1;
        e_backward->num_insertions = (uint64_t)-1;
        e_backward->self_edge_ID = (uint64_t)-1;
        e_backward->reverse_edge_ID = (uint64_t)-1;

        return 1;
    }
    else
    {
        return 0;
    }
}

/********
 * Node *
 ********/

void init_Node_alloc(Node_alloc* list)
{
	memset(list, 0, sizeof(Node_alloc));
	list->size = INIT_NODE_SIZE;
	list->list = (Node*)calloc(list->size, sizeof(Node));
}

void destory_Node_alloc(Node_alloc* list)
{
	uint64_t i;
	for (i = 0; i < list->size; i++) {
		destory_Edge_alloc(&list->list[i].deletion_edges);
		destory_Edge_alloc(&list->list[i].insertion_edges);
		destory_Edge_alloc(&list->list[i].mismatch_edges);
	}
	free(list->list);
	free(list->sort.list);
	free(list->sort.visit);
	free(list->sort.iterative_buffer);
	free(list->sort.iterative_buffer_visit);
}

void clear_Node_alloc(Node_alloc* list)
{
	uint64_t i =0;
	for (i = 0; i < list->length; i++) { // TODO: is this list->size or list->length? The original version is list->length.
		clear_Edge_alloc(&list->list[i].insertion_edges);
		clear_Edge_alloc(&list->list[i].mismatch_edges);
		clear_Edge_alloc(&list->list[i].deletion_edges);
	}
	list->length = 0;
	list->delete_length = 0;
}

uint64_t append_Node_alloc(Node_alloc* list, char base)
{
	if (list->length + 1 > list->size) {
		uint64_t old_size = list->size;
		list->size = list->size * 2;
		list->list = (Node*)realloc(list->list, sizeof(Node) * list->size);
		memset(&list->list[old_size], 0, (list->size - old_size) * sizeof(Node));
	}

	list->list[list->length].ID = list->length;
	list->list[list->length].base = base;
	list->list[list->length].weight = 1;
	list->list[list->length].num_insertions = 0;
	init_Edge_alloc(&list->list[list->length].deletion_edges);
	init_Edge_alloc(&list->list[list->length].insertion_edges);
	init_Edge_alloc(&list->list[list->length].mismatch_edges);

	list->length++;

	return list->length - 1;
}

/*********
 * Graph *
 *********/

void init_Graph(Graph* g)
{
    init_Node_alloc(&g->g_nodes);
    g->g_n_edges = 0;
    g->g_n_nodes = 0;
    g->g_next_nodeID = 0;
    g->s_end_nodeID = 0;
    g->s_start_nodeID = 0;
    g->seq = NULL;
    g->seqID = (uint64_t)-1;

    init_Queue(&(g->node_q));
}

void destory_Graph(Graph* g)
{
    destory_Node_alloc(&g->g_nodes);
    destory_Queue(&(g->node_q));
}

void clear_Graph(Graph* g)
{
    clear_Node_alloc(&g->g_nodes);

    g->g_n_edges = 0;
    g->g_n_nodes = 0;
    g->g_next_nodeID = 0;
    g->s_end_nodeID = 0;
    g->s_start_nodeID = 0;
    g->seq = NULL;
    g->seqID = (uint64_t)-1;

    clear_Queue(&(g->node_q));
}

void addUnmatchedSeqToGraph(Graph* g, char* g_read_seq, long long g_read_length, long long* startID, long long* endID)
{
    long long firstID, lastID, nodeID, i;
    firstID = -1;
    lastID = -1;

    if(g_read_length == 0)
        return;

    ///start node
    nodeID  = add_Node_Graph(g, 'S');
    firstID = nodeID;
    lastID = nodeID;


    for (i = 0; i < g_read_length; i++)
    {
        nodeID = add_Node_Graph(g, g_read_seq[i]);
        
        if (firstID == -1)
        {
            firstID = nodeID;
        }
        if (lastID != -1)
        {
            ///the legnth of match edge is 0, while the length of musmatch is 1
            append_Edge_alloc(&(g->g_nodes.list[lastID].mismatch_edges), lastID, nodeID, 1, 0);
        }

        lastID = nodeID; 
    }

    *startID = firstID;
    *endID = lastID;

    g->s_start_nodeID = firstID;
    g->s_end_nodeID = lastID;
    
}

void addmatchedSeqToGraph(Graph* backbone, long long currentNodeID, char* x_string, long long x_length, 
        char* y_string, long long y_length, CIGAR* cigar, long long backbone_start, long long backbone_end)
{
    
    int x_i, y_i, cigar_i;
    x_i = 0;
    y_i = 0;
    cigar_i = 0;
    int operation;
    int operationLen;
    int i;
    int last_operation = -1;

    
    ///note that node 0 is the start node
    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2 mean y has more bases, while 3 means x has more bases
    while (cigar_i < cigar->length)
    {
        operation = cigar->C_C[cigar_i];
        operationLen = cigar->C_L[cigar_i];

        ///match/mismatch
        if (operation == 0 || operation == 1)
        {
            
            for (i = 0; i < operationLen; i++)
            {
                ///if the previous node is insertion, this node might be mismatch/match
                add_mismatchEdge_weight(backbone, currentNodeID, y_string[y_i], last_operation);    
                x_i++;
                y_i++;
                currentNodeID++;
            }
        }///insertion
        else if (operation == 2)
        {
            ///the begin and end of cigar cannot be 2, so -1 is right here
            ///if (operationLen <= CORRECT_INDEL_LENGTH)
            {
                add_insertionEdge_weight(backbone, currentNodeID, y_string + y_i, operationLen);
                backbone->g_nodes.list[currentNodeID].num_insertions++;
            }
            y_i += operationLen;
        }
        else if (operation == 3)
        {
            ///3 means x has more bases, that means backbone has more bases
            ///like a mismatch (-)
            ///if (operationLen <= CORRECT_INDEL_LENGTH)
            {
                add_deletionEdge_weight(backbone, currentNodeID, operationLen);
            }

            
            currentNodeID += operationLen;
            x_i += operationLen;
        }

        last_operation = operation;
        
        cigar_i++;
    }
}
