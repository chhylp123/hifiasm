#include "POA.h"
#include <stdlib.h>
#define INIT_EDGE_SIZE 50
#define INCREASE_EDGE_SIZE 5
#define INIT_NODE_SIZE 16000


void check_addUnmatchedSeqToGraph(Graph* g, char* g_read_seq, long long g_read_length, long long startID, long long endID)
{
    long long reverse_startID, reverse_endID;
    reverse_startID = startID;
    reverse_endID = endID;
    long long i = 0;
    if (endID - startID + 1 != g_read_length)
    {
        fprintf(stderr, "ERROR Length...\n");
        fprintf(stderr, "startID: %lld, endID: %lld\n", startID, endID);
        fprintf(stderr, "g_read_length: %lld\n", g_read_length);
    }

    if (g_read_length == 0)
    {
        return;
    }
    
    
    while (1)
    {
        if(g->g_nodes.list[startID].base != g_read_seq[i])
        {
            fprintf(stderr, "i: %llu, ERROR Node Base...\n", i);

        }

        if(g->g_nodes.list[startID].outcome_edges.length == 0)
        {
            break;
        }

        startID = g->g_nodes.list[startID].outcome_edges.list[0].out_node;

        i++;
    }

    if (startID != endID)
    {
        fprintf(stderr, "ERROR End Node Base\n");
    }


    i = g_read_length - 1;
    while (1)
    {
        if(g->g_nodes.list[reverse_endID].base != g_read_seq[i])
        {
            fprintf(stderr, "i: %llu, g_read_length: %llu, ERROR Node Base...\n", i, g_read_length);

        }

        if(g->g_nodes.list[reverse_endID].income_edges.length == 0)
        {
            break;
        }

        reverse_endID = g->g_nodes.list[reverse_endID].income_edges.list[0].in_node;

        i--;
    }

    if (reverse_startID != reverse_endID)
    {
        fprintf(stderr, "ERROR Start Node Base\n");
    }

}


void init_Edge_alloc(Edge_alloc* list)
{
    if (list->list == NULL)
    {
        list->size = INIT_EDGE_SIZE;
        list->length = 0;
        list->list = (Edge*)malloc(sizeof(Edge)*list->size);
    }
    else
    {
        list->length = 0;
    }
    
}

void clear_Edge_alloc(Edge_alloc* list)
{
    list->length = 0;
}

void destory_Edge_alloc(Edge_alloc* list)
{
    free(list->list);
}

void append_Edge_alloc(Edge_alloc* list,  uint64_t in_node, uint64_t out_node, uint64_t weight)
{
    if (list->length + 1 > list->size)
    {
        list->size = list->size + INCREASE_EDGE_SIZE;
        list->list = (Edge*)realloc(list->list, sizeof(Edge)*list->size);
    }

    list->list[list->length].in_node = in_node;
    list->list[list->length].out_node = out_node;
    list->list[list->length].weight = weight;

    list->length++;
}








void init_Node_alloc(Node_alloc* list)
{
    list->size = INIT_NODE_SIZE;
    list->length = 0;
    list->list = (Node*)malloc(sizeof(Node)*list->size);
    list->sort.size = 0;
    list->sort.list = NULL;
    list->sort.visit = NULL;

    list->sort.iterative_buffer = NULL;
    list->sort.iterative_buffer_visit = NULL;


    long long i;
    for (i = 0; i < list->size; i++)
    {
        list->list[i].income_edges.list=NULL;
        list->list[i].outcome_edges.list=NULL;
        list->list[i].alignedTo_Nodes.list=NULL;
    }
    
}

void destory_Node_alloc(Node_alloc* list)
{
    uint64_t i =0;
    for (i = 0; i < list->length; i++)
    {
        destory_Edge_alloc(&list->list[i].income_edges);
        destory_Edge_alloc(&list->list[i].outcome_edges);
        destory_Edge_alloc(&list->list[i].alignedTo_Nodes);
    }
    free(list->list);
    free(list->sort.list);
    free(list->sort.visit);
    free(list->sort.iterative_buffer);
    free(list->sort.iterative_buffer_visit);
    ///free(list->topo_order);
}

void clear_Node_alloc(Node_alloc* list)
{
    uint64_t i =0;
    for (i = 0; i < list->length; i++)
    {
        clear_Edge_alloc(&list->list[i].income_edges);
        clear_Edge_alloc(&list->list[i].outcome_edges);
        clear_Edge_alloc(&list->list[i].alignedTo_Nodes);
    }
    list->length = 0;
}


uint64_t append_Node_alloc(Node_alloc* list, char base)
{
   
    if (list->length + 1 > list->size)
    {
        long long i = list->size;

        ///list->topo_order这里用不到，所以不用先分配空间
        ///但是还是一起分配了吧，免得麻烦
        list->size = list->size * 2;
        list->list = (Node*)realloc(list->list, sizeof(Node)*list->size);
        ///list->topo_order = (uint64_t*)realloc(list->topo_order, sizeof(uint64_t)*list->size);
        
        for (; i < list->size; i++)
        {
            list->list[i].income_edges.list=NULL;
            list->list[i].outcome_edges.list=NULL;
            list->list[i].alignedTo_Nodes.list=NULL;
        }
    }
    
    list->list[list->length].ID = list->length;
    list->list[list->length].base = base;
    init_Edge_alloc(&list->list[list->length].income_edges);
    init_Edge_alloc(&list->list[list->length].outcome_edges);
    init_Edge_alloc(&list->list[list->length].alignedTo_Nodes);
    
    list->length++;

    return list->length - 1;
}











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
}

void destory_Graph(Graph* g)
{
    destory_Node_alloc(&g->g_nodes);
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
}





uint64_t inline add_Node_Graph(Graph* g, char base)
{
    return append_Node_alloc(&g->g_nodes, base);
}


void inline add_Edge_Graph(Graph* g, uint64_t start, uint64_t end, uint64_t weight)
{
    if (start >= g->g_nodes.length || end >= g->g_nodes.length)
    {
        fprintf(stderr, "Not existing nodes ...");
        exit(0);
    }
    
    ///对起始节点加出边
    append_Edge_alloc(&g->g_nodes.list[start].outcome_edges, start, end, weight);
    ///对结束节点加入边
    append_Edge_alloc(&g->g_nodes.list[end].income_edges, start, end, weight);
}

void addUnmatchedSeqToGraph(Graph* g, char* g_read_seq, long long g_read_length, long long* startID, long long* endID)
{
    long long firstID, lastID, nodeID, i;
    firstID = -1;
    lastID = -1;

    if(g_read_length == 0)
        return;
    
    for (i = 0; i < g_read_length; i++)
    {
        nodeID = add_Node_Graph(g, g_read_seq[i]);

        ////fprintf(stderr, "nodeID: %llu\n", nodeID);
        
        if (firstID == -1)
        {
            firstID = nodeID;
        }
        if (lastID != -1)
        {
            add_Edge_Graph(g, lastID, nodeID, 1);
        }

        lastID = nodeID; 
    }

    *startID = firstID;
    *endID = lastID;
    
}


void Perform_POA(Graph* g, overlap_region_alloc* overlap_list, All_reads* R_INF, UC_Read* g_read)
{
    long long startNodeID, endNodeID;
    ///第一条序列是read本身，g_read里存的是反向互补，所以首先要恢复回正向
    ///这一步可以优化掉
    reverse_complement(g_read->seq, g_read->length);
    addUnmatchedSeqToGraph(g, g_read->seq, g_read->length, &startNodeID, &endNodeID);

    /**
    if (startNodeID!=0||endNodeID!=g_read->length-1)
    {
        fprintf(stderr, "Error startNodeID or endNodeID ...\n");
    }
    check_addUnmatchedSeqToGraph(g, g_read->seq, g_read->length, startNodeID, endNodeID);
    **/

   get_Topo_Sort_Order(&g->g_nodes, 0);

   uint64_t i = 0;
   for ( i = 0; i < overlap_list->length; i++)
   {
       /**
       if(overlap_list->list[i].x_id == overlap_list->list[i].y_id)
       {
           fprintf(stderr, "Error x_id or y_id ...\n");
       }
       **/

       if (overlap_list->list[i].y_pos_strand)
       {
           recover_UC_Read_RC(g_read, R_INF, overlap_list->list[i].y_id);
       }
       else
       {
           recover_UC_Read(g_read, R_INF, overlap_list->list[i].y_id);
       }

   }

}

void topologicalSortDFS(Node_alloc* list, uint64_t nodeID)
{
    list->sort.visit[nodeID] = 1;

    long long i;
    uint64_t out_nodeID;

    for (i = 0; i < list->list[nodeID].outcome_edges.length; i++)
    {
        out_nodeID = list->list[nodeID].outcome_edges.list[i].out_node;
        if (list->sort.visit[out_nodeID] == 0)
        {
            topologicalSortDFS(list, out_nodeID);
        }
    }

    list->sort.length--;
    list->sort.list[list->sort.length] = nodeID;
}

#define INIT_STACK(stack) stack.iterative_i = 0;
#define PUSH(stack, nodeID, time) stack.iterative_buffer[stack.iterative_i]=nodeID;\
stack.iterative_buffer_visit[stack.iterative_i++]=time;
#define IF_EMPTY(stack) (stack.iterative_i == 0)
#define POP(stack, nodeID, time) --stack.iterative_i;nodeID = stack.iterative_buffer[stack.iterative_i];\
time = stack.iterative_buffer_visit[stack.iterative_i];


void topologicalSortDFS_Iterative(Node_alloc* list, uint64_t nodeID)
{
    long long i;
    uint64_t out_nodeID;
    int flag;

    INIT_STACK(list->sort);
    PUSH(list->sort, nodeID, 0);

    while (!IF_EMPTY(list->sort))
    {
        POP(list->sort, nodeID, flag);
        ///flag == 1说明是第二次访问; flag == 0说明是第一次访问
        if (flag)
        {
            list->sort.length--;
            list->sort.list[list->sort.length] = nodeID;
            continue;
        }
        
        list->sort.visit[nodeID] = 1;
        PUSH(list->sort, nodeID, 1);

        for (i = 0; i < list->list[nodeID].outcome_edges.length; i++)
        {
            out_nodeID = list->list[nodeID].outcome_edges.list[i].out_node;
            if (list->sort.visit[out_nodeID] == 0)
            {
                PUSH(list->sort, out_nodeID, 0);
            }
        }
        
    }
}

uint64_t* get_Topo_Sort_Order(Node_alloc* list, int need_sort)
{
    long long i = 0;
    list->sort.length = list->length;
    if (list->length > list->sort.size)
    {
        list->sort.size = list->length;
        list->sort.list = (uint64_t*)realloc(list->sort.list, sizeof(uint64_t)*list->sort.size);
        list->sort.visit = (uint8_t*)realloc(list->sort.visit, sizeof(uint8_t)*list->sort.size);

        list->sort.iterative_buffer 
            = (uint64_t*)realloc(list->sort.iterative_buffer, sizeof(uint64_t)*list->sort.size);
        list->sort.iterative_buffer_visit 
            = (uint8_t*)realloc(list->sort.iterative_buffer_visit, sizeof(uint8_t)*list->sort.size);
    }
    
    if (!need_sort)
    {
        ///可以循环展开, 作用微乎其微
        for (i = 0; i < list->length; i++)
        {
            list->sort.list[i] = i;
        }
               
    }
    else
    {
        memset(list->sort.visit, 0 , list->length);
        for (i = 0; i < list->length; i++)
        {
            if(list->sort.visit[i] == 0)
            {
                ///topologicalSortDFS(list, i);
                topologicalSortDFS_Iterative(list, i);
            }
        }


        /**
        for (i = 0; i < list->length; i++)
        {
            if (list->sort.list[i] != i)
            {
                fprintf(stderr, "ERROR Sort ....\n");
            }
        }
        **/
        

    }
    
    return list->sort.list;
    
}
