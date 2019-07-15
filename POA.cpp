#include "POA.h"
#include <stdlib.h>
#include "Correct.h"
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

    list->total_start.income_edges.list = NULL;
    list->total_start.outcome_edges.list = NULL;
    list->total_start.alignedTo_Nodes.list = NULL;
    
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
    destory_Edge_alloc(&list->total_start.income_edges);
    destory_Edge_alloc(&list->total_start.outcome_edges);
    destory_Edge_alloc(&list->total_start.alignedTo_Nodes);

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
    clear_Edge_alloc(&list->total_start.income_edges);
    clear_Edge_alloc(&list->total_start.outcome_edges);
    clear_Edge_alloc(&list->total_start.alignedTo_Nodes);
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
    list->list[list->length].weight = 1;
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
            ///0是match边
            add_Edge_Graph(g, lastID, nodeID, 0);
        }

        lastID = nodeID; 
    }

    *startID = firstID;
    *endID = lastID;
    
}

inline long long get_alignToNode(Graph* backbone, long long currentNodeID, char base)
{
    if(backbone->g_nodes.list[currentNodeID].alignedTo_Nodes.length == 0)
    {
        return -1;
    }

    long long i = 0;
    long long nodeID;
    for (i = 0; i < backbone->g_nodes.list[currentNodeID].alignedTo_Nodes.length; i++)
    {
        nodeID = backbone->g_nodes.list[currentNodeID].alignedTo_Nodes.list[i].out_node;
        if(backbone->g_nodes.list[nodeID].base == base)
        {
            return nodeID;
        }
    }

    return -1;
    
}


inline void add_mismatch_to_backbone(Graph* backbone, long long* alignNodeID, char* mis_base, long long mis_base_length)
{
    long long i;
    long long mismatch_nodeID;
    char base;
    
    for (i = 0; i < mis_base_length; i++, (*alignNodeID)++)
    {
        base = mis_base[i];
        mismatch_nodeID = get_alignToNode(backbone, *alignNodeID, base);

        ///如果已经存在一个误配节点，那么给误配节点的权重+1
        if (mismatch_nodeID != -1)
        {
            backbone->g_nodes.list[mismatch_nodeID].weight++;
        }
        else
        {
            mismatch_nodeID = add_Node_Graph(backbone, base);
            ///1代表是mismatch边
            append_Edge_alloc(&backbone->g_nodes.list[*alignNodeID].alignedTo_Nodes, 
                    *alignNodeID, mismatch_nodeID, 1);
        }
    }

}

inline void add_deletion_to_backbone(Graph* backbone, long long* alignNodeID, long long deletion_length)
{
    long long i;
    long long mismatch_nodeID;
    char base;
    
    for (i = 0; i < deletion_length; i++, (*alignNodeID)++)
    {
        base = 'D';
        mismatch_nodeID = get_alignToNode(backbone, *alignNodeID, base);

        ///如果已经存在一个误配节点，那么给误配节点的权重+1
        if (mismatch_nodeID != -1)
        {
            backbone->g_nodes.list[mismatch_nodeID].weight++;
        }
        else
        {
            mismatch_nodeID = add_Node_Graph(backbone, base);
            ///3代表是deletion边
            append_Edge_alloc(&backbone->g_nodes.list[*alignNodeID].alignedTo_Nodes, 
                    *alignNodeID, mismatch_nodeID, 3);
        }
    }

}



///注意这里返回的有可能是新加的节点，也有可能返回的是backbone上的节点
inline long long get_insertion_Node(Graph* backbone, long long currentNodeID, char base)
{
    ///看出边数量
    if(backbone->g_nodes.list[currentNodeID].outcome_edges.length == 0)
    {
        return -1;
    }

    long long i = 0;
    long long nodeID;
    int type;
    ///遍历所有出边
    for (i = 0; i < backbone->g_nodes.list[currentNodeID].outcome_edges.length; i++)
    {
        ///出边类型
        type = backbone->g_nodes.list[currentNodeID].outcome_edges.list[i].weight;

        ///为2的时候才是deletion边
        ///这个边有可能是match边，也就是type = 0
        ///这个似乎不需要...，加上反而坏事
        ///也不一定
        if (type == 2)
        {
            nodeID = backbone->g_nodes.list[currentNodeID].outcome_edges.list[i].out_node;
            if(backbone->g_nodes.list[nodeID].base == base)
            {
                return nodeID;
            }
        }
    }

    return -1;
}



///注意这里返回的有可能是新加的节点，也有可能返回的是backbone上的节点
inline void link_insertion_Node(Graph* backbone, long long currentNodeID, long long backboneNodeID)
{
    ///如果出边数量为0，那这就是个新节点
    if(backbone->g_nodes.list[currentNodeID].outcome_edges.length == 0)
    {

        ///2代表是insertion边
        add_Edge_Graph(backbone, currentNodeID, backboneNodeID, 2);
    }
    else  ////如果不为0，backboneNodeID应该一定在出边中
    {
        long long i = 0;
        long long nodeID;
        int type;
        ///遍历所有出边
        for (i = 0; i < backbone->g_nodes.list[currentNodeID].outcome_edges.length; i++)
        {
            ///出边类型
            type = backbone->g_nodes.list[currentNodeID].outcome_edges.list[i].weight;

            nodeID = backbone->g_nodes.list[currentNodeID].outcome_edges.list[i].out_node;
            if(backboneNodeID == nodeID)
            {
                break;
            }
            
        }

        ///如果这个节点没有被连到backboneNodeID上，就要处理
        if (i >= backbone->g_nodes.list[currentNodeID].outcome_edges.length)
        {
            ///2代表是insertion边
            add_Edge_Graph(backbone, currentNodeID, backboneNodeID, 2);
        }
        
    }
    
}


inline void add_insertion_to_backbone(Graph* backbone, long long alignNodeID, char* insertion_base, long long insertion_length, 
long long backbone_start, long long backbone_end)
{
    long long i;
    long long insertion_nodeID;
    long long backboneID = alignNodeID + 1;
    char base;
    
    for (i = 0; i < insertion_length; i++)
    {
        base = insertion_base[i];
        insertion_nodeID = get_insertion_Node(backbone, alignNodeID, base);

        ///如果已经存在一个insertion节点，那么给insertion节点的权重+1
        ///注意这里返回的有可能是新加的节点，也有可能返回的是backbone上的节点
        ///不可能，这里要是返回了backbone上的节点就错了，最后要验证下
        if (insertion_nodeID != -1)
        {
            /**
            if (insertion_nodeID >= backbone_start && insertion_nodeID <= backbone_end)
            {
                fprintf(stderr, "error\n");
            }
            **/
            
            backbone->g_nodes.list[insertion_nodeID].weight++;
        }
        else
        {
            insertion_nodeID = add_Node_Graph(backbone, base);
            ///2代表是insertion边
            add_Edge_Graph(backbone, alignNodeID, insertion_nodeID, 2);
        }

        alignNodeID = insertion_nodeID;   
    }

    ///最后要把节点接回到backbone上去
    link_insertion_Node(backbone, alignNodeID, backboneID);
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

    
    

    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2是x缺字符（y多字符），而3是y缺字符（x多字符）
    ///while (x_i < x_len && y_i < y_len && cigar_i < cigar->length)
    while (cigar_i < cigar->length)
    {
        operation = cigar->C_C[cigar_i];
        operationLen = cigar->C_L[cigar_i];

        ///这种情况代表匹配
        if (operation == 0)
        {
            
            for (i = 0; i < operationLen; i++)
            {
                backbone->g_nodes.list[currentNodeID].weight++;

                x_i++;
                y_i++;
                currentNodeID++;
            }
        }
        else if (operation == 1)
        {
            add_mismatch_to_backbone(backbone, &currentNodeID, y_string + y_i, operationLen);
            x_i = x_i + operationLen;
            y_i = y_i + operationLen;   
        }
        else if (operation == 2)
        {
            ///记住要传currentNodeID - 1而不是currentNodeID
            ///cigar的起始和结尾不可能是2，所以这里-1没问题
            add_insertion_to_backbone(backbone, currentNodeID - 1, y_string + y_i, operationLen, backbone_start, backbone_end);
            y_i += operationLen;
        }
        else if (operation == 3)
        {
            ///3是y缺字符（x多字符），也就是backbone多字符
            ///这个相当于在backbone对应字符处变成了‘——’
            ///因此可以用mismatch类似的方法处理
            add_deletion_to_backbone(backbone, &currentNodeID, operationLen);
            x_i += operationLen;
        }
        
        cigar_i++;
    }


    /**
    ///cigar的起始和结尾不可能是2
    if (cigar->C_C[0] == 2 || cigar->C_C[cigar->length - 1] == 2)
    {
        fprintf(stderr, "error\n");
    }
    

    if (x_i != x_length)
    {
        fprintf(stderr, "x_i: %d, x_length: %d\n", x_i, x_length);
    }

    if (y_i != y_length)
    {
        fprintf(stderr, "y_i: %d, y_length: %d\n", y_i, y_length);
    }
    **/
    
}

void Graph_debug(Graph* backbone, long long currentNodeID, char* x_string, long long x_length, 
        char* y_string, long long y_length, CIGAR* cigar, long long backbone_start, long long backbone_end)
{
    int x_i, y_i, cigar_i;
    x_i = 0;
    y_i = 0;
    cigar_i = 0;
    int operation;
    int operationLen;
    int i;

    
    

    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2是x缺字符（y多字符），而3是y缺字符（x多字符）
    ///while (x_i < x_len && y_i < y_len && cigar_i < cigar->length)
    while (cigar_i < cigar->length)
    {
        operation = cigar->C_C[cigar_i];
        operationLen = cigar->C_L[cigar_i];

        ///这种情况代表匹配
        if (operation == 0)
        {
            
            for (i = 0; i < operationLen; i++)
            {
                if (backbone->g_nodes.list[currentNodeID].base != y_string[y_i])
                {
                    fprintf(stderr, "error match\n");
                }

                backbone->g_nodes.list[currentNodeID].weight--;
            
                x_i++;
                y_i++;
                currentNodeID++;
            }
        }
        else if (operation == 1)
        {
            for (i = 0; i < operationLen; i++)
            {
                if (backbone->g_nodes.list[currentNodeID].base == y_string[y_i])
                {
                    fprintf(stderr, "error mismatch 1\n");
                }

                long long mismatchID = get_alignToNode(backbone, currentNodeID, y_string[y_i]);



                if(mismatchID == -1)
                {
                    fprintf(stderr, "error mismatch 2\n");
                }
                else
                {
                    backbone->g_nodes.list[mismatchID].weight--;
                }
                

                x_i++;
                y_i++;
                currentNodeID++;
            } 
        }
        else if (operation == 2)
        {
            long long nodeID = currentNodeID - 1;
            long long mismatchID;

            for (i = 0; i < operationLen; i++)
            {
                mismatchID = get_insertion_Node(backbone, nodeID, y_string[y_i]);

                if (mismatchID == -1)
                {
                    fprintf(stderr, "error insertion 1, i: %d\n", i);
                }
                else
                {
                    backbone->g_nodes.list[mismatchID].weight--;
                }

                nodeID = mismatchID;
                
                y_i++;
            }
            ///注意这里是x_string[x_i]而不是x_string[currentNodeID]
            mismatchID = get_insertion_Node(backbone, nodeID, x_string[x_i]);
            if (mismatchID == -1)
            {
                fprintf(stderr, "error insertion 2, i: %d, x_i: %d\n", i, x_i);
            }

            
            if (mismatchID != currentNodeID)
            {
                fprintf(stderr, "error insertion 3, i: mismatchID: %d, currentNodeID: %d\n", mismatchID, currentNodeID);
            }
            
            


        }
        else if (operation == 3)
        {
            for (i = 0; i < operationLen; i++)
            {
 
                long long mismatchID = get_alignToNode(backbone, currentNodeID, 'D');

                if(mismatchID == -1)
                {
                    fprintf(stderr, "error deletion 2\n");
                }
                else
                {
                    backbone->g_nodes.list[mismatchID].weight--;
                }

                x_i++;
                currentNodeID++;
            } 
        }
        
        cigar_i++;
    }


    if (cigar->C_C[0] == 2 || cigar->C_C[cigar->length - 1] == 2)
    {
        fprintf(stderr, "error\n");
    }
    

    if (x_i != x_length)
    {
        fprintf(stderr, "x_i: %d, x_length: %d\n", x_i, x_length);
    }

    if (y_i != y_length)
    {
        fprintf(stderr, "y_i: %d, y_length: %d\n", y_i, y_length);
    }
    
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
