#include <stdio.h>
#include <stdlib.h>
#include "CommandLines.h"
#include "Process_Read.h"
#include "Assembly.h"

/********************************for debug***************************************/
///使用这个函数的时候，必须把Counting_multiple_thr()里的destory_Total_Count_Table(&TCB)注释掉
void debug_Counting()
{
    init_kseq(read_file_name);
    Verify_Counting();
    fprintf(stderr, "debug over!\n");
    destory_kseq();
}


int main(int argc, char *argv[])
{

    if (!CommandLine_process(argc, argv))
        return 1;


    fprintf(stdout, "k-mer length: %d\n",k_mer_length);

    if (load_index_from_disk && load_pre_cauculated_index())
    {
        ;
    }
    else
    {
        Counting_multiple_thr();

        Build_hash_table_multiple_thr();
    }
    

    

    ///verify_Position_hash_table();




    return 1;
}
