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
/********************************for debug***************************************/


int main(int argc, char *argv[])
{

    if (!CommandLine_process(argc, argv))
        return 1;

    init_kseq(read_file_name);

    fprintf(stdout, "k-mer length: %d\n",k_mer_length);

    ///Counting();
    Counting_multiple_thr();

    destory_kseq();

    ///debug_Counting();


    return 1;
}
