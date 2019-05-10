#include <stdio.h>
#include <stdlib.h>
#include "CommandLines.h"
#include "Process_Read.h"
#include "Assembly.h"


int main(int argc, char *argv[])
{

    if (!CommandLine_process(argc, argv))
        return 1;

    init_kseq(read_file_name);



    Counting();



    destory_kseq();


    return 1;
}
