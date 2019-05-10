#include <stdio.h>
#include <stdlib.h>
#include "CommandLines.h"
#include "Process_Read.h"
#include "Assembly.h"


int main(int argc, char *argv[])
{

    if (!CommandLine_process(argc, argv))
        return 1;


    if (!initiReadAllReads(read_file_name))
    {
        fprintf(stdout, "Cannot open read files. \n");
        return 1;
    }


    Counting();

    


    return 1;
}
