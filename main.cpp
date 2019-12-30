#include <stdio.h>
#include <stdlib.h>
#include "CommandLines.h"
#include "Process_Read.h"
#include "Assembly.h"
#include "Levenshtein_distance.h"

int main(int argc, char *argv[])
{
    init_opt(&asm_opt);

    if (!CommandLine_process(argc, argv, &asm_opt)) return 1;

    Correct_Reads(asm_opt.number_of_round);

    destory_opt(&asm_opt);
    
    return 0;
}
