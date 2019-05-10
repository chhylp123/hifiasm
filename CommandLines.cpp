#include "CommandLines.h"
#include <stdlib.h>
#include <stdio.h>
#include "ketopt.h"

char* read_file_name = NULL;
char* output_file_name = NULL;
int thread_num = 1;

void Print_H()
{
    fprintf(stderr, "Incorrect options.\n");
}




int CommandLine_process (int argc, char *argv[])
{

    static ko_longopt_t longopts[] = {
        { "help", ko_no_argument,       100},
        { "seq", ko_required_argument, 101},
        { "output", ko_required_argument, 102},
        { "thread", ko_required_argument, 103},
        { NULL, 0, 0 }
    };

    ketopt_t opt = KETOPT_INIT;

    int i, c;
    while ((c = ketopt(&opt, argc, argv, 1, "ht:o:q:", longopts)) >= 0) {
        if (c == 100 || c == 'h') Print_H();
        else if (c == 103 || c == 't') thread_num = atoi(opt.arg); 
        else if (c == 102 || c == 'o') output_file_name = opt.arg;
        else if (c == 101 || c == 'q') read_file_name = opt.arg;
        else if (c == '?') printf("unknown opt: -%c\n", opt.opt? opt.opt : ':');
        else if (c == ':') printf("missing arg: -%c\n", opt.opt? opt.opt : ':');
    }


    if (argc == 1)
    {
        Print_H();
        return 0;
    }


    return 1;
}