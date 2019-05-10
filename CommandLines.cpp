#include "CommandLines.h"
#include <stdlib.h>
#include <getopt.h>
#include <stdio.h>

char* read_file_name = NULL;
char* output_file_name = NULL;
int thread_num = 1;

void Print_H()
{
    fprintf(stderr, "Incorrect options.\n");
}




int CommandLine_process (int argc, char *argv[])
{
    int o; 
    int index;

    static struct option longOptions[] =
    {
        {"seq",		required_argument, 0, 'q'},
        {"threads",	required_argument, 0, 't'},
        {"threads",	required_argument, 0, 'o'},
        {0,  0,  0, 0},
    };


    while ( (o = getopt_long ( argc, argv, "q:t:o:", longOptions, &index)) != -1 )
    {
        switch (o)
        {
            case 'q':
                read_file_name = optarg;
                break;
            case 'o':
                output_file_name = optarg;
                break;
            case 't':
                thread_num = atoi(optarg);
                break;
            case '?':
                fprintf(stderr, "Unrecognized option!\n" );
                return 0;   
                break;
            default:
                Print_H();
                return 0;   
        }
    }

    if (argc == 1)
    {
        Print_H();
        return 0;
    }


    return 1;
}