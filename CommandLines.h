#ifndef __COMMAND_LINE_PARSER__
#define __COMMAND_LINE_PARSER__

#include <pthread.h>


extern char* read_file_name;
extern char* output_file_name;
extern int thread_num;
extern int k_mer_length;
extern int k_mer_min_freq;
extern int k_mer_max_freq;



int CommandLine_process (int argc, char *argv[]);
double Get_T(void);

#endif