#ifndef __COMMAND_LINE_PARSER__
#define __COMMAND_LINE_PARSER__

#include <pthread.h>

#define HA_VERSION "0.3.0-dirty-r198"

#define VERBOSE 0

typedef struct {
    int num_reads;
    char** read_file_names;
    char* output_file_name;
    char* required_read_name;
    char* pat_index;
    char* mat_index;
    int thread_num;
    int k_mer_length;
	int mz_win;
	int bf_shift;
	float high_factor;
	int no_HPC;
	int no_kmer_flt;
	int max_n_chain;
    int k_mer_min_freq;
    int k_mer_max_freq;
    int load_index_from_disk;
    int write_index_to_disk;
    int number_of_round;
    int adapterLen;
    int clean_round;
    int roundID;
    int max_hang_Len;
    int gap_fuzz;
    int min_overlap_Len;
    int min_overlap_coverage;
    int max_short_tip;
    int min_cnt;
    int mid_cnt;
	int verbose_gfa;

    float max_hang_rate;
    float min_drop_rate;
    float max_drop_rate;

    long long small_pop_bubble_size;
    long long large_pop_bubble_size;
    long long num_bases;
    long long num_corrected_bases;
    long long num_recorrected_bases;
	long long mem_buf;
    long long coverage;
} hifiasm_opt_t;

extern hifiasm_opt_t asm_opt;

void init_opt(hifiasm_opt_t* asm_opt);
void destory_opt(hifiasm_opt_t* asm_opt);
void clear_opt(hifiasm_opt_t* asm_opt, int round);
int CommandLine_process(int argc, char *argv[], hifiasm_opt_t* asm_opt);
double Get_T(void);

#endif
