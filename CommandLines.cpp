#include "CommandLines.h"
#include <zlib.h>
#include <stdlib.h>
#include <stdio.h>
#include "ketopt.h"
#include <sys/time.h>

#define VERSION "0.2.0"
#define DEFAULT_OUTPUT "hifiasm.asm"

hifiasm_opt_t asm_opt;

double Get_T(void)
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec+t.tv_usec/1000000.0;
}

void Print_H(hifiasm_opt_t* asm_opt)
{
    fprintf(stderr, "Usage: hifiasm [options] <in_1.fq> <in_2.fq> <...>\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  Assembly:\n");
    fprintf(stderr, "    -o FILE       prefix of output files [%s]\n", asm_opt->output_file_name);
    ///fprintf(stderr, "    -c FILE       file including trio information\n");
    fprintf(stderr, "    -t INT        number of threads [%d]\n", asm_opt->thread_num);
    fprintf(stderr, "    -r INT        round of correction [%d]\n", asm_opt->number_of_round);
    fprintf(stderr, "    -a INT        round of assembly cleaning [%d]\n", asm_opt->clean_round);
    fprintf(stderr, "    -k INT        k-mer length [%d] (must be < 64)\n", asm_opt->k_mer_length);
    ///fprintf(stderr, "    -w            write all overlaps to disk, can accelerate assembly next time [%d]\n", asm_opt->write_index_to_disk);
    ///fprintf(stderr, "    -l            load all overlaps from disk, can avoid overlap calculation [%d]\n", asm_opt->load_index_from_disk);
    ///fprintf(stderr, "    -i            ignore saved overlaps in *.ovlp*.bin files\n");
    fprintf(stderr, "    -i            ignore saved overlaps in *.ovlp* files\n");
    fprintf(stderr, "    -z INT        length of adapters that should be removed [%d]\n", asm_opt->adapterLen);
    fprintf(stderr, "    -m INT        size of popped large bubbles for contig graph [%lld]\n", 
    asm_opt->large_pop_bubble_size);
    fprintf(stderr, "    -p INT        size of popped small bubbles for haplotype-resolved unitig graph [%lld]\n", 
    asm_opt->small_pop_bubble_size);
    fprintf(stderr, "    -n INT        small removed unitig threshold [%d]\n", asm_opt->max_short_tip);
    fprintf(stderr, "    -x FLOAT      max overlap drop ratio [%.2g]\n", asm_opt->max_drop_rate);
    fprintf(stderr, "    -y FLOAT      min overlap drop ratio [%.2g]\n", asm_opt->min_drop_rate);
    fprintf(stderr, "    -v            show version number\n");
    fprintf(stderr, "    -h            show help information\n");

    fprintf(stderr, "  Trio-partition:\n");
    fprintf(stderr, "    -P FILE       paternal trio index [NULL]\n");
    fprintf(stderr, "    -M FILE       Maternal trio index [NULL]\n");
    fprintf(stderr, "    -c INT        lower bound of the binned k-mer's frequency [%d]\n", asm_opt->min_cnt);
    fprintf(stderr, "    -d INT        upper bound of the binned k-mer's frequency [%d]\n", asm_opt->mid_cnt);

    fprintf(stderr, "Example: ./hifiasm -o NA12878.asm -t 32 NA12878.fq.gz\n");
    fprintf(stderr, "See `man ./hifiasm.1' for detailed description of these command-line options.\n");
}

void init_opt(hifiasm_opt_t* asm_opt)
{
    asm_opt->coverage = -1;
    asm_opt->num_reads = 0;
    asm_opt->read_file_names = NULL;
    asm_opt->output_file_name = (char*)(DEFAULT_OUTPUT);
    asm_opt->required_read_name = NULL;
    asm_opt->pat_index = NULL;
    asm_opt->mat_index = NULL;
    asm_opt->thread_num = 1;
    asm_opt->k_mer_length = 40;
    asm_opt->k_mer_min_freq = 3;
    asm_opt->k_mer_max_freq = 66;
    asm_opt->load_index_from_disk = 1;
    asm_opt->write_index_to_disk = 1;
    asm_opt->number_of_round = 2;
    asm_opt->adapterLen = 0;
    asm_opt->clean_round = 4;
    asm_opt->complete_threads = 0;
    asm_opt->small_pop_bubble_size = 100000;
    asm_opt->large_pop_bubble_size = 10000000;
    asm_opt->min_drop_rate = 0.2;
    asm_opt->max_drop_rate = 0.8;
    asm_opt->max_hang_Len = 1000;
    asm_opt->max_hang_rate = 0.8;
    asm_opt->gap_fuzz = 1000;
    asm_opt->min_overlap_Len = 50;
    asm_opt->min_overlap_coverage = 0;
    asm_opt->max_short_tip = 3;
    asm_opt->min_cnt = 2;
    asm_opt->mid_cnt = 5;
}

void destory_opt(hifiasm_opt_t* asm_opt)
{
    if(asm_opt->read_file_names != NULL)
    {
        free(asm_opt->read_file_names);
    }
}

void clear_opt(hifiasm_opt_t* asm_opt, int last_round)
{
    asm_opt->complete_threads = 0;
    asm_opt->num_bases = 0;
    asm_opt->num_corrected_bases = 0;
    asm_opt->num_recorrected_bases = 0;
    asm_opt->roundID = asm_opt->number_of_round - last_round;
}

int check_file(char* name, const char* opt)
{
    if(!name)
    {
        fprintf(stderr, "[ERROR] file does not exist (-%s)\n", opt);
        return 0;
    } 
    FILE* is_exist = NULL;
    is_exist = fopen(name,"r");
    if(!is_exist)
    {
        fprintf(stderr, "[ERROR] %s does not exist (-%s)\n", name, opt);
        return 0;
    } 

    fclose(is_exist);
    return 1;
}

int check_option(hifiasm_opt_t* asm_opt)
{
    if(asm_opt->read_file_names == NULL || asm_opt->num_reads == 0)
    {
        fprintf(stderr, "[ERROR] missing input: please specify a read file\n");
        return 0;
    }

    if(asm_opt->output_file_name == NULL)
    {
        fprintf(stderr, "[ERROR] missing output: please specify the output name (-o)\n");
        return 0;
    }

    if(asm_opt->thread_num < 1)
    {
        fprintf(stderr, "[ERROR] the number of threads must be > 0 (-t)\n");
        return 0;
    }


    if(asm_opt->number_of_round < 1)
    {
        fprintf(stderr, "[ERROR] the number of rounds for correction must be > 0 (-r)\n");
        return 0;
    }

    if(asm_opt->clean_round < 1)
    {
        fprintf(stderr, "[ERROR] the number of rounds for assembly cleaning must be > 0 (-a)\n");
        return 0;
    }

    if(asm_opt->adapterLen < 0)
    {
        fprintf(stderr, "[ERROR] the length of removed adapters must be >= 0 (-z)\n");
        return 0;
    }


    if(asm_opt->k_mer_length >= 64)
    {
        fprintf(stderr, "[ERROR] the length of k_mer must be < 64 (-k)\n");
        return 0;
    }


    if(asm_opt->max_drop_rate < 0 || asm_opt->max_drop_rate >= 1 ) 
    {
        fprintf(stderr, "[ERROR] max overlap drop ratio must be [0.0, 1.0) (-x)\n");
        return 0;
    }

    
    if(asm_opt->min_drop_rate < 0 || asm_opt->min_drop_rate >= 1)
    {
        fprintf(stderr, "[ERROR] min overlap drop ratio must be [0.0, 1.0) (-y)\n");
        return 0;
    }

    if(asm_opt->max_drop_rate <= asm_opt->min_drop_rate)
    {
        fprintf(stderr, "[ERROR] min overlap drop ratio must be less than max overlap drop ratio (-x/-y)\n");
        return 0;
    }

    if(asm_opt->small_pop_bubble_size < 0)
    {
        fprintf(stderr, "[ERROR] the size of popped small bubbles must be >= 0 (-p)\n");
        return 0;
    }

    if(asm_opt->large_pop_bubble_size < 0)
    {
        fprintf(stderr, "[ERROR] the size of popped large bubbles must be >= 0 (-m)\n");
        return 0;
    }

    if(asm_opt->max_hang_Len < 0)
    {
        fprintf(stderr, "[ERROR] max_hang_Len must be >= 0\n");
        return 0;
    }

    if(asm_opt->max_hang_rate < 0)
    {
        fprintf(stderr, "[ERROR] max_hang_rate must be >= 0\n");
        return 0;
    }

    if(asm_opt->gap_fuzz < 0)
    {
        fprintf(stderr, "[ERROR] gap_fuzz must be >= 0\n");
        return 0;
    }

    if(asm_opt->min_overlap_Len < 0)
    {
        fprintf(stderr, "[ERROR] min_overlap_Len must be >= 0\n");
        return 0;
    }

    if(asm_opt->min_overlap_coverage < 0)
    {
        fprintf(stderr, "[ERROR] min_overlap_coverage must be >= 0\n");
        return 0;
    }


    if(asm_opt->max_short_tip < 0)
    {
        fprintf(stderr, "[ERROR] the length of removal tips must be >= 0 (-n)\n");
        return 0;
    }


    if(asm_opt->pat_index != NULL && check_file(asm_opt->pat_index, "P") == 0) return 0;
    if(asm_opt->mat_index != NULL && check_file(asm_opt->mat_index, "M") == 0) return 0;

    // fprintf(stderr, "input file num: %d\n", asm_opt->num_reads);
    // fprintf(stderr, "output file: %s\n", asm_opt->output_file_name);
    // fprintf(stderr, "number of threads: %d\n", asm_opt->thread_num);
    // fprintf(stderr, "number of rounds for correction: %d\n", asm_opt->number_of_round);
    // fprintf(stderr, "number of rounds for assembly cleaning: %d\n", asm_opt->clean_round);
    // fprintf(stderr, "length of removed adapters: %d\n", asm_opt->adapterLen);
    // fprintf(stderr, "length of k_mer: %d\n", asm_opt->k_mer_length);
    // fprintf(stderr, "min overlap drop ratio: %.2g\n", asm_opt->min_drop_rate);
    // fprintf(stderr, "max overlap drop ratio: %.2g\n", asm_opt->max_drop_rate);
    // fprintf(stderr, "size of popped small bubbles: %lld\n", asm_opt->small_pop_bubble_size);
    // fprintf(stderr, "size of popped large bubbles: %lld\n", asm_opt->large_pop_bubble_size);
    // fprintf(stderr, "small removed unitig threshold: %d\n", asm_opt->max_short_tip);

    return 1;
}

void get_queries(int argc, char *argv[], ketopt_t* opt, hifiasm_opt_t* asm_opt)
{
    if(opt->ind == argc)
    {
        return;
    }

    asm_opt->num_reads = argc - opt->ind;
    asm_opt->read_file_names = (char**)malloc(sizeof(char*)*asm_opt->num_reads);
    
    long long i;
    gzFile dfp;
    for (i = 0; i < asm_opt->num_reads; i++)
    {
        asm_opt->read_file_names[i] = argv[i + opt->ind];
        dfp = gzopen(asm_opt->read_file_names[i], "r");
        if (dfp == 0)
        {
            fprintf(stderr, "[ERROR] Cannot find the input read file: %s\n", 
                    asm_opt->read_file_names[i]);
		    exit(0);
        }
        gzclose(dfp);
    }
}


int CommandLine_process(int argc, char *argv[], hifiasm_opt_t* asm_opt)
{
    ketopt_t opt = KETOPT_INIT;

    int c;

    while ((c = ketopt(&opt, argc, argv, 1, "hvt:o:k:lwm:n:r:a:b:z:x:y:p:c:d:M:P:i", 0)) >= 0) {
        if (c == 'h')
        {
            Print_H(asm_opt);
            return 0;
        } 
        else if (c == 'v')
        {
            fprintf(stderr, "[Version] %s\n", VERSION);
            return 0;
        } 
        else if (c == 't') asm_opt->thread_num = atoi(opt.arg); 
        else if (c == 'o') asm_opt->output_file_name = opt.arg;
        else if (c == 'r') asm_opt->number_of_round = atoi(opt.arg);
        else if (c == 'k') asm_opt->k_mer_length = atoi(opt.arg);
        else if (c == 'i') asm_opt->load_index_from_disk = 0; 
        else if (c == 'l') asm_opt->load_index_from_disk = 1; 
        else if (c == 'w') asm_opt->write_index_to_disk = 1;
        else if (c == 'a') asm_opt->clean_round = atoi(opt.arg); 
        else if (c == 'z') asm_opt->adapterLen = atoi(opt.arg);
        else if (c == 'b') asm_opt->required_read_name = opt.arg;
        else if (c == 'c') asm_opt->min_cnt = atoi(opt.arg);
        else if (c == 'd') asm_opt->mid_cnt = atoi(opt.arg);
        else if (c == 'P') asm_opt->pat_index = opt.arg;
        else if (c == 'M') asm_opt->mat_index = opt.arg;
        else if (c == 'x') asm_opt->max_drop_rate = atof(opt.arg);
        else if (c == 'y') asm_opt->min_drop_rate = atof(opt.arg);
        else if (c == 'p') asm_opt->small_pop_bubble_size = atoll(opt.arg);
        else if (c == 'm') asm_opt->large_pop_bubble_size = atoll(opt.arg);
        else if (c == 'n') asm_opt->max_short_tip = atoll(opt.arg);
        else if (c == ':') 
        {
			fprintf(stderr, "[ERROR] missing option argument in \"%s\"\n", argv[opt.i - 1]);
			return 0;
		} 
        else if (c == '?') 
        {
			fprintf(stderr, "[ERROR] unknown option in \"%s\"\n", argv[opt.i - 1]);
			return 0;
		}
    }

    if (argc == 1)
    {
        Print_H(asm_opt);
        return 0;
    }

    get_queries(argc, argv, &opt, asm_opt);

    return check_option(asm_opt);
}
