#include "CommandLines.h"
#include <stdlib.h>
#include <stdio.h>
#include "ketopt.h"
#include <sys/time.h>

#define VERSION "0.0.0.1"

hifiasm_opt_t asm_opt;

double Get_T(void)
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec+t.tv_usec/1000000.0;
}

void Print_H(hifiasm_opt_t* asm_opt)
{
    fprintf(stderr, "Usage: hifiasm [options] -q <input.fa> -o <output_asm>\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -q FILE       input in the fastq(.gz)/fasta(.gz) formats\n");
    fprintf(stderr, "    -k FILE       output assembly (in gfa format) and corrected reads (in fasta format)\n");
    fprintf(stderr, "    -t INT        number of threads [%d]\n", asm_opt->thread_num);
    fprintf(stderr, "    -r INT        round of correction [%d]\n", asm_opt->number_of_round);
    fprintf(stderr, "    -a INT        round of assembly cleaning [%d]\n", asm_opt->clean_round);
    fprintf(stderr, "    -k INT        k-mer length [%d] (must be < 64)\n", asm_opt->k_mer_length);
    fprintf(stderr, "    -w            write all overlaps to disk, can accelerate assembly next time\n");
    fprintf(stderr, "    -l            load all overlaps from disk, can avoid overlap calculation\n");
    fprintf(stderr, "    -z INT        length of adapters that should be removed [%d]\n", asm_opt->adapterLen);
    fprintf(stderr, "    -p INT        size of popped bubbles [%lld]\n", asm_opt->pop_bubble_size);
    fprintf(stderr, "    -x FLOAT      max overlap drop ratio [%.2g]\n", asm_opt->max_drop_rate);
    fprintf(stderr, "    -y FLOAT      min overlap drop ratio [%.2g]\n", asm_opt->min_drop_rate);
    fprintf(stderr, "    -v            show version number\n");
    fprintf(stderr, "    -h            show help information\n");
    fprintf(stderr, "Example: ./hifiasm -w -l -q NA12878.fq.gz -o NA12878.asm -k 40 -t 32 -r 2 -a 4 -z 0\n");

    
}

void init_opt(hifiasm_opt_t* asm_opt)
{
    asm_opt->read_file_name = NULL;
    asm_opt->output_file_name = NULL;
    asm_opt->required_read_name = NULL;
    asm_opt->thread_num = 1;
    asm_opt->k_mer_length = 40;
    asm_opt->k_mer_min_freq = 3;
    asm_opt->k_mer_max_freq = 66;
    asm_opt->load_index_from_disk = 0;
    asm_opt->write_index_to_disk = 0;
    asm_opt->number_of_round = 2;
    asm_opt->adapterLen = 0;
    asm_opt->clean_round = 4;
    asm_opt->complete_threads = 0;
    asm_opt->pop_bubble_size = 100000;
    asm_opt->min_drop_rate = 0.2;
    asm_opt->max_drop_rate = 0.8;
}

void clear_opt(hifiasm_opt_t* asm_opt, int last_round)
{
    asm_opt->complete_threads = 0;
    asm_opt->num_bases = 0;
    asm_opt->num_corrected_bases = 0;
    asm_opt->num_recorrected_bases = 0;
    asm_opt->roundID = asm_opt->number_of_round - last_round;
}


int check_option(hifiasm_opt_t* asm_opt)
{
    if(asm_opt->read_file_name == NULL)
    {
        fprintf(stderr, "[ERROR] missing input: please specify a read file\n");
        return 0;
    }

    if(asm_opt->output_file_name == NULL)
    {
        fprintf(stderr, "[ERROR] missing output: please specify the output name\n");
        return 0;
    }

    if(asm_opt->thread_num < 1)
    {
        fprintf(stderr, "[ERROR] the number of threads must be > 0\n");
        return 0;
    }


    if(asm_opt->number_of_round < 1)
    {
        fprintf(stderr, "[ERROR] the number of rounds for correction must be > 0\n");
        return 0;
    }

    if(asm_opt->clean_round < 1)
    {
        fprintf(stderr, "[ERROR] the number of rounds for assembly cleaning must be > 0\n");
        return 0;
    }

    if(asm_opt->adapterLen < 0)
    {
        fprintf(stderr, "[ERROR] the length of removed adapters must be >= 0\n");
        return 0;
    }


    if(asm_opt->k_mer_length >= 64)
    {
        fprintf(stderr, "[ERROR] the length of k_mer must be < 64\n");
        return 0;
    }


    if(asm_opt->max_drop_rate < 0 || asm_opt->max_drop_rate >= 1 ) 
    {
        fprintf(stderr, "[ERROR] max overlap drop ratio must be [0.0, 1.0)\n");
        return 0;
    }

    
    if(asm_opt->min_drop_rate < 0 || asm_opt->min_drop_rate >= 1)
    {
        fprintf(stderr, "[ERROR] min overlap drop ratio must be [0.0, 1.0)\n");
        return 0;
    }

    if(asm_opt->max_drop_rate <= asm_opt->min_drop_rate)
    {
        fprintf(stderr, "[ERROR] min overlap drop ratio must be less than max overlap drop ratio\n");
        return 0;
    }

    if(asm_opt->pop_bubble_size < 0)
    {
        fprintf(stderr, "[ERROR] the size of popped bubbles must be >= 0\n");
        return 0;
    }


    // fprintf(stderr, "input file: %s\n", asm_opt->read_file_name);
    // fprintf(stderr, "output file: %s\n", asm_opt->output_file_name);
    // fprintf(stderr, "number of threads: %d\n", asm_opt->thread_num);
    // fprintf(stderr, "number of rounds for correction: %d\n", asm_opt->number_of_round);
    // fprintf(stderr, "number of rounds for assembly cleaning: %d\n", asm_opt->clean_round);
    // fprintf(stderr, "length of removed adapters: %d\n", asm_opt->adapterLen);
    // fprintf(stderr, "length of k_mer: %d\n", asm_opt->k_mer_length);
    // fprintf(stderr, "min overlap drop ratio: %.2g\n", asm_opt->min_drop_rate);
    // fprintf(stderr, "max overlap drop ratio: %.2g\n", asm_opt->max_drop_rate);
    // fprintf(stderr, "size of popped bubbles: %lld\n", asm_opt->pop_bubble_size);

    return 1;
}


int CommandLine_process(int argc, char *argv[], hifiasm_opt_t* asm_opt)
{
    ketopt_t opt = KETOPT_INIT;

    int c;

    while ((c = ketopt(&opt, argc, argv, 1, "hvt:o:q:k:lwm:n:r:a:b:z:x:y:p:", 0)) >= 0) {
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
        else if (c == 'q') asm_opt->read_file_name = opt.arg;
        else if (c == 'n') asm_opt->k_mer_min_freq = atoi(opt.arg);
        else if (c == 'm') asm_opt->k_mer_max_freq = atoi(opt.arg);
        else if (c == 'r') asm_opt->number_of_round = atoi(opt.arg);
        else if (c == 'k') asm_opt->k_mer_length = atoi(opt.arg); 
        else if (c == 'l') asm_opt->load_index_from_disk = 1; 
        else if (c == 'w') asm_opt->write_index_to_disk = 1;
        else if (c == 'a') asm_opt->clean_round = atoi(opt.arg); 
        else if (c == 'z') asm_opt->adapterLen = atoi(opt.arg);
        else if (c == 'b') asm_opt->required_read_name = opt.arg;
        else if (c == 'x') asm_opt->max_drop_rate = atof(opt.arg);
        else if (c == 'y') asm_opt->min_drop_rate = atof(opt.arg);
        else if (c == 'p') asm_opt->pop_bubble_size = atoll(opt.arg);
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


    return check_option(asm_opt);
}