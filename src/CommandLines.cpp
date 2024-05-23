#define __STDC_LIMIT_MACROS
#include <zlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <sys/time.h>
#include "CommandLines.h"
#include "ketopt.h"

#define DEFAULT_OUTPUT "hifiasm.asm"

hifiasm_opt_t asm_opt;

static ko_longopt_t long_options[] = {
	{ "version",       ko_no_argument, 300 },
	{ "dbg-gfa",       ko_no_argument, 301 },
	{ "write-paf",     ko_no_argument, 302 },
	{ "write-ec",      ko_no_argument, 303 },
	{ "skip-triobin",  ko_no_argument, 304 },
	{ "max-od-ec",     ko_no_argument, 305 },
	{ "max-od-final",  ko_no_argument, 306 },
	{ "ex-list",       ko_required_argument, 307 },
	{ "ex-iter",       ko_required_argument, 308 },
    { "hom-cov",     ko_required_argument, 309 },
    { "pri-range",     ko_required_argument, 310 },
    { "lowQ",          ko_required_argument, 312 },
	{ "min-hist-cnt",  ko_required_argument, 313 },
    { "h1",      ko_required_argument, 314 },
    { "h2",      ko_required_argument, 315 },
    { "enzyme",      ko_required_argument, 316 },
    { "b-cov",      ko_required_argument, 317 },
    { "h-cov",      ko_required_argument, 318 },
    { "m-rate",      ko_required_argument, 319 },
    { "primary",      ko_no_argument, 320 },
    { "t-occ",      ko_required_argument, 321 },
    { "seed",      ko_required_argument, 322 },
    { "n-perturb",      ko_required_argument, 323 },
    { "f-perturb",      ko_required_argument, 324 },
    { "n-hap",      ko_required_argument, 325 },
    { "n-weight",      ko_required_argument, 326 },
    { "l-msjoin",      ko_required_argument, 327 },
    { "purge-max",     ko_required_argument, 328 },
    { "fast",     ko_no_argument, 329 },
    { "dp-er",     ko_required_argument, 330},
    { "max-kocc",     ko_required_argument, 331},
    { "hg-size",     ko_required_argument, 332},
    { "ul",     ko_required_argument, 333},
    { "unskew",     ko_no_argument, 334},
    { "kpt-rate",     ko_required_argument, 335},
    { "ul-rate",     ko_required_argument, 336},
    { "dbg-het-cnt",     ko_no_argument, 337},
    { "ul-tip",     ko_required_argument, 338},
    { "low-het",     ko_no_argument, 339},
    { "s-base",     ko_required_argument, 340},
    { "bin-only",     ko_no_argument, 341},
    { "ul-round",     ko_required_argument, 342},
    { "prt-raw",     ko_no_argument, 343},
    { "integer-correct",     ko_required_argument, 344},
    { "dbg-ovec",     ko_no_argument, 345},
    { "path-max",     ko_required_argument, 346},
    { "path-min",     ko_required_argument, 347},
    { "trio-dual",     ko_no_argument, 348},
    { "ul-cut",     ko_required_argument, 349},
    { "dual-scaf",     ko_no_argument, 350},
    { "scaf-gap",     ko_required_argument, 351},
    { "sec-in",     ko_required_argument, 352},
    { "somatic-cov",     ko_required_argument, 353},
    { "telo-m",     ko_required_argument, 354},
    { "telo-p",     ko_required_argument, 355},
    { "telo-d",     ko_required_argument, 356},
    { "telo-s",     ko_required_argument, 357},
    { "ctg-n",     ko_required_argument, 358},
    // { "path-round",     ko_required_argument, 348},
	{ 0, 0, 0 }
};

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
	fprintf(stderr, "  Input/Output:\n");
    fprintf(stderr, "    -o STR       prefix of output files [%s]\n", asm_opt->output_file_name);
    fprintf(stderr, "    -t INT       number of threads [%d]\n", asm_opt->thread_num);
    fprintf(stderr, "    -h           show help information\n");
    fprintf(stderr, "    --version    show version number\n");
	fprintf(stderr, "  Overlap/Error correction:\n");
    fprintf(stderr, "    -k INT       k-mer length (must be <64) [%d]\n", asm_opt->k_mer_length);
	fprintf(stderr, "    -w INT       minimizer window size [%d]\n", asm_opt->mz_win);
	fprintf(stderr, "    -f INT       number of bits for bloom filter; 0 to disable [%d]\n", asm_opt->bf_shift);
	fprintf(stderr, "    -D FLOAT     drop k-mers occurring >FLOAT*coverage times [%.1f]\n", asm_opt->high_factor);
	fprintf(stderr, "    -N INT       consider up to max(-D*coverage,-N) overlaps for each oriented read [%d]\n", asm_opt->max_n_chain);
    fprintf(stderr, "    -r INT       round of correction [%d]\n", asm_opt->number_of_round);
    fprintf(stderr, "    -z INT       length of adapters that should be removed [%d]\n", asm_opt->adapterLen);
    fprintf(stderr, "    --max-kocc   INT\n");
    fprintf(stderr, "                 employ k-mers occurring <INT times to rescue repetitive overlaps [%d]\n", asm_opt->max_kmer_cnt);
    fprintf(stderr, "    --hg-size    INT(k, m or g)\n");
    fprintf(stderr, "                 estimated haploid genome size used for inferring read coverage [auto]\n");
    fprintf(stderr, "  Assembly:\n");
    fprintf(stderr, "    -a INT       round of assembly cleaning [%d]\n", asm_opt->clean_round);
    fprintf(stderr, "    -m INT       pop bubbles of <INT in size in contig graphs [%lld]\n", asm_opt->large_pop_bubble_size);
    fprintf(stderr, "    -p INT       pop bubbles of <INT in size in unitig graphs [%lld]\n", asm_opt->small_pop_bubble_size);
    fprintf(stderr, "    -n INT       remove tip unitigs composed of <=INT reads [%d]\n", asm_opt->max_short_tip);
    fprintf(stderr, "    -x FLOAT     max overlap drop ratio [%.2g]\n", asm_opt->max_drop_rate);
    fprintf(stderr, "    -y FLOAT     min overlap drop ratio [%.2g]\n", asm_opt->min_drop_rate);
    fprintf(stderr, "    -i           ignore saved read correction and overlaps\n");
    fprintf(stderr, "    -u           post-join step for contigs which may improve N50; 0 to disable; 1 to enable\n");
    fprintf(stderr, "                 [%u] and [%u] in default for the UL+HiFi assembly and the HiFi assembly, respectively\n",
                                      asm_opt->ul_pst_join, asm_opt->hifi_pst_join);
    fprintf(stderr, "    --hom-cov    INT\n");
    fprintf(stderr, "                 homozygous read coverage [auto]\n");
    fprintf(stderr, "    --lowQ       INT\n");
    fprintf(stderr, "                 output contig regions with >=INT%% inconsistency in BED format; 0 to disable [%d]\n", asm_opt->bed_inconsist_rate);
    fprintf(stderr, "    --b-cov      INT\n");
    fprintf(stderr, "                 break contigs at positions with <INT-fold coverage; work with '--m-rate'; 0 to disable [%d]\n", asm_opt->b_low_cov);
    fprintf(stderr, "    --h-cov      INT\n");
    fprintf(stderr, "                 break contigs at positions with >INT-fold coverage; work with '--m-rate'; -1 to disable [%d]\n", asm_opt->b_high_cov);
    fprintf(stderr, "    --m-rate     FLOAT\n");
    fprintf(stderr, "                 break contigs at positions with <=FLOAT*coverage exact overlaps;\n");
    fprintf(stderr, "                 only work with '--b-cov' or '--h-cov'[%.2f]\n", asm_opt->m_rate);
    fprintf(stderr, "    --primary    output a primary assembly and an alternate assembly\n");
    fprintf(stderr, "    --ctg-n      INT\n");
    fprintf(stderr, "                 remove tip contigs composed of <=INT reads [%d]\n", asm_opt->max_contig_tip);
    
    
//	fprintf(stderr, "    --pri-range INT1[,INT2]\n");
//	fprintf(stderr, "                keep contigs with coverage in this range in p_ctg.gfa; -1 to disable [auto,inf]\n");

    fprintf(stderr, "  Trio-partition:\n");
    fprintf(stderr, "    -1 FILE      hap1/paternal k-mer dump generated by \"yak count\" []\n");
    fprintf(stderr, "    -2 FILE      hap2/maternal k-mer dump generated by \"yak count\" []\n");
    fprintf(stderr, "    -3 FILE      list of hap1/paternal read names []\n");
	fprintf(stderr, "    -4 FILE      list of hap2/maternal read names []\n");
    fprintf(stderr, "    -c INT       lower bound of the binned k-mer's frequency [%d]\n", asm_opt->min_cnt);
    fprintf(stderr, "    -d INT       upper bound of the binned k-mer's frequency [%d]\n", asm_opt->mid_cnt);
    fprintf(stderr, "    --t-occ      INT\n");
    fprintf(stderr, "                 forcedly remove unitigs with >INT unexpected haplotype-specific reads;\n");
    fprintf(stderr, "                 ignore graph topology; [%d]\n", asm_opt->trio_flag_occ_thres);
    fprintf(stderr, "    --trio-dual  utilize homology information to correct trio phasing errors\n");


    fprintf(stderr, "  Purge-dups:\n");
    fprintf(stderr, "    -l INT       purge level. 0: no purging; 1: light; 2/3: aggressive [0 for trio; 3 for unzip]\n");
    fprintf(stderr, "    -s FLOAT     similarity threshold for duplicate haplotigs in read-level [%g for -l1/-l2, %g for -l3]\n", 
                                      asm_opt->purge_simi_rate_l2, asm_opt->purge_simi_rate_l3);
    fprintf(stderr, "    -O INT       min number of overlapped reads for duplicate haplotigs [%d]\n", 
                                      asm_opt->purge_overlap_len);
    fprintf(stderr, "    --purge-max  INT\n");
    fprintf(stderr, "                 coverage upper bound of Purge-dups [auto]\n");
    fprintf(stderr, "    --n-hap      INT\n");
    fprintf(stderr, "                 number of haplotypes [%d]\n", asm_opt->polyploidy);
    
    // fprintf(stderr, "  Hi-C-partition [experimental, not stable]:\n");
    fprintf(stderr, "  Hi-C-partition:\n");
    fprintf(stderr, "    --h1 FILEs   file names of Hi-C R1  [r1_1.fq,r1_2.fq,...]\n");
    fprintf(stderr, "    --h2 FILEs   file names of Hi-C R2  [r2_1.fq,r2_2.fq,...]\n");
    fprintf(stderr, "    --seed INT   RNG seed [%lu]\n", asm_opt->seed);
    fprintf(stderr, "    --s-base     FLOAT\n");
	fprintf(stderr, "                 similarity threshold for homology detection in base-level;\n"); 
    fprintf(stderr, "                 -1 to disable [%.3g]; -s for read-level (see <Purge-dups>)\n", asm_opt->trans_base_rate_sec);
    
    fprintf(stderr, "    --n-weight   INT\n");
    fprintf(stderr, "                 rounds of reweighting Hi-C links [%d]\n", asm_opt->n_weight);
	fprintf(stderr, "    --n-perturb  INT\n");
    fprintf(stderr, "                 rounds of perturbation [%d]\n", asm_opt->n_perturb);
    fprintf(stderr, "    --f-perturb  FLOAT\n");
	fprintf(stderr, "                 fraction to flip for perturbation [%.3g]\n", asm_opt->f_perturb);
    fprintf(stderr, "    --l-msjoin   INT\n");
    fprintf(stderr, "                 detect misjoined unitigs of >=INT in size; 0 to disable [%lu]\n", asm_opt->misjoin_len);

    fprintf(stderr, "  Ultra-Long-integration:\n");
    fprintf(stderr, "    --ul FILEs   file names of Ultra-Long reads [r1.fq,r2.fq,...]\n");
    fprintf(stderr, "    --ul-rate    FLOAT\n");
    fprintf(stderr, "                 error rate of Ultra-Long reads [%.3g]\n", asm_opt->ul_error_rate);
    fprintf(stderr, "    --ul-tip     INT\n");
    fprintf(stderr, "                 remove tip unitigs composed of <=INT reads for the UL assembly [%d]\n", asm_opt->max_short_ul_tip);
    fprintf(stderr, "    --path-max   FLOAT\n");
    fprintf(stderr, "                 max path drop ratio [%.2g]; higher number may make the assembly cleaner\n", asm_opt->max_path_drop_rate);
    fprintf(stderr, "                 but may lead to more misassemblies\n");
    fprintf(stderr, "    --path-min   FLOAT\n");
    fprintf(stderr, "                 min path drop ratio [%.2g]; higher number may make the assembly cleaner\n", asm_opt->min_path_drop_rate);
    fprintf(stderr, "                 but may lead to more misassemblies\n");
    fprintf(stderr, "    --ul-cut     INT\n");
    fprintf(stderr, "                 filter out <INT UL reads during the UL assembly [%d]\n", asm_opt->ul_min_base);
    // fprintf(stderr, "    --low-het    enable it for genomes with very low het heterozygosity rate (<0.0001%%)\n");

    fprintf(stderr, "  Dual-Scaffolding:\n");
    fprintf(stderr, "    --dual-scaf  output scaffolding\n");
    fprintf(stderr, "    --scaf-gap   INT\n");
    fprintf(stderr, "                 max gap size for scaffolding [%ld]\n", asm_opt->self_scaf_gap_max);

    fprintf(stderr, "  Telomere-identification:\n");
    fprintf(stderr, "    --telo-m     STR\n");
    fprintf(stderr, "                 telomere motif at 5'-end; CCCTAA for human [%s]\n", ((asm_opt->telo_motif)?(asm_opt->telo_motif):("NULL")));///5'-end, check CCCTAA
    fprintf(stderr, "    --telo-p     INT\n");
    fprintf(stderr, "                 non-telomeric penalty [%ld]\n", asm_opt->telo_pen);
    fprintf(stderr, "    --telo-d     INT\n");
    fprintf(stderr, "                 max drop [%ld]\n", asm_opt->telo_drop);
    fprintf(stderr, "    --telo-s     INT\n");
    fprintf(stderr, "                 min score for telomere reads [%ld]\n", asm_opt->telo_mic_sc);


    fprintf(stderr, "Example: ./hifiasm -o NA12878.asm -t 32 NA12878.fq.gz\n");
    fprintf(stderr, "See `https://hifiasm.readthedocs.io/en/latest/' or `man ./hifiasm.1' for complete documentation.\n");
}

void init_opt(hifiasm_opt_t* asm_opt)
{
	memset(asm_opt, 0, sizeof(hifiasm_opt_t));
	///asm_opt->flag = 0;
    asm_opt->flag = HA_F_PARTITION;
    asm_opt->coverage = -1;
    asm_opt->num_reads = 0;
    asm_opt->read_file_names = NULL;
    asm_opt->output_file_name = (char*)(DEFAULT_OUTPUT);
    asm_opt->required_read_name = NULL;
    asm_opt->hic_enzymes = NULL;
    asm_opt->hic_reads[0] = NULL;
    asm_opt->hic_reads[1] = NULL;
    asm_opt->fn_bin_poy = NULL;
    asm_opt->ar = NULL;
    asm_opt->thread_num = 1;
    asm_opt->k_mer_length = 51;
    asm_opt->hic_mer_length = 31;
    asm_opt->ul_mer_length = 19;
    asm_opt->trans_mer_length = 31;
	asm_opt->mz_win = 51;
    asm_opt->ul_mz_win = 19;
    asm_opt->trans_win = 31;
    asm_opt->mz_rewin = 1000;
    asm_opt->ul_mz_rewin = 360;
	asm_opt->mz_sample_dist = 500;
	asm_opt->bf_shift = 37;
	asm_opt->max_kmer_cnt = 2000;
	asm_opt->high_factor = 5.0;
	asm_opt->max_ov_diff_ec = 0.04;
	asm_opt->max_ov_diff_final = 0.03;
	asm_opt->hom_cov = 20;
    asm_opt->het_cov = -1024;
	asm_opt->max_n_chain = MIN_N_CHAIN;
	asm_opt->min_hist_kmer_cnt = 5;
    asm_opt->load_index_from_disk = 1;
    asm_opt->write_index_to_disk = 1;
    asm_opt->number_of_round = 3;
    asm_opt->adapterLen = 0;
    asm_opt->clean_round = 4;
    ///asm_opt->small_pop_bubble_size = 100000;
    asm_opt->small_pop_bubble_size = 0;
    asm_opt->large_pop_bubble_size = 10000000;
    asm_opt->min_drop_rate = 0.2;
    asm_opt->max_drop_rate = 0.8;
    asm_opt->max_hang_Len = 1000;
    asm_opt->max_hang_rate = 0.8;
    asm_opt->gap_fuzz = 1000;
    asm_opt->min_overlap_Len = 50;
    asm_opt->min_overlap_coverage = 0;
    asm_opt->max_short_tip = 3;
    asm_opt->max_short_ul_tip = 6;
    asm_opt->max_contig_tip = 3;
    asm_opt->min_cnt = 2;
    asm_opt->mid_cnt = 5;
    asm_opt->purge_level_primary = 3;
    asm_opt->purge_level_trio = 0;
    asm_opt->purge_simi_rate_l2 = 0.75;
    asm_opt->purge_simi_rate_l3 = 0.55;
    asm_opt->trans_base_rate = 0.93;
    asm_opt->trans_base_rate_sec = 0.5;
    asm_opt->purge_overlap_len = 1;
    ///asm_opt->purge_overlap_len_hic = 50;
    asm_opt->recover_atg_cov_min = -1024;
    asm_opt->recover_atg_cov_max = INT_MAX;
    asm_opt->hom_global_coverage = -1;
    asm_opt->hom_global_coverage_set = 0;
    asm_opt->pur_global_coverage = -1;
    asm_opt->bed_inconsist_rate = 70;
    asm_opt->hic_inconsist_rate = 30;
    ///asm_opt->bub_mer_length = 3;
    asm_opt->bub_mer_length = 1000000;
    asm_opt->b_low_cov = 0;
    asm_opt->b_high_cov = -1;
    asm_opt->m_rate = 0.75;
    asm_opt->hap_occ = 1;
    asm_opt->polyploidy = 2;
    asm_opt->trio_flag_occ_thres = 60;
    asm_opt->seed = 11;
    asm_opt->n_perturb = 10000;
    asm_opt->f_perturb = 0.1;
    asm_opt->n_weight = 3;
    asm_opt->is_alt = 0;
    asm_opt->misjoin_len = 500000;
    asm_opt->scffold = 0;
    asm_opt->dp_min_len = 2000;
    asm_opt->dp_e = 0.0025;
    asm_opt->hg_size = -1;
    asm_opt->kpt_rate = -1;
    asm_opt->infor_cov = 3;
    asm_opt->s_hap_cov = 3;
    asm_opt->ul_error_rate = 0.2/**0.15**/;
    asm_opt->ul_error_rate_low = 0.1;
    asm_opt->ul_error_rate_hpc = 0.2;
    asm_opt->ul_ec_round = 3;
    asm_opt->is_dbg_het_cnt = 0;
    asm_opt->is_low_het_ul = 0;
    asm_opt->is_base_trans = 1;
    asm_opt->is_read_trans = 1;
    asm_opt->is_topo_trans = 1;
    asm_opt->is_bub_trans = 1;
    asm_opt->bin_only = 0;
    asm_opt->ul_clean_round = 1;
    asm_opt->prt_dbg_gfa = 0;
    asm_opt->integer_correct_round = 0;
    asm_opt->dbg_ovec_cal = 0;
    asm_opt->min_path_drop_rate = 0.2;
    asm_opt->max_path_drop_rate = 0.6;
    asm_opt->hifi_pst_join = 1;
    asm_opt->ul_pst_join = 1;
    asm_opt->trio_cov_het_ovlp = -1;
    asm_opt->ul_min_base = 0;
    asm_opt->self_scaf = 0;
    asm_opt->self_scaf_min = 250000;
    asm_opt->self_scaf_reliable_min = 5000000;
    asm_opt->self_scaf_gap_max = 3000000;
    asm_opt->sec_in = NULL;
    asm_opt->somatic_cov = -1;

    asm_opt->telo_motif = NULL;
    asm_opt->telo_pen = 1;
    asm_opt->telo_drop = 2000;
    asm_opt->telo_mic_sc = 500;
}

void destory_enzyme(enzyme* f)
{
    int i;
    if(f != NULL)
    {
        for (i = 0; i < f->n; i++)
        {
            free(f->a[i]);
        }
        free(f->a);
        free(f->l);
        free(f);
    }
}

void destory_opt(hifiasm_opt_t* asm_opt)
{
    if(asm_opt->read_file_names != NULL) free(asm_opt->read_file_names);
    if(asm_opt->hic_enzymes != NULL) destory_enzyme(asm_opt->hic_enzymes);
    if(asm_opt->hic_reads[0] != NULL) destory_enzyme(asm_opt->hic_reads[0]);
    if(asm_opt->hic_reads[1] != NULL) destory_enzyme(asm_opt->hic_reads[1]);
    if(asm_opt->ar != NULL) destory_enzyme(asm_opt->ar);
}

void ha_opt_reset_to_round(hifiasm_opt_t* asm_opt, int round)
{
    asm_opt->num_bases = 0;
    asm_opt->num_corrected_bases = 0;
    asm_opt->num_recorrected_bases = 0;
	asm_opt->mem_buf = 0;
    asm_opt->roundID = round;
}

void ha_opt_update_cov(hifiasm_opt_t *opt, int hom_cov)
{
	int max_n_chain = (int)(hom_cov * opt->high_factor + .499);
	opt->hom_cov = hom_cov;
	if (opt->max_n_chain < max_n_chain)
		opt->max_n_chain = max_n_chain;
	fprintf(stderr, "[M::%s] updated max_n_chain to %d\n", __func__, opt->max_n_chain);
}

void ha_opt_update_cov_min(hifiasm_opt_t *opt, int hom_cov, int min_chain)
{
	int max_n_chain = (int)(hom_cov * opt->high_factor + .499);
	opt->hom_cov = hom_cov; opt->max_n_chain = max_n_chain;
    if(opt->max_n_chain < min_chain) opt->max_n_chain = min_chain;
	fprintf(stderr, "[M::%s] updated max_n_chain to %d\n", __func__, opt->max_n_chain);
}

static int check_file(char* name, const char* opt)
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

static int check_hic_reads(enzyme* f, const char* opt)
{
    int i;
    for (i = 0; i < f->n; i++)
    {
        if(check_file(f->a[i], opt) == 0) return 0;
    }
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

	if (asm_opt->max_ov_diff_ec < asm_opt->max_ov_diff_final) {
		fprintf(stderr, "[ERROR] max_ov_diff_ec shouldn't be smaller than max_ov_diff_final\n");
		return 0;
	}

	if (asm_opt->max_ov_diff_ec < HA_MIN_OV_DIFF) {
		fprintf(stderr, "[ERROR] max_ov_diff_ec shouldn't be smaller than %g\n", HA_MIN_OV_DIFF);
		return 0;
	}

    if(asm_opt->max_short_tip < 0)
    {
        fprintf(stderr, "[ERROR] the length of removal tips must be >= 0 (-n)\n");
        return 0;
    }

    if(asm_opt->purge_level_primary < 0 || asm_opt->purge_level_primary > 3)
    {
        fprintf(stderr, "[ERROR] the level of purge-dup should be [0, 3] (-l)\n");
        return 0;
    }

    if(ha_opt_triobin(asm_opt) && ((asm_opt->purge_level_trio < 0 || asm_opt->purge_level_trio > 1)))
    {
        fprintf(stderr, "[ERROR] the level of purge-dup for trio should be [0, 1] (-l)\n");
        return 0;
    }

    if(asm_opt->hom_global_coverage < 0 && asm_opt->hom_global_coverage != -1)
    {
        fprintf(stderr, "[ERROR] homozygous read coverage should be >= 0 (--hom-cov)\n");
        return 0;
    }

    if(asm_opt->pur_global_coverage < 0 && asm_opt->pur_global_coverage != -1)
    {
        fprintf(stderr, "[ERROR] purge duplication coverage threshold should be >= 0 (--purge-max)\n");
        return 0;
    }

    if(asm_opt->bed_inconsist_rate < 0 || asm_opt->bed_inconsist_rate > 100)
    {
        fprintf(stderr, "[ERROR] inconsistency rate should be [0, 100] (--lowQ)\n");
        return 0;
    }


    if(asm_opt->fn_bin_yak[0] != NULL && check_file(asm_opt->fn_bin_yak[0], "YAK1") == 0) return 0;
    if(asm_opt->fn_bin_yak[1] != NULL && check_file(asm_opt->fn_bin_yak[1], "YAK2") == 0) return 0;
    if(asm_opt->fn_bin_list[0] != NULL && check_file(asm_opt->fn_bin_list[0], "LIST1") == 0) return 0;
    if(asm_opt->fn_bin_list[1] != NULL && check_file(asm_opt->fn_bin_list[1], "LIST2") == 0) return 0;
    if(asm_opt->required_read_name != NULL && check_file(asm_opt->required_read_name, "b") == 0) return 0;
    
    if(asm_opt->hic_reads[0] != NULL && check_hic_reads(asm_opt->hic_reads[0], "HIC1") == 0) return 0;
    if(asm_opt->hic_reads[1] != NULL && check_hic_reads(asm_opt->hic_reads[1], "HIC2") == 0) return 0;
    if(asm_opt->hic_reads[0] != NULL && asm_opt->hic_reads[1] == NULL)
    {
        fprintf(stderr, "[ERROR] lack r2 of HiC reads (--h2)\n");
        return 0;
    }
    if(asm_opt->hic_reads[1] != NULL && asm_opt->hic_reads[0] == NULL)
    {
        fprintf(stderr, "[ERROR] lack r1 of HiC reads (--h1)\n");
        return 0;
    }

    if(asm_opt->hic_reads[0] != NULL && asm_opt->hic_reads[1] != NULL && 
                            asm_opt->hic_reads[0]->n != asm_opt->hic_reads[1]->n)
    {
        fprintf(stderr, "[ERROR] wrong r1 and r2 of HiC reads (--h1 && --h2)\n");
        return 0;
    }

    if(asm_opt->hic_enzymes != NULL && asm_opt->hic_enzymes->n == 0)
    {
        fprintf(stderr, "[ERROR] wrong HiC enzymes (--enzyme)\n");
        return 0;
    }

    if(asm_opt->hic_reads[0] != NULL && asm_opt->hic_reads[0]->n == 0)
    {
        fprintf(stderr, "[ERROR] wrong r1 of HiC reads (--h1)\n");
        return 0;
    }

    if(asm_opt->hic_reads[1] != NULL && asm_opt->hic_reads[1]->n == 0)
    {
        fprintf(stderr, "[ERROR] wrong r2 of HiC reads (--h2)\n");
        return 0;
    }

    if(asm_opt->ar != NULL && check_hic_reads(asm_opt->ar, "UL") == 0) return 0;
    if(asm_opt->ar != NULL && asm_opt->ar->n == 0)
    {
        fprintf(stderr, "[ERROR] wrong UL reads (--ul)\n");
        return 0;
    }

    if(asm_opt->b_low_cov < 0)
    {
        fprintf(stderr, "[ERROR] must >= 0 (--b-cov)\n");
        return 0;
    }

    if(asm_opt->b_high_cov != -1 && asm_opt->b_high_cov < 0)
    {
        fprintf(stderr, "[ERROR] must >= 0 (--h-cov)\n");
        return 0;
    }

    if(asm_opt->m_rate < 0)
    {
        fprintf(stderr, "[ERROR] must >= 0 (--m-rate)\n");
        return 0;
    }

    if(asm_opt->b_high_cov != -1 && asm_opt->b_high_cov <= asm_opt->b_low_cov)
    {
        fprintf(stderr, "[ERROR] [--h-cov] must >= [--b-cov]\n");
        return 0;
    }

    if(asm_opt->purge_simi_thres < 0)
    {
        fprintf(stderr, "[ERROR] [-s] must >= 0\n");
        return 0;
    }

    if(asm_opt->max_kmer_cnt < 0)
    {
        fprintf(stderr, "[ERROR] [--max-kocc] must >= 0\n");
        return 0;
    }

    if(asm_opt->hg_size < -1)
    {
        fprintf(stderr, "[ERROR] [--hg-size] wrong genome size\n");
        return 0;
    }

    if(asm_opt->telo_motif) {
        uint64_t k, tlen = strlen((asm_opt->telo_motif)); char c;
        if(tlen > 32) {
            fprintf(stderr, "[ERROR] [--telo-m] must be no longer than 32\n");
            return 0;
        }
        for (k = 0; k < tlen; k++) {
            c = asm_opt->telo_motif[k];
            if(c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
                c != 'a' && c != 'c' && c != 'g' && c != 't') {
                fprintf(stderr, "[ERROR] [--telo-m] must be A/C/G/T\n");
                return 0;
            }
        }
    }

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

void get_hic_enzymes(char *argv, enzyme** x, int check_name)
{
    int i, k, pre_i, len = strlen(argv);
    (*x) = (enzyme*)calloc(1, sizeof(enzyme));
    if(len == 0)
    {
        (*x)->n = 0; (*x)->l = NULL; (*x)->a = NULL;
        return;
    }


    (*x)->n = 1;
    for (i = pre_i = 0; i < len; i++)
    {
        if(argv[i] == ',')
        {
            (*x)->n++;
            continue;
        } 
        
        if(check_name)
        {
            if(argv[i] != 'A' && argv[i] != 'C' && argv[i] != 'G' && argv[i] != 'T' &&
                argv[i] != 'a' && argv[i] != 'c' && argv[i] != 'g' && argv[i] != 't' && 
                argv[i] != 'N' && argv[i] != 'n')
            {
                (*x)->n = 0;
                (*x)->l = NULL;
                (*x)->a = NULL;
                return;
            }
        }
        
    }
    (*x)->l = (int*)calloc((*x)->n, sizeof(int));
    (*x)->a = (char**)calloc((*x)->n, sizeof(char*));

    for (i = pre_i = k = 0; i < len; i++)
    {
        if(argv[i] == ',')
        {
            (*x)->l[k] = i - pre_i;
            (*x)->a[k] = (char*)malloc(sizeof(char)*((*x)->l[k]+1));
            memcpy((*x)->a[k], argv + pre_i, (*x)->l[k]);
            (*x)->a[k][(*x)->l[k]] = '\0';
            pre_i = i + 1;
            k++;
        }
    }

    (*x)->l[k] = i - pre_i;
    (*x)->a[k] = (char*)malloc(sizeof(char)*((*x)->l[k]+1));
    memcpy((*x)->a[k], argv + pre_i, (*x)->l[k]);
    (*x)->a[k][(*x)->l[k]] = '\0';
}

int64_t inter_gsize(char *argv)
{
    int64_t len = strlen(argv);
    double s;
    if(len <= 1) return -2;
    char t = argv[len-1];
    if(t != 'k' && t != 'K' && t != 'm' && t != 'M' && t != 'g' && t != 'G') return -2;
    char *ss=(char*)malloc(len);
    memcpy(ss, argv, len-1); ss[len-1] = '\0';
    s = atof(ss);
    free(ss);
    if(t == 'k' || t == 'K') return s*1000;
    if(t == 'm' || t == 'M') return s*1000000;
    if(t == 'g' || t == 'G') return s*1000000000;
    return s;
}

int CommandLine_process(int argc, char *argv[], hifiasm_opt_t* asm_opt)
{
    ketopt_t opt = KETOPT_INIT;

    int c;

    while ((c = ketopt(&opt, argc, argv, 1, "hvt:o:k:w:m:n:r:a:b:z:x:y:p:c:d:M:P:if:D:FN:1:2:3:4:5:l:s:O:eu:", long_options)) >= 0) {
        if (c == 'h')
        {
            Print_H(asm_opt);
            return 0;
        } 
        else if (c == 'v' || c == 300)
        {
			puts(HA_VERSION);
            return 0;
        }
		else if (c == 'f') asm_opt->bf_shift = atoi(opt.arg);
        else if (c == 't') asm_opt->thread_num = atoi(opt.arg); 
        else if (c == 'o') asm_opt->output_file_name = opt.arg;
        else if (c == 'r') asm_opt->number_of_round = atoi(opt.arg);
        else if (c == 'k') asm_opt->k_mer_length = atoi(opt.arg);
        else if (c == 'i') asm_opt->load_index_from_disk = 0; 
        else if (c == 'w') asm_opt->mz_win = atoi(opt.arg);
		else if (c == 'D') asm_opt->high_factor = atof(opt.arg);
		else if (c == 'F') asm_opt->flag |= HA_F_NO_KMER_FLT;
		else if (c == 'N') asm_opt->max_n_chain = atoi(opt.arg);
        else if (c == 'a') asm_opt->clean_round = atoi(opt.arg); 
        else if (c == 'z') asm_opt->adapterLen = atoi(opt.arg);
        else if (c == 'b') asm_opt->required_read_name = opt.arg;
        else if (c == 'c') asm_opt->min_cnt = atoi(opt.arg);
        else if (c == 'd') asm_opt->mid_cnt = atoi(opt.arg);
        else if (c == '1' || c == 'P') asm_opt->fn_bin_yak[0] = opt.arg; // -P/-M reserved for backward compatibility
        else if (c == '2' || c == 'M') asm_opt->fn_bin_yak[1] = opt.arg;
        else if (c == '3') asm_opt->fn_bin_list[0] = opt.arg;
        else if (c == '4') asm_opt->fn_bin_list[1] = opt.arg;
        else if (c == '5') asm_opt->fn_bin_poy = opt.arg;
        else if (c == 'x') asm_opt->max_drop_rate = atof(opt.arg);
        else if (c == 'y') asm_opt->min_drop_rate = atof(opt.arg);
        else if (c == 'p') asm_opt->small_pop_bubble_size = atoll(opt.arg);
        else if (c == 'm') asm_opt->large_pop_bubble_size = atoll(opt.arg);
        else if (c == 'n') asm_opt->max_short_tip = atoll(opt.arg);
        else if (c == 'e') asm_opt->flag |= HA_F_BAN_ASSEMBLY;
        else if (c == 'u') {
            if(atoll(opt.arg)) {
                asm_opt->hifi_pst_join = asm_opt->ul_pst_join = 1;
            } else {
                asm_opt->hifi_pst_join = asm_opt->ul_pst_join = 0;
            }
        }
        else if (c == 301) asm_opt->flag |= HA_F_VERBOSE_GFA;
		else if (c == 302) asm_opt->flag |= HA_F_WRITE_PAF;
		else if (c == 303) asm_opt->flag |= HA_F_WRITE_EC;
		else if (c == 304) asm_opt->flag |= HA_F_SKIP_TRIOBIN;
		else if (c == 305) asm_opt->max_ov_diff_ec = atof(opt.arg);
		else if (c == 306) asm_opt->max_ov_diff_final = atof(opt.arg);
		else if (c == 307) asm_opt->extract_list = opt.arg;
		else if (c == 308) asm_opt->extract_iter = atoi(opt.arg);
        else if (c == 309)
        {
            asm_opt->hom_global_coverage = atoi(opt.arg);
            asm_opt->hom_global_coverage_set = 1;
        } 
        else if (c == 310)
        {
            char* s = NULL;
            asm_opt->recover_atg_cov_min = strtol(opt.arg, &s, 10);
			if (*s == ',') asm_opt->recover_atg_cov_max = strtol(s + 1, &s, 10);
            if(asm_opt->recover_atg_cov_min == -1 || asm_opt->recover_atg_cov_max == -1)
            {
                asm_opt->recover_atg_cov_min = asm_opt->recover_atg_cov_max = -1;
            }
        }
        ///else if (c == 311) asm_opt->flag |= HA_F_HIGH_HET;
        else if (c == 312) asm_opt->bed_inconsist_rate = atoi(opt.arg);
        else if (c == 313) asm_opt->min_hist_kmer_cnt = atoi(opt.arg);
        else if (c == 314) get_hic_enzymes(opt.arg, &(asm_opt->hic_reads[0]), 0);
        else if (c == 315) get_hic_enzymes(opt.arg, &(asm_opt->hic_reads[1]), 0);
        else if (c == 316) get_hic_enzymes(opt.arg, &(asm_opt->hic_enzymes), 1);
        else if (c == 317) asm_opt->b_low_cov = atoi(opt.arg);
        else if (c == 318) asm_opt->b_high_cov = atoi(opt.arg);
        else if (c == 319) asm_opt->m_rate = atof(opt.arg);
        else if (c == 320) asm_opt->flag -= HA_F_PARTITION, asm_opt->is_alt = 1;
        else if (c == 321) asm_opt->trio_flag_occ_thres = atoi(opt.arg);
        else if (c == 322) asm_opt->seed = atol(opt.arg);
        else if (c == 323) asm_opt->n_perturb = atoi(opt.arg);
        else if (c == 324) asm_opt->f_perturb = atof(opt.arg);
        else if (c == 325) asm_opt->polyploidy = atoi(opt.arg);
        else if (c == 326) asm_opt->n_weight = atoi(opt.arg);
        else if (c == 327) asm_opt->misjoin_len = atol(opt.arg);
        else if (c == 328) asm_opt->pur_global_coverage = atoi(opt.arg);
        else if (c == 329) asm_opt->flag |= HA_F_FAST;
        else if (c == 330) asm_opt->dp_e = atof(opt.arg);
        else if (c == 331) asm_opt->max_kmer_cnt = atol(opt.arg);
        else if (c == 332) asm_opt->hg_size = inter_gsize(opt.arg);      
        else if (c == 333) get_hic_enzymes(opt.arg, &(asm_opt->ar), 0);
        else if (c == 334) asm_opt->flag |= HA_F_USKEW;
        else if (c == 335) asm_opt->kpt_rate = atof(opt.arg);
        else if (c == 336) asm_opt->ul_error_rate = atof(opt.arg);
        else if (c == 337) asm_opt->is_dbg_het_cnt = 1;
        else if (c == 338) asm_opt->max_short_ul_tip = atol(opt.arg);
        else if (c == 339) asm_opt->is_low_het_ul = 1;
        else if (c == 340) {
            asm_opt->trans_base_rate_sec = atof(opt.arg);
            if(asm_opt->trans_base_rate_sec < 0) asm_opt->is_base_trans = 0;
        } 
        else if (c == 341) asm_opt->bin_only = 1;
        else if (c == 342) asm_opt->ul_clean_round = atol(opt.arg);
        else if (c == 343) asm_opt->prt_dbg_gfa = 1;
        else if (c == 344) asm_opt->integer_correct_round = atol(opt.arg);
        else if (c == 345) asm_opt->dbg_ovec_cal = 1;
        else if (c == 346) asm_opt->max_path_drop_rate = atof(opt.arg);
        else if (c == 347) asm_opt->min_path_drop_rate = atof(opt.arg);
        else if (c == 348) asm_opt->trio_cov_het_ovlp = 1;
        else if (c == 349) asm_opt->ul_min_base = atol(opt.arg);
        else if (c == 350) asm_opt->self_scaf = 1;
        else if (c == 351) asm_opt->self_scaf_gap_max = atol(opt.arg);
        else if (c == 352) get_hic_enzymes(opt.arg, &(asm_opt->sec_in), 0);
        else if (c == 353) asm_opt->somatic_cov = atol(opt.arg);
        else if (c == 354) asm_opt->telo_motif = opt.arg;
        else if (c == 355) asm_opt->telo_pen = atol(opt.arg);
        else if (c == 356) asm_opt->telo_drop = atol(opt.arg);
        else if (c == 357) asm_opt->telo_mic_sc = atol(opt.arg);
        else if (c == 358) asm_opt->max_contig_tip = atol(opt.arg);
        else if (c == 'l') {   ///0: disable purge_dup; 1: purge containment; 2: purge overlap
            asm_opt->purge_level_primary = asm_opt->purge_level_trio = atoi(opt.arg);
        }
        else if (c == 's') asm_opt->purge_simi_rate_l2 = asm_opt->purge_simi_rate_l3 = atof(opt.arg);
        else if (c == 'O') asm_opt->purge_overlap_len = atoll(opt.arg);
        else if (c == ':') 
        {
			fprintf(stderr, "[ERROR] missing option argument in \"%s\"\n", argv[opt.i - 1]);
			return 1;
		} 
        else if (c == '?') 
        {
			fprintf(stderr, "[ERROR] unknown option in \"%s\"\n", argv[opt.i - 1]);
			return 1;
		}
    }
    
    if(asm_opt->purge_level_primary > 2) asm_opt->purge_simi_thres = asm_opt->purge_simi_rate_l3;
    else asm_opt->purge_simi_thres = asm_opt->purge_simi_rate_l2;
    
    if (argc == opt.ind)
    {
        Print_H(asm_opt);
        return 0;
    }

    get_queries(argc, argv, &opt, asm_opt);

    c = ((asm_opt->ar)?(asm_opt->ul_pst_join):(asm_opt->hifi_pst_join));
    if(c) {
        if((asm_opt->flag&HA_F_BAN_POST_JOIN)) asm_opt->flag-=HA_F_BAN_POST_JOIN;
    } else {
        asm_opt->flag |= HA_F_BAN_POST_JOIN;
    }

    // fprintf(stderr, "[M::%s::] post join::%u\n", __func__, (uint32_t)(!(asm_opt->flag & HA_F_BAN_POST_JOIN)));
    // exit(1);

    return check_option(asm_opt);
}
