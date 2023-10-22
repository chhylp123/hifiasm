#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include "Process_Read.h"
#include "htab.h"
#include "Correct.h"
#include "kalloc.h"
#include <assert.h>

#define UL_FLANK 512
uint8_t seq_nt6_table[256] = {
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 0, 5, 1,  5, 5, 5, 2,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  3, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 0, 5, 1,  5, 5, 5, 2,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  3, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

char bit_t_seq_table[256][4] = {{0}};
char bit_t_seq_table_rc[256][4] = {{0}};
char s_H[5] = {'A', 'C', 'G', 'T', 'N'};
char rc_Table[6] = {'T', 'G', 'C', 'A', 'N', 'N'};

void init_All_reads(All_reads* r)
{
	memset(r, 0, sizeof(All_reads));
	r->index_size = READ_INIT_NUMBER;
	r->read_length = (uint64_t*)malloc(sizeof(uint64_t)*r->index_size);
	r->name_index_size = READ_INIT_NUMBER;
	r->name_index = (uint64_t*)malloc(sizeof(uint64_t)*r->name_index_size);
	r->name_index[0] = 0;
}

void destory_All_reads(All_reads* r)
{
	uint64_t i = 0;
	for (i = 0; i < r->total_reads; i++) {
		if (r->N_site[i]) free(r->N_site[i]);
		if (r->read_sperate[i]) free(r->read_sperate[i]);
		if (r->paf && r->paf[i].buffer) free(r->paf[i].buffer);
		if (r->reverse_paf && r->reverse_paf[i].buffer) free(r->reverse_paf[i].buffer);
		///if (r->pb_regions) kv_destroy(r->pb_regions[i].a);
	}
	free(r->paf);
	free(r->reverse_paf);
	free(r->N_site);
	free(r->read_sperate);
	free(r->name);
	free(r->name_index);
	free(r->read_length);
	free(r->trio_flag);
	///if (r->pb_regions) free(r->pb_regions);
}

void write_All_reads(All_reads* r, char* read_file_name)
{
    fprintf(stderr, "Writing reads to disk... \n");
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "w");
	fwrite(&asm_opt.adapterLen, sizeof(asm_opt.adapterLen), 1, fp);
    fwrite(&r->index_size, sizeof(r->index_size), 1, fp);
	fwrite(&r->name_index_size, sizeof(r->name_index_size), 1, fp);
	fwrite(&r->total_reads, sizeof(r->total_reads), 1, fp);
	fwrite(&r->total_reads_bases, sizeof(r->total_reads_bases), 1, fp);
	fwrite(&r->total_name_length, sizeof(r->total_name_length), 1, fp);

	uint64_t i = 0;
	uint64_t zero = 0;
	for (i = 0; i < r->total_reads; i++)
	{
		if (r->N_site[i] != NULL)
		{
			///number of Ns
			fwrite(&r->N_site[i][0], sizeof(r->N_site[i][0]), 1, fp);
			if (r->N_site[i][0])
			{
				fwrite(r->N_site[i]+1, sizeof(r->N_site[i][0]), r->N_site[i][0], fp);
			}
		}
		else
		{
			fwrite(&zero, sizeof(zero), 1, fp);
		}
	}

	fwrite(r->read_length, sizeof(uint64_t), r->total_reads, fp);
	for (i = 0; i < r->total_reads; i++)
	{
		fwrite(r->read_sperate[i], sizeof(uint8_t), r->read_length[i]/4+1, fp);
	}
	
	fwrite(r->name, sizeof(char), r->total_name_length, fp);
	fwrite(r->name_index, sizeof(uint64_t), r->name_index_size, fp);
	fwrite(r->trio_flag, sizeof(uint8_t), r->total_reads, fp);
	fwrite(&(asm_opt.hom_cov), sizeof(asm_opt.hom_cov), 1, fp);
	fwrite(&(asm_opt.het_cov), sizeof(asm_opt.het_cov), 1, fp);

    free(index_name);    
	fflush(fp);
    fclose(fp);
    fprintf(stderr, "Reads has been written.\n");
}

int load_All_reads(All_reads* r, char* read_file_name)
{
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "r");
	if (!fp) {
		free(index_name);
        return 0;
    }
	int local_adapterLen;
	int f_flag;
    f_flag = fread(&local_adapterLen, sizeof(local_adapterLen), 1, fp);
    if(local_adapterLen != asm_opt.adapterLen)
    {
        fprintf(stderr, "the adapterLen of index is: %d, but the adapterLen set by user is: %d\n", 
        local_adapterLen, asm_opt.adapterLen);
		exit(1);
    }
    f_flag += fread(&r->index_size, sizeof(r->index_size), 1, fp);
	f_flag += fread(&r->name_index_size, sizeof(r->name_index_size), 1, fp);
	f_flag += fread(&r->total_reads, sizeof(r->total_reads), 1, fp);
	f_flag += fread(&r->total_reads_bases, sizeof(r->total_reads_bases), 1, fp);
	f_flag += fread(&r->total_name_length, sizeof(r->total_name_length), 1, fp);

	uint64_t i = 0;
	uint64_t zero = 0;
	r->N_site = (uint64_t**)malloc(sizeof(uint64_t*)*r->total_reads);
	for (i = 0; i < r->total_reads; i++)
	{
		f_flag += fread(&zero, sizeof(zero), 1, fp);

		if (zero)
		{
			r->N_site[i] = (uint64_t*)malloc(sizeof(uint64_t)*(zero + 1));
			r->N_site[i][0] = zero;
			if (r->N_site[i][0])
			{
				f_flag += fread(r->N_site[i]+1, sizeof(r->N_site[i][0]), r->N_site[i][0], fp);
			}
		}
		else
		{
			r->N_site[i] = NULL;
		}
	}

	r->read_length = (uint64_t*)malloc(sizeof(uint64_t)*r->total_reads);
	f_flag += fread(r->read_length, sizeof(uint64_t), r->total_reads, fp);

	r->read_size = (uint64_t*)malloc(sizeof(uint64_t)*r->total_reads);
	memcpy (r->read_size, r->read_length, sizeof(uint64_t)*r->total_reads);

	r->read_sperate = (uint8_t**)malloc(sizeof(uint8_t*)*r->total_reads);
	for (i = 0; i < r->total_reads; i++)
    {
		r->read_sperate[i] = (uint8_t*)malloc(sizeof(uint8_t)*(r->read_length[i]/4+1));
		f_flag += fread(r->read_sperate[i], sizeof(uint8_t), r->read_length[i]/4+1, fp);
	}


	r->name = (char*)malloc(sizeof(char)*r->total_name_length);
	f_flag += fread(r->name, sizeof(char), r->total_name_length, fp);

	r->name_index = (uint64_t*)malloc(sizeof(uint64_t)*r->name_index_size);
	f_flag += fread(r->name_index, sizeof(uint64_t), r->name_index_size, fp);

	/****************************may have bugs********************************/
	r->trio_flag = (uint8_t*)malloc(sizeof(uint8_t)*r->total_reads);
	f_flag += fread(r->trio_flag, sizeof(uint8_t), r->total_reads, fp);
	f_flag += fread(&(asm_opt.hom_cov), sizeof(asm_opt.hom_cov), 1, fp);
    f_flag += fread(&(asm_opt.het_cov), sizeof(asm_opt.het_cov), 1, fp);
	/****************************may have bugs********************************/

	r->cigars = (Compressed_Cigar_record*)malloc(sizeof(Compressed_Cigar_record)*r->total_reads);
	r->second_round_cigar = (Compressed_Cigar_record*)malloc(sizeof(Compressed_Cigar_record)*r->total_reads);
	for (i = 0; i < r->total_reads; i++)
	{
		r->second_round_cigar[i].size = r->cigars[i].size = 0;
		r->second_round_cigar[i].length = r->cigars[i].length = 0;
		r->second_round_cigar[i].record = r->cigars[i].record = NULL;

		r->second_round_cigar[i].lost_base_size = r->cigars[i].lost_base_size = 0;
		r->second_round_cigar[i].lost_base_length = r->cigars[i].lost_base_length = 0;
		r->second_round_cigar[i].lost_base = r->cigars[i].lost_base = NULL;
	}
	///r->pb_regions = NULL;

    free(index_name);    
    fclose(fp);
    fprintf(stderr, "Reads has been loaded.\n");

	return 1;
}

void read_ma(ma_hit_t* x, FILE* fp)
{
    int f_flag;
    f_flag = fread(&(x->qns), sizeof(x->qns), 1, fp);
    f_flag += fread(&(x->qe), sizeof(x->qe), 1, fp);
    f_flag += fread(&(x->tn), sizeof(x->tn), 1, fp);
    f_flag += fread(&(x->ts), sizeof(x->ts), 1, fp);
    f_flag += fread(&(x->te), sizeof(x->te), 1, fp);
    f_flag += fread(&(x->el), sizeof(x->el), 1, fp);
    f_flag += fread(&(x->no_l_indel), sizeof(x->no_l_indel), 1, fp);

    uint32_t t;
    f_flag += fread(&(t), sizeof(t), 1, fp);
    x->ml = t;

    f_flag += fread(&(t), sizeof(t), 1, fp);
    x->rev = t;
    
    f_flag += fread(&(t), sizeof(t), 1, fp);
    x->bl = t;

    f_flag += fread(&(t), sizeof(t), 1, fp);
    x->del = t;
}


int append_All_reads(All_reads* r, char *idx, uint32_t id)
{
	char *gfa_name = NULL; MALLOC(gfa_name, (strlen(idx)+100));
	FILE *fp = NULL, *fpo = NULL;

	sprintf(gfa_name, "%s.ec%u.bin", idx, id);
    fp = fopen(gfa_name, "r");
	if (!fp) {free(gfa_name); return 0;}

	sprintf(gfa_name, "%s.ovlp%u.source.bin", idx, id);
	fpo = fopen(gfa_name, "r");
	if (!fpo) {free(gfa_name); fclose(fp); return 0;}
	free(gfa_name);

	int local_adapterLen, f_flag;
    f_flag = fread(&local_adapterLen, sizeof(local_adapterLen), 1, fp);
    if(local_adapterLen != asm_opt.adapterLen) {
        fprintf(stderr, "the adapterLen of index is: %d, but the adapterLen set by user is: %d\n", 
        local_adapterLen, asm_opt.adapterLen);
		exit(1);
    }
	uint64_t index_size0, name_index_size0, total_reads0, total_reads_bases0, total_name_length0;
	index_size0 = r->index_size; 
	total_reads0 = r->total_reads; 
	total_reads_bases0 = r->total_reads_bases;
	total_name_length0 = r->total_name_length;
	name_index_size0 = ((total_reads0)?(total_reads0+1):(0));///r->name_index_size;

    f_flag += fread(&r->index_size, sizeof(r->index_size), 1, fp);
	f_flag += fread(&r->name_index_size, sizeof(r->name_index_size), 1, fp);
	f_flag += fread(&r->total_reads, sizeof(r->total_reads), 1, fp);
	f_flag += fread(&r->total_reads_bases, sizeof(r->total_reads_bases), 1, fp);
	f_flag += fread(&r->total_name_length, sizeof(r->total_name_length), 1, fp);

	r->index_size += index_size0;
	r->name_index_size += name_index_size0; if(name_index_size0) r->name_index_size--;
	r->total_reads += total_reads0;
	r->total_reads_bases += total_reads_bases0;
	r->total_name_length += total_name_length0;


	uint64_t i = 0, zero = 0, k;
	REALLOC(r->N_site, r->total_reads);

	for (i = total_reads0; i < r->total_reads; i++) {
		f_flag += fread(&zero, sizeof(zero), 1, fp);

		if (zero) {
			r->N_site[i] = (uint64_t*)malloc(sizeof(uint64_t)*(zero + 1));
			r->N_site[i][0] = zero;
			if (r->N_site[i][0]) {
				f_flag += fread(r->N_site[i]+1, sizeof(r->N_site[i][0]), r->N_site[i][0], fp);
			}
		} else {
			r->N_site[i] = NULL;
		}
	}

	REALLOC(r->read_length, r->total_reads);
	f_flag += fread(r->read_length + total_reads0, sizeof(uint64_t), r->total_reads-total_reads0, fp);

	REALLOC(r->read_size, r->total_reads);
	memcpy (r->read_size + total_reads0, r->read_length + total_reads0, sizeof(uint64_t)*(r->total_reads-total_reads0));

	REALLOC(r->read_sperate, r->total_reads);
	for (i = total_reads0; i < r->total_reads; i++) {
		r->read_sperate[i] = (uint8_t*)malloc(sizeof(uint8_t)*(r->read_length[i]/4+1));
		f_flag += fread(r->read_sperate[i], sizeof(uint8_t), r->read_length[i]/4+1, fp);
	}


	REALLOC(r->name, r->total_name_length);
	f_flag += fread(r->name + total_name_length0, sizeof(char), r->total_name_length - total_name_length0, fp);

	REALLOC(r->name_index, r->name_index_size); uint64_t sft = 0;
	if(name_index_size0) sft = r->name_index[--name_index_size0];
	// fprintf(stderr, "name_index_size0::%lu, total_reads0::%lu, sft::%lu\n", name_index_size0, total_reads0, sft);
	f_flag += fread(r->name_index + name_index_size0, sizeof(uint64_t), r->name_index_size-name_index_size0, fp);
	if(sft) {
		for (i = name_index_size0; i < r->name_index_size; i++) r->name_index[i] += sft;
	}

	/****************************may have bugs********************************/
	REALLOC(r->trio_flag, r->total_reads);
	f_flag += fread(r->trio_flag + total_reads0, sizeof(uint8_t), r->total_reads - total_reads0, fp);

	int hom_cov0 = asm_opt.hom_cov, het_cov0 = asm_opt.het_cov;
	f_flag += fread(&(asm_opt.hom_cov), sizeof(asm_opt.hom_cov), 1, fp); asm_opt.hom_cov += hom_cov0;
    f_flag += fread(&(asm_opt.het_cov), sizeof(asm_opt.het_cov), 1, fp); asm_opt.het_cov += het_cov0;
	/****************************may have bugs********************************/

	REALLOC(r->cigars, r->total_reads);
	REALLOC(r->second_round_cigar, r->total_reads);
	for (i = total_reads0; i < r->total_reads; i++) {
		r->second_round_cigar[i].size = r->cigars[i].size = 0;
		r->second_round_cigar[i].length = r->cigars[i].length = 0;
		r->second_round_cigar[i].record = r->cigars[i].record = NULL;

		r->second_round_cigar[i].lost_base_size = r->cigars[i].lost_base_size = 0;
		r->second_round_cigar[i].lost_base_length = r->cigars[i].lost_base_length = 0;
		r->second_round_cigar[i].lost_base = r->cigars[i].lost_base = NULL;
	}
	///r->pb_regions = NULL;

	REALLOC(r->paf, r->total_reads);
	memset(r->paf+total_reads0, 0, sizeof((*(r->paf)))*(r->total_reads - total_reads0));
	long long n_read; f_flag = fread(&n_read, sizeof(n_read), 1, fpo); ma_hit_t t;
	for (i = total_reads0; i < r->total_reads; i++) {
		f_flag += fread(&(r->paf[i].is_fully_corrected), sizeof(r->paf[i].is_fully_corrected), 1, fpo);
        f_flag += fread(&(r->paf[i].is_abnormal), sizeof(r->paf[i].is_abnormal), 1, fpo);
        f_flag += fread(&(r->paf[i].length), sizeof(r->paf[i].length), 1, fpo);

        if(r->paf[i].length == 0) continue;
        for (k = 0; k < r->paf[i].length; k++) read_ma(&t, fpo);
		r->paf[i].length = 0;
	}

	REALLOC(r->reverse_paf, r->total_reads); 
	memset(r->reverse_paf+total_reads0, 0, sizeof((*(r->reverse_paf)))*(r->total_reads - total_reads0));

    fclose(fp); fclose(fpo);
    fprintf(stderr, "Reads has been loaded.\n");
	return 1;
}


int destory_read_bin(All_reads* r)
{

	uint64_t i = 0;
	for (i = 0; i < r->total_reads; i++)
	{
		if (r->N_site[i]) free(r->N_site[i]);
		if (r->read_sperate[i]) free(r->read_sperate[i]);
		if (r->cigars[i].record) free(r->cigars[i].record);
		if (r->cigars[i].lost_base) free(r->cigars[i].lost_base);
		if (r->second_round_cigar[i].record) free(r->second_round_cigar[i].record);
		if (r->second_round_cigar[i].lost_base) free(r->second_round_cigar[i].lost_base);
	}

	free(r->N_site);
	free(r->read_length);
	free(r->read_size);
	free(r->read_sperate);
	free(r->name);
	free(r->name_index);
	free(r->trio_flag);
	free(r->cigars);
	free(r->second_round_cigar);
	return 1;
}



void ha_insert_read_len(All_reads *r, int read_len, int name_len)
{
	r->total_reads++;
	r->total_reads_bases += (uint64_t)read_len;
	r->total_name_length += (uint64_t)name_len;

	// must +1
	if (r->index_size < r->total_reads + 2) {
		r->index_size = r->index_size * 2 + 2;
		r->read_length = (uint64_t*)realloc(r->read_length, sizeof(uint64_t) * r->index_size);
		r->name_index_size = r->name_index_size * 2 + 2;
		r->name_index = (uint64_t*)realloc(r->name_index, sizeof(uint64_t) * r->name_index_size);
	}

	r->read_length[r->total_reads - 1] = read_len;
	r->name_index[r->total_reads] = r->name_index[r->total_reads - 1] + name_len;
}

void malloc_All_reads(All_reads* r)
{
	r->read_size = (uint64_t*)malloc(sizeof(uint64_t)*r->total_reads);
	memcpy(r->read_size, r->read_length, sizeof(uint64_t)*r->total_reads);

	r->read_sperate = (uint8_t**)malloc(sizeof(uint8_t*)*r->total_reads);
	long long i = 0;
	for (i = 0; i < (long long)r->total_reads; i++)
	{
		r->read_sperate[i] = (uint8_t*)malloc(sizeof(uint8_t)*(r->read_length[i]/4+1));
	}

	r->cigars = (Compressed_Cigar_record*)malloc(sizeof(Compressed_Cigar_record)*r->total_reads);
	r->second_round_cigar = (Compressed_Cigar_record*)malloc(sizeof(Compressed_Cigar_record)*r->total_reads);
	r->paf = (ma_hit_t_alloc*)malloc(sizeof(ma_hit_t_alloc)*r->total_reads);
	r->reverse_paf = (ma_hit_t_alloc*)malloc(sizeof(ma_hit_t_alloc)*r->total_reads);
	///r->pb_regions = (kvec_t_u64_warp*)malloc(r->total_reads*sizeof(kvec_t_u64_warp));

	for (i = 0; i < (long long)r->total_reads; i++)
	{
		r->second_round_cigar[i].size = r->cigars[i].size = 0;
		r->second_round_cigar[i].length = r->cigars[i].length = 0;
		r->second_round_cigar[i].record = r->cigars[i].record = NULL;

		r->second_round_cigar[i].lost_base_size = r->cigars[i].lost_base_size = 0;
		r->second_round_cigar[i].lost_base_length = r->cigars[i].lost_base_length = 0;
		r->second_round_cigar[i].lost_base = r->cigars[i].lost_base = NULL;
		init_ma_hit_t_alloc(&(r->paf[i]));
		init_ma_hit_t_alloc(&(r->reverse_paf[i]));
		///kv_init(r->pb_regions[i].a);
	}

	r->name = (char*)malloc(sizeof(char)*r->total_name_length);
	r->N_site = (uint64_t**)calloc(r->total_reads, sizeof(uint64_t*));
	r->trio_flag = (uint8_t*)malloc(r->total_reads*sizeof(uint8_t));
	memset(r->trio_flag, AMBIGU, r->total_reads*sizeof(uint8_t));
}

void destory_UC_Read(UC_Read* r)
{
	free(r->seq);
}

void init_aux_table()
{
	if (bit_t_seq_table[0][0] == 0)
	{
		uint64_t i = 0;

		for (i = 0; i < 256; i++)
		{
			bit_t_seq_table[i][0] = s_H[((i >> 6)&(uint64_t)3)];
			bit_t_seq_table[i][1] = s_H[((i >> 4)&(uint64_t)3)];
			bit_t_seq_table[i][2] = s_H[((i >> 2)&(uint64_t)3)];
			bit_t_seq_table[i][3] = s_H[(i&(uint64_t)3)];

			bit_t_seq_table_rc[i][0] = RC_CHAR(bit_t_seq_table[i][3]);
			bit_t_seq_table_rc[i][1] = RC_CHAR(bit_t_seq_table[i][2]);
			bit_t_seq_table_rc[i][2] = RC_CHAR(bit_t_seq_table[i][1]);
			bit_t_seq_table_rc[i][3] = RC_CHAR(bit_t_seq_table[i][0]);
		}
	}
}

void init_UC_Read(UC_Read* r)
{
	r->length = 0;
	r->size = 0;
	r->seq = NULL;
	if (bit_t_seq_table[0][0] == 0)
	{
		uint64_t i = 0;

		for (i = 0; i < 256; i++)
		{
			bit_t_seq_table[i][0] = s_H[((i >> 6)&(uint64_t)3)];
			bit_t_seq_table[i][1] = s_H[((i >> 4)&(uint64_t)3)];
			bit_t_seq_table[i][2] = s_H[((i >> 2)&(uint64_t)3)];
			bit_t_seq_table[i][3] = s_H[(i&(uint64_t)3)];

			bit_t_seq_table_rc[i][0] = RC_CHAR(bit_t_seq_table[i][3]);
			bit_t_seq_table_rc[i][1] = RC_CHAR(bit_t_seq_table[i][2]);
			bit_t_seq_table_rc[i][2] = RC_CHAR(bit_t_seq_table[i][1]);
			bit_t_seq_table_rc[i][3] = RC_CHAR(bit_t_seq_table[i][0]);
		}
	}
}


void recover_UC_Read_sub_region(char* r, int64_t start_pos, int64_t length, uint8_t strand, All_reads* R_INF, int64_t ID)
{
	int64_t readLen = Get_READ_LENGTH((*R_INF), ID);
	int64_t end_pos = start_pos + length - 1, begLen, tailLen, offset, mn, src_i, des_i, i;
	uint8_t* src = Get_READ((*R_INF), ID);

	if(strand == 0) {
		offset = start_pos&3;
		begLen = 4-offset;
		if(begLen > length) begLen = length;
		tailLen = (length-begLen)&3;
		mn = (length - begLen - tailLen)>>2;
		src_i = start_pos; des_i = 0; i = 0;

		if(begLen > 0) {
			memcpy(r+des_i, bit_t_seq_table[src[src_i>>2]]+offset, begLen);
			des_i += begLen; src_i += begLen;
		}

		for (i = 0; i < mn; i++) {
			memcpy(r+des_i, bit_t_seq_table[src[src_i>>2]], 4);
			des_i += 4; src_i += 4;
		}

		if(tailLen > 0) {
			memcpy(r+des_i, bit_t_seq_table[src[src_i>>2]], tailLen);
			des_i += tailLen; src_i += tailLen;
		}

		if (R_INF->N_site[ID]) {
			for (i = 1; i <= (long long)R_INF->N_site[ID][0]; i++) {
				if ((long long)R_INF->N_site[ID][i] >= start_pos && (long long)R_INF->N_site[ID][i] <= end_pos) {
					r[R_INF->N_site[ID][i] - start_pos] = 'N';
				}
				else if((long long)R_INF->N_site[ID][i] > end_pos) {
					break;
				}
			}
		}
	}
	else {
		start_pos = readLen - start_pos - 1;
		end_pos = readLen - end_pos - 1;

		begLen = (start_pos+1)&3;
		offset = 4 - begLen;
		if(begLen > length) begLen = length;
		tailLen = (length-begLen)&3;
		mn = (length - begLen - tailLen)>>2;
		src_i = start_pos; des_i = 0; i = 0;

		if(begLen > 0) {
			memcpy(r+des_i, bit_t_seq_table_rc[src[src_i>>2]]+offset, begLen);
			des_i += begLen; src_i -= begLen;
		}

		for (i = 0; i < mn; i++) {
			memcpy(r+des_i, bit_t_seq_table_rc[src[src_i>>2]], 4);
			des_i += 4; src_i -= 4;
		}

		if(tailLen > 0) {
			memcpy(r+des_i, bit_t_seq_table_rc[src[src_i>>2]], tailLen);
			des_i += tailLen; src_i -= tailLen;
		}

		/**
		if (R_INF->N_site[ID]) {
			offset = readLen - start_pos - 1;
			for (i = 1; i <= (long long)R_INF->N_site[ID][0]; i++) {
				if ((long long)R_INF->N_site[ID][i] >= end_pos && (long long)R_INF->N_site[ID][i] <= start_pos) {
					r[readLen - R_INF->N_site[ID][i] - 1 - offset] = 'N';
				}
				else if((long long)R_INF->N_site[ID][i] > start_pos) {
					break;
				}
			}
		}
		**/
		if (R_INF->N_site[ID]) {
			start_pos = readLen - start_pos - 1; end_pos = readLen - end_pos - 1;
			for (i = 1; i <= (long long)R_INF->N_site[ID][0]; i++) {
				offset = readLen - R_INF->N_site[ID][i] - 1;
				if(offset >= start_pos && offset <= end_pos) {
					r[offset - start_pos] = 'N';
				} else if(offset < start_pos) {
					break;
				}
			}
		}
	}
}


void recover_UC_sub_Read(UC_Read* i_r, long long start_pos, long long length, uint8_t strand, All_reads* R_INF, long long ID)
{
	i_r->length = length;i_r->RID = ID;
	if (i_r->length + 8 > i_r->size)
	{
		i_r->size = i_r->length + 4;
		i_r->seq = (char*)realloc(i_r->seq,sizeof(char)*(i_r->size));
	}
	char* r = i_r->seq;
	long long readLen = Get_READ_LENGTH((*R_INF), ID);
	uint8_t* src = Get_READ((*R_INF), ID);

	long long i;
	long long copyLen;
	long long end_pos = start_pos + length - 1;

	if (strand == 0)
	{
		i = start_pos;
		copyLen = 0;

		long long initLen = start_pos % 4;

		if (initLen != 0)
		{
			memcpy(r, bit_t_seq_table[src[i>>2]] + initLen, 4 - initLen);
			copyLen = copyLen + 4 - initLen;
			i = i + copyLen;
		}

		while (copyLen < length)
		{
			memcpy(r+copyLen, bit_t_seq_table[src[i>>2]], 4);
			copyLen = copyLen + 4;
			i = i + 4;
		}

		if (R_INF->N_site[ID])
		{
			for (i = 1; i <= (long long)R_INF->N_site[ID][0]; i++)
			{
				if ((long long)R_INF->N_site[ID][i] >= start_pos && (long long)R_INF->N_site[ID][i] <= end_pos)
				{
					r[R_INF->N_site[ID][i] - start_pos] = 'N';
				}
				else if((long long)R_INF->N_site[ID][i] > end_pos)
				{
					break;
				}
			}
		}
	}
	else
	{
		start_pos = readLen - start_pos - 1;
		end_pos = readLen - end_pos - 1;

		///start_pos > end_pos
		i = start_pos;
		copyLen = 0;
		long long initLen = (start_pos + 1) % 4;

		if (initLen != 0)
		{
			memcpy(r, bit_t_seq_table_rc[src[i>>2]] + 4 - initLen, initLen);
			copyLen = copyLen + initLen;
			i = i - initLen;
		}

		while (copyLen < length)
		{
			memcpy(r+copyLen, bit_t_seq_table_rc[src[i>>2]], 4);
			copyLen = copyLen + 4;
			i = i - 4;
		}

		if (R_INF->N_site[ID])
		{
			long long offset = readLen - start_pos - 1;

			for (i = 1; i <= (long long)R_INF->N_site[ID][0]; i++)
			{
				if ((long long)R_INF->N_site[ID][i] >= end_pos && (long long)R_INF->N_site[ID][i] <= start_pos)
				{
					r[readLen - R_INF->N_site[ID][i] - 1 - offset] = 'N';
				}
				else if((long long)R_INF->N_site[ID][i] > start_pos)
				{
					break;
				}
			}
		}
	}
}




void recover_UC_Read(UC_Read* r, const All_reads *R_INF, uint64_t ID)
{
	r->length = Get_READ_LENGTH((*R_INF), ID);
	uint8_t* src = Get_READ((*R_INF), ID);

	if (r->length + 4 > r->size) {
		r->size = r->length + 4;
		r->seq = (char*)realloc(r->seq,sizeof(char)*(r->size));
	}

	uint64_t i = 0;

	if(src) {
		while ((long long)i < r->length) {
			memcpy(r->seq+i, bit_t_seq_table[src[i>>2]], 4);
			i = i + 4;
		}

		if (R_INF->N_site[ID]) {
			for (i = 1; i <= R_INF->N_site[ID][0]; i++) r->seq[R_INF->N_site[ID][i]] = 'N';
		}
	} else {///N
		memset(r->seq, 'N', r->length);
	}

	r->RID = ID;	
}

void recover_UC_Read_RC(UC_Read* r, All_reads* R_INF, uint64_t ID)
{
	r->length = Get_READ_LENGTH((*R_INF), ID);
	uint8_t* src = Get_READ((*R_INF), ID);

	if (r->length + 4 > r->size) {
		r->size = r->length + 4;
		r->seq = (char*)realloc(r->seq,sizeof(char)*(r->size));
	}

	long long last_chr = r->length % 4;
	long long i = r->length / 4 - 1 + (last_chr != 0);
	long long index = 0;

	if(src) {
		if(last_chr!=0) {
			memcpy(r->seq + index, bit_t_seq_table_rc[src[i]] + 4 - last_chr, last_chr);
			index = last_chr;
			i--;
		}

		while (i >= 0) {
			memcpy(r->seq + index, bit_t_seq_table_rc[src[i]], 4);
			i--;
			index = index + 4;
		}

		if (R_INF->N_site[ID]) {
			for (i = 1; i <= (long long)R_INF->N_site[ID][0]; i++) {
				r->seq[r->length - R_INF->N_site[ID][i] - 1] = 'N';
			}
		}
	} else {///N
		memset(r->seq, 'N', r->length);
	}

	
}

#define COMPRESS_BASE {c = seq_nt6_table[(uint8_t)src[i]];\
		if (c >= 4)\
		{\
			c = 0;\
			(*N_site_lis)[N_site_i] = i;\
			N_site_i++;\
		}\
		i++;}\

void ha_compress_base(uint8_t* dest, char* src, uint64_t src_l, uint64_t** N_site_lis, uint64_t N_site_occ)
{
	///N_site_lis saves the pos of all Ns in this read
	///N_site_lis[0] is the number of Ns
	free((*N_site_lis));
	if (N_site_occ)
	{
		(*N_site_lis) = (uint64_t*)malloc(sizeof(uint64_t)*(N_site_occ + 1));
		(*N_site_lis)[0] = N_site_occ;
	}
	else
	{
		(*N_site_lis) = NULL;
	}

	uint64_t i = 0;
	uint64_t N_site_i = 1;
	uint64_t dest_i = 0;
	uint8_t tmp = 0;
	uint8_t c = 0;

	while (i + 4 <= src_l)
	{
		tmp = 0;

		COMPRESS_BASE;
		tmp = tmp | (c<<6);

		COMPRESS_BASE;
		tmp = tmp | (c<<4);

		COMPRESS_BASE;
		tmp = tmp | (c<<2);

		COMPRESS_BASE;
		tmp = tmp | c;

		dest[dest_i] = tmp;

		dest_i++;
	}

	//at most 3 bases here
	uint64_t shift = 6;
	if (i < src_l)
	{
		tmp = 0;

		while (i < src_l)
		{
			COMPRESS_BASE;
			tmp = tmp | (c << shift);
			shift = shift -2;
		}
		
		dest[dest_i] = tmp;
		dest_i++;
	}
}

void reverse_complement(char* pattern, uint64_t length)
{
	uint64_t i = 0;
	uint64_t end = length / 2;
	char k;
	uint64_t index;

	for (i = 0; i < end; i++)
	{
		
		index = length - i - 1;
		k = pattern[index];
		pattern[index] = RC_CHAR(pattern[i]);
		pattern[i] = RC_CHAR(k);
	}

	if(length&(uint64_t)1)
	{
		pattern[end] = RC_CHAR(pattern[end]);
	}
}


void init_Debug_reads(Debug_reads* x, const char* file)
{
	int nameLen, i, bufLen = 1000;
	if((uint64_t)(bufLen) < strlen(file) + 50) bufLen = strlen(file) + 50;
	char* Name_Buffer = (char*)malloc(sizeof(char)*bufLen);
	fprintf(stderr, "Queried debugging reads at: %s\n", file);

	x->fp = fopen(file,"r");
	x->query_num = 0;
	
	while(fgets(Name_Buffer, bufLen, x->fp))
	{
		x->query_num++;
	}
	x->read_name = (char**)malloc(sizeof(char*)*x->query_num);
	x->candidate_count = (kvec_t_u64_warp*)malloc(sizeof(kvec_t_u64_warp)*x->query_num);
	x->read_id = NULL; MALLOC(x->read_id, x->query_num); 
	memset(x->read_id, -1, sizeof((*(x->read_id)))*x->query_num);
	fseek(x->fp, 0, SEEK_SET);

	i = 0;
	while(fgets(Name_Buffer, bufLen, x->fp))
	{
		nameLen = strlen(Name_Buffer) - 1;
		x->read_name[i] = (char*)malloc(sizeof(char)*(nameLen+1));
		memcpy(x->read_name[i], Name_Buffer, sizeof(char)*nameLen);
		x->read_name[i][nameLen] = '\0';
		kv_init(x->candidate_count[i].a);
		i++;
	}
	
    fclose(x->fp);

	sprintf(Name_Buffer, "%s.debug.stdout", file);
	x->fp = fopen(Name_Buffer,"w");
	fprintf(stderr, "Print debugging information to: %s\n", Name_Buffer);

	sprintf(Name_Buffer, "%s.debug.r0.fa", file);
	x->fp_r0 = fopen(Name_Buffer,"w");
	fprintf(stderr, "Print raw reads to: %s\n", Name_Buffer);

	sprintf(Name_Buffer, "%s.debug.r1.fa", file);
	x->fp_r1 = fopen(Name_Buffer,"w");
	fprintf(stderr, "Print corrected reads to: %s\n", Name_Buffer);

	free(Name_Buffer);
}

void destory_Debug_reads(Debug_reads* x)
{
	uint64_t i;
	for (i = 0; i < x->query_num; i++)
	{
		free(x->read_name[i]);
		kv_destroy(x->candidate_count[i].a);
	}
	
	free(x->read_name); free(x->read_id);
	fclose(x->fp); fclose(x->fp_r0); fclose(x->fp_r1);
}


void init_all_ul_t(all_ul_t *x, All_reads *hR) {
	memset(x, 0, sizeof(*x));
	x->hR = hR;
	init_aux_table();
}
void destory_all_ul_t(all_ul_t *x) {
	uint64_t i;
	for (i = 0; i < x->n; i++) {
		free(x->a[i].N_site.a); free(x->a[i].r_base.a); free(x->a[i].bb.a);
	}
	free(x->a);

	for (i = 0; i < x->nid.n; i++) free(x->nid.a[i].a);
	free(x->nid.a);
	free(x->ridx.idx.a); free(x->ridx.occ.a);

	// if(x->mm) {
	// 	for (i = 0; i < x->mm->n; i++) free(x->mm->a[i].a);
	// 	free(x->mm->a); free(x->mm); x->mm = NULL;
	// }
}

void ha_encode_base(uint8_t* dest, char* src, uint64_t src_l, N_t *nn, uint64_t nn_offset)
{
    uint64_t i = 0;
    uint64_t dest_i = 0;
    uint8_t tmp = 0;
    uint8_t c = 0;

    while (i + 4 <= src_l) {
        tmp = 0;
		// fprintf(stderr, "i->%lu, src_l->%ld, src[i]->%c, (uint8_t)src[i]->%u\n", i, src_l, src[i], (uint8_t)src[i]);
        c = seq_nt6_table[(uint8_t)src[i]];
		if (c >= 4) {
			c = 0; kv_push(uint32_t, *nn, i+nn_offset);
		}
		i++;
        tmp = tmp | (c<<6);

        c = seq_nt6_table[(uint8_t)src[i]];
		if (c >= 4) {
			c = 0; kv_push(uint32_t, *nn, i+nn_offset);
		}
		i++;
        tmp = tmp | (c<<4);

        c = seq_nt6_table[(uint8_t)src[i]];
		if (c >= 4) {
			c = 0; kv_push(uint32_t, *nn, i+nn_offset);
		}
		i++;
        tmp = tmp | (c<<2);

        c = seq_nt6_table[(uint8_t)src[i]];
		if (c >= 4) {
			c = 0; kv_push(uint32_t, *nn, i+nn_offset);
		}
		i++;
        tmp = tmp | c;

        dest[dest_i] = tmp;

        dest_i++;
    }

    //at most 3 bases here
    uint64_t shift = 6;
    if (i < src_l) {
        tmp = 0;
        while (i < src_l) {
            c = seq_nt6_table[(uint8_t)src[i]];
			if (c >= 4) {
				c = 0; kv_push(uint32_t, *nn, i+nn_offset);
			}
			i++;
            tmp = tmp | (c << shift);
            shift = shift -2;
        }
        
        dest[dest_i] = tmp;
        dest_i++;
    }
}

#define B4L(x) (((x)>>2)+(((x)&3)?1:0))
void push_subblock_original_bases(char* str, all_ul_t *x, ul_vec_t *p, uint32_t s, uint32_t e, uint32_t subLen)//for debug
{
	uc_block_t *b = NULL;
	uint32_t qs = s, qe = e;
	while (qs < e) {
		qe = qs + subLen; if(qe > e) qe = e;
		kv_pushp(uc_block_t, p->bb, &b);
		b->hid = (uint32_t)-1; b->rev = 0;
		b->qs = qs; b->qe = qe;
		b->ts = p->r_base.n; b->te = b->ts + B4L(b->qe - b->qs);
		kv_resize(uint8_t, p->r_base, b->te); p->r_base.n = b->te;
		ha_encode_base(p->r_base.a+b->ts, str+b->qs, b->qe-b->qs, &(p->N_site), b->qs);

		qs = qe;
	}	
}

void append_ul_t_compress_ovlp(all_ul_t *x, uint64_t *rid, char* id, int64_t id_l, char* str, int64_t str_l, ul_ov_t *o, int64_t on) {
	int64_t i, mine, maxs, ovlp, end;
	ul_vec_t *p = NULL;
	nid_t *np = NULL;
	ul_ov_t *z = NULL, *zp = NULL;
	uc_block_t *b = NULL;

	if(id) {
		kv_pushp(nid_t, x->nid, &np);
		np->n = id_l; MALLOC(np->a, np->n+1); memcpy(np->a, id, id_l); np->a[id_l] = '\0';
	}

	if(str) {
		if(rid == NULL) {
			kv_pushp(ul_vec_t, *x, &p);
			memset(p, 0, sizeof(*p));
		} else {
			if((*rid) >= x->m) kv_resize(ul_vec_t, *x, (*rid) + 1);
			if((*rid) >= x->n) {
				memset(x->a+x->n, 0, sizeof(*p)*((*rid) + 1 - x->n));
				x->n = (*rid) + 1;
			}
			p = &(x->a[(*rid)]);
		}
		

		p->bb.n = p->N_site.n = p->r_base.n = 0;
		p->rlen = str_l; 
		
		if(o == NULL || on == 0) on = 0;
		for (i = end = 0, zp = NULL; i < on; i++) {
			z = &(o[i]);
			if(!z->el) continue;
			if(zp) {
				mine = MIN(zp->qe, z->qe); maxs = MAX(zp->qs, z->qs);
				ovlp = mine - maxs; 
				if(zp->qe >= z->qe && zp->qs <= z->qs) continue;
			} else {
				ovlp = -z->qs; 
			}
			if(ovlp < 0) {///push original bases
				kv_pushp(uc_block_t, p->bb, &b);
				b->hid = (uint32_t)-1/**x->mm**/; b->rev = 0; b->base = 1; b->pchain = 0;
				b->qs = end; b->qe = b->qs - ovlp;
				b->ts = p->r_base.n; b->te = b->ts + B4L(-ovlp);
				kv_resize(uint8_t, p->r_base, b->te); p->r_base.n = b->te;
				ha_encode_base(p->r_base.a+b->ts, str+b->qs, b->qe-b->qs, &(p->N_site), b->qs);
			} 

			///push ovlp bases
			kv_pushp(uc_block_t, p->bb, &b);
			b->hid = (z->tn<<1)>>1; b->rev = z->rev; b->base = 0; b->pchain = 0;
			b->qs = z->qs + (ovlp>0?ovlp:0); b->qe = z->qe;
			if(z->rev) {
				b->ts = z->ts; b->te = z->ts + (b->qe - b->qs);
			} else {
				b->ts = z->te - (b->qe - b->qs); b->te = z->te;
			}

			end = MAX(zp?zp->qe:0, z->qe);
		}
		
		if(end < str_l) {///push original bases
			kv_pushp(uc_block_t, p->bb, &b);
			b->hid = (uint32_t)-1/**x->mm**/; b->rev = 0; b->base = 1; b->pchain = 0;
			b->qs = end; b->qe = str_l;
			b->ts = p->r_base.n; b->te = b->ts + B4L(b->qe - b->qs);
			kv_resize(uint8_t, p->r_base, b->te); p->r_base.n = b->te;
			ha_encode_base(p->r_base.a+b->ts, str+b->qs, b->qe-b->qs, &(p->N_site), b->qs);		
			// push_subblock_original_bases(str, x, p, end, str_l, 321);//for debug
		}
		// char *sst = NULL; CALLOC(sst, str_l);//for debug
		// retrieve_ul_t(NULL, sst, x, rid?*rid:x->n-1, 0, 0, -1);
		// if(memcmp(sst, str, str_l)) {
		// 	fprintf(stderr, "ap-Wrong read, id: %ld, [%d, %ld)\n", (int64_t)(rid?*rid:x->n-1), 0, str_l);
		// 	for (i = 0; i < str_l; i++) {
		// 		if(sst[i] != str[i]) fprintf(stderr, "[%ld] input:%c, decompress:%c\n", i, str[i], sst[i]);
		// 	}
		// }
		// free(sst);
	}
}

void debug_append_ul_t(ul_ov_t *o, int64_t on, ul_vec_t *p)
{
	int64_t k, l = 0;
    uint32_t sp = (uint32_t)-1, ep = (uint32_t)-1, qss, qee, m;
    for (k = on-1; k >= 0; k--) {
        if(sp == (uint32_t)-1 || o[k].qe <= sp) {
            if(sp != (uint32_t)-1) l += ep - sp;

			if(sp == (uint32_t)-1) qss = o[k].qe, qee = p->rlen;
			else qss = o[k].qe, qee = sp;
			if(qee > qss) {
				for (m = 0; m < p->bb.n; m++) {
					if(qss == (p->bb.a[m].qs+((p->bb.a[m].hid>>15)&(0x7fffU))) && 
							qee == (p->bb.a[m].qe-(p->bb.a[m].hid&(0x7fffU)))) {
						break;
					}
				}
				if(m >= p->bb.n) fprintf(stderr, "ERROR\n");
			}

            sp = o[k].qs;
            ep = o[k].qe;
        } else {
            sp = MIN(sp, o[k].qs);
        }
    }
    if(sp != (uint32_t)-1) l += ep - sp;

	if(sp == (uint32_t)-1) qss = 0, qee = p->rlen;
	else qss = 0, qee = sp;
	if(qee > qss) {
		for (m = 0; m < p->bb.n; m++) {
			if(qss == (p->bb.a[m].qs+((p->bb.a[m].hid>>15)&(0x7fffU))) && 
					qee == (p->bb.a[m].qe-(p->bb.a[m].hid&(0x7fffU)))) {
				break;
			}
		}
		if(m >= p->bb.n) fprintf(stderr, "ERROR\n");
	}
}

void append_ul_t_back(all_ul_t *x, uint64_t *rid, char* id, int64_t id_l, char* str, int64_t str_l, ul_ov_t *o, int64_t on, float p_chain_rate) {
	int64_t i, mine, maxs, ovlp, st, et, bl = 0, pc = 0;
	uint32_t o_l, o_r;
	ul_vec_t *p = NULL;
	nid_t *np = NULL;
	ul_ov_t *z = NULL;
	uc_block_t *b = NULL, tt;

	if(id) {
		kv_pushp(nid_t, x->nid, &np);
		np->n = id_l; MALLOC(np->a, np->n+1); memcpy(np->a, id, id_l); np->a[id_l] = '\0';
	}

	if(str||str_l) {
		if(rid == NULL) {
			kv_pushp(ul_vec_t, *x, &p);
			memset(p, 0, sizeof(*p));
		} else {
			if((*rid) >= x->m) kv_resize(ul_vec_t, *x, (*rid) + 1);
			if((*rid) >= x->n) {
				memset(x->a+x->n, 0, sizeof(*p)*((*rid) + 1 - x->n));
				x->n = (*rid) + 1;
			}
			p = &(x->a[(*rid)]);
		}
		// if((*rid) == 23) fprintf(stderr, "#rid->%lu, on->%ld\n", *rid, on);

		p->bb.n = p->N_site.n = p->r_base.n = 0; p->dd = 0;
		p->rlen = str_l; 
		
		if(o == NULL || on == 0) on = 0;
		for (i = on-1, st = et = str_l; i >= 0; i--) {
			z = &(o[i]);
			if(z->el) {
				mine = MIN(et, ((int64_t)z->qe)); maxs = MAX(st, ((int64_t)z->qs));
				ovlp = mine - maxs; 
				
				if(ovlp < 0) {///push original bases
					kv_pushp(uc_block_t, p->bb, &b);
					b->hid = 0/**x->mm**/; b->rev = 0; b->base = 1; b->pchain = 0; b->el = 0;
					b->qe = maxs; b->qs = b->qe + ovlp; bl += (b->qe-b->qs);
					o_l = (b->qs >= UL_FLANK?UL_FLANK:b->qs);
					o_r = ((str_l-b->qe)>=UL_FLANK?UL_FLANK:(str_l-b->qe));
					b->hid |= (o_l<<15); b->hid |= o_r;
					b->qs -= o_l; b->qe += o_r;
					b->ts = p->r_base.n; b->te = b->ts + B4L(b->qe-b->qs);
					kv_resize(uint8_t, p->r_base, b->te); p->r_base.n = b->te;
					// if(!str) fprintf(stderr, "+rid->%lu\n", *rid);
					ha_encode_base(p->r_base.a+b->ts, str+b->qs, b->qe-b->qs, &(p->N_site), b->qs);
				} 

				st = MIN(st, z->qs);
			}

			///push ovlp bases
			kv_pushp(uc_block_t, p->bb, &b);
			b->hid = (z->tn<<1)>>1; b->rev = z->rev; b->base = 0; b->el = z->el;
			b->pchain = ((z->tn&((uint32_t)(0x80000000)))?1:0);
			b->qs = z->qs; b->qe = z->qe;
			b->ts = z->ts; b->te = z->te;
			if(b->pchain) pc++;
			
			// st = MIN(st, z->qs);
		}
		
		if(st > 0) {///push original bases
			kv_pushp(uc_block_t, p->bb, &b);
			b->hid = 0/**x->mm**/; b->rev = 0; b->base = 1; b->pchain = 0; b->el = 0;
			b->qe = st; b->qs = 0; bl += (b->qe-b->qs);
			o_l = (b->qs >= UL_FLANK?UL_FLANK:b->qs);
			o_r = ((str_l-b->qe)>=UL_FLANK?UL_FLANK:(str_l-b->qe));
			b->hid |= (o_l<<15); b->hid |= o_r;
			b->qs -= o_l; b->qe += o_r;
			b->ts = p->r_base.n; b->te = b->ts + B4L(b->qe-b->qs);
			kv_resize(uint8_t, p->r_base, b->te); p->r_base.n = b->te;
			// if(!str) fprintf(stderr, "-rid->%lu, st->%ld, str_l->%ld\n", *rid, st, str_l);
			ha_encode_base(p->r_base.a+b->ts, str+b->qs, b->qe-b->qs, &(p->N_site), b->qs);	
			// push_subblock_original_bases(str, x, p, end, str_l, 321);//for debug
		}

		if(pc > 0) p->dd = 3;
		if((pc == on) && ((str_l-bl) > (str_l*p_chain_rate))) p->dd = 2;
		if((pc == on) && (bl == 0)) p->dd = 1;
		// debug_append_ul_t(o, on, p);
		// char *sst = NULL; CALLOC(sst, str_l);//for debug
		// retrieve_ul_t(NULL, sst, x, rid?*rid:x->n-1, 0, 0, -1);
		// if(memcmp(sst, str, str_l)) {
		// 	fprintf(stderr, "ap-Wrong read, id: %ld, [%d, %ld)\n", (int64_t)(rid?*rid:x->n-1), 0, str_l);
		// 	for (i = 0; i < str_l; i++) {
		// 		if(sst[i] != str[i]) fprintf(stderr, "[%ld] input:%c, decompress:%c\n", i, str[i], sst[i]);
		// 	}
		// }
		// free(sst);
		ovlp = p->bb.n>>1;
		for (i = 0; i < ovlp; i++) {
			tt = p->bb.a[i]; 
			p->bb.a[i] = p->bb.a[p->bb.n-i-1];
			p->bb.a[p->bb.n-i-1] = tt;
		}
		
	}
}


void determine_chain_distance(ul_ov_t *o, int64_t on, ul_vec_t *p, ma_hit_t_alloc *src, int64_t max_hang, int64_t min_ovlp, int64_t rid)
{
	int64_t k, i, m, l = 0, r, last_i, last_dis; ul_ov_t *z = NULL; 
	uint32_t li_v, lj_v, t, qn, tn; ma_hit_t_alloc *x = NULL; asg_arc_t te;

	for (k = on-1; k >= 0; --k) {
		z = &(o[k]);
		if((!z->el) || (z->sec == SEC_MODE)) continue;
		///if z->el = 1, z->qn must work
		if(p->bb.a[z->qn].pidx != (uint32_t)-1) continue;
		last_i = -1; last_dis = 0;
		for (i = k, l = 0; i >= 0;) {
			assert((o[i].tn&((uint32_t)(0x80000000))));
			if(o[i].el) {
				last_i = o[i].qn; last_dis = l;
			}
			m = i; i = o[i].sec; if(o[m].sec == SEC_MODE) i = -1;
			// i = ((o[i].sec == SEC_MODE)?-1:o[i].sec); 
			if(i < 0) break;
			// if(i >= on || i < 0) fprintf(stderr, "m->%ld, i->%ld, rid->%ld, on->%ld, o[m].sec->%u\n", m, i, rid, on, o[m].sec);
			li_v = (o[m].tn<<1)|o[m].rev; li_v ^= 1;
			lj_v = (o[i].tn<<1)|o[i].rev; lj_v ^= 1;

			x = &(src[li_v>>1]);
			for (t = 0; t < x->length; t++) {
				qn = Get_qn(x->buffer[t]);
        		tn = Get_tn(x->buffer[t]);
				if(qn == (li_v>>1) && tn == (lj_v>>1)) {
					r = ma_hit2arc(&(x->buffer[t]), Get_READ_LENGTH(R_INF, qn), Get_READ_LENGTH(R_INF, tn),
						max_hang, asm_opt.max_hang_rate, min_ovlp, &te);
            		if(r < 0) continue;
            		if((te.ul>>32) != li_v || te.v != lj_v) continue;
					l += (uint32_t)te.ul;
					break;
				}
			}
			// if(t>=x->length) {
			// 	fprintf(stderr, "m->%ld(id::%u), i->%ld(id::%u), rid->%ld, on->%ld\n", 
			// 	m, (o[m].tn<<1)>>1, i, (o[i].tn<<1)>>1, rid, on);
			// 	// fprintf(stderr, "[i::%ld]id::%u\t\tq::[%u, %u)\tt::[%u, %u)\n", 
			// 	// 	i, o[i].tn, o[i].qs, o[i].qe, o[i].ts, o[i].te);
			// 	// fprintf(stderr, "[m::%ld]id::%u\t\tq::[%u, %u)\tt::[%u, %u)\n", 
			// 	// 	m, o[m].tn, o[m].qs, o[m].qe, o[m].ts, o[m].te);
			// 	fprintf(stderr, "[i::%ld]id::%u\t%.*s\t%c\tq::[%u, %u)\tt::[%u, %u)\n", 
			// 		i, (o[i].tn<<1)>>1, (int)Get_NAME_LENGTH(R_INF, ((o[i].tn<<1)>>1)), Get_NAME(R_INF, ((o[i].tn<<1)>>1)), 
			// 		"+-"[o[i].rev], o[i].qs, o[i].qe, o[i].ts, o[i].te);
			// 	fprintf(stderr, "[m::%ld]id::%u\t%.*s\t%c\tq::[%u, %u)\tt::[%u, %u)\n", 
			// 		m, (o[m].tn<<1)>>1, (int)Get_NAME_LENGTH(R_INF, ((o[m].tn<<1)>>1)), Get_NAME(R_INF, ((o[m].tn<<1)>>1)), 
			// 		"+-"[o[m].rev], o[m].qs, o[m].qe, o[m].ts, o[m].te);
			// 	exit(1);
			// }
			assert(t<x->length);
			if(!(o[i].el)) continue;
			assert(last_i>=0);
			p->bb.a[last_i].pidx = o[i].qn;
			p->bb.a[last_i].pdis = l - last_dis;///TODO: enable pdis
			p->bb.a[p->bb.a[last_i].pidx].aidx = last_i;
		}
	}
}

void append_ul_t(all_ul_t *x, uint64_t *rid, char* id, int64_t id_l, char* str, int64_t str_l, ul_ov_t *o, int64_t on, float p_chain_rate, const ug_opt_t *uopt, uint32_t save_bases) {
	int64_t i, mine, maxs, ovlp, st, et, bl = 0, pc = 0, en = 0;
	uint32_t o_l, o_r;
	ul_vec_t *p = NULL;
	nid_t *np = NULL;
	ul_ov_t *z = NULL;
	uc_block_t *b = NULL, tt;

	if(id) {
		kv_pushp(nid_t, x->nid, &np);
		np->n = id_l; MALLOC(np->a, np->n+1); memcpy(np->a, id, id_l); np->a[id_l] = '\0';
	}

	if(str||str_l) {
		if(rid == NULL) {
			kv_pushp(ul_vec_t, *x, &p);
			memset(p, 0, sizeof(*p));
		} else {
			if((*rid) >= x->m) kv_resize(ul_vec_t, *x, (*rid) + 1);
			if((*rid) >= x->n) {
				memset(x->a+x->n, 0, sizeof(*p)*((*rid) + 1 - x->n));
				x->n = (*rid) + 1;
			}
			p = &(x->a[(*rid)]);
		}
		// if((*rid) == 23) fprintf(stderr, "#rid->%lu, on->%ld\n", *rid, on);

		p->bb.n = p->N_site.n = p->r_base.n = 0; p->dd = 0;
		p->rlen = str_l; 
		// fprintf(stderr, "str_l->%ld, str->%u\n", str_l, str?1:0);
		
		if(o == NULL || on == 0) on = 0; en = 0;
		for (i = on-1, st = et = str_l; i >= 0; i--) {
			z = &(o[i]);
			if(z->el) {
				if((z->tn&((uint32_t)(0x80000000)))) {
					mine = MIN(et, ((int64_t)z->qe)); maxs = MAX(st, ((int64_t)z->qs));
					ovlp = mine - maxs; 
					if(ovlp < 0) bl -= ovlp;
					if(save_bases && ovlp < 0) {///push original bases
						kv_pushp(uc_block_t, p->bb, &b);
						b->hid = 0/**x->mm**/; b->rev = 0; b->base = 1; b->pchain = 0; b->el = 0;
						b->qe = maxs; b->qs = b->qe + ovlp; //bl += (b->qe-b->qs);
						o_l = (b->qs >= UL_FLANK?UL_FLANK:b->qs);
						o_r = ((str_l-b->qe)>=UL_FLANK?UL_FLANK:(str_l-b->qe));
						b->pidx = b->pdis = b->aidx = (uint32_t)-1; 
						b->hid |= (o_l<<15); b->hid |= o_r;
						b->qs -= o_l; b->qe += o_r;
						b->ts = p->r_base.n; b->te = b->ts + B4L(b->qe-b->qs);
						kv_resize(uint8_t, p->r_base, b->te); p->r_base.n = b->te;
						// fprintf(stderr, "\n+rid->%lu, str_l->%ld, b->qs->%u, b->qe->%u\n", *rid, str_l, b->qs, b->qe);
						ha_encode_base(p->r_base.a+b->ts, str+b->qs, b->qe-b->qs, &(p->N_site), b->qs);
					} 

					st = MIN(st, z->qs);
				}

				///push ovlp bases
				kv_pushp(uc_block_t, p->bb, &b);
				b->hid = (z->tn<<1)>>1; b->rev = z->rev; b->base = 0; b->el = z->el;
				b->pchain = ((z->tn&((uint32_t)(0x80000000)))?1:0);
				b->qs = z->qs; b->qe = z->qe;
				b->ts = z->ts; b->te = z->te;
				if(b->pchain) pc++;				
				b->pidx = i; b->pdis = b->aidx = (uint32_t)-1; 
				en++; z->qn = p->bb.n - 1;
			}
		}

		if(st > 0) bl += st;
		if(save_bases && st > 0) {///push original bases
			kv_pushp(uc_block_t, p->bb, &b);
			b->hid = 0/**x->mm**/; b->rev = 0; b->base = 1; b->pchain = 0; b->el = 0;
			b->qe = st; b->qs = 0; //bl += (b->qe-b->qs);
			o_l = (b->qs >= UL_FLANK?UL_FLANK:b->qs);
			o_r = ((str_l-b->qe)>=UL_FLANK?UL_FLANK:(str_l-b->qe));
			b->pidx = b->pdis = b->aidx = (uint32_t)-1; 
			b->hid |= (o_l<<15); b->hid |= o_r;
			b->qs -= o_l; b->qe += o_r;
			b->ts = p->r_base.n; b->te = b->ts + B4L(b->qe-b->qs);
			kv_resize(uint8_t, p->r_base, b->te); p->r_base.n = b->te;
			// if(!str) fprintf(stderr, "-rid->%lu, st->%ld, str_l->%ld\n", *rid, st, str_l);
			ha_encode_base(p->r_base.a+b->ts, str+b->qs, b->qe-b->qs, &(p->N_site), b->qs);	
			// push_subblock_original_bases(str, x, p, end, str_l, 321);//for debug
		}

		if(pc > 0) p->dd = 3;
		if((pc == en) && ((str_l-bl) > (str_l*p_chain_rate))) p->dd = 2;
		if((pc == en) && (bl == 0)) p->dd = 1;
		// debug_append_ul_t(o, on, p);
		// char *sst = NULL; CALLOC(sst, str_l);//for debug
		// retrieve_ul_t(NULL, sst, x, rid?*rid:x->n-1, 0, 0, -1);
		// if(memcmp(sst, str, str_l)) {
		// 	fprintf(stderr, "ap-Wrong read, id: %ld, [%d, %ld)\n", (int64_t)(rid?*rid:x->n-1), 0, str_l);
		// 	for (i = 0; i < str_l; i++) {
		// 		if(sst[i] != str[i]) fprintf(stderr, "[%ld] input:%c, decompress:%c\n", i, str[i], sst[i]);
		// 	}
		// }
		// free(sst);
		ovlp = p->bb.n>>1;
		for (i = 0; i < ovlp; i++) {
			tt = p->bb.a[i]; 
			p->bb.a[i] = p->bb.a[p->bb.n-i-1];
			p->bb.a[p->bb.n-i-1] = tt;
			if(p->bb.a[i].pidx!=(uint32_t)-1) {
				o[p->bb.a[i].pidx].qn = p->bb.n-o[p->bb.a[i].pidx].qn-1;
				p->bb.a[i].pidx = (uint32_t)-1;
			}
			if(p->bb.a[p->bb.n-i-1].pidx!=(uint32_t)-1) {
				o[p->bb.a[p->bb.n-i-1].pidx].qn = p->bb.n-o[p->bb.a[p->bb.n-i-1].pidx].qn-1;
				p->bb.a[p->bb.n-i-1].pidx = (uint32_t)-1;
			}
		}
		if(((uint32_t)p->bb.n)&1) {
			o[p->bb.a[i].pidx].qn = p->bb.n-o[p->bb.a[i].pidx].qn-1;
			p->bb.a[i].pidx = (uint32_t)-1;
		}
		
		determine_chain_distance(o, on, p, uopt->sources, uopt->max_hang, uopt->min_ovlp, *rid);
	}
}

void retrieve_ul_t(UC_Read* i_r, char *i_s, all_ul_t *ref, uint64_t ID, uint8_t strand, int64_t s, int64_t l) {
	ul_vec_t *p = &(ref->a[ID]); 
	if(l < 0) l = p->rlen;
	uc_block_t *b = NULL;
	char *r = NULL;
	uint8_t *src = NULL;
	int64_t k, i, a_n, e = s + l, ssp, sep, sl, rts, rte;
	int64_t offset, begLen, tailLen, src_i, des_i;
	if(i_r) {
		i_r->length = l; i_r->RID = ID;
		if(i_r->length > i_r->size) {
			i_r->size = i_r->length;
			i_r->seq = (char*)realloc(i_r->seq,sizeof(char)*(i_r->size));
		}
		r = i_r->seq;
	}
	if(i_s) r = i_s;
	
	

	if(strand == 0) {
		for (k = 0, des_i = 0; k < (int64_t)p->bb.n; k++) {
			b = &(p->bb.a[k]); 
			if(b->qe <= s) continue;
			if(b->qs >= e) break;
			src = p->r_base.a + b->ts; 
			ssp = MAX(s, b->qs) - b->qs; 
			sep = MIN(e, b->qe) - b->qs;
			sl = sep - ssp;
			
			if(b->base/**b->hid&ref->mm**/){///original bases
				offset = ssp&3; 
				begLen = 4-offset;
				if(begLen > sl) begLen = sl;
				tailLen = (sl-begLen)&3;
				a_n = (sl - begLen - tailLen)>>2;

				src_i = ssp; i = 0;
				if(begLen > 0) {
					memcpy(r+des_i, bit_t_seq_table[src[src_i>>2]]+offset, begLen);
					des_i += begLen; src_i += begLen;
				}

				for (i = 0; i < a_n; i++) {
					memcpy(r+des_i, bit_t_seq_table[src[src_i>>2]], 4);
					des_i += 4; src_i += 4;
				}

				if(tailLen > 0) {
					memcpy(r+des_i, bit_t_seq_table[src[src_i>>2]], tailLen);
					des_i += tailLen; src_i += tailLen;
				}
			} else {///ovlps
				if(b->rev == 0) {
					recover_UC_Read_sub_region(r+des_i, b->ts + ssp, sep - ssp, b->rev, ref->hR, b->hid);
				} else{
					rts = b->ts + (b->qe - sep); rte = rts + sep - ssp;
					recover_UC_Read_sub_region(r+des_i, Get_READ_LENGTH((*ref->hR), b->hid) - rte, 
					sep - ssp, b->rev, ref->hR, b->hid);
				}
				des_i += sep - ssp;
			}
		}
		for (k = 0; k < (int64_t)p->N_site.n; k++) {
			if(p->N_site.a[k] >= s && p->N_site.a[k] < e){
				r[p->N_site.a[k]-s] = 'N';
			}
			else if(p->N_site.a[k] >= e) {
				break;
			}
		}
	} else {
		sep = p->rlen - s; 
		ssp = p->rlen - e;
		s = ssp; e = sep; 
		///[s, e)
		for (k = ((int64_t)p->bb.n)-1, des_i = 0; k >= 0; k--) {
			b = &(p->bb.a[k]); 
			if(b->qe <= s) break;
			if(b->qs >= e) continue;
            src = p->r_base.a + b->ts; 
            ssp = MAX(s, b->qs) - b->qs; 
            sep = MIN(e, b->qe) - b->qs;
            sl = sep - ssp;
			///[ssp, sep)
			
			if(b->base/**b->hid&ref->mm**/){///original bases
				begLen = sep&3;
        		offset = 4 - begLen;
				if(begLen > sl) begLen = sl;
				tailLen = (sl-begLen)&3;
				a_n = (sl - begLen - tailLen)>>2;
				src_i = sep-1; i = 0;

				if(begLen > 0) {
					memcpy(r+des_i, bit_t_seq_table_rc[src[src_i>>2]]+offset, begLen);
					des_i += begLen; src_i -= begLen;
				}

				for (i = 0; i < a_n; i++) {
					memcpy(r+des_i, bit_t_seq_table_rc[src[src_i>>2]], 4);
					des_i += 4; src_i -= 4;
				}

				if(tailLen > 0) {
					memcpy(r+des_i, bit_t_seq_table_rc[src[src_i>>2]], tailLen);
					des_i += tailLen; src_i -= tailLen;
				}
			} else {///ovlps
				if(b->rev == 0) {///b->rev != strand
					rts = b->ts + ssp; rte = rts + sep - ssp;
					recover_UC_Read_sub_region(r+des_i, Get_READ_LENGTH((*ref->hR), b->hid) - rte, sep - ssp, 1, ref->hR, b->hid);
				} else {///b->rev == strand
					rts = b->ts + (b->qe - sep); rte = rts + sep - ssp;
					recover_UC_Read_sub_region(r+des_i, rts, sep - ssp, 0, ref->hR, b->hid);
				}
				des_i += sep - ssp;
			}
		}

		sep = p->rlen - s; 
		ssp = p->rlen - e;
		s = ssp; e = sep; 
		for (k = 0; k < (int64_t)p->N_site.n; k++) {
			sl = p->rlen - p->N_site.a[k] - 1;
			if(sl >= s && sl < e) r[sl-s] = 'N';
			else if(sl < s) {
				break;
			}
		}		
	}
}


void retrieve_u_seq(UC_Read* i_r, char* i_s, ma_utg_t *u, uint8_t strand, int64_t s, int64_t l, void *km)
{
	if(u->m == 0 || u->n == 0) return;
	if(l < 0) l = u->len;
	char *r = NULL, *a = NULL;
	int64_t e = s + l, ssp, sep, rs, re, des_i;
	uint64_t k, rId, ori, r_l;
	if(i_r) {
        i_r->length = l; i_r->RID = 0;
        if(i_r->length > i_r->size) {
            i_r->size = i_r->length;
			if(!km) REALLOC(i_r->seq, i_r->size);
			else KREALLOC(km, i_r->seq, i_r->size);
            // i_r->seq = (char*)realloc(i_r->seq,sizeof(char)*(i_r->size));
        }
        r = i_r->seq;
    }
    if(i_s) r = i_s;

	if(strand == 1) {
		sep = u->len - s; 
        ssp = u->len - e;
        s = ssp; e = sep; 
	}
	for (k = l = des_i = 0; k < u->n; k++) {
		rId = u->a[k]>>33;
		ori = u->a[k]>>32&1;
		r_l = (uint32_t)u->a[k];
		if(r_l == 0) continue;
		ssp = l; sep = l + r_l;
		l += r_l;
		if(sep <= s) continue;
		if(ssp >= e) break;
		rs = MAX(ssp, s); re = MIN(sep, e); 
		a = r + des_i; des_i += re - rs;
		recover_UC_Read_sub_region(a, rs-ssp, re-rs, ori, &R_INF, rId);
	}
	if(strand == 1) {
		char t;
		re = (e - s);
		l = re>>1; 
		for (k = 0; k < (uint64_t)l; k++) {
			des_i = re - k - 1;
			t = r[des_i];
			r[des_i] = RC_CHAR(r[k]);
			r[k] = RC_CHAR(t);
		}
		if(re&1) r[l] = RC_CHAR(r[l]);
	}
}

uint32_t retrieve_u_cov(const ul_idx_t *ul, uint64_t id, uint8_t strand, uint64_t pos, uint8_t dir, int64_t *pi)
{
	uint64_t *a = ul->cc->interval.a + ul->cc->idx[id], cc = 0, ff = 0;
	int64_t a_n = ul->cc->idx[id+1]-ul->cc->idx[id], k = 0, cc_i = pi? *pi:0;
	if(a_n == 0) return 0;
	if(cc_i + 1 >= a_n || cc_i < 0) cc_i = 0;
	if(strand) pos = ul->ug->u.a[id].len - pos - 1;
	if(dir == 0) {
		for (k = cc_i; k + 1 < a_n; k++) {
			if(pos>=(a[k]>>32) && pos<(a[k+1]>>32)) {
				cc = (uint32_t)a[k];
				ff = 1;
				break;
			}
		}

		if(ff == 0) {
			for (k = 0; k < cc_i; k++) {
					if(pos>=(a[k]>>32) && pos<(a[k+1]>>32)) {
					cc = (uint32_t)a[k];
					ff = 1;
					break;
				}
			}
		}
	} else {
		for (k = cc_i; k >= 0; k--) {
			if(pos>=(a[k]>>32) && pos<(a[k+1]>>32)) {
				cc = (uint32_t)a[k];
				ff = 1;
				break;
			}
		}

		if(ff == 0) {
			for (k = cc_i+1; k + 1 < a_n; k++) {
				if(pos>=(a[k]>>32) && pos<(a[k+1]>>32)) {
					cc = (uint32_t)a[k];
					ff = 1;
					break;
				}
			}
		}
	}
	
	if(pi) *pi = ff?k:0;
	return cc;
}

uint64_t retrieve_u_cov_region(const ul_idx_t *ul, uint64_t id, uint8_t strand, uint64_t s, uint64_t e, int64_t *pi)
{
    uint64_t *a = ul->cc->interval.a + ul->cc->idx[id], cc = 0, o = 0, tk, ts, te, tcc = 0;
    int64_t a_n = ul->cc->idx[id+1]-ul->cc->idx[id], k = 0, cc_i = pi? *pi:0;
    if(a_n == 0) return 0;
    if(cc_i + 1 >= a_n || cc_i < 0) cc_i = 0;
    if(strand) {
		tk = s;
		s = ul->ug->u.a[id].len - e;
		e = ul->ug->u.a[id].len - tk;
	}
	// fprintf(stderr,"\nul->ug->u.a[id].len:%u, fs:%lu, fe:%lu\n", ul->ug->u.a[id].len, a[k]>>32, a[k+1]>>32);
	k = cc_i; tk = s;
	if(tk < (a[k]>>32)) {
		for (; k >= 0; k--) {
            if(tk>=(a[k]>>32) && tk<(a[k+1]>>32)) break;
        }
	} else if(tk >= (a[k+1]>>32)) {
		for (; k + 1 < a_n; k++) {
            if(tk>=(a[k]>>32) && tk<(a[k+1]>>32)) break;
        }
	}
	if(pi) *pi = k;

	

	for (; k + 1 < a_n; k++) {
		ts = a[k]>>32; te = a[k+1]>>32; cc = (uint32_t)a[k];
		o = (MIN(e, te) > MAX(s, ts))?(MIN(e, te)-MAX(s, ts)):0;
		tcc += o*cc;
		// fprintf(stderr, ">>k:%ld, s:%lu, e:%lu, ts:%lu, te:%lu, o:%lu, cc:%lu\n", k, s, e, ts, te, o, cc);
		if(e>=(a[k]>>32) && e<(a[k+1]>>32)) break;
	}
    
    // if(s == 54201 && e == 58376 && id == 492) fprintf(stderr, "tcc:%lu\n", tcc);
    return tcc;
}


uint64_t retrieve_r_cov_region(const ul_idx_t *ul, uint64_t id, uint8_t strand, uint64_t s, uint64_t e, int64_t *pi)
{
    uint64_t *a = ul->cr->interval.a + ul->cr->idx[id], cc = 0, o = 0, tk, ts, te, tcc = 0;
    int64_t a_n = ul->cr->idx[id+1]-ul->cr->idx[id], k = 0, cc_i = pi? *pi:0;
    if(a_n == 0) return e>=s?e-s:0;
    if(cc_i + 1 >= a_n || cc_i < 0) cc_i = 0;
    if(strand) {
		tk = s;
		s = ul->ug->u.a[id].len - e;
		e = ul->ug->u.a[id].len - tk;
	}
	// fprintf(stderr,"\nul->ug->u.a[id].len:%u, fs:%lu, fe:%lu\n", ul->ug->u.a[id].len, a[k]>>32, a[k+1]>>32);
	k = cc_i; tk = s;
	if(tk < (a[k]>>32)) {
		for (; k >= 0; k--) {
            if(tk>=(a[k]>>32) && tk<(a[k+1]>>32)) break;
        }
	} else if(tk >= (a[k+1]>>32)) {
		for (; k + 1 < a_n; k++) {
            if(tk>=(a[k]>>32) && tk<(a[k+1]>>32)) break;
        }
	}
	if(pi) *pi = k;
	if(k + 1 >= a_n) return e>=s?e-s:0;
	if(k < 0) k = 0;
	
	///note: for unitig coverage, all 
	for (; k + 1 < a_n; k++) {
		// fprintf(stderr, "[M::%s] k:%ld, a_n:%ld\n", __func__, k, a_n);
		ts = a[k]>>32; te = a[k+1]>>32; cc = (uint32_t)a[k];
		o = (MIN(e, te) > MAX(s, ts))?(MIN(e, te)-MAX(s, ts)):0;
		tcc += o*cc;
		// fprintf(stderr, ">>k:%ld, s:%lu, e:%lu, ts:%lu, te:%lu, o:%lu, cc:%lu\n", k, s, e, ts, te, o, cc);
		if(e>=(a[k]>>32) && e<(a[k+1]>>32)) break;
		if(e<=(a[k]>>32)) break;///if may happend when [s, e) does not overlap with the first region
	}
    
    
    return tcc + (e>=s?e-s:0);
}

uint32_t produce_u_cov(ul_idx_t *ul, uint64_t id, uint8_t strand, uint64_t pos, ma_hit_t_alloc* src, int64_t min_ovlp, int64_t max_hang, int64_t gap_fuzz,
uint8_t *sset, kvec_t_u64_warp *buf)
{
	uint64_t k, l, i, z, s = 0, e = 0, qn, tn, qs, qe;
	ma_utg_t *u = NULL;
	int64_t dp, r; 
	asg_arc_t t;

	if(strand == (uint8_t)-1 || pos == (uint64_t)-1) {
		u = &(ul->ug->u.a[id]);
		buf->a.n = 0; kv_resize(uint64_t, buf->a, u->n*2);
		for (k = l = 0; k < u->n; k++) {
			kv_push(uint64_t, buf->a, l<<1);
			kv_push(uint64_t, buf->a, ((l + Get_READ_LENGTH(R_INF, u->a[k]>>33))<<1)|1);


			i = u->a[k]>>33;///rid
			for (z = 0; z < src[i].length; z++) {
				if(!src[i].buffer[z].el) continue;
				qn = Get_qn(src[i].buffer[z]); tn = Get_tn(src[i].buffer[z]);
				if(sset[tn]) continue;
				if((Get_qe(src[i].buffer[z]) - Get_qs(src[i].buffer[z])) < min_ovlp) continue;
				if((Get_te(src[i].buffer[z]) - Get_ts(src[i].buffer[z])) < min_ovlp) continue;
				r = ma_hit2arc(&(src[i].buffer[z]), Get_READ_LENGTH(R_INF, qn), Get_READ_LENGTH(R_INF, tn), 
				max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
				if(r != MA_HT_TCONT) continue;///tn is contained
				if(((u->a[k]>>32)&1) == 0) {
					qs = Get_qs(src[i].buffer[z]); qe = Get_qe(src[i].buffer[z]);
				} else {
					qs = (Get_READ_LENGTH(R_INF, i)) - Get_qe(src[i].buffer[z]);
					qe = (Get_READ_LENGTH(R_INF, i)) - Get_qs(src[i].buffer[z]);
				}
				kv_push(uint64_t, buf->a, (l+qs)<<1);
				kv_push(uint64_t, buf->a, ((l+qe)<<1)|1);
			}

			l += (uint32_t)u->a[k];
		}
		sort_kvec_t_u64_warp(buf, 0);
		return 0;
	}

	if(strand) pos = ul->ug->u.a[id].len - pos - 1;
	for (k = 0, dp = 0, s = e = 0; k < buf->a.n; k++) {        
		e = buf->a.a[k]>>1;
		// fprintf(stderr, "[M::%s::k:%lu] [s, e)->[%lu, %lu), dp->%ld\n", __func__, k, s, e, dp);
		if(pos >= s && pos < e) break;
		s = buf->a.a[k]>>1;
        if (buf->a.a[k]&1) --dp;
        else ++dp;
	}

	return dp;
}



void produce_u_seq(char* r, ma_utg_t *u, UC_Read *buf)
{
	if(u->m == 0 || u->n == 0) return;
	uint32_t j, k, l = 0;
    uint32_t rId, ori, start, eLen, readLen;
	char *readS = NULL;
	memset(r, 'N', u->len);
	for (j = 0; j < u->n; ++j) {
		rId = u->a[j]>>33;
		///uId = i;
		ori = u->a[j]>>32&1;
		start = l;
		eLen = (uint32_t)u->a[j];
		l += eLen;

		if(eLen == 0) continue;
		recover_UC_Read(buf, &R_INF, rId);

		readS = buf->seq;
		readLen = Get_READ_LENGTH(R_INF, rId);
		
		if (!ori) // forward strand
		{
			for (k = 0; k < eLen; k++)
			{
				r[start + k] = readS[k];
			}
		}
		else
		{
			for (k = 0; k < eLen; k++)
			{
				uint8_t c = (uint8_t)readS[readLen - 1 - k];
				r[start + k] = c >= 128? 'N' : RC_CHAR(c);
			}
		}
	}
}

void debug_retrieve_rc_sub(const ug_opt_t *uopt, all_ul_t *ref, const All_reads *R_INF, ul_idx_t *ul, uint32_t n_step)
{
	uint64_t i, step, s, e, occ, qc, rc;
	UC_Read f, r; 
	init_UC_Read(&f); init_UC_Read(&r);
	kvec_t(char) ss; kv_init(ss);
	if(ref) {
		for (i = 0, occ = 0; i < ref->n; i++) {
			retrieve_ul_t(&f, NULL, ref, i, 0, 0, -1);
			retrieve_ul_t(&r, NULL, ref, i, 1, 0, -1);

			kv_resize(char, ss, ref->a[i].rlen); 
			ss.n = ref->a[i].rlen;
			memcpy(ss.a, r.seq, ss.n); 
			reverse_complement(ss.a, ss.n);
			if(memcmp(ss.a, f.seq, ss.n)) {
				fprintf(stderr, "1-Wrong whole reverse-read, id: %lu\n", i);
				// for (s = 0; s < ss.n; s++) {
				// 	if(ss.a[s] != f.seq[s]) {
				// 		fprintf(stderr,"s:%lu, ss.a[s]:%c, f.seq[s]:%c, r.seq[rs]:%c\n", s, ss.a[s], f.seq[s], r.seq[ref->a[i].rlen - s - 1]);
				// 	}
				// }
			}

			step = ss.n/n_step;
			if(step <= 0) step = 1;

			for (s = 0; s < ss.n; s += step) {
				e = MIN(s+step, ss.n);
				retrieve_ul_t(NULL, ss.a, ref, i, 0, s, e-s);
				if(memcmp(ss.a, f.seq + s, e - s)) {
					fprintf(stderr, "1-Wrong sub forward-read, id: %lu, [%lu, %lu)\n", i, s, e);
				}

				retrieve_ul_t(NULL, ss.a, ref, i, 1, s, e-s);
				if(memcmp(ss.a, r.seq + s, e - s)) {
					fprintf(stderr, "1-Wrong sub reverse-read, id: %lu, [%lu, %lu)\n", i, s, e);
					// uint64_t dk; char *tf = ss.a; char *rf = r.seq + s;
					// for (dk = 0; dk < e - s; dk++) {
					// 	if(tf[dk] != rf[dk]) {
					// 		fprintf(stderr,"dk:%lu, tf[dk]:%c, rf[dk]:%c\n", 
					// 					dk, tf[dk], rf[dk]);
					// 	}
					// }
				}
				occ++;
			}
		}
		fprintf(stderr, "[M::%s::# checking: %lu] ==> all_ul_t\n", __func__, occ);
	}


	if(R_INF) {
		for (i = 0, occ = 0; i < R_INF->total_reads; i++) {
			recover_UC_Read(&f, R_INF, i);
			recover_UC_Read_RC(&r, (All_reads*)R_INF, i);
			kv_resize(char, ss, Get_READ_LENGTH((*R_INF), i)); 
			ss.n = Get_READ_LENGTH((*R_INF), i);
			memcpy(ss.a, r.seq, ss.n); 
			reverse_complement(ss.a, ss.n);
			if(memcmp(ss.a, f.seq, ss.n)) fprintf(stderr, "2-Wrong whole reverse-read, id: %lu\n", i);

			step = ss.n/n_step;
			if(step <= 0) step = 1;
			for (s = 0; s < ss.n; s += step) {
				e = MIN(s+step, ss.n);
				recover_UC_Read_sub_region(ss.a, s, e-s, 0, (All_reads*)R_INF, i);
				if(memcmp(ss.a, f.seq + s, e - s)) fprintf(stderr, "2-Wrong sub forward-read, id: %lu, [%lu, %lu)\n", i, s, e);

				recover_UC_Read_sub_region(ss.a, s, e-s, 1, (All_reads*)R_INF, i);
				if(memcmp(ss.a, r.seq + s, e - s)) fprintf(stderr, "2-Wrong sub reverse-read, id: %lu, [%lu, %lu)\n", i, s, e);
				occ++;
			}
		}
		fprintf(stderr, "[M::%s::# checking: %lu] ==> All_reads\n", __func__, occ);
	}

	if(ul) {

		uint8_t *sset = NULL; CALLOC(sset, R_INF->total_reads);
		for (i = 0; i < ul->ug->u.n; i++) {
			for (s = 0; s < ul->ug->u.a[i].n; s++){
				sset[ul->ug->u.a[i].a[s]>>33] = 1;
			}
		}

		kvec_t_u64_warp buf; memset(&buf, 0, sizeof(buf)); int64_t pi = 0;
		// produce_u_cov(ul, 0, (uint8_t)-1, (uint64_t)-1, &buf);
		// produce_u_cov(ul, 0, 0, ul->ug->u.a[0].len>>1, &buf);

		for (i = 0, occ = 0; i < ul->ug->u.n; i++) {
			kv_resize(char, ss, ul->ug->u.a[i].len); ss.n = ul->ug->u.a[i].len;
			retrieve_u_seq(&f, NULL, &(ul->ug->u.a[i]), 0, 0, -1, NULL);
			produce_u_seq(ss.a, &(ul->ug->u.a[i]), &r);
			if(memcmp(ss.a, f.seq, ss.n)) fprintf(stderr, "4-Wrong whole reverse-read, id: %lu\n", i);
			retrieve_u_seq(&r, NULL, &(ul->ug->u.a[i]), 1, 0, -1, NULL);
			
			memcpy(ss.a, r.seq, ss.n); 
			reverse_complement(ss.a, ss.n);
			if(memcmp(ss.a, f.seq, ss.n)) fprintf(stderr, "3-Wrong whole reverse-read, id: %lu\n", i);

			produce_u_cov(ul, i, (uint8_t)-1, (uint64_t)-1, uopt->sources, uopt->min_ovlp, uopt->max_hang, uopt->gap_fuzz, sset, &buf);
			step = ss.n/n_step;
    		if(step <= 0) step = 1;

			for (s = 0; s < ss.n; s += step) {
				e = MIN(s+step, ss.n);
				retrieve_u_seq(NULL, ss.a, &(ul->ug->u.a[i]), 0, s, e-s, NULL);
				if(memcmp(ss.a, f.seq + s, e - s)) fprintf(stderr, "3-Wrong sub forward-read, id: %lu, [%lu, %lu)\n", i, s, e);

				retrieve_u_seq(NULL, ss.a, &(ul->ug->u.a[i]), 1, s, e-s, NULL);
				if(memcmp(ss.a, r.seq + s, e - s)) fprintf(stderr, "3-Wrong sub reverse-read, id: %lu, [%lu, %lu)\n", i, s, e);
				
				/**if(ul->ug->u.a[i].n > 1)**/ {
					rc = produce_u_cov(ul, i, 0, s, uopt->sources, uopt->min_ovlp, uopt->max_hang, uopt->gap_fuzz, sset, &buf);

					qc = retrieve_u_cov(ul, i, 0, s, 0, &pi);
					if(rc != qc) fprintf(stderr, "4-Wrong coverage, id: %lu, rc:%lu, qc:%lu, pos:%lu\n", i, rc, qc, s);
					pi++;
					qc = retrieve_u_cov(ul, i, 0, s, 0, &pi);
					if(rc != qc) fprintf(stderr, "4-Wrong coverage, id: %lu, rc:%lu, qc:%lu, pos:%lu\n", i, rc, qc, s);
					pi--;
					qc = retrieve_u_cov(ul, i, 0, s, 0, &pi);
					if(rc != qc) fprintf(stderr, "4-Wrong coverage, id: %lu, rc:%lu, qc:%lu, pos:%lu\n", i, rc, qc, s);
					pi++;
					qc = retrieve_u_cov(ul, i, 0, s, 1, &pi);
					if(rc != qc) fprintf(stderr, "4-Wrong coverage, id: %lu, rc:%lu, qc:%lu, pos:%lu\n", i, rc, qc, s);
					pi--;
					qc = retrieve_u_cov(ul, i, 0, s, 1, &pi);
					if(rc != qc) fprintf(stderr, "4-Wrong coverage, id: %lu, rc:%lu, qc:%lu, pos:%lu\n", i, rc, qc, s);




					rc = produce_u_cov(ul, i, 1, s, uopt->sources, uopt->min_ovlp, uopt->max_hang, uopt->gap_fuzz, sset, &buf);
					qc = retrieve_u_cov(ul, i, 1, s, 0, &pi);
					if(rc != qc) fprintf(stderr, "4-Wrong coverage, id: %lu, rc:%lu, qc:%lu, pos:%lu\n", i, rc, qc, s);
					pi++;
					qc = retrieve_u_cov(ul, i, 1, s, 0, &pi);
					if(rc != qc) fprintf(stderr, "4-Wrong coverage, id: %lu, rc:%lu, qc:%lu, pos:%lu\n", i, rc, qc, s);
					pi--;
					qc = retrieve_u_cov(ul, i, 1, s, 0, &pi);
					if(rc != qc) fprintf(stderr, "4-Wrong coverage, id: %lu, rc:%lu, qc:%lu, pos:%lu\n", i, rc, qc, s);
					pi++;
					qc = retrieve_u_cov(ul, i, 1, s, 1, &pi);
					if(rc != qc) fprintf(stderr, "4-Wrong coverage, id: %lu, rc:%lu, qc:%lu, pos:%lu\n", i, rc, qc, s);
					pi--;
					qc = retrieve_u_cov(ul, i, 1, s, 1, &pi);
					if(rc != qc) fprintf(stderr, "4-Wrong coverage, id: %lu, rc:%lu, qc:%lu, pos:%lu\n", i, rc, qc, s);
				}
				occ++;
			}
		}
	
		kv_destroy(buf.a); free(sset);
		fprintf(stderr, "[M::%s::# checking: %lu] ==> ma_utg_v\n", __func__, occ);
	}

	destory_UC_Read(&f); destory_UC_Read(&r);
	kv_destroy(ss);
}


void write_compress_base_disk(FILE *fp, uint64_t ul_rid, char *str, uint32_t len, ul_vec_t *buf)
{
	buf->N_site.n = buf->bb.n = 0; buf->rlen = len; buf->r_base.n = B4L(len); 
	kv_resize(uint8_t, buf->r_base, buf->r_base.n); 
	ha_encode_base(buf->r_base.a, str, len, &(buf->N_site), 0);

	fwrite(&ul_rid, sizeof(ul_rid), 1, fp);
	fwrite(&len, sizeof(len), 1, fp);
	fwrite(&buf->N_site.n, sizeof(buf->N_site.n), 1, fp);
	fwrite(buf->N_site.a, sizeof((*buf->N_site.a)), buf->N_site.n, fp);
	fwrite(buf->r_base.a, sizeof((*buf->r_base.a)), buf->r_base.n, fp);
}

int64_t load_compress_base_disk(FILE *fp, uint64_t *ul_rid, char *dest, uint32_t *len, ul_vec_t *buf)
{
	fread(ul_rid, sizeof((*ul_rid)), 1, fp);
	if(feof(fp)) return 0;
	fread(len, sizeof((*len)), 1, fp);
	fread(&buf->N_site.n, sizeof(buf->N_site.n), 1, fp);
	kv_resize(uint32_t, buf->N_site, buf->N_site.n);
	fread(buf->N_site.a, sizeof((*buf->N_site.a)), buf->N_site.n, fp);
	buf->r_base.n = B4L((*len)); kv_resize(uint8_t, buf->r_base, buf->r_base.n); 
	fread(buf->r_base.a, sizeof((*buf->r_base.a)), buf->r_base.n, fp);

	int64_t ssp, sep, sl, offset, begLen, tailLen, a_n, src_i, des_i, i;
	ssp = 0; sep = (*len); sl = sep - ssp;
	offset = ssp&3; begLen = 4-offset;
	if(begLen > sl) begLen = sl;
	tailLen = (sl-begLen)&3;
	a_n = (sl - begLen - tailLen)>>2;

	src_i = ssp; i = 0; des_i = 0;
	if(begLen > 0) {
		memcpy(dest+des_i, bit_t_seq_table[buf->r_base.a[src_i>>2]]+offset, begLen);
		des_i += begLen; src_i += begLen;
	}

	for (i = 0; i < a_n; i++) {
		memcpy(dest+des_i, bit_t_seq_table[buf->r_base.a[src_i>>2]], 4);
		des_i += 4; src_i += 4;
	}

	if(tailLen > 0) {
		memcpy(dest+des_i, bit_t_seq_table[buf->r_base.a[src_i>>2]], tailLen);
		des_i += tailLen; src_i += tailLen;
	}
	return 1;
}


scaf_res_t *init_scaf_res_t(uint32_t n)
{
	scaf_res_t *p = NULL; CALLOC(p, 1);
	p->n = p->m = n; CALLOC(p->a, n);
	return p;
}

void destroy_scaf_res_t(scaf_res_t *p)
{
	if(p) {
		uint32_t k;
		for (k = 0; k < p->m; k++) {
			free(p->a[k].N_site.a); free(p->a[k].r_base.a); free(p->a[k].bb.a);	
		}
		free(p->a);
		free(p);
	}
}