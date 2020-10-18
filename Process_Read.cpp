#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include "Process_Read.h"

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
char rc_Table[5] = {'T', 'G', 'C', 'A', 'N'};

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


void recover_UC_Read_sub_region(char* r, long long start_pos, long long length, uint8_t strand, All_reads* R_INF, long long ID)
{
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

	if (r->length + 4 > r->size)
	{
		r->size = r->length + 4;
		r->seq = (char*)realloc(r->seq,sizeof(char)*(r->size));
	}

	uint64_t i = 0;

	while ((long long)i < r->length)
	{
		memcpy(r->seq+i, bit_t_seq_table[src[i>>2]], 4);
		i = i + 4;
	}


	if (R_INF->N_site[ID])
	{
		for (i = 1; i <= R_INF->N_site[ID][0]; i++)
		{
			r->seq[R_INF->N_site[ID][i]] = 'N';
		}
	}

	r->RID = ID;
		
}

void recover_UC_Read_RC(UC_Read* r, All_reads* R_INF, uint64_t ID)
{
	r->length = Get_READ_LENGTH((*R_INF), ID);
	uint8_t* src = Get_READ((*R_INF), ID);

	if (r->length + 4 > r->size)
	{
		r->size = r->length + 4;
		r->seq = (char*)realloc(r->seq,sizeof(char)*(r->size));
	}

	long long last_chr = r->length % 4;
	long long i = r->length / 4 - 1 + (last_chr != 0);
	long long index = 0;

	if(last_chr!=0)
	{
		memcpy(r->seq + index, bit_t_seq_table_rc[src[i]] + 4 - last_chr, last_chr);
		index = last_chr;
		i--;
	}

	while (i >= 0)
	{
		memcpy(r->seq + index, bit_t_seq_table_rc[src[i]], 4);
		i--;
		index = index + 4;
	}

	if (R_INF->N_site[ID])
	{
		for (i = 1; i <= (long long)R_INF->N_site[ID][0]; i++)
		{
			r->seq[r->length - R_INF->N_site[ID][i] - 1] = 'N';
		}
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
	
	free(x->read_name);
	fclose(x->fp);
}