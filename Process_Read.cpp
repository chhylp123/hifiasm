#include "Process_Read.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <fcntl.h>


int format = -1;
FILE *_r_fp = NULL;
gzFile g_r_fp;


int check_file_format(char *fileName)
{
	int i = 0;

	char* file_extention;
	i = strlen(fileName) - 1;
	while (i >= 0 && fileName[i] != '.')
	{
		i--;
	}

	if (i < 0)
	{
		fprintf(stderr, 
			"ERROR: the input read files must be FASTQ format (file extension: .fq or .fastq) or compressed FASTQ format (file extension: .fq.gz or .fastq.gz)!\n");
		exit(0);
	}

	if ((strcmp(".fq", fileName+i) == 0)
		||
		strcmp(".fastq", fileName + i) == 0)
	{
		return FASTQ;
	}
	else if (strcmp(".gz", fileName + i) == 0)
	{
		i--;
		while (i >= 0 && fileName[i] != '.')
		{
			i--;
		}


		if (i < 0)
		{
			fprintf(stderr,
				"ERROR: the input read files must be FASTQ format (file extension: .fq or .fastq) or compressed FASTQ format (file extension: .fq.gz or .fastq.gz)!\n");
			exit(0);
		}

		if ((strcmp(".fq.gz", fileName + i) == 0)
			||
			strcmp(".fastq.gz", fileName + i) == 0)
		{
			return FASTQGZ;
		}
		else
		{
			fprintf(stderr,
				"ERROR: the input read files must be FASTQ format (file extension: .fq or .fastq) or compressed FASTQ format (file extension: .fq.gz or .fastq.gz)!\n");
			exit(0);
		}

	}


	fprintf(stderr,
		"ERROR: the input read files must be FASTQ format (file extension: .fq or .fastq) or compressed FASTQ format (file extension: .fq.gz or .fastq.gz)!\n");
	exit(0);
}


int initiReadAllReads(char *fileName)
{

	format = check_file_format(fileName);

    if (format == FASTQ)
	{
		_r_fp = fopen(fileName, "r");

		if (_r_fp == NULL)
		{
			return 0;
		}

        fprintf(stdout,
				" Read files are in FASTQ format...\n");
	}
	else if (format == FASTQGZ)
	{
		int fd = open(fileName, O_CREAT | O_RDONLY, 0666);
		g_r_fp = gzdopen(fd, "rb");

        fprintf(stdout,
				" Read files are in compressed FASTQ format...\n");
	}
    else
    {
        fprintf(stderr,
			"ERROR: the input read files must be FASTQ format (file extension: .fq or .fastq) or compressed FASTQ format (file extension: .fq.gz or .fastq.gz)!\n");
		exit(0);
    }
    
}

#define BUFFER_SIZE 20000
char Read_Buffer[BUFFER_SIZE];
char true_flag[2];


inline char *Input_line_from_file(char *seq)
{
	if (format == FASTQ)
	{
		return fgets(Read_Buffer, BUFFER_SIZE, _r_fp);
	}
	else
	{
		seq = gzgets(g_r_fp, seq, BUFFER_SIZE);

		return (!gzeof(g_r_fp)) ? true_flag : NULL;

	}
  
}

inline int get_whole_line(str_v* line)
{
	while(Input_line_from_file(Read_Buffer))
	{
		append_str(line, Read_Buffer);

		if (line->centext[line->length - 1] == '\n')
		{
			line->length--;
			line->centext[line->length] = '\0';
			return 1;
		}
		
	}

	return 0;
}


int inputRead(Read *read)
{
	int buffer_length;

	if (get_whole_line(&read->name))
	{
		get_whole_line(&read->seq);
		Input_line_from_file(Read_Buffer);
		get_whole_line(&read->qual);
		return 1;
	}
	
	return 0;

}