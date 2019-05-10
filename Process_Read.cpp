#include "Process_Read.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>


gzFile fp;
kseq_t *seq;

void init_kseq(char* file)
{
	fp = gzopen(file, "r");
  	seq = kseq_init(fp);
}

void destory_kseq()
{
	kseq_destroy(seq);
  	gzclose(fp);
}


inline void exchage_kstring_t(kstring_t* a, kstring_t* b)
{
	kstring_t tmp;
	tmp = *a;
	*a = *b;
	*b = tmp;
}

int get_read(kseq_t *s)
{
	int l;

	if ((l = kseq_read(seq)) >= 0)
	{
		
		exchage_kstring_t(&seq->comment, &s->comment);
		exchage_kstring_t(&seq->name, &s->name);
		exchage_kstring_t(&seq->qual, &s->qual);
		exchage_kstring_t(&seq->seq, &s->seq);
		
		return 1;
	}
	else
	{
		return 0;
	}
	
	
}







