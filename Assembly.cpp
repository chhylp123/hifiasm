#include "Assembly.h"
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "Process_Read.h"





void Counting()
{


    kseq_t *seq = (kseq_t*)calloc(1, sizeof(kseq_t));	


    long long read_number = 0;

    while (get_read(seq))
    {

        fprintf(stderr,"@%s\n",seq->name.s);
        fprintf(stderr,"%s\n",seq->seq.s);
        fprintf(stderr,"+\n");
        fprintf(stderr,"%s\n",seq->qual.s);

        read_number++;
    }
    
    fprintf(stdout, "read_number: %lld\n",read_number);



}
