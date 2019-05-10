#include "Assembly.h"
#include <stdio.h>
#include <stdlib.h>
#include "Process_Read.h"

void Counting()
{
    Read r;
    init_Read(&r);

    long long read_number = 0;

    while (inputRead(&r))
    {
        
        fprintf(stderr,"%s\n",r.name);
        
        fprintf(stderr,"%s\n",r.seq);
        fprintf(stderr,"+\n");
        fprintf(stderr,"%s\n",r.qual);
        

        read_number++;

        clear_read(&r);
    }
    
    fprintf(stdout, "read_number: %lld\n",read_number);



}
