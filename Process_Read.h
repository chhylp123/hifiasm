#ifndef __READ__
#define __READ__

#include<stdint.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
#include "kseq.h"


KSEQ_INIT(gzFile, gzread)

void init_kseq(char* file);
void destory_kseq();
int get_read(kseq_t *s);


#endif
