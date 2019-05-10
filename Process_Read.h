#ifndef __READ__
#define __READ__

#include<stdint.h>
#include <string.h>
#include <stdlib.h>


#define FASTQ 1
#define FASTQGZ 2

#define NAME_SZIE 25
#define SEQ_SZIE 20000


typedef struct
{
	char* centext;
	uint32_t length;
	uint32_t size;

}str_v;


inline void init_str_v(str_v* buf, uint32_t length)
{
	buf->size = length + 1;
	buf->length = 0;
	buf->centext = (char*)malloc(sizeof(char)*buf->size);
}


inline void append_str(str_v* buf, char* input)
{
	int length = strlen(input);
	buf->length = buf->length + length;
	if (buf->length > buf->size)
	{
		buf->centext = (char*)realloc(buf->centext, buf->length + 1);
	}
	memcpy(buf->centext + buf->length - length, input, length);
	buf->centext[buf->length] = '\0';
}

inline void clear_str(str_v* buf)
{
	buf->length = 0;
}

typedef struct
{
	str_v name;
	str_v seq;
	str_v qual;

} Read;

inline void init_Read(Read* r)
{
	init_str_v(&r->name, NAME_SZIE);
	init_str_v(&r->seq, SEQ_SZIE);
	init_str_v(&r->qual, SEQ_SZIE);
}

inline void clear_read(Read* r)
{
	clear_str(&r->name);
	clear_str(&r->seq);
	clear_str(&r->qual);
}


int initiReadAllReads(char *fileName);
int inputRead(Read *read);

#endif
