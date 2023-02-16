#ifndef __LEVENSHTEIN__
#define __LEVENSHTEIN__

#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include "emmintrin.h"
#include "nmmintrin.h"
#include "smmintrin.h"
#include <immintrin.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "kvec.h"

extern const unsigned char seq_nt4_table[256];
typedef uint64_t Word;
typedef uint32_t Word_32;
typedef struct {size_t n, m; uint16_t *a; } asg16_v;

inline void get_error(int t_length, int errthold, int init_err, Word VP, Word VN, 
unsigned int* return_err, int* back_site)
{
	(*return_err) = (unsigned int)-1;
	int site = t_length - 1;
	int return_site = -1;
	///in most cases, p_length should be t_length + 2 * errthold
	///int available_i = p_length - t_length;
	int available_i = 2 * errthold;
	

	if ((init_err <= errthold) && ((unsigned int)init_err <= (*return_err)))
	{
		(*return_err) = init_err;
		return_site = site;
	}


	int i = 0;
	unsigned int ungap_error = (unsigned int)-1;

	while (i < available_i)
	{
		init_err = init_err + ((VP >> i)&(Word)1);
		init_err = init_err - ((VN >> i)&(Word)1);
		++i;

		if ((init_err <= errthold) && ((unsigned int)init_err <= *return_err))
		{
			*return_err = init_err;
			return_site = site + i;
		}

		/****************************may have bugs********************************/
		if(i == errthold)
		{
			ungap_error = init_err;
		}
		/****************************may have bugs********************************/
	}

	/****************************may have bugs********************************/
	if((ungap_error<=(unsigned int)errthold) && (ungap_error == (*return_err)))
	{
		return_site = site + errthold;
	}
	/****************************may have bugs********************************/

	(*back_site) = return_site;
}

inline int Reserve_Banded_BPM_Extension
(char *pattern, int p_length, char *text, int t_length, unsigned short errthold, 
unsigned int* return_err, int* return_p_end, int* return_t_end)
{
	(*return_err) = (unsigned int)-1;
	(*return_p_end) = -1;
	(*return_t_end) = -1;

	Word Peq[256];

	unsigned int line_error = (unsigned int)-1;
	int return_site;
	int band_length = (errthold << 1) + 1;
	int i = 0;
	Word tmp_Peq_1 = (Word)1;

	Peq[(uint8_t)'A'] = (Word)0;
	Peq[(uint8_t)'T'] = (Word)0;
	Peq[(uint8_t)'G'] = (Word)0;
	Peq[(uint8_t)'C'] = (Word)0;


	Word Peq_A;
	Word Peq_T;
	Word Peq_C;
	Word Peq_G;

	///band_length = 2k + 1
	for (i = 0; i<band_length; i++)
	{
		Peq[(uint8_t)pattern[i]] = Peq[(uint8_t)pattern[i]] | tmp_Peq_1;
		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	Peq_A = Peq[(uint8_t)'A'];
	Peq_C = Peq[(uint8_t)'C'];
	Peq_T = Peq[(uint8_t)'T'];
	Peq_G = Peq[(uint8_t)'G'];


	memset(Peq, 0, sizeof(Word)* 256);


	Peq[(uint8_t)'A'] = Peq_A;
	Peq[(uint8_t)'C'] = Peq_C;
	Peq[(uint8_t)'T'] = Peq_T;
	Peq[(uint8_t)'G'] = Peq_G;


	

	Word Mask = ((Word)1 << (errthold << 1));

	Word VP = 0;
	Word VN = 0;
	Word X = 0;
	Word D0 = 0;
	Word HN = 0;
	Word HP = 0;


	i = 0;

	int err = 0;

	Word err_mask = (Word)1;

	int i_bd = (errthold << 1);


	int last_high = (errthold << 1);

	int t_length_1 = t_length - 1;

	while (i<t_length_1)
	{
		X = Peq[(uint8_t)text[i]] | VN;

		D0 = ((VP + (X&VP)) ^ VP) | X;

		HN = VP&D0;
		HP = VN | ~(VP | D0);

		X = D0 >> 1;
		VN = X&HP;
		VP = HN | ~(X | HP);

		if (!(D0&err_mask))
		{
			++err;
			if ((err - last_high)>errthold)
            {
                return (*return_t_end);
            }
		}
		get_error(i + 1, errthold, err, VP, VN, &line_error, &return_site);
		if(line_error != (unsigned int)-1)
		{
			(*return_t_end) = i;
			(*return_p_end) = return_site;
			(*return_err) = line_error;
		}

		Peq[(uint8_t)'A'] = Peq[(uint8_t)'A'] >> 1;
		Peq[(uint8_t)'C'] = Peq[(uint8_t)'C'] >> 1;
		Peq[(uint8_t)'G'] = Peq[(uint8_t)'G'] >> 1;
		Peq[(uint8_t)'T'] = Peq[(uint8_t)'T'] >> 1;


		++i;
		++i_bd;
		Peq[(uint8_t)pattern[i_bd]] = Peq[(uint8_t)pattern[i_bd]] | Mask;
	}





	X = Peq[(uint8_t)text[i]] | VN;
	D0 = ((VP + (X&VP)) ^ VP) | X;
	HN = VP&D0;
	HP = VN | ~(VP | D0);
	X = D0 >> 1;
	VN = X&HP;
	VP = HN | ~(X | HP);
	if (!(D0&err_mask))
	{
		++err;
		if ((err - last_high)>errthold)
		{
			return (*return_t_end);
		}
	}
	///i = t_length - 1
	get_error(i + 1, errthold, err, VP, VN, &line_error, &return_site);
	if(line_error != (unsigned int)-1)
	{
		(*return_t_end) = i;
		(*return_p_end) = return_site;
		(*return_err) = line_error;
	}

	return (*return_t_end);
}

inline int Reserve_Banded_BPM_Extension_REV
(char *pattern, int p_length, char *text, int t_length, unsigned short errthold, 
unsigned int* return_err, int* return_p_end, int* return_t_end)
{
    (*return_err) = (unsigned int)-1;
    (*return_p_end) = -1;
    (*return_t_end) = -1;

    Word Peq[256];

    unsigned int line_error = (unsigned int)-1;
    int return_site;
    int band_length = (errthold << 1) + 1;
    int i = 0;
    Word tmp_Peq_1 = (Word)1;

    Peq[(uint8_t)'A'] = (Word)0;
    Peq[(uint8_t)'T'] = (Word)0;
    Peq[(uint8_t)'G'] = (Word)0;
    Peq[(uint8_t)'C'] = (Word)0;


    Word Peq_A;
    Word Peq_T;
    Word Peq_C;
    Word Peq_G;

    ///band_length = 2k + 1
    for (i = 0; i<band_length; i++)
    {
        Peq[(uint8_t)pattern[p_length-i-1]] = Peq[(uint8_t)pattern[p_length-i-1]] | tmp_Peq_1;
        tmp_Peq_1 = tmp_Peq_1 << 1;
    }

    Peq_A = Peq[(uint8_t)'A'];
    Peq_C = Peq[(uint8_t)'C'];
    Peq_T = Peq[(uint8_t)'T'];
    Peq_G = Peq[(uint8_t)'G'];


    memset(Peq, 0, sizeof(Word)* 256);


    Peq[(uint8_t)'A'] = Peq_A;
    Peq[(uint8_t)'C'] = Peq_C;
    Peq[(uint8_t)'T'] = Peq_T;
    Peq[(uint8_t)'G'] = Peq_G;


    

    Word Mask = ((Word)1 << (errthold << 1));

    Word VP = 0;
    Word VN = 0;
    Word X = 0;
    Word D0 = 0;
    Word HN = 0;
    Word HP = 0;


    i = 0;

    int err = 0;

    Word err_mask = (Word)1;

    int i_bd = (errthold << 1);


    int last_high = (errthold << 1);

    int t_length_1 = t_length - 1;

    while (i<t_length_1)
    {
        X = Peq[(uint8_t)text[t_length-i-1]] | VN;

        D0 = ((VP + (X&VP)) ^ VP) | X;

        HN = VP&D0;
        HP = VN | ~(VP | D0);

        X = D0 >> 1;
        VN = X&HP;
        VP = HN | ~(X | HP);

        if (!(D0&err_mask))
        {
            ++err;
            if ((err - last_high)>errthold)
            {
                return (*return_t_end);
            }
        }
        get_error(i + 1, errthold, err, VP, VN, &line_error, &return_site);
        if(line_error != (unsigned int)-1)
        {
            (*return_t_end) = t_length-i-1;
            (*return_p_end) = p_length-return_site-1;
            (*return_err) = line_error;
        }

        Peq[(uint8_t)'A'] = Peq[(uint8_t)'A'] >> 1;
        Peq[(uint8_t)'C'] = Peq[(uint8_t)'C'] >> 1;
        Peq[(uint8_t)'G'] = Peq[(uint8_t)'G'] >> 1;
        Peq[(uint8_t)'T'] = Peq[(uint8_t)'T'] >> 1;


        ++i;
        ++i_bd;
        Peq[(uint8_t)pattern[p_length-i_bd-1]] = Peq[(uint8_t)pattern[p_length-i_bd-1]] | Mask;
    }





    X = Peq[(uint8_t)text[t_length-i-1]] | VN;
    D0 = ((VP + (X&VP)) ^ VP) | X;
    HN = VP&D0;
    HP = VN | ~(VP | D0);
    X = D0 >> 1;
    VN = X&HP;
    VP = HN | ~(X | HP);
    if (!(D0&err_mask))
    {
        ++err;
        if ((err - last_high)>errthold)
        {
            return (*return_t_end);
        }
    }
    ///i = t_length - 1
    get_error(i + 1, errthold, err, VP, VN, &line_error, &return_site);
    if(line_error != (unsigned int)-1)
    {
        (*return_t_end) = t_length-i-1;
        (*return_p_end) = p_length-return_site-1;
        (*return_err) = line_error;
    }

    return (*return_t_end);
}

inline void reverse_string(char* str, int strLen)
{
	int i, Len;
	char k;
	Len = strLen / 2;
	for (i = 0; i < Len; i++)
	{
		k = str[i];
		str[i] = str[strLen - i - 1];
		str[strLen - i - 1] = k;
	}
}

inline int alignment_extension(char *pattern, int p_length, char *text, int t_length, 
unsigned short errthold, int direction, unsigned int* return_err, int* return_p_end, 
int* return_t_end, int* return_aligned_t_len)
{
	(*return_aligned_t_len) = 0;
	
	if(direction == 0)
	{
		Reserve_Banded_BPM_Extension(pattern, p_length, text, t_length, errthold, return_err, 
		return_p_end, return_t_end);
		if((*return_p_end) != -1 && (*return_t_end) != -1)
		{
			(*return_aligned_t_len) = (*return_t_end) + 1;
			return 1;
		}
		else
		{
			return -1;
		}
		
	}
	else
	{
		reverse_string(pattern, p_length);
		reverse_string(text, t_length);

		Reserve_Banded_BPM_Extension(pattern, p_length, text, t_length, errthold, return_err, 
		return_p_end, return_t_end);

		reverse_string(pattern, p_length);
		reverse_string(text, t_length);

		if((*return_p_end) != -1 && (*return_t_end) != -1)
		{
			(*return_aligned_t_len) = (*return_t_end) + 1;
			(*return_p_end) = p_length - (*return_p_end);
			(*return_t_end) = t_length - (*return_t_end);
			return 1;
		}
		else
		{
			return -1;
		}
	}
}

inline void print_bit(Word z, int64_t w, const char *cmd)
{
	int64_t k;//, w = (sizeof(Word)<<3);
	fprintf(stderr, "%s\t", cmd);
	for (k = w-1; k >= 0; k--) fprintf(stderr, "%llu", (z>>k)&(1ULL));
	fprintf(stderr, "\n");
}

inline int32_t ed_band_cal_global(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre)
{
	if((pn > tn + thre) || (tn > pn + thre)) return INT32_MAX;
	if((pn < thre + 1) || (tn < thre + 1)) return INT32_MAX;
	Word Peq[5] = {0}, mm, VP = 0, VN = 0, X = 0, D0 = 0, HN = 0, HP = 0; 
	int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd = thre+1, i_bd = thre;
	// fprintf(stderr, "\n[M::%s::]\n", __func__);
    for (i = 0, mm = (((Word)1)<<thre); i < bd; i++) {
        Peq[seq_nt4_table[(uint8_t)pstr[i]]] |= mm; mm <<= 1;
    }
	Peq[4] = 0;
	err = thre;
	VN = (((Word)1)<<(thre))-1; 
	VP = (((Word)1)<<((thre<<1)+1))-1; VP ^= VN;

	// print_bit(Peq[0], (thre<<1)+1, "Peq[A]");
	// print_bit(Peq[1], (thre<<1)+1, "Peq[C]");
	// print_bit(Peq[2], (thre<<1)+1, "Peq[G]");
	// print_bit(Peq[3], (thre<<1)+1, "Peq[T]");
	// print_bit(Peq[4], (thre<<1)+1);
	// print_bit(VN, (thre<<1)+1, "VN");
	// print_bit(VP, (thre<<1)+1, "VP");

	///should make Peq[4] = 0 if N is always an error
	i = 0; mm = ((Word)1 << (thre<<1));///for the incoming char/last char

    while (i < tn0) {
        X = Peq[seq_nt4_table[(uint8_t)tstr[i]]] | VN;

        D0 = ((VP + (X&VP)) ^ VP) | X;

        HN = VP&D0;
        HP = VN | ~(VP | D0);

        X = D0 >> 1;
        VN = X&HP;
        VP = HN | ~(X | HP);
		// fprintf(stderr, "\n[M::%s::i->%d]\n", __func__, i);
		// print_bit(VN, (thre<<1)+1, "VN");
		// print_bit(VP, (thre<<1)+1, "VP");
		// print_bit(HN, (thre<<1)+1, "HN");
		// print_bit(HP, (thre<<1)+1, "HP");
		// print_bit(D0, (thre<<1)+1, "D0");

        if (!(D0&(1ULL))) {
            ++err;
            if (err>cut) return INT32_MAX;
        }
		// fprintf(stderr, "[M::%s::i->%d] err->%d\n", __func__, i, err);

		Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1; ///Peq[4] >>= 1;    

        ++i; ++i_bd;
        if(i_bd < pn) {
			Peq[seq_nt4_table[(uint8_t)pstr[i_bd]]] |= mm; Peq[4] = 0;
		}
        // if(i < pn) Peq[seq_nt4_table[(uint8_t)pstr[i]]] |= mm;
    }

    X = Peq[seq_nt4_table[(uint8_t)tstr[i]]] | VN;
    D0 = ((VP + (X&VP)) ^ VP) | X;
    HN = VP&D0;
    HP = VN | ~(VP | D0);
    X = D0 >> 1;
    VN = X&HP;
    VP = HN | ~(X | HP);
	// fprintf(stderr, "\n[M::%s::i->%d]\n", __func__, i);
	// print_bit(VN, (thre<<1)+1, "VN");
	// print_bit(VP, (thre<<1)+1, "VP");
	// print_bit(HN, (thre<<1)+1, "HN");
	// print_bit(HP, (thre<<1)+1, "HP");
	// print_bit(D0, (thre<<1)+1, "D0");
    if (!(D0&(1ULL))) {
		++err;
		if (err>cut) return INT32_MAX;
	}
	// fprintf(stderr, "[M::%s::i->%d] err->%d\n", __func__, i, err);

    int32_t site = tn - 1 - thre;///up bound
	for (cut = pn - 1, i = 0; site < cut; site++, i++) {
		// fprintf(stderr, "+[M::%s::site->%d] err->%d\n", __func__, site, err);
		err += ((VP >> i)&(1ULL)); err -= ((VN >> i)&(1ULL));
		// fprintf(stderr, "-[M::%s::site->%d] err->%d\n", __func__, site, err);
	}

	if (site == cut && err <= thre) return err;
    return INT32_MAX;
}

#define EAC_M 0
#define MIS_M 1 
#define MOR_YP 2
#define MOR_XT 3

inline void push_trace(asg16_v *res, uint16_t c, uint32_t len)
{
    uint16_t p; c <<= 14; 
	while (len >= (0x3fff)) {
		p = (c + (0x3fff)); kv_push(uint16_t, *res, p); len -= (0x3fff);
	}
	if(len) {
		p = (c + len); kv_push(uint16_t, *res, p);
	}
}

inline uint32_t pop_trace(asg16_v *res, uint32_t i, uint16_t *c, uint32_t *len)
{
	(*c) = (res->a[i]>>14); (*len) = (res->a[i]&(0x3fff));
	for (i++; (i < res->n) && ((*c) == (res->a[i]>>14)); i++) {
		(*len) += (res->a[i]&(0x3fff));
	}
	return i;
}

inline int32_t pop_trace_back(asg16_v *res, int32_t i, uint16_t *c, uint32_t *len)
{
	(*c) = (res->a[i]>>14); (*len) = (res->a[i]&(0x3fff));
	for (i--; (i >= 0) && ((*c) == (res->a[i]>>14)); i--) {
		(*len) += (res->a[i]&(0x3fff));
	}
	return i;
}

///511 -> 16 64-bits
// #define MAX_E 511
// #define MAX_L 2500

///511 -> 32 64-bits
#define MAX_CNS_E 1023
#define MAX_CNS_L 3072
#define FORCE_CNS_L 256


#define MAX_SIN_E 2047
#define MAX_SIN_L 10000
#define FORCE_SIN_L 512

typedef uint64_t w_sig;
#define bitw (6)
#define bitwbit (64)
#define bitz (63)
// typedef uint32_t w_sig;
// #define bitw (5)
// #define bitwbit (32)
// #define bitz (31)
// typedef uint16_t w_sig;
// #define bitw (4)
// #define bitwbit (16)
// #define bitz (15)
// typedef uint8_t w_sig;
// #define bitw (3)
// #define bitwbit (8)
// #define bitz (7)
typedef struct {w_sig *a;} w128_t;
typedef struct {size_t n, m; w_sig *a;} w64_trace_t;
typedef struct {
	int32_t done_cigar, cigar_n, done_path, path_n;
	int32_t ps, pe, pl, ts, te, tl;
	int32_t thre, err, nword, mword;
	uint32_t m, mm_thres; w_sig *a;
	w128_t Peq[5], mm, VP, VN, X, D0, HN, HP;
	asg16_v cigar; w64_trace_t path; 
} bit_extz_t;

#define is_align(exz) ((exz).err<=(exz).thre)
#define clear_align(exz) ((exz).err=INT32_MAX)

inline uint32_t cigar_check(char *pstr, char *tstr, bit_extz_t *ez)
{
	int32_t pi = ez->ps, ti = ez->ts, err = 0; uint32_t ci = 0, cl, k; uint16_t c;
	while (ci < ez->cigar.n) {
		ci = pop_trace(&(ez->cigar), ci, &c, &cl);
		// fprintf(stderr, "# %u = %u\n", c, cl); 
		if(c == 0) {
			for (k=0;(k<cl)&&(pstr[pi]==tstr[ti]);k++,pi++,ti++);
			if(k!=cl) {
				fprintf(stderr, "ERROR-d-0\n"); 
				return 0;
			}
		} else {
			err += cl;
			if(c == 1) {
				for (k=0;(k<cl)&&(pstr[pi]!=tstr[ti]);k++,pi++,ti++);
				if(k!=cl) {
					fprintf(stderr, "ERROR-d-1\n"); 
					return 0;
				}
			} else if(c == 2) {///more p
				pi+=cl;
			} else if(c == 3) {
				ti+=cl;
			}
		}
	}
	if(err != ez->err) {
		fprintf(stderr, "ERROR-err\n"); 
		return 0;
	}
	return 1;
}

inline int32_t dbg_ext_err(int32_t i, int32_t thre, int32_t err, int32_t pe, w128_t *VP, w128_t *VN)
{
	int32_t poff = i-thre, tmp_e = INT32_MAX, k, bd, k_bd;
	if((poff) + (thre<<1) >= (pe)) { 
		tmp_e = err;
		for ((k) = 0; (poff) < (pe); (poff)++) {
			bd = (k>>bitw); k_bd = (k&bitz);
			(tmp_e) += ((VP->a[bd]>>k_bd)&((w_sig)1)); 
			(tmp_e) -= ((VN->a[bd]>>k_bd)&((w_sig)1));
			(k)++;
		}
	}
	return tmp_e;
}


inline void print_bits(w_sig *az, int64_t w, const char *cmd)
{
	int64_t k, m, s = (sizeof(*az)<<3), sw = (w/s) + (!!(w%s)), ks;
	fprintf(stderr, "%s\t", cmd);
	for (m = sw - 1, k = w-1; m >= 0 && k >= 0; m--) {
		for (ks = k%s; ks >= 0 && k >= 0; ks--, k--) fprintf(stderr, "%llu", (az[m]>>ks)&(1ULL));
	}
	fprintf(stderr, "\n");
}

#define prt_bit_extz_t(ez, w) do {	\
	print_bits((ez).Peq[0].a, w, "Peq[0]");\
	print_bits((ez).Peq[1].a, w, "Peq[1]");\
	print_bits((ez).Peq[2].a, w, "Peq[2]");\
	print_bits((ez).Peq[3].a, w, "Peq[3]");\
	print_bits((ez).Peq[4].a, w, "Peq[4]");\
	print_bits((ez).VP.a, w, "VP");\
	print_bits((ez).VN.a, w, "VN");\
	print_bits((ez).X.a, w, "X");\
	print_bits((ez).D0.a, w, "D0");\
	print_bits((ez).HN.a, w, "HN");\
	print_bits((ez).HP.a, w, "HP");\
} while (0)		
	


#define resize_bit_extz_t(ex, thres) do {	\
		if(((int32_t)(thres)) > ((int32_t)(ex).mm_thres)) {\
			(ex).mm_thres = (((thres)<<1)+1);\
			(ex).mm_thres = (((ex).mm_thres>>bitw)+(!!((ex).mm_thres&bitz)))*12;\
			if((ex).mm_thres > (ex).m) {\
				(ex).a = (w_sig *)realloc((ex).a, (ex).mm_thres * sizeof(*((ex).a)));\
				(ex).mm_thres/=12; (ex).m=0;\
				(ex).Peq[0].a = (ex).a; (ex).m+=(ex).mm_thres;\
				(ex).Peq[1].a = (ex).a+(ex).m; (ex).m+=(ex).mm_thres;\
				(ex).Peq[2].a = (ex).a+(ex).m; (ex).m+=(ex).mm_thres;\
				(ex).Peq[3].a = (ex).a+(ex).m; (ex).m+=(ex).mm_thres;\
				(ex).Peq[4].a = (ex).a+(ex).m; (ex).m+=(ex).mm_thres;\
				(ex).mm.a = (ex).a+(ex).m; (ex).m+=(ex).mm_thres;\
				(ex).VP.a = (ex).a+(ex).m; (ex).m+=(ex).mm_thres;\
				(ex).VN.a = (ex).a+(ex).m; (ex).m+=(ex).mm_thres;\
				(ex).X.a = (ex).a+(ex).m; (ex).m+=(ex).mm_thres;\
				(ex).D0.a = (ex).a+(ex).m; (ex).m+=(ex).mm_thres;\
				(ex).HN.a = (ex).a+(ex).m; (ex).m+=(ex).mm_thres;\
				(ex).HP.a = (ex).a+(ex).m; (ex).m+=(ex).mm_thres;\
				(ex).mword = (ex).mm_thres;\
			}\
			(ex).mm_thres = (thres);\
		}\
	} while (0)	

inline void init_bit_extz_t(bit_extz_t *ex, uint64_t thres) {
	memset(ex, 0, sizeof((*ex)));
	///(bitwbit<<2) >= (((thres)<<1)+1);->at least 4 cells for each w128_t
	if(((uint64_t)thres)<(((((uint64_t)bitwbit)<<2)-1)>>1)) thres=(((bitwbit<<2)-1)>>1);
	resize_bit_extz_t((*ex), (thres));
}

inline void destroy_bit_extz_t(bit_extz_t *ex) {
	free((*ex).a); free((*ex).path.a); free((*ex).cigar.a);
}


inline void gen_trace(bit_extz_t *ez, int32_t ptrim, int32_t reverse)///ptrim = thre for global and extension; = 0 for semi
{
	if(ez->err > ez->thre) return;
	ez->cigar.n = 0;
	int32_t V, H, D, min, cur, tn = (ez->te+1-ez->ts), pn = tn + (ez->thre<<1), bd = (ez->thre<<1)+1;
	int32_t bs = ez->path.n/tn, bbs = bs/5, poff = ez->pe, sft = bd - (pn - ez->pe - ptrim);
	int32_t i = tn, low = bd-1, wi, wm, d = 0, pd = -1, pdn = 0; cur = ez->err;
	w_sig *D0, *VP, *VN, *HP, *HN;

	// fprintf(stderr, "\n[M::%s::] bs::%d, bbs::%d, cur::%d, poff::%d, sft::%d\n", __func__, bs, bbs, cur, poff, sft);

	while (i > 0 && cur > 0) {
		D0 = ez->path.a + ((i-1)*bs); VP = D0 + bbs; VN = VP + bbs; HP = VN + bbs; HN = HP + bbs;

		wi = (sft>>bitw); wm = sft&bitz;
		D = cur - ((~(D0[wi]>>wm))&(1ULL)); d = 0; min = D;
		H = V = INT32_MAX;
		if(sft!=low) {
            H = cur + ((HN[wi]>> wm)&(1ULL)) - ((HP[wi]>> wm)&(1ULL)); 
            if ((H+1) == cur && H <= min) {//prefer indels
                min = H; d = 3;
            }
		}
		if(sft!=0) {
			wi = ((sft-1)>>bitw); wm = (sft-1)&bitz;
			V = cur + ((VN[wi]>> wm)&(1ULL)) - ((VP[wi]>> wm)&(1ULL)); 
			if ((V+1) == cur && V <= min) {//prefer indels
				min = V; d = 2;
			}
		}
		// fprintf(stderr, "[M::%s::] cur::%d, D::%d, V::%d, H::%d, poff::%d, toff::%d, d::%d\n", __func__, cur, D, V, H, poff, i, d);
		if(d == 0) {
			if(D != cur) d = 1;
			i--; poff--;
		} else if(d == 2) {//more pstr
			sft--; poff--;
		} else if(d == 3) {///more tstr
			i--; sft++;
		}

		if(d == pd) {
			pdn++;
		} else {
			if(pdn > 0) push_trace(&(ez->cigar), pd, pdn);
			pd = d; pdn = 1;
		}
		cur = min;
	}

	if (i > 0) {
        d = 0; poff -= i;
		if(d == pd) {
			pdn += i;
		} else {
			if(pdn > 0) push_trace(&(ez->cigar), pd, pdn);
			pd = d; pdn = i;
		}
    }

	poff++;
	// fprintf(stderr, "[M::%s::] poff::%d, ez->ps::%d\n", __func__, poff, ez->ps);
	if(ez->ps < 0 || ez->ps >= ez->pl) {//ps is unavailable
		ez->ps = poff;
	} else if(poff > ez->ps){
		d = 2; i = poff - ez->ps;
		if(d == pd) {
			pdn += i;
		} else {
			if(pdn > 0) push_trace(&(ez->cigar), pd, pdn);
			pd = d; pdn = i;
		}
	}
	if(pdn > 0) push_trace(&(ez->cigar), pd, pdn);

	if(reverse) {
		uint16_t t; pdn = ez->cigar.n>>1; 
		for (i = 0; i < pdn; i++) {
			t = ez->cigar.a[i]; 
			ez->cigar.a[i] = ez->cigar.a[ez->cigar.n-i-1]; 
			ez->cigar.a[ez->cigar.n-i-1] = t;
		}
	}
}



#define init_base_ed(ez, thre, pn, tn) {\
	(ez).thre = (thre), (ez).err = INT32_MAX, (ez).pl = pn, (ez).tl = tn;\
	(ez).done_cigar = (ez).done_path = (ez).cigar_n = (ez).path_n = 0;\
}

#define w_bit(x, b) ((x).a[((b)>>bitw)]|=(((w_sig)1)<<((b)&bitz)))

#define w_get_bit(x, b) (((x).a[((b)>>bitw)]>>((b)&bitz))&((w_sig)1))

/***********************2 words***********************/
#define w_128_clear(x) ((x).a[0]=(x).a[1]=0)

#define w_128_word (2)

#define w_128_self_not(x) ((x).a[0]=~(x).a[0], \
									(x).a[1]=~(x).a[1])

#define w_128_self_or(x, y) ((x).a[0]|=(y).a[0], \
									(x).a[1]|=(y).a[1])

#define w_128_or(r, x, y) ((r).a[0] = (x).a[0]|(y).a[0], \
									(r).a[1] = (x).a[1]|(y).a[1])

#define w_128_and(r, x, y) ((r).a[0] = (x).a[0]&(y).a[0], \
									(r).a[1] = (x).a[1]&(y).a[1])

#define w_128_self_xor(x, y) ((x).a[0]^=(y).a[0], \
									(x).a[1]^=(y).a[1])

#define w_128_self_lsft_1(x) ((x).a[1] = ((x).a[1]<<1)|((x).a[0]>>bitz), \
																		(x).a[0] <<= 1)

#define w_128_self_rsft_1(x) ((x).a[0] = ((x).a[0]>>1)|((x).a[1]<<bitz), \
																		(x).a[1] >>= 1)

#define w_128_rsft_1(x, y) ((x).a[0] = ((y).a[0]>>1)|((y).a[1]<<bitz), \
																		(x).a[1] = (y).a[1]>>1)
#define w_128_self_add(x, y, c) ((x).a[0]+=(y).a[0], \
								(x).a[1]+=(y).a[1]+((x).a[0]<(y).a[0]))

#define w_128_set_bit_lsub(x, l) do {	\
		(x).a[0] = (w_sig)-1, (x).a[1] = 0;	\
		if((l) <= bitwbit) (x).a[0] = (((w_sig)1)<<(l))-1; \
		else (x).a[1] = (((w_sig)1)<<((l)-bitwbit))-1;\
	} while (0)												\

#define w_128_copy(x, y) ((x).a[0]=(y).a[0], (x).a[1]=(y).a[1])

/***********************3 words***********************/
#define w_192_clear(x) ((x).a[0]=(x).a[1]=(x).a[2]=0)

#define w_192_word (3)

#define w_192_self_not(x) ((x).a[0]=~(x).a[0],\
									(x).a[1]=~(x).a[1],\
												(x).a[2]=~(x).a[2])

#define w_192_self_or(x, y) ((x).a[0]|=(y).a[0],\
									(x).a[1]|=(y).a[1],\
												(x).a[2]|=(y).a[2])

#define w_192_or(r, x, y) ((r).a[0]=(x).a[0]|(y).a[0],\
									(r).a[1] = (x).a[1]|(y).a[1],\
												(r).a[2] = (x).a[2]|(y).a[2])

#define w_192_and(r, x, y) ((r).a[0] = (x).a[0]&(y).a[0], \
									(r).a[1] = (x).a[1]&(y).a[1], \
											(r).a[2] = (x).a[2]&(y).a[2])

#define w_192_self_xor(x, y) ((x).a[0]^=(y).a[0], \
									(x).a[1]^=(y).a[1], \
											(x).a[2]^=(y).a[2])

#define w_192_self_lsft_1(x) ((x).a[2] = ((x).a[2]<<1)|((x).a[1]>>bitz), \
										(x).a[1] = ((x).a[1]<<1)|((x).a[0]>>bitz), \
																		(x).a[0] <<= 1)

#define w_192_self_rsft_1(x) ((x).a[0] = ((x).a[0]>>1)|((x).a[1]<<bitz), \
										(x).a[1] = ((x).a[1]>>1)|((x).a[2]<<bitz), \
																		(x).a[2] >>= 1)

#define w_192_rsft_1(x, y) ((x).a[0] = ((y).a[0]>>1)|((y).a[1]<<bitz), \
										(x).a[1] = ((y).a[1]>>1)|((y).a[2]<<bitz), \
													(x).a[2] = (y).a[2]>>1)

#define w_192_self_add(x, y, c) ((x).a[0]+=(y).a[0], c=((x).a[0]<(y).a[0]),\
								(x).a[1]+=c, c=((x).a[1]<c), (x).a[1]+=(y).a[1], c|=((x).a[1]<(y).a[1]),\
								(x).a[2]+=c, (x).a[2]+=(y).a[2])

#define w_192_set_bit_lsub(x, l) do {	\
		(x).a[0] = (x).a[1] = (x).a[2] = 0;	\
		if((l)>>bitw) memset((x).a, -1, sizeof(*((x).a))*((l)>>bitw));\
		if((l)&bitz) (x).a[(l)>>bitw] = (((w_sig)1)<<((l)&bitz))-1; \
	} while (0)												\

#define w_192_copy(x, y) ((x).a[0]=(y).a[0], (x).a[1]=(y).a[1], (x).a[2]=(y).a[2])

/***********************4 words***********************/
#define w_256_clear(x) ((x).a[0]=(x).a[1]=(x).a[2]=(x).a[3]=0)

#define w_256_word (4)

#define w_256_self_not(x) ((x).a[0]=~(x).a[0],\
									(x).a[1]=~(x).a[1],\
												(x).a[2]=~(x).a[2],\
													(x).a[3]=~(x).a[3])

#define w_256_self_or(x, y) ((x).a[0]|=(y).a[0],\
									(x).a[1]|=(y).a[1],\
												(x).a[2]|=(y).a[2],\
													(x).a[3]|=(y).a[3])

#define w_256_or(r, x, y) ((r).a[0]=(x).a[0]|(y).a[0],\
									(r).a[1] = (x).a[1]|(y).a[1],\
												(r).a[2] = (x).a[2]|(y).a[2],\
													(r).a[3] = (x).a[3]|(y).a[3])

#define w_256_and(r, x, y) ((r).a[0] = (x).a[0]&(y).a[0], \
									(r).a[1] = (x).a[1]&(y).a[1], \
											(r).a[2] = (x).a[2]&(y).a[2], \
												(r).a[3] = (x).a[3]&(y).a[3])

#define w_256_self_xor(x, y) ((x).a[0]^=(y).a[0], \
									(x).a[1]^=(y).a[1], \
											(x).a[2]^=(y).a[2],\
												(x).a[3]^=(y).a[3])

#define w_256_self_lsft_1(x) ((x).a[3] = ((x).a[3]<<1)|((x).a[2]>>bitz), \
								(x).a[2] = ((x).a[2]<<1)|((x).a[1]>>bitz), \
										(x).a[1] = ((x).a[1]<<1)|((x).a[0]>>bitz), \
																		(x).a[0] <<= 1)

#define w_256_self_rsft_1(x) ((x).a[0] = ((x).a[0]>>1)|((x).a[1]<<bitz), \
										(x).a[1] = ((x).a[1]>>1)|((x).a[2]<<bitz), \
											(x).a[2] = ((x).a[2]>>1)|((x).a[3]<<bitz), \
																			(x).a[3] >>= 1)

#define w_256_rsft_1(x, y) ((x).a[0] = ((y).a[0]>>1)|((y).a[1]<<bitz), \
								(x).a[1] = ((y).a[1]>>1)|((y).a[2]<<bitz), \
									(x).a[2] = ((y).a[2]>>1)|((y).a[3]<<bitz), \
											(x).a[3] = (y).a[3]>>1)

#define w_256_self_add(x, y, c) ((x).a[0]+=(y).a[0], c=((x).a[0]<(y).a[0]),\
								(x).a[1]+=c, c=((x).a[1]<c), (x).a[1]+=(y).a[1], c|=((x).a[1]<(y).a[1]),\
								(x).a[2]+=c, c=((x).a[2]<c), (x).a[2]+=(y).a[2], c|=((x).a[2]<(y).a[2]),\
								(x).a[3]+=c, (x).a[3]+=(y).a[3])

#define w_256_set_bit_lsub(x, l) do {	\
		(x).a[0] = (x).a[1] = (x).a[2] = (x).a[3] = 0;	\
		if((l)>>bitw) memset((x).a, -1, sizeof(*((x).a))*((l)>>bitw));\
		if((l)&bitz) (x).a[(l)>>bitw] = (((w_sig)1)<<((l)&bitz))-1; \
	} while (0)	

#define w_256_copy(x, y) ((x).a[0]=(y).a[0], (x).a[1]=(y).a[1], (x).a[2]=(y).a[2], (x).a[3]=(y).a[3])

#define cmp_Word(des, src, dd, ws) {\
	for ((dd)=0;(dd)<64;dd+=bitwbit){\
		if((((des)>>dd)&((w_sig)-1))!=(src).a[dd>>bitw]) break;\
	}\
	if((dd)<64) {\
		print_bit((des), (ws), "des");\
		print_bits((src).a, (ws), "src");\
		exit(0);\
	}\
}

#define dump_Word(des, src, dd, ws) {\
	for ((dd)=bitwbit,(des)=((Word)(src).a[0]);(dd)<64;dd+=bitwbit) {\
		(des) |= ((Word)(src).a[dd>>bitw])<<dd;\
	}\
	cmp_Word((des), (src), (dd), (ws));\
}

#define prt_address(Peq, VP, VN, X, D0, HN, HP) {\
	fprintf(stderr, "Peq[0].a::%p, Peq[1].a::%p, Peq[2].a::%p, Peq[3].a::%p, Peq[4].a::%p, VP.a::%p, VN.a::%p, X.a::%p, D0.a::%p, HN.a::%p, HP.a::%p\n",\
	Peq[0].a, Peq[1].a, Peq[2].a, Peq[3].a, Peq[4].a, VP.a, VN.a, X.a, D0.a, HN.a, HP.a);\
}

#define prt_vector(Peq, VP, VN, X, D0, HN, HP, ws) {\
	print_bit((Peq)[0], (ws), "Peq[0]");\
	print_bit((Peq)[1], (ws), "Peq[1]");\
	print_bit((Peq)[2], (ws), "Peq[2]");\
	print_bit((Peq)[3], (ws), "Peq[3]");\
	print_bit((Peq)[4], (ws), "Peq[4]");\
	print_bit((VP), (ws), "VP");\
	print_bit((VN), (ws), "VN");\
	print_bit((X), (ws), "X");\
	print_bit((D0), (ws), "D0");\
	print_bit((HN), (ws), "HN");\
	print_bit((HP), (ws), "HP");\
}

#define ed_core_debug(sf, Peq, VP, VN, X, D0, HN, HP, z, c, ws) {	\
		Word Peq_d[5], VP_d, VN_d, X_d, D0_d, HN_d, HP_d, dd;\
		dump_Word(Peq_d[0], Peq[0], dd, ws);\
		dump_Word(Peq_d[1], Peq[1], dd, ws);\
		dump_Word(Peq_d[2], Peq[2], dd, ws);\
		dump_Word(Peq_d[3], Peq[3], dd, ws);\
		dump_Word(Peq_d[4], Peq[4], dd, ws);\
		dump_Word(VP_d, VP, dd, ws);\
		dump_Word(VN_d, VN, dd, ws);\
		dump_Word(X_d, X, dd, ws);\
		dump_Word(D0_d, D0, dd, ws);\
		dump_Word(HN_d, HN, dd, ws);\
		dump_Word(HP_d, HP, dd, ws);\
		/**X = Peq[seq_nt4_table[(uint8_t)tstr[i]]] | VN;**/\
		prt_address(Peq, VP, VN, X, D0, HN, HP);\
		(c) = seq_nt4_table[(z)];\
		print_bit(Peq_d[(c)], 64, "raw-Peq[]");print_bits(Peq[(c)].a, 64, "bit-Peq[]");\
		print_bit(VN_d, 64, "raw-VN");print_bits((VN).a, 64, "bit-VN");\
		w##sf##or((X), (Peq)[(c)], (VN));\
		print_bits((X).a, 64, "+bit-X");\
		X_d = Peq_d[(c)]|VN_d; fprintf(stderr, "X=Peq[]|VN\n"); \
		print_bit(X_d, 64, "raw-X");print_bits((X).a, 64, "-bit-X");\
		cmp_Word(X_d, (X), dd, ws);\
		/**D0 = ((VP + (X&VP)) ^ VP) | X;**/\
		w##sf##and((D0), (X), (VP));\
		prt_address(Peq, VP, VN, X, D0, HN, HP);\
		D0_d = X_d&VP_d; fprintf(stderr, "D0=X&VP\n"); cmp_Word(D0_d, (D0), dd, ws);\
		w##sf##self_add((D0), (VP), (c));\
		prt_address(Peq, VP, VN, X, D0, HN, HP);\
		D0_d += VP_d; fprintf(stderr, "D0+=VP\n"); cmp_Word(D0_d, (D0), dd, ws);\
		w##sf##self_xor((D0), (VP));\
		prt_address(Peq, VP, VN, X, D0, HN, HP);\
		D0_d ^= VP_d; fprintf(stderr, "D0^=VP\n"); cmp_Word(D0_d, (D0), dd, ws);\
		/**print_bit(D0_d, (ws), "raw-D0"); print_bit(X_d, (ws), "raw-X");**/\
		/**print_bits((D0).a, (ws), "bit-D0"); print_bits((X).a, (ws), "bit-X");**/\
		w##sf##self_or((D0), (X));\
		prt_address(Peq, VP, VN, X, D0, HN, HP);\
		D0_d |= X_d; fprintf(stderr, "D0|=X\n"); cmp_Word(D0_d, (D0), dd, ws);\
        /**HN = VP&D0;**/\
		w##sf##and((HN), (VP), (D0));\
		prt_address(Peq, VP, VN, X, D0, HN, HP);\
		HN_d = VP_d&D0_d; fprintf(stderr, "HN=VP&D0\n"); cmp_Word(HN_d, (HN), dd, ws);\
        /**HP = VN | ~(VP | D0);**/\
		w##sf##or((HP), (VP), (D0));\
		prt_address(Peq, VP, VN, X, D0, HN, HP);\
		HP_d = VP_d|D0_d; fprintf(stderr, "HP=VP|D0\n"); cmp_Word(HP_d, (HP), dd, ws);\
		w##sf##self_not((HP));\
		prt_address(Peq, VP, VN, X, D0, HN, HP);\
		HP_d = ~HP_d; fprintf(stderr, "HP~=HP\n"); cmp_Word(HP_d, (HP), dd, ws);\
		w##sf##self_or((HP), (VN));\
		prt_address(Peq, VP, VN, X, D0, HN, HP);\
		HP_d |= VN_d; fprintf(stderr, "HP|=VN\n"); cmp_Word(HP_d, (HP), dd, ws);\
        /**X = D0 >> 1;**/\
		w##sf##copy((X), (D0));\
		prt_address(Peq, VP, VN, X, D0, HN, HP);\
		X_d = D0_d; fprintf(stderr, "X=D0\n"); cmp_Word(X_d, (X), dd, ws);\
		w##sf##self_rsft_1((X)); \
		prt_address(Peq, VP, VN, X, D0, HN, HP);\
		X_d >>= 1; fprintf(stderr, "X>>=1\n"); cmp_Word(X_d, (X), dd, ws);\
        /**VN = X&HP;**/\
		w##sf##and((VN), (X), (HP));\
		prt_address(Peq, VP, VN, X, D0, HN, HP);\
		VN_d = X_d&HP_d; fprintf(stderr, "VN=X&HP\n"); cmp_Word(VN_d, (VN), dd, ws);\
        /**VP = HN | ~(X | HP);**/\
		w##sf##or((VP), (X), (HP));\
		prt_address(Peq, VP, VN, X, D0, HN, HP);\
		VP_d = X_d|HP_d; fprintf(stderr, "VP=X|HP\n"); cmp_Word(VP_d, (VP), dd, ws);\
		w##sf##self_not((VP));\
		prt_address(Peq, VP, VN, X, D0, HN, HP);\
		VP_d = ~VP_d; fprintf(stderr, "VP=~VP\n"); cmp_Word(VP_d, (VP), dd, ws);\
		w##sf##self_or((VP), (HN));\
		prt_address(Peq, VP, VN, X, D0, HN, HP);\
		VP_d |= HN_d; fprintf(stderr, "VP|=HN\n"); cmp_Word(VP_d, (VP), dd, ws);\
}

#define ed_core(sf, Peq, VP, VN, X, D0, HN, HP, z, c) {	\
		/**X = Peq[seq_nt4_table[(uint8_t)tstr[i]]] | VN;**/\
		(c) = seq_nt4_table[(z)];\
		w##sf##or((X), (Peq)[(c)], (VN));\
		/**D0 = ((VP + (X&VP)) ^ VP) | X;**/\
		w##sf##and((D0), (X), (VP));\
		w##sf##self_add((D0), (VP), (c));\
		w##sf##self_xor((D0), (VP));\
		w##sf##self_or((D0), (X));\
        /**HN = VP&D0;**/\
		w##sf##and((HN), (VP), (D0));\
        /**HP = VN | ~(VP | D0);**/\
		w##sf##or((HP), (VP), (D0));\
		w##sf##self_not((HP));\
		w##sf##self_or((HP), (VN));\
        /**X = D0 >> 1;**/\
		/**w##sf##copy((X), (D0));w##sf##self_rsft_1((X));**/\
		w##sf##rsft_1((X), (D0));\
        /**VN = X&HP;**/\
		w##sf##and((VN), (X), (HP));\
        /**VP = HN | ~(X | HP);**/\
		w##sf##or((VP), (X), (HP));\
		w##sf##self_not((VP));\
		w##sf##self_or((VP), (HN));\
}

#define ed_init_core(i, bd, i_bd, S, Peq) {	\
	for ((i) = 0; (i) < (bd); (i)++) {\
		w_bit((Peq)[seq_nt4_table[(uint8_t)(S)[(i)]]], (i_bd));i_bd++;\
    }\
}

#define HA_ED_INIT(sf)\
inline void ed_band_cal_global_##sf##_w(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)\
{\
	init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = 0;\
	if((pn > tn + thre) || (tn > pn + thre)) return;\
	int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd, i_bd, Peq_i; w_sig c, Peq_m;\
	w_##sf##_clear(ez->Peq[0]); w_##sf##_clear(ez->Peq[1]); w_##sf##_clear(ez->Peq[2]); w_##sf##_clear(ez->Peq[3]); w_##sf##_clear(ez->Peq[4]);\
	\
	bd = thre+1; bd = ((bd<=pn)?bd:pn); i_bd = thre;\
	ed_init_core(i, bd, i_bd, pstr, ez->Peq);\
	bd = thre+1, i_bd = thre;\
	\
	w_##sf##_clear(ez->Peq[4]);\
	err = thre;\
	w_##sf##_set_bit_lsub(ez->VN, thre); /**VN = (((Word)1)<<(thre))-1;**/\ 
	w_##sf##_set_bit_lsub(ez->VP, (thre<<1)+1); /**VP = (((Word)1)<<((thre<<1)+1))-1;**/\
	w_##sf##_self_xor(ez->VP, ez->VN); /**VP ^= VN;**/\
	\
	/** print_bits(VP.a, (thre<<1)+1, "-VP");**/\
	/**should make Peq[4] = 0 if N is always an error**/\
	i = 0; \
	/**for the incoming char/last char**/\
	Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));\
    while (i < tn0) {\
		ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c);\
        if (!(ez->D0.a[0]&(1ULL))) {\
            ++err; if (err>cut) return;\
        }\
		/** Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;**/\
		w_##sf##_self_rsft_1(ez->Peq[0]); w_##sf##_self_rsft_1(ez->Peq[1]);\
		w_##sf##_self_rsft_1(ez->Peq[2]); w_##sf##_self_rsft_1(ez->Peq[3]);\
		\
		++i; ++i_bd; c = 4;\
        if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];\
        if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;\
    }\
	ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c);\
	if (!(ez->D0.a[0]&(1ULL))) {\
		++err; if (err>cut) return;\
	}\
	\
    int32_t site = tn - 1 - thre;/**up bound**/\
	for (cut = pn - 1, i = 0; site < cut; site++, i++) {\
		bd = (i>>bitw); i_bd = (i&bitz);\
        err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));\
        err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));\
	}\
	\
	if (site == cut && err <= thre) {\
		ez->err = err; \
		ez->pe = pn-1; ez->te = tn-1;\
	}\
    return;\
}\
inline void ed_band_cal_semi_##sf##_w(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)\
{\
	init_base_ed(*ez, thre, pn, tn); ez->ps = ez->pe = -1; ez->ts = 0; ez->te = tn-1;\
	int32_t bd, i, err = 0, i_bd, cut = thre+(thre<<1), tn0 = tn - 1, Peq_i; w_sig c, Peq_m;\
	if((pn > tn + cut) || (tn > pn + cut)) return;\
	\
	w_##sf##_clear(ez->Peq[0]);\
	w_##sf##_clear(ez->Peq[1]);\
	w_##sf##_clear(ez->Peq[2]);\
	w_##sf##_clear(ez->Peq[3]);\
	w_##sf##_clear(ez->Peq[4]);\
	w_##sf##_clear(ez->VP);\
	w_##sf##_clear(ez->VN);\
	\
	bd = (thre<<1)+1; bd = ((bd<=pn)?bd:pn); i_bd = 0;\
	ed_init_core(i, bd, i_bd, pstr, ez->Peq);\
	bd = (thre<<1)+1, i_bd = (thre<<1);\
	\
	w_##sf##_clear(ez->Peq[4]);\
	i = 0;\
	/**for the incoming char/last char; mm = ((Word)1 << (thre<<1))**/\
    Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));\
    while (i < tn0) {\
		ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c);\
		if (!(ez->D0.a[0]&(1ULL))) {\
            ++err; if (err>cut) return;\
        }\
		/** Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;**/\
		w_##sf##_self_rsft_1(ez->Peq[0]); w_##sf##_self_rsft_1(ez->Peq[1]);\
		w_##sf##_self_rsft_1(ez->Peq[2]); w_##sf##_self_rsft_1(ez->Peq[3]);\
        ++i; ++i_bd; c = 4;\
        if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];\
        if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;\
    }\
	ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c);\
	if (!(ez->D0.a[0]&(1ULL))) {\
		++err; if (err>cut) return;\
	}\
	\
    int32_t site = tn - 1;/**up bound**/\
    /**in most cases, ai = (thre<<1)**/\
    int32_t ai = pn - tn, uge = INT32_MAX;\
    if ((err <= thre) && (err <= ez->err)) {\
        ez->err = err; ez->pe = site;\
    }\
	\
    i = 0;\
    while (i < ai) {\
        bd = (i>>bitw); i_bd = (i&bitz);\
        err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));\
        err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));\
		++i;\
        if ((err <= thre) && (err <= ez->err)) {\
            ez->err = err; ez->pe = site + i;\
        }\
        if(i == thre) uge = err;\
    }\
	\
    if((uge <= thre) && (uge == ez->err)) ez->pe = site + thre;\
}\
/**require:: (pn >= tn - thre && pn <= tn + thre)**/\
inline void ed_band_cal_extension_##sf##_0_w(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)\
{\
	init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = 0; ez->pe = ez->te = -1;\
	if(pn > tn + thre) pn = tn + thre;\
	else if(tn > pn + thre) tn = pn + thre;\
	int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd, i_bd, k_bd, Peq_i; w_sig c, Peq_m;\
	int32_t poff, pe = pn-1, tmp_e = INT32_MAX, k;\
	\
	w_##sf##_clear(ez->Peq[0]);\
	w_##sf##_clear(ez->Peq[1]);\
	w_##sf##_clear(ez->Peq[2]);\
	w_##sf##_clear(ez->Peq[3]);\
	w_##sf##_clear(ez->Peq[4]);\
	\
	bd = thre+1; bd = ((bd<=pn)?bd:pn); i_bd = thre;\
	ed_init_core(i, bd, i_bd, pstr, ez->Peq);\
	bd = thre+1, i_bd = thre;\
	\
	w_##sf##_clear(ez->Peq[4]);\
	err = thre;\
	w_##sf##_set_bit_lsub(ez->VN, thre); /**VN = (((Word)1)<<(thre))-1; **/\
	w_##sf##_set_bit_lsub(ez->VP, (thre<<1)+1); /**VP = (((Word)1)<<((thre<<1)+1))-1;**/\
	w_##sf##_self_xor(ez->VP, ez->VN); /**VP ^= VN;**/\
	\
	/** print_bits(ez->VP.a, (thre<<1)+1, "-VP");**/\
	/**should make Peq[4] = 0 if N is always an error**/\
	i = 0; \
	/**for the incoming char/last char; mm = ((Word)1 << (thre<<1))**/\
    Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));\
    while (i < tn0) {\
		ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c);\
		if (!(ez->D0.a[0]&(1ULL))) {\
            ++err; if (err>cut) return;\
        }\
		\
        {\
            poff = i-thre; k = i+thre-pe;/**poff:[i-thre, i+thre]**/\
            if(k >= 0) {\
                if(tmp_e == INT32_MAX) {\
                    tmp_e = err;\
                    for ((k) = 0; (poff) < (pe); (poff)++, (k)++) {\
                        bd = (k>>bitw); k_bd = (k&bitz);\
                        (tmp_e) += (((*ez).VP.a[bd]>>k_bd)&((w_sig)1));\
                        (tmp_e) -= (((*ez).VN.a[bd]>>k_bd)&((w_sig)1));\
                    }\
                } else {\
                    k = (thre<<1) - k;\
                    if(k >= 0) {\
                        bd = (k>>bitw); k_bd = (k&bitz);\
                        (tmp_e) += (((*ez).HP.a[bd]>>k_bd)&((w_sig)1));\
                        (tmp_e) -= (((*ez).HN.a[bd]>>k_bd)&((w_sig)1));\
                    }\
                }\
                if((tmp_e) <= (*ez).thre && (tmp_e) < (*ez).err) {\
                    (*ez).err = tmp_e; (*ez).pe = pe; (*ez).te = i;\
                }\
            }\
            /**if(dbg_ext_err(i, thre, err, pe, &((*ez).VP), &((*ez).VN))!=tmp_e)**/\
        }\
		\
		/** Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;**/\
		w_##sf##_self_rsft_1(ez->Peq[0]); w_##sf##_self_rsft_1(ez->Peq[1]);\
		w_##sf##_self_rsft_1(ez->Peq[2]); w_##sf##_self_rsft_1(ez->Peq[3]);\
        ++i; ++i_bd; c = 4;\
        if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];\
        if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;\
    }\
	ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c);\
	if (!(ez->D0.a[0]&(1ULL))) {\
		++err; if (err>cut) return;\
	}\
	\
    int32_t site = tn - 1 - thre;/**up bound; site:[tn - 1 - thre, tn - 1 + thre]**/\
	for (cut = pn - 1, i = 0; site < cut; i++) {\
		bd = (i>>bitw); i_bd = (i&bitz);\
        err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));\
        err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));\
		site++;\
		if(err <= thre && err < ez->err) {\
			ez->err = err; ez->pe = site; ez->te = tn-1;\
		}\
	}\
	if(err <= thre && err < ez->err) {\
		ez->err = err; ez->pe = site; ez->te = tn-1;\
	}\
    return;\
}\
/**require:: (pn >= tn - thre && pn <= tn + thre)**/\
inline void ed_band_cal_extension_##sf##_1_w(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)\
{\
	init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = INT32_MAX; ez->pe = pn-1; ez->te = tn-1;\
	if(pn > tn + thre) pn = tn + thre;\
	else if(tn > pn + thre) tn = pn + thre;\
	int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd, i_bd, k_bd, Peq_i; w_sig c, Peq_m;\
	int32_t poff, pe = pn-1, tmp_e = INT32_MAX, k, pidx = ez->pe, tidx = ez->te;\
	\
	w_##sf##_clear(ez->Peq[0]);\
	w_##sf##_clear(ez->Peq[1]);\
	w_##sf##_clear(ez->Peq[2]);\
	w_##sf##_clear(ez->Peq[3]);\
	w_##sf##_clear(ez->Peq[4]);\
	\
	bd = thre+1; bd = ((bd<=pn)?bd:pn); i_bd = thre;\
	for (i = 0; i < bd; i++, i_bd++) {\
        w_bit((ez->Peq)[seq_nt4_table[(uint8_t)(pstr)[pidx-i]]], (i_bd));\
    }\
	bd = thre+1, i_bd = thre;\
	\
	w_##sf##_clear(ez->Peq[4]);\
	err = thre;\
	w_##sf##_set_bit_lsub(ez->VN, thre); /**VN = (((Word)1)<<(thre))-1; **/\
	w_##sf##_set_bit_lsub(ez->VP, (thre<<1)+1); /**VP = (((Word)1)<<((thre<<1)+1))-1;**/\
	w_##sf##_self_xor(ez->VP, ez->VN); /**VP ^= VN;**/\
	\
	/** print_bits(ez->VP.a, (thre<<1)+1, "-VP");**/\
	/**should make Peq[4] = 0 if N is always an error**/\
	i = 0; \
	/**for the incoming char/last char; mm = ((Word)1 << (thre<<1))**/\
    Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));\
    while (i < tn0) {\
		ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[tidx-i], c);\
		if (!(ez->D0.a[0]&(1ULL))) {\
            ++err; if (err>cut) return;\
        }\
		\
        {\
            poff = i-thre; k = i+thre-pe;/**poff:[i-thre, i+thre]**/\
            if(k >= 0) {\
                if(tmp_e == INT32_MAX) {\
                    tmp_e = err;\
                    for ((k) = 0; (poff) < (pe); (poff)++, (k)++) {\
                        bd = (k>>bitw); k_bd = (k&bitz);\
                        (tmp_e) += (((*ez).VP.a[bd]>>k_bd)&((w_sig)1));\
                        (tmp_e) -= (((*ez).VN.a[bd]>>k_bd)&((w_sig)1));\
                    }\
                } else {\
                    k = (thre<<1) - k;\
                    if(k >= 0) {\
                        bd = (k>>bitw); k_bd = (k&bitz);\
                        (tmp_e) += (((*ez).HP.a[bd]>>k_bd)&((w_sig)1));\
                        (tmp_e) -= (((*ez).HN.a[bd]>>k_bd)&((w_sig)1));\
                    }\
                }\
                if((tmp_e) <= (*ez).thre && (tmp_e) < (*ez).err) {\
                    (*ez).err = tmp_e; (*ez).ps = pidx - pe; (*ez).ts = tidx-i;\
                }\
            }\
            /**if(dbg_ext_err(i, thre, err, pe, &((*ez).VP), &((*ez).VN))!=tmp_e)**/\
        }\
		\
		/** Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;**/\
		w_##sf##_self_rsft_1(ez->Peq[0]); w_##sf##_self_rsft_1(ez->Peq[1]);\
		w_##sf##_self_rsft_1(ez->Peq[2]); w_##sf##_self_rsft_1(ez->Peq[3]);\
        ++i; ++i_bd; c = 4;\
        if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[pidx-i_bd]];\
        if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;\
    }\
	ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[tidx-i], c);\
	if (!(ez->D0.a[0]&(1ULL))) {\
		++err; if (err>cut) return;\
	}\
	\
    int32_t site = tn - 1 - thre;/**up bound; site:[tn - 1 - thre, tn - 1 + thre]**/\
	for (cut = pn - 1, i = 0; site < cut; i++) {\
		bd = (i>>bitw); i_bd = (i&bitz);\
        err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));\
        err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));\
		site++;\
		if(err <= thre && err < ez->err) {\
			ez->err = err; ez->ps = pidx-site; ez->ts = tidx+1-tn;\
		}\
	}\
	if(err <= thre && err < ez->err) {\
		ez->err = err; ez->ps = pidx-site; ez->ts = tidx+1-tn;\
	}\
    return;\
}\
inline void ed_band_cal_global_##sf##_w_trace(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)\
{\
	ez->cigar.n = 0; ez->nword = w_##sf##_word;\
	if(ez->err > thre) {\
		init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = 0;\
	} else if(ez->err == 0) {\
        push_trace(&(ez->cigar), 0, ez->te+1-ez->ts); return;\
    }\
	if((pn > tn + thre) || (tn > pn + thre)) return;\
	int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd, i_bd, ws, Peq_i; w_sig c, Peq_m;\
	w_##sf##_clear(ez->Peq[0]); w_##sf##_clear(ez->Peq[1]); w_##sf##_clear(ez->Peq[2]); w_##sf##_clear(ez->Peq[3]); w_##sf##_clear(ez->Peq[4]);\
	\
	bd = thre+1; bd = ((bd<=pn)?bd:pn); i_bd = thre;\
	ed_init_core(i, bd, i_bd, pstr, ez->Peq);\
	bd = thre+1, i_bd = thre;\
	\
	w_##sf##_clear(ez->Peq[4]);\
	err = thre;\
	w_##sf##_set_bit_lsub(ez->VN, thre); /**VN = (((Word)1)<<(thre))-1;**/\ 
	w_##sf##_set_bit_lsub(ez->VP, (thre<<1)+1); /**VP = (((Word)1)<<((thre<<1)+1))-1;**/\
	w_##sf##_self_xor(ez->VP, ez->VN); /**VP ^= VN;**/\
	\
	/** print_bits(VP.a, (thre<<1)+1, "-VP");**/\
	/**should make Peq[4] = 0 if N is always an error**/\
	i = 0; ws = sizeof(*(ez->a))*(ez->nword);\
	ez->path.n=(ez->nword*tn*5);\
    kv_resize(w_sig, ez->path, ez->path.n); ez->path.n=0;\
	\
	/**for the incoming char/last char**/\
	Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));\
    while (i < tn0) {\
		ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c);\
        if (!(ez->D0.a[0]&(1ULL))) {\
            ++err; if (err>cut) return;\
        }\
		/** Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;**/\
		w_##sf##_self_rsft_1(ez->Peq[0]); w_##sf##_self_rsft_1(ez->Peq[1]);\
		w_##sf##_self_rsft_1(ez->Peq[2]); w_##sf##_self_rsft_1(ez->Peq[3]);\
		\
		++i; ++i_bd; c = 4;\
        if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];\
        if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;\
		\
		memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;\
        memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;\
        memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;\
        memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;\
        memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;\
    }\
	ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c);\
	if (!(ez->D0.a[0]&(1ULL))) {\
		++err; if (err>cut) return;\
	}\
	memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;\
	\
    int32_t site = tn - 1 - thre;/**up bound**/\
	if(ez->err > thre) {\
		for (cut = pn - 1, i = 0; site < cut; site++, i++) {\
			bd = (i>>bitw); i_bd = (i&bitz);\
			err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));\
			err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));\
		}\
		\
		if (site == cut && err <= thre) {\
			ez->err = err; \
			ez->pe = pn-1; ez->te = tn-1;\
		}\
	}\
    gen_trace(ez, thre, 1);\
    return;\
}\
inline void ed_band_cal_semi_##sf##_w_trace(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)\
{\
	ez->cigar.n = 0; ez->nword = w_##sf##_word;\
	if(ez->err > thre) {\
		init_base_ed(*ez, thre, pn, tn); ez->ps = ez->pe = -1; ez->ts = 0; ez->te = tn-1;\
	}	else if(ez->err == 0) {\
        push_trace(&(ez->cigar), 0, ez->te+1-ez->ts);\
        ez->ps = ez->pe - (ez->te-ez->ts);\
        return;\
    }\
	int32_t bd, i, err = 0, i_bd, cut = thre+(thre<<1), tn0 = tn - 1, Peq_i, ws; w_sig c, Peq_m;\
	if((pn > tn + cut) || (tn > pn + cut)) return;\
	\
	w_##sf##_clear(ez->Peq[0]);\
	w_##sf##_clear(ez->Peq[1]);\
	w_##sf##_clear(ez->Peq[2]);\
	w_##sf##_clear(ez->Peq[3]);\
	w_##sf##_clear(ez->Peq[4]);\
	w_##sf##_clear(ez->VP);\
	w_##sf##_clear(ez->VN);\
	\
	bd = (thre<<1)+1; bd = ((bd<=pn)?bd:pn); i_bd = 0;\
	ed_init_core(i, bd, i_bd, pstr, ez->Peq);\
	bd = (thre<<1)+1, i_bd = (thre<<1);\
	\
	w_##sf##_clear(ez->Peq[4]);\
	i = 0; ws = sizeof(*(ez->a))*(ez->nword);\
	ez->path.n=(ez->nword*tn*5);\
    kv_resize(w_sig, ez->path, ez->path.n); ez->path.n=0;\
	/**for the incoming char/last char; mm = ((Word)1 << (thre<<1))**/\
    Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));\
    while (i < tn0) {\
		ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c);\
		if (!(ez->D0.a[0]&(1ULL))) {\
            ++err; if (err>cut) return;\
        }\
		/** Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;**/\
		w_##sf##_self_rsft_1(ez->Peq[0]); w_##sf##_self_rsft_1(ez->Peq[1]);\
		w_##sf##_self_rsft_1(ez->Peq[2]); w_##sf##_self_rsft_1(ez->Peq[3]);\
        ++i; ++i_bd; c = 4;\
        if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];\
        if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;\
		\
		memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;\
		memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;\
		memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;\
		memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;\
		memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;\
    }\
	ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c);\
	if (!(ez->D0.a[0]&(1ULL))) {\
		++err; if (err>cut) return;\
	}\
	memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;\
	\
    int32_t site = tn - 1;/**up bound**/\
    /**in most cases, ai = (thre<<1)**/\
    int32_t ai = pn - tn, uge = INT32_MAX;\
	if(ez->err > thre) {\
		if ((err <= thre) && (err <= ez->err)) {\
			ez->err = err; ez->pe = site;\
		}\
		\
		i = 0;\
		while (i < ai) {\
			bd = (i>>bitw); i_bd = (i&bitz);\
			err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));\
			err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));\
			++i;\
			if ((err <= thre) && (err <= ez->err)) {\
				ez->err = err; ez->pe = site + i;\
			}\
			if(i == thre) uge = err;\
		}\
		\
		if((uge <= thre) && (uge == ez->err)) ez->pe = site + thre;\
	}\
    gen_trace(ez, 0, 1);\
}\
inline void ed_band_cal_extension_##sf##_0_w_trace(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)\
{\
	int32_t done = 1;\
	ez->cigar.n = 0; ez->nword = w_##sf##_word;\
	if(ez->err > thre) {\
		init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = 0; ez->pe = ez->te = -1;\
		if(pn > tn + thre) pn = tn + thre;\
		else if(tn > pn + thre) tn = pn + thre;\
		done = 0;\
	} else if(ez->err == 0) {\
        push_trace(&(ez->cigar), 0, ez->te+1-ez->ts); return;\
    } else {\
        pn = ez->pe + 1 - ez->ps; tn = ez->te + 1 - ez->ts;\
    }\
	int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd, i_bd, k_bd, Peq_i, ws; w_sig c, Peq_m;\
	int32_t poff, pe = pn-1, tmp_e = INT32_MAX, k;\
	\
	w_##sf##_clear(ez->Peq[0]);\
	w_##sf##_clear(ez->Peq[1]);\
	w_##sf##_clear(ez->Peq[2]);\
	w_##sf##_clear(ez->Peq[3]);\
	w_##sf##_clear(ez->Peq[4]);\
	\
	bd = thre+1; bd = ((bd<=pn)?bd:pn); i_bd = thre;\
	ed_init_core(i, bd, i_bd, pstr, ez->Peq);\
	bd = thre+1, i_bd = thre;\
	\
	w_##sf##_clear(ez->Peq[4]);\
	err = thre;\
	w_##sf##_set_bit_lsub(ez->VN, thre); /**VN = (((Word)1)<<(thre))-1; **/\
	w_##sf##_set_bit_lsub(ez->VP, (thre<<1)+1); /**VP = (((Word)1)<<((thre<<1)+1))-1;**/\
	w_##sf##_self_xor(ez->VP, ez->VN); /**VP ^= VN;**/\
	\
	/** print_bits(ez->VP.a, (thre<<1)+1, "-VP");**/\
	/**should make Peq[4] = 0 if N is always an error**/\
	i = 0; ws = sizeof(*(ez->a))*(ez->nword);\
	ez->path.n=(ez->nword*tn*5);\
    kv_resize(w_sig, ez->path, ez->path.n); ez->path.n=0;\
	/**for the incoming char/last char; mm = ((Word)1 << (thre<<1))**/\
    Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));\
    while (i < tn0) {\
		ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c);\
		if (!(ez->D0.a[0]&(1ULL))) {\
            ++err; if (err>cut) return;\
        }\
		\
        if(!done) {\
            poff = i-thre; k = i+thre-pe;/**poff:[i-thre, i+thre]**/\
            if(k >= 0) {\
                if(tmp_e == INT32_MAX) {\
                    tmp_e = err;\
                    for ((k) = 0; (poff) < (pe); (poff)++, (k)++) {\
                        bd = (k>>bitw); k_bd = (k&bitz);\
                        (tmp_e) += (((*ez).VP.a[bd]>>k_bd)&((w_sig)1));\
                        (tmp_e) -= (((*ez).VN.a[bd]>>k_bd)&((w_sig)1));\
                    }\
                } else {\
                    k = (thre<<1) - k;\
                    if(k >= 0) {\
                        bd = (k>>bitw); k_bd = (k&bitz);\
                        (tmp_e) += (((*ez).HP.a[bd]>>k_bd)&((w_sig)1));\
                        (tmp_e) -= (((*ez).HN.a[bd]>>k_bd)&((w_sig)1));\
                    }\
                }\
                if((tmp_e) <= (*ez).thre && (tmp_e) < (*ez).err) {\
                    (*ez).err = tmp_e; (*ez).pe = pe; (*ez).te = i;\
                }\
            }\
            /**if(dbg_ext_err(i, thre, err, pe, &((*ez).VP), &((*ez).VN))!=tmp_e)**/\
        }\
		\
		/** Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;**/\
		w_##sf##_self_rsft_1(ez->Peq[0]); w_##sf##_self_rsft_1(ez->Peq[1]);\
		w_##sf##_self_rsft_1(ez->Peq[2]); w_##sf##_self_rsft_1(ez->Peq[3]);\
        ++i; ++i_bd; c = 4;\
        if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];\
        if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;\
		\
		memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;\
		memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;\
		memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;\
		memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;\
		memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;\
    }\
	ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c);\
	if (!(ez->D0.a[0]&(1ULL))) {\
		++err; if (err>cut) return;\
	}\
	memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;\
	\
    int32_t site = tn - 1 - thre;/**up bound; site:[tn - 1 - thre, tn - 1 + thre]**/\
	if(!done) {\
		for (cut = pn - 1, i = 0; site < cut; i++) {\
			bd = (i>>bitw); i_bd = (i&bitz);\
			err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));\
			err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));\
			site++;\
			if(err <= thre && err < ez->err) {\
				ez->err = err; ez->pe = site; ez->te = tn-1;\
			}\
		}\
		if(err <= thre && err < ez->err) {\
			ez->err = err; ez->pe = site; ez->te = tn-1;\
		}\
	}\
	\
	if((ez->te-ez->ts+1) != tn) {\
        ez->path.n /= tn; ez->path.n *= (ez->te+1-ez->ts);\
    }\
    gen_trace(ez, thre, 1);\
    return;\
}\
inline void ed_band_cal_extension_##sf##_1_w_trace(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)\
{\
	int32_t done = 1;\
	ez->cigar.n = 0; ez->nword = w_##sf##_word;\
	if(ez->err > thre) {\
		init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = INT32_MAX; ez->pe = pn-1; ez->te = tn-1;\
		if(pn > tn + thre) pn = tn + thre;\
		else if(tn > pn + thre) tn = pn + thre;\
		done = 0;\
    } else if(ez->err == 0) {\
        push_trace(&(ez->cigar), 0, ez->te+1-ez->ts); return;\
    } else {\
        pn = ez->pe + 1 - ez->ps; tn = ez->te + 1 - ez->ts;\
    }\
	int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd, i_bd, k_bd, Peq_i, ws; w_sig c, Peq_m;\
	int32_t poff, pe = pn-1, tmp_e = INT32_MAX, k, pidx = ez->pe, tidx = ez->te;\
	\
	w_##sf##_clear(ez->Peq[0]);\
	w_##sf##_clear(ez->Peq[1]);\
	w_##sf##_clear(ez->Peq[2]);\
	w_##sf##_clear(ez->Peq[3]);\
	w_##sf##_clear(ez->Peq[4]);\
	\
	bd = thre+1; bd = ((bd<=pn)?bd:pn); i_bd = thre;\
	for (i = 0; i < bd; i++, i_bd++) {\
        w_bit((ez->Peq)[seq_nt4_table[(uint8_t)(pstr)[pidx-i]]], (i_bd));\
    }\
	bd = thre+1, i_bd = thre;\
	\
	w_##sf##_clear(ez->Peq[4]);\
	err = thre;\
	w_##sf##_set_bit_lsub(ez->VN, thre); /**VN = (((Word)1)<<(thre))-1; **/\
	w_##sf##_set_bit_lsub(ez->VP, (thre<<1)+1); /**VP = (((Word)1)<<((thre<<1)+1))-1;**/\
	w_##sf##_self_xor(ez->VP, ez->VN); /**VP ^= VN;**/\
	\
	/** print_bits(ez->VP.a, (thre<<1)+1, "-VP");**/\
	/**should make Peq[4] = 0 if N is always an error**/\
	i = 0; ws = sizeof(*(ez->a))*(ez->nword);\
	ez->path.n=(ez->nword*tn*5);\
    kv_resize(w_sig, ez->path, ez->path.n); ez->path.n=0;\
	/**for the incoming char/last char; mm = ((Word)1 << (thre<<1))**/\
    Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));\
    while (i < tn0) {\
		ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[tidx-i], c);\
		if (!(ez->D0.a[0]&(1ULL))) {\
            ++err; if (err>cut) return;\
        }\
		\
        if(!done) {\
            poff = i-thre; k = i+thre-pe;/**poff:[i-thre, i+thre]**/\
            if(k >= 0) {\
                if(tmp_e == INT32_MAX) {\
                    tmp_e = err;\
                    for ((k) = 0; (poff) < (pe); (poff)++, (k)++) {\
                        bd = (k>>bitw); k_bd = (k&bitz);\
                        (tmp_e) += (((*ez).VP.a[bd]>>k_bd)&((w_sig)1));\
                        (tmp_e) -= (((*ez).VN.a[bd]>>k_bd)&((w_sig)1));\
                    }\
                } else {\
                    k = (thre<<1) - k;\
                    if(k >= 0) {\
                        bd = (k>>bitw); k_bd = (k&bitz);\
                        (tmp_e) += (((*ez).HP.a[bd]>>k_bd)&((w_sig)1));\
                        (tmp_e) -= (((*ez).HN.a[bd]>>k_bd)&((w_sig)1));\
                    }\
                }\
                if((tmp_e) <= (*ez).thre && (tmp_e) < (*ez).err) {\
                    (*ez).err = tmp_e; (*ez).ps = pidx - pe; (*ez).ts = tidx-i;\
                }\
            }\
            /**if(dbg_ext_err(i, thre, err, pe, &((*ez).VP), &((*ez).VN))!=tmp_e)**/\
        }\
		\
		/** Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;**/\
		w_##sf##_self_rsft_1(ez->Peq[0]); w_##sf##_self_rsft_1(ez->Peq[1]);\
		w_##sf##_self_rsft_1(ez->Peq[2]); w_##sf##_self_rsft_1(ez->Peq[3]);\
        ++i; ++i_bd; c = 4;\
        if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[pidx-i_bd]];\
        if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;\
		\
		memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;\
		memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;\
		memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;\
		memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;\
		memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;\
    }\
	ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[tidx-i], c);\
	if (!(ez->D0.a[0]&(1ULL))) {\
		++err; if (err>cut) return;\
	}\
	memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;\
	\
    int32_t site = tn - 1 - thre;/**up bound; site:[tn - 1 - thre, tn - 1 + thre]**/\
	if(!done) {\
		for (cut = pn - 1, i = 0; site < cut; i++) {\
			bd = (i>>bitw); i_bd = (i&bitz);\
			err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));\
			err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));\
			site++;\
			if(err <= thre && err < ez->err) {\
				ez->err = err; ez->ps = pidx-site; ez->ts = tidx+1-tn;\
			}\
		}\
		if(err <= thre && err < ez->err) {\
			ez->err = err; ez->ps = pidx-site; ez->ts = tidx+1-tn;\
		}\
	}\
	\
	if((ez->te-ez->ts+1) != tn) {\
        ez->path.n /= tn; ez->path.n *= (ez->te+1-ez->ts);\
    }\
	\
    poff = ez->ps; ez->ps = pidx - ez->pe; ez->pe = pidx - poff;\
    gen_trace(ez, thre, 0);\
    poff = ez->ps; ez->ps = pidx - ez->pe; ez->pe = pidx - poff;\
    return;\
}\
inline void ed_band_cal_semi_##sf##_w_absent_diag(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, int32_t abs_diag, bit_extz_t *ez)\
{\
	init_base_ed(*ez, thre, pn, tn); ez->ps = ez->pe = -1; ez->ts = 0; ez->te = tn-1;\
	int32_t bd, i, err = abs_diag, i_bd, cut = thre+(thre<<1), tn0 = tn - 1, Peq_i; w_sig c, Peq_m;\
	if((pn > tn + cut) || (tn > pn + cut)) return;\
	\
	w_##sf##_clear(ez->Peq[0]);\
	w_##sf##_clear(ez->Peq[1]);\
	w_##sf##_clear(ez->Peq[2]);\
	w_##sf##_clear(ez->Peq[3]);\
	w_##sf##_clear(ez->Peq[4]);\
	w_##sf##_clear(ez->VP);\
	w_##sf##_set_bit_lsub(ez->VN, abs_diag);\
	\
	bd = ((thre<<1)+1)-abs_diag; bd = ((bd<=pn)?bd:pn); i_bd = abs_diag; \
	ed_init_core(i, bd, i_bd, pstr, ez->Peq);\
	i_bd = (thre<<1)-abs_diag;\
	\
	w_##sf##_clear(ez->Peq[4]);\
	i = 0;\
	/**for the incoming char/last char; mm = ((Word)1 << (thre<<1))**/\
    Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));\
    while (i < tn0) {\
		ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c);\
		if (!(ez->D0.a[0]&(1ULL))) {\
            ++err; if (err>cut) return;\
        }\
		/** Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;**/\
		w_##sf##_self_rsft_1(ez->Peq[0]); w_##sf##_self_rsft_1(ez->Peq[1]);\
		w_##sf##_self_rsft_1(ez->Peq[2]); w_##sf##_self_rsft_1(ez->Peq[3]);\
        ++i; ++i_bd; c = 4;\
        if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];\
        if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;\
    }\
	ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c);\
	if (!(ez->D0.a[0]&(1ULL))) {\
		++err; if (err>cut) return;\
	}\
	\
    int32_t site = tn - 1 - abs_diag;/**up bound**/\
    /**in most cases, ai = (thre<<1)**/\
    int32_t ai = pn - tn + abs_diag, uge = INT32_MAX; i = 0;\
	for (i = 0; site < 0 && i < ai; i++, site++) {\
        bd = (i>>bitw); i_bd = (i&bitz);\
        err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));\
        err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));\
    }\
    if ((err <= thre) && (err <= ez->err)) {\
        ez->err = err; ez->pe = site;\
    }\
	\
    site -= i;\
    while (i < ai) {\
        bd = (i>>bitw); i_bd = (i&bitz);\
        err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));\
        err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));\
		++i;\
        if ((err <= thre) && (err <= ez->err)) {\
            ez->err = err; ez->pe = site + i;\
        }\
        if(i == thre) uge = err;\
    }\
	\
    if((uge <= thre) && (uge == ez->err)) ez->pe = site + thre;\
}\
inline void ed_band_cal_semi_##sf##_w_absent_diag_trace(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, int32_t abs_diag, bit_extz_t *ez)\
{\
	ez->cigar.n = 0; ez->nword = w_##sf##_word;\
	if(ez->err > thre) {\
		init_base_ed(*ez, thre, pn, tn); ez->ps = ez->pe = -1; ez->ts = 0; ez->te = tn-1;\
	}	else if(ez->err == 0) {\
        push_trace(&(ez->cigar), 0, ez->te+1-ez->ts);\
        ez->ps = ez->pe - (ez->te-ez->ts);\
        return;\
    }\
	int32_t bd, i, err = abs_diag, i_bd, cut = thre+(thre<<1), tn0 = tn - 1, Peq_i, ws; w_sig c, Peq_m;\
	if((pn > tn + cut) || (tn > pn + cut)) return;\
	\
	w_##sf##_clear(ez->Peq[0]);\
	w_##sf##_clear(ez->Peq[1]);\
	w_##sf##_clear(ez->Peq[2]);\
	w_##sf##_clear(ez->Peq[3]);\
	w_##sf##_clear(ez->Peq[4]);\
	w_##sf##_clear(ez->VP);\
	w_##sf##_set_bit_lsub(ez->VN, abs_diag);\
	\
	bd = ((thre<<1)+1)-abs_diag; bd = ((bd<=pn)?bd:pn); i_bd = abs_diag;\
    ed_init_core(i, bd, i_bd, pstr, ez->Peq);\
    i_bd = (thre<<1)-abs_diag;\
	\
	w_##sf##_clear(ez->Peq[4]);\
	i = 0; ws = sizeof(*(ez->a))*(ez->nword);\
	ez->path.n=(ez->nword*tn*5);\
    kv_resize(w_sig, ez->path, ez->path.n); ez->path.n=0;\
	/**for the incoming char/last char; mm = ((Word)1 << (thre<<1))**/\
    Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));\
    while (i < tn0) {\
		ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c);\
		if (!(ez->D0.a[0]&(1ULL))) {\
            ++err; if (err>cut) return;\
        }\
		/** Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;**/\
		w_##sf##_self_rsft_1(ez->Peq[0]); w_##sf##_self_rsft_1(ez->Peq[1]);\
		w_##sf##_self_rsft_1(ez->Peq[2]); w_##sf##_self_rsft_1(ez->Peq[3]);\
        ++i; ++i_bd; c = 4;\
        if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];\
        if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;\
		\
		memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;\
		memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;\
		memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;\
		memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;\
		memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;\
    }\
	ed_core(_##sf##_, ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c);\
	if (!(ez->D0.a[0]&(1ULL))) {\
		++err; if (err>cut) return;\
	}\
	memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;\
	memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;\
	\
    int32_t site = tn - 1 - abs_diag;/**up bound**/\
    /**in most cases, ai = (thre<<1)**/\
    int32_t ai = pn - tn + abs_diag, uge = INT32_MAX; i = 0;\
	if(ez->err > thre) {\
		for (i = 0; site < 0 && i < ai; i++, site++) {\
            bd = (i>>bitw); i_bd = (i&bitz);\
            err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));\
            err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));\
        }\
		if ((err <= thre) && (err <= ez->err)) {\
			ez->err = err; ez->pe = site;\
		}\
		\
		site -= i;\
		while (i < ai) {\
			bd = (i>>bitw); i_bd = (i&bitz);\
			err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));\
			err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));\
			++i;\
			if ((err <= thre) && (err <= ez->err)) {\
				ez->err = err; ez->pe = site + i;\
			}\
			if(i == thre) uge = err;\
		}\
		\
		if((uge <= thre) && (uge == ez->err)) ez->pe = site + thre;\
	}\
    gen_trace(ez, abs_diag, 1);\
}\

HA_ED_INIT(128)
HA_ED_INIT(192)
HA_ED_INIT(256)


#define w_infi_clear(x, nw) (memset((x).a, 0, sizeof(*((x).a))*(nw)))

#define w_infi_set_bit_lsub(x, l, nw) do {	\
		w_infi_clear(x, nw);\
		if((l)>>bitw) memset((x).a, -1, sizeof(*((x).a))*((l)>>bitw));\
		if((l)&bitz) (x).a[(l)>>bitw] = (((w_sig)1)<<((l)&bitz))-1; \
	} while (0)	

#define w_infi_self_xor(x, y, nw, w_z) {\
		for((w_z)=0; (w_z)<(nw); (w_z)++) (x).a[w_z]^=(y).a[w_z];}

#define w_self_add_sin(x, y, w_z, ad) {\
	(x).a[(w_z)]+=(ad), (ad)=((x).a[(w_z)]<(ad)), (x).a[(w_z)]+=(y).a[(w_z)], (ad)|=((x).a[(w_z)]<(y).a[(w_z)]);}

#define ed_infi_core(Peq, VP, VN, X, D0, HN, HP, z, c, ad, w_z, nw) {\
	(c) = seq_nt4_table[(z)];\
	for((w_z)=(ad)=0; (w_z)<(nw); (w_z)++){\
		/**X = Peq[seq_nt4_table[(uint8_t)tstr[i]]] | VN;**/\
		(X).a[(w_z)]=(Peq)[(c)].a[(w_z)]|(VN).a[(w_z)];\
		/**D0 = ((VP + (X&VP)) ^ VP) | X;**/\
		(D0).a[(w_z)]=(X).a[(w_z)]&(VP).a[(w_z)];\
		w_self_add_sin((D0), (VP), (w_z), (ad));\
		(D0).a[(w_z)]^=(VP).a[(w_z)];\
		(D0).a[(w_z)]|=(X).a[(w_z)];\
		/**HN = VP&D0;**/\
		(HN).a[(w_z)]=(VP).a[(w_z)]&(D0).a[(w_z)];\
		/**HP = VN | ~(VP | D0);**/\
		(HP).a[(w_z)]=~((VP).a[(w_z)]|(D0).a[(w_z)]);\
		(HP).a[(w_z)]|=(VN).a[(w_z)];\
	}\
	for((w_z)=(nw)-1,(ad)=0; (w_z)>=0; (w_z)--){\
		/**X = D0 >> 1;**/\
		(X).a[(w_z)]=((D0).a[(w_z)]>>1)|(ad),(ad)=(D0).a[(w_z)]<<bitz;\
		/**VN = X&HP;**/\
		(VN).a[(w_z)]=(X).a[(w_z)]&(HP).a[(w_z)];\
		/**VP = HN | ~(X | HP);**/\
		(VP).a[(w_z)]=~((X).a[(w_z)]|(HP).a[(w_z)]);\
		(VP).a[(w_z)]|=(HN).a[(w_z)];\
	}\
}


#define ed_infi_post_Peq(Peq, w_z, w_z1, nw) {\
	/** Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;**/\
	for((w_z)=0, (w_z1)=1; (w_z1)<(nw); (w_z)++,(w_z1)++){\
		Peq[0].a[(w_z)]=(Peq[0].a[(w_z)]>>1)|(Peq[0].a[(w_z1)]<<bitz);\
		Peq[1].a[(w_z)]=(Peq[1].a[(w_z)]>>1)|(Peq[1].a[(w_z1)]<<bitz);\
		Peq[2].a[(w_z)]=(Peq[2].a[(w_z)]>>1)|(Peq[2].a[(w_z1)]<<bitz);\
		Peq[3].a[(w_z)]=(Peq[3].a[(w_z)]>>1)|(Peq[3].a[(w_z1)]<<bitz);\
	}\
	(w_z)=(nw)-1;\
	Peq[0].a[(w_z)]>>=1;\
	Peq[1].a[(w_z)]>>=1;\
	Peq[2].a[(w_z)]>>=1;\
	Peq[3].a[(w_z)]>>=1;\
}

inline void ed_band_cal_global_infi_w(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, int32_t *nword, bit_extz_t *ez)
{
	init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = 0;
	if((pn > tn + thre) || (tn > pn + thre)) return;
	int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd, i_bd, wz, wz1, Peq_i; w_sig c, ad, Peq_m;
	if(nword) {
		ez->nword = (*nword);
	} else {
		bd = (((thre)<<1)+1); ez->nword = ((bd>>bitw)+(!!(bd&bitz)));
	}
	resize_bit_extz_t((*ez), (thre)); 
	wz = sizeof(*(ez->a))*(ez->nword);
	memset((ez->Peq[0]).a, 0, wz);
	memset((ez->Peq[1]).a, 0, wz);
	memset((ez->Peq[2]).a, 0, wz);
	memset((ez->Peq[3]).a, 0, wz);
	memset((ez->Peq[4]).a, 0, wz);

	bd = thre+1; bd = ((bd<=pn)?bd:pn); i_bd = thre;
	ed_init_core(i, bd, i_bd, pstr, ez->Peq);
	bd = thre+1, i_bd = thre;
	
	memset((ez->Peq[4]).a, 0, wz);
	err = thre;
	w_infi_set_bit_lsub(ez->VN, thre, ez->nword); /**VN = (((Word)1)<<(thre))-1;**/ 
	w_infi_set_bit_lsub(ez->VP, (thre<<1)+1, ez->nword); /**VP = (((Word)1)<<((thre<<1)+1))-1;**/
	w_infi_self_xor(ez->VP, ez->VN, ez->nword, wz) /**VP ^= VN;**/
	
	i = 0; 
	/**for the incoming char/last char; mm = ((Word)1 << (thre<<1))**/
	Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));
    while (i < tn0) {
		// fprintf(stderr, "\ni::%d\n", i);
		// prt_bit_extz_t((*ez), (((thre<<1))+1));
		ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c, ad, wz, ez->nword);
        if (!(ez->D0.a[0]&(1ULL))) {
            ++err; if (err>cut) return;
        }

		ed_infi_post_Peq(ez->Peq, wz, wz1, ez->nword);
		++i; ++i_bd; c = 4;
		if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];
		if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;
		// fprintf(stderr, "c::%u\n", c);
		// print_bits(ez->Peq[c].a, (((thre<<1))+1), "Peq[c]");
    }
	ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c, ad, wz, ez->nword);
	if (!(ez->D0.a[0]&(1ULL))) {
		++err; if (err>cut) return;
	}
	
    int32_t site = tn - 1 - thre;/**up bound**/
	for (cut = pn - 1, i = 0; site < cut; site++, i++) {
		bd = (i>>bitw); i_bd = (i&bitz);
		err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));
		err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));
	}
	
	if (site == cut && err <= thre) {
		ez->err = err; 
		ez->pe = pn-1; ez->te = tn-1;
	}
    return;
}

inline void ed_band_cal_semi_infi_w(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, int32_t *nword, bit_extz_t *ez)
{
	init_base_ed(*ez, thre, pn, tn); ez->ps = ez->pe = -1; ez->ts = 0; ez->te = tn-1;
	int32_t bd, i, err = 0, i_bd, cut = thre+(thre<<1), tn0 = tn - 1, wz, wz1, Peq_i; w_sig c, ad, Peq_m;
	if((pn > tn + cut) || (tn > pn + cut)) return;

	if(nword) {
        ez->nword = (*nword);
    } else {
        bd = (((thre)<<1)+1); ez->nword = ((bd>>bitw)+(!!(bd&bitz)));
    }
    resize_bit_extz_t((*ez), (thre)); 
	wz = sizeof(*(ez->a))*(ez->nword);
	memset((ez->Peq[0]).a, 0, wz);
	memset((ez->Peq[1]).a, 0, wz);
	memset((ez->Peq[2]).a, 0, wz);
	memset((ez->Peq[3]).a, 0, wz);
	memset((ez->Peq[4]).a, 0, wz);
	memset((ez->VP).a, 0, wz);
	memset((ez->VN).a, 0, wz);

	bd = (thre<<1)+1; bd = ((bd<=pn)?bd:pn); i_bd = 0;
    ed_init_core(i, bd, i_bd, pstr, ez->Peq);
    bd = (thre<<1)+1, i_bd = (thre<<1);

	memset((ez->Peq[4]).a, 0, wz);
	i = 0; 
    /**for the incoming char/last char; mm = ((Word)1 << (thre<<1))**/
    Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));
	while (i < tn0) {
        ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c, ad, wz, ez->nword);
        if (!(ez->D0.a[0]&(1ULL))) {
            ++err; if (err>cut) return;
        }
        
        ed_infi_post_Peq(ez->Peq, wz, wz1, ez->nword);
        ++i; ++i_bd; c = 4;
        if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];
        if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;
    }
	ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c, ad, wz, ez->nword);
    if (!(ez->D0.a[0]&(1ULL))) {
        ++err; if (err>cut) return;
    }

	int32_t site = tn - 1;/**up bound**/
	/**in most cases, ai = (thre<<1)**/
    int32_t ai = pn - tn, uge = INT32_MAX;
    if ((err <= thre) && (err <= ez->err)) {
        ez->err = err; ez->pe = site;
    }

	i = 0;
	while (i < ai) {
		bd = (i>>bitw); i_bd = (i&bitz);
		err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));
		err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));
        ++i;
        if ((err <= thre) && (err <= ez->err)) {
            ez->err = err; ez->pe = site + i;
        }
        if(i == thre) uge = err;
    }
    
    if((uge <= thre) && (uge == ez->err)) ez->pe = site + thre;
}

inline void ed_band_cal_extension_infi_0_w(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, int32_t *nword, bit_extz_t *ez)
{
    init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = 0; ez->pe = ez->te = -1;
	if(pn > tn + thre) pn = tn + thre;
    else if(tn > pn + thre) tn = pn + thre;
    int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd, i_bd, k_bd, wz, wz1, Peq_i; w_sig c, ad, Peq_m;
    int32_t poff, pe = pn-1, tmp_e = INT32_MAX, k;
	if(nword) {
        ez->nword = (*nword);
    } else {
        bd = (((thre)<<1)+1); ez->nword = ((bd>>bitw)+(!!(bd&bitz)));
    }
    resize_bit_extz_t((*ez), (thre)); 
	wz = sizeof(*(ez->a))*(ez->nword);
    memset((ez->Peq[0]).a, 0, wz);
    memset((ez->Peq[1]).a, 0, wz);
    memset((ez->Peq[2]).a, 0, wz);
    memset((ez->Peq[3]).a, 0, wz);
    memset((ez->Peq[4]).a, 0, wz);

	bd = thre+1; bd = ((bd<=pn)?bd:pn); i_bd = thre;
    ed_init_core(i, bd, i_bd, pstr, ez->Peq);
    bd = thre+1, i_bd = thre;

	memset((ez->Peq[4]).a, 0, wz);
    err = thre;
	w_infi_set_bit_lsub(ez->VN, thre, ez->nword); /**VN = (((Word)1)<<(thre))-1;**/ 
    w_infi_set_bit_lsub(ez->VP, (thre<<1)+1, ez->nword); /**VP = (((Word)1)<<((thre<<1)+1))-1;**/
    w_infi_self_xor(ez->VP, ez->VN, ez->nword, wz) /**VP ^= VN;**/

	i = 0;
	/**for the incoming char/last char; mm = ((Word)1 << (thre<<1))**/
    Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));
	//for debug
	// memset((ez->X).a, 0, sizeof(*(ez->a))*(ez->nword));
	// memset((ez->D0).a, 0, sizeof(*(ez->a))*(ez->nword));
	// memset((ez->HN).a, 0, sizeof(*(ez->a))*(ez->nword));
	// memset((ez->HP).a, 0, sizeof(*(ez->a))*(ez->nword));
	while (i < tn0) {
		// fprintf(stderr, "i::%d\n", i);
		// prt_bit_extz_t((*ez), ((thre<<1)+1));
        ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c, ad, wz, ez->nword);
        if (!(ez->D0.a[0]&(1ULL))) {
            ++err; if (err>cut) return;
        }

		{
            poff = i-thre; k = i+thre-pe;/**poff:[i-thre, i+thre]**/
			if(k >= 0) {
				if(tmp_e == INT32_MAX) {
					tmp_e = err;
					for ((k) = 0; (poff) < (pe); (poff)++, (k)++) {
						bd = (k>>bitw); k_bd = (k&bitz);
						(tmp_e) += (((*ez).VP.a[bd]>>k_bd)&((w_sig)1)); 
						(tmp_e) -= (((*ez).VN.a[bd]>>k_bd)&((w_sig)1));
					} 
				} else {
					k = (thre<<1) - k; 
					if(k >= 0) {
						bd = (k>>bitw); k_bd = (k&bitz);
						(tmp_e) += (((*ez).HP.a[bd]>>k_bd)&((w_sig)1)); 
						(tmp_e) -= (((*ez).HN.a[bd]>>k_bd)&((w_sig)1));
					}
				}
				if((tmp_e) <= (*ez).thre && (tmp_e) < (*ez).err) {
                    (*ez).err = tmp_e; (*ez).pe = pe; (*ez).te = i;
				}
			}
			// if(dbg_ext_err(i, thre, err, pe, &((*ez).VP), &((*ez).VN))!=tmp_e) {
			// 	fprintf(stderr, "i::%d, tmp_e::%d, dbg_ext_err::%d\n", i, tmp_e, dbg_ext_err(i, thre, err, pe, &((*ez).VP), &((*ez).VN)));
			// }
        }

        ed_infi_post_Peq(ez->Peq, wz, wz1, ez->nword);
        ++i; ++i_bd; c = 4;
        if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];
        if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;
    }
    ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c, ad, wz, ez->nword);
    if (!(ez->D0.a[0]&(1ULL))) {
        ++err; if (err>cut) return;
    }

    int32_t site = tn - 1 - thre;/**up bound; site:[tn - 1 - thre, tn - 1 + thre]**/
    for (cut = pn - 1, i = 0; site < cut; i++) {
		bd = (i>>bitw); i_bd = (i&bitz);
        err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));
        err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));
        site++;
        if(err <= thre && err < ez->err) {
            ez->err = err; ez->pe = site; ez->te = tn-1;
        }
    }
    if(err <= thre && err < ez->err) {
        ez->err = err; ez->pe = site; ez->te = tn-1;
    }
    return;
}

inline void ed_band_cal_extension_infi_1_w(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, int32_t *nword, bit_extz_t *ez)
{
	init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = INT32_MAX; ez->pe = pn-1; ez->te = tn-1;
	if(pn > tn + thre) pn = tn + thre;
    else if(tn > pn + thre) tn = pn + thre;
	int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd, i_bd, k_bd, wz, wz1, Peq_i; 
    int32_t poff, pe = pn-1, tmp_e = INT32_MAX, k, pidx = ez->pe, tidx = ez->te; w_sig c, ad, Peq_m;
	if(nword) {
        ez->nword = (*nword);
    } else {
        bd = (((thre)<<1)+1); ez->nword = ((bd>>bitw)+(!!(bd&bitz)));
    }
    resize_bit_extz_t((*ez), (thre)); 
	wz = sizeof(*(ez->a))*(ez->nword);
    memset((ez->Peq[0]).a, 0, wz);
    memset((ez->Peq[1]).a, 0, wz);
    memset((ez->Peq[2]).a, 0, wz);
    memset((ez->Peq[3]).a, 0, wz);
    memset((ez->Peq[4]).a, 0, wz);

	bd = thre+1; bd = ((bd<=pn)?bd:pn); i_bd = thre;
	for (i = 0; i < bd; i++, i_bd++) { 
		w_bit((ez->Peq)[seq_nt4_table[(uint8_t)(pstr)[pidx-i]]], (i_bd));
	}
    bd = thre+1, i_bd = thre;

	memset((ez->Peq[4]).a, 0, wz);
    err = thre;
	w_infi_set_bit_lsub(ez->VN, thre, ez->nword); /**VN = (((Word)1)<<(thre))-1;**/ 
    w_infi_set_bit_lsub(ez->VP, (thre<<1)+1, ez->nword); /**VP = (((Word)1)<<((thre<<1)+1))-1;**/
    w_infi_self_xor(ez->VP, ez->VN, ez->nword, wz) /**VP ^= VN;**/

	i = 0;
	/**for the incoming char/last char; mm = ((Word)1 << (thre<<1))**/
    Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));
	while (i < tn0) {
        ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[tidx-i], c, ad, wz, ez->nword);
        if (!(ez->D0.a[0]&(1ULL))) {
            ++err; if (err>cut) return;
        }

		{
            poff = i-thre; k = i+thre-pe;/**poff:[i-thre, i+thre]**/
			if(k >= 0) {
				if(tmp_e == INT32_MAX) {
					tmp_e = err;
					for ((k) = 0; (poff) < (pe); (poff)++, (k)++) {
						bd = (k>>bitw); k_bd = (k&bitz);
						(tmp_e) += (((*ez).VP.a[bd]>>k_bd)&((w_sig)1)); 
						(tmp_e) -= (((*ez).VN.a[bd]>>k_bd)&((w_sig)1));
					} 
				} else {
					k = (thre<<1) - k; 
					if(k >= 0) {
						bd = (k>>bitw); k_bd = (k&bitz);
						(tmp_e) += (((*ez).HP.a[bd]>>k_bd)&((w_sig)1)); 
						(tmp_e) -= (((*ez).HN.a[bd]>>k_bd)&((w_sig)1));
					}
				}
				if((tmp_e) <= (*ez).thre && (tmp_e) < (*ez).err) {
                    (*ez).err = tmp_e; (*ez).ps = pidx - pe; (*ez).ts = tidx-i;
				}
			}
        }
		// if(dbg_ext_err(i, thre, err, pe, &((*ez).VP), &((*ez).VN))!=tmp_e) {
		// 	fprintf(stderr, "i::%d, tmp_e::%d, dbg_ext_err::%d\n", i, tmp_e, dbg_ext_err(i, thre, err, pe, &((*ez).VP), &((*ez).VN)));
		// }

        ed_infi_post_Peq(ez->Peq, wz, wz1, ez->nword);
        ++i; ++i_bd; c = 4;
        if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[pidx-i_bd]];
        if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;
    }
    ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[tidx-i], c, ad, wz, ez->nword);
    if (!(ez->D0.a[0]&(1ULL))) {
        ++err; if (err>cut) return;
    }

    int32_t site = tn - 1 - thre;/**up bound; site:[tn - 1 - thre, tn - 1 + thre]**/
    for (cut = pn - 1, i = 0; site < cut; i++) {
		bd = (i>>bitw); i_bd = (i&bitz);
        err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));
        err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));
        site++;
        if(err <= thre && err < ez->err) {
            ez->err = err; ez->ps = pidx-site; ez->ts = tidx+1-tn;
        }
    }
    if(err <= thre && err < ez->err) {
        ez->err = err; ez->ps = pidx-site; ez->ts = tidx+1-tn;
    }
    return;
}

inline void ed_band_cal_global_infi_w_trace(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, int32_t *nword, bit_extz_t *ez)
{
	ez->cigar.n = 0;//diff
	if(ez->err > thre) {//diff
		init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = 0;
	} else if(ez->err == 0) {//diff
        push_trace(&(ez->cigar), 0, ez->te+1-ez->ts); return;//diff
    }
	if((pn > tn + thre) || (tn > pn + thre)) return;
	int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd, i_bd, wz, wz1, Peq_i, ws; w_sig c, ad, Peq_m;
	if(nword) {
		ez->nword = (*nword);
	} else {
		bd = (((thre)<<1)+1); ez->nword = ((bd>>bitw)+(!!(bd&bitz)));
	}
	resize_bit_extz_t((*ez), (thre)); 
	wz = ws = sizeof(*(ez->a))*(ez->nword);//diff
	memset((ez->Peq[0]).a, 0, wz);
	memset((ez->Peq[1]).a, 0, wz);
	memset((ez->Peq[2]).a, 0, wz);
	memset((ez->Peq[3]).a, 0, wz);
	memset((ez->Peq[4]).a, 0, wz);

	bd = thre+1; bd = ((bd<=pn)?bd:pn); i_bd = thre;
	ed_init_core(i, bd, i_bd, pstr, ez->Peq);
	bd = thre+1, i_bd = thre;
	
	memset((ez->Peq[4]).a, 0, wz);
	err = thre;
	w_infi_set_bit_lsub(ez->VN, thre, ez->nword); /**VN = (((Word)1)<<(thre))-1;**/ 
	w_infi_set_bit_lsub(ez->VP, (thre<<1)+1, ez->nword); /**VP = (((Word)1)<<((thre<<1)+1))-1;**/
	w_infi_self_xor(ez->VP, ez->VN, ez->nword, wz) /**VP ^= VN;**/

	ez->path.n=(ez->nword*tn*5);//diff
    kv_resize(w_sig, ez->path, ez->path.n); ez->path.n=0;//diff
	
	i = 0; 
	/**for the incoming char/last char; mm = ((Word)1 << (thre<<1))**/
	Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));
    while (i < tn0) {
		// fprintf(stderr, "\ni::%d\n", i);
		// prt_bit_extz_t((*ez), (((thre<<1))+1));
		ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c, ad, wz, ez->nword);
        if (!(ez->D0.a[0]&(1ULL))) {
            ++err; if (err>cut) return;
        }

		ed_infi_post_Peq(ez->Peq, wz, wz1, ez->nword);
		++i; ++i_bd; c = 4;
		if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];
		if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;
		// fprintf(stderr, "c::%u\n", c);
		// print_bits(ez->Peq[c].a, (((thre<<1))+1), "Peq[c]");

		memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;//diff
		memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;//diff
		memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;//diff
		memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;//diff
		memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;//diff
    }
	ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c, ad, wz, ez->nword);
	if (!(ez->D0.a[0]&(1ULL))) {
		++err; if (err>cut) return;
	}
	memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;//diff
	
    int32_t site = tn - 1 - thre;/**up bound**/
	if(ez->err > thre) {//diff
		for (cut = pn - 1, i = 0; site < cut; site++, i++) {
			bd = (i>>bitw); i_bd = (i&bitz);
			err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));
			err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));
		}
		
		if (site == cut && err <= thre) {
			ez->err = err; 
			ez->pe = pn-1; ez->te = tn-1;
		}
	}
	gen_trace(ez, thre, 1);//diff
    return;
}

inline void ed_band_cal_semi_infi_w_trace(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, int32_t *nword, bit_extz_t *ez)
{
	ez->cigar.n = 0;//diff
	if(ez->err > thre) {//diff
		init_base_ed(*ez, thre, pn, tn); ez->ps = ez->pe = -1; ez->ts = 0; ez->te = tn-1;
	} else if(ez->err == 0) {//diff
        push_trace(&(ez->cigar), 0, ez->te+1-ez->ts); //diff
        ez->ps = ez->pe - (ez->te-ez->ts);//diff
        return;//diff
    }
	int32_t bd, i, err = 0, i_bd, cut = thre+(thre<<1), tn0 = tn - 1, wz, wz1, ws, Peq_i; w_sig c, ad, Peq_m;
	if((pn > tn + cut) || (tn > pn + cut)) return;

	if(nword) {
        ez->nword = (*nword);
    } else {
        bd = (((thre)<<1)+1); ez->nword = ((bd>>bitw)+(!!(bd&bitz)));
    }
    resize_bit_extz_t((*ez), (thre)); 
	wz = ws = sizeof(*(ez->a))*(ez->nword);//diff
	memset((ez->Peq[0]).a, 0, wz);
	memset((ez->Peq[1]).a, 0, wz);
	memset((ez->Peq[2]).a, 0, wz);
	memset((ez->Peq[3]).a, 0, wz);
	memset((ez->Peq[4]).a, 0, wz);
	memset((ez->VP).a, 0, wz);
	memset((ez->VN).a, 0, wz);

	bd = (thre<<1)+1; bd = ((bd<=pn)?bd:pn); i_bd = 0;
    ed_init_core(i, bd, i_bd, pstr, ez->Peq);
    bd = (thre<<1)+1, i_bd = (thre<<1);

	memset((ez->Peq[4]).a, 0, wz);

	ez->path.n=(ez->nword*tn*5);//diff
    kv_resize(w_sig, ez->path, ez->path.n); ez->path.n=0;//diff

	i = 0; 
    /**for the incoming char/last char; mm = ((Word)1 << (thre<<1))**/
    Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));
	while (i < tn0) {
        ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c, ad, wz, ez->nword);
        if (!(ez->D0.a[0]&(1ULL))) {
            ++err; if (err>cut) return;
        }
        
        ed_infi_post_Peq(ez->Peq, wz, wz1, ez->nword);
        ++i; ++i_bd; c = 4;
        if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];
        if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;

		memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;//diff
        memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;//diff
        memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;//diff
        memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;//diff
        memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;//diff
    }
	ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c, ad, wz, ez->nword);
    if (!(ez->D0.a[0]&(1ULL))) {
        ++err; if (err>cut) return;
    }
	memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;//diff

	int32_t site = tn - 1;/**up bound**/
	/**in most cases, ai = (thre<<1)**/
    int32_t ai = pn - tn, uge = INT32_MAX;
	if(ez->err > thre) {//diff
		if ((err <= thre) && (err <= ez->err)) {
			ez->err = err; ez->pe = site;
		}

		i = 0;
		while (i < ai) {
			bd = (i>>bitw); i_bd = (i&bitz);
			err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));
			err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));
			++i;
			if ((err <= thre) && (err <= ez->err)) {
				ez->err = err; ez->pe = site + i;
			}
			if(i == thre) uge = err;
		}
		if((uge <= thre) && (uge == ez->err)) ez->pe = site + thre;
	}
	gen_trace(ez, 0, 1);//diff
}

inline void ed_band_cal_extension_infi_0_w_trace(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, int32_t *nword, bit_extz_t *ez)
{
	int32_t done = 1;//diff
	ez->cigar.n = 0; //diff
	if(ez->err > thre) {//diff
		init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = 0; ez->pe = ez->te = -1;
		if(pn > tn + thre) pn = tn + thre;
		else if(tn > pn + thre) tn = pn + thre;
		done = 0;//diff
	} else if(ez->err == 0) {//diff
        push_trace(&(ez->cigar), 0, ez->te+1-ez->ts); return;//diff
    } else {
        pn = ez->pe + 1 - ez->ps; tn = ez->te + 1 - ez->ts; //diff
    }
    int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd, i_bd, k_bd, wz, wz1, ws, Peq_i; w_sig c, ad, Peq_m;
    int32_t poff, pe = pn-1, tmp_e = INT32_MAX, k;
	if(nword) {
        ez->nword = (*nword);
    } else {
        bd = (((thre)<<1)+1); ez->nword = ((bd>>bitw)+(!!(bd&bitz)));
    }
    resize_bit_extz_t((*ez), (thre)); 
	wz = ws = sizeof(*(ez->a))*(ez->nword);//diff
    memset((ez->Peq[0]).a, 0, wz);
    memset((ez->Peq[1]).a, 0, wz);
    memset((ez->Peq[2]).a, 0, wz);
    memset((ez->Peq[3]).a, 0, wz);
    memset((ez->Peq[4]).a, 0, wz);

	bd = thre+1; bd = ((bd<=pn)?bd:pn); i_bd = thre;
    ed_init_core(i, bd, i_bd, pstr, ez->Peq);
    bd = thre+1, i_bd = thre;

	memset((ez->Peq[4]).a, 0, wz);
    err = thre;
	w_infi_set_bit_lsub(ez->VN, thre, ez->nword); /**VN = (((Word)1)<<(thre))-1;**/ 
    w_infi_set_bit_lsub(ez->VP, (thre<<1)+1, ez->nword); /**VP = (((Word)1)<<((thre<<1)+1))-1;**/
    w_infi_self_xor(ez->VP, ez->VN, ez->nword, wz) /**VP ^= VN;**/

	ez->path.n=(ez->nword*tn*5);//diff
    kv_resize(w_sig, ez->path, ez->path.n); ez->path.n=0;//diff

	i = 0;
	/**for the incoming char/last char; mm = ((Word)1 << (thre<<1))**/
    Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));
	//for debug
	// memset((ez->X).a, 0, sizeof(*(ez->a))*(ez->nword));
	// memset((ez->D0).a, 0, sizeof(*(ez->a))*(ez->nword));
	// memset((ez->HN).a, 0, sizeof(*(ez->a))*(ez->nword));
	// memset((ez->HP).a, 0, sizeof(*(ez->a))*(ez->nword));
	while (i < tn0) {
		// fprintf(stderr, "i::%d\n", i);
		// prt_bit_extz_t((*ez), ((thre<<1)+1));
        ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c, ad, wz, ez->nword);
        if (!(ez->D0.a[0]&(1ULL))) {
            ++err; if (err>cut) return;
        }

		if(!done) {//diff
            poff = i-thre; k = i+thre-pe;/**poff:[i-thre, i+thre]**/
			if(k >= 0) {
				if(tmp_e == INT32_MAX) {
					tmp_e = err;
					for ((k) = 0; (poff) < (pe); (poff)++, (k)++) {
						bd = (k>>bitw); k_bd = (k&bitz);
						(tmp_e) += (((*ez).VP.a[bd]>>k_bd)&((w_sig)1)); 
						(tmp_e) -= (((*ez).VN.a[bd]>>k_bd)&((w_sig)1));
					} 
				} else {
					k = (thre<<1) - k; 
					if(k >= 0) {
						bd = (k>>bitw); k_bd = (k&bitz);
						(tmp_e) += (((*ez).HP.a[bd]>>k_bd)&((w_sig)1)); 
						(tmp_e) -= (((*ez).HN.a[bd]>>k_bd)&((w_sig)1));
					}
				}
				if((tmp_e) <= (*ez).thre && (tmp_e) < (*ez).err) {
                    (*ez).err = tmp_e; (*ez).pe = pe; (*ez).te = i;
				}
			}
			// if(dbg_ext_err(i, thre, err, pe, &((*ez).VP), &((*ez).VN))!=tmp_e) {
			// 	fprintf(stderr, "i::%d, tmp_e::%d, dbg_ext_err::%d\n", i, tmp_e, dbg_ext_err(i, thre, err, pe, &((*ez).VP), &((*ez).VN)));
			// }
        }

        ed_infi_post_Peq(ez->Peq, wz, wz1, ez->nword);
        ++i; ++i_bd; c = 4;
        if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];
        if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;

		memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;//diff
        memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;//diff
        memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;//diff
        memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;//diff
        memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;//diff
    }
    ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c, ad, wz, ez->nword);
    if (!(ez->D0.a[0]&(1ULL))) {
        ++err; if (err>cut) return;
    }
	memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;//diff

    int32_t site = tn - 1 - thre;/**up bound; site:[tn - 1 - thre, tn - 1 + thre]**/
	if(!done) {//diff
		for (cut = pn - 1, i = 0; site < cut; i++) {
			bd = (i>>bitw); i_bd = (i&bitz);
			err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));
			err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));
			site++;
			if(err <= thre && err < ez->err) {
				ez->err = err; ez->pe = site; ez->te = tn-1;
			}
		}
		if(err <= thre && err < ez->err) {
			ez->err = err; ez->pe = site; ez->te = tn-1;
		}
	}

	if((ez->te-ez->ts+1) != tn) {//diff
        ez->path.n /= tn; ez->path.n *= (ez->te+1-ez->ts);//diff
    }//diff
    gen_trace(ez, thre, 1);//diff
    return;
}

inline void ed_band_cal_extension_infi_1_w_trace(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, int32_t *nword, bit_extz_t *ez)
{
	int32_t done = 1;//diff
    ez->cigar.n = 0;//diff
	if(ez->err > thre) {//diff
		init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = INT32_MAX; ez->pe = pn-1; ez->te = tn-1;
		if(pn > tn + thre) pn = tn + thre;
		else if(tn > pn + thre) tn = pn + thre;
		done = 0;//diff
    } else if(ez->err == 0) {//diff
        push_trace(&(ez->cigar), 0, ez->te+1-ez->ts); return;//diff
    } else {
        pn = ez->pe + 1 - ez->ps; tn = ez->te + 1 - ez->ts; //diff
    }
	int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd, i_bd, k_bd, wz, wz1, ws, Peq_i; 
    int32_t poff, pe = pn-1, tmp_e = INT32_MAX, k, pidx = ez->pe, tidx = ez->te; w_sig c, ad, Peq_m;
	if(nword) {
        ez->nword = (*nword);
    } else {
        bd = (((thre)<<1)+1); ez->nword = ((bd>>bitw)+(!!(bd&bitz)));
    }
    resize_bit_extz_t((*ez), (thre)); 
	wz = ws = sizeof(*(ez->a))*(ez->nword);//diff
    memset((ez->Peq[0]).a, 0, wz);
    memset((ez->Peq[1]).a, 0, wz);
    memset((ez->Peq[2]).a, 0, wz);
    memset((ez->Peq[3]).a, 0, wz);
    memset((ez->Peq[4]).a, 0, wz);

	bd = thre+1; bd = ((bd<=pn)?bd:pn); i_bd = thre;
	for (i = 0; i < bd; i++, i_bd++) { 
		w_bit((ez->Peq)[seq_nt4_table[(uint8_t)(pstr)[pidx-i]]], (i_bd));
	}
    bd = thre+1, i_bd = thre;

	memset((ez->Peq[4]).a, 0, wz);
    err = thre;
	w_infi_set_bit_lsub(ez->VN, thre, ez->nword); /**VN = (((Word)1)<<(thre))-1;**/ 
    w_infi_set_bit_lsub(ez->VP, (thre<<1)+1, ez->nword); /**VP = (((Word)1)<<((thre<<1)+1))-1;**/
    w_infi_self_xor(ez->VP, ez->VN, ez->nword, wz) /**VP ^= VN;**/

	ez->path.n=(ez->nword*tn*5);//diff
    kv_resize(w_sig, ez->path, ez->path.n); ez->path.n=0;//diff

	i = 0;
	/**for the incoming char/last char; mm = ((Word)1 << (thre<<1))**/
    Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));
	while (i < tn0) {
        ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[tidx-i], c, ad, wz, ez->nword);
        if (!(ez->D0.a[0]&(1ULL))) {
            ++err; if (err>cut) return;
        }

		if(!done) {//diff
            poff = i-thre; k = i+thre-pe;/**poff:[i-thre, i+thre]**/
			if(k >= 0) {
				if(tmp_e == INT32_MAX) {
					tmp_e = err;
					for ((k) = 0; (poff) < (pe); (poff)++, (k)++) {
						bd = (k>>bitw); k_bd = (k&bitz);
						(tmp_e) += (((*ez).VP.a[bd]>>k_bd)&((w_sig)1)); 
						(tmp_e) -= (((*ez).VN.a[bd]>>k_bd)&((w_sig)1));
					} 
				} else {
					k = (thre<<1) - k; 
					if(k >= 0) {
						bd = (k>>bitw); k_bd = (k&bitz);
						(tmp_e) += (((*ez).HP.a[bd]>>k_bd)&((w_sig)1)); 
						(tmp_e) -= (((*ez).HN.a[bd]>>k_bd)&((w_sig)1));
					}
				}
				if((tmp_e) <= (*ez).thre && (tmp_e) < (*ez).err) {
                    (*ez).err = tmp_e; (*ez).ps = pidx - pe; (*ez).ts = tidx-i;
				}
			}
        }
		// if(dbg_ext_err(i, thre, err, pe, &((*ez).VP), &((*ez).VN))!=tmp_e) {
		// 	fprintf(stderr, "i::%d, tmp_e::%d, dbg_ext_err::%d\n", i, tmp_e, dbg_ext_err(i, thre, err, pe, &((*ez).VP), &((*ez).VN)));
		// }

        ed_infi_post_Peq(ez->Peq, wz, wz1, ez->nword);
        ++i; ++i_bd; c = 4;
        if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[pidx-i_bd]];
        if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;

		memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;//diff
        memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;//diff
        memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;//diff
        memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;//diff
        memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;//diff
    }
    ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[tidx-i], c, ad, wz, ez->nword);
    if (!(ez->D0.a[0]&(1ULL))) {
        ++err; if (err>cut) return;
    }
	memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;//diff

    int32_t site = tn - 1 - thre;/**up bound; site:[tn - 1 - thre, tn - 1 + thre]**/
	if(!done) {//diff
		for (cut = pn - 1, i = 0; site < cut; i++) {
			bd = (i>>bitw); i_bd = (i&bitz);
			err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));
			err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));
			site++;
			if(err <= thre && err < ez->err) {
				ez->err = err; ez->ps = pidx-site; ez->ts = tidx+1-tn;
			}
		}
		if(err <= thre && err < ez->err) {
			ez->err = err; ez->ps = pidx-site; ez->ts = tidx+1-tn;
		}
	}

    if((ez->te-ez->ts+1) != tn) {//diff
        ez->path.n /= tn; ez->path.n *= (ez->te+1-ez->ts);//diff
    }//diff

    poff = ez->ps; ez->ps = pidx - ez->pe; ez->pe = pidx - poff;//diff
    gen_trace(ez, thre, 0);//diff
    poff = ez->ps; ez->ps = pidx - ez->pe; ez->pe = pidx - poff;//diff
    return;
}

inline void ed_band_cal_semi_infi_w_absent_diag(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, int32_t abs_diag, int32_t *nword, bit_extz_t *ez)
{
	init_base_ed(*ez, thre, pn, tn); ez->ps = ez->pe = -1; ez->ts = 0; ez->te = tn-1;
	int32_t bd, i, err = abs_diag, i_bd, cut = thre+(thre<<1), tn0 = tn - 1, wz, wz1, Peq_i; w_sig c, ad, Peq_m;
	if((pn > tn + cut) || (tn > pn + cut)) return;

	if(nword) {
        ez->nword = (*nword);
    } else {
        bd = (((thre)<<1)+1); ez->nword = ((bd>>bitw)+(!!(bd&bitz)));
    }
    resize_bit_extz_t((*ez), (thre)); 
	wz = sizeof(*(ez->a))*(ez->nword);
	memset((ez->Peq[0]).a, 0, wz);
	memset((ez->Peq[1]).a, 0, wz);
	memset((ez->Peq[2]).a, 0, wz);
	memset((ez->Peq[3]).a, 0, wz);
	memset((ez->Peq[4]).a, 0, wz);
	memset((ez->VP).a, 0, wz);
	w_infi_set_bit_lsub(ez->VN, abs_diag, ez->nword); /**VN = (((Word)1)<<(abs_diag))-1; ;**/ 

	bd = ((thre<<1)+1)-abs_diag; bd = ((bd<=pn)?bd:pn); i_bd = abs_diag; 
    ed_init_core(i, bd, i_bd, pstr, ez->Peq);
    i_bd = (thre<<1)-abs_diag;

	memset((ez->Peq[4]).a, 0, wz);
	i = 0; 
    /**for the incoming char/last char; mm = ((Word)1 << (thre<<1))**/
    Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));
	while (i < tn0) {
        ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c, ad, wz, ez->nword);
        if (!(ez->D0.a[0]&(1ULL))) {
            ++err; if (err>cut) return;
        }
        
        ed_infi_post_Peq(ez->Peq, wz, wz1, ez->nword);
        ++i; ++i_bd; c = 4;
        if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];
        if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;
    }
	ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c, ad, wz, ez->nword);
    if (!(ez->D0.a[0]&(1ULL))) {
        ++err; if (err>cut) return;
    }

	int32_t site = tn - 1 - abs_diag;/**up bound**/
    /**in most cases, ai = (thre<<1)**/
    int32_t ai = pn - tn + abs_diag, uge = INT32_MAX; i = 0;
	for (i = 0; site < 0 && i < ai; i++, site++) {
		bd = (i>>bitw); i_bd = (i&bitz);
		err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));
		err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));
	}
    if ((err <= thre) && (err <= ez->err)) {
        ez->err = err; ez->pe = site;
    }
	site -= i;
	while (i < ai) {
		bd = (i>>bitw); i_bd = (i&bitz);
		err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));
		err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));
        ++i;
        if ((err <= thre) && (err <= ez->err)) {
            ez->err = err; ez->pe = site + i;
        }
        if(i == thre) uge = err;
    }
    if((uge <= thre) && (uge == ez->err)) ez->pe = site + thre;
}

inline void ed_band_cal_semi_infi_w_absent_diag_trace(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, int32_t abs_diag, int32_t *nword, bit_extz_t *ez)
{
	ez->cigar.n = 0;//diff
	if(ez->err > thre) {//diff
		init_base_ed(*ez, thre, pn, tn); ez->ps = ez->pe = -1; ez->ts = 0; ez->te = tn-1;
	} else if(ez->err == 0) {//diff
        push_trace(&(ez->cigar), 0, ez->te+1-ez->ts); //diff
        ez->ps = ez->pe - (ez->te-ez->ts);//diff
        return;//diff
    }
	int32_t bd, i, err = abs_diag, i_bd, cut = thre+(thre<<1), tn0 = tn - 1, wz, wz1, ws, Peq_i; w_sig c, ad, Peq_m;
	if((pn > tn + cut) || (tn > pn + cut)) return;

	if(nword) {
        ez->nword = (*nword);
    } else {
        bd = (((thre)<<1)+1); ez->nword = ((bd>>bitw)+(!!(bd&bitz)));
    }
    resize_bit_extz_t((*ez), (thre)); 
	wz = ws = sizeof(*(ez->a))*(ez->nword);//diff
	memset((ez->Peq[0]).a, 0, wz);
	memset((ez->Peq[1]).a, 0, wz);
	memset((ez->Peq[2]).a, 0, wz);
	memset((ez->Peq[3]).a, 0, wz);
	memset((ez->Peq[4]).a, 0, wz);
	memset((ez->VP).a, 0, wz);
	w_infi_set_bit_lsub(ez->VN, abs_diag, ez->nword); /**VN = (((Word)1)<<(abs_diag))-1; ;**/ 

	bd = ((thre<<1)+1)-abs_diag; bd = ((bd<=pn)?bd:pn); i_bd = abs_diag; 
    ed_init_core(i, bd, i_bd, pstr, ez->Peq);
    i_bd = (thre<<1)-abs_diag;

	memset((ez->Peq[4]).a, 0, wz);

	ez->path.n=(ez->nword*tn*5);//diff
    kv_resize(w_sig, ez->path, ez->path.n); ez->path.n=0;//diff

	i = 0; 
    /**for the incoming char/last char; mm = ((Word)1 << (thre<<1))**/
    Peq_i = (((thre<<1))>>bitw); Peq_m = (((w_sig)1)<<(((thre<<1))&bitz));
	while (i < tn0) {
        ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c, ad, wz, ez->nword);
        if (!(ez->D0.a[0]&(1ULL))) {
            ++err; if (err>cut) return;
        }
        
        ed_infi_post_Peq(ez->Peq, wz, wz1, ez->nword);
        ++i; ++i_bd; c = 4;
        if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];
        if(c < 4) ez->Peq[c].a[Peq_i]|=Peq_m;

		memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;//diff
        memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;//diff
        memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;//diff
        memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;//diff
        memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;//diff
    }
	ed_infi_core(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP, (uint8_t)tstr[i], c, ad, wz, ez->nword);
    if (!(ez->D0.a[0]&(1ULL))) {
        ++err; if (err>cut) return;
    }
	memcpy(ez->path.a+ez->path.n, ez->D0.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->VP.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->VN.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->HP.a, ws); ez->path.n += ez->nword;//diff
	memcpy(ez->path.a+ez->path.n, ez->HN.a, ws); ez->path.n += ez->nword;//diff

	int32_t site = tn - 1 - abs_diag;/**up bound**/
    /**in most cases, ai = (thre<<1)**/
    int32_t ai = pn - tn + abs_diag, uge = INT32_MAX; i = 0;
	if(ez->err > thre) {//diff
		for (i = 0; site < 0 && i < ai; i++, site++) {
			bd = (i>>bitw); i_bd = (i&bitz);
			err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));
			err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));
		}
		if ((err <= thre) && (err <= ez->err)) {
			ez->err = err; ez->pe = site;
		}
		site -= i;
		while (i < ai) {
			bd = (i>>bitw); i_bd = (i&bitz);
			err += (((*ez).VP.a[bd]>>i_bd)&((w_sig)1));
			err -= (((*ez).VN.a[bd]>>i_bd)&((w_sig)1));
			++i;
			if ((err <= thre) && (err <= ez->err)) {
				ez->err = err; ez->pe = site + i;
			}
			if(i == thre) uge = err;
		}
		if((uge <= thre) && (uge == ez->err)) ez->pe = site + thre;
	}
	gen_trace(ez, abs_diag, 1);//diff
}


#define ed_core_64(Peq, VP, VN, X, D0, HN, HP, z) {	\
		/**X = Peq[seq_nt4_table[(uint8_t)tstr[i]]] | VN;**/\
		(X) = (Peq)[(seq_nt4_table[(z)])]|(VN);\
		(D0) = (((VP) + ((X)&(VP))) ^ (VP)) | (X);\
        (HN) = (VP)&(D0);\
        (HP) = (VN) | ~((VP) | (D0));\
        (X) = (D0) >> 1;\
        (VN) = (X)&(HP);\
        (VP) = (HN) | ~((X) | (HP));\
}

inline void ed_band_cal_global_64_w(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{
	init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = 0;
    if((pn > tn + thre) || (tn > pn + thre)) return;
	int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd, i_bd; 
	Word c, Peq[5] = {0}, VP, VN, X, D0, HN, HP, mm;

	bd = thre+1; bd = ((bd<=pn)?bd:pn); i_bd = thre;
	for (i = 0, mm = (((Word)1)<<i_bd); i < bd; i++) {
        Peq[seq_nt4_table[(uint8_t)pstr[i]]] |= mm; mm <<= 1;
    }	
	bd = thre+1, i_bd = thre; Peq[4] = 0; err = thre;
	VN = (((Word)1)<<(thre))-1; VP = (((Word)1)<<((thre<<1)+1))-1; VP ^= VN;

	i = 0; mm = ((Word)1 << (thre<<1));///for the incoming char/last char
	while (i < tn0) {
		// fprintf(stderr, "\ni::%d\n", i);
		// prt_vector(Peq, VP, VN, X, D0, HN, HP, (((thre<<1))+1));
		ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[i]);
		if (!(D0&(1ULL))) {
			++err; if (err>cut) return;
		}

		Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;
		++i; ++i_bd;
		if(i_bd < pn) {
			c = seq_nt4_table[(uint8_t)pstr[i_bd]];
			if(c < 4) Peq[c] |= mm;
			// fprintf(stderr, "c::%lu\n", c);
			// print_bit(Peq[c], (((thre<<1))+1), "Peq[c]");
		}
	}
	ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[i]);
	if (!(D0&(1ULL))) {
		++err; if (err>cut) return;
	}

	int32_t site = tn - 1 - thre;/**up bound**/\
    for (cut = pn - 1; site < cut; site++) {
        err += VP&(1ULL); VP >>= 1; err -= VN&(1ULL); VN >>= 1; 
    }

    if (site == cut && err <= thre) {
        ez->err = err; 
        ez->pe = pn-1; ez->te = tn-1;
    }
    return;
}

inline void ed_band_cal_semi_64_w(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{
	init_base_ed(*ez, thre, pn, tn); ez->ps = ez->pe = -1; ez->ts = 0; ez->te = tn-1;
	Word c, Peq[5] = {0}, VP = 0, VN = 0, X, D0, HN, HP, mm;
	int32_t bd, i, err = 0, i_bd, last_high = (thre<<1), tn0 = tn - 1, cut = thre+last_high;
	if((pn > tn + cut) || (tn > pn + cut)) return;

	bd = (thre<<1)+1; bd = ((bd<=pn)?bd:pn); 
	for (i = 0, mm = 1; i < bd; i++) {
        Peq[seq_nt4_table[(uint8_t)pstr[i]]] |= mm; mm <<= 1;
    }
	bd = (thre<<1)+1, i_bd = (thre<<1);

	i = 0; Peq[4] = 0; mm = ((Word)1 << (thre<<1));///for the incoming char/last char**
	while (i < tn0) {
		ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[i]);
		if (!(D0&(1ULL))) {
            ++err; if (err>cut) return;
        }

		Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;
		++i; ++i_bd; c = 4;
		if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];
		if(c < 4) Peq[c] |= mm; 
	}
	ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[i]);
	if (!(D0&(1ULL))) {
		++err; if (err>cut) return;
	}

	int32_t site = tn - 1;/**up bound**/
    /**in most cases, ai = (thre<<1)**/
    int32_t ai = pn - tn, uge = INT32_MAX;
    if ((err <= thre) && (err <= ez->err)) {
        ez->err = err; ez->pe = site;
    }

    i = 0;
    while (i < ai) {\
        err += ((VP >> i)&(1ULL)); err -= ((VN >> i)&(1ULL)); ++i;
        if ((err <= thre) && (err <= ez->err)) {
            ez->err = err; ez->pe = site + i;
        }
        if(i == thre) uge = err;
    }
    
    if((uge <= thre) && (uge == ez->err)) ez->pe = site + thre;
}

inline void ed_band_cal_extension_64_0_w(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{
	init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = 0; ez->pe = ez->te = -1;
	if(pn > tn + thre) pn = tn + thre;
    else if(tn > pn + thre) tn = pn + thre;
	int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd, i_bd, poff, pe = pn-1, tmp_e = INT32_MAX, k;
	Word c, Peq[5] = {0}, VP, VN, X, D0, HN, HP, mm;

	bd = thre+1; bd = ((bd<=pn)?bd:pn); i_bd = thre;
	for (i = 0, mm = (((Word)1)<<i_bd); i < bd; i++) {
        Peq[seq_nt4_table[(uint8_t)pstr[i]]] |= mm; mm <<= 1;
    }	
	bd = thre+1, i_bd = thre; Peq[4] = 0; err = thre;
	VN = (((Word)1)<<(thre))-1; VP = (((Word)1)<<((thre<<1)+1))-1; VP ^= VN;

	i = 0; mm = ((Word)1 << (thre<<1));
	// X = D0 = HN = HP = mm = 0;//for debug
	while (i < tn0) {
		// fprintf(stderr, "i::%d\n", i);
		// prt_vector(Peq, VP, VN, X, D0, HN, HP, ((thre<<1)+1));
		ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[i]);
		if (!(D0&(1ULL))) {
			++err; if (err>cut) return;
		}

		{
			poff = i-thre; k = i+thre-pe;/**poff:[i-thre, i+thre]**/
            if(k >= 0) {
                if(tmp_e == INT32_MAX) {
                    tmp_e = err;
                    for (k = 0; poff < pe; poff++, (k)++) {
                        tmp_e += ((VP>>k)&(1ULL)); tmp_e -= ((VN>>k)&(1ULL));
                    } 
                } else {
                    k = (thre<<1) - k; 
                    if(k >= 0) {
                        tmp_e += ((HP>>k)&(1ULL)); tmp_e -= ((HN>>k)&(1ULL)); 
                    }
                }
                if(tmp_e <= (*ez).thre && tmp_e < (*ez).err) {
					(*ez).err = tmp_e; (*ez).pe = pe; (*ez).te = i;
                }
            }
        }
		// fprintf(stderr, "i::%d, tmp_e::%d, err::%d\n", i, tmp_e, err);

		Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;
		++i; ++i_bd;
        if(i_bd < pn) {
            c = seq_nt4_table[(uint8_t)pstr[i_bd]]; 
            if(c < 4) Peq[c] |= mm;
        }
	}
	ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[i]);
	if (!(D0&(1ULL))) {
		++err; if (err>cut) return;
	}

	int32_t site = tn - 1 - thre;/**up bound; site:[tn - 1 - thre, tn - 1 + thre]**/\
    for (cut = pn - 1; site < cut; ) {
		err += VP&(1ULL); VP >>= 1; err -= VN&(1ULL); VN >>= 1; 
        site++;
        if(err <= thre && err < ez->err) {
            ez->err = err; ez->pe = site; ez->te = tn-1;
        }
    }
    if(err <= thre && err < ez->err) {
        ez->err = err; ez->pe = site; ez->te = tn-1;
    }
    return;
}

inline void ed_band_cal_extension_64_1_w(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{
	init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = INT32_MAX; ez->pe = pn-1; ez->te = tn-1;
    if(pn > tn + thre) pn = tn + thre;
    else if(tn > pn + thre) tn = pn + thre;
	int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd, i_bd;
	int32_t poff, pe = pn-1, tmp_e = INT32_MAX, k, pidx = ez->pe, tidx = ez->te;
	Word c, Peq[5] = {0}, VP, VN, X, D0, HN, HP, mm;

	bd = thre+1; bd = ((bd<=pn)?bd:pn); i_bd = thre;
	for (i = 0, mm = (((Word)1)<<i_bd); i < bd; i++) {
        Peq[seq_nt4_table[(uint8_t)pstr[pidx-i]]] |= mm; mm <<= 1;
    }	
	bd = thre+1, i_bd = thre; Peq[4] = 0; err = thre;
	VN = (((Word)1)<<(thre))-1; VP = (((Word)1)<<((thre<<1)+1))-1; VP ^= VN;

	i = 0; mm = ((Word)1 << (thre<<1));
	// X = D0 = HN = HP = mm = 0;//for debug
	while (i < tn0) {
		// fprintf(stderr, "i::%d\n", i);
		// prt_vector(Peq, VP, VN, X, D0, HN, HP, ((thre<<1)+1));
		ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[tidx-i]);
		if (!(D0&(1ULL))) {
			++err; if (err>cut) return;
		}

		{
			poff = i-thre; k = i+thre-pe;/**poff:[i-thre, i+thre]**/
            if(k >= 0) {
                if(tmp_e == INT32_MAX) {
                    tmp_e = err;
                    for (k = 0; poff < pe; poff++, (k)++) {
                        tmp_e += ((VP>>k)&(1ULL)); tmp_e -= ((VN>>k)&(1ULL));
                    } 
                } else {
                    k = (thre<<1) - k; 
                    if(k >= 0) {
                        tmp_e += ((HP>>k)&(1ULL)); tmp_e -= ((HN>>k)&(1ULL)); 
                    }
                }
                if(tmp_e <= (*ez).thre && tmp_e < (*ez).err) {
					(*ez).err = tmp_e; (*ez).ps = pidx - pe; (*ez).ts = tidx-i;
                }
            }
        }
		// fprintf(stderr, "i::%d, tmp_e::%d, err::%d\n", i, tmp_e, err);

		Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;
		++i; ++i_bd;
        if(i_bd < pn) {
            c = seq_nt4_table[(uint8_t)pstr[pidx-i_bd]]; 
            if(c < 4) Peq[c] |= mm;
        }
	}
	ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[tidx-i]);
	if (!(D0&(1ULL))) {
		++err; if (err>cut) return;
	}

	int32_t site = tn - 1 - thre;/**up bound; site:[tn - 1 - thre, tn - 1 + thre]**/\
    for (cut = pn - 1; site < cut; ) {
		err += VP&(1ULL); VP >>= 1; err -= VN&(1ULL); VN >>= 1; 
        site++;
        if(err <= thre && err < ez->err) {
            ez->err = err; ez->ps = pidx-site; ez->ts = tidx+1-tn;
        }
    }
    if(err <= thre && err < ez->err) {
        ez->err = err; ez->ps = pidx-site; ez->ts = tidx+1-tn;
    }
    return;
}

inline void ed_band_cal_global_64_w_trace(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{	
	ez->cigar.n = 0; ez->nword = 1;//diff
	if(ez->err > thre) {//diff
		init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = 0;//diff
	} else if(ez->err == 0) {//diff
		push_trace(&(ez->cigar), 0, ez->te+1-ez->ts); return;//diff
	}
	if((pn > tn + thre) || (tn > pn + thre)) return;
	int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd, i_bd; 
	w_sig c, Peq[5] = {0}, VP, VN, X, D0, HN, HP, mm;

	bd = thre+1; bd = ((bd<=pn)?bd:pn); i_bd = thre;
	for (i = 0, mm = (((w_sig)1)<<i_bd); i < bd; i++) {
        Peq[seq_nt4_table[(uint8_t)pstr[i]]] |= mm; mm <<= 1;
    }	
	bd = thre+1, i_bd = thre; Peq[4] = 0; err = thre;
	VN = (((w_sig)1)<<(thre))-1; VP = (((w_sig)1)<<((thre<<1)+1))-1; VP ^= VN;

	ez->path.n=(ez->nword*tn*5);//diff
	kv_resize(w_sig, ez->path, ez->path.n); ez->path.n=0;//diff

	i = 0; mm = ((w_sig)1 << (thre<<1));///for the incoming char/last char
	while (i < tn0) {
		// fprintf(stderr, "\ni::%d\n", i);
		// prt_vector(Peq, VP, VN, X, D0, HN, HP, (((thre<<1))+1));
		ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[i]);
		if (!(D0&(1ULL))) {
			++err; if (err>cut) return;
		}

		Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;
		++i; ++i_bd;
		if(i_bd < pn) {
			c = seq_nt4_table[(uint8_t)pstr[i_bd]];
			if(c < 4) Peq[c] |= mm;
		}
		ez->path.a[ez->path.n++] = D0;//diff
		ez->path.a[ez->path.n++] = VP;//diff
		ez->path.a[ez->path.n++] = VN;//diff
		ez->path.a[ez->path.n++] = HP;//diff
		ez->path.a[ez->path.n++] = HN;//diff
	}
	ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[i]);
	if (!(D0&(1ULL))) {
		++err; if (err>cut) return;
	}
	ez->path.a[ez->path.n++] = D0;//diff
	ez->path.a[ez->path.n++] = VP;//diff
	ez->path.a[ez->path.n++] = VN;//diff
	ez->path.a[ez->path.n++] = HP;//diff
	ez->path.a[ez->path.n++] = HN;//diff

	int32_t site = tn - 1 - thre;/**up bound**/\
	if(ez->err > thre) {//diff
		for (cut = pn - 1; site < cut; site++) {
			err += VP&(1ULL); VP >>= 1; err -= VN&(1ULL); VN >>= 1; 
		}

		if (site == cut && err <= thre) {
			ez->err = err; 
			ez->pe = pn-1; ez->te = tn-1;
		}
	}
	///should update ez->path.n for extension
	gen_trace(ez, thre, 1);//diff
    return;
}

inline void ed_band_cal_semi_64_w_trace(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{
	ez->cigar.n = 0; ez->nword = 1;//diff
	if(ez->err > thre) {//diff
		init_base_ed(*ez, thre, pn, tn); ez->ps = ez->pe = -1; ez->ts = 0; ez->te = tn-1;//diff
	} else if(ez->err == 0) {//diff
		push_trace(&(ez->cigar), 0, ez->te+1-ez->ts); //diff
		ez->ps = ez->pe - (ez->te-ez->ts);//diff
		return;//diff
	}
	Word c, Peq[5] = {0}, VP = 0, VN = 0, X, D0, HN, HP, mm;
	int32_t bd, i, err = 0, i_bd, last_high = (thre<<1), tn0 = tn - 1, cut = thre+last_high;
	if((pn > tn + cut) || (tn > pn + cut)) return;

	bd = (thre<<1)+1; bd = ((bd<=pn)?bd:pn); 
	for (i = 0, mm = 1; i < bd; i++) {
        Peq[seq_nt4_table[(uint8_t)pstr[i]]] |= mm; mm <<= 1;
    }
	bd = (thre<<1)+1, i_bd = (thre<<1);

	ez->path.n=(ez->nword*tn*5);//diff
    kv_resize(w_sig, ez->path, ez->path.n); ez->path.n=0;//diff

	i = 0; Peq[4] = 0; mm = ((Word)1 << (thre<<1));///for the incoming char/last char**
	while (i < tn0) {
		ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[i]);
		if (!(D0&(1ULL))) {
            ++err; if (err>cut) return;
        }

		Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;
		++i; ++i_bd; c = 4;
		if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];
		if(c < 4) Peq[c] |= mm; 

		ez->path.a[ez->path.n++] = D0;//diff
        ez->path.a[ez->path.n++] = VP;//diff
        ez->path.a[ez->path.n++] = VN;//diff
        ez->path.a[ez->path.n++] = HP;//diff
        ez->path.a[ez->path.n++] = HN;//diff
	}
	ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[i]);
	if (!(D0&(1ULL))) {
		++err; if (err>cut) return;
	}
	ez->path.a[ez->path.n++] = D0;//diff
	ez->path.a[ez->path.n++] = VP;//diff
	ez->path.a[ez->path.n++] = VN;//diff
	ez->path.a[ez->path.n++] = HP;//diff
	ez->path.a[ez->path.n++] = HN;//diff

	int32_t site = tn - 1;/**up bound**/
    /**in most cases, ai = (thre<<1)**/
    int32_t ai = pn - tn, uge = INT32_MAX;
	if(ez->err > thre) {//diff
		if ((err <= thre) && (err <= ez->err)) {
			ez->err = err; ez->pe = site;
		}

		i = 0;
		while (i < ai) {
			err += ((VP >> i)&(1ULL)); err -= ((VN >> i)&(1ULL)); ++i;
			if ((err <= thre) && (err <= ez->err)) {
				ez->err = err; ez->pe = site + i;
			}
			if(i == thre) uge = err;
		}
		if((uge <= thre) && (uge == ez->err)) ez->pe = site + thre;
	}
    ///should update ez->path.n for extension
    gen_trace(ez, 0, 1);//diff
}

inline void ed_band_cal_extension_64_0_w_trace(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{
	int32_t done = 1;//diff
	ez->cigar.n = 0; ez->nword = 1;//diff
    if(ez->err > thre) {//diff
		init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = 0; ez->pe = ez->te = -1;
		if(pn > tn + thre) pn = tn + thre;
    	else if(tn > pn + thre) tn = pn + thre;
		done = 0;//diff
	} else if(ez->err == 0) {//diff
        push_trace(&(ez->cigar), 0, ez->te+1-ez->ts); return;//diff
    } else {
		pn = ez->pe + 1 - ez->ps; tn = ez->te + 1 - ez->ts; //diff
	}
	int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd, i_bd, poff, pe = pn-1, tmp_e = INT32_MAX, k;
	Word c, Peq[5] = {0}, VP, VN, X, D0, HN, HP, mm;

	bd = thre+1; bd = ((bd<=pn)?bd:pn); i_bd = thre;
	for (i = 0, mm = (((Word)1)<<i_bd); i < bd; i++) {
        Peq[seq_nt4_table[(uint8_t)pstr[i]]] |= mm; mm <<= 1;
    }	
	bd = thre+1, i_bd = thre; Peq[4] = 0; err = thre;
	VN = (((Word)1)<<(thre))-1; VP = (((Word)1)<<((thre<<1)+1))-1; VP ^= VN;

	ez->path.n=(ez->nword*tn*5);//diff
    kv_resize(w_sig, ez->path, ez->path.n); ez->path.n=0;//diff

	i = 0; mm = ((Word)1 << (thre<<1));
	// X = D0 = HN = HP = mm = 0;//for debug
	while (i < tn0) {
		// fprintf(stderr, "i::%d\n", i);
		// prt_vector(Peq, VP, VN, X, D0, HN, HP, ((thre<<1)+1));
		ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[i]);
		if (!(D0&(1ULL))) {
			++err; if (err>cut) return;
		}

		if(!done) {//diff
			poff = i-thre; k = i+thre-pe;/**poff:[i-thre, i+thre]**/
            if(k >= 0) {
                if(tmp_e == INT32_MAX) {
                    tmp_e = err;
                    for (k = 0; poff < pe; poff++, (k)++) {
                        tmp_e += ((VP>>k)&(1ULL)); tmp_e -= ((VN>>k)&(1ULL));
                    } 
                } else {
                    k = (thre<<1) - k; 
                    if(k >= 0) {
                        tmp_e += ((HP>>k)&(1ULL)); tmp_e -= ((HN>>k)&(1ULL)); 
                    }
                }
                if(tmp_e <= (*ez).thre && tmp_e < (*ez).err) {
					(*ez).err = tmp_e; (*ez).pe = pe; (*ez).te = i;
                }
				// fprintf(stderr, "i::%d, (*ez).err::%d, (*ez).pe::%d, (*ez).te::%d, thre::%d, err::%d, tmp_e::%d, poff::%d, k::%d, pe::%d\n", 
				// i, (*ez).err, (*ez).pe, (*ez).te, thre, err, tmp_e, i-thre, i+thre-pe, pe);
            }
        }

		Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;
		++i; ++i_bd;
        if(i_bd < pn) {
            c = seq_nt4_table[(uint8_t)pstr[i_bd]]; 
            if(c < 4) Peq[c] |= mm;
        }
		// if(((ez->path.n+5)>(ez->nword*tn*5))||((ez->path.n+5)>(ez->path.m))) {
		// 	fprintf(stderr, "[M::%s::] pn::%d, tn::%d\n", __func__, pn, tn);
		// }
		ez->path.a[ez->path.n++] = D0;//diff
        ez->path.a[ez->path.n++] = VP;//diff
        ez->path.a[ez->path.n++] = VN;//diff
        ez->path.a[ez->path.n++] = HP;//diff
        ez->path.a[ez->path.n++] = HN;//diff
	}
	ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[i]);
	if (!(D0&(1ULL))) {
		++err; if (err>cut) return;
	}
	// if(((ez->path.n+5)>(ez->nword*tn*5))||((ez->path.n+5)>(ez->path.m))) {
	// 	fprintf(stderr, "[M::%s::] pn::%d, tn::%d\n", __func__, pn, tn);
	// }
	ez->path.a[ez->path.n++] = D0;//diff
	ez->path.a[ez->path.n++] = VP;//diff
	ez->path.a[ez->path.n++] = VN;//diff
	ez->path.a[ez->path.n++] = HP;//diff
	ez->path.a[ez->path.n++] = HN;//diff

	int32_t site = tn - 1 - thre;/**up bound; site:[tn - 1 - thre, tn - 1 + thre]**/
	if(!done) {//diff
		for (cut = pn - 1; site < cut; ) {
			err += VP&(1ULL); VP >>= 1; err -= VN&(1ULL); VN >>= 1; 
			site++;
			if(err <= thre && err < ez->err) {
				ez->err = err; ez->pe = site; ez->te = tn-1;
			}
		}
		if(err <= thre && err < ez->err) {
			ez->err = err; ez->pe = site; ez->te = tn-1;
		}
	}

	if((ez->te-ez->ts+1) != tn) {//diff
		ez->path.n /= tn; ez->path.n *= (ez->te+1-ez->ts);//diff
	}//diff
	gen_trace(ez, thre, 1);//diff
    return;
}

inline void ed_band_cal_extension_64_1_w_trace(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{
	int32_t done = 1;//diff
	ez->cigar.n = 0; ez->nword = 1;//diff
	if(ez->err > thre) {//diff
		init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = INT32_MAX; ez->pe = pn-1; ez->te = tn-1;
		if(pn > tn + thre) pn = tn + thre;
		else if(tn > pn + thre) tn = pn + thre;
		done = 0;//diff
	} else if(ez->err == 0) {//diff
        push_trace(&(ez->cigar), 0, ez->te+1-ez->ts); return;//diff
    } else {
        pn = ez->pe + 1 - ez->ps; tn = ez->te + 1 - ez->ts; //diff
    }

	int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd, i_bd;
	int32_t poff, pe = pn-1, tmp_e = INT32_MAX, k, pidx = ez->pe, tidx = ez->te;
	Word c, Peq[5] = {0}, VP, VN, X, D0, HN, HP, mm;

	bd = thre+1; bd = ((bd<=pn)?bd:pn); i_bd = thre;
	for (i = 0, mm = (((Word)1)<<i_bd); i < bd; i++) {
        Peq[seq_nt4_table[(uint8_t)pstr[pidx-i]]] |= mm; mm <<= 1;
    }	
	bd = thre+1, i_bd = thre; Peq[4] = 0; err = thre;
	VN = (((Word)1)<<(thre))-1; VP = (((Word)1)<<((thre<<1)+1))-1; VP ^= VN;

	ez->path.n=(ez->nword*tn*5);//diff
    kv_resize(w_sig, ez->path, ez->path.n); ez->path.n=0;//diff

	i = 0; mm = ((Word)1 << (thre<<1));
	// X = D0 = HN = HP = mm = 0;//for debug
	while (i < tn0) {
		// fprintf(stderr, "i::%d\n", i);
		// prt_vector(Peq, VP, VN, X, D0, HN, HP, ((thre<<1)+1));
		ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[tidx-i]);
		if (!(D0&(1ULL))) {
			++err; if (err>cut) return;
		}

		if(!done) {//diff
			poff = i-thre; k = i+thre-pe;/**poff:[i-thre, i+thre]**/
            if(k >= 0) {
                if(tmp_e == INT32_MAX) {
                    tmp_e = err;
                    for (k = 0; poff < pe; poff++, (k)++) {
                        tmp_e += ((VP>>k)&(1ULL)); tmp_e -= ((VN>>k)&(1ULL));
                    } 
                } else {
                    k = (thre<<1) - k; 
                    if(k >= 0) {
                        tmp_e += ((HP>>k)&(1ULL)); tmp_e -= ((HN>>k)&(1ULL)); 
                    }
                }
                if(tmp_e <= (*ez).thre && tmp_e < (*ez).err) {
					(*ez).err = tmp_e; (*ez).ps = pidx - pe; (*ez).ts = tidx-i;
                }
            }
        }
		// fprintf(stderr, "i::%d, tmp_e::%d, err::%d\n", i, tmp_e, err);

		Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;
		++i; ++i_bd;
        if(i_bd < pn) {
            c = seq_nt4_table[(uint8_t)pstr[pidx-i_bd]]; 
            if(c < 4) Peq[c] |= mm;
        }

		ez->path.a[ez->path.n++] = D0;//diff
        ez->path.a[ez->path.n++] = VP;//diff
        ez->path.a[ez->path.n++] = VN;//diff
        ez->path.a[ez->path.n++] = HP;//diff
        ez->path.a[ez->path.n++] = HN;//diff
	}
	ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[tidx-i]);
	if (!(D0&(1ULL))) {
		++err; if (err>cut) return;
	}
	ez->path.a[ez->path.n++] = D0;//diff
	ez->path.a[ez->path.n++] = VP;//diff
	ez->path.a[ez->path.n++] = VN;//diff
	ez->path.a[ez->path.n++] = HP;//diff
	ez->path.a[ez->path.n++] = HN;//diff

	int32_t site = tn - 1 - thre;/**up bound; site:[tn - 1 - thre, tn - 1 + thre]**/\
	if(!done) {//diff
		for (cut = pn - 1; site < cut; ) {
			err += VP&(1ULL); VP >>= 1; err -= VN&(1ULL); VN >>= 1; 
			site++;
			if(err <= thre && err < ez->err) {
				ez->err = err; ez->ps = pidx-site; ez->ts = tidx+1-tn;
			}
		}
		if(err <= thre && err < ez->err) {
			ez->err = err; ez->ps = pidx-site; ez->ts = tidx+1-tn;
		}
	}

	if((ez->te-ez->ts+1) != tn) {//diff
        ez->path.n /= tn; ez->path.n *= (ez->te+1-ez->ts);//diff
    }//diff

	poff = ez->ps; ez->ps = pidx - ez->pe; ez->pe = pidx - poff;//diff
	gen_trace(ez, thre, 0);//diff
	poff = ez->ps; ez->ps = pidx - ez->pe; ez->pe = pidx - poff;//diff
    return;
}

inline void ed_band_cal_semi_64_w_absent_diag(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, int32_t abs_diag, bit_extz_t *ez)
{
	init_base_ed(*ez, thre, pn, tn); ez->ps = ez->pe = -1; ez->ts = 0; ez->te = tn-1;
	Word c, Peq[5] = {0}, VP = 0, VN, X, D0, HN, HP, mm;
	int32_t bd, i, err = abs_diag, i_bd, last_high = (thre<<1), tn0 = tn - 1, cut = thre+last_high;
	if((pn > tn + cut) || (tn > pn + cut)) return;

	bd = ((thre<<1)+1)-abs_diag; bd = ((bd<=pn)?bd:pn); i_bd = abs_diag; 
	for (i = 0, mm = (((Word)1)<<i_bd); i < bd; i++) {
        Peq[seq_nt4_table[(uint8_t)pstr[i]]] |= mm; mm <<= 1;
    }
	i_bd = (thre<<1)-abs_diag; VN = (((Word)1)<<(abs_diag))-1; 

	i = 0; Peq[4] = 0; mm = ((Word)1 << (thre<<1));///for the incoming char/last char**
	while (i < tn0) {
		ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[i]);
		if (!(D0&(1ULL))) {
            ++err; if (err>cut) return;
        }

		Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;
		++i; ++i_bd; c = 4;
		if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];
		if(c < 4) Peq[c] |= mm; 
	}
	ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[i]);
	if (!(D0&(1ULL))) {
		++err; if (err>cut) return;
	}

	int32_t site = tn - 1 - abs_diag;/**up bound**/
    /**in most cases, ai = (thre<<1)**/
    int32_t ai = pn - tn + abs_diag, uge = INT32_MAX; i = 0;
	for (i = 0; site < 0 && i < ai; i++, site++) {
		err += ((VP >> i)&(1ULL)); err -= ((VN >> i)&(1ULL));
	}
    if ((err <= thre) && (err <= ez->err)) {
        ez->err = err; ez->pe = site;
    }
    site -= i;
    while (i < ai) {
        err += ((VP >> i)&(1ULL)); err -= ((VN >> i)&(1ULL)); ++i;
        if ((err <= thre) && (err <= ez->err)) {
            ez->err = err; ez->pe = site + i;
        }
        if(i == thre) uge = err;
    }
    
    if((uge <= thre) && (uge == ez->err)) ez->pe = site + thre;
}

inline void ed_band_cal_semi_64_w_absent_diag_trace(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, int32_t abs_diag, bit_extz_t *ez)
{
	ez->cigar.n = 0; ez->nword = 1;//diff
	if(ez->err > thre) {//diff
		init_base_ed(*ez, thre, pn, tn); ez->ps = ez->pe = -1; ez->ts = 0; ez->te = tn-1;//diff
	} else if(ez->err == 0) {//diff
		push_trace(&(ez->cigar), 0, ez->te+1-ez->ts); //diff
		ez->ps = ez->pe - (ez->te-ez->ts);//diff
		return;//diff
	}
	Word c, Peq[5] = {0}, VP = 0, VN, X, D0, HN, HP, mm;
	int32_t bd, i, err = abs_diag, i_bd, last_high = (thre<<1), tn0 = tn - 1, cut = thre+last_high;
	if((pn > tn + cut) || (tn > pn + cut)) return;

	bd = ((thre<<1)+1)-abs_diag; bd = ((bd<=pn)?bd:pn); i_bd = abs_diag; 
	for (i = 0, mm = (((Word)1)<<i_bd); i < bd; i++) {
        Peq[seq_nt4_table[(uint8_t)pstr[i]]] |= mm; mm <<= 1;
    }
	i_bd = (thre<<1)-abs_diag; VN = (((Word)1)<<(abs_diag))-1;

	ez->path.n=(ez->nword*tn*5);//diff
    kv_resize(w_sig, ez->path, ez->path.n); ez->path.n=0;//diff

	i = 0; Peq[4] = 0; mm = ((Word)1 << (thre<<1));///for the incoming char/last char**
	while (i < tn0) {
		ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[i]);
		if (!(D0&(1ULL))) {
            ++err; if (err>cut) return;
        }

		Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;
		++i; ++i_bd; c = 4; 
		if(i_bd < pn) c = seq_nt4_table[(uint8_t)pstr[i_bd]];
		if(c < 4) Peq[c] |= mm; 

		ez->path.a[ez->path.n++] = D0;//diff
        ez->path.a[ez->path.n++] = VP;//diff
        ez->path.a[ez->path.n++] = VN;//diff
        ez->path.a[ez->path.n++] = HP;//diff
        ez->path.a[ez->path.n++] = HN;//diff
	}
	ed_core_64(Peq, VP, VN, X, D0, HN, HP, (uint8_t)tstr[i]);
	if (!(D0&(1ULL))) {
		++err; if (err>cut) return;
	}
	ez->path.a[ez->path.n++] = D0;//diff
	ez->path.a[ez->path.n++] = VP;//diff
	ez->path.a[ez->path.n++] = VN;//diff
	ez->path.a[ez->path.n++] = HP;//diff
	ez->path.a[ez->path.n++] = HN;//diff

	int32_t site = tn - 1 - abs_diag;/**up bound**/
    /**in most cases, ai = (thre<<1)**/
    int32_t ai = pn - tn + abs_diag, uge = INT32_MAX; i = 0;
	if(ez->err > thre) {//diff
		for (i = 0; site < 0 && i < ai; i++, site++) {
			err += ((VP >> i)&(1ULL)); err -= ((VN >> i)&(1ULL));
		}
		if ((err <= thre) && (err <= ez->err)) {
			ez->err = err; ez->pe = site;
		}
		site -= i;
		while (i < ai) {
			err += ((VP >> i)&(1ULL)); err -= ((VN >> i)&(1ULL)); ++i;
			if ((err <= thre) && (err <= ez->err)) {
				ez->err = err; ez->pe = site + i;
			}
			if(i == thre) uge = err;
		}
		
		if((uge <= thre) && (uge == ez->err)) ez->pe = site + thre;
	}
    ///should update ez->path.n for extension
    gen_trace(ez, abs_diag, 1);//diff
}

/**
 pattern is the longer one, while text is the shorter one
 **/
inline int Reserve_Banded_BPM
(char *pattern, int p_length, char *text, int t_length, unsigned short errthold, unsigned int* return_err)
{
	(*return_err) = (unsigned int)-1;

	// bit_extz_t exz; ed_band_cal_semi_128bit(pattern, p_length, text, t_length, errthold, &exz);
	// if(exz.err <= exz.thre) (*return_err) = exz.err;
	// return exz.pe;

	Word Peq[256];

	int band_length = (errthold << 1) + 1;
	int i = 0;
	Word tmp_Peq_1 = (Word)1;

	Peq[(uint8_t)'A'] = (Word)0;
	Peq[(uint8_t)'T'] = (Word)0;
	Peq[(uint8_t)'G'] = (Word)0;
	Peq[(uint8_t)'C'] = (Word)0;


	Word Peq_A;
	Word Peq_T;
	Word Peq_C;
	Word Peq_G;

	///band_length = 2k + 1
	for (i = 0; i<band_length; i++)
	{
		Peq[(uint8_t)pattern[i]] = Peq[(uint8_t)pattern[i]] | tmp_Peq_1;
		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	///Peq[(uint8_t)'T'] = Peq[(uint8_t)'T'] | Peq[(uint8_t)'C'];

	Peq_A = Peq[(uint8_t)'A'];
	Peq_C = Peq[(uint8_t)'C'];
	Peq_T = Peq[(uint8_t)'T'];
	Peq_G = Peq[(uint8_t)'G'];


	memset(Peq, 0, sizeof(Word)* 256);


	Peq[(uint8_t)'A'] = Peq_A;
	Peq[(uint8_t)'C'] = Peq_C;
	Peq[(uint8_t)'T'] = Peq_T;
	Peq[(uint8_t)'G'] = Peq_G;


	

	Word Mask = ((Word)1 << (errthold << 1));

	Word VP = 0;
	Word VN = 0;
	Word X = 0;
	Word D0 = 0;
	Word HN = 0;
	Word HP = 0;


	i = 0;

	int err = 0;

	Word err_mask = (Word)1;


	///band_down = 2k
	///i_bd = 2k
	///int i_bd = i + band_down;
	int i_bd = (errthold << 1);


	int last_high = (errthold << 1);


	/// t_length_1 = SEQ_LENGTH - 1
	int t_length_1 = t_length - 1;
	//while(i<t_length)

	while (i<t_length_1)
	{
		X = Peq[(uint8_t)text[i]] | VN;

		D0 = ((VP + (X&VP)) ^ VP) | X;

		HN = VP&D0;
		HP = VN | ~(VP | D0);

		X = D0 >> 1;
		VN = X&HP;
		VP = HN | ~(X | HP);

		if (!(D0&err_mask))
		{
			++err;

			if ((err - last_high)>(int)errthold)
			{
				return -1;
			}
				
		}


		Peq[(uint8_t)'A'] = Peq[(uint8_t)'A'] >> 1;
		Peq[(uint8_t)'C'] = Peq[(uint8_t)'C'] >> 1;
		Peq[(uint8_t)'G'] = Peq[(uint8_t)'G'] >> 1;
		Peq[(uint8_t)'T'] = Peq[(uint8_t)'T'] >> 1;


		++i;
		++i_bd;
		Peq[(uint8_t)pattern[i_bd]] = Peq[(uint8_t)pattern[i_bd]] | Mask;
	}





	X = Peq[(uint8_t)text[i]] | VN;
	D0 = ((VP + (X&VP)) ^ VP) | X;
	HN = VP&D0;
	HP = VN | ~(VP | D0);
	X = D0 >> 1;
	VN = X&HP;
	VP = HN | ~(X | HP);
	if (!(D0&err_mask))
	{
		++err;
		if ((err - last_high)>errthold)
			return -1;
	}


	////fprintf(stderr, "sucess(2)\n");

	/// last_high = 2k
	/// site = (SEQ_LENGTH + 2k) - 2k -1
	/// site = SEQ_LENGTH - 1
	///int site = p_length - last_high - 1;
	int site = t_length - 1;
	int return_site = -1;
	///in most cases, p_lengthshould be t_length + 2 * errthold
	int available_i = p_length - t_length;
	if ((err <= errthold) && ((unsigned int)err<=*return_err))
	{
		*return_err = err;
		return_site = site;
	}
	i = 0;

	/****************************may have bugs********************************/
	unsigned int ungap_error = (unsigned int)-1;
	/****************************may have bugs********************************/

	while (i < available_i)
	{
		err = err + ((VP >> i)&(Word)1);
		err = err - ((VN >> i)&(Word)1);
		++i;

		if ((err <= (int)errthold) && ((unsigned int)err <= *return_err))
		{
			*return_err = err;
			return_site = site + i;
		}

		/****************************may have bugs********************************/
		if(i == (int)errthold)
		{
			ungap_error = err;
		}
		/****************************may have bugs********************************/
	}

	/****************************may have bugs********************************/
	if((ungap_error<=errthold) && (ungap_error == (*return_err)))
	{
		return_site = site + errthold;
	}
	/****************************may have bugs********************************/

	return return_site;

}



inline int try_cigar(char *pattern, int p_length,
	char *text, int t_length, int end_site, char* path,
	int error, 
	int* return_start_site, 
	int* return_path_length)
{
	int i = 0;
	int tmp_err = 0;
	///start pos of y
	int start_site = end_site - t_length + 1;

	if (start_site >= 0)
	{

		for (i = 0; i < t_length; i++)
		{
			///path[i] = 0;
			///path is saved backwards
			path[t_length - i - 1] = 0;
			if (text[i] != pattern[i + start_site])
			{
				path[t_length - i - 1] = 1;
				tmp_err++;

				if (tmp_err > error)
				{
					return 0;
				}
			}
		}

		if (tmp_err == error)
		{
			(*return_path_length) = t_length;

			(*return_start_site) = start_site;

			return 6;
		}
	}


	return 0;
}



///p_length might be samller than t_length + 2 * errthold
inline int Reserve_Banded_BPM_PATH
(char *pattern, int p_length, char *text, int t_length, unsigned short errthold, 
		unsigned int* return_err, int* return_start_site, int* return_path_length, Word* matrix_bit, char* path, 
		int old_error, int old_end_site)
{
	if (old_error != -1 && old_end_site != -1)
	{
		if (old_error == 0)
		{
			(*return_err) = old_error;
			(*return_start_site) = old_end_site - t_length + 1;
			return old_end_site;
		}

		if (try_cigar(pattern, p_length, text, t_length, old_end_site, path,
		old_error, return_start_site, return_path_length))
		{
			(*return_err) = old_error;
			return old_end_site;
		}
		
	}
	
	(*return_err) = (unsigned int)-1;

	Word Peq[256];

	int band_length = (errthold << 1) + 1;
	int i = 0;
	Word tmp_Peq_1 = (Word)1;

	Peq[(uint8_t)'A'] = (Word)0;
	Peq[(uint8_t)'T'] = (Word)0;
	Peq[(uint8_t)'G'] = (Word)0;
	Peq[(uint8_t)'C'] = (Word)0;


	Word Peq_A;
	Word Peq_T;
	Word Peq_C;
	Word Peq_G;

	///band_length = 2k + 1
	for (i = 0; i<band_length; i++)
	{
		Peq[(uint8_t)pattern[i]] = Peq[(uint8_t)pattern[i]] | tmp_Peq_1;
		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	///Peq[(uint8_t)'T'] = Peq[(uint8_t)'T'] | Peq[(uint8_t)'C'];

	Peq_A = Peq[(uint8_t)'A'];
	Peq_C = Peq[(uint8_t)'C'];
	Peq_T = Peq[(uint8_t)'T'];
	Peq_G = Peq[(uint8_t)'G'];


	memset(Peq, 0, sizeof(Word)* 256);


	Peq[(uint8_t)'A'] = Peq_A;
	Peq[(uint8_t)'C'] = Peq_C;
	Peq[(uint8_t)'T'] = Peq_T;
	Peq[(uint8_t)'G'] = Peq_G;


	

	Word Mask = ((Word)1 << (errthold << 1));

	Word VP = 0;
	Word VN = 0;
	Word X = 0;
	Word D0 = 0;
	Word HN = 0;
	Word HP = 0;


	i = 0;

	Word column_start;

	int err = 0;

	Word err_mask = (Word)1;


	///band_down = 2k
	///i_bd = 2k
	///int i_bd = i + band_down;
	int i_bd = (errthold << 1);


	int last_high = (errthold << 1);


	/// t_length_1 = SEQ_LENGTH - 1
	int t_length_1 = t_length - 1;
	//while(i<t_length)

	while (i<t_length_1)
	{
		///pattern[0]ÔÚPeq[2k], ¶øpattern[2k]ÔÚPeq[0]
		X = Peq[(uint8_t)text[i]] | VN;

		D0 = ((VP + (X&VP)) ^ VP) | X;

		HN = VP&D0;
		HP = VN | ~(VP | D0);

		X = D0 >> 1;
		VN = X&HP;
		VP = HN | ~(X | HP);

		if (!(D0&err_mask))
		{
			++err;

			if ((err - last_high)>(int)errthold)
			{
				return -1;
			}
				
		}


		Peq[(uint8_t)'A'] = Peq[(uint8_t)'A'] >> 1;
		Peq[(uint8_t)'C'] = Peq[(uint8_t)'C'] >> 1;
		Peq[(uint8_t)'G'] = Peq[(uint8_t)'G'] >> 1;
		Peq[(uint8_t)'T'] = Peq[(uint8_t)'T'] >> 1;


		++i;
		++i_bd;
		Peq[(uint8_t)pattern[i_bd]] = Peq[(uint8_t)pattern[i_bd]] | Mask;


		///Peq[(uint8_t)'T'] = Peq[(uint8_t)'T'] | Peq[(uint8_t)'C'];

		column_start = i << 3;
		matrix_bit[column_start] = D0;
		matrix_bit[column_start + 1] = VP;
		matrix_bit[column_start + 2] = VN;
		matrix_bit[column_start + 3] = HP;
		matrix_bit[column_start + 4] = HN;
	}





	X = Peq[(uint8_t)text[i]] | VN;
	D0 = ((VP + (X&VP)) ^ VP) | X;
	HN = VP&D0;
	HP = VN | ~(VP | D0);
	X = D0 >> 1;
	VN = X&HP;
	VP = HN | ~(X | HP);
	if (!(D0&err_mask))
	{
		++err;
		if ((err - last_high)>(int)errthold)
			return -1;
	}


	column_start = (i + 1) << 3;
	matrix_bit[column_start] = D0;
	matrix_bit[column_start + 1] = VP;
	matrix_bit[column_start + 2] = VN;
	matrix_bit[column_start + 3] = HP;
	matrix_bit[column_start + 4] = HN;


	////fprintf(stderr, "sucess(2)\n");

	/// last_high = 2k
	/// site = (SEQ_LENGTH + 2k) - 2k -1
	/// site = SEQ_LENGTH - 1
	///int site = p_length - last_high - 1;
	int site = t_length - 1;
	int return_site = -1;

	/****************************may have bugs********************************/
	unsigned int ungap_error = (unsigned int)-1;
	/****************************may have bugs********************************/

	///in most cases, p_length should be t_length + 2 * errthold
	int available_i = p_length - t_length;
	if ((err <= (int)errthold) && ((unsigned int)err<=*return_err))
	{
		*return_err = err;
		return_site = site;
	}
	i = 0;

	while (i < available_i)
	{
		err = err + ((VP >> i)&(Word)1);
		err = err - ((VN >> i)&(Word)1);
		++i;

		if ((err <= (int)errthold) && ((unsigned int)err <= *return_err))
		{
			*return_err = err;
			return_site = site + i;
		}

		/****************************may have bugs********************************/
		if(i == (int)errthold)
		{
			ungap_error = err;
		}
		/****************************may have bugs********************************/
	}

	


	if ((*return_err) == (unsigned int)-1)
	{
		return return_site;
	}


	/****************************may have bugs********************************/
	if((ungap_error<=errthold) && (ungap_error == (*return_err)))
	{
		return_site = site + errthold;
	}
	/****************************may have bugs********************************/
	

	///need to correct p_length here, since p_length might be smaller than t_length + 2* err_threashlod
	p_length = t_length + 2 * errthold;
	///end_site is always correct
	int end_site = return_site;
	int start_site = end_site;
	int back_track_site = band_length - (p_length - end_site);

	Word v_value, h_value, delta_value, min_value, current_value;
	///Word direction; ///0 is match, 1 is mismatch, 2 is up, 3 is left
	Word direction = 0; ///0 is match, 1 is mismatch, 2 is up, 3 is left
	i = t_length;
	int path_length = 0;
	current_value = *return_err;


	int low_bound = band_length - 1;


	while (i>0)
	{
		if (current_value == 0)
		{
			break;
		}

		column_start = i << 3;

		delta_value = current_value -
			((~(matrix_bit[column_start] >> back_track_site))&err_mask);


		if (back_track_site == 0)
		{
			///HP
			h_value = current_value - ((matrix_bit[column_start + 3] >> back_track_site)&err_mask);
			//HN
			h_value = h_value + ((matrix_bit[column_start + 4] >> back_track_site)&err_mask);


			min_value = delta_value;
			direction = 0;

			if (h_value < min_value)
			{
				min_value = h_value;
				direction = 3;
			}

		}
		else if (back_track_site == low_bound)
		{
			v_value = current_value - ((matrix_bit[column_start + 1] >> (back_track_site - 1))&err_mask);
			v_value = v_value + ((matrix_bit[column_start + 2] >> (back_track_site - 1))&err_mask);

			min_value = delta_value;
			direction = 0;

			if (v_value < min_value)
			{
				min_value = v_value;
				direction = 2;
			}

		}
		else
		{

			h_value = current_value - ((matrix_bit[column_start + 3] >> back_track_site)&(Word)1);


			h_value = h_value + ((matrix_bit[column_start + 4] >> back_track_site)&(Word)1);


			v_value = current_value - ((matrix_bit[column_start + 1] >> (back_track_site - 1))&err_mask);
			v_value = v_value + ((matrix_bit[column_start + 2] >> (back_track_site - 1))&err_mask);


			min_value = delta_value;
			direction = 0;

			if (v_value < min_value)
			{
				min_value = v_value;
				direction = 2;
			}


			if (h_value < min_value)
			{
				min_value = h_value;
				direction = 3;
			}

		}


		if (direction == 0)
		{

			if (delta_value != current_value)
			{
				direction = 1;
			}

			i--;

			start_site--;

		}
		if (direction == 2)///ru guo xiang shang yi dong, bing bu huan lie
		{
			back_track_site--;
			start_site--;
		}
		else if (direction == 3)///ru guo xiang zuo yi dong
		{
			i--;
			back_track_site++;
		}


		path[path_length++] = direction;


		current_value = min_value;

	}


	if (i > 0)
	{
		memset(path + path_length, 0, i);
		start_site = start_site - i;
		direction = 0;
		path_length = path_length + i;
	}

	if (direction != 3)
	{
		start_site++;
	}
	(*return_start_site) = start_site;
	(*return_path_length) = path_length;

	return return_site;
}

////four patterns have the same p_length
inline int Reserve_Banded_BPM_4_SSE_only(char *pattern1, char *pattern2, char *pattern3, char *pattern4, int p_length, char *text, int t_length,
	int* return_sites, unsigned int* return_sites_error, unsigned short errthold, __m128i* Peq_SSE)

{
	memset(return_sites, -1, sizeof(int)* 4);
	memset(return_sites_error, -1, sizeof(unsigned int)* 4);

	Word_32 Peq[256][4];
	int band_length = (errthold << 1) + 1;


	int i;

	Word_32 tmp_Peq_1 = 1;


	memset(Peq[(uint8_t)'A'], 0, sizeof(Word_32)* 4);
	memset(Peq[(uint8_t)'C'], 0, sizeof(Word_32)* 4);
	memset(Peq[(uint8_t)'G'], 0, sizeof(Word_32)* 4);
	memset(Peq[(uint8_t)'T'], 0, sizeof(Word_32)* 4);

	for (i = 0; i<band_length; i++)
	{
		Peq[(uint8_t)pattern1[i]][0] = Peq[(uint8_t)pattern1[i]][0] | tmp_Peq_1;
		Peq[(uint8_t)pattern2[i]][1] = Peq[(uint8_t)pattern2[i]][1] | tmp_Peq_1;
		Peq[(uint8_t)pattern3[i]][2] = Peq[(uint8_t)pattern3[i]][2] | tmp_Peq_1;
		Peq[(uint8_t)pattern4[i]][3] = Peq[(uint8_t)pattern4[i]][3] | tmp_Peq_1;

		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	Peq_SSE[(uint8_t)'A'] = _mm_set_epi32(Peq[(uint8_t)'A'][3], Peq[(uint8_t)'A'][2], Peq[(uint8_t)'A'][1], Peq[(uint8_t)'A'][0]);
	Peq_SSE[(uint8_t)'C'] = _mm_set_epi32(Peq[(uint8_t)'C'][3], Peq[(uint8_t)'C'][2], Peq[(uint8_t)'C'][1], Peq[(uint8_t)'C'][0]);
	Peq_SSE[(uint8_t)'G'] = _mm_set_epi32(Peq[(uint8_t)'G'][3], Peq[(uint8_t)'G'][2], Peq[(uint8_t)'G'][1], Peq[(uint8_t)'G'][0]);
	Peq_SSE[(uint8_t)'T'] = _mm_set_epi32(Peq[(uint8_t)'T'][3], Peq[(uint8_t)'T'][2], Peq[(uint8_t)'T'][1], Peq[(uint8_t)'T'][0]);
	///Peq_SSE[(uint8_t)'T'] = _mm_or_si128(Peq_SSE[(uint8_t)'T'], Peq_SSE[(uint8_t)'C']);

	Word_32 Mask = ((Word_32)1 << (errthold << 1));

	__m128i Mask1 = _mm_set_epi32(0, 0, 0, Mask);
	__m128i Mask2 = _mm_set_epi32(0, 0, Mask, 0);
	__m128i Mask3 = _mm_set_epi32(0, Mask, 0, 0);
	__m128i Mask4 = _mm_set_epi32(Mask, 0, 0, 0);


	__m128i VP = _mm_setzero_si128();
	__m128i VN = _mm_setzero_si128();
	__m128i X = _mm_setzero_si128();
	__m128i D0 = _mm_setzero_si128();
	__m128i HN = _mm_setzero_si128();
	__m128i HP = _mm_setzero_si128();
	__m128i tmp_process;
	__m128i tmp_process1;



	__m128i Err_4 = _mm_setzero_si128();
	__m128i err_mask = _mm_set_epi32(1, 1, 1, 1);
	///for_not li quan shi 1
	__m128i for_not = _mm_set1_epi32(-1);
	__m128i err_arry;
	__m128i cmp_result;


	int i_bd = (errthold << 1);
	int last_high = (errthold << 1);
	int t_length_1 = t_length - 1;
	int err1;
	int err2;
	int err3;
	int err4;

	__m128i pre_end = _mm_set_epi32
		(last_high + errthold, last_high + errthold, last_high + errthold, last_high + errthold);
	i = 0;

	while (i<t_length_1)
	{
		///X = Peq[text[i]] | VN;
		X = _mm_or_si128(Peq_SSE[(uint8_t)text[i]], VN);



		/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/
		///X&VP
		tmp_process1 = _mm_and_si128(X, VP);
		///(VP + (X&VP))
		tmp_process = _mm_add_epi32(tmp_process1, VP);
		///((VP + (X&VP)) ^ VP)
		tmp_process = _mm_xor_si128(tmp_process, VP);
		///((VP + (X&VP)) ^ VP) | X
		D0 = _mm_or_si128(tmp_process, X);
		/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/

		///HN = VP&D0;
		HN = _mm_and_si128(D0, VP);

		///HP = VN | ~(VP | D0);
		tmp_process = _mm_or_si128(D0, VP);
		tmp_process = _mm_andnot_si128(tmp_process, for_not);
		HP = _mm_or_si128(tmp_process, VN);

		///X = D0 >> 1;
		X = _mm_srli_epi32(D0, 1);
		///VN = X&HP;
		VN = _mm_and_si128(X, HP);
		///VP = HN | ~(X | HP);
		tmp_process = _mm_or_si128(X, HP);
		tmp_process = _mm_andnot_si128(tmp_process, for_not);
		VP = _mm_or_si128(HN, tmp_process);

		///D0&err_mask
		err_arry = _mm_and_si128(D0, err_mask);
		Err_4 = _mm_add_epi32(Err_4, err_mask);
		Err_4 = _mm_sub_epi32(Err_4, err_arry);

		/**************** */
		///shi ji shang zhe ge zhi hen xiao d
		cmp_result = _mm_cmpgt_epi32(Err_4, pre_end);

		///jian zhi
		if (_mm_extract_epi32(cmp_result, 0) && _mm_extract_epi32(cmp_result, 1)
			&& _mm_extract_epi32(cmp_result, 2) && _mm_extract_epi32(cmp_result, 3))
			return 1;
		/**************** */


		Peq_SSE[(uint8_t)'A'] = _mm_srli_epi32(Peq_SSE[(uint8_t)'A'], 1);
		Peq_SSE[(uint8_t)'T'] = _mm_srli_epi32(Peq_SSE[(uint8_t)'T'], 1);
		Peq_SSE[(uint8_t)'G'] = _mm_srli_epi32(Peq_SSE[(uint8_t)'G'], 1);
		Peq_SSE[(uint8_t)'C'] = _mm_srli_epi32(Peq_SSE[(uint8_t)'C'], 1);

		++i;
		++i_bd;

		Peq_SSE[(uint8_t)pattern1[i_bd]] = _mm_or_si128(Mask1, Peq_SSE[(uint8_t)pattern1[i_bd]]);
		Peq_SSE[(uint8_t)pattern2[i_bd]] = _mm_or_si128(Mask2, Peq_SSE[(uint8_t)pattern2[i_bd]]);
		Peq_SSE[(uint8_t)pattern3[i_bd]] = _mm_or_si128(Mask3, Peq_SSE[(uint8_t)pattern3[i_bd]]);
		Peq_SSE[(uint8_t)pattern4[i_bd]] = _mm_or_si128(Mask4, Peq_SSE[(uint8_t)pattern4[i_bd]]);
		///Peq_SSE[(uint8_t)'T'] = _mm_or_si128(Peq_SSE[(uint8_t)'T'], Peq_SSE[(uint8_t)'C']);
	}



		///X = Peq[text[i]] | VN;
	X = _mm_or_si128(Peq_SSE[(uint8_t)text[i]], VN);

	/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/
	///X&VP
	tmp_process1 = _mm_and_si128(X, VP);
	///(VP + (X&VP))
	tmp_process = _mm_add_epi32(tmp_process1, VP);
	///((VP + (X&VP)) ^ VP)
	tmp_process = _mm_xor_si128(tmp_process, VP);
	///((VP + (X&VP)) ^ VP) | X
	D0 = _mm_or_si128(tmp_process, X);
	/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/

	///HN = VP&D0;
	HN = _mm_and_si128(D0, VP);

	///HP = VN | ~(VP | D0);
	tmp_process = _mm_or_si128(D0, VP);
	tmp_process = _mm_andnot_si128(tmp_process, for_not);
	HP = _mm_or_si128(tmp_process, VN);


	///X = D0 >> 1;
	X = _mm_srli_epi32(D0, 1);
	///VN = X&HP;
	VN = _mm_and_si128(X, HP);
	///VP = HN | ~(X | HP);
	tmp_process = _mm_or_si128(X, HP);
	tmp_process = _mm_andnot_si128(tmp_process, for_not);
	VP = _mm_or_si128(HN, tmp_process);

	///D0&err_mask
	err_arry = _mm_and_si128(D0, err_mask);
	Err_4 = _mm_add_epi32(Err_4, err_mask);
	Err_4 = _mm_sub_epi32(Err_4, err_arry);

	///shi ji shang zhe ge zhi hen xiao d
	cmp_result = _mm_cmpgt_epi32(Err_4, pre_end);

	///jian zhi
	if (_mm_extract_epi32(cmp_result, 0) && _mm_extract_epi32(cmp_result, 1)
		&& _mm_extract_epi32(cmp_result, 2) && _mm_extract_epi32(cmp_result, 3))
		return 1;

	int site = t_length - 1;
	err1 = _mm_extract_epi32(Err_4, 0);
	err2 = _mm_extract_epi32(Err_4, 1);
	err3 = _mm_extract_epi32(Err_4, 2);
	err4 = _mm_extract_epi32(Err_4, 3);


	if ((err1 <= (int)errthold) && ((unsigned int)err1 <= return_sites_error[0]))
	{
		return_sites[0] = site;
		return_sites_error[0] = err1;
	}
	if ((err2 <= (int)errthold) && ((unsigned int)err2 <= return_sites_error[1]))
	{
		return_sites[1] = site;
		return_sites_error[1] = err2;
	}
	if ((err3 <= (int)errthold) && ((unsigned int)err3 <= return_sites_error[2]))
	{
		return_sites[2] = site;
		return_sites_error[2] = err3;
	}
	if ((err4 <= (int)errthold) && ((unsigned int)err4 <= return_sites_error[3]))
	{
		return_sites[3] = site;
		return_sites_error[3] = err4;
	}


	i = 0;

	/****************************may have bugs********************************/
	unsigned int ungap_error1 = (unsigned int)-1;
	unsigned int ungap_error2 = (unsigned int)-1;
	unsigned int ungap_error3 = (unsigned int)-1;
	unsigned int ungap_error4 = (unsigned int)-1;
	/****************************may have bugs********************************/


	///in most cases, p_length should be t_length + 2 * errthold
	int available_i = p_length - t_length;

	while (i < available_i)
	{
		///err = err + ((VP >> i)&(Word_32)1);
		tmp_process = _mm_srli_epi32(VP, i);
		tmp_process = _mm_and_si128(tmp_process, err_mask);
		Err_4 = _mm_add_epi32(Err_4, tmp_process);

		///err = err - ((VN >> i)&(Word_32)1);
		tmp_process1 = _mm_srli_epi32(VN, i);
		tmp_process1 = _mm_and_si128(tmp_process1, err_mask);
		Err_4 = _mm_sub_epi32(Err_4, tmp_process1);
		++i;

		err1 = _mm_extract_epi32(Err_4, 0);
		err2 = _mm_extract_epi32(Err_4, 1);
		err3 = _mm_extract_epi32(Err_4, 2);
		err4 = _mm_extract_epi32(Err_4, 3);


		if ((err1 <= (int)errthold) && ((unsigned int)err1 <= return_sites_error[0]))
		{
			return_sites[0] = site + i;
			return_sites_error[0] = err1;
		}
		if ((err2 <= (int)errthold) && ((unsigned int)err2 <= return_sites_error[1]))
		{
			return_sites[1] = site + i;
			return_sites_error[1] = err2;
		}
		if ((err3 <= (int)errthold) && ((unsigned int)err3 <= return_sites_error[2]))
		{
			return_sites[2] = site + i;
			return_sites_error[2] = err3;
		}
		if ((err4 <= (int)errthold) && ((unsigned int)err4 <= return_sites_error[3]))
		{
			return_sites[3] = site + i;
			return_sites_error[3] = err4;
		}

		/****************************may have bugs********************************/
		if(i == (int)errthold)
		{
			ungap_error1 = err1;
			ungap_error2 = err2;
			ungap_error3 = err3;
			ungap_error4 = err4;
		}
		/****************************may have bugs********************************/
	}

	/****************************may have bugs********************************/
	if((ungap_error1<=errthold) && (ungap_error1 == return_sites_error[0]))
	{
		return_sites[0] = site + errthold;
	}

	if((ungap_error2<=errthold) && (ungap_error2 == return_sites_error[1]))
	{
		return_sites[1] = site + errthold;
	}

	if((ungap_error3<=errthold) && (ungap_error3 == return_sites_error[2]))
	{
		return_sites[2] = site + errthold;
	}

	if((ungap_error4<=errthold) && (ungap_error4 == return_sites_error[3]))
	{
		return_sites[3] = site + errthold;
	}
	/****************************may have bugs********************************/

	return 1;
}


// void move_trace_gap(uint16_t *trace, int32_t trace_n, int32_t trace_i, 
// char *pstr, int32_t pi, char *tstr, int32_t ti, int32_t *err)
// {
// 	uint16_t c = trace[trace_i]>>14, l = (trace[trace_i]<<2)>>2; 
// 	if(c != 3 && c != 2) return;
// 	trace_i--;
// 	if(c == 3) pi--;
// 	else if(c == 2) ti--;
	

// 	if()

// }

// void adjust_trace(uint16_t *trace, int32_t *trace_n, int32_t *p_beg, int32_t *p_end, int32_t *err, char *pstr, char *tstr)
// {
// 	if((*err) == 0) return;
// 	int32_t i, pi, ti; uint16_t c, l;
// 	for (i = 0; i < (*trace_n) && (trace[i]>>14) == 1; i++) {
// 		trace[i] <<= 2; trace[i] >>= 2; trace[i] += (((uint16_t)3)<<14); 
// 		l = (trace[i]<<2)>>2; (*p_beg) += l;
// 	}
// 	for (i = (*trace_n) - 1; i >= 0 && (trace[i]>>14) == 1; i--) {
// 		trace[i] <<= 2; trace[i] >>= 2; trace[i] += (((uint16_t)3)<<14); 
// 		l = (trace[i]<<2)>>2; (*p_end) -= l;
// 	}

// 	i = 0; pi = (*p_beg); ti = 0;
// 	for (i = 0; i < (*trace_n); i++) {
// 		c = trace[i]>>14; l = (trace[i]<<2)>>2;
// 		if(c == 0 || c == 1) {
//             pi += l; ti += l;
//         } else if(c == 2) {
// 			// move_trace_gap(trace, *trace_n, i, pstr, pi, tstr, ti, err);
//             pi += l;
//         } else if(c == 3) {
//             // move_trace_gap(trace, *trace_n, i, pstr, pi, tstr, ti, err);
//             ti += l;
//         }
// 	}
// }

inline int32_t ungap_trace(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t know_err, int32_t know_end, 
int32_t *r_err, int32_t *r_beg, asg16_v *cigar, int32_t *cigar_l)
{
	int32_t cn = cigar->n, pk, tk, e, l;
	(*r_err) = (*r_beg) = INT32_MAX; (*cigar_l) = 0;
	if(know_err < 0 || know_end < 0) return -1;
	if(know_err == 0) {
		push_trace(cigar, EAC_M, tn);
		(*r_err) = know_err; (*r_beg) = know_end + 1 - tn; (*cigar_l) = cigar->n - cn;
		return know_end;
	}

	pk = know_end+1-tn; tk = 0; e = 0;
	for (l = 0; tk < tn; tk++, pk++) {
		if(pstr[pk]!=tstr[tk]) {
			e++; if(e > know_err) break;
			if(tk > l) push_trace(cigar, EAC_M, tk-l);
			push_trace(cigar, MIS_M, 1); l = tk + 1;
		}
	}
	if(tk == tn) {
		if(tk > l) push_trace(cigar, EAC_M, tk-l);
		(*r_err) = know_err; (*r_beg) = know_end + 1 - tn; (*cigar_l) = cigar->n - cn;
		return know_end;
	}

	cigar->n = cn;
	return -1;
}


// ///p_length might be samller than t_length + 2 * errthold
inline int32_t ed_band_cal_semi_trace(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, 
		int32_t know_err, int32_t know_end, int32_t *r_err, int32_t *r_beg, Word *buf, asg16_v *cigar, int32_t *cigar_l) {
	int32_t cn = cigar->n; (*r_err) = (*r_beg) = INT32_MAX; (*cigar_l) = 0; 
	Word Peq[5] = {0}, mm = (Word)1, VP = 0, VN = 0, X = 0, D0 = 0, HN = 0, HP = 0, i_col, i_col_dux; 
    int32_t bd = (thre<<1)+1, i, err = 0, i_bd = (thre<<1), last_high = (thre<<1), tn0 = tn - 1;
    int32_t cut = thre+last_high; ///kv_resize(uint16_t, *cigar, cigar->n+(uint32_t)know_err+2);//pre-alloc

    if(ungap_trace(pstr, pn, tstr, tn, know_err, know_end, r_err, r_beg, cigar, cigar_l) >= 0) {
		return know_end;
	}

    for (i = 0; i < bd; i++) {
        Peq[seq_nt4_table[(uint8_t)pstr[i]]] |= mm; mm <<= 1;
    }
	Peq[4] = 0;
	///should make Peq[4] = 0 if N is always an error
    i = i_col = 0; mm = ((Word)1 << (thre<<1));///for the incoming char/last char

    while (i < tn0) {
        X = Peq[seq_nt4_table[(uint8_t)tstr[i]]] | VN;

        D0 = ((VP + (X&VP)) ^ VP) | X;

        HN = VP&D0;
        HP = VN | ~(VP | D0);

        X = D0 >> 1;
        VN = X&HP;
        VP = HN | ~(X | HP);

        if (!(D0&(1ULL))) {
            ++err;
			if (err>cut) return -1;    
        }

		Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1; ///Peq[4] >>= 1;

        ++i; ++i_bd; 
        Peq[seq_nt4_table[(uint8_t)pstr[i_bd]]] |= mm; Peq[4] = 0;

		buf[i_col++] = D0; buf[i_col++] = VP; buf[i_col++] = VN; buf[i_col++] = HP; buf[i_col++] = HN;
    }

    X = Peq[seq_nt4_table[(uint8_t)tstr[i]]] | VN;
    D0 = ((VP + (X&VP)) ^ VP) | X;
    HN = VP&D0;
    HP = VN | ~(VP | D0);
    X = D0 >> 1;
    VN = X&HP;
    VP = HN | ~(X | HP);
    if (!(D0&(1ULL))) {
		++err;
		if (err>cut) return -1;    
	}

    buf[i_col++] = D0; buf[i_col++] = VP; buf[i_col++] = VN; buf[i_col++] = HP; buf[i_col++] = HN;

	i_col_dux = i_col/tn;

    int32_t site = tn - 1, end = -1;///up bound
	///in most cases, ai = (thre<<1)
    int32_t ai = pn - tn, uge = INT32_MAX;
    if ((err <= thre) && (err<=(*r_err))) {
        *r_err = err; end = site;
    }
    i = 0;

    while (i < ai) {
        err += ((VP >> i)&(1ULL)); err -= ((VN >> i)&(1ULL)); ++i;
        if ((err <= thre) && (err <= (*r_err))) {
            *r_err = err; end = site + i;
        }
        if(i == thre) uge = err;
    }
	if((uge<=thre) && (uge == (*r_err))) end = site + thre;
    if ((*r_err) > thre) return end;
    

    ///need to correct pn here, since pn might be smaller than tn + 2* thre
    pn = tn + (thre<<1);
    int32_t beg = end, back_track_site = bd - (pn - end);

    Word v_value, h_value, delta_value, min_value, current_value;
    ///Word direction; ///0 is match, 1 is mismatch, 2 is up, 3 is left
    Word direction = 0, *ba, pd, pdl; ///0 is match, 1 is mismatch, 2 is up, 3 is left
    i = tn; pd = (Word)-1; pdl = 0; 
    current_value = *r_err;
    int32_t low_bound = bd - 1;

    while (i > 0) {
        if (current_value == 0) break;
		ba = buf + ((i*i_col_dux) - i_col_dux);
        delta_value = current_value - ((~(ba[0]>>back_track_site))&(1ULL));

        if (back_track_site == 0) {
            ///HP
            h_value = current_value - ((ba[3] >> back_track_site)&(1ULL));
            //HN
            h_value = h_value + ((ba[4] >> back_track_site)&1ULL);
            min_value = delta_value; direction = 0;
            if (h_value < min_value) {
                min_value = h_value;
                direction = 3;
            }
        } else if (back_track_site == low_bound) {
            v_value = current_value - ((ba[1]>>(back_track_site-1))&(1ULL));
            v_value = v_value + ((ba[2]>>(back_track_site-1))&(1ULL));

            min_value = delta_value; direction = 0;
            if (v_value < min_value) {
                min_value = v_value;
                direction = 2;
            }
        }
        else {
            h_value = current_value-((ba[3]>>back_track_site)&(1ULL));
            h_value = h_value+((ba[4]>>back_track_site)&(1ULL));

            v_value = current_value - ((ba[1]>>(back_track_site-1))&(1ULL));
            v_value = v_value + ((ba[2]>>(back_track_site-1))&(1ULL));

            min_value = delta_value; direction = 0;
            if (v_value < min_value) {
                min_value = v_value;
                direction = 2;
            }

            if (h_value < min_value) {
                min_value = h_value;
                direction = 3;
            }
        }


        if (direction == 0) {
            if (delta_value != current_value) {
                direction = 1;
            }
            i--; beg--;
        }
        if (direction == 2) {///ru guo xiang shang yi dong, bing bu huan lie
            back_track_site--; beg--;
        }
        else if (direction == 3) {///ru guo xiang zuo yi dong
            i--;
            back_track_site++;
        }
		
		if(direction != pd) {
			if(pdl > 0) push_trace(cigar, pd, pdl);
			pd = direction; pdl = 1;
		} else {
			pdl++;
		}
        // path[path_length++] = direction;
        current_value = min_value;
    }


    if (i > 0) {
		direction = 0; beg -= i;
		if(direction != pd) {
			if(pdl > 0) push_trace(cigar, pd, pdl);
			pd = direction; pdl = i;
		} else {
			pdl += i;
		}
    }

	if(pdl > 0) push_trace(cigar, pd, pdl);
    if (direction != 3) beg++;

	uint16_t *trac = cigar->a + cn, tt; int32_t trac_n = cigar->n - cn; ai = trac_n>>1;
	for (i = 0; i < ai; i++) {
		tt = trac[i]; trac[i] = trac[trac_n-i-1]; trac[trac_n-i-1] = tt;
	}
	(*cigar_l) = cigar->n - cn;

    (*r_beg) = beg; (*cigar_l) = cigar->n - cn;
    return end;
}


#endif
