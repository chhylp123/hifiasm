#ifndef __LEVENSHTEIN__
#define __LEVENSHTEIN__
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

///p_length might be samller than t_length + 2 * errthold
/// pattern is longer than text
inline int32_t ed_band_cal_semi(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, int32_t *re_err)
{
    (*re_err) = INT32_MAX;
	Word Peq[5] = {0}, mm = (Word)1, VP = 0, VN = 0, X = 0, D0 = 0, HN = 0, HP = 0; 
	int32_t bd = (thre<<1)+1, i, err = 0, i_bd = (thre<<1), last_high = (thre<<1), tn0 = tn - 1;
	int32_t cut = thre+last_high; 

    for (i = 0; i < bd; i++) {
        Peq[seq_nt4_table[(uint8_t)pstr[i]]] |= mm; mm <<= 1;
    }
	///should make Peq[4] = 0 if N is always an error
	Peq[4] = 0;
	i = 0; mm = ((Word)1 << (thre<<1));///for the incoming char/last char

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

    int32_t site = tn - 1, end = -1;///up bound
    ///in most cases, ai = (thre<<1)
    int32_t ai = pn - tn, uge = INT32_MAX;
    if ((err <= thre) && (err<=(*re_err))) {
        *re_err = err; end = site;
    }
    i = 0;

    while (i < ai) {
        err += ((VP >> i)&(1ULL)); err -= ((VN >> i)&(1ULL)); ++i;
        if ((err <= thre) && (err <= (*re_err))) {
            *re_err = err; end = site + i;
        }
        if(i == thre) uge = err;
    }

    if((uge<=thre) && (uge == (*re_err))) end = site + thre;
    return end;
}

inline void print_bit(Word z, int64_t w, const char *cmd)
{
	int64_t k;//, w = (sizeof(Word)<<3);
	fprintf(stderr, "%s\t", cmd);
	for (k = 0; k < w; k++) fprintf(stderr, "%llu", (z>>k)&(1ULL));
	fprintf(stderr, "\n");
}

inline void print_bits(Word *az, int64_t w, const char *cmd)
{
	int64_t k, m, s = (sizeof(*az)<<3), sw = (w/s) + (!!(w%s)), ks;
	fprintf(stderr, "%s\t", cmd);
	for (m = k = 0; m < sw && k < w; m++) {
		for (ks = 0; ks < s && k < w; ks++, k++) fprintf(stderr, "%llu", (az[m]>>ks)&(1ULL));
	}
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

typedef uint64_t w_sig;
typedef struct {w_sig a[2];} w128_t;
#define bitw (6)
#define bitwbit (64)
#define bitz (63)
// typedef uint32_t w_sig;
// typedef struct {w_sig a[2];} w128_t;
// #define bitw (5)
// #define bitwbit (32)
// #define bitz (31)
typedef struct {size_t n, m; w_sig *a;} w64_trace_t;
typedef struct {
	int32_t done_cigar, cigar_n, done_path, path_n;
	int32_t ps, pe, pl, ts, te, tl;
	int32_t thre, err, nword, mword;
	w128_t Peq[5], mm, VP, VN, X, D0, HN, HP;
	// asg16_v cigar; w64_trace_t path; 
} bit_extz_t;

#define w128_bit(x, b) ((x).a[((b)>>bitw)]|=(((w_sig)1)<<((b)&bitz)))

#define w128_get_bit(x, b) (((x).a[((b)>>bitw)]>>((b)&bitz))&((w_sig)1))

#define w128_clear(x) ((x).a[0]=(x).a[1]=0)

#define w128_self_not(x) ((x).a[0]=~(x).a[0], \
									(x).a[1]=~(x).a[1])

#define w128_self_or(x, y) ((x).a[0]|=(y).a[0], \
									(x).a[1]|=(y).a[1])

#define w128_or(r, x, y) ((r).a[0] = (x).a[0]|(y).a[0], \
									(r).a[1] = (x).a[1]|(y).a[1])

#define w128_and(r, x, y) ((r).a[0] = (x).a[0]&(y).a[0], \
									(r).a[1] = (x).a[1]&(y).a[1])

#define w128_self_xor(x, y) ((x).a[0]^=(y).a[0], \
									(x).a[1]^=(y).a[1])

// #define w128_self_lsft_l(x, l) ((x).a[1] = ((x).a[1]<<(l))|((x).a[0]>>(bitwbit-(l))), (x).a[0] <<= (l))

#define w128_self_lsft_1(x) ((x).a[1] = ((x).a[1]<<1)|((x).a[0]>>bitz), \
																		(x).a[0] <<= 1)

#define w128_self_rsft_1(x) ((x).a[0] = ((x).a[0]>>1)|((x).a[1]<<bitz), \
																		(x).a[1] >>= 1)

#define w128_self_add(x, y) ((x).a[0]+=(y).a[0], \
								(x).a[1]+=(y).a[1]+((x).a[0]<(y).a[0]))

#define w128_set_bit_lsub(x, l) do {	\
		(x).a[0] = (w_sig)-1, (x).a[1] = 0;	\
		if((l) <= bitwbit) (x).a[0] = (((w_sig)1)<<(l))-1; \
		else (x).a[1] = (((w_sig)1)<<((l)-bitwbit))-1;\
	} while (0)												\

#define ed_core_w128(Peq, VP, VN, X, D0, HN, HP) {	\
		/**X = Peq[seq_nt4_table[(uint8_t)tstr[i]]] | VN;**/\
		c = seq_nt4_table[(uint8_t)tstr[i]]; w128_or(X, Peq[c], VN);\
		/**D0 = ((VP + (X&VP)) ^ VP) | X;**/\
		w128_and(D0, X, VP);\
		w128_self_add(D0, VP);\
		w128_self_xor(D0, VP);\
		w128_self_or(D0, X);\
        /**HN = VP&D0;**/\
		w128_and(HN, VP, D0);\
        /**HP = VN | ~(VP | D0);**/\
		w128_or(HP, VP, D0);\
		w128_self_not(HP);\
		w128_self_or(HP, VN);\
        /**X = D0 >> 1;**/\
		X = D0; w128_self_rsft_1(X);\
        /**VN = X&HP;**/\
		w128_and(VN, X, HP);\
        /**VP = HN | ~(X | HP);**/\
		w128_or(VP, X, HP);\
		w128_self_not(VP);\
		w128_self_or(VP, HN);\
	} 

#define ed_core_w128_reshift(Peq) {\
		/** Peq[0] >>= 1; Peq[1] >>= 1; Peq[2] >>= 1; Peq[3] >>= 1;**/\
		w128_self_rsft_1(Peq[0]); w128_self_rsft_1(Peq[1]);\
		w128_self_rsft_1(Peq[2]); w128_self_rsft_1(Peq[3]);\
}

#define init_base_ed(ez, thre, pn, tn) {\
	(ez).thre = (thre), (ez).err = INT32_MAX, (ez).pl = pn, (ez).tl = tn;\
	(ez).done_cigar = (ez).done_path = (ez).cigar_n = (ez).path_n = 0;\
}

#define tst_band_err(p_off, t_off, p_end, dif, ez, k) {\
	if((p_off) + (((ez).thre)<<1) >= (p_end)) { \
		for ((k) = 0; (p_off) < (p_end); (p_off)++) {\
			(dif) += w128_get_bit((ez).VP, (k));\
			(dif) -= w128_get_bit((ez).VN, (k));\
			(k)++;}\
		if((dif) <= (ez).thre && (dif) < (ez).err) {\
			(ez).err = dif; (ez).pe = p_off; (ez).te = t_off;}\
	}\
}


inline void ed_band_cal_global_128bit(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{
	init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = 0;
	if((pn > tn + thre) || (tn > pn + thre)) return;
	// if((pn < thre + 1) || (tn < thre + 1)) return;
	int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd = thre+1, i_bd = thre; uint8_t c;
	w128_clear(ez->Peq[0]); w128_clear(ez->Peq[1]); w128_clear(ez->Peq[2]); w128_clear(ez->Peq[3]); w128_clear(ez->Peq[4]);

	w128_clear(ez->mm); w128_bit(ez->mm, thre); ///mm = (((Word)1)<<thre)
    for (i = 0; i < bd && i < pn; i++) {
		w128_self_or(ez->Peq[seq_nt4_table[(uint8_t)pstr[i]]], ez->mm); w128_self_lsft_1(ez->mm);
        // Peq[seq_nt4_table[(uint8_t)pstr[i]]] |= mm; mm <<= 1;
    }
	w128_clear(ez->Peq[4]);
	err = thre;
	w128_set_bit_lsub(ez->VN, thre); ///VN = (((Word)1)<<(thre))-1; 
	w128_set_bit_lsub(ez->VP, (thre<<1)+1); ///VP = (((Word)1)<<((thre<<1)+1))-1;
	w128_self_xor(ez->VP, ez->VN); ///VP ^= VN;

	// print_bits(Peq[0].a, (thre<<1)+1, "-Peq[A]");
	// print_bits(Peq[1].a, (thre<<1)+1, "-Peq[C]");
	// print_bits(Peq[2].a, (thre<<1)+1, "-Peq[G]");
	// print_bits(Peq[3].a, (thre<<1)+1, "-Peq[T]");
	// print_bits(VN.a, (thre<<1)+1, "-VN");
	// print_bits(VP.a, (thre<<1)+1, "-VP");

	///should make Peq[4] = 0 if N is always an error
	i = 0; 
	///for the incoming char/last char
	w128_clear(ez->mm); w128_bit(ez->mm, (thre<<1)); ///mm = ((Word)1 << (thre<<1));
	//VP + ((Peq|VN)&VP)
    while (i < tn0) {
		ed_core_w128(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP);
        if (!(ez->D0.a[0]&(1ULL))) {
            ++err; if (err>cut) return;
        }
		ed_core_w128_reshift(ez->Peq);
        ++i; ++i_bd;
        if(i_bd < pn) {
			c = seq_nt4_table[(uint8_t)pstr[i_bd]]; 
			///if(c < 4) Peq[c] |= mm;
			if(c < 4) w128_self_or(ez->Peq[c], ez->mm);
		}
    }
	ed_core_w128(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP);
	if (!(ez->D0.a[0]&(1ULL))) {
		++err; if (err>cut) return;
	}

    int32_t site = tn - 1 - thre;///up bound
	for (cut = pn - 1; site < cut; site++) {
		// err += ((VP >> i)&(1ULL)); 
		err += ez->VP.a[0]&(1ULL); w128_self_rsft_1(ez->VP); 
		// err -= ((VN >> i)&(1ULL));
		err -= ez->VN.a[0]&(1ULL); w128_self_rsft_1(ez->VN); 
	}

	if (site == cut && err <= thre) {
		ez->err = err; 
		ez->pe = pn-1; ez->te = tn-1;
	}
    return;
}

inline void ed_band_cal_semi_128bit(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{
    // (*re_err) = INT32_MAX;
	init_base_ed(*ez, thre, pn, tn); ez->ps = ez->pe = -1; ez->ts = 0; ez->te = tn-1;
	w128_clear(ez->VP); w128_clear(ez->VN); w128_clear(ez->mm); w128_bit(ez->mm, 0);
	w128_clear(ez->Peq[0]); w128_clear(ez->Peq[1]); w128_clear(ez->Peq[2]); w128_clear(ez->Peq[3]); w128_clear(ez->Peq[4]);
	int32_t bd = (thre<<1)+1, i, err = 0, i_bd = (thre<<1), last_high = (thre<<1), tn0 = tn - 1;
	int32_t cut = thre+last_high; uint8_t c;

    for (i = 0; i < bd; i++) {
		w128_self_or(ez->Peq[seq_nt4_table[(uint8_t)pstr[i]]], ez->mm); w128_self_lsft_1(ez->mm);
        // Peq[seq_nt4_table[(uint8_t)pstr[i]]] |= mm; mm <<= 1;
    }
	///should make Peq[4] = 0 if N is always an error
	// Peq[4] = 0;
	w128_clear(ez->Peq[4]);
	//mm = ((Word)1 << (thre<<1));///for the incoming char/last char
	w128_clear(ez->mm); w128_bit(ez->mm, (thre<<1));

	i = 0;
    while (i < tn0) {
		ed_core_w128(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP);
		if (!(ez->D0.a[0]&(1ULL))) {
            ++err; if (err>cut) return;
        }
		ed_core_w128_reshift(ez->Peq);

        ++i; ++i_bd; 
		// Peq[seq_nt4_table[(uint8_t)pstr[i_bd]]] |= mm; Peq[4] = 0;
		c = seq_nt4_table[(uint8_t)pstr[i_bd]]; 
		if(c < 4) w128_self_or(ez->Peq[c], ez->mm);
    }
    ed_core_w128(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP);
	if (!(ez->D0.a[0]&(1ULL))) {
		++err; if (err>cut) return;
	}

    int32_t site = tn - 1;///up bound
    ///in most cases, ai = (thre<<1)
    int32_t ai = pn - tn, uge = INT32_MAX;
    if ((err <= thre) && (err <= ez->err)) {
        ez->err = err; ez->pe = site;
    }
    i = 0;

    while (i < ai) {
        // err += ((VP >> i)&(1ULL)); 
		err += ez->VP.a[0]&(1ULL); w128_self_rsft_1(ez->VP); 
		// err -= ((VN >> i)&(1ULL)); 
		err -= ez->VN.a[0]&(1ULL); w128_self_rsft_1(ez->VN); 
		++i;
        if ((err <= thre) && (err <= ez->err)) {
            ez->err = err; ez->pe = site + i;
        }
        if(i == thre) uge = err;
    }

    if((uge <= thre) && (uge == ez->err)) ez->pe = site + thre;
}

///require:: (pn >= tn - thre && pn <= tn + thre)
inline void ed_band_cal_extension_128bit(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{
	// fprintf(stderr, "\n[M::%s::] pn::%d, tn::%d, thre::%d\n", __func__, pn, tn, thre);
	if(pn > tn + thre) pn = tn + thre;
	else if(tn > pn + thre) tn = pn + thre;
	// fprintf(stderr, "[M::%s::] pn::%d, tn::%d, thre::%d\n", __func__, pn, tn, thre);
	init_base_ed(*ez, thre, pn, tn); ez->ps = ez->ts = 0; ez->pe = ez->te = -1;
	// if((pn > tn + thre) || (tn > pn + thre)) return;
	// if((pn < thre + 1) || (tn < thre + 1)) return;
	int32_t i, err, tn0 = tn - 1, cut = thre+(thre<<1), bd = thre+1, i_bd = thre; uint8_t c;
	int32_t poff, pe = pn-1, tmp_e, k;
	w128_clear(ez->Peq[0]); w128_clear(ez->Peq[1]); w128_clear(ez->Peq[2]); w128_clear(ez->Peq[3]); w128_clear(ez->Peq[4]);

	w128_clear(ez->mm); w128_bit(ez->mm, thre); ///mm = (((Word)1)<<thre)
    for (i = 0; i < bd && i < pn; i++) {
		w128_self_or(ez->Peq[seq_nt4_table[(uint8_t)pstr[i]]], ez->mm); w128_self_lsft_1(ez->mm);
        // Peq[seq_nt4_table[(uint8_t)pstr[i]]] |= mm; mm <<= 1;
    }
	w128_clear(ez->Peq[4]);
	err = thre;
	w128_set_bit_lsub(ez->VN, thre); ///VN = (((Word)1)<<(thre))-1; 
	w128_set_bit_lsub(ez->VP, (thre<<1)+1); ///VP = (((Word)1)<<((thre<<1)+1))-1;
	w128_self_xor(ez->VP, ez->VN); ///VP ^= VN;

	// print_bits(ez->Peq[0].a, (thre<<1)+1, "-Peq[A]");
	// print_bits(ez->Peq[1].a, (thre<<1)+1, "-Peq[C]");
	// print_bits(ez->Peq[2].a, (thre<<1)+1, "-Peq[G]");
	// print_bits(ez->Peq[3].a, (thre<<1)+1, "-Peq[T]");
	// print_bits(ez->VN.a, (thre<<1)+1, "-VN");
	// print_bits(ez->VP.a, (thre<<1)+1, "-VP");

	///should make Peq[4] = 0 if N is always an error
	i = 0; 
	///for the incoming char/last char
	w128_clear(ez->mm); w128_bit(ez->mm, (thre<<1)); ///mm = ((Word)1 << (thre<<1));
	//VP + ((Peq|VN)&VP)
    while (i < tn0) {
		ed_core_w128(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP);
		if (!(ez->D0.a[0]&(1ULL))) {
            ++err; if (err>cut) return;
        }
		poff = i-thre; tmp_e = err; //poff:[i-thre, i+thre]
		tst_band_err(poff, i, pe, tmp_e, *ez, k);

		ed_core_w128_reshift(ez->Peq);
        ++i; ++i_bd;
        if(i_bd < pn) {
			c = seq_nt4_table[(uint8_t)pstr[i_bd]]; 
			if(c < 4) w128_self_or(ez->Peq[c], ez->mm);
		}
    }
	ed_core_w128(ez->Peq, ez->VP, ez->VN, ez->X, ez->D0, ez->HN, ez->HP);
	if (!(ez->D0.a[0]&(1ULL))) {
		++err; if (err>cut) return;
	}
	///i = tn - 1
    int32_t site = tn - 1 - thre;///up bound; site:[tn - 1 - thre, tn - 1 + thre]
	for (cut = pn - 1; site < cut; ) {
		// err += ((VP >> i)&(1ULL)); 
		err += ez->VP.a[0]&(1ULL); w128_self_rsft_1(ez->VP); 
		// err -= ((VN >> i)&(1ULL));
		err -= ez->VN.a[0]&(1ULL); w128_self_rsft_1(ez->VN); 
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
