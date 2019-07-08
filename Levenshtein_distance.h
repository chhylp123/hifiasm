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

typedef uint64_t Word;
typedef uint32_t Word_32;





/**
 pattern是长的那个，是y
 p_length是长的那个的长度,  p_length实际没用
 text是短的那个，是x
 t_length是短的那个的长度
 errthold是阈值
 return_err是编辑距离
 返回值是结束位置
 **/
inline int Reserve_Banded_BPM
(char *pattern, int p_length, char *text, int t_length, unsigned short errthold, unsigned int* return_err)
{
	(*return_err) = (unsigned int)-1;

	Word Peq[256];

	int band_length = (errthold << 1) + 1;
	int i = 0;
	Word tmp_Peq_1 = (Word)1;

	Peq['A'] = (Word)0;
	Peq['T'] = (Word)0;
	Peq['G'] = (Word)0;
	Peq['C'] = (Word)0;


	Word Peq_A;
	Word Peq_T;
	Word Peq_C;
	Word Peq_G;

	///band_length = 2k + 1
	for (i = 0; i<band_length; i++)
	{
		Peq[pattern[i]] = Peq[pattern[i]] | tmp_Peq_1;
		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	///Peq['T'] = Peq['T'] | Peq['C'];

	Peq_A = Peq['A'];
	Peq_C = Peq['C'];
	Peq_T = Peq['T'];
	Peq_G = Peq['G'];


	memset(Peq, 0, sizeof(Word)* 256);


	Peq['A'] = Peq_A;
	Peq['C'] = Peq_C;
	Peq['T'] = Peq_T;
	Peq['G'] = Peq_G;


	

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
		///pattern[0]ÔÚPeq[2k], ¶øpattern[2k]ÔÚPeq[0]
		X = Peq[text[i]] | VN;

		D0 = ((VP + (X&VP)) ^ VP) | X;

		HN = VP&D0;
		HP = VN | ~(VP | D0);

		X = D0 >> 1;
		VN = X&HP;
		VP = HN | ~(X | HP);

		if (!(D0&err_mask))
		{
			++err;

			///¼´Ê¹È«²¿µÝ¼õ£¬Ò²¾Í¼õ2k
			if ((err - last_high)>errthold)
			{
				///fprintf(stderr, "0 ######, i: %u\n", i);
				return -1;
			}
				
		}


		Peq['A'] = Peq['A'] >> 1;
		Peq['C'] = Peq['C'] >> 1;
		Peq['G'] = Peq['G'] >> 1;
		Peq['T'] = Peq['T'] >> 1;


		++i;
		++i_bd;
		Peq[pattern[i_bd]] = Peq[pattern[i_bd]] | Mask;


		///Peq['T'] = Peq['T'] | Peq['C'];
	}





	X = Peq[text[i]] | VN;
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
	if ((err <= errthold) && (err<=*return_err))
	{
		*return_err = err;
		return_site = site;
	}
	int i_last = i;
	i = 0;




	while (i<errthold)
	{
		err = err + ((VP >> i)&(Word)1);
		err = err - ((VN >> i)&(Word)1);
		++i;

		if ((err <= errthold) && (err <= *return_err))
		{
			*return_err = err;
			return_site = site + i;
		}
	}


	unsigned int ungap_err;
	ungap_err = err;


	while (i<last_high)
	{
		err = err + ((VP >> i)&(Word)1);
		err = err - ((VN >> i)&(Word)1);
		++i;

		if ((err <= errthold) && (err<=*return_err))
		{
			*return_err = err;
			return_site = site + i;
		}




	}


	if ((ungap_err <= errthold) && (ungap_err == *return_err))
	{
		return_site = site + errthold;
	}

	return return_site;

}

inline int Reserve_Banded_BPM_4_SSE_only(char *pattern1, char *pattern2, char *pattern3, char *pattern4, int p_length, char *text, int t_length,
	int* return_sites, unsigned int* return_sites_error, unsigned short errthold, __m128i* Peq_SSE)

{
	memset(return_sites, -1, sizeof(int)* 4);
	memset(return_sites_error, -1, sizeof(unsigned int)* 4);

	Word_32 Peq[256][4];
	int band_length = (errthold << 1) + 1;


	int i;

	Word_32 tmp_Peq_1 = 1;


	memset(Peq['A'], 0, sizeof(Word_32)* 4);
	memset(Peq['C'], 0, sizeof(Word_32)* 4);
	memset(Peq['G'], 0, sizeof(Word_32)* 4);
	memset(Peq['T'], 0, sizeof(Word_32)* 4);

	for (i = 0; i<band_length; i++)
	{
		Peq[pattern1[i]][0] = Peq[pattern1[i]][0] | tmp_Peq_1;
		Peq[pattern2[i]][1] = Peq[pattern2[i]][1] | tmp_Peq_1;
		Peq[pattern3[i]][2] = Peq[pattern3[i]][2] | tmp_Peq_1;
		Peq[pattern4[i]][3] = Peq[pattern4[i]][3] | tmp_Peq_1;

		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	Peq_SSE['A'] = _mm_set_epi32(Peq['A'][3], Peq['A'][2], Peq['A'][1], Peq['A'][0]);
	Peq_SSE['C'] = _mm_set_epi32(Peq['C'][3], Peq['C'][2], Peq['C'][1], Peq['C'][0]);
	Peq_SSE['G'] = _mm_set_epi32(Peq['G'][3], Peq['G'][2], Peq['G'][1], Peq['G'][0]);
	Peq_SSE['T'] = _mm_set_epi32(Peq['T'][3], Peq['T'][2], Peq['T'][1], Peq['T'][0]);
	///Peq_SSE['T'] = _mm_or_si128(Peq_SSE['T'], Peq_SSE['C']);

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
		X = _mm_or_si128(Peq_SSE[text[i]], VN);



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


		Peq_SSE['A'] = _mm_srli_epi32(Peq_SSE['A'], 1);
		Peq_SSE['T'] = _mm_srli_epi32(Peq_SSE['T'], 1);
		Peq_SSE['G'] = _mm_srli_epi32(Peq_SSE['G'], 1);
		Peq_SSE['C'] = _mm_srli_epi32(Peq_SSE['C'], 1);

		++i;
		++i_bd;

		Peq_SSE[pattern1[i_bd]] = _mm_or_si128(Mask1, Peq_SSE[pattern1[i_bd]]);
		Peq_SSE[pattern2[i_bd]] = _mm_or_si128(Mask2, Peq_SSE[pattern2[i_bd]]);
		Peq_SSE[pattern3[i_bd]] = _mm_or_si128(Mask3, Peq_SSE[pattern3[i_bd]]);
		Peq_SSE[pattern4[i_bd]] = _mm_or_si128(Mask4, Peq_SSE[pattern4[i_bd]]);
		///Peq_SSE['T'] = _mm_or_si128(Peq_SSE['T'], Peq_SSE['C']);
	}



		///X = Peq[text[i]] | VN;
	X = _mm_or_si128(Peq_SSE[text[i]], VN);

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


	if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
	{
		return_sites[0] = site;
		return_sites_error[0] = err1;
	}
	if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
	{
		return_sites[1] = site;
		return_sites_error[1] = err2;
	}
	if ((err3 <= errthold) && (err3 <= return_sites_error[2]))
	{
		return_sites[2] = site;
		return_sites_error[2] = err3;
	}
	if ((err4 <= errthold) && (err4 <= return_sites_error[3]))
	{
		return_sites[3] = site;
		return_sites_error[3] = err4;
	}


	i = 0;



	while (i<errthold)
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


		if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
		{
			return_sites[0] = site + i;
			return_sites_error[0] = err1;
		}
		if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
		{
			return_sites[1] = site + i;
			return_sites_error[1] = err2;
		}
		if ((err3 <= errthold) && (err3 <= return_sites_error[2]))
		{
			return_sites[2] = site + i;
			return_sites_error[2] = err3;
		}
		if ((err4 <= errthold) && (err4 <= return_sites_error[3]))
		{
			return_sites[3] = site + i;
			return_sites_error[3] = err4;
		}
	}


	unsigned int ungap_err1;
	unsigned int ungap_err2;
	unsigned int ungap_err3;
	unsigned int ungap_err4;
	ungap_err1 = err1;
	ungap_err2 = err2;
	ungap_err3 = err3;
	ungap_err4 = err4;




	while (i<last_high)
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


		if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
		{
			return_sites[0] = site + i;
			return_sites_error[0] = err1;
		}
		if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
		{
			return_sites[1] = site + i;
			return_sites_error[1] = err2;
		}
		if ((err3 <= errthold) && (err3 <= return_sites_error[2]))
		{
			return_sites[2] = site + i;
			return_sites_error[2] = err3;
		}
		if ((err4 <= errthold) && (err4 <= return_sites_error[3]))
		{
			return_sites[3] = site + i;
			return_sites_error[3] = err4;
		}
	}


	if ((ungap_err1 <= errthold) && ungap_err1 == return_sites_error[0])
	{
		return_sites[0] = site + errthold;
	}


	if ((ungap_err2 <= errthold) && ungap_err2 == return_sites_error[1])
	{
		return_sites[1] = site + errthold;
	}


	if ((ungap_err3 <= errthold) && ungap_err3 == return_sites_error[2])
	{
		return_sites[2] = site + errthold;
	}


	if ((ungap_err4 <= errthold) && ungap_err4 == return_sites_error[3])
	{
		return_sites[3] = site + errthold;
	}

	return 1;
	

}


void output_bit_myers(Word x, int length);
void prase_vertical(Word VP, Word VN, int length, int matrix[1000][1000], int i);
void prase_D0(Word D0, int length, int matrix[1000][1000], int i);
void prase_H(Word HP, Word HN, int length, int matrix[1000][1000], int i);
int Reserve_Banded_BPM_debug(char *pattern, int p_length, char *text, int t_length, unsigned short errthold, 
unsigned int* return_err, int matrix[1000][1000]);
int Reserve_Banded_BPM_new(char *pattern,int p_length,char *text,int t_length,unsigned short errthold,
unsigned short band_down,unsigned short band_below,unsigned short band_length,int* return_err, int thread_id);
int BS_Reserve_Banded_BPM
(char *pattern, int p_length, char *text, int t_length, unsigned short errthold, unsigned int* return_err);

#endif