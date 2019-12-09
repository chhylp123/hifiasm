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
	///p_length大部分情况下应该是t_length + 2 * errthold，这是i要小于last_high = 2 * errthold
	///也就是p_length - t_length
	///那么当p_length < t_length + 2 * errthold, available_i也应该是这个值
	int available_i = p_length - t_length;
	if ((err <= errthold) && (err<=*return_err))
	{
		*return_err = err;
		return_site = site;
	}
	int i_last = i;
	i = 0;

	/****************************may have bugs********************************/
	unsigned int ungap_error = (unsigned int)-1;
	/****************************may have bugs********************************/

	while (i < available_i)
	{
		err = err + ((VP >> i)&(Word)1);
		err = err - ((VN >> i)&(Word)1);
		++i;

		if ((err <= errthold) && (err <= *return_err))
		{
			*return_err = err;
			return_site = site + i;
		}

		/****************************may have bugs********************************/
		if(i == errthold)
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
	///y上的起始位置
	int start_site = end_site - t_length + 1;

	if (start_site >= 0)
	{

		for (i = 0; i < t_length; i++)
		{
			///path[i] = 0;
			///path倒着存
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


///p_length有可能不够，但是t_length总是够的
///就是p_length有可能小于t_length + 2 * errthold
inline int Reserve_Banded_BPM_PATH
(char *pattern, int p_length, char *text, int t_length, unsigned short errthold, 
		unsigned int* return_err, int* return_start_site, int* return_path_length, Word* matrix_bit, char* path, 
		int old_error, int old_end_site)
{
	if (old_error != -1 && old_end_site != -1)
	{
		if (old_error == 0)
		{
			///fprintf(stderr, "0 error\n");
			(*return_err) = old_error;
			(*return_start_site) = old_end_site - t_length + 1;
			return old_end_site;
		}

		if (try_cigar(pattern, p_length, text, t_length, old_end_site, path,
		old_error, return_start_site, return_path_length))
		{
			///fprintf(stderr, "no gap error\n");
			(*return_err) = old_error;
			return old_end_site;
		}
		
	}
	
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
			{
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

		column_start = i << 3;
		matrix_bit[column_start] = D0;
		matrix_bit[column_start + 1] = VP;
		matrix_bit[column_start + 2] = VN;
		matrix_bit[column_start + 3] = HP;
		matrix_bit[column_start + 4] = HN;
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

	///p_length大部分情况下应该是t_length + 2 * errthold，这是i要小于last_high = 2 * errthold
	///也就是p_length - t_length
	///那么当p_length < t_length + 2 * errthold, available_i也应该是这个值
	int available_i = p_length - t_length;
	if ((err <= errthold) && (err<=*return_err))
	{
		*return_err = err;
		return_site = site;
	}
	int i_last = i;
	i = 0;

	while (i < available_i)
	{
		err = err + ((VP >> i)&(Word)1);
		err = err - ((VN >> i)&(Word)1);
		++i;

		if ((err <= errthold) && (err <= *return_err))
		{
			*return_err = err;
			return_site = site + i;
		}

		/****************************may have bugs********************************/
		if(i == errthold)
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
	
    ////注意，这里p_length要矫正啊啊
	///不矫正会出错
	///因为p_length有可能不够
	p_length = t_length + 2 * errthold;
	///end_site是正确的
	int end_site = return_site;
	int start_site = end_site;
	///这个是各个bit-vector里面，end_site对应bit所在的位置
	int back_track_site = band_length - (p_length - end_site);

	Word v_value, h_value, delta_value, min_value, current_value;
	Word direction, is_mismatch; ///0 is match, 1 is mismatch, 2 is up, 3 is left

	///代表pattern到哪了，就是短的那个到哪了
	i = t_length;
	int path_length = 0;

	///到0就结束了，后面的路径可以直接match
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


////这个p_length四个是一样的
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

	/****************************may have bugs********************************/
	unsigned int ungap_error1 = (unsigned int)-1;
	unsigned int ungap_error2 = (unsigned int)-1;
	unsigned int ungap_error3 = (unsigned int)-1;
	unsigned int ungap_error4 = (unsigned int)-1;
	/****************************may have bugs********************************/


	///p_length大部分情况下应该是t_length + 2 * errthold，这是i要小于last_high = 2 * errthold
	///也就是p_length - t_length
	///那么当p_length < t_length + 2 * errthold, available_i也应该是这个值
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

		/****************************may have bugs********************************/
		if(i == errthold)
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