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




inline void output_bit_myers(Word x, int length)
{
	int i = 0;
	while (i < length)
	{
		fprintf(stderr, "%u", (x >> i) & ((Word)1));
		i++;
	}
	fprintf(stderr, "\n");
}


inline void prase_vertical(Word VP, Word VN, int length, int matrix[1000][1000], int i)
{

	i++;
	int k = 0;
	int j = i;
	int diff;
	while (k < length)
	{
		
		int x_p = (VP >> k) & ((Word)1);
		int x_n = (VN >> k) & ((Word)1);
		if (x_p == 1 && x_n == 1)
		{
			fprintf(stderr, "error\n");
		}
		
		
		
		if (x_p == 1)
		{
			diff = 1;
			///fprintf(stderr, "[+1]");
		}
		
		if (x_n == 1)
		{
			diff = -1;
			///fprintf(stderr, "[-1]");
		}

		if (x_p == 0 && x_n == 0)
		{
			diff = 0;
			///fprintf(stderr, "[+0]");
		}
		j++;
		if (matrix[i][j] - matrix[i][j - 1] != diff)
		{
			fprintf(stderr, "*************V(k): %u\n", k);
			///return;
		}
		k++;
	}
	///fprintf(stderr, "\n");
}




inline void prase_D0(Word D0, int length, int matrix[1000][1000], int i)
{

	i++;
	int k = 0;
	int j = i;
	int diff;
	while (k < length)
	{
		
		int diff = (D0 >> k) & ((Word)1);

		if(diff==matrix[i][j] - matrix[i-1][j-1])
		{
			fprintf(stderr, "*************D(k): %u\n", k);
			fprintf(stderr, "diff: %u, matrix[i][j]: %u, matrix[i-1][j-1]: %u\n", diff, matrix[i][j], matrix[i-1][j-1]);
			///return;
		}
		j++;
		k++;
	}
	///fprintf(stderr, "\n");
}

inline void prase_H(Word HP, Word HN, int length, int matrix[1000][1000], int i)
{

	i++;
	int k = 0;
	int j = i;
	int diff;
	while (k < length)
	{
		
		int x_p = (HP >> k) & ((Word)1);
		int x_n = (HN >> k) & ((Word)1);
		if (x_p == 1 && x_n == 1)
		{
			fprintf(stderr, "error\n");
		}
		if (x_p == 1)
		{
			diff = 1;
			///fprintf(stderr, "[+1]");
		}
		if (x_n == 1)
		{
			diff = -1;
			///fprintf(stderr, "[-1]");
		}
		if (x_p == 0 && x_n == 0)
		{
			diff = 0;
			///fprintf(stderr, "[+0]");
		}

		if(diff!=matrix[i][j] - matrix[i-1][j])
		{
			fprintf(stderr, "*************H(k): %u\n", k);
			///return;
		}

		j++;
		k++;
	}
	///fprintf(stderr, "\n");
}


/**
 pattern是长的那个，是y
 p_length是长的那个的长度,  p_length实际没用
 text是短的那个，是x
 t_length是短的那个的长度
 errthold是阈值
 return_err是编辑距离
 返回值是结束位置
 **/
inline int Reserve_Banded_BPM_debug
(char *pattern, int p_length, char *text, int t_length, unsigned short errthold, unsigned int* return_err, int matrix[1000][1000])
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
		fprintf(stderr, "i: %u, j_s: %u, j_e: %u\n", i + 1, i + 1, i + band_length + 1);
		///if (i >= 6)
		if (i >= 0)
		{
			/**
			fprintf(stderr, "VP:\n");
			output_bit_myers(VP, band_length);

			fprintf(stderr, "VN:\n");
			output_bit_myers(VN, band_length);

			fprintf(stderr, "HP:\n");
			output_bit_myers(HP, band_length);

			fprintf(stderr, "HN:\n");
			output_bit_myers(HN, band_length);

			fprintf(stderr, "D0:\n");
			output_bit_myers(D0, band_length);

			fprintf(stderr, "text[i]: %c\n", text[i]);

			fprintf(stderr, "Peq[text[i]]:\n");
			output_bit_myers(Peq[text[i]], band_length);

			for (size_t j = i; j < i + band_length; j++)
			{
				fprintf(stderr, "%c", pattern[j]);
			}
			fprintf(stderr, "\n");
			
			fprintf(stderr, "Previous begin.\n");
			prase_vertical(VP, VN, band_length, matrix, i-1);
			prase_D0(D0, band_length, matrix, i-1);
			prase_H(HP, HN, band_length, matrix, i-1);
			fprintf(stderr, "Previous test done.\n");
			**/

			X = Peq[text[i]] | VN;
			/**
			fprintf(stderr, "#X:\n");
			output_bit_myers(X, band_length);
			**/

			D0 = ((VP + (X&VP)) ^ VP) | X;
			/**
			fprintf(stderr, "#(X&VP):\n");
			output_bit_myers((X&VP), band_length);

			fprintf(stderr, "#(VP + (X&VP)):\n");
			output_bit_myers((VP + (X&VP)), band_length);

			fprintf(stderr, "#((VP + (X&VP)) ^ VP):\n");
			output_bit_myers(((VP + (X&VP)) ^ VP), band_length);

			fprintf(stderr, "#D0:\n");
			output_bit_myers(D0, band_length);
			**/
			HN = VP&D0;
			HP = VN | ~(VP | D0);

			X = D0 >> 1;
			VN = X&HP;
			VP = HN | ~(X | HP);
		}
		else
		{
			///pattern[0]ÔÚPeq[2k], ¶øpattern[2k]ÔÚPeq[0]
			X = Peq[text[i]] | VN;

			D0 = ((VP + (X&VP)) ^ VP) | X;

			HN = VP&D0;
			HP = VN | ~(VP | D0);

			X = D0 >> 1;
			VN = X&HP;
			VP = HN | ~(X | HP);
		}
		

		

		
		
		
		/**
		for (size_t j = i + 1; j <= i + band_length + 1; j++)
		{
			fprintf(stderr, "[%u]", matrix[i + 1][j]);
		}
		fprintf(stderr, "\n");
		**/

		prase_vertical(VP, VN, band_length, matrix, i);
		prase_D0(D0, band_length, matrix, i);
		prase_H(HP, HN, band_length, matrix, i);
		/**
		fprintf(stderr, "VP:\n");
		output_bit_myers(VP, band_length);

		fprintf(stderr, "VN:\n");
		output_bit_myers(VN, band_length);
		**/
		

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

	prase_vertical(VP, VN, band_length, matrix, i);
	prase_D0(D0, band_length, matrix, i);
	prase_H(HP, HN, band_length, matrix, i);


	fprintf(stderr, "err: %d, matrix[][]: %d\n", err, matrix[i+1][i+1]);
	fprintf(stderr, "VP:\n");
	output_bit_myers(VP, band_length);

	fprintf(stderr, "VN:\n");
	output_bit_myers(VN, band_length);

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

		fprintf(stderr, "*i: %u, err: %d\n", i, err);

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

		fprintf(stderr, "*i: %u, err: %d\n", i, err);

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








inline int BS_Reserve_Banded_BPM
(char *pattern, int p_length, char *text, int t_length, unsigned short errthold, unsigned int* return_err)
{
	(*return_err) = (unsigned int)-1;

	///Õâ¸öÊÇÄÇ¸öÐèÒªÔ¤´¦ÀíµÄÏòÁ¿
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
	///ÕâÊÇ°ÑpatternµÄÇ°2k + 1¸ö×Ö·ûÔ¤´¦Àí
	///pattern[0]¶ÔÓ¦Peq[0]
	///pattern[2k]¶ÔÓ¦Peq[2k]
	for (i = 0; i<band_length; i++)
	{
		Peq[pattern[i]] = Peq[pattern[i]] | tmp_Peq_1;
		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	Peq['T'] = Peq['T'] | Peq['C'];

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
		///Èç¹ûÐ±¶Ô½ÇÏß·½ÏòÆ¥ÅäÔòD0ÊÇ1
		///Èç¹û²»Æ¥ÅäÔòD0ÊÇ0
		///Õâ¸öÒâË¼ÊÇÈç¹û×îÉÏÃæÄÇÌõ¶Ô½ÇÏßÉÏµÄÐ±¶Ô½ÇÏß·½Ïò·¢ÉúÎóÅä,ÔòÖ´ÐÐÄÚ²¿³ÌÐò
		///
		if (!(D0&err_mask))
		{
			++err;

			///¼´Ê¹È«²¿µÝ¼õ£¬Ò²¾Í¼õ2k
			if ((err - last_high)>errthold)
				return -1;
		}

		///pattern[0]ÔÚPeq[2k], ¶øpattern[2k]ÔÚPeq[0]
		//ÓÒÒÆÊµ¼ÊÉÏÊÇ°Ñpattern[0]ÒÆµôÁË
		Peq['A'] = Peq['A'] >> 1;
		Peq['C'] = Peq['C'] >> 1;
		Peq['G'] = Peq['G'] >> 1;
		Peq['T'] = Peq['T'] >> 1;


		++i;
		++i_bd;
		///ÕâÊÇ°ÑÐÂµÄpattern[2k]¼Ó½øÀ´, ÕâÃ²ËÆÊÇ¼Óµ½Peq[2k]ÉÏÁË
		Peq[pattern[i_bd]] = Peq[pattern[i_bd]] | Mask;


		Peq['T'] = Peq['T'] | Peq['C'];
	}




	///fprintf(stderr, "sucess(1)\n");


	///Õâ¸öÑ­»·ÄÃ³öÀ´ÊÇÎªÁË·ÀÖ¹ÄÚ´æÐ¹Â¶
	///ÆäÊµÒ²¾ÍÊÇÑ­»·ÀïµÄ×îºóÒ»ÐÐÓï¾ä°É
	///ÍêÈ«¿ÉÒÔ°ÑpatternÔö´óÒ»Î»
	///²»¹ýÕâÑùÒ²ºÃ£¬¿ÉÒÔ¼õÉÙ¼ÆËã¿ªÏú
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
	///´ËÊ±Õâ¸ösiteÃ²ËÆÊÇ×îÉÏÃæÄÇÌõ¶Ô½ÇÏßµÄÎ»ÖÃ
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



inline int Reserve_Banded_BPM_new(char *pattern,int p_length,char *text,int t_length,
                unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length,int* return_err, int thread_id)
{

    Word Peq[128];
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;
    Word tmp_Peq_1=(Word)1;


    Peq['A']=(Word)0;
    Peq['T']=(Word)0;
    Peq['G']=(Word)0;
    Peq['C']=(Word)0;
    Word Peq_A;
    Word Peq_T;
    Word Peq_C;
    Word Peq_G;

    for (r =0; r<band_length; r++)
    {
            Peq[pattern[r]]=Peq[pattern[r]]|tmp_Peq_1;
        tmp_Peq_1=tmp_Peq_1<<1;
    }
    Peq_A=Peq['A'];
    Peq_C=Peq['C'];
    Peq_T=Peq['T'];
    Peq_G=Peq['G'];

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq[symbol]=(Word)0;
        //jump_c[symbol]=128;
    }
    Peq['A']=Peq_A;
    Peq['C']=Peq_C;
    Peq['T']=Peq_T;
    Peq['G']=Peq_G;


    Word Mask_Pre=(Word)1<<(band_length-2);
    Word Mask=(Word)1<<(band_length-1);
    Word VP=0;
    Word VN=0;
    Word X=0;
    Word D0=0;
    Word HN=0;
    Word HP=0;


    int s=0;
    int i = 0;
    int j=0;

    int bound=band_length-2-band_down;
    int err=0;

    Word err_mask=(Word)1;
    int s1=band_length-2;
    int i_bd=i+band_down;
    int last_high=band_length-t_length+p_length-band_down-1;
   int t_length_1=t_length-1;
    //while(i<t_length)
    while(i<t_length_1)
    {

        X=Peq[text[i]]|VN;

        D0=((VP+(X&VP))^VP)|X;

        HN=VP&D0;
        HP=VN|~(VP|D0);

        X=D0>>1;
        VN=X&HP;
        VP=HN|~(X|HP);
        if(!(D0&err_mask))
         {
            ++err;
            if((err-last_high)>errthold)
                return -1;
        }

        Peq['A']=Peq['A']>>1;
        Peq['C']=Peq['C']>>1;
        Peq['G']=Peq['G']>>1;
        Peq['T']=Peq['T']>>1;


        ++i;
        ++i_bd;
        Peq[pattern[i_bd]]=Peq[pattern[i_bd]]|Mask;
    }


    ///这个循环拿出来是为了防止内存泄露

       X=Peq[text[i]]|VN;
        D0=((VP+(X&VP))^VP)|X;
        HN=VP&D0;
        HP=VN|~(VP|D0);
        X=D0>>1;
        VN=X&HP;
        VP=HN|~(X|HP);
        if(!(D0&err_mask))
         {
            ++err;
            if((err-last_high)>errthold)
                return -1;
        }


    int site=p_length-last_high-1;
    int return_site=-1;
    if((err<=errthold)&&(err<*return_err))
    {
            *return_err=err;
            return_site=site;
    }
    int i_last=i;
    i=0;
    while(i<last_high)
    {
        err=err+((VP>>i)&(Word)1);
        err=err-((VN>>i)&(Word)1);
        ++i;

        if((err<=errthold)&&(err<*return_err))
        {
                *return_err=err;
                return_site=site+i;
        }


    }
    return return_site;

}


#endif