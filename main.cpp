#include <stdio.h>
#include <stdlib.h>
#include "CommandLines.h"
#include "Process_Read.h"
#include "Assembly.h"
#include "Levenshtein_distance.h"
#include "edlib.h"
/********************************for debug***************************************/
///使用这个函数的时候，必须把Counting_multiple_thr()里的destory_Total_Count_Table(&TCB)注释掉
void debug_Counting()
{
    init_kseq(read_file_name);
    Verify_Counting();
    fprintf(stderr, "debug over!\n");
    destory_kseq();
}


int matrix[1000][1000] = {0};
///y_length > x_length
int edit_distance_normal(char* y, int y_length, char* x, int x_length)
{   memset(matrix, 0, sizeof(matrix));

    int i, j;
    for (i = 0; i <= x_length; i++)
    {
        matrix[i][0] = i;
    }

    int digonal, up, left, min;

    ///一列列算的
    for (i = 0; i < x_length; i++)
    {
        for (j = 0; j < y_length; j++)
        {
            ///matrix[i + 1][j + 1]
            digonal = matrix[i][j] + (x[i] !=  y[j]);
            up = matrix[i + 1][j] + 1;
            left = matrix[i][j + 1] + 1;
            min = digonal;
            if (up < min)
            {
                min = up;
            }

            if (left< min)
            {
                min = left;
            }

            matrix[i + 1][j + 1] = min;
        }
    }

    min = 999999;
    for (j = 0; j <= y_length; j++)
    {
        if (matrix[i][j] < min)
        {
            min = matrix[i][j];
        }
    }


    return min;
}


///y_length > x_length
int edit_distance_normal_banded(char* y, int y_length, char* x, int x_length, int error)
{   memset(matrix, 0, sizeof(matrix));

    int i, j;
    for (i = 0; i <= x_length; i++)
    {
        for (j = 0; j <= y_length; j++)
        {
            matrix[i][j] = 1000000;
        }
    }

    for (i = 0; i <= x_length; i++)
    {
        matrix[i][0] = i;
    }

    for (i = 0; i <= y_length; i++)
    {
        matrix[0][i] = 0;
    }

    int banded_length = error*2 + 1;

    int digonal, up, left, min;


    for (i = 0; i < x_length; i++)
    {
        ///for (j = 0; j < y_length; j++)
        for (j = i; j < banded_length + i; j++)
        {
            ///matrix[i + 1][j + 1]
            digonal = matrix[i][j] + (x[i] !=  y[j]);
            up = matrix[i + 1][j] + 1;
            left = matrix[i][j + 1] + 1;
            min = digonal;
            if (up < min)
            {
                min = up;
            }

            if (left< min)
            {
                min = left;
            }

            matrix[i + 1][j + 1] = min;
        }
    }
    min = 999999;
    for (j = 0; j <= y_length; j++)
    ///for (j = i; j < banded_length + i; j++)
    {
        if (matrix[i][j] < min)
        {
            min = matrix[i][j];
        }
    }


    return min;
}

void debug_edit_distance()
{
    /**
    char* x = "TTCCATACGATTCCATTCAATTCGAGACCATTCTATTCCTGTCCATTCCTTGTGGTTCGATTCCATTTCACTCTAGTCCATTCCATTCCATTCAATTCCATTCGACTCTATTCCGTTCCACTCAATTCCATTCCATTCGATTCCATTTTTTTCGAGAACCTTCCATTACACTCCCTTCCATTCCAGTGCATTCCATTCCAGTCTCTTCAGTTCGATTCCATTCCATTCGTTTCGATTCCTTTCCATTCCAGCCCATTCCATTCCATTCCATTCCTTTCCTTTCCGTTTCATTAGATTCCATTGCATTCGATTCCATTCAAATCAATTCCGTTCTATTCAATTTGATTCAT";
    char* y = "CCATACGATTCCATTCAATTCGAGACCATTCTATTCCTGTCCATTCCTTGTGGTTCGATTCCATTTCACTCTAGTCCATTCCATTCCATTCAATTCCATTCGACTCTATTCCGTTCCATTCAATTCCATTCCATTCGATTCCATTTTTTTCGAGAACCTTCCATTACACTCCCTTCCATTCCAGTGCATTCCATTCCAGTCTCTTCACTTCGATTCCATTCCATTCGTTTCGATTCCTTTCCATTCCAGCCCATTCCATTCCATTCCATTCCTTTCCTTTCCGTTTCATTAGATTCCATTGCATTCCATTCCATTCAATTCAATTCCGTGCTATTCAATTTGATTCATTTCCATTTAATTCCATTCCATTAGATTCCATT";
    **/
   unsigned short toold = 15;
   char* x 
   = "GAAAGAGAATCAAATGGAATTGAATCGAATGGAATCGAATGGATTGGAAAGGAATAGAATGGAATGGAATGGAATTGACTCAAATGGAATGGACTAGAATGGAATGGATTCGAATGGAAGGCAAAGGAATGGAATCTATCGGAATGGACTGTAATGGAATGGAATGGAAGGGATTGGAATGGATTCGAATGTAATGGACTGCAATAGAAAGGATTCGAATGGAATGAAAAAGAATTGAATGGAATAGAACAGAATGGAATCAAATCGAAGGAAATGGAATGGAATAGAAAGGAATGGAATGAAATGGAATGGAAAGGATTCGAATGGAATGCAATCGAATGGAATGGAATCGAACGGAATGGAATAAAATGGAAG";
    char* y = 
     "GAAAGAGAATCAAATGGAATTGAATCGAATGGAATCGAATGGATTGGAAAGGAATAGAATGGAATGGAATGGAATGGACTCAAATGGAATGTACTAGAATGGAATGGATTCGAATGGAAGGCAAAGGAATGGAATCTATTGGAATGGACTGTAATGGAATGGAATGGAAGGGATTGGAATGGACTCGAATGGAATGGACTGCAATAGAAAGGATTCGAATGGAATGAAAAAGAATTGAATGGAATAGAACAGAATGGAATCAAATCGAATGAAATGGAATGGAATAGAAAGGAATGGAATGAAATGGAATGGAAAGGATTCGAATGGAATGCAATCGAATGGAATGGAATCGAACGGAATGGAATAAATTTTCTG";
    fprintf(stderr, "x_length: %u\n", strlen(x));
    fprintf(stderr, "y_length: %u\n", strlen(y));


    EdlibAlignResult result = edlibAlign(x, strlen(x), y, strlen(y), 
            edlibNewAlignConfig(toold, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

    if (result.status == EDLIB_STATUS_OK) {

        fprintf(stderr, "****\nedlib: %d, alignmentLength: %d, startLocations: %d, endLocations: %d\n", 
        result.editDistance, result.alignmentLength, result.startLocations[0], result.endLocations[0]);
        char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
        fprintf(stderr,"%s\n", cigar);
        free(cigar);
    }
    edlibFreeAlignResult(result);


    
    unsigned int error;
    int end_site = Reserve_Banded_BPM(y, strlen(y), x, strlen(x), toold, &error);

    fprintf(stderr, "BPM: error: %u, end_site: %u\n", error, end_site);


    unsigned short band_length=(toold+1)*3-1-1-toold;
    unsigned short band_down=toold-1;
    unsigned short band_blew=2*(toold+1)-1-1;



    int return_err = 99999;
    ///注意pattern/text和band_down/band_blew是反的
    Reserve_Banded_BPM_new(y, strlen(y), x, strlen(x),
                toold,band_blew,band_down,band_length, &return_err, 0);

    fprintf(stderr, "new BPM: error: %u\n", return_err);


    return_err = edit_distance_normal(y, strlen(y), x, strlen(x));
    fprintf(stderr, "edit_distance_normal: error: %u\n", return_err);

    return_err = edit_distance_normal_banded(y, strlen(y), x, strlen(x), toold);
    fprintf(stderr, "edit_distance_normal_banded: error: %u\n", return_err);

    end_site = Reserve_Banded_BPM_debug(y, strlen(y), x, strlen(x), toold, &error, matrix);

    fprintf(stderr, "BPM debug: error: %u, end_site: %u\n", error, end_site);
}


int main(int argc, char *argv[])
{

    if (!CommandLine_process(argc, argv))
        return 1;


    /**
    debug_edit_distance();

    return 1;
    **/
    



    fprintf(stdout, "k-mer length: %d\n",k_mer_length);

    if (load_index_from_disk && load_pre_cauculated_index())
    {
        ;
    }
    else
    {
        Counting_multiple_thr();

        Build_hash_table_multiple_thr();
    }
    
    ///verify_Position_hash_table();

    Overlap_calculate_multipe_thr();



    return 1;
}
