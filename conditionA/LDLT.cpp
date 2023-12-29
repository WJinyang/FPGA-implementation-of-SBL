//
// Created by imvia on 2023/5/26.
//

/*
 * File: choleskyInverse.c
 *
 * MATLAB Coder version            : 5.3
 * C/C++ source code generated on  : 12-May-2023 12:35:52
 */

/* Include Files */
#include <iostream>
#include "Types.h"
#include <cmath>
#include <string.h>

/* Function Definitions */
/*
 * Arguments    : const complex_fixed matrix[64]
 *                complex_fixed inverse[64]
 * Return Type  : void
 */
void LDLT( complex_fixed matrix[M][M], complex_fixed L1[M][M], real_fixed D_inv[M], real_fixed alpha_sigma2_copy1[N+1], real_fixed alpha_sigma2_copy11[N+1], complex_fixed Y_copy1[M][T], complex_fixed Y_copy11[M][T])
{

    creal_T L[M][M];
#pragma HLS ARRAY_PARTITION variable=L type=complete dim=2
    double D[M];
    real_fixed array_re[M];
    real_fixed array_im[M];
    real_fixed sum_im;
    real_fixed sum_re;
    real_fixed ar;
    real_fixed bi;
    real_fixed bim;
    real_fixed im;
    real_fixed re;
    real_fixed sgnbi;
    real_fixed const1 = 1;
    int i;
    int j;
    int k;
    int m;


//#pragma HLS ARRAY_PARTITION variable=array_re type=complete dim=1
//#pragma HLS ARRAY_PARTITION variable=array_im type=complete dim=1
//#pragma HLS ARRAY_PARTITION variable=D type=complete dim=1

//    memset(&L[0][0], 0, M * M * sizeof(complex_fixed));
    /*  LDL*/
    D[0] = matrix[0][0].re;
    D_inv[0] = 1/D[0];
//    D_inv = const1/matrix[0][0].re;
//    L[0][0].re = 1;
    LDLT_col1:
    for (i = 0; i < M; i++) {
#pragma HLS PIPELINE
        L[i][0].re = double(matrix[i][0].re) / D[0];
        L[i][0].im = double(matrix[i][0].im) / D[0];
        L[i][i].re = 1;
        L[i][i].im = 0;
        alpha_sigma2_copy11[i] = alpha_sigma2_copy1[i];
    }


    LDLT_row:
    for (i = 1; i < M; i++) {
        LDLT_col:
        for (j = 1; j <= M-1; j++) {
#pragma HLS LOOP_FLATTEN
#pragma HLS PIPELINE
            sum_re = 0.0;
            sum_im = 0.0;
            memset(&array_re[0], 0, M * sizeof(real_fixed));
            memset(&array_im[0], 0, M * sizeof(real_fixed));

//            for(m = 0; m < M; m++)
//                std::cout << array[m] << " ";
            if (j == i) {
                switch (j) {
                    case 7 :
                        array_re[6] = (L[j][6].re * L[j][6].re + L[j][6].im * L[j][6].im) * D[6];
                    case 6 :
                        array_re[5] = (L[j][5].re * L[j][5].re + L[j][5].im * L[j][5].im) * D[5];
                    case 5 :
                        array_re[4] = (L[j][4].re * L[j][4].re + L[j][4].im * L[j][4].im) * D[4];
                    case 4 :
                        array_re[3] = (L[j][3].re * L[j][3].re + L[j][3].im * L[j][3].im) * D[3];
                    case 3 :
                        array_re[2] = (L[j][2].re * L[j][2].re + L[j][2].im * L[j][2].im) * D[2];
                    case 2 :
                        array_re[1] = (L[j][1].re * L[j][1].re + L[j][1].im * L[j][1].im) * D[1];
                    default :
                        array_re[0] = (L[j][0].re * L[j][0].re + L[j][0].im * L[j][0].im) * D[0];
                }

//                for (k = 0; k <= j - 1; k++) {
//                    array_re[k] = L[j][k].re * L[j][k].re + L[j][k].im * L[j][k].im;
//                }
                sum:
                for(m = 0; m < M-1; m++) {
                    sum_re += array_re[m];
                }
                D[j] = matrix[j][j].re - sum_re;
                D_inv[j] = 1/D[j];
            } else if(j < i){
                switch (j) {
                    case 7 :
                        array_re[6] = (L[i][6].re * L[j][6].re + L[i][6].im * L[j][6].im) * D[6];
                        array_im[6] = (L[i][6].im * L[j][6].re - L[i][6].re * L[j][6].im) * D[6];
                    case 6 :
                        array_re[5] = (L[i][5].re * L[j][5].re + L[i][5].im * L[j][5].im) * D[5];
                        array_im[5] = (L[i][5].im * L[j][5].re - L[i][5].re * L[j][5].im) * D[5];
                    case 5 :
                        array_re[4] = (L[i][4].re * L[j][4].re + L[i][4].im * L[j][4].im) * D[4];
                        array_im[4] = (L[i][4].im * L[j][4].re - L[i][4].re * L[j][4].im) * D[4];
                    case 4 :
                        array_re[3] = (L[i][3].re * L[j][3].re + L[i][3].im * L[j][3].im) * D[3];
                        array_im[3] = (L[i][3].im * L[j][3].re - L[i][3].re * L[j][3].im) * D[3];
                    case 3 :
                        array_re[2] = (L[i][2].re * L[j][2].re + L[i][2].im * L[j][2].im) * D[2];
                        array_im[2] = (L[i][2].im * L[j][2].re - L[i][2].re * L[j][2].im) * D[2];
                    case 2 :
                        array_re[1] = (L[i][1].re * L[j][1].re + L[i][1].im * L[j][1].im) * D[1];
                        array_im[1] = (L[i][1].im * L[j][1].re - L[i][1].re * L[j][1].im) * D[1];
                    case 1 :
                        array_re[0] = (L[i][0].re * L[j][0].re + L[i][0].im * L[j][0].im) * D[0];
                        array_im[0] = (L[i][0].im * L[j][0].re - L[i][0].re * L[j][0].im) * D[0];
                }


//                for (k = 0; k <= j - 1; k++) {
//                    array_re[k] = L[i][k].re * L[j][k].re + L[i][k].im * L[j][k].im;
//                    array_im[k] = L[i][k].im * L[j][k].re - L[i][k].re * L[j][k].im;
//                }
                sum2:
                for(m = 0; m < M; m++) {
                    sum_re += array_re[m];
                    sum_im += array_im[m];
                }
//                temp.re = matrix[i][j].re - sum_re;
//                temp.im = matrix[i][j].im - sum_im;
//                complex_AdivideB( temp, L[j][j], &L[i][j]);
                re = matrix[i][j].re - sum_re;
                im = matrix[i][j].im - sum_im;
//                no complex divide
                L[i][j].re = re * D_inv[j];
                L[i][j].im = im * D_inv[j];
            }
        }
    }


    copy_Yrow:
    for (i = 0; i < M; i++) {
        copy_Ycol:
        for (j = 0; j < T; j++) {
#pragma HLS LOOP_FLATTEN
#pragma HLS PIPELINE
            Y_copy11[i][j] = Y_copy1[i][j];
            L1[i][j].re = L[i][j].re;
            L1[i][j].im = L[i][j].im;
        }
    }
    copy_as:
    for(k = M; k < N+1; k++) {
#pragma HLS PIPELINE
        alpha_sigma2_copy11[k] = alpha_sigma2_copy1[k];
    }
}



/*
 * File trailer for choleskyInverse.c
 *
 * [EOF]
 */
