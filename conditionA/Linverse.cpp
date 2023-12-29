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
void Linverse( complex_fixed L[M][M], real_fixed D_inv[M], complex_fixed L_inv[M][M], real_fixed D_inv1[M], real_fixed alpha_sigma2_copy11[N+1], real_fixed alpha_sigma2_copy12[N+1], complex_fixed Y_copy11[M][T], complex_fixed Y_copy12[M][T])
{

    real_fixed array_re[M-1];
    real_fixed array_im[M-1];
    real_fixed ar;
    real_fixed bi;
    real_fixed bim;
    real_fixed im;
    real_fixed re;
    real_fixed sgnbi;
    real_fixed sum_im;
    real_fixed const1 = 1;
    int i;
    int i2;
    int j;
    int k;
    real_fixed a,b,c,d;
    real_fixed x,y,z;

    /*  L inverse */
    /*      L_inv = forwardSubstitutionInverse(L); */
    memset(&L_inv[0][0], 0, M*M * sizeof(complex_fixed));

    Linverse_diag:
    for (i = 0; i < M; i++) {
#pragma HLS PIPELINE
        D_inv1[i] = D_inv[i];
        L_inv[i][i].re = 1;
        L_inv[i][i].im = 0;
        alpha_sigma2_copy12[i] = alpha_sigma2_copy11[i];
    }
    Linverse_col:
    for (i = 0; i < M-1; i++) {
        Linverse_row:
        for (j = 1; j < M; j++) {
#pragma HLS LOOP_FLATTEN
#pragma HLS PIPELINE
//            if (j > i) {
//                ar = 0.0;
//                re = 0.0;
//                memset(&array_re[0], 0, M * sizeof(double));
//                memset(&array_im[0], 0, M * sizeof(double));

//                switch (j) {
//                    case 7 :
//                        if (i >= 7) break;
//                        array_re[6] = L[j][6].im * L_inv[6][i].im - L[j][6].re * L_inv[6][i].re;
//                        array_im[6] = L[j][6].re * L_inv[6][i].im + L[j][6].im * L_inv[6][i].re;
//                    case 6 :
//                        if (i >= 6) break;
//                        array_re[5] = L[j][5].im * L_inv[5][i].im - L[j][5].re * L_inv[5][i].re;
//                        array_im[5] = L[j][5].re * L_inv[5][i].im + L[j][5].im * L_inv[5][i].re;
//                    case 5 :
//                        if (i >= 5) break;
//                        array_re[4] = L[j][4].im * L_inv[4][i].im - L[j][4].re * L_inv[4][i].re;
//                        array_im[4] = L[j][4].re * L_inv[4][i].im + L[j][4].im * L_inv[4][i].re;
//                    case 4 :
//                        if (i >= 4) break;
//                        array_re[3] = L[j][3].im * L_inv[3][i].im - L[j][3].re * L_inv[3][i].re;
//                        array_im[3] = L[j][3].re * L_inv[3][i].im + L[j][3].im * L_inv[3][i].re;
//                    case 3 :
//                        if (i >= 3) break;
//                        array_re[2] = L[j][2].im * L_inv[2][i].im - L[j][2].re * L_inv[2][i].re;
//                        array_im[2] = L[j][2].re * L_inv[2][i].im + L[j][2].im * L_inv[2][i].re;
//                    case 2 :
//                        if (i >= 2) break;
//                        array_re[1] = L[j][1].im * L_inv[1][i].im - L[j][1].re * L_inv[1][i].re;
//                        array_im[1] = L[j][1].re * L_inv[1][i].im + L[j][1].im * L_inv[1][i].re;
//                    default :
//                        if (i >= 1) break;
//                        array_re[0] = L[j][0].im * L_inv[0][i].im - L[j][0].re * L_inv[0][i].re;
//                        array_im[0] = L[j][0].re * L_inv[0][i].im + L[j][0].im * L_inv[0][i].re;
//                }
//                Linverse_product:
//                for (i2 = i; i2 < j; i2++) {
//#pragma HLS LOOP_TRIPCOUNT max=7
//                    ar += L[j][i2].im * L_inv[i2][i].im - L[j][i2].re * L_inv[i2][i].re;
//                    re -= L[j][i2].re * L_inv[i2][i].im + L[j][i2].im * L_inv[i2][i].re;
//                }
                Linverse_product:
                for (i2 = 0; i2 < M; i2++) {
//                	#pragma HLS LOOP_FLATTEN
//                	#pragma HLS PIPELINE
                	if (i2 >= i && i2 < j){
                    	a = L[j][i2].re;
                    	b = L[j][i2].im;
                    	c = L_inv[i2][i].re;
                    	d = L_inv[i2][i].im;

                    	x = (a + b) * d;
                    	y = (b - a) * c;
                    	z = (c + d) * a;
        //
                    	L_inv[j][i].re += x - z;
                    	L_inv[j][i].im -= y + z;
//                		L_inv[j][i].re += L[j][i2].im * L_inv[i2][i].im - L[j][i2].re * L_inv[i2][i].re;
//                		L_inv[j][i].im -= L[j][i2].re * L_inv[i2][i].im + L[j][i2].im * L_inv[i2][i].re;
                	}
                }
//                sum2:
//                for(i2 = 0; i2 < M; i2++) {
//                    ar += array_re[i2];
//                    re -= array_im[i2];
//                }

//                L_inv[j][i].re = ar;
//                L_inv[j][i].im = re;

            }
//        }
    }

    copy_Yrow:
    for (i = 0; i < M; i++) {
        copy_Ycol:
        for (j = 0; j < T; j++) {
#pragma HLS LOOP_FLATTEN
#pragma HLS PIPELINE
            Y_copy12[i][j] = Y_copy11[i][j];
        }
    }
    Linverse_copy_as:
    for (k = M; k < N+1; k++) {
#pragma HLS PIPELINE
        alpha_sigma2_copy12[k] = alpha_sigma2_copy11[k];
    }
}


/*
 * File trailer for choleskyInverse.c
 *
 * [EOF]
 */
