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
void Linv_mul( complex_fixed L_inv[M][M], real_fixed D_inv[M], complex_fixed inverse[M][M], real_fixed alpha_sigma2_copy12[N+1], real_fixed alpha_sigma2_copy2[N+1], complex_fixed Y_copy12[M][T], complex_fixed Y_copy2[M][T])
{

    real_fixed bim;
    real_fixed sgnbi;
    real_fixed sum_im;
    real_fixed sum_re;
    int b_i;
    int i;
    int j;
    int i2;
    int k;

    memset(&inverse[0][0], 0, M*M * sizeof(complex_fixed));
    /*  matrix multiply */
    Linv_mul_row:
    for (i = 0; i < M; i++) {
        Linv_mul_col:
        for (j = 0; j < M; j++) {
            Linv_mul_product:
            for (k = 0; k < M; k++) {
#pragma HLS LOOP_FLATTEN
#pragma HLS PIPELINE
                inverse[i][j].re += ( L_inv[k][i].re * L_inv[k][j].re + L_inv[k][i].im * L_inv[k][j].im ) * D_inv[k] ;
                inverse[i][j].im += ( L_inv[k][i].re * L_inv[k][j].im - L_inv[k][i].im * L_inv[k][j].re ) * D_inv[k] ;
            }
        }
    }

    Linv_mul_copy_Yrow:
    for (i = 0; i < M; i++) {
        Linv_mul_copy_Ycol:
        for (j = 0; j < T; j++) {
#pragma HLS LOOP_FLATTEN
#pragma HLS PIPELINE
            Y_copy2[i][j] = Y_copy12[i][j];
        }
    }

    Linv_mul_copy_as:
    for (k = 0; k < N+1; k++) {
#pragma HLS LOOP_FLATTEN
#pragma HLS PIPELINE
        alpha_sigma2_copy2[k] = alpha_sigma2_copy12[k];
    }
}


/*
 * File trailer for choleskyInverse.c
 *
 * [EOF]
 */
