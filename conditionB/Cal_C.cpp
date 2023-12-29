//
// Created by imvia on 2023/5/26.
//

//
// Created by imvia on 2023/5/24.
//

#include "Types.h"
#include <iostream>

void cal_C(real_fixed alpha_sigma2[N+1], complex_fixed C[M][M], real_fixed alpha_sigma2_copy1[N+1], complex_fixed Y[M][T], complex_fixed Y_copy1[M][T])
{
    int i;
    int j;
    int k;
    int m;
    int n;
//    real_fixed re;
//    real_fixed im;
    real_fixed a,b,c,d;
    real_fixed x,y,z;

    /*      C = sigma2 * eye(M) + Phi * diag(alpha) * Phi';  */
    // Initialize Phi, alpha, and C with zeros
    memset(&C[0][0], 0, M*M * sizeof(complex_fixed));

    /* compute (C)i,j */
//    for (m = 0; m < (M>>2); m++) {
        product:
        for(k = 0; k < N; k++) {
            row:
            for(i = 0; i < M; i++) {
#pragma HLS LOOP_FLATTEN
#pragma HLS PIPELINE
            	a = Phi[i][k].re;
            	b = Phi[i][k].im;
                col:
                for(j = 0; j < M; j++) {
//                    j = (m<<2) + n;
//                    j = n;

                	c = Phi[j][k].re;
                	d = -Phi[j][k].im;

                	x = (a + b) * d;
                	y = (b - a) * c;
                	z = (c + d) * a;
    //
                	C[i][j].re += (z - x) * alpha_sigma2[k];
                	C[i][j].im += (y + z) * alpha_sigma2[k];
//                    C[i][j].re += ( Phi[i][k].re * Phi[j][k].re + Phi[i][k].im * Phi[j][k].im ) * alpha_sigma2[k] ;
//                    C[i][j].im += ( Phi[i][k].im * Phi[j][k].re - Phi[i][k].re * Phi[j][k].im ) * alpha_sigma2[k] ;
                }
            }
        }
//    }
    addsigma2:
    for (i = 0; i < M; i++) {
#pragma HLS PIPELINE
        C[i][i].re += alpha_sigma2[N];
        alpha_sigma2_copy1[i] = alpha_sigma2[i];
    }
    copy_Yrow:
    for (i = 0; i < M; i++) {
        copy_Ycol:
        for (j = 0; j < T; j++) {
#pragma HLS LOOP_FLATTEN
#pragma HLS PIPELINE
            Y_copy1[i][j] = Y[i][j];
        }
    }
    copy_as:
    for(k = M; k < N+1; k++) {
#pragma HLS PIPELINE
        alpha_sigma2_copy1[k] = alpha_sigma2[k];
    }
}

