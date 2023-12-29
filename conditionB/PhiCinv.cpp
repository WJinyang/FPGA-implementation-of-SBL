//
// Created by imvia on 2023/5/24.
//

#include "Types.h"

void PhiCinv(complex_fixed C_inv[M][M], complex_fixed b_Phi[N][M], real_fixed alpha_sigma2_copy2[N+1], real_fixed alpha_sigma2_copy3[N+1], complex_fixed Y_copy2[M][T], complex_fixed Y_copy3[M][T])
{

    int i;
    int j;
    int k;
    int m;
    int n;
    real_fixed a,b,c,d;
    real_fixed x,y,z;

    copy_as:
    for (k = 0; k < N+1; k++) {
#pragma HLS PIPELINE
        alpha_sigma2_copy3[k] = alpha_sigma2_copy2[k];
    }
    copy_Yrow:
    for (i = 0; i < M; i++) {
        copy_Ycol:
        for (j = 0; j < T; j++) {
#pragma HLS LOOP_FLATTEN
#pragma HLS PIPELINE
            Y_copy3[i][j] = Y_copy2[i][j];
        }
    }

//    memset(&b_Phi[0][0], 0, 500 * sizeof(complex_fixed));
    /*      b_Phi = Phi' * Cinv */
//    for (m = 0; m < (M>>2); m++) {
        row:
        for (i = 0; i < N; i++) {
            col:
            for (j = 0; j < M; j++) {
#pragma HLS LOOP_FLATTEN
#pragma HLS PIPELINE
            	b_Phi[i][j].re = 0;
            	b_Phi[i][j].im = 0;
                product:
                for (k = 0; k < M; k++) {
                	a = Phi[k][i].re;
                	b = -Phi[k][i].im;
                	c = C_inv[k][j].re;
                	d = C_inv[k][j].im;

                	x = (a + b) * d;
                	y = (b - a) * c;
                	z = (c + d) * a;
    //
                	b_Phi[i][j].re += z - x;
                	b_Phi[i][j].im += y + z;

//                    b_Phi[i][j].re += Phi[k][i].re * C_inv[k][j].re + Phi[k][i].im * C_inv[k][j].im;
//                    b_Phi[i][j].im += Phi[k][i].re * C_inv[k][j].im - Phi[k][i].im * C_inv[k][j].re;
                }
            }
        }
//    }
//        for (k = 0; k < 500; k++) {
//    #pragma HLS PIPELINE
//            n = k;
//        }

}
