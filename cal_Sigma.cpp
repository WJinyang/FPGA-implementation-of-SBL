//
// Created by imvia on 2023/5/24.
//

#include "Types.h"

void cal_Sigma(complex_fixed b_Phi[N][M], real_fixed alpha_sigma2_copy3[N+1], real_fixed diag_Sigma[N], complex_fixed b_Phi_copy1[M][N], real_fixed alpha_sigma2_copy4[N+1], complex_fixed Y_copy3[M][T], complex_fixed Y_copy4[M][T])
{
    real_fixed absxk;
    real_fixed tmp;
    real_fixed scale;
    real_fixed t;
    int i;
    int j;

    /*      Sigma = diag(alpha) - diag(alpha) * Phi' * Cinv * Phi * diag(alpha) */

    memset(&diag_Sigma[0], 0, N * sizeof(real_fixed));
//    add delay
    memset(&b_Phi_copy1[0][0], 0, N*5 * sizeof(complex_fixed));


    row:
    for (i = 0; i < N; i++) {
        product:
        for (j = 0; j < M; j++) {
#pragma HLS LOOP_FLATTEN
#pragma HLS PIPELINE
            diag_Sigma[i] += b_Phi[i][j].re * Phi[j][i].re - b_Phi[i][j].im * Phi[j][i].im;
            b_Phi_copy1[j][i] = b_Phi[i][j];
        }
    }
    a_a:
    for (i = 0; i < N; i++) {
#pragma HLS PIPELINE
        scale = -alpha_sigma2_copy3[i] * alpha_sigma2_copy3[i];
        diag_Sigma[i] *= scale;
        diag_Sigma[i] += alpha_sigma2_copy3[i];

        alpha_sigma2_copy4[i] = alpha_sigma2_copy3[i];
    }
    copy_Yrow:
    for (i = 0; i < M; i++) {
        copy_Ycol:
        for (j = 0; j < T; j++) {
#pragma HLS LOOP_FLATTEN
#pragma HLS PIPELINE
            Y_copy4[i][j] = Y_copy3[i][j];
        }
    }
    alpha_sigma2_copy4[N] = alpha_sigma2_copy3[N];

}
