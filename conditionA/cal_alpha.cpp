//
// Created by imvia on 2023/5/24.
//

#include "Types.h"

void cal_alpha(complex_fixed b_Phi[M][N], real_fixed diag_Sigma[N], complex_fixed Y[M][T], real_fixed alpha_sigma2_copy4[N+1], complex_fixed mu[N][T], real_fixed diag_Sigma_copy1[N], complex_fixed Y_copy5[M][T], real_fixed alpha_sigma2_copy5[N+1])
{
	#pragma HLS ARRAY_PARTITION variable=b_Phi type=complete dim=1
    real_fixed Phi_re;
    real_fixed absxk;
    real_fixed tmp;
    real_fixed re;
    real_fixed im;
    real_fixed musq[N];
    int j;
    int i;
    int k;
    real_fixed a,b,c,d;
    real_fixed x,y,z;


    /*      mu = diag(alpha) * Phi' * C_inv * Y */
    /*      musq = mean(abs(mu).^2,2); */
    /*      alpha = musq + real(diag(Sigma)); */
//    memset(&mu[0][0], 0, N*T * sizeof(complex_fixed));

    memset(&musq[0], 0, N * sizeof(real_fixed));
//    add delay
    memset(&mu[0][0], 0, N*5 * sizeof(complex_fixed));


    row:
    for (i = 0; i < N; i++) {
        col:
        for (j = 0; j < T; j++) {
#pragma HLS LOOP_FLATTEN
#pragma HLS PIPELINE
            re = 0;
            im = 0;
            product:
            for (k = 0; k < M; k++) {
#pragma HLS UNROLL
            	a = b_Phi[k][i].re;
            	b = b_Phi[k][i].im;
            	c = Y[k][j].re;
            	d = Y[k][j].im;

            	x = (a + b) * d;
            	y = (b - a) * c;
            	z = (c + d) * a;
//
				re += z - x;
				im += y + z;
//                re += b_Phi[k][i].re * Y[k][j].re - b_Phi[k][i].im * Y[k][j].im;
//                im += b_Phi[k][i].re * Y[k][j].im + b_Phi[k][i].im * Y[k][j].re;
            }
            mu[i][j].re = alpha_sigma2_copy4[i] * re;
            mu[i][j].im = alpha_sigma2_copy4[i] * im;
            musq[i] += mu[i][j].re * mu[i][j].re + mu[i][j].im * mu[i][j].im;
        }
    }

    copy_Yrow:
    for (i = 0; i < M; i++) {
        copy_Ycol:
        for (j = 0; j < T; j++) {
#pragma HLS LOOP_FLATTEN
#pragma HLS PIPELINE
            Y_copy5[i][j] = Y[i][j];
        }
    }


    copy_as:
    for (i = 0; i < N; i++) {
#pragma HLS PIPELINE
        alpha_sigma2_copy5[i] = (musq[i]>>3) + diag_Sigma[i];
        diag_Sigma_copy1[i] = diag_Sigma[i];
    }

    alpha_sigma2_copy5[N] = alpha_sigma2_copy4[N];

}
