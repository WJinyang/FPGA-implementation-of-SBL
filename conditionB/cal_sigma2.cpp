//
// Created by imvia on 2023/5/24.
//

#include "Types.h"

void cal_sigma2(real_fixed alpha_sigma2_copy5[N+1], complex_fixed Y[M][T], real_fixed diag_Sigma[N], complex_fixed mu[N][T], real_fixed alpha_sigma2_new[N+1])
{
#pragma HLS ARRAY_PARTITION variable=mu type=complete dim=2
    complex_fixed resid[M][T];
#pragma HLS ARRAY_PARTITION variable=resid type=complete dim=2
    real_fixed gamma1[N];
    real_fixed Phi_im;
    real_fixed Phi_re;
    real_fixed absxk;
    real_fixed tmp;
    real_fixed scale;
    real_fixed t;
    real_fixed const1 = 1;
    int b_i;
    int i1;
    int i;
    int j;
    int k;
    real_fixed a,b,c,d;
    real_fixed x,y,z;

    /*      resid = Y - Phi * mu;     */
    memset(&resid[0][0], 0, M*T * sizeof(complex_fixed));
    //    add delay
//    memset(&gamma1[0], 0, N * sizeof(real_fixed));

        product:
        for(k = 0; k < N; k++) {
            row:
            for(i = 0; i < M; i++) {
#pragma HLS LOOP_FLATTEN
#pragma HLS PIPELINE
            	a = Phi[i][k].re;
            	b = Phi[i][k].im;
                col:
                for(j = 0; j < T; j++) {

                	c = mu[k][j].re;
                	d = mu[k][j].im;

                	x = (a + b) * d;
                	y = (b - a) * c;
                	z = (c + d) * a;
    //
                	resid[i][j].re += (z - x);
                	resid[i][j].im += (y + z);
                }
            }
        }

    Phi_re = 0.0;
    norm_row:
    for (b_i = 0; b_i < M; b_i++) {
        norm_col:
        for (i = 0; i < T; i++) {
#pragma HLS LOOP_FLATTEN
#pragma HLS PIPELINE
            resid[b_i][i].re = Y[b_i][i].re - resid[b_i][i].re;
            resid[b_i][i].im = Y[b_i][i].im - resid[b_i][i].im;
            Phi_re += resid[b_i][i].re * resid[b_i][i].re + resid[b_i][i].im * resid[b_i][i].im;
        }
    }

/*      gamma1 = 1 - real(diag(Sigma)) ./ (alpha + eps); */
    Phi_im = 0.0;
    gamma:
    for (i = 0; i < N; i++) {
#pragma HLS LOOP_FLATTEN
#pragma HLS PIPELINE
        gamma1[i] = const1 - diag_Sigma[i] / alpha_sigma2_copy5[i];
        Phi_im += gamma1[i];
        alpha_sigma2_new[i] = alpha_sigma2_copy5[i];

    }
/*      sigma2 = (norm(resid, 'fro')^2 + sigma2 * sum(gamma1) ) / ( M );  */
    alpha_sigma2_new[N] = ((Phi_re>>3) + alpha_sigma2_copy5[N] * Phi_im)>>4;


}
