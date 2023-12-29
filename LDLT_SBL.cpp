//
// Created by imvia on 2023/4/26.
//
//
#include <iostream>
#include <cmath>
#include "Types.h"

void LDLT_SBL(complex_fixed Y[M][T], real_fixed alpha_sigma2[N+1], real_fixed alpha_sigma2_new[N+1])
{
#pragma HLS DATAFLOW
    //  ************  Parameters:   ************
    //    alpha = alpha_sigma2[0:N-1];
    //    sigma2 = alpha_sigma2[N];
    real_fixed diag_Sigma[N];
    real_fixed D_inv1[M];
    real_fixed D_inv[M];
    complex_fixed b_Phi[N][M];
    complex_fixed mu[N][T];
    complex_fixed C[M][M];
    complex_fixed C_inv[M][M];
    complex_fixed L[M][M];
    complex_fixed L_inv[M][M];
    //  ************  Parameters Copy:   ************
    real_fixed alpha_sigma2_copy1[N+1], alpha_sigma2_copy2[N+1], alpha_sigma2_copy3[N+1], alpha_sigma2_copy4[N+1], alpha_sigma2_copy5[N+1], alpha_sigma2_copy11[N+1], alpha_sigma2_copy12[N+1];
    real_fixed diag_Sigma_copy1[N];
    complex_fixed b_Phi_copy1[M][N];
    complex_fixed Y_copy1[M][T], Y_copy2[M][T], Y_copy3[M][T], Y_copy4[M][T], Y_copy5[M][T], Y_copy11[M][T], Y_copy12[M][T];

    /*  Learning loop  */
    /*  Golden reference */

    //    START PROCESSING
    //    ************************************************************************

    /*      C = sigma2 * eye(M) + Phi * diag(alpha) * Phi';  */
    cal_C(alpha_sigma2, C, alpha_sigma2_copy1, Y, Y_copy1);

    int i, j;
//    for(i = 0; i < M; i++)
//        for(j = 0; j < M; j++)
//            std::cout << C[i][j].re << " ";
////
//    std::cout <<'\n';
//    std::cout <<'\n';
//    std::cout <<'\n';

    //    ************************************************************************

    /*      C_inv = cholesky inverse(C); */
//    choleskyinverse(C, C_inv, alpha_sigma2_copy1, alpha_sigma2_copy2, Y_copy1, Y_copy2);

//    cholesky( C, L, alpha_sigma2_copy1, alpha_sigma2_copy11, Y_copy1, Y_copy11);
//
//
//    Linverse( L, L_inv, alpha_sigma2_copy11, alpha_sigma2_copy12, Y_copy11, Y_copy12);
//
//
//
//    Linv_mul( L_inv, C_inv, alpha_sigma2_copy12, alpha_sigma2_copy2, Y_copy12, Y_copy2);
    LDLT( C, L, D_inv, alpha_sigma2_copy1, alpha_sigma2_copy11, Y_copy1, Y_copy11);

//    for(i = 0; i < M; i++)
//        for(j = 0; j < M; j++)
//            std::cout << L[i][j].re << " ";
////
//    std::cout <<'\n';
//    std::cout <<'\n';
//    std::cout << D[0];

    Linverse( L, D_inv, L_inv, D_inv1, alpha_sigma2_copy11, alpha_sigma2_copy12, Y_copy11, Y_copy12);

//    for(i = 0; i < M; i++)
//        for(j = 0; j < M; j++)
//            std::cout << L_inv[i][j].re << " ";
//        std::cout <<'\n';
//    std::cout <<'\n';
//
//    for(i = 0; i < M; i++)
//            std::cout << D_inv[i] << " ";
////
//    std::cout <<'\n';
//    std::cout <<'\n';

    Linv_mul( L_inv, D_inv1, C_inv, alpha_sigma2_copy12, alpha_sigma2_copy2, Y_copy12, Y_copy2);

//    for(i = 0; i < M; i++)
//        for(j = 0; j < M; j++)
//            std::cout << C_inv[i][j].re << " ";
////
//    std::cout <<'\n';
//    std::cout <<'\n';
//
//

    //    ************************************************************************

    /*      b_Phi = Phi' * Cinv */
    PhiCinv(C_inv, b_Phi, alpha_sigma2_copy2, alpha_sigma2_copy3, Y_copy2, Y_copy3);
//    std::cout <<'\n';
//    std::cout <<'\n';
//    for(i = 0; i < N; i++)
//        for(j = 0; j < M; j++)
//            std::cout << b_Phi[i][j].im << " ";
    //    ************************************************************************

    /*      Sigma = diag(alpha) - diag(alpha) * Phi' * Cinv * Phi * diag(alpha) */
    cal_Sigma(b_Phi, alpha_sigma2_copy3, diag_Sigma, b_Phi_copy1, alpha_sigma2_copy4, Y_copy3, Y_copy4);
//    for(j = 0; j < N; j++)
//        std::cout << diag_Sigma[j] << " ";

    //    ************************************************************************

    /*  =========== update alpha ===========   */
    /*      mu = diag(alpha) * Phi' * C_inv * Y */
    /*      musq = mean(abs(mu).^2,2); */
    /*      alpha = musq + real(diag(Sigma)); */
    cal_alpha(b_Phi_copy1, diag_Sigma, Y_copy4, alpha_sigma2_copy4,  mu, diag_Sigma_copy1, Y_copy5, alpha_sigma2_copy5);

    //    ************************************************************************

    /*  =========== update sigma2 ===========   */
    /*      resid = Y - Phi * mu;     */
    /*      gamma1 = 1 - real(diag(Sigma)) ./ (alpha + eps); */
    /*      sigma2 = (norm(resid, 'fro')^2 + sigma2 * sum(gamma1) ) / ( M );  */
    cal_sigma2(alpha_sigma2_copy5, Y_copy5, diag_Sigma_copy1, mu, alpha_sigma2_new);


//        for(j = 0; j < N+1; j++)
//            std::cout << alpha_sigma2_new[j] << " ";
//    std::cout <<'\n';
//    std::cout <<'\n';
    //    ************************************************************************

}

