#include <iostream>
#include <cmath>
#include <cstring>
#include "Types.h"
#include<iomanip>

//#include <E:/ImVIA/FPGA/HLS/HLS_arbitrary_Precision_Types-master/include/ap_fixed.h>
//#include "E:/ImVIA/FPGA/HLS/HLS_arbitrary_Precision_Types-master/include/ap_int.h"
//#include "ap_fixed.h"
//#include "ap_int.h"


int main()
{

    const real_fixed PI = 3.1415926535897932384626433832795;
    int rx = 2;
//    int N = (int)(180 - rx) / rx + 1;
    real_fixed grid[N];
    int i;
    int j;
    for (i = 0; i < N; i++) {
        grid[i] = i * rx;
    };

//    real_fixed theta[] = {60, 88};
//    real_fixed amp[] = {1, 1};
//    int M = 8;
//    real_fixed f = 300;
//    real_fixed c = 1500;
//    real_fixed lamda = c / f;
//    real_fixed pos[M];
//    for (i = 0; i < M; i++) {
//        pos[i] = (i - (M - 1) / 2.0) * lamda / 2.0;
//    };

    /* Input Signal: Y */
//    complex_fixed Y1[M] = {{
//                            0.242062785739483, /* re */
//                            -1.64258702333431  /* im */
//                    },
//                    {
//                            1.8530859463325511, /* re */
//                            -1.244276681746201  /* im */
//                    },
//                    {
//                            1.2819853992513011, /* re */
//                            0.575893648299133   /* im */
//                    },
//                    {
//                            0.064182228782395, /* re */
//                            0.000709539704952  /* im */
//                    },
//                    {
//                            0.556514982050625, /* re */
//                            -1.174466584241536 /* im */
//                    },
//                    {
//                            1.874618905966339, /* re */
//                            -0.601904898143015 /* im */
//                    },
//                    {
//                            1.5530193004753841, /* re */
//                            0.8109754457841     /* im */
//                    },
//                    {
//                            0.148180519193266, /* re */
//                            0.396541881226616  /* im */
//                    }};
    complex_fixed Y1[M] =   {
        {
            0.695494227842183, -0.0920048368751170
        },
        {
            -0.391653354138112, -0.505161695163649
        },
        {
            0.00451381212092238, -1.63083253012561
        },
        {
            1.34412157104494, -1.30084692085581
        },
        {
            1.11606334080248, +0.145750682754139
        },
        {
            -0.165807244341588, -0.137307327912181
        },
        {
            0.333289637362201, -1.48088421452496
        },
        {
            1.76544906247274, -0.949902219070189
        }
    };
    complex_fixed Y[M][T];
    for (i = 0; i < M; i++) {
        for (j = 0; j < T; j++){
            Y[i][j] = Y1[i];
        }
    }
    /* Dictionary matrix: Phi */
//    complex_fixed Phi[720];
//    memset(&Phi[0], 0, 720U * sizeof(complex_fixed));
    int m;
    int n;
    real_fixed V;
    real_fixed pift;
//    for (m = 0; m < M; m++) {
//        for (n = 0; n < N; n++) {
//            V = cos(grid[n] / 180 * PI);
//            pift = 2 * PI * f * pos[m] * V / c;
//            Phi[m+(n<<3)].re = cos(pift);
//            Phi[m+(n<<3)].im = sin(pift);
//        }
//    }
    /* Output : Power spectrum: alpha*/
    int Y_win = 8;
    real_fixed alpha_sigma2_MM[N+1][Y_win];

    real_fixed alpha_sigma2[N+1];
    real_fixed alpha_sigma2_new[N+1];
    for (i = 0; i < 90; i++) {
        for (n = 0; n < Y_win; n++) {
            alpha_sigma2_MM[i][n] = 1.0;
        }
    }
    for (n = 0; n < Y_win; n++){
        alpha_sigma2_MM[N][n] = 0.001;
    }


    int t;

    int counter;
    int iter;

    complex_fixed data[M][T*Y_win];
    for (n = 0; n < Y_win; n++){
        for (i = 0; i < M; i++) {
            for (j = 0; j < T; j++){
                data[i][n*Y_win+j] = Y[i][j];
            }
        }
    }

    double fenzi = 0;
    double fenmu = 0;
    double epsilon[M];


//    printf("******************Begin Test!******************\n");
    std::cout << "******************Begin Test!******************\n" << std::endl;
    for(iter = 0; iter < 3; iter++) {
        for (n = 0; n < Y_win; n++){

//            if ( epsilon[n] < 0.001 ) continue;

            //    	Y = data;
            for (i = 0; i < M; i++) {
                for (j = 0; j < T; j++){
                    Y[i][j] = data[i][n*Y_win+j];
                }
            }
            for (m = 0; m < N+1; m++){
                alpha_sigma2[m] = alpha_sigma2_MM[m][n];
                alpha_sigma2_new[m] = alpha_sigma2_MM[m][n];
            }


            LDLT_SBL(Y, alpha_sigma2, alpha_sigma2_new);

            fenzi = 0;
            fenmu = 0;

            for (m = 0; m < N; m++){
                fenzi += fabs( double(alpha_sigma2_new[m] - alpha_sigma2_MM[m][n] ));
                fenmu += fabs( double(alpha_sigma2_MM[m][n] ));
            }
            epsilon[n] = fenzi/fenmu;



            for (m = 0; m < N+1; m++){
                alpha_sigma2_MM[m][n] = alpha_sigma2_new[m];
            }

//            if (iter == 136){
//                for(i = 0; i < 90; i++)
//                    std::cout << std::setprecision(16)<< alpha_sigma2_new[i] << " ";
//                std::cout <<'\n';
//                std::cout <<'\n';
//            }


        }

    }


    for(i = 0; i < 91; i++)
        std::cout << std::setprecision(16) << alpha_sigma2_new[i] << " ";
    std::cout <<'\n';
    std::cout << "******************Test Pass!******************\n" << std::endl;

//    ap_uint<3> a = 6;
//    std::cout << a <<'\n';
//
//    ap_fixed<32,5> b = 1.3354;
//    std::cout << b;
    return 0;
}

