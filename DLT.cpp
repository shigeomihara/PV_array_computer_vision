/*********************************************************************
 *      DLT.cpp
 ********************************************************************/
// #include <iostream>
// #include "main.hpp"
// #include "DLT.hpp"
// #include "svd.hpp"
// #include "ludcmp.h"

// Constructor
DLT::DLT(KizyunTen kt)
    :m_kt(kt), m_n(kt.get_n()){
}

/** Calculate camera matrix P ***/
void DLT::calcP(VecDoub_O &P){
    double x[m_n][3];
    double X[m_n][4];
    for(int i=0; i<m_n; i++){
        for(int j=0; j<2; j++){
            x[i][j] = m_kt.m_xTilder[i][j];
        }
        x[i][2] = 1.0;
    
        for(int j=0; j<3; j++){
            X[i][j] = m_kt.m_Xtilder[i][j];
        }
        X[i][3] = 1.0;
    }
  
    double A[2*m_n][12];
    for(int i=0; i<2*m_n; i++){
        int ix = i/2;
        for(int j=0; j<12; j++){
            if(i%2==0){  // i is even
                if(j<4) A[i][j] = 0.0;
                else if(j<8) A[i][j] = -x[ix][2]*X[ix][j-4];
                else A[i][j] = x[ix][1]*X[ix][j-8];
            }else{   // i is odd
                if(j<4) A[i][j] = x[ix][2]*X[ix][j];
                else if(j<8) A[i][j] = 0.0;
                else A[i][j] = -x[ix][0]*X[ix][j-8];
            }
        }
    }

    // std::cout << "matrix A\n";//////////////////////////t
    // for(int i=1; i<=2*m_n; i++){
    //     for(int j=1; j<=12; j++){
    //         printf("%10f ", A[i][j]);
    //     }
    //     printf("\n");
    // }

    MatDoub a(2*m_n, 12);
    for(int i=0; i<2*m_n; i++){
        for(int j=0; j<12; j++){
            a[i][j] = A[i][j];
        }
    }

    SVD svd(a);
    
    std::cout << "SVD diag(D)\n";////////////////////////t
    for(int i=0; i<svd.w.size(); i++){
        printf("%10f ", svd.w[i]);
    }
    printf("\n");

    // double P[13];
    for(int i=0; i<12; i++) P[i] = svd.v[i][11];
    
    // std::cout << "hello DLT calc\n";
    // std::cout << "n=" << m_kt.get_n() << std::endl;

    printf("Normalized P\n");///////////////////
    for(int i=0; i<12; i++){//////////////////
        printf("%10f ", P[i]);
        if(i%4==3) printf("\n");
    }
}

/***     *****/
void DLT::denormalize(VecDoub_IO &P, VecDoub_O &Ctilder, MatDoub_O &Minv){
    // Denormalization
    double matP[3][4];
    for(int i=0; i<3; i++){
        for(int j=0; j<4; j++){
            matP[i][j] = P[j+ i*4];
        }
    }

    // std::cout << "matP\n";
    // for(int i=0; i<3; i++){//////////////////t
    //     for(int j=0; j<4; j++){
    //         printf("%15f ",matP[i][j]); 
    //         // std::cout << matP[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
	
    double TinvP[3][4];
    for(int i=0; i<3; i++){
        for(int j=0; j<4; j++){
            TinvP[i][j] = 0.0;
            for(int k=0; k<3; k++){
                TinvP[i][j] += m_kt.m_Tinv[i][k]*matP[k][j];
            }
        }
    }

    double TinvPU[3][4];
    for(int i=0; i<3; i++){
        for(int j=0; j<4; j++){
            TinvPU[i][j] = 0.0;
            for(int k=0; k<4; k++){
                TinvPU[i][j] += TinvP[i][k]*m_kt.m_U[k][j];
            }
        }
    }

    // denormalize P
    for(int i=0; i<3; i++){
        for(int j=0; j<4; j++){
            P[4*i +j] = TinvPU[i][j];
        }
    }
    
    std::cout << "DLT::denormalize(): 'Denormalized P'\n";////////////////////t
    for(int i=0; i<3; i++){//////////////////
        for(int j=0; j<4; j++){
            // printf("%12f ", TinvPU[i][j]/TinvPU[3][4]);
            printf("%12g ", P[i*4+j]/P[11]);
            // std::cout << TinvPU[i][j] << " ";
        }
        std::cout << std::endl;
    }
    
    MatDoub M(3, 3);
    for(int i=0; i<=2; i++){
        for(int j=0; j<=2; j++){
            M[i][j] = TinvPU[i][j];
        }
    }
  
    LUdcmp lud(M);
  
    lud.inverse(Minv);

    // std::cout << "Minv\n";
    // for(int i=0; i<=2; i++){//////////////////t
    //     for(int j=0; j<=2; j++){
    //         printf("%10f ", Minv[i][j]);
    //         // std::cout << Minv[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
  
    // std::cout << "M*Minv\n";
    // for(int i=0; i<=2; i++){//////////////////t
    //     for(int j=0; j<=2; j++){
    //         double x {0.0};
    //         for(int k=0; k<3; k++){
    //             x += M[i][k]*Minv[k][j];
    //         }
    //         printf("%10f ", x);
    //         // std::cout << x << " ";
    //     }
    //     std::cout << std::endl;
    // }

    for(int i=0; i<3; i++){
        Ctilder[i] = 0.0;
        for(int k=0; k<3; k++){
            Ctilder[i] += -Minv[i][k]*TinvPU[k][3];
        }
    }

    std::cout << "DLT::denormalize(): Ctilder\n";
    for(int i=0; i<3; i++){//////////////////t
    	std::cout << "DLT::denormalize(): Ctilder[" << i << "]= " << Ctilder[i] << std::endl;
    }
}

/***** Linear triangulation **************************
 * inputs: P[0..11], Pdash[0..11], m[0,1], mDash[0,1] 
 *         m[2]=1, mDash[2]=1  not given 
 * output: M[0...3]                       ************/
void DLT::triangulation(VecDoub_I &P, VecDoub_I &Pdash,
                        double m[], double mDash[],
                        VecDoub_O &M){
    MatDoub A(4, 4);
    // printf("m[1]=%g, P[9]=%g\n", m[1], P[9]);///////////t
    // exit(0);
    for(int j=0; j<4; j++){
        A[0][j] = m[0]*P[8+j] - P[j];
        A[1][j] = m[1]*P[8+j] - P[4+j];
        A[2][j] = mDash[0]*Pdash[8+j] - Pdash[j];
        A[3][j] = mDash[1]*Pdash[8+j] - Pdash[4+j];
    }

    SVD svd(A);
    
    std::cout << "SVD diag(D)\n";////////////////////////t
    for(int i=0; i<svd.w.size(); i++){
        printf("%10f ", svd.w[i]);
    }
    printf("\n");

    for(int i=0; i<4; i++) M[i] = svd.v[i][3];

}
