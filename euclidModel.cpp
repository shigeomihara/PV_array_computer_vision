/****************************************************************
 *      euclidModel.cpp
 ***************************************************************/
#include <iostream>
#include "main.hpp"
// #include "nr3.h"
// #include "ludcmp.h"
#include <cmath>

/*** Constructor doing nothing *********************************/
EuclidModel::EuclidModel()
    :m_P(3, 4){
}

/*** Constructor *********************************/
EuclidModel::EuclidModel(KizyunTen kt)
    :m_kizyunTen(kt), m_P(3, 4){
    m_w2 = 0.5*m_kizyunTen.m_width;
    m_h2 = 0.5*m_kizyunTen.m_height;
}

/*****************************************************/
KizyunTen EuclidModel::getKizyunTen(){
    return m_kizyunTen;
}

/******************************************************/
void EuclidModel::printImageError(){
    double error {0.0};

    for(int i=0; i<m_kizyunTen.get_n(); i++){
	VecDoub Xi(4);
	for(int j=0; j<3; j++) Xi[j] = m_kizyunTen.m_XtilderOrig[i][j];
	Xi[3] = 1.0;

	VecDoub xi(3, 0.0);
	for(int j=0; j<3; j++){
	    for(int k=0; k<4; k++) xi[j] += m_P[j][k]*Xi[k];
	}

	for(int j=0; j<=1; j++){
	    double del {xi[j]/xi[2]-m_kizyunTen.m_xTilderOrig[i][j]};
	    error += del*del;
	}

	// printf("ref:x[%d]=%g, cal:x[%d]=%g\n", i, m_kizyunTen.m_xTilderOrig[i][0],
	//        i, xi[0]/xi[2]);
	// printf("ref:y[%d]=%g, cal:y[%d]=%g\n", i, m_kizyunTen.m_xTilderOrig[i][1],
	//        i, xi[1]/xi[2]);
    }
    printf("error= %g\n", error);
}

/****** Return matrix product
 * Assume A.ncols()==B.nrows(), A.nrows()<=C.nrows(), B.ncols()<=C.ncols() *******/
void EuclidModel::prod(MatDoub_I &A, MatDoub_I &B, MatDoub_O &C){
    for(int i=0; i<A.nrows(); i++){
	for(int j=0; j<B.ncols(); j++){
	    C[i][j] = 0.0;
	    for(int k=0; k<A.ncols(); k++){
		C[i][j] += A[i][k]*B[k][j];
	    }
	}
    }
}

/****** Print matrix ********************************/
void EuclidModel::printMatrix(MatDoub_I &A, std::string str){
    std::cout << str << std::endl;
    for(int i=0; i<A.nrows(); i++){
	for(int j=0; j<A.ncols(); j++){
	    printf("[%d][%d]=%10.4f ", i, j, A[i][j]);
	}
	printf("\n");
    }
}
    
/*****************************************************/
void EuclidModel::getJ(VecDoub_I &vecAlpha, MatDoub_O &Jval){
    // printf("EuclidModel::getJ\n");///////////////
    m_alpha = vecAlpha[0];
    m_t[0] = vecAlpha[4];
    m_t[1] = vecAlpha[5];
    m_t[2] = vecAlpha[6];
    // for(int i=1; i<=3; i++) printf("vecAlpha[%d]= %g\n", i, vecAlpha[i]);//////////////
    double cx {std::cos(vecAlpha[1])};
    double sx {std::sin(vecAlpha[1])};
    double cy {std::cos(vecAlpha[2])};
    double sy {std::sin(vecAlpha[2])};
    double cz {std::cos(vecAlpha[3])};
    double sz {std::sin(vecAlpha[3])};

    MatDoub K(3, 3, 0.0);
    MatDoub Qx(3, 3, 0.0);
    MatDoub Qy(3, 3, 0.0);
    MatDoub Qz(3, 3, 0.0);
    MatDoub Qxd(3, 3, 0.0);
    MatDoub Qyd(3, 3, 0.0);
    MatDoub Qzd(3, 3, 0.0);
    // MatDoub R(3, 3);
    MatDoub Rt(3, 4);
    MatDoub P(3, 4);
    MatDoub X(m_kizyunTen.get_n(), 4);   // world homogenious coordinates
                                         // of the reference points
    
    MatDoub A(3, 3);  // temporaly matrix used when computing matrix product
    MatDoub B(3, 3);  // temporaly matrix used when computing matrix product

    K[0][0] = m_alpha;
    K[0][2] = m_w2;
    K[1][1] = m_alpha;
    K[1][2] = m_h2;
    K[2][2] = 1.0;
    
    Qx[0][0] = 1.0;
    Qx[1][1] = cx;
    Qx[1][2] = sx;
    Qx[2][1] = -sx;
    Qx[2][2] = cx;
    
    Qxd[1][1] = -sx;
    Qxd[1][2] = cx;
    Qxd[2][1] = -cx;
    Qxd[2][2] = -sx;
    
    Qy[0][0] = cy;
    Qy[0][2] = -sy;
    Qy[1][1] = 1.0;
    Qy[2][0] = sy;
    Qy[2][2] = cy;
    
    Qyd[0][0] = -sy;
    Qyd[0][2] = -cy;
    Qyd[2][0] = cy;
    Qyd[2][2] = -sy;
    
    Qz[0][0] = cz;
    Qz[0][1] = sz;
    Qz[1][0] = -sz;
    Qz[1][1] = cz;
    Qz[2][2] = 1.0;
    
    Qzd[0][0] = -sz;
    Qzd[0][1] = cz;
    Qzd[1][0] = -cz;
    Qzd[1][1] = -sz;

    // Rt[0][0] = cy*cz;
    // Rt[0][1] = sx*sy*cz+cx*sz;
    // Rt[0][2] = -cx*sy*cz+sx*sz;
    // Rt[1][0] = -cy*sz;
    // Rt[1][1] = -sx*sy*sz+cx*cz;
    // Rt[1][2] = cx*sy*sz+sx*cz;
    // Rt[2][0] = sy;
    // Rt[2][1] = -sx*cy;
    // Rt[2][2] = cx*cy;
    // Change the product order to Qx Qz Qy
    MatDoub matTmp(3, 3);
    prod(Qx, Qz, matTmp);
    prod(matTmp, Qy, Rt);

    for(int i=0; i<3; i++) Rt[i][3] = m_t[i];
    // printMatrix(Rt, "Rt");//////////////////

    prod(K, Rt, P);
    // printMatrix(P, "P");///////////////////

    for(int i=0; i<m_kizyunTen.get_n(); i++){
	for(int j=0; j<3; j++){
	    X[i][j] = m_kizyunTen.m_XtilderOrig[i][j];
	}
	X[i][3] = 1.0;
    }
    // printMatrix(X, "X");//////////////////////

    // // check
    // MatDoub RR(3, 3);
    // prod(Qy, Qx, A);
    // prod(Qz, A, RR);
    // printMatrix(R, "R");
    // printMatrix(RR, "RR");

    MatDoub KQxdQzQy(3, 3);
    prod(Qz, Qy, A);
    prod(Qxd, A, B);
    prod(K, B, KQxdQzQy);
    
    MatDoub KQxQzQyd(3, 3);
    prod(Qz, Qyd, A);
    prod(Qx, A, B);
    prod(K, B, KQxQzQyd);
    
    MatDoub KQxQzdQy(3, 3);
    prod(Qzd, Qy, A);
    prod(Qx, A, B);
    prod(K, B, KQxQzdQy);

    // printMatrix(KQzQyQxd, "KQzQyQxd");///////////////
    // printMatrix(KQzQydQx, "KQzQydQx");///////////////
    // printMatrix(KQzdQyQx, "KQzdQyQx");///////////////
    
    for(int i=0; i<m_kizyunTen.get_n(); i++){
	// double PXi2 {0.0};
	// for(int k=0; k<3; k++) PXi2 += Rt[2][k]*m_kizyunTen.m_XtilderOrig[i+1][k+1];
	// PXi2 += m_t[2];

	VecDoub PXi(3, 0.0);
	for(int j=0; j<3; j++){
	    for(int k=0; k<4; k++) PXi[j] += P[j][k]*X[i][k];
	}
	// printf("i=%d, PXi2=%g, PXi[2]=%g\n", i, PXi2, PXi[2]);//////////////////

	double invPXi2square {1.0/(PXi[2]*PXi[2])};

	VecDoub dPda1Xi(3, 0.0);
	VecDoub dPda2Xi(3, 0.0);
	VecDoub dPda3Xi(3, 0.0);
	for(int j=0; j<3; j++){
	    for(int k=0; k<3; k++){
		dPda1Xi[j] += KQxdQzQy[j][k]*X[i][k];
		dPda2Xi[j] += KQxQzQyd[j][k]*X[i][k];
		dPda3Xi[j] += KQxQzdQy[j][k]*X[i][k];
	    }
	}
	// for(int j=0; j<3; j++) printf("j=%d, dPda1Xi=%g, dPda2Xi=%g, dPda3Xi=%g\n",
	// 			      j, dPda1Xi[j], dPda2Xi[j], dPda3Xi[j]);///////////
		
	for(int l=0; l<=1; l++){
	    int j {i*2+l};

	    double numer {0.0};
	    for(int k=0; k<4; k++) numer += Rt[l][k]*X[i][k];

	    Jval[j][0] = numer/PXi[2];

	    Jval[j][1] = invPXi2square*( dPda1Xi[l]*PXi[2] -PXi[l]*dPda1Xi[2] );
	    Jval[j][2] = invPXi2square*( dPda2Xi[l]*PXi[2] -PXi[l]*dPda2Xi[2] );
	    Jval[j][3] = invPXi2square*( dPda3Xi[l]*PXi[2] -PXi[l]*dPda3Xi[2] );

	    if(l==0){
		Jval[j][4] = m_alpha/PXi[2];
		Jval[j][5] = 0.0;
		Jval[j][6] = invPXi2square*(m_w2*PXi[2]-PXi[0]);
	    }else{
		Jval[j][4] = 0.0;
		Jval[j][5] = m_alpha/PXi[2];
		Jval[j][6] = invPXi2square*(m_h2*PXi[2]-PXi[1]);
	    }
	}
    }

    // printMatrix(Jval, "EuclidModel::getJ Jval");///////////////////////
}

/*****************************************************/
void EuclidModel::getJ0(VecDoub_I &P, MatDoub_O &Jval){
    double cx {std::cos(P[1])};
    double sx {std::sin(P[1])};
    double cy {std::cos(P[2])};
    double sy {std::sin(P[2])};
    double cz {std::cos(P[3])};
    double sz {std::sin(P[3])};
    
    for(int i=0; i<m_kizyunTen.get_n(); i++){
	double Xi {m_kizyunTen.m_XtilderOrig[i+1][1]};
	double Yi {m_kizyunTen.m_XtilderOrig[i+1][2]};
	double Zi {m_kizyunTen.m_XtilderOrig[i+1][3]};
	double xi {m_kizyunTen.m_xTilderOrig[i+1][1]};
	double yi {m_kizyunTen.m_xTilderOrig[i+1][2]};
	
	double alpha {P[0]};
	
	double p {cy*cz*Xi+(sx*sy*cz+cx*sz)*Yi+(-cx*sy*cz+sx*sz)*Zi+P[4]};
	double q {-cy*sz*Xi+(-sx*sy*sz+cx*cz)*Yi+(cx*sy*sz+sx*sz)*Zi+P[5]};
	// double r {sy*Xi-sx*cy*Yi+cx*cy*Zi+P[6]};

	double dp_dthetax {(cx*sy*cz-sx*sz)*Yi+(sx*sy*cz+cx*sz)*Zi};
	double dp_dthetay {-sy*cz*Xi+sx*cy*cz*Yi-cx*cy*cz*Zi};
	double dp_dthetaz {-cy*sz*Xi+(-sx*sy*sz+cx*cz)*Yi+(cx*sy*sz+sx*cz)*Zi};

	double dq_dthetax {(-cx*sy*sz-sx*cz)*Yi+(-sx*sy*sz+cx*cz)*Zi};
	double dq_dthetay {sy*sz*Xi-sx*cy*sz*Yi+cx*cy*sz*Zi};
	double dq_dthetaz {-cy*cz*Xi-(sx*sy*cz+cx*sz)*Yi+(cx*sy*cz-sx*sz)*Zi};

	double dr_dthetax {-cx*cy*Yi-sx*cy*Zi};
	double dr_dthetay {cy*Xi+sx*sy*Yi-cx*sy*Zi};
	// double dr_dthetaz {0.0};

	Jval[2*i][0] = p;  // d alpha
	Jval[2*i+1][0] = q;
	Jval[2*i][1] = alpha*dp_dthetax-(xi-m_w2)*dr_dthetax;  // d theta_x
	Jval[2*i+1][1] = alpha*dq_dthetax-(yi-m_h2)*dr_dthetax;
	Jval[2*i][2] = alpha*dp_dthetay-(xi-m_w2)*dr_dthetay;  // d theta_y
	Jval[2*i+1][2] = alpha*dq_dthetay-(yi-m_h2)*dr_dthetay;
	Jval[2*i][3] = alpha*dp_dthetaz;  // d theta_z
	Jval[2*i+1][3] = alpha*dq_dthetaz;
	Jval[2*i][4] = alpha;  // d t1
	Jval[2*i+1][4] = 0.0;
	Jval[2*i][5] = 0.0;  // d t2
	Jval[2*i+1][5] = alpha;
	Jval[2*i][6] = m_w2-xi;  // d t3
	Jval[2*i+1][6] = m_h2-yi;
    }
}

/** for LM algorithm. Calculate Jacobian numerically. ***/
void EuclidModel::getJnum(VecDoub &vecAlpha, MatDoub_O &Jval){
    int fdim {Jval.nrows()};
    VecDoub fval0(fdim);
    VecDoub fval1(fdim);
    
    get_F(vecAlpha, fval0);
	    
    for(int j=0; j<vecAlpha.size(); j++){
        // double delta {std::max(std::abs(1.0E-4*P[j]), 1.0E-6)};
        double delta {std::max(std::abs(1.0E-8*vecAlpha[j]), 1.0E-10)};

        double bak {vecAlpha[j]};
	
        vecAlpha[j] += delta;
        get_F(vecAlpha, fval1);
	    
        for(int i=0; i<fdim; i++){
            Jval[i][j] = (fval1[i] - fval0[i])/delta;
	    // if(i==1 && j==1){
	    // 	printf("fval1=%20f, fval0=%20f, delta=%g\n", fval1[i], fval0[i], delta);//////
	    // }
        }

        vecAlpha[j] = bak;
    }
}

/** for LM algorithm. Calculate Jacobian numerically. ***/
void EuclidModel::getJnum0(VecDoub &P, MatDoub_O &Jval){
    int fdim {Jval.nrows()};
    VecDoub fval0(fdim);
    VecDoub fval1(fdim);
    
    get_f0(P, fval0);
	    
    for(int j=0; j<P.size(); j++){
        // double delta {std::max(std::abs(1.0E-4*P[j]), 1.0E-6)};
        double delta {std::max(std::abs(1.0E-8*P[j]), 1.0E-10)};

        double bak {P[j]};
	
        P[j] += delta;
        get_f0(P, fval1);
	    
        for(int i=0; i<fdim; i++){
            Jval[i][j] = (fval1[i] - fval0[i])/delta;
	    if(i==1 && j==1){
		printf("fval1=%20f, fval0=%20f, delta=%g\n", fval1[i], fval0[i], delta);//////
	    }
        }

        P[j] = bak;
    }
}

/************ My notebook page 22.6-10 ************************************/
void EuclidModel::get_F(VecDoub_I &vecAlpha, VecDoub_O &fval){
    // printf("EuclidModel::get_f\n");///////////////////

    calcPgene(vecAlpha, false);   // calculate m_P[3][4]
    
    for(int i=0; i<m_kizyunTen.get_n(); i++){
	double PXi[3];
	for(int j=0; j<3; j++){
	    PXi[j] = 0.0;
	    for(int k=0; k<3; k++){
		PXi[j] += m_P[j][k]*m_kizyunTen.m_XtilderOrig[i][k];
	    }
	    PXi[j] += m_P[j][3];
	}

	fval[i*2] = PXi[0]/PXi[2];
	fval[i*2+1] = PXi[1]/PXi[2];
    }
}

/**** Calculate the Pprime.  MyNote 22.10-4...11  **************************/
void EuclidModel::calcPprime2(KizyunTen &ktRight, VecDoub_O &vecAlpha){
    MatDoub A(4, 2);
    // MatDoub A(2, 2);
    VecDoub b(A.nrows());
    // printf("m_alpha=%g, m_w2=%g, m_h2=%g\n", m_alpha, m_w2, m_h2);/////////////////////

    setAb(ktRight, A, b, 2, 1, 0); //[2:he, 1:hw]
    setAb(ktRight, A, b, 4, 3, 1); //[4:le, 3:lw]
    setAb(ktRight, A, b, 6, 5, 2); //[6:me, 5:mw]
    setAb(ktRight, A, b, 8, 7, 3); //[8:mhe, 7:mhw]

    // for(int i=0; i<A.nrows(); i++){//////////////////////
    // 	for(int j=0; j<A.ncols(); j++){
    // 	    printf("A[%d][%d]=%g ", i, j, A[i][j]);
    // 	}
    // 	printf("\n");
    // }
    // for(int i=0; i<b.size(); i++){//////////////////////
    // 	printf("b[%d]=%g\n", i, b[i]);
    // }

    MatDoub Apinv(A.ncols(), A.nrows());
    calcPseudInv(A, Apinv);
    // for(int i=0; i<Apinv.nrows(); i++){//////////////////////
    // 	for(int j=0; j<Apinv.ncols(); j++){
    // 	    printf("Apinv[%d][%d]=%g ", i, j, Apinv[i][j]);
    // 	}
    // 	printf("\n");
    // }
    // for(int i=0; i<Apinv.nrows(); i++){//////////////////////
    // 	for(int j=0; j<A.ncols(); j++){
    // 	    double sum {0.0};
    // 	    for(int k=0; k<Apinv.ncols(); k++) sum += Apinv[i][k]*A[k][j];
    // 	    printf("Apinv*A[%d][%d]=%g ", i, j, sum);
    // 	}
    // 	printf("\n");
    // }

    VecDoub Apinv_b(Apinv.nrows());
    for(int i=0; i<Apinv.nrows(); i++){
	double sum {0.0};
	for(int j=0; j<Apinv.ncols(); j++) sum += Apinv[i][j]*b[j];
	Apinv_b[i] = sum;
    }

    // for(int i=0; i<Apinv_b.size(); i++){/////////////////////////
    // 	printf("Apinv_b[%d] = %g\n", i, Apinv_b[i]);
    // }

    // double theta_x {std::atan(Apinv_b[0])};
    m_thetax = std::atan(Apinv_b[0]);
    GlobalConstants gc;
    // printf("theta_x=%g (deg)\n", theta_x*180.0/gc.pi);//////////////////

    // double theta_y {std::atan( 1.0/(Apinv_b[1]*std::cos(theta_x)) )};
    m_thetay = std::atan( 1.0/(Apinv_b[1]*std::cos(m_thetax)) );
    // printf("theta_y=%g (deg)\n", theta_y*180.0/gc.pi);//////////////////

    // Next, estimate the vector t
    MatDoub A2(ktRight.get_n()*2, 3);
    VecDoub b2(ktRight.get_n()*2);

    double cx {std::cos(m_thetax)};
    double sx {std::sin(m_thetax)};
    double cy {std::cos(m_thetay)};
    double sy {std::sin(m_thetay)};
    double w2 {ktRight.m_width*0.5};
    double h2 {ktRight.m_height*0.5};
    for(int i=0; i<ktRight.get_n(); i++){
	double xi {ktRight.m_xTilderOrig[i][1]};
	double yi {ktRight.m_xTilderOrig[i][2]};
	double Xi {ktRight.m_XtilderOrig[i][1]};
	double Yi {ktRight.m_XtilderOrig[i][2]};
	double Zi {ktRight.m_XtilderOrig[i][3]};
	double ri0 {cy*Xi-sy*Zi};
	double ri1 {sx*sy*Xi+cx*Yi+sx*cy*Zi};
	double ri2 {cx*sy*Xi-sx*Yi+cx*cy*Zi};
	
	A2[i*2][0] = 0.0; A2[i*2][1] = -m_alpha; A2[i*2][2] = yi-h2;
	A2[i*2+1][0] = m_alpha; A2[i*2+1][1] = 0.0; A2[i*2+1][2] = w2-xi;

	b2[i*2] = m_alpha*ri1-(yi-h2)*ri2;
	b2[i*2+1] = -m_alpha*ri0+(xi-w2)*ri2;
    }

    MatDoub A2pinv(A2.ncols(), A2.nrows());
    calcPseudInv(A2, A2pinv);
    // for(int i=0; i<A2pinv.nrows(); i++){//////////////////////
    // 	for(int j=0; j<A2pinv.ncols(); j++){
    // 	    printf("A2pinv[%d][%d]=%g ", i, j, A2pinv[i][j]);
    // 	}
    // 	printf("\n");
    // }

    // VecDoub t(A2pinv.nrows());
    for(int i=0; i<A2pinv.nrows(); i++){
	double sum {0.0};
	for(int j=0; j<A2pinv.ncols(); j++) sum += A2pinv[i][j]*b2[j];
	m_t[i] = sum;
    }

    // for(int i=0; i<t.size(); i++){/////////////////////////
    // 	printf("t[%d] = %g\n", i, t[i]);
    // }

    m_thetaz = 0.0;
    
    vecAlpha[0] = m_alpha;
    vecAlpha[1] = m_thetax;
    vecAlpha[2] = m_thetay;
    vecAlpha[3] = m_thetaz;
    vecAlpha[4] = m_t[0];
    vecAlpha[5] = m_t[1];
    vecAlpha[6] = m_t[2];
    for(int i=0; i<vecAlpha.size(); i++){///////////////////////////
	if(i>=1 && i<=3) printf("vecAlpha[%d]=%g (deg)\n", i, vecAlpha[i]*180.0/gc.pi);
	else printf("vecAlpha[%d]=%g\n", i, vecAlpha[i]);
    }
    
    // Calculate a camera matrix Pprime
    MatDoub K(3, 3, 0.0);
    K[0][0] = m_alpha; K[0][2] = w2;
    K[1][1] = m_alpha; K[1][2] = h2;
    K[2][2] = 1.0;
    
    MatDoub Rt(3, 4);
    Rt[0][0]=cy; Rt[0][1]=0.0; Rt[0][2]=-sy; Rt[0][3]=m_t[0];
    Rt[1][0]=sx*sy; Rt[1][1]=cx; Rt[1][2]=sx*cy; Rt[1][3]=m_t[1];
    Rt[2][0]=cx*sy; Rt[2][1]=-sx; Rt[2][2]=cx*cy; Rt[2][3]=m_t[2];
    
    // MatDoub Pprime(3, 4);
    for(int i=0; i<3; i++){
	for(int j=0; j<4; j++){
	    double sum {0.0};
	    for(int k=0; k<3; k++) sum += K[i][k]*Rt[k][j];
	    m_P[i][j] = sum;
	}
    }

    // for(int i=0; i<3; i++){////////////////////////////////////
    // 	for(int j=0; j<4; j++){
    // 	    printf("Pp[%d][%d]=%g ", i, j, m_Pprime[i][j]/m_Pprime[2][3]);////////////////
    // 	}
    // 	printf("\n");//////////////////////
    // }
    
    // // For check, calculate the camera center
    // MatDoub M(3, 3);
    // for(int i=0; i<3; i++){
    // 	for(int j=0; j<3; j++) M[i][j] = m_Pprime[i][j];
    // }

    // LUdcmp lud(M);
    // MatDoub MInv(3, 3);
    // lud.inverse(MInv);
    
    // VecDoub CC(3, 0.0);
    // for(int i=0; i<3; i++){
    // 	for(int k=0; k<3; k++) CC[i] += -MInv[i][k]*m_Pprime[k][3];
    // }
    
    // for(int i=0; i<3; i++){////////////////
    // 	printf("CC[%d]=%g\n", i, CC[i]);
    // }

}

/*******************************************/
void EuclidModel::setAb(KizyunTen &ktRight, MatDoub_IO &A, VecDoub_IO &b,
			int i, int j, int Ai){
    double xi {ktRight.m_xTilderOrig[i][0]};
    double yi {ktRight.m_xTilderOrig[i][1]};
    double xj {ktRight.m_xTilderOrig[j][0]};
    double yj {ktRight.m_xTilderOrig[j][1]};
    double w2 {ktRight.m_width*0.5};
    double h2 {ktRight.m_height*0.5};
    
    A[Ai][0] = xi - xj;
    A[Ai][1] = yj - yi;
    b[Ai] = ( (xi-xj)*(yj-h2) - (xj-w2)*(yi-yj) )/m_alpha;
}

/*********************************************/
void EuclidModel::calcPseudInv(MatDoub_I &A, MatDoub_O &Apinv){
    int m {A.nrows()};
    int n {A.ncols()};
    MatDoub ATA(n, n);

    for(int i=0; i<n; i++){
	for(int j=0; j<n; j++){
	    double sum {0.0};
	    for(int k=0; k<m; k++) sum += A[k][i]*A[k][j];
	    ATA[i][j] = sum;
	}
    }

    LUdcmp lud(ATA);
    MatDoub ATAinv(n, n);
    lud.inverse(ATAinv);

    for(int i=0; i<n; i++){
	for(int j=0; j<m; j++){
	    double sum {0.0};
	    for(int k=0; k<n; k++) sum += ATAinv[i][k]*A[j][k];
	    Apinv[i][j] = sum;
	}
    }
}

/**** Calculate the Pprime.  MyNote 22.10-1,2  **************************/
void EuclidModel::calcPprime(KizyunTen &ktRight){
    // int i {2};   // lw
    // int j {3};   // le
    // int i {0};   // hw
    // int j {1};   // he
    int i {4};   // mw
    int j {5};   // me

    double xi {ktRight.m_xTilderOrig[i][0]};
    double yi {ktRight.m_xTilderOrig[i][1]};
    double xj {ktRight.m_xTilderOrig[j][0]};
    double yj {ktRight.m_xTilderOrig[j][1]};
    
    // printf("xi=%g, xj=%g\n", xi, xj);/////////////////
    // printf("yi=%g, yj=%g\n", yi, yj);////////////////

    double w2 {ktRight.m_width*0.5};
    double h2 {ktRight.m_height*0.5};
    // printf("w2=%g, h2=%g\n", w2, h2);///////////////

    double tanTheta_y {m_alpha*(yi-yj)/( (xj-w2)*(yi-yj)-(yj-h2)*(xi-xj) )};
    printf("tanTheta_y = %g\n", tanTheta_y);//////////////////////

    double theta_y {std::atan(tanTheta_y)};
    printf("theta_y=%g (rad)\n", theta_y);

    GlobalConstants gc;
    printf("theta_y=%g (deg)\n", theta_y*180.0/gc.pi);
}

/************ Not used. linear ************************************/
void EuclidModel::get_f0(VecDoub_I &P, VecDoub_O &fval){
    double cx {std::cos(P[1])};
    double sx {std::sin(P[1])};
    double cy {std::cos(P[2])};
    double sy {std::sin(P[2])};
    double cz {std::cos(P[3])};
    double sz {std::sin(P[3])};

    for(int i=0; i<m_kizyunTen.get_n(); i++){
	double Xi {m_kizyunTen.m_XtilderOrig[i+1][1]};
	double Yi {m_kizyunTen.m_XtilderOrig[i+1][2]};
	double Zi {m_kizyunTen.m_XtilderOrig[i+1][3]};
	double xi {m_kizyunTen.m_xTilderOrig[i+1][1]};
	double yi {m_kizyunTen.m_xTilderOrig[i+1][2]};
	
	double p {cy*cz*Xi+(sx*sy*cz+cx*sz)*Yi+(-cx*sy*cz+sx*sz)*Zi+P[4]};
	double q {-cy*sz*Xi+(-sx*sy*sz+cx*cz)*Yi+(cx*sy*sz+sx*sz)*Zi+P[5]};
	double r {sy*Xi-sx*cy*Yi+cx*cy*Zi+P[6]};

	fval[2*i] = P[0]*p-(xi-m_w2)*r;
	fval[2*i+1] = P[0]*q-(yi-m_h2)*r;
    }
}

/*** calcP 
 *  my note: Section 22.7  *****************************************/
void EuclidModel::calcSimple(){
    // printf("width=%d\n", m_kizyunTen.m_width);//////////////////
    double ZcHalpha = (m_kizyunTen.m_Hw[0]-m_kizyunTen.m_He[0])
	/(m_kizyunTen.m_hw[0]-m_kizyunTen.m_he[0]);
    double ZcLalpha = (m_kizyunTen.m_Lw[0]-m_kizyunTen.m_Le[0])
	/(m_kizyunTen.m_lw[0]-m_kizyunTen.m_le[0]);

    // printf("ZcHalpha=%g, ZcLalpha=%g\n", ZcHalpha, ZcLalpha);////////////////

    GlobalConstants gc;
    double d {ZcHalpha*(m_kizyunTen.m_he[1]-0.5*m_kizyunTen.m_height)
    	      -ZcLalpha*(m_kizyunTen.m_le[1]-0.5*m_kizyunTen.m_height)};
    // double d {ZcHalpha*(m_kizyunTen.m_hw[1]-0.5*m_kizyunTen.m_height)
    // 	      -ZcLalpha*(m_kizyunTen.m_lw[1]-0.5*m_kizyunTen.m_height)};
    // printf("d=%g\n", d);////////////////////

    double theta_rad;   // radian
    double dr {d/gc.mY};
    if(dr >1.0) theta_rad = (90.0-gc.theta_p)*gc.pi/180.0;
    else if(dr < -1.0) theta_rad = (-90.0-gc.theta_p)*gc.pi/180.0;
    else theta_rad = std::asin(dr) -gc.theta_p*gc.pi/180.0;
    // printf("EuclidModel::calcSimple(): theta= %g (deg)\n", theta_rad*180.0/gc.pi);//////////////////////////

    m_c = std::cos(theta_rad);
    m_s = std::sin(theta_rad);
    m_theta_rad = theta_rad;
    // printf("m_c=%g, m_s=%g\n", m_c, m_s);//////////////////

    // Compute t1, t2 (m_t[0], m_t[1])
    double wp2 {m_kizyunTen.m_width*0.5};
    double hp2 {m_kizyunTen.m_height*0.5};
    MatDoub A(m_kizyunTen.get_n(), 2);

    for(int i=0; i<m_kizyunTen.get_n(); i++){
	A[i][0] = m_kizyunTen.m_xTilderOrig[i][1]-hp2;
	A[i][1] = wp2 -m_kizyunTen.m_xTilderOrig[i][0];
    }
    
    VecDoub b(m_kizyunTen.get_n());
    for(int i=0; i<m_kizyunTen.get_n(); i++){
	b[i] = m_kizyunTen.m_XtilderOrig[i][1]
	    *(m_kizyunTen.m_xTilderOrig[i][0]-wp2)*m_c
	    +m_kizyunTen.m_XtilderOrig[i][2]
	    *(m_kizyunTen.m_xTilderOrig[i][0]-wp2)*m_s
	    -(m_kizyunTen.m_xTilderOrig[i][1]-hp2)
	    *m_kizyunTen.m_XtilderOrig[i][0];
    }
	
    SVD svd(A);
    // for(int i=0; i<svd.w.size(); i++){ ////////////////////////
    // 	printf("w[%d]=%g\n", i, svd.w[i]);
    // }
    
    MatDoub Aplus(A.ncols(), A.nrows());  // pseudo inverse
    for(int i=0; i<Aplus.nrows(); i++){
	for(int j=0; j<Aplus.ncols(); j++){
	    Aplus[i][j] = 0.0;
	    for(int k=0; k<svd.w.size(); k++){
		if(svd.w[k] != 0.0){
		    Aplus[i][j] += svd.v[i][k]*(1.0/svd.w[k])*svd.u[j][k];
		}
	    }
	}
    }
    
    // printf("Aplus x A\n");////////////////////
    // for(int i=0; i<Aplus.nrows(); i++){//////////////////
    // 	for(int j=0; j<A.ncols(); j++){
    // 	    double x {0.0};
    // 	    for(int k=0; k<Aplus.ncols(); k++){
    // 		x += Aplus[i][k]*A[k][j];
    // 	    }
    // 	    printf("%8f ", x);
    // 	}
    // 	printf("\n");
    // }

    VecDoub Aplusb(Aplus.nrows());
    for(int i=0; i<Aplusb.size(); i++){
	Aplusb[i] = 0.0;
	for(int j=0; j<Aplus.ncols(); j++){
	    Aplusb[i] += Aplus[i][j]*b[j];
	}
    }
    // printf("Aplusb\n");////////////////////
    // for(int i=0; i<Aplusb.size(); i++){//////////////////
    // 	printf("Aplusb[%d]=%g\n", i, Aplusb[i]);
    // }

    m_t[0] = Aplusb[0];
    m_t[1] = Aplusb[1];

    // printf("EuclidModel::calcSimple(): t1=%g\n", m_t[0]);//////////////////
    // printf("EuclidModel::calcSimple(): t2=%g\n", m_t[1]);//////////////////

    // Get the initial guess of alpha and t3==m_t[2]
    MatDoub A2(m_kizyunTen.get_n()*2, 2);
    VecDoub b2(m_kizyunTen.get_n()*2);

    for(int i=0; i<m_kizyunTen.get_n(); i++){
	A2[2*i][0] = m_kizyunTen.m_xTilderOrig[i][1]-hp2;
	A2[2*i][1] = -(m_c*m_kizyunTen.m_XtilderOrig[i][1]
		       +m_s*m_kizyunTen.m_XtilderOrig[i][2] +m_t[1]);
	A2[2*i+1][0] = wp2-m_kizyunTen.m_xTilderOrig[i][0];
	A2[2*i+1][1] = m_kizyunTen.m_XtilderOrig[i][0] +m_t[0];

	b2[2*i] = (m_kizyunTen.m_xTilderOrig[i][1] -hp2)
	    *(m_s*m_kizyunTen.m_XtilderOrig[i][1]
	      -m_c*m_kizyunTen.m_XtilderOrig[i][2]);
	b2[2*i+1] = (m_kizyunTen.m_xTilderOrig[i][0] -wp2)
	    *(-m_s*m_kizyunTen.m_XtilderOrig[i][1]
	      +m_c*m_kizyunTen.m_XtilderOrig[i][2]);
    }

    // printf("A2 b\n");
    // for(int i=0; i<A2.nrows(); i++){//////////////////
    // 	for(int j=0; j<A2.ncols(); j++){
    // 	    printf("A2[%d][%d]=%8f ", i, j, A2[i][j]);
    // 	}
    // 	printf("b2[%d]=%8f", i, b2[i]);
    // 	printf("\n");
    // }
    
    SVD svd2(A2);
    // for(int i=0; i<svd2.w.size(); i++){ ////////////////////////
    // 	printf("w[%d]=%g\n", i, svd2.w[i]);
    // }

    MatDoub A2plus(A2.ncols(), A2.nrows());  // pseudo inverse
    for(int i=0; i<A2plus.nrows(); i++){
	for(int j=0; j<A2plus.ncols(); j++){
	    A2plus[i][j] = 0.0;
	    for(int k=0; k<2; k++){
		if(svd2.w[k] != 0.0){
		    A2plus[i][j] += svd2.v[i][k]*(1.0/svd2.w[k])*svd2.u[j][k];
		}
	    }
	}
    }

    // printf("A2plus\n");////////////////////
    // for(int i=0; i<A2plus.nrows(); i++){//////////////////
    // 	for(int j=0; j<A2plus.ncols(); j++){
    // 	    printf("%8f ", A2plus[i][j]);
    // 	}
    // 	printf("\n");
    // }    
    // printf("A2plus x A2\n");////////////////////
    // for(int i=0; i<A2plus.nrows(); i++){//////////////////
    // 	for(int j=0; j<A2.ncols(); j++){
    // 	    double x {0.0};
    // 	    for(int k=0; k<A2plus.ncols(); k++){
    // 		x += A2plus[i][k]*A2[k][j];
    // 	    }
    // 	    printf("%8f ", x);
    // 	}
    // 	printf("\n");
    // }

    VecDoub A2plusb2(A2plus.nrows());
    for(int i=0; i<A2plusb2.size(); i++){
	A2plusb2[i] = 0.0;
	for(int j=0; j<A2plus.ncols(); j++){
	    A2plusb2[i] += A2plus[i][j]*b2[j];
	}
    }
    // printf("A2plusb2\n");////////////////////
    // for(int i=0; i<A2plusb2.size(); i++){//////////////////
    // 	printf("A2plusb2[%d]=%g\n", i, A2plusb2[i]);
    // }

    m_t[2] = A2plusb2[0];
    m_alpha = A2plusb2[1];

    // printf("EuclidModel::calcSimple(): t3=%g\n", m_t[2]);/////////////////
    // printf("EuclidModel::calcSimple(): alpha=%g\n", m_alpha);/////////////////

}

/*** calcP 
 *  Calculate the camera matrix P
 *  Use: m_alpha, m_kizyunTen.m_width, m_kizyunTen.m_height, m_c, m_s, m_t[3]
 *  Output: m_P[3][4]
 *  my note: 22.6-7a page  *****************************************/
void EuclidModel::calcP(){
    double w2 {m_kizyunTen.m_width*0.5};
    double h2 {m_kizyunTen.m_height*0.5};
    
    m_P[0][0] = m_alpha;
    m_P[0][1] = -w2*m_s;
    m_P[0][2] = w2*m_c;
    m_P[0][3] = m_alpha*m_t[0]+w2*m_t[2];
    
    m_P[1][0] = 0.0;
    m_P[1][1] = m_alpha*m_c-h2*m_s;
    m_P[1][2] = m_alpha*m_s+h2*m_c;
    m_P[1][3] = m_alpha*m_t[1]+h2*m_t[2];
    
    m_P[2][0] = 0.0;
    m_P[2][1] = -m_s;
    m_P[2][2] = m_c;
    m_P[2][3] = m_t[2];

    // // for check
    // printf("EuclidModel::calcP():  P[3][4]\n");///////////////////
    // for(int i=0; i<=2; i++){
    // 	for(int j=0; j<=3; j++){
    // 	    printf("%12f ", m_P[i][j]/m_P[2][3]);
    // 	}
    // 	printf("\n");
    // }

    // Calculate the camera center
    MatDoub M(3, 3);
    for(int i=0; i<=2; i++){
	for(int j=0; j<=2; j++){
	    M[i][j] = m_P[i][j];
	}
    }

    LUdcmp lud(M);
    MatDoub Minv(3, 3);
    lud.inverse(Minv);

    for(int i=0; i<=2; i++){
	m_Ctilder[i] = 0.0;
	for(int j=0; j<=2; j++){
	    m_Ctilder[i] -= Minv[i][j]*m_P[j][3];
	}
    }
    
    // // for check
    // printf("EuclidModel::calcP(): Ctilder[3], camera center\n");///////////////////
    // for(int i=0; i<=2; i++){
    // 	printf("EuclidModel::calcP(): Ctilder[%d]= %12f\n", i, m_Ctilder[i]);
    // }
}

/*** calcPgene 
 *  Calculate the camera matrix P in general
 *  Input: P[0]:alpha, P[1]:thetax, P[2]:thetay, P[3]:thetaz, P[4]:t1, P[5]:t2, P[6]:t3
 *  Use: m_kizyunTen.m_width, m_kizyunTen.m_height
 *  Output: m_P[3][4], m_alpha, m_theta[x,y,z], m_t[3]
 *  my note: 22.6-10 page  *****************************************/
void EuclidModel::calcPgene(VecDoub_I &P, bool boolPrint){
    double m_w2 {m_kizyunTen.m_width*0.5};
    double m_h2 {m_kizyunTen.m_height*0.5};

    m_alpha = P[0];
    m_thetax = P[1];
    m_thetay = P[2];
    m_thetaz = P[3];
    m_t[0] = P[4];
    m_t[1] = P[5];
    m_t[2] = P[6];
    double cx {std::cos(P[1])};
    double sx {std::sin(P[1])};
    double cy {std::cos(P[2])};
    double sy {std::sin(P[2])};
    double cz {std::cos(P[3])};
    double sz {std::sin(P[3])};
    
    // MatDoub Rt(3, 4);
    // Rt[0][0] = cy*cz;
    // Rt[0][1] = sx*sy*cz+cx*sz;
    // Rt[0][2] = -cx*sy*cz+sx*sz;
    // Rt[0][3] = m_t[0];
    // Rt[1][0] = -cy*sz;
    // Rt[1][1] = -sx*sy*sz+cx*cz;
    // Rt[1][2] = cx*sy*sz+sx*cz;
    // Rt[1][3] = m_t[1];
    // Rt[2][0] = sy;
    // Rt[2][1] = -sx*cy;
    // Rt[2][2] = cx*cy;
    // Rt[2][3] = m_t[2];

    MatDoub K(3, 3, 0.0);
    MatDoub Qx(3, 3, 0.0);
    MatDoub Qy(3, 3, 0.0);
    MatDoub Qz(3, 3, 0.0);
    MatDoub Rt(3, 4);
    
    MatDoub A(3, 3);  // temporaly matrix used when computing matrix product

    K[0][0] = m_alpha;
    K[0][2] = m_w2;
    K[1][1] = m_alpha;
    K[1][2] = m_h2;
    K[2][2] = 1.0;
    
    Qx[0][0] = 1.0;
    Qx[1][1] = cx;
    Qx[1][2] = sx;
    Qx[2][1] = -sx;
    Qx[2][2] = cx;
    
    Qy[0][0] = cy;
    Qy[0][2] = -sy;
    Qy[1][1] = 1.0;
    Qy[2][0] = sy;
    Qy[2][2] = cy;
    
    Qz[0][0] = cz;
    Qz[0][1] = sz;
    Qz[1][0] = -sz;
    Qz[1][1] = cz;
    Qz[2][2] = 1.0;
    
    prod(Qx, Qz, A);
    prod(A, Qy, Rt);

    for(int i=0; i<3; i++) Rt[i][3] = m_t[i];
    // printMatrix(Rt, "Rt");//////////////////

    prod(K, Rt, m_P);
    
    // for(int i=0; i<3; i++){
    // 	for(int j=0; j<4; j++){
    // 	    double sum {0.0};
    // 	    for(int k=0; k<3; k++) sum += K[i][k]*Rt[k][j];
    // 	    m_P[i][j] = sum;
    // 	}
    // }
    
    // for check
    if(boolPrint==true){
	printf("EuclidModel::calcPgene():  P[3][4]\n");////////////////////
	for(int i=0; i<=2; i++){
	    for(int j=0; j<=3; j++){
		printf("%12f ", m_P[i][j]/m_P[2][3]);
	    }
	    printf("\n");
	}
    }
    
    // Calculate the camera center
    MatDoub M(3, 3);
    for(int i=0; i<=2; i++){
	for(int j=0; j<=2; j++){
	    M[i][j] = m_P[i][j];
	}
    }

    LUdcmp lud(M);
    MatDoub Minv(3, 3);
    lud.inverse(Minv);

    for(int i=0; i<=2; i++){
	m_Ctilder[i] = 0.0;
	for(int j=0; j<=2; j++){
	    m_Ctilder[i] -= Minv[i][j]*m_P[j][3];
	}
    }
    
    // for check
    if(boolPrint==true){
	printf("EuclidModel::calcPgene(): Ctilder[3], camera center\n");///////////////////
	for(int i=0; i<=2; i++){
	    printf("EuclidModel::calcPgene(): Ctilder[%d]= %12f\n", i, m_Ctilder[i]);
	}
    }
}

/*** calc *****************************************/
void EuclidModel::calc(){
    // printf("EuclicModel::calc, width=%d, height=%d\n", m_kizyunTen.m_width, m_kizyunTen.m_height);//////////////////

    double wp2 {m_kizyunTen.m_width*0.5};
    double hp2 {m_kizyunTen.m_height*0.5};
    MatDoub A(m_kizyunTen.get_n(), 4);

    for(int i=0; i<m_kizyunTen.get_n(); i++){
	A[i][0] = (m_kizyunTen.m_xTilderOrig[i][0]-wp2)*m_kizyunTen.m_XtilderOrig[i][1];
	A[i][1] = (m_kizyunTen.m_xTilderOrig[i][0]-wp2)*m_kizyunTen.m_XtilderOrig[i][2];
	A[i][2] = hp2-m_kizyunTen.m_xTilderOrig[i][1];
	A[i][3] = m_kizyunTen.m_xTilderOrig[i][0]-wp2;
    }

    // for(int i=0; i<m_kizyunTen.get_n(); i++){//////////////////
    // 	for(int j=0; j<4; j++){
    // 	    printf("%8f ", A[i][j]);
    // 	}
    // 	printf("\n");
    // }    

    VecDoub b(m_kizyunTen.get_n());
    for(int i=0; i<m_kizyunTen.get_n(); i++){
	b[i] = (m_kizyunTen.m_xTilderOrig[i][1]-hp2)*m_kizyunTen.m_XtilderOrig[i][0];
    }

    SVD svd(A);
    // for(int i=0; i<svd.w.size(); i++){ ////////////////////////
    // 	printf("w[%d]=%g\n", i, svd.w[i]);
    // }
    // exit(0);///////////////////

    MatDoub Aplus(A.ncols(), A.nrows());  // pseudo inverse
    for(int i=0; i<Aplus.nrows(); i++){
	for(int j=0; j<Aplus.ncols(); j++){
	    Aplus[i][j] = 0.0;
	    for(int k=0; k<4; k++){
		if(svd.w[k] != 0.0){
		    Aplus[i][j] += svd.v[i][k]*(1.0/svd.w[k])*svd.u[j][k];
		}
	    }
	}
    }

    // printf("Aplus\n");////////////////////
    // for(int i=0; i<Aplus.nrows(); i++){//////////////////
    // 	for(int j=0; j<Aplus.ncols(); j++){
    // 	    printf("%8f ", A[i][j]);
    // 	}
    // 	printf("\n");
    // }    
    // printf("Aplus x A\n");////////////////////
    // for(int i=0; i<Aplus.nrows(); i++){//////////////////
    // 	for(int j=0; j<A.ncols(); j++){
    // 	    double x {0.0};
    // 	    for(int k=0; k<Aplus.ncols(); k++){
    // 		x += Aplus[i][k]*A[k][j];
    // 	    }
    // 	    printf("%8f ", x);
    // 	}
    // 	printf("\n");
    // }

    VecDoub Aplusb(Aplus.nrows());
    for(int i=0; i<Aplusb.size(); i++){
	Aplusb[i] = 0.0;
	for(int j=0; j<Aplus.ncols(); j++){
	    Aplusb[i] += Aplus[i][j]*b[j];
	}
    }
    // printf("Aplusb\n");////////////////////
    // for(int i=0; i<Aplusb.size(); i++){//////////////////
    // 	printf("Aplusb[%d]=%g\n", i, Aplusb[i]);
    // }

    m_c = Aplusb[0];
    m_s = Aplusb[1];
    m_t[0] = Aplusb[2];
    m_t[1] = Aplusb[3];
    // printf("EuclidModel::calc(): c=%g\n", m_c);//////////////////////
    // printf("EuclidModel::calc(): s=%g\n", m_s);//////////////////////
    // printf("EuclidModel::calc(): c^2+s^2=%g\n", m_c*m_c+m_s*m_s);//////////////////////
    // printf("EuclidModel::calc(): t1=%g\n", m_t[0]);//////////////////////
    // printf("EuclidModel::calc(): t2=%g\n", m_t[1]);//////////////////////

    // Get the initial guess of alpha and t3==m_t[2]
    MatDoub A2(m_kizyunTen.get_n()*2, 2);
    VecDoub b2(m_kizyunTen.get_n()*2);

    for(int i=0; i<m_kizyunTen.get_n(); i++){
	A2[2*i][0] = m_kizyunTen.m_xTilderOrig[i][1]-hp2;
	A2[2*i][1] = -(m_c*m_kizyunTen.m_XtilderOrig[i][1]
		       +m_s*m_kizyunTen.m_XtilderOrig[i][2] +m_t[1]);
	A2[2*i+1][0] = wp2-m_kizyunTen.m_xTilderOrig[i][0];
	A2[2*i+1][1] = m_kizyunTen.m_XtilderOrig[i][0] +m_t[0];

	b2[2*i] = (m_kizyunTen.m_xTilderOrig[i][1] -hp2)
	    *(m_s*m_kizyunTen.m_XtilderOrig[i][1]
	      -m_c*m_kizyunTen.m_XtilderOrig[i][2]);
	b2[2*i+1] = (m_kizyunTen.m_xTilderOrig[i][0] -wp2)
	    *(-m_s*m_kizyunTen.m_XtilderOrig[i][1]
	      +m_c*m_kizyunTen.m_XtilderOrig[i][2]);
    }

    // printf("A2 b\n");
    // for(int i=0; i<A2.nrows(); i++){//////////////////
    // 	for(int j=0; j<A2.ncols(); j++){
    // 	    printf("A2[%d][%d]=%8f ", i, j, A2[i][j]);
    // 	}
    // 	printf("b2[%d]=%8f", i, b2[i]);
    // 	printf("\n");
    // }
    
    SVD svd2(A2);
    // for(int i=0; i<svd2.w.size(); i++){ ////////////////////////
    // 	printf("w[%d]=%g\n", i, svd2.w[i]);
    // }

    MatDoub A2plus(A2.ncols(), A2.nrows());  // pseudo inverse
    for(int i=0; i<A2plus.nrows(); i++){
	for(int j=0; j<A2plus.ncols(); j++){
	    A2plus[i][j] = 0.0;
	    for(int k=0; k<2; k++){
		if(svd2.w[k] != 0.0){
		    A2plus[i][j] += svd2.v[i][k]*(1.0/svd2.w[k])*svd2.u[j][k];
		}
	    }
	}
    }

    // printf("A2plus\n");////////////////////
    // for(int i=0; i<A2plus.nrows(); i++){//////////////////
    // 	for(int j=0; j<A2plus.ncols(); j++){
    // 	    printf("%8f ", A2plus[i][j]);
    // 	}
    // 	printf("\n");
    // }    
    // printf("A2plus x A2\n");////////////////////
    // for(int i=0; i<A2plus.nrows(); i++){//////////////////
    // 	for(int j=0; j<A2.ncols(); j++){
    // 	    double x {0.0};
    // 	    for(int k=0; k<A2plus.ncols(); k++){
    // 		x += A2plus[i][k]*A2[k][j];
    // 	    }
    // 	    printf("%8f ", x);
    // 	}
    // 	printf("\n");
    // }

    VecDoub A2plusb2(A2plus.nrows());
    for(int i=0; i<A2plusb2.size(); i++){
	A2plusb2[i] = 0.0;
	for(int j=0; j<A2plus.ncols(); j++){
	    A2plusb2[i] += A2plus[i][j]*b2[j];
	}
    }
    // printf("A2plusb2\n");////////////////////
    // for(int i=0; i<A2plusb2.size(); i++){//////////////////
    // 	printf("A2plusb2[%d]=%g\n", i, A2plusb2[i]);
    // }

    m_t[2] = A2plusb2[0];
    m_alpha = A2plusb2[1];

    // printf("EuclidModel::calc(): t3=%g\n", m_t[2]);/////////////////
    // printf("EuclidModel::calc(): alpha=%g\n", m_alpha);/////////////////
}

/*****  getP  P[7] ***************/
void EuclidModel::getParams(VecDoub &P){
    P[0] = m_alpha;
    P[1] = m_theta_rad;
    P[2] = 0.0;
    P[3] = 0.0;
    P[4] = m_t[0];
    P[5] = m_t[1];
    P[6] = m_t[2];
}

/*****   ***************/
double EuclidModel::getGeometricError(){
    
    return m_kizyunTen.getGeometricError(m_P);
}

/***** calcCameraParams *****************************************************
 *  ノート23.4節    2025.10.16--
 *  posiveAngleIdxは正と予想される回転角のインデックス（ないときは負の値にする）
 *          *************************************************/
void EuclidModel::calcCameraParams(MatDoub_I &P, int positiveAngleIdx,
				   MatDoub_O &K, MatDoub_O &R, VecDoub_O &theta,
				   VecDoub_O &thetaPrime, VecDoub_O &t){
    MatDoub A(3, 3);
    VecDoub a(3);
    for(int i=0; i<3; i++){
	for(int j=0; j<3; j++){
	    A[i][j] = P[i][j];
	}
    }
    
    for(int i=0; i<3; i++) a[i] = P[i][3];

    LinAlg la;
    MatDoub B(3, 3);
    int c;
    la.RQdcmp(A, c, B, R, theta, thetaPrime);

    if(positiveAngleIdx>=0 && theta[positiveAngleIdx]<0.0 && thetaPrime[positiveAngleIdx]>=0.0){
	for(int i=0; i<t.size(); i++){
	    double x {theta[i]};
	    theta[i] = thetaPrime[i];
	    thetaPrime[i] = x;
	}
    }
    
    for(int i=0; i<K.nrows(); i++){
	for(int j=0; j<K.ncols(); j++){
	    K[i][j] = B[i][j]/B[2][2];
	}
    }

    LUdcmp lud(B);
    MatDoub Binv(3, 3);
    lud.inverse(Binv);
    
    for(int i=0; i<t.size(); i++){
	t[i] = 0.0;
	for(int j=0; j<a.size(); j++) t[i] += Binv[i][j]*a[j]/c;
    }
}

/********************************************************/
void EuclidModel::calcCameraParamsYXZ(MatDoub_I &P, int positiveAngleIdx){
    MatDoub K(3, 3);
    calcCameraParamsYXZ(P, positiveAngleIdx, K);
}

/********************************************************/
void EuclidModel::calcCameraParamsYXZ(MatDoub_I &P, int positiveAngleIdx, MatDoub_O &K){
    MatDoub R(3, 3);
    VecDoub theta(3), thetaPrime(3), t(3);

    calcCameraParamsYXZ(P, positiveAngleIdx, K, R, theta, thetaPrime, t);

    LinAlg la;
    la.printMatrix(K, "K");///////////////////
    la.printMatrix(R, "R");///////////////////
    la.printVector(theta, "theta [deg]", 180.0/GlobalConstants::pi);
    la.printVector(thetaPrime, "thetaPrime [deg]", 180.0/GlobalConstants::pi);//////
    la.printVector(t, "t");////////////////

    MatDoub RT(3, 3);
    la.transpose(R, RT);
    VecDoub RTt(3);
    la.prod(RT, t, RTt);
    printf("Camera center: X=%g, Y=%g, Z=%g\n", -RTt[0], -RTt[1], -RTt[2]);
}


/***** calcCameraParamsYZX *****************************************************
 *  2025.10.19  ノート23.4節終わりの方
 *  positiveAngleIdxはXYZの順
 *  出力の回転は、YXZの順 (theta[1]-->theta[0]-->theta[2])    *******************/
void EuclidModel::calcCameraParamsYXZ(MatDoub_I &P, int positiveAngleIdx, MatDoub_O &K,
				      MatDoub_O &R, VecDoub_O &theta, VecDoub_O &thetaPrime,
				      VecDoub_O &t){
    LinAlg la;

    MatDoub P1(3, 4);
    MatDoub Pprime(3, 4);
    la.exchangeRows(P, 0, 1, P1);
    la.exchangeCols(P1, 0, 1, Pprime);
    
    MatDoub A(3, 3);
    VecDoub a(3);
    for(int i=0; i<3; i++){
	for(int j=0; j<3; j++){
	    A[i][j] = Pprime[i][j];
	}
    }
    
    for(int i=0; i<3; i++) a[i] = Pprime[i][3];

    MatDoub B(3, 3);
    int c;
    la.RQdcmp(A, c, B, R, theta, thetaPrime);

    LUdcmp lud(B);
    MatDoub Binv(3, 3);
    lud.inverse(Binv);
    
    for(int i=0; i<t.size(); i++){
	t[i] = 0.0;
	for(int j=0; j<a.size(); j++) t[i] += Binv[i][j]*a[j]/c;
    }

    double x {t[0]};
    t[0] = t[1];
    t[1] = x;

    MatDoub B1(3, 3);
    la.exchangeRows(B, 0, 1, B1);
    la.exchangeCols(B1, 0, 1, B);
    la.scalarMul(1.0/B[2][2], B, K);

    x = theta[0];
    theta[0] = -theta[1];
    theta[1] = -x;
    theta[2] = -theta[2];
    x = thetaPrime[0];
    thetaPrime[0] = -thetaPrime[1];
    thetaPrime[1] = -x;
    thetaPrime[2] = -thetaPrime[2];

    MatDoub R1(3, 3);
    la.exchangeRows(R, 0, 1, R1);
    la.exchangeCols(R1, 0, 1, R);

    if(positiveAngleIdx>=0 && theta[positiveAngleIdx]<0.0 && thetaPrime[positiveAngleIdx]>=0.0){
	for(int i=0; i<t.size(); i++){
	    double x {theta[i]};
	    theta[i] = thetaPrime[i];
	    thetaPrime[i] = x;
	}
    }
}

/***** My note section 22.11 ***/
void EuclidModel::calc7paramsFromKizyunten(VecDoub_O &params, MatDoub_O &P){
    double X[6];
    double Y[6];
    double Z[6];
    int x[6];
    int y[6];

    X[0] = m_kizyunTen.m_Le[0]; Y[0] = m_kizyunTen.m_Le[1]; Z[0] = m_kizyunTen.m_Le[2]; 
    X[1] = m_kizyunTen.m_Lw[0]; Y[1] = m_kizyunTen.m_Lw[1]; Z[1] = m_kizyunTen.m_Lw[2]; 
    X[2] = m_kizyunTen.m_Me[0]; Y[2] = m_kizyunTen.m_Me[1]; Z[2] = m_kizyunTen.m_Me[2]; 
    X[3] = m_kizyunTen.m_Mw[0]; Y[3] = m_kizyunTen.m_Mw[1]; Z[3] = m_kizyunTen.m_Mw[2]; 
    X[4] = m_kizyunTen.m_He[0]; Y[4] = m_kizyunTen.m_He[1]; Z[4] = m_kizyunTen.m_He[2]; 
    X[5] = m_kizyunTen.m_Hw[0]; Y[5] = m_kizyunTen.m_Hw[1]; Z[5] = m_kizyunTen.m_Hw[2];
    
    x[0] = m_kizyunTen.m_le[0]; y[0] = m_kizyunTen.m_le[1];
    x[1] = m_kizyunTen.m_lw[0]; y[1] = m_kizyunTen.m_lw[1];
    x[2] = m_kizyunTen.m_me[0]; y[2] = m_kizyunTen.m_me[1];
    x[3] = m_kizyunTen.m_mw[0]; y[3] = m_kizyunTen.m_mw[1];
    x[4] = m_kizyunTen.m_he[0]; y[4] = m_kizyunTen.m_he[1];
    x[5] = m_kizyunTen.m_hw[0]; y[5] = m_kizyunTen.m_hw[1];

    for(int i=0; i<6; i++){/////////////////////////////
	printf("X[%d]=%f, Y[%d]=%f, Z[%d]=%f\n", i, X[i], i, Y[i], i, Z[i]);
    }
    for(int i=0; i<6; i++){/////////////////////////////
	printf("x[%d]=%d, y[%d]=%d\n", i, x[i], i, y[i]);
    }
    
    double q[5];
    for(int i=0; i<=4; i+=2){
	q[i] = (X[i]-X[i+1])/(x[i]-x[i+1]);
	// printf("q[%d]=%g\n", i, q[i]);///////////////////
    }
    for(int i=1; i<=5; i+=2){
	q[i] = q[i-1];
    }

    double qx[5], qy[5], qq[5];
    for(int i=2; i<=4; i+=2){
	qx[i] = q[i]*x[i]-q[0]*x[0];
	qy[i] = q[i]*y[i]-q[0]*y[0];
	qq[i] = q[i]-q[0];
    }
    for(int i=3; i<=5; i+=2){
	qx[i] = q[i]*x[i] - q[1]*x[1];
 	qy[i] = q[i]*y[i] - q[1]*y[1];
	qq[i] = q[i]-q[1];
    }

    double bunsi {0.0};
    double bunbo {0.0};
    for(int i=2; i<=5; i++){
	bunsi += qq[i]*qx[i];
	bunbo += qq[i]*qq[i];
    }
    
    double px;
    px = bunsi/bunbo;;
    printf("px=%g\n", px);/////////////////////
    
    double r[5];
    for(int i=2; i<=5; i++){
	r[i] = std::sqrt( SQR(Y[i]-Y[0])+SQR(Z[i]-Z[0]) );
	// printf("r[%d]=%g\n", i, r[i]);/////////////////////////
    }
    
    double theta_a {std::atan2(Y[4]-Y[0], Z[4]-Z[0])};
    printf("theta_a=%g (deg)\n", theta_a*180.0/GlobalConstants::pi);////////////////

    MatDoub A(4, 2);
    A[0][0] = qq[2]; A[0][1] = r[2];
    A[1][0] = qq[3]; A[1][1] = r[3];
    A[2][0] = qq[4]; A[2][1] = r[4];
    A[3][0] = qq[5]; A[3][1] = r[5];

    VecDoub b(4);
    for(int i=0; i<4; i++){
	b[i] = qy[i+2];
    }

    SVD svd(A);
    for(int i=0; i<svd.w.size(); i++){///////////////////////////
        printf("SVD: w[%d]=%10f\n", i, svd.w[i]);
    }

    VecDoub pys(2);
    svd.solve(b, pys, -1.0);

    // for(int i=0; i<2; i++) printf("pys[%d]=%g\n", i, pys[i]);//////////////////
    
    double py, theta_x;
    py = pys[0];
    theta_x = std::asin(pys[1]) - theta_a;
    printf("py=%g, theta_x=%g (deg)\n", py, theta_x*180.0/GlobalConstants::pi);

    bunsi = 0.0;
    bunbo = 0.0;
    for(int i=2; i<=5; i++){
	bunsi += r[i]*qq[i];
	bunbo += qq[i]*qq[i];
    }
    double alpha {bunsi/bunbo*std::cos(theta_x+theta_a)};
    printf("alpha=%g\n", alpha);//////////////////////////

    MatDoub R(3, 3, 0.0);
    R[0][0] = 1.0;
    R[1][1] = std::cos(theta_x);  R[1][2] = std::sin(theta_x);
    R[2][1] = -std::sin(theta_x); R[2][2] = std::cos(theta_x);

    VecDoub t(3, 0.0);
    VecDoub XYZ(3), RXYZ(3);
    LinAlg la;
    for(int i=0; i<=5; i++){
	XYZ[0] = X[i]; XYZ[1] = Y[i]; XYZ[2] = Z[i];
	la.prod(R, XYZ, RXYZ);
	t[0] += q[i]*(x[i]-px)-RXYZ[0];
	t[1] += q[i]*(y[i]-py)-RXYZ[1];
	t[2] += q[i]*alpha-RXYZ[2];
    }
    for(int i=0; i<3; i++) t[i] /= 6.0;
    la.printVector(t, "t");

    MatDoub K(3, 3, 0.0), Rt(3, 4);
    K[0][0] = alpha, K[0][2] = px;
    K[1][1] = alpha, K[1][2] = py;
    K[2][2] = 1.0;
    la.concat(R, t, Rt);
    la.prod(K, Rt, P);
    
    params[0] = px;
    params[1] = py;
    params[2] = alpha;
    params[3] = theta_x;
    params[4] = t[0];
    params[5] = t[1];
    params[6] = t[2];

    la.printMatrix(P, "P");/////////////////////t
}

/******** Calculate K, t from the normalized K', t' (2025.11.4)  *****************/
void EuclidModel::denormalizeKt(KizyunTen &kizyunTen, MatDoub_IO &K, MatDoub_I &R,
				VecDoub_IO &t){
    LinAlg la;

    K[0][0] = K[0][0]/kizyunTen.m_Ts;
    K[0][1] = K[0][1]/kizyunTen.m_Ts;
    K[0][2] = (K[0][2]-kizyunTen.m_Tt[0])/kizyunTen.m_Ts;
    K[1][0] = K[1][0]/kizyunTen.m_Ts;
    K[1][1] = K[1][1]/kizyunTen.m_Ts;
    K[1][2] = (K[1][2]-kizyunTen.m_Tt[1])/kizyunTen.m_Ts;

    VecDoub tU(3);
    for(int i=0; i<3; i++) tU[i] = kizyunTen.m_Ut[i];
    
    VecDoub Rt(3);
    la.prod(R, tU, Rt);

    for(int i=0; i<3; i++) Rt[i] += t[i];
    
    for(int i=0; i<3; i++) t[i] = Rt[i]/kizyunTen.m_Us;
}

/******** 2026.1.6,  My note section 28.7 *************************************/
void EuclidModel::calcSecondCameraP(MatDoub_I &K, MatDoub_O &P, std::string strRL){
    double X[6];
    double Y[6];
    double Z[6];
    int x[6];
    int y[6];

    X[0] = m_kizyunTen.m_Le[0]; Y[0] = m_kizyunTen.m_Le[1]; Z[0] = m_kizyunTen.m_Le[2]; 
    X[1] = m_kizyunTen.m_Lw[0]; Y[1] = m_kizyunTen.m_Lw[1]; Z[1] = m_kizyunTen.m_Lw[2]; 
    X[2] = m_kizyunTen.m_Me[0]; Y[2] = m_kizyunTen.m_Me[1]; Z[2] = m_kizyunTen.m_Me[2]; 
    X[3] = m_kizyunTen.m_Mw[0]; Y[3] = m_kizyunTen.m_Mw[1]; Z[3] = m_kizyunTen.m_Mw[2]; 
    X[4] = m_kizyunTen.m_He[0]; Y[4] = m_kizyunTen.m_He[1]; Z[4] = m_kizyunTen.m_He[2]; 
    X[5] = m_kizyunTen.m_Hw[0]; Y[5] = m_kizyunTen.m_Hw[1]; Z[5] = m_kizyunTen.m_Hw[2];
    
    x[0] = m_kizyunTen.m_le[0]; y[0] = m_kizyunTen.m_le[1];
    x[1] = m_kizyunTen.m_lw[0]; y[1] = m_kizyunTen.m_lw[1];
    x[2] = m_kizyunTen.m_me[0]; y[2] = m_kizyunTen.m_me[1];
    x[3] = m_kizyunTen.m_mw[0]; y[3] = m_kizyunTen.m_mw[1];
    x[4] = m_kizyunTen.m_he[0]; y[4] = m_kizyunTen.m_he[1];
    x[5] = m_kizyunTen.m_hw[0]; y[5] = m_kizyunTen.m_hw[1];

    for(int i=0; i<6; i++){/////////////////////////////
	printf("X[%d]=%f, Y[%d]=%f, Z[%d]=%f\n", i, X[i], i, Y[i], i, Z[i]);
    }
    for(int i=0; i<6; i++){/////////////////////////////
	printf("x[%d]=%d, y[%d]=%d\n", i, x[i], i, y[i]);
    }

    double alpha {(K[0][0]+K[1][1])*0.5};
    double px {K[0][2]};
    double py {K[1][2]};
    printf("alpha=%g, px=%g, py=%g\n", alpha, px, py);//////////////////////////
    double alphaInv {1.0/alpha};
    
    double a[6], b[6];
    for(int i=0; i<6; i++){
	a[i] = alphaInv*(x[i]-px);
	b[i] = alphaInv*(y[i]-py);
    }

    MatDoub A(9, 9, 0.0);
    double mx {X[1]-X[0]};
    printf("mx=%g\n", mx);////////////////////////

    A[0][0] = a[0]; A[0][1] = -a[1]; A[0][6] = mx;
    A[1][0] = b[0]; A[1][1] = -b[1]; A[1][7] = mx;
    A[2][0] = 1.0;  A[2][1] = -1.0;  A[2][8] = mx;
    A[3][2] = a[2]; A[3][3] = -a[3]; A[3][6] = mx;
    A[4][2] = b[2]; A[4][3] = -b[3]; A[4][7] = mx;
    A[5][2] = 1.0;  A[5][3] = -1.0;  A[5][8] = mx;
    A[6][4] = a[4]; A[6][5] = -a[5]; A[6][6] = mx;
    A[7][4] = b[4]; A[7][5] = -b[5]; A[7][7] = mx;
    A[8][4] = 1.0;  A[8][5] = -1.0;  A[8][8] = mx;
    
    SVD svd(A);
    // for(int i=0; i<9; i++) printf("w[%d]=%g\n", i, svd.w[i]);/////////////////////

    VecDoub xx(9);
    for(int i=0; i<9; i++) xx[i] = svd.v[i][8];
    // for(int i=0; i<9; i++) printf("xx[%d]=%g\n", i, xx[i]);///////////////////////

    double c;
    c = std::sqrt(xx[6]*xx[6] +xx[7]*xx[7]+xx[8]*xx[8]);

    double cY;

    if(xx[6]<0.0) c = -c;

    cY=xx[6]/c;

    double thetaY;
    if(strRL=="right") thetaY = std::acos(cY);
    else if(strRL=="left") thetaY = -std::acos(cY);
    
    printf("thetaY=%g (deg)\n", thetaY*180.0/GlobalConstants::pi);

    double sY {std::sin(thetaY)};

    double sX {xx[7]/(c*sY)};
    double cX {xx[8]/(c*sY)};
    
    double thetaX {std::atan2(sX, cX)};
    printf("thetaX=%g (deg)\n", thetaX*180.0/GlobalConstants::pi);

    double Zc[6];
    for(int i=0; i<=5; i++) Zc[i] = xx[i]/c;
    
    // for(int i=0; i<=5; i++) printf("Zc[%d]=%g\n", i, Zc[i]);//////////////////

    MatDoub RX(3, 3, 0.), RY(3, 3, 0.), R(3, 3);
    RX[0][0] = 1.0;
    RX[1][1] = cX;  RX[1][2] = sX;
    RX[2][1] = -sX; RX[2][2] = cX;
    RY[0][0] = cY; RY[0][2] = -sY;
    RY[1][1] = 1.0;
    RY[2][0] = sY; RY[2][2] = cY;

    LinAlg la;
    la.prod(RX, RY, R);
    
    VecDoub t(3, 0.0), XYZi(3), RXYZi(3);
    for(int i=0; i<6; i++){
	XYZi[0]=X[i]; XYZi[1]=Y[i]; XYZi[2]=Z[i];
	la.prod(R, XYZi, RXYZi);
	t[0] += Zc[i]*a[i]-RXYZi[0];
	t[1] += Zc[i]*b[i]-RXYZi[1];
	t[2] += Zc[i]-RXYZi[2];
    }
    for(int j=0; j<3; j++) t[j] /= 6.0;
    // la.printVector(t, "t");///////////////////////

    MatDoub Rt(3, 4);
    la.concat(R, t, Rt);

    la.prod(K, Rt, P);
    la.printMatrix(P, "secondCamera P");////////////////////
    
    MatDoub Rtrans(3, 3);
    la.transpose(R, Rtrans);
    VecDoub mCC(3);
    la.prod(Rtrans, t, mCC);
    la.printVector(mCC, "Camera Center", -1.0);
}

/******** 2026.2.3,  My note section 29.6 for paper *************************************/
void EuclidModel::calcSecondCameraPfromK(MatDoub_I &K, MatDoub_O &P, std::string strRL){
    double X[6];
    double Y[6];
    double Z[6];
    int x[6];
    int y[6];

    X[0] = m_kizyunTen.m_Le[0]; Y[0] = m_kizyunTen.m_Le[1]; Z[0] = m_kizyunTen.m_Le[2]; 
    X[1] = m_kizyunTen.m_Lw[0]; Y[1] = m_kizyunTen.m_Lw[1]; Z[1] = m_kizyunTen.m_Lw[2]; 
    X[2] = m_kizyunTen.m_Me[0]; Y[2] = m_kizyunTen.m_Me[1]; Z[2] = m_kizyunTen.m_Me[2]; 
    X[3] = m_kizyunTen.m_Mw[0]; Y[3] = m_kizyunTen.m_Mw[1]; Z[3] = m_kizyunTen.m_Mw[2]; 
    X[4] = m_kizyunTen.m_He[0]; Y[4] = m_kizyunTen.m_He[1]; Z[4] = m_kizyunTen.m_He[2]; 
    X[5] = m_kizyunTen.m_Hw[0]; Y[5] = m_kizyunTen.m_Hw[1]; Z[5] = m_kizyunTen.m_Hw[2];
    
    x[0] = m_kizyunTen.m_le[0]; y[0] = m_kizyunTen.m_le[1];
    x[1] = m_kizyunTen.m_lw[0]; y[1] = m_kizyunTen.m_lw[1];
    x[2] = m_kizyunTen.m_me[0]; y[2] = m_kizyunTen.m_me[1];
    x[3] = m_kizyunTen.m_mw[0]; y[3] = m_kizyunTen.m_mw[1];
    x[4] = m_kizyunTen.m_he[0]; y[4] = m_kizyunTen.m_he[1];
    x[5] = m_kizyunTen.m_hw[0]; y[5] = m_kizyunTen.m_hw[1];

    for(int i=0; i<6; i++){/////////////////////////////
	printf("X[%d]=%f, Y[%d]=%f, Z[%d]=%f\n", i, X[i], i, Y[i], i, Z[i]);
    }
    for(int i=0; i<6; i++){/////////////////////////////
	printf("x[%d]=%d, y[%d]=%d\n", i, x[i], i, y[i]);
    }

    LUdcmp ludcmp(K);
    MatDoub Kinv(3, 3);
    ludcmp.inverse(Kinv);
    
    LinAlg la;
    // la.printMatrix(Kinv, "Kinv");////////////////////////////
    
    double a[6][3];
    for(int i=0; i<6; i++){
	double xx[3], yy[3];
	xx[0] = x[i]; xx[1] = y[i]; xx[2] = 1.0;
	la.prod(Kinv, xx, yy);
	for(int j=0; j<3; j++) a[i][j] = yy[j];
    }

    MatDoub A(9, 9, 0.0);

    A[0][0] = a[0][0]; A[0][1] = -a[1][0]; A[0][6] = X[1]-X[0];
    A[1][0] = a[0][1]; A[1][1] = -a[1][1]; A[1][7] = X[1]-X[0];
    A[2][0] = a[0][2]; A[2][1] = -a[1][2]; A[2][8] = X[1]-X[0];
    A[3][2] = a[2][0]; A[3][3] = -a[3][0]; A[3][6] = X[3]-X[2];
    A[4][2] = a[2][1]; A[4][3] = -a[3][1]; A[4][7] = X[3]-X[2];
    A[5][2] = a[2][2]; A[5][3] = -a[3][2]; A[5][8] = X[3]-X[2];
    A[6][4] = a[4][0]; A[6][5] = -a[5][0]; A[6][6] = X[5]-X[4];
    A[7][4] = a[4][1]; A[7][5] = -a[5][1]; A[7][7] = X[5]-X[4];
    A[8][4] = a[4][2]; A[8][5] = -a[5][2]; A[8][8] = X[5]-X[4];
    
    SVD svd(A);
    // for(int i=0; i<9; i++) printf("w[%d]=%g\n", i, svd.w[i]);/////////////////////

    VecDoub xx(9);
    for(int i=0; i<9; i++) xx[i] = svd.v[i][8];
    // for(int i=0; i<9; i++) printf("xx[%d]=%g\n", i, xx[i]);///////////////////////

    double b;
    b = std::sqrt(xx[6]*xx[6] +xx[7]*xx[7]+xx[8]*xx[8]);

    double k[6], c1, s0s1, c0s1;
    if(xx[6]>=0.0){
	for(int i=0; i<6; i++) k[i] = xx[i]/b;
	c1 = xx[6]/b;
	s0s1 = xx[7]/b;
	c0s1 = xx[8]/b;
    }else{
	for(int i=0; i<6; i++) k[i] = -xx[i]/b;
	c1 = -xx[6]/b;
	s0s1 = -xx[7]/b;
	c0s1 = -xx[8]/b;
    }
	
    double theta1;
    if(strRL=="right") theta1 = std::acos(c1);
    else if(strRL=="left") theta1 = -std::acos(c1);
    
    printf("theta1=%g (deg)\n", theta1*180.0/GlobalConstants::pi);////////////////

    double s1 {std::sin(theta1)};

    double theta0 {std::atan2(s0s1/s1, c0s1/s1)};
    printf("theta0=%g (deg)\n", theta0*180.0/GlobalConstants::pi);

    double c0 {std::cos(theta0)};
    double s0 {std::sin(theta0)};

    MatDoub R(3, 3);
    R[0][0] = c1; R[0][1] = 0.0; R[0][2] = -s1;
    R[1][0] = s0*s1; R[1][1] = c0; R[1][2] = s0*c1;
    R[2][0] = c0*s1; R[2][1] = -s0; R[2][2] = c0*c1; 

    VecDoub t(3, 0.0), kiai(3), Xi(3), RXi(3);
    for(int i=0; i<6; i++){
	Xi[0]=X[i]; Xi[1]=Y[i]; Xi[2]=Z[i];
	la.prod(R, Xi, RXi);
	t[0] += k[i]*a[i][0]-RXi[0];
	t[1] += k[i]*a[i][1]-RXi[1];
	t[2] += k[i]*a[i][2]-RXi[2];
    }
    for(int j=0; j<3; j++) t[j] /= 6.0;
    la.printVector(t, "t");///////////////////////

    MatDoub Rt(3, 4);
    la.concat(R, t, Rt);

    la.prod(K, Rt, P);
    la.printMatrix(P, "secondCamera P");////////////////////

    MatDoub Rtrans(3, 3);
    la.transpose(R, Rtrans);
    VecDoub mCC(3);
    la.prod(Rtrans, t, mCC);
    la.printVector(mCC, "Camera Center", -1.0);
}

/******** 2026.2.4,  My note 2602041602-
    params[0] = px;
    params[1] = py;
    params[2] = alpha;
*************************************/
void EuclidModel::calcSecondCameraP(VecDoub_I &params, MatDoub_O &P, std::string strRL){
    double X[6];
    double Y[6];
    double Z[6];
    int x[6];
    int y[6];

    X[0] = m_kizyunTen.m_Le[0]; Y[0] = m_kizyunTen.m_Le[1]; Z[0] = m_kizyunTen.m_Le[2]; 
    X[1] = m_kizyunTen.m_Lw[0]; Y[1] = m_kizyunTen.m_Lw[1]; Z[1] = m_kizyunTen.m_Lw[2]; 
    X[2] = m_kizyunTen.m_Me[0]; Y[2] = m_kizyunTen.m_Me[1]; Z[2] = m_kizyunTen.m_Me[2]; 
    X[3] = m_kizyunTen.m_Mw[0]; Y[3] = m_kizyunTen.m_Mw[1]; Z[3] = m_kizyunTen.m_Mw[2]; 
    X[4] = m_kizyunTen.m_He[0]; Y[4] = m_kizyunTen.m_He[1]; Z[4] = m_kizyunTen.m_He[2]; 
    X[5] = m_kizyunTen.m_Hw[0]; Y[5] = m_kizyunTen.m_Hw[1]; Z[5] = m_kizyunTen.m_Hw[2];
    
    x[0] = m_kizyunTen.m_le[0]; y[0] = m_kizyunTen.m_le[1];
    x[1] = m_kizyunTen.m_lw[0]; y[1] = m_kizyunTen.m_lw[1];
    x[2] = m_kizyunTen.m_me[0]; y[2] = m_kizyunTen.m_me[1];
    x[3] = m_kizyunTen.m_mw[0]; y[3] = m_kizyunTen.m_mw[1];
    x[4] = m_kizyunTen.m_he[0]; y[4] = m_kizyunTen.m_he[1];
    x[5] = m_kizyunTen.m_hw[0]; y[5] = m_kizyunTen.m_hw[1];

    for(int i=0; i<6; i++){/////////////////////////////
	printf("X[%d]=%f, Y[%d]=%f, Z[%d]=%f\n", i, X[i], i, Y[i], i, Z[i]);
    }
    for(int i=0; i<6; i++){/////////////////////////////
	printf("x[%d]=%d, y[%d]=%d\n", i, x[i], i, y[i]);
    }

    double alpha {params[2]};
    double px {params[0]};
    double py {params[1]};
    printf("alpha=%g, px=%g, py=%g\n", alpha, px, py);//////////////////////////
    
    double a[6], b[6];
    for(int i=0; i<6; i++){
	a[i] = x[i]-px;
	b[i] = y[i]-py;
    }

    MatDoub A(9, 9, 0.0);
    double mx {X[1]-X[0]};
    printf("mx=%g\n", mx);////////////////////////

    A[0][0] = a[0]; A[0][1] = -a[1]; A[0][6] = mx;
    A[1][0] = b[0]; A[1][1] = -b[1]; A[1][7] = mx;
    A[2][0] = alpha;A[2][1] = -alpha;A[2][8] = mx;
    A[3][2] = a[2]; A[3][3] = -a[3]; A[3][6] = mx;
    A[4][2] = b[2]; A[4][3] = -b[3]; A[4][7] = mx;
    A[5][2] = alpha;A[5][3] = -alpha;A[5][8] = mx;
    A[6][4] = a[4]; A[6][5] = -a[5]; A[6][6] = mx;
    A[7][4] = b[4]; A[7][5] = -b[5]; A[7][7] = mx;
    A[8][4] = alpha;A[8][5] = -alpha;A[8][8] = mx;
    
    SVD svd(A);
    // for(int i=0; i<9; i++) printf("w[%d]=%g\n", i, svd.w[i]);/////////////////////

    VecDoub xx(9);
    for(int i=0; i<9; i++) xx[i] = svd.v[i][8];
    // for(int i=0; i<9; i++) printf("xx[%d]=%g\n", i, xx[i]);///////////////////////

    double W;
    W = std::sqrt(xx[6]*xx[6] +xx[7]*xx[7]+xx[8]*xx[8]);

    if(xx[6]<0.0) W = -W;

    double q[6], c1, s0s1, c0s1;
    for(int i=0; i<6; i++) q[i] = xx[i]/W;
    c1 = xx[6]/W;
    s0s1 = xx[7]/W;
    c0s1 = xx[8]/W;
	
    double theta1;
    if(strRL=="right") theta1 = std::acos(c1);
    else if(strRL=="left") theta1 = -std::acos(c1);
    
    printf("theta1=%g (deg)\n", theta1*180.0/GlobalConstants::pi);////////////////

    double s1 {std::sin(theta1)};

    double theta0 {std::atan2(s0s1/s1, c0s1/s1)};
    printf("theta0=%g (deg)\n", theta0*180.0/GlobalConstants::pi);

    double c0 {std::cos(theta0)};
    double s0 {std::sin(theta0)};

    MatDoub R(3, 3);
    R[0][0] = c1; R[0][1] = 0.0; R[0][2] = -s1;
    R[1][0] = s0*s1; R[1][1] = c0; R[1][2] = s0*c1;
    R[2][0] = c0*s1; R[2][1] = -s0; R[2][2] = c0*c1; 

    LinAlg la;
    
    VecDoub t(3, 0.0), qiai(3), Xi(3), RXi(3);
    for(int i=0; i<6; i++){
	Xi[0]=X[i]; Xi[1]=Y[i]; Xi[2]=Z[i];
	la.prod(R, Xi, RXi);
	t[0] += q[i]*a[i]-RXi[0];
	t[1] += q[i]*b[i]-RXi[1];
	t[2] += q[i]*alpha-RXi[2];
    }
    for(int j=0; j<3; j++) t[j] /= 6.0;
    la.printVector(t, "t");///////////////////////

    MatDoub Rt(3, 4);
    la.concat(R, t, Rt);

    MatDoub K(3, 3, 0.0);
    K[0][0] = alpha; K[0][2] = px;
    K[1][1] = alpha; K[1][2] = py;
    K[2][2] = 1.0;

    la.prod(K, Rt, P);
    la.printMatrix(P, "secondCamera P");////////////////////

    MatDoub Rtrans(3, 3);
    la.transpose(R, Rtrans);
    VecDoub mCC(3);
    la.prod(Rtrans, t, mCC);
    la.printVector(mCC, "Camera Center", -1.0);
}

/******** 2026.3.12--,  My note, Section 29.6, 2603111752-
    params[0] = px;
    params[1] = py;
    params[2] = alpha;
    Calculate P from 4 ref. points
*************************************/
void EuclidModel::calcSecondCameraPHomography(KizyunTen kt, MatDoub_I &K, MatDoub_O &P){
    LinAlg la;
    GlobalConstants gc;
    double X[5], Yprime[5], Zprime[5];

    X[1] = 0.0; Yprime[1] = 0.0; Zprime[1] = 0.0;
    X[2] = gc.mX; Yprime[2] = 0.0; Zprime[2] = 0.0;
    X[3] = 0.0; Yprime[3] = gc.mY*4.0+gc.deltaY*3.0; Zprime[3] = 0.0;
    X[4] = gc.mX; Yprime[4] = gc.mY*4.0+gc.deltaY*3.0; Zprime[4] = 0.0;

    // for(int i=1; i<=4; i++){///////////////////////////////
    // 	printf("i=%d, X=%g, Yprime=%g, Zprime=%g\n", i, X[i], Yprime[i], Zprime[i]);
    // }

    double x[5], y[5];
    x[1] = kt.m_le[0]; y[1] = kt.m_le[1];
    x[2] = kt.m_lw[0]; y[2] = kt.m_lw[1];
    x[3] = kt.m_he[0]; y[3] = kt.m_he[1];
    x[4] = kt.m_hw[0]; y[4] = kt.m_hw[1];

    // for(int i=1; i<=4; i++){///////////////////////////////
    // 	printf("i=%d, x=%g, y=%g\n", i, x[i], y[i]);
    // }

    VecDoub u(5), v(5);
    MatDoub T(3, 3);
    kt.normalize(x, y, u, v, T);

    // la.printVector(u, "u");/////////////////////
    // la.printVector(v, "v");
    // la.printMatrix(T, "T");

    MatDoub Tinv(3, 3);
    la.calcInv(T, Tinv);
    
    // for(int i=1; i<=4; i++){//////////////////////////////
    // 	VecDoub ui(3), xi(3);
    // 	ui[0] = u[i]; ui[1] = v[i]; ui[2] = 1.0;
    // 	la.prod(Tinv, ui, xi);
    // 	printf("i=%d\n", i);
    // 	la.printVector(xi, "xi");
    // }
    
    VecDoub U(5), V(5);
    MatDoub S(3, 3);
    kt.normalize(X, Yprime, U, V, S);
    
    // la.printVector(U, "U");/////////////////////
    // la.printVector(V, "V");
    // la.printMatrix(S, "S");

    MatDoub A(8, 8, 0.0);
    VecDoub b(8);
    for(int i=1; i<=4; i++){
	int r {2*(i-1)+1};
	A[r-1][0]=U[i]; A[r-1][1]=V[i]; A[r-1][2]=1.0; A[r-1][6]=-u[i]*U[i]; A[r-1][7]=-u[i]*V[i]; 
	A[r][3]=U[i]; A[r][4]=V[i]; A[r][5]=1.0; A[r][6]=-v[i]*U[i]; A[r][7]=-v[i]*V[i];
	b[r-1] = u[i];
	b[r] = v[i];
    }

    // la.printMatrix(A, "A");////////////////////////
    // la.printVector(b, "b");

    VecDoub w(8);
    LUdcmp lud(A);
    lud.solve(b, w);
    // la.printVector(w, "w");///////////////////////

    MatDoub H(3, 3);
    H[0][0] = w[0]; H[0][1] = w[1]; H[0][2] = w[2]; 
    H[1][0] = w[3]; H[1][1] = w[4]; H[1][2] = w[5]; 
    H[2][0] = w[6]; H[2][1] = w[7]; H[2][2] = 1.0;
    // la.printMatrix(H, "H");//////////////////////

    MatDoub Hprime(3, 3), Xtmp(3, 3);
    la.prod(Tinv, H, Xtmp);
    la.prod(Xtmp, S, Hprime);
    // la.printMatrix(Hprime, "Hprime");/////////////////////

    MatDoub Kinv(3, 3), B(3, 3);
    la.calcInv(K, Kinv);
    la.prod(Kinv, Hprime, B);
    // la.printMatrix(B, "B");////////////////////////

    VecDoub b1(3), b2(3);
    b1[0] = B[0][0]; b2[0] = B[0][1];
    b1[1] = B[1][0]; b2[1] = B[1][1]; 
    b1[2] = B[2][0]; b2[2] = B[2][1];

    VecDoub e1(3), b2prime(3), e2(3), xTmp(3);
    la.scalarMul(1.0/la.getNorm(b1), b1, e1);
    // la.printVector(e1, "e1");//////////////////
    double b2e1 {la.getScalarProd(b2, e1)};
    la.scalarMul(-b2e1, e1, xTmp);
    la.add(b2, xTmp, b2prime);
    la.scalarMul(1.0/la.getNorm(b2prime), b2prime, e2);
    // la.printVector(e2, "e2");//////////////////

    double x1 {la.getNorm(b1)};
    double x2 {la.getScalarProd(b2,e1)};
    double y2 {la.getScalarProd(b2,e2)};
    // printf("x2=%g, y2=%g\n", x2, y2);///////////////////////
    double bunbo {std::sqrt((x1+y2)*(x1+y2)+x2*x2)};
    double xx {(x1+y2)/bunbo};
    double yy {x2/bunbo};

    VecDoub r1(3), r2(3), r3(3), xTmp2(3);
    la.scalarMul(xx, e1, xTmp);
    la.scalarMul(yy, e2, xTmp2);
    la.add(xTmp, xTmp2, r1);
    la.scalarMul(-yy, e1, xTmp);
    la.scalarMul(xx, e2, xTmp2);
    la.add(xTmp, xTmp2, r2);
    // la.printVector(r1, "r1");//////////////////////
    // la.printVector(r2, "r2");//////////////////////
    
    la.calcVectorProd(r1, r2, r3);
    la.printVector(r3, "r3");//////////////////////

    MatDoub R(3, 3);
    R[0][0] = r1[0]; R[0][1] = r2[0]; R[0][2] = r3[0]; 
    R[1][0] = r1[1]; R[1][1] = r2[1]; R[1][2] = r3[1]; 
    R[2][0] = r1[2]; R[2][1] = r2[2]; R[2][2] = r3[2];
    // la.printMatrix(R, "R");/////////////////////

    VecDoub xTilder(3);
    MatDoub a(5, 4);
    for(int i=1; i<=4; i++){
	xTilder[0] = x[i];
	xTilder[1] = y[i];
	xTilder[2] = 1.0;

	la.prod(Kinv, xTilder, xTmp);
	for(int j=1; j<=3; j++) a[i][j] = xTmp[j-1];
    }
    // la.printMatrix(a, "a");//////////////////////
    MatDoub c(5, 4);
    VecDoub mai(3), xTmp3(3), ci(3), bb(12);
    for(int i=1; i<=4; i++){
	for(int j=1; j<=3; j++) mai[j-1] = -a[i][j];
	la.scalarMul(X[i], r1, xTmp);
	la.scalarMul(Yprime[i], r2, xTmp2);
	la.add(xTmp, xTmp2, xTmp3);
	la.calcVectorProd(mai, xTmp3, ci);
	int i0 {3*(i-1)};
	for(int j=0; j<=2; j++) bb[i0+j] = ci[j];
    }
    // la.printVector(bb, "bb");//////////////////////////

    MatDoub AA(12, 3);
    for(int i=0; i<4; i++){
	int i0 {3*i};
	AA[i0][0] = 0.0; AA[i0][1] = -a[i+1][3]; AA[i0][2] = a[i+1][2]; 
	AA[i0+1][0] = a[i+1][3]; AA[i0+1][1] = 0.0; AA[i0+1][2] = -a[i+1][1]; 
	AA[i0+2][0] = -a[i+1][2]; AA[i0+2][1] = a[i+1][1]; AA[i0+2][2] = 0.0;
    }
    // la.printMatrix(AA, "AA");////////////////////////////

    VecDoub t(3), ATb(3);
    MatDoub AT(3, 12), ATA(3, 3), ATAinv(3, 3);
    la.transpose(AA, AT);
    la.prod(AT, AA, ATA);
    la.calcInv(ATA, ATAinv);
    la.prod(AT, bb, ATb);
    la.prod(ATAinv, ATb, t);
    // la.printVector(t, "t");/////////////////////

    MatDoub Q(3, 3, 0.0);
    double co {std::cos(gc.pi*80.0/180.0)};
    double si {std::sin(gc.pi*80.0/180.0)};
    Q[0][0] = 1.0;
    Q[1][1] = co; Q[1][2] = si;
    Q[2][1] = -si; Q[2][2] = co;

    MatDoub RQ(3, 3);
    VecDoub t0(3, 0.0);
    t0[2] = gc.hL;
    la.prod(R, Q, RQ);
    la.prod(RQ, t0, xTmp);
    la.subtract(t, xTmp, xTmp2);
    MatDoub RQt(3, 4);
    la.concat(RQ, xTmp2, RQt);

    la.prod(K, RQt, P);
    la.printMatrix(P, "P");////////////////////
}

