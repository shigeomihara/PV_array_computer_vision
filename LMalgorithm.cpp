/*********************************************************************
 *      LMalgorithm.cpp
 ********************************************************************/
// #include "main.hpp"
// #include "LMalgorithm.hpp"
// #include "ludcmp.h"

LMalgorithm::LMalgorithm(VecDoub_I &X, VecDoub_I &P0, KizyunTen &kt)
    :m_X(X), m_P0(P0), m_kt(kt){
}

LMalgorithm::LMalgorithm(VecDoub_I &X, VecDoub_I &P0, EuclidModel &em)
    :m_X(X), m_P0(P0), m_kt(em.getKizyunTen()), m_em(em){
}

// LMalgorithm::LMalgorithm(Homography &hm)
//     :m_hom(hm){
// }

/**** P0: not normalized camera matrix *****/
LMalgorithm::LMalgorithm(MatDoub_I &P0, KizyunTen &kt)
    :m_X(kt.get_n()*2), m_P0(P0.nrows()*P0.ncols()), m_kt(kt){

    LinAlg la;
    kt.get_xTilder(m_X);
    // la.printVector(m_X, "m_X");//////////////////

    MatDoub P0_nm(3, 4);
    kt.normalize_P(P0, P0_nm);
    // la.printMatrix(P0_nm, "P0_nm");//////////////////

    la.flatten(P0_nm, m_P0);
    // la.printVector(m_P0, "m_P0");////////////////////////
}


// /******************************************************/
// void LMalgorithm::calcHomography(){
//     // printf("hello calcHomography\n");//////////////////////
//     // for(int i=0; i<m_hom.m_H.nrows(); i++){////////////////////
//     // 	for(int j=0; j<m_hom.m_H.ncols(); j++){
//     // 	    printf("H[%d][%d]=%8.5f ", i, j, m_hom.m_H[i][j]);
//     // 	}
//     // 	printf("\n");
//     // }
//     // for(int i=0; i<m_hom.m_uHat.nrows(); i++){////////////////////
//     // 	for(int j=0; j<m_hom.m_uHat.ncols(); j++){
//     // 	    printf("uHat[%d][%d]=%8.5f ", i, j, m_hom.m_uHat[i][j]);
//     // 	}
//     // 	printf("\n");
//     // }

//     VecDoub fval(m_hom.m_n*4);
//     // m_hom.get_f(fval);
//     // for(int i=0; i<fval.size(); i++){/////////////////////////
//     // 	printf("fval[%d]=%g\n", i, fval[i]);
//     // }

//     MatDoub Jval(m_hom.m_n*4, m_hom.m_n*2+9);
//     // MatDoub Jvalnum(m_hom.m_n*4, m_hom.m_n*2+9);
//     // m_hom.getJ(Jval);
//     // m_hom.getJnum(Jvalnum);
//     // for(int i=0; i<Jval.nrows(); i++){/////////////////////
//     // 	for(int j=0; j<Jval.ncols(); j++){
//     // 	    printf("[%d][%d]=%g ", i, j, Jval[i][j]-Jvalnum[i][j]);
//     // 	}
//     // 	printf("\n");
//     // }

    
//     // 初期値設定
//     double lambda {0.0};
//     m_hom.getJ(Jval);
//     for(int j=0; j<Jval.ncols(); j++){
// 	double aii {0.0};
// 	for(int k=0; k<Jval.nrows(); k++) aii += Jval[k][j]*Jval[k][j];
// 	lambda += aii;
//     }
//     lambda /= Jval.ncols();
//     lambda *= 1.0E-3;
//     // printf("initial lambda=%f\n", lambda);/////////////////

//     int Fdim {m_hom.m_n*4};
//     VecDoub epsilon(Fdim);
//     VecDoub Y(Fdim);  // Objective vector
//     m_hom.getY(Y);
//     // for(int i=0; i<m_hom.m_n*4; i++){//////////////////////
//     // 	printf("Y[%d]=%g\n", i, Y[i]);
//     // }

//     double lambdaMax {std::numeric_limits<double>::max()*1.0E-10};
    
//     while(true){
// 	m_hom.get_f(fval);
	
// 	double Ecurr {0.0};
// 	for(int i=0; i<Fdim; i++){
// 	    epsilon[i] = fval[i]-Y[i];
// 	    Ecurr += epsilon[i]*epsilon[i];
// 	}
// 	// printf("Ecurr=%.10f\n", Ecurr);///////////////////
	
// 	m_hom.getJ(Jval);

// 	MatDoub a(Jval.ncols(), Jval.ncols());
// 	VecDoub Delta(Jval.ncols());
// 	VecDoub b(Jval.ncols());
	
// 	while(true){
// 	    for(int i=0; i<Jval.ncols(); i++){
// 		for(int j=0; j<Jval.ncols(); j++){
// 		    double sum {0.0};
// 		    for(int k=0; k<Jval.nrows(); k++) sum += Jval[k][i]*Jval[k][j];
// 		    a[i][j] = sum;
// 		}
// 	    }
// 	    for(int i=0; i<a.nrows(); i++) a[i][i] += lambda;
	
// 	    for(int i=0; i<Jval.ncols(); i++){
// 		double sum {0.0};
// 		for(int k=0; k<Fdim; k++) sum += Jval[k][i]*epsilon[k];
// 		b[i] = -sum;
// 	    }

// 	    LUdcmp lud(a);

// 	    lud.solve(b, Delta);
	
// 	    //for(int i=0; i<vecAlphaDim; i++) vecAlphaNew[i] = vecAlpha[i]+Delta[i];
// 	    for(int i=0; i<m_hom.m_n; i++){
// 		m_hom.m_uHat[i][0] += Delta[i*2];
// 		m_hom.m_uHat[i][1] += Delta[i*2+1];
// 	    }
// 	    for(int i=0; i<3; i++){
// 		for(int j=0; j<3; j++){
// 		    m_hom.m_H[i][j] += Delta[m_hom.m_n*2+i*3+j];
// 		}
// 	    }
		
// 	    m_hom.get_f(fval);
	    
// 	    double Enew {0.0};
// 	    for(int i=0; i<Fdim; i++){
// 		epsilon[i] = fval[i]-Y[i];
// 		Enew += epsilon[i]*epsilon[i];
// 	    }
// 	    // printf("Enew= %.10f, lambda=%g\n", Enew, lambda);///////////////////

// 	    if(Enew < Ecurr){
// 		lambda *= 0.1;
// 		break;
// 	    }else{
// 		lambda *= 10.0;
// 		if(lambda > lambdaMax){  // added !!!
// 		    break;
// 		}
// 		for(int i=0; i<m_hom.m_n; i++){
// 		    m_hom.m_uHat[i][0] -= Delta[i*2];
// 		    m_hom.m_uHat[i][1] -= Delta[i*2+1];
// 		}
// 		for(int i=0; i<3; i++){
// 		    for(int j=0; j<3; j++){
// 			m_hom.m_H[i][j] -= Delta[m_hom.m_n*2+i*3+j];
// 		    }
// 		}
// 	    }
// 	}

// 	double norX2 {0.0};
// 	double norDelta2 {0.0};
	
// 	for(int i=0; i<m_hom.m_n; i++){
// 	    norX2 += m_hom.m_uHat[i][0]*m_hom.m_uHat[i][0];
// 	    norX2 += m_hom.m_uHat[i][1]*m_hom.m_uHat[i][1];
// 	}
// 	for(int i=0; i<3; i++){
// 	    for(int j=0; j<3; j++){
// 		norX2 += m_hom.m_H[i][j]*m_hom.m_H[i][j];
// 	    }
// 	}
	
// 	for(int i=0; i<Delta.size(); i++){
// 	    norDelta2 += Delta[i]*Delta[i];
// 	}

// 	// Finish condition
// 	// printf("E=%g, norDelta2=%g, norX2=%g, lambda=%g\n",
// 	//        Ecurr, norDelta2, norX2, lambda);//////////t
// 	// if(norDelta2/norX2 < 1.0E-10) break;
// 	if(norDelta2/norX2 < 1.0E-16 || lambda > lambdaMax) break;
//     }

// }

/*****************************************************/
void LMalgorithm::calcEuclidModel(VecDoub_O &vecAlpha){
    // VecDoub F(m_kt.get_n()*2);
    // m_em.get_F(vecAlpha, F);

    // MatDoub J(m_kt.get_n()*2, vecAlpha.size());
    // m_em.getJ(vecAlpha, J);
    // m_em.printMatrix(J, "J");

    // MatDoub Jnum(m_kt.get_n()*2, vecAlpha.size());
    // m_em.getJnum(vecAlpha, Jnum);
    // m_em.printMatrix(Jnum, "Jnum");

    // printf("J-Jnum\n");/////////////////////
    // for(int i=0; i<J.nrows(); i++){/////////////////////
    // 	for(int j=0; j<J.ncols(); j++){
    // 	    printf("[%d][%d]=%.5f ", i, j, J[i][j]-Jnum[i][j]);
    // 	}
    // 	printf("\n");
    // }
    
    int Fdim {m_kt.get_n()*2};
    int vecAlphaDim {vecAlpha.size()};
    VecDoub F(Fdim);
    MatDoub J(Fdim, vecAlphaDim);
    double lambda {0.0};
    VecDoub epsilon(Fdim);
    MatDoub a(vecAlphaDim, vecAlphaDim);
    VecDoub b(vecAlphaDim);
    VecDoub Delta(vecAlphaDim);
    VecDoub vecAlphaNew(vecAlphaDim);

    // 初期値設定
    vecAlpha = m_P0;
    m_em.getJ(vecAlpha, J);
    for(int j=0; j<vecAlphaDim; j++){
	double aii {0.0};
	for(int k=0; k<J.nrows(); k++) aii += J[k][j]*J[k][j];
	lambda += aii;
    }
    lambda /= vecAlphaDim;
    lambda *= 1.0E-3;
    // printf("initial lambda=%f\n", lambda);/////////////////
    
    while(true){
	m_em.get_F(vecAlpha, F);
	
	double Ecurr {0.0};
	for(int i=0; i<Fdim; i++){
	    epsilon[i] = F[i]-m_X[i];
	    Ecurr += epsilon[i]*epsilon[i];
	}
	// printf("Ecurr=%.10f\n", Ecurr);///////////////////

	m_em.getJ(vecAlpha, J);

	while(true){
	    for(int i=0; i<vecAlphaDim; i++){
		for(int j=0; j<vecAlphaDim; j++){
		    double sum {0.0};
		    for(int k=0; k<Fdim; k++) sum += J[k][i]*J[k][j];
		    a[i][j] = sum;
		}
	    }
	    for(int i=0; i<vecAlphaDim; i++) a[i][i] += lambda;

	    for(int i=0; i<vecAlphaDim; i++){
		double sum {0.0};
		for(int k=0; k<Fdim; k++) sum += J[k][i]*epsilon[k];
		b[i] = -sum;
	    }

	    LUdcmp *lud = new LUdcmp(a);

	    lud->solve(b, Delta);
	    
	    for(int i=0; i<vecAlphaDim; i++) vecAlphaNew[i] = vecAlpha[i]+Delta[i];
	    
	    m_em.get_F(vecAlphaNew, F);
	    
	    double Enew {0.0};
	    for(int i=0; i<Fdim; i++){
		epsilon[i] = F[i]-m_X[i];
		Enew += epsilon[i]*epsilon[i];
	    }
	    // printf("Enew= %.10f, lambda=%g\n", Enew, lambda);///////////////////
	    if(Enew < Ecurr){
		lambda *= 0.1;
		break;
	    }else{
		lambda *= 10.0;
	    }

	    delete lud;
	}

	vecAlpha = vecAlphaNew;

	double norvecAlpha2 {0.0};
	double norDelta2 {0.0};
	for(int i=0; i<vecAlphaDim; i++){
	    norvecAlpha2 += vecAlpha[i]*vecAlpha[i];
	    norDelta2 += Delta[i]*Delta[i];
	}

	// Finish condition
	// if(norDelta2/norvecAlpha2 < 1.0E-10) break;
	if(norDelta2/norvecAlpha2 < 1.0E-16) break;
    }
}

/********* Not used. linear ******************************/
void LMalgorithm::calcEuclidModel0(VecDoub_O &P){
    // printf("calcEuclidModel\n");///////////////
    int fdim {m_kt.get_n()*2};
    int Pdim {P.size()};
    VecDoub f(fdim);
    MatDoub J(fdim, Pdim);
    double lambda {0.0};
    // VecDoub epsilon(fdim);
    MatDoub a(Pdim, Pdim);
    VecDoub b(Pdim);
    VecDoub Delta(Pdim);
    VecDoub Pnew(Pdim);

    // 初期値設定
    P = m_P0;
    m_em.getJ0(P, J);
    for(int j=0; j<Pdim; j++){
	double aii {0.0};
	for(int k=0; k<J.nrows(); k++) aii += J[k][j]*J[k][j];
	lambda += aii;
    }
    lambda /= Pdim;
    lambda *= 1.0E-3;
    // printf("initial lambda=%f\n", lambda);/////////////////
    
    while(true){
	m_em.get_f0(P, f);
	
	double Ecurr {0.0};
	for(int i=0; i<fdim; i++){
	    // epsilon[i] = f[i];
	    // Ecurr += epsilon[i]*epsilon[i];
	    Ecurr += f[i]*f[i];
	}
	// printf("Ecurr=%.10f\n", Ecurr);///////////////////

	m_em.getJ0(P, J);

	while(true){
	    for(int i=0; i<Pdim; i++){
		for(int j=0; j<Pdim; j++){
		    double sum {0.0};
		    for(int k=0; k<fdim; k++) sum += J[k][i]*J[k][j];
		    a[i][j] = sum;
		}
	    }
	    for(int i=0; i<Pdim; i++) a[i][i] += lambda;

	    for(int i=0; i<Pdim; i++){
		double sum {0.0};
		// for(int k=0; k<fdim; k++) sum += J[k][i]*epsilon[k];
		for(int k=0; k<fdim; k++) sum += J[k][i]*f[k];
		b[i] = -sum;
	    }

	    LUdcmp *lud = new LUdcmp(a);

	    lud->solve(b, Delta);
	    
	    for(int i=0; i<Pdim; i++) Pnew[i] = P[i]+Delta[i];
	    
	    m_em.get_f0(Pnew, f);
	    
	    double Enew {0.0};
	    for(int i=0; i<fdim; i++){
		Enew += f[i] * f[i];
	    }
	    // printf("Enew= %.10f, lambda=%g\n", Enew, lambda);///////////////////
	    if(Enew < Ecurr){
		lambda *= 0.1;
		break;
	    }else{
		lambda *= 10.0;
	    }

	    delete lud;
	}

	P = Pnew;

	double norP2 {0.0};
	double norDelta2 {0.0};
	for(int i=0; i<Pdim; i++){
	    norP2 += P[i]*P[i];
	    norDelta2 += Delta[i]*Delta[i];
	}

	// Finish condition
	// if(norDelta2/norP2 < 1.0E-10) break;
	if(norDelta2/norP2 < 1.0E-16) break;

    }
}

/************** calc ***********************/
/*** Usage example
    VecDoub X(ktSyoumen.get_n()*2);
    ktSyoumen.get_xTilder(X);
    VecDoub P0_flat(12);
    for(int i=0; i<3; i++) for(int j=0; j<4; j++) P0_flat[i*4+j] = P_nm[i][j];
    LMalgorithm lm(X, P0_flat, ktSyoumen);
    MatDoub P_LM(3, 4);
    lm.calc(P_LM);
******************************************/     
void LMalgorithm::calc(MatDoub_O &P){
    VecDoub P_flat(12);
    calc(P_flat);
    for(int i=0; i<3; i++){
	for(int j=0; j<4; j++) P[i][j] = P_flat[i*4+j];
    }
}

/************************************************/
void LMalgorithm::calcDenormalized(MatDoub_O &P){
    MatDoub P_nm(3, 4);
    calc(P_nm);
    m_kt.denormalize_P(P_nm, P);
}

/************** calc ***********************/
void LMalgorithm::calc(VecDoub_O &P){
    int n {m_kt.get_n()};
    int n2 {2*n};
    MatDoub J(n2, 12);
    VecDoub f(n2);
    VecDoub epsilon(n2);
    // VecDoub Pnew(n2+1);
    VecDoub Pnew(P.size());  // 2025.3.14
    double lambda {0.0};
    MatDoub a(12, 12);
    VecDoub b(12);
    VecDoub Delta(12);
    
    // 初期値設定
    P = m_P0;
    m_kt.J(P, J);
    for(int i=0; i<12; i++){
	double aii {0.0};
	for(int k=0; k<n2; k++) aii += J[k][i]*J[k][i];
	lambda += aii;
    }
    lambda /= 12.0;
    lambda *= 1.0E-3;
    // printf("initial lambda=%f\n", lambda);/////////////////

    while(true){
	m_kt.f(P, f);

	double Ecurr {0.0};
	for(int i=0; i<n2; i++){
	    epsilon[i] = f[i] - m_X[i];
	    Ecurr += epsilon[i]*epsilon[i];
	}
	
	// printf("Ecurr= %.10f\n", std::sqrt(Ecurr/12.0));///////////////////
	printf("Ecurr= %.10f\n", std::sqrt(Ecurr/n));// 2026.3.14

	m_kt.J(P, J);

	while(true){
	    for(int i=0; i<12; i++){
		for(int j=0; j<12; j++){
		    double sum {0.0};
		    for(int k=0; k<n2; k++) sum += J[k][i]*J[k][j];
		    a[i][j] = sum;
		}
	    }
	    for(int i=0; i<12; i++) a[i][i] += lambda;

	    for(int i=0; i<12; i++){
		double sum {0.0};
		for(int k=0; k<n2; k++) sum += J[k][i]*epsilon[k];
		b[i] = -sum;
	    }

	    LUdcmp *lud = new LUdcmp(a);

	    lud->solve(b, Delta);
	    
	    for(int i=0; i<12; i++) Pnew[i] = P[i]+Delta[i];
	    
	    m_kt.f(Pnew, f);

	    double Enew {0.0};
	    for(int i=0; i<n2; i++){
		double x {f[i] - m_X[i]};
		Enew += x * x;
	    }

	    if(Enew < Ecurr){
		lambda *= 0.1;
		Ecurr = Enew; // 2026.1.24
		break;
	    }else{
		lambda *= 10.0;
	    }

	    delete lud;
	}


	P = Pnew;

	double norP2 {0.0};
	double norDelta2 {0.0};
	for(int i=0; i<12; i++){
	    norP2 += P[i]*P[i];
	    norDelta2 += Delta[i]*Delta[i];
	}

	// Finish condition
	// if(norDelta2/norP2 < 1.0E-10){
	if(norDelta2/norP2 < 1.0E-20){
	    // printf("Final E= %.10f\n", std::sqrt(Ecurr/12.0));///////////////////
	    printf("Final E= %.10f\n", std::sqrt(Ecurr/n));// 2026.3.14
	    break;
	}
	// if(norDelta2/norP2 < 1.0E-16) break;
    }
}
