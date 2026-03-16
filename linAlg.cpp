/*********************************************************************
 *      linAlg.cpp
 ********************************************************************/
// #include <iostream>
// #include "main.hpp"
// #include "DLT.hpp"
// #include "svd.hpp"
// #include "ludcmp.h"

// Constructor

/****** Return matrix product
 * Assume A.ncols()==B.nrows(), A.nrows()<=C.nrows(), B.ncols()<=C.ncols() *******/
void LinAlg::prod(MatDoub_I &A, MatDoub_I &B, MatDoub_O &C){
    for(int i=0; i<A.nrows(); i++){
	for(int j=0; j<B.ncols(); j++){
	    C[i][j] = 0.0;
	    for(int k=0; k<A.ncols(); k++){
		C[i][j] += A[i][k]*B[k][j];
	    }
	}
    }
}

/****************************************************/
void LinAlg::prod(MatDoub_I &A, VecDoub_I &x, VecDoub_O &y){
    for(int i=0; i<y.size(); i++){
	y[i] = 0.0;
	for(int j=0; j<A.ncols(); j++){
	    y[i] += A[i][j]*x[j];
	}
    }
}

/****************************************************/
void LinAlg::prod(MatDoub_I &A, double x[], double y[]){
    for(int i=0; i<A.nrows(); i++){
	y[i] = 0.0;
	for(int j=0; j<A.ncols(); j++){
	    y[i] += A[i][j]*x[j];
	}
    }
}

/***********************************************************/
void LinAlg::concat(MatDoub_I &A, VecDoub_I &a, MatDoub_O &B){
    for(int i=0; i<A.nrows(); i++){
	for(int j=0; j<A.ncols(); j++){
	    B[i][j] = A[i][j];
	}
	B[i][B.ncols()-1] = a[i];
    }
}

/********************************************************/
void LinAlg::scalarMul(double c, MatDoub_I &A, MatDoub_O &cA){
    for(int i=0; i<A.nrows(); i++){
	for(int j=0; j<A.ncols(); j++){
	    cA[i][j] = c*A[i][j];
	}
    }
}

/****************************************************/
void LinAlg::exchangeRows(MatDoub_I &A, int i0, int i1, MatDoub_O &B){
    for(int j=0; j<A.ncols(); j++){
	for(int i=0; i<A.nrows(); i++){
	    if(i==i0) B[i][j] = A[i1][j];
	    else if(i==i1) B[i][j] = A[i0][j];
	    else B[i][j] = A[i][j];
	}
    }
}

/****************************************************/
void LinAlg::exchangeCols(MatDoub_I &A, int j0, int j1, MatDoub_O &B){
    for(int i=0; i<A.nrows(); i++){
	for(int j=0; j<A.ncols(); j++){
	    if(j==j0) B[i][j] = A[i][j1];
	    else if(j==j1) B[i][j] = A[i][j0];
	    else B[i][j] = A[i][j];
	}
    }
}
    
/********************************************************/
void LinAlg::rotMatrix(int rotAxis, double theta, MatDoub_O &R){
    for(int i=0; i<R.nrows(); i++){
	for(int j=0; j<R.ncols(); j++){
	    R[i][j] = 0.0;
	}
    }

    double c {std::cos(theta)};
    double s {std::sin(theta)};
    if(rotAxis==0){
	R[0][0] = 1.0;
	R[1][1] = c; R[1][2] = s;
	R[2][1] = -s; R[2][2] = c;
    }else if(rotAxis==1){
	R[0][0] = c; R[0][2] = -s;
	R[1][1] = 1.0;
	R[2][0] = s; R[2][2] = c;
    }else if(rotAxis==2){
	R[0][0] = c; R[0][1] = s;
	R[1][0] = -s; R[1][1] = c;
	R[2][2] = 1.0;
    }
}
	
/****** Print matrix ********************************/
void LinAlg::printMatrix(MatDoub_I &A, std::string str){
    std::cout << str << std::endl;
    for(int i=0; i<A.nrows(); i++){
	for(int j=0; j<A.ncols(); j++){
	    printf("[%d][%d]=%12.6f ", i, j, A[i][j]);
	}
	printf("\n");
    }
}

/****** Print matrix
 * Usage: la.printMatrix(&(A[0][0]), "A", 3, 4)   *********************/
void LinAlg::printMatrix(double* A, std::string str, int nRows, int nCols){
    std::cout << str << std::endl;
    for(int i=0; i<nRows; i++){
	for(int j=0; j<nCols; j++){
	    printf("[%d][%d]=%12.6f ", i, j, A[i*nRows+j]);
	}
	printf("\n");
    }
}

/****** Print Vector ********************************/
void LinAlg::printVector(VecDoub_I &v, std::string str){
    std::cout << str << ": ";
    for(int i=0; i<v.size(); i++){
	printf("[%d]=%12.6f ", i, v[i]);
    }
    printf("\n");
}

/****** Print Vector ********************************/
void LinAlg::printVector(VecDoub_I &v, std::string str, double keisu){
    std::cout << str << std::endl;
    for(int i=0; i<v.size(); i++){
	printf("[%d]=%16.10f ", i, v[i]*keisu);
    }
    printf("\n");
}

/****** Print Vector
 * Usage la.printVector(&X, "X", 0.1, 3) *************************/
void LinAlg::printVector(double *v, std::string str, double keisu, int size){
    std::cout << str << ": ";
    for(int i=0; i<size; i++){
	printf("[%d]=%12.6f ", i, v[i]*keisu);
    }
    printf("\n");
}

/*********************************************/
void LinAlg::calcPseudInv(MatDoub_I &A, MatDoub_O &Apinv){
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

/*********************************************/
void LinAlg::calcInv(MatDoub_I &A, MatDoub_O &Ainv){
    LUdcmp lud(A);
    lud.inverse(Ainv);
}

/***** RQ decomposition,  theta [rad]
 * ƒm[ƒg23.4ß 2025.10.16--
*************************/
void LinAlg::RQdcmp(MatDoub_I &A, int &cc, MatDoub_O &B, MatDoub_O &R,
		    VecDoub_O &theta, VecDoub_O &thetaPrime){

    double c, s;
    if(A[2][1]==0.0 && A[2][2]==0.0) theta[0] = 0.0;
    else{
	double bunbo {A[2][1]*A[2][1]+A[2][2]*A[2][2]};
	c = A[2][2]/bunbo;
	s = -A[2][1]/bunbo;
	theta[0] = std::atan2(s, c);
    }

    c = std::cos(theta[0]);
    s = std::sin(theta[0]);

    MatDoub R0T(3, 3, 0.0);
    R0T[0][0] = 1.0;
    R0T[1][1] = c; R0T[1][2] = -s;
    R0T[2][1] = s; R0T[2][2] = c;

    MatDoub A1(3, 3);
    prod(A, R0T, A1);

    if(A1[2][0]==0.0 && A1[2][2]==0.0) theta[1] = 0.0;
    else{
	double bunbo {A1[2][0]*A1[2][0]+A1[2][2]*A1[2][2]};
	c = A1[2][2]/bunbo;
	s = A1[2][0]/bunbo;
	theta[1] = std::atan2(s, c);
    }

    c = std::cos(theta[1]);
    s = std::sin(theta[1]);

    MatDoub R1T(3, 3, 0.0);
    R1T[0][0] = c; R1T[0][2] = s;
    R1T[1][1] = 1.0;
    R1T[2][0] = -s; R1T[2][2] = c;
    
    MatDoub A2(3, 3);
    prod(A1, R1T, A2);

    if(A2[1][0]==0.0 && A2[1][1]==0.0) theta[2] = 0.0;
    else{
	double bunbo {A2[1][0]*A2[1][0]+A2[1][1]*A2[1][1]};
	c = A2[1][1]/bunbo;
	s = -A2[1][0]/bunbo;
	theta[2] = std::atan2(s, c);
    }

    c = std::cos(theta[2]);
    s = std::sin(theta[2]);

    MatDoub R2T(3, 3, 0.0);
    R2T[0][0] = c; R2T[0][1] = -s;
    R2T[1][0] = s; R2T[1][1] = c;
    R2T[2][2] = 1.0;
    
    prod(A2, R2T, B);

    if(B[2][2] >= 0.0) cc = 1;
    else{
	cc = -1;
	for(int i=0; i<3; i++){
	    for(int j=0; j<3; j++){
		B[i][j] = -B[i][j];
	    }
	}
    }

    if(B[0][0]>0.0 && B[1][1]>0.0){
	if(theta[2] >= 0.0) theta[2] -= GlobalConstants::pi;
	else theta[2] += GlobalConstants::pi;

	B[0][0] = -B[0][0]; B[0][1] = -B[0][1];
	B[1][1] = -B[1][1];
    }else if(B[0][0]>0.0 && B[1][1]<=0.0){
	if(theta[2] >= 0.0) theta[2] = GlobalConstants::pi-theta[2];
	else theta[2] = -GlobalConstants::pi-theta[2];
	
	if(theta[1] >= 0.0) theta[1] -= GlobalConstants::pi;
	else theta[1] += GlobalConstants::pi;

	cc = -cc;
	B[0][0] = -B[0][0];
    }else if(B[0][0]<=0.0 && B[1][1]>0.0){
	theta[2] = -theta[2];
	
	if(theta[1] >= 0.0) theta[1] -= GlobalConstants::pi;
	else theta[1] += GlobalConstants::pi;

	cc = -cc;
	B[0][1] = -B[0][1];
	B[1][1] = -B[1][1];
    }

    c = std::cos(theta[0]);
    s = std::sin(theta[0]);
    MatDoub R0(3, 3, 0.0);
    R0[0][0] = 1.0;
    R0[1][1] =  c; R0[1][2] = s;
    R0[2][1] = -s; R0[2][2] = c;

    c = std::cos(theta[1]);
    s = std::sin(theta[1]);
    MatDoub R1(3, 3, 0.0);
    R1[0][0] = c; R1[0][2] = -s;
    R1[1][1] = 1.0;
    R1[2][0] = s; R1[2][2] = c;
    
    c = std::cos(theta[2]);
    s = std::sin(theta[2]);
    MatDoub R2(3, 3, 0.0);
    R2[0][0] =  c; R2[0][1] = s;
    R2[1][0] = -s; R2[1][1] = c;
    R2[2][2] = 1.0;

    MatDoub X(3, 3);
    prod(R2, R1, X);
    prod(X, R0, R);
    
    double pi2 {GlobalConstants::pi*0.5};
    if(theta[1] != pi2 && theta[1] != -pi2){
	if(theta[1] >= 0.0) thetaPrime[1] = GlobalConstants::pi-theta[1];
	else thetaPrime[1] = -GlobalConstants::pi-theta[1];
	
	if(theta[0] >= 0.0) thetaPrime[0] = theta[0]-GlobalConstants::pi;
	else thetaPrime[0] = theta[0]+GlobalConstants::pi;
	
	if(theta[2] >= 0.0) thetaPrime[2] = theta[2]-GlobalConstants::pi;
	else thetaPrime[2] = theta[2]+GlobalConstants::pi;
    }else{
	for(int i=0; i<3; i++) thetaPrime[i] = theta[i];
    }
}

/***** RQ decomposition,  theta [rad]  *************************/
void LinAlg::RQdcmp(MatDoub_I &A, MatDoub_O &R, MatDoub_O &Q, VecDoub_O &theta){
    // printMatrix(A, "RQdcmp");////////////
    
    if(A[2][1]==0.0 && A[2][2]==0.0){
	theta[0] = 0.0;
    }else if(A[2][2] >= 0.0){
	theta[0] = std::asin( -A[2][1]/std::sqrt(A[2][1]*A[2][1]+A[2][2]*A[2][2]) );
    }else{
	theta[0] = std::asin( A[2][1]/std::sqrt(A[2][1]*A[2][1]+A[2][2]*A[2][2]) );
    }
    // printf("theta[0]=%g\n", theta[0]);/////////////////////t

    double c {std::cos(theta[0])};
    double s {std::sin(theta[0])};
    MatDoub Q1(3, 3, 0.0);
    Q1[0][0] = 1.0; Q1[1][1] = c; Q1[1][2] = -s;
    Q1[2][1] = s; Q1[2][2] = c;
    
    MatDoub A1(3, 3);
    prod(A, Q1, A1);
    // printMatrix(A1, "A1");

    if(A1[2][0]==0.0 && A1[2][2]==0.0){
	theta[1] = 0.0;
    }else if(A1[2][2] >= 0.0){
	theta[1] = std::asin( A1[2][0]/std::sqrt(A1[2][0]*A1[2][0]+A1[2][2]*A1[2][2]) );
    }else{
	theta[1] = std::asin( -A1[2][0]/std::sqrt(A1[2][0]*A1[2][0]+A1[2][2]*A1[2][2]) );
    }
    // printf("theta[1]=%g\n", theta[1]);/////////////////////
    
    c = std::cos(theta[1]);
    s = std::sin(theta[1]);
    MatDoub Q2(3, 3, 0.0);
    Q2[0][0] = c; Q2[0][2] = s; Q2[1][1] = 1.0;
    Q2[2][0] = -s; Q2[2][2] = c;

    MatDoub A2(3, 3);
    prod(A1, Q2, A2);
    // printMatrix(A2, "A2");//////////////////////////

    if(A2[1][0]==0.0 && A2[1][1]==0.0){
	theta[2] = 0.0;
    }else if(A2[1][1] >= 0.0){
	theta[2] = std::asin( -A2[1][0]/std::sqrt(A2[1][0]*A2[1][0]+A2[1][1]*A2[1][1]) );
    }else{
	theta[2] = std::asin( A2[1][0]/std::sqrt(A2[1][0]*A2[1][0]+A2[1][1]*A2[1][1]) );
    }
    // printf("theta[2]=%g\n", theta[2]);/////////////////////
    
    c = std::cos(theta[2]);
    s = std::sin(theta[2]);
    MatDoub Q3(3, 3, 0.0);
    Q3[0][0] = c; Q3[0][1] = -s; Q3[1][0] = s;
    Q3[1][1] = c; Q3[2][2] = 1.0;

    prod(A2, Q3, R);
    // printMatrix(R, "R");////////////////////////

    // MatDoub Q1T(3, 3);
    // MatDoub Q2T(3, 3);
    // MatDoub Q3T(3, 3);
    // transpose(Q1, Q1T);
    // transpose(Q2, Q2T);
    // transpose(Q3, Q3T);

    // MatDoub X(3, 3);
    // prod(Q3T, Q2T, X);
    // prod(X, Q1T, Q);

    // printMatrix(Q, "Q");/////////////////////////

    // R[0][0]‚ÆR[1][1]‚Í“¯•„†‚ÅAR[2][2]‚Æ‚ÍˆÙ•„†‚É‚·‚éB note pp.15.3-18ˆÈ~
    bool sonomama {false};
    
    if(R[0][0] > 0.0){
	if(R[1][1] > 0.0){
	    if(R[2][2] > 0.0){  // (i) R00>0, R11>0, R22>0
		printf("(i)\n");////////////////////////
		if(theta[2] >=0.0) theta[2] -= GlobalConstants::pi;
		else theta[2] += GlobalConstants::pi;
	    }else{  // (iii) R00>0, R11>0, R22<0
		printf("(iii)\n");////////////////////////
		sonomama = true;
	    }
	}else{  // R00>0, R11<0
	    if(R[2][2] > 0.0){  // (iv) R00>0, R11<0, R22>0
		printf("(iv)\n");////////////////////////
		if(theta[0] >= 0.0) theta[0]  -= GlobalConstants::pi;
		else theta[0] += GlobalConstants::pi;
		theta[1] = -theta[1];
		theta[2] = -theta[2];
	    }else{  // (v) R00>0, R11<0, R22<0
		printf("(v)\n");////////////////////////
		if(theta[0] >= 0.0) theta[0]  -= GlobalConstants::pi;
		else theta[0] += GlobalConstants::pi;
		theta[1] = -theta[1];
		if(theta[2] >= 0.0) theta[2] = GlobalConstants::pi-theta[2];
		else theta[2] = -GlobalConstants::pi-theta[2];
	    }
	}
    }else{   // R00<0
	if(R[1][1] > 0.0){   // R00<0, R11>0
	    if(R[2][2] > 0.0){   // (vi) R00<0, R11>0, R22>0
		printf("(vi)\n");////////////////////////
		if(theta[0] >= 0.0) theta[0]  -= GlobalConstants::pi;
		else theta[0] += GlobalConstants::pi;
		theta[1] = -theta[1];
		if(theta[2] >= 0.0) theta[2] = GlobalConstants::pi-theta[2];
		else theta[2] = -GlobalConstants::pi-theta[2];
	    }else{   // (vii) R00<0, R11>0, R22<0
		printf("(vii)\n");////////////////////////
		if(theta[0] >= 0.0) theta[0]  -= GlobalConstants::pi;
		else theta[0] += GlobalConstants::pi;
		theta[1] = -theta[1];
		theta[2] = -theta[2];
	    }
	}else{   // R00<0, R11<0
	    if(R[2][2] > 0.0){   // (viii) R00<0, R11<0, R22>0
		printf("(viii)\n");////////////////////////
		sonomama = true;
	    }else{   // (ii) R00<0, R11<0, R22<0
		printf("(ii)\n");////////////////////////
		if(theta[2] >= 0.0) theta[2] -= GlobalConstants::pi;
		else theta[2] += GlobalConstants::pi;
	    }
	}
    }

    MatDoub X(3, 3);
    
    if(sonomama == false){
	c = std::cos(theta[0]);
	s = std::sin(theta[0]);
	Q1[1][1] = c; Q1[1][2] = -s;
	Q1[2][1] = s; Q1[2][2] = c;
	
	c = std::cos(theta[1]);
	s = std::sin(theta[1]);
	Q2[0][0] = c; Q2[0][2] = s;
	Q2[2][0] = -s; Q2[2][2] = c;
	
	c = std::cos(theta[2]);
	s = std::sin(theta[2]);
	Q3[0][0] = c; Q3[0][1] = -s;
	Q3[1][0] = s; Q3[1][1] = c;

	prod(A, Q1, X);
	MatDoub Y(3, 3);
	prod(X, Q2, Y);
	prod(Y, Q3, R);
    }
	
    MatDoub Q1T(3, 3);
    MatDoub Q2T(3, 3);
    MatDoub Q3T(3, 3);
    
    transpose(Q1, Q1T);
    transpose(Q2, Q2T);
    transpose(Q3, Q3T);

    prod(Q3T, Q2T, X);
    prod(X, Q1T, Q);
	
}

/************ Take the transpose of matirx *******************/
void LinAlg::transpose(MatDoub_I &A, MatDoub_O &AT){
    for(int i=0; i<A.nrows(); i++){
	for(int j=0; j<A.ncols(); j++){
	    AT[j][i] = A[i][j];
	}
    }
}

/***************** 2025.9.17 ****************************/
void LinAlg::testRQdcmp(){
    MatDoub R(3, 3, 0.0);
    R[0][0] = 1.0; R[0][1] = 0.1; R[0][2] = 3.0;
    R[1][1] = -0.5; R[1][2] = -2.0;
    R[2][2] = 1.0;

    printMatrix(R, "R orig");

    double theta_x {50.0/180.0*GlobalConstants::pi};
    double theta_y {-40.0/180.0*GlobalConstants::pi};
    double theta_z {70.0/180.0*GlobalConstants::pi};

    VecDoub theta(3);
    theta[0] = theta_x; theta[1] = theta_y; theta[2] = theta_z;

    printVector(theta, "theta orig [deg]", 180.0/GlobalConstants::pi);

    MatDoub QxT(3, 3, 0.0);
    MatDoub QyT(3, 3, 0.0);
    MatDoub QzT(3, 3, 0.0);

    double c {std::cos(theta_x)};
    double s {std::sin(theta_x)};
    QxT[0][0] = 1.0;
    QxT[1][1] = c; QxT[1][2] = s;
    QxT[2][1] = -s; QxT[2][2] = c;
    
    c = std::cos(theta_y);
    s = std::sin(theta_y);
    QyT[0][0] = c; QyT[0][2] = -s;
    QyT[1][1] = 1.0;
    QyT[2][0] = s; QyT[2][2] = c;
    
    c = std::cos(theta_z);
    s = std::sin(theta_z);
    QzT[0][0] = c; QzT[0][1] = s;
    QzT[1][0] = -s; QzT[1][1] = c;
    QzT[2][2] = 1.0;

    MatDoub A(3, 3);
    MatDoub Q(3, 3);
    MatDoub X(3, 3);
    MatDoub Y(3, 3);
    prod(QzT, QyT, X);
    prod(X, QxT, Q);
    prod(R, Q, A);

    printMatrix(Q, "Q orig");
    
    RQdcmp(A, R, Q, theta);
    printMatrix(R, "R");
    printMatrix(Q, "Q");
    printVector(theta, "theta [deg]", 180.0/GlobalConstants::pi);
}

/******************************************/
void LinAlg::flatten(MatDoub_I &P, VecDoub_O &Pf){
    int k {0};
    for(int i=0; i<P.nrows(); i++){
	for(int j=0; j<P.ncols(); j++){
	    Pf[k] = P[i][j];
	    k++;
	}
    }
}

/********************************************/
double LinAlg::getNorm(VecDoub_I &x){
    double ret {0.0};
    for(int i=0; i<x.size(); i++) ret += x[i]*x[i];
    return std::sqrt(ret);
}

/**********************************************/
double LinAlg::getScalarProd(VecDoub_I &x, VecDoub_I &y){
    double ret {0.0};
    for(int i=0; i<x.size(); i++) ret += x[i]*y[i];
    return ret;
}

/*****************************************************/
void LinAlg::scalarMul(double c, VecDoub_I &X, VecDoub_O &Y){
    for(int i=0; i<X.size(); i++) Y[i] = c*X[i];
}

/**************************************************/
void LinAlg::add(VecDoub_I &x, VecDoub_I &y, VecDoub_O &z){
    for(int i=0; i<x.size(); i++) z[i] = x[i]+y[i];
}

/**************************************************/
void LinAlg::subtract(VecDoub_I &x, VecDoub_I &y, VecDoub_O &z){
    for(int i=0; i<x.size(); i++) z[i] = x[i]-y[i];
}

/**************************************************/
void LinAlg::calcVectorProd(VecDoub_I &x, VecDoub_I &y, VecDoub_O &z){
    z[0] = x[1]*y[2]-x[2]*y[1];
    z[1] = x[2]*y[0]-x[0]*y[2];
    z[2] = x[0]*y[1]-x[1]*y[0];
}
