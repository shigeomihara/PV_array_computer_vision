/****************************************************************
 *      main.hpp
 ***************************************************************/
#if !defined(MAIN_HPP)
#define MAIN_HPP
#include <string>
#include <vector>
#include "nr3.h"
#include "ludcmp.h"
#include "svd.h"
// #include <opencv2/core.hpp>

/****** KizyunTen in Japanese == reference point **************/
class KizyunTen{
private:
    char m_ABCD;
    const std::string m_txtFilePath;
    int m_n;   // number of reference points
    void readTxtFile(std::string txtFilePath);
    void setKizyun3D_31();
    void setKizyun3D_410();
    void setKizyun3D_box();
    void setKizyun3D_array();
    
public:
    // Constructors
    KizyunTen();
    KizyunTen(int n, std::string photoPath, std::string photoType);

    // name list of reference points
    std::string m_kizName[20];
    
    // world coordinates
    double m_Hw[3];
    double m_He[3];  // X, Y, Z
    double m_Lw[3];
    double m_Le[3];
    double m_Mw[3];
    double m_Me[3];
    double m_MHw[3];
    double m_MHe[3];
    double m_Hcw[3];
    double m_Hce[3];
    double m_Lcw[3];
    double m_Lce[3];
    double m_MHcw[3];
    double m_MHce[3];
    double m_Mcw[3];
    double m_Mce[3];
    
    double m_Phw[3];  // middle fence pole top
    double m_Phe[3];  // middle fence pole top
 
    double m_Plw[3];  // middle fence pole bottom
    double m_Ple[3];  // middle fence pole bottom
    
    double m_Bpw[3];  // Back array highest panel west corner
    double m_Bpe[3];  // Back array highest panel east corner
    
    double m_Mh[3];
    double m_Hc[3];
    double m_Lc[3];
    double m_Mww[3];
    double m_Mee[3];

    // double m_xTilder[m_n+1][3];
    // double m_Xtilder[m_n+1][4];
    // double m_xTilderOrig[m_n+1][3]; // original before normalization
    // double m_XtilderOrig[m_n+1][4]; // original before normalization

    // Take arrays suffieiently large
    double m_xTilder[20][2];
    double m_Xtilder[20][3];
    double m_xTilderOrig[20][2]; // original before normalization
    double m_XtilderOrig[20][3]; // original before normalization

    // image coordinates
    int m_width;   // image width
    int m_height;  // image height
    int m_hw[2];  // x, y
    int m_he[2];  // x, y
    int m_lw[2];  // x, y
    int m_le[2];  // x, y
    int m_mw[2];  // x, y
    int m_me[2];  // x, y
    int m_mhw[2];  // x, y
    int m_mhe[2];  // x, y
    int m_hcw[2];  // x, y
    int m_hce[2];  // x, y
    int m_lcw[2];  // x, y
    int m_lce[2];  // x, y
    int m_mhcw[2];  // x, y
    int m_mhce[2];  // x, y
    int m_mcw[2];  // x, y
    int m_mce[2];  // x, y
    
    int m_phw[2];
    int m_phe[2];
    
    int m_plw[2];
    int m_ple[2];
    
    int m_bpw[2];
    int m_bpe[2];
    
    int m_mh[2];
  
    // normalization
    double m_Ts;
    double m_Tt[2];   // t[0], t[1]
    double m_Us;
    double m_Ut[3];   // t[0], t[1], t[2]
    double m_Tinv[3][3];   // 3x3 matrix
    double m_U[4][4];   // 4x4 matrix  UX = normalized X
    double m_T[3][3];   // Tx = normalized x
    double m_Uinv[4][4];
  
    int get_n() const{ return m_n; }
    void get_xTilderOrig(VecDoub_O &xTilderOrig)const;
    void get_xTilder(VecDoub_O &X)const;

    // for LM algorithm
    void f(VecDoub_I &P, VecDoub_O &fval);
    void J(VecDoub_I &P, MatDoub_O &Jval);
    void Jnum(VecDoub &P, MatDoub_O &Jval);
    void fOrig(VecDoub_I &P, VecDoub_O &fval);
    
    static std::vector<std::string> splitString(std::string str, char c);

    void denormalize_x(double x[4]);
    void denormalize_M(VecDoub_IO &M);
    void normalize_m(double m[3]);
    void normalize(const double x[], const double y[], VecDoub_O &u, VecDoub_O &v,
		   MatDoub &T);

    VecDoub m_Pbar {VecDoub(13)};   // normalized camera matrix
    void calcPbar();

    double getGeometricError(VecDoub_I &P);
    double getGeometricError(MatDoub_I &P);

    void normalize_P(MatDoub_I &P, MatDoub_O &P_nm);
    void denormalize_P(MatDoub_I &P_nm, MatDoub_O &P);

    void printWorld2D();    // for OpenCV camera calib.
    
    double getReprojectionErrorSquare(VecDoub_I &P);
    double getReprojectionErrorSquare(MatDoub_I &P);
};

/************************************************************/
class GlobalConstants{
public:
    constexpr static double pi = 3.141592653589793;
    const double deltaX {19.385};   // module gap in east-west direction (mm)
    const double deltaY {20.0};   // module gap in north-south direction  (mm)
    const double mX {1652.0};   // module length (mm)
    const double mY {994.0};   // module width (mm)
    const double theta_p {10.0};  // array tilt angle (deg)
    const double hL {550.0};    // height of module south edge (mm)
    const double dX {32.5};    // (mm)
    const double dY {22.5};    // (mm)
    const double aC {156.0};   // (mm)
    const double epsilonX {3.0};  // (mm)
    const double epsilonY {2.6};  // (mm)
 };

/************************************************************/
class EuclidModel{
    KizyunTen m_kizyunTen;
    double m_alpha;
    double m_c;
    double m_s;
    double m_t[3];
    double m_Ctilder[3];  // camera center
    double m_theta_rad;
    
    double m_w2;   // 0.5*m_kizyunTen.m_width
    double m_h2;   // 0.5*m_kizyunTen.m_height

    // general model
    double m_thetax;   // rad
    double m_thetay;   // rad
    double m_thetaz;   // rad

    void setAb(KizyunTen &ktRight, MatDoub_IO &A, VecDoub_IO &b, int i, int j, int Ai);
public:
    EuclidModel();
    EuclidModel(KizyunTen kt);
    // double m_P[3][4];
    MatDoub m_P;
    // double m_Pprime[3][4];
    void calc();
    void calcP();
    void calcSimple();
    void getParams(VecDoub &P);
    KizyunTen getKizyunTen();
    void get_f0(VecDoub_I &P, VecDoub_O &fval);
    void getJ0(VecDoub_I &P, MatDoub_O &Jval);
    void getJnum0(VecDoub &P, MatDoub_O &Jval);
    void calcPgene(VecDoub_I &P, bool boolPrint);
    void get_F(VecDoub_I &vecAlpha, VecDoub_O &fval);
    void getJ(VecDoub_I &vecAlpha, MatDoub_O &Jval);
    void getJnum(VecDoub &vecAlpha, MatDoub_O &Jval);
    void calcPprime(KizyunTen &ktRight);
    void calcPprime2(KizyunTen &ktRight, VecDoub_O &vecAlpha);
    void printImageError();
    double getAlpha()const{ return m_alpha; };
    void setAlpha(double alpha_val){ m_alpha = alpha_val; };
    void calcCameraParams(MatDoub_I &P, int positiveAngleIdx, MatDoub_O &K,
			  MatDoub_O &R, VecDoub_O &theta,
			  VecDoub_O &thetaPrime, VecDoub_O &t);
    void calcCameraParamsYXZ(MatDoub_I &P, int positiveAngleIdx,
			     MatDoub_O &K, MatDoub_O &R,
			     VecDoub_O &theta, VecDoub_O &thetaPrime,
			     VecDoub_O &t);
    void calcCameraParamsYXZ(MatDoub_I &P, int positiveAngleIdx, MatDoub_O &K);
    void calcCameraParamsYXZ(MatDoub_I &P, int positiveAngleIdx);
    double getGeometricError();
    void calc7paramsFromKizyunten(VecDoub_O &params, MatDoub_O &P);
    void denormalizeKt(KizyunTen &kizyunTen, MatDoub_IO &K, MatDoub_I &R, VecDoub_IO &t);
    void calcSecondCameraP(MatDoub_I &K, MatDoub_O &P, std::string strRL);
    void calcSecondCameraPfromK(MatDoub_I &K, MatDoub_O &P, std::string strRL);
    void calcSecondCameraP(VecDoub_I &params, MatDoub_O &P, std::string strRL);
    void calcSecondCameraPHomography(KizyunTen kt, MatDoub_I &K, MatDoub_O &P);

    // Linear Algebra
    void prod(MatDoub_I &A, MatDoub_I &B, MatDoub_O &C);   // C=AB
    void printMatrix(MatDoub_I &A, std::string str);
    void calcPseudInv(MatDoub_I &A, MatDoub_O &Apinv);
};

/************************************************************/
class DLT{
private:
    KizyunTen m_kt;
    const int m_n;
  
public:
    DLT(KizyunTen kt);  // Constructor
  
    void calcP(VecDoub_O &P);
    void denormalize(VecDoub_IO &P, VecDoub_O &Ctilder, MatDoub_O &Minv);
    void triangulation(VecDoub_I &P, VecDoub_I &Pdash,
		       double m[], double mDash[], VecDoub_O &M);
};

/************************************************************/
class LMalgorithm{
private:
    VecDoub m_X;  // f(P)=X
    VecDoub m_P0; // initial P
    KizyunTen m_kt;
    EuclidModel m_em;
    
public:
    LMalgorithm(VecDoub_I &X, VecDoub_I &P0, KizyunTen &kt);
    LMalgorithm(VecDoub_I &X, VecDoub_I &P0, EuclidModel &em);
    // LMalgorithm(Homography &hm);
    LMalgorithm(MatDoub_I &P0, KizyunTen &kt);
    // Homography m_hom;
    void calc(VecDoub_O &P);
    void calc(MatDoub_O &P);
    void calcDenormalized(MatDoub_O &P);
    void calcEuclidModel0(VecDoub_O &P);  // Not used. linear
    void calcEuclidModel(VecDoub_O &vecAlpha);
    void calcHomography();
};

/******************************************************/
class LinAlg{
public:
    void prod(MatDoub_I &A, MatDoub_I &B, MatDoub_O &C);   // C=AB
    void prod(MatDoub_I &A, VecDoub_I &x, VecDoub_O &y);
    void prod(MatDoub_I &A, double x[], double y[]);
    void concat(MatDoub_I &A, VecDoub_I &a, MatDoub_O &B);
    void scalarMul(double c, MatDoub_I &A, MatDoub_O &cA);
    void scalarMul(double c, VecDoub_I &X, VecDoub_O &Y);
    void exchangeRows(MatDoub_I &A, int i0, int i1, MatDoub_O &B);
    void exchangeCols(MatDoub_I &A, int j0, int j1, MatDoub_O &B);
    void rotMatrix(int rotAxis, double theta, MatDoub_O &R);
    void printVector(VecDoub_I &v, std::string str);
    void printVector(VecDoub_I &v, std::string str, double keisu);
    void printVector(double *v, std::string str, double keisu, int size);
    void printMatrix(MatDoub_I &A, std::string str);
    void printMatrix(double* A, std::string str, int nRows, int nCols);
    void calcPseudInv(MatDoub_I &A, MatDoub_O &Apinv);
    void calcInv(MatDoub_I &A, MatDoub_O &Ainv);
    void RQdcmp(MatDoub_I &A, MatDoub_O &R, MatDoub_O &Q, VecDoub_O &theta);
    void RQdcmp(MatDoub_I &A, int &c, MatDoub_O &B, MatDoub_O &R,
		VecDoub_O &theta, VecDoub_O &thetaPrime);
    void transpose(MatDoub_I &A, MatDoub_O &AT);
    void testRQdcmp();
    void flatten(MatDoub_I &P, VecDoub_O &Pf);
    double getNorm(VecDoub_I &x);
    double getScalarProd(VecDoub_I &x, VecDoub_I &y);
    void add(VecDoub_I &x, VecDoub_I &y, VecDoub_O &z);
    void subtract(VecDoub_I &x, VecDoub_I &y, VecDoub_O &z);
    void calcVectorProd(VecDoub_I &x, VecDoub_I &y, VecDoub_O &z);
};

#endif
