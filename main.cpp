#include <iostream>
#include "main.hpp"

#include "kizyunTen.cpp"
#include "euclidModel.cpp"
#include "DLT.cpp"
#include "LMalgorithm.cpp"
#include "linAlg.cpp"

/***********************************************/
int main(int argc, char** argv){
    KizyunTen ktSyoumen(6, "Data/IMG_0752array_A.txt", "array");  // ippan array
    
    EuclidModel em(ktSyoumen);

    VecDoub params(7);
    MatDoub P(3, 4);
    em.calc7paramsFromKizyunten(params, P);

    LinAlg la;

    LMalgorithm lm(P, ktSyoumen);
    lm.calcDenormalized(P);
    la.printMatrix(P, "P");///////////////////////

    // Calculate K from P
    MatDoub K(3, 3);
    em.calcCameraParamsYXZ(P, 0, K);

    double E2syoumen = ktSyoumen.getReprojectionErrorSquare(P);
    printf("E2syoumen= %g, sqrt(E2syoumen)=%g\n", E2syoumen, std::sqrt(E2syoumen));

    KizyunTen ktSecond(4, "Data/IMG_0753right_array_A.txt", "array");  // ippan array
    std::string strRL {"right"}; 
    // std::string strRL {"left"};
    
    EuclidModel em2(ktSecond);
    MatDoub P2(3, 4);
    em2.calcSecondCameraPHomography(ktSecond, K, P2);

    em2.calcCameraParamsYXZ(P2, 0);
    double E2second = ktSecond.getReprojectionErrorSquare(P2);
    printf("E2second= %g, sqrt(E2second)=%g\n", E2second, std::sqrt(E2second));

    LMalgorithm lm2(P2, ktSecond);
    lm2.calcDenormalized(P2);
    la.printMatrix(P2, "LM applied P2");////////////////////////
    
    em2.calcCameraParamsYXZ(P2, 0);
    
    E2second = ktSecond.getReprojectionErrorSquare(P2);
    printf("E2second= %g, sqrt(E2second)=%g\n", E2second, std::sqrt(E2second));

    double Eoverall {std::sqrt(0.5*(E2syoumen+E2second))};
    printf("Eoverall= %g\n", Eoverall);

    // triangulation
    DLT dlt(ktSyoumen);
    VecDoub PF(12), P2F(12);
    la.flatten(P, PF);
    la.flatten(P2, P2F);
    double m[2] {1213.0, 2825.0};  // IMG_0830.JPG stake tip, left
    double mDash[2] {901.0, 2596.0};  // IMG_0831.JPG stake tip, left

    VecDoub M(4);
    dlt.triangulation(PF, P2F, m, mDash, M); 
    la.printVector(M, "M", 1.0/M[3]);
}
