// #include <iostream>
// #include <cmath>
// #include <fstream>
// #include <vector>
// #include "main.hpp"
// #include "kizyunTen.hpp"
// #include "globalConstants.hpp"
// #include "DLT.hpp"

/*** Constructor doing nothing ****/
KizyunTen::KizyunTen(){
}

/*** Constructor ***/
KizyunTen::KizyunTen(int n, std::string txtFilePath, std::string photoType)
    :m_n(n), m_txtFilePath(txtFilePath){

    // std::vector<std::string> vs = splitString(txtFilePath, '.');
    // std::string secondLastString = vs[vs.size()-2];
    // m_ABCD = secondLastString[secondLastString.size()-1];
    m_ABCD = 'A';

    if(photoType=="410") setKizyun3D_410();
    else if(photoType=="box") setKizyun3D_box();
    else if(photoType=="array") setKizyun3D_array();
    else setKizyun3D_31();

    for(int i=0; i<m_n; i++){///////////////////
      for(int j=0; j<3; j++){
	  std::cout << "Xtilder[" << i << "][" << j << "]= " << m_Xtilder[i][j] << " "
		    << m_kizName[i] << " ";
      }
      std::cout << std::endl;
    }

    for(int i=0; i<m_n; i++){///////////////////
        for(int j=0; j<2; j++){
            std::cout << "xTilder[" << i << "][" << j << "]=" << m_xTilder[i][j] << " "
		<< m_kizName[i] << " ";
        }
        std::cout << std::endl;
    }
    // exit(0);//////////////////////////

    // normalization
    double m[2];  // m[j], j=0,1
    for(int j=0; j<2; j++){
        m[j] = 0.0;
        for(int i=0; i<m_n; i++){
            m[j] += m_xTilder[i][j];
        }
        m[j] /= m_n;
    }

    double sum2 {0.0};
    for(int i=0; i<m_n; i++){
        for(int j=0; j<2; j++){
            sum2 += (m_xTilder[i][j] -m[j])*(m_xTilder[i][j] -m[j]);
        }
    }
    
    double s = sqrt(2.0)/sqrt( sum2/m_n );

    double t[3];
    for(int j=0; j<2; j++){
        t[j] = -s * m[j];
    }

    m_Ts = s;
    for(int j=0; j<2; j++) m_Tt[j] = t[j];

    // Matrix m_T[3][3]
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            m_T[i][j] = 0.0;
        }
    }
    m_T[0][0] = m_Ts;
    m_T[1][1] = m_Ts;
    m_T[0][2] = m_Tt[0];
    m_T[1][2] = m_Tt[1];
    m_T[2][2] = 1.0;

    // save original
    for(int i=0; i<m_n; i++){
        for(int j=0; j<2; j++){
            m_xTilderOrig[i][j] = m_xTilder[i][j];
        }
    }

    // normalization
    for(int i=0; i<m_n; i++){
        for(int j=0; j<2; j++){
            m_xTilder[i][j] = m_Ts*m_xTilder[i][j] + m_Tt[j];
        }
    }

    // normalize world coordinates
    double M[3];
    for(int j=0; j<3; j++){
        M[j] = 0.0;
        for(int i=0; i<m_n; i++){
            M[j] += m_Xtilder[i][j];
        }
        M[j] /= m_n;
    }

    sum2 = 0.0;
    for(int i=0; i<m_n; i++){
        for(int j=0; j<3; j++){
            sum2 += (m_Xtilder[i][j] -M[j])*(m_Xtilder[i][j] -M[j]);
        }
    }

    s = sqrt(3.0)/sqrt( sum2/m_n );

    for(int j=0; j<3; j++){
        t[j] = -s * M[j];
    }

    m_Us = s;
    for(int j=0; j<3; j++) m_Ut[j] = t[j];
    
    // Matrix m_U[4][4]
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            m_U[i][j] = 0.0;
        }
    }
    m_U[0][0] = m_Us;
    m_U[1][1] = m_Us;
    m_U[2][2] = m_Us;
    m_U[0][3] = m_Ut[0];
    m_U[1][3] = m_Ut[1];
    m_U[2][3] = m_Ut[2];
    m_U[3][3] = 1.0;

    // Matrix m_Uinv[4][4]
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            m_Uinv[i][j] = 0.0;
        }
    }
    double UsInv {1.0/m_Us};
    m_Uinv[0][0] = UsInv;
    m_Uinv[1][1] = UsInv;
    m_Uinv[2][2] = UsInv;
    m_Uinv[0][3] = -UsInv*m_Ut[0];
    m_Uinv[1][3] = -UsInv*m_Ut[1];
    m_Uinv[2][3] = -UsInv*m_Ut[2];
    m_Uinv[3][3] = 1.0;

    // save original
    for(int i=0; i<m_n; i++){
        for(int j=0; j<3; j++){
            m_XtilderOrig[i][j] = m_Xtilder[i][j];
        }
    }

    // normalization
    for(int i=0; i<m_n; i++){
        for(int j=0; j<3; j++){
            m_Xtilder[i][j] = m_Us*m_Xtilder[i][j] + m_Ut[j];
        }
    }

    // std::cout << "normazlized xTilder\n";
    // for(int i=0; i<m_n; i++){///////////////////
    //     for(int j=0; j<2; j++){
    //         std::cout << "xTilder[" << i << "][" << j << "]=" << m_xTilder[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "normalized Xtilder\n";
    // for(int i=0; i<m_n; i++){///////////////////
    //     for(int j=0; j<3; j++){
    //         std::cout << "Xtilder[" << i << "][" << j << "]=" << m_Xtilder[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // exit(0);////////////////////

    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            m_Tinv[i][j] = 0.0;
        }
    }
    m_Tinv[0][0] = 1.0/m_Ts;
    m_Tinv[1][1] = 1.0/m_Ts;
    m_Tinv[0][2] = -m_Tt[0]/m_Ts;
    m_Tinv[1][2] = -m_Tt[1]/m_Ts;
    m_Tinv[2][2] = 1.0;
  
    // std::cout << "T*Tinv\n";////////////////////
    // for(int i=0; i<3; i++){
    //     for(int j=0; j<3; j++){
    //         double x {0.0};
    //         for(int k=0; k<3; k++){
    //             x += m_T[i][k]*m_Tinv[k][j];
    //         }
    //         printf("%10f ", x);
    //     }
    //     printf("\n");
    // }
}

/*** Print 2D world coordinate. For Opencv camera calib.  *****************************/
void KizyunTen::printWorld2D(){
    GlobalConstants gc;
    
    double X {0.0};
    double Y {0.0};
    printf("%g, %g\n", X, Y);
    X = gc.mX;
    printf("%g, %g\n", X, Y);
    X = 0.0; Y = gc.mY;
    printf("%g, %g\n", X, Y);
    X = gc.mX; Y = gc.mY;
    printf("%g, %g\n", X, Y);
    X = 0.0; Y = gc.mY*4.0+gc.deltaY*3.0;
    printf("%g, %g\n", X, Y);
    X = gc.mX; Y = gc.mY*4.0+gc.deltaY*3.0;
    printf("%g, %g\n", X, Y);
}

/*** Compute the geometric erro for the camera matrix P  *****************************/
double KizyunTen::getGeometricError(VecDoub_I &P){
    double E {0.0};
    VecDoub fval(2*m_n);
    VecDoub X(2*m_n);

    fOrig(P, fval);
    get_xTilderOrig(X);

    for(int i=0; i<2*m_n; i++){
	double epsilon = fval[i] - X[i];
	E += epsilon*epsilon;
    }
    
    return E;
}

/*** Compute the geometric erro for the camera matrix P  *****************************/
double KizyunTen::getGeometricError(MatDoub_I &P){
    VecDoub vecP(P.nrows()*P.ncols());
    for(int i=0; i<P.nrows(); i++){
	for(int j=0; j<P.ncols(); j++){
	    vecP[i*P.ncols()+j] = P[i][j];
	}
    }

    return getGeometricError(vecP);
}

/*** Compute the re-projection error for the camera matrix P, Mar.10,2026 ****************/
double KizyunTen::getReprojectionErrorSquare(VecDoub_I &P){
    double E {0.0};
    VecDoub fval(2*m_n);
    VecDoub X(2*m_n);

    fOrig(P, fval);
    get_xTilderOrig(X);

    for(int i=0; i<2*m_n; i++){
	double epsilon = fval[i] - X[i];
	E += epsilon*epsilon;
    }

    //////////////////////t
    for(int i=0; i<2*m_n; i++){
	printf("fval[%d]=%g, X[%d]=%g\n", i, fval[i], i, X[i]);
    }

    double ave {0.0};//////////////////
    for(int i=0; i<2*m_n; i++) ave += X[i];
    printf("average image point coordinate %g, 2*m_n=%d\n", ave/(2*m_n), 2*m_n);
    
    return E/m_n;
}

/*** Compute the re-projection error for the camera matrix P, Mar.10,2026  ****************/
double KizyunTen::getReprojectionErrorSquare(MatDoub_I &P){
    VecDoub vecP(P.nrows()*P.ncols());
    for(int i=0; i<P.nrows(); i++){
	for(int j=0; j<P.ncols(); j++){
	    vecP[i*P.ncols()+j] = P[i][j];
	}
    }

    return getReprojectionErrorSquare(vecP);
}

/*** 1 kouku, the southest array reference points  *****************************/
void KizyunTen::setKizyun3D_31(){
    GlobalConstants gc;
    
}

/*** Ippan array (410) reference points  *****************************/
void KizyunTen::setKizyun3D_410(){
    GlobalConstants gc;
  
    // std::cout << "ABCD=" << m_ABCD << std::endl;  /////////////////////t

    // Set world coordinates
    if(m_ABCD == 'A'){
        m_He[0] = 0.0;
        m_He[1] = gc.hL + gc.mY*std::sin(gc.theta_p*gc.pi/180.0);
        m_He[2] = gc.mY*std::cos(gc.theta_p*gc.pi/180.0); 

        m_Hw[0] = gc.mX;
        m_Hw[1] = m_He[1];
        m_Hw[2] = m_He[2];

        m_Me[0] = 0.0;
        // m_Me[0] = m_MHe[0];
        // m_Me[0] = m_Hce[0];   // Be carefull !!!
        m_Me[1] = gc.hL + (2.0*gc.aC +gc.dY +1.5*gc.epsilonY)*std::sin(gc.theta_p*gc.pi/180.0);
        m_Me[2] = (2.0*gc.aC +gc.dY +1.5*gc.epsilonY)*std::cos(gc.theta_p*gc.pi/180.0);

        m_Mw[0] = m_Hw[0];
        // m_Mw[0] = m_MHw[0];
        // m_Mw[0] = m_Hcw[0];   // Be carefull !!!
        m_Mw[1] = m_Me[1];
        m_Mw[2] = m_Me[2];
    
        m_Le[0] = 0.0;
        m_Le[1] = gc.hL;
        m_Le[2] = 0.0;

        m_Lw[0] = gc.mX;
        m_Lw[1] = m_Le[1];
        m_Lw[2] = m_Le[2];

        m_Hce[0] = gc.dX +3.0*gc.aC +2.5*gc.epsilonX;
        m_Hce[1] = gc.hL + (gc.mY-gc.dY)*std::sin(gc.theta_p*gc.pi/180.0);
        m_Hce[2] = (gc.mY-gc.dY)*std::cos(gc.theta_p*gc.pi/180.0); 

        m_Hcw[0] = gc.dX +7.0*gc.aC +6.5*gc.epsilonX;
        m_Hcw[1] = m_Hce[1];
        m_Hcw[2] = m_Hce[2];

        m_MHe[0] = gc.dX;
        m_MHe[1] = gc.hL + (4.0*gc.aC +gc.dY +3.5*gc.epsilonY)*std::sin(gc.theta_p*gc.pi/180.0);
        m_MHe[2] = (4.0*gc.aC +gc.dY +3.5*gc.epsilonY)*std::cos(gc.theta_p*gc.pi/180.0);

        m_MHw[0] = gc.mX -gc.dX;
        m_MHw[1] = m_MHe[1];
        m_MHw[2] = m_MHe[2];
	
        m_MHce[0] = m_Hce[0];
        m_MHce[1] = m_MHe[1];
        m_MHce[2] = m_MHe[2];

        m_MHcw[0] = m_Hcw[0];
        m_MHcw[1] = m_MHe[1];
        m_MHcw[2] = m_MHe[2];

        m_Mce[0] = m_Hce[0];
        m_Mce[1] = m_Me[1];
        m_Mce[2] = m_Me[2];

        m_Mcw[0] = m_Hcw[0];
        m_Mcw[1] = m_Me[1];
        m_Mcw[2] = m_Me[2];

        m_Lce[0] = m_Hce[0];
        m_Lce[1] = gc.hL + gc.dY*std::sin(gc.theta_p*gc.pi/180.0);
        m_Lce[2] = gc.dY*std::cos(gc.theta_p*gc.pi/180.0);

        m_Lcw[0] = m_Hcw[0];
        m_Lcw[1] = m_Lce[1];
        m_Lcw[2] = m_Lce[2];

        m_Hc[0] = 0.5*gc.mX;
        m_Hc[1] = m_Hce[1];
        m_Hc[2] = m_Hce[2];
	
        // m_Lc[0] = 0.5*gc.mX;
        // // m_Lc[1] = m_Ph[1];
        // // m_Lc[2] = m_Ph[2];
        // m_Lc[1] = gc.hp;
        // m_Lc[2] = -gc.dpL;
	
        // m_Mww[0] = gc.mX;
        // m_Mww[1] = gc.hL + 0.5*gc.mY*std::sin(gc.theta_p*gc.pi/180.0);
        // m_Mww[2] = 0.5*gc.mY*std::cos(gc.theta_p*gc.pi/180.0);
	
        // m_Mee[0] = 0.0;
        // m_Mee[1] = m_Mww[1];
        // m_Mee[2] = m_Mww[2];
	
        // m_Phw[0] = gc.mX -gc.dX -2.0*gc.aC -1.5*gc.epsilonX;
        // m_Phw[1] = gc.hL +gc.dY*std::sin(gc.theta_p*gc.pi/180.0);
        // m_Phw[2] = gc.dY*std::cos(gc.theta_p*gc.pi/180.0);
	
        // m_Phe[0] = gc.dX +2.0*gc.aC +1.5*gc.epsilonX;
        // m_Phe[1] = gc.hL +gc.dY*std::sin(gc.theta_p*gc.pi/180.0);
        // m_Phe[2] = gc.dY*std::cos(gc.theta_p*gc.pi/180.0);
	
        // m_Plw[0] = gc.mX -gc.dX -2.0*gc.aC -1.5*gc.epsilonX;
        // m_Plw[1] = gc.hL -gc.ironBarBottom*std::cos(gc.theta_p*gc.pi/180.0);
        // m_Plw[2] = gc.ironBarBottom*std::sin(gc.theta_p*gc.pi/180.0);
	
        // m_Ple[0] = gc.dX +2.0*gc.aC +1.5*gc.epsilonX;
        // m_Ple[1] = gc.hL -gc.ironBarBottom*std::cos(gc.theta_p*gc.pi/180.0);
        // m_Ple[2] = gc.ironBarBottom*std::sin(gc.theta_p*gc.pi/180.0);
	
        // m_Bpe[0] = 0.0;
        // m_Bpe[1] = gc.hL +gc.AY*std::sin(gc.theta_p*gc.pi/180.0);
        // m_Bpe[2] = -gc.arrayGap;
	
        // m_Bpw[0] = gc.mX;
        // m_Bpw[1] = gc.hL +gc.AY*std::sin(gc.theta_p*gc.pi/180.0);
        // m_Bpw[2] = -gc.arrayGap;
    }
  
    // std::cout << "Hw.X=" << m_Hw[0] << " Hw.Y=" << m_Hw[1] << " Hw.Z=" << m_Hw[2] << std::endl;//////////t
    // std::cout << "He.X=" << m_He[0] << " He.Y=" << m_He[1] << " He.Z=" << m_He[2] << std::endl;
    // std::cout << "Mw.X=" << m_Mw[0] << " Mw.Y=" << m_Mw[1] << " Mw.Z=" << m_Mw[2] << std::endl;
    // std::cout << "Me.X=" << m_Me[0] << " Me.Y=" << m_Me[1] << " Me.Z=" << m_Me[2] << std::endl;
    // std::cout << "Lw.X=" << m_Lw[0] << " Lw.Y=" << m_Lw[1] << " Lw.Z=" << m_Lw[2] << std::endl;
    // std::cout << "Le.X=" << m_Le[0] << " Le.Y=" << m_Le[1] << " Le.Z=" << m_Le[2] << std::endl;
    // std::cout << "Phe.X=" << m_Phe[0] << " Phe.Y=" << m_Phe[1] << " Phe.Z=" << m_Phe[2] << std::endl;
    // std::cout << "Ple.X=" << m_Ple[0] << " Ple.Y=" << m_Ple[1] << " Ple.Z=" << m_Ple[2] << std::endl;
    // std::cout << "Phw.X=" << m_Phw[0] << " Phw.Y=" << m_Phw[1] << " Phw.Z=" << m_Phw[2] << std::endl;
    // std::cout << "Plw.X=" << m_Plw[0] << " Plw.Y=" << m_Plw[1] << " Plw.Z=" << m_Plw[2] << std::endl;
    // std::cout << "Bpe.X=" << m_Bpe[0] << " Bpe.Y=" << m_Bpe[1] << " Bpe.Z=" << m_Bpe[2] << std::endl;
    // std::cout << "Bpw.X=" << m_Bpw[0] << " Bpw.Y=" << m_Bpw[1] << " Bpw.Z=" << m_Bpw[2] << std::endl;
    // exit(0);

    readTxtFile(m_txtFilePath);

    // for(int i=0; i<m_n; i++){///////////////////
    // 	std::cout << m_kizName[i];
    // 	for(int j=0; j<3; j++){
    // 	  std::cout << " Xtilder[" << i << "][" << j << "]= " << m_Xtilder[i][j];
    // 	}
    // 	std::cout << std::endl;
    // 	std::cout << m_kizName[i];
    //     for(int j=0; j<2; j++){
    //         std::cout << " xTilder[" << i << "][" << j << "]= " << m_xTilder[i][j];
    //     }
    // 	std::cout << std::endl;
    // }
    // exit(0);
}

/*** Box in lab reference points  *****************************/
void KizyunTen::setKizyun3D_box(){
    GlobalConstants gc;
  
    // Set world coordinates
    m_He[0] = 0.0;
    m_He[1] = 320.0;
    m_He[2] = 0.0; 

    m_Hw[0] = 337.0;
    m_Hw[1] = m_He[1];
    m_Hw[2] = m_He[2];

    m_Le[0] = 0.0;
    m_Le[1] = 0.0;
    m_Le[2] = 0.0;

    m_Lw[0] = m_Hw[0];
    m_Lw[1] = 0.0;
    m_Lw[2] = 0.0;

    m_Phw[0] = m_Hw[0];
    m_Phw[1] = m_Hw[1];
    m_Phw[2] = 379.0;
	
    m_Phe[0] = 0.0;
    m_Phe[1] = m_Phw[1];
    m_Phe[2] = m_Phw[2];
	
    m_Ple[0] = 0.0;
    m_Ple[1] = 0.0;
    m_Ple[2] = m_Phe[2];
	
    readTxtFile(m_txtFilePath);

}

/*** Ippan array (410) reference points  *****************************/
void KizyunTen::setKizyun3D_array(){
    GlobalConstants gc;
  
    // std::cout << "ABCD=" << m_ABCD << std::endl;  /////////////////////t

    // Set world coordinates
    if(m_ABCD == 'A'){
        m_He[0] = 0.0;
        m_He[1] = gc.hL + (gc.mY*4+gc.deltaY*3)*std::sin(gc.theta_p*gc.pi/180.0);
        m_He[2] = (gc.mY*4+gc.deltaY*3)*std::cos(gc.theta_p*gc.pi/180.0); 

        m_Hw[0] = gc.mX;
        m_Hw[1] = m_He[1];
        m_Hw[2] = m_He[2];

        m_Me[0] = 0.0;
        // m_Me[0] = m_MHe[0];
        // m_Me[0] = m_Hce[0];   // Be carefull !!!
        m_Me[1] = gc.hL + gc.mY*std::sin(gc.theta_p*gc.pi/180.0);
        m_Me[2] = gc.mY*std::cos(gc.theta_p*gc.pi/180.0);

        m_Mw[0] = m_Hw[0];
        // m_Mw[0] = m_MHw[0];
        // m_Mw[0] = m_Hcw[0];   // Be carefull !!!
        m_Mw[1] = m_Me[1];
        m_Mw[2] = m_Me[2];
    
        m_Le[0] = 0.0;
        m_Le[1] = gc.hL;
        m_Le[2] = 0.0;

        m_Lw[0] = gc.mX;
        m_Lw[1] = m_Le[1];
        m_Lw[2] = m_Le[2];
    }

    readTxtFile(m_txtFilePath);
}

/*** Read image coordinates ***/
void KizyunTen::readTxtFile(std::string txtFileName){
    std::cout << "Read file '" << txtFileName << "'\n";
    
    std::fstream ifs {txtFileName, std::ios_base::in};
    if(!ifs){
        std::cout << "couldn't open '" << txtFileName << "' for reading" << std::endl;
        exit(-1);
    }

    int kizNum {0};   // number of image coordinates
    std::string str;
    while(std::getline(ifs, str)){
	if(str[0]=='/' && str[1]=='/') continue;  // comment out lines
	
        std::vector<std::string> vStr = splitString(str, ',');

        for(int i=0; i<vStr.size(); i++){
	    if(vStr[i]=="Width"){
                m_width = std::stoi(vStr[i+1]);  // image width
	    }else if(vStr[i]=="Height"){
                m_height = std::stoi(vStr[i+1]);  // image height
            }else if(vStr[i]=="Hw"){
                m_hw[0] = std::stoi(vStr[i+1]);  // x 
                m_hw[1] = std::stoi(vStr[i+2]);  // y
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_Hw[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_hw[k]);
		m_kizName[kizNum] = "Hw";
		kizNum++;
            }else if(vStr[i]=="He"){
                m_he[0] = std::stoi(vStr[i+1]);  // x 
                m_he[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_He[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_he[k]);
		m_kizName[kizNum] = "He";
		kizNum++;
            }else if(vStr[i]=="Mw"){
                m_mw[0] = std::stoi(vStr[i+1]);  // x 
                m_mw[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_Mw[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_mw[k]);
		m_kizName[kizNum] = "Mw";
		kizNum++;
            }else if(vStr[i]=="Me"){
                m_me[0] = std::stoi(vStr[i+1]);  // x 
                m_me[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_Me[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_me[k]);
		m_kizName[kizNum] = "Me";
		kizNum++;
            }else if(vStr[i]=="Lw"){
                m_lw[0] = std::stoi(vStr[i+1]);  // x 
                m_lw[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_Lw[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_lw[k]);
		m_kizName[kizNum] = "Lw";
		kizNum++;
            }else if(vStr[i]=="Le"){
                m_le[0] = std::stoi(vStr[i+1]);  // x 
                m_le[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_Le[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_le[k]);
		m_kizName[kizNum] = "Le";
		kizNum++;
            }else if(vStr[i]=="MHw"){
                m_mhw[0] = std::stoi(vStr[i+1]);  // x 
                m_mhw[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_MHw[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_mhw[k]);
		m_kizName[kizNum] = "MHw";
		kizNum++;
            }else if(vStr[i]=="MHe"){
                m_mhe[0] = std::stoi(vStr[i+1]);  // x 
                m_mhe[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_MHe[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_mhe[k]);
		m_kizName[kizNum] = "MHe";
		kizNum++;
            }else if(vStr[i]=="Hcw"){
                m_hcw[0] = std::stoi(vStr[i+1]);  // x 
                m_hcw[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_Hcw[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_hcw[k]);
		m_kizName[kizNum] = "Hcw";
		kizNum++;
            }else if(vStr[i]=="Hce"){
                m_hce[0] = std::stoi(vStr[i+1]);  // x 
                m_hce[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_Hce[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_hce[k]);
		m_kizName[kizNum] = "Hce";
		kizNum++;
            }else if(vStr[i]=="Lcw"){
                m_lcw[0] = std::stoi(vStr[i+1]);  // x 
                m_lcw[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_Lcw[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_lcw[k]);
		m_kizName[kizNum] = "Lcw";
		kizNum++;
            }else if(vStr[i]=="Lce"){
                m_lce[0] = std::stoi(vStr[i+1]);  // x 
                m_lce[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_Lce[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_lce[k]);
		m_kizName[kizNum] = "Lce";
		kizNum++;
            }else if(vStr[i]=="MHcw"){
                m_mhcw[0] = std::stoi(vStr[i+1]);  // x 
                m_mhcw[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_MHcw[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_mhcw[k]);
		m_kizName[kizNum] = "MHcw";
		kizNum++;
            }else if(vStr[i]=="MHce"){
                m_mhce[0] = std::stoi(vStr[i+1]);  // x 
                m_mhce[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_MHce[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_mhce[k]);
		m_kizName[kizNum] = "MHce";
		kizNum++;
            }else if(vStr[i]=="Mcw"){
                m_mcw[0] = std::stoi(vStr[i+1]);  // x 
                m_mcw[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_Mcw[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_mcw[k]);
		m_kizName[kizNum] = "Mcw";
		kizNum++;
            }else if(vStr[i]=="Mce"){
                m_mce[0] = std::stoi(vStr[i+1]);  // x 
                m_mce[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_Mce[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_mce[k]);
		m_kizName[kizNum] = "Mce";
		kizNum++;
            }else if(vStr[i]=="pHe" || vStr[i]=="Phe"){
                m_phe[0] = std::stoi(vStr[i+1]);  // x 
                m_phe[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_Phe[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_phe[k]);
		m_kizName[kizNum] = "Phe";
		kizNum++;
            }else if(vStr[i]=="pLe" || vStr[i]=="Ple" ){
                m_ple[0] = std::stoi(vStr[i+1]);  // x 
                m_ple[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_Ple[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_ple[k]);
		m_kizName[kizNum] = "Ple";
		kizNum++;
            }else if(vStr[i]=="pHw" || vStr[i]=="Phw"){
                m_phw[0] = std::stoi(vStr[i+1]);  // x 
                m_phw[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_Phw[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_phw[k]);
		m_kizName[kizNum] = "Phw";
		kizNum++;
            }else if(vStr[i]=="pLw" || vStr[i]=="Plw" ){
                m_plw[0] = std::stoi(vStr[i+1]);  // x 
                m_plw[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_Plw[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_plw[k]);
		m_kizName[kizNum] = "Plw";
		kizNum++;
            }else if(vStr[i]=="mH"){
                m_mh[0] = std::stoi(vStr[i+1]);  // x 
                m_mh[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_Mh[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_mh[k]);
		m_kizName[kizNum] = "mH";
		kizNum++;
            }else if(vStr[i]=="bpw" ){
                m_bpw[0] = std::stoi(vStr[i+1]);  // x 
                m_bpw[1] = std::stoi(vStr[i+2]);  // y 
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_Bpw[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_bpw[k]);
		m_kizName[kizNum] = "Bpw";
		kizNum++;
            }else if(vStr[i]=="bpe"){
                m_bpe[0] = std::stoi(vStr[i+1]);  // x 
                m_bpe[1] = std::stoi(vStr[i+2]);  // y
		for(int k=0; k<3; k++) m_Xtilder[kizNum][k] = m_Bpe[k];
 		for(int k=0; k<2; k++) m_xTilder[kizNum][k] = static_cast<double>(m_bpe[k]);
		m_kizName[kizNum] = "Bpe";
		kizNum++;
            }
        }

	if(kizNum >= m_n) break;
    }
    // std::cout << "hw.x=" << m_hw[0] << " hw.y=" << m_hw[1] << std::endl;////////////////////
    // std::cout << "he.x=" << m_he[0] << " he.y=" << m_he[1] << std::endl;////////////////////
    // std::cout << "mw.x=" << m_mw[0] << " mw.y=" << m_mw[1] << std::endl;////////////////////
    // std::cout << "me.x=" << m_me[0] << " me.y=" << m_me[1] << std::endl;////////////////////
    // std::cout << "lw.x=" << m_lw[0] << " lw.y=" << m_lw[1] << std::endl;////////////////////
    // std::cout << "le.x=" << m_le[0] << " le.y=" << m_le[1] << std::endl;////////////////////
    // std::cout << "mhw.x=" << m_mhw[0] << " mhw.y=" << m_mhw[1] << std::endl;////////////////////
    // std::cout << "mhe.x=" << m_mhe[0] << " mhe.y=" << m_mhe[1] << std::endl;////////////////////
    // std::cout << "hcw.x=" << m_hcw[0] << " hcw.y=" << m_hcw[1] << std::endl;////////////////////
    // std::cout << "hce.x=" << m_hce[0] << " hce.y=" << m_hce[1] << std::endl;////////////////////
    // std::cout << "lcw.x=" << m_lcw[0] << " lcw.y=" << m_lcw[1] << std::endl;////////////////////
    // std::cout << "lce.x=" << m_lce[0] << " lce.y=" << m_lce[1] << std::endl;////////////////////
    // std::cout << "pt.x=" << m_pt[0] << " pt.y=" << m_pt[1] << std::endl;////////////////////
    // std::cout << "mfpt.x=" << m_mfpt[0] << " mfpt.y=" << m_mfpt[1] << std::endl;////////////////////
    // std::cout << "mfpb.x=" << m_mfpb[0] << " mfpb.y=" << m_mfpb[1] << std::endl;////////////////////
    // exit(0);///////////////////////
}

/********   Split string   ***************/
std::vector<std::string> KizyunTen::splitString(std::string str, char c){
    std::vector<std::string> ret;

    unsigned int i0 = 0;
    bool final = false;
    
    while(true){
        int i1 = str.find_first_of(c, i0);

        if(i0==i1){   // If the first character is space
            i0 = i1+1;
            continue;
        }else if(i1 == std::string::npos){
            i1 = str.size();
            final = true;
        }
	
        std::string subStr {str, i0, i1-i0};
        ret.push_back(subStr);
	
        // std::cout << "i0= " << i0 << " i1= " << i1 << " " << subStr << std::endl;

        i0 = i1+1;

        if(final) break;
    }

    return ret;
}

/** for LM algorithm. f[2*m_n] ***/
void KizyunTen::f(VecDoub_I &P, VecDoub_O &fval){
    for(int k=0; k<m_n; k++){
        double P3Xk {0.0};
        for(int m=0; m<3; m++){
            P3Xk += P[8+m]*m_Xtilder[k][m];
        }
        P3Xk += P[11];  // m=4

        for(int l=0; l<2; l++){
            int l1 = (l==0)?0:4;
            double PlXk {0.0};
            for(int m=0; m<3; m++){
                PlXk += P[l1+m]*m_Xtilder[k][m];
                // printf("l1+m=%d, l=%d, k=%d, m=%d\n", l1+m, l, k, m); /////////t
            }
            PlXk += P[l1+3];  // m=4

            fval[2*k +l] = PlXk/P3Xk;
            // printf("fval[%d]=%f\n", 2*k +l, fval[2*k +l]);////////////////t
        }
    }
}

/** for LM algorithm. f[2*m_n] ***/
void KizyunTen::fOrig(VecDoub_I &P, VecDoub_O &fval){
    for(int k=0; k<m_n; k++){
        double P3Xk {0.0};
        for(int m=0; m<3; m++){
            P3Xk += P[8+m]*m_XtilderOrig[k][m];
        }
        P3Xk += P[11];  // m=4

        for(int l=0; l<2; l++){
            int l1 = (l==0)?0:4;
            double PlXk {0.0};
            for(int m=0; m<3; m++){
                PlXk += P[l1+m]*m_XtilderOrig[k][m];
                // printf("l1+m=%d, l=%d, k=%d, m=%d\n", l1+m, l, k, m); /////////t
            }
            PlXk += P[l1+3];  // m=4

            fval[2*k +l] = PlXk/P3Xk;
            // printf("fval[%d]=%f\n", 2*k +l, fval[2*k +l]);////////////////t
        }
    }
}

/** for LM algorithm. Calculate Jacobian. ***/
void KizyunTen::J(VecDoub_I &P, MatDoub_O &Jval){
    VecDoub P2X(m_n);
    for(int k=0; k<m_n; k++){
        double sum {0.0};
        for(int m=0; m<3; m++){
            sum += P[8+m]*m_Xtilder[k][m];
        }
        P2X[k] = sum + P[11];
    }
    
    for(int i=0; i<2*m_n; i++){
        int k {i/2};
        int l {i%2};
        for(int j=0; j<12; j++){
            int a {j/4};
            int b {j%4};
            if(a<2){
                if(l==a){
                    if(b<3) Jval[i][j] = m_Xtilder[k][b]/P2X[k];
                    else Jval[i][j] = 1.0/P2X[k];
                }else{
                    Jval[i][j] = 0.0;
                }
            }else{  // a==2
                double PlXk {0.0};
                for(int m=0; m<3; m++){
                    PlXk += P[4*l +m]*m_Xtilder[k][m];
                }
                PlXk += P[4*l+3];
		
                if(b<3) Jval[i][j] = -PlXk/(P2X[k]*P2X[k])*m_Xtilder[k][b];
                else Jval[i][j] = -PlXk/(P2X[k]*P2X[k]);
            }
        }
    }
}

/** for LM algorithm. Calculate Jacobian numerically. ***/
void KizyunTen::Jnum(VecDoub &P, MatDoub_O &Jval){
    VecDoub fval0(2*m_n+1);
    VecDoub fval1(2*m_n+1);
    
    f(P, fval0);
	    
    for(int j=1; j<=12; j++){
        double delta {std::max(std::abs(1.0E-4*P[j]), 1.0E-6)};

        double bak {P[j]};
	
        P[j] += delta;
        f(P, fval1);
	    
        for(int i=1; i<=2*m_n; i++){
            Jval[i][j] = (fval1[i] - fval0[i])/delta;
        }

        P[j] = bak;
    }
}

/**** Normalize image coordinates m[0..1], m[2]=1.0 not given  *****/
void KizyunTen::normalize_m(double m[2]){
    for(int i=0; i<2; i++){
        m[i] = m[i]*m_Ts + m_Tt[i];
    }
}

/***** x[1..], y[1..], u[1..], v[1..], T[0..][0..] ****************/
void KizyunTen::normalize(const double x[], const double y[], VecDoub_O &u, VecDoub_O &v,
			  MatDoub &T){
    int n {u.size()-1};
    
    double mx {0.0};
    double my {0.0};
    for(int i=1; i<=n; i++){
	mx += x[i];
	my += y[i];
    }
    mx /= n;  my /= n;

    double vx {0.0};
    double vy {0.0};
    for(int i=1; i<=n; i++){
	vx += (x[i]-mx)*(x[i]-mx);
	vy += (y[i]-my)*(y[i]-my);
    }
    double ax {1.0/std::sqrt(vx/n)};
    double ay {1.0/std::sqrt(vy/n)};

    for(int i=1; i<=n; i++){
	u[i] = ax*(x[i]-mx);
	v[i] = ay*(y[i]-my);
    }
    
    T[0][0] = ax; T[0][1] = 0.0; T[0][2] = -ax*mx;
    T[1][0] = 0.0; T[1][1] = ay; T[1][2] = -ay*my;
    T[2][0] = 0.0; T[2][1] = 0.0; T[2][2] = 1.0;
}

/****** Denormalize image coordinates  x[0..2]   *****/
void KizyunTen::denormalize_x(double x[3]){
    for(int i=0; i<2; i++){
        x[i] = (x[i]-m_Tt[i]*x[2])/m_Ts;
    }
}

/******** Denormalize world coordinates,  M[0...3]    ******/
void KizyunTen::denormalize_M(VecDoub_IO &M){
    for(int i=0; i<=2; i++){
        // M[i] = (M[i]-m_Ut[i+1]*M[3])/m_Us;
        M[i] = (M[i]-m_Ut[i]*M[3])/m_Us;  // modify 2026.1.5
    }
}

/***** Calculate the normalized camera matrix m_Pbar *****/
void KizyunTen::calcPbar(){
    DLT dlt(*this);
    dlt.calcP(m_Pbar);

    // for check. Note: m_Pbar is overwritten. Note!!!!!
    VecDoub C3(3);
    MatDoub P3inv(3, 3);
    dlt.denormalize(m_Pbar, C3, P3inv);
}

/*****************************************************/
void KizyunTen::get_xTilderOrig(VecDoub_O &xTilderOrig)const{
    for(int i=0; i<m_n; i++){
	xTilderOrig[i*2]=m_xTilderOrig[i][0];
	xTilderOrig[i*2+1]=m_xTilderOrig[i][1];
    }
}

/*** for LM algorithm.  X[2*m_n] ***/
void KizyunTen::get_xTilder(VecDoub_O &X)const{
    for(int k=0; k<X.size(); k++){
        int i = k/2;
        int j = k -i*2;

        X[k] = m_xTilder[i][j];
        // printf("k=%d, i=%d, j=%d, X=%f\n", k, i, j, X[k]); ////////////////////
    }
}

/******** 2026.1.5  *****************************************/
void KizyunTen::normalize_P(MatDoub_I &P, MatDoub_O &P_nm){
    MatDoub A(3, 4, 0.0);
    for(int i=0; i<3; i++){
	for(int j=0; j<4; j++){
	    for(int k=0; k<3; k++) A[i][j] += m_T[i][k]*P[k][j];
	}
    }
    
    for(int i=0; i<3; i++){
	for(int j=0; j<4; j++){
	    P_nm[i][j] = 0.0;
	    for(int k=0; k<4; k++) P_nm[i][j] += A[i][k]*m_Uinv[k][j];
	}
    }
}

/******** 2026.1.6  *****************************************/
void KizyunTen::denormalize_P(MatDoub_I &P_nm, MatDoub_O &P){
    MatDoub A(3, 4, 0.0);
    for(int i=0; i<3; i++){
	for(int j=0; j<4; j++){
	    for(int k=0; k<3; k++) A[i][j] += m_Tinv[i][k]*P_nm[k][j];
	}
    }
    
    for(int i=0; i<3; i++){
	for(int j=0; j<4; j++){
	    P[i][j] = 0.0;
	    for(int k=0; k<4; k++) P[i][j] += A[i][k]*m_U[k][j];
	}
    }
}
