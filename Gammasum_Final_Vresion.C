#include "TRint.h"
#include <string.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TObject.h"
#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <TVector3.h>
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TSpectrum.h"
#include "TH1I.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TPaveStats.h"
#include "TVector3.h"
#include <stdlib.h>
#include <ctime>
#include <time.h>
#include <iostream>
#define N 5 //  N : NUMBER OF LEVELS TO CONSIDER ( N=0 -> G.S. )//EUA =14, EUB =15 , Ba = 5
using namespace std;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PROGRAM GAMMASUM
//   PROGRAM TO CALCULATE COINCIDENCE-SUMMING CORRECTIONS FOR Ge GAMMA
//   DETECTOR SPECTRA.  SPECIFICALLY, THE PHOTO-PEAK EFFICIENCY IS
//   COMPUTED GIVEN THE OBSERVED INTENSITIES, THE SOURCE STRENGTH,
//   THE TOTAL EFFICIENCY, AND THE DECAY SCHEME.  THE TOTAL EFFICIENCY
//   MAY BE GIVEN EITHER AS THE TOTAL EFFICIENCY (ETA), OR THE PEAK-TO-
//   TOTAL RATIO (WHICH IS DIRECTLY MEASUREABLE EXPERIMENTALLY).
//   THE EFFECT OF ANGULAR CORRELATIONS IS IGNORED IN THE PROCEDURE.
//   THE DECAY SCHEME IS ARBITRARY, AND IS SPECIFIED BY AN INPUT FILE.

//  TECHNIQUE IS BASED ON MATRIX ALGEBRA; MOST OF ALGORYTHMS IN
//  THIS PROGRAM ARE FROM A PAPER BY T. M. SEMKOW, ET AL.,
//  NUCLEAR INSTRUMENTS AND METHODS A290 (1990) 437-444

//  THE FILE GAMMASUM.INP IS READ TO START THINGS OFF. IT CONTAINS:
//  F_LEV : THE FILE THAT CONTAINS THE DECAY SCHEME
//  F_DAT : THE FILE THAT CONTAINS THE OBSERVED INTENSITIES & SOURCE STRENGTH
//  F_OUT : THE OUTPUT FILE FOR THE CORRECTED EFFICIENCIES
//  F_RAW : THE OUTPUT FILE FOR THE UNCORRECTED EFFICIENCIES
//  F_CHK : THE OUTPUT FILE FOR THE DEBUGGING INFO
//   MODE  : MODE=1 -> USE TOTAL EFFICIENCY
//            MODE=2 -> USE PEAK-TO-TOAL RATIO
//   NGEO  : GEOMETRY CODE
//
// DESCRIPTION OF SOME IMPORTANT VARIABLES:
//
//   F(I)    :  THE FEEDING FRACTION OF THE NUCLEAR STATES (I)
//              f IN T. M. SEMKOW ET AL. PAPER
//   X(J,I)  :  THE FRACTION OF TIMES STATE J DECAYS TO STATE I
//              Xji IN T. M. SEMKOW ET AL. PAPER
//   A(J,I)  :  MATRIX aji IN T. M. SEMKOW ET AL. PAPER
//   E(J,I)  :  MATRIX eji IN T. M. SEMKOW ET AL. PAPER
//   B(J,I)  :  MATRIX bji IN T. M. SEMKOW ET AL. PAPER
//   A1(J,I) :  MATRIX Aji IN T. M. SEMKOW ET AL. PAPER
//   B1(J,I) :  MATRIX Bji IN T. M. SEMKOW ET AL. PAPER
//   QN(J,I) :  MATRIX Nji IN T. M. SEMKOW ET AL. PAPER
//   QM(J,I) :  MATRIX Mji IN T. M. SEMKOW ET AL. PAPER
//   EN(I)   :  THE ENERGIES (IN MeV) OF THE LEVELS (I)
//   C(J,I)  :  MATRIX Cji IN T. M. SEMKOW ET AL. PAPER
//   TMP1(J,I), TMP2(J,I) :  TEMPORARY MATRICES
//   FB(I)   :  THE ROW MATRIX f*B IN T. M. SEMKOW ET AL. PAPER
//   QI(J,I) :  THE KNOWN GAMMA RAY INTENSITIES, NORAMLIZED TO 100 DECAYS
//              OF PARENT.  THIS INFORMATION IS IN PRINCIPLE CONTAINED IN
//              THE DECAY SCHEME, BUT THIS QUANTITY IS OFTEN KNOWN TO GREATER
//              PRECISION EXPERIMENTALLY.
//   S(J,I)  :  OBSERVED GAMMA INTENSITY (IN COUNTS) OF TRANSITION J -> I
//   DS(J,I) :  ERROR ASSOCIATED WITH S(J,I)
//   B0(J,I) :  MATRIX B(0)ji IN T. M. SEMKOW ET AL. PAPER
//   EP(J,I) :  PHOTO-PEAK EFFICIENCY FOR TRANSITION J -> I
//   ET(J,I) :  TOTAL EFFICIENCY FOR TRANSITION J -> I
//   QN0(J,I) : THE MATRIX N(0)ji IN T. M. SEMKOW ET AL. PAPER
//   QN0A(J,I) : THE MATRIX N(0)*a IN T. M. SEMKOW ET AL. PAPER
//   FB0(I)  :  THE ROW MATRIX f*B(0) IN T. M. SEMKOW ET AL. PAPER
//   D(J,I)  :  THE MATRIX Dji IN T. M. SEMKOW ET AL. PAPER
//   QN0C(J,I) : THE MATRIX N(0)*c IN T. M. SEMKOW ET AL. PAPER
//   AL(J,I) :  THE GAMMA-RAY CONVERSION COEFFICIENTS OF TRANSION J->I
//              ALPHA IN T. M. SEMKOW ET AL. PAPER
//   QIE(J,I) : PERCENT ERRORS ASSOCIATED WITH QI(J,I)

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//  PHOTO-FRACTION GOES HERE
double PHI(double EG, double NGEO);

// YOUR TOTAL EFFICIENCY GOES HERE
double ETA(double EG, double NGEO);

// PHOTO-EFFICIENCY FUNCTION GOES HERE
double  EPS(double EG, double NGEO);

// PRINT OUT MATRIX A TO SCREEN; USED FOR DE-BUGGING
void Print(double Matrix[][N], int row, int column);

// PRINT OUT vector A TO SCREEN; USED FOR DE-BUGGING
void display(double mult[N], int rowThird);

// SET MATRIX A EQUAL TO THE IDENTITY MATRI
void Identity ( double Matrix[][N], int row, int column);

// ADD MATRIX A INTO MATRIX B
void Subtract(double firstMatrix[][N], double secondMatrix[][N], double thirdMatrix[][N], int row, int column);

// Multiplying matrix firstMatrix and secondMatrix and storing in array mult.
void MMmultiplyMatrices(double firstMatrix[][N], double secondMatrix[][N], double mult[][N], int rowFirst, int columnFirst, int rowSecond, int columnSecond);

// Multiplying 1D array and 2D array and storing in array mult.
void MmultiplyMatrices(double firstMatrix[N], double secondMatrix[][N], double mult[N], int rowFirst, int rowSecond , int columnSecond);

//ADD MATRIX A INTO MATRIX B
void Addition(double firstMatrix[][N], double secondMatrix[][N], double thirdMatrix[][N], int row, int column);

//PUT MATRIX A INTO MATRIX B
void Put (double firstMatrix[][N], double secondMatrix[][N], int row, int column);

// Initializing elements of matrix A to 0.
void Zero( double Matrix[][N], int row, int column);

// NUMBER OF PARENT NUCLEI DECAYS.
double NUMBER_OF_PARENT_NUCLEI_DECAYS(double Initial_Activity, double time, double half_life, double T);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                               The main program                                                            //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MODE = 2 : MODE = 1 USE TOTAL EFFICIENCY WHILE MODE = 2 USE PHOTO-FRACTION.
// NGEO = 1 : GEOMETRY CODE for close = 1 WHILE for far = 2
//  FEED : FRACTION OF 100 DECAYS OF PARENT LEADING TO THIS GROUP OF STATES. USUALLY 100, BUT SOME NUCLEI (EG. 152 Eu) DECAY INTO TWO DIFFERENT FINAL NUCLEI. (e capture) = 72.08 WHILE BETA MINUS = 27.92
// ================================================================
// The following variables are changing based on the experiment.
// Initial_Activity = 3933.1;// in Bq A0= 0.1063 MicoCi
// time = 520560000; // in Sec -----  May,1,2003 - Oct,29,2019.
// half_life = 427186459;// in Sec ----- 13.516 years
// T = 32202;  // time of the run
// ============================================================
void gammasum(Double_t MODE, Double_t NGEO, Double_t FEED, Double_t Initial_Activity, Double_t time, Double_t half_life, Double_t T ){
    if (MODE == 1){
    cout <<"MODE = 1 USE TOTAL EFFICIENCY"<<endl;
    }
    else if  (MODE == 2){
    cout <<"MODE = 2 USE PHOTO-FRACTION"<<endl;
    }
    else{
    cout <<"MODE ERROR !!!!"<<endl;
    }
    cout << "GEOMETRY CODE = " << NGEO << endl;
    cout <<NUMBER_OF_PARENT_NUCLEI_DECAYS(Initial_Activity,time,half_life,T)<< endl;
double F[N],X[N][N],A[N][N],E[N][N],B[N][N];
double A1[N][N],B1[N][N],QN[N][N],QM[N][N];
double EN[N],C[N][N],TMP1[N][N],TMP2[N][N];
double FB[N],QI[N][N],S[N][N],DS[N][N];
double B0[N][N],EP[N][N],ET[N][N];
double QN0[N][N],QN0A[N][N],FB0[N];
double D[N][N],QN0C[N][N],AL[N][N],QIE[N][N];
cout <<"PHOTO-PEAK EFFICIENCY CALCULATION INCLUDING SUMMING"<<endl;
ifstream EUA01LEV;
EUA01LEV.open(Form("Ba_Elev_1.txt"));// needs to be changed
for (int i =0;i < N ;i++){
EUA01LEV>>EN[i]>>F[i];
}
EUA01LEV.close();
    
int I,J;
double QT1,QT2,QT3,QT4,QS,DQS;
Zero(QI,N,N);
Zero(QIE,N,N);
Zero(X,N,N);
Zero(AL,N,N);
ifstream EUA02LEV;
EUA02LEV.open(Form("Ba_Elev_2.txt")); // needs to be changed
Int_t i=0;
while(!EUA02LEV.eof())
{
EUA02LEV>>J>>I>>QT1>>QT2>>QT3>>QT4;
QI[J][I] = QT1;
QIE[J][I] = QT2;
X[J][I] = QT3;
AL[J][I] = QT4;
i++;
}
EUA02LEV.close();
    
// RENORMALIZE FEEDING FRACTIONS TO THIS SET OF STATES
for (int i =0;i < N;i++){
F[i] = F[i]*( 100.0 / FEED );
}

// CALCULATE MATRIX C
Zero(C,N,N);
    
// Initializing elements of matrix A to 0.
for(int i = 0; i < N-1; ++i){
    for(int j = i+1 ; j < N ; ++j){
        C[j][i]= X[j][i] / ( AL[j][i] + 1.0 ) ;
    }
}
    

// READ THE INPUT INTENSITIES
// R : THE NUMBER OF PARENT NUCLEI DECAYS DURING LIVE TIME REQUIRED TO
//     OBTAIN SPECTRUM
// S(J,I)  : OBSERVED GAMMA INTENSITY (IN COUNTS) OF TRANSITION J -> I
// DS(J,I) : ERROR ASSOCIATED WITH S(J,I)
Zero(S,N,N);
Zero(DS,N,N);
ifstream EUAINP;
EUAINP.open(Form("Ba_input.txt"));// needs to be changed
int k= 0;
while(!EUAINP.eof()){
EUAINP>>J>>I>>QS>>DQS;
S[J][I] = QS;
DS[J][I] = DQS;
k++;
}
EUAINP.close();

// CALL ZERO(X) SETS MATRIX X = 0
// EP : PHOTO - PEAK EFFICIENCY
// ET : TOTAL EFFICIENCY
Zero(EP,N,N);
Zero(ET,N,N);

double EN1;
    
// GO THROUGH INPUT INTENSITIES; CALCULATE THE 0TH ORDER
// PHOTO-PEAK EFFICIENCIES; USE THE EFFICIENCY FUNCTION
// SPECIFIED BY NGEO IF THE OBSERVED INTENSITY ISN'T PROVIDED
    
for(int i = 0; i < N-1; ++i){
    for(int j = i+1; j < N ; ++j){
        if(S[j][i] < 1.0){
            EN1 = EN[j] - EN[i];
            if(EN1 > 0.0){
                EP[j][i] = EPS( EN1 , NGEO );
            }
            else{
                EP[j][i] = 0.0;
            }
        }
        else{
            EP[j][i] = S[j][i] / ( NUMBER_OF_PARENT_NUCLEI_DECAYS(Initial_Activity,time,half_life,T) * QI[j][i] / 100.0 );
        }
    }
}
    
//        CALCULATE B0,FB0,QN0,QN0C PER T. M. SEMKOW ET AL.
//        THESE QUANTITIES DO NOT DEPEND ON THE TOTAL EFFICIENCY OR
//        THE PHOTO-PEAK EFFICIENCY.
//
//        Identity     : SET EQUAL TO IDENTITY MATRIX
//        PUT          : PUT THE FIRST MATRIX IN THE SECOND MATRIX
//        MMmultiplyMatrices(A,B,C)  : C = A * B   FOR MATRICES A,B,C
//        Addition(A,B,C)  : C = A + B   FOR MATRICES A,B,C
//        MmultiplyMatrices(A,B,C) : C = A * B   LIKE MULT, BUT A IS A ROW MATRIX
    
Identity(B0,N,N);
Identity(TMP1,N,N);
for(int i = 1 ; i < N ; ++i){
    MMmultiplyMatrices(X,TMP1,TMP2,N,N,N,N);
    Put(TMP2,TMP1,N,N);
    Addition(B0,TMP1,TMP2,N,N);
    Put(TMP2,B0,N,N);
}
    
MmultiplyMatrices(F,B0,FB0,N,N,N);
Zero(QN0,N,N);
Zero(QN0C,N,N);
for(int i = 0 ; i < N ; ++i){
    QN0[i][i] = FB0[i];
}
MMmultiplyMatrices(QN0,C,QN0C,N,N,N,N);

///////////////////////////////////////////////////////////////////////////////////////////////
// PRINT TO THE SCREEN :
// THE GAMMA RAY ENERGY
// THE KNOWN INTENSITY FROM QI
// THE ERROR IN THE KNOWN INTENSITY
// THE KNOWN INTENSITY FROM THEN LEVEL SCHEME IN QN0C
// WRITE TO THE FILE F_CHK :
// THE GAMMA RAY ENERGY
// THE KNOW INTENSITY FROM QI
// THE ERROR IN THE KNOWN INTENSITY
// THE KNOW INTENSITY FOR THE LEVEL SCHEME IN QN0C
// WRITE TO THE FILE F_RAW :
// THE GAMMA ENERGY
// THE RAW (I. E. UNCORRECTED FOR SUMMING) EFFICIENCY
// THE ERROR IN EFFICIENCY

 ofstream F_RAW;
 F_RAW.open(Form("/home/gula/Desktop/8MeV_analysis/F_RAW.txt"));
 ofstream F_CHK;
 F_CHK.open(Form("/home/gula/Desktop/8MeV_analysis/F_CHK.txt"));
 for(int i = 0 ; i < N-1 ; ++i){
     for(int j = i+1 ; j < N ; ++j){
         F_CHK << j <<"  "<< i <<"  "<< EN[j]-EN[i] <<"  "<< QI[j][i] <<"  "<< QIE[j][i] <<"  " << QN0C[j][i]*FEED<<endl;
             if (S[j][i] >= 0.5){
                 F_RAW<<EN[j]-EN[i]<<" "<< EP[j][i] <<""<< EP[j][i]*DS[j][i]/S[j][i] <<endl;
         }
     }
 }
 F_RAW.close();
 F_CHK.close();

//  THE ITERATIVE LOOP COMES BACK TO THIS POINT.
//  ALL QUANTITIES THAT DEPEND UPON THE PHOTO-EFFICIENCY MUST
//  BE CALCULATED IN THIS LOOP.  ALL QUANTITIES THAT DEPEND UPON
//  THE TOTAL EFFICIENCY ARE ALSO CALCULATED IN THIS LOOP BECAUSE
//  THEY DEPEND OF THE PHOTO-EFFICIENCY IF MODE=2.
//  SET ITERATION NUMBER EQUAL TO 1
//  cout <<" ITERATION NUMBER "<<"      "<< " CONVERGENCE " <<endl;
//  cout << endl;
double DEFF;
int NIT = 1;
do{
    
//  CALCULATE ETA, E, AND B
//  THESE QUANTITIES CLEARLY DEPEND ON THE TOTAL EFFICIENCY
//  ETA DEPENDS ON THE PHOTO-EFFICIENCY IF MODE=2.
    
Zero(E,N,N);
Zero(B,N,N);
for(int i = 0 ; i < N-1 ; ++i){
    for(int j = i+1 ; j < N ; ++j){
        if (MODE == 1){
            ET[j][i] = ETA( EN[j] - EN[i] , NGEO );
        }
        else if (MODE == 2){
            ET[j][i] = EP[j][i] / PHI( EN[j] - EN[i] , NGEO );
        }
        E[j][i] = C[j][i] * ET[j][i];
        B[j][i] = X[j][i] - E[j][i];
    }
}
    
//   CALCULATE B1,QN,QM
//   THESE QUANTITIES DO DEPEND ON THE TOTAL EFFICIENCY
    
    
Identity(B1,N,N);
Identity(TMP1,N,N);
for (int i =1 ; i < N ; i++){
    MMmultiplyMatrices(B,TMP1,TMP2,N,N,N,N);
    Put(TMP2,TMP1,N,N);
    Addition(B1,TMP1,TMP2,N,N);
    Put(TMP2,B1,N,N);
    }

MmultiplyMatrices(F,B1,FB,N,N,N);
Zero(QN,N,N);
Zero(QM,N,N);
for(int i = 0 ; i < N ; ++i){
    QN[i][i] = FB[i];
    QM[i][i] = B1[i][0];
    }
   
//  CALCULATE A
Zero(A,N,N);
for(int i = 0 ; i < N-1 ; ++i){
    for(int j = i+1 ; j < N ; ++j){
        A[j][i] = C[j][i] * EP[j][i];
    }
}

//  CALCULATE A1
Zero(A1,N,N);
Identity(TMP1,N,N);
for (int i =1 ; i < N ; i++){
    MMmultiplyMatrices(A,TMP1,TMP2,N,N,N,N);
    Put(TMP2,TMP1,N,N);
    Addition(A1,TMP1,TMP2,N,N);
    Put(TMP2,A1,N,N);
}

//  CALCULATE D

MMmultiplyMatrices(A1,QM,TMP1,N,N,N,N);
MMmultiplyMatrices(QN,TMP1,TMP2,N,N,N,N);
MMmultiplyMatrices(QN0,A,QN0A,N,N,N,N);
Subtract(TMP2,QN0A,D,N,N);
    
//  CALCULATE THE NEW PHOTO-PEAK EFFICIENCIES.
//  ALSO CALCULATE THE "CONVERGENCE" DEFINED AS THE RMS VALUE
//  OF THE FRACTIONAL CHANGE IN EFFICIENCY FROM THE PRECEDING ITERATION.
    
int NEFF = 0;
double SEFF= 0.0;
double EFF1;
for(int i = 0 ; i < N-1 ; ++i){
    for(int j = i+1 ; j < N ; ++j){
        if(S[j][i] >= 0.5){
            EFF1 = EP[j][i];
            EP[j][i]= (S[j][i]/(NUMBER_OF_PARENT_NUCLEI_DECAYS(Initial_Activity,time,half_life,T)*QI[j][i]/100.0)) - (D[j][i]/QN0C[j][i]);
            SEFF = SEFF + pow( ( (EFF1-EP[j][i]) / EP[j][i] ),2);
            NEFF = NEFF +1;
        }
            
        }
}
DEFF = TMath::Sqrt(SEFF/float(NEFF));
    cout << NIT << "   "<< DEFF << endl;
// INCREMENT ITERATION NUMBER AND LOOP BACK TO STAR
NIT = NIT + 1;
} while(NIT < 10 || (NIT < 20 && DEFF >= 1.0*pow(10,-6)));
   
    
double ERR;
ofstream F_OUT;
F_OUT.open(Form("/home/gula/Desktop/8MeV_analysis/F_OUT.txt"));
    cout << endl;
cout <<"============================================================================="<<endl;
cout <<"============================================================================="<<endl;
cout <<"GAMMA ENERGY"<<"       "<<"PHOTO-EFFICIENCY"<<"       "<<"PHOTO-EFFICIENC"<<endl;
cout <<"   (MeV)    "<<"       "<<"  (CORRECTED)   "<<"       "<<"(NOT CORRECTED)"<<endl;
cout <<"============================================================================="<<endl;
cout <<"============================================================================="<<endl;
for(int i = 0 ; i < N-1 ; ++i){
    for(int j = i ; j < N ; ++j){
     if (S[j][i] >= 0.5){
    ERR = TMath::Sqrt( pow( (QIE[j][i]/QI[j][i]) ,2) + pow( (DS[j][i]/S[j][i]) ,2) );
    //ERR = ( (  QIE[j][i] / QI[j][i] ) + ( DS[j][i] / S[j][i] ) );
    F_OUT<< EN[j]-EN[i] <<"      "<< EP[j][i] <<"      "<< EP[j][i]*ERR <<"      " << S[j][i]/(NUMBER_OF_PARENT_NUCLEI_DECAYS(Initial_Activity,time,half_life,T)*QI[j][i]/100.0)<<endl;
    cout <<"   "<<EN[j]-EN[i] <<"            "<< EP[j][i] <<"                "<<  S[j][i]/(NUMBER_OF_PARENT_NUCLEI_DECAYS(Initial_Activity,time,half_life,T)*QI[j][i]/100.0) <<endl;
          }
    }
}
F_OUT.close();

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    YOUR PHOTO-FRACTION GOES HERE. USED ONLY IF MODE=2.
//    REMEMBER THAT ALL ENERGIES ARE IN MeV.
//    NGEO = "GEOMETRY CODE"

double PHI(double EG, double NGEO){
    double E = 1000. * EG;
double PHI0;
        if (E > 0.0) {
            if(NGEO == 1){
                //                PHI0 = (0.0353744 + 0.60996 * TMath::Exp(-7.57933 * pow(10,-3) * E) + 0.33549 * TMath::Exp(-5.54383 * pow(10,-4) * E));
                PHI0 = (0.0307834 + 0.646614 * TMath::Exp(-7.7344 * pow(10,-3) * E) + 0.331807 * TMath::Exp(-5.17157 * pow(10,-4) * E));
                
            }
            else{
                PHI0 =  (0.0353744 + 0.60996 * TMath::Exp(-7.57933 * pow(10,-3) * E) +
                        0.33549 * TMath::Exp(-5.54383 * pow(10,-4) * E))  * ( 1.0554 + 8.14703 * pow(10,-8) * pow((E-1500.),2) );
            }
        }
        else
            PHI0 = 1.0;

return PHI0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   YOUR TOTAL EFFICIENCY GOES HERE. USED ONLY IF MODE=1.
//   REMEMBER THAT ALL ENERGIES ARE IN MeV.
//   NGEO = "GEOMETRY CODE"
//   APPROXIMATE TOTAL EFFICIENCY FOR 85% HPGe

double ETA(double EG, double NGEO){
    double ETA0;
    ETA0 = 0.085 + 0.05261 * pow((2.341/2.),2) / ( pow(EG,2) + pow((2.341/2.),2) );

    return ETA0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double NUMBER_OF_PARENT_NUCLEI_DECAYS(double Initial_Activity, double time, double half_life, double T){
     double final_Activity = Initial_Activity * TMath::Exp(- ( TMath::Log(2) / half_life ) * time);
     double N0 = final_Activity / ( TMath::Log(2) / half_life);
     double Nf = N0 * exp (- ( TMath::Log(2) / half_life ) * T);// 1080 : is the time during theexperiment
     double M = N0-Nf;
     return M;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    PHOTO-EFFICIENCY FUNCTION GOES HERE; USED FOR CALCULATING THE
//    EFFECTS OF TRANSITIONS WHOSE OBSERVED INTENSITIES ARE NOT GIVEN
//    IN THE INPUT FILE. NGEO = "GEOMETRY CODE"
//    NGEO = 1  --> Pb SHIELDED CLOSE GEOMETRY, 85% DETECTOR
//    NGEO = 2  --> Pb SHIELDED FAR GEOMETRY, 85% DETECTOR

double  EPS(double EG, double NGEO){
    double EFF;
    if (NGEO ==1){
        EFF = TMath::Exp( -3.4986 - 0.64347 * TMath::Log(EG) -2.59583 * pow(10,-2) * TMath::Log(pow(EG,2)));
    } else if(NGEO == 2) {
        EFF = TMath::Exp( -5.4827 -0.57871 * TMath::Log(EG) -3.45118 * pow(10,-2) * TMath::Log(pow(EG,2)));
    }else{
        EFF = 0.01/EG;
    }
    return EFF;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//PRINT OUT MATRIX A TO SCREEN; USED FOR DE-BUGGING

void Print(double Matrix[][N], int row, int column){
    int i, j;
    cout << "Output Matrix:" << endl;
    for(i = 0; i < row; ++i)
    {
        for(j = 0; j < column; ++j)
        {
            cout << Matrix[i][j] << " ";
            if(j == column - 1)
                cout << endl << endl;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//PRINT OUT vector A TO SCREEN; USED FOR DE-BUGGING

void display(double mult[N], int rowThird)
{
    cout << "Output Matrix:" << endl;
    for(int i = 0; i < rowThird; i++)
    {
        cout << mult[i] << endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//SET MATRIX A EQUAL TO THE IDENTITY MATRI

void Identity ( double Matrix[][N], int row, int column){
    
    for(int i = 0; i < row; ++i)
    {
        for(int j = 0; j < column; ++j)
        {
            if (i == j){
                Matrix[i][j] = 1;
            }
            else {
                Matrix[i][j] = 0;
            }
        }
    
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//ADD MATRIX A INTO MATRIX B

void Subtract(double firstMatrix[][N], double secondMatrix[][N], double thirdMatrix[][N], int row, int column){
    int i, j;
    
    for(i = 0; i < row; ++i){
        for(j = 0; j < column; ++j){
            thirdMatrix[i][j] = firstMatrix[i][j] - secondMatrix[i][j];
            }
        }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Multiplying matrix firstMatrix and secondMatrix and storing in array mult.

void MMmultiplyMatrices(double firstMatrix[][N], double secondMatrix[][N], double mult[][N], int rowFirst, int columnFirst, int rowSecond, int columnSecond)
{
    
    // Initializing elements of matrix mult to 0.
    for(int i = 0; i < rowFirst; ++i)
    {
        for(int j = 0; j < columnSecond; ++j)
        {
            mult[i][j] = 0;
        }
    }
    // Multiplying matrix firstMatrix and secondMatrix and storing in array mult.
    for(int i = 0; i < rowFirst; ++i)
    {
        for(int j = 0; j < columnSecond; ++j)
        {
            for(int k=0; k < columnFirst; ++k)
            {
                mult[i][j] += firstMatrix[i][k] * secondMatrix[k][j];
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Multiplying 1D array and 2D array and storing in array mult.

void MmultiplyMatrices(double firstMatrix[N], double secondMatrix[][N], double mult[N], int rowFirst, int rowSecond , int columnSecond)
{
    
    // Initializing elements of matrix mult to 0.
    for(int i = 0; i < rowSecond; i++)
    {
        mult[i] = 0;
    }
    
    // Multiplying matrix firstMatrix and secondMatrix and storing in array mult.
    for(int i = 0; i < rowFirst; i++)
    {
        for(int j = 0; j < columnSecond; j++)
        {
            mult[i] += firstMatrix[j] * secondMatrix[j][i];
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//ADD MATRIX A INTO MATRIX B

void Addition(double firstMatrix[][N], double secondMatrix[][N], double thirdMatrix[][N], int row, int column){
    int i, j;
    for(i = 0; i < row; ++i){
        for(j = 0; j < column; ++j){
            thirdMatrix[i][j] = firstMatrix[i][j] + secondMatrix[i][j];
            }
        }
}

 // Initializing elements of matrix A to 0.

void Zero( double Matrix[][N], int row, int column){
    int i, j;
    // Initializing elements of matrix A to 0.
    for(i = 0; i < row; ++i){
        for(j = 0; j < column; ++j){
            Matrix[i][j]=0;
            }
        }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//PUT MATRIX A INTO MATRIX B

void Put (double firstMatrix[][N], double secondMatrix[][N], int row, int column){
    int i, j;
//PUT MATRIX A INTO MATRIX B
    for(i = 0; i < row; ++i){
        for(j = 0; j < column; ++j){
        secondMatrix[i][j] = firstMatrix[i][j];
        }
    }
}










//0.3443      0.0103621      0.00257689      0.00992655
//0.4111      0.00975472      0.00383823      0.00853133
//0.7789      0.00509471      0.00114824      0.00498666


//0.3443      0.0103621      0.00257689      0.00992655
//0.4111      0.00975472      0.00383823      0.00853133
//0.7789      0.00509471      0.00114824      0.00498666

