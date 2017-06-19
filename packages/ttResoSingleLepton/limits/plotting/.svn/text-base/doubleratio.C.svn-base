//////////////////////////////////////////////////
// Plot the limits with error bands
// ebergeas@cern.ch Jan, Feb 2011
// Updated for 2011 data, May 2011
// Updated for 5 fb-1 data, June 2012
//////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <map>
#include "TString.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "functions.h"
#include "TLine.h"
#include "TArrow.h"
using namespace std;

void doubleratio() {

  bool doKKg=false;
  bool doZprime=!doKKg;
  //7TeV result
  double Zmass7[10]={500,600,700,800,1000,1300,1600,2000,2500,3000};
  double Zobs7[10]={5.05191, 7.10921, 4.63765, 1.61036, 0.433908, 0.117448, 0.0559945, 0.0382263, 0.0342212, 0.0312128};
  double Zexp7[10]={6.70909, 4.15,    2.33659, 1.44608, 0.491346, 0.148025, 0.0802695, 0.0422983, 0.0325333, 0.0276696};
  double Zm1s7[10]={3.67887, 2.41957, 1.46,   0.984545, 0.314286, 0.0904255,0.0485922, 0.0268544, 0.0216652, 0.0189536};
  double Zp1s7[10]={10.1538, 6.18088, 3.67292,1.8885,   0.736047, 0.212821, 0.115429,  0.0644112, 0.0477187, 0.0441};
  double Zm2s7[10]={1.8625,  1.22763, 0.705769,0.45875, 0.155357, 0.0529348,0.0303947, 0.0180769, 0.0152198, 0.014822};
  double Zp2s7[10]={14.8375, 9.1125,  5.20556, 3.03462, 1.095,    0.311071, 0.158553,  0.0859,    0.0661,    0.0721071};
  double Zxsec7[10]={19.58, 10.26,  5.62, 3.22, 1.16, 0.30,  0.09, 0.02,  0.003,  0.001};

  double Kmass7[10]={700,800,1000,1150,1300,1600,1800,2000,2250,2500};
  double Kobs7[10]={4.96192, 2.6416, 0.662837, 0.293317, 0.195234, 0.105768, 0.0716345, 0.0774572, 0.0751281, 0.0765797};
  double Kexp7[10]={3.49348, 1.86093,0.764027, 0.381304, 0.242973, 0.140235,  0.10493,   0.0891852, 0.083838,  0.0781481};
  double Km1s7[10] ={2.16223, 1.28925,0.513978, 0.242128, 0.150207, 0.0822302, 0.0660326, 0.0563386, 0.0541739, 0.0502564};
  double Km2s7[10]={1.02422, 0.622414, 0.210625, 0.121628, 0.0841622, 0.0476042, 0.0408929, 0.0363889, 0.0354592, 0.0334756};
  double Kp1s7[10] ={5.52105, 3.08679,1.14103,  0.584889, 0.367,    0.198115,  0.159432,  0.129412,  0.126328,  0.118692};
  double Kp2s7[10]={7.92917, 4.16818, 1.73088, 0.846, 0.523091, 0.275, 0.231389, 0.19375, 0.179861, 0.172292};
  double Kxsec7[10]={20.76, 11.59, 4.11, 2.07, 1.09, 0.35, 0.18, 0.10, 0.05, 0.03};
  //8TeV result
  TString path="../8TeVresults/combined/130407/";
  const double NZ=10;
  const double NKKg=13;
  double Zxsec8[NZ]={23.17,5.60,1.61, 0.57,0.21,0.0883521, 0.04, 0.0127335,0.007, 0.002};
  double Kxsec8[NKKg]={81.9270, 45.0410, 25.1950, 14.6300, 8.80670, 5.4726, 2.8164, 1.5233, 0.5004, 0.2556, 0.1369, 0.0670, 0.0351};
  double Zmass8[NZ]={500,750,1000,1250,1500,1750,2000,2250,2500,3000}; 
  double Kmass8[NKKg]={500,600,700,800,900,1000,1150,1300,1600,1800,2000,2250,2500};
  double Zobs8[NZ], Zexp8[NZ], Zm1s8[NZ], Zp1s8[NZ],Zm2s8[NZ],Zp2s8[NZ];
  double Kobs8[NKKg], Kexp8[NKKg], Km1s8[NKKg], Kp1s8[NKKg],Km2s8[NKKg],Kp2s8[NKKg];
  TString obs_filename = "Zprime-lim_pretag_syst_obs.dat";
  TString exp_filename = "Zprime-lim_pretag_syst_exp.dat";
  if (doKKg) {
    obs_filename="KKg-lim_pretag_syst_obs.dat";
    exp_filename="KKg-lim_pretag_syst_exp.dat";
  }

 ifstream str_obs(string(path+obs_filename).c_str());
  if(!str_obs.is_open()) {
    cout << "ERROR reading file " << path+obs_filename << endl;
  } 
 ifstream str_exp(string(path+exp_filename).c_str());
 if(!str_exp.is_open()) {
    cout << "ERROR reading file " << path+exp_filename << endl;
  }
 double ee, m95, m68, p68, p95; 
 int j=0, Mindex;

 if(doZprime) {
   cout << "Observed: ";
   j=0;
   while(!str_obs.eof() && j < NZ ) {
     str_obs >> Mindex >> ee;
     Zobs8[j]=ee;
     j++;
     cout << ee << " ";
   }
   cout << endl;
   cout << "Expected mean: "<<endl;
   j=0;
   while(!str_exp.eof() && j < NZ ) {
     str_exp >> Mindex >> ee >> m68 >> p68 >> m95 >> p95;
     Zexp8[j]=ee;
     Zm1s8[j]=m68; Zp1s8[j]=p68; 
     Zm2s8[j]=m95; Zp2s8[j]=p95;
     j++;
     cout << ee << " "<<m68<<" "<<p68<<" "<<m95<<" "<<p95<<endl;
   }
  cout << endl;
 }
 else if(doKKg) {
    cout << "Observed: ";
   j=0;
   while(!str_obs.eof() && j < NKKg ) {
     str_obs >> Mindex >> ee;
     Kobs8[j]=ee;
     j++;
     cout << ee << " ";
   }
   cout << endl;
   cout << "Expected mean: "<<endl;;
   j=0;
   while(!str_exp.eof() && j < NKKg ) {
     str_exp >> Mindex >> ee >> m68 >> p68 >> m95 >> p95;
     Kexp8[j]=ee;
     Km1s8[j]=m68; Kp1s8[j]=p68; 
     Km2s8[j]=m95; Kp2s8[j]=p95;
     j++;
     cout << ee << " "<<m68<<" "<<p68<<" "<<m95<<" "<<p95<<endl;
   }
 cout << endl;
 }
  // double Zobs8[10]={5.56354, 2.00278,  0.352268, 0.15641,   0.149204,  0.0852381, 0.0588308, 0.0465069, 0.0385295, 0.032064};
  // double Zexp8[10]={1.43099, 0.247619, 0.101618, 0.0647407, 0.0416842, 0.0295833, 0.0252632, 0.02144,   0.0163134, 0.00950162 };
 

  // double Kobs8[13]={9.23774, 5.26994, 3.46765, 2.21853, 1.23862, 0.706627, 0.297692, 0.251392, 0.185512, 0.112114, 0.096428, 0.0998155, 0.110219};
  // double Kexp8[13]={5.00727, 3.15, 1.88, 1.16696, 0.740625, 0.578495, 0.33082, 0.215316, 0.127755, 0.104052, 0.0845775, 0.0767647, 0.0799773};
  // double Km1s8[13]={1.72, 0.927632, 0.448101, 0.25875, 0.171983, 0.143226, 0.0972131, 0.0808647, 0.0505645, 0.0437129, 0.0369583, 0.0361786, 0.03824};
  // double Km2s8[13]={0.334034, 0.225375, 0.145982, 0.0702247, 0.0614162, 0.0539344, 0.0431452, 0.0387329, 0.0259091, 0.0230288, 0.0206686, 0.0211339, 0.0229865};
  // double Kp1s8[13]={10.4234, 6.19773, 3.93478, 2.25833, 1.65106, 1.27778, 0.705385, 0.4275, 0.253333, 0.2125, 0.173167, 0.149591, 0.165375}; 
  // double Kp2s8[13]={16.455, 9.76875, 6.35, 3.61154, 2.74688, 2.16222, 1.12143, 0.680156, 0.397308, 0.329583, 0.2705, 0.24525, 0.2765};
 
 
  map<double, double> Zobsmap7;
  map<double, double> Zexpmap7;
  map<double, double> Zm1smap7;
  map<double, double> Zm2smap7;
  map<double, double> Zp1smap7;
  map<double, double> Zp2smap7;
  map<double, double> Zobsmap8;
  map<double, double> Zexpmap8;
  map<double, double> Zm1smap8;
  map<double, double> Zm2smap8;
  map<double, double> Zp1smap8;
  map<double, double> Zp2smap8;

  map<double, double> Kobsmap7;
  map<double, double> Kexpmap7;
  map<double, double> Km1smap7;
  map<double, double> Km2smap7;
  map<double, double> Kp1smap7;
  map<double, double> Kp2smap7;
  map<double, double> Kobsmap8;
  map<double, double> Kexpmap8;
  map<double, double> Km1smap8;
  map<double, double> Km2smap8;
  map<double, double> Kp1smap8;
  map<double, double> Kp2smap8;
  for (int i=0; i<10;i++) {
    Zobsmap7[Zmass7[i]]=Zobs7[i]/Zxsec7[i];
    Zexpmap7[Zmass7[i]]=Zexp7[i]/Zxsec7[i];
    Zm1smap7[Zmass7[i]]=Zm1s7[i]/Zxsec7[i];
    Zp1smap7[Zmass7[i]]=Zp1s7[i]/Zxsec7[i];
    Zm2smap7[Zmass7[i]]=Zm2s7[i]/Zxsec7[i];
    Zp2smap7[Zmass7[i]]=Zp2s7[i]/Zxsec7[i];
  }
  if (doZprime) {
    for (int i=0; i<NZ;i++) {
      Zobsmap8[Zmass8[i]]=Zobs8[i]/Zxsec8[i];
      Zexpmap8[Zmass8[i]]=Zexp8[i]/Zxsec8[i];
      Zm1smap8[Zmass8[i]]=Zm1s8[i]/Zxsec8[i];
      Zp1smap8[Zmass8[i]]=Zp1s8[i]/Zxsec8[i];
      Zm2smap8[Zmass8[i]]=Zm2s8[i]/Zxsec8[i];
      Zp2smap8[Zmass8[i]]=Zp2s8[i]/Zxsec8[i];
    }
  }
  for (int i=0; i<10;i++) {
    Kobsmap7[Kmass7[i]]=Kobs7[i]/Kxsec7[i];
    Kexpmap7[Kmass7[i]]=Kexp7[i]/Kxsec7[i];
    Km1smap7[Kmass7[i]]=Km1s7[i]/Kxsec7[i];
    Kp1smap7[Kmass7[i]]=Kp1s7[i]/Kxsec7[i];
    Km2smap7[Kmass7[i]]=Km2s7[i]/Kxsec7[i];
    Kp2smap7[Kmass7[i]]=Kp2s7[i]/Kxsec7[i];
  }
  if (doKKg) {
  for (int i=0; i<NKKg;i++) {
    Kobsmap8[Kmass8[i]]=Kobs8[i]/Kxsec8[i];
    Kexpmap8[Kmass8[i]]=Kexp8[i]/Kxsec8[i];
    Km1smap8[Kmass8[i]]=Km1s8[i]/Kxsec8[i];
    Kp1smap8[Kmass8[i]]=Kp1s8[i]/Kxsec8[i];
    Km2smap8[Kmass8[i]]=Km2s8[i]/Kxsec8[i];
    Kp2smap8[Kmass8[i]]=Kp2s8[i]/Kxsec8[i];
  } 
  }
 //common point
  double Zmass[5]={500,1000,2000,2500,3000};
  double Kmass[10]={700,800,1000,1150,1300,1600,1800,2000,2250,2500};
 
  double Zobs[5];
  double Zexp[5];
  double Zm1s[5];
  double Zp1s[5];
  double Zm2s[5];
  double Zp2s[5];
  double Kobs[10];
  double Kexp[10];
  double Kxsec[10];
  double Km1s[10];
  double Kp1s[10];
  double Km2s[10];
  double Kp2s[10];
  double err1a, err1b;
  double err2a, err2b;
  double err1,err2;
  double errp1a, errp1b;
  double errm1a, errm1b;
  double errp2a, errp2b;
  double errm2a, errm2b;
 double errp2a, errp2b;
 double errm1,errp1,errm2,errp2;
  if (doZprime) {
  for (int i=0; i<5; i++) {
    if (Zobsmap7.find(Zmass[i])!=Zobsmap7.end() && Zobsmap8.find(Zmass[i])!=Zobsmap8.end()) Zobs[i]=Zobsmap8[Zmass[i]]/Zobsmap7[Zmass[i]];
    if (Zexpmap7.find(Zmass[i])!=Zexpmap7.end() && Zexpmap8.find(Zmass[i])!=Zexpmap8.end()) Zexp[i]=Zexpmap8[Zmass[i]]/Zexpmap7[Zmass[i]];
    if (Zm1smap7.find(Zmass[i])!=Zm1smap7.end() && Zm1smap8.find(Zmass[i])!=Zm1smap8.end()) Zm1s[i]=Zm1smap8[Zmass[i]]/Zm1smap7[Zmass[i]];
    if (Zp1smap7.find(Zmass[i])!=Zp1smap7.end() && Zp1smap8.find(Zmass[i])!=Zp1smap8.end()) Zp1s[i]=Zp1smap8[Zmass[i]]/Zp1smap7[Zmass[i]];
    if (Zm2smap7.find(Zmass[i])!=Zm2smap7.end() && Zm2smap8.find(Zmass[i])!=Zm2smap8.end()) Zm2s[i]=Zm2smap8[Zmass[i]]/Zm2smap7[Zmass[i]];
    if (Zp2smap7.find(Zmass[i])!=Zp2smap7.end() && Zp2smap8.find(Zmass[i])!=Zp2smap8.end()) Zp2s[i]=Zp2smap8[Zmass[i]]/Zp2smap7[Zmass[i]];
    err1a=max(fabs(Zexpmap8[Zmass[i]]-Zm1smap8[Zmass[i]]), fabs(Zexpmap8[Zmass[i]]-Zp1smap8[Zmass[i]]));
    err1b=max(fabs(Zexpmap7[Zmass[i]]-Zm1smap7[Zmass[i]]), fabs(Zexpmap7[Zmass[i]]-Zp1smap7[Zmass[i]]));
    err2a=max(fabs(Zexpmap8[Zmass[i]]-Zm2smap8[Zmass[i]]), fabs(Zexpmap8[Zmass[i]]-Zp2smap8[Zmass[i]]));
    err2b=max(fabs(Zexpmap7[Zmass[i]]-Zm2smap7[Zmass[i]]), fabs(Zexpmap7[Zmass[i]]-Zp2smap7[Zmass[i]]));
    err1=sqrt(pow(err1a/Zexpmap8[Zmass[i]],2)+pow(err1b/Zexpmap7[Zmass[i]],2))*Zexp[i];
    err2=sqrt(pow(err2a/Zexpmap8[Zmass[i]],2)+pow(err2b/Zexpmap7[Zmass[i]],2))*Zexp[i];
    Zm1s[i]= (Zexp[i]-err1)>0?(Zexp[i]-err1):0;
    Zm2s[i]= (Zexp[i]-err2)>0?(Zexp[i]-err2):0;
    Zp1s[i]= Zexp[i]+err1;
    Zp2s[i]= Zexp[i]+err2;
    cout<<Zmass[i]<<" "<<err1a<<" "<<err1b<<" "<<err1<<endl;
     cout<<Zmass[i]<<" "<<err2a<<" "<<err2b<<" "<<err2<<endl;
    // if (Zm1smap7.find(Zmass[i])!=Zm1smap7.end() && Zm1smap8.find(Zmass[i])!=Zm1smap8.end()) Zm1s[i]=sqrt(pow(Zm1smap8[Zmass[i]]/Zexpmap8[Zmass[i]],2)+pow(Zm1smap7[Zmass[i]]/Zexpmap7[Zmass[i]],2))*Zexp[i];
    // Zp1s[i]=2*Zexp[i]-Zm1s[i];
    // if (Zp1smap7.find(Zmass[i])!=Zp1smap7.end() && Zp1smap8.find(Zmass[i])!=Zp1smap8.end()) Zp1s[i]=sqrt(pow(Zp1smap8[Zmass[i]]/Zexpmap8[Zmass[i]],2)+pow(Zp1smap7[Zmass[i]]/Zexpmap7[Zmass[i]],2))*Zexp[i];
    //if (Zm2smap7.find(Zmass[i])!=Zm2smap7.end() && Zm2smap8.find(Zmass[i])!=Zm2smap8.end()) Zm2s[i]=sqrt(pow(Zm2smap8[Zmass[i]]/Zexpmap8[Zmass[i]],2)+pow(Zm2smap7[Zmass[i]]/Zexpmap7[Zmass[i]],2))*Zexp[i];
   //if (Zp2smap7.find(Zmass[i])!=Zp2smap7.end() && Zp2smap8.find(Zmass[i])!=Zp2smap8.end()) Zp2s[i]=sqrt(pow(Zp2smap8[Zmass[i]]/Zexpmap8[Zmass[i]],2)+pow(Zp2smap7[Zmass[i]]/Zexpmap7[Zmass[i]],2))*Zexp[i];
    //Zp2s[i]=2*Zexp[i]-Zm2s[i];
    
    Zmass[i]=Zmass[i]/1000;
    cout<<Zmass[i]<<" "<<Zobs[i]<<" "<<Zexp[i]<<" "<<Zm1s[i]<<" "<<Zp1s[i]<<" "<<Zm2s[i]<<" "<<Zp2s[i]<<endl;
  }
  }
  if (doKKg) {
  for (int i=0; i<10; i++) {
    if (Kobsmap7.find(Kmass[i])!=Kobsmap7.end() && Kobsmap8.find(Kmass[i])!=Kobsmap8.end()) Kobs[i]=Kobsmap8[Kmass[i]]/Kobsmap7[Kmass[i]];
    if (Kexpmap7.find(Kmass[i])!=Kexpmap7.end() && Kexpmap8.find(Kmass[i])!=Kexpmap8.end()) Kexp[i]=Kexpmap8[Kmass[i]]/Kexpmap7[Kmass[i]];
    // err1a=max(fabs(Kexpmap8[Kmass[i]]-Km1smap8[Kmass[i]]), fabs(Kexpmap8[Kmass[i]]-Kp1smap8[Kmass[i]]));
    // err1b=max(fabs(Kexpmap7[Kmass[i]]-Km1smap7[Kmass[i]]), fabs(Kexpmap7[Kmass[i]]-Kp1smap7[Kmass[i]]));
    // err2a=max(fabs(Kexpmap8[Kmass[i]]-Km2smap8[Kmass[i]]), fabs(Kexpmap8[Kmass[i]]-Kp2smap8[Kmass[i]]));
    // err2b=max(fabs(Kexpmap7[Kmass[i]]-Km2smap7[Kmass[i]]), fabs(Kexpmap7[Kmass[i]]-Kp2smap7[Kmass[i]]));
    // err1a=min(fabs(Kexpmap8[Kmass[i]]-Km1smap8[Kmass[i]]), fabs(Kexpmap8[Kmass[i]]-Kp1smap8[Kmass[i]]));
    // err1b=min(fabs(Kexpmap7[Kmass[i]]-Km1smap7[Kmass[i]]), fabs(Kexpmap7[Kmass[i]]-Kp1smap7[Kmass[i]]));
    // err2a=min(fabs(Kexpmap8[Kmass[i]]-Km2smap8[Kmass[i]]), fabs(Kexpmap8[Kmass[i]]-Kp2smap8[Kmass[i]]));
    // err2b=min(fabs(Kexpmap7[Kmass[i]]-Km2smap7[Kmass[i]]), fabs(Kexpmap7[Kmass[i]]-Kp2smap7[Kmass[i]]));
    // err1=sqrt(pow(err1a/Kexpmap8[Kmass[i]],2)+pow(err1b/Kexpmap7[Kmass[i]],2))*Kexp[i];
    // err2=sqrt(pow(err2a/Kexpmap8[Kmass[i]],2)+pow(err2b/Kexpmap7[Kmass[i]],2))*Kexp[i];
    errm1a=fabs(Kexpmap8[Kmass[i]]-Km1smap8[Kmass[i]]);
    errp1a=fabs(Kexpmap8[Kmass[i]]-Kp1smap8[Kmass[i]]);
    errm2a=fabs(Kexpmap8[Kmass[i]]-Km2smap8[Kmass[i]]);
    errp2a=fabs(Kexpmap8[Kmass[i]]-Kp2smap8[Kmass[i]]);
    errm1b=fabs(Kexpmap7[Kmass[i]]-Km1smap7[Kmass[i]]);
    errp1b=fabs(Kexpmap7[Kmass[i]]-Kp1smap7[Kmass[i]]);
    errm2b=fabs(Kexpmap7[Kmass[i]]-Km2smap7[Kmass[i]]);
    errp2b=fabs(Kexpmap7[Kmass[i]]-Kp2smap7[Kmass[i]]);
    errm1=sqrt(pow(errm1a/Kexpmap8[Kmass[i]],2)+pow(errm1b/Kexpmap7[Kmass[i]],2))*Kexp[i];
    errm2=sqrt(pow(errm2a/Kexpmap8[Kmass[i]],2)+pow(errm2b/Kexpmap7[Kmass[i]],2))*Kexp[i];
    errp1=sqrt(pow(errp1a/Kexpmap8[Kmass[i]],2)+pow(errp1b/Kexpmap7[Kmass[i]],2))*Kexp[i];
    errp2=sqrt(pow(errp2a/Kexpmap8[Kmass[i]],2)+pow(errp2b/Kexpmap7[Kmass[i]],2))*Kexp[i];
    Km1s[i]= (Kexp[i]-errm1)>0?(Kexp[i]-errm1):0;
    Km2s[i]= (Kexp[i]-errm2)>0?(Kexp[i]-errm2):0;
    Kp1s[i]= Kexp[i]+errp1;
    Kp2s[i]= Kexp[i]+errp2;

    // if (Km1smap7.find(Kmass[i])!=Km1smap7.end() && Km1smap8.find(Kmass[i])!=Km1smap8.end()) Km1s[i]=Km1smap8[Kmass[i]]/Km1smap7[Kmass[i]];
    // if (Kp1smap7.find(Kmass[i])!=Kp1smap7.end() && Kp1smap8.find(Kmass[i])!=Kp1smap8.end()) Kp1s[i]=Kp1smap8[Kmass[i]]/Kp1smap7[Kmass[i]];
    // if (Km2smap7.find(Kmass[i])!=Km2smap7.end() && Km2smap8.find(Kmass[i])!=Km2smap8.end()) Km2s[i]=Km2smap8[Kmass[i]]/Km2smap7[Kmass[i]];
    // if (Kp2smap7.find(Kmass[i])!=Kp2smap7.end() && Kp2smap8.find(Kmass[i])!=Kp2smap8.end()) Kp2s[i]=Kp2smap8[Kmass[i]]/Kp2smap7[Kmass[i]];
    // Km1s[i]=sqrt(pow(Km1smap8[Kmass[i]]/Kexpmap8[Kmass[i]],2)+pow(Km1smap7[Kmass[i]]/Kexpmap7[Kmass[i]],2))*Kexp[i];
    // Km2s[i]=sqrt(pow(Km2smap8[Kmass[i]]/Kexpmap8[Kmass[i]],2)+pow(Km2smap7[Kmass[i]]/Kexpmap7[Kmass[i]],2))*Kexp[i];
    // Kp1s[i]=2*Kexp[i]-Km1s[i];
    //Kp2s[i]=2*Kexp[i]-Km2s[i];
 
    cout<<Kmass[i]<<" "<<Kobs[i]<<" "<<Kexp[i]<<" "<<Kexpmap8[Kmass[i]]<<" "<<Kexpmap7[Kmass[i]]<<endl;
   Kmass[i]=Kmass[i]/1000;
  }
  }
  
  bool plotinTeV=true;
  bool doLogy=true;
  //Wuppertal style
  gStyle->SetOptStat(0);

  gStyle->SetCanvasBorderMode(0); //frame color of canvas
  gStyle->SetCanvasColor(0);  //bkrd color of canvas
  gStyle->SetStatBorderSize(0); //frame style of stat-box 1

  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);

  gStyle->SetLineWidth(3); // width of ticks
  gStyle->SetPadTickX(1); //0//1:ticks on upper,2: ticks+labels on upper xaxis
  gStyle->SetPadTickY(1); //0

  gStyle->SetPadLeftMargin(0.16); // 0.21// 0.18
  gStyle->SetPadRightMargin(0.05); //0.08
  gStyle->SetPadTopMargin(0.05); //0.07
  gStyle->SetPadBottomMargin(0.17);
  gStyle->SetFrameLineWidth(3);
  gStyle->SetFrameFillColor(0);
  gStyle->SetPaintTextFormat(".2f");
  
  //Added by ebergeas
  gStyle->SetTitleOffset(1.3,"Y");
  gStyle->SetTitleOffset(1.3,"X");
  gStyle->SetTitleOffset(1.2,"Z");
  Int_t font=42;
  Double_t tsize=0.05; 
  Double_t lsize=0.04; 
  gStyle->SetTextFont(font);
  gStyle->SetTextSize(tsize);
  gStyle->SetLabelFont(font,"x");
  gStyle->SetTitleFont(font,"x");
  gStyle->SetLabelFont(font,"y");
  gStyle->SetTitleFont(font,"y");
  gStyle->SetLabelFont(font,"z");
  gStyle->SetTitleFont(font,"z");

  gStyle->SetLabelSize(lsize,"x");
  gStyle->SetTitleSize(tsize,"x");
  gStyle->SetLabelSize(lsize,"y");
  gStyle->SetTitleSize(tsize,"y");
  gStyle->SetLabelSize(lsize,"z");
  gStyle->SetTitleSize(tsize,"z");


  const int lwidth=3;

  /////////////////////////////////////////////////
  /// Draw plots
  //Colours
  int col_theoryband=kRed+1, col_theoryline=kRed+2;
  //if(doDrawTheoryUncert) col_theoryline=kRed+4;
  int col_68=kGreen, col_95=kYellow;

  //Style
  int sty_theoryline=7; //1-solid, 7-dashed
  int fill_theory=3001;
  int sty_exp=2;

  TString lt_exp="Exp. 95% CL upper limit";
  TString lt_obs="Obs. 95% CL upper limit";
  TString lt_68="Exp. 1 #sigma uncertainty";
  TString lt_95="Exp. 2 #sigma uncertainty";
  //TString lt_th="Leptophobic Z' (LO x 1.3)";
  //TString lt_th_BH="Quantum Black Hole";
  //TString lt_th_KKg="Kaluza-Klein gluon (LO)";
  //TString lt_th_KKGrav="Kaluza-Klein graviton";

  //axis limits
  double x_ax_xs[2]={510., 2000.};
  if(doKKg) { 
    x_ax_xs[0]=500; //795; //805; //715
    x_ax_xs[1]=2500; //1885; //1893;
    //x_ax_xs[1]=4000;
  }
  if(doZprime) { 
    //x_ax_xs[0]=700; //715; //710;  
    // x_ax_xs[0]=618; //617
    // x_ax_xs[1]=1882; //1885; 
    x_ax_xs[0]=500; //617
    x_ax_xs[1]=3000; //1885; 
    //  x_ax_xs[1]=2100; //1885; 
 
  } 
 
  if(plotinTeV) {
    x_ax_xs[0]=x_ax_xs[0]/1000.;
    x_ax_xs[1]=x_ax_xs[1]/1000.;
  }

  double y_ax_xs[2];
   y_ax_xs[0]=0; y_ax_xs[1]=4.;
TLine* l1;
 l1=new TLine(x_ax_xs[0],1,x_ax_xs[1],1);
   l1->SetLineColor(4);
   l1->SetLineStyle(2);
   
  //if(doStat) {
  // if(doLogy) {
  //   if(doZprime)    { y_ax_xs[0]=5e-3;  y_ax_xs[1]=5e3; } //9e-3
  //   // if(doZprime)    { y_ax_xs[0]=1e-1;  y_ax_xs[1]=1e3; }
    
  //   //  if(doZprime)    { y_ax_xs[0]=5e-3;  y_ax_xs[1]=5e3; }
   
  //   if(doKKg)       { y_ax_xs[0]=1e-2;  y_ax_xs[1]=6e3; } //2e-2
   
  // }
  // else {
  //   if(doZprime)    { y_ax_xs[0]=7e-1;  y_ax_xs[1]=1e3; }
    
  //   if(doKKg)       { y_ax_xs[0]=2e-2;  y_ax_xs[1]=6e2; }
   
  // }



  //axis titles
  TString x_title="Z' mass ";
  TString y_title="#mu_{95} @8 TeV/ #mu_{95} @7 TeV";

  
  if(doKKg) {
    x_title="g_{KK} mass ";
  }

  if(plotinTeV) x_title+="[TeV]";
  else x_title+="[GeV]";

  //Draw
  int can_w=2400, can_h=2400;
  bool doRectangular=true;
  if(doRectangular) {
    //can_h=2100; can_w= 2970;
    can_h=1200; can_w=1600;
  }

  TCanvas* cv_xs = new TCanvas("cv_xs","Cross section",can_w,can_h);
  //   TCanvas* cv_xs = new TCanvas("cv_xs","Cross section",50,50,2400,2400);


  TGraph *gr_ax_xs = new TGraph(2, x_ax_xs, y_ax_xs);
  //   gr_ax_xs->SetMarkerStyle(1);
  gr_ax_xs->SetMarkerColor(kWhite);
  gr_ax_xs->SetTitle("");
  gr_ax_xs->GetXaxis()->SetTitle(x_title);
    gr_ax_xs->GetYaxis()->SetTitle(y_title);
  TGraph *gr_expected_xs = new TGraph(5, Zmass, Zexp);
  if (doKKg) gr_expected_xs = new TGraph(10, Kmass, Kexp);
  //   gr_expected_xs->SetMarkerStyle(1);
  gr_expected_xs->SetLineColor(kBlack);
  gr_expected_xs->SetLineStyle(sty_exp); //2
  gr_expected_xs->SetLineWidth(lwidth);

  TGraph *gr_observed_xs ;
  TGraph *gr_errorBand68_xs; 
  TGraph *gr_errorBand95_xs; 
  if (doKKg) { gr_observed_xs = new TGraph(10, Kmass, Kobs);
       gr_errorBand68_xs=getErrorBand(10, Kmass, Kp1s, Km1s);
       gr_errorBand95_xs=getErrorBand(10, Kmass, Kp2s, Km2s);
       gr_errorBand68_xs->SetFillColor(kGreen);  
       gr_errorBand68_xs->SetLineColor(kGreen);
       gr_errorBand95_xs->SetFillColor(kYellow);    
       gr_errorBand95_xs->SetLineColor(kYellow);
  }
  else {
       gr_observed_xs = new TGraph(5, Zmass, Zobs);
       gr_errorBand68_xs=getErrorBand(5, Zmass, Zp1s, Zm1s);
       gr_errorBand95_xs=getErrorBand(5, Zmass, Zp2s, Zm2s);
       gr_errorBand68_xs->SetFillColor(kGreen);  
       gr_errorBand68_xs->SetLineColor(kGreen);
       gr_errorBand95_xs->SetFillColor(kYellow);    
       gr_errorBand95_xs->SetLineColor(kYellow);
  }
  gr_observed_xs->SetMarkerStyle(20);
  gr_observed_xs->SetLineColor(kBlack);
  gr_observed_xs->SetLineWidth(lwidth);

  // TGraph *gr_errorBand68_xs=getErrorBand(Nmasses, signal_mass, limHigh68, limLow68);
  // TGraph *gr_errorBand95_xs=getErrorBand(Nmasses, signal_mass, limHigh95, limLow95);
  // TGraph *gr_expected_xs_r;
  // TGraph *gr_observed_xs_r;
  // TGraph *gr_errorBand68_xs_r;
  // TGraph *gr_errorBand95_xs_r;
  // TLine* l1;
  gr_ax_xs->Draw("AP");
  gr_ax_xs->GetXaxis()->SetRangeUser(x_ax_xs[0],x_ax_xs[1]);
 gr_ax_xs->GetYaxis()->SetRangeUser(y_ax_xs[0],y_ax_xs[1]);
  cv_xs->Update();
 
  // gr_errorBand95_xs->Draw("FL");
  // gr_errorBand68_xs->Draw("FL");
    gr_expected_xs->Draw("PL");
    gr_observed_xs->Draw("PL");
 l1->Draw();
   double lxlow=0.45, lylow= 0.66, lxup=0.92, lyup=0.98;
  TLegend *leg_limvsmass_xs = new TLegend(lxlow, lylow, lxup, lyup);
 
    leg_limvsmass_xs->AddEntry(gr_observed_xs, lt_obs,"pl");
    leg_limvsmass_xs->Draw();

  leg_limvsmass_xs->AddEntry(gr_expected_xs, lt_exp,"pl");      
  // leg_limvsmass_xs->AddEntry(gr_errorBand68_xs, lt_68,"lf");
  // leg_limvsmass_xs->AddEntry(gr_errorBand95_xs, lt_95,"lf");
  leg_limvsmass_xs->SetHeader("");
  leg_limvsmass_xs->SetFillColor(0);
  leg_limvsmass_xs->SetBorderSize(0);
  leg_limvsmass_xs->Draw();

  //ATLAS labels
  double a_x=lxlow+0.02, a_y=lylow-0.05; //, l_y=lyup-0.11, l_x=0.19;
  //if(doAtlasLabel) ATLASLabel(a_x, a_y, cv_xs,true);
  double l_y=lyup-0.17, l_x=0.19;
  //ATLASCMEIntL(l_x,l_y,LUMI,cv_xs);

  gPad->RedrawAxis();
  // gPad->SetLogy(1);
 

  // TString plotname="limits"+signal+tool+stsy+"_obs_exp"+sel+linlog;
  // if (!doratio) {
  // cv_xs->Print(plotname+".pdf","pdf");
  // 
  // cv_xs->Print(plotname+".png","png");
  // }
  // else  {
  // cv_xs->Print(plotname+"_ratio.pdf","pdf");
  // cv_xs->Print(plotname+"_ratio.eps","eps");
  // cv_xs->Print(plotname+"_ratio.png","png");
  // }

  cv_xs->Print("double_ratio.eps");
  cv_xs->Print("double_ratio.pdf");
}
