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

void plotLimits_r17(TString path="../8TeVresults/combined/130129",TString selection="_1TeV", TString channel="multi", const int Nmasses=5, bool doObserved=true, bool doAtlasLabel=true, bool doLogy=true, bool doRectangular=false);

void plotLimits_r17(TString path, TString selection, TString channel, const int Nmasses, bool doObserved, bool doAtlasLabel, bool doLogy, bool doRectangular) {
  
  bool debug=true;

  if( debug ) std::cout << "In function plotLimits_r17(...) " << std::endl;

  const int doPrintXsec=true;

  bool doZp13factor=true;
  
  bool plotinTeV=true;
  
  const bool doDrawTheoryUncert=false;

  double LUMI = 20.3;


  //bool printTex = false; 
  bool printTex = true; //print limits as a TeX table
  bool print2sigma=false; //prints only 1 sigma limits, not 2 sigma
  //bool print2sigma=true; //prints 1 and 2 sigma limits

  bool do95Lim=true; //Plot 1 and 2 sigma limits

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

  cout<<"selection: " << selection<<endl;
  const int lwidth=3;

  assert(channel == "ele" || channel == "muo" || channel == "multi");

  //Automatic settings -- what kind of limits are we plotting?
  bool doBlackHole=false;
  bool doKKg=false;
  bool doZprime=false;
  bool doHH=false;
  bool doZprime3p=true;  //add 3% theory line
  bool doKKGrav=false;
  bool doratio=false;
  
  bool doKKgWidth=false;

  bool doCMS=false;
  
  if(selection.Contains("BH-lim")  || selection.Contains("qbh")) {
    doBlackHole=true; 
  }
  else if(selection.Contains("KKg-lim") || selection.Contains("kkg")) {
    doKKg=true; 
  }
  else if(selection.Contains("KKg-Width")) {
    doKKgWidth=true; 
  }
  else if(selection.Contains("KKGrav-lim")) {
    doKKGrav=true; 
    cout<<"KKGrav test \n";
  }
  else if(selection.Contains("HH-lim")) {
    doHH=true; 
    cout<<"HH test \n";
  }
  else doZprime=true;

  if(doKKgWidth) plotinTeV=false;

  bool doStat=true; //Stat only limits
  if(selection.Contains("_sys") || path.Contains("sys")) doStat=false;



  //Legend settings
  //double lxlow=0.45, lylow= 0.6, lxup=0.92, lyup=0.92;
  double lxlow=0.45, lylow= 0.66, lxup=0.92, lyup=0.98;
  if(doStat) {
    lylow= 0.6; lyup=0.92;
  }

  //File names for limits
  TString obs_filename = "obs_limit_events_"+channel+"_"+selection+".dat";
  TString exp_filename = "exp_limit_events_68_95_"+channel+"_"+selection+".dat";
  TString e2x_filename = "events2pb_"+channel+"_"+selection+".dat";

  if(selection.Contains("-lim")) { //if limits are given in pb,
                                   //not as number of sign events
    obs_filename = selection+"_obs.dat";
    exp_filename = selection+"_exp.dat";
    e2x_filename = "null";
  }

  //Open .dat files with the limits
  if(!path.EndsWith("/")) path+="/";
  cout << "Files " << path << exp_filename << endl;
  cout << "      " << path << e2x_filename << endl;
  cout << "      " << path << obs_filename << endl;

  ifstream str_obs(string(path+obs_filename).c_str());
  if(!str_obs.is_open()) {
    cout << "ERROR reading file " << path+obs_filename << endl;
  }
  else std::cout << "Successfully openend file " << path+obs_filename << std::endl;

  ifstream str_exp(string(path+exp_filename).c_str());
  if(!str_exp.is_open()) {
    cout << "ERROR reading file " << path+exp_filename << endl;
  }
  else std::cout << "Successfully openend file " << path+exp_filename << std::endl;

  ifstream str_e2x(string(path+e2x_filename).c_str());
  if(!str_e2x.is_open() && e2x_filename!="null") {
    cout << "ERROR reading file " << path+e2x_filename << endl;
  }
  else std::cout << "Successfully openend file " << path+e2x_filename << std::endl;
  
  
  //===========================================
  //Extract limits from files
  //===========================================
  double EventsObs[Nmasses], EventsExp[Nmasses], Events2pb[Nmasses], ee;
  double EventsExp_m95[Nmasses], EventsExp_m68[Nmasses];
  double EventsExp_p68[Nmasses], EventsExp_p95[Nmasses]; 
  double signal_mass[Nmasses];
  double m95, m68, p68, p95;
  int j=0, Mindex;
  
  //Observed limits
  cout << "Observed: ";
  j=0;
  while(!str_obs.eof() && j < Nmasses ) {
    if(selection.Contains("-lim")) {
      str_obs >> Mindex >> ee;
      cout << "Printing read out " << endl;
      cout << "str_exp " << str_exp << endl;
      cout << Mindex << " " << ee << endl;
    }
    else {
      str_obs >> ee;
    }
    EventsObs[j]=ee;
    j++;
    cout << ee << " ";
  }
  cout << endl;

  if(e2x_filename!="null") { //if limits are given as events
    cout << "Events2pb: ";
    j=0;
    while(!str_e2x.eof() && j < Nmasses ) {
      str_e2x >> Mindex >> ee;
      Events2pb[j]=ee;
      j++;
      cout << ee << " ";
    }
    cout << endl;
  }
  else {
    for(int jj=0;jj<Nmasses;++jj) {
      Events2pb[jj]=1.0;
    }
  }

  //Expected limits
  cout << "Expected mean: ";
  j=0;
  if(selection.Contains("-lim_")) { //order: M, m, -1s, +1s, -2s, +2s
    while(!str_exp.eof() && j < Nmasses ) {
      str_exp >> Mindex >> ee >> m68 >> p68 >> m95 >> p95;
      cout << "Printing read out " << endl;
      cout << "str_exp " << str_exp << endl;
      cout<<Mindex<<" "<<ee<<" "<< m68 <<" "<< p68 <<" "<< m95 <<" "<< p95<<endl;
      EventsExp[j]=ee;
      EventsExp_m68[j]=m68; EventsExp_p68[j]=p68; 
      EventsExp_m95[j]=m95; EventsExp_p95[j]=p95;
      if(plotinTeV) signal_mass[j]=Mindex/1000.;
      else signal_mass[j]=Mindex;
      j++;
      //cout << ee << " ";
      // cout << p95 << ", ";
    }
  }
  else { //order: i, -2s, -1s, m, +1s, +2s
    while(!str_exp.eof() && j < Nmasses ) {
      str_exp >> Mindex >> m95 >> m68 >> ee >> p68 >> p95;
      EventsExp[j]=ee;
      EventsExp_m95[j]=m95; EventsExp_m68[j]=m68; 
      EventsExp_p68[j]=p68; EventsExp_p95[j]=p95;
      j++;
      cout << ee << " ";
    }
    getMass(selection, signal_mass, Nmasses, plotinTeV);
  }
  cout << endl;

  //Scale events (only relevant if Events2pb !=1)
  double limitMean[Nmasses], limitObserved[Nmasses];
  double limHigh68[Nmasses], limLow68[Nmasses], limHigh95[Nmasses], limLow95[Nmasses];
 
  for(int i=0;i<Nmasses;++i) {
    // if (i<Nmasses-2)
    limitMean[i]= EventsExp[i]*Events2pb[i];
    limHigh68[i]= EventsExp_p68[i]*Events2pb[i];
    limLow68[i] = EventsExp_m68[i]*Events2pb[i];
    limHigh95[i]= EventsExp_p95[i]*Events2pb[i];
    limLow95[i] = EventsExp_m95[i]*Events2pb[i];
    limitObserved[i]= EventsObs[i]*Events2pb[i];
  }


  //===========================================
  //Draw final plots
  //===========================================
  
  //------------------------------------
  //Colours
  int col_theoryband=kRed+1, col_theoryline=kRed+2;
  if(doDrawTheoryUncert) col_theoryline=kRed+4;
  int col_68=kGreen, col_95=kYellow;

  //------------------------------------
  //Style
  int sty_theoryline=7; //1-solid, 7-dashed
  int fill_theory=3001;
  int sty_exp=2;

  //------------------------------------
  // Make graphs
  
  //------------------------------------
  //Z' theory
  TString theorypath="theory/";
  // const int Ntheory_DW= (doZp13factor? 216:189);
    const int Ntheory_DW= 14;
  TString theoryfile="zprimexsecLHC_cteq6l1_172.5_100M_20110210+20110505.dat";
   if(doZp13factor) theoryfile="zprimexsecLHC2011_cteq6l1_172.5_20120223.dat";
   bool do8TeV = true;
   if (do8TeV) LUMI=20.3;
   if(doZp13factor && do8TeV) { theoryfile="zprimexsecLHC8TeVwidth_0_012.dat"; }
  double Mt_DW[Ntheory_DW];
  double xsect_DW[Ntheory_DW];
  double xsect_errup_DW[Ntheory_DW];
  double xsect_errdown_DW[Ntheory_DW];
  double e0[Ntheory_DW];
  double Mm, xc; //, xu, xd;
  ifstream str_theory(string(theorypath+theoryfile).c_str());
  if(!str_theory.is_open()) {
    cout << "ERROR reading file " << theoryfile << endl;
  } 
  j=0;
  while(!str_theory.eof() && j < Ntheory_DW ) {
    str_theory >> Mm >> xc;
    cout<<Mm<<" "<<xc<<endl;
    if(plotinTeV) Mt_DW[j]=Mm/1000.;
    else Mt_DW[j]=Mm;  
    //Mt_DW[j]=Mm;
    xsect_DW[j]=xc;
    //Errors +8.3% and -7.4%. Gustaaf Brooijmans on ttbarres list 9 Feb
    xsect_errup_DW[j]=xc*0.083;
    xsect_errdown_DW[j]=xc*0.074;
    if(doZp13factor) {
      xsect_DW[j]=1.3*xc;
      xsect_errup_DW[j]=1.3*xc*0.083;
      xsect_errdown_DW[j]=1.3*xc*0.074;
    }
    cout<<xsect_DW[j]<<endl;
    e0[j]=0;
    j++;
  }

  //error bands  
  TGraph *gr_theory_DW=getErrorBand(Ntheory_DW, Mt_DW, xsect_DW, xsect_errdown_DW, xsect_errup_DW);
  gr_theory_DW->SetLineStyle(sty_theoryline); //1
  gr_theory_DW->SetLineWidth(lwidth);
  gr_theory_DW->SetLineColor(col_theoryline); //kRed+2
  gr_theory_DW->SetFillColor(col_theoryline); //kRed+1
  gr_theory_DW->SetFillStyle(fill_theory);  //3001
  
  //Central value
  TGraph *gr_theory_DWc=new TGraph(Ntheory_DW, Mt_DW, xsect_DW);
  // gr_theory_DWc->SetLineStyle(1);
  // gr_theory_DWc->SetLineWidth(lwidth);
  // gr_theory_DWc->SetLineColor(kRed+4);
  gr_theory_DWc->SetLineStyle(sty_theoryline);
  gr_theory_DWc->SetMarkerStyle(1);
  gr_theory_DWc->SetLineWidth(lwidth);
  gr_theory_DWc->SetLineColor(col_theoryline);
  gr_theory_DWc->SetMarkerColor(col_theoryline);

  //for ratio
  //limits over theory
  double limitMean_r[Nmasses], limitObserved_r[Nmasses];
  double limHigh68_r[Nmasses], limLow68_r[Nmasses], limHigh95_r[Nmasses], limLow95_r[Nmasses];
  map<double, double> xsecmap;
  if(doZprime) {
  
      for(int i=0;i<Ntheory_DW; i++) xsecmap[Mt_DW[i]]=xsect_DW[i];
      if (do8TeV) {
      xsecmap[1.75]=gr_theory_DWc->Eval(1.75,0,"S");
      if (doratio) cout<<"Zprime 1.75 xsec "<< gr_theory_DWc->Eval(1.75,0,"S")<<endl;
      xsecmap[2.25]=gr_theory_DWc->Eval(2.25,0,"S");
      if (doratio) cout<<"Zprime 2.25 xsec "<< gr_theory_DWc->Eval(2.25,0,"S")<<endl;
      }
      else {
	xsecmap[3.]=0.002;
	cout<<"Zprime 3 xsec "<< 0.002<<endl;
      }
  }
  ///zprime 3% width
  TString theorypath_Z3p="theory/";
  TString theoryfile_Z3p="zprimexsecLHC8TeVwidth_0_03.dat";

  double Mt_Z3p[Ntheory_DW];
  double xsect_Z3p[Ntheory_DW];
  double e0_Z3p[Ntheory_DW];
  ifstream str_theory_Z3p(string(theorypath_Z3p+theoryfile_Z3p).c_str());
  if(!str_theory_Z3p.is_open()) {
    cout << "ERROR reading file " << theoryfile_Z3p << endl;
  }
  j=0;

  while(!str_theory_Z3p.eof() && j < Ntheory_DW ) {
    str_theory_Z3p >> Mm >> xc;
    if(plotinTeV) Mt_Z3p[j]=Mm/1000.;
    else Mt_Z3p[j]=Mm;
    //Mt_BH[j]=Mm;
    xsect_Z3p[j]=xc;  e0_Z3p[j]=0;
    
    if(doZp13factor) {
      xsect_Z3p[j]=1.3*xc;
    }


    j++;
  }

  TGraph *gr_theory_Z3p=new TGraph(Ntheory_DW, Mt_Z3p, xsect_Z3p);
  gr_theory_Z3p->SetLineStyle(sty_theoryline+1);
  gr_theory_Z3p->SetMarkerStyle(2);
  gr_theory_Z3p->SetLineWidth(lwidth);
  gr_theory_Z3p->SetLineColor(4);
  gr_theory_Z3p->SetMarkerColor(4);



  ///zprime 2% width
  TString theorypath_Z2p="theory/";
  TString theoryfile_Z2p="zprimexsecLHC8TeVwidth_0_02.dat";

  double Mt_Z2p[Ntheory_DW];
  double xsect_Z2p[Ntheory_DW];
  double e0_Z2p[Ntheory_DW];
  ifstream str_theory_Z2p(string(theorypath_Z2p+theoryfile_Z2p).c_str());
  if(!str_theory_Z2p.is_open()) {
    cout << "ERROR reading file " << theoryfile_Z2p << endl;
  }
  j=0;

  while(!str_theory_Z2p.eof() && j < Ntheory_DW ) {
    str_theory_Z2p >> Mm >> xc;
    if(plotinTeV) Mt_Z2p[j]=Mm/1000.;
    else Mt_Z2p[j]=Mm;
    //Mt_BH[j]=Mm;
    xsect_Z2p[j]=xc;  e0_Z2p[j]=0;
    
    if(doZp13factor) {
      xsect_Z2p[j]=1.3*xc;
    }


    j++;
  }

  TGraph *gr_theory_Z2p=new TGraph(Ntheory_DW, Mt_Z2p, xsect_Z2p);
  gr_theory_Z2p->SetLineStyle(sty_theoryline+2);
  gr_theory_Z2p->SetMarkerStyle(3);
  gr_theory_Z2p->SetLineWidth(lwidth);
  gr_theory_Z2p->SetLineColor(28);
  gr_theory_Z2p->SetMarkerColor(28);




///CMS obs
  TString theorypath_CMSobs="theory/";
  TString theoryfile_CMSobs="zprimeCMSobserved.dat";
  int nCMS=26;

  double Mt_CMSobs[nCMS];
  double xsect_CMSobs[nCMS];
  double e0_CMSobs[nCMS];
  ifstream str_theory_CMSobs(string(theorypath_CMSobs+theoryfile_CMSobs).c_str());
  if(!str_theory_CMSobs.is_open()) {
    cout << "ERROR reading file " << theoryfile_CMSobs << endl;
  }
  j=0;

  while(!str_theory_CMSobs.eof() && j < nCMS ) {
    str_theory_CMSobs >> Mm >> xc;
    if(plotinTeV) Mt_CMSobs[j]=Mm/1000.;
    else Mt_CMSobs[j]=Mm;
    //Mt_BH[j]=Mm;
    xsect_CMSobs[j]=xc;  e0_CMSobs[j]=0;
    j++;
  }

  TGraph *gr_theory_CMSobs=new TGraph(nCMS, Mt_CMSobs, xsect_CMSobs);
  gr_theory_CMSobs->SetLineStyle(sty_theoryline);
  gr_theory_CMSobs->SetMarkerStyle(3);
  gr_theory_CMSobs->SetLineWidth(lwidth);
  gr_theory_CMSobs->SetLineColor(6);
  gr_theory_CMSobs->SetMarkerColor(6);




///CMS exp
  TString theorypath_CMSexp="theory/";
  TString theoryfile_CMSexp="zprimeCMSexpected.dat";


  double Mt_CMSexp[nCMS];
  double xsect_CMSexp[nCMS];
  double e0_CMSexp[nCMS];
  ifstream str_theory_CMSexp(string(theorypath_CMSexp+theoryfile_CMSexp).c_str());
  if(!str_theory_CMSexp.is_open()) {
    cout << "ERROR reading file " << theoryfile_CMSexp << endl;
  }
  j=0;

  while(!str_theory_CMSexp.eof() && j < nCMS ) {
    str_theory_CMSexp >> Mm >> xc;
    if(plotinTeV) Mt_CMSexp[j]=Mm/1000.;
    else Mt_CMSexp[j]=Mm;
    //Mt_BH[j]=Mm;
    xsect_CMSexp[j]=xc;  e0_CMSexp[j]=0;
    j++;
  }

  TGraph *gr_theory_CMSexp=new TGraph(nCMS, Mt_CMSexp, xsect_CMSexp);
  gr_theory_CMSexp->SetLineStyle(sty_theoryline);
  gr_theory_CMSexp->SetMarkerStyle(3);
  gr_theory_CMSexp->SetLineWidth(lwidth);
  gr_theory_CMSexp->SetLineColor(7);
  gr_theory_CMSexp->SetMarkerColor(7);





  //////Black hole theory
  const int Ntheory_BH=7;
  TString theorypath_BH="theory/";
  TString theoryfile_BH="blackholecrossection.dat";
  double Mt_BH[Ntheory_BH];
  double xsect_BH[Ntheory_BH];
  double e0_BH[Ntheory_BH];
  ifstream str_theory_BH(string(theorypath_BH+theoryfile_BH).c_str());
  if(!str_theory_BH.is_open()) {
    cout << "ERROR reading file " << theoryfile_BH << endl;
  } 
  j=0;
  while(!str_theory_BH.eof() && j < Ntheory_BH ) {
    str_theory_BH >> Mm >> xc;
    if(plotinTeV) Mt_BH[j]=Mm/1000.;
    else Mt_BH[j]=Mm;  
    //Mt_BH[j]=Mm;  
    xsect_BH[j]=xc;  e0_BH[j]=0;
    j++;
  }
  
  TGraph *gr_theory_BHc=new TGraph(Ntheory_BH, Mt_BH, xsect_BH);
  gr_theory_BHc->SetLineStyle(1);
  gr_theory_BHc->SetLineWidth(lwidth); 
  gr_theory_BHc->SetLineColor(kRed);

  ////////KKg theory
  TString theorypath_KKg="theory/";
  const int Ntheory_KKg=15;
 
  TString theoryfile_KKg="KKGluon6L1XS.dat";
  if (do8TeV) theoryfile_KKg="KKgtheory8TeV.dat";
  double Mt_KKg[Ntheory_KKg];
  double xsect_KKg[Ntheory_KKg];
  double xsect_errup_KKg[Ntheory_KKg];
  double xsect_errdown_KKg[Ntheory_KKg];
  double e0_KKg[Ntheory_KKg];
  ifstream str_theory_KKg(string(theorypath_KKg+theoryfile_KKg).c_str());
  if(!str_theory_KKg.is_open()) {
    cout << "ERROR reading file " << theoryfile_KKg << endl;
  } 
  j=0;
  while(!str_theory_KKg.eof() && j < Ntheory_KKg ) {
    str_theory_KKg >> Mm >> xc;
    if(plotinTeV) Mt_KKg[j]=Mm/1000.;
    else Mt_KKg[j]=Mm;  
    //Mt_KKg[j]=Mm;  
    xsect_KKg[j]=xc;  e0_KKg[j]=0;
    //Errors +5% and -5%. Gustaaf Brooijmans on ttbarres list 3 Aug 2011
    xsect_errup_KKg[j]=xc*0.05;
    xsect_errdown_KKg[j]=xc*0.05;    
    j++;
  }
 if (doKKg) {
 for(int i=0; i<Ntheory_KKg;++i) xsecmap[Mt_KKg[i]]=xsect_KKg[i];
  }
  TGraph *gr_theory_KKgc=new TGraph(Ntheory_KKg, Mt_KKg, xsect_KKg);
  // gr_theory_KKgc->SetLineStyle(1);
  // gr_theory_KKgc->SetLineWidth(lwidth); //3
  // gr_theory_KKgc->SetLineColor(kRed+4);
  gr_theory_KKgc->SetLineStyle(sty_theoryline);
  gr_theory_KKgc->SetMarkerStyle(1);
  gr_theory_KKgc->SetLineWidth(lwidth); //3
  gr_theory_KKgc->SetLineColor(col_theoryline);
  gr_theory_KKgc->SetMarkerColor(col_theoryline);

  //error bands  
  TGraph *gr_theory_KKg=getErrorBand(Ntheory_KKg, Mt_KKg, xsect_KKg, xsect_errdown_KKg, xsect_errup_KKg);
  // gr_theory_KKg->SetLineStyle(1);
  // gr_theory_KKg->SetLineWidth(lwidth);
  // gr_theory_KKg->SetLineColor(kRed+2);
  // gr_theory_KKg->SetFillColor(kRed+1);
  // gr_theory_KKg->SetFillStyle(3001);
  gr_theory_KKg->SetLineStyle(sty_theoryline);
  gr_theory_KKg->SetLineWidth(lwidth);
  gr_theory_KKg->SetLineColor(col_theoryline); //kRed+2);
  gr_theory_KKg->SetFillColor(col_theoryband);
  gr_theory_KKg->SetFillStyle(fill_theory);
  
  //////KKGWidth theory
  
  const int Ntheory_KKg_width=15;
  TString theoryfile_KKgWidth;
  if(selection.Contains("KKg-Width-1TeV")) 
    theoryfile_KKgWidth="KKGluon_Width_1TeV_theory.dat";
  else if(selection.Contains("KKg-Width-2TeV"))
    theoryfile_KKgWidth="KKGluon_Width_2TeV_theory.dat";
  else if(selection.Contains("KKg-Width-3TeV"))
    theoryfile_KKgWidth="KKGluon_Width_3TeV_theory.dat";
 
  else theoryfile_KKgWidth= "blah"; //KKGluon_Width_1TeV_theory.dat";

  cout<<"selection :"<<  selection<<endl;
  cout<<"theoryfile_KKgWidth :"<<  theoryfile_KKgWidth<<endl;


  double Mt_KKg_width[Ntheory_KKg_width];
  double xsect_KKg_width[Ntheory_KKg_width];

  ifstream str_theory_KKg_width(string(theorypath_KKg+theoryfile_KKgWidth).c_str());
  if(!str_theory_KKg_width.is_open()) {
    cout << "ERROR reading file " << theoryfile_KKg << endl;
  }
  j=0;
  while(!str_theory_KKg.eof() && j < Ntheory_KKg_width ) {
    str_theory_KKg_width >> Mm >> xc;
     Mt_KKg_width[j]=Mm;
    //Mt_KKg[j]=Mm;
    xsect_KKg_width[j]=xc;  
    j++;
  }

  if (doKKgWidth) {
    for(int i=0; i<Ntheory_KKg_width;++i) xsecmap[Mt_KKg_width[i]]=xsect_KKg_width[i];
  }
  TGraph *gr_theory_KKg_widthc=new TGraph(Ntheory_KKg_width, Mt_KKg_width, xsect_KKg_width);
  // gr_theory_KKgc->SetLineStyle(1);
  // gr_theory_KKgc->SetLineWidth(lwidth); //3
  // gr_theory_KKgc->SetLineColor(kRed+4);
  gr_theory_KKg_widthc->SetLineStyle(sty_theoryline);
  gr_theory_KKg_widthc->SetMarkerStyle(1);
  gr_theory_KKg_widthc->SetLineWidth(lwidth); //3
  gr_theory_KKg_widthc->SetLineColor(col_theoryline);
  gr_theory_KKg_widthc->SetMarkerColor(col_theoryline);

  cout<<"col_theoryline: "<<col_theoryline<<endl;


  ////////KKGrav theory
  TString theorypath_KKGrav="theory/";
  const int Ntheory_KKGrav=13;
  TString theoryfile_KKGrav="KKGrav_8TeV.dat";
  double Mt_KKGrav[Ntheory_KKGrav];
  double xsect_KKGrav[Ntheory_KKGrav];
  double xsect_errup_KKGrav[Ntheory_KKGrav];
  double xsect_errdown_KKGrav[Ntheory_KKGrav];
  double e0_KKGrav[Ntheory_KKGrav];
  ifstream str_theory_KKGrav(string(theorypath_KKGrav+theoryfile_KKGrav).c_str());
  if(!str_theory_KKGrav.is_open()) {
    cout << "ERROR reading file " << theoryfile_KKGrav << endl;
  } 
  j=0;
  cout<<"kkgrav theory: \n";
  while(!str_theory_KKGrav.eof() && j < Ntheory_KKGrav ) {
    str_theory_KKGrav >> Mm >> xc;
    cout<<Mm<<" "<<xc<<endl;
    if(plotinTeV) Mt_KKGrav[j]=Mm/1000.;
    else Mt_KKGrav[j]=Mm;  
    xsect_KKGrav[j]=xc;  e0_KKGrav[j]=0;
    //Errors +5% and -5%. Assumption -- needs to be checked
    xsect_errup_KKGrav[j]=xc*0.05;
    xsect_errdown_KKGrav[j]=xc*0.05;    
    j++;
  }

  TGraph *gr_theory_KKGravc=new TGraph(Ntheory_KKGrav, Mt_KKGrav, xsect_KKGrav);
  // gr_theory_KKGravc->SetLineStyle(1);
  // gr_theory_KKGravc->SetLineWidth(lwidth); //3
  // gr_theory_KKGravc->SetLineColor(kRed+4);
  gr_theory_KKGravc->SetLineStyle(sty_theoryline);
  gr_theory_KKGravc->SetMarkerStyle(1);
  gr_theory_KKGravc->SetLineWidth(lwidth); //3
  gr_theory_KKGravc->SetLineColor(col_theoryline);
  gr_theory_KKGravc->SetMarkerColor(col_theoryline);
  //error bands  
  TGraph *gr_theory_KKGrav=getErrorBand(Ntheory_KKGrav, Mt_KKGrav, xsect_KKGrav, xsect_errdown_KKGrav, xsect_errup_KKGrav);
  // gr_theory_KKGrav->SetLineStyle(1);
  // gr_theory_KKGrav->SetLineWidth(lwidth);
  // gr_theory_KKGrav->SetLineColor(kRed+2);
  // gr_theory_KKGrav->SetFillColor(kRed+1);
  // gr_theory_KKGrav->SetFillStyle(3001);
  gr_theory_KKGrav->SetLineStyle(sty_theoryline);
  gr_theory_KKGrav->SetMarkerStyle(1);
  gr_theory_KKGrav->SetLineWidth(lwidth); //3
  gr_theory_KKGrav->SetLineColor(col_theoryline);
  gr_theory_KKGrav->SetMarkerColor(col_theoryline);

  //Make TGraphs of limits and draw

  //legends
  TString selname=getName(selection);
  //TString legname=selname;
  TString legname="";
  if(channel == "ele") legname+= ", ele channel";
  else if(channel == "muo") legname+= ", muo channel";

  //if(doStat)  { legname+=".  Stat. only"; }
  //else {  legname+=".  Syst.+stat.";  }
  if(path.Contains("dRorChi2")) {
    if(doStat)  { legname+="dRmin reconstruction"; }
    else {  legname+="#chi ^{2} reconstruction";  }
  }
  else {
    if(doStat)  { legname+="Stat. only errors"; }
    //else {  legname+="Syst.+stat. errors";  }
  }

  TString lt_exp="Exp. 95% CL upper limit";
  TString lt_obs="Obs. 95% CL upper limit";
  TString lt_68="Exp. 1 #sigma uncertainty";
  TString lt_95="Exp. 2 #sigma uncertainty";
  TString lt_th="Leptophobic Z' (LO x 1.3)";
  TString lt_th_2="Leptophobic Z'(2%) (LO x 1.3)";
  TString lt_th_3="Leptophobic Z'(3%) (LO x 1.3)";
  TString lt_th_BH="Quantum Black Hole";
  TString lt_th_KKg="Kaluza-Klein gluon (LO)";
  TString lt_th_KKGrav="Kaluza-Klein graviton";

  //axis limits
  double x_ax_xs[2]={510., 2000.};
  if(doBlackHole) { x_ax_xs[0]=750;  x_ax_xs[1]=2500; } 
  if(doKKg) { 
    x_ax_xs[0]=500; //795; //805; //715
    x_ax_xs[1]=3000; //1885; //1893;
    //x_ax_xs[1]=4000;
  }
  
  if(doKKgWidth) {
    x_ax_xs[0]=10; //795; //805; //715
    x_ax_xs[1]=40; //1885; //1893;
    //x_ax_xs[1]=4000;
  }

  if(doZprime) { 
    //x_ax_xs[0]=700; //715; //710;  
    // x_ax_xs[0]=618; //617
    // x_ax_xs[1]=1882; //1885; 
    x_ax_xs[0]=400; //617
    x_ax_xs[1]=3000; //1885;
    // x_ax_xs[1]=1000; //1885;
    //  x_ax_xs[1]=2100; //1885; 
    if(Nmasses==9 || Nmasses==8) {
      x_ax_xs[0]=707;
      x_ax_xs[1]=2795; 
    }
  } 
  if(doKKGrav) {
    x_ax_xs[0]=400.; //510.;
    x_ax_xs[1]=2500; //1285;
  }
  if(doHH) {
    x_ax_xs[0]=400.; //510.;
    x_ax_xs[1]=3000; //1285;
  }
  
  if(plotinTeV && ! doKKgWidth) {
    x_ax_xs[0]=x_ax_xs[0]/1000.;
    x_ax_xs[1]=x_ax_xs[1]/1000.;
  }

  double y_ax_xs[2];
  y_ax_xs[0]=0; y_ax_xs[1]=101;
  //if(doStat) {
  if(doLogy) {
     if(doZprime)    { y_ax_xs[0]=5e-3;  y_ax_xs[1]=5e3; } //9e-3
     // if(doZprime)    { y_ax_xs[0]=1e-1;  y_ax_xs[1]=1e3; }
    if(doratio && doZprime)     { y_ax_xs[0]=5e-3;  y_ax_xs[1]=5e+3; }
    //  if(doZprime)    { y_ax_xs[0]=5e-3;  y_ax_xs[1]=5e3; }
    //if(doKKGrav)    { y_ax_xs[0]=3e-2;  y_ax_xs[1]=2e3; }
    if(doKKGrav)    { y_ax_xs[0]=3e-4;  y_ax_xs[1]=2e3; }
    if(doHH)    { y_ax_xs[0]=5e-3;  y_ax_xs[1]=5e3; }
    if(doKKg )       { y_ax_xs[0]=1e-2;  y_ax_xs[1]=6e3; } //2e-2
    if(doKKgWidth )       { y_ax_xs[0]=1e-3;  y_ax_xs[1]=6e3; } //2e-2
    if(doratio && doKKg)     { y_ax_xs[0]=5e-3;  y_ax_xs[1]=5e+3; }
    if(doBlackHole) { y_ax_xs[0]=1e-1;  y_ax_xs[1]=4e3; }
  }
  else {
    if(doZprime)    { y_ax_xs[0]=7e-1;  y_ax_xs[1]=1e3; }
    if(doHH)    { y_ax_xs[0]=7e-1;  y_ax_xs[1]=1e3; }
    if(doKKGrav)    { y_ax_xs[0]=3e-2;  y_ax_xs[1]=8e2; }
    if(doKKg )       { y_ax_xs[0]=2e-2;  y_ax_xs[1]=6e2; cout<<"test 123!! \n";}
    if(doBlackHole) { y_ax_xs[0]=1.5e1; y_ax_xs[1]=4e3; }     
  }



  //axis titles
  TString x_title="Z' mass ";
  TString y_title="#sigma_{Z'} #times BR(Z'#rightarrow t#bar{t}) [pb]";
  
  if(doBlackHole) {
    x_title="QBH mass threshold ";
    y_title="#sigma_{QBH} [pb]";
  }
  else if(doKKg) {
    x_title="g_{KK} mass ";
    y_title="#sigma_{g_{KK}} #times BR(g_{KK}#rightarrow t#bar{t}) [pb]";
  }
  else if(doKKgWidth) {
    x_title="g_{KK} Width (%)";
    y_title="#sigma_{g_{KK}} #times BR(g_{KK}#rightarrow t#bar{t}) [pb]";
  }
  else if(doKKGrav) {
    x_title="G_{KK} mass ";
    y_title="#sigma_{G_{KK}} #times BR(G_{KK}#rightarrow t#bar{t}) [pb]";
  }
  else if(doHH) {
    x_title="scalar resonance mass ";
    y_title="#sigma_{scalar res.} #times BR(scalar res.#rightarrow t#bar{t}) [pb]";
  }
  if( !doKKgWidth){
    if(plotinTeV  ) x_title+="[TeV]";
    else x_title+="[GeV]";
  }
  if (doratio) y_title="95% CL Limit on #mu";
  //Draw
  int can_w=2400, can_h=2400;
  if(doRectangular) {
    //can_h=2100; can_w= 2970;
    can_h=1200; can_w=1600;
  }

  TCanvas* cv_xs = new TCanvas("cv_xs","Cross section",can_w,can_h);
  //   TCanvas* cv_xs = new TCanvas("cv_xs","Cross section",50,50,2400,2400);
  if(doBlackHole || doLogy) {
    gPad->SetLogy(1);
  }

  TGraph *gr_ax_xs = new TGraph(2, x_ax_xs, y_ax_xs);
  //   gr_ax_xs->SetMarkerStyle(1);
  gr_ax_xs->SetMarkerColor(kWhite);
  gr_ax_xs->SetTitle("");
  if(doBlackHole) gr_ax_xs->GetXaxis()->SetNdivisions(505, 'x');

  TGraph *gr_expected_xs = new TGraph(Nmasses, signal_mass, limitMean);
  //   gr_expected_xs->SetMarkerStyle(1);
  gr_expected_xs->SetLineColor(kBlack);
  gr_expected_xs->SetLineStyle(sty_exp); //2
  gr_expected_xs->SetLineWidth(lwidth);

  TGraph *gr_observed_xs = new TGraph(Nmasses, signal_mass, limitObserved);
  gr_observed_xs->SetMarkerStyle(20);
  gr_observed_xs->SetLineColor(kBlack);
  gr_observed_xs->SetLineWidth(lwidth);

  TGraph *gr_errorBand68_xs=getErrorBand(Nmasses, signal_mass, limHigh68, limLow68);
  TGraph *gr_errorBand95_xs=getErrorBand(Nmasses, signal_mass, limHigh95, limLow95);
  TGraph *gr_expected_xs_r;
  TGraph *gr_observed_xs_r;
  TGraph *gr_errorBand68_xs_r;
  TGraph *gr_errorBand95_xs_r;
  TLine* l1;
  if (doratio) {
  for(int i=0;i<Nmasses;++i) {
    if (doratio) {
      limitMean_r[i]= limitMean[i]/xsecmap[signal_mass[i]];
      limHigh68_r[i]=limHigh68[i]/xsecmap[signal_mass[i]];
      limLow68_r[i] = limLow68[i]/xsecmap[signal_mass[i]];
      limHigh95_r[i]= limHigh95[i]/xsecmap[signal_mass[i]];
      limLow95_r[i] = limLow95[i]/xsecmap[signal_mass[i]];
      limitObserved_r[i]=limitObserved[i]/xsecmap[signal_mass[i]];
    }
  }
  gr_expected_xs_r = new TGraph(Nmasses, signal_mass, limitMean_r);
  //   gr_expected_xs->SetMarkerStyle(1);
  gr_expected_xs_r->SetLineColor(kBlack);
  gr_expected_xs_r->SetLineStyle(2);
  gr_expected_xs_r->SetLineWidth(lwidth);
  gr_observed_xs_r = new TGraph(Nmasses, signal_mass, limitObserved_r);
  gr_observed_xs_r->SetMarkerStyle(20);
  gr_observed_xs_r->SetLineColor(kBlack);
  gr_observed_xs_r->SetLineWidth(lwidth);
  gr_errorBand68_xs_r=getErrorBand(Nmasses, signal_mass, limHigh68_r, limLow68_r);
  gr_errorBand95_xs_r=getErrorBand(Nmasses, signal_mass, limHigh95_r, limLow95_r);
    gr_errorBand68_xs_r->SetFillColor(kGreen);
    gr_errorBand95_xs_r->SetFillColor(kYellow);     
    gr_errorBand68_xs_r->SetLineColor(kGreen);
    gr_errorBand95_xs_r->SetLineColor(kYellow);
  }

  if(do95Lim) {
    gr_errorBand68_xs->SetFillColor(kGreen);
    gr_errorBand95_xs->SetFillColor(kYellow);     
    gr_errorBand68_xs->SetLineColor(kGreen);
    gr_errorBand95_xs->SetLineColor(kYellow);
  }
  else {
    gr_errorBand68_xs->SetFillColor(kYellow);
    gr_errorBand95_xs->SetFillColor(kWhite);
    gr_errorBand68_xs->SetLineColor(kBlack);
    gr_errorBand95_xs->SetLineColor(kWhite);
  }

  gr_ax_xs->Draw("AP");
  gr_ax_xs->GetXaxis()->SetRangeUser(x_ax_xs[0],x_ax_xs[1]);
  cv_xs->Update();
  if (!doratio) {
    if(do95Lim) gr_errorBand95_xs->Draw("FL");
    gr_errorBand68_xs->Draw("FL");
    if(doBlackHole) {
      gr_theory_BHc->Draw("L");
    }
    else if(doKKg) {
      if(doDrawTheoryUncert)  gr_theory_KKg->Draw("F");    
      gr_theory_KKgc->Draw("L");
    }
    else if(doKKgWidth) {
      //cout<<"plotting KKgWidth theory line \n";
      //gr_theory_KKg_widthc->SetLineColor(col_theoryline);
      gr_theory_KKg_widthc->Draw("L");
    }

    else if(doKKGrav) {
      if(doDrawTheoryUncert)  gr_theory_KKGrav->Draw("F");    
      gr_theory_KKGravc->Draw("L");
      cout<<"drawing kkgrav theory line \n";
    }
    else if( doHH ){
      std::cout << "No theory lines available for HH" << std::endl;
    }
    else {
      if(doDrawTheoryUncert) gr_theory_DW->Draw("F");
      gr_theory_DWc->Draw("L");   

      if(doZprime3p){
	gr_theory_Z3p->Draw("Lsame");
	gr_theory_Z2p->Draw("Lsame");
      }
      
      if(doCMS){
        gr_theory_CMSobs->Draw("Lsame");
        gr_theory_CMSexp->Draw("Lsame");
      }


    }
    gr_expected_xs->Draw("PL");
    gr_ax_xs->GetXaxis()->SetTitle(x_title);
    gr_ax_xs->GetYaxis()->SetTitle(y_title);
 
  }
  else {
   
    if(do95Lim) gr_errorBand95_xs_r->Draw("FL");
    gr_errorBand68_xs_r->Draw("FL");
    gr_expected_xs_r->Draw("PL");
    if (doObserved) gr_observed_xs_r->Draw("PL");  
    gr_ax_xs->GetXaxis()->SetTitle(x_title);
    gr_ax_xs->GetYaxis()->SetTitle(y_title);
    l1=new TLine(x_ax_xs[0],1,x_ax_xs[1],1);
    l1->SetLineColor(4);
    l1->SetLineStyle(2);
    l1->Draw();
  }

  TLegend *leg_limvsmass_xs = new TLegend(lxlow, lylow, lxup, lyup);
  if(doObserved ) {
    if (!doratio) gr_observed_xs->Draw("PL");
    leg_limvsmass_xs->AddEntry(gr_observed_xs, lt_obs,"pl");
    leg_limvsmass_xs->Draw();
  }
  leg_limvsmass_xs->AddEntry(gr_expected_xs, lt_exp,"pl");      
  leg_limvsmass_xs->AddEntry(gr_errorBand68_xs, lt_68,"lf");
  if(do95Lim) leg_limvsmass_xs->AddEntry(gr_errorBand95_xs, lt_95,"lf");
  if (!doratio) {
 if(doDrawTheoryUncert) {
    if(doBlackHole) leg_limvsmass_xs->AddEntry(gr_theory_BHc, lt_th_BH,"lp");
    //else if(doKKg) leg_limvsmass_xs->AddEntry(gr_theory_KKgc, lt_th_KKg,"lp");
    else if(doKKg || doKKgWidth) leg_limvsmass_xs->AddEntry(gr_theory_KKg, lt_th_KKg,"lfp");
    else if(doKKGrav) leg_limvsmass_xs->AddEntry(gr_theory_KKGrav, lt_th_KKGrav,"lfp");
    else leg_limvsmass_xs->AddEntry(gr_theory_DW, lt_th,"lfp");
 }
  else {
    if(doBlackHole) leg_limvsmass_xs->AddEntry(gr_theory_BHc, lt_th_BH,"l");
    else if(doKKg || doKKgWidth) leg_limvsmass_xs->AddEntry(gr_theory_KKg, lt_th_KKg,"l");
    else if(doKKGrav) leg_limvsmass_xs->AddEntry(gr_theory_KKGrav, lt_th_KKGrav,"l");
    else if(doHH) std::cout << "No legend entry to be added for HH theory graph" << std::endl;
    else {
      leg_limvsmass_xs->AddEntry(gr_theory_DW, lt_th,"lp");
      if(doZprime3p){
	leg_limvsmass_xs->AddEntry(gr_theory_Z2p, lt_th_2,"lp");
	leg_limvsmass_xs->AddEntry(gr_theory_Z3p, lt_th_3,"lp");
      }
      if(doCMS){
        leg_limvsmass_xs->AddEntry(gr_theory_CMSobs, "CMS obs","lp");
        leg_limvsmass_xs->AddEntry(gr_theory_CMSexp, "CMS_exp","lp");
      }


    }
  }
  }
  else {
    if(doKKg) leg_limvsmass_xs->AddEntry(l1, lt_th_KKg,"lfp");
    else leg_limvsmass_xs->AddEntry(l1, lt_th,"l");
  }
  leg_limvsmass_xs->SetHeader(legname);
  leg_limvsmass_xs->SetFillColor(0);
  leg_limvsmass_xs->SetBorderSize(0);
  leg_limvsmass_xs->Draw();

  //ATLAS labels
  double a_x=lxlow+0.02, a_y=lylow-0.05; //, l_y=lyup-0.11, l_x=0.19;
  if(doAtlasLabel) ATLASLabel(a_x, a_y, cv_xs,true);
  double l_y=lyup-0.17, l_x=0.19;
  ATLASCMEIntL(l_x,l_y,LUMI,cv_xs);

  gPad->RedrawAxis();

  if(doZprime){
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(kBlue+2); 
    tex->SetTextSize(0.065);
    tex->DrawLatex(600,1,"excluded(CDF)");
  }
  
  //Save plots with appropriate names
  TString signal="";
  TString tool="";
  TString stsy="";
  TString sel="";
  TString linlog="";

  if(selection.Contains("Zprime-lim")) signal="_zprime";
  else if(selection.Contains("HH-lim")) signal="_hh";
  else if(selection.Contains("KKg-lim")) signal="_kkg";
  else if(selection.Contains("KKg-Width-1TeV-lim")) signal="_kkg_width_1tev";
  else if(selection.Contains("KKg-Width-2TeV-lim")) signal="_kkg_width_2tev";
  else if(selection.Contains("KKg-Width-3TeV-lim")) signal="_kkg_width_3tev";
  else if(selection.Contains("KKGrav-lim")) signal="_kkgrav";
  
  if(doKKgWidth){
  
    const char* label="";
    if( signal=="_kkg_width_1tev" ) label="m_{g_{KK}} = 1 TeV";
    else if( signal=="_kkg_width_2tev" ) label="m_{g_{KK}} = 2 TeV";
    else if( signal=="_kkg_width_3tev" ) label="m_{g_{KK}} = 3 TeV";
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(1); 
    tex->SetTextSize(0.045);
    tex->DrawLatex(l_x,l_y-0.1,label);
  }


  if(path.Contains("gustaaf")) tool="_cls";
  else if(path.Contains("wuppertal")) tool="_bayes";
  else if(path.Contains("calvet")) tool="_gcbayes";

  if(selection.Contains("syst")) stsy="_systconv";
  else if(selection.Contains("stat")) stsy="_stat";

  if(selection.Contains("4j")) sel="_4j";
  else if(selection.Contains("dr")) sel="_dr";
  else if(selection.Contains("nc")) sel="_x2nc";
  else if(selection.Contains("chi2cal")) sel="_x2c";  
  else if(selection.Contains("pretag")) {
    if(path.Contains("Combined")||path.Contains("combined")) sel="_combined";
    else if(path.Contains("Boosted")||path.Contains("boosted")) sel="_boosted";
    else if(path.Contains("Resolved")||path.Contains("resolved")) sel="_resolved";
    else sel="_pretag";  
  }

  if(!doLogy) linlog="_lin";

  TString plotname="limits"+signal+tool+stsy+"_obs_exp"+sel+linlog;
  if (!doratio) {
    cv_xs->Print(plotname+".pdf","pdf");
    cv_xs->Print(plotname+".eps","eps");
    cv_xs->Print(plotname+".png","png");
  }
  else  {
    cv_xs->Print(plotname+"_ratio.pdf","pdf");
    cv_xs->Print(plotname+"_ratio.eps","eps");
    cv_xs->Print(plotname+"_ratio.png","png");
  }

  //========================================
  //Print limits in a table
  //========================================
  
  //Define separator
  TString lend="", sep=""; 
  if(printTex) {
    lend="\\\\"; sep="&"; 
  }
  
  //Print column labels
  cout << "The limits (pb)" << endl;
  cout << " Mass " << sep << "  Obs.  " << sep << "  Exp. " << sep 
       << "   -1std " << sep << "  +1std ";
  
  if(print2sigma) cout << sep << "  -2std " << sep << "  +2std";
  cout<< lend << endl;
  
  const int massprec=(plotinTeV ? 2: 0);
  
  //Print limit numbers
  for(int i=0; i<Nmasses; i++) {
    cout << setiosflags(ios::fixed) << setprecision(massprec)
	 << setw(5) << signal_mass[i] 
	 << setiosflags(ios::fixed) << setprecision(3)
	 << sep << setw(8) << limitObserved[i] << sep << setw(8) << limitMean[i]
	 << sep << setw(8) << limLow68[i] << sep << setw(8) << limHigh68[i];
    if(print2sigma) {
      cout<< sep << setw(8) << limLow95[i] << sep << setw(8) << limHigh95[i];
    }
    cout << lend << endl;
  }
    if (doratio) {
      cout<<"limit/theory xsection "<<endl;
      for(int i=0; i<Nmasses; i++) {
	cout << setiosflags(ios::fixed) << setprecision(massprec)
	     << setw(5) << signal_mass[i] 
	     << setiosflags(ios::fixed) << setprecision(2)
	     << sep << setw(8) << limitObserved_r[i]<< sep << setw(8) << limitMean_r[i]
	     << sep << setw(8) << limLow68_r[i] << sep << setw(8) << limHigh68_r[i];
	if(print2sigma) {
	  cout<< sep << setw(8) << limLow95_r[i] << sep << setw(8) << limHigh95_r[i];
	}
	cout << lend << endl;
      }
    }

  //Print theory xsec
  if(doPrintXsec) {
    cout <<"Theory cross section ";
    if(doBlackHole) {
      cout << "black hole" << endl;
      for(int i=0; i<Ntheory_BH;++i) {
	if((Mt_BH[i]/100.) == int(Mt_BH[i])/100) {
	  cout << setw(6) << setiosflags(ios::fixed) << setprecision(massprec) << Mt_BH[i] 
	       << sep << setw(7) << setiosflags(ios::fixed);
	  cout << setprecision(4);
	  cout << xsect_BH[i] << lend << endl;
	}
      }
    }
    if(doKKg) {
      cout << "kkg" << endl;
      for(int i=0; i<Ntheory_KKg;++i) {
	//	if((Mt_KKg[i]/100.) == int(Mt_KKg[i])/100) 
      {
	  cout << "$g_{\\rm KK}$, $m_{g_{\\rm KK}}$ = "
	       << setw(6) << setiosflags(ios::fixed) << setprecision(massprec) 
	       << Mt_KKg[i] << " \\GeV " << sep 
	       << setw(7) << setiosflags(ios::fixed);
	  // if(Mt_KKg[i]<901) cout << setprecision(1);
	  // else if( Mt_KKg[i]<1901) cout << setprecision(2);
	  // else if( Mt_KKg[i]<2500) cout << setprecision(3);
	  // else cout << setprecision(4);
	   cout << xsect_KKg[i] << lend << endl;
	}
      }
    }
    if(doZprime) {
      cout << "Z'" << endl;
      for(int i=0; i<Ntheory_DW;++i) {
	//	if((Mt_DW[i]/100.) == int(Mt_DW[i])/100) 
         {
	  cout << "Topcolor $Z_t'$, $m_{Z_t'}$ = "
	       << setw(6) << setiosflags(ios::fixed) << setprecision(massprec) 
	       << Mt_DW[i] << " \\GeV " << sep 
	       << setw(7) << setiosflags(ios::fixed);
	  // if(Mt_DW[i]<941) cout << setprecision(2);
	  // else if(Mt_DW[i]<1241) cout << setprecision(3);
	  // else cout << setprecision(5);
          if(Mt_DW[i]>2) cout << setprecision(3);
	  cout << xsect_DW[i] << lend << endl;
	}
      }
    }
    if(doHH) {
      cout << "HH" << endl;
      cout << "No HH theory cross section" << endl;
    }
    if(doKKGrav) {
      cout << "KKG" << endl;
      for(int i=0; i<Ntheory_KKGrav;++i) {
	if((Mt_DW[i]/100.) == int(Mt_DW[i])/100) {
	  cout << "Kaluza--Klein Graviton, $m_{GKK}$ = "
	       << setw(6) << setiosflags(ios::fixed) << setprecision(massprec) 
	       << Mt_DW[i] << " \\GeV " << sep 
	       << setw(7) << setiosflags(ios::fixed);
	  if(Mt_KKGrav[i]<941) cout << setprecision(2);
	  else if(Mt_KKGrav[i]<1241) cout << setprecision(3);
	  else cout << setprecision(4);
	  cout << xsect_KKGrav[i] << lend << endl;
	}
      }
    }
  }


}
