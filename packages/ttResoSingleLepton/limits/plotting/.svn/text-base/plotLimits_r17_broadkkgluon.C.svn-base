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


void plotLimits_r17_broadkkgluon(bool doRectangular = false)
{ 
  int can_w=800, can_h=800;
  if(doRectangular) {
    //can_h=2100; can_w= 2970;
    can_h=600; can_w=800;
  }
  gStyle->SetPadLeftMargin(0.16); // 0.21// 0.18
  gStyle->SetPadRightMargin(0.05); //0.08
  gStyle->SetPadTopMargin(0.01); //0.07
  gStyle->SetPadBottomMargin(0.2);
    gStyle->SetLineWidth(3); // width of ticks
  gStyle->SetPadTickX(1); //0//1:ticks on upper,2: ticks+labels on upper xaxis
  gStyle->SetPadTickY(1); //0

  TCanvas* cv_xs = new TCanvas("cv_xs","Cross section",can_w,can_h);
  TPad *pad1 = new TPad("1","1",0,0,1,0.35);
  pad1->Draw(); 
  gStyle->SetPadTopMargin(0.0); 
  gStyle->SetPadBottomMargin(0.0);
  TPad *pad2 = new TPad("2","2",0,0.35,1,0.67); 
  pad2->Draw();

  gStyle->SetPadTopMargin(0.05); 
  gStyle->SetPadBottomMargin(0.01);
  TPad *pad3 = new TPad("3","3",0,0.67,1,1.0);
  pad3->Draw();
  pad1->cd();
  
  const char* dir = "/data/atlas/atlasdata2/behr/TTbarResonanceSearch/LimitSetting/results_v067_06Oct_smoothed_data/Results/Limits/LimitNumbers/Combined";
  
  plotLimits_r17(dir,"KKg-Width-1TeV-lim_pretag_syst","multi","m_{g_{KK}} = 1 TeV",7,true,true,true, false, false);


  pad2->cd();
  plotLimits_r17(dir,"KKg-Width-2TeV-lim_pretag_syst","multi","m_{g_{KK}} = 2 TeV",7,true,true,true, true, false, true);
  pad3->cd();
  plotLimits_r17(dir,"KKg-Width-3TeV-lim_pretag_syst","multi","m_{g_{KK}} = 3 TeV",7,true,true,true, true, true);
  
  cv_xs->SaveAs("AllWidthScanLimits.eps");
  cv_xs->SaveAs("AllWidthScanLimits.png");
}

void plotLimits_r17(TString path, TString selection, TString channel, TString label, const int Nmasses, bool doObserved, bool doAtlasLabel, bool doLogy, bool same, bool first, bool second=false) {
  bool debug=true;


  const int doPrintXsec=true;

  bool doZp13factor=true;
  
  bool plotinTeV=true;
  
  const bool doDrawTheoryUncert=false;

  double LUMI = 20.3;


  //bool printTex = false; 
  bool printTex = true; //print limits as a TeX table
  bool print2sigma= true; //prints only 1 sigma limits, not 2 sigma
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


  cout<<"selection: " << selection<<endl;
  const int lwidth=3;

  assert(channel == "ele" || channel == "muo" ||channel == "multi");

  //Automatic settings -- what kind of limits are we plotting?
  bool doBlackHole=false;
  bool doKKg=false;
  bool doZprime=false;
  bool doZprime3p=true;  //add 3% theory line
  bool doKKGrav=false;
  bool doratio=false;
  
  bool doKKgWidth=false;

  bool doCMS=false;
  
  if(selection.Contains("BH-lim")  || selection.Contains("qbh")) {
    doBlackHole=true; }
  else if(selection.Contains("KKg-lim") || selection.Contains("kkg")) {
    doKKg=true; }
  else if(selection.Contains("KKg-Width")) {
    doKKgWidth=true; }
  else if(selection.Contains("KKGrav-lim")) {
    doKKGrav=true; 
    cout<<"KKGrav test \n";
  }
  else doZprime=true;

  if(doKKgWidth) plotinTeV=false;

  bool doStat=true; //Stat only limits
  if(selection.Contains("_sys") || path.Contains("sys")) doStat=false;



  //Legend settings
  //double lxlow=0.45, lylow= 0.6, lxup=0.92, lyup=0.92;
  double lxlow=0.45, lylow= 0.55, lxup=0.92, lyup=0.8;
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

  ifstream str_exp(string(path+exp_filename).c_str());
  if(!str_exp.is_open()) {
    cout << "ERROR reading file " << path+exp_filename << endl;
  }

  ifstream str_e2x(string(path+e2x_filename).c_str());
  if(!str_e2x.is_open() && e2x_filename!="null") {
    cout << "ERROR reading file " << path+e2x_filename << endl;
  }

  //Extract limits from files
  double EventsObs[Nmasses], EventsExp[Nmasses], Events2pb[Nmasses], ee;
  double EventsExp_m95[Nmasses], EventsExp_m68[Nmasses];
  double EventsExp_p68[Nmasses], EventsExp_p95[Nmasses]; 
  double signal_mass[Nmasses];
  double m95, m68, p68, p95;
  int j=0, Mindex;
  cout << "Observed: ";
  j=0;
  while(!str_obs.eof() && j < Nmasses ) {
    if(selection.Contains("-lim")) {
      str_obs >> Mindex >> ee;
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


  cout << "Expected mean: ";
  j=0;
  if(selection.Contains("-lim_")) { //order: M, m, -1s, +1s, -2s, +2s
    while(!str_exp.eof() && j < Nmasses ) {
      str_exp >> Mindex >> ee >> m68 >> p68 >> m95 >> p95;
      
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


  /////////////////////////////////////////////////
  /// Draw plots
  //Colours
  int col_theoryband=kRed+1, col_theoryline=kRed+2;
  if(doDrawTheoryUncert) col_theoryline=kRed+4;
  int col_68=kGreen, col_95=kYellow;

  //Style
  int sty_theoryline=7; //1-solid, 7-dashed
  int fill_theory=3001;
  int sty_exp=2;

 
 
  //////KKGWidth theory

  const int Ntheory_KKg_width=9;
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


  int Mt_KKg_width[Ntheory_KKg_width];
  double Mt_KKg_width_double[Ntheory_KKg_width];
  double xsect_KKg_width[Ntheory_KKg_width];

  TString theorypath_KKg = "./";

  ifstream str_theory_KKg_width(string(theorypath_KKg+theoryfile_KKgWidth).c_str());
  if(!str_theory_KKg_width.is_open()) {
    cout << "ERROR reading file " << theoryfile_KKg << endl;
  }
  j=0;
  int Mm;
  float xc;
  while(!str_theory_KKg_width.eof() && j < Ntheory_KKg_width ) {
    str_theory_KKg_width >> Mm >> xc;
     Mt_KKg_width[j]=Mm;    
     Mt_KKg_width_double[j]=(double) Mm;
    xsect_KKg_width[j]=xc;  
    j++;
  }

  double xsecmap[Ntheory_KKg_width];
  xsect_KKg_width[0]/=1.3;
  for(int i=0; i<Ntheory_KKg_width;++i) xsecmap[i]=xsect_KKg_width[i];
  TGraph *gr_theory_KKg_widthc2=new TGraph(2, Mt_KKg_width_double, xsect_KKg_width);
 TGraph *gr_theory_KKg_widthc=new TGraph(Ntheory_KKg_width-1, &Mt_KKg_width_double[1], &xsect_KKg_width[1]);

  // gr_theory_KKgc->SetLineStyle(1);
  // gr_theory_KKgc->SetLineWidth(lwidth); //3
  // gr_theory_KKgc->SetLineColor(kRed+4);
  gr_theory_KKg_widthc->SetLineStyle(sty_theoryline);
  gr_theory_KKg_widthc->SetMarkerStyle(1);
  gr_theory_KKg_widthc->SetLineWidth(lwidth); //3
  gr_theory_KKg_widthc->SetLineColor(col_theoryline);
  gr_theory_KKg_widthc->SetMarkerColor(col_theoryline);

  cout<<"col_theoryline: "<<col_theoryline<<endl;


  cout << " done with broad KK gluon " << endl;



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
  
  TString lt_ATLAS_width="ATLAS: width/mass = 15.3%";
  TString lt_CMS_width="CMS: width/mass = 20%";

  //axis limits
  double x_ax_xs[2]={510., 2000.};

  if(doKKgWidth) {
    x_ax_xs[0]=10; //795; //805; //715
    x_ax_xs[1]=40; //1885; //1893;
    //x_ax_xs[1]=4000;
  }


 

  if(plotinTeV && ! doKKgWidth) {
    x_ax_xs[0]=x_ax_xs[0]/1000.;
    x_ax_xs[1]=x_ax_xs[1]/1000.;
  }

  double y_ax_xs[2];
  y_ax_xs[0]=0; y_ax_xs[1]=101;
  //if(doStat) {
  if(doLogy) {
    if(doKKgWidth && !same)       { y_ax_xs[0]=5e-2;  y_ax_xs[1]=20; } //2e-2 
    if(doKKgWidth && same)      { y_ax_xs[0]=5e-3; y_ax_xs[1]=2;}
   
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
    if (!same) x_title="g_{KK} Width/Mass [%]";
    else  x_title="";
    y_title="#sigma_{g_{KK}} #times BR(g_{KK}#rightarrow t#bar{t}) [pb]";
  }

  else if(doKKGrav) {
    x_title="G_{KK} mass ";
    y_title="#sigma_{G_{KK}} #times BR(G_{KK}#rightarrow t#bar{t}) [pb]";
  }
  if( !doKKgWidth)
    if(plotinTeV  ) x_title+="[TeV]";
    else x_title+="[GeV]";
  if (doratio) y_title="95% CL Limit on #mu";
  //Draw
 
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

  if (!same) gr_ax_xs->Draw("AP");
  else gr_ax_xs->Draw("AP");
  gr_ax_xs->GetXaxis()->SetRangeUser(x_ax_xs[0],x_ax_xs[1]);
  gr_ax_xs->GetXaxis()->SetTitleOffset(0.9); 
  if (same) gr_ax_xs->GetXaxis()->SetLabelSize(0);
  else  gr_ax_xs->GetXaxis()->SetLabelSize(0.08);

  if (same) gr_ax_xs->GetYaxis()->SetLabelSize(0.09);
  else  gr_ax_xs->GetYaxis()->SetLabelSize(0.08);
  gr_ax_xs->GetYaxis()->SetTitleOffset(0.76); 
  gr_ax_xs->GetYaxis()->SetTitleSize(0.08);
  gr_ax_xs->GetXaxis()->SetTitleSize(0.08);
  
  gr_ax_xs->GetXaxis()->SetTitleOffset(1.1);
  if (same) { 
   
    
    gr_ax_xs->GetYaxis()->SetTitleSize(0.09);
    gr_ax_xs->GetYaxis()->SetTitleOffset(0.7); 
  }

  cv_xs->Update();
  if (!doratio) {
    if(do95Lim) gr_errorBand95_xs->Draw("FL");
    gr_errorBand68_xs->Draw("FL");
    if(doBlackHole) {
      gr_theory_BHc->Draw("L");
    }
    else if(doKKgWidth) {
      //cout<<"plotting KKgWidth theory line \n";
      //gr_theory_KKg_widthc->SetLineColor(col_theoryline);
      gr_theory_KKg_widthc2->SetMarkerStyle(28); 
      gr_theory_KKg_widthc2->SetMarkerSize(2);
      gr_theory_KKg_widthc2->SetMarkerColor(gr_theory_KKg_widthc->GetMarkerColor());
      
      
      gr_theory_KKg_widthc2->Draw("P");
      gr_theory_KKg_widthc->Draw("L");
      gr_theory_KKg_widthc->Print();
    }

    else if(doKKGrav) {
      if(doDrawTheoryUncert)  gr_theory_KKGrav->Draw("F");    
      gr_theory_KKGravc->Draw("L");
      cout<<"drawing kkgrav theory line \n";
    }
    else {
      if(doDrawTheoryUncert) gr_theory_DW->Draw("F");
      gr_theory_DWc->Draw("L");   

      if(doZprime3p){


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

  TLegend *leg_limvsmass_xs = new TLegend(0.2,0.57,0.5,0.95);
  if(doObserved ) {
    if (!doratio) gr_observed_xs->Draw("PL");
    leg_limvsmass_xs->AddEntry(gr_observed_xs, lt_obs,"pl");
    if (first) leg_limvsmass_xs->Draw();
  }
  leg_limvsmass_xs->AddEntry(gr_expected_xs, lt_exp,"pl");      
  leg_limvsmass_xs->AddEntry(gr_errorBand68_xs, lt_68,"lf");
  if(do95Lim) leg_limvsmass_xs->AddEntry(gr_errorBand95_xs, lt_95,"lf");
  if (!doratio) {
 if(doDrawTheoryUncert) {
    if(doBlackHole) leg_limvsmass_xs->AddEntry(gr_theory_BHc, lt_th_BH,"lp");
    //else if(doKKg) leg_limvsmass_xs->AddEntry(gr_theory_KKgc, lt_th_KKg,"lp");
    else if(doKKgWidth){
	leg_limvsmass_xs->AddEntry(gr_theory_KKg_widthc, lt_th_KKg,"lfp");
    }
    else if(doKKGrav) leg_limvsmass_xs->AddEntry(gr_theory_KKGrav, lt_th_KKGrav,"lfp");
    else leg_limvsmass_xs->AddEntry(gr_theory_DW, lt_th,"lfp");
 }
  else {
    if(doBlackHole) leg_limvsmass_xs->AddEntry(gr_theory_BHc, lt_th_BH,"l");
    else if(doKKGrav) leg_limvsmass_xs->AddEntry(gr_theory_KKGrav, lt_th_KKGrav,"l");
    else {
      if(doZprime3p){

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
  if (first) leg_limvsmass_xs->Draw();
  
  //Legend for ATLAS/CMS theory points
  TLegend* theoryPointsLeg = new TLegend(0.2,0.70,0.5,1.0);
  theoryPointsLeg->AddEntry(gr_theory_KKg_widthc2, lt_ATLAS_width , "p");
  theoryPointsLeg->AddEntry(gr_theory_KKg_widthc2, lt_CMS_width, "p");
  theoryPointsLeg->SetHeader(legname);
  theoryPointsLeg->SetFillColor(0);
  theoryPointsLeg->SetBorderSize(0);
  theoryPointsLeg->SetTextSize(0.07);
  if(second) theoryPointsLeg->Draw();
  
 
 
  
  //ATLAS labels
 
  TLatex *text22 = new TLatex(32,1.5*y_ax_xs[0],label.Data());
  if (same) text22->SetTextSize(0.09);
  else text22->SetTextSize(0.08);
  text22->Draw();
  
    if (first)
    { 
      
      //      ATLASLabel(0.5, 0.7, cv_xs,true);
    
      // ATLASCMEIntL(0.5,0.6,LUMI,cv_xs);
      TLatex l;
      // l.SetNDC();
      l.SetTextFont(42);
      l.SetTextColor(1);
      l.SetTextSize(0.09);
      char charIntL[200];
      //sprintf(charIntL,"%.2f",IntL);
      sprintf(charIntL,"%.1f",LUMI);
      //l.DrawLatex(x,y,TString("#scale[0.65]{#int} Ldt = "+string(charIntL)+" nb^{-1}      #sqrt{s} = 7 TeV"));
      l.DrawLatex(30,1.,TString("#scale[0.6]{#int}_{  }#font[12]{L_{ }dt} = "+string(charIntL)+" fb^{-1}"));
      //l.DrawLatex(x,y+0.2,"#sqrt{s} = 7 TeV");
      l.DrawLatex(30,0.5,"#sqrt{#font[12]{s}} = 8 TeV");
      cout << " drawing label at NDC position " << endl;
      
    }
 

 gPad->RedrawAxis();
  //Save plots with appropriate names
  TString signal="";
  TString tool="";
  TString stsy="";
  TString sel="";
  TString linlog="";

  if(selection.Contains("Zprime-lim")) signal="_zprime";
  else if(selection.Contains("KKg-lim")) signal="_kkg";
  else if(selection.Contains("KKg-Width-1TeV-lim")) signal="_kkg_width_1tev";
  else if(selection.Contains("KKg-Width-2TeV-lim")) signal="_kkg_width_2tev";
  else if(selection.Contains("KKg-Width-3TeV-lim")) signal="_kkg_width_3tev";
  else if(selection.Contains("KKGrav-lim")) signal="_kkgrav";


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

  //Print limits in a table
  TString lend="", sep=""; 
  if(printTex) {
    lend="\\\\"; sep="&"; 
  }
 
  cout << "The limits (pb)" << endl;
  cout << " Mass " << sep << "  Obs.  " << sep << "  Exp. " << sep 
       << "   -1std " << sep << "  +1std ";
  if(print2sigma) cout << sep << "  -2std " << sep << "  +2std";
  cout<< lend << endl;
  const int massprec=(plotinTeV ? 2: 0);
  for(int i=0; i<Nmasses; i++) {
    /*    cout << setiosflags(ios::fixed) << setprecision(massprec)
	 << setw(5) << signal_mass[i] 
	 << setiosflags(ios::fixed) << setprecision(3)
	 << sep << setw(8) << limitObserved[i] << sep << setw(8) << limitMean[i]
	 << sep << setw(8) << limLow68[i] << sep << setw(8) << limHigh68[i];
    */ 
    cout << setiosflags(ios::fixed) << setprecision(massprec)
	 << setw(5) << signal_mass[i] 
	 << setiosflags(ios::fixed) << setprecision(3)
	 << sep << setw(8) << limitMean[i]
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
