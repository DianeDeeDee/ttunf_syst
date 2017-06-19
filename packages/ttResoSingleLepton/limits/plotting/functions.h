#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include "TString.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TPaveText.h"

void getMass(TString selection, double *mass, const int Nmasses, const bool massinTeV=false);
void getMass(TString selection, double *mass, const int Nmasses, const bool massinTeV) {

  double dummy[9]={500,600,700,800,1000,1300,1600};
  if(selection=="syst4j" || selection=="systdr" || selection=="newlum4j" || selection=="newlumdr") {
    dummy[0]=500;
    dummy[1]=600;
    dummy[2]=700;
    dummy[3]=800;
    dummy[4]=1000;
  }
  if(selection.Contains("r164j") ||selection.Contains("r16dr") ) {
    dummy[0]=500;
    dummy[1]=700;
    dummy[2]=1000;
  }
  if(selection.Contains("Zprime-lims")){
    dummy[0]=500;
    dummy[1]=600;
    dummy[2]=700;
    dummy[3]=800;
    dummy[4]=1000;
  }
  if(selection.Contains("Zprime-lim_pretag")) {
    dummy[0]=500;
    dummy[1]=600;
    dummy[2]=700;
    dummy[3]=800;
    dummy[4]=1000;
    dummy[5]=1300;
    dummy[6]=1600;
    dummy[7]=2000;
    dummy[8]=3000;
  }
  if(selection.Contains("KKg-lim_pretag")){
    dummy[0]=700;
    dummy[1]=800;
    dummy[2]=1000;
    dummy[3]=1300;
    dummy[4]=1600;
    //dummy[5]=1800;
    dummy[5]=2000;
  }
  if(selection.Contains("KKGrav-lim_pretag")){
    dummy[0]=500;
    dummy[1]=600;
    dummy[2]=700;
    dummy[3]=800;
    dummy[4]=900;
    dummy[5]=1000;
    dummy[6]=1300;
  }
  if(selection.Contains("BH-lim")){
    //dummy[0]=750;
    //dummy[1]=1500;
    //dummy[2]=2250;
    dummy[0]=750;
    dummy[1]=1500;
    dummy[2]=2000;
    dummy[3]=2250;
    dummy[4]=2500;
    //dummy[5]=2750;
  }
  if(selection.Contains("qbh")){
    //dummy[0]=750;
    //dummy[1]=1500;
    //dummy[2]=2250;
    dummy[0]=750;
    dummy[1]=1500;
    //dummy[2]=2000;
    dummy[2]=2250;
    dummy[3]=2500;
    dummy[4]=2750;
  }
  if(selection.Contains("_sc")){
    dummy[0]=500;
    dummy[1]=600;
    dummy[2]=700;
    dummy[3]=800;
    dummy[4]=1000;
    dummy[5]=1300;
    dummy[6]=1600;
  }
  for(int i=0;i<Nmasses;++i) {
    if(massinTeV) mass[i]=dummy[i]/1000;
    else mass[i]=dummy[i];
  }

}


void getXsec(TString selection, double *xsec, const int Nmasses);
void getXsec(TString selection, double *xsec, const int Nmasses) {

  //double dummy[5]={500,700,1000,1600,2000};
  //double dummy[5]={3.552, 1.373, 0.310, 0.027, 0.0063};
  double dummy[5];
  dummy[0] = 3.552; dummy[1] = 2.246;
  dummy[2] = 1.373; dummy[3] = 0.826;
  dummy[4] = 0.310;

  if(selection=="10tev") {
    dummy[0]=1;
    dummy[1]=1;
    dummy[2]=1;
  }
  if(selection=="venkat") {
    dummy[0]=3.552;
    dummy[1]=1.373;
    dummy[2]=0.310;
  }
  if(selection=="clermond") {
    dummy[0] = 35.52; 
    dummy[1] = 137.3; 
    dummy[2] = 31.0; 
  }
  if(selection=="spano") {
    dummy[0] = 1.373; 
    dummy[1] = 0.310; 
    dummy[2] = 0.027; 
  }
  if(selection=="lucia4jet" || selection=="luciadrcut") {
    dummy[0] = 3.552; 
    dummy[1] = 1.373; 
    dummy[2] = 0.310; 
  }
  if(selection=="wuppertal") {
    dummy[0] = 3.552; 
    dummy[1] = 1.373; 
    dummy[2] = 0.310; 
    dummy[3] = 0.027; 
    dummy[4] = 0.0063; 
  }
  if(selection.Contains("syst4j") || selection.Contains("systdr") || selection.Contains("newlum4j") || selection.Contains("newlumdr")) {
    dummy[0] = 3.552; 
    dummy[1] = 2.246;
    dummy[2] = 1.373;
    dummy[3] = 0.826;
    dummy[4] = 0.310;
  }
  if(selection.Contains("r164j") || selection.Contains("r16dr")) {
    dummy[0] = 3.552; 
    dummy[1] = 1.373;
    dummy[2] = 0.310;
  }

  for(int i=0;i<Nmasses;++i) {
    xsec[i]=dummy[i];
  }


}

TString getName(TString selection);
TString getName(TString selection) {

  TString Name="";
  if(selection=="10tev") Name="10 TeV MC";  
  else if(selection.Contains("-lims")) Name="CLs, dRmin";
  else if(selection.Contains("-lim_dr")) Name="dRmin";
  else if(selection.Contains("-lim_4j")) Name="4 jets";
  else if(selection.Contains("-lim_nc")) Name="#chi^{2} no cal";
  else if(selection.Contains("-lim_chi2cal")) Name="#chi^{2} cal";  
  else if(selection.Contains("-lim_pretag")) Name="Boosted"; 
  else if(selection.Contains("_sc4j")) Name="GCBT, 4 jets";
  else if(selection.Contains("_scdr")) Name="GCBT, dRmin";
  else if(selection.Contains("_scnc")) Name="GCBT, #chi^{2} no cal";
  else if(selection.Contains("_scx2")) Name="GCBT, #chi^{2} rescale";
  else if(selection.Contains("_kkgsc4j")) Name="GCBT, 4 jets";
  else if(selection.Contains("_kkgscdr")) Name="GCBT, dRmin";
  else if(selection.Contains("_kkgscnc")) Name="GCBT, #chi^{2} no cal";
  else if(selection.Contains("_kkgscx2")) Name="GCBT, #chi^{2} rescale";


  return Name;
}

TGraph * getErrorBand(const int N, double* x, double* yUp, double* yDown);
TGraph * getErrorBand(const int N, double* x, double* yUp, double* yDown) {
  // Returns a TGraph with the error band.
  //const int N=x.size;//(x);
  //const int Ne = (N>1 ? 2*N+1 : 2*N);
  assert(N>0);

  const int Ne = 2*N+1;

  double xe[Ne], ye[Ne];
  if(N>1) {
    for(int i=0;i<N;i++) {
      xe[i] = x[i];
      xe[2*N-1-i] = x[i];
      ye[i] = yUp[i];
      ye[2*N-1-i] = yDown[i];
    }
    xe[Ne-1] = xe[0];
    ye[Ne-1] = yUp[0];
  }
  else {
    xe[0] = x[0];     
    xe[1] = x[0];
    xe[2] = x[0];
    ye[0] = yUp[0];
    ye[1] = yDown[0];
    ye[2] = yUp[0];
  }

  TGraph *errorBand = new TGraph(Ne,xe,ye);

  //cout << "Point  x   y" << endl;
  //for(int i=0;i<Ne;++i) {
  //  cout << i << "  " << xe[i] << "  " << ye[i] << endl;
  //}

  return errorBand;

}

TGraph * getErrorBand(const int N, double* x, double* yCentral, double* yPlus, double* yMinus);
TGraph * getErrorBand(const int N, double* x, double* yCentral, double* yPlus, double* yMinus) {
  // Returns a TGraph with the error band.
  //yCental is the central value, with plus and minus limits
  //const int N=x.size;//(x);
  //const int Ne = (N>1 ? 2*N+1 : 2*N);
  assert(N>0);

  const int Ne = 2*N+1;

  double xe[Ne], ye[Ne];
  if(N>1) {
    for(int i=0;i<N;i++) {
      xe[i] = x[i];
      xe[2*N-1-i] = x[i];
      ye[i] = yCentral[i] + yPlus[i];
      ye[2*N-1-i] = yCentral[i] - fabs(yMinus[i]);
    }
    xe[Ne-1] = xe[0];
    ye[Ne-1] = yCentral[0] + yPlus[0];
  }
  else {
    xe[0] = x[0];     
    xe[1] = x[0];
    xe[2] = x[0];
    ye[0] = yCentral[0] + yPlus[0];
    ye[1] = yCentral[0] - fabs(yMinus[0]);
    ye[2] = yCentral[0] + yPlus[0];
  }

  TGraph *errorBand = new TGraph(Ne,xe,ye);

  //cout << "Point  x   y" << endl;
  //for(int i=0;i<Ne;++i) {
  //  cout << i << "  " << xe[i] << "  " << ye[i] << endl;
  //}

  return errorBand;

}

TString getSystName(TString systematics);
TString getSystName(TString systematics) {

  TString Name=systematics;
  if(systematics=="") Name ="Nominal"; 
  if(systematics=="jer") Name ="JER up"; 
  if(systematics=="jesu") Name ="JES up"; 
  if(systematics=="jesd") Name ="JES down"; 
  if(systematics=="bjesu") Name ="BJES up";
  if(systematics=="bjesd") Name ="BJES down";
  if(systematics=="jrec") Name ="Jet reco";
  if(systematics=="btagd") Name ="Btag down"; 
  if(systematics=="btagu") Name ="Btag up"; 
  if(systematics=="midptd") Name ="MuID pt down"; 
  if(systematics=="midptu") Name ="MuID pt up"; 
  if(systematics=="mmsptd") Name ="MuMS pt down"; 
  if(systematics=="mmsptu") Name ="MuMS pt up";
  if(systematics=="musfd") Name ="MuSF down";
  if(systematics=="musfu") Name ="MuSF up";
  if(systematics=="vvxsd") Name ="VV xs down";
  if(systematics=="vvxsu") Name ="VV xs up";
  if(systematics=="stxsd") Name ="Sing top xs down";
  if(systematics=="stxsu") Name ="Sing top xs up";
  if(systematics=="ttxsd") Name ="ttbar xs down";
  if(systematics=="ttxsu") Name ="ttbar xs up";
  if(systematics=="lumd") Name ="Lumi down";
  if(systematics=="lumu") Name ="Lumi up";
  if(systematics=="elsfd") Name ="ElSF down";
  if(systematics=="elsfu") Name ="ElSF up";
  if(systematics=="elptd") Name ="El pt down";
  if(systematics=="elptu") Name ="El pt up";
  if(systematics=="elresu") Name ="El res up";
  if(systematics=="fsrd") Name ="FSR down";
  if(systematics=="fsru") Name ="FSR up";
  if(systematics=="isrd") Name ="ISR down";
  if(systematics=="isru") Name ="ISR up";
  if(systematics=="ifsrd") Name ="I/FSR down";
  if(systematics=="ifsru") Name ="I/FSR up";
  if(systematics=="wjsfd") Name ="Wjet SF down";
  if(systematics=="wjsfu") Name ="Wjet SF up";
  if(systematics=="wjhf") Name ="Wjet HF";
  if(systematics=="wjallu") Name ="Wjet all up";
  if(systematics=="wjalld") Name ="Wjet all down";
  if(systematics=="qcde") Name ="QCD anti-ele";
  if(systematics=="crosscheck") Name ="Nom. (Mainz)";
  /*if(systematics=="lumd1") Name ="Lumi down 1%";
  if(systematics=="lumd2") Name ="Lumi down 2%";
  if(systematics=="lumd3") Name ="Lumi down 3%";
  if(systematics=="lumd4") Name ="Lumi down 4%";
  if(systematics=="lumd5") Name ="Lumi down 5%";
  if(systematics=="lumd6") Name ="Lumi down 6%";
  if(systematics=="lumd7") Name ="Lumi down 7%";
  if(systematics=="lumd8") Name ="Lumi down 8%";
  if(systematics=="lumd9") Name ="Lumi down 9%";
  if(systematics=="lumd10") Name ="Lumi down 10%";
  if(systematics=="lumd11") Name ="Lumi down 11%";
  */
  if(systematics=="lumu2") Name ="Lumi up 2.2%";
  if(systematics=="lumd2") Name ="Lumi down 2.2%";
  if(systematics=="lumu6") Name ="Lumi up 6.6%";
  if(systematics=="lumd6") Name ="Lumi down 6.6%";
  if(systematics=="lumu10") Name ="Lumi up 11%";
  if(systematics=="lumd10") Name ="Lumi down 11%";
  if(systematics=="lumu12") Name ="Lumi up 13.2%";
  if(systematics=="lumd12") Name ="Lumi down 13.2%";
  if(systematics=="lumu16") Name ="Lumi up 17.6%";
  if(systematics=="lumd16") Name ="Lumi down 17.6%";
  if(systematics=="lumu20") Name ="Lumi up 22%";
  if(systematics=="lumd20") Name ="Lumi down 22%";
  if(systematics=="lumu22") Name ="Lumi up 24.2%";
  if(systematics=="lumd22") Name ="Lumi down 24.2%";
  if(systematics=="lumu26") Name ="Lumi up 28.6%";
  if(systematics=="lumd26") Name ="Lumi down 28.6%";
  if(systematics=="lumu30") Name ="Lumi up 33%";
  if(systematics=="lumd30") Name ="Lumi down 33%";


  return Name;

}

void ATLAS_LABEL(Double_t x,Double_t y,Color_t color,TCanvas* cv) 
{
  if (cv==NULL) return;
  cv->cd();
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);
  l.DrawLatex(x,y,"ATLAS");
  return;
}

void ATLASPRELIM_LABEL(Double_t x,Double_t y,Color_t color,TCanvas* cv) 
{
  if (cv==NULL) return;
  cv->cd();
  TLatex* l = new TLatex(x,y,""); //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l->SetNDC();
  l->SetTextFont(72);
  l->SetTextColor(color);
  l->DrawLatex(x,y,"ATLAS");
  TLatex* l2 = new TLatex(x,y,""); //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l2->SetNDC();
  l2->SetTextColor(color);
  l2->DrawLatex(x+0.175,y,"Preliminary");
  return;
}

void ATLASLabel(Double_t x,Double_t y,TCanvas* cv,bool Preliminary=false,Color_t color=1) {
  if (cv==NULL) return;
  cv->cd();
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);

  double delx = 0.115*696*gPad->GetWh()/(472*gPad->GetWw());
  double dely = 0.029*696*gPad->GetWh()/(472*gPad->GetWw());
  TPaveText *bk = new TPaveText(x,y-0.01,x+delx,y+dely,"NDC");
  //cout << "gPad->GetWh()=" << gPad->GetWh() << endl;
  //cout << "gPad->GetWw()=" << gPad->GetWw() << endl;
  //cout << 696.0*gPad->GetWh()/(472.0*gPad->GetWw()) << endl;
  //cout << "Coordinates: " << x << ", " << y << endl;
  //cout << "           : " << x+delx << ", " << y+dely << endl;
  //bk->AddText("ATLAS");
  bk->SetTextColor(kWhite);
  bk->SetFillColor(kWhite);
  bk->SetBorderSize(0);
  bk->Draw();
  l.DrawLatex(x,y,"ATLAS");
  //l.SetTextFont(42);  l.SetTextSize(0.04);  l.DrawLatex(x,y-0.08,"#sqrt{#font[12]{s}} = 7 TeV");  l.SetTextFont(72);  l.SetTextSize(0.05);
  if (Preliminary) {
    TLatex p; 
    p.SetNDC();
    p.SetTextFont(42);
    p.SetTextColor(color);
    double delx2 = 0.115*696*gPad->GetWh()/(472*gPad->GetWw());
    TPaveText *bk2 = new TPaveText(x+delx,y-0.01,x+delx+delx2,y+dely,"NDC");
    bk2->SetTextColor(kWhite);
    bk2->SetFillColor(kWhite);
    bk2->SetBorderSize(0);
    bk2->Draw();
    //p.DrawLatex(x+delx,y,"Preliminary"); //Preliminary
    //p.DrawLatex(x+delx,y,"for approval"); 
    p.DrawLatex(x+delx,y,"Internal"); 
    //    p.DrawLatex(x,y,"#sqrt{s}=900GeV");
  }
  return;
}

void ATLASCMEIntL(Double_t x,Double_t y,Double_t IntL,TCanvas* cv,Color_t color=1) {
  if (cv==NULL) return;
  cv->cd();
  TLatex l;
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextColor(color);
  l.SetTextSize(0.04);
  char charIntL[200];
  //sprintf(charIntL,"%.2f",IntL);
  sprintf(charIntL,"%.1f",IntL);
  //l.DrawLatex(x,y,TString("#scale[0.65]{#int} Ldt = "+string(charIntL)+" nb^{-1}      #sqrt{s} = 7 TeV"));
  l.DrawLatex(x,y,TString("#scale[0.6]{#int}_{  }#font[12]{L_{ }dt} = "+string(charIntL)+" fb^{-1}"));
  //l.DrawLatex(x,y+0.2,"#sqrt{s} = 7 TeV");
  l.DrawLatex(x,y+0.08,"#sqrt{#font[12]{s}} = 8 TeV");
  return;
}

void ATLASCME(Double_t x,Double_t y,TCanvas* cv,Color_t color=1) {
  if (cv==NULL) return;
  cv->cd();
  TLatex l;
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextColor(color);
  l.SetTextSize(0.04);
  l.DrawLatex(x,y,"#sqrt{s} = 8 TeV");
  return;
}



