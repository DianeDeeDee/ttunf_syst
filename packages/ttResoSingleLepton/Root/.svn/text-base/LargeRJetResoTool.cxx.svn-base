#include <string>
#include <iostream>
#include "TFile.h"
#include "TMath.h"
#include "ttResoSingleLepton/LargeRJetResoTool.h"

LargeRJetResoTool::LargeRJetResoTool(const string &fname){
  inputFile=TFile::Open(fname.data());
  JEREta1 = (TH1F *) inputFile->Get("_histo_jer_eta1");
  JEREta2 = (TH1F *) inputFile->Get("_histo_jer_eta2");
  JEREta3 = (TH1F *) inputFile->Get("_histo_jer_eta3");
  JMREta1 = (TH1F *) inputFile->Get("_histo_jmr_eta1");
  JMREta2 = (TH1F *) inputFile->Get("_histo_jmr_eta2");
  JMREta3 = (TH1F *) inputFile->Get("_histo_jmr_eta3");
}

LargeRJetResoTool::~LargeRJetResoTool(){
  inputFile->Close();
}



float LargeRJetResoTool::GetResolution(float pTGeV, float eta, resType iType){
  float abseta = TMath::Abs(eta);

  TH1F* hist;
  if (iType==JER) {
    if ( abseta<0.8) hist=JEREta1;
    else if (abseta<1.2) hist=JEREta2;
    else if (abseta<2.0) hist=JEREta3;
    else {
      cout << Form("Warning: No large-R jet resolution available for eta=%.2f ", eta) << endl;    
      return -999.;
    }
  }
  else if (iType==JMR) {
    if ( abseta<0.8) hist=JMREta1;
    else if (abseta<1.2) hist=JMREta2;
    else if (abseta<2.0) hist=JMREta3;
    else {
      cout << Form("Warning: No large-R jet resolution available for eta=%.2f ", eta) << endl;    
      return -999.;
    }
  }
  else {
    cout<<"ERROR: Unknown iType"<<endl;
    return -999.;
  }

  if(pTGeV<200)  return hist->GetBinContent(1)/100.;
  else if(pTGeV<300)  return hist->GetBinContent(2)/100.;
  else if(pTGeV<400)  return hist->GetBinContent(3)/100.;
  else if(pTGeV<500)  return hist->GetBinContent(4)/100.;
  else if(pTGeV<600)  return hist->GetBinContent(5)/100.;
  else if(pTGeV<700)  return hist->GetBinContent(6)/100.;
  else if(pTGeV<900)  return hist->GetBinContent(7)/100.;
  else if(pTGeV<1800) return hist->GetBinContent(8)/100.;
  else {
    cout << Form("Warning: No large-R jet resolution available for pt=%.1f GeV ", pTGeV) << endl;        
    return hist->GetBinContent(8)/100.;
  }
}


float LargeRJetResoTool::GetSmearFactor(float pTGeV, float eta, resType iType, UInt_t seed) {
  SetSeed(seed);
  float  res=GetResolution(pTGeV, eta, iType);
  double smear=rdn.Gaus(0, 1);
  return 0.66*res*smear;  // Smear to achieve 20% larger resolution (current recommendation of resolution uncertainties)
}


