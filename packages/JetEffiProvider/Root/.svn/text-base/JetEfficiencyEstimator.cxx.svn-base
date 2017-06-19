/*
EXAMPLE of how to use:
    root -l 
    gROOT->LoadMacro("JetEfficiencyEstimator.cxx+");
    JetEfficiencyEstimator* JEE = new JetEfficiencyEstimator();
    for a given jet pT and eta, do:
    bool status = JEE->isGoodJet(pt, eta);
    status is false if the jet is thrown away and true if kept by the rejection method
*/
#include "JetEffiProvider/JetEfficiencyEstimator.h"
#include <iostream>
#include <cmath>

//_________________________________________________________________________________________________
JetEfficiencyEstimator::JetEfficiencyEstimator(unsigned int seed): 
  m_rndmGen(seed) { // default seed
  // these are the efficiencies before averaging over new bins using the tool as provided by Tancredi Carli on  February 1st 2013
  // this map was used to derive smoother efficiencies.
  /*
  std::vector<double> h_in_x, h_in_y, h_in_e;
  double xmin = ( 5.5+5.778)*1.377;
  double xmax = (49.5+5.778)*1.377;
  for (std::map<double,double>::const_iterator it = effiMap.begin(); it != effiMap.end(); ++it){
    double x = ((*it).first +5.778)*1.377;
    double y = (*it).second;
    double e = 0.5 *( effiMapHigh[(*it).first] - effiMapLow[(*it).first]);
    h_in_x.push_back(x);
    h_in_y.push_back(y);
    h_in_e.push_back(TMath::Abs(e));
  }
  TH1D* h_in = new TH1D("h_in","",65,15,80); h_in->Sumw2();
  for (int i=0; i< h_in_x.size(); i++){
    double x = h_in_x[i];
    double y = h_in_y[i];
    double e = h_in_e[i];
    int b = h_in->FindBin(x);
    h_in->SetBinContent(b,y);
    h_in->SetBinError(b,e);
  }
   std::vector<Double_t> bin;
   bin.push_back( 15);
   bin.push_back( 20);
   bin.push_back( 30);
   bin.push_back( 45);
   bin.push_back( 65);
   bin.push_back( 80);
   TH1D* h_out = binAverage(h_in,bin);
   */
  std::cout << "JetEfficiencyEstimator : Instanciating a JetEfficiencyEstimator constructeur with seed " << seed << std::endl;
  std::cout << "JetEfficiencyEstimator : The input jet pT is given in *GeV*" << std::endl;

  m_effiMap[   5.5] = 0.935711; m_effiMapLow[   5.5] = 0.850000; m_effiMapHigh[   5.5] = 1.04160;
  m_effiMap[   6.5] = 0.970928; m_effiMapLow[   6.5] = 0.919684; m_effiMapHigh[   6.5] = 1.02217;
  m_effiMap[   7.5] = 0.981771; m_effiMapLow[   7.5] = 0.948085; m_effiMapHigh[   7.5] = 1.01546;
  m_effiMap[   8.5] = 0.967075; m_effiMapLow[   8.5] = 0.910947; m_effiMapHigh[   8.5] = 1.02320;
  m_effiMap[   9.5] = 0.997874; m_effiMapLow[   9.5] = 0.988423; m_effiMapHigh[   9.5] = 1.00733;
  m_effiMap[  10.5] = 0.999399; m_effiMapLow[  10.5] = 0.992375; m_effiMapHigh[  10.5] = 1.00642;
  m_effiMap[  11.5] = 0.988956; m_effiMapLow[  11.5] = 0.969349; m_effiMapHigh[  11.5] = 1.00856;
  m_effiMap[  12.5] = 0.996535; m_effiMapLow[  12.5] = 0.989615; m_effiMapHigh[  12.5] = 1.00346;
  m_effiMap[  13.5] = 0.997884; m_effiMapLow[  13.5] = 0.992921; m_effiMapHigh[  13.5] = 1.00285;
  m_effiMap[  14.5] = 0.997568; m_effiMapLow[  14.5] = 0.992598; m_effiMapHigh[  14.5] = 1.00254;
  m_effiMap[  15.5] = 0.997912; m_effiMapLow[  15.5] = 0.993604; m_effiMapHigh[  15.5] = 1.00222;
  m_effiMap[  16.5] = 0.996423; m_effiMapLow[  16.5] = 0.989768; m_effiMapHigh[  16.5] = 1.00308;
  m_effiMap[  17.5] = 1.000000; m_effiMapLow[  17.5] = 0.998961; m_effiMapHigh[  17.5] = 1.00104;
  m_effiMap[  18.5] = 0.998934; m_effiMapLow[  18.5] = 0.995993; m_effiMapHigh[  18.5] = 1.00188;
  m_effiMap[  19.5] = 1.000000; m_effiMapLow[  19.5] = 0.998374; m_effiMapHigh[  19.5] = 1.00163;
  m_effiMap[  20.5] = 1.000000; m_effiMapLow[  20.5] = 0.997964; m_effiMapHigh[  20.5] = 1.00204;
  m_effiMap[  21.5] = 1.000000; m_effiMapLow[  21.5] = 0.997716; m_effiMapHigh[  21.5] = 1.00228;
  m_effiMap[  22.5] = 1.000000; m_effiMapLow[  22.5] = 0.997254; m_effiMapHigh[  22.5] = 1.00275;
  m_effiMap[  23.5] = 1.000000; m_effiMapLow[  23.5] = 0.996558; m_effiMapHigh[  23.5] = 1.00344;
  m_effiMap[  24.5] = 0.996154; m_effiMapLow[  24.5] = 0.985610; m_effiMapHigh[  24.5] = 1.00670;
  m_effiMap[  25.5] = 1.000000; m_effiMapLow[  25.5] = 0.994786; m_effiMapHigh[  25.5] = 1.00521;
  m_effiMap[  26.5] = 1.000000; m_effiMapLow[  26.5] = 0.994462; m_effiMapHigh[  26.5] = 1.00554;
  m_effiMap[  27.5] = 1.000000; m_effiMapLow[  27.5] = 0.992183; m_effiMapHigh[  27.5] = 1.00782;
  m_effiMap[  28.5] = 1.000000; m_effiMapLow[  28.5] = 0.992104; m_effiMapHigh[  28.5] = 1.00790;
  m_effiMap[  29.5] = 1.000000; m_effiMapLow[  29.5] = 0.990066; m_effiMapHigh[  29.5] = 1.00993;
  m_effiMap[  30.5] = 1.000000; m_effiMapLow[  30.5] = 0.989332; m_effiMapHigh[  30.5] = 1.01067;
  m_effiMap[  31.5] = 1.000000; m_effiMapLow[  31.5] = 0.989303; m_effiMapHigh[  31.5] = 1.01070;
  m_effiMap[  32.5] = 1.000000; m_effiMapLow[  32.5] = 0.985875; m_effiMapHigh[  32.5] = 1.01413;
  m_effiMap[  33.5] = 1.000000; m_effiMapLow[  33.5] = 0.985657; m_effiMapHigh[  33.5] = 1.01434;
  m_effiMap[  34.5] = 1.000000; m_effiMapLow[  34.5] = 0.978832; m_effiMapHigh[  34.5] = 1.02117;
  m_effiMap[  35.5] = 1.000000; m_effiMapLow[  35.5] = 0.978097; m_effiMapHigh[  35.5] = 1.02190;
  m_effiMap[  36.5] = 1.000000; m_effiMapLow[  36.5] = 0.966820; m_effiMapHigh[  36.5] = 1.03318;
  m_effiMap[  37.5] = 1.000000; m_effiMapLow[  37.5] = 0.970104; m_effiMapHigh[  37.5] = 1.02990;
  m_effiMap[  38.5] = 1.000000; m_effiMapLow[  38.5] = 0.982668; m_effiMapHigh[  38.5] = 1.01733;
  m_effiMap[  39.5] = 1.000000; m_effiMapLow[  39.5] = 0.963181; m_effiMapHigh[  39.5] = 1.03682;
  m_effiMap[  40.5] = 1.000000; m_effiMapLow[  40.5] = 0.968612; m_effiMapHigh[  40.5] = 1.03139;
  m_effiMap[  41.5] = 1.000000; m_effiMapLow[  41.5] = 0.958191; m_effiMapHigh[  41.5] = 1.04181;
  m_effiMap[  42.5] = 1.000000; m_effiMapLow[  42.5] = 0.954372; m_effiMapHigh[  42.5] = 1.04563;
  m_effiMap[  43.5] = 1.000000; m_effiMapLow[  43.5] = 0.931296; m_effiMapHigh[  43.5] = 1.06870;
  m_effiMap[  44.5] = 1.000000; m_effiMapLow[  44.5] = 0.943704; m_effiMapHigh[  44.5] = 1.05630;
  m_effiMap[  45.5] = 1.000000; m_effiMapLow[  45.5] = 0.944162; m_effiMapHigh[  45.5] = 1.05584;
  m_effiMap[  46.5] = 1.000000; m_effiMapLow[  46.5] = 0.943704; m_effiMapHigh[  46.5] = 1.05630;
  m_effiMap[  47.5] = 1.000000; m_effiMapLow[  47.5] = 0.898648; m_effiMapHigh[  47.5] = 1.10135;
  m_effiMap[  48.5] = 1.000000; m_effiMapLow[  48.5] = 0.908510; m_effiMapHigh[  48.5] = 1.09149;
  m_effiMap[  49.5] = 1.000000; m_effiMapLow[  49.5] = 0.946459; m_effiMapHigh[  49.5] = 1.05354;
  m_useICHEPResult = false; // default
}

//_________________________________________________________________________________________________
bool JetEfficiencyEstimator::isGoodJet(double pt, double eta){
  double eff = 1; 
  if (std::fabs(eta) > 2.1) return true; // beyond 2.5-0.4, no way to check.
  if (pt < 15) pt =15; // below 5GeV threshold, keep threshold.
  if (m_useICHEPResult==true){
    eff = 0.98; // 2% for ICHEP flat in pT
  }
  else{
    /*
    // find first the jet pt bin:
    double kt = pt/1.377 - 5.778;
    for (std::map<double, double>::iterator itr = m_effiMap.begin(); itr != m_effiMap.end(); ++itr){
      if (itr->first > kt) {
	eff = m_effiMap[itr->first];
 	break; //
      }
    } // end of loop over bins
    */
    if      (pt >=30         ) return true; // fully efficient
    else if (pt >=15 && pt<20) eff = 0.97359; // 2.7% 
    else if (pt >=20 && pt<30) eff = 0.99772; // 0.23%
    else                       eff = 1;
  } // end of case ICHEP 
  double r = m_rndmGen.Uniform();
  if (r < std::fabs(1-eff)) return false;
  return true;
}
//_________________________________________________________________________________________________
bool JetEfficiencyEstimator::checkEvent(std::vector<double> pt, std::vector<double> eta, int& N){
  std::vector<bool> keep;
  return this->checkEvent(pt,eta,N,keep);
}

//_________________________________________________________________________________________________
bool JetEfficiencyEstimator::checkEvent(std::vector<double> pt, std::vector<double> eta, int& N, std::vector<bool>& keep){
  // check for inconsistent input first
  if (pt.size() != eta.size()) return false;
  N = pt.size(); // default value
  keep.resize(N); for (int i=0; i<N; i++) keep[i] = true;
  for (unsigned int i=0; i< pt.size(); i++){
    if (this->isGoodJet(pt[i], eta[i]) == false){
      N--; // this jet should be discarded
      keep[i] = false;
    }
  } // end of loop over all jets
  return true;
}
