#include <iostream>
#include <string>
#include <vector>
#include <math.h>

#include "TFile.h"
#include "TH1F.h"

//Returns the total muhat uncertainty minus the uncertainty related to the nuisance parameters referred to from the "toSubtract" vector (see the config/*xml file)
void getFiducialCrossSectionMuhatUncertainty(std::string folder) {
  std::string total="total";
  std::string statistical="statistical";
  std::vector<std::string> toSubtract;
  toSubtract.push_back("theoSigInclPlusAccept");
  int nToSubtract = toSubtract.size();

  //Histograms retrieved from the different files
  TH1 * hists[nToSubtract+2];

  TFile * fTotal = new TFile(std::string("root-files/Breakdown_" + folder + "_breakdown_add/" + total + ".root").c_str());
  TFile * fStat = new TFile(std::string("root-files/Breakdown_" + folder + "_breakdown_add/" + statistical + ".root").c_str());
  hists[0] = (TH1*)fTotal->Get(total.c_str());
  hists[1] = (TH1*)fStat->Get(statistical.c_str());

  TFile * fToSubtract[nToSubtract];
  for (int i=0; i<nToSubtract; i++) {
    fToSubtract[i] = new TFile(std::string("root-files/Breakdown_" + folder + "_breakdown_add/" + toSubtract[i] + ".root").c_str());
    hists[i+2] = (TH1*)fToSubtract[i]->Get(toSubtract[i].c_str());
  }

  //Build an output hist
  TH1 * hOutput = (TH1*)hists[0]->Clone("hOutput");
  
  //Now compute the difference in quadrature
  for (int ibin=1; ibin<hists[0]->GetNbinsX()+1; ibin++) {
    float newTot = pow(hists[0]->GetBinContent(ibin),2);
    for (int ih=0; ih<nToSubtract; ih++) 
      newTot -= pow(hists[ih+2]->GetBinContent(ibin),2) - pow(hists[1]->GetBinContent(ibin),2);
    hOutput->SetBinContent(ibin, sqrt(newTot));
  }

  //Dumping info for bin #2 (+ unc. on muhat) and #3 (- unc. on muhat)
  std::cout << "Uncertainty on muhat (excluding ";
  for (int i=0; i<nToSubtract; i++) 
    std::cout << ((i) ? ", " : "") << toSubtract[i].c_str();
  std::cout << ") : ";
  std::cout << "+" << hOutput->GetBinContent(2);
  std::cout << " / -" << hOutput->GetBinContent(3) << std::endl;
  
}
