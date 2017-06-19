#include "TopMuonSFUtils/MuonSF.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include <iostream>

using namespace std;

int main() {
  MuonSF muonsf("data/mc11b/");

  cout << "period2011BtoI = " << TopMuon::period2011BtoI << endl;
  cout << "period2011LtoM = " << TopMuon::period2011LtoM << endl;

  TRandom3 rand(419);

  for(int i=0; i<10; ++i) {
    TLorentzVector testmu1;    
    double pt = -1.0;
    while(pt < 15.0) {
      pt = rand.Gaus(60.0,40.0);
    }
    double eta = -3.0;
    while(eta < -2.5 || eta > 2.5) {
      eta = rand.Gaus(0.0,1.0);
    }    double phi = rand.Uniform(TMath::Pi());
    if(rand.Rndm() > 0.5) phi *= -1.0;    testmu1.SetPtEtaPhiM( pt, eta, phi, 0.106);

    cout << "pt, eta, phi = " << pt << ", " << eta << ", " << phi << endl;
    cout << "\tSF to I = " << muonsf.mu_trigger_SF(testmu1.Eta(), testmu1.Phi(), TopMuon::period2011BtoI ) << endl;
    cout << "\tSFerr to I = " << muonsf.mu_trigger_SF_err(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011BtoI ) << endl;
  }//muon loop

  cout << "Muon ID SF to I = " << muonsf.mu_ID_SF(TopMuon::period2011BtoI) << endl;
  cout << "Muon ID SF J-K  = " << muonsf.mu_ID_SF(TopMuon::period2011JtoK) << endl;
  cout << "Muon ID SF L-M  = " << muonsf.mu_ID_SF(TopMuon::period2011LtoM) << endl;

  return 0;
}

