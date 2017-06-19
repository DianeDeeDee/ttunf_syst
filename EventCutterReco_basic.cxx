#include "EventCutterReco.h"
#include <cmath>
#include <iostream>
using namespace std;
#include "packages/GoodRunsLists/GoodRunsLists/DQHelperFunctions.h"
//Used when we preselect!e
EventCutterReco::EventCutterReco() {
}

EventCutterReco::~EventCutterReco() {
  }
bool EventCutterReco::select(const Event &e, Event &sel) {
int els = 0;
  int mus = 0;
  int fatjets = 0;
  int jets = 0;
  int btags = 0;
  sel.cutFlow().clear();
  // CUT 0
  sel.cutFlow()[0] += 1;
 if (e.isData()) {
    if (!DQ::PassRunLB(e.runNumber(), e.lbn())){return false;}
    sel.cutFlow()[1] += 1;
  }
  if (e.terr() == 2) return false;
  if ((e.cfl() & 0x40000) != 0) return false;
  if (e.lerr() > 1) return false;

  // CUT 1 - Quality
  sel.cutFlow()[2] += 1;
  if(sel.npv() < 1) return false;
  // CUT 2 - Vertex
  sel.cutFlow()[3] += 1;

 // preselect jets for OR
  bool passIsBadLoose = true; //this flag is dor Data only
for (int k = 0; k < e.largeJet().size(); ++k) {
}

  for (int k = 0; k < e.jet().size(); ++k) {
    if (!e.jet()[k].passBadLoose()) passIsBadLoose = false;
    if (e.jet()[k].pass()) {
      sel.jet().push_back(e.jet()[k]);
      jets++;
    }
  }
  sel.cutFlow()[4] += 1;
 if (!passIsBadLoose) {
 return false;
}
  // CUT 3 - Jet  Quality
 sel.cutFlow()[5] += 1;
 std::vector<std::vector<Jet> > selectedJets;
  int elidx = -1;

  for (int k = 0; k < e.electron().size(); ++k) {
    selectedJets.push_back(std::vector<Jet>());
    if (e.electron()[k].pass()) {
      selectedJets[k] = sel.jet();
      float dr = 99.;
      int jidx = -1;
      dr = e.electron()[k].minDeltaR<Jet>(selectedJets[k], jidx);
      if (dr < 0.2) {
        selectedJets[k].erase(selectedJets[k].begin() + jidx);
        dr = e.electron()[k].minDeltaR<Jet>(selectedJets[k], jidx);
      }
      if (dr < 0.4) {
        continue;
      }
      sel.electron().push_back(e.electron()[k]);
      elidx = k;
      els++;
    }
  }
  
  TLorentzVector lepton;
  bool electron = true;
  if (e.triggerElectron() && els == 1) {
    lepton = sel.electron()[0].mom();
    electron = true;
    sel.jet() = selectedJets[elidx];
  }

  for (int k = 0; k < e.muon().size(); ++k) {
      float dr = e.muon()[k].minDeltaR(sel.jet());
      if (dr >= 0.04 + 10e3/e.muon()[k].mom().Perp()) {
	sel.muon().push_back(e.muon()[k]);
        mus++;
      }
    }

  // CUT 4 - Trigger
 // if (e.triggerMuon() || e.triggerElectron())   sel.cutFlow().push_back(1.0);
if (e.triggerMuon() && els == 0 && mus == 1) {
    lepton = sel.muon()[0].mom();
    electron = false;
    sel.cutFlow()[6]+=1;
  } else if (e.triggerElectron() && els == 1 && mus == 0) {
    electron = true;
    sel.cutFlow()[7]+=1;
  } else
    return false;
  sel.cutFlow()[8]+=1;
  // CUT 5 - == 1 lepton
  sel.triggerElectron() = e.triggerElectron();
  sel.triggerMuon() = e.triggerMuon();

  sel.cutFlow()[9]+=1;

//all in MeV
  sel.met(e.met().Px(), e.met().Py());
  float mtw = std::sqrt(2*lepton.Perp()*sel.met().Perp()*(1 - std::cos(lepton.Phi() - sel.met().Phi())) );
 // CUT 6 - MET . 20 GeV
  if (e.met().Perp() > 20e3) sel.cutFlow()[10] += 1 ;

//if (e.met().Perp() < 30e3 || e.met().Perp() + mtw < 60e3)  return false;
if (e.met().Perp() < 20e3 || e.met().Perp() + mtw < 60e3)  return false;

  // CUT7 MET+MTW > 60 GeV
  sel.cutFlow()[11] += 1;

  // require a jet close to the lepton (dR < 1.5) and select the highest p_T one
  // for future cut
  double dr = 1.5;   
  int sel_jet = -1;
  
  double lpt = 0;
  for (int k = 0; k < sel.jet().size(); ++k) {
    double tdr = sel.jet()[k].mom().DeltaR(lepton);
    if (tdr < dr && sel.jet()[k].mom().Perp() > lpt) {
      sel_jet = k;
      lpt = sel.jet()[k].mom().Perp();
    }
  }
  if (sel_jet == -1) return false;

 // CUT9 dR(l,j) < 1.5
  sel.cutFlow()[12] += 1;
  
for (int k = 0; k < e.largeJet().size(); ++k) {

    //if (e.largeJet()[k].pass()) 
    if (e.largeJet()[k].passLoose()) {
      bool passdR = sel.jet()[sel_jet].mom().DeltaR(e.largeJet()[k].mom()) > 1.5; //1.5
      //the 2nd SM top is boosted nut not enough to be a fat jet==> lepton and b are well seperated
      bool passdPhi = true;
      //e.largeJet()[k].mom().DeltaPhi(lepton) > 2.3; //2.3
      bool passM = e.largeJet()[k].mom().M() > 60e3; //100e3
      bool passD12 = e.largeJet()[k].split12() > 20e3;
      if (passdR && passdPhi && passM && passD12) {
        sel.largeJet().push_back(e.largeJet()[k]);
        fatjets++;
      }
    }
  }
  
//  if (fatjets < 1) return false;
if (jets <=3  && fatjets < 1) return false;
 
  // CUT10 1 fatjet
  sel.cutFlow()[13] += 1;

  for (int k = 0; k < sel.jet().size(); ++k) {
    if (sel.jet()[k].btag()) {
      btags++;
    }
  }
  sel.cutFlow()[14] += 1;
  if (btags < 1) {
  	return false;
}
 // CUT 8 : 1 b-tag
 sel.cutFlow()[15] += 1;
  return true;

}
