#include "EventCutterPart14.h"
#include <cmath>

EventCutterPart14::EventCutterPart14() {
}

EventCutterPart14::~EventCutterPart14() {
}

bool EventCutterPart14::select(const Event &e, Event &sel) {
  int els = 0;
  int mus = 0;
  int fatjets = 0;
  int jets = 0;

  for (int k = 0; k < e.partJet().size(); ++k) {
    if (e.partJet()[k].pass()) {
      sel.partJet().push_back(e.partJet()[k]);
      jets++;
    }
  }

  std::vector<std::vector<Jet> > selectedJets;
  int elidx = -1;
  for (int k = 0; k < e.partElectron().size(); ++k) {
    selectedJets.push_back(std::vector<Jet>());
    if (e.partElectron()[k].pass()) {
      selectedJets[k] = sel.partJet();
      float dr = 99;
      int jidx = -1;
      dr = e.partElectron()[k].minDeltaR(selectedJets[k], jidx);
      if (dr < 0.2) {
        selectedJets[k].erase(selectedJets[k].begin() + jidx);
        dr = e.partElectron()[k].minDeltaR(selectedJets[k], jidx);
      }
      if (dr < 0.4) {
        continue;
      }
      sel.partElectron().push_back(e.partElectron()[k]);
      elidx = k;
      els++;
    }
  }

  TLorentzVector lepton;
  bool electron = true;
  if (els == 1) {
    lepton = sel.partElectron()[0].mom();
    electron = true;
    sel.partJet() = selectedJets[elidx];
  }

  for (int k = 0; k < e.partMuon().size(); ++k) {
    if (e.partMuon()[k].pass()) {
      float dr = e.partMuon()[k].minDeltaR(sel.partJet());
      if (dr >= 0.04 + 10e3/e.partMuon()[k].mom().Perp()) {
        sel.partMuon().push_back(e.partMuon()[k]);
        mus++;
      }
    }
  }
  if (els == 0 && mus == 1) {
    lepton = sel.partMuon()[0].mom();
    electron = false;
  } else if (els == 1 && mus == 0) {
    lepton = sel.partElectron()[0].mom();
    electron = true;
  } else
    return false;


  sel.partMet(e.partMet().Px(), e.partMet().Py());
  float mtw = std::sqrt(2*lepton.Perp()*sel.partMet().Perp()*(1 - std::cos(lepton.Phi() - sel.partMet().Phi())) );
  //if (e.partMet().Perp() < 30e3 || e.partMet().Perp() + mtw < 60e3)  return false;


  if (jets < 1) return false;

  // require a jet close to the lepton (dR < 1.5) and select the highest p_T one
  // for future cut
  int sel_jet = -1;
  double dr = 1.5;
  double lpt = 0;
  for (int k = 0; k < sel.partJet().size(); ++k) {
    double tdr = sel.partJet()[k].mom().DeltaR(lepton);
    if (tdr < dr && sel.partJet()[k].mom().Perp() > lpt) {
      sel_jet = k;
      lpt = sel.partJet()[k].mom().Perp();
    }
  }
  if (sel_jet == -1) return false;


  for (int k = 0; k < e.partLargeJet().size(); ++k) {
    if (e.partLargeJet()[k].passLoose()) {
      bool passdR = true;
      bool passdPhi = true;
      bool passM = e.partLargeJet()[k].mom().M() > 70e3;
      bool passD12 = true;//e.largeJet()[k].split12() > 40e3;
      if (passdR && passdPhi && passM && passD12) {
        sel.partLargeJet().push_back(e.partLargeJet()[k]);
        fatjets++;
      }
    }
  }
  if (fatjets < 1) return false;

  int btags = 0;
  for (int k = 0; k < e.partJet().size(); ++k) {
    if (e.partJet()[k].pass() && (e.partJet()[k].trueFlavour() == 5)) {
      btags++;
    }
  }
  if (btags < 1) return false;

  sel.partMom().clear();
  for (int k = 0; k < e.partMom().size(); ++k) {
    sel.partMom().push_back(e.partMom()[k]);
  }

  return true;
}

