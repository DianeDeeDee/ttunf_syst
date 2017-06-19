#include "EventCutterReco14.h"
#include <cmath>

#include "GoodRunsLists/DQHelperFunctions.h"

EventCutterReco14::EventCutterReco14() {
}

EventCutterReco14::~EventCutterReco14() {
}

bool EventCutterReco14::select(const Event &e, Event &sel) {
  int els = 0;
  int mus = 0;
  int fatjets = 0;
  int jets = 0;

  for (int k = 0; k < e.jet().size(); ++k) {
    if (e.jet()[k].pass()) {
      sel.jet().push_back(e.jet()[k]);
      jets++;
    }
  }

  std::vector<std::vector<Jet> > selectedJets;
  int elidx = -1;
  for (int k = 0; k < e.electron().size(); ++k) {
    selectedJets.push_back(std::vector<Jet>());
    if (e.electron()[k].passLoose()) {
      selectedJets[k] = sel.jet();
      float dr = 99;
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
    if (e.muon()[k].passLoose()) {
      float dr = e.muon()[k].minDeltaR(sel.jet());
      if (dr >= 0.04 + 10e3/e.muon()[k].mom().Perp()) {
        sel.muon().push_back(e.muon()[k]);
        mus++;
      }
    }
  }


  if (els == 0 && mus == 1) {
    lepton = sel.muon()[0].mom();
    electron = false;
  } else if (els == 1 && mus == 0) {
    lepton = sel.electron()[0].mom();
    electron = true;
  } else
    return false;

  sel.triggerElectron() = e.triggerElectron();
  sel.triggerMuon() = e.triggerMuon();

  sel.met(e.met().Px(), e.met().Py());
  float mtw = std::sqrt(2*lepton.Perp()*sel.met().Perp()*(1 - std::cos(lepton.Phi() - sel.met().Phi())) );
  //if (e.met().Perp() < 30e3 || e.met().Perp() + mtw < 60e3)  return false;

  // require a jet close to the lepton (dR < 1.5) and select the highest p_T one
  // for future cut
  int sel_jet = -1;
  double dr = 1.5;
  double lpt = 0;
  for (int k = 0; k < sel.jet().size(); ++k) {
    double tdr = sel.jet()[k].mom().DeltaR(lepton);
    if (tdr < dr && sel.jet()[k].mom().Perp() > lpt) {
      sel_jet = k;
      lpt = sel.jet()[k].mom().Perp();
    }
  }
  if (sel_jet == -1) return false;


  for (int k = 0; k < e.largeJet().size(); ++k) {
    if (e.largeJet()[k].passLoose()) {
      // now check topological cuts
      bool passdR = true;
      bool passdPhi = true;
      bool passM = e.largeJet()[k].mom().M() > 100e3;
      bool passD12 = true; //e.largeJet()[k].split12() > 40e3;
      if (passdR && passdPhi && passM && passD12) {
        sel.largeJet().push_back(e.largeJet()[k]);
        fatjets++;
      }
    }
  }
  if (fatjets < 1) return false;



  int btags = 0;
  for (int k = 0; k < sel.jet().size(); ++k) {
    if (sel.jet()[k].pass() && sel.jet()[k].trueFlavour() == 5) {
      btags++;
    }
  }
  if (btags < 1) return false;


  return true;
}

