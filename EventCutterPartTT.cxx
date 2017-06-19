#include "EventCutterPartTT.h"
#include <cmath>

EventCutterPartTT::EventCutterPartTT() {
}

EventCutterPartTT::~EventCutterPartTT() {
}

bool EventCutterPartTT::select(const Event &e, Event &sel) {
  int els = 0;
  int mus = 0;
  int fatjets = 0;
  int jets = 0;

  // preselect jets for OR
  /*
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
      dr = e.partElectron()[k].minDeltaR<Jet>(selectedJets[k], jidx);
      if (dr < 0.2) {
        selectedJets[k].erase(selectedJets[k].begin() + jidx);
        dr = e.partElectron()[k].minDeltaR<Jet>(selectedJets[k], jidx);
      }
      if (dr < 0.4) {
        continue;
      }

      sel.partElectron().push_back(e.partElectron()[k]);
      elidx = k;
      els++;
    }
  }

  //TLorentzVector lepton;
  bool electron = true;
  if (els == 1) {
  //  lepton = sel.partElectron()[0].mom();
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
  }*/

  //if (els == 0 && mus == 1) {
  //  lepton = sel.partMuon()[0].mom();
  //  electron = false;
  //} else if (els == 1 && mus == 0) {
  //  electron = true;
  //} else
  //  return false;

  int elidx = -1;
  std::map<int, int> closeJet;
  std::vector<std::set<int> > jetToEl(e.partJet().size());
  for (int k = 0; k < e.partElectron().size(); ++k) {
    closeJet[k] = -1;
    if (e.partElectron()[k].pass()) {
      int jidx = -1;
      double dr = e.partElectron()[k].minDeltaR<Jet>(e.partJet(), jidx);
      if (jidx < 0) continue;
      closeJet[k] = jidx;
      jetToEl[jidx].insert(k);
    }
  }
  sel.partJet() = e.partJet();
  for (int k = 0; k < sel.partJet().size(); ++k) {
    for (std::set<int>::iterator l = jetToEl[k].begin(); l!= jetToEl[k].end(); ++l) {
      sel.partJet()[k].mom() -= e.partElectron()[*l].mom();
    }
  }
  sel.partElectron().clear();
  for (int k = 0; k < e.partElectron().size(); ++k) {
    if (e.partElectron()[k].pass()) {
      bool good = true;
      for (int l = 0; l < sel.partJet().size(); ++l) {
        if (sel.partJet()[l].mom().Perp() < 25e3)
          continue;
        if (sel.partJet()[l].mom().DeltaR(e.partElectron()[k].mom()) < 0.2) {
          sel.partJet()[closeJet[k]].mom() += e.partElectron()[k].mom();
          good = false;
        }
      }
      if (good) {
        sel.partElectron().push_back(e.partElectron()[k]);
        els++;
      }
    }
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


  sel.partMet(e.partMet().Px(), e.partMet().Py());
  //float mtw = std::sqrt(2*lepton.Perp()*sel.partMet().Perp()*(1 - std::cos(lepton.Phi() - sel.partMet().Phi())) );

  //if (e.partMet().Perp() < 20e3) return false;
  //if (e.partMet().Perp() + mtw < 60e3)  return false;

  // require a jet close to the lepton (dR < 1.5) and select the highest p_T one
  // for future cut
  //int sel_jet = -1;
  //double dr = 1.5;
  //double lpt = 0;
  //for (int k = 0; k < sel.partJet().size(); ++k) {
  //  double tdr = sel.partJet()[k].mom().DeltaR(lepton);
  //  if (tdr < dr && sel.partJet()[k].mom().Perp() > lpt && (sel.partJet()[k].trueflavour() == 5)) {
  //    sel_jet = k;
  //    lpt = sel.partJet()[k].mom().Perp();
  //  }
  //}
  //if (sel_jet == -1) return false;

  for (int k = 0; k < e.partLargeJet().size(); ++k) {
    if (e.partLargeJet()[k].mom().Perp() > 150e3 && std::fabs(e.partLargeJet()[k].mom().Eta()) < 2.0) {
      // now check topological cuts
      sel.partLargeJet().push_back(e.partLargeJet()[k]);
      fatjets++;
    }
  }
  if (fatjets < 1) return false;

  for (int k = 0; k < e.partLargeJetBB().size(); ++k) {
    if (e.partLargeJetBB()[k].mom().Perp() > 150e3 && std::fabs(e.partLargeJetBB()[k].mom().Eta()) < 2.0) {
      // now check topological cuts
      sel.partLargeJetBB().push_back(e.partLargeJetBB()[k]);
    }
  }

  sel.partMom() = e.partMom();

  return true;
}

