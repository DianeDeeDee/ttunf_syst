#include "EventCutterRecoTTSig.h"
#include <cmath>

#include "GoodRunsLists/DQHelperFunctions.h"

EventCutterRecoTTSig::EventCutterRecoTTSig() {
}

EventCutterRecoTTSig::~EventCutterRecoTTSig() {
}

bool EventCutterRecoTTSig::select(const Event &e, Event &sel) {
  int els = 0;
  int mus = 0;
  int fatjets = 0;
  int jets = 0;

  sel.cutFlow().clear();
  // CUT 0
  sel.cutFlow().push_back(1.0);


  //if (e.isData()) {
  //  if (!DQ::PassRunLB(e.runNumber(), e.lbn()))
  //    return false;
  //}
  //if (e.terr() == 2) return false;
  //if ((e.cfl() & 0x40000) != 0) return false;
  //if (e.lerr() > 1) return false;

  // CUT 1 - Quality
  sel.cutFlow().push_back(1.0);

  //if (e.npv() < 1) return false;
  // CUT 2 - Vertex
  sel.cutFlow().push_back(1.0);


  // preselect jets for OR
  bool passIsBadLoose = true;
  for (int k = 0; k < e.jet().size(); ++k) {
    if (!e.jet()[k].passBadLooseMinus())
      passIsBadLoose = false;

    if (e.jet()[k].pass()) {
      sel.jet().push_back(e.jet()[k]);
      jets++;
    }
  }
  //if (!passIsBadLoose)
  //  return false;

  // CUT 3 - Jet Quality
  sel.cutFlow().push_back(1.0);

  for (int k = 0; k < e.electron().size(); ++k) {
    if (e.electron()[k].pass()) {
      sel.electron().push_back(e.electron()[k]);
      els++;
    }
  }

  /*
  std::vector<std::vector<Jet> > selectedJets;
  int elidx = -1;
  for (int k = 0; k < e.electron().size(); ++k) {
    selectedJets.push_back(std::vector<Jet>());
    if (e.electron()[k].pass()) {
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
  }*/

  for (int k = 0; k < e.muon().size(); ++k) {
    if (e.muon()[k].pass()) {
      float dr = e.muon()[k].minDeltaR(sel.jet());
      if (dr >= 0.04 + 10e3/e.muon()[k].mom().Perp()) {
        sel.muon().push_back(e.muon()[k]);
        mus++;
      }
    }
  }

  sel.triggerElectron() = e.triggerElectron();
  sel.triggerMuon() = e.triggerMuon();

  sel.met(e.met().Px(), e.met().Py());


  for (int k = 0; k < e.largeJet().size(); ++k) {
    if (e.largeJet()[k].passLoose()) {
      // now check topological cuts
      sel.largeJet().push_back(e.largeJet()[k]);
      fatjets++;
    }
  }
  if (fatjets < 1) return false;

  // copy ljetbb
  for (int k = 0; k < e.largeJetBB().size(); ++k) {
    if (e.largeJetBB()[k].passLoose()) {
      // now check topological cuts
      sel.largeJetBB().push_back(e.largeJetBB()[k]);
    }
  }

  // CUT10 1 fatjet
  sel.cutFlow().push_back(1.0);

  sel.partMom() = e.partMom();

  //sel.hfor() = e.hfor();

  return true;
}

