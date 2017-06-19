#include "EventCutterRecoTTBkg.h"
#include <cmath>

#include "GoodRunsLists/DQHelperFunctions.h"

EventCutterRecoTTBkg::EventCutterRecoTTBkg() {
}

EventCutterRecoTTBkg::~EventCutterRecoTTBkg() {
}

bool EventCutterRecoTTBkg::select(const Event &e, Event &sel) {
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

  // copy ljetbb
  for (int k = 0; k < e.largeJetBB().size(); ++k) {
    if (e.largeJetBB()[k].passLoose()) {
      // now check topological cuts
      sel.largeJetBB().push_back(e.largeJetBB()[k]);
    }
  }

  if (fatjets < 1) return false;

  // CUT10 1 fatjet
  sel.cutFlow().push_back(1.0);

  sel.partMom() = e.partMom();
  //sel.hfor() = e.hfor();

  return true;
}

