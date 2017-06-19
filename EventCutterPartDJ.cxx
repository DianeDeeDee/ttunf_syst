#include "EventCutterPartDJ.h"
#include <cmath>

EventCutterPartDJ::EventCutterPartDJ() {
}

EventCutterPartDJ::~EventCutterPartDJ() {
}

bool EventCutterPartDJ::select(const Event &e, Event &sel) {
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

  if (jets < 1) return false;

  for (int k = 0; k < e.partLargeJet().size(); ++k) {
    if (e.partLargeJet()[k].pass()) {
      sel.partLargeJet().push_back(e.partLargeJet()[k]);
      fatjets++;
    }
  }
  if (fatjets < 1) return false;

  //Add dijet cut
  if( fatjets > 1 && fabs(sel.partLargeJet()[0].mom().Phi() - sel.partLargeJet()[1].mom().Phi() ) > 2.6 ) return false; 

  int btags = 0;
  for (int k = 0; k < e.partJet().size(); ++k) {
    if (e.partJet()[k].pass() && (e.partJet()[k].trueFlavour() == 5)) {
      btags++;
    }
  }
  if (btags < 1) return false;

  return true;
}

