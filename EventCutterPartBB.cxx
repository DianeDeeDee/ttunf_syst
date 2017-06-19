#include "EventCutterPartBB.h"
#include <cmath>

EventCutterPartBB::EventCutterPartBB() {
}

EventCutterPartBB::~EventCutterPartBB() {
}

bool EventCutterPartBB::select(const Event &e, Event &sel) {

  sel.partLargeJetBB().push_back(e.partLargeJetBB()[0]);

  if( e.largeJetBB()[0].mom().DeltaR(e.partLargeJetBB()[0].mom()) < (0.75*1.2) ) return false;

  return true;
}

