#include "EventCutterRecoDJ.h"
#include <cmath>

EventCutterRecoDJ::EventCutterRecoDJ() {
}

EventCutterRecoDJ::~EventCutterRecoDJ() {
}

bool EventCutterRecoDJ::select(const Event &e, Event &sel) {
  int fatjets = 0;
  int jets = 0;

 

  // preselect jets for OR
  bool passIsBadLoose = true;
  for (int k = 0; k < e.jet().size(); ++k) {
    if (!e.jet()[k].passBadLoose())
      passIsBadLoose = false;

    if (e.jet()[k].pass()) {
      sel.jet().push_back(e.jet()[k]);
      jets++;
    }
  }
  if (!passIsBadLoose)
    return false;




  // Fatjet Trigger
  /*if(e.triggerLargeJet()){

   for (int k = 0; k < e.largeJet().size(); ++k) {
     if (e.largeJet()[k].pass()) {
       sel.largeJet().push_back(e.largeJet()[k]);
       fatjets++;
     }
   }
  } //FJ Trigger
*/
  if (fatjets < 1) return false;

  //Add dijet cut
  if( fatjets > 1 && fabs(sel.largeJet()[0].mom().Phi() - sel.largeJet()[1].mom().Phi() ) > 2.6 ) return false; 

 // sel.triggerLargeJet() = e.triggerLargeJet();

 
  return true;
}

