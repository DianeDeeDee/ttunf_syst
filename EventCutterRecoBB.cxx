#include "EventCutterRecoBB.h"
#include <cmath>

#include "GoodRunsLists/DQHelperFunctions.h"

EventCutterRecoBB::EventCutterRecoBB() {
}

EventCutterRecoBB::~EventCutterRecoBB() {
}

bool EventCutterRecoBB::select(const Event &e, Event &sel) {
  int jets = 0;



  if (e.isData()) {
    if (!DQ::PassRunLB(e.runNumber(), e.lbn()))
      return false;
  }


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



  //event cleaning
  bool passeventcleaning = true;
  if( (e.lerr() > 1) || ( e.terr() == 2) ||  ((e.cfl()&0x40000) != 0 )) 
    passeventcleaning = false;
  if (!passeventcleaning) return false;



  //Npv
  bool passnpv = true;
  if( e.npv() < 2 )
    passnpv = false;
  if (!passnpv) return false;

 
  sel.largeJetBB().push_back(e.largeJetBB()[0]);

  std::cout << " LargeR jet pT: " <<  e.largeJetBB()[0].mom().Perp()/1000 << " Mass: " <<  e.largeJetBB()[0].mom().M()/1000 <<  " Eta: " <<  e.largeJetBB()[0].mom().Eta() <<  std::endl;
  std::cout << " UG Tau1: " <<  e.largeJetBB()[0].ug_tau1() << " UG Tau2: " <<  e.largeJetBB()[0].ug_tau2() << std::endl;
  std::cout << " BDRS jet pT: " <<  e.largeJetBB()[0].bdrs_pt()/1000 << " Mass: " <<  e.largeJetBB()[0].bdrs_m()/1000 <<  " Eta: " <<  e.largeJetBB()[0].bdrs_eta() <<  std::endl;


  std::cout << " After Selection " << std::endl;

  std::cout << " LargeR jet pT: " <<  sel.largeJetBB()[0].mom().Perp()/1000 << " Mass: " << sel.largeJetBB()[0].mom().M()/1000 <<  " Eta: " <<  sel.largeJetBB()[0].mom().Eta() <<  std::endl;
  std::cout << " UG Tau1: " <<  sel.largeJetBB()[0].ug_tau1() << " UG Tau2: " <<  sel.largeJetBB()[0].ug_tau2() << std::endl;
  std::cout << " BDRS jet pT: " <<  sel.largeJetBB()[0].bdrs_pt()/1000 << " Mass: " <<  sel.largeJetBB()[0].bdrs_m()/1000 <<  " Eta: " <<  sel.largeJetBB()[0].bdrs_eta() <<  std::endl;

  std::cout << " --------------------------------- " << std::endl;
 
  return true;
}

