#include "EventCutterRecoTTBkgData.h"
#include <cmath>

#include "GoodRunsLists/DQHelperFunctions.h"

#include "SkimReader.h"

EventCutterRecoTTBkgData::EventCutterRecoTTBkgData() {
}

EventCutterRecoTTBkgData::~EventCutterRecoTTBkgData() {
}

bool EventCutterRecoTTBkgData::select(const Event &e, Event &sel) {
  int els = 0;
  int mus = 0;
  int fatjets = 0;
  int jets = 0;

  float cf_w = 1.0;

  sel.cutFlow().clear();
  for (int z = 0; z < 30; ++z) sel.cutFlow().push_back(0.0);

  // CUT 0
  sel.cutFlow()[0] += cf_w;

  // CUT 1 - GRL
  /*
  const LargeJet *leading = 0;
  for (int k = 0; k < e.largeJet().size(); ++k) {
    if (e.largeJet()[k].passLoose()) {
      leading = &(e.largeJet()[k]);
      break;
    }
  }

  if (e.isData()) {
    if (leading) {
      if (leading->mom().Perp() > 250e3 && leading->mom().Perp() < 400e3) {
	DQ::SetXMLFile("extra/periodI_grl_250GeV_400GeV.xml");
      } else if (leading->mom().Perp() >= 400e3) {
	DQ::SetXMLFile("extra/periodI_grl_gr_400GeV.xml");
      }

      if (!DQ::PassRunLB(e.runNumber(), e.lbn()))
	return false;
    }
  }
  */
  if (e.isData()) {
    if (!DQ::PassRunLB(e.runNumber(), e.lbn()))
      return false;
  }
  sel.cutFlow()[1] += cf_w;

  // CUT 2 - Quality
  if (e.terr() == 2) return false;
  if ((e.cfl() & 0x40000) != 0) return false;
  if (e.lerr() > 1) return false;

  sel.cutFlow()[2] += cf_w;

  // CUT 3 - Vertex
  if (e.npv_good() < 1) return false;
  sel.cutFlow()[3] += cf_w;

  // CUT 4 - trigger
  if (e.isData() && !e.triggerElectron())
    return false;
  sel.cutFlow()[4] += cf_w;

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
  if (!passIsBadLoose)
    return false;

  // CUT 5 - Jet Quality
  sel.cutFlow()[5] += cf_w;

  sel.triggerElectron() = e.triggerElectron();
  sel.triggerMuon() = e.triggerMuon();
  sel.triggerLargeJet() = e.triggerLargeJet();

  sel.met(e.met().Px(), e.met().Py());

  sel.electron().clear();
  for (int k = 0; k < e.electron().size(); ++k) {
    if (e.electron()[k].passLoose()) {
      sel.electron().push_back(e.electron()[k]);
      els++;
    }
  }
  sel.muon().clear();
  for (int k = 0; k < e.muon().size(); ++k) {
    if (e.muon()[k].passLoose()) {
      mus++;
    }
  }

  // veto events that pass the trigger and have offline electrons
  if (e.isData() && (els >= 1))
    return false;

  sel.cutFlow()[6] += cf_w;

  /*
  // CUT 6 - trigger match
  // check trigger
  TLorentzVector triggerJet;
  int iTriggerJet = -1;
  float drTriggerJet = 999.0;
  int iJetMatched = -1;
  for (size_t iTrig = 0; iTrig < e.sr()->trig_EF_jet_n; iTrig++) {
    // check trigger bit
    if (!e.sr()->trig_EF_jet_EF_j360_a4tchad->at(iTrig)) continue;
    // get fourvector
    triggerJet.SetPtEtaPhiE(e.sr()->trig_EF_jet_pt->at(iTrig),
			    e.sr()->trig_EF_jet_eta->at(iTrig),
			    e.sr()->trig_EF_jet_phi->at(iTrig),
			    e.sr()->trig_EF_jet_E->at(iTrig));
    
    // now loop over jets
    for (size_t iJet = 0; iJet < sel.jet().size(); iJet++) {
      if (std::fabs(sel.jet().at(iJet).mom().DeltaR(triggerJet)) < drTriggerJet) {
	iTriggerJet = iTrig;
	drTriggerJet = std::fabs(sel.jet().at(iJet).mom().DeltaR(triggerJet));
        iJetMatched = iJet;
      }
    }
  }
  if (iTriggerJet != -1 && drTriggerJet < 0.4) {
    triggerJet.SetPtEtaPhiE(e.sr()->trig_EF_jet_pt->at(iTriggerJet),
			    e.sr()->trig_EF_jet_eta->at(iTriggerJet),
			    e.sr()->trig_EF_jet_phi->at(iTriggerJet),
			    e.sr()->trig_EF_jet_E->at(iTriggerJet));
  } else {
    return false;
  }
  sel.triggerJet4mom() = triggerJet;
  */

  // CUT 5
  // cut on large-R jets
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
  sel.cutFlow()[7] += cf_w;

  /*
  if (sel.largeJet()[0].mom().DeltaR(sel.jet()[iJetMatched].mom()) < 0.4) return false;
  sel.cutFlow().push_back(1.0);

  if (sel.largeJet()[0].mom().DeltaPhi(sel.jet()[iJetMatched].mom()) < 0.3) return false;
  sel.cutFlow().push_back(1.0);

  if (sel.largeJet()[0].mom().DeltaPhi(sel.jet()[iJetMatched].mom()) > TMath::Pi() - 0.3) return false;
  sel.cutFlow().push_back(1.0);
  */

  sel.partMom() = e.partMom();
  //sel.hfor() = e.hfor();

  // save the trigger electron
  TLorentzVector triggerEl;
  float tPt = 0;
  int iTriggerEl = -1;
  for (size_t iTrig = 0; iTrig < e.sr()->trig_EF_el_n; iTrig++) {
    // check trigger bit
    if (!e.sr()->trig_EF_el_EF_e24vhi_medium1->at(iTrig) && !e.sr()->trig_EF_el_EF_e60_medium1->at(iTrig)) continue;
    // get fourvector
    triggerEl.SetPtEtaPhiE(e.sr()->trig_EF_el_pt->at(iTrig),
                           e.sr()->trig_EF_el_eta->at(iTrig),
                           e.sr()->trig_EF_el_phi->at(iTrig),
                           e.sr()->trig_EF_el_E->at(iTrig));
    if (triggerEl.Perp() > tPt) {
      iTriggerEl = iTrig;
      tPt = triggerEl.Perp();
    }
  }
  if (iTriggerEl != -1) {
    triggerEl.SetPtEtaPhiE(e.sr()->trig_EF_el_pt->at(iTriggerEl),
                           e.sr()->trig_EF_el_eta->at(iTriggerEl),
                           e.sr()->trig_EF_el_phi->at(iTriggerEl),
                           e.sr()->trig_EF_el_E->at(iTriggerEl));
  }
  sel.triggerJet4mom() = triggerEl;

  return true;
}

