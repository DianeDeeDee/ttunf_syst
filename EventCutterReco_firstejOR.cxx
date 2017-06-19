#include "EventCutterReco.h"
#include <cmath>
#include "SkimReader.h"

#include "GoodRunsLists/DQHelperFunctions.h"

EventCutterReco::EventCutterReco(bool loose)
  :m_loose(loose) {
}

EventCutterReco::~EventCutterReco() {
}

bool EventCutterReco:: ef_mu24Hypo(float eta, float pt) {
  if(eta < 1.05) {
    if(pt > 23.34) return true;
  } else if(eta < 1.5) {
    if(pt > 23.19) return true;
  } else if(eta < 2.0) {
    if(pt > 23.14) return true;
  } else {
    if(pt > 23.06) return true;
  }
  return false;
}

bool EventCutterReco:: ef_mu36Hypo(float eta, float pt) {
  if(eta < 1.05) {
    if(pt > 34.96) return true;
  } else if(eta < 1.5) {
    if(pt > 34.78) return true;
  } else if(eta < 2.0) {
    if(pt > 34.69) return true;
  } else {
    if(pt > 34.63) return true;
  }
  return false;
}

bool EventCutterReco::select(const Event &e, Event &sel) {
  int els = 0;
  int mus = 0;
  int fatjets = 0;
  int jets = 0;

  const int cf_mu = 15;
  float cf_w = 1.0;

  sel.cutFlow().clear();
  for (int z = 0; z < 30; ++z) sel.cutFlow().push_back(0.0);

  // CUT 0
  sel.cutFlow()[0] += cf_w;
  sel.cutFlow()[cf_mu] += cf_w;

  // CUT 1 - Trigger
  if (e.terr() == 2) return false;
  if ((e.cfl() & 0x40000) != 0) return false;
  if (e.lerr() > 1) return false;

  if (e.triggerElectron()) sel.cutFlow()[1] += cf_w;
  if (e.triggerMuon()) sel.cutFlow()[cf_mu+1] += cf_w;

  if (!e.triggerElectron() && !e.triggerMuon()) return false;

  if (e.isData()) {
    if (!DQ::PassRunLB(e.runNumber(), e.lbn()))
      return false;
  }

  // CUT 2 - GRL
  if (e.triggerElectron()) sel.cutFlow()[2] += cf_w;
  if (e.triggerMuon()) sel.cutFlow()[cf_mu+2] += cf_w;

  if (e.npv_good() < 1) return false;

  // CUT 3 - Vertex
  if (e.triggerElectron()) sel.cutFlow()[3] += cf_w;
  if (e.triggerMuon()) sel.cutFlow()[cf_mu+3] += cf_w;


  // preselect jets for OR
  bool passIsBadLoose = true;
  sel.jet().clear();
  for (int k = 0; k < e.jet().size(); ++k) {
    if (!e.jet()[k].passBadLooseMinus()){
     //if (!e.jet()[k].passBadLoose())
      passIsBadLoose = false;}

    if (e.jet()[k].pass()) {
      sel.jet().push_back(e.jet()[k]);
      jets++;
    }
  }
/* if (!passIsBadLoose) {
                        return false;
                } 
*/
  int els_loose = 0;
  int els_tight = 0;
  sel.electron().clear();
  for (int k = 0; k < e.electron().size(); ++k) {
    if (e.electron()[k].pass() && !m_loose) {
      sel.electron().push_back(e.electron()[k]);
      els++;
    }
    // only loose for QCD here:
    if (e.electron()[k].passLoose() && m_loose) {
      sel.electron().push_back(e.electron()[k]);
      els++;
      els_loose++;
    }
    if (e.electron()[k].pass() && m_loose) els_tight++;
  }

  int mus_loose = 0;
  int mus_tight = 0;
  sel.muon().clear();
  for (int k = 0; k < e.muon().size(); ++k) {
    if (e.muon()[k].pass() && !m_loose) {
      float dr = e.muon()[k].minDeltaR(sel.jet());
      if (dr >= 0.04 + 10e3/e.muon()[k].mom().Perp()) {
        sel.muon().push_back(e.muon()[k]);
        mus++;
      }
    }
    // only loose selection for QCD here
    if (e.muon()[k].passLoose() && m_loose) {
      float dr = e.muon()[k].minDeltaR(sel.jet());
      if (dr >= 0.04 + 10e3/e.muon()[k].mom().Perp()) {
        sel.muon().push_back(e.muon()[k]);
        mus++;
        mus_loose++;
      }
    }
    if (e.muon()[k].pass() && m_loose) {
      float dr = e.muon()[k].minDeltaR(sel.jet());
      if (dr >= 0.04 + 10e3/e.muon()[k].mom().Perp()) {
        mus_tight++;
      }
    }
  }

  // CUT 4 - >= 1 lepton
  if (!m_loose) {
    if (e.triggerElectron() && els >= 1) sel.cutFlow()[4] += cf_w;
    if (e.triggerMuon() && mus >= 1) sel.cutFlow()[cf_mu+4] += cf_w;
  } else { // loose
    if (e.triggerElectron() && els_loose >= 1) sel.cutFlow()[4] += cf_w;
    if (e.triggerMuon() && mus_loose >= 1) sel.cutFlow()[cf_mu+4] += cf_w;
  }

  // CUT 5 - == 1 lepton
  if (!m_loose) {
    if (e.triggerElectron() && els == 1) sel.cutFlow()[5] += cf_w;
    if (e.triggerMuon() && mus == 1) sel.cutFlow()[cf_mu+5] += cf_w;
  } else {
    if (e.triggerElectron() && els_loose == 1) sel.cutFlow()[5] += cf_w;
    if (e.triggerMuon() && mus_loose == 1) sel.cutFlow()[cf_mu+5] += cf_w;
  }

  // CUT 6 - == 0 other lepton
  if (!m_loose) {
    if (e.triggerElectron() && els == 1 && mus == 0) sel.cutFlow()[6] += cf_w;
    if (e.triggerMuon() && mus == 1 && els == 0) sel.cutFlow()[cf_mu+6] += cf_w;
  } else {
    if (e.triggerElectron() && els_loose == 1 && mus_tight == 0) {
      els = 1;
      mus = 0;
      sel.cutFlow()[6] += cf_w;
    }
    if (e.triggerMuon() && mus_loose == 1 && els_tight == 0){
      els = 0;
      mus = 1;
      sel.cutFlow()[cf_mu+6] += cf_w;
    }
  }

  TLorentzVector lepton;
  bool electron = true;
  if (e.triggerMuon() && els == 0 && mus == 1) {
    lepton = sel.muon()[0].mom();
    electron = false;
  } else if (e.triggerElectron() && els == 1 && mus == 0) {
    lepton = sel.electron()[0].mom();
    electron = true;
    //sel.jet() = selectedJets[elidx];
  } else
    return false;
/// 947 entries/3912

 /// CUT 7 - trigger match (TODO)
  if (electron) {
    TLorentzVector triggerEl;
    float mindr = 999;
    int iTriggerEl = -1;
    for (size_t iTrig = 0; iTrig < e.sr()->trig_EF_el_n; iTrig++) {
      // check trigger bit
      if (!e.sr()->trig_EF_el_EF_e24vhi_medium1->at(iTrig) && !e.sr()->trig_EF_el_EF_e60_medium1->at(iTrig)) continue;
      // get fourvector
      triggerEl.SetPtEtaPhiE(e.sr()->trig_EF_el_pt->at(iTrig),
                             e.sr()->trig_EF_el_eta->at(iTrig),
                             e.sr()->trig_EF_el_phi->at(iTrig),
                             e.sr()->trig_EF_el_E->at(iTrig));
      if (triggerEl.DeltaR(lepton) < mindr) {
        iTriggerEl = iTrig;
        mindr = triggerEl.DeltaR(lepton);
      }
    }
    if (iTriggerEl != -1) {
      triggerEl.SetPtEtaPhiE(e.sr()->trig_EF_el_pt->at(iTriggerEl),
                             e.sr()->trig_EF_el_eta->at(iTriggerEl),
                             e.sr()->trig_EF_el_phi->at(iTriggerEl),
                             e.sr()->trig_EF_el_E->at(iTriggerEl));
    }
    sel.triggerJet4mom() = triggerEl;
    if (mindr >= 0.15)
      return false;
  } else {
    TLorentzVector triggerMu;
    float mindr = 999;
    int iTriggerMu = -1;
    int iTrkMu = -1;
    for (size_t iTrig = 0; iTrig < e.sr()->trig_EF_trigmuonef_n; iTrig++) {
      // check trigger bit
      if (!e.sr()->trig_EF_trigmuonef_EF_mu24i_tight->at(iTrig) && !e.sr()->trig_EF_trigmuonef_EF_mu36_tight->at(iTrig)) continue;

      for (size_t iTrk = 0; iTrk < e.sr()->trig_EF_trigmuonef_track_n->at(iTrig); iTrk++) {
        // get fourvector
        bool passes = false;
        if (e.sr()->trig_EF_trigmuonef_EF_mu24i_tight->at(iTrig) && ef_mu24Hypo(e.sr()->trig_EF_trigmuonef_track_CB_eta->at(iTrig).at(iTrk), e.sr()->trig_EF_trigmuonef_track_CB_pt->at(iTrig).at(iTrk)*1e-3)) passes = true;
        if (e.sr()->trig_EF_trigmuonef_EF_mu36_tight->at(iTrig) && ef_mu36Hypo(e.sr()->trig_EF_trigmuonef_track_CB_eta->at(iTrig).at(iTrk), e.sr()->trig_EF_trigmuonef_track_CB_pt->at(iTrig).at(iTrk)*1e-3)) passes = true;
        if (!passes) continue;

        triggerMu.SetPtEtaPhiM(e.sr()->trig_EF_trigmuonef_track_CB_pt->at(iTrig).at(iTrk),
                               e.sr()->trig_EF_trigmuonef_track_CB_eta->at(iTrig).at(iTrk),
                               e.sr()->trig_EF_trigmuonef_track_CB_phi->at(iTrig).at(iTrk),
                               0);
        if (triggerMu.DeltaR(lepton) < mindr) {
          iTriggerMu = iTrig;
          iTrkMu = iTrk;
          mindr = triggerMu.DeltaR(lepton);
        }
      }
    }
    if (iTriggerMu != -1) {
      triggerMu.SetPtEtaPhiM(e.sr()->trig_EF_trigmuonef_track_CB_pt->at(iTriggerMu).at(iTrkMu),
                             e.sr()->trig_EF_trigmuonef_track_CB_eta->at(iTriggerMu).at(iTrkMu),
                             e.sr()->trig_EF_trigmuonef_track_CB_phi->at(iTriggerMu).at(iTrkMu),
                             0);
    }
    sel.triggerJet4mom() = triggerMu;
    if (mindr >= 0.15)
      return false;
  }
  if (electron) sel.cutFlow()[7] += cf_w;
  else sel.cutFlow()[cf_mu+7] += cf_w;
  
  // CUT 8 - jet cleaning
  if (!passIsBadLoose)   return false;

  if (electron) sel.cutFlow()[8] += cf_w;
  else sel.cutFlow()[cf_mu+8] += cf_w;

  sel.triggerElectron() = e.triggerElectron();
  sel.triggerMuon() = e.triggerMuon();

  sel.met(e.met().Px(), e.met().Py());
  float mtw = std::sqrt(2*lepton.Perp()*sel.met().Perp()*(1 - std::cos(lepton.DeltaPhi(sel.met()))) );

  // CUT 9 - MET > 20 GeV
  if (e.met().Perp() < 20e3) return false;
  if (electron) sel.cutFlow()[9] += cf_w;
  else sel.cutFlow()[cf_mu+9] += cf_w;

  // CUT 10 MET + MTW > 60 GeV
  if (e.met().Perp() + mtw < 60e3)  return false;
  if (electron) sel.cutFlow()[10] += cf_w;
  else sel.cutFlow()[cf_mu+10] += cf_w;
//New: electron OR for better QCD
if (electron){
int sel_jetOR = -1;
  double drOR = 0.4;
  double mydR=1.5;
  for (int k = 0; k < sel.jet().size(); ++k) {
    double tdrOR = sel.jet()[k].mom().DeltaR(lepton);

   if ( tdrOR < mydR) {
      sel_jetOR = k; //mark also the electron....
      //lpt = sel.jet()[k].mom().Perp();
           mydR=tdrOR;
            }
        }
    if(mydR<0.4)  return false;
  }                    
  // require a jet close to the lepton (dR < 1.5) and select the highest p_T one
  // for future cut
  int sel_jet = -1;
  double dr = 1.5;
  double lpt = 0;
  for (int k = 0; k < sel.jet().size(); ++k) {
    double tdr = sel.jet()[k].mom().DeltaR(lepton);
    if (tdr < dr && sel.jet()[k].mom().Perp() > lpt && sel.jet()[k].btag()) {
      sel_jet = k;
      lpt = sel.jet()[k].mom().Perp();
    }
  }
  if (sel_jet == -1) return false;

  // CUT 11 -- jet close to lepton (differs from tt res due to the btag)
  if (electron) sel.cutFlow()[11] += cf_w;
  else sel.cutFlow()[cf_mu+11] += cf_w;

  sel.largeJet().clear();
  for (int k = 0; k < e.largeJet().size(); ++k) {
    if (e.largeJet()[k].passLoose()) {
       bool passdR = true;//sel.jet()[sel_jet].mom().DeltaR(e.largeJet()[k].mom()) > 1.5;
      bool passdPhi =  true;//e.largeJet()[k].mom().DeltaPhi(lepton) > 2.3; //2.3
      bool passM = e.largeJet()[k].mom().M() > 100e3;
      bool passD12 = e.largeJet()[k].split12() > 20e3; //40e3 pour ttbar paper
      if (passdR && passdPhi && passM && passD12) {
      // now check topological cuts
      sel.largeJet().push_back(e.largeJet()[k]);
      fatjets++;
      }
    }
  }
 //if (fatjets < 1) return false;
 if (jets <= 3 && fatjets < 1) return false;
  // CUT 12 -- Akt10 p_T > 200 GeV, |eta| < 2.0 (differs from tt res due to the missing topo cuts and lower p_T cuts)
  sel.cutFlow()[12] += cf_w;
  sel.cutFlow()[cf_mu+12] += cf_w;

  // copy ljetbb
/*  sel.largeJetBB().clear();
  for (int k = 0; k < e.largeJetBB().size(); ++k) {
    if (e.largeJetBB()[k].passLoose()) {
      // now check topological cuts
      sel.largeJetBB().push_back(e.largeJetBB()[k]);
    }
  }

  sel.partMom() = e.partMom();
  sel.hfor() = e.hfor();
*/
  return true;
}

