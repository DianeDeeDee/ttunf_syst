#include "EventCutterReco.h"
#include "SkimReader.h"
#include "KinematicUtils.h"
#include "KinematicUtils.cxx"
using namespace std;
#include "packages/GoodRunsLists/GoodRunsLists/DQHelperFunctions.h"
//#include "TMath.h"
#include <stdio.h>
//#include <stdlib.h> // for "abs"
#include <math.h> // for "fabs"
//#include <cmath> // for "std::abs" and "std::fabs"

#include <iostream>

//Used when we preselect!e
//EventCutterReco::EventCutterReco() {
//}
EventCutterReco::EventCutterReco(bool loose)
  :m_loose(loose) {
}
EventCutterReco::~EventCutterReco() {
}

bool EventCutterReco::ef_mu24Hypo(float eta, float pt) {
    
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

//trig_EF_trigmuonef
bool EventCutterReco::ef_mu36Hypo(float eta, float pt) {
    
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
	//sel.cutFlow().clear();
	using KinematicUtils::deltaR;
	int Events_before_trigMatching = 0;
	int Events_after_trigMatching = 0;
	//std::cout<<"m_loose= "<<m_loose<<std::endl;
	if (m_loose == false) {	//tight selection: use --doLoose 0 for preselect
		//sel.cutFlow().clear();
		int els = 0;
		int mus = 0;
		int fatjets = 0;
		int jets = 0;
		int btags = 0;	
		// CUT 0
		sel.cutFlow()[0] += 1;
		if (e.isData()) {
			if (!DQ::PassRunLB(e.runNumber(), e.lbn())){return false;}
			}
		if (e.terr() == 2) return false;
		if ((e.cfl() & 0x40000) != 0) return false;
		if (e.lerr() > 1) return false;
		
		// CUT 1 - Quality
		sel.cutFlow()[1] +=1;
		
		if(sel.npv() < 1) return false;
		// CUT 2 - Vertex
		sel.cutFlow()[2] += 1;
		
		// preselect jets for OR
		bool passIsBadLoose = true; //this flag is dor Data only
		for (int k = 0; k < e.largeJet().size(); ++k) {
		}
		
		for (int k = 0; k < e.jet().size(); ++k) {
			if (!e.jet()[k].passBadLoose()) passIsBadLoose = false;
			if (e.jet()[k].pass()) {
				sel.jet().push_back(e.jet()[k]);
				jets++;
				}
			}
		if (!passIsBadLoose) {
			return false;
			}
		
		// CUT 3 - Jet  Quality
		sel.cutFlow()[3] += 1;
		
		std::vector<std::vector<Jet> > selectedJets;
		int elidx = -1;
		
		for (int k = 0; k < e.electron().size(); ++k) {
			selectedJets.push_back(std::vector<Jet>());
			if (e.electron()[k].pass()) {
				selectedJets[k] = sel.jet();
				float dr = 99.;
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
			if (!e.muon()[k].pass()) continue;
			float dr = e.muon()[k].minDeltaR(sel.jet());
			if (dr >= 0.04 + 10e3/e.muon()[k].mom().Perp()) {
				sel.muon().push_back(e.muon()[k]);
				mus++;
				}
			}
		
		// CUT 4 - Trigger
		if (e.triggerMuon() || e.triggerElectron())   sel.cutFlow()[4] += 1;
		
		if (e.triggerMuon() && els == 0 && mus == 1) {
			lepton = sel.muon()[0].mom();
			electron = false;
			} 
		else if (e.triggerElectron() && els == 1 && mus == 0) {
			electron = true;
			} else
		return false;
		
		sel.partMom() = e.partMom();
		// CUT 5 - == 1 lepton
		sel.cutFlow()[5] += 1;
		sel.triggerElectron() = e.triggerElectron();
		sel.triggerMuon() = e.triggerMuon();
		// CUT 6 Trigger
		sel.cutFlow()[6]+=1;
		
		
		Events_before_trigMatching = Events_before_trigMatching+1;
		
		
		// Trigger Matching	
		if (e.triggerMuon() && els == 0 && mus == 1) {
			lepton = sel.muon()[0].mom();		
			int itrk;
			// The matching cut.
			double dR_Trig_Mu_cut = 0.15;	
			// Check the trigger matching using the EF trigger object.
			double dR_Trig_Mu = 1000.0;
			//int trigPass = 0;
			bool trigPass_mu24i = 0;
			bool trigPass_mu36i = 0;
			int itrEF = 0;		
			while(dR_Trig_Mu >= dR_Trig_Mu_cut && itrEF < e.sr()->trig_EF_trigmuonef_n) {
			if(e.sr()->trig_EF_trigmuonef_EF_mu24i_tight->at(itrEF) <= 0 && e.sr()->trig_EF_trigmuonef_EF_mu36_tight->at(itrEF) <= 0) {
				itrEF++; continue;
				}
			itrk = 0;
			while(dR_Trig_Mu >= dR_Trig_Mu_cut && itrk < e.sr()->trig_EF_trigmuonef_track_n->at(itrEF)) {
				bool passedHypo=false;
				if (e.sr()->trig_EF_trigmuonef_EF_mu24i_tight->at(itrEF) && ef_mu24Hypo(fabs(e.sr()->trig_EF_trigmuonef_track_CB_eta->at(itrEF).at(itrk)),e.sr()->trig_EF_trigmuonef_track_CB_pt->at(itrEF).at(itrk)/1.e+3)) passedHypo=true;
				if (e.sr()->trig_EF_trigmuonef_EF_mu36_tight->at(itrEF) && ef_mu36Hypo(fabs(e.sr()->trig_EF_trigmuonef_track_CB_eta->at(itrEF).at(itrk)),e.sr()->trig_EF_trigmuonef_track_CB_pt->at(itrEF).at(itrk)/1.e+3)) passedHypo=true;
				if (!passedHypo) { itrk++; continue;}
				dR_Trig_Mu = deltaR(sel.muon()[0].mom().Eta(), e.sr()->trig_EF_trigmuonef_track_CB_eta->at(itrEF).at(itrk), sel.muon()[0].mom().Phi(), e.sr()->trig_EF_trigmuonef_track_CB_phi->at(itrEF).at(itrk));
				if (dR_Trig_Mu < dR_Trig_Mu_cut && e.sr()->trig_EF_trigmuonef_EF_mu24i_tight) trigPass_mu24i =1;//|= (1 << 0); // 1
				if (dR_Trig_Mu < dR_Trig_Mu_cut && e.sr()->trig_EF_trigmuonef_EF_mu36_tight) trigPass_mu36i =1;//|= (1 << 1); // 2
				itrk++;
				}
			itrEF++;
			}
			if(trigPass_mu36i==0 && trigPass_mu24i==0 ) return false;
			if(dR_Trig_Mu >= dR_Trig_Mu_cut) return false;
			Events_after_trigMatching = Events_after_trigMatching+1;
			}
		//End of Trigger Matching
		//CUT 7 Muon Trigg matching
	    sel.cutFlow()[7] += 1;
	    
		//New Electron trigger matching	    
		if (e.triggerElectron() && els == 1 && mus == 0) {
			lepton = sel.electron()[0].mom();
			bool trigPass_e24vhi = 0;
			bool trigPass_e60 = 0;
			double dR_Trig_El, dR_Trig_El_cut;
			int itrEF_el = 0;
			dR_Trig_El_cut = 0.15;
			// Loop over the trigger objects and find the one which is closest to the offline electron.
			dR_Trig_El = 1000.0;
			//int trigPass = 0; 
			// VD: in mc11b/c the MC menu is correctly simulated for the mc so we can use the same triggers in data and MC
			//     at present using and OR of EF_e22vh_medium1 and EF_e45_medium1 to be fully efficienct also for hight Pt analyses
			while(dR_Trig_El >= dR_Trig_El_cut && itrEF_el < e.sr()->trig_EF_el_n) {
				if(e.sr()->trig_EF_el_EF_e24vhi_medium1 <= 0  && e.sr()->trig_EF_el_EF_e60_medium1 <= 0) {// did not trigger this EF alg.
					itrEF_el++; continue;
					}
				dR_Trig_El = deltaR(sel.electron()[0].mom().Eta(), e.sr()->trig_EF_el_eta->at(itrEF_el), sel.electron()[0].mom().Phi(), e.sr()->trig_EF_el_phi->at(itrEF_el));
				if (dR_Trig_El < dR_Trig_El_cut && e.sr()->trig_EF_el_EF_e24vhi_medium1) trigPass_e24vhi = 1;// |= (1 << 0); // 1
				if (dR_Trig_El < dR_Trig_El_cut && e.sr()->trig_EF_el_EF_e60_medium1) trigPass_e60 = 1; //trigPass |= (1 << 1); // 2
					itrEF_el++;
					}
			if(trigPass_e60==0 && trigPass_e24vhi==0 ) return false;
			if(dR_Trig_El >= dR_Trig_El_cut) return 0; // check association cut.
			//return trigPass; 
			}    
        //CUT 8 El Trigg matching
	    sel.cutFlow()[8] += 1;   
	    //End Electron Trigger matching
	
		//all in MeV
		sel.met(e.met().Px(), e.met().Py());
		float mtw = std::sqrt(2*lepton.Perp()*sel.met().Perp()*(1 - std::cos(lepton.Phi() - sel.met().Phi())) );		
		
		// CUT 9 - MET . 20 GeV
	    if (e.met().Perp() > 20e3) sel.cutFlow()[9] += 1;
		if (e.met().Perp() < 20e3 || e.met().Perp() + mtw < 60e3)  return false;
		//if (e.met().Perp() < 30e3 || e.met().Perp() + mtw < 60e3)  return false;
		// CUT10 MET+MTW > 60 GeV
		sel.cutFlow()[10] += 1;
		
		// require a jet close to the lepton (dR < 1.5) and select the highest p_T one
		// for future cut
		double dr = 1.5;
		int sel_jet = -1;
		
		double lpt = 0;
		for (int k = 0; k < sel.jet().size(); ++k) {
			double tdr = sel.jet()[k].mom().DeltaR(lepton);
			if (tdr < dr && sel.jet()[k].mom().Perp() > lpt) {
				sel_jet = k;
				lpt = sel.jet()[k].mom().Perp();
				}
			}
		if (sel_jet == -1) return false;
		
		// CUT11 dR(l,j) < 1.5
		sel.cutFlow()[11] += 1;
		
		for (int k = 0; k < e.largeJet().size(); ++k) {
			
			//if (e.largeJet()[k].pass()) 
			if (e.largeJet()[k].passLoose()) {
				// now check topological cuts
				//the 1st SM Top which is so boosted that it's a fat jet
				//A.DeltaR( B )= DeltaR(A,B) if a and B are TLorentzVector: sel.jet()[sel_jet].mom().DeltaR(e.largeJet()[k].mom()) = DeltaR(jet,fatjet)> 1.5
				bool passdR = sel.jet()[sel_jet].mom().DeltaR(e.largeJet()[k].mom()) > 1.5; //1.5
				//the 2nd SM top is boosted nut not enough to be a fat jet==> lepton and b are well seperated
				bool passdPhi = true;
				//e.largeJet()[k].mom().DeltaPhi(lepton) > 2.3; //2.3
				bool passM = e.largeJet()[k].mom().M() > 60e3; //100e3
				bool passD12 = e.largeJet()[k].split12() > 20e3;
				if (passdR && passdPhi && passM && passD12) {
					sel.largeJet().push_back(e.largeJet()[k]);
					fatjets++;
					}
				}
			}
		
		//  if (fatjets < 1) return false;
		if (jets < 3 && fatjets < 1) return false;
		
		
		// CUT12 1 fatjet
		sel.cutFlow()[12] += 1;
		
		for (int k = 0; k < sel.jet().size(); ++k) {
			if (sel.jet()[k].btag()) {
			btags++;
			}
			}
		if (btags < 1) {
			return false;
			}
		// CUT 13 : 1 b-tag
		sel.cutFlow()[13] += 1;
		return true;
		}
	
	//For QCD
	if (m_loose == true) {
		int els = 0;
		int mus = 0;
		int fatjets = 0;
		int jets = 0;
		int btags = 0;
		//sel.cutFlow().clear();
		// CUT 14
		sel.cutFlow()[20] += 1;
		if (e.isData()) {
			if (!DQ::PassRunLB(e.runNumber(), e.lbn())){return false;}
			}
		if (e.terr() == 2) return false;
		if ((e.cfl() & 0x40000) != 0) return false;
		if (e.lerr() > 1) return false;
		
		// CUT 15 - Quality
		sel.cutFlow()[21] += 1;
		if(sel.npv() < 1) return false;
		// CUT 16 - Vertex
		sel.cutFlow()[22] += 1;
		
		// preselect jets for OR
		bool passIsBadLoose = true; //this flag is dor Data only
		
		for (int k = 0; k < e.jet().size(); ++k) {
			if (!e.jet()[k].passBadLoose()) passIsBadLoose = false;
			if (e.jet()[k].pass()) {
			sel.jet().push_back(e.jet()[k]);
			jets++;
			}
			}
		if (!passIsBadLoose) {
			return false;
			}
		
		// CUT 17 - Jet  Quality
		sel.cutFlow()[23] += 1;
		for (int k = 0; k < e.electron().size(); ++k) {
			if (e.electron()[k].passLoose()) {
				sel.electron().push_back(e.electron()[k]);
				els++;
				}
			}
		
		for (int k = 0; k < e.muon().size(); ++k) {
			if (!e.muon()[k].passLoose()) continue;
			float dr = e.muon()[k].minDeltaR(sel.jet());
			if (dr >= 0.04 + 10e3/e.muon()[k].mom().Perp()) {
				sel.muon().push_back(e.muon()[k]);
				mus++;
				}
			}
		
		// CUT 18 - Trigger
		if (e.triggerMuon() || e.triggerElectron())   sel.cutFlow()[24] += 1;
		TLorentzVector lepton;
		bool electron = true;
		if (e.triggerMuon() && els == 0 && mus == 1) {
			lepton = sel.muon()[0].mom();
			electron = false;
			} else if (e.triggerElectron() && els == 1 && mus == 0) {
			lepton = sel.electron()[0].mom();
			electron = true;
			} else
		return false;
		
 	
		sel.partMom() = e.partMom();
		// CUT 19 - == 1 lepton
		sel.cutFlow()[25] += 1;
		sel.triggerElectron() = e.triggerElectron();
		sel.triggerMuon() = e.triggerMuon();
		
		// CUT 20 - == Trigger
		sel.cutFlow()[26] += 1;
		
		//New Trigger matching
		int itrk;
		// The matching cut.
		double dR_Trig_Mu_cut = 0.15;	
		// Check the trigger matching using the EF trigger object.
		double dR_Trig_Mu = 1000.0;
		//int trigPass = 0;
		bool trigPass_mu24i = 0;
		bool trigPass_mu36i = 0;
		int itrEF = 0;	
		if (e.triggerMuon() && els == 0 && mus == 1){ 	
			while(dR_Trig_Mu >= dR_Trig_Mu_cut && itrEF < e.sr()->trig_EF_trigmuonef_n) {
				if(e.sr()->trig_EF_trigmuonef_EF_mu24i_tight->at(itrEF) <= 0 && e.sr()->trig_EF_trigmuonef_EF_mu36_tight->at(itrEF) <= 0) {
					itrEF++; continue;
					}
				itrk = 0;
				
				while(dR_Trig_Mu >= dR_Trig_Mu_cut && itrk < e.sr()->trig_EF_trigmuonef_track_n->at(itrEF)) {
					bool passedHypo=false;
					
					if (e.sr()->trig_EF_trigmuonef_EF_mu24i_tight->at(itrEF) && ef_mu24Hypo(fabs(e.sr()->trig_EF_trigmuonef_track_CB_eta->at(itrEF).at(itrk)),e.sr()->trig_EF_trigmuonef_track_CB_pt->at(itrEF).at(itrk)/1.e+3)) passedHypo=true;
					
					if (e.sr()->trig_EF_trigmuonef_EF_mu36_tight->at(itrEF) && ef_mu36Hypo(fabs(e.sr()->trig_EF_trigmuonef_track_CB_eta->at(itrEF).at(itrk)),e.sr()->trig_EF_trigmuonef_track_CB_pt->at(itrEF).at(itrk)/1.e+3)) passedHypo=true;
					if (!passedHypo) { itrk++; continue;}
					
					dR_Trig_Mu = deltaR(sel.muon()[0].mom().Eta(), e.sr()->trig_EF_trigmuonef_track_CB_eta->at(itrEF).at(itrk), sel.muon()[0].mom().Phi(), e.sr()->trig_EF_trigmuonef_track_CB_phi->at(itrEF).at(itrk));
					
					if (dR_Trig_Mu < dR_Trig_Mu_cut && e.sr()->trig_EF_trigmuonef_EF_mu24i_tight) trigPass_mu24i =1;//|= (1 << 0); // 1
					if (dR_Trig_Mu < dR_Trig_Mu_cut && e.sr()->trig_EF_trigmuonef_EF_mu36_tight) trigPass_mu36i =1;//|= (1 << 1); // 2
					itrk++;
					}
				itrEF++;
				}
			if(trigPass_mu36i==0 && trigPass_mu24i==0 ) return false;
			
			if(dR_Trig_Mu >= dR_Trig_Mu_cut) return false; //??????? check association cut.
			
			}
		//CUT 21 Mu Trigger matching
		sel.cutFlow()[27] += 1 ;
		// 	//End Trigger matching	
		//New Electron trigger matching	    
		if (e.triggerElectron() && els == 1 && mus == 0) {
			lepton = sel.electron()[0].mom();
			bool trigPass_e24vhi = 0;
			bool trigPass_e60 = 0;
			double dR_Trig_El, dR_Trig_El_cut;
			int itrEF_el = 0;
			dR_Trig_El_cut = 0.15;
			// Loop over the trigger objects and find the one which is closest to the offline electron.
			dR_Trig_El = 1000.0;
			//int trigPass = 0; 
			// VD: in mc11b/c the MC menu is correctly simulated for the mc so we can use the same triggers in data and MC
			//     at present using and OR of EF_e22vh_medium1 and EF_e45_medium1 to be fully efficienct also for hight Pt analyses
			while(dR_Trig_El >= dR_Trig_El_cut && itrEF_el < e.sr()->trig_EF_el_n) {
				if(e.sr()->trig_EF_el_EF_e24vhi_medium1 <= 0  && e.sr()->trig_EF_el_EF_e60_medium1 <= 0) {// did not trigger this EF alg.
					itrEF_el++; continue;
					}
				dR_Trig_El = deltaR(sel.electron()[0].mom().Eta(), e.sr()->trig_EF_el_eta->at(itrEF_el), sel.electron()[0].mom().Phi(), e.sr()->trig_EF_el_phi->at(itrEF_el));
				if (dR_Trig_El < dR_Trig_El_cut && e.sr()->trig_EF_el_EF_e24vhi_medium1) trigPass_e24vhi = 1;// |= (1 << 0); // 1
				if (dR_Trig_El < dR_Trig_El_cut && e.sr()->trig_EF_el_EF_e60_medium1) trigPass_e60 = 1; //trigPass |= (1 << 1); // 2
					itrEF_el++;
					}
			if(trigPass_e60==0 && trigPass_e24vhi==0 ) return false;
			if(dR_Trig_El >= dR_Trig_El_cut) return 0; // check association cut.
			//return trigPass; 
			}    
        //CUT 8 El Trigg matching
	    sel.cutFlow()[28] += 1;   
	    //End Electron Trigger matching

		//all in MeV
		sel.met(e.met().Px(), e.met().Py());
		float mtw = std::sqrt(2*lepton.Perp()*sel.met().Perp()*(1 - std::cos(lepton.Phi() - sel.met().Phi())) );
		
		// CUT 22 - MET . 20 GeV
		if (e.met().Perp() > 20e3) sel.cutFlow()[29] += 1;
	        if (e.met().Perp() < 20e3 || e.met().Perp() + mtw < 60e3)  return false;	
	//	if (e.met().Perp() < 30e3 || e.met().Perp() + mtw < 60e3)  return false;		
		// CUT23 MET+MTW > 60 GeV
		sel.cutFlow()[30] += 1;
		
		// require a jet close to the lepton (dR < 1.5) and select the highest p_T one
		double dr = 1.5;
		int sel_jet = -1;
		
		double lpt = 0;
		for (int k = 0; k < sel.jet().size(); ++k) {
			double tdr = sel.jet()[k].mom().DeltaR(lepton);
			if (tdr < dr && sel.jet()[k].mom().Perp() > lpt) {
				sel_jet = k;
				lpt = sel.jet()[k].mom().Perp();
				}
			}
		if (sel_jet == -1) return false;
		
		// CUT24 dR(l,j) < 1.5
		sel.cutFlow()[31] += 1;
		
		for (int k = 0; k < e.largeJet().size(); ++k) {
			if (e.largeJet()[k].passLoose()) {
			bool passdR = sel.jet()[sel_jet].mom().DeltaR(e.largeJet()[k].mom()) > 1.5; //1.5
			bool passdPhi = true;
			bool passM = e.largeJet()[k].mom().M() > 60e3; //100e3
			bool passD12 = e.largeJet()[k].split12() > 20e3;
			if (passdR && passdPhi && passM && passD12) {
				sel.largeJet().push_back(e.largeJet()[k]);
				fatjets++;
				}
			}
		}
		
		if (jets < 3 && fatjets < 1) return false;
		
		// CUT25 1 fatjet
		sel.cutFlow()[32] += 1;
		
		for (int k = 0; k < sel.jet().size(); ++k) {
			if (sel.jet()[k].btag()) {
				btags++;
				}
			}
		if (btags < 1) {
			return false;
			}
		// CUT 26 : 1 b-tag
		sel.cutFlow()[33] +=1;
		return true;	
		
		}//m_loose == true
		
	//End For QCD
	
	
}
