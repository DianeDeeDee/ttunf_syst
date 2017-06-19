#include "AllCorrections.h"
#include "Event.h"
#include <fastjet/PseudoJet.hh>
#include "SkimReader.h"


void Corrections::ElectronCorr::run(Event &e) {
  tools->m_seedTool.setEventNumber(e.eventNumber());
  tools->m_eRescaler.SetRandomSeed(tools->m_seedTool.getSeed(SeedTool::forEnergyRescaler));
  for (int i = 0; i < e.electron().size(); ++i) {

    Electron &el = e.electron()[i];
    //el.mom().SetPtEtaPhiE(std::sqrt(std::pow(el.mom().caloMom().E(), 2) - std::pow(el.mom().M(), 2))/std::cosh(eta), eta, phi, el.mom().caloMom().E());

    double mass = 0.511;
    //double mass = el.mom().M();
    double eta = el.trkMom().Eta();
    double phi = el.trkMom().Phi();
  
    if (el.nSiHits() < 4) {
      eta = el.caloMom().Eta();
      phi = el.caloMom().Phi();
    }

    double f = 1;
    double fUp = 1;
    double fDown = 1;
    if (!e.isData()) {
      f = tools->m_eRescaler.getSmearingCorrection(el.caloMom().Eta(), el.caloMom().E(), egRescaler::EnergyRescalerUpgrade::NOMINAL);
      if (tools->m_syst) {
        fUp = tools->m_eRescaler.getSmearingCorrection(el.caloMom().Eta(), el.caloMom().E(), egRescaler::EnergyRescalerUpgrade::ERR_UP);
        fDown = tools->m_eRescaler.getSmearingCorrection(el.caloMom().Eta(), el.caloMom().E(), egRescaler::EnergyRescalerUpgrade::ERR_DOWN);
      }
    }

    double r = 0;
    double rUp = 0;
    double rDown = 0;
    if (e.isData()) {
      r = tools->m_eRescaler.applyEnergyCorrection(el.caloMom().Eta(), el.caloMom().E()*f, egRescaler::EnergyRescalerUpgrade::Electron, egRescaler::EnergyRescalerUpgrade::Nominal, 1, e.runNumber())/(el.caloMom().E()*f) - 1;
    }

    if (!e.isData()) {
      double totShiftUp = 0;
      double totShiftDown = 0;

      if (tools->m_syst) {
        rUp = tools->m_eRescaler.applyEnergyCorrection(el.caloMom().Eta(), el.caloMom().E()*f, egRescaler::EnergyRescalerUpgrade::Electron, egRescaler::EnergyRescalerUpgrade::AllUp, 1, e.runNumber())/(el.caloMom().E()*f) - 1;
        rDown = tools->m_eRescaler.applyEnergyCorrection(el.caloMom().Eta(), el.caloMom().E()*f, egRescaler::EnergyRescalerUpgrade::Electron, egRescaler::EnergyRescalerUpgrade::AllDown, 1, e.runNumber())/(el.caloMom().E()*f) - 1;

        double zeeUp = tools->m_eRescaler.applyEnergyCorrection(el.caloMom().Eta(), el.caloMom().E()*f, egRescaler::EnergyRescalerUpgrade::Electron, egRescaler::EnergyRescalerUpgrade::ZeeAllUp, 1, e.runNumber())/(el.caloMom().E()*f) - 1;
        double matUp = tools->m_eRescaler.applyEnergyCorrection(el.caloMom().Eta(), el.caloMom().E()*f, egRescaler::EnergyRescalerUpgrade::Electron, egRescaler::EnergyRescalerUpgrade::R12StatUp, 1, e.runNumber())/(el.caloMom().E()*f) - 1;
        double psUp = tools->m_eRescaler.applyEnergyCorrection(el.caloMom().Eta(), el.caloMom().E()*f, egRescaler::EnergyRescalerUpgrade::Electron, egRescaler::EnergyRescalerUpgrade::PSStatUp, 1, e.runNumber())/(el.caloMom().E()*f) - 1;
        double lowUp = tools->m_eRescaler.applyEnergyCorrection(el.caloMom().Eta(), el.caloMom().E()*f, egRescaler::EnergyRescalerUpgrade::Electron, egRescaler::EnergyRescalerUpgrade::LowPtUp, 1, e.runNumber())/(el.caloMom().E()*f) - 1;

        //double zeeDown = tools->m_eRescaler.applyEnergyCorrection(el.caloMom().Eta(), el.caloMom().E()*f, egRescaler::EnergyRescalerUpgrade::Electron, egRescaler::EnergyRescalerUpgrade::ZeeAllDown, 1, e.runNumber())/(el.caloMom().E()*f) - 1;
        //double matDown = tools->m_eRescaler.applyEnergyCorrection(el.caloMom().Eta(), el.caloMom().E()*f, egRescaler::EnergyRescalerUpgrade::Electron, egRescaler::EnergyRescalerUpgrade::R12StatDown, 1, e.runNumber())/(el.caloMom().E()*f) - 1;
        //double psDown = tools->m_eRescaler.applyEnergyCorrection(el.caloMom().Eta(), el.caloMom().E()*f, egRescaler::EnergyRescalerUpgrade::Electron, egRescaler::EnergyRescalerUpgrade::PSStatDown, 1, e.runNumber())/(el.caloMom().E()*f) - 1;
        //double lowDown = tools->m_eRescaler.applyEnergyCorrection(el.caloMom().Eta(), el.caloMom().E()*f, egRescaler::EnergyRescalerUpgrade::Electron, egRescaler::EnergyRescalerUpgrade::LowPtDown, 1, e.runNumber())/(el.caloMom().E()*f) - 1;

        totShiftUp = std::sqrt( std::pow(zeeUp, 2) + std::pow(matUp, 2) + std::pow(psUp, 2) + std::pow(lowUp, 2) );
        totShiftDown = -totShiftUp;//std::sqrt( std::pow(zeeDown, 2) + std::pow(matDown, 2) + std::pow(psDown, 2) + std::pow(lowDown, 2) );
      }

      double extra = 1;
      if (e.atlFastII()) {
        extra = tools->m_eRescaler.applyAFtoG4(el.caloMom().Eta());
      }
      r = (r+1)*extra - 1;
      rUp = (rUp+1)*extra - 1;
      rDown = (rDown+1)*extra - 1;
    }
  
    TLorentzVector v;
    v.SetPtEtaPhiM(el.caloMom().E()*f*(r+1)/std::cosh(eta), eta, phi, mass);
    TLorentzVector vUp;
    vUp.SetPtEtaPhiM(el.caloMom().E()*fUp*(r+1)/std::cosh(eta), eta, phi, mass);
    TLorentzVector vDown;
    vDown.SetPtEtaPhiM(el.caloMom().E()*fDown*(r+1)/std::cosh(eta), eta, phi, mass);

    el.setCorrection(v); // set it to the default
    el.corr("", true) = v;
    el.corr("eSmearUp", true) = vUp;
    el.corr("eSmearDown", true) = vDown;

    TLorentzVector vrUp;
    vrUp.SetPtEtaPhiM(el.caloMom().E()*f*(rUp+1)/std::cosh(eta), eta, phi, mass);
    TLorentzVector vrDown;
    vrDown.SetPtEtaPhiM(el.caloMom().E()*f*(rDown+1)/std::cosh(eta), eta, phi, mass);

    el.corr("eRescaleUp", true) = vrUp;
    el.corr("eRescaleDown", true) = vrDown;
  }
}

void Corrections::MuonCorr::run(Event &e) {
  if (e.isData()) return;
  for (int i = 0; i < e.muon().size(); ++i) {
    Muon &mu = e.muon()[i];

    //tools->m_mcpSmear.SetSeed((int) (mu.mom().Phi()*10000), i); 
    tools->m_mcpSmear.SetSeed(tools->m_seedTool.getSeed(SeedTool::forMuonSmear), i);

    double m = mu.mom().M();
    double eta = mu.mom().Eta();
    double phi = mu.mom().Phi();
    double ptcb = mu.mom().Perp();
    double newPt(ptcb), newCB(ptcb), newME(mu.momME().Perp()), newID(mu.momID().Perp());
    double newCB_idup = ptcb;
    double newCB_idlow = ptcb;
    double newCB_msup = ptcb;
    double newCB_mslow = ptcb;
    double newCB_scale = ptcb;
    double newID_idup = newID;
    double newID_idlow = newID;
    double newID_msup = newID;
    double newID_mslow = newID;
    double newID_scale = newID;
    double newME_idup = newME;
    double newME_idlow = newME;
    double newME_msup = newME;
    double newME_mslow = newME;
    double newME_scale = newME;

    // Get Smeared Pts 
    double pTCB_smeared; 
    double pTMS_smeared; 
    double pTID_smeared;

    if (std::fabs(eta) <= 2.7) {
      if (mu.st()) {
        if (!e.isData()) {
          tools->m_mcpSmear.Event(mu.momID().Perp(), eta, "ID", mu.charge(), phi);
          newPt = newCB = newCB_idup = newCB_idlow = tools->m_mcpSmear.pTID();
          newME = 1e-20; // should be 0 anyway
          newID = newID_idup = newID_idlow = tools->m_mcpSmear.pTID();
          if (tools->m_syst) {
            tools->m_mcpSmear.PTVar(newCB_idup, "IDUP");
            tools->m_mcpSmear.PTVar(newCB_idlow, "IDLOW");
            tools->m_mcpSmear.PTVar(newID_idup, "IDUP");
            tools->m_mcpSmear.PTVar(newID_idlow, "IDLOW");
          }
        }
      } else if (mu.sa()) {
        if (!e.isData()) {
          tools->m_mcpSmear.Event(mu.momME().Perp(), eta, "MS", mu.charge(), phi);
          newPt = newCB = newCB_msup = newCB_mslow = tools->m_mcpSmear.pTMS();
          newME = newME_msup = newME_mslow = tools->m_mcpSmear.pTMS();
          newID = 1e-20;
          if (tools->m_syst) {
            tools->m_mcpSmear.PTVar(newCB_msup, "MSUP");
            tools->m_mcpSmear.PTVar(newCB_mslow, "MSLOW");
            tools->m_mcpSmear.PTVar(newME_msup, "MSUP");
            tools->m_mcpSmear.PTVar(newME_mslow, "MSLOW");
          }
        }
      } else if (mu.cb()) {
        if (!e.isData()) {
          tools->m_mcpSmear.Event(mu.momME().Perp(), mu.momID().Perp(), ptcb, eta, mu.charge(), phi);
          newPt = newCB = newCB_mslow = newCB_msup = newCB_idup = newCB_idlow = (mu.momME().Perp() < 1e-19)?tools->m_mcpSmear.pTID():tools->m_mcpSmear.pTCB();
          newME = newME_msup = newME_mslow = newME_idup = newME_idlow = tools->m_mcpSmear.pTMS();
          newID = newID_msup = newID_mslow = newID_idup = newID_idlow = tools->m_mcpSmear.pTID();
          if (tools->m_syst) {
            tools->m_mcpSmear.PTVar(newME_msup, newID_msup, newCB_msup, "MSUP");
            tools->m_mcpSmear.PTVar(newME_mslow, newID_mslow, newCB_mslow, "MSLOW");
            tools->m_mcpSmear.PTVar(newME_idup, newID_idup, newCB_idup, "IDUP");
            tools->m_mcpSmear.PTVar(newME_idlow, newID_idlow, newCB_idlow, "IDLOW");
            tools->m_mcpSmear.PTVar(newME_scale, newID_scale, newCB_scale, "SCALEUP");
            if (mu.momME().Perp() < 1e-19) {
              newCB_msup = newID_msup;
              newCB_mslow = newID_mslow;
              newCB_idup = newID_idup;
              newCB_idlow = newID_idlow;
              newCB_scale = newID_scale;
            }
          }
        }
      }
    }
  
    // Use the smeared pT but keep the invariant mass the same
    m = 105.658;
    TLorentzVector v; v.SetPtEtaPhiM(newCB, eta, phi, m);
    mu.setCorrection(v);
    mu.corr("", true) = v;
    mu.momMECorr().SetPtEtaPhiM(newME, eta, phi, m);

    TLorentzVector vIDUP; vIDUP.SetPtEtaPhiM(newCB_idup, eta, phi, m);
    TLorentzVector vIDLOW; vIDLOW.SetPtEtaPhiM(newCB_idlow, eta, phi, m);
    TLorentzVector vMSUP; vMSUP.SetPtEtaPhiM(newCB_msup, eta, phi, m);
    TLorentzVector vMSLOW; vMSLOW.SetPtEtaPhiM(newCB_mslow, eta, phi, m);
    TLorentzVector vSCALEUP; vSCALEUP.SetPtEtaPhiM(newCB_scale, eta, phi, m);
    mu.corr("muSmearIDUP", true) = vIDUP;
    mu.corr("muSmearIDLOW", true) = vIDLOW;
    mu.corr("muSmearMSUP", true) = vMSUP;
    mu.corr("muSmearMSLOW", true) = vMSLOW;
    mu.corr("muSmearSCALEUP", true) = vSCALEUP;

  }
}

// for the e-jet OR with e subtraction from jets
void Corrections::ElJetOR::run(Event &e) {
  if (!e.sr()) return;

  std::vector<float> el_cl_E,el_cl_eta,el_cl_phi,el_trk_eta,el_trk_phi;
  std::vector<int> el_GSF_trk_index;
  std::vector<int> idx;
  for (int j = 0; j < e.electron().size(); ++j) {
    e.electron()[j].passOR() = true;
    if (e.electron()[j].pass()) {
      idx.push_back(j);
      el_cl_E.push_back(e.electron()[j].caloMom().E());
      el_cl_eta.push_back(e.electron()[j].caloMom().Eta());
      el_cl_phi.push_back(e.electron()[j].caloMom().Phi());
      el_trk_eta.push_back(e.electron()[j].trkMom().Eta());
      el_trk_phi.push_back(e.electron()[j].trkMom().Phi());
      el_GSF_trk_index.push_back(e.electron()[j].GSF_trk_index());
    }
  }

  std::vector<float> jet_pt, jet_eta, jet_phi, jet_E;
  std::vector<float> jet_jvf, jet_sumPtTrk;
  std::vector< std::vector<int> > jet_TrackAssoc;
  for (int i = 0; i < e.jet().size(); ++i) {
    Jet &j = e.jet().at(i);
    jet_pt.push_back(j.mom().Perp());
    jet_eta.push_back(j.mom().Eta());
    jet_phi.push_back(j.mom().Phi());
    jet_E.push_back(j.mom().E());
    jet_jvf.push_back(j.jvf());
    jet_sumPtTrk.push_back(j.sumPtTrk());
    jet_TrackAssoc.push_back(j.TrackAssoc_index());
  }

  // all jets from D3PD
  tools->fElJetOverlapTool.LoadJets(
  jet_pt,
  jet_eta,
  jet_phi,
  jet_E,
  jet_jvf,
  jet_sumPtTrk,
  jet_TrackAssoc);
  // ONLY SELECTED ELECTRONS
  tools->fElJetOverlapTool.LoadSelectedElectrons(
  el_cl_E, el_cl_eta, el_cl_phi,
  el_trk_eta, el_trk_phi,
  el_GSF_trk_index);
  // all tracks from D3PD
  tools->fElJetOverlapTool.LoadTracks(
  *(e.sr()->GSF_trk_trk_index),
  *(e.sr()->trk_pt),
  *(e.sr()->trk_eta),
  *(e.sr()->trk_phi_wrtPV),
  *(e.sr()->trk_theta_wrtPV),
  *(e.sr()->trk_z0_wrtPV),
  *(e.sr()->trk_d0_wrtPV));

  // JVF for D3PD each jet
  std::vector<float> JVFs = tools->fElJetOverlapTool.NewJetJVFs();
  // new jet TLVs (after electron subtraction)
  std::vector<TLorentzVector> JetTLVs = tools->fElJetOverlapTool.NewJetTLVs();
  // which jets pass the overlap removal
  // which electrons pass the overlap removal
  std::vector<bool> IsGoodEl = tools->fElJetOverlapTool.GoodEls();
  for (int i = 0; i < e.jet().size(); ++i) {
    Jet &j = e.jet().at(i);
    //j.corr("original", true) = j.mom();
    j.mom() = JetTLVs[i];
    j.jvf() = JVFs[i];
  }
  for (int i = 0; i < idx.size(); ++i) {
    int midx = idx[i];
    Electron &el = e.electron().at(midx);
    el.passOR() = IsGoodEl[i];
  }

}

//New Test
void Corrections::ElJetORLoose::run(Event &e) {
  if (!e.sr()) return;

  std::vector<float> el_cl_E,el_cl_eta,el_cl_phi,el_trk_eta,el_trk_phi;
  std::vector<int> el_GSF_trk_index;
  std::vector<int> idx;
  for (int j = 0; j < e.electron().size(); ++j) {
    e.electron()[j].passOR() = true;
    if (e.electron()[j].passLoose()) {
      idx.push_back(j);
      el_cl_E.push_back(e.electron()[j].caloMom().E());
      el_cl_eta.push_back(e.electron()[j].caloMom().Eta());
      el_cl_phi.push_back(e.electron()[j].caloMom().Phi());
      el_trk_eta.push_back(e.electron()[j].trkMom().Eta());
      el_trk_phi.push_back(e.electron()[j].trkMom().Phi());
      el_GSF_trk_index.push_back(e.electron()[j].GSF_trk_index());
    }
  }

  std::vector<float> jet_pt, jet_eta, jet_phi, jet_E;
  std::vector<float> jet_jvf, jet_sumPtTrk;
  std::vector< std::vector<int> > jet_TrackAssoc;
  for (int i = 0; i < e.jet().size(); ++i) {
    Jet &j = e.jet().at(i);
    jet_pt.push_back(j.mom().Perp());
    jet_eta.push_back(j.mom().Eta());
    jet_phi.push_back(j.mom().Phi());
    jet_E.push_back(j.mom().E());
    jet_jvf.push_back(j.jvf());
    jet_sumPtTrk.push_back(j.sumPtTrk());
    jet_TrackAssoc.push_back(j.TrackAssoc_index());
  }
   tools->fElJetOverlapTool.LoadJets(
  jet_pt,
  jet_eta,
  jet_phi,
  jet_E,
  jet_jvf,
  jet_sumPtTrk,
  jet_TrackAssoc);
  tools->fElJetOverlapTool.LoadSelectedElectrons(
  el_cl_E, el_cl_eta, el_cl_phi,
  el_trk_eta, el_trk_phi,
  el_GSF_trk_index);
   tools->fElJetOverlapTool.LoadTracks(
  *(e.sr()->GSF_trk_trk_index),
  *(e.sr()->trk_pt),
  *(e.sr()->trk_eta),
  *(e.sr()->trk_phi_wrtPV),
  *(e.sr()->trk_theta_wrtPV),
  *(e.sr()->trk_z0_wrtPV),
  *(e.sr()->trk_d0_wrtPV));
   std::vector<float> JVFs = tools->fElJetOverlapTool.NewJetJVFs();
  std::vector<TLorentzVector> JetTLVs = tools->fElJetOverlapTool.NewJetTLVs();
   std::vector<bool> IsGoodEl = tools->fElJetOverlapTool.GoodEls();

  for (int i = 0; i < e.jet().size(); ++i) {
    Jet &j = e.jet().at(i);
    j.mom() = JetTLVs[i];
    j.jvf() = JVFs[i];
  }
  for (int i = 0; i < idx.size(); ++i) {
    int midx = idx[i];
    Electron &el = e.electron().at(midx);
    el.passOR() = IsGoodEl[i];
  }

}
//End New Test 11aout
void Corrections::JetCalib::run(Event &e) {
  double mu      = e.mu();
  double rho     = e.rho();
  for (int i = 0; i < e.jet().size(); ++i) {
    Jet &j = e.jet().at(i);

    double Eraw    = j.detE();
    double eta_det = j.detEta();
    double eta     = j.detEta();
    double phi     = j.detPhi();
    double m       = j.detM();
    double Ax      = j.Ax();
    double Ay      = j.Ay();
    double Az      = j.Az();
    double Ae      = j.Ae();

    double eta_origin = e.sr()->jet_AntiKt4LCTopo_EtaOrigin->at(i);
    double phi_origin = e.sr()->jet_AntiKt4LCTopo_PhiOrigin->at(i);
    double m_origin   = e.sr()->jet_AntiKt4LCTopo_MOrigin->at(i);

    double fEM3    = -999; //EM3 correction is not applied for LC jets;
    double fTile0  = -999; //Tile0 correction is not applied for LC jets;
    double nTrk    = e.sr()->jet_AntiKt4LCTopo_nTrk_pv0_1GeV->at(i);  
    double trackWIDTH = e.sr()->jet_AntiKt4LCTopo_trackWIDTH_pv0_1GeV->at(i); //This variable may also be called jet_AntiKt4LCTopo_trackWIDTH_pv0_1GeV

    double Nsegments = 0;
    if (!tools->m_isAtlFastII) {
      // Match mspn container to jet delR<0.4
      double delR = 100;
      double delta_phi;
      double delta_eta;
      int index_musp=-1;
      for (unsigned int i = 0; i < e.sr()->musp_phi->size(); i++) {
        delta_phi = fabs(phi - e.sr()->musp_phi->at(i)); //calculate the distance in phi
        if (delta_phi > TMath::Pi()) delta_phi = (2 * TMath::Pi()) - delta_phi; // always take the smaller angle (below 180°)
        delta_eta = eta - e.sr()->musp_eta->at(i); // distance in eta
        if (sqrt( pow(delta_phi,2) + pow(delta_eta,2)) < delR) {
          delR = sqrt( pow(delta_phi,2) + pow(delta_eta,2));
          index_musp = i;
        }
      }
      if (delR >= 0.4) {
        index_musp = -1;
      }
      if (index_musp >= 0)
        Nsegments = e.sr()->musp_innerSegments->at(index_musp) + e.sr()->musp_outerSegments->at(index_musp) + e.sr()->musp_middleSegments->at(index_musp);
    }


    TLorentzVector v = tools->m_JESCalibTool->ApplyJetAreaOffsetOriginEtaJESGSC(Eraw,eta_det,phi,m,eta_origin,phi_origin,m_origin,Ax,Ay,Az,Ae,rho,trackWIDTH,nTrk,fTile0, fEM3, Nsegments, mu, e.npv());
    //TLorentzVector v = tools->m_JESCalibTool->ApplyJetAreaOffsetEtaJES(Eraw,eta,phi,m,Ax,Ay,Az,Ae,rho,mu,e.npv());
    j.setCorrection(v);
    j.mom() = v;
    j.corr("", true) = j.mom();
    j.corr("original", true) = j.mom();
    j.corr("originalnom", true) = j.mom();
  }


  for (int i = 0; i < e.largeJet().size(); ++i) {
  
    LargeJet &lj = e.largeJet().at(i);
    if (lj.mom().Perp() < 150e3 ) {
      TLorentzVector v;
      v.SetPtEtaPhiE(2e3,0,0,2e3);
      lj.setCorrection(v);
      lj.corr("", true) = lj.mom();
      continue;
    }

    double Eraw    = lj.detE();
    double eta_det = lj.detEta();
    double eta     = lj.detEta();
    double phi     = lj.detPhi();
    double m       = lj.detM();
    TLorentzVector v = tools->m_fatJESCalibTool->ApplyEtaMassJES(Eraw,eta_det, eta, phi, m);
    float jpt = v.Perp();
    if ( (jpt != jpt) || (jpt < 2e3) ) {
      v.SetPtEtaPhiE(2e3,0,0,2e3);
      continue;
    }
    lj.setCorrection(v);
    lj.corr("", true) = lj.mom();

    // Calibrate sub-jets
    for (int k = 0; k < lj.subjet().size(); ++k) {
      //fastjet::PseudoJet &fsj = lj.subjet()[k];
      //fastjet::PseudoJet fsj_areavec = fsj.area_4vector();
      TLorentzVector &sj = lj.subjet()[k].v;
      TLorentzVector &sj_a = lj.subjet_area()[k];
      TLorentzVector &sjj = lj.subjetJER()[k].v;
      sjj = sj;
      double Eraw    = sj.E();
      double eta     = sj.Eta();
      double phi     = sj.Phi();
      double m       = sj.M();
      double Ax      = sj_a.Px();
      double Ay      = sj_a.Py();
      double Az      = sj_a.Pz();
      double Ae      = sj_a.E();
      TLorentzVector calib_jet = tools->m_subjetCalibTool->ApplyJetAreaOffsetEtaJES(Eraw, eta, phi, m, Ax, Ay, Az, Ae, rho, mu, e.npv());
      double shift = 0;
      //if (calib_jet.M() < 100) shift = 100;
      calib_jet.SetPtEtaPhiE(calib_jet.Perp(), calib_jet.Eta(), calib_jet.Phi(), calib_jet.E()+shift);
      lj.subjet()[k].v = calib_jet;
      //lj.subjet()[k].reset_momentum(calib_jet.Px(), calib_jet.Py(), calib_jet.Pz(), calib_jet.E());

      // apply JER on subjets
      std::pair<float, float> su= tools->GetJERSmearFactors(0.15, sjj.Eta(), sjj.Perp()*1e-3);
      float S = su.first;
      float U = su.second;
      float smearingFactorSyst = std::sqrt(std::pow(S+U, 2) - std::pow(S, 2));
      if (smearingFactorSyst > 1) smearingFactorSyst = 1;
      float smearingGaus = tools->r.Gaus(1, smearingFactorSyst);
      while (smearingGaus < 0) {
        smearingGaus = tools->r.Gaus(1, smearingFactorSyst);
      }
      sjj.SetPtEtaPhiE(sjj.Perp()*smearingGaus, sjj.Eta(), sjj.Phi(), sjj.E()*smearingGaus);

      Eraw    = sjj.E();
      eta     = sjj.Eta();
      phi     = sjj.Phi();
      m       = sjj.M();
      calib_jet = tools->m_subjetCalibTool->ApplyJetAreaOffsetEtaJES(Eraw, eta, phi, m, Ax, Ay, Az, Ae, rho, mu, e.npv());
      calib_jet.SetPtEtaPhiE(calib_jet.Perp(), calib_jet.Eta(), calib_jet.Phi(), calib_jet.E()+shift);
      lj.subjetJER()[k].v = calib_jet;
    }
  }


 for (int i = 0; i < e.largeJetBB().size(); ++i) {
  
    LargeJet &lj = e.largeJetBB().at(i);
  
    if (lj.mom().Perp() < 150e3 ) {
      TLorentzVector v;
      v.SetPtEtaPhiE(2e3,0,0,2e3);
      lj.setCorrection(v);
      lj.corr("", true) = lj.mom();
      continue;
    }

    double Eraw    = lj.detE();
    double eta_det = lj.detEta();
    double eta     = lj.detEta();
    double phi     = lj.detPhi();
    double m       = lj.detM();
    TLorentzVector v = tools->m_fatJESBBCalibTool->ApplyEtaMassJES(Eraw,eta_det, eta, phi, m);
    float jpt = v.Perp();
    if ( (jpt != jpt) || (jpt < 2e3) ) {
      v.SetPtEtaPhiE(2e3,0,0,2e3);
      continue;
    }
    lj.setCorrection(v);

    // Calibrate sub-jets
    for (int k = 0; k < lj.subjet().size(); ++k) {
      //fastjet::PseudoJet &fsj = lj.subjet()[k];
      //fastjet::PseudoJet fsj_areavec = fsj.area_4vector();
      TLorentzVector &sj = lj.subjet()[k].v;
      TLorentzVector &sj_a = lj.subjet_area()[k];
      TLorentzVector &sjj = lj.subjetJER()[k].v;
      sjj = sj;
      double Eraw    = sj.E();
      double eta     = sj.Eta();
      double phi     = sj.Phi();
      double m       = sj.M();
      double Ax      = sj_a.Px();
      double Ay      = sj_a.Py();
      double Az      = sj_a.Pz();
      double Ae      = sj_a.E();
      TLorentzVector calib_jet = tools->m_subjetCalibTool->ApplyJetAreaOffsetEtaJES(Eraw, eta, phi, m, Ax, Ay, Az, Ae, rho, mu, e.npv());
      double shift = 0;
      //if (calib_jet.M() < 100) shift = 100;
      calib_jet.SetPtEtaPhiE(calib_jet.Perp(), calib_jet.Eta(), calib_jet.Phi(), calib_jet.E()+shift);
      lj.subjet()[k].v = calib_jet;
      //lj.subjet()[k].reset_momentum(calib_jet.Px(), calib_jet.Py(), calib_jet.Pz(), calib_jet.E());

      // apply JER on subjets
      std::pair<float, float> su= tools->GetJERSmearFactors(0.15, sjj.Eta(), sjj.Perp()*1e-3);
      float S = su.first;
      float U = su.second;
      float smearingFactorSyst = std::sqrt(std::pow(S+U, 2) - std::pow(S, 2));
      float smearingGaus = tools->r.Gaus(1, smearingFactorSyst);
      if (smearingFactorSyst > 1) smearingFactorSyst = 1;
      while (smearingGaus < 0) {
        smearingGaus = tools->r.Gaus(1, smearingFactorSyst);
      }
      sjj.SetPtEtaPhiE(sjj.Perp()*smearingGaus, sjj.Eta(), sjj.Phi(), sjj.E()*smearingGaus);

      Eraw    = sjj.E();
      eta     = sjj.Eta();
      phi     = sjj.Phi();
      m       = sjj.M();
      calib_jet = tools->m_subjetCalibTool->ApplyJetAreaOffsetEtaJES(Eraw, eta, phi, m, Ax, Ay, Az, Ae, rho, mu, e.npv());
      calib_jet.SetPtEtaPhiE(calib_jet.Perp(), calib_jet.Eta(), calib_jet.Phi(), calib_jet.E()+shift);
      lj.subjetJER()[k].v = calib_jet;
    }
  }



}




void Corrections::BoostedJES::run(Event &e) {
  if (e.isData()) return;
  for (int i = 0; i < e.largeJet().size(); ++i) {
    LargeJet &j = e.largeJet().at(i);

    TLorentzVector jJES = j.mom();

    double pT = jJES.Perp()*1e-3;
    double eta = jJES.Eta();
    double mass = jJES.M()*1e-3;

    j.corr("boostedJER", true) = jJES;
    j.corr("boostedJESUp", true) = jJES;
    j.corr("boostedJESDown", true) = jJES;
    j.corr("boostedJMSUp", true) = jJES;
    j.corr("boostedJMSDown", true) = jJES;

    double jptsu_gj = 0;
    double jptsu_topo = 0;
    double jptsu_inter = 0;
    double shift_NPV = 0;
    double shift_Mu = 0;
    if ( ((pT > 200) && (pT < 700) && (std::fabs(eta) < 2.0)) ||
         ((pT >= 700) && (std::fabs(eta) < 2.0) && (mass/pT >= 0) && (mass/pT <= 1)) ) {
      //double jptsu_comp01 = tools->m_uJES.getRelUncertComponent("Component_dataMC", pT, mass, eta, UJUncertaintyProvider::JPTS);
      //double jptsu_comp02 = tools->m_uJES.getRelUncertComponent("Component_dPhiCut", pT, mass, eta, UJUncertaintyProvider::JPTS);
      //double jptsu_comp03 = tools->m_uJES.getRelUncertComponent("Component_pt2Cut", pT, mass, eta, UJUncertaintyProvider::JPTS);
      //double jptsu_comp04 = tools->m_uJES.getRelUncertComponent("Component_photonPurity", pT, mass, eta, UJUncertaintyProvider::JPTS);
      //double jptsu_comp05 = tools->m_uJES.getRelUncertComponent("Component_PES", pT, mass, eta, UJUncertaintyProvider::JPTS);
      //double jptsu_comp06 = tools->m_uJES.getRelUncertComponent("Component_generator", pT, mass, eta, UJUncertaintyProvider::JPTS);
      //double jptsu_comp07 = tools->m_uJES.getRelUncertComponent("Component_kterm", pT, mass, eta, UJUncertaintyProvider::JPTS);
      //double jptsu_comp08 = tools->m_uJES.getRelUncertComponent("Component_JER", pT, mass, eta, UJUncertaintyProvider::JPTS);
      //double jptsu_comp09 = tools->m_uJES.getRelUncertComponent("Component_akt4insideOutsideLargeR", pT, mass, eta, UJUncertaintyProvider::JPTS);
      //double jptsu_comp10 = tools->m_uJES.getRelUncertComponent("Component_more1smallJetInsideLargeR", pT, mass, eta, UJUncertaintyProvider::JPTS);
      //double jptsu_comp11 = tools->m_uJES.getRelUncertComponent("Component_stats", pT, mass, eta, UJUncertaintyProvider::JPTS);
      //double jptsu_comp12 = tools->m_uJES.getRelUncertComponent("Component_MoverPt", pT, mass, eta, UJUncertaintyProvider::JPTS);

      //if you are interested only in the quadratic sum of those 12 components of the gamma+jet method, you can use the function:
      jptsu_gj = tools->m_uJES.getRelUncert_GammaJet_Apr042014(pT, mass, eta);

      //jptsu_gj += std::pow(jptsu_comp01,2);
      //jptsu_gj += std::pow(jptsu_comp02,2);
      //jptsu_gj += std::pow(jptsu_comp03,2);
      //jptsu_gj += std::pow(jptsu_comp04,2);
      //jptsu_gj += std::pow(jptsu_comp05,2);
      //jptsu_gj += std::pow(jptsu_comp06,2);
      //jptsu_gj += std::pow(jptsu_comp07,2);
      //jptsu_gj += std::pow(jptsu_comp08,2);
      //jptsu_gj += std::pow(jptsu_comp09,2);
      //jptsu_gj += std::pow(jptsu_comp10,2);
      //jptsu_gj += std::pow(jptsu_comp11,2);
      //jptsu_gj += std::pow(jptsu_comp12,2);
      //jptsu_gj = std::sqrt(jptsu_gj);

      //This can be used to quickly check susceptibility to this uncertainty, but the recommendation for the final systematic is to propagate
      // each of the 12 uncertainty components individually
      //
      jptsu_topo = tools->m_uJES.getRelUncertComponent("Component_topology", pT, mass, eta, UJUncertaintyProvider::JPTS);
      jptsu_inter = tools->m_uJES.getRelUncertComponent("Component_DoubleRatioInterpolation700to900", pT, mass, eta, UJUncertaintyProvider::JPTS);
      //Pile-up effects are available separately with getRelUncert_GammaJet_ShiftNPV and getRelUncert_GammaJet_ShiftMu
      shift_NPV = tools->m_uJES.getRelUncert_GammaJet_ShiftNPV(pT, e.npv());
      shift_Mu = tools->m_uJES.getRelUncert_GammaJet_ShiftMu(pT, e.mu());
    }

    double jmsu = 0;
    if ( (pT > 0) && (fabs(eta) < 2.0) && (mass/pT >= 0) && (mass/pT <= 1) ) {
      jmsu = tools->m_uJESDoubleRatio.getRelUncert_DoubleRatio(pT, mass, eta, UJUncertaintyProvider::JMS); // NOTE: mass now passed to all functions
    }

    // uncertainties for taggers
    //double jd12u = tools->m_uJESDoubleRatio.getRelUncert_DoubleRatio(pT, mass, eta, UJUncertaintyProvider::D12);
    //double jd23u = tools->m_uJESDoubleRatio.getRelUncert_DoubleRatio(pT, mass, eta, UJUncertaintyProvider::D23);
    //double jtau21u = tools->m_uJESDoubleRatio.getRelUncert_DoubleRatio(pT, mass, eta, UJUncertaintyProvider::TAU21);
    //double jtau32u = tools->m_uJESDoubleRatio.getRelUncert_DoubleRatio(pT, mass, eta, UJUncertaintyProvider::TAU32);

    double jptsu = std::sqrt(std::pow(jptsu_gj,2) + std::pow(jptsu_topo,2) + std::pow(jptsu_inter,2) + std::pow(shift_NPV, 2) + std::pow(shift_Mu, 2));
    jJES.SetPtEtaPhiM(j.mom().Perp()*(1+jptsu), j.mom().Eta(), j.mom().Phi(), j.mom().M());

    j.corr("boostedJESUp", true) = jJES;

    jJES.SetPtEtaPhiM(j.mom().Perp()*(1-jptsu), j.mom().Eta(), j.mom().Phi(), j.mom().M());
    j.corr("boostedJESDown", true) = jJES;

    jJES.SetPtEtaPhiM(j.mom().Perp(), j.mom().Eta(), j.mom().Phi(), j.mom().M()*(1+jmsu));
    j.corr("boostedJMSUp", true) = jJES;
    jJES.SetPtEtaPhiM(j.mom().Perp(), j.mom().Eta(), j.mom().Phi(), j.mom().M()*(1-jmsu));
    j.corr("boostedJMSDown", true) = jJES;

    jJES.SetPtEtaPhiM(j.mom().Perp(), j.mom().Eta(), j.mom().Phi(), j.mom().M());
    float pt_jer = j.mom().Perp()*1e-3;
    float eta_jer = j.mom().Eta();

    float smearFactor = 0;
    if (std::fabs(eta_jer) < 2.0 && pt_jer < 1800)
      smearFactor = tools->m_largeRes.GetSmearFactor(pt_jer, eta_jer, LargeRJetResoTool::JER);

    jJES.SetPtEtaPhiM(j.mom().Perp()*(1.0+smearFactor), j.mom().Eta(), j.mom().Phi(), j.mom().M());
    j.corr("boostedJER", true) = jJES;
  }
}

void Corrections::JES::run(Event &e) {
  if (e.isData()) return;
  for (int i = 0; i < e.jet().size(); ++i) {
    Jet &j = e.jet().at(i);

    TLorentzVector jJES;
    TLorentzVector jJES_up;
    TLorentzVector jJES_down;

    float pTJES_up = j.mom().Perp();
    float pTJES_down = j.mom().Perp();

    float JESunc = 0;

    float dR = 1000;
    float fCloseBy = 0;
    for (int k = 0; k < e.jet().size(); ++k) {
      if (k == i) continue;

      if (e.jet()[k].mom().Perp() <= 12e3)
        continue;

      // jet EM+JES recalibration
      float mydr = e.jet()[k].mom().DeltaR(e.jet()[i].mom());
      if (mydr > 1.1)
        continue;

      if (mydr < dR) {
        dR = mydr;
      }
      fCloseBy += e.jet()[k].mom().Vect().Dot(e.jet()[i].mom().Vect())/(std::pow(e.jet()[k].mom().P(), 2));
    }

    float jetPt = j.mom().Perp();
    float jetEta = j.mom().Eta();
    if (jetPt < 15e3) jetPt = 15.0001e3;
    if (fabs(jetEta) >= 4.5) jetEta = (jetEta > 0)?4.49:-4.49;

    float relUncertUp   = 0;
    float relUncertDown = 0;

    float jetPhi = j.mom().Phi();

    double delR = 0;
    double delta_phi;
    double delta_eta;

    // Match mspn container to jet delR<0.4
    delR = 100;
    int index_musp=-1;
    for (unsigned int i = 0; i < e.sr()->musp_phi->size(); i++) {
      delta_phi = fabs(jetPhi - e.sr()->musp_phi->at(i)); //calculate the distance in phi
      if (delta_phi > TMath::Pi()) delta_phi = (2 * TMath::Pi()) - delta_phi; // always take the smaller angle (below 180°)
      delta_eta = jetEta - e.sr()->musp_eta->at(i); // distance in eta
      if (sqrt( pow(delta_phi,2) + pow(delta_eta,2)) < delR) {
        delR = sqrt( pow(delta_phi,2) + pow(delta_eta,2));
        index_musp = i;
      }
    }
    if (delR >= 0.4) {
      index_musp = -1;
    }
    double Nsegments = 0;
    if (index_musp >= 0)
      Nsegments = e.sr()->musp_innerSegments->at(index_musp) + e.sr()->musp_outerSegments->at(index_musp) + e.sr()->musp_middleSegments->at(index_musp);
      int Ncomp=tools->m_multiJES->getNUncertaintyComponents();
    //relUncertUp = tools->m_multiJES->getRelUncert(jetPt, jetEta, dR, true, e.npv(), e.mu(), (j.trueFlavour() == 5));
    //relUncertDown = tools->m_multiJES->getRelUncert(jetPt, jetEta, dR, false, e.npv(), e.mu(), (j.trueFlavour() == 5));
    if (jetPt*1e-3 > 20) {
		for (int icomp=1;icomp<Ncomp;++icomp) {
		TLorentzVector jJESunc_up;//new feb2016
		TLorentzVector jJESunc_down; //new feb2016
		float pTJESunc_up = j.mom().Perp();//new feb2016
		float pTJESunc_down = j.mom().Perp();//new feb2016
		TString compName = tools->m_multiJES->getComponentNames().at(icomp);
		TString compDesc = tools->m_multiJES->getComponentDescriptions().at(icomp);
		//compCategory categoryEnum = tools->m_multiJES->getComponentCategories().at(icomp); // returns a enumeration of the category for the nuisance parameter
		//String categoryName = tools->m_multiJES->getCategoryStringFromEnum(categoryEnum); // returns a string for the category name from the enumeration
		double unc = tools->m_multiJES->getRelUncertComponent(icomp, jetPt, jetEta); // nuisance paramter amplitude (with sign) ("relative uncertainty" of component)
		//  tools->m_uncJES->getRelUncertComponent(compName, jetPt, jetEta);
		// now scale your jet (the full 4-vector)
		pTJESunc_up *= 1+unc;
		pTJESunc_down *= 1-unc;
		jJESunc_up.SetPtEtaPhiM(pTJESunc_up, j.mom().Eta(), j.mom().Phi(), j.mom().M());//new feb2016
		jJESunc_down.SetPtEtaPhiM(pTJESunc_down, j.mom().Eta(), j.mom().Phi(), j.mom().M());//new feb2016
		//j.corr(Form(compName,"High") true) = jJESunc_up;//new feb2016
		//j.corr(Form(compName,"Low"), true) = jJESunc_down;//new feb2016
		//j.corr(Form("original%s", compName), true) = jJESunc_up;//new feb2016
		//j.corr(Form("original%s", compName), true) = jJESunc_down;//new feb2016
		if(icomp == 1 ){
			j.corr("EffectiveNP_Statistical1High", true) = jJESunc_up;//new feb2016
			j.corr("EffectiveNP_Statistical1Low", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Statistical1High"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Statistical1Low"), true) = jJESunc_down;//new feb2016
		}
		
		if(icomp == 2 ){
			j.corr("EffectiveNP_Statistical2High", true) = jJESunc_up;//new feb2016
			j.corr("EffectiveNP_Statistical2Low", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Statistical2High"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Statistical2Low"), true) = jJESunc_down;//new feb2016
		}
		
		if(icomp == 3 ){
			j.corr("EffectiveNP_Statistical3High", true) = jJESunc_up;//new feb2016
			j.corr("EffectiveNP_Statistical3Low", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Statistical3High"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Statistical3Low"), true) = jJESunc_down;//new feb2016
		}
		if(icomp == 4 ){
			j.corr("EffectiveNP_Statistical4High", true) = jJESunc_up;//new feb2016
			j.corr("EffectiveNP_Statistical4Low", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Statistical4High"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Statistical4Low"), true) = jJESunc_down;//new feb2016
		}
		if(icomp == 5 ){
			j.corr("EffectiveNP_Modelling1High", true) = jJESunc_up;//new feb2016
			j.corr("EffectiveNP_Modelling1Low", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Modelling1High"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Modelling1Low"), true) = jJESunc_down;//new feb2016
		}
		if(icomp == 6 ){
			j.corr("EffectiveNP_Modelling2High", true) = jJESunc_up;//new feb2016
			j.corr("EffectiveNP_Modelling2Low", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Modelling2High"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Modelling2Low"), true) = jJESunc_down;//new feb2016
		}
		if(icomp == 7 ){
			j.corr("EffectiveNP_Modelling3High", true) = jJESunc_up;//new feb2016
			j.corr("EffectiveNP_Modelling3Low", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Modelling3High"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Modelling3Low"), true) = jJESunc_down;//new feb2016
		}
		if(icomp == 8 ){
			j.corr("EffectiveNP_Modelling4High", true) = jJESunc_up;//new feb2016
			j.corr("EffectiveNP_Modelling4Low", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Modelling4High"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Modelling4Low"), true) = jJESunc_down;//new feb2016
		}
		if(icomp == 9 ){
		    j.corr("EffectiveNP_Detector1High", true) = jJESunc_up;//new feb2016
			j.corr("EffectiveNP_Detector1Low", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Detector1High"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Detector1Low"), true) = jJESunc_down;//new feb2016
		}
		if(icomp == 10 ){
			j.corr("EffectiveNP_Detector2High", true) = jJESunc_up;//new feb2016
			j.corr("EffectiveNP_Detector2Low", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Detector2High"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Detector2Low"), true) = jJESunc_down;//new feb2016
		
		}
		if(icomp == 11 ){
			j.corr("EffectiveNP_Detector3High", true) = jJESunc_up;//new feb2016
			j.corr("EffectiveNP_Detector3Low", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Detector3High"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Detector3Low"), true) = jJESunc_down;//new feb2016
		}
		if(icomp == 12 ){
			j.corr("EffectiveNP_Mixed1High", true) = jJESunc_up;//new feb2016
			j.corr("EffectiveNP_Mixed1Low", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Mixed1High"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Mixed1Low"), true) = jJESunc_down;//new feb2016
		}
		if(icomp == 13 ){
			j.corr("EffectiveNP_Mixed2High", true) = jJESunc_up;//new feb2016
			j.corr("EffectiveNP_Mixed2Low", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Mixed2High"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Mixed2Low"), true) = jJESunc_down;//new feb2016
		}
		if(icomp == 14 ){
			j.corr("EffectiveNP_Mixed3High", true) = jJESunc_up;//new feb2016
			j.corr("EffectiveNP_Mixed3Low", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Mixed3High"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Mixed3Low"), true) = jJESunc_down;//new feb2016
		}
		if(icomp == 15 ){
			j.corr("EffectiveNP_Mixed4High", true) = jJESunc_up;//new feb2016
			j.corr("EffectiveNP_Mixed4Low", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Mixed4High"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "EffectiveNP_Mixed4Low"), true) = jJESunc_down;//new feb2016
		}
		
		//Special components
		if(icomp == 16 ){
		// Eta intercalibration: theory uncertainty
			j.corr("EtaIntercalibration_ModellingHigh", true) = jJESunc_up;//new feb2016
			j.corr("EtaIntercalibration_ModellingLow", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "EtaIntercalibration_ModellingHigh"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "EtaIntercalibration_ModellingLow"), true) = jJESunc_down;//new feb2016
		}
		if(icomp == 17 ){
		// Eta intercalibration: total statistical and method uncertainty
			j.corr("EtaIntercalibration_TotalStatHigh", true) = jJESunc_up;//new feb2016
			j.corr("EtaIntercalibration_TotalStatLow", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "EtaIntercalibration_TotalStatHigh"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "EtaIntercalibration_TotalStatLow"), true) = jJESunc_down;//new feb2016
		}
		if(icomp == 18 ){
		// Pileup: Mu term
			j.corr("Pileup_OffsetMuHigh", true) = jJESunc_up;//new feb2016
			j.corr("Pileup_OffsetMuLow", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "Pileup_OffsetMuHigh"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "Pileup_OffsetMuLow"), true) = jJESunc_down;//new feb2016
		}
		
		if(icomp == 19 ){
		// Pileup: NPV term
			j.corr("Pileup_OffsetNPVHigh", true) = jJESunc_up;//new feb2016
			j.corr("Pileup_OffsetNPVLow", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "Pileup_OffsetNPVHigh"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "Pileup_OffsetNPVLow"), true) = jJESunc_down;//new feb2016
		}
		
		if(icomp == 20 ){
		// Pileup: pT term (mu part)
			j.corr("Pileup_PtTerm_MuHigh", true) = jJESunc_up;//new feb2016
			j.corr("Pileup_PtTerm_MuLow", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "Pileup_PtTerm_MuHigh"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "Pileup_PtTerm_MuLow"), true) = jJESunc_down;//new feb2016
		}
		if(icomp == 21 ){
		// Pileup: pT term (NPV part)
			j.corr("Pileup_PtTerm_NPVHigh", true) = jJESunc_up;//new feb2016
			j.corr("Pileup_PtTerm_NPVLow", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "Pileup_PtTerm_NPVHigh"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "Pileup_PtTerm_NPVLow"), true) = jJESunc_down;//new feb2016
		}
		
		if(icomp == 22 ){
		//Pileup: rho topology
			j.corr("Pileup_RhoTopologyHigh", true) = jJESunc_up;//new feb2016
			j.corr("Pileup_RhoTopologyLow", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "Pileup_RhoTopologyHigh"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "Pileup_RhoTopologyLow"), true) = jJESunc_down;//new feb2016
		}
		if(icomp == 23 ){
		// Calibration closure
			j.corr("SingleParticle_HighPtHigh", true) = jJESunc_up;//new feb2016
			j.corr("SingleParticle_HighPtLow", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "SingleParticle_HighPtHigh"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "SingleParticle_HighPtLow"), true) = jJESunc_down;//new feb2016
		}
		
		if(icomp == 24 ){
		// Calibration closure
			j.corr("RelativeNonClosure_MCTYPEHigh", true) = jJESunc_up;//new feb2016
			j.corr("RelativeNonClosure_MCTYPELow", true) = jJESunc_down;//new feb2016
			j.corr(Form("original%s", "RelativeNonClosure_MCTYPEHigh"), true) = jJESunc_up;//new feb2016
			j.corr(Form("original%s", "RelativeNonClosure_MCTYPELow"), true) = jJESunc_down;//new feb2016
		}
		//
		//            if(icomp == 19 ){
		//                j.corr("RelativeNonClosure_Pythia8High", true) = jJESunc_up;//new feb2016
		//                j.corr("RelativeNonClosure_Pythia8Low", true) = jJESunc_down;//new feb2016
		//                j.corr(Form("original%s", "RelativeNonClosure_Pythia8High"), true) = jJESunc_up;//new feb2016
		//                j.corr(Form("original%s", "RelativeNonClosure_Pythia8Low"), true) = jJESunc_down;//new feb2016
		//             }
		// printf("JES uncertainty source %2d: %s\n",icomp,compName.Data());
  }	
      relUncertUp = tools->m_multiJES->getRelUncert(jetPt, jetEta, fCloseBy, true, e.npv(), e.mu(), (j.mom().Perp() > 20e3 && j.trueFlavour() == 5), 0, Nsegments);
      relUncertDown = tools->m_multiJES->getRelUncert(jetPt, jetEta, fCloseBy, false, e.npv(), e.mu(), (j.mom().Perp() > 20e3 && j.trueFlavour() == 5), 0, Nsegments);
    }

    if (relUncertUp == -1) {
      relUncertUp = 0;
    }
    if (relUncertDown == -1) {
      relUncertDown = 0;
    }
    pTJES_up *= 1 + relUncertUp;
    pTJES_down *= 1 - relUncertDown;

    // if no JES variation happens, this is the same as the unscaled one ...
    jJES_up.SetPtEtaPhiM(pTJES_up, j.mom().Eta(), j.mom().Phi(), j.mom().M());
    jJES_down.SetPtEtaPhiM(pTJES_down, j.mom().Eta(), j.mom().Phi(), j.mom().M());

    j.corr("jesUp", true) = jJES_up;
    j.corr("jesDown", true) = jJES_down;
    j.corr(Form("original%s", "jesUp"), true) = jJES_up;
    j.corr(Form("original%s", "jesDown"), true) = jJES_down;
  }
}

void Corrections::JER::run(Event &e) {
  if (e.isData()) return;
  tools->m_JetSmearingTool->SetSeed(tools->m_seedTool.getSeed(SeedTool::forJetEnergyResolution));
  for (int i = 0; i < e.jet().size(); ++i) {
    Jet &j = e.jet().at(i);

    TLorentzVector v = j.mom();

    if (e.atlFastII()) {
      tools->m_JetSmearingTool->SmearJet_Syst_AFII(v);
    } else {
      tools->m_JetSmearingTool->SmearJet_Syst(v);
    }
    j.corr("jer", true) = v;
    j.corr(Form("original%s", "jer"), true) = v;
  }
}

void Corrections::JEE::run(Event &e) {
  if (e.isData()) return;
  for (int i = 0; i < e.jet().size(); ++i) {
    Jet &j = e.jet().at(i);
    bool keepJet = tools->m_JEE.isGoodJet(j.mom().Perp()*1e-3, j.mom().Eta());
    TLorentzVector v = j.mom();
    if (!keepJet) {
      v.SetPtEtaPhiE(2e3, 0, 0, 2e3); // adds an artificial jet with 2 GeV: this is always rejected in the jet obj. def.
    }
    j.corr("jee", true) = v;
    j.corr(Form("original%s", "jee"), true) = v;
  }
}

void Corrections::METRecalc::run(Event &e) {
  if (!e.sr()) return;

  tools->m_metSys.reset(); // clear all information

  std::vector<float> el_shift, el_zeros;
  el_zeros.resize(e.sr()->el_n);
  for (int k = 0; k < e.sr()->el_n; ++k) {
    el_shift.push_back(e.electron()[k].mom().Perp()/e.sr()->el_pt->at(k) - 1.0);
  }

  std::vector<float> mu_shift, mu_shift_ms, mu_zeros;
  std::vector<float> mu_qoverp, mu_charge, mu_phi, mu_theta;
  mu_zeros.resize(e.sr()->mu_muid_n);

  for (int k = 0; k < e.sr()->mu_muid_n; ++k) {
    mu_theta.push_back(e.sr()->mu_muid_ms_theta->at(k));
    mu_phi.push_back(e.sr()->mu_muid_ms_phi->at(k));
    mu_qoverp.push_back(std::fabs(e.sr()->mu_muid_ms_qoverp->at(k)));
    mu_charge.push_back(e.muon()[k].charge());
    mu_shift.push_back(e.muon()[k].mom().Perp()/e.sr()->mu_muid_pt->at(k) - 1.0);
    float mspt = (e.muon()[k].momME().Perp() < 1e-19)? e.muon()[k].momMS().Perp() : e.muon()[k].momMS().Perp()*(e.muon()[k].momMECorr().Perp()/e.muon()[k].momME().Perp());
    mu_shift_ms.push_back(mspt/e.sr()->mu_muid_pt->at(k) - 1.0);
    if (mu_shift_ms[mu_shift_ms.size()-1] != mu_shift_ms[mu_shift_ms.size()-1])
      mu_shift_ms[mu_shift_ms.size()-1] = 0;

  }

  std::vector<float> jet_shift, jet_zeros;
  jet_zeros.resize(e.sr()->jet_AntiKt4LCTopo_n);
  for (int k = 0; k < e.sr()->jet_AntiKt4LCTopo_n; ++k) {
    jet_shift.push_back(e.jet()[k].corr("original").Perp()/e.sr()->jet_AntiKt4LCTopo_pt->at(k) - 1.0);
  }

  tools->m_metSys.setElectronParameters(e.sr()->el_pt, \
                                        e.sr()->el_eta, \
                                        e.sr()->el_phi, \
                                        e.sr()->el_MET_tightpp_wet, \
                                        e.sr()->el_MET_tightpp_wpx, \
                                        e.sr()->el_MET_tightpp_wpy, \
                                        e.sr()->el_MET_tightpp_statusWord);
  tools->m_metSys.setJetParameters(e.sr()->jet_AntiKt4LCTopo_pt, \
                                   e.sr()->jet_AntiKt4LCTopo_eta, \
                                   e.sr()->jet_AntiKt4LCTopo_phi, \
                                   e.sr()->jet_AntiKt4LCTopo_E, \
                                   e.sr()->jet_AntiKt4LCTopo_MET_tightpp_wet, \
                                   e.sr()->jet_AntiKt4LCTopo_MET_tightpp_wpx, \
                                   e.sr()->jet_AntiKt4LCTopo_MET_tightpp_wpy, \
                                   e.sr()->jet_AntiKt4LCTopo_MET_tightpp_statusWord);
  tools->m_metSys.setOriJetParameters(e.sr()->jet_AntiKt4LCTopo_pt);
  tools->m_metSys.setMuonParameters(e.sr()->mu_muid_pt, \
                                    e.sr()->mu_muid_eta, \
                                    e.sr()->mu_muid_phi, \
                                    e.sr()->mu_muid_MET_tightpp_wet, \
                                    e.sr()->mu_muid_MET_tightpp_wpx, \
                                    e.sr()->mu_muid_MET_tightpp_wpy, \
                                    e.sr()->mu_muid_MET_tightpp_statusWord);
  tools->m_metSys.setExtraMuonParameters(&mu_qoverp, &mu_theta, &mu_phi, &mu_charge);

  tools->m_metSys.setMETTerm(METUtil::CellOut, \
                             e.sr()->MET_CellOut_tightpp_etx, \
                             e.sr()->MET_CellOut_tightpp_ety, \
                             e.sr()->MET_CellOut_tightpp_sumet);
  tools->m_metSys.setMETTerm(METUtil::SoftJets, \
                             e.sr()->MET_SoftJets_etx, \
                             e.sr()->MET_SoftJets_ety, \
                             e.sr()->MET_SoftJets_sumet);
  tools->m_metSys.setMETTerm(METUtil::Truth, \
                             e.sr()->MET_Truth_NonInt_etx, \
                             e.sr()->MET_Truth_NonInt_ety, \
                             e.sr()->MET_Truth_NonInt_sumet);
  tools->m_metSys.setNvtx(e.npv_met());

  tools->m_metSys.setObjectEnergyUncertainties(METUtil::Electrons, \
                                               el_shift, el_zeros);
  tools->m_metSys.setObjectEnergyUncertainties(METUtil::Jets, \
                                               jet_shift, jet_zeros);
  tools->m_metSys.setObjectResolutionShift(METUtil::MuonsComboMS, \
                                           mu_shift, mu_zeros);
  tools->m_metSys.setObjectResolutionShift(METUtil::SpectroMuons, \
                                           mu_shift_ms, mu_zeros);

  METUtil::Systematics metSyst = METUtil::None;
  std::vector<std::string> allSystNames;
  std::vector<METUtil::Systematics> allSyst;
  allSystNames.push_back("");
  allSystNames.push_back("metResoSoftTermsUp");
  allSystNames.push_back("metResoSoftTermsDown");
  allSystNames.push_back("metScaleSoftTermsUp");
  allSystNames.push_back("metScaleSoftTermsDown");
  allSyst.push_back(METUtil::None);
  allSyst.push_back(METUtil::ResoSoftTermsUp);
  allSyst.push_back(METUtil::ResoSoftTermsDown);
  allSyst.push_back(METUtil::ScaleSoftTermsUp);
  allSyst.push_back(METUtil::ScaleSoftTermsDown);

  for (int x = 0; x < allSystNames.size(); ++x) {
    metSyst = allSyst[x];
    
    METObject METcell    = tools->m_metSys.getMissingET(METUtil::CellOut, metSyst);
    METObject METsoft    = tools->m_metSys.getMissingET(METUtil::SoftJets, metSyst);
    METObject METsoft_pt = tools->m_metSys.getMissingET(METUtil::SoftTerms, metSyst);

    METObject METele = tools->m_metSys.getMissingET(METUtil::RefEle, METUtil::EESUp);
    METObject METjet = tools->m_metSys.getMissingET(METUtil::RefJet, METUtil::JESUp);
    METObject METmuon= tools->m_metSys.getMissingET(METUtil::MuonTotal, METUtil::MERMSUp);

    float ex    = METele.etx()+METjet.etx()+METmuon.etx()+METcell.etx();
    float ey    = METele.ety()+METjet.ety()+METmuon.ety()+METcell.ety();
    float sumet = METele.sumet()+METjet.sumet()+METcell.sumet();
    e.metCorr(allSystNames[x], true).SetPxPyPzE(ex, ey, 0, std::sqrt(std::pow(ex, 2) + std::pow(ey, 2)));
  }
  e.met(e.metCorr("").Px(), e.metCorr("").Py());
}

