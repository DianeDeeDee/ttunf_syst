#include "RawReader.h"
#include "SkimReader.h"
#include <fastjet/PseudoJet.hh>
#include "Correction.h"
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/AreaDefinition.hh>
#include <fastjet/contrib/Nsubjettiness.hh>
#include <fastjet/contrib/Njettiness.hh>
#include <fastjet/contrib/NjettinessPlugin.hh>
#include <fastjet/tools/Filter.hh>

#include "HEPTopTagger.hh"

#include "egammaAnalysisUtils/IsEMPlusPlusDefs.h"

using namespace fastjet;
using namespace fastjet::contrib;


RawReader::RawReader(SkimReader &sr)
  : Reader(sr), m_doParticleLevelSelection(false) {
}

RawReader::~RawReader() {
}



Njettiness::AxesMode axisMode1 = Njettiness::onepass_kt_axes;
double beta1 = 1.0;
double R0 = 1.0;
double Rcut = 1.0;    



void RawReader::electron(std::vector<Electron> &v) {
  v.clear();
  for (int k = 0; k < m_sr->el_n; ++k) {
    v.push_back(Electron());
    Electron &e = v[k];
    e.mom().SetPtEtaPhiE(m_sr->el_pt->at(k),
                         m_sr->el_eta->at(k),
                         m_sr->el_phi->at(k),
                         m_sr->el_E->at(k));
    e.caloMom().SetPtEtaPhiE(m_sr->el_cl_pt->at(k),
                             m_sr->el_cl_eta->at(k),
                             m_sr->el_cl_phi->at(k),
                             m_sr->el_cl_E->at(k));
    e.trkMom().SetPtEtaPhiE(m_sr->el_cl_pt->at(k),
                             m_sr->el_tracketa->at(k), // CHANGE
                             m_sr->el_trackphi->at(k), // CHANGE
                             m_sr->el_cl_E->at(k));
    int tpp = m_sr->el_tightPP->at(k);
    
    double DEmaxs1 = (m_sr->el_emaxs1->at(k) - m_sr->el_Emax2->at(k))/(m_sr->el_emaxs1->at(k) + m_sr->el_Emax2->at(k));
    tpp = isTightPlusPlus(m_sr->el_etas2->at(k), m_sr->el_cl_E->at(k)/TMath::CosH(m_sr->el_etas2->at(k)),m_sr->el_f3->at(k),
                          m_sr->el_Ethad->at(k)/(m_sr->el_cl_E->at(k)/TMath::CosH(m_sr->el_etas2->at(k))), // rHad
                          m_sr->el_Ethad1->at(k)/(m_sr->el_cl_E->at(k)/TMath::CosH(m_sr->el_etas2->at(k))), // rHad1
                          m_sr->el_reta->at(k), m_sr->el_weta2->at(k),
                          m_sr->el_f1->at(k), m_sr->el_wstot->at(k),
                          DEmaxs1, //(m_sr->el_emaxs1->at(k) - m_sr->el_Emax2->at(k))/(m_sr->el_emaxs1->at(k) + m_sr->el_Emax2->at(k)),
                          m_sr->el_deltaeta1->at(k), m_sr->el_trackd0_physics->at(k),
                          m_sr->el_TRTHighTOutliersRatio->at(k), m_sr->el_nTRTHits->at(k), m_sr->el_nTRTOutliers->at(k),
                          m_sr->el_nSiHits->at(k), m_sr->el_nSCTOutliers->at(k) + m_sr->el_nPixelOutliers->at(k),
                          m_sr->el_nPixHits->at(k), m_sr->el_nPixelOutliers->at(k),
                          m_sr->el_nBLHits->at(k), m_sr->el_nBLayerOutliers->at(k), m_sr->el_expectHitInBLayer->at(k),
                          m_sr->el_cl_E->at(k) * fabs(m_sr->el_trackqoverp->at(k)), m_sr->el_deltaphi2->at(k),
                          //m_sr->el_isEM->at(k) & (1 << egammaPID::ConversionMatch_Electron),
                          m_sr->el_isEM->at(k) & (1 << 1),
                          egammaMenu::eg2012,false, false);
    e.setTightPP(tpp);
    int mpp = m_sr->el_mediumPP->at(k);
    mpp = isMediumPlusPlus(m_sr->el_etas2->at(k), m_sr->el_cl_E->at(k)/TMath::CosH(m_sr->el_etas2->at(k)),m_sr->el_f3->at(k),
                           m_sr->el_Ethad->at(k)/(m_sr->el_cl_E->at(k)/TMath::CosH(m_sr->el_etas2->at(k))), // rHad
                           m_sr->el_Ethad1->at(k)/(m_sr->el_cl_E->at(k)/TMath::CosH(m_sr->el_etas2->at(k))), // rHad1
                           m_sr->el_reta->at(k), m_sr->el_weta2->at(k),
                           m_sr->el_f1->at(k), m_sr->el_wstot->at(k),
                           (m_sr->el_emaxs1->at(k) - m_sr->el_Emax2->at(k))/(m_sr->el_emaxs1->at(k) + m_sr->el_Emax2->at(k)), // DEmaxs1
                           m_sr->el_deltaeta1->at(k), m_sr->el_trackd0_physics->at(k),
                           m_sr->el_TRTHighTOutliersRatio->at(k), m_sr->el_nTRTHits->at(k), m_sr->el_nTRTOutliers->at(k),
                           m_sr->el_nSiHits->at(k), m_sr->el_nSCTOutliers->at(k) + m_sr->el_nPixelOutliers->at(k),
                           m_sr->el_nPixHits->at(k), m_sr->el_nPixelOutliers->at(k),
                           m_sr->el_nBLHits->at(k), m_sr->el_nBLayerOutliers->at(k), m_sr->el_expectHitInBLayer->at(k),egammaMenu::eg2012,
                           false, false);
    e.setMediumPP(mpp);
    e.setLoosePP(m_sr->el_loosePP->at(k));
    e.passOR() = true;

    e.setMI(m_sr->el_MI10_max40_ptsum->at(k));
    //e.setMI(-1); // CHANGE
    e.z0() = m_sr->el_trackz0pvunbiased->at(k);
    //e.z0() = 0;//m_sr->el_z0->at(k);
    e.author() = m_sr->el_author->at(k);
    e.nSiHits() = m_sr->el_nSiHits->at(k);//+m_sr->el_nPixHits->at(k);
    e.oq() = m_sr->el_OQ->at(k);
    e.isEM() = m_sr->el_isEM->at(k);
    e.GSF_trk_index() = m_sr->el_GSF_trk_index->at(k);
  }
}

void RawReader::muon(std::vector<Muon> &v) {
  bool useMuid = true;
  v.clear();
  if (useMuid) {
    for (int k = 0; k < m_sr->mu_muid_n; ++k) {
      v.push_back(Muon());
      Muon &m = v[k];
      m.mom().SetPtEtaPhiE(m_sr->mu_muid_pt->at(k),
                           m_sr->mu_muid_eta->at(k),
                           m_sr->mu_muid_phi->at(k),
                           m_sr->mu_muid_E->at(k));
      m.setTight(m_sr->mu_muid_tight->at(k));
      m.setMI(m_sr->mu_muid_MI10_max40_ptsum->at(k));
      m.z0() = m_sr->mu_muid_trackz0pvunbiased->at(k);
      m.z0_exPV() = m_sr->mu_muid_id_z0_exPV->at(k);
      m.author() = m_sr->mu_muid_author->at(k);
      m.d0() = m_sr->mu_muid_trackd0pvunbiased->at(k);
      m.sd0() = m_sr->mu_muid_tracksigd0pvunbiased->at(k);
      m.author() = m_sr->mu_muid_author->at(k);
      m.passTrkCuts() = true;
      if (m_sr->mu_muid_nPixHits->at(k) + m_sr->mu_muid_nPixelDeadSensors->at(k) <= 0) m.passTrkCuts() = false;
      if (m_sr->mu_muid_nSCTHits->at(k) + m_sr->mu_muid_nSCTDeadSensors->at(k) <= 4) m.passTrkCuts() = false;
      if (m_sr->mu_muid_nPixHoles->at(k) + m_sr->mu_muid_nSCTHoles->at(k) >= 3) m.passTrkCuts() = false;
      double n = m_sr->mu_muid_nTRTHits->at(k) + m_sr->mu_muid_nTRTOutliers->at(k);
      if (std::fabs(m.mom().Eta()) > 0.1 && std::fabs(m.mom().Eta()) < 1.9) {
        if (! ((n > 5) && (((double) m_sr->mu_muid_nTRTOutliers->at(k))/n < 0.9)) )
          m.passTrkCuts() = false;
      }
      double pME  = (m_sr->mu_muid_me_qoverp->at(k) == 0) ? 0 : std::fabs(1.0/m_sr->mu_muid_me_qoverp->at(k));
      double etaME=-std::log(std::tan(m_sr->mu_muid_me_theta->at(k)/2.0));
      double ptME = (etaME != etaME) ? 1e-20 : pME/std::cosh(std::fabs(etaME));
      double phiME = m_sr->mu_muid_me_phi->at(k);
      m.momME().SetPtEtaPhiM(ptME, etaME, phiME, 0);
      m.momMECorr().SetPtEtaPhiM(ptME, etaME, phiME, 0);
  
      double pMS  = (m_sr->mu_muid_ms_qoverp->at(k) == 0) ? 0 : std::fabs(1.0/m_sr->mu_muid_ms_qoverp->at(k));
      double etaMS=-std::log(std::tan(m_sr->mu_muid_ms_theta->at(k)/2.0));
      double ptMS = (etaMS != etaMS) ? 1e-20 : pMS/std::cosh(std::fabs(etaMS));
      double phiMS = m_sr->mu_muid_ms_phi->at(k);
      m.momMS().SetPtEtaPhiM(ptMS, etaMS, phiMS, 0);
  
      double pID  =(m_sr->mu_muid_id_qoverp->at(k) == 0) ? 0 : std::fabs(1/m_sr->mu_muid_id_qoverp->at(k));
      double etaID=-std::log(std::tan(m_sr->mu_muid_id_theta->at(k)/2.0));
      double ptID =(etaID != etaID) ? 1e-20 : pID/std::cosh(std::fabs(etaID));
      double phiID = m_sr->mu_muid_id_phi->at(k);
      m.momID().SetPtEtaPhiM(ptID, etaID, phiID, 0);

      double pTrack  =(m_sr->mu_muid_trackqoverp->at(k) == 0) ? 0 : std::fabs(1.0/m_sr->mu_muid_trackqoverp->at(k));
      double etaTrack=-std::log(std::tan(m_sr->mu_muid_tracktheta->at(k)/2));
      double ptTrack =(etaTrack != etaTrack) ? 1e-20 : pTrack/cosh(fabs(etaTrack));
      double phiTrack = m_sr->mu_muid_trackphi->at(k);
      m.momTrk().SetPtEtaPhiM(ptTrack, etaTrack, phiTrack, 0);

      m.charge() = m_sr->mu_muid_charge->at(k);
      m.setST(m_sr->mu_muid_isSegmentTaggedMuon->at(k));
      m.setSA(m_sr->mu_muid_isStandAloneMuon->at(k));
      m.setCB(m_sr->mu_muid_isCombinedMuon->at(k));
    }
  } else {
    /*
    for (int k = 0; k < m_sr->mu_n; ++k) {
      v.push_back(Muon());
      Muon &m = v[k];
      m.mom().SetPtEtaPhiE(m_sr->mu_pt->at(k),
                           m_sr->mu_eta->at(k),
                           m_sr->mu_phi->at(k),
                           m_sr->mu_E->at(k));
      m.setTight(m_sr->mu_tight->at(k));
      m.setMI(m_sr->mu_MI10_max40_ptsum->at(k));
      m.z0() = m_sr->mu_trackz0pvunbiased->at(k);
      m.z0_exPV() = m_sr->mu_id_z0_exPV->at(k);
      m.author() = m_sr->mu_author->at(k);
      m.d0() = m_sr->mu_trackd0pvunbiased->at(k);
      m.sd0() = m_sr->mu_tracksigd0pvunbiased->at(k);
      m.author() = m_sr->mu_author->at(k);
      m.passTrkCuts() = true;
      if (m_sr->mu_nPixHits->at(k) + m_sr->mu_nPixelDeadSensors->at(k) <= 0) m.passTrkCuts() = false;
      if (m_sr->mu_nSCTHits->at(k) + m_sr->mu_nSCTDeadSensors->at(k) <= 4) m.passTrkCuts() = false;
      if (m_sr->mu_nPixHoles->at(k) + m_sr->mu_nSCTHoles->at(k) >= 3) m.passTrkCuts() = false;
      double n = m_sr->mu_nTRTHits->at(k) + m_sr->mu_nTRTOutliers->at(k);
      if (std::fabs(m.mom().Eta()) > 0.1 && std::fabs(m.mom().Eta()) < 1.9) {
        if (! ((n > 5) && (((double) m_sr->mu_nTRTOutliers->at(k))/n < 0.9)) )
          m.passTrkCuts() = false;
      }
      double pME  = (m_sr->mu_me_qoverp->at(k) == 0) ? 0 : std::fabs(1.0/m_sr->mu_me_qoverp->at(k));
      double etaME=-std::log(std::tan(m_sr->mu_me_theta->at(k)/2.0));
      double ptME = (etaME != etaME) ? 1e-20 : pME/std::cosh(std::fabs(etaME));
      double phiME = m_sr->mu_me_phi->at(k);
      m.momME().SetPtEtaPhiM(ptME, etaME, phiME, 0);
  
      double pMS  = (m_sr->mu_ms_qoverp->at(k) == 0) ? 0 : std::fabs(1.0/m_sr->mu_ms_qoverp->at(k));
      double etaMS=-std::log(std::tan(m_sr->mu_ms_theta->at(k)/2.0));
      double ptMS = (etaMS != etaMS) ? 1e-20 : pMS/std::cosh(std::fabs(etaMS));
      double phiMS = m_sr->mu_ms_phi->at(k);
      m.momMS().SetPtEtaPhiM(ptMS, etaMS, phiMS, 0);
  
      double pID  =(m_sr->mu_id_qoverp->at(k) == 0) ? 0 : std::fabs(1/m_sr->mu_id_qoverp->at(k));
      double etaID=-std::log(std::tan(m_sr->mu_id_theta->at(k)/2.0));
      double ptID =(etaID != etaID) ? 1e-20 : pID/std::cosh(std::fabs(etaID));
      double phiID = m_sr->mu_id_phi->at(k);
      m.momID().SetPtEtaPhiM(ptID, etaID, phiID, 0);

      double pTrack  =(m_sr->mu_trackqoverp->at(k) == 0) ? 0 : std::fabs(1.0/m_sr->mu_trackqoverp->at(k));
      double etaTrack=-std::log(std::tan(m_sr->mu_tracktheta->at(k)/2));
      double ptTrack =(etaTrack != etaTrack) ? 1e-20 : pTrack/cosh(fabs(etaTrack));
      double phiTrack = m_sr->mu_trackphi->at(k);
      m.momTrk().SetPtEtaPhiM(ptTrack, etaTrack, phiTrack, 0);

      m.charge() = m_sr->mu_charge->at(k);
      m.setST(m_sr->mu_isSegmentTaggedMuon->at(k));
      m.setSA(m_sr->mu_isStandAloneMuon->at(k));
      m.setCB(m_sr->mu_isCombinedMuon->at(k));
    }
    */
  }
}

void RawReader::met(float &met_x, float &met_y) {
  met_x = m_sr->MET_RefFinal_tightpp_etx;
  met_y = m_sr->MET_RefFinal_tightpp_ety;
}

void RawReader::trigger(bool &triggerElectron, bool &triggerMuon, bool &triggerLargeJet){// bool &triggerElectron24, bool &triggerElectron60) {
  triggerElectron = m_sr->EF_e24vhi_medium1 || m_sr->EF_e60_medium1;
  triggerMuon = m_sr->EF_mu24i_tight || m_sr->EF_mu36_tight;
  triggerLargeJet = m_sr -> EF_j360_a4tchad;
 // triggerElectron24 = m_sr->EF_e24vhi_medium1;
 // triggerElectron60 = m_sr->EF_e60_medium1;
}




void RawReader::jet(std::vector<Jet> &v, std::vector<LargeJet> &u){//, std::vector<LargeJet> &w) {
  double minLargeJetPt = 150e3;
  //double minJetPt = 25e3;
  //double minSubjetPt = 20e3;
  double minSubjetPt = 2e3;
  v.clear();
  u.clear();
  //w.clear();


  std::vector<PseudoJet> cluster;
  for (int k = 0; k < m_sr->cl_lc_n; ++k) {
    double px = m_sr->cl_lc_pt->at(k)*std::cos(m_sr->cl_lc_phi->at(k));
    double py = m_sr->cl_lc_pt->at(k)*std::sin(m_sr->cl_lc_phi->at(k));
    double pz = m_sr->cl_lc_pt->at(k)*std::sinh(m_sr->cl_lc_eta->at(k));
    double E  = std::sqrt(std::pow(px, 2) + std::pow(py, 2) + std::pow(pz, 2));
    cluster.push_back(PseudoJet(px, py, pz, E));
  }
  
  JetDefinition jet_def_large = JetDefinition(antikt_algorithm, 1.0, fastjet::E_scheme, fastjet::Best);
  JetDefinition jet_def_largeBB;
  //if (Correction::globalTools->m_useCA12)
    //jet_def_largeBB = JetDefinition(cambridge_algorithm, 1.2, fastjet::E_scheme, fastjet::Best);
  //else
    jet_def_largeBB = JetDefinition(cambridge_algorithm, 1.0, fastjet::E_scheme, fastjet::Best);

  JetDefinition jet_def_small;
 // if (Correction::globalTools->m_useR03)
   // jet_def_small = JetDefinition(cambridge_algorithm, 0.3, fastjet::E_scheme, fastjet::Best);
  //else
    jet_def_small = JetDefinition(cambridge_algorithm, 0.2, fastjet::E_scheme, fastjet::Best);



  for (int k = 0; k < m_sr->jet_AntiKt4LCTopo_n; ++k) {
    v.push_back(Jet());
    v[k].mom().SetPtEtaPhiE(m_sr->jet_AntiKt4LCTopo_pt->at(k),
                            m_sr->jet_AntiKt4LCTopo_eta->at(k),
                            m_sr->jet_AntiKt4LCTopo_phi->at(k),
                            m_sr->jet_AntiKt4LCTopo_E->at(k));
    v[k].corr("original", true) = v[k].mom();
    v[k].corr("originalnom", true) = v[k].mom();
    v[k].corr("originaljesUp", true) = v[k].mom();
    v[k].corr("originaljesDown", true) = v[k].mom();
    v[k].corr("originaljer", true) = v[k].mom();
    v[k].corr("originaljee", true) = v[k].mom();
    if (!Correction::globalTools->m_isData) {
      v[k].trueFlavour() = m_sr->jet_AntiKt4LCTopo_flavor_truth_label->at(k);
    } else {
      v[k].trueFlavour() = -999;
    }
   
    v[k].mv1() = m_sr->jet_AntiKt4LCTopo_flavor_weight_MV1->at(k);
    v[k].jvf() = m_sr->jet_AntiKt4LCTopo_jvtxf->at(k);
    v[k].detE() = m_sr->jet_AntiKt4LCTopo_constscale_E->at(k);
    v[k].detEta() = m_sr->jet_AntiKt4LCTopo_constscale_eta->at(k);
    v[k].detPhi() = m_sr->jet_AntiKt4LCTopo_constscale_phi->at(k);
    v[k].detM() = m_sr->jet_AntiKt4LCTopo_constscale_m->at(k);
    v[k].Ax() = m_sr->jet_AntiKt4LCTopo_ActiveAreaPx->at(k);
    v[k].Ay() = m_sr->jet_AntiKt4LCTopo_ActiveAreaPy->at(k);
    v[k].Az() = m_sr->jet_AntiKt4LCTopo_ActiveAreaPz->at(k);
    v[k].Ae() = m_sr->jet_AntiKt4LCTopo_ActiveAreaE->at(k);
    v[k].isBadLoose() = m_sr->jet_AntiKt4LCTopo_isBadLoose->at(k);
    v[k].isBadLooseMinus() = m_sr->jet_AntiKt4LCTopo_isBadLooseMinus->at(k);
    v[k].TrackAssoc_index() = m_sr->jet_AntiKt4LCTopo_TrackAssoc_index->at(k);
    v[k].sumPtTrk() = m_sr->jet_AntiKt4LCTopo_sumPtTrk_pv0_500MeV->at(k);
    v[k].jfit_deltaEta() = m_sr->jet_AntiKt4LCTopo_flavor_component_jfit_deltaEta->at(k);
    v[k].jfit_deltaPhi() = m_sr->jet_AntiKt4LCTopo_flavor_component_jfit_deltaPhi->at(k);
  }




  fastjet::AreaDefinition areaDefinition(active_area);

  for (int k = 0; k < m_sr->jet_AntiKt10LCTopoTrimmedPtFrac5SmallR30_n; ++k) {
    u.push_back(LargeJet());
    u[k].mom().SetPtEtaPhiE(m_sr->jet_AntiKt10LCTopoTrimmedPtFrac5SmallR30_pt->at(k), \
                            m_sr->jet_AntiKt10LCTopoTrimmedPtFrac5SmallR30_eta->at(k), \
                            m_sr->jet_AntiKt10LCTopoTrimmedPtFrac5SmallR30_phi->at(k), \
                            m_sr->jet_AntiKt10LCTopoTrimmedPtFrac5SmallR30_E->at(k));
    u[k].detE() = m_sr->jet_AntiKt10LCTopoTrimmedPtFrac5SmallR30_constscale_E->at(k);
    u[k].detEta() = m_sr->jet_AntiKt10LCTopoTrimmedPtFrac5SmallR30_constscale_eta->at(k);
    u[k].detPhi() = m_sr->jet_AntiKt10LCTopoTrimmedPtFrac5SmallR30_constscale_phi->at(k);
    u[k].detM() = m_sr->jet_AntiKt10LCTopoTrimmedPtFrac5SmallR30_constscale_m->at(k);
    u[k].Ax() = 0;
    u[k].Ay() = 0;
    u[k].Az() = 0;
    u[k].Ae() = 0;
    u[k].split12() = m_sr->jet_AntiKt10LCTopoTrimmedPtFrac5SmallR30_SPLIT12->at(k);
    u[k].trueFlavour() = 0;

    u[k].trimmed_tau1() = m_sr->jet_AntiKt10LCTopoTrimmedPtFrac5SmallR30_Tau1->at(k);
    u[k].trimmed_tau2() = m_sr->jet_AntiKt10LCTopoTrimmedPtFrac5SmallR30_Tau2->at(k);

    if (!Correction::globalTools->m_isData) {
      for (int z = 0; z < m_sr->mc_n; ++z) {
        if (m_sr->mc_pt->at(z) < 50e3) continue;
        int pdg = (int) std::fabs(m_sr->mc_pdgId->at(z));
        TLorentzVector v;
        float pt = m_sr->mc_pt->at(z);
        float eta = m_sr->mc_eta->at(z);
        float phi = m_sr->mc_phi->at(z);
        float m = m_sr->mc_m->at(z);
        v.SetPtEtaPhiM(pt, eta, phi, m);
        if (v.DeltaR(u[k].mom()) > 0.3) continue;
        if (pdg == 25) {
          u[k].trueFlavour() = 25;
          break;
        }
        if (pdg == 6) {
          u[k].trueFlavour() = 6;
          break;
        }
        if (pdg == 5) {
          u[k].trueFlavour() = 5;
          break;
        }
      }
    }

    std::vector<PseudoJet> cl;
    int pind = m_sr->jet_AntiKt10LCTopoTrimmedPtFrac5SmallR30_Parent_index->at(k)[0];
    for (int z = 0; z < m_sr->jet_AntiKt10LCTopo_constit_n->at(pind); ++z) {
      int idx = m_sr->jet_AntiKt10LCTopo_constit_index->at(pind).at(z);
      cl.push_back(fastjet::PseudoJet());
      TLorentzVector v;
      v.SetPtEtaPhiM(m_sr->cl_lc_pt->at(idx), \
                     m_sr->cl_lc_eta->at(idx), \
                     m_sr->cl_lc_phi->at(idx), \
                     0);
      cl[z].reset_momentum(v.Px(), v.Py(), v.Pz(), v.E());
    }

    for (int z = 0; z < cl.size(); ++z) cl[z] *= 1e-3;
    fastjet::ClusterSequence fj_cl(cl, jet_def_large);
    vector<PseudoJet> fat_jet =  sorted_by_pt(fj_cl.inclusive_jets());
    u[k].htt() = -1;
    if (fat_jet.size() > 0) {
      HEPTopTagger htt(fj_cl, fat_jet[0], 172.3, 80.3);
      htt.set_top_range(172.3 - 25.0, 172.3 + 25.0);
      htt.run_tagger();
      if (!htt.is_masscut_passed()) {
        u[k].htt() = -1;
      } else {
        u[k].htt() = htt.top_candidate().m();
      }
    }
    for (int z = 0; z < cl.size(); ++z) cl[z] *= 1e3;

    if (m_readBB) {
      fastjet::ClusterSequence fj_cl(cl, jet_def_large);
      vector<PseudoJet> fat_jet =  sorted_by_pt(fj_cl.inclusive_jets());
     /*
      Njettiness::AxesMode axisMode1 = Njettiness::onepass_kt_axes;
      double beta1 = 1.0;
      double R0 = 1.0;
      double Rcut = 1.0;
     */
      Nsubjettiness nSub1(1, axisMode1, beta1, R0, Rcut);
      //double tau1 = nSub1(fat_jet[0]);
      Nsubjettiness nSub2(2, axisMode1, beta1, R0, Rcut);
      // double tau2 = nSub2(fat_jet[0]);
      //Nsubjettiness nSub3(3, axisMode1, beta1, R0, Rcut);
      // double tau3 = nSub3(fat_jet[0]);

      u[k].trimmed_am_tau1() =  nSub1(fat_jet[0]);
      u[k].trimmed_am_tau2() =  nSub2(fat_jet[0]);
    }

    fastjet::ClusterSequenceArea cs_cl(cl, jet_def_small, areaDefinition);
    std::vector<fastjet::PseudoJet> subjets = sorted_by_pt(cs_cl.inclusive_jets(minSubjetPt));
    u[k].subjet().clear();
    u[k].subjetJER().clear();
    u[k].subjet_area().clear();
    for (int j = 0; j < subjets.size(); ++j) {
      //u[k].subjet().push_back(subjets[j]);
      //u[k].subjet_area().push_back(subjets[j].area_4vector());
      /*fastjet::PseudoJet a = subjets[j].area_4vector();
      u[k].subjet().push_back(LargeJet::Subjet());
      u[k].subjet()[j].v.SetPxPyPzE(subjets[j].px(), subjets[j].py(), subjets[j].pz(), subjets[j].e());
      u[k].subjet()[j].lcpt = subjets[j].perp();
      u[k].subjetJER().push_back(LargeJet::Subjet());
      u[k].subjetJER()[j].v.SetPxPyPzE(subjets[j].px(), subjets[j].py(), subjets[j].pz(), subjets[j].e());
      u[k].subjetJER()[j].lcpt = subjets[j].perp();
      u[k].subjet_area().push_back(TLorentzVector());
      u[k].subjet_area()[j].SetPxPyPzE(a.px(), a.py(), a.pz(), a.e());
    */}
  }

  if (m_readBB) {
    ClusterSequence cs_largeBB(cluster, jet_def_largeBB);
    //std::vector<PseudoJet> ljetsBB = sorted_by_pt(cs_largeBB.inclusive_jets(minLargeJetPt));
    /*for (int k = 0; k < ljetsBB.size(); ++k) {
      w.push_back(LargeJet());
      //fastjet::PseudoJet alj = ljetsBB[k].area_4vector();
      w[k].mom().SetPxPyPzE(ljetsBB[k].px(), ljetsBB[k].py(), ljetsBB[k].pz(), ljetsBB[k].e());
      w[k].detE() = w[k].mom().E();
      w[k].detEta() = w[k].mom().Eta();
      w[k].detPhi() = w[k].mom().Phi();
      w[k].detM() = w[k].mom().M();
      w[k].Ax() = 0; //alj.px();
      w[k].Ay() = 0; //alj.py();
      w[k].Az() = 0; //alj.pz();
      w[k].Ae() = 0; //alj.e();
      w[k].split12() = 0;
      w[k].trueFlavour() = 0;
  
      std::vector<PseudoJet> cl = cs_largeBB.constituents(ljetsBB[k]);
      fastjet::ClusterSequenceArea cs_cl(cl, jet_def_small, areaDefinition);
      std::vector<fastjet::PseudoJet> subjets = sorted_by_pt(cs_cl.inclusive_jets(minSubjetPt));
      w[k].subjet().clear();
      w[k].subjetJER().clear();
      w[k].subjet_area().clear();
      for (int j = 0; j < subjets.size(); ++j) {
        fastjet::PseudoJet a = subjets[j].area_4vector();
        w[k].subjet().push_back(LargeJet::Subjet());
        w[k].subjet()[j].v.SetPxPyPzE(subjets[j].px(), subjets[j].py(), subjets[j].pz(), subjets[j].e());
        w[k].subjet()[j].lcpt = subjets[j].perp();
        w[k].subjetJER().push_back(LargeJet::Subjet());
        w[k].subjetJER()[j].v.SetPxPyPzE(subjets[j].px(), subjets[j].py(), subjets[j].pz(), subjets[j].e());
        w[k].subjetJER()[j].lcpt = subjets[j].perp();

        //      std::cout << " Subjet pT: " << subjets[j].pt() << " Eta: " <<  subjets[j].eta() << std:: endl;

        w[k].subjet_area().push_back(TLorentzVector());
        w[k].subjet_area()[j].SetPxPyPzE(a.px(), a.py(), a.pz(), a.e());
      }

      Nsubjettiness nSub1_noam(1, Njettiness::kt_axes, beta1, R0, Rcut);
      Nsubjettiness nSub2_noam(2, Njettiness::kt_axes, beta1, R0, Rcut);
      w[k].ug_tau1() =  nSub1_noam(ljetsBB[k]);
      w[k].ug_tau2() =  nSub2_noam(ljetsBB[k]);
 
      Nsubjettiness nSub1(1, axisMode1, beta1, R0, Rcut);
      Nsubjettiness nSub2(2, axisMode1, beta1, R0, Rcut);
      w[k].ug_am_tau1() =  nSub1(ljetsBB[k]);
      w[k].ug_am_tau2() =  nSub2(ljetsBB[k]);

      w[k].trimmed_tau1() = -1;
      w[k].trimmed_tau2() = -1;
      w[k].trimmed_am_tau1() = -1;
      w[k].trimmed_am_tau2() = -1;
      w[k].bdrs() = w[k].mom();
    }*/
  } //else if (Correction::globalTools->m_useCA12) {
    for (int k = 0; k < m_sr->jet_CamKt12LCTopo_n; ++k) {
  /*    w.push_back(LargeJet());
      w[k].mom().SetPtEtaPhiE(m_sr->jet_CamKt12LCTopo_pt->at(k), \
                            m_sr->jet_CamKt12LCTopo_eta->at(k), \
                            m_sr->jet_CamKt12LCTopo_phi->at(k), \
                            m_sr->jet_CamKt12LCTopo_E->at(k));
      w[k].detE() = m_sr->jet_CamKt12LCTopo_constscale_E->at(k);
      w[k].detEta() = m_sr->jet_CamKt12LCTopo_constscale_eta->at(k);
      w[k].detPhi() = m_sr->jet_CamKt12LCTopo_constscale_phi->at(k);
      w[k].detM() = m_sr->jet_CamKt12LCTopo_constscale_m->at(k);
      w[k].Ax() = 0;
      w[k].Ay() = 0;
      w[k].Az() = 0;
      w[k].Ae() = 0;
      w[k].split12() = m_sr->jet_CamKt12LCTopo_SPLIT12->at(k);
      w[k].trueFlavour() = 0;
*/
      std::vector<PseudoJet> cl;
      for (int z = 0; z < m_sr->jet_CamKt12LCTopo_constit_n->at(k); ++z) {
        int idx = m_sr->jet_CamKt12LCTopo_constit_index->at(k).at(z);
        cl.push_back(fastjet::PseudoJet());
        TLorentzVector v;
        v.SetPtEtaPhiM(m_sr->cl_lc_pt->at(idx), \
                       m_sr->cl_lc_eta->at(idx), \
                       m_sr->cl_lc_phi->at(idx), \
                       0);
        cl[z].reset_momentum(v.Px(), v.Py(), v.Pz(), v.E());
      }

      fastjet::ClusterSequenceArea cs_cl(cl, jet_def_small, areaDefinition);
      std::vector<fastjet::PseudoJet> subjets = sorted_by_pt(cs_cl.inclusive_jets(minSubjetPt));
    /*  w[k].subjet().clear();
      w[k].subjetJER().clear();
      w[k].subjet_area().clear();
      */
      for (int j = 0; j < subjets.size(); ++j) {
        /*fastjet::PseudoJet a = subjets[j].area_4vector();
        w[k].subjet().push_back(LargeJet::Subjet());
        w[k].subjet()[j].v.SetPxPyPzE(subjets[j].px(), subjets[j].py(), subjets[j].pz(), subjets[j].e());
        w[k].subjet()[j].lcpt = subjets[j].perp();
        w[k].subjetJER().push_back(LargeJet::Subjet());
        w[k].subjetJER()[j].v.SetPxPyPzE(subjets[j].px(), subjets[j].py(), subjets[j].pz(), subjets[j].e());
        w[k].subjetJER()[j].lcpt = subjets[j].perp();
        w[k].subjet_area().push_back(TLorentzVector());
        w[k].subjet_area()[j].SetPxPyPzE(a.px(), a.py(), a.pz(), a.e());
      */}

      //w[k].htt() = -1;
      for (int z = 0; z < cl.size(); ++z) cl[z] *= 1e-3;
      fastjet::ClusterSequence fj_cl(cl, jet_def_large);
      vector<PseudoJet> fat_jet =  sorted_by_pt(fj_cl.inclusive_jets());
      if (fat_jet.size() > 0) {
        HEPTopTagger htt(fj_cl, fat_jet[0], 172.3, 80.3);
        htt.set_top_range(172.3 - 25.0, 172.3 + 25.0);
        htt.run_tagger();
        if (!htt.is_masscut_passed()) {
          //w[k].htt() = -1;
        } else {
        //  w[k].htt() = htt.top_candidate().m();
        }
      }
      for (int z = 0; z < cl.size(); ++z) cl[z] *= 1e3;
    //}
  }
} 

void RawReader::readBB(const bool readBB) {
  m_readBB = readBB;
}

void RawReader::partElectron(std::vector<Particle> &v) {
  v.clear();
  if (!m_doParticleLevelSelection) return;
  for (int k = 0; k < m_sr->mc_n; ++k) {
    if (m_sr->mc_pt->at(k) == 0) continue;
    if (std::fabs(m_sr->mc_pdgId->at(k)) != 11) continue;
    if (m_sr->mc_status->at(k) != 1 || m_sr->mc_barcode->at(k) >= 200000) continue;
    v.push_back(Particle());
    Particle &e = v[v.size()-1];
    float pt = m_sr->mc_pt->at(k);
    float eta = m_sr->mc_eta->at(k);
    float phi = m_sr->mc_phi->at(k);
    float m = m_sr->mc_m->at(k);
    e.mom().SetPtEtaPhiM(pt, eta, phi, m);
    float lpt = e.mom().Perp();
    float drCut = 10e3/lpt; // mini-isolation
    float isol = 0;
    for (int l = 0; l < m_sr->mc_n; ++l) {
      if (l == k) continue;
      if (m_sr->mc_status->at(l) != 1 || m_sr->mc_barcode->at(l) >= 200000) continue;
      int pdg = std::fabs(m_sr->mc_pdgId->at(l));
      if (pdg == 12 || pdg == 14 || pdg == 16) continue;
      TLorentzVector x;
      x.SetPtEtaPhiM(m_sr->mc_pt->at(l),
                     m_sr->mc_eta->at(l),
                     m_sr->mc_phi->at(l),
                     m_sr->mc_m->at(l));
      float dr = x.DeltaR(e.mom());
      if (dr < 0.1) { // photon
        e.mom() += x;
      } else if (dr < drCut && x.Perp() > 500) {
        isol += x.Perp();
      }
    }
    e.setIsol(isol);
  }
}

void RawReader::partMom(std::vector<Particle> &v) {
  v.clear();
  int pdgs[] = {25, 6, 5, 24};
  std::vector<int> done;
  for (int k = 0; k < m_sr->mc_n; ++k) {
    if (m_sr->mc_pt->at(k) == 0) continue;
    if (std::find(done.begin(), done.end(), k) != done.end()) continue;
    int pdg = (int) std::fabs(m_sr->mc_pdgId->at(k));
    bool isTop = pdg == 6;
    bool isLep = false;
    for (int l = 0; l < sizeof(pdgs)/sizeof(int); ++l) {
      if (pdg == pdgs[l]) {
        int child = k;
        bool last = true;
        do {
          last = true;
          for (int m = 0; m < m_sr->mc_child_index->at(child).size(); ++m) {
            int child_idx = (m_sr->mc_child_index->at(child).at(m));
            int child_pdg = (int) std::fabs(m_sr->mc_pdgId->at(child_idx));
            if (child_pdg == pdg) {
              child = child_idx;
              last = false;
              break;
            }
            if (isTop && child_pdg == 24) {
              int wchild = child_idx;
              bool wlast = true;
              do {
                wlast = true;
                for (int n = 0; n < m_sr->mc_child_index->at(wchild).size(); ++n) {
                  int wchild_idx = (m_sr->mc_child_index->at(wchild).at(n));
                  int wchild_pdg = (int) std::fabs(m_sr->mc_pdgId->at(wchild_idx));
                  if (wchild_pdg == 24) {
                    wchild = wchild_idx;
                    wlast = false;
                    break;
                  }
                  if ((wchild_pdg == 11) || (wchild_pdg == 13) || (wchild_pdg == 15)) {
                    isLep = true;
                  }
                }
              } while (!wlast);
            }
          }
        } while (!last);
        if (std::find(done.begin(), done.end(), child) != done.end()) continue;
        done.push_back(child);
        v.push_back(Particle());
        v[v.size()-1].mom().SetPtEtaPhiM(m_sr->mc_pt->at(child),
                                         m_sr->mc_eta->at(child),
                                         m_sr->mc_phi->at(child),
                                         m_sr->mc_m->at(child));
        v[v.size()-1].setId(m_sr->mc_pdgId->at(child));
        if (isTop && isLep) {
          v.push_back(Particle());
          v[v.size()-1].mom().SetPtEtaPhiM(m_sr->mc_pt->at(child),
                                           m_sr->mc_eta->at(child),
                                           m_sr->mc_phi->at(child),
                                           m_sr->mc_m->at(child));
          v[v.size()-1].setId(m_sr->mc_pdgId->at(child)*100000);
        } else if (isTop && !isLep) {
          v.push_back(Particle());
          v[v.size()-1].mom().SetPtEtaPhiM(m_sr->mc_pt->at(child),
                                           m_sr->mc_eta->at(child),
                                           m_sr->mc_phi->at(child),
                                           m_sr->mc_m->at(child));
          v[v.size()-1].setId(m_sr->mc_pdgId->at(child)*1000000);
        }
        break;
      }
    }
  }
}

void RawReader::partMuon(std::vector<Particle> &v) {
  v.clear();
  if (!m_doParticleLevelSelection) return;
  for (int k = 0; k < m_sr->mc_n; ++k) {
    if (m_sr->mc_pt->at(k) == 0) continue;
    if (std::fabs(m_sr->mc_pdgId->at(k)) != 13) continue;
    if (m_sr->mc_status->at(k) != 1 || m_sr->mc_barcode->at(k) >= 200000) continue;
    v.push_back(Particle());
    Particle &m = v[v.size()-1];
    m.mom().SetPtEtaPhiM(m_sr->mc_pt->at(k),
                         m_sr->mc_eta->at(k),
                         m_sr->mc_phi->at(k),
                         m_sr->mc_m->at(k));
    float lpt = m.mom().Perp();
    float drCut = 10e3/lpt; // mini-isolation
    float isol = 0;
    for (int l = 0; l < m_sr->mc_n; ++l) {
      if (l == k) continue;
      if (m_sr->mc_status->at(l) != 1 || m_sr->mc_barcode->at(l) >= 200000) continue;
      int pdg = std::fabs(m_sr->mc_pdgId->at(l));
      if (pdg == 12 || pdg == 14 || pdg == 16) continue;
      TLorentzVector x;
      x.SetPtEtaPhiM(m_sr->mc_pt->at(l),
                     m_sr->mc_eta->at(l),
                     m_sr->mc_phi->at(l),
                     m_sr->mc_m->at(l));
      float dr = x.DeltaR(m.mom());
      if (dr < 0.3 && x.Perp() > 500) {
        isol += x.Perp();
      }
    }
    m.setIsol(isol);
  }
}

void RawReader::partMet(float &met_x, float &met_y) {
  met_x = m_sr->MET_Truth_NonInt_etx;
  met_y = m_sr->MET_Truth_NonInt_ety;
}


 void RawReader::partJet(std::vector<Jet> &v, std::vector<LargeJet> &u){//, std::vector<LargeJet> &w) 
  if (!m_doParticleLevelSelection) return;
  double minLargeJetPt = 200e3;
  //double minJetPt = 25e3;
  double minSubjetPt = 5e3;
  v.clear();
  u.clear();
 // w.clear();

  std::vector<PseudoJet> cluster;
  for (int k = 0; k < m_sr->mc_n; ++k) {
    if (m_sr->mc_pt->at(k) == 0) continue;
    if (m_sr->mc_status->at(k) != 1 || m_sr->mc_barcode->at(k) >= 200000) continue;
    //if (m_sr->mc_child_index->at(k).size() != 0) continue;
    int pdg = std::fabs(m_sr->mc_pdgId->at(k));
    if (pdg == 12 || pdg == 14 || pdg == 16) continue;
    //if (std::fabs(m_sr->mc_charge->at(k)) < 1e-4) continue;

    TLorentzVector tmp;
    tmp.SetPtEtaPhiM(m_sr->mc_pt->at(k),
                     m_sr->mc_eta->at(k),
                     m_sr->mc_phi->at(k),
                     m_sr->mc_m->at(k));

    double px = tmp.Px();
    double py = tmp.Py();
    double pz = tmp.Pz();
    double E  = tmp.E();
    if (E <= 0) continue;
    cluster.push_back(PseudoJet(px, py, pz, E));
  }
  JetDefinition jet_def_large = JetDefinition(antikt_algorithm, 1.0);
  JetDefinition jet_def_largeBB = JetDefinition(cambridge_algorithm, 1.2);
  JetDefinition jet_def_small;
  //if (Correction::globalTools->m_useR03)
    //jet_def_small = JetDefinition(cambridge_algorithm, 0.3);
  //else
    jet_def_small = JetDefinition(cambridge_algorithm, 0.2);

  for (int k = 0; k < m_sr->AntiKt4Truth_n; ++k) {
    v.push_back(Jet());
    v[k].mom().SetPtEtaPhiE(m_sr->AntiKt4Truth_pt->at(k),
                            m_sr->AntiKt4Truth_eta->at(k),
                            m_sr->AntiKt4Truth_phi->at(k),
                            m_sr->AntiKt4Truth_E->at(k));
    v[k].trueFlavour() = m_sr->AntiKt4Truth_flavor_truth_label->at(k);
    v[k].mv1() = 0;
    v[k].jvf() = 1;
    if (v[k].trueFlavour() == 5)
      v[k].mv1() = 99;
  }

  fastjet::AreaDefinition areaDefinition(active_area);
  ClusterSequence cs_large(cluster, jet_def_large);
  std::vector<PseudoJet> ljets = sorted_by_pt(cs_large.inclusive_jets(minLargeJetPt));
  Filter trimmer(0.3, SelectorPtFractionMin(0.05));
  for (int k = 0; k < ljets.size(); ++k) {
    PseudoJet tjet = trimmer(ljets[k]);
    u.push_back(LargeJet());
    u[k].mom().SetPxPyPzE(tjet.px(), tjet.py(), tjet.pz(), tjet.e());
    u[k].split12() = -1;
    u[k].trueFlavour() = 0;
    for (int z = 0; z < m_sr->mc_n; ++z) {
      if (m_sr->mc_pt->at(z) < 50e3) continue;
      int pdg = (int) std::fabs(m_sr->mc_pdgId->at(z));
      TLorentzVector v;
      float pt = m_sr->mc_pt->at(z);
      float eta = m_sr->mc_eta->at(z);
      float phi = m_sr->mc_phi->at(z);
      float m = m_sr->mc_m->at(z);
      v.SetPtEtaPhiM(pt, eta, phi, m);
      if (v.DeltaR(u[k].mom()) > 0.3) continue;
      if (pdg == 25) {
        u[k].trueFlavour() = 25;
        break;
      }
      if (pdg == 6) {
        u[k].trueFlavour() = 6;
        break;
      }
      if (pdg == 5) {
        u[k].trueFlavour() = 5;
        break;
      }
    }

    std::vector<PseudoJet> cl = ljets[k].constituents();
    fastjet::ClusterSequenceArea cs_cl(cl, jet_def_small, areaDefinition);
    std::vector<fastjet::PseudoJet> subjets = sorted_by_pt(cs_cl.inclusive_jets(minSubjetPt));
    u[k].subjet().clear();
    u[k].subjet_area().clear();
 
    for (int j = 0; j < subjets.size(); ++j) {
      //u[k].subjet().push_back(subjets[j]);
      //u[k].subjet_area().push_back(subjets[j].area_4vector());
      /*fastjet::PseudoJet a = subjets[j].area_4vector();
      u[k].subjet().push_back(LargeJet::Subjet());
      u[k].subjet()[j].v.SetPxPyPzE(subjets[j].px(), subjets[j].py(), subjets[j].pz(), subjets[j].e());
      u[k].subjet()[j].lcpt = subjets[j].perp();
      u[k].subjet_area().push_back(TLorentzVector());
      u[k].subjet_area()[j].SetPxPyPzE(a.px(), a.py(), a.pz(), a.e());
*/
    }
  }
 
  // this is commented out
  // since it is not used
  // and not present in std. NTUP_COMMON
  if (m_readBB){// || Correction::globalTools->m_useCA12) 
    fastjet::AreaDefinition areaDefinitionBB(active_area);
    ClusterSequence cs_largeBB(cluster, jet_def_largeBB);
    //std::vector<PseudoJet> ljetsBB = sorted_by_pt(cs_largeBB.inclusive_jets(minLargeJetPt));
   /* for (int k = 0; k < ljetsBB.size(); ++k) {
      w.push_back(LargeJet());
      w[k].mom().SetPxPyPzE(ljetsBB[k].px(), ljetsBB[k].py(), ljetsBB[k].pz(), ljetsBB[k].e());
      w[k].split12() = -1;
      w[k].trueFlavour() = 0;
      for (int z = 0; z < m_sr->mc_n; ++z) {
        if (m_sr->mc_pt->at(z) < 50e3) continue;
        int pdg = (int) std::fabs(m_sr->mc_pdgId->at(z));
        TLorentzVector v;
        float pt = m_sr->mc_pt->at(z);
        float eta = m_sr->mc_eta->at(z);
        float phi = m_sr->mc_phi->at(z);
        float m = m_sr->mc_m->at(z);
        v.SetPtEtaPhiM(pt, eta, phi, m);
        if (v.DeltaR(w[k].mom()) > 0.3) continue;
        if (pdg == 25) {
          w[k].trueFlavour() = 25;
          break;
        }
        if (pdg == 6) {
          w[k].trueFlavour() = 6;
          break;
        }
        if (pdg == 5) {
          w[k].trueFlavour() = 5;
          break;
        }
      }

      //std::vector<PseudoJet> clBB = ljetsBB[k].constituents();
      fastjet::ClusterSequenceArea cs_clBB(clBB, jet_def_small, areaDefinitionBB);
      std::vector<fastjet::PseudoJet> subjetsBB = sorted_by_pt(cs_clBB.inclusive_jets(minSubjetPt));
      w[k].subjet().clear();
      w[k].subjet_area().clear();
 
      for (int j = 0; j < subjetsBB.size(); ++j) {
        //w[k].subjet().push_back(subjets[j]);
        //w[k].subjet_area().push_back(subjets[j].area_4vector());
        /*fastjet::PseudoJet a = subjetsBB[j].area_4vector();
        w[k].subjet().push_back(LargeJet::Subjet());
        w[k].subjet()[j].v.SetPxPyPzE(subjetsBB[j].px(), subjetsBB[j].py(), subjetsBB[j].pz(), subjetsBB[j].e());
        w[k].subjet()[j].lcpt = subjetsBB[j].perp();
        w[k].subjet_area().push_back(TLorentzVector());
        w[k].subjet_area()[j].SetPxPyPzE(a.px(), a.py(), a.pz(), a.e());
      }
    }*/
  }
}

void RawReader::doParticleLevelSelection(const bool p) {
  m_doParticleLevelSelection = p;
}

Event::DataPeriod RawReader::getPeriodFromDataRun(int run) {
  if (run <  195847) return Event::perUnknown;
  if (run <= 203227) return Event::perAtoB3;
  if (run <= 206247) return Event::perB4toB14;
  if (run <= 208178) return Event::perC6toD3;
  if (run <= 210183) return Event::perD4toX;
  return Event::perXon;
}

void RawReader::general(int &channelNumber, bool &isData, int &runNumber, int &eventNumber, float &mu, int &npv, int &npv_met, float &rho, int &lerr, int &terr, int &cfl, bool &atlFastII, Event::DataPeriod &p, float &vxZ, float &mcWeight, unsigned int &lbn, int &npv_good) {
  channelNumber = (int) m_sr->mc_channel_number;
  isData = Correction::globalTools->m_isData;
  //isData = false;
  //if (channelNumber < 0 || channelNumber > 1000000) {
  //  isData = true;
  //}
  runNumber = m_sr->RunNumber;
  eventNumber = m_sr->EventNumber;
  //hfor = m_sr->top_hfor_type;

  p = Event::perUnknown;
  if (isData) {
    p = getPeriodFromDataRun(runNumber);
  } else {
    int rn = runNumber;
    if (rn == 212399)
      rn = 195848;
    p = (Event::DataPeriod) Correction::globalTools->m_pileupTool.GetRandomPeriodNumber(rn);
  }
  
  mu = m_sr->averageIntPerXing;
  rho = m_sr->Eventshape_rhoKt4LC;
  npv = 0;
  if (m_sr->vxp_nTracks) {
    for (int i = 0; i < m_sr->vxp_nTracks->size(); ++i) {
      if (m_sr->vxp_nTracks->at(i) >= 2) {
        npv++;
      }
    }
  }

  npv_good = 0;
  if (m_sr->vxp_type && m_sr->vxp_nTracks && m_sr->vxp_n > 0) {
    if (m_sr->vxp_type->at(0) == 1 || m_sr->vxp_type->at(0) == 3)
      if (m_sr->vxp_nTracks->at(0) > 4) npv_good++;
  }


  npv_met = 0;
  if (m_sr->vxp_n > 0 && m_sr->vxp_nTracks) {
    bool goodPV = false;
    for(int i = 0; i < m_sr->vxp_n; i++) {
      if (m_sr->vxp_type->at(i) == 1 && \
          m_sr->vxp_nTracks->at(i) > 2 && \
          std::fabs(m_sr->vxp_z->at(i)) < 200) {
        goodPV = true;
      }
    }
    if (goodPV) {
      for (int i = 0; i < m_sr->vxp_n; i++) {
        if (m_sr->vxp_nTracks->at(i) > 2)
          npv_met++;
      }
    }
  }

  if (m_sr->vxp_nTracks) {
    for (int i = 0; i < m_sr->vxp_nTracks->size(); ++i) {
      if (m_sr->vxp_nTracks->at(i) >= 2) {
        npv_met++;
      }
    }
  }

  lerr = m_sr->larError;
  terr = m_sr->tileError;
  cfl = m_sr->coreFlags; 

  

  atlFastII = Correction::globalTools->m_isAtlFastII;

  vxZ = 0;
  mcWeight = 1.0;

  if (!isData){ // added by A.Knue, code crashes for data otherwise!
    if (!Correction::globalTools->m_isAtlFastII) {
      if (m_sr->mc_vx_z) {
        for (int i = 0; i < m_sr->mc_vx_z->size(); ++i) {
          if (m_sr->mc_vx_z->at(i) != 0) {
            vxZ = m_sr->mc_vx_z->at(i);
            break;
          }
        }
      }
      mcWeight = m_sr->mc_event_weight;
    }
  }

  lbn = m_sr->lbn;
}

