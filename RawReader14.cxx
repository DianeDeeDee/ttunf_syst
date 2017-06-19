#include "RawReader14.h"
#include "SkimReader.h"
#include <fastjet/PseudoJet.hh>
#include <fastjet/tools/Filter.hh>
#include "Correction.h"
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/AreaDefinition.hh>

using namespace fastjet;

RawReader14::RawReader14(SkimReader &sr)
  : Reader(sr) {
}

RawReader14::~RawReader14() {
}

void RawReader14::electron(std::vector<Electron> &v) {
  v.clear();
  /*
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
                             m_sr->el_tracketa->at(k),
                             m_sr->el_trackphi->at(k),
                             m_sr->el_cl_E->at(k));
    e.setTightPP(m_sr->el_tightPP->at(k));
    e.setMI(m_sr->el_MI10_max40_ptsum->at(k));
    e.z0() = m_sr->el_trackz0pvunbiased->at(k);
    e.author() = m_sr->el_author->at(k);
    e.nSiHits() = m_sr->el_nSiHits->at(k);
    e.oq() = m_sr->el_OQ->at(k);
    e.GSF_trk_index() = m_sr->el_GSF_trk_index->at(k);
    e.corr("",true) = e.mom();
  }*/
}

void RawReader14::muon(std::vector<Muon> &v) {
  v.clear();
  /*
  for (int k = 0; k < m_sr->mu_muid_n; ++k) {
    v.push_back(Muon());
    Muon &m = v[k];
    m.mom().SetPtEtaPhiE(m_sr->mu_muid_pt->at(k),
                         m_sr->mu_muid_eta->at(k),
                         m_sr->mu_muid_phi->at(k),
                         m_sr->mu_muid_E->at(k));
    m.setTight(m_sr->mu_muid_tight->at(k));
    m.setMI(m_sr->mu_muid_MI10_max40_ptsum->at(k));
    m.z0() = m_sr->mu_muid_trackz0pv->at(k);
    m.author() = m_sr->mu_muid_author->at(k);
    m.d0() = m_sr->mu_muid_trackd0pv->at(k);
    m.sd0() = m_sr->mu_muid_tracksigd0pv->at(k);
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
    m.corr("",true) = m.mom();
  }*/
}

void RawReader14::met(float &met_x, float &met_y) {
  met_x = m_sr->MET_RefFinal_tightpp_etx;
  met_y = m_sr->MET_RefFinal_tightpp_ety;
}

void RawReader14::trigger(bool &triggerElectron, bool &triggerMuon, bool &triggerLargeJet) {
  triggerElectron = false;//m_sr->EF_e24vhi_medium1 || m_sr->EF_e60_medium1;
  triggerMuon = false;//m_sr->EF_mu24i_tight || m_sr->EF_mu36_tight;
  triggerLargeJet = false;//m_sr -> EF_j240_a10tcem;
}




void RawReader14::jet(std::vector<Jet> &v, std::vector<LargeJet> &u, std::vector<LargeJet> &w) {
  double minLargeJetPt = 200e3;
  //double minJetPt = 25e3;
  double minSubjetPt = 20e3;
  v.clear();
  u.clear();
  w.clear();

  /*
  std::vector<PseudoJet> cluster;
  for (int k = 0; k < m_sr->cl_lc_n; ++k) {
    double px = m_sr->cl_lc_pt->at(k)*std::cos(m_sr->cl_lc_phi->at(k));
    double py = m_sr->cl_lc_pt->at(k)*std::sin(m_sr->cl_lc_phi->at(k));
    double pz = m_sr->cl_lc_pt->at(k)*std::sinh(m_sr->cl_lc_eta->at(k));
    double E  = std::sqrt(std::pow(px, 2) + std::pow(py, 2) + std::pow(pz, 2));
    if (E <= 0) continue;
    cluster.push_back(PseudoJet(px, py, pz, E+100));
  }
  
//JetDefinition jet_def_large = JetDefinition(cambridge_algorithm, 1.5);
  JetDefinition jet_def_small = JetDefinition(cambridge_algorithm, 0.2);

  for (int k = 0; k < m_sr->jet_AntiKt4LCTopo_n; ++k) {
    v.push_back(Jet());
    v[k].mom().SetPtEtaPhiE(m_sr->jet_AntiKt4LCTopo_pt->at(k),
                            m_sr->jet_AntiKt4LCTopo_eta->at(k),
                            m_sr->jet_AntiKt4LCTopo_phi->at(k),
                            m_sr->jet_AntiKt4LCTopo_E->at(k));
    v[k].trueFlavour() = m_sr->jet_AntiKt4LCTopo_flavor_truth_label->at(k);
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
    v[k].corr("",true) = v[k].mom();
  }



  fastjet::AreaDefinition areaDefinition(active_area);
  //ClusterSequenceArea cs_large(cluster, jet_def_large, areaDefinition);
  //std::vector<PseudoJet> ljets = sorted_by_pt(cs_large.inclusive_jets(minLargeJetPt));
  //for (int k = 0; k < ljets.size(); ++k) {
  //  u.push_back(LargeJet());
  //  fastjet::PseudoJet alj = ljets[k].area_4vector();
  //  u[k].mom().SetPxPyPzE(ljets[k].px(), ljets[k].py(), ljets[k].pz(), ljets[k].e());
  //  u[k].detE() = u[k].mom().E();
  //  u[k].detEta() = u[k].mom().Eta();
  //  u[k].detPhi() = u[k].mom().Phi();
  //  u[k].detM() = u[k].mom().M();
  //  u[k].Ax() = alj.px();
  //  u[k].Ay() = alj.py();
  //  u[k].Az() = alj.pz();
  //  u[k].Ae() = alj.e();
  //
  //  std::vector<PseudoJet> cl = cs_large.constituents(ljets[k]);

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
    u[k].corr("",true) = u[k].mom();


    std::vector<PseudoJet> cl;
    for (int z = 0; z < m_sr->jet_AntiKt10LCTopoTrimmedPtFrac5SmallR30_constit_n->at(k); ++z) {
      int idx = m_sr->jet_AntiKt10LCTopoTrimmedPtFrac5SmallR30_constit_index->at(k).at(z);
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
    u[k].subjet().clear();
    u[k].subjet_area().clear();
    for (int j = 0; j < subjets.size(); ++j) {
      //u[k].subjet().push_back(subjets[j]);
      //u[k].subjet_area().push_back(subjets[j].area_4vector());
      fastjet::PseudoJet a = subjets[j].area_4vector();
      u[k].subjet().push_back(TLorentzVector());
      u[k].subjet()[j].SetPxPyPzE(subjets[j].px(), subjets[j].py(), subjets[j].pz(), subjets[j].e());
      u[k].subjet_area().push_back(TLorentzVector());
      u[k].subjet_area()[j].SetPxPyPzE(a.px(), a.py(), a.pz(), a.e());
    }
  }


  if (m_sr->jet_CamKt12LCTopo_pt != 0 && m_sr->jet_CamKt12LCTopoSplitFilteredMu67SmallR0YCut9_m != 0) {
    for (int k = 0; k < m_sr->jet_CamKt12LCTopo_n; ++k) {
      w.push_back(LargeJet());
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
      w[k].corr("",true) = w[k].mom();
 
      w[k].ug_tau1() = m_sr->jet_CamKt12LCTopo_Tau1->at(k);
      w[k].ug_tau2() = m_sr->jet_CamKt12LCTopo_Tau2->at(k);
   
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
      w[k].subjet().clear();
      w[k].subjet_area().clear();
      for (int j = 0; j < subjets.size(); ++j) {
        //u[k].subjet().push_back(subjets[j]);
        //u[k].subjet_area().push_back(subjets[j].area_4vector());
        fastjet::PseudoJet a = subjets[j].area_4vector();
        w[k].subjet().push_back(TLorentzVector());
        w[k].subjet()[j].SetPxPyPzE(subjets[j].px(), subjets[j].py(), subjets[j].pz(), subjets[j].e());
        w[k].subjet_area().push_back(TLorentzVector());
        w[k].subjet_area()[j].SetPxPyPzE(a.px(), a.py(), a.pz(), a.e());
      }    
    }

    for (int j = 0; j < m_sr->jet_CamKt12LCTopoSplitFilteredMu67SmallR0YCut9_m->size(); ++j) {
      w[j].bdrs().SetPtEtaPhiM(m_sr->jet_CamKt12LCTopoSplitFilteredMu67SmallR0YCut9_pt->at(j),
                               m_sr->jet_CamKt12LCTopoSplitFilteredMu67SmallR0YCut9_eta->at(j),
                               m_sr->jet_CamKt12LCTopoSplitFilteredMu67SmallR0YCut9_phi->at(j),
                               m_sr->jet_CamKt12LCTopoSplitFilteredMu67SmallR0YCut9_m->at(j));
      w[j].bdrs_tau1() = m_sr->jet_CamKt12LCTopoSplitFilteredMu67SmallR0YCut9_Tau1->at(j);
      w[j].bdrs_tau2() = m_sr->jet_CamKt12LCTopoSplitFilteredMu67SmallR0YCut9_Tau2->at(j);
    }
  }
  */
}

void RawReader14::partElectron(std::vector<Particle> &v) {
  v.clear();
  for (int k = 0; k < m_sr->mc_n; ++k) {
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
      if (pdg == 22 && dr < 0.1) { // photon
        e.mom() += x;
      } else if (dr < drCut && x.Perp() > 500) {
        isol += x.Perp();
      }
    }
    e.setIsol(isol);
    e.corr("",true) = e.mom();
  }
}
void RawReader14::partMom(std::vector<Particle> &v) {
  v.clear();
  int pdgs[] = {25, 6, 5};
  std::vector<int> done;
  for (int k = 0; k < m_sr->mc_n; ++k) {
    if (std::find(done.begin(), done.end(), k) != done.end()) continue;
    int pdg = (int) std::fabs(m_sr->mc_pdgId->at(k));
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
          }
        } while (!last);
        if (std::find(done.begin(), done.end(), child) != done.end()) continue;
        if (m_sr->mc_pt->at(child) < 20e3) continue;
        done.push_back(child);
        v.push_back(Particle());
        v[v.size()-1].mom().SetPtEtaPhiM(m_sr->mc_pt->at(child),
                                         m_sr->mc_eta->at(child),
                                         m_sr->mc_phi->at(child),
                                         m_sr->mc_m->at(child));
        v[v.size()-1].setId( m_sr->mc_pdgId->at(child));
        break;
      }
    }
  }
}

void RawReader14::partMuon(std::vector<Particle> &v) {
  v.clear();
  for (int k = 0; k < m_sr->mc_n; ++k) {
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
    m.corr("",true) = m.mom();
  }
}

void RawReader14::partMet(float &met_x, float &met_y) {
  met_x = m_sr->MET_Truth_NonInt_etx;
  met_y = m_sr->MET_Truth_NonInt_ety;
}


 void RawReader14::partJet(std::vector<Jet> &v, std::vector<LargeJet> &u, std::vector<LargeJet> &w) {
  double minLargeJetPt = 100e3;
  //double minJetPt = 25e3;
  double minSubjetPt = 20e3;
  v.clear();
  u.clear();
  w.clear();

  std::vector<PseudoJet> cluster;
  for (int k = 0; k < m_sr->mc_n; ++k) {
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
  JetDefinition jet_def_large = JetDefinition(antikt_algorithm, 1.5);
  Filter trimmer(JetDefinition(cambridge_algorithm, 0.3), SelectorPtFractionMin(0.05));
  JetDefinition jet_def_small = JetDefinition(cambridge_algorithm, 0.2);


  fastjet::AreaDefinition areaDefinition(active_area);
  JetDefinition jet_def_normal = JetDefinition(antikt_algorithm, 0.4);
  ClusterSequenceArea cs_normal(cluster, jet_def_normal, areaDefinition);
  std::vector<PseudoJet> jets_normal = sorted_by_pt(cs_normal.inclusive_jets(20));
  for (int k = 0; k < jets_normal.size(); ++k) {
    v.push_back(Jet());
    v[k].mom().SetPxPyPzE(jets_normal[k].px(), jets_normal[k].py(), jets_normal[k].pz(), jets_normal[k].e());
  //for (int k = 0; k < m_sr->AntiKt4Truth_n; ++k) {
  //  v.push_back(Jet());
  //  v[k].mom().SetPtEtaPhiE(m_sr->AntiKt4Truth_pt->at(k),
  //                          m_sr->AntiKt4Truth_eta->at(k),
  //                          m_sr->AntiKt4Truth_phi->at(k),
  //                          m_sr->AntiKt4Truth_E->at(k));
  //  //v[k].trueFlavour() = m_sr->AntiKt4Truth_flavor_truth_label->at(k);
    v[k].trueFlavour() = 0;
    for (int z = 0; z < m_sr->mc_n; ++z) {
      if (m_sr->mc_pt->at(z) < 20e3) continue;
      int pdg = (int) std::fabs(m_sr->mc_pdgId->at(z));
      TLorentzVector t;
      float pt = m_sr->mc_pt->at(z);
      float eta = m_sr->mc_eta->at(z);
      float phi = m_sr->mc_phi->at(z);
      float m = m_sr->mc_m->at(z);
      t.SetPtEtaPhiM(pt, eta, phi, m);
      if (t.DeltaR(v[k].mom()) > 0.3) continue;
      if (pdg == 5) {
        v[k].trueFlavour() = 5;
        break;
      }
    }
    v[k].mv1() = 0;
    v[k].jvf() = 1;
    if (v[k].trueFlavour() == 5)
      v[k].mv1() = 99;
    v[k].corr("",true) = v[k].mom();
  }

  ClusterSequenceArea cs_large(cluster, jet_def_large, areaDefinition);
  std::vector<PseudoJet> ljets = sorted_by_pt(cs_large.inclusive_jets(minLargeJetPt));
  for (int k = 0; k < ljets.size(); ++k) {
    u.push_back(LargeJet());
    //PseudoJet trimjet = trimmer(ljets[k]);
    u[k].mom().SetPxPyPzE(ljets[k].px(), ljets[k].py(), ljets[k].pz(), ljets[k].e());
    u[k].split12() = -1;
    u[k].corr("",true) = u[k].mom();
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
    std::vector<PseudoJet> cl = cs_large.constituents(ljets[k]);
    fastjet::ClusterSequenceArea cs_cl(cl, jet_def_small, areaDefinition);
    std::vector<fastjet::PseudoJet> subjets = sorted_by_pt(cs_cl.inclusive_jets(minSubjetPt));
    u[k].subjet().clear();
    u[k].subjet_area().clear();
    for (int j = 0; j < subjets.size(); ++j) {
      //u[k].subjet().push_back(subjets[j]);
      //u[k].subjet_area().push_back(subjets[j].area_4vector());
      fastjet::PseudoJet a = subjets[j].area_4vector();
      u[k].subjet().push_back(TLorentzVector());
      u[k].subjet()[j].SetPxPyPzE(subjets[j].px(), subjets[j].py(), subjets[j].pz(), subjets[j].e());
      u[k].subjet_area().push_back(TLorentzVector());
      u[k].subjet_area()[j].SetPxPyPzE(a.px(), a.py(), a.pz(), a.e());

    }
  }

  /*
   for (int k = 0; k < m_sr->jet_CamKt12Truth_n; ++k) {
    w.push_back(LargeJet());
    w[k].mom().SetPtEtaPhiE(m_sr->jet_CamKt12Truth_pt->at(k), \
                            m_sr->jet_CamKt12Truth_eta->at(k), \
                            m_sr->jet_CamKt12Truth_phi->at(k), \
                            m_sr->jet_CamKt12Truth_E->at(k));
 
 
    w[k].corr("",true) = w[k].mom();
    std::vector<PseudoJet> cl;
    for (int z = 0; z < m_sr->jet_CamKt12Truth_constit_n->at(k); ++z) {
      int idx = m_sr->jet_CamKt12Truth_constit_index->at(k).at(z);
      cl.push_back(fastjet::PseudoJet());
      TLorentzVector v;
      v.SetPtEtaPhiM(m_sr->mc_pt->at(idx), \
                     m_sr->mc_eta->at(idx), \
                     m_sr->mc_phi->at(idx), \
                     m_sr->mc_m->at(idx));
      cl[z].reset_momentum(v.Px(), v.Py(), v.Pz(), v.E());
    }
    fastjet::AreaDefinition areaDefinition(active_area);
    fastjet::ClusterSequenceArea cs_cl(cl, jet_def_small, areaDefinition);
    std::vector<fastjet::PseudoJet> subjets = sorted_by_pt(cs_cl.inclusive_jets(minSubjetPt));
    w[k].subjet().clear();
    w[k].subjet_area().clear();
    for (int j = 0; j < subjets.size(); ++j) {
      fastjet::PseudoJet a = subjets[j].area_4vector();
      w[k].subjet().push_back(TLorentzVector());
      w[k].subjet()[j].SetPxPyPzE(subjets[j].px(), subjets[j].py(), subjets[j].pz(), subjets[j].e());
      w[k].subjet_area().push_back(TLorentzVector());
      w[k].subjet_area()[j].SetPxPyPzE(a.px(), a.py(), a.pz(), a.e());


    }
  }*/
}

void RawReader14::general(int &channelNumber, bool &isData, int &runNumber, int &eventNumber, float &mu, int &npv, int &npv_met, float &rho, int &lerr, int &terr, int &cfl, bool &atlFastII, Event::DataPeriod &p, float &vxZ, float &mcWeight, unsigned int &lbn) {
  channelNumber = (int) m_sr->mc_channel_number;
  isData = false;
  if (channelNumber < 0 || channelNumber > 1000000) {
    isData = true;
  }
  runNumber = m_sr->RunNumber;
  eventNumber = m_sr->EventNumber;

  p = Event::perUnknown;
  
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
  for (int i = 0; i < m_sr->mc_vx_z->size(); ++i) {
    if (m_sr->mc_vx_z->at(i) != 0) {
      vxZ = m_sr->mc_vx_z->at(i);
      break;
    }
  }
  mcWeight = m_sr->mc_event_weight;
  lbn = m_sr->lbn;
}

