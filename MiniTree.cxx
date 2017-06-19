#include "MiniTree.h"
#include "Event.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include <vector>
#include <string>
#include "TChain.h"

MiniTree::MiniTree(bool toWrite, const std::string &file)
  : m_toWrite(toWrite), m_fileToWrite(0), m_chain(0) {
  if (toWrite) {
    m_fileToWrite = new TFile(file.c_str(), "RECREATE");
    m_sfs = new std::vector<std::string>;
    m_systs = new std::vector<std::string>;
    m_chain = new TTree("analysis", "");
    m_num = new TTree("num", "");
  } else {
    m_fileToWrite = TFile::Open(file.c_str());
    std::vector<std::string> *sfs = (std::vector<std::string> *) m_fileToWrite->Get("weightNames");
    std::vector<std::string> *systs = (std::vector<std::string> *) m_fileToWrite->Get("systematicsNames");
    m_sfs = new std::vector<std::string>;
    m_systs = new std::vector<std::string>;
    if (sfs)
      *m_sfs = *sfs;
    if (systs)
      *m_systs = *systs;
    m_chain = new TChain("analysis");
    m_num = new TChain("num");
    ((TChain *) m_chain)->Add(file.c_str());
    ((TChain *) m_num)->Add(file.c_str());
  }
  m_sumWeights = 0;
  m_sumWeights_var = 0;

  m_passCuts = 0;

  prepareBranches();
}

MiniTree::~MiniTree() {
  if (m_toWrite) {
    m_fileToWrite->cd();
    m_chain->Write();
    m_num->Fill();
    m_num->Write();
    gDirectory->WriteObject(m_sfs, "weightNames");
    gDirectory->WriteObject(m_systs, "systematicsNames");
  } else {
    if (m_sfs) delete m_sfs;
    if (m_systs) delete m_systs;
  }
  if (m_chain) delete m_chain;
  if (m_fileToWrite) delete m_fileToWrite;
}

std::vector<std::string> *MiniTree::weightNames() {
  return m_sfs;
}

std::vector<std::string> *MiniTree::systematicsNames() {
  return m_systs;
}

void MiniTree::read(int event, Event &e, int &_systIdx) {
  e.clear();
  m_chain->GetEntry(event);
  //e.hfor() = hfor;
  e.channelNumber() = mc_channel_number;
  _systIdx = systIdx;
  e.passReco() = passReco;
  e.passPart() = passPart;
  e.triggerElectron() = triggerElectron;
  /*e.triggerElectron24() = triggerElectron24;
  e.triggerElectron60() = triggerElectron60;
  */e.triggerMuon() = triggerMuon;
  e.triggerLargeJet() = triggerLargeJet;

  e.isTight() = isTight;

  e.npv() = npv;
  e.mu() = avmu;

  for (int k = 0; k < el_n; ++k) {
    e.electron().push_back(Electron());
    e.electron()[k].mom().SetPtEtaPhiE(el_pt->at(k), el_eta->at(k), el_phi->at(k), el_E->at(k));
    e.electron()[k].setMI(0);
    e.electron()[k].setTightPP(true);
    e.electron()[k].caloMom() = e.electron()[k].mom();
    e.electron()[k].trkMom() = e.electron()[k].mom();
    e.electron()[k].z0() = 0;
    e.electron()[k].author() = 1;
  }
  for (int k = 0; k < mu_n; ++k) {
    e.muon().push_back(Muon());
    e.muon()[k].mom().SetPtEtaPhiE(mu_pt->at(k), mu_eta->at(k), mu_phi->at(k), mu_E->at(k));
    e.muon()[k].setMI(0);
    e.muon()[k].setTight(true);
    e.muon()[k].z0() = 0;
    e.muon()[k].d0() = 0;
    e.muon()[k].sd0() = 0;
    e.muon()[k].author() = 0;
    e.muon()[k].passTrkCuts() = true;
  }

  for (int k = 0; k < jet_n; ++k) {
    e.jet().push_back(Jet());
    e.jet()[k].mom().SetPtEtaPhiE(jet_pt->at(k), jet_eta->at(k), jet_phi->at(k), jet_E->at(k));
    e.jet()[k].trueFlavour() = jet_trueflav==0?-99:jet_trueflav->at(k);
    e.jet()[k].mv1() = jet_mv1==0?-99:jet_mv1->at(k);
  }
  e.triggerJet4mom().SetPtEtaPhiE(trigjet_pt, trigjet_eta, trigjet_phi, trigjet_E);

  for (int k = 0; k < ljet_n; ++k) {
    e.largeJet().push_back(LargeJet());
    e.largeJet()[k].mom().SetPtEtaPhiE(ljet_pt->at(k), ljet_eta->at(k), ljet_phi->at(k), ljet_E->at(k));
    e.largeJet()[k].split12() = ljet_split12->at(k);
    e.largeJet()[k].trueFlavour() = ljet_trueflav==0?-99:ljet_trueflav->at(k);
    /*for (int l = 0; l < ljet_sub_n->at(k); ++l) {
      e.largeJet()[k].subjet().push_back(LargeJet::Subjet());
      e.largeJet()[k].subjet()[l].v.SetPtEtaPhiE(ljet_sub_pt->at(k).at(l), ljet_sub_eta->at(k).at(l), ljet_sub_phi->at(k).at(l), ljet_sub_E->at(k).at(l));
      e.largeJet()[k].subjet()[l].lcpt = ljet_sub_lcpt->at(k).at(l);
      e.largeJet()[k].trimmed_tau1() = ljet_trimmed_tau1->at(k);
      e.largeJet()[k].trimmed_tau2() = ljet_trimmed_tau2->at(k);
      e.largeJet()[k].trimmed_am_tau1() = ljet_trimmed_am_tau1->at(k);
      e.largeJet()[k].trimmed_am_tau2() = ljet_trimmed_am_tau2->at(k);
      if (ljet_htt) e.largeJet()[k].htt() = ljet_htt->at(k);
    }*/
    if (ljet_subunc_n) {
    for (int l = 0; l < ljet_subunc_n->at(k); ++l) {
      e.largeJet()[k].subjetJER().push_back(LargeJet::Subjet());
      e.largeJet()[k].subjetJER()[l].v.SetPtEtaPhiE(ljet_subunc_pt->at(k).at(l), ljet_subunc_eta->at(k).at(l), ljet_subunc_phi->at(k).at(l), ljet_subunc_E->at(k).at(l));
      e.largeJet()[k].subjetJER()[l].lcpt = ljet_subunc_lcpt->at(k).at(l);
    }
    }
  }

  for (int k = 0; k < ljetBB_n; ++k) {
    /*e.largeJetBB().push_back(LargeJet());
    e.largeJetBB()[k].mom().SetPtEtaPhiE(ljetBB_pt->at(k), ljetBB_eta->at(k), ljetBB_phi->at(k), ljetBB_E->at(k));
    e.largeJetBB()[k].split12() = ljetBB_split12->at(k);
    e.largeJetBB()[k].ug_tau1() = ljetBB_ug_tau1->at(k);
    e.largeJetBB()[k].ug_tau2() = ljetBB_ug_tau2->at(k);
    e.largeJetBB()[k].ug_am_tau1() = ljetBB_ug_am_tau1->at(k);
    e.largeJetBB()[k].ug_am_tau2() = ljetBB_ug_am_tau2->at(k);

    for (int l = 0; l < ljetBB_sub_n->at(k); ++l) {
      e.largeJetBB()[k].subjet().push_back(LargeJet::Subjet());
      e.largeJetBB()[k].subjet()[l].v.SetPtEtaPhiE(ljetBB_sub_pt->at(k).at(l), ljetBB_sub_eta->at(k).at(l), ljetBB_sub_phi->at(k).at(l), ljetBB_sub_E->at(k).at(l));
      e.largeJetBB()[k].subjet()[l].lcpt = ljetBB_sub_lcpt->at(k).at(l);
    }
    if (ljetBB_subunc_n) {
    for (int l = 0; l < ljetBB_subunc_n->at(k); ++l) {
      e.largeJetBB()[k].subjetJER().push_back(LargeJet::Subjet());
      e.largeJetBB()[k].subjetJER()[l].v.SetPtEtaPhiE(ljetBB_subunc_pt->at(k).at(l), ljetBB_subunc_eta->at(k).at(l), ljetBB_subunc_phi->at(k).at(l), ljetBB_subunc_E->at(k).at(l));
      e.largeJetBB()[k].subjetJER()[l].lcpt = ljetBB_subunc_lcpt->at(k).at(l);
    }
    }
    e.largeJetBB()[k].bdrs().SetPtEtaPhiM(ljetBB_bdrs_pt->at(k),
                                          //ljetBB_bdrs_eta->at(k), // DANILO TEMPORARY REMOVAL
                                          //ljetBB_bdrs_phi->at(k), // DANILO TEMPORARY REMOVAL
                                          ljetBB_eta->at(k), // DANILO TEMPORARY REMOVAL
                                          ljetBB_phi->at(k), // DANILO TEMPORARY REMOVAL
                                          ljetBB_bdrs_m->at(k));
    e.largeJetBB()[k].ug_bdrs_tau1() = ljetBB_ug_bdrs_tau1->at(k);
    e.largeJetBB()[k].ug_bdrs_tau2() = ljetBB_ug_bdrs_tau2->at(k);
    if (ljetBB_htt) e.largeJetBB()[k].htt() = ljetBB_htt->at(k);
  */}

  e.met(met_etx, met_ety);
  e.triggerElectron() = triggerElectron;
  /*e.triggerElectron24() = triggerElectron24;
  e.triggerElectron60() = triggerElectron60;
 */ e.triggerMuon() = triggerMuon;
  e.triggerLargeJet() = triggerLargeJet;

  for (int k = 0; k < partel_n; ++k) {
    e.partElectron().push_back(Particle());
    e.partElectron()[k].mom().SetPtEtaPhiE(partel_pt->at(k), partel_eta->at(k), partel_phi->at(k), partel_E->at(k));
    e.partElectron()[k].setIsol(0);
    e.partElectron()[k].setId(11);
  }
  for (int k = 0; k < partmu_n; ++k) {
    e.partMuon().push_back(Particle());
    e.partMuon()[k].mom().SetPtEtaPhiE(partmu_pt->at(k), partmu_eta->at(k), partmu_phi->at(k), partmu_E->at(k));
    e.partMuon()[k].setIsol(0);
    e.partMuon()[k].setId(13);
  }

  for (int k = 0; k < partjet_n; ++k) {
    e.partJet().push_back(Jet());
    e.partJet()[k].mom().SetPtEtaPhiE(partjet_pt->at(k), partjet_eta->at(k), partjet_phi->at(k), partjet_E->at(k));
    e.partJet()[k].trueFlavour() = partjet_trueflav->at(k);
    e.partJet()[k].mv1() = 99;
  }

  for (int k = 0; k < partmom_n; ++k) {
    e.partMom().push_back(Particle());
    e.partMom()[k].mom().SetPtEtaPhiE(partmom_pt->at(k), partmom_eta->at(k), partmom_phi->at(k), partmom_E->at(k));
    e.partMom()[k].setId(partmom_trueflav->at(k));
  }

  for (int k = 0; k < partljet_n; ++k) {
    e.partLargeJet().push_back(LargeJet());
    e.partLargeJet()[k].mom().SetPtEtaPhiE(partljet_pt->at(k), partljet_eta->at(k), partljet_phi->at(k), partljet_E->at(k));
    e.partLargeJet()[k].split12() = partljet_split12->at(k);
    e.partLargeJet()[k].trueFlavour() = partljet_trueflav->at(k);
    for (int l = 0; l < partljet_sub_n->at(k); ++l) {
      e.partLargeJet()[k].subjet().push_back(LargeJet::Subjet());
      e.partLargeJet()[k].subjet()[l].v.SetPtEtaPhiE(partljet_sub_pt->at(k).at(l), partljet_sub_eta->at(k).at(l), partljet_sub_phi->at(k).at(l), partljet_sub_E->at(k).at(l));
    }
  }

  for (int k = 0; k < partljetBB_n; ++k) {
    /*e.partLargeJetBB().push_back(LargeJet());
    e.partLargeJetBB()[k].mom().SetPtEtaPhiE(partljetBB_pt->at(k), partljetBB_eta->at(k), partljetBB_phi->at(k), partljetBB_E->at(k));
    e.partLargeJetBB()[k].split12() = partljetBB_split12->at(k);
    for (int l = 0; l < partljetBB_sub_n->at(k); ++l) {
      e.partLargeJetBB()[k].subjet().push_back(LargeJet::Subjet());
      e.partLargeJetBB()[k].subjet()[l].v.SetPtEtaPhiE(partljetBB_sub_pt->at(k).at(l), partljetBB_sub_eta->at(k).at(l), partljetBB_sub_phi->at(k).at(l), partljetBB_sub_E->at(k).at(l));
    }*/
  }


  e.partMet(partmet_etx, partmet_ety);

  for (int k = 0; k < weight->size(); ++k) {
    e.weight(m_sfs->at(k), true) = weight->at(k);
  }
}

double &MiniTree::sumWeights() {
  return m_sumWeights;
}

std::vector<double> &MiniTree::sumWeights_var() {
  return *m_sumWeights_var;
}

std::vector<float> &MiniTree::passCuts() {
  return *m_passCuts;
}

void MiniTree::write(const Event &e, const int _systIdx) {
  // TODO

  //hfor = e.hfor();
  mc_channel_number = e.channelNumber();
  triggerElectron = e.triggerElectron();
  /*triggerElectron24 = e.triggerElectron24();
  triggerElectron60 = e.triggerElectron60();
  */triggerMuon = e.triggerMuon();
  triggerLargeJet = e.triggerLargeJet();

  isTight = e.isTight();

  passReco = e.passReco();
  passPart = e.passPart();

  npv = e.npv();
  avmu = e.mu();


  systIdx = _systIdx;

  el_n = 0;
  el_pt->clear();
  el_eta->clear();
  el_phi->clear();
  el_E->clear();
  for (int k = 0; k < e.electron().size(); ++k) {
    el_pt->push_back(e.electron()[k].mom().Perp());
    el_eta->push_back(e.electron()[k].mom().Eta());
    el_phi->push_back(e.electron()[k].mom().Phi());
    el_E->push_back(e.electron()[k].mom().E());
    el_n++;
  }

  mu_n = 0;
  mu_pt->clear();
  mu_eta->clear();
  mu_phi->clear();
  mu_E->clear();
  for (int k = 0; k < e.muon().size(); ++k) {
    mu_pt->push_back(e.muon()[k].mom().Perp());
    mu_eta->push_back(e.muon()[k].mom().Eta());
    mu_phi->push_back(e.muon()[k].mom().Phi());
    mu_E->push_back(e.muon()[k].mom().E());
    mu_n++;
  }

  jet_n = 0;
  jet_pt->clear();
  jet_eta->clear();
  jet_phi->clear();
  jet_E->clear();
  jet_mv1->clear();
  jet_trueflav->clear();
  for (int k = 0; k < e.jet().size(); ++k) {
    jet_pt->push_back(e.jet()[k].mom().Perp());
    jet_eta->push_back(e.jet()[k].mom().Eta());
    jet_phi->push_back(e.jet()[k].mom().Phi());
    jet_E->push_back(e.jet()[k].mom().E());
    jet_mv1->push_back(e.jet()[k].mv1());
    jet_trueflav->push_back(e.jet()[k].trueFlavour());
    jet_n++;
  }
  trigjet_pt = e.triggerJet4mom().Perp();
  trigjet_eta = e.triggerJet4mom().Eta();
  trigjet_phi = e.triggerJet4mom().Phi();
  trigjet_E = e.triggerJet4mom().E();

  ljet_n = 0;
  ljet_pt->clear();
  ljet_eta->clear();
  ljet_phi->clear();
  ljet_E->clear();
  ljet_split12->clear();
  ljet_trueflav->clear();
  ljet_sub_n->clear();
  ljet_sub_pt->clear();
  ljet_sub_eta->clear();
  ljet_sub_phi->clear();
  ljet_sub_E->clear();
  ljet_sub_lcpt->clear();
  ljet_subunc_n->clear();
  ljet_subunc_pt->clear();
  ljet_subunc_eta->clear();
  ljet_subunc_phi->clear();
  ljet_subunc_E->clear();
  ljet_subunc_lcpt->clear();
  ljet_trimmed_tau1->clear();
  ljet_trimmed_tau2->clear();
  ljet_trimmed_am_tau1->clear();
  ljet_trimmed_am_tau2->clear();
  ljet_htt->clear();

  ljetBB_n = 0;
  ljetBB_pt->clear();
  ljetBB_eta->clear();
  ljetBB_phi->clear();
  ljetBB_E->clear();
  ljetBB_split12->clear();
  ljetBB_sub_n->clear();
  ljetBB_sub_pt->clear();
  ljetBB_sub_eta->clear();
  ljetBB_sub_phi->clear();
  ljetBB_sub_E->clear();
  ljetBB_sub_lcpt->clear();
  ljetBB_subunc_n->clear();
  ljetBB_subunc_pt->clear();
  ljetBB_subunc_eta->clear();
  ljetBB_subunc_phi->clear();
  ljetBB_subunc_E->clear();
  ljetBB_subunc_lcpt->clear();

  ljetBB_ug_tau1->clear();
  ljetBB_ug_tau2->clear();
  ljetBB_bdrs_pt->clear(); 
  ljetBB_bdrs_m->clear();
  ljetBB_bdrs_phi->clear(); 
  ljetBB_bdrs_eta->clear();
  ljetBB_ug_bdrs_tau1->clear();
  ljetBB_ug_bdrs_tau2->clear();
  ljetBB_ug_tau1->clear();
  ljetBB_ug_tau2->clear();
  ljetBB_ug_am_tau1->clear();
  ljetBB_ug_am_tau2->clear();
  ljetBB_htt->clear();



  for (int k = 0; k < e.largeJet().size(); ++k) {
    ljet_pt->push_back(e.largeJet()[k].mom().Perp());
    ljet_eta->push_back(e.largeJet()[k].mom().Eta());
    ljet_phi->push_back(e.largeJet()[k].mom().Phi());
    ljet_E->push_back(e.largeJet()[k].mom().E());
    ljet_split12->push_back(e.largeJet()[k].split12());
    ljet_trueflav->push_back(e.largeJet()[k].trueFlavour());
    ljet_sub_n->push_back(0);
    ljet_sub_pt->push_back(std::vector<float>());
    ljet_sub_eta->push_back(std::vector<float>());
    ljet_sub_phi->push_back(std::vector<float>());
    ljet_sub_E->push_back(std::vector<float>());
    ljet_sub_lcpt->push_back(std::vector<float>());
    ljet_subunc_n->push_back(0);
    ljet_subunc_pt->push_back(std::vector<float>());
    ljet_subunc_eta->push_back(std::vector<float>());
    ljet_subunc_phi->push_back(std::vector<float>());
    ljet_subunc_E->push_back(std::vector<float>());
    ljet_subunc_lcpt->push_back(std::vector<float>());
    ljet_trimmed_tau1->push_back(e.largeJet()[k].trimmed_tau1());
    ljet_trimmed_tau2->push_back(e.largeJet()[k].trimmed_tau2());
    ljet_trimmed_am_tau1->push_back(e.largeJet()[k].trimmed_am_tau1());
    ljet_trimmed_am_tau2->push_back(e.largeJet()[k].trimmed_am_tau2());
    for (int l = 0; l < e.largeJet()[k].subjet().size(); ++l) {
      const TLorentzVector &v = e.largeJet()[k].subjet()[l].v;
      ljet_sub_pt->at(k).push_back(v.Perp());
      ljet_sub_eta->at(k).push_back(v.Eta());
      ljet_sub_phi->at(k).push_back(v.Phi());
      ljet_sub_E->at(k).push_back(v.E());
      ljet_sub_lcpt->at(k).push_back(e.largeJet()[k].subjet()[l].lcpt);
      (*ljet_sub_n)[k]++;
    }
    for (int l = 0; l < e.largeJet()[k].subjetJER().size(); ++l) {
      const TLorentzVector &v = e.largeJet()[k].subjetJER()[l].v;
      ljet_subunc_pt->at(k).push_back(v.Perp());
      ljet_subunc_eta->at(k).push_back(v.Eta());
      ljet_subunc_phi->at(k).push_back(v.Phi());
      ljet_subunc_E->at(k).push_back(v.E());
      ljet_subunc_lcpt->at(k).push_back(e.largeJet()[k].subjetJER()[l].lcpt);
      (*ljet_subunc_n)[k]++;
    }
    ljet_htt->push_back(e.largeJet()[k].htt());
    ljet_n++;
  }

  for (int k = 0; k < e.largeJetBB().size(); ++k) {
    ljetBB_pt->push_back(e.largeJetBB()[k].mom().Perp());
    ljetBB_eta->push_back(e.largeJetBB()[k].mom().Eta());
    ljetBB_phi->push_back(e.largeJetBB()[k].mom().Phi());
    ljetBB_E->push_back(e.largeJetBB()[k].mom().E());
    ljetBB_split12->push_back(e.largeJetBB()[k].split12());
    ljetBB_ug_tau1->push_back(e.largeJetBB()[k].ug_tau1());
    ljetBB_ug_tau2->push_back(e.largeJetBB()[k].ug_tau2());
    ljetBB_ug_am_tau1->push_back(e.largeJetBB()[k].ug_am_tau1());
    ljetBB_ug_am_tau2->push_back(e.largeJetBB()[k].ug_am_tau2());
    ljetBB_sub_n->push_back(0);
    ljetBB_sub_pt->push_back(std::vector<float>());
    ljetBB_sub_eta->push_back(std::vector<float>());
    ljetBB_sub_phi->push_back(std::vector<float>());
    ljetBB_sub_E->push_back(std::vector<float>());
    ljetBB_sub_lcpt->push_back(std::vector<float>());
    for (int l = 0; l < e.largeJetBB()[k].subjet().size(); ++l) {
      const TLorentzVector &v = e.largeJetBB()[k].subjet()[l].v;
      ljetBB_sub_pt->at(k).push_back(v.Perp());
      ljetBB_sub_eta->at(k).push_back(v.Eta());
      ljetBB_sub_phi->at(k).push_back(v.Phi());
      ljetBB_sub_E->at(k).push_back(v.E());
      ljetBB_sub_lcpt->at(k).push_back(e.largeJetBB()[k].subjet()[l].lcpt);
      (*ljetBB_sub_n)[k]++;
    }
    ljetBB_subunc_n->push_back(0);
    ljetBB_subunc_pt->push_back(std::vector<float>());
    ljetBB_subunc_eta->push_back(std::vector<float>());
    ljetBB_subunc_phi->push_back(std::vector<float>());
    ljetBB_subunc_E->push_back(std::vector<float>());
    ljetBB_subunc_lcpt->push_back(std::vector<float>());
    for (int l = 0; l < e.largeJetBB()[k].subjetJER().size(); ++l) {
      const TLorentzVector &v = e.largeJetBB()[k].subjetJER()[l].v;
      ljetBB_subunc_pt->at(k).push_back(v.Perp());
      ljetBB_subunc_eta->at(k).push_back(v.Eta());
      ljetBB_subunc_phi->at(k).push_back(v.Phi());
      ljetBB_subunc_E->at(k).push_back(v.E());
      ljetBB_subunc_lcpt->at(k).push_back(e.largeJetBB()[k].subjetJER()[l].lcpt);
      (*ljetBB_subunc_n)[k]++;
    }
    ljetBB_bdrs_m->push_back(e.largeJetBB()[k].bdrs_m());
    ljetBB_bdrs_pt->push_back(e.largeJetBB()[k].bdrs_pt());
    ljetBB_bdrs_eta->push_back(e.largeJetBB()[k].bdrs_eta());
    ljetBB_bdrs_phi->push_back(e.largeJetBB()[k].bdrs_phi());
    ljetBB_ug_bdrs_tau1->push_back(e.largeJetBB()[k].ug_bdrs_tau1());
    ljetBB_ug_bdrs_tau2->push_back(e.largeJetBB()[k].ug_bdrs_tau2());
    ljetBB_htt->push_back(e.largeJetBB()[k].htt());
    ljetBB_n++;
  }

  met_etx = e.met().Px();
  met_ety = e.met().Py();

  partel_n = 0;
  partel_pt->clear();
  partel_eta->clear();
  partel_phi->clear();
  partel_E->clear();
  for (int k = 0; k < e.partElectron().size(); ++k) {
    partel_pt->push_back(e.partElectron()[k].mom().Perp());
    partel_eta->push_back(e.partElectron()[k].mom().Eta());
    partel_phi->push_back(e.partElectron()[k].mom().Phi());
    partel_E->push_back(e.partElectron()[k].mom().E());
    partel_n++;
  }

  partmu_n = 0;
  partmu_pt->clear();
  partmu_eta->clear();
  partmu_phi->clear();
  partmu_E->clear();
  for (int k = 0; k < e.partMuon().size(); ++k) {
    partmu_pt->push_back(e.partMuon()[k].mom().Perp());
    partmu_eta->push_back(e.partMuon()[k].mom().Eta());
    partmu_phi->push_back(e.partMuon()[k].mom().Phi());
    partmu_E->push_back(e.partMuon()[k].mom().E());
    partmu_n++;
  }

  partjet_n = 0;
  partjet_pt->clear();
  partjet_eta->clear();
  partjet_phi->clear();
  partjet_E->clear();
  partjet_trueflav->clear();
  for (int k = 0; k < e.partJet().size(); ++k) {
    partjet_pt->push_back(e.partJet()[k].mom().Perp());
    partjet_eta->push_back(e.partJet()[k].mom().Eta());
    partjet_phi->push_back(e.partJet()[k].mom().Phi());
    partjet_E->push_back(e.partJet()[k].mom().E());
    partjet_trueflav->push_back(e.partJet()[k].trueFlavour());
    partjet_n++;
  }

  partmom_n = 0;
  partmom_pt->clear();
  partmom_eta->clear();
  partmom_phi->clear();
  partmom_E->clear();
  partmom_trueflav->clear();

  partljet_n = 0;
  partljet_pt->clear();
  partljet_eta->clear();
  partljet_phi->clear();
  partljet_E->clear();
  partljet_split12->clear();
  partljet_trueflav->clear();
  partljet_sub_n->clear();
  partljet_sub_pt->clear();
  partljet_sub_eta->clear();
  partljet_sub_phi->clear();
  partljet_sub_E->clear();

  partljetBB_n = 0;
  partljetBB_pt->clear();
  partljetBB_eta->clear();
  partljetBB_phi->clear();
  partljetBB_E->clear();
  partljetBB_split12->clear();
  partljetBB_sub_n->clear();
  partljetBB_sub_pt->clear();
  partljetBB_sub_eta->clear();
  partljetBB_sub_phi->clear();
  partljetBB_sub_E->clear();


  partmom_n = 0;
  for (int k = 0; k < e.partMom().size(); ++k) {
    partmom_pt->push_back(e.partMom()[k].mom().Perp());
    partmom_eta->push_back(e.partMom()[k].mom().Eta());
    partmom_phi->push_back(e.partMom()[k].mom().Phi());
    partmom_E->push_back(e.partMom()[k].mom().E());
    partmom_trueflav->push_back(e.partMom()[k].id());
    partmom_n++;
  }

  for (int k = 0; k < e.partLargeJet().size(); ++k) {
    partljet_pt->push_back(e.partLargeJet()[k].mom().Perp());
    partljet_eta->push_back(e.partLargeJet()[k].mom().Eta());
    partljet_phi->push_back(e.partLargeJet()[k].mom().Phi());
    partljet_E->push_back(e.partLargeJet()[k].mom().E());
    partljet_split12->push_back(e.partLargeJet()[k].split12());
    partljet_trueflav->push_back(e.partLargeJet()[k].trueFlavour());
    partljet_sub_n->push_back(0);
    partljet_sub_pt->push_back(std::vector<float>());
    partljet_sub_eta->push_back(std::vector<float>());
    partljet_sub_phi->push_back(std::vector<float>());
    partljet_sub_E->push_back(std::vector<float>());
    for (int l = 0; l < e.partLargeJet()[k].subjet().size(); ++l) {
      const TLorentzVector &v = e.partLargeJet()[k].subjet()[l].v;
      partljet_sub_pt->at(k).push_back(v.Perp());
      partljet_sub_eta->at(k).push_back(v.Eta());
      partljet_sub_phi->at(k).push_back(v.Phi());
      partljet_sub_E->at(k).push_back(v.E());
      (*partljet_sub_n)[k]++;
    }
    partljet_n++;
  }





  for (int k = 0; k < e.partLargeJetBB().size(); ++k) {
    /*partljetBB_pt->push_back(e.partLargeJetBB()[k].mom().Perp());
    partljetBB_eta->push_back(e.partLargeJetBB()[k].mom().Eta());
    partljetBB_phi->push_back(e.partLargeJetBB()[k].mom().Phi());
    partljetBB_E->push_back(e.partLargeJetBB()[k].mom().E());
    partljetBB_split12->push_back(e.partLargeJetBB()[k].split12());
    partljetBB_sub_n->push_back(0);
    partljetBB_sub_pt->push_back(std::vector<float>());
    partljetBB_sub_eta->push_back(std::vector<float>());
    partljetBB_sub_phi->push_back(std::vector<float>());
    partljetBB_sub_E->push_back(std::vector<float>());
    for (int l = 0; l < e.partLargeJetBB()[k].subjet().size(); ++l) {
      const TLorentzVector &v = e.partLargeJetBB()[k].subjet()[l].v;
      partljetBB_sub_pt->at(k).push_back(v.Perp());
      partljetBB_sub_eta->at(k).push_back(v.Eta());
      partljetBB_sub_phi->at(k).push_back(v.Phi());
      partljetBB_sub_E->at(k).push_back(v.E());
      (*partljetBB_sub_n)[k]++;
    }
    partljetBB_n++;*/
  }






  partmet_etx = e.partMet().Px();
  partmet_ety = e.partMet().Py();

  weight->clear();
  std::vector<std::string> wnames = e.weightNames();
  for (int k = 0; k < wnames.size(); ++k) {
    weight->push_back(e.weight(wnames[k]));
  }

  m_chain->Fill();
}

void MiniTree::addFileToRead(const std::string &fname) {
  ((TChain *) m_chain)->Add(fname.c_str());
  ((TChain *) m_num)->Add(fname.c_str());
}

double MiniTree::getSumWeights() {
  double sum = 0;
  for (int k = 0; k < m_num->GetEntries(); ++k) {
    m_num->GetEntry(k);
    sum += m_sumWeights;
  }
  return sum;
}

std::vector<float> MiniTree::getPassCuts() {
  std::vector<float> sum;
  for (int l = 0; l < m_passCuts->size(); l++) {
    sum.push_back(0);
  }
  for (int k = 0;k < m_num->GetEntries();++k) {
    m_num->GetEntry(k);
    for (int l = 0; l < m_passCuts->size(); l++) {
      sum[l] += m_passCuts->at(l);
    }
  }
  return sum;
}

int MiniTree::GetEntries() {
  return m_chain->GetEntries();
}

void MiniTree::prepareBranches() {
  systIdx = 0;
  //hfor = 0;

  el_pt = 0;
  el_eta = 0;
  el_phi = 0;
  el_E = 0;

  mu_pt = 0;
  mu_eta = 0;
  mu_phi = 0;
  mu_E = 0;

  jet_pt = 0;
  jet_eta = 0;
  jet_phi = 0;
  jet_E = 0;
  jet_mv1 = 0;
  jet_trueflav = 0;

  ljet_pt = 0;
  ljet_eta = 0;
  ljet_phi = 0;
  ljet_E = 0;
  ljet_split12 = 0;
  ljet_trueflav = 0;
  ljet_sub_n = 0;
  ljet_sub_pt = 0;
  ljet_sub_eta = 0;
  ljet_sub_phi = 0;
  ljet_sub_E = 0;
  ljet_sub_lcpt = 0;
  ljet_subunc_n = 0;
  ljet_subunc_pt = 0;
  ljet_subunc_eta = 0;
  ljet_subunc_phi = 0;
  ljet_subunc_E = 0;
  ljet_subunc_lcpt = 0;
  ljet_trimmed_tau1 = 0;
  ljet_trimmed_tau2 = 0;   
  ljet_trimmed_am_tau1 = 0;
  ljet_trimmed_am_tau2 = 0;   
  ljet_htt = 0;

  ljetBB_pt = 0;
  ljetBB_eta = 0;
  ljetBB_phi = 0;
  ljetBB_E = 0;
  ljetBB_split12 = 0;
  ljetBB_sub_n = 0;
  ljetBB_sub_pt = 0;
  ljetBB_sub_eta = 0;
  ljetBB_sub_phi = 0;
  ljetBB_sub_E = 0;
  ljetBB_sub_lcpt = 0;
  ljetBB_subunc_n = 0;
  ljetBB_subunc_pt = 0;
  ljetBB_subunc_eta = 0;
  ljetBB_subunc_phi = 0;
  ljetBB_subunc_E = 0;
  ljetBB_subunc_lcpt = 0;
  ljetBB_ug_tau1 = 0;
  ljetBB_ug_tau2 = 0;
  ljetBB_ug_am_tau1 = 0;
  ljetBB_ug_am_tau2 = 0;
  ljetBB_bdrs_m = 0;
  ljetBB_bdrs_pt = 0;
  ljetBB_bdrs_eta = 0;
  ljetBB_bdrs_phi = 0;
  ljetBB_ug_bdrs_tau1 = 0;
  ljetBB_ug_bdrs_tau2 = 0;
  ljetBB_htt = 0;



  partel_pt = 0;
  partel_eta = 0;
  partel_phi = 0;
  partel_E = 0;

  partmu_pt = 0;
  partmu_eta = 0;
  partmu_phi = 0;
  partmu_E = 0;

  partjet_pt = 0;
  partjet_eta = 0;
  partjet_phi = 0;
  partjet_E = 0;
  partjet_trueflav = 0;

  partmom_n = 0;
  partmom_pt = 0;
  partmom_eta = 0;
  partmom_phi = 0;
  partmom_E = 0;
  partmom_trueflav = 0;

  partljet_pt = 0;
  partljet_eta = 0;
  partljet_phi = 0;
  partljet_E = 0;
  partljet_split12 = 0;
  partljet_trueflav = 0;
  partljet_sub_n = 0;
  partljet_sub_pt = 0;
  partljet_sub_eta = 0;
  partljet_sub_phi = 0;
  partljet_sub_E = 0;

  partljetBB_pt = 0;
  partljetBB_eta = 0;
  partljetBB_phi = 0;
  partljetBB_E = 0;
  partljetBB_split12 = 0;
  partljetBB_sub_n = 0;
  partljetBB_sub_pt = 0;
  partljetBB_sub_eta = 0;
  partljetBB_sub_phi = 0;
  partljetBB_sub_E = 0;



  weight = 0;

  if (m_toWrite) {
    m_num->Branch("sumWeights", &m_sumWeights);
    m_num->Branch("sumWeights_var", &m_sumWeights_var);
    m_num->Branch("passCuts", &m_passCuts);

    m_chain->Branch("isTight", &isTight);

    //m_chain->Branch("hfor", &hfor);
    m_chain->Branch("mc_channel_number", &mc_channel_number);

    m_chain->Branch("systIdx", &systIdx);
    m_chain->Branch("passReco", &passReco);
    m_chain->Branch("passPart", &passPart);

    m_chain->Branch("triggerElectron", &triggerElectron);
    //m_chain->Branch("triggerElectron24", &triggerElectron24);
    //m_chain->Branch("triggerElectron60", &triggerElectron60);
    m_chain->Branch("triggerMuon", &triggerMuon);
    m_chain->Branch("triggerLargeJet", &triggerLargeJet);

    m_chain->Branch("weight", &weight);

    m_chain->Branch("el_n", &el_n);
    m_chain->Branch("el_pt", &el_pt);
    m_chain->Branch("el_eta", &el_eta);
    m_chain->Branch("el_phi", &el_phi);
    m_chain->Branch("el_E", &el_E);

    m_chain->Branch("mu_n", &mu_n);
    m_chain->Branch("mu_pt", &mu_pt);
    m_chain->Branch("mu_eta", &mu_eta);
    m_chain->Branch("mu_phi", &mu_phi);
    m_chain->Branch("mu_E", &mu_E);

    m_chain->Branch("jet_n", &jet_n);
    m_chain->Branch("jet_pt", &jet_pt);
    m_chain->Branch("jet_eta", &jet_eta);
    m_chain->Branch("jet_phi", &jet_phi);
    m_chain->Branch("jet_E", &jet_E);
    m_chain->Branch("jet_mv1", &jet_mv1);
    m_chain->Branch("jet_trueflav", &jet_trueflav);

    m_chain->Branch("trigjet_pt", &trigjet_pt);
    m_chain->Branch("trigjet_eta", &trigjet_eta);
    m_chain->Branch("trigjet_phi", &trigjet_phi);
    m_chain->Branch("trigjet_E", &trigjet_E);

    m_chain->Branch("ljet_n", &ljet_n);
    m_chain->Branch("ljet_pt", &ljet_pt);
    m_chain->Branch("ljet_eta", &ljet_eta);
    m_chain->Branch("ljet_phi", &ljet_phi);
    m_chain->Branch("ljet_E", &ljet_E);
    m_chain->Branch("ljet_split12", &ljet_split12);
    m_chain->Branch("ljet_trueflav", &ljet_trueflav);
    m_chain->Branch("ljet_sub_n", &ljet_sub_n);
    m_chain->Branch("ljet_sub_pt", &ljet_sub_pt);
    m_chain->Branch("ljet_sub_eta", &ljet_sub_eta);
    m_chain->Branch("ljet_sub_phi", &ljet_sub_phi);
    m_chain->Branch("ljet_sub_E", &ljet_sub_E);
    m_chain->Branch("ljet_sub_lcpt", &ljet_sub_lcpt);
    m_chain->Branch("ljet_subunc_n", &ljet_subunc_n);
    m_chain->Branch("ljet_subunc_pt", &ljet_subunc_pt);
    m_chain->Branch("ljet_subunc_eta", &ljet_subunc_eta);
    m_chain->Branch("ljet_subunc_phi", &ljet_subunc_phi);
    m_chain->Branch("ljet_subunc_E", &ljet_subunc_E);
    m_chain->Branch("ljet_subunc_lcpt", &ljet_subunc_lcpt);
    m_chain->Branch("ljet_trimmed_tau1", &ljet_trimmed_tau1);
    m_chain->Branch("ljet_trimmed_tau2", &ljet_trimmed_tau2);
    m_chain->Branch("ljet_trimmed_am_tau1", &ljet_trimmed_am_tau1);
    m_chain->Branch("ljet_trimmed_am_tau2", &ljet_trimmed_am_tau2);
    m_chain->Branch("ljet_htt", &ljet_htt);

    m_chain->Branch("ljetBB_n", &ljetBB_n);
    m_chain->Branch("ljetBB_pt", &ljetBB_pt);
    m_chain->Branch("ljetBB_eta", &ljetBB_eta);
    m_chain->Branch("ljetBB_phi", &ljetBB_phi);
    m_chain->Branch("ljetBB_E", &ljetBB_E);
    m_chain->Branch("ljetBB_split12", &ljetBB_split12);
    m_chain->Branch("ljetBB_sub_n", &ljetBB_sub_n);
    m_chain->Branch("ljetBB_sub_pt", &ljetBB_sub_pt);
    m_chain->Branch("ljetBB_sub_eta", &ljetBB_sub_eta);
    m_chain->Branch("ljetBB_sub_phi", &ljetBB_sub_phi);
    m_chain->Branch("ljetBB_sub_E", &ljetBB_sub_E);
    m_chain->Branch("ljetBB_sub_lcpt", &ljetBB_sub_lcpt);
    m_chain->Branch("ljetBB_subunc_n", &ljetBB_subunc_n);
    m_chain->Branch("ljetBB_subunc_pt", &ljetBB_subunc_pt);
    m_chain->Branch("ljetBB_subunc_eta", &ljetBB_subunc_eta);
    m_chain->Branch("ljetBB_subunc_phi", &ljetBB_subunc_phi);
    m_chain->Branch("ljetBB_subunc_E", &ljetBB_subunc_E);
    m_chain->Branch("ljetBB_subunc_lcpt", &ljetBB_subunc_lcpt);
    m_chain->Branch("ljetBB_ug_tau1", &ljetBB_ug_tau1);
    m_chain->Branch("ljetBB_ug_tau2", &ljetBB_ug_tau2);
    m_chain->Branch("ljetBB_ug_am_tau1", &ljetBB_ug_am_tau1);
    m_chain->Branch("ljetBB_ug_am_tau2", &ljetBB_ug_am_tau2);
    m_chain->Branch("ljetBB_bdrs_m", &ljetBB_bdrs_m);
    m_chain->Branch("ljetBB_bdrs_pt", &ljetBB_bdrs_pt);
    m_chain->Branch("ljetBB_bdrs_eta", &ljetBB_bdrs_eta);
    m_chain->Branch("ljetBB_bdrs_phi", &ljetBB_bdrs_phi);
    m_chain->Branch("ljetBB_ug_bdrs_tau1", &ljetBB_ug_bdrs_tau1);
    m_chain->Branch("ljetBB_ug_bdrs_tau2", &ljetBB_ug_bdrs_tau2);
    m_chain->Branch("ljetBB_htt", &ljetBB_htt);




    m_chain->Branch("met_etx", &met_etx);
    m_chain->Branch("met_ety", &met_ety);
  
    m_chain->Branch("partel_n", &partel_n);
    m_chain->Branch("partel_pt", &partel_pt);
    m_chain->Branch("partel_eta", &partel_eta);
    m_chain->Branch("partel_phi", &partel_phi);
    m_chain->Branch("partel_E", &partel_E);
 
    m_chain->Branch("partmu_n", &partmu_n);
    m_chain->Branch("partmu_pt", &partmu_pt);
    m_chain->Branch("partmu_eta", &partmu_eta);
    m_chain->Branch("partmu_phi", &partmu_phi);
    m_chain->Branch("partmu_E", &partmu_E);

    m_chain->Branch("partjet_n", &partjet_n);
    m_chain->Branch("partjet_pt", &partjet_pt);
    m_chain->Branch("partjet_eta", &partjet_eta);
    m_chain->Branch("partjet_phi", &partjet_phi);
    m_chain->Branch("partjet_E", &partjet_E);
    m_chain->Branch("partjet_trueflav", &partjet_trueflav);
  
    m_chain->Branch("partmom_n", &partmom_n);
    m_chain->Branch("partmom_pt", &partmom_pt);
    m_chain->Branch("partmom_eta", &partmom_eta);
    m_chain->Branch("partmom_phi", &partmom_phi);
    m_chain->Branch("partmom_E", &partmom_E);
    m_chain->Branch("partmom_trueflav", &partmom_trueflav);

    m_chain->Branch("partljet_n", &partljet_n);
    m_chain->Branch("partljet_pt", &partljet_pt);
    m_chain->Branch("partljet_eta", &partljet_eta);
    m_chain->Branch("partljet_phi", &partljet_phi);
    m_chain->Branch("partljet_E", &partljet_E);
    m_chain->Branch("partljet_split12", &partljet_split12);
    m_chain->Branch("partljet_trueflav", &partljet_trueflav);
    m_chain->Branch("partljet_sub_n", &partljet_sub_n);
    m_chain->Branch("partljet_sub_pt", &partljet_sub_pt);
    m_chain->Branch("partljet_sub_eta", &partljet_sub_eta);
    m_chain->Branch("partljet_sub_phi", &partljet_sub_phi);
    m_chain->Branch("partljet_sub_E", &partljet_sub_E);

    m_chain->Branch("partljetBB_n", &partljetBB_n);
    m_chain->Branch("partljetBB_pt", &partljetBB_pt);
    m_chain->Branch("partljetBB_eta", &partljetBB_eta);
    m_chain->Branch("partljetBB_phi", &partljetBB_phi);
    m_chain->Branch("partljetBB_E", &partljetBB_E);
    m_chain->Branch("partljetBB_split12", &partljetBB_split12);
    m_chain->Branch("partljetBB_sub_n", &partljetBB_sub_n);
    m_chain->Branch("partljetBB_sub_pt", &partljetBB_sub_pt);
    m_chain->Branch("partljetBB_sub_eta", &partljetBB_sub_eta);
    m_chain->Branch("partljetBB_sub_phi", &partljetBB_sub_phi);
    m_chain->Branch("partljetBB_sub_E", &partljetBB_sub_E);



    m_chain->Branch("partmet_etx", &partmet_etx);
    m_chain->Branch("partmet_ety", &partmet_ety);

    m_chain->Branch("avmu", &avmu);
    m_chain->Branch("npv", &npv);
  } else {
    m_num->SetBranchAddress("sumWeights", &m_sumWeights);
    m_num->SetBranchAddress("sumWeights_var", &m_sumWeights_var);
    m_num->SetBranchAddress("passCuts", &m_passCuts);

    m_chain->SetBranchAddress("isTight", &isTight);

    //m_chain->SetBranchAddress("hfor", &hfor);
    m_chain->SetBranchAddress("mc_channel_number", &mc_channel_number);

    m_chain->SetBranchAddress("systIdx", &systIdx);
    m_chain->SetBranchAddress("passReco", &passReco);
    m_chain->SetBranchAddress("passPart", &passPart);

    m_chain->SetBranchAddress("triggerElectron", &triggerElectron);
    //m_chain->SetBranchAddress("triggerElectron24", &triggerElectron24);
  //  m_chain->SetBranchAddress("triggerElectron60", &triggerElectron60);
    m_chain->SetBranchAddress("triggerMuon", &triggerMuon);
    m_chain->SetBranchAddress("triggerLargeJet", &triggerLargeJet);

    m_chain->SetBranchAddress("weight", &weight);

    m_chain->SetBranchAddress("el_n", &el_n);
    m_chain->SetBranchAddress("el_pt", &el_pt);
    m_chain->SetBranchAddress("el_eta", &el_eta);
    m_chain->SetBranchAddress("el_phi", &el_phi);
    m_chain->SetBranchAddress("el_E", &el_E);

    m_chain->SetBranchAddress("mu_n", &mu_n);
    m_chain->SetBranchAddress("mu_pt", &mu_pt);
    m_chain->SetBranchAddress("mu_eta", &mu_eta);
    m_chain->SetBranchAddress("mu_phi", &mu_phi);
    m_chain->SetBranchAddress("mu_E", &mu_E);

    m_chain->SetBranchAddress("jet_n", &jet_n);
    m_chain->SetBranchAddress("jet_pt", &jet_pt);
    m_chain->SetBranchAddress("jet_eta", &jet_eta);
    m_chain->SetBranchAddress("jet_phi", &jet_phi);
    m_chain->SetBranchAddress("jet_E", &jet_E);
    m_chain->SetBranchAddress("jet_mv1", &jet_mv1);
    m_chain->SetBranchAddress("jet_trueflav", &jet_trueflav);

    m_chain->SetBranchAddress("trigjet_pt", &trigjet_pt);
    m_chain->SetBranchAddress("trigjet_eta", &trigjet_eta);
    m_chain->SetBranchAddress("trigjet_phi", &trigjet_phi);
    m_chain->SetBranchAddress("trigjet_E", &trigjet_E);

    m_chain->SetBranchAddress("ljet_n", &ljet_n);
    m_chain->SetBranchAddress("ljet_pt", &ljet_pt);
    m_chain->SetBranchAddress("ljet_eta", &ljet_eta);
    m_chain->SetBranchAddress("ljet_phi", &ljet_phi);
    m_chain->SetBranchAddress("ljet_E", &ljet_E);
    m_chain->SetBranchAddress("ljet_split12", &ljet_split12);
    m_chain->SetBranchAddress("ljet_trueflav", &ljet_trueflav);
    m_chain->SetBranchAddress("ljet_sub_n", &ljet_sub_n);
    m_chain->SetBranchAddress("ljet_sub_pt", &ljet_sub_pt);
    m_chain->SetBranchAddress("ljet_sub_eta", &ljet_sub_eta);
    m_chain->SetBranchAddress("ljet_sub_phi", &ljet_sub_phi);
    m_chain->SetBranchAddress("ljet_sub_E", &ljet_sub_E);
    m_chain->SetBranchAddress("ljet_sub_lcpt", &ljet_sub_lcpt);
    m_chain->SetBranchAddress("ljet_subunc_n", &ljet_subunc_n);
    m_chain->SetBranchAddress("ljet_subunc_pt", &ljet_subunc_pt);
    m_chain->SetBranchAddress("ljet_subunc_eta", &ljet_subunc_eta);
    m_chain->SetBranchAddress("ljet_subunc_phi", &ljet_subunc_phi);
    m_chain->SetBranchAddress("ljet_subunc_E", &ljet_subunc_E);
    m_chain->SetBranchAddress("ljet_subunc_lcpt", &ljet_subunc_lcpt);
    m_chain->SetBranchAddress("ljet_trimmed_tau1", &ljet_trimmed_tau1);
    m_chain->SetBranchAddress("ljet_trimmed_tau2", &ljet_trimmed_tau2);
    m_chain->SetBranchAddress("ljet_trimmed_am_tau1", &ljet_trimmed_am_tau1);
    m_chain->SetBranchAddress("ljet_trimmed_am_tau2", &ljet_trimmed_am_tau2);
    m_chain->SetBranchAddress("ljet_htt", &ljet_htt);



    m_chain->SetBranchAddress("ljetBB_n", &ljetBB_n);
    m_chain->SetBranchAddress("ljetBB_pt", &ljetBB_pt);
    m_chain->SetBranchAddress("ljetBB_eta", &ljetBB_eta);
    m_chain->SetBranchAddress("ljetBB_phi", &ljetBB_phi);
    m_chain->SetBranchAddress("ljetBB_E", &ljetBB_E);
    m_chain->SetBranchAddress("ljetBB_split12", &ljetBB_split12);
    m_chain->SetBranchAddress("ljetBB_sub_n", &ljetBB_sub_n);
    m_chain->SetBranchAddress("ljetBB_sub_pt", &ljetBB_sub_pt);
    m_chain->SetBranchAddress("ljetBB_sub_eta", &ljetBB_sub_eta);
    m_chain->SetBranchAddress("ljetBB_sub_phi", &ljetBB_sub_phi);
    m_chain->SetBranchAddress("ljetBB_sub_E", &ljetBB_sub_E);
    m_chain->SetBranchAddress("ljetBB_sub_lcpt", &ljetBB_sub_lcpt);
    m_chain->SetBranchAddress("ljetBB_subunc_n", &ljetBB_subunc_n);
    m_chain->SetBranchAddress("ljetBB_subunc_pt", &ljetBB_subunc_pt);
    m_chain->SetBranchAddress("ljetBB_subunc_eta", &ljetBB_subunc_eta);
    m_chain->SetBranchAddress("ljetBB_subunc_phi", &ljetBB_subunc_phi);
    m_chain->SetBranchAddress("ljetBB_subunc_E", &ljetBB_subunc_E);
    m_chain->SetBranchAddress("ljetBB_subunc_lcpt", &ljetBB_subunc_lcpt);
    m_chain->SetBranchAddress("ljetBB_ug_tau1", &ljetBB_ug_tau1);
    m_chain->SetBranchAddress("ljetBB_ug_tau2", &ljetBB_ug_tau2);
    m_chain->SetBranchAddress("ljetBB_ug_am_tau1", &ljetBB_ug_am_tau1);
    m_chain->SetBranchAddress("ljetBB_ug_am_tau2", &ljetBB_ug_am_tau2);
    m_chain->SetBranchAddress("ljetBB_bdrs_m", &ljetBB_bdrs_m);
    m_chain->SetBranchAddress("ljetBB_bdrs_phi", &ljetBB_bdrs_phi);
    m_chain->SetBranchAddress("ljetBB_bdrs_eta", &ljetBB_bdrs_eta);
    m_chain->SetBranchAddress("ljetBB_bdrs_pt", &ljetBB_bdrs_pt);
    m_chain->SetBranchAddress("ljetBB_ug_bdrs_tau1", &ljetBB_ug_bdrs_tau1);
    m_chain->SetBranchAddress("ljetBB_ug_bdrs_tau2", &ljetBB_ug_bdrs_tau2);
    m_chain->SetBranchAddress("ljetBB_htt", &ljetBB_htt);


    m_chain->SetBranchAddress("met_etx", &met_etx);
    m_chain->SetBranchAddress("met_ety", &met_ety);
  
    m_chain->SetBranchAddress("partel_n", &partel_n);
    m_chain->SetBranchAddress("partel_pt", &partel_pt);
    m_chain->SetBranchAddress("partel_eta", &partel_eta);
    m_chain->SetBranchAddress("partel_phi", &partel_phi);
    m_chain->SetBranchAddress("partel_E", &partel_E);
 
    m_chain->SetBranchAddress("partmu_n", &partmu_n);
    m_chain->SetBranchAddress("partmu_pt", &partmu_pt);
    m_chain->SetBranchAddress("partmu_eta", &partmu_eta);
    m_chain->SetBranchAddress("partmu_phi", &partmu_phi);
    m_chain->SetBranchAddress("partmu_E", &partmu_E);

    m_chain->SetBranchAddress("partjet_n", &partjet_n);
    m_chain->SetBranchAddress("partjet_pt", &partjet_pt);
    m_chain->SetBranchAddress("partjet_eta", &partjet_eta);
    m_chain->SetBranchAddress("partjet_phi", &partjet_phi);
    m_chain->SetBranchAddress("partjet_E", &partjet_E);
    m_chain->SetBranchAddress("partjet_trueflav", &partjet_trueflav);

    m_chain->SetBranchAddress("partmom_n", &partmom_n);
    m_chain->SetBranchAddress("partmom_pt", &partmom_pt);
    m_chain->SetBranchAddress("partmom_eta", &partmom_eta);
    m_chain->SetBranchAddress("partmom_phi", &partmom_phi);
    m_chain->SetBranchAddress("partmom_E", &partmom_E);
    m_chain->SetBranchAddress("partmom_trueflav", &partmom_trueflav);
  
    m_chain->SetBranchAddress("partljet_n", &partljet_n);
    m_chain->SetBranchAddress("partljet_pt", &partljet_pt);
    m_chain->SetBranchAddress("partljet_eta", &partljet_eta);
    m_chain->SetBranchAddress("partljet_phi", &partljet_phi);
    m_chain->SetBranchAddress("partljet_E", &partljet_E);
    m_chain->SetBranchAddress("partljet_split12", &partljet_split12);
    m_chain->SetBranchAddress("partljet_trueflav", &partljet_trueflav);
    m_chain->SetBranchAddress("partljet_sub_n", &partljet_sub_n);
    m_chain->SetBranchAddress("partljet_sub_pt", &partljet_sub_pt);
    m_chain->SetBranchAddress("partljet_sub_eta", &partljet_sub_eta);
    m_chain->SetBranchAddress("partljet_sub_phi", &partljet_sub_phi);
    m_chain->SetBranchAddress("partljet_sub_E", &partljet_sub_E);

    m_chain->SetBranchAddress("partljetBB_n", &partljetBB_n);
    m_chain->SetBranchAddress("partljetBB_pt", &partljetBB_pt);
    m_chain->SetBranchAddress("partljetBB_eta", &partljetBB_eta);
    m_chain->SetBranchAddress("partljetBB_phi", &partljetBB_phi);
    m_chain->SetBranchAddress("partljetBB_E", &partljetBB_E);
    m_chain->SetBranchAddress("partljetBB_split12", &partljetBB_split12);
    m_chain->SetBranchAddress("partljetBB_sub_n", &partljetBB_sub_n);
    m_chain->SetBranchAddress("partljetBB_sub_pt", &partljetBB_sub_pt);
    m_chain->SetBranchAddress("partljetBB_sub_eta", &partljetBB_sub_eta);
    m_chain->SetBranchAddress("partljetBB_sub_phi", &partljetBB_sub_phi);
    m_chain->SetBranchAddress("partljetBB_sub_E", &partljetBB_sub_E);


    m_chain->SetBranchAddress("partmet_etx", &partmet_etx);
    m_chain->SetBranchAddress("partmet_ety", &partmet_ety);

    m_chain->SetBranchAddress("avmu", &avmu);
    m_chain->SetBranchAddress("npv", &npv);
  }
}

