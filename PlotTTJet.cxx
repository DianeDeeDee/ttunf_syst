#include "Plot.h"
#include "PlotTTJet.h"
#include "Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "HistogramService.h"

#include <fastjet/PseudoJet.hh>

#include <cmath>

#include "SD-tosvn/Exception.h"

using namespace fastjet;
using namespace Deconstruction;

PlotTTJet::PlotTTJet(const std::string &filename, const std::vector<std::string> &systs)
  : Plot(filename, systs),
    signal(0), background(0), isr(0), deconstruct(0), param("input_card.dat"), paramb("input_card.dat") {
  m_hSvc.create1D("callargeJetPt", "; Large jet p_{T} [GeV] ; Events", 28, 200, 1600);
  m_hSvc.create1D("callargeJetEta", "; Large jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("callargeJetPhi", "; Large jet #phi [rad] ; Events", 64, -3.2, 3.2);
  m_hSvc.create1D("callargeJetM", "; Large jet M ; Events", 40, 0, 400);

  m_hSvc.create1D("llargeJetPt", "; Large jet p_{T} [GeV] ; Events", 28, 200, 1600);
  m_hSvc.create1D("llargeJetEta", "; Large jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("llargeJetPhi", "; Large jet #phi [rad] ; Events", 64, -3.2, 3.2);
  m_hSvc.create1D("llargeJetM", "; Large jet M ; Events", 40, 0, 400);

  m_hSvc.create1D("lsubjetN", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 11, -0.5, 10.5);
  m_hSvc.create1D("lsubjet0Pt", "; Anti-k_{t} R=1.0 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("lsubjet1Pt", "; Anti-k_{t} R=1.0 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet2Pt", "; Anti-k_{t} R=1.0 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet3Pt", "; Anti-k_{t} R=1.0 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet4Pt", "; Anti-k_{t} R=1.0 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);

  m_hSvc.create1D("lsubjet0Eta", "; Anti-k_{t} R=1.0 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet1Eta", "; Anti-k_{t} R=1.0 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet2Eta", "; Anti-k_{t} R=1.0 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet3Eta", "; Anti-k_{t} R=1.0 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet4Eta", "; Anti-k_{t} R=1.0 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  m_hSvc.create1D("calsubjetN", "; C/A R=1.2 Sub-jet multiplicity ; Events", 11, -0.5, 10.5);
  m_hSvc.create1D("calsubjet0Pt", "; C/A R=1.2 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("calsubjet1Pt", "; C/A R=1.2 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet2Pt", "; C/A R=1.2 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet3Pt", "; C/A R=1.2 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet4Pt", "; C/A R=1.2 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);

  m_hSvc.create1D("calsubjet0Eta", "; C/A R=1.2 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet1Eta", "; C/A R=1.2 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet2Eta", "; C/A R=1.2 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet3Eta", "; C/A R=1.2 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet4Eta", "; C/A R=1.2 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  m_hSvc.create1D("lsubjetNChi40", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 11, -0.5, 10.5);
  m_hSvc.create1D("lsubjet0PtChi40", "; Anti-k_{t} R=1.0 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("lsubjet1PtChi40", "; Anti-k_{t} R=1.0 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet2PtChi40", "; Anti-k_{t} R=1.0 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet3PtChi40", "; Anti-k_{t} R=1.0 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet4PtChi40", "; Anti-k_{t} R=1.0 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);

  m_hSvc.create1D("calsubjetNChi40", "; C/A R=1.2 Sub-jet multiplicity ; Events", 11, -0.5, 10.5);
  m_hSvc.create1D("calsubjet0PtChi40", "; C/A R=1.2 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("calsubjet1PtChi40", "; C/A R=1.2 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet2PtChi40", "; C/A R=1.2 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet3PtChi40", "; C/A R=1.2 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet4PtChi40", "; C/A R=1.2 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);

  m_hSvc.create1D("lsubjetNChib40", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 11, -0.5, 10.5);
  m_hSvc.create1D("lsubjet0PtChib40", "; Anti-k_{t} R=1.0 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("lsubjet1PtChib40", "; Anti-k_{t} R=1.0 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet2PtChib40", "; Anti-k_{t} R=1.0 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet3PtChib40", "; Anti-k_{t} R=1.0 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet4PtChib40", "; Anti-k_{t} R=1.0 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);

  m_hSvc.create1D("calsubjetNChib40", "; C/A R=1.2 Sub-jet multiplicity ; Events", 11, -0.5, 10.5);
  m_hSvc.create1D("calsubjet0PtChib40", "; C/A R=1.2 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("calsubjet1PtChib40", "; C/A R=1.2 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet2PtChib40", "; C/A R=1.2 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet3PtChib40", "; C/A R=1.2 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet4PtChib40", "; C/A R=1.2 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);

  m_hSvc.create1D("lsubjet0EtaChi40", "; Anti-k_{t} R=1.0 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet1EtaChi40", "; Anti-k_{t} R=1.0 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet2EtaChi40", "; Anti-k_{t} R=1.0 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet3EtaChi40", "; Anti-k_{t} R=1.0 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet4EtaChi40", "; Anti-k_{t} R=1.0 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  m_hSvc.create1D("calsubjet0EtaChi40", "; C/A R=1.2 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet1EtaChi40", "; C/A R=1.2 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet2EtaChi40", "; C/A R=1.2 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet3EtaChi40", "; C/A R=1.2 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet4EtaChi40", "; C/A R=1.2 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  m_hSvc.create1D("lsubjet0EtaChib40", "; Anti-k_{t} R=1.0 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet1EtaChib40", "; Anti-k_{t} R=1.0 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet2EtaChib40", "; Anti-k_{t} R=1.0 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet3EtaChib40", "; Anti-k_{t} R=1.0 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet4EtaChib40", "; Anti-k_{t} R=1.0 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  m_hSvc.create1D("calsubjet0EtaChib40", "; C/A R=1.2 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet1EtaChib40", "; C/A R=1.2 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet2EtaChib40", "; C/A R=1.2 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet3EtaChib40", "; C/A R=1.2 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet4EtaChib40", "; C/A R=1.2 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  m_hSvc.create1D("lchi", "; SD log(#chi) ; Events", 30, -15, 15);
  m_hSvc.create1D("lsd_sig", "; log(SD signal weight) ; Events", 25, -40, -15);
  m_hSvc.create1D("lsd_bkg", "; log(SD background weight) ; Events", 25, -40, -15);

  m_hSvc.create1D("calchi", "; SD log(#chi) ; Events", 30, -15, 15);
  m_hSvc.create1D("calsd_sig", "; log(SD signal weight) ; Events", 25, -40, -15);
  m_hSvc.create1D("calsd_bkg", "; log(SD background weight) ; Events", 25, -40, -15);

  m_hSvc.create1D("lchib", "; SD log(#chi) ; Events", 30, -15, 15);
  m_hSvc.create1D("lsd_sigb", "; log(SD signal weight) ; Events", 25, -40, -15);
  m_hSvc.create1D("lsd_bkgb", "; log(SD background weight) ; Events", 25, -40, -15);

  m_hSvc.create1D("calchib", "; SD log(#chi) ; Events", 30, -15, 15);
  m_hSvc.create1D("calsd_sigb", "; log(SD signal weight) ; Events", 25, -40, -15);
  m_hSvc.create1D("calsd_bkgb", "; log(SD background weight) ; Events", 25, -40, -15);

  m_hSvc.create1D("lchi_pt550", "; SD log(#chi) ; Events", 30, -15, 15);
  m_hSvc.create1D("lchi_350pt550", "; SD log(#chi) ; Events", 30, -15, 15);

  m_hSvc.create1D("lchiM", "; SD log(#chi) ; Events", 30, -15, 15);
  m_hSvc.create1D("lchiMD12", "; SD log(#chi) ; Events", 30, -15, 15);

  m_hSvc.create1D("mu", "; <#mu> ; Events", 40, 0, 40);
  m_hSvc.create1D("npv", "; npv ; Events", 30, 0, 30);

  double llb[] = {350, 400, 450, 500, 1200};
  m_hSvc.create1DVar("llargeJetPtEff", "; Large jet p_{T} [GeV] ; Events", 4, llb);
  m_hSvc.create1DVar("llargeJetPtEffCa", "; Large jet p_{T} [GeV] ; Events", 4, llb);
  m_hSvc.create1DVar("llargeJetPtEffChi40", "; Large jet p_{T} [GeV] ; Events", 4, llb);
  m_hSvc.create1DVar("llargeJetPtEffChib40", "; Large jet p_{T} [GeV] ; Events", 4, llb);
  m_hSvc.create1DVar("llargeJetPtEffCaChi40", "; Large jet p_{T} [GeV] ; Events", 4, llb);
  m_hSvc.create1DVar("llargeJetPtEffCaChib40", "; Large jet p_{T} [GeV] ; Events", 4, llb);

  m_hSvc.create1D("llargeJetPtPos", "; Large jet p_{T} [GeV] ; Events", 8, 350, 1150);
  m_hSvc.create1D("llargeJetPtPosCa", "; Large jet p_{T} [GeV] ; Events", 8, 350, 1150);
  m_hSvc.create1D("llargeJetPtPosChi40", "; Large jet p_{T} [GeV] ; Events", 8, 350, 1150);
  m_hSvc.create1D("llargeJetPtPosChib40", "; Large jet p_{T} [GeV] ; Events", 8, 350, 1150);
  m_hSvc.create1D("llargeJetPtPosCaChi40", "; Large jet p_{T} [GeV] ; Events", 8, 350, 1150);
  m_hSvc.create1D("llargeJetPtPosCaChib40", "; Large jet p_{T} [GeV] ; Events", 8, 350, 1150);

  m_hSvc.create1D("llargeJetMPos", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosCa", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosChi40", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosCaChi40", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosChib40", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosCaChib40", "; Large jet M ; Events", 40, 0, 400);

  m_hSvc.create1D("llargeJetSDMChi40", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMCaChi40", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMChib40", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMCaChib40", "; Large jet M ; Events", 40, 0, 400);

  signal = new TopGluonModel(param);
  background = new BackgroundModel(param);
  isr = new ISRModel(param);
  deconstruct = new Deconstruct(param, *signal, *background, *isr);

  paramb["useBtag"] = 1;

  signalb = new TopGluonModel(paramb);
  backgroundb = new BackgroundModel(paramb);
  isrb = new ISRModel(paramb);
  deconstructb = new Deconstruct(paramb, *signalb, *backgroundb, *isrb);

  m_hSvc.m_tree->Branch("syst", &m_tree_syst);

  m_hSvc.m_tree->Branch("nsub", &m_tree_nsub);
  m_hSvc.m_tree->Branch("chi", &m_tree_chi);
  m_hSvc.m_tree->Branch("chib", &m_tree_chib);
  m_hSvc.m_tree->Branch("bmatch", &m_tree_bmatch);
  m_hSvc.m_tree->Branch("pt", &m_tree_pt);
  m_hSvc.m_tree->Branch("eta", &m_tree_eta);
  m_hSvc.m_tree->Branch("truept", &m_tree_truept);
  m_hSvc.m_tree->Branch("trueeta", &m_tree_trueeta);
  m_hSvc.m_tree->Branch("m", &m_tree_m);
  m_hSvc.m_tree->Branch("d12", &m_tree_d12);

  m_hSvc.m_tree->Branch("cansub", &m_tree_cansub);
  m_hSvc.m_tree->Branch("cachi", &m_tree_cachi);
  m_hSvc.m_tree->Branch("cachib", &m_tree_cachib);
  m_hSvc.m_tree->Branch("cabmatch", &m_tree_cabmatch);
  m_hSvc.m_tree->Branch("capt", &m_tree_capt);
  m_hSvc.m_tree->Branch("caeta", &m_tree_caeta);
  m_hSvc.m_tree->Branch("catruept", &m_tree_catruept);
  m_hSvc.m_tree->Branch("catrueeta", &m_tree_catrueeta);

  m_hSvc.m_tree->Branch("cahtt", &m_tree_cahtt);

  m_hSvc.m_tree->Branch("npv", &m_tree_npv);
  m_hSvc.m_tree->Branch("mu", &m_tree_mu);
  m_hSvc.m_tree->Branch("w", &m_tree_w);
}

PlotTTJet::~PlotTTJet() {
  if (deconstruct) delete deconstruct;
  if (isr) delete isr;
  if (background) delete background;
  if (signal) delete signal;
}


std::vector<LargeJet::Subjet> PlotTTJet::applySjUnc(const std::vector<LargeJet::Subjet> &vsj, const std::vector<LargeJet::Subjet> &vsj_unc, const std::string &s) {
  bool doJer = false;
  bool doJes = false;
  bool isDown = false;
  if (s.find("_subJes") != std::string::npos) {
    doJes = true;
    if (s.find("Down") != std::string::npos) isDown = true;
    size_t pref = std::string("_subJes").size();
  } else if (s.find("_subJer") != std::string::npos) {
    doJer = true;
  }

  std::vector<LargeJet::Subjet> msj;
  for (int k = 0; k < vsj.size(); ++k) {
    LargeJet::Subjet sj = vsj[k];

    if (doJer) {
      sj = vsj_unc[k];
    } else if (doJes) {
      float unc = 0;
      if (std::fabs(sj.v.Eta()) < 0.8) {
        if (sj.v.Perp() < 30e3) {
          unc = 5.7e-2;
        } else if (sj.v.Perp() > 30e3 && sj.v.Perp() < 500e3) {
          unc = 2.8e-2;
        } else if (sj.v.Perp() > 500e3 && sj.v.Perp() < 600e3) {
          unc = 2.9e-2;
        } else if (sj.v.Perp() > 600e3 && sj.v.Perp() < 800e3) {
          unc = 3.3e-2;
        } else if (sj.v.Perp() > 800e3) {
          unc = 4.6e-2;
        }
      } else if (std::fabs(sj.v.Eta()) < 1.4) {
        if (sj.v.Perp() < 30e3) {
          unc = 6.2e-2;
        } else if (sj.v.Perp() > 30e3 && sj.v.Perp() < 500e3) {
          unc = 3.6e-2;
        } else if (sj.v.Perp() > 500e3 && sj.v.Perp() < 600e3) {
          unc = 3.7e-2;
        } else if (sj.v.Perp() > 600e3 && sj.v.Perp() < 800e3) {
          unc = 4.1e-2;
        } else if (sj.v.Perp() > 800e3) {
          unc = 5.2e-2;
        }
      } else if (std::fabs(sj.v.Eta()) < 2.1) {
        if (sj.v.Perp() < 30e3) {
          unc = 7.8e-2;
        } else if (sj.v.Perp() > 30e3 && sj.v.Perp() < 500e3) {
          unc = 5.9e-2;
        } else if (sj.v.Perp() > 500e3 && sj.v.Perp() < 600e3) {
          unc = 6.0e-2;
        } else if (sj.v.Perp() > 600e3 && sj.v.Perp() < 800e3) {
          unc = 6.2e-2;
        } else if (sj.v.Perp() > 800e3) {
          unc = 7.0e-2;
        }
      }
      if (isDown) unc = -unc;
      sj.v *= (1 + unc);
    }
    if (sj.v.Perp() < 20e3 || std::fabs(sj.v.Eta()) > 2.1) continue;
    msj.push_back(sj);
  }
  return msj;
}

void PlotTTJet::run(const Event &e, double weight, double pweight, const std::string &s) {
  HistogramService *h = &m_hSvc;
  if (e.passReco()) {
    // only run SD for large jets passing the following cuts for basic studies
    // pT > 350 GeV makes top be contained in the large jet
    int fji = -1;
    for (int k = 0; k < e.largeJet().size(); ++k) {
      if (e.largeJet()[k].passLoose() && e.largeJet()[k].mom().Perp() > 350e3) {
        fji = k;
        break;
      }
    }
    int cafji = -1;
    for (int k = 0; k < e.largeJetBB().size(); ++k) {
      if (e.largeJetBB()[k].passLoose() && e.largeJetBB()[k].mom().Perp() > 350e3) {
        cafji = k;
        break;
      }
    }
    if (fji != -1 && cafji != -1) {
      const LargeJet &lj = e.largeJet()[fji];
      const LargeJet &calj = e.largeJetBB()[cafji];
      if (lj.mom().DeltaR(calj.mom()) > 0.75) return;

      float mdr = lj.mom().DeltaR(e.triggerJet4mom());
      if (mdr < 1.5) {
        return;
      }

      h->h1D("mu", "", s)->Fill(e.mu(), weight);
      h->h1D("npv", "", s)->Fill(e.npv(), weight);

      h->h1D("llargeJetM", "", s)->Fill(lj.mom().M()*1e-3, weight);
      h->h1D("llargeJetPt", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
      h->h1D("llargeJetEta", "", s)->Fill(lj.mom().Eta(), weight);
      h->h1D("llargeJetPhi", "", s)->Fill(lj.mom().Phi(), weight);

      h->h1D("callargeJetM", "", s)->Fill(calj.mom().M()*1e-3, weight);
      h->h1D("callargeJetPt", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
      h->h1D("callargeJetEta", "", s)->Fill(calj.mom().Eta(), weight);
      h->h1D("callargeJetPhi", "", s)->Fill(calj.mom().Phi(), weight);

      int sjn = lj.subjet().size();
      std::vector<LargeJet::Subjet> vsj = applySjUnc(lj.subjet(), lj.subjetJER(), s);
      sjn = vsj.size();

      int casjn = calj.subjet().size();
      std::vector<LargeJet::Subjet> cavsj = applySjUnc(calj.subjet(), lj.subjetJER(), s);
      casjn = cavsj.size();
       
      h->h1D("lsubjetN", "", s)->Fill(sjn, weight);
      if (sjn >= 1)
         h->h1D("lsubjet0Pt", "", s)->Fill(vsj[0].v.Perp()*1e-3, weight);
      if (sjn >= 2)
        h->h1D("lsubjet1Pt", "", s)->Fill(vsj[1].v.Perp()*1e-3, weight);
      if (sjn >= 3)
        h->h1D("lsubjet2Pt", "", s)->Fill(vsj[2].v.Perp()*1e-3, weight);
      if (sjn >= 4)
        h->h1D("lsubjet3Pt", "", s)->Fill(vsj[3].v.Perp()*1e-3, weight);
      if (sjn >= 5)
        h->h1D("lsubjet4Pt", "", s)->Fill(vsj[4].v.Perp()*1e-3, weight);

      if (sjn >= 1)
        h->h1D("lsubjet0Eta", "", s)->Fill(vsj[0].v.Eta(), weight);
      if (sjn >= 2)
        h->h1D("lsubjet1Eta", "", s)->Fill(vsj[1].v.Eta(), weight);
      if (sjn >= 3)
        h->h1D("lsubjet2Eta", "", s)->Fill(vsj[2].v.Eta(), weight);
      if (sjn >= 4)
        h->h1D("lsubjet3Eta", "", s)->Fill(vsj[3].v.Eta(), weight);
      if (sjn >= 5)
        h->h1D("lsubjet4Eta", "", s)->Fill(vsj[4].v.Eta(), weight);

      h->h1D("calsubjetN", "", s)->Fill(casjn, weight);
      if (casjn >= 1)
        h->h1D("calsubjet0Pt", "", s)->Fill(cavsj[0].v.Perp()*1e-3, weight);
      if (casjn >= 2)
        h->h1D("calsubjet1Pt", "", s)->Fill(cavsj[1].v.Perp()*1e-3, weight);
      if (casjn >= 3)
        h->h1D("calsubjet2Pt", "", s)->Fill(cavsj[2].v.Perp()*1e-3, weight);
      if (casjn >= 4)
        h->h1D("calsubjet3Pt", "", s)->Fill(cavsj[3].v.Perp()*1e-3, weight);
      if (casjn >= 5)
        h->h1D("calsubjet4Pt", "", s)->Fill(cavsj[4].v.Perp()*1e-3, weight);
      if (casjn >= 1)
        h->h1D("calsubjet0Eta", "", s)->Fill(cavsj[0].v.Eta(), weight);
      if (casjn >= 2)
        h->h1D("calsubjet1Eta", "", s)->Fill(cavsj[1].v.Eta(), weight);
      if (casjn >= 3)
        h->h1D("calsubjet2Eta", "", s)->Fill(cavsj[2].v.Eta(), weight);
      if (casjn >= 4)
        h->h1D("calsubjet3Eta", "", s)->Fill(cavsj[3].v.Eta(), weight);
      if (casjn >= 5)
        h->h1D("calsubjet4Eta", "", s)->Fill(cavsj[4].v.Eta(), weight);

      std::vector<fastjet::PseudoJet> sd_subs;
      for (int k = 0; k < sjn; ++k) {
        const TLorentzVector &sj = vsj[k].v;
        fastjet::PseudoJet psj(sj.Px()*1e-3, sj.Py()*1e-3, sj.Pz()*1e-3, sj.E()*1e-3);
        psj.set_user_index(0);
        for (int z = 0; z < e.jet().size(); ++z) {
          double dr = e.jet()[z].mom().DeltaR(sj);
          if (dr < 0.2) {
            if (e.jet()[z].btag())
              psj.set_user_index(1);
            else
              psj.set_user_index(-1);
            break;
          }
        }
        sd_subs.push_back(psj);
      }
      std::vector<fastjet::PseudoJet> sd_sorted = sorted_by_pt(sd_subs);
      if (sd_sorted.size() > 9)
        sd_sorted.erase(sd_sorted.begin() + 9, sd_sorted.begin() + sd_sorted.size());

      std::vector<fastjet::PseudoJet> casd_subs;
      for (int k = 0; k < casjn; ++k) {
        const TLorentzVector &sj = cavsj[k].v;
        fastjet::PseudoJet psj(sj.Px()*1e-3, sj.Py()*1e-3, sj.Pz()*1e-3, sj.E()*1e-3);
        psj.set_user_index(0);
        for (int z = 0; z < e.jet().size(); ++z) {
          double dr = e.jet()[z].mom().DeltaR(sj);
          if (dr < 0.2) {
            if (e.jet()[z].btag())
              psj.set_user_index(1);
            else
              psj.set_user_index(-1);
            break;
          }
        }
        casd_subs.push_back(psj);
      }
      std::vector<fastjet::PseudoJet> casd_sorted = sorted_by_pt(casd_subs);
      if (casd_sorted.size() > 9)
        casd_sorted.erase(casd_sorted.begin() + 9, casd_sorted.begin() + casd_sorted.size());

      double chi = -100;
      double wSig = 0;
      double wBkg = 0;
      fastjet::PseudoJet chi_top(0,0,0,0);
      double chi_signal = 0;
      std::multimap<double, std::vector<fastjet::PseudoJet> > signalWeight;
      try {
        chi = deconstruct->deconstruct(sd_sorted, wSig, wBkg);
        signalWeight = deconstruct->signalWeight();
        for (std::multimap<double, std::vector<fastjet::PseudoJet> >::const_iterator it = signalWeight.begin(); it != signalWeight.end(); ++it) {
          fastjet::PseudoJet thisTop(0,0,0,0);
          for (int z = 0; z < it->second.size(); ++z) thisTop += it->second[z];
          chi_top += it->first*thisTop;
          chi_signal += it->first;
        }
        chi_top /= chi_signal;
      } catch (Deconstruction::Exception &e) {
        std::cout << "Exception while running SD: " << e.what() << std::endl;
      }
      float lchi = -100;
      if (chi > 0) lchi = std::log(chi);
      float lwSig = -100;
      if (wSig > 0) lwSig = std::log(wSig);
      float lwBkg = -100;
      if (wBkg > 0) lwBkg = std::log(wBkg);

      double cachi = -100;
      double cawSig = 0;
      double cawBkg = 0;
      double cachi_signal = 0;
      fastjet::PseudoJet cachi_top(0,0,0,0);
      std::multimap<double, std::vector<fastjet::PseudoJet> > casignalWeight;
      try {
        cachi = deconstruct->deconstruct(casd_sorted, cawSig, cawBkg);
        casignalWeight = deconstruct->signalWeight();
        for (std::multimap<double, std::vector<fastjet::PseudoJet> >::const_iterator it = casignalWeight.begin(); it != casignalWeight.end(); ++it) {
          fastjet::PseudoJet thisTop(0,0,0,0);
          for (int z = 0; z < it->second.size(); ++z) thisTop += it->second[z];
          cachi_top += it->first*thisTop;
          cachi_signal += it->first;
        }
        cachi_top /= cachi_signal;
      } catch (Deconstruction::Exception &e) {
        std::cout << "Exception while running SD: " << e.what() << std::endl;
      }
      float calchi = -100;
      if (cachi > 0) calchi = std::log(cachi);
      float calwSig = -100;
      if (cawSig > 0) calwSig = std::log(cawSig);
      float calwBkg = -100;
      if (cawBkg > 0) calwBkg = std::log(cawBkg);

      double cachib = -100;
      double cawSigb = 0;
      double cawBkgb = 0;
      fastjet::PseudoJet cachib_top(0,0,0,0);
      double cachib_signal = 0;
      std::multimap<double, std::vector<fastjet::PseudoJet> > casignalWeightb;
      try {
        cachib = deconstructb->deconstruct(casd_sorted, cawSigb, cawBkgb);
        casignalWeightb = deconstructb->signalWeight();
        for (std::multimap<double, std::vector<fastjet::PseudoJet> >::const_iterator it = casignalWeightb.begin(); it != casignalWeightb.end(); ++it) {
          fastjet::PseudoJet thisTop(0,0,0,0);
          for (int z = 0; z < it->second.size(); ++z) thisTop += it->second[z];
          cachib_top += it->first*thisTop;
          cachib_signal +=  it->first;
        }
        cachib_top /= cachib_signal;
      } catch (Deconstruction::Exception &e) {
        std::cout << "Exception while running SD: " << e.what() << std::endl;
      }
      float calchib = -100;
      if (cachib > 0) calchib = std::log(cachib);
      float calwSigb = -100;
      if (cawSigb > 0) calwSigb = std::log(cawSigb);
      float calwBkgb = -100;
      if (cawBkgb > 0) calwBkgb = std::log(cawBkgb);

      double chib = -100;
      double wSigb = 0;
      double wBkgb = 0;
      fastjet::PseudoJet chib_top(0,0,0,0);
      double chib_signal = 0;
      std::multimap<double, std::vector<fastjet::PseudoJet> > signalWeightb;
      try {
        chib = deconstructb->deconstruct(sd_sorted, wSigb, wBkgb);
        signalWeightb = deconstructb->signalWeight();
        for (std::multimap<double, std::vector<fastjet::PseudoJet> >::const_iterator it = signalWeightb.begin(); it != signalWeightb.end(); ++it) {
          fastjet::PseudoJet thisTop(0,0,0,0);
          for (int z = 0; z < it->second.size(); ++z) thisTop += it->second[z];
          chib_top += it->first*thisTop;
          chib_signal +=  it->first;
        }
        chib_top /= chib_signal;
      } catch (Deconstruction::Exception &e) {
        std::cout << "Exception while running SD: " << e.what() << std::endl;
      }
      float lchib = -100;
      if (chib > 0) lchib = std::log(chib);
      float lwSigb = -100;
      if (wSigb > 0) lwSigb = std::log(wSig);
      float lwBkgb = -100;
      if (wBkgb > 0) lwBkgb = std::log(wBkg);

      h->h1D("lchi", "", s)->Fill(lchi, weight);
      h->h1D("lsd_sig", "", s)->Fill(lwSig, weight);
      h->h1D("lsd_bkg", "", s)->Fill(lwBkg, weight);

      h->h1D("calchi", "", s)->Fill(calchi, weight);
      h->h1D("calsd_sig", "", s)->Fill(calwSig, weight);
      h->h1D("calsd_bkg", "", s)->Fill(calwBkg, weight);

      h->h1D("lchib", "", s)->Fill(lchib, weight);
      h->h1D("lsd_sigb", "", s)->Fill(lwSigb, weight);
      h->h1D("lsd_bkgb", "", s)->Fill(lwBkgb, weight);

      h->h1D("calchib", "", s)->Fill(calchib, weight);
      h->h1D("calsd_sigb", "", s)->Fill(calwSigb, weight);
      h->h1D("calsd_bkgb", "", s)->Fill(calwBkgb, weight);

      if (e.largeJet()[0].mom().Perp() >= 550e3) {
        h->h1D("lchi_pt550", "", s)->Fill(lchi, weight);
      } else {
        h->h1D("lchi_350pt550", "", s)->Fill(lchi, weight);
      }

      if (e.largeJet()[0].mom().M() > 100e3)
        h->h1D("lchiM", "", s)->Fill(lchi, weight);
      if (e.largeJet()[0].mom().M() > 100e3 && e.largeJet()[0].split12() > 40e3)
        h->h1D("lchiMD12", "", s)->Fill(lchi, weight);

      h->h1D("llargeJetPtEff", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
      h->h1D("llargeJetPtEffCa", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
      if (lchi > 2.7) h->h1D("llargeJetPtEffChi40", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
      if (calchi > 2.8) h->h1D("llargeJetPtEffCaChi40", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
      if (lchib > 2.6) h->h1D("llargeJetPtEffChib40", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
      if (calchib > 2.8) h->h1D("llargeJetPtEffCaChib40", "", s)->Fill(calj.mom().Perp()*1e-3, weight);

      h->h1D("llargeJetPtPos", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
      h->h1D("llargeJetPtPosCa", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
      if (lchi > 2.7) h->h1D("llargeJetPtPosChi40", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
      if (calchi > 2.8) h->h1D("llargeJetPtPosCaChi40", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
      if (lchib > 2.6) h->h1D("llargeJetPtPosChib40", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
      if (calchib > 2.8) h->h1D("llargeJetPtPosCaChib40", "", s)->Fill(calj.mom().Perp()*1e-3, weight);

      h->h1D("llargeJetMPos", "", s)->Fill(lj.mom().M()*1e-3, weight);
      h->h1D("llargeJetMPosCa", "", s)->Fill(calj.mom().M()*1e-3, weight);
      if (lchi > 2.7) h->h1D("llargeJetMPosChi40", "", s)->Fill(lj.mom().M()*1e-3, weight);
      if (calchi > 2.8) h->h1D("llargeJetMPosCaChi40", "", s)->Fill(calj.mom().M()*1e-3, weight);
      if (lchib > 2.6) h->h1D("llargeJetMPosChib40", "", s)->Fill(lj.mom().M()*1e-3, weight);
      if (calchib > 2.8) h->h1D("llargeJetMPosCaChib40", "", s)->Fill(calj.mom().M()*1e-3, weight);

      if (lchi > 2.7) h->h1D("llargeJetSDMChi40", "", s)->Fill(chi_top.m(), weight);
      if (calchi > 2.8) h->h1D("llargeJetSDMCaChi40", "", s)->Fill(cachi_top.m(), weight);
      if (lchib > 2.6) h->h1D("llargeJetSDMChib40", "", s)->Fill(chib_top.m(), weight);
      if (calchib > 2.8) h->h1D("llargeJetSDMCaChib40", "", s)->Fill(cachib_top.m(), weight);

      if (lchi > 2.7) {
        h->h1D("lsubjetNChi40", "", s)->Fill(sjn, weight);
        if (sjn >= 1)
          h->h1D("lsubjet0PtChi40", "", s)->Fill(vsj[0].v.Perp()*1e-3, weight);
        if (sjn >= 2)
          h->h1D("lsubjet1PtChi40", "", s)->Fill(vsj[1].v.Perp()*1e-3, weight);
        if (sjn >= 3)
          h->h1D("lsubjet2PtChi40", "", s)->Fill(vsj[2].v.Perp()*1e-3, weight);
        if (sjn >= 4)
          h->h1D("lsubjet3PtChi40", "", s)->Fill(vsj[3].v.Perp()*1e-3, weight);
        if (sjn >= 5)
          h->h1D("lsubjet4PtChi40", "", s)->Fill(vsj[4].v.Perp()*1e-3, weight);

        if (sjn >= 1)
          h->h1D("lsubjet0EtaChi40", "", s)->Fill(vsj[0].v.Eta(), weight);
        if (sjn >= 2)
          h->h1D("lsubjet1EtaChi40", "", s)->Fill(vsj[1].v.Eta(), weight);
        if (sjn >= 3)
          h->h1D("lsubjet2EtaChi40", "", s)->Fill(vsj[2].v.Eta(), weight);
        if (sjn >= 4)
          h->h1D("lsubjet3EtaChi40", "", s)->Fill(vsj[3].v.Eta(), weight);
        if (sjn >= 5)
          h->h1D("lsubjet4EtaChi40", "", s)->Fill(vsj[4].v.Eta(), weight);

      }

      if (calchi > 2.8) {
        h->h1D("calsubjetNChi40", "", s)->Fill(casjn, weight);
        if (casjn >= 1)
          h->h1D("calsubjet0PtChi40", "", s)->Fill(cavsj[0].v.Perp()*1e-3, weight);
        if (casjn >= 2)
          h->h1D("calsubjet1PtChi40", "", s)->Fill(cavsj[1].v.Perp()*1e-3, weight);
        if (casjn >= 3)
          h->h1D("calsubjet2PtChi40", "", s)->Fill(cavsj[2].v.Perp()*1e-3, weight);
        if (casjn >= 4)
          h->h1D("calsubjet3PtChi40", "", s)->Fill(cavsj[3].v.Perp()*1e-3, weight);
        if (casjn >= 5)
          h->h1D("calsubjet4PtChi40", "", s)->Fill(cavsj[4].v.Perp()*1e-3, weight);
        if (casjn >= 1)
          h->h1D("calsubjet0EtaChi40", "", s)->Fill(cavsj[0].v.Eta(), weight);
        if (casjn >= 2)
          h->h1D("calsubjet1EtaChi40", "", s)->Fill(cavsj[1].v.Eta(), weight);
        if (casjn >= 3)
          h->h1D("calsubjet2EtaChi40", "", s)->Fill(cavsj[2].v.Eta(), weight);
        if (casjn >= 4)
          h->h1D("calsubjet3EtaChi40", "", s)->Fill(cavsj[3].v.Eta(), weight);
        if (casjn >= 5)
          h->h1D("calsubjet4EtaChi40", "", s)->Fill(cavsj[4].v.Eta(), weight);
      }

      if (lchib > 2.6) {
        h->h1D("lsubjetNChib40", "", s)->Fill(sjn, weight);
        if (sjn >= 1)
          h->h1D("lsubjet0PtChib40", "", s)->Fill(vsj[0].v.Perp()*1e-3, weight);
        if (sjn >= 2)
          h->h1D("lsubjet1PtChib40", "", s)->Fill(vsj[1].v.Perp()*1e-3, weight);
        if (sjn >= 3)
          h->h1D("lsubjet2PtChib40", "", s)->Fill(vsj[2].v.Perp()*1e-3, weight);
        if (sjn >= 4)
          h->h1D("lsubjet3PtChib40", "", s)->Fill(vsj[3].v.Perp()*1e-3, weight);
        if (sjn >= 5)
          h->h1D("lsubjet4PtChib40", "", s)->Fill(vsj[4].v.Perp()*1e-3, weight);

        if (sjn >= 1)
          h->h1D("lsubjet0EtaChib40", "", s)->Fill(vsj[0].v.Eta(), weight);
        if (sjn >= 2)
          h->h1D("lsubjet1EtaChib40", "", s)->Fill(vsj[1].v.Eta(), weight);
        if (sjn >= 3)
          h->h1D("lsubjet2EtaChib40", "", s)->Fill(vsj[2].v.Eta(), weight);
        if (sjn >= 4)
          h->h1D("lsubjet3EtaChib40", "", s)->Fill(vsj[3].v.Eta(), weight);
        if (sjn >= 5)
          h->h1D("lsubjet4EtaChib40", "", s)->Fill(vsj[4].v.Eta(), weight);

      }

      if (calchib > 2.8) {
        h->h1D("calsubjetNChib40", "", s)->Fill(casjn, weight);
        if (casjn >= 1)
          h->h1D("calsubjet0PtChib40", "", s)->Fill(cavsj[0].v.Perp()*1e-3, weight);
        if (casjn >= 2)
          h->h1D("calsubjet1PtChib40", "", s)->Fill(cavsj[1].v.Perp()*1e-3, weight);
        if (casjn >= 3)
          h->h1D("calsubjet2PtChib40", "", s)->Fill(cavsj[2].v.Perp()*1e-3, weight);
        if (casjn >= 4)
          h->h1D("calsubjet3PtChib40", "", s)->Fill(cavsj[3].v.Perp()*1e-3, weight);
        if (casjn >= 5)
          h->h1D("calsubjet4PtChib40", "", s)->Fill(cavsj[4].v.Perp()*1e-3, weight);
        if (casjn >= 1)
          h->h1D("calsubjet0EtaChib40", "", s)->Fill(cavsj[0].v.Eta(), weight);
        if (casjn >= 2)
          h->h1D("calsubjet1EtaChib40", "", s)->Fill(cavsj[1].v.Eta(), weight);
        if (casjn >= 3)
          h->h1D("calsubjet2EtaChib40", "", s)->Fill(cavsj[2].v.Eta(), weight);
        if (casjn >= 4)
          h->h1D("calsubjet3EtaChib40", "", s)->Fill(cavsj[3].v.Eta(), weight);
        if (casjn >= 5)
          h->h1D("calsubjet4EtaChib40", "", s)->Fill(cavsj[4].v.Eta(), weight);

      }
      if (m_hSvc.m_tree) {
        int bmatch = 0;
        int cabmatch = 0;
        for (int z = 0; z < e.jet().size(); ++z) {
          if (e.jet()[z].btag() && e.jet()[z].mom().DeltaR(lj.mom()) < 1.0) {
            bmatch = 1;
            break;
          }
        }
        for (int z = 0; z < e.jet().size(); ++z) {
          if (e.jet()[z].btag() && e.jet()[z].mom().DeltaR(calj.mom()) < 1.0) {
            cabmatch = 1;
            break;
          }
        }
        int syst = 0;
        if (s == "_subJesUp")   syst = 1;
        if (s == "_subJesDown") syst = 2;
        m_tree_syst = syst;
        m_tree_nsub = sd_sorted.size();
        m_tree_chi = lchi;
        m_tree_chib = lchib;
        m_tree_bmatch = bmatch;
        m_tree_pt = lj.mom().Perp()*1e-3;
        m_tree_eta = lj.mom().Eta();
        m_tree_truept = 0;//ljp.mom().Perp()*1e-3;
        m_tree_trueeta = 0;//ljp.mom().Eta();
        //m_tree_truept = lj.mom().Perp()*1e-3;
        //m_tree_trueeta = lj.mom().Eta();
        m_tree_m = lj.mom().M()*1e-3;
        m_tree_d12 = lj.split12()*1e-3;

        m_tree_cansub = casd_sorted.size();
        m_tree_cachi = calchi;
        m_tree_cachib = calchib;
        m_tree_cabmatch = cabmatch;
        m_tree_capt = calj.mom().Perp()*1e-3;
        m_tree_caeta = calj.mom().Eta();
        m_tree_catruept = 0;//caljp.mom().Perp()*1e-3;
        m_tree_catrueeta = 0;//caljp.mom().Eta();
        //m_tree_catruept = calj.mom().Perp()*1e-3;
        //m_tree_catrueeta = calj.mom().Eta();

        //m_tree_cahtt = caljp.htt() > 0;
        m_tree_cahtt = 0;

        m_tree_npv = e.npv();
        m_tree_mu = e.mu();
        m_tree_w = pweight;
        m_hSvc.m_tree->Fill();
      }
    }
  }

}

