#include "Plot.h"
#include "PlotSemilepTT.h"
#include "Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "HistogramService.h"

#include <fastjet/PseudoJet.hh>

#include <cmath>

//#include "SD-tosvn/Exception.h"

#include <algorithm>

using namespace fastjet;
//using namespace Deconstruction;

bool SortByPt(const LargeJet::Subjet &a, const LargeJet::Subjet &b) {
  return a.v.Perp() >= b.v.Perp();
}

PlotSemilepTT::PlotSemilepTT(const std::string &filename, int mode, const std::vector<std::string> &systs){

  m_requireCloseBtag = -1.0;

  minPt = 350; //350
  int nbinsPt = 16; //8

  double llb50[] = {350, 400, 450, 500, 600, 700};
  int nllb50 = 5;
  //double llb50[] = {200, 250, 300, 350, 400, 450, 500, 600, 700};
  //int nllb50 = 8;
  //

  m_massCut = -1;


  m_hSvc.create1D("lepPt", "; Lepton p_{T} [GeV] ; Events", 40, 25, 425);
  m_hSvc.create1D("lepEta", "; Lepton #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lepPhi", "; Lepton #phi [rad] ; Events", 64, -3.2, 3.2);
  m_hSvc.create1D("jetPt", "; Jet p_{T} [GeV] ; Events", 100, 25, 1025);
  m_hSvc.create1D("jetEta", "; Jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("jetPhi", "; Jet #phi [rad] ; Events", 64, -3.2, 3.2);

  m_hSvc.create1D("callargeJetPt", "; Large jet p_{T} [GeV] ; Events", 30, 0, 1500);
  m_hSvc.create1D("callargeJetEta", "; Large jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("callargeJetPhi", "; Large jet #phi [rad] ; Events", 64, -3.2, 3.2);
  m_hSvc.create1D("callargeJetM", "; Large jet M ; Events", 40, 0, 400);

  m_hSvc.create1D("llargeJetPt", "; Large jet p_{T} [GeV] ; Events", 30, 0, 1500);
  m_hSvc.create1D("llargeJetEta", "; Large jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("llargeJetPhi", "; Large jet #phi [rad] ; Events", 64, -3.2, 3.2);
  m_hSvc.create1D("llargeJetM", "; Large jet M ; Events", 40, 0, 400);

  m_hSvc.create1D("lsubjetN", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 11, -0.5, 10.5);
  m_hSvc.create1D("lsubjet0at3Pt", "; Anti-k_{t} R=1.0 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("lsubjet1at3Pt", "; Anti-k_{t} R=1.0 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet0Pt", "; Anti-k_{t} R=1.0 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("lsubjet1Pt", "; Anti-k_{t} R=1.0 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet2Pt", "; Anti-k_{t} R=1.0 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet3Pt", "; Anti-k_{t} R=1.0 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet4Pt", "; Anti-k_{t} R=1.0 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);

  m_hSvc.create1D("lsubjetBN", "; Anti-k_{t} R=1.0 b-tagged Sub-jet multiplicity ; Events", 3, -0.5, 2.5);
  m_hSvc.create1D("lsubjet0Bdr", "; Anti-k_{t} R=1.0 Leading sub-jet #Delta R to b-jet ; Events", 20, 0, 1.0);
  m_hSvc.create1D("lsubjet1Bdr", "; Anti-k_{t} R=1.0 Sub-leading sub-jet #Delta R to b-jet ; Events", 20, 0, 1.0);
  m_hSvc.create1D("lsubjet2Bdr", "; Anti-k_{t} R=1.0 Third leading sub-jet #Delta R to b-jet ; Events", 20, 0, 1.0);
  m_hSvc.create1D("lsubjet3Bdr", "; Anti-k_{t} R=1.0 Fourth leading sub-jet #Delta R to b-jet ; Events", 20, 0, 1.0);
  m_hSvc.create1D("lsubjet4Bdr", "; Anti-k_{t} R=1.0 Fifth leading sub-jet #Delta R to b-jet ; Events", 20, 0, 1.0);

  m_hSvc.create1D("lsubjet0Ldr", "; Anti-k_{t} R=1.0 Leading sub-jet #Delta R to b-jet ; Events", 20, 0, 1.0);
  m_hSvc.create1D("lsubjet1Ldr", "; Anti-k_{t} R=1.0 Sub-leading sub-jet #Delta R to b-jet ; Events", 20, 0, 1.0);
  m_hSvc.create1D("lsubjet2Ldr", "; Anti-k_{t} R=1.0 Third leading sub-jet #Delta R to b-jet ; Events", 20, 0, 1.0);
  m_hSvc.create1D("lsubjet3Ldr", "; Anti-k_{t} R=1.0 Fourth leading sub-jet #Delta R to b-jet ; Events", 20, 0, 1.0);
  m_hSvc.create1D("lsubjet4Ldr", "; Anti-k_{t} R=1.0 Fifth leading sub-jet #Delta R to b-jet ; Events", 20, 0, 1.0);

  m_hSvc.create1D("lsubjet0Eta", "; Anti-k_{t} R=1.0 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet1Eta", "; Anti-k_{t} R=1.0 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet2Eta", "; Anti-k_{t} R=1.0 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet3Eta", "; Anti-k_{t} R=1.0 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet4Eta", "; Anti-k_{t} R=1.0 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  m_hSvc.create1D("calsubjetN", "; C/A R=1.2 Sub-jet multiplicity ; Events", 11, -0.5, 10.5);
  m_hSvc.create1D("calsubjet0at3Pt", "; C/A R=1.2 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("calsubjet1at3Pt", "; C/A R=1.2 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet0Pt", "; C/A R=1.2 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("calsubjet1Pt", "; C/A R=1.2 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet2Pt", "; C/A R=1.2 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet3Pt", "; C/A R=1.2 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet4Pt", "; C/A R=1.2 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);

  m_hSvc.create1D("calsubjetBN", "; C/A R=1.2 b-tagged Sub-jet multiplicity ; Events", 3, -0.5, 2.5);
  m_hSvc.create1D("calsubjet0Bdr", "; C/A R=1.2 Leading sub-jet p_{T} [GeV] ; Events", 20, 0, 1.0);
  m_hSvc.create1D("calsubjet1Bdr", "; C/A R=1.2 Sub-leading sub-jet p_{T} [GeV] ; Events", 20, 0, 1.0);
  m_hSvc.create1D("calsubjet2Bdr", "; C/A R=1.2 Third leading sub-jet p_{T} [GeV] ; Events", 20, 0, 1.0);
  m_hSvc.create1D("calsubjet3Bdr", "; C/A R=1.2 Fourth leading sub-jet p_{T} [GeV] ; Events", 20, 0, 1.0);
  m_hSvc.create1D("calsubjet4Bdr", "; C/A R=1.2 Fifth leading sub-jet p_{T} [GeV] ; Events", 20, 0, 1.0);

  m_hSvc.create1D("calsubjet0Ldr", "; C/A R=1.2 Leading sub-jet p_{T} [GeV] ; Events", 20, 0, 1.0);
  m_hSvc.create1D("calsubjet1Ldr", "; C/A R=1.2 Sub-leading sub-jet p_{T} [GeV] ; Events", 20, 0, 1.0);
  m_hSvc.create1D("calsubjet2Ldr", "; C/A R=1.2 Third leading sub-jet p_{T} [GeV] ; Events", 20, 0, 1.0);
  m_hSvc.create1D("calsubjet3Ldr", "; C/A R=1.2 Fourth leading sub-jet p_{T} [GeV] ; Events", 20, 0, 1.0);
  m_hSvc.create1D("calsubjet4Ldr", "; C/A R=1.2 Fifth leading sub-jet p_{T} [GeV] ; Events", 20, 0, 1.0);

  m_hSvc.create1D("calsubjet0Eta", "; C/A R=1.2 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet1Eta", "; C/A R=1.2 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet2Eta", "; C/A R=1.2 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet3Eta", "; C/A R=1.2 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet4Eta", "; C/A R=1.2 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  double llbsj[] = {2.5, 3.5, 4.5, 5.5, 6.5, 10.5};
  m_hSvc.create1DVar("lsubjetNSF", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("calsubjetNSF", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);

  m_hSvc.create1DVar("lsubjetNLE", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("lsubjetNHE", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);

  m_hSvc.create1DVar("calsubjetNLE", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("calsubjetNHE", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);

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

  m_hSvc.create1D("lm12", "; m_{12} [GeV] ; Events", 20, 0, 200);
  m_hSvc.create1D("lm23", "; m_{23} [GeV] ; Events", 20, 0, 200);
  m_hSvc.create1D("lm123", "; m_{123} [GeV] ; Events", 20, 0, 300);

  m_hSvc.create1D("roc_largeJetPt", "; Large jet p_{T} [GeV] ; Events", 30, 0, 1500);
  m_hSvc.create1D("roc_largeJetEta", "; Large jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("roc_largeJetPhi", "; Large jet #phi [rad] ; Events", 64, -3.2, 3.2);
  m_hSvc.create1D("roc_largeJetM", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("roc_chi", "; SD log(#chi) ; Events", 30, -15, 15);
  m_hSvc.create1D("roc_sd_sig", "; log(SD signal weight) ; Events", 25, -40, -15);
  m_hSvc.create1D("roc_sd_bkg", "; log(SD background weight) ; Events", 25, -40, -15);

  m_hSvc.create1D("roc_subjetN", "; Sub-jet multiplicity ; Events", 11, -0.5, 10.5);
  m_hSvc.create1D("roc_subjet0Pt", "; Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("roc_subjet1Pt", "; Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("roc_subjet2Pt", "; Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("roc_subjet3Pt", "; Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("roc_subjet4Pt", "; Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);

  m_hSvc.create1D("roc_m12", "; m_{12} [GeV] ; Events", 20, 0, 200);
  m_hSvc.create1D("roc_m23", "; m_{23} [GeV] ; Events", 20, 0, 200);
  m_hSvc.create1D("roc_m123", "; m_{123} [GeV] ; Events", 20, 0, 300);

  m_hSvc.create1D("mu", "; <#mu> ; Events", 5, 12, 22);
  m_hSvc.create1D("npv", "; npv ; Events", 5, 10, 20);

  m_hSvc.create1D("met", "; Missing E_{T} [GeV] ; Events", 20, 0, 400);
  m_hSvc.create1D("met_phi", "; Missing energy #phi  ; Events", 64, -3.2, 3.2);
  m_hSvc.create1D("mtw", "; m_{T} [GeV] ; Events", 30, 0, 300);

  m_hSvc.create1D("ldrfjb", "; #Delta R(large-R jet, b-jet) ; Events", 45, 0, 4.5);
  m_hSvc.create1D("ldpfjl", "; #Delta #phi(large-R jet, lepton) ; Events", 32, 0, 3.2);
  m_hSvc.create1D("ldrfjl", "; #Delta R(large-R jet, lepton) ; Events", 45, 0, 4.5);
  m_hSvc.create1D("ldrjl", "; min #Delta R(jet, lepton) ; Events", 45, 0, 4.5);
  m_hSvc.create1D("ldpml", "; #Delta #phi(MET, lepton) ; Events", 64, -3.2, 3.2);
  m_hSvc.create1D("ldrml", "; #Delta R(MET, lepton) ; Events", 45, 0, 4.5);

  m_hSvc.create1D("subjetlcpt", "; subjet lc pt [GeV] ; Events", 50, 0, 50);

  m_hSvc.create1D("llargeJetPtPos", "; Large jet p_{T} [GeV] ; Events", nbinsPt, minPt, 1150);
  m_hSvc.create1D("llargeJetPtPosCa", "; Large jet p_{T} [GeV] ; Events", nbinsPt, minPt, 1150);

  m_hSvc.create1D("llargeJetMPos", "; Large jet M [GeV] ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosCa", "; Large jet M [GeV] ; Events", 40, 0, 400);

  // dist. after chi cuts
  //WP1
  m_hSvc.create1DVar("lsubjetNChiWP1", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("lsubjetNLEChiWP1", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("lsubjetNHEChiWP1", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1D("lsubjetNLETestChiWP1", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 11, -0.5, 10.5);

  m_hSvc.create1D("lsubjet0PtChiWP1", "; Anti-k_{t} R=1.0 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("lsubjet1PtChiWP1", "; Anti-k_{t} R=1.0 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet2PtChiWP1", "; Anti-k_{t} R=1.0 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet3PtChiWP1", "; Anti-k_{t} R=1.0 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet4PtChiWP1", "; Anti-k_{t} R=1.0 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);

  m_hSvc.create1DVar("calsubjetNChiWP1", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("calsubjetNLEChiWP1", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("calsubjetNHEChiWP1", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);

  m_hSvc.create1D("calsubjet0PtChiWP1", "; C/A R=1.2 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("calsubjet1PtChiWP1", "; C/A R=1.2 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet2PtChiWP1", "; C/A R=1.2 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet3PtChiWP1", "; C/A R=1.2 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet4PtChiWP1", "; C/A R=1.2 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);


  m_hSvc.create1DVar("lsubjetNChibWP1", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("lsubjetNLEChibWP1", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("lsubjetNHEChibWP1", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);

  m_hSvc.create1D("lsubjet0PtChibWP1", "; Anti-k_{t} R=1.0 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("lsubjet1PtChibWP1", "; Anti-k_{t} R=1.0 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet2PtChibWP1", "; Anti-k_{t} R=1.0 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet3PtChibWP1", "; Anti-k_{t} R=1.0 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet4PtChibWP1", "; Anti-k_{t} R=1.0 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);

  m_hSvc.create1DVar("calsubjetNChibWP1", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("calsubjetNLEChibWP1", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("calsubjetNHEChibWP1", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);

  m_hSvc.create1D("calsubjet0PtChibWP1", "; C/A R=1.2 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("calsubjet1PtChibWP1", "; C/A R=1.2 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet2PtChibWP1", "; C/A R=1.2 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet3PtChibWP1", "; C/A R=1.2 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet4PtChibWP1", "; C/A R=1.2 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);

  m_hSvc.create1D("lsubjet0EtaChiWP1", "; Anti-k_{t} R=1.0 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet1EtaChiWP1", "; Anti-k_{t} R=1.0 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet2EtaChiWP1", "; Anti-k_{t} R=1.0 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet3EtaChiWP1", "; Anti-k_{t} R=1.0 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet4EtaChiWP1", "; Anti-k_{t} R=1.0 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  m_hSvc.create1D("calsubjet0EtaChiWP1", "; C/A R=1.2 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet1EtaChiWP1", "; C/A R=1.2 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet2EtaChiWP1", "; C/A R=1.2 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet3EtaChiWP1", "; C/A R=1.2 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet4EtaChiWP1", "; C/A R=1.2 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  m_hSvc.create1D("lsubjet0EtaChibWP1", "; Anti-k_{t} R=1.0 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet1EtaChibWP1", "; Anti-k_{t} R=1.0 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet2EtaChibWP1", "; Anti-k_{t} R=1.0 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet3EtaChibWP1", "; Anti-k_{t} R=1.0 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet4EtaChibWP1", "; Anti-k_{t} R=1.0 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  m_hSvc.create1D("calsubjet0EtaChibWP1", "; C/A R=1.2 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet1EtaChibWP1", "; C/A R=1.2 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet2EtaChibWP1", "; C/A R=1.2 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet3EtaChibWP1", "; C/A R=1.2 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet4EtaChibWP1", "; C/A R=1.2 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  m_hSvc.create1D("llargeJetPtPosChiWP1", "; Large jet p_{T} [GeV] ; Events", nbinsPt, minPt, 1150);
  m_hSvc.create1D("llargeJetPtPosChibWP1", "; Large jet p_{T} [GeV] ; Events", nbinsPt, minPt, 1150);
  m_hSvc.create1D("llargeJetPtPosCaChiWP1", "; Large jet p_{T} [GeV] ; Events", nbinsPt, minPt, 1150);
  m_hSvc.create1D("llargeJetPtPosCaChibWP1", "; Large jet p_{T} [GeV] ; Events", nbinsPt, minPt, 1150);

  m_hSvc.create1D("llargeJetMPosChiWP1", "; Large jet M [GeV] ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosCaChiWP1", "; Large jet M [GeV] ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosChibWP1", "; Large jet M [GeV] ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosCaChibWP1", "; Large jet M [GeV] ; Events", 40, 0, 400);

  m_hSvc.create1D("llargeJetSDMChiWP1", "; Large jet M [GeV] ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMCaChiWP1", "; Large jet M [GeV] ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMChibWP1", "; Large jet M [GeV] ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMCaChibWP1", "; Large jet M [GeV] ; Events", 40, 0, 400);

  m_hSvc.create1D("lm12ChiWP1", "; m_{12} [GeV] ; Events", 20, 0, 200);
  m_hSvc.create1D("lm23ChiWP1", "; m_{23} [GeV] ; Events", 20, 0, 200);
  m_hSvc.create1D("lm123ChiWP1", "; m_{123} [GeV] ; Events", 20, 0, 300);

  m_hSvc.create1D("lm12ChibWP1", "; m_{12} [GeV] ; Events", 20, 0, 200);
  m_hSvc.create1D("lm23ChibWP1", "; m_{23} [GeV] ; Events", 20, 0, 200);
  m_hSvc.create1D("lm123ChibWP1", "; m_{123} [GeV] ; Events", 20, 0, 300);

  //WP2
  m_hSvc.create1DVar("lsubjetNChiWP2", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("lsubjetNLEChiWP2", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("lsubjetNHEChiWP2", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1D("lsubjetNLETestChiWP2", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 11, -0.5, 10.5);

  m_hSvc.create1D("lsubjet0PtChiWP2", "; Anti-k_{t} R=1.0 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("lsubjet1PtChiWP2", "; Anti-k_{t} R=1.0 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet2PtChiWP2", "; Anti-k_{t} R=1.0 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet3PtChiWP2", "; Anti-k_{t} R=1.0 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet4PtChiWP2", "; Anti-k_{t} R=1.0 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);

  m_hSvc.create1DVar("calsubjetNChiWP2", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("calsubjetNLEChiWP2", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("calsubjetNHEChiWP2", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);

  m_hSvc.create1D("calsubjet0PtChiWP2", "; C/A R=1.2 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("calsubjet1PtChiWP2", "; C/A R=1.2 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet2PtChiWP2", "; C/A R=1.2 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet3PtChiWP2", "; C/A R=1.2 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet4PtChiWP2", "; C/A R=1.2 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);


  m_hSvc.create1DVar("lsubjetNChibWP2", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("lsubjetNLEChibWP2", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("lsubjetNHEChibWP2", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);

  m_hSvc.create1D("lsubjet0PtChibWP2", "; Anti-k_{t} R=1.0 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("lsubjet1PtChibWP2", "; Anti-k_{t} R=1.0 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet2PtChibWP2", "; Anti-k_{t} R=1.0 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet3PtChibWP2", "; Anti-k_{t} R=1.0 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet4PtChibWP2", "; Anti-k_{t} R=1.0 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);

  m_hSvc.create1DVar("calsubjetNChibWP2", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("calsubjetNLEChibWP2", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("calsubjetNHEChibWP2", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);

  m_hSvc.create1D("calsubjet0PtChibWP2", "; C/A R=1.2 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("calsubjet1PtChibWP2", "; C/A R=1.2 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet2PtChibWP2", "; C/A R=1.2 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet3PtChibWP2", "; C/A R=1.2 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet4PtChibWP2", "; C/A R=1.2 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);

  m_hSvc.create1D("lsubjet0EtaChiWP2", "; Anti-k_{t} R=1.0 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet1EtaChiWP2", "; Anti-k_{t} R=1.0 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet2EtaChiWP2", "; Anti-k_{t} R=1.0 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet3EtaChiWP2", "; Anti-k_{t} R=1.0 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet4EtaChiWP2", "; Anti-k_{t} R=1.0 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  m_hSvc.create1D("calsubjet0EtaChiWP2", "; C/A R=1.2 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet1EtaChiWP2", "; C/A R=1.2 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet2EtaChiWP2", "; C/A R=1.2 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet3EtaChiWP2", "; C/A R=1.2 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet4EtaChiWP2", "; C/A R=1.2 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  m_hSvc.create1D("lsubjet0EtaChibWP2", "; Anti-k_{t} R=1.0 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet1EtaChibWP2", "; Anti-k_{t} R=1.0 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet2EtaChibWP2", "; Anti-k_{t} R=1.0 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet3EtaChibWP2", "; Anti-k_{t} R=1.0 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet4EtaChibWP2", "; Anti-k_{t} R=1.0 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  m_hSvc.create1D("calsubjet0EtaChibWP2", "; C/A R=1.2 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet1EtaChibWP2", "; C/A R=1.2 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet2EtaChibWP2", "; C/A R=1.2 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet3EtaChibWP2", "; C/A R=1.2 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet4EtaChibWP2", "; C/A R=1.2 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  m_hSvc.create1D("llargeJetPtPosChiWP2", "; Large jet p_{T} [GeV] ; Events", nbinsPt, minPt, 1150);
  m_hSvc.create1D("llargeJetPtPosChibWP2", "; Large jet p_{T} [GeV] ; Events", nbinsPt, minPt, 1150);
  m_hSvc.create1D("llargeJetPtPosCaChiWP2", "; Large jet p_{T} [GeV] ; Events", nbinsPt, minPt, 1150);
  m_hSvc.create1D("llargeJetPtPosCaChibWP2", "; Large jet p_{T} [GeV] ; Events", nbinsPt, minPt, 1150);

  m_hSvc.create1D("llargeJetMPosChiWP2", "; Large jet M [GeV] ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosCaChiWP2", "; Large jet M [GeV] ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosChibWP2", "; Large jet M [GeV] ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosCaChibWP2", "; Large jet M [GeV] ; Events", 40, 0, 400);

  m_hSvc.create1D("llargeJetSDMChiWP2", "; Large jet M [GeV] ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMCaChiWP2", "; Large jet M [GeV] ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMChibWP2", "; Large jet M [GeV] ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMCaChibWP2", "; Large jet M [GeV] ; Events", 40, 0, 400);

  m_hSvc.create1D("lm12ChiWP2", "; m_{12} [GeV] ; Events", 20, 0, 200);
  m_hSvc.create1D("lm23ChiWP2", "; m_{23} [GeV] ; Events", 20, 0, 200);
  m_hSvc.create1D("lm123ChiWP2", "; m_{123} [GeV] ; Events", 20, 0, 300);

  m_hSvc.create1D("lm12ChibWP2", "; m_{12} [GeV] ; Events", 20, 0, 200);
  m_hSvc.create1D("lm23ChibWP2", "; m_{23} [GeV] ; Events", 20, 0, 200);
  m_hSvc.create1D("lm123ChibWP2", "; m_{123} [GeV] ; Events", 20, 0, 300);

  //WP3
  m_hSvc.create1DVar("lsubjetNChiWP3", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("lsubjetNLEChiWP3", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("lsubjetNHEChiWP3", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1D("lsubjetNLETestChiWP3", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 11, -0.5, 10.5);

  m_hSvc.create1D("lsubjet0PtChiWP3", "; Anti-k_{t} R=1.0 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("lsubjet1PtChiWP3", "; Anti-k_{t} R=1.0 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet2PtChiWP3", "; Anti-k_{t} R=1.0 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet3PtChiWP3", "; Anti-k_{t} R=1.0 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet4PtChiWP3", "; Anti-k_{t} R=1.0 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);

  m_hSvc.create1DVar("calsubjetNChiWP3", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("calsubjetNLEChiWP3", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("calsubjetNHEChiWP3", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);

  m_hSvc.create1D("calsubjet0PtChiWP3", "; C/A R=1.2 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("calsubjet1PtChiWP3", "; C/A R=1.2 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet2PtChiWP3", "; C/A R=1.2 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet3PtChiWP3", "; C/A R=1.2 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet4PtChiWP3", "; C/A R=1.2 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);


  m_hSvc.create1DVar("lsubjetNChibWP3", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("lsubjetNLEChibWP3", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("lsubjetNHEChibWP3", "; Anti-k_{t} R=1.0 Sub-jet multiplicity ; Events", 5, llbsj);

  m_hSvc.create1D("lsubjet0PtChibWP3", "; Anti-k_{t} R=1.0 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("lsubjet1PtChibWP3", "; Anti-k_{t} R=1.0 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet2PtChibWP3", "; Anti-k_{t} R=1.0 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet3PtChibWP3", "; Anti-k_{t} R=1.0 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("lsubjet4PtChibWP3", "; Anti-k_{t} R=1.0 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);

  m_hSvc.create1DVar("calsubjetNChibWP3", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("calsubjetNLEChibWP3", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);
  m_hSvc.create1DVar("calsubjetNHEChibWP3", "; C/A R=1.2 Sub-jet multiplicity ; Events", 5, llbsj);

  m_hSvc.create1D("calsubjet0PtChibWP3", "; C/A R=1.2 Leading sub-jet p_{T} [GeV] ; Events", 50, 0, 1000);
  m_hSvc.create1D("calsubjet1PtChibWP3", "; C/A R=1.2 Sub-leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet2PtChibWP3", "; C/A R=1.2 Third leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet3PtChibWP3", "; C/A R=1.2 Fourth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);
  m_hSvc.create1D("calsubjet4PtChibWP3", "; C/A R=1.2 Fifth leading sub-jet p_{T} [GeV] ; Events", 25, 0, 500);

  m_hSvc.create1D("lsubjet0EtaChiWP3", "; Anti-k_{t} R=1.0 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet1EtaChiWP3", "; Anti-k_{t} R=1.0 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet2EtaChiWP3", "; Anti-k_{t} R=1.0 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet3EtaChiWP3", "; Anti-k_{t} R=1.0 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet4EtaChiWP3", "; Anti-k_{t} R=1.0 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  m_hSvc.create1D("calsubjet0EtaChiWP3", "; C/A R=1.2 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet1EtaChiWP3", "; C/A R=1.2 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet2EtaChiWP3", "; C/A R=1.2 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet3EtaChiWP3", "; C/A R=1.2 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet4EtaChiWP3", "; C/A R=1.2 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  m_hSvc.create1D("lsubjet0EtaChibWP3", "; Anti-k_{t} R=1.0 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet1EtaChibWP3", "; Anti-k_{t} R=1.0 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet2EtaChibWP3", "; Anti-k_{t} R=1.0 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet3EtaChibWP3", "; Anti-k_{t} R=1.0 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lsubjet4EtaChibWP3", "; Anti-k_{t} R=1.0 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  m_hSvc.create1D("calsubjet0EtaChibWP3", "; C/A R=1.2 Leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet1EtaChibWP3", "; C/A R=1.2 Sub-leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet2EtaChibWP3", "; C/A R=1.2 Third leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet3EtaChibWP3", "; C/A R=1.2 Fourth leading sub-jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("calsubjet4EtaChibWP3", "; C/A R=1.2 Fifth leading sub-jet #eta ; Events", 50, -2.5, 2.5);

  m_hSvc.create1D("llargeJetPtPosChiWP3", "; Large jet p_{T} [GeV] ; Events", nbinsPt, minPt, 1150);
  m_hSvc.create1D("llargeJetPtPosChibWP3", "; Large jet p_{T} [GeV] ; Events", nbinsPt, minPt, 1150);
  m_hSvc.create1D("llargeJetPtPosCaChiWP3", "; Large jet p_{T} [GeV] ; Events", nbinsPt, minPt, 1150);
  m_hSvc.create1D("llargeJetPtPosCaChibWP3", "; Large jet p_{T} [GeV] ; Events", nbinsPt, minPt, 1150);

  m_hSvc.create1D("llargeJetMPosChiWP3", "; Large jet M [GeV] ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosCaChiWP3", "; Large jet M [GeV] ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosChibWP3", "; Large jet M [GeV] ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosCaChibWP3", "; Large jet M [GeV] ; Events", 40, 0, 400);

  m_hSvc.create1D("llargeJetSDMChiWP3", "; Large jet M [GeV] ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMCaChiWP3", "; Large jet M [GeV] ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMChibWP3", "; Large jet M [GeV] ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMCaChibWP3", "; Large jet M [GeV] ; Events", 40, 0, 400);

  m_hSvc.create1D("lm12ChiWP3", "; m_{12} [GeV] ; Events", 20, 0, 200);
  m_hSvc.create1D("lm23ChiWP3", "; m_{23} [GeV] ; Events", 20, 0, 200);
  m_hSvc.create1D("lm123ChiWP3", "; m_{123} [GeV] ; Events", 20, 0, 300);

  m_hSvc.create1D("lm12ChibWP3", "; m_{12} [GeV] ; Events", 20, 0, 200);
  m_hSvc.create1D("lm23ChibWP3", "; m_{23} [GeV] ; Events", 20, 0, 200);
  m_hSvc.create1D("lm123ChibWP3", "; m_{123} [GeV] ; Events", 20, 0, 300);

  // efficiencies
  m_hSvc.create1D("muEff", "; <#mu> ; Events", 5, 12, 22);
  m_hSvc.create1D("muEffCa", "; <#mu> ; Events", 5, 12, 22);

  m_hSvc.create1D("muEffChiWP1", "; <#mu> ; Events", 5, 12, 22);
  m_hSvc.create1D("muEffCaChiWP1", "; <#mu> ; Events", 5, 12, 22);
  m_hSvc.create1D("muEffChibWP1", "; <#mu> ; Events", 5, 12, 22);
  m_hSvc.create1D("muEffCaChibWP1", "; <#mu> ; Events", 5, 12, 22);

  m_hSvc.create1D("muEffChiWP2", "; <#mu> ; Events", 5, 12, 22);
  m_hSvc.create1D("muEffCaChiWP2", "; <#mu> ; Events", 5, 12, 22);
  m_hSvc.create1D("muEffChibWP2", "; <#mu> ; Events", 5, 12, 22);
  m_hSvc.create1D("muEffCaChibWP2", "; <#mu> ; Events", 5, 12, 22);

  m_hSvc.create1D("muEffChiWP3", "; <#mu> ; Events", 5, 12, 22);
  m_hSvc.create1D("muEffCaChiWP3", "; <#mu> ; Events", 5, 12, 22);
  m_hSvc.create1D("muEffChibWP3", "; <#mu> ; Events", 5, 12, 22);
  m_hSvc.create1D("muEffCaChibWP3", "; <#mu> ; Events", 5, 12, 22);

  m_hSvc.create1D("npvEff", "; NPV ; Events", 5, 4, 24);
  m_hSvc.create1D("npvEffCa", "; NPV ; Events", 5, 4, 24);

  m_hSvc.create1D("npvEffChiWP1", "; NPV ; Events", 5, 4, 24);
  m_hSvc.create1D("npvEffCaChiWP1", "; NPV ; Events", 5, 4, 24);
  m_hSvc.create1D("npvEffChibWP1", "; NPV ; Events", 5, 4, 24);
  m_hSvc.create1D("npvEffCaChibWP1", "; NPV ; Events", 5, 4, 24);

  m_hSvc.create1D("npvEffChiWP2", "; NPV ; Events", 5, 4, 24);
  m_hSvc.create1D("npvEffCaChiWP2", "; NPV ; Events", 5, 4, 24);
  m_hSvc.create1D("npvEffChibWP2", "; NPV ; Events", 5, 4, 24);
  m_hSvc.create1D("npvEffCaChibWP2", "; NPV ; Events", 5, 4, 24);

  m_hSvc.create1D("npvEffChiWP3", "; NPV ; Events", 5, 4, 24);
  m_hSvc.create1D("npvEffCaChiWP3", "; NPV ; Events", 5, 4, 24);
  m_hSvc.create1D("npvEffChibWP3", "; NPV ; Events", 5, 4, 24);
  m_hSvc.create1D("npvEffCaChibWP3", "; NPV ; Events", 5, 4, 24);

  m_hSvc.create1D("pre", ";  ; Events", 1, 0, 1);
  m_hSvc.create1D("preCa", ";  ; Events", 1, 0, 1);

  m_hSvc.create1D("preChiWP1", ";  ; Events", 1, 0, 1);
  m_hSvc.create1D("preChibWP1", ";  ; Events", 1, 0, 1);
  m_hSvc.create1D("preCaChiWP1", ";  ; Events", 1, 0, 1);
  m_hSvc.create1D("preCaChibWP1", ";  ; Events", 1, 0, 1);

  m_hSvc.create1D("preChiWP2", ";  ; Events", 1, 0, 1);
  m_hSvc.create1D("preChibWP2", ";  ; Events", 1, 0, 1);
  m_hSvc.create1D("preCaChiWP2", ";  ; Events", 1, 0, 1);
  m_hSvc.create1D("preCaChibWP2", ";  ; Events", 1, 0, 1);

  m_hSvc.create1D("preChiWP3", ";  ; Events", 1, 0, 1);
  m_hSvc.create1D("preChibWP3", ";  ; Events", 1, 0, 1);
  m_hSvc.create1D("preCaChiWP3", ";  ; Events", 1, 0, 1);
  m_hSvc.create1D("preCaChibWP3", ";  ; Events", 1, 0, 1);

  m_hSvc.create1DVar("llargeJetPtEff", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtEffCa", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);

  m_hSvc.create1DVar("llargeJetPtEffChiWP1", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtEffChibWP1", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtEffCaChiWP1", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtEffCaChibWP1", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);

  m_hSvc.create1DVar("llargeJetPtEffChiWP2", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtEffChibWP2", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtEffCaChiWP2", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtEffCaChibWP2", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);

  m_hSvc.create1DVar("llargeJetPtEffChiWP3", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtEffChibWP3", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtEffCaChiWP3", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtEffCaChibWP3", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);

  m_hSvc.create1DVar("llargeJetPtLEEff", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtLEEffCa", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);

  m_hSvc.create1DVar("llargeJetPtLEEffChiWP1", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtLEEffChibWP1", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtLEEffCaChiWP1", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtLEEffCaChibWP1", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);

  m_hSvc.create1DVar("llargeJetPtLEEffChiWP2", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtLEEffChibWP2", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtLEEffCaChiWP2", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtLEEffCaChibWP2", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);

  m_hSvc.create1DVar("llargeJetPtLEEffChiWP3", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtLEEffChibWP3", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtLEEffCaChiWP3", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtLEEffCaChibWP3", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);

  m_hSvc.create1DVar("llargeJetPtHEEff", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtHEEffCa", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);

  m_hSvc.create1DVar("llargeJetPtHEEffChiWP1", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtHEEffChibWP1", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtHEEffCaChiWP1", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtHEEffCaChibWP1", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);

  m_hSvc.create1DVar("llargeJetPtHEEffChiWP2", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtHEEffChibWP2", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtHEEffCaChiWP2", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtHEEffCaChibWP2", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);

  m_hSvc.create1DVar("llargeJetPtHEEffChiWP3", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtHEEffChibWP3", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtHEEffCaChiWP3", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);
  m_hSvc.create1DVar("llargeJetPtHEEffCaChibWP3", "; Large jet p_{T} [GeV] ; Events", nllb50, llb50);

i/*  signal = new TopGluonModel(param);
  background = new BackgroundModel(param);
  isr = new ISRModel(param);
  deconstruct = new Deconstruct(param, *signal, *background, *isr);

  paramb["useBtag"] = 1;

  signalb = new TopGluonModel(paramb);
  backgroundb = new BackgroundModel(paramb);
  isrb = new ISRModel(paramb);
  deconstructb = new Deconstruct(paramb, *signalb, *backgroundb, *isrb);

  m_HFsysTool = new HFsys("packages/ttResoSingleLepton/data/HFsys_WjetYields.txt", 
                          "packages/ttResoSingleLepton/data/HFsysboosted_SFs.root", HFSYS::BOOSTED);

  if (m_mode >= 2) { // only for ROC
    m_hSvc.m_tree->Branch("syst", &m_tree_syst);

    m_hSvc.m_tree->Branch("nsub", &m_tree_nsub);
    m_hSvc.m_tree->Branch("chi", &m_tree_chi);
    m_hSvc.m_tree->Branch("chib", &m_tree_chib);
    m_hSvc.m_tree->Branch("bmatch", &m_tree_bmatch);
    m_hSvc.m_tree->Branch("pt", &m_tree_pt);
    m_hSvc.m_tree->Branch("subpt0", &m_tree_subpt0);
    m_hSvc.m_tree->Branch("subpt1", &m_tree_subpt1);
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
    m_hSvc.m_tree->Branch("casubpt0", &m_tree_casubpt0);
    m_hSvc.m_tree->Branch("casubpt1", &m_tree_casubpt1);
    m_hSvc.m_tree->Branch("caeta", &m_tree_caeta);
    m_hSvc.m_tree->Branch("catruept", &m_tree_catruept);
    m_hSvc.m_tree->Branch("catrueeta", &m_tree_catrueeta);

    m_hSvc.m_tree->Branch("cahtt", &m_tree_cahtt);

    m_hSvc.m_tree->Branch("npv", &m_tree_npv);
    m_hSvc.m_tree->Branch("mu", &m_tree_mu);
    m_hSvc.m_tree->Branch("w", &m_tree_w);
  }*/

}

PlotSemilepTT::~PlotSemilepTT() {
 // if (deconstruct) delete deconstruct;
  //if (isr) delete isr;
  //if (background) delete background;
  if (signal) delete signal;
  //if (m_HFsysTool) delete m_HFsysTool;
}

std::vector<LargeJet::Subjet> PlotSemilepTT::applySjUnc(const std::vector<LargeJet::Subjet> &vsj, const std::vector<LargeJet::Subjet> &vsj_unc, const std::string &s) {
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

void PlotSemilepTT::run(const Event &e, double weight, double pweight, const std::string &s) {
  // HFOR for Wjets
  const int datasetsWjets[] = {190001, 190002, 190003, 190004, 190005, 190011, 190012, 190013, 190014, 190015, 190021, 190022, 190023, 190024, 190025, 190050, 190051, 190052, 190053, 190040, 190041, 190042, 190043, 190030, 190031, 190032, 190033, 190034, 110801, 110802, 110803, 110804, 126601, 126602, 126603, 126604, 126605, 126606, 126607, 126608, 126609, 117700, 117701, 117702, 117703, 117704, 117705, 117680, 117681, 117682, 117683, 117684, 117685, 117690, 117691, 117692, 117693, 117694, 117695, 147025, 147026, 147027, 147028, 147029, 147030, 147033, 147034, 147035, 147036, 147037, 147038, 147041, 147042, 147043, 147044, 147045, 147046, 200156, 200157, 200158, 200159, 200056, 200057, 200058, 200059, 200060, 200256, 200257, 200258, 200259};
  for (int k = 0; k < sizeof(datasetsWjets)/sizeof(int); ++k) {
    if (e.channelNumber() == datasetsWjets[k]) {
      if (m_mode < 2) {
        std::string ws = "nominal";
        if (s == "_wjetsnlo") ws = "nlo";
        if (s == "_wjetsqcdUp") ws = "qcdup";
        if (s == "_wjetsqcdDown") ws = "qcddw";
        std::string lepton_channel = "el";
        if (e.electron().size() == 1) lepton_channel = "el";
        if (e.muon().size() == 1) lepton_channel = "mu";
        //double w_wjets = m_HFsysTool->GetWjetSFWeight(ws, lepton_channel, e.hfor(), e.jet().size());
        weight *= w_wjets;
        pweight *= w_wjets;
      }
      break;
    }
  }

  for (int k = 0; k < sizeof(datasetsWjets)/sizeof(int); ++k) {
    if (e.channelNumber() == datasetsWjets[k]) {
      if (e.hfor() < 0 || e.hfor() >= 4)
        return;
    }
  }

  // HFOR for Zjets
  const int datasetsZjets[] = {117650, 117651, 117652, 117653, 117654, 117655, 117660, 117661, 117662, 117663, 117664, 117665, 117670, 117671, 117672, 117673, 117674, 117675, 110817, 110818, 110819, 110820, 110821, 110822, 110823, 110824, 110825, 110826, 110827, 110828, 110805, 110806, 110807, 110808, 110809, 110810, 110811, 110812, 110813, 110814, 110815, 110816, 147105, 147106, 147107, 147108, 147109, 147110, 200332, 200333, 200334, 200335, 200432, 200433, 200434, 200435, 147113, 147114, 147115, 147116, 147117, 147118, 147121, 147122, 147123, 147124, 147125, 147126, 200340, 200341, 200342, 200343, 200348, 200349, 200350, 200351, 200440, 200441, 200442, 200443, 200448, 200449, 200450, 200451};
  for (int k = 0; k < sizeof(datasetsZjets)/sizeof(int); ++k) {
    if (e.channelNumber() == datasetsZjets[k]) {
      if (e.hfor() < 0 || e.hfor() >= 4)
        return;
    }
  }
  const int datasetsWjetsLowPt[] = {110801, 110802, 110803, 110804, 126601, 126602, 126603, 126604, 126605, 126606, 126607, 126608, 126609, 117700, 117701, 117702, 117703, 117704, 117705, 117680, 117681, 117682, 117683, 117684, 117685, 117690, 117691, 117692, 117693, 117694, 117695, 147025, 147026, 147027, 147028, 147029, 147030, 147033, 147034, 147035, 147036, 147037, 147038, 147041, 147042, 147043, 147044, 147045, 147046, 200156, 200157, 200158, 200159, 200056, 200057, 200058, 200059, 200060, 200256, 200257, 200258, 200259};
  for (int k = 0; k < sizeof(datasetsWjetsLowPt)/sizeof(int); ++k) {
    if (e.channelNumber() == datasetsWjetsLowPt[k]) {
      // only available for nominal!
      if (e.partLargeJet().size() >= 1 && e.partLargeJet()[0].mom().Perp() >= 250e3) {
        return;
      }
      break;
    }
  }
  const int datasetsWjetsHighPt[] = {190001, 190002, 190003, 190004, 190005, 190011, 190012, 190013, 190014, 190015, 190021, 190022, 190023, 190024, 190025, 190050, 190051, 190052, 190053, 190040, 190041, 190042, 190043, 190030, 190031, 190032, 190033, 190034};
  for (int k = 0; k < sizeof(datasetsWjetsHighPt)/sizeof(int); ++k) {
    if (e.channelNumber() == datasetsWjetsHighPt[k]) {
      // only available for nominal!
      if (e.partLargeJet().size() >= 1 && e.partLargeJet()[0].mom().Perp() < 250e3) {
        return;
      }
      break;
    }
  }

  // top pt reweighting (only when splitting into matched and non-matched)
  // this allows onlyMatche == 0 for the systematic unc.
  if ((e.channelNumber() == 117050 || e.channelNumber() == 117001) && (m_ttbarPtWeight != 0) && e.partMom().size() != 0) {
    double tw = topPtWeight(e);
    weight *= tw;
    pweight *= tw;
  }
  HistogramService *h = &m_hSvc;
  
  if (e.passReco()) {

    TLorentzVector l;
    if (m_mode == 0) {
      l = e.electron()[0].mom();
    } else if (m_mode == 1) {
      l = e.muon()[0].mom();
    }

    if (m_mode <= 1) {
      // only run SD for large jets passing the following cuts for basic studies
      // pT > 350 GeV makes top be contained in the large jet
      int fji = -1;
      for (int k = 0; k < e.largeJet().size(); ++k) {
        if (e.largeJet()[k].passLoose() && e.largeJet()[k].mom().Perp() > minPt*1e3 && l.DeltaR(e.largeJet()[k].mom()) > 1.5) {
          fji = k;
          break;
        }
      }
    
     
  // for Signal ROC match the particle jet to the top and the reco jet to the particle jet
  // and demand pT > 550 GeV, |eta| < 1.2
  // for background only use the leading jet without extra cuts
  if (m_mode == 2 || m_mode == 3) {
    bool isSig = true;
    if (m_mode == 2) isSig = true;
    else if (m_mode == 3) isSig = false;

    int rocj = -1;
    int rocjp = -1;
    int rocjCA = -1;

    
    if (!isSig && e.passPart() && e.passReco()) {
      rocjp = 0;
      // loop over all reco jets
      for (int m = 0; m < e.largeJet().size(); ++m) {
        // apply pt and eta cuts
        if (e.largeJet()[m].mom().Perp() > 350e3) {
          // calculate dR
          if (e.partLargeJet()[0].mom().Perp() > 150e3 && std::fabs(e.partLargeJet()[0].mom().Eta()) < 2.0) {
            double dr = e.partLargeJet()[0].mom().DeltaR(e.largeJet()[m].mom());
            // if matched ...
            if (dr < 0.75) {
              rocj = m;
              ok = true;
              break;
            }
          }
        }
      }

    }
   
  } //mode 2 ou 3
    }//mode <=1
}// passReco
}
