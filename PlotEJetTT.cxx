#include "Plot.h"
#include "PlotEJetTT.h"
#include "Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "HistogramService.h"

#include <fastjet/PseudoJet.hh>

#include <cmath>

#include "SD-tosvn/Exception.h"

#include <algorithm>

using namespace fastjet;
using namespace Deconstruction;

bool SortByPt2(const LargeJet::Subjet &a, const LargeJet::Subjet &b) {
  return a.v.Perp() >= b.v.Perp();
}

PlotEJetTT::PlotEJetTT(const std::string &filename, int mode, const std::vector<std::string> &systs)
  : Plot(filename, systs), m_mode(mode),
    signal(0), background(0), isr(0), deconstruct(0), param("input_card.dat"), paramb("input_card.dat") {
  m_hSvc.create1D("lepPt", "; Lepton p_{T} [GeV] ; Events", 40, 25, 425);
  m_hSvc.create1D("lepEta", "; Lepton #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("lepPhi", "; Lepton #phi [rad] ; Events", 64, -3.2, 3.2);
  m_hSvc.create1D("jetPt", "; Jet p_{T} [GeV] ; Events", 100, 25, 1025);
  m_hSvc.create1D("jetEta", "; Jet #eta ; Events", 50, -2.5, 2.5);
  m_hSvc.create1D("jetPhi", "; Jet #phi [rad] ; Events", 64, -3.2, 3.2);

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

  m_hSvc.create1D("mu", "; <#mu> ; Events", 5, 12, 22);
  m_hSvc.create1D("npv", "; npv ; Events", 5, 10, 20);

  m_hSvc.create1D("met", "; Missing E_{T} [GeV] ; Events", 20, 0, 400);
  m_hSvc.create1D("mtw", "; m_{T} [GeV] ; Events", 10, 0, 200);

  m_hSvc.create1D("llargeJetPtPos", "; Large jet p_{T} [GeV] ; Events", 16, 350, 1150);
  m_hSvc.create1D("llargeJetPtPosCa", "; Large jet p_{T} [GeV] ; Events", 16, 350, 1150);

  m_hSvc.create1D("llargeJetMPos", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosCa", "; Large jet M ; Events", 40, 0, 400);

  m_hSvc.create1D("subjetlcpt", "; subjet lc pt ; Events", 50, 0, 50);

  // after chi cuts
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

  m_hSvc.create1D("llargeJetPtPosChiWP1", "; Large jet p_{T} [GeV] ; Events", 16, 350, 1150);
  m_hSvc.create1D("llargeJetPtPosChibWP1", "; Large jet p_{T} [GeV] ; Events", 16, 350, 1150);
  m_hSvc.create1D("llargeJetPtPosCaChiWP1", "; Large jet p_{T} [GeV] ; Events", 16, 350, 1150);
  m_hSvc.create1D("llargeJetPtPosCaChibWP1", "; Large jet p_{T} [GeV] ; Events", 16, 350, 1150);

  m_hSvc.create1D("llargeJetMPosChiWP1", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosCaChiWP1", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosChibWP1", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosCaChibWP1", "; Large jet M ; Events", 40, 0, 400);

  m_hSvc.create1D("llargeJetSDMChiWP1", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMCaChiWP1", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMChibWP1", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMCaChibWP1", "; Large jet M ; Events", 40, 0, 400);

  m_hSvc.create1D("lm12ChiWP1", "; m_{12} ; Events", 20, 0, 200);
  m_hSvc.create1D("lm23ChiWP1", "; m_{23} ; Events", 20, 0, 200);
  m_hSvc.create1D("lm123ChiWP1", "; m_{123} ; Events", 20, 0, 300);

  m_hSvc.create1D("lm12ChibWP1", "; m_{12} ; Events", 20, 0, 200);
  m_hSvc.create1D("lm23ChibWP1", "; m_{23} ; Events", 20, 0, 200);
  m_hSvc.create1D("lm123ChibWP1", "; m_{123} ; Events", 20, 0, 300);

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

  m_hSvc.create1D("llargeJetPtPosChiWP2", "; Large jet p_{T} [GeV] ; Events", 16, 350, 1150);
  m_hSvc.create1D("llargeJetPtPosChibWP2", "; Large jet p_{T} [GeV] ; Events", 16, 350, 1150);
  m_hSvc.create1D("llargeJetPtPosCaChiWP2", "; Large jet p_{T} [GeV] ; Events", 16, 350, 1150);
  m_hSvc.create1D("llargeJetPtPosCaChibWP2", "; Large jet p_{T} [GeV] ; Events", 16, 350, 1150);

  m_hSvc.create1D("llargeJetMPosChiWP2", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosCaChiWP2", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosChibWP2", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosCaChibWP2", "; Large jet M ; Events", 40, 0, 400);

  m_hSvc.create1D("llargeJetSDMChiWP2", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMCaChiWP2", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMChibWP2", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMCaChibWP2", "; Large jet M ; Events", 40, 0, 400);

  m_hSvc.create1D("lm12ChiWP2", "; m_{12} ; Events", 20, 0, 200);
  m_hSvc.create1D("lm23ChiWP2", "; m_{23} ; Events", 20, 0, 200);
  m_hSvc.create1D("lm123ChiWP2", "; m_{123} ; Events", 20, 0, 300);

  m_hSvc.create1D("lm12ChibWP2", "; m_{12} ; Events", 20, 0, 200);
  m_hSvc.create1D("lm23ChibWP2", "; m_{23} ; Events", 20, 0, 200);
  m_hSvc.create1D("lm123ChibWP2", "; m_{123} ; Events", 20, 0, 300);

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

  m_hSvc.create1D("llargeJetPtPosChiWP3", "; Large jet p_{T} [GeV] ; Events", 16, 350, 1150);
  m_hSvc.create1D("llargeJetPtPosChibWP3", "; Large jet p_{T} [GeV] ; Events", 16, 350, 1150);
  m_hSvc.create1D("llargeJetPtPosCaChiWP3", "; Large jet p_{T} [GeV] ; Events", 16, 350, 1150);
  m_hSvc.create1D("llargeJetPtPosCaChibWP3", "; Large jet p_{T} [GeV] ; Events", 16, 350, 1150);

  m_hSvc.create1D("llargeJetMPosChiWP3", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosCaChiWP3", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosChibWP3", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetMPosCaChibWP3", "; Large jet M ; Events", 40, 0, 400);

  m_hSvc.create1D("llargeJetSDMChiWP3", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMCaChiWP3", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMChibWP3", "; Large jet M ; Events", 40, 0, 400);
  m_hSvc.create1D("llargeJetSDMCaChibWP3", "; Large jet M ; Events", 40, 0, 400);

  m_hSvc.create1D("lm12ChiWP3", "; m_{12} ; Events", 20, 0, 200);
  m_hSvc.create1D("lm23ChiWP3", "; m_{23} ; Events", 20, 0, 200);
  m_hSvc.create1D("lm123ChiWP3", "; m_{123} ; Events", 20, 0, 300);

  m_hSvc.create1D("lm12ChibWP3", "; m_{12} ; Events", 20, 0, 200);
  m_hSvc.create1D("lm23ChibWP3", "; m_{23} ; Events", 20, 0, 200);
  m_hSvc.create1D("lm123ChibWP3", "; m_{123} ; Events", 20, 0, 300);

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

  m_hSvc.create1D("npvEff", "; NPV ; Events", 5, 10, 20);
  m_hSvc.create1D("npvEffCa", "; NPV ; Events", 5, 10, 20);

  m_hSvc.create1D("npvEffChiWP1", "; NPV ; Events", 5, 10, 20);
  m_hSvc.create1D("npvEffCaChiWP1", "; NPV ; Events", 5, 10, 20);
  m_hSvc.create1D("npvEffChibWP1", "; NPV ; Events", 5, 10, 20);
  m_hSvc.create1D("npvEffCaChibWP1", "; NPV ; Events", 5, 10, 20);

  m_hSvc.create1D("npvEffChiWP2", "; NPV ; Events", 5, 10, 20);
  m_hSvc.create1D("npvEffCaChiWP2", "; NPV ; Events", 5, 10, 20);
  m_hSvc.create1D("npvEffChibWP2", "; NPV ; Events", 5, 10, 20);
  m_hSvc.create1D("npvEffCaChibWP2", "; NPV ; Events", 5, 10, 20);

  m_hSvc.create1D("npvEffChiWP3", "; NPV ; Events", 5, 10, 20);
  m_hSvc.create1D("npvEffCaChiWP3", "; NPV ; Events", 5, 10, 20);
  m_hSvc.create1D("npvEffChibWP3", "; NPV ; Events", 5, 10, 20);
  m_hSvc.create1D("npvEffCaChibWP3", "; NPV ; Events", 5, 10, 20);


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

  double llb50[] = {350, 400, 450, 500, 600, 700};
  int nllb50 = 5;
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

  signal = new TopGluonModel(param);
  background = new BackgroundModel(param);
  isr = new ISRModel(param);
  deconstruct = new Deconstruct(param, *signal, *background, *isr);

  paramb["useBtag"] = 1;

  signalb = new TopGluonModel(paramb);
  backgroundb = new BackgroundModel(paramb);
  isrb = new ISRModel(paramb);
  deconstructb = new Deconstruct(paramb, *signalb, *backgroundb, *isrb);

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
  }

  TFile *f = new TFile("ptFkrWeight.root");
  m_pol = (TF1 *) f->Get("pol1")->Clone("mypol");

}

PlotEJetTT::~PlotEJetTT() {
  if (deconstruct) delete deconstruct;
  if (isr) delete isr;
  if (background) delete background;
  if (signal) delete signal;
}

std::vector<LargeJet::Subjet> PlotEJetTT::applySjUnc(const std::vector<LargeJet::Subjet> &vsj, const std::vector<LargeJet::Subjet> &vsj_unc, const std::string &s) {
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

void PlotEJetTT::run(const Event &e, double weight, double pweight, const std::string &s) {
  HistogramService *h = &m_hSvc;
  if (e.passReco()) {

    TLorentzVector l;
    l = e.triggerJet4mom();
    bool testdR = true;
    bool testdRBB = true;
    

    // only run SD for large jets passing the following cuts for basic studies
    // pT > 350 GeV makes top be contained in the large jet
    int fji = -1;
    for (int k = 0; k < e.largeJet().size(); ++k) {
      testdR = true;
      if (e.isData())
        testdR = (l.DeltaR(e.largeJet()[k].mom()) > 1.5) && ((std::fabs(l.Eta()) < 1.37) || (std::fabs(l.Eta()) > 1.52)) && (std::fabs(l.Eta()) < 2.4);
      if (e.largeJet()[k].passLoose() && e.largeJet()[k].mom().Perp() > 350e3 && testdR) {
      //if (e.largeJet()[k].passLoose() && e.largeJet()[k].mom().Perp() > 200e3 && l.DeltaR(e.largeJet()[k].mom()) > 1.5) {
        fji = k;
        break;
      }
    }
    int cafji = -1;
    for (int k = 0; k < e.largeJetBB().size(); ++k) {
      testdRBB = true;
      if (e.isData())
        testdRBB = (l.DeltaR(e.largeJetBB()[k].mom()) > 1.5) && ((std::fabs(l.Eta()) < 1.37) || (std::fabs(l.Eta()) > 1.52)) && (std::fabs(l.Eta()) < 2.4);
      if (e.largeJetBB()[k].passLoose() && e.largeJetBB()[k].mom().Perp() > 350e3 && testdRBB) {
      //if (e.largeJetBB()[k].passLoose() && e.largeJetBB()[k].mom().Perp() > 200e3 && l.DeltaR(e.largeJetBB()[k].mom()) > 1.5) {
        cafji = k;
        break;
      }
    }
    if (fji < 0 || cafji < 0) return;

    const LargeJet &lj = e.largeJet()[fji];
    const LargeJet &calj = e.largeJet()[cafji];
    if (lj.mom().DeltaR(calj.mom()) > 0.75) return;
    if (e.isData()) {
      if (lj.mom().DeltaR(l) < 1.5) return;
      if (calj.mom().DeltaR(l) < 1.5) return;
    }

    if (!e.isData()) {
      double ptWeight = m_pol->Eval(lj.mom().Perp()*1e-3);
      weight *= ptWeight;
      pweight *= ptWeight;
    }

    if (e.isData()) {
      h->h1D("lepPt", "", s)->Fill(l.Perp()*1e-3, weight);
      h->h1D("lepEta", "", s)->Fill(l.Eta(), weight);
      h->h1D("lepPhi", "", s)->Fill(l.Phi(), weight);
    }

    if (e.jet().size() >= 1) {
      const TLorentzVector &j = e.jet()[0].mom();
      h->h1D("jetPt", "", s)->Fill(j.Perp()*1e-3, weight);
      h->h1D("jetEta", "", s)->Fill(j.Eta()*1e-3, weight);
      h->h1D("jetPhi", "", s)->Fill(j.Phi(), weight);
    }

    h->h1D("met", "", s)->Fill(e.met().Perp()*1e-3, weight);
    if (e.isData()) {
      float mtw = std::sqrt(2*l.Perp()*e.met().Perp()*(1 - std::cos(l.Phi() - e.met().Phi())) )*1e-3;
      h->h1D("mtw", "", s)->Fill(mtw, weight);
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
    std::sort(vsj.begin(), vsj.end(), SortByPt2);
    sjn = vsj.size();

    int casjn = calj.subjet().size();
    std::vector<LargeJet::Subjet> cavsj = applySjUnc(calj.subjet(), calj.subjetJER(), s);
    std::sort(cavsj.begin(), cavsj.end(), SortByPt2);
    casjn = cavsj.size();
    
    for (size_t x = 0; x < sjn; ++x) h->h1D("subjetlcpt","",s)->Fill(vsj[x].lcpt*1e-3, weight);

    h->h1D("lsubjetN", "", s)->Fill(sjn, weight);
    h->h1D("lsubjetNSF", "", s)->Fill(sjn, weight);
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
    h->h1D("calsubjetNSF", "", s)->Fill(casjn, weight);
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
    int bn = 0;
    for (int k = 0; k < sjn; ++k) {
      const TLorentzVector &sj = vsj[k].v;
      fastjet::PseudoJet psj(sj.Px()*1e-3, sj.Py()*1e-3, sj.Pz()*1e-3, sj.E()*1e-3);
      psj.set_user_index(0);

      for (int z = 0; z < e.jet().size(); ++z) {
        double dr = e.jet()[z].mom().DeltaR(sj);
        if (dr < 0.2) {
          if (e.jet()[z].btag()) {
            psj.set_user_index(1);
            bn++;
          } else
            psj.set_user_index(-1);
          break;
        }
      }
      sd_subs.push_back(psj);
    }
    std::vector<fastjet::PseudoJet> sd_sorted = sorted_by_pt(sd_subs);
    if (sd_sorted.size() > 9)
      sd_sorted.erase(sd_sorted.begin() + 9, sd_sorted.begin() + sd_sorted.size());

    std::vector<float> sdb_r;
    std::vector<float> sdl_r;
    for (int k = 0; k < sd_sorted.size(); ++k) {
      const TLorentzVector sj(sd_sorted[k].px(), sd_sorted[k].py(), sd_sorted[k].pz(), sd_sorted[k].e());
      float minDr = 99;
      float minDrL = 99;
      for (int z = 0; z < e.jet().size(); ++z) {
        float mdr = e.jet()[z].mom().DeltaR(sj);
        if (mdr < minDrL) minDrL = mdr;

        if (!e.jet()[z].btag()) continue;
        if (mdr < minDr) minDr = mdr;
      }
      sdb_r.push_back(minDr);
      sdl_r.push_back(minDrL);
    }

    std::vector<fastjet::PseudoJet> casd_subs;
    int cabn = 0;
    for (int k = 0; k < casjn; ++k) {
      const TLorentzVector &sj = cavsj[k].v;
      fastjet::PseudoJet psj(sj.Px()*1e-3, sj.Py()*1e-3, sj.Pz()*1e-3, sj.E()*1e-3);
      psj.set_user_index(0);
      for (int z = 0; z < e.jet().size(); ++z) {
        double dr = e.jet()[z].mom().DeltaR(sj);
        if (dr < 0.2) {
          if (e.jet()[z].btag()) {
            psj.set_user_index(1);
            cabn++;
          } else
            psj.set_user_index(-1);
          break;
        }
      }
      casd_subs.push_back(psj);
    }
    std::vector<fastjet::PseudoJet> casd_sorted = sorted_by_pt(casd_subs);
    if (casd_sorted.size() > 9)
      casd_sorted.erase(casd_sorted.begin() + 9, casd_sorted.begin() + casd_sorted.size());

    std::vector<float> casdb_r;
    std::vector<float> casdl_r;
    for (int k = 0; k < casd_sorted.size(); ++k) {
      const TLorentzVector sj(casd_sorted[k].px(), casd_sorted[k].py(), casd_sorted[k].pz(), casd_sorted[k].e());
      float minDr = 99;
      float minDrL = 99;
      for (int z = 0; z < e.jet().size(); ++z) {
        float mdr = e.jet()[z].mom().DeltaR(sj);
        if (mdr < minDrL) minDrL = mdr;

        if (!e.jet()[z].btag()) continue;
        if (mdr < minDr) minDr = mdr;
      }
      casdb_r.push_back(minDr);
      casdl_r.push_back(minDrL);
    }

    h->h1D("lsubjetBN", "", s)->Fill(bn, weight);
    if (sdb_r.size() >= 1)
      h->h1D("lsubjet0Bdr", "", s)->Fill(sdb_r[0], weight);
    if (sdb_r.size() >= 2)
      h->h1D("lsubjet1Bdr", "", s)->Fill(sdb_r[1], weight);
    if (sdb_r.size() >= 3)
      h->h1D("lsubjet2Bdr", "", s)->Fill(sdb_r[2], weight);
    if (sdb_r.size() >= 4)
      h->h1D("lsubjet3Bdr", "", s)->Fill(sdb_r[3], weight);
    if (sdb_r.size() >= 5)
      h->h1D("lsubjet4Bdr", "", s)->Fill(sdb_r[4], weight);

    if (sdl_r.size() >= 1)
      h->h1D("lsubjet0Ldr", "", s)->Fill(sdl_r[0], weight);
    if (sdl_r.size() >= 2)
      h->h1D("lsubjet1Ldr", "", s)->Fill(sdl_r[1], weight);
    if (sdl_r.size() >= 3)
      h->h1D("lsubjet2Ldr", "", s)->Fill(sdl_r[2], weight);
    if (sdl_r.size() >= 4)
      h->h1D("lsubjet3Ldr", "", s)->Fill(sdl_r[3], weight);
    if (sdl_r.size() >= 5)
      h->h1D("lsubjet4Ldr", "", s)->Fill(sdl_r[4], weight);

    h->h1D("calsubjetBN", "", s)->Fill(cabn, weight);
    if (casdb_r.size() >= 1)
      h->h1D("calsubjet0Bdr", "", s)->Fill(casdb_r[0], weight);
    if (casdb_r.size() >= 2)
      h->h1D("calsubjet1Bdr", "", s)->Fill(casdb_r[1], weight);
    if (casdb_r.size() >= 3)
      h->h1D("calsubjet2Bdr", "", s)->Fill(casdb_r[2], weight);
    if (casdb_r.size() >= 4)
      h->h1D("calsubjet3Bdr", "", s)->Fill(casdb_r[3], weight);
    if (casdb_r.size() >= 5)
      h->h1D("calsubjet4Bdr", "", s)->Fill(casdb_r[4], weight);

    if (casdl_r.size() >= 1)
      h->h1D("calsubjet0Ldr", "", s)->Fill(casdl_r[0], weight);
    if (casdl_r.size() >= 2)
      h->h1D("calsubjet1Ldr", "", s)->Fill(casdl_r[1], weight);
    if (casdl_r.size() >= 3)
      h->h1D("calsubjet2Ldr", "", s)->Fill(casdl_r[2], weight);
    if (casdl_r.size() >= 4)
      h->h1D("calsubjet3Ldr", "", s)->Fill(casdl_r[3], weight);
    if (casdl_r.size() >= 5)
      h->h1D("calsubjet4Ldr", "", s)->Fill(casdl_r[4], weight);

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

    h->h1D("pre", "", s)->Fill(0.5, weight);
    h->h1D("preCa", "", s)->Fill(0.5, weight);
    if (lchi > CHIWP1) h->h1D("preChiWP1", "", s)->Fill(0.5, weight);
    if (lchib > CHIBWP1) h->h1D("preChibWP1", "", s)->Fill(0.5, weight);
    if (calchi > CACHIWP1) h->h1D("preCaChiWP1", "", s)->Fill(0.5, weight);
    if (calchib > CACHIBWP1) h->h1D("preCaChibWP1", "", s)->Fill(0.5, weight);

    if (lchi > CHIWP2) h->h1D("preChiWP2", "", s)->Fill(0.5, weight);
    if (lchib > CHIBWP2) h->h1D("preChibWP2", "", s)->Fill(0.5, weight);
    if (calchi > CACHIWP2) h->h1D("preCaChiWP2", "", s)->Fill(0.5, weight);
    if (calchib > CACHIBWP2) h->h1D("preCaChibWP2", "", s)->Fill(0.5, weight);

    if (lchi > CHIWP3) h->h1D("preChiWP3", "", s)->Fill(0.5, weight);
    if (lchib > CHIBWP3) h->h1D("preChibWP3", "", s)->Fill(0.5, weight);
    if (calchi > CACHIWP3) h->h1D("preCaChiWP3", "", s)->Fill(0.5, weight);
    if (calchib > CACHIBWP3) h->h1D("preCaChibWP3", "", s)->Fill(0.5, weight);

    if (std::fabs(lj.mom().Eta()) < 0.7) {
      h->h1D("preLE", "", s)->Fill(0.5, weight);
      if (lchi > CHIWP1) h->h1D("preLEChiWP1", "", s)->Fill(0.5, weight);
      if (lchib > CHIBWP1) h->h1D("preLEChibWP1", "", s)->Fill(0.5, weight);
      if (lchi > CHIWP2) h->h1D("preLEChiWP2", "", s)->Fill(0.5, weight);
      if (lchib > CHIBWP2) h->h1D("preLEChibWP2", "", s)->Fill(0.5, weight);
      if (lchi > CHIWP3) h->h1D("preLEChiWP3", "", s)->Fill(0.5, weight);
      if (lchib > CHIBWP3) h->h1D("preLEChibWP3", "", s)->Fill(0.5, weight);
    } else {
      h->h1D("preHE", "", s)->Fill(0.5, weight);
      if (lchi > CHIWP1) h->h1D("preHEChiWP1", "", s)->Fill(0.5, weight);
      if (lchib > CHIBWP1) h->h1D("preHEChibWP1", "", s)->Fill(0.5, weight);
      if (lchi > CHIWP2) h->h1D("preHEChiWP2", "", s)->Fill(0.5, weight);
      if (lchib > CHIBWP2) h->h1D("preHEChibWP2", "", s)->Fill(0.5, weight);
      if (lchi > CHIWP3) h->h1D("preHEChiWP3", "", s)->Fill(0.5, weight);
      if (lchib > CHIBWP3) h->h1D("preHEChibWP3", "", s)->Fill(0.5, weight);
    }

    if (std::fabs(calj.mom().Eta()) < 0.7) {
      h->h1D("preLECa", "", s)->Fill(0.5, weight);
      if (calchi > CACHIWP1) h->h1D("preLECaChiWP1", "", s)->Fill(0.5, weight);
      if (calchib > CACHIBWP1) h->h1D("preLECaChibWP1", "", s)->Fill(0.5, weight);
      if (calchi > CACHIWP2) h->h1D("preLECaChiWP2", "", s)->Fill(0.5, weight);
      if (calchib > CACHIBWP2) h->h1D("preLECaChibWP2", "", s)->Fill(0.5, weight);
      if (calchi > CACHIWP3) h->h1D("preLECaChiWP3", "", s)->Fill(0.5, weight);
      if (calchib > CACHIBWP3) h->h1D("preLECaChibWP3", "", s)->Fill(0.5, weight);
    } else {
      h->h1D("preHECa", "", s)->Fill(0.5, weight);
      if (calchi > CACHIWP1) h->h1D("preHECaChiWP1", "", s)->Fill(0.5, weight);
      if (calchib > CACHIBWP1) h->h1D("preHECaChibWP1", "", s)->Fill(0.5, weight);
      if (calchi > CACHIWP2) h->h1D("preHECaChiWP2", "", s)->Fill(0.5, weight);
      if (calchib > CACHIBWP2) h->h1D("preHECaChibWP2", "", s)->Fill(0.5, weight);
      if (calchi > CACHIWP3) h->h1D("preHECaChiWP3", "", s)->Fill(0.5, weight);
      if (calchib > CACHIBWP3) h->h1D("preHECaChibWP3", "", s)->Fill(0.5, weight);
    }

    h->h1D("muEff", "", s)->Fill(e.mu(), weight);
    h->h1D("muEffCa", "", s)->Fill(e.mu(), weight);
    if (lchi > CHIWP1) h->h1D("muEffChiWP1", "", s)->Fill(e.mu(), weight);
    if (calchi > CACHIWP1) h->h1D("muEffCaChiWP1", "", s)->Fill(e.mu(), weight);
    if (lchib > CHIBWP1) h->h1D("muEffChibWP1", "", s)->Fill(e.mu(), weight);
    if (calchib > CACHIBWP1) h->h1D("muEffCaChibWP1", "", s)->Fill(e.mu(), weight);

    if (lchi > CHIWP2) h->h1D("muEffChiWP2", "", s)->Fill(e.mu(), weight);
    if (calchi > CACHIWP2) h->h1D("muEffCaChiWP2", "", s)->Fill(e.mu(), weight);
    if (lchib > CHIBWP2) h->h1D("muEffChibWP2", "", s)->Fill(e.mu(), weight);
    if (calchib > CACHIBWP2) h->h1D("muEffCaChibWP2", "", s)->Fill(e.mu(), weight);

    h->h1D("npvEff", "", s)->Fill(e.npv(), weight);
    h->h1D("npvEffCa", "", s)->Fill(e.npv(), weight);
    if (lchi > CHIWP1) h->h1D("npvEffChiWP1", "", s)->Fill(e.npv(), weight);
    if (calchi > CACHIWP1) h->h1D("npvEffCaChiWP1", "", s)->Fill(e.npv(), weight);
    if (lchib > CHIBWP1) h->h1D("npvEffChibWP1", "", s)->Fill(e.npv(), weight);
    if (calchib > CACHIBWP1) h->h1D("npvEffCaChibWP1", "", s)->Fill(e.npv(), weight);

    if (lchi > CHIWP2) h->h1D("npvEffChiWP2", "", s)->Fill(e.npv(), weight);
    if (calchi > CACHIWP2) h->h1D("npvEffCaChiWP2", "", s)->Fill(e.npv(), weight);
    if (lchib > CHIBWP2) h->h1D("npvEffChibWP2", "", s)->Fill(e.npv(), weight);
    if (calchib > CACHIBWP2) h->h1D("npvEffCaChibWP2", "", s)->Fill(e.npv(), weight);

    if (lchi > CHIWP3) h->h1D("npvEffChiWP3", "", s)->Fill(e.npv(), weight);
    if (calchi > CACHIWP3) h->h1D("npvEffCaChiWP3", "", s)->Fill(e.npv(), weight);
    if (lchib > CHIBWP3) h->h1D("npvEffChibWP3", "", s)->Fill(e.npv(), weight);
    if (calchib > CACHIBWP3) h->h1D("npvEffCaChibWP3", "", s)->Fill(e.npv(), weight);

    h->h1D("llargeJetPtEff", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
    h->h1D("llargeJetPtEffCa", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
    if (lchi > CHIWP1) h->h1D("llargeJetPtEffChiWP1", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
    if (calchi > CACHIWP1) h->h1D("llargeJetPtEffCaChiWP1", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
    if (lchib > CHIBWP1) h->h1D("llargeJetPtEffChibWP1", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
    if (calchib > CACHIBWP1) h->h1D("llargeJetPtEffCaChibWP1", "", s)->Fill(calj.mom().Perp()*1e-3, weight);

    if (lchi > CHIWP2) h->h1D("llargeJetPtEffChiWP2", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
    if (calchi > CACHIWP2) h->h1D("llargeJetPtEffCaChiWP2", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
    if (lchib > CHIBWP2) h->h1D("llargeJetPtEffChibWP2", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
    if (calchib > CACHIBWP2) h->h1D("llargeJetPtEffCaChibWP2", "", s)->Fill(calj.mom().Perp()*1e-3, weight);

    if (lchi > CHIWP3) h->h1D("llargeJetPtEffChiWP3", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
    if (calchi > CACHIWP3) h->h1D("llargeJetPtEffCaChiWP3", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
    if (lchib > CHIBWP3) h->h1D("llargeJetPtEffChibWP3", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
    if (calchib > CACHIBWP3) h->h1D("llargeJetPtEffCaChibWP3", "", s)->Fill(calj.mom().Perp()*1e-3, weight);

    if (std::fabs(lj.mom().Eta()) < 0.7) {
      h->h1D("llargeJetPtLEEff", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
      if (lchi > CHIWP1) h->h1D("llargeJetPtLEEffChiWP1", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
      if (lchib > CHIBWP1) h->h1D("llargeJetPtLEEffChibWP1", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
      if (lchi > CHIWP2) h->h1D("llargeJetPtLEEffChiWP2", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
      if (lchib > CHIBWP2) h->h1D("llargeJetPtLEEffChibWP2", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
      if (lchi > CHIWP3) h->h1D("llargeJetPtLEEffChiWP3", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
      if (lchib > CHIBWP3) h->h1D("llargeJetPtLEEffChibWP3", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
    } else {
      h->h1D("llargeJetPtHEEff", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
      if (lchi > CHIWP1) h->h1D("llargeJetPtHEEffChiWP1", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
      if (lchib > CHIBWP1) h->h1D("llargeJetPtHEEffChibWP1", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
      if (lchi > CHIWP2) h->h1D("llargeJetPtHEEffChiWP2", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
      if (lchib > CHIBWP2) h->h1D("llargeJetPtHEEffChibWP2", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
      if (lchi > CHIWP3) h->h1D("llargeJetPtHEEffChiWP3", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
      if (lchib > CHIBWP3) h->h1D("llargeJetPtHEEffChibWP3", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
    }
    if (std::fabs(calj.mom().Eta()) < 0.7) {
      h->h1D("llargeJetPtLEEffCa", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
      if (calchi > CACHIWP1) h->h1D("llargeJetPtLEEffCaChiWP1", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
      if (calchib > CACHIBWP1) h->h1D("llargeJetPtLEEffCaChibWP1", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
      if (calchi > CACHIWP2) h->h1D("llargeJetPtLEEffCaChiWP2", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
      if (calchib > CACHIBWP2) h->h1D("llargeJetPtLEEffCaChibWP2", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
      if (calchi > CACHIWP3) h->h1D("llargeJetPtLEEffCaChiWP3", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
      if (calchib > CACHIBWP3) h->h1D("llargeJetPtLEEffCaChibWP3", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
    } else {
      h->h1D("llargeJetPtHEEffCa", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
      if (calchi > CACHIWP1) h->h1D("llargeJetPtHEEffCaChiWP1", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
      if (calchib > CACHIBWP1) h->h1D("llargeJetPtHEEffCaChibWP1", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
      if (calchi > CACHIWP2) h->h1D("llargeJetPtHEEffCaChiWP2", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
      if (calchib > CACHIBWP2) h->h1D("llargeJetPtHEEffCaChibWP2", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
      if (calchi > CACHIWP3) h->h1D("llargeJetPtHEEffCaChiWP3", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
      if (calchib > CACHIBWP3) h->h1D("llargeJetPtHEEffCaChibWP3", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
    }


    h->h1D("llargeJetPtPos", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
    h->h1D("llargeJetPtPosCa", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
    if (lchi > CHIWP1) h->h1D("llargeJetPtPosChiWP1", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
    if (calchi > CACHIWP1) h->h1D("llargeJetPtPosCaChiWP1", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
    if (lchib > CHIBWP1) h->h1D("llargeJetPtPosChibWP1", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
    if (calchib > CACHIBWP1) h->h1D("llargeJetPtPosCaChibWP1", "", s)->Fill(calj.mom().Perp()*1e-3, weight);

    if (lchi > CHIWP2) h->h1D("llargeJetPtPosChiWP2", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
    if (calchi > CACHIWP2) h->h1D("llargeJetPtPosCaChiWP2", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
    if (lchib > CHIBWP2) h->h1D("llargeJetPtPosChibWP2", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
    if (calchib > CACHIBWP2) h->h1D("llargeJetPtPosCaChibWP2", "", s)->Fill(calj.mom().Perp()*1e-3, weight);

    if (lchi > CHIWP3) h->h1D("llargeJetPtPosChiWP3", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
    if (calchi > CACHIWP3) h->h1D("llargeJetPtPosCaChiWP3", "", s)->Fill(calj.mom().Perp()*1e-3, weight);
    if (lchib > CHIBWP3) h->h1D("llargeJetPtPosChibWP3", "", s)->Fill(lj.mom().Perp()*1e-3, weight);
    if (calchib > CACHIBWP3) h->h1D("llargeJetPtPosCaChibWP3", "", s)->Fill(calj.mom().Perp()*1e-3, weight);

    h->h1D("llargeJetMPos", "", s)->Fill(lj.mom().M()*1e-3, weight);
    h->h1D("llargeJetMPosCa", "", s)->Fill(calj.mom().M()*1e-3, weight);
    if (lchi > CHIWP1) h->h1D("llargeJetMPosChiWP1", "", s)->Fill(lj.mom().M()*1e-3, weight);
    if (calchi > CACHIWP1) h->h1D("llargeJetMPosCaChiWP1", "", s)->Fill(calj.mom().M()*1e-3, weight);
    if (lchib > CHIBWP1) h->h1D("llargeJetMPosChibWP1", "", s)->Fill(lj.mom().M()*1e-3, weight);
    if (calchib > CACHIBWP1) h->h1D("llargeJetMPosCaChibWP1", "", s)->Fill(calj.mom().M()*1e-3, weight);

    if (lchi > CHIWP2) h->h1D("llargeJetMPosChiWP2", "", s)->Fill(lj.mom().M()*1e-3, weight);
    if (calchi > CACHIWP2) h->h1D("llargeJetMPosCaChiWP2", "", s)->Fill(calj.mom().M()*1e-3, weight);
    if (lchib > CHIBWP2) h->h1D("llargeJetMPosChibWP2", "", s)->Fill(lj.mom().M()*1e-3, weight);
    if (calchib > CACHIBWP2) h->h1D("llargeJetMPosCaChibWP2", "", s)->Fill(calj.mom().M()*1e-3, weight);

    if (lchi > CHIWP3) h->h1D("llargeJetMPosChiWP3", "", s)->Fill(lj.mom().M()*1e-3, weight);
    if (calchi > CACHIWP3) h->h1D("llargeJetMPosCaChiWP3", "", s)->Fill(calj.mom().M()*1e-3, weight);
    if (lchib > CHIBWP3) h->h1D("llargeJetMPosChibWP3", "", s)->Fill(lj.mom().M()*1e-3, weight);
    if (calchib > CACHIBWP3) h->h1D("llargeJetMPosCaChibWP3", "", s)->Fill(calj.mom().M()*1e-3, weight);

    if (lchi > CHIWP1) h->h1D("llargeJetSDMChiWP1", "", s)->Fill(chi_top.m(), weight);
    if (calchi > CACHIWP1) h->h1D("llargeJetSDMCaChiWP1", "", s)->Fill(cachi_top.m(), weight);
    if (lchib > CHIBWP1) h->h1D("llargeJetSDMChibWP1", "", s)->Fill(chib_top.m(), weight);
    if (calchib > CACHIBWP1) h->h1D("llargeJetSDMCaChibWP1", "", s)->Fill(cachib_top.m(), weight);

    if (lchi > CHIWP2) h->h1D("llargeJetSDMChiWP2", "", s)->Fill(chi_top.m(), weight);
    if (calchi > CACHIWP2) h->h1D("llargeJetSDMCaChiWP2", "", s)->Fill(cachi_top.m(), weight);
    if (lchib > CHIBWP2) h->h1D("llargeJetSDMChibWP2", "", s)->Fill(chib_top.m(), weight);
    if (calchib > CACHIBWP2) h->h1D("llargeJetSDMCaChibWP2", "", s)->Fill(cachib_top.m(), weight);

    if (lchi > CHIWP3) h->h1D("llargeJetSDMChiWP3", "", s)->Fill(chi_top.m(), weight);
    if (calchi > CACHIWP3) h->h1D("llargeJetSDMCaChiWP3", "", s)->Fill(cachi_top.m(), weight);
    if (lchib > CHIBWP3) h->h1D("llargeJetSDMChibWP3", "", s)->Fill(chib_top.m(), weight);
    if (calchib > CACHIBWP3) h->h1D("llargeJetSDMCaChibWP3", "", s)->Fill(cachib_top.m(), weight);

    if (lchi > CHIWP1) {
      h->h1D("lsubjetNChiWP1", "", s)->Fill(sjn, weight);
      if (std::fabs(lj.mom().Eta()) < 0.7) {
        h->h1D("lsubjetNLEChiWP1", "", s)->Fill(sjn, weight);
        if (lj.mom().Perp() > 500e3 && lj.mom().Perp() < 550e3)
          h->h1D("lsubjetNLETestChiWP1", "", s)->Fill(sjn, weight);
      } else {
        h->h1D("lsubjetNHEChiWP1", "", s)->Fill(sjn, weight);
      }
      if (sjn >= 1)
        h->h1D("lsubjet0PtChiWP1", "", s)->Fill(vsj[0].v.Perp()*1e-3, weight);
      if (sjn >= 2)
        h->h1D("lsubjet1PtChiWP1", "", s)->Fill(vsj[1].v.Perp()*1e-3, weight);
      if (sjn >= 3)
        h->h1D("lsubjet2PtChiWP1", "", s)->Fill(vsj[2].v.Perp()*1e-3, weight);
      if (sjn >= 4)
        h->h1D("lsubjet3PtChiWP1", "", s)->Fill(vsj[3].v.Perp()*1e-3, weight);
      if (sjn >= 5)
        h->h1D("lsubjet4PtChiWP1", "", s)->Fill(vsj[4].v.Perp()*1e-3, weight);

      if (sjn >= 1)
        h->h1D("lsubjet0EtaChiWP1", "", s)->Fill(vsj[0].v.Eta(), weight);
      if (sjn >= 2)
        h->h1D("lsubjet1EtaChiWP1", "", s)->Fill(vsj[1].v.Eta(), weight);
      if (sjn >= 3)
        h->h1D("lsubjet2EtaChiWP1", "", s)->Fill(vsj[2].v.Eta(), weight);
      if (sjn >= 4)
        h->h1D("lsubjet3EtaChiWP1", "", s)->Fill(vsj[3].v.Eta(), weight);
      if (sjn >= 5)
        h->h1D("lsubjet4EtaChiWP1", "", s)->Fill(vsj[4].v.Eta(), weight);

      if (sjn >= 2)
        h->h1D("lm12ChiWP1", "", s)->Fill((vsj[0].v + vsj[1].v).M()*1e-3, weight);
      if (sjn >= 3) {
        h->h1D("lm23ChiWP1", "", s)->Fill((vsj[1].v + vsj[2].v).M()*1e-3, weight);
        h->h1D("lm123ChiWP1", "", s)->Fill((vsj[0].v + vsj[1].v + vsj[2].v).M()*1e-3, weight);
      }
    }

    if (lchi > CHIWP2) {
      h->h1D("lsubjetNChiWP2", "", s)->Fill(sjn, weight);
      if (std::fabs(lj.mom().Eta()) < 0.7) {
        h->h1D("lsubjetNLEChiWP2", "", s)->Fill(sjn, weight);
        if (lj.mom().Perp() > 500e3 && lj.mom().Perp() < 550e3)
          h->h1D("lsubjetNLETestChiWP2", "", s)->Fill(sjn, weight);
      } else {
        h->h1D("lsubjetNHEChiWP2", "", s)->Fill(sjn, weight);
      }
      if (sjn >= 1)
        h->h1D("lsubjet0PtChiWP2", "", s)->Fill(vsj[0].v.Perp()*1e-3, weight);
      if (sjn >= 2)
        h->h1D("lsubjet1PtChiWP2", "", s)->Fill(vsj[1].v.Perp()*1e-3, weight);
      if (sjn >= 3)
        h->h1D("lsubjet2PtChiWP2", "", s)->Fill(vsj[2].v.Perp()*1e-3, weight);
      if (sjn >= 4)
        h->h1D("lsubjet3PtChiWP2", "", s)->Fill(vsj[3].v.Perp()*1e-3, weight);
      if (sjn >= 5)
        h->h1D("lsubjet4PtChiWP2", "", s)->Fill(vsj[4].v.Perp()*1e-3, weight);

      if (sjn >= 1)
        h->h1D("lsubjet0EtaChiWP2", "", s)->Fill(vsj[0].v.Eta(), weight);
      if (sjn >= 2)
        h->h1D("lsubjet1EtaChiWP2", "", s)->Fill(vsj[1].v.Eta(), weight);
      if (sjn >= 3)
        h->h1D("lsubjet2EtaChiWP2", "", s)->Fill(vsj[2].v.Eta(), weight);
      if (sjn >= 4)
        h->h1D("lsubjet3EtaChiWP2", "", s)->Fill(vsj[3].v.Eta(), weight);
      if (sjn >= 5)
        h->h1D("lsubjet4EtaChiWP2", "", s)->Fill(vsj[4].v.Eta(), weight);

      if (sjn >= 2)
        h->h1D("lm12ChiWP2", "", s)->Fill((vsj[0].v + vsj[1].v).M()*1e-3, weight);
      if (sjn >= 3) {
        h->h1D("lm23ChiWP2", "", s)->Fill((vsj[1].v + vsj[2].v).M()*1e-3, weight);
        h->h1D("lm123ChiWP2", "", s)->Fill((vsj[0].v + vsj[1].v + vsj[2].v).M()*1e-3, weight);
      }
    }

    if (lchi > CHIWP3) {
      h->h1D("lsubjetNChiWP3", "", s)->Fill(sjn, weight);
      if (std::fabs(lj.mom().Eta()) < 0.7) {
        h->h1D("lsubjetNLEChiWP3", "", s)->Fill(sjn, weight);
        if (lj.mom().Perp() > 500e3 && lj.mom().Perp() < 550e3)
          h->h1D("lsubjetNLETestChiWP3", "", s)->Fill(sjn, weight);
      } else {
        h->h1D("lsubjetNHEChiWP3", "", s)->Fill(sjn, weight);
      }
      if (sjn >= 1)
        h->h1D("lsubjet0PtChiWP3", "", s)->Fill(vsj[0].v.Perp()*1e-3, weight);
      if (sjn >= 2)
        h->h1D("lsubjet1PtChiWP3", "", s)->Fill(vsj[1].v.Perp()*1e-3, weight);
      if (sjn >= 3)
        h->h1D("lsubjet2PtChiWP3", "", s)->Fill(vsj[2].v.Perp()*1e-3, weight);
      if (sjn >= 4)
        h->h1D("lsubjet3PtChiWP3", "", s)->Fill(vsj[3].v.Perp()*1e-3, weight);
      if (sjn >= 5)
        h->h1D("lsubjet4PtChiWP3", "", s)->Fill(vsj[4].v.Perp()*1e-3, weight);

      if (sjn >= 1)
        h->h1D("lsubjet0EtaChiWP3", "", s)->Fill(vsj[0].v.Eta(), weight);
      if (sjn >= 2)
        h->h1D("lsubjet1EtaChiWP3", "", s)->Fill(vsj[1].v.Eta(), weight);
      if (sjn >= 3)
        h->h1D("lsubjet2EtaChiWP3", "", s)->Fill(vsj[2].v.Eta(), weight);
      if (sjn >= 4)
        h->h1D("lsubjet3EtaChiWP3", "", s)->Fill(vsj[3].v.Eta(), weight);
      if (sjn >= 5)
        h->h1D("lsubjet4EtaChiWP3", "", s)->Fill(vsj[4].v.Eta(), weight);

      if (sjn >= 2)
        h->h1D("lm12ChiWP3", "", s)->Fill((vsj[0].v + vsj[1].v).M()*1e-3, weight);
      if (sjn >= 3) {
        h->h1D("lm23ChiWP3", "", s)->Fill((vsj[1].v + vsj[2].v).M()*1e-3, weight);
        h->h1D("lm123ChiWP3", "", s)->Fill((vsj[0].v + vsj[1].v + vsj[2].v).M()*1e-3, weight);
      }
    }

    if (calchi > CACHIWP1) {
      h->h1D("calsubjetNChiWP1", "", s)->Fill(casjn, weight);
      if (std::fabs(calj.mom().Eta()) < 0.7) {
        h->h1D("calsubjetNLEChiWP1", "", s)->Fill(casjn, weight);
      } else {
        h->h1D("calsubjetNHEChiWP1", "", s)->Fill(casjn, weight);
      }
      if (casjn >= 1)
        h->h1D("calsubjet0PtChiWP1", "", s)->Fill(cavsj[0].v.Perp()*1e-3, weight);
      if (casjn >= 2)
        h->h1D("calsubjet1PtChiWP1", "", s)->Fill(cavsj[1].v.Perp()*1e-3, weight);
      if (casjn >= 3)
        h->h1D("calsubjet2PtChiWP1", "", s)->Fill(cavsj[2].v.Perp()*1e-3, weight);
      if (casjn >= 4)
        h->h1D("calsubjet3PtChiWP1", "", s)->Fill(cavsj[3].v.Perp()*1e-3, weight);
      if (casjn >= 5)
        h->h1D("calsubjet4PtChiWP1", "", s)->Fill(cavsj[4].v.Perp()*1e-3, weight);
      if (casjn >= 1)
        h->h1D("calsubjet0EtaChiWP1", "", s)->Fill(cavsj[0].v.Eta(), weight);
      if (casjn >= 2)
        h->h1D("calsubjet1EtaChiWP1", "", s)->Fill(cavsj[1].v.Eta(), weight);
      if (casjn >= 3)
        h->h1D("calsubjet2EtaChiWP1", "", s)->Fill(cavsj[2].v.Eta(), weight);
      if (casjn >= 4)
        h->h1D("calsubjet3EtaChiWP1", "", s)->Fill(cavsj[3].v.Eta(), weight);
      if (casjn >= 5)
        h->h1D("calsubjet4EtaChiWP1", "", s)->Fill(cavsj[4].v.Eta(), weight);
    }

    if (calchi > CACHIWP2) {
      h->h1D("calsubjetNChiWP2", "", s)->Fill(casjn, weight);
      if (std::fabs(calj.mom().Eta()) < 0.7) {
        h->h1D("calsubjetNLEChiWP2", "", s)->Fill(casjn, weight);
      } else {
        h->h1D("calsubjetNHEChiWP2", "", s)->Fill(casjn, weight);
      }
      if (casjn >= 1)
        h->h1D("calsubjet0PtChiWP2", "", s)->Fill(cavsj[0].v.Perp()*1e-3, weight);
      if (casjn >= 2)
        h->h1D("calsubjet1PtChiWP2", "", s)->Fill(cavsj[1].v.Perp()*1e-3, weight);
      if (casjn >= 3)
        h->h1D("calsubjet2PtChiWP2", "", s)->Fill(cavsj[2].v.Perp()*1e-3, weight);
      if (casjn >= 4)
        h->h1D("calsubjet3PtChiWP2", "", s)->Fill(cavsj[3].v.Perp()*1e-3, weight);
      if (casjn >= 5)
        h->h1D("calsubjet4PtChiWP2", "", s)->Fill(cavsj[4].v.Perp()*1e-3, weight);
      if (casjn >= 1)
        h->h1D("calsubjet0EtaChiWP2", "", s)->Fill(cavsj[0].v.Eta(), weight);
      if (casjn >= 2)
        h->h1D("calsubjet1EtaChiWP2", "", s)->Fill(cavsj[1].v.Eta(), weight);
      if (casjn >= 3)
        h->h1D("calsubjet2EtaChiWP2", "", s)->Fill(cavsj[2].v.Eta(), weight);
      if (casjn >= 4)
        h->h1D("calsubjet3EtaChiWP2", "", s)->Fill(cavsj[3].v.Eta(), weight);
      if (casjn >= 5)
        h->h1D("calsubjet4EtaChiWP2", "", s)->Fill(cavsj[4].v.Eta(), weight);
    }

    if (calchi > CACHIWP3) {
      h->h1D("calsubjetNChiWP3", "", s)->Fill(casjn, weight);
      if (std::fabs(calj.mom().Eta()) < 0.7) {
        h->h1D("calsubjetNLEChiWP3", "", s)->Fill(casjn, weight);
      } else {
        h->h1D("calsubjetNHEChiWP3", "", s)->Fill(casjn, weight);
      }
      if (casjn >= 1)
        h->h1D("calsubjet0PtChiWP3", "", s)->Fill(cavsj[0].v.Perp()*1e-3, weight);
      if (casjn >= 2)
        h->h1D("calsubjet1PtChiWP3", "", s)->Fill(cavsj[1].v.Perp()*1e-3, weight);
      if (casjn >= 3)
        h->h1D("calsubjet2PtChiWP3", "", s)->Fill(cavsj[2].v.Perp()*1e-3, weight);
      if (casjn >= 4)
        h->h1D("calsubjet3PtChiWP3", "", s)->Fill(cavsj[3].v.Perp()*1e-3, weight);
      if (casjn >= 5)
        h->h1D("calsubjet4PtChiWP3", "", s)->Fill(cavsj[4].v.Perp()*1e-3, weight);
      if (casjn >= 1)
        h->h1D("calsubjet0EtaChiWP3", "", s)->Fill(cavsj[0].v.Eta(), weight);
      if (casjn >= 2)
        h->h1D("calsubjet1EtaChiWP3", "", s)->Fill(cavsj[1].v.Eta(), weight);
      if (casjn >= 3)
        h->h1D("calsubjet2EtaChiWP3", "", s)->Fill(cavsj[2].v.Eta(), weight);
      if (casjn >= 4)
        h->h1D("calsubjet3EtaChiWP3", "", s)->Fill(cavsj[3].v.Eta(), weight);
      if (casjn >= 5)
        h->h1D("calsubjet4EtaChiWP3", "", s)->Fill(cavsj[4].v.Eta(), weight);
    }

    if (lchib > CHIBWP1) {
      h->h1D("lsubjetNChibWP1", "", s)->Fill(sjn, weight);
      if (std::fabs(lj.mom().Eta()) < 0.7) {
        h->h1D("lsubjetNLEChibWP1", "", s)->Fill(sjn, weight);
      } else {
        h->h1D("lsubjetNHEChibWP1", "", s)->Fill(sjn, weight);
      }
      if (sjn >= 1)
        h->h1D("lsubjet0PtChibWP1", "", s)->Fill(vsj[0].v.Perp()*1e-3, weight);
      if (sjn >= 2)
        h->h1D("lsubjet1PtChibWP1", "", s)->Fill(vsj[1].v.Perp()*1e-3, weight);
      if (sjn >= 3)
        h->h1D("lsubjet2PtChibWP1", "", s)->Fill(vsj[2].v.Perp()*1e-3, weight);
      if (sjn >= 4)
        h->h1D("lsubjet3PtChibWP1", "", s)->Fill(vsj[3].v.Perp()*1e-3, weight);
      if (sjn >= 5)
        h->h1D("lsubjet4PtChibWP1", "", s)->Fill(vsj[4].v.Perp()*1e-3, weight);

      if (sjn >= 1)
        h->h1D("lsubjet0EtaChibWP1", "", s)->Fill(vsj[0].v.Eta(), weight);
      if (sjn >= 2)
        h->h1D("lsubjet1EtaChibWP1", "", s)->Fill(vsj[1].v.Eta(), weight);
      if (sjn >= 3)
        h->h1D("lsubjet2EtaChibWP1", "", s)->Fill(vsj[2].v.Eta(), weight);
      if (sjn >= 4)
        h->h1D("lsubjet3EtaChibWP1", "", s)->Fill(vsj[3].v.Eta(), weight);
      if (sjn >= 5)
        h->h1D("lsubjet4EtaChibWP1", "", s)->Fill(vsj[4].v.Eta(), weight);

      if (sjn >= 2)
        h->h1D("lm12ChibWP1", "", s)->Fill((vsj[0].v + vsj[1].v).M()*1e-3, weight);
      if (sjn >= 3) {
        h->h1D("lm23ChibWP1", "", s)->Fill((vsj[1].v + vsj[2].v).M()*1e-3, weight);
        h->h1D("lm123ChibWP1", "", s)->Fill((vsj[0].v + vsj[1].v + vsj[2].v).M()*1e-3, weight);
      }
    }

    if (lchib > CHIBWP2) {
      h->h1D("lsubjetNChibWP2", "", s)->Fill(sjn, weight);
      if (std::fabs(lj.mom().Eta()) < 0.7) {
        h->h1D("lsubjetNLEChibWP2", "", s)->Fill(sjn, weight);
      } else {
        h->h1D("lsubjetNHEChibWP2", "", s)->Fill(sjn, weight);
      }
      if (sjn >= 1)
        h->h1D("lsubjet0PtChibWP2", "", s)->Fill(vsj[0].v.Perp()*1e-3, weight);
      if (sjn >= 2)
        h->h1D("lsubjet1PtChibWP2", "", s)->Fill(vsj[1].v.Perp()*1e-3, weight);
      if (sjn >= 3)
        h->h1D("lsubjet2PtChibWP2", "", s)->Fill(vsj[2].v.Perp()*1e-3, weight);
      if (sjn >= 4)
        h->h1D("lsubjet3PtChibWP2", "", s)->Fill(vsj[3].v.Perp()*1e-3, weight);
      if (sjn >= 5)
        h->h1D("lsubjet4PtChibWP2", "", s)->Fill(vsj[4].v.Perp()*1e-3, weight);

      if (sjn >= 1)
        h->h1D("lsubjet0EtaChibWP2", "", s)->Fill(vsj[0].v.Eta(), weight);
      if (sjn >= 2)
        h->h1D("lsubjet1EtaChibWP2", "", s)->Fill(vsj[1].v.Eta(), weight);
      if (sjn >= 3)
        h->h1D("lsubjet2EtaChibWP2", "", s)->Fill(vsj[2].v.Eta(), weight);
      if (sjn >= 4)
        h->h1D("lsubjet3EtaChibWP2", "", s)->Fill(vsj[3].v.Eta(), weight);
      if (sjn >= 5)
        h->h1D("lsubjet4EtaChibWP2", "", s)->Fill(vsj[4].v.Eta(), weight);

      if (sjn >= 2)
        h->h1D("lm12ChibWP2", "", s)->Fill((vsj[0].v + vsj[1].v).M()*1e-3, weight);
      if (sjn >= 3) {
        h->h1D("lm23ChibWP2", "", s)->Fill((vsj[1].v + vsj[2].v).M()*1e-3, weight);
        h->h1D("lm123ChibWP2", "", s)->Fill((vsj[0].v + vsj[1].v + vsj[2].v).M()*1e-3, weight);
      }
    }

    if (lchib > CHIBWP3) {
      h->h1D("lsubjetNChibWP3", "", s)->Fill(sjn, weight);
      if (std::fabs(lj.mom().Eta()) < 0.7) {
        h->h1D("lsubjetNLEChibWP3", "", s)->Fill(sjn, weight);
      } else {
        h->h1D("lsubjetNHEChibWP3", "", s)->Fill(sjn, weight);
      }
      if (sjn >= 1)
        h->h1D("lsubjet0PtChibWP3", "", s)->Fill(vsj[0].v.Perp()*1e-3, weight);
      if (sjn >= 2)
        h->h1D("lsubjet1PtChibWP3", "", s)->Fill(vsj[1].v.Perp()*1e-3, weight);
      if (sjn >= 3)
        h->h1D("lsubjet2PtChibWP3", "", s)->Fill(vsj[2].v.Perp()*1e-3, weight);
      if (sjn >= 4)
        h->h1D("lsubjet3PtChibWP3", "", s)->Fill(vsj[3].v.Perp()*1e-3, weight);
      if (sjn >= 5)
        h->h1D("lsubjet4PtChibWP3", "", s)->Fill(vsj[4].v.Perp()*1e-3, weight);

      if (sjn >= 1)
        h->h1D("lsubjet0EtaChibWP3", "", s)->Fill(vsj[0].v.Eta(), weight);
      if (sjn >= 2)
        h->h1D("lsubjet1EtaChibWP3", "", s)->Fill(vsj[1].v.Eta(), weight);
      if (sjn >= 3)
        h->h1D("lsubjet2EtaChibWP3", "", s)->Fill(vsj[2].v.Eta(), weight);
      if (sjn >= 4)
        h->h1D("lsubjet3EtaChibWP3", "", s)->Fill(vsj[3].v.Eta(), weight);
      if (sjn >= 5)
        h->h1D("lsubjet4EtaChibWP3", "", s)->Fill(vsj[4].v.Eta(), weight);

      if (sjn >= 2)
        h->h1D("lm12ChibWP3", "", s)->Fill((vsj[0].v + vsj[1].v).M()*1e-3, weight);
      if (sjn >= 3) {
        h->h1D("lm23ChibWP3", "", s)->Fill((vsj[1].v + vsj[2].v).M()*1e-3, weight);
        h->h1D("lm123ChibWP3", "", s)->Fill((vsj[0].v + vsj[1].v + vsj[2].v).M()*1e-3, weight);
      }
    }

    if (calchib > CACHIBWP1) {
      h->h1D("calsubjetNChibWP1", "", s)->Fill(casjn, weight);
      if (std::fabs(calj.mom().Eta()) < 0.7) {
        h->h1D("calsubjetNLEChibWP1", "", s)->Fill(casjn, weight);
      } else {
        h->h1D("calsubjetNHEChibWP1", "", s)->Fill(casjn, weight);
      }
      if (casjn >= 1)
        h->h1D("calsubjet0PtChibWP1", "", s)->Fill(cavsj[0].v.Perp()*1e-3, weight);
      if (casjn >= 2)
        h->h1D("calsubjet1PtChibWP1", "", s)->Fill(cavsj[1].v.Perp()*1e-3, weight);
      if (casjn >= 3)
        h->h1D("calsubjet2PtChibWP1", "", s)->Fill(cavsj[2].v.Perp()*1e-3, weight);
      if (casjn >= 4)
        h->h1D("calsubjet3PtChibWP1", "", s)->Fill(cavsj[3].v.Perp()*1e-3, weight);
      if (casjn >= 5)
        h->h1D("calsubjet4PtChibWP1", "", s)->Fill(cavsj[4].v.Perp()*1e-3, weight);
      if (casjn >= 1)
        h->h1D("calsubjet0EtaChibWP1", "", s)->Fill(cavsj[0].v.Eta(), weight);
      if (casjn >= 2)
        h->h1D("calsubjet1EtaChibWP1", "", s)->Fill(cavsj[1].v.Eta(), weight);
      if (casjn >= 3)
        h->h1D("calsubjet2EtaChibWP1", "", s)->Fill(cavsj[2].v.Eta(), weight);
      if (casjn >= 4)
        h->h1D("calsubjet3EtaChibWP1", "", s)->Fill(cavsj[3].v.Eta(), weight);
      if (casjn >= 5)
        h->h1D("calsubjet4EtaChibWP1", "", s)->Fill(cavsj[4].v.Eta(), weight);
    }

    if (calchib > CACHIBWP2) {
      h->h1D("calsubjetNChibWP2", "", s)->Fill(casjn, weight);
      if (std::fabs(calj.mom().Eta()) < 0.7) {
        h->h1D("calsubjetNLEChibWP2", "", s)->Fill(casjn, weight);
      } else {
        h->h1D("calsubjetNHEChibWP2", "", s)->Fill(casjn, weight);
      }
      if (casjn >= 1)
        h->h1D("calsubjet0PtChibWP2", "", s)->Fill(cavsj[0].v.Perp()*1e-3, weight);
      if (casjn >= 2)
        h->h1D("calsubjet1PtChibWP2", "", s)->Fill(cavsj[1].v.Perp()*1e-3, weight);
      if (casjn >= 3)
        h->h1D("calsubjet2PtChibWP2", "", s)->Fill(cavsj[2].v.Perp()*1e-3, weight);
      if (casjn >= 4)
        h->h1D("calsubjet3PtChibWP2", "", s)->Fill(cavsj[3].v.Perp()*1e-3, weight);
      if (casjn >= 5)
        h->h1D("calsubjet4PtChibWP2", "", s)->Fill(cavsj[4].v.Perp()*1e-3, weight);
      if (casjn >= 1)
        h->h1D("calsubjet0EtaChibWP2", "", s)->Fill(cavsj[0].v.Eta(), weight);
      if (casjn >= 2)
        h->h1D("calsubjet1EtaChibWP2", "", s)->Fill(cavsj[1].v.Eta(), weight);
      if (casjn >= 3)
        h->h1D("calsubjet2EtaChibWP2", "", s)->Fill(cavsj[2].v.Eta(), weight);
      if (casjn >= 4)
        h->h1D("calsubjet3EtaChibWP2", "", s)->Fill(cavsj[3].v.Eta(), weight);
      if (casjn >= 5)
        h->h1D("calsubjet4EtaChibWP2", "", s)->Fill(cavsj[4].v.Eta(), weight);
    }

    if (calchib > CACHIBWP3) {
      h->h1D("calsubjetNChibWP3", "", s)->Fill(casjn, weight);
      if (std::fabs(calj.mom().Eta()) < 0.7) {
        h->h1D("calsubjetNLEChibWP3", "", s)->Fill(casjn, weight);
      } else {
        h->h1D("calsubjetNHEChibWP3", "", s)->Fill(casjn, weight);
      }
      if (casjn >= 1)
        h->h1D("calsubjet0PtChibWP3", "", s)->Fill(cavsj[0].v.Perp()*1e-3, weight);
      if (casjn >= 2)
        h->h1D("calsubjet1PtChibWP3", "", s)->Fill(cavsj[1].v.Perp()*1e-3, weight);
      if (casjn >= 3)
        h->h1D("calsubjet2PtChibWP3", "", s)->Fill(cavsj[2].v.Perp()*1e-3, weight);
      if (casjn >= 4)
        h->h1D("calsubjet3PtChibWP3", "", s)->Fill(cavsj[3].v.Perp()*1e-3, weight);
      if (casjn >= 5)
        h->h1D("calsubjet4PtChibWP3", "", s)->Fill(cavsj[4].v.Perp()*1e-3, weight);
      if (casjn >= 1)
        h->h1D("calsubjet0EtaChibWP3", "", s)->Fill(cavsj[0].v.Eta(), weight);
      if (casjn >= 2)
        h->h1D("calsubjet1EtaChibWP3", "", s)->Fill(cavsj[1].v.Eta(), weight);
      if (casjn >= 3)
        h->h1D("calsubjet2EtaChibWP3", "", s)->Fill(cavsj[2].v.Eta(), weight);
      if (casjn >= 4)
        h->h1D("calsubjet3EtaChibWP3", "", s)->Fill(cavsj[3].v.Eta(), weight);
      if (casjn >= 5)
        h->h1D("calsubjet4EtaChibWP3", "", s)->Fill(cavsj[4].v.Eta(), weight);
    }
  }
}

