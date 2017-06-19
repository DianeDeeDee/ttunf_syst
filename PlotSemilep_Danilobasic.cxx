#include "Plot.h"
#include "PlotSemilep.h"
#include "Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "HistogramService.h"

bool sortFunction(const LargeJet &a, const LargeJet &b){
  return a.mom().M() > b.mom().M();
}

//PlotSemilep::PlotSemilep(const std::string &filename, bool electron, const std::vector<std::string> &systs)
  //: Plot(filename, systs), m_electron(electron) {
 PlotSemilep::PlotSemilep(const std::string &filename, bool electron, const std::vector<std::string> &systs)
: Plot(filename, systs), m_electron(electron) {
  m_hSvc.create1D("lepPt", "; Lepton p_{T} ; Events", 100, 0, 1000);
  m_hSvc.create1D("jetPt", "; Jet p_{T} ; Events", 100, 0, 1000);
  m_hSvc.create1D("largeJetPt", "; Large jet p_{T} ; Events", 100, 0, 1000);
  m_hSvc.create1D("largeJetM", "; Large jet M ; Events", 100, 0, 1000);
  m_hSvc.create1D("NbJets", "; Number of b-jets ; Events", 8,0,8);
  m_hSvc.create1D("NLargeJets", "; Number of Fat jets ; Events", 8, 0, 8);
  m_hSvc.create1D("NSmallJets", "; Number of Small jets ; Events", 12,0,12);
}

PlotSemilep::~PlotSemilep() {
}

void PlotSemilep::run(const Event &e, double weight, double pweight, const std::string &s) {
  HistogramService *h = &m_hSvc;
  if (e.passReco()) {
    TLorentzVector l;
    if (m_electron) {
      l = e.electron()[0].mom();
    } else {
      l = e.muon()[0].mom();
    }

 std::vector<LargeJet> sortedJets = e.largeJet(); //to sort the jets by mass
    std::sort(sortedJets.begin(), sortedJets.end(), sortFunction);

    h->h1D("lepPt", "", s)->Fill(l.Perp()*1e-3, weight);

    const TLorentzVector &j = e.jet()[0].mom();
    h->h1D("jetPt", "", s)->Fill(j.Perp()*1e-3, weight);
    h->h1D("NSmallJets", "", s)->Fill(e.jet().size(), weight);
    if(e.largeJet().size()>0) {
     const TLorentzVector &lj = e.largeJet()[0].mom();
     h->h1D("largeJetPt", "", s)->Fill(lj.Perp()*1e-3, weight);
     }
     h->h1D("NLargeJets", "", s)->Fill(e.largeJet().size(), weight); 
     
}
}

