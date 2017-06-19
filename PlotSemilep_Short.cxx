#include "Plot.h"
#include "PlotSemilep.h"
#include "Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "HistogramService.h"
#include <iostream>
#include "LargeJet.h"
#include "Jet.h"
#include <algorithm>
using namespace std;
#include <math.h>

bool sortFunction(const LargeJet &a, const LargeJet &b){
  return a.mom().M() > b.mom().M();
}
PlotSemilep::PlotSemilep(const std::string &filename, bool electron, const std::vector<std::string> &systs)
  : Plot(filename, systs), m_electron(electron) {
  m_hSvc.create1D("lepPt", "; Lepton p_{T} ; Events", 100, 0, 1000);
  m_hSvc.create1D("jetPt", "; Jet p_{T} ; Events", 100, 0, 1000);
  m_hSvc.create1D("largeJetPt", "; Large jet p_{T} ; Events", 100, 0, 1000);
  m_hSvc.create1D("largeJetM", "; Large jet M ; Events", 100, 0, 1000);

}

PlotSemilep::~PlotSemilep() {
}

void PlotSemilep::run(const Event &e, double weight, double pweight, const std::string &s) {
  HistogramService *h = &m_hSvc;
  float mass = 0;
  /*if ((e.channelNumber() == 117050 || e.channelNumber() == 117001) && (m_ttbarPtWeight != 0) && e.partMom().size() != 0) {
    double tw = topPtWeight(e);
    weight *= tw;
    pweight *= tw;
  }*/
  if (e.passReco()) {
     TLorentzVector l;
std::cout<<"Plot: s= "<<s<<"   pweight= "<<pweight<<"  weight= "<<weight<<" m_electron= "<<m_electron<<std::endl;
    if (m_electron) {
std::cout<<"Plot2: s= "<<s<<"   pweight= "<<pweight<<"  weight= "<<weight<<" m_electron= "<<m_electron<<std::endl;
      l = e.electron()[0].mom();
std::cout<<"Plot2p: s= "<<s<<"   pweight= "<<pweight<<"  weight= "<<weight<<" m_electron= "<<m_electron<<std::endl;

    } else {
      l = e.muon()[0].mom();
std::cout<<"Plot3: s= "<<s<<"   pweight= "<<pweight<<"  weight= "<<weight<<" m_electron= "<<m_electron<<std::endl;

    }	
std::cout<<"Plot3: s= "<<s<<"   pweight= "<<pweight<<"  weight= "<<weight<<" m_electron= "<<m_electron<<std::endl;

    /*TLorentzVector l;
      if (m_electron) {
      l = e.electron()[0].mom();
      mass = 0.510998910; // mass of the electron in MeV
         } else {
      l = e.muon()[0].mom();
      mass = 105.6583668; // mass of the muon in MeV
    } */
    
    h->h1D("lepPt", "", s)->Fill(l.Perp()*1e-3, weight);
    
    const TLorentzVector &j = e.jet()[0].mom();
    h->h1D("jetPt", "", s)->Fill(j.Perp()*1e-3, weight);

    const TLorentzVector &lj = e.largeJet()[0].mom();
    h->h1D("largeJetPt", "", s)->Fill(lj.Perp()*1e-3, weight);
    h->h1D("largeJetM", "", s)->Fill(lj.M()*1e-3, weight);

  
  }

}

