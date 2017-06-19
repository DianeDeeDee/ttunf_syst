#include "Plot.h"
#include "PlotSemilep.h"
#include "Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "HistogramService.h"
#include <iostream>
using namespace std;
bool sortFunction(const LargeJet &a, const LargeJet &b){
  return a.mom().M() > b.mom().M();
}

PlotSemilep::PlotSemilep(const std::string &filename, bool electron, const std::vector<std::string> &systs)
: Plot(filename, systs), m_electron(electron) {
    m_hSvc.create1D("lepPt", "; Lepton p_{T} (GeV) ; Events", 40, 0., 1000.);
    m_hSvc.create1D("MET", "; Missing ET (GeV)", 40, 0., 1000.);
    m_hSvc.create1D("jetPt", "; Jet p_{T} (GeV) ; Events", 40, 0., 1000.);
    m_hSvc.create1D("NLargeJets", "; Number of Fat jets ; Events", 8, 0, 8);
    m_hSvc.create1D("NSmallJets", "; Number of Small jets ; Events", 12,0,12);
    m_hSvc.create1D("NbJets", "; Number of b-jets ; Events", 8,0,8);
    m_hSvc.create1D("seljetPt", "; Selected Jet p_{T} (GeV) ; Events", 40, 0., 1000.);

    m_hSvc.create1D("lepPtt", "; Lepton p_{T} (GeV) ; Events", 40, 0., 1000.);
    m_hSvc.create1D("METt", "; Missing ET (GeV)", 40, 0., 1000.);
    m_hSvc.create1D("jetPtt", "; Jet p_{T} (GeV) ; Events", 40, 0., 1000.);
    m_hSvc.create1D("NLargeJetst", "; Number of Fat jets ; Events", 8, 0, 8);
    m_hSvc.create1D("NSmallJetst", "; Number of Small jets ; Events", 12,0,12);
    m_hSvc.create1D("NbJetst", "; Number of b-jets ; Events", 8,0,8);
    m_hSvc.create1D("seljetPtt", "; Selected Jet p_{T} (GeV) ; Events", 40, 0., 1000.);

    m_hSvc.create1D("dR", "; dR(largejet, sel_jet) ; Events", 14, 1., 8.);
    m_hSvc.create1D("dPhi", "; dPhi(largejet, lepton) ; Events", 6, 1., 4.);
    m_hSvc.create1D("M", "; Large jet mass (GeV) ; Events", 25, 0., 500.);
    m_hSvc.create1D("Pt", "; Large Jet p_{T} (GeV) ; Events", 30, 100., 1000.);
    m_hSvc.create1D("D12", "; D12 (GeV) ; Events", 10, 0., 400.);
    m_hSvc.create1D("dR_seljEl", "; dR(Seljet, lepton) ; Events", 25, 0., 5.);	    
    m_hSvc.create1D("dR_bjEl", "; dR(bjet, lepton) ; Events", 25, 0., 5.);
    m_hSvc.create1D("Eta_Sj_0", "; Eta of the first Smalljet ; Events", 16, -4., 4.); 
    m_hSvc.create1D("dR_SjEl", "; dR(Smalljet, lepton) ; Events", 25, 0., 5.);
    m_hSvc.create1D("Pt_Sj_0", "; first samll Jet p_{T} (GeV) ; Events", 20, 0., 1000.);
    m_hSvc.create1D("Pt_Sj", "; small Jet.pass p_{T} (GeV) ; Events", 20, 0., 1000.);
    m_hSvc.create1D("lepEta", "; lepton Eta ; Events", 120, -3., 3.);
    m_hSvc.create1D("seljetE", "; Selected Jet E_{T} (GeV) ; Events", 21, -100., 2000.);
    m_hSvc.create1D("E_Sj", "; Small Jet.pass E_{T} (GeV) ; Events", 21, -100., 2000.);
    m_hSvc.create1D("E_Sj_0", "; first Small Jet.pass E_{T} (GeV) ; Events", 21, -100., 2000.);
}

PlotSemilep::~PlotSemilep() {
}

void PlotSemilep::run(const Event &e, double weight, double pweight, const std::string &s) {
    HistogramService *h = &m_hSvc;
    double Mmj = 0.;
    double Mej = 0.;
    double Mmjnu = 0.;
    double Mejnu = 0.;
    double MWj=0.;
    // double  MWb = 0;
    int goodjet = 0;
    int sel_jet = -1;
    double dr = 1.5;
    double lpt =0.;
    int fatjets = 0;
    int btags =0;
    int jetsmall = 0;
    int idx = -1;
    //For MET
    float METx = e.met().Px();//[Missing E_T x direction]; // in MeV
    float METy = e.met().Py();//[Missing E_T y direction]; // in MeV
    float mass = 0; // mas of the lepton
    double Menu = 0;
    double Mmnu = 0;
    
    // apply extra cuts
    
    if (e.passReco()) {// events passes cuts from EventCutterReco
        TLorentzVector l;
        TLorentzVector momNu;
        // this will be returned by the function -> it is the neutrino vector including the p_z
        TLorentzVector momW; // this will be returned by the function -> it is the W boson vector after the neutrino p_z is found
        TLorentzVector momLepton; float NuX; float NuY;float NuZ; float lepMass;
        if (m_electron) {
            l = e.electron()[0].mom();
            mass = 0.510998910; // mass of the electron in MeV
        } else {
            l = e.muon()[0].mom();
            mass = 105.6583668; // mass of the muon in MeV
        }
        h->h1D("lepPtt", "", s)->Fill(l.Perp()*1e-3, weight);
         h->h1D("lepEta", "", s)->Fill(l.Eta(), weight);
        getWFromLeptonicDecay(l, METx, METy, momNu, momW, mass);
        std::vector<LargeJet> sortedJets = e.largeJet(); //to sort the jets by mass
        std::sort(sortedJets.begin(), sortedJets.end(), sortFunction);
       
        bool passIsBadLoose = true; 
        double mindr0 = 9999;
        for (int k = 0; k < e.jet().size(); ++k) {
	    if (!e.jet()[k].passBadLooseMinus()) passIsBadLoose = false;
            double tdr = e.jet()[k].mom().DeltaR(l);
            if (e.jet()[k].pass()) {
            	    jetsmall++;
               	    if (tdr < mindr0) mindr0 = tdr;
                    h->h1D("Pt_Sj", "", s)->Fill(e.jet()[k].mom().Perp()*1e-3, weight);
                    h->h1D("E_Sj", "", s)->Fill(e.jet()[k].mom().E()*1e-3, weight);
               
            }
	  }
        h->h1D("jetPtt", "", s)->Fill(e.jet()[0].mom().Perp()*1e-3, weight);
        //if (!passIsBadLoose) return;
        h->h1D("dR_SjEl", "", s)->Fill(mindr0, weight);
        h->h1D("NSmallJetst", "", s)->Fill(jetsmall, weight);
        h->h1D("Eta_Sj_0", "", s)->Fill(e.jet()[0].mom().Eta(), weight);
        h->h1D("Pt_Sj_0", "", s)->Fill(e.jet()[0].mom().Perp()*1e-3, weight);
        //h->h1D("Pt_Sj_1", "", s)->Fill(e.jet()[1].mom().Perp()*1e-3, weight);
        h->h1D("E_Sj_0", "", s)->Fill(e.jet()[0].mom().E()*1e-3, weight);
        //h->h1D("E_Sj_0", "", s)->Fill(e.jet()[0].mom().E()*1e-3, weight);
	double mindr = 9999;
        for (int k = 0; k < e.jet().size(); ++k) { 
	   double tdr = e.jet()[k].mom().DeltaR(l);
            if (e.jet()[k].pass()) {
	   // if (tdr < dr && e.jet()[k].mom().Perp() > lpt && e.jet()[k].btag()) {
            if (tdr < dr && e.jet()[k].mom().Perp() > lpt) {
                sel_jet = k;
                lpt = e.jet()[k].mom().Perp();
                if (tdr < mindr) mindr = tdr;
            }
          }
        }
        if (sel_jet == -1) return;
	const TLorentzVector &selj = e.jet()[sel_jet].mom();
        h->h1D("seljetPtt", "", s)->Fill(selj.Perp()*1e-3, weight);
        h->h1D("seljetE", "", s)->Fill(selj.E()*1e-3, weight);
        h->h1D("dR_seljEl", "", s)->Fill(mindr, weight);
 
	if(m_electron){
            Mej = (e.electron()[0].mom()+e.jet()[sel_jet].mom()).M();
        } else{
            Mmj = (e.muon()[0].mom() + e.jet()[sel_jet].mom()).M();
        }
        MWj = momW.M()+e.jet()[sel_jet].mom().M();
        //h->h1D("lepPt", "", s)->Fill(l.Perp()*1e-3, weight);
        h->h1D("METt","",s)->Fill(momNu.Perp()*1e-3, weight);
        //const TLorentzVector &selj = e.jet()[sel_jet].mom();
        const TLorentzVector &j = e.jet()[0].mom();
        
        float mtw = std::sqrt(2*l.Perp()*momNu.Perp()*(1 - std::cos(l.Phi() - momNu.Phi())) );
  // std::cout<<"momNu.Perp()= "<<momNu.Perp()<<"   momNu.Phi()="<<momNu.Phi()<<"     cos(l.Phi() - e.momNu.Phi())="<<std::cos(l.Phi() - e.momNu.Phi())<<std::endl;
  //    //std::cout<<"l.Perp()="<<l.Perp()<<"     l.Phi()="<<l.Phi()<<std::endl;
        if((momNu.Perp()<20e3)) std::cout<<"MET<20e3!!!  MET="<<(momNu.Perp())/1000.<<std::endl;
            if((mtw+momNu.Perp())<60e3) std::cout<<"mtw+MET<60e3!!!  mtw+MET="<<(momNu.Perp()+mtw)/1000.<<std::endl;

        for (int k = 0; k < e.largeJet().size(); k++) {
            bool passdPhi=false;
            bool passdR = false;
            bool passM = false;
            bool passD12 = false;
            bool passPt = false;
            if (e.largeJet()[k].pass()) {
                passdR = e.jet()[sel_jet].mom().DeltaR(e.largeJet()[k].mom()) > 1.5;
                if(m_electron)passdPhi=e.largeJet()[k].mom().DeltaPhi(e.electron()[0].mom()) > 2.3;
                else passdPhi=e.largeJet()[k].mom().DeltaPhi(e.muon()[0].mom()) > 2.3;
                passM = e.largeJet()[k].mom().M() > 100e3;
                passPt = e.largeJet()[k].mom().Perp() > 300e3;//350e3;130e3; //if true, it lower electron yield
                passD12 = e.largeJet()[k].split12() > 40e3;
                if (passdR && passdPhi && passM && passD12 &&passPt) {
                    fatjets++;
                    idx = k;
                    h->h1D("dR", "", s)->Fill(e.jet()[sel_jet].mom().DeltaR(e.largeJet()[k].mom()), weight);
                    h->h1D("dPhi", "", s)->Fill(e.largeJet()[k].mom().DeltaPhi(l), weight);
                    h->h1D("M", "", s)->Fill(e.largeJet()[k].mom().M()*1e-3, weight);
                    h->h1D("Pt", "", s)->Fill(e.largeJet()[k].mom().Perp()*1e-3, weight);
                    h->h1D("D12", "", s)->Fill(e.largeJet()[k].split12()*1e-3, weight);
            	   }
		}
        }
         if (jetsmall <5  && fatjets < 1) return;
        //if (fatjets < 1) return;
        h->h1D("NLargeJetst", "", s)->Fill(fatjets, weight);
       const TLorentzVector &lj = e.largeJet()[idx].mom(); 
       /* const TLorentzVector &selj = e.jet()[sel_jet].mom();
        
        const TLorentzVector &j = e.jet()[0].mom();
        const TLorentzVector &lj = e.largeJet()[idx].mom();
        */
      
     /*   h->h1D("NbJets", "", s)->Fill(btags, weight);
        h->h1D("lepPt", "", s)->Fill(l.Perp()*1e-3, weight);
        h->h1D("jetPt", "", s)->Fill(j.Perp()*1e-3, weight);*/
      //  h->h1D("seljetPt", "", s)->Fill(selj.Perp()*1e-3, weight);
        double mindr2 = 9999;
        for (int k = 0; k < e.jet().size(); ++k) {
            double tdr = e.jet()[k].mom().DeltaR(l);
            if (e.jet()[k].pass() && e.jet()[k].btag()){//&& (e.jet()[k].trueFlavour() == 5)) {
                btags++;
                if (tdr < mindr2) mindr2 = tdr;
            }
        }
        if (btags < 1) return;
        h->h1D("dR_bjEl", "", s)->Fill(mindr2, weight);
        h->h1D("NSmallJets", "", s)->Fill(jetsmall, weight);
        h->h1D("NLargeJets", "", s)->Fill(fatjets, weight);
        h->h1D("NbJets", "", s)->Fill(btags, weight);
        h->h1D("seljetPt", "", s)->Fill(selj.Perp()*1e-3, weight);
        h->h1D("MET","",s)->Fill(momNu.Perp()*1e-3, weight);
        h->h1D("jetPt", "", s)->Fill(j.Perp()*1e-3, weight);
        h->h1D("lepPt", "", s)->Fill(l.Perp()*1e-3, weight);
    }//e.passReco loop
    
}//Run function

bool PlotSemilep::getWFromLeptonicDecay(TLorentzVector &momLepton, float NuX, float NuY, TLorentzVector &momNu, TLorentzVector &momW, float lepMass) {
    float massW = 80.385*1e+3;//  float A =0.;
    float A = ((massW*massW - lepMass*lepMass)*0.5 + momLepton.Px()*NuX + momLepton.Py()*NuY)/momLepton.E();
    float B = momLepton.Pz()/momLepton.E();
    float discriminator = A*A*B*B - (B*B - 1)*(A*A - NuX*NuX - NuY*NuY);
    // We need sqrt(discriminator) ... if discriminator < 0, reject this lepton
    
    if (discriminator < 0) {
        // no solution
        // let's keep the real part only
        float NuZ = -A*B/(B*B - 1);
        momNu.SetPxPyPzE(NuX, NuY, NuZ, sqrt(std::pow(NuX, 2) + std::pow(NuY, 2) + std::pow(NuZ, 2)));
        momW = momLepton + momNu;
        float dm = fabs(momW.M() - massW);
    } else if (discriminator == 0) {
        float NuZ = 0;
        NuZ = -A*B/(B*B - 1);
        momNu.SetPxPyPzE(NuX, NuY, NuZ, sqrt(std::pow(NuX, 2) + std::pow(NuY, 2) + std::pow(NuZ, 2)));
        momW = momLepton + momNu;
    } else { // two solutions
        float NuZ1 = 0;
        float NuZ2 = 0;
        NuZ1 = (-A*B + std::sqrt(discriminator))/(B*B - 1);
        NuZ2 = (-A*B - std::sqrt(discriminator))/(B*B - 1);
        TLorentzVector momNu1;
        TLorentzVector momNu2;
        momNu1.SetPxPyPzE(NuX, NuY, NuZ1, sqrt(std::pow(NuX, 2) + std::pow(NuY, 2) + std::pow(NuZ1, 2)));
        momNu2.SetPxPyPzE(NuX, NuY, NuZ2, sqrt(std::pow(NuX, 2) + std::pow(NuY, 2) + std::pow(NuZ2, 2)));
        TLorentzVector momW1;
        TLorentzVector momW2;
        momW1 = momLepton + momNu1;
        momW2 = momLepton + momNu2;
        
        
        if (fabs(NuZ2) < fabs(NuZ1)) {
            momNu = momNu2;
            momW = momW2;
        } else {
            momNu = momNu1;
            momW = momW1;
        }
    }
    
    return true;
    
}  


/*TLorentzVector j1;
 TLorentzVector j2;
 int jets_that_pass = 0;
 for (int k = 0; k < e.jets().size(); ++k) {
 if (e.jets()[k].pass()) {
 if (jets_that_pass == 0) j1 = e.jets()[k].mom();
 if (jets_that_pass == 1) j2 = e.jets()[k].mom();
 jets_that_pass++;
 }
 }*/

/* m_hSvc.create2D("subjetN_pr_mig", "; Sub-jet multiplicity (part+reco); Events", 6, 0.5, 6.5, 6, 0.5, 6.5);
 m_hSvc.create1D("subjetN_pr", "; Sub-jet multiplicity (part+reco); Events", 6, 0.5, 6.5);
 m_hSvc.create1D("subjetN_npr", "; Sub-jet multiplicity (!part+reco); Events", 6, 0.5, 6.5);
 m_hSvc.create1D("p_largeJetPt", "; Large jet p_{T}(GeV) ; Events", 100, 0, 1000);
 m_hSvc.create1D("p_largeJetM", "; Large jet M (GeV); Events", 100, 0, 1000);
 m_hSvc.create1D("p_subjetPt", "; Sub-jet p_{T}(GeV) ; Events", 100, 0, 1000);
 m_hSvc.create1D("p_subjetN", "; Sub-jet multiplicity ; Events", 6, 0.5, 6.5);*/

