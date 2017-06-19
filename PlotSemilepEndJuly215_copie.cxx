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
//bool sortFunction2(const Jet &a, const Jet &b){
//  return a.mom().Perp() > b.mom().Perp();
//}
PlotSemilep::PlotSemilep(const std::string &filename, bool electron, const std::vector<std::string> &systs)
  : Plot(filename, systs), m_electron(electron) {
 
  m_hSvc.create1D("MET", " Missing ET (GeV)",81, -10, 800); //80, 0 800

 m_hSvc.create1D("seljetPt", " Selected Jet p_{T} (GeV) ; Selected Jet p_{T} (GeV)", 10, 0, 1400); //old: 10, 0, 1400
  m_hSvc.create1D("jetPt", " Jet p_{T} (GeV) ; Jet p_{T} (GeV)", 30, 0, 3000); //30 0 3000
   m_hSvc.create1D("lepPt", " Lepton p_{T} (GeV) ; Lepton Pt (GeV)", 30, 0, 3000); //old: 30 0 3000
 m_hSvc.create1D("lepPt2", " Lepton p_{T} (GeV) when dR < 0.4; Pt(dR<0.4) (GeV)", 100, 0, 500); //old: 30 0 3000
m_hSvc.create1D("Eta", "Eta(Lepton)", 30, -3., 3); //done
	m_hSvc.create1D("NbJets", "; Number of b-jets ; N_{bjets}", 20,0,20);
//	m_hSvc.create1D("NbJets1", "; Number of b-jets1 ; Events", 13,-0.5,12.5);
//	m_hSvc.create1D("NbJets2", "; Number of b-jets2 ; Events", 13,-0.5,12.5);
//	m_hSvc.create1D("NbJets3", "; Number of b-jets3 ; Events", 13,-0.5,12.5);
//	m_hSvc.create1D("NbJets4", "; Number of b-jets4 ; Events", 13,-0.5,12.5);
//	m_hSvc.create1D("NbJets5", "; Number of b-jets5 ; Events", 13,-0.5,12.5);
	
	m_hSvc.create1D("NSmallJets0", "; Size of Small jets ; Events", 20,0,20);
	m_hSvc.create1D("NSmallJets", "; Number of Small jets ; N_{Small jets}", 20,0,20);  
//	m_hSvc.create1D("NSmallJets1", "; Number of Small jets1 ; Events", 21,-0.5,20.5);
//	m_hSvc.create1D("NSmallJets3", "; Number of Small jets3 ; Events", 21,-0.5,20.5);
//	m_hSvc.create1D("NSmallJets4", "; Number of Small jets4 ; Events", 21,-0.5,20.5);
//	m_hSvc.create1D("NSmallJets5", "; Number of Small jets5 ; Events", 21,-0.5,20.5);
//	
	m_hSvc.create1D("NSelJets", "; Number of Selected jets ; N_{Sel. jet}", 20,0,20);
//	m_hSvc.create1D("NSelJets1", "; Number of Selected jets1 ; Events", 31,-1,30);
//	m_hSvc.create1D("NSelJets2", "; Number of Selected jets2 ; Events", 16,-0.5,15.5);
//	m_hSvc.create1D("NSelJets3", "; Number of Selected jets3 ; Events", 16,-0.5,15.5);
//	m_hSvc.create1D("NSelJets4", "; Number of Selected jets4 ; Events", 16,-0.5,15.5);
//	m_hSvc.create1D("NSelJets5", "; Number of Selected jets5 ; Events", 16,-0.5,15.5);
	
	
	m_hSvc.create1D("NSmaljetsFarLep", "; Number of Small jets far from the lepton ; Events", 20,0,20);   
//	m_hSvc.create1D("NSmaljetsFarLep3", "; Number of SmaljetsFarLep3 ; Events", 16,-0.5,15.5);
//	m_hSvc.create1D("NSmaljetsFarLep4", "; Number of SmaljetsFarLep4 ; Events", 16,-0.5,15.5);
//	m_hSvc.create1D("NSmaljetsFarLep5", "; Number of SmaljetsFarLep5 ; Events", 16,-0.5,15.5);
	
	m_hSvc.create1D("NSmaljetsCloseLep", "; Number of Small jets close to the lepton ; Events", 20,0,20);
//	m_hSvc.create1D("NSmaljetsCloseLep3", "; Number of SmaljetsCloseLep3 ; Events", 16,-0.5,15.5);
//	m_hSvc.create1D("NSmaljetsCloseLep4", "; Number of SmaljetsCloseLep4 ; Events", 16,-0.5,15.5);
//	m_hSvc.create1D("NSmaljetsCloseLep5", "; Number of SmaljetsCloseLep5 ; Events", 16,-0.5,15.5);
//	
//	m_hSvc.create1D("Nfatjets0", "; Size of Fat jets ; Events", 20, 0, 20);
	m_hSvc.create1D("NLargeJets", "; Number of Fat jets ; N_{Fat jets}", 20, 0, 20);
//	m_hSvc.create1D("Nfatjets1", "; Number of Fat jets1 ; Events", 7, -0.5, 6.5);
//	m_hSvc.create1D("Nfatjets2", "; Number of Fat jets2 ; Events", 7, -0.5, 6.5);
//	m_hSvc.create1D("Nfatjets3", "; Number of Fat jets3 ; Events", 7, -0.5, 6.5);
//	m_hSvc.create1D("Nfatjets4", "; Number of Fat jets4 ; Events", 7, -0.5, 6.5);
//	m_hSvc.create1D("Nfatjets5", "; Number of Fat jets5 ; Events", 7, -0.5, 6.5);
//	
//	m_hSvc.create1D("NSmallJets_chan_FatJet0", "; Number of small jets -0 Fat Jets ; Events", 21, -1, 20);
//	m_hSvc.create1D("SmaljetsCloseLep_chan_FatJet0", "; Number of small jets close to lept-0 Fat Jets ; Events", 21, -1, 20);
//	m_hSvc.create1D("SmaljetsFarLep_chan_FatJet0", "; Number of small jets far from lept-0 Fat Jets ; Events", 21, -1, 20);
//	
	m_hSvc.create1D("NSmallJets_chan_FatJet1", "; Number of small jets -1 Fat Jets ; Events", 21, -1, 20);
	m_hSvc.create1D("SmaljetsCloseLep_chan_FatJet1", "; Number of small jets close to lept-1 Fat Jets ; Events", 21, -1, 20);
	m_hSvc.create1D("SmaljetsFarLep_chan_FatJet1", "; Number of small jets far from lept-1 Fat Jets ; Events", 21, -1, 20);
	
	m_hSvc.create1D("NSmallJets_chan_FatJet2", "; Number of small jets -2 Fat Jets ; Events", 21, -1, 20);
	m_hSvc.create1D("SmaljetsCloseLep_chan_FatJet2", "; Number of small jets close to lept-2 Fat Jets ; Events", 21, -1, 20);
	m_hSvc.create1D("SmaljetsFarLep_chan_FatJet2", "; Number of small jets far from lept-2 Fat Jets ; Events", 21, -1, 20);
	
	m_hSvc.create1D("MET_plus_mtW", "; MET+mtW ; Events", 50, 0, 100);
	m_hSvc.create1D("mtW", "; mtW ; Events", 50, 0, 100);

//1 Fat jet
 m_hSvc.create1D("MassGStar1FatJet2jets_Fatjet2jets", " M(Jet_Fat[0]+SjetsCloseLep[0]+SjetsCloseLep[1]) (GeV)", 360, 0, 3600);
 m_hSvc.create1D("MassZprimeFatJet_FatJet0", " ZPrime_1FatJet: M(FatJet[0]) (GeV)", 360, 0, 3600);
 m_hSvc.create1D( "Pt_1FatJet_2jets_FatJet0", " 2Sjets_1FatJet: FatJet[0] p_{T} (GeV); Events", 50, 0, 2000);//done old: 20, 0, 2000
 m_hSvc.create1D( "Pt_1FatJet_2jets_jetsCloseLep0", " 1FatJet_2jets: JetCloseLep[0] p_{T} (GeV); Events", 35, 0, 1400);//done old: 20, 0, 2000
 m_hSvc.create1D( "Pt_1FatJet_2jets_jetsCloseLep1", " 1FatJet_2jets: JetCloseLep[1] p_{T} (GeV); Events", 40, 0, 800);//done old: 20, 0, 2000
 m_hSvc.create1D("Mass_1FatJet_2jets_FatJet0", " 1FatJet_2jets: M(FatJet[0]) (GeV)", 50, 0, 1000); //done old: 36, 0, 3600
 m_hSvc.create1D("Eta_1FatJet_2jets_FatJet0", " 1FatJet_2jets: Eta(FatJet[0])", 30, -3., 3); //done
 m_hSvc.create1D("Eta_1FatJet_2jets_jetsCloseLep0", " 1FatJet_2jets: Eta(JetCloseLep[0])", 30, -3., 3); //done?
 m_hSvc.create1D("Eta_1FatJet_2jets_jetsCloseLep1", " 1FatJet_2jets: Eta(JetCloseLep[1])", 30, -3., 3); //done?
 m_hSvc.create1D( "PtZprimeFatJet_W", " 1FatJet: W p_{T} (GeV); Events", 200, 0, 2000);  
 m_hSvc.create1D( "Pt1FatJet2jets_2jets", " 2jets_1FatJet: Jet[1+2] p_{T} (GeV); Events", 200, 0, 2000);
 m_hSvc.create1D("Pt1FatJet_2jets", " 1FatJet_2jetsCloseToLep p_{T} (GeV) ; Events", 200, 0, 2000);
  // m_hSvc.create1D("MassZprime_FatJet", " Z', T, gkk mass with a fat Jet (GeV)", 360, 0, 3600);


//2 fat jets
 m_hSvc.create1D("MassGStar2FatJet_2FatJets", " M(Jet_Fat[0]+Jet_Fat[1]) (GeV)", 360, 0, 3600);
 m_hSvc.create1D("Pt2FatJet_FatJet1", " FatJet_1 p_{T} (GeV) ; Events", 60, 0, 1200); //done
 m_hSvc.create1D("Eta2FatJet_FatJet0", " 2FatJets: Eta(FatJet[0])", 30, -3., 3); //done
 m_hSvc.create1D("Eta2FatJet_FatJet1", " 2FatJets: Eta(FatJet[1])", 30, -3., 3); //done
 m_hSvc.create1D("Pt2FatJet", " FatJets p_{T} (GeV) ; Events", 200, 0, 2000);
 m_hSvc.create1D("Mass2FatJet_FatJet0", " First Large jet M (GeV); Events", 40, 0, 800); //done old: 100, 0, 1000
 m_hSvc.create1D("PtGStar_2FatJet", " GStar_2FatJet p_{T} (GeV) ; Events", 200, 0, 2000);
 m_hSvc.create1D("Mass2FatJet_FatJet1", " Second Large jet M (GeV); Events", 80, 0, 800); //old: 100, 0, 1000
 m_hSvc.create1D("MassGStar_2FatJet", " GStar Mass (with 2 fat jets + ttbar) (GeV)", 360, 0, 3600);
 m_hSvc.create1D("Pt2FatJet_FatJet0", " FatJet_0 p_{T} (GeV) ; Events", 35, 0, 1400); //done
 m_hSvc.create1D("dR", " dR (lepton,closestJet) ; dR (lepton,closestJet)", 35, 0., 3.5); //done 
              	
}

PlotSemilep::~PlotSemilep() {
}

void PlotSemilep::run(const Event &e, double weight, double pweight, const std::string &s) {
 HistogramService *h = &m_hSvc;
// h->h1D("RunNumber", "", s)->Fill(e.runNumber());  
 double Mmj = 0.;
 double Mej = 0.;
 double Mmjnu = 0.;
 double Mejnu = 0.;
 double MWj=0.;
 // double  MWb = 0;
 int sel_jet = -1;
 int b_jet = -1;
 int SjetsFarL = -1;
 int SjetsCloseL = -1;
 double dr = 1.5;
 double lpt =0.;
 int btags =0;   int btags1 =0;int btags2 =0;int btags3 =0;int btags4 =0;int btags5 =0; 
 int idx = -1;
 float mj0j1j2j3j4 = 0.; //mass of T
 float mj0j1j2 = 0.;
  //For MET
 float METx = e.met().Px();//[Missing E_T x direction]; // in MeV
 float METy = e.met().Py();//[Missing E_T y direction]; // in MeV
  /*std::cout<<""<<std::endl;
  std::cout<<"----In the PlotSemilep.cxx"<<std::endl;
  std::cout<<"METx= "<<METx<<" ;METy= "<<METy<<"  ;METz= "<<e.met().Pz()<<std::endl;*/
 float mass = 0; // mas of the lepton
 double Menu = 0;
 double Mmnu = 0;
 
 int SmaljetsFarLep = 0;
 int SmaljetsCloseLep = 0;
int fatjets = 0; int fatjets1=0; int fatjets2=0;
 int jetsmall = 0; int jetsmall1=0; int jetsmall2=0;

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
    getWFromLeptonicDecay(l, METx, METy, momNu, momW, mass);
  
    std::vector<LargeJet> sortedJets = e.largeJet(); //to sort the jets by mass
    std::sort(sortedJets.begin(), sortedJets.end(), sortFunction);
//  
//  std::vector<Jet> sortedJets2 = e.Jet(); //to sort the jets by mass
//  std::sort(sortedJets2.begin(), sortedJets2.end(), sortFunction2);
//  h->h1D("NSmallJets0", "", s)->Fill(e.jet().size(), weight);
//  h->h1D("Nfatjets0", "", s)->Fill(e.largeJet().size(), weight);
  vector<TLorentzVector> SjetsCloseLep;
  vector<TLorentzVector> SjetsFarLep;
  double mindr = 9999;
  for (int k = 0; k < e.jet().size(); ++k) {
   	double tdr = e.jet()[k].mom().DeltaR(l);
   	if(tdr < dr && e.jet()[k].mom().Perp() > lpt) {
  		    	sel_jet = k;
  		    	lpt = e.jet()[k].mom().Perp();
                       if (tdr < mindr) mindr = tdr;  
		    	}
  }

  
//  if (sel_jet == -1) return; 
  
   for (int k = 0; k < e.jet().size(); ++k) {
  	jetsmall++; //decomment for boosted
  	double tdr = e.jet()[k].mom().DeltaR(l);
   	if(tdr < dr){
   		  SjetsCloseL = k;
   		  if (SjetsCloseL == sel_jet) continue; 
  		  SjetsCloseLep.push_back(e.jet()[k].mom());
  		   SmaljetsCloseLep++;
   	}
    else if(tdr > dr){
    	SjetsFarL = k;
    	if (SjetsFarL == SjetsCloseL) continue;
        SjetsFarLep.push_back(e.jet()[k].mom());
        SmaljetsFarLep++;

          }     
  }

  
   float mtw = std::sqrt(2*l.Perp()*momNu.Perp()*(1 - std::cos(l.Phi() - momNu.Phi())) );
  // std::cout<<"momNu.Perp()= "<<momNu.Perp()<<"   momNu.Phi()="<<momNu.Phi()<<"     cos(l.Phi() - e.momNu.Phi())="<<std::cos(l.Phi() - e.momNu.Phi())<<std::endl;
   //std::cout<<"l.Perp()="<<l.Perp()<<"     l.Phi()="<<l.Phi()<<std::endl;
   if((momNu.Perp()<20e3)) std::cout<<"MET<20e3!!!  MET="<<(momNu.Perp())/1000.<<std::endl;
   if((mtw+momNu.Perp())<60e3) std::cout<<"mtw+MET<60e3!!!  mtw+MET="<<(momNu.Perp()+mtw)/1000.<<std::endl;
   //else std::cout<<"mtw="<<mtw/1000.<<"  MET="<<momNu.Perp()/1000.<<"   mtw+MET="<<(momNu.Perp()+mtw)/1000.<<std::endl;
   h->h1D("mtW", "", s)->Fill(mtw/1000., weight);
   h->h1D("MET_plus_mtW", "", s)->Fill((momNu.Perp()+mtw)/1000., weight);

  h->h1D("MET","",s)->Fill(momNu.Perp()/1000., weight);
  vector<TLorentzVector> Jet_Fat;
  for (int k = 0; k < e.largeJet().size(); k++) {
  	bool passdPhi=false;
    bool passdR = false;
    bool passM = false;
    bool passD12 = false;
    bool passEta = false;
    bool passPt = false;
    if (e.largeJet()[k].passLoose()) {
       passdR = e.jet()[sel_jet].mom().DeltaR(e.largeJet()[k].mom()) > 1.5;
      passM = e.largeJet()[k].mom().M() > 100e3;//60e3;
      passEta = true;//std::fabs(e.largeJet()[k].mom().Eta()) <  2.0;
      passPt = true;//e.largeJet()[k].mom().Perp() > 300e3;//160e3;
       passdPhi = true;//e.largeJet()[k].mom().DeltaPhi(l) > 2.3;
      passD12 = e.largeJet()[k].split12() > 40e3;
      if (passdR && passdPhi && passM  && passEta && passPt &&passD12) {
      	Jet_Fat.push_back(e.largeJet()[k].mom());
        fatjets++;
        idx = k;
      }
    }
  }
 //std::cout <<"PlotSemi: jetsmall= "<<jetsmall<<"  ;fatjets= "<<fatjets<<std::endl; 
 if (jetsmall < 4 && fatjets < 1) return;
//if(fatjets > 0) return;
//if (jetsmall < 4) return;
// if(fatjets < 1) return;
//  if(jetsmall < 4) std::cout <<"PlotSemi2: jetsmall= "<<jetsmall<<"  ;fatjets= "<<fatjets<<std::endl;
  const TLorentzVector &selj = e.jet()[sel_jet].mom();
  const TLorentzVector &j = e.jet()[0].mom();
 
/*  h->h1D("lepPt", "", s)->Fill(l.Perp()*1e-3, weight);
  h->h1D("jetPt", "", s)->Fill(j.Perp()*1e-3, weight);
  h->h1D("seljetPt", "", s)->Fill(selj.Perp()*1e-3, weight);
  h->h1D("NSmallJets", "", s)->Fill(jetsmall, weight);
  h->h1D("NSmaljetsCloseLep", "", s)->Fill(SmaljetsCloseLep, weight);
  h->h1D("NSmaljetsFarLep", "", s)->Fill(SmaljetsFarLep, weight);
  h->h1D("NLargeJets", "", s)->Fill(fatjets, weight);
  */
double mindr2 = 9999;
  for (int k = 0; k < e.jet().size(); ++k) {
double tdr = e.jet()[k].mom().DeltaR(l);
  	if (e.jet()[k].pass() && e.jet()[k].btag()) {
      btags++; b_jet = k;
     //   if(tdr <= 3.){
          if (tdr < mindr2) mindr2 = tdr;
       // }
  }
 }
  

  if (btags < 1) return;
const TLorentzVector &bj = e.jet()[b_jet].mom();
   h->h1D("NbJets", "", s)->Fill(btags, weight);
  h->h1D("lepPt", "", s)->Fill(l.Perp()*1e-3, weight);
  h->h1D("jetPt", "", s)->Fill(j.Perp()*1e-3, weight);
  h->h1D("seljetPt", "", s)->Fill(selj.Perp()*1e-3, weight);
  h->h1D("NSmallJets", "", s)->Fill(jetsmall, weight);
  h->h1D("NSmaljetsCloseLep", "", s)->Fill(SmaljetsCloseLep, weight);
  h->h1D("NSmaljetsFarLep", "", s)->Fill(SmaljetsFarLep, weight);
	  
//h->h1D("dR", "", s)->Fill(e.jet()[sel_jet].mom().DeltaR(l), weight);
h->h1D("dR", "", s)->Fill(mindr2, weight);
   h->h1D("NSelJets", "", s)->Fill(sel_jet, weight);
  if(mindr<0.4) h->h1D("lepPt2", "", s)->Fill(l.Perp()*1e-3, weight);
h->h1D("Eta", "", s)->Fill(l.Eta(), weight);  
//if(m_electron) h->h1D("Eta", "", s)->Fill(e.electron()[0].mom().Eta(), weight); 
  //else h->h1D("Eta", "", s)->Fill(e.muon()[0].mom().Eta(), weight); //done
std::cout <<"PlotSemi2: jetsmall= "<<jetsmall<<"  ;fatjets= "<<fatjets<<std::endl;

   if(m_electron){
  //      h->h1D("Eta", "", s)->Fill(e.electron()[0].mom().Eta(), weight);
       if (fatjets == 1) {	
    	    	  h->h1D("NSmallJets_chan_FatJet1", "", s)->Fill(jetsmall, weight);
    	    	  h->h1D("SmaljetsCloseLep_chan_FatJet1", "", s)->Fill(SmaljetsCloseLep, weight);
    	    	 // h->h1D("SmaljetsFarLep_chan_FatJet1", "", s)->Fill(SmaljetsFarLep, weight);
    	    	 // h->h1D("MassZprime_FatJet", "", s)->Fill((Jet_Fat[0]+e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).M()*1e-3, weight);
//                  h->h1D("PtZprime_FatJet", "", s)->Fill((Jet_Fat[0]+e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).Perp()*1e-3, weight);
//                  h->h1D("PtZprimeFatJet_W", "", s)->Fill((e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).Perp()*1e-3, weight);
//                  h->h1D("MassZprimeFatJet_FatJet0", "", s)->Fill((Jet_Fat[0]).M()*1e-3, weight);

    	    }
    	if(fatjets == 1 && SmaljetsCloseLep ==2){
				h->h1D("MassGStar_1FatJet_2jets", "", s)->Fill((Jet_Fat[0]+SjetsCloseLep[0]+SjetsCloseLep[1]+e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).M()*1e-3, weight);
				h->h1D("PtGStar_1FatJet_2jets", "", s)->Fill((Jet_Fat[0]+SjetsCloseLep[0]+SjetsCloseLep[1]+e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).Perp()*1e-3, weight);
				h->h1D("Pt1FatJet_2jets", "", s)->Fill((Jet_Fat[0]+SjetsCloseLep[0]+SjetsCloseLep[1]).Perp()*1e-3, weight);
				h->h1D("Pt1FatJet2jets_2jets", "", s)->Fill((SjetsCloseLep[0]+SjetsCloseLep[1]).Perp()*1e-3, weight);
				h->h1D("MassGStar1FatJet2jets_Fatjet2jets", "", s)->Fill((Jet_Fat[0]+SjetsCloseLep[0]+SjetsCloseLep[1]).M()*1e-3, weight); 
				h->h1D("Pt_1FatJet_2jets_FatJet0", "", s)->Fill((Jet_Fat[0]).Perp()*1e-3, weight);
				h->h1D("Mass_1FatJet_2jets_FatJet0", "", s)->Fill((Jet_Fat[0]).M()*1e-3, weight);
				h->h1D("Eta_1FatJet_2jets_FatJet0", "", s)->Fill((Jet_Fat[0]).Eta(), weight);
				h->h1D("Pt_1FatJet_2jets_jetsCloseLep0", "", s)->Fill((SjetsCloseLep[0]).Perp()*1e-3, weight); 
				h->h1D("Pt_1FatJet_2jets_jetsCloseLep1", "", s)->Fill((SjetsCloseLep[1]).Perp()*1e-3, weight); 
				
				h->h1D("Eta_1FatJet_2jets_jetsCloseLep0", "", s)->Fill((SjetsCloseLep[0]).Eta(), weight);
				h->h1D("Eta_1FatJet_2jets_jetsCloseLep1", "", s)->Fill((SjetsCloseLep[1]).Eta(), weight);
			}
		else if(fatjets > 1) {
				h->h1D("NSmallJets_chan_FatJet2", "", s)->Fill(jetsmall, weight);
				h->h1D("SmaljetsCloseLep_chan_FatJet2", "", s)->Fill(SmaljetsCloseLep, weight);
				h->h1D("SmaljetsFarLep_chan_FatJet2", "", s)->Fill(SmaljetsFarLep, weight);
				h->h1D("MassGStar_2FatJet", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]+e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).M()*1e-3, weight);
				h->h1D("PtGStar_2FatJet", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]+e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).Perp()*1e-3, weight);
				h->h1D("Pt2FatJet", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]).Perp()*1e-3, weight);
				h->h1D("Mass2FatJet_FatJet0", "", s)->Fill(Jet_Fat[0].M()*1e-3, weight);
				h->h1D("Pt2FatJet_FatJet0", "", s)->Fill(Jet_Fat[0].Perp()*1e-3, weight);
				h->h1D("Mass2FatJet_FatJet1", "", s)->Fill(Jet_Fat[1].M()*1e-3, weight);
				h->h1D("Pt2FatJet_FatJet1", "", s)->Fill(Jet_Fat[1].Perp()*1e-3, weight);
				h->h1D("Eta2FatJet_FatJet0", "", s)->Fill(Jet_Fat[0].Eta(), weight);
				h->h1D("Eta2FatJet_FatJet1", "", s)->Fill(Jet_Fat[1].Eta(), weight);                  	    	     
				h->h1D("MassGStar2FatJet_2FatJets", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]).M()*1e-3, weight);  
			
			}     
    	            	    
//    else if(fatjets == 0){
//    	h->h1D("NSmallJets_chan_FatJet0", "", s)->Fill(jetsmall, weight);
//    	h->h1D("SmaljetsCloseLep_chan_FatJet0", "", s)->Fill(SmaljetsCloseLep, weight);
//    	h->h1D("SmaljetsFarLep_chan_FatJet0", "", s)->Fill(SmaljetsFarLep, weight);    	  
//    }

//   else if (fatjets > 1) {
//            	  h->h1D("NSmallJets_chan_FatJet2", "", s)->Fill(jetsmall, weight);
//    	    	  h->h1D("SmaljetsCloseLep_chan_FatJet2", "", s)->Fill(SmaljetsCloseLep, weight);
//    	    	  h->h1D("SmaljetsFarLep_chan_FatJet2", "", s)->Fill(SmaljetsFarLep, weight);
////    	    	  h->h1D("MassGStar_2FatJet", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]+e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).M()*1e-3, weight);
////    	    	  h->h1D("PtGStar_2FatJet", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]+e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).Perp()*1e-3, weight);
////    	          h->h1D("Pt2FatJet", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]).Perp()*1e-3, weight);
////    	          h->h1D("Mass2FatJet_FatJet0", "", s)->Fill(Jet_Fat[0].M()*1e-3, weight);
////                  h->h1D("Pt2FatJet_FatJet0", "", s)->Fill(Jet_Fat[0].Perp()*1e-3, weight);
////                  h->h1D("Mass2FatJet_FatJet1", "", s)->Fill(Jet_Fat[1].M()*1e-3, weight);
////                  h->h1D("Pt2FatJet_FatJet1", "", s)->Fill(Jet_Fat[1].Perp()*1e-3, weight);
////                  h->h1D("Eta2FatJet_FatJet0", "", s)->Fill(Jet_Fat[0].Eta(), weight);
////                  h->h1D("Eta2FatJet_FatJet1", "", s)->Fill(Jet_Fat[1].Eta(), weight);                  	    	     
////                  h->h1D("MassGStar2FatJet_2FatJets", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]).M()*1e-3, weight); 
//
//    	    }
    }
   else{
//    h->h1D("Eta", "", s)->Fill(e.muon()[0].mom().Eta(), weight); 
    if (fatjets == 1) {	
    	h->h1D("NSmallJets_chan_FatJet1", "", s)->Fill(jetsmall, weight);
    	h->h1D("SmaljetsCloseLep_chan_FatJet1", "", s)->Fill(SmaljetsCloseLep, weight);
    	//h->h1D("SmaljetsFarLep_chan_FatJet1", "", s)->Fill(SmaljetsFarLep, weight);
    	//h->h1D("MassZprime_FatJet", "", s)->Fill((Jet_Fat[0]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).M()*1e-3, weight);
//                  h->h1D("PtZprime_FatJet", "", s)->Fill((Jet_Fat[0]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
//                  h->h1D("PtZprimeFatJet_W", "", s)->Fill((e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
//                  h->h1D("MassZprimeFatJet_FatJet0", "", s)->Fill((Jet_Fat[0]).M()*1e-3, weight);

    	    } 
     if(fatjets == 1 && SmaljetsCloseLep ==2){
				h->h1D("MassGStar_1FatJet_2jets", "", s)->Fill((Jet_Fat[0]+SjetsCloseLep[0]+SjetsCloseLep[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).M()*1e-3, weight);
				h->h1D("PtGStar_1FatJet_2jets", "", s)->Fill((Jet_Fat[0]+SjetsCloseLep[0]+SjetsCloseLep[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
				h->h1D("Pt1FatJet_2jets", "", s)->Fill((Jet_Fat[0]+SjetsCloseLep[0]+SjetsCloseLep[1]).Perp()*1e-3, weight);
				h->h1D("Pt1FatJet2jets_2jets", "", s)->Fill((SjetsCloseLep[0]+SjetsCloseLep[1]).Perp()*1e-3, weight);
				h->h1D("MassGStar1FatJet2jets_Fatjet2jets", "", s)->Fill((Jet_Fat[0]+SjetsCloseLep[0]+SjetsCloseLep[1]).M()*1e-3, weight); 
				h->h1D("Pt_1FatJet_2jets_FatJet0", "", s)->Fill((Jet_Fat[0]).Perp()*1e-3, weight);
				h->h1D("Mass_1FatJet_2jets_FatJet0", "", s)->Fill((Jet_Fat[0]).M()*1e-3, weight);
				h->h1D("Eta_1FatJet_2jets_FatJet0", "", s)->Fill((Jet_Fat[0]).Eta(), weight);
				h->h1D("Pt_1FatJet_2jets_jetsCloseLep0", "", s)->Fill((SjetsCloseLep[0]).Perp()*1e-3, weight); 
				h->h1D("Pt_1FatJet_2jets_jetsCloseLep1", "", s)->Fill((SjetsCloseLep[1]).Perp()*1e-3, weight);
				h->h1D("Eta_1FatJet_2jets_jetsCloseLep0", "", s)->Fill((SjetsCloseLep[0]).Eta(), weight);
				h->h1D("Eta_1FatJet_2jets_jetsCloseLep1", "", s)->Fill((SjetsCloseLep[1]).Eta(), weight);  	    	     
			}
	  if(fatjets > 1) {
				h->h1D("MassGStar_2FatJet", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).M()*1e-3, weight);
				h->h1D("PtGStar_2FatJet", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
				h->h1D("Pt2FatJet", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]).Perp()*1e-3, weight);
				h->h1D("Mass2FatJet_FatJet0", "", s)->Fill(Jet_Fat[0].M()*1e-3, weight);
				h->h1D("Pt2FatJet_FatJet0", "", s)->Fill(Jet_Fat[0].Perp()*1e-3, weight);
				h->h1D("Mass2FatJet_FatJet1", "", s)->Fill(Jet_Fat[1].M()*1e-3, weight);
				h->h1D("Pt2FatJet_FatJet1", "", s)->Fill(Jet_Fat[1].Perp()*1e-3, weight);
				h->h1D("Eta2FatJet_FatJet0", "", s)->Fill(Jet_Fat[0].Eta(), weight);
				h->h1D("Eta2FatJet_FatJet1", "", s)->Fill(Jet_Fat[1].Eta(), weight);                  	    	     
				h->h1D("MassGStar2FatJet_2FatJets", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]).M()*1e-3, weight);  
			
			}  	    
//    else if(fatjets == 0){
//    	h->h1D("NSmallJets_chan_FatJet0", "", s)->Fill(jetsmall, weight);
//    	h->h1D("SmaljetsCloseLep_chan_FatJet0", "", s)->Fill(SmaljetsCloseLep, weight);
//    	h->h1D("SmaljetsFarLep_chan_FatJet0", "", s)->Fill(SmaljetsFarLep, weight);                            	    	  
//    }
 
 
//    else if (fatjets > 1) {
//    	h->h1D("NSmallJets_chan_FatJet2", "", s)->Fill(jetsmall, weight);
//    	h->h1D("SmaljetsCloseLep_chan_FatJet2", "", s)->Fill(SmaljetsCloseLep, weight);
//    	h->h1D("SmaljetsFarLep_chan_FatJet2", "", s)->Fill(SmaljetsFarLep, weight); 
////    	h->h1D("MassGStar_2FatJet", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).M()*1e-3, weight);
////    	    	h->h1D("PtGStar_2FatJet", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
////    	        h->h1D("Pt2FatJet", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]).Perp()*1e-3, weight);
////    	        h->h1D("Mass2FatJet_FatJet0", "", s)->Fill(Jet_Fat[0].M()*1e-3, weight);
////                h->h1D("Pt2FatJet_FatJet0", "", s)->Fill(Jet_Fat[0].Perp()*1e-3, weight);
////                h->h1D("Mass2FatJet_FatJet1", "", s)->Fill(Jet_Fat[1].M()*1e-3, weight);
////                h->h1D("Pt2FatJet_FatJet1", "", s)->Fill(Jet_Fat[1].Perp()*1e-3, weight);
////                h->h1D("Eta2FatJet_FatJet0", "", s)->Fill(Jet_Fat[0].Eta(), weight);
////                h->h1D("Eta2FatJet_FatJet1", "", s)->Fill(Jet_Fat[1].Eta(), weight);                  	    	     
////                h->h1D("MassGStar2FatJet_2FatJets", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]).M()*1e-3, weight);
//    	    }
    } //Muon's Loop
  //std::cout<<"at the end: momNu.Perp()*1e-3= "<<momNu.Perp()*1e-3<<std::endl;
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
 
