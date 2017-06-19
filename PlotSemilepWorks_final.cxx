#include "Plot.h"
#include "PlotSemilep.h"
#include "Event.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include "HistogramService.h"
#include <iostream>
using namespace std;


PlotSemilep::PlotSemilep(const std::string &filename, bool electron, const std::vector<std::string> &systs)
  : Plot(filename, systs), m_electron(electron) {
  m_hSvc.create1D("lepPt", "; Lepton p_{T} (GeV) ; Events", 100, 0, 1000);
  m_hSvc.create1D("MassMuj", "; Invariant Mass M(jet-Muon) (GeV)", 100, 0, 1000);
  m_hSvc.create1D("MassElecj", "; Invariant Mass M(jet-Electron) (GeV)", 100, 0, 1000);
  m_hSvc.create1D("MET", "; Missing ET (MeV)", 100, 0, 1000);
  m_hSvc.create1D("MassElecjNu", "; SM Top Mass (electron channel) (GeV)", 200, 0, 200);
  m_hSvc.create1D("MassSMTop", "; SM Top Mass(GeV)",50, 0, 500);
  m_hSvc.create1D("MassSMTop2", "; SM Top Mass(GeV) where we have Wb", 25, 0, 1000);
  m_hSvc.create1D("MassWj", "; Invariant Mass M(W+j) (GeV)", 50, 0, 500);
  m_hSvc.create1D("MassWb", "; Invariant Mass M(W+b) (GeV)", 25, 0, 1000);
  m_hSvc.create1D("MassTprime", "; VLQ Top Mass (jet channel) (GeV)", 27, 0, 2700);
  m_hSvc.create1D("MassZprime", "; Z' (GeV)", 27, 0, 2700);
  m_hSvc.create1D("jetPt", "; Jet p_{T} (GeV) ; Events", 100, 0, 1000);
  m_hSvc.create1D("largeJetPt", "; Large jet p_{T} (GeV); Events", 100, 0, 1000);
  m_hSvc.create1D("largeJetM", "; Large jet M (GeV); Events", 100, 0, 1000);
  m_hSvc.create1D("subjetPt", "; Sub-jet p_{T} (GeV) ; Events", 100, 0, 1000);
  m_hSvc.create1D("subjetN", "; Sub-jet multiplicity ; Events", 6, 0.5, 6.5);
  m_hSvc.create1D("NLargeJets", "; Number of Fat jets ; Events", 7, -0.5, 6.5);
  m_hSvc.create1D("NSmallJets", "; Number of Small jets ; Events", 21,-0.5,20.5);
  m_hSvc.create1D("NbJets", "; Number of b-jets ; Events", 13,-0.5,12.5);
  
  
  m_hSvc.create1D("MassFatJet", "; Leading Fat Mass (GeV)", 8, 0, 800);
  m_hSvc.create1D("MassGStar_3jets", "; GStar Mass (with 3 jets DR > 1.5) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassGStar_2jets", "; GStar Mass (with 2 jets DR > 1.5) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassZprime_3jets", "; Z',T , gkk mass (with 3 jets DR > 1.5) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassZprime_2jets", "; Z', T , gkk mass (with 2 jets DR > 1.5) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassZprime_FatJet", "; Z', T, gkk mass with a fat Jet (GeV)", 36, 0, 3600);
  
   m_hSvc.create1D("jetPt", "; Jet p_{T} (GeV) ; Events", 100, 0, 1000);
  m_hSvc.create1D("largeJetPt", "; Large jet p_{T} (GeV); Events", 100, 0, 1000);
  m_hSvc.create1D("largeJetM", "; Large jet M (GeV); Events", 100, 0, 1000);
  m_hSvc.create1D("subjetPt", "; Sub-jet p_{T} (GeV) ; Events", 100, 0, 1000);
  m_hSvc.create1D("subjetN", "; Sub-jet multiplicity ; Events", 6, 0.5, 6.5);
  m_hSvc.create1D("NLargeJets", "; Number of Fat jets ; Events", 7, -0.5, 6.5);
  m_hSvc.create1D("NSmallJets", "; Number of Small jets ; Events", 13,-0.5,12.5);
  m_hSvc.create1D("NbJets", "; Number of b-jets ; Events", 13,-0.5,12.5);
  
  m_hSvc.create1D("MassGStar_2FatJet", "; GStar Mass (with 2 fat jets + ttbar) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassGStar_1FatJet_2jets", "; GStar Mass (with 1 fat jet and 2 small jets + ttbar) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassGStar_4jets", "; GStar Mass (with 4 small jets + ttbar) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassGStar_5jets", "; GStar Mass (with 5 small jets + ttbar) (GeV)", 36, 0, 3600);
  
  m_hSvc.create1D("seljetPt", "; Selected Jet p_{T} (GeV) ; Events", 20, 0, 2000);
  m_hSvc.create1D("PtFatJet_0", "; FatJet_0 p_{T} (GeV) ; Events", 20, 0, 2000);
  m_hSvc.create1D("PtZprime_FatJet", "; Zprime_FatJet p_{T} (GeV) ; Events", 20, 0, 2000);
  m_hSvc.create1D("PtZprime_2jets", "; Zprime_2jets p_{T} (GeV) ; Events", 20, 0, 2000);
  m_hSvc.create1D("PtGStar_1FatJet_2jets", "; GStar_1FatJet_2jets p_{T} (GeV) ; Events", 40, 0, 4000);
  m_hSvc.create1D("PtZprime_3jets", "; Zprime_3jets p_{T} (GeV) ; Events", 20, 0, 2000);
  m_hSvc.create1D("PtGStar_4jets", "; GStar_4jets p_{T} (GeV) ; Events", 20, 0, 2000);
  m_hSvc.create1D("PtGStar_5jets", "; GStar_5jets p_{T} (GeV) ; Events", 20, 0, 2000);
  m_hSvc.create1D("PtGStar_2FatJet", "; GStar_2FatJet p_{T} (GeV) ; Events", 20, 0, 2000);
  
  m_hSvc.create1D("PtFatJet", "; FatJet p_{T} (GeV) ; Events", 20, 0, 2000);
  m_hSvc.create1D("Pt2jets", "; 2jets p_{T} (GeV) ; Events", 20, 0, 2000);
  m_hSvc.create1D("Pt1FatJet_2jets", "; 1FatJet_2jets p_{T} (GeV) ; Events", 20, 0, 2000);
  m_hSvc.create1D("Pt3jets", "; 3jets p_{T} (GeV) ; Events", 20, 0, 2000);
  m_hSvc.create1D("Pt4jets", "; 4jets p_{T} (GeV) ; Events", 20, 0, 2000);
  m_hSvc.create1D("Pt5jets", "; jets p_{T} (GeV) ; Events", 20, 0, 2000);
  m_hSvc.create1D("Pt2FatJet", "; FatJet p_{T} (GeV) ; Events", 20, 0, 2000);

  
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
  std::cout<<""<<std::endl;
  std::cout<<"----In the PlotSemilep.cxx"<<std::endl;
  std::cout<<"METx= "<<METx<<" ;METy= "<<METy<<"  ;METz= "<<e.met().Pz()<<std::endl;
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
  getWFromLeptonicDecay(l, METx, METy, momNu, momW, mass);
 std::cout<<"after getWfrom: momNu.Perp()*1e-3= "<<momNu.Perp()*1e-3<<std::endl;
   for (int k = 0; k < e.jet().size(); ++k) {
  	//std::cout<<"In the jet loop k= "<<k<<std::endl;
  	if (e.jet()[k].pass()) {
  		jetsmall++;
  	}
  	//std::cout<<"after e.jet()[k].pass"<<std::endl;
    double tdr = e.jet()[k].mom().DeltaR(l);
    if (tdr < dr && e.jet()[k].mom().Perp() > lpt) {
    	sel_jet = k;
        lpt = e.jet()[k].mom().Perp();
        }
   }   
   //std::cout<<"1: sel_jet= "<<sel_jet<<std::endl; 
   if (sel_jet == -1) return;
  //std::cout<<"sel_jet= "<<sel_jet<<std::endl;
  
     
    if(m_electron){
    	 Mej = (e.electron()[0].mom()+e.jet()[sel_jet].mom()).M();
    	 //Menu = (e.electron()[0].mom()+e.momNu).M();
    }
    else{
    	 Mmj = (e.muon()[0].mom() + e.jet()[sel_jet].mom()).M();
    	 //Mmnu = (e.muon()[0].mom() +e.momNu).M();
    }
    MWj = momW.M()+e.jet()[sel_jet].mom().M();
    std::cout<<"MassW=momW.M()= "<<momW.M()<<std::endl;
    std::cout<<""<<std::endl;
    h->h1D("MassMuj","",s)->Fill(Mmj*1e-3, weight);
    h->h1D("MassElecj","",s)->Fill(Mej*1e-3, weight);
    h->h1D("MET","",s)->Fill(momNu.Perp(), weight);
    h->h1D("MassSMTop","",s)->Fill(MWj*1e-3, weight); //h->h1D("MassSMTop","",s)->Fill((MWj+e.electron()[0].mom().M())*1e-3, weight);
    h->h1D("MassElecjNu","",s)->Fill(Mej*1e-3, weight);
    h->h1D("MassWj", "", s)->Fill(MWj*1e-3, weight);
    
  for (int k = 0; k < e.largeJet().size(); k++) {
  	bool passdPhi=false;
    bool passdR = false;
    bool passM = false;
    bool passD12 = false;
    if (e.largeJet()[k].pass()) {
       passdR = e.jet()[sel_jet].mom().DeltaR(e.largeJet()[k].mom()) > 1.5;
      if(m_electron){
         passdPhi=e.largeJet()[k].mom().DeltaPhi(e.electron()[0].mom()) > 2.3;
      }
      else {
      	passdPhi=e.largeJet()[k].mom().DeltaPhi(e.muon()[0].mom()) > 2.3;
      }
      passM = e.largeJet()[k].mom().M() > 60e3;
      passD12 = e.largeJet()[k].split12() > 20e3;
      if (passdR && passdPhi && passM && passD12) {
        fatjets++;
        idx = k;
      }
    }
  }
 // if (fatjets < 1) return;   
  if (goodjet < 4 && fatjets < 1) return;
  vector<TLorentzVector> goodJets_Fat;
   for (int k = 0; k < e.largeJet().size(); k++) {
  	bool passdPhi=false;
    bool passdR = false;
    bool passM = false;
    bool passD12 = false;
    if (e.largeJet()[k].pass()) {
       passdR = e.jet()[sel_jet].mom().DeltaR(e.largeJet()[k].mom()) > 1.5;
      if(m_electron){
         passdPhi=e.largeJet()[k].mom().DeltaPhi(e.electron()[0].mom()) > 2.3;
      }
      else {
      	passdPhi=e.largeJet()[k].mom().DeltaPhi(e.muon()[0].mom()) > 2.3;
      }
      passM = e.largeJet()[k].mom().M() > 60e3;
      passD12 = e.largeJet()[k].split12() > 20e3;
      if (passdR && passdPhi && passM && passD12) {
        goodJets_Fat.push_back(e.largeJet()[k].mom());
      }
    }
  }

  const TLorentzVector &selj = e.jet()[sel_jet].mom();
  
  const TLorentzVector &j = e.jet()[0].mom();
  const TLorentzVector &lj = e.largeJet()[idx].mom();
  h->h1D("lepPt", "", s)->Fill(l.Perp()*1e-3, weight);
  h->h1D("jetPt", "", s)->Fill(j.Perp()*1e-3, weight);
  h->h1D("seljetPt", "", s)->Fill(selj.Perp()*1e-3, weight);
  h->h1D("NSmallJets", "", s)->Fill(jetsmall, weight);
  h->h1D("largeJetPt", "", s)->Fill(lj.Perp()*1e-3, weight);
  h->h1D("largeJetM", "", s)->Fill(lj.M()*1e-3, weight);
  h->h1D("NLargeJets", "", s)->Fill(fatjets, weight);
 
  vector<TLorentzVector> goodJets;
  for (int k = 0; k < e.jet().size(); ++k) {
  	double tdr = e.jet()[k].mom().DeltaR(l);
  	if (e.jet()[k].pass() && tdr > dr) {
  		goodJets.push_back(e.jet()[k].mom());
  		goodjet++;
  		}
  	}
   float mj0j1j2j3j4 = (goodJets[0] + goodJets[1]+ goodJets[2]+ goodJets[3]+ goodJets[4]).M(); //mass of T
   float mj0j1j2 = (goodJets[0] + goodJets[1]+ goodJets[2]).M(); //mass of T
   if(jetsmall >= 5){
   	std::cout<<"jetsmall>=5: n_jetSamll= "<<jetsmall<<std::endl;
   	 h->h1D("MassTprime", "", s)->Fill(mj0j1j2j3j4*1e-3, weight);
   }
   
  
 

 // if ((MWb+Menu)>150*1e-3) std::cout<<"TopMass2= "<<(MWb+Menu)*1e-3<<std::endl;
  //h->h1D("MassSMTop2","",s)->Fill((goodJets_b[0].M()+Menu)*1e-3, weight);
   
  // if (momW.Perp() != 0) std::cout<<"MassWb= "<<MWb<<" PtW= "<<momW.Perp()<<std::endl;
 
  h->h1D("MassFatJet", "", s)->Fill(goodJets_Fat[0].M()*1e-3, weight);
   h->h1D("PtFatJet_0", "", s)->Fill(goodJets_Fat[0].Perp()*1e-3, weight);
  
  //fatjet + lepton +MET + highest P_T jet 
   std::cout<<"Before Z' calc: momNu.Perp()*1e-3= "<<momNu.Perp()*1e-3<<std::endl;
   if(m_electron){
    	 Mej = (e.electron()[0].mom()+e.jet()[sel_jet].mom()).M();
    	 h->h1D("MassZprime", "", s)->Fill((goodJets_Fat[0]+goodJets[0]+e.electron()[0].mom()).M()*1e-3+momNu.Perp()*1e-3, weight);
    	 //Menu = (e.electron()[0].mom()+e.momNu).M();
    }
   else{
//    	 Mmj = (e.muon()[0].mom() + e.jet()[sel_jet].mom()).M();
//    	 h->h1D("MassZprime", "", s)->Fill((goodJets_Fat[0]+goodJets[0]+e.muon()[0].mom()).M()*1e-3+momNu.Perp()*1e-3, weight);
    	 //Mmnu = (e.muon()[0].mom() +e.momNu).M();
    	 std::cout<<"==============in mu case fatjets= "<<fatjets<<" goodjet= "<<goodjet<<std::endl;
    	    Mmj = (e.muon()[0].mom() + e.jet()[sel_jet].mom()).M(); 
    	    h->h1D("MassZprime_FatJet", "", s)->Fill((goodJets_Fat[0]+e.jet()[sel_jet].mom()+e.muon()[0].mom()).M()*1e-3+momNu.Perp()*1e-3, weight);
            h->h1D("PtZprime_FatJet", "", s)->Fill((goodJets_Fat[0]+e.jet()[sel_jet].mom()+e.muon()[0].mom()).Perp()*1e-3+momNu.Perp()*1e-3, weight);   	    
    	   // h->h1D("PtFatJet", "", s)->Fill((goodJets_Fat[0]).Perp()*1e-3, weight);
    	   
            if (goodjet == 2){
    	    	 h->h1D("MassZprime_2jets", "", s)->Fill((goodJets[0]+goodJets[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()).M()*1e-3+momNu.Perp()*1e-3, weight);
    	    	 h->h1D("PtZprime_2jets", "", s)->Fill((goodJets[0]+goodJets[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()).Perp()*1e-3+momNu.Perp()*1e-3, weight);
    	    	 h->h1D("Pt2jets", "", s)->Fill((goodJets[0]+goodJets[1]).Perp()*1e-3, weight);
    	    	 std::cout<<"Good jet ==1:  sel_fatjets= "<<fatjets<<"  Ngoodjet="<<goodjet<<std::endl;
                 if(fatjets == 1){
    	    	 	 h->h1D("MassGStar_1FatJet_2jets", "", s)->Fill((goodJets_Fat[0]+goodJets[0]+goodJets[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()).M()*1e-3+momNu.Perp()*1e-3, weight);
    	    	 	 h->h1D("PtGStar_1FatJet_2jets", "", s)->Fill((goodJets_Fat[0]+goodJets[0]+goodJets[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()).Perp()*1e-3+momNu.Perp()*1e-3, weight);
    	    	 	 std::cout<<"fatjet==1case: sel_fatjets= "<<fatjets<<" goodjet="<<goodjet<<"  Pt1FatJet_2jets="<<((goodJets_Fat[0]+goodJets[0]+goodJets[1]).Perp()*1e-3)<<std::endl;
    	    	     h->h1D("Pt1FatJet_2jets", "", s)->Fill((goodJets_Fat[0]+goodJets[0]+goodJets[1]).Perp()*1e-3, weight);
                 }
    	    }
    	 	if (goodjet > 2){
    	 		std::cout<<"----Good jet >2 :  fatjets= "<<fatjets<<"  Ngoodjet="<<goodjet<<std::endl;
    	 		if (goodjet == 3){
    	 			 h->h1D("MassZprime_3jets", "", s)->Fill((goodJets[0]+goodJets[1]+goodJets[3]+e.jet()[sel_jet].mom()+e.muon()[0].mom()).M()*1e-3+momNu.Perp()*1e-3, weight);
    	 			 h->h1D("PtZprime_3jets", "", s)->Fill((goodJets[0]+goodJets[1]+goodJets[3]+e.jet()[sel_jet].mom()+e.muon()[0].mom()).Perp()*1e-3+momNu.Perp()*1e-3, weight);
    	 			 std::cout<<"goodjet==3: case: sel_fatjets= "<<fatjets<<" goodjet="<<goodjet<<"  Pt3jets="<<((goodJets[0]+goodJets[1]+goodJets[3]).Perp()*1e-3)<<std::endl;
    	 		     h->h1D("Pt3jets", "", s)->Fill((goodJets[0]+goodJets[1]+goodJets[3]).Perp()*1e-3, weight);
                          }
    	 		 if (goodjet == 4){
    	 		 	 h->h1D("MassGStar_4jets", "", s)->Fill((goodJets[0]+goodJets[1]+goodJets[2]+goodJets[3]+e.jet()[sel_jet].mom()+e.muon()[0].mom()).M()*1e-3+momNu.Perp()*1e-3, weight);
    	 		 	 h->h1D("PtGStar_4jets", "", s)->Fill((goodJets[0]+goodJets[1]+goodJets[2]+goodJets[3]+e.jet()[sel_jet].mom()+e.muon()[0].mom()).Perp()*1e-3+momNu.Perp()*1e-3, weight);
    	 		         h->h1D("Pt4jets", "", s)->Fill((goodJets[0]+goodJets[1]+goodJets[2]+goodJets[3]).Perp()*1e-3, weight);  
                         }
    	 		 if (goodjet > 4) {
    	 		 	h->h1D("MassGStar_5jets", "", s)->Fill((goodJets[0]+goodJets[1]+goodJets[2]+goodJets[3]+goodJets[4]+e.jet()[sel_jet].mom()+e.muon()[0].mom()).M()*1e-3+momNu.Perp()*1e-3, weight);
    	 		 	h->h1D("PtGStar_5jets", "", s)->Fill((goodJets[0]+goodJets[1]+goodJets[2]+goodJets[3]+goodJets[4]+e.jet()[sel_jet].mom()+e.muon()[0].mom()).Perp()*1e-3+momNu.Perp()*1e-3, weight);
    	 		        h->h1D("Pt5jets", "", s)->Fill((goodJets[0]+goodJets[1]+goodJets[2]+goodJets[3]+goodJets[4]).Perp()*1e-3, weight);
                         }
    	 	}
    	    if(fatjets > 1) {
    	    	h->h1D("MassGStar_2FatJet", "", s)->Fill((goodJets_Fat[0]+goodJets_Fat[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()).M()*1e-3+momNu.Perp()*1e-3, weight);
    	    	h->h1D("PtGStar_2FatJet", "", s)->Fill((goodJets_Fat[0]+goodJets_Fat[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()).Perp()*1e-3+momNu.Perp()*1e-3, weight);
    	        h->h1D("PtFatJet", "", s)->Fill((goodJets_Fat[0]+goodJets_Fat[1]).Perp()*1e-3, weight);   
          //std::cout<<"2fatjets= "<<fatjets<<std::endl;
    	    }
    	    if(fatjets == 0 || fatjets == 1) {
    	    	std::cout<<"fatjets_inf1= "<<fatjets<<"    ;goodjet= "<<goodjet<<std::endl;
    	 	   h->h1D("MassGStar_1FatJet_2jets", "", s)->Fill((goodJets_Fat[0]+goodJets[0]+goodJets[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()).M()*1e-3+momNu.Perp()*1e-3, weight);
    	       if (goodjet == 4) h->h1D("MassGStar_4jets", "", s)->Fill((goodJets[0]+goodJets[1]+goodJets[2]+goodJets[3]+e.jet()[sel_jet].mom()+e.muon()[0].mom()).M()*1e-3+momNu.Perp()*1e-3, weight);
    	       if (goodjet > 4) h->h1D("MassGStar_5jets", "", s)->Fill((goodJets[0]+goodJets[1]+goodJets[2]+goodJets[3]+goodJets[4]+e.jet()[sel_jet].mom()+e.muon()[0].mom()).M()*1e-3+momNu.Perp()*1e-3, weight); 
    	    } 
    }
  std::cout<<"at the end: momNu.Perp()*1e-3= "<<momNu.Perp()*1e-3<<std::endl;
   for (int k = 0; k < e.jet().size(); ++k) {
  	if (e.jet()[k].pass() && e.jet()[k].btag()) {
      btags++;
      }
   }
  if (btags < 1) return;
  
  vector<TLorentzVector> goodJets_b;
  for (int k = 0; k < e.jet().size(); ++k) {
  	if (e.jet()[k].pass() && e.jet()[k].btag()) goodJets_b.push_back(e.jet()[k].mom());
  	}
  float MWb = (goodJets_b[0] + momW).M(); //mass of T
   h->h1D("NbJets", "", s)->Fill(btags, weight);
   h->h1D("MassWb", "", s)->Fill(MWb*1e-3, weight);
    h->h1D("MassSMTop2","",s)->Fill(MWb*1e-3, weight);
      if(jetsmall < 5) {
   	std::cout<<"jetsmall<5: n_jetSamll= "<<jetsmall<<std::endl;
   	h->h1D("MassTprime", "", s)->Fill((mj0j1j2+MWb)*1e-3, weight);
   }



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
