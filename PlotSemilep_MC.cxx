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

bool sortFunction(const LargeJet &a, const LargeJet &b){
  return a.mom().M() > b.mom().M();
}
//bool sortFunction2(const Jet &a, const Jet &b){
//  return a.mom().Perp() > b.mom().Perp();
//}
PlotSemilep::PlotSemilep(const std::string &filename, bool electron, const std::vector<std::string> &systs)
  : Plot(filename, systs), m_electron(electron) {
  m_hSvc.create1D("MassMuj", "; Invariant Mass M(jet-Muon) (GeV)", 100, 0, 1000);
  m_hSvc.create1D("MassElecj", "; Invariant Mass M(jet-Electron) (GeV)", 100, 0, 1000);
  m_hSvc.create1D("MassElecjNu", "; SM Top Mass (electron channel) (GeV)", 200, 0, 200);
  m_hSvc.create1D("MassWj", "; Invariant Mass M(W+j) (GeV)", 50, 0, 500);
  m_hSvc.create1D("MassWb", "; Invariant Mass M(W+b) (GeV)", 25, 0, 1000);
  m_hSvc.create1D("MassTprime", "; VLQ Top Mass (jet channel) (GeV)", 27, 0, 2700);

  
  m_hSvc.create1D("lepPt", "; Lepton p_{T} (GeV) ; Events", 100, 0, 1000);
  m_hSvc.create1D("MET", "; Missing ET (MeV)", 100, 0, 1000);
  m_hSvc.create1D("MassSMTop", "; SM Top Mass(GeV)",50, 0, 500);
  m_hSvc.create1D("MassSMTop2", "; SM Top Mass(GeV) where we have Wb", 25, 0, 1000);
  m_hSvc.create1D("jetPt_0", "; First Jet p_{T} (GeV) ; Events", 100, 0, 1000);
  m_hSvc.create1D("largeJetPt_0", "; First Large jet p_{T} (GeV); Events", 1000, 0, 1000);
  m_hSvc.create1D("MassFatJet_0", "; First Large jet M (GeV); Events", 100, 0, 1000);
  m_hSvc.create1D("largeJetM_0", "; First Large jet M (GeV); Events", 1000, 0, 1000);
  m_hSvc.create1D("subjetPt", "; Sub-jet p_{T} (GeV) ; Events", 100, 0, 1000);
  m_hSvc.create1D("subjetN", "; Sub-jet multiplicity ; Events", 6, 0.5, 6.5);
  m_hSvc.create1D("NLargeJets", "; Number of Fat jets ; Events", 7, -0.5, 6.5);
  m_hSvc.create1D("NSmallJets", "; Number of Small jets ; Events", 21,-0.5,20.5);
  m_hSvc.create1D("NbJets", "; Number of b-jets ; Events", 13,-0.5,12.5);
  m_hSvc.create1D("MassFatJet_1", "; Second Large jet M (GeV); Events", 100, 0, 1000);
  
  m_hSvc.create1D("largeJetPt_1", "; Second Large jet p_{T} (GeV); Events", 1000, 0, 1000);
  m_hSvc.create1D("largeJetM_1", "; Second Large jet M (GeV); Events", 1000, 0, 1000);
  
 // m_hSvc.create1D("NSmaljetsCloseLep", "; Number of Small jets close to the lepton ; Events", 21,-0.5,20.5);
  m_hSvc.create1D("NSmaljetsFarLepLJet", "; Number of Small jets Far from the lepton and the fat jet; Events", 21,-0.5,20.5);
  
  m_hSvc.create1D("NSmaljetsFarLep", "; Number of Small jets far from the lepton ; Events", 21,-0.5,20.5);   
  
  m_hSvc.create1D("MassGStar_3jets", "; GStar Mass (with 3 jets DR > 1.5) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassGStar_2jets", "; GStar Mass (with 2 jets DR > 1.5) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassZprime_3jets", "; Z',T , gkk mass (with 3 jets DR > 1.5) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassZprime_2jets", "; Z', T , gkk mass (with 2 jets DR > 1.5) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassZprime_FatJet", "; Z', T, gkk mass with a fat Jet (GeV)", 36, 0, 3600);
    
  m_hSvc.create1D("MassGStar_2FatJet", "; GStar Mass (with 2 fat jets + ttbar) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassGStar_1FatJet_2jets", "; GStar Mass (with 1 fat jet and 2 small jets + ttbar) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassGStar_4jets", "; GStar Mass (with 4 small jets + ttbar) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassGStar_5jets", "; GStar Mass (with 5 small jets + ttbar) (GeV)", 36, 0, 3600);
  
  m_hSvc.create1D("seljetPt", "; Selected Jet p_{T} (GeV) ; Events", 20, 0, 2000);
  m_hSvc.create1D("PtFatJet_0", "; FatJet_0 p_{T} (GeV) ; Events", 20, 0, 2000);
  m_hSvc.create1D("PtFatJet_1", "; FatJet_1 p_{T} (GeV) ; Events", 20, 0, 2000);
   
  m_hSvc.create1D("PtZprime_FatJet", "; Zprime_FatJet p_{T} (GeV) ; Events", 20, 0, 2000);
  m_hSvc.create1D("PtZprime_2jets", "; Zprime_2jets p_{T} (GeV) ; Events", 20, 0, 2000);
  m_hSvc.create1D("PtGStar_1FatJet_2jets", "; GStar_1FatJet_2jets p_{T} (GeV) ; Events", 20, 0, 2000);
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
  
  m_hSvc.create1D("largeJetPt1", "; Large jet p_{T} (GeV); Events", 20, 0, 2000);
  m_hSvc.create1D("largeJetPt2", "; Large jet p_{T} (GeV); Events", 20, 0, 2000);
 
  
  m_hSvc.create1D( "Pt5jets_01", "; 5jets: Jet[0+1] p_{T} (GeV); Events", 20, 0, 2000);
  m_hSvc.create1D( "Pt5jets_234", "; 5jets: Jet[2+3+4] p_{T} (GeV); Events", 20, 0, 2000);
  m_hSvc.create1D( "Pt4jets_01", "; 4jets: Jet[0+1] p_{T} (GeV); Events", 20, 0, 2000);  
  m_hSvc.create1D( "Pt4jets_23", "; 4jets: Jet[2+3] p_{T} (GeV); Events", 20, 0, 2000);
  
  m_hSvc.create1D( "Pt3jets_01", "; 3jets: Jet[0+1] p_{T} (GeV); Events", 20, 0, 2000);
  m_hSvc.create1D( "Pt3jets_3", "; 3jets: Jet[3] p_{T} (GeV); Events", 20, 0, 2000);
   
  m_hSvc.create1D( "Pt1FatJet2jets_1Fat", "; 2jets_1FatJet: FatJet[0] p_{T} (GeV); Events", 20, 0, 2000);
  m_hSvc.create1D( "Pt1FatJet2jets_2jets", "; 2jets_1FatJet: Jet[1+2] p_{T} (GeV); Events", 20, 0, 2000);
  
  
  m_hSvc.create1D( "PtZprime2jets_W", "; 2jets: W p_{T} (GeV); Events", 20, 0, 2000);
  m_hSvc.create1D( "PtZprimeFatJet_W", "; 1FatJet: W p_{T} (GeV); Events", 20, 0, 2000);
  
  
  m_hSvc.create1D("Mass1FatJet_FatJet0", "; 1FatJet_2Closesjets: First Leading Large jet Mass (GeV); Events", 36, 0, 3600);
  m_hSvc.create1D("MassZprimeFatJet_FatJet0", "; ZPrime_1FatJet: M(FatJet[0]) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassZprime2jets_2jets", "; ZPrime_2SamllJets: M(SjetsFarLep[0]+SjetsFarLep[1]) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassGStar1FatJet2jets_Fatjet2jets", "; M(Jet_Fat[0]+SjetsCloseLep[0]+SjetsCloseLep[1]) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassZprime3jets_3Jets", "; M(SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassGStar4jets_4jets", "; M(SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassGStar5jets_5jets", "; M(SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+SjetsFarLep[4]) (GeV)", 36, 0, 3600);
  m_hSvc.create1D("MassGStar2FatJet_2FatJets", "; M(Jet_Fat[0]+Jet_Fat[1]) (GeV)", 36, 0, 3600);
  
  
  
  
  
  //New Mass
 m_hSvc.create1D("MGStar_true_reco_2Sjets", "; 2Sjets: Truth_Reco G* Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MT_true_reco_2Sjets", "; 2Sjets: Truth_Reco T' Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MTop_true_reco_2Sjets", "; 2Sjets: Truth_Reco Leptonic Top Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MW_true_reco_2Sjets", "; 2Sjets: Truth_Reco leptonic W Mass Comparisons; Events", 72, -3600, 3600);
 
 m_hSvc.create1D("MGStar_true_reco_1FatJet_2Sjets", "; 1FatJet2Sjets: Truth_Reco G* Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MT_true_reco_1FatJet_2Sjets", "; 1FatJet2Sjets: Truth_Reco T' Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MZ_true_reco_1FatJet_2Sjets", "; 1FatJet2Sjets: Truth_Reco  Z Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MW_true_reco_1FatJet_2Sjets", "; 1FatJet2Sjets: Truth_Reco leptonic W Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MH_true_reco_1FatJet_2Sjets", "; 1FatJet2Sjets: Truth_Reco H Comparisons; Events", 72, -3600, 3600);

 m_hSvc.create1D("MGStar_true_reco_3Sjets", "; 3Sjets: Truth_Reco G* Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MT_true_reco_3Sjets", "; 3Sjets: Truth_Reco T' Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MTop_true_reco_3Sjets", "; 3Sjets: Truth_Reco Leptonic Top Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MW_true_reco_3Sjets", "; 3Sjets: Truth_Reco leptonic W Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MH_true_reco_3Sjets", "; 3Sjets: Truth_Reco  H Mass Comparisons; Events", 72, -3600, 3600);
 
 
 m_hSvc.create1D("MGStar_true_reco_4Sjets", "; 4Sjets: Truth_Reco G* Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MT_true_reco_4Sjets", "; 4Sjets: Truth_Reco T' Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MTop_true_reco_4Sjets", "; 4Sjets: Truth_Reco Leptonic Top Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MW_true_reco_4Sjets", "; 4Sjets: Truth_Reco leptonic W Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MZ_true_reco_4Sjets", "; 4Sjets: Truth_Reco Z Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MH_true_reco_4Sjets", "; 4Sjets: Truth_Reco H Mass Comparisons; Events", 72, -3600, 3600);
 
 m_hSvc.create1D("MGStar_true_reco_5Sjets", "; 5Sjets: Truth_Reco G* Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MT_true_reco_5Sjets", "; 5Sjets: Truth_Reco T' Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MTop_true_reco_5Sjets", "; 5Sjets: Truth_Reco Leptonic Top Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MW_true_reco_5Sjets", "; 5Sjets: Truth_Reco leptonic W Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MZ_true_reco_5Sjets", "; 5Sjets: Truth_Reco Z Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MH_true_reco_5Sjets", "; 5Sjets: Truth_Reco H Mass Comparisons; Events", 72, -3600, 3600);
 
 m_hSvc.create1D("MGStar_true_reco_2FatJet", "; 2FatJet: Truth_Reco G* Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MT_true_reco_2FatJet", "; 2FatJet: Truth_Reco T' Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MTop_true_reco_2FatJet", "; 2FatJet: Truth_Reco Leptonic Top Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MW_true_reco_2FatJet", "; 2FatJet: Truth_Reco leptonic W Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MZ_true_reco_2FatJet", "; 2FatJet: Truth_Reco Z Mass Comparisons; Events", 72, -3600, 3600);
 m_hSvc.create1D("MH_true_reco_2FatJet", "; 2FatJet: Truth_Reco H Mass Comparisons; Events", 72, -3600, 3600);
  
 m_hSvc.create1D("MGStar_true", "; 2FatJet: Truth G* Mass ; Events", 36, 0, 3600);
 m_hSvc.create1D("MT_true", "; 2FatJet: Truth VLQ Top Mass ; Events", 36, 0, 3600);



 //New Pt
 m_hSvc.create1D("PtGStar_true_reco_2Sjets", "; 2Sjets: Truth_Reco G* Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtT_true_reco_2Sjets", "; 2Sjets: Truth_Reco T' Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtTop_true_reco_2Sjets", "; 2Sjets: Truth_Reco Leptonic Top Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtW_true_reco_2Sjets", "; 2Sjets: Truth_Reco leptonic W Pt Comparisons; Events", 20, 0, 2000);
 
 m_hSvc.create1D("PtGStar_true_reco_1FatJet_2Sjets", "; 1FatJet2Sjets: Truth_Reco G* Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtT_true_reco_1FatJet_2Sjets", "; 1FatJet2Sjets: Truth_Reco T' Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtZ_true_reco_1FatJet_2Sjets", "; 1FatJet2Sjets: Truth_Reco  Z Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtW_true_reco_1FatJet_2Sjets", "; 1FatJet2Sjets: Truth_Reco leptonic W Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtH_true_reco_1FatJet_2Sjets", "; 1FatJet2Sjets: Truth_Reco H Comparisons; Events", 20, 0, 2000);

 m_hSvc.create1D("PtGStar_true_reco_3Sjets", "; 3Sjets: Truth_Reco G* Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtT_true_reco_3Sjets", "; 3Sjets: Truth_Reco T' Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtTop_true_reco_3Sjets", "; 3Sjets: Truth_Reco Leptonic Top Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtW_true_reco_3Sjets", "; 3Sjets: Truth_Reco leptonic W Pt Comparisons; Events", 20, 0, 2000);
 
 
 m_hSvc.create1D("PtGStar_true_reco_4Sjets", "; 4Sjets: Truth_Reco G* Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtT_true_reco_4Sjets", "; 4Sjets: Truth_Reco T' Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtTop_true_reco_4Sjets", "; 4Sjets: Truth_Reco Leptonic Top Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtW_true_reco_4Sjets", "; 4Sjets: Truth_Reco leptonic W Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtZ_true_reco_4Sjets", "; 4Sjets: Truth_Reco Z Pt Comparisons; Events", 20, 0, 2000);
 
 m_hSvc.create1D("PtGStar_true_reco_5Sjets", "; 5Sjets: Truth_Reco G* Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtT_true_reco_5Sjets", "; 5Sjets: Truth_Reco T' Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtTop_true_reco_5Sjets", "; 5Sjets: Truth_Reco Leptonic Top Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtW_true_reco_5Sjets", "; 5Sjets: Truth_Reco leptonic W Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtZ_true_reco_5Sjets", "; 5Sjets: Truth_Reco Z Pt Comparisons; Events", 20, 0, 2000);
 
 m_hSvc.create1D("PtGStar_true_reco_2FatJets", "; 2FatJets: Truth_Reco G* Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtT_true_reco_2FatJets", "; 2FatJets: Truth_Reco T' Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtTop_true_reco_2FatJets", "; 2FatJets: Truth_Reco Leptonic Top Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtW_true_reco_2FatJets", "; 2FatJets: Truth_Reco leptonic W Pt Comparisons; Events", 20, 0, 2000);
 m_hSvc.create1D("PtZ_true_reco_2FatJets", "; 2FatJets: Truth_Reco Z Pt Comparisons; Events", 20, 0, 2000);
 
 m_hSvc.create1D("SjetsFarLepDR", "; DR(Jet_Fat[0],SjetsFarLep); Events", 100, -1, 9);
 
  
//PtTop_had
m_hSvc.create1D("PtTop_had_1FatJet", "; 1FatJet: Pt of the hadronic Top; Events", 20, 0, 2000);
m_hSvc.create1D("PtTop_had_1FatJet2jets", "; 1FatJet2jets: Pt of the hadronic Top; Events", 20, 0, 2000);
m_hSvc.create1D("PtTop_had_4jets", "; 4jets: Pt of the hadronic Top; Events", 20, 0, 2000);
m_hSvc.create1D("PtTop_had_5jets", "; 5jets: Pt of the hadronic Top; Events", 20, 0, 2000);
 m_hSvc.create1D("PtTop_had_2FatJet", "; 2FatJet: Pt of the hadronic Top; Events", 20, 0, 2000);
//PtSMTop_lep
 m_hSvc.create1D("PtSMTop_lep", "; Pt of the leptonic Top; Events", 20, 0, 2000);

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
  int sel_jet = -1;
  int SjetsFarL = -1;
  int SjetsFarLepLJ = -1;
 // int SjetsCloseL = -1;
  double dr = 1.5;
  double lpt =0.;
  int btags =0; 
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
 
 int fatjets = 0;
 int jetsmall = 0; 
 int SmaljetsFarLep = 0;
 //int SmaljetsCloseLep = 0;
 int SmaljetsFarLepLJet = 0;

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


  //vector<TLorentzVector> SjetsCloseLep;
  vector<TLorentzVector> SjetsFarLep;
  vector<TLorentzVector> SjetsFarLepLJet;
  
  for (int k = 0; k < e.jet().size(); ++k) {
  	jetsmall++;
   	double tdr = e.jet()[k].mom().DeltaR(l);
   	if(tdr < dr && e.jet()[k].mom().Perp() > lpt) {
  		    	sel_jet = k;
  		    	lpt = e.jet()[k].mom().Perp();
  		    	}
  }
  if (sel_jet == -1) return;
  
    if(m_electron){
    	 Mej = (e.electron()[0].mom()+e.jet()[sel_jet].mom()).M();
    	 //Menu = (e.electron()[0].mom()+e.momNu).M();
    }
   else{
    	 Mmj = (e.muon()[0].mom() + e.jet()[sel_jet].mom()).M();
    	 //Mmnu = (e.muon()[0].mom() +e.momNu).M();
    }
    MWj = momW.M()+e.jet()[sel_jet].mom().M();
   /* std::cout<<"MassW=momW.M()= "<<momW.M()<<std::endl;
    std::cout<<""<<std::endl;*/
    h->h1D("MassMuj","",s)->Fill(Mmj*1e-3, weight);
    h->h1D("MassElecj","",s)->Fill(Mej*1e-3, weight);
    h->h1D("MET","",s)->Fill(momNu.Perp(), weight);
    if(m_electron) {
    	h->h1D("MassSMTop","",s)->Fill((e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).M()*1e-3, weight);
    	h->h1D("PtSMTop_lep","",s)->Fill((e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).Perp()*1e-3, weight);
    }
    else {
    	h->h1D("MassSMTop","",s)->Fill((e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).M()*1e-3, weight);
    	h->h1D("PtSMTop_lep","",s)->Fill((e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
    }
    h->h1D("MassElecjNu","",s)->Fill(Mej*1e-3, weight);
    h->h1D("MassWj", "", s)->Fill(MWj*1e-3, weight);
    
//   for (int k = 0; k < e.jet().size(); ++k) {
//      jetsmall++;
//     
//    }
  	     
  TLorentzVector ljet;
   for (int k = 0; k < e.largeJet().size(); k++) {
  	bool passdPhi=false;
    bool passdR = false;
    bool passM = false;
    bool passD12 = false;
    bool passEta = false;
    bool passPt = false;
   
    if (e.largeJet()[k].passLoose()) {
    	h->h1D("largeJetPt1", "", s)->Fill(e.largeJet()[k].mom().Perp()*1e-3, weight);
       passdR = e.jet()[sel_jet].mom().DeltaR(e.largeJet()[k].mom()) > 1.5;
       passdPhi =true;
//      if(m_electron) passdPhi=e.largeJet()[k].mom().DeltaPhi(e.electron()[0].mom()) > 2.3;

      passM = e.largeJet()[k].mom().M() > 60e3;
      passEta = std::fabs(e.largeJet()[k].mom().Eta()) <  2.0;
      passPt = e.largeJet()[k].mom().Perp() > 200e3;

      passD12 = e.largeJet()[k].split12() > 20e3;
      if (passdR && passdPhi && passM && passD12 && passEta && passPt) {
        fatjets++;
        idx = k;
        h->h1D("largeJetPt2", "", s)->Fill(e.largeJet()[k].mom().Perp()*1e-3, weight);
        ljet = e.largeJet()[k].mom();
      }
    }
  }
 //if (fatjets < 1) return;   
// if (SmaljetsFarLep < 4 && fatjets < 1) return;
 if (jetsmall < 3 && fatjets < 1) return;
 
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
      passdPhi =true;
      passM = e.largeJet()[k].mom().M() > 60e3;
      passD12 = e.largeJet()[k].split12() > 20e3;
      passEta = std::fabs(e.largeJet()[k].mom().Eta()) <  2.0;
      passPt = e.largeJet()[k].mom().Perp() > 200e3;
      if (passdR && passdPhi && passM && passD12 && passEta && passPt) {
        Jet_Fat.push_back(e.largeJet()[k].mom());
      }
    }
    h->h1D("PtFatJet", "", s)->Fill(e.largeJet()[k].mom().Perp(), weight);
  }
  
  for (int k = 0; k < e.jet().size(); ++k) {
  	double tdr = e.jet()[k].mom().DeltaR(l);
  	double tdr2 = e.jet()[k].mom().DeltaR(ljet);
    if(tdr > dr){
    	SjetsFarL = k;
    	if (SjetsFarL == sel_jet) continue;
        SmaljetsFarLep++;
        h->h1D("SjetsFarLepDR", "", s)->Fill(e.jet()[k].mom().DeltaR(ljet), weight);
     }
   if(tdr > dr && tdr2 > 1.5){
   	    SjetsFarLepLJ = k; 
   	    //if (SjetsFarLepLJ == SjetsFarL) continue; 
   	    SmaljetsFarLepLJet++;
   	}  
   
  } 
  for (int k = 0; k < e.jet().size(); ++k) {
  	/*deltaR(lepton, SmalljetsCloseLep) > 1.5 and deltaR(SmalljetsCloseLep, FatJet) > 1.0 */
  	double tdr = e.jet()[k].mom().DeltaR(l);
  	double tdr2 = e.jet()[k].mom().DeltaR(ljet);
   	//if(tdr < dr) SjetsCloseLep.push_back(e.jet()[k].mom());   	
    if(tdr > dr){
    	 SjetsFarLep.push_back(e.jet()[k].mom());
    }
   if(tdr > dr && tdr2 > 1.5){SjetsFarLepLJet.push_back(e.jet()[k].mom());}
   
  }
   


  const TLorentzVector &selj = e.jet()[sel_jet].mom();
  const TLorentzVector &j = e.jet()[0].mom();
 
  h->h1D("lepPt", "", s)->Fill(l.Perp()*1e-3, weight);
  h->h1D("jetPt", "", s)->Fill(j.Perp()*1e-3, weight);
  h->h1D("seljetPt", "", s)->Fill(selj.Perp()*1e-3, weight);
  h->h1D("NSmallJets", "", s)->Fill(jetsmall, weight);
  
  h->h1D("NSmaljetsFarLepLJet", "", s)->Fill(SmaljetsFarLepLJet, weight);
  //h->h1D("NSmaljetsCloseLep", "", s)->Fill(SmaljetsCloseLep, weight);
  h->h1D("NSmaljetsFarLep", "", s)->Fill(SmaljetsFarLep, weight);
  if (fatjets > 0) {
  	const TLorentzVector &lj = e.largeJet()[idx].mom();
   	h->h1D("largeJetPt_0", "", s)->Fill(lj.Perp()*1e-3, weight);
    h->h1D("largeJetM_0", "", s)->Fill(lj.M()*1e-3, weight); 
    h->h1D("NLargeJets", "", s)->Fill(fatjets, weight);
   }

  if (fatjets > 1) {
    h->h1D("largeJetM_1", "", s)->Fill(Jet_Fat[1].M()*1e-3, weight);
    h->h1D("largeJetPt_1", "", s)->Fill(Jet_Fat[1].Perp()*1e-3, weight); 
  }           
                
   if(SmaljetsFarLep >= 5){
   	mj0j1j2j3j4 = (SjetsFarLep[0] + SjetsFarLep[1]+ SjetsFarLep[2]+ SjetsFarLep[3]+ SjetsFarLep[4]).M(); //mass of T
   	 h->h1D("MassTprime", "", s)->Fill(mj0j1j2j3j4*1e-3, weight);
   }
   
  
  for (int k = 0; k < e.jet().size(); ++k) {
  	if (e.jet()[k].pass() && e.jet()[k].btag()) {
      btags++;
      }
   }
  if (btags < 1) return;
  
  vector<TLorentzVector> SjetsFarLep_b;
  for (int k = 0; k < e.jet().size(); ++k) {
  	if (e.jet()[k].pass() && e.jet()[k].btag()) SjetsFarLep_b.push_back(e.jet()[k].mom());
  	}
  float MWb = (SjetsFarLep_b[0] + momW).M(); //mass of T
  if(SmaljetsFarLep < 5 && SmaljetsFarLep > 2) {
  	   mj0j1j2 = (SjetsFarLep[0] + SjetsFarLep[1]+ SjetsFarLep[2]).M(); //mass of T
   	//std::cout<<"jetsmall<5: n_jetSamll= "<<jetsmall<<std::endl;
   	h->h1D("MassTprime", "", s)->Fill((mj0j1j2+MWb)*1e-3, weight);
   }

  h->h1D("MassSMTop2","",s)->Fill(MWb*1e-3, weight);  
  h->h1D("MassWb", "", s)->Fill(MWb*1e-3, weight);
  // if (fatjets > 0) h->h1D("MassFatJet_0", "", s)->Fill(Jet_Fat[0].M()*1e-3, weight);
  h->h1D("NbJets", "", s)->Fill(btags, weight);
   //if (fatjets > 0) h->h1D("PtFatJet_0", "", s)->Fill(Jet_Fat[0].Perp()*1e-3, weight);


   if(m_electron){
//    	 Mej = (e.electron()[0].mom()+e.jet()[sel_jet].mom()).M();
//    	 h->h1D("MassZprime", "", s)->Fill((Jet_Fat[0]+SjetsFarLep[0]+e.electron()[0].mom()+momNu).M()*1e-3, weight);
    	 //Menu = (e.electron()[0].mom()+e.momNu).M();
    }
   else{
    	    Mmj = (e.muon()[0].mom() + e.jet()[sel_jet].mom()).M(); 
    	    //if (fatjets > 0) 
    	    if (fatjets == 1) {	
    	    	h->h1D("MassZprime_FatJet", "", s)->Fill((Jet_Fat[0]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).M()*1e-3, weight);
                h->h1D("PtZprime_FatJet", "", s)->Fill((Jet_Fat[0]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
                
                h->h1D("PtZprimeFatJet_W", "", s)->Fill((e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
                h->h1D("MassZprimeFatJet_FatJet0", "", s)->Fill((Jet_Fat[0]).M()*1e-3, weight);
                h->h1D("PtTop_had_1FatJet", "", s)->Fill((Jet_Fat[0]).Perp()*1e-3, weight);
    	    }   	    
    	    if(fatjets == 0 && SmaljetsFarLep == 2){
    	    	  	h->h1D("MassZprime_2jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).M()*1e-3, weight);
    	    	    h->h1D("PtZprime_2jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
    	    	    
    	    	    h->h1D("Pt2jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]).Perp()*1e-3, weight);
    	    	    h->h1D("PtZprime2jets_W", "", s)->Fill((e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
    	    	    h->h1D("MassZprime2jets_2jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]).M()*1e-3, weight);
 	 		      
    	   		    	   
    	    	  }
    	    	  
    	    	  if(fatjets == 1 && SmaljetsFarLepLJet==2){
    	    	 	 h->h1D("MassGStar_1FatJet_2jets", "", s)->Fill((Jet_Fat[0]+SjetsFarLepLJet[0]+SjetsFarLepLJet[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).M()*1e-3, weight);
    	    	 	 h->h1D("PtGStar_1FatJet_2jets", "", s)->Fill((Jet_Fat[0]+SjetsFarLepLJet[0]+SjetsFarLepLJet[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
    	    	     h->h1D("Pt1FatJet_2jets", "", s)->Fill((Jet_Fat[0]+SjetsFarLepLJet[0]+SjetsFarLepLJet[1]).Perp()*1e-3, weight);
    	    	     h->h1D("Pt1FatJet2jets_1Fat", "", s)->Fill((Jet_Fat[0]).Perp()*1e-3, weight);
    	    	     h->h1D("Pt1FatJet2jets_2jets", "", s)->Fill((SjetsFarLepLJet[0]+SjetsFarLepLJet[1]).Perp()*1e-3, weight);
    	    	     
    	    	     h->h1D("Mass1FatJet_FatJet0", "", s)->Fill(Jet_Fat[0].M()*1e-3, weight);
    	    	     h->h1D("MassGStar1FatJet2jets_Fatjet2jets", "", s)->Fill((Jet_Fat[0]+SjetsFarLepLJet[0]+SjetsFarLepLJet[1]).M()*1e-3, weight);
    	    	     h->h1D("PtTop_had_1FatJet2jets", "", s)->Fill((Jet_Fat[0]).Perp()*1e-3, weight); 
    	    	        	   	     
                 }


//                 if(fatjets == 1 && SmaljetsCloseLep ==2){
//    	    	 	 h->h1D("MassGStar_1FatJet_2jets", "", s)->Fill((Jet_Fat[0]+SjetsCloseLep[0]+SjetsCloseLep[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).M()*1e-3, weight);
//    	    	 	 h->h1D("PtGStar_1FatJet_2jets", "", s)->Fill((Jet_Fat[0]+SjetsCloseLep[0]+SjetsCloseLep[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
//    	    	     h->h1D("Pt1FatJet_2jets", "", s)->Fill((Jet_Fat[0]+SjetsCloseLep[0]+SjetsCloseLep[1]).Perp()*1e-3, weight);
//    	    	     h->h1D("Pt1FatJet2jets_1Fat", "", s)->Fill((Jet_Fat[0]).Perp()*1e-3, weight);
//    	    	     h->h1D("Pt1FatJet2jets_2jets", "", s)->Fill((SjetsCloseLep[0]+SjetsCloseLep[1]).Perp()*1e-3, weight);
//    	    	     
//    	    	     h->h1D("Mass1FatJet_FatJet0", "", s)->Fill(Jet_Fat[0].M()*1e-3, weight);
//    	    	     h->h1D("MassGStar1FatJet2jets_Fatjet2jets", "", s)->Fill((Jet_Fat[0]+SjetsCloseLep[0]+SjetsCloseLep[1]).M()*1e-3, weight); 	    	   	     
//                 }
    	 	if (SmaljetsFarLep > 2 && fatjets == 0){
    	 		if (SmaljetsFarLep == 3){
    	 			 h->h1D("MassZprime_3jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).M()*1e-3, weight);
    	 			 h->h1D("PtZprime_3jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
    	 		     h->h1D("Pt3jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[3]).Perp()*1e-3, weight);
    	 		     
    	 		     h->h1D("Pt3jets_01", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]).Perp()*1e-3, weight);
    	 		     h->h1D("Pt3jets_3", "", s)->Fill((SjetsFarLep[3]).Perp()*1e-3, weight);
    	 		     h->h1D("MassZprime3jets_3Jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]).M()*1e-3, weight);
    	 		 
    	   			 
                          }
    	 		 if (SmaljetsFarLep == 4){
    	 		 	 h->h1D("MassGStar_4jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).M()*1e-3, weight);
    	 		 	 h->h1D("PtGStar_4jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
    	 		     h->h1D("Pt4jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]).Perp()*1e-3, weight);
    	 		     
    	 		     h->h1D("Pt4jets_01", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]).Perp()*1e-3, weight);  
    	 		    // h->h1D("Pt4jets_1", "", s)->Fill((SjetsFarLep[1]).Perp()*1e-3, weight);  
    	 		     h->h1D("Pt4jets_23", "", s)->Fill((SjetsFarLep[2]+SjetsFarLep[3]).Perp()*1e-3, weight);
                     h->h1D("MassGStar4jets_4jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]).M()*1e-3, weight);
                     h->h1D("PtTop_had_4jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]).Perp()*1e-3, weight);  
                    	   			                         
                     }
    	 		 if (SmaljetsFarLep > 4) {
    	 		 	h->h1D("MassGStar_5jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+SjetsFarLep[4]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).M()*1e-3, weight);
    	 		 	h->h1D("PtGStar_5jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+SjetsFarLep[4]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
    	 		    h->h1D("Pt5jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+SjetsFarLep[4]).Perp()*1e-3, weight);
    	 		    
    	 		    h->h1D("Pt5jets_01", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]).Perp()*1e-3, weight);
    	 		    //h->h1D("Pt5jets_1", "", s)->Fill((SjetsFarLep[1]).Perp()*1e-3, weight);
    	 		    h->h1D("Pt5jets_234", "", s)->Fill((SjetsFarLep[2]+SjetsFarLep[3]+SjetsFarLep[4]).Perp()*1e-3, weight);
    	 		    h->h1D("MassGStar5jets_5jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+SjetsFarLep[4]).M()*1e-3, weight);
    	 		    h->h1D("PtTop_had_5jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]).Perp()*1e-3, weight);
    	 		  	
                         }
    	 	}
    	    if(fatjets > 1) {
    	    	h->h1D("MassGStar_2FatJet", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).M()*1e-3, weight);
    	    	h->h1D("PtGStar_2FatJet", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
    	        h->h1D("Pt2FatJet", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]).Perp()*1e-3, weight);
    	        h->h1D("MassFatJet_0", "", s)->Fill(Jet_Fat[0].M()*1e-3, weight);
                h->h1D("PtFatJet_0", "", s)->Fill(Jet_Fat[0].Perp()*1e-3, weight);
                h->h1D("MassFatJet_1", "", s)->Fill(Jet_Fat[1].M()*1e-3, weight);
                h->h1D("PtFatJet_1", "", s)->Fill(Jet_Fat[1].Perp()*1e-3, weight);
                h->h1D("PtTop_had_2FatJet", "", s)->Fill(Jet_Fat[1].Perp()*1e-3, weight); 
                
                h->h1D("MassGStar2FatJet_2FatJets", "", s)->Fill((Jet_Fat[0]+Jet_Fat[1]).M()*1e-3, weight); 
               

    	    }
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
 
