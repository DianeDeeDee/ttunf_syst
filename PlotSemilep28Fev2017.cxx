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
bool sortFunction2(const Jet &a, const Jet &b){
	return a.mom().Perp() > b.mom().Perp();
	}
//className::MethodeName ->les Methodes sont run () and getWFromLeptonicDecay()
//Main function is in runsd.cxx
PlotSemilep::PlotSemilep(const std::string &filename, bool electron, const std::vector<std::string> &systs)
: Plot(filename, systs), m_electron(electron) { //Utilisation du constructeur et son implementation
	
	m_hSvc.create1D("jetPt", "; Jet p_{T} (GeV) ; Events", 20,25,500);//20, 0, 800); //Faire: 20,25,500
	m_hSvc.create1D("MET", "; Missing ET (GeV); Events",10,-20,200);//41, -20, 800); //Faire: 10,-20,200
	m_hSvc.create1D("seljetPt", "; Selected Jet p_{T} (GeV) ; Events",10,0.,300);// 20, 0, 1000); //Faire 10,0.,300
	m_hSvc.create1D("lepPt", " Lepton p_{T} (GeV) ; Events",15,0, 250);// 15, 0, 500); //Faire 15,0, 250
	m_hSvc.create1D("NbJets", "; $b$-tagged jet multiplicity; Events ",6,0.5,6.5);// 20,0,20);//Faire 6,0.5,6.5
	m_hSvc.create1D("NSmallJets", "; Small-$R$ jets multiplicity ; Events",8,2.5,10.5);// 20,0,20);  //Faire 8,2.5,10.5
	m_hSvc.create1D("NLargeJets", "; Large-$R$ jet multiplicity ; Events", 4,-0.5,3.5);//11, -0.5, 10.5);//Faire 4,-0.5,3.5
	m_hSvc.create1D("NSelJets", "; Selected jets multiplicity; Events", 8,-0.5,7.5);//20,0,20);//Faire: 8,-0.5,7.5
	//0 Fat jet
	m_hSvc.create1D("MassGStar_3SjFarNSeljL", "; GStar Mass (3Sjets + W) (GeV); Events", 12,200,2600);//10,200,1700);//10, 200, 1700);//faire 10,0,1500
	m_hSvc.create1D("MassT_3SjFar", "; T Mass (3Sjets) (GeV); Events",10,0,1000);// 20, 0, 3600);////Faire 10,0,1000
	m_hSvc.create1D("MassGStar_5SjFarNSeljL", "; GStar Mass (5Sjets + W) (GeV); Events", 12,200,2600);//12,300,2100);//20, 0, 3600);//Faire 12,300,2100
	m_hSvc.create1D("MassT_5SjFarNSeljL", "; T Mass (5Sjets) (GeV); Events", 15,0,1500);//20, 0, 3600);//Faire 15,0,1500
	m_hSvc.create1D("MassGStar_4jets", "; GStar Mass (4Sjets + W) (GeV); Events" ,12,200,2600);//12,200,2000);//20, 0, 3600);//Faire 12,200,2000
	m_hSvc.create1D("MassT_4jets", "; T Mass (4Sjets) (GeV); Events" ,12,0,1200);//20, 0, 3600);//Faire 12,0,1200
	
	//1 Fat jet
	m_hSvc.create1D("MassGStar_1Fatj2SjFarLepLJetNSeljL", "; GStar Mass (2Sjets + 1Fat jets + W) (GeV)",12,200,2600);//12, 500, 2600);// 20, 0, 3600); //Faire 12, 500, 2600
	m_hSvc.create1D("MassT_1Fatj2SjFarLepLJetNSeljL", " ; T Mass (2Sjets + 1Fat jets) (GeV)",20,0,1100);// 20, 0, 3600);//Faire 12,0,1100	
	m_hSvc.create1D("PtGStar_3SjFarNSeljL", "; GStar p_{T} (3Sjets + W) (GeV) ; Events", 20, 0, 2000);
	m_hSvc.create1D("PtGStar_4SjFarNSeljL", "; GStar p_{T} (4Sjets + W) (GeV) ; Events", 20, 0, 2000);
	m_hSvc.create1D("PtGStar_5SjFarNSeljL", "; GStar p_{T} (5Sjets + W) (GeV) ; Events", 20, 0, 2000);
	//  m_hSvc.create1D("MassGStar_4SjFarNSeljL", "; GStar Mass (4Sjets + W) (GeV)" ,36, 0, 3600); 
	m_hSvc.create1D("Pt_3SjFar", "; 3Farjets p_{T} (GeV) ; Events", 20, 0, 2000);
	m_hSvc.create1D("Pt_4SjFar", "; 4Farjets p_{T} (GeV) ; Events", 20, 0, 2000);
	m_hSvc.create1D("Pt_5SjFar", "; 5Farjets p_{T} (GeV) ; Events", 20, 0, 2000);
	m_hSvc.create1D("D12", "Large Jet D12 (GeV); Events",30, 0., 300.); //80, 0 800 
	m_hSvc.create1D("D12_0", "Large Jet D12 (GeV)",30, 0., 300.); //80, 0 800 
	m_hSvc.create1D("D12_00", "Large Jet D12 (GeV)",30, 0., 300.); //80, 0 800
	m_hSvc.create1D("D120", "Large Jet D12 (GeV)",30, 0., 300.); //80, 0 800
	m_hSvc.create1D("Eta1", "Eta(Lepton) at low dR", 10, -3., 3);
	m_hSvc.create1D("lepPt1", " Lepton p_{T} (GeV) at low dR; Lepton Pt (GeV)", 20, 0, 1000);
	m_hSvc.create1D("seljetPt1", " Selected Jet p_{T} (GeV) at low dR; Selected Jet p_{T} (GeV)", 20, 0, 1000);
	m_hSvc.create1D("lepPt2", " Lepton p_{T} (GeV) when dR < 0.4; Pt(dR<0.4) (GeV)", 10, 0, 500); //old: 30 0 3000
	m_hSvc.create1D("Eta", "Eta(Lepton)", 10, -3., 3); //done
	m_hSvc.create1D("NSmallJets0", "; Size of Small jets ; Events", 20,0,20);
	m_hSvc.create1D("NSmallJets1", "; Number of Small jets when Nlarge jet>0 ; Events", 21,-0.5,20.5);
	m_hSvc.create1D("NSmaljetsFarLepLJet","; Number of Small jets far from the lepton and large jet; Events", 20,0,20);
	m_hSvc.create1D("NSmaljetsFarLep", "; Number of Small jets far from the lepton ; Events", 20,0,20);   
	m_hSvc.create1D("NSmaljetsCloseLep", "; Number of Small jets close to the lepton ; Events", 20,0,20);
	m_hSvc.create1D("NLargeJets_0", "; Number of Fat jets ; N_{Fat jets}", 11, -0.5, 10.5);
	m_hSvc.create1D("NLargeJets_00", "; Number of Fat jets ; N_{Fat jets}", 11, -0.5, 10.5);
	m_hSvc.create1D("NSmallJets_chan_FatJet0", "; Number of small jets -0 Fat Jets ; Events", 20, 0, 20);
	m_hSvc.create1D("NSjetsFarLepLJet_chan_FatJet0", "; Number of small jets far from to lept-0 Fat Jets ; Events", 20, 0, 20);
	m_hSvc.create1D("NSmaljetsFarLep_chan_FatJet0", "; Number of small jets far from lept-0 Fat Jets ; Events", 20, 0, 20);
	m_hSvc.create1D("NSmallJets_FatJet1", "; Number of small jets -1 Fat Jets ; Events", 20, 0, 20);
	m_hSvc.create1D("NSjetsFarLepLJet_FatJet1", "; Number of small jets far from lept and large jet -1 Fat Jets ; Events", 20, 0, 20);
	m_hSvc.create1D("NSmaljetsFarLep_chan_FatJet1", "; Number of small jets far from lept-1 Fat Jets ; Events", 20, 0, 20);
	m_hSvc.create1D("NSmallJets_FatJet2", "; Number of small jets -2 Fat Jets ; Events", 20, 0, 20);
	m_hSvc.create1D("NSmaljetsFarLepLJet_FatJet2", "; Number of small jets far from lept and large jet-2 Fat Jets ; Events", 20, 0, 20);
	m_hSvc.create1D("NSmaljetsFarLep_chan_FatJet2", "; Number of small jets far from lept-2 Fat Jets ; Events", 20, 0, 20);
	m_hSvc.create1D("MET_plus_mtW", "; MET+mtW ; Events", 20, 0, 100);
	m_hSvc.create1D("mtW", "; mtW ; Events", 20, 0, 100);
	//1 Fat jet
	m_hSvc.create1D("MassZprime_2jets", "; Zprime_2jets Mass (GeV) ; Events", 20, 0, 2000);
	m_hSvc.create1D("PtZprime_2jets", "; Zprime_2jets p_{T} (GeV) ; Events", 20, 0, 2000);
	m_hSvc.create1D("Pt2jets", "; 2jets p_{T} (GeV) ; Events", 20, 0, 2000);
	m_hSvc.create1D("Pt_1Fatj2SjFarLepLJet", " 1FatJet_2jetsCloseToLep p_{T} (GeV) ; Events", 20, 0, 2000);
	m_hSvc.create1D("PtGStar_1Fatj2SjFarLepLJetNSeljL", "; GStar p_{T} (2Sjets + 1Fat jets + W) (GeV) ; Events", 20, 0, 1000);
	m_hSvc.create1D("MassZprimeFatJet_FatJet0", " ZPrime_1FatJet: M(FatJet[0]) (GeV)", 20, 0, 3600);
	m_hSvc.create1D( "Pt_1Fatj", " 2Sjets_1FatJet: FatJet[0] p_{T} (GeV); Events", 20, 0, 2000);//done old: 20, 0, 2000
	m_hSvc.create1D( "Pt_SjFarLepLJet0", " 1FatJet_2jets: JetCloseLep[0] p_{T} (GeV); Events", 20, 0, 1000);//done old: 20, 0, 2000
	m_hSvc.create1D( "Pt_SjFarLepLJet1", " 1FatJet_2jets: JetCloseLep[1] p_{T} (GeV); Events", 20, 0, 1000);//done old: 20, 0, 2000
	m_hSvc.create1D("Mass_1Fatj", " 1FatJet_2jets: M(FatJet[0]) (GeV)", 20, 0, 1000); //done old: 36, 0, 3600
	m_hSvc.create1D("Eta_1Fatj", " 1FatJet_2jets: Eta(FatJet[0])", 10, -3., 3); //done
	m_hSvc.create1D("Eta_SjFarLepLJet0", " 1FatJet_2jets: Eta(JetCloseLep[0])", 10, -3., 3); //done?
	m_hSvc.create1D("Eta_SjFarLepLJet1", " 1FatJet_2jets: Eta(JetCloseLep[1])", 10, -3., 3); //done?
	m_hSvc.create1D( "PtZprimeFatJet_W", " 1FatJet: W p_{T} (GeV); Events", 20, 0, 2000);  
	m_hSvc.create1D( "Pt_2SjFarLepLJet", " 2jets_1FatJet: Jet[1+2] p_{T} (GeV); Events", 20, 0, 2000);
	m_hSvc.create1D("PtLj", "Fat jet Pt  (GeV)", 32, 0, 800);
	m_hSvc.create1D("PtLj_00", "Fat jet Pt  (GeV)", 32, 0, 800);
	m_hSvc.create1D("Mass_00", "M (GeV)", 20, 0, 1000);
	m_hSvc.create1D("Dr_00", " dR (selj,Lj)", 20, 0., 10.);
	m_hSvc.create1D("DPhi_00", " dPhi (lepton,Lj)", 20, 0., 10.);
	m_hSvc.create1D("Eta_00", " Eta (Lj)", 80, -10., 10.);
	m_hSvc.create1D("PtLj_0", "Fat jet Pt  (GeV)", 32, 0, 800);
	m_hSvc.create1D("Mass_0", "M (GeV)", 20, 0, 1000);
	m_hSvc.create1D("Dr_0", " dR (selj,Lj)", 20, 0., 10.);
	m_hSvc.create1D("DPhi_0", " dPhi (lepton,Lj)", 20, 0., 10.);
	m_hSvc.create1D("Eta_0", " Eta (Lj)", 10, -5., 5.);
	m_hSvc.create1D("Mass", "M (GeV)", 20, 0, 1000);
	m_hSvc.create1D("Dr", " dR (selj,Lj)", 20, 0., 10.);
	m_hSvc.create1D("DPhi", " dPhi (lepton,Lj)", 20, 0., 10.);
	m_hSvc.create1D("Eta0", " Eta (Lj)", 80, -10., 10.);
	m_hSvc.create1D("largeJetM", "Large jet Mass (GeV)", 15, 0, 500);
	m_hSvc.create1D("Mtop", "Leptonic Top Mass (GeV)", 20, 0, 1000);
	m_hSvc.create1D("largeJetPt", "Fat jet Pt  (GeV)", 14,150,850);//20, 0, 800);//Faire: 14,150,850
	m_hSvc.create1D("largeJetPt_0", "Fat jet Pt  (GeV)", 20, 0, 800);
	m_hSvc.create1D("largeJetPt_00", "Fat jet Pt  (GeV)", 20, 0, 800);
	m_hSvc.create1D("largeJetPt0", "Fat jet Pt  (GeV)", 20, 0, 800);
	m_hSvc.create1D("Pttop", "Leptonic Top Pt (GeV)", 20, 0, 1000);
	//2 fat jets
	m_hSvc.create1D("Mass_2Fatj", " M(Jet_Fat[0]+Jet_Fat[1]) (GeV)", 200, 0, 3600);
	m_hSvc.create1D("Pt2FatJet_FatJet1", " FatJet_1 p_{T} (GeV) ; Events", 200, 0, 2000); //done
	m_hSvc.create1D("Eta2FatJet_FatJet0", " 2FatJets: Eta(FatJet[0])", 90, -3., 3); //done
	m_hSvc.create1D("Eta2FatJet_FatJet1", " 2FatJets: Eta(FatJet[1])", 90, -3., 3); //done
	m_hSvc.create1D("Pt2FatJet", " FatJets p_{T} (GeV) ; Events", 200, 0, 2000);
	m_hSvc.create1D("Mass2FatJet_FatJet0", " First Large jet M (GeV); Events", 200, 0, 800); //done old: 100, 0, 1000
	m_hSvc.create1D("PtGStar_2FatjNSeljL", " GStar p_{T} (2Fat jets + W) (GeV) ; Events", 200, 0, 2000);
	m_hSvc.create1D("Mass2FatJet_FatJet1", " Second Large jet M (GeV); Events", 200, 0, 800); //old: 100, 0, 1000
	m_hSvc.create1D("MassGStar_2FatjNSeljL", " GStar Mass (2Fat jets + W) (GeV)", 200, 0, 3600);
	m_hSvc.create1D("Pt2FatJet_FatJet0", " FatJet_0 p_{T} (GeV) ; Events", 200, 0, 1000); //done
	m_hSvc.create1D("dR_closejet", " dR (lepton,closestJet) ; dR (lepton,closestJet)", 20, 0., 10.); //done 
	m_hSvc.create1D("dR_seljet", " dR (lepton,SelJet) ; dR (lepton,SelJet)", 10, 0., 3.); //done
	m_hSvc.create1D("dR_bjet", " dR (lepton,bjet) ; dR (lepton,SelJet)", 20, 0., 10.); //done
	}

PlotSemilep::~PlotSemilep() { //le destructeur
/* Rien à mettre ici car on ne fait pas d'allocation dynamique
 *
 *     dans la classe Personnage. Le destructeur est donc inutile mais
 *
 *         je le mets pour montrer à quoi cela ressemble.
 *
 *             En temps normal, un destructeur fait souvent des delete et quelques
 *
 *                 autres vérifications si nécessaire avant la destruction de l'objet. */
	}

void PlotSemilep::run(const Event &e, double weight, double pweight, const std::string &s) {
	HistogramService *h = &m_hSvc;  //*h is the pointer variable; &m_hSvc store address of m_hSvc in pointer variable cout<<h==>adress stored in h; cout<<*h=>gives value of *h variable
	// h->h1D("RunNumber", "", s)->Fill(e.runNumber());  
	//std::cout << "P:148" << std::endl; 
	double Mmj = 0.;
	double Mej = 0.;
	double Mmjnu = 0.;
	double Mejnu = 0.;
	double MWj=0.;
	int SmaljetsFarLep = 0;
    int SmaljetsCloseLep = 0;
    int SjetsFarL = -1;
    int sel_jet = -1;
    int SjetsCloseL = -1;
	double dr = 1.5;
	double lpt = 0;


	int b_jet = -1;
	
	//double dr = 1.5;
	// double lpt =0.;
	int btags =0;  
	int idx = -1;
	float mj0j1j2j3j4 = 0.; //mass of T
	float mj0j1j2 = 0.;
	//For MET
	float METx = e.met().Px();//[Missing E_T x direction]; // in MeV
	float METy = e.met().Py();//[Missing E_T y direction]; // in MeV
	
	float mass = 0; // mas of the lepton
	double Menu = 0;
	double Mmnu = 0;
	
	
	int fatjets = 0; int fatjets1=0; int fatjets2=0;
	int jetsmall = 0; int jetsmall1=0; int jetsmall2=0;
	// apply extra cuts
	
	if (e.passReco()) {// events passes cuts from EventCutterReco
		TLorentzVector l;
		TLorentzVector momNu;
		TLorentzVector momW; // this will be returned by the function -> it is the W boson vector after the neutrino p_z is found
		TLorentzVector momLepton; float NuX; float NuY;float NuZ; float lepMass;
		if (m_electron) {
			l = e.electron()[0].mom();
			mass = 0.510998910; // mass of the electron in MeV
			} 
		else {
			l = e.muon()[0].mom();
			mass = 105.6583668; // mass of the muon in MeV
			} 
		getWFromLeptonicDecay(l, METx, METy, momNu, momW, mass);
		
		std::vector<LargeJet> sortedJets = e.largeJet(); //to sort the jets by mass
		std::sort(sortedJets.begin(), sortedJets.end(), sortFunction);
		
		vector<TLorentzVector> SjetsCloseLep;
		vector<TLorentzVector> Jet_Fat;
		TLorentzVector ljet;
		vector<TLorentzVector> SjetsFarLep;
		vector<TLorentzVector> SjetsFarLepLJet;
		int SjetsFarLepLJ = -1;
		int SmaljetsFarLepLJet =0;
		double mindr = 9999;
		for (int k = 0; k < e.jet().size(); ++k) { 
			double tdr = e.jet()[k].mom().DeltaR(l);
			if(tdr < dr && e.jet()[k].mom().Perp() > lpt) {
				sel_jet = k;
				lpt = e.jet()[k].mom().Perp();
				if (tdr < mindr) mindr = tdr;  
				}
			}
		
		if (sel_jet == -1) return;
		
		for (int k = 0; k < e.jet().size(); ++k) {
			jetsmall++;
			double tdr = e.jet()[k].mom().DeltaR(l);
			if(tdr < dr){
				SjetsCloseL = k;
				if (SjetsCloseL == sel_jet) continue; 
				SmaljetsCloseLep++;
				}
			else if(tdr > dr){
				SjetsFarL = k;
				//if (SjetsFarL == SjetsCloseL) continue;
				if (SjetsFarL == sel_jet) continue;
				SmaljetsFarLep++;
				}     
			}

		for (int k = 0; k < e.jet().size(); ++k) {
			double tdr = e.jet()[k].mom().DeltaR(l);
			if(tdr < dr) SjetsCloseLep.push_back(e.jet()[k].mom());   	
			else if(tdr > dr) SjetsFarLep.push_back(e.jet()[k].mom());
			}	
				
		float mtw = std::sqrt(2*l.Perp()*momNu.Perp()*(1 - std::cos(l.Phi() - momNu.Phi())) );
		// std::cout<<"momNu.Perp()= "<<momNu.Perp()<<"   momNu.Phi()="<<momNu.Phi()<<"     cos(l.Phi() - e.momNu.Phi())="<<std::cos(l.Phi() - e.momNu.Phi())<<std::endl;
		//std::cout<<"l.Perp()="<<l.Perp()<<"     l.Phi()="<<l.Phi()<<std::endl;
		if((momNu.Perp()<20e3)) std::cout<<"MET<20e3!!!  MET="<<(momNu.Perp())*1e-3<<std::endl;
		if((mtw+momNu.Perp())<60e3) std::cout<<"mtw+MET<60e3!!!  mtw+MET="<<(momNu.Perp()+mtw)*1e-3<<std::endl;
		
		for (int k = 0; k < e.largeJet().size(); k++) {
			bool passdPhi=false;
			bool passdR = false;
			bool passM = false;
			bool passD12 = false;
			bool passEta = false;
			bool passPt = false;
			/*h->h1D("PtLj_0", "", s)->Fill((e.largeJet()[k].mom().Perp())*1e-3, weight);
			h->h1D("Mass_0", "", s)->Fill((e.largeJet()[k].mom().M())*1e-3, weight);
			h->h1D("Dr_0", "", s)->Fill(e.jet()[sel_jet].mom().DeltaR(e.largeJet()[k].mom()), weight);
			h->h1D("DPhi_0", "", s)->Fill(e.largeJet()[k].mom().DeltaPhi(l), weight);
			h->h1D("Eta_0", "", s)->Fill(e.largeJet()[k].mom().Eta(),weight);
			*/
			if (e.largeJet()[k].passLoose()) {
				passdR = e.jet()[sel_jet].mom().DeltaR(e.largeJet()[k].mom()) > 1.5;
				passM = e.largeJet()[k].mom().M() > 100e3;//60e3;
				passEta = std::fabs(e.largeJet()[k].mom().Eta()) <  2.0;//true
				passPt = e.largeJet()[k].mom().Perp() > 300e3;//300e3;//300e3;//160e3;//true
				passdPhi = e.largeJet()[k].mom().DeltaPhi(l) > 2.3;//true
				passD12 = e.largeJet()[k].split12() > 40e3;
				if (passPt && passdR && passdPhi && passM  && passEta &&passD12) {
					idx = k;
					ljet = e.largeJet()[k].mom();
//					Jet_Fat.push_back(e.largeJet()[k].mom());
					fatjets++;
					}
				}
			}
		
/*		for (int k = 0; k < e.jet().size(); ++k) {
			double tdr = e.jet()[k].mom().DeltaR(l);
			double tdr2 = e.jet()[k].mom().DeltaR(ljet);
			if(tdr > 1.5){
				//	 std::cout<<"381 k="<<k<<std::endl;
				SjetsFarL = k;
				if (SjetsFarL == sel_jet) continue;
				// std::cout<<"384"<<std::endl;
				SjetsFarLep.push_back(e.jet()[k].mom());
				// std::cout<<"386"<<std::endl;
				SmaljetsFarLep++;
				//   std::cout<<"387"<<std::endl;
				//h->h1D("SjetsFarLepDR", "", s)->Fill(e.jet()[k].mom().DeltaR(ljet), weight);
				}
			if(tdr > 1.5 && tdr2 > 1.5){
				SjetsFarLepLJ = k;
				SjetsFarLepLJet.push_back(e.jet()[k].mom());
				//if (SjetsFarLepLJ == SjetsFarL) continue;
				SmaljetsFarLepLJet++;
				}
			if (e.jet()[k].pass()) jetsmall++;			
			}
		//std::cout<<"405 SmaljetsFarLep= "<<SmaljetsFarLep<<" SmaljetsFarLepLJet="<<SmaljetsFarLepLJet<<"  jetsmall="<<jetsmall<<std::endl;
		
		if (jetsmall <4  && fatjets < 1) return;
		
	*/
if (jetsmall < 4  && fatjets < 1) return;	   	
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
				passEta = std::fabs(e.largeJet()[k].mom().Eta()) <  2.0;//true
				passPt = e.largeJet()[k].mom().Perp() > 300e3;//300e3;//300e3;//160e3;//true
				passdPhi = e.largeJet()[k].mom().DeltaPhi(l) > 2.3;//true
				passD12 = e.largeJet()[k].split12() > 40e3;
				if (passPt && passdR && passdPhi && passM  && passEta &&passD12) {
					Jet_Fat.push_back(e.largeJet()[k].mom());
					}
				}
			}
			
//		if (jetsmall < 4  && fatjets < 1) return;	
		const TLorentzVector &selj = e.jet()[sel_jet].mom();
        const TLorentzVector &j = e.jet()[0].mom();		
		double mindr2 = 9999;
		
		for (int k = 0; k < e.jet().size(); ++k) {
			double tdr = e.jet()[k].mom().DeltaR(l);
			if (e.jet()[k].pass() && e.jet()[k].btag()) {
				btags++; b_jet = k;
				//if (tdr < mindr2) mindr2 = tdr;
				}
			}
		if (btags < 1) return;
		
		//if (fatjets>1) std::cout<<"fatjets= "<<fatjets<<std::endl;
		h->h1D("mtW", "", s)->Fill(mtw*1e-3, weight);
		h->h1D("MET_plus_mtW", "", s)->Fill((momNu.Perp()+mtw)*1e-3, weight);
		h->h1D("MET","",s)->Fill(momNu.Perp()*1e-3, weight);
		h->h1D("NSmallJets", "", s)->Fill(jetsmall, weight);
		h->h1D("NSmaljetsFarLepLJet", "", s)->Fill(SmaljetsFarLepLJet, weight);
		h->h1D("NSmaljetsFarLep", "", s)->Fill(SmaljetsFarLep, weight);
		h->h1D("NLargeJets", "", s)->Fill(fatjets, weight);
		h->h1D("NSelJets", "", s)->Fill(sel_jet, weight);		
		h->h1D("lepPt", "", s)->Fill(l.Perp()*1e-3, weight);
		h->h1D("jetPt", "", s)->Fill(j.Perp()*1e-3, weight);
		h->h1D("seljetPt", "", s)->Fill(selj.Perp()*1e-3, weight);
		h->h1D("NbJets", "", s)->Fill(btags, weight);
		h->h1D("dR_seljet", "", s)->Fill(e.jet()[sel_jet].mom().DeltaR(l), weight);
		h->h1D("dR_bjet", "", s)->Fill(e.jet()[b_jet].mom().DeltaR(l), weight);
		h->h1D("Eta", "", s)->Fill(l.Eta(), weight);  
		/////////////////////// Electron steam
		if(m_electron){
			h->h1D("Mtop", "", s)->Fill((e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).M()*1e-3, weight);
			h->h1D("Pttop", "", s)->Fill((e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).Perp()*1e-3, weight);
			
			if (SmaljetsFarLep > 2){
				h->h1D("MassGStar_3SjFarNSeljL", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).M()*1e-3, weight);
				h->h1D("PtGStar_3SjFarNSeljL", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).Perp()*1e-3, weight);
				h->h1D("Pt_3SjFar", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[3]).Perp()*1e-3, weight);
				h->h1D("MassT_3SjFar", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[3]).M()*1e-3, weight);
				if (SmaljetsFarLep > 3){
					//h->h1D("MassGStar_4SjFarNSeljL", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).M()*1e-3, weight);
					h->h1D("MassGStar_4jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).M()*1e-3, weight);
					h->h1D("MassT_4jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]).M()*1e-3, weight);
					
					h->h1D("PtGStar_4SjFarNSeljL", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).Perp()*1e-3, weight);
					h->h1D("Pt_4SjFar", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]).Perp()*1e-3, weight); 
					} 
				if (SmaljetsFarLep > 4){
					h->h1D("MassGStar_5SjFarNSeljL", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+SjetsFarLep[4]+e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).M()*1e-3, weight);
					h->h1D("MassT_5SjFarNSeljL", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+SjetsFarLep[4]).M()*1e-3, weight);
					h->h1D("PtGStar_5SjFarNSeljL", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+SjetsFarLep[4]+e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).Perp()*1e-3, weight);
					h->h1D("Pt_5SjFar", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+SjetsFarLep[4]).Perp()*1e-3, weight);
					}
				}
			if(fatjets > 0) {
				h->h1D("largeJetM", "", s)->Fill((Jet_Fat[0]).M()*1e-3, weight);
				h->h1D("largeJetPt", "", s)->Fill(Jet_Fat[0].Perp()*1e-3, weight);
				if(SmaljetsFarLep > 1){
					h->h1D("Pt2jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]).Perp()*1e-3, weight);
					h->h1D("MassGStar_1Fatj2SjFarLepLJetNSeljL", "", s)->Fill((Jet_Fat[0]+SjetsFarLep[0]+SjetsFarLep[1]+e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).M()*1e-3, weight);
					h->h1D("MassT_1Fatj2SjFarLepLJetNSeljL", "", s)->Fill((Jet_Fat[0]+SjetsFarLep[0]+SjetsFarLep[1]).M()*1e-3, weight);
					
					h->h1D("PtGStar_1Fatj2SjFarLepLJetNSeljL", "", s)->Fill((Jet_Fat[0]+SjetsFarLep[0]+SjetsFarLep[1]+e.jet()[sel_jet].mom()+e.electron()[0].mom()+momNu).Perp()*1e-3, weight);
					h->h1D("Pt_1Fatj2SjFarLepLJet", "", s)->Fill((Jet_Fat[0]+SjetsFarLep[0]+SjetsFarLep[1]).Perp()*1e-3, weight);
					h->h1D("Pt_2SjFarLepLJet", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]).Perp()*1e-3, weight);
					h->h1D("Pt_1Fatj", "", s)->Fill((Jet_Fat[0]).Perp()*1e-3, weight);
					h->h1D("Mass_1Fatj", "", s)->Fill((Jet_Fat[0]).M()*1e-3, weight);
					h->h1D("Eta_1Fatj", "", s)->Fill((Jet_Fat[0]).Eta(), weight);
					h->h1D("Pt_SjFarLepLJet0", "", s)->Fill((SjetsFarLep[0]).Perp()*1e-3, weight); 
					h->h1D("Pt_SjFarLepLJet1", "", s)->Fill((SjetsFarLep[1]).Perp()*1e-3, weight); 
					h->h1D("Eta_SjFarLepLJet0", "", s)->Fill((SjetsFarLep[0]).Eta(), weight);
					h->h1D("Eta_SjFarLepLJet1", "", s)->Fill((SjetsFarLep[1]).Eta(), weight);
					
					}
				}    
			}//Electron loop
		
		///////////////////////Muon steam
		else{ 
			h->h1D("Mtop", "", s)->Fill((e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).M()*1e-3, weight);
			h->h1D("Pttop", "", s)->Fill((e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
			
			if (SmaljetsFarLep > 2){
				h->h1D("MassGStar_3SjFarNSeljL", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).M()*1e-3, weight);
				h->h1D("Pt_3SjFar", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[3]).Perp()*1e-3, weight);
				h->h1D("MassT_3SjFar", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[3]).M()*1e-3, weight);
				if (SmaljetsFarLep > 3){
					//h->h1D("MassGStar_4SjFarNSeljL", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).M()*1e-3, weight);
					h->h1D("MassGStar_4jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).M()*1e-3, weight);	
					h->h1D("PtGStar_4SjFarNSeljL", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
					h->h1D("MassT_4jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]).M()*1e-3, weight);
					h->h1D("Pt_4SjFar", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]).Perp()*1e-3, weight);
					}	  
				if (jetsmall > 4){
					h->h1D("MassGStar_5SjFarNSeljL", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+SjetsFarLep[4]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).M()*1e-3, weight);
					h->h1D("MassT_5SjFarNSeljL", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+SjetsFarLep[4]).M()*1e-3, weight);
					h->h1D("PtGStar_5SjFarNSeljL", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+SjetsFarLep[4]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
					h->h1D("Pt_5SjFar", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]+SjetsFarLep[2]+SjetsFarLep[3]+SjetsFarLep[4]).Perp()*1e-3, weight);
					}	
				}
			
			if(fatjets > 0) {
				h->h1D("largeJetM", "", s)->Fill((Jet_Fat[0]).M()*1e-3, weight);
				h->h1D("largeJetPt", "", s)->Fill(Jet_Fat[0].Perp()*1e-3, weight);
				if(SmaljetsFarLep > 1){
					h->h1D("Pt2jets", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]).Perp()*1e-3, weight);
					h->h1D("MassGStar_1Fatj2SjFarLepLJetNSeljL", "", s)->Fill((Jet_Fat[0]+SjetsFarLep[0]+SjetsFarLep[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).M()*1e-3, weight);
					h->h1D("PtGStar_1Fatj2SjFarLepLJetNSeljL", "", s)->Fill((Jet_Fat[0]+SjetsFarLep[0]+SjetsFarLep[1]+e.jet()[sel_jet].mom()+e.muon()[0].mom()+momNu).Perp()*1e-3, weight);
					h->h1D("Pt_1Fatj2SjFarLepLJet", "", s)->Fill((Jet_Fat[0]+SjetsFarLep[0]+SjetsFarLep[1]).Perp()*1e-3, weight);
					h->h1D("Pt_2SjFarLepLJet", "", s)->Fill((SjetsFarLep[0]+SjetsFarLep[1]).Perp()*1e-3, weight);
					h->h1D("MassT_1Fatj2SjFarLepLJetNSeljL", "", s)->Fill((Jet_Fat[0]+SjetsFarLep[0]+SjetsFarLep[1]).M()*1e-3, weight); 
					h->h1D("Pt_1Fatj", "", s)->Fill((Jet_Fat[0]).Perp()*1e-3, weight);
					h->h1D("Mass_1Fatj", "", s)->Fill((Jet_Fat[0]).M()*1e-3, weight);
					h->h1D("Eta_1Fatj", "", s)->Fill((Jet_Fat[0]).Eta(), weight);
					h->h1D("Pt_SjFarLepLJet0", "", s)->Fill((SjetsFarLep[0]).Perp()*1e-3, weight); 
					h->h1D("Pt_SjFarLepLJet1", "", s)->Fill((SjetsFarLep[1]).Perp()*1e-3, weight); 
					h->h1D("Eta_SjFarLepLJet0", "", s)->Fill((SjetsFarLep[0]).Eta(), weight);
					h->h1D("Eta_SjFarLepLJet1", "", s)->Fill((SjetsFarLep[1]).Eta(), weight);
					
					}
				}     
			} //Muon's Loop
		}//e.passReco loop
	
	}//PlotSemilep::run function


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
		} 
	else if (discriminator == 0) {
		float NuZ = 0;
		NuZ = -A*B/(B*B - 1);
		momNu.SetPxPyPzE(NuX, NuY, NuZ, sqrt(std::pow(NuX, 2) + std::pow(NuY, 2) + std::pow(NuZ, 2)));
		momW = momLepton + momNu;
		} 
	else { // two solutions
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
			} 
		else {
			momNu = momNu1;
			momW = momW1;
			}
		}
	
	return true;
	
	} 

