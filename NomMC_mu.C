#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cmath>

#include "Rtypes.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TGraphErrors.h"

void NomMC_mu(){
	

	char name[100];
	TFile *file0 = TFile::Open("/afs/cern.ch/work/d/dshoaleh/private/Analysis/ttunf_syst/ttbarhistograms_mu.root","read");
	TFile *file1 = TFile::Open("/afs/cern.ch/work/d/dshoaleh/private/Analysis/ttunf_syst/Dibosonhistograms_mu.root","read");
	TFile *file2 = TFile::Open("/afs/cern.ch/work/d/dshoaleh/private/Analysis/ttunf_syst/qcdDataMuhistograms_mu.root","read");
	TFile *file3 = TFile::Open("/afs/cern.ch/work/d/dshoaleh/private/Analysis/ttunf_syst/SingleTophistograms_mu.root","read");
	TFile *file4 = TFile::Open("/afs/cern.ch/work/d/dshoaleh/private/Analysis/ttunf_syst/ttBosonhistograms_mu.root","read");
	TFile *file5 = TFile::Open("/afs/cern.ch/work/d/dshoaleh/private/Analysis/ttunf_syst/Wjhistograms_mu.root","read");
	TFile *file6 = TFile::Open("/afs/cern.ch/work/d/dshoaleh/private/Analysis/ttunf_syst/Zjhistograms_mu.root","read");
	TFile* outfile = TFile::Open("/afs/cern.ch/work/d/dshoaleh/private/Analysis/ttunf_syst/Nom_mu.root", "recreate");
	
	TH1F *inS2_0= (TH1F*)file0->Get("MassGStar_5SjFarNSeljL");
	inS2_0->SetName("hBkgttbarS2_mu_obs_cuts");
	TH1F *inS2_1= (TH1F*)file1->Get("MassGStar_5SjFarNSeljL");
	inS2_1->SetName("hBkgDibosonS2_mu_obs_cuts");
	TH1F *inS2_2= (TH1F*)file2->Get("MassGStar_5SjFarNSeljL");
	inS2_2->SetName("hBkgQcdS2_mu_obs_cuts");
	TH1F *inS2_3= (TH1F*)file3->Get("MassGStar_5SjFarNSeljL");
	inS2_3->SetName("hBkgSingleTopS2_mu_obs_cuts");
	TH1F *inS2_4= (TH1F*)file4->Get("MassGStar_5SjFarNSeljL");
	inS2_4->SetName("hBkgttBosonS2_mu_obs_cuts");
	TH1F *inS2_5= (TH1F*)file5->Get("MassGStar_5SjFarNSeljL");
	inS2_5->SetName("hBkgWjS2_mu_obs_cuts");
	TH1F *inS2_6= (TH1F*)file6->Get("MassGStar_5SjFarNSeljL");
	inS2_6->SetName("hBkgZjS2_mu_obs_cuts");
//
    inS2_0->SetTitle("MassGStar_5SjFarNSeljL ttbar Nom");
	inS2_1->SetTitle("MassGStar_5SjFarNSeljL DiBoson Nom");
	inS2_2->SetTitle("MassGStar_5SjFarNSeljL qcd Nom");
	inS2_3->SetTitle("MassGStar_5SjFarNSeljL SingleTop Nom");
	inS2_4->SetTitle("MassGStar_5SjFarNSeljL ttBoson Nom");
	inS2_5->SetTitle("MassGStar_5SjFarNSeljL Wj Nom");
	inS2_6->SetTitle("MassGStar_5SjFarNSeljL Zj Nom");
	inS2_0->Integral();
	inS2_1->Integral();
	inS2_2->Integral();
	inS2_3->Integral();
    inS2_4->Integral();
    inS2_5->Integral();
    inS2_6->Integral();
	//std::cout<<<<"inS2_0->Integral()="<<inS2_0->Integral()<<"inS2_1->Integral()="<<inS2_1->Integral()<<std::endl;
        inS2_1->Add(inS2_0);
	//std::cout<<<<"0+1"<<std::endl;
//std::cout<<<<"inS2_1+0->Integral()="<<inS2_1->Integral()<<std::endl;
	inS2_2->Add(inS2_1);
	//std::cout<<<<"1+2"<<std::endl;
	//std::cout<<<<"inS2_1_2->Integral()="<<inS2_2->Integral()<<std::endl;
	inS2_3->Add(inS2_2);
	//std::cout<<<<"2+3"<<std::endl;
	//std::cout<<<<"inS2_23->Integral()="<<inS2_3->Integral()<<std::endl;
	inS2_4->Add(inS2_3);
	//std::cout<<<<"3+4"<<std::endl;
	//std::cout<<<<"inS2_34->Integral()="<<inS2_4->Integral()<<std::endl;
	inS2_5->Add(inS2_4);
	//std::cout<<<<"4+5"<<std::endl;
	//std::cout<<<<"inS2_45->Integral()="<<inS2_5->Integral()<<std::endl;
	inS2_6->Add(inS2_5);
	//std::cout<<<<"5+6"<<std::endl;
	//std::cout<<<<"inS2_56->Integral()="<<inS2_6->Integral()<<std::endl;
  
  //:  inS2_6->Clone();
	TH1F *inS2_=(TH1F*)inS2_6->Clone();
	inS2_->SetName("hBkgMCS2_mu_obs_cuts");
	inS2_->SetTitle("MassGStar_5SjFarNSeljL All Normalized MC Nomminal");
    //std::cout<<<<"6=MCNom"<<std::endl;
    //std::cout<<<<"inS2_MC->Integral()="<<inS2_->Integral()<<std::endl;
//inS2_6->SetDirectory(0);	
	
//	inS2_0->Clone();
//	TH1F *inS2_=(TH1F*)inS2_0->Clone();
//	
//	TH1F *inS2_muSmearIDUP= (TH1F*)file->Get("MassGStar_5SjFarNSeljLmuSmearIDUP");
//	inS2_muSmearIDUP->Clone();
//    TH1F *inS2_muSmearIDHighNorm=(TH1F*)inS2_muSmearIDUP->Clone();             
//	 inS2_muSmearIDHighNorm->SetName("hSigGstartT_MG1000MT600_TtHmuSmearIDHigh_MassGStar_5SjFarNSeljL_mu_obs_cutsNorm");
//    inS2_muSmearIDHighNorm->SetTitle("MassGStar_5SjFarNSeljL GstartT_MG1000MT600_TtH muSmearIDHighNorm");

	TH1F *inS1_0= (TH1F*)file0->Get("MassGStar_3SjFarNSeljL");
	inS1_0->SetName("hBkgttbarS1_mu_obs_cuts");
	TH1F *inS1_1= (TH1F*)file1->Get("MassGStar_3SjFarNSeljL");
	inS1_1->SetName("hBkgDibosonS1_mu_obs_cuts");
	TH1F *inS1_2= (TH1F*)file2->Get("MassGStar_3SjFarNSeljL");
	inS1_2->SetName("hBkgQcdS1_mu_obs_cuts");
	TH1F *inS1_3= (TH1F*)file3->Get("MassGStar_3SjFarNSeljL");
	inS1_3->SetName("hBkgSingleTopS1_mu_obs_cuts");
	TH1F *inS1_4= (TH1F*)file4->Get("MassGStar_3SjFarNSeljL");
	inS1_4->SetName("hBkgttBosonS1_mu_obs_cuts");
	TH1F *inS1_5= (TH1F*)file5->Get("MassGStar_3SjFarNSeljL");
	inS1_5->SetName("hBkgWjS1_mu_obs_cuts");
	TH1F *inS1_6= (TH1F*)file6->Get("MassGStar_3SjFarNSeljL");
	inS1_6->SetName("hBkgZjS1_mu_obs_cuts");
//
    inS1_0->SetTitle("MassGStar_3SjFarNSeljL ttbar Nom");
	inS1_1->SetTitle("MassGStar_3SjFarNSeljL DiBoson Nom");
	inS1_2->SetTitle("MassGStar_3SjFarNSeljL qcd Nom");
	inS1_3->SetTitle("MassGStar_3SjFarNSeljL SingleTop Nom");
	inS1_4->SetTitle("MassGStar_3SjFarNSeljL ttBoson Nom");
	inS1_5->SetTitle("MassGStar_3SjFarNSeljL Wj Nom");
	inS1_6->SetTitle("MassGStar_3SjFarNSeljL Zj Nom");
	inS1_0->Integral();
	inS1_1->Integral();
	inS1_2->Integral();
	inS1_3->Integral();
	inS1_4->Integral();
	inS1_5->Integral();
	inS1_6->Integral();
	//std::cout<<<<"inS1_0->Integral()="<<inS1_0->Integral()<<"inS1_1->Integral()="<<inS1_1->Integral()<<std::endl;
        inS1_1->Add(inS1_0);
	//std::cout<<<<"0+1"<<std::endl;
//std::cout<<<<"inS1_1+0->Integral()="<<inS1_1->Integral()<<std::endl;
	inS1_2->Add(inS1_1);
	//std::cout<<<<"1+2"<<std::endl;
	//std::cout<<<<"inS1_1_2->Integral()="<<inS1_2->Integral()<<std::endl;
	inS1_3->Add(inS1_2);
	//std::cout<<<<"2+3"<<std::endl;
	//std::cout<<<<"inS1_23->Integral()="<<inS1_3->Integral()<<std::endl;
	inS1_4->Add(inS1_3);
	//std::cout<<<<"3+4"<<std::endl;
	//std::cout<<<<"inS1_34->Integral()="<<inS1_4->Integral()<<std::endl;
	inS1_5->Add(inS1_4);
	//std::cout<<<<"4+5"<<std::endl;
	//std::cout<<<<"inS1_45->Integral()="<<inS1_5->Integral()<<std::endl;
	inS1_6->Add(inS1_5);
	//std::cout<<<<"5+6"<<std::endl;
	//std::cout<<<<"inS1_56->Integral()="<<inS1_6->Integral()<<std::endl;
  
  //:  inS1_6->Clone();
	TH1F *inS1_=(TH1F*)inS1_6->Clone();
	inS1_->SetName("hBkgMCS1_mu_obs_cuts");
	inS1_->SetTitle("MassGStar_3SjFarNSeljL All Normalized MC Nomminal");
    //std::cout<<<<"6=MCNom"<<std::endl;
    //std::cout<<<<"inS1_MC->Integral()="<<inS1_->Integral()<<std::endl;
    
    
    
    	
	TH1F *inS3_0= (TH1F*)file0->Get("MassGStar_1Fatj2SjFarLepLJetNSeljL");
	inS3_0->SetName("hBkgttbarS3_mu_obs_cuts");
	TH1F *inS3_1= (TH1F*)file1->Get("MassGStar_1Fatj2SjFarLepLJetNSeljL");
	inS3_1->SetName("hBkgDibosonS3_mu_obs_cuts");
	TH1F *inS3_2= (TH1F*)file2->Get("MassGStar_1Fatj2SjFarLepLJetNSeljL");
	inS3_2->SetName("hBkgQcdS3_mu_obs_cuts");
	TH1F *inS3_3= (TH1F*)file3->Get("MassGStar_1Fatj2SjFarLepLJetNSeljL");
	inS3_3->SetName("hBkgSingleTopS3_mu_obs_cuts");
	TH1F *inS3_4= (TH1F*)file4->Get("MassGStar_1Fatj2SjFarLepLJetNSeljL");
	inS3_4->SetName("hBkgttBosonS3_mu_obs_cuts");
	TH1F *inS3_5= (TH1F*)file5->Get("MassGStar_1Fatj2SjFarLepLJetNSeljL");
	inS3_5->SetName("hBkgWjS3_mu_obs_cuts");
	TH1F *inS3_6= (TH1F*)file6->Get("MassGStar_1Fatj2SjFarLepLJetNSeljL");
	inS3_6->SetName("hBkgZjS3_mu_obs_cuts");
//
    inS3_0->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljL ttbar Nom");
	inS3_1->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljL DiBoson Nom");
	inS3_2->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljL qcd Nom");
	inS3_3->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljL SingleTop Nom");
	inS3_4->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljL ttBoson Nom");
	inS3_5->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljL Wj Nom");
	inS3_6->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljL Zj Nom");
	inS3_0->Integral();
	inS3_1->Integral();
	inS3_2->Integral();
	inS3_3->Integral();
	inS3_4->Integral();
	inS3_5->Integral();
	inS3_6->Integral();
	//std::cout<<<<"inS3_0->Integral()="<<inS3_0->Integral()<<"inS3_1->Integral()="<<inS3_1->Integral()<<std::endl;
        inS3_1->Add(inS3_0);
	//std::cout<<<<"0+1"<<std::endl;
//std::cout<<<<"inS3_1+0->Integral()="<<inS3_1->Integral()<<std::endl;
	inS3_2->Add(inS3_1);
	//std::cout<<<<"1+2"<<std::endl;
	//std::cout<<<<"inS3_1_2->Integral()="<<inS3_2->Integral()<<std::endl;
	inS3_3->Add(inS3_2);
	//std::cout<<<<"2+3"<<std::endl;
	//std::cout<<<<"inS3_23->Integral()="<<inS3_3->Integral()<<std::endl;
	inS3_4->Add(inS3_3);
	//std::cout<<<<"3+4"<<std::endl;
	//std::cout<<<<"inS3_34->Integral()="<<inS3_4->Integral()<<std::endl;
	inS3_5->Add(inS3_4);
	//std::cout<<<<"4+5"<<std::endl;
	//std::cout<<<<"inS3_45->Integral()="<<inS3_5->Integral()<<std::endl;
	inS3_6->Add(inS3_5);
	//std::cout<<<<"5+6"<<std::endl;
	//std::cout<<<<"inS3_56->Integral()="<<inS3_6->Integral()<<std::endl;
  
  //:  inS3_6->Clone();
	TH1F *inS3_=(TH1F*)inS3_6->Clone();
	inS3_->SetName("hBkgMCS3_mu_obs_cuts");
	inS3_->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljL All Normalized MC Nomminal");
    //std::cout<<<<"6=MCNom"<<std::endl;
    //std::cout<<<<"inS3_MC->Integral()="<<inS3_->Integral()<<std::endl;
    
    
    
    	outfile->Write();
	//outfile->ls();
	outfile->Close();

}
