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

void NomSys_mu(){
	

	char name[100];
	TFile *file0 = TFile::Open("/afs/cern.ch/work/d/dshoaleh/private/Analysis/ttunf_syst/ttbarhistograms_mu.root","read");
	TFile *file1 = TFile::Open("/afs/cern.ch/work/d/dshoaleh/private/Analysis/ttunf_syst/Dibosonhistograms_mu.root","read");
	TFile *file2 = TFile::Open("/afs/cern.ch/work/d/dshoaleh/private/Analysis/ttunf_syst/qcdDataMuhistograms_mu.root","read");
	TFile *file3 = TFile::Open("/afs/cern.ch/work/d/dshoaleh/private/Analysis/ttunf_syst/SingleTophistograms_mu.root","read");
	TFile *file4 = TFile::Open("/afs/cern.ch/work/d/dshoaleh/private/Analysis/ttunf_syst/ttBosonhistograms_mu.root","read");
	TFile *file5 = TFile::Open("/afs/cern.ch/work/d/dshoaleh/private/Analysis/ttunf_syst/Wjhistograms_mu.root","read");
	TFile *file6 = TFile::Open("/afs/cern.ch/work/d/dshoaleh/private/Analysis/ttunf_syst/Zjhistograms_mu.root","read");
	TFile* outfile = TFile::Open("/afs/cern.ch/work/d/dshoaleh/private/Analysis/ttunf_syst/Sys_mu.root", "recreate");
	
	TH1F *inUpS2_0= (TH1F*)file0->Get("MassGStar_5SjFarNSeljLmuSmearIDUP");
	inUpS2_0->SetName("hBkgttbarUpS2_mu_muSmearID");
	TH1F *inUpS2_1= (TH1F*)file1->Get("MassGStar_5SjFarNSeljLmuSmearIDUP");
	inUpS2_1->SetName("hBkgDibosonUpS2_mu_muSmearID");
//	TH1F *inUpS2_2= (TH1F*)file2->Get("MassGStar_5SjFarNSeljLmuSmearIDUP");
//	inUpS2_2->SetName("hBkgQcdUpS2_mu_muSmearID");
	TH1F *inUpS2_3= (TH1F*)file3->Get("MassGStar_5SjFarNSeljLmuSmearIDUP");
	inUpS2_3->SetName("hBkgSingleTopUpS2_mu_muSmearID");
	TH1F *inUpS2_4= (TH1F*)file4->Get("MassGStar_5SjFarNSeljLmuSmearIDUP");
	inUpS2_4->SetName("hBkgttBosonUpS2_mu_muSmearID");
	TH1F *inUpS2_5= (TH1F*)file5->Get("MassGStar_5SjFarNSeljLmuSmearIDUP");
	inUpS2_5->SetName("hBkgWjUpS2_mu_muSmearID");
	TH1F *inUpS2_6= (TH1F*)file6->Get("MassGStar_5SjFarNSeljLmuSmearIDUP");
	inUpS2_6->SetName("hBkgZjUpS2_mu_muSmearID");
//
    inUpS2_0->SetTitle("MassGStar_5SjFarNSeljLmuSmearIDUP ttbar muSmearID");
	inUpS2_1->SetTitle("MassGStar_5SjFarNSeljLmuSmearIDUP DiBoson muSmearID");
//	inUpS2_2->SetTitle("MassGStar_5SjFarNSeljLmuSmearIDUP qcd muSmearID");
	inUpS2_3->SetTitle("MassGStar_5SjFarNSeljLmuSmearIDUP SingleTop muSmearID");
	inUpS2_4->SetTitle("MassGStar_5SjFarNSeljLmuSmearIDUP ttBoson muSmearID");
	inUpS2_5->SetTitle("MassGStar_5SjFarNSeljLmuSmearIDUP Wj muSmearID");
	inUpS2_6->SetTitle("MassGStar_5SjFarNSeljLmuSmearIDUP Zj muSmearID");
	inUpS2_0->Integral();
	inUpS2_1->Integral();
//	inUpS2_2->Integral();
	inUpS2_3->Integral();
    inUpS2_4->Integral();
    inUpS2_5->Integral();
    inUpS2_6->Integral();
	//std::cout<<<<"inUpS2_0->Integral()="<<inUpS2_0->Integral()<<"inUpS2_1->Integral()="<<inUpS2_1->Integral()<<std::endl;
        inUpS2_1->Add(inUpS2_0);
	//std::cout<<<<"0+1"<<std::endl;
//std::cout<<<<"inUpS2_1+0->Integral()="<<inUpS2_1->Integral()<<std::endl;
//	inUpS2_2->Add(inUpS2_1);
	//std::cout<<<<"1+2"<<std::endl;
	//std::cout<<<<"inUpS2_1_2->Integral()="<<inUpS2_2->Integral()<<std::endl;
	inUpS2_3->Add(inUpS2_1);
	//std::cout<<<<"2+3"<<std::endl;
	//std::cout<<<<"inUpS2_23->Integral()="<<inUpS2_3->Integral()<<std::endl;
	inUpS2_4->Add(inUpS2_3);
	//std::cout<<<<"3+4"<<std::endl;
	//std::cout<<<<"inUpS2_34->Integral()="<<inUpS2_4->Integral()<<std::endl;
	inUpS2_5->Add(inUpS2_4);
	//std::cout<<<<"4+5"<<std::endl;
	//std::cout<<<<"inUpS2_45->Integral()="<<inUpS2_5->Integral()<<std::endl;
	inUpS2_6->Add(inUpS2_5);
	//std::cout<<<<"5+6"<<std::endl;
	//std::cout<<<<"inUpS2_56->Integral()="<<inUpS2_6->Integral()<<std::endl;
  
  //:  inUpS2_6->Clone();
	TH1F *inUpS2_=(TH1F*)inUpS2_6->Clone();
	inUpS2_->SetName("hBkgMCUpS2_mu_muSmearID");
	inUpS2_->SetTitle("MassGStar_5SjFarNSeljLmuSmearIDUP All Normalized MC muSmearIDUP");
    //std::cout<<<<"6=MCNom"<<std::endl;
    //std::cout<<<<"inUpS2_MC->Integral()="<<inUpS2_->Integral()<<std::endl;

	TH1F *inUpS1_0= (TH1F*)file0->Get("MassGStar_3SjFarNSeljLmuSmearIDUP");
	inUpS1_0->SetName("hBkgttbarUpS1_mu_muSmearID");
	TH1F *inUpS1_1= (TH1F*)file1->Get("MassGStar_3SjFarNSeljLmuSmearIDUP");
	inUpS1_1->SetName("hBkgDibosonUpS1_mu_muSmearID");
//	TH1F *inUpS1_2= (TH1F*)file2->Get("MassGStar_3SjFarNSeljLmuSmearIDUP");
//	inUpS1_2->SetName("hBkgQcdUpS1_mu_muSmearID");
	TH1F *inUpS1_3= (TH1F*)file3->Get("MassGStar_3SjFarNSeljLmuSmearIDUP");
	inUpS1_3->SetName("hBkgSingleTopUpS1_mu_muSmearID");
	TH1F *inUpS1_4= (TH1F*)file4->Get("MassGStar_3SjFarNSeljLmuSmearIDUP");
	inUpS1_4->SetName("hBkgttBosonUpS1_mu_muSmearID");
	TH1F *inUpS1_5= (TH1F*)file5->Get("MassGStar_3SjFarNSeljLmuSmearIDUP");
	inUpS1_5->SetName("hBkgWjUpS1_mu_muSmearID");
	TH1F *inUpS1_6= (TH1F*)file6->Get("MassGStar_3SjFarNSeljLmuSmearIDUP");
	inUpS1_6->SetName("hBkgZjUpS1_mu_muSmearID");
//
    inUpS1_0->SetTitle("MassGStar_3SjFarNSeljLmuSmearIDUP ttbar muSmearID");
	inUpS1_1->SetTitle("MassGStar_3SjFarNSeljLmuSmearIDUP DiBoson muSmearID");
//	inUpS1_2->SetTitle("MassGStar_3SjFarNSeljLmuSmearIDUP qcd muSmearID");
	inUpS1_3->SetTitle("MassGStar_3SjFarNSeljLmuSmearIDUP SingleTop muSmearID");
	inUpS1_4->SetTitle("MassGStar_3SjFarNSeljLmuSmearIDUP ttBoson muSmearID");
	inUpS1_5->SetTitle("MassGStar_3SjFarNSeljLmuSmearIDUP Wj muSmearID");
	inUpS1_6->SetTitle("MassGStar_3SjFarNSeljLmuSmearIDUP Zj muSmearID");
	inUpS1_0->Integral();
	inUpS1_1->Integral();
//	inUpS1_2->Integral();
	inUpS1_3->Integral();
	inUpS1_4->Integral();
	inUpS1_5->Integral();
	inUpS1_6->Integral();
	//std::cout<<<<"inUpS1_0->Integral()="<<inUpS1_0->Integral()<<"inUpS1_1->Integral()="<<inUpS1_1->Integral()<<std::endl;
        inUpS1_1->Add(inUpS1_0);
	//std::cout<<<<"0+1"<<std::endl;
//std::cout<<<<"inUpS1_1+0->Integral()="<<inUpS1_1->Integral()<<std::endl;
//	inUpS1_2->Add(inUpS1_1);
	//std::cout<<<<"1+2"<<std::endl;
	//std::cout<<<<"inUpS1_1_2->Integral()="<<inUpS1_2->Integral()<<std::endl;
	inUpS1_3->Add(inUpS1_1);
	//std::cout<<<<"2+3"<<std::endl;
	//std::cout<<<<"inUpS1_23->Integral()="<<inUpS1_3->Integral()<<std::endl;
	inUpS1_4->Add(inUpS1_3);
	//std::cout<<<<"3+4"<<std::endl;
	//std::cout<<<<"inUpS1_34->Integral()="<<inUpS1_4->Integral()<<std::endl;
	inUpS1_5->Add(inUpS1_4);
	//std::cout<<<<"4+5"<<std::endl;
	//std::cout<<<<"inUpS1_45->Integral()="<<inUpS1_5->Integral()<<std::endl;
	inUpS1_6->Add(inUpS1_5);
	//std::cout<<<<"5+6"<<std::endl;
	//std::cout<<<<"inUpS1_56->Integral()="<<inUpS1_6->Integral()<<std::endl;
  
  //:  inUpS1_6->Clone();
	TH1F *inUpS1_=(TH1F*)inUpS1_6->Clone();
	inUpS1_->SetName("hBkgMCUpS1_mu_muSmearID");
	inUpS1_->SetTitle("MassGStar_3SjFarNSeljLmuSmearIDUP All Normalized MC muSmearIDUP");
    //std::cout<<<<"6=MCNom"<<std::endl;
    //std::cout<<<<"inUpS1_MC->Integral()="<<inUpS1_->Integral()<<std::endl;
    
    
    
    	
	TH1F *inUpS3_0= (TH1F*)file0->Get("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDUP");
	inUpS3_0->SetName("hBkgttbarUpS3_mu_muSmearID");
	TH1F *inUpS3_1= (TH1F*)file1->Get("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDUP");
	inUpS3_1->SetName("hBkgDibosonUpS3_mu_muSmearID");
//	TH1F *inUpS3_2= (TH1F*)file2->Get("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDUP");
//	inUpS3_2->SetName("hBkgQcdUpS3_mu_muSmearID");
	TH1F *inUpS3_3= (TH1F*)file3->Get("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDUP");
	inUpS3_3->SetName("hBkgSingleTopUpS3_mu_muSmearID");
	TH1F *inUpS3_4= (TH1F*)file4->Get("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDUP");
	inUpS3_4->SetName("hBkgttBosonUpS3_mu_muSmearID");
	TH1F *inUpS3_5= (TH1F*)file5->Get("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDUP");
	inUpS3_5->SetName("hBkgWjUpS3_mu_muSmearID");
	TH1F *inUpS3_6= (TH1F*)file6->Get("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDUP");
	inUpS3_6->SetName("hBkgZjUpS3_mu_muSmearID");
//
    inUpS3_0->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDUP ttbar muSmearID");
	inUpS3_1->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDUP DiBoson muSmearID");
//	inUpS3_2->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDUP qcd muSmearID");
	inUpS3_3->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDUP SingleTop muSmearID");
	inUpS3_4->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDUP ttBoson muSmearID");
	inUpS3_5->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDUP Wj muSmearID");
	inUpS3_6->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDUP Zj muSmearID");
	inUpS3_0->Integral();
	inUpS3_1->Integral();
//	inUpS3_2->Integral();
	inUpS3_3->Integral();
	inUpS3_4->Integral();
	inUpS3_5->Integral();
	inUpS3_6->Integral();
	//std::cout<<<<"inUpS3_0->Integral()="<<inUpS3_0->Integral()<<"inUpS3_1->Integral()="<<inUpS3_1->Integral()<<std::endl;
        inUpS3_1->Add(inUpS3_0);
	//std::cout<<<<"0+1"<<std::endl;
//std::cout<<<<"inUpS3_1+0->Integral()="<<inUpS3_1->Integral()<<std::endl;
//	inUpS3_2->Add(inUpS3_1);
	//std::cout<<<<"1+2"<<std::endl;
	//std::cout<<<<"inUpS3_1_2->Integral()="<<inUpS3_2->Integral()<<std::endl;
	inUpS3_3->Add(inUpS3_1);
	//std::cout<<<<"2+3"<<std::endl;
	//std::cout<<<<"inUpS3_23->Integral()="<<inUpS3_3->Integral()<<std::endl;
	inUpS3_4->Add(inUpS3_3);
	//std::cout<<<<"3+4"<<std::endl;
	//std::cout<<<<"inUpS3_34->Integral()="<<inUpS3_4->Integral()<<std::endl;
	inUpS3_5->Add(inUpS3_4);
	//std::cout<<<<"4+5"<<std::endl;
	//std::cout<<<<"inUpS3_45->Integral()="<<inUpS3_5->Integral()<<std::endl;
	inUpS3_6->Add(inUpS3_5);
	//std::cout<<<<"5+6"<<std::endl;
	//std::cout<<<<"inUpS3_56->Integral()="<<inUpS3_6->Integral()<<std::endl;
  
  //:  inUpS3_6->Clone();
	TH1F *inUpS3_=(TH1F*)inUpS3_6->Clone();
	inUpS3_->SetName("hBkgMCUpS3_mu_muSmearID");
	inUpS3_->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDUP All Normalized MC Nomminal");
    //std::cout<<<<"6=MCNom"<<std::endl;
    //std::cout<<<<"inUpS3_MC->Integral()="<<inUpS3_->Integral()<<std::endl;
    
    ////////////////////////DOWN
    TH1F *inLowS2_0= (TH1F*)file0->Get("MassGStar_5SjFarNSeljLmuSmearIDLOW");
	inLowS2_0->SetName("hBkgttbarLowS2_mu_muSmearID");
	TH1F *inLowS2_1= (TH1F*)file1->Get("MassGStar_5SjFarNSeljLmuSmearIDLOW");
	inLowS2_1->SetName("hBkgDibosonLowS2_mu_muSmearID");
//	TH1F *inLowS2_2= (TH1F*)file2->Get("MassGStar_5SjFarNSeljLmuSmearIDLOW");
//	inLowS2_2->SetName("hBkgQcdLowS2_mu_muSmearID");
	TH1F *inLowS2_3= (TH1F*)file3->Get("MassGStar_5SjFarNSeljLmuSmearIDLOW");
	inLowS2_3->SetName("hBkgSingleTopLowS2_mu_muSmearID");
	TH1F *inLowS2_4= (TH1F*)file4->Get("MassGStar_5SjFarNSeljLmuSmearIDLOW");
	inLowS2_4->SetName("hBkgttBosonLowS2_mu_muSmearID");
	TH1F *inLowS2_5= (TH1F*)file5->Get("MassGStar_5SjFarNSeljLmuSmearIDLOW");
	inLowS2_5->SetName("hBkgWjLowS2_mu_muSmearID");
	TH1F *inLowS2_6= (TH1F*)file6->Get("MassGStar_5SjFarNSeljLmuSmearIDLOW");
	inLowS2_6->SetName("hBkgZjLowS2_mu_muSmearID");
//
    inLowS2_0->SetTitle("MassGStar_5SjFarNSeljLmuSmearIDLOW ttbar muSmearID");
	inLowS2_1->SetTitle("MassGStar_5SjFarNSeljLmuSmearIDLOW DiBoson muSmearID");
//	inLowS2_2->SetTitle("MassGStar_5SjFarNSeljLmuSmearIDLOW qcd muSmearID");
	inLowS2_3->SetTitle("MassGStar_5SjFarNSeljLmuSmearIDLOW SingleTop muSmearID");
	inLowS2_4->SetTitle("MassGStar_5SjFarNSeljLmuSmearIDLOW ttBoson muSmearID");
	inLowS2_5->SetTitle("MassGStar_5SjFarNSeljLmuSmearIDLOW Wj muSmearID");
	inLowS2_6->SetTitle("MassGStar_5SjFarNSeljLmuSmearIDLOW Zj muSmearID");
	inLowS2_0->Integral();
	inLowS2_1->Integral();
//	inLowS2_2->Integral();
	inLowS2_3->Integral();
    inLowS2_4->Integral();
    inLowS2_5->Integral();
    inLowS2_6->Integral();
	//std::cout<<<<"inLowS2_0->Integral()="<<inLowS2_0->Integral()<<"inLowS2_1->Integral()="<<inLowS2_1->Integral()<<std::endl;
        inLowS2_1->Add(inLowS2_0);
	//std::cout<<<<"0+1"<<std::endl;
//std::cout<<<<"inLowS2_1+0->Integral()="<<inLowS2_1->Integral()<<std::endl;
//	inLowS2_2->Add(inLowS2_1);
	//std::cout<<<<"1+2"<<std::endl;
	//std::cout<<<<"inLowS2_1_2->Integral()="<<inLowS2_2->Integral()<<std::endl;
	inLowS2_3->Add(inLowS2_1);
	//std::cout<<<<"2+3"<<std::endl;
	//std::cout<<<<"inLowS2_23->Integral()="<<inLowS2_3->Integral()<<std::endl;
	inLowS2_4->Add(inLowS2_3);
	//std::cout<<<<"3+4"<<std::endl;
	//std::cout<<<<"inLowS2_34->Integral()="<<inLowS2_4->Integral()<<std::endl;
	inLowS2_5->Add(inLowS2_4);
	//std::cout<<<<"4+5"<<std::endl;
	//std::cout<<<<"inLowS2_45->Integral()="<<inLowS2_5->Integral()<<std::endl;
	inLowS2_6->Add(inLowS2_5);
	//std::cout<<<<"5+6"<<std::endl;
	//std::cout<<<<"inLowS2_56->Integral()="<<inLowS2_6->Integral()<<std::endl;
  
  //:  inLowS2_6->Clone();
	TH1F *inLowS2_=(TH1F*)inLowS2_6->Clone();
	inLowS2_->SetName("hBkgMCLowS2_mu_muSmearID");
	inLowS2_->SetTitle("MassGStar_5SjFarNSeljLmuSmearIDLOW All Normalized MC muSmearIDLOW");
    //std::cout<<<<"6=MCNom"<<std::endl;
    //std::cout<<<<"inLowS2_MC->Integral()="<<inLowS2_->Integral()<<std::endl;

	TH1F *inLowS1_0= (TH1F*)file0->Get("MassGStar_3SjFarNSeljLmuSmearIDLOW");
	inLowS1_0->SetName("hBkgttbarLowS1_mu_muSmearID");
	TH1F *inLowS1_1= (TH1F*)file1->Get("MassGStar_3SjFarNSeljLmuSmearIDLOW");
	inLowS1_1->SetName("hBkgDibosonLowS1_mu_muSmearID");
//	TH1F *inLowS1_2= (TH1F*)file2->Get("MassGStar_3SjFarNSeljLmuSmearIDLOW");
//	inLowS1_2->SetName("hBkgQcdLowS1_mu_muSmearID");
	TH1F *inLowS1_3= (TH1F*)file3->Get("MassGStar_3SjFarNSeljLmuSmearIDLOW");
	inLowS1_3->SetName("hBkgSingleTopLowS1_mu_muSmearID");
	TH1F *inLowS1_4= (TH1F*)file4->Get("MassGStar_3SjFarNSeljLmuSmearIDLOW");
	inLowS1_4->SetName("hBkgttBosonLowS1_mu_muSmearID");
	TH1F *inLowS1_5= (TH1F*)file5->Get("MassGStar_3SjFarNSeljLmuSmearIDLOW");
	inLowS1_5->SetName("hBkgWjLowS1_mu_muSmearID");
	TH1F *inLowS1_6= (TH1F*)file6->Get("MassGStar_3SjFarNSeljLmuSmearIDLOW");
	inLowS1_6->SetName("hBkgZjLowS1_mu_muSmearID");
//
    inLowS1_0->SetTitle("MassGStar_3SjFarNSeljLmuSmearIDLOW ttbar muSmearID");
	inLowS1_1->SetTitle("MassGStar_3SjFarNSeljLmuSmearIDLOW DiBoson muSmearID");
//	inLowS1_2->SetTitle("MassGStar_3SjFarNSeljLmuSmearIDLOW qcd muSmearID");
	inLowS1_3->SetTitle("MassGStar_3SjFarNSeljLmuSmearIDLOW SingleTop muSmearID");
	inLowS1_4->SetTitle("MassGStar_3SjFarNSeljLmuSmearIDLOW ttBoson muSmearID");
	inLowS1_5->SetTitle("MassGStar_3SjFarNSeljLmuSmearIDLOW Wj muSmearID");
	inLowS1_6->SetTitle("MassGStar_3SjFarNSeljLmuSmearIDLOW Zj muSmearID");
	inLowS1_0->Integral();
	inLowS1_1->Integral();
//	inLowS1_2->Integral();
	inLowS1_3->Integral();
	inLowS1_4->Integral();
	inLowS1_5->Integral();
	inLowS1_6->Integral();
	//std::cout<<<<"inLowS1_0->Integral()="<<inLowS1_0->Integral()<<"inLowS1_1->Integral()="<<inLowS1_1->Integral()<<std::endl;
        inLowS1_1->Add(inLowS1_0);
	//std::cout<<<<"0+1"<<std::endl;
//std::cout<<<<"inLowS1_1+0->Integral()="<<inLowS1_1->Integral()<<std::endl;
//	inLowS1_2->Add(inLowS1_1);
	//std::cout<<<<"1+2"<<std::endl;
	//std::cout<<<<"inLowS1_1_2->Integral()="<<inLowS1_2->Integral()<<std::endl;
	inLowS1_3->Add(inLowS1_1);
	//std::cout<<<<"2+3"<<std::endl;
	//std::cout<<<<"inLowS1_23->Integral()="<<inLowS1_3->Integral()<<std::endl;
	inLowS1_4->Add(inLowS1_3);
	//std::cout<<<<"3+4"<<std::endl;
	//std::cout<<<<"inLowS1_34->Integral()="<<inLowS1_4->Integral()<<std::endl;
	inLowS1_5->Add(inLowS1_4);
	//std::cout<<<<"4+5"<<std::endl;
	//std::cout<<<<"inLowS1_45->Integral()="<<inLowS1_5->Integral()<<std::endl;
	inLowS1_6->Add(inLowS1_5);
	//std::cout<<<<"5+6"<<std::endl;
	//std::cout<<<<"inLowS1_56->Integral()="<<inLowS1_6->Integral()<<std::endl;
  
  //:  inLowS1_6->Clone();
	TH1F *inLowS1_=(TH1F*)inLowS1_6->Clone();
	inLowS1_->SetName("hBkgMCLowS1_mu_muSmearID");
	inLowS1_->SetTitle("MassGStar_3SjFarNSeljLmuSmearIDLOW All Normalized MC muSmearIDLOW");
    //std::cout<<<<"6=MCNom"<<std::endl;
    //std::cout<<<<"inLowS1_MC->Integral()="<<inLowS1_->Integral()<<std::endl;
    
    
    
    	
	TH1F *inLowS3_0= (TH1F*)file0->Get("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDLOW");
	inLowS3_0->SetName("hBkgttbarLowS3_mu_muSmearID");
	TH1F *inLowS3_1= (TH1F*)file1->Get("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDLOW");
	inLowS3_1->SetName("hBkgDibosonLowS3_mu_muSmearID");
//	TH1F *inLowS3_2= (TH1F*)file2->Get("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDLOW");
//	inLowS3_2->SetName("hBkgQcdLowS3_mu_muSmearID");
	TH1F *inLowS3_3= (TH1F*)file3->Get("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDLOW");
	inLowS3_3->SetName("hBkgSingleTopLowS3_mu_muSmearID");
	TH1F *inLowS3_4= (TH1F*)file4->Get("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDLOW");
	inLowS3_4->SetName("hBkgttBosonLowS3_mu_muSmearID");
	TH1F *inLowS3_5= (TH1F*)file5->Get("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDLOW");
	inLowS3_5->SetName("hBkgWjLowS3_mu_muSmearID");
	TH1F *inLowS3_6= (TH1F*)file6->Get("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDLOW");
	inLowS3_6->SetName("hBkgZjLowS3_mu_muSmearID");
//
    inLowS3_0->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDLOW ttbar muSmearID");
	inLowS3_1->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDLOW DiBoson muSmearID");
//	inLowS3_2->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDLOW qcd muSmearID");
	inLowS3_3->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDLOW SingleTop muSmearID");
	inLowS3_4->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDLOW ttBoson muSmearID");
	inLowS3_5->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDLOW Wj muSmearID");
	inLowS3_6->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDLOW Zj muSmearID");
	inLowS3_0->Integral();
	inLowS3_1->Integral();
//	inLowS3_2->Integral();
	inLowS3_3->Integral();
	inLowS3_4->Integral();
	inLowS3_5->Integral();
	inLowS3_6->Integral();
	//std::cout<<<<"inLowS3_0->Integral()="<<inLowS3_0->Integral()<<"inLowS3_1->Integral()="<<inLowS3_1->Integral()<<std::endl;
        inLowS3_1->Add(inLowS3_0);
	//std::cout<<<<"0+1"<<std::endl;
//std::cout<<<<"inLowS3_1+0->Integral()="<<inLowS3_1->Integral()<<std::endl;
//	inLowS3_2->Add(inLowS3_1);
	//std::cout<<<<"1+2"<<std::endl;
	//std::cout<<<<"inLowS3_1_2->Integral()="<<inLowS3_2->Integral()<<std::endl;
	inLowS3_3->Add(inLowS3_1);
	//std::cout<<<<"2+3"<<std::endl;
	//std::cout<<<<"inLowS3_23->Integral()="<<inLowS3_3->Integral()<<std::endl;
	inLowS3_4->Add(inLowS3_3);
	//std::cout<<<<"3+4"<<std::endl;
	//std::cout<<<<"inLowS3_34->Integral()="<<inLowS3_4->Integral()<<std::endl;
	inLowS3_5->Add(inLowS3_4);
	//std::cout<<<<"4+5"<<std::endl;
	//std::cout<<<<"inLowS3_45->Integral()="<<inLowS3_5->Integral()<<std::endl;
	inLowS3_6->Add(inLowS3_5);
	//std::cout<<<<"5+6"<<std::endl;
	//std::cout<<<<"inLowS3_56->Integral()="<<inLowS3_6->Integral()<<std::endl;
  
  //:  inLowS3_6->Clone();
	TH1F *inLowS3_=(TH1F*)inLowS3_6->Clone();
	inLowS3_->SetName("hBkgMCLowS3_mu_muSmearID");
	inLowS3_->SetTitle("MassGStar_1Fatj2SjFarLepLJetNSeljLmuSmearIDLOW All Normalized MC Nomminal");
    //std::cout<<<<"6=MCNom"<<std::endl;
    //std::cout<<<<"inLowS3_MC->Integral()="<<inLowS3_->Integral()<<std::endl;
    
    
    
    	outfile->Write();
	//outfile->ls();
	outfile->Close();

}
