#include <TROOT.h>
#include <TH1D.h>
#include <vector>
#include <iostream>
#include <TFile.h>
#include <map>
#include <TMath.h>
#include <TCanvas.h>
#include "config_list.h"

using namespace std;

void smoothSystematics(TH1D* h_var, TH1D* h_nom);
int  GetMaxStatErrBin(TH1D* h0, TH1D* h1, double thr=1e9, bool skipZero=true);
//int  GetNExtrema(TH1D* h1, int mode=0); //mode: 0=extrema; 1=maxima; 2=minima

void smoothSpectra()
{
  cout<<"smoothSystmatics \n";
  
  //update existing file

  //TFile* spectraFile=TFile::Open( ("test_"+InputSpectraName).data(), "UPDATE");
  TFile* spectraFile=TFile::Open( InputSpectraName.data(), "UPDATE");
  //gROOT->cd();
  
  vector<string> systematics;  
  vector<string> samples;
  vector<string> channels;
  vector<string> leptons; 
  vector<string> directions;

  
  
  channels.push_back("massTTbarChi2LPC");
  channels.push_back("masstT");
  channels.push_back("massTTbarChi2LPC_cat1");
  channels.push_back("massTTbarChi2LPC_cat2");
  channels.push_back("massTTbarChi2LPC_cat3");
  channels.push_back("masstT_cat1");
  channels.push_back("masstT_cat2");
  channels.push_back("masstT_cat3");
  

  leptons.push_back("e");
  leptons.push_back("mu");

  directions.push_back("up");
  directions.push_back("dw");


  buildSystList(&systematics, 3);
  
  samples.push_back("tt");
  samples.push_back("W");
//  samples.push_back("single-top");
//  samples.push_back("Z");
//  samples.push_back("Diboson");
//  samples.push_back("ttV");
//  samples.push_back("QCDe");
//  samples.push_back("QCDmu");
  samples.push_back("smallBgr");
  samples.push_back("Bgr");

  buildSigList(&samples);

  //samples.push_back("HH400");
  //samples.push_back("HH500");
  //samples.push_back("HH750");
  //samples.push_back("HH1000");
  //samples.push_back("HH1250");
  //samples.push_back("HH1500");
  //samples.push_back("HH1750");
  //samples.push_back("HH2000");
  //samples.push_back("HH2250");
  //samples.push_back("HH2500");
  //samples.push_back("HH2750");
  //samples.push_back("HH3000");

  for(unsigned int samp=0; samp<samples.size(); samp++)    {
    for(unsigned int c=0; c<channels.size(); c++)	{
      for(unsigned int l=0; l<leptons.size(); l++)	    {
	for(unsigned int sys=0; sys<systematics.size(); sys++)		{
	  for(unsigned int d=0; d<directions.size(); d++){

	    string sysName=samples.at(samp)+"_"+channels.at(c)+"_"+leptons.at(l)+"_"+systematics.at(sys)+"_"+directions.at(d)+";1";
	    string nomName=samples.at(samp)+"_"+channels.at(c)+"_"+leptons.at(l)+";1";
	    TH1D* h_sys=(TH1D*)(spectraFile->Get(sysName.data()));
	    TH1D* h_nom=(TH1D*)(spectraFile->Get(nomName.data()));
	    cout<<sysName<<"\t";  // followed by merged binning
	    smoothSystematics(h_sys, h_nom);
	    h_sys->Write("",TObject::kOverwrite);
	  }
	}
      }
    }
  }

  cout<<"closing files \n";
  spectraFile->Close();

  cout<<"finished \n";

}

int GetMaxStatErrBin(TH1D* h0, TH1D* h1,  double thr, bool skipZero) {
  int iMin=-999;
  double vStatMax=0;
  for (int i=1; i<=h1->GetNbinsX(); ++i) {   // ignore under/overflow bins for now
    double val=h0->GetBinContent(i);
    if (skipZero && val==0) continue;

    double err0=h0->GetBinError(i);
    double err1=h1->GetBinError(i);
    double errR=fabs(sqrt(err0*err0+err1*err1)/val);
    if (errR<thr && errR>vStatMax) {
      vStatMax=errR;
      iMin=i;
    }
  }
  return iMin;
}

//int  GetNExtrema(TH1D* h1, int mode) { //mode: 0=extrema; 1=maxima; 2=minima
//  int nExtrema=0;
//  for (int i=2; i<=h1->GetNbinsX()-1; ++i) {   // ignore under/overflow bins for now
//    double v0=h1->GetBinContent(i);
//    double v1=h1->GetBinContent(i-1);
//    double v2=h1->GetBinContent(i+1);
//    if ((v0-v1)*(v0-v2)>0) {
//      if (mode==1 && v0<v1) continue;
//      if (mode==2 && v0>v1) continue;
//      nExtrema++;
//    }
//  }
//  return nExtrema;
//}


/*
Start from the bin with lowest stat. 
If the syst variation is smaller than stat uncertainty, merge it with the adjacent bin with lower stat.
Else check the next lowest bin.
Restart after merging is done, until no more merging is needed (either single-bin or all bins with syst larger than stat)
*/
void smoothSystematics(TH1D* h_var, TH1D* h_nom)
{
  // merge if syst variation smaller than n-sigma of stat uncertainty (from both var and nom)
  const float systOverStatThreshold=2.; 


  bool verbose=true;
  verbose=false;
  
  TH1D* h_nom_merged=new TH1D(*h_nom);
  TH1D* h_var_merged=new TH1D(*h_var);
  while (h_nom_merged->GetNbinsX()>1) {    // iterate until single-bin, or nothing to merge

    // Obtain current binning
    int nbins =h_nom_merged->GetNbinsX();
    vector<Double_t> xbins;
    xbins.resize(nbins);
    h_nom_merged->GetXaxis()->GetLowEdge(&xbins[0]);
    xbins.push_back(h_nom_merged->GetXaxis()->GetXmax());

    // Search for bins to merge
    int iMerge=-999;
    double maxErrR=1e9;
    while (iMerge==-999) {  // Loop over all bins, starting from the one with largest stat error. Stop if merge happens.
      int iMin=GetMaxStatErrBin(h_nom_merged, h_var_merged, maxErrR);
      if (verbose) cout<<"Check "<<iMin<<"-th bin out of "<<h_var_merged->GetNbinsX()<<" bins"<<endl;
      if (iMin<0) break;
      double var =h_var_merged->GetBinContent(iMin);
      double nom =h_nom_merged->GetBinContent(iMin);
      double stat=sqrt(pow(h_var_merged->GetBinError(iMin),2) + pow(h_nom_merged->GetBinError(iMin),2));
      maxErrR = stat/fabs(nom);
      if (verbose) cout<<Form("yields=%f(%f)\tRelSyst=%f\tRelStat=%f", nom, h_nom_merged->GetBinError(iMin), (var-nom)/nom, maxErrR)<<endl;

      if (fabs(var-nom)/stat<systOverStatThreshold) { 
	// Check if adjacent bin can be merged: low syst/stat ratio, or (syst-syst0)<stat0
	double nomL =h_nom_merged->GetBinContent(iMin-1);
	double systL=h_var_merged->GetBinContent(iMin-1)-nomL;
	double statL=sqrt(pow(h_var_merged->GetBinError(iMin-1),2) + pow(h_nom_merged->GetBinError(iMin-1),2));
	double nomR =h_nom_merged->GetBinContent(iMin+1);
	double systR=h_var_merged->GetBinContent(iMin+1)-nomR;
	double statR=sqrt(pow(h_var_merged->GetBinError(iMin+1),2) + pow(h_nom_merged->GetBinError(iMin+1),2));

	//if (systL*systR<0) continue; // Skip the bins where systematics turning opposite signs (these bins are often merged by following criteria, which remove the shape information)
	
	bool  mergeL=(iMin!=1 && nomL!=0);
	bool  mergeR=(iMin!=nbins && nomR!=0);
	// The neighbor bin is considered not good for merge, if its systematics is statistically significant enough, and can get this bin's systematics over-estimated after merging.
	if ( fabs(systL)/statL>systOverStatThreshold && fabs(systL+var-nom)/(nomL+nom)*nom/stat>systOverStatThreshold ) mergeL=false;
	if ( fabs(systR)/statR>systOverStatThreshold && fabs(systR+var-nom)/(nomR+nom)*nom/stat>systOverStatThreshold ) mergeR=false;
	
	
	
	if (!mergeL && !mergeR) continue;  // Skip bin if both side not appropriate for merge
	iMerge=iMin; // else this bin will be marked as to be merged
	
	if (mergeL && mergeR) {  // Choose the side with worse statistics to merge, if both sides are appropriate
	  double vL=statL/nomL;
	  double vR=statR/nomR;
	  if (vR>vL) mergeL=false;  
	  else mergeR=false;
	}
	
	if (verbose) {
	  cout<<Form("L:%d\tN=%f\tsyst=%f\tstat=%f\tR:%d\tN=%f\tsyst=%f\tstat=%f", mergeL, nomL, systL/nomL, statL/nomL, mergeR, nomR, systR/nomR, statR/nomR)<<endl;
	}
	
	if (mergeR) xbins.erase(xbins.begin()+iMerge);            
        if (mergeL) xbins.erase(xbins.begin()+iMerge-1);
      }
    }
    if (iMerge==-999) break;

    // Rebin and replace 
    TH1D* h_nom_new=(TH1D*)h_nom_merged->Rebin(nbins-1, h_nom_merged->GetName(), &xbins[0]);
    TH1D* h_var_new=(TH1D*)h_var_merged->Rebin(nbins-1, h_var_merged->GetName(), &xbins[0]);
    delete h_nom_merged; 
    delete h_var_merged;
    h_nom_merged=h_nom_new;
    h_var_merged=h_var_new;
  }

  // Take relative variation from merged hists
  h_var_merged->Add(h_var_merged, h_nom_merged, 1, -1);
  h_var_merged->Divide(h_nom_merged);

  //cout<<Form("Merged %d bins into %d bins. %d Extrema remains.", h_var->GetNbinsX(), h_nom_merged->GetNbinsX(), GetNExtrema(h_var_merged))<<endl;
  cout<<Form("Merged %d bins into %d bins.", h_var->GetNbinsX(), h_nom_merged->GetNbinsX())<<endl;


  for (int i=1; i<=h_var->GetNbinsX(); ++i) {  
    double vN=h_nom->GetBinContent(i);
    double vX=h_var->GetXaxis()->GetBinCenter(i);
    int    iX=h_var_merged->GetXaxis()->FindBin(vX); 
    double vY=h_var_merged->GetBinContent(iX);
    if (verbose) cout<<Form("bin %d (%d-th in merged binning): vN=%.3f vOld=%.3f vNew=%.3f relU=%.5f", i, iX, vN, h_var->GetBinContent(i), (1+vY)*vN, vY  )<<endl;
    h_var->SetBinContent(i, (1+vY)*vN);
  }

  delete h_nom_merged;
  delete h_var_merged;

}


