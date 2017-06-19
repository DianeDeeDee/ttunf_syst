#include <TH1D.h>
#include <vector>
#include <iostream>
#include <TFile.h>
#include <map>
#include <TMath.h>
#include <assert.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "config_list.h"

//Andrew Altheimer



using namespace std;


void mergeSystematics(string outName, vector<string> inName, TFile * infile, TFile * outfile, bool debug=false);


void mergeBackgrounds_full(int mode=1)  // mode: 1=merge Bgr; 2=merge leptons
{
  cout<<"mergeBackgrounds \n";

  TFile * infile=TFile::Open( InputSpectraName.data(), "UPDATE");
  
  
  vector<string> systematics;
  vector<string> smallSamples;
  vector<string> allSamples;
  vector<string> channels;
  vector<string> leptons; 
  vector<string> directions;

  
  channels.push_back("massTTbarChi2LPC");  
  channels.push_back("massTTbarChi2LPC_cat1");
  channels.push_back("massTTbarChi2LPC_cat2");
  channels.push_back("massTTbarChi2LPC_cat3");
  channels.push_back("masstT");
  channels.push_back("masstT_cat1");
  channels.push_back("masstT_cat2");
  channels.push_back("masstT_cat3");

  leptons.push_back("e");
  leptons.push_back("mu");

  directions.push_back("up");
  directions.push_back("dw");

  buildSystList(&systematics);
  


  allSamples.push_back("tt");
  allSamples.push_back("W");
  allSamples.push_back("single-top");
  allSamples.push_back("Z");
  allSamples.push_back("Diboson");
  allSamples.push_back("ttV");
  allSamples.push_back("QCDe");
  allSamples.push_back("QCDmu");

  if (mode==2) {
    allSamples.push_back("Bgr");
    allSamples.push_back("smallBgr");
    allSamples.push_back("data");
    allSamples.push_back("QCD");
  }

  smallSamples.push_back("single-top");
  smallSamples.push_back("Z");
  smallSamples.push_back("Diboson");
  smallSamples.push_back("ttV");
  smallSamples.push_back("QCDe");
  smallSamples.push_back("QCDmu");



  //string prefix="binkkg_";
  string prefix="";

  string outName;
  vector<string> inName;


  for(unsigned int c=0; c<channels.size(); c++)	{

    if (mode==1) {
      cout<<"merge Bgr/smallBgr"<<endl;
      // merge bkgs
      for(unsigned int l=0; l<leptons.size(); l++)	    {
	//cout<<" "<<c<<" "<<l<<endl;
	
	//smallBgr
	outName=prefix+"smallBgr_"+channels.at(c)+"_"+leptons.at(l);
	inName.clear();
	for(unsigned int samp=0; samp<smallSamples.size(); samp++) {
	  inName.push_back(prefix+smallSamples.at(samp)+"_"+channels.at(c)+"_"+leptons.at(l));
	}
	mergeSystematics(outName, inName, infile, infile);
	
	//allBgr
	outName=prefix+"Bgr_"+channels.at(c)+"_"+leptons.at(l);
	inName.clear();
	for(unsigned int samp=0; samp<allSamples.size(); samp++) {
	  inName.push_back(prefix+allSamples.at(samp)+"_"+channels.at(c)+"_"+leptons.at(l));
	}
	mergeSystematics(outName, inName, infile, infile, true);
	
	
	
	for(unsigned int sys=0; sys<systematics.size(); sys++)		{
	  for(unsigned int d=0; d<directions.size(); d++){
	    
	    //cout<<" "<<c<<" "<<l<<" "<<sys<<" "<<d<<endl;
	    
	    //smallBgr
	    outName=prefix+"smallBgr_"+channels.at(c)+"_"+leptons.at(l)+"_"+systematics.at(sys)+"_"+directions.at(d);
	    inName.clear();
	    for(unsigned int samp=0; samp<smallSamples.size(); samp++) {
	      inName.push_back(prefix+smallSamples.at(samp)+"_"+channels.at(c)+"_"+leptons.at(l)+"_"+systematics.at(sys)+"_"+directions.at(d));
	  }
	    
	    mergeSystematics(outName, inName, infile, infile);
	    
	    //allBgr
	    outName=prefix+"Bgr_"+channels.at(c)+"_"+leptons.at(l)+"_"+systematics.at(sys)+"_"+directions.at(d);
	    inName.clear();
	    for(unsigned int samp=0; samp<allSamples.size(); samp++) {
	      inName.push_back(prefix+allSamples.at(samp)+"_"+channels.at(c)+"_"+leptons.at(l)+"_"+systematics.at(sys)+"_"+directions.at(d));
	    }
	    
	    mergeSystematics(outName, inName, infile, infile);
	    
	    
	    //infile->Close();
	    //infile=TFile::Open( (inputFileName).c_str());
	    
	    
	  }//direction
	}//sys
      }//leptons
    }
    if (mode==2) {
      cout<<"merge e/mu"<<endl;
      // merge leptons
      for(unsigned int samp=0; samp<allSamples.size(); samp++) {
	outName=prefix+allSamples.at(samp)+"_"+channels.at(c)+"_"+"lep";
	inName.clear();
	for(unsigned int l=0; l<leptons.size(); l++)	    {
	  if (allSamples.at(samp)=="QCD") {
	    inName.push_back(prefix+allSamples.at(samp)+"e_"+channels.at(c)+"_"+leptons.at(l));
	    inName.push_back(prefix+allSamples.at(samp)+"mu_"+channels.at(c)+"_"+leptons.at(l));
	  }
	  else
	    inName.push_back(prefix+allSamples.at(samp)+"_"+channels.at(c)+"_"+leptons.at(l));
	}
	mergeSystematics(outName, inName, infile, infile, true);
	
	if (allSamples.at(samp)=="data") continue;
	for(unsigned int sys=0; sys<systematics.size(); sys++)		{
	  for(unsigned int d=0; d<directions.size(); d++){
	    outName=prefix+allSamples.at(samp)+"_"+channels.at(c)+"_"+"lep"+"_"+systematics.at(sys)+"_"+directions.at(d);
	    inName.clear();
	    for(unsigned int l=0; l<leptons.size(); l++)	    {
	      if (allSamples.at(samp)=="QCD") {
		inName.push_back(prefix+allSamples.at(samp)+"e_"+channels.at(c)+"_"+leptons.at(l)+"_"+systematics.at(sys)+"_"+directions.at(d));
		inName.push_back(prefix+allSamples.at(samp)+"mu_"+channels.at(c)+"_"+leptons.at(l)+"_"+systematics.at(sys)+"_"+directions.at(d));
	      }
	      else 
		inName.push_back(prefix+allSamples.at(samp)+"_"+channels.at(c)+"_"+leptons.at(l)+"_"+systematics.at(sys)+"_"+directions.at(d));
	    }
	    
	    mergeSystematics(outName, inName, infile, infile);
	    
	  }//direction
	}//sys
      }// samples
    }
  }//channels


  cout<<"closing files \n";
  
  //outfile->Close();
  infile->Close();
  cout<<"finished \n";

}


void mergeSystematics(string outName, vector<string> inName, TFile * infile, TFile * outfile, bool debug)
{
  //cout<<outName<<endl;

  TH1D * spectra=0; 

  for(unsigned int i=0; i<inName.size(); i++){
    //cout<<"   "<<inName.at(i)<<endl;
    TKey *key = infile->FindKey( inName.at(i).c_str());
    if (key==0){
      cout << inName.at(i) << "  histogram does not exist \n";
      continue;
    }


    if(i==0)
      {
	spectra = (TH1D*) ((TH1D*) infile->Get( inName.at(i).c_str()))->Clone(outName.c_str());
   
	spectra->SetTitle(outName.c_str());
	spectra->SetName(outName.c_str());
	if(debug) cout<<"   "<<inName.at(i)<<" "<<spectra->Integral()<<" "<<spectra->Integral()<<endl;


      }
    else
      { 
	TH1D * hist=(TH1D*) ((TH1D*) infile->Get( inName.at(i).c_str()));
	spectra->Add( hist);
	if(debug) cout<<"   "<<inName.at(i)<<" "<<hist->Integral()<<" "<<spectra->Integral()<<endl;


	
	  }
  }

  if(spectra!=0){
    outfile->cd();
    spectra->Write(outName.c_str());

    delete spectra;
  }

}
