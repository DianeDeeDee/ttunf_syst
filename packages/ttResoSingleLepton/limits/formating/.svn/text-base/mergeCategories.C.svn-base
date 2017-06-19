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


void mergeSystematics(string outName, vector<string> inName, TFile * infile, TFile * outfile);


void mergeCategories()
{
  cout<<"mergeCategories \n";
  string inputFileName="PreSpectra_2014Jul29Alldata.root";
  string outputFileName="PreSpectra_2014Jul29Alldata_CatMerge.root";

  //TFile * infile=TFile::Open( (inputFileName).c_str(), "UPDATE");
  //if(inputFileName==outputFileName) infile=TFile::Open( (inputFileName).c_str(), "UPDATE");
  //else infile=TFile::Open( (inputFileName).c_str());


  //TFile * outfile=infile;
  //if(inputFileName!=outputFileName) outfile=TFile::Open( (outputFileName).c_str(),"recreate");
  //else outfile=infile;
  
  TFile * outfile=TFile::Open( (inputFileName).c_str(), "UPDATE");
  TFile * infile=outfile;
  
  vector<string> systematics;
  vector<string> samples;
  vector<string> channels;
  vector<string> leptons; 
  vector<string> directions;

  
  
  channels.push_back("massTTbarChi2LPC");
  channels.push_back("masstT");

  leptons.push_back("e");
  leptons.push_back("mu");

  directions.push_back("up");
  directions.push_back("dw");

  //directions.push_back("smeared_up");
  //directions.push_back("smeared_dw");

  buildSystList(&systematics);

  
  samples.push_back("tt");
  samples.push_back("W");
  samples.push_back("single-top");
  samples.push_back("Z");
  samples.push_back("Diboson");
  samples.push_back("ttV");
  samples.push_back("QCDe");
  samples.push_back("QCDmu");
  samples.push_back("Bgr");
  samples.push_back("smallBgr");
  //samples.push_back("Z1500");
  samples.push_back("data");


  samples.push_back("Z400");
  samples.push_back("Z500");
  samples.push_back("Z750");
  samples.push_back("Z1000");
  samples.push_back("Z1250");
  samples.push_back("Z1500");
  samples.push_back("Z1750");
  samples.push_back("Z2000");
  samples.push_back("Z2250");
  samples.push_back("Z2500");
  samples.push_back("Z3000");

  samples.push_back("KKg400");
  samples.push_back("KKg500");
  samples.push_back("KKg600");
  samples.push_back("KKg700");
  samples.push_back("KKg800");
  samples.push_back("KKg900");
  samples.push_back("KKg1000");
  samples.push_back("KKg1150");
  samples.push_back("KKg1300");
  samples.push_back("KKg1600");
  samples.push_back("KKg1800");
  samples.push_back("KKg2000");
  samples.push_back("KKg2250");
  samples.push_back("KKg2500");
  samples.push_back("KKg2750");
  samples.push_back("KKg3000");
  
  samples.push_back("KKg1000_width10pc");
  samples.push_back("KKg1000_width15pc");
  samples.push_back("KKg1000_width20pc");
  samples.push_back("KKg1000_width25pc");
  samples.push_back("KKg1000_width30pc");
  samples.push_back("KKg1000_width35pc");
  samples.push_back("KKg1000_width40pc");

  samples.push_back("KKg2000_width10pc");
  samples.push_back("KKg2000_width15pc");
  samples.push_back("KKg2000_width20pc");
  samples.push_back("KKg2000_width25pc");
  samples.push_back("KKg2000_width30pc");
  samples.push_back("KKg2000_width35pc");
  samples.push_back("KKg2000_width40pc");

  samples.push_back("KKg3000_width10pc");
  samples.push_back("KKg3000_width15pc");
  samples.push_back("KKg3000_width20pc");
  samples.push_back("KKg3000_width25pc");
  samples.push_back("KKg3000_width30pc");
  samples.push_back("KKg3000_width35pc");
  samples.push_back("KKg3000_width40pc");

  samples.push_back("RSG400");
  samples.push_back("RSG500");
  samples.push_back("RSG600");
  samples.push_back("RSG700");
  samples.push_back("RSG800");
  samples.push_back("RSG900");
  samples.push_back("RSG1000");
  samples.push_back("RSG1200");
  samples.push_back("RSG1400");
  samples.push_back("RSG1600");
  samples.push_back("RSG1800");
  samples.push_back("RSG2000");
  samples.push_back("RSG2500");

  samples.push_back("HH400");
  samples.push_back("HH500");
  samples.push_back("HH750");
  samples.push_back("HH1000");
  samples.push_back("HH1250");
  samples.push_back("HH1500");
  samples.push_back("HH1750");
  samples.push_back("HH2000");
  samples.push_back("HH2250");
  samples.push_back("HH2500");
  samples.push_back("HH2750");
  samples.push_back("HH3000");

  string outName;
  vector<string> inName;

  cout<<samples.size()<<endl;
  for(unsigned int samp=0; samp<samples.size(); samp++)    {
    //infile->Close();
    //infile=TFile::Open( (inputFileName).c_str());
    //cout<<samp<<" "<<samples.at(samp)<<endl;
    for(unsigned int c=0; c<channels.size(); c++)	{
      for(unsigned int l=0; l<leptons.size(); l++)	    {

	//cout<<samp<<" "<<c<<" "<<l<<endl;
	outName="binkkg_"+samples.at(samp)+"_"+channels.at(c)+"_"+leptons.at(l);
	inName.clear();
	inName.push_back("binkkg_"+samples.at(samp)+"_"+channels.at(c)+"_cat1_"+leptons.at(l));
	inName.push_back("binkkg_"+samples.at(samp)+"_"+channels.at(c)+"_cat2_"+leptons.at(l));
	inName.push_back("binkkg_"+samples.at(samp)+"_"+channels.at(c)+"_cat3_"+leptons.at(l));
	
	mergeSystematics(outName, inName, infile, outfile);
	
	if(samples.at(samp)=="data") continue;

	for(unsigned int sys=0; sys<systematics.size(); sys++)		{
	  for(unsigned int d=0; d<directions.size(); d++){
	    //cout<<samp<<" "<<c<<" "<<l<<" "<<sys<<" "<<d<<endl;
	    outName="binkkg_"+samples.at(samp)+"_"+channels.at(c)+"_"+leptons.at(l)+"_"+systematics.at(sys)+"_"+directions.at(d);
	    inName.clear();
	    inName.push_back("binkkg_"+samples.at(samp)+"_"+channels.at(c)+"_cat1_"+leptons.at(l)+"_"+systematics.at(sys)+"_"+directions.at(d));
	    inName.push_back("binkkg_"+samples.at(samp)+"_"+channels.at(c)+"_cat2_"+leptons.at(l)+"_"+systematics.at(sys)+"_"+directions.at(d));
	    inName.push_back("binkkg_"+samples.at(samp)+"_"+channels.at(c)+"_cat3_"+leptons.at(l)+"_"+systematics.at(sys)+"_"+directions.at(d));

	    mergeSystematics(outName, inName, infile, outfile);
	  }
	}
      }
    }
  }

  cout<<"closing files \n";

  outfile->Close();

  cout<<"finished \n";

}


void mergeSystematics(string outName, vector<string> inName, TFile * infile, TFile * outfile)
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


    if(spectra==0)
      {
	spectra = (TH1D*) ((TH1D*) infile->Get( inName.at(i).c_str()))->Clone(outName.c_str());
   
	spectra->SetTitle(outName.c_str());
	spectra->SetName(outName.c_str());



      }
    else spectra->Add( (TH1D*) ((TH1D*) infile->Get( inName.at(i).c_str())));
    
  }

  if(spectra!=0){
    outfile->cd();
    spectra->Write(outName.c_str());

    delete spectra;
  }

}
