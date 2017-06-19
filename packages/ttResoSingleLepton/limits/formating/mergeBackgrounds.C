#include <TH1D.h>
#include <vector>
#include <iostream>
#include <TFile.h>
#include <map>
#include <TMath.h>
#include <assert.h>
#include <TCanvas.h>
#include <TLegend.h>

//Andrew Altheimer



using namespace std;


void mergeSystematics(string outName, vector<string> inName, TFile * infile, TFile * outfile);


void mergeBackgrounds()
{
  cout<<"mergeCategories \n";
  //string outputFileName="outputs_march23/mergedCategories_J1Pt25.root";
  //TFile * outfile=TFile::Open( (outputFileName).c_str(),"recreate");
  
  //string inputFileName="outputs_march23/topStatSpectra_J1Pt25.root";
  //string inputFileName="outputs_april21/topStatSpectra_EWS_merged.root";
  string inputFileName="outputs_may5/lepPreSpectra_merged_full.root";
  TFile * infile=TFile::Open( (inputFileName).c_str(), "UPDATE");
  
  
  vector<string> systematics;
  vector<string> samples;
  vector<string> channels;
  vector<string> leptons; 
  vector<string> directions;

  
  
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

  systematics.clear();
  systematics.push_back("EleSF");
  systematics.push_back("MuSF");
  systematics.push_back("Btag");
  systematics.push_back("BtagC");
  systematics.push_back("BtagL");

  systematics.push_back("MCGen");
  systematics.push_back("PartonShower");
  systematics.push_back("EWS");
  systematics.push_back("IFSR");
  systematics.push_back("topmass");
  systematics.push_back("JetEnerRes");

  systematics.push_back("JES_ALL");
  systematics.push_back("smallJES");

  systematics.push_back("BoostedJES");
  
  systematics.push_back("norm_tt");
  systematics.push_back("luminosity");
  
  systematics.push_back("BoostedJES0");
  systematics.push_back("BoostedJES1");
  systematics.push_back("BoostedJES2");


  systematics.push_back("Btag0");
  systematics.push_back("Btag1");
  systematics.push_back("Btag2");
  systematics.push_back("Btag3");
  systematics.push_back("Btag4");
  systematics.push_back("Btag5");
  systematics.push_back("Btag6");
  systematics.push_back("Btag7");
  systematics.push_back("Btag8");
  systematics.push_back("Btag9");
  systematics.push_back("Btag10");

  systematics.push_back("BtagC0");
  systematics.push_back("BtagC1");
  systematics.push_back("BtagC2");
  systematics.push_back("BtagC3");
  systematics.push_back("BtagC4");
  systematics.push_back("BtagC5");
  systematics.push_back("BtagC6");
  
  systematics.push_back("BtagL0");
  systematics.push_back("BtagL1");
  systematics.push_back("BtagL2");
  systematics.push_back("BtagL3");
  systematics.push_back("BtagL4");
  systematics.push_back("BtagL5");
  systematics.push_back("BtagL6");
  systematics.push_back("BtagL7");
  systematics.push_back("BtagL8");
  systematics.push_back("BtagL9");
  systematics.push_back("BtagL10");
  systematics.push_back("BtagL11");
  systematics.push_back("BtagL12");
  

  
  systematics.push_back("JES0");  
  systematics.push_back("JES1");
  systematics.push_back("JES2");  
  systematics.push_back("JES3");  
  systematics.push_back("JES4");  
  systematics.push_back("JES5");  
  systematics.push_back("JES6");  
  systematics.push_back("JES7");  
  systematics.push_back("JES8");  
  systematics.push_back("JES9");
  systematics.push_back("JES10");
  systematics.push_back("JES11");
  systematics.push_back("JES12");
  systematics.push_back("JES13");
  //systematics.push_back("JES14");
  systematics.push_back("JES15");
  systematics.push_back("JES16");
  systematics.push_back("JES17");
  systematics.push_back("JES18");
  //systematics.push_back("JES19");
  systematics.push_back("JES20");
  systematics.push_back("JES21");
  systematics.push_back("JES22");
  




  


  //samples.push_back("tt");
  //samples.push_back("W");
  samples.push_back("single-top");
  samples.push_back("Z");
  samples.push_back("Diboson");
  samples.push_back("ttV");
  samples.push_back("QCDe");
  samples.push_back("QCDmu");




  

  string outName;
  vector<string> inName;

  //cout<<samples.size()<<endl;
  //for(unsigned int samp=0; samp<samples.size(); samp++)    {


  for(unsigned int c=0; c<channels.size(); c++)	{
    for(unsigned int l=0; l<leptons.size(); l++)	    {
      cout<<" "<<c<<" "<<l<<endl;
      
      outName="binkkg_smallBgr_"+channels.at(c)+"_"+leptons.at(l);
      
      inName.clear();
      for(unsigned int samp=0; samp<samples.size(); samp++) {
	inName.push_back("binkkg_"+samples.at(samp)+"_"+channels.at(c)+"_"+leptons.at(l));
      }
      mergeSystematics(outName, inName, infile, infile);


      for(unsigned int sys=0; sys<systematics.size(); sys++)		{
	for(unsigned int d=0; d<directions.size(); d++){
	    
	  cout<<" "<<c<<" "<<l<<" "<<sys<<" "<<d<<endl;
	  outName="binkkg_smallBgr_"+channels.at(c)+"_"+leptons.at(l)+"_"+systematics.at(sys)+"_"+directions.at(d);
	  
	  inName.clear();
	  for(unsigned int samp=0; samp<samples.size(); samp++) {
	    inName.push_back("binkkg_"+samples.at(samp)+"_"+channels.at(c)+"_"+leptons.at(l)+"_"+systematics.at(sys)+"_"+directions.at(d));
	  }
	
	  mergeSystematics(outName, inName, infile, infile);
	  //infile->Close();
	    //infile=TFile::Open( (inputFileName).c_str());
	  
	    
	}//direction
      }//sys
    }//leptons
  }//channels
  

  cout<<"closing files \n";
  
  //outfile->Close();
  infile->Close();
  cout<<"finished \n";

}


void mergeSystematics(string outName, vector<string> inName, TFile * infile, TFile * outfile)
{
  cout<<outName<<endl;

  TH1D * spectra=0; 

  for(unsigned int i=0; i<inName.size(); i++){
    cout<<"   "<<inName.at(i)<<endl;
    
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
