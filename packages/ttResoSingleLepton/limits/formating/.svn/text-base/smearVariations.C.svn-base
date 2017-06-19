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

void smearVariations();


void smoothSystematics(string outName, string inName, string nominalName, TFile * infile, TFile * outfile, bool copyNominal, bool debug);
void smearSystematics(string outName, string inName, string nominalName, TFile * infile, TFile * outfile, bool copyNominal, bool debug);


void smearVariations()
{
  cout<<"mergeCategories \n";
  
  
  //update existing file

  string outputFileName="outputs_may9_rush/lepPreSpectra_full_merged_wsmear.root";
  TFile * outfile=TFile::Open( (outputFileName).c_str(), "UPDATE");
  //TFile * outfile=TFile::Open( (outputFileName).c_str(),"recreate");
  
  //string inputFileName="outputs_may9_rush/lepPreSpectra_full_merged.root";
  //TFile * infile=TFile::Open( (inputFileName).c_str());
  TFile * infile=outfile;
  
  vector<string> systematics;
  vector<string> samples;
  vector<string> channels;
  vector<string> leptons; 
  vector<string> directions;

  
  
  //channels.push_back("massTTbarChi2LPC");
  //channels.push_back("masstT");

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


  /*systematics.push_back("Btag0");
  systematics.push_back("Btag1");
  systematics.push_back("Btag2");
  systematics.push_back("Btag3");
  systematics.push_back("Btag4");
  systematics.push_back("Btag5");
  systematics.push_back("Btag6");*/
  systematics.push_back("Btag7");
  systematics.push_back("Btag8");
  systematics.push_back("Btag9");
  systematics.push_back("Btag10");
  systematics.push_back("smallBtag");

  /*systematics.push_back("BtagC0");
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
  */
  
  //systematics.push_back("JES0");
  //systematics.push_back("JES1");
  //systematics.push_back("JES2");
  systematics.push_back("JES3");
  //systematics.push_back("JES4");
  //systematics.push_back("JES5");
  //systematics.push_back("JES6");
  systematics.push_back("JES7");
  //systematics.push_back("JES8");
  //systematics.push_back("JES9");
  //systematics.push_back("JES10");
  //systematics.push_back("JES11");
  systematics.push_back("JES12");
  //systematics.push_back("JES13");
  ////systematics.push_back("JES14");
  //systematics.push_back("JES15");
  //systematics.push_back("JES16");
  //systematics.push_back("JES17");
  systematics.push_back("JES18");
  //systematics.push_back("JES19");
  systematics.push_back("JES20");
  systematics.push_back("JES21");
  systematics.push_back("JES22");

  
  
  
  samples.push_back("tt");
  samples.push_back("W");
  //samples.push_back("single-top");
  //samples.push_back("Z");
  //samples.push_back("Diboson");
  //samples.push_back("ttV");
  //samples.push_back("QCDe");
  //samples.push_back("QCDmu");
  samples.push_back("Bgr");
  samples.push_back("smallBgr");
  
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
  

  string outName;
  vector<string> inName;
  bool copyNominal=false;
  bool debug=false;
  cout<<samples.size()<<endl;
  for(unsigned int samp=0; samp<samples.size(); samp++)    {


    for(unsigned int c=0; c<channels.size(); c++)	{
      for(unsigned int l=0; l<leptons.size(); l++)	    {
	if(infile==outfile) copyNominal=true;//do we copy the nominal spectra?
	else copyNominal=false;

	for(unsigned int sys=0; sys<systematics.size(); sys++)		{
	  for(unsigned int d=0; d<directions.size(); d++){
	    
	    if(sys==2) debug=true;
	    else debug=false;
	    string outname;
	    if(infile==outfile) outname="binkkg_"+samples.at(samp)+"_"+channels.at(c)+"_"+leptons.at(l)+"_"+systematics.at(sys)+"_smeared_"+directions.at(d);
	    else outname="binkkg_"+samples.at(samp)+"_"+channels.at(c)+"_"+leptons.at(l)+"_"+systematics.at(sys)+"_"+directions.at(d);
	    //string outname="binkkg_"+samples.at(samp)+"_"+channels.at(c)+"_"+leptons.at(l)+"_"+systematics.at(sys)+"_"+directions.at(d);
	    string inname="binkkg_"+samples.at(samp)+"_"+channels.at(c)+"_"+leptons.at(l)+"_"+systematics.at(sys)+"_"+directions.at(d);
	    string nominalName="binkkg_"+samples.at(samp)+"_"+channels.at(c)+"_"+leptons.at(l);
	    smearSystematics(outname, inname, nominalName, infile, outfile, copyNominal, debug);
	    copyNominal=false;
	    
	  }
	}
      }
    }
  }

  cout<<"closing files \n";

  outfile->Close();

  cout<<"finished \n";

}


void smoothSystematics(string outName, string inName, string nominalName, TFile * infile, TFile * outfile, bool copyNominal, bool debug)
{
  cout<<"smearSystematics \n";
  cout<<outName<<endl;

  TH1D * outSpectra=(TH1D*) ((TH1D*) infile->Get( inName.c_str()))->Clone(outName.c_str());
  TH1D * inSpectra=((TH1D*) infile->Get( inName.c_str()));
  TH1D * nominalSpectra=((TH1D*) infile->Get( nominalName.c_str()));

  int nbins=inSpectra->GetNbinsX();
  
  int firstBin=0;
  while(nominalSpectra->GetBinContent(firstBin)==0)
    firstBin++;

  cout<<"first bin: "<<firstBin<<endl;

  inSpectra->GetXaxis()->SetRange(firstBin,nbins-1);

  outSpectra->SetTitle(outName.c_str());
  outSpectra->SetName(outName.c_str());

  inSpectra->Divide(nominalSpectra);

  //cout<<"inspectra: ";
  //for(int i=0; i<20; i++)
  //  cout<<inSpectra->GetBinContent(i)<<" ";
  //cout<<endl;
  
  inSpectra->Smooth(1,"R");
  //for(int i=0; i<20; i++)
  //  cout<<inSpectra->GetBinContent(i)<<" ";
  //cout<<endl;

  
  

  for(int b=1; b<nbins+1; b++){
    outSpectra->SetBinContent(b,nominalSpectra->GetBinContent(b)*inSpectra->GetBinContent(b));
    outSpectra->SetBinError(b,nominalSpectra->GetBinError(b)*inSpectra->GetBinContent(b));
  }
 

  outfile->cd();
  outSpectra->Write(outName.c_str());
  if(copyNominal) nominalSpectra->Write();
  delete inSpectra;
  delete outSpectra;
  delete nominalSpectra;
}

void smearSystematics(string outName, string inName, string nominalName, TFile * infile, TFile * outfile, bool copyNominal, bool debug)
{
  cout<<"smearSystematics \n";
  cout<<outName<<endl;

  TH1D * outSpectra=(TH1D*) ((TH1D*) infile->Get( inName.c_str()))->Clone(outName.c_str());
  TH1D * inSpectra=((TH1D*) infile->Get( inName.c_str()));
  TH1D * nominalSpectra=((TH1D*) infile->Get( nominalName.c_str()));

  outSpectra->SetTitle(outName.c_str());
  outSpectra->SetName(outName.c_str());

  double weights[] = {.25,.5,.25};
  inSpectra->Divide(nominalSpectra);

  
  int nbins=inSpectra->GetNbinsX();
  
  double low, mid, high, avg;
  double lowWeight, midWeight, highWeight;
  for(int b=1; b<nbins+1; b++){
    if(nominalSpectra->GetBinContent(b)==0) {
      outSpectra->SetBinContent(b,0);
      outSpectra->SetBinError(b,0);
      continue;
    }
    cout<<b<<" low edge: "<<nominalSpectra->GetBinLowEdge(b)<<" nominal count: "<<nominalSpectra->GetBinContent(b)<<" var ratio: "<<inSpectra->GetBinContent(b)<<endl;
    low=inSpectra->GetBinContent(b-1);
    mid=inSpectra->GetBinContent(b);
    high=inSpectra->GetBinContent(b+1);

    lowWeight=0; midWeight=0; highWeight=0;
    if(nominalSpectra->GetBinContent(b-1)!=0)
      lowWeight=nominalSpectra->GetBinContent(b-1)/nominalSpectra->GetBinError(b-1);
    if(nominalSpectra->GetBinContent(b)!=0)
      midWeight=nominalSpectra->GetBinContent(b)/nominalSpectra->GetBinError(b);
    if(nominalSpectra->GetBinContent(b+1)!=0)
      highWeight=nominalSpectra->GetBinContent(b+1)/nominalSpectra->GetBinError(b+1);
    cout<<"weights:    "<<lowWeight<<" "<<midWeight<<" "<<highWeight<<endl;

    if(low==0) low=mid; 
    if(high==0) high=mid;
   
    avg=(weights[0]*low*lowWeight+weights[1]*mid*midWeight+weights[2]*high*highWeight)/(weights[0]*lowWeight+weights[1]*midWeight+weights[2]*highWeight);
    cout<<"     "<<low<<" "<<mid<<" "<<high<< " "<<avg<<endl;

    outSpectra->SetBinContent(b,nominalSpectra->GetBinContent(b)*avg);
    outSpectra->SetBinError(b,nominalSpectra->GetBinError(b)*avg);
    cout<<b<<" "<<nominalSpectra->GetBinLowEdge(b)<<" "<<outSpectra->GetBinContent(b)<<" "<<outSpectra->GetBinContent(b)/nominalSpectra->GetBinContent(b)<<endl;
    
  }
  

  outfile->cd();
  outSpectra->Write(outName.c_str());
  if(copyNominal) nominalSpectra->Write();
  delete inSpectra;
  delete outSpectra;
  delete nominalSpectra;
}
