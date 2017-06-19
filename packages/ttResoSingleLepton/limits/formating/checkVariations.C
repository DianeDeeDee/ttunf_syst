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


void checkVariations()
{
  gErrorIgnoreLevel = kWarning;
  
  cout<<"mergeCategories \n";
    
  string inputFileName="PreSpectra_2014Jul29Alldata.root";
  TFile * infile=TFile::Open( (inputFileName).c_str());
  
  
  vector<string> systematics;
  vector<string> samples;
  vector<string> channels;
  vector<string> leptons; 
  vector<string> categories;
  vector<string> directions;

  string prefix="binkkg_";
  
  channels.push_back("massTTbarChi2LPC");
  channels.push_back("masstT");
  
  categories.push_back("cat1");
  categories.push_back("cat2");
  categories.push_back("cat3");


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
  systematics.push_back("PDF");
  systematics.push_back("JetEnerRes");


  systematics.push_back("norm_tt");
  systematics.push_back("luminosity");

  systematics.push_back("BoostedJES0");
  //systematics.push_back("BoostedJES1");
  //systematics.push_back("BoostedJES2");

  systematics.push_back("BoostedJES13");
  systematics.push_back("BoostedJES14");
  systematics.push_back("BoostedJES15");
  systematics.push_back("BoostedJES16");
  systematics.push_back("BoostedJMS");
  systematics.push_back("BoostedJER");
  systematics.push_back("BoostedJMR");

  systematics.push_back("BoostedJES_ALL");

  systematics.push_back("Btag7");
  systematics.push_back("Btag8");
  systematics.push_back("Btag9");
  systematics.push_back("Btag10");
  systematics.push_back("smallBtag");

  


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
  
  systematics.push_back("JES_ALL");
  systematics.push_back("JESsmall");
  systematics.push_back("JVFCut");

  


  

  samples.push_back("Bgr");
  samples.push_back("tt");
  samples.push_back("W");
  samples.push_back("smallBgr");
    
   /*samples.push_back("single-top");
  samples.push_back("Z");
  samples.push_back("Diboson");
  samples.push_back("ttV");
  samples.push_back("QCDe");
  samples.push_back("QCDmu");
  samples.push_back("Bgr");
  samples.push_back("smallBgr");
  */

  /*  samples.push_back("Z400");
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

  */
  

  
  vector<string> inName;

  for(unsigned int samp=0; samp<samples.size(); samp++) {
    for(unsigned int c=0; c<channels.size(); c++) {
      string significantSystematics="";
      string normSystematics="";
      string insignificantSystematics="";
      for(unsigned int sys=0; sys<systematics.size(); sys++) {
	double maxChi2=0;
	double sumChi2=0;
	double maxdN=0;


	
	for(unsigned int cat=0; cat<categories.size(); cat++) {
	  for(unsigned int l=0; l<leptons.size(); l++) {
	    
	    
	    string nominalName=prefix+samples.at(samp)+"_"+channels.at(c)+"_"+categories.at(cat)+"_"+leptons.at(l);
	    string systNameUp= nominalName+"_"+systematics.at(sys)+"_up";
	    string systNameDw= nominalName+"_"+systematics.at(sys)+"_dw";
	    
	    cout<<nominalName<<endl;

	    
	    TKey *key = infile->FindKey( systNameUp.c_str());
	    if (key==0){
	    cout <<systNameUp  << "  histogram does not exist \n";
	      continue;
	     }
	    
	    TH1D* nominal = ((TH1D*) infile->Get( nominalName.c_str()));
	    TH1D* systVarUp = ((TH1D*) infile->Get( systNameUp.c_str()));
	    TH1D* systVarDw = ((TH1D*) infile->Get( systNameDw.c_str()));
	    
	    double chi2Up, chi2Dw, dN, pUp, pDw;
	    
	    
	    chi2Up=nominal->Chi2Test(systVarUp,"WWChi2");
	    chi2Dw=nominal->Chi2Test(systVarDw,"WWChi2");
	    pUp=nominal->Chi2Test(systVarUp,"WW");
	    pDw=nominal->Chi2Test(systVarDw,"WW");
	    dN=(abs(systVarUp->Integral()-nominal->Integral())+abs(systVarDw->Integral()-nominal->Integral()))/(2*nominal->Integral());
	    
	    
	   
	    
	    if(chi2Up>maxChi2) maxChi2=chi2Up;
	    if(chi2Dw>maxChi2) maxChi2=chi2Dw;

	    //cout<<dN<<" "<<maxdN<<endl;
	    if(dN>maxdN) { maxdN=dN;   }
	    sumChi2+=(chi2Up+chi2Dw);
	    
	   
	    cout<<endl;
	    cout<<" "<<systNameUp<<" "<<chi2Up<<"("<<pUp<<") "<<chi2Dw<<"("<<pDw<<") "<<dN<<endl;
	  }//leptons
	}//cat
	cout<<"summary: "<<samples.at(samp)<<" "<<channels.at(c)<<" "<<systematics.at(sys)<<endl;
	cout<<"         "<<maxChi2<<" "<<sumChi2<<" "<<maxdN<<endl;
	
	if(maxChi2>1 ) significantSystematics+=(" "+systematics.at(sys));
	else if(maxdN>.01 ) normSystematics+=(" "+systematics.at(sys));
	else insignificantSystematics+=(" "+systematics.at(sys));
      }//sys
      cout<<"final summary: "<<samples.at(samp)<<" "<<channels.at(c)<<endl;
      cout<<"significant systematics: "<<significantSystematics<<endl;
      cout<<"norm systematics: "<<normSystematics<<endl;
      cout<<"insignificant systematics: "<<insignificantSystematics<<endl;
    } //channel
  }//sample  
      
      

  cout<<"finished";


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
