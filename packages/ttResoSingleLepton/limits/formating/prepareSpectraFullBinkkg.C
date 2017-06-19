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


void prepareSpectraFullBinkkg();


//main function for merging resolved spectra
void mergeResolvedSpectra();
//void prepareHadronic();

//copies list of systematics from input to output with formating changes, samples given by 'prefixes.'  nominal spectra contained in nominalFile, while spectraFile contains the variations (can be the same file)
void AddSystematicSet(TFile * nominalFile, TFile *spectraFile, TFile * outfile, vector<string> prefixes, vector<string> systematics);

//combine list of systematics then copy to outfile
void mergeSystematicSet(TFile * nominalFile, TFile *spectraFile, TFile * outfile, vector<string> prefixes, vector<string> systematics, string outSysName);

//copy nominal spectra
void Copy(TFile * input, TFile * outfile ,vector<string> prefixes);

//not currently used
//void SumBgr(TFile * input, TFile * outfile ,vector<string> prefixes, string outname);

//helper function used in many places.  automatically does rebinning
TH1D * GetSpectra(TFile * file, string histo);

//can be used to draw variations (slow, and not very pretty)
void Plot(TH1D * Nom, TH1D * up, TH1D * dw, string name);

//used by various scripts to adjust for the case that the nominal spectra is slightly different in two input files
TH1D* GetScaleFactors(TH1D * Nominal, TH1D * NominalShift);


//used for hadronic channel
//void CopySys(TFile * input, TFile * outfile ,vector<string> prefixes, vector<string> systematics);

//had some problems with root...
void assign(std::vector<double> &vec,int len,double what) ;

//create a normalization systematics (lumi, norm_tt, etc) in the output file
void addNormalizationSys(TFile * NominalSpectra,TFile * outfile,string prefix, string systematic,double down, double up);


//normalization of wjets shape systematics 
//normalization uncertainties
//single sided systematic (?)

//double mttbins_boosted[21];
//double mttbins_resolved[19];
//double mttbins[21];

void prepareSpectraFullBinkkg()
{
  //prepareHadronic();
  mergeResolvedSpectra();
  return;
}

void mergeResolvedSpectra()
{
  //OUTPUTFILE
  string outputFileName="PreSpectra_2014Jul29Alldata.root";
  TFile * outfile=TFile::Open( (outputFileName).c_str(),"recreate");

  //mttbins_boosted[0]=0;
  //for(int b=1;b<=10;b++) mttbins_boosted[b] = mttbins_boosted[b-1]+80.;
  //for(int b=11;b<=15;b++) mttbins_boosted[b] = mttbins_boosted[b-1]+120.;
  //for(int b=16;b<=18;b++) mttbins_boosted[b] = mttbins_boosted[b-1]+200.;
  //for(int b=19;b<=19;b++) mttbins_boosted[b] = mttbins_boosted[b-1]+500.;
  //for(int b=20;b<=20;b++) mttbins_boosted[b] = mttbins_boosted[b-1]+1100.;
  //
  //mttbins_resolved[0]=0;
  //for(int b=1;b<=10;b++) mttbins_resolved[b] = mttbins_resolved[b-1]+80.;
  //for(int b=11;b<=15;b++) mttbins_resolved[b] = mttbins_resolved[b-1]+120.;
  //for(int b=16;b<=17;b++) mttbins_resolved[b] = mttbins_resolved[b-1]+200.;
  //mttbins_resolved[18]=3600;
  //
  ////for(int b=19;b<=19;b++) mttbins_boosted[b] = mttbins_boosted[b-1]+500.;
  ////for(int b=20;b<=20;b++) mttbins_boosted[b] = mttbins[b-1]+1100.;

  //string nominalFileName="inputs/Nominal.root";
  
  //INPUTFILE
  string nominalFileName="inputs/MttSpectra_2014Jul29Alldata.root";   

  TFile * sysSpectra;
  TFile * NominalSpectra=TFile::Open( (nominalFileName).c_str());

  


  vector<string> samples;
  vector<string> SigSamples;
  vector<string> dataSamples;
  vector<string> channels;
  
  dataSamples.push_back("data");
  

  samples.push_back("tt");
  samples.push_back("W");
  samples.push_back("single-top");
  samples.push_back("Z");
  samples.push_back("Diboson");
  samples.push_back("ttV");
  samples.push_back("QCDe");
  samples.push_back("QCDmu");

  buildSigList(&samples);

  buildChList(&channels);
  

  vector<string> prefixes;
  vector<string> SigPrefixes;
  vector<string> dataPrefixes;



  for(unsigned int s=0; s<samples.size(); s++){
    for(unsigned int c=0; c<channels.size(); c++){
      prefixes.push_back(samples.at(s)+"_"+channels.at(c));
    }
  }

  for(unsigned int s=0; s<SigSamples.size(); s++){
    for(unsigned int c=0; c<channels.size(); c++){
      //SigPrefixes.push_back(SigSamples.at(s)+"_"+channels.at(c));
    }
  }


  for(unsigned int s=0; s<dataSamples.size(); s++){
    for(unsigned int c=0; c<channels.size(); c++){
      dataPrefixes.push_back(dataSamples.at(s)+"_"+channels.at(c));
    }
  }



  vector<string> systematics;
  buildSystList(&systematics, false);

  //sysSpectra->Close();
  //MAIN SYS LIST
  sysSpectra=NominalSpectra;

//  systematics.clear();
//  systematics.push_back("EleSF");
//  systematics.push_back("MuSF");
//  //systematics.push_back("Btag");
//  systematics.push_back("BtagC");
//  systematics.push_back("BtagL");
//
//  systematics.push_back("MCGen");
//  systematics.push_back("PartonShower");
//  systematics.push_back("EWS");
//  systematics.push_back("IFSR");
//  systematics.push_back("topmass");
//  systematics.push_back("PDF");
//  systematics.push_back("JetEnerRes");
//
//  systematics.push_back("iqopt3");
//  systematics.push_back("ptjmin10");
//  systematics.push_back("Whfsf");
//
//  systematics.push_back("BoostedJES0");
//  /*systematics.push_back("BoostedJES1");
//  systematics.push_back("BoostedJES2");
//  systematics.push_back("BoostedJES3");
//  systematics.push_back("BoostedJES4");
//  systematics.push_back("BoostedJES5");
//  systematics.push_back("BoostedJES6");
//  systematics.push_back("BoostedJES7");
//  systematics.push_back("BoostedJES8");
//  systematics.push_back("BoostedJES9");
//  systematics.push_back("BoostedJES10");
//  systematics.push_back("BoostedJES11");
//  systematics.push_back("BoostedJES12");*/
//  systematics.push_back("BoostedJES13");
//  systematics.push_back("BoostedJES14");
//  systematics.push_back("BoostedJES15");
//  systematics.push_back("BoostedJES16");
//  systematics.push_back("BoostedJMS");
//  systematics.push_back("BoostedJER");
//  systematics.push_back("BoostedJMR");
//  
//  systematics.push_back("BoostedJES_ALL");
//
//  //systematics.push_back("Btag0");
//  //systematics.push_back("Btag1");
//  //systematics.push_back("Btag2");
//  //systematics.push_back("Btag3");
//  //systematics.push_back("Btag4");
//  //systematics.push_back("Btag5");
//  //systematics.push_back("Btag6");
//  systematics.push_back("Btag7");
//  systematics.push_back("Btag8");
//  systematics.push_back("Btag9");
//  systematics.push_back("Btag10");
//  /*
//  systematics.push_back("BtagC0");
//  systematics.push_back("BtagC1");
//  systematics.push_back("BtagC2");
//  systematics.push_back("BtagC3");
//  systematics.push_back("BtagC4");
//  systematics.push_back("BtagC5");
//  systematics.push_back("BtagC6");
//  
//  systematics.push_back("BtagL0");
//  systematics.push_back("BtagL1");
//  systematics.push_back("BtagL2");
//  systematics.push_back("BtagL3");
//  systematics.push_back("BtagL4");
//  systematics.push_back("BtagL5");
//  systematics.push_back("BtagL6");
//  systematics.push_back("BtagL7");
//  systematics.push_back("BtagL8");
//  systematics.push_back("BtagL9");
//  systematics.push_back("BtagL10");
//  systematics.push_back("BtagL11");
//  systematics.push_back("BtagL12");
//
//  */
//
//  
//  //systematics.push_back("JES0");  
//  //systematics.push_back("JES1");
//  //systematics.push_back("JES2");  
//  systematics.push_back("JES3");  
//  //systematics.push_back("JES4");  
//  //systematics.push_back("JES5");  
//  //systematics.push_back("JES6");  
//  systematics.push_back("JES7");  
//  //systematics.push_back("JES8");  
//  //systematics.push_back("JES9");
//  //systematics.push_back("JES10");
//  //systematics.push_back("JES11");
//  systematics.push_back("JES12");
//  //systematics.push_back("JES13");
//  ////systematics.push_back("JES14"); //zero!
//  //systematics.push_back("JES15");
//  //systematics.push_back("JES16");
//  //systematics.push_back("JES17");
//  //systematics.push_back("JES18");
//  ////systematics.push_back("JES19"); //zero!
//  systematics.push_back("JES20");
//  //systematics.push_back("JES21");
//  //systematics.push_back("JES22");
//  
//  systematics.push_back("JES_ALL");
//  systematics.push_back("JESsmall");
//  systematics.push_back("JVFCut");

  
   
  
  cout<<"AddSystematics \n";
   AddSystematicSet(NominalSpectra, sysSpectra, outfile, prefixes, systematics);
  cout<<"systematics added, prepare small systematics \n";

//  vector<string> smallSystematics;
//  smallSystematics.push_back("JES0");
//  smallSystematics.push_back("JES1");
//  smallSystematics.push_back("JES2");
//  //smallSystematics.push_back("JES3");
//  smallSystematics.push_back("JES4");
//  smallSystematics.push_back("JES5");
//  smallSystematics.push_back("JES6");
//  //smallSystematics.push_back("JES7");
//  smallSystematics.push_back("JES8");
//  smallSystematics.push_back("JES9");
//  smallSystematics.push_back("JES10");
//  smallSystematics.push_back("JES11");
//  //smallSystematics.push_back("JES12");
//  smallSystematics.push_back("JES13");
//  //smallSystematics.push_back("JES14");
//  smallSystematics.push_back("JES15");
//  smallSystematics.push_back("JES16");
//  smallSystematics.push_back("JES17");
//  //smallSystematics.push_back("JES18");
//  //smallSystematics.push_back("JES19");
//  //smallSystematics.push_back("JES20");
//  //smallSystematics.push_back("JES21");
//  //smallSystematics.push_back("JES22");
//
//  cout<<"mergeSystematicSet \n";
//  //mergeSystematicSet(NominalSpectra, sysSpectra, outfile, prefixes, smallSystematics, "smallJES");
//
//  smallSystematics.clear();
//  
//  smallSystematics.push_back("Btag0");
//  smallSystematics.push_back("Btag1");
//  smallSystematics.push_back("Btag2");
//  smallSystematics.push_back("Btag3");
//  smallSystematics.push_back("Btag4");
//  smallSystematics.push_back("Btag5");
//  smallSystematics.push_back("Btag6");
//
//  //mergeSystematicSet(NominalSpectra, sysSpectra, outfile, prefixes, smallSystematics, "smallBtag");
  


  vector<string> normSystematics;
  
  

  
  for(unsigned int c=0; c<channels.size(); c++) {
    for(unsigned int s=0; s<samples.size(); s++) {
      string prefix=samples.at(s)+"_"+channels.at(c);

      if(samples.at(s)=="QCDe" || samples.at(s)=="QCDmu" || samples.at(s)=="W") {
	addNormalizationSys(NominalSpectra,outfile,prefix,"luminosity",0,0 );
      }
      else{
	addNormalizationSys(NominalSpectra,outfile,prefix,"luminosity",-.028,.028 );
      }
      
      if(samples.at(s)=="tt"){
	addNormalizationSys(NominalSpectra,outfile,prefix,"norm_tt", -0.0643, .0606 );
      }
      else{
	addNormalizationSys(NominalSpectra,outfile,prefix,"norm_tt", 0,0);
      }


      if(samples.at(s)=="QCDe"){
	if (prefix.find("massTTbarChi2LPC") != std::string::npos) {
	  addNormalizationSys(NominalSpectra,outfile,prefix,"norm_QCDe", -0.2006, .2006);
	}
	else if (prefix.find("masstT") != std::string::npos) {
	  addNormalizationSys(NominalSpectra,outfile,prefix,"norm_QCDe", -0.1939, .1939);
	}
      }
      else{
        addNormalizationSys(NominalSpectra,outfile,prefix,"norm_QCDe", 0,0);
      }

      if(samples.at(s)=="QCDmu"){
	if (prefix.find("massTTbarChi2LPC") != std::string::npos) {
	  addNormalizationSys(NominalSpectra,outfile,prefix,"norm_QCDmu", -0.2262, .2262);
	}
	else if (prefix.find("masstT") != std::string::npos) {
	  addNormalizationSys(NominalSpectra,outfile,prefix,"norm_QCDmu", -0.1885, .1885);
	}
      }
      else{
        addNormalizationSys(NominalSpectra,outfile,prefix,"norm_QCDmu", 0,0);
      }
      
      
      
    }
    

  }

  
  
  cout<<"Copy Nominal Spectra\n";
  //copy nominal backgrounds
  Copy(NominalSpectra, outfile, prefixes);
  cout<<"Copy Signals Spectra\n";
  Copy(NominalSpectra, outfile, SigPrefixes);

  //copy data?
  Copy(NominalSpectra, outfile, dataPrefixes);
  


  cout<<"closing files \n";
  
  //NominalSpectra->Close();
  outfile->Close();
  cout<<"finished \n";



}


//void prepareHadronic()
//{
//  
//  for(int b=1;b<=10;b++) mttbins[b] = mttbins[b-1]+80.;
//  for(int b=11;b<=15;b++) mttbins[b] = mttbins[b-1]+120.;
//  for(int b=16;b<=18;b++) mttbins[b] = mttbins[b-1]+200.;
//  for(int b=19;b<=19;b++) mttbins[b] = mttbins[b-1]+500.;
//  for(int b=20;b<=20;b++) mttbins[b] = mttbins[b-1]+1100.;
//
//  string inputFileName="inputs/Spectra_all_oldSettings.root";
//  //string inputFileName="inputs/Spectra_all_3.root";
//  TFile * inputSpectra=TFile::Open( (inputFileName).c_str());
// 
//
//  string outputFileName="outputs/hadChannel_oldSetting.root";
//  TFile * outfile=TFile::Open( (outputFileName).c_str(),"recreate");
//
//
//  vector<string> samples;
//  vector<string> channels;
//
//  samples.push_back("tt");
//  //samples.push_back("W");
//  //samples.push_back("single-top");
//  //samples.push_back("Z");
//  //samples.push_back("Diboson");
//  samples.push_back("dataqcd");
//  samples.push_back("Z500");
//  samples.push_back("Z750");
//  samples.push_back("Z1000");
//  samples.push_back("Z1250");
//  samples.push_back("Z1500");
//  samples.push_back("Z1750");
//  samples.push_back("Z2000");
//  samples.push_back("Z2250");
//  samples.push_back("Z2500");
//  samples.push_back("Z3000");
//
//  channels.push_back("allhad");
//  
//
//  vector<string> prefixes;
//
//  
//  
//  for(unsigned int s=0; s<samples.size(); s++){
//    for(unsigned int c=0; c<channels.size(); c++){
//      //prefixes.push_back(samples.at(s)+"_"+channels.at(c));
//      prefixes.push_back(samples.at(s)+"_"+channels.at(c));
//    }
//  }
//
//  vector<string> systematics;
//  //systematics.push_back("Btag");
//  //systematics.push_back("BtagC");
//  //systematics.push_back("BtagL");
//  systematics.push_back("FJES");
//  systematics.push_back("JESE0");
//  systematics.push_back("JESE1");
//  systematics.push_back("QCDB");
//  systematics.push_back("QCDD");
//  systematics.push_back("norm_tt");
//  
//  //systematics.push_back("JVF");
//  //systematics.push_back("BoostedJES0");
//  //systematics.push_back("BoostedJES1");
//  //systematics.push_back("BoostedJES2");
//  
//
//
//
//  
//
//  CopySys(inputSpectra, outfile, prefixes, systematics);
//  SumBgr(inputSpectra, outfile, prefixes, "binkkg_Bgr_allhad_had");
//
//  cout<<"closing files \n";
//
//  inputSpectra->Close();
//  outfile->Close();
//  cout<<"finished \n";
//
//  
//}
//
//void SumBgr(TFile * input, TFile * outfile ,vector<string> prefixes, string outname){
//  cout<<"summing backgrounds "<<outname<<endl;
//  
//  TH1D * allBgr=0;
//
//  for(unsigned int p=0; p<prefixes.size(); p++){
//    cout<<prefixes.at(p)<<endl;
//    if(prefixes.at(p).at(0)=='Z'){ cout<<"skipping \n"; continue;}
//
//    
//    TH1D * Nominal= GetSpectra(input, prefixes.at(p));
//    if(p==0){
//      allBgr= (TH1D *) Nominal->Clone(outname.c_str());
//      allBgr->SetTitle(outname.c_str());
//    }
//    else{
//      allBgr->Add(Nominal);
//    }
//    
//  }
//    
//  outfile->cd();
//  allBgr->Write(outname.c_str());
//}
//
//
//void CopySys(TFile * input, TFile * outfile ,vector<string> prefixes, vector<string> systematics){
//  map<string, float> scale;
//  /*scale["Z500_allhad"]=17.82     * 1.3;
//  scale["Z750_allhad"]=4.31      * 1.3;
//  scale["Z1000_allhad"]=1.24      * 1.3;
//  scale["Z1250_allhad"]=0.441     * 1.3 ;
//  scale["Z1500_allhad"]=0.160     * 1.3;
//  scale["Z1750_allhad"]=0.07      * 1.3;
//  scale["Z2000_allhad"]=0.0275    * 1.3;
//  scale["Z2250_allhad"]=0.013     * 1.3;
//  scale["Z2500_allhad"]=0.0053    * 1.3;
//  scale["Z3000_allhad"]=0.00116   * 1.3;
//  */
//
//  for(unsigned int p=0; p<prefixes.size(); p++){
//
//    double xscale=1;
//    if(scale.count(prefixes.at(p))) xscale=scale[prefixes.at(p)];
//    cout<<prefixes.at(p)<<" "<<xscale<<endl;
//
//    TH1D * Nominal= GetSpectra(input, prefixes.at(p));
//    Nominal->SetTitle(("binkkg_"+prefixes.at(p)+"_had").c_str());
//    Nominal->SetName(("binkkg_"+prefixes.at(p)+"_had").c_str());
//    Nominal->Scale(1/xscale);
//    for (unsigned int s=0; s<systematics.size(); s++){
//      cout<<prefixes.at(p)+"_"+systematics.at(s)<<endl;
//      string systematic =systematics.at(s);
//      if(systematic=="Btag") systematic="hBtag";
//      if(systematic=="BtagC") systematic="hBtagC";
//      if(systematic=="BtagL") systematic="hBtagL";
//      if(systematic=="norm_tt") systematic="hnorm_tt";
//
//      TH1D * up;
//      TH1D * dw;
//
//
//      if(systematics.at(s)=="QCDB" ||systematics.at(s)=="QCDD"){
//	cout<<"QCDB and QCDD not implemented \n";
//	continue;
//	
//      }
//      else{
//	up= GetSpectra(input, prefixes.at(p)+"_"+systematics.at(s)+"_up");
//	dw= GetSpectra(input, prefixes.at(p)+"_"+systematics.at(s)+"_dw");
//	
//	
//
//	if(!up && !dw){ 
//	  cout<<"WARNING: systematics does not exist!!! \n";
//	  up=(TH1D * )Nominal->Clone(("binkkg_"+prefixes.at(p)+"_had_"+systematic+"_up").c_str()); 
//	  dw=(TH1D * )Nominal->Clone(("binkkg_"+prefixes.at(p)+"_had_"+systematic+"_dw").c_str());
//	  up->Scale(1/xscale);
//	  dw->Scale(1/xscale);
//	}
//	else if(!dw || !up) {cout<<"warning: no single sided systematics implemented!";}
//      }//double sided systematics
//	
//      up->Scale(1/xscale);
//      dw->Scale(1/xscale);
//      
//      up->SetTitle(("binkkg_"+prefixes.at(p)+"_had_"+systematic+"_up").c_str());
//      up->SetName(("binkkg_"+prefixes.at(p)+"_had_"+systematic+"_up").c_str());
//      dw->SetTitle(("binkkg_"+prefixes.at(p)+"_had_"+systematic+"_dw").c_str());
//      dw->SetName(("binkkg_"+prefixes.at(p)+"_had_"+systematic+"_dw").c_str());
//      
//      cout<<Nominal->GetName()<<" "<<up->GetName()<<" "<<dw->GetName()<<endl;
//      cout<<Nominal->GetTitle()<<" "<<up->GetTitle()<<" "<<dw->GetTitle()<<endl;
//      cout<<Nominal->Integral()<<" "<<up->Integral()<<" "<<dw->Integral()<<endl;
//      
//      outfile->cd();
//      up->Write(("binkkg_"+prefixes.at(p)+"_had_"+systematic+"_up").c_str());
//      dw->Write(("binkkg_"+prefixes.at(p)+"_had_"+systematic+"_dw").c_str());
//      
//      //check that histogram wrote properly
//      //double upInt=up->Integral();
//      //double dwInt=dw->Integral();
//      //up= GetSpectra(outfile, prefixes.at(p)+"_"+systematic+"_up");
//      //dw= GetSpectra(outfile, prefixes.at(p)+"_"+systematic+"_dw");
//      
//      //assert(upInt==up->Integral());
//      //assert(dwInt==dw->Integral());
//      delete up;
//      delete dw;
//
//    
//    }//end loop over systematics
//  outfile->cd();
//  //double NInt=Nominal->Integral();
//  Nominal->Write(("binkkg_"+prefixes.at(p)+"_had").c_str());
//  //TH1D * Nominal= GetSpectra(outfile,("binkkg_"+prefixes.at(p)+"_allHad").c_str());
//  
//  delete Nominal;
//    }// end loop over prefixes
//}



void Copy(TFile * input, TFile * outfile ,vector<string> prefixes){
  for(unsigned int p=0; p<prefixes.size(); p++){
    cout<<prefixes.at(p)<<endl;
    
    TH1D * Nominal= GetSpectra(input, prefixes.at(p));
         
    outfile->cd();
    
    Nominal->Write(("binkkg_"+prefixes.at(p)).c_str());
    //TH1D * Nominal= GetSpectra(outfile,("binkkg_"+prefixes.at(p)+"_allHad").c_str());

    delete Nominal;
    

  }
}

/*
Notes (Jiahang):
Scale the difference between nominal spectra and systematics spectra (a legacy when systematics spectra were produced by a different group failed to match perfectly with nominal. Not needed anymore)
*/
void AddSystematicSet(TFile * nominalFile, TFile *spectraFile, TFile * outfile, vector<string> prefixes, vector<string> systematics)
{
  bool debug=false;
    
  for(unsigned int p=0; p<prefixes.size(); p++){
    cout<<prefixes.at(p)<<endl;
    
    TH1D * Nominal= GetSpectra(nominalFile, prefixes.at(p));
    TH1D * NominalShift= GetSpectra(spectraFile, prefixes.at(p));
    if(!NominalShift) NominalShift=Nominal;
    TH1D * scaleFactors=GetScaleFactors(Nominal, NominalShift);
    if(debug) cout<<"Nominal: "<<Nominal->Integral()<<" shifted: "<<NominalShift->Integral()<<endl;
     
    for (unsigned int s=0; s<systematics.size(); s++){
       cout<<prefixes.at(p)<<"_"<<systematics.at(s)<<endl;
       TH1D * up=GetSpectra(spectraFile, prefixes.at(p)+"_"+systematics.at(s)+"_up");
       TH1D * dw=GetSpectra(spectraFile, prefixes.at(p)+"_"+systematics.at(s)+"_dw");
       if(!up && ! dw) {
	 cout<<"WARNING: up and dw histograms do not exist, using nominal: "<< prefixes.at(p)<<"_"<<systematics.at(s)<<endl;
	 
	 up=(TH1D*) NominalShift->Clone(("binkkg_"+prefixes.at(p)+"_"+systematics.at(s)+"_up").c_str());
         dw=(TH1D*) NominalShift->Clone(("binkkg_"+prefixes.at(p)+"_"+systematics.at(s)+"_dw").c_str());
	 
	 up->SetTitle(("binkkg_"+prefixes.at(p)+"_"+systematics.at(s)+"_up").c_str());
	 up->SetName(("binkkg_"+prefixes.at(p)+"_"+systematics.at(s)+"_up").c_str());
	 dw->SetTitle(("binkkg_"+prefixes.at(p)+"_"+systematics.at(s)+"_dw").c_str());
	 dw->SetName(("binkkg_"+prefixes.at(p)+"_"+systematics.at(s)+"_dw").c_str());
       }

       //Plot(NominalShift, up, dw, "initial_"+prefixes.at(p)+"_"+systematics.at(s));

       assert(up->GetNbinsX()==dw->GetNbinsX());
       assert(up->GetNbinsX()==Nominal->GetNbinsX());

       if(debug) cout<<"before : "<<up->Integral()<<" "<<NominalShift->Integral()<<" "<<dw->Integral()<<endl;
       if(debug) cout<<"before : "<<up->Integral()/NominalShift->Integral()<<" "<<dw->Integral()/NominalShift->Integral()<<endl;
       
       for(int b=0; b<=up->GetNbinsX()+1; b++)
	 {
	   double up_c, up_err, dw_c, dw_err, sf;
	   sf=scaleFactors->GetBinContent(b);
	   
	   //cout<<b<<" "<<sf<<endl;
	   if(sf<0) {
	     up_c=Nominal->GetBinContent(b);
	     up_err=Nominal->GetBinError(b);
	     dw_c=Nominal->GetBinContent(b);
	     dw_err=Nominal->GetBinError(b);

	   }
	   else{
	     up_c=up->GetBinContent(b)*sf;
	     up_err=up->GetBinError(b)*sf;
	     dw_c=dw->GetBinContent(b)*sf;
	     dw_err=dw->GetBinError(b)*sf;
	   }
	  

	   up->SetBinContent(b, up_c);
	   up->SetBinError(b, up_err);

	   dw->SetBinContent(b, dw_c);
	   dw->SetBinError(b, dw_err);
	   

	 }
       
       if(debug) cout<<"after : "<<up->Integral()<<" "<<Nominal->Integral()<<" "<<dw->Integral()<<endl;
       if(debug) cout<<"after : "<<up->Integral()/Nominal->Integral()<<" "<<dw->Integral()/Nominal->Integral()<<endl;
       double upInt=up->Integral();
       double dwInt=dw->Integral();

       outfile->cd();


       up->Write( ("binkkg_"+prefixes.at(p)+"_"+systematics.at(s)+"_up").c_str());
       dw->Write( ("binkkg_"+prefixes.at(p)+"_"+systematics.at(s)+"_dw").c_str());
       
       
       //check that histogram wrote properly
       up= GetSpectra(outfile, "binkkg_"+prefixes.at(p)+"_"+systematics.at(s)+"_up");
       dw= GetSpectra(outfile, "binkkg_"+prefixes.at(p)+"_"+systematics.at(s)+"_dw");

       if(debug)cout<<"check that histogram wrote properly:(addSystematics) "<<upInt<<" "<<up->Integral()<<" "<<upInt-up->Integral()<<endl;
       
       assert(upInt==up->Integral());
       assert(dwInt==dw->Integral());

       delete up;
       delete dw;

       //Plot(Nominal, up, dw, "final_"+prefixes.at(p)+"_"+systematics.at(s));

    }
  
    if(debug) cout<<"debug: deleting nominal histos \n";
    delete Nominal;
    delete scaleFactors;
      
    //known memory leak here...
    //if(NominalShift) delete NominalShift;
  }
  cout<<"finished AddSystematicSet \n";
}

void AddSystematicSetBtag(TFile * nominalFile, TFile *spectraFile, TFile * outfile, vector<string> prefixes, vector<string> systematics)
{
  cout<<"AddSystematicSetBtag \n";
  
  for(unsigned int p=0; p<prefixes.size(); p++){
    cout<<prefixes.at(p)<<endl;

     TH1D * Nominal= GetSpectra(nominalFile, prefixes.at(p));
     TH1D * NominalShift= GetSpectra(spectraFile, prefixes.at(p));
     if(!NominalShift) NominalShift=Nominal;
     TH1D * scaleFactors=GetScaleFactors(Nominal, NominalShift);

     for (unsigned int s=0; s<systematics.size(); s++){
       string theSys=systematics.at(s);
       int delim=theSys.find("_");
       string theSysBase=theSys.substr(0, delim);
       char chan=theSys.at(delim+1);
       string channel;
       if(chan=='r') channel="massTTbarChi2LPC";
       else if(chan=='b') channel="masstT";
       else cout<<"WARNING: bad channel!!!!!!!!!!!!!!!!!!!!!!\n";
       string cat="cat";
       
       cat+=theSys.at(delim+2);
       cout<<prefixes.at(p)<<"_"<<systematics.at(s)<<" "<<delim<<" "<<theSysBase<<" "<<channel<<" "<<cat<<endl;
       
       bool rightChannel=false;
       //cout<<prefixes.at(p).find(cat)<<" "<<prefixes.at(p).find(channel)<<" "<<string::npos<<endl;
       if(prefixes.at(p).find(cat)!=string::npos && prefixes.at(p).find(channel)!=string::npos)rightChannel=true;

       cout<<"is right channel? "<<rightChannel<<endl;
       
       TH1D * up;
       TH1D * dw;

       if(rightChannel){
	 up=GetSpectra(spectraFile, prefixes.at(p)+"_"+theSysBase+"_up");
	 dw=GetSpectra(spectraFile, prefixes.at(p)+"_"+theSysBase+"_dw");
	 if(!up && !dw){
	   cout<<"Warning, up and dw histograms do not exist, using nominal , AddSystematicSetBtag\n";
	   up=(TH1D*) NominalShift->Clone(("binkkg_"+prefixes.at(p)+"_"+theSysBase+"_up").c_str());
	   dw=(TH1D*) NominalShift->Clone(("binkkg_"+prefixes.at(p)+"_"+theSysBase+"_dw").c_str());
	 }
       }
       else{
	 up=(TH1D*) NominalShift->Clone(("binkkg_"+prefixes.at(p)+"_"+theSysBase+"_up").c_str());
	 dw=(TH1D*) NominalShift->Clone(("binkkg_"+prefixes.at(p)+"_"+theSysBase+"_dw").c_str());
       }
       
       up->SetTitle(("binkkg_"+prefixes.at(p)+"_"+theSys+"_up").c_str());
       up->SetName(("binkkg_"+prefixes.at(p)+"_"+theSys+"_up").c_str());
       dw->SetTitle(("binkkg_"+prefixes.at(p)+"_"+theSys+"_dw").c_str());
       dw->SetName(("binkkg_"+prefixes.at(p)+"_"+theSys+"_dw").c_str());

       
       //Plot(NominalShift, up, dw, "initial_"+prefixes.at(p)+"_"+theSys);

       assert(up->GetNbinsX()==dw->GetNbinsX());
       assert(up->GetNbinsX()==Nominal->GetNbinsX());

       cout<<"before : "<<up->Integral()<<" "<<NominalShift->Integral()<<" "<<dw->Integral()<<endl;
       cout<<"before : "<<up->Integral()/NominalShift->Integral()<<" "<<dw->Integral()/NominalShift->Integral()<<endl;
       
       for(int b=0; b<=up->GetNbinsX()+1; b++)
	 {
	   double up_c, up_err, dw_c, dw_err, sf;
	   sf=scaleFactors->GetBinContent(b);
	   
	   //cout<<b<<" "<<sf<<endl;
	   if(sf<0) {
	     up_c=Nominal->GetBinContent(b);
	     up_err=Nominal->GetBinError(b);
	     dw_c=Nominal->GetBinContent(b);
	     dw_err=Nominal->GetBinError(b);

	   }
	   else{
	     up_c=up->GetBinContent(b)*sf;
	     up_err=up->GetBinError(b)*sf;
	     dw_c=dw->GetBinContent(b)*sf;
	     dw_err=dw->GetBinError(b)*sf;
	   }
	  

	   up->SetBinContent(b, up_c);
	   up->SetBinError(b, up_err);

	   dw->SetBinContent(b, dw_c);
	   dw->SetBinError(b, dw_err);
	   

	 }
       
       cout<<"after : "<<up->Integral()<<" "<<Nominal->Integral()<<" "<<dw->Integral()<<endl;
       cout<<"after : "<<up->Integral()/Nominal->Integral()<<" "<<dw->Integral()/Nominal->Integral()<<endl;
       double upInt=up->Integral();
       double dwInt=dw->Integral();

       outfile->cd();
       up->Write(("binkkg_"+prefixes.at(p)+"_"+theSys+"_up").c_str());
       dw->Write(("binkkg_"+prefixes.at(p)+"_"+theSys+"_dw").c_str());
        
       delete up;
       delete dw;
      
       
       //check that histogram wrote properly
       up= GetSpectra(outfile, prefixes.at(p)+"_"+theSys+"_up");
       dw= GetSpectra(outfile, prefixes.at(p)+"_"+theSys+"_dw");

       cout<<upInt<<" "<<up->Integral()<<endl;
       assert(upInt==up->Integral());
       assert(dwInt==dw->Integral());
       

       //Plot(Nominal, up, dw, "final_"+prefixes.at(p)+"_"+theSys);
       cout<<"test \n";
       delete up;
       cout<<"test1 \n";
       delete dw;
       cout<<"test2 \n";
     }
     cout<<"test3 \n";
     delete Nominal;
     cout<<"test4 "<<Nominal<<" "<<NominalShift<<"\n";
     if(NominalShift!=Nominal)
       delete NominalShift;
     cout<<"test5 \n";

  }


}


void mergeSystematicSet(TFile * nominalFile, TFile *spectraFile, TFile * outfile, vector<string> prefixes, vector<string> systematics, string outSysName)
{
  bool debug=false;
  
  for(unsigned int p=0; p<prefixes.size(); p++){
    cout<<prefixes.at(p)<<endl;
    
    TH1D * Nominal= GetSpectra(nominalFile, prefixes.at(p));
    TH1D * NominalShift= GetSpectra(spectraFile, prefixes.at(p));
    if(!NominalShift) NominalShift=Nominal;
    TH1D * scaleFactors=GetScaleFactors(Nominal, NominalShift);
    if(debug) cout<<"Nominal: "<<Nominal->Integral()<<" shifted: "<<NominalShift->Integral()<<endl;
    
    int nbins=Nominal->GetNbinsX();
    vector<double> upVar2, dwVar2;
    assign(upVar2,nbins+2,0.0);
    assign(dwVar2,nbins+2,0.0);
    
    if(!Nominal || !NominalShift || !scaleFactors) cout<<"error, nominal plots not loaded \n";

    for (unsigned int s=0; s<systematics.size(); s++){
       cout<<prefixes.at(p)<<"_"<<systematics.at(s)<<endl;
       TH1D * up=GetSpectra(spectraFile, prefixes.at(p)+"_"+systematics.at(s)+"_up");
       TH1D * dw=GetSpectra(spectraFile, prefixes.at(p)+"_"+systematics.at(s)+"_dw");
       if(!up || ! dw) {
	 cout<<"Warning, up and dw histograms do not exist! skipping";
	 continue;
       }

       //Plot(NominalShift, up, dw, "initial_"+prefixes.at(p)+"_"+systematics.at(s));

       assert(up->GetNbinsX()==dw->GetNbinsX());
       assert(up->GetNbinsX()==Nominal->GetNbinsX());

       double upVar, dwVar;
       if(debug) cout<<"filling variations<<endl;";
       for(int b=0; b<=up->GetNbinsX()+1; b++)
	 {
	   if(debug) cout<<"  bin "<<b<<endl;
	   upVar=up->GetBinContent(b)-Nominal->GetBinContent(b);
	   dwVar=dw->GetBinContent(b)-Nominal->GetBinContent(b);
	   if(debug) cout<<"    "<<upVar<<" "<<dwVar<<endl;
	   
	   //avoid overflows?
	   //if(abs(upVar)<.0001) upVar=0;
	   //if(abs(dwVar)<.0001) dwVar=0;

	   
	   if(upVar < dwVar) { //variations are reversed
	     upVar+=dwVar;
	     dwVar=upVar-dwVar;
	     upVar-=dwVar;
	   }
 	   if(debug) cout<<"    "<<upVar<<" "<<dwVar<<endl;
	   //signs don't actually matter since we are squaring the variation
	   //this will implicitly flip the sign of the smaller variation if both are in the same direction
	   
	   //actually, lets make sure the signs are correct...
	   if(upVar<0) upVar=0;
	   if(dwVar>0) dwVar=0;

	  
	   upVar2[b]+=upVar*upVar;
	   dwVar2[b]+=dwVar*dwVar;
	   
	   if(debug) cout<<"    "<<upVar2[b]<<" "<<dwVar2[b]<<endl;
	 }
       if(debug) cout<<"deleting up and down "<<up<<" "<<dw<<endl;
       delete up;
       delete dw;
     

    }
    if(debug) cout<<"create summed systematic variation \n";

    
    TH1D * upSum=(TH1D*) Nominal->Clone(("binkkg_"+prefixes.at(p)+"_"+outSysName+"_up").c_str());
    TH1D * dwSum=(TH1D*) Nominal->Clone(("binkkg_"+prefixes.at(p)+"_"+outSysName+"_dw").c_str());
    
    upSum->SetTitle(("binkkg_"+prefixes.at(p)+"_"+outSysName+"_up").c_str());
    upSum->SetName(("binkkg_"+prefixes.at(p)+"_"+outSysName+"_up").c_str());

    dwSum->SetTitle(("binkkg_"+prefixes.at(p)+"_"+outSysName+"_dw").c_str());
    dwSum->SetName(("binkkg_"+prefixes.at(p)+"_"+outSysName+"_dw").c_str());

    
    if(debug) cout<<"start loop \n";
    for(int b=0; b<=Nominal->GetNbinsX()+1; b++)
      {
	if(debug) cout<<"bin "<<b<<" "<<Nominal->GetBinLowEdge(b)<<endl;
	double up_c, up_err, dw_c, dw_err, sf;
	sf=scaleFactors->GetBinContent(b);
	if(debug) cout<<"sf "<<sf<<endl;
	if(debug) cout<<"var2 "<<dwVar2[b]<<" "<<upVar2[b]<<endl;

	//cout<<b<<" "<<sf<<endl;
	if(sf<0) {  //if sf<=0, then the sys uncertainty doesn't really make sense for this bin
	  up_c=Nominal->GetBinContent(b);
	  up_err=Nominal->GetBinError(b);
	  dw_c=Nominal->GetBinContent(b);
	  dw_err=Nominal->GetBinError(b);
	  
	}
	else{ //only the error needs to be scalled to the official nominal spectrum
	  up_c=Nominal->GetBinContent(b)+TMath::Sqrt(upVar2[b])*sf;
	  up_err=Nominal->GetBinError(b);//isn't really used anyway
	  dw_c=Nominal->GetBinContent(b)-TMath::Sqrt(upVar2[b])*sf;
	  dw_err=Nominal->GetBinError(b);
	}
	if(debug) cout<<up_c<<" "<<up_err<<" "<<dw_c<<" "<<dw_err<<endl;


	upSum->SetBinContent(b, up_c);
	upSum->SetBinError(b, up_err);
	
	dwSum->SetBinContent(b, dw_c);
	dwSum->SetBinError(b, dw_err);
	
	
      }
    
    if(debug) cout<<"finished loop over bins! \n";
    double upInt=upSum->Integral();
    double dwInt=dwSum->Integral();

    cout<<"writing histograms \n";
    outfile->cd();
    upSum->Write( ("binkkg_"+prefixes.at(p)+"_"+outSysName+"_up").c_str());
    dwSum->Write( ("binkkg_"+ prefixes.at(p)+"_"+outSysName+"_dw").c_str());
    
    cout<<"getting histograms \n";
    //check that histogram wrote properly
    upSum= GetSpectra(outfile, "binkkg_"+prefixes.at(p)+"_"+outSysName+"_up");
    dwSum= GetSpectra(outfile, "binkkg_"+prefixes.at(p)+"_"+outSysName+"_dw");
    
    if(debug) {
      cout<<"check that histogram wrote properly: " <<upSum<<" "<<dwSum<<endl;;
      cout<<upInt<<" "<<upSum->Integral()<<" "<<upInt-upSum->Integral()<<endl;
      cout<<upInt<<" "<<upSum->Integral()<<" "<<upInt-upSum->Integral()<<endl;
    }
    
    if(abs(upInt-upSum->Integral())>.1 || abs(dwInt-dwSum->Integral())>.1)
      {
	cout<<"WARNING: Histograms did not write properly! \n";
	cout<<upInt<<" "<<upSum->Integral()<<" "<<upInt-upSum->Integral()<<endl;
	cout<<upInt<<" "<<upSum->Integral()<<" "<<upInt-upSum->Integral()<<endl;
      }
    assert( abs(upInt-upSum->Integral())<.1 );
    assert( abs(dwInt-dwSum->Integral())<.1 );
    
    
    
    //Plot(Nominal, upSum, dwSum, "final_"+prefixes.at(p)+"_"+outSysName);
    
    
    //delete up;
    //delete dw;
    
    

    
    //delete Nominal;
    //if(NominalShift) delete NominalShift;
  }
}

void addNormalizationSys(TFile * nominalFile,TFile * outfile,string prefix,string systematic, double downSF, double upSF){

  TH1D * Nominal= GetSpectra(nominalFile, prefix);
  
  TH1D * up=(TH1D*) Nominal->Clone(("binkkg_"+prefix+"_"+systematic+"_up").c_str());
  TH1D * dw=(TH1D*) Nominal->Clone(("binkkg_"+prefix+"_"+systematic+"_dw").c_str());

  up->SetTitle(("binkkg_"+prefix+"_"+systematic+"_up").c_str());
  up->SetName(("binkkg_"+prefix+"_"+systematic+"_up").c_str());
  dw->SetTitle(("binkkg_"+prefix+"_"+systematic+"_dw").c_str());
  dw->SetName(("binkkg_"+prefix+"_"+systematic+"_dw").c_str());

  up->Scale(1+upSF);
  dw->Scale(1+downSF);

  outfile->cd();
  up->Write( ("binkkg_"+prefix+"_"+systematic+"_up").c_str());
  dw->Write( ("binkkg_"+ prefix+"_"+systematic+"_dw").c_str());

  delete up;
  delete dw;
  delete Nominal;

}



TH1D * GetSpectra(TFile * file, string histo){
  bool debug=false;
  bool dropFirstBins=true;
  if(debug) cout<<" getting histo "<< histo<<endl;
  
  TKey *key = file->FindKey( (histo).c_str());
  if (key==0){
    cout << histo << "  histogram does not exist \n";
    return 0;
  }

  TH1D * spectra = (TH1D*) ((TH1D*) file->Get( histo.c_str()));
  //if(histo.find("_masstT_")!=std::string::npos)spectra->Rebin(20,(histo).c_str(), mttbins_boosted);
  //else if(histo.find("massTTbarChi2LPC")!=std::string::npos)spectra->Rebin(20,(histo).c_str(), mttbins_boosted);
  ////else if(histo.find("massTTbarChi2LPC")!=std::string::npos)spectra->Rebin(18,(histo).c_str(), mttbins_resolved);
  //else cout<<"warning: no rebining \n";
  //spectra->SetBinContent(21,0);
  //spectra->SetBinError(21,0);
  ////spectra=(TH1D*) spectra
  
  int binsToDrop=0;
  if(histo.find("_masstT_")!=std::string::npos) binsToDrop=5;
  else if(histo.find("massTTbarChi2LPC")!=std::string::npos) binsToDrop=3;


  if(binsToDrop!=0 && dropFirstBins){
    for(int i=0; i<=binsToDrop; i++){//might as well drop underfil bin as well
      spectra->SetBinContent(i,0);
      spectra->SetBinError(i,0);
    }
  }
      
  spectra->SetTitle(("binkkg_"+histo).c_str());
  spectra->SetName(("binkkg_"+histo).c_str());


  return spectra;
  //return 0;
}


TH1D* GetScaleFactors(TH1D * Nominal, TH1D * NominalShift){
  bool debug=false;
  //cout<<Nominal<<" "<<NominalShift<<endl;
  cout<<Nominal->Integral()<<" "<<NominalShift->Integral()<<endl;
  assert(Nominal->GetNbinsX()==NominalShift->GetNbinsX());
  TH1D * SF=(TH1D *) Nominal->Clone("SFs");
  for(int b=0; b<=Nominal->GetNbinsX()+1; b++)
    {
      if(Nominal->GetBinLowEdge(b) !=  NominalShift->GetBinLowEdge(b)) cout<<"WARNING: binning mismatch! "<<Nominal->GetBinLowEdge(b)<<NominalShift->GetBinLowEdge(b)<<endl;
      
      float num=Nominal->GetBinContent(b);
      float denom=NominalShift->GetBinContent(b);

      float sf=num/denom;
      if(denom==0 && num==0) sf=1;
      else if(denom==0) sf=-1;
      else if(num==0) sf=0;

      if(debug)cout<<Nominal->GetBinLowEdge(b)<<" "<<num<<" "<<denom<<" "<<sf<<endl;
      
      SF->SetBinContent(b,sf);


    }
  
  return SF;
}


void Plot(TH1D * Nom, TH1D * up, TH1D * dw, string name)
{
  cout<<"Plot: "<<name<<endl;
  TCanvas c1;
  c1.SetLogy();
  TLegend *leg_hist=new TLegend(.6,.6, .9, .9);
  leg_hist->SetFillColor(0);
  leg_hist->SetLineColor(0);
  leg_hist->SetShadowColor(0);

  leg_hist->AddEntry(Nom,"Nominal");
  leg_hist->AddEntry(up,"up");
  leg_hist->AddEntry(dw,"dw");
  
  cout<<Nom<<" "<<up<<" "<<dw<<endl;
  if(!(Nom && up && dw)) 
  cout<<"in Plot: nom up dw "<< Nom->Integral()<<" "<<up->Integral()<<" "<<dw->Integral()<<" "<<endl;

  up->SetLineColor(2);
  up->Draw();
  Nom->SetMarkerStyle(1);
  Nom->Draw("same");

  dw->SetLineColor(3);
  dw->Draw("same");

  //leg_hist->Draw("same");
  c1.Print((name+".png").c_str());

}

//apparently root does not support the assign method...
void assign(std::vector<double> &vec,int len,double what) 
{
  vec.erase(vec.begin(),vec.end());
  for(int i=0; i<len; ++i) {
    vec.push_back(what);
  }
}
