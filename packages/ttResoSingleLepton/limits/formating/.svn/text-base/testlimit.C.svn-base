// Based on testlimit.C by Tatjana Lenz and Thorsten Kuhl
// Modified for Glasgow spectra by Elin Kuutmann 2012
// 
// To run:
// make
// ./testlimit &>log

#include <assert.h>
#include <iostream>
#include <TH1D.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>

using namespace std;

TRandom *rnd=0;

int nbins = 720;
double mttbins[720];
//void splitTH1D(const TH1D* orig_nom, const TH1D* orig_shift, TH1D* split, double cutoff, int option);

TH1D* splitTH1D(const TH1D* orig_nom, const TH1D* orig_shift, double cutoff, int option);

TH1D* splitin3TH1D(const TH1D* orig_nom, const TH1D* orig_shift, double cutoff, int option);

// specify your binning here if it is different to the input histogram
void setBins(string binning){
  if(binning=="bin40"){
    nbins=90;
    for(int b=0;b<=91;b++) mttbins[b] = b*40.; 
  } 
  else if(binning=="bin80"){
    nbins=45;
    for(int b=0;b<=46;b++) mttbins[b] = b*80.; 
  } 
  else if(binning=="bin100"){
    nbins=36;
    for(int b=0;b<=37;b++) mttbins[b] = b*100.; 
  } 
  else if(binning=="bincrasy"){
    nbins=36;
    mttbins[0]=0;
    for(int b=1;b<=20;b++) mttbins[b] = mttbins[b-1]+40.;
    for(int b=21;b<=25;b++) mttbins[b] = mttbins[b-1]+60.;
    for(int b=26;b<=30;b++) mttbins[b] = mttbins[b-1]+100.;
    for(int b=31;b<=32;b++) mttbins[b] = mttbins[b-1]+200.;
    for(int b=33;b<=35;b++) mttbins[b] = mttbins[b-1]+500.;
    for(int b=36;b<=36;b++) mttbins[b] = mttbins[b-1]+600.;
  } 
  else if(binning=="binkkg"){
    /*
      nbins=32;
      mttbins[0]=0;
      for(int b=1;b<=15;b++) mttbins[b] = mttbins[b-1]+40.;
      for(int b=16;b<=20;b++) mttbins[b] = mttbins[b-1]+60.;
      for(int b=21;b<=25;b++) mttbins[b] = mttbins[b-1]+100.;
      for(int b=26;b<=28;b++) mttbins[b] = mttbins[b-1]+200.;
      for(int b=29;b<=31;b++) mttbins[b] = mttbins[b-1]+500.;
      for(int b=32;b<=32;b++) mttbins[b] = mttbins[b-1]+600.;
    */
    /*
    nbins=21;
    mttbins[0]=0;
    for(int b=1;b<=10;b++) mttbins[b] = mttbins[b-1]+80.;
    for(int b=11;b<=15;b++) mttbins[b] = mttbins[b-1]+120.;
    for(int b=16;b<=18;b++) mttbins[b] = mttbins[b-1]+200.;
    for(int b=19;b<=20;b++) mttbins[b] = mttbins[b-1]+500.;
    for(int b=21;b<=21;b++) mttbins[b] = mttbins[b-1]+600.;
    */
    nbins=20;
    mttbins[0]=0;
    for(int b=1;b<=10;b++) mttbins[b] = mttbins[b-1]+80.;
    for(int b=11;b<=15;b++) mttbins[b] = mttbins[b-1]+120.;
    for(int b=16;b<=18;b++) mttbins[b] = mttbins[b-1]+200.;
    for(int b=19;b<=19;b++) mttbins[b] = mttbins[b-1]+500.;
    for(int b=20;b<=20;b++) mttbins[b] = mttbins[b-1]+1100.;
  } 
}

//______________________________________________________________________
// routine to add the bins in rebinning

TH1D*  rebinTH1D(TH1D* hN, string binning, bool sys){

  string hname(hN->GetName());
  TH1D* testhist=(TH1D*)hN->Clone();
  setBins(binning);

  //std::cout << " Debug: nbins = " << nbins << endl;
  /*std::cout << " Debug: bin edges: ";
  for(int ijk=0;ijk<=nbins;++ijk) {
    std::cout << ijk << ", " << mttbins[ijk] << "; ";
  }
  std::cout << endl;
  */
  TH1D* rebinned= new TH1D(binning+"_"+hname,binning+"_"+hname,nbins,mttbins);
  rebinned->Reset();
  rebinned->Sumw2();
  for (int ibin=0; ibin< testhist->GetNbinsX()+2;ibin++){
    double oldvalue=testhist->GetBinContent(ibin);
    double olderror=testhist->GetBinError(ibin);
    double center=testhist->GetBinCenter(ibin);
    int bin=rebinned->FindBin(center);
    double value=rebinned->GetBinContent(bin);
    value+=oldvalue;
    //    value+=oldvalue/binwidth;
    rebinned->SetBinContent(bin,value);   
    double error=rebinned->GetBinError(bin);
    if (sys) {
      error = error + olderror;      
      rebinned->SetBinError(bin,error);
    } else {  
      error = error*error+olderror*olderror;        
      rebinned->SetBinError(bin,sqrt(error));
    }       
  }
  return rebinned;
}

TH1F RemoveNegativeEntries(TH1* h_input, double eps=1e-6, bool printDebug=false) {
  int allBins = h_input->GetNbinsX();
  assert (allBins <= 10000);

  if(printDebug) cout << "eps = " << eps << endl;
  if(printDebug) cout << "Histogram " << h_input->GetName() << endl;
  double binc, bine;
  assert (allBins <= 10000);
  double binLimits[10000];
  if(printDebug) cout << "Nbins " << allBins << endl;
  for(int i=1;i<allBins+2;++i)  {
    binc=h_input->GetBinContent(i);
    bine=h_input->GetBinError(i);
    binLimits[i-1] = h_input->GetBinLowEdge(i);
    if(printDebug) cout << "bin " << i << ", low edge " << binLimits[i-1] << ", cont=" << binc << ", err=" << bine << endl;
  }
  string tmpname=h_input->GetName(); tmpname+="_";
  TH1F result(tmpname,h_input->GetTitle(),allBins,binLimits);
  result.Sumw2();

  if(printDebug) cout << "  Replace negative bins" << endl;
  for(int i=0;i<allBins+2;++i) {
    binc=h_input->GetBinContent(i);
    bine=h_input->GetBinError(i);
    if(binc<0) {
      result.SetBinContent(i,eps);
      result.SetBinError(i,bine);
      if(printDebug) {
	if(i<allBins+1) cout << "bin " << i << " (" << binLimits[i] << ", " << binLimits[i+1] << ")";
	else cout << "overflow bin " << i;
	cout << ", cont = " << binc << ", err = " << bine << endl;
      }
    }
    else {
      result.SetBinContent(i,binc);
      result.SetBinError(i,bine);
    }
  }
 
  //Average away 0 bins
  
  //cout << "   0 bin averaging..." << endl;
  for(int i=7;i<allBins;++i) {
    binc=result.GetBinContent(i);
    bine=result.GetBinError(i);
    if(binc<0.9*eps) { 
      
      double previousm=result.GetBinContent(i-1);
      double previousp=result.GetBinContent(i+1);
      double mm=(previousp+previousm+binc)/3.;

      double em=result.GetBinError(i-1);
      double ep=result.GetBinError(i+1);
      double ee=sqrt(em*em+ep*ep+bine*bine)/3.;
      result.SetBinContent(i-1,mm);
      result.SetBinContent(i,mm);
      result.SetBinContent(i+1,mm);
      
      result.SetBinError(i-1,ee);
      result.SetBinError(i,ee);
      result.SetBinError(i+1,ee);
      
      if(printDebug) {
	cout << "bin " << i-1 << " (" << binLimits[i-1] << ", " << binLimits[i] << ")";
	cout << ", cont = " << mm << ", err = " << ee << " (averaged) " << endl;
	cout << "bin " << i << " (" << binLimits[i] << ", " << binLimits[i+1] << ")";
	cout << ", cont = " << mm << ", err = " << ee << " (averaged) " << endl;
	cout << "bin " << i+1 << " (" << binLimits[i+1] << ", " << binLimits[i+2] << ")";
	cout << ", cont = " << mm << ", err = " << ee << " (averaged) " << endl;
      }
    }
  }
  result.SetTitle(h_input->GetTitle());
  result.GetXaxis()->SetTitle(h_input->GetXaxis()->GetTitle());
  result.GetYaxis()->SetTitle(h_input->GetYaxis()->GetTitle());

  if(printDebug) {
    cout << "Final histo cross check" << endl;
    int NN=result.GetNbinsX();
    for(int i=0;i<=NN+1;++i) {
      binc=result.GetBinContent(i);
      bine=result.GetBinError(i);
      cout << "bin " << i << " lowedge=" << result.GetBinLowEdge(i);
      cout << ", cont = " << binc << ", err = " << bine << endl;
    }
  }
  return result;
  
}

TH1D* RemoveFirstBins(TH1* h_input, double threshold, bool printDebug=false) {
  //Set the first bins (up to the bin that has upper edge = threshold)
  // to 0.
  //If the threshold falls between bins, the first bin to contain the threshold is kept.

  int allBins = h_input->GetNbinsX();
  const int Nbinsmax=10000;
  assert (allBins <= Nbinsmax);
  double binLimits[Nbinsmax];
  if(printDebug) cout << "Histogram " << h_input->GetName() << endl;
  //double binlowedge=0., binhighedge=0., 
  double binc=0., bine=0.;
  if(printDebug) cout << "Nbins " << allBins << endl;
  if(printDebug) cout << "Cut-off threshold " << threshold << endl;
  for(int i=1;i<allBins+2;++i)  {
    binc=h_input->GetBinContent(i);
    bine=h_input->GetBinError(i);
    binLimits[i-1] = h_input->GetBinLowEdge(i);
    if(printDebug) cout << "bin " << i << ", low edge " << binLimits[i-1] << ", cont=" << binc << ", err=" << bine << endl;
  }
  string tmpname=h_input->GetName(); //tmpname+="_";
  //TH1D* result(tmpname,h_input->GetTitle(),allBins,binLimits);
  TH1D *result=(TH1D*)h_input->Clone();
  result->SetName(tmpname);
  result->Sumw2();

  if(printDebug) cout << "  Set first bins to 0" << endl;
  for(int i=1;i<=allBins+1;++i)  {
    //binlowedge = h_input->GetBinLowEdge(i);
    //binhighedge = h_input->GetBinLowEdge(i+1);
    binc = h_input->GetBinContent(i);
    bine = h_input->GetBinError(i);
    //if(printDebug) cout << "bin " << i << ", low edge " << binlowedge << ", high edge " << binhighedge << ", cont=" << endl;

    if(binLimits[i]*0.99999<threshold) {
      result->SetBinContent(i,0);
      result->SetBinError(i,0);
      if(printDebug) {
	cout << "Set to 0: bin " << i << " (" << binLimits[i-1] << ", " << binLimits[i] << ")";
	cout << ", nominal cont = " << binc << endl;
      }
    }
    else {
      result->SetBinContent(i,binc);
      result->SetBinError(i,bine);
    }
  }
 
  result->SetTitle(h_input->GetTitle());
  result->GetXaxis()->SetTitle(h_input->GetXaxis()->GetTitle());
  result->GetYaxis()->SetTitle(h_input->GetYaxis()->GetTitle());

  if(printDebug) {
    cout << "Final histo cross check" << endl;
    int NN=result->GetNbinsX();
    for(int i=0;i<=NN+1;++i) {
      binc=result->GetBinContent(i);
      bine=result->GetBinError(i);
      cout << "bin " << i << " lowedge=" << result->GetBinLowEdge(i);
      cout << ", cont = " << binc << ", err = " << bine << endl;
    }
  }
  return result;
  
}


TH1D* splitTH1D(const TH1D* orig_nom, const TH1D* orig_shift, double cutoff, int option){

  //option: 
  //0-- lower end is shifted, higher end not
  //1-- higher end is shifted, lower end is not

  string hname(orig_shift->GetName());
  TH1D *split=(TH1D*)orig_nom->Clone();

  split->Sumw2();
  bool lower=false, higher=false;

  for (int ibin=0; ibin< orig_shift->GetNbinsX()+2; ibin++){
    double shiftvalue=orig_shift->GetBinContent(ibin);
    double shifterror=orig_shift->GetBinError(ibin);
    double nomvalue=orig_nom->GetBinContent(ibin);
    double nomerror=orig_nom->GetBinError(ibin);

    double center=orig_shift->GetBinCenter(ibin);
    //cout << "bin " << ibin
    if(center<=cutoff) {
      if(option==0) {
	split->SetBinContent(ibin, shiftvalue);
	split->SetBinError(ibin, shifterror);
      }
      else {
	split->SetBinContent(ibin, nomvalue);
	split->SetBinError(ibin, nomerror);
      }
      if(!lower) lower=true;
    }
    else { //center>cutoff
      if(option==1) {
	split->SetBinContent(ibin, shiftvalue);
	split->SetBinError(ibin, shifterror);
      }
      else {
	split->SetBinContent(ibin, nomvalue);
	split->SetBinError(ibin, nomerror);
      }
      if(!higher) higher=true;
    }

  }

  if(!lower) cout << "WARNING! Histogram split: " << cutoff << " lower than all bin values" << endl;
  if(!higher) cout << "WARNING! Histogram split: " << cutoff << " higher than all bin values" << endl;
  if(!lower || !higher) cout << "Low bin " <<  orig_shift->GetBinCenter(1) << ", high bin " << orig_shift->GetBinCenter(orig_shift->GetNbinsX()+1) << ", nbins " << orig_shift->GetNbinsX() << endl;

  return split;

}

TH1D* splitin3TH1D(const TH1D* orig_nom, const TH1D* orig_shift, double co1, double co2, int option){

  //option: 
  //0-- lower end is shifted, the others not
  //1-- middle part is shifted, the others not
  //2-- higher end is shifted, the others not

  if(fabs(co1-co2)<5) {
    cout << "cutoff values must be different!" << endl;
    exit(1);
  }
  if(!(option==0 || option==1 || option==2)) {
    cout << "splitin3TH1D: Invalid option -- use 0, 1 or 2" << endl;
    exit(1);
  } 

  double cutoff1, cutoff2;
  if(co1 < co2) {
    cutoff1=co1;
    cutoff2=co2;
  }
  else if(co1 > co2) {
    cutoff1=co2;
    cutoff2=co1;
  }

  if(option==0) {
    TH1D *split=splitTH1D(orig_nom,orig_shift,cutoff1,0);
    return split;
  }
  else if(option==2) {
    TH1D *split=splitTH1D(orig_nom,orig_shift,cutoff2,1);
    return split;
  }
  else {//option==1
    string hname(orig_shift->GetName());
    TH1D *split=(TH1D*)orig_nom->Clone();

    split->Sumw2();
    bool lower=false, higher=false;

    for (int ibin=0; ibin< orig_shift->GetNbinsX()+2; ibin++){
      double shiftvalue=orig_shift->GetBinContent(ibin);
      double shifterror=orig_shift->GetBinError(ibin);
      double nomvalue=orig_nom->GetBinContent(ibin);
      double nomerror=orig_nom->GetBinError(ibin);
      
      double center=orig_shift->GetBinCenter(ibin);
      //cout << "bin " << ibin
      if(center<=cutoff1 || center>cutoff2) {
	split->SetBinContent(ibin, nomvalue);
	split->SetBinError(ibin, nomerror);
      }
      else if(center>cutoff1 && center<=cutoff2){ 
	split->SetBinContent(ibin, shiftvalue);
	split->SetBinError(ibin, shifterror);
      }
    }
    return split;
  }
}

//___________________________________________________________________


int main(int argc, char const *argv[]) {

  const bool doDebug=false;
  //const bool doDebug=false;
  //const bool doPartial=true; //saves no nominal histograms
  const bool doPartial=false;
  
  cout<<"doData=false !!!!!!!!!!!!!!!!! \n";
  const bool doData=false;
  //const bool doData=true;

  const bool doSetFirstBinsToZero=true;
  const bool DebugSetFirstBinsToZero=false; //true;

  const double cutoff1=800; //values at which we split the BoostJES systematic
  const double cutoff2=1300; 
  //0-800, 800-1300, 1300+

  TH1::SetDefaultSumw2(true); 

  const int Nbinnings=1; //2; 
  string binname[Nbinnings];
  binname[0]="binkkg"; //binname[1]="bin40"; binname[2]="bin100";
  //binname[1]="bin80";

  //string pdfbin="binkkg"; //binning of the PDF syst

  //signal cross sections and sample names
 
  

  vector<string>  m_SignalDample;
  m_SignalDample.push_back("Z400");
  m_SignalDample.push_back("Z500");
  m_SignalDample.push_back("Z750");
  m_SignalDample.push_back("Z1000");
  m_SignalDample.push_back("Z1250");
  m_SignalDample.push_back("Z1500");
  m_SignalDample.push_back("Z1750");
  m_SignalDample.push_back("Z2000");
  m_SignalDample.push_back("Z2250");
  m_SignalDample.push_back("Z2500");
  m_SignalDample.push_back("Z3000");

  /*m_SignalDample.push_back("KKg400");
  m_SignalDample.push_back("KKg500");
  m_SignalDample.push_back("KKg600");
  m_SignalDample.push_back("KKg700");
  m_SignalDample.push_back("KKg800");
  m_SignalDample.push_back("KKg900");
  m_SignalDample.push_back("KKg1000");
  m_SignalDample.push_back("KKg1150");
  m_SignalDample.push_back("KKg1300");
  m_SignalDample.push_back("KKg1600");
  m_SignalDample.push_back("KKg1800");
  m_SignalDample.push_back("KKg2000");
  m_SignalDample.push_back("KKg2250");
  m_SignalDample.push_back("KKg2500");
  m_SignalDample.push_back("KKg2750");
  m_SignalDample.push_back("KKg3000");
  
  m_SignalDample.push_back("KKg1000_width10pc");
  m_SignalDample.push_back("KKg1000_width15pc");
  m_SignalDample.push_back("KKg1000_width20pc");
  m_SignalDample.push_back("KKg1000_width25pc");
  m_SignalDample.push_back("KKg1000_width30pc");
  m_SignalDample.push_back("KKg1000_width35pc");
  m_SignalDample.push_back("KKg1000_width40pc");

  m_SignalDample.push_back("KKg2000_width10pc");
  m_SignalDample.push_back("KKg2000_width15pc");
  m_SignalDample.push_back("KKg2000_width20pc");
  m_SignalDample.push_back("KKg2000_width25pc");
  m_SignalDample.push_back("KKg2000_width30pc");
  m_SignalDample.push_back("KKg2000_width35pc");
  m_SignalDample.push_back("KKg2000_width40pc");

  m_SignalDample.push_back("KKg3000_width10pc");
  m_SignalDample.push_back("KKg3000_width15pc");
  m_SignalDample.push_back("KKg3000_width20pc");
  m_SignalDample.push_back("KKg3000_width25pc");
  m_SignalDample.push_back("KKg3000_width30pc");
  m_SignalDample.push_back("KKg3000_width35pc");
  m_SignalDample.push_back("KKg3000_width40pc");

  m_SignalDample.push_back("RSG400");
  m_SignalDample.push_back("RSG500");
  m_SignalDample.push_back("RSG600");
  m_SignalDample.push_back("RSG700");
  m_SignalDample.push_back("RSG800");
  m_SignalDample.push_back("RSG900");
  m_SignalDample.push_back("RSG1000");
  m_SignalDample.push_back("RSG1200");
  m_SignalDample.push_back("RSG1400");
  m_SignalDample.push_back("RSG1600");
  m_SignalDample.push_back("RSG1800");
  m_SignalDample.push_back("RSG2000");
  m_SignalDample.push_back("RSG2500");
  */

  vector<string> m_SignalIn=m_SignalDample;

  int NSAMPLE=m_SignalDample.size();
  double m_AMIXSection[NSAMPLE];

  //kkg1000,2000 widths
  for(int nsamp=0; nsamp<NSAMPLE; nsamp++)
    {
      m_AMIXSection[nsamp]=1.0;
    }

  /*m_SignalDample[23] = "KKg1000Width10pc";
  m_SignalDample[24] = "KKg1000Width15pc";
  m_SignalDample[25] = "KKg1000Width20pc";
  m_SignalDample[26] = "KKg1000Width25pc";
  m_SignalDample[27] = "KKg1000Width30pc";
  m_SignalDample[28] = "KKg1000Width35pc";
  m_SignalDample[29] = "KKg1000Width40pc";

  m_SignalDample[30] = "KKg2000Width10pc";
  m_SignalDample[31] = "KKg2000Width15pc";
  m_SignalDample[32] = "KKg2000Width20pc";
  m_SignalDample[33] = "KKg2000Width25pc";
  m_SignalDample[34] = "KKg2000Width30pc";
  m_SignalDample[35] = "KKg2000Width35pc";
  m_SignalDample[36] = "KKg2000Width40pc";
  */

  //KK graviton samples
  // m_SignalDample[19] = "KKGrav500";
  // m_SignalDample[20] = "KKGrav600";
  // m_SignalDample[21] = "KKGrav700";
  // m_SignalDample[22] = "KKGrav800";
  // m_SignalDample[23] = "KKGrav900";
  // m_SignalDample[24] = "KKGrav1000";
  // m_SignalDample[25] = "KKGrav1150";
  // m_SignalDample[26] = "KKGrav1300";



  
  /*m_SignalIn[10] = "kkg500";
  m_SignalIn[11] = "kkg600";
  m_SignalIn[12] = "kkg700";
  m_SignalIn[13] = "kkg800";
  m_SignalIn[14] = "kkg900";
  m_SignalIn[15] = "kkg1000";
  m_SignalIn[16] = "kkg1150";
  m_SignalIn[17] = "kkg1300";
  m_SignalIn[18] = "kkg1600";
  m_SignalIn[19] = "kkg1800";
  m_SignalIn[20] = "kkg2000";
  m_SignalIn[21] = "kkg2250";
  m_SignalIn[22] = "kkg2500";
  */
  /*m_SignalIn[23] = "kkg1000width10pc";
  m_SignalIn[24] = "kkg1000width15pc";
  m_SignalIn[25] = "kkg1000width20pc";
  m_SignalIn[26] = "kkg1000width25pc";
  m_SignalIn[27] = "kkg1000width30pc";
  m_SignalIn[28] = "kkg1000width35pc";
  m_SignalIn[29] = "kkg1000width40pc";

  m_SignalIn[30] = "kkg2000width10pc";
  m_SignalIn[31] = "kkg2000width15pc";
  m_SignalIn[32] = "kkg2000width20pc";
  m_SignalIn[33] = "kkg2000width25pc";
  m_SignalIn[34] = "kkg2000width30pc";
  m_SignalIn[35] = "kkg2000width35pc";
  m_SignalIn[36] = "kkg2000width40pc";
  */
  
  // m_SignalIn[19] = "KKGrav500";
  // m_SignalIn[20] = "KKGrav600";
  // m_SignalIn[21] = "KKGrav700";
  // m_SignalIn[22] = "KKGrav800";
  // m_SignalIn[23] = "KKGrav900";
  // m_SignalIn[24] = "KKGrav1000";
  // m_SignalIn[25] = "KKGrav1150";
  // m_SignalIn[26] = "KKGrav1300";
  
  /*
  //test set-up for speed
  m_SignalIn[3] = "NONE"; //"Zprime800";
  m_SignalIn[4] = "NONE"; //"Zprime1000";
  m_SignalIn[5] = "NONE"; //"Zprime1300";
  m_SignalIn[6] = "NONE"; //"Zprime1600";
  m_SignalIn[7] = "NONE"; //"Zprime2000";
  m_SignalIn[8] = "NONE"; //"Zprime2500";
  m_SignalIn[9] = "NONE"; //"Zprime3000";
  
  m_SignalIn[10] = "NONE"; //"KKG700";
  m_SignalIn[11] = "NONE"; //"KKG800";
  m_SignalIn[12] = "NONE"; //"KKG900";
  m_SignalIn[13] = "NONE"; //"KKG1000";
  m_SignalIn[14] = "NONE"; //"KKG1150";
  m_SignalIn[15] = "NONE"; //"KKG1300";
  m_SignalIn[16] = "NONE"; //"KKG1600";
  m_SignalIn[17] = "NONE"; //"KKg1800";
  m_SignalIn[18] = "NONE"; //"KKG2000";
  */


  int isample =0;

  // mass reconstruction schemes
  /*  const int nreco=3; //1; //2; //3;
  string reco[nreco];
  reco[0]="masstT"; // boosted
  // dRmin
  reco[1]="massTTbarChi2LPC"; // chi2
  //reco[0]="massTTbarChi2LPC"; // chi2
  reco[2]="massTTbarDRMin"; 
  */

  const int nreco=6; //1; //2; //3;
  string reco[nreco];
  reco[0]="masstT_cat1"; // boosted
  reco[1]="masstT_cat2"; // boosted
  reco[2]="masstT_cat3"; // boosted
  //reco[3]="";//"masstT_all"; // boosted
  // dRmin
  reco[3]="massTTbarChi2LPC_cat1"; // chi2
  reco[4]="massTTbarChi2LPC_cat2"; // chi2
  reco[5]="massTTbarChi2LPC_cat3"; // chi2
  //reco[6]="masstT";
  //reco[7]="massTTbarChi2LPC"; // chi2



  const double BinsZeroThreshold[nreco]={400.,400.,400., 240., 240., 240.};
  //const double BinsZeroThreshold[nreco]={400.,400.,400., 240., 240., 240.,400., 240.};
  //const double ThesholdResolved=240.;
  //const double ThesholdBoosted=400;


  int ireco =0;

  // input and output file
  TFile* infile;
  TFile* outfile;

  
  
  //string infilename="Hadronic/outputs_Jiahang/resolvedPreSpectra.root";
  //string outname="Hadronic/outputs_Jiahang/resolvedSpectra";

  string infilename="outputs_april22/lepPreSpectra_EWS.root";
  string outname="outputs_april22/topStatSpectra_EWS_noQCD_minimal";


  infile=TFile::Open(infilename);
  
  //string outname="test";
  if(doPartial) outname+="_partial";
  outfile=TFile::Open(outname+".root","recreate");  

  //  TSysLimit limit(0);
  std::cout << " start of loop " << std::endl;
 for (isample=0; isample < NSAMPLE ; isample++){
   // for (isample=0; isample < 1 ; isample++){
    std::cout << " Signal sample " << m_SignalDample[isample]<< std::endl;
    if (m_SignalDample[isample] == "NONE" || m_SignalIn[isample] == "NONE") continue;
    for (int ichannel=0; ichannel<2;ichannel++){
      string channelname, channel_in;
      if (ichannel==0) {
	channelname="e"; channel_in="e"; }
      if (ichannel==1) {
	channelname="mu"; channel_in="mu"; }

      for (ireco=0; ireco<nreco; ireco++){
	if(reco[ireco]=="") continue;
	string prefix_histname="";
	string suffix="_"+reco[ireco]+"_"+channel_in;
	string this_histname;

	string up="_up", dw="_dw";
	string upIn="_up", dwIn="_dw";

	TH1D *histBgr;
	TH1D *histBgrtmp_up;
	TH1D *histBgrtmp_dw;
	TH1D *histBgrtmp2;
	TH1D *histBgrtmpd;
	TH1D *histSignal;
	TH1D *histSignaltmp_up;
	TH1D *histSignaltmp_dw;
	TH1D *histSum;
	TH1D *htmp;
	TH1D *hbgr;

	//const int NSYS=(doPartial ? 13: 20);
	//const int NSYS=(doPartial ? 60: 60);
	//const int NSYS=(doPartial ? 61: 61);
	//const int NSYS= 74;
 	//string name[NSYS];
 	//string nameIn[NSYS];
	vector<string> name;
	vector<string> nameIn;

	/*for(int i=0; i< NSYS; i++)
	  {
	    name[i]="";
	    nameIn[i]="";
	    }
	*/


	name.push_back("EleSF");
	name.push_back("MuSF");

	name.push_back("Btag");
	name.push_back("BtagC");
	name.push_back("BtagL");

	name.push_back("MCGen");
	name.push_back("PartonShower");
	name.push_back("EWS");
	name.push_back("IFSR");
	name.push_back("topmass");
	name.push_back("JetEnerRes");
	
	name.push_back("JES_ALL");
	name.push_back("smallJES");
	
	name.push_back("BoostedJES");
	name.push_back("BoostedJES0");
	name.push_back("BoostedJES1");
	name.push_back("BoostedJES2");
	
	name.push_back("norm_tt");
	name.push_back("luminosity");
	
	
	/*name.push_back("Btag0");
	name.push_back("Btag1");
	name.push_back("Btag2");
	name.push_back("Btag3");
	name.push_back("Btag4");
	name.push_back("Btag5");
	name.push_back("Btag6");
	name.push_back("Btag7");
	name.push_back("Btag8");
	name.push_back("Btag9");
	name.push_back("Btag10");
	name.push_back("BtagC0");
	name.push_back("BtagC1");
	name.push_back("BtagC2");
	name.push_back("BtagC3");
	name.push_back("BtagC4");
	name.push_back("BtagC5");
	name.push_back("BtagC6");
	name.push_back("BtagL0");
	name.push_back("BtagL1");
	name.push_back("BtagL2");
	name.push_back("BtagL3");
	name.push_back("BtagL4");
	name.push_back("BtagL5");
	name.push_back("BtagL6");
	name.push_back("BtagL7");
	name.push_back("BtagL8");
	name.push_back("BtagL9");
	name.push_back("BtagL10");
	name.push_back("BtagL11");
	name.push_back("BtagL12");

	
	name.push_back("JES0");
	name.push_back("JES1");
	name.push_back("JES2");
	name.push_back("JES3");
	name.push_back("JES4");
	name.push_back("JES5");
	name.push_back("JES6");
	name.push_back("JES7");
	name.push_back("JES8");
	name.push_back("JES9");
	name.push_back("JES10");
	name.push_back("JES11");
	name.push_back("JES12");
	name.push_back("JES13");
	//name.push_back("JES14");
	name.push_back("JES15");
	name.push_back("JES16");
	name.push_back("JES17");
	name.push_back("JES18");
	//name.push_back("JES19");
	name.push_back("JES20");
	name.push_back("JES21");
	name.push_back("JES22");
	*/
	

	
	
	nameIn=name;
	
	int NSYS= name.size();
       
	for(int ns=0; ns<NSYS; ns++)
	  {
	    if(nameIn[ns]=="JetEnerRes") nameIn[ns]="JetEnerRes_up";
	    if(nameIn[ns]=="IFSR") nameIn[ns]="IFSR_up";
	    if(nameIn[ns]=="MCGen") nameIn[ns]="MCGen_up";
	    if(nameIn[ns]=="PartonShower") nameIn[ns]="PartonShower_up";
	    if(nameIn[ns]=="topmass") nameIn[ns]="topmass_up";
	    if(nameIn[ns]=="norm_tt") nameIn[ns]="ttbarnorm";
	  }
	
	
	//name.push_back(10"JVF");	          nameIn.push_back(10"JVF");
	//name.push_back(12"norm_QCDe");	  nameIn.push_back(12"norm_QCD");
	//name.push_back(18"PDF");	          nameIn.push_back(18"pdfrw");
        //name.push_back(19"norm_QCDmu");	  nameIn.push_back(19"norm_QCD");
	//name.push_back(21"NormW");	  nameIn.push_back(21"Wjets5Syst");
	//name.push_back(30"WHFC0");	  nameIn.push_back(30"Wjets0Syst");
	//name.push_back(31"WHFC3");	  nameIn.push_back(31"Wjets3Syst");
	//name.push_back(32"WHFC4");	  nameIn.push_back(32"Wjets4Syst");
	



	//pure scaling systematics
	/*
	if(!doPartial) {
	  const int iNP=2;
	  name[iNP+0]="luminosity";     nameIn[iNP+0]="luminosity";
	  name[iNP+1]="norm_tt";        nameIn[iNP+1]="norm_tt";
	  name[iNP+2]="norm_W";         nameIn[iNP+2]="norm_W";
	  name[iNP+3]="norm_QCD";       nameIn[iNP+3]="norm_QCD";
	  name[iNP+4]="norm_Z";         nameIn[iNP+4]="norm_Z";
	  name[iNP+5]="norm_single-top";nameIn[iNP+5]="norm_single-top";
	  name[iNP+6]="norm_Diboson";   nameIn[iNP+6]="norm_Diboson";
	  //name[iNP+7]="scale_e";        nameIn[iNP+7]="scale_e"; //Trigger + Reco 
	  //name[iNP+8]="scale_mu";       nameIn[iNP+8]="scale_mu"; // Reco
	}
	*/
	TH1D *histBgrSys[NSYS][2];
	TH1D *histSignalSys[NSYS][2];
	TH1D *histData;


	if(doDebug) cout << NSYS << " systematics.\nLast systematic is " << name[NSYS-1] << endl;
        

	//string channelname;
	//const int NBGR=7;
	const int NBGR=8;
	string bgrname[NBGR];
	string bgrIn[NBGR];
		
	for(int z=0; z<NBGR; z++){
	  bgrname[z]="";
	  bgrIn[z]="";
	}
	  


	bgrname[0]="tt";
	bgrname[1]="W";
	//bgrname[2]="WHF";
	bgrname[2]="Z";   
	//bgrname[3]="Zbb";   
	bgrname[3]="single-top";   
	bgrname[4]="Diboson";  
	bgrname[5]=""; //"QCDe";  
	bgrname[6]=""; //"QCDmu";  
	bgrname[7]="ttV";  
	


      	bgrIn[0]="tt";
	bgrIn[1]="W";
	//bgrIn[2]="WHFjets";
	bgrIn[2]="Z";   
	//bgrIn[3]="Zbbjets";   
	bgrIn[3]="single-top";   
	bgrIn[4]="Diboson";  
	bgrIn[5]="";//"QCDe";
	bgrIn[6]="";//"QCDmu";
	bgrIn[7]="ttV";  

	string dataname="Data";

	TH1D *histBgrtmp[NBGR];
	TH1D *histSingleBgrSys[NSYS][2][NBGR];

	//Scale factors
	double Bgr_SF[NBGR][2];
	for(int i=0; i<NBGR;++i) {
	  for (int c=0;c<2;++c) {
	    Bgr_SF[i][c]=1.;
	    if(bgrname[i]=="W") {
	      if(c==0) Bgr_SF[i][c]=1; //0.5; 
	      else if(c==1) Bgr_SF[i][c]=1; //0.5; 
	    }
	    if(bgrname[i]=="tt") {
	      //Bgr_SF[i][c]=176.8/166.8; //Hathor 1.3 change
	      Bgr_SF[i][c]=1.;
	    }
	  }
	}

	// central histograms
	double scale;
	//if (ichannel==0) channelname="el";
	//if (ichannel==1) channelname="mu";

	//Zprime500_massTFrom_e
	this_histname=prefix_histname+m_SignalIn[isample]+suffix;
	
	//std::cout << "channel: " << channelname << std::endl;
	//std::cout << this_histname  << std::endl;
	histSignal= (TH1D*) infile->Get(this_histname)->Clone(m_SignalDample[isample]+"_"+reco[ireco]+"_"+channelname);
	scale =  1./ m_AMIXSection[isample]; 
	histSignal->Scale(scale);
	if (isample==0){
	  //Backgrounds
	   for (int ibgr=0; ibgr<NBGR; ibgr++){
	     if(bgrIn[ibgr]=="") continue;
	   this_histname=prefix_histname+bgrIn[ibgr]+suffix;
	   //std::cout << this_histname << std::endl;
	    if (ibgr==0){
	      
	      histBgrtmp2= (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname);
	      if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
	      histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
	      //cout<<"andrew : "<<histBgrtmp2->Integral()<<endl;
	      histBgrtmp[ibgr]=histBgrtmp2;
	      histBgr= (TH1D*) infile->Get(this_histname)->Clone("Bgr_"+reco[ireco]+"_"+channelname);
	      histBgr->Scale(Bgr_SF[ibgr][ichannel]);	      
	      //if(doDebug) cout << "Added bgr integral " << histBgr->Integral(-1,9999) << endl;
	    }
	    else{
	      histBgrtmp2=  (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname);
	      if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
	      histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
	      histBgrtmp[ibgr]=histBgrtmp2;
	      histBgr->Add(histBgrtmp2);
	      //if(doDebug) cout << "Added bgr integral " << histBgr->Integral(-1,9999) << endl;
	    }
	   }
	}
	std::cout << "background done" << std::endl;

	// systematics
	std::cout << " now load systematics " << std::endl;
	for(int s=0; s< NSYS; s++){
	  //std::cout<<"sys "<<s<<endl;
	  if(name[s]=="") continue;
	  this_histname=prefix_histname+m_SignalIn[isample];
	  std::cout << " Name: " << name[s] << std::endl;
	  ////////-------------------------
	  /*
	  if(name[s]=="BoostedJESlow" || name[s]=="BoostedJESmid"|| name[s]=="BoostedJEShigh") {
	    int splitoption=-1;
	    if(name[s].Contains("low")) {
	      splitoption=0;
	      if(doDebug) cout << "BoostedJES split low: mtt <" << cutoff1 << " shifted." << endl;
	    }
	    else if(name[s].Contains("mid")) {
	      splitoption=1;
	      if(doDebug) cout << "BoostedJES split mid: " << cutoff1 << " < mtt <" << cutoff2 << " shifted." << endl;
	    }
	    else if(name[s].Contains("high")) {
	      splitoption=2;
	      if(doDebug) cout << "BoostedJES split high: mtt >" << cutoff2 << " shifted." << endl;
	    }
	    this_histname=prefix_histname+m_SignalIn[isample]+suffix;
	    string this_dw="_"+nameIn[s]+dwIn, this_up="_"+nameIn[s]+upIn;
	    if(doDebug) std::cout << " " << this_histname+this_up << std::endl;
	    htmp=(TH1D*) infile->Get(this_histname+this_up)->Clone(m_SignalDample[isample]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
	    histSignaltmp_up= splitin3TH1D(histSignal, htmp, cutoff1, cutoff2, splitoption);
	    histSignaltmp_up->SetName(htmp->GetName());
	    cout << htmp->GetName() << endl;
	    histSignaltmp_up->Scale(scale);
	    //	      histSignaltmp_up->Print();
	    if(doDebug) std::cout << this_histname+this_dw << std::endl;
	    htmp=(TH1D*)infile->Get(this_histname+this_dw)->Clone(m_SignalDample[isample]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
	    histSignaltmp_dw= splitin3TH1D(histSignal, htmp, cutoff1, cutoff2, splitoption);
	    histSignaltmp_dw->SetName(htmp->GetName());
	    histSignaltmp_dw->Scale(scale);
	    //	      histSignaltmp_dw->Print();

	    if (isample==0){
	      for (int ibgr=0; ibgr<NBGR; ibgr++){
		this_histname=prefix_histname+bgrIn[ibgr]+suffix;
		if (ibgr==0){
		  if(doDebug) std::cout << this_histname+this_up << std::endl;
		  htmp=  (TH1D*) infile->Get(this_histname+this_up)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		  if(doDebug) std::cout << this_histname << std::endl;
		  hbgr=  (TH1D*) infile->Get(this_histname)->Clone("tmp");
		  if(doDebug) std::cout << "SF " << Bgr_SF[ibgr][ichannel] << std::endl;
		  histBgrtmp2=splitin3TH1D(hbgr, htmp, cutoff1, cutoff2, splitoption);
		  histBgrtmp2->SetName(htmp->GetName());
		  histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		  histSingleBgrSys[s][0][ibgr]=histBgrtmp2;

		  htmp= (TH1D*) infile->Get(this_histname+this_up)->Clone("Bgr_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		  histBgrtmp_up=splitin3TH1D(hbgr, htmp, cutoff1, cutoff2, splitoption);
		  histBgrtmp_up->SetName(htmp->GetName());
		  histBgrtmp_up->Scale(Bgr_SF[ibgr][ichannel]);

		  std::cout << this_histname+this_dw << std::endl;
		  htmp=  (TH1D*) infile->Get(this_histname+this_dw)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		  histBgrtmp2=splitin3TH1D(hbgr, htmp, cutoff1, cutoff2, splitoption);
		  histBgrtmp2->SetName(htmp->GetName());
		  histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		  //histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
		  histSingleBgrSys[s][1][ibgr]=histBgrtmp2;

		  htmp= (TH1D*) infile->Get(this_histname+this_dw)->Clone("Bgr_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		  histBgrtmp_dw=splitin3TH1D(hbgr, htmp, cutoff1, cutoff2, splitoption);
		  histBgrtmp_dw->SetName(htmp->GetName());
		  histBgrtmp_dw->Scale(Bgr_SF[ibgr][ichannel]);
		  if(doDebug) cout << "Added bgr integrals (up, dw) " << histBgrtmp_up->Integral(-1,9999) << ", " << histBgrtmp_dw->Integral(-1,9999) << endl;
		}
		else {
		  if (!(bgrname[ibgr]=="QCDe") && !(bgrname[ibgr]=="QCDmu")&& !(bgrname[ibgr]=="QCD")){
		    //this_histname=prefix_histname+bgrIn[ibgr]+"_"+nameIn[s];
		    std::cout << this_histname << std::endl;
		    hbgr=  (TH1D*) infile->Get(this_histname)->Clone("tmp");
		    this_histname=prefix_histname+bgrIn[ibgr]+suffix;
		    std::cout << this_histname+this_up << std::endl;
		    htmp=  (TH1D*) infile->Get(this_histname+this_up)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		    histBgrtmp2=splitin3TH1D(hbgr, htmp, cutoff1, cutoff2, splitoption);
		    histBgrtmp2->SetName(htmp->GetName());
if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
		    histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		    histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
		    histBgrtmp_up->Add(histBgrtmp2);

		    std::cout << this_histname+this_dw << std::endl;
		    htmp=  (TH1D*) infile->Get(this_histname+this_dw)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		    histBgrtmp2=splitin3TH1D(hbgr, htmp, cutoff1, cutoff2, splitoption);
		    histBgrtmp2->SetName(htmp->GetName());
		    if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
		    histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		    histSingleBgrSys[s][1][ibgr]=histBgrtmp2;
		    histBgrtmp_dw->Add(histBgrtmp2);
		    if(doDebug) cout << "Added bgr integrals (up, dw) " << histBgrtmp_up->Integral(-1,9999) << ", " << histBgrtmp_dw->Integral(-1,9999) << endl;
		  }
		  else{//if (bgrname[ibgr]=="QCD")
		    this_histname=prefix_histname+bgrIn[ibgr]+suffix;
		    std::cout << this_histname << std::endl;
		    histBgrtmp2=  (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up );
		    if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
		    histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		    histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
		    histBgrtmp_up->Add(histBgrtmp2);
		    histBgrtmp2=  (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		    if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
		    histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		    histSingleBgrSys[s][1][ibgr]=histBgrtmp2;
		    histBgrtmp_dw->Add(histBgrtmp2);
		    if(doDebug) cout << "Added bgr integrals (up, dw) " << histBgrtmp_up->Integral(-1,9999) << ", " << histBgrtmp_dw->Integral(-1,9999) << endl;
		  }
		}
	      }
	    }
	  }
	  */
	  ////////-------------------------
	  //else if
	  //name[s].SubString("JES")!=0
	  if( name[s].SubString("JES")!="" || name[s]=="BoostedJES0" || name[s]=="BoostedJES1"|| name[s]=="BoostedJES2" || name[s]=="Btag"|| name[s]=="BtagC"|| name[s]=="BtagL"|| name[s]=="JES" || name[s]=="JES0" || name[s]=="JES1" || name[s]=="JES2" || name[s]=="JES3" || name[s]=="JES4" || name[s]=="JES5" || name[s]=="JES6" || name[s]=="JES7" || name[s]=="JES8" || name[s]=="JES9" || name[s]=="JESP" || name[s]=="JESm" || name[s]=="JESb" || name[s]=="JESCloseby" || name[s]=="JESFlavComp" || name[s]=="JESFlavResp" || name[s]=="BoostedJES" || name[s]=="ScaleTtbar" || name[s]=="MetCellOut" || name[s]=="MetPileUp" || name[s]=="Mureco" || name[s]=="MuResId" || name[s]=="MuResMs" || name[s]=="EnerScale" || name[s]=="EleReco" || name[s]=="EleEnerRes" || name[s]=="EleID" || name[s]=="EleSF" || name[s]=="EleTrig" || name[s]=="LArHole" || name[s]=="MuTrig" || name[s]=="MuTrigMatch" || name[s]=="MuID" || name[s]=="Mureco" || name[s]=="MuEnerScale" || name[s]=="MuSF" || name[s]=="PU" || name[s]=="EWS"|| name[s]=="WHFC0"|| name[s]=="WHFC1"|| name[s]=="WHFC2"|| name[s]=="WHFC3"|| name[s]=="WHFC4"|| name[s]=="NormW"|| name[s]=="MI" ||name[s]=="JVF"|| name[s]=="NormZ"||name[s]=="MetScale"||name[s]=="MetReso" ||name[s]=="norm_W" || name[s].Contains("Btag")) {//name[s]=="norm_W" || name[s]=="norm_Z" || 	
	    //if (!((channelname=="e"  && (name[s]=="MuResId" || name[s]=="MuResMs" || name[s]=="Mureco"|| name[s]=="MuTrig" || name[s]=="MuID" || name[s]=="MuEnerScale")) || (channelname=="mu" && (name[s]=="EnerScale" || name[s]=="EleEnerRes" || name[s]=="EleID")))) { 
	    cout<<"name[s].SubString(JES) _"<<name[s].SubString("JES")<<"_"<<endl;
	    if(true) {
	      // double sided
	      
	      this_histname=prefix_histname+m_SignalIn[isample]+suffix;
	      string this_dw="_"+nameIn[s]+dwIn, this_up="_"+nameIn[s]+upIn;
	      if(doDebug) std::cout << " " << this_histname+this_up << std::endl;
 	     //  if (name[s]=="WHFC0" || name[s]=="WHFC3" || name[s]=="WHFC4") histSignaltmp_up=(TH1D*) infile->Get(this_histname)->Clone(m_SignalDample[isample]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);	
//   	      else
		
	      histSignaltmp_up=(TH1D*) infile->Get(this_histname+this_up)->Clone(m_SignalDample[isample]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
	      histSignaltmp_up->Scale(scale);
	      //	      histSignaltmp_up->Print();
	      if(doDebug) std::cout << this_histname+this_dw << std::endl;
	    
 	      // if (name[s]=="WHFC0" || name[s]=="WHFC3" || name[s]=="WHFC4") histSignaltmp_dw=(TH1D*) infile->Get(this_histname)->Clone(m_SignalDample[isample]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
//                 else 
		  
	      histSignaltmp_dw=(TH1D*)infile->Get(this_histname+this_dw)->Clone(m_SignalDample[isample]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
	      histSignaltmp_dw->Scale(scale);     

	      //	      histSignaltmp_dw->Print();
	      if (isample==0){
		for (int ibgr=0; ibgr<NBGR; ibgr++){
		  if(bgrIn[ibgr]=="") continue;
		  //this_histname=prefix_histname+bgrIn[ibgr]+"_"+nameIn[s];
		  this_histname=prefix_histname+bgrIn[ibgr]+suffix;
		  //if (!( name[s]=="IFSR"|| name[s]=="FSR" ||name[s]=="ISR")){
		  if(true) {
		    if (ibgr==0){//ttbar background
		      //std::cout << this_histname+this_up << std::endl;		      
		      // if (name[s]=="QCDMMFake"|| name[s]=="QCDMMEff")  histBgrtmp2=(TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up); 
// 		      else 
			histBgrtmp2=  (TH1D*) infile->Get(this_histname+this_up)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		      if(doDebug) std::cout << "SF " << Bgr_SF[ibgr][ichannel] << std::endl;
		      histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		      histSingleBgrSys[s][0][ibgr]=histBgrtmp2;		      
                    //  if (name[s]=="QCDMMFake"|| name[s]=="QCDMMEff")  histBgrtmp_up= (TH1D*) infile->Get(this_histname)->Clone("Bgr_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
// 		     else 
		       histBgrtmp_up= (TH1D*) infile->Get(this_histname+this_up)->Clone("Bgr_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		      histBgrtmp_up->Scale(Bgr_SF[ibgr][ichannel]);
		      //std::cout << this_histname+this_dw << std::endl;		     
		     //  if (name[s]=="QCDMMFake"|| name[s]=="QCDMMEff") histBgrtmp2=  (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
// 		      else  
			histBgrtmp2=  (TH1D*) infile->Get(this_histname+this_dw)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		      histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		      //histSingleBgrSys[s][0][ibgr]=histBgrtmp2;

		      histSingleBgrSys[s][1][ibgr]=histBgrtmp2;
		     
		     //  if (name[s]=="QCDMMFake"|| name[s]=="QCDMMEff") histBgrtmp_dw= (TH1D*) infile->Get(this_histname)->Clone("Bgr_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
// 		      else 
			histBgrtmp_dw= (TH1D*) infile->Get(this_histname+this_dw)->Clone("Bgr_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		      histBgrtmp_dw->Scale(Bgr_SF[ibgr][ichannel]);
		      if(doDebug) cout << "Added bgr integrals (up, dw) " << histBgrtmp_up->Integral(-1,9999) << ", " << histBgrtmp_dw->Integral(-1,9999) << endl;
		    }
		    else {
		      if (!(bgrname[ibgr]=="QCDe") && !(bgrname[ibgr]=="QCDmu") && !(bgrname[ibgr]=="QCD"))
			//|| (bgrname[ibgr]=="QCD" && (name[s]=="QCDMMFake" || name[s]=="QCDMMEff")
			{
			//this_histname=prefix_histname+bgrIn[ibgr]+"_"+nameIn[s];
			this_histname=prefix_histname+bgrIn[ibgr]+suffix;
			//std::cout << this_histname+this_up <<" take up/dw from origin "<< std::endl;
			histBgrtmp2=  (TH1D*) infile->Get(this_histname+this_up)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
			if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
			histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);

			histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
			histBgrtmp_up->Add(histBgrtmp2);
			//std::cout << this_histname+this_dw <<" take up/dw from origin "<< std::endl;
			histBgrtmp2=  (TH1D*) infile->Get(this_histname+this_dw)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
			if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
			histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
			histSingleBgrSys[s][1][ibgr]=histBgrtmp2;
			histBgrtmp_dw->Add(histBgrtmp2);
			if(doDebug) cout << "Added bgr integrals (up, dw) " << histBgrtmp_up->Integral(-1,9999) << ", " << histBgrtmp_dw->Integral(-1,9999) << endl;
		      }
		      else{//if (bgrname[ibgr]=="QCD")
			this_histname=prefix_histname+bgrIn[ibgr]+suffix;
			//std::cout << this_histname <<std::endl;
			//std::cout << bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up<<" take nominal from origin as up/dw"<< std::endl;
			histBgrtmp2=  (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up );
			if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
			histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
			histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
			histBgrtmp_up->Add(histBgrtmp2);
			histBgrtmp2=  (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
			if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
			histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
			histSingleBgrSys[s][1][ibgr]=histBgrtmp2;
			histBgrtmp_dw->Add(histBgrtmp2);
			if(doDebug) cout << "Added bgr integrals (up, dw) " << histBgrtmp_up->Integral(-1,9999) << ", " << histBgrtmp_dw->Integral(-1,9999) << endl;
		      }
		    }
		    // only signal double sidedsystematics
		  }
		  /*else if ( name[s]=="IFSR"|| name[s]=="FSR" ||name[s]=="ISR"){
		    this_histname=prefix_histname+bgrIn[ibgr]; //+"_"+nameIn[s];
		    if (ibgr==0){
		      histBgrtmp2=(TH1D*) infile->Get(this_histname+this_up)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		      if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
		      histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		      histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
		      histBgrtmp_up= (TH1D*) infile->Get(this_histname+this_up)->Clone("Bgr_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		      histBgrtmp_up->Scale(Bgr_SF[ibgr][ichannel]);

		      histBgrtmp2=(TH1D*) infile->Get(this_histname+this_dw)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		      histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		      histSingleBgrSys[s][1][ibgr]=histBgrtmp2;
		      histBgrtmp_dw= (TH1D*) infile->Get(this_histname+this_dw)->Clone("Bgr_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		      histBgrtmp_dw->Scale(Bgr_SF[ibgr][ichannel]);
		    }
		    else{
		      this_histname=prefix_histname+bgrIn[ibgr]+suffix;
		      histBgrtmp2=  (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		      if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
		      histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		      histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
		      histBgrtmp_up->Add(histBgrtmp2);
		      histBgrtmp2=  (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		      histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		      histSingleBgrSys[s][1][ibgr]=histBgrtmp2;
		      histBgrtmp_dw->Add(histBgrtmp2);
		    }
		  }*/	      
		}  
	      } 
	    }
	    /*else {
	      // for e or mu only
	      std::cout << " and not existing for this channel-> use dummy";
	      this_histname=prefix_histname+m_SignalIn[isample]+suffix;
	      histSignaltmp_up=(TH1D*) infile->Get(this_histname)->Clone(m_SignalDample[isample]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
	      histSignaltmp_up->Scale(scale);
	      histSignaltmp_dw=(TH1D*)infile->Get(this_histname)->Clone(m_SignalDample[isample]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
	      histSignaltmp_dw->Scale(scale);
	      if(isample==0){
		for (int ibgr=0; ibgr<NBGR; ibgr++){
		  this_histname=prefix_histname+bgrIn[ibgr]+suffix;
		  if (ibgr==0){
		    histBgrtmp2 = (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		    if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
		    histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		    histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
		    histBgrtmp_up= (TH1D*) infile->Get(this_histname)->Clone("Bgr_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		    histBgrtmp_up->Scale(Bgr_SF[ibgr][ichannel]);
		    histBgrtmp2= (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		    histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		    histSingleBgrSys[s][1][ibgr]=histBgrtmp2;
		    histBgrtmp_dw= (TH1D*) infile->Get(this_histname)->Clone("Bgr_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		    histBgrtmp_dw->Scale(Bgr_SF[ibgr][ichannel]);
		  }
		  else{
		    histBgrtmp2=  (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		    histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		    histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
		    histBgrtmp_up->Add(histBgrtmp2);
		    histBgrtmp2=  (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		    histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		    histSingleBgrSys[s][1][ibgr]=histBgrtmp2;
		    histBgrtmp_dw->Add(histBgrtmp2);
		  }
		}
	      }
	    }*/
	    std::cout <<  name[s] << " ----> done " << std::endl;
	  }
	  //else if ( name[s]=="PDF") {
	    //Special binning. Rebin QCD here.
	    // double sided
	  //}
	  else if ( name[s]=="JetEffMC" || name[s]=="PileUp" || name[s]=="Fragmentation" ||  name[s]=="JetEnerRes" || name[s]=="PDF" || name[s]=="WRWiqopt2" || name[s]=="WRWiqopt3" || name[s]=="WRWptjmin10" || name[s]=="WRWptjmin20" || name[s]=="WRWqfacKtfac" ||  name[s]=="shape_QCD" || name[s]=="MuScale" || name[s]=="JEE" ||name[s]=="Norm_Diboson" || name[s]=="Norm_single-top" || name[s]=="NormTtbar" ||  name[s]=="JER0" || name[s]=="JER1" || name[s]=="JER2"||name[s]=="BoostedTrig"||name[s]=="PartonShower"|| name[s]=="IFSR" ||name[s]=="topmass" ||name[s]=="MCGen"){//name[s]=="NormZ"  ||name[s]=="PartonShower"|| name[s]=="IFSR" || 
	    std::cout << " shape systematics " << name[s] <<  " single sided";
	    this_histname=prefix_histname+m_SignalIn[isample]+suffix;
	    //if(doDebug) std::cout << std::endl << " " << this_histname << std::endl;
	    if ( name[s]=="WRWiqopt2" ||name[s]=="WRWiqopt3" ||name[s]=="WRWptjmin10"||name[s]=="WRWptjmin20"){
	    //if(false) {
	      if(doDebug) std::cout << " " << nameIn[s] << std::endl;
	      if(doDebug) std::cout << this_histname+"_"+nameIn[s] << std::endl;
	      histSignaltmp_up=(TH1D*) infile->Get(this_histname+"_"+nameIn[s])->Clone(m_SignalDample[isample]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
	      histSignaltmp_up->Scale(scale);
	      
	      if(isample==0){
		for (int ibgr=0; ibgr<NBGR; ibgr++){
		  if(bgrIn[ibgr]=="") continue;
		  this_histname=prefix_histname+bgrIn[ibgr]+suffix;
		  
		  if (ibgr==0){//no shift for ttbar
		    if(doDebug) std::cout << this_histname+"_"+nameIn[s] << std::endl;
		    histBgrtmp_up= (TH1D*) infile->Get(this_histname+"_"+nameIn[s])->Clone("Bgr_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		    if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
		    histBgrtmp_up->Scale(Bgr_SF[ibgr][ichannel]);

		    histBgrtmp2 = (TH1D*) infile->Get(this_histname+"_"+nameIn[s])->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		    histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		    histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
		    if(doDebug) std::cout << this_histname << std::endl;
		    histBgrtmp_up= (TH1D*) infile->Get(this_histname)->Clone("Bgr_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		    histBgrtmp_up->Scale(Bgr_SF[ibgr][ichannel]);
		    if(doDebug) std::cout << this_histname << std::endl;
		    histBgrtmpd = (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		    histBgrtmpd->Scale(Bgr_SF[ibgr][ichannel]);
		    histBgrtmpd->Scale(2.);
		    histBgrtmpd->Add(histBgrtmp2,-1.);
		    histSingleBgrSys[s][1][ibgr]=histBgrtmpd; 
		  }
		  else if (bgrname[ibgr]=="W" ||bgrname[ibgr]=="WHF") { 
		    if(doDebug) std::cout << this_histname+"_"+nameIn[s] << std::endl;
		    histBgrtmp2=  (TH1D*) infile->Get(this_histname+"_"+nameIn[s])->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		    if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
		    histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);

		    //Normalise to nominal yield
		    if(doDebug) std::cout << this_histname << std::endl;
		    hbgr=  (TH1D*) infile->Get(this_histname)->Clone("tmp");
		    hbgr->Scale(Bgr_SF[ibgr][ichannel]);
		    double wscale_nom = hbgr->Integral(-1,9999);
		    if(doDebug) std::cout << "Nom integral " << wscale_nom << std::endl;
		    double wscale = wscale_nom/histBgrtmp2->Integral(-1,9999);
		    if(doDebug) std::cout << "Up integral " << histBgrtmp2->Integral(-1,9999) << std::endl;
		    if(doDebug) std::cout << "Norm scale (up): " << wscale << std::endl;
		    histBgrtmp2->Scale(wscale);

		    //Get dw histo
		    if(doDebug) std::cout << this_histname << std::endl;
		    histBgrtmpd = (TH1D*)hbgr->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		    if(doDebug) std::cout << "Before symmetrising: Dw integral " << histBgrtmpd->Integral(-1,9999) << std::endl;

		    //symmetrise
		    histBgrtmpd->Scale(2.);
		    histBgrtmpd->Add(histBgrtmp2,-1.);

		    //Cross check normalisation to nominal yield
		    if(doDebug) std::cout << "Crosscheck: Up integral " << histBgrtmp2->Integral(-1,9999) << std::endl;
		    if(doDebug) std::cout << "Crosscheck: Dw integral " << histBgrtmpd->Integral(-1,9999) << std::endl;

		    //add to hist vector and tot Bgr
		    histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
		    histBgrtmp_up->Add(histBgrtmp2);
		    histSingleBgrSys[s][1][ibgr]=histBgrtmpd;
		  }
		  else {//no shift for other backgrounds
		    if(doDebug) std::cout << this_histname<< std::endl;
		    histBgrtmp2=  (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		    if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
		    histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		    histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
		    //		      histBgrtmp2->Write();
		    histBgrtmp_up->Add(histBgrtmp2);
		    if(doDebug) std::cout << this_histname<< std::endl;
		    histBgrtmpd=  (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		    histBgrtmpd->Scale(Bgr_SF[ibgr][ichannel]);
		    histBgrtmpd->Scale(2.);
		    histBgrtmpd->Add(histBgrtmp2,-1.);
		    histSingleBgrSys[s][1][ibgr]=histBgrtmpd;
		  }
		}//end loop over backgrounds
	      }//end if(isample==0)
	    }//end if iqopt, ptjmin
	    else {//all other single-sided backgrounds
	      if(doDebug) std::cout << nameIn[s] << std::endl;
	      if(doDebug) std::cout << this_histname+"_"+nameIn[s] << std::endl;
	      histSignaltmp_up=(TH1D*) infile->Get(this_histname+"_"+nameIn[s])->Clone(m_SignalDample[isample]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
	      histSignaltmp_up->Scale(scale);
	      if(isample==0){
		for (int ibgr=0; ibgr<NBGR; ibgr++){
		  if(bgrIn[ibgr]=="") continue;
		  this_histname=prefix_histname+bgrIn[ibgr]+suffix;
		  if(true) {
		    this_histname=prefix_histname+bgrIn[ibgr]+suffix;
		    if (ibgr==0){
		      if(doDebug) std::cout << this_histname+"_"+nameIn[s] << std::endl;
		      histBgrtmp_up= (TH1D*) infile->Get(this_histname+"_"+nameIn[s])->Clone("Bgr_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		      if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
		      histBgrtmp_up->Scale(Bgr_SF[ibgr][ichannel]);
		      histBgrtmp2= (TH1D*) infile->Get(this_histname+"_"+nameIn[s])->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		      histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		      histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
		      if(doDebug) std::cout << this_histname << std::endl;
		      histBgrtmpd=  (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		      histBgrtmpd->Scale(Bgr_SF[ibgr][ichannel]);
		      histBgrtmpd->Scale(2.);
		      histBgrtmpd->Add(histBgrtmp2,-1.);
		      histSingleBgrSys[s][1][ibgr]=histBgrtmpd;
		    }//end if (ibgr==0)
		    else{ //ibgr>0
		      if (!(bgrname[ibgr]=="QCDe") && !(bgrname[ibgr]=="QCDmu")&& !(bgrname[ibgr]=="QCD")){
			if(doDebug) std::cout << this_histname+"_"+nameIn[s] << std::endl;
			histBgrtmp2=  (TH1D*) infile->Get(this_histname+"_"+nameIn[s])->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
			if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
			histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
			histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
			histBgrtmp_up->Add(histBgrtmp2);

			//down histo
			if(doDebug) std::cout << this_histname<< std::endl;
			histBgrtmpd=  (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
			histBgrtmpd->Scale(Bgr_SF[ibgr][ichannel]);
			histBgrtmpd->Scale(2.);
			histBgrtmpd->Add(histBgrtmp2,-1.);
			histSingleBgrSys[s][1][ibgr]=histBgrtmpd;
		      }
		      else{ //bgrname[ibgr]=="QCD"
			//No change for QCD
			if(doDebug) std::cout << this_histname<< std::endl;
			histBgrtmp2=  (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
			if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
			histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
			histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
			//		      histBgrtmp2->Write();
			histBgrtmp_up->Add(histBgrtmp2);
			if(doDebug) std::cout << this_histname<< std::endl;
			histBgrtmpd=  (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
			histBgrtmpd->Scale(Bgr_SF[ibgr][ichannel]);
			histBgrtmpd->Scale(2.);
			histBgrtmpd->Add(histBgrtmp2,-1.);
			histSingleBgrSys[s][1][ibgr]=histBgrtmpd;
			//		      histBgrtmpd->Write();
		      }//end else bgrname[ibgr]=="QCD"
		    } //end else ibgr>0		  
		  }//////end if(true)
		  else if(name[s]=="shape_QCD") {
		    this_histname=prefix_histname+bgrIn[ibgr]+suffix;
		    if (ibgr==0){
		      if(doDebug) std::cout << this_histname << std::endl;
		      histBgrtmp2 = (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		      if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
		      histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		      histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
		      if(doDebug) std::cout << this_histname << std::endl;
		      histBgrtmp_up= (TH1D*) infile->Get(this_histname)->Clone("Bgr_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		      histBgrtmp_up->Scale(Bgr_SF[ibgr][ichannel]);
		      if(doDebug) std::cout << this_histname << std::endl;
		      histBgrtmpd = (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		      histBgrtmpd->Scale(Bgr_SF[ibgr][ichannel]);
		      histBgrtmpd->Scale(2.);
		      histBgrtmpd->Add(histBgrtmp2,-1.);
		      histSingleBgrSys[s][1][ibgr]=histBgrtmpd; 
		    }
		    else if (bgrname[ibgr]== "QCDe" || bgrname[ibgr]== "QCDmu"|| bgrname[ibgr]== "QCD") { //change ONLY for QCD
		      if(doDebug) std::cout << this_histname+"_"+nameIn[s] << std::endl;
		      histBgrtmp2=  (TH1D*) infile->Get(this_histname+"_"+nameIn[s])->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		      if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
		      histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		      histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
		      histBgrtmp_up->Add(histBgrtmp2);
		      //down histo
		      if(doDebug) std::cout << this_histname << std::endl;
		      histBgrtmpd=  (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		      histBgrtmpd->Scale(Bgr_SF[ibgr][ichannel]);
		      histBgrtmpd->Scale(2.);
		      histBgrtmpd->Add(histBgrtmp2,-1.);
		      histSingleBgrSys[s][1][ibgr]=histBgrtmpd;
		    } 
		    else { //all other backgrounds
		      if(doDebug) std::cout << this_histname << std::endl;
		      histBgrtmp2=  (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		      if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
		      histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		      histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
		      histBgrtmp_up->Add(histBgrtmp2);
		      if(doDebug) std::cout << this_histname << std::endl;
		      histBgrtmpd = (TH1D*) infile->Get(this_histname)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		      histBgrtmpd->Scale(Bgr_SF[ibgr][ichannel]);
		      histBgrtmpd->Scale(2.);
		      histBgrtmpd->Add(histBgrtmp2,-1.);
		      histSingleBgrSys[s][1][ibgr]=histBgrtmpd;
		    }
		  }//if syst = shape_QCD
		}//end loop over backgrounds
	      }//end if(isample==0) 
	    }//all other single-sided backgrounds

	    //Symmetrise signal and tot background
	    histSignaltmp_dw= (TH1D*) histSignal->Clone(string(histSignal->GetName())+"_"+name[s]+dw);
	    histSignaltmp_dw->Scale(2.);
	    histSignaltmp_dw->Add(histSignaltmp_up,-1.);
	    if(isample==0) {
	      histBgrtmp_dw=(TH1D*) histBgr->Clone(string(histBgr->GetName())+"_"+name[s]+dw); 
	      histBgrtmp_dw->Scale(2.);
	      histBgrtmp_dw->Add(histBgrtmp_up,-1.);
	    }
	    std::cout <<  name[s] << " ----> done " << std::endl;
	  }
	  //else if ( name[s]=="luminosity" ||  name[s]=="scale_e" ||name[s]=="scale_mu"|| name[s]=="norm_tt" || name[s]=="norm_W" || name[s]=="norm_QCD" || name[s]=="norm_Z" || name[s]=="norm_single-top" || name[s]=="norm_Diboson"){
	  else if (name[s]=="luminosity" ||  name[s]=="scale_e" ||  name[s]=="scale_mu" || name[s]=="norm_tt" || name[s]=="norm_QCD" || name[s]=="norm_QCDe" || name[s]=="norm_QCDmu" || name[s]=="norm_Z" || name[s]=="norm_single-top" || name[s]=="norm_Diboson"){
	    // QCD and W+jets are data driven -> no systematics
	    //Lumi uncert 3.9% (Kevin Einsweiler, Sep 27 2012)
	    //Previous uncert: 1.8%
	    double lumisys_signal_e_up=0.028;
	    double lumisys_signal_e_dw=-0.028;
	    double lumisys_signal_mu_up=0.028;
	    double lumisys_signal_mu_dw=-0.028;
	                               //ttbar, W,  Z,     s-top, VV,    QCDe,  QCDmu, ttv
	    double lumisys_e_up[NBGR] ={ 0.028, 0., 0.028, 0.028, 0.028, 0.0, 0.0 , 0.28};
	    double lumisys_e_dw[NBGR] ={-0.028, 0.,-0.028,-0.028,-0.028, 0.0, 0.0 , -0.28};
	    double lumisys_mu_up[NBGR]={ 0.028, 0., 0.028, 0.028, 0.028, 0.0, 0.0 , 0.28};
	    double lumisys_mu_dw[NBGR]={-0.028, 0.,-0.028,-0.028,-0.028, 0.0, 0.0 , -0.28};	    
	    //Scal systematic
	    double scalesys_signal_e_up=0.018;
	    double scalesys_signal_e_dw=-0.018;
	    double scalesys_signal_mu_up=0.010;
	    double scalesys_signal_mu_dw=-0.010;
	                                //ttbar, W,  Z,     s-top, VV,    QCDe,QCDmu 
	    double scalesys_e_up[NBGR] ={ 0.018, 0., 0.018, 0.018, 0.018, 0.0, 0.0, 0.018};
	    double scalesys_e_dw[NBGR] ={-0.018, 0.,-0.018,-0.018,-0.018, 0.0, 0.0, -0.018};
	    double scalesys_mu_up[NBGR]={ 0.018, 0., 0.018, 0.018, 0.018, 0.0, 0.0, 0.018};
	    double scalesys_mu_dw[NBGR]={-0.018, 0.,-0.018,-0.018,-0.018, 0.0, 0.0, -0.018};

	    //Norm uncert
	                               //ttbar, W,  Z,    s-top, VV,    QCDe, QCDmu 
	    //double normsys_e_up[NBGR] ={ 0.098, 0., 0., 0.077, 0., 0.5,  0.5};
	    //double normsys_e_dw[NBGR] ={-0.106, 0.,-0.,-0.077,-0.,-0.5, -0.5};
	    //double normsys_mu_up[NBGR]={ 0.098, 0., 0., 0.077, 0., 0.5,  0.5};
	    //double normsys_mu_dw[NBGR]={-0.106, 0.,-0.,-0.077,-0.,-0.5, -0.5};
	    double normsys_e_up[NBGR] ={ .0606, 0., 0., 0.077, 0., 0.5,  0.5, 0};
	    double normsys_e_dw[NBGR] ={-0.0643, 0.,-0.,-0.077,-0.,-0.5, -0.5, 0};
	    double normsys_mu_up[NBGR]={ 0.0606, 0., 0., 0.077, 0., 0.5,  0.5, 0};
	    double normsys_mu_dw[NBGR]={-0.0643, 0.,-0.,-0.077,-0.,-0.5, -0.5, 0};

	    //Set values
	    double sys_up;
	    double sys_dw;

	    
	    if (name[s]=="luminosity"){
	      if (channelname=="e"){
		sys_up = lumisys_signal_e_up;  sys_dw = lumisys_signal_e_dw;}
	      else { 
		sys_up = lumisys_signal_mu_up; sys_dw = lumisys_signal_mu_dw;}
	    } 
	    else if (name[s]=="scale_e"){
	      if (channelname=="e") {
		sys_up = scalesys_signal_e_up;sys_dw = scalesys_signal_e_dw;}
	      else  {
		sys_up = 0.; sys_dw = 0.;}
	    } 
	    else if (name[s]=="scale_mu"){
	      if (channelname=="e") {
		sys_up = 0.;sys_dw = 0.;}
	      else {
		sys_up = scalesys_signal_mu_up; sys_dw = scalesys_signal_mu_dw;}
	    }  
	    else if (name[s]=="norm_tt" ||name[s]=="norm_W" || name[s]=="norm_QCD"|| name[s]=="norm_QCDe" || name[s]=="norm_QCDmu" ||name[s]=="norm_Z" ||name[s]=="norm_single-top"|| name[s]=="norm_Diboson"){
	      sys_up =0; sys_dw =0;
	    }
	    std::cout << " scale systematics: " << name[s] << std::endl;
	    //std::cout << "   " << histSignal->GetName() << "_" << name[s] << "_up" << std::endl;
	    histSignaltmp_up = (TH1D*) histSignal->Clone(string(histSignal->GetName())+"_"+name[s]+up);
	    histSignaltmp_up->Scale(1.+sys_up);
 	    //std::cout << "   " << histSignal->GetName() << "_" << name[s] << "_dw" << std::endl;
	    histSignaltmp_dw = (TH1D*) histSignal->Clone(string(histSignal->GetName())+"_"+name[s]+dw); 
	    histSignaltmp_dw->Scale(1.+sys_dw); 
	    for (int ibgr=0; ibgr<NBGR; ibgr++){
	      if(bgrIn[ibgr]=="") continue;
	      sys_up = 0;sys_dw = 0;//bugfix from andrew

	      //cout<<"andrew: "<< name[s]<<" "<<bgrname[ibgr]<<" "<<channelname<<endl;
	      if (name[s]=="luminosity"){
		if (channelname=="e") {
		  sys_up = lumisys_e_up[ibgr]; sys_dw = lumisys_e_dw[ibgr];
		}
		else {
		  sys_up = lumisys_mu_up[ibgr]; sys_dw = lumisys_mu_dw[ibgr];
		}
	      } 
	      else if (name[s]=="scale_e"){
		if (channelname=="e") {
		  sys_up = scalesys_e_up[ibgr];sys_dw = scalesys_e_dw[ibgr];
		}
		else {
		  sys_up = 0.;sys_dw = 0.;}
	      } 
	      else if (name[s]=="scale_mu"){
		if (channelname=="e") {sys_up = 0;sys_dw = 0;}
		else {
		  sys_up = scalesys_mu_up[ibgr];sys_dw = scalesys_mu_dw[ibgr];}
	      }
	      else if (name[s]=="norm_tt") {
		if (bgrname[ibgr]=="tt"){
		  if (channelname=="e") {
		    sys_up = normsys_e_up[ibgr];sys_dw = normsys_e_dw[ibgr];}
		  else {
		    sys_up = normsys_mu_up[ibgr];sys_dw = normsys_mu_dw[ibgr];}
		} 
		else {
		  sys_up = 0;sys_dw = 0;
		}
	      } 
	      else if (name[s]=="norm_W" ){
		if (bgrname[ibgr]=="W"){
		  if (channelname=="e") {
		    sys_up = normsys_e_up[ibgr];sys_dw = normsys_e_dw[ibgr];}
		  else {
		    sys_up = normsys_mu_up[ibgr];sys_dw = normsys_mu_dw[ibgr];}
		} 
		else {
		  sys_up = 0;sys_dw = 0;
		}
	      } 
	      else if (name[s]=="norm_QCD"){
		if (bgrname[ibgr]=="QCD" || bgrname[ibgr]=="QCDe" || bgrname[ibgr]=="QCDmu"){
		  if (channelname=="e") {
		    sys_up = normsys_e_up[ibgr];sys_dw = normsys_e_dw[ibgr];}
		  else {
		    sys_up = normsys_mu_up[ibgr];sys_dw = normsys_mu_dw[ibgr];}
		} 
		else {
		  sys_up = 0;sys_dw = 0;  //bugfix from andrew
		}
	      }
	      else if (name[s]=="norm_QCDe"){
		if (bgrname[ibgr]=="QCD" || bgrname[ibgr]=="QCDe" || bgrname[ibgr]=="QCDmu"){
		  //cout<<"andrew: test0 \n";
		  if (channelname=="e"){//ele channel shift only
		    sys_up = normsys_e_up[ibgr];sys_dw = normsys_e_dw[ibgr];
		  } 
		  else {
		    sys_up = 0;sys_dw = 0;
		  }
		}
		else {
		  sys_up = 0;sys_dw = 0; //bugfix from andrew
		}
	      }
	      else if (name[s]=="norm_QCDmu"){
		if (bgrname[ibgr]=="QCD" || bgrname[ibgr]=="QCDe" || bgrname[ibgr]=="QCDmu"){
		  if (channelname=="mu"){//muo channel shift only
		    //cout<<"andrew: test0 \n";
		    sys_up = normsys_mu_up[ibgr];sys_dw = normsys_mu_dw[ibgr];
		  } 
		  else {
		    sys_up = 0;sys_dw = 0;
		  }
		}
		else {
		  sys_up = 0;sys_dw = 0;
		}
	      }
	      else if (name[s]=="norm_Z"){
		if (bgrname[ibgr]=="Z"){
		  if (channelname=="e") {
		    sys_up = normsys_e_up[ibgr];sys_dw = normsys_e_dw[ibgr];}
		  else {
		    sys_up = normsys_mu_up[ibgr];sys_dw = normsys_mu_dw[ibgr];}
		} 
		else {
		  sys_up = 0;sys_dw = 0;
		}
	      } 
	      else if (name[s]=="norm_single-top"){
		if (bgrname[ibgr]=="single-top"){
		  if (channelname=="e") {
		    sys_up = normsys_e_up[ibgr];sys_dw = normsys_e_dw[ibgr];}
		  else {
		    sys_up = normsys_mu_up[ibgr];sys_dw = normsys_mu_dw[ibgr];}
		} 
		else {
		  sys_up = 0;sys_dw = 0;
		}
	      } 
	      else if (name[s]=="norm_Diboson"){
		if (bgrname[ibgr]=="Diboson"){
		  if (channelname=="e") {
		    sys_up = normsys_e_up[ibgr];sys_dw = normsys_e_dw[ibgr];}
		  else {
		    sys_up = normsys_mu_up[ibgr];sys_dw = normsys_mu_dw[ibgr];}
		} 
		else {
		  sys_up = 0;sys_dw = 0;
		}
	      }
	      //cout<<"andrew: "<<sys_up<<" "<<sys_dw<<endl;
	      if(isample==0){
		this_histname=prefix_histname+bgrIn[ibgr];
		if (ibgr==0){
		  cout << this_histname << std::endl;
		  histBgrtmp2= (TH1D*) infile->Get(this_histname+suffix)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		  if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
		  histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		  histBgrtmp2->Scale(1.+sys_up); 
		  histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
		  histBgrtmp_up= (TH1D*) infile->Get(this_histname+suffix)->Clone("Bgr_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		  histBgrtmp_up->Scale(Bgr_SF[ibgr][ichannel]);
		  histBgrtmp_up->Scale(1.+sys_up);
		  histBgrtmp2= (TH1D*) infile->Get(this_histname+suffix)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		  histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		  histBgrtmp2->Scale(1.+sys_dw); 
		  histSingleBgrSys[s][1][ibgr]=histBgrtmp2;
		  histBgrtmp_dw= (TH1D*) infile->Get(this_histname+suffix)->Clone("Bgr_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		  histBgrtmp_dw->Scale(Bgr_SF[ibgr][ichannel]);
		  histBgrtmp_dw->Scale(1.+sys_dw);
		}
		else{
		  histBgrtmp2=  (TH1D*) infile->Get(this_histname+suffix)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+up);
		  if(doDebug) cout << "SF " << Bgr_SF[ibgr][ichannel] << endl;
		  histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		  histBgrtmp2->Scale(1.+sys_up);
		  histSingleBgrSys[s][0][ibgr]=histBgrtmp2;
		  histBgrtmp_up->Add(histBgrtmp2);
		  histBgrtmp2=  (TH1D*) infile->Get(this_histname+suffix)->Clone(bgrname[ibgr]+"_"+reco[ireco]+"_"+channelname+"_"+name[s]+dw);
		  histBgrtmp2->Scale(Bgr_SF[ibgr][ichannel]);
		  histBgrtmp2->Scale(1.+sys_dw); 
		  histSingleBgrSys[s][1][ibgr]=histBgrtmp2;
		  histBgrtmp_dw->Add(histBgrtmp2);
		}
	      }
	    }
	  }
	  std::cout<< " fill sys array " << std::endl;//--
	  histSignalSys[s][0]=histSignaltmp_up;
	  histSignalSys[s][1]=histSignaltmp_dw;
	  if(isample==0){
	    histBgrSys[s][0]=histBgrtmp_up;
	    histBgrSys[s][1]=histBgrtmp_dw;
	  }
	  std::cout <<  name[s] << " ----> done " <<std::endl;
	}//end for(int s=0; s< NSYS; s++): loop over systematics
	
	if(isample==0 && doData){
	  std::cout << " finally data " ;
	  
	  histData= (TH1D*) infile->Get(prefix_histname+dataname+suffix)->Clone("datahisto_"+reco[ireco]+"_"+channelname);
	  std::cout <<  " ----> done " << std::endl;
	}
	

	//Rebinning histograms + set first bins to 0
	if(!doPartial) {
	  histSignal->Print();
	  //andrew no unbinned histSignal->Write();

	  //for(int abin=0; abin<22; abin++){
	  //  cout<<"error -> 0 " <<abin<<endl;
	    //histSignal->SetBinError(abin, 0);
	  //}
	  
	  for(int b=0; b<Nbinnings;++b) {
	    TH1D* newhist= rebinTH1D(histSignal,binname[b],false);
	    if(doSetFirstBinsToZero) newhist = RemoveFirstBins(newhist, BinsZeroThreshold[ireco], DebugSetFirstBinsToZero);
	    newhist->Write();
	  }
	  
	  

	  
	  if (isample==0){
	    histBgr->Print();
	    //if(!do80GeVbinning) 
	    //andrew no unbinned histBgr->Write();

	    if(doData) {
	      histData->Print();
	      //if(!do80GeVbinning) 
	      //andrew no unbinned histData->Write();
	    }
	    for(int b=0; b<Nbinnings;++b) {
	      TH1D* newhist= rebinTH1D(histBgr,binname[b],false);
	      if(doSetFirstBinsToZero) newhist = RemoveFirstBins(newhist, BinsZeroThreshold[ireco], DebugSetFirstBinsToZero);
	      newhist->Write();
	      if(doData) {
		newhist= rebinTH1D(histData,binname[b],false);
		if(doSetFirstBinsToZero) newhist = RemoveFirstBins(newhist, BinsZeroThreshold[ireco], DebugSetFirstBinsToZero);
		newhist->Write();
	      }
	    }
	  }
	}//end if(!doPartial)
	for(Int_t s=0;s<NSYS;s++) {
	  if(name[s]=="") continue;
	  std::cout << " XXXX summary for " << m_SignalDample[isample] << std::endl;

	  for(int i=0; i<2;++i) {//loop over up/down
	    histSignalSys[s][i]->Print();
	    //andrew no unbinned histSignalSys[s][i]->Write();
	    for(int b=0; b<Nbinnings;++b) {
	      TH1D* newhist= rebinTH1D(histSignalSys[s][i],binname[b],false);
	      if(doSetFirstBinsToZero) newhist = RemoveFirstBins(newhist, BinsZeroThreshold[ireco], DebugSetFirstBinsToZero);
	      newhist->Write();
	    }
	    
	    if(isample==0){
	      histBgrSys[s][i]->Print();
	      //if(!do80GeVbinning) 
	      //andrew no unbinned histBgrSys[s][i]->Write();
	      for(int b=0; b<Nbinnings;++b) {
		TH1D* newhist= rebinTH1D(histBgrSys[s][i],binname[b],false);
		if(doSetFirstBinsToZero) newhist = RemoveFirstBins(newhist, BinsZeroThreshold[ireco], DebugSetFirstBinsToZero);
		newhist->Write();
	      }
	    }
	  }//end loop over up/down
	}//end loop over systematics
	if (isample==0) {
	  for (int ibgr=0; ibgr<NBGR; ibgr++){
	    if(bgrIn[ibgr]=="") continue;
	    if(!doPartial) {
	      histBgrtmp[ibgr]->Print();
	      //andrew no unbinned histBgrtmp[ibgr]->Write();
	      for(int b=0; b<Nbinnings;++b) {
		//cout<<"andrew : "<<histBgrtmp[ibgr]->Integral()<<endl;
		TH1D* newhist=rebinTH1D(histBgrtmp[ibgr],binname[b],false);
		//cout<<"andrew : "<<newhist->Integral()<<endl;
		if(doSetFirstBinsToZero) newhist = RemoveFirstBins(newhist, BinsZeroThreshold[ireco], DebugSetFirstBinsToZero);
		//cout<<"andrew : "<<BinsZeroThreshold[ireco]<<" "<<newhist->Integral()<<endl;
		newhist->Write();
	      }
	    }//end if(!doPartial)
	    
	    //Loop over bgr syst
	    for(Int_t s=0;s<NSYS;s++) {
	      if(name[s]=="") continue;
	      for(int i=0; i<2;++i) {//loop over up/down
		histSingleBgrSys[s][i][ibgr]->Print();
	        //andrew no unbinned histSingleBgrSys[s][i][ibgr]->Write();
		for(int b=0; b<Nbinnings;++b) {
		  TH1D* newhist=rebinTH1D(histSingleBgrSys[s][i][ibgr],binname[b],false);
		  if(doSetFirstBinsToZero) newhist = RemoveFirstBins(newhist, BinsZeroThreshold[ireco], DebugSetFirstBinsToZero);
		  newhist->Write();
		}
	      }//end loop over up/down
	    }//end loop over systs
	    std::cout << " done, next " << std::endl;
	  }//end loop over backgrounds
	}//end if isample==0
      }
    }
  }
 
 cout<<"closing files \n";
 cout<<infile<<endl;

 //outfile->ls();
 //infile->Close();
 //cout<<"infile closed \n";
 outfile->Close();
 cout<<"outfile closed \n";
 if(doSetFirstBinsToZero) {
   std::cout << std::endl << " INFO: the first bins set to 0" << std::endl;
   for(int nm=0; nm<nreco; nm++)
     {
       cout<< " low bin threshold = " << BinsZeroThreshold[nm] << " GeV  "<<reco[nm]<<endl;
     }
     //<< " low bin threshold = " << BinsZeroThreshold[1] << " GeV (resolved), "
     //<< BinsZeroThreshold[0] << " GeV (boosted)." << std::endl;
  }
  std::cout << std::endl << " ==== All done, ending ====" << std::endl;
  return 0;
}



