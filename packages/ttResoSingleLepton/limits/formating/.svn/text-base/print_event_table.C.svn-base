#include<iostream>
#include<iomanip>
#include "assert.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "config_list.h"
#include "TString.h"

using namespace std;
enum{e,mu,comb};

//void print_event_table(string dir=".", string filename="Systematics_all", string srebin="binkkg");
void print_event_table(TString dir=".", TString filename="MttSpectra_Oct06_smoothed.root", string srebin="binkkg");


void print_event_table(TString dir, TString filename, string srebin) {
  //Print number of events
  //ebergeas@cern.ch 2012

  const int specsign=4; //index of signal to look at. 4=Z' 1.5

  //const bool doData=false;
  const bool doData=true;

  //const bool doQCD=false;
  const bool doQCD=true;

  bool doDebug=false;
  //bool doDebug=true; 

  //string prefix = srebin+"_"; // no binkkg_ prefix anymore
  string prefix = ""; 

  const int Nch=2;
  string channelname[Nch]={"_e","_mu"};
  string addout="_add";

  
  const int Nsel=2; 
  string selname[Nsel] = {"massTTbarChi2LPC", "masstT"};
  string selprint[Nsel]={"Resolved", "Boosted"};
  //string selname[Nsel] = {"massTTbarDRMin", "massTTbarChi2LPC", "masstT"};
  //string selprint[Nsel]={"Resolved", "Resolved chi2", "Boosted"};
  //string selname[Nsel] = {"massTTbarChi2LPC"};
  //string selprint[Nsel]={"Resolved chi2"};

  const int Ntop=3;
  enum{tt,st,ttv};
  string topname[Ntop]={"tt_", "single-top_","ttV_"};
  string topprint[Ntop]={"\\ttbar", "Single top","ttV"};
  string tout="top";
  double Evt_top[Ntop][Nsel][Nch];

  
  const int Ndd=3;
  enum{QCDe,QCDmu,W};
  string ddname[Ndd]={"QCDe_", "QCDmu_", "W_"};
  string ddprint[Ndd]={"QCD e", "QCD mu", "$W$+jets"};
  string dout="dd";
  double Evt_dd[Ndd][Nsel][Nch];

  const int Nvv=2;
  enum{Z,Db};
  string vvname[Nvv]={"Z_", "Diboson_"};
  string vvprint[Nvv]={"$Z$+jets", "Di-bosons"};
  string vout="vv";
  double Evt_vv[Nvv][Nsel][Nch];

  string mcout="mc";

  double Evt_sign[Nsel][Nch];

  string bkgname="Bgr_";
  double Evt_totele[Nsel], Evt_totmuo[Nsel], Evt_totadd[Nsel];
  double Evt_signele[Nsel], Evt_signmuo[Nsel], Evt_signadd[Nsel];
  for(int nn=0;nn<Nsel;++nn) {
    for(int cc=0;cc<Nch;++cc) {
      Evt_totele[nn]=0; 
      Evt_totmuo[nn]=0;
      Evt_totadd[nn]=0;
      Evt_signele[nn]=0; 
      Evt_signmuo[nn]=0; 
      Evt_signadd[nn]=0;
    }
  }
  
  string dataname="data_";
  string dataprint="Data";
  double Evt_data[Nsel][Nch];
  
  //const int Nsign=26;
  const int Nsign=16;
     
  string signalname[Nsign];
  double signalxsec[Nsign];
  
  signalname[0] = "tt_";   signalxsec[0]=23.17;
  signalname[1] = "Z750_";   signalxsec[1]=5.60;
  signalname[2] = "Z1000_";  signalxsec[2]=1.61;
  signalname[3] = "Z1250_";  signalxsec[3]=0.57;
  signalname[4] = "Z1500_"; signalxsec[4]=0.21;
  signalname[5] = "Z1750_";  signalxsec[5]=-1;
  signalname[6] = "Z2000_";  signalxsec[6]=0.04;
  signalname[7] = "Z2250_";  signalxsec[7]=-1;
  signalname[8] = "Z2500_";  signalxsec[8]=0.007;
  signalname[9] = "Z3000_";  signalxsec[9]=0.002;
  signalname[10] = "KKg700_";   signalxsec[10]=25.20;
  signalname[11] = "KKg800_";   signalxsec[11]=14.63;
  signalname[12] = "KKg1000_";  signalxsec[12]=5.473;
  signalname[13] = "KKg1300_";  signalxsec[13]=1.523;
  signalname[14] = "KKg1600_";  signalxsec[14]=0.5;
  signalname[15] = "KKg1800_";  signalxsec[15]=0.255;


  /*signalname[0] = "Z500_"; signalxsec[0]=19.6;
  signalname[1] = "Z600_"; signalxsec[1]=10.3;
  signalname[2] = "Z700_"; signalxsec[2]=5.6;
  signalname[3] = "Z800_"; signalxsec[3]=3.2;
  signalname[4] = "Z1000_"; signalxsec[4]=1.2;
  signalname[5] = "Z1300_"; signalxsec[5]=0.30;
  signalname[6] = "Z1600_"; signalxsec[6]=0.086;
  signalname[7] = "Z2000_"; signalxsec[7]=0.019;
  signalname[8] = "Z2500_"; signalxsec[8]=0.0030;
  signalname[9] = "Z3000_"; signalxsec[9]=0.0;
  
  signalname[10] = "KKg700_"; signalxsec[10]=20.8;
  signalname[11] = "KKg800_"; signalxsec[11]=11.6;
  signalname[12] = "KKg900_"; signalxsec[12]=6.8;
  signalname[13] = "KKg1000_"; signalxsec[13]=4.1;
  signalname[14] = "KKg1150_"; signalxsec[14]=2.1;
  signalname[15] = "KKg1300_"; signalxsec[15]=1.1;
  signalname[16] = "KKg1600_"; signalxsec[16]=0.35;
  signalname[17] = "KKg1800_"; signalxsec[17]=0.18;
  signalname[18] = "KKg2000_"; signalxsec[18]=0.095;
  
  //KK graviton samples
  signalname[19] = "KKGrav600_"; signalxsec[19]=2.2;
  signalname[20] = "KKGrav700_"; signalxsec[20]=0.93;
  signalname[21] = "KKGrav800_"; signalxsec[21]=0.42;
  signalname[22] = "KKGrav900_"; signalxsec[22]=0.0;
  signalname[23] = "KKGrav1000_"; signalxsec[23]=0.099;
  signalname[24] = "KKGrav1150_"; signalxsec[24]=0.0;
  signalname[25] = "KKGrav1300_"; signalxsec[25]=0.015;
  */




  string hname;

  TH1D *d_mcall;
  TH1D *d_data;

  string signalfile;
  TFile *x;
  if(!dir.EndsWith("/")) dir+="/";
  if(!filename.EndsWith(".root")) filename+=".root";
  signalfile=dir+filename; 
  //cout << "Signalfile: " << signalfile << endl;
  x=0;
  x = TFile::Open(signalfile.data());
  if (!x) cout << signalfile << " does not exist" << endl;
  //else cout << "File opened!" << endl;

  hname=prefix+topname[0]+selname[0]+channelname[0];
  if(doDebug) cout << hname << endl;

  if (x){
    for(int n=0;n<Nsel;++n) { //loop over selections (azpr,azbt)

      for(int c=0;c<Nch;c++) {  

	//data
	if(doData) {
	  hname=prefix+dataname+selname[n]+channelname[c];
	  cout << hname << endl;
	  d_data=(TH1D*)x->Get(hname.data());
	  Evt_data[n][c]=d_data->Integral(-1,9999);
	}
	else {
	  if(doDebug) cout << "No data stored" << endl;
	  Evt_data[n][c]=0;
	}

	//--------------------
	//MC, top
	for(int i=0; i<Ntop;++i) {
	  hname=prefix+topname[i]+selname[n]+channelname[c];
	  if(doDebug) cout << hname << endl;
	  d_mcall=(TH1D*)x->Get(hname.data());
	  Evt_top[i][n][c]=d_mcall->Integral(-1,9999);
	}
	 
	//MC, vv
	for(int i=0; i<Nvv;++i) {
	  hname=prefix+vvname[i]+selname[n]+channelname[c];
	  if(doDebug) cout << hname << endl;
	  d_mcall=(TH1D*)x->Get(hname.data());
	  Evt_vv[i][n][c]=d_mcall->Integral(-1,9999);  
	}
	
	//DD
	for(int i=0; i<Ndd;++i) {
	  if(doQCD || !TString(ddname[i]).Contains("QCD")) {
	    hname=prefix+ddname[i]+selname[n]+channelname[c];
	    if(doDebug) cout << hname << endl;
	    d_mcall=(TH1D*)x->Get(hname.data());
	    if(!d_mcall) d_mcall=new TH1D("dummy","dummy",10,0,10);
	    Evt_dd[i][n][c]=d_mcall->Integral(-1,9999);  
	  }
	  else Evt_dd[i][n][c]=0;
	}

	//Signal
	double thisevt;
	hname=prefix+signalname[specsign]+selname[n]+channelname[c];
	if(doDebug) cout << hname << endl;
	d_mcall=(TH1D*)x->Get(hname.data());
	thisevt=d_mcall->Integral(-1,9999);
	Evt_sign[n][c]=thisevt*signalxsec[specsign]; 
	//if(doDebug) cout<<n<<" "<<c<<endl;
	//if(doDebug) cout<<d_mcall->Integral(-1,9999)<<" "<<d_mcall->Integral(-1,9999)<<" "<<d_mcall->GetEntries()<<" "<<signalxsec[specsign]<<" "<<Evt_sign[n][c]<<endl;
		
	cout << channelname[c] << " done!" << endl<< endl;
      }//end loop over channels

      

      ///Print the table

      //absolute error
      if(doDebug) cout<<"getting absolute errors \n";
       vector<double> tt_err=getAbsErr(x,  prefix,"tt_", selname[n]);
       vector<double> singleTop_err=getAbsErr(x,  prefix,"single-top_", selname[n]);
       vector<double> ttV_err=getAbsErr(x,  prefix,"ttV_", selname[n]);
       vector<double> QCDe_err=getAbsErr(x,  prefix,"QCDe_", selname[n]);
       vector<double> QCDmu_err=getAbsErr(x,  prefix,"QCDmu_", selname[n]);
       vector<double> Wjet_err=getAbsErr(x,  prefix,"W_", selname[n]);
       vector<double> Z_err=getAbsErr(x,  prefix,"Z_", selname[n]);
       vector<double> Diboson_err=getAbsErr(x,  prefix,"Diboson_", selname[n]);
       vector<double> Bgr_err=getAbsErr(x,  prefix,"Bgr_", selname[n]);
       vector<double> sign_err=getAbsErr(x,  prefix,signalname[specsign], selname[n]);
    
       double err_top[Ntop][3];
       double err_dd[Ndd][3];
       double err_vv[Nvv][3];
       double err_Bgr[3];
       double err_sign[3];

       for(int c=0;c<3;c++)
	 {
	   err_top[tt][c]=tt_err.at(c);
	   err_top[st][c]=singleTop_err.at(c);
	   err_top[ttv][c]=ttV_err.at(c);
	   err_dd[QCDe][c]=QCDe_err.at(c);
	   err_dd[QCDmu][c]=QCDmu_err.at(c);
	   err_dd[W][c]=Wjet_err.at(c);
	   err_vv[Z][c]=Z_err.at(c);
	   err_vv[Db][c]=Diboson_err.at(c);
	   err_Bgr[c]=Bgr_err.at(c);
	   err_sign[c]=sign_err.at(c);
	 }

         

       
       //Norm errors
      double norm_top[Ntop][2];
      double norm_dd[Ndd][3];
      double norm_vv[Nvv][2];
      
      double norm_sign[2];

      //ele
      //ttbar, single top
      norm_top[0][0]=0.107;  //ttbar 
      norm_top[1][0]=0.077;  //single top
      
      //QCD e, QCD mu, W
      norm_dd[0][0]=0.6;    //qcd e
      norm_dd[1][0]=0.;     //qcd mu
      norm_dd[2][0]=0.19;   // w+jet
      if(selname[n]=="massTTbarChi2LPC" || selname[n]=="massTTbarDRMin") {
	norm_dd[2][0]=0.105;
      }
      
      //Z+j, di-bos
      norm_vv[0][0]=0.48; //z+jet
      norm_vv[1][0]=0.344; //diboson
      
      //muo
      norm_top[0][1]=0.107;
      norm_top[1][1]=0.077;
      
      //QCD e, QCD mu, W
      norm_dd[0][1]=0.0;
      norm_dd[1][1]=0.6;
      norm_dd[2][1]=0.18;
      if(selname[n]=="massTTbarChi2LPC" || selname[n]=="massTTbarDRMin") {
	norm_dd[2][1]=0.105;
      }
     
      //Z+j, di-bos
      norm_vv[0][1]=0.48;
      norm_vv[1][1]=0.344;

      norm_sign[0]=0.05;
      norm_sign[1]=0.05;

      double err_ele=0, err_muo=0, err_add=0, toterr_ele=0, toterr_muo=0, toterr_add=0;

      bool fullErrors=true;
      bool pdgRound=true;

      cout << "Summary table, " << selprint[n]<< endl;
      //top
      for(int b=0; b<Ntop;++b) {
	if(fullErrors){
	  err_ele=err_top[b][e];	 
	  err_muo=err_top[b][mu];
	  err_add=err_top[b][comb];
	}
	else{
	  err_ele=Evt_top[b][n][0]*norm_top[b][0];
	  err_muo=Evt_top[b][n][1]*norm_top[b][1];
	  err_add=err_ele+err_muo;
	}

	if(!pdgRound){
	if(Evt_top[b][n][0]>5) cout<<setiosflags(ios::fixed)<<setprecision(0);
	else if(Evt_top[b][n][0]>1) cout<<setiosflags(ios::fixed)<<setprecision(1);
	else cout<<setiosflags(ios::fixed)<<setprecision(2);
	}
	if(pdgRound) SetPrecisionPDG(err_ele);
	cout << topprint[b] 
	     << " & " << Evt_top[b][n][0] << " &  " << err_ele;
	if(pdgRound) SetPrecisionPDG(err_muo);
	cout<< " & " << Evt_top[b][n][1] << " &  " << err_muo;
	if(pdgRound) SetPrecisionPDG(err_add);
	cout<< " & " << Evt_top[b][n][0]+Evt_top[b][n][1] << " &  "
	     << err_add
	     << " \\\\ \\hline" << endl;
	toterr_ele+= err_ele*err_ele;
	toterr_muo+= err_muo*err_muo;
	Evt_totele[n]+=Evt_top[b][n][0];
	Evt_totmuo[n]+=Evt_top[b][n][1];
	Evt_totadd[n]+=Evt_top[b][n][0]+Evt_top[b][n][1];
      }
      //dd
      for(int b=0; b<Ndd;++b) {
	if(fullErrors){
          err_ele=err_dd[b][e];
          err_muo=err_dd[b][mu];
          err_add=err_dd[b][comb];
        }
        else{
	  err_ele=Evt_dd[b][n][0]*norm_dd[b][0];
	  err_muo=Evt_dd[b][n][1]*norm_dd[b][1];
	  err_add=err_ele+err_muo;
	}
	if(!pdgRound){
	  if(Evt_dd[b][n][0]>5) cout<<setiosflags(ios::fixed)<<setprecision(0);
	  else if(Evt_dd[b][n][0]>1) cout<<setiosflags(ios::fixed)<<setprecision(1);
	  else cout<<setiosflags(ios::fixed)<<setprecision(2);
	}
	if(pdgRound) SetPrecisionPDG(err_ele);
	cout << ddprint[b] 
	     << " & " << Evt_dd[b][n][0] << " &  " << err_ele;
	if(pdgRound) SetPrecisionPDG(err_muo);
	cout<< " & " << Evt_dd[b][n][1] << " &  " << err_muo;
	if(pdgRound) SetPrecisionPDG(err_add);
	cout<< " & " << Evt_dd[b][n][0]+Evt_dd[b][n][1] << " &  "
	    << err_add
	    << " \\\\ \\hline" << endl;
	toterr_ele+= err_ele*err_ele;
	toterr_muo+= err_muo*err_muo;
	Evt_totele[n]+=Evt_dd[b][n][0];
	Evt_totmuo[n]+=Evt_dd[b][n][1];
	Evt_totadd[n]+=Evt_dd[b][n][0]+Evt_dd[b][n][1];

      }
      //vv
      for(int b=0; b<Nvv;++b) {
	if(fullErrors){
          err_ele=err_vv[b][e];
          err_muo=err_vv[b][mu];
          err_add=err_vv[b][comb];
        }
        else{
	  err_ele=Evt_vv[b][n][0]*norm_vv[b][0];
	  err_muo=Evt_vv[b][n][1]*norm_vv[b][1];
	  err_add=err_ele+err_muo;
	}
	if(!pdgRound){
	  if(Evt_vv[b][n][0]>5) cout<<setiosflags(ios::fixed)<<setprecision(0);
	  else if(Evt_vv[b][n][0]>1) cout<<setiosflags(ios::fixed)<<setprecision(1);
	  else cout<<setiosflags(ios::fixed)<<setprecision(2);
	}
	if(pdgRound) SetPrecisionPDG(err_ele);
	cout << vvprint[b] 
	     << " & " << Evt_vv[b][n][0] << " &  " << err_ele;
	if(pdgRound) SetPrecisionPDG(err_muo);
	cout<< " & " << Evt_vv[b][n][1] << " &  " << err_muo;
	if(pdgRound) SetPrecisionPDG(err_add);
	cout<< " & " << Evt_vv[b][n][0]+Evt_vv[b][n][1] << " &  "
	     << err_add
	     << " \\\\ \\hline" << endl;
	toterr_ele+= err_ele*err_ele;
	toterr_muo+= err_muo*err_muo;
	Evt_totele[n]+=Evt_vv[b][n][0];
	Evt_totmuo[n]+=Evt_vv[b][n][1];
	Evt_totadd[n]+=Evt_vv[b][n][0]+Evt_vv[b][n][1];

      }
      //total
      if(fullErrors){
	toterr_ele=err_Bgr[e];
	toterr_muo=err_Bgr[mu];
	toterr_add=err_Bgr[comb];
      }
      else{
	toterr_ele=sqrt(toterr_ele);
	toterr_muo=sqrt(toterr_muo);
	toterr_add=toterr_ele+toterr_muo;
      }
      
      cout<<setiosflags(ios::fixed)<<setprecision(0);
      

      cout << "Total" 
	   << " & " << Evt_totele[n] << " &  " << toterr_ele 
	   << " & " << Evt_totmuo[n] << " &  " << toterr_muo
	   << " & " << Evt_totadd[n] << " &  " << toterr_add
	   << " \\\\ \\hline" << endl;
      cout<<setiosflags(ios::fixed)<<setprecision(0);
      cout << "Data  & " << Evt_data[n][0] << " & " << Evt_data[n][1] 
	   << "  & " << Evt_data[n][0]+Evt_data[n][1] << " \\\\ \\hline" << endl;
      
      
      if(fullErrors){
        err_ele=err_sign[e];
        err_muo=err_sign[mu];
        err_add=err_sign[comb];
      }
      else{
        err_ele=norm_sign[0]*Evt_sign[n][0];
        err_muo= norm_sign[1]*Evt_sign[n][1];
        err_add=norm_sign[0]*Evt_sign[n][0]+norm_sign[1]*Evt_sign[n][1];
      }


      if(!pdgRound)cout<<setiosflags(ios::fixed)<<setprecision(2);

      if(pdgRound) SetPrecisionPDG(err_ele);
      cout << "Signal " << signalname[specsign] << "&" 
	   << Evt_sign[n][0] << " & " << err_ele<< " & " ;
      if(pdgRound) SetPrecisionPDG(err_muo);
      cout<< Evt_sign[n][1] << " & " << err_muo << " & " ;
      if(pdgRound) SetPrecisionPDG(err_add);
      cout<< Evt_sign[n][0]+Evt_sign[n][1] << " & " << err_add 
	   << " \\\\ \\hline" << endl;
      cout << "Cross section " << signalxsec[specsign] << endl;
      cout << endl;

      cout << "Selection " << selname[n] << " done" << endl;
    }
  }

  if(srebin!="") cout << "Histograms rebinned by " << srebin << endl;
  cout << "All done!" << endl;

}

TH1F* TH1DtoTH1F(TH1D* histD) {
  int allBins = histD->GetNbinsX();
  int nBins = allBins;
  assert (nBins <= 10000);
  double binLimits[10000];
  for (int bin=0; bin <= nBins+2; bin++) { //bugfix
    binLimits[bin] = histD->GetBinLowEdge(bin+1); }
    //binLimits[bin] = histD->GetBinLowEdge(bin); }
    //binLimits[nBins] = histD->GetBinLowEdge(allBins) + histD->GetBinWidth(allBins);

  string hhname=histD->GetName(); hhname+="_f";
  string hhtitle=histD->GetTitle(); hhtitle+="_f";
  //TH1F *histF = new TH1F(histD->GetName(),histD->GetTitle(),nBins,binLimits);
  TH1F *histF = new TH1F(hhname,hhtitle,nBins,binLimits);
  if (histD->GetSumw2N()) histF->Sumw2();

  for (int bin=0; bin <= nBins+1; bin++) { //including the over/underflow
    histF->SetBinContent(bin,histD->GetBinContent(bin));
    histF->SetBinError(bin,histD->GetBinError(bin));
  }
  
  //cross check
  /*
  cout << "histD    histF" << endl;
  for (int bin=0; bin <= nBins+1; bin++) { //including the over/underflow
    histF->SetBinContent(bin,histD->GetBinContent(bin));
    histF->SetBinError(bin,histD->GetBinError(bin));
  }
  */
  return histF;
}


void FillByBinFluctuation(const TH1* h1, TH1* h2, TRandom3* randomEngine) {
  //Get a histogram representing expected events in each bin, and fluctuate the content
  //of its bins according to Poisson, to simulate observed data.
  //G. Choudalakis 2011

  unsigned int nBins1 = h1->GetNbinsX();
  unsigned int nBins2 = h2->GetNbinsX();
  assert(nBins1==nBins2);
  for (unsigned int bin=0; bin <= nBins1+1; bin++) {
    double expected = h1->GetBinContent(bin);
    double observed;
    observed = randomEngine->PoissonD(expected);
    h2->SetBinContent(bin,observed);
    h2->SetBinError(bin,sqrt(observed));
  }
  return;
}


//number of digits to display after the decimal point
int getPrecision(double num)
{
  double n=fabs(num); //just in case...
  if(n==0) return 1;

  
  int digits=2;
  if(n>=.100 && n<.355) return 2;
  if(n>=.355 && n<.950) return 1;
  if(n>=.950 && n<1.000) return 1;
  
  if(n>=1.00)return getPrecision(n/10)-1;
  if(n<.100) return getPrecision(n/10)+1;

  cout<<"something has gone wrong in getPrecision()! \n";
  return 2;
}

int SetPrecisionPDG(double num)
{
  int prec= getPrecision(num);
  if(prec<0) prec=0;
  cout<<setiosflags(ios::fixed)<<setprecision(prec);
    //cout<<setiosflags(ios::fixed)<<setprecision(5);
  return prec;

}




void print(double num)
{
  cout<<setiosflags(ios::fixed)<<setprecision(getPrecision(num));
  cout<<num<<endl;
}


//  binkkg_tt_massTTbarChi2LPC_e


vector<double> getAbsErr(string filename="07_02_2013_glgw_Systematics.root", string prefix="binkkg_", string mcName="Bgr_", string selName="masstT_")
{
  
  TFile * file = TFile::Open(filename);
  if(!file)cout<<"failed to open file! \n";
  //cout<<"test\n";
  return getAbsErr(file);
}  
vector<double> getAbsErr(TFile *file, string prefix="binkkg_", string mcName="Bgr_", string selName="masstT_")
//vector<double> getAbsErr(TFile *file, string prefix, string mcName, string selName)
{
  bool debug=false;
  bool AllSyst=true;
  
  //if(mcName="Bgr_") debug=true;

  //const int Nsysts=31;
  //string syst[Nsysts];
  vector<string> syst;

  buildSystList(&syst);

//  syst.clear();
//  syst.push_back("EleSF");
//  syst.push_back("MuSF");
// 
//  syst.push_back("MCGen");
//  syst.push_back("PartonShower");
//  syst.push_back("EWS");
//  syst.push_back("IFSR");
//  syst.push_back("topmass");
//  syst.push_back("PDF");
//
//  syst.push_back("luminosity");
//  syst.push_back("norm_tt");
//  syst.push_back("mcStat");
//
//  syst.push_back("BoostedJES0");
//  syst.push_back("BoostedJES13");
//  syst.push_back("BoostedJES14");
//  syst.push_back("BoostedJES15");
//  syst.push_back("BoostedJES16");
//  syst.push_back("BoostedJMS");
//  syst.push_back("BoostedJER");
//  syst.push_back("BoostedJMR");
//
//
//  /*syst.push_back("Btag0");
//  syst.push_back("Btag1");
//  syst.push_back("Btag2");
//  syst.push_back("Btag3");
//  syst.push_back("Btag4");
//  syst.push_back("Btag5");
//  syst.push_back("Btag6");*/
//  syst.push_back("Btag7");
//  syst.push_back("Btag8");
//  syst.push_back("Btag9");
//  syst.push_back("Btag10");
//  //syst.push_back("smallBtag");
//  /*  syst.push_back("BtagC0");
//  syst.push_back("BtagC1");
//  syst.push_back("BtagC2");
//  syst.push_back("BtagC3");
//  syst.push_back("BtagC4");
//  syst.push_back("BtagC5");
//  syst.push_back("BtagC6");
//  syst.push_back("BtagL0");
//  syst.push_back("BtagL1");
//  syst.push_back("BtagL2");
//  syst.push_back("BtagL3");
//  syst.push_back("BtagL4");
//  syst.push_back("BtagL5");
//  syst.push_back("BtagL6");
//  syst.push_back("BtagL7");
//  syst.push_back("BtagL8");
//  syst.push_back("BtagL9");
//  syst.push_back("BtagL10");
//  syst.push_back("BtagL11");
//  syst.push_back("BtagL12");*/
//  //syst.push_back("Btag");
//  syst.push_back("BtagC");
//  syst.push_back("BtagL");
// 
//  
//  //syst.push_back("JES0");
//  //syst.push_back("JES1");
//  //syst.push_back("JES2");
//  syst.push_back("JES3");
//  //syst.push_back("JES4");
//  //syst.push_back("JES5");
//  //syst.push_back("JES6");
//  syst.push_back("JES7");
//  //syst.push_back("JES8");
//  //syst.push_back("JES9");
//  //syst.push_back("JES10");
//  //syst.push_back("JES11");
//  syst.push_back("JES12");
//  //syst.push_back("JES13");
//  //syst.push_back("JES14");
//  //syst.push_back("JES15");
//  //syst.push_back("JES16");
//  //syst.push_back("JES17");
//  //syst.push_back("JES18");
//  //syst.push_back("JES19");
//  syst.push_back("JES20");
//  //syst.push_back("JES21");
//  //syst.push_back("JES22");
//
//  syst.push_back("JESsmall");
//  syst.push_back("JetEnerRes");
//
//
//  syst.push_back("iqopt3");
//  syst.push_back("ptjmin10");
//  syst.push_back("Whfsf");
//
//  syst.push_back("norm_QCDe");
//  syst.push_back("norm_QCDmu");
  




  int Nsysts=syst.size();
 
  //cout<<"Nsysts="<<syst.size()<<endl;

  
  string hname_el=prefix+mcName+selName+"_e";
  string hname_mu=prefix+mcName+selName+"_mu";
  
  if(debug)cout<<hname_el<<endl;
  
  
  TH1D *d_mcall_el=(TH1D*)file->Get(hname_el.data());//->Clone();
  TH1D *d_mcall_mu=(TH1D*)file->Get(hname_mu.data());//->Clone();
  if(!d_mcall_el && !d_mcall_mu) {
    vector<double> results(3,0);
    return results;
  }


  unsigned int NBINS=d_mcall_el->GetNbinsX();
  if(debug) cout<<"NBINS:" <<NBINS<<endl;

  vector<double> diffsum2_mx_el, centralvalue_el; //, relerr_mx;
  vector<double> diffsum2_mx_mu, centralvalue_mu; //, relerr_mx;
  //vector< vector<double> > totalerr_mx_el;
  //vector< vector<double> > totalerr_mx_mu;

  double totalerr_mx_el[100][100], totalerr_mx_mu[100][100];
  
  //diffsum2_mx.assign(NBINS+2,0);
  //centralvalue.assign(NBINS+2,0);
  diffsum2_mx_el.resize(NBINS+2);
  diffsum2_mx_mu.resize(NBINS+2);
  centralvalue_el.resize(NBINS+2);
  centralvalue_mu.resize(NBINS+2);
  
  //totalerr_mx_el.resize(Nsysts);
  //totalerr_mx_mu.resize(Nsysts);
  
  
  for(int s=0; s<Nsysts; ++s) {
    //totalerr_mx_el.at(s).resize(NBINS+2);
    //totalerr_mx_mu.at(s).resize(NBINS+2);
    for(int b=0;b<NBINS+2;++b) {
 
      totalerr_mx_el[s][b]=0;
      totalerr_mx_mu[s][b]=0;
      
    }    
  }

  for(int b=0;b<NBINS+2;++b) {
    diffsum_mx_el.at(b)=0;
    diffsum_mx_mu.at(b)=0;
    diffsum2_mx_el.at(b)=0;
    diffsum2_mx_mu.at(b)=0;
    centralvalue_el.at(b)=0;
    centralvalue_mu.at(b)=0;
  }

  double bin0, binu, bind, diff2, du, dd, maxud; //, sdiff;
    //double nom, shift, diff, reldiff;

  double events_el, eventsUp_el, eventsDw_el, events_mu, eventsUp_mu, eventsDw_mu,totErr2_el, totErr2_mu, totErr2_comb;
  totErr2_el=0;totErr2_mu=0;totErr2_comb=0;

  for(int s=0; s<Nsysts; ++s) {
    if(syst[s]=="NONE") continue;
    if(debug) cout<<"systematic: "<<syst[s]<<endl; 
   
    if(syst[s]=="mcStat"){
      double statErr_mu, statErr_el;
      events_el=d_mcall_el->IntegralAndError(-1,9999, statErr_el);
      events_mu=d_mcall_mu->IntegralAndError(-1,9999, statErr_mu);

      if(debug) cout<<"  "<< statErr_el<<" "<<statErr_mu<<" "<<totErr2_el<<" "<<totErr2_mu<<endl;
      totErr2_el+=statErr_el*statErr_el;
      totErr2_mu+=statErr_mu*statErr_mu;
      totErr2_comb+=statErr_el*statErr_el+statErr_mu*statErr_mu;
      continue;
    }

    string hname_el_up=prefix+mcName+selName+"_e_"+syst[s]+"_up";
    string hname_el_dw=prefix+mcName+selName+"_e_"+syst[s]+"_dw";
    string hname_mu_up=prefix+mcName+selName+"_mu_"+syst[s]+"_up";
    string hname_mu_dw=prefix+mcName+selName+"_mu_"+syst[s]+"_dw";
    if(debug) cout<< hname_el_up<<endl;
    TH1D *d_el_up=(TH1D*)file->Get(hname_el_up.data())->Clone();
    TH1D *d_el_dw=(TH1D*)file->Get(hname_el_dw.data())->Clone();
    TH1D *d_mu_up=(TH1D*)file->Get(hname_mu_up.data())->Clone();
    TH1D *d_mu_dw=(TH1D*)file->Get(hname_mu_dw.data())->Clone();


    events_el=d_mcall_el->Integral(-1,9999);
    eventsUp_el=d_el_up->Integral(-1,9999);
    eventsDw_el=d_el_dw->Integral(-1,9999);
    
    events_mu=d_mcall_mu->Integral(-1,9999);
    eventsUp_mu=d_mu_up->Integral(-1,9999);
    eventsDw_mu=d_mu_dw->Integral(-1,9999);
    
    du=fabs(events_el-eventsUp_el);
    dd=fabs(events_el-eventsDw_el);
    

    maxud=TMath::Max(du,dd);
    totErr2_el+=maxud*maxud;

    if(debug) cout<<" "<<du<<" "<<dd<<" "<<maxud<<" "<<totErr2_el<<endl;
    
    du=fabs(events_mu-eventsUp_mu);
    dd=fabs(events_mu-eventsDw_mu);

    maxud=TMath::Max(du,dd);
    totErr2_mu+=maxud*maxud;
    if(debug) cout<<" "<<du<<" "<<dd<<" "<<maxud<<" "<<totErr2_mu<<endl;
    du=fabs(events_el+events_mu-eventsUp_el-eventsUp_mu);
    dd=fabs(events_el+events_mu-eventsDw_el-eventsDw_mu);

    maxud=TMath::Max(du,dd);
    totErr2_comb+=maxud*maxud;
    

    for(int b=0;b<NBINS+2;++b) {
      bin0=d_mcall_el->GetBinContent(b);
      centralvalue_el[b]=bin0;
      binu=d_el_up->GetBinContent(b);
      bind=d_el_dw->GetBinContent(b);
      du=fabs(bin0-binu);
      dd=fabs(bin0-bind);

      maxud=TMath::Max(du,dd);
      diffsum2_mx_el[b]+=maxud*maxud;
      totalerr_mx_el[s][b]=maxud;
      

      //cout<<b<<" "<<bin0<<" "<<binu<<" "<<bind<<endl;
    }// end loop over bins
    
    for(int b=0;b<NBINS+2;++b) {
      bin0=d_mcall_mu->GetBinContent(b);
      centralvalue_mu[b]=bin0;
      binu=d_mu_up->GetBinContent(b);
      bind=d_mu_dw->GetBinContent(b);
      du=fabs(bin0-binu);
      dd=fabs(bin0-bind);

      maxud=TMath::Max(du,dd);
      diffsum2_mx_mu[b]+=maxud*maxud;
      totalerr_mx_mu[s][b]=maxud;

 
    }// end loop over bins

   
  }//end loop over syst

  for(int b=0;b<NBINS+2;++b) {
    diffsum2_mx_el[b]=sqrt(diffsum2_mx_el[b]);
    diffsum2_mx_mu[b]=sqrt(diffsum2_mx_mu[b]);
  }

  double err2_el[50];
  double err2_mu[50];
  double err2_comb[50];
  double err_el_tot=0;
  double err_mu_tot=0;
  double err_comb_tot=0;

  for(int b=0;b<NBINS+2;++b) {
    err2_el[b]=0;
    err2_mu[b]=0;
    err2_comb[b]=0;
    for(int s=0; s<Nsysts; ++s) {
      if(syst[s]=="NONE" || syst[s]=="mcStat") continue;
      //cout<<b<<" "<<s<<" "<< totalerr_mx_mu[s][b]<<" "<<totalerr_mx_el[s][b]<<endl;
      
      double err_comb=totalerr_mx_el[s][b]+totalerr_mx_mu[s][b];
      err2_el[b]+=totalerr_mx_el[s][b]*totalerr_mx_el[s][b];
      err2_mu[b]+=totalerr_mx_mu[s][b]*totalerr_mx_mu[s][b];
      err2_comb[b]+=err_comb*err_comb;
      //cout<<b<<" "<<s<<" "<<err_comb*err_comb<<endl;

	}
    //cout<<b<<" "<<sqrt(err2_comb[b])<<endl;

    err_el_tot+=sqrt(err2_el[b]);
    err_mu_tot+=sqrt(err2_mu[b]);
    err_comb_tot+=sqrt(err2_comb[b]);
    
  }

  

  

  vector<double> resultsBin, results;
  resultsBin.push_back(err_el_tot);
  resultsBin.push_back(err_mu_tot);
  resultsBin.push_back(err_comb_tot);
  results.push_back(sqrt(totErr2_el));
  results.push_back(sqrt(totErr2_mu));
  results.push_back(sqrt(totErr2_comb));

  
  if(debug) cout<<resultsBin.at(0)<<" "<<resultsBin.at(1)<<" "<<resultsBin.at(2)<<endl;
  if(debug) cout<<results.at(0)<<" "<<results.at(1)<<" "<<results.at(2)<<endl;
  return results;

}
