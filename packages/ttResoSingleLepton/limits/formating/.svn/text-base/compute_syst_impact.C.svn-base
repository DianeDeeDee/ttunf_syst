#include<vector>
#include<utility>
#include<iostream> //for cout etc
#include<iomanip>
#include<TCanvas.h>
#include<TF1.h>
#include<TF2.h>
#include<TFile.h>
#include<TGraph.h>
#include<TGraphErrors.h>
#include<TH1.h>
#include<TH2.h>
#include<TLatex.h>
#include "TLegend.h"
#include "TMath.h"
#include<TMatrixD.h>
#include<TMinuit.h>
#include<TROOT.h>
#include<TStyle.h>
#include "functions.h"
#include "config_list.h"

using namespace std;

void compute_syst_impact(string file , string syst="JES", string sshort="", string reco="masstT", bool doSubBkg=false, bool doDebug=false);
void run_all(string file = "glgw_Systematics.root", string reco="masstT", bool doDebug=false);

void  compute_syst_impact()
{

  string file="/afs/cern.ch/user/j/jzhong/WP/public/CommonD3PD/20141006/MttSpectra_Oct06_raw.root";
  
  cout<<"\n\n resolved, \n ";
  run_all(file.c_str(),"massTTbarChi2LPC",false);

  cout<<"\n \n Boosted  \n";
  run_all(file.c_str(),"masstT",false);

}

void run_all(string file, string reco, bool doDebug) {
  //const bool doDebug=false;

  //doDebug=true;
  //run_all("Systematics_all.root", "massTTbarDRMin");
  //run_all("Systematics_all.root", "masstT");
  //run_all("Systematics_all.root", "massTTbarChi2LPC");

  if(doDebug) cout << "Selection " << reco << endl;

  vector<string> syst;
  buildSystList(&syst);
  syst.push_back("mcStat");


//  syst.push_back("mcStat");
//  syst.push_back("luminosity");
//  
//
//  syst.push_back("EleSF");
//  syst.push_back("MuSF");
//
//  //syst.push_back("JES_ALL");
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
//  ////syst.push_back("JES14");
//  //syst.push_back("JES15");
//  //syst.push_back("JES16");
//  //syst.push_back("JES17");
//  //syst.push_back("JES18");
//  ////syst.push_back("JES19");
//  syst.push_back("JES20");
//  //syst.push_back("JES21");
//  //syst.push_back("JES22");
//  
//  syst.push_back("JESsmall");
//  syst.push_back("JetEnerRes");
//  syst.push_back("JVFCut");
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
//  /*syst.push_back("BtagC0");
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
//  syst.push_back("BtagL12");  */
//  //syst.push_back("Btag");
//  syst.push_back("BtagC");
//  syst.push_back("BtagL");
//
//  syst.push_back("norm_tt");
//  syst.push_back("MCGen");
//  syst.push_back("PartonShower");
//  syst.push_back("topmass");
//  syst.push_back("IFSR");
//  syst.push_back("EWS");
//
//  //syst.push_back("iqopt3");
//  //syst.push_back("ptjmin10");
//  syst.push_back("Whfsf");
//
//  syst.push_back("norm_QCDe");
//  syst.push_back("norm_QCDmu");
//
//  syst.push_back("PDF");

  const int Nsysts=syst.size();

 //  cout<<"Systematic effect & "<<"tot.bgr & "<<"\ttbar\ & "<<"sing.top & "<<"$W$+jets &"<<" multi-jet &"<<" $Z$+jets &"<<" Di-bosons &"<<" \Zprime\ 1.5 \TeV\ \\hline"<<endl;
  for(int i=0;i<Nsysts;++i) {
    if(syst[i]=="NONE") continue;
    compute_syst_impact(file, syst[i], syst[i], reco, false, doDebug);
  }
  
}

void compute_syst_impact(string file, string syst, string sshort, string reco, bool doSubBkg, bool doDebug) {

  bool doPrintRelLim=false;
  //if(sshort!="") doPrintRelLim=true;
  string rellim_path="./";

  const bool doData=false; //false

  const int spec_sign=4; //index of signal in question,

  gStyle->SetPalette(1,0);
  gStyle->SetOptStat(0);
  //gROOT->SetStyle("ATLAS");
  //gROOT->ForceStyle();

  if(doDebug) cout << "Selection " << reco << endl;

  //File source
  TFile *_unif = new TFile(file.data(),"read");

  const int NBGR=8;
  string bgrname[NBGR];
  string binname="";
  bgrname[0]=binname+"tt";
  bgrname[1]=binname+"W";
  bgrname[2]=binname+"Z";   
  bgrname[3]=binname+"single-top";   
  bgrname[4]=binname+"Diboson";  
  bgrname[5]=binname+"QCDe";
  bgrname[6]="QCDmu";
  bgrname[7]="ttV";

  const int NSIGN=16, NZpMASSES=9;
  string signame[NSIGN];
  signame[0] = binname+"Z500";
  signame[1] = binname+"Z750";
  signame[2] = binname+"Z1000";
  signame[3] = binname+"Z1250";
  signame[4] = binname+"Z1500";
  signame[5] = binname+"Z1750";
  signame[6] = binname+"Z2000";
  signame[7] = binname+"Z2250";
  signame[8] = binname+"Z2500";  
  signame[9] = binname+"KKg700";
  signame[10] = binname+"KKg800";
  signame[11] = binname+"KKg1000";
  signame[12] = binname+"KKg1300";
  signame[13] = binname+"KKg1600";
  signame[14] = binname+"KKg1800";
  signame[15] = binname+"tt";
  //signame[4] = binname+"HH3000";  // hack

  string dataname="datahisto";

  string name;
  name=syst;
  string prefix="";

  string systType="normal";

  const int Nch=2;
  string channelname[Nch]={"e", "mu"};
  string hname;
  //const int Nreco=2;
  //string reco[Nreco]={"massTTbarDRMin", "masstT"};
  //string reco="masstT";
  //string reco="massTTbarDRMin";
  string canname, cantitle;
  string up="_up", dw="_dw";

  double nEvt[Nch], nEvt_up[Nch], nEvt_dw[Nch];
  double sig_nEvt[Nch], sig_nEvt_up[Nch], sig_nEvt_dw[Nch];
  double nBgrEvt[Nch][NBGR], nBgrEvt_up[Nch][NBGR], nBgrEvt_dw[Nch][NBGR];
  double tot_nEvt=0, tot_nEvt_up=0, tot_nEvt_dw=0;
  double tot_nBgrEvt[NBGR], tot_nBgrEvt_up[NBGR], tot_nBgrEvt_dw[NBGR];
  double totsig_nEvt=0, totsig_nEvt_up=0, totsig_nEvt_dw=0;
  double n, n_up, n_dw;

  for(int c=0;c<Nch;++c) {
    nEvt[c]=0;
    nEvt_up[c]=0;
    nEvt_dw[c]=0;
  }
  for(int b=0;b<NBGR;++b) {
    tot_nBgrEvt[b]=0;
    tot_nBgrEvt_up[b]=0;
    tot_nBgrEvt_dw[b]=0;
  }

  for(int c=0;c<Nch;++c) {

    //Data comparison
    if(doData) {    
      TH1F *h_data;
      hname=dataname+"_"+reco+"_"+channelname[c];
      if(doDebug) cout << hname << endl;
      h_data=(TH1F*)_unif->Get(hname.data());
      
      //over-and underflow check
      int nbins=h_data->GetNbinsX();
      if(doDebug) cout << "Data, nbins " << nbins << endl;
    }

    //Backgrounds
    for(int b=0; b<NBGR;++b) {
      TH1F *h_bgr;
      hname=prefix+bgrname[b]+"_"+reco+"_"+channelname[c];
      if(doDebug) cout << hname << endl;
      h_bgr=(TH1F*)_unif->Get(hname.data());
     
      //over-and underflow check
      n=h_bgr->Integral(-1,9999);
      nBgrEvt[c][b]=n;
      nEvt[c]+=n;
      if(doDebug) cout << "Background " << bgrname[b] << ", nEvt " << n << endl;

      if(syst=="") { cout<<"syst=empty"<<endl; }
      else if (syst=="mcStat") {
	double statErr;
	h_bgr->IntegralAndError(-1,9999, statErr);
	n_up=nBgrEvt[c][b]+statErr;
	nBgrEvt_up[c][b]=n_up;
	nEvt_up[c]+=n_up;
	n_dw=nBgrEvt[c][b]-statErr;
	nBgrEvt_dw[c][b]=n_dw;
	nEvt_dw[c]+=n_dw;
      }
      else {
	if(doDebug) cout << syst << " " << bgrname[b] << endl;
	TH1F *h_bgr_up, *h_bgr_dw;
	hname=prefix+bgrname[b]+"_"+reco+"_"+channelname[c]+"_"+name+up;
	/*if(syst=="PDF" && bgrname[b]!="tt") {
	  hname=prefix+bgrname[b]+"_"+reco+"_"+channelname[c];
	}
	else if(syst.Contains("QCD") && bgrname[b]!="QCD") {
	  hname=prefix+bgrname[b]+"_"+reco+"_"+channelname[c];
	  }*/

	if(doDebug) cout << hname << endl;
	h_bgr_up=(TH1F*)_unif->Get(hname.data());
	
	hname=prefix+bgrname[b]+"_"+reco+"_"+channelname[c]+"_"+name+dw;
	/*if(syst=="PDF" && bgrname[b]!="tt") {
	  hname=prefix+bgrname[b]+"_"+reco+"_"+channelname[c];
	}
	else if(syst.Contains("QCD") && bgrname[b]!="QCD") {
	  hname=prefix+bgrname[b]+"_"+reco+"_"+channelname[c];
	  }*/

	if(doDebug) cout << hname << endl;
	h_bgr_dw=(TH1F*)_unif->Get(hname.data());
	
	//over-and underflow check
	n_up=h_bgr_up->Integral(-1,9999);
	nBgrEvt_up[c][b]=n_up;
	nEvt_up[c]+=n_up;
	n_dw=h_bgr_dw->Integral(-1,9999);
	nBgrEvt_dw[c][b]=n_dw;
	nEvt_dw[c]+=n_dw;
	if(doDebug) {
	  cout << "Syst " << syst<< " up" << ", nEvt " << n_up << endl;
	  cout << "Syst " << syst<< " dw" << ", nEvt " << n_dw << endl;
	}	
      }
      tot_nBgrEvt[b]+=nBgrEvt[c][b];
      tot_nBgrEvt_up[b]+=nBgrEvt_up[c][b];
      tot_nBgrEvt_dw[b]+=nBgrEvt_dw[c][b];
    }
    tot_nEvt+=nEvt[c]; 
    tot_nEvt_up+=nEvt_up[c];
    tot_nEvt_dw+=nEvt_dw[c]; 


    //Signal
    TH1F *h_sig;
    hname=prefix+signame[spec_sign]+"_"+reco+"_"+channelname[c];
    if(doDebug) cout << hname << endl;
    h_sig=(TH1F*)_unif->Get(hname.data());
    
    //over-and underflow check
    n=h_sig->Integral(-1,9999);
    sig_nEvt[c]=n;
    totsig_nEvt+=n;
    if(doDebug) cout << "Signal " << signame[spec_sign] << ", nEvt " << n << endl;

    if(syst=="") { cout<<"syst=empty"<<endl; }
    else if (syst=="mcStat") {
      double statErr;
      h_sig->IntegralAndError(-1,9999, statErr);
      n_up=sig_nEvt[c]+statErr;
      sig_nEvt_up[c]=n_up;
      totsig_nEvt_up+=n_up;
      n_dw=sig_nEvt[c]-statErr;
      sig_nEvt_dw[c]=n_dw;
      totsig_nEvt_dw+=n_dw;
    }
    else {
      if(doDebug) cout << syst << " " << signame[spec_sign] << endl;
      TH1F *h_sig_up, *h_sig_dw;
      hname=prefix+signame[spec_sign]+"_"+reco+"_"+channelname[c]+"_"+name+up;
      
      if(doDebug) cout << hname << endl;
      h_sig_up=(TH1F*)_unif->Get(hname.data());
      
      hname=prefix+signame[spec_sign]+"_"+reco+"_"+channelname[c]+"_"+name+dw;
      
      if(doDebug) cout << hname << endl;
      h_sig_dw=(TH1F*)_unif->Get(hname.data());
      
      //Event count
      n_up=h_sig_up->Integral(-1,9999);
      sig_nEvt_up[c]=n_up;
      totsig_nEvt_up+=n_up;
      n_dw=h_sig_dw->Integral(-1,9999);
      sig_nEvt_dw[c]=n_dw;
      totsig_nEvt_dw+=n_dw;
      if(doDebug) {
	cout << "Syst " << syst<< " up" << ", nEvt " << n_up << endl;
	cout << "Syst " << syst<< " dw" << ", nEvt " << n_dw << endl;
      }	
    }
  }
  
  //Print summary

  double tot_impact=(fabs(1-tot_nEvt_up/tot_nEvt) + fabs(1-tot_nEvt_dw/tot_nEvt))/2.;
  double sig_impact=(fabs(1-totsig_nEvt_up/totsig_nEvt) + fabs(1-totsig_nEvt_dw/totsig_nEvt))/2.;

  double bgr_impact[NBGR];
  for(int ib=0;ib<NBGR;++ib) {
    bgr_impact[ib]=(fabs(1-(nBgrEvt_up[0][ib]+nBgrEvt_up[1][ib])/tot_nBgrEvt[ib]) + fabs(1-(nBgrEvt_dw[0][ib]+nBgrEvt_dw[1][ib])/tot_nBgrEvt[ib]))/2.;
  }
  double tt_impact=bgr_impact[0];
  double stop_impact=bgr_impact[3];
  double W_impact=bgr_impact[1];
  double  QCD_impact=0;
  if(NBGR>5)  QCD_impact=(fabs(1-(nBgrEvt_up[0][6]+nBgrEvt_up[1][6]+nBgrEvt_up[0][5]+nBgrEvt_up[1][5])/(tot_nBgrEvt[6]+tot_nBgrEvt[5])) + fabs(1-(nBgrEvt_dw[0][6]+nBgrEvt_dw[1][6]+nBgrEvt_dw[0][5]+nBgrEvt_dw[1][5])/(tot_nBgrEvt[6]+tot_nBgrEvt[5])))/2.; //sumQCDe+QCDmu
  // double QCD_impact=bgr_impact[5];
  double Z_impact=bgr_impact[2]; 
  double Dibo_impact=bgr_impact[4];
  double ttV_impact=bgr_impact[7];
  //double EWK_impact=(fabs(1-(nBgrEvt_up[0][3]+nBgrEvt_up[1][3]+nBgrEvt_up[0][4]+nBgrEvt_up[1][4]+nBgrEvt_up[0][5]+nBgrEvt_up[1][5])/(tot_nBgrEvt[3]+tot_nBgrEvt[4]+tot_nBgrEvt[5])) + fabs(1-(nBgrEvt_dw[0][3]+nBgrEvt_dw[1][3]+nBgrEvt_dw[0][4]+nBgrEvt_dw[1][4]+nBgrEvt_dw[0][5]+nBgrEvt_dw[1][5])/(tot_nBgrEvt[3]+tot_nBgrEvt[4]+tot_nBgrEvt[5])))/2.;


  if(doDebug) {
    cout << "Impact  of " << syst << " on background" << endl;
    cout << "Channel  nom       up       down     av. impact (%)" << endl;
    cout << "ele      " << nEvt[0] << "   " << nEvt_up[0]<< "  " 
	 << nEvt_dw[0] 
	 << "  " << 100.*(fabs(1-nEvt_up[0]/nEvt[0]) + fabs(1-nEvt_dw[0]/nEvt[0]))/2. 
	 << endl;
    cout << "muo      " << nEvt[1] << "   " << nEvt_up[1]<< "  " 
	 << nEvt_dw[1] 
	 << "  " << 100.*(fabs(1-nEvt_up[1]/nEvt[1]) + fabs(1-nEvt_dw[1]/nEvt[1]))/2. 
	 << endl;
    cout << "tot      " << nEvt[0]+nEvt[1] << "   " << nEvt_up[0]+nEvt_up[1]
	 << "  "  << nEvt_dw[0]+nEvt_dw[1] << "  " 
	 << 100*tot_impact  << endl<<endl;

    cout << "Impact  of " << syst << " on signal " << signame[spec_sign] << endl;
    cout << "Channel  nom       up       down     av. impact (%)" << endl;
    cout << "ele      " << sig_nEvt[0] << "   " << sig_nEvt_up[0]<< "  " 
	 << sig_nEvt_dw[0] 
	 << "  " << 100.*(fabs(1-sig_nEvt_up[0]/sig_nEvt[0]) + fabs(1-sig_nEvt_dw[0]/sig_nEvt[0]))/2. 
	 << endl;
    cout << "muo      " << sig_nEvt[1] << "   " << sig_nEvt_up[1]<< "  " 
	 << sig_nEvt_dw[1] 
	 << "  " << 100.*(fabs(1-sig_nEvt_up[1]/sig_nEvt[1]) + fabs(1-sig_nEvt_dw[1]/sig_nEvt[1]))/2. 
	 << endl;
    cout << "tot      " << sig_nEvt[0]+sig_nEvt[1] << "   " << sig_nEvt_up[0]+sig_nEvt_up[1]
	 << "  "  << sig_nEvt_dw[0]+sig_nEvt_dw[1] << "  " 
	 << 100*sig_impact  << endl;
  }
  else {
    //Print the short version for the note table
    cout << setiosflags(ios::fixed) << setprecision(1);
    //cout << setiosflags(ios::fixed) << setprecision(5);
    TString longlabel=getSystLongName(sshort);
    cout << longlabel << "& $" 
	 << 100*tot_impact << "$& $"
	 << 100*tt_impact << "$& $" 
	 << 100*stop_impact << "$& $" 
	 << 100*W_impact << "$& $" 
	 << 100*QCD_impact << "$& $" 
      //<< 100*EWK_impact << "$& $" 
	 << 100*Z_impact << "$& $"
	 << 100*Dibo_impact << "$& $" 
	 << 100*ttV_impact<< "$& $"
      	 << 100*sig_impact << "$";
    if(doPrintRelLim) {
      string filename=rellim_path+sshort+"u.dat";
      if(doDebug) cout << filename << endl;
      ifstream str_rel(string(filename).c_str());
      if(!str_rel.is_open()) cout << "ERROR reading file" << endl;
      int j=0;
      double ii, rel_up, rel_dw;
      double all_rel_up[NZpMASSES], all_rel_dw[NZpMASSES];
      while(!str_rel.eof() && j < NZpMASSES ) {
	str_rel >> ii >> rel_up >> rel_dw;
	all_rel_up[j]=rel_up;
	all_rel_dw[j]=rel_dw;
	++j;
      }
      cout << "& $" << all_rel_up[spec_sign] << "$&  $" << all_rel_dw[spec_sign] << "$ ";
    }
    cout << "\\\\";
    cout << endl;
  }
  if(doSubBkg) {
    double tot_Bgrimpact;
    for(int b=0;b<NBGR;++b) {
      tot_Bgrimpact=(fabs(1-tot_nBgrEvt_up[b]/tot_nBgrEvt[b]) + fabs(1-tot_nBgrEvt_dw[b]/tot_nBgrEvt[b]))/2.;

      if(doDebug) {
	cout << "Impact  of " << syst << " on " << bgrname[b] << " background" << endl;
	cout << "Channel  nom       up       down     av. impact (%)" << endl;
	cout << "ele      " << nBgrEvt[0][b] << "   " << nBgrEvt_up[0][b]<< "  " 
	     << nBgrEvt_dw[0][b] 
	     << "  " << 100.*(fabs(1-nBgrEvt_up[0][b]/nBgrEvt[0][b]) + fabs(1-nBgrEvt_dw[0][b]/nBgrEvt[0][b]))/2. 
	     << endl;
	cout << "muo      " << nBgrEvt[1][b] << "   " << nBgrEvt_up[1][b]<< "  " 
	     << nBgrEvt_dw[1][b] 
	     << "  " << 100.*(fabs(1-nBgrEvt_up[1][b]/nBgrEvt[1][b]) + fabs(1-nBgrEvt_dw[1][b]/nBgrEvt[1][b]))/2. 
	     << endl;
	cout << "tot      " << nBgrEvt[0][b]+nBgrEvt[1][b] << "   " << nBgrEvt_up[0][b]+nBgrEvt_up[1][b]
	     << "  "  << nBgrEvt_dw[0][b]+nBgrEvt_dw[1][b] << "  " 
	     << 100*tot_Bgrimpact  << endl;
      }//end if(doDebug)
    }//end lop over backgrounds
  }//end if(doSubBgr)
}
