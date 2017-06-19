//##############################################################################
//##                                                                          ##
//## TTBarSystematics.C script to plot systematics histograms                 ##
//##                                                                          ##
//## author  : Tatjana Lenz (tlenz@cern.ch)                                   ##
//##           Tobias Heck (theck@cern.ch)                                    ##
//## version : 1.2                                                            ##
//## date    : 2012-10-02                                                     ##
//##                                                                          ##
//##############################################################################


#include <iostream>
#include <iomanip>
#include "TROOT.h"
#include "TPad.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include <vector>

void SetATLASStyle();

using namespace std;

//===========================================================================================================
// Main function
//
int main(int argc, char**argv) {

  cout << endl;
  cout << " =================================== " << endl;
  cout << "    Plot TTBar Systematics Histos    " << endl;
  cout << " =================================== " << endl;  
  cout << endl;

  SetATLASStyle();


  // ================================================================
  //                       In/Output settings
  // ================================================================

  // TString inputdir  = "/afs/physik.uni-mainz.de/home/weina/workdir/TopPhys/top_statistics_rs/TopResonance/LPC";
//   TString inputdir  = "/data/atlas/atlasdata2/behr/TTbarResonanceSearch/LimitSetting/outputs/Spectra15July2014/";
  TString inputdir  = "/data/atlas/atlasdata2/behr/TTbarResonanceSearch/LimitSetting/outputs/Spectra06October2014/";
  TString outputdir = "plots";

  TString rootfile_smoothed  = "MttSpectra_Oct06_smoothed.root";
  TString rootfile_unsmoothed  = "MttSpectra_Oct06_unsmoothed.root";
//    TString rootfile  = "PreSpectra_2014Jul14.root";
  //  TString rootfile  = "23_02_2013_glgw_Systematics.root";
//   TString outputdir = "/afs/physik.uni-mainz.de/home/theck/public_html/TTBarSystematics/range/2012_11_24";
//   TString rootfile  = "2012_11_24_glgw_Systematics.root";

  // ================================================================
  //                         Other settings
  // ================================================================
   bool plotSmooth=true;

  
  TString Lumi     = "20.4";
  Bool_t  doPrelim = false;
  // int     rebin    = 20;
  bool    doEPS    = true;
  bool    doPNG    = true;
  bool    doCheck  = false;
  //if (doEPS) outputdir = "/afs/physik.uni-mainz.de/home/weina/public_html/TTBarSystematics/20130407/eps";
  gStyle->SetLabelFont(42, "xyz");

  // ================================================================
  //                           Systematics
  // ================================================================

  //Int_t nSys = 47;//61;
  //TString systematics[nSys];

  // systematics[0] = "PDF";
  vector<TString>  systematics;
  
//   systematics.push_back("ElectronSF");
//   systematics.push_back("MuonSF");
  systematics.push_back("EleSF");
  systematics.push_back("MuSF");
  systematics.push_back("Btag"); 
  systematics.push_back("BtagC");
  systematics.push_back("BtagL");

  systematics.push_back("norm_tt");
//   systematics.push_back("Norm_ttbar");
  systematics.push_back("PDF");
  systematics.push_back("JVFCut");
  systematics.push_back("luminosity");
  systematics.push_back("IFSR");
  systematics.push_back("MCGen");
//   systematics.push_back("MC_Gen");
  systematics.push_back("PartonShower");
//   systematics.push_back("PS");
  systematics.push_back("EWS");
  systematics.push_back("topmass");
//   systematics.push_back("top_mass");
  systematics.push_back("hdamp");

  systematics.push_back("iqopt3");
  systematics.push_back("ptjmin10");
  systematics.push_back("Whfsf");
  
  systematics.push_back("norm_QCDe");
  systematics.push_back("norm_QCDmu");

  systematics.push_back("BoostedJES_ALL");
  systematics.push_back("BoostedJESothers");
  systematics.push_back("BoostedJES0");
  systematics.push_back("BoostedJES1");
  systematics.push_back("BoostedJES2");
  systematics.push_back("BoostedJES3");
  systematics.push_back("BoostedJES4");
  systematics.push_back("BoostedJES5");
  systematics.push_back("BoostedJES6");
  systematics.push_back("BoostedJES7");
  systematics.push_back("BoostedJES8");
  systematics.push_back("BoostedJES9");
  systematics.push_back("BoostedJES10");
  systematics.push_back("BoostedJES11");
  systematics.push_back("BoostedJES12");
  systematics.push_back("BoostedJES13");
  systematics.push_back("BoostedJES14");
  systematics.push_back("BoostedJES15");
  systematics.push_back("BoostedJES16");
  systematics.push_back("BoostedJMS");
  systematics.push_back("BoostedJMR");
  systematics.push_back("BoostedJER");

  systematics.push_back("JES_ALL");  
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
  systematics.push_back("JES14");
  systematics.push_back("JES15");
  systematics.push_back("JES16");
  systematics.push_back("JES17");
  systematics.push_back("JES18");
//   systematics.push_back("JES19");
  systematics.push_back("JES20");
  systematics.push_back("JES21");
  systematics.push_back("JES22");
  systematics.push_back("JESsmall");
//   systematics.push_back("smallJES");
//   systematics.push_back("JER");
//   systematics.push_back("JVF");
  systematics.push_back("JetEnerRes");

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
//   systematics.push_back("BtagC0");
//   systematics.push_back("BtagC1");
//   systematics.push_back("BtagC2");
//   systematics.push_back("BtagC3");
//   systematics.push_back("BtagC4");
//   systematics.push_back("BtagC5");
//   systematics.push_back("BtagL0");
//   systematics.push_back("BtagL1");
//   systematics.push_back("BtagL2");
//   systematics.push_back("BtagL3");
//   systematics.push_back("BtagL4");
//   systematics.push_back("BtagL5");
//   systematics.push_back("BtagL6");
//   systematics.push_back("BtagL7");
//   systematics.push_back("BtagL8");
//   systematics.push_back("BtagL9");
//   systematics.push_back("BtagL10");
//   systematics.push_back("BtagL11");
//   systematics.push_back("smallBtag");
  
//   systematics.push_back("btag6");
//   systematics.push_back("btag7");
//   systematics.push_back("btag8");
//   systematics.push_back("btag9");
//   systematics.push_back("btag10");
//   systematics.push_back("ctag");
//   systematics.push_back("mistag");

  Int_t nSys=systematics.size();
    
 

  // ================================================================
  //                             Channels
  // ================================================================
  
  vector<TString> channels;
  /*channels.push_back("masstT_cat1_mu");
  channels.push_back("masstT_cat1_e");
  channels.push_back("masstT_cat2_mu");
  channels.push_back("masstT_cat2_e");
  channels.push_back("masstT_cat3_mu");
  channels.push_back("masstT_cat3_e");
  channels.push_back("massTTbarChi2LPC_cat1_mu");
  channels.push_back("massTTbarChi2LPC_cat1_e");
  channels.push_back("massTTbarChi2LPC_cat2_mu");
  channels.push_back("massTTbarChi2LPC_cat2_e");
  channels.push_back("massTTbarChi2LPC_cat3_mu");
  channels.push_back("massTTbarChi2LPC_cat3_e");*/

  channels.push_back("masstT_mu");
  channels.push_back("masstT_e");
  channels.push_back("massTTbarChi2LPC_mu");
  channels.push_back("massTTbarChi2LPC_e");

  Int_t nChan = channels.size();
  
  // ================================================================
  //                             Samples
  // ================================================================

  vector<TString> samples;
  samples.push_back("Bgr");
  samples.push_back("tt");
  samples.push_back("W");
  samples.push_back("smallBgr");

//   samples.push_back("binkkg_Bgr");
//   samples.push_back("binkkg_tt");
//   samples.push_back("binkkg_W");
//   samples.push_back("binkkg_smallBgr");
//samples.push_back("binkkg_Z");  
  //samples.push_back("binkkg_Diboson");
  //samples.push_back("binkkg_single-top");
  //samples.push_back("binkkg_QCDe");  
  //samples.push_back("binkkg_QCDmu");
  //samples.push_back("binkkg_ttV");
  //samples.push_back("binkkg_Z750");
  //samples.push_back("binkkg_Z1000");
  //samples.push_back("binkkg_Z1500");  
  //samples.push_back("binkkg_Z2000");

  Int_t nSamp = samples.size();
  
  
  // ================================================================
  //                          Read histograms
  // ================================================================

  TString filename_smoothed = inputdir+rootfile_smoothed;
  cout<<"filename_smoothed "<<filename_smoothed<<endl;
  TFile *file_smoothed = TFile::Open(filename_smoothed);
  
  TString filename_unsmoothed = inputdir+rootfile_unsmoothed;
  cout<<"filename_unsmoothed "<<filename_unsmoothed<<endl;
  TFile *file_unsmoothed = TFile::Open(filename_unsmoothed);
  
  if(!file_smoothed) return 0;
  if(!file_unsmoothed) return 0;
  
  //loop samples  nSamp
  for(unsigned int isamp = 0; isamp <nSamp; isamp++) {

    cout << samples[isamp] << endl;

    //loop channels
    for(unsigned int ichan = 0; ichan < nChan; ichan++) {

      cout << "  " << channels[ichan] << endl;

      //get nominal histogram
      TH1F *h1p_hist = (TH1F*)file_smoothed->Get(samples[isamp]+"_"+channels[ichan]);
      if(!h1p_hist) {
	cout << " WARNING :: Histogram " << samples[isamp]+"_"+channels[ichan] << " not available. Skipping..." << endl;
	continue;
      }

      //set nominal histogram parameters
      TH1F h1_hist = *h1p_hist;
      h1_hist.SetTitle("");
      h1_hist.SetLineStyle(1);
      h1_hist.SetLineWidth(2);
      h1_hist.SetLineColor(kBlack); 
      h1_hist.SetFillColor(kWhite);       
      
      //create nominal error histogram
      TH1F h1_hist_err = *((TH1F*)h1_hist.Clone());
      //h1_hist_err.SetFillStyle(3244);
      h1_hist_err.SetFillColor(33);
      h1_hist_err.SetLineColor(33);

      //loop systematics
       for(unsigned int isys = 0; isys < nSys; isys++) {

	 // for(unsigned int isys = 35; isys < 36; isys++) {
	//cout << "    " << systematics[isys] << endl;

	//get variation histograms
	TH1F *h1p_var_up   = (TH1F*)file_smoothed->Get(samples[isamp]+"_"+channels[ichan]+"_"+systematics[isys]+"_up");
	TH1F *h1p_var_down = (TH1F*)file_smoothed->Get(samples[isamp]+"_"+channels[ichan]+"_"+systematics[isys]+"_dw");

	TH1F *h1p_unsmoothed_var_up   = (TH1F*)file_unsmoothed->Get(samples[isamp]+"_"+channels[ichan]+"_"+systematics[isys]+"_up");
	TH1F *h1p_unsmoothed_var_down = (TH1F*)file_unsmoothed->Get(samples[isamp]+"_"+channels[ichan]+"_"+systematics[isys]+"_dw");
	
	//Check availability
	if(!h1p_var_up || !h1p_var_down) {
	  cout << " WARNING :: Missing variation histograms for " << samples[isamp]+"_"+channels[ichan]+"_"+systematics[isys]+". Skipping..." << endl;
	  continue;
	}

	if(!h1p_unsmoothed_var_up || !h1p_unsmoothed_var_down) {
	  cout << " WARNING :: Missing variation histograms for " << samples[isamp]+"_"+channels[ichan]+"_"+systematics[isys]+". Skipping..." << endl;
	  continue;
	}
	
	//set var up histogram parameters
	TH1F h1_var_up = *h1p_var_up;
	h1_var_up.SetTitle("");
	h1_var_up.SetLineStyle(1);
	h1_var_up.SetLineWidth(2);
	h1_var_up.SetLineColor(kRed+2);
	h1_var_up.SetLineStyle(1);
	h1_var_up.SetFillColor(kWhite);

	TH1F h1_unsmoothed_var_up = *h1p_unsmoothed_var_up;
	h1_unsmoothed_var_up.SetTitle("");
	h1_unsmoothed_var_up.SetLineStyle(1);
	h1_unsmoothed_var_up.SetLineWidth(2);
	h1_unsmoothed_var_up.SetLineColor(kRed+2);
	h1_unsmoothed_var_up.SetLineStyle(2);
	h1_unsmoothed_var_up.SetFillColor(kWhite);
	
	//set var down histogram parameters
	TH1F h1_var_down = *h1p_var_down;
	h1_var_down.SetTitle("");
	h1_var_down.SetLineStyle(1);
	h1_var_down.SetLineWidth(2);
	h1_var_down.SetLineColor(kBlue+2);
	h1_var_down.SetLineStyle(1); 
	h1_var_down.SetFillColor(kWhite);

	TH1F h1_unsmoothed_var_down = *h1p_unsmoothed_var_down;
	h1_unsmoothed_var_down.SetTitle("");
	h1_unsmoothed_var_down.SetLineStyle(1);
	h1_unsmoothed_var_down.SetLineWidth(2);
	h1_unsmoothed_var_down.SetLineColor(kBlue+2);
	h1_unsmoothed_var_down.SetLineStyle(3); 
	h1_unsmoothed_var_down.SetFillColor(kWhite);
	
	//creatre relative difference var up histogram
	TH1F h1_var_up_ratio = *((TH1F*)h1_var_up.Clone());
	h1_var_up_ratio.Add(&h1_hist, -1); 
	h1_var_up_ratio.Divide(&h1_hist);
	h1_var_up_ratio.SetLineColor(kRed+2);
	h1_var_up_ratio.SetLineStyle(1);
	h1_var_up_ratio.SetFillColor(kWhite);

	TH1F h1_unsmoothed_var_up_ratio = *((TH1F*)h1_unsmoothed_var_up.Clone());
	h1_unsmoothed_var_up_ratio.Add(&h1_hist, -1); 
	h1_unsmoothed_var_up_ratio.Divide(&h1_hist);
	h1_unsmoothed_var_up_ratio.SetLineColor(kRed+2);
	h1_unsmoothed_var_up_ratio.SetLineStyle(2);
	h1_unsmoothed_var_up_ratio.SetFillColor(kWhite);
	
	//creatre relative difference var down histogram	
	TH1F h1_var_down_ratio = *((TH1F*)h1_var_down.Clone());
	h1_var_down_ratio.Add(&h1_hist, -1);
	h1_var_down_ratio.Divide(&h1_hist); 
	h1_var_down_ratio.SetLineColor(kBlue+2);
	h1_var_down_ratio.SetLineStyle(1);   
	h1_var_down_ratio.SetFillColor(kWhite);

	TH1F h1_unsmoothed_var_down_ratio = *((TH1F*)h1_unsmoothed_var_down.Clone());
	h1_unsmoothed_var_down_ratio.Add(&h1_hist, -1);
	h1_unsmoothed_var_down_ratio.Divide(&h1_hist); 
	h1_unsmoothed_var_down_ratio.SetLineColor(kBlue+2);
	h1_unsmoothed_var_down_ratio.SetLineStyle(3);   
	h1_unsmoothed_var_down_ratio.SetFillColor(kWhite);
	
	//Some checks if histograms are properly filled
	if(doCheck) {
	  if(h1_hist.GetMaximum()     == 0) cout << " WARNING :: Probably non filled reference histogram (max = 0)"      << endl;
	  if(h1_var_up.GetMaximum()   == 0) cout << " WARNING :: Probably non filled up-variation histogram (max = 0)"   << endl;
	  if(h1_var_down.GetMaximum() == 0) cout << " WARNING :: Probably non filled down-variation histogram (max = 0)" << endl;
	  if(h1_var_up_ratio.GetMaximum()   == 0 || h1_var_up_ratio.GetMinimum()   == 0) cout << " WARNING :: Reference and up-variation histograms probably equal" << endl;
	  if(h1_var_down_ratio.GetMaximum() == 0 || h1_var_down_ratio.GetMinimum() == 0) cout << " WARNING :: Reference and down-variation histograms probably equal" << endl;
	}

	//create error ratio histogram
	double maxerr = 0.0;	
	TH1F h1_hist_err_ratio = *((TH1F*)h1_hist.Clone());
	//h1_hist_err_ratio.SetFillStyle(3244);
	h1_hist_err_ratio.SetFillColor(33);

	//loop all bins
	for(unsigned int ibin = 0; ibin <= h1_hist_err_ratio.GetNbinsX(); ibin++) {
	  double binerr = h1_hist_err_ratio.GetBinError(ibin);
	  double bincon = h1_hist_err_ratio.GetBinContent(ibin);
	  double newbinerr = 0.0;
	  if(bincon > 0) newbinerr = binerr/bincon;
	  h1_hist_err_ratio.SetBinContent(ibin, 0);
	  h1_hist_err_ratio.SetBinError(ibin, newbinerr);
	  if(newbinerr > maxerr) maxerr = newbinerr;
	}

	
	//Create template histogram
	TH1F h1_tpl = *((TH1F*)h1_hist.Clone());
	h1_tpl.Reset();
	h1_tpl.SetTitle("; ; Entries; ");
	h1_tpl.GetXaxis()->SetLabelSize(0.0000001);
	h1_tpl.GetYaxis()->SetLabelSize(0.04);
	h1_tpl.GetYaxis()->SetTitleOffset(1.5);
	h1_tpl.SetMaximum(1.3*h1_var_up.GetMaximum());
	//h1_tpl.SetMaximum(8500);
	//h1_tpl.SetMaximum(220);

	//Create template ratio histogram
	TH1F h1_tpl_ratio = *((TH1F*)h1_hist.Clone());
	h1_tpl_ratio.Reset();
	h1_tpl_ratio.SetTitle("; t#bar{t} mass [GeV] ; relative difference ; ");
	Float_t diffmax = max(h1_var_up_ratio.GetMaximum(), h1_var_down_ratio.GetMaximum());
	Float_t diffmin = min(h1_var_up_ratio.GetMinimum(), h1_var_down_ratio.GetMinimum());

	//calculate y-range of ratio histogram
 	double range = 1.0;
	if(fabs(diffmax) < 0.08 && fabs(diffmin) < 0.08) range = 0.1;
	else if(fabs(diffmax) < 0.18 && fabs(diffmin) < 0.18) range = 0.2;
	else if(fabs(diffmax) < 0.28 && fabs(diffmin) < 0.28) range = 0.3;
	else if(fabs(diffmax) < 0.38 && fabs(diffmin) < 0.38) range = 0.4;
	else if(fabs(diffmax) < 0.48 && fabs(diffmin) < 0.48) range = 0.5;
	else if(fabs(diffmax) < 0.58 && fabs(diffmin) < 0.58) range = 0.6;
	else if(fabs(diffmax) < 0.68 && fabs(diffmin) < 0.68) range = 0.7;
 	else range = 0.8;
	//	range = 0.5;
	//range = 0.2;

	//Set ratio template parameters
	h1_tpl_ratio.SetMaximum(range);
	h1_tpl_ratio.SetMinimum(-range);
	h1_tpl_ratio.GetXaxis()->SetTitleSize(0.11);
	h1_tpl_ratio.GetXaxis()->SetLabelSize(0.11);
	h1_tpl_ratio.GetYaxis()->SetTitleSize(0.11);
	h1_tpl_ratio.GetYaxis()->SetTitleOffset(0.45);
	h1_tpl_ratio.GetYaxis()->SetLabelSize(0.11);
	h1_tpl_ratio.SetTickLength(0.03,"y");
	h1_tpl_ratio.SetTickLength(0.10,"x");
	h1_tpl_ratio.SetTickLength(0.025, "Y");
	h1_tpl_ratio.SetNdivisions(305, "Y");
	h1_tpl.GetYaxis()->SetLabelSize(0.04);
	h1_tpl.GetYaxis()->SetTitleOffset(1.5);


	//insets
	TLatex lumi;
	lumi.SetNDC();
	lumi.SetTextFont(42);
	lumi.SetTextSize(0.032);
	lumi.SetTextColor(1);
	
	TLatex chan;
	chan.SetNDC();
	chan.SetTextFont(42);
	chan.SetTextSize(0.05);
	chan.SetTextColor(1);
	
	TLatex prelim;
	prelim.SetNDC();
	prelim.SetTextFont(42);
	prelim.SetTextColor(1);
      	
	//legend
	TLegend leg(0.70, 0.70, 0.95, 0.95);
	leg.SetBorderSize(0);
	leg.SetTextSize(0.032);
	leg.SetFillColor(0);
	leg.AddEntry(&h1_hist,     " nominal ",      "l");	
	leg.AddEntry(&h1_var_up,   " up (smoothed)",           "l");
	leg.AddEntry(&h1_var_down, " down (smoothed)",         "l");
	leg.AddEntry(&h1_unsmoothed_var_up,   " up (unsmoothed)",           "l");
	leg.AddEntry(&h1_unsmoothed_var_down, " down (unsmoothed)",         "l");
	//	leg.AddEntry(&h1_hist_err, " nominal error", "f");
	leg.AddEntry(&h1_hist_err, " statistical uncertainty", "f");

	
	//canvas and pads

	TCanvas canv("canvas", "", 800, 800);
	
	TPad padTop("padTop", "", 0.0, 0.2, 0.98, 0.98);
	padTop.SetTopMargin(0.01);
	padTop.SetRightMargin(0.01);
	padTop.SetLeftMargin(0.12);
	padTop.SetFillColor(0); 

	TPad padBot("padBot", "", 0.0, 0.03, 0.98, 0.3);
	padBot.SetTopMargin(0.02);
	padBot.SetRightMargin(0.01);
	padBot.SetBottomMargin(0.15);
	padBot.SetLeftMargin(0.12);
	padBot.SetFillColor(0); 
	padBot.SetGridy();		
	
	//Draw
	padTop.Draw();
	padTop.cd();

	h1_tpl.Draw("");
	h1_hist_err.Draw("E2 SAME");
	h1_var_up.Draw("SAME HIST");
	h1_var_down.Draw("SAME HIST");
	h1_unsmoothed_var_up.Draw("SAME HIST");
	h1_unsmoothed_var_down.Draw("SAME HIST");
	h1_hist.Draw("SAME HIST");

	leg.Draw();
	//prelim.DrawLatex(0.27, 0.88, "#font[72]{ATLAS} Preliminary");	
	if(channels[ichan].Contains("_e")) chan.DrawLatex(0.50, 0.90, "e+jets");
	else if(channels[ichan].Contains("_mu")) chan.DrawLatex(0.50, 0.90, "#mu+jets");
	
	lumi.DrawLatex(0.70, 0.60, "#sqrt{s} = 8 TeV");
	lumi.DrawLatex(0.70, 0.52, "#int L dt = "+Lumi+" fb^{-1}");
	lumi.DrawLatex(0.70, 0.38, systematics[isys]);
	lumi.DrawLatex(0.70, 0.33, samples[isamp]);

	padTop.Update();
	padTop.RedrawAxis();

	canv.cd();

	padBot.Draw();
	padBot.cd();
	
	h1_tpl_ratio.Draw();
	h1_hist_err_ratio.Draw("SAME E2");
	h1_var_up_ratio.Draw("SAME HIST");
	h1_var_down_ratio.Draw("SAME HIST");
	h1_unsmoothed_var_up_ratio.Draw("SAME HIST");
	h1_unsmoothed_var_down_ratio.Draw("SAME HIST");
	h1_hist_err_ratio.Draw("SAME HIST"); //redraw without error bars to have 0-line

	padBot.Update();
	padBot.RedrawAxis();     //Redraw axis 
	padBot.RedrawAxis("g");  //Redraw gridlines


	TString out = outputdir+"/"+samples[isamp]+"/"+channels[ichan];
	system("if [ ! -d "+ out + " ]; then mkdir -p " + out + "; fi");

	if(doEPS) canv.Print(out+"/"+samples[isamp]+"_"+channels[ichan]+"_"+systematics[isys]+".eps");
        if(doPNG) canv.Print(out+"/"+samples[isamp]+"_"+channels[ichan]+"_"+systematics[isys]+".png");
      }
    }  
  }

  cout << endl;
  cout << " =================================== " << endl;
  cout << "             Skript done             " << endl;
  cout << " =================================== " << endl;  
  cout << endl;

  return 1;
}


//===========================================================================================================
// Set common ATLAS Style for plots
//
void SetATLASStyle() {

  gStyle->SetOptStat(0);

  gStyle->SetCanvasBorderMode(0); //frame color of canvas
  gStyle->SetCanvasColor(0);  //bkrd color of canvas
  gStyle->SetStatBorderSize(0); //frame style of stat-box 1

  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);

  gStyle->SetLineWidth(2); // width of ticks
  gStyle->SetPadTickX(1); //1:ticks on upper,2: ticks+labels on upper xaxis
  gStyle->SetPadTickY(1);

  gStyle->SetPadLeftMargin(0.18); // 0.18
  gStyle->SetPadRightMargin(0.08);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.17);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetFrameFillColor(0);
  gStyle->SetPaintTextFormat(".2f");
  
  gStyle->SetTitleOffset(1.1,"X");
  gStyle->SetTitleOffset(1.1,"Y");
  gStyle->SetTextFont(42);  
  
}
 


