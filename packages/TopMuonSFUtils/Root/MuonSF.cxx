#include "TopMuonSFUtils/MuonSF.h"
    
using namespace std;

using namespace TopMuon;

MuonSF::MuonSF(const TString& root_data_directory):
  PeriodAndIndex_muID(), PeriodAndIndex_muTrigEffSF(), TH1DVector_muID(), TH1DVector_muID_staterr(), TH1DVector_muID_systerr(),
  TH1DVector_muID_toterr(), TH2DVector_muTrigEffSF_barrel1(), TH2DVector_muTrigEffSF_barrel1_staterr(),
  TH2DVector_muTrigEffSF_barrel1_systerr(), TH2DVector_muTrigEffSF_barrel1_toterr(), TH2DVector_muTrigEffData_barrel1(),
  TH2DVector_muTrigEffMC_barrel1(), TH2DVector_muTrigEffSF_endcap1(), TH2DVector_muTrigEffSF_endcap1_staterr(),
  TH2DVector_muTrigEffSF_endcap1_systerr(), TH2DVector_muTrigEffSF_endcap1_toterr(), TH2DVector_muTrigEffData_endcap1(),
  TH2DVector_muTrigEffMC_endcap1(), TH2DVector_muTrigEffSF_barrel2(), TH2DVector_muTrigEffSF_barrel2_staterr(),
  TH2DVector_muTrigEffSF_barrel2_systerr(), TH2DVector_muTrigEffSF_barrel2_toterr(), TH2DVector_muTrigEffData_barrel2(),
  TH2DVector_muTrigEffMC_barrel2(), TH2DVector_muTrigEffSF_endcap2(), TH2DVector_muTrigEffSF_endcap2_staterr(),
  TH2DVector_muTrigEffSF_endcap2_systerr(), TH2DVector_muTrigEffSF_endcap2_toterr(), TH2DVector_muTrigEffData_endcap2(),
  TH2DVector_muTrigEffMC_endcap2(), ifile_muID(0), ifile_muTrigEff_barrel1(0), ifile_muTrigEff_endcap1(0) {
  std::cout << "MuonSF : Start reading muon scale maps from " << root_data_directory << std::endl;
  ifile_muID = TFile::Open(root_data_directory+"/muID_Summary.root","read");
  if(!ifile_muID) cerr << "Did not open " << root_data_directory+"/muID_Summary.root" << endl;
  assert(ifile_muID);
  for(int i = 0 ; i < ifile_muID->GetListOfKeys()->GetEntries() ; i++){
    TString objname = ifile_muID->GetListOfKeys()->At(i)->GetName();
    if(!objname.Contains("sfOverall_period_"))continue;
    objname.Remove(0,objname.Index("period")+7);
    TH1D* hsf         = (TH1D*)ifile_muID->Get("sfOverall_period_"+objname);
    TH1D* hsf_staterr = (TH1D*)ifile_muID->Get("sfOverall_stat_period_"+objname);
    TH1D* hsf_systerr = (TH1D*)ifile_muID->Get("sfOverall_syst_period_"+objname);
    TH1D* hsf_toterr  = (TH1D*)ifile_muID->Get("sfOverall_toterr_period_"+objname);
    assert(hsf        );
    assert(hsf_staterr);
    assert(hsf_systerr);
    assert(hsf_toterr );
    TH1DVector_muID.push_back(hsf);
    TH1DVector_muID_staterr.push_back(hsf_staterr);
    TH1DVector_muID_systerr.push_back(hsf_systerr);
    TH1DVector_muID_toterr.push_back(hsf_toterr);
    if( objname.CompareTo("BtoD")==0 ) {
      PeriodAndIndex_muID[period2011BtoD] = TH1DVector_muID.size() - 1;
    } else if( objname.CompareTo("EtoH")==0 ) {
      PeriodAndIndex_muID[period2011EtoH] = TH1DVector_muID.size() - 1;
    } else if( objname.CompareTo("I")==0 ) {
      PeriodAndIndex_muID[period2011I] = TH1DVector_muID.size() - 1;
    } else if( objname.CompareTo("JtoK")==0 ) {
      PeriodAndIndex_muID[period2011JtoK] = TH1DVector_muID.size() - 1;
    } else if( objname.CompareTo("L_badRPC")==0 ) {
      PeriodAndIndex_muID[period2011L_badrpc] = TH1DVector_muID.size() - 1;
    } else if( objname.CompareTo("L_goodRPC")==0 ) {
      PeriodAndIndex_muID[period2011L_goodrpc] = TH1DVector_muID.size() - 1;
    } else if( objname.CompareTo("M")==0 ) {
      PeriodAndIndex_muID[period2011M] = TH1DVector_muID.size() - 1;
    } else if( objname.CompareTo("BtoI")==0 ) {
      PeriodAndIndex_muID[period2011BtoI] = TH1DVector_muID.size() - 1;
    } else if( objname.CompareTo("LtoM")==0 ) {
      PeriodAndIndex_muID[period2011LtoM] = TH1DVector_muID.size() - 1;
    } else {
      cerr << "Unrecognised ID histogram " << objname << endl;
    }
    
  }//loop on mu ID file objects
  
  ifile_muTrigEff_barrel1 = TFile::Open(root_data_directory+"/muTrigEff_Summary_Matsushita2_barrelBinning.root","read"); 
  ifile_muTrigEff_endcap1 = TFile::Open(root_data_directory+"/muTrigEff_Summary_Matsushita2_endcapBinning.root","read");
  if(!ifile_muTrigEff_barrel1) cerr << "Did not open " << root_data_directory+"/muTrigEff_Summary_Matsushita2_barrelBinning.root" << endl;
  if(!ifile_muTrigEff_endcap1) cerr << "Did not open " << root_data_directory+"/muTrigEff_Summary_Matsushita2_endcapBinning.root" << endl;
  assert(ifile_muTrigEff_barrel1);
  assert(ifile_muTrigEff_endcap1);
  for(int i = 0 ; i < ifile_muTrigEff_barrel1->GetListOfKeys()->GetEntries() ; i++){
    TString objname = ifile_muTrigEff_barrel1->GetListOfKeys()->At(i)->GetName();
    if(!objname.Contains("sf_period_"))continue;
    objname.Remove(0,objname.Index("period")+7);
    if(1){
      TH2D* hsf         = (TH2D*)ifile_muTrigEff_barrel1->Get("sf_period_"+objname);
      TH2D* hsf_staterr = (TH2D*)ifile_muTrigEff_barrel1->Get("sf_stat_period_"+objname);
      TH2D* hsf_systerr = (TH2D*)ifile_muTrigEff_barrel1->Get("sf_syst_period_"+objname);
      TH2D* hsf_toterr  = (TH2D*)ifile_muTrigEff_barrel1->Get("sf_toterr_period_"+objname);
      assert(hsf        );
      assert(hsf_staterr);
      assert(hsf_systerr);
      assert(hsf_toterr );
      TH2DVector_muTrigEffSF_barrel1.push_back(hsf);
      TH2DVector_muTrigEffSF_barrel1_staterr.push_back(hsf_staterr);
      TH2DVector_muTrigEffSF_barrel1_systerr.push_back(hsf_systerr);
      TH2DVector_muTrigEffSF_barrel1_toterr.push_back(hsf_toterr);
      TH2D* heffdata = (TH2D*)ifile_muTrigEff_barrel1->Get("eff_data_period_"+objname);
      TH2D* heffmc   = (TH2D*)ifile_muTrigEff_barrel1->Get("eff_mc_period_"+objname);
      assert(heffdata);
      assert(heffmc  );
      TH2DVector_muTrigEffData_barrel1.push_back(heffdata);
      TH2DVector_muTrigEffMC_barrel1.push_back(heffmc);
    }
    if(1){
      TH2D* hsf         = (TH2D*)ifile_muTrigEff_endcap1->Get("sf_period_"+objname);
      TH2D* hsf_staterr = (TH2D*)ifile_muTrigEff_endcap1->Get("sf_stat_period_"+objname);
      TH2D* hsf_systerr = (TH2D*)ifile_muTrigEff_endcap1->Get("sf_syst_period_"+objname);
      TH2D* hsf_toterr  = (TH2D*)ifile_muTrigEff_endcap1->Get("sf_toterr_period_"+objname);
      assert(hsf        );
      assert(hsf_staterr);
      assert(hsf_systerr);
      assert(hsf_toterr );
      TH2DVector_muTrigEffSF_endcap1.push_back(hsf);
      TH2DVector_muTrigEffSF_endcap1_staterr.push_back(hsf_staterr);
      TH2DVector_muTrigEffSF_endcap1_systerr.push_back(hsf_systerr);
      TH2DVector_muTrigEffSF_endcap1_toterr.push_back(hsf_toterr);
      TH2D* heffdata = (TH2D*)ifile_muTrigEff_endcap1->Get("eff_data_period_"+objname);
      TH2D* heffmc   = (TH2D*)ifile_muTrigEff_endcap1->Get("eff_mc_period_"+objname);
      assert(heffdata);
      assert(heffmc  );
      TH2DVector_muTrigEffData_endcap1.push_back(heffdata);
      TH2DVector_muTrigEffMC_endcap1.push_back(heffmc);
    }
    if( objname.CompareTo("BtoD")==0 ) {
      PeriodAndIndex_muTrigEffSF[period2011BtoD] = TH2DVector_muTrigEffSF_barrel1.size() - 1;
    } else if( objname.CompareTo("EtoH")==0) {
      PeriodAndIndex_muTrigEffSF[period2011EtoH] = TH2DVector_muTrigEffSF_barrel1.size() - 1;
    } else if( objname.CompareTo("I")==0) {
      PeriodAndIndex_muTrigEffSF[period2011I] = TH2DVector_muTrigEffSF_barrel1.size() - 1;
    } else if( objname.CompareTo("JtoK")==0) {
      PeriodAndIndex_muTrigEffSF[period2011JtoK] = TH2DVector_muTrigEffSF_barrel1.size() - 1;
    } else if( objname.CompareTo("L_badRPC")==0) {
      PeriodAndIndex_muTrigEffSF[period2011L_badrpc] = TH2DVector_muTrigEffSF_barrel1.size() - 1;
    } else if( objname.CompareTo("L_goodRPC")==0) {
      PeriodAndIndex_muTrigEffSF[period2011L_goodrpc] = TH2DVector_muTrigEffSF_barrel1.size() - 1;
    } else if( objname.CompareTo("M")==0) {
      PeriodAndIndex_muTrigEffSF[period2011M] = TH2DVector_muTrigEffSF_barrel1.size() - 1;
    } else if( objname.CompareTo("BtoI")==0) {
      PeriodAndIndex_muTrigEffSF[period2011BtoI] = TH2DVector_muTrigEffSF_barrel1.size() - 1;
    } else if( objname.CompareTo("LtoM")==0) {
      PeriodAndIndex_muTrigEffSF[period2011LtoM] = TH2DVector_muTrigEffSF_barrel1.size() - 1;
    } else {
      cerr << "Unrecognised trigger histogram " << objname << endl;
    }
  }
  assert(PeriodAndIndex_muID.size()>0 && PeriodAndIndex_muTrigEffSF.size()>0);
  std::cout << "MuonSF : Successful load : ID " << PeriodAndIndex_muID.size() << ", Trig " << PeriodAndIndex_muTrigEffSF.size() << std::endl;

}
DataPeriods MuonSF::getPeriod( const int& runNumber ) {
  if( runNumber >= 177986 && runNumber <= 186493 ) return period2011BtoI;
  if( runNumber >= 186516 && runNumber <= 187815 ) return period2011JtoK;
  if( runNumber >= 188902 && runNumber <= 191933 ) return period2011LtoM;
  
  
  /* for more detailed breakdown:
  if( runNumber >= 177986 && runNumber <= 180481 ) return period2011BtoD;
  if( runNumber >= 180614 && runNumber <= 184169 ) return period2011EtoH;
  if( runNumber >= 185353 && runNumber <= 186493 ) return period2011I;
  if( runNumber >= 186516 && runNumber <= 187815 ) return period2011JtoK;
  if( runNumber >= 190503 && runNumber <= 191933 ) return period2011M;
  if( runNumber >= 188902 && runNumber <= 190343 ) {
    if( runNumber >= 189207 && runNumber <= 189610 ) return period2011L_badRPC;
    else                                             return period2011L_goodRPC;
  }
  */
  
  std::cout << "MuonSF::getPeriod : unknown run " << runNumber << std::endl;
  assert(0);
  return periodUnknown;
}
double MuonSF::mu_ID_SF(const TopMuon::DataPeriods& period) {
  assert(PeriodAndIndex_muID.count(period)==1);
  double sf = TH1DVector_muID.at( PeriodAndIndex_muID[ period ] )->GetBinContent(1);
  return sf;
}
double MuonSF::mu_ID_SF_err(const TopMuon::DataPeriods& period) {
  assert( PeriodAndIndex_muID.count(period)==1);
  return TH1DVector_muID_toterr.at( PeriodAndIndex_muID[ period ] )->GetBinContent(1);
}

double MuonSF::mu_trigger_SF( const double& eta,  const double& phi,  const TopMuon::DataPeriods& period) {
  assert(PeriodAndIndex_muTrigEffSF.count(period)==1);
  TH2D* h = (fabs(eta)<1.05)
    ? TH2DVector_muTrigEffSF_barrel1.at( PeriodAndIndex_muTrigEffSF[ period ] )
    : TH2DVector_muTrigEffSF_endcap1.at( PeriodAndIndex_muTrigEffSF[ period ] );
  double sf = h->GetBinContent( h->GetXaxis()->FindBin(eta),h->GetYaxis()->FindBin(phi));
  return sf;
}
double MuonSF::mu_trigger_SF_err( const double& eta, const double& phi, const TopMuon::DataPeriods& period ) {
  assert(PeriodAndIndex_muTrigEffSF.count(period)==1);
  TH2D* htoterr = (fabs(eta)<1.05)
    ? TH2DVector_muTrigEffSF_barrel1_toterr.at( PeriodAndIndex_muTrigEffSF[ period ] )
    : TH2DVector_muTrigEffSF_endcap1_toterr.at( PeriodAndIndex_muTrigEffSF[ period ] );
  double sf_toterr = htoterr->GetBinContent( htoterr->GetXaxis()->FindBin(eta),htoterr->GetYaxis()->FindBin(phi));
  return sf_toterr;
}

double MuonSF::mu_trigger_eff_data( const double& eta,  const double& phi, const TopMuon::DataPeriods& period) {
  assert(PeriodAndIndex_muTrigEffSF.count(period)==1);
  TH2D* h = (fabs(eta)<1.05)
    ? TH2DVector_muTrigEffData_barrel1.at( PeriodAndIndex_muTrigEffSF[ period ] )
    : TH2DVector_muTrigEffData_endcap1.at( PeriodAndIndex_muTrigEffSF[ period ] );
  return h->GetBinContent( h->GetXaxis()->FindBin(eta),h->GetYaxis()->FindBin(phi) );
}
double MuonSF::mu_trigger_eff_mc( const double& eta,  const double& phi, const TopMuon::DataPeriods& period) {
  assert(PeriodAndIndex_muTrigEffSF.count(period)==1);
  TH2D* h = (fabs(eta)<1.05)
    ? TH2DVector_muTrigEffMC_barrel1.at( PeriodAndIndex_muTrigEffSF[ period ] )
    : TH2DVector_muTrigEffMC_endcap1.at( PeriodAndIndex_muTrigEffSF[ period ] );
  return h->GetBinContent( h->GetXaxis()->FindBin(eta),h->GetYaxis()->FindBin(phi) );
}

/**
 * Destructor - close the input files
 */
MuonSF::~MuonSF() {
  if(ifile_muID) ifile_muID->Close();
  if(ifile_muTrigEff_barrel1) ifile_muTrigEff_barrel1->Close();
  if(ifile_muTrigEff_endcap1) ifile_muTrigEff_endcap1->Close();
}
