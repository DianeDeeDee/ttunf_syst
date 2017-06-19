#ifndef _TOPPHYS_MUON_SF_H_
#define _TOPPHYS_MUON_SF_H_
#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <cassert>
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TFile.h"

/**
 * Enumeration for the different data
 * periods we have evaluated efficiencies for.
 * 
 * The current recommendation for mc11a is to use the following grouped efficiency periods:
 * period2011BtoI
 * period2011JtoK
 * period2011LtoM
 */
namespace TopMuon {
  enum DataPeriods{
    periodUnknown=0,
    period2011BtoD=1,
    period2011EtoH=2,
    period2011I=4,
    period2011JtoK=8,
    period2011L_badrpc=16,
    period2011L_goodrpc=32,
    period2011M=64,

    period2011BtoH = (period2011BtoD | period2011EtoH),
    period2011BtoI = (period2011BtoD | period2011EtoH | period2011I),
    period2011LtoM = (period2011L_badrpc | period2011L_goodrpc | period2011M)

  
  };
}//namespace TopMuon

/**
 * Class to access the muon ID & trigger scale factors / efficiencies.
 *
 * authors: Masato Aoki (masato.aoki@cern.ch) and Mark Owen (markowen@cern.ch)
 * 
 * The class should be instantiated like:
 * MuonSF muonsftool("/your/path/MuonSF/TopMuonSFUtils/data/");
 *
 * The scale factors are accessed like:
 * double sf = muonsftool.mu_trigger_SF(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011BtoI );
 * The possible values for the period splitting are detailed in the TopMuon::DataPeriods enum (above).
 *
 * The current recommendation for mc11a is to use the following grouped efficiency periods:
 * period2011BtoI
 * period2011JtoK
 * period2011LtoM
 */
class MuonSF{
 public:
  ~MuonSF();
  MuonSF( const TString& root_data_directory );

  double mu_trigger_SF( const double& eta,  const double& phi,  const TopMuon::DataPeriods& period) ;
  double mu_trigger_SF_err( const double& eta,  const double& phi,  const TopMuon::DataPeriods& period) ;
  double mu_trigger_eff_data( const double& eta,  const double& phi,  const TopMuon::DataPeriods& period) ;
  double mu_trigger_eff_mc( const double& eta,  const double& phi,  const TopMuon::DataPeriods& period) ;  

  double mu_ID_SF( const TopMuon::DataPeriods& period ) ;
  double mu_ID_SF_err( const TopMuon::DataPeriods& period ) ;

  TopMuon::DataPeriods getPeriod(  const int& runNumber ) ;
 private:
  std::map<TopMuon::DataPeriods,int> PeriodAndIndex_muID;
  std::map<TopMuon::DataPeriods,int> PeriodAndIndex_muTrigEffSF;
  std::vector<TH1D*> TH1DVector_muID;
  std::vector<TH1D*> TH1DVector_muID_staterr;
  std::vector<TH1D*> TH1DVector_muID_systerr;
  std::vector<TH1D*> TH1DVector_muID_toterr;

  std::vector<TH2D*> TH2DVector_muTrigEffSF_barrel1;
  std::vector<TH2D*> TH2DVector_muTrigEffSF_barrel1_staterr;
  std::vector<TH2D*> TH2DVector_muTrigEffSF_barrel1_systerr;
  std::vector<TH2D*> TH2DVector_muTrigEffSF_barrel1_toterr;
  std::vector<TH2D*> TH2DVector_muTrigEffData_barrel1;
  std::vector<TH2D*> TH2DVector_muTrigEffMC_barrel1;

  std::vector<TH2D*> TH2DVector_muTrigEffSF_endcap1;
  std::vector<TH2D*> TH2DVector_muTrigEffSF_endcap1_staterr;
  std::vector<TH2D*> TH2DVector_muTrigEffSF_endcap1_systerr;
  std::vector<TH2D*> TH2DVector_muTrigEffSF_endcap1_toterr;
  std::vector<TH2D*> TH2DVector_muTrigEffData_endcap1;
  std::vector<TH2D*> TH2DVector_muTrigEffMC_endcap1;


  std::vector<TH2D*> TH2DVector_muTrigEffSF_barrel2;
  std::vector<TH2D*> TH2DVector_muTrigEffSF_barrel2_staterr;
  std::vector<TH2D*> TH2DVector_muTrigEffSF_barrel2_systerr;
  std::vector<TH2D*> TH2DVector_muTrigEffSF_barrel2_toterr;
  std::vector<TH2D*> TH2DVector_muTrigEffData_barrel2;
  std::vector<TH2D*> TH2DVector_muTrigEffMC_barrel2;

  std::vector<TH2D*> TH2DVector_muTrigEffSF_endcap2;
  std::vector<TH2D*> TH2DVector_muTrigEffSF_endcap2_staterr;
  std::vector<TH2D*> TH2DVector_muTrigEffSF_endcap2_systerr;
  std::vector<TH2D*> TH2DVector_muTrigEffSF_endcap2_toterr;
  std::vector<TH2D*> TH2DVector_muTrigEffData_endcap2;
  std::vector<TH2D*> TH2DVector_muTrigEffMC_endcap2;

  TFile* ifile_muID;
  TFile* ifile_muTrigEff_barrel1;
  TFile* ifile_muTrigEff_endcap1;
};
#endif
