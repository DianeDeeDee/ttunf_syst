/**
 * ########################################## WARNING #################################################################
 * #                                                                                                                  #
 * #  This script is depcrecated and currently only kept for histogrical reasons, HFsys for boosted are implementd    #
 * #  in the HFsys.cxx/.h files aswell, please use those for boosted Wjet SFs !! (deprecate warning added 2014/11/04) #
 * #                                                                                                                  #
 * ########################################## WARNING #################################################################
 */

#ifndef HFSYSBOOSTED_H
#define HFSYSBOOSTED_H

#include "TString.h"

// FROM ELECTRON FILE ON TWIKI https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/WplusJetsBackgroundsforTopAnalyses#Scale_factors_for_final_7_TeV_20 [iwatson]

// contact: M.Vreeswijk, h73@nikhef.nl, version 12dec -->  NEW qcd 6GeV iso elec eta/Pt , >4jet tagged  no iso, loose,  DR/Pt?
// updated with now HF extrapolation uncertainties and including PDF uncertenaties
// pre-declaration of helper function (dont call this function directly):
void GetFFactors_boosted_elec(int idsys, TString sysname, double Fbb[9],double Fcc[9],double Fc[9], double Fll[9], double cafac[9],
		      double fbb[9],double fcc[9],double fc[9], double fll[9],
		      double winpretag[9],double wintag[9],double woutpretag[9], double wouttag[9]);
void GetFFactors8TeV_boosted_elec(int idsys, TString sysname, double Fbb[9],double Fcc[9],double Fc[9], double Fll[9], double cafac[9],
		      double fbb[9],double fcc[9],double fc[9], double fll[9],
		      double winpretag[9],double wintag[9],double woutpretag[9], double wouttag[9]);

// This is the function to use:
void SetWflavors_boosted_elec(int idsys, TString sysname, int ijet, double Wjets_in[5],double Wjets_out[5], double& canorm, int _mode=0);

// FROM MUON FILE ON TWIKI [iwatson]

// contact: M.Vreeswijk, h73@nikhef.nl, version 12dec2012-->  qcd6_update (updating old numbers)
// updated with now HF extrapolation uncertainties and including PDF uncertenaties
// pre-declaration of helper function (dont call this function directly):
void GetFFactors_boosted_muon(int idsys, TString sysname, double Fbb[9],double Fcc[9],double Fc[9], double Fll[9], double cafac[9],
		      double fbb[9],double fcc[9],double fc[9], double fll[9],
		      double winpretag[9],double wintag[9],double woutpretag[9], double wouttag[9]);
void GetFFactors8TeV_boosted_muon(int idsys, TString sysname, double Fbb[9],double Fcc[9],double Fc[9], double Fll[9], double cafac[9],
		      double fbb[9],double fcc[9],double fc[9], double fll[9],
		      double winpretag[9],double wintag[9],double woutpretag[9], double wouttag[9]);

// This is the function to use:
void SetWflavors_boosted_muon(int idsys, TString sysname, int ijet, double Wjets_in[5],double Wjets_out[5], double& canorm, int _mode=0);

#endif

