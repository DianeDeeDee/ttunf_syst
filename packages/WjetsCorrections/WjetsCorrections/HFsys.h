// from electron file
// contact: M.Vreeswijk, h73@nikhef.nl, version 201201605 --> muon_hfxCA_elec
// pre-declaration of helper function (dont call this function directly):

// from muon file
// contact: M.Vreeswijk, h73@nikhef.nl, version 20120905, /project/atlas/users/h73/kfac/furious/log.qcdsym
// pre-declaration of helper function (dont call this function directly):

#ifndef HFSYS_H
#define HFSYS_H

#include "TString.h"

void GetFFactors_muon(int idsys, TString sysname, double Fbb[9],double Fcc[9],double Fc[9], double Fll[9], double cafac[9],
	double fbb[9],double fcc[9],double fc[9], double fll[9],
	double winpretag[9],double wintag[9],double woutpretag[9], double wouttag[9]);


// This is the function to use:
void SetWflavors_muon(int idsys, TString sysname, int ijet, double Wjets_in[4],double Wjets_out[4], double& canorm, int _mode=0);

void GetFFactors_elec(int idsys, TString sysname, double Fbb[9],double Fcc[9],double Fc[9], double Fll[9], double cafac[9],
	double fbb[9],double fcc[9],double fc[9], double fll[9],
	double winpretag[9],double wintag[9],double woutpretag[9], double wouttag[9]);


// This is the function to use:
void SetWflavors_elec(int idsys, TString sysname, int ijet, double Wjets_in[5],double Wjets_out[5], double& canorm, int _mode=0);
  


#endif

