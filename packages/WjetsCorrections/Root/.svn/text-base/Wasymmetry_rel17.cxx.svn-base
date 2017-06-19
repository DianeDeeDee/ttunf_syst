// ------------------------------
// 
// -- W charge asymmetry macro --
// Trieste, 06/10/2011
// updates:
//  * 19/10/2011 : 1fb-1 final version
//  * 24/10/2011 : central values fix and more syst splitting
//  * 28/10/2011 : JES uncertainty updated and +/- split; mcst corrected
//  * 22/11/2011 : bugfix in JES and MC non W backgrounds. To be considered as the main version of the macro for 1 fb-1.
//  * 05/12/2011 : PDF fixed and QCD uncertainty quoted separately
//  * 21/12/2011 : 2 fb-1 numbers
// Udine/ICTP group:
// B. Acharya, U. De Sanctis, M. Pinamonti, K. Shaw, R. Soualah
// contact: michele.pinamonti@gmail.com
// 
// ------------------------------

// ------- //
// CONTENT //
// ------- //

// The macro provides two sets of functions:
// - the functions in the firts part ("INDEPENDENT ESTIMATE")
//   are supposed to be used to re-do the W charge asymmetry estimate 
//   with different data-set or different event selection
//   (currently pre-tag estimates (pure charge asymmetry method) and tag estimates (with the method described below), syst uncertainties not included...)
// - the functions in the second part ("DEFAULT ESTIMATE") 
//   can be used to get the SFs to apply to W+jets MC as calculated using 2fb-1 of data; 
//   this part includes all the syst uncertainties for both pre-tag and b-tag estimate with the Top Group Baseline Selection

#include "WjetsCorrections/Wasymmetry_rel17.h"

// ---------------------------------------------------------------------------------------------------------------------------------

// -------------------- //
// INDEPENDENT ESTIMATE //
// -------------------- //

// ---------------------------------------------------------------------------------------------------------------------------------

// ------------------------------
// function to get the W estimate
// 
// Inputs: PRETAG ESTIMATES
// - DataPlus: number of events in data with a positively charged lepton
// - DataMinus: number of events in data with a negatively charged lepton
// - SumBkgPlus: sum of background events (non-W, non-ttbar MC + QCD) with a positively charged lepton
// - SumBkgMinus: sum of background events (non-W, non-ttbar MC + QCD) with a negatively charged lepton
// - WmcPlus: number of MC W+jets events with a positively charged lepton
// - WmcMinus: number of MC W+jets events with a negatively charged lepton

// Inputs: TAG ESTIMATES

// DETERMINATION of ftag_2j ratio
// - DataPre2j: number of events in data pretag 2j bin
// - DataTag2j: number of events in data tag 2j bin
// - SumBkgPre2j: sum of background events (non-W, non-ttbar MC + QCD) in 2j bin pretag
// - SumBkgTag2j: sum of background events (non-W, non-ttbar MC + QCD) in 2j bin tag
// - ttPre2j: number of ttbar events in 2j bin pretag;
// - ttTag2j: number of ttbar events in 2j bin tag;

// DETERMINATION of ftag_2j-->ftag_(ith)j factor
// - WPre2j: number of events in W+jets MC pretag 2j bin
// - WTag2j: number of events in W+jets MC tag 2j bin
// - WPreNj: number of events in W+jets MC pretag i-th bin
// - WTagNj: number of events in W+jets MC tag i-th bin

// ABSOLUTE ESTIMATES

// estimated number of W+jets events PRETAG(dd)
double GetWestimatePretag( int DataPlus, int DataMinus, double SumBkgPlus, double SumBkgMinus, double WmcPlus, double WmcMinus ){
 double R = WmcPlus / WmcMinus;
 double C = (R+1.)/(R-1.);
 double DataCorrPlus = DataPlus - SumBkgPlus;
 double DataCorrMinus = DataMinus - SumBkgMinus;
 return ( DataCorrPlus - DataCorrMinus )*C;
}

// Computation of the f_2j tag/pretag ratio in DATA in 2j bin after the subtraction of all non W-background and ttbar contributions
double Getf2j( int DataPre2j, int DataTag2j, double SumBkgPre2j, double SumBkgTag2j, double ttPre2j, double ttTag2j ){
   return ( (DataTag2j-SumBkgTag2j-ttTag2j)/(DataPre2j-SumBkgPre2j-ttPre2j));
}

// Computation of the ftag_2j-->ftag_(ith)j factor in W+jets MC 
double Getf2toN( double WPre2j, double WTag2j, double WPreNj, double WTagNj){
   double f4j = WTagNj/WPreNj;
   double f2j = WTag2j/WPre2j;
   return f4j/f2j ;
}

// THE FINAL TAG ESTIMATE IS THE product of the pretag estimate TIMES the f_2j tag ratio TIMES the f_2j to i-th bin extrapolation factor.

// SCALE FACTORS


// dd/MC scale factor for PRETAG
double GetWscaleFactorPretag( int DataPlus, int DataMinus, double SumBkgPlus, double SumBkgMinus, double WmcPlus, double WmcMinus ){
 double R = WmcPlus / WmcMinus;
 double C = (R+1.)/(R-1.);
 double DataCorrPlus = DataPlus - SumBkgPlus;
 double DataCorrMinus = DataMinus - SumBkgMinus;
 double estimate = ( DataCorrPlus - DataCorrMinus )*C;
 return estimate / ( WmcPlus + WmcMinus );
}

// dd/MC scale factor for TAG: choose your favorite bin and put in WTagNj the number of events in W+jets MC in that tagged bin
double GetWscaleFactorTag( int DataPlus, int DataMinus, double SumBkgPlus, double SumBkgMinus, double WmcPlus, double WmcMinus,
int DataPre2j, int DataTag2j, double SumBkgPre2j, double SumBkgTag2j, double ttPre2j, double ttTag2j,
double WPre2j, double WTag2j, double WPreNj, double WTagNj){
   double pretag  = GetWestimatePretag( DataPlus, DataMinus, SumBkgPlus, SumBkgMinus, WmcPlus, WmcMinus) ;
   double f2j     = Getf2j( DataPre2j, DataTag2j, SumBkgPre2j, SumBkgTag2j, ttPre2j, ttTag2j ) ;
   double f2jtoNj = Getf2toN( WPre2j, WTag2j, WPreNj, WTagNj) ;
   double estimate = pretag*f2j*f2jtoNj ;
   return estimate/WTagNj ;

}


// STATISTICAL UNCERTAINTIES

// The statistical error on the PRETAG estimate. 

// relative statistical uncertainty on the W pretag estimate (or on the W dd/MC scale factor)
double GetWstatRelPretag( int DataPlus, int DataMinus, double SumBkgPlus, double SumBkgMinus ){
 double DataCorrPlus = DataPlus - SumBkgPlus;
 double DataCorrMinus = DataMinus - SumBkgMinus;
 return sqrt( DataPlus + DataMinus )/( DataCorrPlus - DataCorrMinus );
}

// The statistical error on the f2j ratio. 


// relative statistical uncertainty
double GetWstatRelf2j( int DataPre2j, int DataTag2j, double SumBkgPre2j, double SumBkgTag2j, double ErrSumBkgPre2j, double ErrSumBkgTag2j, double ttPre2j, double ttTag2j, double ErrttPre2j, double ErrttTag2j ){
 double errnum = sqrt( DataTag2j + ErrSumBkgTag2j*ErrSumBkgTag2j + ErrttTag2j*ErrttTag2j);
 double errden = sqrt( DataPre2j + ErrSumBkgPre2j*ErrSumBkgPre2j + ErrttPre2j*ErrttPre2j);
 double num = DataTag2j-SumBkgTag2j-ttTag2j;
 double den = DataPre2j-SumBkgPre2j-ttPre2j;
 return sqrt((errnum/num)*(errnum/num)+(errden/den)*(errden/den));
}

// The statistical error on the extrapolation factor 2j -> i-th j bin


// relative statistical uncertainty 
double GetWstatRelf2toN( double WPre2j, double WTag2j, double WPreNj, double WTagNj, double ErrWPre2j, double ErrWTag2j, double ErrWPreNj, double ErrWTagNj){
 double errnum = (WTagNj/WPreNj)*sqrt((ErrWTagNj/WTagNj)*(ErrWTagNj/WTagNj) + (ErrWPreNj/WPreNj)*(ErrWPreNj/WPreNj)); 
 double errden = (WTag2j/WPre2j)*sqrt((ErrWTag2j/WTag2j)*(ErrWTag2j/WTag2j) + (ErrWPre2j/WPre2j)*(ErrWPre2j/WPre2j));
 double num = WTagNj/WPreNj;
 double den = WTag2j/WPre2j;
 return sqrt((errnum/num)*(errnum/num)+(errden/den)*(errden/den));

}

// The statistical error on the TAG estimate is from the error propagation of the RELATIVE or ABSOLUTE UNCERTANTIES as you prefer.


// ------------------------------
// Functions to get a weight for each event
// 
// Inputs:
// - charge: the charge of the lepton in the event
// - R: is the ratio WmcPlus / WmcMinus
// - isData: whether the event is a data event
// - isBkg: whether the event is a background (non-W, non-ttbar MC, QCD) event

// weight: apply it to each data and background (non-W, non-ttbar MC, QCD) event, setting the two bool accordingly
double GetWweight( int charge, double R, bool isData, bool isBkg ){
 if ( charge > 0 ) charge = 1;
 if ( charge < 0 ) charge = -1;
 double C = (R+1.)/(R-1.);
 double result = 0.;
 if(isData) result += charge*C;
 if(isBkg) result -= charge*C;
 return result;
}

// apply it to each data event to have the squared sum of the statistical uncertainty: the uncertainty is then the square root of the intergal
double GetWweightSqErr( int charge, double R, bool isData ){
 if ( charge > 0 ) charge = 1;
 if ( charge < 0 ) charge = -1;
 double C = (R+1.)/(R-1.);
 double result = 0.;
 if(isData) result += C*C;
 return result;
}

// ---------------------------------------------------------------------------------------------------------------------------------

// ---------------- //
// SYSTEMATICS FOR THE DEFAULT ESTIMATE //
// ---------------- //

// ---------------------------------------------------------------------------------------------------------------------------------

// Inputs:
//  channel: 0 = e, 1 = mu
//  njet: 4.1 = 4 or more, 4 = exactly 4
//  btag: 0 = pretag, 1 = jfcNN70

// Example of usage:
// // W+jets SF in mu channel, 3j bin, after tagging:
//  cout << GetSFvalue(1,3,1);
// // Statistical uncertainty:
//  cout << " +/-" << GetSFunc("stat",1,3,1) << " (stat)";
// // Systematic uncertainty:
//  cout << " +/-" << GetSFunc("syst",1,3,1) << " (syst)" << endl;
// // One should get:
//  0.813 +/-0.0590828 (stat) +/-0.15925 (syst)

// to specify the type for the split of the syst:
//  "all" : stat + syst
//   "stat" : statistical uncertainty only, including:
//    "stat1" : pretag stat only (uncorrelated between different jet bins)
//    "stat2" : ftag2j stat only (always correlated between the different jet bins, = 0 for pretag estimate)
//   "syst" : systematic uncertainty only; it includes:
//    "misid" : charge mis-id
//    "mcgen" : MC generator (Alpgen vs Sherpa)
//    "pdf" : PDF uncertainty
//    "mcst1" : MC statistic in the pretag estimate (on R)
//    "mcst2" : MC statistic in the f2->N factor for b-tag estimate (= 0 for pretag)
//    "xsec" : Cross-section systematics. Negligible in pretag, reported for ftag2j and then propagates to btag estimate
//    "hf1" : heavy flavour rescaling in W+jets MC (Wbb/cc and Wc symultaneous shift up and down)
//    "hf2" : heavy flavour rescaling in W+jets MC (the 25% for extrapolation to !=2jet bins for Wbb/cc)
//    "hf3" : heavy flavour rescaling in W+jets MC (the 25% for extrapolation to !=2jet bins for Wc)
//    "jes" : JES
//    "bes" : b-JES
//    "btag" : b/c-tag efficiency uncertainty (= 0 for pretag)
//    "ltag" : light-tag efficiency uncertainty (= 0 for pretag)

// ---------------------------------------------------------------------------------------------------------------------------------

// function declaration

double GetRvalue(int channel, float njet, int btag);     // R = W+/W- value from MC
double GetWvalue(int channel, float njet, int btag);    // W(data-driven) value, from 2fb-1 2011 data
double GetSFvalue(int channel, float njet, int btag);    // SF = W(data-driven)/W(MC) value, from 2fb-1 2011 data

// relative uncertainty
double GetSFuncRel(string type, int channel, float njet, int btag);         // uncertainties are summed in quadrature
double GetSFuncRelSingle(string type, int channel, float njet, int btag); // single uncertainty, with sign

// ---------------------------------------------------------------------------------------------------------------------------------

// Tabulated values
// 1j,2j,3j,4+j,4j,5+j

// R=W+/W- from MC
const double R_el_pretag[6]  = { 1.4869, 1.5053, 1.5310, 1.5548, 1.5477, 1.5765 };
const double R_mu_pretag[6]  = { 1.4990, 1.5122, 1.5569, 1.6463, 1.6391, 1.6727 };
const double R_el_jfcNN70[6] = { 1.2608, 1.3067, 1.4726, 1.4897, 1.4512, 1.5935 };
const double R_mu_jfcNN70[6] = { 1.1983, 1.3284, 1.4153, 1.5838, 1.5486, 1.6796 };

// estimaetd W from 2fb-1 2011 data
const double W_el_pretag[6] = {  468521,  126996,  33125,  11476,  8639,  2844   };
const double W_mu_pretag[6] = {  1150948,  286850,  63357,  19282,  15173,  4152   };
const double W_el_jfcNN70[6] = {  12848,  8330,  3567,  1829,  1293,  536  };
const double W_mu_jfcNN70[6] = {  31910,  19188,  6862,  3036,  2225,  802   };
// const double W_el_pretag[6]  = {  468401.1,  127009.7,  33108.2,  11473.7,  8637.3,  2844.5   };// obsolete
// const double W_mu_pretag[6]  = {  1150366.5,  286710.7,  63331.3,  19272.7,  15165.8,  4150   };// obsolete
// const double W_el_jfcNN70[6] = {  12293,  7978.3,  3415.4,  1752.2,  1239.3,  513.5   };// obsolete
// const double W_mu_jfcNN70[6] = {  31951.3,  19226.8,  6880.2,  3043.9,  2231,  804.6   };// obsolete

// SF = W(dd)/W(MC) from 2fb-1 2011 data
const double SF_el_pretag[6] = {  0.912,  0.871,  0.862,  0.899,  0.899,  0.901   };
const double SF_mu_pretag[6] = {  0.980,  0.939,  0.856,  0.827,  0.844,  0.776   };
const double SF_el_jfcNN70[6] = {  0.846,  0.808,  0.800,  0.834,  0.834,  0.836   };
const double SF_mu_jfcNN70[6] = {  0.950,  0.91,  0.83,  0.801,  0.819,  0.752   };
// const double SF_el_pretag[6]  = { 0.912, 0.871, 0.862, 0.899, 0.899, 0.901 };// obsolete
// const double SF_mu_pretag[6]  = { 0.980, 0.939, 0.856, 0.827, 0.844, 0.776 };// obsolete
// const double SF_el_jfcNN70[6] = { 0.810, 0.774, 0.766, 0.799, 0.799, 0.800 }; // obsolete
// const double SF_mu_jfcNN70[6] = { 0.952, 0.913, 0.832, 0.803, 0.821, 0.754 }; // obsolete

// Uncertainties

// statistical uncertainty on pretag estimate - UNCORRELATED btw jet bins and channels
const double SFstat1_el_pretag[6]  = { 0.0081, 0.0161, 0.0333, 0.0671, 0.0738, 0.1492 }; 
const double SFstat1_mu_pretag[6]  = { 0.005, 0.01, 0.0221, 0.0430, 0.0465, 0.1029 };
const double SFstat1_el_jfcNN70[6] = { 0.0214, 0.0227, 0.0389, 0.0709, 0.078, 0.1533 };
const double SFstat1_mu_jfcNN70[6] = { 0.0146, 0.0152, 0.0266, 0.0461, 0.0499, 0.1057 };

// statistical uncertainty on 2j bin tag rate - CORRELATED btw jet bins, UNCORRELATED btw channels
const double SFstat2_el_pretag[6]  = { 0, 0, 0, 0, 0, 0 };
const double SFstat2_mu_pretag[6]  = { 0, 0, 0, 0, 0, 0 };
const double SFstat2_el_jfcNN70[6] = { 0.0126,0.0126,0.0126,0.0126,0.0126,0.0126};
const double SFstat2_mu_jfcNN70[6] = { 0.0092, 0.0092, 0.0092, 0.0092, 0.0092, 0.0092};

// systematic uncertainties on pretag estimate

// charge mis-id
const double SFmisid_el_pretag[6]  = { 0.0182, 0.0199, 0.0187, 0.0226, 0.0222, 0.0232 };
const double SFmisid_mu_pretag[6]  = { 0, 0, 0, 0, 0, 0 };
const double SFmisid_el_jfcNN70[6] = { 0.0182, 0.0199, 0.0187, 0.0226, 0.0222, 0.0232 };
const double SFmisid_mu_jfcNN70[6] = { 0, 0, 0, 0, 0, 0 };

// Alpgen vs Sherpa
const double SFmcgen_el_pretag[6]  = { 0.0826, 0.0637, 0.1340, 0.0517, 0.0396, 0.0894 };
const double SFmcgen_mu_pretag[6]  = { 0.0811, 0.0675, 0.0575, 0.0806, 0.0883, 0.0554 };
const double SFmcgen_el_jfcNN70[6] = { 0.0826, 0.0637, 0.1340, 0.0517, 0.0396, 0.0894 };
const double SFmcgen_mu_jfcNN70[6] = { 0.0008, 0.0675, 0.0575, 0.0806, 0.0883, 0.0554 };

// PDFs
const double SFpdf_el_pretag[6]  = { 0.0550, 0.0574, 0.0557, 0.0553, 0.0550, 0.0556 };
const double SFpdf_mu_pretag[6]  = { 0.0557, 0.0588, 0.0545, 0.0498, 0.0496, 0.0489 };
const double SFpdf_el_jfcNN70[6] = { 0.0550, 0.0574, 0.0557, 0.0553, 0.0550, 0.0556 };
const double SFpdf_mu_jfcNN70[6] = { 0.0557, 0.0588, 0.0545, 0.0498, 0.0496, 0.0489 };

// MC statistics (in pretag R value)
const double SFmcst1_el_pretag[6]  = { 0.0029, 0.0037, 0.0058, 0.0083, 0.0097, 0.0159 };
const double SFmcst1_mu_pretag[6]  = { 0.0018, 0.0024, 0.0044, 0.0069, 0.0079, 0.0135 };
const double SFmcst1_el_jfcNN70[6] = { 0.0029, 0.0037, 0.0058, 0.0083, 0.0097, 0.0159 };
const double SFmcst1_mu_jfcNN70[6] = { 0.0018, 0.0024, 0.0044, 0.0069, 0.0079, 0.0135  };

// systematic uncertainties on both pretag and btag

// bb/cc and c Up and Down - CORRELATED btw jet bins and channels
const double SFhf1_el_pretag[6]  = {  0.0434,  0.0555,  0.0521,  0.0354,  0.0311,  0.0509 };
const double SFhf1_mu_pretag[6]  = {  0.0376,  0.0508,  0.0629,  0.0505,  0.0605,  0.0158 };
const double SFhf1_el_jfcNN70[6] = {  0.0095,  0.0555,  0.0786,  0.0635,  0.0477,  0.1079 };
const double SFhf1_mu_jfcNN70[6] = {  0.0087,  0.0508,  0.0915,  0.0884,  0.0963,  0.0578 };

// bb/cc Up and Down by 25% - UNCORRELATED btw jet bins, CORRELATED btw channels (?)
const double SFhf2_el_pretag[6]  = {  0.0001,  0., -0.0011, -0.0114, -0.0143, -0.0016 };
const double SFhf2_mu_pretag[6]  = { -0.0023,  0.,  0.0028, -0.0040, -0.0011, -0.0146 };
const double SFhf2_el_jfcNN70[6] = {  0.0558,  0.,  0.1208,  0.1226,  0.1125,  0.1508 };
const double SFhf2_mu_jfcNN70[6] = {  0.0467,  0.,  0.1260,  0.1317,  0.1325,  0.1258 };

// c Up and Down by 25% - UNCORRELATED btw jet bins, CORRELATED btw channels (?)
const double SFhf3_el_pretag[6]  = {  0.0346,  0.,  0.0432,  0.0453,  0.0462,  0.0430 };
const double SFhf3_mu_pretag[6]  = {  0.0335,  0.,  0.0461,  0.0464,  0.0500,  0.0344 };
const double SFhf3_el_jfcNN70[6] = {  0.1349,  0.,  0.0915,  0.0767,  0.0796,  0.0690 };
const double SFhf3_mu_jfcNN70[6] = {  0.1420,  0.,  0.0908,  0.0796,  0.0848,  0.0635 };

// JES
const double SFjesUp_el_pretag[6]  = { -0.0038, 0.0106,  0.0088, -0.0219, -0.0750, 0.1450 };
const double SFjesUp_mu_pretag[6]  = { -0.0200, 0.0048, -0.0451,  0.0361,  0.0296, 0.0596 };
const double SFjesUp_el_jfcNN70[6] = {  0.0179, 0.0219,  0.0082, -0.0199, -0.0851, 0.1604 };
const double SFjesUp_mu_jfcNN70[6] = {  0.0078, 0.0110, -0.0383,  0.0291,  0.0119, 0.0479 };

const double SFjesDown_el_pretag[6]  = { -0.0193,  0.0158, -0.0490, -0.0895, -0.0919, -0.0859 };
const double SFjesDown_mu_pretag[6]  = {  0.0184,  0.0025, -0.0081,  0.0001, -0.0097,  0.0227 };
const double SFjesDown_el_jfcNN70[6] = { -0.0171,  0.0065, -0.0298, -0.1070, -0.1101, -0.0831 };
const double SFjesDown_mu_jfcNN70[6] = {  0.0203, -0.0052,  0.0016,  0.0148,  0.0124,  0.0371 };

// bJES
const double SFbesUp_el_pretag[6]  = {  0.0003, -0.0007,  0.0021, -0.0051, -0.0084,  0.0048 };
const double SFbesUp_mu_pretag[6]  = { -0.0001, -0.0004,  0.0001, -0.0002,  0.0003, -0.0018 };
const double SFbesUp_el_jfcNN70[6] = { -0.0056, -0.0092, -0.0091, -0.0156, -0.0194, -0.0047 };
const double SFbesUp_mu_jfcNN70[6] = { -0.0017, -0.0007,  0.0017, -0.0008,  0.0000, -0.0032 };

const double SFbesDown_el_pretag[6]  = {  0.0004, -0.0004,  0.0029, -0.0019, -0.0013, -0.0034 };
const double SFbesDown_mu_pretag[6]  = {  0.0003, -0.0004, -0.0008, -0.0019, -0.0016, -0.0028 };
const double SFbesDown_el_jfcNN70[6] = { -0.0086, -0.0069,  0.0011, -0.0023, -0.0010, -0.0055 };
const double SFbesDown_mu_jfcNN70[6] = {  0.0008,  0.0002,  0.0002,  0.0015,  0.0018,  0.0006 };

// QCD uncertainty systematic
const double SFqcd_el_pretag[6] = {0., 0., 0., 0., 0., 0. };
const double SFqcd_mu_pretag[6] = { 0., 0., 0., 0., 0., 0. };
const double SFqcd_el_jfcNN70[6] = { -0.0423, -0.0423, -0.0423, -0.0423, -0.0423, -0.0423 };
const double SFqcd_mu_jfcNN70[6] = { -0.0464, -0.0464, -0.0464, -0.0464, -0.0464, -0.0464 };

// X-sec systematics
const double SFxsec_el_pretag[6]  = { 0., 0., 0., 0., 0., 0. };
const double SFxsec_mu_pretag[6]  = { 0., 0., 0., 0., 0., 0. };
const double SFxsec_el_jfcNN70[6] = { -0.0271, -0.0271, -0.0271, -0.0271, -0.0271, -0.0271 };
const double SFxsec_mu_jfcNN70[6] = { -0.0168, -0.0168, -0.0168, -0.0168, -0.0168, -0.0168 };

// systematic uncertainties on btag only

// bTag eff
const double SFbtag_el_pretag[6]  = { 0, 0, 0, 0, 0, 0 };
const double SFbtag_mu_pretag[6]  = { 0, 0, 0, 0, 0, 0 };
const double SFbtag_el_jfcNN70[6] = { -0.0337, -0.0380, -0.0371, -0.0463, -0.0488, -0.0401 };
const double SFbtag_mu_jfcNN70[6] = { -0.0207, -0.0264, -0.0299, -0.0347, -0.0326, -0.0401 };

// lightTag eff
const double SFltag_el_pretag[6]  = { 0, 0, 0, 0, 0, 0 };
const double SFltag_mu_pretag[6]  = { 0, 0, 0, 0, 0, 0 };
const double SFltag_el_jfcNN70[6] = { -0.0150, -0.0171, -0.0245, -0.0207, -0.0198, -0.0222 };
const double SFltag_mu_jfcNN70[6] = { -0.0037, -0.0082, -0.0090, -0.0126, -0.0120, -0.0140 };

// statistical uncertainty on 2-> N factor - CORRELATED btw jet bins, UNCORRELATED btw channels
const double SFmcst2_el_pretag[6]  = { 0, 0, 0, 0, 0, 0 };
const double SFmcst2_mu_pretag[6]  = { 0, 0, 0, 0, 0, 0 };
const double SFmcst2_el_jfcNN70[6] = { 0.0152, 0.0099, 0.0157, 0.0194, 0.0221, 0.0326 };
const double SFmcst2_mu_jfcNN70[6] = { 0.0102, 0.0069, 0.0116, 0.0136, 0.0155, 0.0225 };

// ---------------------------------------------------------------------------------------------------------------------------------

// function implementation

// R = W+/W- value from MC
double GetRvalue(int channel, float njet, int btag){
 int nj = -1;
 if(njet<4) nj = (int)njet-1;
 if(((int)njet)==4 && njet>4) nj = 3;
 if(njet == 4) nj = 4;
 if(njet >= 5) nj = 5;
 if(channel == 0 && btag == 0) return R_el_pretag[nj];
 if(channel == 1 && btag == 0) return R_mu_pretag[nj];
 if(channel == 0 && btag == 1) return R_el_jfcNN70[nj];
 if(channel == 1 && btag == 1) return R_mu_jfcNN70[nj];
 return 0.;
}

// W(data-driven) value, from 2fb-1 2011 data
double GetWvalue(int channel, float njet, int btag){
 int nj = -1;
 if(njet<4) nj = (int)njet-1;
 if(((int)njet)==4 && njet>4) nj = 3;
 if(njet == 4) nj = 4;
 if(njet >= 5) nj = 5;
 if(channel == 0 && btag == 0) return W_el_pretag[nj];
 if(channel == 1 && btag == 0) return W_mu_pretag[nj];
 if(channel == 0 && btag == 1) return W_el_jfcNN70[nj];
 if(channel == 1 && btag == 1) return W_mu_jfcNN70[nj];
 return 0.;
}

// SF = W(data-driven)/W(MC) value, from 2fb-1 2011 data
double GetSFvalue(int channel, float njet, int btag){
 int nj = -1;
 if(njet<4) nj = (int)njet-1;
 if(((int)njet)==4 && njet>4) nj = 3;
 if(njet == 4) nj = 4;
 if(njet >= 5) nj = 5;
 if(channel == 0 && btag == 0) return SF_el_pretag[nj];
 if(channel == 1 && btag == 0) return SF_mu_pretag[nj];
 if(channel == 0 && btag == 1) return SF_el_jfcNN70[nj];
 if(channel == 1 && btag == 1) return SF_mu_jfcNN70[nj];
 return 0.;
}

// relative uncertainty
//
// uncertainties are summed in quadrature
double GetSFuncRel(string type, int channel, float njet, int btag){
 if(type=="stat") type = "stat1,stat2";
 if(type=="syst") type = "misid,mcgen,pdf,mcst1,mcst2,xsec,qcd,hf1,hf2,hf3,btag,ltag,jes,bes";
 if(type=="all" || type=="tot") type = "stat1,stat2,misid,mcgen,pdf,mcst1,mcst2,xsec,qcd,hf1,hf2,hf3,btag,ltag,jes,bes";
 if(type=="hf") type = "hf1,hf2,hf3";
 if(type=="mcst") type = "mcst1,mcst2";
 if(type=="other") type = "mcst1,mcst2,misid,mcgen,xsec,qcd";
 double SFunc = 0.;
 if(type.find("stat1")!=string::npos) SFunc = sqrt(SFunc*SFunc + pow( GetSFuncRelSingle("stat1",channel,njet,btag),2 ));
 if(type.find("stat2")!=string::npos) SFunc = sqrt(SFunc*SFunc + pow( GetSFuncRelSingle("stat2",channel,njet,btag),2 ));
 if(type.find("misid")!=string::npos) SFunc = sqrt(SFunc*SFunc + pow( GetSFuncRelSingle("misid",channel,njet,btag),2 ));
 if(type.find("mcgen")!=string::npos) SFunc = sqrt(SFunc*SFunc + pow( GetSFuncRelSingle("mcgen",channel,njet,btag),2 ));
 if(type.find("pdf")!=string::npos) SFunc = sqrt(SFunc*SFunc + pow( GetSFuncRelSingle("pdf",channel,njet,btag),2 ));
 if(type.find("mcst1")!=string::npos) SFunc = sqrt(SFunc*SFunc + pow( GetSFuncRelSingle("mcst1",channel,njet,btag),2 ));
 if(type.find("mcst2")!=string::npos) SFunc = sqrt(SFunc*SFunc + pow( GetSFuncRelSingle("mcst2",channel,njet,btag),2 ));
 if(type.find("xsec")!=string::npos) SFunc = sqrt(SFunc*SFunc + pow( GetSFuncRelSingle("xsec",channel,njet,btag),2 ));
 if(type.find("qcd")!=string::npos) SFunc = sqrt(SFunc*SFunc + pow( GetSFuncRelSingle("qcd",channel,njet,btag),2 ));
 if(type.find("hf1")!=string::npos) SFunc = sqrt(SFunc*SFunc + pow( GetSFuncRelSingle("hf1",channel,njet,btag),2 ));
 if(type.find("hf2")!=string::npos) SFunc = sqrt(SFunc*SFunc + pow( GetSFuncRelSingle("hf2",channel,njet,btag),2 ));
 if(type.find("hf3")!=string::npos) SFunc = sqrt(SFunc*SFunc + pow( GetSFuncRelSingle("hf3",channel,njet,btag),2 ));
 if(type.find("jes")!=string::npos) SFunc = sqrt(SFunc*SFunc + pow( GetSFuncRelSingle("jes",channel,njet,btag),2 ));
 if(type.find("bes")!=string::npos) SFunc = sqrt(SFunc*SFunc + pow( GetSFuncRelSingle("bes",channel,njet,btag),2 ));
 if(type.find("btag")!=string::npos) SFunc = sqrt(SFunc*SFunc + pow( GetSFuncRelSingle("btag",channel,njet,btag),2 ));
 if(type.find("ltag")!=string::npos) SFunc = sqrt(SFunc*SFunc + pow( GetSFuncRelSingle("ltag",channel,njet,btag),2 ));
 return SFunc;
}
//
// single uncertainty, with sign
// - this is for "up" uncertainties; for "down" use (-1)*this
// - for JES and bJES only, one can specify "jesup" or "jesdown" to have to different effect; "jes" alone is giving a positive number only
double GetSFuncRelSingle(string type, int channel, float njet, int btag){
 int nj;
 if(njet<4) nj = (int)njet-1;
 if(((int)njet)==4 && njet>4) nj = 3;
 if(njet == 4) nj = 4;
 if(njet >= 5) nj = 5;
 if(type=="stat1"){
  if(channel == 0 && btag == 0) return SFstat1_el_pretag[nj];
  if(channel == 1 && btag == 0) return SFstat1_mu_pretag[nj];
  if(channel == 0 && btag == 1) return SFstat1_el_jfcNN70[nj];
  if(channel == 1 && btag == 1) return SFstat1_mu_jfcNN70[nj];
 }
 if(type=="stat2"){
  if(channel == 0 && btag == 0) return SFstat2_el_pretag[nj];
  if(channel == 1 && btag == 0) return SFstat2_mu_pretag[nj];
  if(channel == 0 && btag == 1) return SFstat2_el_jfcNN70[nj];
  if(channel == 1 && btag == 1) return SFstat2_mu_jfcNN70[nj];
 }
 if(type=="misid"){
  if(channel == 0 && btag == 0) return SFmisid_el_pretag[nj];
  if(channel == 1 && btag == 0) return SFmisid_mu_pretag[nj];
  if(channel == 0 && btag == 1) return SFmisid_el_jfcNN70[nj];
  if(channel == 1 && btag == 1) return SFmisid_mu_jfcNN70[nj];
 }
 if(type=="mcgen"){
  if(channel == 0 && btag == 0) return SFmcgen_el_pretag[nj];
  if(channel == 1 && btag == 0) return SFmcgen_mu_pretag[nj];
  if(channel == 0 && btag == 1) return SFmcgen_el_jfcNN70[nj];
  if(channel == 1 && btag == 1) return SFmcgen_mu_jfcNN70[nj];
 }
 if(type=="pdf"){
  if(channel == 0 && btag == 0) return SFpdf_el_pretag[nj];
  if(channel == 1 && btag == 0) return SFpdf_mu_pretag[nj];
  if(channel == 0 && btag == 1) return SFpdf_el_jfcNN70[nj];
  if(channel == 1 && btag == 1) return SFpdf_mu_jfcNN70[nj];
 }
 if(type=="mcst1"){
  if(channel == 0 && btag == 0) return SFmcst1_el_pretag[nj];
  if(channel == 1 && btag == 0) return SFmcst1_mu_pretag[nj];
  if(channel == 0 && btag == 1) return SFmcst1_el_jfcNN70[nj];
  if(channel == 1 && btag == 1) return SFmcst1_mu_jfcNN70[nj];
 }
 if(type=="mcst2"){
  if(channel == 0 && btag == 0) return SFmcst2_el_pretag[nj];
  if(channel == 1 && btag == 0) return SFmcst2_mu_pretag[nj];
  if(channel == 0 && btag == 1) return SFmcst2_el_jfcNN70[nj];
  if(channel == 1 && btag == 1) return SFmcst2_mu_jfcNN70[nj];
 }
 if(type=="xsec"){
  if(channel == 0 && btag == 0) return SFxsec_el_pretag[nj];
  if(channel == 1 && btag == 0) return SFxsec_mu_pretag[nj];
  if(channel == 0 && btag == 1) return SFxsec_el_jfcNN70[nj];
  if(channel == 1 && btag == 1) return SFxsec_mu_jfcNN70[nj];
 }
 if(type=="qcd"){
  if(channel == 0 && btag == 0) return SFqcd_el_pretag[nj];
  if(channel == 1 && btag == 0) return SFqcd_mu_pretag[nj];
  if(channel == 0 && btag == 1) return SFqcd_el_jfcNN70[nj];
  if(channel == 1 && btag == 1) return SFqcd_mu_jfcNN70[nj];
 }
 if(type=="hf1"){
  if(channel == 0 && btag == 0) return SFhf1_el_pretag[nj];
  if(channel == 1 && btag == 0) return SFhf1_mu_pretag[nj];
  if(channel == 0 && btag == 1) return SFhf1_el_jfcNN70[nj];
  if(channel == 1 && btag == 1) return SFhf1_mu_jfcNN70[nj];
 }
 if(type=="hf2"){
  if(channel == 0 && btag == 0) return SFhf2_el_pretag[nj];
  if(channel == 1 && btag == 0) return SFhf2_mu_pretag[nj];
  if(channel == 0 && btag == 1) return SFhf2_el_jfcNN70[nj];
  if(channel == 1 && btag == 1) return SFhf2_mu_jfcNN70[nj];
 }
 if(type=="hf3"){
  if(channel == 0 && btag == 0) return SFhf3_el_pretag[nj];
  if(channel == 1 && btag == 0) return SFhf3_mu_pretag[nj];
  if(channel == 0 && btag == 1) return SFhf3_el_jfcNN70[nj];
  if(channel == 1 && btag == 1) return SFhf3_mu_jfcNN70[nj];
 }
 if(type=="jesup"){
  if(channel == 0 && btag == 0) return SFjesUp_el_pretag[nj];
  if(channel == 1 && btag == 0) return SFjesUp_mu_pretag[nj];
  if(channel == 0 && btag == 1) return SFjesUp_el_jfcNN70[nj];
  if(channel == 1 && btag == 1) return SFjesUp_mu_jfcNN70[nj];
 }
 if(type=="jesdown"){
  if(channel == 0 && btag == 0) return SFjesDown_el_pretag[nj];
  if(channel == 1 && btag == 0) return SFjesDown_mu_pretag[nj];
  if(channel == 0 && btag == 1) return SFjesDown_el_jfcNN70[nj];
  if(channel == 1 && btag == 1) return SFjesDown_mu_jfcNN70[nj];
 }
 if(type=="jes"){
  if(channel == 0 && btag == 0) return (fabs(SFjesUp_el_pretag[nj])+fabs(SFjesDown_el_pretag[nj]))/2.;
  if(channel == 1 && btag == 0) return (fabs(SFjesUp_mu_pretag[nj])+fabs(SFjesDown_mu_pretag[nj]))/2.;
  if(channel == 0 && btag == 1) return (fabs(SFjesUp_el_jfcNN70[nj])+fabs(SFjesDown_el_jfcNN70[nj]))/2.;
  if(channel == 1 && btag == 1) return (fabs(SFjesUp_mu_jfcNN70[nj])+fabs(SFjesDown_mu_jfcNN70[nj]))/2.;
 }
 if(type=="besup"){
  if(channel == 0 && btag == 0) return SFbesUp_el_pretag[nj];
  if(channel == 1 && btag == 0) return SFbesUp_mu_pretag[nj];
  if(channel == 0 && btag == 1) return SFbesUp_el_jfcNN70[nj];
  if(channel == 1 && btag == 1) return SFbesUp_mu_jfcNN70[nj];
 }
 if(type=="besdown"){
  if(channel == 0 && btag == 0) return SFbesDown_el_pretag[nj];
  if(channel == 1 && btag == 0) return SFbesDown_mu_pretag[nj];
  if(channel == 0 && btag == 1) return SFbesDown_el_jfcNN70[nj];
  if(channel == 1 && btag == 1) return SFbesDown_mu_jfcNN70[nj];
 }
 if(type=="bes"){
  if(channel == 0 && btag == 0) return (fabs(SFbesUp_el_pretag[nj])+fabs(SFbesDown_el_pretag[nj]))/2.;
  if(channel == 1 && btag == 0) return (fabs(SFbesUp_mu_pretag[nj])+fabs(SFbesDown_mu_pretag[nj]))/2.;
  if(channel == 0 && btag == 1) return (fabs(SFbesUp_el_jfcNN70[nj])+fabs(SFbesDown_el_jfcNN70[nj]))/2.;
  if(channel == 1 && btag == 1) return (fabs(SFbesUp_mu_jfcNN70[nj])+fabs(SFbesDown_mu_jfcNN70[nj]))/2.;
 }
 if(type=="btag"){
  if(channel == 0 && btag == 0) return SFbtag_el_pretag[nj];
  if(channel == 1 && btag == 0) return SFbtag_mu_pretag[nj];
  if(channel == 0 && btag == 1) return SFbtag_el_jfcNN70[nj];
  if(channel == 1 && btag == 1) return SFbtag_mu_jfcNN70[nj];
 }
 if(type=="ltag"){
  if(channel == 0 && btag == 0) return SFltag_el_pretag[nj];
  if(channel == 1 && btag == 0) return SFltag_mu_pretag[nj];
  if(channel == 0 && btag == 1) return SFltag_el_jfcNN70[nj];
  if(channel == 1 && btag == 1) return SFltag_mu_jfcNN70[nj];
 }
 return 0.;
}

