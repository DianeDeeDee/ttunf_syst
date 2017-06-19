#include "WjetsCorrections/HFsys_factor_ttbar_emu.h"
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

void SetWflavors(int idsys, TString sysname, int ijet, double Wjets_in[4],double Wjets_out[4]) {
  //-
  // INPUTS:
  // for idsys there are two input modes: the individul systematics (long list), or the combined systematics (4 checks).
  // Long list:
  // input idsys>0..999 or sysname are passed to the function GetFFactors (see below) that returns the relevant HF factors.
  // This allows to use the full list of systematics.
  // Short list:
  // most analysers will use the total sum of all systematics. In that case only four checks up and two down are needed.
  // Complete set: idsys=1000 (up), -1000 (down) ---> checks Fbb and anticorrelated Fc. Fl used to renormalize.
  //                     1001 (up), -1001 (down) ---> checks HF versus LF, thus change Fl and use (Fbb+Fc) to renormalize
  //                     1002 (up), -1002 (down) ---> check Fbb+25% in 3jet bin and Fbb+sqrt(2)*25% in 4jet bin (and renormalize all).
  //                     1003 (up), -1003 (down) ---> check Fc+25% in 3jet bin and Fc+sqrt(2)*25% in 4jet bin (and renormalize all).
  //
  // input ijet = jetbin you want to set.
  // input Wjets_in are the PRETAGGED event numbers (for ijet) for 0:Wbb, 1:Wcc, 2:Wc, 3:Wlight
  // The flavor fractions defined as bb, cc, c and light flavor jets correspond to the HFOR flag equals to 0, 1, 2 and 3 respectively
  //--
  // OUTPUTS:
  // output Wjets_out are the event numbers (for ijet) for 0:Wbb, 1:Wcc, 2:Wc, 3:Wlight 
  //
  //
  // Typical use for the full systematics list:
  //      .....
  //      int idsys=0; // for nominal, but you have to loop over many idsys, see in GetFFactors + explanation on www.nikhef.nl/~h73/factors.html
  //      int ijet=3; // in case you input the 3jetbin
  //      double Wjets_in[4]={ 6888, 14281, 22937, 102347} ; 
  //      double Wjets_out[4]; 
  //      SetWflavors(idsys, "" , ijet, Wjets_in, Wjets_out);
  //      .....
  // Typical use when the combined uncertainties (short list)  are used:
  //      .....
  //      int ijet=3; // in case you input the 3jetbin
  //      int idsys=1000; // 
  //      double Wjets_in[4]={ 6888, 14281, 22937, 102347} ;
  //      double Wjets_out[4];
  //      SetWflavors(idsys, "" , ijet, Wjets_in, Wjets_out);
  //      do this also for idsys=1001, 1002, 1003 and all negative: -1000, -1001, -1002, -1003. 
  //      all positive effects in quadrature, and all negative as usual
  //      (use nominal HF factors) as of your systematic study. Same for Wc.
  //      .....
  
  // get HF Kfactors
  double Fbb[5], Fcc[5], Fc[5],  Fll[5];
  if (idsys*idsys<999*999) { // go for the full list of systematics
    GetFFactors( idsys, sysname,  Fbb,Fcc, Fc,  Fll);
  } else {                  // go for the reduced errors.
    GetFFactors( 0, "",  Fbb,Fcc, Fc,  Fll); // get nominal
    double sigma_bb[2]={ 0.85 ,  - 0.88};
    double sigma_c[2] ={ 0.49 ,  - 0.32};
    double sigma_ll[2] ={0.073,  - 0.087};
    
    if (idsys==1000) {
      Fbb[ijet]*=1+sigma_bb[0];
      Fcc[ijet]*=1+sigma_bb[0];
      Fc[ijet]*=1+sigma_c[1];
    } else if (idsys==-1000) {
      Fbb[ijet]*=1+sigma_bb[1];
      Fcc[ijet]*=1+sigma_bb[1];
      Fc[ijet]*=1+sigma_c[0];
    }
    if (idsys==1001) {
      Fll[ijet]*=1+sigma_ll[0];
    } else if (idsys==-1001) {
      Fll[ijet]*=1+sigma_ll[1];
    }
    if (idsys==1002) { // change bb in 3jet and 4jet bin, later all get renormalized.
      Fbb[3]*=1+0.25;
      Fcc[3]*=1+0.25;
      Fbb[4]*=1+0.25*sqrt(2);
      Fcc[4]*=1+0.25*sqrt(2);
    } else if (idsys==-1002) {
      Fbb[3]*=1-0.25;
      Fcc[3]*=1-0.25;
      Fbb[4]*=1-0.25*sqrt(2);
      Fcc[4]*=1-0.25*sqrt(2);
    }
    if (idsys==1003) { // change c in 3jet and 4jet bin, later all get renormalized.
      Fc[3]*=1+0.25;
      Fc[4]*=1+0.25*sqrt(2);
    } else if (idsys==-1003) {
      Fc[3]*=1-0.25;
      Fc[4]*=1-0.25*sqrt(2);
    }
    if (abs(idsys)>1003) cout << " FATAL: unknown systematic check for HF factors " ;
  }
  // calculate new Wjets_out
  Wjets_out[0]=Fbb[ijet]*Wjets_in[0];
  Wjets_out[1]=Fcc[ijet]*Wjets_in[1];
  Wjets_out[2]=Fc[ijet] *Wjets_in[2];
  Wjets_out[3]=Fll[ijet]*Wjets_in[3];
  //  renormalize  to ensure that the overal normalization remains constant
  double nw_in=0;
  double nw_raw=0;
  for (int iw=0;iw<4;iw++) {
    nw_in+=Wjets_in[iw];
    nw_raw+=Wjets_out[iw];
  }
  if (nw_raw==0) { nw_in=0; nw_raw=1;} // some protection
  for (int iw=0;iw<4;iw++) {
    Wjets_out[iw]*=nw_in/nw_raw;
  }
}




void GetFFactors(int idsys, TString sysname, double Fbb[5],double Fcc[5],double Fc[5], double Fll[5]) {
  //--
  // Ffactors for r17. ttbar cuts 28/1/2012.  ELEC+MUON
  // the Ffactors Fbb, Fcc, Fc, Fll contain weight factors for Wbb events, Wcc events, Wc event and Wlight events. 
  // These factors do NOT contains the normalization of MC Wjets events to data.
  //--
  // input: int idsys --> the id of the systematic variation, then call with sysname="" 
  //        TString sysname --> the name of the systematic variation, then call with idsys=-1 (NOT: idsys=0, because that is nominal).
  // output: Fbb[5],Fcc[5],Fc[5],Fll[5]
  //         position Fxx[0] is not used, position 1,2,3,4 are for the jetbins 1,2,3,4
  // typical usage case:
  // You call this function once for nominal and then again for each systematic variation under study.
  // The input is either the integer idsys or the TString sysname. For the first mode with idsys, you typically may use:
  //      .....
  //      int idsys; 
  //      double Fbb[5], Fcc[5], Fc[5],  Fll[5];
  //      idsys=0; // for nominal
  //      GetFFactors(idsys, "" , Fbb,Fcc,Fc,Fll); 
  //      .....
  // IF you prefer the mode with the TString sysname, you typically may use:
  //      .....
  //      TString sysname;  
  //      double Fbb[5], Fcc[5], Fc[5],  Fll[5];
  //      sysname="nominal";
  //      GetFFactors(0, sysname , Fbb,Fcc,Fc,Fll); 
  //      .....
  //
  
  Fbb[1]= 1.11542;   Fcc[1]= 1.11542;   Fc[1]= 1.10602;   Fll[1]= 0.976007;  
  Fbb[2]= 1.10391;   Fcc[2]= 1.10391;   Fc[2]= 1.09461;   Fll[2]= 0.965936;  
  Fbb[3]= 1.09729;   Fcc[3]= 1.09729;   Fc[3]= 1.08804;   Fll[3]= 0.960144;  
  Fbb[4]= 1.09115;   Fcc[4]= 1.09115;   Fc[4]= 1.08195;   Fll[4]= 0.954767;  
  if ( idsys==0   || sysname=="norminal" || sysname=="nominal" ) { 
    // OK
  } 
  else if ( idsys==1   || sysname=="tt_up" ) { 
    Fbb[1]= 0.93297;   Fcc[1]= 0.93297;   Fc[1]= 1.15281;   Fll[1]= 0.978473;  
    Fbb[2]= 0.931613;   Fcc[2]= 0.931613;   Fc[2]= 1.15114;   Fll[2]= 0.97705;  
    Fbb[3]= 0.933717;   Fcc[3]= 0.933717;   Fc[3]= 1.15373;   Fll[3]= 0.979256;  
    Fbb[4]= 0.937467;   Fcc[4]= 0.937467;   Fc[4]= 1.15837;   Fll[4]= 0.98319;  
  } 
  else if ( idsys==2   || sysname=="tt_down" ) { 
    Fbb[1]= 1.30008;   Fcc[1]= 1.30008;   Fc[1]= 1.05866;   Fll[1]= 0.973511;  
    Fbb[2]= 1.27521;   Fcc[2]= 1.27521;   Fc[2]= 1.03841;   Fll[2]= 0.954886;  
    Fbb[3]= 1.25727;   Fcc[3]= 1.25727;   Fc[3]= 1.0238;   Fll[3]= 0.941452;  
    Fbb[4]= 1.23862;   Fcc[4]= 1.23862;   Fc[4]= 1.00862;   Fll[4]= 0.92749;  
  } 
  else if ( idsys==27   || sysname=="WbbWcc_up" ) { 
    Fbb[1]= 0.583259;   Fcc[1]= 0.583259;   Fc[1]= 1.15667;   Fll[1]= 1.02071;  
    Fbb[2]= 0.606981;   Fcc[2]= 0.606981;   Fc[2]= 1.20371;   Fll[2]= 1.06222;  
    Fbb[3]= 0.627959;   Fcc[3]= 0.627959;   Fc[3]= 1.24531;   Fll[3]= 1.09893;  
    Fbb[4]= 0.651843;   Fcc[4]= 0.651843;   Fc[4]= 1.29268;   Fll[4]= 1.14073;  
  } 
  else if ( idsys==28   || sysname=="WbbWcc_down" ) { 
    Fbb[1]= 2.17979;   Fcc[1]= 2.17979;   Fc[1]= 1.08069;   Fll[1]= 0.953658;  
    Fbb[2]= 2.09781;   Fcc[2]= 2.09781;   Fc[2]= 1.04004;   Fll[2]= 0.917794;  
    Fbb[3]= 2.03599;   Fcc[3]= 2.03599;   Fc[3]= 1.0094;   Fll[3]= 0.890749;  
    Fbb[4]= 1.96979;   Fcc[4]= 1.96979;   Fc[4]= 0.976574;   Fll[4]= 0.861783;  
  } 
  else if ( idsys==19   || sysname=="Kbb_up" ) { 
    Fbb[1]= 1.15128;   Fcc[1]= 1.15128;   Fc[1]= 1.1041;   Fll[1]= 0.974318;  
    Fbb[2]= 1.13713;   Fcc[2]= 1.13713;   Fc[2]= 1.09053;   Fll[2]= 0.962343;  
    Fbb[3]= 1.12846;   Fcc[3]= 1.12846;   Fc[3]= 1.08222;   Fll[3]= 0.955004;  
    Fbb[4]= 1.1201;   Fcc[4]= 1.1201;   Fc[4]= 1.0742;   Fll[4]= 0.94793;  
  } 
  else if ( idsys==20   || sysname=="Kbb_down" ) { 
    Fbb[1]= 1.07944;   Fcc[1]= 1.07944;   Fc[1]= 1.10794;   Fll[1]= 0.977702;  
    Fbb[2]= 1.07045;   Fcc[2]= 1.07045;   Fc[2]= 1.09871;   Fll[2]= 0.969557;  
    Fbb[3]= 1.06579;   Fcc[3]= 1.06579;   Fc[3]= 1.09393;   Fll[3]= 0.96534;  
    Fbb[4]= 1.06178;   Fcc[4]= 1.06178;   Fc[4]= 1.08981;   Fll[4]= 0.961703;  
  } 
  else if ( idsys==21   || sysname=="Kc_up" ) { 
    Fbb[1]= 1.11278;   Fcc[1]= 1.11278;   Fc[1]= 1.12089;   Fll[1]= 0.973695;  
    Fbb[2]= 1.1009;   Fcc[2]= 1.1009;   Fc[2]= 1.10892;   Fll[2]= 0.963297;  
    Fbb[3]= 1.09434;   Fcc[3]= 1.09434;   Fc[3]= 1.10232;   Fll[3]= 0.957557;  
    Fbb[4]= 1.08841;   Fcc[4]= 1.08841;   Fc[4]= 1.09634;   Fll[4]= 0.952368;  
  } 
  else if ( idsys==22   || sysname=="Kc_down" ) { 
    Fbb[1]= 1.11808;   Fcc[1]= 1.11808;   Fc[1]= 1.09107;   Fll[1]= 0.978331;  
    Fbb[2]= 1.10695;   Fcc[2]= 1.10695;   Fc[2]= 1.08021;   Fll[2]= 0.96859;  
    Fbb[3]= 1.10027;   Fcc[3]= 1.10027;   Fc[3]= 1.07369;   Fll[3]= 0.962746;  
    Fbb[4]= 1.0939;   Fcc[4]= 1.0939;   Fc[4]= 1.06748;   Fll[4]= 0.957177;  
  } 
  else if ( idsys==23   || sysname=="Kll_up" ) { 
    Fbb[1]= 1.10967;   Fcc[1]= 1.10967;   Fc[1]= 1.10031;   Fll[1]= 0.977274;  
    Fbb[2]= 1.0988;   Fcc[2]= 1.0988;   Fc[2]= 1.08953;   Fll[2]= 0.967702;  
    Fbb[3]= 1.09253;   Fcc[3]= 1.09253;   Fc[3]= 1.08332;   Fll[3]= 0.962186;  
    Fbb[4]= 1.08671;   Fcc[4]= 1.08671;   Fc[4]= 1.07755;   Fll[4]= 0.957059;  
  } 
  else if ( idsys==24   || sysname=="Kll_down" ) { 
    Fbb[1]= 1.12124;   Fcc[1]= 1.12124;   Fc[1]= 1.11179;   Fll[1]= 0.974728;  
    Fbb[2]= 1.10908;   Fcc[2]= 1.10908;   Fc[2]= 1.09973;   Fll[2]= 0.964155;  
    Fbb[3]= 1.10209;   Fcc[3]= 1.10209;   Fc[3]= 1.0928;   Fll[3]= 0.958084;  
    Fbb[4]= 1.09562;   Fcc[4]= 1.09562;   Fc[4]= 1.08638;   Fll[4]= 0.952456;  
  } 
  else if ( idsys==31   || sysname=="WbbWccjet_up" ) { 
    Fbb[1]= 1.02319;   Fcc[1]= 1.02319;   Fc[1]= 1.18268;   Fll[1]= 0.96881;  
    Fbb[2]= 1.01508;   Fcc[2]= 1.01508;   Fc[2]= 1.17331;   Fll[2]= 0.961131;  
    Fbb[3]= 1.01244;   Fcc[3]= 1.01244;   Fc[3]= 1.17025;   Fll[3]= 0.958628;  
    Fbb[4]= 1.01147;   Fcc[4]= 1.01147;   Fc[4]= 1.16913;   Fll[4]= 0.957706;  
  } 
  else if ( idsys==32   || sysname=="WbbWccjet_down" ) { 
    Fbb[1]= 1.219;   Fcc[1]= 1.219;   Fc[1]= 1.01359;   Fll[1]= 0.98244;  
    Fbb[2]= 1.20581;   Fcc[2]= 1.20581;   Fc[2]= 1.00262;   Fll[2]= 0.971811;  
    Fbb[3]= 1.20204;   Fcc[3]= 1.20204;   Fc[3]= 0.999489;   Fll[3]= 0.968772;  
    Fbb[4]= 1.19582;   Fcc[4]= 1.19582;   Fc[4]= 0.994318;   Fll[4]= 0.963759;  
  } 
  else if ( idsys==39   || sysname=="Wcjet_up" ) { 
    Fbb[1]= 0.565277;   Fcc[1]= 0.565277;   Fc[1]= 1.59919;   Fll[1]= 0.949983;  
    Fbb[2]= 0.557411;   Fcc[2]= 0.557411;   Fc[2]= 1.57694;   Fll[2]= 0.936765;  
    Fbb[3]= 0.554859;   Fcc[3]= 0.554859;   Fc[3]= 1.56971;   Fll[3]= 0.932475;  
    Fbb[4]= 0.564547;   Fcc[4]= 0.564547;   Fc[4]= 1.59712;   Fll[4]= 0.948757;  
  } 
  else if ( idsys==40   || sysname=="Wcjet_down" ) { 
    Fbb[1]= 1.44369;   Fcc[1]= 1.44369;   Fc[1]= 0.853686;   Fll[1]= 1.00543;  
    Fbb[2]= 1.40824;   Fcc[2]= 1.40824;   Fc[2]= 0.832729;   Fll[2]= 0.980753;  
    Fbb[3]= 1.37126;   Fcc[3]= 1.37126;   Fc[3]= 0.810861;   Fll[3]= 0.954997;  
    Fbb[4]= 1.33645;   Fcc[4]= 1.33645;   Fc[4]= 0.790273;   Fll[4]= 0.93075;  
  } 
  else if ( idsys==3   || sysname=="wt_up" ) { 
    Fbb[1]= 0.919397;   Fcc[1]= 0.919397;   Fc[1]= 1.14606;   Fll[1]= 0.980348;  
    Fbb[2]= 0.919002;   Fcc[2]= 0.919002;   Fc[2]= 1.14557;   Fll[2]= 0.979926;  
    Fbb[3]= 0.921709;   Fcc[3]= 0.921709;   Fc[3]= 1.14895;   Fll[3]= 0.982813;  
    Fbb[4]= 0.926062;   Fcc[4]= 0.926062;   Fc[4]= 1.15437;   Fll[4]= 0.987455;  
  } 
  else if ( idsys==4   || sysname=="wt_down" ) { 
    Fbb[1]= 1.31379;   Fcc[1]= 1.31379;   Fc[1]= 1.06533;   Fll[1]= 0.971643;  
    Fbb[2]= 1.28735;   Fcc[2]= 1.28735;   Fc[2]= 1.04389;   Fll[2]= 0.952091;  
    Fbb[3]= 1.26841;   Fcc[3]= 1.26841;   Fc[3]= 1.02853;   Fll[3]= 0.938085;  
    Fbb[4]= 1.2488;   Fcc[4]= 1.2488;   Fc[4]= 1.01263;   Fll[4]= 0.923583;  
  } 
  else if ( idsys==5   || sysname=="st_up" ) { 
    Fbb[1]= 1.03345;   Fcc[1]= 1.03345;   Fc[1]= 1.11807;   Fll[1]= 0.9786;  
    Fbb[2]= 1.0272;   Fcc[2]= 1.0272;   Fc[2]= 1.11131;   Fll[2]= 0.972684;  
    Fbb[3]= 1.02487;   Fcc[3]= 1.02487;   Fc[3]= 1.10879;   Fll[3]= 0.970475;  
    Fbb[4]= 1.02346;   Fcc[4]= 1.02346;   Fc[4]= 1.10726;   Fll[4]= 0.969144;  
  } 
  else if ( idsys==6   || sysname=="st_down" ) { 
    Fbb[1]= 1.19793;   Fcc[1]= 1.19793;   Fc[1]= 1.09382;   Fll[1]= 0.97341;  
    Fbb[2]= 1.18046;   Fcc[2]= 1.18046;   Fc[2]= 1.07787;   Fll[2]= 0.959217;  
    Fbb[3]= 1.16903;   Fcc[3]= 1.16903;   Fc[3]= 1.06743;   Fll[3]= 0.949925;  
    Fbb[4]= 1.15763;   Fcc[4]= 1.15763;   Fc[4]= 1.05702;   Fll[4]= 0.940659;  
  } 
  else if ( idsys==7   || sysname=="z_up" ) { 
    Fbb[1]= 1.17024;   Fcc[1]= 1.17024;   Fc[1]= 1.12565;   Fll[1]= 0.969694;  
    Fbb[2]= 1.15374;   Fcc[2]= 1.15374;   Fc[2]= 1.10978;   Fll[2]= 0.956021;  
    Fbb[3]= 1.14381;   Fcc[3]= 1.14381;   Fc[3]= 1.10023;   Fll[3]= 0.947792;  
    Fbb[4]= 1.13434;   Fcc[4]= 1.13434;   Fc[4]= 1.09112;   Fll[4]= 0.939946;  
  } 
  else if ( idsys==8   || sysname=="z_down" ) { 
    Fbb[1]= 1.06729;   Fcc[1]= 1.06729;   Fc[1]= 1.08656;   Fll[1]= 0.981919;  
    Fbb[2]= 1.05993;   Fcc[2]= 1.05993;   Fc[2]= 1.07906;   Fll[2]= 0.975146;  
    Fbb[3]= 1.05605;   Fcc[3]= 1.05605;   Fc[3]= 1.07511;   Fll[3]= 0.971574;  
    Fbb[4]= 1.05265;   Fcc[4]= 1.05265;   Fc[4]= 1.07165;   Fll[4]= 0.96845;  
  } 
  else if ( idsys==9   || sysname=="di_up" ) { 
    Fbb[1]= 1.11377;   Fcc[1]= 1.11377;   Fc[1]= 1.10619;   Fll[1]= 0.976071;  
    Fbb[2]= 1.10237;   Fcc[2]= 1.10237;   Fc[2]= 1.09488;   Fll[2]= 0.966086;  
    Fbb[3]= 1.09584;   Fcc[3]= 1.09584;   Fc[3]= 1.08839;   Fll[3]= 0.960365;  
    Fbb[4]= 1.0898;   Fcc[4]= 1.0898;   Fc[4]= 1.08239;   Fll[4]= 0.955067;  
  } 
  else if ( idsys==10   || sysname=="di_down" ) { 
    Fbb[1]= 1.11708;   Fcc[1]= 1.11708;   Fc[1]= 1.10584;   Fll[1]= 0.975944;  
    Fbb[2]= 1.10545;   Fcc[2]= 1.10545;   Fc[2]= 1.09433;   Fll[2]= 0.965787;  
    Fbb[3]= 1.09874;   Fcc[3]= 1.09874;   Fc[3]= 1.08769;   Fll[3]= 0.959923;  
    Fbb[4]= 1.09249;   Fcc[4]= 1.09249;   Fc[4]= 1.08151;   Fll[4]= 0.954467;  
  } 
  else if ( idsys==11   || sysname=="qcd_up" ) { 
    Fbb[1]= 1.02494;   Fcc[1]= 1.02494;   Fc[1]= 1.08547;   Fll[1]= 0.984468;  
    Fbb[2]= 1.02041;   Fcc[2]= 1.02041;   Fc[2]= 1.08067;   Fll[2]= 0.980111;  
    Fbb[3]= 1.01869;   Fcc[3]= 1.01869;   Fc[3]= 1.07885;   Fll[3]= 0.978463;  
    Fbb[4]= 1.01764;   Fcc[4]= 1.01764;   Fc[4]= 1.07773;   Fll[4]= 0.977452;  
  } 
  else if ( idsys==12   || sysname=="qcd_down" ) { 
    Fbb[1]= 1.1986;   Fcc[1]= 1.1986;   Fc[1]= 1.12654;   Fll[1]= 0.967961;  
    Fbb[2]= 1.17974;   Fcc[2]= 1.17974;   Fc[2]= 1.10882;   Fll[2]= 0.952734;  
    Fbb[3]= 1.16807;   Fcc[3]= 1.16807;   Fc[3]= 1.09785;   Fll[3]= 0.943307;  
    Fbb[4]= 1.15677;   Fcc[4]= 1.15677;   Fc[4]= 1.08722;   Fll[4]= 0.934179;  
  } 
  else if ( idsys==25   || sysname=="preqcd_up" ) { 
    Fbb[1]= 1.26526;   Fcc[1]= 1.26526;   Fc[1]= 1.11217;   Fll[1]= 0.966608;  
    Fbb[2]= 1.24124;   Fcc[2]= 1.24124;   Fc[2]= 1.09106;   Fll[2]= 0.948254;  
    Fbb[3]= 1.22535;   Fcc[3]= 1.22535;   Fc[3]= 1.07709;   Fll[3]= 0.936112;  
    Fbb[4]= 1.20945;   Fcc[4]= 1.20945;   Fc[4]= 1.06312;   Fll[4]= 0.923969;  
  } 
  else if ( idsys==26   || sysname=="preqcd_down" ) { 
    Fbb[1]= 0.982929;   Fcc[1]= 0.982929;   Fc[1]= 1.09721;   Fll[1]= 0.984876;  
    Fbb[2]= 0.980563;   Fcc[2]= 0.980563;   Fc[2]= 1.09457;   Fll[2]= 0.982505;  
    Fbb[3]= 0.980783;   Fcc[3]= 0.980783;   Fc[3]= 1.09481;   Fll[3]= 0.982726;  
    Fbb[4]= 0.98196;   Fcc[4]= 0.98196;   Fc[4]= 1.09613;   Fll[4]= 0.983905;  
  } 
  else if ( idsys==13   || sysname=="jes_up" ) { 
    Fbb[1]= 1.35611;   Fcc[1]= 1.35611;   Fc[1]= 1.1252;   Fll[1]= 0.962165;  
    Fbb[2]= 1.32345;   Fcc[2]= 1.32345;   Fc[2]= 1.0981;   Fll[2]= 0.938988;  
    Fbb[3]= 1.30106;   Fcc[3]= 1.30106;   Fc[3]= 1.07952;   Fll[3]= 0.923107;  
    Fbb[4]= 1.27913;   Fcc[4]= 1.27913;   Fc[4]= 1.06133;   Fll[4]= 0.907549;  
  } 
  else if ( idsys==14   || sysname=="jes_down" ) { 
    Fbb[1]= 0.934617;   Fcc[1]= 0.934617;   Fc[1]= 1.06388;   Fll[1]= 0.992893;  
    Fbb[2]= 0.936451;   Fcc[2]= 0.936451;   Fc[2]= 1.06596;   Fll[2]= 0.994841;  
    Fbb[3]= 0.939091;   Fcc[3]= 0.939091;   Fc[3]= 1.06897;   Fll[3]= 0.997645;  
    Fbb[4]= 0.942797;   Fcc[4]= 0.942797;   Fc[4]= 1.07319;   Fll[4]= 1.00158;  
  } 
  else if ( idsys==67   || sysname=="eer_up" ) { 
    Fbb[1]= 1.11913;   Fcc[1]= 1.11913;   Fc[1]= 1.10335;   Fll[1]= 0.97623;  
    Fbb[2]= 1.10745;   Fcc[2]= 1.10745;   Fc[2]= 1.09184;   Fll[2]= 0.96604;  
    Fbb[3]= 1.10063;   Fcc[3]= 1.10063;   Fc[3]= 1.08512;   Fll[3]= 0.960096;  
    Fbb[4]= 1.09428;   Fcc[4]= 1.09428;   Fc[4]= 1.07885;   Fll[4]= 0.954553;  
  } 
  else if ( idsys==68   || sysname=="eer_down" ) { 
    Fbb[1]= 1.11509;   Fcc[1]= 1.11509;   Fc[1]= 1.10499;   Fll[1]= 0.976189;  
    Fbb[2]= 1.10364;   Fcc[2]= 1.10364;   Fc[2]= 1.09364;   Fll[2]= 0.966168;  
    Fbb[3]= 1.09703;   Fcc[3]= 1.09703;   Fc[3]= 1.08709;   Fll[3]= 0.960382;  
    Fbb[4]= 1.09094;   Fcc[4]= 1.09094;   Fc[4]= 1.08106;   Fll[4]= 0.955049;  
  } 
  else if ( idsys==65   || sysname=="ees_up" ) { 
    Fbb[1]= 1.11315;   Fcc[1]= 1.11315;   Fc[1]= 1.10589;   Fll[1]= 0.976151;  
    Fbb[2]= 1.1018;   Fcc[2]= 1.1018;   Fc[2]= 1.09462;   Fll[2]= 0.966206;  
    Fbb[3]= 1.09529;   Fcc[3]= 1.09529;   Fc[3]= 1.08815;   Fll[3]= 0.960493;  
    Fbb[4]= 1.08933;   Fcc[4]= 1.08933;   Fc[4]= 1.08223;   Fll[4]= 0.955263;  
  } 
  else if ( idsys==66   || sysname=="ees_down" ) { 
    Fbb[1]= 1.12176;   Fcc[1]= 1.12176;   Fc[1]= 1.10218;   Fll[1]= 0.976283;  
    Fbb[2]= 1.10994;   Fcc[2]= 1.10994;   Fc[2]= 1.09056;   Fll[2]= 0.965988;  
    Fbb[3]= 1.10297;   Fcc[3]= 1.10297;   Fc[3]= 1.08371;   Fll[3]= 0.959924;  
    Fbb[4]= 1.09645;   Fcc[4]= 1.09645;   Fc[4]= 1.07731;   Fll[4]= 0.954253;  
  } 
  else if ( idsys==63   || sysname=="jer_one" ) { 
    Fbb[1]= 1.12666;   Fcc[1]= 1.12666;   Fc[1]= 1.21951;   Fll[1]= 0.959133;  
    Fbb[2]= 1.10883;   Fcc[2]= 1.10883;   Fc[2]= 1.20022;   Fll[2]= 0.94396;  
    Fbb[3]= 1.10043;   Fcc[3]= 1.10043;   Fc[3]= 1.19113;   Fll[3]= 0.93681;  
    Fbb[4]= 1.09424;   Fcc[4]= 1.09424;   Fc[4]= 1.18443;   Fll[4]= 0.931542;  
  } 
  else if ( idsys==64   || sysname=="jef_one" ) { 
    Fbb[1]= 1.11448;   Fcc[1]= 1.11448;   Fc[1]= 1.10517;   Fll[1]= 0.976192;  
    Fbb[2]= 1.10307;   Fcc[2]= 1.10307;   Fc[2]= 1.09385;   Fll[2]= 0.966193;  
    Fbb[3]= 1.09652;   Fcc[3]= 1.09652;   Fc[3]= 1.08736;   Fll[3]= 0.960461;  
    Fbb[4]= 1.09043;   Fcc[4]= 1.09043;   Fc[4]= 1.08131;   Fll[4]= 0.955122;  
  } 
  else if ( idsys==61   || sysname=="musms_up" ) { 
    Fbb[1]= 1.11463;   Fcc[1]= 1.11463;   Fc[1]= 1.10419;   Fll[1]= 0.976347;  
    Fbb[2]= 1.10325;   Fcc[2]= 1.10325;   Fc[2]= 1.09292;   Fll[2]= 0.966378;  
    Fbb[3]= 1.09668;   Fcc[3]= 1.09668;   Fc[3]= 1.08641;   Fll[3]= 0.960621;  
    Fbb[4]= 1.09062;   Fcc[4]= 1.09062;   Fc[4]= 1.08041;   Fll[4]= 0.955319;  
  } 
  else if ( idsys==62   || sysname=="musms_down" ) { 
    Fbb[1]= 1.12135;   Fcc[1]= 1.12135;   Fc[1]= 1.10343;   Fll[1]= 0.976099;  
    Fbb[2]= 1.10951;   Fcc[2]= 1.10951;   Fc[2]= 1.09177;   Fll[2]= 0.965785;  
    Fbb[3]= 1.10255;   Fcc[3]= 1.10255;   Fc[3]= 1.08493;   Fll[3]= 0.959731;  
    Fbb[4]= 1.09604;   Fcc[4]= 1.09604;   Fc[4]= 1.07852;   Fll[4]= 0.954065;  
  } 
  else if ( idsys==57   || sysname=="mscale_up" ) { 
    Fbb[1]= 1.11293;   Fcc[1]= 1.11293;   Fc[1]= 1.10548;   Fll[1]= 0.97623;  
    Fbb[2]= 1.10162;   Fcc[2]= 1.10162;   Fc[2]= 1.09425;   Fll[2]= 0.966308;  
    Fbb[3]= 1.09513;   Fcc[3]= 1.09513;   Fc[3]= 1.08781;   Fll[3]= 0.960622;  
    Fbb[4]= 1.08912;   Fcc[4]= 1.08912;   Fc[4]= 1.08184;   Fll[4]= 0.955347;  
  } 
  else if ( idsys==58   || sysname=="mscale_down" ) { 
    Fbb[1]= 1.11173;   Fcc[1]= 1.11173;   Fc[1]= 1.10659;   Fll[1]= 0.976122;  
    Fbb[2]= 1.10046;   Fcc[2]= 1.10046;   Fc[2]= 1.09537;   Fll[2]= 0.966225;  
    Fbb[3]= 1.09403;   Fcc[3]= 1.09403;   Fc[3]= 1.08897;   Fll[3]= 0.960582;  
    Fbb[4]= 1.08806;   Fcc[4]= 1.08806;   Fc[4]= 1.08304;   Fll[4]= 0.955343;  
  } 
  else if ( idsys==55   || sysname=="metpileup_up" ) { 
    Fbb[1]= 1.12248;   Fcc[1]= 1.12248;   Fc[1]= 1.10648;   Fll[1]= 0.975538;  
    Fbb[2]= 1.11046;   Fcc[2]= 1.11046;   Fc[2]= 1.09463;   Fll[2]= 0.965085;  
    Fbb[3]= 1.10339;   Fcc[3]= 1.10339;   Fc[3]= 1.08766;   Fll[3]= 0.958943;  
    Fbb[4]= 1.09682;   Fcc[4]= 1.09682;   Fc[4]= 1.08118;   Fll[4]= 0.953233;  
  } 
  else if ( idsys==56   || sysname=="metpileup_down" ) { 
    Fbb[1]= 1.11925;   Fcc[1]= 1.11925;   Fc[1]= 1.10431;   Fll[1]= 0.976068;  
    Fbb[2]= 1.1075;   Fcc[2]= 1.1075;   Fc[2]= 1.09272;   Fll[2]= 0.965821;  
    Fbb[3]= 1.10067;   Fcc[3]= 1.10067;   Fc[3]= 1.08598;   Fll[3]= 0.959864;  
    Fbb[4]= 1.09432;   Fcc[4]= 1.09432;   Fc[4]= 1.07971;   Fll[4]= 0.954325;  
  } 
  else if ( idsys==53   || sysname=="met_up" ) { 
    Fbb[1]= 1.11866;   Fcc[1]= 1.11866;   Fc[1]= 1.10423;   Fll[1]= 0.976115;  
    Fbb[2]= 1.10698;   Fcc[2]= 1.10698;   Fc[2]= 1.0927;   Fll[2]= 0.965923;  
    Fbb[3]= 1.10018;   Fcc[3]= 1.10018;   Fc[3]= 1.08599;   Fll[3]= 0.959988;  
    Fbb[4]= 1.09386;   Fcc[4]= 1.09386;   Fc[4]= 1.07975;   Fll[4]= 0.954471;  
  } 
  else if ( idsys==54   || sysname=="met_down" ) { 
    Fbb[1]= 1.12529;   Fcc[1]= 1.12529;   Fc[1]= 1.1018;   Fll[1]= 0.976141;  
    Fbb[2]= 1.11319;   Fcc[2]= 1.11319;   Fc[2]= 1.08996;   Fll[2]= 0.965651;  
    Fbb[3]= 1.10602;   Fcc[3]= 1.10602;   Fc[3]= 1.08294;   Fll[3]= 0.959427;  
    Fbb[4]= 1.09928;   Fcc[4]= 1.09928;   Fc[4]= 1.07634;   Fll[4]= 0.95358;  
  } 
  else if ( idsys==69   || sysname=="powhe_one" ) { 
    Fbb[1]= 1.07165;   Fcc[1]= 1.07165;   Fc[1]= 1.11918;   Fll[1]= 0.97628;  
    Fbb[2]= 1.06278;   Fcc[2]= 1.06278;   Fc[2]= 1.10992;   Fll[2]= 0.968203;  
    Fbb[3]= 1.05847;   Fcc[3]= 1.05847;   Fc[3]= 1.10542;   Fll[3]= 0.964279;  
    Fbb[4]= 1.05495;   Fcc[4]= 1.05495;   Fc[4]= 1.10173;   Fll[4]= 0.961065;  
  } 
  else if ( idsys==70   || sysname=="powpy_one" ) { 
    Fbb[1]= 1.11086;   Fcc[1]= 1.11086;   Fc[1]= 1.10206;   Fll[1]= 0.976918;  
    Fbb[2]= 1.09983;   Fcc[2]= 1.09983;   Fc[2]= 1.09112;   Fll[2]= 0.967224;  
    Fbb[3]= 1.09349;   Fcc[3]= 1.09349;   Fc[3]= 1.08483;   Fll[3]= 0.961649;  
    Fbb[4]= 1.08761;   Fcc[4]= 1.08761;   Fc[4]= 1.07899;   Fll[4]= 0.956472;  
  } 
  else if ( idsys==15   || sysname=="btag_up" ) { 
    Fbb[1]= 0.371277;   Fcc[1]= 0.371277;   Fc[1]= 1.24659;   Fll[1]= 0.994378;  
    Fbb[2]= 0.381941;   Fcc[2]= 0.381941;   Fc[2]= 1.2824;   Fll[2]= 1.02294;  
    Fbb[3]= 0.393381;   Fcc[3]= 0.393381;   Fc[3]= 1.32081;   Fll[3]= 1.05358;  
    Fbb[4]= 0.407979;   Fcc[4]= 0.407979;   Fc[4]= 1.36982;   Fll[4]= 1.09268;  
  } 
  else if ( idsys==16   || sysname=="btag_down" ) { 
    Fbb[1]= 2.03578;   Fcc[1]= 2.03578;   Fc[1]= 0.956742;   Fll[1]= 0.949219;  
    Fbb[2]= 1.92288;   Fcc[2]= 1.92288;   Fc[2]= 0.903682;   Fll[2]= 0.896576;  
    Fbb[3]= 1.83825;   Fcc[3]= 1.83825;   Fc[3]= 0.863906;   Fll[3]= 0.857113;  
    Fbb[4]= 1.75203;   Fcc[4]= 1.75203;   Fc[4]= 0.823388;   Fll[4]= 0.816913;  
  } 
  else if ( idsys==17   || sysname=="bmtag_up" ) { 
    Fbb[1]= 1.14646;   Fcc[1]= 1.14646;   Fc[1]= 1.09572;   Fll[1]= 0.975975;  
    Fbb[2]= 1.13301;   Fcc[2]= 1.13301;   Fc[2]= 1.08286;   Fll[2]= 0.964526;  
    Fbb[3]= 1.12467;   Fcc[3]= 1.12467;   Fc[3]= 1.07489;   Fll[3]= 0.957427;  
    Fbb[4]= 1.11658;   Fcc[4]= 1.11658;   Fc[4]= 1.06716;   Fll[4]= 0.950539;  
  } 
  else if ( idsys==18   || sysname=="bmtag_down" ) { 
    Fbb[1]= 1.08464;   Fcc[1]= 1.08464;   Fc[1]= 1.11625;   Fll[1]= 0.976037;  
    Fbb[2]= 1.07497;   Fcc[2]= 1.07497;   Fc[2]= 1.1063;   Fll[2]= 0.967337;  
    Fbb[3]= 1.06998;   Fcc[3]= 1.06998;   Fc[3]= 1.10117;   Fll[3]= 0.962851;  
    Fbb[4]= 1.0657;   Fcc[4]= 1.0657;   Fc[4]= 1.09676;   Fll[4]= 0.958994;  
  } 
  
  else { 
    cout << " WARNING: Ffactors request for unknown variation. Set to Nominal! " << endl; 
  }
  return;
}
